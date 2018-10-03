/*
Copyright (c) 2006-2018 Elmar Pruesse <elmar.pruesse@ucdenver.edu>

This file is part of SINA.
SINA is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your
option) any later version.

SINA is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with SINA.  If not, see <http://www.gnu.org/licenses/>.

Additional permission under GNU GPL version 3 section 7

If you modify SINA, or any covered work, by linking or combining it
with components of ARB (or a modified version of that software),
containing parts covered by the terms of the
ARB-public-library-license, the licensors of SINA grant you additional
permission to convey the resulting work. Corresponding Source for a
non-source form of such a combination shall include the source code
for the parts of ARB used as well as that of the covered work.
*/

#include "kmer_search.h"
#include "kmer.h"
#include "idset.h"
#include "query_arb.h"
#include "helpers.h"
#include "timer.h"
#include "log.h"

#include <vector>
using std::vector;

#include <string>
using std::string;

#include <algorithm>
using std::sort;

#include <iostream>
#include <iomanip>
using std::endl;
using std::uppercase;
using std::hex;

#include <unordered_map>
#include <unordered_set>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <boost/progress.hpp>
using progress = boost::progress_display;

#include <boost/thread/mutex.hpp>

#include <stdio.h>
#include <sys/stat.h>

#include "zlib.h"

using namespace sina;

static const char* module_name = "kmer_search";
static auto logger = Log::create_logger(module_name);

std::map<string, kmer_search*> kmer_search::indices;
static boost::mutex indices_access;

kmer_search*
kmer_search::get_kmer_search(std::string filename, int k) {
    boost::mutex::scoped_lock lock(indices_access);
    std::stringstream str;
    str << filename << "%%k=" << k;
    if (indices.size() == 0) {
        atexit(kmer_search::destroy_indices);
    }
    if (not indices.count(str.str())) {
        indices[str.str()] = new kmer_search(query_arb::getARBDB(filename), k);
    }
    return indices[str.str()];
}

void
kmer_search::destroy_indices() {
    for (const std::pair<std::string, kmer_search*>& pair: indices) {
        delete pair.second;
    }
}

class kmer_search::index {
    friend class kmer_search;
    int k;
    int n_kmers;
    int n_sequences;

    std::vector<std::string> sequence_names;
    std::vector<idset*> kmer_idx;

    query_arb* arbdb;
    timer timeit;

public:
    index(int k_, query_arb* arbdb_) :
        k(k_),
        n_kmers(1<<(k_*2)),
        n_sequences(0),
        sequence_names(),
        kmer_idx(1<<(k_*2), NULL),
        arbdb(arbdb_),
        timeit()
    {
    }
    ~index() {
        logger->info("Timings for Kmer Search: {}", timeit);
    }
};

const uint64_t idx_header_magic = 0x53494e414b494458; // SINAKIDX
const uint16_t idx_header_vers  = 0;
struct idx_header {
    uint64_t magic;
    uint16_t vers;  // 0
    uint16_t k;
    uint32_t n_sequences;
    uint64_t kmer_idx_size;
    uint64_t total_size;
};

const uint64_t chunk_size = 1024*1024;

kmer_search::kmer_search(query_arb* arbdb, int k)
    : data(* new index(k, arbdb))
{
    init();
}

kmer_search::~kmer_search() {
    delete &data;
}

void
kmer_search::init() {
    build_index();
    /*
    logger->info("Trying to load index from disk");
    if (data.try_load_index()) {
        logger->info("Loaded index");
        return;
    } else {
        build_index();
        data.save_index();
    }
    */
}

void
kmer_search::build_index() {
    data.arbdb->loadCache(); // parallel load - 10% faster overall
    data.sequence_names = data.arbdb->getSequenceNames();
    data.n_sequences = data.sequence_names.size();

    logger->info("Building index...");
    boost::progress_display p(data.n_sequences, std::cerr);
    timer t;
    t.start();
    data.kmer_idx.clear();
    data.kmer_idx.reserve(data.n_kmers);
    for (int i=0; i < data.n_kmers; i++) {
        data.kmer_idx.push_back(new vlimap(data.n_sequences));
    }
    t.stop("alloc");

    std::unordered_set<unsigned int> seen;
    for (int i=0; i < data.n_sequences; i++) {
        const cseq& c = data.arbdb->getCseq(data.sequence_names[i]);
        const auto& bases = c.const_getAlignedBases();
        t.stop("load");
        for (const auto& kmer: unique_kmers(bases, seen, data.k)) {
            data.kmer_idx[kmer]->push_back(i);
        }
        t.stop("hash");
        t.end_loop(2);
        ++p;
    }
    t.stop("count");

    logger->info("Indexed {} sequences", data.n_sequences);
    logger->info("Timings: {}", t);
}

struct score {
    unsigned int val;
    unsigned int id;
    bool operator<(const score& rhs) const {
        return val > rhs.val;
    }
};

double
kmer_search::match(std::vector<cseq>& results,
                   const cseq& query,
                   int min_match,
                   int max_match,
                   float min_score,
                   float max_score,
                   query_arb* arb,
                   bool noid,
                   int min_len,
                   int num_full,
                   int full_min_len,
                   int range_cover,
                   bool leave_query_out) {
    find(query, results, max_match);
    if (results.size() == 0) {
        return 0;
    } else {
        return results[0].getScore();
    }
}

void
kmer_search::find(const cseq& query, std::vector<cseq>& results, int max) {
    if (data.n_sequences == 0) {
        return;
    }
    data.timeit.start();
    const vector<aligned_base>& bases = query.const_getAlignedBases();
    vector<uint16_t> scores(data.n_sequences, 0);
    data.timeit.stop("load");

    std::unordered_set<unsigned int> seen;

    for (unsigned int kmer: unique_kmers(bases, seen, data.k)) {
        data.kmer_idx[kmer]->increment(scores);
    }
    data.timeit.stop("sum");

    std::vector<std::pair<int, string> > scored_names;
    scored_names.reserve(data.n_sequences);
    std::transform(scores.begin(), scores.end(),
                   data.sequence_names.begin(),
                   std::back_inserter(scored_names),
                   [](int score, string name) {
                       return std::make_pair(score, name);
                   }
        );
    data.timeit.stop("assemble");

    std::partial_sort(scored_names.begin(), scored_names.begin()+max, scored_names.end(),
                      std::greater<std::pair<int,string> >());
    data.timeit.stop("sort");

    results.clear();
    results.reserve(max);
    std::transform(scored_names.begin(), scored_names.begin()+max,
                   std::back_inserter(results),
                   [&] (std::pair<int, string> hit) {
                       cseq c = data.arbdb->getCseq(hit.second);
                       c.setScore(hit.first);
                       return c;
                   });
    data.timeit.stop("output");
}

/*
  Local Variables:
  mode:c++
  c-file-style:"stroustrup"
  c-file-offsets:((innamespace . 0)(inline-open . 0)(case-label . +))
  indent-tabs-mode:nil
  fill-column:99
  End:
*/
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:encoding=utf-8:textwidth=99 :
