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

#include <iostream>

#include <unordered_map>
#include <unordered_set>

#include <boost/progress.hpp>
using progress = boost::progress_display;

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

#include <boost/thread/mutex.hpp>

#include <tbb/tbb.h>
#include <tbb/cache_aligned_allocator.h>

#include <cstdio>
#include <sys/stat.h>

#include "zlib.h"

using namespace sina;

static const char* module_name = "kmer_search";
static auto logger = Log::create_logger(module_name);

std::map<string, kmer_search*> kmer_search::indices;
static boost::mutex indices_access;

kmer_search*
kmer_search::get_kmer_search(const fs::path& filename, int k) {
    boost::mutex::scoped_lock lock(indices_access);
    std::stringstream str;
    str << filename << "%%k=" << k;
    if (indices.empty()) {
        atexit(kmer_search::destroy_indices);
    }
    if (indices.count(str.str()) == 0u) {
        indices[str.str()] = new kmer_search(query_arb::getARBDB(filename), k);
    }
    return indices[str.str()];
}

void
kmer_search::destroy_indices() {
    for (auto& pair: indices) {
        delete pair.second;
    }
}

class kmer_search::index {
public:
    friend class kmer_search;
    friend class BuildIndex;
    int k;
    int n_kmers;
    int n_sequences{0};

    std::vector<std::string> sequence_names;
    std::vector<vlimap*> kmer_idx;

    query_arb* arbdb;
    timer timeit;


    index(int k_, query_arb* arbdb_)
        : k(k_),
          n_kmers(1<<(k_*2)),
          kmer_idx(1<<(k_*2), nullptr),
          arbdb(arbdb_)
    {}
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
}

class BuildIndex {
    kmer_search::index *data;
    boost::progress_display *p;
public:
    std::vector<vlimap*> kmer_idx;

    void operator()(const tbb::blocked_range<size_t>& r) {
        size_t end = r.end();
        std::unordered_set<unsigned int> seen;
        for (size_t i = r.begin(); i < end; ++i) {
            const cseq& c = data->arbdb->getCseqUncached(data->sequence_names[i]);
            const auto& bases = c.const_getAlignedBases();
            for (unsigned int kmer: unique_prefix_kmers(bases, seen, data->k, 1, BASE_A)) {
                //for (const auto& kmer: unique_kmers(bases, seen, data->k)) {
                if (kmer_idx[kmer] == nullptr) {
                    kmer_idx[kmer] = new vlimap(data->n_sequences);
                }
                kmer_idx[kmer]->push_back(i);
            }
            ++(*p);
        }
    }

    void join(BuildIndex& other) {
        for (int i = 0; i < data->n_kmers; ++i) {
            if (other.kmer_idx[i] != nullptr) {
                if (kmer_idx[i] != nullptr) {
                    kmer_idx[i]->append(*other.kmer_idx[i]);
                } else {
                    kmer_idx[i] = other.kmer_idx[i];
                    other.kmer_idx[i] = nullptr;
                }
            }
        }
    }

    ~BuildIndex() {
        for (int i = 0; i < data->n_kmers; ++i) {
            delete kmer_idx[i];
        }
    }

    BuildIndex(BuildIndex& x, tbb::split)
        : data(x.data), p(x.p), kmer_idx(x.data->n_kmers, nullptr)
    {
    }

    BuildIndex(kmer_search::index *data_, boost::progress_display *p_)
        : data(data_), p(p_), kmer_idx(data->n_kmers, nullptr)
    {
    }
};

void
kmer_search::build_index() {
    timer t;
    t.start();

    logger->info("Loading names...");
    data.sequence_names = data.arbdb->getSequenceNames();
    data.n_sequences = data.sequence_names.size();

    logger->info("Building index...");
    {
        boost::progress_display p(data.n_sequences, std::cerr);
        BuildIndex bi(&data, &p);
        tbb::parallel_reduce(tbb::blocked_range<size_t>(0, data.n_sequences), bi);
        data.kmer_idx.clear();
        data.kmer_idx.reserve(data.n_kmers);
        for (int i=0; i < data.n_kmers; i++) {
            if (bi.kmer_idx[i] != nullptr) {
                data.kmer_idx.push_back(new vlimap(*bi.kmer_idx[i]));
            } else {
                data.kmer_idx.push_back(new vlimap(data.n_sequences));
            }
        }
    }
    t.stop("build");

    for (int i=0; i < data.n_kmers; i++) {
        if (data.kmer_idx[i]->size() > data.n_sequences / 2) {
            data.kmer_idx[i]->invert();
        }
    }
    t.stop("shrink");

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
                   int  /*min_match*/,
                   int max_match,
                   float  /*min_score*/,
                   float  /*max_score*/,
                   query_arb*  /*arb*/,
                   bool  /*noid*/,
                   int  /*min_len*/,
                   int  /*num_full*/,
                   int  /*full_min_len*/,
                   int  /*range_cover*/,
                   bool  /*leave_query_out*/) {
    find(query, results, max_match);
    if (results.empty()) {
        return 0;
    }
    return results[0].getScore();
}

void
kmer_search::find(const cseq& query, std::vector<cseq>& results, int max) {
    if (data.n_sequences == 0) {
        return;
    }
    data.timeit.start();
    const vector<aligned_base>& bases = query.const_getAlignedBases();
    idset::inc_t scores(data.n_sequences, 0);
    data.timeit.stop("load");

    std::unordered_set<unsigned int> seen(query.size()*2-1);

    int offset = 0;
    // for (unsigned int kmer: all_kmers(bases, data.k)) {
    // for (unsigned int kmer: unique_kmers(bases, seen, data.k)) {
    // for (unsigned int kmer: unique_prefix_kmers(bases, seen, data.k, 1, BASE_A)) {
    for (unsigned int kmer: prefix_kmers(bases, data.k, 1, BASE_A)) {
        if (data.kmer_idx[kmer]->size() < data.n_sequences) {
            offset += data.kmer_idx[kmer]->increment(scores);
        }
    }
    data.timeit.stop("sum");

    using pair = std::pair<idset::inc_t::value_type, int>;
    std::vector<pair> ranks;
    ranks.reserve(data.n_sequences);
    int n = 0;
    for (auto score: scores) {
        ranks.emplace_back(score + offset, n++);
    }
    data.timeit.stop("assemble");

    std::partial_sort(ranks.begin(), ranks.begin()+max, ranks.end(),
                      std::greater<pair>());
    data.timeit.stop("sort");

    results.clear();
    results.reserve(max);
    for (int i=0; i<max; i++) {
        results.emplace_back(data.arbdb->getCseq(data.sequence_names[ranks[i].second]));
        results.back().setScore(ranks[i].first);
    }
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
