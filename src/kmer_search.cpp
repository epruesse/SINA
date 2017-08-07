/*
Copyright (c) 2006-2017 Elmar Pruesse <elmar.pruesse@ucdenver.edu>

This file is part of SINA.
123456789012345678901234567890123456789012345678901234567890123456789012
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
#include "query_arb.h"
#include "helpers.h"

#include <vector>
using std::vector;

#include <string>
using std::string;

#include <algorithm>
using std::sort;

#include <iostream>
#include <iomanip>
using std::cerr;
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


template<typename T>
class const_range {
public:
    typedef typename T::const_iterator const_iterator;
    const_range(const const_iterator & begin, const const_iterator & end) :
        _begin(begin), _end(end) {}
    const_iterator begin() const { return _begin; }
    const_iterator end() const { return _end; }
private:
    const_iterator _begin, _end;
};
template<typename T>
const_range<T>
make_const_range(const T &t, size_t offset, size_t length) {
    return const_range<T>(t.begin() + offset, t.begin() + offset + length);
}

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

    /* database kmer index:
     *
     * sequence_names: list of sequence IDs
     * kmer_list:      references into sequence_names
     * kmer_offsets:   references into list
     * kmer_counts:    number of occurence of kmer
     */

    typedef typename std::vector<unsigned int> idx_type;
    std::vector<std::string> sequence_names;
    idx_type kmer_idx;
    idx_type kmer_counts;
    idx_type kmer_offsets;

    query_arb* arbdb;

public:
    index(int k_, query_arb* arbdb_) :
        k(k_),
        n_kmers(1<<(k_*2)),
        n_sequences(0),
        sequence_names(),
        kmer_counts(1<<(k_*2), 0),
        kmer_offsets(1<<(k*2)+1),
        kmer_idx(),
        arbdb(arbdb_)
    {
    }

    const_range<idx_type> get_matching_seqnum(unsigned int kmer) {
        return make_const_range(kmer_idx, kmer_offsets[kmer], kmer_counts[kmer]);
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
    cerr << "Trying to load index from disk..." << endl;
    if (data.try_load_index()) {
        cerr << "Loaded index" << endl;
        return;
    } else {
        build_index();
        data.save_index();
    }
    */
}

void
kmer_search::build_index() {
    data.sequence_names = data.arbdb->getSequenceNames();
    
    // count in how many sequences each kmer occurs
    cerr << "Loading sequences and counting kmers..." << endl;
    {
        boost::progress_display p(data.sequence_names.size(), cerr);
        for (string& name : data.sequence_names) {
            const cseq& c = data.arbdb->getCseq(name);
            const vector<aligned_base>& bases = c.const_getAlignedBases();
            for (unsigned int kmer: unique_kmers(bases, data.k)) {
                ++data.kmer_counts[kmer];
            }
            ++data.n_sequences;
            ++p;
        }
    }

    cerr << "Calculating offsets..." << endl;
    // calculate offsets
    unsigned long total = 0;
    for (unsigned int i = 0; i < data.n_kmers; ++i) {
        data.kmer_offsets[i] = total;
        total += data.kmer_counts[i];
    }
    data.kmer_offsets[data.n_kmers] = total;
    
    // reset counts
    std::fill(data.kmer_counts.begin(), data.kmer_counts.end(), 0);

    cerr << "Filling index..." << endl;
    // fill indices
    data.kmer_idx.resize(total);
    boost::progress_display p(data.sequence_names.size(), cerr);
    for (int i=0; i<data.n_sequences; ++i) {
        const cseq& c = data.arbdb->getCseq(data.sequence_names[i]);
        const vector<aligned_base>& bases = c.const_getAlignedBases();
        for (unsigned int kmer : unique_kmers(bases, data.k)) {
            data.kmer_idx[data.kmer_offsets[kmer] +
                          data.kmer_counts[kmer]++] = i;
        }
        ++p;
    }

    cerr << "Index done" << endl;
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
    kmer_generator kmer(data.k);
    const vector<aligned_base>& bases = query.const_getAlignedBases();
    vector<int> scores(data.n_sequences, 0);

    for (unsigned int kmer: unique_kmers(bases, data.k)) {
        for (int idx : data.get_matching_seqnum(kmer)) {
            scores[idx]++;
        }
    }
    std::vector<std::pair<int, string> > scored_names;
    scored_names.reserve(data.n_sequences);
    std::transform(scores.begin(), scores.end(),
                   data.sequence_names.begin(),
                   std::back_inserter(scored_names),
                   [](int score, string name) {
                       return std::make_pair(score, name);
                   }
        );
    std::partial_sort(scored_names.begin(), scored_names.begin()+max, scored_names.end(),
                      std::greater<std::pair<int,string> >());
    results.clear();
    results.reserve(max);
    std::transform(scored_names.begin(), scored_names.begin()+max,
                   std::back_inserter(results),
                   [&] (std::pair<int, string> hit) {
                       cseq c = data.arbdb->getCseq(hit.second);
                       c.setScore(hit.first);
                       return c;
                   });
    
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
