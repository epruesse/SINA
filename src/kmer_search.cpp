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

using namespace sina;

class kmer_generator {
private:
    const unsigned int k;
    const unsigned int mask;
    unsigned int val;
    unsigned int good_count;
public:
    kmer_generator(int k)
        : val(0), good_count(0), k(10), mask((1UL<<(2*k))-1)
    {
        if (sizeof(val)*8 < 2*k) {
            throw std::runtime_error("K too large!");
        }
    }
    
    void push(const base_iupac& b) {
        if (unlikely(b.is_ambig())) {
            good_count = 0;
        } else {
            good_count++;
            val <<= 2;
            val &= mask;
            val += b.getBaseType();
        }
    }
    
    bool good() const {
        return good_count >= k;
    }
    
    operator unsigned int() const {
        return val;
    }
};

class unique_kmers {
    typedef const std::vector<aligned_base> bases;
    typedef bases::const_iterator base_iterator;
public:
    unique_kmers(bases &v, int k) :
        _begin(v.begin()), _end(v.end()), _k(k) {}

    class const_iterator {
    public:
        const_iterator() : _k(10) {}
        const_iterator(base_iterator& begin, base_iterator& end, int k) :
            _k(k), _begin(begin), _end(end) {
            seen.reserve(_begin - _end);
            this->operator++();
        }
        unsigned int operator*() const {
            return _k;
        }
        const_iterator& operator++() {
            _k.push(*_begin++);
            while (likely(_begin != _end)
                   &&
                   not (_k.good()
                        &&
                        likely(seen.insert(_k).second)
                       )) {
                _k.push(*_begin++);
            }
        }
        bool operator!=(const const_iterator& rhs) const {
            return _begin != _end;
        }
    private:
        kmer_generator _k;
        base_iterator _begin, _end;
        std::unordered_set<unsigned int> seen;
    };

    const_iterator begin() {
        return const_iterator(_begin, _end, _k);
    }
    const_iterator end() {
        return const_iterator();
    }


private:
    base_iterator _begin, _end;
    unsigned int _k;
};

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
    index(int k_) :
        k(k_),
        n_kmers(1<<k*2),
        n_sequences(0),
        sequence_names(),
        kmer_counts(n_kmers, 0),
        kmer_offsets(n_kmers+1),
        kmer_idx(),
        arbdb(0)
    {
    }

    const_range<idx_type> get_matching_seqnum(unsigned int kmer) {
        return make_const_range(kmer_idx, kmer_offsets[kmer], kmer_counts[kmer]);
    }
};

kmer_search::kmer_search(int k)
    : data(* new index(k))
{
}

kmer_search::~kmer_search() {
    delete &data;
}

void
kmer_search::build_index(query_arb& arbdb) {
    data.arbdb = &arbdb;
    data.sequence_names = arbdb.getSequenceNames();
    
    // count in how many sequences each kmer occurs
    cerr << "Loading sequences and counting kmers..." << endl;
    {
        boost::progress_display p(data.sequence_names.size(), cerr);
        for (string& name : data.sequence_names) {
            const cseq& c = arbdb.getCseq(name);
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
        const cseq& c = arbdb.getCseq(data.sequence_names[i]);
        const vector<aligned_base>& bases = c.const_getAlignedBases();
        for (unsigned int kmer : unique_kmers(bases, data.k)) {
            data.kmer_idx[data.kmer_offsets[kmer] +
                          data.kmer_counts[kmer]++] = i;
        }
        ++p;
    }
}

struct score {
    unsigned int val;
    unsigned int id;
    bool operator<(const score& rhs) const {
        return val > rhs.val;
    }
};

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
    query_arb* arbdb = data.arbdb;
    std::transform(scored_names.begin(), scored_names.begin()+max,
                   std::back_inserter(results),
                   [&] (std::pair<int, string> hit) {
                       cseq c = arbdb->getCseq(hit.second);
                       c.setScore(hit.first);
                       return c;
                   });
    
}


kmer_search::result_iterator::result_iterator() {
}

kmer_search::result_iterator::~result_iterator() {
}

cseq&
kmer_search::result_iterator::operator*() {
}

kmer_search::result_iterator&
kmer_search::result_iterator::operator++() {
}

bool
kmer_search::result_iterator::operator==(const result_iterator& /*rhs*/ ) const {
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
