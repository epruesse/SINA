/*
Copyright (c) 2006-2017 Elmar Pruesse <elmar.pruesse@ucdenver.edu>

This file is part of SINA.

SINA is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "kmer_search.h"
#include "query_arb.h"

#include <vector>
using std::vector;

#include <string>
using std::string;

#include <algorithm>
using std::sort;

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH



using namespace sina;

class hash_fourbase {
public:
    static const unsigned short K = 10;
    hash_fourbase() 
        : val(0), good_count(0)
    {
    }
    
    void
    push(const aligned_base& ab) {
        if (ab.is_ambig()) {
            good_count = 0;
        } else {
            good_count++;
            val <<= 2;
            val &= 0xFFFFF;
            val += ab;
        }
    }
    
    bool good() {
        return good_count >= K;
    }
    
    unsigned int
    value() const {
        return val;
    }
    
    unsigned int 
    max() const {
        return 0xFFFFF;
    }
private:
    unsigned int val;
    unsigned int good_count;
};

kmer_search::kmer_search() {
}

kmer_search::~kmer_search() {
}

void
kmer_search::build_index(query_arb& arbdb) {
    seqNames = arbdb.getSequenceNames();
    hash_fourbase hash;
    map.resize(hash.max()+1);
    
    int seqno = 0;
    foreach(string name, seqNames) {
        const cseq& c = arbdb.getCseq(name);
        const vector<aligned_base>& bases = c.const_getAlignedBases();
        foreach (const aligned_base& ab, bases) {
            hash.push(ab);
            if (hash.good()) {
                map[hash.value()].push_back(seqno);
            }
        }
        seqno++;
    }
    
    foreach(vector<unsigned int>& l, map) {
        std::cout << l.size() << "," << l.capacity() << std::endl;
    }
    
}

struct score {
    unsigned int val;
    unsigned int id;
    bool operator<(const score& rhs) const {
        return val > rhs.val;
    }
};

std::pair<kmer_search::result_iterator, 
	  kmer_search::result_iterator>
kmer_search::find(const cseq& query, unsigned int /*max*/) {

    vector<score> scores;
    const int num_seqs = seqNames.size();
    scores.resize(num_seqs);
    for (int i=0; i < num_seqs; i++) {
        scores[i].id=i;
    }
    
    hash_fourbase hash;
    const vector<aligned_base>& bases = query.const_getAlignedBases();
    foreach (aligned_base ab, bases) {
        hash.push(ab);
        if (hash.good()) {
            foreach(int j, map[hash.value()]) {
                scores[j].val++;
            }
        }
    }
    
    std::partial_sort(scores.begin(), scores.begin()+40, scores.end());
    
    /*
      std::cerr << query.getName() << ": ";
      for (int i=0; i<40; i++) {
      std::cerr << seqNames[scores[i].id] << " ";
      }
      std::cerr << std::endl;
    */
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
