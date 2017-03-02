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

#ifndef _KMER_SEARCH_H_
#define _KMER_SEARCH_H_

#include <utility> // std::pair
#include <iterator>
#include <vector>
#include <string>

namespace sina {

class query_arb;
class cseq;

class kmer_search {
public:
    class result_iterator;

    kmer_search();
    ~kmer_search();
    
    void 
    build_index(query_arb& db);
    
    std::pair<result_iterator, result_iterator>
    find(const cseq& query, unsigned int max);

private:
    std::vector<std::string> seqNames;
    std::vector<std::vector<unsigned int> > map;
};

class kmer_search::result_iterator {
public:
    result_iterator();
    ~result_iterator();

    typedef cseq value_type;
    typedef cseq* pointer;
    typedef cseq& reference;
    typedef std::forward_iterator_tag iterator_category;
    typedef int difference_type;

    cseq& operator*();
    result_iterator& operator++();
    bool operator==(const result_iterator& rhs) const;
};


} // namespace sina

#endif // _KMER_SEARCH_H_


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
