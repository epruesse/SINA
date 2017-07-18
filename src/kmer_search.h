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

    kmer_search(int k=8);
    ~kmer_search();
    
    void 
    build_index(query_arb& db);

    void
    find(const cseq& query, std::vector<cseq>& family, int max);

private:
    class index;
    index &data;
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
