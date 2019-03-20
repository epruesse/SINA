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

#ifndef _SEARCH_H_
#define _SEARCH_H_

#include <boost/core/noncopyable.hpp>
#include <boost/program_options.hpp>
#include <string>
#include <vector>
#include "cseq.h"

namespace sina {

class query_arb;

enum ENGINE_TYPE {
    ENGINE_ARB_PT=0,
    ENGINE_SINA_KMER=1
};
std::ostream& operator<<(std::ostream& out, const sina::ENGINE_TYPE& t);
void validate(boost::any& v, const std::vector<std::string>& values,
              sina::ENGINE_TYPE* /*unused*/,int /*unused*/);


class search : private boost::noncopyable {
protected:
    search() = default;
public:
    virtual ~search() = default;

    /**
     * match runs a word search using the PT server
     *
     * arguments:
     *  family:    will contain scored results
     *  query:     query sequence
     *  min_match: minimum number of results required
     *  max_match: maximum number of results desired
     *  min_score: minimum relative score
     *  max_score: maximum relative score
     *  arb:       pointer to matching arb database
     *  noid:      skip matches containing query
     *  min_len:   skip matches shorter
     *  num_full:  minimum "full length", ignoring score
     *  minlen_full: length to be considered "full"
     *  range_cover: minimum sequences touching alignment edge
     *  leave_query_out: drop sequence with matching id
     */
    virtual double match(std::vector<cseq> &family,
                         const cseq& query,
                         int min_match,
                         int max_match,
                         float min_score,
                         float max_score,
                         query_arb *arb,
                         bool noid,
                         int minlen,
                         int num_full,
                         int minlen_full,
                         int range_cover,
                         bool leave_query_out) = 0;

    virtual void find(const cseq& query, std::vector<cseq>& results, int max) = 0;

    virtual unsigned int size() const = 0;
};

} // namespace sina

#endif // _SEARCH_H_

/*
  Local Variables:
  mode:c++
  c-file-style:"stroustrup"
  c-file-offsets:((innamespace . 0)(inline-open . 0)(case-label . 0))
  indent-tabs-mode:nil
  fill-column:99
  End:
*/
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:encoding=utf-8:textwidth=99 :
