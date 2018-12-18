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

#ifndef _KMER_SEARCH_H_
#define _KMER_SEARCH_H_

#include <map>
#include <boost/filesystem.hpp>

#include "search.h"

namespace sina {

class kmer_search : public search {
public:
    class result_iterator;

    static kmer_search* get_kmer_search(const boost::filesystem::path& filename, int k=10);

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
    double match(std::vector<cseq> &results,
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
                         bool leave_query_out) override;

    double match(std::vector<cseq> &family,
                         const cseq& sequence,
                         int min_match,
                         int max_match,
                         float min_score) override {
        return match(family, sequence, min_match, max_match, min_score, 2.0,
                     nullptr, false, 0, 0, 0, 0, false);
    };
    
    void build_index();
    void init();
    void find(const cseq& query, std::vector<cseq>& results, int max) override;

    class index;
private:
    kmer_search(query_arb* arbdb, int k=8);
    ~kmer_search() override;


    index &data;
    static std::map<std::string, kmer_search*> indices;
    static void destroy_indices();
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
