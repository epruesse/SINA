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

#ifndef _QUERY_PT_H_
#define _QUERY_PT_H_

#include "search.h"

#include <exception>

#include <boost/program_options.hpp>

namespace sina {

class query_pt_exception : public std::exception {
    std::string message;
public:
    query_pt_exception(std::string _message) throw();
    ~query_pt_exception() throw();
    virtual const char* what() const throw();
};


class query_pt : public search {
private:
    void init();
    void exit();
    void restart();
public:
    query_pt(const char* portname,
             const char* dbname,
             bool fast=true,
             int k=10,
             int mk=0,
             bool norel=false);
    ~query_pt();

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
                         bool leave_query_out) override;

    virtual double match(std::vector<cseq> &family,
                         const cseq& sequence,
                         int min_match,
                         int max_match,
                         float min_score) {
        return match(family, sequence, min_match, max_match, min_score, 2.0,
                     NULL, false, 0, 0, 0, 0, false);
    };

    int turn_check(const cseq& query, bool all);

    void set_find_type_fast(bool fast);
    void set_probe_len(int len);
    void set_mismatches(int len);
    void set_sort_type(bool absolute);
    void set_range(int startpos, int stoppos);
    void unset_range();

    static void get_options_description(boost::program_options::options_description& all,
                                        boost::program_options::options_description& adv);
    static void validate_vm(boost::program_options::variables_map&,
                            boost::program_options::options_description&);

private:
    struct priv_data;
    priv_data &data;
    struct options;
    static struct options *opts;
};

} // namespace sina

#endif // _QUERY_PT_H_

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
