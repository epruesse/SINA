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

#ifndef _QUERY_PT_H_
#define _QUERY_PT_H_

#include "search.h"

#include <exception>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace sina {

class query_pt_exception : public std::exception {
    std::string message;
public:
    query_pt_exception(std::string  msg) noexcept;
    ~query_pt_exception() noexcept override;
    const char* what() const noexcept override;
};

class query_pt_pool : public search {
    struct pimpl;
    std::shared_ptr<pimpl> impl;
public:
    static query_pt_pool* get_pool(boost::filesystem::path filename,
                            int k=10, bool fast=true, bool norel=false, int mk=0,
                            std::string portname="");
    query_pt_pool(std::shared_ptr<pimpl>);
    ~query_pt_pool() override;
private:
    query_pt_pool() = delete;
    query_pt_pool(const query_pt_pool&) = delete;


    void find(const cseq& query, result_vector& results, unsigned int max) override;

    double match(result_vector &family,
                 const cseq& queryc,
                 int min_match,
                 int max_match,
                 float min_score,
                 float max_score,
                 query_arb *arb,
                 bool noid,
                 int min_len,
                 int num_full,
                 int full_min_len,
                 int range_cover,
                 bool leave_query_out) override;

    unsigned int size() const override;
};

class query_pt : public search {
public:
    static query_pt* get_pt_search(const boost::filesystem::path& filename,
                                   int k=10,
                                   bool fast=true,
                                   bool norel=false,
                                   int mk=0,
                                   std::string portname="");

    query_pt(const char* portname,
             const char* dbname,
             bool fast=true,
             int k=10,
             int mk=0,
             bool norel=false);
    ~query_pt() override;

    void find(const cseq& query, result_vector& results, unsigned int max) override;
    unsigned int size() const override;

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
    double match(result_vector &family,
                 const cseq& queryc,
                 int min_match,
                 int max_match,
                 float min_score,
                 float max_score,
                 query_arb *arb,
                 bool noid,
                 int min_len,
                 int num_full,
                 int full_min_len,
                 int range_cover,
                 bool leave_query_out) override;

    void set_find_type_fast(bool fast);
    void set_probe_len(int len);
    void set_mismatches(int len);
    void set_sort_type(bool absolute);
    void set_range(int startpos=-1, int stoppos=-1);
    void unset_range();

    static void get_options_description(boost::program_options::options_description& all,
                                        boost::program_options::options_description& adv);
    static void validate_vm(boost::program_options::variables_map& /*unused*/,
                            boost::program_options::options_description& /*unused*/);

private:
    struct priv_data;
    priv_data *data;
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
