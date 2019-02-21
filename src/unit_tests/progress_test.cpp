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


# include "../progress.h"

#define BOOST_TEST_MODULE progress

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
namespace utf = boost::unit_test;
namespace bdata = utf::data;

#include <string>
#include <iostream>
#include "spdlog/spdlog.h"

BOOST_AUTO_TEST_SUITE(progress_test);

using namespace sina;


BOOST_TEST_DECORATOR(* utf::timeout(1))
BOOST_DATA_TEST_CASE(
    format_bar_to,
    bdata::make({0, 1, 2, 3, 100, 10000}) * // totals
    bdata::make({0u, 1u, 2u, 3u, 10u, 100u}) *    // widths
    bdata::make({0., 0.001, 0.01, .1, .5, .51, .9, .99, .999, 1.}), // fills
    total, width, fill)
{
    int expected_full = fill * width + 1./12;
    int expected_remain = (fill * width - expected_full) * 12;

    FILE *devnull = fopen("/dev/null", "w");
    Progress p("description", total, /*ascii = */true, /*file = */devnull);
    
    fmt::memory_buffer buf;
    p.format_bar_to(buf, width, fill);
    BOOST_CHECK_EQUAL(buf.size(), width);
    BOOST_CHECK_EQUAL(std::count(buf.begin(), buf.end(), '#'), expected_full);
    BOOST_CHECK_EQUAL(std::count(buf.begin(), buf.end(), char('0' + expected_remain - 1))
                      , expected_remain?1:0);
    fclose(devnull);
}

BOOST_DATA_TEST_CASE(
    show_progress,
    bdata::make({60u,43u,42u}) *
    bdata::make({10,50,99}),
    width, n)
{
    fmt::memory_buffer buf;
    buf.resize(500);
    FILE * filebuf = fmemopen(buf.begin(), buf.size(), "rw");

    int total = 100;
    std::fill_n(buf.begin(), buf.size(), 0);
    Progress p("text", total, /*ascii=*/true, filebuf, width);
    p.update(n);
    fseek(filebuf, 0, SEEK_SET);
    p.show_progress();
    fflush(filebuf);
    std::string result(buf.begin());
    BOOST_TEST_INFO("result: '" << result << "'");
    std::string beg = fmt::format("text:  {}% |", n);
    std::string end = fmt::format("| {}/100 [00:00:00 / 00:00:00]\n", n);
    BOOST_TEST_INFO("expected begin: '" << beg << "'");
    BOOST_TEST_INFO("expected end: '" << end << "'");
    BOOST_REQUIRE_LE(result.size(), width+1);
    BOOST_REQUIRE_GE(result.size(), beg.size() + end.size());
    BOOST_CHECK_EQUAL(result.substr(0, beg.size()), beg);
    BOOST_CHECK_EQUAL(result.substr(result.size()-end.size(), end.size()), end);

    fclose(filebuf);
}

BOOST_AUTO_TEST_SUITE_END(); // progress_test


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
