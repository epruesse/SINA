/*
  Copyright (c) 2006-2017 Elmar Pruesse <elmar.pruesse@ucdenver.edu>

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

#include "../kmer_search.h"

#define BOOST_TEST_MODULE kmer_search

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
namespace bdata = boost::unit_test::data;

#include <random>
#include <set>
#include <algorithm>

#include "../query_arb.h"

BOOST_AUTO_TEST_SUITE(kmer_search_test);

const char* DATABASE = "test_data/ltp_reduced.arb";
const int N = 100;
const int M = 10;

using namespace sina;

/** shuffles only the first n items **/
template<class ITER>
ITER shuffle_n(ITER begin, ITER end, size_t n) {
  size_t items = std::distance(begin, end);
  while (n--) {
    ITER it = begin;
    int ran = std::rand() % items;
    std::advance(it, ran);
    std::swap(*begin, *it);
    --items, ++begin;
  }
  return begin;
}

BOOST_AUTO_TEST_CASE(simple, *boost::unit_test::tolerance(0.0001)) {
  kmer_search *search_index = kmer_search::get_kmer_search(DATABASE);

  query_arb *arbdb = query_arb::getARBDB(DATABASE);
  std::vector<std::string> ids = arbdb->getSequenceNames();
  BOOST_TEST(ids.size() > N);
  shuffle_n(ids.begin(), ids.end(), N);

  std::vector<cseq> family;
  for (int i=0; i<N; i++) {
    cseq query = arbdb->getCseq(ids[i]);
    search_index->find(query, family, M);
    float max_score = family[0].getScore();
    std::vector<cseq>::iterator self;
    self = std::find_if(family.begin(), family.end(),
			[&](const cseq &c) {
			  return c.getName() == query.getName();}
			);
    BOOST_TEST((self != family.end()));
    BOOST_TEST(self->getScore() == max_score);
  }
}

BOOST_AUTO_TEST_SUITE_END(); // kmer_search_test

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
