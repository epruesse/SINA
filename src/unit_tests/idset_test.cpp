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

#include "../idset.h"

#define BOOST_TEST_MODULE bitmap
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
namespace bdata = boost::unit_test::data;

#include <random>
#include <set>

BOOST_TEST_SPECIALIZED_COLLECTION_COMPARE(std::vector<idset::inc_t>);

BOOST_AUTO_TEST_SUITE(bitmap_test);


struct test_set {
    test_set() = default;
    void init (unsigned int size, int fill, int seed) {
        // fill data with random sequence of monotonically rising numbers
        n = size * fill / 100;
        std::mt19937 engine(seed);
        std::uniform_real_distribution<> dist;
        std::set<unsigned int> d;
        for (int i=0; i<n ; i++) {
            unsigned int val = dist(engine) * (size-1);
            while (not d.insert(val).second)
            {
                val++;
                if (val >= size) {
                    val=0;
                }
            }
        }
        data.resize(n);
        std::copy(d.begin(), d.end(), data.begin());
        std::sort(data.begin(), data.end());

        // fill expected_count
        expected_counts.resize(size, 0);
        for (auto i : data) {
            expected_counts[i]++;
        }
    }
    std::vector<unsigned int> data;
    idset::inc_t expected_counts;
    int n;
};

#if 1
int map_sizes[] = {0, 255,256,257, 10000};
int map_fill[]  = {0, 10, 50, 100};
int map_seed[]  = {132456, 54321, 242424};
#else
int map_sizes[] = {10};
int map_fill[]  = {0, 50, 100};
int map_seed[]  = {132456};
#endif
idset* map_type[] =  {new bitmap(0), new imap_abs(0), new vlimap_abs(0), new vlimap(0) };


BOOST_DATA_TEST_CASE_F(test_set,
		       bitmap_test,
		       bdata::make(map_sizes) *
                       bdata::make(map_fill) *
                       bdata::make(map_seed),
		       map_size, map_fill, map_seed) {
    // init random set
    init(map_size, map_fill, map_seed);
  
    // init bitmap
    bitmap b(map_size);
    for (auto i : data) {
        b.set(i);
    }
  
    // check get()
    int matching_set = 0;
    for (auto i: data) {
        if (b.get(i)) {
            matching_set++;
        }
    }
    BOOST_CHECK_EQUAL(n, matching_set);
  
    // check count()
    BOOST_CHECK_EQUAL(n, b.count());

    // check increment()
    idset::inc_t count(map_size, 0);
    b.increment(count);
    matching_set = 0;
    for (auto i: data) {
        matching_set += count[i];
    }
    BOOST_CHECK_EQUAL(n, matching_set);
}


BOOST_DATA_TEST_CASE_F(test_set,
		       idset_test,
		       bdata::make(map_sizes) *
		       bdata::make(map_fill) *
		       bdata::make(map_seed) *
		       bdata::make(map_type),
		       map_size, map_fill, map_seed, type) {
    // init random set
    init(map_size, map_fill, map_seed);
  
    // init bitmap
    idset* b = type->make_new(map_size);
    BOOST_REQUIRE_EQUAL(b->size(), 0);
    for (auto i : data) {
        b->push_back(i);
    }
    BOOST_CHECK_EQUAL(data.size(), b->size());

    // check increment()
    idset::inc_t count(map_size, 0);
    b->increment(count);
    int matching_set = 0;
    for (auto i: data) {
        matching_set += count[i];
    }
    BOOST_CHECK_EQUAL(n, matching_set);
    if (n != matching_set) {
        for (auto i: data) {
            BOOST_CHECK_EQUAL(count[i]*i, i);
        }
    }
}


BOOST_DATA_TEST_CASE_F(test_set,
		       vlipmap_abs_test,
		       bdata::make(map_sizes) *
                       bdata::make(map_fill) *
                       bdata::make(map_seed),
		       map_size, map_fill, map_seed) {
    // init random set
    init(map_size, map_fill, map_seed);

    // init vlimap
    vlimap_abs a, b;

    int mid = data.size() / 2;
    for (int i = 0; i < mid; i++) {
        a.push_back(data[i]);
    }
    for (int i = mid; i < data.size(); i++) {
        b.push_back(data[i]);
    }

    idset::inc_t count(map_size, 0);
    a.increment(count);
    b.increment(count);
    BOOST_TEST(count == expected_counts, boost::test_tools::per_element()) ;
}


BOOST_DATA_TEST_CASE_F(test_set,
		       vlipmap_test,
		       bdata::make(map_sizes) *
                       bdata::make(map_fill) *
                       bdata::make(map_seed),
		       map_size, map_fill, map_seed) {
    // init random set
    init(map_size, map_fill, map_seed);

    // init vlimap
    vlimap a(map_size), b(map_size);

    int mid = data.size() / 2;
    for (int i = 0; i < mid; i++) {
        a.push_back(data[i]);
    }
    for (int i = mid; i < data.size(); i++) {
        b.push_back(data[i]);
    }

    {
        idset::inc_t count(map_size, 0);
        a.increment(count);
        b.increment(count);
        BOOST_TEST(count == expected_counts, boost::test_tools::per_element()) ;
    }

    BOOST_TEST_CONTEXT("split = " << mid <<
                       "(" << (data.size()?data[mid-1]:-1) <<
                       "/" << (data.size()?data[mid]:-1) << ")") {
        idset::inc_t count(map_size, 0);
        a.append(b);
        a.increment(count);
        BOOST_TEST(count == expected_counts, boost::test_tools::per_element()) ;
    }

    {
        idset::inc_t count(map_size, 1);
        a.invert();
        int res = a.increment(count);
        BOOST_CHECK_EQUAL(res, 1);
        BOOST_TEST(count == expected_counts, boost::test_tools::per_element()) ;
    }
    {
        idset::inc_t count(map_size, 1);
        vlimap *tmp = new vlimap(a);
        int res = tmp->increment(count);
        BOOST_CHECK_EQUAL(res, 1);
        BOOST_TEST(count == expected_counts, boost::test_tools::per_element()) ;
        delete tmp;
    }
}



BOOST_AUTO_TEST_SUITE_END(); // cseq_test

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
