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

#include "../kmer_search.h"
#include "../query_pt.h"

#define BOOST_TEST_MODULE kmer_search

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

#include <random>
#include <set>
#include <algorithm>

#include "../query_arb.h"


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

struct Fixture {
    query_arb *arbdb;
    std::vector<std::string> ids;
    const unsigned int N{1000};
    const unsigned int M{50};

    Fixture() {
        std::srand(1234);
        int argc = boost::unit_test::framework::master_test_suite().argc;
        char** argv = boost::unit_test::framework::master_test_suite().argv;
        BOOST_REQUIRE(argc>1);
        fs::path database = argv[1];
        BOOST_TEST_MESSAGE("Loading test database " <<  database << " for testing");
        arbdb = query_arb::getARBDB(database);
        ids = arbdb->getSequenceNames();
        BOOST_REQUIRE(ids.size() > N);
        shuffle_n(ids.begin(), ids.end(), N);
    }

    ~Fixture() {
        BOOST_TEST_MESSAGE("Destroying fixture");
    }
};

BOOST_AUTO_TEST_SUITE(search_tests);

BOOST_FIXTURE_TEST_CASE(kmer_simple1, Fixture) {
    fs::path dbname = arbdb->getFileName();
    fs::path idxname = fs::path(dbname).replace_extension("sidx");
    if (fs::exists(idxname)) {
        fs::remove(idxname);
    }
    // test runs twice, once with index cache absent (removed above if present),
    // and once with index cache generated during first run.
    for (int run = 0; run < 2; ++run) {
        kmer_search::release_kmer_search(dbname);
        kmer_search *search_index = kmer_search::get_kmer_search(dbname);
        search::result_vector family;
        for (unsigned int i = 0; i < N; i++) {
            if (i % (N/50) == 0) {
                std::cerr << ".";
            }
            cseq query = arbdb->getCseq(ids[i]);
            search_index->find(query, family, M);
            BOOST_TEST(family.size() == M);
            float max_score = family[0].score;
            auto self = std::find_if(family.begin(), family.end(),
                                     [&](const search::result_item &i) {
                                         return i.sequence->getName() == query.getName();}
                );
            BOOST_REQUIRE((self != family.end()));
            if (self != family.end()) {
                BOOST_TEST(self->score == max_score);
            }
        }
        std::cerr << std::endl;
        delete search_index;
    }
}


BOOST_FIXTURE_TEST_CASE(pt_simple, Fixture, *boost::unit_test::tolerance(0.0001)) {
    search *search_index = query_pt::get_pt_search(arbdb->getFileName());
    search::result_vector family;
    for (unsigned int i = 0; i < N; i++) {
        if (i % (N/50) == 0) {
            std::cerr << ".";
        }
        cseq query = arbdb->getCseq(ids[i]);
        search_index->find(query, family, M);
        BOOST_TEST(family.size() == M);
        float max_score = family[0].score;
        auto self = std::find_if(family.begin(), family.end(),
                                 [&](const search::result_item &i) {
                                return i.sequence->getName() == query.getName();}
            );
        BOOST_TEST((self != family.end()));
        // PT server counts duplicate kmers twice, allow for some discrepancy
        BOOST_TEST(self->score > max_score - 5); // FIXME: now relative?
    }
    std::cerr << std::endl;
    delete search_index;
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
