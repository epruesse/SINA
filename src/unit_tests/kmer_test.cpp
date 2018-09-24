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

#include "../kmer.h"

#define BOOST_TEST_MODULE kmer_search

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
namespace bdata = boost::unit_test::data;

#include <string>

BOOST_AUTO_TEST_SUITE(kmer_test);

using namespace sina;

BOOST_AUTO_TEST_CASE(kmer_generator_test_bounds) {
    BOOST_CHECK_THROW(kmer_generator kmer(0), std::runtime_error);
    BOOST_CHECK_THROW(kmer_generator kmer(17), std::runtime_error);
}


// data for BOOST_DATA_TEST_CASE(kmer_generator_test...):

const char test_sequence[] =
    "AGCTN"
    "AGCTAGCTN"
    "AGCTAGCTAGCTN";

const int test_sequence_len = std::extent<decltype(test_sequence)>::value - 1;

const int try_k[] = {1,2,3,4,8};
const int try_k_len = std::extent<decltype(try_k)>::value;

const bool valid_k[try_k_len][test_sequence_len] = {
    { // k = 1
        1, 1, 1, 1, 0,
        1, 1, 1, 1,  1, 1, 1, 1, 0,
        1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1, 0
    },
    { // k = 2
        0, 1, 1, 1, 0,
        0, 1, 1, 1,  1, 1, 1, 1, 0,
        0, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1, 0
    },
    { // k = 3
        0, 0, 1, 1, 0,
        0, 0, 1, 1,  1, 1, 1, 1, 0,
        0, 0, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1, 0
    },
    { // k = 4
        0, 0, 0, 1, 0,
        0, 0, 0, 1,  1, 1, 1, 1, 0,
        0, 0, 0, 1,  1, 1, 1, 1,  1, 1, 1, 1, 0
    },
    { // k = 8
        0, 0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 1, 0,
        0, 0, 0, 0,  0, 0, 0, 1,  1, 1, 1, 1, 0
    }
};

const bool first_k[try_k_len][test_sequence_len] = {
    { // k = 1
        1, 1, 1, 1, 0,
        0, 0, 0, 0,  0, 0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0
    },
    { // k = 2
        0, 1, 1, 1, 0,
        0, 0, 0, 0,  1, 0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0
    },
    { // k = 3
        0, 0, 1, 1, 0,
        0, 0, 0, 0,  1, 1, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0
    },
    { // k = 4
        0, 0, 0, 1, 0,
        0, 0, 0, 0,  1, 1, 1, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0, 0
    },
    { // k = 8
        0, 0, 0, 0, 0,
        0, 0, 0, 0,  0, 0, 0, 1, 0,
        0, 0, 0, 0,  0, 0, 0, 0,  1, 1, 1, 0, 0
    }
};

const unsigned int A = BASE_A;
const unsigned int G = BASE_G;
const unsigned int C = BASE_C;
const unsigned int T = BASE_TU;

const unsigned int kmers_k[try_k_len][test_sequence_len] = {
    { // k = 1
        A, G, C, T, 0,
        A, G, C, T,  A, G, C, T, 0,
        A, G, C, T,  A, G, C, T,  A, G, C, T, 0,
    },
    { // k = 2
        0, A<<2|G, G<<2|C, C<<2|T, 0,
        0, A<<2|G, G<<2|C, C<<2|T,
        T<<2|A, A<<2|G, G<<2|C, C<<2|T, 0,
        0, A<<2|G, G<<2|C, C<<2|T,
        T<<2|A, A<<2|G, G<<2|C, C<<2|T,
        T<<2|A, A<<2|G, G<<2|C, C<<2|T, 0
    },
    { // k = 3
        0, 0, A<<4|G<<2|C, G<<4|C<<2|T, 0,
        0, 0, A<<4|G<<2|C, G<<4|C<<2|T,
        C<<4|T<<2|A, T<<4|A<<2|G, A<<4|G<<2|C, G<<4|C<<2|T, 0,
        0, 0, A<<4|G<<2|C, G<<4|C<<2|T,
        C<<4|T<<2|A, T<<4|A<<2|G, A<<4|G<<2|C, G<<4|C<<2|T,
        C<<4|T<<2|A, T<<4|A<<2|G, A<<4|G<<2|C, G<<4|C<<2|T, 0
    },
    { // k = 4
        0, 0, 0, A<<6|G<<4|C<<2|T, 0,
        0, 0, 0, A<<6|G<<4|C<<2|T,
        G<<6|C<<4|T<<2|A, C<<6|T<<4|A<<2|G, T<<6|A<<4|G<<2|C, A<<6|G<<4|C<<2|T, 0,
        0, 0, 0, A<<6|G<<4|C<<2|T,
        G<<6|C<<4|T<<2|A, C<<6|T<<4|A<<2|G, T<<6|A<<4|G<<2|C, A<<6|G<<4|C<<2|T,
        G<<6|C<<4|T<<2|A, C<<6|T<<4|A<<2|G, T<<6|A<<4|G<<2|C, A<<6|G<<4|C<<2|T, 0
    },
    { // k = 8
        0, 0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, A<<14|G<<12|C<<10|T<<8|A<<6|G<<4|C<<2|T, 0,
        0, 0, 0, 0,
        0, 0, 0, A<<14|G<<12|C<<10|T<<8|A<<6|G<<4|C<<2|T,
        G<<14|C<<12|T<<10|A<<8|G<<6|C<<4|T<<2|A,
        C<<14|T<<12|A<<10|G<<8|C<<6|T<<4|A<<2|G,
        T<<14|A<<12|G<<10|C<<8|T<<6|A<<4|G<<2|C,
        A<<14|G<<12|C<<10|T<<8|A<<6|G<<4|C<<2|T, 0
    }
};

BOOST_DATA_TEST_CASE(kmer_generator_test,
		     bdata::xrange(try_k_len) ^ bdata::make(try_k),
		     n, k
    ) {
    const auto& valid = valid_k[n];
    const auto& kmers = kmers_k[n];

    kmer_generator kmer(k);

    for (int i = 0; i < test_sequence_len; ++i) {
        kmer.push(test_sequence[i]);
        BOOST_CHECK_EQUAL(kmer.good(), valid[i]);
        BOOST_CHECK_MESSAGE(kmer.good() == valid[i], ""
                            << " i=" << i
                            << " kmer=" << kmer
                            << " b=" << test_sequence[i]
                            << " gc=" << kmer.get_good_count());
        if (kmer.good()) {
            BOOST_CHECK_EQUAL(kmer, kmers[i]);
            BOOST_CHECK_MESSAGE(kmer == kmers[i], ""
                                << " i=" << i
                                << " b=" << test_sequence[i]
                                << " kmer=" << kmer
                                << " test_kmer=" << kmer_generator(k, kmers[i])
                );
        }
    }
}

BOOST_DATA_TEST_CASE(kmer_iterable_generator_test,
		     bdata::xrange(try_k_len) ^ bdata::make(try_k),
		     n, k
    ) {
    const auto& valid = valid_k[n];
    const auto& kmers = kmers_k[n];

    int i = 0;
    std::string test_string(test_sequence);
    for (const auto& kmer : all_kmers(test_string, k)) {
        while (not valid[i]) { ++i; }
        BOOST_CHECK_EQUAL(kmer, kmers[i]);
        BOOST_CHECK_MESSAGE(kmer == kmers[i], ""
                            << " i=" << i
                            << " b=" << test_sequence[i]
                            << " kmer=" << kmer
                            << " test_kmer=" << kmer_generator(k, kmers[i])
            );
        ++i;
    }
}

BOOST_DATA_TEST_CASE(kmer_unique_generator_test,
		     bdata::xrange(try_k_len) ^ bdata::make(try_k),
		     n, k
    ) {
    const auto& valid = first_k[n];
    const auto& kmers = kmers_k[n];

    std::unordered_set<unsigned int> seen;

    unique_kmer_generator kmer(seen, k);

    for (int i = 0; i < test_sequence_len; ++i) {
        kmer.push(test_sequence[i]);
        BOOST_CHECK_EQUAL(kmer.good(), valid[i]);
        BOOST_CHECK_MESSAGE(kmer.good() == valid[i], ""
                            << " i=" << i
                            << " kmer=" << kmer
                            << " b=" << test_sequence[i]
                            << " gc=" << kmer.get_good_count());
        if (kmer.good()) {
            BOOST_CHECK_EQUAL(kmer, kmers[i]);
            BOOST_CHECK_MESSAGE(kmer == kmers[i], ""
                                << " i=" << i
                                << " b=" << test_sequence[i]
                                << " kmer=" << kmer
                                << " test_kmer=" << kmer_generator(k, kmers[i])
                );
        }
    }
}

BOOST_DATA_TEST_CASE(kmer_iterable_unique_generator_test,
		     bdata::xrange(try_k_len) ^ bdata::make(try_k),
		     n, k
    ) {
    const auto& valid = first_k[n];
    const auto& kmers = kmers_k[n];

    int i = 0;
    std::string test_string(test_sequence);
    std::unordered_set<unsigned int> seen;
    for (const auto& kmer : unique_kmers(test_string, seen, k)) {
        while (not valid[i]) { ++i; }
        BOOST_CHECK_EQUAL(kmer, kmers[i]);
        BOOST_CHECK_MESSAGE(kmer == kmers[i], ""
                            << " i=" << i
                            << " b=" << test_sequence[i]
                            << " kmer=" << kmer
                            << " test_kmer=" << kmer_generator(k, kmers[i])
            );
        ++i;
    }
}


// data for BOOST_DATA_TEST_CASE(prefix_kmer_generator_test,...):

const int kmer_prefix_n = 2;
const unsigned int kmer_prefix_lens[try_k_len][kmer_prefix_n] = {
    { 1, 1 }, // k=1
    { 1, 2 }, // k=2
    { 2, 3 }, // k=3
    { 2, 4 }, // k=4
    { 2, 4 }, // k=8
};

const unsigned int kmer_prefixes[try_k_len][kmer_prefix_n] = {
    { A, T }, // k=1
    { A, G<<2|C }, // k=2
    { G<<2|C, T<<4|A<<2|G }, // k=3
    { C<<2|T, G<<6|C<<4|T<<2|A }, // k=4
    { T<<2|A, A<<6|G<<4|C<<2|T } // k=8
};

const unsigned int kmer_prefix_counts[try_k_len][kmer_prefix_n] = {
    { 6, 6 }, // k=1
    { 6, 6 }, // k=2
    { 6, 3 }, // k=3
    { 3, 3 }, // k=4
    { 1, 3 }, // k=8
};

const unsigned int kmer_unique_prefix_counts[try_k_len][kmer_prefix_n] = {
    { 1, 1 }, // k=1
    { 1, 1 }, // k=2
    { 1, 1 }, // k=3
    { 1, 1 }, // k=4
    { 1, 1 }, // k=8
};

BOOST_DATA_TEST_CASE(prefix_kmer_generator_test,
		     (bdata::xrange(try_k_len) ^ bdata::make(try_k)) *
		     bdata::xrange(kmer_prefix_n),
		     n, k, p
    ) {
    const auto& valid = valid_k[n];
    const auto& kmers = kmers_k[n];
    const auto& prefix_len = kmer_prefix_lens[n][p];
    const auto& prefix = kmer_prefixes[n][p];
    const auto& prefix_count = kmer_prefix_counts[n][p];

    prefix_kmer_generator kmer(k, prefix_len, prefix);

    int count = 0;
    for (int i = 0; i < test_sequence_len; ++i) {
        kmer.push(test_sequence[i]);
        unsigned int current_prefix = kmers[i] >> (k-prefix_len)*2;
        bool should_be_valid = valid[i] && current_prefix == prefix;
        BOOST_CHECK_EQUAL(kmer.good(), should_be_valid);
        BOOST_CHECK_MESSAGE(kmer.good() == should_be_valid, ""
                            << " i=" << i
                            << " kmer=" << kmer
                            << " b=" << test_sequence[i]
                            << " prefix_len=" << prefix_len
                            << " prefix=" << kmer_generator(prefix_len,prefix)
                            << " cur_prefix=" << kmer_generator(prefix_len,current_prefix)

                            << " gc=" << kmer.get_good_count());
        if (kmer.good()) {
            BOOST_CHECK_EQUAL(kmer, kmers[i]);
            BOOST_CHECK_MESSAGE(kmer == kmers[i], ""
                                << " i=" << i
                                << " b=" << test_sequence[i]
                                << " kmer=" << kmer
                                << " test_kmer=" << kmer_generator(k, kmers[i])
                );
            ++count;
        }
    }
    BOOST_CHECK_EQUAL(count, prefix_count);
}


BOOST_DATA_TEST_CASE(kmer_iterable_prefix_generator_test,
		     (bdata::xrange(try_k_len) ^ bdata::make(try_k)) *
		     bdata::xrange(kmer_prefix_n),
		     n, k, p
    ) {
    const auto& valid = valid_k[n];
    const auto& kmers = kmers_k[n];
    const auto& prefix_len = kmer_prefix_lens[n][p];
    const auto& prefix = kmer_prefixes[n][p];
    const auto& prefix_count = kmer_prefix_counts[n][p];

    int i = 0;
    int count = 0;
    std::string test_string(test_sequence);
    for (const auto& kmer : prefix_kmers(test_string, k, prefix_len, prefix)) {
        bool should_be_valid = false;
        unsigned int current_prefix;
        do {
            current_prefix = kmers[i] >> (k-prefix_len)*2;
            should_be_valid = valid[i] && current_prefix == prefix;
        } while (should_be_valid == false && ++i);
        BOOST_CHECK_EQUAL(kmer, kmers[i]);
        BOOST_CHECK_MESSAGE(kmer == kmers[i], ""
                            << " i=" << i
                            << " b=" << test_sequence[i]
                            << " kmer=" << kmer
                            << " test_kmer=" << kmer_generator(k, kmers[i])
                            << " b=" << test_sequence[i]
                            << " prefix_len=" << prefix_len
                            << " prefix=" << kmer_generator(prefix_len,prefix)
                            << " cur_prefix=" << kmer_generator(prefix_len,current_prefix)
            );
        ++i;
        ++count;
    }
    BOOST_CHECK_EQUAL(count, prefix_count);
}


BOOST_DATA_TEST_CASE(kmer_iterable_unique_prefix_generator_test,
		     (bdata::xrange(try_k_len) ^ bdata::make(try_k)) *
		     bdata::xrange(kmer_prefix_n),
		     n, k, p
    ) {
    const auto& valid = first_k[n];
    const auto& kmers = kmers_k[n];
    const auto& prefix_len = kmer_prefix_lens[n][p];
    const auto& prefix = kmer_prefixes[n][p];
    const auto& prefix_count = kmer_unique_prefix_counts[n][p];

    int i = 0;
    int count = 0;
    std::unordered_set<unsigned int> seen;

    std::string test_string(test_sequence);
    for (const auto& kmer : unique_prefix_kmers(test_string, seen, k, prefix_len, prefix)) {
        bool should_be_valid = false;
        unsigned int current_prefix;
        do {
            current_prefix = kmers[i] >> (k-prefix_len)*2;
            should_be_valid = valid[i] && current_prefix == prefix;
        } while (should_be_valid == false && ++i);
        BOOST_CHECK_EQUAL(kmer, kmers[i]);
        BOOST_CHECK_MESSAGE(kmer == kmers[i], ""
                            << " i=" << i
                            << " b=" << test_sequence[i]
                            << " kmer=" << kmer
                            << " test_kmer=" << kmer_generator(k, kmers[i])
                            << " b=" << test_sequence[i]
                            << " prefix_len=" << prefix_len
                            << " prefix=" << kmer_generator(prefix_len,prefix)
                            << " cur_prefix=" << kmer_generator(prefix_len,current_prefix)
            );
        ++i;
        ++count;
    }
    BOOST_CHECK_EQUAL(count, prefix_count);
}

BOOST_AUTO_TEST_SUITE_END(); // kmer_test


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
