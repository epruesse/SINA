#include "../famfinder.h"

#define BOOST_TEST_MODULE famfinder

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
namespace bdata = boost::unit_test::data;

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <random>
#include <set>
#include <algorithm>
#include <initializer_list>

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
    const unsigned int N{20};
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


void configure(std::initializer_list<const char*> l) {
    const char* cmd[l.size()+1];
    cmd[0] = "sina";
    int i = 1;
    for (auto lv : l) {
        cmd[i++] = lv;
    }
    po::variables_map vm;
    po::options_description od, adv_od;
    famfinder::get_options_description(od, adv_od);
    od.add(adv_od);
    po::store(
        po::parse_command_line(l.size()+1, cmd, od),
        vm);
    try {
        //BOOST_TEST(vm["db"].as<fs::path>() == "test_data/ltp_reduced.arb");
        po::notify(vm);
        famfinder::validate_vm(vm, od);
    } catch (std::logic_error &e) {
        BOOST_TEST(e.what() == "");
    }
}

BOOST_AUTO_TEST_SUITE(famfinder_tests)


BOOST_DATA_TEST_CASE_F(Fixture, turn, bdata::make({"internal", "pt-server"}), engine) {
    configure({"--db", arbdb->getFileName().c_str(), "--fs-engine", engine});

    famfinder finder;

    for (unsigned int i = 0; i < N; i++) {
        cseq query = arbdb->getCseq(ids[i]);
        BOOST_TEST_CONTEXT("Test " << i << "(sequence " << query.getName() << ")") {
            BOOST_TEST_INFO("Unturned; just revcomp");
            BOOST_TEST(finder.turn_check(query, false) == 0);
            BOOST_TEST_INFO("Unturned; all orientations");
            BOOST_TEST(finder.turn_check(query, true) == 0);
            query.reverse();
            BOOST_TEST_INFO("Reversed; all orientations");
            BOOST_TEST(finder.turn_check(query, true) == 1);
            query.complement();
            BOOST_TEST_INFO("Revcomp'ed; just revcomp");
            BOOST_TEST(finder.turn_check(query, false) == 3);
            BOOST_TEST_INFO("Revcomp'ed; all orientations");
            BOOST_TEST(finder.turn_check(query, true) == 3);
            query.reverse();
            BOOST_TEST_INFO("Comp'ed; all orientations");
            BOOST_TEST(finder.turn_check(query, true) == 2);
        }
        std::cerr << ".";
    }
    std::cerr << std::endl;
}



BOOST_AUTO_TEST_SUITE_END(); // famfinder_tests

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
