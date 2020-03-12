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

#include "../query_arb.h"
#include "../tempfile.h"
#include "../log.h"

#define BOOST_TEST_MODULE query_arb_test

#include <boost/test/unit_test.hpp>
namespace utf = boost::unit_test;
#include <boost/test/data/test_case.hpp>
namespace bdata = utf::data;
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

#include <spdlog/sinks/ostream_sink.h>

using namespace sina;


struct GlobalFixture {
    static GlobalFixture*& instance() {
        static GlobalFixture* inst{nullptr};
        return inst;
    }

    GlobalFixture() {
        instance() = this;

        int argc = boost::unit_test::framework::master_test_suite().argc;
        if (argc != 2) {
            throw std::runtime_error("Need exactly 1 argument: the reference arb file");
        }

        char** argv = boost::unit_test::framework::master_test_suite().argv;
        query_arb* ref = query_arb::getARBDB(argv[1]);
        query_arb* tmp = query_arb::getARBDB(small_arb_path());
        std::vector<std::string> ids = ref->getSequenceNames();

        std::srand(1234);
        std::random_shuffle(ids.begin(), ids.end());
        for (size_t i = 0; i < n_seq; ++i) {
            cseq c = ref->getCseq(ids[i]);
            tmp->putCseq(c);
        }
        tmp->save();
        query_arb::closeOpenARBDBs();
    }
    ~GlobalFixture() {
        instance() = nullptr;
    }

    static TempFile& small_arb_path() {
        return (*instance())._small_arb_path;
    }

    TempFile _small_arb_path;
    const size_t n_seq{50};
};

struct CaseFixture {
    CaseFixture() {}
    ~CaseFixture() {
        query_arb::closeOpenARBDBs();
    }

    query_arb* tmparb() {
        if (_tmparb == nullptr) {
            _tmparb = query_arb::getARBDB(tmpfile);
        }
        return _tmparb;
    }
    query_arb* smallarb() {
        if (_smallarb == nullptr) {
            _smallarb = query_arb::getARBDB(GlobalFixture::small_arb_path());
        }
        return _smallarb;
    }
    TempFile tmpfile;
    query_arb *_tmparb{nullptr}, *_smallarb{nullptr};
};


struct what_starts_with {
    explicit what_starts_with(const std::string& prefix) : _prefix(prefix) {}
    bool operator()(const std::exception& e) const {
        return std::string(e.what()).rfind(_prefix, 0) == 0;
    }
    std::string _prefix;
};

BOOST_GLOBAL_FIXTURE(GlobalFixture);

BOOST_FIXTURE_TEST_SUITE(query_arb_test, CaseFixture);

BOOST_AUTO_TEST_CASE(new_file) {
    query_arb* arb = tmparb();
    arb->save();
    fs::exists(tmpfile);
    query_arb::closeOpenARBDBs();
    fs::remove(tmpfile);
    arb = query_arb::getARBDB(tmpfile);
    arb->save();
    fs::exists(tmpfile);
}


BOOST_AUTO_TEST_CASE(cache) {
    query_arb *arb = smallarb();
    arb->loadCache();
    std::vector<cseq*> cache = arb->getCacheContents();
    BOOST_CHECK_EQUAL(cache.size(), GlobalFixture::instance()->n_seq);
}

BOOST_AUTO_TEST_CASE(empty_filename) {
    BOOST_CHECK_EXCEPTION(
        query_arb::getARBDB(""),
        query_arb_exception,
        what_starts_with("Empty ARB database")
        );
}

BOOST_AUTO_TEST_CASE(unable_to_open) {
    TempFile notAnArbFile;
    {
        boost::filesystem::ofstream out(notAnArbFile);
        out << "Nothing" << std::endl;
    }

    BOOST_CHECK_EXCEPTION(
        query_arb::getARBDB(notAnArbFile),
        query_arb_exception,
        what_starts_with("Unable to open")
        );
}

BOOST_AUTO_TEST_CASE(unable_to_save) {
    query_arb *arb = tmparb();
    TempFile write_protected;
    {
        boost::filesystem::ofstream out(write_protected);
    }
    boost::filesystem::permissions(
        write_protected,
        boost::filesystem::owner_read);

    std::ostringstream oss;
    auto ostream_sink = std::make_shared<spdlog::sinks::ostream_sink_mt> (oss);
    ostream_sink->set_level(spdlog::level::err);
    Log::add_sink(ostream_sink);
    arb->saveAs(write_protected);
    Log::remove_sink(ostream_sink);

    BOOST_CHECK(oss.str().find("Error while trying to save") != std::string::npos);
}



BOOST_AUTO_TEST_SUITE_END(); // rw_csv_test

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
