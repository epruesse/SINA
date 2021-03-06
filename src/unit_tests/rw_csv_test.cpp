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

#include "../rw_csv.h"
#include "../tempfile.h"

#define BOOST_TEST_MODULE rw

#include <boost/test/unit_test.hpp>
namespace utf = boost::unit_test;
#include <boost/test/data/test_case.hpp>
namespace bdata = utf::data;
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#include <boost/program_options.hpp>
namespace po = boost::program_options;


using namespace sina;


struct Fixture {
    fs::path model{"sina-%%%%-%%%%.csv"};
    std::unique_ptr<TempFile> outfile;
    rw_csv::writer* writer{nullptr};
    int copy_relatives{0};
    std::vector<std::string> fields;

    std::string _output;
    int _nlines{0};

    Fixture()
        : outfile(new TempFile(model))
    {
        configure({});
    }

    void set_ext(const std::string& ext) {
        outfile = std::unique_ptr<TempFile>(
            new TempFile(fs::path("sina-%%%%-%%%%.") += ext)
            );
    }

    void write(const cseq& c) {
        if (!writer) {
            writer = new rw_csv::writer(*outfile, copy_relatives, fields);
        }
        tray t;
        t.aligned_sequence = const_cast<cseq*>(&c);
        (*writer)(t);
    }

    std::string& output() {
        if (writer) {
            delete writer;
            _output = outfile->load();
            writer = nullptr;
        }
        return _output;
    }

    int nlines() {
        auto txt = output();
        return std::count(txt.begin(), txt.end(), '\n');
    }

    void configure(std::initializer_list<const char*> l) {
        const char* cmd[l.size()+1];
        cmd[0] = "sina";
        int i = 1;
        for (auto lv : l) {
            cmd[i++] = lv;
        }
        po::variables_map vm;
        po::options_description od, adv_od;
        rw_csv::get_options_description(od, adv_od);
        od.add(adv_od);
        po::store(
            po::parse_command_line(l.size()+1, cmd, od),
            vm);
        try {
            po::notify(vm);
            rw_csv::validate_vm(vm, od);
        } catch (std::logic_error &e) {
            BOOST_TEST(e.what() == "");
        }
    }
};

BOOST_FIXTURE_TEST_SUITE(rw_csv_test, Fixture);

BOOST_AUTO_TEST_CASE(one_seq_no_data) {
    write(cseq("test_sequence"));
    BOOST_CHECK_EQUAL(nlines(), 2);
    BOOST_CHECK_EQUAL(
        output(),
        "name\n"
        "test_sequence\n"
        );
}

BOOST_AUTO_TEST_CASE(two_seq_no_data) {
    write(cseq("test_sequence1"));
    write(cseq("test_sequence2"));
    BOOST_CHECK_EQUAL(nlines(), 3);
    BOOST_CHECK_EQUAL(
        output(),
        "name\n"
        "test_sequence1\n"
        "test_sequence2\n");
}

BOOST_AUTO_TEST_CASE(one_seq_one_string) {
    cseq c("test_sequence");
    c.set_attr("col1", "Some test data");
    write(c);
    BOOST_CHECK_EQUAL(nlines(), 2);
    BOOST_CHECK_EQUAL(
        output(),
        "name,col1\n"
        "test_sequence,Some test data\n"
        );
}

BOOST_AUTO_TEST_CASE(one_seq_one_int) {
    cseq c("test_sequence");
    c.set_attr("col1", 123);
    write(c);
    BOOST_CHECK_EQUAL(nlines(), 2);
    BOOST_CHECK_EQUAL(
        output(),
        "name,col1\n"
        "test_sequence,123\n"
        );
}

BOOST_AUTO_TEST_CASE(escape_seqname) {
    cseq c("test_se,quence");
    write(c);
    BOOST_CHECK_EQUAL(
        output(),
        "name\n"
        "\"test_se,quence\"\n"
        );
}

BOOST_AUTO_TEST_CASE(escape_comma) {
    cseq c("seq");
    c.set_attr("data", "d1,d2,d2");
    write(c);
    BOOST_CHECK_EQUAL(
        output(),
        "name,data\n"
        "seq,\"d1,d2,d2\"\n"
        );
}

BOOST_AUTO_TEST_CASE(escape_quote) {
    cseq c("seq1");
    c.set_attr("data", "a quote '\"'");
    write(c);
    cseq d("seq2");
    d.set_attr("data", "more \"quotes\" '\"\"'");
    write(d);
    BOOST_CHECK_EQUAL(
        output(),
        "name,data\n"
        "seq1,\"a quote '\"\"'\"\n"
        "seq2,\"more \"\"quotes\"\" '\"\"\"\"'\"\n"
        );
}

BOOST_AUTO_TEST_CASE(escape_multilines) {
    cseq c("seq");
    c.set_attr("data", "multiple\nlines");
    write(c);
    BOOST_CHECK_EQUAL(
        output(),
        "name,data\n"
        "seq,\"multiple\nlines\"\n"
        );
}


BOOST_AUTO_TEST_CASE(cannot_write) {
    auto check_msg =
        [](const std::runtime_error &e) -> bool {
            return std::string(e.what()).rfind("Unable to open", 0) == 0;
        };
    BOOST_CHECK_EXCEPTION(
        rw_csv::writer("/proc/version", copy_relatives, fields),
        std::runtime_error,
        check_msg
        );
}


BOOST_AUTO_TEST_CASE(sep_tsv) {
    set_ext("tsv");
    cseq c("test_sequence");
    c.set_attr("col1", "Some test data");
    write(c);
    BOOST_CHECK_EQUAL(nlines(), 2);
    BOOST_CHECK_EQUAL(
        output(),
        "name\tcol1\n"
        "test_sequence\tSome test data\n"
        );
}


BOOST_AUTO_TEST_CASE(sep_tsv_override) {
    set_ext("tsv");
    configure({"--csv-sep", ";;"});
    cseq c("test_sequence");
    c.set_attr("col1", "Some test data");
    write(c);
    BOOST_CHECK_EQUAL(nlines(), 2);
    BOOST_CHECK_EQUAL(
        output(),
        "name;;col1\n"
        "test_sequence;;Some test data\n"
        );
}

BOOST_AUTO_TEST_CASE(nullptr_seq) {
    writer = new rw_csv::writer(*outfile, copy_relatives, fields);
    tray t;
    t.aligned_sequence = nullptr;
    (*writer)(t);
    cseq c("seq");
    write(c);
    (*writer)(t);
    BOOST_CHECK_EQUAL(
        output(),
        "name\nseq\n"
        );
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
