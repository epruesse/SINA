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

#define BOOST_TEST_MODULE rw

#include <boost/test/unit_test.hpp>
namespace utf = boost::unit_test;
#include <boost/test/data/test_case.hpp>
namespace bdata = utf::data;
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

BOOST_AUTO_TEST_SUITE(rw_csv_test);

using namespace sina;


class TempFile : boost::noncopyable {
public:
    TempFile(const fs::path model = "sina-%%%%-%%%%-%%%%")
        : _path(fs::temp_directory_path() / fs::unique_path(model))
    {}
    ~TempFile() {
        fs::remove(_path);
    }
    operator const fs::path() const {
        return _path;
    }
    const fs::path& path() const {
        return _path;
    }
    void dump(std::ostream& out = std::cerr) const {
        out << "Dumping Tempfile " << this << std::endl;
        std::string sep = "----------------------";
        out << sep << std::endl;
        fs::ifstream file(*this);
        std::string line;
        while (getline(file, line).good()) {
            std::cerr << line << std::endl;
        }
        out << sep << std::endl;
    }
    std::string load() const {
        auto closer = [](FILE* fp){ if(fp) std::fclose(fp);};
        auto fp = std::unique_ptr<FILE, decltype(closer)>(
            std::fopen(_path.native().c_str(), "r"), closer);
        std::string res;
        res.resize(fs::file_size(_path));
        auto len = std::fread(const_cast<char*>(res.data()), 1, res.size(), fp.get());
        return res;
    }
private:
    fs::path _path;
};
std::ostream& operator<<(std::ostream& out, const TempFile& t) { return out << t.path(); }


struct F {
    TempFile outfile;
    rw_csv::writer* writer{nullptr};
    int copy_relatives{0};
    std::vector<std::string> fields;

    std::string _output;
    int _nlines;

    void write(const cseq& c) {
        if (!writer) {
            writer = new rw_csv::writer(outfile, copy_relatives, fields);
        }
        tray t;
        t.aligned_sequence = const_cast<cseq*>(&c);
        (*writer)(t);
    }

    std::string& output() {
        if (writer) {
            delete writer;
            _output = outfile.load();
            writer = nullptr;
        }
        return _output;
    }

    int nlines() {
        auto txt = output();
        return std::count(txt.begin(), txt.end(), '\n');
    }
};

BOOST_FIXTURE_TEST_CASE(one_seq_no_data, F) {
    write(cseq("test_sequence"));
    BOOST_CHECK_EQUAL(nlines(), 2);
    BOOST_CHECK_EQUAL(
        output(),
        "name\n"
        "test_sequence\n"
        );
}

BOOST_FIXTURE_TEST_CASE(two_seq_no_data, F) {
    write(cseq("test_sequence1"));
    write(cseq("test_sequence2"));
    BOOST_CHECK_EQUAL(nlines(), 3);
    BOOST_CHECK_EQUAL(
        output(),
        "name\n"
        "test_sequence1\n"
        "test_sequence2\n");
}

BOOST_FIXTURE_TEST_CASE(one_seq_one_string, F) {
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

BOOST_FIXTURE_TEST_CASE(one_seq_one_int, F) {
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

BOOST_FIXTURE_TEST_CASE(escape_seqname, F) {
    cseq c("test_se,quence");
    write(c);
    BOOST_CHECK_EQUAL(
        output(),
        "name\n"
        "\"test_se,quence\"\n"
        );
}

BOOST_FIXTURE_TEST_CASE(escape_comma, F) {
    cseq c("seq");
    c.set_attr("data", "d1,d2,d2");
    write(c);
    BOOST_CHECK_EQUAL(
        output(),
        "name,data\n"
        "seq,\"d1,d2,d2\"\n"
        );
}

BOOST_FIXTURE_TEST_CASE(escape_quote, F) {
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

BOOST_FIXTURE_TEST_CASE(escape_multilines, F) {
    cseq c("seq");
    c.set_attr("data", "multiple\nlines");
    write(c);
    BOOST_CHECK_EQUAL(
        output(),
        "name,data\n"
        "seq,\"multiple\nlines\"\n"
        );
}


// TODO:
// - test trying to write to unwriteable path
// - test writing gzip


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
