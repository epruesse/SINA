/*
Copyright (c) 2006-2017 Elmar Pruesse <elmar.pruesse@ucdenver.edu>

This file is part of SINA.
123456789012345678901234567890123456789012345678901234567890123456789012
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

#include "rw_arb.h"

#include <iostream>
using std::cerr;
using std::endl;

#include <fstream>
using std::istream;
using std::ifstream;

#include <string>
using std::string;

#include <sstream>
using std::stringstream;

#include <vector>
using std::vector;

#include <map>
using std::map;

#include <boost/thread/thread.hpp>
using boost::thread;

#include <boost/shared_ptr.hpp>
using boost::shared_ptr;

#include <boost/bind.hpp>
using boost::bind;
using boost::ref;

#include <boost/algorithm/string.hpp>
using boost::split;
using boost::is_any_of;

#include "query_arb.h"
#include "align.h"


/** Writes newly aligned sequence to arb database
 * \param in input family stream, last entry of each family is stored
 * \param out output family stream, in is passed untouched
 * \param arb arb interface object, needed to access arb database
 */

using namespace sina;
namespace po = boost::program_options;


/** Section:  Configuration
 */

struct rw_arb::options {
    bool markcopied;
    bool markaligned;
    int prot_lvl;
    string select_file;
    int select_step;
    int select_skip;
    string extra_fields;
    vector<string> v_extra_fields;
    string pt_database;
};
struct rw_arb::options *rw_arb::opts;


void
rw_arb::get_options_description(po::options_description& /*main*/,
                                po::options_description& adv) {
    opts = new options();

    po::options_description od("ARB I/O");
    od.add_options()
        ("markcopied",
         po::bool_switch(&opts->markcopied),
         "mark copied references")

        ("markaligned",
         po::bool_switch(&opts->markaligned),
         "mark aligned sequences")

        ("prot-level",
         po::value<int>(&opts->prot_lvl)->default_value(4, ""),
         "arb export protection level (4)")

        ("select-file",
         po::value<string>(&opts->select_file)->default_value(""),
         "file containting arb names to be used ('-' for STDIN)")

        ("select-step",
         po::value<int>(&opts->select_step)->default_value(1, ""),
         "use every n-th sequence (1)" )

        ("select-skip",
         po::value<int>(&opts->select_skip)->default_value(0,""),
         "skip the first n seuqences (0)")

        ("extra-fields",
         po::value<string>(&opts->extra_fields)->default_value(""),
         "load additional fields, colon separated")
        ;

    adv.add(od);
}

void
rw_arb::validate_vm(po::variables_map& /*vm*/,
                    po::options_description& /*desc*/) {
//    opts->pt_database = vm["ptdb"].as<string>();
}


/** Section: Common
 */
string
make_datetime() {
    time_t  t;
    tm      tm;
    char   buf[50];

    time(&t);
    gmtime_r(&t, &tm);
    strftime(buf, 50, "%F %T", &tm);

    return string(buf);
}


/** Section: Reader
 */

class rw_arb::reader : public PipeElement<void, tray> {
    friend class rw_arb;
    reader(string infile);
    ~reader();
    query_arb* arb;
    istream *in;
public:
    tray operator()(void);
    std::string getName() const {return "arb::reader";}
};


PipeElement<void,tray>*
rw_arb::make_reader(po::variables_map& vm) {
    return new reader(vm["in"].as<string>());
}

rw_arb::reader::reader(string infile)
{
    arb = query_arb::getARBDB(infile);

    if (opts->select_file=="-") {
        in = &std::cin;
    } else if (!opts->select_file.empty()) {
        in = new ifstream(opts->select_file.c_str());
    } else {
        stringstream *tmp = new stringstream();
        vector<string> cache = arb->getSequenceNames();
        for (vector<string>::iterator it = cache.begin();
              it != cache.end(); ++it) {
            *tmp << *it << endl;
        }
        in = tmp;
    }

    // ignore first <select_skip> names
    string tmp;
    for (int i = 0; i < opts->select_skip; i++) {
        (*in) >> tmp;
        //cerr << "skipping " << tmp << endl;
    }

    split(opts->v_extra_fields, opts->extra_fields, is_any_of(":,"));
}

rw_arb::reader::~reader() {
}

tray
rw_arb::reader::operator()(void) {
    tray t;
    string name;
    t.logstream=new stringstream();
    
retry:
    if (in->bad())
        throw PipeEOF();
    for (int i = 2; i <= opts->select_step; i++) {
        (*in) >> name;
        //cerr << "skipping " << name << endl;
    }
    (*in) >> name;
    if (name.empty())
        throw PipeEOF();
    try {
        t.input_sequence = new cseq(arb->getCseq(name));
    } catch (base_iupac::bad_character_exception& e) {
        std::cerr << "ERROR: Bad character " << e.character << " in " << name
                  << ". Skipping sequence." << std::endl;
        goto retry;
    }
    if (!opts->extra_fields.empty()) {
        BOOST_FOREACH(string f, opts->v_extra_fields) {
            arb->loadKey(*t.input_sequence, f);
        }
    }
    return t;
}


/** Section: Writer
 */
class rw_arb::writer :  public PipeElement<tray, void> {
    friend class rw_arb;
    writer(string outfile);
    ~writer();
    query_arb *arb, *ptarb;
public:
    string arb_fname;
    void operator()(tray);
    std::string getName() const {return "arb::writer";}
private:
    int copyref;
};


PipeElement<tray,void>*
rw_arb::make_writer(po::variables_map& vm) {
    return new writer(vm["out"].as<string>());
}

rw_arb::writer::writer(string outfile)
    :  arb_fname(outfile)
{
    arb = query_arb::getARBDB(outfile);
    arb->setProtectionLevel(opts->prot_lvl);
    
    if (copyref) {
        ptarb = query_arb::getARBDB(opts->pt_database);
    }
}

rw_arb::writer::~writer() {
    if (arb_fname != ":") // cannot save remote db
        arb->save();
}

void
rw_arb::writer::operator()(tray t) {
    if (t.aligned_sequence == 0) {
    } else {
        arb->putCseq(*t.aligned_sequence);
    }
    t.destroy();
}

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
