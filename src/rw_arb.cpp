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

#include "rw_arb.h"

#include <iostream>

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
#include "log.h"


using namespace sina;
namespace po = boost::program_options;

static const char* module_name = "ARB I/O";
static auto logger = Log::create_logger(module_name);


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

    po::options_description od(module_name);
    od.add_options()
        // write
        ("markcopied",
         po::bool_switch(&opts->markcopied),
         "mark copied references")
        ("markaligned",
         po::bool_switch(&opts->markaligned),
         "mark aligned sequences")
        ("prot-level",
         po::value<int>(&opts->prot_lvl)->default_value(4, ""),
         "arb export protection level (4)")

        // read
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
    split(opts->v_extra_fields, opts->extra_fields, is_any_of(":,"));
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

struct rw_arb::reader::priv_data {
    query_arb* arb;
    istream *in;

    priv_data() : arb(NULL), in(NULL) {}
};

rw_arb::reader::reader() {}
rw_arb::reader::reader(const reader& o) : data(o.data) {}
rw_arb::reader::~reader() {}
rw_arb::reader& rw_arb::reader::operator=(const reader& o) {
    data = o.data;
    return *this;
}


rw_arb::reader::reader(std::string infile)
    : data(new priv_data())
{
    data->arb = query_arb::getARBDB(infile);

    if (opts->select_file ==  "-") {
        data->in = &std::cin;
    } else if (!opts->select_file.empty()) {
        data->in = new ifstream(opts->select_file.c_str());
    } else {
        stringstream *tmp = new stringstream();
        vector<string> cache = data->arb->getSequenceNames();
        for (vector<string>::iterator it = cache.begin();
              it != cache.end(); ++it) {
            *tmp << *it << std::endl;
        }
        data->in = tmp;
    }

    // ignore first <select_skip> names
    string tmp;
    for (int i = 0; i < opts->select_skip; i++) {
        (*data->in) >> tmp;
    }
}


bool
rw_arb::reader::operator()(tray& t) {
    string name;
    t.input_sequence = 0; // we may get initialized tray

    while (not t.input_sequence) {
        if (data->in->bad()) {
            return false;
        }

        for (int i = 2; i <= opts->select_step; i++) {
            (*data->in) >> name;
        }
        (*data->in) >> name;
    
        if (name.empty()) {
            return false;
        }

        try {
            t.input_sequence = new cseq(data->arb->getCseq(name));
        } catch (base_iupac::bad_character_exception& e) {
            logger->error("Bad character {} in sequence {}",
                          e.character, name);
        }
    }

    if (!opts->extra_fields.empty()) {
        for (string& f: opts->v_extra_fields) {
            data->arb->loadKey(*t.input_sequence, f);
        }
    }

    logger->debug("loaded sequence {}", t.input_sequence->getName());
    return true;
}


/** Section: Writer
 */
struct rw_arb::writer::priv_data {
    query_arb *arb;
    string arb_fname;
    ~priv_data() {
        logger->info("Saving...");
        if (arb_fname != ":") {
            arb->save();
        }
    }    
};

rw_arb::writer::writer() {}
rw_arb::writer::writer(const writer& o) : data(o.data) {}
rw_arb::writer::~writer() {}
rw_arb::writer& rw_arb::writer::operator=(const writer& o) {
    data = o.data;
    return *this;
}

rw_arb::writer::writer(string outfile)
    :  data(new priv_data)
{
    data->arb_fname = outfile;
    data->arb = query_arb::getARBDB(outfile);
    data->arb->setProtectionLevel(opts->prot_lvl);
}


tray
rw_arb::writer::operator()(tray t) {
    if (t.aligned_sequence) {
        data->arb->putCseq(*t.aligned_sequence);
    }
    return t;
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
