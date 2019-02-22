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

#include <unordered_set>

#include <boost/algorithm/string.hpp>
using boost::split;
using boost::is_any_of;

#include "query_arb.h"
#include "log.h"
#include "progress.h"

using namespace sina;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

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
         "skip the first n sequences (0)")
        ;

    adv.add(od);
}

void
rw_arb::validate_vm(po::variables_map& /*vm*/,
                    po::options_description& /*desc*/) {
}


/** Section: Reader
 */

struct rw_arb::reader::priv_data {
    query_arb* arb{nullptr};
    istream *in{nullptr};
    int seqno{0};
    int total_expected_sequences{0};
    vector<string>& v_fields;
    Progress *p{nullptr};
    explicit priv_data(vector<string>& fields)
        : v_fields(fields)
    {
    }

    ~priv_data() {
        logger->info("read {} sequences", seqno);
    }
};

rw_arb::reader::reader() = default;
rw_arb::reader::reader(const reader& o) = default;
rw_arb::reader::~reader() = default;
rw_arb::reader& rw_arb::reader::operator=(const reader& o) = default;

rw_arb::reader::reader(fs::path infile,
                       vector<string>& fields)
    : data(new priv_data(fields))
{
    data->arb = query_arb::getARBDB(infile);

    int n_seqs_db = data->arb->getSeqCount();
    int n_seqs_sel = 0;

    if (opts->select_file ==  "-") {
        data->in = &std::cin;
    } else if (!opts->select_file.empty()) {
        data->in = new ifstream(opts->select_file.c_str());
    } else {
        auto *tmp = new stringstream();
        vector<string> cache = data->arb->getSequenceNames();
        for (auto & it : cache) {
            *tmp << it << std::endl;
        }
        data->in = tmp;
        n_seqs_sel = n_seqs_db;
    }

    // ignore first <select_skip> names
    if (opts->select_skip) {
        logger->info("Skipping first {} sequences", opts->select_skip);
        string tmp;
        for (int i = 0; i < opts->select_skip; i++) {
            if (data->in->bad()) {
                logger->error("After skipping {} sequences, none where left", i);
                break;
            }
            (*data->in) >> tmp;
        }
        n_seqs_sel -= opts->select_skip;
    }

    if (opts->select_step > 1) {
        logger->info("Processing only every {}th sequence", opts->select_step);
        n_seqs_sel = 1 + (n_seqs_sel - 1) / opts->select_step;
    }

    if (n_seqs_sel > 0 && n_seqs_sel < n_seqs_db) {
        logger->info("Processing {} sequences out of {} in the input database",
                     n_seqs_sel, n_seqs_db);
    }   
    data->total_expected_sequences = (n_seqs_sel > 0) ? n_seqs_sel : 0;
}

void
rw_arb::reader::set_progress(Progress &p) {
    data->p = &p;
    data->p->set_total(data->total_expected_sequences);
}


bool
rw_arb::reader::operator()(tray& t) {
    string name;
    t.seqno = ++data->seqno;

    t.input_sequence = nullptr; // we may get initialized tray

    while (t.input_sequence == nullptr) {
        if (data->in->bad()) {
            data->total_expected_sequences = data->seqno - 1;
            data->p->set_total(data->total_expected_sequences);
            return false;
        }

        for (int i = 2; i <= opts->select_step; i++) {
            (*data->in) >> name;
        }
        (*data->in) >> name;
    
        if (name.empty()) {
            data->total_expected_sequences = data->seqno - 1;
            data->p->set_total(data->total_expected_sequences);
            return false;
        }

        try {
            t.input_sequence = new cseq(data->arb->getCseq(name));
        } catch (base_iupac::bad_character_exception& e) {
            --data->total_expected_sequences;
            data->p->set_total(data->total_expected_sequences);
            logger->error("Bad character {} in sequence {}",
                          e.character, name);
        }
    }
    for (const auto& f: data->v_fields) {
        data->arb->loadKey(*t.input_sequence, f);
    }

    logger->debug("loaded sequence {}", t.input_sequence->getName());
    return true;
}


/** Section: Writer
 */
struct rw_arb::writer::priv_data {
    query_arb *arb;
    fs::path arb_fname;
    int count;
    int excluded;
    std::unordered_set<string> relatives_written;
    unsigned long copy_relatives;
    vector<string>& v_fields;
    priv_data(fs::path& arb_fname_,
              unsigned long copy_relatives_,
              vector<string>& fields)
        : arb(nullptr),
          arb_fname(arb_fname_),
          count(0),
          excluded(0),
          copy_relatives(copy_relatives_),
          v_fields(fields)
    {
    }
    ~priv_data() {
        if (arb == nullptr) { // might never have been initialized
            return;
        }
        logger->info("wrote {} sequences ({} excluded, {} relatives)",
                     count, excluded, relatives_written.size());
        if (arb_fname.native() != ":") {
            arb->save();
        }
    }    
};

rw_arb::writer::writer() = default;
rw_arb::writer::writer(const writer& o) = default;
rw_arb::writer::~writer() = default;
rw_arb::writer& rw_arb::writer::operator=(const writer& o) = default;

rw_arb::writer::writer(fs::path outfile, unsigned int copy_relatives,
                       vector<string>& fields)
    :  data(new priv_data(outfile, copy_relatives, fields))
{
    data->arb = query_arb::getARBDB(outfile); // might throw
    data->arb->setProtectionLevel(opts->prot_lvl);
}


tray
rw_arb::writer::operator()(tray t) {
    if (t.aligned_sequence == nullptr) {
        logger->info("Not writing sequence {} (>{}): not aligned",
                     t.seqno, t.input_sequence->getName());
        ++data->excluded;
        return t;
    }
    cseq &c = *t.aligned_sequence;

    data->arb->putCseq(c);

    if (data->copy_relatives != 0u) {
        // FIXME: we should copy if reference is an arb database
        auto* relatives = t.search_result != nullptr ? t.search_result : t.alignment_reference;
        if (relatives != nullptr) {
            int i = data->copy_relatives;
            for (auto& seq : *relatives) {
                if (data->relatives_written.insert(seq.getName()).second) {
                    data->arb->putCseq(seq);
                    data->count++;
                }
                if (--i == 0) {
                    break;
                }
            }
        }
    }
    data->count++;
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
