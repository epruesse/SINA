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

#include "config.h"
#include "query_pt.h"

#include <string.h>

#include <iostream>
#include <sstream>
#include <cmath>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

using std::endl;
using std::cerr;
using std::stringstream;
using std::string;
using std::vector;

#include "cseq.h"

#ifndef TEST

#include "query_arb.h"

#include <arbdb.h>
#include <PT_com.h>
#include <client.h>
#include <boost/thread/mutex.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include "../include/pstream.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;


boost::mutex arb_pt_start;

namespace sina {

struct query_pt::options {
};
struct query_pt::options *query_pt::opts;

po::options_description
query_pt::get_options_description() {
    po::options_description od("ARB PT-server search");
    opts =  new struct query_pt::options();
    od.add_options()
        ;
}

void query_pt::validate_vm(po::variables_map& vm) {
}


struct query_pt::priv_data {
    priv_data(const char* port, const char* db)
        : portname(port),
          dbname(db),
          arb_pt_server(NULL),
          range_begin(-1),
          range_end(-1),
          find_type_fast(false)
    {}
    aisc_com *link;
    T_PT_MAIN com;
    T_PT_LOCS locs;
    T_PT_FAMILYFINDER ffinder;
    boost::mutex arb_pt_access;
    string portname;
    string dbname;
    redi::pstream *arb_pt_server;
    int range_begin;
    int range_end;
    bool find_type_fast;

    void create_context();
};

void
query_pt::priv_data::create_context() {
    if (!link) {
        throw exception("Could not register connection context with PT server");
    }

    boost::mutex::scoped_lock lock(arb_pt_start);
    if (aisc_create(link,
                    PT_MAIN, com,
                    MAIN_LOCS, PT_LOCS, locs,
                    NULL)) {
        throw exception("Unable to connect to PT server! (code 02)");
    }

    if (aisc_create(link,
                    PT_LOCS, locs,
                    LOCS_FFINDER, PT_FAMILYFINDER, ffinder,
                    NULL)) {
        throw exception("Unable to connect to PT server! (code 03)");
    }
}


void
query_pt::init() {
    // Try to connect to a PT server on the given port. Maybe it's already
    // running
    {  
        boost::mutex::scoped_lock lock(arb_pt_start);
        GB_ERROR error = NULL;
        data.link = aisc_open(data.portname.c_str(), data.com, AISC_MAGIC_NUMBER, &error);
        if (error) {
           throw exception(error);
        }
    }

    // If that worked, create context and return.
    if (data.link) { 
        data.create_context();
        return;
    }

    // Check if we have a database file.
    if (data.dbname.empty()) {
        // no chance to go on without one
        throw exception("Missing reference database");
    }

    struct stat ptindex_stat, arbdb_stat;
    string ptindex = data.dbname + ".index.arb.pt";
    if (stat(data.dbname.c_str(), &arbdb_stat)) {
        perror("Error accessing ARB database file");
        throw exception("Failed to launch PT server.");
    }
    if (stat(ptindex.c_str(), &ptindex_stat) 
        || arbdb_stat.st_mtime > ptindex_stat.st_mtime) {
        if (arbdb_stat.st_mtime > ptindex_stat.st_mtime) {
            cerr << "PT server index missing. Building... " << endl;
        } else {
            cerr << "PT server index out of date. Rebuilding..." << endl;
        }
        
        string cmd = string("cp ") + data.dbname + " " + data.dbname + ".index.arb &&"
            + "arb_pt_server -build_clean -D" + data.dbname + ".index.arb &&"
            + "arb_pt_server -build -D" + data.dbname + ".index.arb";

        system(cmd.c_str());

        if (stat(ptindex.c_str(), &ptindex_stat) 
            || arbdb_stat.st_mtime > ptindex_stat.st_mtime) {
            throw exception("Failed to (re)build PT server index! (out of memory?)");
        }
    }

    data.dbname = data.dbname + ".index.arb";

    int split = data.portname.find(":");
    string host = data.portname.substr(0, split);
    string port = data.portname.substr(split+1);
    if (!host.empty() && host != "localhost") {
        throw exception("Starting a PT server on hosts other than localhost not supported");
    }

    string cmd = string("arb_pt_server -D") + data.dbname + " -T" + data.portname;
    cerr << "Launching PT server..."
         << cmd << endl;
    data.arb_pt_server = new redi::pstream(cmd);

    // read the pt server output. once it says "ok"
    // we can be sure that it's ready for connections
    // (connecting at wrong time causes lockup)
    // FIXME: abort if waiting for too long
    // FIXME: the lockup should be fixed in ARB
    string line;
    while (std::getline(*data.arb_pt_server, line)) {
        cerr << line << endl;
        if (line == "ok, server is running.") break;
    }

    cerr << "Launched PT server. Connecting... ";
    {
        boost::mutex::scoped_lock lock(arb_pt_start);
        GB_ERROR error = NULL;
        data.link = aisc_open(data.portname.c_str(), 
                              data.com, AISC_MAGIC_NUMBER, &error);
        if (error) {
           throw exception(error);
        }
    }

    if (!data.link) {
        cerr << "[FAIL]" << endl;
        throw("Failed to start PT server. Do you have enough memory?");
    }
    cerr << "[OK]" << endl;

    data.create_context();
}

void
query_pt::exit() {
    if (data.arb_pt_server) { // we started our own pt server.
        cerr << "Terminating PT server..." << endl;
        aisc_nput(data.link, PT_MAIN, data.com, MAIN_SHUTDOWN, "47@#34543df43%&3667gh", NULL);
        // probably should wait here. delete below /should/ cause SIGHUP to child...
        delete data.arb_pt_server; 
    } 
    aisc_close(data.link, data.com);
    delete &data;
}

void
query_pt::restart() {
    cerr << "Trying to restart pt server connection..." << endl;
    exit();
    cerr << "Terminated aisc. Sleeping for 5." << endl;
    sleep(5);
    init();
    cerr << "Done. Hopefully." << endl;
    // FIXME: settings need to be restored!
}

query_pt::query_pt(const char* portname, const char* dbname)
    : data(*(new priv_data(portname,dbname)))
{
    init();
}

query_pt::~query_pt() {
    exit();
}

void
query_pt::set_find_type_fast(bool fast) {
    boost::mutex::scoped_lock lock(data.arb_pt_access);
    int err = aisc_put(data.link, 
                       PT_FAMILYFINDER, data.ffinder,
                       FAMILYFINDER_FIND_TYPE, fast?1:0,
                       NULL);
    if (err) cerr << "Unable to set find_type" << endl;
    else data.find_type_fast = fast;
}

void
query_pt::set_probe_len(int len) {
    boost::mutex::scoped_lock lock(data.arb_pt_access);
    int err = aisc_put(data.link,
                       PT_FAMILYFINDER, data.ffinder,
                       FAMILYFINDER_PROBE_LEN, len,
                       NULL);
    if (err) cerr << "Unable to set probe len" << endl;
}

void
query_pt::set_mismatches(int len) {
    boost::mutex::scoped_lock lock(data.arb_pt_access);
    int err = aisc_put(data.link,
                       PT_FAMILYFINDER, data.ffinder,
                       FAMILYFINDER_MISMATCH_NUMBER, len,
                       NULL);

    if (err) cerr << "Unable to set mismatch number" << endl;
}

void
query_pt::set_sort_type(bool absolute) {
    boost::mutex::scoped_lock lock(data.arb_pt_access);
    int err = aisc_put(data.link,
                       PT_FAMILYFINDER, data.ffinder,
                       FAMILYFINDER_SORT_TYPE, absolute?0:1,
                       NULL);
    if (err) cerr << "Unable to set sort type" << endl;
}

void
query_pt::set_range(int startpos, int stoppos) {
    boost::mutex::scoped_lock lock(data.arb_pt_access);
/*
    int err = aisc_put(data.link,
                       PT_FAMILYFINDER, data.ffinder,
                       FAMILYFINDER_RANGE_STARTPOS, startpos,
                       FAMILYFINDER_RANGE_ENDPOS, stoppos,
                       NULL);
    if (err) cerr << "Unable to set matching range" << endl;
*/

    data.range_begin = startpos;
    data.range_end = stoppos;
}

void
query_pt::unset_range() {
    set_range(-1,-1);
}

double
query_pt::match(std::vector<cseq> &family, const cseq& queryc,
                int min_match, int max_match, float min_score, float max_score,
                query_arb *arb, bool noid, int min_len,
                int num_full, int full_min_len, int range_cover, bool leave_query_out
                ) {
    string query_str = queryc.getBases();
    const char* query = query_str.c_str();
    int maxfail = 5;
    int range_cover_left = range_cover;
    int range_cover_right = range_cover;

    if (!query || strlen(query) < 20) {
        cerr << "Sequence too short (" << strlen(query)
             << " bases) for PT server search!" << endl;
        return 0;
    }

    T_PT_FAMILYLIST f_list;

    bytestring bs;
    bs.data = const_cast<char*>(query);
    bs.size = strlen(query)+1;

match_retry:
    family.reserve(max_match);
    boost::mutex::scoped_lock lock(data.arb_pt_access);
    
/*
    if (max_score > 1 && !noid) {
        // if we really only want the top max_match, 
        // limit sorting to that amount. if we have a max_score,
        // we might as well have the pt server sort it all.
        aisc_put(data.link, 
                       PT_FAMILYFINDER, data.ffinder,
                 FAMILYFINDER_SORT_MAX, max_match,
                 NULL);
    } else {
        aisc_put(data.link,
                       PT_FAMILYFINDER, data.ffinder,
                 FAMILYFINDER_SORT_MAX, 0,
                 NULL);
    }
*/

    int err = aisc_put(data.link,
                       PT_FAMILYFINDER, data.ffinder,
                       FAMILYFINDER_FIND_FAMILY, &bs,
                       NULL);
    if (err) {
        cerr << "Unable to execute find_family command on pt-server" 
             << endl;
        if (--maxfail==0) {
            cerr << "Retried too often. Aborting." << endl;
            return 0;
        }
        restart();
        goto match_retry;
    }

    err = aisc_get(data.link,
                   PT_FAMILYFINDER, data.ffinder,
                   FAMILYFINDER_FAMILY_LIST, f_list.as_result_param(),
                   NULL);
    if (err || !f_list.exists() ) {
        return 0;
    }

    double f_relscore = 0.f;
    int matches = 0, skipped = 0;
    char *f_name;

    do {
        err = aisc_get(data.link, PT_FAMILYLIST, f_list,
                       FAMILYLIST_NAME, &f_name,
                       FAMILYLIST_REL_MATCHES, &f_relscore,
                       FAMILYLIST_NEXT, f_list.as_result_param(),
                       NULL);
        if (err) {
            cerr << "Unable to get next item in family list" << endl;
            break;
        }

        if (leave_query_out && f_name == queryc.getName()) {
            cerr << "Leaving out query." << endl;
            free(f_name);
            continue;
        }

        if (find_type_fast) {
            f_relscore *= 4;
        }

        // correct relscore according to formula based on edgar2004
        f_relscore = 1 - log(f_relscore + 1.0/bs.size)/log(1.0/bs.size);

        if (matches <= min_match || f_relscore >= min_score) {
            if (arb) {
                bool sequence_broken=false;
                cseq seq(f_name);
                try {
                    seq = arb->getCseq(f_name);
                } catch (base_iupac::bad_character_exception& e) {
                    cerr << "Reference Sequence " << f_name 
                         << " contains invalid character " << e.character 
                         << ". Skipping." << endl;
                    sequence_broken=true;
                }
                seq.setScore(f_relscore);
                if (max_score >= 2 /*|| queryc.identity_with(seq) <= max_scoreFIXME*/) {
                    if (!sequence_broken && (long)seq.size() > min_len &&
                        (!noid || !boost::algorithm::icontains(seq.getBases(),query)) ) {
                        family.push_back(seq);

                        ++matches;
                        if (num_full && (long)seq.size() > full_min_len)
                            num_full--;
                        if (range_cover_right && 
                            seq.getById(seq.size()-1).getPosition() >= data.range_end)
                            range_cover_right--;
                        if (range_cover_left &&
                            seq.begin()->getPosition() <= data.range_begin)
                            range_cover_left--;
                    } else {
                        skipped++;
                    }
                }
            } else {
                if (f_relscore <= max_score) {
                    family.push_back(cseq(f_name, f_relscore));
                    ++matches;
                }
            }
        }

        free(f_name);
    } while (matches < max_match
             && (matches <= min_match || f_relscore >= min_score)
             && f_list.exists());

    // get full length sequence
    if (arb) {
        while (f_list.exists() && num_full + range_cover_right + range_cover_left) {
            err = aisc_get(data.link, PT_FAMILYLIST, f_list,
                           FAMILYLIST_NAME, &f_name,
                           FAMILYLIST_MATCHES, &f_relscore,
                           FAMILYLIST_NEXT, f_list.as_result_param(),
                           NULL);

            if (err) {
                cerr << "Unable to get next item in family list" << endl;
                break;
            }
            if (find_type_fast) {
                f_relscore *= 4;
            }

            cseq seq = arb->getCseq(f_name);
            seq.setScore(f_relscore);
            if (max_score >= 2 /*|| queryc.identity_with(seq) <= max_score FIXME*/) {
                bool keep = false;
                if (num_full && (long)seq.size() > full_min_len) {
                    num_full--;
                    keep = true;
                }

                if (range_cover_right && (long)seq.size() > min_len &&
                    seq.getById(seq.size()-1).getPosition() >= data.range_end) {
                    range_cover_right--;
                    keep = true;
                }

                if (range_cover_left && (long)seq.size() > min_len &&
                    seq.begin()->getPosition() <= data.range_begin) {
                    range_cover_left--;
                    keep = true;
                } 
                if (keep) {
                    family.push_back(seq);
                } 
            }
            free(f_name);
        }
    }

    if (skipped>0) {
        std::cerr << "ptmatch skipped " << skipped << " sequences" << endl;
    }

    return f_relscore;
}

int
query_pt::turn_check(const cseq& query, bool all) {
    std::vector<cseq> matches;
    double score[4];

    score[0] = match(matches, query, 1, 1, 0.0f);

    cseq turn(query);
    turn.reverse();
    if (all) {
        score[1] = match(matches, turn, 1, 1, 0.0f);

        cseq comp(query);
        comp.complement();
        score[2] = match(matches, comp, 1, 1, 0.0f);
    } else {
        score[1] = score[2] = 0;
    }

    turn.complement();
    score[3] = match(matches, turn, 1, 1, 0.0f);

    double max = 0;
    int best = 0;
    for (int i = 0; i < 4; i++)
        if (max < score[i])
            max = score[i], best = i;

    return best;
}


query_pt::exception::exception(std::string msg) throw()
    : message(msg)
{
}

query_pt::exception::~exception() throw() {
}

const char*
query_pt::exception::what() const throw() {
    return message.c_str();
}

} // namespace sina

#else // TEST

#include <string>
#include <fstream>
#include <algorithm>
#include <iterator>

void
read_fasta(istream& ifs, vector<cseq>& seqs) {
  while(!ifs.eof()) {
    char buf[1024];
    ifs.getline(buf,1024);
    if (*buf == '>') {
      string s(buf+1);
      while (ifs.fail() && !ifs.eof()) {
        ifs.clear();
        ifs.getline(buf,1024);
        s.append(buf);
      }
      seqs.push_back(cseq(s.c_str()));

    } else {
        if (seqs.empty()) continue; // malformed input, data before ">"
      cseq& s=seqs.back();
      s.append(buf);
      while (ifs.fail() && !ifs.eof()) {
        ifs.clear();
        ifs.getline(buf,1024);
        s.append(buf);
      }
    }
  }
}


int main(int argc, char **argv) {
    if (argc!=4)  exit(1);
    query_pt *pt = new query_pt(argv[2], argv[1]);

    std::ifstream ifs(argv[3]);
    std::vector<cseq> query, family;

    read_fasta(ifs, query);

    std::string seq = query[0].getBases();
    pt->match(family, seq.c_str(), 5, 40, .7);

    std::copy(family.begin(), family.end(),
              std::ostream_iterator<cseq>(std::cout, " "));

    delete pt;
    return 0;
}

#endif // TEST

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
