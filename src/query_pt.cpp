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

#include "config.h"
#include "query_pt.h"

#include <cstring>

#include <iostream>
#include <sstream>
#include <cmath>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

using std::endl;
using std::stringstream;
using std::string;
using std::vector;

#include "cseq.h"
#include "log.h"

#ifndef TEST

#include "query_arb.h"

#include <dlfcn.h>
#include <cstdlib>

#include <arbdb.h>
#include <PT_com.h>
#include <client.h>
#include <boost/thread/mutex.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/dll.hpp>
#include <pstream.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

#include <boost/system/error_code.hpp>
#include <utility>
namespace sys = boost::system;

namespace sina {

static auto logger = Log::create_logger("Search (ARB PT)");
static auto pt_logger = Log::create_logger("ARB_PT_SERVER");


class managed_pt_server {
    redi::ipstream* process;
    string dbname;
    string portname;
public:
    managed_pt_server(string  dbname, string  portname);
    managed_pt_server(const managed_pt_server&);
    ~managed_pt_server();
};


managed_pt_server::managed_pt_server(string  dbname_, string  portname_)
    : dbname(std::move(dbname_)), portname(std::move(portname_))
{
    // Check that database specified and file accessible
    if (dbname.empty()) {
        throw query_pt_exception("Missing reference database");
    }

    struct stat arbdb_stat;
    if (stat(dbname.c_str(), &arbdb_stat)) {
        perror("Error accessing ARB database file");
        throw query_pt_exception("Failed to launch PT server.");
    }

    // Make sure ARBHOME is set; guess if possible
    fs::path ARBHOME;
    if (const char* arbhome = std::getenv("ARBHOME")) {
        ARBHOME = arbhome;
        logger->info("Using ARBHOME={}", ARBHOME);
    } else {
        ARBHOME = boost::dll::symbol_location(GB_open).parent_path().parent_path();
        logger->info("Setting ARBHOME={}", ARBHOME);
        setenv("ARBHOME", ARBHOME.c_str(), 1);  // no setenv in C++/STL
    }

    // Locate PT server binary
    fs::path arb_pt_server("arb_pt_server");
    fs::path arb_pt_server_path = fs::system_complete(arb_pt_server);
    if (!fs::exists(arb_pt_server_path)) { // not in PATH
        logger->debug("{} not found in PATH", arb_pt_server);
        arb_pt_server_path = ARBHOME / "bin" / arb_pt_server;
        if (!fs::exists(arb_pt_server_path)) { // not in ARBHOME
            logger->debug("{} not found in ARBHOME", arb_pt_server);
            // 3. Next to us
            fs::path self_dir = boost::dll::program_location().parent_path();
            arb_pt_server_path = self_dir / arb_pt_server;
            if (!fs::exists(arb_pt_server_path)) { // not next to us either?
                logger->debug("{} not found in {}", arb_pt_server, self_dir);
                throw query_pt_exception("Failed to locate 'arb_pt_server'");
            }
        }
    }

    // (Re)build index if missing or older than database
    struct stat ptindex_stat;
    string ptindex = dbname + ".index.arb.pt";
    if (stat(ptindex.c_str(), &ptindex_stat) 
        || arbdb_stat.st_mtime > ptindex_stat.st_mtime) {
        if (arbdb_stat.st_mtime > ptindex_stat.st_mtime) {
            logger->info("PT server index missing for {}. Building:", dbname);
        } else {
            logger->info("PT server index out of date for {}. Rebuilding:", dbname);
        }

        vector<string> cmds;
        cmds.push_back(string("cp ") + dbname + " " + dbname + ".index.arb");
        cmds.push_back(arb_pt_server_path.native() + " -build_clean -D" + dbname + ".index.arb");
        cmds.push_back(arb_pt_server_path.native() + " -build -D" + dbname + ".index.arb");
        for (auto& cmd :  cmds) {
            logger->debug("Executing '{}'", cmd);
            system(cmd.c_str());
            logger->debug("Command finished");
        }

        if (stat(ptindex.c_str(), &ptindex_stat) 
            || arbdb_stat.st_mtime > ptindex_stat.st_mtime) {
            throw query_pt_exception("Failed to (re)build PT server index! (out of memory?)");
        }
    }

    // Reset database name to the index version created during build
    dbname = dbname + ".index.arb";

    // Check portname: allowed are localhost:PORT, :PORT and :SOCKETFILE
    int split = portname.find(":");
    string host = portname.substr(0, split);
    string port = portname.substr(split+1);
    if (!host.empty() && host != "localhost") {
        throw query_pt_exception("Starting a PT server on hosts other than localhost not supported");
    }

    // Actually launch the server now:
    vector<string> cmd{arb_pt_server_path.native(), "-D"+dbname, "-T"+portname};
    logger->info("Launching background PT server for {} on {}", dbname, portname);
    logger->debug("Executing ['{}']", fmt::join(cmd, "', '"));
    process = new redi::ipstream(cmd, redi::pstreams::pstdout|redi::pstreams::pstderr);

    // read the pt server output. once it says "ok"
    // we can be sure that it's ready for connections
    // (connecting at wrong time causes lockup)
    // FIXME: abort if waiting for too long
    string line;
    while (std::getline(*process, line)) {
        pt_logger->debug(line);
        if (line == "ok, server is running.") {
            break;
        }
    }

    logger->info("Launched PT server ({} on {}).", dbname, portname);
}

managed_pt_server::~managed_pt_server() {
    logger->info("Terminating PT server ({} on {})", dbname, portname);
    process->rdbuf()->kill();
}

struct query_pt::options {
};
struct query_pt::options *query_pt::opts;

void
query_pt::get_options_description(po::options_description& /*main*/,
                                  po::options_description& /*adv*/) {
}

void
query_pt::validate_vm(po::variables_map& /* vm */,
                      po::options_description& /*desc*/) {
}


struct query_pt::priv_data {
    aisc_com         *link{nullptr};
    T_PT_MAIN         com;
    T_PT_LOCS         locs;
    T_PT_FAMILYFINDER ffinder;

    boost::mutex arb_pt_access;

    int  range_begin{-1};
    int  range_end{-1};
    bool find_type_fast{false};
    int  kmer_len;
    int  num_mismatch;
    bool relative_sort;

    static std::map<string, std::weak_ptr<managed_pt_server>> servers;
    std::shared_ptr<managed_pt_server> server;

    bool            connect_server(string portname);
    void            disconnect_server();
};


bool
query_pt::priv_data::connect_server(string portname) {
    boost::mutex::scoped_lock lock(arb_pt_access);
    GB_ERROR error = nullptr;
    link = aisc_open(portname.c_str(), com, AISC_MAGIC_NUMBER, &error);
    if (error) {
        throw query_pt_exception(error);
    }
    if (!link) {
        return false;
    }

    if (aisc_create(link,
                    PT_MAIN, com,
                    MAIN_LOCS, PT_LOCS, locs,
                    NULL)) {
        throw query_pt_exception("Unable to connect to PT server! (code 02)");
    }

    if (aisc_create(link,
                    PT_LOCS, locs,
                    LOCS_FFINDER, PT_FAMILYFINDER, ffinder,
                    NULL)) {
        throw query_pt_exception("Unable to connect to PT server! (code 03)");
    }

    return true;
}

void
query_pt::priv_data::disconnect_server() {
    boost::mutex::scoped_lock lock(arb_pt_access);
    aisc_close(link, com);
}



std::map<string, std::weak_ptr<managed_pt_server>> query_pt::priv_data::servers;


query_pt::query_pt(const char* portname, const char* dbname,
                   bool fast, int k, int mk, bool norel)
    : data(new priv_data())
{
    if (data->servers.count(portname)) {
        data->server = data->servers[portname].lock();
    }

    if (!data->connect_server(portname)) {
        data->server = std::make_shared<managed_pt_server>(dbname, portname);
        if (!data->connect_server(portname)) {
            throw query_pt_exception("Failed to start PT server. Do you have enough memory?");
        }
        data->servers[portname] = data->server;
    }

    set_find_type_fast(fast);
    set_probe_len(k);
    set_mismatches(mk);
    set_sort_type(norel);
}

query_pt::~query_pt() {
    delete data;
}

#if 0
void
query_pt::restart() {
    logger->info("Trying to restart pt server connection...");
    exit();
    logger->info("Terminated PT server. Sleeping for 5 seconds.");
    sleep(5);
    logger->info("Restarting PT server");
    init();
    logger->info("Done.");
    // FIXME: settings need to be restored!
}
#endif

void
query_pt::set_find_type_fast(bool fast) {
    boost::mutex::scoped_lock lock(data->arb_pt_access);
    int err = aisc_put(data->link,
                       PT_FAMILYFINDER, data->ffinder,
                       FAMILYFINDER_FIND_TYPE, fast?1:0,
                       NULL);
    if (err) {
        logger->warn("Unable to set find_type = {}", fast ? "fast" : "normal");
    } else {
        data->find_type_fast = fast;
    }
}

void
query_pt::set_probe_len(int len) {
    boost::mutex::scoped_lock lock(data->arb_pt_access);
    int err = aisc_put(data->link,
                       PT_FAMILYFINDER, data->ffinder,
                       FAMILYFINDER_PROBE_LEN, len,
                       NULL);
    if (err) {
        logger->warn("Unable to set k = {}", len);
    } else {
        data->kmer_len = len;
    }
}

void
query_pt::set_mismatches(int len) {
    boost::mutex::scoped_lock lock(data->arb_pt_access);
    int err = aisc_put(data->link,
                       PT_FAMILYFINDER, data->ffinder,
                       FAMILYFINDER_MISMATCH_NUMBER, len,
                       NULL);

    if (err) {
        logger->warn("Unable to set allowable mismatches to {}", len);
    } else {
        data->num_mismatch = len;
    }
}

void
query_pt::set_sort_type(bool absolute) {
    boost::mutex::scoped_lock lock(data->arb_pt_access);
    int err = aisc_put(data->link,
                       PT_FAMILYFINDER, data->ffinder,
                       FAMILYFINDER_SORT_TYPE, absolute?0:1,
                       NULL);
    if (err) {
        logger->warn("Unable to set sort type = {}", absolute ? "absolute" : "relative");
    } else {
        data->relative_sort = !absolute;
    }
}

void
query_pt::set_range(int startpos, int stoppos) {
    boost::mutex::scoped_lock lock(data->arb_pt_access);
#if 0
    int err = aisc_put(data->link,
                       PT_FAMILYFINDER, data->ffinder,
                       FAMILYFINDER_RANGE_STARTPOS, startpos,
                       FAMILYFINDER_RANGE_ENDPOS, stoppos,
                       NULL);
#else
    int err = 0;
#endif
    if (err) {
        logger->warn("Unable to constain matching to {}-{}", startpos, stoppos);
    } else {
        data->range_begin = startpos;
        data->range_end = stoppos;
    }
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
    int matches = 0;
    int skipped_max_score = 0;
    int skipped_broken = 0;
    int skipped_min_len = 0;
    int skipped_noid = 0;

    if (query_str.size() < 20) {
        logger->warn("Sequence {} too short ({} bp) for PT server search",
                     queryc.getName(), strlen(query));
        return 0;
    }

    T_PT_FAMILYLIST f_list;

    bytestring bs;
    bs.data = const_cast<char*>(query);
    bs.size = strlen(query)+1;

match_retry:
    family.reserve(max_match);
    boost::mutex::scoped_lock lock(data->arb_pt_access);
    
/*
    if (max_score > 1 && !noid) {
        // if we really only want the top max_match, 
        // limit sorting to that amount. if we have a max_score,
        // we might as well have the pt server sort it all.
        aisc_put(data->link,
                       PT_FAMILYFINDER, data->ffinder,
                 FAMILYFINDER_SORT_MAX, max_match,
                 NULL);
    } else {
        aisc_put(data->link,
                       PT_FAMILYFINDER, data->ffinder,
                 FAMILYFINDER_SORT_MAX, 0,
                 NULL);
    }
*/

    int err = aisc_put(data->link,
                       PT_FAMILYFINDER, data->ffinder,
                       FAMILYFINDER_FIND_FAMILY, &bs,
                       NULL);
    if (err) {
        logger->error("Unable to execute find_family command on pt-server");
        if (--maxfail==0) {
            logger->error("No retries left; aborting.");
            return 0;
        } else {
            logger->error("Retrying...");
        }
        //FIXME restart();
        goto match_retry;
    }

    err = aisc_get(data->link,
                   PT_FAMILYFINDER, data->ffinder,
                   FAMILYFINDER_FAMILY_LIST, f_list.as_result_param(),
                   NULL);
    if (err) {
        logger->error("Unable to get results for search");
        return 0;
    }
    if (!f_list.exists()) {
        return 0;
    }


    char   *f_name;
    double  f_relscore = 0.f;
    do {
        err = aisc_get(data->link, PT_FAMILYLIST, f_list,
                       FAMILYLIST_NAME, &f_name,
                       FAMILYLIST_REL_MATCHES, &f_relscore,
                       FAMILYLIST_NEXT, f_list.as_result_param(),
                       NULL);
        if (err) {
            logger->error("Unable to get next item in family list");
            break;
        }

        if (leave_query_out && f_name == queryc.getName()) {
            logger->info("Omitting sequence with same name as query from result");
            free(f_name);
            continue;
        }

        if (data->find_type_fast) {
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
                    logger->error("Sequence {} contains invalid character{}. Skipping",
                                  f_name, e.character);
                    sequence_broken=true;
                }
                seq.setScore(f_relscore);
                /*if (max_score <= 2 && queryc.identity_with(seq) <= max_score) {
                    skipped_max_score ++;
                    } else*/
                if (sequence_broken) {
                    skipped_broken ++;
                } else if ((long)seq.size() < min_len) {
                    skipped_min_len ++;
                } else if (noid && boost::algorithm::icontains(seq.getBases(), query)) {
                    skipped_noid ++;
                } else {
                    matches ++;

                    family.push_back(seq);

                    if (num_full && (long)seq.size() > full_min_len)
                        num_full--;
                    if (range_cover_right &&
                        seq.getById(seq.size()-1).getPosition() >= data->range_end)
                        range_cover_right--;
                    if (range_cover_left &&
                        seq.begin()->getPosition() <= data->range_begin)
                        range_cover_left--;
                }
            } else {
                if (f_relscore <= max_score) {
                    family.emplace_back(f_name, f_relscore);
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
            err = aisc_get(data->link, PT_FAMILYLIST, f_list,
                           FAMILYLIST_NAME, &f_name,
                           FAMILYLIST_MATCHES, &f_relscore,
                           FAMILYLIST_NEXT, f_list.as_result_param(),
                           NULL);

            if (err) {
                logger->warn("Unable to get next item in family list");
                break;
            }
            if (data->find_type_fast) {
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
                    seq.getById(seq.size()-1).getPosition() >= data->range_end) {
                    range_cover_right--;
                    keep = true;
                }

                if (range_cover_left && (long)seq.size() > min_len &&
                    seq.begin()->getPosition() <= data->range_begin) {
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

    if (skipped_max_score || skipped_broken || skipped_min_len || skipped_noid) {
        logger->warn("Skipped {} sequences ({} id < {}, {} broken, {} len < {}, {} noid)",
                     skipped_max_score + skipped_broken + skipped_min_len + skipped_noid,
                     skipped_max_score, max_score,
                     skipped_broken,
                     skipped_min_len, min_len,
                     skipped_noid);
    }

    return f_relscore;
}

query_pt_exception::query_pt_exception(std::string  msg) noexcept
    : message(std::move(msg))
{
}

query_pt_exception::~query_pt_exception() noexcept {
}

const char*
query_pt_exception::what() const noexcept {
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
