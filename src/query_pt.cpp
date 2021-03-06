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
#include "timer.h"

#include <cstring>

#include <iostream>
#include <sstream>
#include <cmath>
#include <mutex>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

using std::stringstream;
using std::string;
using std::vector;

#include "cseq.h"
#include "cseq_comparator.h"
#include "log.h"

#include "query_arb.h"

#include <dlfcn.h>
#include <cstdlib>

#include <arbdb.h>
#include <PT_com.h>
#include <client.h>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/dll.hpp>
#include <pstream.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

#include <utility>

namespace sina {

static auto logger = Log::create_logger("Search (ARB PT)");
static auto pt_logger = Log::create_logger("ARB_PT_SERVER");


class managed_pt_server {
    redi::ipstream* process;
    const fs::path dbname;
    const string portname;
public:
    managed_pt_server(fs::path dbname_, string portname_);
    managed_pt_server(const managed_pt_server&);
    ~managed_pt_server();
    fs::path ensure_env_ARBHOME();
    fs::path get_pt_server_path();
    fs::path ensure_index_exists();
    bool build_index(fs::path index_arb_file);
};

static std::mutex build_lock, launch_lock;

managed_pt_server::managed_pt_server(fs::path dbname_, string  portname_)
    : dbname(std::move(dbname_)), portname(std::move(portname_))
{
    // Check that database specified and file accessible
    if (dbname.empty() or not fs::exists(dbname)) {
        throw query_pt_exception("Missing reference database");
    }

    ensure_env_ARBHOME();
    fs::path arb_pt_server = get_pt_server_path();
    fs::path index_arb_file = ensure_index_exists();

    // Check portname: allowed are localhost:PORT, :PORT and :SOCKETFILE
    int split = portname.find(':');
    string host = portname.substr(0, split);
    string port = portname.substr(split+1);
    if (!host.empty() && host != "localhost") {
        throw query_pt_exception("Starting a PT server on hosts other than localhost not supported");
    }

    // Actually launch the server now:
    vector<string> cmd{arb_pt_server.native(), "-D" + index_arb_file.native(), "-T" +portname};
    logger->info("Launching background PT server for {} on {}", dbname, portname);
    logger->debug("Executing ['{}']", fmt::join(cmd, "', '"));
    {
        std::lock_guard<std::mutex> lock(launch_lock);
        // something in here appears to not be thread safe
        process = new redi::ipstream(cmd, redi::pstreams::pstdout|redi::pstreams::pstderr);
    }

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
    if (line.empty()) {
        process->rdbuf()->kill();
        throw query_pt_exception("PT server failed to respond. Do you have enough memory?");
    }
    if (process->rdbuf()->exited()) {
        throw query_pt_exception(
            fmt::format("PT server exited immediately. Exit status was {}",
                        process->rdbuf()->status()));
    }

    logger->warn("Launched PT server ({} on {}).", dbname, portname);
}

managed_pt_server::~managed_pt_server() {
    logger->warn("Terminating PT server ({} on {})", dbname, portname);
    process->rdbuf()->kill();
}

fs::path
managed_pt_server::ensure_index_exists() {
    std::lock_guard<std::mutex> lock(build_lock);

    fs::path index_arb_file = dbname;
    index_arb_file += ".index.arb";
    fs::path index_pt_file = index_arb_file;
    index_pt_file += ".pt";


    if (fs::exists(index_pt_file) && fs::exists(index_arb_file)) {
        if (fs::last_write_time(index_arb_file) >= fs::last_write_time(dbname)
            && fs::last_write_time(index_pt_file) >= fs::last_write_time(index_arb_file)) {
            // some.arb older-than some.arb.index.arb older-than some.arb.index.arb.pt
            // -> good index exists
            return index_arb_file;
        } else {
            logger->info("PT server index out of date for {}. Rebuilding:", dbname);
        }
    } else {
        logger->info("PT server index missing or incomplete for {}. Building:", dbname);
    }

    build_index(index_arb_file);

    if (not fs::exists(index_pt_file) ||
        fs::last_write_time(index_pt_file) < fs::last_write_time(dbname)) {
        throw query_pt_exception("Failed to (re)build PT server index! (out of memory?)");
    }
    return index_arb_file;
}

bool managed_pt_server::build_index(fs::path index_arb_file) {
    fs::path arb_pt_server = get_pt_server_path();
    vector<string> cmds;
    // copy database file
    cmds.push_back(string("cp ") + dbname.native() + " " + index_arb_file.native());
    // shrink database (convert to index db)
    cmds.push_back(arb_pt_server.native() + " -build_clean -D" + index_arb_file.native());
    // build index (build .arb.pt)
    cmds.push_back(arb_pt_server.native() + " -build -D" + index_arb_file.native());
    for (auto& cmd : cmds) {
        logger->debug("Executing '{}'", cmd);
        int rval = system(cmd.c_str());
        logger->debug("Command finished");
        if (rval) {
            logger->error("Command {} failed with exit code {}", cmd, rval);
            return false;
        }
    }
    return true;
}


fs::path
managed_pt_server::ensure_env_ARBHOME() {
    fs::path ARBHOME;
    // Make sure ARBHOME is set; guess if possible
    if (const char* arbhome = std::getenv("ARBHOME")) {
        ARBHOME = arbhome;
        logger->info("Using ARBHOME={}", ARBHOME);
    } else {
        ARBHOME = boost::dll::symbol_location(GB_open).parent_path().parent_path();
        logger->info("Setting ARBHOME={}", ARBHOME);
        setenv("ARBHOME", ARBHOME.c_str(), 1);  // no setenv in C++/STL
    }
    return ARBHOME;
}

fs::path
managed_pt_server::get_pt_server_path() {
    // Locate PT server binary
    fs::path arb_pt_server("arb_pt_server");
    fs::path arb_pt_server_path = fs::system_complete(arb_pt_server);
    if (!fs::exists(arb_pt_server_path)) { // not in PATH
        logger->debug("{} not found in PATH", arb_pt_server);
        arb_pt_server_path = ensure_env_ARBHOME() / "bin" / arb_pt_server;
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
    return arb_pt_server_path;
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

    std::mutex arb_pt_access;

    unsigned int  range_begin{0};
    unsigned int  range_end{INT_MAX};
    bool find_type_fast{false};
    int  kmer_len;
    int  num_mismatch;
    bool relative_sort;

    query_arb *arbdb;

    timer timeit;

    static std::map<string, std::weak_ptr<managed_pt_server>> servers;
    std::shared_ptr<managed_pt_server> server;

    bool            connect_server(const string& portname);
    void            disconnect_server();
};


bool
query_pt::priv_data::connect_server(const string& portname) {
    std::lock_guard<std::mutex> lock(arb_pt_access);
    GB_ERROR error = nullptr;
    link = aisc_open(portname.c_str(), com, AISC_MAGIC_NUMBER, &error);
    if (error != nullptr) {
        throw query_pt_exception(error);
    }
    if (link == nullptr) {
        return false;
    }

    if (aisc_create(link,
                    PT_MAIN, com,
                    MAIN_LOCS, PT_LOCS, locs,
                    NULL) != 0) {
        throw query_pt_exception("Unable to connect to PT server! (code 02)");
    }

    if (aisc_create(link,
                    PT_LOCS, locs,
                    LOCS_FFINDER, PT_FAMILYFINDER, ffinder,
                    NULL) != 0) {
        throw query_pt_exception("Unable to connect to PT server! (code 03)");
    }

    return true;
}

void
query_pt::priv_data::disconnect_server() {
    std::lock_guard<std::mutex> lock(arb_pt_access);
    aisc_close(link, com);
}


std::map<string, std::weak_ptr<managed_pt_server>> query_pt::priv_data::servers;


query_pt*
query_pt::get_pt_search(const fs::path& filename, int k,
                        bool fast,
                        bool norel,
                        int mk,
                        std::string portname) {
    if (portname.empty()) {
        portname = ":" + (fs::temp_directory_path() / fs::unique_path()).native();
    }

    return new query_pt(portname.c_str(), filename.native().c_str(),
                        fast, k, mk, norel);
}


query_pt::query_pt(const char* portname, const char* dbname,
                   bool fast, int k, int mk, bool norel)
    : data(new priv_data())
{
    if (priv_data::servers.count(portname) != 0u) {
        data->server = priv_data::servers[portname].lock();
    }

    if (!data->connect_server(portname)) {
        data->server = std::make_shared<managed_pt_server>(dbname, portname);
        if (!data->connect_server(portname)) {
            delete data;
            throw query_pt_exception("Failed to start PT server. Do you have enough memory?");
        }
        priv_data::servers[portname] = data->server;
    }

    data->arbdb = query_arb::getARBDB(dbname);

    set_find_type_fast(fast);
    set_probe_len(k);
    set_mismatches(mk);
    set_sort_type(norel);
}

query_pt::~query_pt() {
    logger->info("Timings for PT Search: {}", data->timeit);
    delete data;
}


void
query_pt::set_find_type_fast(bool fast) {
    std::lock_guard<std::mutex> lock(data->arb_pt_access);
    int err = aisc_put(data->link,
                       PT_FAMILYFINDER, data->ffinder,
                       FAMILYFINDER_FIND_TYPE, fast?1L:0L,
                       NULL);
    if (err != 0) {
        logger->warn("Unable to set find_type = {}", fast ? "fast" : "normal");
    } else {
        data->find_type_fast = fast;
    }
}

void
query_pt::set_probe_len(int len) {
    std::lock_guard<std::mutex> lock(data->arb_pt_access);
    int err = aisc_put(data->link,
                       PT_FAMILYFINDER, data->ffinder,
                       FAMILYFINDER_PROBE_LEN, long(len),
                       NULL);
    if (err != 0) {
        logger->warn("Unable to set k = {}", len);
    } else {
        data->kmer_len = len;
    }
}

void
query_pt::set_mismatches(int len) {
    std::lock_guard<std::mutex> lock(data->arb_pt_access);
    int err = aisc_put(data->link,
                       PT_FAMILYFINDER, data->ffinder,
                       FAMILYFINDER_MISMATCH_NUMBER, long(len),
                       NULL);

    if (err != 0) {
        logger->warn("Unable to set allowable mismatches to {}", len);
    } else {
        data->num_mismatch = len;
    }
}

void
query_pt::set_sort_type(bool absolute) {
    std::lock_guard<std::mutex> lock(data->arb_pt_access);
    int err = aisc_put(data->link,
                       PT_FAMILYFINDER, data->ffinder,
                       FAMILYFINDER_SORT_TYPE, absolute?0L:1L,
                       NULL);
    if (err != 0) {
        logger->warn("Unable to set sort type = {}", absolute ? "absolute" : "relative");
    } else {
        data->relative_sort = !absolute;
    }
}

void
query_pt::set_range(int startpos, int stoppos) {
    std::lock_guard<std::mutex> lock(data->arb_pt_access);

    int err = aisc_put(data->link,
                       PT_FAMILYFINDER, data->ffinder,
                       FAMILYFINDER_RANGE_STARTPOS, long(startpos),
                       FAMILYFINDER_RANGE_ENDPOS, long(stoppos),
                       NULL);
    if (err != 0) {
        logger->warn("Unable to constain matching to {}-{}", startpos, stoppos);
    } else {
        data->range_begin = startpos;
        data->range_end = stoppos;
    }
}

double
query_pt::match(search::result_vector &family, const cseq& queryc,
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
    cseq_comparator cmp(CMP_IUPAC_OPTIMISTIC, CMP_DIST_NONE, CMP_COVER_QUERY, false);

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
    std::lock_guard<std::mutex> lock(data->arb_pt_access);
    
    int err = aisc_put(data->link,
                       PT_FAMILYFINDER, data->ffinder,
                       FAMILYFINDER_FIND_FAMILY, &bs,
                       NULL);
    if (err != 0) {
        logger->error("Unable to execute find_family command on pt-server");
        if (--maxfail==0) {
            logger->error("No retries left; aborting.");
            return 0;
        }
        logger->error("Retrying...");

        //FIXME restart();
        goto match_retry;
    }

    err = aisc_get(data->link,
                   PT_FAMILYFINDER, data->ffinder,
                   FAMILYFINDER_FAMILY_LIST, f_list.as_result_param(),
                   NULL);
    if (err != 0) {
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
        if (err != 0) {
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
            try {
                const cseq& seq = arb->getCseq(f_name);

                if (max_score <= 2 && cmp(queryc, seq) > max_score) {
                    skipped_max_score ++;
                } else if ((long)seq.size() < min_len) {
                    skipped_min_len ++;
                } else if (noid && boost::algorithm::icontains(seq.getBases(), query)) {
                    skipped_noid ++;
                } else {
                    matches ++;

                    family.emplace_back(f_relscore, &seq);

                    if ((num_full != 0) && (long)seq.size() > full_min_len) {
                        num_full--;
                    }
                    if ((range_cover_right != 0) &&
                        seq.getById(seq.size()-1).getPosition() >= data->range_end) {
                        range_cover_right--;
                    }
                    if ((range_cover_left != 0) &&
                        seq.begin()->getPosition() <= data->range_begin) {
                        range_cover_left--;
                    }
                }
            } catch (base_iupac::bad_character_exception& e) {
                logger->error("Sequence {} contains invalid character{}. Skipping",
                              f_name, e.character);
                skipped_broken++;
            }
        }

        free(f_name);
    } while (matches < max_match
             && (matches <= min_match || f_relscore >= min_score)
             && f_list.exists());

    // get full length sequence
    while (f_list.exists() && ((num_full + range_cover_right + range_cover_left) != 0)) {
        err = aisc_get(data->link, PT_FAMILYLIST, f_list,
                       FAMILYLIST_NAME, &f_name,
                       FAMILYLIST_MATCHES, &f_relscore,
                       FAMILYLIST_NEXT, f_list.as_result_param(),
                       NULL);

        if (err != 0) {
            logger->warn("Unable to get next item in family list");
            break;
        }
        if (data->find_type_fast) {
            f_relscore *= 4;
        }

        const cseq& seq = arb->getCseq(f_name);
        if (max_score >= 2 || cmp(queryc, seq) <= max_score) {
            bool keep = false;
            if ((num_full != 0) && (long)seq.size() > full_min_len) {
                num_full--;
                keep = true;
            }

            if ((range_cover_right != 0) && (long)seq.size() > min_len &&
                seq.getById(seq.size()-1).getPosition() >= data->range_end) {
                range_cover_right--;
                keep = true;
            }

            if ((range_cover_left != 0) && (long)seq.size() > min_len &&
                seq.begin()->getPosition() <= data->range_begin) {
                range_cover_left--;
                keep = true;
            }
            if (keep) {
                family.emplace_back(f_relscore, &seq);
            }
        }
        free(f_name);
    }

    if ((skipped_max_score != 0) || (skipped_broken != 0) || (skipped_min_len != 0) || (skipped_noid != 0)) {
        logger->warn("Skipped {} sequences ({}x id > {}, {}x broken, {}x len < {}, {}x noid)",
                     skipped_max_score + skipped_broken + skipped_min_len + skipped_noid,
                     skipped_max_score, max_score,
                     skipped_broken,
                     skipped_min_len, min_len,
                     skipped_noid);
    }

    return f_relscore;
}

void
query_pt::find(const cseq& query, search::result_vector& results, unsigned int max) {
    data->timeit.start();
    char *error = nullptr;
    results.clear();
    results.reserve(max);

    std::lock_guard<std::mutex> lock(data->arb_pt_access);
    data->timeit.stop("acquire lock");

    bytestring bs;
    bs.data = strdup(query.getBases().c_str());
    bs.size = query.size()+1;
    if (aisc_put(data->link,
                 PT_FAMILYFINDER, data->ffinder,
                 FAMILYFINDER_SORT_MAX, long(max),
                 FAMILYFINDER_FIND_FAMILY, &bs,
                 NULL) != 0) {
        logger->error("Unable to execute find_family command on pt-server");
    }
    free(bs.data);
    data->timeit.stop("send query");

    T_PT_FAMILYLIST f_list;
    aisc_get(data->link,
             PT_FAMILYFINDER, data->ffinder,
             FAMILYFINDER_FAMILY_LIST, f_list.as_result_param(),
             FAMILYFINDER_ERROR, &error,
             NULL);
    if (error && error[0] != 0) {
        logger->error("Unable to get results for search: {}", error);
        return;
    }
    data->timeit.stop("get first");

    char* f_name;
    double f_rel_matches = 0.f;
    long f_matches = 0;
    int mult = (data->find_type_fast) ? 4:1;

    std::vector<std::pair<float, string> > scored_names;
    while (max-- && f_list.exists()) {
        aisc_get(data->link, PT_FAMILYLIST, f_list,
                 FAMILYLIST_NAME, &f_name,
                 FAMILYLIST_REL_MATCHES, &f_rel_matches,
                 FAMILYLIST_MATCHES, &f_matches,
                 FAMILYLIST_NEXT, f_list.as_result_param(),
                 NULL);
        f_rel_matches = 1 - log((f_rel_matches*mult) + 1.0/bs.size)/log(1.0/bs.size);
        scored_names.emplace_back(f_rel_matches, f_name);
        free(f_name);
    }
    data->timeit.stop("get all");

    for (auto& res : scored_names) {
        results.emplace_back(res.first, &data->arbdb->getCseq(res.second));
    }

    data->timeit.stop("load seqs");
}

unsigned int query_pt::size() const {
    return data->arbdb->getSeqCount();
}

struct query_pt_pool::pimpl {
    pimpl(fs::path& filename, int k, bool fast, bool norel, int km, string portname)
        : _filename(filename), _k(k), _fast(fast), _norel(norel), _km(km), _portname(portname)
    {}
    ~pimpl() {
        for (auto pt : _pts) {
            delete pt;
        }
    }

    using pool_map = std::map<fs::path, std::shared_ptr<query_pt_pool::pimpl>>;
    static pool_map _pools;

    // parameters for creating new query_pts
    fs::path _filename;
    int _k;
    bool _fast;
    bool _norel;
    int _km;
    string _portname;

    std::mutex _access_pts;
    std::list<query_pt*> _pts;
    int _count{0};
    query_pt* borrow() {
        int n = 0;
        {
            std::lock_guard<std::mutex> lock(_access_pts);
            if (!_pts.empty()) {
                query_pt* pt = _pts.front();
                _pts.pop_front();
                return pt;
            }
            n = _count++;
        }
        std::string portname;
        if (n > 0) {
            // FIXME: manage the port better. This works for unix sockets, but not
            // for TCP ports.
            portname = fmt::format("{}_{}", _portname, n);
        } else {
            portname = _portname;
        }
        query_pt *pt = query_pt::get_pt_search(
            _filename, _k, _fast, _norel, _km, portname
            );
        return pt;
    }
    void giveback(query_pt* pt) {
        std::lock_guard<std::mutex> lock(_access_pts);
        _pts.push_front(pt);
    }
};

query_pt_pool::pimpl::pool_map query_pt_pool::pimpl::_pools;

query_pt_pool* query_pt_pool::get_pool(fs::path filename,
                                       int k, bool fast, bool norel, int mk,
                                       std::string portname) {
    if (query_pt_pool::pimpl::_pools.count(filename) == 0u) {
        query_pt_pool::pimpl::_pools.emplace(
            filename, std::make_shared<query_pt_pool::pimpl>(filename, k, fast, norel, mk, portname)
            );
    }
    return new query_pt_pool(pimpl::_pools.at(filename));
}

query_pt_pool::query_pt_pool(std::shared_ptr<query_pt_pool::pimpl> p)
    : impl(p) {
}

query_pt_pool::~query_pt_pool() {}

void
query_pt_pool::find(const cseq& query, result_vector& results, unsigned int max) {
    query_pt *pt = impl->borrow();
    pt->find(query, results, max);
    impl->giveback(pt);
}

double
query_pt_pool::match(result_vector& family, const cseq& queryc, int min_match,
                     int max_match, float min_score, float max_score, query_arb *arb,
                     bool noid, int min_len, int num_full, int full_min_len,
                     int range_cover, bool leave_query_out) {
    query_pt *pt = impl->borrow();
    double res = pt->match(family, queryc, min_match, max_match, min_score, max_score, arb,
                           noid, min_len, num_full, full_min_len, range_cover, leave_query_out);
    impl->giveback(pt);
    return res;
}

unsigned int
query_pt_pool::size() const {
    query_pt *pt = impl->borrow();
    unsigned int res = pt->size();
    impl->giveback(pt);
    return res;
}


query_pt_exception::query_pt_exception(std::string msg) noexcept
    : message(std::move(msg))
{
}

query_pt_exception::~query_pt_exception() noexcept = default;

const char*
query_pt_exception::what() const noexcept {
    return message.c_str();
}

} // namespace sina


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
