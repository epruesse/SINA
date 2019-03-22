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

/* this file contains the code that communicates with arb
 * => get data from db, query pt-server */

#include "config.h"

#include "query_arb.h"
#include "timer.h"
#include "alignment_stats.h"
#include "log.h"
#include "progress.h"

#include <iostream>
using std::uppercase;
using std::hex;

#include <string>
using std::string;

#include <vector>
using std::vector;

#include <list>
using std::list;

#include <fstream>
using std::ifstream;
using std::ofstream;

#include <cmath>
#include <exception>
#include <mutex>

#include <sstream>
using std::stringstream;

#include <map>
using std::map;
using std::pair;

#include <unordered_map>
#include <tbb/concurrent_unordered_map.h>

#include <cstdio>

using std::runtime_error;

#ifdef HAVE_TBB
#  include "tbb/pipeline.h"
#endif

// ARB needs either DEBUG or NDEBUG defined
#ifdef DEBUG
#  define typeof __typeof__
#else
# ifndef NDEBUG
#  define NDEBUG 1
# endif
#endif

#include <arbdbt.h>
#include <BI_helix.hxx>
#include <arb_handlers.h>

#ifndef HAVE_GBT_FIND_SEQUENCE
inline GBDATA* GBT_find_sequence(GBDATA* gbd, const char* ali) {
    return GBT_read_sequence(gbd, ali);
}
#endif

#include <boost/functional/hash/hash.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;

using namespace sina;

// const fieldnames for arb export
const char* query_arb::fn_turn       = "turn";
const char* query_arb::fn_acc        = "acc";
const char* query_arb::fn_start      = "start";
const char* query_arb::fn_used_rels  = "used_rels";
const char* query_arb::fn_fullname   = "full_name";
const char* query_arb::fn_nuc        = "nuc";

const char* query_arb::fn_qual       = "align_quality_slv";
const char* query_arb::fn_head       = "align_cutoff_head_slv";
const char* query_arb::fn_tail       = "align_cutoff_tail_slv";
const char* query_arb::fn_date       = "aligned_slv";
const char* query_arb::fn_astart     = "align_startpos_slv";
const char* query_arb::fn_astop      = "align_stoppos_slv";
const char* query_arb::fn_idty       = "align_ident_slv";
const char* query_arb::fn_nuc_gene   = "nuc_gene_slv";
const char* query_arb::fn_bpscore    = "align_bp_score_slv";
const char* query_arb::fn_family_str = "align_family_slv";

const char* query_arb::fn_family     = "NONE";
const char* query_arb::fn_align_log  = "align_log_slv";


// Global lock -- ARB database access is not thread safe! Not even between
//                open datases.
static std::mutex arb_db_access;

// List of opened ARB databases
map<fs::path, query_arb*> query_arb::open_arb_dbs;

static auto arb_logger = Log::create_logger("libARBDB");
static auto logger = Log::create_logger("ARB I/O");

struct query_arb::priv_data {
    ~priv_data() {
        if (default_alignment != nullptr) {
            free(const_cast<char*>(default_alignment));
        }
        if (gbmain != nullptr) {
            logger->warn("Closing ARB database '{}' ...", filename);
            GB_close(gbmain);
        }
    }

    static GB_shell* the_arb_shell;

    using sequence_cache_type = tbb::concurrent_unordered_map<std::string, cseq>;
    using gbdata_cache_type = std::unordered_map<string, GBDATA*, boost::hash<string>>;
    using error_list_type = list<std::string>;

    sequence_cache_type sequence_cache;
    gbdata_cache_type gbdata_cache;
    error_list_type write_errors;
    bool have_cache{false};
    const char* default_alignment{nullptr};
    int alignment_length{0};
    fs::path filename;
    GBDATA *gbmain{nullptr}, *gblast{nullptr}, *gbspec{nullptr};
    int count{0};

    GBDATA* getGBDATA(const string& name);
    string getSequence(const char* name, const char* ali);
    void get_weights_nolock(vector<float>& res,
                            const char *name, const char *ali) ;
};

GB_shell *query_arb::priv_data::the_arb_shell = nullptr;

GBDATA*
query_arb::priv_data::getGBDATA(const string& name) {
    if (have_cache) {
        return gbdata_cache[name];
    }

    auto it = gbdata_cache.find(name);
    if (it != gbdata_cache.end()) {
        return it->second;
    }

    GBDATA* gbd = GBT_find_species(gbmain, name.c_str());
    gbdata_cache[name] =  gbd;

    return gbd;
}

string
query_arb::priv_data::getSequence(const char *name, const char *ali) {
    // if there is a preloaded cache, just hand out sequence
    if (have_cache) {
        return sequence_cache[name].getAligned();
    }

    if (ali == nullptr) {
        ali = default_alignment;
    }

    std::lock_guard<std::mutex> lock(arb_db_access);
    GB_transaction trans(gbmain);

    // get sequence root entry ("species")
    GBDATA *gbdata;
    const char *res;

    if (((gbdata = getGBDATA(name)) != nullptr) &&
        ((gbdata = GBT_find_sequence(gbdata, ali)) != nullptr) &&
        ((res = GB_read_char_pntr(gbdata)) != nullptr)) {
        string out(res);
        GB_flush_cache(gbdata);
        return out;
    }
    return "";
}

query_arb::query_arb(const fs::path& arbfile)
    : data(new priv_data()) {
    data->filename = arbfile;
    if (arbfile.empty()) {
        throw runtime_error("Empty ARB database name?!");
    }

    data->gbmain = GB_open(arbfile.c_str(), "rwc");
    if (data->gbmain == nullptr) {
        throw runtime_error(fmt::format("Unable to open ARB database {}.", arbfile));
    }

    setProtectionLevel(6); // drop privileges

    GB_transaction trans(data->gbmain);

    data->default_alignment = GBT_get_default_alignment(data->gbmain);

    if (data->default_alignment == nullptr) {
        GBT_create_alignment(data->gbmain, "ali_16s", 2000, 0, 4, "rna");
        GBT_set_default_alignment(data->gbmain, "ali_16s");
        data->default_alignment=strdup("ali_16s");
        logger->warn("Created new alignment ali_16s in '{}'", data->filename);
    }

    data->alignment_length =  GBT_get_alignment_len(data->gbmain,
                                                   data->default_alignment);
    if (data->alignment_length < 0) {
        throw runtime_error(
            fmt::format(
                "Width of default alignment \"{}\" in {} is <0",
                data->default_alignment, data->filename
                )
            );
    }

    data->gbspec = GB_search(data->gbmain, "species_data", GB_CREATE_CONTAINER);

    int spec_count = 0;
    for ( GBDATA *gbspec = GBT_first_species(data->gbmain);
          gbspec != nullptr; gbspec = GBT_next_species(gbspec)) {
        spec_count++;
    }
    data->count=spec_count;

    logger->info("Loading names map... (for {})", data->filename);

    logger_progress q(logger, "Scanning", spec_count);

    Progress p("Scanning", spec_count);
    for ( GBDATA* gbspec = GBT_first_species(data->gbmain);
          gbspec != nullptr; gbspec = GBT_next_species(gbspec)) {
        data->gbdata_cache[GBT_read_name(gbspec)] = gbspec;
        ++p;
    }
}

query_arb::~query_arb() {
}


void
query_arb::setProtectionLevel(int p) {
    GB_transaction trans(data->gbmain);
    GB_change_my_security(data->gbmain, p);
}


void
query_arb::closeOpenARBDBs() {
    // atexit registered method
    std::lock_guard<std::mutex> lock(arb_db_access);
    for (auto& it: open_arb_dbs) {
        if(it.second->hasErrors()){
            it.second->printErrors(std::cerr);
        }

        delete it.second;
    }
}


static void log_arb_err(const char *msg) {
    arb_logger->error(msg);
};
static void log_arb_warn(const char *msg) {
    arb_logger->warn(msg);
};
static void log_arb_info(const char *msg) {
    arb_logger->info(msg);
};

logger_progress *arb_progress{nullptr};
static void arb_openstatus(const char* title) {
    if (arb_progress)
        delete arb_progress;
    arb_progress = new logger_progress(logger, title, 100);
}
static void arb_closestatus() {
    delete arb_progress;
    arb_progress = nullptr;
}

static ARB_STATUS_RETURN_TYPE arb_set_title(const char* title) {
    arb_logger->info("Progress title: {}", title);
    return ARB_STATUS_RETURN_VALUE;
}
static ARB_STATUS_RETURN_TYPE arb_set_subtitle(const char* title) {
    arb_logger->info("Progress subttitle: {}", title);
    return ARB_STATUS_RETURN_VALUE;
}
static ARB_STATUS_RETURN_TYPE arb_set_gauge(double gauge) {
    if (arb_progress) {
        auto cur = arb_progress->count();
        auto set_to = arb_progress->size() * gauge;
        arb_progress->update(set_to - cur);
    }
    return ARB_STATUS_RETURN_VALUE;
}
static bool arb_user_abort() {
    return false;
}

static arb_status_implementation log_arb_status {
    AST_RANDOM, arb_openstatus, arb_closestatus, arb_set_title,
        arb_set_subtitle, arb_set_gauge, arb_user_abort
};

static arb_handlers arb_log_handlers = {
    log_arb_err, log_arb_warn, log_arb_info, log_arb_status
};


query_arb*
query_arb::getARBDB(const fs::path& file_name) {
    std::lock_guard<std::mutex> lock(arb_db_access);
    if (query_arb::priv_data::the_arb_shell == nullptr) {
        query_arb::priv_data::the_arb_shell = new GB_shell();

        ARB_install_handlers(arb_log_handlers);
        ARB_redirect_handlers_to(stderr, stderr);
    }
    if (open_arb_dbs.empty()) {
        atexit(query_arb::closeOpenARBDBs);
    }
    if (open_arb_dbs.count(file_name) == 0u) {
        open_arb_dbs[file_name] = new query_arb(file_name);
    }
    return open_arb_dbs[file_name];
}


void
query_arb::save() {
    saveAs(data->filename);
}

const fs::path&
query_arb::getFileName() const {
    return data->filename;
}

void
query_arb::saveAs(const fs::path& fname, const char* type) {
    logger->info("Saving database {}", fname);
    {
        GB_transaction trans(data->gbmain);
        logger->info("Checking alignment...");
        GB_ERROR err = GBT_check_data(data->gbmain, data->default_alignment);
        if (err != nullptr) {
            logger->error(err);
        }
    }

    if (GB_ERROR err = GB_save_as(data->gbmain, fname.c_str(), type)) {
        logger->error("Error '{}' while trying to save {}", err, fname);
    }
}

bool
query_arb::good() const {
    return (data->gbmain != nullptr) && (data->default_alignment != nullptr);
}

static void
loadKey(cseq& c, const string& key, GBDATA* gbspec) {
    GBDATA *gbd = GB_find(gbspec, key.c_str(), SEARCH_CHILD);
    if (gbd != nullptr) {
        switch(GB_read_type(gbd)) {
            case GB_STRING:
                c.set_attr(key, (const char*)GB_read_pntr(gbd));
                return;
            case GB_BYTE:
                c.set_attr(key, (char)GB_read_byte(gbd));
                return;
            case GB_INT:
                c.set_attr(key, (int)GB_read_int(gbd));
                return;
            case GB_FLOAT:
                c.set_attr(key, (float)GB_read_float(gbd));
                return;
            case GB_BITS:
            default:
                logger->error("loadKey failed: type unsupported");
                return;
        }
    }
}

void
query_arb::loadKey(cseq& c, const string& key) {
    std::lock_guard<std::mutex> lock(arb_db_access);
    GB_transaction trans(data->gbmain);
    ::loadKey(c, key, data->getGBDATA(c.getName()));
}

struct query_arb::storeKey_visitor
    : public boost::static_visitor<> {
    /**
     * @param gbmain
     * @param gbspec
     * @param key
     * @param queryArb Pointer to the instance of query_arb. Needed to call the write methods
     */
    storeKey_visitor(GBDATA* gbmain, GBDATA* gbspec, const string& key, query_arb& queryArb)
        : _gbmain(gbmain), _gbspec(gbspec), _key(key), _query_arb(queryArb)  {}

    GBDATA* _gbmain;
    GBDATA* _gbspec;
    const string& _key;
    query_arb& _query_arb;

    void operator()(const int& i) {
        GBT_add_new_changekey(_gbmain, _key.c_str(), GB_INT);
        GBDATA *gbd = GB_entry(_gbspec, _key.c_str());
        if ((gbd != nullptr) && GB_read_type(gbd) != GB_INT) {
            GB_delete(gbd);
            gbd = nullptr;
        }

        if (gbd == nullptr) {
            gbd = GB_create(_gbspec, _key.c_str(), GB_INT);
        }
        _query_arb.write(gbd, i);


    }
    void operator()(const unsigned int& ui) {
        operator()((int)ui);
    }
    void operator()(const bool& b) {
        operator()((int)b);
    }
    void operator()(const float& f) {
        GBT_add_new_changekey(_gbmain, _key.c_str(), GB_FLOAT);
        GBDATA *gbd = GB_entry(_gbspec, _key.c_str());
        if ((gbd != nullptr) && GB_read_type(gbd) != GB_FLOAT) {
            GB_delete(gbd);
            gbd = nullptr;
        }
        if (gbd == nullptr) {
            gbd = GB_create(_gbspec, _key.c_str(), GB_FLOAT);
        }

        _query_arb.write(gbd, f);
    }
    void operator()(const string& str) {
        GBT_add_new_changekey(_gbmain, _key.c_str(), GB_STRING);
        GBDATA *gbd = GB_entry(_gbspec, _key.c_str());
        if ((gbd != nullptr) && GB_read_type(gbd) != GB_STRING) {
            GB_delete(gbd);
            gbd = nullptr;
        }
        if (gbd == nullptr) {
            gbd = GB_create(_gbspec, _key.c_str(), GB_STRING);
        }

        _query_arb.write(gbd, str.c_str());

    }

    void operator()(const vector<cseq>& /*vc*/) {
    }
};

inline void
query_arb::storeKey(GBDATA* gbmain, GBDATA* gbspec, const std::string& key,
                    cseq::variant var) {
    storeKey_visitor vis(gbmain, gbspec, key,*this);
    boost::apply_visitor(vis, var);
}

void
query_arb::storeKey(cseq& c, const std::string& key) {
    std::lock_guard<std::mutex> lock(arb_db_access);
    GB_transaction trans(data->gbmain);
    storeKey(data->gbmain, data->getGBDATA(c.getName()), key,
               c.get_attr<cseq::variant>(key));
}

void
query_arb::loadCache(std::vector<std::string>& keys) {
    GBDATA *gbspec;

    const char *ali = data->default_alignment;

    std::lock_guard<std::mutex> lock(arb_db_access);
    GB_transaction trans(data->gbmain);

    logger->info("Loading {} sequences...", data->count);
    logger_progress p(logger, "Loading sequences", data->count);

    // re-init sequence_cache with size data->count?

#undef HAVE_TBB
#ifndef HAVE_TBB // serial implementation
    timer t;
    for ( gbspec = GBT_first_species(data->gbmain);
          gbspec;
          gbspec = GBT_next_species(gbspec)) {
        t.start();
        const char* name_str = GBT_read_name(gbspec);
        string name(name_str);
        t.stop("name");
        cseq sequence;
        if (not data->sequence_cache.count(name)) {
            GBDATA *gb_data = GBT_find_sequence(gbspec,ali);
            if (not gb_data) continue;
            const char* ptr =  GB_read_char_pntr(gb_data);
            t.stop("arb load");
            sequence = cseq(name.c_str(), 0.f, ptr);
            t.stop("cseq convert");
            GB_flush_cache(gb_data);
        } else {
            sequence = data->sequence_cache[name];
        }

        for(vector<string>::iterator it = keys.begin();
            it != keys.end(); ++it) {
            if (not sequence.get_attrs().count(*it)) {
                ::loadKey(sequence, *it, gbspec);
            }
        }


        data->sequence_cache[sequence.getName()] = sequence;
        ++p;
    }
    logger->info("Timings for Cache Load: {}", t);

#else // HAVE_TBB -- parallel implementation

    gbspec = GBT_first_species(data->gbmain);
    if (gbspec == nullptr) {
        logger->error("Failed to load sequences -- database empty?");
        return;
    }

    timer arb, tkeys, all;
    all.start();
    int n=0;
    tbb::parallel_pipeline(
        /*max_number_of_live_token=*/ 1024,
        tbb::make_filter<void, std::pair<const char*, const char*>>(
            tbb::filter::serial,
            [&](tbb::flow_control& fc) -> std::pair<const char*, const char*>
            {
                arb.start();
                const char* name_str;
                GBDATA *gb_sequence = nullptr;
                while ((gbspec != nullptr) && (gb_sequence == nullptr)) {
                    name_str = GBT_read_name(gbspec);
                    gb_sequence = GBT_find_sequence(gbspec, ali);
                    gbspec = GBT_next_species(gbspec);
                }
                arb.stop("arb find");
                if ((gbspec == nullptr) && (gb_sequence == nullptr)) {
                    fc.stop();
                    return std::make_pair("","");
                }
                const char* seq_str = GB_read_char_pntr(gb_sequence);
                arb.stop("arb load");
                n++;
                return std::make_pair(name_str, seq_str);
            }) &
        tbb::make_filter<std::pair<const char*, const char*>, cseq>(
            tbb::filter::parallel,
            [&](std::pair<const char*, const char*> p) -> cseq
            {
                return cseq(p.first, 0.f, p.second);
            }) &
        tbb::make_filter<cseq, void>(
            tbb::filter::serial,
            [&](cseq sequence) -> void
            {
                tkeys.start();
                for(auto & key : keys) {
                    if (sequence.get_attrs().count(key) == 0u) {
                        ::loadKey(sequence, key, gbspec); //FIXME!!!
                    }
                }
                data->sequence_cache[sequence.getName()] = sequence;
                ++p;
                tkeys.stop("keys");
            })
        );
    all.stop("all");
    logger->info("Timings for cache load: {} {} {}", all, tkeys, arb);
#endif // have TBB

    logger->info("Loaded {} sequences", data->sequence_cache.size());

    data->have_cache = true;
}

vector<cseq*>
query_arb::getCacheContents() {
    vector<cseq*> tmp;
    tmp.reserve(data->sequence_cache.size());
    for (auto & it : data->sequence_cache) {
        tmp.push_back(&it.second);
    }
    return tmp;
}

long
query_arb::getAlignmentWidth() {
    return data->alignment_length;
}

vector<string>
query_arb::getSequenceNames() {
    vector<string> tmp;
    tmp.reserve(data->gbdata_cache.size());
    for (auto & it : data->gbdata_cache) {
        tmp.push_back(it.first);
    }
    return tmp;
}

cseq&
query_arb::getCseq(const string& name) { //, bool nocache) {
    // if there is a preloaded cache, just hand out sequence
    if (data->have_cache) {
        return data->sequence_cache[name];
    }
    
    // if not, check whether we already loaded it
    auto it = data->sequence_cache.find(name);
    if (it != data->sequence_cache.end()) {
        return it->second;
    }

    // if all fails, fetch sequence from arb and cache
    cseq tmp(name.c_str(), 0.f, 
             data->getSequence(name.c_str(), data->default_alignment).c_str()
        );
    data->sequence_cache[name]=tmp;
#if defined(DEBUG)
    int n_cached = data->sequence_cache.size();
    if (n_cached % 1000 == 0) {
        logger->error("Cache size: {}", n_cached);
        long total = 0;
        for (auto& c : data->sequence_cache) {
            total += c.second.getAlignedBases().size();
        }
        logger->error("  {} bases of {} bytes = {} MB",
                      total, sizeof(aligned_base), total*sizeof(aligned_base) /1024/1024);
    }
#endif
    return data->sequence_cache[name];
}

cseq
query_arb::getCseqUncached(const string& name) { //, bool nocache) {
    return {name.c_str(), 0.f,
            data->getSequence(name.c_str(), data->default_alignment).c_str()};
}

int
query_arb::getSeqCount() const {
    return data->count;
}

void
query_arb::putCseq(const cseq& seq) {
    putSequence(seq);

    std::lock_guard<std::mutex> lock(arb_db_access);
    GB_transaction trans(data->gbmain);
    GBDATA* gbspec = data->getGBDATA(seq.getName());
    for (auto& ap: seq.get_attrs()) {
        storeKey(data->gbmain, gbspec, ap.first, ap.second);
    }
}

void
query_arb::putSequence(const cseq& seq) {
    std::lock_guard<std::mutex> lock(arb_db_access);
    GB_transaction trans(data->gbmain);

    const char *ali = data->default_alignment;
    GBDATA *gbdata;

    string aseq_string = seq.getAligned();
    const char *aseq = aseq_string.c_str();

    gbdata = data->getGBDATA(seq.getName());
    if (gbdata == nullptr) {
        gbdata =  GB_create_container(data->gbspec, "species");
        GBDATA *gbname  = GB_create(gbdata, "name", GB_STRING);

        if (seq.getName().empty()) {
            stringstream tmp;
            tmp << "slv_" << data->count;
            write(gbname, tmp.str().c_str());
        } else {
            write(gbname, seq.getName().c_str());
        }

        GBDATA *gbseq = GBT_create_sequence_data(gbdata, ali, "data", GB_STRING, 0);
        write(gbseq, aseq);

        GBDATA *gbfname = GB_create(gbdata, "full_name", GB_STRING);
        write(gbfname, seq.getName().c_str());

        GBDATA *gbacc = GB_create(gbdata, "acc", GB_STRING);
        stringstream tmp;
        tmp.str("");
        tmp << "ARB_" << uppercase << hex
            << GB_checksum(aseq, strlen(aseq), 1, ".-");
        write(gbacc, tmp.str().c_str());
        data->gbdata_cache[seq.getName()]=gbdata;
    }
    data->gblast = gbdata;

    GBDATA *gbseq = GBT_find_sequence(gbdata, ali);
    if (gbseq == nullptr) {
        gbseq = GBT_create_sequence_data(gbdata, ali, "data", GB_STRING, 0);
    }
    write(gbseq, aseq);
}

void
query_arb::copySequence(query_arb& other, const std::string& name, bool mark) {
    std::lock_guard<std::mutex> lock(arb_db_access);

    // lock underlying arb database
    GB_transaction t1(data->gbmain);
    GB_transaction t2(other.data->gbmain);

    GBDATA *gbdest;

    // don't copy if sequence with identical name exists
    if ( (gbdest=GBT_find_species(data->gbmain, name.c_str())) != nullptr ) {
        logger->error("Species \"{}\" already in target db. Not copying.", name);
        if (mark) {
            write_flag(gbdest, 1l);
        }
        return;
    }

    GBDATA *gbsource = other.data->getGBDATA(name);
    gbdest = GB_create_container(data->gbspec, "species");
    if ((gbsource != nullptr) && (gbdest != nullptr)) {
        GB_copy(gbdest,gbsource);
        logger->info("Copied species {}", name);
        data->gblast = gbdest;
        if (mark) {
            write_flag(gbdest, 1l);
        }
        return;
    }
    logger->error("Error while copying species \"{}\".", name);

    data->gblast = nullptr;
}

void
query_arb::setMark() {
    std::lock_guard<std::mutex> lock(arb_db_access);
    if (data->gblast != nullptr) {
        GB_transaction trans(data->gbmain);
        write_flag(data->gblast,1l);
    }
}

void
query_arb::setMark(const std::string& name) {
    std::lock_guard<std::mutex> lock(arb_db_access);
    GB_transaction trans(data->gbmain);
    GBDATA *gbdata = data->getGBDATA(name);
    if (gbdata != nullptr) {
        write_flag(gbdata,1);
        data->gblast = gbdata;
    } else {
        logger->error("Failed to mark species {} - name not found", name);
        data->gblast = nullptr;
    }
}


/////// filter stuff

string
query_arb::getFilter(const string& name) {
    std::lock_guard<std::mutex> lock(arb_db_access);
    GB_transaction trans(data->gbmain);

    // structure of most sai entries:
    // name -> <name>
    // <alignment> -> {   // (alignment name begins with 'ali'
    //   data -> <chararray>
    // }

    // find SAI container named <name>
    GBDATA *gbsai = GBT_find_SAI(data->gbmain, name.c_str());
    if (gbsai == nullptr) {
        return "";
    }

    // descend into alignment tag
    gbsai = GB_find(gbsai, data->default_alignment, SEARCH_CHILD);
    if (gbsai == nullptr) {
        return "";
    }

    // find data entry
    gbsai = GB_find(gbsai, "data", SEARCH_CHILD);
    if (gbsai == nullptr) {
        return "";
    }

    // read data entry and return as string
    return string(GB_read_char_pntr(gbsai));
}

vector<alignment_stats>
query_arb::getAlignmentStats() {
    vector<alignment_stats> res;
    vector<int> pairs = getPairs();
    
    // What ARB calls "Transitions" is really "Mutations"
    // Actual transitions could be computed by subtracting transversions
    const char *pvp_names[] =
        { "NA", "NC", "NG", "NU", "TRANSITIONS", "TRANSVERSIONS", nullptr };
    enum {
        NA = 0, NC = 1, NG = 2, NU = 3, TRNS = 4, TRVRS = 5
    };

    std::lock_guard<std::mutex> lock(arb_db_access);
    GB_transaction trans(data->gbmain);

    for (GBDATA *gbsai = GBT_first_SAI(data->gbmain); gbsai != nullptr;
         gbsai = GBT_next_SAI(gbsai)) {
        // get sai data for current alignment
        GBDATA *gbname = GB_find(gbsai, "name", SEARCH_CHILD);
        if (gbname == nullptr) {
            logger->error("SAI without name? Broken DB!");
            continue;
        }
        string name = string(GB_read_char_pntr(gbname));

        GBDATA *gbali = GB_find(gbsai, data->default_alignment, SEARCH_CHILD);
        if (gbali == nullptr) {
            continue; // no data for this alignment
        }
        
        

        // get contents of type field
        GBDATA *gbtype = GB_find(gbali, "_TYPE", SEARCH_CHILD);
        if (gbtype == nullptr) {
            continue; // no _TYPE field
        }
        string type = string(GB_read_char_pntr(gbtype));

        // check if this is a "positional variability by parsimony" sai
        if (type.substr(0,4) != "PVP:") {
            continue;
        }
        
        // extract number of taxa from the type field 
        // (yeah, weird place, hope the value will keep being 
        // printed there...)
        size_t ntaxa_end = type.rfind("ntaxa ") + 6;
        int ntaxa = atoi(type.substr(ntaxa_end).c_str());
        
        // get frequencies container
        GBDATA *gbfreq = GB_find(gbali, "FREQUENCIES", SEARCH_CHILD);
        if (gbfreq == nullptr) {
            logger->error("ERROR: SAI '{}' is of type PVP but lacks"
                          "contained 'FREQUENCIES'. Your DB might be corrupted!",
                          name);
            continue;
        }
    
        // load frequencies
        unsigned int *pvp_data[6]; 
        for ( int i=0; pvp_names[i] != nullptr; i++) {
            GBDATA *gbdata = GB_find(gbfreq, pvp_names[i], SEARCH_CHILD);
            if (gbdata == nullptr) {
                logger->error("unable to find PVP data {}", pvp_names[i]);
                continue;
            }
            pvp_data[i] = GB_read_ints(gbdata);
        }

        // make an alignment_stats from the frequencies
        alignment_stats a(name,
                          ntaxa, data->alignment_length,
                          pvp_data[NA], pvp_data[NG], pvp_data[NC],
                          pvp_data[NU], pvp_data[TRNS], pvp_data[TRVRS],
                          pairs);
        res.push_back(a);
    }
    return res;
}

vector<int>
query_arb::getPairs() {
    vector<int> pairs;
    std::lock_guard<std::mutex> lock(arb_db_access);
    GB_transaction trans(data->gbmain);

    const char *ali = data->default_alignment;
    BI_helix helix;

    pairs.resize(data->alignment_length);

    const char* error = helix.init(data->gbmain, ali);
    if (error == nullptr) {
        for (int i=0; i<data->alignment_length; ++i) {
            pairs[i]=helix.entry(i).pair_pos;
        }
    } else {
        logger->error("No HELIX filter found in ARB file. "
                      "Disabling secondary structure features.");
        for (int i=0; i<data->alignment_length; ++i) {
            pairs[i]=0;
        }
    }
    return pairs;
}

void
query_arb::write(GBDATA *pData, double value) {
    if (pData == nullptr) {
        throw query_arb_exception("GB_write_float pData is null");
    }

    GB_ERROR err = GB_write_float(pData, value);
    if (err != nullptr) { //GB_write_int returns NULL if no error occurred.
        std::stringstream ss;
        ss << "GB_write_float(value = " << value << ") failed. Reason: " << err;
        addError(ss.str());
    }
}

void
query_arb::write(GBDATA *pData, int value) {
    if (pData == nullptr) {
        throw query_arb_exception("GB_write_int pData is null");
    }

    GB_ERROR err = GB_write_int(pData, value);
    if (err != nullptr) {//GB_write_int returns NULL if no error occurred.
        std::stringstream ss;
        ss << "GB_write_int(value = " << value << ") failed. Reason: " << err;
        addError(ss.str());
    }

}



void
query_arb::write(GBDATA *pData, const char* pValue) {
    if (pData == nullptr) {
        throw query_arb_exception("GB_write_string pData is null");
    }

    GB_ERROR err = GB_write_string(pData, pValue);
    if (err != nullptr) {  //GB_write_int returns NULL if no error occurred.
        std::stringstream ss;
        ss << "GB_write_string(value = " << pValue << ") failed. Reason: " << err;
        addError(ss.str());
    }
}

void
query_arb::write_flag(GBDATA *pData, long value) {
    if(pData == nullptr) {
        throw query_arb_exception("GB_write_flag pData is null");
    }

    // GB_write_flag kills arb if an error occurs.
    // => no error handling possible
    GB_write_flag(pData, value);
}

void
query_arb::addError(const std::string& message) {
    data->write_errors.push_back(message);
}

bool
sina::query_arb::hasErrors() const {
    return !data->write_errors.empty();
}

void
query_arb::printErrors(std::ostream& stream){
    if(hasErrors()) {
        stream << "Following errors occurred while querying arb:" << std::endl;
        for (auto& msg: data->write_errors) {
            stream << msg.substr(0,70) << std::endl;
        }
    }
    else{
        stream << "No errors" << std::endl;
    }
}

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
