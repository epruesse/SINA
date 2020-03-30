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
const char* query_arb::fn_family     = "align_family_slv";
const char* query_arb::fn_align_log  = "align_log_slv";
const char* query_arb::fn_filter     = "align_filter_slv";
const char* query_arb::fn_nearest    = "nearest_slv";

// Global lock -- ARB database access is not thread safe! Not even between
//                open datases.
static std::mutex arb_db_access;

// List of opened ARB databases
map<fs::path, query_arb*> query_arb::open_arb_dbs;

static auto arb_logger = Log::create_logger("libARBDB");
static auto logger = Log::create_logger("ARB I/O");

template<typename... Args>
query_arb_exception make_exception(const char *fmt, const Args&... args) {
    return query_arb_exception(fmt::format(fmt, args...));
}


/////////////////////   Hidden Implementation Class ///////////////

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
    const char* default_alignment{nullptr};
    int alignment_length{0};
    fs::path filename;
    GBDATA *gbmain{nullptr}, *gblast{nullptr}, *gbspec{nullptr};
    int count{0};

    void write(GBDATA* pData, const string& str);
    GBDATA* getGBDATA(const string& name, bool create = false);
    const char* getSequence(GBDATA*);
    void putSequence(const cseq& seq);
    void get_weights_nolock(vector<float>& res,
                            const char *name, const char *ali);
};

GB_shell *query_arb::priv_data::the_arb_shell = nullptr;


void
query_arb::priv_data::write(GBDATA *pData, const string& str) {
    if (pData == nullptr) {
        throw make_exception("GB_write_string pData is null");
    }

    GB_ERROR err = GB_write_string(pData, str.c_str());
    if (err != nullptr) {
        logger->error("GB_write_string(value = {}) failed. Reason: {}",
                      str, err);
    }
}


GBDATA*
query_arb::priv_data::getGBDATA(const string& name, bool create) {
    logger->info("getGBDATA {} {}", name, create);
    auto it = gbdata_cache.find(name);
    if (it != gbdata_cache.end()) {
        return it->second;
    }
    GBDATA* gbd = GBT_find_species(gbmain, name.c_str());
    if (gbd == nullptr) {
        if (!create) {
            throw make_exception("No sequence \"{}\" in {}",
                                 name, filename.filename());
        }
        logger->info("Creating new sequence in {}", filename.filename());
        gbd = GB_create_container(gbspec, "species");
        GBDATA *gbname = GB_create(gbd, "name", GB_STRING);
        ++count;
        if (name.empty()) {
            write(gbname, fmt::format("sina_{}", count));
        } else {
            write(gbname, name);
        }
    }
    gbdata_cache[name] =  gbd;

    return gbd;
}

const char*
query_arb::priv_data::getSequence(GBDATA* gbdata) {
    gbdata = GBT_find_sequence(gbdata, default_alignment);
    if (gbdata != nullptr) {
        const char *res = GB_read_char_pntr(gbdata);
        if (res != nullptr) {
            return res;
        }
    }
    return nullptr;
}

void
query_arb::priv_data::putSequence(const cseq& seq) {
    string aseq = seq.getAligned();

    std::lock_guard<std::mutex> lock(arb_db_access);
    GB_transaction trans(gbmain);
    GBDATA *gbdata = getGBDATA(seq.getName(), true);
    GBDATA *gbseq = GBT_find_sequence(gbdata, default_alignment);
    if (gbseq == nullptr) {
        gbseq = GBT_create_sequence_data(gbdata, default_alignment, "data", GB_STRING, 0);
        GBDATA *gbacc = GB_create(gbdata, "acc", GB_STRING);
        string acc = fmt::format(
            "ARB_{:X}", GB_checksum(aseq.c_str(), aseq.size(), 1, ".-")
            );
        write(gbacc, acc);
    }
    write(gbseq, aseq);
}

/////////////////////  Reading / Writing attributes  ///////////////

/* Loads a metadata attribute from ARB into cseq
 */
void loadKey(cseq& c, const string& key, GBDATA* gbspec) {
    GBDATA *gbd = GB_find(gbspec, key.c_str(), SEARCH_CHILD);
    if (gbd != nullptr) {
        switch(GB_read_type(gbd)) {
            case GB_STRING:
                c.set_attr(key, (const char*)GB_read_pntr(gbd));
                return;
            case GB_INT:
                c.set_attr(key, (int)GB_read_int(gbd));
                return;
            case GB_FLOAT:
                c.set_attr(key, (float)GB_read_float(gbd));
                return;
            case GB_BYTE:
            case GB_BITS:
            default:
                logger->error("loadKey failed: type unsupported");
                return;
        }
    } else {
        logger->error("loadKey failed: sequence not found");
    }
}


struct storeKey_visitor : public boost::static_visitor<> {
    storeKey_visitor(GBDATA* gbmain, GBDATA* gbspec, const string& key)
        : _gbmain(gbmain), _gbspec(gbspec), _key(key)
    {}

    GBDATA* _gbmain;
    GBDATA* _gbspec;
    const string& _key;

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

        GB_ERROR err = GB_write_int(gbd, i);
        if (err != nullptr) {
            logger->error("GB_write_int(,{}) failed. Reason: {}",
                          i, err);
        }
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
        GB_ERROR err = GB_write_float(gbd, f);
        if (err != nullptr) {
            logger->error("GB_write_float(,{}) failed. Reason: {}",
                          f, err);
        }
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
        GB_ERROR err = GB_write_string(gbd, str.c_str());

        if (err != nullptr) {
            logger->error("GB_write_string(,{}) failed. Reason: {}",
                          str, err);
        }
    }

    void operator()(const vector<cseq>& /*vc*/) {
    }
};


/////////////////////   Outside Class   ///////////////

void
query_arb::storeKey(cseq& c, const std::string& key) {
    std::lock_guard<std::mutex> lock(arb_db_access);
    GB_transaction trans(data->gbmain);
    GBDATA *gbspec = data->getGBDATA(c.getName());
    storeKey_visitor vis(data->gbmain, gbspec, key);
    boost::apply_visitor(vis, c.get_attr<cseq::variant>(key));
}

void
query_arb::loadKey(const cseq& c, const string& key, bool reload) {
    if (!reload && c.has_attr(key)) return;
    std::lock_guard<std::mutex> lock(arb_db_access);
    GB_transaction trans(data->gbmain);
    // FIXME: should we check if the cseq is actually ours?
    ::loadKey(const_cast<cseq&>(c), key, data->getGBDATA(c.getName()));
}


query_arb::query_arb(const fs::path& arbfile)
    : data(new priv_data()) {
    data->filename = arbfile;
    if (arbfile.empty()) {
        throw make_exception("Empty ARB database name?!");
    }

    data->gbmain = GB_open(arbfile.c_str(), "rwc");
    if (data->gbmain == nullptr) {
        throw make_exception("Unable to open ARB database {}.", arbfile);
    }

    setProtectionLevel(6); // drop privileges

    GB_transaction trans(data->gbmain);

    GBDATA *gbd;
    if ((gbd = GB_entry(data->gbmain, "ptserver")) != nullptr &&
        (gbd = GB_entry(gbd, "dbstate")) != nullptr &&
        GB_read_int(gbd) > 0) {
        throw make_exception(
            "{} has been compressed for use by the ARB PT server"
            "and cannot be accessed by SINA.",
            arbfile.filename()
            );
    }

    data->default_alignment = GBT_get_default_alignment(data->gbmain);

    if (data->default_alignment == nullptr) {
        GBT_create_alignment(data->gbmain, "ali_16s", 2000, 0, 4, "rna");
        GBT_set_default_alignment(data->gbmain, "ali_16s");
        data->default_alignment=strdup("ali_16s");
        logger->warn("Created new alignment ali_16s in '{}'", data->filename);
    }

    data->alignment_length = GBT_get_alignment_len(data->gbmain,
                                                   data->default_alignment);
    if (data->alignment_length < 0) {
        // This should not actually be possible. LCOV_EXCL_START
        throw make_exception(
            "Width of default alignment \"{}\" in {} is <0 ?!?!",
            data->default_alignment, data->filename
            );
        // LCOV_EXCL_STOP
    }

    data->gbspec = GB_search(data->gbmain, "species_data", GB_CREATE_CONTAINER);

    int spec_count = 0;
    for ( GBDATA *gbspec = GBT_first_species(data->gbmain);
          gbspec != nullptr; gbspec = GBT_next_species(gbspec)) {
        spec_count++;
    }
    data->count = spec_count;

    logger->info("Loading names map... (for {})", data->filename);

    logger_progress p(logger, "Scanning", spec_count);
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



/////////////////////   ARB logging callbacks  ///////////////

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


/////////////////////   Static Instance Acessors  ///////////////

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
query_arb::closeOpenARBDBs() {
    // atexit registered method
    std::lock_guard<std::mutex> lock(arb_db_access);
    for (auto arb : open_arb_dbs) {
        delete arb.second;
    }
    open_arb_dbs.clear();
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
            // LCOV_EXCL_START
            logger->error("Error '{}' while checking ARB database alignment");
            // LCOV_EXCL_STOP
        }
    }

    if (GB_ERROR err = GB_save_as(data->gbmain, fname.c_str(), type)) {
        logger->error("Error while trying to save {}", fname);
        logger->error("  ARB said: \n{}\n", err);
    }
}


void
query_arb::loadCache(std::vector<std::string>& keys) {
    GBDATA *gbspec;

    const char *ali = data->default_alignment;

    std::lock_guard<std::mutex> lock(arb_db_access);
    GB_transaction trans(data->gbmain);

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
            sequence = cseq(name.c_str(), ptr);
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
query_arb::getAlignmentWidth() const {
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

const cseq&
query_arb::getCseq(const string& name) { //, bool nocache) {
    // if not, check whether we already loaded it
    auto it = data->sequence_cache.find(name);
    if (it != data->sequence_cache.end()) {
        return it->second;
    }

    // not checked, load cache and return item

    auto res = data->sequence_cache.emplace(
        name, std::move(getCseqUncached(name)));
    return res.first->second;
}

cseq
query_arb::getCseqUncached(const string& name) { //, bool nocache) {
    std::lock_guard<std::mutex> lock(arb_db_access);
    GB_transaction trans(data->gbmain);
    GBDATA *spec = data->getGBDATA(name);
    const char* seq = data->getSequence(spec);
    if (seq == nullptr) {
        throw make_exception("No alignment for sequence \"{}\" in {}",
                             name, data->filename.filename());
    }
    cseq c = {name.c_str(), seq};
    GB_flush_cache(spec);
    return c;
}

int
query_arb::getSeqCount() const {
    return data->count;
}

void
query_arb::putCseq(const cseq& seq) {
    data->putSequence(seq);

    std::lock_guard<std::mutex> lock(arb_db_access);
    GB_transaction trans(data->gbmain);
    GBDATA* gbspec = data->getGBDATA(seq.getName());
    for (auto& ap: seq.get_attrs()) {
        storeKey_visitor vis(data->gbmain, gbspec, ap.first);
        boost::apply_visitor(vis, ap.second);
    }
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
            GB_write_flag(gbdest, 1);
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
            GB_write_flag(gbdest, 1);
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
        GB_write_flag(data->gblast, 1);
    }
}

void
query_arb::setMark(const std::string& name) {
    std::lock_guard<std::mutex> lock(arb_db_access);
    GB_transaction trans(data->gbmain);
    GBDATA *gbdata = data->getGBDATA(name);
    if (gbdata != nullptr) {
        GB_write_flag(gbdata, 1);
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
