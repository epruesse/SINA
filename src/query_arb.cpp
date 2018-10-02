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

#include <sstream>
using std::stringstream;

#include <map>
using std::map;
using std::pair;

#include <unordered_map>
using std::unordered_map;

#include <stdio.h>

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

#ifndef HAVE_GBT_FIND_SEQUENCE
inline GBDATA* GBT_find_sequence(GBDATA* gbd, const char* ali) {
    return GBT_read_sequence(gbd, ali);
}
#endif

#include <boost/lambda/lambda.hpp>
using boost::lambda::_1;

#include <boost/lambda/bind.hpp>
using boost::lambda::bind;

#include <boost/thread/mutex.hpp>
#include <boost/functional/hash/hash.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/progress.hpp>
#include <boost/foreach.hpp>

using namespace sina;

// const fieldnames for arb export
const char* query_arb::fn_turn       = "turn";
const char* query_arb::fn_acc        = "acc";
const char* query_arb::fn_start      = "start";
const char* query_arb::fn_qual       = "align_quality_slv";
const char* query_arb::fn_head       = "align_cutoff_head_slv";
const char* query_arb::fn_tail       = "align_cutoff_tail_slv";
const char* query_arb::fn_date       = "aligned_slv";
const char* query_arb::fn_astart     = "align_startpos_slv";
const char* query_arb::fn_astop      = "align_stoppos_slv";
const char* query_arb::fn_idty       = "align_ident_slv";
const char* query_arb::fn_nuc        = "nuc";
const char* query_arb::fn_nuc_gene   = "nuc_gene_slv";
const char* query_arb::fn_bpscore    = "align_bp_score_slv";
const char* query_arb::fn_family_str = "align_family_slv";
const char* query_arb::fn_used_rels  = "used_rels";
const char* query_arb::fn_family     = "NONE";
const char* query_arb::fn_fullname   = "full_name";


// Global lock -- ARB database access is not guaranteed to be thread safe.
// Note: if this ever changes and lock is moved into query ARB, locking
//       order must be considered when accessing two databases at once
//       to prevent chance of deadlock.
static boost::mutex arb_db_access;

// List of opened ARB databases
map<string, query_arb*> query_arb::open_arb_dbs;

static auto logger = Log::create_logger("ARB I/O");

struct query_arb::priv_data {
    priv_data() :
        have_cache(false),
        default_alignment(0L),
        alignment_length(0),
        gbmain(0L),
        gblast(0L),
        gbspec(0L),
        count(0)
    {}

    ~priv_data() {
        if (default_alignment) {
            free(const_cast<char*>(default_alignment));
        }
        if (gbmain) {
            logger->info("Closing ARB database '{}'", filename);
            GB_close(gbmain);
        }
    }

    static GB_shell* the_arb_shell;

    typedef unordered_map<std::string, cseq> sequence_cache_type;
    typedef list<std::string> error_list_type;
    typedef std::unordered_map<string, GBDATA*,
                               boost::hash<string>
                               > gbdata_cache_type;

    sequence_cache_type sequence_cache;
    gbdata_cache_type gbdata_cache;
    error_list_type write_errors;
    bool have_cache;
    const char* default_alignment;
    int alignment_length;
    string filename;
    GBDATA *gbmain, *gblast, *gbspec;
    int count;

    GBDATA* getGBDATA(string name);
    void loadCache();
    void storeCache();
    string getSequence(const char* name, const char* ali);
    void get_weights_nolock(vector<float>& res,
                            const char *name, const char *ali) ;
};

GB_shell *query_arb::priv_data::the_arb_shell = NULL;

GBDATA*
query_arb::priv_data::getGBDATA(string name) {
    if (have_cache) {
        return gbdata_cache[name];
    }

    gbdata_cache_type::iterator it = gbdata_cache.find(name);
    if (it != gbdata_cache.end()) {
        return it->second;
    }

    GBDATA* gbd = GBT_find_species(gbmain, name.c_str());
    gbdata_cache[name] =  gbd;

    return gbd;
}

void
query_arb::priv_data::loadCache() {
    // implementation removed; serialize was too slow to be useful
}

void
query_arb::priv_data::storeCache() {
    // implementation removed; serialize was too slow to be useful
}

string
query_arb::priv_data::getSequence(const char *name, const char *ali) {
    // if there is a preloaded cache, just hand out sequence
    if (have_cache) {
        return sequence_cache[name].getAligned();
    }

    boost::mutex::scoped_lock lock(arb_db_access);
    GB_transaction trans(gbmain);

    if (not ali) {
        ali = default_alignment;
    }

    // get sequence root entry ("species")
    GBDATA *gbdata;
    const char *res;

    if ((gbdata = getGBDATA(name)) &&
        (gbdata = GBT_find_sequence(gbdata, ali)) &&
        (res = GB_read_char_pntr(gbdata))) {
        return res;
    } else {
        return "";
    }
}

void
query_arb::init(const char *db_server) {
    if (not db_server) {
        throw std::runtime_error(string("NULL passed to query_arb::init!"));
    }

    data.gbmain = GB_open(db_server, "rwc");
    if (not data.gbmain) {
        throw std::runtime_error(string("Unable to open ARB database \"") +
                                 db_server + "\".");
    }

    setProtectionLevel(6); // drop privileges

    GB_transaction trans(data.gbmain);

    data.default_alignment = GBT_get_default_alignment(data.gbmain);

    if (not data.default_alignment) {
        GBT_create_alignment(data.gbmain, "ali_16s", 2000, 0, 4, "rna");
        GBT_set_default_alignment(data.gbmain, "ali_16s");
        data.default_alignment=strdup("ali_16s");
        logger->warn("Created new alignment ali_16s in '{}'", data.filename);
    }

    data.alignment_length =  GBT_get_alignment_len(data.gbmain,
                                                   data.default_alignment);
    if (data.alignment_length < 0) {
        throw std::runtime_error(string("Width of default alignment \"") +
                                 data.default_alignment  +
                                 "\" is smaller than zero.");
    }

    data.gbspec = GB_search(data.gbmain, "species_data", GB_CREATE_CONTAINER);
    data.filename = db_server;

    int spec_count = 0;
    for ( GBDATA *gbspec = GBT_first_species(data.gbmain);
          gbspec; gbspec = GBT_next_species(gbspec)) {
        spec_count++;
    }
    data.count=spec_count;

    logger->info("Loading names map... (for {})", data.filename);
    boost::progress_display p(spec_count, std::cerr);
    for ( GBDATA* gbspec = GBT_first_species(data.gbmain);
          gbspec; gbspec = GBT_next_species(gbspec)) {
        data.gbdata_cache[GBT_read_name(gbspec)] = gbspec;
        ++p;
    }
}


query_arb::query_arb(std::string arbfile)
    : data(* new(priv_data))
{
    init(arbfile.c_str());
}

query_arb::~query_arb() {
    delete &data;
}


void
query_arb::setProtectionLevel(int p) {
    GB_transaction trans(data.gbmain);
    GB_change_my_security(data.gbmain, p);
}


void
query_arb::closeOpenARBDBs() {
    for (std::map<string, query_arb*>::iterator it = open_arb_dbs.begin();
         it != open_arb_dbs.end(); ++it) {

        if(it->second->hasErrors()){
            it->second->printErrors(std::cerr);
        }

        delete it->second;
    }
}

query_arb*
query_arb::getARBDB(std::string file_name) {
    boost::mutex::scoped_lock lock(arb_db_access);
    if (not query_arb::priv_data::the_arb_shell) {
        query_arb::priv_data::the_arb_shell = new GB_shell();
    }
    if (open_arb_dbs.size() == 0) {
        atexit(query_arb::closeOpenARBDBs);
    }
    if (not open_arb_dbs.count(file_name)) {
        open_arb_dbs[file_name] = new query_arb(file_name);
    }
    return open_arb_dbs[file_name];
}


void
query_arb::save() {
    saveAs(data.filename.c_str());
}

std::string
query_arb::getFileName() const {
    return data.filename;
}

void
query_arb::saveAs(const char* fname, const char* type) {
    logger->info("Saving database {}", fname);
    {
        GB_transaction trans(data.gbmain);
        logger->info("Checking alignment...");
        GB_ERROR err = GBT_check_data(data.gbmain, data.default_alignment);
        if (err) {
            logger->error(err);
        }
    }

    if (GB_ERROR err = GB_save_as(data.gbmain, fname, type)) {
        logger->error("Error '{}' while trying to save {}", err, fname);
    }
}

bool
query_arb::good() const {
    return data.gbmain && data.default_alignment;
}

static void
loadKey(cseq& c, const string key, GBDATA* gbspec) {
    GBDATA *gbd = GB_find(gbspec, key.c_str(), SEARCH_CHILD);
    if (gbd) {
        switch(GB_read_type(gbd)) {
            case GB_STRING:
                c.set_attr(key, (const char*)GB_read_pntr(gbd));
                return;
            case GB_BYTE:
                c.set_attr(key, (const char)GB_read_byte(gbd));
                return;
            case GB_INT:
                c.set_attr(key, (const int)GB_read_int(gbd));
                return;
            case GB_FLOAT:
                c.set_attr(key, (const float)GB_read_float(gbd));
                return;
            case GB_BITS:
            default:
                logger->error("loadKey failed: type unsupported");
                return;
        }
    }
}

void
query_arb::loadKey(cseq& c, const string key) {
    boost::mutex::scoped_lock lock(arb_db_access);
    GB_transaction trans(data.gbmain);
    ::loadKey(c, key, data.getGBDATA(c.getName()));
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
        if (gbd && GB_read_type(gbd) != GB_INT) {
            GB_delete(gbd);
            gbd = 0;
        }

        if (not gbd) {
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
        if (gbd && GB_read_type(gbd) != GB_FLOAT) {
            GB_delete(gbd);
            gbd = 0;
        }
        if (not gbd) {
            gbd = GB_create(_gbspec, _key.c_str(), GB_FLOAT);
        }

        _query_arb.write(gbd, f);
    }
    void operator()(const string& str) {
        GBT_add_new_changekey(_gbmain, _key.c_str(), GB_STRING);
        GBDATA *gbd = GB_entry(_gbspec, _key.c_str());
        if (gbd && GB_read_type(gbd) != GB_STRING) {
            GB_delete(gbd);
            gbd = 0;
        }
        if (not gbd) {
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
query_arb::storeKey(cseq& c, const std::string key) {
    boost::mutex::scoped_lock lock(arb_db_access);
    GB_transaction trans(data.gbmain);
    storeKey(data.gbmain, data.getGBDATA(c.getName()), key,
               c.get_attr<cseq::variant>(key));
}

void
query_arb::loadCache(std::vector<std::string>& keys) {
    GBDATA *gbspec;

    const char *ali = data.default_alignment;

    boost::mutex::scoped_lock lock(arb_db_access);
    GB_transaction trans(data.gbmain);
    GB_set_cache_size(data.gbmain, 0);

    data.loadCache();
    unsigned int scache_size = data.sequence_cache.size();

    logger->info("Loading {} sequences...", data.count);
    boost::progress_display p(data.count, std::cerr);

    data.sequence_cache.reserve(data.count);

#ifndef HAVE_TBB // serial implementation
    timer t;
    for ( gbspec = GBT_first_species(data.gbmain);
          gbspec;
          gbspec = GBT_next_species(gbspec)) {
        const char* name_str = GBT_read_name(gbspec);
        string name(name_str);
        t.stop("name");
        cseq sequence;
        if (not data.sequence_cache.count(name)) {
            t.start();
            GBDATA *gb_data = GBT_find_sequence(gbspec,ali);
            if (not gb_data) continue;
            const char* ptr =  GB_read_char_pntr(gb_data);
            t.stop("arb load");
            sequence = cseq(name.c_str(), 0.f, ptr);
            t.stop("cseq convert");
        } else {
            sequence = data.sequence_cache[name];
        }

        for(vector<string>::iterator it = keys.begin();
            it != keys.end(); ++it) {
            if (not sequence.get_attrs().count(*it)) {
                ::loadKey(sequence, *it, gbspec);
            }
        }

        data.sequence_cache[sequence.getName()] = sequence;
        ++p;
    }
    logger.info("Timings for Cache Load: {}", t);

#else // HAVE_TBB -- parallel implementation

    gbspec = GBT_first_species(data.gbmain);
    if (not gbspec) {
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
                GBDATA *gb_sequence = NULL;
                while (gbspec && not gb_sequence) {
                    name_str = GBT_read_name(gbspec);
                    gb_sequence = GBT_find_sequence(gbspec, ali);
                    gbspec = GBT_next_species(gbspec);
                }
                arb.stop("arb find");
                if (not gbspec && not gb_sequence) {
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
                for(vector<string>::iterator it = keys.begin();
                    it != keys.end(); ++it) {
                    if (not sequence.get_attrs().count(*it)) {
                        ::loadKey(sequence, *it, gbspec);
                    }
                }
                data.sequence_cache[sequence.getName()] = sequence;
                ++p;
                tkeys.stop("keys");
            })
        );
    all.stop("all");
    logger->info("Timings for cache load: {} {} {}", all, tkeys, arb);
#endif // have TBB

    logger->info("Loaded {} sequences", data.sequence_cache.size());
    if (data.sequence_cache.size() > scache_size)
        data.storeCache();

    data.have_cache = true;
}

vector<cseq*>
query_arb::getCacheContents() {
    vector<cseq*> tmp;
    tmp.reserve(data.sequence_cache.size());
    for (priv_data::sequence_cache_type::iterator it = data.sequence_cache.begin();
         it != data.sequence_cache.end(); ++it) {
        tmp.push_back(&it->second);
    }
    return tmp;
}

long
query_arb::getAlignmentWidth() {
    return data.alignment_length;
}

vector<string>
query_arb::getSequenceNames() {
    vector<string> tmp;
    tmp.reserve(data.gbdata_cache.size());
    for (priv_data::gbdata_cache_type::iterator it = data.gbdata_cache.begin();
         it != data.gbdata_cache.end(); ++it) {
        tmp.push_back(it->first);
    }
    return tmp;
}

cseq&
query_arb::getCseq(string name) {
    // if there is a preloaded cache, just hand out sequence
    if (data.have_cache)
        return data.sequence_cache[name];

    // if not, check whether we already loaded it
    priv_data::sequence_cache_type::iterator it = data.sequence_cache.find(name);
    if (it != data.sequence_cache.end())
        return it->second;

    // if all fails, fetch sequence from arb and cache
    cseq tmp(name.c_str(), 0.f, 
             data.getSequence(name.c_str(), data.default_alignment).c_str()
        );
    data.sequence_cache[name]=tmp;
    return data.sequence_cache[name];
}

int
query_arb::getSeqCount() const {
    return data.count;
}

void
query_arb::putCseq(const cseq& seq) {
    putSequence(seq);

    boost::mutex::scoped_lock lock(arb_db_access);
    GB_transaction trans(data.gbmain);
    pair<string,cseq::variant> ap;
    GBDATA* gbspec = data.getGBDATA(seq.getName());
    BOOST_FOREACH(ap, seq.get_attrs()) {
        storeKey(data.gbmain, gbspec, ap.first, ap.second);
    }
}

void
query_arb::putSequence(const cseq& seq) {
    boost::mutex::scoped_lock lock(arb_db_access);
    GB_transaction trans(data.gbmain);

    const char *ali = data.default_alignment;
    GBDATA *gbdata;

    string aseq_string = seq.getAligned();
    const char *aseq = aseq_string.c_str();

    gbdata = data.getGBDATA(seq.getName());
    if (not gbdata) {
        gbdata =  GB_create_container(data.gbspec, "species");
        GBDATA *gbname  = GB_create(gbdata, "name", GB_STRING);

        if (seq.getName().empty()) {
            stringstream tmp;
            tmp << "slv_" << data.count;
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
            << GB_checksum(aseq, strlen(aseq), true, ".-");
        write(gbacc, tmp.str().c_str());
        data.gbdata_cache[seq.getName()]=gbdata;
    }
    data.gblast = gbdata;

    GBDATA *gbseq = GBT_find_sequence(gbdata, ali);
    if (not gbseq) {
        gbseq = GBT_create_sequence_data(gbdata, ali, "data", GB_STRING, 0);
    }
    write(gbseq, aseq);
}

void
query_arb::copySequence(query_arb& other, std::string name, bool mark) {
    boost::mutex::scoped_lock lock(arb_db_access);

    // lock underlying arb database
    GB_transaction t1(data.gbmain);
    GB_transaction t2(other.data.gbmain);

    GBDATA *gbdest;

    // don't copy if sequence with identical name exists
    if ( (gbdest=GBT_find_species(data.gbmain, name.c_str())) ) {
        logger->error("Species \"{}\" already in target db. Not copying.", name);
        if (mark) write_flag(gbdest,1l);
        return;
    }

    GBDATA *gbsource = other.data.getGBDATA(name);
    gbdest = GB_create_container(data.gbspec, "species");
    if (gbsource && gbdest) {
        GB_copy(gbdest,gbsource);
        logger->info("Copied species {}", name);
        data.gblast = gbdest;
        if (mark) write_flag(gbdest,1l);
        return;
    } else {
        logger->error("Error while copying species \"{}\".", name);
    }

    data.gblast = 0;
}

void
query_arb::setMark() {
    boost::mutex::scoped_lock lock(arb_db_access);
    if (data.gblast) {
        GB_transaction trans(data.gbmain);
        write_flag(data.gblast,1l);
    }
}

void
query_arb::setMark(const std::string& name) {
    boost::mutex::scoped_lock lock(arb_db_access);
    GB_transaction trans(data.gbmain);
    GBDATA *gbdata = data.getGBDATA(name);
    if (gbdata) {
        write_flag(gbdata,1);
        data.gblast = gbdata;
    } else {
        logger->error("Failed to mark species {} - name not found", name);
        data.gblast = 0;
    }
}


/////// filter stuff

string
query_arb::getFilter(string name) {
    boost::mutex::scoped_lock lock(arb_db_access);
    GB_transaction trans(data.gbmain);

    // structure of most sai entries:
    // name -> <name>
    // <alignment> -> {   // (alignment name begins with 'ali'
    //   data -> <chararray>
    // }

    // find SAI container named <name>
    GBDATA *gbsai = GBT_find_SAI(data.gbmain, name.c_str());
    if (not gbsai)
        return "";

    // descend into alignment tag
    gbsai = GB_find(gbsai, data.default_alignment, SEARCH_CHILD);
    if (not gbsai)
        return "";

    // find data entry
    gbsai = GB_find(gbsai, "data", SEARCH_CHILD);
    if (not gbsai)
        return "";

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
        { "NA", "NC", "NG", "NU", "TRANSITIONS", "TRANSVERSIONS", 0L };
    enum {
        NA = 0, NC = 1, NG = 2, NU = 3, TRNS = 4, TRVRS = 5
    };

    boost::mutex::scoped_lock lock(arb_db_access);
    GB_transaction trans(data.gbmain);

    for (GBDATA *gbsai = GBT_first_SAI(data.gbmain); gbsai; 
         gbsai = GBT_next_SAI(gbsai)) {
        // get sai data for current alignment
        GBDATA *gbname = GB_find(gbsai, "name", SEARCH_CHILD);
        if (not gbname) {
            logger->error("SAI without name? Broken DB!");
            continue;
        }
        string name = string(GB_read_char_pntr(gbname));

        GBDATA *gbali = GB_find(gbsai, data.default_alignment, SEARCH_CHILD);
        if (not gbali) {
            continue; // no data for this alignment
        }
        
        

        // get contents of type field
        GBDATA *gbtype = GB_find(gbali, "_TYPE", SEARCH_CHILD);
        if (not gbtype) {
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
        if (not gbfreq) {
            logger->error("ERROR: SAI '{}' is of type PVP but lacks"
                          "contained 'FREQUENCIES'. Your DB might be corrupted!",
                          name);
            continue;
        }
    
        // load frequencies
        unsigned int *pvp_data[6]; 
        for ( int i=0; pvp_names[i]; i++) {
            GBDATA *gbdata = GB_find(gbfreq, pvp_names[i], SEARCH_CHILD);
            if (not gbdata) {
                logger->error("unable to find PVP data {}", pvp_names[i]);
                continue;
            }
            pvp_data[i] = GB_read_ints(gbdata);
        }

        // make an alignment_stats from the frequencies
        alignment_stats a(name,
                          ntaxa, data.alignment_length, 
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
    boost::mutex::scoped_lock lock(arb_db_access);
    GB_transaction trans(data.gbmain);

    const char *ali = data.default_alignment;
    BI_helix helix;

    pairs.resize(data.alignment_length);

    const char* error = helix.init(data.gbmain, ali);
    if (not error) {
        for (int i=0; i<data.alignment_length; ++i)
            pairs[i]=helix.entry(i).pair_pos;
    } else {
        logger->error("No HELIX filter found in ARB file. "
                      "Disabling secondary structure features.");
        for (int i=0; i<data.alignment_length; ++i) {
            pairs[i]=0;
        }
    }
    return pairs;
}

void
query_arb::write(GBDATA *pData, double value) {
    if (not pData) {
        throw query_arb_exception("GB_write_float pData is null");
    }

    GB_ERROR err = GB_write_float(pData, value);
    if (err) { //GB_write_int returns NULL if no error occurred.
        std::stringstream ss;
        ss << "GB_write_float(value = " << value << ") failed. Reason: " << err;
        addError(ss.str());
    }
}

void
query_arb::write(GBDATA *pData, int value) {
    if (not pData) {
        throw query_arb_exception("GB_write_int pData is null");
    }

    GB_ERROR err = GB_write_int(pData, value);
    if (err) {//GB_write_int returns NULL if no error occurred.
        std::stringstream ss;
        ss << "GB_write_int(value = " << value << ") failed. Reason: " << err;
        addError(ss.str());
    }

}



void
query_arb::write(GBDATA *pData, const char* pValue) {
    if (not pData) {
        throw query_arb_exception("GB_write_string pData is null");
    }

    GB_ERROR err = GB_write_string(pData, pValue);
    if (err) {  //GB_write_int returns NULL if no error occurred.
        std::stringstream ss;
        ss << "GB_write_string(value = " << pValue << ") failed. Reason: " << err;
        addError(ss.str());
    }
}

void
query_arb::write_flag(GBDATA *pData, long value) {
    if(not pData) {
        throw query_arb_exception("GB_write_flag pData is null");
    }

    // GB_write_flag kills arb if an error occurs.
    // => no error handling possible
    GB_write_flag(pData, value);
}

void
query_arb::addError(const std::string& message) {
    data.write_errors.push_back(message);
}

bool
sina::query_arb::hasErrors() const {
    return data.write_errors.size() > 0;
}

void
query_arb::printErrors(std::ostream& stream){
    if(hasErrors()) {
        stream << "Following errors occurred while querying arb:" << std::endl;
        BOOST_FOREACH(std::string msg,data.write_errors) {
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
