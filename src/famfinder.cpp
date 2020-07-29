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

#include "famfinder.h"
#include "config.h"
#include "query_pt.h"
#include "query_arb.h"
#include "kmer_search.h"
#include "log.h"
#include "cseq_comparator.h"

#include <iostream>
using std::ostream;

#include <sstream>
using std::stringstream;

#include <iomanip>
using std::setprecision;

#include <string>
using std::string;

#include <vector>
using std::vector;

using std::pair;

#include <exception>
using std::logic_error;


#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

#include <boost/algorithm/string/predicate.hpp>
using boost::algorithm::istarts_with;
using boost::algorithm::iequals;



namespace sina {

static auto logger = Log::create_logger("famfinder");



struct options {
    TURN_TYPE turn_which;
    ENGINE_TYPE engine;

    int gene_start;
    int gene_end;

    string posvar_filter;
    string posvar_autofilter_field;
    float  posvar_autofilter_thres;

    unsigned int   fs_min;
    unsigned int   fs_max;
    float          fs_msc;
    float          fs_msc_max;
    bool           fs_leave_query_out;
    unsigned int   fs_req;
    unsigned int   fs_req_full;
    unsigned int   fs_full_len;
    unsigned int   fs_req_gaps;
    bool           fs_no_fast;
    unsigned int   fs_kmer_len;
    unsigned int   fs_kmer_mm;
    bool           fs_kmer_norel;
    unsigned int   fs_min_len;
    unsigned int   fs_cover_gene;

    fs::path database;
    string   pt_port;

    bool oldmatch;
};

static options opts;



void validate(boost::any& v,
              const std::vector<std::string>& values,
              TURN_TYPE* /*tt*/, int /*unused*/) {
    using namespace boost::program_options;
    validators::check_first_occurrence(v);
    const std::string& s = validators::get_single_string(values);
    if (iequals(s, "none")) {
        v = TURN_NONE;
    } else if (iequals(s, "revcomp")) {
        v = TURN_REVCOMP;
    } else if (iequals(s, "all")) {
        v = TURN_ALL;
    } else {
        throw po::invalid_option_value(s);
    }
}

std::ostream& operator<<(std::ostream& out, const TURN_TYPE& t) {
    switch(t) {
    case TURN_NONE: out << "none"; break;
    case TURN_REVCOMP: out << "revcomp"; break;
    case TURN_ALL: out << "all"; break;
    default: out << "[UNKNOWN!]";
    }
    return out;
}

void
famfinder::get_options_description(po::options_description& main,
                                   po::options_description& adv) {

    main.add_options()
        ("db,r", po::value<fs::path>(&opts.database), "reference database")
        ("turn,t", po::value<TURN_TYPE>(&opts.turn_which)
         ->default_value(TURN_NONE, "")
         ->implicit_value(TURN_REVCOMP, ""),
         "check other strand as well\n"
         "'all' checks all four frames")
        ;

    po::options_description mid("Reference Selection");
    mid.add_options()
        ("fs-engine", po::value<ENGINE_TYPE>(&opts.engine)->default_value(ENGINE_SINA_KMER, ""),
         "search engine to use for reference selection "
         "[*pt-server*|internal]")
        ("fs-kmer-len", po::value<unsigned int>(&opts.fs_kmer_len)->default_value(10u,""),
         "length of k-mers (10)")
        ("fs-req", po::value<unsigned int>(&opts.fs_req)->default_value(1u,""),
         "required number of reference sequences (1)\n"
         "queries with less matches will be dropped")
        ("fs-min", po::value<unsigned int>(&opts.fs_min)->default_value(40u,""),
         "number of references used regardless of shared fraction (40)")
        ("fs-max", po::value<unsigned int>(&opts.fs_max)->default_value(40u,""),
         "number of references used at most (40)")
        ("fs-msc", po::value<float>(&opts.fs_msc)->default_value(.7, ""),
         "required fractional identity of references (0.7)")
        ("fs-req-full", po::value<unsigned int>(&opts.fs_req_full)->default_value(1u, ""),
         "required number of full length references (1)")
        ("fs-full-len", po::value<unsigned int>(&opts.fs_full_len)->default_value(1400u, ""),
         "minimum length of full length reference (1400)")
        ("fs-req-gaps", po::value<unsigned int>(&opts.fs_req_gaps)->default_value(10u, ""),
         "ignore references with less internal gaps (10)")
        ("fs-min-len", po::value<unsigned int>(&opts.fs_min_len)->default_value(150u, ""),
         "minimal reference length (150)")
        ;
    main.add(mid);

    po::options_description od("Advanced Reference Selection");
    od.add_options()
        ("ptdb", po::value<fs::path>(&opts.database),
         "PT server database (old name)")
        ("ptport", po::value<string>(&opts.pt_port)
         ->default_value(fmt::format(":/tmp/sina_pt_{}", getpid()), ""),
         "PT server port (:/tmp/sina_pt_<pid>)")
        ("fs-kmer-no-fast", po::bool_switch(&opts.fs_no_fast),
         "don't use fast family search")
        ("fs-kmer-mm", po::value<unsigned int>(&opts.fs_kmer_mm)->default_value(0,""),
         "allowed mismatches per k-mer (0)")
        ("fs-kmer-norel", po::bool_switch(&opts.fs_kmer_norel),
         "don't score k-mer distance relative to target length")
        ("fs-msc-max", po::value<float>(&opts.fs_msc_max)->default_value(2, ""),
         "max identity of used references (for evaluation)")
        ("fs-leave-query-out", po::bool_switch(&opts.fs_leave_query_out),
         "ignore candidate if found in reference (for evaluation)")
        ("gene-start", po::value<int>(&opts.gene_start)->default_value(0,""),
         "alignment position of first base of gene (0)")
        ("gene-end", po::value<int>(&opts.gene_end)->default_value(0,""),
         "alignment position of last base of gene (0)")
        ("fs-cover-gene", po::value<unsigned int>(&opts.fs_cover_gene)->default_value(0,""),
         "required number of references covering each gene end (0)")
        ("filter", po::value<string>(&opts.posvar_filter)->default_value(""),
         "select posvar filter")
        ("auto-filter-field", po::value<string>(&opts.posvar_autofilter_field)
         ->default_value(""), "select field for auto filter selection")
        ("auto-filter-threshold",  po::value<float>(&opts.posvar_autofilter_thres)
         ->default_value(0.8, ""), "quorum for auto filter selection (0.8)")
        ("fs-oldmatch", po::bool_switch(&opts.oldmatch),
         "use legacy family composition algorithm")
        ;
    adv.add(od);
}

void famfinder::validate_vm(po::variables_map& vm,
                            po::options_description&  /*desc*/) {
    if (vm["db"].empty() && vm["ptdb"].empty()) {
        throw logic_error("Family Finder: Must have reference database (--db/-r)");
    }
    if (not vm["ptdb"].empty()) {
        logger->warn("Option --ptdb deprecated; please use --db/-r instead");
    }
    if (not vm["ptdb"].empty() && not vm["db"].empty()) {
        throw logic_error("Family Finder: please use only new --db/-r option");
    }
    if (not fs::exists(opts.database)) {
        if (opts.database.compare(":") == 0) {
            logger->warn("Loading reference database from running ARB DB server");
        } else {
            throw logic_error(fmt::format("Reference database file {} does not exist",
                                          opts.database));
        }
    }
    if (vm["fs-req"].as<unsigned int>() < 1) {
        throw logic_error("Family Finder: fs-req must be >= 1");
    }
    if (vm["fs-oldmatch"].as<bool>() && vm["fs-engine"].as<ENGINE_TYPE>() != ENGINE_ARB_PT) {
        throw logic_error("Legacy family composition only available for pt-server engine");
    }
}

ENGINE_TYPE famfinder::get_engine() {
    return opts.engine;
}

class famfinder::impl {
public:
    search *index{nullptr};
    query_arb *arb{nullptr};
    vector<alignment_stats> vastats;

    void do_turn_check(cseq& /*c*/);
    int turn_check(const cseq& /*query*/, bool /*all*/);
    void select_astats(tray &t);
    void match(search::result_vector& results, const cseq& query);

    impl();
    impl(const impl&);
    ~impl();
    tray operator()(tray /*t*/);
};

// pimpl wrappers
famfinder::famfinder() : pimpl(new impl()) {}
famfinder::famfinder(const famfinder& o) = default;
famfinder& famfinder::operator=(const famfinder& o) = default;
famfinder::~famfinder() = default;
tray famfinder::operator()(const tray& t) { return (*pimpl)(t); }

int famfinder::turn_check(const cseq& query, bool all) {
    return pimpl->turn_check(query, all);
}

famfinder::impl::impl()
    : arb(query_arb::getARBDB(opts.database))
{
    switch(opts.engine) {
    case ENGINE_ARB_PT:
        index = query_pt_pool::get_pool(
            opts.database,
            opts.fs_kmer_len,
            not opts.fs_no_fast,
            opts.fs_kmer_norel,
            opts.fs_kmer_mm,
            opts.pt_port);
        logger->warn("Using ARB PT server for reference search");
        break;
    case ENGINE_SINA_KMER:
        index = kmer_search::get_kmer_search(
            opts.database,
            opts.fs_kmer_len,
            opts.fs_no_fast);
        logger->warn("Using internal engine for reference search");
        break;
    default:
        throw std::runtime_error("Unknown sequence search engine type");
    }
    vastats = arb->getAlignmentStats();
    //index->set_range(opts.gene_start, opts.gene_end);

    //posvar_filter
    //readonly
}


famfinder::impl::~impl() {
    delete index;
}


void
famfinder::impl::do_turn_check(cseq &c) {
    // Turning sequence.
    // Strictly, this could be considered a "modification" of the input
    // sequence. However, we're really only correcting its representation.
    // The purpose of keeping the original is so that we can compare
    // changed made to the alignment at the end. This is way easier if we
    // don't have to worry about sequence orientation.
    if (opts.turn_which != TURN_NONE) {
        switch(turn_check(c, opts.turn_which==TURN_ALL)) {
        case 0:
            c.set_attr(query_arb::fn_turn, "none");
            break;
        case 1:
            c.set_attr(query_arb::fn_turn, "reversed");
            c.reverse();
            break;
        case 2:
            c.set_attr(query_arb::fn_turn, "complemented");
            c.complement();
            break;
        case 3:
            c.set_attr(query_arb::fn_turn, "reversed and complemented");
            c.reverse();
            c.complement();
            break;
        }
    } else {
        c.set_attr(query_arb::fn_turn, "turn-check disabled");
    }
}


int
famfinder::impl::turn_check(const cseq& query, bool all) {
    search::result_vector matches;
    double score[4];
    index->find(query, matches, 1);
    score[0] = matches.empty() ? 0 : matches[0].score;

    cseq turn(query);
    turn.reverse();
    if (all) {
        index->find(turn, matches, 1);
        score[1] = matches.empty() ? 0 : matches[0].score;

        cseq comp(query);
        comp.complement();

        index->find(comp, matches, 1);
        score[2] = matches.empty() ? 0 : matches[0].score;
    } else {
        score[1] = score[2] = 0;
    }

    turn.complement();
    index->find(turn, matches, 1);
    score[3] = matches.empty() ? 0 : matches[0].score;

    double max = 0;
    int best = 0;
    for (int i = 0; i < 4; i++) {
        if (max < score[i]) {
            max = score[i], best = i;
        }
    }
    return best;
}


void
famfinder::impl::select_astats(tray& t) {
    alignment_stats *astats = nullptr;

    // load default as per --filter
    if (!opts.posvar_filter.empty()) {
        for (alignment_stats &as: vastats) {
            if (as.getName() == opts.posvar_filter
                || as.getName() == opts.posvar_filter + ":ALL"
                || as.getName() == opts.posvar_filter + ":all"
                ) {
                astats = &as;
            }
        }
    }

    // do autoselection if --auto-filter-field given
    // this uses a quorum of the alignment reference:
    // at least the fraction given by posvar_autofilter_thres
    // must agree on the filter chosen
    // fixme: prefers higher level filters, should prefer most specific
    // filter, i.e. bacteria will always beat bacteria;proteobacteria
    if (opts.posvar_autofilter_field.length() > 0) {
        auto &vc = *t.alignment_reference;
        int best_count = 0;
        using vastats_t = pair<string, alignment_stats>;
        alignment_stats *best = nullptr;
        for (alignment_stats& p: vastats) {
            string filter_name = p.getName();
            int n = 0;
            for (auto &r: vc) {
                string f = opts.posvar_filter + ":" + r.sequence->get_attr<string>(opts.posvar_autofilter_field);
                if (boost::algorithm::istarts_with(f, filter_name)) {
                    ++n;
                }
            }

            if (n > best_count) {
                best_count = n;
                best = &p;
            }
        }
        if (best_count > vc.size() * opts.posvar_autofilter_thres) {
            t.log << "autofilter: " << best->getName() << ";";
            astats = best;
        } else {
            t.log << "autofilter: no match;";
        }
    }

    if (astats == nullptr) {
        astats = new alignment_stats();
    }

    t.astats = astats;
}


tray
famfinder::impl::operator()(tray t) {
    t.alignment_reference = new search::result_vector();
    auto &vc = *t.alignment_reference;
    cseq &c = *t.input_sequence;

    do_turn_check(c);

    // FIXME: int noid = opts.realign
    bool noid = false;
    if (opts.oldmatch) {
        index->match(vc, c, opts.fs_min, opts.fs_max, opts.fs_msc, opts.fs_msc_max,
                     arb, noid, opts.fs_min_len, opts.fs_req_full,
                     opts.fs_full_len, opts.fs_cover_gene, opts.fs_leave_query_out);
    } else {
        match(vc, c);
    }

    // prepare log string for alignment reference
    fmt::memory_buffer tmp;
    for (auto &r: vc) {
        if (opts.posvar_autofilter_field.length() > 0) {
            arb->loadKey(*r.sequence, opts.posvar_autofilter_field);
        }
        arb->loadKey(*r.sequence, query_arb::fn_acc);
        arb->loadKey(*r.sequence, query_arb::fn_start);
        fmt::format_to(tmp, "{}.{}:{:.2f} ",
                       r.sequence->get_attr<string>(query_arb::fn_acc),
                       r.sequence->get_attr<string>(query_arb::fn_start, "0"),
                       r.score);
    }
    c.set_attr(query_arb::fn_family, string(tmp.data(), tmp.size()));

    // remove sequences having too few gaps
    // FIXME: this should be done in match()
    if (opts.fs_req_gaps != 0) {
        auto too_few_gaps = [&](search::result_item& i) {
            return 0 == i.sequence->size()
            || i.sequence->rbegin()->getPosition() - i.sequence->size() + 1 < opts.fs_req_gaps;
        };
        vc.erase(std::remove_if(vc.begin(), vc.end(), too_few_gaps), vc.end());
    }

    // load apropriate alignment statistics into tray
    select_astats(t);
    
    // no reference => no alignment
    if (vc.size() < opts.fs_req) {
        t.log << "unable to align: too few relatives (" << vc.size() << ");";
        delete t.alignment_reference;
        t.alignment_reference = nullptr;
        return t;
    } 

    return t;
}


void
famfinder::impl::match(search::result_vector& results, const cseq& query) {
    unsigned int min_match = opts.fs_min;
    unsigned int max_match = opts.fs_max;
    float min_score = opts.fs_msc;
    float max_score = opts.fs_msc_max;
    bool noid = false;
    unsigned int min_len = opts.fs_min_len;
    unsigned int num_full = opts.fs_req_full;
    unsigned int full_min_len = opts.fs_full_len;
    unsigned int range_cover = opts.fs_cover_gene;
    bool leave_query_out = opts.fs_leave_query_out;

    using item_t = search::result_item;

    size_t range_begin = 0, range_end = 0;
    auto is_full = [full_min_len](const item_t& result) {
        return result.sequence->size() >= full_min_len;
    };
    auto is_range_left = [range_begin](const item_t& result) {
        return result.sequence->begin()->getPosition() <= range_begin;
    };
    auto is_range_right = [range_end](const item_t& result) {
        return result.sequence->getById(result.sequence->size()-1).getPosition() >= range_end;
    };

    size_t have = 0, have_full = 0, have_cover_left = 0, have_cover_right = 0;
    auto count_good = [&](const item_t& result) {
        ++have;
        if (num_full && is_full(result)) {
            ++have_full;
        }
        if (range_cover && is_range_right(result)) {
            ++have_cover_right;
        }
        if (range_cover && is_range_left(result)) {
            ++have_cover_left;
        }
        return false;
    };

    // matches results shorter than min_len
    auto remove_short = [min_len](const item_t& result) {
        return result.sequence->size() < min_len;
    };

    // matches results sharing name with query
    auto remove_query = [&, leave_query_out](const item_t& result) {
        return leave_query_out && query.getName() == result.sequence->getName();
    };

    // matches results containing query
    auto remove_superstring = [&, noid](const item_t& result) {
        return noid && boost::algorithm::icontains(result.sequence->getBases(), query.getBases());
    };

    // matches results too similar to query
    cseq_comparator cmp(CMP_IUPAC_OPTIMISTIC, CMP_DIST_NONE, CMP_COVER_QUERY, false);
    auto remove_similar = [&, max_score](const item_t& result) {
        return max_score <= 2 && cmp(query, *result.sequence) > max_score;
    };

    auto min_reached = [&](const item_t&) {
        return have >= min_match;
    };
    auto max_reached = [&](const item_t&) {
        return have >= max_match;
    };
    auto score_good = [&](const item_t& result) {
        return result.score < min_score;
    };
    auto adds_to_full = [&](const item_t& result) {
        return num_full && have_full < num_full && is_full(result);
    };
    auto adds_to_range = [&](const item_t& result) {
        return
        (range_cover && have_cover_right < range_cover && is_range_right(result))
        || (range_cover && have_cover_left < range_cover && is_range_left(result))
        ;
    };

    auto remove = [&](const item_t& result) {
        return
        remove_short(result) ||
        remove_query(result) ||
        remove_superstring(result) ||
        remove_similar(result) || (
            min_reached(result) &&
            (max_reached(result) || !score_good(result)) &&
            !adds_to_full(result) &&
            !adds_to_range(result) ) ||
        count_good(result);
    };

    size_t max_results = max_match + 1;
    search::result_vector::iterator from;
    while (have < max_match || have_full < num_full ||
           have_cover_left < range_cover || have_cover_right < range_cover) {

        results.clear();
        index->find(query, results, max_results);
        if (results.empty()) {
            return;
        }

        have = 0, have_full = 0, have_cover_left = 0, have_cover_right = 0;
        from = std::remove_if(results.begin(), results.end(), remove);
        if (max_results >= index->size()) {
            break;
        }
        max_results *= 10;
    }

    results.erase(from, results.end());
    return;
}



} // namespace sina


#if 0
void fixme() {
    int termini_begin = -1, termini_end = -1;
    string termini = arb->getFilter("termini");
    if (!termini.empty()) {
        termini_begin = termini.find_first_of('x')+1 ;
        termini_end = termini.find_last_of('x')+1;
        logger->info("Found TERMINI filter: {} - {}",
                     termini_begin, termini_end);
    }

    // FIXME: find a good way to do this with program_options
    if (opts.gene_start < 1) {
        if (termini_begin == -1) {
            opts.gene_start = 0;
        } else {
            opts.gene_start = termini_begin;
        }
    }
    if (opts.gene_end < 1 || opts.gene_end > arb->getAlignmentWidth()) {
        if (termini_end == -1) {
            opts.gene_end = arb->getAlignmentWidth();
        } else {
            opts.gene_end = termini_end;
        }
    }
    log->info("Range of gene within alignment: {} - {}",
              opts.gene_start, opts.gene_end);
    // decrement range ... we start at 0 around here
    --opts.gene_start;
    --opts.gene_end;
}
#endif


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
