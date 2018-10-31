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

#include <list>
using std::list;
using std::pair;

#include <exception>
using std::exception;
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

enum ENGINE_TYPE {
    ENGINE_ARB_PT=0,
    ENGINE_SINA_KMER=1
};


struct famfinder::options {
    TURN_TYPE turn_which;
    ENGINE_TYPE engine;

    int gene_start;
    int gene_end;

    string posvar_filter;
    string posvar_autofilter_field;
    float  posvar_autofilter_thres;

    int   fs_min;
    int   fs_max;
    float fs_msc;
    float fs_msc_max;
    bool  fs_leave_query_out;
    int   fs_req;
    int   fs_req_full;
    int   fs_full_len;
    int   fs_req_gaps;
    bool  fs_no_fast;
    int   fs_kmer_len;
    int   fs_kmer_mm;
    bool  fs_kmer_norel;
    int   fs_min_len;
    int   fs_cover_gene;

    fs::path database;
    string   pt_port;
};
struct famfinder::options *famfinder::opts;


void validate(boost::any& v,
              const std::vector<std::string>& values,
              ENGINE_TYPE* /*tt*/, int) {
    using namespace boost::program_options;
    validators::check_first_occurrence(v);
    const std::string& s = validators::get_single_string(values);
    if (iequals(s, "pt-server")) v = ENGINE_ARB_PT;
    else if (iequals(s, "internal")) v = ENGINE_SINA_KMER;
    else throw po::invalid_option_value(s);
}

std::ostream& operator<<(std::ostream& out, const ENGINE_TYPE& t) {
    switch(t) {
    case ENGINE_ARB_PT: out << "pt-server"; break;
    case ENGINE_SINA_KMER: out << "internal"; break;
    default: out << "[UNKNOWN!]";
    }
    return out;
}

void validate(boost::any& v,
              const std::vector<std::string>& values,
              TURN_TYPE* /*tt*/, int) {
    using namespace boost::program_options;
    validators::check_first_occurrence(v);
    const std::string& s = validators::get_single_string(values);
    if (iequals(s, "none")) v = TURN_NONE;
    else if (iequals(s, "revcomp")) v = TURN_REVCOMP;
    else if (iequals(s, "all")) v = TURN_ALL;
    else throw po::invalid_option_value(s);
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
    opts = new struct famfinder::options();

    main.add_options()
        ("db,r", po::value<fs::path>(&opts->database), "reference database")
        ("turn,t", po::value<TURN_TYPE>(&opts->turn_which)
         ->default_value(TURN_NONE, "")
         ->implicit_value(TURN_REVCOMP, ""),
         "check other strand as well\n"
         "'all' checks all four frames")
        ;

    po::options_description mid("Reference Selection");
    mid.add_options()
        ("fs-engine", po::value<ENGINE_TYPE>(&opts->engine),
         "search engine to use for reference selection "
         "[*pt-server*|internal]")
        ("fs-kmer-len", po::value<int>(&opts->fs_kmer_len)->default_value(10,""),
         "length of k-mers (10)")
        ("fs-req", po::value<int>(&opts->fs_req)->default_value(1,""),
         "required number of reference sequences (1)\n"
         "queries with less matches will be dropped")
        ("fs-min", po::value<int>(&opts->fs_min)->default_value(40,""),
         "number of references used regardless of shared fraction (40)")
        ("fs-max", po::value<int>(&opts->fs_max)->default_value(40,""),
         "number of references used at most (40)")
        ("fs-msc", po::value<float>(&opts->fs_msc)->default_value(.7, ""),
         "required fractional identity of references (0.7)")
        ("fs-req-full", po::value<int>(&opts->fs_req_full)->default_value(1, ""),
         "required number of full length references (1)")
        ("fs-full-len", po::value<int>(&opts->fs_full_len)->default_value(1400, ""),
         "minimum length of full length reference (1400)")
        ("fs-req-gaps", po::value<int>(&opts->fs_req_gaps)->default_value(10, ""),
         "ignore references with less internal gaps (10)")
        ("fs-min-len", po::value<int>(&opts->fs_min_len)->default_value(150, ""),
         "minimal reference length (150)")
        ;
    main.add(mid);

    po::options_description od("Advanced Reference Selection");
    od.add_options()
        ("ptdb", po::value<fs::path>(&opts->database),
         "PT server database (old name)")
        ("ptport", po::value<string>(&opts->pt_port)
         ->default_value(fmt::format(":/tmp/sina_pt_{}", getpid()), ""),
         "PT server port (:/tmp/sina_pt_<pid>)")
        ("fs-kmer-no-fast", po::bool_switch(&opts->fs_no_fast),
         "don't use fast family search")
        ("fs-kmer-mm", po::value<int>(&opts->fs_kmer_mm)->default_value(0,""),
         "allowed mismatches per k-mer (0)")
        ("fs-kmer-norel", po::bool_switch(&opts->fs_kmer_norel),
         "don't score k-mer distance relative to target length")
        ("fs-msc-max", po::value<float>(&opts->fs_msc_max)->default_value(2, ""),
         "max identity of used references (for evaluation)")
        ("fs-leave-query-out", po::bool_switch(&opts->fs_leave_query_out),
         "ignore candidate if found in reference (for evaluation)")
        ("gene-start", po::value<int>(&opts->gene_start)->default_value(0,""),
         "alignment position of first base of gene (0)")
        ("gene-end", po::value<int>(&opts->gene_end)->default_value(0,""),
         "alignment position of last base of gene (0)")
        ("fs-cover-gene", po::value<int>(&opts->fs_cover_gene)->default_value(0,""),
         "required number of references covering each gene end (0)")
        ("filter", po::value<string>(&opts->posvar_filter)->default_value(""),
         "select posvar filter")
        ("auto-filter-field", po::value<string>(&opts->posvar_autofilter_field)
         ->default_value(""), "select field for auto filter selection")
        ("auto-filter-threshold",  po::value<float>(&opts->posvar_autofilter_thres)
         ->default_value(0.8, ""), "quorum for auto filter selection (0.8)")
        ;
    adv.add(od);
}

void famfinder::validate_vm(po::variables_map& vm,
                            po::options_description& desc) {
    if (vm["db"].empty() && vm["ptdb"].empty()) {
        throw logic_error("Family Finder: PT server database not set");
    }
    if (not vm["ptdb"].empty()) {
        logger->warn("Option --ptdb deprecated; please use --db instead");
    }
    if (not vm["ptdb"].empty() && not vm["db"].empty()) {
        throw logic_error("Family Finder: please use only new --db option");
    }
    if (vm["fs-req"].as<int>() < 1) {
        throw logic_error("Family Finder: fs-req must be >= 1");
    }
}

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
    if (opts->gene_start < 1) {
        if (termini_begin == -1) {
            opts->gene_start = 0;
        } else {
            opts->gene_start = termini_begin;
        }
    }
    if (opts->gene_end < 1 || opts->gene_end > arb->getAlignmentWidth()) {
        if (termini_end == -1) {
            opts->gene_end = arb->getAlignmentWidth();
        } else {
            opts->gene_end = termini_end;
        }
    }
    log->info("Range of gene within alignment: {} - {}",
              opts->gene_start, opts->gene_end);
    // decrement range ... we start at 0 around here
    --opts->gene_start;
    --opts->gene_end;
}
#endif

class famfinder::_famfinder {
    search *index;
    query_arb *arb;
    vector<alignment_stats> vastats;
    
    void do_turn_check(cseq&);
    int turn_check(const cseq&, bool);
    void select_astats(tray &t);
    
public:
    explicit _famfinder(int n);
    _famfinder(const _famfinder&);
    ~_famfinder();
    tray operator()(tray);
    std::string getName() const {return "famfinder";}
};


famfinder::finder::finder(int n)
    : data(new _famfinder(n))
{
}

famfinder::finder::finder(const finder& o)
    : data(o.data)
{
}

famfinder::finder&
famfinder::finder::operator=(const finder& o) {
    data = o.data;
    return *this;
}

famfinder::finder::~finder() {
}

tray
famfinder::finder::operator()(tray t) {
    return (*data)(t);
}

famfinder::_famfinder::_famfinder(int n)
    : arb(query_arb::getARBDB(opts->database))
{
    string pt_port = opts->pt_port;
    // FIXME: manage the port better. This works for unix sockets, but not
    // for TCP ports.
    if (n != 0) {
        pt_port +=  boost::lexical_cast<std::string>(n);
    }
    switch(opts->engine) {
    case ENGINE_ARB_PT:
        index = new query_pt(pt_port.c_str(), opts->database.c_str(),
                             not opts->fs_no_fast,
                             opts->fs_kmer_len,
                             opts->fs_kmer_mm,
                             opts->fs_kmer_norel);
        break;
    case ENGINE_SINA_KMER:
        index = kmer_search::get_kmer_search(opts->database, opts->fs_kmer_len);
        break;
    default:
        throw std::runtime_error("Unknown sequence search engine type");
    }
    vastats = arb->getAlignmentStats();
    //index->set_range(opts->gene_start, opts->gene_end);

    //posvar_filter
    //readonly
}


famfinder::_famfinder::~_famfinder() {
    delete index;
}


void
famfinder::_famfinder::do_turn_check(cseq &c) {
    // Turning sequence.
    // Strictly, this could be considered a "modification" of the input
    // sequence. However, we're really only correcting its representation.
    // The purpose of keeping the original is so that we can compare
    // changed made to the alignment at the end. This is way easier if we
    // don't have to worry about sequence orientation.
    if (opts->turn_which != TURN_NONE) {
        switch(turn_check(c, opts->turn_which==TURN_ALL)) {
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
famfinder::_famfinder::turn_check(const cseq& query, bool all) {
    std::vector<cseq> matches;
    double score[4];

    score[0] = index->match(matches, query, 1, 1, 0.0f);

    cseq turn(query);
    turn.reverse();
    if (all) {
        score[1] = index->match(matches, turn, 1, 1, 0.0f);

        cseq comp(query);
        comp.complement();
        score[2] = index->match(matches, comp, 1, 1, 0.0f);
    } else {
        score[1] = score[2] = 0;
    }

    turn.complement();
    score[3] = index->match(matches, turn, 1, 1, 0.0f);

    double max = 0;
    int best = 0;
    for (int i = 0; i < 4; i++)
        if (max < score[i])
            max = score[i], best = i;

    return best;
}


void
famfinder::_famfinder::select_astats(tray& t) {
    alignment_stats *astats = nullptr;

    // load default as per --filter
    if (opts->posvar_filter != "") {
        for (alignment_stats &as: vastats) {
            if (as.getName() == opts->posvar_filter
                || as.getName() == opts->posvar_filter + ":ALL"
                || as.getName() == opts->posvar_filter + ":all"
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
    if (opts->posvar_autofilter_field.length() > 0) {
        vector<cseq> &vc = *t.alignment_reference;
        int best_count = 0;
        typedef pair<string, alignment_stats> vastats_t;
        alignment_stats *best = nullptr;
        for (alignment_stats& p: vastats) {
            string filter_name = p.getName();
            int n = 0;
            for (cseq &r: vc) {
                string f = opts->posvar_filter + ":" + r.get_attr<string>(opts->posvar_autofilter_field);
                if (boost::algorithm::istarts_with(f, filter_name)) {
                    ++n;
                }
            }

            if (n > best_count) {
                best_count = n;
                best = &p;
            }
        }
        if (best_count > vc.size() * opts->posvar_autofilter_thres) {
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

/* tests if cseq has less than n gaps before last base */
struct has_max_n_gaps {
    typedef bool result_type;
    const int n_gaps;
    explicit has_max_n_gaps(int _n_gaps) : n_gaps(_n_gaps) {}
    bool operator()(const cseq& c) {
        return 0 == c.size() // safety, must have bases
            || c.rbegin()->getPosition() - c.size() +1 < n_gaps; 
    }
};

tray
famfinder::_famfinder::operator()(tray t) {
    t.alignment_reference = new vector<cseq>();
    vector<cseq> &vc = *t.alignment_reference;
    cseq &c = *t.input_sequence;

    do_turn_check(c);

    // FIXME: int noid = opts->realign
    int noid = false;
    index->match(vc, c, opts->fs_min, opts->fs_max, opts->fs_msc, opts->fs_msc_max,
                 arb, noid, opts->fs_min_len, opts->fs_req_full,
                 opts->fs_full_len, opts->fs_cover_gene, opts->fs_leave_query_out);
    
    stringstream tmp;
    for (cseq &r: vc) {
        if (opts->posvar_autofilter_field.length() > 0) {
            arb->loadKey(r,opts->posvar_autofilter_field);
        }
        arb->loadKey(r, query_arb::fn_acc);
        arb->loadKey(r, query_arb::fn_start);
        tmp << r.get_attr<string>(query_arb::fn_acc) << "."
            << r.get_attr<string>(query_arb::fn_start) << ":"
            << setprecision(2) << r.getScore() << " ";
    }
    c.set_attr(query_arb::fn_family_str, tmp.str());

    // remove sequences having too few gaps
    if (opts->fs_req_gaps) {
        vc.erase(std::remove_if(vc.begin(), vc.end(), 
                                has_max_n_gaps(opts->fs_req_gaps)), 
                 vc.end());
    }

    // load apropriate alignment statistics into tray
    select_astats(t);

    
    // no reference => no alignment
    if (vc.size() < opts->fs_req) {
        t.log << "unable to align: too few relatives (" << vc.size() << ");";
        delete t.alignment_reference;
        t.alignment_reference = nullptr;
        return t;
    } 

    return t;
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
