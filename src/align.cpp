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


#include "align.h"
#include "config.h"

#include <string>
using std::string;

#include <vector>
using std::vector;

#include <list>
using std::list;
using std::pair;

#include <iostream>
using std::endl;
using std::clog;
using std::ostream;

#include <iomanip>
using std::setprecision;

#include <fstream>
using std::ofstream;

#include <iterator>
using std::ostream_iterator;

#include <sstream>
using std::stringstream;

#include <exception>
using std::exception;
using std::logic_error;

#include <algorithm>
using std::find_if;

#include <boost/bind.hpp>
using boost::bind;

#include <boost/shared_ptr.hpp>
using boost::shared_ptr;

#include <boost/thread.hpp>
using boost::thread;

#include <boost/lexical_cast.hpp>
using boost::lexical_cast;

#include <boost/algorithm/string/predicate.hpp>
using boost::algorithm::istarts_with;
using boost::algorithm::iequals;

#include <boost/assert.hpp>
#include <boost/algorithm/string/find.hpp>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <sys/types.h>
#include <unistd.h> //for getpid()


#include "query_pt.h"
#include "mesh.h"
#include "mesh_debug.h"
#include "mseq.h"
#include "pseq.h"
#include "cseq_comparator.h"

namespace sina {

// const fieldnames for arb export
const char* aligner::fn_turn       = "turn_slv";
const char* aligner::fn_acc        = "acc";
const char* aligner::fn_start      = "start";
const char* aligner::fn_qual       = "align_quality_slv";
const char* aligner::fn_head       = "align_cutoff_head_slv";
const char* aligner::fn_tail       = "align_cutoff_tail_slv";
const char* aligner::fn_date       = "aligned_slv";
const char* aligner::fn_astart     = "align_startpos_slv";
const char* aligner::fn_astop      = "align_stoppos_slv";
const char* aligner::fn_idty       = "align_ident_slv";
const char* aligner::fn_nuc        = "nuc";
const char* aligner::fn_nuc_gene   = "nuc_gene_slv";
const char* aligner::fn_bpscore    = "align_bp_score_slv";
const char* aligner::fn_family_str ="align_family_slv";
const char* aligner::fn_used_rels  = "used_rels";
const char* aligner::fn_family     = "NONE";
const char* aligner::fn_fullname   = "full_name";

/// option stuff ///

struct aligner::options {
    bool realign;
    TURN_TYPE turn_which;
    OVERHANG_TYPE overhang;
    LOWERCASE_TYPE lowercase;
    INSERTION_TYPE insertion;
    bool calc_idty;

    string posvar_filter;
    string posvar_autofilter_field;
    float  posvar_autofilter_thres;

    int gene_start;
    int gene_end;

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
    bool  fs_no_graph;
    float fs_weight;
    int   fs_cover_gene;

    float match_score;
    float mismatch_score;
    float gap_penalty;
    float gap_ext_penalty;

    string pt_database;
    string pt_port;

    bool debug_graph;
    bool write_used_rels;

    bool use_subst_matrix;
};
struct aligner::options *aligner::opts;

void validate(boost::any& v,
              const std::vector<std::string>& values,
              TURN_TYPE* /*tt*/, int) {
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
      throw po::invalid_option_value("must be one of 'none','revcomp' or 'all'");
  }
}

std::ostream& operator<<(std::ostream& out, const TURN_TYPE& t) {
  switch(t) {
  case TURN_NONE:
    out << "none";
    break;
  case TURN_REVCOMP:
    out << "revcomp";
    break;
  case TURN_ALL:
    out << "all";
    break;
  default:
    out << "[UNKNOWN!]";
  }
  return out;
}


void validate(boost::any& v,
              const std::vector<std::string>& values,
              OVERHANG_TYPE* /*tt*/, int) {
  using namespace boost::program_options;
  validators::check_first_occurrence(v);
  const std::string& s = validators::get_single_string(values);
  if (iequals(s, "attach")) {
      v = OVERHANG_ATTACH;
  } else if (iequals(s, "remove")) {
      v = OVERHANG_REMOVE;
  } else if (iequals(s, "edge")) {
      v = OVERHANG_EDGE;
  } else {
      throw po::invalid_option_value("must be one of 'attach', 'remove' or 'edge'");
  }
}
std::ostream& operator<<(std::ostream& out, const OVERHANG_TYPE& t) {
  switch(t) {
  case OVERHANG_ATTACH:
    out << "attach";
    break;
  case OVERHANG_REMOVE:
    out << "remove";
    break;
  case OVERHANG_EDGE:
    out << "edge";
    break;
  default:
    out << "[UNKNOWN!]";
  }
  return out;
}

void validate(boost::any& v,
              const std::vector<std::string>& values,
              LOWERCASE_TYPE* /*tt*/, int) {
  using namespace boost::program_options;
  validators::check_first_occurrence(v);
  const std::string& s = validators::get_single_string(values);
  if (iequals(s, "none")) {
      v = LOWERCASE_NONE;
  } else if (iequals(s, "original")) {
      v = LOWERCASE_ORIGINAL;
  } else if (iequals(s, "unaligned")) {
      v = LOWERCASE_UNALIGNED;
  } else {
      throw po::invalid_option_value("must be one of 'none', 'original' or 'unaligned'");
  }
}
std::ostream& operator<<(std::ostream& out, const LOWERCASE_TYPE& t) {
  switch(t) {
  case LOWERCASE_NONE:
    out << "none";
    break;
  case LOWERCASE_ORIGINAL:
    out << "original";
    break;
  case LOWERCASE_UNALIGNED:
    out << "unaligned";
    break;
  default:
    out << "[UNKNOWN!]";
  }
  return out;
}

void validate(boost::any& v,
              const std::vector<std::string>& values,
              INSERTION_TYPE* /*tt*/, int) {
  using namespace boost::program_options;
  validators::check_first_occurrence(v);
  const std::string& s = validators::get_single_string(values);
  if (iequals(s, "shift")) {
      v = INSERTION_SHIFT;
  } else if (iequals(s, "forbid")) {
      v = INSERTION_FORBID;
  } else if (iequals(s, "remove")) {
      v = INSERTION_REMOVE;
  } else {
      throw po::invalid_option_value("must be one of 'shift', 'forbid' or 'remove'");
  }
}
std::ostream& operator<<(std::ostream& out, const INSERTION_TYPE& t) {
  switch(t) {
  case INSERTION_SHIFT:
    out << "shift";
    break;
  case INSERTION_FORBID:
    out << "forbid";
    break;
  case INSERTION_REMOVE:
    out << "remove";
    break;
  default:
    out << "[UNKNOWN!]";
  }
  return out;
}


po::options_description
aligner::get_options_description() {
    po::options_description od("Aligner");

    opts = new struct aligner::options();
    od.add_options()
        ("ptdb",
         po::value<string>(&opts->pt_database),
         "PT server database")

        ("ptport",
#ifdef HAVE_GETPID
         po::value<string>(&opts->pt_port)->default_value(":/tmp/sina_pt_"
                           + boost::lexical_cast<std::string>(getpid())),
#else
         po::value<string>(&opts->pt_port)->default_value("localhost:4040"),
#endif
         "PT server port")

        ("turn",
         po::value<TURN_TYPE>(&opts->turn_which)->default_value(TURN_ALL),
         "(none|revcomp|all) try reversing/complementing sequence")

        ("realign",
         po::bool_switch(&opts->realign),
         "do not copy alignment from reference")

        ("overhang",
         po::value<OVERHANG_TYPE>(&opts->overhang)->default_value(OVERHANG_ATTACH),
         "(attach|remove|edge) select type of overhang placement")

        ("lowercase",
         po::value<LOWERCASE_TYPE>(&opts->lowercase)->default_value(LOWERCASE_NONE),
         "(none|original|unaligned) select which bases to put in lower case")

        ("insertion",
         po::value<INSERTION_TYPE>(&opts->insertion)->default_value(INSERTION_SHIFT),
         "(shift|forbid|remove) handling of insertions not accomodatable by reference alignment")

        ("filter",
         po::value<string>(&opts->posvar_filter)->default_value("none"),
         "select posvar filter")

        ("auto-filter-field",
         po::value<string>(&opts->posvar_autofilter_field)->default_value(""),
         "select field for auto filter selection")

        ("auto-filter-threshold",
         po::value<float>(&opts->posvar_autofilter_thres)->default_value(0.8, "0.8"),
         "quorum for auto filter selection")

        ("fs-min",
         po::value<int>(&opts->fs_min)->default_value(40),
         "min number of used references")

        ("fs-max",
         po::value<int>(&opts->fs_max)->default_value(40),
         "max number of used references")

        ("fs-msc",
         po::value<float>(&opts->fs_msc)->default_value(.7, "0.7"),
         "min identity of used references")

        ("fs-msc-max",
         po::value<float>(&opts->fs_msc_max)->default_value(2, "none"),
         "max identity of used references (for evaluation)")

        ("fs-leave-query-out",
         po::bool_switch(&opts->fs_leave_query_out),
         "ignore candidate if found in reference (for evaluation)")

        ("fs-req",
         po::value<int>(&opts->fs_req)->default_value(1),
         "required number of references")

        ("fs-req-full",
         po::value<int>(&opts->fs_req_full)->default_value(1),
         "required number of full length references")

        ("fs-full-len",
         po::value<int>(&opts->fs_full_len)->default_value(1400),
         "minimum length of full length reference")

        ("fs-req-gaps",
         po::value<int>(&opts->fs_req_gaps)->default_value(10),
         "required number of gaps before last base")

        ("fs-kmer-no-fast",
         po::bool_switch(&opts->fs_no_fast),
         "don't use fast family search")

        ("fs-kmer-len",
         po::value<int>(&opts->fs_kmer_len)->default_value(10),
         "length of k-mers")

        ("fs-kmer-mm",
         po::value<int>(&opts->fs_kmer_mm)->default_value(0),
         "allowed mismatches per k-mer")

        ("fs-kmer-norel",
         po::bool_switch(&opts->fs_kmer_norel),
         "don't score k-mer distance relative to target length")

        ("fs-min-len",
         po::value<int>(&opts->fs_min_len)->default_value(150),
         "min reference length")

        ("fs-no-graph",
         po::bool_switch(&opts->fs_no_graph),
         "use profile vector instead of DAG to as template")

        ("fs-weight",
         po::value<float>(&opts->fs_weight)->default_value(1),
         "scales weight derived from fs base freq")

        ("gene-start",
         po::value<int>(&opts->gene_start)->default_value(0),
         "alignment position of first base of gene")

        ("gene-end",
         po::value<int>(&opts->gene_end)->default_value(0),
         "alignment position of last base of gene")

        ("fs-cover-gene",
         po::value<int>(&opts->fs_cover_gene)->default_value(0),
         "required number of references covering each gene end")

        ("match-score", 
         po::value<float>(&opts->match_score)->default_value(2),
         "score awarded for a match")

        ("mismatch-score",
         po::value<float>(&opts->mismatch_score)->default_value(-1),
         "score awarded for a mismatch")

        ("pen-gap",
         po::value<float>(&opts->gap_penalty)->default_value(5.0),
         "affine gap penalty (open)")

        ("pen-gapext",
         po::value<float>(&opts->gap_ext_penalty)->default_value(2.0),
         "affine gap penalty (extend)")

        ("debug-graph",
         po::bool_switch(&opts->debug_graph),
         "dump reference graphs to disk")

        ("use-subst-matrix",
         po::bool_switch(&opts->use_subst_matrix),
         "use experimental scoring system (slow)")

        ("write-used-rels",
         po::bool_switch(&opts->write_used_rels),
         "write used reference sequences to field 'used_rels'")

        ("calc-idty",
         po::bool_switch(&opts->calc_idty),
         "calculate highest identity of aligned sequence with any reference")
        ;


  return od;
}

void aligner::validate_vm(boost::program_options::variables_map& vm) {
    if (vm.count("no-align")) {
        return;
    }

    if (!opts) {
        throw logic_error("aligner options not parsed?!");
    }
    if (vm["ptdb"].empty()) {
        throw logic_error(string("PT server database not set"));
    }
    if (opts->fs_req < 1) {
        throw logic_error("fs-req must be >= 1");
    }

    query_arb *arb = query_arb::getARBDB(opts->pt_database);

    int termini_begin = -1, termini_end = -1;
    string termini = arb->getFilter("termini");
    if (!termini.empty()) {
        termini_begin = termini.find_first_of('x')+1 ;
        termini_end = termini.find_last_of('x')+1;
        std::cerr << "Found TERMINI filter: "
                  << termini_begin << " - " << termini_end
                  << endl;
    }

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
    std::cerr << "Range of gene within alignment: "
              << opts->gene_start << " - " << opts->gene_end
              << endl;
    // decrement range ... we start at 0 around here
    --opts->gene_start;
    --opts->gene_end;


}


} // namespace sina;

using namespace sina;

static string
make_datetime() {
        time_t  t;
        tm      tm;
        char   buf[50];

        time(&t);
        gmtime_r(&t, &tm);
        strftime(buf, 50, "%F %T", &tm);

        return string(buf);
}



class aligner::famfinder
    : public PipeElement<tray, tray > {
    friend class aligner;
    friend class galigner;
    query_pt pt;
    query_arb *arb;

    void do_turn_check(cseq&);
public:
    famfinder();
    ~famfinder();
    tray operator()(tray);
    std::string getName() const {return "famfinder";}
};

class aligner::galigner
    : public PipeElement<tray, tray> {
    friend class aligner;
    galigner();
    vector<alignment_stats> vastats;
    alignment_stats *nullstats;
    std::vector<int> pairs;
    cseq_comparator comparator;

    void select_astats(tray &t);
    template<typename SCORING_SCHEME, typename MASTER>
    void choose_transition(cseq&, cseq&, MASTER&, SCORING_SCHEME&, ostream&);
    template<typename TRANSITION, typename MASTER>
    void do_align(cseq&, cseq&, MASTER&, TRANSITION&, ostream&);
    void write_attrs(cseq &c);
public:
    tray operator()(tray);
    std::string getName() const {return "galigner";}
};


void
aligner::make_aligner(PipeElement<tray,tray> **f,
                      PipeElement<tray,tray> **g) {
    *f = new famfinder();
    *g = new galigner();
}


aligner::famfinder::famfinder()
    : pt(opts->pt_port.c_str(), opts->pt_database.c_str()),
      arb(query_arb::getARBDB(opts->pt_database))
{
    pt.set_find_type_fast(!opts->fs_no_fast);
    pt.set_probe_len(opts->fs_kmer_len);
    pt.set_mismatches(opts->fs_kmer_mm);
    pt.set_sort_type(opts->fs_kmer_norel);
    pt.set_range(opts->gene_start, opts->gene_end);

    //posvar_filter
    //readonly
}

aligner::galigner::galigner()
{
    query_arb *arb = query_arb::getARBDB(opts->pt_database);
    vastats = arb->getAlignmentStats();
    pairs = arb->getPairs();

    nullstats = new alignment_stats();
    // yeah, the above is technically a leak, as it will live
    // until program termination

    comparator = cseq_comparator(CMP_IUPAC_OPTIMISTIC,
                                 CMP_DIST_NONE,
                                 CMP_COVER_OVERLAP,
                                 false);
}
aligner::famfinder::~famfinder() {
}

void
aligner::famfinder::do_turn_check(cseq &c) {
    // Turning sequence.
    // Strictly, this could be considered a "modification" of the input
    // sequence. However, we're really only correcting its representation.
    // The purpose of keeping the original is so that we can compare
    // changed made to the alignment at the end. This is way easier if we
    // don't have to worry about sequence orientation.
    if (opts->turn_which != TURN_NONE) {
        switch(pt.turn_check(c, opts->turn_which==TURN_ALL)) {
        case 0:
            c.set_attr(fn_turn, "none");
            break;
        case 1:
            c.set_attr(fn_turn, "reversed");
            c.reverse();
            break;
        case 2:
            c.set_attr(fn_turn, "complemented");
            c.complement();
            break;
        case 3:
            c.set_attr(fn_turn, "reversed and complemented");
            c.reverse();
            c.complement();
            break;
        }
    } else {
        c.set_attr(fn_turn, "turn-check disabled");
    }
}

/* tests if cseq has less than n gaps before last base */
struct has_max_n_gaps {
    typedef bool result_type;
    const int n_gaps;
    has_max_n_gaps(int _n_gaps) : n_gaps(_n_gaps) {}
    bool operator()(const cseq& c) {
        return 0 == c.size() // safety, must have bases
            || c.rbegin()->getPosition() - c.size() +1 < n_gaps; 
    }
};


tray
aligner::famfinder::operator()(tray t) {
    t.alignment_reference = new vector<cseq>();
    vector<cseq> &vc = *t.alignment_reference;
    cseq &c = *t.input_sequence;

    do_turn_check(c);

    pt.match(vc, c, opts->fs_min, opts->fs_max, opts->fs_msc, opts->fs_msc_max,
             arb, opts->realign, opts->fs_min_len, opts->fs_req_full,
             opts->fs_full_len, opts->fs_cover_gene, opts->fs_leave_query_out);

    stringstream tmp;
    BOOST_FOREACH(cseq &r, vc) {
        if (opts->posvar_autofilter_field.length() > 0) {
            arb->loadKey(r,opts->posvar_autofilter_field);
        }
        arb->loadKey(r, fn_acc);
        arb->loadKey(r, fn_start);
        tmp << r.get_attr<string>(fn_acc) << "."
            << r.get_attr<string>(fn_start) << ":"
            << setprecision(2) << r.getScore() << " ";
    }
    c.set_attr(fn_family_str, tmp.str());

    // remove sequences having too few gaps
    if (opts->fs_req_gaps) {
        vc.erase(std::remove_if(vc.begin(), vc.end(), 
                                has_max_n_gaps(opts->fs_req_gaps)), 
                 vc.end());
    }

    // compute statistics?

    return t;
}

static int calc_nuc_term(unsigned int term_begin, unsigned int term_end, cseq& c) {
    int n = 0;
    cseq::iterator it = c.begin();
    cseq::iterator end = c.end();

    while (it != end && it->getPosition() < term_begin) ++it;
    while (it != end && it->getPosition() < term_end) { ++it, ++n; }

    return n;
}

void
aligner::galigner::select_astats(tray& t) {
    alignment_stats *astats = 0;

    // load default as per --filter
    if (opts->posvar_filter != "") {
        BOOST_FOREACH(alignment_stats &as, vastats) {
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
        alignment_stats *best=0;
        BOOST_FOREACH(alignment_stats& p, vastats) {
            string filter_name = p.getName();
            int n = 0;
            BOOST_FOREACH(cseq &r, vc) {
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
            t.log() << "autofilter: " << best->getName() << ";";
            astats = best;
        } else {
            t.log() << "autofilter: no match;";
        }
    }

    if (astats ==  0) {
        astats = new alignment_stats();
    }

    t.astats = astats;
}

struct not_icontains {
    typedef bool result_type;
    const string bases;
    not_icontains(const string& _bases) : bases(_bases) {}
    bool operator()(const cseq& c) {
        return !boost::algorithm::icontains(c.getBases(), bases);
    }
};

struct iequals_cmp {
    typedef bool result_type;
    const string bases;
    iequals_cmp(const string& _bases) : bases(_bases) {}
    bool operator()(const cseq& c) {
        return iequals(bases, c.getBases());
    }
};


tray
aligner::galigner::operator()(tray t) {
    // load apropriate alignment statistics into tray
    select_astats(t);

    // prepare variables
    cseq &c = *(new cseq(*t.input_sequence)); // working copy
    vector<cseq> &vc = *t.alignment_reference;
    string bases = c.getBases(); // unaligned sequence

    if (opts->lowercase != LOWERCASE_ORIGINAL) {
        c.upperCaseAll();
    }

    // sort reference sequences containing candidate to end of family
    vector<cseq>::iterator it;
    it = partition(vc.begin(), vc.end(), not_icontains(bases));

    // if there are such sequences...
    if (it != vc.end()) {
        if (opts->realign) { // ...either realign (throw them out)
            t.log() << "sequences ";
            for (vector<cseq>::iterator it2 = it;
                 it2 != vc.end(); ++it2) {
                t.log() << it->get_attr<string>(fn_acc) << " ";
            }
            t.log() << "containing exact candidate removed from family;";
            vc.erase(it, vc.end());
        } else { // ...or steal their alignment
            vector<cseq>::iterator exact_match = find_if(it,vc.end(),iequals_cmp(bases));
            if (exact_match != vc.end()) {
                c.setAlignedBases(exact_match->getAlignedBases());
                t.log() << "copied alignment from identical template sequence "
                        << exact_match->get_attr<string>(fn_acc) << ":"
                        << exact_match->get_attr<string>(fn_start)
                        << "; ";
            } else {
                vector<aligned_base> subalignment, refalignment;
                string refsequence = it->getBases();
                boost::iterator_range<string::iterator> substr;
                refalignment = it->getAlignedBases();

                substr = boost::ifind_first(refsequence,bases);
                subalignment.reserve(substr.size());
                std::copy( refalignment.begin() + std::distance(refsequence.begin(), substr.begin()),
                           refalignment.begin() + std::distance(refsequence.begin(), substr.end()),
                           std::back_inserter(subalignment) );


                c.setAlignedBases(subalignment);
                t.log() << "copied alignment from (longer) template sequence "
                        << it->get_attr<string>(fn_acc) << ":"
                        << it->get_attr<string>(fn_start)
                        << "; ";
                BOOST_ASSERT(bases == c.getBases());
           }
            c.setWidth(it->getWidth());
            write_attrs(c);
            c.set_attr(fn_qual, 100);
            c.set_attr(fn_idty, 100.f);
            c.set_attr(fn_head, 0);
            c.set_attr(fn_tail, 0);
            c.set_attr("align_filter_slv", "");
            t.aligned_sequence = &c;
            return t;
        }
    }

    // no reference => no alignment
    if (vc.size() < opts->fs_req) {
        t.log() << "unable to align: too few relatives (" << vc.size() << ");";
        return t;
    } 

    if (!opts->fs_no_graph) {
        // prepare reference
        mseq m(vc.begin(), vc.end(), opts->fs_weight);
        // (remove duplicate edges:)
        m.sort();
        m.reduce_edges();

        if (!opts->use_subst_matrix) {
            if (t.astats->getWidth() == 0) {
                scoring_scheme_simple s(-opts->match_score, -opts->mismatch_score,
                                        opts->gap_penalty, opts->gap_ext_penalty);
                choose_transition(c, *t.input_sequence, m, s, t.log());
            } else {
                vector<float> weights = t.astats->getWeights();
                scoring_scheme_weighted s(-opts->match_score, -opts->mismatch_score,
                                          opts->gap_penalty, opts->gap_ext_penalty,
                                          weights);
                choose_transition(c, *t.input_sequence, m, s, t.log());
            }
        } else {
            vector<float> weights(vc.begin()->getWidth(), 1.f);
            if (t.astats->getWidth() == 0) {
                weights = t.astats->getWeights();
            }
            float dist = vc.begin()->getScore();
            t.log() << "using dist: " << dist << endl;
            scoring_scheme_matrix<aligned_base::matrix_type>
                s(opts->gap_penalty, opts->gap_ext_penalty, weights,
                  t.astats->getSubstMatrix(dist));
            choose_transition(c, *t.input_sequence, m, s, t.log());
        }
    } else {
        pseq p(vc.begin(), vc.end());
        scoring_scheme_profile s(-opts->match_score, -opts->mismatch_score,
                                 opts->gap_penalty, opts->gap_ext_penalty);
        choose_transition(c, *t.input_sequence, p, s, t.log());
    }

    if (opts->write_used_rels) {
        stringstream tmp;
        BOOST_FOREACH(cseq &s, vc) {
            tmp << s.getName() << " ";
        }
        c.set_attr(fn_used_rels, tmp.str());
    }

    if (opts->calc_idty) {
        float idty = 0;
        BOOST_FOREACH(cseq &s, vc) {
            idty = std::max(idty, comparator(c, s));
        }
        c.set_attr(fn_idty, 100.f*idty);
    }

    write_attrs(c);
    c.set_attr("align_filter_slv", t.astats->getName());
    t.aligned_sequence = &c;

    return t;
}

void
aligner::galigner::write_attrs(cseq &c) {
    c.set_attr(fn_date, make_datetime());
    c.set_attr(fn_nuc_gene, (int)(calc_nuc_term(opts->gene_start, opts->gene_end, c)));
}

template<typename SCORING_SCHEME, typename MASTER>
void
aligner::galigner::choose_transition(cseq& c, cseq& orig, MASTER& m,
                                     SCORING_SCHEME& s, ostream& log) {
    if (opts->insertion == INSERTION_FORBID) {
        typedef transition_aspace_aware<SCORING_SCHEME, MASTER, cseq> transition;
        transition tr(s);
        do_align(c, orig, m, tr, log);
    } else {
        typedef transition_simple<SCORING_SCHEME, MASTER, cseq> transition;
        transition tr(s);
        do_align(c, orig, m, tr, log);
    }
}

template<typename transition, typename MASTER>
void
aligner::galigner::do_align(cseq& c, cseq& orig, MASTER &m,
                            transition &tr, ostream& log) {

    typedef compute_node_simple<transition> cnsts_type;
    typedef typename cnsts_type::data_type data_type;
    cnsts_type cns(tr);

    // create the alignment "mesh" (not quite a matrix)
    mesh<MASTER, cseq, data_type> A(m, c);

    int oh_head, oh_tail;
#ifdef DEBUG
    log << "refsize: " << m.size() << "; ";
#endif

    // compute values of mesh nodes
    compute(A, cns);

    // remove bases from sequence container
    c.clearSequence();

    // run backtracking on the mesh
    backtrack(A, c, tr, opts->overhang, opts->lowercase, opts->insertion,
              oh_head, oh_tail, log);
    // alignment done :-)
    c.set_attr(fn_head, oh_head);
    c.set_attr(fn_tail, oh_tail);
    c.set_attr(fn_qual, (int)std::min(100.f, std::max(0.f, 100.f * c.getScore())));

    if (opts->debug_graph) {
        stringstream tmp;
        tmp <<"mseq_" << c.getName() << ".dot";
        ofstream out(tmp.str().c_str());
        m.print_graphviz(out,"reference");

        list<unsigned int> bad_parts = orig.find_differing_parts(c);
        for (list<unsigned int>::iterator it=bad_parts.begin();
             it!=bad_parts.end(); ++it) {
            stringstream tmp;
            list<unsigned int>::iterator begin=it++;
            tmp << "mesh_" << c.getName() << "_" << *begin
                << "_" << *it << ".dot";
            mesh_to_svg(A, *begin, *it, tmp.str().c_str());
        }
    }
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
