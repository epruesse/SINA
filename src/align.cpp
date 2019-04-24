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


#include "align.h"
#include "config.h"

#include <string>
using std::string;

#include <utility>
#include <vector>
using std::vector;

#include <list>
using std::list;
using std::pair;

#include <iostream>
using std::endl;
using std::ostream;

#include <fstream>
using std::ofstream;

#include <iterator>
using std::ostream_iterator;

#include <sstream>
using std::stringstream;

#include <exception>
using std::exception;

#include <algorithm>
using std::find_if;

#ifdef HAVE_TBB
#  include "tbb/tbb_allocator.h"
#endif

#include <boost/shared_ptr.hpp>
using boost::shared_ptr;

#include <boost/algorithm/string/predicate.hpp>
using boost::algorithm::iequals;

#include <boost/assert.hpp>
#include <boost/algorithm/string/find.hpp>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <sys/types.h>
#include <unistd.h> //for getpid()

#include "query_pt.h"
#include "pseq.h"
#include "cseq_comparator.h"
#include "query_arb.h"

#include "log.h"
auto logger = sina::Log::create_logger("align");

#include "mesh.h"
#include "mesh_debug.h"
#include "mseq.h"

namespace sina {

template<typename SCORING_SCHEME, typename MASTER>
void choose_transition(cseq& c, const cseq& orig, MASTER& m, SCORING_SCHEME& s, ostream& log);

template<typename transition, typename MASTER>
void do_align(cseq& c, const cseq& orig, MASTER& m, transition& tr, ostream& log);

struct aligner::options {
    bool realign;
    OVERHANG_TYPE overhang;
    LOWERCASE_TYPE lowercase;
    INSERTION_TYPE insertion;
    bool calc_idty;

    bool fs_no_graph;
    float fs_weight;

    float match_score;
    float mismatch_score;
    float gap_penalty;
    float gap_ext_penalty;

    bool debug_graph;
    bool write_used_rels;

    bool use_subst_matrix;
};
struct aligner::options *aligner::opts;


void validate(boost::any& v,
              const std::vector<std::string>& values,
              OVERHANG_TYPE* /*tt*/, int /*unused*/) {
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
              LOWERCASE_TYPE* /*tt*/, int /*unused*/) {
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
              INSERTION_TYPE* /*tt*/, int /*unused*/) {
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


void
aligner::get_options_description(po::options_description&  /*main*/,
                                 po::options_description& adv) {
    opts = new struct aligner::options();

    po::options_description od("Aligner");
    od.add_options()
        ("realign",
         po::bool_switch(&opts->realign),
         "do not copy alignment from reference")
        ("overhang",
         po::value<OVERHANG_TYPE>(&opts->overhang)->default_value(OVERHANG_ATTACH,""),
         "select type of overhang placement [*attach*|remove|edge] ")
        ("lowercase",
         po::value<LOWERCASE_TYPE>(&opts->lowercase)->default_value(LOWERCASE_NONE,""),
         "select which bases to put in lower case [*none*|original|unaligned] ")
        ("insertion",
         po::value<INSERTION_TYPE>(&opts->insertion)->default_value(INSERTION_SHIFT,""),
         "handling of insertions not accomodatable by reference alignment [*shift*|forbid|remove]")
        ("fs-no-graph",
         po::bool_switch(&opts->fs_no_graph),
         "use profile vector instead of DAG to as template")
        ("fs-weight",
         po::value<float>(&opts->fs_weight)->default_value(1,""),
         "scales weight derived from fs base freq (1)")
        ("match-score", 
         po::value<float>(&opts->match_score)->default_value(2,""),
         "score awarded for a match (2)")
        ("mismatch-score",
         po::value<float>(&opts->mismatch_score)->default_value(-1,""),
         "score awarded for a mismatch (-1)")
        ("pen-gap",
         po::value<float>(&opts->gap_penalty)->default_value(5.0,""),
         "gap open penalty (5)")
        ("pen-gapext",
         po::value<float>(&opts->gap_ext_penalty)->default_value(2.0, ""),
         "gap extend penalty (2)")
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

    adv.add(od);
}

void aligner::validate_vm(boost::program_options::variables_map& /*vm*/,
                          po::options_description& /*desc*/) {
}

} // namespace sina

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


aligner::aligner() = default;
aligner::~aligner() = default;
aligner::aligner(const aligner&) = default;
aligner& aligner::operator=(const aligner&  /*a*/) = default;


tray
aligner::operator()(tray t) {
    // skip if requirements missing
    if ((t.input_sequence == nullptr) ||
        (t.alignment_reference == nullptr) ||
        (t.astats == nullptr) ) {
        logger->error("Internal error - incomplete data for alignment");
        return t;
    }
    // prepare variables
    cseq &c = *(new cseq(*t.input_sequence)); // working copy
    search::result_vector &vc = *t.alignment_reference;
    const string bases = c.getBases(); // unaligned sequence

    if (opts->lowercase != LOWERCASE_ORIGINAL) {
        c.upperCaseAll();
    }

    // sort reference sequences containing candidate to end of family
    auto not_contains_query = [&](search::result_item& item) {
        // FIXME: we can do this w/o converting to string
        return !boost::algorithm::icontains(item.sequence->getBases(), bases);
    };
    auto begin_containing = partition(vc.begin(), vc.end(), not_contains_query);

    // if there are such sequences...
    if (begin_containing != vc.end()) {
        if (opts->realign) { // realign means ignore those sequences
            // FIXME: this should be done in famfinder to re-fill family
            t.log << "sequences ";
            for (auto it = begin_containing; it != vc.end(); ++it) {
                t.log << it->sequence->get_attr<string>(query_arb::fn_acc) << " ";
            }
            t.log << "containing exact candidate removed from family;";
            vc.erase(begin_containing, vc.end());
            if (vc.empty()) {
                t.log << "that's ALL of them. skipping sequence;";
                return t;
            }
        } else { // otherwise, we steal the alignment
            auto same_as_query = [&](search::result_item& item) {
                return iequals(bases, item.sequence->getBases());
            };
            auto exact_match = find_if(begin_containing, vc.end(), same_as_query);
            if (exact_match != vc.end()) {
                c.setAlignedBases(exact_match->sequence->getAlignedBases());
                t.log << "copied alignment from identical template sequence "
                      << exact_match->sequence->get_attr<string>(query_arb::fn_acc) << ":"
                      << exact_match->sequence->get_attr<string>(query_arb::fn_start, "0")
                      << "; ";
            } else {
                vector<aligned_base> subalignment;
                const vector<aligned_base>& refalignment = begin_containing->sequence->getAlignedBases();
                string refsequence = begin_containing->sequence->getBases();
                boost::iterator_range<string::iterator> substr = boost::ifind_first(refsequence, bases);
                subalignment.reserve(substr.size());
                std::copy(refalignment.begin() + std::distance(refsequence.begin(), substr.begin()),
                          refalignment.begin() + std::distance(refsequence.begin(), substr.end()),
                          std::back_inserter(subalignment));
                c.setAlignedBases(subalignment);
                t.log << "copied alignment from (longer) template sequence "
                      << begin_containing->sequence->get_attr<string>(query_arb::fn_acc) << ":"
                      << begin_containing->sequence->get_attr<string>(query_arb::fn_start, "0")
                      << "; ";
                BOOST_ASSERT(bases == c.getBases());
            }

            c.setWidth(begin_containing->sequence->getWidth());
            c.set_attr(query_arb::fn_date, make_datetime());
            c.set_attr(query_arb::fn_qual, 100);
            if (opts->calc_idty) {
                c.set_attr(query_arb::fn_idty, 100.f);
            }
            c.set_attr(query_arb::fn_head, 0);
            c.set_attr(query_arb::fn_tail, 0);
            c.set_attr(query_arb::fn_filter, "");
            t.aligned_sequence = &c;
            return t;
        }
    }

    std::vector<const cseq*> vcp;
    vcp.reserve(vc.size());
    for (auto& r : vc) {
        vcp.push_back(r.sequence);
    }

    if (!opts->fs_no_graph) {
        // prepare reference
        mseq m(vcp.begin(), vcp.end(), opts->fs_weight);
        // (remove duplicate edges:)
        m.sort();
        m.reduce_edges();

        if (not opts->use_subst_matrix) {
            if (t.astats->getWidth() == 0) {
                scoring_scheme_simple s(-opts->match_score, -opts->mismatch_score,
                                        opts->gap_penalty, opts->gap_ext_penalty);
                choose_transition(c, *t.input_sequence, m, s, t.log);
            } else {
                vector<float> weights = t.astats->getWeights();
                scoring_scheme_weighted s(-opts->match_score, -opts->mismatch_score,
                                          opts->gap_penalty, opts->gap_ext_penalty,
                                          weights);
                choose_transition(c, *t.input_sequence, m, s, t.log);
            }
        } else {
            vector<float> weights(vc.begin()->sequence->getWidth(), 1.f);
            if (t.astats->getWidth() == 0) { // FIXME: this looks broken
                weights = t.astats->getWeights();
            }
            float dist = vc.begin()->score;
            t.log << "using dist: " << dist << endl;
            scoring_scheme_matrix<aligned_base::matrix_type>
                s(opts->gap_penalty, opts->gap_ext_penalty, weights,
                  t.astats->getSubstMatrix(dist));
            choose_transition(c, *t.input_sequence, m, s, t.log);
        }
    } else {
        pseq p(vcp.begin(), vcp.end());
        scoring_scheme_profile s(-opts->match_score, -opts->mismatch_score,
                                 opts->gap_penalty, opts->gap_ext_penalty);
        choose_transition(c, *t.input_sequence, p, s, t.log);
    }

    if (opts->write_used_rels) {
        stringstream tmp;
        for (auto &s: vc) {
            tmp << s.sequence->getName() << " ";
        }
        c.set_attr(query_arb::fn_used_rels, tmp.str());
    }

    if (opts->calc_idty) {
        cseq_comparator calc_id(CMP_IUPAC_OPTIMISTIC,
                                CMP_DIST_NONE,
                                CMP_COVER_OVERLAP,
                                false);
        float idty = 0;
        for (auto &s: vc) {
            idty = std::max(idty, calc_id(c, *s.sequence));
        }
        c.set_attr(query_arb::fn_idty, 100.f*idty);
    }

    c.set_attr(query_arb::fn_date, make_datetime());
    c.set_attr(query_arb::fn_filter, t.astats->getName());
    t.aligned_sequence = &c;

    return t;
}

template<typename SCORING_SCHEME, typename MASTER>
void
sina::choose_transition(cseq& c, const cseq& orig, MASTER& m,
                  SCORING_SCHEME& s, ostream& log) {
    if (aligner::opts->insertion == INSERTION_FORBID) {
        transition_aspace_aware<SCORING_SCHEME, MASTER, cseq> tr(s);
        do_align(c, orig, m, tr, log);
    } else {
        transition_simple<SCORING_SCHEME, MASTER, cseq> tr(s);
        do_align(c, orig, m, tr, log);
    }
}

template<typename transition, typename MASTER>
void
sina::do_align(cseq& c, const cseq& orig, MASTER &m,
               transition &tr, ostream& log) {
    using cnsts_type = compute_node_simple<transition>;
    using data_type = typename cnsts_type::data_type;
    cnsts_type cns(tr);

    // create the alignment "mesh" (not quite a matrix)
    using mesh_t = mesh<MASTER, cseq, data_type, tbb::tbb_allocator<data_type>>;
    logger->debug("Allocating {}MB for alignment matrix", mesh_t::guessMem(m, c)/1024/1024);
    mesh_t A(m, c);

    int oh_head, oh_tail;

    // compute values of mesh nodes
    compute(A, cns);

    // remove bases from sequence container
    c.clearSequence();

    // run backtracking on the mesh
    float score = backtrack(A, c, tr,
                            aligner::opts->overhang,
                            aligner::opts->lowercase,
                            aligner::opts->insertion,
                            oh_head, oh_tail, log);
    // alignment done :-)
    c.set_attr(query_arb::fn_head, oh_head);
    c.set_attr(query_arb::fn_tail, oh_tail);
    c.set_attr(query_arb::fn_qual, (int)std::min(100.f, std::max(0.f, 100.f * score)));

    if (aligner::opts->debug_graph) {
        ofstream out(fmt::format("mseq_{}.dot", c.getName()));
        m.print_graphviz(out, "reference");

        for (auto part : orig.find_differing_parts(c)) {
            mesh_to_svg(A, part.first, part.second,
                        fmt::format("mesh_{}_{}_{}.dot",
                                    c.getName(), part.first, part.second).c_str());
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
