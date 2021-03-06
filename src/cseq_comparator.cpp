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

#include "cseq_comparator.h"

#include <boost/algorithm/string/predicate.hpp>
using boost::algorithm::istarts_with;
using boost::algorithm::iequals;

#include <boost/foreach.hpp>

#include <cmath>

namespace po = boost::program_options;

namespace sina {

static float jukes_cantor(float in) {
    return -3.0/4 * log( 1.0 - 4.0/3*in);
}


cseq_comparator::cseq_comparator(CMP_IUPAC_TYPE iupac, CMP_DIST_TYPE dist, 
				 CMP_COVER_TYPE cover, bool filter_lc) 
    : iupac_rule(iupac), dist_rule(dist), cover_rule(cover), 
      filter_lc_rule(filter_lc)
{
}

cseq_comparator::cseq_comparator() = default;

template<typename FUNC> 
void
traverse(const cseq_base& A, const cseq_base& B, FUNC F) {
    auto a = A.bases.begin();
    auto a_end = A.bases.end();
    auto b = B.bases.begin();
    auto b_end = B.bases.end();

    // skip filtered bases at beginning
    while (a != a_end && F->filtered(*a)) {
        ++a;
    }
    while (b != b_end && F->filtered(*b)) {
        ++b;
    }

    // skip filtered bases at end
    while (a != a_end && F->filtered(*(a_end-1)) ) {
        --a_end;
    }
    while (b != b_end && F->filtered(*(b_end-1)) ) {
        --b_end;
    }
    
    // count left overhang
    if (a->getPosition() < b->getPosition()) {
        // overhang in sequence A
        while (a != a_end && 
               a->getPosition() < b->getPosition()) {
            F->overhangA(*a);
            ++a;
        }
    } else {
        // overhang in sequence B or no overhang
        while (b != b_end && 
               a->getPosition() > b->getPosition()) {
            F->overhangB(*b);
            ++b;
        }
    }

    // count overlaping zone
    while (a != a_end && b != b_end) {
        int diff = a->getPosition() - b->getPosition();
        if (diff > 0) { // a>b
            F->onlyB(*b);
            ++b;
        } else if (diff < 0) { // a<b
            F->onlyA(*a);
            ++a;
        } else { // a==b
            F->both(*a,*b);
            ++a;
            ++b;
        }
    }

    // count right overhang
    while (a != a_end) { F->overhangA(*a); ++a; }
    while (b != b_end) { F->overhangB(*b); ++b; }
}

struct base_comp_optimistic {
    bool operator()(const base_iupac& a, const base_iupac& b) const {
        return a.comp(b);
    }
};

struct base_comp_pessimistic {
    bool operator()(const base_iupac& a, const base_iupac& b) const {
        return a.comp_pessimistic(b);
    }
};

struct base_comp_exact {
    bool operator()(const base_iupac& a, const base_iupac& b) const {
        return a.comp_exact(b);
    }
};

struct filter_none {
    bool filtered(const aligned_base& /*unused*/) const {
        return false;
    }
};

struct filter_lowercase {
    bool filtered(const aligned_base& b) const {
        return b.isLowerCase();
    }
};

  

struct match_counter {
    int only_a_overhang{0};
    int only_b_overhang{0};
    int only_a{0};
    int only_b{0};
    int match{0};
    int mismatch{0};

    template<typename BCOMP, typename FILTER>     
    struct counter;

    template<typename BCOMP, typename FILTER>
    counter<BCOMP, FILTER>* func() {
        return static_cast<counter<BCOMP, FILTER>* >(this);
    }
};

template<typename BCOMP, typename FILTER>     
struct match_counter::counter : public match_counter, public FILTER {
    counter(const counter&);
    void overhangA(const aligned_base& b) {
        if (!FILTER::filtered(b)) {
            only_a_overhang++;
        }
    } 
    void overhangB(const aligned_base& b) {
        if (!FILTER::filtered(b)) {
            only_b_overhang++;
        }
    } 
    void onlyA(const aligned_base& b) {
        if (!FILTER::filtered(b)) {
            only_a++;
        }
    }
    void onlyB(const aligned_base& b) {
        if (!FILTER::filtered(b)) {
            only_b++;
        }
    }
    void both(const aligned_base& b1, 
               const aligned_base& b2) {
        if (!FILTER::filtered(b1) && !FILTER::filtered(b2)) {
            if (cmp(b1,b2)) {
                match++;
            } else {
                mismatch++;
            }
        } else if (!FILTER::filtered(b1)) {
            only_a++;
        } else if (!FILTER::filtered(b2)) {
            only_b++;
        }   
    }
    BCOMP cmp;
};



float 
cseq_comparator::operator()(const cseq& query, const cseq& target) {
    match_counter m;
    int base;
   
    switch(iupac_rule) {
    case CMP_IUPAC_OPTIMISTIC:
        if (filter_lc_rule) {
            traverse(query, target, m.func<base_comp_optimistic, filter_lowercase>());
        } else {
            traverse(query, target, m.func<base_comp_optimistic, filter_none>());
        }
        break;
    case CMP_IUPAC_PESSIMISTIC:
        if (filter_lc_rule) {
            traverse(query, target, m.func<base_comp_pessimistic, filter_lowercase>());
        } else {
            traverse(query, target, m.func<base_comp_pessimistic, filter_none>());
        }
        break;
    case CMP_IUPAC_EXACT:
        if (filter_lc_rule) {
            traverse(query, target, m.func<base_comp_exact, filter_lowercase>());
        } else {
            traverse(query, target, m.func<base_comp_exact, filter_none>());
        }
        break;
    default:
        throw std::logic_error("unknown iupac rule");
    }
    
    switch(cover_rule) {
    case CMP_COVER_ABS:
        base = 1;
        break;
    case CMP_COVER_QUERY:
        base = m.match + m.mismatch + m.only_a + m.only_a_overhang;
        break;
    case CMP_COVER_TARGET:
        base = m.match + m.mismatch + m.only_b + m.only_b_overhang;
        break;
    case CMP_COVER_OVERLAP:
        base = m.match + m.mismatch + m.only_a + m.only_b;
        break;
    case CMP_COVER_ALL:
        base = m.match + m.mismatch + m.only_a + m.only_b 
            + m.only_a_overhang + m.only_b_overhang;
        break;
    case CMP_COVER_AVERAGE:
        base = m.match + m.mismatch + 
            (m.only_a + m.only_b 
             + m.only_a_overhang + m.only_b_overhang)/2;
        break;
    case CMP_COVER_MIN:
        base = m.match + m.mismatch 
            + std::min(m.only_a + m.only_a_overhang, 
                       m.only_b + m.only_b_overhang);
        break;
    case CMP_COVER_MAX:
        base = m.match + m.mismatch   
            + std::max(m.only_a + m.only_a_overhang, 
                       m.only_b + m.only_b_overhang);
        break;
    case CMP_COVER_NOGAP:
        base = m.match + m.mismatch;
        break;
    default:
        throw std::logic_error("unknown cover rule");
    }

    float dist = (float)m.match / base;

    switch(dist_rule) {
    case CMP_DIST_NONE:
        break;
    case CMP_DIST_JC:
        dist = jukes_cantor(dist);
        break;
    default:
        throw std::logic_error("unknown dist rule");
    }

    return dist;
}


////////////////////////// OPTION PARSING //////////////////////////

void
validate(boost::any& v,
	 const std::vector<std::string>& values,
	 CMP_IUPAC_TYPE* /*unused*/, int /*unused*/) {
    po::validators::check_first_occurrence(v);
    const std::string& s = po::validators::get_single_string(values);
    if (istarts_with("optimistic", s)) {
        v = CMP_IUPAC_OPTIMISTIC;
    } else if (istarts_with("pessimistic", s)) {
        v = CMP_IUPAC_PESSIMISTIC;
    } else if (istarts_with("exact", s)) {
        v = CMP_IUPAC_EXACT;
    } else {
        throw po::invalid_option_value
            ("iupac matching must be either optimistic or pessimistic");
    }
}

std::ostream& 
operator<<(std::ostream& out, const CMP_IUPAC_TYPE& t) {
    switch(t) {
    case CMP_IUPAC_OPTIMISTIC:
        out << "optimistic";
        break;
    case CMP_IUPAC_PESSIMISTIC:
        out << "pessimistic";
        break;
    case CMP_IUPAC_EXACT:
        out << "exact";
        break;
    default:
        out << "UNDEFINED!";
    }
    return out;
}

void
validate(boost::any& v,
	 const std::vector<std::string>& values,
	 CMP_DIST_TYPE* /*unused*/, int /*unused*/) {
    po::validators::check_first_occurrence(v);
    const std::string& s = po::validators::get_single_string(values);
    if (iequals(s, "none")) {
        v = CMP_DIST_NONE;
    } else if (iequals(s, "jc")) {
        v = CMP_DIST_JC;
    } else {
        throw po::invalid_option_value
            ("distance correction must be either none or jc");
    }
}

std::ostream&
operator<<(std::ostream& out, const CMP_DIST_TYPE& t) {
    switch(t) {
    case CMP_DIST_NONE:
        out << "none";
        break;
    case CMP_DIST_JC:
        out << "jc";
        break;
    default:
        out << "UNDEFINED!";
    }
    return out;
}

void
validate(boost::any& v,
	 const std::vector<std::string>& values,
	 CMP_COVER_TYPE* /*unused*/, int /*unused*/) {
    po::validators::check_first_occurrence(v);
    const std::string& s = po::validators::get_single_string(values);
    if (iequals(s, "abs")) {
        v = CMP_COVER_ABS;
    } else if (iequals(s, "query")) {
        v = CMP_COVER_QUERY;
    } else if (iequals(s, "target")) {
        v = CMP_COVER_TARGET;
    } else if (iequals(s, "overlap")) {
        v = CMP_COVER_OVERLAP;
    } else if (iequals(s, "all")) {
        v = CMP_COVER_ALL;
    } else if (iequals(s, "average")) {
        v = CMP_COVER_AVERAGE;
    } else if (iequals(s, "min")) {
        v = CMP_COVER_MIN;
    } else if (iequals(s, "max")) {
        v = CMP_COVER_MAX;
    } else if (iequals(s, "nogap")) {
        v = CMP_COVER_NOGAP;
    } else {
        throw po::invalid_option_value
            ("coverage type must be one of abs, query, target, overlap,"
             "average, nogap, min or max");
    }
}

std::ostream&
operator<<(std::ostream& out, const CMP_COVER_TYPE& t) {
    switch(t) {
    case CMP_COVER_ABS:
        out << "abs";
        break;
    case CMP_COVER_QUERY:
        out << "query";
        break;
    case CMP_COVER_TARGET:
        out << "target";
        break;
    case CMP_COVER_OVERLAP:
        out << "overlap";
        break;
    case CMP_COVER_ALL:
        out << "all";
        break;
    case CMP_COVER_AVERAGE:
        out << "average";
        break;
    case CMP_COVER_MIN:
        out << "min";
        break;
    case CMP_COVER_MAX:
        out << "max";
        break;
    case CMP_COVER_NOGAP:
        out << "nogap";
        break;
    default:
        out << "UNDEFINED!";
    }
    return out;
}


po::options_description
cseq_comparator::get_options_description(const char* prefix) {
    po::options_description od;
    std::string p(prefix);
    od.add_options()
        ((p + "iupac").c_str(),
         po::value<CMP_IUPAC_TYPE>()->default_value(CMP_IUPAC_OPTIMISTIC, ""),
         "strategy for comparing ambiguous bases [pessimistic|*optimistic*|exact]")
        
        ((p + "correction").c_str(),
         po::value<CMP_DIST_TYPE>()->default_value(CMP_DIST_NONE,""),
         "apply distance correction. [*none*|jc]")
      
        ((p + "cover").c_str(),
         po::value<CMP_COVER_TYPE>()->default_value(CMP_COVER_QUERY,""),
         "compute comparative measure relative to\n"
         "  abs: 1\n"
         " *query: query sequence length\n"
         "  target: target sequence length\n"
         "  min: length of shorter sequence\n"
         "  max: length of longer sequence\n"
         "  avg: average length\n"
         "  overlap: length of overlap\n"
         "  all: columns with bases in either\n"
         "  nogap: columns with bases in both\n")
        
        ((p + "filter-lowercase").c_str(), 
         po::bool_switch(),
         "ignore bases in lowercase when comparing sequences")
      ;
    return od;
}

cseq_comparator
cseq_comparator::make_from_variables_map(po::variables_map& vm,
                                         const char* prefix){
    std::string p(prefix);
    
    CMP_IUPAC_TYPE iupac = vm[p + "iupac"].as<CMP_IUPAC_TYPE>();
    CMP_DIST_TYPE dist = vm[p + "correction"].as<CMP_DIST_TYPE>();
    CMP_COVER_TYPE cover =  vm[p + "cover"].as<CMP_COVER_TYPE>();
    bool filter_lc = vm[p + "filter-lowercase"].as<bool>();
    
    if (CMP_COVER_ABS == cover && CMP_DIST_NONE != dist) {
        throw std::logic_error
            ("only fractional identity can be distance corrected");
    }

    return {iupac, dist, cover, filter_lc};
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
