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

#include "log.h"
 
#include <list>
using std::list;

#include <string>
using std::string;

#include <sstream>
using std::stringstream;

#include <iostream>
using std::cerr;
using std::endl;

#include <map>
using std::map;
using std::pair;

#include <fstream>

#include <boost/program_options.hpp>
#include <boost/foreach.hpp>

#include "cseq.h"
#include "query_arb.h"
#include "query_arb.h"

using namespace sina;
namespace po = boost::program_options;

/// options parsing ///

struct Log::options {
    int verbosity;
    bool show_diff;
    bool show_dist;
    bool colors;
    string origdb;
    string logfile;
};

struct Log::options *Log::opts;

void
Log::get_options_description(po::options_description& main,
                             po::options_description& adv) {
    opts = new struct options;

    main.add_options()
        ("log-file",
         po::value<string>(&opts->logfile)->default_value("/dev/stderr", ""),
         "write log to <arg> (stderr)")
        ;

    po::options_description od("Logging");
    od.add_options()
        ("verbosity",
         po::value<int>(&opts->verbosity)->default_value(3,""),
         "verbosity level (3)")
        ("show-diff",
         po::bool_switch(&opts->show_diff),
         "show difference to original alignment")
        ("show-dist",
         po::bool_switch(&opts->show_dist),
         "show distance to original alignment")
        ("orig-db",
         po::value<string>(&opts->origdb),
         "ARB DB containing original alignment")
        ("colors",
         po::bool_switch(&opts->colors),
         "distinguish printed bases using colors")
        ;

    adv.add(od);
}


void
Log::validate_vm(po::variables_map& vm,
                 po::options_description& /*desc*/) {
    if (vm["orig-db"].empty()) {
        if (!vm["ptdb"].empty()) {
            opts->origdb = vm["ptdb"].as<string>();
        }
    }
}


/// pipeline stuff ///

class Log::printer : public PipeElement<tray, tray> {
private:
    friend class Log;
    printer();
    ~printer();
public:
    tray operator()(tray);
    std::string getName() const {return "log::printer";}
private:
    int sequence_num;

    // stats
    double total_sps;
    double total_error;
    double total_cpm;
    double total_idty;
    double total_bps;
    double total_score;
    std::ofstream out;
    query_arb *arb;

    std::vector<int> helix_pairs;
};



PipeElement<tray, tray>*
Log::make_printer() {
    return new printer();
}


Log::printer::printer()
    : sequence_num(0),
      total_sps(0), total_error(0), total_cpm(0), total_idty(0), total_bps(0),
      total_score(0), arb(0)
{
    if (!opts->origdb.empty()) {
        arb = query_arb::getARBDB(opts->origdb);
    }
    if (arb) {
        helix_pairs = arb->getPairs();
    } 

    out.open(opts->logfile.c_str(),std::ios_base::app);
    if (!out) {
        stringstream tmp; 
        tmp << "Unable to open file \"" << opts->logfile << "\" for writing.";
        throw std::runtime_error(tmp.str());
    }
}

Log::printer::~printer() {
    if (opts->show_dist) {
        out << "avg_sps: " << total_sps / sequence_num << endl
            << "avg_cpm: " << total_cpm / sequence_num << endl
            << "avg_idty: " << total_idty / sequence_num << endl           
            << "avg_error: " << total_error / sequence_num << endl
            << "avg_bps: " << total_bps / sequence_num << endl
            << "avg_score: " << total_score / sequence_num << endl
            ;
    }
}

static int calc_nuc_term(unsigned int term_begin, unsigned int term_end, cseq& c) {
    int n = 0;
    cseq::iterator it = c.begin();
    cseq::iterator end = c.end();

    while (it != end && it->getPosition() < term_begin) ++it;
    while (it != end && it->getPosition() < term_end) { ++it, ++n; }

    return n;
}


/// actual "filter" ///

tray
Log::printer::operator()(tray t) {
    stringstream tmp;
/*
    c.set_attr(fn_qual, std::min(100, std::max(0, (int)(-100 * c.getScore()))));
*/
    if (t.input_sequence == 0) {
        throw std::runtime_error("Received broken tray in " __FILE__);
    }

    tmp << "sequence_number: " << ++sequence_num << endl
        << "sequence_identifier: " << t.input_sequence->getName() << endl;

    if (!t.aligned_sequence) {
        out << tmp.str() 
            << "align_log_slv:" << t.logstream->str() << endl
            << query_arb::fn_fullname << ":" << t.input_sequence->get_attr<string>(query_arb::fn_fullname) << endl
            << "alignment failed!" << endl << endl;
        return t;
    }
    t.aligned_sequence->set_attr("align_log_slv", 
                                 t.logstream->str());

    cseq& c = *t.aligned_sequence;

    float bps = c.calcPairScore(helix_pairs);

    c.set_attr(query_arb::fn_nuc, (int)c.size());
    c.set_attr(query_arb::fn_bpscore, (int)(100 * bps));
    if (c.size()) {
        c.set_attr(query_arb::fn_astart, (int)c.begin()->getPosition());
        c.set_attr(query_arb::fn_astop, (int)((--c.end())->getPosition()));
    }

    
    tmp << "sequence_score: " << c.getScore() << endl;

    const std::map<string,cseq::variant>& attrs = c.get_attrs();
    pair<string,cseq::variant> ap;
    BOOST_FOREACH(ap, attrs) {
        tmp << ap.first << ": "
            << boost::apply_visitor(lexical_cast_visitor<string>(), 
                                    ap.second) << endl;
    }

    bool tmp_show_diff = false;
    if (opts->show_dist) {
        cseq o = *t.input_sequence;
        if (arb) {
            string name = o.getName();
            name = name.substr(0,name.find_first_of(' '));
            o = arb->getCseq(name);
            tmp << "len-orig: " << o.size() << endl 
                << "len-alig: " << c.size() << endl;
        }
        /*
        boost::tuple<int,int,int> p = c.compare_simple(o);
        double sps = (double) p.get<0>() / o.size();
        double error = (double) (p.get<1>()+p.get<2>()) / o.size();
        tmp << "sps: " << sps << endl
            << "error: " << error << endl
            << "matches: " << p.get<0>() << endl
            << "mismatches: " << p.get<1>() << endl
            << "overhang: " << p.get<2>() << endl;
            ;
        if (((float)c.size() / o.size()) > 0.5) {
            total_sps += sps;
            total_error += error;
        } else {
            tmp << "more than 50% of sequence unaligned. not counting towards sps total." << endl;
            }*/

//        if (bps > 0) {
        tmp << "bps: " << bps << endl;
        total_bps += bps;
//        }

        /*
        if (t.alignment_reference || t.search_result) {
            std::vector<cseq> *ref = t.alignment_reference;
            if (!ref) ref = t.search_result;
            if (ref->size() > 0) {
                BOOST_FOREACH(cseq &r, *ref) {
                    r.setScore(o.identity_with(r));
                }
                std::sort(ref->begin(), ref->end());

                double idty = ref->rbegin()->getScore();
                total_idty += idty;

                double achieved_idty = c.identity_with(*ref->rbegin());
                
                boost::tuple<int,int,int> q = ref->rbegin()->compare_simple(o);
                double cpm = (c.size()-q.get<0>()>0)?(double)(p.get<0>() - q.get<0>())/(c.size() - q.get<0>()):1;
                total_cpm += cpm;
                tmp 
                    << "cpm: " << cpm << endl
                    << "idty: "  << idty << endl
                    << "achieved_idty: " << achieved_idty << endl
                    ;
            } else {
                tmp << "reference / search result empty?" << endl;
            }
        }
        */
        total_score += c.getScore();
    }


    if ((t.alignment_reference || t.search_result) 
        && (opts->show_diff || tmp_show_diff)) {

        std::vector<cseq> *ref = t.alignment_reference;
        cseq *orig = t.input_sequence;
        if (!ref) {
            ref = t.search_result;
            orig = &*ref->rbegin();
            ref->pop_back();
        } 
        
        list<unsigned int> bad_parts = orig->find_differing_parts(*t.aligned_sequence);
        list<unsigned int>::iterator it = bad_parts.begin();
        list<unsigned int>::iterator it_end = bad_parts.end();
        ref->push_back(*orig);
        ref->push_back(*t.aligned_sequence);
        while (it != it_end) {
            cseq::idx_type begin = *it++;
            cseq::idx_type end   = *it++;
            cseq::write_alignment(tmp,*ref, begin, end, opts->colors);
        }
        ref->pop_back();
        ref->pop_back();
        tmp << endl << endl;
    }
    
    out << tmp.str() << endl;

    return t;
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
