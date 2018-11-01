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

#include "log.h"
 
#include <list>
using std::list;

#include <string>
using std::string;

#include <sstream>
using std::stringstream;

#include <iostream>
using std::endl;

#include <map>
using std::map;
using std::pair;

#include <fstream>
#include <memory>

#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/filesystem.hpp>

#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/sinks/basic_file_sink.h"

#include "cseq.h"
#include "query_arb.h"
#include "query_arb.h"

using namespace sina;
namespace po = boost::program_options;
namespace fs = boost::filesystem;
using spdlog::level::level_enum;

static auto logger = Log::create_logger("log");

// having a counter in boost::program_options is difficult:

/* Type to be handled as counter */
template<typename T>
struct counting_type {
    T val;
    counting_type(const T& t) : val(t) {}
    static T initial() { return 0; };
    static T increment(const T& t) { return t+1; }
};

/* Validator handling options of counting_type<T> */
template<typename T>
void validate(boost::any& v,
	      const std::vector<std::string>& xs,
	      counting_type<T>*, long) {
    if (v.empty()) {
        v = counting_type<T>::increment(counting_type<T>::initial());
    } else {
        v = counting_type<T>::increment(boost::any_cast<T>(v));
    }
}

/* Specialization of typed_value for counting_type
 *
 * This is necessary so the value store can contain T while handling
 * is done as counting_type<T>
 */
template<typename T>
class counting_value : public po::typed_value<counting_type<T>, char> {
public:
    /* The store is of type T, but typed_value doesn't know that. */
    using super = po::typed_value<counting_type<T>, char>;
    counting_value(T* store_to)
        : super(reinterpret_cast<counting_type<T>*>(store_to))
    {
        super::zero_tokens();
        //default_value(counting_type<T>::initial(), "");
    }

    counting_value* default_value(const T& v, const std::string& x) {
        super::default_value(counting_type<T>(v),x);
    }

    /* Same here, need turn any(int) to any(counting<int>) before notify */
    void notify(const boost::any& value_store) const override {
        const auto* value = boost::any_cast<T>(&value_store);
        if (value) {
            boost::any vs(*reinterpret_cast<const counting_type<T>*>(value));
            super::notify(vs);
        } else {
            super::notify(value_store);
        }
    }
};


template<typename T>
po::typed_value<counting_type<T> >* counter(T* t) {
  return new counting_value<T>(t);
}

template<typename T>
po::typed_value<counting_type<T> >* counter() {
  return new counting_value<T>(0);
}



struct Log::options {
    int quiet_count{0};
    int verbose_count{0};
    level_enum verbosity{spdlog::level::warn};
    bool show_diff{false};
    bool show_dist{false};
    bool colors{false};
    fs::path origdb;
    fs::path logfile;
};

std::unique_ptr<Log::options> Log::opts;

static std::vector<spdlog::sink_ptr> sinks;

namespace spdlog { namespace level {
std::ostream& operator<<(std::ostream& out, const level_enum& i) {
    return out << to_c_str(i);
}
}}

void
Log::get_options_description(po::options_description& main,
                             po::options_description& adv) {
    opts = std::unique_ptr<options>(new options);

    main.add_options()
        ("verbose,v", counter<int>(&opts->verbose_count), "increase verbosity")
        ("quiet,q", counter<int>(&opts->quiet_count), "decrease verbosity")
        ("log-file", po::value<fs::path>(&opts->logfile), "file to write log to")
        ;

    po::options_description od("Logging");
    od.add_options()
        ("show-diff", po::bool_switch(&opts->show_diff), "show difference to original alignment")
        ("show-dist", po::bool_switch(&opts->show_dist), "show distance to original alignment")
        ("orig-db", po::value<fs::path>(&opts->origdb), "ARB DB containing original alignment")
        ("colors", po::bool_switch(&opts->colors), "distinguish printed bases using colors")
        ;

    adv.add(od);
}

void
Log::validate_vm(po::variables_map& vm,
                 po::options_description& /*desc*/) {
    // calculate effective log level
    auto verbosity = static_cast<int>(opts->verbosity);
    verbosity += opts->quiet_count;
    verbosity -= opts->verbose_count;
    opts->verbosity = static_cast<level_enum>(verbosity);
    opts->verbosity = std::min(
        std::max(opts->verbosity, spdlog::level::trace),
        spdlog::level::off);

    // create logging sinks
    auto console_sink = sinks[0];
    console_sink->set_level(opts->verbosity);
    console_sink->set_pattern("%T [%n] %^%v%$");

    logger->info("Loglevel set to {}", opts->verbosity);

    if (vm.count("log-file")) {
        auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(
            opts->logfile.native(), true);
        file_sink->set_level(std::min(spdlog::level::info, opts->verbosity));
        sinks.push_back(file_sink);
    }

    // update sinks of pre-existing loggers
    spdlog::apply_all([&](std::shared_ptr<spdlog::logger> l) {
            l->sinks() = sinks;
            l->set_level(spdlog::level::trace);
        });

    // database for computing distance to test case
    if (vm["orig-db"].empty()) {
        if (!vm["db"].empty()) {
            opts->origdb = vm["db"].as<fs::path>();
        }
    }
}

std::shared_ptr<spdlog::logger>
Log::create_logger(std::string name) {
    auto logger = spdlog::get(name);
    if (logger) {
        return logger;
    }
    if (sinks.empty()) {
        sinks.push_back(std::make_shared<spdlog::sinks::stderr_color_sink_mt>());
    }
    logger = std::make_shared<spdlog::logger>(name, sinks.begin(), sinks.end());
    spdlog::register_logger(logger);
    return logger;
}


/// pipeline stuff ///

struct Log::printer::priv_data {
    ~priv_data();

    int sequence_num{0};

    // stats
    double total_sps{0};
    double total_error{0};
    double total_cpm{0};
    double total_idty{0};
    double total_bps{0};
    double total_score{0};

    std::ofstream out;

    query_arb *arb{nullptr};

    std::vector<int> helix_pairs;
};


Log::printer::printer()
    : data(new priv_data())
{
    if (!opts->origdb.empty()) {
        data->arb = query_arb::getARBDB(opts->origdb);
    }
    if (data->arb) {
        data->helix_pairs = data->arb->getPairs();
    } 

    /*
    data->out.open(opts->logfile.c_str(), std::ios_base::app);

    if (!data->out) {
        stringstream tmp; 
        tmp << "Unable to open file \"" << opts->logfile << "\" for writing.";
        throw std::runtime_error(tmp.str());
    }
    */
}

Log::printer::printer(const printer& o) = default;
Log::printer& Log::printer::operator=(const printer& o) = default;
Log::printer::~printer() = default;

Log::printer::priv_data::~priv_data() {
    if (Log::opts->show_dist) {
        logger->info("avg_sps: {}", total_sps / sequence_num);
        logger->info("avg_cpm: {}", total_cpm / sequence_num);
        logger->info("avg_idty: {}", total_idty / sequence_num);
        logger->info("avg_error: {}", total_error / sequence_num);
        logger->info("avg_bps: {}", total_bps / sequence_num);
        logger->info("avg_score: {}", total_score / sequence_num);
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
    if (t.input_sequence == nullptr) {
        throw std::runtime_error("Received broken tray in " __FILE__);
    }

    tmp << "sequence_number: " << t.seqno << endl;
    tmp << "sequence_identifier: " << t.input_sequence->getName() << endl;

    if (!t.aligned_sequence) {
        data->out << tmp.str()
                  << "align_log_slv:" << t.log.str() << endl
                  << query_arb::fn_fullname << ":"
                  << t.input_sequence->get_attr<string>(query_arb::fn_fullname) << endl
                  << "alignment failed!" << endl << endl;
        return t;
    }
    t.aligned_sequence->set_attr("align_log_slv", 
                                 t.log.str());

    cseq& c = *t.aligned_sequence;

    float bps = c.calcPairScore(data->helix_pairs);

    c.set_attr(query_arb::fn_nuc, (int)c.size());
    c.set_attr(query_arb::fn_bpscore, (int)(100 * bps));
    if (c.size()) {
        c.set_attr(query_arb::fn_astart, (int)c.begin()->getPosition());
        c.set_attr(query_arb::fn_astop, (int)((--c.end())->getPosition()));
    }

    tmp << "sequence_score: " << c.getScore() << endl;

    const std::map<string,cseq::variant>& attrs = c.get_attrs();
    for (auto& ap: attrs) {
        tmp << ap.first << ": "
            << boost::apply_visitor(lexical_cast_visitor<string>(), 
                                    ap.second) << endl;
    }

    bool tmp_show_diff = false;
    if (opts->show_dist) {
        cseq o = *t.input_sequence;
        if (data->arb) {
            string name = o.getName();
            name = name.substr(0,name.find_first_of(' '));
            o = data->arb->getCseq(name);
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
        data->total_bps += bps;
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
        data->total_score += c.getScore();
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
        auto it = bad_parts.begin();
        auto it_end = bad_parts.end();
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

    logger->info(tmp.str());

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
