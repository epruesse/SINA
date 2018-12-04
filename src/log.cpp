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

#include <fstream>
#include <memory>

#include <boost/program_options.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/filesystem.hpp>

#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/sinks/basic_file_sink.h"

#include "cseq.h"
#include "cseq_comparator.h"
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
	      const std::vector<std::string>&  /*xs*/,
	      counting_type<T>* /*unused*/, long /*unused*/) {
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
} // namespace level
} // namespace spdlog

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

    if (vm.count("log-file") != 0u) {
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
    if (data->arb != nullptr) {
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
    logger->info("avg_score: {}", total_score / sequence_num);
    if (Log::opts->show_dist) {
        logger->info("avg_sps: {}", total_sps / sequence_num);
        logger->info("avg_cpm: {}", total_cpm / sequence_num);
        logger->info("avg_idty: {}", total_idty / sequence_num);
        logger->info("avg_error: {}", total_error / sequence_num);
        logger->info("avg_bps: {}", total_bps / sequence_num);
    }
}

static int calc_nuc_term(unsigned int term_begin, unsigned int term_end, cseq& c) {
    int n = 0;
    cseq::iterator it = c.begin();
    cseq::iterator end = c.end();

    while (it != end && it->getPosition() < term_begin) {
        ++it;
    }
    while (it != end && it->getPosition() < term_end) {
        ++it, ++n;
    }

    return n;
}


/// actual "filter" ///

tray
Log::printer::operator()(tray t) {
    stringstream tmp;
    if (t.input_sequence == nullptr) {
        throw std::runtime_error("Received broken tray in " __FILE__);
    }

    logger->info("sequence_number: {}", t.seqno);
    logger->info("sequence_identifier: {}", t.input_sequence->getName());

    if (t.aligned_sequence == nullptr) {
        logger->info("{}: {}", query_arb::fn_align_log, t.log.str());
        logger->info("{}: {}", query_arb::fn_fullname,
                     t.input_sequence->get_attr<string>(query_arb::fn_fullname));
        logger->info("alignment failed!");
        return t;
    }
    ++data->sequence_num;
    cseq& aligned = *t.aligned_sequence;

    float bps = aligned.calcPairScore(data->helix_pairs);
    data->total_bps += bps;
    aligned.set_attr(query_arb::fn_bpscore, (int)(100 * bps));
    aligned.set_attr(query_arb::fn_align_log, t.log.str());
    aligned.set_attr(query_arb::fn_nuc, (int)aligned.size());

    if (aligned.size() != 0u) {
        aligned.set_attr(query_arb::fn_astart, (int)aligned.begin()->getPosition());
        aligned.set_attr(query_arb::fn_astop, (int)((--aligned.end())->getPosition()));
    } else {  // shouldn't happen, but let's be careful
        aligned.set_attr(query_arb::fn_astart, 0);
        aligned.set_attr(query_arb::fn_astop, 0);
    }

    logger->info("sequence_score: {}", aligned.getScore());
    data->total_score += aligned.getScore();

    for (auto& ap: aligned.get_attrs()) {
        string val = boost::apply_visitor(lexical_cast_visitor<string>(),
                                          ap.second);
        logger->info("{}: {}", ap.first, val);
    }

    std::vector<cseq> ref;
    if (t.alignment_reference != nullptr) {
        ref = *t.alignment_reference;
    } else if (t.search_result != nullptr) {
        ref = *t.search_result;
    }

    if (opts->show_dist) {
        cseq& orig = *t.input_sequence;
        if (data->arb != nullptr) {  // we have a comparison db
            string name = orig.getName();
            orig = data->arb->getCseq(name);
            logger->info("len-orig: {}", orig.size());
            logger->info("len-alig: {}", aligned.size());
        }
        cseq_comparator cmp(CMP_IUPAC_EXACT, CMP_DIST_NONE,
                            CMP_COVER_QUERY, false);
        float sps = cmp(orig, aligned);

        logger->log((sps > 0.9999) ? spdlog::level::info : spdlog::level::warn,
                    "orig_idty: {}", sps);
        data->total_sps += sps;

        if (!ref.empty()) {
            cseq_comparator cmp(CMP_IUPAC_OPTIMISTIC, CMP_DIST_NONE,
                                CMP_COVER_QUERY, false);
            for (auto& r : ref) {
                r.setScore(cmp(orig, r));
            }
            std::sort(ref.begin(), ref.end());
            cseq &closest = *ref.rbegin();

            float orig_idty = closest.getScore();
            data->total_idty += orig_idty;
            logger->info("orig_closest_idty: {}", orig_idty);

            float aligned_idty = cmp(aligned, closest);
            logger->info("closest_idty: {}", aligned_idty);
            float cpm = orig_idty - aligned_idty;
            logger->info("cpm: {}", cpm);

            data->total_cpm += cpm;
        } else {
            tmp << "reference / search result empty?" << endl;
        }
    }

    if (opts->show_diff) {
        cseq& orig = *t.input_sequence;
        ref.push_back(orig);
        ref.push_back(aligned);
        for (auto part : orig.find_differing_parts(aligned)) {
            cseq::write_alignment(tmp, ref, part.first, part.second, opts->colors);
        }
        ref.pop_back();
        ref.pop_back();
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
