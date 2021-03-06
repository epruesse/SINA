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

#include "config.h"
#include <iostream>

#include <fstream>
using std::ifstream;

#include <string>
using std::string;

#include <vector>
using std::vector;

#include <tuple>
using std::tuple;
using std::get;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>
using boost::algorithm::iequals;

#include <boost/core/demangle.hpp>
using boost::core::demangle;

using std::exception;
using std::logic_error;

#define TBB_PREVIEW_FLOW_GRAPH_FEATURES 1
#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include <tbb/flow_graph.h>
namespace tf = tbb::flow;

#include "famfinder.h"
#include "align.h"
#include "rw_arb.h"
#include "rw_fasta.h"
#include "rw_csv.h"
#include "log.h"
#include "search_filter.h"
#include "timer.h"
#include "cseq_comparator.h"
#include "progress.h"
#include "search.h"

using namespace sina;

static auto logger = Log::create_logger("SINA");

// define new type of configuration selection of input/output type
enum SEQUENCE_DB_TYPE {
    SEQUENCE_DB_NONE,
    SEQUENCE_DB_AUTO,
    SEQUENCE_DB_ARB,
    SEQUENCE_DB_FASTA,
    SEQUENCE_DB_CSV,
};

// make above type printable
std::ostream& operator<<(std::ostream& out, const SEQUENCE_DB_TYPE& db) {
    switch(db) {
    case SEQUENCE_DB_NONE:  out << "NONE";  break;
    case SEQUENCE_DB_AUTO:  out << "AUTO";  break;
    case SEQUENCE_DB_ARB:   out << "ARB";   break;
    case SEQUENCE_DB_FASTA: out << "FASTA"; break;
    case SEQUENCE_DB_CSV:   out << "CSV"; break;
    default:                out << "Undef!";
    }
    return out;
}

template<typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& vs) {
    bool first = true;
    for (auto& v : vs) {
        if (first) {
            first = false;
        } else {
            out << " ";
        }
        out << v;
    }
    return out;
}

// make above type parseable by boost::program_options
void validate(boost::any& v,
              const vector<string>& values,
              SEQUENCE_DB_TYPE* /*db*/, int /*unused*/) {
    //po::validators::check_first_occurrence(v);
    const std::string& s = po::validators::get_single_string(values);
    if (iequals(s, "NONE")) {
        v = SEQUENCE_DB_NONE;
    } else if (iequals(s, "AUTO")) {
        v = SEQUENCE_DB_AUTO;
    } else if (iequals(s, "ARB")) {
        v = SEQUENCE_DB_ARB;
    } else if (iequals (s, "FASTA")) {
        v = SEQUENCE_DB_FASTA;
    } else if (iequals (s, "CSV")) {
        v = SEQUENCE_DB_CSV;
    } else {
        throw po::invalid_option_value(s);
    }
}


// make known any<> types printable
template <class T,
          typename = typename std::enable_if<std::is_same<T, boost::any>::value>::type >
std::ostream& operator<<(std::ostream& out,
                         const T& a) {
    using boost::any_cast;
    if (any_cast<bool>(&a) != nullptr) {
        out << any_cast<bool>(a);
    } else if (any_cast<int>(&a) != nullptr) {
        out << any_cast<int>(a);
    } else if (any_cast<unsigned int>(&a) != nullptr) {
        out << any_cast<unsigned int>(a);
    } else if (any_cast<long>(&a) != nullptr) {
        out << any_cast<long>(a);
    } else if (any_cast<float>(&a) != nullptr) {
        out << any_cast<float>(a);
    } else if (any_cast<string>(&a) != nullptr) {
        out << any_cast<string>(a);
    } else if (any_cast<TURN_TYPE>(&a) != nullptr) {
        out << any_cast<TURN_TYPE>(a);
    } else if (any_cast<OVERHANG_TYPE>(&a) != nullptr) {
        out << any_cast<OVERHANG_TYPE>(a);
    } else if (any_cast<INSERTION_TYPE>(&a) != nullptr) {
        out << any_cast<INSERTION_TYPE>(a);
    } else if (any_cast<LOWERCASE_TYPE>(&a) != nullptr) {
        out << any_cast<LOWERCASE_TYPE>(a);
    } else if (any_cast<FASTA_META_TYPE>(&a) != nullptr) {
        out << any_cast<FASTA_META_TYPE>(a);
    } else if (any_cast<SEQUENCE_DB_TYPE>(&a) != nullptr) {
        out << any_cast<SEQUENCE_DB_TYPE>(a);
    } else if (any_cast<std::vector<SEQUENCE_DB_TYPE>>(&a) != nullptr) {
        out << any_cast<std::vector<SEQUENCE_DB_TYPE>>(a);
    } else if (any_cast<CMP_IUPAC_TYPE>(&a) != nullptr) {
        out << any_cast<CMP_IUPAC_TYPE>(a);
    } else if (any_cast<CMP_DIST_TYPE>(&a) != nullptr) {
        out << any_cast<CMP_DIST_TYPE>(a);
    } else if (any_cast<CMP_COVER_TYPE>(&a) != nullptr) {
        out << any_cast<CMP_COVER_TYPE>(a);
    } else if (any_cast<ENGINE_TYPE>(&a) != nullptr) {
        out << any_cast<ENGINE_TYPE>(a);
    } else if (any_cast<fs::path>(&a) != nullptr) {
        out << any_cast<fs::path>(a);
    } else if (any_cast<std::vector<fs::path>>(&a) != nullptr) {
        out << any_cast<std::vector<fs::path>>(a);
    } else {
        out << "UNKNOWN TYPE: '" << a.type().name()<<"'";
    }
    return out;
}

void show_conf(po::variables_map& vm) {
    std::cerr << "Effective parameters:" << std::endl;
    for (auto& pv: vm) {
        std::cerr << pv.first << " = ";
        try {
            std::cerr << pv.second.value() << std::endl;
        } catch (boost::bad_any_cast &e) {
            std::cerr << "UNKNOWN TYPE" << std::endl;
        }
    }
    std::cerr << std::endl;
}

struct options {
    SEQUENCE_DB_TYPE intype;
    fs::path in;
    std::vector<SEQUENCE_DB_TYPE> outtype;
    std::vector<fs::path> out;
    std::vector<std::pair<SEQUENCE_DB_TYPE, fs::path>> out_merged;
    unsigned int copy_relatives;
    bool noalign;
    bool skip_align;
    bool do_search;
    bool inorder;
    unsigned int threads;
    unsigned int num_pt_servers;
    unsigned int max_trays;
    string has_cli_vers;
    string fields;
    vector<string> v_fields;
};

static options opts;

// define hidden options
void get_options_description(po::options_description& main,
                             po::options_description& adv) {
    int tbb_threads = tbb::task_scheduler_init::default_num_threads();
    unsigned int tbb_automatic = tbb::task_scheduler_init::automatic;

    adv.add_options()
        ("show-conf", "show effective configuration")
        ("intype",
         po::value<SEQUENCE_DB_TYPE>(&opts.intype)->default_value(SEQUENCE_DB_AUTO),
         "override input file type [*auto*|none|arb|fasta|csv]")
        ("outtype",
         po::value<std::vector<SEQUENCE_DB_TYPE>>(&opts.outtype),
         "override output file type for next output file [*auto*|none|arb|fasta|csv]")
        ("preserve-order", po::bool_switch(&opts.inorder),
         "maintain order of sequences")
        ("max-in-flight", po::value<unsigned int>(&opts.max_trays)
         ->default_value(tbb_threads*2),
         "max number of sequences processed at a time")
        ("has-cli-vers", po::value<string>(&opts.has_cli_vers), "verify support of cli version")
        ("no-align", po::bool_switch(&opts.noalign),
         "disable alignment stage (same as prealigned)")
        ("fields,f", po::value<string>(&opts.fields),
         "select fields to write")
        ;

    main.add_options()
        ("help,h", "show short help")
        ("help-all,H", "show full help (long)")
        ("in,i", po::value<fs::path>(&opts.in)->default_value("-"),
         "input file (arb or fasta)")
        ("out,o", po::value<std::vector<fs::path>>(&opts.out)->multitoken(),
         "output file (arb, fasta or csv; may be specified multiple times)")
        ("add-relatives", po::value<unsigned int>(&opts.copy_relatives)->default_value(0, ""),
         "add the ARG nearest relatives for each sequence to output")
        ("search,S", po::bool_switch(&opts.do_search), "enable search stage")
        ("prealigned,P", po::bool_switch(&opts.skip_align), "skip alignment stage")
        ("threads,p", po::value<unsigned int>(&opts.threads)->default_value(tbb_automatic, ""),
         "limit number of threads (automatic)")
        ("num-pts", po::value<unsigned int>(&opts.num_pt_servers)->default_value(tbb_threads),
         "number of PT servers to start")
        ("version,V", "show version")
        ;
}

void validate_vm(po::variables_map& vm, const po::options_description&  /*all_od*/,
                 po::parsed_options& parsed_options) {
    if (vm.count("has-cli-vers") != 0u) {
        std::cerr << "** SINA (SILVA Incremental Aligner) " << PACKAGE_VERSION
                  << " present" << std::endl;
        const char* supported_versions[]{"1", "2", "ARB5.99"};
        for (auto& supported : supported_versions) {
            if (opts.has_cli_vers == supported) {
                exit(EXIT_SUCCESS);
            }
        }

        std::cerr << "** Error: requested CLI version not supported!" << std::endl;
        exit(EXIT_FAILURE);
    }

    if (vm.count("version") != 0u) {
        std::cerr << PACKAGE_STRING
#ifdef PACKAGE_BUILDINFO
                  << " (" << PACKAGE_BUILDINFO << ")"
#endif
                  << std::endl;
        exit(EXIT_SUCCESS);
    }

    // Autodetect / validate intype selection
    if (opts.intype == SEQUENCE_DB_AUTO) {
        if (opts.in.extension() == ".arb" || opts.in.native() == ":") {
            opts.intype = SEQUENCE_DB_ARB;
        } else if ((opts.in.extension() == ".csv") ||
                   (opts.in.extension() == ".gz" &&
                    opts.in.stem().extension() == ".csv")) {
            opts.intype = SEQUENCE_DB_CSV;
        } else {
            opts.intype = SEQUENCE_DB_FASTA;
        }
    }

    if (opts.intype == SEQUENCE_DB_NONE) {
        throw logic_error("Input type NONE invalid - need something to process");
    }
    if (opts.intype == SEQUENCE_DB_CSV) {
        throw logic_error("Input type CSV invalid - can't parse sequences from that");
    }

    SEQUENCE_DB_TYPE type_val = SEQUENCE_DB_AUTO;
    int type_idx = 0, out_idx = 0;
    for (auto &opt : parsed_options.options) {
        if (opt.string_key == "outtype") {
            type_val = opts.outtype[type_idx++];
        } else if (opt.string_key == "out") {
            for (size_t i = 0; i < opt.value.size(); i++) {
                fs::path out = opts.out[out_idx++];
                SEQUENCE_DB_TYPE outtype = type_val;
                if (outtype == SEQUENCE_DB_AUTO) {
                    if (out.extension() == ".arb" || out.native() == ":") {
                        outtype = SEQUENCE_DB_ARB;
                    } else if (out == "/dev/null") {
                        continue;
                    } else if (
                        (out.extension() == ".csv") ||
                        (out.extension() == ".gz" && out.stem().extension() == ".csv")
                        ) {
                        outtype = SEQUENCE_DB_CSV;
                    } else {
                        outtype = SEQUENCE_DB_FASTA;
                    }
                }
                opts.out_merged.emplace_back(outtype, out);
            }
            type_val = SEQUENCE_DB_AUTO;
        }
    }

    if (out_idx == 0) { // no --out specified
        if (opts.intype == SEQUENCE_DB_ARB) {
            opts.out_merged.push_back(std::make_pair(SEQUENCE_DB_ARB, opts.in));
            logger->warn("No explicit output file provided. "
                         "Reading and writing to same ARB database.");
        } else if (type_val != SEQUENCE_DB_NONE) {
            opts.out_merged.push_back(std::make_pair(SEQUENCE_DB_FASTA, "-"));
        }
    }

    // Split copy_fields
    boost::split(opts.v_fields, opts.fields, boost::is_any_of(":,"));
    if (opts.v_fields.back().empty()) {
        opts.v_fields.pop_back();
    }
    // Add full_name if no fields are specified
    if (opts.v_fields.empty()) {
        opts.v_fields.push_back(query_arb::fn_fullname);
    }
}

void show_help(po::options_description* od,
               po::options_description* adv = nullptr) {
    std::cerr << "Usage:" << std::endl
              << " sina -i input [-o output] [--prealigned|--db reference] [--search] "
              << "[--search-db search.arb] [options]"
              << std::endl << std::endl
              << *od << std::endl;
    if (adv != nullptr) {
        std::cerr << *adv << std::endl;
    }
    exit(EXIT_SUCCESS);
}

// do the messy option parsing/validating
po::variables_map
parse_options(int argc, char** argv) {
    po::variables_map vm;

    string infile, outfile;
    po::options_description od("Options"), adv_od("Advanced Options");

    get_options_description(od, adv_od);
    Log::get_options_description(od, adv_od);
    rw_arb::get_options_description(od, adv_od);
    rw_fasta::get_options_description(od, adv_od);
    rw_csv::get_options_description(od, adv_od);
    aligner::get_options_description(od, adv_od);
    famfinder::get_options_description(od, adv_od);
    search_filter::get_options_description(od, adv_od);
    query_pt::get_options_description(od, adv_od);

    po::options_description all_od(od);
    all_od.add(adv_od);

    po::positional_options_description no_positional;
    try {
        po::command_line_parser parser(argc, argv);
        parser.options(all_od);
        parser.positional(no_positional);
        po::parsed_options parsed_options = parser.run();
        po::store(parsed_options, vm);

        if (vm.count("help") != 0u) {
            show_help(&od);
        }
        if (vm.count("help-all") != 0u) {
            show_help(&od, &adv_od);
        }

        po::notify(vm);

        validate_vm(vm, all_od, parsed_options);
        Log::validate_vm(vm, all_od);
        rw_arb::validate_vm(vm, all_od);
        rw_fasta::validate_vm(vm, all_od);
        rw_csv::validate_vm(vm, all_od);
        if (!opts.skip_align && !opts.noalign) {
            famfinder::validate_vm(vm, all_od);
            aligner::validate_vm(vm, all_od);
        }
        if (opts.do_search) {
            search_filter::validate_vm(vm, all_od);
        }
        query_pt::validate_vm(vm, all_od);
    } catch (std::logic_error &e) {
        std::cerr << "Configuration error:" << std::endl
                  << e.what() << std::endl
                  << "Use \"--help\" to show options" << std::endl
                  << std::endl;
        if (vm.count("show-conf") != 0u) {
            show_conf(vm);
        }
        exit(EXIT_FAILURE);
    }
    return vm;
}


int real_main(int argc, char** argv) {
    po::variables_map vm = parse_options(argc, argv);
    logger->warn("This is {}.", PACKAGE_STRING);
    if (vm.count("show-conf") != 0u) {
         show_conf(vm);
    }

    tbb::task_scheduler_init init(vm["threads"].as<unsigned int>());

    tf::graph g;  // Main data flow graph (pipeline)
    logger_progress p(logger, "Processing");

    vector<std::unique_ptr<tf::graph_node>> nodes; // Nodes (for cleanup)
    tf::sender<tray> *last_node; // Last tray producing node

    using source_node = tf::source_node<tray>;
    using filter_node = tf::function_node<tray, tray>;
    using limiter_node = tf::limiter_node<tray>;
    filter_node *node;

    // Make source node reading sequences
    source_node *source; // will be activated once graph complete
    switch (opts.intype) {
    case SEQUENCE_DB_ARB: {
        auto arbreader = rw_arb::reader(opts.in, opts.v_fields);
        arbreader.set_progress(p);
        source = new source_node(g, arbreader, false);
    }
        break;
    case SEQUENCE_DB_FASTA: {
        auto fastareader = rw_fasta::reader(opts.in, opts.v_fields);
        fastareader.set_progress(p);
        source = new source_node(g, fastareader, false);
    }
        break;
    default:
        throw logic_error("input type undefined");
    }
    nodes.emplace_back(source);
    last_node = source;

    // Make node limiting in-flight sequence trays
    auto *limiter = new limiter_node(g, opts.max_trays);
    tf::make_edge(*last_node, *limiter);
    nodes.emplace_back(limiter);
    last_node = limiter;

    // determine number of pt servers for search and align
    if (not opts.skip_align && not opts.noalign && opts.do_search) {
        opts.num_pt_servers /= 2;
    }
    if (opts.num_pt_servers < 1) {
        opts.num_pt_servers = 1;
    }
    if (famfinder::get_engine() == ENGINE_SINA_KMER) {
        opts.num_pt_servers = tf::unlimited;
    }

    if (opts.skip_align || opts.noalign) {
        // Just copy alignment over
        node = new filter_node(g, tf::unlimited, [](tray t) -> tray {
                t.aligned_sequence = new cseq(*t.input_sequence);
                return t;
            });
        tf::make_edge(*last_node, *node);
        nodes.emplace_back(node);
        last_node = node;
    } else {
        node = new filter_node(g, opts.num_pt_servers, famfinder());
        tf::make_edge(*last_node, *node);
        nodes.emplace_back(node);
        last_node = node;

        node = new filter_node(g, tf::unlimited, aligner());
        tf::make_edge(*last_node, *node);
        nodes.emplace_back(node);
        last_node = node;
    }

    if (opts.do_search) {
        node = new filter_node(g, opts.num_pt_servers, search_filter());
        tf::make_edge(*last_node, *node);
        nodes.emplace_back(node);
        last_node = node;
    }

    if (opts.inorder) {
        using sequencer_node = tf::sequencer_node<tray>;
        sequencer_node *node = new sequencer_node(
            g, [](const tray& t) -> int {
                return t.seqno - 1;
            });
        tf::make_edge(*last_node, *node);
        nodes.emplace_back(node);
        last_node = node;
    }

    // Make node writing sequences
    for (auto& out : opts.out_merged) {
        switch(out.first) {
        case SEQUENCE_DB_ARB:
            node = new filter_node(g, 1, rw_arb::writer(out.second,
                                                        opts.copy_relatives,
                                                        opts.v_fields));
            break;
        case SEQUENCE_DB_FASTA:
            node = new filter_node(g, 1, rw_fasta::writer(out.second,
                                                          opts.copy_relatives,
                                                          opts.v_fields));
            break;
        case SEQUENCE_DB_CSV:
            node = new filter_node(g, 1, rw_csv::writer(out.second,
                                                        opts.copy_relatives,
                                                        opts.v_fields));
            break;
        default:
            throw logic_error("output type undefined");
        }
        tf::make_edge(*last_node, *node);
        nodes.emplace_back(node);
        last_node = node;
    }

    node = new filter_node(g, 1, Log::printer());
    tf::make_edge(*last_node, *node);
    nodes.emplace_back(node);
    last_node = node;

    int count = 0;
    // needs very new tbb: tf::function_node<tray, tf::continue_msg, tf::lightweight>
    tf::function_node<tray, tf::continue_msg>
        sink(g, 1, [&](tray t) -> tf::continue_msg {
                count++;
                ++p;
                t.destroy();
                return tf::continue_msg();
            });
    tf::make_edge(sink, limiter->decrement);
    tf::make_edge(*last_node, sink);

    logger->warn("Aligner ready. Processing sequences");
    timestamp before;
    source->activate();
    g.wait_for_all();
    timestamp after;
    logger->warn("Took {} to align {} sequences ({} sequences/s)",
                 after-before, count, count/(after-before));
    nodes.clear();
    logger->warn("SINA finished.");
    return 0;
}

int main(int argc, char** argv) {
    try {
        return real_main(argc, argv);
    } catch (query_arb_exception &e) {
        logger->error(e.what());
        logger->error("The ARB database you were trying to use is likely corrupted.");
        return EXIT_FAILURE;
    } catch (std::exception &e) {
        logger->error("Error during program execution: {} {}",
                         demangle(typeid(e).name()),
                         e.what());
        return EXIT_FAILURE;
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
