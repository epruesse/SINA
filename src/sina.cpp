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
#include "log.h"
#include "search_filter.h"
#include "timer.h"
#include "cseq_comparator.h"

using namespace sina;

static auto logger = Log::create_logger("SINA");

// define new type of configuration selection of input/output type
enum SEQUENCE_DB_TYPE {
    SEQUENCE_DB_NONE,
    SEQUENCE_DB_AUTO,
    SEQUENCE_DB_ARB,
    SEQUENCE_DB_FASTA,
};

// make above type printable
std::ostream& operator<<(std::ostream& out,
                         const SEQUENCE_DB_TYPE& db) {
    switch(db) {
    case SEQUENCE_DB_NONE:  out << "NONE";  break;
    case SEQUENCE_DB_AUTO:  out << "AUTO";  break;
    case SEQUENCE_DB_ARB:   out << "ARB";   break;
    case SEQUENCE_DB_FASTA: out << "FASTA"; break;
    default:                out << "Undef!";
    }
    return out;
}

// make above type parseable by boost::program_options
void validate(boost::any& v,
              const vector<string>& values,
              SEQUENCE_DB_TYPE* /*db*/, int /*unused*/) {
    po::validators::check_first_occurrence(v);
    const std::string& s = po::validators::get_single_string(values);
    if (iequals(s, "NONE")) {
        v = SEQUENCE_DB_NONE;
    } else if (iequals(s, "AUTO")) {
        v = SEQUENCE_DB_AUTO;
    } else if (iequals(s, "ARB")) {
        v = SEQUENCE_DB_ARB;
    } else if (iequals (s, "FASTA")) {
        v = SEQUENCE_DB_FASTA;
    } else { throw po::invalid_option_value(s);
}
}

// make known any<> types printable
std::ostream& operator<<(std::ostream& out,
                         const boost::any& a) {
    using boost::any_cast;
    if (any_cast<bool>(&a) != nullptr) {
        out << any_cast<bool>(a);
    } else if (any_cast<int>(&a) != nullptr) {
        out << any_cast<int>(a);
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
    } else if (any_cast<CMP_IUPAC_TYPE>(&a) != nullptr) {
        out << any_cast<CMP_IUPAC_TYPE>(a);
    } else if (any_cast<CMP_DIST_TYPE>(&a) != nullptr) {
        out << any_cast<CMP_DIST_TYPE>(a);
    } else if (any_cast<CMP_COVER_TYPE>(&a) != nullptr) {
        out << any_cast<CMP_COVER_TYPE>(a);
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
    SEQUENCE_DB_TYPE outtype;
    fs::path in;
    fs::path out;
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
         "override input file type")
        ("outtype",
         po::value<SEQUENCE_DB_TYPE>(&opts.outtype)->default_value(SEQUENCE_DB_AUTO),
         "override output file type")
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
        ("out,o", po::value<fs::path>(&opts.out)->default_value(""),
         "output file (arb or fasta)")
        ("add-relatives", po::value<unsigned int>(&opts.copy_relatives)->default_value(0, ""),
         "add the ARG nearest relatives for each sequence to output")
        ("search,S", po::bool_switch(&opts.do_search), "enable search stage")
        ("prealigned,P", po::bool_switch(&opts.skip_align), "skip alignment stage")
        ("threads,p", po::value<unsigned int>(&opts.threads)->default_value(tbb_automatic, ""),
         "limit number of threads (automatic)")
        ("num-pts", po::value<unsigned int>(&opts.num_pt_servers)->default_value(1, ""),
         "number of PT servers to start (1)")
        ("version,V", "show version")
        ;
}

void validate_vm(po::variables_map& vm, const po::options_description&  /*all_od*/) {
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
        } else {
            opts.intype = SEQUENCE_DB_FASTA;
        }
    }

    if (opts.intype == SEQUENCE_DB_NONE) {
        throw logic_error("Input type NONE invalid - need something to process");
    }

    // Output default is infile if that is ARB, else "-"
    if (opts.out == "") {
        if (opts.intype == SEQUENCE_DB_ARB) {
            opts.out = opts.in;
        } else {
            opts.out = "-";
        }
    }

    // Autodetect / validate outtype selection
    if (opts.outtype == SEQUENCE_DB_AUTO) {
        if (opts.out.extension() == ".arb" || opts.out.native() == ":") {
            opts.outtype = SEQUENCE_DB_ARB;
        } else if (opts.out == "/dev/null") {
            opts.outtype = SEQUENCE_DB_NONE;
        } else {
            opts.outtype = SEQUENCE_DB_FASTA;
        }
    }

    // Split copy_fields
    boost::split(opts.v_fields, opts.fields, boost::is_any_of(":,"));
    if (opts.v_fields.back().empty()) {
        opts.v_fields.pop_back();
    }
    if (opts.fields.find(query_arb::fn_fullname) == string::npos) {
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
    aligner::get_options_description(od, adv_od);
    famfinder::get_options_description(od, adv_od);
    search_filter::get_options_description(od, adv_od);
    query_pt::get_options_description(od, adv_od);

    po::options_description all_od(od);
    all_od.add(adv_od);

    try {
        po::store(po::parse_command_line(argc,argv,all_od),vm);

        if (vm.count("help") != 0u) {
            show_help(&od);
        }
        if (vm.count("help-all") != 0u) {
            show_help(&od, &adv_od);
        }

        po::notify(vm);

        validate_vm(vm, all_od);
        Log::validate_vm(vm, all_od);
        rw_arb::validate_vm(vm, all_od);
        rw_fasta::validate_vm(vm, all_od);
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
    vector<std::unique_ptr<tf::graph_node>> nodes; // Nodes (for cleanup)
    tf::sender<tray> *last_node; // Last tray producing node

    using source_node = tf::source_node<tray>;
    using filter_node = tf::function_node<tray, tray>;
    using limiter_node = tf::limiter_node<tray>;
    filter_node *node;

    // Make source node reading sequences
    source_node *source; // will be activated once graph complete
    switch (opts.intype) {
    case SEQUENCE_DB_ARB:
        source = new source_node(g, rw_arb::reader(opts.in,
                                                   opts.v_fields),
                                 false);
        break;
    case SEQUENCE_DB_FASTA:
        source = new source_node(g, rw_fasta::reader(opts.in,
                                                     opts.v_fields),
                                 false);
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
        // Make node(s) finding reference set
        using tray_and_finder = tuple<tray, famfinder::finder>;
        using finder_node = tf::multifunction_node<tray_and_finder, tray_and_finder>;
        using finder_node_out = finder_node::output_ports_type;
        using finder_buffer_node = tf::buffer_node<famfinder::finder>;
        using tray_and_finder_join_node = tf::join_node<tray_and_finder>;

        auto *buffer = new finder_buffer_node(g);
        auto *join = new tray_and_finder_join_node(g);

        finder_node *family_find = new finder_node(
            g, opts.num_pt_servers,
            [&](const tray_and_finder &in, finder_node_out &out) -> void {
                const tray& t = get<0>(in);
                famfinder::finder finder(get<1>(in));
                tray t_out = finder(t);
                get<0>(out).try_put(t_out);
                get<1>(out).try_put(finder);
            });


        // ---tray--> join ---------> family_find --> tray
        //            ^                         |
        //            |                         |
        //            +-finder- buffer <-finder-+
        tf::make_edge(*buffer, tf::input_port<1>(*join));
        tf::make_edge(tf::output_port<1>(*family_find), *buffer);
        nodes.emplace_back(join);
        nodes.emplace_back(buffer);

        tf::make_edge(*last_node, tf::input_port<0>(*join));
        tf::make_edge(*join, *family_find);
        nodes.emplace_back(family_find);

        buffer->try_put(famfinder::finder(0));
        tbb::parallel_for(1U, opts.num_pt_servers, [&](unsigned int i) {
                logger->warn("Launching PT server no {}", i);
                buffer->try_put(famfinder::finder(i));
            });

        node = new filter_node(g, tf::unlimited, aligner());
        tf::make_edge(tf::output_port<0>(*family_find), *node);
        nodes.emplace_back(node);
        last_node = node;
    }

    if (opts.do_search) {
        node = new filter_node(g, 1, search_filter());
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
    switch(opts.outtype) {
    case SEQUENCE_DB_ARB:
        node = new filter_node(g, 1, rw_arb::writer(opts.out,
                                                    opts.copy_relatives,
                                                    opts.v_fields));
        break;
    case SEQUENCE_DB_FASTA:
        node = new filter_node(g, 1, rw_fasta::writer(opts.out,
                                                      opts.copy_relatives,
                                                      opts.v_fields));
        break;
    case SEQUENCE_DB_NONE:
        node = nullptr;
        break;
    default:
        throw logic_error("output type undefined");
    }
    if (node != nullptr) {
        tf::make_edge(*last_node, *node);
        nodes.emplace_back(node);
        last_node = node;
    }

    node = new filter_node(g, 1, Log::printer());
    tf::make_edge(*last_node, *node);
    nodes.emplace_back(node);
    last_node = node;

    int count = 0;
    tf::function_node<tray, tf::continue_msg, tf::lightweight>
        sink(g, 1, [&](tray t) -> tf::continue_msg {
                count++;
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
    } catch (std::exception &e) {
        logger->critical("Error during program execution: {}", e.what());
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
