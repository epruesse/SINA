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

#include <tuple>
using std::tuple;
using std::get;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

#include <boost/algorithm/string/predicate.hpp>
using boost::algorithm::iends_with;
using boost::algorithm::iequals;
using boost::algorithm::equals;

using std::exception;
using std::logic_error;

#include <tbb/task_scheduler_init.h>
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
    SEQUENCE_DB_ARB,
    SEQUENCE_DB_SILVA,
    SEQUENCE_DB_FASTA,
};

// make above type printable
std::ostream& operator<<(std::ostream& out,
                         const SEQUENCE_DB_TYPE& db) {
    switch(db) {
    case SEQUENCE_DB_NONE:
        out << "NONE";
        break;
    case SEQUENCE_DB_ARB:
        out << "ARB";
        break;
    case SEQUENCE_DB_SILVA:
        out << "SILVA";
        break;
    case SEQUENCE_DB_FASTA:
        out << "FASTA";
        break;
    default:
        out << "Undef!";
    }
    return out;
}

// make above type parseable by boost::program_options
void validate(boost::any& v,
              const std::vector<std::string>& values,
              SEQUENCE_DB_TYPE* /*db*/, int) {
    po::validators::check_first_occurrence(v);
    const std::string& s = po::validators::get_single_string(values);
    if (iequals(s, "ARB")) {
        v = SEQUENCE_DB_ARB;
    } else if (iequals (s, "FASTA")) {
        v = SEQUENCE_DB_FASTA;
    } else if (iequals (s, "SILVA")) {
        v = SEQUENCE_DB_SILVA;
    } else {
        throw po::invalid_option_value("unknown input/output type");
    }
}

// make known any<> types printable
std::ostream& operator<<(std::ostream& out,
                         const boost::any& a) {
    using boost::any_cast;
    if (any_cast<bool>(&a)) out << any_cast<bool>(a);
    else if (any_cast<int>(&a)) out << any_cast<int>(a);
    else if (any_cast<long>(&a)) out << any_cast<long>(a);
    else if (any_cast<float>(&a)) out << any_cast<float>(a);
    else if (any_cast<string>(&a)) out << any_cast<string>(a);
    else if (any_cast<TURN_TYPE>(&a)) out << any_cast<TURN_TYPE>(a);
    else if (any_cast<OVERHANG_TYPE>(&a)) out << any_cast<OVERHANG_TYPE>(a);
    else if (any_cast<INSERTION_TYPE>(&a)) out << any_cast<INSERTION_TYPE>(a);
    else if (any_cast<LOWERCASE_TYPE>(&a)) out << any_cast<LOWERCASE_TYPE>(a);
    else if (any_cast<FASTA_META_TYPE>(&a)) out << any_cast<FASTA_META_TYPE>(a);
    else if (any_cast<SEQUENCE_DB_TYPE>(&a)) out << any_cast<SEQUENCE_DB_TYPE>(a);
    else if (any_cast<CMP_IUPAC_TYPE>(&a)) out << any_cast<CMP_IUPAC_TYPE>(a);
    else if (any_cast<CMP_DIST_TYPE>(&a)) out << any_cast<CMP_DIST_TYPE>(a);
    else if (any_cast<CMP_COVER_TYPE>(&a)) out << any_cast<CMP_COVER_TYPE>(a);
    else out << "UNKNOWN TYPE: '" << a.type().name()<<"'";
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

// define hidden options
po::options_description
global_get_hidden_options_description() {
    po::options_description hidden("Advanced Options");
    hidden.add_options()
        ("show-conf",
         "show effective configuration")
        ("intype",
         po::value<SEQUENCE_DB_TYPE>()->default_value(SEQUENCE_DB_NONE,""),
         "override input file type")
        ("outtype",
         po::value<SEQUENCE_DB_TYPE>()->default_value(SEQUENCE_DB_NONE,""),
         "override output file type")

        ("has-cli-vers",
         po::value<string>(),
         "verify support of cli version")
        ("no-align",
         po::bool_switch(),
         "disable alignment stage (same as prealigned)")

        ;
    return hidden;
}

// define global options
po::options_description
global_get_options_description() {
    po::options_description glob("Options");
    glob.add_options()
        ("help,h",
         "show short help")
        ("help-all,H",
         "show full help (long)")
        ("in,i",
         po::value<fs::path>(),
         "input file (arb or fasta)")
        ("out,o",
         po::value<fs::path>(),
         "output file (arb or fasta)")
        ("search,S",
         po::bool_switch(),
         "enable search stage")
        ("prealigned,P",
         po::bool_switch(),
         "skip alignment stage")
        ("threads,p",
         po::value<unsigned int>()->default_value(
             tbb::task_scheduler_init::automatic, ""),
         "limit number of threads (automatic)")
        ("version,V",
         "show version")

        ;
    return glob;
}

// do the messy option parsing/validating
po::variables_map
parse_options(int argc, char** argv) {
    po::variables_map vm;

    string infile, outfile;
    po::options_description
        opts = global_get_options_description(),
        adv_opts = global_get_hidden_options_description();

    Log::get_options_description(opts, adv_opts);
    rw_arb::get_options_description(opts, adv_opts);
    rw_fasta::get_options_description(opts, adv_opts);
    aligner::get_options_description(opts, adv_opts);
    famfinder::get_options_description(opts, adv_opts);
    search_filter::get_options_description(opts, adv_opts);
    query_pt::get_options_description(opts, adv_opts);

    po::options_description all_opts(opts);
    all_opts.add(adv_opts);

    try {
        po::store(po::parse_command_line(argc,argv,all_opts),vm);

        if (vm.count("help")) {
            std::cerr << opts << std::endl;
            exit(EXIT_SUCCESS);
        }
        if (vm.count("help-all")) {
            std::cerr << opts << std::endl
                      << adv_opts << std::endl;
            exit(EXIT_SUCCESS);
        }
        if (vm.count("has-cli-vers")) {
            std::cerr << "** SINA (SILVA Incremental Aligner) " << PACKAGE_VERSION
                      << " present" << std::endl;
            string requested =  vm["has-cli-vers"].as<string>();
            if (requested == "1" || requested == "2" || requested == "ARB5.99") {
                exit(EXIT_SUCCESS);
            }

            std::cerr << "** Error: requested CLI version not supported!" << std::endl;
            exit(EXIT_FAILURE);
        }

        if (vm.count("version")) {
            std::cerr << PACKAGE_STRING
#ifdef PACKAGE_BUILDINFO
                      << " (" << PACKAGE_BUILDINFO << ")"
#endif
                      << std::endl;
            exit(EXIT_SUCCESS);
        }


        // Autodetect / validate intype selection
        if (vm["intype"].defaulted() && vm.count("in")) {
            const fs::path in = vm["in"].as<fs::path>();
            if (in.extension() == ".arb" || in.native() == ":") {
                std::vector<string> cmd(2);
                cmd[0]="--intype";
                cmd[1]="ARB";
                po::store(po::command_line_parser(cmd).options(all_opts).run(),
                          vm);
            } else {
                std::vector<string> cmd(2);
                cmd[0]="--intype";
                cmd[1]="FASTA";
                po::store(po::command_line_parser(cmd).options(all_opts).run(),
                          vm);
            } 
        }

        // Pick suitable output if no output given
        if (not vm.count("out") && vm.count("in")) {
            std::vector<string> cmd(4);
            switch(vm["intype"].as<SEQUENCE_DB_TYPE>()) {
            case SEQUENCE_DB_ARB: 
                // ARB files can be used for input and output
                cmd[0]="-o";
                cmd[1]=vm["in"].as<fs::path>().native();
                po::store(po::command_line_parser(cmd).options(all_opts).run(), vm);
                break;
            case SEQUENCE_DB_FASTA: 
                // Use "input.fasta.aligned"
                cmd[0]="-o";
                if (vm["in"].as<fs::path>() == "/dev/stdin") {
                    cmd[1] = "/dev/stdout";
                } else if (vm["in"].as<fs::path>() == "-") {
                    cmd[1] = "-";
                } else {
                    cmd[1] = vm["in"].as<fs::path>().native()+".aligned";
                }
                po::store(po::command_line_parser(cmd).options(all_opts).run(), vm);
                break;
            default:
                throw logic_error("broken output type");
            }
        }

        // Autodetect / validate outtype selection
        if (vm["outtype"].defaulted() && vm.count("out")) {
            const fs::path out = vm["out"].as<fs::path>();
            if (out.extension() == ".arb" || out.native() == ":") {
                std::vector<string> cmd(2);
                cmd[0]="--outtype";
                cmd[1]="ARB";
                po::store(po::command_line_parser(cmd).options(all_opts).run(),
                          vm);
            } else {
                std::vector<string> cmd(2);
                cmd[0]="--outtype";
                cmd[1]="FASTA";
                po::store(po::command_line_parser(cmd).options(all_opts).run(),
                          vm);
            }
        }

        if (vm["no-align"].as<bool>() || vm["prealigned"].as<bool>()) {
            if (vm["db"].empty() && vm["ptdb"].empty()) {
                std::vector<string> cmd(2);
                cmd[0]="--db";
                cmd[1]="NONE";
                po::store(po::command_line_parser(cmd).options(all_opts).run(),
                          vm);
            }
        }
        
        // enable calc-idty if min-idty requested
        if (!vm["min-idty"].empty()) {
            std::vector<string> cmd(1);
            cmd[0]="--calc-idty";
            po::store(po::command_line_parser(cmd).options(all_opts).run(), vm);
        }

        Log::validate_vm(vm, all_opts);
        rw_arb::validate_vm(vm, all_opts);
        rw_fasta::validate_vm(vm, all_opts);
        aligner::validate_vm(vm, all_opts);
        famfinder::validate_vm(vm, all_opts);
        search_filter::validate_vm(vm, all_opts);
        query_pt::validate_vm(vm, all_opts);

        if (vm.count("in") == 0) {
            throw logic_error("need input file");
        }

        po::notify(vm);


    } catch (std::logic_error &e) {
        std::cerr << "Configuration error:" << std::endl
                  << e.what() << std::endl
                  << "Use \"--help\" to show options" << std::endl
                  << std::endl;
        if (vm.count("show-conf")) {
            show_conf(vm);
        }
        exit(EXIT_FAILURE);
    }
    return vm;
}


int real_main(int argc, char** argv) {
    po::variables_map vm = parse_options(argc, argv);
    logger->warn("This is {}.", PACKAGE_STRING);
    if (vm.count("show-conf")) {
         show_conf(vm);
    }

    tbb::task_scheduler_init init(vm["threads"].as<unsigned int>());

    SEQUENCE_DB_TYPE intype = vm["intype"].as<SEQUENCE_DB_TYPE>();
    SEQUENCE_DB_TYPE outtype = vm["outtype"].as<SEQUENCE_DB_TYPE>();
    fs::path input_file = vm["in"].as<fs::path>();
    fs::path output_file = vm["out"].as<fs::path>();
    bool do_align = not (vm["no-align"].as<bool>() || vm["prealigned"].as<bool>());
    bool do_search = vm["search"].as<bool>();
    int max_trays = 10;

    tf::graph g;  // Main data flow graph (pipeline)
    std::vector<std::unique_ptr<tf::graph_node>> nodes; // Nodes (for cleanup)
    tf::sender<tray> *last_node; // Last tray producing node

    typedef tf::source_node<tray> source_node;
    typedef tf::function_node<tray, tray> filter_node;
    typedef tf::limiter_node<tray> limiter_node;
    filter_node *node;

    // Make source node reading sequences
    source_node *source; // will be activated once graph complete
    switch (intype) {
    case SEQUENCE_DB_ARB:
        source = new source_node(g, rw_arb::reader(input_file), false);
        break;
    case SEQUENCE_DB_FASTA:
        source = new source_node(g, rw_fasta::reader(input_file), false);
        break;
    default:
        throw logic_error("input type undefined");
    }
    nodes.emplace_back(source);
    last_node = source;

    // Make node limiting in-flight sequence trays
    limiter_node *limiter = new limiter_node(g, max_trays);
    tf::make_edge(*last_node, *limiter);
    nodes.emplace_back(limiter);
    last_node = limiter;

    if (!do_align) {
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
        typedef tuple<tray, famfinder::finder> tray_and_finder;
        typedef tf::multifunction_node<tray_and_finder, tray_and_finder> finder_node;
        typedef finder_node::output_ports_type finder_node_out;
        typedef tf::buffer_node<famfinder::finder> finder_buffer_node;
        typedef tf::join_node<tray_and_finder> tray_and_finder_join_node;

        finder_buffer_node *buffer = new finder_buffer_node(g);
        tray_and_finder_join_node *join = new tray_and_finder_join_node(g);

        finder_node *family_find = new finder_node(
            g, 1,
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

        for (int i=0; i<3; i++) {
            buffer->try_put(famfinder::finder(i));
        }

        node = new filter_node(g, tf::unlimited, aligner());
        tf::make_edge(tf::output_port<0>(*family_find), *node);
        nodes.emplace_back(node);
        last_node = node;
    }

    if (do_search) {
        node = new filter_node(g, 1, search_filter());
        tf::make_edge(*last_node, *node);
        nodes.emplace_back(node);
        last_node = node;
    }

    // Make node writing sequences
    switch(outtype) {
    case SEQUENCE_DB_ARB:
        node = new filter_node(g, 1, rw_arb::writer(output_file));
        break;
    case SEQUENCE_DB_FASTA:
        node = new filter_node(g, 1, rw_fasta::writer(output_file));
        break;
    default:
        throw logic_error("input type undefined");
    }
    tf::make_edge(*last_node, *node);
    nodes.emplace_back(node);
    last_node = node;

    node = new filter_node(g, 1, Log::printer());
    tf::make_edge(*last_node, *node);
    nodes.emplace_back(node);
    last_node = node;

    int count = 0;
    tf::function_node<tray, tf::continue_msg>
        sink(g, 1, [&](tray t) -> tf::continue_msg {
                count++;
                t.destroy();
                return tf::continue_msg();
            });
    tf::make_edge(sink, limiter->decrement);
    tf::make_edge(*last_node, sink);

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
        logger->critical("uncaught exception: {}", e.what());
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
