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

#include "config.h"
#include <iostream>
using std::cerr;
using std::endl;

#include <fstream>
using std::ifstream;

#include <string>
using std::string;

#include <tuple>
using std::tuple;
using std::get;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <boost/algorithm/string/predicate.hpp>
using boost::algorithm::iends_with;
using boost::algorithm::iequals;
using boost::algorithm::equals;

#include <boost/thread.hpp>
using boost::thread_group;

using std::exception;
using std::logic_error;

#ifdef HAVE_TBB
#include "tbb/flow_graph.h"
namespace tf = tbb::flow;
#endif

#include "famfinder.h"
#include "align.h"
#include "rw_arb.h"
#include "rw_fasta.h"
#include "pipe.h"
#include "log.h"
#include "search_filter.h"
#include "null_filter.h"
#include "copy_alignment.h"
#include "timer.h"
#include "cseq_comparator.h"

using namespace sina;

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
    std::cerr << "Effective parameters:" << endl;
    std::pair<std::string,po::variable_value> pv;
    BOOST_FOREACH(pv, vm) {
        std::cerr << pv.first << " = ";
        try {
            std::cerr << pv.second.value() << endl;
        } catch (boost::bad_any_cast &e) {
            std::cerr << "UNKNOWN TYPE" << endl;
        }
    }
    std::cerr << endl;
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

        ("prealigned",
         po::bool_switch(),
         "skip alignment stage")

        ("has-cli-vers",
         po::value<string>(),
         "verify support of cli version")
        ("no-align",
         po::bool_switch(),
         "disable alignment stage (deprecated)")

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
        ("help-all",
         "show full help (long)")
        ("in,i",
         po::value<string>(),
         "input file (arb or fasta)")
        ("out,o",
         po::value<string>(),
         "output file (arb or fasta)")

        ("version",
         "show version")
        ("search",
         po::bool_switch(),
         "enable search stage")
        ;
    return glob;
}

// do the messy option parsing/validating
po::variables_map
parse_options(int argc, char** argv) {
    po::variables_map vm;

    string infile, outfile;
    thread_group threads;
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
            cerr << opts << endl;
            exit(EXIT_SUCCESS);
        }
        if (vm.count("help-all")) {
            cerr << opts << endl << adv_opts << endl;
            exit(EXIT_SUCCESS);
        }
        if (vm.count("has-cli-vers")) {
            cerr << "** SINA (SILVA Incremental Aligner) " << PACKAGE_VERSION
                 << " present" << endl;
            string requested =  vm["has-cli-vers"].as<string>();
            if (requested == "1" || requested == "2" || requested == "ARB5.99") {
                exit(EXIT_SUCCESS);
            }

            cerr << "** Error: requested CLI version not supported!" << endl;
            exit(EXIT_FAILURE);
        }

        if (vm.count("version")) {
            cerr << PACKAGE_STRING
#ifdef PACKAGE_BUILDINFO
                 << " (" << PACKAGE_BUILDINFO << ")"
#endif
                 << endl;
            exit(EXIT_SUCCESS);
        }


        // Autodetect / validate intype selection
        if (vm["intype"].defaulted() && vm.count("in")) {
            const string in = vm["in"].as<string>();
            if (iends_with(in, ".arb")
                || iequals(in, ":")) {
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
                cmd[1]=vm["in"].as<string>();
                po::store(po::command_line_parser(cmd).options(all_opts).run(), vm);
                break;
            case SEQUENCE_DB_FASTA: 
                // Use "input.fasta.aligned"
                cmd[0]="-o";
                if (equals(vm["in"].as<string>(), "/dev/stdin")) {
                    cmd[1]="/dev/stdout";
                } else {
                    cmd[1]=vm["in"].as<string>()+".aligned";
                }
                po::store(po::command_line_parser(cmd).options(all_opts).run(), vm);
                break;
            default:
                throw logic_error("broken output type");
            }
        }

        // Autodetect / validate outtype selection
        if (vm["outtype"].defaulted() && vm.count("out")) {
            const string in = vm["out"].as<string>();
            if (iends_with(in, ".arb") || iequals(in, ":")) {
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
        cerr << "Configuration error:" << endl
             << e.what() << endl
             << "Use \"--help\" to show options" << endl << endl;
        if (vm.count("show-conf")) {
            show_conf(vm);
        }
        exit(EXIT_FAILURE);
    }
    return vm;
}

int main(int argc, char** argv) {
    po::variables_map vm = parse_options(argc, argv);
    cerr << "This is " << PACKAGE_STRING << "." << endl;
    if (vm.count("show-conf")) {
         show_conf(vm);
    }

    typedef PipeElement<tray,tray>* kernel;

    try {
        PipeElement<void,tray> *source;
        switch (vm["intype"].as<SEQUENCE_DB_TYPE>()) {
        case SEQUENCE_DB_ARB:
            source = rw_arb::make_reader(vm);
            break;
        case SEQUENCE_DB_FASTA:
            source = rw_fasta::make_reader(vm);
            break;
        default:
            throw logic_error("input type undefined");
        }

        PipeElement<tray,void> *sink;
        switch(vm["outtype"].as<SEQUENCE_DB_TYPE>()) {
        case SEQUENCE_DB_ARB:
            sink = rw_arb::make_writer(vm);
            break;
        case SEQUENCE_DB_FASTA:
            sink = rw_fasta::make_writer(vm);
            break;
        default:
            throw logic_error("output type outdefined");
        }

        PipeElement<tray, tray> *famfinder = 0;
        PipeElement<tray, tray> *aligner;
        if (vm["no-align"].as<bool>()
            ||
            vm["prealigned"].as<bool>()) {
            aligner = copy_alignment::make_copy_alignment();
        } else {
            aligner = aligner::make_aligner();
            famfinder = famfinder::make_famfinder();
        }

        PipeElement<tray, tray> *search;
        if (vm["search"].as<bool>()) {
            search = search_filter::make_search_filter();
        } else {
            search = null_filter::make_null_filter();
        }

        PipeElement<tray, tray> *printer;
        printer = Log::make_printer();

#ifdef HAVE_TBB
        int max_trays = 10;
        tf::graph g;

        // source node
        tf::source_node<tray> node_src(g, [&](tray &t) -> bool {
                try {
                    t = source->operator()();
//                    cerr << "R";
                    return true;
                } catch (PipeEOF &peof) {
                    return false;
                }
            }, /*is_active=*/false);
        tf::limiter_node<tray> node_limit(g, max_trays);
        tf::make_edge(node_src, node_limit);

        // famfinder with pt server
        tf::buffer_node<kernel> node_famfinder_buffer(g);
        tf::join_node< tuple<tray, kernel>, tf::reserving > node_famfinder_join(g);
        typedef tf::multifunction_node< tuple<tray,kernel>, tuple<tray, kernel> > ff_node_t;
        typedef ff_node_t::output_ports_type ff_output_t;
        ff_node_t node_famfinder(g, tf::unlimited, [&](const tuple<tray, kernel> &in, ff_output_t &out) -> void {
//                cerr << "P";
                const tray &t = get<0>(in);
                kernel famfinder = get<1>(in);
                get<0>(out).try_put((*famfinder)(t));
                get<1>(out).try_put(famfinder);
            });

        tf::make_edge(node_famfinder_buffer, tf::input_port<1>(node_famfinder_join));
        tf::make_edge(node_famfinder_join, node_famfinder);
        tf::make_edge(tf::output_port<1>(node_famfinder), node_famfinder_buffer);

        tf::function_node<tray,tray> node_aligner(g, 3, [&](tray t) -> tray {
                //              cerr << "A";
                return (*aligner)(t);
            });
        if (famfinder != 0) {
            tf::make_edge(node_src, tf::input_port<0>(node_famfinder_join));
            tf::make_edge(tf::output_port<0>(node_famfinder), node_aligner);
            node_famfinder_buffer.try_put(famfinder);
            node_famfinder_buffer.try_put(famfinder::make_famfinder(1));
            node_famfinder_buffer.try_put(famfinder::make_famfinder(2));
            node_famfinder_buffer.try_put(famfinder::make_famfinder(3));
            node_famfinder_buffer.try_put(famfinder::make_famfinder(4));
            node_famfinder_buffer.try_put(famfinder::make_famfinder(5));
            node_famfinder_buffer.try_put(famfinder::make_famfinder(6));

        } else {
            tf::make_edge(node_limit, node_aligner);
        }

        tf::function_node<tray,tray> node_search(g, 1, [&](tray t) -> tray {
//                cerr << "S";
                return (*search)(t);
            });
        tf::make_edge(node_aligner, node_search);
        tf::function_node<tray,tray> node_printer(g, 1, [&](tray t) -> tray {
                //              cerr << "p";
                return (*printer)(t);
            });
        tf::make_edge(node_search, node_printer);
        tf::function_node<tray, tf::continue_msg>
            node_sink(g, 1, [&](const tray t) -> tf::continue_msg {
//                    cerr << "W";
                    (*sink)(t);
                return tf::continue_msg();
                });
        tf::make_edge(node_printer, node_sink);
        tf::make_edge(node_sink, node_limit.decrement);

        timestamp before;
        node_src.activate();
        g.wait_for_all();
        timestamp after;
        cerr << "Time for alignment phase: " << after-before << "s" << endl;
        delete sink;
#else
        if (famfinder) {
            aligner = new PipeSerialSegment<tray, tray>(*famfinder | *aligner);
        }
        Pipe p = *source | *aligner | *search | *printer | *sink;

        // workaround for double reference because aligner is pss itself
        if (famfinder) {
            delete aligner;
        }


        timestamp before;
        p.run();
        timestamp after;
        cerr << "Time for alignment phase: " << after-before << "s" << endl;

        p.destroy();
#endif
  
        cerr << "SINA finished." << endl;
    } catch (std::exception &e) {
        cerr << "Fatal error: " << e.what() << endl;
        exit(EXIT_FAILURE);
    }
    return 0;
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
