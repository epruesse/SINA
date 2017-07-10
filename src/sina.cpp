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
    po::options_description hidden("");
    hidden.add_options()
        ("has-cli-vers",
         po::value<string>(),
         "verify support of cli version")
    ;
    return hidden;
}

// define global options
po::options_description
global_get_options_description() {
    po::options_description glob("SINA global options");
    glob.add_options()
        ("help,h",
         "print help message")
        ("version",
         "print version string")
        ("show-conf",
         "show effective configuration")

        // source (none => FAIL)
        ("in,i",
         po::value<string>(),
         "input file")
        ("intype",
         po::value<SEQUENCE_DB_TYPE>()->default_value(SEQUENCE_DB_NONE),
         "input file type")

        // sink (none => source)
        ("out,o",
         po::value<string>(),
         "output file")
        ("outtype",
         po::value<SEQUENCE_DB_TYPE>()->default_value(SEQUENCE_DB_NONE),
         "output file type")

        // pipe form
        ("no-align",
         po::bool_switch(),
         "disable alignment stage")
        ("prealigned",
         po::bool_switch(),
         "keep input alignment")
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
    po::options_description desc;

    desc.add(global_get_options_description());
    desc.add(Log::get_options_description());
    desc.add(rw_arb::get_options_description());
    desc.add(rw_fasta::get_options_description());
    desc.add(aligner::get_options_description());
    desc.add(search_filter::get_options_description());

    po::options_description all_opts(desc);
    all_opts.add(global_get_hidden_options_description());

    try {
        po::store(po::parse_command_line(argc,argv,all_opts),vm);

        if (vm.count("help")) {
            cerr << desc << endl;
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

        if (vm.count("in") == 0) {
            throw logic_error("need input file");
        }

        // Autodetect / validate intype selection
        if (vm["intype"].defaulted()) {
            const string in = vm["in"].as<string>();
            if (iends_with(in, ".arb")
                || iequals(in, ":")) {
                std::vector<string> cmd(2);
                cmd[0]="--intype";
                cmd[1]="ARB";
                po::store(po::command_line_parser(cmd).options(desc).run(),
                          vm);
            } else {
                std::vector<string> cmd(2);
                cmd[0]="--intype";
                cmd[1]="FASTA";
                po::store(po::command_line_parser(cmd).options(desc).run(),
                          vm);
            } 
        }

        // Pick suitable output if no output given
        if (vm.count("out") == 0) {
            std::vector<string> cmd(4);
            switch(vm["intype"].as<SEQUENCE_DB_TYPE>()) {
            case SEQUENCE_DB_ARB: 
                // ARB files can be used for input and output
                cmd[0]="-o";
                cmd[1]=vm["in"].as<string>();
                po::store(po::command_line_parser(cmd).options(desc).run(), vm);
                break;
            case SEQUENCE_DB_FASTA: 
                // Use "input.fasta.aligned"
                cmd[0]="-o";
                if (equals(vm["in"].as<string>(), "/dev/stdin")) {
                    cmd[1]="/dev/stdout";
                } else {
                    cmd[1]=vm["in"].as<string>()+".aligned";
                }
                po::store(po::command_line_parser(cmd).options(desc).run(), vm);
                break;
            default:
                throw logic_error("broken output type");
            }
        }

        // Autodetect / validate outtype selection
        if (vm["outtype"].defaulted()) {
            const string in = vm["out"].as<string>();
            if (iends_with(in, ".arb")
                || iequals(in, ":")) {
                std::vector<string> cmd(2);
                cmd[0]="--outtype";
                cmd[1]="ARB";
                po::store(po::command_line_parser(cmd).options(desc).run(),
                          vm);
            } else {
                std::vector<string> cmd(2);
                cmd[0]="--outtype";
                cmd[1]="FASTA";
                po::store(po::command_line_parser(cmd).options(desc).run(),
                          vm);
            }
        }
        
        // enable calc-idty if min-idty requested
        if (!vm["min-idty"].empty()) {
            std::vector<string> cmd(1);
            cmd[0]="--calc-idty";
            po::store(po::command_line_parser(cmd).options(desc).run(), vm);
        }

        Log::validate_vm(vm);
        aligner::validate_vm(vm);
        rw_arb::validate_vm(vm);
        rw_fasta::validate_vm(vm);
        search_filter::validate_vm(vm);

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

    try {
        typed_PipeElement<void,tray> *source;
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

        typed_PipeElement<tray,void> *sink;
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

        typed_PipeElement<tray, tray> *aligner;
        if (vm["no-align"].as<bool>()) {
            aligner = null_filter::make_null_filter();
        } else if (vm["prealigned"].as<bool>()) {
            aligner = copy_alignment::make_copy_alignment();
        } else {
            aligner = aligner::make_aligner();        
        }

        typed_PipeElement<tray, tray> *search;
        if (vm["search"].as<bool>()) {
            search = search_filter::make_search_filter();
        } else {
            search = null_filter::make_null_filter();
        }

        typed_PipeElement<tray, tray> *printer;
        printer = Log::make_printer();

        Pipe p = *source | *aligner | *search | *printer | *sink;

        // workaround for double reference because aligner is pss itself
        if (vm["no-align"].as<bool>() == false
            && vm["prealigned"].as<bool>() == false) {
            delete aligner;
        }


        timestamp before;
        p.run();
        timestamp after;
        cerr << "Time for alignment phase: " << after-before << "s" << endl;

        p.destroy();
  
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
