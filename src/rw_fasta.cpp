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

#include "rw_fasta.h"

#include <fstream>
#include <utility>
#include <vector>
#include <iostream>
#include <unordered_set>
#include <string>
#include <sstream>

#include <functional>

#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/counter.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>

#include "query_arb.h"
#include "log.h"
#include "progress.h"

using std::stringstream;
using std::vector;
using std::string;
using boost::algorithm::iequals;

namespace bi = boost::iostreams;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

namespace sina {

static const char* module_name = "FASTA I/O";
static auto logger = Log::create_logger(module_name);

// define extra datatype for metadata output format
std::ostream& operator<<(std::ostream& out, const sina::FASTA_META_TYPE& m) {
    switch(m) {
    case FASTA_META_NONE:
        out << "none";
    break;
    case FASTA_META_HEADER:
        out << "header";
    break;
    case FASTA_META_COMMENT:
        out << "comment";
    break;
    case FASTA_META_CSV:
        out << "csv";
    break;
    default:
        out << "[UNKNOWN!] (value=" << (int)m << ")";
    }
    return out;
}
void validate(boost::any& v,
              const std::vector<std::string>& values,
              sina::FASTA_META_TYPE* /*m*/, int /*unused*/) {
    po::validators::check_first_occurrence(v);
    const std::string &s = po::validators::get_single_string(values);
    if (iequals(s, "none")) {
        v = FASTA_META_NONE;
    } else if (iequals(s, "header")) {
        v = FASTA_META_HEADER;
    } else if (iequals(s, "comment")) {
        v = FASTA_META_COMMENT;
    } else if (iequals(s, "csv")) {
        v = FASTA_META_CSV;
    } else {
        throw po::invalid_option_value("must be one of 'none', 'header', 'comment' or 'cvs'");
    }
}

// struct holding configuration options
struct rw_fasta::options {
    FASTA_META_TYPE fastameta;
    int line_length;
    float min_idty;
    long fasta_block;
    long fasta_idx;
    bool out_dots;
    bool out_dna;
};
struct rw_fasta::options *rw_fasta::opts;

void
rw_fasta::get_options_description(po::options_description& main,
                                  po::options_description& adv) {
    opts = new options;

    main.add_options()
        ("meta-fmt",
         po::value<FASTA_META_TYPE>(&opts->fastameta)->default_value(FASTA_META_NONE,""),
         "meta data in (*none*|header|comment|csv)")
        ;

    po::options_description od(module_name);
    od.add_options()
        // write
        ("line-length",
         po::value<int>(&opts->line_length)->default_value(0, ""),
         "wrap output sequence (unlimited)")

        ("min-idty",
         po::value<float>(&opts->min_idty)->default_value(0.f, ""),
         "only write sequences with align_idty_slv > X, implies calc-idty")
        ("fasta-write-dna",
         po::bool_switch(&opts->out_dna),
         "Write DNA sequences (default: RNA)")
        ("fasta-write-dots",
         po::bool_switch(&opts->out_dots),
         "Use dots instead of dashes to distinguish unknown sequence data from indels")

        // read
        ("fasta-idx",
         po::value<long>(&opts->fasta_idx)->default_value(0, ""),
         "process only sequences beginning in block <arg>")
        ("fasta-block",
         po::value<long>(&opts->fasta_block)->default_value(0, ""),
         "length of blocks")
        ;
    adv.add(od);
}


void
rw_fasta::validate_vm(po::variables_map& /*vm*/,
                      po::options_description& /*desc*/) {
}


////////// reader

struct rw_fasta::reader::priv_data {
    bi::file_descriptor_source file;
    bi::filtering_istream in;
    bi::counter *counter{0};
    fs::path filename;
    size_t file_size{0};
    int lineno{0};
    int seqno{0};

    vector<string>& v_fields;
    logger_progress* p{nullptr};

    priv_data(fs::path filename_, vector<string>& fields)
        : filename(std::move(filename_)),
          v_fields(fields)
    {}
    ~priv_data() {
        logger->info("read {} sequences from {} lines", seqno-1, lineno-1);
    }
};

rw_fasta::reader::reader(const fs::path& infile, vector<string>& fields)
    : data(new priv_data(infile, fields))
{
    if (infile == "-") {
        data->file.open(STDIN_FILENO, bi::never_close_handle);
    } else {
        data->file.open(infile.c_str(), std::ios_base::binary);
        if (fs::is_regular_file(infile)) {
            data->file_size = fs::file_size(infile);
        }
    }
    if (!data->file.is_open()) {
        stringstream msg; 
        msg << "Unable to open file \"" << infile << "\" for reading.";
        throw std::runtime_error(msg.str());
    }
    if (infile.extension() == ".gz") {
        data->in.push(bi::gzip_decompressor());
    }
    if (data->file_size > 0) {
        data->in.push(bi::counter());
        data->counter = data->in.component<bi::counter>(data->in.size()-1);
    }
    data->in.push(data->file);

    // if fasta blocking enabled, seek to selected block
    if (opts->fasta_block > 0) {
        if (infile == "-") {
            throw std::logic_error("Cannot use --fasta-idx when input is piped");
        }

        data->in.seekg(opts->fasta_block * opts->fasta_idx);
    }
}

rw_fasta::reader::reader(const reader& r) = default;
rw_fasta::reader& rw_fasta::reader::operator=(const reader& r) = default;
rw_fasta::reader::~reader() = default;

void
rw_fasta::reader::set_progress(logger_progress& p) {
    data->p = &p;
}

bool
rw_fasta::reader::operator()(tray& t) {
    t.seqno = ++data->seqno;
    t.input_sequence = new cseq();
    cseq &c = *t.input_sequence;
    if (data->in.fail()) {
        return false;
    }
    
    // if fasta blocking enabled, check if we've passed block
    // boundary in last sequence
    if (opts->fasta_block > 0
        && data->in.tellg() > opts->fasta_block * (opts->fasta_idx + 1)) {
        return false;
    }

    string line;
    // skip lines not beginning with '>'
    while (data->in.peek() != '>' && getline(data->in, line).good()) {
        data->lineno++;
    }

    // parse title
    data->lineno++;
    if (getline(data->in, line).good()) {
        if (line[line.size()-1] == '\r') {  // "\r" at end
            line.resize(line.size()-1);
        }
        
        // set name to text between first '>' and first ' '
        unsigned int blank = line.find_first_of(" \t");
        if (blank == 0) {
            blank = line.size();
        }
        c.setName(line.substr(1, blank-1));
        if (blank < line.size()) {
            c.set_attr<string>(query_arb::fn_fullname, line.substr(blank+1));
        }
    } else { // didn't get a title
        return false;
    }

    // handle comments
    while (data->in.peek() == ';' && getline(data->in, line).good()) {
        data->lineno++;

        // if the comment contains an attribute: add it. 
        // Otherwise ignore the comment

        size_t equalsign = line.find_first_of('=');
        if (equalsign != string::npos) {
            string key = line.substr(1, equalsign-1);
            boost::trim(key);
            string value = line.substr(equalsign+1);
            boost::trim(value);
            c.set_attr(key,value);
        }
    }

    try {
        // all lines until eof or next /^>/ are data
        while (data->in.peek() != '>' && data->in.good()) {
            getline(data->in, line);
            data->lineno++;
            c.append(line);
        }
    } catch (base_iupac::bad_character_exception& e) {
        logger->error("Skipping sequence {} (>{}) at {}:{} (contains character '{}')",
                      data->seqno, c.getName(),
                      data->filename, data->lineno,
                      e.character);
        while (data->in.peek() != '>' && getline(data->in, line).good()) {
            data->lineno++;
        }
        delete t.input_sequence;
        return (*this)(t); // FIXME: stack size?
    }

    if (data->p && data->counter) {
        auto bytes_read = data->counter->characters();
        if (bytes_read) {
            data->p->set_total(data->seqno * data->file_size / bytes_read);
        }
    }

    logger->debug("loaded sequence {}", t.input_sequence->getName());
    return true;
}


/////////////// writer

struct rw_fasta::writer::priv_data {
    bi::file_descriptor_sink file;
    bi::filtering_ostream out;
    std::ofstream out_csv;
    int count;
    int excluded;
    std::unordered_set<string> relatives_written;
    unsigned long copy_relatives;
    vector<string>& v_fields;
    priv_data(unsigned int copy_relatives_,
              vector<string>& fields)
        : count(0),
          excluded(0),
          copy_relatives(copy_relatives_),
          v_fields(fields)
    {}
    ~priv_data() {
        logger->info("wrote {} sequences ({} excluded, {} relatives)",
                     count, excluded, relatives_written.size());
    }
    void write(const cseq& c);
};

rw_fasta::writer::writer(const fs::path& outfile,
                         unsigned int copy_relatives,
                         vector<string>& fields)
    : data(new priv_data(copy_relatives, fields))
{
    if (outfile == "-") {
        data->file.open(STDOUT_FILENO, bi::never_close_handle);
    } else {
        data->file.open(outfile.c_str(), std::ios_base::binary);
    }
    if (!data->file.is_open()) {
        stringstream msg; 
        msg << "Unable to open file \"" << outfile << "\" for writing.";
        throw std::runtime_error(msg.str());
    }
    if (outfile.extension() == ".gz") {
        data->out.push(bi::gzip_compressor());
    }
    data->out.push(data->file);

    if (opts->fastameta == FASTA_META_CSV) {
        fs::path outcsv(outfile);
        outcsv.replace_extension(".csv");
        data->out_csv.open(outcsv.native());
        if (data->out_csv.fail()) {
            stringstream msg; 
            msg << "Unable to open file \"" << outfile << ".csv\" for writing.";
            throw std::runtime_error(msg.str());
        }
    }
}

rw_fasta::writer::writer(const writer& o) = default;
rw_fasta::writer& rw_fasta::writer::operator=(const writer& o) = default;
rw_fasta::writer::~writer() = default;

string escape_string(const string& in) {
    if (in.find_first_of("\",\r\n") == string::npos) {
        return in;
    }
    stringstream tmp;
    tmp << "\"";
    size_t j = 0;
    for (size_t i = in.find('"'); i != string::npos; 
         j=i+1, i = in.find('"',i+1)) {
        tmp << in.substr(j, i-j) << "\"\"";
    }
    tmp << in.substr(j) << "\"";
    return tmp.str();
}

tray
rw_fasta::writer::operator()(tray t) {
    if (t.input_sequence == nullptr) {
        throw std::runtime_error("Received broken tray in " __FILE__);
    }
    if (t.aligned_sequence == nullptr) {
        logger->info("Not writing sequence {} (>{}): not aligned",
                     t.seqno, t.input_sequence->getName());
        ++data->excluded;
        return t;
    }
    if (opts->min_idty > 0) {
        auto idty = t.aligned_sequence->get_attr<float>(query_arb::fn_idty, 0);
        if (opts->min_idty > idty) {
            logger->info("Not writing sequence {} (>{}): below identity threshold ({}<={})",
                         t.seqno, t.input_sequence->getName(),
                         idty, opts->min_idty);
            ++data->excluded;
            return t;
        }
    }
    cseq &c = *t.aligned_sequence;

    data->write(c);

    if (data->copy_relatives != 0u) {
        auto *relatives = t.search_result != nullptr ? t.search_result : t.alignment_reference;
        if (relatives != nullptr) {
            int i = data->copy_relatives;
            for (auto& item : *relatives) {
                if (data->relatives_written.insert( item.sequence->getName() ).second) {
                    data->write(*item.sequence);
                }
                if (--i == 0) {
                    break;
                }
            }
        }
    }

    return t;
}

void
rw_fasta::writer::priv_data::write(const cseq& c) {
    const auto& attrs = c.get_attrs();

    out << ">" << c.getName();
    string fname = c.get_attr<string>(query_arb::fn_fullname, "");
    if (!fname.empty()) {
        out << " " << fname;
    }

    switch (opts->fastameta) {
    case FASTA_META_NONE:
        out << std::endl;
        break;
    case FASTA_META_HEADER:
        for (auto& ap: attrs) {
            if (ap.first == query_arb::fn_family) {
                continue; // alignment family is too much
            }
            if (ap.first == query_arb::fn_fullname) {
                continue; // already written as description in header
            }
            string val = boost::apply_visitor(lexical_cast_visitor<string>(),
                                              ap.second);
            if (!val.empty()) {
                out << " [" << ap.first << "="
                    << val
                    << "]";
            }
        }
        out << std::endl;
        break;
    case FASTA_META_COMMENT:
        out << std::endl;

        for (auto& ap: attrs) {
            if (ap.first == query_arb::fn_family) {
                continue; // alignment family is too much
            }
            if (ap.first == query_arb::fn_fullname) {
                continue; // already written as description in header
            }
            out << "; " << ap.first << "="
                << boost::apply_visitor(lexical_cast_visitor<string>(),
                                        ap.second)
                << std::endl;
        }
        break;
    case FASTA_META_CSV:
        out << std::endl;

        // print header
        if (count == 0) {
            out_csv << "name";
            for (auto& ap: attrs) {
              if (ap.first == query_arb::fn_family) {
                  continue; // alignment family is too much
              }
              out_csv << "," << escape_string(ap.first);
            }
            out_csv << "\r\n";
        }

        out_csv << c.getName();
        for (auto& ap: attrs) {
            if (ap.first == query_arb::fn_family) {
                continue;
            }
            out_csv << ","
                    << escape_string(
                        boost::apply_visitor(
                            lexical_cast_visitor<string>(),
                            ap.second
                            )
                        );
        }

        out_csv << "\r\n";
        break;
    default:
        throw std::runtime_error("Unknown meta-fmt output option");
    }

    string seq  = c.getAligned(!opts->out_dots, opts->out_dna);
    int len = seq.size();
    if (opts->line_length > 0) {
        for (int i=0; i<len; i+=opts->line_length) {
            out << seq.substr(i, opts->line_length) << std::endl;
        }
    } else {
        out << seq << std::endl;
    }
    count++;
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
