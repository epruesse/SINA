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

#include "rw_csv.h"
#include "log.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
namespace bi = boost::iostreams;

namespace sina {
namespace rw_csv {

static const char* module_name = "CSV I/O";
static auto logger = Log::create_logger(module_name);

struct options {
};
static options opts;

void get_options_description(po::options_description &main,
                             po::options_description &adv) {
}

void validate_vm(po::variables_map &, po::options_description &) {
}

struct writer::priv_data {
    bi::file_descriptor_sink file;
    bi::filtering_ostream out;
    unsigned long copy_relatives;
    std::vector<std::string> v_fields;
    std::vector<std::string> headers;
};

writer::writer(const fs::path& outfile, unsigned int copy_relatives,
               std::vector<std::string>& fields)
    : data(new priv_data())
{
    data->copy_relatives = copy_relatives;
    data->v_fields = fields;

    if (outfile == "-") {
        data->file.open(STDOUT_FILENO, bi::never_close_handle);
    } else {
        data->file.open(outfile.c_str(), std::ios_base::binary);
    }
    if (!data->file.is_open()) {
        auto msg = "Unable to open file {} for writing.";
        throw std::runtime_error(fmt::format(msg, outfile));
    }
    if (outfile.extension() == ".gz") {
        data->out.push(bi::gzip_compressor());
    }
    data->out.push(data->file);
}
writer::writer(const writer&) = default;
writer& writer::operator=(const writer&) = default;
writer::~writer() = default;


std::string escape_string(const std::string& in) {
    if (in.find_first_of("\",\r\n") == std::string::npos) {
        return in;
    }
    std::stringstream tmp;
    tmp << "\"";
    size_t j = 0;
    for (size_t i = in.find('"'); i != std::string::npos;
         j=i+1, i = in.find('"',i+1)) {
        tmp << in.substr(j, i-j) << "\"\"";
    }
    tmp << in.substr(j) << "\"";
    return tmp.str();
}

static inline void append(fmt::memory_buffer& buf, const std::string& str) {
    buf.append(str.data(), str.data() + str.size());
}


tray writer::operator()(tray t) {
    const char sep[] = ",";
    const char crlf[] = "\r\n";
    const char id[] = "ID";
    fmt::memory_buffer buf;

    if (t.aligned_sequence == nullptr) {
        return t;
    }
    if (data->headers.empty()) {
        auto attrs = t.aligned_sequence->get_attrs();
        data->headers.reserve(attrs.size()+1);
        buf.append(id, id + sizeof(id));
        for (auto& ap : attrs) {
            data->headers.push_back(ap.first);
            buf.append(sep, sep + sizeof(sep));
            append(buf, ap.first);
        }
        buf.append(crlf, crlf+sizeof(crlf));
    }
    const std::string& name = t.aligned_sequence->getName();
    append(buf, name);
    for (auto& key : data->headers) {
        buf.append(sep, sep + sizeof(sep));
        append(buf, t.aligned_sequence->get_attr<std::string>(key));
    }
    buf.append(crlf, crlf+sizeof(crlf));

    fmt::internal::write(data->out, buf);

    return t;
}

} // namespace rw_csv
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
