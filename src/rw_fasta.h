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

#ifndef _RW_FASTA_H_
#define _RW_FASTA_H_

#include "tray.h"

#include <fstream>
#include <memory>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace sina {

enum FASTA_META_TYPE {
    FASTA_META_NONE=0,
    FASTA_META_HEADER=1,
    FASTA_META_COMMENT=2,
    FASTA_META_CSV=3
};
std::ostream& operator<<(std::ostream&, const sina::FASTA_META_TYPE&);
void validate(boost::any&, const std::vector<std::string>&,
              sina::FASTA_META_TYPE*, int);

class rw_fasta {
public:
    class reader {
        struct priv_data;
        std::shared_ptr<priv_data> data;
    public:
        explicit reader(const boost::filesystem::path&);
        reader(const reader&);
        reader& operator=(const reader&);
        ~reader();
        bool operator()(tray&);
    };

    class writer {
        struct priv_data;
        std::shared_ptr<priv_data> data;
    public:
        explicit writer(const boost::filesystem::path&);
        writer(const writer&);
        writer& operator=(const writer&);
        ~writer();
        tray operator()(tray);
    };

    static void get_options_description(boost::program_options::options_description& all,
                                        boost::program_options::options_description& adv);
    static void validate_vm(boost::program_options::variables_map&,
                            boost::program_options::options_description&);

private:
    struct options;
    static struct options *opts;
};

} // namespace sina

#endif // _RW_FASTA_H

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
