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


#ifndef _ALIGN_H_
#define _ALIGN_H_

#include "tray.h"

#include <vector>
#include <boost/program_options.hpp>

namespace sina {

/** Parameter Types 
 **/
enum OVERHANG_TYPE {
    OVERHANG_ATTACH,
    OVERHANG_REMOVE,
    OVERHANG_EDGE
};
std::ostream& operator<<(std::ostream& out, const sina::OVERHANG_TYPE& t);
void validate(boost::any& v, const std::vector<std::string>& values,
              sina::OVERHANG_TYPE* /*unused*/, int /*unused*/);

enum LOWERCASE_TYPE {
    LOWERCASE_NONE,
    LOWERCASE_ORIGINAL,
    LOWERCASE_UNALIGNED
};
std::ostream& operator<<(std::ostream& out, const sina::LOWERCASE_TYPE& t);
void validate(boost::any& v, const std::vector<std::string>& values,
              sina::LOWERCASE_TYPE* /*unused*/, int /*unused*/);

enum INSERTION_TYPE {
  INSERTION_SHIFT,
  INSERTION_FORBID,
  INSERTION_REMOVE
};
std::ostream& operator<<(std::ostream& out, const sina::INSERTION_TYPE& t);
void validate(boost::any& v, const std::vector<std::string>& values,
              sina::INSERTION_TYPE* /*unused*/, int /*unused*/);


class aligner {
    struct options;
public:
    static struct options *opts;
    aligner();
    aligner(const aligner& rhs);
    ~aligner();
    aligner& operator=(const aligner& rhs);
    tray operator()(tray t);

    static void get_options_description(boost::program_options::options_description& all,
                                        boost::program_options::options_description& adv);
    static void validate_vm(boost::program_options::variables_map& /*unused*/,
                            boost::program_options::options_description& /*unused*/);
};


} // namespace sina


#endif // _ALIGN_H

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
