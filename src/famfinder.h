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

#ifndef _FAMFINDER_H_
#define _FAMFINDER_H_

#include "tray.h"

#include <memory>
#include <boost/program_options.hpp>


namespace sina {

enum TURN_TYPE {
    TURN_NONE=0,
    TURN_REVCOMP=1,
    TURN_ALL=2
};
std::ostream& operator<<(std::ostream& out, const sina::TURN_TYPE& t);
void validate(boost::any& v, const std::vector<std::string>& values,
              sina::TURN_TYPE* /*unused*/,int /*unused*/);


class famfinder {
private:
    struct options;
    static struct options *opts;
    std::vector<float> weights;
public:
    class _famfinder;

    class finder {
        std::shared_ptr<_famfinder> data;
    public:
        finder() = default;
        explicit finder(int n);
        finder(const finder& o);
        finder& operator=(const finder& o);
        ~finder();
        tray operator()(tray t);
    };

    static void get_options_description(boost::program_options::options_description& main,
                                        boost::program_options::options_description& adv);
    static void validate_vm(boost::program_options::variables_map& vm,
                            boost::program_options::options_description& desc);
};


} // namespace sina



#endif // _FAMFINDER_H_

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
