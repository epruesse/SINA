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

#ifndef _LOG_H_
#define _LOG_H_

#include <boost/program_options.hpp>
#include "spdlog/spdlog.h"
#include "spdlog/fmt/ostr.h" // logging from ostreamable

#include "tray.h"

namespace sina {

class Log {
private:
    struct options;
    static std::unique_ptr<options> opts;
public:
    class printer {
        struct priv_data;
        std::shared_ptr<priv_data> data;
    public:
        printer();
        ~printer();
        printer(const printer& /*o*/);
        printer& operator=(const printer& /*o*/);
        tray operator()(tray /*t*/);
    };

    static std::shared_ptr<spdlog::logger> create_logger(std::string name);
    static void get_options_description(boost::program_options::options_description& main,
                                        boost::program_options::options_description& adv);
    static void validate_vm(boost::program_options::variables_map& /*vm*/,
                            boost::program_options::options_description& /*unused*/);

};

} // namespace sina


#endif // _LOG_H_
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
