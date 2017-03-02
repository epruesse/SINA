/*
Copyright (c) 2006-2017 Elmar Pruesse <elmar.pruesse@ucdenver.edu>

This file is part of SINA.

SINA is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _RW_FASTA_H_
#define _RW_FASTA_H_

#include "pipe.h"
#include "tray.h"

#include <fstream>
#include <boost/program_options.hpp>

namespace sina {

class query_arb;

enum FASTA_META_TYPE {
    FASTA_META_NONE=0,
    FASTA_META_HEADER=1,
    FASTA_META_COMMENT=2,
    FASTA_META_CSV=3
};

class rw_fasta {
public:
    class reader;
    class writer;

    static typed_PipeElement<void,tray>* make_reader(boost::program_options::variables_map&);
    static typed_PipeElement<tray,void>* make_writer(boost::program_options::variables_map&);

    static boost::program_options::options_description get_options_description();
    static void validate_vm(boost::program_options::variables_map&);    
private:
    struct options;
    static struct options *opts;
};

std::ostream& operator<<(std::ostream&, const sina::FASTA_META_TYPE&);
void validate(boost::any&, const std::vector<std::string>&,
              sina::FASTA_META_TYPE*,int);


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
