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


#ifndef _ALIGN_H_
#define _ALIGN_H_

#include "query_arb.h"
#include "query_pt.h"
#include "pipe.h"
#include "tray.h"

#include <vector>
#include <boost/program_options.hpp>

namespace sina {

enum TURN_TYPE {
    TURN_NONE=0,
    TURN_REVCOMP=1,
    TURN_ALL=2
};

enum OVERHANG_TYPE {
    OVERHANG_ATTACH,
    OVERHANG_REMOVE,
    OVERHANG_EDGE
};

enum LOWERCASE_TYPE {
    LOWERCASE_NONE,
    LOWERCASE_ORIGINAL,
    LOWERCASE_UNALIGNED
};

enum INSERTION_TYPE {
  INSERTION_SHIFT,
  INSERTION_FORBID,
  INSERTION_REMOVE
};

class aligner {
    query_arb *arb;

    std::vector<float> weights;
    std::vector<int> pairs;

public:
    class famfinder;
    class galigner;

    static typed_PipeElement<tray,tray>* make_aligner();

    static const char* fn_turn;
    static const char* fn_acc;
    static const char* fn_start;
    static const char* fn_qual;
    static const char* fn_head;
    static const char* fn_tail;
    static const char* fn_date;
    static const char* fn_astart;
    static const char* fn_astop;
    static const char* fn_idty;
    static const char* fn_family;
    static const char* fn_family_str;
    static const char* fn_nuc;
    static const char* fn_nuc_gene;
    static const char* fn_bpscore;
    static const char* fn_used_rels;
    static const char* fn_fullname;

    static boost::program_options::options_description get_options_description();
    static void validate_vm(boost::program_options::variables_map&);
private:
    struct options;
    static struct options *opts;
};


std::ostream& operator<<(std::ostream&, const sina::TURN_TYPE&);
std::ostream& operator<<(std::ostream&, const sina::OVERHANG_TYPE&);
std::ostream& operator<<(std::ostream&, const sina::LOWERCASE_TYPE&);
std::ostream& operator<<(std::ostream&, const sina::INSERTION_TYPE&);

void validate(boost::any&, const std::vector<std::string>&,
              sina::TURN_TYPE*,int);
void validate(boost::any&, const std::vector<std::string>&,
              sina::OVERHANG_TYPE*,int);
void validate(boost::any&, const std::vector<std::string>&,
              sina::LOWERCASE_TYPE*,int);
void validate(boost::any&, const std::vector<std::string>&,
              sina::INSERTION_TYPE*,int);


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
