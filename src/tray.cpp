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

#include "tray.h"
#include <sstream>

namespace sina {

tray::tray() 
    : input_sequence(0), 
      aligned_sequence(0), 
      alignment_reference(0), 
      search_result(0), 
      logstream(0)
{
}

tray::tray(const tray& o) 
    : input_sequence(o.input_sequence), 
      aligned_sequence(o.aligned_sequence), 
      alignment_reference(o.alignment_reference), 
      search_result(o.search_result), 
      logstream(o.logstream)
{
}

tray&
tray::operator=(const tray& o) {
    input_sequence=o.input_sequence;
    aligned_sequence=o.aligned_sequence;
    alignment_reference=o.alignment_reference;
    search_result=o.search_result;
    logstream=o.logstream;

    return *this;
}

tray::~tray() {
}

void
tray::destroy() {
    if (input_sequence) delete input_sequence;
    if (aligned_sequence) delete aligned_sequence;
    if (alignment_reference) delete alignment_reference;
    if (search_result) delete search_result;
    if (logstream) delete logstream;
}

};


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
