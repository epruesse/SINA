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

#ifndef _NULL_FILTER_H_
#define _NULL_FILTER_H_

#include "pipe.h"
#include "cseq.h"

namespace sina {

class null_filter 
  : public PipeElement<tray,tray>
{
public:
  static typed_PipeElement<tray,tray>* make_null_filter() {
        return new null_filter();
    }
    tray operator()(tray c) {
        return c;
    }
    std::string getName() const { return "null-filter"; }
    null_filter() {};  
    ~null_filter() {};
private:
};

} // namespace sina


#endif // _SEARCH_H_
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
