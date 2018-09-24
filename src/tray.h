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

#ifndef _SINA_H_
#define _SINA_H_

#include "cseq.h"
#include "alignment_stats.h"
#include <boost/shared_ptr.hpp>
#include <iostream>

namespace sina {


class tray {
public:
    cseq *input_sequence;
    cseq *aligned_sequence;
    std::vector<cseq> *alignment_reference;
    std::vector<cseq> *search_result;
    std::stringstream *logstream;
    alignment_stats *astats;
    
    std::ostream& log() { 
        if (logstream) 
            return *logstream; 
        else 
            return std::cerr; 
    }
    
    tray();
    tray(const tray&);
    tray& operator=(const tray&);
    ~tray();
    void destroy();
};


} // namespace sina

#endif // _SINA_H_

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
