/*
Copyright (c) 2006-2017 Elmar Pruesse <elmar.pruesse@ucdenver.edu>

This file is part of SINA.
123456789012345678901234567890123456789012345678901234567890123456789012
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

#include "graph.h"

#include <iostream>

using namespace std;
struct aligned_base {
    aligned_base(int pos=0, char base='-')
        : _position(pos), _base(base) {}
    int _position;
    char _base;
};

ostream& operator<<(ostream& out, aligned_base& b) {
    out << b._base << "(" << b._position << ")";
}

int main(int argc, char** argv) {
    typedef dag<aligned_base> tdag;
    tdag mydag;
    tdag::node_ref a,b;

    a = mydag.insert(aligned_base(1,'G'));
    b = mydag.insert(aligned_base(2,'A'));
    mydag.link(a,b);
    a = mydag.insert(aligned_base(3,'T'));
    mydag.link(b,a);
    a = mydag.insert(aligned_base(3,'C'));
    mydag.link(b,a);

    mydag.print_graphviz(cout,"test");
}

/*
  Local Variables:
  mode:c++
  c-file-style:"stroustrup"
  c-file-offsets:((innamespace . 0)(inline-open . 0)(case-label . +))
  indent-tabs-mode:nil
  fill-column:99
  End:
*/
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:encoding=utf-8:textwidth=99 :

