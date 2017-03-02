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

#include "pseq.h"

#include <string>

namespace sina {
std::ostream&
operator<<(std::ostream& out, const aligned_base_profile& ab)
{
    std::stringstream tmp;
    tmp <<  ab.getBase().getString() << "(" << ab.getPosition() << ")";
    return out << tmp.str();
}
}

using namespace sina;
pseq::pseq(std::vector<cseq>::iterator seqs_it,
           std::vector<cseq>::iterator seqs_end) {

    width = seqs_it->getWidth();
    int height = 0;

    std::vector<cseq::iterator> base_iterators, base_ends;
    for (; seqs_it != seqs_end; ++seqs_it) {
        base_iterators.push_back(seqs_it->begin());
        base_ends.push_back(seqs_it->end());
        ++height;
    }

    bool gap[height];
    for (int i=0; i<height; ++i) {
        gap[i] = true;
    }

    typedef aligned_base::idx_type idx_type;

    idx_type current_column = 0;
    while (current_column < width) {
        idx_type next_column = width;
        int A=0, G=0, C=0, T=0;
        int gapOpen=0, gapExtend=0;
        for (int row = 0; row < height; ++row) {
            cseq::iterator it = base_iterators[row];
            if (it != base_ends[row] &&
                it->getPosition() == current_column) {
                if (it->ambig_order() > 0) {
                    int points = 12 / it->ambig_order();
                    if (it->has_A()) A+= points;
                    if (it->has_G()) G+= points;
                    if (it->has_C()) C+= points;
                    if (it->has_TU()) T+= points;
                    gap[row] = false;
                }
                ++base_iterators[row];
            } else {
                if (gap[row]) {
                    ++gapExtend;
                } else {
                    gap[row] = true;
                    ++gapOpen;
                }
            }

            if (base_iterators[row] != base_ends[row]) {
                next_column = std::min(
                    next_column,
                    base_iterators[row]->getPosition()
                    );
            }
        }

        base_profile bp(A, G, C, T, gapOpen*12, gapExtend*12);
        profile.push_back(aligned_base_profile(current_column, bp));

        current_column = next_column;
    }
}

void
pseq::print_graphviz(std::ostream& out, std::string /*name*/) {
    for (iterator it = begin(); it != end(); ++it) {
        out << it->getBase().getString() << std::endl;
    }
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

