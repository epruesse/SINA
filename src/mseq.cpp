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

#include "mseq.h"
#include <set>
#include <algorithm>
#include <limits>
#include "timer.h"
#include <cmath>

using std::vector;
using std::list;
using std::cerr;
using std::endl;
using std::set;
using std::min;

#ifndef TEST

using namespace sina;

/* turns a number of sequences into a graph */
mseq::mseq(vector<cseq>::iterator seqs_begin,
           vector<cseq>::iterator seqs_end,
           float weight)
    : num_seqs(seqs_end-seqs_begin), bases_width(seqs_begin->getWidth())
{
    // Sanity check input
    for (vector<cseq>::iterator it = seqs_begin;
         it != seqs_end; ++it) {
        if (bases_width != it->getWidth()) {
            cerr <<
                "Aligned sequences to be stored in mseq of differ in length!"
                 << endl;
            cerr << "Sequence " << it->getName() << "(" << it-seqs_begin  
                 << "/" << seqs_end - seqs_begin << ")"
                 << ": length = " << it->getWidth() << " expected " << 
                bases_width << endl;
            exit(1);
        }
    }

    vector<cseq::iterator> citv;
    citv.reserve(num_seqs);
    for (vector<cseq>::iterator seq_it = seqs_begin;
         seq_it != seqs_end; ++seq_it) {
        citv.push_back(seq_it->begin());
    }
    aligned_base::idx_type min_next=0;
    vector<iterator>::size_type nodes_size = std::numeric_limits<value_type::base_type>().max();
    vector<iterator> nodes(nodes_size);

    vector<iterator> last(num_seqs);
    for (unsigned int i=0; i < bases_width; i++) {
        if (min_next > i) continue;
        min_next = std::numeric_limits<int>().max();
        nodes.assign(nodes_size,iterator());

        for (unsigned int j = 0; j < num_seqs; j++) {
            if ( citv[j]->getPosition() == i) {
                iterator newnode;
                unsigned char base = citv[j]->getBase();
                if (nodes[base].isNull()) {
                    nodes[base] = newnode = insert(*citv[j]);
                } else {
                    newnode = nodes[base];
                    newnode->weight += 1.f;
                }
                if (!last[j].isNull())
                    link(last[j],newnode);
                last[j]=newnode;

                ++citv[j];
            }
            min_next=min(min_next,citv[j]->getPosition());
        }

        for (vector<iterator>::iterator it = nodes.begin();
             it != nodes.end(); ++it) {
            if (!it->isNull()) {
                (*it)->weight = 1.0/(weight+1) + weight * ((*it)->weight/num_seqs);
                    //std::min(20.0, -log(1 - (*it)->weight/num_seqs));
            }
        }
    }
}

#else // TEST
int main(int, char**) {
    std::cout << "mseq: no tests yet" << endl;
};
#endif
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

