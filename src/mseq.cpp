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

#include "log.h"
static const char* module_name = "mseq";
static auto logger = sina::Log::create_logger(module_name);

#include "mseq.h"
#include "timer.h"
#include "log.h"
#include <set>
#include <algorithm>
#include <limits>
#include <cmath>

using std::vector;
using std::list;
using std::set;
using std::min;


using namespace sina;


/* turns a number of sequences into a graph */
mseq::mseq(std::vector<cseq>::iterator seqs_begin,
           std::vector<cseq>::iterator seqs_end,
           float weight)
    : num_seqs(seqs_end-seqs_begin), bases_width(seqs_begin->getWidth())
{
    // Sanity check input
    for (vector<cseq>::iterator it = seqs_begin;
         it != seqs_end; ++it) {
        if (bases_width != it->getWidth()) {
            logger->critical("Sequence {} ({}/{}): length = {} expected {}",
                             it->getName(), it-seqs_begin, seqs_end - seqs_begin,
                             it->getWidth(), bases_width);
            throw std::runtime_error(
                "Aligned sequences to be stored in mseq of differ in length!"
                );
        }
    }

    vector<cseq::iterator> citv;
    citv.reserve(num_seqs);
    vector<cseq::iterator> citv_end;
    citv_end.reserve(num_seqs);
    for (vector<cseq>::iterator seq_it = seqs_begin;
         seq_it != seqs_end; ++seq_it) {
        citv.push_back(seq_it->begin());
        citv_end.push_back(seq_it->end());
    }
    aligned_base::idx_type min_next=0;
    vector<iterator>::size_type nodes_size = std::numeric_limits<value_type::base_type>().max();
    vector<iterator> nodes(nodes_size);

    vector<iterator> last(num_seqs);
    // iterate over alignment columns left to right
    for (unsigned int i=0; i < bases_width; i++) {
        if (min_next > i) continue;
        min_next = std::numeric_limits<int>().max();
        nodes.assign(nodes_size,iterator());

        // check all sequences for that column
        for (unsigned int j = 0; j < num_seqs; j++) {
            if (citv[j] != citv_end[j] && citv[j]->getPosition() == i) {
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
            if (citv[j] != citv_end[j]) {
                min_next=min(min_next,citv[j]->getPosition());
            }
        }

        for (auto & node : nodes) {
            if (!node.isNull()) {
                node->weight = 1.0/(weight+1) + weight * (node->weight/num_seqs);
                    //std::min(20.0, -log(1 - (*it)->weight/num_seqs));
            }
        }
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

