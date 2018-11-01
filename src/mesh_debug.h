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

// this file contains functions to prettyprint the contents of
// the alignment matrix for debugging purposes

#include <iostream>
#include <fstream>
#include "mesh.h"

/* writes dot format commands to out that represent the given sequence
 * care is taken to display the nodes either vertically or horizontally
 */
template<typename SEQUENCE>
void
draw_axis(SEQUENCE &s,
	  typename SEQUENCE::iterator &begin,
	  typename SEQUENCE::iterator &end,
	  unsigned int from, unsigned int to, std::ostream& out,
	  bool horizontal) {
  const char *nname = (horizontal?"h":"v");

  typename SEQUENCE::iterator it = s.begin();
  begin = s.begin();
  end = s.end();

  // find parts to show;
  while (it != end && it->getPosition() < from) {
    ++it;
  }
  begin=it;
  while (it != end && it->getPosition() < to) {
    ++it;
  }
  end=it;

  // print node labels
  for (it = begin; it != end && it->getPosition() < to; ++it) {
    int mnode_id = get_node_id(s,it);
    out << nname << mnode_id << " [label=\""
	<< *it << "\",style=solid]; " << std::endl;
  }

  // cluster master nodes, link with edges
  out << "{ edge [style=invis]; " << std::endl;
  if (horizontal) {
    out << " rank=min;" << std:: endl;
  }
  out << "origin -> ";
  for (it = begin; it != end && it->getPosition() < to-1; ++it) {
    out << nname << get_node_id(s,it) << " -> ";
  }
  out << nname << get_node_id(s,it) << std::endl
      << "}" << std::endl;

  // show master edges
  out << "{ edge [style=solid, constraint=true]; "
      << std::endl;
  for (it = begin; it != end; ++it) {
    int mnode_id = get_node_id(s,it);

    typename SEQUENCE::pn_iterator
      pit = it.next_begin(),
      pit_end = it.next_end();

    for (; pit != pit_end; ++pit) {
      if (pit->getPosition() < to) {
	out << nname << mnode_id <<
	  " -> " << nname << get_node_id(s,pit) << std::endl;
      }
    }
  }
  out << "}" << std::endl;
}

template<typename MESH>
void
mesh_to_svg(MESH& mesh, unsigned int from, unsigned int to, std::ostream& out) {
  out << "digraph {" << std::endl
    //      << "graph [rankdir=LR]; " << std::endl
      << "node [style=invis]; " << std::endl
      << "origin [style=invis]; " << std::endl
    ;

  using midx_type = typename MESH::master_idx_type;
  using sidx_type = typename MESH::slave_idx_type;

  typename MESH::master_type::iterator  mit, mit_end, mit_begin;
  draw_axis(mesh._master, mit_begin, mit_end, from, to, out, false);

  typename MESH::slave_type::iterator   sit, sit_end, sit_begin;
  draw_axis(mesh._slave, sit_begin, sit_end, from, to, out, true);

  // print node labels
  for (mit = mit_begin; mit != mit_end; ++mit) {
    midx_type midx = get_node_id(mesh._master, mit);
    for (sit = sit_begin; sit != sit_end; ++sit) {
      sidx_type sidx = get_node_id(mesh._slave, sit);
      out << "f_" << midx << "_" << sidx << " [label=<<TABLE BORDER=\"0\""
	  << R"( CELLBORDER="1" CELLSPACING="0">)"
	  << "<TR><TD>" <<  -mesh(midx,sidx).value
	  << " (" << mesh(midx,sidx).value - 
	mesh(mesh(midx,sidx).value_midx,
	     mesh(midx,sidx).value_sidx).value << ")"
	  << "</TD></TR><TR><TD>" << -mesh(midx,sidx).gapm_val
	  << "/" << -mesh(midx,sidx).gaps_val
	  << "|" << mesh(midx,sidx).gaps_max
	  << "</TD></TR><TR><TD>" << *mit
	  << "/" << *sit
	  << "</TD></TR></TABLE>>];" << std:: endl;
    }
  }

  // link horizontally
  for (mit = mit_begin; mit != mit_end; ++mit) {
    midx_type midx = get_node_id(mesh._master, mit);
    out << "{ rank=same; edge [style=invis]; v" << midx;
    for (sit = sit_begin; sit != sit_end; ++sit) {
      sidx_type sidx = get_node_id(mesh._slave, sit);
      out << " -> " << "f_" << midx << "_" << sidx;
    }
    out << "}" << std::endl;
  }

  // link vertically
  for (sit = sit_begin; sit != sit_end; ++sit) {
    sidx_type sidx = get_node_id(mesh._slave, sit);
    out << "{ edge [style=invis]; h" << sidx;
    for (mit = mit_begin; mit != mit_end; ++mit) {
      midx_type midx = get_node_id(mesh._master, mit);

      out << " -> " << "f_" << midx << "_" << sidx;
    }
    out << "}" << std::endl;
  }

  out << "edge [style=solid,constraint=true]; " << endl;

  for (sit = sit_begin; sit != sit_end; ++sit) {
    sidx_type sidx = get_node_id(mesh._slave, sit);
    for (mit = mit_begin; mit != mit_end; ++mit) {
      midx_type midx = get_node_id(mesh._master, mit);
      midx_type tgt_midx = mesh(midx,sidx).value_midx;
      sidx_type tgt_sidx = mesh(midx,sidx).value_sidx;
      if (mesh._master.getById(tgt_midx).getPosition() >= from  &&
	  mesh._slave.getById(tgt_sidx).getPosition() >= from) {
	out << "f_" << midx << "_" << sidx
	    << " -> "
	    << "f_" << mesh(midx,sidx).value_midx
	     << "_" <<  mesh(midx,sidx).value_sidx
	     << ";" << endl;
	}
    }
  }
  out << "}" << std::endl;
}


template<typename MESH>
void
mesh_to_svg(MESH& mesh, unsigned int from, unsigned int to, const char* outfile) {
  std::ofstream out(outfile);
  if (out.bad()) {
    std::cerr << "unable to open file " << outfile << std::endl;
    return;
  }

  mesh_to_svg(mesh, from, to, out);
}




struct default_weight {
    int operator[](int) { return 1; }
};

template<typename L, typename R>
std::pair<float,float>
seq_compare(L& left, R& right) {
    default_weight w;
    return seq_compare(left,right,w);
}


template<typename L, typename R, typename W>
std::pair<float,float>
seq_compare(L& left, R& right, W& weight) {
    typename L::iterator l_it = left.begin(), l_end = left.end();
    typename R::iterator r_it = right.begin(), r_end = right.end();

    float mismatches = 0;
    float matches = 0;

    if (left.size() != right.size()) {
        std::cout << "seq_compare: length differs by "
                  << left.size()-right.size() << std::endl;
    }
    if (left.getBases() != right.getBases()) {
        std::cout << "seq_compare: bases differ" << std::endl;
        for(;l_it != l_end && r_it != r_end; ++l_it, ++r_it) {
            if(l_it->getBase() != r_it->getBase()) {
                std::cout << *l_it << " " << *r_it << std::endl;
            }
        }
    }

    while(l_it != l_end && r_it != r_end) {
        int lpos = l_it->getPosition();
        int rpos = r_it->getPosition();
        if (lpos < rpos) {
            mismatches += weight[rpos];
            ++l_it;
        } else if ( rpos < lpos ) {
            mismatches += weight[rpos];
            ++r_it;
        } else { // rpos <=> lops
            if (l_it->getBase() != r_it->getBase()) {
                mismatches += weight[rpos];
            } else {
                matches += weight[rpos];
	    }
            ++r_it;
            ++l_it;
        }
    }
    return std::make_pair(mismatches,matches);
}



/*

  digraph G {
  graph [rankdir=TD]
  edge [style=invis];
  o [style=invis];
  subgraph cluserHead{
      rank=min;
      o x1 x2  x3  x4  x5
  }

  subgraph cluserLeft{
     o ->  y1 -> y2 ->  y3 -> y4 -> y5
  }

  f11 [label=<<TABLE BORDER="0"><TR><TD COLSPAN="2">6</TD></TR><TR><TD>3</TD><TD>2</TD></TR></TABLE>>]; 

  { rank=same; y1 -> f11 -> f21 -> f31 -> f41 -> f51 }
  { rank=same; y2 -> f12 -> f22 -> f32 -> f42 -> f52 }
  { rank=same; y3 -> f13 -> f23 -> f33 -> f43 -> f53 }
  { rank=same; y4 -> f14 -> f24 -> f34 -> f44 -> f54 }
  { rank=same; y5 -> f15 -> f25 -> f35 -> f45 -> f55 }

  { x1 -> f11 -> f12 -> f13 -> f14 -> f15 }
  { x2 -> f21 -> f22 -> f23 -> f24 -> f25 }
  { x3 -> f31 -> f32 -> f33 -> f34 -> f35 }
  { x4 -> f41 -> f42 -> f43 -> f44 -> f45 }
  { x5 -> f51 -> f52 -> f53 -> f54 -> f55 }

  edge [style=solid];
  { edge [constraint=false];
    f22 -> f11;
    f45 -> f22;
}
*/
