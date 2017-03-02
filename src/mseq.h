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

#ifndef _MSEQ_H_
#define _MSEQ_H_

#include "graph.h"
#include "cseq.h"

#include <vector>

namespace sina {

class mseq;

class mseq_node : public aligned_base 
{
  mseq_node();
public:
  mseq_node(const aligned_base& b) 
    : aligned_base(b), weight(1.f) 
  {}
  
  mseq_node(int i, char c) 
  : aligned_base(i,c), weight(1.f) 
  {}
  
  float getWeight() const { 
    return weight; 
  }


  bool operator<(const mseq_node &rhs) const { 
    return aligned_base::operator<(rhs); 
  }
private:
    float weight;
    friend class mseq;
};

class mseq : public dag<mseq_node>
{
public:
  mseq() 
  : num_seqs(0), bases_width(0) 
  {}

  mseq(std::vector<cseq>::iterator,
       std::vector<cseq>::iterator,
       float);

  unsigned int getWidth() { 
    return bases_width; 
  }

private:
  unsigned int num_seqs;
  unsigned int bases_width;
};

} // namespace sina

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

