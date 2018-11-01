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

#ifndef _ALIGNMENT_STATS_H_
#define _ALIGNMENT_STATS_H_

#include "aligned_base.h"

#include <vector>
#include <string>

namespace sina {
  
class alignment_stats {
public:
  struct freqs {
    unsigned int num_a{0};
    unsigned int num_g{0};
    unsigned int num_c{0};
    unsigned int num_u{0};
    unsigned int num_mutations{0};
    unsigned int num_transversions{0};
  };
  alignment_stats();
  alignment_stats(std::string  name,
                  unsigned int ntaxa, unsigned int alen,
		  unsigned int *na, unsigned int *nc, unsigned int *ng, 
		  unsigned int *nu, unsigned int *nM, unsigned int *nT,
                  std::vector<int>  _pairs);
  
  const std::vector<float>& getWeights() const { return weights; }
  const std::vector<int>& getPairs() const { return pairs; }
  const aligned_base::matrix_type getSubstMatrix(double identity) const;
  int getWidth() const {return width;} 
  std::string getName() const {return name;}
private:
  std::string name;
  unsigned int num_taxa{0};
  unsigned int width{0};
  
  freqs global_freqs;
  std::vector<freqs> column_freqs;

  std::vector<int> pairs;

  std::vector<float> weights;
  float maxweight{0}, minweight{0}, sumweight{0};
  int weighted_columns{0};
};
  
}
#endif //_ALIGNMENT_STATS_H_

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
