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

#include "alignment_stats.h"
#include <cmath>
#include <algorithm>
#include <iostream>
using std::cerr;
using std::endl;
using std::string;
using std::vector;
namespace sina {

float jukes_cantor(float in) {
    return -3.0/4 * log( 1.0 - 4.0/3*in);
}

alignment_stats::alignment_stats() 
    : num_taxa(0), width(0), global_freqs(), maxweight(0), 
      minweight(0), sumweight(0), weighted_columns(0)
{

    global_freqs.num_a=1000;
    global_freqs.num_g=1000;
    global_freqs.num_c=1000;
    global_freqs.num_u=1000;
    global_freqs.num_mutations=20;
    global_freqs.num_transversions=10;
}

alignment_stats::alignment_stats(
    const std::string& _name,
    unsigned int ntaxa, unsigned int alen,
    unsigned int *na, unsigned int *ng, unsigned int *nc, 
    unsigned int *nu, unsigned int *nM, unsigned int *nT,
    const std::vector<int>& _pairs
    ) 
    : name(_name),
      num_taxa(ntaxa), width(alen), global_freqs(), 
      pairs(_pairs),
      maxweight(0), 
      minweight(9999999), sumweight(0), weighted_columns(0)
{
    column_freqs.resize(width);
    weights.resize(width);
    int first_weighted = width, last_weighted=0;
    for (int i=0; i<width; i++) {
        freqs &f = column_freqs[i];
        f.num_a = na[i];
        f.num_c = nc[i];
        f.num_g = ng[i];
        f.num_u = nu[i];
        f.num_mutations = nM[i];
        f.num_transversions = nT[i];

        global_freqs.num_a += f.num_a;
        global_freqs.num_c += f.num_c;
        global_freqs.num_g += f.num_g;
        global_freqs.num_u += f.num_u;
        global_freqs.num_mutations += f.num_mutations;
        global_freqs.num_transversions += f.num_transversions;

        int sum = f.num_a + f.num_c + f.num_g + f.num_u;
        if (sum > ntaxa * 0.2) {
            double rate = std::min((double)f.num_mutations / sum, .95 * .75);
            rate = std::min(jukes_cantor(rate), 1.f);
            double weight = .5 - log(rate);
            if (weight > 20 ) {
                std:: cerr << "extreme weight '" << weight << "'for column " << i 
                           << " clamped to 20" << std::endl;
                weight = 20;
            }
            weights[i] = weight;
            sumweight += weight;
            maxweight = std::max(maxweight, (float)weight);
            minweight = std::min(minweight, (float)weight);
            weighted_columns++;
            if (i<first_weighted) first_weighted=i;
            if (i>last_weighted) last_weighted=i;
        } else {
            weights[i] = 1;
        }
    }

    int total_bases = global_freqs.num_a + global_freqs.num_c + 
        global_freqs.num_g + global_freqs.num_u;
  
    cerr << "alignment stats for subset "<< name << endl 
         << "weighted/unweighted columns = " << weighted_columns 
         << "/" << width - weighted_columns << endl
         << "average weight = " << sumweight / weighted_columns << endl
         << "minimum weight = " << minweight << endl
         << "maximum weight = " << maxweight << endl
         << "ntaxa = " << ntaxa << endl
         << "base frequencies: na=" << (double)global_freqs.num_a/total_bases
         <<" nc=" << (double)global_freqs.num_c/total_bases
         <<" ng=" << (double)global_freqs.num_g/total_bases
         <<" nu=" << (double)global_freqs.num_u/total_bases << endl
         << "mutation frequencies: any=" << (double)global_freqs.num_mutations/total_bases
         <<" transversions=" << (double)global_freqs.num_transversions/total_bases << endl
         <<" first/last weighted column=" << first_weighted <<"/"<< last_weighted << endl
        ;
}

const aligned_base::matrix_type 
alignment_stats::getSubstMatrix(double identity) const {
    aligned_base::matrix_type m;
    int total_bases = global_freqs.num_a + global_freqs.num_c + 
        global_freqs.num_g + global_freqs.num_u;
    double f[BASE_MAX];
    f[BASE_A] = (double)global_freqs.num_a/total_bases;
    f[BASE_C] = (double)global_freqs.num_c/total_bases;
    f[BASE_G] = (double)global_freqs.num_g/total_bases;
    f[BASE_TU] = (double)global_freqs.num_u/total_bases;

    double avgmm=0;
    for (int i=0; i < BASE_MAX; i++) {
        for (int j=0; j < BASE_MAX; j++) {
            double p;
            if (i==j) {
                p = identity / 4;
            } else {
                p = (1-identity) / 12;
            }
            avgmm += m.v[i*BASE_MAX+j] = -log( p / (f[i] * f[j]) );
        }
    }
    avgmm /= 12;

    std::cerr << "avgmm:"<<avgmm<<std::endl;
    return m;
}


/*
const aligned_base::matrix_type
alignment_stats::getHky85Matrix(double identity) {
    aligned_base::matrix_type m;

    int total_bases = global_freqs.num_a + global_freqs.num_c + 
        global_freqs.num_g + global_freqs.num_u;
    double f[aligned_base::BASE_MAX], fa, fc, fg, fu;
    fa = f[aligned_base::BASE_A] = (double)global_freqs.num_a/total_bases;
    fc = f[aligned_base::BASE_C] = (double)global_freqs.num_c/total_bases;
    fg = f[aligned_base::BASE_G] = (double)global_freqs.num_g/total_bases;
    fu = f[aligned_base::BASE_TU] = (double)global_freqs.num_u/total_bases;

    double k = global_freqs.num_transversions / 
        (global_freqs.num_mutations-global_freqs.num_transversions);
    double beta = 1 / 
        (2 * (fa+fg)*(fc+fu) + 2*k*(fa*fg + fc*ft));

}

*/

} // namespace sina

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
