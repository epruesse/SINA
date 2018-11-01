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
#include "log.h"

#include <cmath>
#include <algorithm>
#include <iostream>
#include <utility>
using std::string;
using std::vector;
namespace sina {

float jukes_cantor(float in) {
    return -3.0/4 * log( 1.0 - 4.0/3*in);
}

alignment_stats::alignment_stats() {

    global_freqs.num_a=1000;
    global_freqs.num_g=1000;
    global_freqs.num_c=1000;
    global_freqs.num_u=1000;
    global_freqs.num_mutations=20;
    global_freqs.num_transversions=10;
}

alignment_stats::alignment_stats(
    std::string  name_,
    unsigned int ntaxa, unsigned int alen,
    const unsigned int *na, const unsigned int *nc, const unsigned int *ng,
    const unsigned int *nu, const unsigned int *nM, const unsigned int *nT,
    std::vector<int>  pairs_
    ) 
    : name(std::move(name_)),
      num_taxa(ntaxa), width(alen),
      pairs(std::move(pairs_)),
      minweight(9999999)
{
    auto console = Log::create_logger("alignment_stats");

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
                console->info("extreme weight {} for column {} clamped to 20", weight, i);
                weight = 20;
            }
            weights[i] = weight;
            sumweight += weight;
            maxweight = std::max(maxweight, (float)weight);
            minweight = std::min(minweight, (float)weight);
            weighted_columns++;
            if (i < first_weighted) {
                first_weighted = i;
            }
            if (i > last_weighted) {
                last_weighted = i;
            }
        } else {
            weights[i] = 1;
        }
    }

    int total_bases = global_freqs.num_a + global_freqs.num_c + 
        global_freqs.num_g + global_freqs.num_u;

    console->info("alignment stats for subset {}", name);
    console->info("weighted/unweighted columns = {}/{}",
                  weighted_columns, width - weighted_columns);
    console->info("average weight = {}", sumweight / weighted_columns);
    console->info("minimum weight = {}", minweight);
    console->info("maximum weight = {}", maxweight);
    console->info("ntaxa = {}", ntaxa);
    console->info("base frequencies: na={} nu={} nc={} ng={}",
                  (double)global_freqs.num_a/total_bases,
                  (double)global_freqs.num_c/total_bases,
                  (double)global_freqs.num_g/total_bases,
                  (double)global_freqs.num_u/total_bases);
    console->info("mutation frequencies: any={} transversions={}",
                  (double)global_freqs.num_mutations/total_bases,
                  (double)global_freqs.num_transversions/total_bases);
    console->info("first/last weighted column={}/{}",
                  first_weighted, last_weighted);
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

#if 0
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
# endif

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
