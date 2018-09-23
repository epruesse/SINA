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

#include "aligned_base.h"
#include <sstream>

#ifndef TEST

using namespace sina;

std::ostream&
operator<<(std::ostream& out, aligned_base ab)
{
    std::stringstream tmp;
    tmp <<  ab.getBase() << "(" << ab.getPosition() << ")";
    return out << tmp.str();
}


/* === IUPAC ==
 *
 * A = adenine
 * C = cytosine
 * G = guanine
 * T= thymine
 * U = uracil
 * R = G A (purine)
 * Y = T C (pyrimidine)
 * K = G T (keto)
 * M = A C (amino)
 * S = G C
 * W = A T
 * B = G T C
 * D = G A T
 * H = A C T
 * V = G C A
 * N = A G C T (any)
 *
 * Complements:
 * http://www.animalgenome.org/edu/gene/genetic-code.html
 *
 *
 */

const base_iupac::value_type
base_iupac::iupac_char_to_bmask[] = {
    0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  0, 0, 0, 0, // 0-7
    0, 0, 0, 0,  0, 0, 0, 0, // 8-?

    0, BASEM_A, BASEM_G|BASEM_TU|BASEM_C, BASEM_C,  BASEM_G|BASEM_A|BASEM_TU, 0, 0, BASEM_G, // @-G
    BASEM_A|BASEM_C|BASEM_TU, 0, 0, BASEM_G|BASEM_TU,  0, BASEM_A|BASEM_C, BASEM_A|BASEM_G|BASEM_C|BASEM_TU, 0, // H-O
    0, 0, BASEM_G|BASEM_A, BASEM_G|BASEM_C,  BASEM_TU, BASEM_TU, BASEM_G|BASEM_C|BASEM_A, BASEM_A|BASEM_TU, // P-W
    0, BASEM_TU|BASEM_C, 0, 0,  0, 0, 0, 0, // X-_
    0, BASEM_LC|BASEM_A, BASEM_LC|BASEM_G|BASEM_TU|BASEM_C, BASEM_LC|BASEM_C,  BASEM_LC|BASEM_G|BASEM_A|BASEM_TU, 0, 0, BASEM_LC|BASEM_G, // `-g
    BASEM_LC|BASEM_A|BASEM_C|BASEM_TU, 0, 0, BASEM_LC|BASEM_G|BASEM_TU,  0, BASEM_LC|BASEM_A|BASEM_C, BASEM_LC|BASEM_A|BASEM_G|BASEM_C|BASEM_TU, 0, // h-o
    0, 0, BASEM_LC|BASEM_G|BASEM_A, BASEM_LC|BASEM_G|BASEM_C,  BASEM_LC|BASEM_TU, BASEM_LC|BASEM_TU, BASEM_LC|BASEM_G|BASEM_C|BASEM_A, BASEM_LC|BASEM_A|BASEM_TU, // p-w
    0, BASEM_LC|BASEM_TU|BASEM_C, 0, 0,  0, 0, 0, 0, // x-

    0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  0, 0, 0, 0,

    0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  0, 0, 0, 0
};

const unsigned char
base_iupac::bmask_to_iupac_rna_char[] = {
 //0   1   10   11  100  101  110  111  1000 1001 1010 1011 1100 1101 1110 1111
 '.', 'A', 'G', 'R', 'C', 'M', 'S', 'V', 'U', 'W', 'K', 'D', 'Y', 'H', 'B', 'N',
 '.', 'a', 'g', 'r', 'c', 'm', 's', 'v', 'u', 'w', 'k', 'd', 'y', 'h', 'b', 'n'
};

const unsigned char
base_iupac::bmask_to_iupac_dna_char[] = {
 //0   1   10   11  100  101  110  111  1000 1001 1010 1011 1100 1101 1110 1111
 '.', 'A', 'G', 'R', 'C', 'M', 'S', 'V', 'T', 'W', 'K', 'D', 'Y', 'H', 'B', 'N',
 '.', 'a', 'g', 'r', 'c', 'm', 's', 'v', 't', 'w', 'k', 'd', 'y', 'h', 'b', 'n'
};

const float
base_iupac::base_pairings[] = {
//  .    A   G   R   C   M   S   V   U   W   K   D   Y   H   B   N
    0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f, // .
    0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,1.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f, // A
    0.f,0.f,0.f,0.f,1.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f, // G
    0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f, // R
    0.f,0.f,1.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f, // C
    0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f, // M
    0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f, // S
    0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f, // V
    0.f,1.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f, // U
    0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f, // W
    0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f, // K
    0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f, // D
    0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f, // Y
    0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f, // H
    0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f,0.f, // B
};

#else // TEST

int main(int, char**) {
    return 0;
}

#endif // TEST


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

