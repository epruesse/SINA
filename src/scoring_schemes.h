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

#ifndef _SCORING_SCHEMES_H_
#define _SCORING_SCHEMES_H_

#include <vector>

namespace sina {


class scoring_scheme_profile {
public:
    using value_type = float;

    scoring_scheme_profile(value_type _m, value_type _mm,
                           value_type _gp, value_type _gpe)
        : match_score(_m), mismatch_score(_mm),
          gap_penalty(_gp), gap_extension_penalty(_gpe)
    {}

    template<typename base_type_a, typename base_type_b>
    value_type
    insertion(value_type prev,
              const base_type_a& /*b1*/,
              const base_type_b& /*b2*/) const
    {
        return prev + gap_penalty;
    }

    template<typename base_type_a, typename base_type_b>
    value_type
    insertion_ext(value_type prev,
                  const base_type_a& /*b1*/,
                  const base_type_b& /*b2*/,
                  int /*offset*/) const
    {
        return prev + gap_extension_penalty;
    }

    template<typename base_type_a, typename base_type_b>
    value_type
    deletion(value_type prev,
             const base_type_a& /*b1*/,
             const base_type_b& /*b2*/) const
    {
        return prev + gap_penalty;
    }

    template<typename base_type_a, typename base_type_b>
    value_type
    deletion_ext(value_type prev,
                 const base_type_a& /*b1*/,
                 const base_type_b& /*b2*/,
                 int /*offset*/) const
    {
        return prev + gap_extension_penalty;
    }

    template<typename base_type_a, typename base_type_b>
    value_type
    match(value_type prev,
          const base_type_a& b1,
          const base_type_b& b2) const
    {
        return prev + (b1.comp(b2, match_score, mismatch_score,
                               gap_penalty, gap_extension_penalty));
    }

private:
    const value_type match_score;
    const value_type mismatch_score;
    const value_type gap_penalty;
    const value_type gap_extension_penalty;
};

class scoring_scheme_simple {
public:
    using value_type = float;

    scoring_scheme_simple(value_type _m, value_type _mm,
                          value_type _gp, value_type _gpe)
        : match_score(_m), mismatch_score(_mm),
          gap_penalty(_gp), gap_extension_penalty(_gpe)
    {}

    template<typename base_type_a, typename base_type_b>
    value_type
    insertion(value_type prev,
              const base_type_a& /*b1*/,
              const base_type_b& /*b2*/) const
    {
        return prev + gap_penalty;
    }

    template<typename base_type_a, typename base_type_b>
    value_type
    insertion_ext(value_type prev,
                  const base_type_a& /*b1*/,
                  const base_type_b& /*b2*/,
                  int /*offset*/) const
    {
        return prev + gap_extension_penalty;
    }

    template<typename base_type_a, typename base_type_b>
    value_type
    deletion(value_type prev,
             const base_type_a& /*b1*/,
             const base_type_b& /*b2*/) const
    {
        return prev + gap_penalty;// * b1.getWeight();
    }

    template<typename base_type_a, typename base_type_b>
    value_type
    deletion_ext(value_type prev,
                 const base_type_a& /*b1*/,
                 const base_type_b& /*b2*/,
                 int /*offset*/) const
    {
        return prev + gap_extension_penalty;// * b1.getWeight();
    }

    template<typename base_type_a,typename base_type_b>
    value_type
    match(value_type prev,
          const base_type_a& b1,
          const base_type_b& b2) const
    {
        return prev + (b1.comp(b2)?match_score:mismatch_score) * b1.getWeight();
    }

private:
    const value_type match_score;
    const value_type mismatch_score;
    const value_type gap_penalty;
    const value_type gap_extension_penalty;
};

class scoring_scheme_weighted {
public:
    using value_type = float;

    scoring_scheme_weighted(value_type _m, value_type _mm,
                            value_type _gp, value_type _gpe,
                            std::vector<value_type>& _weights)
        : match_score(_m), mismatch_score(_mm),
          gap_penalty(_gp), gap_extension_penalty(_gpe),
          weights(_weights)
    {}

    template<typename base_type_a, typename base_type_b>
    value_type
    insertion(value_type prev,
              const base_type_a& b1,
              const base_type_b& /*b2*/) const
    {
        // the insertion will be placed in the column
        // immediately following the current column. use
        // the weight from that column.
        return prev + gap_penalty * weights[b1.getPosition() +1];
    }

    template<typename base_type_a, typename base_type_b>
    value_type
    insertion_ext(value_type prev,
                  const base_type_a& b1,
                  const base_type_b& /*b2*/,
                  int offset) const
    {
        // the insertion will be placed in the column
        // immediately following the current column plus
        // offset (number of previous insertions)
        return prev + gap_extension_penalty
            * weights[b1.getPosition() + 1 + offset];
    }

    template<typename base_type_a, typename base_type_b>
    value_type
    deletion(value_type prev,
             const base_type_a& b1,
             const base_type_b& /*b2*/) const
    {
        // deletions always happen in current column
        return prev + gap_penalty * weights[b1.getPosition()];// * b1.getWeight();
    }

    template<typename base_type_a, typename base_type_b>
    value_type
    deletion_ext(value_type prev,
                 const base_type_a& b1,
                 const base_type_b& /*b2*/,
                 int /*offset*/) const
    {
        // deletions always happen in current column
        return prev + gap_extension_penalty
            * weights[b1.getPosition()] ;//* b1.getWeight();
    }

    template<typename base_type_a, typename base_type_b>
    value_type
    match(value_type prev,
          const base_type_a& b1,
          const base_type_b& b2) const
    {
        return prev + (b1.comp(b2)?match_score:mismatch_score) * weights[b1.getPosition()] * b1.getWeight();
    }

private:
    const value_type match_score;
    const value_type mismatch_score;
    const value_type gap_penalty;
    const value_type gap_extension_penalty;
    std::vector<value_type>& weights;
};


template<typename MATRIX_TYPE>
class scoring_scheme_matrix {
public:
    using value_type = float;
    using matrix_type = MATRIX_TYPE;

    scoring_scheme_matrix(value_type _gp, value_type _gpe,
                          std::vector<value_type>& _weights,
                          const matrix_type& _matrix)
        : gap_penalty(_gp),
          gap_extension_penalty(_gpe),
          weights(_weights),
          matrix(_matrix)
    {}

    template<typename base_type_a, typename base_type_b>
    value_type
    insertion(value_type prev,
              const base_type_a& b1,
              const base_type_b& /*b2*/) const
    {
        return prev + gap_penalty * weights[b1.getPosition()];
    }

    template<typename base_type_a, typename base_type_b>
    value_type
    insertion_ext(value_type prev,
                  const base_type_a& b1,
                  const base_type_b& /*b2*/,
                  int /*offset*/) const
    {
        return prev + gap_extension_penalty * weights[b1.getPosition()];
    }

    template<typename base_type_a, typename base_type_b>
    value_type
    deletion(value_type prev,
             const base_type_a& b1,
             const base_type_b& b2) const
    {
        return insertion(prev, b1, b2);
    }

    template<typename base_type_a, typename base_type_b>
    value_type
    deletion_ext(value_type prev,
                 const base_type_a& b1,
                 const base_type_b& b2,
                 int offset) const
    {
        return insertion_ext(prev, b1, b2, offset);
    }

    template<typename base_type_a, typename base_type_b>
    value_type
    match(value_type prev,
          const base_type_a& b1,
          const base_type_b& b2) const
    {
        return prev + b1.comp(b2,matrix) * weights[b1.getPosition()];
    }

private:
    const value_type gap_penalty;
    const value_type gap_extension_penalty;
    const std::vector<value_type>& weights;
    const matrix_type& matrix;
};


} // namespace sina

#endif // _MESH_H_

/*
  Local Variables:
  mode:c++
  c-file-style:"stroustrup"
  c-file-offsets:((innamespace . 0)(inline-open . 0)(case-label . 0))
  indent-tabs-mode:nil
  fill-column:99
  End:
*/
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:encoding=utf-8:textwidth=99 :
