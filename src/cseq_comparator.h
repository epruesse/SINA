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

#ifndef _CSEQ_COMPARATOR_H_
#define _CSEQ_COMPARATOR_H_

#include "cseq.h"
#include <boost/program_options.hpp>

namespace sina {

/**
 * Type defining how IUPAC encoded bases should
 * be compared. 
 * OPTIMISTIC: A=N
 * PESSIMISTIC: A!=N
 */
enum CMP_IUPAC_TYPE {
    CMP_IUPAC_OPTIMISTIC,
    CMP_IUPAC_PESSIMISTIC
};

/**
 * Type defining which distance correction should
 * be applied. 
 * NONE: do not correct
 * JC: use Jukes Cantor
 */
enum CMP_DIST_TYPE {
    CMP_DIST_NONE,
    CMP_DIST_JC
};

/**
 * Type defining relative to what the distance/identity/similarity
 * should be computet.
 * ABS: absolute (relative to 1)
 * QUERY: relative to query length
 * TARGET: relative to reference sequence length
 * OVERLAP: relative to length of overlapping part
 * ALL: relative to total alignment positions
 *      e.g.:  s1:------------AGCUAGCU-----
 *             s2:----------------AGCUAGCU-
 *             positions:     ^^^^^^^^^^^^ = 12
 *             equal:             ^^^^     = 4
 *             distance:          4/12     = 0.33333...
 * AVERAGE: relative to average of query and target length
 * MIN, MAX: relative to length of shorter/longer of query/target
 * NOGAP: relative to match/mismatch columns
 */
enum CMP_COVER_TYPE {
    CMP_COVER_ABS,
    CMP_COVER_QUERY,
    CMP_COVER_TARGET,
    CMP_COVER_OVERLAP,
    CMP_COVER_ALL,
    CMP_COVER_AVERAGE,
    CMP_COVER_MIN,
    CMP_COVER_MAX,
    CMP_COVER_NOGAP
};

/**
 * validator interpreting CLI encoding of CMP_IUPAC_TYPE
 * for boost_program_options
 */
void 
validate(boost::any&, const std::vector<std::string>&,
         CMP_IUPAC_TYPE*, int);

/**
 * validator interpreting CLI encoding of CMP_DIST_TYPE
 * for boost_program_options
 */
void 
validate(boost::any&, const std::vector<std::string>&,
         CMP_DIST_TYPE*, int);

/**
 * validator interpreting CLI encoding of CMP_COVER_TYPE
 * for boost_program_options
 */
void 
validate(boost::any&, const std::vector<std::string>&,
         CMP_COVER_TYPE*, int);


class cseq_comparator {
public:
    /**
     * Builds options_description to pass class parameters
     * to CLI interface.
     * @param prefix   prefix for option names 
     */
    static boost::program_options::options_description
    get_options_description(const char* prefix="");
    
    /**
     * Factory method building a comparator from the variables
     * map retrieved via boost program_options.
     * @param vm       variables map containing parsed cmdline
     * @param prefix   prefix for option names 
     */
    static cseq_comparator
    make_from_variables_map(boost::program_options::variables_map& vm,
                            const char* prefix="");
    
    /** 
     * Constructor taking explicit configuration.
     * See type definition for an explanation of the values.
     * @param filter_lc if true lower case bases at the beginning
     *                  or end are ignored.
     */
    cseq_comparator(CMP_IUPAC_TYPE iupac, CMP_DIST_TYPE dist, 
                    CMP_COVER_TYPE cover, bool filter_lc);
    
    /**
     * Default constructor
     */
    cseq_comparator();
    
    /**
     * Operator provided by the comparator object.
     * @param query  query sequence 
     * @param target target sequence 
     */
    float 
    operator()(const cseq& query, const cseq& target);
    
private:
    CMP_IUPAC_TYPE iupac_rule;
    CMP_DIST_TYPE dist_rule;
    CMP_COVER_TYPE cover_rule;
    bool filter_lc_rule;
};

/**
 * ostream output operator for CMP_IUPAC_TYPE
 */
std::ostream& 
operator<<(std::ostream&, const CMP_IUPAC_TYPE&);

/**
 * ostream output operator for CMP_IUPAC_TYPE
 */
std::ostream& 
operator<<(std::ostream&, const CMP_DIST_TYPE&);

/**
 * ostream output operator for CMP_IUPAC_TYPE
 */
std::ostream& 
operator<<(std::ostream&, const CMP_COVER_TYPE&);

} // namespace sina


#endif // _CSEQ_COMPARATOR_H_

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
