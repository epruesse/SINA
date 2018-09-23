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

/*
 * queryarbadapter.h
 *
 *  Created on: Mar 20, 2012
 *      Author: aboeckma
 */

#ifndef QUERYARBADAPTER_H_
#define QUERYARBADAPTER_H_

#include "query_sequence.h"

namespace sina
{

/**
 * Adapts query_arb to the query_sequence interface
 */
class query_arb_adapter : public query_sequence
{

    class cseq;
    /**
     * TODO
     *  * reference documentation from interface
     */


public:

    query_arb_adapter(const string& filename);

    virtual ~query_arb_adapter();

    virtual size_t getAlignmentWidth() throw();

    virtual cseq& getSequence(const std::string& name) throw(query_exception);

    virtual void putSequence(const cseq& sequence) throw(query_exception)

    virtual size_t getSequenceCount() throw(query_exception);

    virtual std::vector<std::string> getSequenceNames() throw(query_exception);

    virtual void save() throw(query_exception);

    virtual void saveAs(const std::string& fileName) throw(query_exception);

    virtual std::string getFilter(const std::string& name);

private:

    query_arb* m_pQuery_arb;


};

} /* namespace sina */
#endif /* QUERYARBADAPTER_H_ */
