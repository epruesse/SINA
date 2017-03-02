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
