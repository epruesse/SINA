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
 * queryarbadapter.cpp
 *
 *  Created on: Mar 20, 2012
 *      Author: aboeckma
 */

#include "query_arb_adapter.h"
using namespace std;

namespace sina
{

query_arb_adapter::query_arb_adapter(const string& filename) :
        m_pQuery_arb(query_arb::getARBDB(filename))
{
    //caching should be transparent to the user.
    // make set of fields to load
    std::vector<std::string> fields;
    fields.push_back("acc");
    fields.push_back("start");
    fields.push_back("stop");
    // cache sequences and meta data
    m_pQuery_arb->loadCache(fields);

}

query_arb_adapter::~query_arb_adapter()
{
    delete m_pQuery_arb;
}



size_t query_arb_adapter::getAlignmentWidth()
{
    return m_pQuery_arb->getAlignmentWidth();
}



cseq & query_arb_adapter::getSequence(const std::string & name)
{
    return m_pQuery_arb->getCseq(name);
}



void query_arb_adapter::putSequence(const cseq & sequence)
{
    m_pQuery_arb->putCseq(sequence);
}




size_t query_arb_adapter::getSequenceCount()
{
    return m_pQuery_arb->getSeqCount();
}



std::vector<std::string> query_arb_adapter::getSequenceNames()
{
    return m_pQuery_arb->getSequenceNames();
}



void query_arb_adapter::save()
{
    m_pQuery_arb->save();

}



void query_arb_adapter::saveAs(const std::string & fileName)
{
    //save as binary file
    m_pQuery_arb->saveAs(fileName.c_str());

}

std::string query_arb_adapter::getFilter(const std::string& name)
{
    return m_pQuery_arb->getFilter(name);

}

 /* namespace sina */
