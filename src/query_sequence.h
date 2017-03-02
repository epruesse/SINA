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
 * data_interface.h
 * A basic interface to query sequences.
 *  Created on: Mar 12, 2012
 *      Author: aboeckma
 */

#ifndef QUERY_SEQUENCE_H_
#define QUERY_SEQUENCE_H_

#include <map>
#include <string>
#include <vector>


class query_exception : public std::exception
{
public:
    query_exception(std::string& what): m_what(what){}
    query_exception(const std::string& what): m_what(what){}
    query_exception(const char* what): m_what(what){}
    query_exception(){}
    virtual const char* what() const throw(){ return m_what.c_str();}
    virtual ~query_exception() throw(){}

private:
    std::string m_what;
};

/**
 * Interface for different sequence providers.
 */
class query_sequence{

    class cseq;

    /**
     * TODO error handling
     * TODO exception documentation
     */

    /**
     * TODO what is the alignment width? The number of bases?
     * The length of the string? no idea... yet :)
     */
    virtual size_t getAlignmentWidth() throw() = 0;


    /**
     * TODO
     */
    virtual std::vector<alignment_stats> getAlignmentStats() throw() = 0;

    /**
     *  Returns the sequence identified by the specified name.
     *  @throw query_exception if no matching sequence was found
     */
    virtual cseq& getSequence(const std::string& name) throw(query_exception) = 0;

    /**
     * Stores the specified sequence.
     */
    virtual void putSequence(const cseq& sequence) throw(query_exception) = 0;


    virtual size_t getSequenceCount() throw(query_exception) = 0;

    /**
     * Returns a list of sequence names.
     * The name format may differ depending on the underlying database.
     */
    virtual std::vector<std::string> getSequenceNames() throw(query_exception) = 0;

    /**
     * save using default name
     */
    virtual void save() throw(query_exception) = 0;

    /**
     * Save using the specified filename
     */
    virtual void saveAs(const std::string& fileName) throw(query_exception) = 0;

    /**
     * TODO no idea what this method does.
     * TODO Determine whether this method should be part of the interface or not.
     */
    virtual std::string getFilter(const std::string& name) = 0;

};



#endif /* QUERY_SEQUENCE_H_ */
