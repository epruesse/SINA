/*
Copyright (c) 2006-2017 Elmar Pruesse <elmar.pruesse@ucdenver.edu>

This file is part of SINA.
123456789012345678901234567890123456789012345678901234567890123456789012
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

#ifndef _QUERY_ARB_H_
#define _QUERY_ARB_H_

#include <string>
#include <vector>
#include <list>
#include "cseq.h"
#include "alignment_stats.h"

struct GBDATA;

namespace sina {

/**
 * Base exception thrown by query_arb
 */
class query_arb_exception : public std::exception
{
public:
    query_arb_exception(std::string& what): m_what(what){}
    query_arb_exception(const std::string& what): m_what(what){}
    query_arb_exception(const char* what): m_what(what){}
    query_arb_exception(){}
    virtual const char* what() const throw(){ return m_what.c_str();}
    virtual ~query_arb_exception() throw(){}

private:
    std::string m_what;
};


class query_arb{

    query_arb(std::string);
    ~query_arb();



 public:
    /**
     * Factory method.
     * @param file_name ...
     * @returns one instance per database file. If the method is called multiple times on the same file the same instance is returned.
     */
    static query_arb* getARBDB(std::string file_name);

    static const char* fn_turn;
    static const char* fn_acc;
    static const char* fn_start;
    static const char* fn_qual;
    static const char* fn_head;
    static const char* fn_tail;
    static const char* fn_date;
    static const char* fn_astart;
    static const char* fn_astop;
    static const char* fn_idty;
    static const char* fn_family;
    static const char* fn_family_str;
    static const char* fn_nuc;
    static const char* fn_nuc_gene;
    static const char* fn_bpscore;
    static const char* fn_used_rels;
    static const char* fn_fullname;


    /**
     * Saves the database using the default filename in binary format.
     * @note The filename is specified in init().
     *       Do not call this method before init() has been called.
     *
     * @see saveAs(const char* fname, const char* type="b")
     */
    void save();

    /**
     * Saves the database using the specified filename.
     * @param fname the filename.
     * @param type Defines how the database should be saved.
     *        The value is a combination of one or more of the following
     *        values:
     *        * 'a' ascii
     *        * 'b' binary
     *        * 'm' save mapfile (only together with binary)
     *        * 'f' force saving even in disabled path to a different directory (out of order save)
     *        * 'S' save to stdout (for debugging)
     *        @note Either 'a' or 'b' must be defined.
     */
    void saveAs(const char* fname, const char* type="b");

    /**
     * Returns filename of underlying ARB database file
     */
    std::string getFileName() const;

    void setProtectionLevel(int);

    int getSeqCount() const;
    std::vector<std::string> getSequenceNames();


    cseq& getCseq(std::string name);
    void putCseq(const cseq&);
    void putSequence(const cseq&);//calls write
    void loadKey(cseq&, std::string);
    void storeKey(cseq&, std::string);




    long getAlignmentWidth();

    std::string getFilter(std::string name);
    std::vector<int> getPairs();
    std::vector<alignment_stats> getAlignmentStats();


    void loadCache();
    void loadCache(std::vector<std::string>&);
    std::vector<cseq*> getCacheContents();


private:

    /**
     * Prints an error statistic to the specified stream.
     */
    void printErrors(std::ostream& stream);

    bool bad() const { return !good(); }

    bool good() const;

    /**
     * @return True if errors occurred while accessing the arb-db. False if no errors occurred.
     */
    bool hasErrors() const;

    void setMark();
    void setMark(const std::string&);
    void setMark(const cseq& cs);


    void copySequence(query_arb& qa, const cseq& cs, bool m); //calls write
    void copySequence(query_arb& qa, const std::string s, bool m);

    // make query_arb non-copyable
    query_arb(const query_arb&);
    query_arb& operator=(const query_arb&);


    void init(const char*);
    static void closeOpenARBDBs();

    inline void
    storeKey(GBDATA* gbmain, GBDATA* gbspec, const std::string& key,
             cseq::variant var);

    /**
     * Wrapper for arbs GB_write_float method.
     * @throw query_arb_exception if pData is null.
     */
    void write(GBDATA* pData, double value);
    /**
     * Wrapper for arbs GB_write_float method.
     * @throw query_arb_exception if pData is null.
     */
    void write(GBDATA* pData, int value);
    /**
     * Wrapper for arbs GB_write_string method.
     * @throw query_arb_exception if pData is null.
     */
    void write(GBDATA* pData, const char* pValue);
    /**
     * Wrapper for arbs GB_write_flag method.
     * @throw query_arb_exception if pData is null.
     */
    void write_flag(GBDATA* pData, long value);

    /**
     * Adds an error to the internal error list
     */
    void addError(const std::string& message);

    struct priv_data;
    priv_data& data;
    struct storeKey_visitor;
    struct putKeyVal_visitor;
    static std::map<std::string, query_arb*> open_arb_dbs;

};

// inline implementations
inline void
query_arb::copySequence(query_arb& qa, const cseq& cs, bool m) {
    copySequence(qa,cs.getName(), m);
}

inline void
query_arb::setMark(const cseq& cs) {
    setMark(cs.getName());
}

inline void
query_arb::loadCache() {
  std::vector<std::string> vs; loadCache(vs);
}


} // namespace sina

#endif // _QUERY_ARB_H_

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
