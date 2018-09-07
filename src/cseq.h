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

#ifndef _CSEQ_H_
#define _CSEQ_H_

#include <string>
#include <list>
#include <vector>
#include <map>
#include <sstream>
#include <iostream>
#include <algorithm>
#include "aligned_base.h"

#include <boost/variant.hpp>
#include <boost/lexical_cast.hpp>
/*
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/level.hpp>
#include <boost/serialization/tracking.hpp>
#include <boost/serialization/access.hpp>
*/
#include <boost/tuple/tuple.hpp>

namespace sina {

class cseq;
template<typename T> class lexical_cast_visitor;

/* compressed sequence. instead of gapping the sequence, all characters
   have positions */
class cseq {
public:
    typedef unsigned int idx_type;
    typedef aligned_base::idx_type vidx_type;
    typedef aligned_base value_type;
    class iterator;
    class const_iterator;
    class const_reverse_iterator;
    typedef iterator pn_iterator;
    typedef const_iterator const_pn_iterator;
    typedef boost::variant<std::string, char, int, float, std::vector<cseq> > variant;

    // Constructors / assignment operator

    cseq(const char* name, float score = 0.f, const char* data = NULL);
    cseq();
    cseq& operator=(const cseq& rhs); 
    cseq(const cseq& orig);

    // writing methods
    void clearSequence();

    cseq& append(const char* aligned_rna_sequence);
    cseq& append(const std::string aligned_rna_sequence);
    cseq& append(const aligned_base& a);

    //FIXME what does assign do?
    cseq& assign(std::vector<unsigned char> &dat);
  void assignFromCompressed(const void *buffer, size_t len);

    void setWidth(vidx_type pos);
    //FIXME what does fix_duplicate_positions do?
    void fix_duplicate_positions(std::ostream&, bool lowercase, bool remove);
    //FIXME cseq does not inherit anything. why does the empty sort method exist?
    void sort() {} // does nothing, cseq is always sorted by positions
    void reverse();
    void complement();
    void upperCaseAll();


    // eval methods

    std::vector<aligned_base> getAlignedBases() { return bases; }
    const std::vector<aligned_base>& const_getAlignedBases() const { return bases; }
    void setAlignedBases(const std::vector<aligned_base>& vab) { bases = vab; }
    std::string getAligned(bool nodots=false, bool dna=false) const;
    std::string getAlignedNoDots() const {return getAligned(true);}
    std::string getBases() const;
    std::string getName() const { return name; }
    void setName(std::string n) { name=n; }
    std::string getNameScore() const;
    // fixme: handle "-" and "." correctly
    //FIXME where is the difference between size and width? why is there a difference?
    vidx_type size() const { return bases.size(); }
    vidx_type getWidth() const { return alignment_width; }

    void compressAligned(std::vector<unsigned char> &out);

    //Fixme how is the score calculated?
    //      what should the content of pairs be?
    float calcPairScore(const std::vector<int>& pairs);

    iterator begin();
    const_iterator begin() const;
    const_reverse_iterator rbegin() const;
    iterator end();
    const_iterator end() const;
    const_reverse_iterator rend() const;

    iterator getIterator(cseq::vidx_type i);
    const_iterator getIterator(cseq::vidx_type i) const;

    char operator [](vidx_type);
    const aligned_base& getById(idx_type i) const { return bases[i]; } 

    // meta-data

    float getScore() const { return score; }
    void setScore(float f) { score = f; }

    // io operations dealing with vector<*cseq>
    static void write_alignment(std::ostream& ofs, std::vector<cseq>& seqs,
                                cseq::idx_type start, cseq::idx_type stop, 
                                bool color_code = false);    
    static void write_alignment(std::ostream& ofs, std::vector<cseq>& seqs,
                                bool color_code = false);

    template<typename T>
    void set_attr(std::string key, T val) {
        attributes[key]=val;
    }

    template<typename T>
    T get_attr(std::string attr) {
        return boost::apply_visitor(lexical_cast_visitor<T>(), attributes[attr]); 
    }

    const std::map<std::string,variant>& get_attrs() const { return attributes; }


    friend std::ostream& operator<<(std::ostream& out, const cseq&);

    bool operator==(const cseq& rhs) const {
        return name == rhs.name && bases == rhs.bases &&
            attributes == rhs.attributes;// && score-rhs.score<0.0001;
    }
    bool operator!=(const cseq& rhs) const {
        return name != rhs.name || bases != rhs.bases ||
            attributes != rhs.attributes; // || score != rhs.score;
    }
    bool operator<(const cseq& rhs) const { return score < rhs.score; }
    bool operator>(const cseq& rhs) const { return score > rhs.score; }
    std::list<unsigned int> find_differing_parts(const cseq&) const;
protected:

private:
    std::string name;
    std::vector<aligned_base> bases;
    unsigned int alignment_width;
    std::map<std::string,variant> attributes;
    float score;
    
    template <typename FUNC>
    friend void traverse(const cseq& A, const cseq& B, FUNC F);
};

template<>
inline cseq::variant cseq::get_attr<cseq::variant>(std::string attr) { return attributes[attr]; }

std::ostream& operator<<(std::ostream& out, const cseq&);



class cseq::iterator
    : public std::vector<aligned_base>::iterator
{
public:
    iterator(std::vector<aligned_base>::iterator it)
        : std::vector<aligned_base>::iterator(it) {}
    iterator() : std::vector<aligned_base>::iterator() {}

    typedef iterator pn_iterator;
    iterator prev_begin() const { iterator n(*this); return --n; }
    iterator prev_end() const { return (*this); }
    iterator next_begin() const { iterator n(*this); return ++n; }
    iterator next_end() const { iterator n(*this); return n+2; }
protected:

private:

};

class cseq::const_iterator
    : public std::vector<aligned_base>::const_iterator
{
 public:
    const_iterator(std::vector<aligned_base>::const_iterator it)
        : std::vector<aligned_base>::const_iterator(it) {}
    const_iterator() : std::vector<aligned_base>::const_iterator() {}

    typedef const_iterator const_pn_iterator;
    const_iterator prev_begin() const {
        const_iterator n(*this); return --n; }
    const_iterator prev_end() const { return (*this); }
    const_iterator next_begin() const {
        const_iterator n(*this); return ++n; }
    const_iterator next_end() const {
        const_iterator n(*this); return n+2; }
 protected:

 private:

};

class cseq::const_reverse_iterator
    : public std::vector<aligned_base>::const_reverse_iterator
{
 public:
    const_reverse_iterator(std::vector<aligned_base>::const_reverse_iterator it)
      : std::vector<aligned_base>::const_reverse_iterator(it) {}
    const_reverse_iterator() : std::vector<aligned_base>::const_reverse_iterator() {}

  /*
    typedef const_iterator const_pn_iterator;
    const_iterator prev_begin() const {
        const_iterator n(*this); return --n; }
    const_iterator prev_end() const { return (*this); }
    const_iterator next_begin() const {
        const_iterator n(*this); return ++n; }
    const_iterator next_end() const {
        const_iterator n(*this); return n+2; }
  */
 protected:

 private:

};


} // namespace sina

#include "cseq_impl.h"

#if 0
// elminate serialization overhead at the cost of
// never being able to increase the version.
BOOST_CLASS_IMPLEMENTATION(sina::cseq,
                           boost::serialization::object_serializable);

// eliminate object tracking (even if serialized through a pointer)
// at the risk of a programming error creating duplicate objects.
BOOST_CLASS_TRACKING(sina::cseq,
                     boost::serialization::track_never);
#endif


#endif

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

