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

#ifndef _CSEQ_H_
#define _CSEQ_H_

#include <string>
#include <list>
#include <utility>
#include <vector>
#include <map>
#include <sstream>
#include <iostream>
#include <algorithm>
#include "aligned_base.h"

#include <boost/variant.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tuple/tuple.hpp>

namespace sina {

/* compressed sequence. instead of gapping the sequence, all characters
   have positions */
class cseq_base {
public:
    using idx_type = unsigned int;
    using vidx_type = aligned_base::idx_type;
    using value_type = aligned_base;
    class iterator;
    class const_iterator;
    class const_reverse_iterator;
    using pn_iterator = iterator;
    using const_pn_iterator = const_iterator;

    // Constructors / assignment operator

    cseq_base(const char* _name, const char* _data = nullptr);
    cseq_base() = default;
    cseq_base& operator=(const cseq_base& rhs) = default;
    cseq_base(const cseq_base& orig) = default;

    /* remove sequence */
    void clearSequence();

    /* append bases to alignment */
    cseq_base& append(const char* str);
    cseq_base& append(const std::string& str);
    cseq_base& append(const aligned_base& ab);

    /* get size in bases */
    vidx_type size() const { return bases.size(); }

    /* get aligned base vector */
    const std::vector<aligned_base>& getAlignedBases() const { return bases; }

    /* get size in columns */
    vidx_type getWidth() const { return alignment_width; }

    /* set number of columns
     *  - shifts bases inwards from the right as needed
     *  - throws if newWidth smaller than size()
     */
    void setWidth(vidx_type newWidth);

    /* handle insertions (multiple bases with identical position)
     * created during alignment */
    void fix_duplicate_positions(std::ostream& /*log*/, bool lowercase, bool remove);

    /* does nothing, cseq_base is always sorted by positions */
    static void sort() {}

    /* reverse the sequence */
    void reverse();

    /* complement all bases */
    void complement();

    /* convert all bases to upper case */
    void upperCaseAll();

    // eval methods

    void setAlignedBases(const std::vector<aligned_base>& vab) { bases = vab; }
    std::string getAligned(bool nodots=false, bool dna=false) const;
    std::string getBases() const;

    /* get sequence name */
    std::string getName() const { return name; }

    /* set sequence name */
    void setName(std::string n) { name=std::move(n); }

    void assignFromCompressed(const void *data, size_t len);
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

    iterator getIterator(cseq_base::vidx_type i);
    const_iterator getIterator(cseq_base::vidx_type i) const;

    char operator [](vidx_type i) const;
    const aligned_base& getById(idx_type i) const { return bases[i]; } 


    // io operations dealing with vector<*cseq_base>
    static void write_alignment(std::ostream& ofs, std::vector<const cseq_base*>& seqs,
                                cseq_base::idx_type from_pos, cseq_base::idx_type to_pos,
                                bool colors = false);



    friend std::ostream& operator<<(std::ostream& out, const cseq_base& c);

    bool operator==(const cseq_base& rhs) const {
        return name == rhs.name && bases == rhs.bases; // &&
        //attributes == rhs.attributes;
           // && score-rhs.score<0.0001;
    }
    bool operator!=(const cseq_base& rhs) const {
        return name != rhs.name || bases != rhs.bases; // ||
        //attributes != rhs.attributes;
            // || score != rhs.score;
    }
    bool operator<(const cseq_base& rhs) const { return name < rhs.name; }
    bool operator>(const cseq_base& rhs) const { return name > rhs.name; }
    std::vector<std::pair<unsigned int, unsigned int>> find_differing_parts(const cseq_base& right) const;

private:
    std::string name;
    std::vector<aligned_base> bases;
    unsigned int alignment_width{0};

    template <typename FUNC>
    friend void traverse(const cseq_base& A, const cseq_base& B, FUNC F);
};


std::ostream& operator<<(std::ostream& out, const cseq_base& c);


class cseq_base::iterator
    : public std::vector<aligned_base>::iterator
{
public:
    iterator(std::vector<aligned_base>::iterator it)
        : std::vector<aligned_base>::iterator(it) {}
    iterator() = default;

    using pn_iterator = iterator;
    iterator prev_begin() const { iterator n(*this); return --n; }
    iterator prev_end() const { return (*this); }
    iterator next_begin() const { iterator n(*this); return ++n; }
    iterator next_end() const { iterator n(*this); return n+2; }
protected:

private:

};

class cseq_base::const_iterator
    : public std::vector<aligned_base>::const_iterator
{
 public:
    const_iterator(std::vector<aligned_base>::const_iterator it)
        : std::vector<aligned_base>::const_iterator(it) {}
    const_iterator() = default;

    using const_pn_iterator = const_iterator;
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

class cseq_base::const_reverse_iterator
    : public std::vector<aligned_base>::const_reverse_iterator
{
 public:
    const_reverse_iterator(const std::vector<aligned_base>::const_reverse_iterator& it)
      : std::vector<aligned_base>::const_reverse_iterator(it) {}
    const_reverse_iterator() = default;

 protected:

 private:

};

template<typename T> class lexical_cast_visitor;

class annotated_cseq : public cseq_base {
public:
    using variant = boost::variant<std::string, char, int, float>;

    annotated_cseq(const char* _name, const char* _data = nullptr)
        : cseq_base(_name, _data) {}
    annotated_cseq() : cseq_base() {}

    template<typename T>
    void set_attr(const std::string& key, T val) {
        attributes[key] = val;
    }

    bool has_attr(const std::string& key) const {
        return attributes.find(key) != attributes.end();
    }

    template<typename T>
    T get_attr(const std::string& attr) const {
        const auto it = attributes.find(attr);
        if (it != attributes.end()) {
            return boost::apply_visitor(lexical_cast_visitor<T>(), it->second);
        } else {
            return T();
        }
    }
    template<typename T>
    T get_attr(const std::string& attr, T value) const {
        const auto it = attributes.find(attr);
        if (it != attributes.end()) {
            return boost::apply_visitor(lexical_cast_visitor<T>(), it->second);
        } else {
            return value;
        }
    }

    const std::map<std::string, variant>& get_attrs() const { return attributes; }
private:
    std::map<std::string,variant> attributes;
};

/* specialization returning variant directly */
template<> inline annotated_cseq::variant
annotated_cseq::get_attr<annotated_cseq::variant>(const std::string& attr) const {
    return attributes.at(attr);
}


typedef annotated_cseq cseq;

} // namespace sina

#include "cseq_impl.h"

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

