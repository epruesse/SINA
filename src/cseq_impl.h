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

#ifndef _CSEQ_IMPL_H_
#define _CSEQ_IMPL_H_

namespace sina {

inline cseq& 
cseq::append(const std::string s) {
  return append(s.c_str());
}

template<typename T>
class lexical_cast_visitor
    : public boost::static_visitor<T> {
public:
    template<typename S>
    const T operator()(const S& s) const {
        try {
            return boost::lexical_cast<T>(s);
        } catch (...) {
            return T();
        }
    }

    // no convert for identical type

    const T operator()(const T& t) const {
        return t;
    }

    // return default value for vector (no sense to convert)
    const T operator()(const std::vector<cseq>& /*v*/) const {
        return T();
    }
};

inline cseq::iterator 
prev_begin(const cseq& c, const cseq::iterator& it) {
  if (c.begin() != it) {
    return it-1;
  }
  return it;
}
inline cseq::iterator 
prev_end(const cseq& /*c*/, const cseq::iterator& it) {
  return it;
}
inline bool
has_prev(const cseq& c, const cseq::iterator& it) {
  return c.begin() != it;
}

inline cseq::iterator
next_begin(const cseq& /*c*/, const cseq::iterator& it) {
  return it+1;
}
inline cseq::iterator
next_end(const cseq& c, const cseq::iterator& it) {
  if (it+1 == c.end()) {
    return it + 1;
  }
  return it + 2;
}
inline bool
has_next(const cseq& c, const cseq::iterator& it) {
  return (it+1) != c.end();
}


inline cseq::idx_type
get_node_id(cseq& c, const cseq::iterator& it){
    return it - c.begin();
}

inline cseq::iterator
cseq::begin() {
    return {bases.begin()};
}

inline cseq::const_iterator
cseq::begin() const {
    return {bases.begin()};
}

inline cseq::const_reverse_iterator
cseq::rbegin() const {
    return {bases.rbegin()};
}


inline cseq::iterator
cseq::end() {
    return {bases.end()};
}

inline cseq::const_iterator
cseq::end() const {
    return {bases.end()};
}
inline cseq::const_reverse_iterator
cseq::rend() const {
    return {bases.rend()};
}

inline cseq::iterator
cseq::getIterator(cseq::vidx_type i) {
  // this is weird. FIXME
    return {std::lower_bound(bases.begin(),
                                bases.end(),aligned_base(i,'.'))};
}

inline cseq::const_iterator
cseq::getIterator(cseq::vidx_type i) const {
    return {lower_bound(bases.begin(),
                                      bases.end(),aligned_base(i,'.'))};
}


template<typename T>
void
for_each_prev(cseq::iterator cit, T t) {
    t(*(cit-1));
}

} // namespace sina

#endif // _CSEQ_IMPL_H_

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
