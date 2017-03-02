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

// Copyright 2011 Elmar Pruesse

#ifndef _PSEQ_H_
#define _PSEQ_H_

#include "cseq.h"
#include "aligned_base.h"

#include <vector>
#include <sstream>

namespace sina {

class base_profile {
public:
  enum base_types {
    BASE_A=0,
    BASE_G=1,
    BASE_C=2,
    BASE_TU=3,
    BASE_MAX=4,
  };

  base_profile(); // no default

  base_profile(int a, int g, int c, int t,
               int open, int extend) {
    int sum = a + g + c + t + open + extend;
    bases[BASE_A] = (float)a/sum;
    bases[BASE_G] = (float)g/sum;
    bases[BASE_C] = (float)c/sum;
    bases[BASE_TU]= (float)t/sum;
    gapOpen = (float)open/sum;
    gapExtend = (float)extend/sum;
  }

  base_profile(const base_iupac& b)
    : gapOpen(0), gapExtend(0)
  {
    if (b.ambig_order() > 0) {
      float val = 1.f / b.ambig_order();
      for (int i=0; i<BASE_MAX; i++) {
        bases[i]=0.f; // c++0x would come in handy here...
      }
      if (b.has_A())  bases[BASE_A] = val;
      if (b.has_G())  bases[BASE_G] = val;
      if (b.has_C())  bases[BASE_C] = val;
      if (b.has_TU()) bases[BASE_TU] = val;
    }
  }

  void complement() {
    std::swap(bases[BASE_A],bases[BASE_TU]);
    std::swap(bases[BASE_G],bases[BASE_C]);
  }

  /*
  float getA() const { return bases[BASE_A]; }
  float getG() const { return bases[BASE_G]; }
  float getC() const { return bases[BASE_C]; }
  float getT() const { return bases[BASE_TU]; }
  */

  float comp(const base_profile& rhs, float match, float mismatch,
             float gap, float gap_ext) const {
    float res = 0;
    for (int i=0; i<BASE_MAX; i++) {
      for (int j=0; j<BASE_MAX; j++) {
        if (i==j) {
          res += match * bases[i] * rhs.bases[j];
        } else {
          res += mismatch * bases[i] * rhs.bases[j];
        }
      }
    }
    return res +  gap * gapOpen + gap_ext * gapExtend;
  }

  float comp(const base_iupac& base, float match, float mismatch,
             float gap, float gap_ext) const {
    return comp(base_profile(base), match, mismatch,
                gap, gap_ext);
  }

  std::string getString() {
    std::stringstream out;
    const char l[4]  = {'A','G','C','T'};
    for (int i=0; i<BASE_MAX; i++) {
        out << l[i] << bases[i] << ":";
    }

    return out.str();
  }
 private:
  float bases[4];
  float gapOpen;
  float gapExtend;
};

typedef aligned<base_profile> aligned_base_profile;

class pseq {
public:
  typedef unsigned int idx_type;
  class iterator;
  class const_iterator;
  typedef iterator pn_iterator;
  typedef aligned_base_profile value_type;

  pseq(std::vector<cseq>::iterator, std::vector<cseq>::iterator);

  idx_type size() const { return profile.size(); }
  idx_type getWidth() const { return width; }
  iterator begin();
  const_iterator begin() const;
  iterator end();
  const_iterator end() const;
  pn_iterator pn_first_begin();
  pn_iterator pn_first_end();
  pn_iterator pn_last_begin();
  pn_iterator pn_last_end();
  const aligned_base_profile& getById(idx_type i) const {return profile[i];}
  void sort() {}

  void print_graphviz(std::ostream& out, std::string name);

private:
  idx_type width;
  std::vector<aligned_base_profile> profile;
  std::vector<idx_type> positions;

  friend pseq::idx_type get_node_id(const pseq& p, const pseq::iterator& i);
};

class pseq::iterator
  : public std::vector<aligned_base_profile>::iterator
{
public:
  iterator(std::vector<aligned_base_profile>::iterator it)
    : std::vector<aligned_base_profile>::iterator(it) {}
  iterator() : std::vector<aligned_base_profile>::iterator() {}

  typedef iterator pn_iterator;
  iterator prev_begin() const { iterator n(*this); return --n; }
  iterator prev_end() const { return (*this); }
  iterator next_begin() const { iterator n(*this); return ++n; }
  iterator next_end() const { iterator n(*this); return n+2; }
};

class pseq::const_iterator
  : public std::vector<aligned_base_profile>::const_iterator
{
public:
  const_iterator(std::vector<aligned_base_profile>::const_iterator it)
    : std::vector<aligned_base_profile>::const_iterator(it) {}
  const_iterator() : std::vector<aligned_base_profile>::const_iterator() {}

  typedef const_iterator const_pn_iterator;
  const_iterator prev_begin() const { const_iterator n(*this); return --n; }
  const_iterator prev_end() const { return (*this); }
  const_iterator next_begin() const { const_iterator n(*this); return ++n; }
  const_iterator next_end() const { const_iterator n(*this); return n+2; }
};

inline pseq::idx_type
get_node_id(pseq& p, const pseq::iterator& i)  {
  return i - p.begin();
}

inline pseq::iterator
pseq::begin() {
  return iterator(profile.begin());
}

inline pseq::iterator
pseq::end() {
  return iterator(profile.end());
}

inline pseq::const_iterator
pseq::begin() const {
  return const_iterator(profile.begin());
}

inline pseq::const_iterator
pseq::end() const {
  return const_iterator(profile.end());
}

inline  pseq::pn_iterator
pseq::pn_first_begin() {
  return begin();
}
inline  pseq::pn_iterator
pseq::pn_first_end() {
  return begin()+1;
}
inline  pseq::pn_iterator
pseq::pn_last_begin() {
  return end()-1;
}
inline  pseq::pn_iterator
pseq::pn_last_end() {
  return end(); 
}

inline pseq::iterator
prev_begin(const pseq& p, const pseq::iterator& it) {
  if (p.begin() != it) {
    return it-1;
  } else {
    return it;
  }
}
inline pseq::iterator
prev_end(const pseq& /*c*/, const pseq::iterator& it) {
  return it;
}
inline pseq::iterator
next_begin(const pseq& /*c*/, const pseq::iterator& it) {
  return it+1;
}
inline pseq::iterator
next_end(const pseq& p, const pseq::iterator& it) {
  if (it+1 == p.end()) {
    return it + 1;
  } else {
    return it + 2;
  }
}

std::ostream& operator<<(std::ostream& out, const sina::aligned_base_profile& ab);
} // namespace sina



#endif // _PSEQ_H_

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

