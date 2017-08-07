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

#ifndef _KMER_H_
#define _KMER_H_

#include <stdexcept> // for runtime_error
#include <unordered_set>
#include <vector>
#include <iostream>
#include "aligned_base.h"
#include "helpers.h" // for unlikely()


namespace sina {

namespace kmer {

/** converts base_iupac sequence to kmer
 */
class generator {
protected:
    const unsigned int _k;
    const unsigned int _mask;
    unsigned int _val;
    unsigned int _good_count;

public:
    generator(unsigned int k, unsigned int val=0)
        :  _k(k), _mask((1UL<<(2*k))-1), _val(val), _good_count(0)
    {
        if (_k < 1) {
            throw std::runtime_error("K must be at least 1");
        }
        if (sizeof(_val)*8 < 2UL*k) {
            throw std::runtime_error("K too large!");
        }
    }

    /** shift kmer by adding a base to the right
     *
     * use good() to check if the new window is a valid kmer
     */
    void push(const base_iupac& b) {
        if (unlikely(b.is_ambig())) {
            _good_count = 0;
        } else {
            _good_count++;
            _val <<= 2;
            _val &= _mask;
            _val += b.getBaseType();
        }
    }

    /** checks if the current sequence window is a valid kmer */
    bool good() const {
        return _good_count >= _k;
    }

    /** implicit conversion to int representation of kmer */
    operator unsigned int() const {
        return _val;
    }

    /** explicit converto to int representation of kmer */
    unsigned int val() const {
        return _val;
    }

    /** access good count **/
    unsigned int get_good_count() const {
        return _good_count;
    };

    /** pretty print self **/
    std::ostream& print_to(std::ostream& out) const {
        for (int i = _k-1; i >= 0; --i) {
            out << base_iupac( base_types( (_val >> (i*2)) & 3 ));
        }
        return out;
    }
};

/** filter kmers by prefix */
template<class generator>
class prefix_filter : public generator {
private:
    const unsigned int _p_mask;
    const unsigned int _p_val;
public:
    template<typename... ARGS>
    prefix_filter(unsigned int k, unsigned int p_len, unsigned int p_val, ARGS&&... data)
        : generator(k, std::forward<ARGS>(data)...),
          _p_mask(((1<<p_len*2)-1)<<((k-p_len)*2)),
          _p_val(p_val<<((k-p_len)*2)) {
    }
    bool good() const {
        return generator::good() && (generator::val() & _p_mask) == _p_val;
    }
};

/** filter all but first occurrence */
template<class generator>
class unique_filter : public generator {
private:
    std::unordered_set<unsigned int> _seen;
    bool is_good;
public:
    template<typename... ARGS>
    unique_filter(ARGS&&... data)
        : generator(std::forward<ARGS>(data)...),
          is_good(false) {
        // seen.reserve(?)
    }

    void push(const base_iupac& b) {
        generator::push(b);
        is_good = generator::good() && _seen.insert(generator::val()).second;
    }

    bool good() const {
        return is_good;
    }
};

/** provide container-like access for for-loops */
template<typename generator, typename bases>
class iterable : generator {
private:
    typedef typename bases::const_iterator bases_iterator;
    bases_iterator _begin, _end;
public:
    template<typename... ARGS>
    iterable(bases &v, ARGS&&... data)
        : generator(std::forward<ARGS>(data)...),
          _begin(v.begin()), _end(v.end())
    {}

    class iterator;
    iterator begin() {
        return iterator(*this, _begin, _end);
    }
    iterator end() {
        return iterator(*this, _end, _end);
    }

    class iterator {
    private:
        bases_iterator _begin, _end;
        iterable _generator;
    public:
        iterator(iterable& gen,
                 bases_iterator& begin, bases_iterator& end)
            : _generator(gen), _begin(begin), _end(end)
        {
            if (begin != end) {
                ++*this;
            }
        }

        iterator& operator++() {
            do {
                _generator.push(*_begin++);
            } while (not _generator.good() &&_begin != _end);
            return *this;
        }

        unsigned int operator*() const {
            return _generator;
        }

        bool operator!=(const iterator& rhs) const {
            return _begin != rhs._begin;
        }
    };
};


} // namespace sina::kmer

typedef kmer::generator kmer_generator;
typedef kmer::unique_filter<kmer_generator> unique_kmer_generator;
typedef kmer::prefix_filter<kmer_generator> prefix_kmer_generator;
typedef kmer::unique_filter<prefix_kmer_generator> unique_prefix_kmer_generator;

/**
 * iterate over all kmers in container v
 *
 * params:
 *  v   container to iterate over
 *  k   size of k
 *  val supply initial value (optional)
 */
template<typename bases, typename... ARGS>
kmer::iterable<kmer_generator, bases>
all_kmers(bases &v, ARGS&&... data) {
    return kmer::iterable<kmer_generator, bases>(v, std::forward<ARGS>(data)...);
}

/**
 * iterate over unique kmers in container v
 *
 * params:
 *  v   container to iterate over
 *  k   size of k
 *  val supply initial value (optional)
 */
template<typename bases, typename... ARGS>
kmer::iterable<unique_kmer_generator, bases>
unique_kmers(bases &v, ARGS&&... data) {
    return kmer::iterable<unique_kmer_generator, bases>(v, std::forward<ARGS>(data)...);
}

/**
 * iterate over kmers in container v with specific prefix
 *
 * params:
 *  v     container to iterate over
 *  k     size of k
 *  p_len prefix length
 *  p_val prefix value
 *  val   supply initial value (optional)
 */
template<typename bases, typename... ARGS>
kmer::iterable<prefix_kmer_generator, bases>
prefix_kmers(bases &v, ARGS&&... data) {
    return kmer::iterable<prefix_kmer_generator, bases>(v, std::forward<ARGS>(data)...);
}


/**
 * iterate over unique kmers in container v with specific prefix
 *
 * params:
 *  v     container to iterate over
 *  k     size of k
 *  p_len prefix length
 *  p_val prefix value
 *  val   supply initial value (optional)
 */
template<typename bases, typename... ARGS>
kmer::iterable<unique_prefix_kmer_generator, bases>
unique_prefix_kmers(bases &v, ARGS&&... data) {
    return kmer::iterable<unique_prefix_kmer_generator, bases>(v, std::forward<ARGS>(data)...);
}



} // namespace sina
#endif // _KMER_H_

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
