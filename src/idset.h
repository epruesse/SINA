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

#include <cstdint>
#include <cstddef>
#include <vector>

#include "helpers.h"

/** Abstract base class for container of ids
 *
 */
class idset {
public:
    typedef uint32_t value_type;
    typedef std::vector<uint16_t> inc_t;
    typedef std::vector<uint8_t> data_t;

    /* virtual destructor */
    virtual ~idset() {}


    size_t size() const {
        return data.size();
    }
    

    /* add id to container; MUST be monotonically rising! */
    virtual void push_back(value_type n) = 0;
    
    /* increment vector at offsets included in set */
    virtual void increment(inc_t& data) const = 0;

    ///// for unit testing: ////

    /* create a new instance */
    virtual idset* make_new(value_type size) const = 0;

    /* get name of subclass */
    virtual const char* name() const { return "idset"; }

protected:
    mutable data_t data;
};


/** idset implementation as bitmap **/
class bitmap : public idset {
public:
    typedef typename data_t::value_type block_type;
    const size_t bits_per_block = sizeof(block_type) * 8;

    bitmap(value_type maxid) {
        int blocks_needed = (maxid + bits_per_block - 1 )/ bits_per_block;
        data.resize(blocks_needed, 0);
    }

    size_t block_index(value_type id) const {
        return id / bits_per_block;
    }

    size_t block_offset(value_type id) const {
        return id % bits_per_block;
    }

    /* set bit at @id */
    void set(value_type id) {
        data[block_index(id)] |= 1 << (block_offset(id));
    }

    /* get value of bit at @id */
    bool get(value_type id) const {
        return data[block_index(id)] & (1 << (block_offset(id)));
    }

    /* count total number of set bits */
    value_type count() const {
        value_type total = 0;
        for (const auto& block : data) {
            total += __builtin_popcount(block);
        }
        return total;
    }

    virtual idset* make_new(value_type size) const override {
        return new bitmap(size);
    }
    
    virtual const char* name() const override {
        return "bitmap";
    }
    
    virtual void push_back(value_type id) override {
        set(id);
    }


    virtual void increment(inc_t& t) const override {
        for (value_type i=0; i < data.size(); i++) {
            block_type block = data[i];
            while (block != 0) {
                auto j = __builtin_ctz(block);
                ++t[i*bits_per_block+j];
                block ^= 1<<j;
            }
        }
    }
};

/* idset implementation using variable length integers
 * (encoded using most significant bit to indicate need for another byte)
 */
class vlimap_abs : public idset {
public:
    vlimap_abs(int) {}

    /* once-forward iterator over contents */
    class const_iterator {
        data_t::const_iterator _it;
    public:
        const_iterator(const data_t::const_iterator& it) : _it(it) {}
        inline value_type operator*()  __attribute__((always_inline)) { // DO NOT CALL TWICE
            value_type val;
            // first byte
            uint8_t byte = *_it;
#define SINA_UNROLL
#ifdef SINA_UNROLL
            if (!(byte & 0x80)) {
                return byte;
            }
            // second byte
            val = byte - 0x80;
            byte = *(++_it);
            val += byte << 7;
            if (!(byte & 0x80)) {
                return val;
            }
            // third byte
            val -= 0x80 << 7;
            byte = *(++_it);
            val += byte << 14;
            if (!(byte & 0x80)) {
                return val;
            }
            // fourth byte
            val -= 0x80 << 14;
            byte = *(++_it);
            val += byte << 21;
            if (!(byte & 0x80)) {
                return val;
            }
            // fifth byte
            val -= 0x80 << 21;
            byte = *(++_it);
            val += byte << 28;
            //if (unlikely(byte & 0x80)) {
            //    throw something
            //}

            return val;
#else // SINA_UNROLL undefined
            if (likely(byte < 128)) {
                val = byte;
            } else {
                val = byte & 0x7f;
                unsigned int shift = 7;
                do {
                    byte = *(++_it);
                    val |= value_type(byte & 0x7f) << shift;
                    shift += 7;
                } while (byte >= 128);
            }
            return val;
#endif // SINA_UNROLL
        }
        const_iterator& operator++() {
            _it++; // step to next encoded value
            return *this;
        }
        bool operator!=(const const_iterator& rhs) const {
            return _it != rhs._it;
        }
    };
    
    const_iterator begin() const {
        return const_iterator(data.begin());
    }
    
    const_iterator end() const {
        return const_iterator(data.end());
    }
    
    size_t size() {
        return data.size();
    }

    virtual idset* make_new(value_type size) const override {
        return new vlimap_abs(size);
    }
    
    virtual const char* name() const override {
        return "vlimap_abs";
    }
    
    virtual void push_back(value_type n) override {
        while (n > 127) {
            data.push_back(n | 0x80);
            n >>= 7;
        }
        data.push_back(n);
    }
    
    virtual void increment(inc_t& t) const override {
        for (const_iterator it = begin(); it != end(); ++it) {
            t[*it]++;
        }
    }
};

/* idset implementation using variable length integer on distances */
class vlimap : public vlimap_abs {
public:
    vlimap(value_type size) : vlimap_abs(size), last(0) {}
    
    virtual idset* make_new(value_type size) const override {
        return new vlimap(size);
    }
    
    virtual const char* name() const override {
        return "vlimap";
    }
    
    virtual void push_back(value_type n) override {
        vlimap_abs::push_back(n-last);
        last = n;
    }
    
    virtual void increment(inc_t& t) const override {
        value_type last = 0;
        for (const_iterator it = begin(); it != end(); ++it) {
            last += *it;
            ++t[last];
        }
    }
private:
    value_type last;
};

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
