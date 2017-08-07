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

#include <cstdint>
#include <cstddef>
#include <vector>

#include "helpers.h"

/** Abstract base class for container of ids
 *
 */
class idset {
public:
    typedef uint32_t id_t;
    typedef std::vector<uint32_t> inc_t;
    typedef std::vector<uint8_t> data_t;

    /* virtual destructor */
    virtual ~idset() {}
    
    /* create a new instance */
    virtual idset* make_new(id_t size) const = 0;

    /* get name of subclass */
    virtual const char* name() const { return "idset"; }
    
    /* add id to container; MUST be monotonically rising! */
    virtual void push_back(id_t n) = 0;
    
    /* increment vector at offsets included in set */
    virtual void increment(inc_t& data) const = 0;
protected:
    data_t data;
};

/** idset implementation as bitmap **/
class bitmap : public idset {
public:
    typedef typename data_t::value_type block_t;
    const size_t block_size = sizeof(block_t);

    bitmap(id_t maxid) {
        data.resize((maxid+block_size-1)/block_size, 0);
    }

    /* set bit at @id */
    void set(id_t id) {
        data[id / block_size] |= 1 << (id % block_size);
    }

    /* get value of bit at @id */
    bool get(id_t id) const {
        return data[id / block_size] & (1 << (id % block_size));
    }

    /* count total number of set bits */
    id_t count() const {
        id_t total = 0;
        for (const auto& block : data) {
            total += __builtin_popcount(block);
        }
        return total;
    }

    virtual idset* make_new(id_t size) const override {
        return new bitmap(size);
    }
    
    virtual const char* name() const override {
        return "bitmap";
    }
    
    virtual void push_back(id_t id) override {
        set(id);
    }
    
    virtual void increment(inc_t& t) const override {
        for (id_t i=0; i < data.size(); i++) {
            block_t block = data[i];
            while (block != 0) {
                auto j = __builtin_ctz(block);
                ++t[i*block_size+j];
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
        id_t operator*() { // DO NOT CALL TWICE
            id_t val;
            uint8_t byte = *_it;
            if (likely(byte < 128)) {
                val = byte;
            } else {
                val = byte & 0x7f;
                unsigned int shift = 7;
                do {
                    byte = *(++_it);
                    val |= id_t(byte & 0x7f) << shift;
                    shift += 7;
                } while (byte >= 128);
            }
            return val;
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

    virtual idset* make_new(id_t size) const override {
        return new vlimap_abs(size);
    }
    
    virtual const char* name() const override {
        return "vlimap_abs";
    }
    
    virtual void push_back(id_t n) override {
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
    vlimap(id_t size) : vlimap_abs(size), last(0) {}
    
    virtual idset* make_new(id_t size) const override {
        return new vlimap(size);
    }
    
    virtual const char* name() const override {
        return "vlimap";
    }
    
    virtual void push_back(id_t n) override {
        vlimap_abs::push_back(n-last);
        last = n;
    }
    
    virtual void increment(inc_t& t) const override {
        id_t last = 0;
        for (const_iterator it = begin(); it != end(); ++it) {
            last += *it;
            t[last]++;
        }
    }
private:
    id_t last;
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
