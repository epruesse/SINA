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
#include <cassert>
#include <vector>
#include <iostream>

#include <tbb/scalable_allocator.h>
#include <tbb/cache_aligned_allocator.h>

#include "helpers.h"

/** Abstract base class for container of ids
 *
 */
class idset {
public:
    using value_type = uint32_t;
    // using inc_t = std::vector<int16_t, tbb::cache_aligned_allocator<int16_t>>;
    //using data_t = std::vector<uint8_t, tbb::cache_aligned_allocator<uint8_t>>;
    using inc_t = std::vector<int16_t>;
    using data_t = std::vector<uint8_t>;

    /* virtual destructor */
    virtual ~idset() = default;

    size_t size() const {
        return _size;
    }

    /* add id to container; MUST be monotonically rising! */
    virtual void push_back(value_type n) = 0;
    
    /* increment vector at offsets included in set */
    virtual int increment(inc_t& data) const = 0;

    ///// for unit testing: ////

    /* create a new instance */
    virtual idset* make_new(value_type size) const = 0;

    void shrink_to_fit() { data.shrink_to_fit(); }

    virtual std::ostream& write(std::ostream& o) const {return o;}

protected:
    mutable data_t data;
    size_t _size{0};
};


/** idset implementation as bitmap **/
class bitmap : public idset {
public:
    using block_type = typename data_t::value_type;
    const size_t bits_per_block = sizeof(block_type) * 8;

    explicit bitmap(value_type maxid) : idset() {
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
        return (data[block_index(id)] & (1 << (block_offset(id)))) != 0;
    }

    /* count total number of set bits */
    value_type count() const {
        value_type total = 0;
        for (const auto& block : data) {
            total += __builtin_popcount(block);
        }
        return total;
    }

    idset* make_new(value_type size) const override {
        return new bitmap(size);
    }
    
    void push_back(value_type id) override {
        if (!get(id)) ++_size;
        set(id);
    }


    int increment(inc_t& t) const override {
        for (value_type i=0; i < data.size(); i++) {
            block_type block = data[i];
            while (block != 0) {
                auto j = __builtin_ctz(block);
                ++t[i*bits_per_block+j];
                block ^= 1<<j;
            }
        }
        return 0;
    }
};



/* idset implementation using plain integers
 */
class imap_abs : public idset {
public:
    explicit imap_abs(int /*unused*/) : idset() {};

    /* once-forward iterator over contents */
    class const_iterator {
        data_t::const_iterator _it;
    public:
        explicit const_iterator(const data_t::const_iterator& it) : _it(it) {}
        value_type operator*() {
            return *(value_type*)(&*_it);
        }
        const_iterator& operator++() {
            _it += sizeof(value_type);
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

    idset* make_new(value_type size) const override {
        return new imap_abs(size);
    }

    void push_back(value_type n) override {
        size_t offset = data.size();
        data.resize(data.size() + sizeof(value_type));
        *(value_type*)(data.data() + offset) = n;
        ++_size;
    }

    int increment(inc_t& t) const override {
        for (idset::value_type it : *this) {
            ++t[it];
        }
        return 0;
    }
};



/* idset implementation using variable length integers
 * (encoded using most significant bit to indicate need for another byte)
 */
class vlimap_abs : public idset {
public:
    explicit vlimap_abs(int /*maxsize*/=0) : idset() {}

    /* once-forward iterator over contents */
    class const_iterator {
        data_t::const_iterator _it;
    public:
        data_t::const_iterator iter() {
            return _it;
        }

        explicit const_iterator(const data_t::const_iterator& it) : _it(it) {}
        inline value_type operator*()  __attribute__((always_inline)) { // DO NOT CALL TWICE
            value_type val;
            // first byte
            uint8_t byte = *_it;
//#define SINA_UNROLL
#ifdef SINA_UNROLL
            if ((byte & 0x80) == 0) {
                return byte;
            }
            // second byte
            val = byte - 0x80;
            byte = *(++_it);
            val += byte << 7;
            if ((byte & 0x80) == 0) {
                return val;
            }
            // third byte
            val -= 0x80 << 7;
            byte = *(++_it);
            val += byte << 14;
            if ((byte & 0x80) == 0) {
                return val;
            }
            // fourth byte
            val -= 0x80 << 14;
            byte = *(++_it);
            val += byte << 21;
            if ((byte & 0x80) == 0) {
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
        inline const_iterator& operator++()  __attribute__((always_inline))  {
            ++_it; // step to next encoded value
            return *this;
        }
        inline bool operator!=(const const_iterator& rhs) const  __attribute__((always_inline))  {
            return _it != rhs._it;
        }
    };
    
    const_iterator begin() const {
        return const_iterator(data.begin());
    }
    
    const_iterator end() const {
        return const_iterator(data.end());
    }
    
    idset* make_new(value_type size) const override {
        return new vlimap_abs(size);
    }
    
    void push_back(value_type n) override {
        while (n > 127) {
            data.push_back(n | 0x80);
            n >>= 7;
        }
        data.push_back(n);
        ++_size;
    }
    
    int increment(inc_t& t) const override {
        for (idset::value_type it : *this) {
            ++t[it];
        }
        return 0;
    }
};

/* idset implementation using variable length integer on distances */
class vlimap final : public vlimap_abs {
public:
    vlimap(value_type maxsize, inc_t::value_type inc=1)
        : vlimap_abs(maxsize),
          _inc(inc),
          _last(0),
          _maxsize(maxsize)
    {}
    
    idset* make_new(value_type maxsize) const override {
        return new vlimap(maxsize);
    }
    
    void push_back(value_type n) override {
        vlimap_abs::push_back(n-_last);
        _last = n;
    }

    inline int increment(inc_t& t) const override {
        auto it = data.begin();
        auto end = data.end();
        value_type last = 0;
        while (it != end) {
            uint8_t byte = *it;
            if (likely(byte < 128)) {
                last += byte;
            } else {
                value_type val = byte & 0x7f;
                unsigned int shift = 7;
                do {
                    byte = *(++it);
                    val |= value_type(byte & 0x7f) << shift;
                    shift += 7;
                } while (byte >= 128);
                last += val;
            }
            t[last] += _inc;
            ++it;
        }
        return 1-(_inc+1)/2;
    }

    /* append contents of another vlimap
     * - largest value in LHS must be smaller than smallest in RHS
     * - invertion status must be equal
     */
    void append(const vlimap& other) {
        if (other.data.size() == 0) {
            return;
        }
        if (data.size() == 0) {
            data = other.data;
            _last = other._last;
            return;
        }

        auto it = other.begin();
        auto val = *it;
        push_back(val);

        data_t::const_iterator data_it = (++it).iter();
        data_t::const_iterator end = other.data.end();
        data.insert(data.end(), data_it, end);
        _last = other._last;
    }

    /* invert map
     * replaces "in map" with "not in map"
     * increment() on inverted object will actually decrement and return 1 instead of 0
     */
    void invert() {
        vlimap res(0, -1);
        value_type next = 0, last = 0;
        for (auto inc : *this) {
            next += inc;
            while (last < next) {
                res.push_back(last++);
            }
            ++last;
        }
        while (last < _maxsize) {
            res.push_back(last++);
        }
        std::swap(data, res.data);
        std::swap(_inc, res._inc);
        std::swap(_last, res._last);
        std::swap(_maxsize, res._maxsize);
    }

    struct file_header {
        uint32_t inc, last, bytesize, size;
    };

    std::ostream& write(std::ostream& out) const override {
        file_header head;
        head.inc = _inc;
        head.last = _last;
        head.bytesize = data.size();
        head.size = size();
        out.write((char*)&head, sizeof(file_header));
        out.write((char*)data.data(), sizeof(data_t::value_type) * data.size());
        return out;
    }

    std::istream& read(std::istream& in) {
        file_header head;
        in.read((char*)&head, sizeof(file_header));
        _inc = head.inc;
        _last = head.last;
        _size = head.size;
        data.resize(head.bytesize);
        in.read((char*)data.data(), sizeof(data_t::value_type) * data.size());
        return in;
    }
private:
    inc_t::value_type _inc{1};
    value_type _last, _maxsize;
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
