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

#ifndef _BUFFER_H_
#define _BUFFER_H_

#include <tbb/scalable_allocator.h>

namespace sina {

template<typename T>
class buffer {
public:
    using size_type = size_t;
    explicit buffer(size_type size) {
        _start = (T*)scalable_malloc(sizeof(T) * size);
    }
    ~buffer() {
        scalable_free(_start);
    }
    T& operator[](size_type idx) {
        return _start[idx];
    }
private:
    T* _start;
};

template<typename T, size_t ALIGN=64>
class aligned_buffer {
public:
    using size_type = size_t;
    explicit aligned_buffer(size_type size) {
        size_type alloc = sizeof(T) * size + sizeof(offset_t) + ALIGN - 1;
        void *ptr = scalable_malloc(alloc);
        if (!ptr) { throw std::bad_alloc(); }
        size_t res = ((size_t)ptr + sizeof(offset_t) + ALIGN - 1) & ~(ALIGN -1);
        *((offset_t*)res - 1) = (offset_t)((size_t)res - (size_t)ptr);
        _start = (T*)res;
    }
    ~aligned_buffer() {
        offset_t offset = *((offset_t*)_start - 1);
        scalable_free((char*)_start - offset);
    }
    T& operator[](size_type idx) {
        return _start[idx];
    }
private:
    using offset_t = uint16_t;
    T* _start;
};

} // namespace sina

#endif // _BUFFER_H_

/*
  Local Variables:
  mode:c++
  c-file-style:"stroustrup"
  c-file-offsets:((innamespace . 0)(inline-open . 0)(case-label . 0))
  indent-tabs-mode:nil
  fill-column:99
  End:
*/
// vim: filetype=cpp:expandtab:shiftwidth=4:tabstop=8:softtabstop=4:encoding=utf-8:textwidth=99 :

