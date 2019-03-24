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

#ifndef _CACHE_H_
#define _CACHE_H_

#include <mutex>
#include <list>
#include <unordered_map>

namespace sina {

template<typename KEY, typename VALUE>
class fifo_cache {
public:
    explicit fifo_cache(size_t size) : _size(size) {}
    void store(KEY key, VALUE&& value) {
        std::lock_guard<std::mutex> lock(_mutex);
        auto it = _keys.find(key);
        if (it != _keys.end()) {
            _items.erase(it->second);
            _keys.erase(it);
        }
        _items.emplace_front(std::make_pair(KEY(key), value));
        _keys.emplace(key, _items.begin());
        if (_keys.size() > _size) {
            _keys.erase(_items.back().first);
            _items.pop_back();
        }
    }

    bool try_get(const KEY& key, VALUE& val) {
        std::lock_guard<std::mutex> lock(_mutex);
        auto it = _keys.find(key);
        if (it == _keys.end()) {
            return false;
        }
        val = std::move(it->second->second);
        _items.erase(it->second);
        _keys.erase(it);
        return true;
    }

private:
    size_t _size;
    using list_type = std::list<std::pair<KEY, VALUE>>;
    list_type _items;
    std::unordered_map<KEY, typename list_type::iterator> _keys;
    std::mutex _mutex;
};



} // namespace SINA

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

