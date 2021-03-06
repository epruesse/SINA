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

#include <sys/time.h>
#include <ctime>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <iterator>
#include <string>
#include <iostream>
#include <sstream>
#include <numeric>
#include <thread>

#include <tbb/concurrent_unordered_map.h>

#ifndef _TIMER_H_
#define _TIMER_H_

/* defined in time.h 
struct timeval {
  time_t         tv_sec;   // seconds
  suseconds_t    tv_usec;  // microseconds
};
*/

namespace sina {

struct timestamp : private timeval {
    timestamp(int /*unused*/) {
        tv_sec = 0;
        tv_usec = 0;
    }

    void get() {
        // Adjusting the system clock breaks things
        // should use clock_gettime(CLOCK_MONOTONIC, ...)
        // (but that needs macos 10.12)
        gettimeofday(this, nullptr);
    }

    timestamp() {
        get();
    }

    operator float() {
        return tv_sec + float(tv_usec)/1000000;
    }
    
    timestamp operator-(const timestamp& rval) {
        timestamp result(0);
        if (tv_usec < rval.tv_usec) {
            result.tv_sec = tv_sec - rval.tv_sec - 1;
            result.tv_usec = 1000000 + tv_usec - rval.tv_usec;
        } else {
            result.tv_sec =  tv_sec - rval.tv_sec;
            result.tv_usec = tv_usec - rval.tv_usec;
        }
        return result;
    }

    timestamp operator+(const timestamp& rval) {
        timestamp result(0);
        if (tv_usec + rval.tv_usec >= 1000000) {
            result.tv_sec = tv_sec + rval.tv_sec + 1;
            result.tv_usec = tv_usec + rval.tv_usec - 1000000;
        } else {
            result.tv_sec = tv_sec + rval.tv_sec;
            result.tv_usec = tv_usec +  rval.tv_usec;
        }
        return result;
    }

    timestamp& operator+=(const timestamp& rval) {
        if (tv_usec + rval.tv_usec >= 1000000) {
            tv_sec += rval.tv_sec + 1;
            tv_usec = tv_usec + rval.tv_usec - 1000000;
        } else {
            tv_sec += rval.tv_sec;
            tv_usec += rval.tv_usec;
        }
        return *this;
    }

    friend std::ostream& operator<<(std::ostream& out, const timestamp& t) {
        return out << t.tv_sec
                   << "." << std::setfill('0') << std::setw(3) << t.tv_usec/1000
                   << "s";
    }
};



class timer {
    std::vector<timestamp> timings;
    std::vector<const char*> names;
    std::vector<timestamp>::iterator time_it;
    timestamp t_last;
    unsigned int calls{0};
public:
    timer() : timings(1, 0), t_last(0) {}

    void start() {
        time_it = timings.begin();
        t_last.get();
        ++calls;
    }

    void stop(const char* name=nullptr) {
        timestamp t_now;
        if (++time_it == timings.end()) {
            names.push_back(name);
            timings.emplace_back(0);
            time_it = timings.end() - 1;
        }
        *time_it +=  t_now - t_last;
        t_last.get();
    }

    void end_loop(int i) {
        time_it-=i;
        t_last.get();
    }

    timer& operator+=(const timer& o) {
        if (timings.size() != o.timings.size()) {
            throw std::runtime_error("Tried to add incompatible timers");
        }
        for (size_t i = 0; i < timings.size(); i++) {
            timings[i] += o.timings[i];
        }
        calls += o.calls;
        return *this;
    }

    friend std::ostream& operator<<(std::ostream& out , const timer& t) {
        out << std::accumulate(t.timings.begin(), t.timings.end(), timestamp(0))
            << " (" << t.calls << " calls, ";
        for (size_t i = 0; i < t.names.size(); ++i) {
            if (i > 0) {
                out << ", ";
            }
            if (t.names[i] != nullptr) {
                out << t.names[i];
            } else {
                out << i;
            }
            out << ": " << t.timings[i];
        }
        out << ")";
        return out;
    }
};

class timer_mt {
    tbb::concurrent_unordered_map<std::thread::id, timer, std::hash<std::thread::id>> timers;
public:
    timer& get_timer() {
        return timers[std::this_thread::get_id()];
    }

    friend std::ostream& operator<<(std::ostream& out , const timer_mt& mt) {
        auto it = mt.timers.begin();
        auto end = mt.timers.end();
        if (it == end) {
            out << "never called";
        } else {
            timer sum((it++)->second);
            for (; it!=end; ++it) {
                sum += it->second;
            }
            out << sum;
        }
        return out;
    }
};

} // namespace sina
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

