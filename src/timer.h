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

#include <iostream>
#include <sys/time.h>
#include <time.h>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <iterator>

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
    timestamp(int) {
       tv_sec = 0;
       tv_usec = 0;
    }

    void get() {
        gettimeofday(this,0);
    }

    timestamp() {
        get();
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
        return out << t.tv_sec << "." << std::setfill('0') << std::setw(6) << t.tv_usec;
    }
};



class timer {
    bool first_run;
    std::vector<timestamp> timestamps;
    std::vector<timestamp>::iterator time_it;
    timestamp t_last;

public:
    timer() : first_run(true), timestamps(1,0), t_last(0) {}

    void start() {
        if (timestamps.size() > 1 )  first_run=false;
        time_it = timestamps.begin();
        t_last.get();
    }

    void stop() {
        timestamp t_now;
        *time_it +=  t_now - t_last;
        ++time_it;
        if (time_it == timestamps.end()) {
            timestamps.push_back(timestamp(0));
            time_it = timestamps.end()-1;
        }
        t_last.get();
    }

    void end_loop(int i) {
        time_it-=i;
        t_last.get();
    }
    friend std::ostream& operator<<(std::ostream& out , const timer& t) {
        out << "timed " << t.timestamps.size() << ":";
        std::copy(t.timestamps.begin(), --t.timestamps.end(), 
                  std::ostream_iterator<timestamp>(out," "));
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

