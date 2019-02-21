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

#ifndef _PROGRESS_H_
#define _PROGRESS_H_

#include <mutex>
#include <unistd.h>
#include <sys/ioctl.h>

#include "spdlog/spdlog.h"
#include "spdlog/fmt/bundled/chrono.h"


namespace sina {

struct print_n {
    print_n(size_t n, const char *sym) : _n(n), _sym(sym) {}
    friend std::ostream& operator<<(std::ostream& out, const print_n& p) {
        for (size_t i=0; i<p._n; i++) out << p._sym;
        return out;
    }
    size_t _n;
    const char *_sym;
};

static const char* bar_syms_unicode[] = {
    " ",
    "\xE2\x96\x8F", "\xE2\x96\x8E", "\xE2\x96\x8D", "\xE2\x96\x8C",
    "\xE2\x96\x8B", "\xE2\x96\x8A", "\xE2\x96\x89", "\xE2\x96\x88"
};
static const char* bar_syms_ascii[] = {
    " ", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "#"
};


class Progress {
    using clock_t = std::chrono::steady_clock;
    using timepoint_t = clock_t::time_point;
    using duration_t = clock_t::duration;

public:
    Progress(std::string desc="", unsigned int total=0, bool ascii=false,
             FILE* file=stderr, unsigned int width=0)
        : _file(file),
          _ncols(width),
          _total(total),
          _desc(desc),
          _bar_syms(ascii?bar_syms_ascii:bar_syms_unicode),
          _nsyms(ascii?std::extent<decltype(bar_syms_ascii)>::value
                 :std::extent<decltype(bar_syms_unicode)>::value)
    {
        if (_ncols == 0) {
            update_term_width();
        }
        show_progress();
        print(term_move_up);
    }

    void restart(std::string desc="", unsigned int total=0) {
        _n = 0;
        _total = total;
        _desc = desc;
        show_progress();
        print(term_move_up);
    }

    unsigned int count() {
        return _n;
    }

    void update_term_width() {
        int fd = fileno(_file);
        struct winsize size;
        if (ioctl(fd, TIOCGWINSZ, &size) == 0) {
            _ncols = size.ws_col;
        }
    }

    void format_bar_to(fmt::memory_buffer& buf, unsigned int width, float frac) {
        if (width == 0) {
            return;
        }
        buf.reserve(buf.size() + width * 3);
        auto it = std::back_inserter(buf);

        size_t complete    = frac * width * _nsyms;
        size_t full_blocks = complete / _nsyms;
        size_t frac_block  = complete % _nsyms;
        size_t fill_length = width - full_blocks;

        auto append_n = [&](size_t n, unsigned int idx) {
            size_t len_sym = strlen(_bar_syms[idx]);
            if (len_sym == 1) {
                std::fill_n(it, n, _bar_syms[idx][0]);
            } else {
                for (size_t i=0; i<n; i++)
                    std::copy_n(_bar_syms[idx], len_sym, it);
            }
        };

        append_n(full_blocks, _nsyms-1);
        if (frac_block) {
            append_n(1, frac_block);
            --fill_length;
        }
        append_n(fill_length, 0);
    }

    void update(unsigned int n=1) {
        _n += n;
        if (_n == _total) {
            std::lock_guard<std::mutex> lock(_mutex);
            show_progress();
            fflush(_file);
            return;
        }
        if (_n >= _last_print_n + _miniterations) {
            timepoint_t now = clock_t::now();
            duration_t delta_time = now - _last_update;
            if (delta_time > _mininterval) {
                std::lock_guard<std::mutex> lock(_mutex);
                _last_update = now;
                _miniterations = (_n - _last_print_n) * _mininterval / delta_time;
                show_progress(now);
                print(term_move_up);
                fflush(_file);
            }
        }
    }

    Progress& operator++() {
        update();
        return *this;
    }

    void operator+=(unsigned int n) {
        update(n);
    }

    void show_progress(timepoint_t now=clock_t::now()) {
        _last_print_n = _n;
        float frac = (float) _n / _total;
        auto elapsed = now - _started_at;
        auto eta =  elapsed * (1/frac -1);
        auto remaining = (frac > 0) ? elapsed * (1/frac - 1) : duration_t(0);

        fmt::memory_buffer left, right;

        auto arg_desc    = fmt::arg("desc", _desc);
        float percent    = frac * 100;
        auto arg_frac    = fmt::arg("frac", percent);
        auto arg_n       = fmt::arg("n", _n);
        auto arg_total   = fmt::arg("total", _total);
        auto arg_elapsed = fmt::arg("elapsed", elapsed);
        auto arg_remain  = fmt::arg("remaining", remaining);
        auto arg_eol     = fmt::arg("eol", term_eol);
        auto args = fmt::make_format_args(arg_desc, arg_frac, arg_n, arg_total,
                                          arg_elapsed, arg_remain, arg_eol);

        fmt::vformat_to(left, lbar_fmt, args);
        fmt::vformat_to(right, rbar_fmt, args);

        int space_for_bar = _ncols - left.size() - right.size() + term_eol.size();
        if (space_for_bar > 0) {
            format_bar_to(left, space_for_bar, frac);
        }
        left.reserve(left.size() + right.size());
        std::copy(right.begin(), right.end(), std::back_inserter(left));

        fwrite(left.data(), 1, left.size(), _file);
    }

    void print(const std::string &str) {
        fwrite(str.data(), sizeof(char), str.size(), _file);
    }

private:
    FILE* _file;
    unsigned int _ncols;
    std::atomic<unsigned int> _n{0};
    unsigned int _last_print_n{0};
    unsigned int _total;
    std::string _desc;
    const char **_bar_syms;
    unsigned int _nsyms;
    timepoint_t _started_at{clock_t::now()};
    timepoint_t _last_update{std::chrono::seconds(0)};
    duration_t _mininterval{std::chrono::milliseconds(10)};
    unsigned int _miniterations{1};
    std::string _bar_tpl;
    std::mutex _mutex;

    const std::string term_move_up = "\x1B[A";
    const std::string term_erase_line = "\x1B[0K";
    const std::string term_eol = "\n";
    const std::string lbar_fmt = "{desc}: {frac:3.0f}% |";
    const std::string rbar_fmt = "| {n}/{total} [{elapsed:%T} / {remaining:%T}]{eol}";
};


} // namespace sina


#endif // _LOG_H_
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
