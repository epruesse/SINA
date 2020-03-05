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
#include <list>
#include <unistd.h>
#include <sys/ioctl.h>
#include <signal.h>

#include "spdlog/spdlog.h"
#include "spdlog/fmt/bundled/chrono.h"
#include "spdlog/sinks/ansicolor_sink.h"

namespace sina {

static const char* bar_syms_unicode[] = {
    " ",
    "\xE2\x96\x8F", "\xE2\x96\x8E", "\xE2\x96\x8D", "\xE2\x96\x8C",
    "\xE2\x96\x8B", "\xE2\x96\x8A", "\xE2\x96\x89", "\xE2\x96\x88"
};
static const char* bar_syms_ascii[] = {
    " ", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "#"
};


/* abstract base rendering progress bar
 * needs show_progress() filled in */
class base_progress {
public:
    using clock_t = std::chrono::steady_clock;
    using timepoint_t = clock_t::time_point;
    using duration_t = clock_t::duration;

    base_progress(std::string desc="", unsigned int total=0, bool ascii=false)
        : _total(total),
          _desc(desc),
          _bar_syms(ascii?bar_syms_ascii:bar_syms_unicode),
          _nsyms(ascii?std::extent<decltype(bar_syms_ascii)>::value
                 :std::extent<decltype(bar_syms_unicode)>::value)
    {}

    void restart(std::string desc="", unsigned int total=0) {
        _n = 0;
        _total = total;
        _desc = desc;
        show_progress();
    }

    unsigned int count() {
        return _n;
    }

    base_progress& operator++() {
        update();
        return *this;
    }

    void operator+=(unsigned int n) {
        update(n);
    }

    void set_total(unsigned int n) {
        _total = n;
    }

    unsigned int size() {
        return _total;
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
            }
        }
    }

    void render_progress(timepoint_t now, unsigned int width, fmt::memory_buffer& buf) {
        _last_print_n = _n;
        auto arg_desc    = fmt::arg("desc", _desc);
        auto elapsed = now - _started_at;
        auto arg_elapsed = fmt::arg("elapsed", elapsed);
        auto arg_eol     = fmt::arg("eol", term_eol);
        auto arg_n       = fmt::arg("n", _n);

        if (_total == 0) {
            fmt::format_to(buf, nototal_fmt, arg_desc, arg_elapsed, arg_eol, arg_n);
            return;
        }

        float frac = (float) _n / _total;
        auto eta =  elapsed * (1/frac -1);
        auto remaining = (frac > 0) ? elapsed * (1/frac - 1) : duration_t(0);

        fmt::memory_buffer right;

        float percent    = frac * 100;
        auto arg_frac    = fmt::arg("frac", percent);
        auto arg_total   = fmt::arg("total", _total);
        auto arg_remain  = fmt::arg("remaining", remaining);
        auto args = fmt::make_format_args(arg_desc, arg_frac, arg_n, arg_total,
                                          arg_elapsed, arg_remain, arg_eol);

        fmt::vformat_to(buf, lbar_fmt, args);
        fmt::vformat_to(right, rbar_fmt, args);

        int space_for_bar = width - buf.size() - right.size() + term_eol.size();
        if (space_for_bar > 0) {
            format_bar_to(buf, space_for_bar, frac);
        }
        buf.reserve(buf.size() + right.size());
        std::copy(right.begin(), right.end(), std::back_inserter(buf));
    }

    virtual void show_progress(timepoint_t now=clock_t::now()) = 0;
private:

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

    const std::string term_erase_line = "\x1B[0K";
    const std::string term_eol = "\n";
    const std::string lbar_fmt = "{desc}: {frac:3.0f}% |";
    const std::string rbar_fmt = "| {n}/{total} [{elapsed:%T} / {remaining:%T}]{eol}";
    const std::string nototal_fmt = "{desc}: {n} [{elapsed:%T}]{eol}";
};

/* Progress monitor writing directly to a file */
class Progress final : public base_progress {
public:
    Progress(std::string desc="", unsigned int total=0, bool ascii=false,
             FILE* file=stderr, unsigned int width=0)
        : base_progress(desc, total, ascii),
          _width(width),
          _file(file)

    {
        if (_width == 0) {
            update_term_width();
        }
    }

    void update_term_width() {
        int fd = fileno(_file);
        struct winsize size;
        if (ioctl(fd, TIOCGWINSZ, &size) == 0) {
            _width = size.ws_col;
        }
    }

    void show_progress(timepoint_t now=clock_t::now()) override final {
        // called from base_progress when it decided the progress bar needs
        // to be printed again
        fmt::memory_buffer buf;
        render_progress(now, _width, buf);
        std::copy(term_move_up.begin(), term_move_up.end(), std::back_inserter(buf));
        fwrite(buf.data(), 1, buf.size(), _file);
        fflush(_file);
    }
    const std::string term_move_up = "\x1B[A";
private:
    unsigned int _width;
    FILE* _file;
};

/******  spdlog integration ****/

class status_line;

class status_line_registry {
public:
    virtual ~status_line_registry() {}

    void add_status_line(status_line *msg) {
        _status_lines.push_back(msg);
        ++_change_since_last_print;
    }
    void remove_status_line(status_line *msg) {
        _status_lines.erase(
            std::remove(_status_lines.begin(), _status_lines.end(), msg),
            _status_lines.end());
        --_change_since_last_print;
    }
    const std::vector<status_line*>& get_status_lines() {
        return _status_lines;
    }
    int get_change_since_last_print() {
        return _change_since_last_print;
    }
    void clear_change_since_last_print() {
        _change_since_last_print = 0;
    }

    virtual void print_status_line_registry() = 0;

private:
    std::vector<status_line*> _status_lines;
    int _change_since_last_print{0};
};


class status_line {
public:
    status_line(std::shared_ptr<spdlog::logger> logger,
               spdlog::level::level_enum level)
        : _logger(logger), _level(level)
    {
        /* on creation, register it with all sinks that understand about us */
        for (auto &sink : _logger->sinks()) {
            auto ptr = dynamic_cast<status_line_registry*>(sink.get());
            if (ptr) {
                ptr->add_status_line(this);
            }
        }
    }

    ~status_line() {
        /* on deletion, deregister from all sinks */
        for (auto &sink : _logger->sinks()) {
            auto ptr = dynamic_cast<status_line_registry*>(sink.get());
            if (ptr) {
                ptr->remove_status_line(this);
            }
        }
    }

    /* trigger (re)print of this status line */
    void trigger_print_status_lines() {
        for (auto &sink : _logger->sinks()) {
            auto ptr = dynamic_cast<status_line_registry*>(sink.get());
            if (ptr) {
                ptr->print_status_line_registry();
            }
        }
    }

    virtual void render_status_line(fmt::memory_buffer &buf, unsigned int width) = 0;

    static const char* magic_filename() {
        static const char* file = "PROGRESS MONITOR";
        return file;
    }
    bool should_log() {
        return _logger->should_log(_level);
    }

    void logx(fmt::memory_buffer &buf) {
        using spdlog::details::fmt_helper::to_string_view;
        spdlog::source_loc loc;
        loc.filename = magic_filename();
        spdlog::details::log_msg msg{loc, &_logger->name(), _level, to_string_view(buf)};
        for (auto &sink : _logger->sinks()) {
            if (sink->should_log(_level)) {
                sink->log(msg);
            }
        }
    }

private:
    std::shared_ptr<spdlog::logger> _logger;
    spdlog::level::level_enum _level;
};


class sigwinch_mixin {
protected:
    sigwinch_mixin(bool install) {
        if (install) {
            install_handler();
        }
        instances().push_back(this);
    }
    virtual ~sigwinch_mixin() {
        instances().remove(this);
    }
    inline bool check_got_sigwinch() {
        if (got_sigwinch() != 0) {
            got_sigwinch() = 0;
            notify_instances();
        }
        if (_needs_update) {
            _needs_update = 0;
            return true;
        }
        return false;
    }
private:
    static sig_atomic_t& got_sigwinch() {
        static sig_atomic_t flag;
        return flag;
    }
    static std::list<sigwinch_mixin*>& instances() {
        static std::list<sigwinch_mixin*> list;
        return list;
    }
    static void handle_sigwinch(int) {
        got_sigwinch() = 1;
    }
    static void install_handler() {
        struct sigaction sa;
        memset(&sa, 0, sizeof(sa));
        if (sigaction(SIGWINCH, nullptr, &sa)) {
            return; // failed
        }
        if (sa.sa_handler == handle_sigwinch) {
            return; // already installed
        }
        memset(&sa, 0, sizeof(sa));
        sa.sa_handler = handle_sigwinch;
        sigfillset(&sa.sa_mask); // block other signals during handler
        if (sigaction(SIGWINCH, &sa, nullptr)) {
            return; // failed
        }
    }
    static void notify_instances() {
        //std::lock_guard<mutex_t> lock(super::mutex_);
        for (auto& instance : instances()) {
            instance->_needs_update = true;
        }
    }

    bool _needs_update{false};
};


template<typename TargetStream, class ConsoleMutex>
class terminal_sink final
    : public spdlog::sinks::ansicolor_sink<TargetStream, ConsoleMutex>,
      public status_line_registry, public sigwinch_mixin {
public:
    using super = spdlog::sinks::ansicolor_sink<TargetStream, ConsoleMutex>;
    using mutex_t = typename ConsoleMutex::mutex_t;

    terminal_sink() : super(), sigwinch_mixin(super::should_do_colors_) {
        if (super::should_do_colors_) {
            update_term_width();
        }
    }
    ~terminal_sink() override = default;

    // normal log message coming in
    void log(const spdlog::details::log_msg &msg) override final {
        _print(&msg);
    }

    // update to status coming in
    void print_status_line_registry() override final {
        _print(nullptr);
    }

    void _print(const spdlog::details::log_msg *msg) {
        if (not super::should_do_colors_) { // not a tty
            if (msg != nullptr) { // got a msg
                super::log(*msg); // pass it upstream
            }
            return;
        }

        if (check_got_sigwinch()) {
            update_term_width();
        }

        // render status lines
        fmt::memory_buffer status_text;
        for (auto line : get_status_lines()) {
            line->render_status_line(status_text, _ncols);
        }
        int nlines = get_status_lines().size() - get_change_since_last_print();
        clear_change_since_last_print();

        // move back up to last log line
        std::lock_guard<mutex_t> lock(super::mutex_);
        for (int i=0; i < nlines; ++i) {
            fwrite(term_move_up.data(), 1, term_move_up.size(), super::target_file_);
        }

        // print message if we have one and it's not an update for other sinks
        if (msg != nullptr && msg->source.filename != status_line::magic_filename()) {
            fwrite(term_clear_line.data(), 1, term_clear_line.size(), super::target_file_);
            super::log(*msg);
        }

        fwrite(status_text.data(), 1, status_text.size(), super::target_file_);
        fflush(super::target_file_);
    }

    void update_term_width() {
        int fd = fileno(super::target_file_);
        struct winsize size;
        if (ioctl(fd, TIOCGWINSZ, &size) == 0) {
            std::lock_guard<mutex_t> lock(super::mutex_);
            _ncols = size.ws_col;
        }
    }

    const std::string term_move_up = "\x1B[A";
    const std::string term_clear_line = "\x1B[K";
private:
    unsigned int _ncols;
};

using terminal_stdout_sink_mt = terminal_sink<spdlog::details::console_stdout, spdlog::details::console_mutex>;
using terminal_stdout_sink_st = terminal_sink<spdlog::details::console_stdout, spdlog::details::console_nullmutex>;
using terminal_stderr_sink_mt = terminal_sink<spdlog::details::console_stderr, spdlog::details::console_mutex>;
using terminal_stderr_sink_st = terminal_sink<spdlog::details::console_stderr, spdlog::details::console_nullmutex>;

template<typename Factory = spdlog::default_factory>
inline std::shared_ptr<spdlog::logger> stdout_terminal_mt(const std::string &logger_name)
{
    return Factory::template create<terminal_stdout_sink_mt>(logger_name);
}
template<typename Factory = spdlog::default_factory>
inline std::shared_ptr<spdlog::logger> stderr_terminal_mt(const std::string &logger_name)
{
    return Factory::template create<terminal_stderr_sink_mt>(logger_name);
}


class logger_progress final : public  base_progress, status_line {
public:
    logger_progress(std::shared_ptr<spdlog::logger> logger,
                    std::string desc="", unsigned int total=0, bool ascii=false,
                    unsigned int /*width*/=0,
                    spdlog::level::level_enum level=spdlog::level::warn
        )
        : base_progress(desc, total, ascii),
          status_line(logger, level)
    {
    }
    ~logger_progress() {
    }

    void show_progress(timepoint_t now=clock_t::now()) override final {
        if (!should_log()) return;
        duration_t delta_time = now - _last_update;
        status_line::trigger_print_status_lines();
        if (delta_time > _mininterval) {
            _last_update = now;

            fmt::memory_buffer buf;
            base_progress::render_progress(now, 60, buf);
            status_line::logx(buf);
        }
    }

    void render_status_line(fmt::memory_buffer &buf, unsigned int width) override final {
        base_progress::render_progress(clock_t::now(), width, buf);
    }

private:
    timepoint_t _last_update{std::chrono::seconds(0)};
    duration_t _mininterval{std::chrono::milliseconds(500)};

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
