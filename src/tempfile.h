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

#ifndef _TEMPFILE_H_
#define _TEMPFILE_H_

#include <boost/filesystem.hpp>

namespace sina {

class TempFile : boost::noncopyable {
public:
    TempFile(const boost::filesystem::path model = "sina-%%%%-%%%%-%%%%")
        : _path(boost::filesystem::temp_directory_path()
                / boost::filesystem::unique_path(model))
    {}
    ~TempFile() {
        boost::filesystem::remove(_path);
    }
    operator const boost::filesystem::path() const {
        return _path;
    }
    const boost::filesystem::path& path() const {
        return _path;
    }
    void dump(std::ostream& out) const {
        out << "Dumping Tempfile " << this << std::endl;
        std::string sep = "----------------------";
        out << sep << std::endl;
        boost::filesystem::ifstream file(*this);
        std::string line;
        while (getline(file, line).good()) {
            out << line << std::endl;
        }
        out << sep << std::endl;
    }
    std::string load() const {
        auto closer = [](FILE* fp){ if(fp) std::fclose(fp);};
        auto fp = std::unique_ptr<FILE, decltype(closer)>(
            std::fopen(_path.native().c_str(), "r"), closer);
        std::string res;
        res.resize(boost::filesystem::file_size(_path));
        auto len = std::fread(const_cast<char*>(res.data()), 1, res.size(), fp.get());
        return res;
    }
private:
    boost::filesystem::path _path;
};
std::ostream& operator<<(std::ostream& out, const TempFile& t) { return out << t.path(); }

} // namespace sina

#endif // _TEMPFILE_H_

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
