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

#include "cseq.h"
#include "log.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <functional>
#include <cmath>
#include <zlib.h>

using namespace std;
using namespace sina;

static auto logger = Log::create_logger("cseq");

cseq::cseq(const char *_name, float _score, const char *_data) 
    : name(_name), score(_score)
{
    if (_data != nullptr) {
        append(_data);
    }
}

  
void
cseq::clearSequence() { 
    bases.clear(); 
    alignment_width = 0; 
}

cseq&
cseq::append(const char *str) {
    //FIXME: keep internal '.'s
    while(*str != 0) {
        if (*str != ' ' && *str != '\t' && *str != '\n' &&  *str != '\r') {
            if (*str != '-' && *str != '.') {
                bases.emplace_back(alignment_width, *str);
            }
            alignment_width++;
        }
        str++;
    }

    return *this;
}

cseq&
cseq::append(const aligned_base& ab) {
    if (ab.getPosition() >= alignment_width) {
        // it's allowed to have more than one base on the
        // same position. we fix them in one go later on.
        bases.push_back(ab);
        alignment_width = ab.getPosition();
    } else {
        // but the new base must not come _before_ the last one 
        logger->error("$ cseq::append(): wrong order! {}({}<{})",
                      ab.getBase(), ab.getPosition(), alignment_width);

        bases.emplace_back(alignment_width, ab.getBase());
    } 

    return *this;
}


void
cseq::setWidth(vidx_type newWidth) {
    if (bases.empty() || newWidth >= bases.back().getPosition() + 1) {
        // modify at will if changing only number of trailing gaps
        alignment_width = newWidth;
        return;
    }

    // we can't shrink to less than 0 gaps
    if (newWidth < size()) {
        logger->critical("Cannot shrink '{}' aligment width to {} - got {} bases",
                         getName(), newWidth, size());
        throw std::runtime_error(
                "Attempted to shrink alignment width below base count"
            );
    }

    // find the number of bases from the right where
    // <position-of-base> + <number-of-following-bases>
    // is at most <width-of-alignment>
    int skip;
    for(skip = 0; skip < size(); skip++) {
        if (bases[size() - skip - 1].getPosition() + skip < newWidth) {
            break;
        }
    }

    for (int i = skip; i > 0; --i) {
        bases[size() - i].setPosition(newWidth - i);
    }
    alignment_width = newWidth;

    logger->warn("moved last {} bases to shrink alignment to {} columns",
                 skip, alignment_width);
}


string
cseq::getAligned(bool nodots, bool dna) const {
    string aligned;
    aligned.reserve(alignment_width);

    char dot='.';
    if (nodots) {
        dot='-';
    }

    auto it = bases.begin(), it_end = bases.end();

    unsigned int cursor = 0;
    for (; it != it_end; ++it) {
        unsigned int pos = it->getPosition();

        // advance cursor by filling with gap character
        aligned.append(pos - cursor, dot);
        // (only the first "gap" is dots)
        dot = '-';
        cursor = pos;

        // print base
        if (dna) {
            aligned.append(1, it->getBase().iupac_dna());
        } else {
            aligned.append(1, it->getBase().iupac_rna());
        }
        cursor++;
    }

    if (cursor < alignment_width) {
        if (!nodots) {
            dot='.';
        }
        aligned.append(alignment_width - cursor, dot);
    }

    return aligned;
}

string
cseq::getBases() const {
    string basestr;
    basestr.reserve(bases.size());

    for (auto base : bases) {
        basestr += base.getBase();
    }

    return basestr;
}

struct compressed_data {
    char id;
    unsigned short size;
    unsigned char data[];
};

void
cseq::compressAligned(std::vector<unsigned char> &out) {
    vector<unsigned char> buf;
    using uint = unsigned int;

    bases.emplace_back(alignment_width);
    const uint bas = bases.size();

    const uint orig_size = sizeof(aligned_base) * bas;
    buf.resize(orig_size);

    for (uint i=0; i<bas; ++i) {
        buf[i]=bases[i].getBase();
    }

    idx_type last=0;
    for (uint i=0; i < bas; i++) {
        idx_type idx = bases[i].getPosition();
        idx_type diff = idx - last;
        for (uint j = 0; j < sizeof(idx_type); ++j) {
            buf[(j+1)*bas+i]  = (unsigned char)(diff &0xFF);
            diff >>=8;
        }
        last = idx;
    }

    unsigned long compr_size = compressBound(orig_size);
    out.resize(compr_size);
    auto *cd = reinterpret_cast<compressed_data*>(&out.front());

    cd->id='#';
    cd->size=(unsigned short)orig_size;

    compress2(cd->data, &compr_size, &buf.front(), orig_size,9);

    out.resize(compr_size+sizeof(compressed_data));
}

void
cseq::assignFromCompressed(const void* data, size_t len) {
    vector<unsigned char> buf;
    using uint = unsigned int;
    const auto *cd = reinterpret_cast<const compressed_data*>(data);
    buf.resize(cd->size);

    unsigned long compr_size = len - sizeof(compressed_data);
    unsigned long orig_size = cd->size;

    uncompress(&buf.front(), &orig_size, cd->data, compr_size);

    const uint bas = orig_size / sizeof(aligned_base);

    bases.clear();
    bases.reserve(bas);


    idx_type last = 0;
    for (uint i=0; i<bas; ++i) {
        idx_type diff = 0;

        for (int j = sizeof(idx_type); j != -1; --j) {
            diff <<= 8;
            diff |= buf[(j+1)*bas+i];
        }
        last+=diff;
        bases.emplace_back(last, buf[i]);
    }
    alignment_width = bases.back().getPosition();
    bases.pop_back();
}



char
cseq::operator[](cseq::vidx_type i) {
    vector<aligned_base>::const_iterator it = getIterator(i);
    if (it != bases.end() && i == it->getPosition()) {
        return it->getBase();
    }
    return '-';
}


std::ostream&
sina::operator<<(std::ostream& out, const cseq& c) {
    out << c.getName();
    return out;
}

void
cseq::reverse() {
    std::reverse(bases.begin(), bases.end());
    std::for_each(bases.begin(), bases.end(),
                  aligned_base_reverse_position(alignment_width-1));
}

void
cseq::complement() {
    std::for_each(bases.begin(), bases.end(),
                  std::mem_fun_ref(&aligned_base::complement));
}

void
cseq::upperCaseAll() {
    std::for_each(bases.begin(), bases.end(),
                  std::mem_fun_ref(&aligned_base::setUpperCase));
}


template<typename RandIt, typename StrictWeakOrdering>
RandIt group_by(RandIt& a, RandIt& b,
                StrictWeakOrdering & /*cmp*/) {
    sort(a,b);
    RandIt end = unique(a,b);
    return end;
}

template<typename RandIt>
typename RandIt::value_type group_by(RandIt& a, RandIt& b) {
    return group_by(a,b, less<typename RandIt::value_type>());
}

string color_code(const string& in) {
    string::const_iterator in_end = in.end();
    stringstream tmp;
    bool colored = false;
    for (string::const_iterator it = in.begin(); it != in_end; ++it) {
        switch(*it) {
                case 'a':
                case 'A':
                    tmp << "\033[34m";
                    colored = true;
                    break;
                case 'g':
                case 'G':
                    tmp << "\033[35m";
                    colored = true;
                    break;
                case 'c':
                case 'C':
                    tmp << "\033[32m";
                    colored = true;
                    break;
                case 't':
                case 'T':
                case 'u':
                case 'U':
                    tmp << "\033[33m";
                    colored = true;
                    break;
                default:
                    if (colored) {
                        tmp << "\033[0m";
                        colored = false;
                    }
        }
        tmp << *it;
    }
    if (colored) {
        tmp << "\033[0m";
    }
    return tmp.str();
}

void
cseq::write_alignment(std::ostream& ofs, std::vector<cseq>& seqs,
                      cseq::idx_type from_pos,
                      cseq::idx_type to_pos,
                      bool colors
                      ) {
    if (seqs.empty()) {
        ofs << "cseq::write_alignment(): no sequences?" << endl;
        return;
    }
    if (from_pos > to_pos || to_pos >= seqs[0].getWidth()) {
        ofs << "cseq::write_alignment(): range out of bounds!" << endl;
        return;
    }

    vector<string> out(seqs.size());

    char outchar[seqs.size()];
    auto jmax = seqs.size();

    for (auto i = from_pos; i <= to_pos; ++i) {
        bool gap = true;
        for (auto j = 0; j < jmax; ++j) {
            outchar[j] = seqs[j][i];
            if (outchar[j] != '-') {
                gap = false;
            }
        }

        if (!gap || i == to_pos-1 ) {
            for (auto j = 0; j < jmax; ++j) {
                out[j].append(1, outchar[j]);
            }
        }
    }

    std::map<string,list<int> > mymap;

    size_t maxlen = 0;
    for (auto it = out.begin(); it != out.end(); ++it) {
        maxlen = std::max(maxlen, it->size());
        mymap[*it].push_back(it - out.begin());
    }

    ofs << "Dumping pos " << from_pos << " through " << to_pos << ":" << endl;
    for (unsigned int i = 0; i < maxlen; i+=100) {
        int len = maxlen - i;
        len = std::min(100,len);
        for (auto & it : mymap) {
            if (colors) {
                ofs << color_code(it.first.substr(i, len)) << " ";
            } else {
                ofs << it.first.substr(i, len) << " ";
            }
            bool range=false, is_last=false, is_secondlast=false;
            int last=-2;
            for (int jt : it.second) {
                if (range) {
                    if (jt != last+1) {
                        ofs << last;
                        range=false;
                        ofs << " " << jt;
                    }
                } else {
                    if (jt == last+1) {
                        range=true;
                        ofs << "-";
                    } else {
                        ofs << " " << jt;
                    }
                }
                last = jt;
                if (jt +1 == (int)seqs.size()) {
                    is_last = true;
                }
                if (jt +2 == (int)seqs.size()) {
                    is_secondlast = true;
                }
            }
            if (range) {
                ofs << last;
            }
            if (is_last) {
                ofs << " <---(## NEW ##) ";
            }
            if (is_secondlast) {
                ofs << " <---(%% ORIG %%) ";
            }

            ofs << endl;
        }
        ofs << endl;
    }
}

void
cseq::fix_duplicate_positions(std::ostream& log, bool lowercase, bool remove) {
    idx_type total_inserts = 0, longest_insert = 0, orig_inserts = 0;

    auto last_it = bases.begin();
    auto bases_end = bases.end();
    auto curr_it = last_it+1;
    idx_type last_idx = last_it->getPosition();


    if (remove) {
        log << "insertion=remove not implemented, using shift; ";
    }

    for (; curr_it < bases_end; ++curr_it) {
        idx_type curr_idx = curr_it->getPosition();
        // check for insertions
        if (last_idx == curr_idx) { 
            // no move -> curr is insertion
            if (curr_it+1 != bases_end) {
                // not at end of sequence -> da capo
                continue;
            }
            // curr is the last base.
            // move iterator to end and fall through to placement
            ++curr_it;
        }
        idx_type num_inserts = curr_it - last_it - 1;
        if (num_inserts == 0) {
            // no insertions. leave base untouched. 
            // remember last good position.
            last_it = curr_it;
            last_idx = curr_idx;
            continue;
        }

        // WE HAVE AN INSERTION :)

        // Determine the range where it may be placed:
        // first free position
        idx_type range_begin = last_it->getPosition() + 1; 
        // first occupied/invalid position
        idx_type range_end   = (curr_it == bases_end) ?
            alignment_width : curr_it->getPosition();

        // make last_it first base to be repositioned
        ++last_it; 
        // make curr_it last base to be repositioned
        --curr_it;

        orig_inserts = num_inserts;
        if (range_end - range_begin < num_inserts) { // range too small
            log << "shifting bases to fit in "<<num_inserts<<" bases at pos "<<range_begin<<" to "<<range_end<<";";
            // last_it and curr_it are the bases enclosing our insertion
            while (range_end - range_begin < num_inserts) {
                int next_left_gap;
                int next_right_gap;
                auto left = last_it;
                auto right = curr_it;
                
                // find first free gap to left of range
                if (left == bases.begin()) {
                    if (range_begin > 0) {
                        next_left_gap = range_begin - 1;
                    } else {
                        next_left_gap = -1;
                    }
                } else {
                    if ((left-1)->getPosition() + 1 < range_begin) {
                        next_left_gap = range_begin - 1;
                    } else {
                        --left;
                        while(left != bases.begin() && 
                              (left-1)->getPosition() + 1 >= left->getPosition()) {
                            --left;
                        }
                        next_left_gap = left->getPosition() - 1;
                    }
                }
                    
                // find first free gap to right of range

                if (right + 1 == bases.end()) {
                    if (range_end < alignment_width) {
                        next_right_gap = range_end;
                    } else {
                        next_right_gap = -1;
                    }
                } else {
                    if ((right+1)->getPosition() > range_end) {
                        next_right_gap = range_end;
                    } else {
                        ++right;
                        while(right+1 != bases.end() && 
                              (right)->getPosition() +1 >= (right+1)->getPosition()) {
                            ++right;
                        }
                        next_right_gap = right->getPosition() + 1;
                    }
                }


                if (next_right_gap == -1 || 
                    (next_left_gap != -1 && 
                     range_begin - next_left_gap <= next_right_gap - (range_end - 1)) ) {
                    if (next_left_gap == -1) {
                        throw runtime_error("ERROR: no space to left and right?? "
                                            "sequence longe than alignment?!");
                    }
                    num_inserts += last_it - left;
                    range_begin = next_left_gap;
                    last_it = left;
                } else {
                    num_inserts += right - curr_it;
                    range_end = next_right_gap + 1;
                    curr_it = right;
                }

            }
        } else {
            range_begin = range_end - num_inserts;
        }
        // make curr_it first base not to be positioned
        ++curr_it;

        for (; last_it != curr_it; ++last_it) {
            last_it->setPosition(range_begin++);
            if (lowercase) {
                last_it->setLowerCase();
            }
        }

        total_inserts+=num_inserts;
        longest_insert = std::max(longest_insert,num_inserts);

        last_it = curr_it;
        last_idx = curr_it->getPosition();
    }
    if (total_inserts > 0) {
        log << "total inserted bases=" << total_inserts << ";"
            << "longest insertion=" << longest_insert << ";"
            << "total inserted bases before shifting=" << orig_inserts << ";";
    }
}

std::vector<std::pair<unsigned int, unsigned int>>
cseq::find_differing_parts(const cseq& right) const {
    using bases_it = std::vector<aligned_base>::const_iterator;
    auto l_it = bases.begin(), l_end = bases.end();
    auto r_it = right.bases.begin(), r_end = right.bases.end();
    
    std::vector<std::pair<unsigned int, unsigned int>> result;
    int score = 0;
    bool bad = false;
    unsigned int start = 0;

    int lpos = l_it->getPosition();
    int rpos = r_it->getPosition();

    while(l_it != l_end && r_it != r_end) {
        if (lpos < rpos) {
            score = 4;
            ++l_it;
        } else if ( rpos < lpos ) {
            score = 4;
            ++r_it;
        } else { // rpos <=> lops
            if ((char)l_it->getBase() != (char)r_it->getBase()) {
                score = 4;
            }
            ++r_it;
            ++l_it;
        }
        if (l_it != l_end) {
            lpos = l_it->getPosition();
        }
        if (r_it != r_end) {
            rpos = r_it->getPosition();
        }
        if (score > 0) {
            if (!bad) {
                int rpos = std::max(right.bases.begin(), r_it - 2)->getPosition();
                start = std::min(lpos, rpos);
                bad=true;
            } else {
                if (--score <= 0 && lpos == rpos) {
                    result.push_back({start, lpos});
                    bad=false;
                }
            }
        }
    }
    if (bad) {
        result.push_back({start, std::min(lpos, rpos)});
    }
    
    return result;
}


float
cseq::calcPairScore(const std::vector<int>& pairs) {
    // create array to hold counts for base combinations
    std::vector<int> count;
    count.resize(65536);
    for (int i=0; i<65536; ++i) {
        count[i] = 0;
    }

    int num=0;

    for(unsigned int i=0; i<pairs.size(); ++i) {
        if (pairs[i] != 0) { // alignment position has "helix-parter"
            unsigned char left, right;
            left = operator[](i);
            right = operator[](pairs[i]);

            if ((left != '.' and right != '.')
                and
                (left != '-' or right != '-')) {
                num++;
                if (left < right) {
                    count[((int)(left)<<8) + right]++;
                } else {
                    count[((int)(right)<<8) + left]++;
                }
            }
        }
    }

#ifdef DEBUG
    stringstream detail;
    for (int i=0; i<65536; ++i) {
       if (count[i]>0) {
           detail << (char)(i>>8) << (char)(i&0xFF) << ": "
                  << count[i]/2 << "  ";
       }
    }
    if (detail.str() != "") {
        logger->info("bp detail: {}", detail.str());
    }
#endif


#if 1
    float score =
        count[((int)('-')<<8) + 'A'] * 0 +
        count[((int)('-')<<8) + 'C'] * 0 +
        count[((int)('-')<<8) + 'G'] * 0 +
        count[((int)('-')<<8) + 'U'] * 0 +
        count[((int)('A')<<8) + 'A'] * 0 +
        count[((int)('A')<<8) + 'C'] * 0 +
        count[((int)('A')<<8) + 'G'] * .5 +
        count[((int)('A')<<8) + 'U'] * 1.1 +
        count[((int)('C')<<8) + 'C'] * 0 +
        count[((int)('C')<<8) + 'G'] * 1.5 +
        count[((int)('C')<<8) + 'U'] * 0 +
        count[((int)('G')<<8) + 'G'] * .4 +
        count[((int)('G')<<8) + 'U'] * .9 +
        count[((int)('U')<<8) + 'U'] * 0;

#else
  float score =
        count[((int)('-')<<8) + 'A'] * -1 +
        count[((int)('-')<<8) + 'C'] * -1 +
        count[((int)('-')<<8) + 'G'] * -1 +
        count[((int)('-')<<8) + 'U'] * -1 +
        count[((int)('A')<<8) + 'A'] * -1 +
        count[((int)('A')<<8) + 'C'] * -1 +
        count[((int)('A')<<8) + 'G'] * 0 +
        count[((int)('A')<<8) + 'U'] * 1 +
        count[((int)('C')<<8) + 'C'] * -1 +
        count[((int)('C')<<8) + 'G'] * 1 +
        count[((int)('C')<<8) + 'U'] * -1 +
        count[((int)('G')<<8) + 'G'] * 0 +
        count[((int)('G')<<8) + 'U'] * 1 +
        count[((int)('U')<<8) + 'U'] * -1;
#endif

    score /= num;

    return score;
}

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

