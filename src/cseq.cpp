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

#include "cseq.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/foreach.hpp>
#include <functional>
#include <cmath>
#include <zlib.h>

using namespace std;
using boost::lambda::bind;

using namespace sina;


cseq::cseq() 
    : name(), bases(), alignment_width(0), attributes(), score(0.f)
{
}

cseq::cseq(const char *_name, float _score, const char *_data) 
    : alignment_width(0), score(_score)
{
    name.assign(_name);
    if (_data)
        append(_data);
}

cseq::cseq(const cseq& orig) 
    : name(orig.name), bases(orig.bases),
      alignment_width(orig.alignment_width),
      attributes(orig.attributes),
      score(orig.score) 
{
}

cseq& 
cseq::operator=(const cseq& rhs) 
{
  name=rhs.name;
  bases=rhs.bases;
  alignment_width=rhs.alignment_width;
  attributes=rhs.attributes;
  score=rhs.score;
  return *this;
}
  




void
cseq::clearSequence() { 
    bases.clear(); 
    alignment_width = 0; 
}

cseq&
cseq::append(const char *str) {
    //FIXME: keep internal '.'s
    while(*str) {
        if (*str != ' ' && *str != '\n' &&  *str != '\r') {
            if (*str != '-' && *str != '.') {
                bases.push_back(aligned_base(alignment_width,*str));
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
        std::cerr << "$ cseq::append(): wrong order! "
                  << ab.getBase() << "(" << ab.getPosition() << "<" 
                  << alignment_width << ")"
                  << std::endl;

        bases.push_back(aligned_base(alignment_width, ab.getBase()));
        //FIXME shoudln't this increase the alignment_width as well?
    } 

    return *this;
}

cseq&
cseq::assign(vector<unsigned char>& dat) {
    if (dat[0]=='#')
        assignFromCompressed(&dat.front(), dat.size());
    else {
        if (dat.back() != '0')
            dat.push_back('0');
//        append((char*)&dat.front());
    #warning fix cseq::assign
    }
    return *this;
}


void
cseq::setWidth(vidx_type newWidth) {
    if (bases.empty() || newWidth >= bases.back().getPosition()) {
        // modify at will if changing only number of trailing gaps
        alignment_width = newWidth;
        return;
    }

    int bases_size=bases.size();
    int skip=1;
    
    // find the number of bases from the right where
    // <position-of-base> + <number-of-following-bases>
    // is at most <width-of-alignment>
    for(; skip<bases_size; skip++) {
        if (bases[bases_size - skip].getPosition() + skip
            <= alignment_width) {
            break;
        }
    }

    // make sure we can safely go on
    if (skip > bases_size) {
        std::clog << "cseq: cannot shrink aligment width below #bases" 
                  << std::endl;
        return;
    }

    for (; skip>0; skip--) {
        bases[bases_size - skip].setPosition(alignment_width-skip);
    }

    std::clog << "cseq: had to move last" << skip
              << "bases to shrink alignment to " << alignment_width << "." << std::endl;
}


string
cseq::getNameScore() const {
    stringstream tmp;
    tmp << name << "(" << score << ")";
    return tmp.str();
}

string
cseq::getAligned(bool nodots) const {
    string aligned;
    aligned.reserve(alignment_width);

    char dot='.';
    if (nodots) dot='-';

    vector<aligned_base>::const_iterator 
        it = bases.begin(),
        it_end = bases.end();

    unsigned int cursor = 0;
    for (; it != it_end; ++it) {
        unsigned int pos = it->getPosition();
        if (cursor > pos) {
            stringstream recent;
            int i=0;
            int left = (it_end - it) - 1;
            for (vector<aligned_base>::const_iterator jt = it + std::min(2,left); 
                 jt != bases.begin() && i<10; jt--, i++) {
                recent << " " << *jt;
            }

            std::cerr << "ERROR: broken sequence!"
                      << " C=" << cursor << " P=" << pos << " L=" << left
                      << " B=" << *it 
                      << " P=\"" << recent.str() << "\""
                      << std::endl;

        } else {
            aligned.append(pos - cursor, dot);
            cursor = pos;
        }
        dot = '-'; // (only the first "gap" is dots)
        aligned.append(1,it->getBase());
        cursor++;
        if (cursor != aligned.size()) {
            cerr << "c=" << cursor << "a=" << aligned.size()<<endl;
        }
    }

    if (cursor < alignment_width) {
        if (!nodots) dot='.';
        aligned.append(alignment_width - cursor, dot);
    }

    return aligned;
}

string
cseq::getBases() const {
    string basestr;
    basestr.reserve(bases.size());

// for_each(bases.begin(), bases.end(),
    //         basestr += bind<char>(&aligned_base::getBase,boost::lambda::_1));
    for (uint i=0; i < bases.size(); i++) {
        basestr += bases[i].getBase();
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
    typedef unsigned int uint;

    bases.push_back(aligned_base(alignment_width));
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
    compressed_data *cd = (compressed_data*)&out.front();

    cd->id='#';
    cd->size=(unsigned short)orig_size;

    compress2(cd->data, &compr_size, &buf.front(), orig_size,9);

    out.resize(compr_size+sizeof(compressed_data));
}

void
cseq::assignFromCompressed(const void* data, size_t len) {
    vector<unsigned char> buf;
    typedef unsigned int uint;
    compressed_data *cd = (compressed_data*)data;
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
        bases.push_back(aligned_base(last,buf[i]));
    }
    alignment_width = bases.back().getPosition();
    bases.pop_back();
}



char
cseq::operator[](cseq::vidx_type i) {
    vector<aligned_base>::const_iterator it = getIterator(i);
    if (it != bases.end() && i == it->getPosition())
        return it->getBase();
    else
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
                StrictWeakOrdering &cmp) {
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
    for (string::const_iterator it = in.begin(); it != in_end; ++it) {
        switch(*it) {
                case 'a':
                case 'A':
                    tmp << "\033[34m";
                    break;
                case 'g':
                case 'G':
                    tmp << "\033[35m";
                    break;
                case 'c':
                case 'C':
                    tmp << "\033[32m";
                    break;
                case 't':
                case 'T':
                case 'u':
                case 'U':
                    tmp << "\033[33m";
                    break;
                default:
                    tmp << "\033[0m";
        }
        tmp << *it;
    }
    tmp << "\033[0m";
    return tmp.str();
}

void
cseq::write_alignment(ostream& ofs, vector<cseq>& seqs,
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
    vector<cseq*>::size_type jmax = seqs.size();

    for (cseq::idx_type i = from_pos; i < to_pos; ++i) {
        bool gap = true;
        for (vector<cseq*>::size_type j = 0; j < jmax; ++j) {
            outchar[j] = seqs[j].operator[](i);
            if (outchar[j] != '-') gap=false;
        }

        if (!gap || i == to_pos-1 ) {
            for (vector<cseq*>::size_type j = 0; j < jmax; ++j) {
                out[j].append(1,outchar[j]);
            }
        }
    }

    std::map<string,list<int> > mymap;

    size_t maxlen = 0;
    for (vector<string>::iterator it = out.begin(); it != out.end(); ++it) {
        maxlen = std::max(maxlen, it->size());
        mymap[*it].push_back(it-out.begin());
    }

    ofs << "Dumping pos " << from_pos << " through " << to_pos << ":" << endl;
    for (unsigned int i = 0; i < maxlen; i+=100) {
        int len = maxlen - i;
        len = std::min(100,len);
        for (std::map<string,list<int> >::iterator it = mymap.begin(); 
             it != mymap.end(); ++it) {
            if (colors) {
                ofs << color_code(it->first.substr(i, len)) << " ";
            } else {
                ofs << it->first.substr(i, len) << " ";
            }
            bool range=false, is_last=false, is_secondlast=false;
            int last=-2;
            for (std::list<int>::iterator jt = it->second.begin();
                 jt != it->second.end(); ++jt) {
                if (range) {
                    if (*jt != last+1) {
                        ofs << last;
                        range=false;
                        ofs << " " << *jt;
                    }
                } else {
                    if (*jt == last+1) {
                        range=true;
                        ofs << "-";
                    } else
                        ofs << " " << *jt;
                }
                last = *jt;
                if (*jt +1 == (int)seqs.size()) is_last=true;
                if (*jt +2 == (int)seqs.size()) is_secondlast=true;
            }
            if (range) ofs << last;
            if (is_last) ofs << " <---(## NEW ##) ";
            if (is_secondlast) ofs << " <---(%% ORIG %%) ";

            ofs << endl;
        }
        ofs << endl;
    }
}

void
cseq::write_alignment(ostream& ofs, vector<cseq>& seqs, bool color_code) {
    cseq::idx_type imax = seqs[0].getWidth();
    cseq::write_alignment(ofs, seqs, 0, imax, color_code);
}

void
cseq::fix_duplicate_positions(std::ostream& log, bool lowercase, bool remove) {
    idx_type total_inserts = 0, longest_insert = 0, orig_inserts = 0;

    vector<aligned_base>::iterator last_it = bases.begin();
    vector<aligned_base>::iterator bases_end = bases.end();
    vector<aligned_base>::iterator curr_it = last_it+1;
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
            } else {
                // curr is the last base.
                // move iterator to end and fall through to placement
                ++curr_it;
            }
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
                vector<aligned_base>::iterator left = last_it;
                vector<aligned_base>::iterator right = curr_it;
                
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
                        std::cerr << "ERROR: no space to left and right?? sequence longe than alignment?!" << endl;
                        break;
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
            if (lowercase) last_it->setLowerCase();
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

std::list<unsigned int>
cseq::find_differing_parts(const cseq& right) const {
    typedef std::vector<aligned_base>::const_iterator bases_it;
    bases_it l_it = bases.begin(), l_end = bases.end();
    bases_it r_it = right.bases.begin(), r_end = right.bases.end();
    
    std::list<unsigned int> start_stops;
    int score = 0;
    bool bad = false;

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
            if ((char)l_it->getBase() != (char)r_it->getBase())
                score = 4;
            ++r_it;
            ++l_it;
        }
        if (l_it != l_end) 
            lpos = l_it->getPosition();
        if (r_it != r_end)
            rpos = r_it->getPosition();
        if (score > 0) {
            if (!bad) {
                int rpos = (r_it-2)->getPosition();
                start_stops.push_back(std::min(lpos,rpos));
                bad=true;
            } else {
                if (--score <= 0 && lpos == rpos) {
                    start_stops.push_back(lpos);
                    bad=false;
                }
            }
        }
    }
    if (bad==true)
        start_stops.push_back(std::min(lpos,rpos));
    
    return start_stops;
}


float
cseq::calcPairScore(const std::vector<int>& pairs) {
    // create array to hold counts for base combinations
    std::vector<int> count;
    count.resize(65536);
    for (int i=0; i<65536; ++i) count[i]=0;

    int num=0;

    for(unsigned int i=0; i<pairs.size(); ++i) {
        unsigned char left, right;
        if (pairs[i]) { // alignment position has "helix-parter"
            left = operator[](i);
            right = operator[](pairs[i]);

            if ((left != '.' and right != '.')
                and
                (left != '-' or right != '-')) {
                num++;
                if (left < right)
                    count[((int)(left)<<8) + right]++;
                else
                    count[((int)(right)<<8) + left]++;
            }
        }
    }

#ifdef DEBUG
    std::cerr << "bp_detail: ";
    for (int i=0; i<65536; ++i) {
       if (count[i]>0) {
           std::cerr << (char)(i>>8) << (char)(i&0xFF) << ": " << count[i]/2 <<"  ";
        }
    }
    std::cerr << endl;
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

struct counter {
    int agap, bgap, match, mismatch;
    counter() : agap(0), bgap(0), match(0), mismatch(0) {}
    void agap(const aligned_base&) { agap++; }
    void bgap(const aligned_base&) { bgap++; }
    void match(const aligned_base&) { match++; }
    void mismatch(const aligned_base&, const aligned_base&) { mismatch++; }
};

struct iupac_compare {
    typedef bool value_type;
    bool operator()(const aligned_base& a, const aliged_base& b) const {
        return a.comp(b);
    }
}



cseq_compare(const cseq& a, const cseq& b) {
    bases_it a_it = a.bases.begin();
    bases_it a_end = a.bases.end();
    bases_it b_it = b.bases.begin();
    bases_it b_end = b.bases.end();
    
    int a_pos = a_it->getPosition();
    int b_pos = b_it->getPosition();
    while (a_it != a_end || b_it != b_end) {
        while (a_pos > b_pos && b_it != b_end) {
            BGAP(*b_it);
            b_pos = ++b_it;
        }
        while (b_pos > a_pos && a_it != a_end) {
            AGAP(*a_it);
            a_pos = ++a_it;
        }
        while (a_pos == b_pos && a_it != a_end && b_it != b_end) {
            if (COMPARE(*a_it,*b_it)) {
                MATCH(*a_it);
            } else {
                MISMATCH(*a_it, *b_it);
            }
            a_pos = ++a_it;
            b_pos = ++b_it;
        }
    }
}
             
*/


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

