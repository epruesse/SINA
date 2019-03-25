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

// template classes for alignment matrix and alignment functions
#ifndef _MESH_H_
#define _MESH_H_

#include <vector>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <cmath>
#include <set>
#include <memory>
#include <string>


#include "align.h"
#include "graph.h"
#include "scoring_schemes.h"
#include "buffer.h"

namespace sina {

template<typename SEQ_MASTER, // mseq
         typename SEQ_SLAVE,  // cseq
         typename DATA,       // data_type
         typename ALLOCATOR = std::allocator<DATA> >
class mesh
{
public:
    using master_type = SEQ_MASTER;
    using master_iterator_type = typename master_type::iterator;
    using master_idx_type = typename master_type::idx_type;

    using slave_type = SEQ_SLAVE;
    using slave_iterator_type = typename slave_type::iterator;
    using slave_idx_type = typename slave_type::idx_type;

    using value_type = DATA;
    using vector_type = aligned_buffer<value_type>;
    using size_type = typename vector_type::size_type;

    class iterator;

    mesh() = delete;
    mesh(SEQ_MASTER &m, SEQ_SLAVE &s)
        : _master(m),
          _slave(s),
          _data(m.size() * s.size()),
          _master_begin(_master.begin()),
          _slave_begin(_slave.begin())
    {}


    template<typename IT_M, typename IT_S>
    value_type&
    operator()(const IT_M& itm, const IT_S& its) {
        return (*this)(get_node_id(_master,itm),
                       get_node_id(_slave,its));
    }

    value_type&
    operator()(const master_idx_type& mi, const slave_idx_type& si) {
        return _data[mi * _slave.size() + si];
    }

    value_type&
    operator()(size_type s) {
        return _data[s];
    }

    iterator
    begin() { 
      return iterator(*this, _master.begin(), _slave.begin()); 
    }

    iterator
    end() { 
      return iterator(*this, _master.end(), _slave.end()); 
    }

    // private:
    SEQ_MASTER _master;
    SEQ_SLAVE _slave;
    vector_type _data;

    master_iterator_type _master_begin;
    slave_iterator_type _slave_begin;

    template<typename T> friend void compute(T&);

    static size_type guessMem(SEQ_MASTER &m, SEQ_SLAVE &s) {
        return sizeof(value_type) * m.size() * s.size();
    }

#if 0
    // unused code

    /* == missing in mseq
    bool
    operator==(const mesh& rhs) {
        return ( (_data == rhs._data) && (_master == rhs._master) && (_slave == rhs._master) );
    }
    */


    void
    reuse(SEQ_MASTER &m, SEQ_SLAVE &s) {
        _master=m;
        _slave=s;
        unsigned int needed = m.size() * s.size();
        if (_data.size() < needed) _data.resize(needed);
        _master_begin=m.begin();
        _slave_begin=s.begin();
    }
#endif
};

template<typename SEQ_MASTER, // mseq, pseq
         typename SEQ_SLAVE,  // cseq
         typename DATA,
         typename ALLOCATOR>
class mesh<SEQ_MASTER,
           SEQ_SLAVE,
           DATA,
           ALLOCATOR>::iterator {
public:
    using mesh_type = mesh<SEQ_MASTER, SEQ_SLAVE, DATA, ALLOCATOR>;
    using master_type = SEQ_MASTER;
    using miterator_type = typename master_type::iterator;
    using midx_type = typename master_type::idx_type;
    using slave_type = SEQ_SLAVE;
    using siterator_type = typename slave_type::iterator;
    using sidx_type = typename slave_type::idx_type;
    using value_type = DATA;
    using size_type = typename mesh_type::size_type;

    iterator(mesh_type& _mesh,
             const miterator_type& _mit,
             const siterator_type& _sit)
        : mymesh(_mesh), mit(_mit), sit(_sit)
    {
        sidx = get_node_id(mymesh._slave, sit);
        midx = get_node_id(mymesh._master, mit);
    }

    iterator&
    operator++() {
        ++sit;
        sidx = get_node_id(mymesh._slave, sit);

        if (sit != mymesh._slave.end()) {
            return *this;
        }

        ++mit;
        midx = get_node_id(mymesh._master, mit);

        if (mit == mymesh._master.end()) {
            return *this;
        }

        sit = mymesh._slave.begin();
        sidx = get_node_id(mymesh._slave,sit);
        return *this;
    }

    bool
    operator==(const iterator& rhs) const {
        return (
//  == missing in mseq and therefor in mesh
//            (mymesh == rhs.mymesh) && ==
            (mit == rhs.mit) &&
            (sit == rhs.sit)
            );
    }

    bool operator!=(const iterator& rhs) const { return !(*this == rhs); }

    value_type&
    operator*() {
        return mymesh(*this);
    }

    operator size_type() {
        return (midx * mymesh._slave.size() + sidx);
    }

    midx_type getMidx() const { return midx; }
    sidx_type getSidx() const { return sidx; }

    mesh_type      &mymesh;
    miterator_type mit;
    midx_type      midx;
    siterator_type sit;
    sidx_type      sidx;
};


template<typename SEQ_MASTER, typename SEQ_SLAVE, typename DATA>
std::ostream&
operator<<(std::ostream& out, mesh<SEQ_MASTER,SEQ_SLAVE,DATA> m) {
    typename SEQ_MASTER::iterator mit     = m._master.begin();
    typename SEQ_MASTER::iterator mit_end = m._master.end();
    typename SEQ_SLAVE ::iterator sit_end = m._slave.end();

    out << "BEGIN MESH DATA" << std::endl;
    out << "SEQ_MASTER: " << m._master << std::endl;
    out << "SEQ_SLAVE: " << m._slave << std::endl;

    {
        typename SEQ_SLAVE::iterator sit = m._slave.begin();
        out << "         ";
        for (; sit!=sit_end; ++sit) {
            out << std::setw(5) << *sit << " ";
        }
        out << std::endl;
    }

    for (; mit!=mit_end; ++mit) {
        out << std::setw(2) << get_node_id(m._master,mit)  << ":"
            << std::setw(5) << *mit << " ";
        typename SEQ_SLAVE::iterator sit=m._slave.begin();
        for (; sit!=sit_end; ++sit) {
            out << std::setw(5) << m(mit,sit).value << " ";
        }
        out << std::endl;
    }
    out << "END MESH DATA" << std::endl;
    return out;
}

template<typename SCORING_SCHEME,
         typename SEQ_MASTER, // mseq, pseq
         typename SEQ_SLAVE>  // cseq
class transition_simple {
public:
    using master_type = SEQ_MASTER;
    using slave_type = SEQ_SLAVE;
    using value_type = typename SCORING_SCHEME::value_type;
    struct data_type;
    using mesh_type = mesh<master_type, slave_type, data_type>;
    using mesh_iterator_type = typename mesh_type::iterator;

  //private:
    const SCORING_SCHEME& s;

public:

    explicit transition_simple(const SCORING_SCHEME& _s)
        : s(_s)
    {
    }

    struct data_type {
        typename master_type::idx_type value_midx;
        typename slave_type::idx_type value_sidx;
        typename master_type::idx_type gapm_idx;
        typename slave_type::idx_type gaps_idx;

        value_type value;
        value_type gapm_val;
        value_type gaps_val;

        bool operator<(const data_type& o) const { return value<o.value; }

        void init_edge() {
            value = gapm_val = gaps_val = 1,
                value_midx = value_sidx = gapm_idx = gaps_idx = 0;
        }
        void init() {
            value = gapm_val = gaps_val = 1000000,
                value_midx = value_sidx = gapm_idx = gaps_idx = 0;
        }
    };


    template<typename BASE_TYPE_A, typename BASE_TYPE_B>
    void
    deletion(const data_type &src, data_type &dest,
             const BASE_TYPE_A& b1, const BASE_TYPE_B& b2,
             typename SEQ_MASTER::idx_type midx,
             typename SEQ_MASTER::idx_type sidx) {

        value_type value = s.deletion(src.value, b1, b2);
        value_type gap_val = s.deletion_ext(src.gapm_val, b1, b2, 0);

        if (value < gap_val) {
            dest.gapm_val = value;
            dest.gapm_idx = midx;
        } else {
            dest.gapm_val = gap_val;
            dest.gapm_idx = src.gapm_idx;
            value = gap_val;
            midx = src.gapm_idx;
        }

        if ( value < dest.value ) {
            dest.value = value;
            dest.value_midx = midx;
            dest.value_sidx = sidx;
        }
    }

    template<typename BASE_TYPE_A, typename BASE_TYPE_B>
    void
    insertion(const data_type &src, data_type &dest,
              const BASE_TYPE_A& b1, const BASE_TYPE_B& b2,
              typename SEQ_MASTER::idx_type midx,
              typename SEQ_SLAVE::idx_type sidx,
              typename SEQ_SLAVE::idx_type /*smax*/) {

        if (src.gaps_val != src.value ) {
            // opening gap
            dest.gaps_val = s.insertion(src.value, b1, b2);
            dest.gaps_idx = sidx;
        } else {
            // extending gap
            dest.gaps_val = s.insertion_ext(src.gaps_val, b1, b2,
                                            sidx-src.gaps_idx);
            dest.gaps_idx = src.gaps_idx;
        }

        if ( dest.gaps_val <= dest.value ) {
            // if gap open/extend is currently the best
            // option, remember this
            dest.value = dest.gaps_val;
            dest.value_sidx = dest.gaps_idx;
            dest.value_midx = midx;
        }
    }

    template<typename BASE_TYPE_A, typename BASE_TYPE_B>
    void
    match(const data_type &src, data_type &dest,
          const BASE_TYPE_A& b1, const BASE_TYPE_B& b2,
          typename SEQ_MASTER::idx_type midx,
          typename SEQ_SLAVE::idx_type sidx) {

        value_type value = s.match(src.value, b1, b2) ;

        if (value < dest.value) {
            dest.value = value;
            dest.value_midx = midx;
            dest.value_sidx = sidx;
        }
    }
};

template<typename SCORING_SCHEME,
         typename SEQ_MASTER_,
         typename SEQ_SLAVE_>
class transition_aspace_aware
    : public transition_simple<SCORING_SCHEME, SEQ_MASTER_, SEQ_SLAVE_>  {
public:
    typedef transition_simple<SCORING_SCHEME,
                              SEQ_MASTER_, SEQ_SLAVE_> transition;
    explicit transition_aspace_aware(const SCORING_SCHEME& _s)
        : transition(_s)
    {
    }

    struct data_type
        : public transition::data_type {
        typename transition::slave_type::idx_type gaps_max;
        void init_edge() {
            transition::data_type::init_edge();
            gaps_max = 0;
        }
        void init() {
            transition::data_type::init();
            gaps_max = 0;
        }
    };

    template<typename BASE_TYPE_A, typename BASE_TYPE_B>
    void
    insertion(const data_type &src, data_type &dest,
              const BASE_TYPE_A& b1, const BASE_TYPE_B& b2,
              typename transition::master_type::idx_type midx,
              typename transition::slave_type::idx_type sidx,
              typename transition::slave_type::idx_type smax) {

        if (smax < 1) {
            return; // can't insert
        }

        if (src.gaps_val != src.value ) {
            // opening gap
            dest.gaps_val = transition::s.insertion(src.value, b1, b2);
            dest.gaps_idx = sidx;
            dest.gaps_max = smax-1;
        } else if (src.gaps_max > 0) {
            // extending gap
            dest.gaps_val = transition::s.insertion_ext(src.gaps_val, b1, b2,
                                                        sidx-src.gaps_idx);
            dest.gaps_idx = src.gaps_idx;
            dest.gaps_max = src.gaps_max-1;
        } else {
            return;
        }

        if ( dest.gaps_val <= dest.value ) {
            // if gap open/extend is currently the best
            // option, remember this
            dest.value = dest.gaps_val;
            dest.value_sidx = dest.gaps_idx;
            dest.value_midx = midx;
        }
    }
};

/* computes all possible transitions to this node and choses the best */
template<typename TRANSITION>
struct compute_node_simple {
    using master_type = typename TRANSITION::master_type;
    using slave_type = typename TRANSITION::slave_type;
    using data_type = typename TRANSITION::data_type;
    using midx_type = typename master_type::idx_type;
    using sidx_type = typename slave_type::idx_type;
    using mpit_type = typename master_type::pn_iterator;
    using spit_type = typename slave_type::pn_iterator;

    explicit compute_node_simple(TRANSITION _t) : t(_t) {}

    template<typename T>
    void
    calc(T& it)
    {
        using mesh_type = typename T::mesh_type;
        mesh_type &mesh = it.mymesh;

        typename master_type::iterator m = it.mit;
        typename slave_type::iterator s = it.sit;

        mpit_type mpit_end = prev_end(mesh._master, m);
        spit_type spit_end = prev_end(mesh._slave, s);
        midx_type midx = it.getMidx();
        sidx_type sidx = it.getSidx();

        data_type d;
        if (prev_begin(mesh._master,m) == mpit_end || s == mesh._slave.begin()) {
            d.init_edge();
        } else {
            d.init();
        }

        for (mpit_type mpit = prev_begin(mesh._master,m); mpit != mpit_end; ++mpit) {
            midx_type mi = get_node_id(mesh._master,mpit);
            t.deletion(mesh(mi,sidx),d,*m,*s,mi,sidx);
        }

        unsigned int min_mpos = 1000000;
        for (mpit_type mpit = next_begin(mesh._master,m); mpit != next_end(mesh._master,m); ++mpit) {
          min_mpos = std::min(min_mpos, mpit->getPosition());
        }
        int max_insert = min_mpos - m->getPosition() - 1;

        for (spit_type spit = prev_begin(mesh._slave,s);
             spit != spit_end; ++spit) {
            sidx_type si = get_node_id(mesh._slave,spit);
            t.insertion(mesh(midx,si),d,*m,*s,midx,si,max_insert);
        }

        for (mpit_type mpit = prev_begin(mesh._master,m); mpit != mpit_end; ++mpit) {
            midx_type mi = get_node_id(mesh._master,mpit);
            for (spit_type spit = prev_begin(mesh._slave,s);
                 spit != spit_end; ++spit) {
                sidx_type si = get_node_id(mesh._slave,spit);

                t.match(mesh(mi,si),d,*m,*s,mi,si);
            }
        }
        *it = d;
    }
    TRANSITION t;
};


/* computes values in a mesh by executing NODE_COMPUTOR for every
   position */
template<typename MESH_TYPE,
         typename NODE_COMPUTOR>
void
compute(MESH_TYPE& mesh, NODE_COMPUTOR& nc) {
    using master_iterator_type = typename MESH_TYPE::master_iterator_type;
    using slave_iterator_type = typename MESH_TYPE::slave_iterator_type;

    mesh._master.sort();
    mesh._slave.sort();

    master_iterator_type m = mesh._master.begin();
    master_iterator_type mend   = mesh._master.end();
    slave_iterator_type  send   = mesh._slave.end();

    typename MESH_TYPE::iterator it = mesh.begin(), end = mesh.end();

    for (;it != end ; ++it) {
        nc.calc(it);
    }
}


/** backtrack creates the cseq result sequence from a mesh that has
 * been computed. returns float score.
 */
template<typename MESH_TYPE, typename TRANSITION>
float
backtrack(MESH_TYPE& mesh, cseq& out, TRANSITION &tr,
          OVERHANG_TYPE overhang_pos, LOWERCASE_TYPE lowercase,
          INSERTION_TYPE insertion,
          int& cutoff_head, int& cutoff_tail,
          std::ostream& log) {
    using master_type = typename MESH_TYPE::master_type;
    using slave_type = typename MESH_TYPE::slave_type;
    using master_idx_type = typename MESH_TYPE::master_idx_type;
    using slave_idx_type = typename MESH_TYPE::slave_idx_type;
    using m_pnit_type = typename master_type::pn_iterator;
    using s_pnit_type = typename slave_type::pn_iterator;
    using master_base_type = typename master_type::value_type;
    using slave_base_type = typename slave_type::value_type;

    master_idx_type alig_width = mesh._master.getWidth();

    // find outer nodes in slave sequence (hack... :/ )
    slave_idx_type sbegin = get_node_id(mesh._slave, mesh._slave.begin());
    slave_idx_type send   = get_node_id(mesh._slave, --mesh._slave.end());

    // find outer nodes in reference graph
    std::set<master_idx_type> set_mbegin, set_mend;
    for (m_pnit_type pit = mesh._master.pn_first_begin();
         pit != mesh._master.pn_first_end(); ++pit) {
        set_mbegin.insert(get_node_id(mesh._master, pit));
    }
    for (m_pnit_type pit = mesh._master.pn_last_begin();
         pit != mesh._master.pn_last_end(); ++pit) {
        set_mend.insert(get_node_id(mesh._master, pit));
    }

    // find starting point for alignment
    // check right side
    master_idx_type m = get_node_id(mesh._master,
                                    mesh._master.pn_last_begin());
    // match last slave base with each graph node
    for (typename master_type::iterator pit = mesh._master.begin();
         pit != mesh._master.end(); ++pit) {
        master_idx_type tmp = get_node_id(mesh._master, pit);
        if (mesh(tmp,send)<mesh(m,send)) {
            m=tmp;
        }
    }
    slave_idx_type s = send;
    // for each graph end-node
    for (m_pnit_type pit = mesh._master.pn_last_begin();
         pit != mesh._master.pn_last_end(); ++pit) {
        master_idx_type mtmp = get_node_id(mesh._master, pit);
        // for each
        for (s_pnit_type sit = mesh._slave.begin();
             sit != mesh._slave.end(); ++sit) {
            slave_idx_type stmp = get_node_id(mesh._slave, sit);
            if (mesh(mtmp, stmp)<mesh(m,s)) {
                m=mtmp, s=stmp;
            }
        }
    }

    // handle right hand overhang
    cutoff_tail = send - s;
    if (cutoff_tail && overhang_pos != OVERHANG_REMOVE) {
        int pos;
        if (overhang_pos == OVERHANG_ATTACH) {
            pos = alig_width -1 - mesh._master.getById(m).getPosition() - cutoff_tail;
        } else { // OVERHANG_EDGE
            pos = 0;
        }
        auto it = mesh._slave.rbegin();
        auto end = it + cutoff_tail;

        if (lowercase == LOWERCASE_UNALIGNED) {
            for (; it != end; ++it) {
                out.append(aligned_base(std::max(0, pos++), it->getBase().setLowerCase()));
            }
        } else {
            for (; it != end; ++it) {
                out.append(aligned_base(std::max(0, pos++), it->getBase()));
            }
        }
    }

    // calculate score
    auto rval = (float)mesh(m,s).value;

    unsigned int pos=alig_width -1 -mesh._master.getById(m).getPosition();
    float sum_weight=0;
    int aligned_bases=0;

    // add pair slave base with position from reference base and
    // append the newly aligned base to result vector.
    slave_base_type ab1(pos,mesh._slave.getById(s).getBase());
    out.append(ab1);
    aligned_bases++;

    // count used weights 
    {
      // create a copy of the master base object
      master_base_type ab2(mesh._master.getById(m));
      // set the copy to be a match to the slave
      ab2.setBase(mesh._slave.getById(s).getBase());
      // count score "had there been a match"
      sum_weight = tr.s.match(sum_weight,ab2,ab1);
    }

    // while we have neither reached the leftmost base of the slave
    // nor one of the leftmost bases of the reference
    while (s!=sbegin && set_mbegin.find(m) == set_mbegin.end()) {
        // get the slave/reference base pair that has
        // yielded the best score to get here
        slave_idx_type snew = mesh(m,s).value_sidx;
        m = mesh(m,s).value_midx;

        // if we had a deletion in the slave, the next cell
        // will point to the same slave idx. we want the second
        // cell, as that is the one that was reached via a "match"
        // => if moving purely along master (and not at end of slave)
        //    skip first cell (the one indicating the deletion)
        if (snew == mesh(m,snew).value_sidx && snew != 0) {
          m = mesh(m,snew).value_midx;
        }

        // get the position in the reference alignment
        pos = alig_width - 1 - mesh._master.getById(m).getPosition();

        // step left in the slave sequence until we reach the base
        // the score was based on (only once for match, multiple
        // iterations of the loop only if slave sequence has had
        // insertions relative to reference)
        while (s != snew) {
            --s;

#ifdef DEBUG
            if (s != snew) {
              log << "multi insertion: " << alig_width - pos << "; ";
            }
#endif

            // create aligned base from position/base pairing           
            slave_base_type ab1(pos, mesh._slave.getById(s).getBase());

            // remember
            out.append(ab1);
            aligned_bases++;

            // count used weights (see above for explanation)
            master_base_type ab2(mesh._master.getById(m));
            ab2.setBase(mesh._slave.getById(s).getBase());
            sum_weight = tr.s.match(sum_weight,ab2,ab1);
        }
    }

    // if and while we have not yet reached the leftmost base of the slave
    // sequence, that is if we have overhang of the slave sequence on the
    // left side
    if (s != sbegin) {
        cutoff_head = s - sbegin;

        switch(overhang_pos) {
        case OVERHANG_ATTACH:
            while (s-- != sbegin) {
                char b = mesh._slave.getById(s).getBase();
                aligned_base ab(std::min(alig_width-1,++pos), b);
                if (lowercase == LOWERCASE_UNALIGNED) {
                    ab.setLowerCase();
                }
                out.append(ab);
            }
            break;
        case OVERHANG_REMOVE:
            // do nothing
            break;
        case OVERHANG_EDGE:
            int n = s-sbegin;
            while (n--) {
                char b = mesh._slave.getById(n).getBase();
                aligned_base ab(alig_width - n - 1,b);
                if (lowercase == LOWERCASE_UNALIGNED) {
                    ab.setLowerCase();
                }
                out.append(ab);
            }
            break;
        }
    } else {
        cutoff_head = 0;
    }

    out.setWidth(alig_width);
    out.reverse();
    out.fix_duplicate_positions(log, lowercase==LOWERCASE_UNALIGNED,
                                insertion==INSERTION_REMOVE);

    if (out.getWidth() > alig_width) {
        const char* wrn = "warning: result sequence too wide!";
        log << wrn;
    }

    log << "scoring: raw=" << rval << ", weight=" << sum_weight
        << ", query-len=" << mesh._slave.size()
        << ", aligned-bases=" << aligned_bases
        << ", score=" << rval/sum_weight <<"; ";

    return rval/sum_weight;
}

} // namespace sina

#endif // _MESH_H_

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

