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

// template classes for bidirectional graph as simple container
// nodes know previous and next nodes, edges are implicit
// insert() adds node, link() creates an edge

#ifndef _GRAPH_H_
#define _GRAPH_H_

#include <vector>
#include <list>
#include <iostream>
#include <iterator>
#include <algorithm>

template<typename T> class dag_node;

template<typename T>
class dag
{
public:
    class iterator;
    class pn_iterator;

    using value_type = T;
    using pointer = T *;
    using reference = T &;

    using node_type = dag_node<T>;
    using node_container = std::list<node_type>;
    using node_ref = typename node_container::iterator;
    using const_node_ref = typename node_container::const_iterator;
    using idx_type = unsigned int;

    dag() : _nodes()
    {
         // ewww ... this is hacky! FIXME
        _nodes.push_back(node_type(T(0,'.'),-1));
        _begin=_nodes.begin();
    }

    void link(node_ref a, node_ref b);
    void link(iterator a, iterator b);

    iterator insert(T a);
    iterator begin();
    iterator end();

    pn_iterator pn_first_begin();
    pn_iterator pn_first_end();
    pn_iterator pn_last_begin();
    pn_iterator pn_last_end();


    template<typename C> void sort(C&);
    void sort(); 
    void reduce_edges();

    const T& getById(idx_type idx) const {
        auto it=_nodes.begin(), end=_nodes.end();
        while (it != end && it->_id != idx ) ++it;
#ifdef DEBUG
        if (it == end) {
            logger->critical("ARG234: {}", idx);
            exit(1);
        }
#endif
        return it->data;
    }

    unsigned int size() { return _nodes_size; };

    void print_graphviz(std::ostream& out, const char* name);

//protected:
    node_container _nodes;
    node_ref _begin;
    unsigned int _nodes_size{0};
};

template<typename T>
std::ostream&
operator<<(std::ostream& out, dag<T> d) {
    std::copy(d.begin(), d.end(), std::ostream_iterator<T>(out,","));
    return out;
}

template<typename T> struct fix_idx;

template<typename T>
class dag_node
{
public:
    using node_type = typename dag<T>::node_type;
    using node_ref = typename dag<T>::node_ref;
    using node_ref_container = std::list<node_ref>;

    dag_node(T t, unsigned int id) : data(t),_id(id) {}

    T data;

    bool operator<(const dag_node& other) const {
      // FIXME
      //return data < other.data;
      return _id < other._id;
    }
    bool operator==(const dag_node& other) const {
        return data == other.data;
    }

    unsigned int id() const {return _id;}

    node_ref_container _previous;
    node_ref_container _next;
    unsigned int _id;
};

template<typename T>
std::ostream&
operator<<(std::ostream& out, dag_node<T> dag_node) {
    out << dag_node.data;
    return out;
}

// Iterator

template<typename T>
class dag<T>::pn_iterator
{
public:
    using reference = T &;
    using value_type = T;
    using pointer = T *;
    using _container = dag<T>;
    using _self = typename dag<T>::pn_iterator;

    using node_ref_it = typename dag_node<T>::node_ref_container::iterator;

    pn_iterator(const node_ref_it& nri) : _iter(nri) {}

    bool operator!=(const pn_iterator& i) { return _iter != i._iter; }
    bool operator==(const pn_iterator& i) { return !(*this != i); }

    dag_node<T>& get_node() const { return *(*_iter); }
    _self& operator++() { ++_iter; return *this; }

    reference operator*() { return get_node().data; }
    pointer operator->() { return &(get_node().data); }

private:
    node_ref_it _iter;
};



template<typename T>
class dag<T>::iterator
{
public:
    using reference = T &;
    using value_type = T;
    using pointer = T *;
    using iterator_category = typename node_ref::iterator_category;
    using difference_type = int;


    using dag_t = dag<T>;
    friend class dag<T>;

    iterator(node_ref idx) : _idx(idx), _isNull(false) {}
    iterator(const iterator& orig) : _idx(orig._idx), _isNull(false) {}
    iterator() {}

    dag_node<T>& get_node() const { return *_idx; }
    node_ref get_node_ref() const { return _idx; }


    iterator& operator++() { ++_idx; return *this; }
    iterator& operator--() { --_idx; return *this; }

    bool operator!=(const iterator& i) { return _idx!=i._idx; }
    bool operator==(const iterator& i) { return _idx==i._idx; }
    bool operator<(const iterator& i) { return _idx<i._idx; }
    bool isNull() { return _isNull; }

    using pn_iterator = typename dag<T>::pn_iterator;

    pn_iterator prev_begin() const {
        return pn_iterator(get_node()._previous.begin());
    }
    pn_iterator prev_end() const {
        return pn_iterator(get_node()._previous.end());
    }
    pn_iterator next_begin() const {
        return pn_iterator(get_node()._next.begin());
    }
    pn_iterator next_end() const {
        return pn_iterator(get_node()._next.end());
    }

    reference operator*() { return get_node().data; }
    pointer operator->() { return &(get_node().data); }
//protected:
    node_ref _idx;
    bool     _isNull{true};
};

template<typename T>
inline typename dag<T>::pn_iterator 
prev_begin(const dag<T>&, const typename dag<T>::iterator& i)
{
  return i.prev_begin();
}

template<typename T>
inline typename dag<T>::pn_iterator 
prev_end(const dag<T>&, const typename dag<T>::iterator& i)
{
  return i.prev_end();
}

template<typename T>
inline typename dag<T>::pn_iterator 
next_begin(const dag<T>&, const typename dag<T>::iterator& i)
{
  return i.next_begin();
}

template<typename T>
inline typename dag<T>::pn_iterator 
next_end(const dag<T>&, const typename dag<T>::iterator& i)
{
  return i.next_end();
}

template<typename T>
inline const typename dag<T>::idx_type
get_node_id(const dag<T>&, typename dag<T>::iterator i) {
    return i.get_node().id();
}

template<typename T>
inline const typename dag<T>::idx_type
get_node_id(const dag<T>&, typename dag<T>::pn_iterator i) {
    return i.get_node().id();
}

template<typename S, typename T>
S
for_each_prev(T& git, S s) {
    using node_ref = typename T::dag_t::node_ref;
    using pn_iterator = typename std::list<node_ref>::iterator;
    pn_iterator it = git.get_node()._previous.begin();
    pn_iterator end = git.get_node()._previous.end();
    for(; it != end; ++it) s(*it);
    return s;
}


template<typename T>
typename dag<T>::iterator
dag<T>::begin()
{
    return iterator(++_nodes.begin());
}

template<typename T>
typename dag<T>::iterator
dag<T>::end()
{
    return iterator(_nodes.end());
}

template<typename T>
typename dag<T>::pn_iterator
dag<T>::pn_first_begin()
{
    return pn_iterator(_nodes.begin()->_next.begin());
}

template<typename T>
typename dag<T>::pn_iterator
dag<T>::pn_first_end()
{
    return pn_iterator(_nodes.begin()->_next.end());
}

template<typename T>
typename dag<T>::pn_iterator
dag<T>::pn_last_begin()
{
    return pn_iterator(_nodes.begin()->_previous.begin());
}

template<typename T>
typename dag<T>::pn_iterator
dag<T>::pn_last_end()
{
    return pn_iterator(_nodes.begin()->_previous.end());
}


// Implementations:

template<typename T>
void
dag<T>::link(node_ref a, node_ref b)
{
    a->_next.push_back(b);
    b->_previous.push_back(a);
    _nodes.begin()->_previous.remove(a);
    _nodes.begin()->_next.remove(b);
}

template<typename T>
void
dag<T>::link(iterator a, iterator b)
{
    link(a._idx,b._idx);
}

template<typename T>
typename dag<T>::iterator
dag<T>::insert(T a) {
    auto it = _nodes.insert(_nodes.end(), node_type(a, _nodes_size++));
    _nodes.begin()->_next.push_back(it);
    _nodes.begin()->_previous.push_back(it);

    return iterator(it);
}

template<typename T>
void
dag<T>::print_graphviz(std::ostream& out, const char* name)
{
    out << "digraph " << name << " { " << std::endl;
    out << "rotate=90" << std::endl;

    auto it   = _nodes.begin();
    auto iend = _nodes.end();

    for (;it != iend; ++it) {
        out << "n" << it->id() << " [ label = \""
            << it->data
            << "\" ]; ";

        {
            auto jt   = it->_next.begin();
            auto jend = it->_next.end();
            for (;  jt != it->_next.end(); ++jt) {
                out << "n" << it->id() << " -> n" << (*jt)->id() << "; ";
            }
        }
        /*
        out << " // ";
        {
            auto jt   = it->_previous.begin();
            auto jend = it->_previous.end();
            for (;  jt != it->_previous.end(); ++jt) {
                out << "n" << it->id() << " -> n" << (*jt)->id() << "; ";
            }
        }
        */
        out << std::endl;
    }
    out << "}" << std::endl;
}

// *** Sorting (reorders nodes in vector) ***

template<typename F>
struct dereference {
    dereference(F f) : _f(f){}
    dereference() : _f(){}
    using result_type = typename F::result_type;

    template<typename A, typename B>
    result_type operator()(const A& a,
                           const B& b) {
        return _f(*a, *b);
    }
    F _f;
};

// default ordering: use operator< on T
template<typename T, typename F>
struct by_data {
    by_data(F f) : _f(f){}
    bool operator()(const dag_node<T>& a, const dag_node<T>& b) {
        return _f(a.data,b.data);
    }
    F _f;
};


// maps old idx to new idx using lookup vector
template<typename node_ref>
struct lookup {
    lookup(std::vector<node_ref> &nrv) : _nrv(nrv) {}
    node_ref operator()(node_ref nr) {
        return _nrv[nr];
    }
    std::vector<node_ref>& _nrv;
};

// fixes node_references of given node using lookup vector
template<typename T>
struct fix_idx {
    using node_ref = typename dag<T>::node_ref;
    fix_idx(std::vector<node_ref> &nrv) : _nrv(nrv) {}
    void operator()(dag_node<T>& dn) {
        // fix incoming edges
        std::transform(dn._previous.begin(),dn._previous.end(),
                       dn._previous.begin(),lookup<node_ref>(_nrv));
        // fix outgoing edges
        std::transform(dn._next.begin(),dn._next.end(),
                       dn._next.begin(),lookup<node_ref>(_nrv));
        // fix self "pointer"
        dn._self = _nrv[dn._self];
    }
    std::vector<node_ref>& _nrv;
};

template<typename T> template<typename C>
void
dag<T>::sort(C& comp) {
  // sort nodes using given comparator
  _nodes.sort(by_data<T,C>(comp));
}

template<typename T> 
void 
dag<T>::sort() { 
  std::less<T> c; 
  sort(c);  
}

// reduce_edges() [ removes duplicate edges from nodes ]
template<typename T>
struct reduce_edge {
  using node_type = typename dag<T>::node_type;
  using node_ref = typename dag<T>::node_ref;
  void
  operator()(node_type& node) {
    node._previous.sort(dereference<std::less<node_type> >());
    node._previous.erase(unique(node._previous.begin(),
                                node._previous.end()),
                         node._previous.end());
    
    node._next.sort(dereference<std::less<node_type> >());
    node._next.erase(unique(node._next.begin(),
                            node._next.end()),
                     node._next.end());
  }
};

template<typename T>
void
dag<T>::reduce_edges() {
    for_each(_nodes.begin(),_nodes.end(), reduce_edge<T>());
}


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




