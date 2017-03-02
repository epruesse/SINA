/*
Copyright (c) 2006-2017 Elmar Pruesse <elmar.pruesse@ucdenver.edu>

This file is part of SINA.

SINA is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/* copyright Elmar Pruesse (epruesse@mpi-bremen.de) 2006-2011
 *
 * This file contains a series of classes implementing thread safe
 * copyable reader and writer objects to an implicitly created pipe.
 * Once all writers have been close()d or destroyed, all blocking
 * readers will be awoken.
 *
 * Usage example:
 * void producer(ObjectOStream<int> out) {
 *   for (int i=0; i<200; i++) {
 *     out.push(i); // may block if queue full!
 *   }
 * }
 *
 * void consumer(ObjectIStream<int> in) {
 *   int i;
 *   while (in.pop(i)) {
 *      cout << "got an " << i << endl;
 *   }
 * }
 *
 * void start_threads() {
 *   ObjectOStream<int> out;
 *
 *   thread consumer_thread(bind(&consumer, out.newIstream()));
 *   new thread(bind(&producer, out.newOStream()));
 *   consumer_thread.wait();
 * }
 */

#ifndef _OBJECT_PIPE_H_
#define _OBJECT_PIPE_H_

#include <execinfo.h>

#include <boost/thread/mutex.hpp>
#include <boost/thread/condition.hpp>
#include <boost/shared_ptr.hpp>

#include <queue>
#include <iostream>
#include <exception>
#include <string>

//#define DEBUG_LOCAL 1

template<typename T> class ObjectStream;
template<typename T> class ObjectOStream;
template<typename T> class ObjectIStream;

/**
 *  baisc_ObjectPipe
 *
 *  base class for ObjectPipie containing type-independent members
 *  (debug, reader/writer counting, mutex, condition, size)
 */

class basic_ObjectPipe {
protected: // all of class is protected
  basic_ObjectPipe(unsigned int size, std::string name)
        : max_size(size), size(0),  not_empty(), not_full(), mutex(),
          m_good(true), m_readers(0), m_writers(0),
          m_name(name)
    {
        debug("new pipe");
    }

    ~basic_ObjectPipe() {
         debug("delete pipe");
    }

  /**
   * increases writer count
   */
  void add_writer() {
    boost::mutex::scoped_lock lock(mutex);
    m_writers++;
    debug("add_writer", m_writers);
  }

  /**
   * decreases writer count
   */
  void remove_writer() {
    boost::mutex::scoped_lock lock(mutex);
    debug("remove_writer", m_writers);
    if (--m_writers>0) return;
    debug("remove_writer, removed last writer, waking readers");
    m_good=false;
    not_empty.notify_all();
  }

  /**
   * inreases reader count
   */
  void add_reader() {
    boost::mutex::scoped_lock lock(mutex);
    m_readers++;
    debug("add_reader", m_readers);
  }

  /**
   * decreases reader count
   */
  void remove_reader() {
    boost::mutex::scoped_lock lock(mutex);
    m_readers--;
    debug("remove_reader", m_readers);
  }

#ifdef DEBUG_LOCAL
  void debug(std::string s, int t) {
    if (!m_name.empty()) {
      std::cerr << m_name << ": " << s << " " << t << std::endl;
    } else {
      std::cerr << this << ": " << s << " " << t << std::endl;
    }
  }
  void debug(std::string s) {
    if (!m_name.empty()) {
      std::cerr << m_name << ": " << s << " " << std::endl;
    } else {
      std::cerr << this << ": " << s << " " << std::endl;
    }
  }
#else
  void debug(std::string, int) {}
  void debug(std::string) {}
#endif

  unsigned int max_size, size;
  boost::condition not_empty, not_full;
  boost::mutex mutex;
  bool m_good;
  int m_readers;
  int m_writers;
  std::string m_name;
};

/**
 * ObjectPipe<T>
 *
 * This class contains the type specific stuff. The STL queue and the
 * push/pop functions to access it safely. Access is not allowed directly.
 * Instead, use the factory methods newIStream and newOStream to create
 * ObjectIStream and ObjectOStream objects.
 */
template<typename T>
class ObjectPipe : private boost::noncopyable, public basic_ObjectPipe {
public:
    typedef T value_type;
    typedef T& reference;
    typedef T* pointer;

    typedef ObjectOStream<value_type> OStream;
    typedef ObjectIStream<value_type> IStream;

    ObjectPipe(unsigned int size, std::string name="")
        : basic_ObjectPipe(size,name)
    {
    }

    ObjectOStream<T> newOStream() {
        return ObjectOStream<T>(this);
    }

    ObjectIStream<T> newIStream() {
        return ObjectIStream<T>(this);
    }

private:
    void push(const T& data) {
        boost::mutex::scoped_lock lock(mutex);
        if (!m_good) {
            debug("push on closed pipe");
            return;
        }
        if (size == max_size) {
            debug("push, pipe full, waiting");
            not_full.wait(lock);
            debug("push unlocked");
        }
        q.push(data);
        size++;
        debug("push, size",size);
        not_empty.notify_one();
    }

    bool try_pop(T& data) {
        return pop(data,false);
    }

    bool pop(T& data, bool block=true) {
        boost::mutex::scoped_lock lock(mutex);
        while(q.empty()) {
            if (!m_good) {
                debug("pop, queue closed and empty");
                return false;
            }
            if (block) {
                debug("pop, pipe empty, waiting");
                not_empty.wait(lock);
            }
        }
        data = q.front();
        q.pop();
        --size;
        debug("pop, size", size);
        not_full.notify_one();
        return true;
    }

    friend class ObjectStream<T>;
    friend class ObjectIStream<T>;
    friend class ObjectOStream<T>;

    void print() {
        debug("print called");
    }

    std::queue<T> q;
};

template<>
class ObjectPipe<void> {
};

/**
 * ObjectStream
 *
 * Baseclass for ObjectIStream and ObjectOStream. 
 */
template<typename T>
class ObjectStream {
public:
    bool good() { return m_pipe != 0 && m_pipe->m_good; }

    ObjectOStream<T> newOStream() {
        return ObjectOStream<T>(m_pipe);
    }

    ObjectIStream<T> newIStream() {
        return ObjectIStream<T>(m_pipe);
    }

    void print() { m_pipe->print(); }

protected:
    typedef T  value_type;
    typedef T& reference;
    typedef T* pointer;

    boost::shared_ptr<ObjectPipe<T> > m_pipe;

    ObjectStream(boost::shared_ptr<ObjectPipe<T> > p)
        : m_pipe(p)
    {}
    ObjectStream(const ObjectStream<T>& other)
        : m_pipe(other.m_pipe)
    {}


protected:
private:
    ObjectStream& operator=(ObjectStream&);
};

/**
 * ObjectOStream<T>
 *
 * Objects of this class can be used to push into the ObjectPipe.
 * If all ObjectOStream objects on an ObjectPipe have been closed
 * or destroyed, ObjectIStreams blocking in pop will be woken
 * with FALSE as return value.
 * Class contains factory methods to create more ObjectIStreams
 * and ObjectOStreams.
 */
template<typename T>
class ObjectOStream : public ObjectStream<T> {
public:
  typedef ObjectStream<T> base;
  typedef T  value_type;
  typedef T& reference;
  typedef T* pointer;

  /* create ObjectPipe with custom size implicitly */
  ObjectOStream(int size)
    : base(boost::shared_ptr<ObjectPipe<T> >(new ObjectPipe<T>(size)))
  {
    base::m_pipe->add_writer();
  }

  /* create ObjectPipe with size=200 implicitly */
  ObjectOStream()
    : base(boost::shared_ptr<ObjectPipe<T> >(new ObjectPipe<T>(200)))
  {
    base::m_pipe->add_writer();
  }

  /* take ObjectPipe-shared_ptr as argument */
  ObjectOStream(boost::shared_ptr<ObjectPipe<T> > p)
    : base(p)
  {
    p->add_writer();
  }

  /* get ObjectPipe from other ObjectOStream */
  ObjectOStream(const ObjectOStream& other)
    : base(other)
  {
    if (base::m_pipe)
      base::m_pipe->add_writer();
    base::m_pipe->debug("copy construct ostream");
  }

  /* destruct implicitly closes */
  ~ObjectOStream() {
    close();
  }

  /* assign changes our objectPipe */
  ObjectOStream& operator=(const ObjectOStream& other) {
    // fixme: shouldn't we close first?
    base::m_pipe = other.base::m_pipe;
    if (base::m_pipe)
      base::m_pipe->add_writer();
    base::m_pipe->debug("assign ostream");
    return *this;
  }

  /* close = deregister with object pipe */
  void close() {
    if (base::m_pipe) {
      base::m_pipe->remove_writer();
      base::m_pipe.reset();
    }
  }

  /* push is handled by pipe. */
  void push(reference t) {
    if (base::m_pipe)     // fixme: only check this if debug?
      base::m_pipe->push(t);
    else
      std::cerr << this << " pushing stuff nowhere?" << std::endl;
    }

  /*
    ObjectOStream& operator<<(const reference t) {
    base::m_pipe.push(t);
    return *this;
    }*/
private:
};

/**
 * ObjectIStream<T>
 *
 * Objects of this class allow read access on their ObjectStream<T>.
 * bool pop(T&) will block until either an object has been pushed or
 * the last reader on the ObjectPipe has been closed. If no more
 * readers are live, pop will return false.
 */
template<typename T>
class ObjectIStream : public ObjectStream<T> {
public:
  typedef ObjectStream<T> base;
  typedef T  value_type;
  typedef T& reference;
  typedef T* pointer;

  ObjectIStream(int size)
    : base(boost::shared_ptr<ObjectPipe<T> >(new ObjectPipe<T>(size)))
  {
    base::m_pipe->add_reader();
  }

  ObjectIStream()
    : base(boost::shared_ptr<ObjectPipe<T> >(new ObjectPipe<T>(200)))
  {
    base::m_pipe->add_reader();
  }

  ObjectIStream(boost::shared_ptr<ObjectPipe<T> > p) : base(p) {
    p->add_reader();
  }

  ObjectIStream(const ObjectIStream& other)
    : base(other) {
    if (base::m_pipe)
      base::m_pipe->add_reader();
  }

  ObjectIStream& operator=(const ObjectIStream& other){
    base::m_pipe = other.base::m_pipe;
    if (base::m_pipe)
      base::m_pipe->add_reader();

    return *this;
  }

  void close() {
    if (base::m_pipe) {
      base::m_pipe->remove_reader();
      base::m_pipe.reset();
    }
  }

  ~ObjectIStream() {
    close();
  }

  bool try_pop(reference t) {
    if (base::m_pipe)
      return base::m_pipe->try_pop(t);
    return false;
  }

  bool pop(reference t) {
    if (base::m_pipe)
      return base::m_pipe->pop(t);
    std::cerr << "popping stuff from nowhere?" << std::endl;
    return false;
  }

  /*
    ObjectIStream& operator>>(reference t) {
    if (!base::m_pipe.pop(t)) {
    return *this;
    }
  */
};

#endif // _OBJECT_PIPE_H_

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
