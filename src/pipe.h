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

#ifndef _PIPE_H_
#define _PIPE_H_

#include <list>
#include <string>

#include <boost/foreach.hpp>

#include "object_pipe.h"

//#define DEBUG_PIPE
#ifdef DEBUG_PIPE
#define DBG(x) std::cerr << x << std::endl
#else
#define DBG(x)
#endif

template<typename IN, typename OUT> class PipeSerialSegment;

struct PipeEOF {};

class basic_PipeElement
{
protected:
public:
    basic_PipeElement() {}

    virtual void run(void) = 0;
    virtual ~basic_PipeElement() {}
    virtual std::string getName() const { return "unkown"; }
    virtual std::list<const basic_PipeElement*> getPEList() {
        std::list<const basic_PipeElement*> l;
        l.push_back(this);
        return l;
    }
};



/////////////// typed_PipeElement ////////////////

template<typename IN>
class input_PipeElement {
public:
    input_PipeElement(ObjectIStream<IN> i) : in(i) {}
    input_PipeElement() : in() {}
    ObjectIStream<IN> in;
    virtual void setIn(ObjectIStream<IN> i) {
        DBG("setIn 1 {");
        in = i;
        DBG("} // set IN 1");
    }

    virtual ~input_PipeElement() {}
};

template<typename OUT>
class output_PipeElement {
public:
    output_PipeElement(ObjectOStream<OUT> o) : out(o) {}
    output_PipeElement() : out() {}
    ObjectOStream<OUT> out;
    virtual void setOut(ObjectOStream<OUT> o) { out = o; }
    virtual ~output_PipeElement() {}
};

template <typename IN, typename OUT>
class typed_PipeElement
    : public basic_PipeElement,
      public input_PipeElement<IN>,
      public output_PipeElement<OUT>
{
public:
    typed_PipeElement(const typed_PipeElement<IN, OUT>& tpe)
        :  input_PipeElement<IN>(tpe.in),
           output_PipeElement<OUT>(tpe.out)
    {}
    typed_PipeElement()
        : input_PipeElement<IN>(),
          output_PipeElement<OUT>()
    {}

    template<typename NEW_OUT>
    PipeSerialSegment<IN, NEW_OUT>
    operator|(typed_PipeElement<OUT, NEW_OUT>& rhs) {
        DBG("tpe<IN,OUT>.op|:" << getName() << "|" << rhs.getName());
        PipeSerialSegment<IN,OUT> pss(*this);
        return pss | rhs;
    }
};

template <typename IN>
class typed_PipeElement<IN, void>
    : public basic_PipeElement,
      public input_PipeElement<IN>
{
public:
    typed_PipeElement(const typed_PipeElement<IN, void>& tpe)
        : input_PipeElement<IN>(tpe.in) {}
    typed_PipeElement() {}

};


template <typename OUT>
class typed_PipeElement<void, OUT>
    : public basic_PipeElement,
      public output_PipeElement<OUT>
{
public:
    typed_PipeElement(const typed_PipeElement<void, OUT>& tpe)
        : output_PipeElement<OUT>(tpe.out) {}
    typed_PipeElement() {}

    template<typename NEW_OUT>
    PipeSerialSegment<void, NEW_OUT>
    operator|(typed_PipeElement<OUT, NEW_OUT>& rhs) {
        DBG("tpe<void,OUT>.op|:" << getName() << "|" << rhs.getName());
        PipeSerialSegment<void,OUT> pss(*this);
        return pss | rhs;
    }
};

template <>
class typed_PipeElement<void,void>
    : public basic_PipeElement
{
public:
    typed_PipeElement() {}
    typed_PipeElement(const typed_PipeElement<void, void>& /*tpe*/) {}
};



//////////////// PipeElement ////////////

// pass-through pipe (middle piece)
template<typename IN, typename OUT>
class PipeElement
    : public typed_PipeElement<IN, OUT>
{
public:
    typedef const PipeElement<IN, OUT>& const_reference;


    PipeElement(const typed_PipeElement<IN, OUT>& tpe)
        : typed_PipeElement<IN, OUT>(tpe) {}
    PipeElement()
        : typed_PipeElement<IN,OUT>() {}

    virtual OUT operator()(IN) = 0;

    void run(void) {
        IN i;
        if (!this->in.pop(i)) {
            this->out.close();
            throw PipeEOF();
        }

        OUT o = this->operator()(i);
        this->out.push(o);
    }

};

// output pipe (source)
template<typename OUT>
class PipeElement<void, OUT>
    : public typed_PipeElement<void, OUT>
{
    typedef typed_PipeElement<void, OUT> base_t;

public:
    typedef const PipeElement<void, OUT>& const_reference;

    PipeElement() {}
    template<typename S>
    PipeElement(const PipeElement<void, OUT>& pe)
        : base_t(pe) {}
    PipeElement(const base_t& tpe)
        : base_t(tpe) {}

    virtual OUT operator()(void) = 0;


    void run(void) {
        try {
            OUT o = this->operator()();
            this->out.push(o);
        } catch (PipeEOF &peof) {
            this->out.close();
            throw peof;
        }
    }
};

// input pipe (sink)
template<typename IN>
class PipeElement<IN, void>
    : public typed_PipeElement<IN, void>
{
    typedef typed_PipeElement<IN, void> base_t;
public:
    typedef const PipeElement<IN, void>& const_reference;

    PipeElement() {}
    template<typename S>
    PipeElement(const PipeElement<IN, void>& pe)
        : base_t(pe)  {}
    PipeElement(const typed_PipeElement<IN, void>& tpe)
        : base_t(tpe) {}

    virtual void operator()(IN) = 0;

    void run(void) {
        IN i;
        if (!this->in.pop(i)) throw PipeEOF();
        this->operator()(i);
    }
};

// completed pipe
template <>
class PipeElement<void,void>
    :  public typed_PipeElement<void, void>
{
    typedef typed_PipeElement<void, void> base_t;
public:
    PipeElement() {}
    template<typename S>
    PipeElement(const PipeElement<void, void>& pe)
        : base_t(pe) {}
    PipeElement(const typed_PipeElement<void, void>& tpe)
        : base_t(tpe) {}

    virtual void operator()(void) = 0;
    typedef const PipeElement<void, void>& const_reference;
    void run() {
        this->operator()();
    }
};



//////////////// PipeSerialSegment ///////////////////

typedef PipeSerialSegment<void,void> Pipe;

template<typename IN, typename OUT>
class basic_PipeSerialSegment
    : public PipeElement<IN, OUT>
{
    typedef PipeElement<IN, OUT> base;
public:
    // not really public, but there is no "access only from template"
    std::list<const basic_PipeElement*> pipeElems;

    basic_PipeSerialSegment(const typed_PipeElement<IN, OUT>& pe)
        : base(pe)
    {
        pipeElems.push_back(&pe);
    }

    basic_PipeSerialSegment(const basic_PipeSerialSegment<IN, OUT>& pss)
        : base(pss)
    {
        std::copy(pss.pipeElems.begin(), pss.pipeElems.end(),
                  std::back_inserter(pipeElems));
    }

    basic_PipeSerialSegment() {}

    void destroy() {
      for (std::list<const basic_PipeElement*>::iterator it = pipeElems.begin();
           it != pipeElems.end(); ++it) {
           delete *it;
      }
    }

    template<typename NEW_OUT>
    PipeSerialSegment<IN, NEW_OUT>
    operator|(typed_PipeElement<OUT, NEW_OUT>& rhs) {
        DBG("construct PSS from " << getName() << "|" << rhs.getName() << "{");
        PipeSerialSegment<IN, NEW_OUT> pss;
        std::copy(pipeElems.begin(), pipeElems.end(),
                  std::back_inserter(pss.pipeElems));
        std::list<const basic_PipeElement*> pel = rhs.getPEList();
        rhs.setIn(base::out.newIStream());
        std::copy(pel.begin(), pel.end(),
                  std::back_inserter(pss.pipeElems));
        pss.out = rhs.out;
        DBG("} // construct PSS from...");
        return pss;
    }

    PipeSerialSegment<IN, void>
    operator|(typed_PipeElement<OUT, void>& rhs) {
        DBG("construct(,void) PSS from " << getName() << "|" << rhs.getName() << "{");
        PipeSerialSegment<IN, void> pss;
        std::copy(pipeElems.begin(), pipeElems.end(),
                  std::back_inserter(pss.pipeElems));
        std::list<const basic_PipeElement*> pel = rhs.getPEList();
        rhs.setIn(base::out.newIStream());
        std::copy(pel.begin(), pel.end(),
                  std::back_inserter(pss.pipeElems));
        //pss.out = rhs.out;
        DBG("} // construct(,void) PSS from...");
        return pss;
    }

    std::string getName() const {
        std::string s;
        BOOST_FOREACH(const basic_PipeElement* pe, pipeElems) {
            s += " | " + pe->getName();
            //std::cerr << pe << std::endl;
        }
        return "PipeSerialSegment(" + s + ")";
    }

    std::list<const basic_PipeElement*> getPEList() {
        return pipeElems;
    }

};

template<typename IN, typename OUT>
class PipeSerialSegment
    : public basic_PipeSerialSegment<IN, OUT>  {
public:
    PipeSerialSegment(const typed_PipeElement<IN, OUT>& pe)
        : basic_PipeSerialSegment<IN,OUT>(pe) {}
    PipeSerialSegment() {}

    void setIn(ObjectIStream<IN> i) {
        DBG("setIn 2 " << this-> getName() << "{");
        this->in = i;
        basic_PipeElement *bpe = const_cast<basic_PipeElement*>(*this->pipeElems.begin());
        dynamic_cast<input_PipeElement<IN>* >(bpe)->setIn(i);
        DBG("} // setIN 2");
    }
    virtual void setOut(ObjectOStream<OUT> o) {
        this->out = o;
        basic_PipeElement *bpe = const_cast<basic_PipeElement*>(this->pipeElems.back());
        dynamic_cast<output_PipeElement<OUT>* >(bpe)->setOut(o);
    }

    OUT operator()(IN) {
        std::cerr << "error in line " << __LINE__ << std::endl;
        return OUT();
    }
};

template<typename IN>
class PipeSerialSegment<IN,void>
    : public basic_PipeSerialSegment<IN,void>  {
public:
    PipeSerialSegment(const typed_PipeElement<IN, void>& pe)
        : basic_PipeSerialSegment<IN, void>(pe) {}
    PipeSerialSegment() {}

    void operator()(IN) {
        std::cerr << "error in line " << __LINE__ << std::endl;
    }
};

template<typename OUT>
class PipeSerialSegment<void,OUT>
    : public basic_PipeSerialSegment<void,OUT>  {
public:
    PipeSerialSegment(const typed_PipeElement<void, OUT>& pe)
        : basic_PipeSerialSegment<void, OUT>(pe) {}
    PipeSerialSegment() {}
    OUT operator()() {
        std::cerr << "error in line " << __LINE__ << std::endl;
        return OUT();
    }
};

template<>
class PipeSerialSegment<void,void>
    : public basic_PipeSerialSegment<void,void>  {
public:
    PipeSerialSegment(const typed_PipeElement<void, void>& pe)
        : basic_PipeSerialSegment<void, void>(pe) {}
    PipeSerialSegment(const PipeSerialSegment<void,void>& pss)
        : basic_PipeSerialSegment<void, void>(pss) {}
    PipeSerialSegment() {}
    void operator()() {
        std::list<const basic_PipeElement*> myElems;
        std::copy(pipeElems.begin(), pipeElems.end(),
                  std::back_inserter(myElems));
        while(!myElems.empty()) {
            std::list<const basic_PipeElement*>::iterator i, end = myElems.end();
            for (i=myElems.begin(); i != end; ++i) {
                try {
                    basic_PipeElement* j = const_cast<basic_PipeElement*>(*i);
                    DBG("exec " << this->getName() << " " << j->getName());
                    j->run();
                } catch (PipeEOF &peof) {
                    DBG("remove " << *i);
                    myElems.erase(i--);
                }
            }
        }
    }
};

/*
class PipeManager {
    int m_threads;
    int m_queue_size;
public:
    PipeManager(int threads, int queue_size);

    void execute(Pipe& p);
};
*/

#endif //_PIPE_H_

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
