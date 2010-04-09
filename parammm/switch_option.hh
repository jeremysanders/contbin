// Interface to switches using virtual fns
// for parammm library

// This class allows parameters to be returned for switches

//      Copyright (C) 2000 Jeremy Sanders
//      Contact: jss@ast.cam.ac.uk
//               Institute of Astronomy, Madingley Road,
//               Cambridge, CB3 0HA, UK.

//      See the file COPYING for full licence details.

//      This program is free software; you can redistribute it and/or modify
//      it under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 2 of the License, or
//      (at your option) any later version.

//      This program is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.

//      You should have received a copy of the GNU General Public License
//      along with this program; if not, write to the Free Software
//      Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

#ifndef SWITCH_OPTION_PARAMMM_HH
#define SWITCH_OPTION_PARAMMM_HH

#include <string>
#include <sstream>

namespace parammm {

  // abstract base class for switch_opt
  class switch_opt
  {
  public:
    virtual ~switch_opt();

    virtual void interpret(const std::string &) const;
    virtual bool takesoption() const; // default true

    virtual switch_opt* makecopy() const = 0;

  private:
    virtual void setfromstream(std::istream *) const = 0;
  };

  // read an integer
  class pint_opt : public switch_opt {
  public:
    pint_opt(int *p) : m_p(p) {}

    pint_opt* makecopy() const;
  private:
    void setfromstream(std::istream *stream) const;
    int *m_p;
  };

  // if the switch exists, set the bool to true
  class pbool_noopt : public switch_opt {
  public:
    pbool_noopt(bool *p) : m_p(p) {}
    pbool_noopt* makecopy() const;
  private:
    void setfromstream(std::istream *stream) const;
    bool takesoption() const;  // no option
    bool *m_p;
  };

  // take a 1,0 or y,n answer
  class pbool_opt : public switch_opt {
  public:
    pbool_opt(bool *p) : m_p(p) {}
    pbool_opt* makecopy() const;
  private:
    void setfromstream(std::istream *stream) const;
    bool *m_p;
  };

  // read a c++ string
  class pstring_opt : public switch_opt {
  public:
    pstring_opt(std::string *p) : m_p(p) {}
    pstring_opt* makecopy() const;
  private:
    void setfromstream(std::istream *stream) const;
    std::string *m_p;
  };

  // read a char str, option size parameter specifies max len
  class pchar_opt : public switch_opt {
  public:
    pchar_opt(char *p, int len)  : m_p(p),m_len(len) {}
    // set len to -1 to not give a length
    pchar_opt* makecopy() const;
  private:
    void setfromstream(std::istream *stream) const;
    char *m_p;
    int m_len;
  };

  // read a floating point double
  class pdouble_opt : public switch_opt {
  public:
    pdouble_opt(double *p) : m_p(p) {}
    pdouble_opt* makecopy() const;
  private:
    void setfromstream(std::istream *stream) const;
    double *m_p;
  };

  //////////////////////////////////////////////////////
  // callbacks - callback class is run if switch exists
  // probably shouldn't be used, as deriving class
  // from switch_opt is a much better idea

  // you'll have to typecast to put your own callback
  // functions in...

  class callback {
  public:
    virtual ~callback();
  };

  typedef void (callback::* fn_callback_noopt)();
  typedef void (callback::* fn_callback_opt)(const std::string &);

  // functions that don't require switch options
  class pcallback_noopt : public switch_opt {
  public:
    pcallback_noopt(callback *cb, fn_callback_noopt fn)
      : m_cb(cb), m_fn(fn) {}
    pcallback_noopt* makecopy() const;
  private:
    void setfromstream(std::istream *stream) const;
    bool takesoption() const;  // no option
    callback *m_cb;
    fn_callback_noopt m_fn;
  };

  // functions that do require switch options
  class pcallback_opt : public switch_opt {
  public:
    pcallback_opt(callback *cb, fn_callback_opt fn)
      : m_cb(cb), m_fn(fn) {}
    pcallback_opt* makecopy() const;
  private:
    void setfromstream(std::istream *stream) const;
    callback *m_cb;
    fn_callback_opt m_fn;
  };

}  // namespace

#endif
