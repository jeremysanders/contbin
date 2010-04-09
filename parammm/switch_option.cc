// implementation of switch options for parammm

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

#include <cstring>
#include <iostream>
#include "switch_option.hh"
#include "exceptions.hh"

using std::string;

namespace parammm
{

  ////////////////////////////////////////////////////
  // Abstract base class

  switch_opt::~switch_opt()
  {
  }

  void switch_opt::interpret(const std::string &str) const
  {
    std::istringstream stream(str);
    setfromstream(&stream);
    if(!stream)
      throw except_inv_switch_option(str);
  }

  bool switch_opt::takesoption() const
  {
    return true;
  }

  ///////////////////////////////////////////
  // Now type definitions

  void pint_opt::setfromstream(std::istream *str) const
  {
    *str >> *m_p;
  }

  pint_opt* pint_opt::makecopy() const
  {
    return new pint_opt(m_p);
  }

  void pbool_noopt::setfromstream(std::istream *str) const
  {
    *m_p = true;
  }

  pbool_noopt* pbool_noopt::makecopy() const
  {
    return new pbool_noopt(m_p);
  }

  bool pbool_noopt::takesoption() const
  {
    return false;
  }

  void pbool_opt::setfromstream(std::istream *str) const
  {
    char c;
    *str >> c;
    switch(c) {
    case 'y':
    case 'Y':
    case '1':
      *m_p=true; break;
    case 'n':
    case 'N':
    case '0':
      *m_p=false; break;
    default:
      str->setstate(std::ios::failbit);
    }
  }

  pbool_opt* pbool_opt::makecopy() const
  {
    return new pbool_opt(m_p);
  }

  void pstring_opt::setfromstream(std::istream *str) const
  {
    getline(*str, *m_p, char(1));
  }

  pstring_opt* pstring_opt::makecopy() const
  {
    return new pstring_opt(m_p);
  }

  void pchar_opt::setfromstream(std::istream *str) const
  {
    string s;
    getline(*str, s, char(1));
    if(m_len<0)
      strcpy(m_p, s.c_str() );
    else {
      strncpy(m_p, s.c_str(), m_len-1);
      m_p[m_len-1] = 0;
    }
  }

  pchar_opt* pchar_opt::makecopy() const
  {
    return new pchar_opt(m_p, m_len);
  }

  void pdouble_opt::setfromstream(std::istream *str) const
  {
    *str >> *m_p;
  }

  pdouble_opt* pdouble_opt::makecopy() const
  {
    return new pdouble_opt(m_p);
  }

  void pcallback_noopt::setfromstream(std::istream *str) const
  {
    (m_cb->*m_fn)();  // call obj->fn()
  }

  pcallback_noopt* pcallback_noopt::makecopy() const
  {
    return new pcallback_noopt(m_cb, m_fn);
  }

  bool pcallback_noopt::takesoption() const
  {
    return false;
  }

  void pcallback_opt::setfromstream(std::istream *str) const
  {
    string param;
    getline(*str, param, char(1));
    (m_cb->*m_fn)(param);  // call obj->fn()
  }

  pcallback_opt* pcallback_opt::makecopy() const
  {
    return new pcallback_opt(m_cb, m_fn);
  }

}
