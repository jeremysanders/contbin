// definitions for parammm library.

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

#include "pswitch.hh"

using std::string;

namespace parammm
{

  pswitch::pswitch(const string &longopt,
		   char shortopt,
		   const switch_opt &sopt,
		   const string &switchdescr,
		   const string &switchoptdescr)
    : m_longopt(longopt), m_shortopt(shortopt),
      m_sopt(sopt.makecopy()), m_switchdescr(switchdescr),
      m_switchoptdescr(switchoptdescr)
  {
  }

  pswitch::pswitch(const pswitch &s)
    : m_longopt(s.m_longopt), m_shortopt(s.m_shortopt),
      m_sopt(s.m_sopt->makecopy()), m_switchdescr(s.m_switchdescr),
      m_switchoptdescr(s.m_switchoptdescr)
  {
  }

  pswitch::~pswitch()
  {
    delete m_sopt;
  }

  pswitch& pswitch::operator=(const pswitch &s)
  {
    delete m_sopt;

    m_longopt = s.m_longopt;
    m_shortopt = s.m_shortopt;
    m_sopt = s.m_sopt->makecopy();
    m_switchdescr = s.m_switchdescr;

    return *this;
  }


}
