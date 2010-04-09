// Autoversion implementation
// 

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

#include <iostream>
#include <cstdlib>
#include "autoversion.hh"

using std::string;
using std::cout;
using std::cerr;
using std::clog;
using std::endl;
using std::exit;

namespace parammm {

  autoversion_opt::autoversion_opt(const string &progver,
				   const string &authors,
				   const string &licence)
    : m_progver(progver), m_authors(authors),
      m_licence(licence)
  {
  }

  bool autoversion_opt::takesoption() const
  {
    return false;
  }

  autoversion_opt* autoversion_opt::makecopy() const
  {
    return new autoversion_opt(m_progver, m_authors, m_licence);
  }

  void autoversion_opt::setfromstream(std::istream *s) const
  {
    clog << "Version " << m_progver << endl
	 << "Copyright (C) " << m_authors << endl;
    if(! m_licence.empty() )
      clog << m_licence << endl;
    
    exit(1);
  }

}
