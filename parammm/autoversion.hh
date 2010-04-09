// Auto-version class
// UNUSED CODE at the moment

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

#ifndef PARAMMM_AUTOVERSION_HH
#define PARAMMM_AUTOVERSION_HH

#include <string>
#include <iostream>
#include "pswitch.hh"

namespace parammm {

  class autoversion_opt : public switch_opt
  {
  public:
    autoversion_opt(const std::string &progver,
		    const std::string &authors,
		    const std::string &licence);

    bool takesoption() const;
    autoversion_opt* makecopy() const;

  private:
    void setfromstream(std::istream *) const;

    std::string m_progver, m_authors, m_licence;
  };

}

#endif
