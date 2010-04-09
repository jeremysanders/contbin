// Auto-help class
//  Used internally to make autohelp message

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


#ifndef PARAMMM_AUTOHELP_HH
#define PARAMMM_AUTOHELP_HH

#include <string>
#include <vector>
#include <iostream>
#include "pswitch.hh"

namespace parammm {

  class autohelp_opt : public switch_opt
  {
  public:
    autohelp_opt(const std::vector<pswitch> &switches,
		 const std::string &progdescr,  // descript at top
		 const std::string &prognotes); // notes at bottom

    bool takesoption() const;
    autohelp_opt* makecopy() const;

  private:
    void setfromstream(std::istream *) const;

    const std::vector<pswitch> &m_switches;
    std::string m_progdescr, m_prognotes;
  };

}

#endif
