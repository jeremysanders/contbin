// Switch definitions
// for parammm library

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

#ifndef PARAMMM_PSWITCH_HH
#define PARAMMM_PSWITCH_HH

#include <string>
#include "switch_option.hh"

namespace parammm
{

  class pswitch
  {
  public:
    pswitch(const std::string &longopt,      // long option (without --)
	    char shortopt,              // short option (or \0)
	    const switch_opt &sopt,     // return val of opt
	    const std::string &switchdescr,  // description of switch
	    const std::string &switchoptdescr);  // opt descr (eg FILE)
    pswitch(const pswitch &);
    ~pswitch();
    bool takes_switch_option() const;
    const std::string& long_option() const;
    char short_option() const;
    const std::string& description() const;
    const std::string& switchopt_description() const;
    void interpret_switch_option(const std::string &swopt) const;

    // check to see whether selected
    bool operator==(const std::string &opt) const;
    bool operator!=(const std::string &opt) const;
    bool operator==(char opt) const;
    bool operator!=(char opt) const;
    pswitch& operator=(const pswitch &);

  private:
    std::string m_longopt;
    char m_shortopt;
    switch_opt *m_sopt;
    std::string m_switchdescr;
    std::string m_switchoptdescr;
  };

  /////////////////////////////////////////////////////////

  inline bool pswitch::takes_switch_option() const
  {
    return m_sopt->takesoption();
  }

  inline const std::string& pswitch::long_option() const
  {
    return m_longopt;
  }

  inline char pswitch::short_option() const
  {
    return m_shortopt;
  }

  inline const std::string& pswitch::description() const
  {
    return m_switchdescr;
  }

  inline const std::string& pswitch::switchopt_description() const
  {
    return m_switchoptdescr;
  }

  inline void pswitch::interpret_switch_option(const std::string &swopt) const
  {
    m_sopt->interpret(swopt);
  }

  inline bool pswitch::operator==(const std::string &opt) const
  {
    return(m_longopt == opt);
  }

  inline bool pswitch::operator!=(const std::string &opt) const
  {
    return(m_longopt != opt);
  }

  inline bool pswitch::operator==(char opt) const
  {
    return(m_shortopt == opt);
  }

  inline bool pswitch::operator!=(char opt) const
  {
    return(m_shortopt != opt);
  }

}

#endif
