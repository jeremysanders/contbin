// exceptions for parammm library

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

#ifndef PARAMMM_EXCEPTIONS_HH
#define PARAMMM_EXCEPTIONS_HH

#include <string>

namespace parammm
{
  class except
  {
  public:
  };

  class except_returnstr : public except
  {
  public:
    except_returnstr(const std::string& uswitch) : m_switch(uswitch) {}
    const std::string &operator()() const { return m_switch; }
  private:
    const std::string m_switch;
  };

  class except_undef_switch : public except_returnstr
  {
  public:
    except_undef_switch(const std::string &uswitch)
      : except_returnstr(uswitch) {}
  };

  class except_xs_switch_param : public except_returnstr
  {
  public:
    except_xs_switch_param(const std::string &uswitch)
      : except_returnstr(uswitch) {}
  };

  class except_no_switch_param : public except_returnstr
  {
  public:
    except_no_switch_param(const std::string &uswitch)
      : except_returnstr(uswitch) {}
  };

  class except_inv_switch_option : public except_returnstr
  {
  public:
    except_inv_switch_option(const std::string &switcho)
      : except_returnstr(switcho) {}
  };

  class except_invalid_at_file : public except_returnstr
  {
  public:
    except_invalid_at_file(const std::string &filename)
      : except_returnstr(filename) {}
  };

}

#endif
