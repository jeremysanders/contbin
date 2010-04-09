// Parameter library interface
//  Main object defn for parammm library

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

#ifndef PARAMMM_PARAM_HH
#define PARAMMM_PARAM_HH

#include <vector>
#include "pswitch.hh"
#include "exceptions.hh"
#include "switch_option.hh"
#include "autohelp.hh"
#include "autoversion.hh"

namespace parammm
{
  typedef std::vector<std::string> str_vec;

  class param
  {
  public:
    param(int argc, const char * const *argv);
    param(const str_vec &vec);
    ~param();

    void add_switch(const pswitch &);

    void set_autohelp(const std::string &program_description,
		      const std::string &program_notes);
    void enable_autohelp();
    void enable_at_expansion();
    void show_autohelp();
    void enable_autoversion(const std::string &progver,
			    const std::string &authors,
			    const std::string &licence);

    void interpret();
    void interpret_and_catch();  // catch exception to print autohelp
    const str_vec& args();

  private:
    void addarg(const std::string &);
    void addlongopt(const std::string &opt, const std::string &next,
		    bool *moveon);
    void addshortopt(const std::string &opt, const std::string &next,
		     bool *moveon);

    void expand_at_file(const std::string &filename, int where);

  private:
    //    int m_argc;
    //    const char * const *m_argv;

    std::vector<std::string> m_argv;

    bool m_gotdoubledash;
    std::string m_progdescr;
    std::string m_prognotes;
    bool m_at_expansion;

    std::vector<pswitch> m_switches;
    str_vec m_args;
  };

  //////////////////////////////////////

  inline const str_vec& param::args()
  {
    return m_args;
  }

}

#endif
