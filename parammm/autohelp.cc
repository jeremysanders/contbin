// Autohelp implementation for parammm

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

#include <algorithm>
#include <iostream>
#include "autohelp.hh"

using namespace std;

namespace parammm {

  autohelp_opt::autohelp_opt(const vector<pswitch> &switches,
			     const string &progdescr,
			     const string &prognotes)
    : m_switches(switches), m_progdescr(progdescr),
      m_prognotes(prognotes)
  {
  }

  bool autohelp_opt::takesoption() const
  {
    return false;
  }

  autohelp_opt* autohelp_opt::makecopy() const
  {
    return new autohelp_opt(m_switches, m_progdescr, m_prognotes);
  }

  void displayswitch(const pswitch &s)
  {
    cerr << "  ";
    if(s.short_option() != 0)
      cerr << '-' << s.short_option();
    else
      cerr << "  ";
    if(s.short_option() != 0  && !s.long_option().empty() )
      cerr << ", ";
    else
      cerr << "  ";

    int space = 20;
    if(! s.long_option().empty() ) {
      cerr << "--" << s.long_option();
      space -= (2+s.long_option().length());
      if( ! s.switchopt_description().empty() ) {
	cerr << '=' << s.switchopt_description();
	space -= (1+s.switchopt_description().length());
      }
    }

    for(int i=space-1; i>=0; --i) cerr << ' ';

    cerr << " " << s.description() << '\n';
  }

  void autohelp_opt::setfromstream(std::istream *s) const
  {
    cerr << m_progdescr << "\n\n";
    for_each(m_switches.begin(), m_switches.end(), displayswitch);
    cerr << endl << m_prognotes << '\n';
    exit(1);
  }

}
