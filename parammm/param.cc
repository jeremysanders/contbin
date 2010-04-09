// Param implementation for parammm library

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
#include <fstream>
#include <algorithm>

#include "exceptions.hh"
#include "param.hh"

using std::string;
using std::cerr;
using std::endl;
using std::vector;
using std::ifstream;

namespace parammm
{

  param::param(int argc, const char * const *argv)
    : m_at_expansion(false)
  {
    // copy parameters into m_argv
    for(int i=1; i<argc; ++i)
      m_argv.push_back( argv[i] );
  }

  param::param(const str_vec &vec)
    : m_argv(vec),
      m_at_expansion(false)
  {
  }

  param::~param()
  {
  }

  void param::add_switch(const pswitch &s)
  {
    m_switches.push_back(s);
  }

  void param::set_autohelp(const string &progdescr,
			   const string &prognotes)
  {
    m_progdescr = progdescr;
    m_prognotes = prognotes;
  }

  void param::enable_autohelp()
  {
    add_switch( pswitch("help", '\0', autohelp_opt(m_switches, m_progdescr,
						   m_prognotes),
			"display this help message", "") );
  }
  
  void param::enable_autoversion(const string &progver,
				 const string &authors,
				 const string &licence)
  {
    add_switch( pswitch("version", 'V', autoversion_opt(progver,
							authors, licence),
			"display the program version", "") );
  }

  void param::show_autohelp()
  {
    autohelp_opt help(m_switches, m_progdescr, m_prognotes);
    help.interpret("");
  }

  void param::interpret_and_catch()
  {
    try {
      interpret();
    }
    catch(const except_invalid_at_file e) {
      cerr << "Cannot open parameter file " << e() << "\n\n";
      show_autohelp();
    }
    catch(const except_returnstr e) {
      cerr << "Option " << e() << " invalid\n\n";
      show_autohelp();
    }
  }

  void param::interpret()
  {
    m_gotdoubledash = false;

    int argcnt = 0;
    while( argcnt < int(m_argv.size()) ) {
      const string currarg = m_argv[argcnt];

      argcnt++;              // set to next number (first!!)
      const string nextarg = ( argcnt<int(m_argv.size()) ) ?
	m_argv[argcnt] : "";

      // single character, or had double-dash before
      if(currarg.length() <= 1 || m_gotdoubledash) {
	addarg(currarg);
	continue;
      }

      // long option
      if(currarg[0] == '-' && currarg[1] == '-') {
	bool moveon = false;
	addlongopt(currarg.substr(2), nextarg, &moveon);
	if(moveon) argcnt++;
	continue;
      }

      // short option
      if(currarg[0] == '-') {
	bool moveon = false;
	addshortopt(currarg.substr(1), nextarg, &moveon);
	if(moveon) argcnt++;
	continue;
      }

      // at expansion
      if(currarg[0] == '@' && m_at_expansion) {
	expand_at_file( currarg.substr(1), argcnt );
	continue;
      }

      addarg(currarg);
    } // go to next argument...

  }

  void param::addarg(const string &arg)
  {
    m_args.push_back(arg);
  }

  void param::addlongopt(const string &opt, const string &next,
			 bool *moveon)
  {
    if(opt.empty()) { // option
      m_gotdoubledash = true;
      return;
    }

    // want option before '=' sign, if any
    const size_t posn = opt.find('=');

    const string strippedopt = (posn == string::npos) ? opt :
      opt.substr(0, posn);

    const vector<pswitch>::iterator switchp = find(m_switches.begin(),
						   m_switches.end(),
						   strippedopt);

    if( switchp == m_switches.end() )
      {
	throw except_undef_switch( strippedopt );
      }

    const bool takesoption = switchp -> takes_switch_option();

    // it doesn't take an option and one's present
    if( !takesoption && posn != string::npos )
      throw except_xs_switch_param( strippedopt );

    // it takes one and there doesn't seem to be one (or next
    // param starts with a '-')
    if( takesoption && posn == string::npos && (next.empty() || next[0]=='-') )
      throw except_no_switch_param( strippedopt );

    // we have tested it for rough conformance...
    // now we 'run' the option

    string optparam;
    if( takesoption ) {
      if( posn == string::npos ) {
	optparam = next;
	*moveon = true;
      } else
	optparam = opt.substr(posn+1);
    }

    // now it's all up to the option to handle
    switchp->interpret_switch_option(optparam);
  }

  void param::addshortopt(const string &opt, const string &next,
			  bool *moveon)
  {
    const unsigned posn = opt.find('=');
    const unsigned len = opt.length();

    for(unsigned i=0; i < len && opt[i] != '='; i++) {
      // string version of option
      const string optstr( string("")+opt[i] );

      const vector<pswitch>::iterator switchp = find(m_switches.begin(),
						     m_switches.end(), opt[i]);
      if(switchp==m_switches.end())
	throw except_undef_switch( optstr );

      const bool takesoption = switchp -> takes_switch_option();

      // doesn't take option and there's an '=' next
      if( !takesoption && posn==i+1 )
	throw except_xs_switch_param( optstr );

      // takes an option but we're not at the last char and there's no '='
      if( takesoption && i<len-1 && posn!=i+1 ) 
	throw except_no_switch_param( optstr );

      // takes an option, we're at the last char, but there's no next arg,
      // or starts with a '-'
      if( takesoption && i==len-1 && (next.empty() || next[0]=='-') )
	throw except_no_switch_param( optstr );

      string optparam;
      if( takesoption ) {
	if( posn == i+1 )    // an equals is next
	  optparam=opt.substr(posn + 1);
	else {
	  optparam = next;
	  *moveon = true;
	}
      }
      switchp->interpret_switch_option(optparam);

    } // end loop

  } // end fn

  void param::enable_at_expansion()
  {
    m_at_expansion = true;
  }

  void param::expand_at_file(const string &filename, int where)
  {
    ifstream at_file(filename.c_str());
    vector<string> args;

    if( ! at_file )
      throw except_invalid_at_file( filename );

    while( ! at_file.eof() ) {
      string line;
      getline(at_file, line);

      if( line.empty() ) continue;
      if( line[0] == '#' ) continue;

      // break line up into WS
      // bad algorithm, needs rewriting
      const int len = line.size();
      bool in_quote = false;

      int i = 0;
      string temp;
      while(i < len) {
	const char c = line[i];

	if( c == '"' ) {
	  if( i > 0 && line[i-1] == '\\' ) {
	    temp[ temp.size() - 1 ] = c;
	  } else {
	    in_quote = ! in_quote;
	  }
	  i++;
	  continue;
	}

	if( c == '#' && !in_quote ) break;  // ignore comments

	if( (c == ' ' || c == '\t') && !in_quote ) {
	  if( ! temp.empty() )
	    args.push_back(temp);
	  temp.erase();
	} else {
	  temp += c;
	}
	i++;
      }
      if(!temp.empty()) args.push_back(temp);

    }

    m_argv.insert(m_argv.begin()+where, args.begin(), args.end());
  }

}
