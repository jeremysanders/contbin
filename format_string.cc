#include <vector>
#include <stdexcept>
#include "format_string.hh"

using namespace std;

static vector<string> do_splitting(const string& s)
{
  string build;
  vector<string> items;

  // keep adding characters and building strings
  const string::const_iterator e = s.end();
  for( string::const_iterator i = s.begin(); i != e; ++i )
    {
      const char c = *i;
      // end of current item
      if( c == 0 )
	{
	  items.push_back(build);
	  build.clear();
	}
      else
	{
	  build += c;
	}
    }

  // add remainder if reqd
  if( ! build.empty() )
    items.push_back(build);

  return items;
}

static string do_formatting( const string& format,
			     const vector<string>& items )
{
  string out;

  // iterate over chars in format string
  string::const_iterator i = format.begin();
  const string::const_iterator  e = format.end();
  while( i != e )
    {
      const char c = *i;
      if( c == '%' )
	{
	  // turn %% into %
	  if( i+1 != e && *(i+1) == '%' )
	    {
	      out.push_back('%');
	      i += 2;
	      continue;
	    }

	  // get a number
	  string no;
	  i++;
	  while( i != e && *i >= '0' && *i <= '9' )
	    {
	      no.push_back(*i);
	      i++;
	    }

	  // didn't find one
	  if( no.empty() )
	    throw invalid_argument("Util::FormatString:"
				   " invalid character after '%'");

	  // convert number to int
	  istringstream s(no);
	  int index = -1;
	  s >> index;
	  index -= 1;

	  // add on the string to the build
	  if( index >= 0 && index < int(items.size()) )
	    {
	      out += items[index];
	    }
	  else
	    {
	      throw out_of_range("Util::FormatString: index out of range");
	    }
	}
      else
	{
	  out.push_back(c);
	  i++;
	}
    }

  return out;
}

void Util::FormatString::format()
{
  // split up items
  const vector<string> items = do_splitting( this->str() );

  _output = do_formatting( _format, items );
  _done_format = true;
}
