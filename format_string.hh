#ifndef UTIL_FORMAT_STRING_HH
#define UTIL_FORMAT_STRING_HH

#include <string>
#include <sstream>

namespace Util
{
  // class allows strings to be formatted like:
  // Util::FormatString s("%1 there %2");
  // s << "Hi" << ends
  //   << "Fred" << ends;
  // std::cout << s << '\n';

  class FormatString : public std::ostringstream
  {
  public:
    explicit FormatString(const std::string& format)
      : _format(format), _done_format(false)
    {}
    
    // get the formatted string
    const std::string& get()
    {
      if( ! _done_format )
	format();
      return _output;
    }
    
    // other convenient ways to get the string
    operator const std::string&() { return this->get(); }
    operator const char*() { return this->get().c_str(); }

  private:
    // do the formatting
    void format();

  private:
    const std::string _format;
    
    bool _done_format;
    std::string _output;
  };

}

#endif
