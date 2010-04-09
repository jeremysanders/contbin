#ifndef MISC_BINACC_HH
#define MISC_BINACC_HH

#include <vector>
#include <string>
#include <algorithm>

#include "memimage.hh"

typedef dm::memimage<double> image_dbl;
typedef dm::memimage<float> image_float;
typedef dm::memimage<long> image_long;
typedef dm::memimage<short> image_short;

typedef std::vector<unsigned> vec_unsigned;
typedef std::vector<double> vec_dbl;

// split a string into parts separated by delim
inline std::vector<std::string> split_string(const std::string& in, const char delim)
{
  std::vector<std::string> out;

  std::string::const_iterator i = in.begin();

  while(true)
    {
      std::string::const_iterator e = std::find(i, in.end(), delim);
      if(e == in.end())
	break;

      out.push_back( std::string(i, e) );
      i = e+1;
    }

  return out;
}

#endif
