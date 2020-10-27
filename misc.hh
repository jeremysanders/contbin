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
typedef dm::memimage<bool> image_bool;

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

// A simple RIAA pointer class to ensure pointers are deleted at end
// of scope.  This should be replaced by something standard for C++11.
template<class T> class delete_ptr
{
public:
  delete_ptr(T* p=0) : _ptr(p) {}
  ~delete_ptr() { delete _ptr; }

  T& operator*() const { return *_ptr; }
  T* operator->() const { return _ptr; }
  void operator=(T* p) { _ptr = p; }  // note this doesn't delete!

  T* ptr() const { return _ptr; }
  T** pptr() { return &_ptr; }

private:
  T* _ptr;
};

#endif
