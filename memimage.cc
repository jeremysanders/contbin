#include <iostream>
#include <fstream>
#include "memimage.hh"

template<class T> dm::memimage<T>::memimage(const std::string& filename)
{
  std::ifstream file(filename.c_str());

  file >> m_xw >> m_yw;

  m_data.resize(m_xw*m_yw);
  if( file )
    {
      for(unsigned y=0; y<m_yw; ++y)
	for(unsigned x=0; x<m_xw; ++x)
	  {
	    file >> operator()(x, y);
	  }
    }

  if( ! file )
    {
      std::cerr << "Error reading from file\n";
      std::exit(1);
    }

}

template<class T> void dm::memimage<T>::dump_to_file
(const std::string& filename) const
{
  std::ofstream file(filename.c_str());

  file << xw() << ' ' << yw() << '\n';
  for(unsigned y=0; y<yw(); ++y)
    {
      for(unsigned x=0; x<xw(); ++x)
	{
	  file << operator()(x, y) << ' ';
	}
      file << '\n';
    }
}

// checked pixel access
template<class T> T& dm::memimage<T>::pixel(unsigned x, unsigned y)
{
  if( x >= m_xw || y >= m_yw)
    throw out_of_range_exception();
  return operator() (x, y);
}

// checked const pixel access
template<class T> T dm::memimage<T>::pixel(unsigned x, unsigned y) const
{
  if( x >= m_xw || y >= m_yw)
    throw out_of_range_exception();

  return operator() (x, y);
}

#define DM_DEFINE_TEMPL(TYPE) \
   template class dm::memimage<TYPE>;

DM_DEFINE_TEMPL(short)
DM_DEFINE_TEMPL(long)
DM_DEFINE_TEMPL(float)
DM_DEFINE_TEMPL(double)
DM_DEFINE_TEMPL(unsigned char)
DM_DEFINE_TEMPL(unsigned short)
DM_DEFINE_TEMPL(unsigned long)

#undef DM_DEFINE_TEMPL
