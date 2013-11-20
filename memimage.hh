#ifndef DM_MEMIMAGE_HH
#define DM_MEMIMAGE_HH

#include <valarray>
#include <string>
#include <limits>
#include <algorithm>

namespace dm {

  template<class T> class memimage
  {
  public:
    // blank image
    memimage(const unsigned xw, const unsigned yw, const T val = 0)
      : m_xw(xw), m_yw(yw), m_data( val, xw*yw )
    {}

    // copy image from another
    memimage(const memimage<T>& other)
      : m_xw( other.m_xw ), m_yw( other.m_yw ), m_data( other.m_data )
    {}

    // initialise from C-style array
    memimage(const unsigned xw, const unsigned yw, const T* data)
      : m_xw( xw ), m_yw( yw ), m_data( data, xw*yw )
    {}

    // initialise from other datatypes
    template<class T2> explicit memimage(const memimage<T2>& other);
    template<class T2> memimage(const unsigned xw,
				const unsigned yw, const T2* const data);

    // read from dumpdata file
    memimage(const std::string& filename);
    // dump to file
    void dump_to_file(const std::string &filename) const;

    // set all the pixels
    void set_all(const T val = 0)  { m_data = val; }

    // get access to pixels
    T& operator() (const unsigned x, const unsigned y)
     { return m_data[x+y*m_xw]; }
    T operator() (const unsigned x, const unsigned y) const
     { return m_data[x+y*m_xw]; }

    // get flat access to pixels
    T& flatdata(const unsigned i)
    { return m_data[i]; }
    T flatdata(const unsigned i) const
    { return m_data[i]; }

    // checked access to pixels
    class out_of_range_exception {};
    T& pixel(const unsigned x, const unsigned y);
    T pixel(const unsigned x, const unsigned y) const;

    class size_mismatch_exception {};

    // get maximum value in image
    T max() const
    {
      // get minimum value (not smallest value)
      T maxval = std::numeric_limits<T>::is_integer
	? std::numeric_limits<T>::min()
	: -std::numeric_limits<T>::max();

      const size_t len = m_xw*m_yw;
      for( size_t i = 0; i != len; ++i )
	maxval = std::max( maxval, m_data[i] );
      return maxval;
    }

    // get minimum value of image
    T min() const
    {
      T minval = std::numeric_limits<T>::max();

      const size_t len = m_xw*m_yw;
      for( size_t i = 0; i != len; ++i )
	minval = std::min( minval, m_data[i] );
      return minval;
    }

    T sum() const
    {
      T tot = 0;
      const size_t len = m_xw*m_yw;
      for( size_t i = 0; i != len; ++i )
	tot += m_data[i];
      return tot;
    }

    // make all values <= upperval
    void trim_down(const T upperval)
    {
      const size_t len = m_xw*m_yw;
      for( size_t i = 0; i != len; ++i )
	m_data[i] = std::min( m_data[i], upperval );
    }

    // make all values >= lowerval
    void trim_up(const T lowerval)
    {
      const size_t len = m_xw*m_yw;
      for( size_t i = 0; i != len; ++i )
	m_data[i] = std::max( m_data[i], lowerval );
    }

  private:
    void assert_size_other(const memimage<T>& other) const
    {
      if( m_xw != other.xw() || m_yw != other.yw() )
	throw size_mismatch_exception();
    }

  public:
    // various operations with images
    //  multiply image by another
    const memimage<T>& operator *= (const memimage<T>& other)
    {
      assert_size_other(other); m_data *= other.data(); return *this;
    }
    const memimage<T>& operator /= (const memimage<T>& other)
    {
      assert_size_other(other); m_data /= other.data(); return *this;
    }
    const memimage<T>& operator -= (const memimage<T>& other)
    {
      assert_size_other(other); m_data -= other.data(); return *this;
    }
    const memimage<T>& operator += (const memimage<T>& other)
    {
      assert_size_other(other); m_data += other.data(); return *this;
    }
    
    // with constants
    const memimage<T>& operator *= (const T other)
    {
      m_data *= other; return *this;
    }
    const memimage<T>& operator /= (const T other)
    {
      m_data /= other; return *this;
    }
    const memimage<T>& operator -= (const T other)
    {
      m_data -= other; return *this;
    }
    const memimage<T>& operator += (const T other)
    {
      m_data += other; return *this;
    }

    // return information about the image
    unsigned xw() const { return m_xw; }  // return width
    unsigned yw() const { return m_yw; }  // return height
    unsigned nelem() const { return m_xw*m_yw; } // no elements
    const std::valarray<T>& data() const { return m_data; } // return data

    // these make temporaries
    memimage<T> operator * (const memimage<T>& other) const
    {
      memimage<T> t = *this; t *= other; return t;
    }
    memimage<T> operator / (const memimage<T>& other) const
    {
      memimage<T> t = *this; t /= other; return t;
    }
    memimage<T> operator - (const memimage<T>& other) const
    {
      memimage<T> t = *this; t -= other; return t;
    }
    memimage<T> operator + (const memimage<T>& other) const
    {
      memimage<T> t = *this; t += other; return t;
    }

    // with constants
    memimage<T> operator * (const T other) const
    {
      memimage<T> t = *this; t *= other; return t;
    }
    memimage<T> operator / (const T other) const
    {
      memimage<T> t = *this; t /= other; return t;
    }
    memimage<T> operator + (const T other) const
    {
      memimage<T> t = *this; t += other; return t;
    }
    memimage<T> operator - (const T other) const
    {
      memimage<T> t = *this; t -= other; return t;
    }
    
  private:
    unsigned m_xw, m_yw;
    std::valarray<T> m_data;
  };


  } // namespace

// convert types
// put this in the header, as we'd generate 10^6 versions of this
// otherwise.
template<class T> template<class T2>
dm::memimage<T>::memimage(const unsigned xw, const unsigned yw,
			  const T2* const data)
  : m_xw(xw), m_yw(yw),
    m_data( data, nelem() )
{
  const unsigned size = nelem();
  for(unsigned i=0; i != size; ++i)
    m_data[i] = static_cast<T>( data[i] );
}

// copy image of different type
template<class T> template<class T2>
dm::memimage<T>::memimage(const memimage<T2>& other)
  : m_xw( other.xw() ), m_yw( other.yw() ),
    m_data( static_cast<T>(0), nelem() )
{
  const unsigned size = nelem();
  const std::valarray<T2>& otherdata = other.data();

  for(unsigned i=0; i != size; ++i)
    m_data[i] = static_cast<T>(otherdata[i]);
}

#endif
