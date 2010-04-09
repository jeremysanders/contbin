#include <cassert>
#include <iostream>
#include <algorithm>
#include <iomanip>

#include "flux_estimator.hh"

using namespace std;

// work out integerised radius
inline static unsigned unsigned_radius(int x, int y)
{
  return unsigned( sqrt( double(x*x + y*y) ) );
}

// simple squaring function
template<class T> inline static T square(T v)
{
  return v*v;
}

// estimate the error squared on c counts
// uses formula of Gehrels 1986 ApJ, 303, 336) eqn 7
inline static double error_sqd_est(double c)
{
  return square( 1. + sqrt(c + 0.75) );
}

/////////////////////////////////////////////////////////////////////////

flux_estimator::flux_estimator( const image_float* const in_image,
				const image_float* const back_image,
				const image_short* const mask_image,
				const image_float* const expmap_image,
				const image_float* const bg_expmap_image,
				const image_float* const noisemap_image,
				
				const double minsn = 10  )
  : _xw( in_image->xw() ), _yw( in_image->yw() ),
    _minsn( minsn ),
    _in_image(in_image), _back_image(back_image), _mask_image(mask_image),
    _expmap_image(expmap_image), _bg_expmap_image(bg_expmap_image),
    _noisemap_image(noisemap_image),
    _max_annuli( unsigned_radius(_xw, _yw)+1 ),
    _annuli_points( _max_annuli ),
    _done( false ),
    _iteration_image( _xw, _yw ),
    _estimated_errors( _xw, _yw )
{
  // check background image is the same size and input image
  assert( back_image == 0 ||
	  (back_image->xw() == _xw && back_image->yw() == _yw ) );
  assert( mask_image->xw() == _xw && mask_image->yw() == _yw );
}

void flux_estimator::precalculate_annuli()
{
  for(int y=-(_yw-1); y<int(_yw); ++y)
    {
      for(int x=-(_xw-1); x<int(_xw); ++x)
	{
	  const unsigned r = unsigned_radius(x, y);
	  
	  // add pixel to appropriate radius
	  _annuli_points[r].push_back( _point(x, y) );
	}
    }
}

const image_float& flux_estimator::operator()()
{
  if( ! _done )
    {
      precalculate_annuli();
      do_estimation();
      _done = true;
    }

  return _iteration_image;
}

void flux_estimator::do_estimation()
{
  smooth();
}

struct _pixel_trim
{
  unsigned x, y;
  double foreground, background;

  _pixel_trim(unsigned _x, unsigned _y, double _foreground,
	      double _background)
    : x(_x), y(_y), foreground(_foreground), background(_background)
  {}
};

class _compare_pixel_trim
{
  const double back_factor;

public:
  _compare_pixel_trim(double bf) : back_factor(bf)
  {};

  int operator()(const _pixel_trim& p1, const _pixel_trim& p2)
  {
    return (p1.foreground - p2.foreground -
      (p1.background - p2.background)*back_factor) < 0;
  }
};

static inline double _square(const double d)
{
  return d*d;
}

void flux_estimator::smooth()
{
  static int c = 0;

  const double min_sn_2 = _minsn*_minsn;

  // iterate over each pixel
  for(unsigned y=0; y != _yw; ++y)
    {
      // write out percentage on each y iteration
      std::cout << '\r'
		<< std::setprecision(3)
		<< std::setw(8)
		<< std::showpoint
		<< y * 100. / _yw
		<< "%";
      std::cout.flush();

      for(unsigned x=0; x != _xw; ++x)
	{
	  // skip masked pixels
	  if( (*_mask_image)(x, y) < 1 )
	    continue;

	  double fg_sum = 0;
	  double bg_sum = 0;
	  double bg_sum_weight = 0;
	  double expratio_sum_2 = 0;
	  double sn_2 = 0;
	  double noise_2 = 0;

	  // if a noise map is supplied
	  double noise_2_total = 0.;

	  // keep track of number of pixels and radius
	  unsigned count = 0;
	  unsigned radius = 0;

	  // loop over pixels until signal to noise >= _minsn
	  while ( radius < _max_annuli && sn_2 < min_sn_2 )
	    {
	      // iterate over points in radius
	      const _point_vec::const_iterator e =
		_annuli_points[radius].end();
	      for( _point_vec::const_iterator i =
		     _annuli_points[radius].begin(); i != e; ++i )
		{
		  const int xp = int(x) + i->x;
		  const int yp = int(y) + i->y;
		  // skip pixels we don't have
		  if( xp < 0 || yp < 0 || xp >= int(_xw) || yp >= int(_yw) )
		    continue;
		  
		  // skip masked pixels
		  if( (*_mask_image)(xp, yp) < 1 )
		    continue;

		  const double in = (*_in_image)(xp, yp);

		  // count up background if any
		  if( _back_image != 0 )
		    {
		      const double bg = (*_back_image)(xp, yp);
		      const double expratio = (*_expmap_image)(xp, yp) /
			(*_bg_expmap_image)(xp, yp);
		      bg_sum += bg;
		      bg_sum_weight += bg*expratio;
		      expratio_sum_2 += expratio*expratio;
		    }		      

		  // add up noise if supplied
		  if( _noisemap_image != 0 )
		    {
		      noise_2_total += square((*_noisemap_image)(xp, yp));
		    }

		  fg_sum += in;
		  count++;
		}
	      
	      // next shell

	      if ( _noisemap_image != 0 )
		{
		  // we have a noise map
		  noise_2 = noise_2_total;
		}
	      else
		{
		  // calculate noise
		  noise_2 = error_sqd_est(fg_sum);

		  if( _back_image != 0 )
		    noise_2 += (expratio_sum_2 / count) * error_sqd_est(bg_sum);
		}

	      sn_2 = square(fg_sum - bg_sum_weight) / noise_2;
	      radius++;
	    }
	  
	  _iteration_image(x, y) = (fg_sum - bg_sum_weight)/count;
	  _estimated_errors(x, y) = sqrt( noise_2 );
	}
    }

  std::cout << '\n'; c++;
}
