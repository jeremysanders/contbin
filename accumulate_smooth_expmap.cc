// accumulate-smoothing with exposure map
// Jeremy Sanders 2002-2005
// Released under the GNU Public License

#include <vector>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <string>

#include "parammm/parammm.hh"

#include "misc.hh"
#include "fitsio_simple.hh"

using namespace std;

class flux_estimator
{
public:
  flux_estimator(const image_float* const in_image,
		 const image_float* const expmap_image,
		 const image_short* const mask_image,
		 const double minsn );

  const image_float& operator()();

  struct _point
  {
    _point(int xp, int yp) : x(xp), y(yp) {}
    int x, y;
  };
  typedef std::vector<_point> _point_vec;
  typedef std::vector<_point_vec> _point_vec_vec;

private:
  void do_estimation();

  // work out which points are in which annuli
  void precalculate_annuli();
  void smooth();

  void bin();

private:

  const unsigned _xw, _yw; // dimensions of images
  const double _minsn; // minimum signal:noise

  const image_float* const _in_image;
  const image_float* const _expmap_image;
  const image_short* const _mask_image;

  // precalculated list of which points are in which annuli
  const unsigned _max_annuli;
  _point_vec_vec _annuli_points;

  bool _done;

  image_float _out_image; // output image
  image_float _estimated_errors; // errors on iteration
};

/////////////////////////////////////////////////////////////////////////

// work out integerised radius
inline static unsigned unsigned_radius(int x, int y)
{
  return unsigned( sqrt( double(x*x + y*y) ) );
}

flux_estimator::flux_estimator( const image_float* const in_image,
				const image_float* const expmap_image,
				const image_short* const mask_image,
				const double minsn = 10  )
  : _xw( in_image->xw() ), _yw( in_image->yw() ),
    _minsn( minsn ),
    _in_image(in_image), _expmap_image(expmap_image),
    _mask_image(mask_image),
    _max_annuli( unsigned_radius(_xw, _yw)+1 ),
    _annuli_points( _max_annuli ),
    _done( false ),
    _out_image( _xw, _yw ),
    _estimated_errors( _xw, _yw )
{
  // check background image is the same size and input image
  assert( expmap_image->xw() == _xw && expmap_image->yw() == _yw );
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

  return _out_image;
}

void flux_estimator::do_estimation()
{
  smooth();
}

static inline double _square(const double d)
{
  return d*d;
}

void flux_estimator::smooth()
{
  static int c = 0;

  const double SN_2 = _minsn*_minsn;

  // iterate over each pixel
  for(unsigned y=0; y != _yw; ++y)
    {
      if( y % (_yw/10) == 0 )
	{
	  std::cout << y / (_yw/10) << ' ';
	  std::cout.flush();
	}

      for(unsigned x=0; x != _xw; ++x)
	{
	  // skip masked pixels
	  if( (*_mask_image)(x, y) < 1 )
	    continue;

	  double sum_corrected = 0;

	  double foreground = 0;
	  double noise_2 = 0;
	  unsigned count = 0;
	  unsigned lastcount = 0;

	  unsigned radius = 0;

	  // loop over pixels until signal to noise >= _minsn
	  while ( radius < _max_annuli &&
		  (noise_2 == 0. ||
		   ( _square(foreground) / noise_2
		     < SN_2 )) )
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

		  const double corrected = (*_in_image)(xp, yp);
		  const double expmap = (*_expmap_image)(xp, yp);

		  const double in = corrected*expmap;
		  foreground += in;
		  noise_2 += in;
		  sum_corrected += corrected;
		  count++;
		}
	      
// 	      if( radius > 20 && lastcount == count )
// 		break;

	      // next shell
	      lastcount = count;
	      radius++;
	    }
	  
	  _out_image(x, y) = sum_corrected/count;
	  _estimated_errors(x, y) = sqrt( foreground ) / count;

	}
    }

  std::cout << '\n'; c++;
}

////////////////////////////////////////////////////////////////////////////

// accumulate smoothing

template<class T> void load_image(const string& filename,
				  T** image)
{
  FITSFile ds(filename);
  ds.readImage(image);
}

void write_image(const string& filename, const image_float& img)
{
  FITSFile ds(filename, FITSFile::Create);
  ds.writeImage(img);
}

int main(int argc, char* argv[])
{
  string back_file, mask_file;
  string out_file = "acsmooth.fits";
  double sn = 15;

  parammm::param params(argc, argv);
  params.add_switch( parammm::pswitch( "mask", 'm',
				       parammm::pstring_opt(&mask_file),
				       "set mask file",
				       "FILE"));
  params.add_switch( parammm::pswitch( "out", 'o',
				       parammm::pstring_opt(&out_file),
				       "set output file (def acsmooth.fits)",
				       "FILE"));
  params.add_switch( parammm::pswitch("sn", 's',
				      parammm::pdouble_opt(&sn),
				      "set signal:noise threshold (def 15)",
				      "VAL"));
  params.set_autohelp("Usage: accumulate_smooth_expmap [OPTIONS] expcorrect.fits expmap.fits\n"
		      "Accumulate smoothing program (exposure map).\n"
		      "Written by Jeremy Sanders 2004.",
		      "Report bugs to <jss@ast.cam.ac.uk>");
  params.enable_autohelp();
  params.enable_autoversion("0.1",
			    "Jeremy Sanders",
			    "Licenced under the GPL - see the file COPYING");
  params.enable_at_expansion();
  params.interpret_and_catch();

  if(params.args().size() != 2)
    {
      params.show_autohelp();
    }

  image_float* in_image;
  load_image( params.args()[0], &in_image);

  image_float* expmap_image;
  load_image( params.args()[1], &expmap_image );

  image_short* mask_image;
  if( mask_file.empty() )
    {
      mask_image = new image_short( in_image->xw(), in_image->yw(), 1 );
    }
  else
    {
      load_image( mask_file, &mask_image );
    }

  flux_estimator fe( in_image, expmap_image, mask_image, sn);
  image_float out = fe();

  write_image(out_file, out);

  delete in_image;
  delete expmap_image;
  delete mask_image;

  return 0;
}
