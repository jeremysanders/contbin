// accumulate-smoothing with exposure map
// Jeremy Sanders 2002-2014
// Released under the GNU Public License

#include <vector>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <string>
#include <limits>

#include "parammm/parammm.hh"

#include "misc.hh"
#include "fitsio_simple.hh"
#include "image_disk_access.hh"

namespace
{
  // square function
  template <typename T> T sqr(T v) { return v*v; }

  // structure to collect x,y values
  struct point
  {
    point(int xp, int yp) : x(xp), y(yp) {}
    int x, y;
  };
  typedef std::vector<point> point_vec;
  typedef std::vector<point_vec> point_vec_vec;

  // retains more precision that a standard summation
  class KahanSum
  {
  public:
    KahanSum()
      : _sum(0), _comp(0)
    {}

    void reset() { _sum = _comp = 0.; }

    void operator +=(double val)
    {
      double y = val - _comp;
      double t = _sum + y;
      _comp = (t - _sum) - y;
      _sum = t;
    }

    double sum() const { return _sum; }

  private:
    double _sum, _comp;
  };

  // for keeping track of the smoothing process
  // adds the total cts, background and exposure
  struct State
  {
    State(double fgexp, double bgexp)
      : _fgexp(fgexp), _bgexp(bgexp)
    {
      reset();
    }

    void reset()
    {
      radius = 0;
      counts.reset();
      back.reset();
      exposure.reset();
    }

    // get square of the signal to noise
    double sn2() const
    {
      if( counts.sum() == 0. and back.sum() == 0. )
        return 0;
      else
        // note, note that the max prevents the S/N being -ve and giving a large
        // postive sn2
        return sqr(std::max((counts.sum()/_fgexp - back.sum()/_bgexp), 0.)) /
          (counts.sum()/sqr(_fgexp) + back.sum()/sqr(_bgexp));
    }

    // return background and exposure-corrected surface brightness
    double pixval() const
    {
      return (counts.sum() - back.sum()*_fgexp/_bgexp) / exposure.sum();
    }

    int radius;
    KahanSum counts;
    KahanSum back;
    KahanSum exposure;

  private:
    double _fgexp;
    double _bgexp;
  };

  // class for doing smoothing
  class Smoother
  {
  public:
    Smoother(const image_float& ct_image,
             const image_float& bg_image,
             const image_float& expmap_image,
             const image_short& mask_image,
             double exptimefg, double exptimebg,
             double sn)
      : _ct_image(ct_image),
        _bg_image(bg_image),
        _expmap_image(expmap_image),
        _mask_image(mask_image),
        _xw(ct_image.xw()), _yw(ct_image.yw()),
        _target_sn2(sqr(sn)),
        _state(exptimefg, exptimebg),
        out_image(ct_image.xw(), ct_image.yw(),
                  std::numeric_limits<double>::quiet_NaN())
    {
      precalculate_circles();
      precalculate_shifts();
      reset_state();
    }

    void smooth_all();

  private:
    void precalculate_circles();
    void precalculate_shifts();
    void reset_state();
    void add_shift(int x, int y, int r, int sign, bool doiny, bool mirror);
    void new_pixel(int x, int y);

    // this is a template to speed up the inner loop as VAL only equals -1 and 1
    template <int VAL> void add_or_remove_circle(int x, int y, int r)
    {
      const point_vec& circ(_circles[r]);
      for(auto i = circ.begin(); i != circ.end(); ++i)
        {
          auto nx = x + i->x;
          auto ny = y + i->y;
        
          if( nx >= 0 and nx < _xw and ny >= 0 and ny < _yw and _mask_image(nx, ny) )
            {
              _state.counts   += VAL*_ct_image(nx, ny);
              _state.back     += VAL*_bg_image(nx, ny);
              _state.exposure += VAL*_expmap_image(nx, ny);
            }
        }
    }

  private:
    const image_float& _ct_image;
    const image_float& _bg_image;
    const image_float& _expmap_image;
    const image_short& _mask_image;
    const int _xw;
    const int _yw;
    const double _target_sn2;

    State _state;

    point_vec_vec _circles;
    point_vec_vec _shift_incl;

    int _last_x, _last_y;

  public:
    image_float out_image;
  };

}

// precalculate pixels included as a function of radius in circle
void Smoother::precalculate_circles()
{
  auto maxrad = int(sqrt(double( sqr(_xw)+sqr(_yw) ))) + 1;
  _circles.resize(maxrad);

  for(auto y=-(_yw-1); y<_yw; ++y)
    for(auto x=-(_xw-1); x<_xw; ++x)
      {
        auto r = int( sqrt(double( sqr(x)+sqr(y) )) );
        _circles[r].push_back( point(x, y) );
      }
}

// precalculate pixels included when a circle is shifted to the right
void Smoother::precalculate_shifts()
{
  auto maxrad = int(sqrt(double( sqr(_xw)+sqr(_yw) ))) + 1;
  _shift_incl.resize(maxrad);

  for(auto y=-(_yw-1); y<_yw; ++y)
    for(auto x=-(_xw-1); x<_xw; ++x)
      {
        auto r1 = int( sqrt(double(sqr(x) + sqr(y))) );
        auto r2 = int( sqrt(double(sqr(x+1) + sqr(y))) );

        if( r1 < r2 )
          {
            _shift_incl[r1].push_back( point(x, y) );
          }
      }
}

void Smoother::reset_state()
{
  _state.reset();
  _last_x = _last_y = -9999;
}

// this includes the pixels in the sum when shifting a circle
// sign: 1 or -1, depending on whether to add or subtract the shift
// doiny: shift in y, not x
// mirror: shift to left, not right
void Smoother::add_shift(int x, int y, int r, int sign, bool doiny, bool mirror)
{
  const point_vec& shift(_shift_incl[r]);
  for(auto i = shift.begin(); i != shift.end(); ++i)
    {
      int dx = i->x;
      int dy = i->y;

      if(mirror)
        dx = -dx;

      if(doiny)
        std::swap(dx, dy);

      int nx = x + dx;
      int ny = y + dy;

      if( nx >= 0 and nx < _xw and ny >= 0 and ny < _yw and _mask_image(nx, ny) )
        {
          _state.counts   += sign*_ct_image(nx, ny);
          _state.back     += sign*_bg_image(nx, ny);
          _state.exposure += sign*_expmap_image(nx, ny);
        }
    }
}

void Smoother::new_pixel(int x, int y)
{
  bool reverse;

  // debugging
  const bool forcereset = false;

  if( ((abs(x-_last_x) == 1 and _last_y == y) or
       (_last_x == x and abs(y-_last_y) == 1)) and not forcereset )
    {
      // shift along the summed values to the next pixel
      const bool iny = _last_y != y;
      const bool mirror = _last_x > x or _last_y > y;

      // remove points from previous circle
      add_shift(_last_x, _last_y, _state.radius, -1, iny, not mirror);
      // add points from new circle
      add_shift(x, y, _state.radius, 1, iny, mirror);

      // move backwards or forwards in radius depending on current S/N
      reverse = _state.sn2() >= _target_sn2;
    }
  else
    {
      reset_state();
      // we have to add the initial pixel, as this is checked below
      add_or_remove_circle<1>(x, y, _state.radius);
      reverse = false;
    }

  // fixed radius test case
  // reset_state();
  // while(_radius < 10)
  //   {
  //     _radius++;
  //     add_circle(x, y, _radius);
  //   }

  if(not reverse)
    {
      // Standard forward increase of radius of the bin

      while(_state.sn2() < _target_sn2)
        {
          _state.radius++;
          add_or_remove_circle<1>(x, y, _state.radius);
        }
    }
  else
    {
      // Reverse direction.
      // We check whether subtracting this radius
      // pushes us under the threshold.

      auto oldsn2 = _state.sn2();
      while(_state.radius >= 0)
        {
          auto oldstate = _state;
          add_or_remove_circle<-1>(x, y, _state.radius);
          auto newsn2 = _state.sn2();

          if( oldsn2 >= _target_sn2 and newsn2 < _target_sn2 )
            {
              // was fine, last time, so reset
              _state = oldstate;
              break;
            }

          _state.radius--;
          oldsn2 = newsn2;
        }
    }

  out_image(x, y) = _state.pixval();

  _last_x = x;
  _last_y = y;
}

void Smoother::smooth_all()
{
  int x = 0;
  int y = 0;
  int xdir = +1;

  const int ydelt = _yw / 10;

  // this is a lambda function to show y value every few steps
  auto showy = [&y, &ydelt]()
    {
          if(y % ydelt == 0)
            {
              std::cout << y/ydelt << ' ';
              std::cout.flush();
            }
    };

  while(y<_yw)
    {
      if(_mask_image(x, y))
        new_pixel(x, y);

      x += xdir;
      if(x == -1)
        {
          xdir = +1;
          x++;
          y++;
          showy();
        }
      if(x == _xw)
        {
          xdir = -1;
          x--;
          y++;
          showy();
        }
    }

  std::cout << '\n';
}

////////////////////////////////////////////////////////////////////////////

// accumulate smoothing

int main(int argc, char* argv[])
{
  double sn = 15;
  std::string back_file, mask_file, expmap_file;
  std::string out_file = "expsmooth.fits";

  parammm::param params(argc, argv);
  params.add_switch( parammm::pswitch( "bg", 'b',
				       parammm::pstring_opt(&back_file),
				       "set background file",
				       "FILE"));
  params.add_switch( parammm::pswitch( "mask", 'm',
				       parammm::pstring_opt(&mask_file),
				       "set mask file",
				       "FILE"));
  params.add_switch( parammm::pswitch( "expmap", 'e',
				       parammm::pstring_opt(&expmap_file),
				       "set exposure map file",
				       "FILE"));
  params.add_switch( parammm::pswitch( "out", 'o',
				       parammm::pstring_opt(&out_file),
				       "set output file (def expsmooth.fits)",
				       "FILE"));
  params.add_switch( parammm::pswitch("sn", 's',
				      parammm::pdouble_opt(&sn),
				      "set signal:noise threshold (def 15)",
				      "VAL"));

  params.set_autohelp("Usage: exposure_smooth [OPTIONS] infile.fits\n"
		      "Accumulative smoothing program with exposure map.\n"
		      "Copyright Jeremy Sanders 2009-2014",
		      "Report bugs to <jsanders@mpe.mpg.de>");
  params.enable_autohelp();
  params.enable_autoversion("0.2",
			    "Jeremy Sanders",
			    "Licenced under the GPL - see the file COPYING");
  params.enable_at_expansion();
  params.interpret_and_catch();

  if(params.args().size() != 1)
    {
      params.show_autohelp();
    }

  const std::string in_filename = params.args()[0];

  // load fg image
  double in_exposure = 1.;
  image_float* in_image = 0;
  load_image(in_filename, &in_exposure, &in_image);

  // load bg image
  double bg_exposure = 1.;
  image_float* bg_image = 0;

  if( back_file.empty() )
    {
      std::cout << "Using blank background\n";
      bg_image = new image_float(in_image->xw(), in_image->yw(), 0.);
      bg_exposure = in_exposure;
    }
  else
    {
      load_image(back_file, &bg_exposure, &bg_image);
    }

  // load expmap image
  image_float* expmap_image = 0;
  if( expmap_file.empty() )
    {
      std::cout << "Using blank exposure map\n";
      expmap_image = new image_float(in_image->xw(), in_image->yw(), 1.);
     }
  else
    {
      load_image(expmap_file, 0, &expmap_image);
    }

  // load mask
  image_short* mask_image = 0;
  if( mask_file.empty() )
    {
      std::cout << "Using blank mask\n";
      mask_image = new image_short( in_image->xw(), in_image->yw(), 1 );

      if( ! expmap_file.empty() )
	{
	  std::cout << "Using exposure map to create mask\n";
	  for(unsigned x = 0; x != in_image->xw(); ++x)
	    for(unsigned y = 0; y != in_image->yw(); ++y)
	      if( (*expmap_image)(x, y) < 1. )
		(*mask_image)(x, y) = 0;
	}
    }
  else
    {
      load_image( mask_file, 0, &mask_image );
    }


  // actually do the work
  Smoother smoother(*in_image, *bg_image, *expmap_image, *mask_image,
                    in_exposure, bg_exposure, sn);
  smoother.smooth_all();

  // write output image
  write_image(out_file, smoother.out_image);

  // clean up
  delete in_image;
  delete bg_image;
  delete mask_image;
  delete expmap_image;

  return 0;
}
