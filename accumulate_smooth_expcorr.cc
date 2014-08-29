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

namespace
{

  template <typename T> T sqr(T v) { return v*v; }

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

  class Smoother
  {
  public:
    Smoother(const image_short& ct_image,
             const image_float& expcorr_image,
             const image_short& mask_image,
             double sn)
      : _ct_image(ct_image),
        _expcorr_image(expcorr_image),
        _mask_image(mask_image),
        _xw(ct_image.xw()), _yw(ct_image.yw()),
        _target_sn2(sn*sn),
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
              _tot_ct += VAL*_ct_image(nx, ny);
              _tot_expcorr += VAL*_expcorr_image(nx, ny);
              _tot_pix += VAL;
            }
        }
    }

    // calculate signal to noise squared
    double sn2() const { return _tot_ct; }

  private:
    const image_short& _ct_image;
    const image_float& _expcorr_image;
    const image_short& _mask_image;
    const int _xw;
    const int _yw;
    const double _target_sn2;

    point_vec_vec _circles;
    point_vec_vec _shift_incl;

    int _radius;
    int _tot_ct;
    KahanSum _tot_expcorr;
    int _tot_pix;

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
  _radius = 0;
  _tot_ct = 0;
  _tot_expcorr.reset();
  _tot_pix = 0;

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
          _tot_ct += sign*_ct_image(nx, ny);
          _tot_expcorr += sign*_expcorr_image(nx, ny);
          _tot_pix += sign;
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
      bool iny, mirror;
      iny = _last_y != y;
      mirror = _last_x > x or _last_y > y;

      add_shift(_last_x, _last_y, _radius, -1, iny, not mirror);
      add_shift(x, y, _radius, 1, iny, mirror);

      // move backwards or forwards in radius depending on current S/N
      reverse = sn2() >= _target_sn2;
    }
  else
    {
      reset_state();
      // we have to add the initial pixel, as this is checked below
      add_or_remove_circle<1>(x, y, _radius);
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
      while(sn2() < _target_sn2)
        {
          _radius++;
          add_or_remove_circle<1>(x, y, _radius);
        }
    }
  else
    {
      // Reverse direction.
      // We check whether subtracting this radius
      // pushes us under the threshold.

      while(_radius >= 0)
        {
          auto oldct = _tot_ct;
          auto oldexpcorr = _tot_expcorr;
          auto oldpix = _tot_pix;
          auto oldsn2 = sn2();

          add_or_remove_circle<-1>(x, y, _radius);
          auto newsn2 = sn2();

          if( oldsn2 >= _target_sn2 and newsn2 < _target_sn2 )
            {
              // was fine, last time, so reset
              _tot_ct = oldct;
              _tot_expcorr = oldexpcorr;
              _tot_pix = oldpix;
              break;
            }

          _radius--;
        }
    }

  out_image(x, y) = _tot_expcorr.sum() / _tot_pix;

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

namespace
{
  template<class T> void load_image(const std::string& filename,
                                    T** image)
  {
    FITSFile ds(filename);
    ds.readImage(image);
  }

  void write_image(const std::string& filename, const image_float& img,
                   FITSFile* indataset)
  {
    FITSFile ds(filename, FITSFile::Create);
    ds.writeImage(img);
    indataset->copyHeaderTo(ds);
  }
}

int main(int argc, char* argv[])
{
  std::string back_file, mask_file;
  std::string out_file = "acsmooth.fits";
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
  params.set_autohelp("Usage: accumulate_smooth_expcorr [OPTIONS] img.fits expcorr.fits\n"
		      "Accumulate smoothing program (using exposure corrected image).\n"
		      "Written by Jeremy Sanders 2014.",
		      "Report bugs to <jsanders@mpe.mpg.de>");
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

  image_short* in_image;
  FITSFile indataset(params.args()[0]);
  indataset.readImage(&in_image);

  image_float* expcorr_image;
  load_image(params.args()[1], &expcorr_image);

  image_short* mask_image;
  if( mask_file.empty() )
    {
      mask_image = new image_short(in_image->xw(), in_image->yw(), 1);
    }
  else
    {
      load_image(mask_file, &mask_image);
    }

  Smoother smoother(*in_image, *expcorr_image, *mask_image, sn);
  smoother.smooth_all();

  write_image(out_file, smoother.out_image, &indataset);

  delete in_image;
  delete expcorr_image;
  delete mask_image;

  return 0;
}
