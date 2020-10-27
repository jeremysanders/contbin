// accumulate smoothing
// Jeremy Sanders 2002-2005
// Released under the GNU Public License

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>

#include "parammm/parammm.hh"
#include "misc.hh"
#include "image_disk_access.hh"

using std::string;
using std::cout;

struct Point
{
  Point(int _x, int _y) : x(_x), y(_y) {}
  int x, y;
};

typedef std::vector<Point> PointVec;
typedef std::vector<PointVec> PointVecVec;

PointVecVec cachePVV()
{
  // cache delta x, y for each r-squared
  const int maxrad = 100;
  PointVecVec pvv;
  for(int y=-maxrad; y<=maxrad; ++y)
    for(int x=-maxrad; x<=maxrad; ++x)
      {
        unsigned r2 = unsigned(x*x + y*y);
        if(r2 <= maxrad*maxrad)
          {
            while(pvv.size() <= r2)
              pvv.push_back(PointVec());
            pvv[r2].push_back(Point(x,y));
          }
      }
  return pvv;
}

void construct_scale(const image_float& inimg, const image_short& maskimg, double sn,
                     image_short& scaleimg)
{
  const PointVecVec pvv(cachePVV());
  const double maxcts = sn*sn;

  // now iterate over each pixel and measure binning scale
  for(unsigned y=0; y<inimg.yw(); ++y)
    {
      if(y%100 == 0)
        std::cout << y << '\n';
      for(unsigned x=0; x<inimg.xw(); ++x)
        {
          if(maskimg(x,y) < 1 && maskimg(x,y) != -2)
            continue;

          double sum = 0;
          unsigned r2 = 0;
          while( r2 < pvv.size() )
            {
              for(auto pt : pvv[r2])
                {
                  int xi = int(x)+pt.x;
                  int yi = int(y)+pt.y;
                  if(xi>=0 && yi>=0 && xi<int(inimg.xw()) && yi<int(inimg.yw()) && maskimg(xi,yi)>0)
                    sum += double(inimg(xi,yi));
                }
              if(sum >= maxcts)
                break;
              r2 += 1;
            }

          scaleimg(x,y) = r2;
        }
    }
}

void apply_scale(const image_float& inimg, const image_short& maskimg,
                 const image_short& scaleimg, image_float& outimg)
{
  const PointVecVec pvv(cachePVV());

  for(unsigned y=0; y<inimg.yw(); ++y)
    {
      if(y%100 == 0)
        std::cout << y << '\n';
      for(unsigned x=0; x<inimg.xw(); ++x)
        {
          if(maskimg(x,y) < 1 && maskimg(x,y) != -2)
            continue;

          double sum = 0;
          unsigned npix = 0;
          for(int r2=0; r2<=scaleimg(x,y) && r2<int(pvv.size()); ++r2)
            {
              for(auto pt : pvv[r2])
                {
                  int xi = int(x)+pt.x;
                  int yi = int(y)+pt.y;
                  if(xi>=0 && yi>=0 && xi<int(inimg.xw()) && yi<int(inimg.yw()) && maskimg(xi,yi)>0)
                    {
                      ++npix;
                      sum += double(inimg(xi,yi));
                    }
                }
            }
          outimg(x,y) = sum / npix;
        }
    }
}

#define MAXEXP 12.0
#define EXPNSTEPS 1024
#define STEPSIZE (MAXEXP/EXPNSTEPS)
void make_exp_cache(std::vector<float>& cache)
{
  for(size_t i=0; i<EXPNSTEPS; ++i)
    {
      double v = i*-STEPSIZE;
      cache.push_back(std::exp(v));
    }
}

inline float quick_exp(const std::vector<float>& cache, float val)
{
  float fidx = val * float(-1/STEPSIZE);
  int iidx = int(fidx);
  if(iidx<0 or iidx+1>=int(cache.size()))
    return 0;
  return (fidx-iidx)*cache[iidx+1] + (1+iidx-fidx)*cache[iidx];
}

void apply_scale_gaussian(const image_float& inimg, const image_short& maskimg,
                          const image_short& scaleimg, image_float& outimg)
{
  std::vector<float> expcache;
  make_exp_cache(expcache);

  for(unsigned y=0; y<inimg.yw(); ++y)
    {
      if(y%10 == 0)
        std::cout << y << '\n';
      for(unsigned x=0; x<inimg.xw(); ++x)
        {
          if((maskimg(x,y)<1 && maskimg(x,y)!=-2) || scaleimg(x,y)<0)
            continue;

          float sum = 0;
          float sum_weights = 0;
          float sigma = std::max(1.f, std::sqrt(float(scaleimg(x,y))));
          float nh_invsigma2 = -0.5f/(sigma*sigma);

          int rng = int(sigma*4);
          for(int dy=-rng; dy<=rng; ++dy)
            for(int dx=-rng; dx<=rng; ++dx)
              {
                int nx = x+dx;
                int ny = y+dy;
                if( nx>=0 && ny>=0 && nx<int(inimg.xw()) && ny<int(inimg.yw()) &&
                    maskimg(nx,ny)>0 )
                  {
                    int rad2 = dx*dx+dy*dy;
                    //float weight = std::exp(nh_invsigma2*rad2);
                    float weight = quick_exp(expcache, nh_invsigma2*rad2);
                    sum_weights += weight;
                    sum += weight*inimg(nx,ny);
                  }
              }

          outimg(x,y) = sum / sum_weights;
        }
    }

}

int main(int argc, char* argv[])
{
  string back_file, mask_file;
  string scale_file = "acscale.fits";
  string app_file = "applied.fits";
  double sn = 15;
  bool apply_mode = false;
  bool apply_gaussian = false;

  parammm::param params(argc, argv);
  params.add_switch( parammm::pswitch( "apply", 'a',
                                       parammm::pbool_noopt(&apply_mode),
                                       "apply scales to existing data",
                                       ""));
  params.add_switch( parammm::pswitch( "gaussian", 'g',
                                       parammm::pbool_noopt(&apply_gaussian),
                                       "apply scales in gaussian mode",
                                       ""));
  params.add_switch( parammm::pswitch( "mask", 'm',
				       parammm::pstring_opt(&mask_file),
				       "set mask file",
				       "FILE"));
  params.add_switch( parammm::pswitch( "applied", 'a',
				       parammm::pstring_opt(&app_file),
				       "set output file (def out.fits)",
				       "FILE"));
  params.add_switch( parammm::pswitch( "scale", 'c',
				       parammm::pstring_opt(&scale_file),
				       "set output file (def acscale.fits)",
				       "FILE"));
  params.add_switch( parammm::pswitch("sn", 's',
				      parammm::pdouble_opt(&sn),
				      "set signal:noise threshold (def 15)",
				      "VAL"));
  params.set_autohelp("Usage: accumulate_counts [OPTIONS] file.fits\n"
		      "Measure smoothing scale from count data, to be applied later to other data.\n"
		      "Written by Jeremy Sanders 2020.",
		      "Report bugs to <jeremy@jeremysanders.net>");
  params.enable_autohelp();
  params.enable_autoversion("0.1",
			    "Jeremy Sanders",
			    "Licenced under the GPL - see the file COPYING");
  params.enable_at_expansion();
  params.interpret_and_catch();

  if(params.args().size() != 1)
    {
      params.show_autohelp();
    }

  const string filename = params.args()[0];
  image_float* in_image;

  load_image( filename, 0, &in_image);

  image_short* mask_image;
  if( mask_file.empty() )
    {
      mask_image = new image_short( in_image->xw(), in_image->yw(), 1 );
    }
  else
    {
      load_image( mask_file, 0, &mask_image );
    }

  if( ! apply_mode )
    {
      image_short* scale_img = new image_short(in_image->xw(), in_image->yw(), -1);

      construct_scale(*in_image, *mask_image, sn, *scale_img);

      write_image(scale_file, *scale_img);
    }
  else
    {
      image_float* out_img = new image_float(in_image->xw(), in_image->yw(),
                                             std::numeric_limits<float>::quiet_NaN());
      image_short* scale_img;

      load_image( scale_file, 0, &scale_img );

      if(apply_gaussian)
        apply_scale_gaussian(*in_image, *mask_image, *scale_img, *out_img);
      else
        apply_scale(*in_image, *mask_image, *scale_img, *out_img);

      write_image(app_file, *out_img);
    }

  return 0;
}
