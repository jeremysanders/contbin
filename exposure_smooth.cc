#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>

#include "parammm/parammm.hh"
#include "fitsio_simple.hh"
#include "misc.hh"

#include "image_disk_access.hh"

// this is a program to accumulatively smooth an X-ray image
// with an optional background image and exposure map image

// it assumes that the background image has the same relative exposure map
// as the fg image, but uses the EXPOSURE key to calculate the fg/bg exposure
// time ratio

using std::string;
using std::cout;
using std::cerr;
using std::vector;
using std::sqrt;
using std::max;

// basic stuct to store integer x, y points
struct point
{
  point(int _x, int _y) : x(_x), y(_y) {};

  int x, y;
};

typedef vector<point> point_vec;
typedef vector<point_vec> point_vec_vec;

// get the points as a function of integer radius
point_vec_vec collectRadii(const int maxrad)
{
  point_vec_vec retn(maxrad+1);

  for(int y = -maxrad; y <= maxrad; ++y)
    for(int x = -maxrad; x <= maxrad; ++x)
      {
	int r = int( sqrt(x*x+y*y) );
        if(r<=maxrad)
          retn[r].push_back( point(x, y) );
      }
  return retn;
}

// square value
template<class T> T sqd(T v)
{
  return v*v;
}

// calculate signal to noise ratio squared given fg and bg counts
// and exposure times fgtime and bgtime (note inv=1/time)
float SNratio2(float fg, float bg, float invfgtime, float invbgtime)
{
  if(fg == 0 && bg == 0)
    return 0;

  return sqd(fg*invfgtime - bg*invbgtime) /
    (fg*sqd(invfgtime) + bg*sqd(invbgtime));
}

// do the actual smoothing
void smoothImage(const image_float& inimage, const image_float& bgimage,
		 const image_float& expmapimage,
		 float sn, int maxrad,
                 float exptimefg, float exptimebg,
		 image_float& outimage)
{
  const int xw = inimage.xw();
  const int yw = inimage.yw();

  if(maxrad <= 0)
    // use diagonal of image as maximum radius if not specified
    maxrad = int(std::sqrt(xw*xw+yw*yw))+1;

  point_vec_vec ptsatradii(collectRadii(maxrad));

  const float invexptimefg = 1/exptimefg;
  const float invexptimebg = 1/exptimebg;
  const float sn2 = sqd(sn);

  for(int y=0; y<yw; ++y)
    {
      if(y%20==0)
        cout << "y=" << y << '/' << yw << '\n';

      for(int x=0; x<xw; ++x)
        {
          if(expmapimage(x, y) <= 0)
            continue;

          float totalfg = 0;
          float totalbg = 0;
          float totalexp = 0;

          for(int radius = 0;
              (SNratio2(totalfg, totalbg, invexptimefg, invexptimebg)<sn2) &&
                (radius<=maxrad);
               ++radius )
            {
              for(auto const& p : ptsatradii[radius])
                {
                  const int nx = x + p.x;
                  const int ny = y + p.y;

                  if(nx >= 0 && ny >= 0 && nx < xw && ny < yw)
                    {
                      float expos = expmapimage(nx, ny);
                      if(expos > 0)
                        {
                          totalexp += expos;
                          totalfg += inimage(nx, ny);
                          totalbg += bgimage(nx, ny);
                        }
                    }
                }
            }
          outimage(x, y) = (totalfg - totalbg * exptimefg / exptimebg) / totalexp;
        }
    }

  cout << '\n';
}

int main(int argc, char* argv[])
{
  double sn = 15;
  int maxrad = -1;
  string back_file, mask_file, expmap_file;
  string out_file = "expsmooth.fits";

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
  params.add_switch( parammm::pswitch("maxrad", 'r',
				      parammm::pint_opt(&maxrad),
				      "maximum radius (def -1 or infinite)",
				      "VAL"));

  params.set_autohelp("Usage: exposure_smooth [OPTIONS] infile.fits\n"
		      "Accumulative smoothing program with exposure map.\n"
		      "Copyright Jeremy Sanders 2009-2018",
		      "Report bugs to <jeremy@jeremysanders.net>");
  params.enable_autohelp();
  params.enable_autoversion("0.3",
			    "Jeremy Sanders",
			    "Licenced under the GPL - see the file COPYING");
  params.enable_at_expansion();
  params.interpret_and_catch();

  if(params.args().size() != 1)
    {
      params.show_autohelp();
    }

  const string in_filename = params.args()[0];

  // load fg image
  double in_exposure = 1.;
  image_float* in_image = 0;
  load_image(in_filename, &in_exposure, &in_image);

  // load bg image
  double bg_exposure = 1.;
  image_float* bg_image = 0;

  if( back_file.empty() )
    {
      cout << "Using blank background\n";
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
      cout << "Using blank exposure map\n";
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
      cout << "Using blank mask\n";
    }
  else
    {
      load_image( mask_file, 0, &mask_image );
    }

  // check image dimensions
  if(in_image->xw() != expmap_image->xw() || in_image->yw() != expmap_image->yw() ||
     (mask_image != 0 && (in_image->xw() != mask_image->xw() ||
                          in_image->yw() != mask_image->yw())))
    {
      cerr << "Input images have different dimensions\n";
      return 1;
    }

  // set exposure map to zero where mask is zero
  if(mask_image != 0)
    for(unsigned y=0; y<in_image->yw(); ++y)
      for(unsigned x=0; x<in_image->xw(); ++x)
        if((*mask_image)(x, y) == 0)
          (*expmap_image)(x, y) = 0;

  // make a new image full of NaNs to use as output
  image_float* out_image = new image_float(in_image->xw(), in_image->yw(),
					   std::numeric_limits<float>::quiet_NaN());

  // actually do the work
  smoothImage(*in_image, *bg_image, *expmap_image,
	      sn, maxrad, in_exposure, bg_exposure, *out_image);

  // write output image
  write_image(out_file, *out_image);

  // clean up
  delete in_image;
  delete bg_image;
  delete mask_image;
  delete expmap_image;
  delete out_image;
}
