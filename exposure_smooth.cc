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
void collectRadii(const int maxrad, point_vec_vec& retn)
{
  retn.clear();
  retn.resize( int( M_SQRT2*maxrad )+1 );

  for( int x = -maxrad; x <= maxrad; ++x )
    for( int y = -maxrad; y <= maxrad; ++y )
      {
	const size_t r = size_t( sqrt(x*x+y*y) );
	retn[r].push_back( point(x, y) );
      }
}

// calculate signal to noise ratio given fg and bg counts
// and exposure times fgtime and bgtime
double SNratio(double fg, double bg, double fgtime, double bgtime)
{
  if( fg == 0. and bg == 0. )
    return 0;

  return (fg / fgtime - bg / bgtime) / sqrt( fg / (fgtime*fgtime) +
					     bg / (bgtime*bgtime) );

}

// do the actual smoothing
void smoothImage(const image_float& inimage, const image_float& bgimage,
		 const image_float& expmapimage,
		 const image_short& maskimage,
		 double sn, double exptimefg, double exptimebg,
		 image_float& outimage)
{
  point_vec_vec ptsatradii;
  collectRadii( max(inimage.xw(), inimage.yw()), ptsatradii);

  const int xw = inimage.xw();
  const int yw = inimage.yw();

  for(int y = 0; y < yw; ++y )
    for(int x = 0; x < xw; ++x )
      {
	if( x == 0 and y % (yw/10) == 0 )
	  {
	    cout << y / (yw/10) << ' ';
	    cout.flush();
	  }

	if( ! maskimage(x, y) )
	  continue;

	double totalfg = 0.;
	double totalbg = 0.;
	double totalexp = 0.;
	int pixels = 0;
	
	for( size_t radius = 0;
	     SNratio(totalfg, totalbg, exptimefg, exptimebg) < sn;
	     ++radius )
	  {
	    const point_vec& pts = ptsatradii[radius];
	    for(point_vec::const_iterator p = pts.begin(); p != pts.end(); ++p)
	      {
		const int nx = x + p->x;
		const int ny = y + p->y;
		
		if( nx < 0 or ny < 0 or nx >= xw or ny >= yw or !maskimage(nx, ny) )
		  {
		    continue;
		  }		    

		totalfg += inimage(nx, ny);
		totalbg += bgimage(nx, ny);
		totalexp += expmapimage(nx, ny);
		pixels++;
	      }
	  }

	outimage(x, y) = (totalfg / exptimefg - totalbg / exptimebg) / totalexp;
      }

  cout << '\n';
}

int main(int argc, char* argv[])
{
  double sn = 15;
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

  params.set_autohelp("Usage: exposure_smooth [OPTIONS] infile.fits\n"
		      "Accumulative smoothing program with exposure map.\n"
		      "Copyright Jeremy Sanders 2009",
		      "Report bugs to <jss@ast.cam.ac.uk>");
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

  // load mask
  image_short* mask_image = 0;
  if( mask_file.empty() )
    {
      cout << "Using blank mask\n";
      mask_image = new image_short( in_image->xw(), in_image->yw(), 1 );
    }
  else
    {
      load_image( mask_file, 0, &mask_image );
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

  // make a new image full of NaNs to use as output
  image_float* out_image = new image_float(in_image->xw(), in_image->yw(),
					   std::numeric_limits<float>::quiet_NaN());

  // actually do the work
  smoothImage(*in_image, *bg_image, *expmap_image, *mask_image,
	      sn, in_exposure, bg_exposure, *out_image);

  // write output image
  write_image(out_file, *out_image);

  // clean up
  delete in_image;
  delete bg_image;
  delete mask_image;
  delete expmap_image;
  delete out_image;
}
