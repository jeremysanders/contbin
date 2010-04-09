// Contour binning program
// See Sanders (2005)

// Copyright 2002-2005 Jeremy Sanders
// Released under the GNU Public Licence (GPL)

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

#include <cmath>
#include <cassert>
#include <cstdlib>

#include "parammm/parammm.hh"

#include "binner.hh"
#include "flux_estimator.hh"
#include "fitsio_simple.hh"

const char* const CONTBIN_VERSION = "1.3";

using std::string;
using std::vector;
using std::find;
using std::ostringstream;
using std::cout;

class program
{
public:
  program(int argc, char **argv);

  void run();

private:
  void auto_mask(const image_float& in_data, image_short* mask);
  
  template<class T> void save_image(const string& filename, const T& image,
				    FITSFile* indataset);
  template<class T> void load_image(const string& filename, double *exposure,
				    T** image);

private:
  string _out_fname, _sn_fname, _binmap_fname;
  string _bg_fname;
  string _in_fname;
  string _mask_fname;
  string _smoothed_fname;
  string _expmap_fname, _bg_expmap_fname;
  string _noisemap_fname;
  double _sn_threshold;
  double _smooth_sn;
  bool _do_automask;
  bool _constrain_fill;
  double _constrain_val;
  bool _noscrub;
  bool _binup;
  double _scrub_large;
};

program::program(int argc, char **argv)
  : _out_fname("contbin_out.fits"),
    _sn_fname("contbin_sn.fits"),
    _binmap_fname("contbin_binmap.fits"),
    _sn_threshold(15.),
    _smooth_sn(15.),
    _do_automask(false),
    _constrain_fill(false),
    _constrain_val(3),
    _noscrub(false),
    _binup(false),
    _scrub_large(-1)
{
  parammm::param params(argc, argv);
  params.add_switch( parammm::pswitch("out", 'o',
				      parammm::pstring_opt(&_out_fname),
				      "set out file (def contbin_out.fits)",
				      "FILE"));
  params.add_switch( parammm::pswitch("outsn", 'e',
				      parammm::pstring_opt(&_sn_fname),
				      "set signal:noise out file (def contbin_sn.fits)",
				      "FILE"));
  params.add_switch( parammm::pswitch("outbinmap", 'n',
				      parammm::pstring_opt(&_binmap_fname),
				      "set binmap out file (def contbin_binmap.fits)",
				      "FILE"));
  params.add_switch( parammm::pswitch("bg", 'b',
				      parammm::pstring_opt(&_bg_fname),
				      "Set background image file (def none)",
				      "FILE"));
  params.add_switch( parammm::pswitch("mask", 'm',
				      parammm::pstring_opt(&_mask_fname),
				      "Set mask image file (def none)",
				      "FILE"));
  params.add_switch( parammm::pswitch("smoothed", 'm',
				      parammm::pstring_opt(&_smoothed_fname),
				      "Set smoothed image file (def none)",
				      "FILE"));
  params.add_switch( parammm::pswitch("expmap", 0,
				      parammm::pstring_opt(&_expmap_fname),
				      "Set exposure map (fg) (def none)",
				      "FILE"));
  params.add_switch( parammm::pswitch("bgexpmap", 0,
				      parammm::pstring_opt(&_bg_expmap_fname),
				      "Set exposure map (bg) (def none)",
				      "FILE"));
  params.add_switch( parammm::pswitch("noisemap", 0,
				      parammm::pstring_opt(&_noisemap_fname),
				      "Set noise map (def none)",
				      "FILE"));

  params.add_switch( parammm::pswitch("sn", 's',
				      parammm::pdouble_opt(&_sn_threshold),
				      "set signal:noise threshold (def 15)",
				      "VAL"));

  params.add_switch( parammm::pswitch("automask", 0,
				      parammm::pbool_noopt(&_do_automask),
				      "automatically mask image",
				      ""));

  params.add_switch( parammm::pswitch("constrainfill", 0,
				      parammm::pbool_noopt(&_constrain_fill),
				      "constrain filling-factor",
				      ""));
  params.add_switch( parammm::pswitch("constrainval", 0,
				      parammm::pdouble_opt(&_constrain_val),
				      "set constrain ratio (def 3)",
				      "VAL"));
  params.add_switch( parammm::pswitch("smoothsn", 0,
				      parammm::pdouble_opt(&_smooth_sn),
				      "set smoothing signal:noise (def 15)",
				      "VAL"));

  params.add_switch( parammm::pswitch("noscrub", 0,
				      parammm::pbool_noopt(&_noscrub),
				      "don't scrub low S/N bins",
				      ""));
  params.add_switch( parammm::pswitch("binup", 0,
				      parammm::pbool_noopt(&_binup),
				      "start binning from bottom",
				      ""));

  params.add_switch( parammm::pswitch("scrublarge", 0,
				      parammm::pdouble_opt(&_scrub_large),
				      "Scrub bins with area frac > this",
				      "VAL"));

  params.set_autohelp("Usage: contbin [OPTIONS] file.fits\n"
		      "Contour binning program\n"
		      "Written by Jeremy Sanders 2002-2010",
		      "Report bugs to <jss@ast.cam.ac.uk>");
  params.enable_autohelp();
  params.enable_autoversion(CONTBIN_VERSION,
			    "Jeremy Sanders",
			    "Licenced under the GPL - see the file COPYING");
  params.enable_at_expansion();

  params.interpret_and_catch();

  if(params.args().size() != 1)
    {
      params.show_autohelp();
    }
  else
    {
      _in_fname = params.args()[0];
    }
}

void program::auto_mask(const image_float& in_data, image_short* mask)
{
  cout << "(i) Automasking... ";
  cout.flush();

  const unsigned blocksize = 8;

  const unsigned xw = in_data.xw();
  const unsigned yw = in_data.yw();
  assert( mask->xw() == xw && mask->yw() == yw );

  const unsigned xw_bl = xw / blocksize + 1;
  const unsigned yw_bl = yw / blocksize + 1;

  // default, we see all pixels
  mask->set_all(1);

  for(unsigned yb = 0; yb < yw_bl; ++yb)
    for(unsigned xb = 0; xb < xw_bl; ++xb)
      {
	const unsigned sx = xb*blocksize;
	const unsigned sy = yb*blocksize;
	double sum = 0.;

	// count counts in blocksize x blocksize region
	for(unsigned dy=0; dy<blocksize; dy++)
	  for(unsigned dx=0; dx<blocksize; dx++)
	    {
	      const unsigned x = sx+dx, y = sy+dy;
	      if( x >= xw || y >= yw)
		continue;

	      sum += in_data(x, y);
	    }

	// no counts, therefore throw away
	if( fabs(sum) < 1e-5 )
	  {
	    for(unsigned dy=0; dy<blocksize; dy++)
	      for(unsigned dx=0; dx<blocksize; dx++)
		{
		  const unsigned x = sx+dx, y = sy+dy;
		  if( x >= xw || y >= yw)
		    continue;

		  mask->pixel(x, y) = 0;
		}
	  }
      }

  cout << "Done\n";
}

// handy template to load an image
template<class T> void program::load_image(const string& filename, double *exposure,
					   T** image)
{
  FITSFile dataset(filename);
  
  if(exposure != 0)
    {
      double defval = 1.;
      dataset.readKey("EXPOSURE", exposure, &defval);
    }
  
  dataset.readImage(image);
}

// write an image with some sensible headers
// copies keywords from indataset
template<class T> void program::save_image(const string& filename, const T& image,
					   FITSFile* indataset)
{
  FITSFile dataset(filename, FITSFile::Create);
  dataset.writeImage(image);

  indataset->copyHeaderTo(dataset);

  dataset.writeDatestamp("contbin");

  ostringstream o;
  o << "Generated by contbin (Jeremy Sanders 2005)\n"
    << "This filename: " << filename << '\n'
    << "Input image: " << _in_fname << '\n'
    << "Back image: " << _bg_fname << '\n'
    << "Mask image: " << _mask_fname << '\n'
    << "Smoothed image: " << _smoothed_fname << '\n'
    << "Expmap image: " << _expmap_fname << '\n'
    << "Back expmap image: " << _bg_expmap_fname << '\n'
    << "Noise map image: " << _noisemap_fname << '\n'
    << "SN threshold: " << _sn_threshold << '\n'
    << "Smooth SN: " << _smooth_sn << '\n'
    << "Automask: " << _do_automask << '\n'
    << "Constrain fill: " << _constrain_fill << '\n'
    << "Constrain val: " << _constrain_val << '\n'
    << "No scrub: " << _noscrub << '\n'
    << "Bin up: " << _binup << '\n'
    << "Scrub large: " << _scrub_large << '\n';

  // split output text and write as lines of history
  vector<string> items = split_string(o.str(), '\n');
  for(vector<string>::const_iterator i = items.begin();
      i != items.end(); ++i)
    {
      dataset.writeHistory(*i);
    }
}

void program::run()
{
  /////////////////////////////////////////////////////////////////
  // load input images

  image_float* in_image;
  double in_exposure;

  // load in the main dataset
  // (keep open to copy header)
  cout << "(i) Loading image " << _in_fname << '\n';

  FITSFile indataset(_in_fname);
  double defval = 1.;
  indataset.readKey("EXPOSURE", &in_exposure, &defval);
  indataset.readImage(&in_image);

  // do automasking (if any)
  image_short mask(in_image->xw(), in_image->yw(), 1);
  if( _do_automask )
    auto_mask(*in_image, &mask);

  // load mask (if any)
  if( ! _mask_fname.empty() )
    {
      cout << "(i) Loading masking image " << _mask_fname << '\n';

      image_short* maskim;
      load_image( _mask_fname, 0, &maskim );
      mask = *maskim;
      delete maskim;
    }

  image_float* expmap;
  if( ! _expmap_fname.empty() )
    {
      cout << "(i) Loading foreground exposure map "
	   << _expmap_fname << '\n';
      load_image( _expmap_fname, 0, &expmap );

      // mask out pixels with no exposure
      for(unsigned y = 0; y != expmap->yw(); ++y)
	for(unsigned x = 0; x != expmap->xw(); ++x)
	  {
	    if( (*expmap)(x, y) < 1. )
	      {
		mask(x, y) = 0;
	      }
	  }
    }
  else
    {
      cout << "(i) Using blank foreground exposure (exp="
	   << in_exposure << ")\n";
      expmap = new image_float( in_image->xw(), in_image->yw(),
				in_exposure );
    }

  image_float* bg_image;
  double bg_exposure;
      
  // load in background file (if any)
  if( ! _bg_fname.empty() )
    {
      cout << "(i) Loading background image " << _bg_fname << '\n';
      load_image( _bg_fname, &bg_exposure, &bg_image );
    }
  else
    {
      bg_image = 0;
      bg_exposure = 1;
    }

  // load in background exposure map
  image_float* bg_expmap;
  if( ! _bg_expmap_fname.empty() )
    {
      cout << "(i) Loading background exposure map "
	   << _bg_expmap_fname << '\n';
      load_image( _bg_expmap_fname, 0, &bg_expmap );
    }
  else
    {
      cout << "(i) Using blank background exposure (exp="
	   << bg_exposure << ")\n";
      bg_expmap = new image_float( in_image->xw(), in_image->yw(),
				   bg_exposure );
    }

  // we avoid division by zero by doing this hack
  bg_expmap->trim_up(1e-7);
  expmap->trim_up(1e-7);

  // load in noise map if passed
  image_float* noisemap = 0;
  if( ! _noisemap_fname.empty() )
    {
      cout << "(i) Loading noise map " << _noisemap_fname << '\n';
      load_image( _noisemap_fname, 0, &noisemap);
    }

  // smooth data, or use passed file
  image_float* smoothed_image;
  if( _smoothed_fname.empty() )
    {
      cout << "(i) Smoothing data (S/N = "
	   << _smooth_sn << ")\n";
      // smooth data
      flux_estimator fe( in_image, bg_image, &mask, expmap,
			 bg_expmap, noisemap, _smooth_sn );
      smoothed_image = new image_float( fe() );
    }
  else
    {
      cout << "(i) Loading smoothed image " << _smoothed_fname
	   << '\n';
      load_image( _smoothed_fname, 0, &smoothed_image );
    }

  //////////////////////////////////////////////////////////////////
  // actually do the binning
  binner the_binner(in_image, smoothed_image, _sn_threshold);
  the_binner.set_back_image( bg_image, expmap, bg_expmap );
  the_binner.set_noisemap_image(noisemap);
  the_binner.set_mask_image(&mask);
  the_binner.set_constrain_fill(_constrain_fill, _constrain_val);
  the_binner.set_scrub_large_bins(_scrub_large);

  the_binner.do_binning(!_binup);
  if( ! _noscrub )
    the_binner.do_scrub();
  the_binner.calc_outputs();

  ///////////////////////////////////////////////////////////////////
  // write output images
  save_image(_out_fname, the_binner.get_output_image(), &indataset);
  save_image(_sn_fname, the_binner.get_sn_image(), &indataset);
  save_image(_binmap_fname, the_binner.get_binmap_image(), &indataset);
  save_image("contbin_mask.fits", mask, &indataset);

  delete in_image;
  delete bg_image;
}

int main(int argc, char *argv[])
{
  program prog(argc, argv);
  prog.run();

  return 0;
}
