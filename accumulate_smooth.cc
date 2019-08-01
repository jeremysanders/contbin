// accumulate smoothing
// Jeremy Sanders 2002-2005
// Released under the GNU Public License

#include <iostream>
#include <string>

#include "parammm/parammm.hh"
#include "flux_estimator.hh"
#include "misc.hh"
#include "image_disk_access.hh"

using std::string;
using std::cout;

int main(int argc, char* argv[])
{
  string back_file, mask_file;
  string out_file = "acsmooth.fits";
  double sn = 15;

  parammm::param params(argc, argv);
  params.add_switch( parammm::pswitch( "bg", 'b',
				       parammm::pstring_opt(&back_file),
				       "set background file",
				       "FILE"));
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
  params.set_autohelp("Usage: accumulate_smooth [OPTIONS] file.fits\n"
		      "Accumulate smoothing program.\n"
		      "Written by Jeremy Sanders 2004.",
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

  double in_exposure = 1.;
  image_float* in_image;

  load_image( filename, &in_exposure, &in_image);

  double bg_exposure = 1.;
  image_float* bg_image;

  if( back_file.empty() )
    {
      bg_image = 0;
    }
  else
    {
      load_image( back_file, &bg_exposure, &bg_image);
    }

  image_short* mask_image;
  if( mask_file.empty() )
    {
      mask_image = new image_short( in_image->xw(), in_image->yw(), 1 );
    }
  else
    {
      load_image( mask_file, 0, &mask_image );
    }

  const image_float fg_exp(in_image->xw(), in_image->yw(), in_exposure);
  const image_float bg_exp(in_image->xw(), in_image->yw(), bg_exposure);

  flux_estimator fe( in_image, bg_image, mask_image,
		     &fg_exp, &bg_exp, 0, sn);
  image_float out = fe();

  write_image(out_file, out);

  return 0;
}
