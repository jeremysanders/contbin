// code takes a binmap and extracts a set of region files

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>

#include "parammm/parammm.hh"

#include "misc.hh"
#include "fitsio_simple.hh"

using namespace std;

const char* c_prog_version = "0.4";

class extractor
{
public:
  extractor(const string& binmapfile, const string& outdir);
  ~extractor();

  void set_binning(double minx, double miny, double binning);
  void extract();

private:
  // get region file for bin, returns false if error
  bool extract_bin(int no, const string& outfname);

private:
  double _minx, _miny, _binning;

  image_long* _binmap;
  string _outdir;
};

extractor::extractor(const string& binmapfile, const string& outdir)
  : _minx(0), _miny(0), _binning(1),
    _outdir(outdir)
{
  FITSFile file(binmapfile);

  file.readImage(&_binmap);
}

extractor::~extractor()
{
  delete _binmap;
}

void extractor::set_binning(double minx, double miny, double binning)
{
  _minx = minx; _miny = miny; _binning = binning;
}

void extractor::extract()
{
  int no = 0;
  for(;;) {

    ostringstream o;
    o << _outdir << "/xaf_" << no << ".reg";

    const string fname = o.str();

    cout << "Bin " << no << " (" << fname << ")\n";

    if( ! extract_bin(no, fname) )
      break;

    no++;
  }

}

bool extractor::extract_bin(int no, const string& outfname)
{
  const int ixw = _binmap->xw();
  const int iyw = _binmap->yw();

  int found = false;
  dm::memimage<short> inbinimage(ixw, iyw);
  for(int y=0; y<iyw; ++y)
    for(int x=0; x<ixw; ++x) 
      if( (*_binmap)(x, y) == no )
	{
	  inbinimage(x, y) = 1;
	  found = true;
	}

  if(!found) return false;

  ofstream file( outfname.c_str() );
  if( ! file ) {
    cerr << "Unable writing to file " << outfname << endl;
    exit(1);
  }

  file << "# Region file format: CIAO version 1.0\n";

  // single pass rectangular bin identification
  for(int y=0; y<iyw; ++y)
    for(int x=0; x<ixw; ++x) {

      // ignore if bin isn't set
      if( not inbinimage(x, y) )
	continue;

      // try and extend pixel...
      int xw=1, yw=1;

      for(;;)
	{
	  // check for y extend
	  bool yextend = true;
	  for(int xi = 0; xi < xw; xi++)
	    {
	      if( (y+yw) >= iyw || not inbinimage(x+xi, y+yw) )
		{
		  yextend = false;
		  break;
		}
	    }
	  // check for x extend
	  bool xextend = true;
	  for(int yi = 0; yi < yw; yi++)
	    {
	      if( (x+xw) >= ixw || not inbinimage(x+xw, y+yi) )
		{
		  xextend = false;
		  break;
		}
	    }
	
	  bool xandyextend = false;
	  if( (x+xw) < ixw && (y+yw) < iyw )
	    xandyextend = inbinimage(x+xw, y+yw);

	  if( ! xextend && ! yextend ) break;
	  if( xextend && yextend && !xandyextend ) xw++;
	  else
	    {
	      if( xextend ) xw++;
	      if( yextend ) yw++;
	    }
	} // extend loop

      // zero pixels in "bin"
      for(int yi = 0; yi < yw; yi++)
	for(int xi = 0; xi < xw; xi++)
	  inbinimage(x+xi, y+yi) = 0;

      file << "rotbox("
	   << _minx + _binning*(x + xw*0.5) << ','
	   << _miny + _binning*(y + yw*0.5) << ','
	   << _binning*xw << ','
	   << _binning*yw << ','
	   << "0)\n" ;
    }


  if( ! file ) {
    cerr << "Error writing to " << outfname << endl;
    exit(1);
  }
  return true;
}

int main(int argc, char *argv[])
{
  parammm::param params(argc, argv);

  params.set_autohelp("Usage: xaf_make_region_files [OPTION] "
		      "--minx=val --miny=val --bin=val "
		      "--outdir=dir/ binmap.fits\n"
		      "Create region files from binmap.\n"
		      "Written by Jeremy Sanders, 2002.",
		      "Report bugs to <jss@ast.cam.ac.uk>");

  double minx = 0, miny = 0, bin = 1;
  string outdir = ".";

  params.add_switch( parammm::pswitch("minx", 'x',
				      parammm::pdouble_opt(&minx),
				      "Set minimum sky x coord for "
				      "bin map (req)", "PIX") );
  params.add_switch( parammm::pswitch("miny", 'y',
				      parammm::pdouble_opt(&miny),
				      "Set minimum sky y coord for "
				      "bin map (req)", "PIX") );
  params.add_switch( parammm::pswitch("bin", 'b',
				      parammm::pdouble_opt(&bin),
				      "Set sky binning factor for "
				      "bin map (def 1)", "PIX") );
  params.add_switch( parammm::pswitch("outdir", 'o',
				      parammm::pstring_opt(&outdir),
				      "Set output directory (def .)",
				      "DIR") );

  params.enable_autohelp();
  params.enable_at_expansion();
  params.enable_autoversion(c_prog_version, "Jeremy Sanders",
			    "Released under the GPL");
  params.interpret_and_catch();

  if(params.args().size() != 1)
    params.show_autohelp();

  cout << "Input binmap: " << params.args()[0] << endl
       << "Output directory: " << outdir << endl
       << "Minimum x: " << minx << endl
       << "Minimum y: " << miny << endl
       << "Binning factor: " << bin << endl;

  extractor e(params.args()[0], outdir);
  e.set_binning(minx, miny, bin);
  e.extract();

  return 0;
}
