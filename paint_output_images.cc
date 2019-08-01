// program takes output files from xaf_fit_spectra
// and creates output images

// Copyright (C) 2002-2006 Jeremy Sanders <jeremy@jeremysanders.net>
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License VERSION 2 as
// published by the Free Software Foundation.  You are not allowed to
// use any other version of the license; unless you have the explicit
// permission from the author to do so.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#include <iostream>
#include <map>
#include <string>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <cmath>
#include <vector>
#include <limits>

#include "parammm/parammm.hh"
#include "fitsio_simple.hh"

#include "misc.hh"
#include "format_string.hh"

// changes
// 1.1 - modified algorithm to read all bins in one go,
//       and output them in one go -
//       much faster
// 1.2 - now uses maps so input params don't have to be in same order
// 2.0 - rewrite to use memimage as output images

using namespace std;

const string c_prog_version = "2.0";

typedef map<unsigned, double> DataMap;
typedef map<string, DataMap> AllDataMap;

struct Bin
{
  Bin(const string& f, int n) : filename(f), number(n) {};

  string filename;
  int number;
};
typedef vector<Bin> BinList;

class painter
{
public:
  painter(int argc, char **argv);

  void run();

private:
  void read_bin_list();
  // get the list of bins to paint from region_list.txt

  void read_variables(unsigned bin, istream& file);
  // paint variables to images and write to disk
  void paint_variables();
  // loop over the bins, processing them
  void iterate_bins();

private:
  string m_binmap_filename;
  string m_input_dir;
  string m_output_dir;
  bool m_gzip;

  AllDataMap m_output_data;
  BinList m_binlist;
};

painter::painter(int argc, char **argv)
  : m_binmap_filename("binmap.fits"),
    m_input_dir("."),
    m_output_dir("."),
    m_gzip(false)
{

  parammm::param params(argc, argv);

  params.add_switch( parammm::pswitch("binmap", 'n',
				      parammm::pstring_opt(&m_binmap_filename),
				      "Set binmap file"
				      " (def. binmap.fits)",
				      "FILE") );
  params.add_switch( parammm::pswitch("input_dir", 'i',
				      parammm::pstring_opt(&m_input_dir),
				      "Input directory containing fit results"
				      " (def. '.')",
				      "DIR") );
  params.add_switch( parammm::pswitch("output_dir", 'o',
				      parammm::pstring_opt(&m_output_dir),
				      "Output directory for FITS images"
				      " (def. '.')",
				      "DIR") );
  params.add_switch( parammm::pswitch("gzip", 0,
				      parammm::pbool_noopt(&m_gzip),
				      "Gzip output files",
				      ""));

  params.set_autohelp("Usage: paint_output_images [OPTION]\n"
		      "Make FITS images from fit results\n"
		      "Written by Jeremy Sanders, 2002-2006.",
		      "Report bugs to <jeremy@jeremysanders.net>");
  params.enable_autohelp();
  params.enable_autoversion(c_prog_version, "Jeremy Sanders",
			    "Released under the GPL");
  params.enable_at_expansion();

  params.interpret_and_catch();
}

// read the variables from file
void painter::read_variables(unsigned bin, istream& file)
{
  while( file )
    {
      string line;
      getline(file, line);
      if( ! file )
	break;
      
      // get parameter and value from line
      istringstream i(line);
      string param;
      string value;
      double val;
      i >> param >> value;
      
      // special case nan
      if( value == "nan" )
	{
	  val = numeric_limits<double>::quiet_NaN();
	} else {
	  istringstream t(value);
	  t >> val;
	}

      if( ! i )
	{
	  cerr << "Error in line '" << line << "'\n";
	  exit(1);
	}
      
      m_output_data[param][bin] = val;
    }

}

void painter::paint_variables()
{
  FITSFile binmapds(m_binmap_filename);
  image_long* binmap_image;
  binmapds.readImage(&binmap_image);

  const unsigned xw = binmap_image->xw();
  const unsigned yw = binmap_image->yw();

  image_dbl im(xw, yw);

  cout << "Painting output...\n";

  // iterate over fit parameters
  for( AllDataMap::const_iterator param = m_output_data.begin();
       param != m_output_data.end(); ++param )
    {
      cout << " Parameter " << param->first << '\n';
      const DataMap& data = param->second;

      // blank image (filled with NaN)
      im.set_all(numeric_limits<double>::quiet_NaN());

      // iterate over pixels
      for(unsigned y=0; y<yw; ++y)
	for(unsigned x=0; x<xw; ++x)
	  {
	    // get the bin
	    const long bin = (*binmap_image)(x, y);
	    if( bin < 0 )
	      continue;

	    // paint the fit parameter
	    const DataMap::const_iterator val = data.find(bin);
	    if( val != data.end() )
	      {
		im(x, y) = val->second;
	      }
	  }

      // write output image
      Util::FormatString filename("%1/%2_out.fits%3");
      filename << m_output_dir << ends << param->first << ends;
      filename << (m_gzip ? ".gz" : "") << ends;

      FITSFile file(filename.get(), FITSFile::Create);
      file.writeImage(im);
      binmapds.copyHeaderTo(file);
      file.writeDatestamp("paint_output_images");

      // write some header items
      ostringstream o;
      o << "Generated by paint_output_images (Jeremy Sanders 2006-2014)\n"
	<< "This filename: " << filename.get() << '\n'
	<< "Input binmap: " << m_binmap_filename << '\n'
	<< "Variable: " << param->first << '\n';

      // split output text and write as lines of history
      vector<string> items = split_string(o.str(), '\n');
      for(vector<string>::const_iterator i = items.begin();
	  i != items.end(); ++i)
	{
	  file.writeHistory(*i);
	}
    }

  cout << "Done\n";
}

// get a sequence of numbers from somewhere in a string
int get_number(const string& s)
{
  string nums;
  for(string::const_iterator i = s.begin();
      i != s.end();
      ++i)
    {
      if( *i >= '0' and *i <= '9' )
	{
	  // add the character
	  nums.push_back(*i);
	}
      else
	{
	  // exit if we run out of numbers
	  if( ! nums.empty() )
	    break;
	}
    }

  istringstream i(nums);
  int num;
  i >> num;
  if(!i)
    {
      cerr << "String '" << s << "' does not contain a number\n";
      exit(1);
    }

  return num;
}

void painter::read_bin_list()
{
  Util::FormatString filename("%1/region_list.txt");
  filename << m_input_dir << ends;

  ifstream infile(filename.get().c_str());

  if( ! infile )
    {
      cerr << "Cannot open list of region file "
	   << filename.get() << '\n';
      exit(1);
    }

  while(infile)
    {
      string line;
      getline(infile, line);
      istringstream i(line);

      string name;
      string filename;
      i >> name >> filename;

      if(infile && i)
	{
	  Util::FormatString resfile("%1/%2_fit_out.txt");
	  resfile << m_input_dir << ends
		  << name << ends;
	  const int num = get_number(name);
	  m_binlist.push_back( Bin(resfile.get(), num) );
	}
    }
}

void painter::iterate_bins()
{
  cout << "Reading bin data from files... " << flush;

  unsigned count = 0;
  for( BinList::const_iterator i = m_binlist.begin();
       i != m_binlist.end(); ++i )
    {
      const string& filename = i->filename;
      const int number = i->number;

      ifstream file(filename.c_str());
      if( ! file ) continue;

      cout << "Reading bin " << number << '\n';
      read_variables(number, file);
      ++count;
    }

  cout << count << " bins read\n";
}

void painter::run()
{
  read_bin_list();
  iterate_bins();
  paint_variables();
}

////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
  painter prog(argc, argv);
  prog.run();

  return 0;
}
