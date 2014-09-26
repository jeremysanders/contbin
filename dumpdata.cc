#include <vector>
#include <string>
#include <cassert>
#include <fstream>
#include <string>
#include <cmath>
#include <iostream>

#include "misc.hh"
#include "fitsio_simple.hh"

class point_holder
{
public:
  //  point_holder();
  point_holder(const image_long& binmap,
	       const image_dbl& vals);

  unsigned no_bins() { return _no_bins; }
  void get_bin_info(unsigned bin,
		    double* ret_x, double* ret_y, double* ret_val,
		    double* ret_rms_dx, double* ret_rms_dy,
		    double* ret_total,
		    unsigned* pix_count);

private:
  unsigned count_no_bins(const image_long& binmap);
  void calc_points(const image_long& binmap,
		   const image_dbl& vals);

private:
  struct bininfo
  {
    double x, y;
    double val;
    double total;
    double rms_dx, rms_dy;
    unsigned pix_count;

    bininfo() { x = y = val = total = rms_dx = rms_dy = 0; pix_count = 0; }
  };
  typedef std::vector<bininfo> vec_bininfo;

private:
  const unsigned _no_bins;
  vec_bininfo _bininfos;
};


point_holder::point_holder(const image_long& binmap,
			   const image_dbl& vals)
  : _no_bins( count_no_bins( binmap ) ),
    _bininfos( _no_bins )

{
  calc_points(binmap, vals);
  std::cout << "(i)  Read " << _no_bins << " bins\n";
}

unsigned point_holder::count_no_bins(const image_long& binmap)
{
  int no = -1;
  for(unsigned y=0; y<binmap.yw(); ++y)
    for(unsigned x=0; x<binmap.xw(); ++x)
      {
	if( binmap(x, y) > int(no) )
	  no = binmap(x, y);
      }
  return no+1;
}

void point_holder::calc_points(const image_long& binmap,
			       const image_dbl& vals)
{
  std::vector<unsigned> counts(_no_bins);
  vec_dbl totx(_no_bins), toty(_no_bins);
  vec_dbl bval(_no_bins), total(_no_bins);

  // calculate centroid
  for(unsigned y=0; y<binmap.yw(); ++y)
    for(unsigned x=0; x<binmap.xw(); ++x)
      {
	long bin = binmap(x, y);
	if( bin >= 0 )
	  {
	    totx[bin] += x;
	    toty[bin] += y;
	    bval[bin] = vals(x, y);
	    total[bin] += vals(x, y);
	    counts[bin] ++;
	  }
      }

  for(unsigned bin=0; bin<_no_bins; ++bin)
    {
      assert( counts[bin] != 0 );
      _bininfos[bin].x = totx[bin]/counts[bin];
      _bininfos[bin].y = toty[bin]/counts[bin];
      _bininfos[bin].val = bval[bin];
      _bininfos[bin].total = total[bin];
      _bininfos[bin].pix_count = counts[bin];
    }

  // using centroid, we can get the widths
  vec_dbl totx_2(_no_bins), toty_2(_no_bins);

  for(unsigned y=0; y<binmap.yw(); ++y)
    for(unsigned x=0; x<binmap.xw(); ++x)
      {
	long bin = binmap(x, y);
	if( bin >= 0 )
	  {
	    // subtract mean x (preserves accuracy)
	    const double dx = x - _bininfos[bin].x;
	    const double dy = y - _bininfos[bin].y;

	    totx_2[bin] += dx*dx;
	    toty_2[bin] += dy*dy;
	  }
      }

  for(unsigned bin=0; bin<_no_bins; ++bin)
    {
      _bininfos[bin].rms_dx = std::sqrt( totx_2[bin]/counts[bin] );
      _bininfos[bin].rms_dy = std::sqrt( toty_2[bin]/counts[bin] );
    }
}

void point_holder::get_bin_info(unsigned bin,
				double* ret_x, double* ret_y,
				double* ret_val,
				double* ret_rms_dx, double* ret_rms_dy,
				double* ret_total,
				unsigned* ret_pix_count)
{
  assert( bin < _no_bins );

  const bininfo& bi = _bininfos[bin];

  *ret_x = bi.x;
  *ret_y = bi.y;
  *ret_val = bi.val;
  *ret_rms_dx = bi.rms_dx;
  *ret_rms_dy = bi.rms_dy;
  *ret_total = bi.total;
  *ret_pix_count = bi.pix_count;
}


int main(int argc, char* argv[])
{
  if( argc != 6 )
    {
      std::cerr << 
	"Usage:\n"
	" dumpdata infile.fits nerr_infile.fits perr_infile.fits "
	"binmap.fits out.dat\n";
      return 1;
    }

  const std::string indata = argv[1];
  const std::string indata_nerr = argv[2];
  const std::string indata_perr = argv[3];
  const std::string inbinmap = argv[4];
  const std::string outdata = argv[5];

  image_dbl* in_image;
  image_dbl* in_image_nerr;
  image_dbl* in_image_perr;
  image_long* in_binmap;

  // load in the main dataset
  {
    std::cout << "(i) Loading image " << indata << '\n';
    FITSFile dataset(indata);
    dataset.readImage(&in_image);
  }
  // load in the nerr dataset
  {
    std::cout << "(i) Loading image " << indata_nerr << '\n';
    FITSFile dataset(indata_nerr);
    dataset.readImage(&in_image_nerr);
  }
  // load in the perr dataset
  {
    std::cout << "(i) Loading image " << indata_perr << '\n';
    FITSFile dataset(indata_perr);
    dataset.readImage(&in_image_perr);
  }

  // load binmap
  {
    std::cout << "(i) Loading image " << inbinmap << '\n';
    FITSFile dataset(inbinmap);
    dataset.readImage(&in_binmap);
  }

  point_holder ph(*in_binmap, *in_image);
  point_holder ph_nerr(*in_binmap, *in_image_nerr);
  point_holder ph_perr(*in_binmap, *in_image_perr);

  std::ofstream fileout(outdata.c_str());
  assert( fileout );

  const unsigned no_bins = ph.no_bins();
  for(unsigned bin=0; bin<no_bins; ++bin)
    {
      double x, y, val, nerr, perr;
      double rms_dx, rms_dy, total;
      unsigned pixcount;

      // get data for bin
      ph_nerr.get_bin_info(bin, &x, &y, &nerr, &rms_dx, &rms_dy, &total,
			   &pixcount);
      ph_perr.get_bin_info(bin, &x, &y, &perr, &rms_dx, &rms_dy, &total,
			   &pixcount);
      ph.get_bin_info(bin, &x, &y, &val, &rms_dx, &rms_dy, &total,
		      &pixcount);

      // subtract
      const double err_n = nerr - val;
      const double err_p = perr - val;

      // kludge to make symmetric errors
      const double error = std::sqrt(0.5*(err_n*err_n + err_p*err_p));

      fileout << x << '\t' << y << '\t'
	      << val << '\t' << error << '\t'
	      << rms_dx << '\t' << rms_dy << '\t' << total
	      << '\t' << pixcount << '\t' << bin << '\t'
	      << err_n << '\t' << err_p << '\n';

    }
}
