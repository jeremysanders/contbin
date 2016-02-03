#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <iostream>

#include "binner.hh"
#include "scrubber.hh"
#include "terminal.hh"

using namespace std;

// blank images
binner::binner( const image_float* in_image,
		const image_float* smoothed_image,
		double threshold )
  : _xw(in_image->xw()), _yw(in_image->yw()),
    _bins_image( _xw, _yw, -1 ),
    _binned_image( _xw, _yw ),
    _sn_image( _xw, _yw ),

    _bin_helper( in_image, smoothed_image, &_bins_image, threshold ),
    _bin_counter( 0 )
{
}

// class sorts pixels in reverse order (or order)
// based on their value in an image
class _flux_sort_points
{
public:
  _flux_sort_points(const image_float& im, bool bindown = true)
    : _im(im), _bindown(bindown)
  {
  }

  bool operator() ( const point_ushort p1,
		    const point_ushort p2 ) const
  {
    const float a = _im(p1.x(), p1.y());
    const float b = _im(p2.x(), p2.y());

    // order according to binning direction
    if( _bindown )
      return a > b;
    else
      return a < b;
  }

private:
  const image_float& _im;
  const bool _bindown;
};

// sort pixels into reverse flux order
void binner::sort_pixels(const bool bin_down)
{
  std::cout << "(i) Sorting pixels, binning from ";
  if(bin_down)
    std::cout << "top";
  else
    std::cout << "bottom";
  std::cout << "...";
  std::cout.flush();

  // add all pixels in image
  const image_short& in_mask = *_bin_helper.mask_image();
  _sorted_pixels.clear();

  for(unsigned y=0; y != _yw; ++y)
    for(unsigned x=0; x != _xw; ++x)
      {
	if( in_mask(x, y) >= 1 )
	  _sorted_pixels.push_back( point_ushort(x, y) );
      }

  // sort in reverse flux order
  std::sort( _sorted_pixels.begin(), _sorted_pixels.end(),
	     _flux_sort_points(*_bin_helper.smoothed_image(), bin_down) );

  // this is where we are in the list of sorted pixels
  _sorted_pix_posn = _sorted_pixels.begin();

  std::cout << " Done.\n";
}

// get number of unmasked pixels
unsigned binner::no_unmasked_pixels() const
{
  const image_short& in_mask = *_bin_helper.mask_image();

  unsigned no_unmasked = 0;
  for(unsigned y=0; y<_yw; ++y)
    for(unsigned x=0; x<_xw; ++x)
      {
	if( in_mask(x, y) >= 1. )
	  no_unmasked++;
      }
  return no_unmasked;
}

// find the highest flux pixel
point_int binner::find_next_pixel()
{
  const image_long& in_bins = *_bin_helper.bins_image();

  // iterate through sorted list until there are no pixels
  while( _sorted_pix_posn != _sorted_pixels.end() )
    {
      point_ushort p( *_sorted_pix_posn );

      // is pixel unbinned?
      if( in_bins(p.x(), p.y()) < 0 )
	{
	  return point_int( p.x(), p.y() );
	}

      ++ _sorted_pix_posn;
    }

  return point_int(-1, -1);
}

void binner::do_binning(const bool bin_down)
{
  // so we can interrupt binning
  terminal T;

  // sort pixels into flux order to find starting pixels
  sort_pixels(bin_down);

  const image_float* in_image = _bin_helper.in_image();
  const image_float* in_back = _bin_helper.back_image();

  // safety check for dimensions of images
  assert( _bins_image.xw()==_xw && _bins_image.yw()==_yw );
  assert( in_back == 0 || (in_back->xw()==_xw && in_back->yw() == _yw) );
  assert( in_image->xw()==_xw && in_image->yw()==_yw );
  assert( _sn_image.xw()==_xw && _sn_image.yw()==_yw );
  assert( _binned_image.xw()==_xw && _binned_image.yw()==_yw );

  std::cout << "(i) Starting binning\n";

  if( T.is_terminal() )
    {
      std::cout << "(i)  Press Esc to abort binning\n";
    }

  unsigned pix_counter = 0; // how many pixels processed
  const unsigned no_unmasked = no_unmasked_pixels();

  // get next pixel
  point_int nextpoint = find_next_pixel();
  assert( nextpoint.x() >= 0 && nextpoint.y() >= 0 );

  // repeat binnings, adding centroids and weights of bins
  // to above variables
  while( nextpoint.x() >= 0 && nextpoint.y() >= 0 )
    {
      // ESC pressed
      if( T.get_char() == 27 )
	{
	  cerr << "\nEsc pressed: aborting binning\n";
	  break;
	}

      // progress counter
      const long counter = _bin_helper.no_bins();
      if( counter % 10 == 0 && counter > 0)
	{
	  std::cout << std::setw(5) << counter << ' ';
	  std::cout.flush();
	  if( counter % 100 == 0 )
	    {
	      std::cout.setf(std::ios::fixed);
	      std::cout << " [" << std::setprecision(1)
			<< pix_counter*100./no_unmasked
			<< "%]\n";
	    }
	}

      // make the new bin and do the binning
      bin newbin( &_bin_helper );
      newbin.do_binning( nextpoint.x(), nextpoint.y() );
      _bins.push_back( newbin );

      // keep track of all the pixels binned
      pix_counter += newbin.count();

      // find the next pixel
      nextpoint = find_next_pixel();
    }

  _bin_counter = _bin_helper.no_bins();

  std::cout << " [100.0%]\n";

  std::cout << "(i) Done binning (" << _bin_counter << " bins)\n";
}

void binner::do_scrub()
{
  scrubber scrub( _bin_helper, _bins );
  scrub.scrub();

  if( _bin_helper.scrub_large_bins() > 0. )
    scrub.scrub_large_bins( _bin_helper.scrub_large_bins() );

  scrub.renumber();
}

// create output images and make histograms of signal/noise
void binner::calc_outputs()
{
  const size_t no_bins = _bins.size();
  std::vector<double> signal(no_bins);
  std::vector<double> noise_2(no_bins);
  std::vector<unsigned> pixcounts(no_bins);
  std::vector<double> sn(no_bins);

  double min_sn = 1e100, max_sn = -1e100;
  double min_signal = 1e100, max_signal = -1e100;

  // iterate over bins & collect info
  for(size_t i = 0; i != no_bins; ++i)
    {
      const bin& b = _bins[i];
      const long no = b.bin_no();
      if( no < 0 )
	continue;

      assert( no < long(no_bins) );

      signal[no] = b.signal();
      max_signal = std::max( signal[no], max_signal );
      min_signal = std::min( signal[no], min_signal );

      noise_2[no] = b.noise_2();
      pixcounts[no] = b.count();

      sn[no] = std::sqrt( b.sn_2() );

      if( ! std::isfinite(sn[no]) || sn[no] < 0 )
        {
          std::cerr << "WARNING: invalid value in signal to noise. "
            "This can be caused by a negative input image.\n";
        }

      max_sn = std::max( sn[no], max_sn );
      min_sn = std::min( sn[no], min_sn );
    }

  // now make output images
  _sn_image.set_all(-1);
  _binned_image.set_all(-1);
  for(unsigned y = 0; y != _yw; ++y)
    for(unsigned x = 0; x != _xw; ++x)
      {
	const int bin = _bins_image(x, y);
	if( bin >= 0 )
	  {
	    _sn_image(x, y) = sn[bin];
	    _binned_image(x, y) = signal[bin] / pixcounts[bin];
	  }
      }

  // build histogram of signal to noises
  {
    const unsigned no_hbins = 30;
    const double delta_sn = (max_sn-min_sn+0.0001)/no_hbins;
    const double delta_signal = (max_signal-min_signal+0.0001)/no_hbins;
    std::vector<unsigned> histo_sn(no_hbins);
    std::vector<unsigned> histo_signal(no_hbins);

    for(size_t bin = 0; bin != no_bins; ++bin)
      {
        if( _bins[bin].bin_no() < 0 )
          continue;

	const unsigned index_sn = unsigned( (sn[bin]-min_sn) / delta_sn );
	const unsigned index_signal = unsigned( (signal[bin]-min_signal) /
						delta_signal );

	assert( index_sn < no_hbins && index_signal < no_hbins );

	++histo_sn[index_sn];
	++histo_signal[index_signal];
      }

    // output histogram data in file
    std::ofstream stream_sn("bin_sn_stats.qdp");
    std::ofstream stream_signal("bin_signal_stats.qdp");
    stream_sn <<
      "label x Signal:Noise\n"
      "label y Number\n"
      "line step\n";
    stream_signal <<
      "label x Counts\n"
      "label y Number\n"
      "line step\n";

    // write out histograms
    for(unsigned h = 0; h != no_hbins; ++h)
      {
	stream_sn << min_sn + (h+0.5)*delta_sn << '\t'
		  << histo_sn[h] << '\n';
	stream_signal << min_signal + (h+0.5)*delta_signal << '\t'
		      << histo_signal[h] << '\n';
      }
  }
}
