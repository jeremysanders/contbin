#ifndef BINNER_BINACC_HH
#define BINNER_BINACC_HH

#include "misc.hh"
#include "point.hh"
#include "bin.hh"

#include <list>


class binner
{
public:
  // set up binner with input image (poisson statistics)
  binner( const image_float* in_image,
	  const image_float* smoothed_image,
	  double threshold );

  // set background image (poisson statistics)
  // correction factor is multiplies error on back_image
  // (if exposure of back_image != input image)
  void set_back_image( const image_float* back_image,
		       const image_float* expmap_image,
		       const image_float* bg_expmap_image )
  {
    _bin_helper.set_back( back_image, expmap_image, bg_expmap_image );
  }

  void set_noisemap_image( const image_float* noisemap_image )
  {
    _bin_helper.set_noisemap( noisemap_image );
  }

  void set_mask_image( const image_short* mask_image )
  {
    _bin_helper.set_mask( mask_image );
  }

  void set_constrain_fill( bool constrain_fill, double constrain_val )
  {
    _bin_helper.set_constrain_fill( constrain_fill, constrain_val );
  }

  void set_scrub_large_bins( double fraction )
  {
    _bin_helper.set_scrub_large_bins( fraction );
  }

  // do the binning
  // bin_down true if start at highest pixels, false at lowest
  void do_binning(const bool bin_down);

  // scrub bins
  void do_scrub();

  // calculate output images (returned below)
  void calc_outputs();

  // get output image, binmap, and signal:noise image
  const image_float& get_output_image() const { return _binned_image; };
  const image_long& get_binmap_image() const { return _bins_image; };
  const image_float& get_sn_image() const { return _sn_image; };

private:
  // find the pixel with the highest smoothed flux
  point_int find_next_pixel();

  // no unmasked pixels
  unsigned no_unmasked_pixels() const;

  // sort smoothed pixels
  void sort_pixels(const bool bin_down);

private:
  const unsigned _xw, _yw;  // size of input images

  image_long _bins_image;  // currently set bins
  image_float _binned_image;  // output image
  image_float _sn_image;  // signal:noise image

  bin_helper _bin_helper;
  unsigned _bin_counter;

  bin_vector _bins; // keep all of the bins

  typedef std::vector< point_ushort >  _Pt_sorted_vec;

  _Pt_sorted_vec _sorted_pixels;
  _Pt_sorted_vec::const_iterator _sorted_pix_posn;
};

#endif
