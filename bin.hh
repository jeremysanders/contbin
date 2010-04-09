#ifndef BIN_HH
#define BIN_HH

#include <list>
#include <vector>

#include "misc.hh"
#include "point.hh"

// const size_t bin_no_neigh = 8;
// const int bin_neigh_x[bin_no_neigh] = {  0, -1, 1, 0, 1, 1, -1, -1 };
// const int bin_neigh_y[bin_no_neigh] = { -1,  0, 0, 1, 1,-1, 1, -1  };

const size_t bin_no_neigh = 4;
const int bin_neigh_x[bin_no_neigh] = {  0, -1, 1, 0 };
const int bin_neigh_y[bin_no_neigh] = { -1,  0, 0, 1 };

// simple squaring function
template<class T> inline T square(T v)
{
  return v*v;
}

// estimate the error squared on c counts
// uses formula of Gehrels 1986 ApJ, 303, 336) eqn 7
inline double error_sqd_est(double c)
{
  return square( 1. + sqrt(c + 0.75) );
}

// keep track of the parameters for the bin class
class bin_helper
{
public:
  bin_helper( const image_float* in_image,
	      const image_float* smoothed_image,
	      image_long* bins_image,
	      double threshold );

  void set_back( const image_float* back_image,
		 const image_float* expmap_image,
		 const image_float* bg_expmap_image )
  {
    _back_image = back_image;
    _expmap_image = expmap_image;
    _bg_expmap_image = bg_expmap_image;
  }

  void set_noisemap( const image_float* noisemap_image )
  {
    _noisemap_image = noisemap_image;
  }

  void set_mask( const image_short* mask_image )
  {
    _mask_image = *mask_image;
  }

  void set_constrain_fill( bool constrain_fill, double constrain_val )
  {
    _constrain_fill = constrain_fill;
    _constrain_val = constrain_val;
  }

  void set_scrub_large_bins( double fraction )
  {
    _scrub_large_bins = fraction;
  }

public:
  typedef std::vector<point_int> point_vec;
  typedef std::vector<point_vec> point_vec_vec;

public:
  // accessors
  const image_float* in_image() const { return _in_image; }
  const image_float* back_image() const { return _back_image; }

  const image_float* expmap_image() const { return _expmap_image; }
  const image_float* bg_expmap_image() const { return _bg_expmap_image; }
  const image_float* noisemap_image() const { return _noisemap_image; }

  const image_float* smoothed_image() const { return _smoothed_image; }
  const image_short* mask_image() const { return &_mask_image; }
  image_long* bins_image() const { return _bins_image; }

  double threshold() const { return _threshold; }
  unsigned xw() const { return _xw; }
  unsigned yw() const { return _yw; }
  unsigned max_annuli() const { return _max_annuli; }
  bool constrain_fill() const { return _constrain_fill; }
  double constrain_val() const { return _constrain_val; }
  double scrub_large_bins() const { return _scrub_large_bins; }

  // get the pixels at a particular radius
  const point_vec& points_at_annuli(const unsigned r) const
  { return _annuli_points[r]; }

  // return the next number for a bin
  long bin_counter() { const long t = _bin_counter; ++_bin_counter; return t; }

  // return how many bins have been processed
  long no_bins() const { return _bin_counter; }

  unsigned get_radius_for_area(unsigned area) const
  {
    return std::upper_bound(_areas.begin(), _areas.end(), area) -
      _areas.begin();
  }

private:
  // store whichs pixels line in an annuli with a particular radius
  void precalculate_annuli();
  // areas corresponding to each radius
  void precalculate_areas();

private:
  const image_float* _in_image;
  const image_float* _smoothed_image;
  image_long* const _bins_image;
  const double _threshold;

  const unsigned _xw;
  const unsigned _yw;

  const image_float* _back_image;
  const image_float* _expmap_image;
  const image_float* _bg_expmap_image;
  const image_float* _noisemap_image;
  image_short _mask_image;

  const unsigned _max_annuli;
  point_vec_vec _annuli_points;
  std::vector<unsigned> _areas;

  long _bin_counter;

  bool _constrain_fill;
  double _constrain_val;
  double _scrub_large_bins;
};

////////////////////////////////////////////////////////////////////////////
// class for constructing a bin
class bin
{
public:
  bin( bin_helper* helper );

  // delete all points in bin
  void drop_bin();

  // start bin with the specified pixel
  void do_binning(const unsigned x, const unsigned y);

  // return number of counts binned
  unsigned count() const { return _count; }

  // get signal in bin
  double signal() const
  {
    return _fg_sum - _bg_sum_weight;
  }

  // get noise in bin
  double noise_2() const
  {
    if( _helper->noisemap_image() == 0 )
      {
	// using background image
	double n = error_sqd_est(_fg_sum);

	if( _helper->back_image() != 0 )
	  n += (_expratio_sum_2 / _count) * error_sqd_est(_bg_sum);

	return n;
      }
    else
      {
	// using noisemap
	return _noisemap_2_sum;
      }
  }

  // get signal : noise squared
  double sn_2() const
  {
    const double csignal = signal();
    const double cnoise_2 = noise_2();

    if( cnoise_2 < 1e-7 )
      return 1e-7;
    else
      return csignal*csignal / cnoise_2;
  }

  // calculate ratio of edge length / a circle of same area
  //  bool check_constraint() const;
  bool check_constraint(const unsigned x, const unsigned y) const;

  typedef std::list<point_int> _point_list;
  typedef std::vector<point_int> _point_vector;
  typedef _point_vector _Pt_container;

  // get list of all pojnts in bin
  _Pt_container& get_all_points() { return _all_points; }

  // get list of all points on edge of bin
  _Pt_container& get_edge_points() { return _edge_points; }

  // get bin number
  long bin_no() const { return _bin_no; }
  // set number
  void set_bin_no(const long num) { _bin_no = num; }

  // add or remove point from the bin
  void add_point( const int x, const int y );
  void remove_point( const int x, const int y );

  // paint bin onto bins image
  void paint_bins_image() const;

private:
  bool add_next_pixel();

private:
  // helper things for binning
  bin_helper* _helper;

  // number in binmap of bin
  long _bin_no;

  // keep track of points on the edge of the bin
  _Pt_container _edge_points;
  // keep track of all the points
  _Pt_container _all_points;

  // pixel value we try to aim for
  double _aimval;

  // add up these
  double _fg_sum;              // foreground counts
  double _bg_sum;              // background counts
  double _bg_sum_weight;       // sum of the background*expratio
  double _noisemap_2_sum;      // sum of square of values from noisemap
  double _expratio_sum_2;      // sum of expratio^2

  // centroid
  point_dbl _centroid_sum;
  double _centroid_weight;
  unsigned _count;
};

typedef std::vector<bin> bin_vector;

#endif
