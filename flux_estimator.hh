#ifndef FLUX_ESTIMATOR_HH
#define FLUX_ESTIMATOR_HH

#include <vector>
#include "misc.hh"

class flux_estimator
{
public:
  flux_estimator(const image_float* const in_image,
		 const image_float* const back_image,
		 const image_short* const mask_image,

		 const image_float* const expmap_image,
		 const image_float* const bg_expmap_image,

		 const image_float* const noisemap_image,

		 const double minsn );

  const image_float& operator()();

  struct _point
  {
    _point(int xp, int yp) : x(xp), y(yp) {}
    int x, y;
  };
  typedef std::vector<_point> _point_vec;
  typedef std::vector<_point_vec> _point_vec_vec;

private:
  void do_estimation();

  // work out which points are in which annuli
  void precalculate_annuli();
  void smooth();

  void bin();

private:

  const unsigned _xw, _yw; // dimensions of images
  const double _minsn; // minimum signal:noise

  const image_float* const _in_image;
  const image_float* const _back_image;
  const image_short* const _mask_image;
  const image_float* const _expmap_image;
  const image_float* const _bg_expmap_image;
  const image_float* const _noisemap_image;

  // precalculated list of which points are in which annuli
  const unsigned _max_annuli;
  _point_vec_vec _annuli_points;

  bool _done;

  image_float _iteration_image; // output image
  image_float _estimated_errors; // errors on iteration
};

#endif
