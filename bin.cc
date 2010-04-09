#include <cmath>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <math.h>

#include "point.hh"
#include "misc.hh"
#include "bin.hh"

namespace
{
  // work out integerised radius
  inline unsigned unsigned_radius(const int x, const int y)
  {
    return unsigned( std::sqrt( double(x*x + y*y) ) );
  }
}

////////////////////////////////////////////////////////////////////////////

bin_helper::bin_helper( const image_float* in_image,
			const image_float* smoothed_image,
			image_long* bins_image,
			double threshold )
  
  : _in_image(in_image),
    _smoothed_image(smoothed_image),
    _bins_image(bins_image),
    _threshold(threshold),
    
    _xw( in_image->xw() ), _yw( in_image->yw() ),

    _back_image( 0 ),
    _expmap_image( 0 ),
    _bg_expmap_image( 0 ),
    _noisemap_image( 0 ),
    _mask_image( _xw, _yw, 1 ),

    _max_annuli( unsigned_radius(_xw, _yw) + 1 ),
    _bin_counter( 0 ),

    _constrain_fill( false ),
    _constrain_val( 4 ),
    
    _scrub_large_bins( -1 )

{
  precalculate_annuli();
  precalculate_areas();
}

void bin_helper::precalculate_annuli()
{
  _annuli_points.clear();
  _annuli_points.resize( _max_annuli );

  // identify which pixels lie at a particular radius
  for(int y=-int(_yw-1); y<int(_yw); ++y)
    {
      for(int x=-int(_xw-1); x<int(_xw); ++x)
        {
          const unsigned r = unsigned_radius(x, y);

          // add pixel to appropriate radius
          _annuli_points[r].push_back( point_int(x, y) );
        }
    }
}

void bin_helper::precalculate_areas()
{
  _areas.resize( _max_annuli );

  unsigned total = 0;
  for(unsigned radius=0; radius<_max_annuli; ++radius)
    {
      const unsigned area = _annuli_points[radius].size();
      total += area;
      _areas[radius] = total;
    }
}

///////////////////////////////////////////////////////////////////////

// find an initial pixel for the bin
// search for the nearest pixel with the highest flux

bin::bin( bin_helper* helper )
  : _helper(helper),
    _bin_no( _helper->bin_counter() ),
    _aimval( -1 ),
    _fg_sum(0), _bg_sum(0),
    _bg_sum_weight(0), _noisemap_2_sum(0),_expratio_sum_2(0),
    _centroid_sum(0, 0), _centroid_weight(0), _count(0)
{
}

void bin::drop_bin( )
{
  // drop all points from bin
  _fg_sum = 0;
  _bg_sum = 0;
  _bg_sum_weight = 0;
  _noisemap_2_sum = 0;
  _expratio_sum_2 = 0;
  _centroid_sum.x() = 0;
  _centroid_sum.y() = 0;
  _centroid_weight = 0;
  _count = 0;
  _all_points.clear();
}

void bin::remove_point( const int x, const int y )
{
  // get rid of point from lists
  {
    const point_int P(x,y);

    _Pt_container::iterator pt = std::find( _all_points.begin(),
					    _all_points.end(),
					    P );
    assert( pt != _all_points.end() );
    _all_points.erase( pt );

    {
      _Pt_container::iterator ept = std::find( _edge_points.begin(),
					       _edge_points.end(),
					       P );
      if( ept != _edge_points.end() )
	  _edge_points.erase( ept );
    }
  }

  image_long* const bins_image = _helper->bins_image();

  // now get rid of the counts
  _fg_sum -= (*_helper->in_image())(x, y);
  _count--;
  (*bins_image) (x, y) = -1;

  if( _helper->back_image() != 0 )
    {
      const double bs = (*_helper->expmap_image())(x, y) /
	(*_helper->bg_expmap_image())(x, y);
      const double bg = (*_helper->back_image())(x, y);

      _bg_sum -= bg;
      _bg_sum_weight -= bg*bs;
      _expratio_sum_2 -= bs*bs; 
    }

  // remove from noisemap sum (if any supplied)
  if( _helper->noisemap_image() != 0 )
    {
      _noisemap_2_sum -= square( (*_helper->noisemap_image())(x, y) );
    }

  // add points that are on the edge of this point in this bin into the
  // edge list

  const int xw = _helper->xw();
  const int yw = _helper->yw();

  for(size_t n = 0; n != bin_no_neigh; ++n)
    {
      const int xp = x + bin_neigh_x[n];
      const int yp = y + bin_neigh_y[n];

      // if neighbour is in this bin, mark as edge
      // if not already marked so
      if( xp >= 0 && yp >= 0 && xp < xw && yp < yw &&
	  (*bins_image)(xp, yp) == _bin_no )
	{
	  if( std::find( _edge_points.begin(), _edge_points.end(),
			 point_int(xp, yp) ) == _edge_points.end() )
	    {
	      _edge_points.push_back( point_int(xp, yp) );
	    }
	}
    } // loop over neighbours
}

void bin::add_point(const int x, const int y)
{
  _all_points.push_back( point_int(x, y) );

  double signal = (*_helper->in_image())(x, y);
  _fg_sum += signal;
  _count++;
  (*_helper->bins_image()) (x, y) = _bin_no;

  if(  _helper->back_image() != 0 )
    {
      const double bs = (*_helper->expmap_image())(x, y) /
	(*_helper->bg_expmap_image())(x, y);
      const double back = (*_helper->back_image())(x, y);
      _bg_sum += back;
      _bg_sum_weight += back*bs;
      _expratio_sum_2 += bs*bs;

      signal -= back*bs;
    }

  // add to noisemap sum (if any supplied)
  if( _helper->noisemap_image() != 0 )
    {
      _noisemap_2_sum += square( (*_helper->noisemap_image())(x, y) );
    }

  // update centroid
  {
    const double cs = signal < 1e-7 ? 1e-7 : signal;
    _centroid_sum += point_dbl( x, y ) * cs;
    _centroid_weight += cs;
  }

  // put into edge (it might not be, but it will get flushed out)
  if( std::find( _edge_points.begin(), _edge_points.end(),
		 point_int(x, y) ) == _edge_points.end() )
    {
      _edge_points.push_back( point_int(x, y) );
    }
}

// paint bin onto bins_image
void bin::paint_bins_image() const
{
  image_long& bins_image = * _helper->bins_image();

  typedef _Pt_container::const_iterator CI; 
  const CI e = _all_points.end();
  for( CI pix = _all_points.begin(); pix != e; ++pix )
    {
      bins_image( pix->x(), pix->y() ) = _bin_no;
    }
}

// adds the next pixel to the bin!
bool bin::add_next_pixel()
{
  // easier access to images
  const int xw = _helper->xw();
  const int yw = _helper->yw();
  const image_short& mask_image = *_helper->mask_image();
  const image_long& bins_image = *_helper->bins_image();
  const image_float& smoothed_image = *_helper->smoothed_image();
  const bool constrain_fill = _helper->constrain_fill();

  // find pixel nearest in value to aimval, looking around the edge pixels
  // of the bin
  double delta = 1e99;
  int bestx = -1;
  int besty = -1;

  // iterate over the proper edge points
  //const bool toolarge = calc_length_ratio() > constrain_val;
  _Pt_container::iterator pix = _edge_points.begin();
  while( pix != _edge_points.end() )
    {
      // where the pixel is
      const int x = pix->x();
      const int y = pix->y();

      // check whether point is edge
      bool is_edge = false;

      // iterate over neighbours, looking for any better than we have
      for( size_t n = 0; n != bin_no_neigh; ++n )
	{
	  const int xp = x + bin_neigh_x[n];
	  const int yp = y + bin_neigh_y[n];

	  if( xp >= 0 && yp >= 0 && xp < xw && yp < yw )
	    {
	      const long bin = bins_image(xp, yp);
	      if( bin != _bin_no )
		is_edge = true;

	      // this pixel isn't taken
	      if( bin < 0 && mask_image(xp, yp) == 1 )
		{
		  if ( ! constrain_fill || check_constraint(xp, yp) )
		    {
		      const double newdelta = fabs( smoothed_image(xp, yp) -
						    _aimval );
		      if( newdelta < delta )
			{
			  delta = newdelta;
			  bestx = xp; besty = yp;
			}
		    }
		} // pixel not taken
	    } // if in image

	} // loop over neighbours

      // go to next point (removing current point if not edge)
      if( is_edge )
	pix++;
      else
	pix = _edge_points.erase(pix);

    } // loop over edge pixels in bin

  // we didn't find any nice neighbours
  if( bestx == -1 )
    {
      return false;
    }

  // update stuff
  add_point( bestx, besty );

  return true;
}

// do the binning until the threshold reached
void bin::do_binning(const unsigned x, const unsigned y)
{
  _aimval = (*_helper->smoothed_image())(x, y);
  add_point(x, y);

  const double sn_threshold_2 = _helper->threshold()*_helper->threshold();

  // keep adding pixels until add_next_pixel complains, or s/n reached
  while( sn_2() < sn_threshold_2 )
    {
      if( ! add_next_pixel() )
	break;
    }
}

// is the constraint still satisfied if we add this pixel?
bool bin::check_constraint(const unsigned x, const unsigned y) const
{


//   if( _count < 20 )
//     {
//       return true;
//     }

  const point_dbl c = _centroid_sum / _centroid_weight;
  const double dx = c.x() - x;
  const double dy = c.y() - y;
  const double r2 = dx*dx + dy*dy;

  // get radius for area
  const unsigned circradius = _helper->get_radius_for_area(_count)+1;

//   std::cout << _count << ' ' << circradius << ' ' <<
//     _count / (circradius*circradius*M_PI) << '\n';

  // estimate of radius of circle sqd, A = PI r^2
  // const double circ_radius_2 =  _count / M_PI;

  return (r2 / (circradius*circradius)) < square(_helper->constrain_val());
}
