#include <iostream>
#include <algorithm>
#include <cassert>
#include <iomanip>

#include "scrubber.hh"

using namespace std;

// square number
inline double square(const double d)
{
  return d*d;
}

scrubber::scrubber( bin_helper& helper, bin_vector& bins )
  : _helper(helper),
    _bins( bins ),
    _no_bins( _bins.size() ),

    _scrub_sn_2( square(helper.threshold() ) ),

    _cannot_dissolve( _no_bins ),

    _xw( helper.xw() ), _yw( helper.yw() )
{
}

void scrubber::find_best_neighbour(bin* thebin, bool allow_unconstrained,
				   int* bestx, int* besty,
				   int* bestbin )
{
  const image_float& smoothed_image = *_helper.smoothed_image();
  const image_long& bins_image = * _helper.bins_image();
  const long binno = thebin->bin_no();

  double bestdelta = 1e99;

  bin::_Pt_container& edgepoints = thebin->get_edge_points();

  size_t pt = 0;
  while( pt < edgepoints.size() )
    {
      const int x = edgepoints[pt].x();
      const int y = edgepoints[pt].y();
      const double v = smoothed_image(x, y);

      // loop over neighbours of edge point
      bool anyneighbours = false;
      for( size_t n = 0; n != bin_no_neigh; ++n )
	{
	  const int xp = x + bin_neigh_x[n];
	  const int yp = y + bin_neigh_y[n];

	  // select pixels in neighbouring bins
	  if( xp >= 0 && yp >= 0 && xp < int(_xw) && yp < int(_yw) )
	    {
	      // if the neighbour is a real and different bin
	      const long nbin = bins_image(xp, yp);
	      if( nbin != -1 && nbin != binno )
		{
		  anyneighbours = true;

		  // we skip neighbours with too long a constraint if reqstd
		  if( _helper.constrain_fill() &&
		      ! allow_unconstrained &&
		      ! _bins[nbin].check_constraint(xp, yp) )
		    continue;

		  const double delta = fabs( v - smoothed_image(xp, yp) );
		  if( delta < bestdelta )
		    {
		      bestdelta = delta;
		      *bestx = x;
		      *besty = y;
		      *bestbin = bins_image(xp, yp);
		    }
		} // valid neighbour
	    } // in image range
	} // loop over neighbours

      // remove edge pixels without any neighbours
      if( ! anyneighbours )
        edgepoints.erase(edgepoints.begin() + pt);
      else
	++pt;
    }

}

void scrubber::dissolve_bin( bin* thebin )
{
  // loop until no pixels remaining
  while( thebin->count() != 0 )
    {
      int bestx = -1;
      int besty = -1;
      int bestbin = -1;

      // find best neighbour within constraints
      find_best_neighbour( thebin, false, &bestx, &besty, &bestbin );

      // if none, then ignore constraints
      if( bestx == -1 && _helper.constrain_fill() )
	find_best_neighbour( thebin, true, &bestx, &besty, &bestbin );

      // stop dissolving bin if we have no neigbours for our remaining
      // pixels
      if( bestbin == -1 )
	{
	  const long binno = thebin->bin_no();
	  cout << "WARNING: Could not dissolve bin "
	       << binno << " into surroundings\n";
	  _cannot_dissolve[ binno ] = true;
	  return;
	}

      // reassign pixel
      thebin->remove_point(bestx, besty);
      _bins[ bestbin ].add_point(bestx, besty);
    }
}

void scrubber::scrub()
{
  std::cout << "(i) Starting scrubbing...\n";

  // put bins into a pointer array so we can discard them quickly when
  // we don't need to consider them
  typedef std::vector<bin*> BPV;
  BPV bin_ptrs;
  bin_ptrs.reserve( _no_bins );
  for( unsigned i = 0; i != _no_bins; ++i )
    {
      if( _bins[i].sn_2() < _scrub_sn_2 )
	bin_ptrs.push_back( &_bins[i] );
    }

  // we keep looping until the lowest S/N bin is removed
  for( ;; )
    {
      double lowest_SN_2 = 1e99;

      // iterate over bins, find those with the lowest S/N

      BPV::iterator i = bin_ptrs.begin();
      BPV::iterator lowest_bin = i-1;
      while( i != bin_ptrs.end() )
	{
	  const double SN_2 = (*i)->sn_2();
	  // if this bin has a larger S/N than threshold, remove it
	  // iterators below this one should still be valid
	  if( SN_2 >= _scrub_sn_2 )
	    {
	      i = bin_ptrs.erase(i);
	    }
	  else
	    {
	      // if this is lower than before, store it
	      if( SN_2 < lowest_SN_2 )
		{
		  lowest_SN_2 = SN_2;
		  lowest_bin = i;
		}
	      ++i;
	    }
	} // loop over bins

      // exit if no more bins remaining
      if( lowest_bin < bin_ptrs.begin() || lowest_SN_2 >= _scrub_sn_2 )
	break;

      // get rid of that bin (if it cannot be disolved, it doesn't matter
      dissolve_bin( *lowest_bin );
      bin_ptrs.erase( lowest_bin );

      // show progress to user
      if( bin_ptrs.size() % 10 == 0 )
        {
          std::cout << std::setw(5) << bin_ptrs.size() << ' ';
          std::cout.flush();
          if( bin_ptrs.size() % 100 == 0 )
            std::cout << '\n';
        }
    }

  std::cout << "(i) Done\n";
}

void scrubber::scrub_large_bins(double fraction)
{
  std::cout << "(i) Scrubbing bins with fraction of area > "
	    << fraction << "...\n";

  typedef bin_vector::iterator BI;
  const BI e = _bins.end();

  // get total number of pixels in bins
  unsigned totct = 0;
  for( BI i = _bins.begin(); i != e; ++i )
    totct += i->count();

  // no get rid of large bins
  for( BI i = _bins.begin(); i != e; ++i )
    {
      const double thisfrac = double(i->count()) / totct;
      if( thisfrac >= fraction )
	{
	  std::cout << " Scrubbing bin " << i->bin_no()
		    << '\n';

	  i->drop_bin();
	}
    }
}

void scrubber::renumber()
{
  std::cout << "(i) Starting renumbering...\n";

  {
    // split bins into those with counts and those without
    bin_vector::iterator e = std::partition( _bins.begin(), _bins.end(),
					     std::mem_fun_ref(&bin::count) );
    _bins.erase(e, _bins.end());
  }

  // now clear bin image, and repaint everything (doing renumber)
  _helper.bins_image()->set_all( -1 );

  long number = 0;
  typedef bin_vector::iterator BI;
  const BI e = _bins.end();
  for( BI i = _bins.begin(); i != e; ++i )
    {
      i -> set_bin_no( number );
      i -> paint_bins_image();
      number++;
    }

  std::cout << "(i)  " << number << " bins when finished\n"
	    << "(i) Done\n";
}
