#ifndef SCRUBBER_HH
#define SCRUBBER_HH

#include "bin.hh"

// go over image assigning stray bins to other bins
class scrubber
{
public:
  scrubber( bin_helper& helper,
	    bin_vector& bins );

  void scrub();

  // renumber bins, throwing away bins with zero counts
  void renumber();

  // get rid of bins with fraction of pixels >= fraction
  void scrub_large_bins( double fraction );

private:
  void dissolve_bin( bin* diss_bin );
  void find_best_neighbour(bin* thebin, bool allow_unconstrained,
			   int* bestx, int* besty, int* bestbin);

private:
  bin_helper& _helper;
  bin_vector& _bins;

  const unsigned _no_bins;
  const double _scrub_sn_2;   // minimum signal:noise squared allowed

  std::vector<bool> _cannot_dissolve; // set if cannot dissolve a bin

  const unsigned _xw, _yw;
};


#endif
