#ifndef IMAGE_DISK_ACCESS_HH
#define IMAGE_DISK_ACCESS_HH

#include <iostream>
#include <string>

#include "misc.hh"
#include "fitsio_simple.hh"

template<class T> void load_image(const std::string& filename, double *exposure,
				  T** image)
{
  std::cout << "(i) Loading " << filename << '\n';

  FITSFile ds(filename);
  ds.readImage(image);

  if( exposure != 0 )
    {
      double defval = 1;
      ds.readKey("EXPOSURE", exposure, &defval);
      std::cout << "(i)  exposure = " << *exposure << '\n';
    }
}

template<class T> void write_image(const std::string& filename, const T& img)
{
  std::cout << "(i) Writing " << filename << '\n';

  FITSFile ds(filename, FITSFile::Create);
  ds.writeImage(img);
}

#endif
