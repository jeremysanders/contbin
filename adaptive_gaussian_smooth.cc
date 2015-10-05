#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>
#include <string>

#include "parammm/parammm.hh"
#include "memimage.hh"
#include "fitsio_simple.hh"
#include "image_disk_access.hh"

template<class T> T sqd(T v)
{
  return v*v;
}

class Kernels
{
public:
  const image_float* getKernel(unsigned idx, float sigma);
  ~Kernels();

private:
  std::vector<image_float*> kernels;
};

const image_float* Kernels::getKernel(unsigned idx, float sigma)
{
  while(idx >= kernels.size())
    kernels.push_back(0);

  if(kernels[idx] != 0)
    return kernels[idx];

  // construct kernel
  unsigned nsigma = 3;
  unsigned width = unsigned(std::ceil(sigma*nsigma))*2+1;
  image_float* kern = new image_float(width, width);

  float invsigma2 = -0.5f/sqd(sigma);

  for(unsigned y=0; y<width; ++y)
    {
      int dy2 = sqd(int(y)-int(width/2));
      for(unsigned x=0; x<width; ++x)
        {
          int dx2 = sqd(int(x)-int(width/2));
          (*kern)(x, y) = std::exp((dx2+dy2)*invsigma2);
        }
    }

  kernels[idx] = kern;
  return kern;
}

Kernels::~Kernels()
{
  for(std::vector<image_float*>::iterator i=kernels.begin();
      i!=kernels.end(); ++i)
    delete (*i);
}

struct KernResult
{
  float avexpcorr;
  float avexpmap;
};

KernResult getKernApplied(unsigned x, unsigned y,
                          const image_float& kern,
                          const image_float& expcorrimg,
                          const image_float& expmapimg,
                          const image_float& maskimg)
{
  unsigned kernsize = kern.xw();

  // clip edges of convolution to edge of input image
  // this is faster than manually checking in the loop
  unsigned xw = expcorrimg.xw();
  unsigned yw = expcorrimg.yw();
  unsigned kx0 = x>kernsize/2 ? 0 : kernsize/2-x;
  unsigned kx1 = x+kernsize/2<xw ? kernsize : kernsize/2+xw-x;
  unsigned ky0 = y>kernsize/2 ? 0 : kernsize/2-y;
  unsigned ky1 = y+kernsize/2<yw ? kernsize : kernsize/2+yw-y;

  // 2d convolution
  float sum = 0;
  float sumweight = 0;
  float sumexpmap = 0;

  for(unsigned ky=ky0; ky<ky1; ++ky)
    {
      unsigned cy = y-kernsize/2+ky;
      for(unsigned kx=kx0; kx<kx1; ++kx)
        {
          unsigned cx = x-kernsize/2+kx;
          float k = kern(kx, ky) * maskimg(cx, cy);
          sum += expcorrimg(cx, cy)*k;
          sumexpmap += expmapimg(cx, cy)*k;
          sumweight += k;
        }
    }

  KernResult result = {sum*(1.f/sumweight), sumexpmap*(1.f/sumweight)};
  return result;
}

void applySmoothing(const image_float& expcorrimg,
                    const image_float& expmapimg,
                    const image_float& maskimg,
                    float snthresh,
                    image_float* outimg)
{
  Kernels kernels;

  const unsigned xw = expcorrimg.xw();
  const unsigned yw = expcorrimg.yw();
  const float sn2thresh = sqd(snthresh);

  for(unsigned y=0; y<yw; ++y)
    {
      if(y % 10 == 0)
        std::cout << y << '\n';
      for(unsigned x=0; x<xw; ++x)
        {
          if(maskimg(x, y) <= 0)
            continue;

          for(unsigned sidx=1; sidx<2000; ++sidx)
            {
              float sigma = sidx*0.25f;
              const image_float* kern = kernels.getKernel(sidx, sigma);
              KernResult res = getKernApplied(x, y, *kern, expcorrimg,
                                              expmapimg, maskimg);
              float cts = (res.avexpcorr*res.avexpmap)*float(M_PI)*sqd(2*sigma);
              float sn2 = cts;

              if(sn2 >= sn2thresh)
                {
                  (*outimg)(x, y) = res.avexpcorr;
                  // exit increasing sigma
                  break;
                }

            } // loop sigma

        } // loop x
    } // loop y
}

// make mask in float (so it can be easily multiplied)
// (also converts masked pixels to 0 in input image)
void makeFloatMask(const image_short& maskimg, image_float& expcorrimg,
                   image_float& maskflt)
{
  const unsigned xw = expcorrimg.xw();
  const unsigned yw = expcorrimg.yw();

  for(unsigned y=0; y<yw; ++y)
    for(unsigned x=0; x<xw; ++x)
      {
        int mask = maskimg(x,y) != 0 && !std::isnan(expcorrimg(x,y));

        maskflt(x, y) = mask;
        if(mask == 0)
          expcorrimg(x, y) = 0;
      }
}

int main(int argc, char *argv[])
{
  double sn=15;
  std::string maskfile;
  std::string outfile = "ags.fits";

  parammm::param params(argc, argv);
  params.add_switch( parammm::pswitch( "mask", 'm',
				       parammm::pstring_opt(&maskfile),
				       "set mask file",
				       "FILE"));
  params.add_switch( parammm::pswitch( "out", 'o',
				       parammm::pstring_opt(&outfile),
				       "set output file (def ags.fits)",
				       "FILE"));

  params.add_switch( parammm::pswitch("sn", 's',
				      parammm::pdouble_opt(&sn),
				      "set signal:noise threshold (def 15)",
				      "VAL"));

  params.set_autohelp("Usage: adaptive_gaussian_smooth "
                      "[OPTIONS] expcorr.fits expmap.fits\n"
		      "Adaptive Gaussian smoothing program.\n"
		      "Written by Jeremy Sanders 2015.",
		      "Report bugs to <jeremy@jeremysanders.net>");
 
  params.enable_autohelp();
  params.enable_autoversion("0.1",
			    "Jeremy Sanders",
			    "Licenced under the GPL - see the file COPYING");
  params.enable_at_expansion();
  params.interpret_and_catch();


  if(params.args().size() != 2)
    {
      params.show_autohelp();
    }

  std::string expcorrfile = params.args()[0];
  std::string expmapfile = params.args()[1];

  image_float* expcorrimg;
  load_image(expcorrfile, 0, &expcorrimg);
  image_float* expmapimg;
  load_image(expmapfile, 0, &expmapimg);
  image_short* maskimg;
  if( !maskfile.empty() )
    {
      load_image(maskfile, 0, &maskimg);
    }
  else
    {
      maskimg = new image_short(expcorrimg->xw(), expcorrimg->yw());
      maskimg->set_all(1);
    }

  image_float maskflt(maskimg->xw(), maskimg->yw());
  makeFloatMask(*maskimg, *expcorrimg, maskflt);

  image_float outimg(expcorrimg->xw(), expcorrimg->yw());
  outimg.set_all(std::numeric_limits<float>::quiet_NaN());
  applySmoothing(*expcorrimg, *expmapimg, maskflt, sn, &outimg);

  write_image(outfile, outimg);

  delete expcorrimg;
  delete expmapimg;
  delete maskimg;

  return 0;
};
