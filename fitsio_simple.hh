// quick routines for doing fitsio reading and writing of files
// Jeremy Sanders (2005)
// Released under the GPL

#ifndef FITSIO_SIMPLE_JSS__HH
#define FITSIO_SIMPLE_JSS__HH

#include <string>
#include <sstream>
#include <iostream>
//#include <cfitsio/fitsio.h>
#include <fitsio.h>

#include "memimage.hh"

//////////////////////////////////////////////////////////////////////
// FITS file interface

class FITSFile
{
public:
  enum OpenMode {RO, RW, Create};
    // modes: existing file, read-only
    //        existing file, read-write
    //        create new file

public:
  FITSFile();
  FITSFile(const std::string& filename, OpenMode mode = RO);
  ~FITSFile();

  // open and close files
  void open(const std::string& filename, OpenMode mode = RO);
  void close();

  // header keyword functions
  template<class T> void readKey(const std::string& key,
				 T* val,
				 const T* defaultval = 0,
				 std::string* comment = 0);
  void readKey(const std::string& key, std::string* val,
	       const std::string* defaultval = 0,
	       std::string* comment = 0);
  template<class T> void updateKey(const std::string& key,
				   const T val,
				   const std::string* comment = 0);
  void updateKey(const std::string& key, const std::string& val,
		 const std::string* comment = 0);
  void updateKey(const std::string& key, const char* val,
		 const std::string* comment = 0)
  {
    updateKey(key, std::string(val), comment);
  }
  void writeDate();
  void writeComment(const std::string& comment);
  void writeHistory(const std::string& history);

  void copyHeaderTo(FITSFile& other);

  void writeDatestamp(const std::string& program);

  // image reading and writing
  template<class T> void readImage(dm::memimage<T>** image);
  template<class T> void writeImage(const dm::memimage<T>& image);

private:
  void _checkStatus(const std::string& operation);

private:
  fitsfile* _file;
  int _status;
  bool _verbose;

  std::string _filename;
};


////////////////////////////////////////////////////////////////////////
// implementation below (for template routines)

// overloaded functions to convert c datatypes to fitsio datatypes
#define _DEFINE_CNVT(CTYPE, FITSIMGTYPE, FITSTYPE) \
inline int _FITSImg_Datatype(const dm::memimage<CTYPE>* image) \
{ return FITSIMGTYPE; } \
inline int _FITSVal_Datatype(const CTYPE* val) \
{ return FITSTYPE; }
// saves lots of space
_DEFINE_CNVT(unsigned char, BYTE_IMG, TBYTE);
_DEFINE_CNVT(signed char, SBYTE_IMG, TSBYTE);
_DEFINE_CNVT(short, SHORT_IMG, TSHORT);
_DEFINE_CNVT(unsigned short, USHORT_IMG, TUSHORT);
_DEFINE_CNVT(int, LONG_IMG, TINT);
_DEFINE_CNVT(long, LONG_IMG, TLONG);
_DEFINE_CNVT(unsigned long, ULONG_IMG, TULONG);
_DEFINE_CNVT(float, FLOAT_IMG, TFLOAT);
_DEFINE_CNVT(double, DOUBLE_IMG, TDOUBLE);
#undef _DEFINE_CNVT

template<class T> void FITSFile::readKey(const std::string& key, T* val,
					 const T* defaultval,
					 std::string* comment)
{
  const int fits_datatype = _FITSVal_Datatype(val);
  // yes, a stupid buffer overflow possibility here, but that's the
  // cfitsio api
  char commentbuf[512];

  const int retval = fits_read_key(_file, fits_datatype,
				   const_cast<char*>(key.c_str()),
				   val, commentbuf, &_status);

  if( retval != 0 )
    {
      if(defaultval == 0)
	{
	  // else print error message
	  std::ostringstream o;
	  o << "Reading header keyword " << key;
	  _checkStatus(o.str());
	  // doesn't get here...
	}

      *val = *defaultval;
      fits_clear_errmsg();
      _status = 0;
    }

  if( comment != 0 )
    *comment = commentbuf;
}

template<class T> void FITSFile::updateKey(const std::string& key,
					   const T val,
					   const std::string* comment)
{
  const int fits_datatype = _FITSVal_Datatype(&val);

  if(comment != 0)
    {
      fits_update_key(_file, fits_datatype,
		      const_cast<char*>(key.c_str()),
		      const_cast<T*>(&val),
		      const_cast<char*>(comment->c_str()), &_status);
    }
  else
    {
      fits_update_key(_file, fits_datatype,
		      const_cast<char*>(key.c_str()),
		      const_cast<T*>(&val),
		      0, &_status);
    }

  std::ostringstream o;
  o << "Updating header keyword " << key;
  _checkStatus(o.str());
}

template<class T> void FITSFile::readImage(dm::memimage<T>** image)
{
  // work out fitsio datatype for datatype
  const int fits_datatype = _FITSVal_Datatype( static_cast<T*>(0) );

  // get axis dimensions and create image
  int xw, yw; 
  readKey("NAXIS1", &xw);
  readKey("NAXIS2", &yw);
  *image = new dm::memimage<T>(xw, yw);

  if(_verbose)
    std::cout << "Reading image (" << xw << "x" << yw << ")\n";

  // actually read image
  fits_read_img(_file, fits_datatype, 1, xw*yw, 0,
		&((*image)->flatdata(0)),
		0, &_status);

  _checkStatus("Read image");
}

template<class T> void FITSFile::writeImage(const dm::memimage<T>& image)
{
  const int fits_imagetype = _FITSImg_Datatype(&image);
  const int fits_datatype = _FITSVal_Datatype( static_cast<T*>(0) );
  long axes[2];
  axes[0] = image.xw(); axes[1] = image.yw();

  if(_verbose)
    std::cout << "Writing image (" << image.xw()
	      << "x" << image.yw() << ")\n";

  // try to read existing axis to determine whether an image header exists
  int a1;
  const int def = -999;
  readKey("AXIS1", &a1, &def);

  if( a1 == def )
    {
      // create a new image header if required
      fits_create_img(_file, fits_imagetype, 2, axes, &_status);
      _checkStatus("Writing image header");
    } else {
      // just update the keywords otherwise
      fits_resize_img(_file, fits_imagetype, 2, axes, &_status);
      _checkStatus("Resizing image");
    }

  // write the data
  // const_cast is evil, but cannot get address of real data unless
  // variable is non-const
  // data are not changed
  dm::memimage<T>* img_no_const = const_cast<dm::memimage<T>*>(&image);
  fits_write_img(_file, fits_datatype, 1, image.xw()*image.yw(),
		 &(img_no_const->flatdata(0)), &_status);
  _checkStatus("Writing image");
}

#endif
