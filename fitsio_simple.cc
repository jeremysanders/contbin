#include <iostream>
#include <set>

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include "fitsio_simple.hh"

FITSFile::FITSFile()
  : _file(0), _status(0), _verbose(true)
{
}

FITSFile::FITSFile(const std::string& filename, OpenMode mode)
  : _file(0), _status(0), _verbose(true)
{
  open(filename, mode);
}

FITSFile::~FITSFile()
{
  if(_file != 0)
    close();
}

void FITSFile::_checkStatus(const std::string& operation)
{
  if(_status != 0)
    {
      // dump out fitsio errors if there is one
      std::cerr << "FITSIO operation failed when:\n "
		<< operation
		<< "\nFilename:\n "
		<< _filename
		<< '\n';
      std::cerr.flush();

      fits_report_error(stderr, _status);

      // quit the program
      exit(1);
    }
}

void FITSFile::open(const std::string& filename, OpenMode mode)
{
  _file = 0;
  _status = 0;
  _filename = filename;

  switch(mode)
    {
    case RO:
      if(_verbose)
	std::cout << "Opening " << filename
		  << " (RO)" << std::endl;
      fits_open_file(&_file, filename.c_str(), READONLY, &_status);
      _checkStatus("Opening file (RO)");
      break;
      
    case RW:
      if(_verbose)
	std::cout << "Opening " << filename
		  << " (RW)" << std::endl;
      fits_open_file(&_file, filename.c_str(), READWRITE, &_status);
      _checkStatus("Opening file (RW)");
      break;
      
    case Create:
      if(_verbose)
	std::cout << "Creating " << filename
		  << std::endl;
      
      unlink(filename.c_str());
      fits_create_file(&_file, filename.c_str(), &_status);
      _checkStatus("Creating file");

      break;
    }

}

void FITSFile::close()
{
  if(_file != 0)
    {
      if(_verbose)
	std::cout << "Closing " << _filename
		  << std::endl;

      fits_close_file(_file, &_status);
      _checkStatus("Closing file");

      _file = 0;
      _status = 0;
      _filename.clear();
    };
}

void FITSFile::readKey(const std::string& key, std::string* val,
		       const std::string* defaultval,
		       std::string* comment)
{
  char commentbuf[512];
  char* buffer = 0;
  const int retval = fits_read_key_longstr(_file,
					   const_cast<char*>(key.c_str()),
					   &buffer,
					   commentbuf, &_status);

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
  else
    {
      *val = buffer;
    }
  
  if(buffer != 0)
    free(buffer);

  // keep comment if wanted
  if(comment != 0)
    {
      *comment = commentbuf;
    }
}

void FITSFile::updateKey(const std::string& key, const std::string& val,
			 const std::string* comment)
{
  if(comment != 0)
    {
      fits_update_key_longstr(_file,
			      const_cast<char*>(key.c_str()),
			      const_cast<char*>(val.c_str()),
			      const_cast<char*>(comment->c_str()), &_status);
    }
  else
    {
      fits_update_key_longstr(_file,
			      const_cast<char*>(key.c_str()),
			      const_cast<char*>(val.c_str()),
			      0, &_status);
    }

  std::ostringstream o;
  o << "Updating header keyword " << key;
  _checkStatus(o.str());
}

void FITSFile::writeComment(const std::string& comment)
{
  fits_write_comment(_file, comment.c_str(), &_status);
  _checkStatus("Writing comment");
}

void FITSFile::writeHistory(const std::string& history)
{
  fits_write_history(_file, history.c_str(), &_status);
  _checkStatus("Writing history");
}

void FITSFile::writeDate()
{
  fits_write_date(_file, &_status);
  _checkStatus("Writing date");
}

// list of exclusions for copying headers
static const int _hdr_exclude_num = 26;
static const char* const _hdr_exclude[_hdr_exclude_num] =
  {
    "NAXIS1", "NAXIS2", "NAXIS3", "NAXIS4", "NAXIS5",
    "BZERO", "BSCALE", "BUNIT",
    "SIMPLE", "BITPIX", "NAXIS", "EXTEND", "XTENSION", "PCOUNT", 
    "GCOUNT", "TFIELDS", "TTYPE", "TBCOL", "TFORM", "TUNIT", "THEAP",
    "TDIM", "GROUPS", "DATASUM", "CHECKSUM", "END"
  };

// copy fits header to another file, avoiding keywords already set
void FITSFile::copyHeaderTo(FITSFile& other)
{
  // make a set to test for existing headers
  std::set<std::string> dontcopy(_hdr_exclude, _hdr_exclude+_hdr_exclude_num);

  // get number of keywords in this header
  int keysexist, morekeys;
  fits_get_hdrspace(_file, &keysexist, &morekeys, &_status);
  _checkStatus("Get number of keywords");

  // iterate over each of the keywords in the header
  for(int i=1; i<=keysexist; ++i)
    {
      // read the record
      char card[81];
      fits_read_record(_file, i, card, &_status);
      _checkStatus("Read keyword");

      // get the keyname from the record
      char keyname[81];
      int length;
      fits_get_keyname(card, keyname, &length, &_status);

      // copy keyword if allowed
      if( dontcopy.find( std::string(keyname) ) == dontcopy.end() )
	{
	  fits_write_record(other._file, card, &other._status);
	  other._checkStatus("Writing keyword record");
	}
    }
}

void FITSFile::writeDatestamp(const std::string& program)
{
  char date[64];
  int timeref = 0;
  fits_get_system_time(date, &timeref, &_status);
  _checkStatus("Get date");

  std::ostringstream o;
  o << program << ' ' << date;
  writeHistory(o.str());
}

// int main()
// {
//   FITSFile a("/data/jss/chandra/per_mega/binning/temp_bg_binned.fits",
// 	     FITSFile::RW);

//   const int defval = 101;
//   std::string i;
//   std::string comment;
//   a.readKey("NAXIS1", &i);
//   std::cout << i << '\n';
//   std::cout << comment << '\n';

//   a.updateKey("FOO", "foo bar", &comment);

//   dm::memimage<unsigned char>* img;
//   a.readImage(&img);

//   std::cout << (*img)(100, 100) << ' '
// 	    << (*img)(1460, 1619) << ' '
// 	    << '\n';

//   *img += 10;

//   FITSFile b("foo.fits", FITSFile::RW);
//   b.writeImage(*img);
//   //a.copyHeaderTo(b);
// }
