/*
Arroyo - software for the simulation of electromagnetic wave propagation
through turbulence and optics.

Copyright (c) 2000-2004 California Institute of Technology.  Written by
Dr. Matthew Britton.  For comments or questions about this software,
please contact the author at mbritton@astro.caltech.edu.

This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as  published by the
Free Software Foundation; either version 2 of the License, or (at your
option) any later version.

This program is provided "as is" and distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  In no
event shall California Institute of Technology be liable to any party
for direct, indirect, special, incidental or consequential damages,
including lost profits, arising out of the use of this software and its
documentation, even if the California Institute of Technology has been
advised of the possibility of such damage.   The California Institute of
Technology has no obligation to provide maintenance, support, updates,
enhancements or modifications.  See the GNU General Public License for
more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
*/

#include <iostream>
#include <iomanip>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>
#include "fitsio.h"
#include "iofits.h"

using namespace std;

namespace Arroyo {

  int iofits::verbose_level = 0;

  iofits::iofits(const char * fname, int iomode){
    fp=NULL; status = 0;
    this->open(fname, iomode);
  } 

  iofits::~iofits(){
    if(fp!=NULL){
      if (fits_close_file(fp, &status)){
	cerr << "iofits error - failed to close file " << filename << endl;
	fits_report_error(stderr, status);
	throw(string("iofits::~iofits"));
      }
    }
  }

  iofits & iofits::operator = (const iofits & iof){
    if(this == &iof)
      return(*this);
    filename = iof.filename; 
    status = iof.status;
    if (fits_open_file(&fp, filename.c_str(), READONLY, &status)) {
      cerr << "iofits op = error - failed to open file " << filename << endl;
      fits_report_error(stderr, status);
      fp = NULL;
      throw(string("iofits op ="));
    }
    // Move to same hdu number in the fits file
    int hdunum, hdutype;
    if(fits_get_hdu_num(iof.fp, &hdunum) || 
       fits_movabs_hdu(fp, hdunum, &hdutype, &status)){
      cerr << "iofit s op = error - could not " << endl;
      fits_report_error(stderr, status);
      throw(string("iofits op ="));
    }
    return(*this);
  }

  void iofits::open(const char * fname, int iomode){
    if(fp!=NULL) close();
    struct stat filestat;
    if(stat(fname, &filestat)!=0){
      cerr << "iofits::open error - file " << fname << " does not exist\n";
      throw(string("iofits::open"));
    }
    if (fits_open_file(&fp, fname, iomode, &status)) {
      cerr << "iofits error - failed to open file " << fname << endl;
      fits_report_error(stderr, status);
      fp = NULL; status = 0;
      throw(string("iofits::open"));
    }
    filename = fname;
  }

  void iofits::close(){
    if(fp!=NULL){
      if (fits_close_file(fp, &status)){
	cerr << "iofits error - failed to close file " << filename << endl;
	fits_report_error(stderr, status);
	throw(string("iofits::close"));
      }
    }
    fp = NULL; status = 0;
  }

  void iofits::create(const char * fname){
    if(fp!=NULL) close();
    struct stat filestat;
    if(stat(fname, &filestat)==0){
      if(verbose_level) cerr << "iofits::create - removing old file " << fname << endl;
      remove(fname);
    }
  
    if(fits_create_file(&fp, fname, &status)){
      cerr << "iofits::create - error creating file " << fname << endl;
      fits_report_error(stderr, status);
      throw(string("iofits::create"));
    }
    filename = fname;
  }

  string iofits::record_head(int record) const {
    string rec = read_record(record);
    if(iofits::verbose_level) 
      cerr << "iofits::record_head: parsing record head from " << rec << endl;
    string record_head(rec, 0, rec.find(' '));
    if(record_head.length() == 9)
      record_head.resize(record_head.length()-1); // clip = sign
    return(record_head);
  }

  string iofits::record_val(int record) const {
    string record_val = read_record(record);
    record_val = string(record_val, 10);
    if(iofits::verbose_level) 
      cerr << "Parsing record value from " << read_record(record) << endl;

    if(record_val[0] == '\''){ 
      record_val.erase(record_val.begin(),record_val.begin()+1);
      record_val.erase(record_val.rfind('\''));
    } else {
      record_val.erase(record_val.rfind('/'));
    }

    if(iofits::verbose_level) 
      cerr << "iofits::record_val - initial value of record_val " 
	   << record << "  x" << record_val << "x" << endl;

    // clear any whitespace at the beginning or end of the string
    while(record_val.length()>0 && record_val[record_val.length()-1] == ' '){
      record_val.erase(record_val.end()-1);
      if(iofits::verbose_level) 
	cerr << "iofits::record_val - trimmed version of record_val " 
	     << record << "  x" << record_val << "x" << endl;
    } 
    while(record_val.length()>0 && record_val[0] == ' '){
      record_val.erase(record_val.begin(), record_val.begin()+1);
      if(iofits::verbose_level) 
	cerr << "iofits::record_val - trimmed version of record_val " 
	     << record << "  x" << record_val << "x" << endl;
    } 
    
    return(record_val);
  }

  int iofits::get_img_dim() const {
    int dimen;
    if(fits_get_img_dim(fp, &dimen, &status)){
      cerr << "iofits::get_img_dim error\n";
      fits_report_error(stderr, status);
      throw(string("iofits::get_img_dim"));
    }
    return(dimen);
  }

  vector<long> iofits::get_img_size() const {
    long ndimen = get_img_dim();
    long * dimens = new long[ndimen];
    if(fits_get_img_size(fp, ndimen, dimens, &status)){
      cerr << "iofits::get_img_size error\n";
      fits_report_error(stderr, status);
      throw(string("iofits::get_img_size"));
    }
    vector<long> vdimens(ndimen);
    for(int i=0; i<ndimen; i++) 
      vdimens[i] = dimens[i];
    delete [] dimens;
    return(vdimens);
  }

  bool iofits::key_exists(const char * keyname) const {
    string value, comment;
    try{read_key(const_cast<char*>(keyname), value, comment);}
    catch(...){
      return false;
    }
    return true;
  }

  void iofits::read_record(int record, string & value) const {
    if(verbose_level) cerr << "iofits::read_record - reading record " 
			   << record << endl;
    if(fits_read_record(fp, record, const_cast<char*>(val), &status)!=0){
      if(verbose_level) 
	cerr << "iofits::read_record - error reading record " << record << endl;
      status = 0;
      throw(string("iofits::read_record"));
    }
    if(verbose_level) cerr << "iofits::read_record - read x" << val << "x" << endl;
    value = val;
  }

  void iofits::delete_key(const char * keyname){
    if(fits_delete_key(fp, const_cast<char*>(keyname), &status)){
      cerr << "iofits::delete_key error - could not delete key " << keyname << endl;
      fits_report_error(stderr, status);
      throw(string("iofits::delete_key"));
    }
  }

  void iofits::delete_record(int record){
    if(fits_delete_record(fp,record,&status)){
      cerr << "iofits::delete_record error - could not delete record " << record << endl;
      fits_report_error(stderr, status);
      throw(string("iofits::delete_record"));
    }
  }

  int iofits::movabs_hdu(const int & hdunum) const {
    int hdutype;
    if(fits_movabs_hdu(fp, hdunum, &hdutype, &status)!=0){
      if(verbose_level){
	cerr << "iofits::movabs_hdu error\n";
	fits_report_error(stderr,status);
      }
      status = 0;
      throw(string("iofits::movabs_hdu"));
    }
    return(hdutype);
  }

  int iofits::movrel_hdu(const int & hdunum) const {
    int hdutype;
    if(fits_movrel_hdu(fp, hdunum, &hdutype, &status)!=0){
      cerr << "iofits::movrel_hdu error\n";
      status = 0;
      throw(string("iofits::movrel_hdu"));
    }
    return(hdutype);
  }

  int iofits::get_num_hdus() const {
    int nhdu;
    fits_get_num_hdus(fp,&nhdu, &status);
    return(nhdu);
  }

  int iofits::get_hdu_num() const {
    int nhdu;
    fits_get_hdu_num(fp,&nhdu);
    return(nhdu);
  }

  void iofits::create_hdu(){
    if(fits_create_hdu(fp,&status)){
      cerr << "iofits::create_hdu error - could not create hdu.  status " << status << endl;
      fits_report_error(stderr, status);
      throw(string("iofits::create_hdu"));
    }
  }

  void iofits::create_image(const vector<long> & axes, imagetype imgtp) {

    if(imgtp==iofits::BYTEIMG) {
      char * t;
      this->create_image(axes, t);
    } else if(imgtp==iofits::USHORTIMG) {
      unsigned short * t;
      create_image(axes, t);
    } else if(imgtp==iofits::SHORTIMG) {
      short * t;
      create_image(axes, t);
    } else if(imgtp==iofits::LONGIMG) {
      long * t;
      create_image(axes, t);
    } else if(imgtp==iofits::FLOATIMG) {
      float * t;
      create_image(axes, t);
    } else if(imgtp==iofits::DOUBLEIMG) {
      double * t;
      create_image(axes, t);
    } else {
      cerr << "iofits::create_image error - creation of images of type "
	   << imgtp << " not currently supported\n";
      throw(string("iofits::create_image"));
    }
  }

  void iofits::read_image_header(int & bitpix, vector<long> & axes) const {
    int maxdim = 32;
    int naxes;
    long * laxes = new long[maxdim];
    if(verbose_level) cerr << "iofits::read_image_header - reading image header\n";
    if(fits_read_imghdr(fp, maxdim, NULL, &bitpix, &naxes, laxes, 
			NULL, NULL, NULL, const_cast<int*>(&status))){
      cerr << "iofits::read_image_header error - could not read image header\n";
      fits_report_error(stderr, status);
      throw(string("iofits::read_image_header"));
    }
  
    axes.resize(naxes);
    for(int i=0; i<naxes; i++)
      axes[i] = laxes[i];

    if(verbose_level) cerr << "iofits::read_image_header - finished reading image header\n";
    delete [] laxes;
  }

  void iofits::write_image_header(int bitpix, const vector<long> & axes) {

    int naxes = axes.size();
    long * laxes = new long[naxes];
    for(int i=0; i<naxes; i++) {
      laxes[i] = axes[i];
    }

    if(verbose_level) cerr << "iofits::write_image_header - writing image header\n";
    if(fits_write_imghdr(fp, bitpix, naxes, laxes, const_cast<int*>(&status))){
      cerr << "iofits::write_image_header error - could not write image header\n";
      fits_report_error(stderr, status);

      cerr << endl;
      this->print_header(cerr, "iofits::write_image_header error");
      cerr << endl << endl;

      throw(string("iofits::write_image_header"));
    }
    if(verbose_level) cerr << "iofits::write_image_header - finished writing image header\n";
    delete [] laxes;
  }

  void iofits::print_header(ostream & os, const char * prefix) const {
    int nkeys;
    string junk;

    fits_get_hdrspace(fp, &nkeys, NULL, &status);

    int vlspc = 30;
  
    for (int i = 1; i <= nkeys; i++)  { 
      this->read_record(i, junk);
      os << setw(vlspc) << prefix << junk << endl;
    }
    os << setw(vlspc) << prefix << "END\n";
  }


  template<>
    void iofits::read_key(const char * keyname, bool & value, string & comment) const {
    int v;
    if(fits_read_key(fp, 
		     TLOGICAL,
		     const_cast<char*>(keyname), 
		     &v,
		     const_cast<char*>(charcomm), 
		     &status)!=0){
      if(verbose_level){
	cerr << "iofits::read_key(bool) - error reading key " << keyname << endl;
	fits_report_error(stderr, status);
      }
      status = 0;
      throw(string("iofits::read_key(bool)"));
    }
    if(v==0) value = false;
    else value = true;
    comment = charcomm;
  }

  template<>
    void iofits::read_key(const char * keyname, string & value, string & comment) const {
    if(verbose_level) cerr << "iofits::read_key(string) - reading keyword " 
			   << keyname << endl;
    if(fits_read_key(fp, 
		     TSTRING,	
		     const_cast<char*>(keyname), 
		     const_cast<char*>(val), 
		     const_cast<char*>(charcomm), 
		     &status)!=0){
      if(verbose_level){
	cerr << "iofits::read_key(string) - error reading key " << keyname << endl;
	fits_report_error(stderr, status);
      }
      status = 0;
      throw(string("iofits::read_key"));
    }
    if(verbose_level) cerr << "iofits::read_key(string) - read x" << val << "x" 
			   << "\tcomment x" << charcomm << "x" << endl;
    comment = charcomm;
    value = val;
  }

  template<>
    void iofits::read_key(const char * keyname, short & value, string & comment) const {
    if(verbose_level) cerr << "iofits::read_key(short) - reading keyword " << 
			keyname << endl;
    if(fits_read_key(fp, 
		     TSHORT,
		     const_cast<char*>(keyname), 
		     &value, 
		     const_cast<char*>(charcomm), 
		     &status)!=0){
      if(verbose_level){
	cerr << "iofits::read_key(long) - error reading key " << keyname << endl;
	fits_report_error(stderr, status);
      }
      status = 0;
      throw(string("iofits::read_key"));
    }
    if(verbose_level) cerr << "iofits::read_key(short) - read: x" << value << "x" 
			   << "\tcomment: x" << charcomm << "x" << endl;
    comment = charcomm;
  }

  template<>
    void iofits::read_key(const char * keyname, unsigned short & value, string & comment) const {
    if(verbose_level) cerr << "iofits::read_key(unsigned short) - reading keyword " << 
			keyname << endl;
    if(fits_read_key(fp, 
		     TUSHORT,
		     const_cast<char*>(keyname), 
		     &value, 
		     const_cast<char*>(charcomm), 
		     &status)!=0){
      if(verbose_level){
	cerr << "iofits::read_key(unsigned short) - error reading key " << keyname << endl;
	fits_report_error(stderr, status);
      }
      status = 0;
      throw(string("iofits::read_key"));
    }
    if(verbose_level) cerr << "iofits::read_key(unsigned short) - read: x" << value << "x" 
			   << "\tcomment: x" << charcomm << "x" << endl;
    comment = charcomm;
  }

  template<>
    void iofits::read_key(const char * keyname, int & value, string & comment) const {
    if(verbose_level) cerr << "iofits::read_key(int) - reading keyword "
			   << keyname << endl;
    if(fits_read_key(fp, 
		     TINT,
		     const_cast<char*>(keyname), 
		     &value, 
		     const_cast<char*>(charcomm), 
		     &status)!=0){
      if(verbose_level){
	cerr << "iofits::read_key(long) - error reading key " << keyname << endl;
	fits_report_error(stderr, status);
      }
      status = 0;
      throw(string("iofits::read_key"));
    }
    if(verbose_level) cerr << "iofits::read_key(int) - read: x" << value << "x" 
			   << "\tcomment: x" << charcomm << "x" << endl;
    comment = charcomm;
  }

  template<>
    void iofits::read_key(const char * keyname, unsigned int & value, string & comment) const {
    if(verbose_level) cerr << "iofits::read_key(unsigned int) - reading keyword " 
			   << keyname << endl;
    if(fits_read_key(fp, 
		     TUINT,
		     const_cast<char*>(keyname), 
		     &value, 
		     const_cast<char*>(charcomm), 
		     &status)!=0){
      if(verbose_level){
	cerr << "iofits::read_key(unsigned int) - error reading key " << keyname << endl;
	fits_report_error(stderr, status);
      }
      status = 0;
      throw(string("iofits::read_key"));
    }
    if(verbose_level) cerr << "iofits::read_key(unsigned int) - read: x" << value << "x" 
			   << "\tcomment: x" << charcomm << "x" << endl;
    comment = charcomm;
  }

  template<>
    void iofits::read_key(const char * keyname, long & value, string & comment) const {
    if(verbose_level) cerr << "iofits::read_key(long) - reading keyword " << 
			keyname << endl;
    if(fits_read_key(fp,
		     TLONG,
		     const_cast<char*>(keyname), 
		     &value, 
		     const_cast<char*>(charcomm), 
		     &status)!=0){
      if(verbose_level){
	cerr << "iofits::read_key(long) - error reading key " << keyname << endl;
	fits_report_error(stderr, status);
      }
      status = 0;
      throw(string("iofits::read_key"));
    }
    if(verbose_level) cerr << "iofits::read_key(long) - read: x" << value << "x" 
			   << "\tcomment: x" << charcomm << "x" << endl;
    comment = charcomm;
  }

  template<>
    void iofits::read_key(const char * keyname, float & value, string & comment) const {
    if(verbose_level) cerr << "iofits::read_key(float) - reading keyword " 
			   << keyname << endl;
    if(fits_read_key(fp, 
		     TFLOAT,
		     const_cast<char*>(keyname), 
		     &value, 
		     const_cast<char*>(charcomm), 
		     &status)!=0){
      if(verbose_level){
	cerr << "iofits::read_key(float) - error reading key " << keyname << endl;
	fits_report_error(stderr, status);
      }
      status = 0;
      throw(string("iofits::read_key"));
    }
    if(verbose_level) cerr << "iofits::read_key(double) - read x" << value << "x"
			   << "\tcomment x" << charcomm << "x" << endl;
    comment = charcomm;
  }

  template<>
    void iofits::read_key(const char * keyname, double & value, string & comment) const {
    if(verbose_level) cerr << "iofits::read_key(double) - reading keyword " 
			   << keyname << endl;
    if(fits_read_key(fp, 
		     TDOUBLE,
		     const_cast<char*>(keyname), 
		     &value, 
		     const_cast<char*>(charcomm), 
		     &status)!=0){
      if(verbose_level){
	cerr << "iofits::read_key(double) - error reading key " << keyname << endl;
	fits_report_error(stderr, status);
      }
      status = 0;
      throw(string("iofits::read_key"));
    }
    if(verbose_level) cerr << "iofits::read_key(double) - read x" << value << "x"
			   << "\tcomment x" << charcomm << "x" << endl;
    comment = charcomm;
  }

  template<>
    void iofits::write_key(const char * keyname, const bool & value, const string & comment) {

    int v = 1;
    if(value==false) v = 0;

    if(verbose_level) cerr << "iofits::write_key(logical) - writing keyword " 
			   << keyname << " value x" << v << "x comment x" << comment << "x" << endl;
    if(fits_write_key(fp, 
		      TLOGICAL,
		      const_cast<char*>(keyname), 
		      (void *)&v,
		      const_cast<char*>(comment.c_str()), 
		      const_cast<int*>(&status))!=0){
      if(verbose_level){
	cerr << "iofits::write_key(bool) - error writing key " << keyname << endl;
	fits_report_error(stderr, status);
      }
      throw(string("iofits::write_key"));
    }
    if(verbose_level) cerr << "iofits::write_key(logical) - wrote x" << value << "x"
			   << "\tcomment x" << comment << "x"<< endl;
  }

  template<>
    void iofits::write_key(const char * keyname, const string & value, const string & comment) {
    if(verbose_level) cerr << "iofits::write_key(string) - writing keyword " 
			   << keyname << " value x" << value << "x comment x" << comment << "x" << endl;

    if(fits_write_key(fp, 
		      TSTRING,
		      const_cast<char*>(keyname), 
		      (void *)(const_cast<char*>(value.c_str())),
		      const_cast<char*>(comment.c_str()), 
		      const_cast<int*>(&status))!=0){
      if(verbose_level){
	cerr << "iofits::write_key(string) - error writing key " << keyname << endl;
	fits_report_error(stderr, status);
      }
      throw(string("iofits::write_key"));
    }
    if(verbose_level) cerr << "iofits::write_key(string) - wrote x" << value << "x"
			   << "\tcomment x" << comment << "x"<< endl;
  }

  template<>
    void iofits::write_key(const char * keyname, const short & value, const string & comment) {
    if(verbose_level) cerr << "iofits::write_key(short) - writing keyword " 
			   << keyname << " value x" << value << "x comment x" << comment << "x" << endl;
    if(fits_write_key(fp, 
		      TSHORT,
		      const_cast<char*>(keyname), 
		      (void *)&value,
		      const_cast<char*>(comment.c_str()), 
		      const_cast<int*>(&status))!=0){
      if(verbose_level){
	cerr << "iofits::write_key(short) - error writing key " << keyname << endl;
	fits_report_error(stderr, status);
      }
      throw(string("iofits::write_key"));
    }
    if(verbose_level) cerr << "iofits::write_key(short) - wrote x" << value << "x"
			   << "\tcomment x" << comment << "x" << endl;
  }

  template<>
    void iofits::write_key(const char * keyname, const unsigned short & value, const string & comment) {
    if(verbose_level) cerr << "iofits::write_key(unsigned short) - writing keyword " 
			   << keyname << " value x" << value << "x comment x" << comment << "x" << endl;
    if(fits_write_key(fp, 
		      TUSHORT,
		      const_cast<char*>(keyname), 
		      (void *)&value,
		      const_cast<char*>(comment.c_str()), 
		      const_cast<int*>(&status))!=0){
      if(verbose_level){
	cerr << "iofits::write_key(unsigned short) - error writing key " << keyname << endl;
	fits_report_error(stderr, status);
      }
      throw(string("iofits::write_key"));
    }
    if(verbose_level) cerr << "iofits::write_key(unsigned short) - wrote x" << value << "x"
			   << "\tcomment x" << comment << "x" << endl;
  }

  template<>
    void iofits::write_key(const char * keyname, const int & value, const string & comment) {
    if(verbose_level) cerr << "iofits::write_key(int) - writing keyword " 
			   << keyname << " value x" << value << "x comment x" << comment << "x" << endl;
    if(fits_write_key(fp, 
		      TINT,
		      const_cast<char*>(keyname), 
		      (void *)&value,
		      const_cast<char*>(comment.c_str()), 
		      const_cast<int*>(&status))!=0){
      if(verbose_level){
	cerr << "iofits::write_key(int) - error writing key " << keyname << endl;
	fits_report_error(stderr, status);
      }
      throw(string("iofits::write_key"));
    }
    if(verbose_level) cerr << "iofits::write_key(int) - wrote x" << value << "x"
			   << "\tcomment x" << comment << "x" << endl;
  }

  template<>
    void iofits::write_key(const char * keyname, const unsigned int & value, const string & comment) {
    if(verbose_level) cerr << "iofits::write_key(unsigned int) - writing keyword " 
			   << keyname << " value x" << value << "x comment x" << comment << "x" << endl;
    if(fits_write_key(fp, 
		      TUINT,
		      const_cast<char*>(keyname), 
		      (void *)&value,
		      const_cast<char*>(comment.c_str()), 
		      const_cast<int*>(&status))!=0){
      if(verbose_level){
	cerr << "iofits::write_key(unsigned int) - error writing key " << keyname << endl;
	fits_report_error(stderr, status);
      }
      throw(string("iofits::write_key"));
    }
    if(verbose_level) cerr << "iofits::write_key(unsigned int) - wrote x" << value << "x"
			   << "\tcomment x" << comment << "x" << endl;
  }

  template<>
    void iofits::write_key(const char * keyname, const long & value, const string & comment) {
    if(verbose_level) cerr << "iofits::write_key(long) - writing keyword " 
			   << keyname << " value x" << value << "x comment x" << comment << "x" << endl;
    if(fits_write_key(fp, 
		      TLONG,
		      const_cast<char*>(keyname), 
		      (void *)&value,
		      const_cast<char*>(comment.c_str()), 
		      const_cast<int*>(&status))!=0){
      if(verbose_level){
	cerr << "iofits::write_key(long) - error writing key " << keyname << endl;
	fits_report_error(stderr, status);
      }
      throw(string("iofits::write_key"));
    }
    if(verbose_level) cerr << "iofits::write_key(long) - wrote x" << value << "x"
			   << "\tcomment x" << comment << "x" << endl;
  }

  template<>
    void iofits::write_key(const char * keyname, const float & value, const string & comment) {
    if(verbose_level) cerr << "iofits::write_key(float) - writing keyword " 
			   << keyname << " value x" << value << "x comment x" << comment << "x" << endl;
    if(fits_write_key(fp, 
		      TFLOAT,
		      const_cast<char*>(keyname), 
		      (void *)&value, 
		      const_cast<char*>(comment.c_str()), 
		      const_cast<int*>(&status))!=0){
      if(verbose_level){
	cerr << "iofits::write_key(float) - error writing key " << keyname << endl;
	fits_report_error(stderr, status);
      }
      throw(string("iofits::write_key"));
    }
    if(verbose_level) cerr << "iofits::write_key(float) - wrote x" << value << "x"
			   << "\tcomment x" << comment << "x" << endl;
  }

  template<>
    void iofits::write_key(const char * keyname, const double & value, const string & comment) {
    if(verbose_level) cerr << "iofits::write_key(double) - writing keyword " 
			   << keyname << " value x" << value << "x comment x" << comment << "x" << endl;
    if(fits_write_key(fp, 
		      TDOUBLE,
		      const_cast<char*>(keyname), 
		      (void *)&value, 
		      const_cast<char*>(comment.c_str()), 
		      const_cast<int*>(&status))!=0){
      if(verbose_level){
	cerr << "iofits::write_key(double) - error writing key " << keyname << endl;
	fits_report_error(stderr, status);
      }
      throw(string("iofits::write_key"));
    }
    if(verbose_level) cerr << "iofits::write_key(double) - wrote x" << value << "x"
			   << "\tcomment x" << comment << "x" << endl;
  }

  template<>
    void iofits::update_key(const char * keyname, const bool & value, const string & comment) {

    int v = 1;
    if(value==false) v = 0;

    if(verbose_level) cerr << "iofits::update_key(logical) - writing keyword " 
			   << keyname << " value x" << v << "x comment x" << comment << "x" << endl;
    if(fits_update_key_log(fp, const_cast<char*>(keyname), v,
			  const_cast<char*>(comment.c_str()), const_cast<int*>(&status))!=0){
      if(verbose_level){
	cerr << "iofits::update_key(bool) - error writing key " << keyname << endl;
	fits_report_error(stderr, status);
      }
      throw(string("iofits::update_key"));
    }
    if(verbose_level) cerr << "iofits::update_key(logical) - wrote x" << value << "x"
			   << "\tcomment x" << comment << "x"<< endl;
  }

  template<>
    void iofits::update_key(const char * keyname, const string & value, const string & comment) {
    if(verbose_level) cerr << "iofits::update_key(string) - writing keyword " 
			   << keyname << " value x" << value << "x comment x" << comment << "x" << endl;
    if(fits_update_key_str(fp, const_cast<char*>(keyname), const_cast<char*>(value.c_str()),
			  const_cast<char*>(comment.c_str()), const_cast<int*>(&status))!=0){
      if(verbose_level){
	cerr << "iofits::update_key(string) - error writing key " << keyname << endl;
	fits_report_error(stderr, status);
      }
      throw(string("iofits::update_key"));
    }
    if(verbose_level) cerr << "iofits::update_key(string) - wrote x" << value << "x"
			   << "\tcomment x" << comment << "x"<< endl;
  }

  template<>
    void iofits::update_key(const char * keyname, const long & value, const string & comment) {
    if(verbose_level) cerr << "iofits::update_key(long) - writing keyword " 
			   << keyname << " value x" << value << "x comment x" << comment << "x" << endl;
    if(fits_update_key_lng(fp, const_cast<char*>(keyname), value,
			  const_cast<char*>(comment.c_str()), const_cast<int*>(&status))!=0){
      if(verbose_level){
	cerr << "iofits::update_key(long) - error writing key " << keyname << endl;
	fits_report_error(stderr, status);
      }
      throw(string("iofits::update_key"));
    }
    if(verbose_level) cerr << "iofits::update_key(long) - wrote x" << value << "x"
			   << "\tcomment x" << comment << "x" << endl;
  }

  template<>
    void iofits::update_key(const char * keyname, const double & value, const string & comment) {
    if(verbose_level) cerr << "iofits::update_key(double) - writing keyword " 
			   << keyname << " value x" << value << "x comment x" << comment << "x" << endl;
    if(fits_update_key_dbl(fp, const_cast<char*>(keyname), value, 8,
			  const_cast<char*>(comment.c_str()), const_cast<int*>(&status))!=0){
      if(verbose_level){
	cerr << "iofits::update_key(double) - error writing key " << keyname << endl;
	fits_report_error(stderr, status);
      }
      throw(string("iofits::update_key"));
    }
    if(verbose_level) cerr << "iofits::update_key(double) - wrote x" << value << "x"
			   << "\tcomment x" << comment << "x" << endl;
  }

  template<>
    void iofits::create_image(const vector<long> & axes, const char * pixframe) {
    
    int anynull;
    if(verbose_level) cerr << "iofits::create_image(char) - creating image with " 
			   << axes.size() << " axes\n";
    long * longaxes = new long[axes.size()];
    for(int i=0; i<axes.size(); i++) longaxes[i] = axes[i];
    if(iofits::verbose_level){
      cerr << "iofits::create_image - creating image with " << axes.size() << endl;
      for(int i=0; i<axes.size(); i++)
	cerr << "\taxis " << i << "\t" << longaxes[i] << endl;
    }
    if(fits_create_img(fp, BYTEIMG, axes.size(), longaxes, const_cast<int*>(&status))){
      cerr << "iofits::create_image(char) error - could not create image\n";
      fits_report_error(stderr, status);
      throw(string("iofits::create_image"));
    }
    delete [] longaxes;
    if(verbose_level) cerr << "iofits::create_image(char) - finished creating image\n";
  }

  template<>
    void iofits::create_image(const vector<long> & axes, const unsigned short * pixframe) {

    int anynull;
    if(verbose_level) cerr << "iofits::create_image(unsigned short) - creating image with " 
			   << axes.size() << " axes\n";
    long * longaxes = new long[axes.size()];
    for(int i=0; i<axes.size(); i++) longaxes[i] = axes[i];
    if(iofits::verbose_level){
      cerr << "iofits::create_image - creating image with " << axes.size() << endl;
      for(int i=0; i<axes.size(); i++)
	cerr << "\taxis " << i << "\t" << longaxes[i] << endl;
    }
    if(fits_create_img(fp, USHORT_IMG, axes.size(), longaxes, const_cast<int*>(&status))){
      cerr << "iofits::create_image(short) error - could not create image\n";
      fits_report_error(stderr, status);
      throw(string("iofits::create_image"));
    }
    delete [] longaxes;
    if(verbose_level) cerr << "iofits::create_image(short) - finished creating image\n";
  }

  template<>
    void iofits::create_image(const vector<long> & axes, const short * pixframe) {

    int anynull;
    if(verbose_level) cerr << "iofits::create_image(short) - creating image with " 
			   << axes.size() << " axes\n";
    long * longaxes = new long[axes.size()];
    for(int i=0; i<axes.size(); i++) longaxes[i] = axes[i];
    if(iofits::verbose_level){
      cerr << "iofits::create_image - creating image with " << axes.size() << endl;
      for(int i=0; i<axes.size(); i++)
	cerr << "\taxis " << i << "\t" << longaxes[i] << endl;
    }
    if(fits_create_img(fp, SHORTIMG, axes.size(), longaxes, const_cast<int*>(&status))){
      cerr << "iofits::create_image(short) error - could not create image\n";
      fits_report_error(stderr, status);
      throw(string("iofits::create_image"));
    }
    delete [] longaxes;
    if(verbose_level) cerr << "iofits::create_image(short) - finished creating image\n";
  }

  template<>
    void iofits::create_image(const vector<long> & axes, const long * pixframe) {

    int anynull;
    if(verbose_level) cerr << "iofits::create_image(long) - creating image with "
			   << axes.size() << " axes\n";
    long * longaxes = new long[axes.size()];
    for(int i=0; i<axes.size(); i++) longaxes[i] = axes[i];
    if(iofits::verbose_level){
      cerr << "iofits::create_image - creating image with " << axes.size() << endl;
      for(int i=0; i<axes.size(); i++)
	cerr << "\taxis " << i << "\t" << longaxes[i] << endl;
    }
  
    if(fits_create_img(fp, LONGIMG, axes.size(), longaxes, const_cast<int*>(&status))){
      cerr << "iofits::create_image(long) error - could not create image\n";
      fits_report_error(stderr, status);
      throw(string("iofits::create_image"));
    }
    delete [] longaxes;  
    if(verbose_level) cerr << "iofits::create_image(long) - finished creating image\n";
  }

  template<>
    void iofits::create_image(const vector<long> & axes, const float * pixframe) {

    int anynull;
    if(verbose_level) cerr << "iofits::create_image(float) - creating image with "
			   << axes.size() << " axes\n";
    long * longaxes = new long[axes.size()];
    for(int i=0; i<axes.size(); i++) longaxes[i] = axes[i];
    if(iofits::verbose_level){
      cerr << "iofits::create_image - creating image with " << axes.size() << endl;
      for(int i=0; i<axes.size(); i++)
	cerr << "\taxis " << i << "\t" << longaxes[i] << endl;
    }
    if(fits_create_img(fp, FLOATIMG, axes.size(), longaxes, const_cast<int*>(&status))){
      cerr << "iofits::create_image(float) error - could not create image\n";
      fits_report_error(stderr, status);
      throw(string("iofits::create_image"));
    }
    delete [] longaxes;
    if(verbose_level) cerr << "iofits::create_image(float) - finished creating image\n";
  }

  template<>
    void iofits::create_image(const vector<long> & axes, const double * pixframe) {

    int anynull;
    if(verbose_level) cerr << "iofits::create_image(double) - creating image with " 
			   << axes.size() << " axes\n";
    long * longaxes = new long[axes.size()];
    for(int i=0; i<axes.size(); i++) longaxes[i] = axes[i];
    if(iofits::verbose_level){
      cerr << "iofits::create_image - creating image with " << axes.size() << endl;
      for(int i=0; i<axes.size(); i++)
	cerr << "\taxis " << i << "\t" << longaxes[i] << endl;
    }
    if(fits_create_img(fp, DOUBLEIMG, axes.size(), longaxes, const_cast<int*>(&status))){
      cerr << "iofits::create_image(double) error - could not create image\n";
      fits_report_error(stderr, status);
      throw(string("iofits::create_image"));
    }
    delete [] longaxes;
    if(verbose_level) cerr << "iofits::create_image(double) - finished creating image\n";
  }

  template<>
    void iofits::read_image(int first, int last, char * pixframe) const {
    int anynull;
    if(verbose_level) cerr << "iofits::read_image(char) - reading image\n";
    if(fits_read_img_byt(fp, 1, first+1,last+1, 0,
			 (unsigned char*)pixframe, &anynull, &status)){
      cerr << "iofits::read_image(char) error - could not read image\n";
      fits_report_error(stderr, status);
      throw(string("iofits::read_image"));
    }
    if(verbose_level) cerr << "iofits::read_image(char) - finished reading image\n";
  }

  template<>
    void iofits::read_image(int first, int last, unsigned short * pixframe) const {
    int anynull;
    if(verbose_level) cerr << "iofits::read_image(unsigned short) - reading image\n";
    if(fits_read_img_usht(fp, 1, first+1,last+1, 0,
			  pixframe, &anynull, &status)){
      cerr << "iofits::read_image(unsigned short) error - could not read image\n";
      fits_report_error(stderr, status);
      throw(string("iofits::read_image"));
    }
    if(verbose_level) cerr << "iofits::read_image(unsigned short) - finished reading image\n";
  }

  template<>
    void iofits::read_image(int first, int last, short * pixframe) const {
    int anynull;
    if(verbose_level) cerr << "iofits::read_image(short) - reading image\n";
    if(fits_read_img_sht(fp, 1, first+1,last+1, 0,
			 pixframe, &anynull, &status)){
      cerr << "iofits::read_image(short) error - could not read image\n";
      fits_report_error(stderr, status);
      throw(string("iofits::read_image"));
    }
    if(verbose_level) cerr << "iofits::read_image(short) - finished reading image\n";
  }

  template<>
    void iofits::read_image(int first, int last, long * pixframe) const {
    int anynull;
    if(verbose_level) cerr << "iofits::read_image(long) - reading image\n";
    if(fits_read_img_lng(fp, 1, first+1,last+1, 0,
			 pixframe, &anynull, &status)){
      cerr << "iofits::read_image(long) error - could not read image\n";
      fits_report_error(stderr, status);
      throw(string("iofits::read_image"));
    }
    if(verbose_level) cerr << "iofits::read_image(long) - finished reading image\n";
  }

  template<>
    void iofits::read_image(int first, int last, float * pixframe) const {
    int anynull;
    if(verbose_level) cerr << "iofits::read_image(float) - reading image\n";
    if(fits_read_img_flt(fp, 1, first+1,last+1, 0, 
			 pixframe, &anynull, &status)){
      cerr << "iofits::read_image(float) error - could not read image\n";
      fits_report_error(stderr, status);
      throw(string("iofits::read_image"));
    }
    if(verbose_level) cerr << "iofits::read_image(float) - finished reading image\n";
  }

  template<>
    void iofits::read_image(int first, int last, double * pixframe) const {
    int anynull;
    if(verbose_level) cerr << "iofits::read_image(double) - reading image\n";
    if(fits_read_img_dbl(fp, 1, first+1,last+1, 0,
			 pixframe, &anynull, &status)){
      cerr << "iofits::read_image(double) error - could not read image\n";
      fits_report_error(stderr, status);
      throw(string("iofits::read_image"));
    }
    if(verbose_level) cerr << "iofits::read_image(double) - finished reading image\n";
  }

  template<>
    void iofits::write_image(int first, int last, const char * pixframe) {
    if(verbose_level) cerr << "iofits::write_image(byte) - writing image with "
			   << last-first+1 << " elements\n";
    if(fits_write_img_byt(fp, 1, first+1, last+1, (unsigned char*) (const_cast<char*>(pixframe)), const_cast<int*>(&status))){
      cerr << "iofits::write_image(byte) error - could not write image\n";
      fits_report_error(stderr, status);
      throw(string("iofits::write_image"));
    }
    if(verbose_level) cerr << "iofits::write_image(byte) - finished writing image\n";
  }

  template<>
    void iofits::write_image(int first, int last, const unsigned short * pixframe) {
    if(verbose_level) cerr << "iofits::write_image(unsigned short) - writing image with "
			   << last-first+1 << " elements\n";
    if(fits_write_img_usht(fp, 1, first+1, last+1, const_cast<unsigned short*>(pixframe), const_cast<int*>(&status))){
      cerr << "iofits::write_image(unsigned short) error - could not write image\n";
      fits_report_error(stderr, status);
      throw(string("iofits::write_image"));
    }
    if(verbose_level) cerr << "iofits::write_image(unsigned short) - finished writing image\n";
  }

  template<>
    void iofits::write_image(int first, int last, const short * pixframe) {
    if(verbose_level) cerr << "iofits::write_image(short) - writing image with "
			   << last-first+1 << " elements\n";
    if(fits_write_img_sht(fp, 1, first+1, last+1, const_cast<short*>(pixframe), const_cast<int*>(&status))){
      cerr << "iofits::write_image(short) error - could not write image\n";
      fits_report_error(stderr, status);
      throw(string("iofits::write_image"));
    }
    if(verbose_level) cerr << "iofits::write_image(short) - finished writing image\n";
  }

  template<>
    void iofits::write_image(int first, int last, const long * pixframe) {
    if(verbose_level) cerr << "iofits::write_image(long) - writing image\n"
			   << last-first+1 << " elements\n";
    if(fits_write_img_lng(fp, 1, first+1, last+1, const_cast<long*>(pixframe), const_cast<int*>(&status))){
      cerr << "iofits::write_image(long) error - could not write image\n";
      fits_report_error(stderr, status);
      throw(string("iofits::write_image"));
    }
    if(verbose_level) cerr << "iofits::write_image(long) - finished writing image\n";
  }

  template<>
    void iofits::write_image(int first, int last, const float * pixframe) {
    if(verbose_level) cerr << "iofits::write_image(float) - writing image\n"
			   << last-first+1 << " elements\n";
    if(fits_write_img_flt(fp, 1, first+1, last+1, const_cast<float*>(pixframe), const_cast<int*>(&status))){
      cerr << "iofits::write_image(float) error - could not write image\n";
      fits_report_error(stderr, status);
      throw(string("iofits::write_image"));
    }
    if(verbose_level) cerr << "iofits::write_image(float) - finished writing image\n";
  }

  template<>
    void iofits::write_image(int first, int last, const double * pixframe) {
    if(verbose_level) cerr << "iofits::write_image(double) - writing image with "
			   << last-first+1 << " elements\n";
    if(fits_write_img_dbl(fp, 1, first+1, last+1, const_cast<double*>(pixframe), const_cast<int*>(&status))){
      cerr << "iofits::write_image(double) error - could not write image\n";
      fits_report_error(stderr, status);
      throw(string("iofits::write_image"));
    }
    if(verbose_level) cerr << "iofits::write_image(double) - finished writing image\n";
  }



}

