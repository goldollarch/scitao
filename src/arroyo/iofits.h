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

#ifndef IOFITS_H
#define IOFITS_H

#include <iostream>
#include <string>
#include <vector>
#include "AO_cpp.h"
#include "fitsio.h"

namespace Arroyo {

  using std::string;
  using std::vector;
  using std::ostream;
  using std::cerr;
  using std::endl;

  ///
  /// A class to wrap cfitsio function calls
  ///

  class iofits {

  protected:

    /// pointer to the fits file underlying this instance
    fitsfile * fp;

    /// name of the fits file
    string filename;

    /// mutable status variable used in cfitsio function calls
    mutable int status;

  private:

    /// char array used to store temporarily comments from header keywords
    char charcomm[256];

    /// char array used to store temporarily header keywords parsed as chars
    char val[256];

  public:

    enum imagetype { FLOATIMG=-32, USHORTIMG=20, DOUBLEIMG=-64, BYTEIMG=8,
    					SHORTIMG=16, LONGIMG=32 };

    ///////////////////////////////////////////
    ///  Null constructor
    iofits(){fp = NULL; status = 0;};

    ///////////////////////////////////////////
    ///  Construct iofits object from file
    iofits(const char * infile, int iomode = READONLY);

    ///////////////////////////////////////////
    ///  Copy constructor
    iofits(const iofits & iof){iofits::operator=(iof);};

    ///////////////////////////////////////////
    ///  Destructor that closes the file
    virtual ~iofits();

    ///////////////////////////////////////////
    ///  Operator =  : creates another fitsfile pointer to file
    iofits & operator = (const iofits & iof);

    ///////////////////////////////////////////
    ///  Open file
    void open(const char * infile, int iomode = READONLY);

    ///////////////////////////////////////////
    ///  Close file
    void close();

    ///////////////////////////////////////////
    ///  Create file
    void create(const char * infile);

    ///////////////////////////////////////////
    ///  Get name
    string name() const {return(filename);}

    ///////////////////////////////////////////
    ///  Function to check whether keyword exists
    bool key_exists(const char * keyname) const;

    ///////////////////////////////////////////
    ///  fits_read_record
    void read_record(int record, 
		     string & value) const;

    ///////////////////////////////////////////
    ///  templatized fits_read_key 
    template<class T>
      void read_key(const char * keyname, 
		    T & value, 
		    string & comment) const;
    
    ///////////////////////////////////////////
    ///  templatized fits_write_key 
    template<class T>
      void write_key(const char * keyname, 
		     const T & value,
		     const string & comment); 

    ///////////////////////////////////////////
    ///  templatized fits_update_key 
    template<class T>
      void update_key(const char * keyname, 
		      const T & value,
		      const string & comment); 

   ///////////////////////////////////////////
    ///  fits_delete_key 
    void delete_key(const char * keyname);  

    ///////////////////////////////////////////
    ///  fits_delete_record
    void delete_record(int record);  

    ///////////////////////////////////////////
    ///  Returns the record of entry record
    string read_record(int record) const {
      char card[FLEN_CARD]; 
      fits_read_record(fp, record, card, &status);
      return(string(card));
    };

    ///////////////////////////////////////////
    ///  Function to get the record head of entry record
    string record_head(int record) const;

    ///////////////////////////////////////////
    ///  Function to get the record value of entry record
    string record_val(int record) const;
  
    ///////////////////////////////////////////
    ///  fits_get_img_dim
    int get_img_dim() const;

    ///////////////////////////////////////////
    ///  Return the status
    int get_status() const {return(status);};

    ///////////////////////////////////////////
    ///  fits_get_img_size
    vector<long> get_img_size() const;

    ///////////////////////////////////////////
    ///  fits_movabs_hdu
    int movabs_hdu(const int & hdunum) const;

    ///////////////////////////////////////////
    ///  fits_movrel_hdu
    int movrel_hdu(const int & hdunum) const;

    ///////////////////////////////////////////
    ///  fits_get_num_hdus
    int get_num_hdus() const;

    ///////////////////////////////////////////
    ///  fits_get_hdu_num
    int get_hdu_num() const;

    ///////////////////////////////////////////
    ///  Create header
    void create_hdu();

    ///////////////////////////////////////////
    ///  fits_create_image for doubles
    void create_image(const vector<long> & axes, imagetype imgtp);

    ///////////////////////////////////////////
    ///  templatized fits_create_image 
    template<class T> 
      void create_image(const vector<long> & axes, const T * t);

    ///////////////////////////////////////////
    ///  templatized fits_read_img 
    template<class T>
      void read_image(int first, int last, T * pixframe) const;

    ///////////////////////////////////////////
    ///  templatized fits_read_img 
    template<class T>
      void write_image(int first, int last, const T * pixframe);

    ///////////////////////////////////////////
    ///  fits_read_imghdr
    void read_image_header(int & bitpix, vector<long> & axes) const;

    ///////////////////////////////////////////
    ///  fits_write_imghdr
    void write_image_header(int bitpix, const vector<long> & axes);

    ///////////////////////////////////////////
    ///  Write out all the items in the header
    void print_header(ostream & os, const char * prefix = "") const;

    /// A verbose_level for printing messages
    static int verbose_level;
  
  };

  template<> void iofits::read_key(const char * keyname,
				bool & value, string & comment) const;
  template<> void iofits::read_key(const char * keyname,
				string & value, string & comment) const;
  template<> void iofits::read_key(const char * keyname,
				long & value, string & comment) const;
  template<> void iofits::read_key(const char * keyname,
				float & value, string & comment) const;
  template<> void iofits::read_key(const char * keyname,
				double & value, string & comment) const;

  template<> void iofits::write_key(const char * keyname,
				const bool & value, const string & comment);
  template<> void iofits::write_key(const char * keyname,
				const string & value, const string & comment);
  template<> void iofits::write_key(const char * keyname,
				const long & value, const string & comment);
  template<> void iofits::write_key(const char * keyname,
				const float & value, const string & comment);
  template<> void iofits::write_key(const char * keyname,
				const double & value, const string & comment);

  template<> void iofits::create_image(const vector<long> & axes,
				const char * pixframe);
  template<> void iofits::create_image(const vector<long> & axes,
				const unsigned short * pixframe);
  template<> void iofits::create_image(const vector<long> & axes,
				const short * pixframe);
  template<> void iofits::create_image(const vector<long> & axes,
				const long * pixframe);
  template<> void iofits::create_image(const vector<long> & axes,
				const float * pixframe);
  template<> void iofits::create_image(const vector<long> & axes,
				const double * pixframe);

  template<> void iofits::read_image(int first, int last, char * pixframe) const;
  template<> void iofits::read_image(int first, int last, unsigned short * pixframe) const;
  template<> void iofits::read_image(int first, int last, short * pixframe) const;
  template<> void iofits::read_image(int first, int last, long * pixframe) const;
  template<> void iofits::read_image(int first, int last, float * pixframe) const;
  template<> void iofits::read_image(int first, int last, double * pixframe) const;

  template<> void iofits::write_image(int first, int last, const char * pixframe);
  template<> void iofits::write_image(int first, int last, const unsigned short * pixframe);
  template<> void iofits::write_image(int first, int last, const short * pixframe);
  template<> void iofits::write_image(int first, int last, const long * pixframe);
  template<> void iofits::write_image(int first, int last, const float * pixframe);
  template<> void iofits::write_image(int first, int last, const double * pixframe);

}

#endif
