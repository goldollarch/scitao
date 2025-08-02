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

#ifndef FITS_HEADER_DATA_H
#define FITS_HEADER_DATA_H

#include <vector>
#include <string>
#include <iomanip>
#include "iofits.h"

namespace Arroyo {

  using std::string;
  using std::vector;
  using std::ostream;

  ///
  /// A class to hold fits scale factors bscale and bzero
  ///

  class fits_scale_factor {

  protected:
    /// The fits scale factor
    double bscale;      
 
    /// The scale factor comment
    string bscale_comment;    
   
    /// The fits zero point
    double bzero;    
   
    /// The zero point comment
    string bzero_comment;        

  public:

    ///////////////////////////////////////////
    ///  Null constructor
    fits_scale_factor();

    ///////////////////////////////////////////
    ///  Construct from iofits object
    fits_scale_factor(const iofits & iof);

    ///////////////////////////////////////////
    ///  Copy constructor
    fits_scale_factor(const fits_scale_factor & fsf){
      fits_scale_factor::operator=(fsf);
    };

    ///////////////////////////////////////////
    ///  Trivial destructor
    ~fits_scale_factor(){};

    ///////////////////////////////////////////
    ///  Operator=
    fits_scale_factor & operator= (const fits_scale_factor & fsf);
 
    ///////////////////////////////////////////
    ///  Read from file
    void read(const iofits & iof){
      operator=(fits_scale_factor(iof));
    };

    ///////////////////////////////////////////
    ///  Write to file
    void write(iofits & iof) const;

    ///////////////////////////////////////////
    ///  Print
    void print(ostream & os, const char * prefix = "") const;
  };

  ///
  /// A class to hold mandatory fits header data
  ///

  template<class T>
    class fits_header_data {

    private:

    ///////////////////////////////////////////
    ///  Get the bitpix
    iofits::imagetype get_image_type() const;
    
    protected:

    /// The bitpix comment
    string bitpix_comment;       

    /// The naxis comment
    string naxis_comment;        

    /// The axes
    vector <long> axes;  

    /// The axes comments
    vector<string> axes_comment; 

    public:
    ///////////////////////////////////////////
    ///  Null constructor
    fits_header_data();

    ///////////////////////////////////////////
    ///  Construct from iofits object
    fits_header_data(const iofits & iof);

    ///////////////////////////////////////////
    ///  Copy constructor
    fits_header_data(const fits_header_data & fhd){
      fits_header_data::operator=(fhd);
    };

    ///////////////////////////////////////////
    ///  Construct from the bits
    fits_header_data(const vector<long> & in_axes, 
		     const vector<string> & in_axes_comment = vector<string>());

    ///////////////////////////////////////////
    ///  Trivial destructor
    virtual ~fits_header_data(){};

    ///////////////////////////////////////////
    ///  Operator=
    fits_header_data & operator= (const fits_header_data & fhd);
 
    ///////////////////////////////////////////
    ///  Read from iofits object
    void read(const iofits & iof);

    ///////////////////////////////////////////
    ///  Write to iofits object
    void write(iofits & iof) const;

    ///////////////////////////////////////////
    ///  Print
    virtual void print(ostream & os, const char * prefix = "") const;

    ///////////////////////////////////////////
    ///  Function to get the axes
    vector<long> get_axes() const {
      return(axes);
    };

    ///////////////////////////////////////////
    ///  Function to print the axes
    void print_axes(ostream & os) const;

    ///////////////////////////////////////////
    ///  Function to set the axes
    virtual void set_axes(const vector<long> & in_axes, 
			  const vector<long> & in_axes_comment = vector<long>(0));

    ///////////////////////////////////////////
    ///  Function to return the total number of elements
    long total_space() const;

    ///////////////////////////////////////////
    ///  Operator == for fits_header_data
    template<class U, class V>
    friend bool operator==(const fits_header_data<U> & fhd1,
				const fits_header_data<V> & fhd2);

  };

  ///////////////////////////////////////////
  ///  Operator == for fits_header_data
  template<class U, class V>
  bool operator==(const fits_header_data<U> & fhd1,
				const fits_header_data<V> & fhd2) {
    if(fhd1.get_image_type()!=fhd2.get_image_type()) return(false);
    if(fhd1.axes!=fhd2.axes) return(false);
    return(true);
  }

  ///////////////////////////////////////////
  ///  Operator !=  for fits_header_data
  template<class U, class V>
    bool operator !=(const fits_header_data<U> &fhd1,
			const fits_header_data<V> &fhd2){
    return(!operator==(fhd1,fhd2));
  }

  template<class T>
    fits_header_data<T>::fits_header_data(){
    bitpix_comment = "number of bits per data pixel";
    naxis_comment  = "number of data axes";
  }

  template<class T>
    fits_header_data<T>::fits_header_data(const iofits & iof){
    this->read(iof);
  }

  template<class T>
    fits_header_data<T>::fits_header_data(const vector<long> & in_axes,
			       const vector<string> & in_axes_comment){
    bitpix_comment = "number of bits per data pixel";
    naxis_comment  = "number of data axes";
    axes = in_axes;
    if(in_axes_comment.size()!=in_axes.size()){
      char tmp[64];
      axes_comment.resize(in_axes.size());
      for(unsigned int i=0; i<in_axes.size(); i++){
	sprintf(tmp, "length of data axis %d", i+1);
	axes_comment[i] = tmp;
      }
    } else 
      axes_comment = in_axes_comment;
  }

  template<class T>
    fits_header_data<T> & fits_header_data<T>::operator= (
    				const fits_header_data & fhd){
    if(this == &fhd) 
      return(*this);
    int i;
    bitpix_comment     = fhd.bitpix_comment;    
    naxis_comment      = fhd.naxis_comment;     
    axes               = fhd.axes;        
    axes_comment       = fhd.axes_comment;
    return(*this);
  }

  template<class T>
    void fits_header_data<T>::read(const iofits & iof) {
    int bitpix;
    iof.read_image_header(bitpix, axes);
    if(bitpix!=this->get_image_type()) {
      if(this->get_image_type()==20 && bitpix==16){
      } else {
	cerr << "fits_header_data::read error - mismatch between native bitpix "
	     << this->get_image_type() << " and iofits bitpix " << bitpix << endl;
	throw(string("fits_header_data::read"));
      }
    }
    axes_comment.resize(axes.size());
    char tmp[64];
    for(unsigned int i=0; i<axes.size(); i++){
      sprintf(tmp, "length of data axis %d", i+1);
      axes_comment[i] = tmp;
    }  
  }

  template<class T>
    void fits_header_data<T>::write(iofits & iof) const {
    iof.create_image(axes, (Arroyo::iofits::imagetype)(this->get_image_type()));
  }

  template<class T>
    void fits_header_data<T>::print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    int i;
    os.setf(std::ios::left, std::ios::adjustfield); 
    os << prefix << "BITPIX     = " << std::setw(vlspc) << this->get_image_type()
       << "/" << bitpix_comment     << std::endl;
    os << prefix << "NAXIS      = " << std::setw(vlspc) << axes.size()      
       << "/" << naxis_comment      << std::endl;
    for(i=0; i<(int)axes.size(); i++)
      os << prefix << "NAXIS" << i+1 << "     = " << std::setw(vlspc) << axes[i] 
	 << "/" << axes_comment[i] << std::endl;
  }

  template<class T>
    void fits_header_data<T>::print_axes(ostream & os) const {
    os << "fits_header_data::print_axes - axes size " << axes.size() << endl;
    for(int i=0; i<axes.size(); i++) {
      os << "fits_header_data::print_axes - axis "
         << i << " size " << axes[i] << endl;
    }
  }

  template<class T>
    long fits_header_data<T>::total_space() const {
    if(axes.size()==0) return(0);
    long nelem = 1;
    for(int i=0; i<axes.size(); i++) nelem *= axes[i];
    return(nelem);
  }

  template<class T>
    void fits_header_data<T>::set_axes(const vector<long> & in_axes, 
			       const vector<long> & in_axes_comment) {
    unsigned int naxis = in_axes.size();
    axes.resize(naxis);
    axes_comment.resize(naxis);
    for(unsigned int i=0; i<naxis; i++)
      axes[i] = in_axes[i];

    if(in_axes_comment.size()!=naxis){
      char comment[64];
      for(unsigned int i=0; i<naxis; i++){
	sprintf(comment, "length of data axis %d", i);
	axes_comment[i] = comment;
      }
    } else {
      for(unsigned int i=0; i<naxis; i++)
	axes_comment[i] = in_axes_comment[i];
    }
  }
}

#endif
