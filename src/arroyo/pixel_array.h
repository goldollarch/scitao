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

#ifndef PIXEL_ARRAY_H
#define PIXEL_ARRAY_H

#include <vector>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include <cfloat>
#include <iostream>
#include <fstream>
#include "AO_algo.h"
#include "colormap.h"
#include "iofits.h"
#include "fits_header_data.h"
#include "fft_manager.h"

using namespace std;

namespace Arroyo {

  using std::string;
  using std::vector;
  using std::cerr;
  using std::endl;
  using std::ostream;
    
  /* forward declarations */
  template <class T> class pixel_array;
  template <class T> pixel_array<T> operator + (const pixel_array<T> &p1,
  						const pixel_array<T> &p2);
  template <class T> pixel_array<T> operator - (const pixel_array<T> &p1,
  						const pixel_array<T> &p2);
  template <class T> pixel_array<T> operator * (const pixel_array<T> &p1,
  						const pixel_array<T> &p2);
  template <class T> pixel_array<T> operator / (const pixel_array<T> &p1,
  						const pixel_array<T> &p2);
  template <class T> pixel_array<T> operator + (const pixel_array<T> &p1,
  						double & fac);
  template <class T> pixel_array<T> operator - (const pixel_array<T> &p1,
  						double & fac);
  template <class T> pixel_array<T> operator * (const pixel_array<T> &p1,
  						double & fac);
  template <class T> pixel_array<T> operator / (const pixel_array<T> &p1,
  						double & fac);
  template <class T> bool operator == (const pixel_array<T> &p1,
				       const pixel_array<T> &p2);
  template <class T> bool operator != (const pixel_array<T> &p1,
				       const pixel_array<T> &p2);

  ///
  /// A class to hold and manipulate rectangular images.
  /// Template parameter specifies storage type of 
  /// data - e.g. double, float, int. 

  template<class T> 
    class pixel_array {

    private:

    ///////////////////////////////////////////
    ///  Implementation of the rotate_and_shift_by_fft 
    ///  function.  Rotation is performed through three
    ///  shearing operations.  These shearing operations 
    ///  are performed using the shift theorem.  The first 
    ///  and last require a 1D forward and backward fft of 
    ///  each row in the image.  The middle requires a 1D 
    ///  forward and backward fft of each column in the image.
    void simple_rotate_and_shift_by_fft(double dx, double dy,
    					double angle, bool window);

    mutable int private_nelem;

    protected:

    /// axes of underlying array
    vector<long> axes;
  
    /// pointer to the data
    T * pixeldata;

    /// pointer to the weights
    float * pixelwts;

    ///////////////////////////////////////////
    ///  Reallocate memory for this pixel_array
    ///  pixeldata is allocated with in_axes elements
    ///  and initialized to zero
    ///  if pixelwts is null it remains so, otherwise
    ///  it is allocated and initialized to zero
    void set_axes(const vector<long> & in_axes);

    ///////////////////////////////////////////
    ///  Allocate weights.  If weights have already
    ///  been allocated, they are reinitialized to wt
    ///  If pixeldata == NULL, pixelwts is set to NULL
    void allocate_weights(double wt);

    ///////////////////////////////////////////
    ///  Deallocate weights 
    void deallocate_weights();

    ///////////////////////////////////////////
    ///  Function to normalize data by weights
    double normalize_by_wts();

    ///////////////////////////////////////////
    ///  Function to bin together n amplitudes
    ///  So far this has only been coded for 2d 
    ///  pixel_arrays
    void decimate(int nadd);

    ///////////////////////////////////////////
    ///  Function to flag weights that are zero in arg 
    template<class U>
      void flag_zero_wts(const pixel_array<U> & pixarr);
    
    ///////////////////////////////////////////
    ///  Function to make this pixel_array
    ///  into a mask of 1s and 0s based
    ///  on the value of the pixeldata
    void mask();

    ///////////////////////////////////////////
    ///  Function to apply a mask
    template<class U>
      void mask(const pixel_array<U> * pixarr);

    ///////////////////////////////////////////
    ///  Construct from iofits object
    pixel_array(const iofits & iof);

    public:

    ///////////////////////////////////////////
    ///  Null constructor  
    pixel_array();

    ///////////////////////////////////////////
    ///  Copy constructor
    ///  template<class U>  - there is not yet an
    ///  operator = for different template types
    pixel_array(const pixel_array<T> & pixarr);

    ///////////////////////////////////////////
    ///  Construct an instance with pixel limits given by pixel_limits
    ///  These limits must be contained by pixarr's arrays
    ///
    ///  The array pixel_limits has 4 elements:
    ///  xmin, xmax, ymin, ymax
    ///
    ///  So far this has only been coded for 2d 
    ///  pixel_arrays
    template<class U>
      pixel_array(const pixel_array<U> & pixarr,
		  const vector<long> & pixel_limits);

    ///////////////////////////////////////////
    ///  Construct from arrays 
    ///  If wts == NULL, weights will not be allocated
    ///  If data == NULL && wts == NULL data will be allocated
    ///    and initialized to zero
    pixel_array(const vector<long> & in_axes, 
		const T * data = NULL, 
		const float * wts = NULL);

    ///////////////////////////////////////////
    ///  Destructor
    virtual ~pixel_array();  

    ///////////////////////////////////////////
    ///  Operator =
    pixel_array<T> & operator = (const pixel_array<T> & pixarr);

    ///////////////////////////////////////////
    ///  Function to read data from iofits
    virtual void read(const iofits & iof);

    ///////////////////////////////////////////
    ///  Function to write data to iofits
    virtual void write(iofits & iof) const;

    ///////////////////////////////////////////
    ///  Function to write data to a ppm file
    void write_to_ppm(double min, double max,
		      bool logscale,
		      bool colorbar,
		      colormap * cmap, 
		      const char * filename,
		      long min_dimen = -1) const;

    ///////////////////////////////////////////
    ///  Function to compute total number of elements in array
    inline int total_space() const {

      if(this->private_nelem>0){
	return(this->private_nelem);
      }
      //cout << "\tATS " << this->private_nelem << "\t";
      int naxes = axes.size();
      if(naxes == 0) 
	this->private_nelem=0;
      else {
	this->private_nelem = 1;
	for(int i=0; i<naxes; i++)
	  this->private_nelem *= axes[i];
      }
      //cout << "\tBTS " << this->private_nelem << endl;
      return(this->private_nelem);
    };

    ///////////////////////////////////////////
    ///  Function to report whether weights are allocated
    bool weights_allocated() const;

    ///////////////////////////////////////////
    ///  Function to set the axes  
    vector<long> get_axes() const {return(axes);};

    ///////////////////////////////////////////
    ///  Function to copy pixel array 
    ///  This function performs template conversions
    ///  through implicit casts.
    template<class U>
      void copyfrom(const pixel_array<U> & pixarr);

    ///////////////////////////////////////////
    ///  Function to print the axes
    void print_axes(ostream & os) const;

    ///////////////////////////////////////////
    ///  Get nth data element
    /// 
    ///  The data is indexed as n = i*axes[0]+j, where
    ///
    ///  0 <= i < axes[1]
    ///
    ///  0 <= j < axes[0]
    T data(long n) const {
      if(n<0 || n>=this->total_space()){
	cerr << "pixel_array::data error - element " << n
	     << " out of range\n";
	throw(string("pixel_array::data"));
      }
      return(pixeldata[n]);
    };

    ///////////////////////////////////////////
    ///  Set nth data element to val
    /// 
    ///  The data is indexed as n = i*axes[0]+j, where
    ///
    ///  0 <= i < axes[1]
    ///
    ///  0 <= j < axes[0]
    void set_data(long n, T val) const {
      if(n<0 || n>=this->total_space()){
	cerr << "pixel_array::set_data error - element " << n
	     << " out of range\n";
	throw(string("pixel_array::set_data"));
      }
      pixeldata[n] = val;
    };

    ///////////////////////////////////////////
    ///  Function to return weight element
    float wt(long elem) const{
      if(pixelwts!=NULL) return(*(pixelwts+elem));
      else {
	cerr << "pixel_array::wt error - weights not allocated\n";
	throw(string("pixel_array::wt"));
      }
    };
  
    ///////////////////////////////////////////
    ///  Find min and max of pixel array
    void min_and_max(double & min, double & max) const;

    ///////////////////////////////////////////
    ///  Find min and max of pixel array
    void min_and_max(double & min, vector<int> & minpixel, 
		     double & max, vector<int> & maxpixel) const;

    ///////////////////////////////////////////
    ///  Find min and max within limits
    void min_and_max(double & min, vector<int> & minpixel, 
		     double & max, vector<int> & maxpixel,
		     vector<int> axis_0_limits, 
		     vector<int> axis_1_limits) const;

    ///////////////////////////////////////////
    ///  Flip data about x axis
    void flip_x();
  
    ///////////////////////////////////////////
    ///  Flip data about y axis
    void flip_y();
  
    ///////////////////////////////////////////
    ///  Flip data about x and y axis
    void flip_xy();

    ///////////////////////////////////////////
    ///  Flip data about 45 degree line
    void flip_45();

    ///////////////////////////////////////////
    ///  Pad each edge of the array by npad pixels
    ///  and initialize to the specified value
    void pad_array(int npad, double value = 0);

    ///////////////////////////////////////////
    /// Clip each edge of the array by nclip pixels 
    void clip_array(int nclip);

    ///////////////////////////////////////////
    ///  Shift a pixel_array by dx and dy
    void shift_by_fft(double dx, double dy);

    ///////////////////////////////////////////
    ///  Rotate pixel_array by an angle 
    /// 
    ///  The resulting pixel array may be larger than the original, so as
    ///  to be able to include the corners that have rotated outside the
    ///  original array bounds.
    void rotate_by_fft(double angle, bool window=true);

    ///////////////////////////////////////////
    ///  Rotate pixel_array by an angle 
    ///  and shift by dx and dy
    /// 
    ///  The resulting pixel array may be larger than the original, so as
    ///  to be able to include the corners that have rotated outside the
    ///  original array bounds.
    ///
    ///  This function does not yet support non-zero shifts
    void rotate_and_shift_by_fft(double dx, 
				 double dy,
				 double angle, 
				 bool window=true);

    ///////////////////////////////////////////
    ///  Implementation of the rotate_and_shift_by_fft 
    ///  function, based on an incomprehensible technical
    ///  report by owen and makedon.
    void owen_makedon_rotate_and_shift_by_fft(double dx, 
					      double dy, 
					      double angle);

    ///////////////////////////////////////////
    ///  Function to cross-correlate with another pixel_array
    ///  Cross correlation performed via fft
    pixel_array<T> cross_correlate(const pixel_array<T> & pixarr) const;

    ///////////////////////////////////////////
    ///  Function to fit the offset 
    ///  between the two pixel arrays.  The computation
    ///  is performed using a shift by fft
    ///  This function returns the number of pixels you have to shift 
    ///  "this" by in order to optimally align it with "pixarr"
    void offset(const pixel_array<T> & pixarr,
		vector<double> & offsets, 
		long range = -1) const;

    ///////////////////////////////////////////
    ///  Friend declaration for operator +=  for pixel_arrays
    template<class U, class V> 
      friend pixel_array<U> & operator+=(pixel_array<U> & lhs,
					 const pixel_array<V> & rhs);

    ///////////////////////////////////////////
    ///  Friend declaration for operator -=  for pixel_arrays
    template<class U, class V> 
      friend pixel_array<U> & operator-=(pixel_array<U> & lhs,
					 const pixel_array<V> & rhs);

    ///////////////////////////////////////////
    ///  Friend declaration for operator *=  for pixel_arrays
    template<class U, class V>  
      friend pixel_array<U> & operator*=(pixel_array<U> & lhs,
					 const pixel_array<V> & rhs);

    ///////////////////////////////////////////
    ///  Friend declaration for operator /=  for pixel_arrays
    template<class U, class V>
      friend pixel_array<U> & operator/=(pixel_array<U> & lhs,
					 const pixel_array<V> & rhs);

    ///////////////////////////////////////////
    ///  Operator +=  for doubles
    virtual pixel_array<T> & operator += (const double & fac);

    ///////////////////////////////////////////
    ///  Operator -=  for doubles 
    virtual pixel_array<T> & operator -= (const double & fac);

    ///////////////////////////////////////////
    ///  Operator *=  for doubles 
    virtual pixel_array<T> & operator *= (const double & fac);

    ///////////////////////////////////////////
    ///  Operator /=  for doubles 
    virtual pixel_array<T> & operator /= (const double & fac);

    ///////////////////////////////////////////
    ///  Friend operator ==  for pixel_array
    friend bool operator ==(const pixel_array<T> &p1, const pixel_array<T> &p2){

      if(p1.axes.size()!=p2.axes.size()) return(0);

      if((p1.pixelwts!=NULL && p2.pixelwts==NULL) ||
	 (p1.pixelwts==NULL && p2.pixelwts!=NULL)){
	return(0);
      }

      int nelem = 1;
      for(uint i=0; i<p1.axes.size(); i++){
	if(p1.axes[i]!=p2.axes[i]) 
	  return(0);
	nelem *= p1.axes[i];
      }

      for(int i=0; i<nelem; i++)
	if(p1.pixeldata[i]!=p2.pixeldata[i]) return(0);

      if(p1.pixelwts!=NULL)
	for(int i=0; i<nelem; i++)
	  if(p1.pixelwts[i]!=p2.pixelwts[i]) return(0);

      return(1);
    };

    ///  Verbose level
    static int verbose_level;

  };

  template<class T>
    int pixel_array<T>::verbose_level = 0;

  template<class T>
    void pixel_array<T>::set_axes(const vector<long> & in_axes){
    
    int nelem = in_axes.size()==0 ? 0 : 1;
    for(unsigned int i=0; i<in_axes.size(); i++)
      nelem *= in_axes[i];

    if(nelem<0){
      cerr << "pixel_array::set_axes error - "
	   << "total number of elements " << nelem
	   << " less than zero\n";
      throw(string("pixel_array::set_axes"));
    }

    if(pixel_array<T>::verbose_level) 
      cout << "pixel_array::set_axes - have " << total_space() 
	   << " and need " << nelem << endl;

    // Check to see if we can avoid memory reallocation
    if(total_space() != nelem){
      
      axes = in_axes;

      // reset this flag so that this will be recomputed after
      // axes have been set to in_axes
      this->private_nelem = -1;

      if(pixel_array<T>::verbose_level) 
	cout << "pixel_array::set_axes - reallocating memory\n";
      
      bool wts_initialized = this->weights_allocated();

      if(pixeldata!=NULL) 
	delete [] pixeldata;
      if(pixelwts!=NULL){
	delete [] pixelwts;
	pixelwts = NULL;
      }

      if(nelem==0){
	pixeldata=NULL;
	pixelwts=NULL;
	axes = in_axes;

	return;
      }

      try{
	if(pixel_array<T>::verbose_level) 
	  cout << "pixel_array::set_axes - allocating data: nelems " 
	       << nelem << endl;
	pixeldata = new T[nelem];
      } catch(...){
	cerr << "pixel_array::set_axes - error reallocating space for data array\n";
	throw(string("pixel_array::set_axes"));
      }
      
      for(int i=0; i<nelem; i++)
	pixeldata[i] = 0;

      if(wts_initialized) 
	this->allocate_weights(0);
    } else {
      // This is the case where the number of elements
      // remain the same, but their dimensional distribution
      // may change.  This includes the important case when
      // there are no elements, but the number of dimensions
      // is specified through in_axes.size()
      axes = in_axes;
      for(int i=0; i<nelem; i++)
	pixeldata[i] = 0;
      if(this->weights_allocated())
	for(int i=0; i<nelem; i++)
	  pixelwts[i] = 0;
    }
  }

  template<class T>
    void pixel_array<T>::allocate_weights(double wt){

    int nelem = total_space();
    if(pixeldata==NULL || nelem == 0){
      pixelwts=NULL;
      return;
    }

    if(pixelwts!=NULL){
      for(int i=0; i<nelem; i++)
	pixelwts[i] = wt;
      return;
    }

    try{
      if(pixel_array<T>::verbose_level) 
	cout << "pixel_array::allocate_weights - " << total_space() 
	     << " weights being initialized to " << wt << endl;
      pixelwts = new float[nelem];}
    catch(...){
      cerr << "pixel_array::allocate_weights - could not allocate memory\n";
      throw(string("pixel_array::allocate_weights"));
    }
    for(int i=0; i<nelem; i++) pixelwts[i] = wt;
  }

  template<class T>
    void pixel_array<T>::deallocate_weights(){
    if(pixelwts==NULL) return;
    delete [] pixelwts;
    pixelwts = NULL;
  }

  template<> double pixel_array<long>::normalize_by_wts();

  template<class T>
    double pixel_array<T>::normalize_by_wts(){
    double maxwt = 0;
    if(!weights_allocated()) 
      return(0);
    int nelem = total_space();
    for(int i=0; i<nelem; i++){
      if(pixelwts[i]>maxwt && pixelwts[i]!=0)
	maxwt = pixelwts[i];
    }
    for(int i=0; i<nelem; i++){
      if(pixelwts[i]!=maxwt && pixelwts[i]!=0){
	pixeldata[i] *= maxwt/(float)pixelwts[i];
	pixelwts[i] = maxwt;
      } else if(pixelwts[i]==0) 
	pixeldata[i] = 0;
    }
    return(maxwt);
  }

  template<class T>
    void pixel_array<T>::decimate(int nadd){
    if(axes.size()!=2){
      cerr << "pixel_array::decimate - cannot decimate "
           << "pixel_array with number of dimensions " 
	   << axes.size() << endl;
      throw(string("pixel_array::decimate"));
    }

    if(nadd<0 || nadd>axes[0] || nadd>axes[1]){
      cerr << "pixel_array::decimate - error decimating by a factor of "
           << nadd << endl;
      throw(string("pixel_array::decimate"));
    }

    if(nadd==0 || nadd==1) return;

    vector<long> newaxes(2);
    for(uint i=0; i<newaxes.size(); i++) 
      newaxes[i] = axes[i]/nadd;

    T * olddata = pixeldata;
    float * oldwts;
    
    pixeldata = new T[newaxes[0]*newaxes[1]];

    if(weights_allocated()){
      oldwts = pixelwts;
      pixelwts = new float[newaxes[0]*newaxes[1]];
    }

    int wtsum;
    for(int i=0; i<newaxes[1]; i++){
      for(int j=0; j<newaxes[0]; j++){
	pixeldata[i*newaxes[0]+j] = 0;
	if(weights_allocated()){
	  pixelwts[i*newaxes[0]+j] = 0;
	  wtsum = 0;
	}
	for(int k=0; k<nadd; k++){
	  for(int l=0; l<nadd; l++){
	    if(weights_allocated()){
	      pixeldata[i*newaxes[0]+j] += 
		olddata[(i*nadd+k)*axes[0]+j*nadd+l]*
		oldwts[(i*nadd+k)*axes[0]+j*nadd+l];
	      pixelwts[i*newaxes[0]+j] += oldwts[(i*nadd+k)*axes[0]+j*nadd+l];
	    } else {
	      pixeldata[i*newaxes[0]+j] += olddata[(i*nadd+k)*axes[0]+j*nadd+l];
	    }
	  }
	}
      }
    }

    axes = newaxes;
    delete [] olddata;
  }


  template<class T>
    void pixel_array<T>::pad_array(int npad, double value) {
    if(npad<0){
      cerr << "pixel_array::pad_array error - cannot pad by "
           << npad << " pixels\n";
      throw(string("pixel_array::pad_array"));
    }
    if(npad==0 || axes.size()==0) return;
    
    if(axes.size()!=2){
      cerr << "pixel_array::pad_array error - cannot pad by "
           << axes.size() << " dimensions, "
	   << " as this generalization has not yet been coded\n";
      throw(string("pixel_array::pad_array"));
    }    

    vector<long> new_dimen = axes;
    int nelem = 1;
    for(uint i=0; i<new_dimen.size(); i++){
      new_dimen[i] += 2*npad;
      nelem *= new_dimen[i];
    }

    // allocate new arrays
    T * newdata = new T[nelem];
    float * newwts = NULL;
    if(weights_allocated())
      newwts = new float[nelem];

    // initialize arrays to the value
    for(int i=0; i<new_dimen[1]; i++){
      for(int j=0; j<new_dimen[0]; j++){
	newdata[i*new_dimen[0]+j] = value;
	if(weights_allocated())
	  newwts[i*new_dimen[0]+j] = 0;
      }
    }

    // copy the old data into the new array
    for(int i=0; i<axes[1]; i++){
      for(int j=0; j<axes[0]; j++){
	newdata[(i+npad)*new_dimen[0]+npad+j] = 
	  pixeldata[i*axes[0]+j];
	if(weights_allocated())
	  newwts[(i+npad)*new_dimen[0]+npad+j] = 
	    pixelwts[i*axes[0]+j];
      }
    }
  
    delete [] pixeldata;
    pixeldata = newdata;

    if(weights_allocated()){
      delete pixelwts;
      pixelwts = newwts;
    }

    axes = new_dimen;
  }

  template<class T>
    void pixel_array<T>::clip_array(int nclip) {

    if(nclip==0) return;
    if(nclip<0){
      cerr << "pixel_array::clip_array error - cannot clip by "
           << nclip << " pixels\n";
      throw(string("pixel_array::clip_array"));
    }

    vector<long> new_axes = axes;
    int nelem = 1;
    for(uint i=0; i<new_axes.size(); i++){
      new_axes[i] -= 2*nclip;
      if(new_axes[i]<=0){
	cerr << "pixel_array::clip_array error - clipping original array by "
	     << nclip << " pixels leaves non-positive array size\n";
	throw(string("pixel_array::clip_array"));
      }
      nelem *= new_axes[i];
    }

    // allocate new arrays
    T * newdata = new T[nelem];
    float * newwts = NULL;
    if(weights_allocated())
      newwts = new float[nelem];

    // copy the old data into the new array
    for(int i=0; i<new_axes[1]; i++)
      for(int j=0; j<new_axes[0]; j++){
	newdata[i*new_axes[0]+j] = 
	  pixeldata[(i+nclip)*axes[0]+nclip+j];
	if(weights_allocated())
	  newwts[i*new_axes[0]+j+1] = 
	    pixelwts[(i+nclip)*axes[0]+nclip+j];
      }

    // make the switcheroo
    delete [] pixeldata;
    pixeldata = newdata;
    if(weights_allocated()){
      delete [] pixelwts;
      pixelwts = newwts;
    }
    axes = new_axes;
  } 

  template<class T>
    template<class U>
    void pixel_array<T>::flag_zero_wts(const pixel_array<U> & pixarr){

    if(!pixarr.weights_allocated()){
      cerr << "pixel_array::flag_zero_wts error - "
           << "reference array has no weights\n";
      throw(string("pixel_array::flag_zero_wts"));
    }
    if(!this->weights_allocated()) 
      this->allocate_weights(1);
    int nelem = total_space();
    for(int i=0; i<nelem; i++){
      if(pixarr.wt(i)==0){
	pixeldata[i] = 0;
	pixelwts[i] = 0;
      }
    }
  }
    
  template<class T>
    void pixel_array<T>::mask(){
    if(!this->weights_allocated())
      this->allocate_weights(1);
    int nelem = total_space();
    for(int i=0; i<nelem; i++){
      if(pixeldata[i]==0){
	pixelwts[i] = 0;
      } else {
	pixeldata[i] = 1;
	pixelwts[i] = 1;
      }
    }
  }

  template<class T>
    template<class U>
    void pixel_array<T>::mask(const pixel_array<U> * pixarr){
    if(axes!=pixarr->get_axes()){
      cerr << "pixel_array::mask error - mismatched axes\n";
      throw(string("pixel_array::mask"));
    }
    if(this->weights_allocated()==0)
      this->allocate_weights(1);
    int count = 0;
    int nelem = total_space();
    for(int i=0; i<nelem; i++){
      if(pixarr->data(i)==0){
	pixeldata[i] = 0;
	pixelwts[i] = 0;
	count++;
      } 
    }
  }

  template<class T>
    void pixel_array<T>::read(const iofits & iof){
  
    //int ndimen = iof.get_img_dim();
    vector<long> tmpaxes = iof.get_img_size();
    this->set_axes(tmpaxes);

    int nelems = this->total_space();

    try{
      if(pixel_array<T>::verbose_level) 
	cout << "pixel_array::pixel_array - reading data\n";
      iof.read_image(0, nelems-1, pixeldata);
    } catch(...){
      cerr << "pixel_array::pixel_array - error reading image\n";
      throw(string("pixel_array::pixel_array"));
    }

    // INCOMPATIBILITY!!!!
    // this is a bad idea, as it relies on a global call,
    // so that pixel array cannot be written as part of
    // an aggregate class
    // Possible solution - write a weights keyword to the
    // header, but at the same time add something to the
    // AO_observation loader that can compensate for the
    // fact that this keyword doesn't exist in the present
    // header format.
    /*
      if(iof.get_num_hdus()==2){
      this->allocate_weights(0);
      if(pixel_array<T>::verbose_level) 
      cout << "pixel_array::pixel_array - reading weights\n";
      iof.movabs_hdu(2);
      iof.read_image(0,nelems-1,pixelwts);
      } else pixelwts = NULL;
    */

    // Backwards compatibility hack:

    // Before types other than AO_observations were
    // supported, I implicitly assumed that all fits 
    // files contained a primary array and zero or one
    // image extensions, depending on whether the 
    // weights were defined for a particular 
    // AO_observation.  One could check the total number 
    // of hdus to determine which case this was.

    // With the addition of other classes that required
    // fits file storage, this simple assumption became
    // too restrictive.  In particular, I wanted to be
    // able to write many pixel_arrays to the same fits
    // file, so that I could store aggregate objects.
    // Thus additional information was required to 
    // resolve whether the next image extension contained
    // weights from this pixel array or a new pixel array 
    // altogether - or some other form of data.

    // In order to fix this while maintaining backwards
    // compatibility, I added a weights keyword to the 
    // headers of the new classes.  If the weights keyword
    // is not present then this pixel array is assumed
    // to belong to an AO_observation, and the number
    // of hdu's is checked as before.  If the weights
    // key is present and is true, then we advance to
    // the next image extension and read the weights.
    // If the weights key is false, we stop.

    // Here we check whether there is a TYPE key,
    // which indicates that this pixel array is
    // part of an object other than an AO_observation.

    if(!iof.key_exists("WEIGHTS")){
      if(iof.get_num_hdus()==2){
	if(pixel_array<T>::verbose_level) 
	  cout << "pixel_array::pixel_array - reading weights\n";
	this->allocate_weights(0);
	iof.movabs_hdu(2);
	iof.read_image(0,nelems-1,pixelwts);
      } else pixelwts = NULL;
    } else {
      bool weights;
      string comment;
      iof.read_key("WEIGHTS", weights, comment);
      if(weights){
	if(pixel_array<T>::verbose_level) 
	  cout << "pixel_array::pixel_array - reading weights\n";
	this->allocate_weights(0);
	iof.movrel_hdu(1);
	iof.read_image(0, nelems-1, pixelwts);
      } 
      
      if(iof.get_hdu_num()<iof.get_num_hdus())
	iof.movrel_hdu(1);
    }
  }

  template<class T>
    void pixel_array<T>::write(iofits & iof) const {

    int space = total_space();
    if(space!=0) iof.write_image(0,space-1,pixeldata);

    string comment = "weights present";
    if(pixelwts!=NULL){
      iof.write_key("WEIGHTS", true, comment);
      Arroyo::fits_header_data<float> fhd(this->get_axes());
      fhd.write(iof);
      if(space!=0) iof.write_image(0,space-1,pixelwts); 
    } else       
      iof.write_key("WEIGHTS", false, comment);
 
  }

  template<class T>
    void pixel_array<T>::write_to_ppm(double min, double max,
				      bool logscale,
				      bool colorbar,
				      colormap * cmap, 
				      const char * filename, 
				      long min_dimen) const {
  
    if(axes.size()!=2){
      cerr << "pixel_array::write_to_ppm error - array has "
	   << axes.size() << " axes, rather than two\n";
      throw(string("pixel_array::write_to_ppm"));
    }

    if(min>=max){
      cerr << "pixel_array::write_to_ppm error - min " 
	   << min << " and max " << max 
	   << " supplied to this function are inconsistent\n";
      throw(string("pixel_array::write_to_ppm"));
    }

    long ppm_dimen;
    long ppm_pix_per_elem;
    if(min_dimen==-1){
      ppm_dimen = axes[0]>axes[1]?axes[0] : axes[1];
      ppm_pix_per_elem = 1;
    } else {
      if(min_dimen<axes[0] || min_dimen<axes[1]){
	cerr << "pixel_array::write_to_ppm error - minimum dimension " 
	     << min_dimen << " less than array dimensions "
	     << axes[0] << "x" << axes[1] << endl;
	throw(string("pixel_array::write_to_ppm"));
      }      
      ppm_dimen = min_dimen;
      ppm_pix_per_elem = axes[0]>axes[1] ? min_dimen/axes[0] : min_dimen/axes[1];
    }
    double ppm_min = min;

    long axis_zero_min =
      axes[0]>axes[1] ? 0 : ((axes[1]-axes[0])*ppm_pix_per_elem)/2;
    long axis_zero_max =
      axes[0]>axes[1] ? ppm_dimen : axes[0]*ppm_pix_per_elem+axis_zero_min;
    long axis_one_min =
      axes[1]>axes[0] ? 0 : ((axes[0]-axes[1])*ppm_pix_per_elem)/2;
    long axis_one_max =
      axes[1]>axes[0] ? ppm_dimen : axes[1]*ppm_pix_per_elem+axis_one_min;

    /*    
	  cout << min_dimen 
	  << "\t" << ppm_dimen 
	  << "\t" << ppm_pix_per_elem
	  << "\t" << axes[0] << "\t" << axes[1]
	  << "\t" << axis_zero_min << "\t" << axis_zero_max
	  << "\t" << axis_one_min << "\t" << axis_one_max
	  << endl;
    */

    std::ofstream fs;
    fs.open(filename);
    fs << "P6\n";

    long colorbar_padding = 0;
    if(colorbar)
      colorbar_padding = .1*ppm_dimen>60 ? (long)(.1*ppm_dimen) : 60;
    fs << ppm_dimen << " " << ppm_dimen+colorbar_padding << endl;
    fs << "255\n";
    fs.close();

    fs.open(filename,
	    std::ios::app | std::ios::out | std::ios::binary | std::ios::ate);

    // Here I inverted the write order on axes[1]
    // to force display to show the results in the
    // same orientation as does ds9
    int index;
    int nelem = axes[0]*axes[1];
    for(int i=ppm_dimen-1; i>=0; i--){
      for(int j=0; j<ppm_dimen; j++){
	if(i<axis_one_min || i>=axis_one_max || j<axis_zero_min || j>=axis_zero_max)
	  fs << cmap->get_R(ppm_min, min, max, logscale)
	     << cmap->get_G(ppm_min, min, max, logscale)
	     << cmap->get_B(ppm_min, min, max, logscale);
	else {
	  index = (((i-axis_one_min)/ppm_pix_per_elem)*axes[0])
	    +((j-axis_zero_min)/ppm_pix_per_elem);
	  fs << cmap->get_R(pixeldata[index], min, max, logscale)
	     << cmap->get_G(pixeldata[index], min, max, logscale)
	     << cmap->get_B(pixeldata[index], min, max, logscale);
	}
      }
    }

    if(colorbar){
      int first_limit = (long)(.3*colorbar_padding);
      int second_limit = (long)(.5*colorbar_padding);

      for(int i=0; i<first_limit; i++)
	for(int j=0; j<ppm_dimen; j++)	
	  fs << cmap->get_R(ppm_min, min, max)
	     << cmap->get_G(ppm_min, min, max)
	     << cmap->get_B(ppm_min, min, max);

      for(int i=first_limit; i<second_limit; i++){
	int first = (long)(.1*ppm_dimen);
	int second = (long)(.9*ppm_dimen);
	for(int j=0; j<first-1; j++)
	  fs << cmap->get_R(ppm_min, min, max)
	     << cmap->get_G(ppm_min, min, max)
	     << cmap->get_B(ppm_min, min, max);

	fs << '0' << '0' << '0';

	for(int j=first; j<second; j++){
	  if(i==first_limit || i==second_limit-1)
	    fs << '0' << '0' << '0';
	  else 
	    fs << cmap->get_R(min+(j-first)*(max-min)/(double)(second-first), min, max)
	       << cmap->get_G(min+(j-first)*(max-min)/(double)(second-first), min, max)
	       << cmap->get_B(min+(j-first)*(max-min)/(double)(second-first), min, max);
	}

	fs << '0' << '0' << '0';

	for(int j=second+1; j<ppm_dimen; j++)
	  fs << cmap->get_R(ppm_min, min, max)
	     << cmap->get_G(ppm_min, min, max)
	     << cmap->get_B(ppm_min, min, max);
      }

      for(int i=second_limit; i<colorbar_padding; i++)
	for(int j=0; j<ppm_dimen; j++)	
	  fs << cmap->get_R(ppm_min, min, max)
	     << cmap->get_G(ppm_min, min, max)
	     << cmap->get_B(ppm_min, min, max);
    }
  }

  template<class T>
    pixel_array<T>::pixel_array() {
    private_nelem = -1;
    pixeldata = NULL; 
    pixelwts = NULL;
  }

  template<class T>
    pixel_array<T>::pixel_array(const pixel_array<T> & pixarr){
    private_nelem = -1;
    pixeldata = NULL;
    pixelwts = NULL;
    pixel_array::operator=(pixarr);
  }

  template<class T>
    pixel_array<T>::pixel_array(const iofits & iof){
    private_nelem = -1;
    pixeldata = NULL;
    pixelwts = NULL;
    this->read(iof);
  }

  template<class T>
    template<class U>
    pixel_array<T>::pixel_array(const pixel_array<U> & pixarr,
    				const vector<long> & pixel_limits){


    private_nelem = -1;
    vector<long> paxes = pixarr.get_axes();
    if(paxes.size()!=2){
      cerr << "pixel_array::pixel_array error - "
           << "cannot construct through extraction from " 
	   << axes.size() << " dimensions, "
	   << " as this generalization has not yet been coded\n";
      throw(string("pixel_array::pixel_array"));
    }    

    if(pixel_limits.size()!=2*paxes.size()){
      cerr << "pixel_array::pixel_array error - "
	   << "cannot construct from a pixel array of other than 2 axes\n";
      throw(string("pixel_array::pixel_array"));
    }      
    
    if(pixel_array<T>::verbose_level)
      cout << "pixel_array::pixel_array - original axes " 
	   << paxes[0] << "\t" << paxes[1] << endl;

    vector<long> tmpaxes(paxes.size());

    for(int i=0; i<pixel_limits.size(); i+=2){

      if(pixel_limits[i+1] >= paxes[i/2] || pixel_limits[i] < 0){
	cerr << "pixel_array::pixel_array - "
	     << "requested pixels out of ";
	if(i==0) cerr << "horizontal";
	else cerr << "vertical";
	cerr << " range of input pixel_array:\n\t" 
	     << "requested range\t" << pixel_limits[i]
	     << " - " << pixel_limits[i+1]
	     << "\n\tinput range\t" << 0 << " - " << paxes[i/2]-1 << endl;
	throw(string("pixel_array::pixel_array"));
      }	

      tmpaxes[i/2] = pixel_limits[i+1] - pixel_limits[i] + 1;
      if(tmpaxes[i/2]<=0){
	cerr << "pixel_array::pixel_array - undefined pixel limits " 
	     << pixel_limits[i+1] << " - " << pixel_limits[i] << endl;
	throw(string("pixel_array::pixel_array"));
      }
    }

    pixeldata = NULL;
    pixelwts = NULL;
    this->set_axes(tmpaxes);
    if(pixarr.weights_allocated())
      this->allocate_weights(0);

    // Note:  fits files are stored in fortran array format rather than C format
    // Hence the loops are reversed 
    if(pixel_array<T>::verbose_level)
      cout << "pixel_array::pixel_array - extracting horizontal " 
	   << pixel_limits[0] << " - " << pixel_limits[1]
	   << " vertical " << pixel_limits[2] << " - " << pixel_limits[3] << endl;
    int tmpa, tmpb;
    int axes_one = axes[1];
    int axes_zero = axes[0];
    int paxes_zero = paxes[0];
    long pixel_limits_zero = pixel_limits[0];
    long pixel_limits_two = pixel_limits[2];
    for(int i=0; i<axes_one; i++){
      for(int j=0; j<axes_zero; j++){
	tmpa = i*axes_zero+j;
	tmpb = (pixel_limits_two+i)*paxes_zero+pixel_limits_zero+j;
	pixeldata[tmpa] = pixarr.data(tmpb);
	if(pixelwts!=NULL) pixelwts[tmpa] = pixarr.wt(tmpb);
      }
    }
  }

  template<class T>
    pixel_array<T>::pixel_array(const vector<long> & in_axes, 
				const T * data, 
				const float * wts){

    if(data == NULL && wts != NULL){
      cerr << "pixel_array::pixel_array error - data == NULL, wts != NULL" << endl;
      throw(string("pixel_array::pixel_array"));
    }

    private_nelem = -1;
    pixeldata = NULL;
    pixelwts = NULL;
    this->set_axes(in_axes);

    int nelems = total_space();
    
    if(data!=NULL){
      for(int i=0; i<nelems; i++) 
	pixeldata[i] = data[i];

      if(wts!=NULL){
	this->allocate_weights(0);
	for(int i=0; i<nelems; i++) 
	  pixelwts[i] = wts[i];
      }
    } else 
      for(int i=0; i<nelems; i++) 
	pixeldata[i] = 0;
  }
  
  template<class T>
    pixel_array<T>::~pixel_array(){
    if(pixeldata!=NULL)
      delete [] pixeldata;
    if(pixelwts!=NULL)
      delete [] pixelwts;
  }

  template<class T>
    pixel_array<T> & pixel_array<T>::operator = (const pixel_array<T> & pixarr){

    if(this == &pixarr){
      return(*this);
    }
    this->copyfrom(pixarr);
    return(*this);
  }

  template<class T>
    bool pixel_array<T>::weights_allocated() const {
    if(pixelwts == NULL) return(0);
    return(1);
  }

  template<class T>
    template<class U>
    void pixel_array<T>::copyfrom(const pixel_array<U> & pixarr){

    if((void*)this == (void*)&pixarr) return;

    this->set_axes(pixarr.get_axes());

    // Copy over the data
    int nelem = pixarr.total_space();
    for(int i=0; i<nelem; i++)
      pixeldata[i] = static_cast<T>(pixarr.data(i));

    // Copy over the weights if they exist
    if(pixarr.weights_allocated()){
      if(!this->weights_allocated())
	this->allocate_weights(0);
      for(int i=0; i<nelem; i++)
	pixelwts[i] = pixarr.wt(i);      
    }
  }

  template<class T>
    void pixel_array<T>::print_axes(ostream & os) const {
    os << "pixel_array::print_axes - axes size " << axes.size() << endl;
    for(uint i=0; i<axes.size(); i++)
      os << "pixel_array::print_axes - axis " << i << " size " << axes[i] << endl;
  }

  template<class T>
    void pixel_array<T>::min_and_max(double & min, double & max) const {
    vector<int> minpix(2,0), maxpix(2,0);
    this->min_and_max(min, minpix, max, maxpix);
  }
  
  template<class T>
    void pixel_array<T>::min_and_max(double & min, vector<int> & minpixel, 
				     double & max, vector<int> & maxpixel) const {
  
    vector<int> axis_0_limits(2,0);
    vector<int> axis_1_limits(2,0);
    axis_0_limits[1] = axes[0]-1;
    axis_1_limits[1] = axes[1]-1;
    this->min_and_max(min, minpixel,
		      max, maxpixel,
		      axis_0_limits, axis_1_limits);
  }

  template<class T>
    void pixel_array<T>::min_and_max(double & min, vector<int> & minpixel, 
				     double & max, vector<int> & maxpixel,
				     vector<int> axis_0_limits, 
				     vector<int> axis_1_limits) const {
  
    if(axis_1_limits.size()!=2 || axis_1_limits.size()!=2){
      cerr << "pixel_array<T>::min_and_max error - axis limits wrong dimension\n";
      throw(string("pixel_array<T>::min_and_max"));
    }
  
    if(axis_0_limits[0] < 0 || axis_0_limits[1] >= this->axes[0]){
      cerr << "pixel_array<T>::min_and_max error - axis 0 limits out of range\n";
      throw(string("pixel_array<T>::min_and_max"));
    }
    if(axis_1_limits[0] < 0 || axis_1_limits[1] >= this->axes[1]){
      cerr << "pixel_array<T>::min_and_max error - axis 1 limits out of range\n";
      throw(string("pixel_array<T>::min_and_max"));
    }
  
    if(verbose_level) 
      cout << "pixel_array<T>::min_and_max - "
           << "searching for min and max within limits " 
	   << "(" << axis_0_limits[0] << "," << axis_0_limits[1] << ")" 
	   << "(" << axis_1_limits[0] << "," << axis_1_limits[1] << ")" << endl;
  
    int index;
    int axes_0_limits_zero = axis_0_limits[0];
    int axes_0_limits_one = axis_0_limits[1];
    int axes_1_limits_zero = axis_1_limits[0];
    int axes_1_limits_one = axis_1_limits[1];
    int axes_zero = axes[0];
  
    minpixel.resize(2);
    maxpixel.resize(2);
    min = pixeldata[0];
    max = pixeldata[0];
    if(pixelwts!=NULL){
      for(int i=axes_1_limits_zero; i<axes_1_limits_one; i++){
	for(int j=axes_0_limits_zero; j<axes_0_limits_one; j++){
	  index = i*axes_zero+j;
	  if(min>pixeldata[index] && pixelwts[index]!=0){
	    min = pixeldata[index];
	    minpixel[0] = j; minpixel[1] = i;
	  }
	  if(max<pixeldata[index] && pixelwts[index]!=0){
	    max = pixeldata[index];
	    maxpixel[0] = j; maxpixel[1] = i;
	  }
	}
      }
    } else {
      for(int i=axes_1_limits_zero; i<axes_1_limits_one; i++){
	for(int j=axes_0_limits_zero; j<axes_0_limits_one; j++){
	  index = i*axes_zero+j;
	  if(min>pixeldata[index]){
	    min = pixeldata[index];
	    minpixel[0] = j; minpixel[1] = i;
	  }
	  if(max<pixeldata[index]){
	    max = pixeldata[index];
	    maxpixel[0] = j; maxpixel[1] = i;
	  }
	}
      }
    }
  }

  template<class T>
    void pixel_array<T>::flip_x(){
    double tmp;
    for(int i=0; i<axes[1]; i++){
      for(int j=0; j<axes[0]/2; j++){
	tmp = pixeldata[i*axes[0]+j];
	pixeldata[i*axes[0]+j] = pixeldata[i*axes[0]+axes[0]-j-1];
	pixeldata[i*axes[0]+axes[0]-j-1] = tmp;
      }
    }
  }
  
  template<class T>
    void pixel_array<T>::flip_y(){
    double tmp;
    for(int i=0; i<axes[1]/2; i++){
      for(int j=0; j<axes[0]; j++){
	tmp = pixeldata[i*axes[0]+j];
	pixeldata[i*axes[0]+j] = pixeldata[(axes[1]-i-1)*axes[0]+j];
	pixeldata[(axes[1]-i-1)*axes[0]+j] = tmp;
      }
    }
  }
  
  template<class T>
    void pixel_array<T>::flip_xy(){
    double tmp;
    for(int i=0; i<axes[1]/2; i++){
      for(int j=0; j<axes[0]/2; j++){
	tmp = pixeldata[i*axes[0]+j];
	pixeldata[i*axes[0]+j] = pixeldata[(axes[1]-i-1)*axes[0] + (axes[0]-j-1)];
	pixeldata[(axes[1]-i-1)*axes[0] + (axes[0]-j-1)] = tmp;
      
	tmp = pixeldata[i*axes[0]+(axes[0]-j-1)];
	pixeldata[i*axes[0]+(axes[0]-j-1)] = pixeldata[(axes[1]-i-1)*axes[0]+j];
	pixeldata[(axes[1]-i-1)*axes[0]+j] = tmp;
      }
    }
  }

  template<class T>
    void pixel_array<T>::flip_45(){
    double tmp;

    if(axes[1]!=axes[0]){
      cerr << "pixel_array::flip_45 error - cannot flip an array that isn't square\n"
	   << "dimens "
	   << axes[0] << "x" << axes[1] << endl;
      throw(string("pixel_array::flip_45"));
    }

    for(int i=0; i<axes[1]; i++){
      for(int j=i; j<axes[0]; j++){
	tmp = pixeldata[i*axes[0]+j];
	pixeldata[i*axes[0]+j] = pixeldata[j*axes[0]+i];
	pixeldata[j*axes[0]+i] = tmp;
      }
    }
  }
  
  template<class T>
    void pixel_array<T>::shift_by_fft(double dx, double dy){ 

    if(dx==0 && dy==0) return;
 
    if(axes.size()!=2){
      cerr << "pixel_array::shift_by_fft error - "
           << "array does not have 2 dimensions\n"; 
      throw(string("pixel_array::shift_by_fft")); 
    }  
    if(axes[0]<=1 || axes[1]<=1){
      cerr << "pixel_array::shift_by_fft error - "
           << "array doesn't contain enough elements\n";
      throw(string("pixel_array::shift_by_fft")); 
    }  
      
    if(fabs(dx)>=1 || fabs(dy)>=1){
      cerr << "pixel_array::shift_by_fft error - "
           << "cannot shift by more than 1 pixel\n"; 
      throw(string("pixel_array::shift_by_fft")); 
    }

    if(pixel_array<T>::verbose_level)
      cout << "pixel_array::shift_by_fft - shifting " << dx << " pixels in x, "
	   << dy << " pixels in y\n";

    // array must be padded by an extra complex element  
    // for the real to complex fft 
    int nelem = axes[1]*2*(axes[0]/2+1); 
  
    T * thisarr;
    T * shiftarr;
    try{ 
      thisarr = new T[nelem]; 
      shiftarr = new T[2*axes[0]*axes[1]];
    } catch(...) { 
      cerr << "pixel_array::shift_by_fft error - error allocating memory\n"; 
      throw(string("pixel_array::shift_by_fft")); 
    }  
    Arroyo::rfft_manager<T> rfft_mgr;
    Arroyo::fft_manager<T> fft_mgr;

    int axes_zero = axes[0];
    int axes_one = axes[1];
    for(int i=0; i<axes_one; i++)
      for(int j=0; j<axes_zero; j++)
	thisarr[i*2*(axes_zero/2+1)+j] = pixeldata[i*axes_zero+j]; 

    if(pixel_array<T>::verbose_level)
      cout << "pixel_array::shift_by_fft - performing forward fft\n";

    //  SWITCHING FFT SCHEMES TO FFTW_MANAGER
    // fftrcf2d(axes[1], axes[0], thisarr); 
    // for the fortran ordering...
    bool fftw_estimate = true;
    bool fftw_in_place = true;
    vector<long> flipped_axes(2);
    flipped_axes[0] = axes[1];
    flipped_axes[1] = axes[0];
    rfft_mgr.real_to_complex_fft(flipped_axes,
				 fftw_estimate, fftw_in_place,
				 thisarr);

    double dxs = 2*M_PI*dx/(float)axes[0];
    double dys = 2*M_PI*dy/(float)axes[1];
    double phase, tmp, cp, sp;
    // positive frequencies
    int index, indexa, indexb;
    for(int i=0; i<axes_one; i++){
      for(int j=0; j<axes_zero/2+1; j++){
	indexa = i*axes_zero + j;
	indexb = i*(axes_zero/2+1) + j;
	shiftarr[2*indexa] = thisarr[2*indexb];
	shiftarr[2*indexa+1] = thisarr[2*indexb+1];
      }
    }

    // negative frequencies
    for(int j=axes_zero/2+1; j<axes_zero; j++){
      indexa = j;
      indexb = axes_zero-j;
      shiftarr[2*indexa] = thisarr[2*indexb];
      shiftarr[2*indexa+1] = -thisarr[2*indexb+1];
    }      

    for(int i=1; i<axes_one; i++){
      for(int j=axes_zero/2+1; j<axes_zero; j++){
	indexa = i*axes_zero + j;
	indexb = (axes_one-i)*(axes_zero/2+1) + axes_zero - j;
	shiftarr[2*indexa] = thisarr[2*indexb];
	shiftarr[2*indexa+1] = -thisarr[2*indexb+1];
      }
    }

    for(int i=0; i<axes_one; i++){
      for(int j=0; j<axes_zero; j++){
	if(j<axes_zero/2) phase = dxs*j;  
	else phase = -dxs*(axes_zero-j);
	if(i<axes_one/2) phase += dys*i;  
	else phase -= dys*(axes_one-i);  
	cp = cos(-phase);
	sp = sin(-phase);
	index = 2*(i*axes_zero+j);
	tmp = shiftarr[index]*cp-shiftarr[index+1]*sp;
	shiftarr[index+1] = shiftarr[index]*sp + shiftarr[index+1]*cp;
	shiftarr[index] = tmp;
      }
    }

    if(pixel_array<T>::verbose_level)
      cout << "pixel_array::shift_by_fft - performing backwards fft\n";
    //  SWITCHING FFT SCHEMES TO FFTW_MANAGER
    //  fftccb2d(axes[1], axes[0], shiftarr);
    //  original line - fftccb2d replaced it: fftcrb2dip(axes[1], axes[0], thisarr);
    fft_mgr.backward_fft(flipped_axes, fftw_estimate, fftw_in_place, shiftarr);

    double scale = 1/(double)(axes[0]*axes[1]);
    for(int i=0; i<axes_one; i++)
      for(int j=0; j<axes_zero; j++){
	pixeldata[i*axes_zero+j] =
	  static_cast<T>(shiftarr[2*(i*axes_zero+j)]*scale);
      }

    delete [] thisarr;
    delete [] shiftarr;
  }

  template<class T>
    void pixel_array<T>::rotate_by_fft(double angle, bool window){
    if(angle==0) return;
    this->rotate_and_shift_by_fft(0, 0, angle, window);
  }

  template<class T>
    void pixel_array<T>::rotate_and_shift_by_fft(
						 double dx, double dy, double angle, bool window){
    this->simple_rotate_and_shift_by_fft(dx, dy, angle, window);
  }

  template<class T>
    void pixel_array<T>::simple_rotate_and_shift_by_fft(
							double dx, double dy, double angle, bool window){

    // not yet supported
    if(dx!=0 || dy!=0){
      cerr << "pixel_array::simple_rotate_and_shift_by_fft error - "
	   << "shifting by fft not yet supported by this function\n";
      throw(string("pixel_array::simple_rotate_and_shift_by_fft"));
    }
  
    // probably can never be supported
    if(this->weights_allocated()){
      cerr << "pixel_array::simple_rotate_and_shift_by_fft error - "
	   << "this member function does not support rotation of "
	   << "arrays with weights defined\n";
      throw(string("pixel_array::simple_rotate_and_shift_by_fft"));
    }

    angle = fmod(angle, 2*M_PI);
    if(angle<0) angle += 2*M_PI;

    if(pixel_array<T>::verbose_level)
      cout << "pixel_array::simple_rotate_and_shift_by_fft - "
           << "rotating through angle " 
	   << angle << endl;

    // mod the angle into the range 0 < angle < M_PI
    if(angle>=3*M_PI_2){
      // here the pixel_array has the correct amount of memory
      // allocated, but we have to reshuffle to switch the indexing
      int nelem = axes[0]*axes[1];
      T * tmparray = new T[nelem];
      for(int i=0; i<axes[1]; i++)
	for(int j=0; j<axes[0]; j++)
	  tmparray[(axes[0]-j-1)*axes[1]+i] = this->pixeldata[i*axes[0]+j];
      long tmp = axes[0];
      axes[0] = axes[1];
      axes[1] = tmp;
      for(int i=0; i<nelem; i++)
	this->pixeldata[i] = tmparray[i];
      delete [] tmparray;
      angle -= 3*M_PI_2;
    } else if(angle>=M_PI) {
      this->flip_xy();
      angle -= M_PI;
    } else if(angle>=M_PI_2) {
      // as above for angle>=3*M_PI_2
      int nelem = axes[0]*axes[1];
      T * tmparray = new T[nelem];
      for(int i=0; i<axes[1]; i++)
	for(int j=0; j<axes[0]; j++)
	  tmparray[j*axes[1]+(axes[1]-i-1)] = this->pixeldata[i*axes[0]+j];
      long tmp = axes[0];
      axes[0] = axes[1];
      axes[1] = tmp;
      for(int i=0; i<nelem; i++)
	this->pixeldata[i] = tmparray[i];
      delete [] tmparray;
      angle -= M_PI_2;
    }

    // if we're down to zero angle, we can just return
    if(dx==0 && dy==0 && angle==0) return;

    if(pixel_array<T>::verbose_level)
      cout << "pixel_array::simple_rotate_and_shift_by_fft - "
	   << "performing rotation by fft through angle " << angle << endl;

    // define new axes as the minimal size required to contain the rotated image
    vector<long> new_axes(2);
    new_axes[0] = (int)(ceil(axes[0] + axes[1]*tan(angle/2.0)));
    new_axes[1] = (int)(2*(ceil(axes[1] + axes[0]*sin(angle))/2));

    // do this again to account for the two row shears
    new_axes[0] = (int)(2*(ceil(new_axes[0] + axes[1]*tan(angle/2.0))/2));

    if(pixel_array<T>::verbose_level)
      cout << "pixel_array::simple_rotate_and_shift_by_fft - "
	   << "old axes " << axes[0] << "\t" << axes[1]
	   << " new axes " << new_axes[0] << "\t" << new_axes[1] << endl;

    int row_pad = new_axes[0]-axes[0];
    int column_pad = new_axes[1]-axes[1];

    int nelem = 2*(new_axes[0]/2+1)*new_axes[1];
    int column_stride = 2*(new_axes[0]/2+1);

    double * tmparr;
    try{tmparr = new double[nelem];} 
    catch(...) {
      cerr << "pixel_array::simple_rotate_and_shift_by_fft error - "
           << "could not allocate memory\n";
      throw(string("pixel_array::simple_rotate_and_shift_by_fft"));
    }

    // Copy the pixeldata into the temporary array
    for(int i=0; i<nelem; i++)
      tmparr[i] = 0;

    for(int i=0; i<axes[1]; i++)
      for(int j=0; j<axes[0]; j++)
	tmparr[(i+column_pad/2)*column_stride + j + row_pad/2]
	  = pixeldata[i*axes[0]+j];

    bool fftw_estimate = false;
    bool fftw_in_place = true;

    Arroyo::rfft_manager<double> rfft_row_mgr;

    int index;
    double dyslope, shear, phase, cp, sp, tmp;
    double amp, tphase;
    for(int i=0; i<new_axes[1]; i++){
      rfft_row_mgr.real_to_complex_fft(vector<long>(1,new_axes[0]),
				       fftw_estimate, fftw_in_place,
				       &(tmparr[i*column_stride]));
      //shear = -tan(angle/2.0)*(i-new_axes[1]/2);
      shear = tan(angle/2.0)*(i-new_axes[1]/2);
      dyslope = 2*M_PI*shear/(float)new_axes[0];
      for(int j=0; j<column_stride/2; j++){
	phase = dyslope*j;
	if(window){
	  cp = (1-(j/(double)(column_stride/2+1)))*cos(phase)/(double)column_stride;
	  sp = (1-(j/(double)(column_stride/2+1)))*sin(phase)/(double)column_stride;
	} else {
	  cp = cos(phase)/(double)column_stride;
	  sp = sin(phase)/(double)column_stride;
	}
	index = i*column_stride+2*j;
	tmp = tmparr[index]*cp-tmparr[index+1]*sp;      
	tmparr[index+1] = tmparr[index]*sp + tmparr[index+1]*cp;
	tmparr[index] = tmp;
      }
      rfft_row_mgr.complex_to_real_fft(vector<long>(1,new_axes[0]),
				       fftw_estimate, fftw_in_place,
				       &(tmparr[i*column_stride]));
    } 

    Arroyo::rfft_manager<double> rfft_column_mgr;

    double * tmp_1d_array = new double[2*(new_axes[1]/2+1)];  
    for(int i=0; i<new_axes[0]; i++){
      for(int j=0; j<new_axes[1]; j++)
	tmp_1d_array[j] = tmparr[j*column_stride+i];
      rfft_column_mgr.real_to_complex_fft(vector<long>(1,new_axes[1]),
					  fftw_estimate, fftw_in_place,
					  tmp_1d_array);

      //shear = sin(angle)*(i-new_axes[0]/2);
      shear = -sin(angle)*(i-new_axes[0]/2);
      dyslope = 2*M_PI*shear/(float)new_axes[1];
      for(int j=0; j<new_axes[1]/2+1; j++){
	phase = dyslope*j;
	if(window){
	  cp = (1-(j/(double)(new_axes[1]/2+1)))*cos(phase)/(double)new_axes[1];
	  sp = (1-(j/(double)(new_axes[1]/2+1)))*sin(phase)/(double)new_axes[1];
	} else {
	  cp = cos(phase)/(double)new_axes[1];
	  sp = sin(phase)/(double)new_axes[1];
	} 
	tmp = tmp_1d_array[2*j]*cp-tmp_1d_array[2*j+1]*sp;      
	tmp_1d_array[2*j+1] = tmp_1d_array[2*j]*sp + tmp_1d_array[2*j+1]*cp;
	tmp_1d_array[2*j] = tmp;
      }

      rfft_column_mgr.complex_to_real_fft(vector<long>(1,new_axes[1]),
					  fftw_estimate, fftw_in_place,
					  tmp_1d_array);
      for(int j=0; j<new_axes[1]; j++)
	tmparr[j*column_stride+i] = tmp_1d_array[j];
    }
    delete [] tmp_1d_array;

    for(int i=0; i<new_axes[1]; i++){
      rfft_row_mgr.real_to_complex_fft(vector<long>(1,new_axes[0]),
				       fftw_estimate, fftw_in_place,
				       &(tmparr[i*column_stride]));
      //shear = -tan(angle/2.0)*(i-new_axes[1]/2);
      shear = tan(angle/2.0)*(i-new_axes[1]/2);
      dyslope = 2*M_PI*shear/(float)new_axes[0];
      for(int j=0; j<column_stride/2; j++){
	phase = dyslope*j;
	if(window){
	  cp = (1-(j/(double)(column_stride/2+1)))*cos(phase)/(double)column_stride;
	  sp = (1-(j/(double)(column_stride/2+1)))*sin(phase)/(double)column_stride;
	} else {
	  cp = cos(phase)/(double)column_stride;
	  sp = sin(phase)/(double)column_stride;
	}
	index = i*column_stride+2*j;
	tmp = tmparr[index]*cp-tmparr[index+1]*sp;      
	tmparr[index+1] = tmparr[index]*sp + tmparr[index+1]*cp;
	tmparr[index] = tmp;
      }
      rfft_row_mgr.complex_to_real_fft(vector<long>(1,new_axes[0]),
				       fftw_estimate, fftw_in_place,
				       &(tmparr[i*column_stride]));
    } 

    this->set_axes(new_axes);
    for(int i=0; i<new_axes[1]; i++)
      for(int j=0; j<new_axes[0]; j++){
	this->pixeldata[i*new_axes[0]+j] = tmparr[i*column_stride+j];
      }
    delete [] tmparr;  
  }

  template<class T>
    void pixel_array<T>::owen_makedon_rotate_and_shift_by_fft(
							      double dx, double dy, double angle){

    // Change this to account for different angles by flipping in x or y
    if(angle < -M_PI_2 || angle > M_PI_2){
      cerr << "pixel_array::owen_makedon_rotate_and_shift_by_fft error - "
	   << "requested rotation angle "
	   << angle << " outside the range -PI/2 < angle < PI/2\n";
      throw(string("pixel_array::rotate"));
    }

    // Compute the max shear
    double ta = tan(angle/2.0);
    long max_shear = (long)(ceil(axes[1]*ta));
    if(max_shear%1!=0) max_shear++;

    cout << "pixel_array::owen_makedon_rotate_and_shift_by_fft - max shear "
    	 << max_shear 
	 << " for angle " << angle << " axes * tan angle / 2 "
	 << axes[1]*ta << endl;

    // allocate a temporary array
    vector<long> new_axes(2);
    new_axes[0] = (axes[0]+2*max_shear);
    new_axes[1] = 2*axes[1];

    cout << "pixel_array::owen_makedon_rotate_and_shift_by_fft - new axes "
	 << new_axes[0] << "\t" << new_axes[1] << endl;

    int nelem = 2*(new_axes[0]/2+1)*new_axes[1];
    int column_stride = 2*(new_axes[0]/2+1);

    cout << "pixel_array::owen_makedon_rotate_and_shift_by_fft - allocating "
         << nelem << " elements\n";

    double * tmparr;
    try{tmparr = new double[nelem];} 
    catch(...) {
      cerr << "pixel_array::owen_makedon_rotate_and_shift_by_fft error - "
           << "could not allocate memory\n";
      throw(string("pixel_array::owen_makedon_rotate_and_shift_by_fft"));
    }

    // Copy the pixeldata into the temporary array
    for(int i=0; i<nelem; i++)
      tmparr[i] = 0;

    for(int i=0; i<axes[1]; i++)
      for(int j=0; j<axes[0]; j++)
	tmparr[(i+axes[1]/2)*column_stride+j+max_shear] = pixeldata[i*axes[0]+j];

    Arroyo::rfft_manager<double> rfft_row_mgr;
    bool fftw_estimate = false;
    bool fftw_in_place = true;

    // Following the paper "High quality alias free image rotation" by
    // Owen & Makedon.  Dartmouth Dept. of Computer Science PCS-TR96-301

    // Step 1
    // Here we follow a slightly different strategy, by allocating all the space
    // we will need up front, and then picking through the array to do the correct
    // transforms

    // STEPS 2 and 3
    // shift each row by the appropriate sheared amount
    int index;
    double dyslope, shear, phase, cp, sp, tmp;
    double amp, tphase;
    for(int i=axes[1]/2; i<3*axes[1]/2; i++){
      rfft_row_mgr.real_to_complex_fft(vector<long>(1,new_axes[0]),
				       fftw_estimate, fftw_in_place,
				       &(tmparr[i*column_stride]));
      shear = -tan(angle/2.0)*(i-axes[1]);
      dyslope = 2*M_PI*shear/(float)new_axes[0];
      for(int j=0; j<column_stride/2; j++){
	phase = dyslope*j;
	cp = cos(phase);
	sp = sin(phase);
	index = i*column_stride+2*j;
	tmp = tmparr[index]*cp-tmparr[index+1]*sp;      
	tmparr[index+1] = tmparr[index]*sp + tmparr[index+1]*cp;
	tmparr[index] = tmp;
      }
      /*
      //rfft_row_mgr.complex_to_real_fft(vector<long>(1,new_axes[0]),
      fftw_estimate, fftw_in_place,
      &(tmparr[i*column_stride]));
      */
    } 

    Arroyo::fft_manager<double> fft_column_mgr;

    // STEP 4 
    // Transform columns
    double * tmp_1d_array = new double[2*new_axes[1]];  
    for(int i=0; i<2*axes[1]; i++) 
      tmp_1d_array[i] = 0;

    for(int i=0; i<new_axes[0]/2; i++){
      for(int j=0; j<axes[1]; j++){
	tmp_1d_array[2*j] = tmparr[(axes[1]/2+j)*column_stride+2*i];
	tmp_1d_array[2*j+1] = tmparr[(axes[1]/2+j)*column_stride+2*i+1];
      }
      fft_column_mgr.forward_fft(vector<long>(1,axes[1]), false, true, tmp_1d_array);
      for(int j=0; j<axes[1]; j++){
	tmparr[(axes[1]/2+j)*column_stride+2*i] = tmp_1d_array[2*j];
	tmparr[(axes[1]/2+j)*column_stride+2*i+1] = tmp_1d_array[2*j+1];
      }
    }

    // Intermediate step - do the shift by fft by adding phase

    // STEP 5
    // duplicate rows out of band
    for(int i=0; i<axes[1]/2; i++)
      for(int j=0; j<2*(new_axes[0]/2+1); j++){
	tmparr[i*column_stride+j] = 
	  tmparr[(axes[1]+i)*column_stride+j];
	tmparr[(3*axes[1]/2+i)*column_stride+j] = 
	  tmparr[(axes[1]/2+i)*column_stride+j];
      }

      
    // STEPS 6, 7, and 8
    // inverse fft rows, shift to effect column shear, and fft back
    for(int i=0; i<new_axes[1]; i++){ 
      rfft_row_mgr.complex_to_real_fft(vector<long>(1,new_axes[0]),
				       fftw_estimate, fftw_in_place,
				       &(tmparr[i*column_stride]));
      shear = -2*sin(angle/2.0)*(i-new_axes[1]/2);
      dyslope = 2*M_PI*shear/(float)new_axes[0];
      for(int j=0; j<column_stride/2; j++){
	phase = dyslope*j;
	cp = cos(phase);
	sp = sin(phase);
	index = i*column_stride+2*j;
	tmp = tmparr[index]*cp-tmparr[index+1]*sp;      
	tmparr[index+1] = tmparr[index]*sp + tmparr[index+1]*cp;
	tmparr[index] = tmp;
      }
      rfft_row_mgr.real_to_complex_fft(vector<long>(1,new_axes[0]),
				       fftw_estimate, fftw_in_place,
				       &(tmparr[i*column_stride]));
    }


    // STEP 9 Window - here we want to determine whether the coordinates
    // of each spectral point lies within the original rectangular
    // region.  To do so, we apply the inverse transformation on the
    // coordinates to determine the original coordinates, and window if
    // these lie outside the rectangle.
    // 
    // Recall that we got to this stage by applying a spatial row shear
    // of -tan theta/2 - equivalent to a spectral column shear of tan
    // theta/2 - followed by a spectral row shear of -2 sin theta.
    // So the inverse transformation is
    // 
    //  ( x' )  =  (1             0) (1   2 sin theta) ( x )
    //  ( y' )     (-tan theta/2  1) (0       1      ) ( y )
    double x_halfpix=0, y_halfpix=0;
    int x_extrapix=1, y_extrapix=1;
    if(axes[1]%2==0){
      x_halfpix = .5;
      x_extrapix = 0;
    }
    if(axes[0]%2==0){
      y_halfpix = .5;
      y_extrapix = 0;
    }
    for(int i=-new_axes[1]/2; i<new_axes[1]/2+x_extrapix; i++){
      for(int j=0; j<new_axes[0]/2+y_extrapix; j++){
	index = (i+new_axes[1]/2)*column_stride+2*j;
	cout << i << "\t" << j << "\t";
	if(fabs( (i+x_halfpix) + 2*sin(angle)*(j+y_halfpix) )>axes[1]/2)
	  tmparr[index] = tmparr[index+1] = 0;
	if(fabs( -tan(angle/2.0) * (i+x_halfpix)
		 + (1-2*sin(angle)*tan(angle/2.0))*(j+y_halfpix) )>axes[0]/2)
	  tmparr[index] = tmparr[index+1] = 0;
	cout << fabs( (i+x_halfpix) + 2*sin(angle)*(j+y_halfpix))
	     << "\t" 
	     << fabs( -tan(angle/2.0) * (i+x_halfpix)
		      + (1-2*sin(angle)*tan(angle/2.0))*(j+y_halfpix))
	     << "\t";
	cout << tmparr[index] << "\t" << tmparr[index+1] << endl;
      }
    }


    // STEPS 10 and 11
    // inverse fft columns, shift to effect row shear
    Arroyo::fft_manager<double> new_column_mgr;

    for(int i=0; i<new_axes[0]/2; i++){
      for(int j=0; j<new_axes[1]; j++){
	tmp_1d_array[2*j] = tmparr[j*column_stride+2*i];
	tmp_1d_array[2*j+1] = tmparr[j*column_stride+2*i+1];
      }

      new_column_mgr.backward_fft(vector<long>(1,new_axes[1]),
				  fftw_estimate, fftw_in_place,
				  tmp_1d_array);

      shear = .5*tan(angle/2.0)*i;
      dyslope = 2*M_PI*shear/(float)new_axes[1];
      for(int j=0; j<new_axes[1]; j++){
	phase = dyslope*(j-new_axes[1]/2);
	cp = cos(phase);
	sp = sin(phase);
	index = 2*j;
	tmp = tmp_1d_array[index]*cp-tmp_1d_array[index+1]*sp;      
	tmp_1d_array[index+1] = tmp_1d_array[index]*sp + tmp_1d_array[index+1]*cp;
	tmp_1d_array[index] = tmp;
      }

      for(int j=0; j<new_axes[1]; j++){
	tmparr[j*column_stride+2*i] = tmp_1d_array[2*j];
	tmparr[j*column_stride+2*i+1] = tmp_1d_array[2*j+1];
      }
    }

    // STEP 12
    // inverse fft rows
    for(int i=0; i<new_axes[1]; i++)
      rfft_row_mgr.complex_to_real_fft(vector<long>(1,new_axes[0]),
				       fftw_estimate, fftw_in_place,
				       &(tmparr[i*column_stride]));

    delete [] tmp_1d_array;

    this->set_axes(new_axes);
    for(int i=0; i<new_axes[1]; i++)
      for(int j=0; j<new_axes[0]; j++){
	this->pixeldata[i*new_axes[0]+j] = tmparr[i*2*(new_axes[0]/2+1)+j];
      }
    delete [] tmparr;  

  }

  template<class T>
    pixel_array<T> pixel_array<T>::cross_correlate(
						   const pixel_array<T> & pixarr) const {
    if(pixarr.axes.size()!=2 || axes.size()!=2
       || axes[0]!=pixarr.axes[0] || axes[1]!=pixarr.axes[1]){
      cerr << "pixel_array::cross_correlate error - array mismatch\n";
      throw(string("pixel_array::cross_correlate"));
    }
    if(axes[0]<=1 || axes[1]<=1 || pixarr.axes[0]<=1 || pixarr.axes[1]<=1){
      cerr << "pixel_array::cross_correlate error - "
	   << "array doesn't contain enough elements\n";
      throw(string("pixel_array::cross_correlate")); 
    }
    
    // arrays must be padded by an extra complex element 
    // for the real to complex fft
    int nelem = axes[1]*2*(axes[0]/2+1);

    T * fft_xcorr;
    T * xcorr;
    T * thisarr;
    T * inarr;
    try{
      if(pixel_array<T>::verbose_level)
	cout << "pixel_array::cross_correlate - allocating data nelem "
	     << nelem << endl;
      thisarr = new T[nelem];
      inarr = new T[nelem];
      fft_xcorr = new T[nelem];
      if(pixel_array<T>::verbose_level)
	cout << "pixel_array::cross_correlate - allocating data nelem "
	     << axes[0]*axes[1] << endl;
      xcorr = new T[axes[0]*axes[1]];
    } catch(...) {
      cerr << "pixel_array::cross_correlate error - error allocating memory\n";
      throw(string("pixel_array::cross_correlate"));
    }
    Arroyo::rfft_manager<T> rfft_mgr;

    int axes_one = axes[1];
    int axes_zero = axes[0];
    for(int i=0; i<axes_one; i++){
      for(int j=0; j<axes_zero; j++){
	thisarr[i*2*(axes_zero/2+1)+j] = this->pixeldata[i*axes_zero+j];
	inarr[i*2*(axes_zero/2+1)+j] = pixarr.pixeldata[i*axes_zero+j];
      }
    }

    //  SWITCHING FFT SCHEMES TO FFTW_MANAGER
    //fftrcf2d(axes[1], axes[0], thisarr);
    //fftrcf2d(axes[1], axes[0], inarr);
    bool fftw_estimate = true;
    bool fftw_in_place = true;
    vector<long> flipped_axes(2);
    flipped_axes[0] = axes[1];
    flipped_axes[1] = axes[0];
    if(pixel_array<T>::verbose_level)
      cout << "pixel_array::cross_correlate - transforming this pixel_array\n";
    rfft_mgr.real_to_complex_fft(flipped_axes,
				 fftw_estimate, fftw_in_place,
				 thisarr); 
    if(pixel_array<T>::verbose_level)
      cout << "pixel_array::cross_correlate - transforming arg pixel_array\n";
    rfft_mgr.real_to_complex_fft(flipped_axes,
				 fftw_estimate, fftw_in_place,
				 inarr);

    int indexa, indexb;
    double scale = 1.0/(double)(axes[0]*axes[1]);
    for(int i=0; i<axes_one; i++){
      for(int j=0; j<axes_zero/2+1; j++){
	indexa = i*2*(axes_zero/2+1) + j;
	indexb = i*2*(axes_zero/2+1) + 2*j;
	fft_xcorr[indexb] = (thisarr[indexb]*inarr[indexb] + 
			     thisarr[indexb+1]*inarr[indexb+1])*scale;
	fft_xcorr[indexb+1] = (thisarr[indexb]*inarr[indexb+1] - 
			       thisarr[indexb+1]*inarr[indexb])*scale;
      }
    }
    delete [] thisarr;
    delete [] inarr;

    // Perform inverse fft to obtain autocorrelation function
    //  SWITCHING FFT SCHEMES TO FFTW_MANAGER
    //fftcrb2d(axes[1],axes[0],fft_xcorr, xcorr);
    fftw_in_place = false;
    if(pixel_array<T>::verbose_level)
      cout << "pixel_array::cross_correlate - "
	   << "back-transforming to retrieve cross-correlated array\n";
    rfft_mgr.complex_to_real_fft(flipped_axes,
				 fftw_estimate, fftw_in_place,
				 fft_xcorr, xcorr);
  
    delete [] fft_xcorr;

    pixel_array<T> xcorr_pixarr = pixel_array<T>(axes, xcorr);
    delete [] xcorr;
    return(xcorr_pixarr);
  }
  
  template<class T>
    void pixel_array<T>::offset(const pixel_array<T> & pixarr, 
				vector<double> & offsets, 
				long range) const {
    
    if(pixarr.axes.size()!=2 || axes.size()!=2
       || axes[0]!=pixarr.axes[0] || axes[1]!=pixarr.axes[1]){
      cerr << "pixel_array::offset error - array mismatch\n";
      throw(string("pixel_array::offset"));
    }
    if(axes[0]<=1 || axes[1]<=1 || pixarr.axes[0]<=1 || pixarr.axes[1]<=1){
      cerr << "pixel_array::offset error - "
	   << "array doesn't contain enough elements\n";
      throw(string("pixel_array::offset")); 
    } 

    if(pixel_array<T>::verbose_level)
      cout << "pixel_array::offset - seeking offset between arrays of size\n"
	   << "\t(" << this->axes[0] << ", " << this->axes[1]  
	   << ") and (" << pixarr.axes[0] << "," << pixarr.axes[1] 
	   << ")" << endl;

    pixel_array<T> xcorr_pixarr = this->cross_correlate(pixarr);
    vector<int> xcorr_minpixel(2,0), xcorr_maxpixel(2,0);

    double min, max;
    xcorr_pixarr.min_and_max(min, xcorr_minpixel, max, xcorr_maxpixel);

    // make the offset (-axis/2,axis/2] rather than [0,axis-1] 
    for(int i=0; i<2; i++){
      if(xcorr_maxpixel[i]>axes[i]/2){
	if(pixel_array<T>::verbose_level) 
	  cout << "pixel_array<T>::offset - correcting dec " << i << endl;
	xcorr_maxpixel[i] -= axes[i];
      }
      if(xcorr_maxpixel[i]<-axes[i]/2){
	if(pixel_array<T>::verbose_level) 
	  cout << "pixel_array<T>::offset - correcting axis " << i << endl;
	xcorr_maxpixel[i] += axes[i];
      }
    }

    // Find the maximum of this pixel array within limits
    // These limits ensure that the maximum lies at least
    // <range points> from the edges of the array

    vector<int> tmp_minpixel, tmp_maxpixel;
    if(range!=-1){
      vector<int> axis_0_range(2,range), axis_1_range(2,range);
      axis_0_range[1] = axes[0]-range-1;
      axis_1_range[1] = axes[1]-range-1;
      
      this->min_and_max(min, tmp_minpixel, max, tmp_maxpixel, axis_0_range, axis_1_range);
    } else 
      this->min_and_max(min, tmp_minpixel, max, tmp_maxpixel);
    
    // Ensure that the maximum is not shifted off the edge by xcorr_maxpixel
    for(int i=0; i<2; i++){
      if(tmp_maxpixel[i]+xcorr_maxpixel[i]>axes[i])
	xcorr_maxpixel[i] -= axes[i];
      if(tmp_maxpixel[i]-xcorr_maxpixel[i]<0)
	xcorr_maxpixel[i] += axes[i];
    }

    if(pixel_array<T>::verbose_level) {
      cout << "pixel_array<T>::offset - xcorr max " << xcorr_maxpixel[0] << "\t" 
	   << xcorr_maxpixel[1] << "\t" << max << endl;
      cout << "pixel_array<T>::offset - this max " << max << " at " 
	   << tmp_maxpixel[0] << "," << tmp_maxpixel[1] << endl;
    }

    pixel_array<T> tpix, apix;
    if(range!=-1){
      // Extract a region from this pixel array around the 
      // maximum set by the maxrange argument
	
      long maxrange = range;
      vector<long> this_pixlimits(4), arg_pixlimits(4);
	
      this_pixlimits[0] = tmp_maxpixel[0] - maxrange;
      this_pixlimits[1] = tmp_maxpixel[0] + maxrange-1;
      this_pixlimits[2] = tmp_maxpixel[1] - maxrange;
      this_pixlimits[3] = tmp_maxpixel[1] + maxrange-1;
	
      if(pixel_array<T>::verbose_level) 
	cout << "pixel_array<T>::offset - extracting pixels " 
	     << this_pixlimits[0] << " - " << this_pixlimits[1] << " in x " 
	     << this_pixlimits[2] << " - " << this_pixlimits[3] << " in y "
	     << this_pixlimits[1] - this_pixlimits[0] + 1 << " x " 
	     << this_pixlimits[3] - this_pixlimits[2] + 1 << endl;
      tpix = pixel_array<T>(*this, this_pixlimits);
	
      // Extract a region from the argument pixel array around the 
      // maximum set by the maxrange argument
      /*
	arg_pixlimits[0] = (this_pixlimits[0]+xcorr_maxpixel[0]+axes[0])%axes[0];
	arg_pixlimits[1] = (this_pixlimits[1]+xcorr_maxpixel[0]+axes[0])%axes[0];
	arg_pixlimits[2] = (this_pixlimits[2]+xcorr_maxpixel[1]+axes[1])%axes[1];
	arg_pixlimits[3] = (this_pixlimits[3]+xcorr_maxpixel[1]+axes[1])%axes[1];
      */
	
      arg_pixlimits[0] = this_pixlimits[0]+xcorr_maxpixel[0];
      arg_pixlimits[1] = this_pixlimits[1]+xcorr_maxpixel[0];
      arg_pixlimits[2] = this_pixlimits[2]+xcorr_maxpixel[1];
      arg_pixlimits[3] = this_pixlimits[3]+xcorr_maxpixel[1];
	
      if(pixel_array<T>::verbose_level) 
	cout << "pixel_array<T>::offset - extracting pixels " 
	     << arg_pixlimits[0] << " - " << arg_pixlimits[1] << " in x " 
	     << arg_pixlimits[2] << " - " << arg_pixlimits[3] << " in y " 
	     << arg_pixlimits[1] - arg_pixlimits[0] + 1 << " x " 
	     << arg_pixlimits[3] - arg_pixlimits[2] + 1 << endl;
      apix = pixel_array<T>(pixarr, arg_pixlimits);
    } else {
      tpix = *this;
      apix = pixarr;
    }
	
    tpix.min_and_max(min, tmp_minpixel, max, tmp_maxpixel);

    if(pixel_array<T>::verbose_level)
      cout << "pixel_array<T>::offset - tpix max " << max << " at " 
	   << tmp_maxpixel[0] << "," << tmp_maxpixel[1] << endl;
    apix.min_and_max(min, tmp_minpixel, max, tmp_maxpixel);
    if(pixel_array<T>::verbose_level)
      cout << "pixel_array<T>::offset - apix max " << max << " at " 
	   << tmp_maxpixel[0] << "," << tmp_maxpixel[1] << endl;

    // Convert from amp phase to product of amplitudes and phase differences
    // At the same time, resort the array to reconstruct the correct ordering
    // of elements.  For a real to complex fft, we have the symmetry
    // X[a,b] = X*[n_a - a, n_b - b]
    // Therefore, we must duplicate the array (with the exception of DC and,
    // for n even, nyquist)
    // by flipping it about the horizontal (symmetry above) and the vertical
    // (to get it in order of increasing frequency) and append it to the front of 
    // the original array
    // We do this on the arrays with the zero padding defined above, neglecting the
    // last complex element (the zero padding for the fft

    // Treat the rows in pairs, first/last, first-1/last-1...
    // First form the amplitude product and the phase difference
    // of the positive frequencies
    // Then form these for the negative frequencies stored in the
    // other array of the pair
    // This stage must take account of whether the array has an
    // even or odd number of elements
    // Finally, copy them into the original array
    int nelem = tpix.axes[1]*2*(tpix.axes[0]/2+1);

    T * tarr;
    T * aarr;
    T * product_amps;
    T * phase_difference;
    try{
      if(pixel_array<T>::verbose_level)
	cout << "pixel_array::offset - allocating data nelem " << nelem << endl;
      tarr = new T[nelem];
      aarr = new T[nelem];
      product_amps = new T[nelem];
      phase_difference = new T[nelem];
    } catch(...) {
      cerr << "pixel_array::offset error - error allocating memory\n";
      throw(string("pixel_array::offset"));
    }
    Arroyo::rfft_manager<T> rfft_mgr;  

    int scale = tpix.axes[0]*tpix.axes[1];
    int tpix_axes_one = tpix.axes[1];
    int tpix_axes_zero = tpix.axes[0];
    for(int i=0; i<tpix_axes_one; i++){
      for(int j=0; j<tpix_axes_zero; j++){
	tarr[i*2*(tpix_axes_zero/2+1)+j] = tpix.pixeldata[i*tpix_axes_zero+j];
	aarr[i*2*(tpix_axes_zero/2+1)+j] = apix.pixeldata[i*tpix_axes_zero+j];
      }
    }

    //  SWITCHING FFT SCHEMES TO FFTW_MANAGER
    //fftrcf2d(tpix.axes[1], tpix.axes[0], tarr);
    //fftrcf2d(apix.axes[1], apix.axes[0], aarr);
    bool fftw_estimate = true;
    bool fftw_in_place = true;
    vector<long> flipped_axes(2);
    flipped_axes[0] = tpix.axes[1];
    flipped_axes[1] = tpix.axes[0];
    rfft_mgr.real_to_complex_fft(flipped_axes,
				 fftw_estimate, fftw_in_place,
				 tarr);
    rfft_mgr.real_to_complex_fft(flipped_axes,
				 fftw_estimate, fftw_in_place,
				 aarr);

    double tmpamp, tmpphase;
    int indexa, indexb;
    // positive frequencies
    for(int i=0; i<tpix_axes_one; i++){
      for(int j=0; j<tpix_axes_zero/2+1; j++){
	indexa = i*2*(tpix_axes_zero/2+1) + j;
	indexb = i*2*(tpix_axes_zero/2+1) + 2*j;
	product_amps[indexa] = 
	  sqrt(tarr[indexb]*tarr[indexb] + tarr[indexb+1]*tarr[indexb+1])*
	  sqrt(aarr[indexb]*aarr[indexb] + aarr[indexb+1]*aarr[indexb+1]);
	phase_difference[indexa] = 
	  atan2(tarr[indexb+1],tarr[indexb]) - 
	  atan2(aarr[indexb+1],aarr[indexb]);
      }
    }
    delete [] tarr;
    delete [] aarr;

    // negative frequencies
    for(int j=tpix_axes_zero/2+1; j<tpix_axes_zero; j++){
      indexa = j;
      indexb = tpix_axes_zero - j;
      product_amps[indexa] = product_amps[indexb];
      phase_difference[indexa] = -phase_difference[indexb];
    }

    for(int i=1; i<tpix_axes_one; i++){
      for(int j=tpix_axes_zero/2+1; j<tpix_axes_zero; j++){
	indexa = i*2*(tpix_axes_zero/2+1) + j;
	indexb = (tpix_axes_one-i)*2*(tpix_axes_zero/2+1) + tpix_axes_zero - j;
	product_amps[indexa] = product_amps[indexb];
	phase_difference[indexa] = -phase_difference[indexb];
      }
    }

    // Perform a search for the best fitting fractional offset.
    double faca = 2*M_PI/tpix.axes[0];
    double facb = 2*M_PI/tpix.axes[1];
    double tmp, tmpa, tmpb;
    double fa, fb, b, chisq;
    double dfa_du, dfa_dv, dfb_du, dfb_dv;
    double min_chisq = 0, min_fa, min_fb;
    double iindex, jindex;
    double last_chisq = 0;
    double frac_axis_0=0, frac_axis_1=0;
    int index; 

    // Solve for best offsets
    while(1){
      fa=fb=dfa_du=dfa_dv=dfb_dv=chisq=0;
      for(int i=0; i<tpix_axes_one; i++){
	for(int j=0; j<tpix_axes_zero; j++){
	  index = 2*i*(tpix_axes_zero/2+1)+j;
	  if(i<=tpix_axes_one/2) iindex = i;
	  else iindex = i-tpix_axes_one;
	  if(j<=tpix_axes_zero/2) jindex = j;
	  else jindex = j-tpix_axes_zero;

	  tmpa = product_amps[index]*sin(phase_difference[index] - 
					 iindex*frac_axis_1*facb -
					 jindex*frac_axis_0*faca);
	  fa += iindex*tmpa;
	  fb += jindex*tmpa;
	  tmpb = product_amps[index]*cos(phase_difference[index] - 
					 iindex*frac_axis_1*facb -
					 jindex*frac_axis_0*faca); 

	  dfa_du += iindex*iindex*tmpb;
	  dfb_dv += jindex*jindex*tmpb;
	  dfa_dv += iindex*jindex*tmpb;
	  chisq -= tmpb;
	}
      }
      
      dfb_du = dfa_dv;
      tmpa = (fa*dfb_dv - fb*dfa_dv)/(dfa_du*dfb_dv-dfa_dv*dfb_du);
      tmpb = (-fa*dfb_du + fb*dfa_du)/(dfa_du*dfb_dv-dfa_dv*dfb_du);
      frac_axis_0+=tmpb;
      frac_axis_1+=tmpa;

      if(fabs(tmpa)<1e-6 && fabs(tmpb)<1e-6)
	break;
      last_chisq = chisq;
    }
    if(pixel_array<T>::verbose_level) 
      cout << "pixel_array<T>::offset - " << endl
	   << " dfa_du " << dfa_du << " dfa_dv " << dfa_dv
	   << " dfb_du " << dfb_du << " dfb_dv " << dfb_dv << endl
	   << " fa " << setw(10) << setprecision(3) << fa 
	   << " fb " << setw(10) << setprecision(3) << fb << endl 
	   << " tmpa " << tmpa << " tmpb " << tmpb << endl
	   << " chisq " << setw(10) << setprecision(6) << chisq
	   << setprecision(3) << endl 
	   << " dra " << setw(10) << frac_axis_0 << endl
	   << " ddec " << setw(10) << frac_axis_1 << endl    
	   << " axis 0 " << setw(10) << setprecision(6)
	   << xcorr_maxpixel[0] + frac_axis_0 << endl
	   << " axis 1 " << setw(10) << xcorr_maxpixel[1] + frac_axis_1 << endl;
    


    /*
      if(fabs(frac_axis_0)>1 || fabs(frac_axis_1)>1){
      cerr << "pixel_array::offset error - value for fractional pixel shift "
      << frac_axis_0 << "\t" << frac_axis_1 
      << " has magnitude greater than one\n";
      throw(string("pixel_array::offset"));
      }
    */

    offsets.resize(2);
    offsets[0] = xcorr_maxpixel[0] + frac_axis_0;
    offsets[1] = xcorr_maxpixel[1] + frac_axis_1;

    if(pixel_array<T>::verbose_level) 
      cout << "pixel_array<T>::offset - final offsets " 
	   << offsets[0] << "\t" << offsets[1] << endl;

    // Clean up
    delete [] product_amps;
    delete [] phase_difference;
  }

  // template <> pixel_array<long> & operator+=<long, float>(
	// 						  pixel_array<long> & lhs, const pixel_array<float> & rhs);
  // template <> pixel_array<long> & operator+=<long, double>(
	// 						   pixel_array<long> & lhs, const pixel_array<double> & rhs);

  template<class T, class U>
    pixel_array<T> & operator += (pixel_array<T> & lhs,
				  const pixel_array<U> & rhs){
    if(lhs.get_axes()!=rhs.get_axes()){
      cerr << "pixel_array::operator+= error - " 
	   << "mismatched array sizes\n";
      for(uint i=0; i<lhs.get_axes().size(); i++){
	cerr << setw(8) << lhs.get_axes()[i];
      }
      cerr << endl;
      for(uint i=0; i<rhs.get_axes().size(); i++){
	cerr << setw(8) << rhs.get_axes()[i];
      }
      cerr << endl;
      throw(string("pixel_array::operator+="));
    }

    if((lhs.pixelwts!=NULL && rhs.weights_allocated()==0) ||
       (lhs.pixelwts==NULL && rhs.weights_allocated()!=0)){
      cerr << "pixel_array<T> & operator += error - "
	   << "weights not defined for both objects\n";
      throw(string("pixel_array<T> & operator +="));
    }

    int nelem = lhs.total_space();
    if(pixel_array<T>::verbose_level)
      cout << "pixel_array::operator+= - adding data\n";
    for(int i=0; i<nelem; i++){
      if(lhs.pixelwts!=NULL && lhs.pixelwts[i] == 0 && rhs.wt(i) == 0){
	lhs.pixeldata[i] = 0;	
	lhs.pixelwts[i] = 0;
      } else if(lhs.pixelwts!=NULL && lhs.pixelwts[i] == 0 && rhs.wt(i) != 0){ 
	lhs.pixeldata[i] = rhs.pixeldata[i];
	lhs.pixelwts[i] = rhs.pixelwts[i];
      } else if(lhs.pixelwts!=NULL){ 
	lhs.pixeldata[i] += rhs.pixeldata[i];
	lhs.pixelwts[i] += rhs.pixelwts[i];
      } else {
	lhs.pixeldata[i] += rhs.pixeldata[i];
      }
    }
    return(lhs);
  }

  // template <> pixel_array<long> & operator-=<long, float>(
	// 						  pixel_array<long> & lhs, const pixel_array<float> & rhs);
  // template <> pixel_array<long> & operator-=<long, double>(
	// 						   pixel_array<long> & lhs, const pixel_array<double> & rhs);

  template<class T, class U>
    pixel_array<T> & operator -= (
				  pixel_array<T> & lhs, 
				  const pixel_array<U> & rhs){
    if(lhs.get_axes()!=rhs.get_axes()){
      cerr << "pixel_array::operator-= error - " 
	   << "mismatched array sizes\n";
      for(uint i=0; i<lhs.get_axes().size(); i++){
	cerr << setw(8) << lhs.get_axes()[i];
      }
      cerr << endl;
      for(uint i=0; i<rhs.get_axes().size(); i++){
	cerr << setw(8) << rhs.get_axes()[i];
      }
      cerr << endl;
      throw(string("pixel_array::operator-="));
    }

    if((lhs.pixelwts!=NULL && rhs.weights_allocated()==0) ||
       (lhs.pixelwts==NULL && rhs.weights_allocated()!=0)){
      cerr << "pixel_array<T> & operator-= error - "
	   << "weights not defined for both objects\n";
      throw(string("pixel_array<T> & operator-="));
    }

    int nelem = lhs.total_space();
    if(pixel_array<T>::verbose_level)
      cout << "pixel_array::operator-= - subtracting data\n";
    for(int i=0; i<nelem; i++){
      if(lhs.pixelwts!=NULL && lhs.pixelwts[i] == 0 && rhs.pixelwts[i] == 0){
	lhs.pixeldata[i] = 0;	
	lhs.pixelwts[i] = 0;
      } else if(lhs.pixelwts!=NULL && lhs.pixelwts[i] == 0 && rhs.pixelwts[i] != 0){ 
	lhs.pixeldata[i] = -rhs.pixeldata[i];
	lhs.pixelwts[i] = rhs.pixelwts[i];
      } else if(lhs.pixelwts!=NULL){ 
	lhs.pixeldata[i] -= rhs.pixeldata[i];
	lhs.pixelwts[i] += rhs.pixelwts[i];
      } else {
	lhs.pixeldata[i] -= rhs.pixeldata[i];
      }
    }
    return(lhs);
  }

  // template <> pixel_array<long> & operator*=<long, float>(
	// 						  pixel_array<long> & lhs, const pixel_array<float> & rhs);
  // template <> pixel_array<long> & operator*=<long, double>(
	// 						   pixel_array<long> & lhs, const pixel_array<double> & rhs);

  template<class T, class U>
    pixel_array<T> & operator *= (
				  pixel_array<T> & lhs, 
				  const pixel_array<U> & rhs){
    if(lhs.get_axes()!=rhs.get_axes()){
      cerr << "pixel_array::operator*= error - " 
	   << "mismatched array sizes\n";
      for(uint i=0; i<lhs.get_axes().size(); i++){
	cerr << setw(8) << lhs.get_axes()[i];
      }
      cerr << endl;
      for(uint i=0; i<rhs.get_axes().size(); i++){
	cerr << setw(8) << rhs.get_axes()[i];
      }
      cerr << endl;
      throw(string("pixel_array::operator*="));
    }

    if((lhs.pixelwts!=NULL && rhs.weights_allocated()==0) ||
       (lhs.pixelwts==NULL && rhs.weights_allocated()!=0)){
      cerr << "pixel_array<T> & operator*= error - "
	   << "weights not defined for both objects\n";
      throw(string("pixel_array<T> & operator*="));
    }

    int nelem = lhs.total_space();
    if(pixel_array<T>::verbose_level)
      cout << "pixel_array::operator*= - multiplying data\n";
    for(int i=0; i<nelem; i++){
      if(lhs.pixelwts!=NULL && lhs.pixelwts[i] == 0 && rhs.pixelwts[i] == 0){
	lhs.pixeldata[i] = 0;	
	lhs.pixelwts[i] = 0;
      } else if(lhs.pixelwts!=NULL && lhs.pixelwts[i] == 0
		&& rhs.pixelwts[i] != 0){ 
	lhs.pixeldata[i] = rhs.pixeldata[i];
	lhs.pixelwts[i] = rhs.pixelwts[i];
      } else if(lhs.pixelwts!=NULL){ 
	lhs.pixeldata[i] *= rhs.pixeldata[i];
	lhs.pixelwts[i] = sqrt(1/((rhs.pixeldata[i]*rhs.pixeldata[i])/
				  (lhs.pixelwts[i]*lhs.pixelwts[i])
				  + (lhs.pixeldata[i]*lhs.pixeldata[i])/
				  (rhs.pixelwts[i]*rhs.pixelwts[i])));
      } else {
	lhs.pixeldata[i] *= rhs.pixeldata[i];
      }
    }
    return(lhs);
  }

  // template <> pixel_array<long> & operator/=<long, float>(
	// 						  pixel_array<long> & lhs, const pixel_array<float> & rhs);
  // template <> pixel_array<long> & operator/=<long, double>(
	// 						   pixel_array<long> & lhs, const pixel_array<double> & rhs);

  template<class T, class U>
    pixel_array<T> & operator /= (
				  pixel_array<T> & lhs, 
				  const pixel_array<U> & rhs){
    if(lhs.get_axes()!=rhs.get_axes()){
      cerr << "pixel_array::operator/= error - " 
	   << "mismatched array sizes\n";
      for(uint i=0; i<lhs.get_axes().size(); i++){
	cerr << setw(8) << lhs.get_axes()[i];
      }
      cerr << endl;
      for(uint i=0; i<rhs.get_axes().size(); i++){
	cerr << setw(8) << rhs.get_axes()[i];
      }
      cerr << endl;
      throw(string("pixel_array::operator/="));
    }

    if((lhs.pixelwts!=NULL && rhs.weights_allocated()==0) ||
       (lhs.pixelwts==NULL && rhs.weights_allocated()!=0)){
      cerr << "pixel_array<T> & operator/= error - "
	   << "weights not defined for both objects\n";
      throw(string("pixel_array<T> & operator/="));
    }

    int nelem = lhs.total_space();
    if(pixel_array<T>::verbose_level)
      cout << "pixel_array::operator/= - dividing data\n";
    for(int i=0; i<nelem; i++){
      if(rhs.pixeldata[i]==0){
	lhs.pixeldata[i] = 0;	
	if(lhs.pixelwts!=NULL)
	  lhs.pixelwts[i] = 0;
      } else if(lhs.pixelwts!=NULL){ 
	lhs.pixeldata[i] /= rhs.pixeldata[i];
	lhs.pixelwts[i] = 
	  sqrt(1/(1/(lhs.pixelwts[i]*lhs.pixelwts[i]*
		     rhs.pixeldata[i]*rhs.pixeldata[i]) +
		  (lhs.pixeldata[i]*lhs.pixeldata[i])/
		  (rhs.pixeldata[i]*rhs.pixeldata[i]*
		   rhs.pixeldata[i]*rhs.pixeldata[i]
		   *lhs.pixelwts[i]*lhs.pixelwts[i])));
      } else lhs.pixeldata[i] /= rhs.pixeldata[i];
    }
    return(lhs);
  }

  template <> pixel_array<long> & pixel_array<long>::operator += (
								  const double & fac);

  template<class T>
    pixel_array<T> & pixel_array<T>::operator += (const double & fac){
    if(fac==0) return(*this);
    int nelem = total_space();
    for(int i=0; i<nelem; i++)
      this->pixeldata[i] += fac;
    return(*this);
  }

  template <> pixel_array<long> & pixel_array<long>::operator -= (
								  const double & fac);

  template<class T>
    pixel_array<T> & pixel_array<T>::operator -= (const double & fac){
    if(fac==0) return(*this);
    int nelem = total_space();
    for(int i=0; i<nelem; i++)
      this->pixeldata[i] -= fac;
    return(*this);
  }

  template <> pixel_array<long> & pixel_array<long>::operator *= (
								  const double & fac);

  template<class T>
    pixel_array<T> & pixel_array<T>::operator *= (const double & fac){
    if(fac==1) return(*this);
    int nelem = total_space();
    for(int i=0; i<nelem; i++){
      this->pixeldata[i] *= fac;
    }   
    if(pixelwts!=NULL) 
      for(int i=0; i<nelem; i++)
	this->pixelwts[i] *= fac;
    return(*this);
  }

  template <> pixel_array<long> & pixel_array<long>::operator /= (
								  const double & fac);

  template<class T>
    pixel_array<T> & pixel_array<T>::operator /= (const double & fac){
    if(fac==1) return(*this);
    if(fac==0){
      cerr << "pixel_array::operator /= division by zero error\n";
      throw(string("pixel_array::operator /="));
    }
    int nelem = total_space();
    for(int i=0; i<nelem; i++)
      this->pixeldata[i] /= fac;
    if(pixelwts!=NULL)
      for(int i=0; i<nelem; i++)
	this->pixelwts[i] /= fac;

    return(*this);
  }

  template<class T> 
    pixel_array<T> operator + (const pixel_array<T> &p1,
			       const pixel_array<T> &p2){
    if(p1.get_axes()!=p2.get_axes()){
      cerr << "pixel_array::operator+ error - " 
	   << "mismatched array sizes\n";
      for(uint i=0; i<p1.get_axes().size(); i++){
	cerr << setw(8) << p1.get_axes()[i];
      }
      cerr << endl;
      for(uint i=0; i<p2.get_axes().size(); i++){
	cerr << setw(8) << p2.get_axes()[i];
      }
      cerr << endl;
      throw(string("pixel_array::operator+"));
    }
    pixel_array<T> pixarr(p1);
    pixarr += p2;
    return(pixarr);
  }

  template <class T> 
    pixel_array<T> operator - (const pixel_array<T> &p1,
			       const pixel_array<T> &p2){
    if(p1.get_axes()!=p2.get_axes()){
      cerr << "pixel_array::operator- error - " 
	   << "mismatched array sizes\n";
      for(uint i=0; i<p1.get_axes().size(); i++){
	cerr << setw(8) << p1.get_axes()[i];
      }
      cerr << endl;
      for(uint i=0; i<p2.get_axes().size(); i++){
	cerr << setw(8) << p2.get_axes()[i];
      }
      cerr << endl;
      throw(string("pixel_array::operator-"));
    }
    pixel_array<T> pixarr(p1);
    pixarr -= p2;
    return(pixarr);
  }

  template <class T> 
    pixel_array<T> operator * (const pixel_array<T> &p1,
			       const pixel_array<T> &p2){
    if(p1.get_axes()!=p2.get_axes()){
      cerr << "pixel_array::operator* error - " 
	   << "mismatched array sizes\n";
      for(uint i=0; i<p1.get_axes().size(); i++){
	cerr << setw(8) << p1.get_axes()[i];
      }
      cerr << endl;
      for(uint i=0; i<p2.get_axes().size(); i++){
	cerr << setw(8) << p2.get_axes()[i];
      }
      cerr << endl;
      throw(string("pixel_array::operator/"));
    }
    pixel_array<T> pixarr(p1);
    pixarr *= p2;
    return(pixarr);
  }

  template <class T> 
    pixel_array<T> operator / (const pixel_array<T> &p1,
			       const pixel_array<T> &p2){
    if(p1.get_axes()!=p2.get_axes()){
      cerr << "pixel_array::operator/ error - " 
	   << "mismatched array sizes\n";
      for(uint i=0; i<p1.get_axes().size(); i++){
	cerr << setw(8) << p1.get_axes()[i];
      }
      cerr << endl;
      for(uint i=0; i<p2.get_axes().size(); i++){
	cerr << setw(8) << p2.get_axes()[i];
      }
      cerr << endl;
      throw(string("pixel_array::operator/"));
    }
    pixel_array<T> pixarr(p1);
    pixarr /= p2;
    return(pixarr);
  }

  template <class T> 
    pixel_array<T> operator + (const pixel_array<T> &p1, double & fac){
    pixel_array<T> pixarr(p1);
    pixarr += fac;
    return(pixarr);
  }

  template <class T> 
    pixel_array<T> operator - (const pixel_array<T> &p1, double & fac){
    pixel_array<T> pixarr(p1);
    pixarr -= fac;
    return(pixarr);
  }

  template <class T> 
    pixel_array<T> operator * (const pixel_array<T> &p1, double & fac){
    pixel_array<T> pixarr(p1);
    pixarr *= fac;
    return(pixarr);
  }

  template <class T> 
    pixel_array<T> operator / (const pixel_array<T> &p1, double & fac){
    pixel_array<T> pixarr(p1);
    pixarr /= fac;
    return(pixarr);
  }

  ///////////////////////////////////////////
  ///  Operator !=  for pixel_amp_array
  template <class T> 
    bool operator != (const pixel_array<T> &p1, const pixel_array<T> &p2){
    return(!(p1==p2));
  }
}

#endif 
