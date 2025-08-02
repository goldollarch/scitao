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

#ifndef PIXEL_PHASE_ARRAY_H
#define PIXEL_PHASE_ARRAY_H

#include <iostream>
#include "pixel_array.h"
#include "zernike.h"

namespace Arroyo {

  using std::vector;

  ///
  /// A class to hold and manipulate square phase images
  ///

  template <class T>
    class pixel_phase_array :
    public pixel_array<T> {

    public:

    ///////////////////////////////////////////
    ///  Null constructor  
    pixel_phase_array(){};

    ///////////////////////////////////////////
    ///  Copy constructor
    pixel_phase_array(const pixel_phase_array<T> & pixpharr){
      operator= (pixpharr);
    };

    ///////////////////////////////////////////
    ///  Construct from iofits object
    pixel_phase_array(const iofits & iof) : pixel_array<T>(iof) {};

    ///////////////////////////////////////////
    ///  Construct from the bits.
    /// 
    ///  The data array should hold the phases
    ///  measured in radians.
    pixel_phase_array(const vector<long> & in_axes, 
		      const T * data = NULL, 
		      const float * wts = NULL) :
      pixel_array<T>(in_axes, data, wts){};

    ///////////////////////////////////////////
    ///  Construct an instance with pixel limits given by pixel_limits
    ///  These limits must be contained by pixarr's arrays
    template<class U>
      pixel_phase_array(const pixel_phase_array<U> & pixarr,
			const std::vector<long> & pixel_limits) 
        : pixel_array<T>(pixarr, pixel_limits) {}

    ///////////////////////////////////////////
    ///  Destructor
    ~pixel_phase_array(){};  

    ///////////////////////////////////////////
    ///  Operator =
    pixel_phase_array & operator = (const pixel_phase_array<T> & pixpharr){
      if(this == &pixpharr)
	return(*this);
      pixel_array<T>::operator=(pixpharr);
      return(*this);
    };

    ///////////////////////////////////////////
    ///  Wrap phase array into [-M_PI,M_PI)
    void wrap(){
      int nelem = this->total_space();
      double twopi = 2*M_PI;
      for(int i=0; i<nelem; i++){
	this->pixeldata[i] = fmod(this->pixeldata[i],twopi);
	if(this->pixeldata[i]<-M_PI)
	  this->pixeldata[i] += twopi;
	if(this->pixeldata[i]>=M_PI)
	  this->pixeldata[i] -= twopi;
      }
    };

    ///////////////////////////////////////////
    ///  Decimate resolution of pixel phase array
    void decimate(int nadd){
      if(nadd<0 || nadd>this->axes[0] || nadd>this->axes[1]){
	cerr << "pixel_phase_array::decimate - error decimating by a factor of " << nadd << endl;
	throw(string("pixel_phase_array::decimate"));
      }
      if(this->axes.size()!=2){
	cerr << "pixel_phase_array::decimate - cannot decimate pixel_array with number of dimensions " 
	     << this->axes.size() << endl;
	throw(string("pixel_phase_array::decimate"));
      }

      if(nadd==0 || nadd==1) return;

      vector<long> newaxes(2);
      for(int i=0; i<newaxes.size(); i++) 
	newaxes[i] = this->axes[i]/nadd;

      float * olddata = this->pixeldata;
      this->pixeldata = new T[newaxes[0]*newaxes[1]];

      float sumcos, sumsin;
      for(int i=0; i<newaxes[1]; i++){
	for(int j=0; j<newaxes[0]; j++){
	  sumcos = 0;
	  sumsin = 0;
	  for(int k=0; k<nadd; k++){
	    for(int l=0; l<nadd; l++){
	      sumcos += cos(olddata[(i*nadd+k)*this->axes[0]+j*nadd+l]);
	      sumsin += sin(olddata[(i*nadd+k)*this->axes[0]+j*nadd+l]);
	    }
	  }
	  this->pixeldata[i*newaxes[0]+j] = atan2(sumsin,sumcos);
	}
      }

      this->axes = newaxes;
      delete [] olddata;
    };

    ///////////////////////////////////////////
    ///  Friend operator ==  for pixel_phase_array
    friend int operator ==(const pixel_phase_array<T> &p1, const pixel_phase_array<T> &p2){
      if(p1!=p2) return(0);
      for(int i=0; i<p1.axes.size(); i++)
	for(int j=0; j<p1.axes[i]; j++)
	  if(p1.pixeldata[i]!=p2.pixeldata[i]) return(0);
    
      if((p1.pixelwts!=NULL && p2.pixelwts==NULL) ||
	 (p1.pixelwts==NULL && p2.pixelwts!=NULL)){
	cerr << "pixel_phase_array & operator += error - weights not defined for both objects\n";
	throw(string("pixel_phase_array & operator +="));
      }
    
      if(p1.pixelwts!=NULL)
	for(int i=0; i<p1.axes.size(); i++)
	  for(int j=0; j<p1.axes[i]; j++)
	    if(p1.pixelwts[i]!=p2.pixelwts[i]) return(0);
    
      return(1);
    }; 
  }; 

  template<class T>
  pixel_phase_array<T> operator + (const pixel_phase_array<T> &p1, const pixel_phase_array<T> &p2){
    if(p1.get_axes()!=p2.get_axes()){
      cerr << "pixel_phase_array::operator+= error - " 
	   << "mismatched array sizes:\n";
      throw(string("pixel_array::operator/="));
    }
    
    pixel_phase_array<T> pixpharr(p1);
    pixpharr += p2;
    return(pixpharr);
  }

  template<class T>
  pixel_phase_array<T> operator - (const pixel_phase_array<T> &p1, const pixel_phase_array<T> &p2){
    if(p1.get_axes()!=p2.get_axes()){
      cerr << "pixel_phase_array::operator-= error - " 
	   << "mismatched array sizes:\n";
      throw(string("pixel_array::operator/="));
    }

    pixel_phase_array<T> pixpharr(p1);
    pixpharr -= p2;
    return(pixpharr);
  }

  template<class T>
  pixel_phase_array<T> operator * (const pixel_phase_array<T> &p1, const pixel_phase_array<T> &p2){
    cerr << "operator* for pixel_phase_arrays has been intentionally disabled\n";
    throw(string("operator*"));
  }

  template<class T>
  pixel_phase_array<T> operator / (const pixel_phase_array<T> &p1, const pixel_phase_array<T> &p2){
    cerr << "operator/ for pixel_phase_arrays has been intentionally disabled\n";
    throw(string("operator/"));
  }

  template<class T>
  pixel_phase_array<T> operator + (const pixel_phase_array<T> &p1, double & fac){
    pixel_phase_array<T> pixpharr(p1);
    pixpharr += fac;
    return(pixpharr);
  }

  template<class T>
  pixel_phase_array<T> operator - (const pixel_phase_array<T> &p1, double & fac){
    pixel_phase_array<T> pixpharr(p1);
    pixpharr -= fac;
    return(pixpharr);
  }

  template<class T>
  pixel_phase_array<T> operator * (const pixel_phase_array<T> &p1, double & fac){
    pixel_phase_array<T> pixpharr(p1);
    pixpharr *= fac;
    return(pixpharr);
  }

  template<class T>
  pixel_phase_array<T> operator / (const pixel_phase_array<T> &p1, double & fac){
    pixel_phase_array<T> pixpharr(p1);
    pixpharr /= fac;
    return(pixpharr);
  }

  template<class T>
  int operator != (const pixel_phase_array<T> &p1, const pixel_phase_array<T> &p2){
    return(!(p1==p2));
  }

}

#endif
