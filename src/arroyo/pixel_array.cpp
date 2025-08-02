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

#include "pixel_array.h"

namespace Arroyo {
  template<>
    void pixel_array<long>::shift_by_fft(double dx, double dy){
    cerr << "pixel_array<long>::shift_by_fft error - cannot shift "
	 << "array of longs using fft\n";
    throw(string("pixel_array::shift_by_fft"));
  }

  template<>
    double pixel_array<long>::normalize_by_wts(){
    cerr << "pixel_array<long>::normalize_by_wts error - "
	 << "cannot normalize long instantiation of pixel_array by weights\n";
    throw(string("pixel_array<long>::normalize_by_wts"));
  }

  template<>
  pixel_array<long> & pixel_array<long>::operator += (const double & fac){
    if(fmod(fac, 1.0)!=0){
      cerr << "pixel_array<long>::operator+= error - cannot add a double to a long pixel array "
	   << " to get a long pixel array\n";
      throw(string("pixel_array<long>::operator+="));
    }
    long lfac = (long)(floor(fac));
    if(lfac==0) return(*this);
    int nelem = total_space();
    for(int i=0; i<nelem; i++)
      this->pixeldata[i] += lfac;
    return(*this);
  }

  template<>
  pixel_array<long> & pixel_array<long>::operator -= (const double & fac){
    if(fmod(fac, 1.0)!=0){
      cerr << "pixel_array<long>::operator-= error - cannot subtract a double from a long pixel array "
	   << " to get a long pixel array\n";
      throw(string("pixel_array<long>::operator-="));
    }
    long lfac = (long)(floor(fac));
    if(lfac==0) return(*this);
    int nelem = total_space();
    for(int i=0; i<nelem; i++)
      this->pixeldata[i] -= lfac;
    return(*this);
  }

  template<>
  pixel_array<long> & pixel_array<long>::operator *= (const double & fac){
    if(fmod(fac, 1.0)!=0){
      cerr << "pixel_array<long>::operator*= error - cannot multiply a long pixel array by a double "
	   << " to get a long pixel array\n";
      throw(string("pixel_array<long>::operator*="));
    }
    long lfac = (long)(floor(fac));
    if(lfac==1) return(*this);
    int nelem = total_space();
    for(int i=0; i<nelem; i++){
      this->pixeldata[i] *= lfac;
    }   
    if(pixelwts!=NULL) 
      for(int i=0; i<nelem; i++)
	this->pixelwts[i] *= lfac;
    return(*this);
  }

  template<>
  pixel_array<long> & pixel_array<long>::operator /= (const double & fac){
    cerr << "pixel_array<long>::operator/= error - cannot divide a long pixel array by a double "
	 << " to get a long pixel array\n";
    throw(string("pixel_array<long>::operator/="));    
  }


  template<>
  pixel_array<long> & operator+=<long, float>(pixel_array<long> & lhs, const pixel_array<float> & rhs){
    cerr << "operator+=<long, float> error - cannot add a float pixel array to a long pixel array "
	 << " to get a long pixel array\n";
    throw(string("operator+=<long, float>"));
  }
  template<>
  pixel_array<long> & operator+=<long, double>(pixel_array<long> & lhs, const pixel_array<double> & rhs){
    cerr << "operator+=<long, double> error - cannot add a double pixel array to a long pixel array "
	 << " to get a long pixel array\n";
    throw(string("operator+=<long, double>"));
  }

  template<>
  pixel_array<long> & operator-=<long, float>(pixel_array<long> & lhs, const pixel_array<float> & rhs){
    cerr << "operator-=<long, float> error - cannot subtract a float pixel array from a long pixel array "
	 << " to get a long pixel array\n";
    throw(string("operator-=<long, float>"));
  }
  template<>
  pixel_array<long> & operator-=<long, double>(pixel_array<long> & lhs, const pixel_array<double> & rhs){
    cerr << "operator-=<long, double> error - cannot subtract a double pixel array from a long pixel array "
	 << " to get a long pixel array\n";
    throw(string("operator-=<long, double>"));
  }

  template<>
  pixel_array<long> & operator*=<long, float>(pixel_array<long> & lhs, const pixel_array<float> & rhs){
    cerr << "operator*=<long, float> error - cannot multiply a float pixel array by a long pixel array "
	 << " to get a long pixel array\n";
    throw(string("operator*=<long, float>"));
  }
  template<>
  pixel_array<long> & operator*=<long, double>(pixel_array<long> & lhs, const pixel_array<double> & rhs){
    cerr << "operator*=<long, double> error - cannot multiply a double pixel array by a long pixel array "
	 << " to get a long pixel array\n";
    throw(string("operator*=<long, double>"));
  }

  template<>
  pixel_array<long> & operator/=<long, float>(pixel_array<long> & lhs, const pixel_array<float> & rhs){
    cerr << "operator/=<long, float> error - cannot divide a float pixel array by a long pixel array "
	 << " to get a long pixel array\n";
    throw(string("operator/=<long, float>"));
  }
  template<>
  pixel_array<long> & operator/=<long, double>(pixel_array<long> & lhs, const pixel_array<double> & rhs){
    cerr << "operator/=<long, double> error - cannot divide a double pixel array by a long pixel array "
	 << " to get a long pixel array\n";
    throw(string("operator/=<long, double>"));
  }

}
