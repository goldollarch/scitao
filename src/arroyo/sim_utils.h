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

#ifndef SIM_UTILS_H
#define SIM_UTILS_H

#include <cmath>
#include <iostream>
#include <cassert>
#include "three_point.h"
#include "three_vector.h"
#include "three_frame.h"
#include "fits_header_data.h"

namespace Arroyo {

  using std::vector;
  using std::string;
  using std::ostream;
  using std::bad_alloc;

  class alloc_size {
  public:
    alloc_size(size_t dim1, size_t dim2 = 1, size_t dim3 = 1) {
      dsz = (double)dim1 * (double)dim2 * (double)dim3;
      isz = dim1 * dim2 * dim3;
    }
    operator size_t() {
      if (dsz >= max_size_t) {
        throw (new bad_alloc);
      }
      return (isz);
    }
    friend ostream& operator<<(ostream&, const alloc_size&);
  private:
    static double max_size_t;
    double dsz;
    size_t isz;
  };

  ostream& operator<<(ostream&, const alloc_size&);


  /// forward declaration
  template <typename T> class diffractive_wavefront;


  ///////////////////////////////////////////////
  ///  Returns the fresnel sin and cos integrals
  ///  for the value x
  ///  Modified from the numerical recipes in C 
  ///  routine.  Warning - the text version of the
  ///  NRC routine contains bugs.
  void fresnel_integral(double x, double & cf, double & sf);


  ///////////////////////////////////////////////
  /// Returns a wavefront constructed via the 
  /// analytic solution for wave propagation
  /// through a rectangular aperture.
  ///
  /// The initial physical dimension is that
  /// of the aperture, in meters
  /// The final physical dimension specifies
  /// the dimensions of the final wavefront, in 
  /// meters
  /// distance is measured from the aperture,in 
  /// meters
  /// wavelength is in meters
  /// pixel scale of the final wavefront is in 
  /// meters/pixel
  diffractive_wavefront<double> 
    propagation_from_rectangular_aperture(vector<double> aperture_physical_dimensions,
					  vector<long> final_axes,
					  double dist, 
					  double wavelength,
					  double pixel_scale);


  ///////////////////////////////////////////////
  /// This function effects the Goertzel-Reinsch
  /// recursion relation for computing the values 
  /// data*cos(angle) and data*sin(angle) summed 
  /// over dimen data elements.  
  ///
  /// The algorithm is described on page 84 of
  /// Stoer & Burlisch, "Introduction to Numerical 
  /// Analysis" 2nd edition (Springer-Verlag 1993)
  ///
  /// Time scale:
  ///
  /// This routine requires order 3*dimen additions and dimen
  /// multiplications to compute these sums.
  ///
  /// Arguments:
  ///
  /// dimen is the dimensionality of the data
  ///
  /// data is a pointer to the beginning of the data
  ///
  /// Angle is in radians
  template<class T>
    void goertzel_reinsch_recursor(int dimen, T * data, double angle,
				   T & real, T & imag) {
    double ca = cos(angle);
    double sa = sin(angle);
    double delta_real, delta_imag;
    double value_real, value_imag;
    double lambda;
    int dc_index = dimen/2;

    if(ca>0){
      lambda = -2*(1-ca);  // -4*sin(angle/2)*sin(angle/2)

      // traverse the positive frequency elements
      delta_real = delta_imag = value_real = value_imag = 0;
      for(int i=dimen-1; i>=dc_index; i--){
	value_real += delta_real;
	value_imag += delta_imag;
	delta_real = lambda*value_real + delta_real + data[2*i];
	delta_imag = lambda*value_imag + delta_imag + data[2*i+1];
      }
      real = delta_real - .5*lambda*value_real - sa*value_imag;
      imag = delta_imag - .5*lambda*value_imag + sa*value_real;

      // traverse the negative frequency elements
      delta_real = delta_imag = value_real = value_imag = 0;
      for(int i=0; i<=dc_index; i++){
	value_real += delta_real;
	value_imag += delta_imag;
	delta_real = lambda*value_real + delta_real + data[2*i];
	delta_imag = lambda*value_imag + delta_imag + data[2*i+1];
      }
      real += delta_real - .5*lambda*value_real + sa*value_imag;
      imag += delta_imag - .5*lambda*value_imag - sa*value_real;

      // subtract out dc, which has been double counted
      real -= data[2*dc_index];
      imag -= data[2*dc_index+1];

    } else {
      lambda = 2*(1+ca);   // 4*cos(angle/2)*cos(angle/2)

      // traverse the positive frequency elements
      delta_real = delta_imag = value_real = value_imag = 0;
      for(int i=dimen-1; i>=dc_index; i--){
	value_real = delta_real - value_real;
	value_imag = delta_imag - value_imag;
	delta_real = lambda*value_real - delta_real + data[2*i];
	delta_imag = lambda*value_imag - delta_imag + data[2*i+1];
      }
      real = delta_real - .5*lambda*value_real - sa*value_imag;
      imag = delta_imag - .5*lambda*value_imag + sa*value_real;

      // traverse the negative frequency elements
      delta_real = delta_imag = value_real = value_imag = 0;
      for(int i=0; i<=dc_index; i++){
	value_real = delta_real - value_real;
	value_imag = delta_imag - value_imag;
	delta_real = lambda*value_real - delta_real + data[2*i];
	delta_imag = lambda*value_imag - delta_imag + data[2*i+1];
      }
      real += delta_real - .5*lambda*value_real + sa*value_imag;
      imag += delta_imag - .5*lambda*value_imag - sa*value_real;

      // subtract out dc, which has been double counted
      real -= data[2*dc_index];
      imag -= data[2*dc_index+1];

    }
  }

  ///////////////////////////////////////////////
  ///  The user must allocate memory for the 
  ///  final data array before calling this 
  ///  function
  template<class T>
    void goertzel_reinsch_transform(vector<long> initial_axes, vector<long> final_axes, 
				    double sampling_factor, 
				    T * initial_data, T * final_data) {

    if(initial_axes.size()!=2){
      cerr << "goertzel_reinsch_transform error - "
	   << "initial axes of dimension "
	   << initial_axes.size() 
	   << " rather than dimension 2\n"; 
      throw(string("goertzel_reinsch_transform"));
    }

    if(initial_axes[0]<=0 || initial_axes[1]<=0){
      cerr << "goertzel_reinsch_transform error - "
	   << "dimensions of initial axes "
	   << initial_axes[0] << "x" << initial_axes[1]
	   << " are invalid\n";
      throw(string("goertzel_reinsch_transform"));
    }

    if(final_axes.size()!=2){
      cerr << "goertzel_reinsch_transform error - "
	   << "final axes of dimension "
	   << final_axes.size() 
	   << " rather than dimension 2\n"; 
      throw(string("goertzel_reinsch_transform"));
    }

    if(final_axes[0]<=0 || final_axes[1]<=0){
      cerr << "goertzel_reinsch_transform error - "
	   << "dimensions of final axes "
	   << final_axes[0] << "x" << final_axes[1]
	   << " are invalid\n";
      throw(string("goertzel_reinsch_transform"));
    }

    double real, imag;
    T * tmp;
    alloc_size sz(initial_axes[1], 2);
    //try{tmp = new T[2*initial_axes[1]];}
    try{tmp = new T[sz];}
    catch(...){
      cerr << "goertzel_reinsch_transform error - "
	   << "unable to allocate memory\n";
      throw;
    }

    double tmp_real, tmp_imag;
    for(int i=0; i<final_axes[1]; i++){ 
      for(int j=0; j<initial_axes[1]; j++){
	goertzel_reinsch_recursor(initial_axes[0], 
				  &(initial_data[2*j*initial_axes[0]]),
				  2*M_PI*sampling_factor*(i-final_axes[1]/2),
				  tmp[2*j], tmp[2*j+1]);
      }

      for(int j=0; j<final_axes[0]; j++){
	goertzel_reinsch_recursor(initial_axes[1], 
				  tmp, 
				  2*M_PI*sampling_factor*(j-final_axes[0]/2),
				  final_data[2*(j*final_axes[0]+i)],
				  final_data[2*(j*final_axes[0]+i)+1]);
      } 
    }

    delete [] tmp;
  }
}

#endif
