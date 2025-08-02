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

#include <complex>
#include <iostream>
#include "computational_geometry.h"
#include "sim_utils.h"
#include "diffractive_wavefront.h"
#include "three_frame.h"

using namespace std;

namespace Arroyo {

// Detect overflow of integer operations for memory allocation
double alloc_size::max_size_t = rint(pow(2.0, (double)(((sizeof(size_t)*8)) - 1)));
ostream& operator<<(ostream& s, const alloc_size& as) {
  return s << fixed << setprecision(0) << as.dsz;
}

void fresnel_integral(double x, double & cf, double & sf) {
  
  if(x==0) {
    cf = sf = 0; 
    return;
  }

  double crossover_computation_limit = 1.5;
    
  double precision_limit=1e-15;
  double test;
  
  if(fabs(x)<crossover_computation_limit){
    double sum = 0;
    int sign = 1;
      double coeff = x*x*M_PI_2;
      int n=3, count=1;
      double term = fabs(x);
	cf = fabs(x);
	sf = 0;
      do {
      term *= coeff/count;
	sum += sign*term/n;
	test=fabs(sum)*precision_limit;
	if(count%2==1){
	  sf = sum;
	  sum = cf;
	  sign*=-1;
	} else {
	cf = sum;
	  sum = sf;
	  }
	count++;
	n += 2;
      }  while(test<term);
    } else {
      double fac = M_PI*x*x;
      complex<double> one(1,0);
      complex<double> four(4,0);
      complex<double> tiny(1e-30,0);

      complex<double> a;
      complex<double> b;
      complex<double> f(tiny);
      complex<double> C(f);
      complex<double> D(0,0);
      complex<double> val, lastval;

      int n = -3;
      int count =0;

      do {      
	lastval = val;
	count ++;

	// here a and b define the continued fraction 
	// That is:
	//  
	//                 a_1
	//  c.f. = b_0 + ---------------------------
	//                        a_2
	//                 b_1 + -------------------
	//                               a_3
	//                        b_2 + ------------
	//                               b_4 + ...
	// 
	//  For this case we have b_0 = 0,
	//  b_n = 2z^{2} + 1 + (n-1)*4 and
	//  a_1 = 1, a_n = -n(n+1)
	n += 2;
	if(n==-1) a = complex<double>(1,0);
	else a = complex<double>(-n*(n+1),0);

	if(n==-1) b = complex<double>(1.0,-fac);
	else b = b + four;

	C = b + a / C;
	if(abs(C)==0) C = tiny;

	D = b + a * D;
	if(abs(D)==0) D = tiny;
	D = one / D;

	f = f * C * D;

	val = f;

	// with this multiply we have now gotten e^{z^{2}} erfc z
	val *= complex<double>(fabs(x),-fabs(x));
    
	// remove the e^{z^{2}}
	val *= complex<double>(cos(.5*fac),sin(.5*fac));

	// finally, we convert from erfc z to erf z
	// and remove the prefactor to get cf + i sf
	val = one - val;
	val *= complex<double>(.5,.5);

	cf = val.real();
	sf = val.imag();      
      
      } while(fabs(abs(val)-abs(lastval)) > precision_limit);
    }
  
  if(x<0) {
    cf = -cf;
    sf = -sf;
  }

}

  // physical dimensions in meters
  // wavelength in meters
  // distance in meters
  // pixel scale in meters per pixel
  diffractive_wavefront<double> propagation_from_rectangular_aperture(vector<double> aperture_physical_dimensions,
								      vector<long> final_axes,
								      double distance, 
								      double wavelength,
								      double pixel_scale) {
  

    if(wavelength<=0){
      cerr << "propagation_from_rectangular_aperture error - "
	   << "pixel scale " << wavelength << " isn't positive\n";
      throw(string("propagation_from_rectangular_aperture"));
    }

    if(pixel_scale<=0){
      cerr << "propagation_from_rectangular_aperture error - "
	   << "pixel scale " << pixel_scale << " isn't positive\n";
      throw(string("propagation_from_rectangular_aperture"));
    }

    if(aperture_physical_dimensions.size()!=2 ||
       aperture_physical_dimensions[0] <= 0 ||
       aperture_physical_dimensions[1] <= 0){

      cerr << "propagation_from_rectangular_aperture error - "
	   << "invalid aperture dimensions:\n"
	   << "\t" << aperture_physical_dimensions[0] 
	   << " x " << aperture_physical_dimensions[1] << endl;
      throw(string("propagation_from_rectangular_aperture"));
    }

    if(final_axes.size()!=2 ||
       final_axes[0] <= 0 ||
       final_axes[1] <= 0){

      cerr << "propagation_from_rectangular_aperture error - "
	   << "invalid final axes:\n"
	   << "\t" << final_axes[0] 
	   << " x " << final_axes[1] << endl;
      throw(string("propagation_from_rectangular_aperture"));
    }

    int nelem=final_axes[0]*final_axes[1];

    double * data;
    double constant_phase = 2*M_PI*fmod(distance/wavelength,1.0) - M_PI_2;
    double fac;
    double a_1, a_2, b_1, b_2;
    double cos_fresnel, sin_fresnel;
    double real_1, real_2, imag_1, imag_2;
    double real, imag;
    double twopi = 2*M_PI;

    double x_halfpix=0, y_halfpix=0;
    int x_extrapix=1, y_extrapix=1;
    if(final_axes[1]%2==0){
      x_halfpix = .5;
      x_extrapix = 0;
    }
    if(final_axes[0]%2==0){
      y_halfpix = .5;
      y_extrapix = 0;
    }

    try{data = new double[2*nelem];}
    catch(...){
      cerr << "propagation_from_rectangular_aperture error - "
	   << "unable to allocate memory for wavefront data\n";
      throw;
    }

    if(wavelength == 0 || distance == 0){
      // Here what you should do is fill in axes according to whether
      // they are less than or greater than the aperture dimensions
      cerr << "propagation_from_rectangular_aperture error - wavelength or distance zero\n";
      throw(string("propagation_from_rectangular_aperture"));
    } else {
  
      fac = sqrt(2/wavelength/distance);
    
      // Here we only calculate one quadrant, since the others
      // are the same by symmetry
      for(int i=0; i<final_axes[1]/2+x_extrapix; i++){
	for(int j=0; j<final_axes[0]/2+y_extrapix; j++){
	
	  a_1 = -(aperture_physical_dimensions[1]/2.0+(i+x_halfpix)*pixel_scale)*fac;
	  a_2 = (aperture_physical_dimensions[1]/2.0-(i+x_halfpix)*pixel_scale)*fac;
	  b_1 = -(aperture_physical_dimensions[0]/2.0+(j+y_halfpix)*pixel_scale)*fac;
	  b_2 = (aperture_physical_dimensions[0]/2.0-(j+y_halfpix)*pixel_scale)*fac;
	
	  fresnel_integral(a_1, cos_fresnel, sin_fresnel);
	  real_1 = -cos_fresnel;  imag_1 = -sin_fresnel;
	
	  fresnel_integral(a_2, cos_fresnel, sin_fresnel);
	  real_1 += cos_fresnel;  imag_1 += sin_fresnel;

	  fresnel_integral(b_1, cos_fresnel, sin_fresnel);
	  real_2 = -cos_fresnel;  imag_2 = -sin_fresnel;
	
	  fresnel_integral(b_2, cos_fresnel, sin_fresnel);
	  real_2 += cos_fresnel;  imag_2 += sin_fresnel;
	
	  real = .5*(real_1*real_2 - imag_1*imag_2);
	  imag = .5*(real_1*imag_2 + real_2*imag_1);
	
	  data[(i+final_axes[1]/2)*final_axes[0]+j+final_axes[0]/2] = sqrt(real*real + imag*imag);
	  data[(i+final_axes[1]/2)*final_axes[0]+j+final_axes[0]/2+nelem] = atan2(imag,real) + constant_phase;
	  while(data[(i+final_axes[1]/2)*final_axes[0]+j+final_axes[0]/2+nelem] < -M_PI)
	    data[(i+final_axes[1]/2)*final_axes[0]+j+final_axes[0]/2+nelem] += twopi;
	  while(data[(i+final_axes[1]/2)*final_axes[0]+j+final_axes[0]/2+nelem] > M_PI)
	    data[(i+final_axes[1]/2)*final_axes[0]+j+final_axes[0]/2+nelem] -= twopi;
	
	  data[(final_axes[1]/2-i-1+x_extrapix)*final_axes[0]+j+final_axes[0]/2] = 
	    data[(final_axes[1]/2-i-1+x_extrapix)*final_axes[0]+final_axes[0]/2-j-1+y_extrapix] = 
	    data[(i+final_axes[1]/2)*final_axes[0]+final_axes[0]/2-j-1+y_extrapix] = 
	    data[(i+final_axes[1]/2)*final_axes[0]+j+final_axes[0]/2];
	  data[(final_axes[1]/2-i-1+x_extrapix)*final_axes[0]+j+final_axes[0]/2+nelem] = 
	    data[(final_axes[1]/2-i-1+x_extrapix)*final_axes[0]+final_axes[0]/2-j-1+nelem+y_extrapix] = 
	    data[(i+final_axes[1]/2)*final_axes[0]+final_axes[0]/2-j-1+nelem+y_extrapix] = 
	    data[(i+final_axes[1]/2)*final_axes[0]+j+final_axes[0]/2+nelem];
	}
      } 
    }

    three_frame tf;
    three_translation ttrans(distance*tf.z());
    ttrans.transform(tf);
    diffractive_wavefront_header<double> wfh(final_axes, tf, wavelength, pixel_scale);
    diffractive_wavefront<double> wf(wfh);
  
    wf.install(pixel_amp_array<double>(final_axes, data));
    wf.install(pixel_phase_array<double>(final_axes, &(data[nelem])));

    wf.set_timestamp(distance/299792458);

    delete [] data;

    return(wf);
  }
}

  
