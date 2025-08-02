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

#ifndef DIFFRACTIVE_WAVEFRONT_H
#define DIFFRACTIVE_WAVEFRONT_H

#include <string.h>
#include <complex>
#include <fstream>
#include "colormap.h"
#include "fft_manager.h"
#include "pixel_amp_array.h"
#include "pixel_phase_array.h"
#include "sim_utils.h"
#include "wavefront.h"

namespace Arroyo {

  using std::string;
  using std::vector;
  using std::ostream;
  using std::cerr;
  using std::endl;

  /* forward declarations */
  template <class T> class diffractive_wavefront;
  template <class T> bool operator == (const diffractive_wavefront<T> &p1, const diffractive_wavefront<T> &p2);
  template <class T> bool operator != (const diffractive_wavefront<T> &p1, const diffractive_wavefront<T> &p2);
  class optic;

  ///
  /// A class to represent a diffractive diffractive_wavefront.
  ///

  template <class T>
    class diffractive_wavefront : 
    virtual public wavefront,
    virtual public diffractive_wavefront_header<T>, 
    protected fft_manager<T> {

    private:

    /// The optics are going to transform the diffractive_wavefront,
    /// and any practical, efficient implementation requires access to
    /// the protected wfdata array.  So we make the base class of the optic
    /// inheritance hierarchy a friend of diffractive_wavefront, and then write 
    /// protected member functions in optic that return protected data
    /// from diffractive_wavefront.  
    friend class one_to_one_optic;

    /// The optics are going to transform the diffractive_wavefront,
    /// and any practical, efficient implementation requires access to
    /// the protected wfdata array.  So we make the base class of the optic
    /// inheritance hierarchy a friend of diffractive_wavefront, and then write 
    /// protected member functions in optic that return protected data
    /// from diffractive_wavefront.  
    friend class one_to_many_optic;

    ///////////////////////////////////////////
    ///  Multiply amplitudes by those in a pixel_amp_array
    template<class U>
      void multiply(const pixel_amp_array<U> & pixamparr);

    ///////////////////////////////////////////
    ///  Add to phases those in a pixel_phase_array
    template<class U>  
      void add(const pixel_phase_array<U> & pixpharr);

    protected:

    /// The array of complex elements
    mutable T * wfdata;  

    /// The following functions non-destructively 
    /// modify the order and representation of the 
    /// data.  These conversions help to modify
    /// the data so that it may be transformed efficiently 
    /// by different fft libraries and may be modified
    /// by optics.
    /// e.g. FFTW requires interleaved real-imag, 
    /// while Intel MKL requires non-interleaved real-imag
  
    /// A flag to indicate whether data is real-imag or amp-phase
    mutable bool real_imag;

    /// A flag to indicate whether data is interleaved
    mutable bool interleaved;

    ///////////////////////////////////////////
    ///  Convert array to real/imag storage.
    ///  This conversion is idempotent.
    void real_imag_conversion() const ;
  
    ///////////////////////////////////////////
    ///  Convert array to amp/phase storage.
    ///  This conversion is idempotent.
    void amp_phase_conversion() const ;

    ///////////////////////////////////////////
    ///  Convert array to interleaved storage
    ///
    ///  That is, when real_imag == TRUE, we have
    ///  real[0],imag[0],real[1],imag[1]... 
    ///  and when real_imag == FALSE, we have
    ///  amp[0],phase[0],amp[1],phase[1]... 
    ///
    ///  This conversion is idempotent.
    void interleaved_conversion() const;

    ///////////////////////////////////////////
    ///  Convert array to contiguous storage.
    ///
    ///  That is, when real_imag == TRUE, we have
    ///  real[0],...real[N-1],imag[0]...imag[N-1] 
    ///  and when real_imag == FALSE, we have
    ///  amp[0],...amp[N-1],phase[0]...phase[N-1] 
    ///
    ///  This conversion is idempotent.
    void non_interleaved_conversion() const;

    private:

    static const bool factory_registration;

    ///////////////////////////////////////////
    ///  Cyclically permutes the
    ///  wfdata, shifting in x and y by a number
    ///  of pixels xshift and yshift.  
    void cyclic_permutation(long xshift, long yshift) const;

    ///////////////////////////////////////////
    ///  This function handles the near field
    ///  propagation for the exact angular
    ///  spectrum propagator and the near field 
    ///  fresnel (paraxial) propagator.
    ///  Choose fresnel propagation by setting
    ///  the fresnel argument to true.  Choose
    ///  propagation using the exact angular 
    ///  spectrum by setting it to false
    ///
    ///  Arguments:
    /// 
    ///  Distance is in meters
    void near_field_propagator(double distance, bool fresnel);
  
    ///////////////////////////////////////////
    ///  This function handles all three far field
    ///  propagators: fresnel, fraunhoffer, and
    ///  one based on the Goertzel-Reinsch algorithm.
    ///
    ///  Choose fresnel propagation by setting the fresnel argument to
    ///  true.
    ///
    ///  Choose fraunhoffer propagation by setting the fresnel argument
    ///  to false
    ///
    ///  Choose Goertzel-Reinsch propagator by supplying a new pixel
    ///  scale and axes to this function The new pixel scale must be
    ///  greater than zero and the new axes must have two positive
    ///  dimensions.
    ///
    ///  Arguments:
    /// 
    ///  Distance is in meters
    ///
    ///  Optional final_pixel_scale is in meters per pixel
    /// 
    ///  Optional axes must have 2 dimensions with positive
    ///  definite values.  These values may be even or odd
    void far_field_propagator(double distance, bool fresnel,
			      double final_pixel_scale=-1,
			      vector<long> final_axes = vector<long>());

    ///////////////////////////////////////////
    ///  Update the timestamp based on a propagation
    ///  distance dist
    ///
    ///  Speed of light defined as 299,792,458 m/s
    void update_timestamp(double dist) {timestamp += dist/299792458;};

    ///////////////////////////////////////////
    ///  Perform a forward fft on the data
    ///
    ///  This function does the real_imag and
    ///  interleaved conversions before calling
    ///  fft_manager<T>::forward_fft on the 
    ///  wavefront data
    void forward_fft();

    ///////////////////////////////////////////
    ///  Perform a backward fft on the data
    ///
    ///  This function does the real_imag and
    ///  interleaved conversions before calling
    ///  fft_manager<T>::backward_fft on the
    ///  wavefront data
    void backward_fft();
  
    public:

    ///////////////////////////////////////////
    ///  Null constructor
    diffractive_wavefront();

    ///////////////////////////////////////////
    ///  Copy constructor
    diffractive_wavefront(const diffractive_wavefront<T> & dwf);

    ///////////////////////////////////////////
    ///  Construct from a file
    diffractive_wavefront(const char * filename);

    ///////////////////////////////////////////
    ///  Construct from an iofits object
    diffractive_wavefront(const iofits & iof);

    ///////////////////////////////////////////
    ///  Construct from the bits
    diffractive_wavefront(const diffractive_wavefront_header<T> & dwfh, T * data = NULL, 
			  bool rl_img = true, bool intrlvd = true);

    ///////////////////////////////////////////
    ///  Null destructor
    ~diffractive_wavefront();

    ///////////////////////////////////////////
    ///  Operator = 
    diffractive_wavefront & operator=(const diffractive_wavefront<T> & dwf);

    ///////////////////////////////////////////
    ///  Resize the array of data 
    void resize(const vector<long> & axes);

    ///////////////////////////////////////////
    ///  Perform a pixel by pixel check, zeroing the wavefront phases
    ///  if the corresponding wavefront amplitudes are zero.
    void mask();

    ///////////////////////////////////////////
    ///  Install amplitudes from a pixel_amp_array
    template<class U>
      void install(const pixel_amp_array<U> & pixamparr);

    ///////////////////////////////////////////
    ///  Install phases from a pixel_phase_array
    template<class U>
      void install(const pixel_phase_array<U> & pixpharr);

    ///////////////////////////////////////////
    ///  Extract amplitudes from the diffractive_wavefront
    pixel_amp_array<T> extract_amps() const;

    ///////////////////////////////////////////
    ///  Extract phases from the diffractive_wavefront
    pixel_phase_array<T> extract_phases() const;

    ///////////////////////////////////////////
    ///  Wrap phases in the diffractive_wavefront
    ///  into the interval [-M_PI,M_PI)
    void wrap_phases(){
      // Optimize later
      pixel_phase_array<T> pixpharr = this->extract_phases();
      pixpharr.wrap();
      this->install(pixpharr);
    };

    ///////////////////////////////////////////
    ///  Get the wigner transform of the array
    /// wigner_transform get_vigner_transform();

    ///////////////////////////////////////////
    ///  Read diffractive_wavefront from a file
    void read(const char * filename);

    ///////////////////////////////////////////
    ///  Read diffractive_wavefront from an iofits object
    void read(const iofits & iof);

    ///////////////////////////////////////////
    ///  Write diffractive_wavefront to a file
    ///  
    ///  In general, details of the fits file 
    ///  format should be invisible to the library
    ///  user.  However, for diffractive wavefronts
    ///  it is useful to know that you can display
    ///  the amplitudes using 
    ///  ds9 filename
    ///  and the phases using 
    ///  ds9 filename[1]
    ///
    ///  Amplitudes mean the amplitude
    ///  of the complex phasor - not the amplitude
    ///  squared.
    void write(const char * filename) const;

    ///////////////////////////////////////////
    ///  Write diffractive_wavefront to an iofits object
    ///  
    ///  In general, details of the fits file 
    ///  format should be invisible to the library
    ///  user.  However, for diffractive wavefronts
    ///  it is useful to know that you can display
    ///  the amplitudes using 
    ///  ds9 filename
    ///  and the phases using 
    ///  ds9 filename[1]
    void write(iofits & iof) const;

    ///////////////////////////////////////////
    ///  Write diffractive_wavefront amplitudes
    ///  to a ppm file.  The data is discretized
    ///  into 255 levels per RGB color in this 
    ///  operation, using the colormap provided.
    ///  If the amplitudes fall outside the range
    ///  [min, max], this function throws an error.
    ///
    ///  The format of a ppm file may be found at
    ///  www.dcs.ed.ac.uk/home/mxr/gfx/
    ///
    ///  This function is meant to facilitate the 
    ///  creation of mpegs using ppmtompeg.
    void write_amps_to_ppm(double min, double max, 
			   bool logscale, 
			   bool colorbar, 
			   colormap * cmap, 
			   const char * filename, 
			   long min_dimen=-1) const;

    ///////////////////////////////////////////
    ///  Write diffractive_wavefront phases
    ///  to a ppm file.  The data is discretized
    ///  into 255 levels per RGB color in this 
    ///  operation, using the colormap provided.
    ///  If the phases fall outside the range
    ///  [min, max], this function throws an error.
    ///
    ///  The format of a ppm file may be found at
    ///  www.dcs.ed.ac.uk/home/mxr/gfx/
    ///
    ///  This function is meant to facilitate the 
    ///  creation of mpegs using ppmtompeg.
    void write_phases_to_ppm(double min, double max, 
			     bool colorbar, 
			     colormap * cmap, 
			     const char * filename, 
			     long min_dimen=-1) const;

    ///////////////////////////////////////////
    ///  Print information about the diffractive_wavefront 
    void print(ostream & os, const char * prefix) const;

    ///////////////////////////////////////////
    ///  Return nth data element as a complex number.
    /// 
    ///  The data is indexed as n - i*axes[0]+j, where
    ///
    ///  0 <= i < axes[1]
    ///
    ///  0 <= j < axes[0]
    std::complex<T> data(int n) const;

    ///////////////////////////////////////////
    ///  Set nth data element to a complex number
    /// 
    ///  The data is indexed as n - i*axes[0]+j, where
    ///
    ///  0 <= i < axes[1]
    ///
    ///  0 <= j < axes[0]
    void set_data(int n, std::complex<T> & dat);

    ///////////////////////////////////////////
    ///  Set the propagation direction to lie
    ///  along the direction prop_dir.  This member
    ///  function modifies the phase of the 
    ///  wavefront to account for the change in 
    ///  direction.
    void set_propagation_direction(const three_vector & prop_dir);

    ///////////////////////////////////////////
    ///  Sets wavefront curvature in the header and removes
    ///  this curvature from the phase
    void set_wavefront_curvature(double curvature);

    ///////////////////////////////////////////
    ///  Propagate without approximation by direct
    ///  computation of the Rayleigh-Sommerfeld
    ///  diffraction formula
    ///
    ///  This propagator permits a change in the
    ///  dimensionality and pixel scale of the 
    ///  final array.  These parameters may be 
    ///  reset using the default arguments, whose
    ///  default values leave these parameters 
    ///  unchanged.
    ///
    ///  Time scale:
    ///
    ///  for NxN array, this requires N^{4} operations
    ///
    ///  Arguments:
    ///
    ///  distance is in meters
    /// 
    ///  optional final_pixel_scale must be positive definite
    ///
    ///  optional axes must have 2 dimensions with positive definite
    ///  values.  These values may be even or odd.
    void exact_propagator(double distance, double final_pixel_scale=-1, 
			  vector<long> final_axes=vector<long>());


    ///////////////////////////////////////////
    ///  Propagate the wavefront geometrically.
    ///  The effect of this propagator is to
    ///  translate the wavefront along the z axis
    ///  of its three frame by the distance provided.
    void geometric_propagator(double distance);

    ///////////////////////////////////////////
    ///  Propagate using angular spectrum
    ///
    ///  This function will later handle verging beams
    ///  using GLAD's formulation
    ///
    ///  Time scale:
    ///
    ///  Requires 2 ffts
    ///
    ///  Arguments:
    ///
    ///  Distance is in meters
    void near_field_angular_propagator(double distance);

    ///////////////////////////////////////////
    ///  Propagate using near field fresnel propagator
    ///  (the paraxial approximation for the angular spectrum)
    ///  This function is an approximate form for
    ///  the near_field_angular_propagator
    ///
    ///  This function will later handle verging beams
    ///  using GLAD's formulation
    ///
    ///  Time scale:
    ///
    ///  Requires 2 ffts
    ///
    ///  Arguments:
    ///
    ///  Distance is in meters
    void near_field_fresnel_propagator(double distance);

    ///////////////////////////////////////////
    ///  Propagate using far field frenel propagator
    ///  (the paraxial approximation)
    ///
    ///  Time scale:
    ///
    ///  Requires 1 fft
    ///
    ///  Arguments:
    ///
    ///  Distance is in meters
    void far_field_fresnel_propagator(double distance);

    ///////////////////////////////////////////
    ///  Propagate using far field fraunhoffer propagator
    ///  This function is an approximate form for
    ///  the far_field_fresnel_propagator
    ///
    ///  Time scale:
    ///
    ///  Requires 1 fft
    ///
    ///  Arguments
    ///
    ///  Distance is in meters
    void far_field_fraunhoffer_propagator(double distance);

    ///////////////////////////////////////////
    ///  Propagate using far field fresnel propagator 
    ///  and the Goertzel-Reinsch algorithm
    ///
    ///  This propagator permits a change in the
    ///  dimensionality and pixel scale of the 
    ///  final array.
    ///
    ///  Time scale: 
    ///
    ///  unknown
    ///
    ///  Arguments
    ///
    ///  Distance is in meters
    ///
    ///  Optional final_pixel_scale is in meters per pixel
    ///
    ///  Optional axes must have 2 dimensions with positive definite
    ///  values.  These values may be even or odd.
    void far_field_fresnel_goertzel_reinsch_propagator(double distance, 
						       double final_pixel_scale, 
						       vector<long> final_axes);

    ///////////////////////////////////////////
    ///  Propagate using far field fraunhoffer propagator 
    ///  and the Goertzel-Reinsch algorithm
    ///
    ///  This propagator permits a change in the
    ///  dimensionality and pixel scale of the 
    ///  final array.
    ///
    ///  Time scale: 
    ///
    ///  Unknown
    ///
    ///  Arguments
    ///
    ///  Distance is in meters
    ///
    ///  Optional final_pixel_scale is in meters per pixel
    ///
    ///  Optional axes must have 2 dimensions with positive definite
    ///  values.  These values may be even or odd.
    void far_field_fraunhoffer_goertzel_reinsch_propagator(double distance, 
							   double final_pixel_scale, 
							   vector<long> final_axes);

    ///////////////////////////////////////////
    ///  Propagate using finite difference method
    ///
    ///  Arguments:
    ///
    ///  Distance is in meters
    void finite_difference_method_propagator(double distance);

    ///////////////////////////////////////////
    ///  Rotate array by angle 
    ///
    ///  Arguments:
    ///
    ///  Angle is in radians
    void rotate(double angle);

    ///////////////////////////////////////////
    ///  Return total power in the image
    double total_power() const;

    ///////////////////////////////////////////
    ///  Pad each edge of the array by npad pixels
    ///  and initialize to the specified value
    template<class U>
      void pad_array(int npad, std::complex<U> value);

    ///////////////////////////////////////////
    /// Clip each edge of the array by nclip pixels 
    void clip_array(int nclip);
  
    ///////////////////////////////////////////
    ///  Friend declaration for operator +=  for diffractive_wavefronts
    template<class U, class V> 
      friend diffractive_wavefront<U> & operator+=(diffractive_wavefront<U> & lhs_dwf, const diffractive_wavefront<V> & rhs_dwf);

    ///////////////////////////////////////////
    ///  Friend declaration for operator -=  for diffractive_wavefronts
    template<class U, class V> 
      friend diffractive_wavefront<U> & operator-=(diffractive_wavefront<U> & lhs_dwf, const diffractive_wavefront<V> & rhs_dwf);

    ///////////////////////////////////////////
    ///  Friend declaration for operator *=  for diffractive_wavefronts
    template<class U, class V>  
      friend diffractive_wavefront<U> & operator*=(diffractive_wavefront<U> & lhs_dwf, const diffractive_wavefront<V> & rhs_dwf);

    ///////////////////////////////////////////
    ///  Friend declaration for operator /=  for diffractive_wavefronts
    template<class U, class V>
      friend diffractive_wavefront<U> & operator/=(diffractive_wavefront<U> & lhs_dwf, const diffractive_wavefront<V> & rhs_dwf);

    ///////////////////////////////////////////
    ///  Operator +=  for doubles
    template<class U>
      diffractive_wavefront<T> & operator +=(std::complex<U> c);

    ///////////////////////////////////////////
    ///  Operator -=  for doubles
    template<class U>
      diffractive_wavefront<T> & operator -=(std::complex<U> c);

    ///////////////////////////////////////////
    ///  Operator *=  for doubles
    template<class U>
      diffractive_wavefront<T> & operator *=(std::complex<U> c);

    ///////////////////////////////////////////
    ///  Operator /=  for doubles
    template<class U>
      diffractive_wavefront<T> & operator /=(std::complex<U> c);

    ///////////////////////////////////////////
    ///  Friend operator ==  for diffractive_wavefronts
    friend bool operator ==(const diffractive_wavefront<T> & dwf1, const diffractive_wavefront<T> & dwf2) {
      if(operator!=(static_cast<const diffractive_wavefront_header<T> >(dwf1), 
		    static_cast<const diffractive_wavefront_header<T> >(dwf2)))
	return(false);
      int nelem = dwf1.total_space();
      for(int i=0; i<nelem; i++)
	if(dwf1.wfdata[i]!=dwf2.wfdata[i]) return(false);
      return(true);
    };

  };


  template<class T>
    void diffractive_wavefront<T>::resize(const vector<long> & in_axes) {

    if(in_axes.size()!=2){
      cerr << "diffractive_wavefront::resize error - axes argument does not have 2 dimensions\n";
      throw(string("diffractive_wavefront::resize"));
    }
    int nelem = 1;
    for(uint i=0; i<in_axes.size(); i++) {
      if(in_axes[i]<0){
	cerr << "diffractive_wavefront::resize error - cannot resize using an axis dimension " << in_axes[i] 
	     << " less than or equal to zero\n";
	throw(string("diffractive_wavefront::resize"));
      }
      nelem *= in_axes[i];
    }

    if(nelem != this->diffractive_wavefront_header<T>::total_space()){
      if(wavefront::verbose_level>1){
	cerr << "diffractive_wavefront::resize - resizing diffractive_wavefront from "
	     << this->axes[0] << "x" << this->axes[1] << " to " 
	     << in_axes[0] << "x" << in_axes[1] 
	     << " for a total of " << nelem << " elements\n";
      }

      if(wfdata!=NULL) delete [] wfdata;
      try{wfdata = new T[2*nelem];}
      catch(...) {
	cerr << "diffractive_wavefront::resize - unable to allocate memory\n";
	throw(string("diffractive_wavefront::resize"));
      }
    }
 
    int tmp = 2*nelem;
    for(int i=0; i<tmp; i++) 
      wfdata[i] = 0;

    this->axes = in_axes;
  }

  template<class T>
    void diffractive_wavefront<T>::cyclic_permutation(long xshift, long yshift) const {
    this->interleaved_conversion(); 
    Arroyo::complex_cyclic_permutation(this->axes, xshift, yshift, wfdata);
  }

  template<class T>
    void diffractive_wavefront<T>::real_imag_conversion() const {
  
    if(real_imag) return;

    int nelem = this->diffractive_wavefront_header<T>::total_space();
    int index;
    T tmpamp, tmpphase;
    if(interleaved){
      for(int i=0; i<nelem; i++){
	index = 2*i;
	tmpamp = wfdata[index];
	tmpphase = wfdata[index+1];
	wfdata[index] = tmpamp * cos(tmpphase);
	wfdata[index+1] = tmpamp * sin(tmpphase);
      }
    } else {
      for(int i=0; i<nelem; i++){
	tmpamp = wfdata[i];
	tmpphase = wfdata[i+nelem];
	wfdata[i] = tmpamp * cos(tmpphase);
	wfdata[i+nelem] = tmpamp * sin(tmpphase);
      }
    }

    real_imag = true;
  }
 
  template<class T>
    void diffractive_wavefront<T>::amp_phase_conversion() const {

    if(!real_imag) return;

    int nelem = this->diffractive_wavefront_header<T>::total_space();
    int index;
    T tmpreal, tmpimag;
    if(interleaved){
      for(int i=0; i<nelem; i++){
	index = 2*i;
	tmpreal = wfdata[index];
	tmpimag = wfdata[index+1];
	wfdata[index] = sqrt(tmpreal*tmpreal + tmpimag*tmpimag);
	wfdata[index+1] = atan2(tmpimag,tmpreal);
      }
    } else {
      for(int i=0; i<nelem; i++){
	tmpreal = wfdata[i];
	tmpimag = wfdata[i+nelem];
	wfdata[i] = sqrt(tmpreal*tmpreal + tmpimag*tmpimag);
	wfdata[i+nelem] = atan2(tmpimag,tmpreal);
      }
    }

    real_imag = false;
  }

  template<class T>
    void diffractive_wavefront<T>::interleaved_conversion() const {

    // here we put wfdata[i] => wfdata[2*i]
    // and wfdata[i+nelem] => wfdata[2*i+1]

    if(interleaved) return;

    int nelem = diffractive_wavefront_header<T>::total_space();
    T * tmp;

    try{tmp=new T[2*nelem];}
    catch(...){
      cerr << "diffractive_wavefront::interleaved_conversion error - "
	   << "unable to allocate memory\n";
      throw(string("diffractive_wavefront::interleaved_conversion"));
    }

    memcpy(tmp, wfdata, sizeof(T)*nelem*2);
    for(int i=0; i<nelem; i++){
      wfdata[2*i] = tmp[i];
      wfdata[2*i+1] = tmp[i+nelem];
    }

    interleaved = true;
    delete [] tmp;
  }

  template<class T>
    void diffractive_wavefront<T>::non_interleaved_conversion() const {

    // here we put wfdata[2*i] => wfdata[i]
    // and wfdata[2*i+1] => wfdata[i+nelem]

    if(!interleaved) return;

    int nelem = diffractive_wavefront_header<T>::total_space();
    T * tmp;

    try{tmp= new T[nelem];}
    catch(...){
      cerr << "diffractive_wavefront::non_interleaved_conversion error - "
	   << "unable to allocate memory\n";
      throw(string("diffractive_wavefront::non_interleaved_conversion"));
    }

    for(int i=0; i<nelem; i++)
      tmp[i] = wfdata[2*i+1];
  
    for(int i=0; i<nelem; i++)
      wfdata[i] = wfdata[2*i];

    for(int i=0; i<nelem; i++)
      wfdata[i+nelem] = tmp[i];

    interleaved = false;
    delete [] tmp;
  }
  
 template<class T>
    diffractive_wavefront<T>::diffractive_wavefront(){
    wfdata = NULL;
    real_imag = true;
    interleaved = true;
  }

  template<class T>
    diffractive_wavefront<T>::diffractive_wavefront(const diffractive_wavefront<T> & dwf){
    wfdata = NULL;
    this->operator=(dwf);
  }

  template<class T>
    diffractive_wavefront<T>::diffractive_wavefront(const char * filename){
    wfdata = NULL;
    interleaved = false;
    this->read(filename);
  }

  template<class T>
    diffractive_wavefront<T>::diffractive_wavefront(const iofits & iof){

    wfdata = NULL;
    interleaved = false;
    this->read(iof);
  }

  template<class T>
    diffractive_wavefront<T>::diffractive_wavefront(const diffractive_wavefront_header<T> & dwfh, 
						    T * data, bool rl_img, bool intrlvd) {

    real_imag = rl_img;
    interleaved = intrlvd;

    wfdata = NULL;
    vector<long> tmp = dwfh.get_axes();
    this->resize(dwfh.get_axes());
    if(data!=NULL){
      int nelem = 2*dwfh.get_axes()[0]*dwfh.get_axes()[1];
      for(int i=0; i<nelem; i++) wfdata[i] = data[i];
    }
    this->diffractive_wavefront_header<T>::operator=(dwfh);
  }

  template<class T>
    diffractive_wavefront<T>::~diffractive_wavefront(){
    delete [] wfdata;
  }

  template<class T>
    diffractive_wavefront<T> & diffractive_wavefront<T>::operator=(const diffractive_wavefront<T> & dwf){

    if(this==&dwf)
      return(*this);

    real_imag = dwf.real_imag;
    interleaved = dwf.interleaved;
    if(this->diffractive_wavefront_header<T>::total_space()!=dwf.diffractive_wavefront_header<T>::total_space())
      this->resize(dwf.get_axes());

    int nelem = dwf.diffractive_wavefront_header<T>::total_space();
    for(int i=0; i<2*nelem; i++)
      wfdata[i] = dwf.wfdata[i];

    this->diffractive_wavefront_header<T>::operator=(dwf); 
    this->fft_manager<T>::operator=(dwf);
    return(*this);
  }

  template<class T>
    void diffractive_wavefront<T>::mask(){
 
    int nelem = diffractive_wavefront_header<T>::total_space();
    if(interleaved){
      for(int i=0; i<nelem; i++)
	if(wfdata[2*i]==0) wfdata[2*i+1] = 0;
    } else {
      for(int i=0; i<nelem; i++)
	if(wfdata[i]==0) wfdata[i+nelem] = 0;
    }
  }

  template<class T>
    template<class U>
    void diffractive_wavefront<T>::install(const pixel_amp_array<U> & pixamparr){

    int nelem = diffractive_wavefront_header<T>::total_space();
    if(diffractive_wavefront_header<T>::axes!=pixamparr.get_axes()){
      cerr << "diffractive_wavefront::install error - " 
	   << "mismatched axes\n";
      cerr << "diffractive wavefront axes " << this->axes[0] << ", " << this->axes[1] << endl;
      cerr << "pixel amp array axes " << pixamparr.get_axes()[0] << ", " << pixamparr.get_axes()[1] << endl;
      throw(string("diffractive_wavefront::install"));
    }

    this->amp_phase_conversion();

    if(interleaved){
      for(int i=0; i<nelem; i++)
	wfdata[2*i] = pixamparr.data(i);
    } else {
      for(int i=0; i<nelem; i++)
	wfdata[i] = pixamparr.data(i);
    }
  }

  template<class T>
    template<class U>
    void diffractive_wavefront<T>::install(const pixel_phase_array<U> & pixpharr){

    int nelem = diffractive_wavefront_header<T>::total_space();
    if(diffractive_wavefront_header<T>::axes!=pixpharr.get_axes()){
      cerr << "diffractive_wavefront::install error - " 
	   << "mismatched axes\n";
      cerr << "diffractive wavefront axes " << this->axes[0] << ", " << this->axes[1] << endl;
      cerr << "pixel phase array axes " << pixpharr.get_axes()[0] << ", " << pixpharr.get_axes()[1] << endl;
      throw(string("diffractive_wavefront::install"));
    }

    this->amp_phase_conversion();

    if(interleaved){
      for(int i=0; i<nelem; i++) {
	if(wfdata[2*i]==0) wfdata[2*i+1]=0;
	else wfdata[2*i+1] = pixpharr.data(i);
      }
    } else {
      for(int i=0; i<nelem; i++){
	if(wfdata[i]==0) wfdata[i+nelem]=0;
	else wfdata[i+nelem] = pixpharr.data(i);
      }
    }
  }

  template<class T>
    pixel_amp_array<T> diffractive_wavefront<T>::extract_amps() const {

    int nelem = this->diffractive_wavefront_header<T>::total_space();
    T * amps;
    try{amps = new T[nelem];}
    catch(...){
      cerr << "diffractive_wavefront::extract_amps error - "
	   << "unable to allocate memory\n";
      throw(string("diffractive_wavefront::extract_amps"));
    }
    this->amp_phase_conversion();

    if(interleaved){
      for(int i=0; i<nelem; i++)
	amps[i] = wfdata[2*i];
    } else {
      for(int i=0; i<nelem; i++)
	amps[i] = wfdata[i];
    }

    pixel_amp_array<T> pixamparr(this->diffractive_wavefront_header<T>::axes, amps);
    delete [] amps;
    return(pixamparr);
  }

  template<class T>
    pixel_phase_array<T> diffractive_wavefront<T>::extract_phases() const {

    int nelem = this->diffractive_wavefront_header<T>::total_space();
    T * phases;
    try{phases = new T[nelem];}
    catch(...){
      cerr << "diffractive_wavefront::extract_phases error - "
	   << "unable to allocate memory\n";
      throw(string("diffractive_wavefront::extract_phases"));
    }

    this->amp_phase_conversion();

    if(interleaved){
      for(int i=0; i<nelem; i++)
	phases[i] = wfdata[2*i+1];
    } else {
      for(int i=0; i<nelem; i++)
	phases[i] = wfdata[i+nelem];
    }

    pixel_phase_array<T> pixpharr(this->diffractive_wavefront_header<T>::axes, phases);
    delete [] phases;
    return(pixpharr);
  }

  template<class T>
    template<class U>
    void diffractive_wavefront<T>::multiply(const pixel_amp_array<U> & pixamparr) {

    int nelem = diffractive_wavefront_header<T>::total_space();
    if(diffractive_wavefront_header<T>::axes!=pixamparr.get_axes()){
      cerr << "diffractive_wavefront::multiply error - " 
	   << "mismatched axes\n";
      throw(string("diffractive_wavefront::multiply"));
    }

    this->amp_phase_conversion();

    if(interleaved){
      for(int i=0; i<nelem; i++) 
	wfdata[2*i] *= pixamparr.data(i);
    } else {
      for(int i=0; i<nelem; i++)
	wfdata[i] *= pixamparr.data(i);
    }
  }

  template<class T>
    template<class U>
    void diffractive_wavefront<T>::add(const pixel_phase_array<U> & pixpharr) {

    int nelem = diffractive_wavefront_header<T>::total_space();
    if(diffractive_wavefront_header<T>::axes!=pixpharr.get_axes()){
      cerr << "diffractive_wavefront::add error - " 
	   << "mismatched axes\n";
      throw(string("diffractive_wavefront::add"));
    }

    this->amp_phase_conversion();

    if(interleaved){
      for(int i=0; i<nelem; i++) {
	if(wfdata[2*i]==0) wfdata[2*i+1]=0;
	else wfdata[2*i+1] += pixpharr.data(i);
      }
    } else {
      for(int i=0; i<nelem; i++){
	if(wfdata[i]==0) wfdata[i+nelem]=0;
	else wfdata[i+nelem] += pixpharr.data(i);
      }
    }
  }

  template<class T>
    void diffractive_wavefront<T>::read(const char * filename) {
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "diffractive_wavefront::read - "
	   << "error opening file " << filename << endl;
      throw(string("diffractive_wavefront::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "diffractive_wavefront::read - "
	   << "error reading diffractive wavefront from file " 
	   << filename << endl;
      throw(string("diffractive_wavefront::read"));
    }
  }

  template<class T>
    void diffractive_wavefront<T>::read(const iofits & iof) {

    // read into temporary header
    // since resize uses the data
    // residing in the diffractive_wavefront_header
    diffractive_wavefront_header<T> wfhdr(iof);
    this->resize(wfhdr.get_axes());
    this->operator=(wfhdr);
  
    if(wavefront::verbose_level>1)
      cout << "diffractive_wavefront::read - reading amplitudes\n";
    pixel_amp_array<T> pixamparr(iof);
    this->install(pixamparr);

    if(wavefront::verbose_level>1)
      cout << "diffractive_wavefront::read - reading phases\n";
    pixel_phase_array<T> pixpharr(iof);
    this->install(pixpharr);

  }

  template<class T>
    void diffractive_wavefront<T>::write(const char * filename) const {
    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "diffractive_wavefront::write - "
	   << "error opening file " << filename << endl;
      throw(string("diffractive_wavefront::write"));
    }
    try{this->write(iof);}
    catch(...){
      cerr << "diffractive_wavefront::write - "
	   << "error writing diffractive wavefront to file " 
	   << filename << endl;
      iof.print_header(cerr, "error");
      throw(string("diffractive_wavefront::write"));
    }
  }

  template<class T>
    void diffractive_wavefront<T>::write(iofits & iof) const {

    try{this->diffractive_wavefront_header<T>::write(iof);}
    catch(...){
      cerr << "diffractive_wavefront::write error - could not write header\n";
      throw(string("diffractive_wavefront::write"));
    }

    try{this->extract_amps().write(iof);}
    catch(...){
      cerr << "diffractive_wavefront::write error - could not write amps\n";
      throw(string("diffractive_wavefront::write"));
    }
    

    // Here we have to start a new hdu, since 
    // pixel_array will not start one for us
    fits_header_data<T> fhd(this->get_axes());
    fhd.write(iof);

    try{this->extract_phases().write(iof);}
    catch(...){
      cerr << "diffractive_wavefront::write error - could not write phases\n";
      throw(string("diffractive_wavefront::write"));
    }
  }

  template<class T> 
    void diffractive_wavefront<T>::print(ostream & os, const char * prefix) const {
    this->diffractive_wavefront_header<T>::print(os, prefix); 
  }

  template<class T>
    std::complex<T> diffractive_wavefront<T>::data(int n) const {

    if(n<0 || n>=this->diffractive_wavefront_header<T>::total_space()){
      cerr << "diffractive_wavefront<T>::data error - index " << n 
	   << " supplied to this function is out of range 0 - "
	   << this->diffractive_wavefront_header<T>::total_space()-1 << endl;
      throw(string("diffractive_wavefront<T>::data"));
    }

    // NOTE:  NEED TO FIX THIS UP TO ACCOUNT FOR CURVATURE, IF PRESENT
    // CHANGE angle definitions below to add curvature term

    std::complex<T> c;
    if(real_imag && interleaved) {
      c = std::complex<T>(wfdata[2*n], wfdata[2*n+1]);
    } else if(real_imag && !interleaved) {
      c = std::complex<T>(wfdata[n], wfdata[n+this->diffractive_wavefront_header<T>::total_space()]);
    } else if(!real_imag && interleaved) {
      double tmp = wfdata[2*n];
      double angle = wfdata[2*n+1];
      c = std::complex<T>(tmp*cos(angle), tmp*sin(angle));
    } else if(!real_imag && !interleaved) {
      double tmp = wfdata[n];
      double angle = wfdata[n+this->diffractive_wavefront_header<T>::total_space()];
      c = std::complex<T>(tmp*cos(angle), tmp*sin(angle));
    }
  
    return(c);
  }

  template<class T>
    void diffractive_wavefront<T>::set_data(int n, std::complex<T> & dat) {

    if(n<0 || n>=this->diffractive_wavefront_header<T>::total_space()){
      cerr << "diffractive_wavefront<T>::set_data error - index " << n 
	   << " supplied to this function is out of range 0 - "
	   << this->diffractive_wavefront_header<T>::total_space()-1 << endl;
      throw(string("diffractive_wavefront<T>::set_data"));
    }

    // NOTE:  NEED TO FIX THIS UP TO ACCOUNT FOR CURVATURE, IF PRESENT
    // CHANGE angle definitions below to add curvature term

    if(real_imag && interleaved) {
      wfdata[2*n] = dat.real();
      wfdata[2*n+1] = dat.imag();
    } else if(real_imag && !interleaved) {
      wfdata[n] = dat.real();
      wfdata[n+this->diffractive_wavefront_header<T>::total_space()] = dat.imag();
    } else if(!real_imag && interleaved) {
      wfdata[2*n] = abs(dat);
      wfdata[2*n+1] = arg(dat);
    } else if(!real_imag && !interleaved) {
      wfdata[n] = abs(dat);
      wfdata[n+this->diffractive_wavefront_header<T>::total_space()] = arg(dat);
    }
  }
  
  template<class T>
    void diffractive_wavefront<T>::set_propagation_direction(const three_vector & prop_dir){

    if(prop_dir.length()<three_frame::precision){
      cerr << "diffractive_wavefront::update error - "
	   << "propagation direction vector provided to this function is null\n";
      throw(string("diffractive_wavefront::update"));
    } 

    three_vector normalized_prop_dir = prop_dir * (1/prop_dir.length());

    if(cross_product(normalized_prop_dir, this->z()).length()<three_frame::precision)
      return;

    this->amp_phase_conversion();

    // modify wf phase to absorb tilt introduced by changing
    // propagation direction
    double x_halfpix=0, y_halfpix=0;
    int x_extrapix=1, y_extrapix=1;
    if(this->axes[1]%2==0){
      x_halfpix = .5;
      x_extrapix = 0;
    }
    if(this->axes[0]%2==0){
      y_halfpix = .5;
      y_extrapix = 0;
    }

    // Slopes are in radians of phase delay per pixel
    double xslope = 2*M_PI*this->get_pixel_scale()/this->get_wavelength()*
      dot_product(normalized_prop_dir, this->x());
    double yslope = 2*M_PI*this->get_pixel_scale()/this->get_wavelength()*
      dot_product(normalized_prop_dir, this->y());
    /*
    cout << "dot products " << dot_product(normalized_prop_dir, this->x())
	 << "\t" << dot_product(normalized_prop_dir, this->y()) << endl;
    cout << "xslope " << xslope << " yslope " << yslope << endl;
    */
    if(interleaved) {
      for(int i=-this->axes[1]/2; i<this->axes[1]/2+x_extrapix; i++)
	for(int j=-this->axes[0]/2; j<this->axes[0]/2+y_extrapix; j++)
	  if(wfdata[2*((i+this->axes[1]/2)*this->axes[0]+j+this->axes[0]/2)]!=0)
	    wfdata[2*((i+this->axes[1]/2)*this->axes[0]+j+this->axes[0]/2)+1] += 
	      (i+x_halfpix)*xslope + (j+y_halfpix)*yslope;
    } else {
      int nelem = this->axes[0]*this->axes[1];
      for(int i=-this->axes[1]/2; i<this->axes[1]/2+x_extrapix; i++)
	for(int j=-this->axes[0]/2; j<this->axes[0]/2+y_extrapix; j++)
	  if(wfdata[(i+this->axes[1]/2)*this->axes[0]+j+this->axes[0]/2]!=0)
	    wfdata[(i+this->axes[1]/2)*this->axes[0]+j+this->axes[0]/2+nelem+1] += 
	      (i+x_halfpix)*xslope + (j+y_halfpix)*yslope;
    }

    // Rotate the wf three_frame so that its z axis lies along prop_dir
    three_vector axis_of_rotation = cross_product(this->z(), normalized_prop_dir);
    double rotation_angle = asin(axis_of_rotation.length());
    //Arroyo::three_rotation trot(*this, axis_of_rotation, -rotation_angle);
    Arroyo::three_rotation trot(*this, axis_of_rotation, rotation_angle);
    trot.transform(*this);

  }

  template<class T>
    void diffractive_wavefront<T>::set_wavefront_curvature(double curvature){

    this->amp_phase_conversion();
    this->interleaved_conversion();
    
    // modify wf phase to absorb curvature in header

    double x_halfpix=0, y_halfpix=0;
    int x_extrapix=1, y_extrapix=1;
    if(this->axes[1]%2==0){
      x_halfpix = .5;
      x_extrapix = 0;
    }
    if(this->axes[0]%2==0){
      y_halfpix = .5;
      y_extrapix = 0;
    }

    // NOTE: positive curvature indicates diverging beam. 

    double curvature_term = 
      (this->get_curvature() - curvature)*M_PI*this->get_pixel_scale()*this->get_pixel_scale()/this->wavelength;
    
    int index;
    double twopi = 2*M_PI;

    for(int i=-this->axes[1]/2; i<this->axes[1]/2+x_extrapix; i++){
      for(int j=-this->axes[0]/2; j<this->axes[0]/2+y_extrapix; j++){
	index = 2*((i+this->axes[1]/2)*this->axes[0]+j+this->axes[0]/2);
	if(wfdata[index]!=0)
	  wfdata[index+1] =
	    fmod(wfdata[index+1] + 
		 curvature_term*((i+x_halfpix)*(i+x_halfpix)+(j+y_halfpix)*(j+y_halfpix)),
		 twopi);
      }
    }
  
    this->set_curvature(curvature);
  }

  template<class T>
    void diffractive_wavefront<T>::write_amps_to_ppm(double min, double max,
						     bool logscale,
						     bool colorbar,
						     colormap * cmap, 
						     const char * filename,
						     long min_dimen) const {

    this->extract_amps().write_to_ppm(min, max, logscale, colorbar, cmap, filename, min_dimen);
  }

  template<class T>
    void diffractive_wavefront<T>::write_phases_to_ppm(double min, double max, 
						       bool colorbar,
						       colormap * cmap, 
						       const char * filename,
						       long min_dimen) const {

    this->extract_phases().write_to_ppm(min, max, false, colorbar, cmap, filename, min_dimen);
  }

  template<class T>
    void diffractive_wavefront<T>::geometric_propagator(double distance) {
    this->three_point::operator+=(distance*this->z());
    if(this->curvature!=0){
      double initial_position = 1/this->curvature;
      this->pixel_scale *= (initial_position+distance)/initial_position;
      this->curvature = 1/(initial_position+distance);
    }
  } 

  template<class T>
    void diffractive_wavefront<T>::exact_propagator(double distance, 
						    double final_pixel_scale, 
						    vector<long> final_axes) {

    int nelem = this->diffractive_wavefront_header<T>::total_space();
    if(nelem==0){
      cerr << "diffractive_wavefront::exact_propagator error - "
	   << "cannot propagate this diffractive_wavefront, as it has no data\n";
      throw(string("diffractive_wavefront::exact_propagator"));
    }

    if(final_pixel_scale==-1){
      final_pixel_scale = this->pixel_scale;
      final_axes = this->axes;
    }

    if(this->pixel_scale<0){
      cerr << "diffractive_wavefront::exact_propagator error - "
	   << "requested pixel scale " << this->pixel_scale 
	   << " less than zero\n";
      throw(string("diffractive_wavefront::exact_propagator"));
    }

    if(final_axes.size()!=2 || final_axes[0]<=0 || final_axes[1]<=0){
      cerr << "diffractive_wavefront::exact_propagator error - "
	   << "requested axes invalid\n";
      for(int i=0; i<this->axes.size(); i++)
	cerr << "\taxis " << i << "\t" << final_axes[i] << endl;
      throw(string("diffractive_wavefront::exact_propagator"));
    }
   
    // this array holds the new pixel data
    // in non-interleaved format.  
    int new_nelem = final_axes[0]*final_axes[1];
    T * newdata;
    try{newdata = new T[2*new_nelem];}
    catch(...){
      cerr << "diffractive_wavefront::exact_propagator error - "
	   << "could not allocate memory\n";
      throw(string("diffractive_wavefront::exact_propagator"));
    }
      
    double initial_power;

    if(wavefront::verbose_level) initial_power = this->total_power();

    // the half pixel definitions - so that
    // when you have even axes the pixel
    // centroids lie on half pixel values.
    double old_x_halfpix=0, old_y_halfpix=0;
    double new_x_halfpix=0, new_y_halfpix=0;
    int old_x_extrapix=1, old_y_extrapix=1;
    int new_x_extrapix=1, new_y_extrapix=1;
    if(this->axes[1]%2==0){
      old_x_halfpix = .5;
      old_x_extrapix = 0;
    }
    if(this->axes[0]%2==0){
      old_y_halfpix = .5;
      old_y_extrapix = 0;
    }

    if(final_axes[1]%2==0){
      new_x_halfpix = .5;
      new_x_extrapix = 0;
    }
    if(final_axes[0]%2==0){
      new_y_halfpix = .5;
      new_y_extrapix = 0;
    }

    this->real_imag_conversion();
    this->interleaved_conversion();

    // Here we're unnecessarily computing the same sqrt 4 times - once
    // for each quadrant.  We should really do the loop over just one
    // quadrant to cut down on compute time.  However, the original data
    // is not necessarily symmetric - just the sqrt.  The minor difficulty
    // is that if we do this, we need to determine whether axes are even
    // or odd, and if odd treat the center row and column separately
    double pixel_to_pixel_distance;
    double phase, cos_phase, sin_phase;
    double distance_squared = distance*distance;
    double real, imag;
    double wavenumber = 2*M_PI/this->wavelength;
    double fac, obliquity;
    double twopi = M_PI * 2;
    int curvature_sign = this->curvature>0 ? 1 : -1;

    for(int i=-final_axes[1]/2; i<final_axes[1]/2+new_x_extrapix; i++){
      for(int j=-final_axes[0]/2; j<final_axes[0]/2+new_y_extrapix; j++){

	real = imag = 0;

	for(int k=-this->axes[1]/2; k<this->axes[1]/2+old_x_extrapix; k++){
	  for(int l=-this->axes[0]/2; l<this->axes[0]/2+old_y_extrapix; l++){
	  
	    if(wfdata[2*((k+this->axes[1]/2)*this->axes[0]+l+this->axes[0]/2)]==0 && 
	       wfdata[2*((k+this->axes[1]/2)*this->axes[0]+l+this->axes[0]/2)+1]==0)
	      continue;

	    pixel_to_pixel_distance = 
	      sqrt(distance_squared +
		   ((i+new_x_halfpix)*final_pixel_scale-(k+old_x_halfpix)*this->pixel_scale)*
		   ((i+new_x_halfpix)*final_pixel_scale-(k+old_x_halfpix)*this->pixel_scale) +
		   ((j+new_y_halfpix)*final_pixel_scale-(l+old_y_halfpix)*this->pixel_scale)*
		   ((j+new_y_halfpix)*final_pixel_scale-(l+old_y_halfpix)*this->pixel_scale));
	  
	    obliquity = distance/pixel_to_pixel_distance;
	    phase = pixel_to_pixel_distance*wavenumber - M_PI_2;

	    // Term to account for the extra propagation distance if
	    // wavefront has finite curvature.  The exact quantity is
	    // sqrt(lateral_displacement^{2}+1/(curvature^{2})) -
	    // 1/curvature.
	    if(this->curvature)
	      phase += curvature_sign*
		(sqrt(this->pixel_scale*this->pixel_scale*((k+old_x_halfpix)*(k+old_x_halfpix)+
							   (l+old_y_halfpix)*(l+old_y_halfpix)) + 
		      1/(this->curvature*this->curvature)) - 1/this->curvature)/this->wavelength;
	    
	    cos_phase = cos(phase);
	    sin_phase = sin(phase);
	    fac = obliquity*this->pixel_scale*this->pixel_scale/pixel_to_pixel_distance/this->wavelength;
	  
	    real += fac *
	      (wfdata[2*((k+this->axes[1]/2)*this->axes[0]+l+this->axes[0]/2)]*cos_phase - 
	       wfdata[2*((k+this->axes[1]/2)*this->axes[0]+l+this->axes[0]/2)+1]*sin_phase);
	    imag += fac *
	      (wfdata[2*((k+this->axes[1]/2)*this->axes[0]+l+this->axes[0]/2)]*sin_phase +
	       wfdata[2*((k+this->axes[1]/2)*this->axes[0]+l+this->axes[0]/2)+1]*cos_phase);

	    /*
	      if(i==0&&j==0)
	      cout << setw(5) << k << setw(5) << l 
	      << setw(12) << obliquity 
	      << setw(12) << fmod(phase,twopi)
	      << setw(12) << fac 
	      << setw(12) << wfdata[2*((k+axes[1]/2)*axes[0]+l+axes[0]/2)]
	      << setw(12) << wfdata[2*((k+axes[1]/2)*axes[0]+l+axes[0]/2)+1]
	      << setw(12) << real 
	      << setw(12) << imag << endl;
	    */
	  }
	}	

	newdata[2*((i+final_axes[1]/2)*final_axes[0]+j+final_axes[0]/2)] = real;
	newdata[2*((i+final_axes[1]/2)*final_axes[0]+j+final_axes[0]/2)+1] = imag;

      } 
    }

    delete [] this->wfdata;
    this->wfdata = newdata;
    this->axes = final_axes;

    if(wavefront::verbose_level) {
      double final_power = this->total_power();
      cout << "initial power " << initial_power << " final power " << final_power
	   << " ratio "  << final_power / initial_power << endl;
    }

    this->pixel_scale = final_pixel_scale;
    this->curvature = 0;

    three_translation tt(distance*this->three_frame::z());
    tt.transform(*this);
    update_timestamp(distance);
  }

  template<class T>
    void diffractive_wavefront<T>::forward_fft() {
    vector<long> flipped_axes(2, this->axes[0]);
    flipped_axes[0] = this->axes[1];
    //this->mk_fftw_plan(flipped_axes, false, true, FFTW_FORWARD);
    this->real_imag_conversion();
    this->interleaved_conversion();  
    this->fft_manager<T>::forward_fft(flipped_axes, false, true, wfdata);
  }

  template<class T>
    void diffractive_wavefront<T>::backward_fft(){
    vector<long> flipped_axes(2, this->axes[0]);
    flipped_axes[0] = this->axes[1];
    //this->mk_fftw_plan(flipped_axes, false, true, FFTW_BACKWARD);
    this->real_imag_conversion();
    this->interleaved_conversion();  
    this->fft_manager<T>::backward_fft(flipped_axes, false, true, wfdata);
  }

  namespace {

    // A convenience function to apply the near field transfer function
    // to the data
    //
    // This function requires that the input data be interleaved
    template<class T>
      void apply_transfer_function(double axial_phase, vector<long> axes, 
				   vector<double> phase_coefficient, bool fresnel, T * wfdata){

      double twopi = 2*M_PI;
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

      int index;
      double tmp_val, xfer_fn_phase;
      Arroyo::complex_cyclic_permutation(axes, axes[0]/2, axes[1]/2, wfdata);
      for(int i=-axes[1]/2; i<axes[1]/2+x_extrapix; i++){
	for(int j=-axes[0]/2; j<axes[0]/2+y_extrapix; j++){
	  index = (i+axes[1]/2)*axes[0]+j+axes[0]/2;
	  tmp_val = phase_coefficient[1]*(i+x_halfpix)*(i+x_halfpix) + 
	    phase_coefficient[0]*(j+y_halfpix)*(j+y_halfpix);
	  
	  if(tmp_val>1){
	    wfdata[2*i*axes[0]] = 0;
	    wfdata[2*(axes[1]-i)*axes[0]] = 0;       
	  } else {
	    if(fresnel) xfer_fn_phase = fmod((axial_phase*(1 - .5*tmp_val)),twopi);
	    else xfer_fn_phase = fmod((axial_phase*sqrt(1 - tmp_val)),twopi);
	    wfdata[2*index+1] += xfer_fn_phase;
	  }
	}
      }
      Arroyo::complex_cyclic_permutation(axes, -axes[0]/2, -axes[1]/2, wfdata);
    }

    // A convenience function to apply the near field transfer function
    // to the data
    //
    // This function requires that the input data be interleaved
    template<class T>
      void optimized_apply_transfer_function(double axial_phase, vector<long> axes, 
					     vector<double> phase_coefficient, bool fresnel, T * wfdata){

      // In the following application of the transfer function, we check
      // to see whether the sqrt is going to have a negative argument.
      // This condition indicates evanescent waves, so we set the transfer
      // function to zero.  See Goodman "Introduction to Fourier Optics"
      // eq. 3-47 for more explanation
      
      // First, do the transfer function for the dc elements
      // if one or both axes are even, then at the same time we do nyquist
      bool odd_zero = axes[0]%2;
      int lim_one = axes[1]/2+axes[1]%2;
      double xfer_fn_phase;
      double tmp_val;
      double twopi = 2*M_PI;

      for(int i=1; i<lim_one; i++){

	tmp_val = phase_coefficient[1]*i*i;
	if(tmp_val>1){
	  wfdata[2*i*axes[0]] = 0;
	  wfdata[2*(axes[1]-i)*axes[0]] = 0;       
	} else {
	  if(fresnel) xfer_fn_phase = fmod((axial_phase*(1 - .5*tmp_val)), twopi);
	  else xfer_fn_phase = fmod((axial_phase*sqrt(1 - tmp_val)), twopi);
	  wfdata[2*i*axes[0]+1] += xfer_fn_phase;
	  wfdata[2*(axes[1]-i)*axes[0]+1] += xfer_fn_phase;
	}

	if(!odd_zero){
	  tmp_val = phase_coefficient[1]*i*i + phase_coefficient[0]*axes[0]/2*axes[0]/2;
	  if(tmp_val>1){
	    wfdata[2*(i*axes[0]+axes[0]/2)] = 0; 
	    wfdata[2*((axes[1]-i)*axes[0]+axes[0]/2)] = 0;
	  } else {
	    if(fresnel) xfer_fn_phase = fmod((axial_phase*(1 - .5*tmp_val)), twopi);
	    else xfer_fn_phase = fmod((axial_phase*sqrt(1 - tmp_val)),twopi);
	    wfdata[2*(i*axes[0]+axes[0]/2)+1] += xfer_fn_phase; 
	    wfdata[2*((axes[1]-i)*axes[0]+axes[0]/2)+1] += xfer_fn_phase; 
	  }
	}
      }    
      
      bool odd_one = axes[1]%2;
      int lim_zero = axes[0]/2+axes[0]%2;
      for(int j=1; j<lim_zero; j++){
	
	tmp_val = phase_coefficient[0]*j*j;
	if(tmp_val>1){
	  wfdata[2*j] = 0;
	  wfdata[2*(axes[0]-j)] = 0;
	} else {
	  if(fresnel) xfer_fn_phase = fmod((axial_phase*(1 - .5*tmp_val)),twopi);
	  else xfer_fn_phase = fmod((axial_phase*sqrt(1 - tmp_val)),twopi);
	  wfdata[2*j+1] += xfer_fn_phase;
	  wfdata[2*(axes[0]-j)+1] += xfer_fn_phase;
	}
	
	if(!odd_one){
	  tmp_val = phase_coefficient[0]*j*j + phase_coefficient[1]*axes[1]/2*axes[1]/2;
	  if(tmp_val>1){
	    wfdata[2*(axes[1]/2*axes[0]+j)] = 0;
	    wfdata[2*(axes[1]/2*axes[0]+axes[0]-j)] = 0;
	  } else {
	    if(fresnel) xfer_fn_phase = fmod((axial_phase*(1 - .5*tmp_val)),twopi);
	    else xfer_fn_phase = fmod((axial_phase*sqrt(1 - tmp_val)),twopi);
	    wfdata[2*(axes[1]/2*axes[0]+j)+1] += xfer_fn_phase;
	    wfdata[2*(axes[1]/2*axes[0]+axes[0]-j)+1] += xfer_fn_phase;
	  }
	}
      }
      // The last 4 points we missed by indexing 
      // loops from one rather than from zero
      // This indexing was required because DC and nyquist
      // must be treated differently by the transfer function,
      // whereas all other points have symmetric properties
      // when i=>-i or j=>-j
      wfdata[1] += axial_phase;
      if(!odd_zero){
	tmp_val = phase_coefficient[0]*axes[0]/2*axes[0]/2;
	if(tmp_val>1) wfdata[axes[0]] = 0;
	else {
	  if(fresnel) wfdata[axes[0]+1] += fmod((axial_phase*(1 - .5*tmp_val)),twopi);
	  else wfdata[axes[0]+1] += fmod((axial_phase*sqrt(1 - tmp_val)),twopi);
	}
      }
      if(!odd_one){
	tmp_val = phase_coefficient[1]*axes[1]/2*axes[1]/2;
	if(tmp_val>1) wfdata[2*(axes[1]/2*axes[0])] = 0;
	else {
	  if(fresnel) wfdata[2*(axes[1]/2*axes[0])+1] += fmod((axial_phase*(1 - .5*tmp_val)),twopi);
	  else wfdata[2*(axes[1]/2*axes[0])+1] += fmod((axial_phase*sqrt(1 - tmp_val)),twopi);
	}
      }
      if(!odd_zero && !odd_one){
	tmp_val = phase_coefficient[0]*axes[0]/2*axes[0]/2 + phase_coefficient[1]*axes[1]/2*axes[1]/2;
	if(tmp_val>1) wfdata[2*(axes[1]/2*axes[0]+axes[0]/2)] = 0;
	else {
	  if(fresnel) wfdata[2*(axes[1]/2*axes[0]+axes[0]/2)+1] += fmod((axial_phase*(1 - .5*tmp_val)),twopi);
	  else wfdata[2*(axes[1]/2*axes[0]+axes[0]/2)+1] += fmod((axial_phase*sqrt(1 - tmp_val)),twopi);
	}
      }
      
      // Now do the rest of the data - easy since
      // the transfer function in the 4 quadrants 
      // is symmetric
      for(int i=1; i<lim_one; i++){
	for(int j=1; j<lim_zero; j++){
	  tmp_val = phase_coefficient[1]*i*i + phase_coefficient[0]*j*j;
	  if(tmp_val>1){
	    wfdata[2*(i*axes[0]+j)] = 0;
	    wfdata[2*((axes[1]-i)*axes[0]+j)] = 0;
	    wfdata[2*(i*axes[0]+axes[0]-j)] = 0;
	    wfdata[2*((axes[1]-i)*axes[0]+axes[0]-j)] = 0;
	  } else {
	    if(fresnel) xfer_fn_phase = fmod((axial_phase*(1 - .5*tmp_val)),twopi);
	    else xfer_fn_phase = fmod((axial_phase*sqrt(1 - tmp_val)),twopi);
	    wfdata[2*(i*axes[0]+j)+1] += xfer_fn_phase;
	    wfdata[2*((axes[1]-i)*axes[0]+j)+1] += xfer_fn_phase;
	    wfdata[2*(i*axes[0]+axes[0]-j)+1] += xfer_fn_phase;
	    wfdata[2*((axes[1]-i)*axes[0]+axes[0]-j)+1] += xfer_fn_phase;
	  }
	}
      }
    }
  }

  template<class T>
    void diffractive_wavefront<T>::near_field_propagator(double distance, bool fresnel) {

    // NOTES:
    //   For verging beams, remember that you have to operate
    //   on a diffractive_wavefront from which the curvature has been removed.
    //   In this case, you do the regular thing but use a different
    //   distance.  Then you fix things up at the end.

    if(wavefront::verbose_level>2){
      string fname = "angular_pre_backwards_fft_diffractive_wavefront.fits";
      if(fresnel)
	fname = "fresnel_pre_backwards_fft_diffractive_wavefront.fits";
      cout << "diffractive_wavefront::near_field_propagator - writing initial diffractive_wavefront to file " 
	   << fname << endl;
      this->write(fname.c_str());
    }

    // Perform inverse fourier transform
    if(diffractive_wavefront::wavefront::verbose_level)
      cout << "diffractive_wavefront::near_field_propagator - performing backwards fft\n";

    this->forward_fft();

    this->amp_phase_conversion();
    this->interleaved_conversion();

    // This factor will multiply the indices to produce the correct
    // value in the exponential phase term.  It's effect includes
    // conversion from pixels to frequency interval, plus the prefactor
    // of M_PI*dist*wavelength.  The conversion to frequency interval
    // depends on the dimensionality of the axis, so this factor may
    // differ between the two axes See GLAD manual page 4.15
    vector<double> phase_coefficient(2);
    phase_coefficient[0] = this->wavelength*this->wavelength/(this->pixel_scale*this->pixel_scale*this->axes[0]*this->axes[0]); 
    phase_coefficient[1] = this->wavelength*this->wavelength/(this->pixel_scale*this->pixel_scale*this->axes[1]*this->axes[1]);

    // this is the phase of light travelling directly down the optical axis
    double axial_phase = 2*M_PI*distance/this->wavelength;
    double effective_distance = distance/(this->curvature*distance+1);

    if(this->curvature)
      axial_phase = 2*M_PI*effective_distance/this->wavelength;

    int nelem = this->axes[0]*this->axes[1];

    //apply_transfer_function(axial_phase, axes, phase_coefficients, fresnel, wfdata);
    optimized_apply_transfer_function(axial_phase, this->axes, phase_coefficient, fresnel, wfdata);

    // Perform fourier transform
    this->backward_fft();

    if(wavefront::verbose_level>2){
      string fname = "angular_forwards_fft_xfer_diffractive_wavefront.fits";
      if(fresnel) fname = "fresnel_post_backwards_fft_diffractive_wavefront.fits";
      cout << "diffractive_wavefront::near_field_propagator - "
	   << "writing back-transformed diffractive_wavefront with transfer function applied to file " 
	   << fname << endl;
      this->write(fname.c_str());
    }

    // FFTW leaves a factor of N*M 
    // after forward and backward fft
    double inv = 1/(1+distance*this->curvature)/(double)(this->axes[0]*this->axes[1]);
    for(int i=0; i<2*nelem; i++)
      wfdata[i] *= inv;

    if(this->curvature){
      this->pixel_scale *= (1+distance*this->curvature);
      this->curvature = 1/(distance+1/this->curvature);
    }

    this->three_point::operator+=(distance*this->three_frame::z());
    update_timestamp(distance);
  }

  template<class T>
    void diffractive_wavefront<T>::near_field_angular_propagator(double distance) {
    this->near_field_propagator(distance, false);
  }

  template<class T>
    void diffractive_wavefront<T>::near_field_fresnel_propagator(double distance) {
    this->near_field_propagator(distance, true);
  }


  // Be very careful about changing this function - four different
  // propagators rely on it.
  // NOTE:  when you add in the curvature data member
  //        to diffractive_wavefront, you will have to include 
  //        operations using this term in this function
  template<class T>
    void diffractive_wavefront<T>::far_field_propagator(double distance, 
							bool fresnel, 
							double final_pixel_scale, 
							vector<long> final_axes) {


    bool goertzel_reinsch = false;
    if(final_axes.size()==2 && final_pixel_scale>0)
      goertzel_reinsch = true;
    else if((final_axes.size()!=0 && final_axes.size()!=2) ||
	    final_pixel_scale!=-1){
      cout << "diffractive_wavefront::far_field_propagator error - " 
	   << "this function was invoked by implicitly requesting " << endl
	   << "propagation using the Goertzel-Reinsch algorithm - " << endl
	   << "but this invocation was corrupt."  << endl
	   << "Specifically, the requested final pixel scale was " 
	   << final_pixel_scale << " and the requested axes were\n";
      for(int i=0; i<final_axes.size(); i++)
	cout << final_axes[i] << endl;
      throw(string("diffractive_wavefront::far_field_propagator"));
    }

    if(wavefront::verbose_level){
      if(fresnel && goertzel_reinsch)
	cout << "diffractive_wavefront::far_field_propagator - applying fresnel goertzel reinsch propagator\n";
      else if(!fresnel && goertzel_reinsch)
	cout << "diffractive_wavefront::far_field_propagator - applying fraunhoffer goertzel reinsch propagator\n";
      else if(fresnel && !goertzel_reinsch)
	cout << "diffractive_wavefront::far_field_propagator - applying fresnel propagator\n";
      else
	cout << "diffractive_wavefront::far_field_propagator - applying fraunhoffer propagator\n";
    }

    if(!goertzel_reinsch && (this->axes[0]!=this->axes[1])){
      cerr << "diffractive_wavefront::far_field_propagator error - "
	   << "cannot perform a far field propagation via fft " << endl
	   << "on an input array with axes of different size, "
	   << "because the software can only handle square pixels\n";
      cerr << "If you need to do this, try the Goertzel-Reinsch propagators\n";
      throw(string("diffractive_wavefront::far_field_propagator"));
    }

    this->interleaved_conversion();
    this->amp_phase_conversion();

    // Correct for halfpixel values At the same time, adjust the phase
    // slope to generate half a pixel shift in the image plane, so that
    // DC is not off-center For fresnel propagation, we also apply the
    // quadratic phase term.  Finally, we renormalize the amplitudes
    //
    // There's a subtlety here with half pixels DC is at 0,0.  We need
    // to put the center of the diffractive_wavefront there.  If the
    // diffractive_wavefront has an even number of pixels, we end up
    // with a phase slope in the image plane corresponding to half a
    // pixel shift in the pupil plane.  Likewise, after the
    // transformation we end up with a DC term at 0,0, which corresponds
    // to the center of the array.  Again if there are an even number of
    // pixels across the array we end up with a half pixel shift after
    // reordering the quadrants.  So we use the shift theorem to fix
    // things up, taking out a phase slope in the pupil plane to account
    // for the image plane half pixel shift, and subtracting a phase
    // slope in the image plane to account for the half pixel shift in
    // the pupil plane.  These effects cancel out in the near field
    // propagator case because we do both a forward and backward fft.

    vector<long> initial_axes(this->axes);

    double initial_x_halfpix=0, initial_y_halfpix=0;
    int initial_x_extrapix=1, initial_y_extrapix=1;
    if(initial_axes[1]%2==0){
      initial_x_halfpix = .5;
      initial_x_extrapix = 0;
    }
    if(initial_axes[0]%2==0){
      initial_y_halfpix = .5;
      initial_y_extrapix = 0;
    }

    if(final_axes.size()==0) final_axes = initial_axes;

    double final_x_halfpix=0, final_y_halfpix=0;
    int final_x_extrapix=1, final_y_extrapix=1;
    if(final_axes[1]%2==0){
      final_x_halfpix = .5;
      final_x_extrapix = 0;
    }
    if(final_axes[0]%2==0){
      final_y_halfpix = .5;
      final_y_extrapix = 0;
    }

    // definitions to account for halfpixel shifts for even final array dimensions
    double initial_pixel_scale = this->get_pixel_scale();
    double nyquist_pixel_scale = this->wavelength * distance / (double) initial_axes[0] / initial_pixel_scale;
    if(final_pixel_scale==-1) final_pixel_scale = nyquist_pixel_scale;
    double xslope = final_axes[1]%2==1 ? 0 : M_PI*final_pixel_scale/nyquist_pixel_scale/(double)initial_axes[1];
    double yslope = final_axes[0]%2==1 ? 0 : M_PI*final_pixel_scale/nyquist_pixel_scale/(double)initial_axes[0];

    // Term to account for the extra propagation distance if
    // wavefront has finite curvature.  The exact quantity is
    // sqrt(lateral_displacement^{2}+1/(curvature^{2})) -
    // 1/curvature, but here we approximate this by
    // .5*lateral_displacement^{2}*curvature^{2}
    int index;
    double twopi = 2*M_PI;
    double initial_curvature_term = M_PI*initial_pixel_scale*initial_pixel_scale*this->curvature/this->wavelength;
    double quadratic_coefficient = 0;
    if(fresnel)
      quadratic_coefficient = M_PI*initial_pixel_scale*initial_pixel_scale/this->wavelength/distance;

    for(int i=-initial_axes[1]/2; i<initial_axes[1]/2+initial_x_extrapix; i++){
      for(int j=-initial_axes[0]/2; j<initial_axes[0]/2+initial_y_extrapix; j++){
	index = (i+initial_axes[1]/2)*initial_axes[0]+j+initial_axes[0]/2;

	wfdata[2*index+1] += 
	  fmod(((quadratic_coefficient + initial_curvature_term)*
		((i+initial_x_halfpix)*(i+initial_x_halfpix)+(j+initial_y_halfpix)*(j+initial_y_halfpix))
		- xslope*(i+initial_x_halfpix) - yslope*(j+initial_y_halfpix)), twopi);
      }
    }

    if(wavefront::verbose_level>2){
      string fname = "far_field_fraunhoffer_xfer_diffractive_wavefront.fits";
      if(fresnel && goertzel_reinsch) fname = "far_field_fresnel_goertzel_reinsch_xfer_diffractive_wavefront.fits";
      else if(!fresnel && goertzel_reinsch) fname = "far_field_fraunhoffer_goertzel_reinsch_xfer_diffractive_wavefront.fits";
      else if(fresnel && !goertzel_reinsch) fname = "far_field_fresnel_xfer_diffractive_wavefront.fits";
      cout << "diffractive_wavefront::far_field_propagator - "
	   << "writing diffractive_wavefront with transfer function applied to file " 
	   << fname << endl;
      this->write(fname.c_str());
    }

    // for the normalization factor, we have one factor of
    // 1/wavelength/distance for the prefactor, two more for changing
    // variables in the fft.  Then we have two factors of the pixel
    // scale for each of the two dimensions
    double normalization_factor = 
      fabs((double)(initial_pixel_scale*initial_pixel_scale/this->wavelength/distance));

    if(goertzel_reinsch){
      long dimen = 2*final_axes[0]*final_axes[1];
      T * newdata;
      try{newdata = new T[dimen];}
      catch(...){
	cerr << "diffractive_wavefront::far_field_propagator error - "
	     << "unable to allocate memory\n";
	throw(string("diffractive_wavefront::far_field_propagator"));
      }

      this->real_imag_conversion();
      // negative sampling factor for forward fft
      double sampling_factor = -initial_pixel_scale*final_pixel_scale / this->wavelength / distance;
      goertzel_reinsch_transform(initial_axes, final_axes, sampling_factor, wfdata, newdata);
      delete [] wfdata;
      wfdata = newdata;
      this->axes = final_axes;
    } else {
      this->cyclic_permutation(initial_axes[0]/2, initial_axes[1]/2);
      // Perform inverse fourier transform
      if(wavefront::verbose_level)
	cout << "diffractive_wavefront::far_field_propagator - performing backwards fft\n";
      this->forward_fft();
      this->cyclic_permutation(-initial_axes[0]/2, -initial_axes[1]/2);
    }

    if(wavefront::verbose_level>2){
      string fname = "far_field_fraunhoffer_xformed_diffractive_wavefront.fits";
      if(fresnel && goertzel_reinsch) fname = "far_field_fresnel_goertzel_reinsch_xformed_diffractive_wavefront.fits";
      else if(!fresnel && goertzel_reinsch) fname = "far_field_fraunhoffer_goertzel_reinsch_xformed_diffractive_wavefront.fits";
      else if(fresnel && !goertzel_reinsch) fname = "far_field_fresnel_xformed_diffractive_wavefront.fits";
      cout << "diffractive_wavefront::far_field_propagator - "
	   << "writing transformed diffractive_wavefront to file " 
	   << fname << endl;
      this->write(fname.c_str());
    }


    // redefinitions to account for halfpixel shifts for even initial array dimensions
    nyquist_pixel_scale = this->wavelength * distance / (double) final_axes[0] / final_pixel_scale;
    xslope = initial_axes[1]%2==1 ? 0 : M_PI*initial_pixel_scale/nyquist_pixel_scale/(double)final_axes[1];
    yslope = initial_axes[0]%2==1 ? 0 : M_PI*initial_pixel_scale/nyquist_pixel_scale/(double)final_axes[0];

    // Apply the quadratic phase term to the diffractive_wavefront
    // this time, the normalization factor is used to 
    // normalize away the dimensional factor that arises
    // in a forwards to backwards fft
    this->amp_phase_conversion();
    quadratic_coefficient = final_pixel_scale*final_pixel_scale*M_PI/this->wavelength/distance;
    double constant_phase = 2*M_PI*fmod(distance/this->wavelength,1.0) - M_PI_2;
    for(int i=-final_axes[1]/2; i<final_axes[1]/2+final_x_extrapix; i++){
      for(int j=-final_axes[0]/2; j<final_axes[0]/2+final_y_extrapix; j++){
	index = (i+final_axes[1]/2)*final_axes[0]+j+final_axes[0]/2;
	wfdata[2*index] *= normalization_factor;
	wfdata[2*index+1] += constant_phase - xslope*(i + final_x_halfpix) - yslope*(j + final_y_halfpix);
	if(fresnel)
	  wfdata[2*index+1] = 
	    fmod(wfdata[2*index+1] + M_PI + quadratic_coefficient*
		 ((i+final_x_halfpix)*(i+final_x_halfpix)+(j+final_y_halfpix)*(j+final_y_halfpix)), twopi) - M_PI_2;
      }
    }
    
    this->pixel_scale = final_pixel_scale;
    this->three_point::operator+=(distance*this->three_frame::z());
    this->curvature = 0;
    update_timestamp(distance);
  }

  template<class T>
    void diffractive_wavefront<T>::far_field_fresnel_propagator(double distance) {
    this->far_field_propagator(distance, true);
  }

  template<class T>
    void diffractive_wavefront<T>::far_field_fraunhoffer_propagator(double distance) {
    this->far_field_propagator(distance, false);
  }

  template<class T>
    void diffractive_wavefront<T>::far_field_fresnel_goertzel_reinsch_propagator(double distance, 
										 double final_pixel_scale, 
										 vector<long> final_axes) {
    this->far_field_propagator(distance, true, final_pixel_scale, final_axes);
  }

  template<class T>
    void diffractive_wavefront<T>::far_field_fraunhoffer_goertzel_reinsch_propagator(double distance, 
										     double final_pixel_scale, 
										     vector<long> final_axes) {
    this->far_field_propagator(distance, false, final_pixel_scale, final_axes);
  }

  template<class T>
    void diffractive_wavefront<T>::finite_difference_method_propagator(double distance) {

    // fast algorithm for short distance propagation
    cerr << "diffractive_wavefront::finite_difference_method_propagator error - "
	 << "this function has not yet been coded\n";
    throw(string("diffractive_wavefront::finite_difference_method_propagator"));

  }

  template<class T>
    void diffractive_wavefront<T>::rotate(double angle) {

    //   To rotate a 2d data set through an angle a,
    //   you want to multiply each (x,y) pixel value
    //   by the matrix
    //  
    //   [  cos(a)   -sin(a)  ]
    //   [  sin(a)    cos(a)  ]
    //
    //   This is equivalent to multiplying by the 3 matrices
    // 
    //   [   1   -tan(a/2)  ] [   1      0][   1     -tan(a/2)  ]
    //   [   0      1       ] [ sin(a)   1][   0         1      ]
    //
    //  This multiplication may be performed in separate passes.
    //  It performs the rotation in place, and requires no
    //  reordering of the array.  However, it does rotate straight
    //  lines into bumpy, pixel resolved diagonal lines
    cerr << "diffractive_wavefront::rotate not yet coded\n";
    throw(string("diffractive_wavefront::rotate"));

  }

  template<class T>
    template<class U>
    void diffractive_wavefront<T>::pad_array(int npad, std::complex<U> value) {

    // NOTES:
    // This function still works if the diffractive_wavefront starts with no elements

    if(npad==0) return;
    if(npad<0){
      cerr << "diffractive_wavefront::pad_array error - cannot pad by " << npad << " pixels\n";
      throw(string("diffractive_wavefront::pad_array"));
    }

    vector<long> new_axes = this->axes;
    int old_nelem = diffractive_wavefront_header<T>::total_space();
    int nelem=1;
    for(int i=0; i<new_axes.size(); i++){
      new_axes[i] += 2*npad;
      nelem *= new_axes[i];
    }

    if(wavefront::verbose_level)
      cout << "diffractive_wavefront::pad_array - padding by " << npad 
	   << " from " << this->axes[0] << "x" << this->axes[1]
	   << " to " << new_axes[0] << "x" << new_axes[1] << endl;

    this->interleaved_conversion();

    T * olddata = wfdata;
    try{wfdata = new T[2*nelem];}
    catch(...){
      cerr << "diffractive_wavefront::pad_array error - "
	   << "unable to allocate memory\n";
      throw(string("diffractive_wavefront::pad_array"));
    }

    if(this->real_imag){
      U real = value.real();
      U imag = value.imag();
      for(int i=0; i<nelem; i++){
	wfdata[2*i] = real;
	wfdata[2*i+1] = imag;
      }
    } else {
      U amp = abs(value);
      U phase = arg(value);
      for(int i=0; i<nelem; i++){
	wfdata[2*i] = amp;
	wfdata[2*i+1] = phase;
      }
    }      

    // copy the old data into the new array
    for(int i=0; i<this->axes[1]; i++)
      for(int j=0; j<this->axes[0]; j++){
	wfdata[2*((i+npad)*new_axes[0]+npad+j)] = 
	  olddata[2*(i*this->axes[0]+j)];
	wfdata[2*((i+npad)*new_axes[0]+npad+j)+1] = 
	  olddata[2*(i*this->axes[0]+j)+1];
      }

    delete [] olddata;
    this->axes = new_axes;
  }

  template<class T>
    void diffractive_wavefront<T>::clip_array(int nclip) {

    if(nclip==0) return;
    if(nclip<0){
      cerr << "diffractive_wavefront::clip_array error - cannot clip by " << nclip << " pixels\n";
      throw(string("diffractive_wavefront::clip_array"));
    }

    vector<long> new_axes = this->axes;
    int nelem = 1;
    for(int i=0; i<new_axes.size(); i++){
      new_axes[i] -= 2*nclip;
      if(new_axes[i]<=0){
	cerr << "diffractive_wavefront::clip_array error - clipping original array by "
	     << nclip << " pixels leaves non-positive array size\n";
	throw(string("diffractive_wavefront::clip_array"));
      }
      nelem *= new_axes[i];
    }
    
    if(wavefront::verbose_level)
      cout << "diffractive_wavefront::clip_array - clipping by " << nclip
	   << " from " << this->axes[0] << "x" << this->axes[1]
	   << " to " << new_axes[0] << "x" << new_axes[1] << endl;
  
    T * newdata;
    try{newdata = new T[2*nelem];}
    catch(...){
      cerr << "diffractive_wavefront::pad_array error - "
	   << "unable to allocate memory\n";
      throw(string("diffractive_wavefront::clip_array"));
    }

    this->interleaved_conversion();

    // copy the old data into the new array
    for(int i=0; i<new_axes[1]; i++)
      for(int j=0; j<new_axes[0]; j++){
	newdata[2*(i*new_axes[0]+j)] = 
	  wfdata[2*((i+nclip)*this->axes[0]+nclip+j)];
	newdata[2*(i*new_axes[0]+j)+1] = 
	  wfdata[2*((i+nclip)*this->axes[0]+nclip+j)+1];
      }

    // make the switcheroo
    delete [] wfdata;
    wfdata = newdata;

    this->axes = new_axes;
  } 

  template<class T>
    double diffractive_wavefront<T>::total_power() const {

    double total_power = 0;
    int nelem = diffractive_wavefront_header<T>::total_space();
    if(real_imag){
      for(int i=0; i<2*nelem; i++)
	total_power += wfdata[i]*wfdata[i];
    } else {
      if(interleaved)
	for(int i=0; i<nelem; i++)
	  total_power += wfdata[2*i]*wfdata[2*i];
      else 
	for(int i=0; i<nelem; i++)
	  total_power += wfdata[i]*wfdata[i];
    }
    return(total_power);
  }

  template<class T, class U>
    diffractive_wavefront<T> & operator+=(diffractive_wavefront<T> & lhs_wf, const diffractive_wavefront<U> & rhs_wf) {

    if(lhs_wf.get_axes()!=rhs_wf.get_axes()){
      cerr << "operator+= error - mismatched arrays in diffractive_wavefronts\n";
      throw(string("diffractive_wavefront::operator+="));
    }

    lhs_wf.real_imag_conversion();
    rhs_wf.real_imag_conversion();
    if(lhs_wf.interleaved && !rhs_wf.interleaved) rhs_wf.interleaved_conversion();
    if(!lhs_wf.interleaved && rhs_wf.interleaved) rhs_wf.non_interleaved_conversion();

    int nelem = lhs_wf.diffractive_wavefront_header<T>::total_space();
    for(int i=0; i<2*nelem; i++)
      lhs_wf.wfdata[i] += rhs_wf.wfdata[i];
  
    return(lhs_wf);
  }

  template<class T, class U>
    diffractive_wavefront<T> & operator-=(diffractive_wavefront<T> & lhs_wf, const diffractive_wavefront<U> & rhs_wf) {

    if(lhs_wf.get_axes()!=rhs_wf.get_axes()){
      cerr << "operator-= error - mismatched arrays in diffractive_wavefronts\n";
      throw(string("diffractive_wavefront::operator-="));
    }

    lhs_wf.real_imag_conversion();
    rhs_wf.real_imag_conversion();
    if(lhs_wf.interleaved && !rhs_wf.interleaved) rhs_wf.interleaved_conversion();
    if(!lhs_wf.interleaved && rhs_wf.interleaved) rhs_wf.non_interleaved_conversion();

    int nelem = lhs_wf.diffractive_wavefront_header<T>::total_space();
    for(int i=0; i<2*nelem; i++) 
      lhs_wf.wfdata[i] -= rhs_wf.wfdata[i];

    return(lhs_wf);
  }

  template<class T, class U>
    diffractive_wavefront<T> & operator*=(diffractive_wavefront<T> & lhs_wf, const diffractive_wavefront<U> & rhs_wf) {

    if(lhs_wf.get_axes()!=rhs_wf.get_axes()){
      cerr << "operator*= error - mismatched arrays in diffractive_wavefronts\n";
      throw(string("diffractive_wavefront::operator*="));
    }

    lhs_wf.amp_phase_conversion();
    rhs_wf.amp_phase_conversion();
    if(lhs_wf.interleaved && !rhs_wf.interleaved) rhs_wf.interleaved_conversion();
    if(!lhs_wf.interleaved && rhs_wf.interleaved) rhs_wf.non_interleaved_conversion();

    // make sure that the diffractive_wavefronts are 
    // both interleaved or both not interleaved
    if(lhs_wf.interleaved) rhs_wf.interleaved_conversion();
    if(!lhs_wf.interleaved) rhs_wf.non_interleaved_conversion();

    int nelem = lhs_wf.diffractive_wavefront_header<T>::total_space();
    if(lhs_wf.interleaved){
      for(int i=0; i<nelem; i++){
	lhs_wf.wfdata[2*i] *= rhs_wf.wfdata[2*i];
	lhs_wf.wfdata[2*i+1] += rhs_wf.wfdata[2*i+1];
      }
    } else {
      for(int i=0; i<nelem; i++){
	lhs_wf.wfdata[i] *= rhs_wf.wfdata[i];
	lhs_wf.wfdata[i+nelem] += rhs_wf.wfdata[i+nelem];
      }
    }

    return(lhs_wf);
  }

  template<class T, class U>
    diffractive_wavefront<T> & operator/=(diffractive_wavefront<T> & lhs_wf, const diffractive_wavefront<U> & rhs_wf) {

    if(lhs_wf.get_axes()!=rhs_wf.get_axes()){
      cerr << "operator/= error - mismatched arrays in diffractive_wavefronts\n";
      throw(string("diffractive_wavefront::operator/="));
    }

    lhs_wf.amp_phase_conversion();
    rhs_wf.amp_phase_conversion();
    if(lhs_wf.interleaved) rhs_wf.interleaved_conversion();
    if(!lhs_wf.interleaved) rhs_wf.non_interleaved_conversion();

    int nelem = lhs_wf.diffractive_wavefront_header<T>::total_space();
    if(lhs_wf.interleaved){
      for(int i=0; i<nelem; i++){
	if(rhs_wf.wfdata[2*i] == 0) lhs_wf.wfdata[2*i] = 0;
	else lhs_wf.wfdata[2*i] /= rhs_wf.wfdata[2*i];
	lhs_wf.wfdata[2*i+1] += rhs_wf.wfdata[2*i+1];
      }
    } else {
      for(int i=0; i<nelem; i++){
	if(rhs_wf.wfdata[i] == 0) lhs_wf.wfdata[i] = 0;
	else lhs_wf.wfdata[i] /= rhs_wf.wfdata[i];
	lhs_wf.wfdata[i+nelem] -= rhs_wf.wfdata[i+nelem];
      }
    }

    return(lhs_wf);
  }

  template<class T>
    template<class U>
    diffractive_wavefront<T> & diffractive_wavefront<T>::operator+=(std::complex<U> c) {

    this->real_imag_conversion();

    int nelem = this->diffractive_wavefront_header<T>::total_space();
    double real = c.real();
    double imag = c.imag();
    if(interleaved){
      for(int i=0; i<nelem; i++){
	this->wfdata[2*i] += real;
	this->wfdata[2*i+1] += imag;
      }
    } else {
      for(int i=0; i<nelem; i++){
	this->wfdata[i] += real;
	this->wfdata[i+nelem] += imag;
      }
    }

    return(*this);
  }

  template<class T>
    template<class U>
    diffractive_wavefront<T> & diffractive_wavefront<T>::operator-=(std::complex<U> c) {

    this->real_imag_conversion();

    int nelem = this->diffractive_wavefront_header<T>::total_space();
    double real = c.real();
    double imag = c.imag();
    if(interleaved){
      for(int i=0; i<nelem; i++){
	this->wfdata[2*i] -= real;
	this->wfdata[2*i+1] -= imag;
      }
    } else {
      for(int i=0; i<nelem; i++){
	this->wfdata[i] -= real;
	this->wfdata[i+nelem] -= imag;
      }
    }

    return(*this);
  }

  template<class T>
    template<class U>
    diffractive_wavefront<T> & diffractive_wavefront<T>::operator*=(std::complex<U> c) {

    this->amp_phase_conversion();

    int nelem = this->diffractive_wavefront_header<T>::total_space();
    double abs_val = abs(c);
    double arg_val = arg(c);
    if(interleaved){
      for(int i=0; i<nelem; i++){
	this->wfdata[2*i] *= abs_val;
	this->wfdata[2*i+1] += arg_val;
      }
    } else {
      for(int i=0; i<nelem; i++){
	this->wfdata[i] *= abs_val;
	this->wfdata[i+nelem] += arg_val;
      }
    }

    return(*this);
  }

  template<class T>
    template<class U>
    diffractive_wavefront<T> & diffractive_wavefront<T>::operator/=(std::complex<U> c) {

    this->amp_phase_conversion();

    int nelem = this->diffractive_wavefront_header<T>::total_space();
    double abs = c.abs();
    double arg = c.arg();
    if(interleaved){
      for(int i=0; i<nelem; i++){
	if(abs==0) this->wfdata[2*i] = 0;
	else this->wfdata[2*i] /= abs;
	this->wfdata[2*i+1] -= arg;
      }
    } else {
      for(int i=0; i<nelem; i++){
	if(abs==0) this->wfdata[i] = 0;
	else this->wfdata[i] /= abs;
	this->wfdata[i+nelem] -= arg;
      }
    }

    return(*this);
  }

  template<class T, class U>
    diffractive_wavefront<T> operator+(const diffractive_wavefront<T> & wf1, const diffractive_wavefront<U> & wf2) {
    diffractive_wavefront<T> tmpwf(wf1);
    tmpwf+=wf2;
    return(tmpwf);
  }

  template<class T, class U>
    diffractive_wavefront<T> operator-(const diffractive_wavefront<T> & wf1, const diffractive_wavefront<U> & wf2) {
    diffractive_wavefront<T> tmpwf(wf1);
    tmpwf-=wf2;
    return(tmpwf);
  } 

  template<class T, class U>
    diffractive_wavefront<T> operator*(const diffractive_wavefront<T> & wf1, const diffractive_wavefront<U> & wf2) {
    diffractive_wavefront<T> tmpwf(wf1);
    tmpwf*=wf2;
    return(tmpwf);
  }

  template<class T, class U>
    diffractive_wavefront<T> operator/(const diffractive_wavefront<T> & wf1, const diffractive_wavefront<U> & wf2) {
    diffractive_wavefront<T> tmpwf(wf1);
    tmpwf/=wf2;
    return(tmpwf);
  }

  template<class T, class U>
    diffractive_wavefront<T> operator+(const diffractive_wavefront<T> & wf1, std::complex<U> c) {  
    diffractive_wavefront<T> tmpwf(wf1);
    tmpwf+=c;
    return(tmpwf);
  }

  template<class T, class U>
    diffractive_wavefront<T> operator-(const diffractive_wavefront<T> & wf1, std::complex<U> c) {  
    diffractive_wavefront<T> tmpwf(wf1);
    tmpwf-=c;
    return(tmpwf);
  }

  template<class T, class U>
    diffractive_wavefront<T> operator*(const diffractive_wavefront<T> & wf1, std::complex<U> c) {  
    diffractive_wavefront<T> tmpwf(wf1);
    tmpwf*=c;
    return(tmpwf);
  }

  template<class T, class U>
    diffractive_wavefront<T> operator/(const diffractive_wavefront<T> & wf1, std::complex<U> c) {  
    diffractive_wavefront<T> tmpwf(wf1);
    tmpwf/=c;
    return(tmpwf);
  }

  template<class T>
    bool operator !=(const diffractive_wavefront<T> & wf1, const diffractive_wavefront<T> & wf2) {
    return(!(wf1==wf2));
  }

}

#endif
