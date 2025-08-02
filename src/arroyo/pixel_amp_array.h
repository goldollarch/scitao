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

#ifndef PIXEL_AMP_ARRAY_H
#define PIXEL_AMP_ARRAY_H

#include <vector>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include "linear_algebra.h"
#include "fft_manager.h"
#include "pixel_array.h"

namespace Arroyo {

  using std::string;
  using std::vector;
  using std::ostream;

  /* forward declarations */
  template <class precision> class pixel_amp_array;
  template <class precision> bool operator == (const pixel_amp_array<precision> &p1,
					const pixel_amp_array<precision> &p2);
  template <class precision> bool operator != (const pixel_amp_array<precision> &p1,
					const pixel_amp_array<precision> &p2);

  ///
  ///  A class to represent pixel arrays in which entries should be 
  ///  interpreted as amplitudes
  ///

  template <class precision>
    class pixel_amp_array :
    public pixel_array<precision> {

    public:

    ///////////////////////////////////////////
    ///  Null constructor  
    pixel_amp_array() : pixel_array<precision>() {};

    ///////////////////////////////////////////
    ///  Copy constructor
    template<class U>
      pixel_amp_array(const pixel_amp_array<U> & pixamparr){
      pixel_amp_array::operator=(pixamparr);
    }

    ///////////////////////////////////////////
    ///  Copy from pixel_array
    template<class U>
      pixel_amp_array(const pixel_array<U> & pixarr){
      pixel_array<precision>::operator=(pixarr);
    }

    ///////////////////////////////////////////
    ///  Construct from iofits object
    pixel_amp_array(const iofits & iof) : pixel_array<precision>(iof) {};

    ///////////////////////////////////////////
    ///  Construct from arrays 
    pixel_amp_array(const std::vector<long> & in_axes, 
		    const precision * data = NULL, 
		    const float * wts = NULL) :
      pixel_array<precision>(in_axes, data, wts) {};

    ///////////////////////////////////////////
    ///  Construct an instance with pixel limits given by pixel_limits
    ///  These limits must be contained by pixarr's arrays
    template<class U>
      pixel_amp_array(const pixel_amp_array<U> & pixarr,
			const std::vector<long> & pixel_limits) 
        : pixel_array<precision>(pixarr, pixel_limits) {}

    ///////////////////////////////////////////
    ///  Construct from median of a vector of pixel_amp_arrays
    template <class U>
      pixel_amp_array(std::vector<pixel_amp_array<U> *> & paas);

    ///////////////////////////////////////////
    ///  Destructor
    ~pixel_amp_array(){};  

    ///////////////////////////////////////////
    ///  Operator =
    pixel_amp_array & operator = (const pixel_amp_array<precision> & pixamparr){
      if(this == &pixamparr)
	return(*this);
      pixel_array<precision>::operator=(pixamparr);
      return(*this);
    };

    ///////////////////////////////////////////
    ///  Find the rms of a pixel_array
    ///  The function sorts the array and discards 
    ///  discard_fraction of the points at either end
    ///  of the resulting array
    void mean_and_rms(double & mean, double & rms) const;

    ///////////////////////////////////////////
    ///  Sigma clip a pixel_amp_array
    ///  (i.e. set the weight to zero if value>sigma_clip * sigma
    void sigma_clip(const int & niter, const double & sigma_clip);

    ///////////////////////////////////////////
    ///  Amp clip the pixel_amp_array
    void amp_clip(const double & min, const double & max);

    ///////////////////////////////////////////
    ///  Function to interpolate over gaps in the 
    ///  array.  Only points surrounded on at least 
    ///  3 sides by points with nonzero weights are fixed.
    void interpolate();
    
    ///////////////////////////////////////////
    ///  Function to model this data as a single
    ///  shifted and scaled copy of the array
    ///  psf.  The best fitting model is returned
    ///  in fitted_psf, and the residuals are 
    ///  returned in fit_residual.  Relative
    ///  offsets are measured in pixels.
    /// 
    ///  Arrays must be the same dimension.
    double fit(const pixel_amp_array<precision> & psf, 
	       double & differential_amplitude,
	       vector<double> & relative_offsets,
	       pixel_amp_array<precision> & fitted_psf,
	       pixel_amp_array<precision> & fit_residual,
	       pixel_amp_array<precision> & orig) const;

    ///////////////////////////////////////////
    ///  Function to model this data as a single shifted and scaled
    ///  copy of the array psf, plus a DC offset.  The best fitting
    ///  model is returned in fitted_psf, and the residuals are
    ///  returned in fit_residual.  Relative offsets are measured in
    ///  pixels.
    /// 
    ///  Arrays must be the same dimension.
    double fit(const pixel_amp_array<precision> & psf, 
	       double & differential_amplitude,
	       double & offset,
	       bool fit_for_offset,
	       vector<double> & relative_offsets,
	       pixel_amp_array<precision> & fitted_psf,
	       pixel_amp_array<precision> & fit_residual,
	       pixel_amp_array<precision> & orig) const;

    ///////////////////////////////////////////
    ///  Function to model this data as npsfs
    ///  shifted and scaled copies of the array
    ///  psf.  The best fitting model is returned
    ///  in fitted_psf, and the residuals are 
    ///  returned in fit_residual.  Relative
    ///  offsets are measured in pixels.
    /// 
    ///  Arrays must be the same dimension.
    double fit(int npsfs,
	       const pixel_amp_array<precision> & psf, 
	       vector<double> & differential_amplitudes,
	       vector<vector<double> > & relative_offsets,
	       pixel_amp_array<precision> & fitted_model,
	       pixel_amp_array<precision> & fit_residual) const;

    ///////////////////////////////////////////
    ///  Function to perform aperture photometry
    ///  All arguments in pixels
    double aperture_photometry(double x, double y, double inner_radius,
						double outer_radius) const;

    ///////////////////////////////////////////
    /// This function compares this pixel_amp_array
    /// to the argument and deletes points in this instance for which the
    /// difference between the two instances is greater than sigma_clip*rms
    /// and for which the point in this instance is greater than sigma_clip*rms
    template<class U>
      void flag_outliers(const pixel_amp_array<U> & pixamparr,
						const double & threshold);

    ///////////////////////////////////////////
    ///  Function to perform an azimuthal average around a central point
    std::vector<float> azav(int xcenter, int ycenter);

    ///////////////////////////////////////////
    ///  Friend operator ==  for pixel_amp_array
    friend bool operator ==(const pixel_amp_array<precision> &p1,
				const pixel_amp_array<precision> &p2) {
      if(p1.get_axes()!=p2.get_axes()) return(0);
      for(int i=0; i<p1.axes.size(); i++)
	for(int j=0; j<p1.axes[i]; j++)
	  if(p1.pixeldata[i]!=p2.pixeldata[i]) return(0);
    
      if((p1.pixelwts!=NULL && p2.pixelwts==NULL) ||
	 (p1.pixelwts==NULL && p2.pixelwts!=NULL)){
	cerr << "pixel_amp_array<precision> & operator == error - "
	     << "weights not defined for both objects\n";
	throw(std::string("pixel_amp_array<precision> & operator =="));
      }
    
      if(p1.pixelwts!=NULL)
	for(int i=0; i<p1.axes.size(); i++)
	  for(int j=0; j<p1.axes[i]; j++)
	    if(p1.pixelwts[i]!=p2.pixelwts[i]) return(0);
    
      return(1);
    };
  };  


  template<class precision>
    template <class U>
    pixel_amp_array<precision>::pixel_amp_array(std::vector<pixel_amp_array<U> *> & paas){
  
    if(pixel_array<precision>::verbose_level) 
      cout << "pixel_amp_array::pixel_amp_array - forming array from median of " 
	   << paas.size() << " images\n";
  
    if(paas.size()==0){
      pixel_amp_array<U> pixarr;
      operator=(pixarr);
      return;
    }
    if(paas.size()==1 || paas.size()==2){
      operator=(*(paas[0]));
      return;
    }
  
    for(int i=0; i<paas.size(); i++){
      if(paas[i]->axes!=paas[0]->axes){
	cerr << "pixel_amp_array::pixel_amp_array error - "
	     << "pixel_amp_arrays have mismatched axes\n";  
	paas[i]->print_axes(cerr);
	paas[0]->print_axes(cerr);
	throw(std::string("pixel_amp_array::pixel_amp_array"));
      }
    }    
  
    int narrays = paas.size();
    int medelem = narrays/2;
    int nelem = paas[0]->total_space();
    std::vector<U> array(paas.size());
    pixel_array<precision> pa(paas[0]->axes, NULL, NULL);
    pixel_array<precision>::operator=(pa);
    for(int i=0; i<nelem; i++){
      for(int j=0; j<narrays; j++)
	array[j] = (precision)(paas[j]->data(i));
      sort(array.begin(),array.end());
      pixel_array<precision>::pixeldata[i] = array[medelem];
    }
  }

  template<class precision>
    void pixel_amp_array<precision>::mean_and_rms(double & mean, double & rms) const {
  
    int nelems = pixel_array<precision>::total_space();
    if(nelems < 2){
      cerr << "pixel_amp_array<precision>::mean_and_rms error - "
	   << "cannot compute rms with less than 2 points\n";
      throw(std::string("pixel_amp_array<precision>::mean_and_rms"));
    }
  
    mean = 0;
    double meansq = 0, wts = 0;
    if(pixel_array<precision>::pixelwts==NULL){
      for(int i=0; i<nelems; i++)
	mean += pixel_array<precision>::pixeldata[i];
      mean = mean/((double)nelems);
      for(int i=0; i<nelems; i++)
	meansq += pixel_array<precision>::pixeldata[i]*pixel_array<precision>::pixeldata[i] - mean*mean;
      meansq = meansq/((double)(nelems));
    } else {
      for(int i=0; i<nelems; i++){
	mean += pixel_array<precision>::pixeldata[i]*pixel_array<precision>::pixelwts[i];
	wts += pixel_array<precision>::pixelwts[i];
      }
      if(wts==0){
	cerr << "pixel_amp_array<precision>::mean_and_rms error - weights are all zero\n";
	throw(std::string("pixel_amp_array<precision>::mean_and_rms"));
      }
      mean /= wts;
      for(int i=0; i<nelems; i++)
	meansq += (pixel_array<precision>::pixeldata[i] - mean)*(pixel_array<precision>::pixeldata[i]-mean)*pixel_array<precision>::pixelwts[i];
      meansq /= wts;
    }
  
    if(meansq<0){
      cerr << "pixel_amp_array<precision>::mean_and_rms error - mean square < 0\n";
      cerr << "mean " << mean << "\twts " << wts << "\tmeansq " << meansq << endl;
      throw(std::string("pixel_amp_array<precision>::mean_and_rms"));
    }
  
    if(pixel_array<precision>::verbose_level)
      cerr << "pixel_amp_array<precision>::mean_and_rms - returning meansq "
           << meansq << " mean " 
	   << mean << "\trms " << sqrt(meansq) << endl;
    rms = sqrt(meansq);
  }

  template<class precision>
    void pixel_amp_array<precision>::sigma_clip(const int & niter,
    					const double & sigma_clip){

    if(pixel_array<precision>::pixelwts==NULL) pixel_array<precision>::allocate_weights(1);
  
    int bad_pixel_count = 0;
    double mean, rms=0, last_rms=-1, wtsq;
    int nelem = this->total_space();
  
    for(int i=0; i<niter; i++){
      mean_and_rms(mean, rms);
      if(rms == 0)
	cout << "pixel_amp_array<precision>::sigma_clip warning:"
	     << " rms = 0 - setting all weights to zero\n";
      if(pixel_array<precision>::verbose_level) 
	cout << "pixel_amp_array<precision>::sigma_clip - "
	     << "clipping all points with value less than " 
	     << mean + sigma_clip*rms << endl;

      if(last_rms==rms) break;

      for(int i=0; i<nelem; i++){
	if(pixel_array<precision>::pixelwts[i]!=0 && fabs(pixel_array<precision>::pixeldata[i]-mean) > sigma_clip*rms){
	  pixel_array<precision>::pixeldata[i] = 0;
	  pixel_array<precision>::pixelwts[i] = 0;
	  bad_pixel_count++;
	}
      }
      last_rms = rms;
    }
    if(pixel_array<precision>::verbose_level) {
      cerr << "pixel_amp_array<precision>::sigma_clip - clipped " 
	   << bad_pixel_count << " out of " << this->total_space() << " points ("
	   << ((double)(bad_pixel_count)/((double)this->total_space()*100.0)) << "%)\n";
      cerr << "pixel_amp_array<precision>::sigma_clip - mean " << mean
           << "\trms " << rms << endl;
    }
  }

  template<class precision>
    void pixel_amp_array<precision>::amp_clip(const double & min, const double & max){

    int nelem = this->total_space();
    if(this->pixelwts==NULL) this->allocate_weights(1);
    for(int i=0; i<nelem; i++){
      if(this->pixeldata[i]>max || this->pixeldata[i]<min){
	if(this->verbose_level) 
	  cout << "pixel_amp_array::amp_clip - clipping ("
	       << i/this->axes[0] << "," << i%this->axes[0] << ")\t"
	       << this->pixeldata[i] << endl;
	this->pixeldata[i] = 0;
	this->pixelwts[i] = 0;
      }
    }
  }


  template<class precision>
    void pixel_amp_array<precision>::interpolate(){
  
    int nelem = this->total_space();
    int nnonzero;
    double sum, wsum;
    if(!this->weights_allocated()){
      cerr << "pixel_amp_array<precision>::interpolate error - weights are not allocated\n";
      throw(std::string("pixel_amp_array<precision>::interpolate"));
    }
    for(int i=0; i<nelem; i++){
      if(this->pixelwts[i] == 0){
      
	nnonzero = 0;
	sum = wsum = 0;
	if(i-this->axes[0]>0 && this->pixelwts[i-this->axes[0]]!=0){
	  nnonzero++; 
	  sum += this->pixeldata[i-this->axes[0]]*this->pixelwts[i-this->axes[0]];
	  wsum += this->pixelwts[i-this->axes[0]];
	}
	if(i+this->axes[0]>0 && this->pixelwts[i+this->axes[0]]!=0){
	  nnonzero++; 
	  sum += this->pixeldata[i+this->axes[0]]*this->pixelwts[i+this->axes[0]];
	  wsum += this->pixelwts[i+this->axes[0]];
	}
	if(i%this->axes[0]>0 && this->pixelwts[i-1]!=0){
	  nnonzero++; 
	  sum += this->pixeldata[i-1]*this->pixelwts[i-1];
	  wsum += this->pixelwts[i-1];
	}
	if(i%this->axes[0]<this->axes[0]-1 && this->pixelwts[i+1]!=0){
	  nnonzero++; 
	  sum += this->pixeldata[i+1]*this->pixelwts[i+1];
	  wsum += this->pixelwts[i+1];
	}
      
	if(nnonzero>=3){
	  if(pixel_array<precision>::verbose_level)
	    cout << "pixel_amp_array<precision>::interpolate - correcting point " 
		 << i/this->axes[0] << ", " << i%this->axes[0] << endl;
	  this->pixeldata[i] = sum/wsum;
	  this->pixelwts[i] = wsum;
	}
      }
    }
  }
    
  template<class T>
    double pixel_amp_array<T>::fit(const pixel_amp_array<T> & psf, 
				   double & differential_amplitude,
				   vector<double> & relative_offsets,
				   pixel_amp_array<T> & fitted_psf,
				   pixel_amp_array<T> & fit_residual,
				   pixel_amp_array<T> & orig) const {
    
    if(psf.axes.size()!=2 || this->axes!=psf.axes){
      cerr << "pixel_amp_array::fit error - array mismatch\n";
      cerr << "this axes " << this->axes[0] << " x " << this->axes[1] << endl;
      cerr << "psf axes " << psf.get_axes()[0] << " x " << psf.get_axes()[1] << endl;
      throw(string("pixel_amp_array::fit"));
    }
    if(this->axes[0]<=1 || this->axes[1]<=1 || psf.axes[0]<=1 || psf.axes[1]<=1){
      cerr << "pixel_amp_array::fit error - "
	   << "array doesn't contain enough elements\n";
      throw(string("pixel_amp_array::fit")); 
    } 

    if(pixel_array<T>::verbose_level)
      cout << "pixel_amp_array::fit - seeking offset between arrays of dimension \n"
	   << this->axes[0] << "x" << this->axes[1]  
	   << endl;

    pixel_array<T> xcorr_psf = this->cross_correlate(psf);
    vector<int> xcorr_minpixel(2,0), xcorr_maxpixel(2,0);

    double min, max;
    xcorr_psf.min_and_max(min, xcorr_minpixel, max, xcorr_maxpixel);

    if(pixel_array<T>::verbose_level) {
      cout << "pixel_array<T>::fit - xcorr max " 
	   << xcorr_maxpixel[0] << "\t" 
	   << xcorr_maxpixel[1] 
	   << "\t" << max 
	   << endl;
    }

    double xslope = this->axes[1]%2==1 ? 0 : M_PI/2.;
    double yslope = this->axes[0]%2==1 ? 0 : M_PI/2.;

    cerr << "xslope " << xslope << " yslope " << yslope << endl;

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


    if(pixel_array<T>::verbose_level)
      cerr << "pixel_amp_array::fit - allocating memory\n";
    int nelem = this->axes[0]*this->axes[1];
    T * this_arr;
    T * psf_arr;
    T * diff_arr;
    try{
      this_arr = new T[2*nelem];
      psf_arr = new T[2*nelem];
      diff_arr = new T[2*nelem];
    } catch(...) {
      cerr << "pixel_amp_array::fit error - error allocating memory\n";
      throw(string("pixel_amp_array::fit"));
    }

    if(pixel_array<T>::verbose_level)
      cerr << "pixel_amp_array::fit - filling in arrays\n";
    int index;
    double amp, phase;
    for(int i=-this->axes[1]/2; i<this->axes[1]/2+x_extrapix; i++){
      for(int j=-this->axes[0]/2; j<this->axes[0]/2+y_extrapix; j++){
	index = (i+this->axes[1]/2)*this->axes[0]+j+this->axes[0]/2;

	phase = -xslope*(i+x_halfpix) - yslope*(j+y_halfpix);

	amp = this->pixeldata[index];
	this_arr[2*index] = amp*cos(phase);
	this_arr[2*index+1] = amp*sin(phase);

	amp = psf.pixeldata[index];
	psf_arr[2*index] = amp*cos(phase);
	psf_arr[2*index+1] = amp*sin(phase);
      }
    }

    if(pixel_array<T>::verbose_level)
      cerr << "pixel_amp_array::fit - performing transform\n";
    Arroyo::fft_manager<T> fft_mgr;  

    vector<long> flipped_axes(2, this->axes[0]);
    flipped_axes[0] = this->axes[1];
    fft_mgr.forward_fft(flipped_axes, false, true, this_arr);
    fft_mgr.forward_fft(flipped_axes, false, true, psf_arr);

    Arroyo::complex_cyclic_permutation(this->axes, 
				       this->axes[0]/2,
				       this->axes[1]/2, 
				       this_arr);

    Arroyo::complex_cyclic_permutation(psf.get_axes(), 
				       psf.get_axes()[0]/2,
				       psf.get_axes()[1]/2, 
				       psf_arr);


    if(pixel_array<T>::verbose_level)
      cerr << "pixel_amp_array::fit - allocating memory\n";
    T *real_product, *imag_product;
    try{
      real_product = new T[nelem];
      imag_product = new T[nelem];
    } catch(...) {
      cerr << "pixel_amp_array::fit error - error allocating memory\n";
      throw(string("pixel_amp_array::fit"));
    }

    if(pixel_array<T>::verbose_level)
      cerr << "pixel_amp_array::fit - computing intermediate arrays\n";
    double this_total_power=0, psf_total_power=0, cross_power=0;
    for(int i=-this->axes[1]/2; i<this->axes[1]/2+x_extrapix; i++){
      for(int j=-this->axes[0]/2; j<this->axes[0]/2+y_extrapix; j++){
	
	index = (i+this->axes[1]/2)*this->axes[0]+j+this->axes[0]/2;

	this_total_power += 
	  this_arr[2*index]*this_arr[2*index] + this_arr[2*index+1]*this_arr[2*index+1];
	psf_total_power += 
	  psf_arr[2*index]*psf_arr[2*index] + psf_arr[2*index+1]*psf_arr[2*index+1];
	
	real_product[index] = 
	  this_arr[2*index]*psf_arr[2*index] + this_arr[2*index+1]*psf_arr[2*index+1];
	imag_product[index] = 
	  this_arr[2*index+1]*psf_arr[2*index] - this_arr[2*index]*psf_arr[2*index+1];
	
	cross_power += real_product[index];
      }
    }

    if(pixel_array<T>::verbose_level)
      cout << "This total power "
	   << this_total_power 
	   << " PSF total power "
	   << psf_total_power
	   << " cross power "
	   << cross_power
	   << endl;

    // Perform a search for the best fitting fractional offset.
    double faca = 2*M_PI/(double)(this->axes[0]);
    double facb = 2*M_PI/(double)(this->axes[1]);
    double dcross_du, dcross_dv;
    double dsqcross_dusq, dsqcross_dvsq, dsqcross_dudv;
    //double x_offset=-xcorr_maxpixel[0], y_offset=-xcorr_maxpixel[1];
    double x_offset=0, y_offset=0;
    double delta_x_offset, delta_y_offset;
    double real_part_of_cross_term;
    double deriv_real_part_of_cross_term;
    double sqderiv_real_part_of_cross_term;
    double sum_real_part_of_cross_term;
    double chisq, last_chisq=DBL_MAX;
    while(1){
      dcross_du = dcross_dv = dsqcross_dusq = dsqcross_dvsq = dsqcross_dudv = 0;
      sum_real_part_of_cross_term = 0;
      for(int i=-this->axes[1]/2; i<this->axes[1]/2+x_extrapix; i++){
	for(int j=-this->axes[0]/2; j<this->axes[0]/2+y_extrapix; j++){

	  index = (i+this->axes[1]/2)*this->axes[0]+j+this->axes[0]/2;

	  real_part_of_cross_term = 
	    real_product[index] * cos(i*y_offset*facb + j*x_offset*faca) +
	    imag_product[index] * sin(i*y_offset*facb + j*x_offset*faca);

	  sum_real_part_of_cross_term += 
	    real_part_of_cross_term;

	  deriv_real_part_of_cross_term = 
	    -1 * real_product[index] * sin(i*y_offset*facb + j*x_offset*faca) +
	    imag_product[index] * cos(i*y_offset*facb + j*x_offset*faca);

	  dcross_du += i*facb*deriv_real_part_of_cross_term;
	  dcross_dv += j*faca*deriv_real_part_of_cross_term;

	  sqderiv_real_part_of_cross_term = 
	    -1 * real_product[index] * cos(i*y_offset*facb + j*x_offset*faca) -
	    imag_product[index] * sin(i*y_offset*facb + j*x_offset*faca);
	  
	  dsqcross_dusq += i*i*facb*facb*sqderiv_real_part_of_cross_term;
	  dsqcross_dudv += i*j*facb*faca*sqderiv_real_part_of_cross_term;
	  dsqcross_dvsq += j*j*faca*faca*sqderiv_real_part_of_cross_term;
	}
      }

      differential_amplitude = fabs(cross_power / psf_total_power);

      chisq = .5*(this_total_power + 
		  differential_amplitude*differential_amplitude*psf_total_power - 
		  2*differential_amplitude*sum_real_part_of_cross_term);

      //////////////////////
      // Compute chi squared
      double angle, real, imag, real_arg=0, imag_arg=0, xchisq=0, xxchisq=0;
      for(int u=-this->axes[1]/2; u<this->axes[1]/2+x_extrapix; u++){
	for(int v=-this->axes[0]/2; v<this->axes[0]/2+y_extrapix; v++){
	  index = (u+this->axes[1]/2)*this->axes[0]+v+this->axes[0]/2;
		
	  real = imag = 0;
	  angle = -(y_offset*facb*u + x_offset*faca*v);
	  real += differential_amplitude*cos(angle);
	  imag += differential_amplitude*sin(angle);
	
	  real_arg = this_arr[2*index] - psf_arr[2*index]*real - psf_arr[2*index+1]*imag;
	  imag_arg = this_arr[2*index+1] - psf_arr[2*index+1]*real + psf_arr[2*index]*imag;

	  xchisq += .5*(real_arg*real_arg + imag_arg*imag_arg);

	  xxchisq += 
	    real_product[index] * cos(angle) -
	    imag_product[index] * sin(angle);

	}
      }
      double tmp_xxchisq = .5*(this_total_power + 
			       differential_amplitude*differential_amplitude*psf_total_power - 
			       2*differential_amplitude*xxchisq);
      ///////////////////////////

      if(pixel_array<T>::verbose_level)
	cout << " chisq " << chisq
	     << " xchisq " << xchisq
	     << " xxchisq " << tmp_xxchisq
	     << " amp " << differential_amplitude
	     << "\t";


      // Find the zero in two dimensions
      double denom = dsqcross_dusq*dsqcross_dvsq - dsqcross_dudv*dsqcross_dudv;
      delta_y_offset = (dcross_du*dsqcross_dvsq - dcross_dv*dsqcross_dudv)/denom;
      delta_x_offset = (-dcross_du*dsqcross_dudv + dcross_dv*dsqcross_dusq)/denom;

      if(pixel_array<T>::verbose_level)
	cerr << "dcross_du " << dcross_du 
	     << " dcross_dv " << dcross_dv
	     << " dsqcross_dusq " << dsqcross_dusq
	     << " dsqcross_dudv " << dsqcross_dudv
	     << " dsqcross_dvsq " << dsqcross_dvsq
	     << " delta_x " << delta_x_offset
	     << " delta_y " << delta_y_offset
	     << endl;

      x_offset-=delta_x_offset;
      y_offset-=delta_y_offset;

      differential_amplitude = sum_real_part_of_cross_term / psf_total_power;
      
      chisq = .5*(this_total_power + 
		  differential_amplitude*differential_amplitude*psf_total_power - 
		  2*differential_amplitude*sum_real_part_of_cross_term);

      //////////////////////
      // Compute chi squared
      real_arg = imag_arg = xchisq = 0;
      double ychisq = 0;
      for(int u=-this->axes[1]/2; u<this->axes[1]/2+x_extrapix; u++){
	for(int v=-this->axes[0]/2; v<this->axes[0]/2+y_extrapix; v++){
	  index = (u+this->axes[1]/2)*this->axes[0]+v+this->axes[0]/2;
		
	  real = imag = 0;
	  angle = -(y_offset*facb*u + x_offset*faca*v);
	  real += differential_amplitude*cos(angle);
	  imag += differential_amplitude*sin(angle);
	  
	  real_arg = this_arr[2*index] - psf_arr[2*index]*real - psf_arr[2*index+1]*imag;
	  imag_arg = this_arr[2*index+1] - psf_arr[2*index+1]*real + psf_arr[2*index]*imag;

	  ychisq += .5*(real_arg*real_arg + imag_arg*imag_arg);

	}
      }

      ///////////////////////////

      if(pixel_array<T>::verbose_level)
	cout << " chisq " << chisq
	     << " xchisq " << ychisq
	     << " xt " << sum_real_part_of_cross_term
	     << " amp " << differential_amplitude
	     << " off " << x_offset 
	     << " " << delta_x_offset 
	     << " " << y_offset 
	     << " " << delta_y_offset 
	     << endl;

      /*
      if(last_chisq<chisq)
	break;
      else 
	last_chisq = chisq;
      */
      if(fabs(delta_x_offset)<1e-6 &&
	 fabs(delta_y_offset)<1e-6)
	break;

    }

    relative_offsets.resize(2);
    relative_offsets[0] = x_offset;
    relative_offsets[1] = y_offset;

    // Transform back to get residual from fit
    if(pixel_array<T>::verbose_level)
      cerr << "pixel_amp_array::fit - computing residuals\n";

    double this_real, this_imag;
    for(int i=-this->axes[1]/2; i<this->axes[1]/2+x_extrapix; i++){
      for(int j=-this->axes[0]/2; j<this->axes[0]/2+y_extrapix; j++){
	index = (i+this->axes[1]/2)*this->axes[0]+j+this->axes[0]/2;

	amp = 
	  differential_amplitude*
	  sqrt(psf_arr[2*index]*psf_arr[2*index] + 
	       psf_arr[2*index+1]*psf_arr[2*index+1]);

	/*
	phase = 
	  atan2(psf_arr[2*index+1],psf_arr[2*index]) - 
	  i*y_offset*facb -
	  j*x_offset*faca - 
	  xslope*(i+x_halfpix) -
	  yslope*(j+y_halfpix);
	*/

	phase = 
	  atan2(psf_arr[2*index+1],psf_arr[2*index]) + 
	  i*y_offset*facb + j*x_offset*faca;


	psf_arr[2*index] = amp*cos(phase);
	psf_arr[2*index+1] = amp*sin(phase);

	diff_arr[2*index] = this_arr[2*index] - amp*cos(phase);
	diff_arr[2*index+1] = this_arr[2*index+1] - amp*sin(phase);

      }
    }


    Arroyo::complex_cyclic_permutation(this->get_axes(), 
				       -this->get_axes()[0]/2,
				       -this->get_axes()[1]/2, 
				       psf_arr);
    
    Arroyo::complex_cyclic_permutation(this->get_axes(), 
				       -this->get_axes()[0]/2,
				       -this->get_axes()[1]/2, 
				       this_arr);

    Arroyo::complex_cyclic_permutation(this->get_axes(), 
				       -this->get_axes()[0]/2,
				       -this->get_axes()[1]/2, 
				       diff_arr);

    fft_mgr.backward_fft(flipped_axes, false, true, this_arr);
    fft_mgr.backward_fft(flipped_axes, false, true, psf_arr);
    fft_mgr.backward_fft(flipped_axes, false, true, diff_arr);

    double norm_fac = 1/(double)(this->total_space());
    //double norm_fac = 1;
    for(int i=0; i<nelem; i++){
      psf_arr[i] = norm_fac*psf_arr[2*i];
      this_arr[i] = norm_fac*this_arr[2*i];

      /*
      psf_arr[i] = norm_fac*
	sqrt(psf_arr[2*i]*psf_arr[2*i]+
	     psf_arr[2*i+1]*psf_arr[2*i+1]);
      this_arr[i] = norm_fac*
	sqrt(this_arr[2*i]*this_arr[2*i]+
	     this_arr[2*i+1]*this_arr[2*i+1]);
      */
    }

    fitted_psf = pixel_array<T>(this->axes, psf_arr);
    orig = pixel_array<T>(this->axes, this_arr);

    fit_residual = orig - fitted_psf;

    // Clean up
    if(pixel_array<T>::verbose_level)
      cerr << "pixel_amp_array::fit - cleaning up\n";

    delete [] this_arr;
    delete [] psf_arr;
    delete [] diff_arr;
    delete [] real_product;
    delete [] imag_product;

    return(chisq);
  }


  template<class T>
    double pixel_amp_array<T>::fit(const pixel_amp_array<T> & psf, 
				   double & differential_amplitude,
				   double & offset,
				   bool fit_for_offset,
				   vector<double> & relative_offsets,
				   pixel_amp_array<T> & fitted_psf,
				   pixel_amp_array<T> & fit_residual,
				   pixel_amp_array<T> & orig) const {
    
    if(psf.axes.size()!=2 || this->axes!=psf.axes){
      cerr << "pixel_amp_array::fit error - array mismatch\n";
      cerr << "this axes " << this->axes[0] << " x " << this->axes[1] << endl;
      cerr << "psf axes " << psf.get_axes()[0] << " x " << psf.get_axes()[1] << endl;
      throw(string("pixel_amp_array::fit"));
    }
    if(this->axes[0]<=1 || this->axes[1]<=1 || psf.axes[0]<=1 || psf.axes[1]<=1){
      cerr << "pixel_amp_array::fit error - "
	   << "array doesn't contain enough elements\n";
      throw(string("pixel_amp_array::fit")); 
    } 

    if(pixel_array<T>::verbose_level)
      cout << "pixel_amp_array::fit - seeking offset between arrays of dimension \n"
	   << this->axes[0] << "x" << this->axes[1]  
	   << endl;

    pixel_array<T> xcorr_psf = this->cross_correlate(psf);
    vector<int> xcorr_minpixel(2,0), xcorr_maxpixel(2,0);

    double min, max;
    xcorr_psf.min_and_max(min, xcorr_minpixel, max, xcorr_maxpixel);

    if(pixel_array<T>::verbose_level) {
      cout << "pixel_array<T>::fit - xcorr max " 
	   << xcorr_maxpixel[0] << "\t" 
	   << xcorr_maxpixel[1] 
	   << "\t" << max 
	   << endl;
    }

    double xslope = this->axes[1]%2==1 ? 0 : M_PI/2.;
    double yslope = this->axes[0]%2==1 ? 0 : M_PI/2.;

    
    if(pixel_array<T>::verbose_level)
      cerr << "pixel_amp_array::fit - xslope " << xslope << " yslope " << yslope << endl;

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


    if(pixel_array<T>::verbose_level)
      cerr << "pixel_amp_array::fit - allocating memory\n";
    int nelem = this->axes[0]*this->axes[1];
    T * this_arr;
    T * psf_arr;
    try{
      this_arr = new T[2*nelem];
      psf_arr = new T[2*nelem];
    } catch(...) {
      cerr << "pixel_amp_array::fit error - error allocating memory\n";
      throw(string("pixel_amp_array::fit"));
    }

    if(pixel_array<T>::verbose_level)
      cerr << "pixel_amp_array::fit - filling in arrays\n";
    int index;
    double amp, phase;
    for(int i=-this->axes[1]/2; i<this->axes[1]/2+x_extrapix; i++){
      for(int j=-this->axes[0]/2; j<this->axes[0]/2+y_extrapix; j++){
	index = (i+this->axes[1]/2)*this->axes[0]+j+this->axes[0]/2;

	phase = -xslope*(i+x_halfpix) - yslope*(j+y_halfpix);

	amp = this->pixeldata[index];
	this_arr[2*index] = amp*cos(phase);
	this_arr[2*index+1] = amp*sin(phase);

	amp = psf.pixeldata[index];
	psf_arr[2*index] = amp*cos(phase);
	psf_arr[2*index+1] = amp*sin(phase);
      }
    }

    if(pixel_array<T>::verbose_level)
      cerr << "pixel_amp_array::fit - performing transform\n";
    Arroyo::fft_manager<T> fft_mgr;  

    vector<long> flipped_axes(2, this->axes[0]);
    flipped_axes[0] = this->axes[1];
    fft_mgr.forward_fft(flipped_axes, false, true, this_arr);
    fft_mgr.forward_fft(flipped_axes, false, true, psf_arr);

    double this_dc = this_arr[0];
    double psf_dc = psf_arr[0];
    double this_dc_0 = this_arr[0];
    double this_dc_1 = this_arr[1];
    double psf_dc_0 = psf_arr[0];
    double psf_dc_1 = psf_arr[1];

    Arroyo::complex_cyclic_permutation(this->axes, 
				       this->axes[0]/2,
				       this->axes[1]/2, 
				       this_arr);

    Arroyo::complex_cyclic_permutation(psf.get_axes(), 
				       psf.get_axes()[0]/2,
				       psf.get_axes()[1]/2, 
				       psf_arr);


    if(pixel_array<T>::verbose_level)
      cerr << "pixel_amp_array::fit - allocating memory\n";
    T *real_product, *imag_product;
    try{
      real_product = new T[nelem];
      imag_product = new T[nelem];
    } catch(...) {
      cerr << "pixel_amp_array::fit error - error allocating memory\n";
      throw(string("pixel_amp_array::fit"));
    }

    if(pixel_array<T>::verbose_level)
      cerr << "pixel_amp_array::fit - computing intermediate arrays\n";
    double this_total_power=0, psf_total_power=0, cross_power=0;
    for(int i=-this->axes[1]/2; i<this->axes[1]/2+x_extrapix; i++){
      for(int j=-this->axes[0]/2; j<this->axes[0]/2+y_extrapix; j++){
	
	index = (i+this->axes[1]/2)*this->axes[0]+j+this->axes[0]/2;

	this_total_power += 
	  this_arr[2*index]*this_arr[2*index] + this_arr[2*index+1]*this_arr[2*index+1];
	psf_total_power += 
	  psf_arr[2*index]*psf_arr[2*index] + psf_arr[2*index+1]*psf_arr[2*index+1];
	
	real_product[index] = 
	  this_arr[2*index]*psf_arr[2*index] + this_arr[2*index+1]*psf_arr[2*index+1];
	imag_product[index] = 
	  this_arr[2*index+1]*psf_arr[2*index] - this_arr[2*index]*psf_arr[2*index+1];
	
	cross_power += real_product[index];
      }
    }

    if(pixel_array<T>::verbose_level)
      cout << "This total power "
	   << this_total_power 
	   << " PSF total power "
	   << psf_total_power
	   << " cross power "
	   << cross_power << endl
	   << " this dc " << this_arr[nelem]
	   << "\t" << this_arr[nelem+1] 
	   << "\t" << this_dc_0 
	   << "\t" << this_dc_1 
	   << endl
	   << " psf dc " << psf_arr[nelem] 
	   << "\t" << psf_arr[nelem+1] 
	   << "\t" << psf_dc_0 
	   << "\t" << psf_dc_1 
	   << endl;

    // Perform a search for the best fitting fractional offset.
    double faca = 2*M_PI/(double)(this->axes[0]);
    double facb = 2*M_PI/(double)(this->axes[1]);
    double dcross_du, dcross_dv;
    double dsqcross_dusq, dsqcross_dvsq, dsqcross_dudv;
    double x_offset=0, y_offset=0;
    double delta_x_offset, delta_y_offset;
    double real_part_of_cross_term;
    double deriv_real_part_of_cross_term;
    double sqderiv_real_part_of_cross_term;
    double sum_real_part_of_cross_term;
    double chisq, last_chisq=DBL_MAX;
    while(1){
      dcross_du = dcross_dv = dsqcross_dusq = dsqcross_dvsq = dsqcross_dudv = 0;
      sum_real_part_of_cross_term = 0;
      for(int i=-this->axes[1]/2; i<this->axes[1]/2+x_extrapix; i++){
	for(int j=-this->axes[0]/2; j<this->axes[0]/2+y_extrapix; j++){

	  index = (i+this->axes[1]/2)*this->axes[0]+j+this->axes[0]/2;

	  real_part_of_cross_term = 
	    real_product[index] * cos(i*y_offset*facb + j*x_offset*faca) +
	    imag_product[index] * sin(i*y_offset*facb + j*x_offset*faca);

	  sum_real_part_of_cross_term += 
	    real_part_of_cross_term;

	  deriv_real_part_of_cross_term = 
	    -1 * real_product[index] * sin(i*y_offset*facb + j*x_offset*faca) +
	    imag_product[index] * cos(i*y_offset*facb + j*x_offset*faca);

	  dcross_du += i*facb*deriv_real_part_of_cross_term;
	  dcross_dv += j*faca*deriv_real_part_of_cross_term;

	  sqderiv_real_part_of_cross_term = 
	    -1 * real_product[index] * cos(i*y_offset*facb + j*x_offset*faca) -
	    imag_product[index] * sin(i*y_offset*facb + j*x_offset*faca);
	  
	  dsqcross_dusq += i*i*facb*facb*sqderiv_real_part_of_cross_term;
	  dsqcross_dudv += i*j*facb*faca*sqderiv_real_part_of_cross_term;
	  dsqcross_dvsq += j*j*faca*faca*sqderiv_real_part_of_cross_term;
	}
      }

      // Find the zero in two dimensions
      double denom = dsqcross_dusq*dsqcross_dvsq - dsqcross_dudv*dsqcross_dudv;
      delta_y_offset = (dcross_du*dsqcross_dvsq - dcross_dv*dsqcross_dudv)/denom;
      delta_x_offset = (-dcross_du*dsqcross_dudv + dcross_dv*dsqcross_dusq)/denom;
      x_offset-=delta_x_offset;
      y_offset-=delta_y_offset;

      if(pixel_array<T>::verbose_level)
	cerr << "dcross_du " << dcross_du 
	     << " dcross_dv " << dcross_dv
	     << " dsqcross_dusq " << dsqcross_dusq
	     << " dsqcross_dudv " << dsqcross_dudv
	     << " dsqcross_dvsq " << dsqcross_dvsq
	     << " delta_x " << delta_x_offset
	     << " delta_y " << delta_y_offset
	     << endl;

      if(!fit_for_offset){
	differential_amplitude = sum_real_part_of_cross_term / psf_total_power;
	offset = 0;
      } else {
	denom = psf_total_power - psf_dc*psf_dc;
	differential_amplitude = (sum_real_part_of_cross_term - psf_dc*this_dc)/denom;
	offset = (-sum_real_part_of_cross_term*psf_dc + psf_total_power*this_dc)/denom;
      }
      chisq = .5*(this_total_power + 
		  differential_amplitude*differential_amplitude*psf_total_power - 
		  2*differential_amplitude*sum_real_part_of_cross_term +
		  2*differential_amplitude*offset*psf_dc -
		  2*offset*this_dc +
		  offset*offset);

      //////////////////////
      // Compute chi squared
      double angle, real, imag, real_arg=0, imag_arg=0, ychisq=0;
      for(int u=-this->axes[1]/2; u<this->axes[1]/2+x_extrapix; u++){
	for(int v=-this->axes[0]/2; v<this->axes[0]/2+y_extrapix; v++){
	  index = (u+this->axes[1]/2)*this->axes[0]+v+this->axes[0]/2;
		
	  real = imag = 0;
	  angle = -(y_offset*facb*u + x_offset*faca*v);
	  real += differential_amplitude*cos(angle);
	  imag += differential_amplitude*sin(angle);
	  
	  real_arg = this_arr[2*index] - psf_arr[2*index]*real - psf_arr[2*index+1]*imag;
	  imag_arg = this_arr[2*index+1] - psf_arr[2*index+1]*real + psf_arr[2*index]*imag;

	  ychisq += .5*(real_arg*real_arg + imag_arg*imag_arg);

	}
      }

      ///////////////////////////

      if(pixel_array<T>::verbose_level){
	cout << " chisq " << chisq
	     << " ychisq " << ychisq
	     << " xt " << sum_real_part_of_cross_term
	     << " amp " << differential_amplitude;
	if(fit_for_offset)
	  cout << " offset " << offset;
	
	cout << " shift " << x_offset 
	     << " " << delta_x_offset 
	     << " " << y_offset 
	     << " " << delta_y_offset 
	     << endl;
      }

      /*
      if(last_chisq<chisq)
	break;
      else 
	last_chisq = chisq;
      */
      if(fabs(delta_x_offset)<1e-6 &&
	 fabs(delta_y_offset)<1e-6)
	break;

    }

    relative_offsets.resize(2);
    relative_offsets[0] = x_offset;
    relative_offsets[1] = y_offset;

    // Transform back to get residual from fit
    if(pixel_array<T>::verbose_level)
      cerr << "pixel_amp_array::fit - computing residuals\n";

    double this_real, this_imag;
    for(int i=-this->axes[1]/2; i<this->axes[1]/2+x_extrapix; i++){
      for(int j=-this->axes[0]/2; j<this->axes[0]/2+y_extrapix; j++){
	index = (i+this->axes[1]/2)*this->axes[0]+j+this->axes[0]/2;

	amp = 
	  differential_amplitude*
	  sqrt(psf_arr[2*index]*psf_arr[2*index] + 
	       psf_arr[2*index+1]*psf_arr[2*index+1]);

	/*
	phase = 
	  atan2(psf_arr[2*index+1],psf_arr[2*index]) - 
	  i*y_offset*facb -
	  j*x_offset*faca - 
	  xslope*(i+x_halfpix) -
	  yslope*(j+y_halfpix);
	*/

	phase = 
	  atan2(psf_arr[2*index+1],psf_arr[2*index]) + 
	  i*y_offset*facb + j*x_offset*faca;


	psf_arr[2*index] = amp*cos(phase);
	psf_arr[2*index+1] = amp*sin(phase);

      }
    }


    Arroyo::complex_cyclic_permutation(this->get_axes(), 
				       -this->get_axes()[0]/2,
				       -this->get_axes()[1]/2, 
				       psf_arr);
    
    Arroyo::complex_cyclic_permutation(this->get_axes(), 
				       -this->get_axes()[0]/2,
				       -this->get_axes()[1]/2, 
				       this_arr);

    psf_arr[0] += offset;

    fft_mgr.backward_fft(flipped_axes, false, true, this_arr);
    fft_mgr.backward_fft(flipped_axes, false, true, psf_arr);

    double norm_fac = 1/(double)(this->total_space());
    //double norm_fac = 1;
    for(int i=0; i<nelem; i++){
      psf_arr[i] = norm_fac*psf_arr[2*i];
      this_arr[i] = norm_fac*this_arr[2*i];

      /*
      psf_arr[i] = norm_fac*
	sqrt(psf_arr[2*i]*psf_arr[2*i]+
	     psf_arr[2*i+1]*psf_arr[2*i+1]);
      this_arr[i] = norm_fac*
	sqrt(this_arr[2*i]*this_arr[2*i]+
	     this_arr[2*i+1]*this_arr[2*i+1]);
      */
    }

    fitted_psf = pixel_array<T>(this->axes, psf_arr);
    orig = pixel_array<T>(this->axes, this_arr);

    fit_residual = orig - fitted_psf;

    double mean, rms;
    fit_residual.mean_and_rms(mean, rms);

    // Clean up
    if(pixel_array<T>::verbose_level)
      cerr << "pixel_amp_array::fit - cleaning up\n";

    delete [] this_arr;
    delete [] psf_arr;
    delete [] real_product;
    delete [] imag_product;

    return(chisq);
  }

  template<class T>
    double pixel_amp_array<T>::fit(int npsfs,
				   const pixel_amp_array<T> & psf, 
				   vector<double> & differential_amplitudes,
				   vector<vector<double> > & relative_offsets,
				   pixel_amp_array<T> & fitted_model,
				   pixel_amp_array<T> & fit_residual) const {
    

    
    if(npsfs<=0){
      cerr << "pixel_amp_array::fit error - cannot perform a fit to "
	   << npsfs
	   << " PSF's\n";
      throw(string("pixel_amp_array::fit"));
    }

    if(psf.axes.size()!=2 || this->axes!=psf.axes){
      cerr << "pixel_amp_array::fit error - array mismatch\n";
      cerr << "this axes " << this->axes[0] << " x " << this->axes[1] << endl;
      cerr << "psf axes " << psf.get_axes()[0] << " x " << psf.get_axes()[1] << endl;
      throw(string("pixel_amp_array::fit"));
    }
    if(this->axes[0]<=1 || this->axes[1]<=1 || psf.axes[0]<=1 || psf.axes[1]<=1){
      cerr << "pixel_amp_array::fit error - "
	   << "array doesn't contain enough elements\n";
      throw(string("pixel_amp_array::fit")); 
    } 


    if(pixel_array<T>::verbose_level)
      cout << "pixel_amp_array::fit - seeking offset between arrays of dimension \n"
	   << this->axes[0] << "x" << this->axes[1]  
	   << endl;

    pixel_array<T> xcorr_psf = this->cross_correlate(psf);
    vector<int> xcorr_minpixel(2,0), xcorr_maxpixel(2,0);

    double min, max;
    xcorr_psf.min_and_max(min, xcorr_minpixel, max, xcorr_maxpixel);

    if(pixel_array<T>::verbose_level) {
      cout << "pixel_array<T>::fit - xcorr max " 
	   << xcorr_maxpixel[0] << "\t" 
	   << xcorr_maxpixel[1] 
	   << "\t" << max 
	   << endl;
    }

    double xslope = this->axes[1]%2==1 ? 0 : M_PI/2.;
    double yslope = this->axes[0]%2==1 ? 0 : M_PI/2.;

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


    if(pixel_array<T>::verbose_level)
      cerr << "pixel_amp_array::fit - allocating memory\n";
    int nelem = this->total_space();
    T * this_arr;
    T * psf_arr;
    T * fit_arr;
    try{
      this_arr = new T[2*nelem];
      psf_arr = new T[2*nelem];
      fit_arr = new T[2*nelem];
    } catch(...) {
      cerr << "pixel_amp_array::fit error - error allocating memory\n";
      throw(string("pixel_amp_array::fit"));
    }

    if(pixel_array<T>::verbose_level)
      cerr << "pixel_amp_array::fit - filling in arrays\n";
    int index;
    double amp, phase;
    for(int i=-this->axes[1]/2; i<this->axes[1]/2+x_extrapix; i++){
      for(int j=-this->axes[0]/2; j<this->axes[0]/2+y_extrapix; j++){
	index = (i+this->axes[1]/2)*this->axes[0]+j+this->axes[0]/2;

	phase = -xslope*(i+x_halfpix) - yslope*(j+y_halfpix);

	amp = this->pixeldata[index]/(double)nelem;
	this_arr[2*index] = amp*cos(phase);
	this_arr[2*index+1] = amp*sin(phase);

	amp = psf.pixeldata[index]/(double)nelem;
	psf_arr[2*index] = amp*cos(phase);
	psf_arr[2*index+1] = amp*sin(phase);
      }
    }

    if(pixel_array<T>::verbose_level)
      cerr << "pixel_amp_array::fit - performing transform\n";
    Arroyo::fft_manager<T> fft_mgr;  

    vector<long> flipped_axes(2, this->axes[0]);
    flipped_axes[0] = this->axes[1];
    fft_mgr.forward_fft(flipped_axes, false, true, this_arr);
    fft_mgr.forward_fft(flipped_axes, false, true, psf_arr);

    Arroyo::complex_cyclic_permutation(this->axes, 
				       this->axes[0]/2,
				       this->axes[1]/2, 
				       this_arr);

    Arroyo::complex_cyclic_permutation(psf.get_axes(), 
				       psf.get_axes()[0]/2,
				       psf.get_axes()[1]/2, 
				       psf_arr);

    if(pixel_array<T>::verbose_level)
      cerr << "pixel_amp_array::fit - allocating memory\n";
    T *real_cross_product, *imag_cross_product, *psf_power;
    try{
      real_cross_product = new T[nelem];
      imag_cross_product = new T[nelem];
      psf_power = new T[nelem];
    } catch(...) {
      cerr << "pixel_amp_array::fit error - error allocating memory\n";
      throw(string("pixel_amp_array::fit"));
    }

    if(pixel_array<T>::verbose_level)
      cerr << "pixel_amp_array::fit - computing intermediate arrays\n";
    double this_total_power=0, psf_total_power=0, cross_power=0;
    for(int i=-this->axes[1]/2; i<this->axes[1]/2+x_extrapix; i++){
      for(int j=-this->axes[0]/2; j<this->axes[0]/2+y_extrapix; j++){
	
	index = (i+this->axes[1]/2)*this->axes[0]+j+this->axes[0]/2;

	this_total_power += 
	  this_arr[2*index]*this_arr[2*index] + this_arr[2*index+1]*this_arr[2*index+1];

	psf_power[index] = 
	  (psf_arr[2*index]*psf_arr[2*index] + psf_arr[2*index+1]*psf_arr[2*index+1]);
	psf_total_power += psf_power[index];
	
	real_cross_product[index] = 
	  (this_arr[2*index]*psf_arr[2*index] + this_arr[2*index+1]*psf_arr[2*index+1]);
	imag_cross_product[index] = 
	  (this_arr[2*index]*psf_arr[2*index+1] - this_arr[2*index+1]*psf_arr[2*index]);
	
	cross_power += real_cross_product[index];
      }
    }

    if(pixel_array<T>::verbose_level)
      cout << "This total power "
	   << this_total_power 
	   << " PSF total power "
	   << psf_total_power
	   << " cross power "
	   << cross_power
	   << endl;

    // Perform a search for the best fitting fractional offset.
    double facx = 2*M_PI/(double)(this->axes[0]);
    double facy = 2*M_PI/(double)(this->axes[1]);
    double chisq, last_chisq=DBL_MAX;

    differential_amplitudes = 
      vector<double>(npsfs, fabs(cross_power/psf_total_power));
    relative_offsets.resize(npsfs);
    for(int i=0; i<npsfs; i++)
      relative_offsets[i].resize(2);

    vector<double> incr_differential_amplitudes(npsfs);
    vector<vector<double> > incr_relative_offsets(npsfs);
    for(int i=0; i<npsfs; i++)
      incr_relative_offsets[i].resize(2);



    // Setup stuff for the SVD
    if(pixel_array<T>::verbose_level)
      cout << "svd setup\n";

    char jobu[1];
    char jobvt[1];
    bool store_eigenmodes = true;
    if(store_eigenmodes)
      jobu[0] = jobvt[0] = 'A';
    else 
      jobu[0] = jobvt[0] = 'S';
    int lwork=-1, info;
    T optimal_workspace_size;
    T *singular_values, *g_gtranspose_eigenmode_data, *gtranspose_g_eigenmode_data, *work;


    T * matrix;
    int matrix_dimen = 3*npsfs;
    try{
      int nelem = matrix_dimen*matrix_dimen;
      matrix = new T[nelem];
      for(int i=0; i<nelem; i++)
	matrix[i] = 0;
    } catch(...){
      cerr << "pixel_amp_array::fit error - "
	   << "unable to allocate " 
	   << nelem
	   << " elements of memory\n";
      throw(string("pixel_amp_array::fit"));
    }

    int axes_0 = matrix_dimen;
    int axes_1 = matrix_dimen;
    
    singular_value_decomposition<T>(jobu, 
				    jobvt, 
				    axes_0,
				    axes_1, 
				    matrix,
				    axes_0, 
				    singular_values, 
				    g_gtranspose_eigenmode_data, 
				    axes_0, 
				    gtranspose_g_eigenmode_data, 
				    axes_1,
				    &optimal_workspace_size, 
				    lwork, 
				    info);
    
    lwork = (int)optimal_workspace_size;
    
    try{
      singular_values = new T[this->axes[1]];
      if(store_eigenmodes)
	g_gtranspose_eigenmode_data = new T[this->axes[0]*this->axes[0]];
      else
	g_gtranspose_eigenmode_data = new T[this->axes[0]*this->axes[1]];
      gtranspose_g_eigenmode_data = new T[this->axes[1]*this->axes[1]];
      work = new T[(int)(optimal_workspace_size)];
    } catch(...) {
      cerr << "pixel_amp_array::fit error - "
	   << "unable to allocate memory to perform the singular value decomposition\n";
      throw(string("pixel_amp_array::fit"));
    }


    double real, imag;
    double tmp, tmp2, tmp3;
    double angle, diff_angle, real_sum_exps, imag_sum_exps;
    double omega, xi;
    vector<double> dchisq_dx(npsfs), dchisq_dy(npsfs), dchisq_db(npsfs);

    while(1){

      // Compute matrix elements and first derivatives
      
      for(int r=0; r<npsfs; r++){
	for(int s=r; s<npsfs; s++){

	  dchisq_dx[r] = dchisq_dy[r] = dchisq_db[r] = 0;

	  for(int u=-this->axes[1]/2; u<this->axes[1]/2+x_extrapix; u++){
	    for(int v=-this->axes[0]/2; v<this->axes[0]/2+y_extrapix; v++){
	      index = (u+this->axes[1]/2)*this->axes[0]+v+this->axes[0]/2;

	      diff_angle = 
		-(relative_offsets[r][0]-relative_offsets[s][0])*facx*v + 
		(relative_offsets[r][1]-relative_offsets[s][1])*facy*u;

	      real_sum_exps = imag_sum_exps = 0;
	      for(int k=0; k<npsfs; k++){
		real_sum_exps += 
		  differential_amplitudes[k]*cos(relative_offsets[k][0]*facx*v+relative_offsets[k][1]*facy*u);
		imag_sum_exps += 
		  differential_amplitudes[k]*sin(relative_offsets[k][0]*facx*v+relative_offsets[k][1]*facy*u);
	      }


	      tmp = psf_power[index]*cos(diff_angle);
	      if(r==s){
		
		angle = -(relative_offsets[r][0]*facx*v + relative_offsets[r][1]*facy*u);

		tmp2 = 
		  real_cross_product[index]*sin(angle) +
		  imag_cross_product[index]*cos(angle) +
		  psf_power[index]*(cos(angle)*imag_sum_exps +
				    sin(angle)*real_sum_exps);

		tmp3 = 
		  real_cross_product[index]*cos(angle)-
		  imag_cross_product[index]*sin(angle)+
		  psf_power[index]*(cos(angle)*real_sum_exps -
				    sin(angle)*imag_sum_exps);
		  
		
		dchisq_dx[r] += differential_amplitudes[r]*v*differential_amplitudes[r]*tmp2;
		dchisq_dy[r] += differential_amplitudes[r]*u*differential_amplitudes[r]*tmp2;
		dchisq_db[r] += 
		  -(real_cross_product[index]*cos(angle)-
		  imag_cross_product[index]*sin(angle)) +
		  psf_power[index]*(cos(angle)*real_sum_exps -
				    sin(angle)*imag_sum_exps);
		
		omega = 
		  differential_amplitudes[r]*differential_amplitudes[r]*psf_power[index] -
		  differential_amplitudes[r]*tmp3;

		xi = tmp2;

	      } else {
		omega = differential_amplitudes[r]*differential_amplitudes[s]*tmp;
		xi = differential_amplitudes[s]*psf_power[index]*sin(diff_angle);
	      }

	      // d^{2}chi^{2}/dx_{r}dx_{s}
	      matrix[matrix_dimen*r+s] += v*v*omega;
	      if(r!=s)
		matrix[matrix_dimen*s+r] += v*v*omega;

	      // d^{2}chi^{2}/dx_{r}dy_{s}
	      matrix[matrix_dimen*r+s+npsfs] += u*v*omega;
	      matrix[matrix_dimen*(s+npsfs)+r] += u*v*omega;

	      // d^{2}chi^{2}/dx_{r}db_{s}
	      matrix[matrix_dimen*r+s+2*npsfs] += v*xi;
	      matrix[matrix_dimen*(s+2*npsfs)+r] += -v*xi;

	      // d^{2}chi^{2}/dy_{r}dy_{s}
	      matrix[matrix_dimen*(r+npsfs)+s+npsfs] += u*u*omega;
	      if(r!=s)
		matrix[matrix_dimen*(s+npsfs)+r+npsfs] += u*u*omega;

	      // d^{2}chi^{2}/dy_{r}db_{s}
	      matrix[matrix_dimen*(r+npsfs)+s+2*npsfs] += u*xi;
	      matrix[matrix_dimen*(s+2*npsfs)+r+npsfs] += -u*xi;

	      // d^{2}chi^{2}/db_{r}db_{s}
	      matrix[matrix_dimen*(r+2*npsfs)+s+2*npsfs] += tmp;
	      if(r!=s)
		matrix[matrix_dimen*(s+2*npsfs)+r+2*npsfs] += tmp;

	    }
	  }
	}
      }


      vector<vector<double> > tmp_matrix(matrix_dimen);
      for(int r=0; r<matrix_dimen; r++){
	tmp_matrix[r].resize(matrix_dimen);
	for(int s=0; s<matrix_dimen; s++){
	  tmp_matrix[r][s] = matrix[r*matrix_dimen+s];
	}
      }



      // Perform SVD    
      singular_value_decomposition<T>(jobu, 
				      jobvt, 
				      axes_0, 
				      axes_1, 
				      matrix,
				      axes_0, 
				      singular_values, 
				      g_gtranspose_eigenmode_data, 
				      axes_0, 
				      gtranspose_g_eigenmode_data, 
				      axes_1,
				      work, 
				      lwork, 
				      info);
    
      /*
      if(pixel_array<T>::verbose_level){
	cout << "GGtrans\n";
	for(int r=0; r<matrix_dimen; r++){
	  for(int s=0; s<matrix_dimen; s++){
	    cout << setw(15) << g_gtranspose_eigenmode_data[r*matrix_dimen+s]; 
	  }
	  cout << endl;
	}
      }

      if(pixel_array<T>::verbose_level){
	cout << "Singular vals\n";
	for(int r=0; r<matrix_dimen; r++){
	  cout << setw(15) << singular_values[r] << endl;
	}
      }

      if(pixel_array<T>::verbose_level){
	cout << "GtransG\n";
	for(int r=0; r<matrix_dimen; r++){
	  for(int s=0; s<matrix_dimen; s++){
	    cout << setw(15) << gtranspose_g_eigenmode_data[r*matrix_dimen+s]; 
	  }
	  cout << endl;
	}
      }

      if(pixel_array<T>::verbose_level){
	cout << "Reconstructed Matrix\n";
	for(int r=0; r<matrix_dimen; r++){
	  for(int s=0; s<matrix_dimen; s++){
	    tmp = 0;
	    for(int k=0; k<matrix_dimen; k++){
	      tmp += gtranspose_g_eigenmode_data[r*matrix_dimen+k]*
		g_gtranspose_eigenmode_data[s+k*matrix_dimen]*singular_values[k];
	      if(r==0 && s==0)
		cout << " r==s==0 " 
		     << k 
		     << "\t" << gtranspose_g_eigenmode_data[r*matrix_dimen+k]
		     << "\t" << singular_values[k] 
		     << "\t" << g_gtranspose_eigenmode_data[s+k*matrix_dimen] 
		     << "\t" << tmp << endl;
	    }
	    cout << setw(15) << tmp;
	  }
	  cout << endl;
	}
      }
      */

      // Form the inverse from V S^{+} U^{T}
      for(int r=0; r<matrix_dimen; r++){
	for(int s=0; s<matrix_dimen; s++){
	  tmp = 0;
	  for(int k=0; k<matrix_dimen; k++){
	    tmp += gtranspose_g_eigenmode_data[r*matrix_dimen+k]*
	      g_gtranspose_eigenmode_data[s+k*matrix_dimen]/singular_values[k];
	  }
	  matrix[s*matrix_dimen+r]=tmp;
	}
      }

      if(pixel_array<T>::verbose_level){
	cout << "Identity\n";
	for(int r=0; r<matrix_dimen; r++){
	  for(int s=0; s<matrix_dimen; s++){
	    double tmp = 0;
	    for(int t=0; t<matrix_dimen; t++)
	      tmp += tmp_matrix[r][t]*matrix[t*matrix_dimen+s];
	    cout << "\t" << tmp;
	  }
	  cout << endl;
	}
      }

      if(pixel_array<T>::verbose_level){
	cout << "Inverse\n";
	for(int r=0; r<matrix_dimen; r++){
	  for(int s=0; s<matrix_dimen; s++){
	    cout << setw(15) << matrix[r*matrix_dimen+s]; 
	  }
	  if(r<npsfs)
	    cout << "\t\t" << dchisq_dx[r] << endl;
	  else if(r>=npsfs && r<2*npsfs)
	    cout << "\t\t" << dchisq_dy[r-npsfs] << endl;
	  else
	    cout << "\t\t" << dchisq_db[r-2*npsfs] << endl;
	}
      }



      // Multiply initial estimate by inverse
      for(int r=0; r<npsfs; r++){
	incr_differential_amplitudes[r] = 0;
	incr_relative_offsets[r][0] = 0;
	incr_relative_offsets[r][1] = 0;
	for(int s=0; s<npsfs; s++){

	  incr_relative_offsets[r][0] += 
	    matrix[matrix_dimen*r+s]*dchisq_dx[s] + 
	    matrix[matrix_dimen*r+s+npsfs]*dchisq_dy[s] + 
	    matrix[matrix_dimen*r+s+2*npsfs]*dchisq_db[s];

	  incr_relative_offsets[r][1] += 
	    matrix[matrix_dimen*(r+npsfs)+s]*dchisq_dx[s] + 
	    matrix[matrix_dimen*(r+npsfs)+s+npsfs]*dchisq_dy[s] + 
	    matrix[matrix_dimen*(r+npsfs)+s+2*npsfs]*dchisq_db[s];

	  incr_differential_amplitudes[r] += 
	    matrix[matrix_dimen*(r+2*npsfs)+s]*dchisq_dx[s] + 
	    matrix[matrix_dimen*(r+2*npsfs)+s+npsfs]*dchisq_dy[s] + 
	    matrix[matrix_dimen*(r+2*npsfs)+s+2*npsfs]*dchisq_db[s];

	  /*
	  incr_relative_offsets[r][0] += 
	    matrix[matrix_dimen*s+r]*dchisq_dx[s] + 
	    matrix[matrix_dimen*s+r+npsfs]*dchisq_dy[s] + 
	    matrix[matrix_dimen*s+r+2*npsfs]*dchisq_db[s];

	  incr_relative_offsets[r][1] += 
	    matrix[matrix_dimen*(s+npsfs)+r]*dchisq_dx[s] + 
	    matrix[matrix_dimen*(s+npsfs)+r+npsfs]*dchisq_dy[s] + 
	    matrix[matrix_dimen*(s+npsfs)+r+2*npsfs]*dchisq_db[s];

	  incr_differential_amplitudes[r] += 
	    matrix[matrix_dimen*(s+2*npsfs)+r]*dchisq_dx[s] + 
	    matrix[matrix_dimen*(s+2*npsfs)+r+npsfs]*dchisq_dy[s] + 
	    matrix[matrix_dimen*(s+2*npsfs)+r+2*npsfs]*dchisq_db[s];
	  */

	}
      }

      // Compute chi squared
      double real_arg=0, imag_arg=0;
      chisq = 0;
      for(int u=-this->axes[1]/2; u<this->axes[1]/2+x_extrapix; u++){
	for(int v=-this->axes[0]/2; v<this->axes[0]/2+y_extrapix; v++){
	  index = (u+this->axes[1]/2)*this->axes[0]+v+this->axes[0]/2;
		
	  real = imag = 0;
	  for(int k=0; k<npsfs; k++){
	    angle = -(relative_offsets[k][0]*facx*v + relative_offsets[k][1]*facy*u);
	    real += differential_amplitudes[k]*cos(angle);
	    imag += differential_amplitudes[k]*sin(angle);
	  }
	  
	  real_arg = this_arr[2*index] - psf_arr[2*index]*real + psf_arr[2*index+1]*imag;
	  imag_arg = this_arr[2*index+1] - psf_arr[2*index+1]*real - psf_arr[2*index]*imag;

	  chisq += .5*(real_arg*real_arg + imag_arg*imag_arg);

	}
      }

      if(pixel_array<T>::verbose_level)
	cout << " chisq " << chisq << "\t";


      // Increment estimates
      for(int r=0; r<npsfs; r++){
	relative_offsets[r][0] += incr_relative_offsets[r][0];
	relative_offsets[r][1] += incr_relative_offsets[r][1];
	differential_amplitudes[r] += incr_differential_amplitudes[r];
      }

      // Compute chi squared
      chisq = 0;
      for(int u=-this->axes[1]/2; u<this->axes[1]/2+x_extrapix; u++){
	for(int v=-this->axes[0]/2; v<this->axes[0]/2+y_extrapix; v++){
	  index = (u+this->axes[1]/2)*this->axes[0]+v+this->axes[0]/2;
		
	  real = imag = 0;
	  for(int k=0; k<npsfs; k++){
	    angle = -(relative_offsets[k][0]*facx*v + relative_offsets[k][1]*facy*u);
	    real += differential_amplitudes[k]*cos(angle);
	    imag += differential_amplitudes[k]*sin(angle);
	  }
	  
	  real_arg = this_arr[2*index] - psf_arr[2*index]*real + psf_arr[2*index+1]*imag;
	  imag_arg = this_arr[2*index+1] - psf_arr[2*index+1]*real - psf_arr[2*index]*imag;

	  chisq += .5*(real_arg*real_arg + imag_arg*imag_arg);

	}
      }

      if(pixel_array<T>::verbose_level){
	cout << " chisq " << chisq;
	for(int k=0; k<npsfs; k++)
	  cout << " " << relative_offsets[k][0]
	       << " " << incr_relative_offsets[k][0]
	       << " " << relative_offsets[k][1]
	       << " " << incr_relative_offsets[k][1]
	       << " " << differential_amplitudes[k] 
	       << " " << incr_differential_amplitudes[k];
	cout << endl;
      }

      if(chisq > last_chisq)
	break;
    }

    // Transform back to get residual from fit
    if(pixel_array<T>::verbose_level)
      cerr << "pixel_amp_array::fit - computing residuals\n";

    for(int i=-this->axes[1]/2; i<this->axes[1]/2+x_extrapix; i++){
      for(int j=-this->axes[0]/2; j<this->axes[0]/2+y_extrapix; j++){
	index = (i+this->axes[1]/2)*this->axes[0]+j+this->axes[0]/2;
	real = imag = 0;

	for(int k=0; k<npsfs; k++){

	  amp = sqrt(psf_arr[2*index]*psf_arr[2*index] + 
		     psf_arr[2*index+1]*psf_arr[2*index+1]);
	  phase = atan2(psf_arr[2*index+1],psf_arr[2*index]);

	  real += amp*cos(phase);
	  imag += amp*sin(phase);
	  
	}
	
	fit_arr[2*index] = real;
	fit_arr[2*index+1] = imag;
	
	this_arr[2*index] -= real;
	this_arr[2*index+1] -= imag;
      }
    }


    Arroyo::complex_cyclic_permutation(this->get_axes(), 
				       -this->get_axes()[0]/2,
				       -this->get_axes()[1]/2, 
				       fit_arr);
    
    Arroyo::complex_cyclic_permutation(this->get_axes(), 
				       -this->get_axes()[0]/2,
				       -this->get_axes()[1]/2, 
				       this_arr);
    
    fft_mgr.backward_fft(flipped_axes, false, true, fit_arr);
    fft_mgr.backward_fft(flipped_axes, false, true, this_arr);

    double norm_fac = 1/(double)(this->total_space());
    for(int i=0; i<nelem; i++){
      fit_arr[i] = fit_arr[2*i]*norm_fac;
      this_arr[i] = this_arr[2*i]*norm_fac;
    }

    fitted_model = pixel_array<T>(this->axes, fit_arr);
    fit_residual = pixel_array<T>(this->axes, this_arr);

    // Clean up
    if(pixel_array<T>::verbose_level)
      cerr << "pixel_amp_array::fit - cleaning up\n";

    // SVD stuff
    delete [] work;
    delete [] singular_values;
    delete [] gtranspose_g_eigenmode_data;
    delete [] g_gtranspose_eigenmode_data;
    delete [] matrix;

    delete [] this_arr;
    delete [] psf_arr;
    delete [] fit_arr;
    delete [] real_cross_product;
    delete [] imag_cross_product;
    delete [] psf_power;

    return(chisq);
  }

  template<class precision>
    double pixel_amp_array<precision>::aperture_photometry(
				double x, double y, 
				double inner_radius, double outer_radius) const {

    // Ensure that specified region lies within the image
    if(x < 0 || x > this->axes[1] || y < 0 || y > this->axes[0]){
      cerr << "pixel_amp_array::aperture_photometry error - coordinates "
           << x << ", " << y
	   << " are out of range\n";
      this->print_axes(cerr);
      throw(std::string("pixel_amp_array::aperture_photometry"));
    }

    if(x-outer_radius < 0 || x+outer_radius > this->axes[1] ||
    			y-outer_radius < 0 || y + outer_radius > this->axes[0]){
      cerr << "pixel_amp_array::aperture_photometry error - coordinates "
           << x << ", " << y
	   << " place an aperture of radius " << outer_radius << " out of range\n";
      this->print_axes(cerr);
      throw(std::string("pixel_amp_array::aperture_photometry"));
    }

    // Do the calculation
    double inner_flux = 0, outer_flux = 0;
    int inner_pixel_count = 0, outer_pixel_count = 0;
    int nelems = this->total_space();
    double distance;
    int axes_zero = this->axes[0];
    for(int i=0; i<nelems; i++){
      distance = sqrt((i/axes_zero-x)*(i/axes_zero-x) + 
		      (i%axes_zero-y)*(i%axes_zero-y)); 
      if(distance < inner_radius && this->pixelwts[i]!=0){
	inner_flux += this->pixeldata[i];
	inner_pixel_count++;
      } else if ((distance > inner_radius)
      			&& (distance < outer_radius) && (this->pixelwts[i]!=0)) {
	outer_flux += this->pixeldata[i];
	outer_pixel_count++;
      }
    }
    
    outer_flux *= inner_pixel_count/(double)outer_pixel_count;
    inner_flux -= outer_flux;
    return(inner_flux);
  }

  template<class precision>
    template<class U>
    void pixel_amp_array<precision>::flag_outliers(
    			const pixel_amp_array<U> & pixamparr,
			const double & threshold){

    if(!pixamparr.weights_allocated()){
      cerr << "pixel_amp_array::flag_outliers error - "
           << "reference array has no weights\n";
      throw(std::string("pixel_amp_array::flag_outliers"));
    }

    if(!this->weights_allocated()) this->allocate_weights(1);

    std::vector<long> paxes = pixamparr.get_axes();
    if(paxes.size()!=this->axes.size()){
      cerr << "pixel_amp_array<precision>::flag_outliers error - "
           << "mismatched number of axes\n";
      throw(std::string("pixel_amp_array<precision>::flag_outliers"));
    }
    for(int i=0; i<this->axes.size(); i++){
      if(paxes[i]!=this->axes[i]){
	cerr << "pixel_amp_array<precision>::flag_outliers error - "
	     << "mismatched number of elements for axis "
	     << i << "\t" << paxes[i] << "\t" << this->axes[i] << endl;
	throw(std::string("pixel_amp_array<precision>::flag_outliers"));
      }
    }

    int nelems = this->total_space();
    double diff;
    for(int i=0; i<nelems; i++){
      if(this->pixelwts[i]!=0 && pixamparr.wt(i)!=0){
	diff = this->pixeldata[i] - pixamparr.data(i);
	// Here we have 3 tests
	// The first tests whether a pixel in this array is much
	// greater than the one in the arg array
	// The second tests whether the pixel in this array
	// is much greater than the threshold.
	// The third tests whether the pixel in the arg array
	// is less than the threshold.
	// The last two test attempt to ensure that we are not
	// flagging points on the star itself.
	if(fabs(diff) > threshold && 
	   this->pixeldata[i] > threshold && 
	   pixamparr.data(i) < threshold){
	  if(pixel_array<precision>::verbose_level)
	    cout << "pixel_amp_array::flag_outliers - flagging point " 
		 << i/this->axes[0] << "," << i%this->axes[0] << "\t" 
		 << this->pixeldata[i] << "\t" << pixamparr.data(i) << endl;
	  this->pixeldata[i] = 0;
	  this->pixelwts[i] = 0;
	}
      }
    }
  }

  template<class precision>
    std::vector<float> pixel_amp_array<precision>::azav(int xcenter, int ycenter){

    if(this->axes.size()!=2){
      cerr << "pixel_array::decimate - cannot decimate "
           << "pixel_array with number of dimensions " 
	   << this->axes.size() << endl;
      throw(std::string("pixel_array::decimate"));
    }
    if(xcenter<0 || xcenter>this->axes[0] || ycenter<0 || ycenter>this->axes[1]){
      cerr << "pixel_array::azav error - cannot average around point " 
	   << xcenter << ", " << ycenter << endl;
      throw(std::string("pixel_array::azav"));
    }

    int radius;
    if(xcenter<this->axes[1]-xcenter) radius = xcenter;
    else radius = this->axes[1]-xcenter;
    if(ycenter<radius) radius = ycenter;
    if(this->axes[0]-ycenter<radius) radius = this->axes[0] - ycenter;

    std::vector<float> azave(radius, 0), real(radius, 0), imag(radius,0);
    
    double rad, irad;
    for(int i=0; i<this->axes[1]; i++){
      for(int j=0; j<this->axes[0]; j++){
	rad = sqrt((i-xcenter)*(i-xcenter) + (j-ycenter)*(j-ycenter));
	irad = (int)(rad + .5);
	if(irad>=radius) continue;
	if(this->weights_allocated()){
	  if(this->pixelwts[i*(this->axes[0])+j] != 0){
	    real[irad] += cos((double)(this->pixeldata[i*(this->axes[0])+j]));
	    imag[irad] += sin((double)(this->pixeldata[i*(this->axes[0])+j]));
	  }
	} else {
	  real[irad] += cos((double)(this->pixeldata[i*(this->axes[0])+j]));
	  imag[irad] += sin((double)(this->pixeldata[i*(this->axes[0])+j]));
	}	  
      }
    }
    
    for(int i=0; i<azave.size(); i++)
      azave[i] = atan2(imag[i],real[i]);

    return(azave);
  }

  template <class precision>
    pixel_amp_array<precision> operator + (const pixel_amp_array<precision> &p1,
    					const pixel_amp_array<precision> &p2){
    if(p1.get_axes()!=p2.get_axes()){
      cerr << "pixel_amp_array<precision>::operator+ error - " 
	   << "mismatched array sizes:\n";
      throw(std::string("pixel_amp_array<precision>::operator+"));
    }
    pixel_amp_array<precision> pixamparr(p1);
    pixamparr += p2;
    return(pixamparr);
  }

  template <class precision>
    pixel_amp_array<precision> operator - (const pixel_amp_array<precision> &p1,
    					const pixel_amp_array<precision> &p2){
    if(p1.get_axes()!=p2.get_axes()){
      cerr << "pixel_amp_array<precision>::operator- error - " 
	   << "mismatched array sizes:\n";
      throw(std::string("pixel_amp_array<precision>::operator-"));
    }
    pixel_amp_array<precision> pixamparr(p1);
    pixamparr -= p2;
    return(pixamparr);
  }

  template <class precision>
    pixel_amp_array<precision> operator * (const pixel_amp_array<precision> &p1,
    					const pixel_amp_array<precision> &p2){
    if(p1.get_axes()!=p2.get_axes()){
      cerr << "pixel_amp_array<precision>::operator* error - " 
	   << "mismatched array sizes:\n";
      throw(std::string("pixel_amp_array<precision>::operator/"));
    }
    pixel_amp_array<precision> pixamparr(p1);
    pixamparr *= p2;
    return(pixamparr);
  }

  template <class precision>
    pixel_amp_array<precision> operator / (const pixel_amp_array<precision> &p1,
    					const pixel_amp_array<precision> &p2){
    if(p1.get_axes()!=p2.get_axes()){
      cerr << "pixel_amp_array<precision>::operator/ error - " 
	   << "mismatched array sizes:\n";
      throw(std::string("pixel_amp_array<precision>::operator/"));
    }
    pixel_amp_array<precision> pixamparr(p1);
    pixamparr /= p2;
    return(pixamparr);
  }

  template <class precision>
    pixel_amp_array<precision> operator + (const pixel_amp_array<precision> &p1, double & fac){
    pixel_amp_array<precision> pixamparr(p1);
    pixamparr += fac;
    return(pixamparr);
  }

  template <class precision>
    pixel_amp_array<precision> operator - (const pixel_amp_array<precision> &p1, double & fac){
    pixel_amp_array<precision> pixamparr(p1);
    pixamparr -= fac;
    return(pixamparr);
  }

  template <class precision>
    pixel_amp_array<precision> operator * (const pixel_amp_array<precision> &p1, double & fac){
    pixel_amp_array<precision> pixamparr(p1);
    pixamparr *= fac;
    return(pixamparr);
  }

  template <class precision>
    pixel_amp_array<precision> operator / (const pixel_amp_array<precision> &p1, double & fac){
    pixel_amp_array<precision> pixamparr(p1);
    pixamparr /= fac;
    return(pixamparr);
  }

  ///////////////////////////////////////////
  ///  Operator ==  for pixel_amp_array
  template<class precision>
    bool operator != (const pixel_amp_array<precision> &p1, const pixel_amp_array<precision> &p2){
    return(!(p1==p2));
  }

}

#endif
