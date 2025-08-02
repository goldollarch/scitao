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
#include <algorithm>
#include "special_functions.h"
#include "iofits.h"
#include "fits_factory.h"
#include "fits_header_data.h"
#include "pixel_array.h"
#include "power_spectrum.h"
#include "structure_function.h"
#include "refractive_atmospheric_layer.h"

using namespace std;

namespace Arroyo {

  namespace factory_register {
    const fits_keyval_set & get_structure_function_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "structure function"));
      return *fkvs;
    }
    
    AO_sim_base * create_structure_function(const iofits & iof) {
      return new structure_function(iof);
    }
  } 

  const bool structure_function::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_structure_function_keyval_set(), 
				      factory_register::create_structure_function);

  structure_function::structure_function() : pixel_array<double>(vector<long>(1,0)) {
    pixel_scale = 0;
  }

  structure_function::
    structure_function(const structure_function & struct_fn){
    this->operator=(struct_fn);
  }

  structure_function::structure_function(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "structure_function::read - "
	   << "error opening file " << filename << endl;
      throw(string("structure_function::read"));
    }
    this->read(iof);
  }

  structure_function::structure_function(const iofits & iof){
    this->read(iof);
  }

  structure_function::
    structure_function(const pixel_array<double> & pixarr, double pixscale) :
      pixel_array<double>(pixarr) {
    if(pixarr.get_axes().size()!=1 && pixarr.get_axes().size()!=2){
      cerr << "structure_function::structure_function error - "
	   << "cannot allocate structure function with "
	   << pixarr.get_axes().size() << " axes\n";
      throw(string("structure_function::structure_function"));
    }
    pixel_scale = pixscale;
  }
  
  structure_function::structure_function(const power_spectrum * pspectrum,
					 vector<long> axes, 
					 double pixscale) {

    // Check arguments

    if(axes.size()!=1 && axes.size()!=2){
      cerr << "structure_function::structure_function error - "
	   << "cannot get a structure function with " << axes.size() << " axes\n";
      throw(string("structure_function::structure_function"));
    }

    for(int i=0; i<axes.size(); i++){
      if(axes[i]<=0){
	cerr << "structure_function::structure_function error - "
	     << "axis " << i+1 << " has value " << axes[i] << endl;
	throw(string("structure_function::structure_function"));
      }    
    }

    if(pixscale<=0){
      cerr << "structure_function::structure_function error - "
	   << "cannot construct structure function with pixel scale " << pixscale << endl; 
      throw(string("structure_function::structure_function"));
    }

    
    // Check that the power spectrum is an isotropic power law with no
    // outer scale or a von Karman outer scale, and with exponent -11/3

    const isotropic_power_law_spectrum<power_law, null_inner_scale> * no_outer_isops =
      dynamic_cast<const isotropic_power_law_spectrum<power_law, null_inner_scale> *>(pspectrum);

    const isotropic_power_law_spectrum<von_karman_power_law, null_inner_scale> * vk_outer_isops =
      dynamic_cast<const isotropic_power_law_spectrum<von_karman_power_law, null_inner_scale> *>(pspectrum);

    double exponent;

    if(no_outer_isops==NULL && vk_outer_isops==NULL){
	cerr << "structure_function::structure_function error - "
	     << "do not yet know how to construct a structure function "
	     << "for power laws other than those with no outer scale or a von Karman outer scale\n";
	throw(string("structure_function::structure_function"));
    } else if(no_outer_isops!=NULL) exponent = no_outer_isops->get_power_law().get_exponent();
    else if(vk_outer_isops!=NULL) exponent = vk_outer_isops->get_power_law().get_exponent();
      
    if(fabs(exponent+11/3.0)>1e-14*11/3.0){
      cerr << "structure_function::structure_function error - "
	   << "do not yet know how to construct structure function "
	   << " with exponent " << exponent << endl;
      throw(string("structure_function::structure_function"));
    }
      
    // initialize the data members and fill in the pixel data
    
    this->pixel_scale = pixscale;
    this->pixel_array<double>::set_axes(axes);

    if(vk_outer_isops!=NULL){
      // see Johannson & Gavel
      double tmp = 2.61*vk_outer_isops->get_power_law().get_coefficient()/2.0/.033/M_PI;
      double constant_term = tmp*.6*pow((vk_outer_isops->get_power_law().get_outer_scale()/2/M_PI), 5/3.0);
      double bessel_prefac = tmp*pow((vk_outer_isops->get_power_law().get_outer_scale()/4/M_PI), 5/6.0)/gamma_function(11/6.0);

      if(axes.size()==1){
	this->pixeldata[0] = 0;
	for(int i=1; i<axes[0]; i++)
	  this->pixeldata[i] = constant_term - 
	    bessel_prefac*pow(pixscale*i, 5/6.0)*
	    bessel_Knu(5/6.0, 2*M_PI*pixscale*i/vk_outer_isops->get_power_law().get_outer_scale());
      } else {
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
	double pixel_location;
	for(int i=-axes[1]/2; i<axes[1]/2+x_extrapix; i++){
	  for(int j=-axes[0]/2; j<axes[0]/2+y_extrapix; j++){
	    pixel_location = pixscale*sqrt((i+x_halfpix)*(i+x_halfpix)+
					   (j+y_halfpix)*(j+y_halfpix));
	    
	    
	    if(pixel_location==0) this->pixeldata[(i+axes[1]/2)*axes[0]+j+axes[0]/2] = 0;
	    else 
	      this->pixeldata[(i+axes[1]/2)*axes[0]+j+axes[0]/2] = constant_term - 
		bessel_prefac*pow(pixel_location, 5/6.0)*
		bessel_Knu(5/6.0, 2*M_PI*pixel_location/vk_outer_isops->get_power_law().get_outer_scale());
	  }
	}
      }
    } else {
      double coefficient = 2.91*no_outer_isops->get_power_law().get_coefficient()/2.0/.033/M_PI;      

      if(axes.size()==1){
	for(int i=0; i<axes[0]; i++)
	  this->pixeldata[i] = coefficient*pow(pixscale*i,5/3.0);
      } else {
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
	double pixel_location;
	for(int i=-axes[1]/2; i<axes[1]/2+x_extrapix; i++){
	  for(int j=-axes[0]/2; j<axes[0]/2+y_extrapix; j++){
	    pixel_location = pixscale*sqrt((i+x_halfpix)*(i+x_halfpix)+
					   (j+y_halfpix)*(j+y_halfpix));	    
	    this->pixeldata[(i+axes[1]/2)*axes[0]+j+axes[0]/2] = coefficient*pow(pixel_location, 5/3.0);
	  }
	}
      }
    }
  }


  structure_function::structure_function(const power_spectrum * pspectrum,
					 const subharmonic_method & subm,
					 vector<long> axes, 
					 double pixscale){

    // Check arguments

    if(axes.size()!=1 && axes.size()!=2){
      cerr << "isotropic_power_law_spectrum::get_expected_structure_function error - "
	   << "cannot get a structure function with " << axes.size() << " axes\n";
      throw(string("isotropic_power_law_spectrum::get_expected_structure_function"));
    }

    for(int i=0; i<axes.size(); i++){
      if(axes[i]<=0){
	cerr << "isotropic_power_law_spectrum::get_expected_structure_function error - "
	     << "axis " << i+1 << " has value " << axes[i] << endl;
	throw(string("isotropic_power_law_spectrum::get_expected_structure_function"));
      }    
    }

    if(pixscale<=0){
      cerr << "isotropic_power_law_spectrum::get_expected_structure_function error - "
	   << "cannot construct structure function with pixel scale " << pixscale << endl; 
      throw(string("isotropic_power_law_spectrum::get_expected_structure_function"));
    }

    vector<long> two_d_axes;

    // Here we fix up the axes so that they correspond to the
    // intrinsic dimensionality of the subharmonic method.
    if(axes.size()==1){
      long dimen = 2*axes[0];
      if(subm.intrinsic_dimensionality()!=-1){
	if(dimen%2!=subm.intrinsic_dimensionality()) 
	  dimen++;
      } else 
	if(axes[0]%2 != dimen%2) dimen++;
      two_d_axes = vector<long>(2,dimen);
    } else {
      two_d_axes = axes;
      if(subm.intrinsic_dimensionality()!=-1){
	for(int i=0; i<two_d_axes.size(); i++){
	  if(two_d_axes[i]%2 != subm.intrinsic_dimensionality())
	    two_d_axes[i]++;
	}
      }
    }

    // initialize data members    
    this->pixel_scale = pixscale;
    this->pixel_array<double>::set_axes(two_d_axes);



    // FIX THIS UP LATER

    double * data;
    //try{data = new double[2*two_d_axes[0]*two_d_axes[1]];}
    alloc_size sz(two_d_axes[0], two_d_axes[1], 2);
    try{
      data = new double[sz];
    }
    catch(...){
      cerr << "isotropic_power_law_spectrum::get_expected_structure_function error - "
	   << "could not allocate memory for structure function\n";
      throw;
    }    

    // Initialize a frequency space array.
    // This array has real values equal to the sqrt of the power
    // spectral density, and imaginary values equal to zero.
    initialize_frequency_array<double>(data, pspectrum, two_d_axes, pixscale, false, subm);

    // Construct the autocorrelation function by fourier transforming
    // the power spectral density.  

    fft_manager<double> fft_mgr;
    halfpixel_fft<double>(two_d_axes, data, fft_mgr);

    for(int i=0; i<two_d_axes[1]; i++)
      for(int j=0; j<two_d_axes[0]; j++)
	data[i*two_d_axes[0]+j] = data[2*(i*two_d_axes[0]+j)];

    pixel_array<double> pixarr;
    if(two_d_axes.size()==1){
      int dimen = two_d_axes[0];
      // At this point the real correlation function is stored in the
      // first axes[0]*axes[0] elements of the data array.  Now we 
      // convert to a 1d strucure function with axes[0]/2+1 elements
      int central_index;
      if(dimen%2) central_index = (dimen*dimen)/2;
      else central_index = (dimen+1)*(dimen/2)-1;
      double central_element = data[central_index];
      for(int i=0; i<axes[0]; i++)
	data[i] = 2*(central_element - data[central_index+i]);
      
      pixarr = pixel_array<double>(vector<long>(1, axes[0]), data);
    } else {
      double central_element = data[two_d_axes[1]/2*two_d_axes[0]+two_d_axes[0]/2];
      if(two_d_axes[0]%2) central_element = data[two_d_axes[1]/2*two_d_axes[0]+two_d_axes[0]/2 + 1];

      for(int i=0; i<two_d_axes[1]; i++)
	for(int j=0; j<two_d_axes[0]; j++)
	  data[i*two_d_axes[0]+j] = 2*(central_element - data[i*two_d_axes[0]+j]);
      this->pixel_array<double>::operator=(pixel_array<double>(two_d_axes, data));
    }
  }


  structure_function & structure_function::
    operator=(const structure_function & struct_fn){
    if(this==&struct_fn)
      return(*this);
    this->pixel_array<double>::operator=(struct_fn);
    pixel_scale = struct_fn.pixel_scale;
    return(*this);
  }

  void structure_function::set_axes(const vector<long> & in_axes){
    if(in_axes.size()!=1 && in_axes.size()!=2){
      cerr << "structure_function::set_axes error - "
	   << "cannot set axes to dimension "
	   << in_axes.size() << " - must be 1 or 2 dimensional\n";
      throw(string("structure_function::set_axes"));
    }
    pixel_array<double>::set_axes(in_axes);
  }

  void structure_function::read(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "structure_function::read - "
	   << "error opening file " << filename << endl;
      throw(string("structure_function::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "structure_function::read - "
	   << "error reading "
	   << this->unique_name() << " from file " 
	   << filename << endl;
      throw(string("structure_function::read"));
    }
  }

  void structure_function::read(const iofits & iof) {
    if(!iof.key_exists("TYPE")){
      cerr << "structure_function::read error - "
	   << "unrecognized type of file\n";
      throw(string("structure_function::read"));
    }
    string type, comment;
    iof.read_key("TYPE", type, comment);
    if(type!=this->unique_name()){
      cerr << "structure_function::read error - file of type " 
	   << type << " rather than of type " 
	   << this->unique_name() << endl;
      throw(string("structure_function::read"));
    }
    iof.read_key("PIXSCALE", pixel_scale, comment);
    this->pixel_array<double>::read(iof);
  
    if(iof.get_hdu_num()<iof.get_num_hdus())
      iof.movrel_hdu(1);

  }

  void structure_function::write(const char * filename) const {
    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "structure_function::write - "
	   << "error opening file " << filename << endl;
      throw(string("structure_function::write"));
    }
    try{this->write(iof);}
    catch(...){
      cerr << "structure_function::write - "
	   << "error writing "
	   << this->unique_name() << " to file "
	   << filename << endl;
      throw(string("structure_function::write"));
    }
  }

  void structure_function::write(iofits & iof) const {

    fits_header_data<double> fhd(this->axes);
    fhd.write(iof);

    string type = this->unique_name();
    string comment = "object type";
    iof.write_key("TYPE", type, comment);
    iof.write_key("PIXSCALE", pixel_scale, "pixel scale (meters)");

    // This is a bit of an inefficient hack to get
    // access to protected pixel_array::normalize_by_weights,
    // while at the same time avoiding the const declaration 
    // in the write member function.
    const structure_function tmpstrfn(*this);
    const_cast<structure_function *>(this)->pixel_array<double>::normalize_by_wts();
    this->pixel_array<double>::write(iof);
    const_cast<structure_function *>(this)->operator=(tmpstrfn);
  }

  void structure_function::
  print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << this->unique_name()
       << "/" << "object type" << endl;
    fits_header_data<double> fhd(this->axes);
    fhd.print(os, prefix);
    os << prefix << "PIXSCALE   = " << setw(vlspc) << pixel_scale
       << "/" << "pixel scale (meters)" << endl;
  }

#if 0
  template<class T>
  void structure_function::
  add_statistics(const refractive_atmospheric_layer<T> & ref_atm_layer){

    vector<long> ref_layer_axes = 
      ref_atm_layer.get_axes();    
    double ref_layer_pixel_scale = ref_atm_layer.get_pixel_scale();

    if(this->get_axes().size()==1){

      // If this is the first time, and weights haven't been allocated,
      // we reinitialize the pixel array and set the weights
      // Also, if the weights have been allocated but the ref_atm_layer
      // is of a different size, we reinitialize as above
      vector<long> sorted_axes = ref_layer_axes;
      sort(sorted_axes.begin(), sorted_axes.end());

      // if the long axis is even, subtract one element.
      long truncated_axis = sorted_axes[0];
      if(sorted_axes[0]%2==0) truncated_axis--;
      long half_truncated_axis = truncated_axis/2;

      if(!this->weights_allocated() || 
	 this->get_axes()[0]!=half_truncated_axis+1 ||
	 ref_layer_pixel_scale!=this->pixel_scale){
	this->set_axes(vector<long>(1,half_truncated_axis+1));
	this->allocate_weights(0);
	this->pixel_scale = ref_layer_pixel_scale;
      }

      double tmp, centerval;
      if(ref_layer_axes[0]>=ref_layer_axes[1]){
	for(int i=0; i<ref_layer_axes[1]; i++){
	  centerval = ref_atm_layer.data(i*ref_layer_axes[0]+half_truncated_axis);
	  for(int j=0; j<truncated_axis; j++){
	    tmp = ref_atm_layer.data(i*ref_layer_axes[0] + j) - centerval;
	    this->pixeldata[abs(j-half_truncated_axis)] += tmp*tmp;
	    this->pixelwts[abs(j-half_truncated_axis)]++;
	  }
	}
      } else {
	for(int j=0; j<ref_layer_axes[0]; j++){
	  centerval = ref_atm_layer.data((half_truncated_axis)*ref_layer_axes[0]+j);
	  for(int i=0; i<truncated_axis; i++){	
	    tmp = ref_atm_layer.data(i*ref_layer_axes[0] + j) - centerval;
	    this->pixeldata[abs(i-half_truncated_axis)] += tmp*tmp;
	    this->pixelwts[abs(i-half_truncated_axis)]++;
	  }
	}
      }

      /*
	double tmp;
	for(int i=0; i<ref_layer_axes[1] && i<sorted_axes[0]; i++){

	for(int j=i; j<ref_layer_axes[0]; j++){
	tmp = 
	ref_atm_layer.pixel_array<double>::data(i*ref_layer_axes[0] + j) -
	ref_atm_layer.pixel_array<double>::data(i*ref_layer_axes[0] + i);
	this->pixel_array<double>::pixeldata[j-i] += tmp*tmp;
	this->pixel_array<double>::pixelwts[j-i]++;
	}

	for(int j=i; j<ref_layer_axes[1]; j++){
	tmp = 
	ref_atm_layer.pixel_array<double>::data(i*ref_layer_axes[0] + j) -
	ref_atm_layer.pixel_array<double>::data(i*ref_layer_axes[0] + i);
	this->pixel_array<double>::pixeldata[j-i] += tmp*tmp;
	this->pixel_array<double>::pixelwts[j-i]++;
	}
	}
      */

    } else {

      // As above, check for dimensional consistency with the ref_atm_layer
      vector<long> midpoint(2);
      if(!this->weights_allocated() || 
	 this->get_axes()!=ref_layer_axes){
	for(int i=0; i<ref_layer_axes.size(); i++){
	  if(ref_layer_axes[i]%2==0) midpoint[i] = ref_layer_axes[i]/2-1;
	  else midpoint[i] = ref_layer_axes[i]/2;
	}
	this->set_axes(ref_layer_axes);
	this->pixel_array<double>::allocate_weights(0);
      }

      double tmp;
      double midpoint_val = 
	ref_atm_layer.pixel_array<double>::data(midpoint[1]*ref_layer_axes[0] + 
						midpoint[0]);
      for(int i=0; i<ref_layer_axes[1]; i++){
	for(int j=0; j<ref_layer_axes[0]; j++){
	  tmp = 
	    ref_atm_layer.pixel_array<double>::data(i*ref_layer_axes[0] + j) -
	    midpoint_val;
	  this->pixel_array<double>::pixeldata[i*ref_layer_axes[0] + j] +=
	    tmp*tmp;
	  this->pixel_array<double>::pixelwts[i*ref_layer_axes[0] + j]++;
	}
      }
    }
  }
#endif

}
