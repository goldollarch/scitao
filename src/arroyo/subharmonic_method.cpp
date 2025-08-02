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

#include <typeinfo>
#include <cmath>
#include <iostream>
#include <iomanip>
#include "AO_algo.h"
#include "three_frame.h"
#include "fits_factory.h"
#include "power_spectrum.h"
#include "subharmonic_method.h"

using namespace std;

namespace Arroyo {

  bool operator==(const subharmonic_method & subm1, const subharmonic_method & subm2) {
    return typeid(subm1)==typeid(subm2) && subm1.equal(subm2);
  }

  bool operator!=(const subharmonic_method & subm1, const subharmonic_method & subm2) {
    return !operator==(subm1, subm2);
  }

  int subharmonic_method::verbose_level = 0;

  subharmonic_method * subharmonic_method::subharmonic_method_factory(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...) {
      cerr << "subharmonic_method::subharmonic_method_factory - "
	   << "error opening file " << filename << endl;
      throw(string("subharmonic_method::subharmonic_method_factory"));
    }
    return(subharmonic_method::subharmonic_method_factory(iof));
  }

  subharmonic_method * subharmonic_method::subharmonic_method_factory(const iofits & iof){
    AO_sim_base * aosb = AO_sim_base::AO_sim_base_factory(iof);
    subharmonic_method * subm = dynamic_cast<subharmonic_method *>(aosb);
    if(subm==NULL)
      throw(string("subharmonic_method::subharmonic_method_factory"));
     return(subm);    
  }

  void subharmonic_method::write(iofits & iof) const {
    fits_header_data<char> tmphdr;
    tmphdr.write(iof);
  }

  namespace factory_register {
    const fits_keyval_set & get_null_subharmonic_method_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "null subharmonic method"));
      return *fkvs;
    }
    
    AO_sim_base * create_null_subharmonic_method(const iofits & iof) {
      return new null_subharmonic_method(iof);
    }
  }

  const bool null_subharmonic_method::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_null_subharmonic_method_keyval_set(), 
				      factory_register::create_null_subharmonic_method);

  void null_subharmonic_method::read(const iofits & iof){
    if(!iof.key_exists("TYPE")){
      cerr << "null_subharmonic_method::read error - "
	   << "unrecognized type of file\n";
      throw(string("null_subharmonic_method::read"));
    }
    // leave iof pointing to the next header, if it exists
    if(iof.get_hdu_num()!=iof.get_num_hdus())
      iof.movrel_hdu(1);
    return;
  }

  void null_subharmonic_method::write(const char * filename) const {
    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "null_subharmonic_method::write - "
	   << "error opening file " << filename << endl;
      throw(string("null_subharmonic_method::write"));
    }
    
    try{this->write(iof);}
    catch(...){
      cerr << "null_subharmonic_method::write - "
	   << "error writing "
	   << this->unique_name() << " to file "
	   << filename << endl;
      throw(string("null_subharmonic_method::write"));
    }
  }
  
  void null_subharmonic_method::write(iofits & iof) const {
    this->subharmonic_method::write(iof);
    string type = "null subharmonic method";
    string comment = "object type";
    iof.write_key("TYPE", type, comment);
  }

  void null_subharmonic_method::print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << "Null subharmonic method"
       << "/" << "object type" << endl;
  }

  namespace factory_register {
    const fits_keyval_set & get_generalized_subharmonic_method_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "generalized subharmonic method"));
      return *fkvs;
    }
    
    AO_sim_base * create_generalized_subharmonic_method(const iofits & iof) {
      return new generalized_subharmonic_method(iof);
    }
  }

  const bool generalized_subharmonic_method::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_generalized_subharmonic_method_keyval_set(), 
				      factory_register::create_generalized_subharmonic_method);

  generalized_subharmonic_method::generalized_subharmonic_method(long sbdpth, 
								 long nsubpix_per_lvl,
								 long nsubpix_per_pix){
    if(sbdpth<0){
      cerr << "generalized_subharmonic_method::generalized_subharmonic_method error - "
	   << "can't construct with subharmonic_level " << sbdpth 
	   << " less than zero\n";
      throw(string("generalized_subharmonic_method::generalized_subharmonic_method"));
    }
    
    if(sbdpth==0){
      this->operator=(generalized_subharmonic_method());
      return;
    }

    if(nsubpix_per_pix<=1 || nsubpix_per_pix%2==0){
      cerr << "generalized_subharmonic_method::generalized_subharmonic_method error - "
	   << "this function requires nsubpix_per_pix odd and >= 3.\n"
	   << "Instead, " << nsubpix_per_pix << " subpixels per pixel supplied to this function\n";
      throw(string("generalized_subharmonic_method::generalized_subharmonic_method"));
    }
    
    if(nsubpix_per_lvl<nsubpix_per_pix || nsubpix_per_lvl>=3*nsubpix_per_pix || nsubpix_per_lvl%2==0){
      cerr << "generalized_subharmonic_method::generalized_subharmonic_method error -\n"
	   << "this function requires an odd number of subpixels per level greater than "
	   << "the number of subpixels per pixel and less than three times the number of subpixels per pixel\n" 
	   << "Instead, " << nsubpix_per_lvl << " subpixels per level supplied to this function\n";
      throw(string("generalized_subharmonic_method::generalized_subharmonic_method"));
    }
    
    subharmonic_depth = sbdpth;
    nsubpixels_per_pixel = nsubpix_per_pix;
    nsubpixels_per_level = nsubpix_per_lvl;
  }

  void generalized_subharmonic_method::read(const char * filename) {
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "generalized_subharmonic_method::read - "
	   << "error opening file " << filename << endl;
      throw(string("generalized_subharmonic_method::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "generalized_subharmonic_method::read - "
	   << "error reading "
	   << this->unique_name() << " from file " 
	   << filename << endl;
      throw(string("generalized_subharmonic_method::read"));
    }
  }

  void generalized_subharmonic_method::read(const iofits & iof) {
    if(!iof.key_exists("TYPE")){
      cerr << "generalized_subharmonic_method::read error - "
	   << "unrecognized type of file\n";
      throw(string("generalized_subharmonic_method::read"));
    }
    string type, comment;
    iof.read_key("TYPE", type, comment);
    if(type!="generalized subharmonic method"){
      cerr << "generalized_subharmonic_method::read error - file of type " 
	   << type << " rather than of type circular aperture\n";
      throw(string("generalized_subharmonic_method::read"));
    }
    iof.read_key("SBHMDPTH", subharmonic_depth, comment);
    iof.read_key("NPXPERPX", nsubpixels_per_pixel, comment);
    iof.read_key("NPXPERLV", nsubpixels_per_level, comment);
    // leave iof pointing to the next header, if it exists
    if(iof.get_hdu_num()!=iof.get_num_hdus())
      iof.movrel_hdu(1);

  }

  void generalized_subharmonic_method::write(const char * filename) const {
    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "generalized_subharmonic_method::write - "
	   << "error opening file " << filename << endl;
      throw(string("generalized_subharmonic_method::write"));
    }

    try{this->write(iof);}
    catch(...){
      cerr << "generalized_subharmonic_method::write - "
	   << "error writing "
	   << this->unique_name() << " to file "
	   << filename << endl;
      throw(string("generalized_subharmonic_method::write"));
    }
  }

  void generalized_subharmonic_method::write(iofits & iof) const {
    this->subharmonic_method::write(iof);
    string type = "generalized subharmonic method";
    string comment = "object type";
    iof.write_key("TYPE", type, comment);
    comment = "subharmonic depth";
    iof.write_key("SBHMDPTH", subharmonic_depth, comment);
    comment = "number of subpixels per pixel";
    iof.write_key("NPXPERPX", nsubpixels_per_pixel, comment);
    comment = "number of subpixels per level";
    iof.write_key("NPXPERLV", nsubpixels_per_level, comment);
  }

  void generalized_subharmonic_method::print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << "generalized subharmonic method"
       << "/" << "object type" << endl;
    os << prefix << "SBHMDPTH   = " << setw(vlspc) << subharmonic_depth
       << "/" << "subharmonic depth" << endl;
    os << prefix << "NPXPERPX   = " << setw(vlspc) << nsubpixels_per_pixel
       << "/" << "number of subpixels per pixel" << endl;
    os << prefix << "NPXPERLV   = " << setw(vlspc) << nsubpixels_per_level
       << "/" << "number of subpixels per level" << endl;
  };
  
  template<class T>
  void generalized_subharmonic_method::apply_subharmonic_correction_impl(
  				const power_spectrum & pspec, 
				const vector<long> & axes, 
				double pixscale, 
				bool random,
				T * data) const {

    if(this->subharmonic_depth == 0) return;
    if(axes.size()!=2 || axes[0]<=0 || axes[1]<=0){
      cerr << "generalized_subharmonic_method::apply_subharmonic_correction error - "
	   << "axes supplied to this function are incomprehensible\n";
      cerr << "size " << axes.size() << " first " << axes[0] << " second " << axes[1] << endl;
      throw(string("generalized_subharmonic_method::apply_subharmonic_correction"));
    }

    if(axes[0]%2==0 || axes[1]%2==0){
      cerr << "generalized_subharmonic_method::apply_subharmonic_correction error - "
	   << "cannot perform correction on array with even dimensional axes: "
	   << axes[0] << " x " << axes[1] << endl;
      throw(string("generalized_subharmonic_method::apply_subharmonic_correction"));
    }
      
    if(pixscale <= 0){
      cerr << "generalized_subharmonic_method::apply_subharmonic_correction error - "
	   << "nonpositive pixel scale " << pixscale << " supplied to this function\n";
      throw(string("generalized_subharmonic_method::apply_subharmonic_correction"));
    }

    if(axes[0]<nsubpixels_per_level*nsubpixels_per_pixel ||
       axes[1]<nsubpixels_per_level*nsubpixels_per_pixel){
      cerr << "generalized_subharmonic_method::apply_subharmonic_correction error - "
	   << "axes " << axes[0] << "x" << axes[1] << " too small to be corrected using "
	   << nsubpixels_per_level << " subpixels per level, with " << nsubpixels_per_pixel
	   << " subpixels per pixel\n";
      throw(string("generalized_subharmonic_method::apply_subharmonic_correction"));
    }

    // The halfpixel information - this is always
    // the same because both the axes and the number
    // of subpixels per layer are always odd in this 
    // function.
    double x_halfpix=0, y_halfpix=0;
    int x_extrapix=1, y_extrapix=1;

    // Set up the weights for the subharmonic levels.  These account
    // for potential overlap from subsequent subharmonic levels, as in
    // the Johansson and Gavel method. 
    vector<vector<double> > subharmonic_pixel_weights(nsubpixels_per_level, 
						      vector<double>(nsubpixels_per_level,1));
    double subharmonic_layer_boundary_in_pixels = nsubpixels_per_level/(double)nsubpixels_per_pixel/2.0;
    double tmp, weight;
    for(int l=-nsubpixels_per_level/2; l<nsubpixels_per_level/2+x_extrapix; l++){
      for(int m=-nsubpixels_per_level/2; m<nsubpixels_per_level/2+y_extrapix; m++){

	if((abs(l)+.5)<=subharmonic_layer_boundary_in_pixels && 
	   (abs(m)+.5)<=subharmonic_layer_boundary_in_pixels) {
	  subharmonic_pixel_weights[l+nsubpixels_per_level/2][m+nsubpixels_per_level/2] = 0;

	  /*
	  if(subharmonic_method::verbose_level)
	    cout << setw(15) << l << setw(15) << m 
		 << setw(15) << (abs(l)+.5)
		 << setw(15) << (abs(l)+.5)  
		 << setw(15) << subharmonic_layer_boundary_in_pixels << "\tzero weight\n";
	  */

       	} else if((abs(l)-.5)>=subharmonic_layer_boundary_in_pixels || 
		  (abs(m)-.5)>=subharmonic_layer_boundary_in_pixels){
	  subharmonic_pixel_weights[l+nsubpixels_per_level/2][m+nsubpixels_per_level/2] = 1;
	  /*
	  if(subharmonic_method::verbose_level)
	    cout << setw(15) << l << setw(15) << m 
		 << setw(15) << (abs(l)-.5)
		 << setw(15) << (abs(m)-.5)  
		 << setw(15) << subharmonic_layer_boundary_in_pixels << "\tunit weight\n";
	  */
	} else {
	  weight = 1;
	  tmp = subharmonic_layer_boundary_in_pixels - (abs(l) - .5);
	  if(tmp>=1) tmp=1;
	  if(tmp>=0) weight*=tmp;

	  tmp = subharmonic_layer_boundary_in_pixels - (abs(m) - .5);
	  if(tmp>=1) tmp=1;
	  if(tmp>=0) weight*=tmp;

	  subharmonic_pixel_weights[l+nsubpixels_per_level/2][m+nsubpixels_per_level/2] = (1-weight);
	  /*
	  if(subharmonic_method::verbose_level)
	    cout << setw(15) << l << setw(15) << m
		 << setw(15) << (abs(l)+.5)
		 << setw(15) << (abs(m)+.5)  
		 << setw(15) << subharmonic_layer_boundary_in_pixels 
		 << setw(15) << subharmonic_pixel_weights[l+nsubpixels_per_level/2][m+nsubpixels_per_level/2] << endl;
	  */
	}
      }
    }

    if(subharmonic_method::verbose_level){
      cout << "generalized_subharmonic_method::apply_subharmonic_correction weighting matrix:\n";
      for(int i=0; i<nsubpixels_per_level; i++){
	for(int j=0; j<nsubpixels_per_level; j++){
	  cout << setw(15) << subharmonic_pixel_weights[i][j];
	}
	cout << endl;
      }
    }


    // Downweight or zero pixels near zero spatial frequency that will
    // be partially or fully corrected using the subharmonic method.
    // In this subharmonic method, the zero spatial frequency pixel
    // always receives zero weight
    for(int i=axes[1]/2-nsubpixels_per_level/2; i<=axes[1]/2+nsubpixels_per_level/2; i++){
      for(int j=axes[0]/2-nsubpixels_per_level/2; j<=axes[0]/2+nsubpixels_per_level/2; j++){
	weight = subharmonic_pixel_weights[i-axes[1]/2+nsubpixels_per_level/2][j-axes[0]/2+nsubpixels_per_level/2];
	if(weight==1) continue;
	else if(weight==0)
	  data[2*(i*axes[0]+j)] = data[2*(i*axes[0]+j)+1] = 0;
	else {
	  data[2*(i*axes[0]+j)] *= weight;
	  data[2*(i*axes[0]+j)+1] *= weight;
	}
      }
    }

    // Now we do the subharmonics
    // NOTE: 
    // It seems that this works because the contribution
    // from each subharmonic level is decreasing, albeit slowly,
    // for a komolgorov spectrum.  More specifically, note that
    // a power law with exponent 11/6 is increasing slightly more
    // slowly than quadratically as the wavenumber decreases.
    // The subharmonic technique decreases the areal contribution
    // at each layer by a factor of 4 while halving the frequency.
    // Thus we get a decreasing contribution of order 2^{-1/6}
    // Consequently, this may not work for power laws with exponents
    // steeper than 2

    double pi_times_pixel_location[2], subpixel_location[2];
    double xfreq_interval = 2*M_PI/(pixscale*axes[1]);
    double yfreq_interval = 2*M_PI/(pixscale*axes[0]);
    double subpixel_weight = 1;
    double val, r1, r2, fac;
 
    // The factor to fix the units on the delta function
    // contribution to make it look like a smooth probability
    // distribution.
    // See Johansson & Gavel eq. 
    double norm;

    // Here if the zero frequency value of the power spectrum
    // is infinite, we will not fill it in on the final subharmonic
    // iteration.
    bool skip_zero_frequency = false;
    if(pspec.pole_at_zero_spatial_frequency()) skip_zero_frequency = true;


    for(int k=0; k<subharmonic_depth; k++){

      subpixel_weight /= (double) nsubpixels_per_pixel;

      norm = sqrt(xfreq_interval*yfreq_interval*subpixel_weight*subpixel_weight);

      for(int l=-nsubpixels_per_level/2; l<nsubpixels_per_level/2+x_extrapix; l++){
	for(int m=-nsubpixels_per_level/2; m<nsubpixels_per_level/2+y_extrapix; m++){

	  // Here we weight the contribution according to the subharmonic pixel weights
	  val = 1;
	  if(k<subharmonic_depth-1)
	    val = subharmonic_pixel_weights[l+nsubpixels_per_level/2][m+nsubpixels_per_level/2];

	  if(val==0 || (l==0 && m==0 && skip_zero_frequency)) continue; 

	  subpixel_location[0] = (l+x_halfpix)*xfreq_interval*subpixel_weight;
	  subpixel_location[1] = (m+y_halfpix)*yfreq_interval*subpixel_weight;

	  val = val*norm*
	    sqrt(pspec.value(sqrt(subpixel_location[0]*subpixel_location[0]+
				  subpixel_location[1]*subpixel_location[1])));


	  if(random) box_mueller(r1,r2);

	  int index;
	  for(int i=-axes[1]/2; i<axes[1]/2+x_extrapix; i++){
	    for(int j=-axes[0]/2; j<axes[0]/2+y_extrapix; j++){
	    
	      pi_times_pixel_location[0] = M_PI*(i + x_halfpix - (l+x_halfpix)*subpixel_weight);
	      pi_times_pixel_location[1] = M_PI*(j + y_halfpix - (m+y_halfpix)*subpixel_weight);

	      fac = val;

	      if(pi_times_pixel_location[0]!=0) fac *= sin(pi_times_pixel_location[0])/(pi_times_pixel_location[0]);
	      if(pi_times_pixel_location[1]!=0) fac *= sin(pi_times_pixel_location[1])/(pi_times_pixel_location[1]);
	    
	      index = (i+axes[1]/2)*axes[0] + j + axes[0]/2;
	      if(random){
		data[2*index] += r1*fac;
		data[2*index+1] += r2*fac;
	      } else {
		data[2*index] += fac*val;
		data[2*index+1] = 0;
	      } 
	    }
	  }
	}
      }
    }
  }

  void generalized_subharmonic_method::apply_subharmonic_correction(
  				const power_spectrum & pspec, 
				const vector<long> & axes, 
				double pixscale, 
				bool random,
				float * data) const {
    apply_subharmonic_correction_impl<float>(pspec, axes, pixscale, random, data);
  }

  void generalized_subharmonic_method::apply_subharmonic_correction(
  				const power_spectrum & pspec, 
				const vector<long> & axes, 
				double pixscale, 
				bool random,
				double * data) const {
    apply_subharmonic_correction_impl<double>(pspec, axes, pixscale, random, data);
  }

  namespace factory_register {
    const fits_keyval_set & get_Lane_subharmonic_method_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "Lane subharmonic method"));
      return *fkvs;
    }
    
    AO_sim_base * create_Lane_subharmonic_method(const iofits & iof) {
      return new Lane_subharmonic_method(iof);
    }
  }

  const bool Lane_subharmonic_method::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_Lane_subharmonic_method_keyval_set(), 
				      factory_register::create_Lane_subharmonic_method);

  void Lane_subharmonic_method::read(const char * filename) {
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "Lane_subharmonic_method::read - "
	   << "error opening file " << filename << endl;
      throw(string("Lane_subharmonic_method::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "Lane_subharmonic_method::read - "
	   << "error reading "
	   << this->unique_name() << " from file " 
	   << filename << endl;
      throw(string("Lane_subharmonic_method::read"));
    }
  }

  void Lane_subharmonic_method::read(const iofits & iof) {
    if(!iof.key_exists("TYPE")){
      cerr << "Lane_subharmonic_method::read error - "
	   << "unrecognized type of file\n";
      throw(string("Lane_subharmonic_method::read"));
    }
    string type, comment;
    iof.read_key("TYPE", type, comment);
    if(type!="Lane subharmonic method"){
      cerr << "Lane_subharmonic_method::read error - file of type " 
	   << type << " rather than of type circular aperture\n";
      throw(string("Lane_subharmonic_method::read"));
    }
    iof.read_key("SBHMDPTH", subharmonic_depth, comment);
    // leave iof pointing to the next header, if it exists
    if(iof.get_hdu_num()!=iof.get_num_hdus())
      iof.movrel_hdu(1);

  }

  void Lane_subharmonic_method::write(const char * filename) const {
    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "Lane_subharmonic_method::write - "
	   << "error opening file " << filename << endl;
      throw(string("Lane_subharmonic_method::write"));
    }

    try{this->write(iof);}
    catch(...){
      cerr << "Lane_subharmonic_method::write - "
	   << "error writing "
	   << this->unique_name() << " to file "
	   << filename << endl;
      throw(string("Lane_subharmonic_method::write"));
    }
  }

  void Lane_subharmonic_method::write(iofits & iof) const {
    this->subharmonic_method::write(iof);
    string type = "Lane subharmonic method";
    string comment = "object type";
    iof.write_key("TYPE", type, comment);
    comment = "subharmonic depth";
    iof.write_key("SBHMDPTH", subharmonic_depth, comment);
  }

  void Lane_subharmonic_method::print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << "Lane subharmonic method"
       << "/" << "object type" << endl;
    os << prefix << "SBHMDPTH   = " << setw(vlspc) << subharmonic_depth
       << "/" << "subharmonic depth" << endl;
  };
  
  void Johansson_Gavel_subharmonic_method::read(const char * filename) {
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "Johansson_Gavel_subharmonic_method::read - "
	   << "error opening file " << filename << endl;
      throw(string("Johansson_Gavel_subharmonic_method::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "Johansson_Gavel_subharmonic_method::read - "
	   << "error reading "
	   << this->unique_name() << " from file " 
	   << filename << endl;
      throw(string("Johansson_Gavel_subharmonic_method::read"));
    }
  }

  namespace factory_register {
    const fits_keyval_set & get_Johansson_Gavel_subharmonic_method_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "Johansson and Gavel subharmonic method"));
      return *fkvs;
    }
    
    AO_sim_base * create_Johansson_Gavel_subharmonic_method(const iofits & iof) {
      return new Johansson_Gavel_subharmonic_method(iof);
    }
  }

  const bool Johansson_Gavel_subharmonic_method::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_Johansson_Gavel_subharmonic_method_keyval_set(), 
				      factory_register::create_Johansson_Gavel_subharmonic_method);
  
  void Johansson_Gavel_subharmonic_method::read(const iofits & iof) {
    if(!iof.key_exists("TYPE")){
      cerr << "Johansson_Gavel_subharmonic_method::read error - "
	   << "unrecognized type of file\n";
      throw(string("Johansson_Gavel_subharmonic_method::read"));
    }
    string type, comment;
    iof.read_key("TYPE", type, comment);
    if(type!="Johansson and Gavel subharmonic method"){
      cerr << "Johansson_Gavel_subharmonic_method::read error - file of type " 
	   << type << " rather than of type circular aperture\n";
      throw(string("Johansson_Gavel_subharmonic_method::read"));
    }
    iof.read_key("SBHMDPTH", subharmonic_depth, comment);
    // leave iof pointing to the next header, if it exists
    if(iof.get_hdu_num()!=iof.get_num_hdus())
      iof.movrel_hdu(1);

  }

  void Johansson_Gavel_subharmonic_method::write(const char * filename) const {
    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "Johansson_Gavel_subharmonic_method::write - "
	   << "error opening file " << filename << endl;
      throw(string("Johansson_Gavel_subharmonic_method::write"));
    }

    try{this->write(iof);}
    catch(...){
      cerr << "Johansson_Gavel_subharmonic_method::write - "
	   << "error writing "
	   << this->unique_name() << " to file "
	   << filename << endl;
      throw(string("Johansson_Gavel_subharmonic_method::write"));
    }
  }

  void Johansson_Gavel_subharmonic_method::write(iofits & iof) const {
    this->subharmonic_method::write(iof);
    string type = "Johansson and Gavel subharmonic method";
    string comment = "object type";
    iof.write_key("TYPE", type, comment);
    comment = "subharmonic depth";
    iof.write_key("SBHMDPTH", subharmonic_depth, comment);
  }

  void Johansson_Gavel_subharmonic_method::print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << "Johansson and Gavel subharmonic method"
       << "/" << "object type" << endl;
    os << prefix << "SBHMDPTH   = " << setw(vlspc) << subharmonic_depth
       << "/" << "subharmonic depth" << endl;
  };
  
  namespace factory_register {
    const fits_keyval_set & get_quad_pixel_subharmonic_method_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "quad pixel subharmonic method"));
      return *fkvs;
    }
    
    AO_sim_base * create_quad_pixel_subharmonic_method(const iofits & iof) {
      return new quad_pixel_subharmonic_method(iof);
    }
  }

  const bool quad_pixel_subharmonic_method::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_quad_pixel_subharmonic_method_keyval_set(), 
				      factory_register::create_quad_pixel_subharmonic_method);

  void quad_pixel_subharmonic_method::read(const char * filename) {
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "quad_pixel_subharmonic_method::read - "
	   << "error opening file " << filename << endl;
      throw(string("quad_pixel_subharmonic_method::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "quad_pixel_subharmonic_method::read - "
	   << "error reading "
	   << this->unique_name() << " from file " 
	   << filename << endl;
      throw(string("quad_pixel_subharmonic_method::read"));
    }
  }

  void quad_pixel_subharmonic_method::read(const iofits & iof) {
    if(!iof.key_exists("TYPE")){
      cerr << "quad_pixel_subharmonic_method::read error - "
	   << "unrecognized type of file\n";
      throw(string("quad_pixel_subharmonic_method::read"));
    }
    string type, comment;
    iof.read_key("TYPE", type, comment);
    if(type!="quad pixel subharmonic method"){
      cerr << "quad_pixel_subharmonic_method::read error - file of type " 
	   << type << " rather than of type circular aperture\n";
      throw(string("quad_pixel_subharmonic_method::read"));
    }
    iof.read_key("SBHMDPTH", subharmonic_depth, comment);
    // leave iof pointing to the next header, if it exists
    if(iof.get_hdu_num()!=iof.get_num_hdus())
      iof.movrel_hdu(1);

  }

  void quad_pixel_subharmonic_method::write(const char * filename) const {
    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "quad_pixel_subharmonic_method::write - "
	   << "error opening file " << filename << endl;
      throw(string("quad_pixel_subharmonic_method::write"));
    }

    try{this->write(iof);}
    catch(...){
      cerr << "quad_pixel_subharmonic_method::write - "
	   << "error writing "
	   << this->unique_name() << " to file "
	   << filename << endl;
      throw(string("quad_pixel_subharmonic_method::write"));
    }
  }

  void quad_pixel_subharmonic_method::write(iofits & iof) const {
    this->subharmonic_method::write(iof);
    string type = "quad pixel subharmonic method";
    string comment = "object type";
    iof.write_key("TYPE", type, comment);
    comment = "subharmonic depth";
    iof.write_key("SBHMDPTH", subharmonic_depth, comment);
  }

  void quad_pixel_subharmonic_method::print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << "quad pixel subharmonic method"
       << "/" << "object type" << endl;
    os << prefix << "SBHMDPTH   = " << setw(vlspc) << subharmonic_depth
       << "/" << "subharmonic depth" << endl;
  };

  template<class T>
  void quad_pixel_subharmonic_method::apply_subharmonic_correction_impl(
  			const power_spectrum & pspec, 
			const vector<long> & axes, 
			double pixscale, 
			bool random,
			T * data) const {

    if(this->subharmonic_depth == 0) return;
    if(axes.size()!=2 || axes[0]<=0 || axes[1]<=0){
      cerr << "quad_pixel_subharmonic_method::apply_subharmonic_correction error - "
	   << "axes supplied to this function are poorly formed\n";
      cerr << "size " << axes.size() << " first " << axes[0] << " second " << axes[1] << endl;
      throw(string("quad_pixel_subharmonic_method::apply_subharmonic_correction"));
    }

    if(axes[0]%2==1 || axes[1]%2==1){
      cerr << "quad_pixel_subharmonic_method::apply_subharmonic_correction error - "
	   << "cannot perform correction on array with even dimensional axes: "
	   << axes[0] << " x " << axes[1] << endl;
      throw(string("quad_pixel_subharmonic_method::apply_subharmonic_correction"));
    }
      
    if(pixscale <= 0){
      cerr << "quad_pixel_subharmonic_method::apply_subharmonic_correction error - "
	   << "bad pixel scale " << pixscale << " supplied to this function\n";
      throw(string("quad_pixel_subharmonic_method::apply_subharmonic_correction"));
    }

    // Zero out the 4 central pixels, which will be treated instead by
    // the subharmonics.  
    int index;
    index = (axes[1]/2-1)*(axes[0]/2-1) + axes[0]/2 - 1;
    data[2*index] = data[2*index+1] = 0;
    index = (axes[1]/2-1)*(axes[0]/2) + axes[0]/2;
    data[2*index] = data[2*index+1] = 0;
    index = (axes[1]/2)*(axes[0]/2-1) + axes[0]/2 - 1;
    data[2*index] = data[2*index+1] = 0;
    index = (axes[1]/2)*(axes[0]/2) + axes[0]/2;
    data[2*index] = data[2*index+1] = 0;
    
    // Now we do the subharmonics
    // NOTE: 
    // It seems that this works because the contribution
    // from each subharmonic level is decreasing, albeit slowly,
    // for a komolgorov spectrum.  More specifically, note that
    // a power law with exponent 11/6 is increasing slightly more
    // slowly than quadratically as the wavenumber decreases.
    // The subharmonic technique decreases the areal contribution
    // at each layer by a factor of 4 while halving the frequency.
    // Thus we get a decreasing contribution of order 2^{-1/6}
    // Consequently, this may not work for power laws with exponents
    // steeper than 2

    // For even array dimensions, we will divide 2 pixels into subpixels.
    // In this case, we divide each subpixel in half.

    int n_x_pixels = 2, n_y_pixels = 2;
    int x_subpixels = 2, y_subpixels = 2;

    double pixel_location[2], subpixel_location[2];
    double xfreq_interval = 2*M_PI/(pixscale*axes[1]);
    double yfreq_interval = 2*M_PI/(pixscale*axes[0]);
    double xsubscale = 1;
    double ysubscale = 1;
    double val, r1, r2, fac;
 
    // the halfpixel information
    double x_halfpix=.5, y_halfpix=.5;
    int x_extrapix=0, y_extrapix=0;

    // The factor to fix the units on the delta function
    // contribution to make it look like a smooth probability
    // distribution.
    // See Johansson & Gavel 
    double norm;

    for(int k=0; k<subharmonic_depth; k++){
      xsubscale /= (double)(x_subpixels);
      ysubscale /= (double)(y_subpixels);
      norm = sqrt(xfreq_interval*yfreq_interval*xsubscale*ysubscale);

      for(int l=-(n_x_pixels*x_subpixels)/2; l<(n_x_pixels*x_subpixels)/2+x_extrapix; l++){
	for(int m=-(n_y_pixels*y_subpixels)/2; m<(n_y_pixels*y_subpixels)/2+y_extrapix; m++){

	  // If this is not the last subharmonic, we 
	  // don't want to fill in the central subharmonics
	  // because they will be filled in at the next subharmonic
	  // level.

	  if(k<subharmonic_depth-1 && (l==-1||l==0) && (m==-1||m==0)) continue;

	  // these are {-3/4, 3/4}, {-3/8, 3/8}, {-3/16, 3/16} ... (axis even)
	  subpixel_location[0] = (l+x_halfpix)*xfreq_interval*xsubscale;
	  subpixel_location[1] = (m+y_halfpix)*yfreq_interval*ysubscale;

	  val = norm*
	    sqrt(pspec.value(sqrt(subpixel_location[0]*subpixel_location[0]+
				  subpixel_location[1]*subpixel_location[1])));

	  if(random) box_mueller(r1,r2);

	  for(int i=-axes[1]/2; i<axes[1]/2+x_extrapix; i++){
	    for(int j=-axes[0]/2; j<axes[0]/2+x_extrapix; j++){
	    
	      pixel_location[0] = i + x_halfpix - (l+x_halfpix)*xsubscale;
	      pixel_location[1] = j + y_halfpix - (m+y_halfpix)*ysubscale;

	      fac = val;

	      if(pixel_location[0]!=0) fac *= sin(M_PI*pixel_location[0])/(M_PI*pixel_location[0]);
	      if(pixel_location[1]!=0) fac *= sin(M_PI*pixel_location[1])/(M_PI*pixel_location[1]);
	    
	      index = (i+axes[1]/2)*axes[0] + j + axes[0]/2;
	      if(random){
		data[2*index] += r1*fac;
		data[2*index+1] += r2*fac;
	      } else {
		data[2*index] += fac*val;
		data[2*index+1] = 0;
	      } 
	    }
	  }
	}
      }
    }
  }
  void quad_pixel_subharmonic_method::apply_subharmonic_correction(
  			const power_spectrum & pspec, 
			const vector<long> & axes, 
			double pixscale, 
			bool random,
			float * data) const {
    apply_subharmonic_correction_impl<float>(pspec, axes, pixscale, random, data);
  }
  void quad_pixel_subharmonic_method::apply_subharmonic_correction(
  			const power_spectrum & pspec, 
			const vector<long> & axes, 
			double pixscale, 
			bool random,
			double * data) const {
    apply_subharmonic_correction_impl<double>(pspec, axes, pixscale, random, data);
  }


  // The attempt to cut out most of the sinc calculation duplication
  /*
  void Lane_subharmonic_method::apply_subharmonic_correction(
  					const power_spectrum & pspec, 
					vector<long> axes, 
					double * data, 
					double pixscale) const {

    if(this->subharmonic_depth == 0) return;
    if(axes.size()!=2 || axes[0]<=0 || axes[1]<=0){
      cerr << "Lane_subharmonic_method::apply_subharmonic_correction error - "
	   << "axes supplied to this function are poorly formed\n";
      cerr << "size " << axes.size() << " first " << axes[0]
           << " second " << axes[1] << endl;
      throw(string("Lane_subharmonic_method::apply_subharmonic_correction"));
    }

    if(pixscale <= 0){
      cerr << "Lane_subharmonic_method::apply_subharmonic_correction error - "
	   << "bad pixel scale " << pixscale << " supplied to this function\n";
      throw(string("Lane_subharmonic_method::apply_subharmonic_correction"));
    }

    // Zero out the central pixels, which will be treated instead by
    // the subharmonics.  There are 1, 2, or 4 central pixels,
    // depending respectively on whether both axes are odd, one of
    // them is odd, or neither is odd.
    int index;
    if(axes[1]%2){
      if(axes[0]%2){
	index = (axes[1]/2)*(axes[0]/2) + axes[0]/2;
	data[2*index] = data[2*index+1] = 0;
      } else {
	index = (axes[1]/2)*(axes[0]/2-1) + axes[0]/2 - 1;
	data[2*index] = data[2*index+1] = 0;
	index = (axes[1]/2)*(axes[0]/2) + axes[0]/2;
	data[2*index] = data[2*index+1] = 0;
      }
    } else {
      if(axes[0]%2){
	index = (axes[1]/2-1)*(axes[0]/2) + axes[0]/2;
	data[2*index] = data[2*index+1] = 0;
	index = (axes[1]/2)*(axes[0]/2) + axes[0]/2;
	data[2*index] = data[2*index+1] = 0;
      } else {
	index = (axes[1]/2-1)*(axes[0]/2-1) + axes[0]/2 - 1;
	data[2*index] = data[2*index+1] = 0;
	index = (axes[1]/2-1)*(axes[0]/2) + axes[0]/2;
	data[2*index] = data[2*index+1] = 0;
	index = (axes[1]/2)*(axes[0]/2-1) + axes[0]/2 - 1;
	data[2*index] = data[2*index+1] = 0;
	index = (axes[1]/2)*(axes[0]/2) + axes[0]/2;
	data[2*index] = data[2*index+1] = 0;
      }
    }
    
    // Now we do the subharmonics
    // NOTE: 
    // It seems that this works because the contribution
    // from each subharmonic level is decreasing, albeit slowly,
    // for a komolgorov spectrum.  More specifically, note that
    // a power law with exponent 11/6 is increasing slightly more
    // slowly than quadratically as the wavenumber decreases.
    // The subharmonic technique decreases the areal contribution
    // at each layer by a factor of 4 while halving the frequency.
    // Thus we get a decreasing contribution of order 2^{-1/6}
    // Consequently, this may not work for power laws with exponents
    // steeper than 2

    // For even array dimensions, we will divide 2 pixels into subpixels.
    // In this case, we divide each subpixel in half.
    // For odd array dimensions, we will divide 1 pixel into subpixels
    // In this case, we divide each subpixel into thirds.
    // 
    // Thus, there are three possibilities:
    // both axes even results in dividing 4 pixels into 16 subpixels
    // both axes odd results in dividing 1 pixel into 9 subpixels
    // one even and one odd axis results in dividing 2 pixels into 12 subpixels

    int x_subpixels = 4, y_subpixels = 4;
    if(axes[1]%2) x_subpixels = 3;
    if(axes[0]%2) y_subpixels = 3;

    double pixel_location[2], subpixel_location[2];
    double subpixel_weight = 1;
    double x_subpixel_scale = 1/(double)x_subpixels;
    double y_subpixel_scale = 1/(double)y_subpixels;
    double xfreq_interval = 1/(pixscale*axes[1]);
    double yfreq_interval = 1/(pixscale*axes[0]);
    double r1, r2, fac;
 
    // the halfpixel information
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

    vector<double> vals(1, x_subpixels*y_subpixels);
    double * precomputed_sinc;
    alloc_size sz(axes[0], axes[1]);
    try{
      //precomputed_sinc = new double[axes[0]*axes[1]];
      precomputed_sinc = new double[sz];
    } catch(...) {
      cerr << "subharmonic_method::apply_subharmonic_correction error - "
	   << "unable to allocate memory for sinc arrays\n";
      throw;
    }

    // The factor to fix the units on the delta function
    // contribution to make it look like a smooth probability
    // distribution.
    // See Johansson & Gavel
    double norm = 1/pow(2*M_PI, 5/6.0)/pixscale/sqrt((float)(axes[0]*axes[1]));

    for(int k=0; k<subharmonic_depth; k++){

      cout << "subharmonic level " << k << endl;

      subpixel_weight /= (x_subpixels*y_subpixels);

      // Fill in the vals array, placing zeroes in the subharmonic locations
      // that we don't want to compute - either because they are part of a
      // subsequent subharmonic level or because the nominal power law value
      // is infinite.  This is the ugly part
      int elements_remaining = 0;
      for(int l=0; l<x_subpixels; l++){
	for(int m=0; m<y_subpixels; m++){
	  if(x_subpixels==3 && y_subpixels==3 && l==1 && m==1) 
	    vals[l*y_subpixels+m] = 0; 
	  else if(k<subharmonic_depth-1 && 
		  x_subpixels==4 && y_subpixels==4 && 
		  (l==1||l==2) && (m==1||m==2)) 
	    vals[l*y_subpixels+m] = 0;
	  else if(k<subharmonic_depth-1 && 
		  x_subpixels==4 && y_subpixels==3 && 
		  (l==1||l==2) && m==1) 
	    vals[l*y_subpixels+m] = 0;
	  else if(k<subharmonic_depth-1 && 
		  x_subpixels==3 && y_subpixels==4 && 
		  l==1 && (m==1||m==2)) 
	    vals[l*y_subpixels+m] = 0;
	  else {
	    elements_remaining++;
	    vals[l*y_subpixels+m] = 1;
	  }
	}
      }

      cout << "total elements " << elements_remaining << endl;

      while(elements_remaining){

      cout << "elements remaining " << elements_remaining << endl;

	int xindex, yindex;
	for(int l=0; l<vals.size(); l++)
	  if(vals[l]!=0){
	    xindex = l/y_subpixels;
	    yindex = l%y_subpixels;
	    break;
	  }

	// precalculate the sinc array
	cout << "precomputing sinc\n";
	vector<double> subpixel_location(2);
	subpixel_location[0] = (xindex-x_subpixels/2+x_halfpix)*x_subpixel_scale*xfreq_interval;
	subpixel_location[1] = (yindex-y_subpixels/2+y_halfpix)*y_subpixel_scale*yfreq_interval;

	for(int i=-axes[1]/2; i<axes[1]/2+x_extrapix; i++){
	  for(int j=-axes[0]/2; j<axes[0]/2+x_extrapix; j++){
	    pixel_location[0] = i + x_halfpix - subpixel_location[0];
	    pixel_location[1] = j + y_halfpix - subpixel_location[1];
	    precomputed_sinc[(i+axes[1]/2)*axes[0] + j + axes[0]/2] = 
	      sin(M_PI*pixel_location[0])/(M_PI*pixel_location[0])*
	      sin(M_PI*pixel_location[1])/(M_PI*pixel_location[1]);
	  }
	}

	fac = sqrt(pspec.value(sqrt(subpixel_location[0]*subpixel_location[0]+
				    subpixel_location[1]*subpixel_location[1]))) * subpixel_weight*norm;

	// run over subharmonics looking for ones that have some
	// permutation of the subpixel locations we have prepared
	vector<double> tmp_subpixel_location(2);
	long index2;
	for(int l=0; l<x_subpixels; l++){
	  for(int m=0; m<y_subpixels; m++){
	    tmp_subpixel_location[0] = (l-x_subpixels/2+x_halfpix)*x_subpixel_scale*xfreq_interval;
	    tmp_subpixel_location[1] = (m-y_subpixels/2+y_halfpix)*y_subpixel_scale*yfreq_interval;


	    if(fabs(tmp_subpixel_location[0]-subpixel_location[0]) < three_frame::precision &&
	       fabs(tmp_subpixel_location[1]-subpixel_location[1]) < three_frame::precision){
	      cout << "found match - a\n";
	      box_mueller(r1,r2);
	      for(int i=0; i<axes[1]; i++){
		for(int j=0; j<axes[0]; j++){
		  index = (i*axes[0]+j);
		  data[2*index] += r1*fac*precomputed_sinc[index];
		  data[2*index+1] += r2*fac*precomputed_sinc[index];
		}
	      }
	      vals[l*y_subpixels+m] = 0;
	      elements_remaining--;
	    } else if(fabs(tmp_subpixel_location[0]-subpixel_location[0]) < three_frame::precision &&
		      fabs(tmp_subpixel_location[1]+subpixel_location[1]) < three_frame::precision){
	      cout << "found match - b\n";
	      box_mueller(r1,r2);
	      for(int i=0; i<axes[1]; i++){
		for(int j=0; j<axes[0]; j++){
		  index = (i*axes[0]+j);
		  index2 = (axes[1]-i)*axes[0]+j;
		  data[2*index] += r1*fac*precomputed_sinc[index2];
		  data[2*index+1] += r2*fac*precomputed_sinc[index2];
		}
	      }
	      vals[l*y_subpixels+m] = 0;
	      elements_remaining--;
	    } else if(fabs(tmp_subpixel_location[0]+subpixel_location[0]) < three_frame::precision &&
		      fabs(tmp_subpixel_location[1]+subpixel_location[1]) < three_frame::precision){
	      cout << "found match - c\n";
	      box_mueller(r1,r2);
	      for(int i=0; i<axes[1]; i++){
		for(int j=0; j<axes[0]; j++){
		  index = (i*axes[0]+j);
		  index2 = (axes[1]-i)*axes[0]+(axes[0]-j);
		  data[2*index] += r1*fac*precomputed_sinc[index2];
		  data[2*index+1] += r2*fac*precomputed_sinc[index2];
		}
	      }
	      vals[l*y_subpixels+m] = 0;
	      elements_remaining--;
	    } else if(fabs(tmp_subpixel_location[0]+subpixel_location[0]) < three_frame::precision &&
		      fabs(tmp_subpixel_location[1]-subpixel_location[1]) < three_frame::precision){
	      cout << "found match - d\n";
	      box_mueller(r1,r2);
	      for(int i=0; i<axes[1]; i++){
		for(int j=0; j<axes[0]; j++){
		  index = (i*axes[0]+j);
		  index2 = i*axes[0]+(axes[0]-j);
		  data[2*index] += r1*fac*precomputed_sinc[index2];
		  data[2*index+1] += r2*fac*precomputed_sinc[index2];
		}
	      }
	      vals[l*y_subpixels+m] = 0;
	      elements_remaining--;
	    } else if(fabs(tmp_subpixel_location[0]-subpixel_location[1]) < three_frame::precision &&
		      fabs(tmp_subpixel_location[1]-subpixel_location[0]) < three_frame::precision){
	      cout << "found match - e\n";
	      box_mueller(r1,r2);
	      for(int i=0; i<axes[1]; i++){
		for(int j=0; j<axes[0]; j++){
		  index = (i*axes[0]+j);
		  index2 = j*axes[1]+i;
		  data[2*index] += r1*fac*precomputed_sinc[index2];
		  data[2*index+1] += r2*fac*precomputed_sinc[index2];
		}
	      }
	      vals[l*y_subpixels+m] = 0;
	      elements_remaining--;
	    } else if(fabs(tmp_subpixel_location[0]-subpixel_location[1]) < three_frame::precision &&
		      fabs(tmp_subpixel_location[1]+subpixel_location[0]) < three_frame::precision){
	      box_mueller(r1,r2);
	      cout << "found match - f\n";
	      for(int i=0; i<axes[1]; i++){
		for(int j=0; j<axes[0]; j++){
		  index = (i*axes[0]+j);
		  index2 = (axes[0]-j)*axes[1]+i;
		  data[2*index] += r1*fac*precomputed_sinc[index2];
		  data[2*index+1] += r2*fac*precomputed_sinc[index2];
		}
	      }
	      vals[l*y_subpixels+m] = 0;
	      elements_remaining--;
	    } else if(fabs(tmp_subpixel_location[0]+subpixel_location[1]) < three_frame::precision &&
		      fabs(tmp_subpixel_location[1]+subpixel_location[0]) < three_frame::precision){
	      cout << "found match - g\n";
	      box_mueller(r1,r2);
	      for(int i=0; i<axes[1]; i++){
		for(int j=0; j<axes[0]; j++){
		  index = (i*axes[0]+j);
		  index2 = j*axes[1]+(axes[0]-i);
		  data[2*index] += r1*fac*precomputed_sinc[index2];
		  data[2*index+1] += r2*fac*precomputed_sinc[index2];
		}
	      }
	      vals[l*y_subpixels+m] = 0;
	      elements_remaining--;
	    } else if(fabs(tmp_subpixel_location[0]+subpixel_location[1]) < three_frame::precision &&
		      fabs(tmp_subpixel_location[1]-subpixel_location[0]) < three_frame::precision){
	      cout << "found match - h\n";
	      box_mueller(r1,r2);
	      for(int i=0; i<axes[1]; i++){
		for(int j=0; j<axes[0]; j++){
		  index = (i*axes[0]+j);
		  index2 = (axes[0]-j)*axes[1]+(axes[0]-i);
		  data[2*index] += r1*fac*precomputed_sinc[index2];
		  data[2*index+1] += r2*fac*precomputed_sinc[index2];
		}
	      }
	      vals[l*y_subpixels+m] = 0;
	      elements_remaining--;
	    }
	  }
	}
      }
      if(axes[1]%2) x_subpixel_scale /= 3.0;
      else x_subpixel_scale /= 2.0;
      if(axes[0]%2) y_subpixel_scale /= 3.0;
      else y_subpixel_scale /= 2.0;
    }
    delete [] precomputed_sinc;
  }
  */

  // Here I tried to convert to a complex to real transform, and failed.
  /*
  void Lane_subharmonic_method::apply_subharmonic_correction(const power_spectrum & pspec, 
							     vector<long> complex_axes, 
							     double * data, 
							     double pixscale) const {

    if(this->subharmonic_depth == 0) return;
    if(complex_axes.size()!=2 || complex_axes[0]<=0 || complex_axes[1]<=0){
      cerr << "Lane_subharmonic_method::apply_subharmonic_correction error - "
	   << "complex_axes supplied to this function are poorly formed\n";
      cerr << "size " << complex_axes.size() << " first " << complex_axes[0] << " second " << complex_axes[1] << endl;
      throw(string("Lane_subharmonic_method::apply_subharmonic_correction"));
    }

    if(pixscale <= 0){
      cerr << "Lane_subharmonic_method::apply_subharmonic_correction error - "
	   << "bad pixel scale " << pixscale << " supplied to this function\n";
      throw(string("Lane_subharmonic_method::apply_subharmonic_correction"));
    }

    // the halfpixel information
    double x_halfpix=0, y_halfpix=0;
    int x_extrapix=1, y_extrapix=1;
    if(complex_axes[1]%2==0){
      x_halfpix = .5;
      x_extrapix = 0;
    }
    if(complex_axes[0]%2==0){
      y_halfpix = .5;
      y_extrapix = 0;
    }

    // Zero out the central pixels, which will be treated instead by
    // the subharmonics.  There are 1 or 2 central pixels.  The former
    // case occurs when the non-shortened axis of the halfcomplex array
    // is odd.  The latter case occurs when this axis is even
    int index;
    if(complex_axes[1]%2){
      index = (complex_axes[1]/2)*(complex_axes[0]/2+1);
      data[2*index] = data[2*index+1] = 0;
    } else {
      index = (complex_axes[1]/2)*(complex_axes[0]/2+1);
      data[2*index] = data[2*index+1] = 0;
      index = (complex_axes[1]/2+1)*(complex_axes[0]/2+1);
      data[2*index] = data[2*index+1] = 0;
    }

    // now we do the subharmonics
    // NOTE: 
    // It seems that this works because the contribution
    // from each subharmonic level is decreasing, albeit slowly,
    // for a komolgorov spectrum.  More specifically, note that
    // a power law with exponent 11/6 is increasing slightly more
    // slowly than quadratically as the wavenumber decreases.
    // The subharmonic technique decreases the areal contribution
    // at each layer by a factor of 4 while halving the frequency.
    // Thus we get a decreasing contribution of order 2^{-1/6}
    // Consequently, this may not work for power laws with exponents
    // steeper than 2

    // For even array dimensions, we will divide 2 pixels into subpixels.
    // In this case, we divide each subpixel in half.
    // For odd array dimensions, we will divide 1 pixel into subpixels
    // In this case, we divide each subpixel into thirds.
    // 
    // Thus, there are three possibilities:
    // both complex_axes even results in dividing 4 pixels into 16 subpixels
    // both complex_axes odd results in dividing 1 pixel into 9 subpixels
    // one even and one odd axis results in dividing 2 pixels into 12 subpixels

    int x_subpixels = 4, y_subpixels = 4;
    if(complex_axes[1]%2) x_subpixels = 3;
    if(complex_axes[0]%2) y_subpixels = 3;

    double pixel_location[2], subpixel_location[2];
    double subpixel_weight = 1;
    double x_subpixel_scale = 1/(double)x_subpixels;
    double y_subpixel_scale = 1/(double)y_subpixels;
    double xfreq_interval = 1/(pixscale*complex_axes[1]);
    double yfreq_interval = 1/(pixscale*complex_axes[0]);
    double val, r1, r2, fac;
 
    // The factor to fix the units on the delta function
    // contribution to make it look like a smooth probability
    // distribution.
    // See Johansson & Gavel
    double norm = 1/pow(2*M_PI, 5/6.0)/pixscale/sqrt((float)(complex_axes[0]*complex_axes[1]));

    for(int k=0; k<subharmonic_depth; k++){
      subpixel_weight /= (x_subpixels*y_subpixels);

      for(int l=0; l<x_subpixels; l++){
	for(int m=0; m<y_subpixels; m++){

	  // There are two circumstances under which we want
	  // to skip this particular subpixel.
	  // First, if this is not the last subharmonic, we 
	  // don't want to fill in the central subharmonics
	  // because they will be filled in at the next subharmonic
	  // level.
	  // Second, we never want to fill in the central pixel
	  // if we have 3 subpixels in both the x and y dimension,
	  // since the nominal value is infinite.

	  if(k<subharmonic_depth-1){
	    if(x_subpixels==4 && y_subpixels==4 && (l==1||l==2) && (m==1||m==2)) continue;
	    if(x_subpixels==4 && y_subpixels==3 && (l==1||l==2) && m==1) continue;
	    if(x_subpixels==3 && y_subpixels==4 && l==1 && (m==1||m==2)) continue;
	  }
	  if(x_subpixels==3 && y_subpixels==3 && l==1 && m==1) continue; 

	  // these are {-3/4, 3/4}, {-3/8, 3/8}, {-3/16, 3/16} ... (axis even)
	  // or {-1/3, 1/3}, {-1/9, 1/9}, {-1/27, 1/27} ... (axis odd)
	  subpixel_location[0] = (l-x_subpixels/2+x_halfpix)*x_subpixel_scale*xfreq_interval;
	  subpixel_location[1] = (m-y_subpixels/2+y_halfpix)*y_subpixel_scale*yfreq_interval;

	  val = sqrt(pspec.value(sqrt(subpixel_location[0]*subpixel_location[0]+
				      subpixel_location[1]*subpixel_location[1])));

	  cout << "depth " << k+1 
	       << " subharmonic " << l << ", " << m 
	       << "\tlocation " << setprecision(5) << setw(8) << subpixel_location[0] << ", " 
	       << setprecision(5) << setw(8) << subpixel_location[1] 
	       << "\tval " << setprecision(5) << setw(8) << val 
	       << "\tnorm " << setprecision(5) << setw(8) << norm 
	       << "\tweight " << setprecision(5) << setw(8) << subpixel_weight
	       << "\twtd val "<< setprecision(5)  << setw(8) << val*norm*subpixel_weight
	       << endl;
	  
	  val *= subpixel_weight*norm;

	  box_mueller(r1,r2);

	  for(int i=-complex_axes[1]/2; i<complex_axes[1]/2+x_extrapix; i++){
	    for(int j=0; j<complex_axes[0]/2+1; j++){
	    
	      pixel_location[0] = i + x_halfpix - subpixel_location[0];
	      pixel_location[1] = j + y_halfpix - subpixel_location[1];

	      fac = val;

	      if(pixel_location[0]!=0) fac *= sin(M_PI*pixel_location[0])/(M_PI*pixel_location[0]);
	      if(pixel_location[1]!=0) fac *= sin(M_PI*pixel_location[1])/(M_PI*pixel_location[1]);
	    
	      index = (i+complex_axes[1]/2)*(complex_axes[0]/2+1) + j;
	      data[2*index] += r1*fac;
	      data[2*index+1] += r2*fac;
	    }
	  }
	}
      }
      if(complex_axes[1]%2) x_subpixel_scale /= 3.0;
      else x_subpixel_scale /= 2.0;
      if(complex_axes[0]%2) y_subpixel_scale /= 3.0;
      else y_subpixel_scale /= 2.0;
    }
  }
  */

}   
