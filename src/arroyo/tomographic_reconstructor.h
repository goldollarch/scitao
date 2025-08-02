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

#ifndef TOMOGRAPHIC_RECONSTRUCTOR_H
#define TOMOGRAPHIC_RECONSTRUCTOR_H

#include "AO_sim_base.h"
#include "pixel_array.h"
#include <map>
#include <vector>
#include "aperture.h"
#include "covariance.h"
#include "observation.h"

namespace Arroyo {

  ///
  /// A class to hold a pseudo open loop tomographic reconstructor
  ///

  template<class precision>
    class pseudo_open_loop_tomographic_reconstructor :
    public AO_sim_base {
    
    private:

    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("pseudo open loop tomographic reconstructor"));};

    protected:

    // guide star emitters
    vector<emitter*> guide_star_emitters;

    // emitter for which tomography is estimated
    emitter * target_emitter;

    // emitter for which the final covariance matrix below has
    // been computed.  This covariance matrix must be recomputed
    // when a different emitter is 
    mutable emitter * arg_emitter;

    // The four matrices below store OPD variances
    // (i.e. OPD^{2})

    // Reconstruction matrix of dimension N^{2}xN^{2}xngs
    // This matrix takes measured phases in the direction
    // of the guide stars and computes a phase in the 
    // direction of the target emitter
    precision * sigma_ba_inv_sigma_aa;

    // Covariance matrix of dimension N^{2}xN^{2}
    // This is the covariance matrix between the 
    // estimated phases for the target emitter
    precision * sigma_ba_inv_sigma_aa_sigma_ab;

    // Covariance matrix of dimension N^{2}xN^{2}
    // This is the covariance matrix between the 
    // estimated phases for the target emitter
    // and for the arg emitter 
    mutable precision * sigma_ba_inv_sigma_aa_sigma_ac;

    // Covariance matrix of dimension N^{2}xN^{2}
    // This is the covariance matrix for the
    // arg emitter 
    mutable precision * sigma_cc;

    // mask for pixel_array <=> matrix conversion
    mutable short * aperture_mask;
 
    // The circular aperture
    circular_aperture circ_ap;

    // The refractive atmospheric model
    refractive_atmospheric_model ref_atm_model_at_zenith;

    double pupil_plane_pixel_scale_meters;

    long nsteps_in_integration;

    mutable long nactuators;

    
    ///////////////////////////////////////////
    ///  Null constructor
    pseudo_open_loop_tomographic_reconstructor(){
      target_emitter = arg_emitter = NULL;
      sigma_ba_inv_sigma_aa = NULL;
      sigma_ba_inv_sigma_aa_sigma_ab = NULL;
      sigma_ba_inv_sigma_aa_sigma_ac = NULL;
      sigma_cc = NULL;
      aperture_mask = NULL;
    };


    ///////////////////////////////////////////
    ///  Get an array corresponding to the aperture mask
    ///  Memory must be freed by the calling routine
    void initialize_aperture_mask() const;


    ///////////////////////////////////////////
    ///  Initialize a correlation matrix
    void initialize_auto_correlation_matrix(const vector<emitter *> & emtrs,
					    precision * array) const;

    ///////////////////////////////////////////
    ///  Initialize a correlation matrix
    void initialize_cross_correlation_matrix(const emitter * emtr,
					     const vector<emitter *> & emtrs,
					     precision * array) const;

    ///////////////////////////////////////////
    ///  Invert sigma_aa by SVD
    void invert_sigma_aa_via_SVD(precision * sigma_aa);

    ///////////////////////////////////////////
    ///  Emitter initialization
    void initialize_emitter(const emitter & emtr) const;

    ///////////////////////////////////////////
    ///  private function supporting two member functions
    double private_aperture_averaged_phase_variance(emitter & emtr,
						    double wavelength_meters,
						    bool differential) const;

    ///////////////////////////////////////////
    ///  private function supporting two member functions
    pixel_array<precision> private_phase_variance(emitter & emtr,
						  double wavelength_meters,
						  bool differential) const;

    ///////////////////////////////////////////
    ///  private function supporting two member functions
    pixel_array<precision> private_phase_covariance(emitter & emtr,
						    int xindex,
						    int yindex,
						    double wavelength_meters,
						    bool differential) const;

    public:

    ///////////////////////////////////////////
    ///  Construct from the bits
    pseudo_open_loop_tomographic_reconstructor(const emitter & target,
					       const vector<emitter*> & guide_star_emitters,
					       const refractive_atmospheric_model & ref_atm_model_at_zenith,
					       double pupil_plane_pixel_scale_meters,
					       int nsteps_in_integration,
					       const circular_aperture & circ_ap,
					       bool perform_SVD = true);

    ///////////////////////////////////////////
    ///  Copy constructor
    pseudo_open_loop_tomographic_reconstructor(const pseudo_open_loop_tomographic_reconstructor & polc_tomo_recon);

    ///////////////////////////////////////////
    ///  Construct from file
    pseudo_open_loop_tomographic_reconstructor(const char * filename);

    ///////////////////////////////////////////
    ///  Construct from iofits object
    pseudo_open_loop_tomographic_reconstructor(const Arroyo::iofits & iof);

    ///////////////////////////////////////////
    ///  Virtual destructor
    ~pseudo_open_loop_tomographic_reconstructor();
    
    ///////////////////////////////////////////
    ///  Operator = 
    pseudo_open_loop_tomographic_reconstructor & 
      operator=(const pseudo_open_loop_tomographic_reconstructor & polc_tomo_recon);

    ///////////////////////////////////////////
    ///  Read from file
    void read(const char * filename);

    ///////////////////////////////////////////
    ///  Read from iofits
    void read(const Arroyo::iofits & iof);

    ///////////////////////////////////////////
    ///  Write to file
    void write(const char * filename) const;
 
    ///////////////////////////////////////////
    ///  Write to iofits
    void write(Arroyo::iofits & iof) const;
 
    ///////////////////////////////////////////
    ///  Print
    void print(std::ostream & os, const char * prefix="") const;

    ///////////////////////////////////////////
    ///  Get aperture averaged variance
    double aperture_averaged_phase_variance(emitter & emtr,
					    double wavelength_meters) const;

    ///////////////////////////////////////////
    ///  Get aperture averaged differential variance
    double aperture_averaged_differential_phase_variance(emitter & emtr,
							 double wavelength_meters) const;

    ///////////////////////////////////////////
    ///  Get variance
    pixel_array<precision> phase_variance(emitter & emtr,
					  double wavelength_meters) const;

    ///////////////////////////////////////////
    ///  Get variance
    pixel_array<precision> differential_phase_variance(emitter & emtr,
						       double wavelength_meters) const;

    ///////////////////////////////////////////
    ///  Get variance
    pixel_array<precision> phase_covariance(emitter & emtr,
					    int xindex,
					    int yindex,
					    double wavelength_meters) const;

    ///////////////////////////////////////////
    ///  Get variance
    pixel_array<precision> differential_phase_covariance(emitter & emtr,
							 int xindex,
							 int yindex,
							 double wavelength_meters) const;

    //////////////////////////////////////////////////////////////////
    ///  Returns a pixel array containing the value of the structure
    ///  function computed betweeen a point tp1 and all other points
    ///  in the aperture.  The structure function is computed on a
    ///  grid that covers the aperture with sampling set by the
    ///  argument pixel_scale_meters.
    ///
    ///  The point tp must lie within the aperture, or this function
    ///  throws an error
    pixel_array<precision> phase_structure_function(const emitter & emtr,
						    double wavelength_meters,
						    int xindex,
						    int yindex) const;
    
    //////////////////////////////////////////////////////////////////
    ///  Returns a pixel array containing the long exposure optical
    ///  transfer function
    basic_otf<precision> optical_transfer_function(const emitter & emtr,
						   double wavelength_meters) const;
    

    //////////////////////////////////////////////////////////////////
    ///  Returns a pixel array containing the long exposure point
    ///  spread function
    basic_observation<precision> point_spread_function(const emitter & emtr,
						       double wavelength_meters,
						       double field_size_arcsecs,
						       double oversampling_factor) const;
    ///////////////////////////////////////////
    ///  Reconstruct the zernike and zonal residuals from a
    ///  Shack Hartmann centroid class instance
    void reconstruct(const vector<Arroyo::pixel_array<double> > & measured_phases,
		     Arroyo::pixel_array<double> & commands) const;

    ///  Verbose level
    static int verbose_level;

  };
  
  template<class precision>
    int pseudo_open_loop_tomographic_reconstructor<precision>::verbose_level = 0;

  template<class precision>
    void pseudo_open_loop_tomographic_reconstructor<precision>::
    initialize_aperture_mask() const {

    if(aperture_mask!=NULL) return;

    vector<long> axes(2,(long)ceil(this->circ_ap.get_diameter()/this->pupil_plane_pixel_scale_meters));
    double normalization_factor = 2*this->pupil_plane_pixel_scale_meters/circ_ap.get_diameter();
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

    try{
      this->aperture_mask = new short[axes[0]*axes[1]];
    } catch(...) {
      cerr << "pseudo_open_loop_tomographic_reconstructor::initialize_aperture_mask - error allocating memory\n";
      throw(string("pseudo_open_loop_tomographic_reconstructor::initialize_aperture_mask"));
    }

    this->nactuators = 0;
    three_vector pixel_vector;
    int index;
    for(int m=-axes[1]/2; m<axes[1]/2+x_extrapix; m++){
      for(int n=-axes[0]/2; n<axes[0]/2+y_extrapix; n++){
	
	index = (m+axes[1]/2)*axes[0]+n+axes[0]/2;
	
	pixel_vector = three_vector((m+x_halfpix)*normalization_factor,
				    (n+y_halfpix)*normalization_factor,
				    0,
				    this->circ_ap);
	
	if(pixel_vector.length()<=1){
	  this->nactuators++;
	  this->aperture_mask[index] = 1;
	} else
	  this->aperture_mask[index] = 0;

      }
    }
  }


  template<class precision>
    void pseudo_open_loop_tomographic_reconstructor<precision>::
    initialize_auto_correlation_matrix(const vector<emitter *> & emtrs,
				       precision * array) const {
    

    this->initialize_aperture_mask();

    int index;
    int block_index, subblock_index;
    int inner_count, outer_count;
    int nstars = emtrs.size();
    double ref_wavelength_meters = 2*M_PI; // with this choice, we store OPD

    vector<long> axes(2,(long)ceil(this->circ_ap.get_diameter()/this->pupil_plane_pixel_scale_meters));

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

    pixel_array<precision> pixarr;
    for(int i=0; i<nstars; i++){
      for(int j=0; j<=i; j++){
	
	Arroyo::phase_covariance<precision, circular_aperture> phase_cvnce(*(emtrs[i]),
									   *(emtrs[j]),
									   this->ref_atm_model_at_zenith,
									   this->circ_ap);
      
	block_index = this->nactuators*this->nactuators*nstars*i+j*this->nactuators;
      
	outer_count = 0;
	for(int m=-axes[1]/2; m<axes[1]/2+x_extrapix; m++){
	  for(int n=-axes[0]/2; n<axes[0]/2+y_extrapix; n++){
	  
	    index = (m+axes[1]/2)*axes[0]+n+axes[0]/2;
	  
	    if(this->aperture_mask[index]==1){
	    
	      pixarr = phase_cvnce.covariance(this->pupil_plane_pixel_scale_meters,
					      ref_wavelength_meters,
					      this->nsteps_in_integration,
					      m,
					      n);
	    
	      if(pseudo_open_loop_tomographic_reconstructor<precision>::verbose_level>=3){
		iofits iof;
		stringstream ss;
		ss << "phase_covariance_aa_" << i << "_" << j << "_" << m << "_" << n << ".fits";
		iof.create(ss.str().c_str());
		iof.create_image(pixarr.get_axes(), array);
		pixarr.write(iof);
	      }
	    
	    
	      inner_count = 0;
	      for(int p=-axes[1]/2; p<axes[1]/2+x_extrapix; p++){
		for(int q=-axes[0]/2; q<axes[0]/2+y_extrapix; q++){
		
		  index = (p+axes[1]/2)*axes[0]+q+axes[0]/2;
		
		  if(this->aperture_mask[index]==1){
		    subblock_index = block_index + outer_count*nactuators*nstars + inner_count;
		    array[subblock_index] = pixarr.data(index);
		    if(i!=j){
		      subblock_index = nactuators*nactuators*nstars*j+i*nactuators +
			inner_count*nactuators*nstars + outer_count;
		      array[subblock_index] = pixarr.data(index);
		    }
		    inner_count++;
		  }
		}
	      }
	      outer_count++;
	    }
	  }
	}
      }
    }
  }

  template<class precision>
    void pseudo_open_loop_tomographic_reconstructor<precision>::
    initialize_cross_correlation_matrix(const emitter * emtr,
					const vector<emitter *> & emtrs,
					precision * array) const {
    
    this->initialize_aperture_mask();

    int block_index, subblock_index;
    int inner_count, outer_count;
    int index;
    double ref_wavelength_meters = 2*M_PI; // with this choice, we store OPD

    vector<long> axes(2,(long)ceil(this->circ_ap.get_diameter()/this->pupil_plane_pixel_scale_meters));
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

    pixel_array<precision> pixarr;
    for(int i=0; i<emtrs.size(); i++){
      Arroyo::phase_covariance<precision, circular_aperture> phase_cvnce(*(emtr),
									 *(emtrs[i]),
									 this->ref_atm_model_at_zenith,
									 this->circ_ap);
      
      block_index = this->nactuators*this->nactuators*i;
      
      outer_count = 0;
      for(int m=-axes[1]/2; m<axes[1]/2+x_extrapix; m++){
	for(int n=-axes[0]/2; n<axes[0]/2+y_extrapix; n++){
	  
	  index = (m+axes[1]/2)*axes[0]+n+axes[0]/2;
	  
	  if(this->aperture_mask[index]==1){

	    pixarr = phase_cvnce.covariance(this->pupil_plane_pixel_scale_meters,
					    ref_wavelength_meters,
					    this->nsteps_in_integration,
					    m,
					    n);

	    if(pseudo_open_loop_tomographic_reconstructor<precision>::verbose_level>=3){
	      iofits iof;
	      stringstream ss;
	      ss << "phase_covariance_ba_" << i << "_" << m << "_" << n << ".fits";
	      iof.create(ss.str().c_str());
	      iof.create_image(pixarr.get_axes(), array);
	      pixarr.write(iof);
	    }


	    inner_count = 0;
	    for(int p=-axes[1]/2; p<axes[1]/2+x_extrapix; p++){
	      for(int q=-axes[0]/2; q<axes[0]/2+y_extrapix; q++){

		index = (p+axes[1]/2)*axes[0]+q+axes[0]/2;
		    
		if(this->aperture_mask[index]==1){
		  subblock_index = block_index + outer_count*nactuators + inner_count;
		  array[subblock_index] = pixarr.data(index);
		  inner_count++;
		}
	      }
	    }
	    outer_count++;
	  }
	}
      }      
    }
  }

  
  template<class precision>
    void pseudo_open_loop_tomographic_reconstructor<precision>::
    invert_sigma_aa_via_SVD(precision * sigma_aa) {
    
    
    bool store_eigenmodes = false;
    if(pseudo_open_loop_tomographic_reconstructor<precision>::verbose_level>=3)
      store_eigenmodes = true;
    
    char jobu[1];
    char jobvt[1];
    if(store_eigenmodes)
      jobu[0] = jobvt[0] = 'A';
    else 
      jobu[0] = jobvt[0] = 'S';
    int lwork=-1, info;
    precision optimal_workspace_size;
    precision *singular_values, *g_gtranspose_eigenmode_data, *gtranspose_g_eigenmode_data, *work;

    vector<long> sigma_aa_axes(2,this->nactuators*this->guide_star_emitters.size());
    int sigma_aa_axes_0 = sigma_aa_axes[0];
    int sigma_aa_axes_1 = sigma_aa_axes[1];

    singular_value_decomposition<precision>(jobu, 
					    jobvt, 
					    sigma_aa_axes_0,
					    sigma_aa_axes_1, 
					    sigma_aa, 
					    sigma_aa_axes_0, 
					    singular_values, 
					    g_gtranspose_eigenmode_data, 
					    sigma_aa_axes_0, 
					    gtranspose_g_eigenmode_data, 
					    sigma_aa_axes_1,
					    &optimal_workspace_size, 
					    lwork, 
					    info);
    
    lwork = (int)optimal_workspace_size;

    try{
      singular_values = new precision[sigma_aa_axes[1]];
      if(store_eigenmodes)
	g_gtranspose_eigenmode_data = new precision[sigma_aa_axes[0]*sigma_aa_axes[0]];
      else
	g_gtranspose_eigenmode_data = new precision[sigma_aa_axes[0]*sigma_aa_axes[1]];
      gtranspose_g_eigenmode_data = new precision[sigma_aa_axes[1]*sigma_aa_axes[1]];
      work = new precision[(int)(optimal_workspace_size)];
    } catch(...) {
      cerr << "pseudo_open_loop_tomographic_reconstructor::pseudo_open_loop_tomographic_reconstructor error - "
	   << "unable to allocate memory to perform the singular value decomposition\n";
      throw(string("pseudo_open_loop_tomographic_reconstructor::pseudo_open_loop_tomographic_reconstructor"));
    }

    singular_value_decomposition<precision>(jobu, 
					    jobvt, 
					    sigma_aa_axes_0, 
					    sigma_aa_axes_1, 
					    sigma_aa, 
					    sigma_aa_axes_0, 
					    singular_values, 
					    g_gtranspose_eigenmode_data, 
					    sigma_aa_axes_0, 
					    gtranspose_g_eigenmode_data, 
					    sigma_aa_axes_1,
					    work, 
					    lwork, 
					    info);
    
    delete [] work;

    if(pseudo_open_loop_tomographic_reconstructor<precision>::verbose_level>=3){
      {
	iofits iof;
	iof.create("sigma_aa_singular_values.fits");
	iof.create_image(vector<long>(1,sigma_aa_axes[1]), 
			 singular_values);
	iof.write_image(0, 
			sigma_aa_axes[1]-1, 
			singular_values);
      }
      {
	iofits iof;
	iof.create("sigma_aa_g_gtranspose.fits");
	iof.create_image(sigma_aa_axes, 
			 g_gtranspose_eigenmode_data);
	iof.write_image(0, 
			sigma_aa_axes[0]*sigma_aa_axes[1]-1, 
			g_gtranspose_eigenmode_data);
      }
      {
	iofits iof;
	iof.create("work_sigma_aa.fits");
	iof.create_image(sigma_aa_axes, 
			 sigma_aa);
	iof.write_image(0, 
			sigma_aa_axes[0]*sigma_aa_axes[1]-1, 
			sigma_aa);
      }
    }

    // Form the pseudoinverse from V S^{+} U^{T}
    double tmp;
    int index;
    for(int i=0; i<sigma_aa_axes[1]; i++){
      for(int j=0; j<sigma_aa_axes[0]; j++){
	tmp = 0;
	index = i*sigma_aa_axes[1];
	for(int k=0; k<sigma_aa_axes[1]; k++)
	  tmp += (gtranspose_g_eigenmode_data[index+k])*
	    (g_gtranspose_eigenmode_data[j+k*sigma_aa_axes[0]])/(singular_values[k]);

	sigma_aa[i*sigma_aa_axes[0]+j]=tmp;
      }
    }

    if(pseudo_open_loop_tomographic_reconstructor<precision>::verbose_level){
      {
	iofits iof;
	iof.create("inv_sigma_aa.fits");
	iof.create_image(sigma_aa_axes, 
			 sigma_aa);
	iof.write_image(0, 
			sigma_aa_axes[0]*sigma_aa_axes[1]-1, 
			sigma_aa);
      }
    }
      
    delete [] singular_values;
    delete [] g_gtranspose_eigenmode_data;
    delete [] gtranspose_g_eigenmode_data;
  }

  template<class precision>
    void pseudo_open_loop_tomographic_reconstructor<precision>::
    initialize_emitter(const emitter & emtr) const {
    
    const plane_wave_emitter *tmp_e1, *tmp_e2;

    if(this->arg_emitter!=NULL &&
       (tmp_e1=dynamic_cast<const plane_wave_emitter *>(this->arg_emitter)) &&
       (tmp_e2=dynamic_cast<const plane_wave_emitter *>(&emtr))){

      if(operator==(*tmp_e1,*tmp_e2)){
	return;
      }
    }

    if(this->arg_emitter!=NULL)
      delete this->arg_emitter;

    this->arg_emitter = emitter::emitter_factory(&emtr);
    
    if(this->sigma_ba_inv_sigma_aa_sigma_ac==NULL){
      try{
	this->sigma_ba_inv_sigma_aa_sigma_ac = 
	  new precision[this->nactuators*this->nactuators];
      } catch(...) {
	cerr << "pseudo_open_loop_tomographic_reconstructor::initialize_emitter "
	     << " error allocating memory\n"; 
	throw(string("pseudo_open_loop_tomographic_reconstructor::initialize_emitter"));
      }
    }  

    if(this->sigma_cc==NULL){
      try{
	this->sigma_cc = 
	  new precision[this->nactuators*this->nactuators];
      } catch(...) {
	cerr << "pseudo_open_loop_tomographic_reconstructor::initialize_emitter "
	     << " error allocating memory\n"; 
	throw(string("pseudo_open_loop_tomographic_reconstructor::initialize_emitter"));
      }
    }  


    //////////////////////////////////
    // Allocate memory for sigma_ca //
    //////////////////////////////////
    precision *sigma_ca;
    vector<long> sigma_ca_axes(2, this->nactuators*this->guide_star_emitters.size());
    sigma_ca_axes[0] = this->nactuators;
    try{
      sigma_ca = new precision[sigma_ca_axes[0]*sigma_ca_axes[1]];
    } catch(...) {
      cerr << "pseudo_open_loop_tomographic_reconstructor::initialize_emitter\n"
	   << " error allocating memory for covariance arrays\n";
      throw(string("pseudo_open_loop_tomographic_reconstructor::initialize_emitter"));
    }


    /////////////////////////
    // Initialize sigma_ca //
    /////////////////////////
    try{

      clock_t ci = clock();

      this->initialize_cross_correlation_matrix(this->arg_emitter,
						this->guide_star_emitters,
						sigma_ca);

      if(pseudo_open_loop_tomographic_reconstructor<precision>::verbose_level){
	cout << "pseudo_open_loop_tomographic_reconstructor::initialize_emitter\n"
	     << "\tformed sigma_ca in " 
	     << (clock()-ci)/(double)CLOCKS_PER_SEC 
	     << endl;
	iofits iof;
	iof.create("sigma_ca.fits");
	iof.create_image(sigma_ca_axes, 
			 sigma_ca);
	iof.write_image(0, 
			sigma_ca_axes[0]*sigma_ca_axes[1]-1, 
			sigma_ca);
      }
    } catch(...) {
      cerr << "pseudo_open_loop_tomographic_reconstructor::initialize_emitter\n"
	   << "error forming sigma_ca\n";
      throw(string("pseudo_open_loop_tomographic_reconstructor::initialize_emitter"));
    }

    ///////////////////////////////////////////////
    // Initialize sigma_ba_inv_sigma_aa_sigma_ac //
    ///////////////////////////////////////////////
    try{
      clock_t ci = clock();
      double tmp;
      for(int i=0; i<this->nactuators; i++){
	for(int j=0; j<this->nactuators; j++){
	  tmp = 0;
	  for(int k=0; k<sigma_ca_axes[1]; k++){
	    tmp += this->sigma_ba_inv_sigma_aa[k*sigma_ca_axes[0]+i]*
	      sigma_ca[k*sigma_ca_axes[0]+j];
	  }

	  this->sigma_ba_inv_sigma_aa_sigma_ac[i*this->nactuators+j]=tmp;
	}
      }

      if(pseudo_open_loop_tomographic_reconstructor<precision>::verbose_level){
	cout << "pseudo_open_loop_tomographic_reconstructor::initialize_emitter\n"
	     << "\tformed sigma_ba_inv_sigma_aa_sigma_ac in " 
	     << (clock()-ci)/(double)CLOCKS_PER_SEC 
	     << endl;
	iofits iof;
	iof.create("sigma_ba_inv_sigma_aa_sigma_ac.fits");
	iof.create_image(vector<long>(2,this->nactuators), 
			 this->sigma_ba_inv_sigma_aa_sigma_ac);
	iof.write_image(0,
			this->nactuators*this->nactuators-1, 
			this->sigma_ba_inv_sigma_aa_sigma_ac);
      }
    } catch(...) {
      cerr << "pseudo_open_loop_tomographic_reconstructor::initialize_emitter\n"
	   << "error forming sigma_ba_inv_sigma_aa_sigma_ac\n";
      throw(string("pseudo_open_loop_tomographic_reconstructor::initialize_emitter"));
    }

    delete [] sigma_ca;

    /////////////////////////
    // Initialize sigma_cc //
    /////////////////////////
    try{
      clock_t ci = clock();

      vector<long> sigma_cc_axes(2, this->nactuators);

      this->initialize_auto_correlation_matrix(vector<emitter *>(1,this->arg_emitter),
					       sigma_cc);

      if(pseudo_open_loop_tomographic_reconstructor<precision>::verbose_level){
	cout << "pseudo_open_loop_tomographic_reconstructor::initialize_emitter\n"
	     << "\tformed sigma_cc in " 
	     << (clock()-ci)/(double)CLOCKS_PER_SEC 
	     << endl;
	iofits iof;
	iof.create("sigma_cc.fits");
	iof.create_image(vector<long>(2,sigma_cc_axes[0]), 
			 this->sigma_cc);
	iof.write_image(0,
			sigma_cc_axes[0]*sigma_cc_axes[0]-1, 
			this->sigma_cc);
      }
    } catch(...) {
      cerr << "pseudo_open_loop_tomographic_reconstructor::initialize_emitter\n"
	   << "error forming sigma_cc\n";
      throw(string("pseudo_open_loop_tomographic_reconstructor::initialize_emitter"));
    }

    if(pseudo_open_loop_tomographic_reconstructor<precision>::verbose_level){
      vector<long> sigma_cc_axes(2, this->nactuators);
      int nelem = sigma_cc_axes[0]*sigma_cc_axes[1];
      precision * sigma_delta_delta = new precision[nelem];
      for(int i=0; i<nelem; i++){
	sigma_delta_delta[i] = this->sigma_cc[i] +
	  this->sigma_ba_inv_sigma_aa_sigma_ab[i] -
	  2*this->sigma_ba_inv_sigma_aa_sigma_ac[i];
      }
      iofits iof;
      iof.create("sigma_delta_delta.fits");
      iof.create_image(vector<long>(2,sigma_cc_axes[0]), 
		       sigma_delta_delta);
      iof.write_image(0,
		      sigma_cc_axes[0]*sigma_cc_axes[0]-1, 
		      sigma_delta_delta);
      delete [] sigma_delta_delta;
    }
  }

  template<typename precision>
    pseudo_open_loop_tomographic_reconstructor<precision>::
    pseudo_open_loop_tomographic_reconstructor(const pseudo_open_loop_tomographic_reconstructor & polc_tomo_recon){
    this->target_emitter = arg_emitter = NULL;
    this->sigma_ba_inv_sigma_aa = NULL;
    this->sigma_ba_inv_sigma_aa_sigma_ab = NULL;
    this->sigma_ba_inv_sigma_aa_sigma_ac = NULL;
    this->sigma_cc = NULL;
    this->aperture_mask = NULL;
    this->operator=(polc_tomo_recon);
  }

  template<typename precision>
    pseudo_open_loop_tomographic_reconstructor<precision>::
    pseudo_open_loop_tomographic_reconstructor(const char * filename){
    this->target_emitter = arg_emitter = NULL;
    this->sigma_ba_inv_sigma_aa = NULL;
    this->sigma_ba_inv_sigma_aa_sigma_ab = NULL;
    this->sigma_ba_inv_sigma_aa_sigma_ac = NULL;
    this->sigma_cc = NULL;
    this->aperture_mask = NULL;
    this->read(filename);
  }

  template<typename precision>
    pseudo_open_loop_tomographic_reconstructor<precision>::
    pseudo_open_loop_tomographic_reconstructor(const Arroyo::iofits & iof){
    this->target_emitter = arg_emitter = NULL;
    this->sigma_ba_inv_sigma_aa = NULL;
    this->sigma_ba_inv_sigma_aa_sigma_ab = NULL;
    this->sigma_ba_inv_sigma_aa_sigma_ac = NULL;
    this->sigma_cc = NULL;
    this->aperture_mask = NULL;
    this->read(iof);
  }
  
  template<typename precision>
    pseudo_open_loop_tomographic_reconstructor<precision>::
    pseudo_open_loop_tomographic_reconstructor(const emitter & target,
					       const vector<emitter*> & guide_star_emitters,
					       const refractive_atmospheric_model & ref_atm_model_at_zenith,
					       double pupil_plane_pixel_scale_meters,
					       int nsteps_in_integration,
					       const circular_aperture & circ_ap,
					       bool perform_SVD){

    /////////////////////
    // Check arguments //
    /////////////////////
    try{
      dynamic_cast<const plane_wave_emitter &>(target);
    } catch(...) {
      cerr << "pseudo_open_loop_tomographic_reconstructor::pseudo_open_loop_tomographic_reconstructor error\n"
	   << "target emitter must be a plane wave emitter\n";
      throw(string("pseudo_open_loop_tomographic_reconstructor::pseudo_open_loop_tomographic_reconstructor"));
    }
 
    if(guide_star_emitters.size()==0){
      cerr << "pseudo_open_loop_tomographic_reconstructor::pseudo_open_loop_tomographic_reconstructor error\n"
	   << "no guide stars passed to the constructor\n";
      throw(string("pseudo_open_loop_tomographic_reconstructor::pseudo_open_loop_tomographic_reconstructor"));
    }
   
    for(int i=0; i<guide_star_emitters.size(); i++){
      if(dynamic_cast<plane_wave_emitter *>(guide_star_emitters[i])==NULL){
	cerr << "pseudo_open_loop_tomographic_reconstructor::pseudo_open_loop_tomographic_reconstructor error\n"
	     << "guide star emitters must be plane wave emitters\n";
	cerr << "Emitter " << i << " is not a plane wave emitter\n";
	throw(string("pseudo_open_loop_tomographic_reconstructor::pseudo_open_loop_tomographic_reconstructor"));
      }
    }
    // Here you should make sure the guide stars are not degenerate!!!


    if(pupil_plane_pixel_scale_meters<=0){
      cerr << "pseudo_open_loop_tomographic_reconstructor::pseudo_open_loop_tomographic_reconstructor error\n"
	   << "pupil plane pixel scale "
	   << pupil_plane_pixel_scale_meters
	   << " must be greater than zero\n";
      throw(string("pseudo_open_loop_tomographic_reconstructor::pseudo_open_loop_tomographic_reconstructor"));
    }
   
    if(nsteps_in_integration<=0){
      cerr << "pseudo_open_loop_tomographic_reconstructor::pseudo_open_loop_tomographic_reconstructor error\n"
	   << "nsteps in integration "
	   << nsteps_in_integration
	   << " must be greater than zero\n";
      throw(string("pseudo_open_loop_tomographic_reconstructor::pseudo_open_loop_tomographic_reconstructor"));
    }
    
    if(circ_ap.get_diameter()<=0){
      cerr << "pseudo_open_loop_tomographic_reconstructor::pseudo_open_loop_tomographic_reconstructor error\n"
	   << "circular aperture diameter "
	   << circ_ap.get_diameter()
	   << " must be greater than zero\n";
      throw(string("pseudo_open_loop_tomographic_reconstructor::pseudo_open_loop_tomographic_reconstructor"));
    }

    /////////////////////////////
    // Initialize data members //
    /////////////////////////////
    this->ref_atm_model_at_zenith = ref_atm_model_at_zenith;

    this->pupil_plane_pixel_scale_meters = pupil_plane_pixel_scale_meters;
    
    this->nsteps_in_integration = nsteps_in_integration;

    this->circ_ap = circ_ap;

    this->target_emitter = emitter::emitter_factory(&target);
    this->arg_emitter = NULL;
    for(int i=0; i<guide_star_emitters.size(); i++){
      this->guide_star_emitters.push_back(emitter::emitter_factory(guide_star_emitters[i]));
    }

    //////////////////////////////////////////////////
    // Compute the number of actuators in the pupil //
    //////////////////////////////////////////////////
    vector<long> axes(2,(long)ceil(this->circ_ap.get_diameter()/this->pupil_plane_pixel_scale_meters));
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

    this->aperture_mask = NULL;
    this->initialize_aperture_mask();
    if(pseudo_open_loop_tomographic_reconstructor::verbose_level)
      cout << this->nactuators << " out of " << axes[0]*axes[1] << " actuators lie in the pupil plane\n";

    //////////////////////////////////
    // Allocate memory for sigma_aa //
    //////////////////////////////////
    precision *sigma_aa;
    vector<long> sigma_aa_axes(2, this->guide_star_emitters.size()*nactuators);
    try{
      sigma_aa = new precision[sigma_aa_axes[0]*sigma_aa_axes[1]];
    } catch(...) {
      cerr << "pseudo_open_loop_tomographic_reconstructor::pseudo_open_loop_tomographic_reconstructor\n"
	   << " error allocating memory for sigma_aa\n";
      throw(string("pseudo_open_loop_tomographic_reconstructor::pseudo_open_loop_tomographic_reconstructor"));
    }
    
    /////////////////////////
    // Initialize sigma_aa //
    /////////////////////////
    clock_t ci = clock();

    try{
      this->initialize_auto_correlation_matrix(this->guide_star_emitters,
					       sigma_aa);
    } catch(...) {
      cerr << "pseudo_open_loop_tomographic_reconstructor::pseudo_open_loop_tomographic_reconstructor\n"
	   << "error initializing sigma_aa\n";
      throw(string("pseudo_open_loop_tomographic_reconstructor::pseudo_open_loop_tomographic_reconstructor"));
    }

    if(pseudo_open_loop_tomographic_reconstructor<precision>::verbose_level){
      cout << "pseudo_open_loop_tomographic_reconstructor::pseudo_open_loop_tomographic_reconstructor\n"
	   << "\tformed sigma_aa in " 
	   << (clock()-ci)/(double)CLOCKS_PER_SEC 
	   << endl;
      iofits iof;
      iof.create("sigma_aa.fits");
      iof.create_image(sigma_aa_axes, sigma_aa);
      iof.write_image(0, 
		      sigma_aa_axes[0]*sigma_aa_axes[1]-1, 
		      sigma_aa);
    }

    /////////////////////////////
    // Invert sigma_aa via SVD //
    // (consider using LAPACK  //
    // ssysv.f/dsysv.f         //
    // sspsf.f/dspsv.f)        //
    /////////////////////////////

    try{
      ci = clock();

      this->invert_sigma_aa_via_SVD(sigma_aa);

      if(pseudo_open_loop_tomographic_reconstructor<precision>::verbose_level){
	cout << "pseudo_open_loop_tomographic_reconstructor::pseudo_open_loop_tomographic_reconstructor\n"
	     << "\tinverted sigma_aa in " 
	     << (clock()-ci)/(double)CLOCKS_PER_SEC 
	     << endl;
      }
    } catch(...) {
      cerr << "pseudo_open_loop_tomographic_reconstructor::pseudo_open_loop_tomographic_reconstructor\n"
	   << "error inverting sigma_aa\n";
      throw(string("pseudo_open_loop_tomographic_reconstructor::pseudo_open_loop_tomographic_reconstructor"));
    }


    ////////////////////////////////////////
    // Allocate memory for other matrices //
    ////////////////////////////////////////
    precision *sigma_ba;
    this->sigma_cc = NULL;
    this->sigma_ba_inv_sigma_aa_sigma_ac = NULL;
    vector<long> sigma_ba_axes(2, this->guide_star_emitters.size()*nactuators);
    sigma_ba_axes[0] = nactuators;
    try{
      this->sigma_ba_inv_sigma_aa_sigma_ab = new precision[nactuators*nactuators];
      this->sigma_ba_inv_sigma_aa = new precision[sigma_ba_axes[0]*sigma_ba_axes[1]];
      sigma_ba = new precision[sigma_ba_axes[0]*sigma_ba_axes[1]];
    } catch(...) {
      cerr << "pseudo_open_loop_tomographic_reconstructor::pseudo_open_loop_tomographic_reconstructor\n"
	   << " error allocating memory for covariance arrays\n";
      throw(string("pseudo_open_loop_tomographic_reconstructor::pseudo_open_loop_tomographic_reconstructor"));
    }

    
    /////////////////////////
    // Initialize sigma_ba //
    /////////////////////////
    try{

      ci = clock();

      this->initialize_cross_correlation_matrix(this->target_emitter,
						this->guide_star_emitters,
						sigma_ba);

      if(pseudo_open_loop_tomographic_reconstructor<precision>::verbose_level){
	cout << "pseudo_open_loop_tomographic_reconstructor::pseudo_open_loop_tomographic_reconstructor\n"
	     << "\tformed sigma_ba in " 
	     << (clock()-ci)/(double)CLOCKS_PER_SEC 
	     << endl;
	iofits iof;
	iof.create("sigma_ba.fits");
	iof.create_image(sigma_ba_axes, 
			 sigma_ba);
	iof.write_image(0, 
			sigma_ba_axes[0]*sigma_ba_axes[1]-1, 
			sigma_ba);
      }
    } catch(...) {
      cerr << "pseudo_open_loop_tomographic_reconstructor::pseudo_open_loop_tomographic_reconstructor\n"
	   << "error forming sigma_ba\n";
      throw(string("pseudo_open_loop_tomographic_reconstructor::pseudo_open_loop_tomographic_reconstructor"));
    }
      

    ////////////////////////////////////////
    // Initialize sigma_ba_inv_sigma_aa   //
    ////////////////////////////////////////

    try{

      ci = clock();
      double tmp;
      int index;

      for(int i=0; i<sigma_ba_axes[0]; i++){
	for(int j=0; j<sigma_ba_axes[1]; j++){
	  tmp = 0;
	  index = j*sigma_ba_axes[1];
	  for(int k=0; k<sigma_ba_axes[1]; k++){
	    tmp += sigma_ba[k*sigma_ba_axes[0]+i]*sigma_aa[index+k];
	  }
	  this->sigma_ba_inv_sigma_aa[j*sigma_ba_axes[0]+i]=tmp;
	}
      }

      if(pseudo_open_loop_tomographic_reconstructor<precision>::verbose_level){
	cout << "pseudo_open_loop_tomographic_reconstructor::pseudo_open_loop_tomographic_reconstructor\n"
	     << "\tformed sigma_ba_inv_sigma_aa in " 
	     << (clock()-ci)/(double)CLOCKS_PER_SEC 
	     << endl;
	iofits iof;
	iof.create("sigma_ba_inv_sigma_aa.fits");
	iof.create_image(sigma_ba_axes, 
			 this->sigma_ba_inv_sigma_aa);
	iof.write_image(0, 
			sigma_ba_axes[0]*sigma_ba_axes[1]-1, 
			this->sigma_ba_inv_sigma_aa);
      }
    } catch(...) {
      cerr << "pseudo_open_loop_tomographic_reconstructor::pseudo_open_loop_tomographic_reconstructor\n"
	   << "error forming sigma_ba_inv_sigma_aa\n";
      throw(string("pseudo_open_loop_tomographic_reconstructor::pseudo_open_loop_tomographic_reconstructor"));
    }
      
    delete [] sigma_aa;

    ///////////////////////////////////////////////
    // Initialize sigma_ba_inv_sigma_aa_sigma_ab //
    ///////////////////////////////////////////////
    try{
      ci = clock();
      double tmp;
      for(int i=0; i<this->nactuators; i++){
	for(int j=0; j<this->nactuators; j++){
	  tmp = 0;
	  for(int k=0; k<sigma_ba_axes[1]; k++){
	    tmp += this->sigma_ba_inv_sigma_aa[k*sigma_ba_axes[0]+i]*
	      sigma_ba[k*sigma_ba_axes[0]+j];
	  }

	  this->sigma_ba_inv_sigma_aa_sigma_ab[i*this->nactuators+j]=tmp;
	}
      }

      if(pseudo_open_loop_tomographic_reconstructor<precision>::verbose_level){
	cout << "pseudo_open_loop_tomographic_reconstructor::pseudo_open_loop_tomographic_reconstructor\n"
	     << "\tformed sigma_ba_inv_sigma_aa_sigma_ab in " 
	     << (clock()-ci)/(double)CLOCKS_PER_SEC 
	     << endl;
	iofits iof;
	iof.create("sigma_ba_inv_sigma_aa_sigma_ab.fits");
	iof.create_image(vector<long>(2,this->nactuators), 
			 this->sigma_ba_inv_sigma_aa_sigma_ab);
	iof.write_image(0,
			this->nactuators*this->nactuators-1, 
			this->sigma_ba_inv_sigma_aa_sigma_ab);
      }
    } catch(...) {
      cerr << "pseudo_open_loop_tomographic_reconstructor::pseudo_open_loop_tomographic_reconstructor\n"
	   << "error forming sigma_ba_inv_sigma_aa_sigma_ab\n";
      throw(string("pseudo_open_loop_tomographic_reconstructor::pseudo_open_loop_tomographic_reconstructor"));
    }
      
    
    delete [] sigma_ba;

  }

  template<typename precision>
    pseudo_open_loop_tomographic_reconstructor<precision>::
    ~pseudo_open_loop_tomographic_reconstructor(){

    if(target_emitter!=NULL) delete target_emitter;
    for(int i=0; i<this->guide_star_emitters.size(); i++)
      delete this->guide_star_emitters[i];

    if(arg_emitter!=NULL) delete arg_emitter;
    if(sigma_ba_inv_sigma_aa!=NULL) delete [] sigma_ba_inv_sigma_aa;
    if(sigma_ba_inv_sigma_aa_sigma_ab!=NULL){
      delete [] sigma_ba_inv_sigma_aa_sigma_ab;
      // do this in case we've set the pointer
      // sigma_ba_inv_sigma_aa_sigma_ac = sigma_ba_inv_sigma_aa_sigma_ac
      sigma_ba_inv_sigma_aa_sigma_ab = NULL;
    }
    if(sigma_ba_inv_sigma_aa_sigma_ac!=NULL) delete [] sigma_ba_inv_sigma_aa_sigma_ac;
    if(sigma_cc!=NULL) delete [] sigma_cc;
    if(aperture_mask!=NULL) delete [] aperture_mask;
  }

  template<typename precision>
    pseudo_open_loop_tomographic_reconstructor<precision> & 
    pseudo_open_loop_tomographic_reconstructor<precision>::
    operator=(const pseudo_open_loop_tomographic_reconstructor<precision> & polc_tomo_recon){

    if(this==&polc_tomo_recon)
      return(*this);

    if(this->target_emitter!=NULL)
      delete this->target_emitter;
    target_emitter = emitter::emitter_factory(polc_tomo_recon.target_emitter);

    if(this->arg_emitter!=NULL)
      delete this->arg_emitter;
    arg_emitter = NULL;

    for(int i=0; i<this->guide_star_emitters.size();i++){
      if(this->guide_star_emitters[i]!=NULL)
	delete this->guide_star_emitters[i];
    }

    this->guide_star_emitters.resize(polc_tomo_recon.guide_star_emitters.size());
    for(int i=0; i<this->guide_star_emitters.size();i++){
      this->guide_star_emitters[i] = 
	emitter::emitter_factory(polc_tomo_recon.guide_star_emitters[i]);
    }

    this->nactuators = polc_tomo_recon.nactuators;
    this->nsteps_in_integration = polc_tomo_recon.nsteps_in_integration;
    this->pupil_plane_pixel_scale_meters = polc_tomo_recon.pupil_plane_pixel_scale_meters;
    this->circ_ap = polc_tomo_recon.circ_ap;
    this->ref_atm_model_at_zenith = polc_tomo_recon.ref_atm_model_at_zenith;

    if(this->sigma_ba_inv_sigma_aa!=NULL)
      delete [] this->sigma_ba_inv_sigma_aa;

    int nelem = this->nactuators*this->nactuators*this->guide_star_emitters.size();
    try{
      this->sigma_ba_inv_sigma_aa = 
      new precision[nelem];
    } catch(...) {
      cerr << "pseudo_open_loop_tomographic_reconstructor::operator=\n"
	   << "error allocating memory\n";
      throw(string("pseudo_open_loop_tomographic_reconstructor::pseudo_open_loop_tomographic_reconstructor"));
    }

    for(int i=0; i<nelem; i++)
      this->sigma_ba_inv_sigma_aa[i] = 
	polc_tomo_recon.sigma_ba_inv_sigma_aa[i];


    if(this->sigma_ba_inv_sigma_aa_sigma_ab!=NULL)
      delete [] this->sigma_ba_inv_sigma_aa_sigma_ab;

    nelem = this->nactuators*this->nactuators;
    try{
      this->sigma_ba_inv_sigma_aa_sigma_ab = 
      new precision[nelem];
    } catch(...) {
      cerr << "pseudo_open_loop_tomographic_reconstructor::operator=\n"
	   << "error allocating memory\n";
      throw(string("pseudo_open_loop_tomographic_reconstructor::pseudo_open_loop_tomographic_reconstructor"));
    }

    for(int i=0; i<nelem; i++)
      this->sigma_ba_inv_sigma_aa_sigma_ab[i] = 
	polc_tomo_recon.sigma_ba_inv_sigma_aa_sigma_ab[i];


    if(this->sigma_ba_inv_sigma_aa_sigma_ac!=NULL){
      delete [] this->sigma_ba_inv_sigma_aa_sigma_ac;
      this->sigma_ba_inv_sigma_aa_sigma_ac = NULL;
    }

    if(this->sigma_cc!=NULL){
      delete [] this->sigma_cc;
      this->sigma_cc = NULL;
    }

    if(this->aperture_mask!=NULL)
      delete [] this->aperture_mask;
    this->aperture_mask = NULL;

    return(*this);
  }

  template<typename precision>
    void pseudo_open_loop_tomographic_reconstructor<precision>::
    read(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "pseudo_open_loop_tomographic_reconstructor::read - "
	   << "error opening file " << filename << endl;
      throw(string("pseudo_open_loop_tomographic_reconstructor::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "pseudo_open_loop_tomographic_reconstructor::read - "
	   << "error reading "
	   << this->unique_name() << " from file " 
	   << filename << endl;
      throw(string("pseudo_open_loop_tomographic_reconstructor::read"));
    }
  }
  
  template<typename precision>
    void pseudo_open_loop_tomographic_reconstructor<precision>::
    read(const Arroyo::iofits & iof){

    if(!iof.key_exists("TYPE")){
      cerr << "pseudo_open_loop_tomographic_reconstructor::read error - "
	   << "unrecognized type of file\n";
      throw(string("pseudo_open_loop_tomographic_reconstructor::read"));
    }

    string type, comment;
    iof.read_key("TYPE", type, comment);
    if(type!=this->unique_name()){
      cerr << "pseudo_open_loop_tomographic_reconstructor::read error - file of type " 
	   << type << " rather than of type "
	   << this->unique_name() << endl;
      throw(string("pseudo_open_loop_tomographic_reconstructor::read"));
    }

    int nguidestars;
    iof.read_key("NGUIDEST", nguidestars, comment);

    iof.read_key("NACTUATR", this->nactuators, comment);

    iof.read_key("NSTEPINT", this->nsteps_in_integration, comment);

    iof.read_key("PUPIXSCL", this->pupil_plane_pixel_scale_meters, comment);

    iof.movrel_hdu(1);
    int nelem = this->nactuators*this->nactuators*this->guide_star_emitters.size();
    try{
      if(this->sigma_ba_inv_sigma_aa!=NULL)
	delete [] this->sigma_ba_inv_sigma_aa;

      this->sigma_ba_inv_sigma_aa = 
	new precision[nelem];
    } catch(...) {
      cerr << "pseudo_open_loop_tomographic_reconstructor::read - error allocating memory\n";
      throw(string("pseudo_open_loop_tomographic_reconstructor::read"));
    }
    iof.read_image(0,nelem-1,this->sigma_ba_inv_sigma_aa);
    iof.movrel_hdu(1);


    nelem = this->nactuators*this->nactuators;
    try{
      if(this->sigma_ba_inv_sigma_aa_sigma_ab!=NULL)
	delete [] this->sigma_ba_inv_sigma_aa_sigma_ab;
      
      this->sigma_ba_inv_sigma_aa_sigma_ab =
	new precision[nelem];
    } catch(...) {
      cerr << "pseudo_open_loop_tomographic_reconstructor::read - error allocating memory\n";
      throw(string("pseudo_open_loop_tomographic_reconstructor::read"));
    }
    iof.read_image(0,nelem-1,this->sigma_ba_inv_sigma_aa_sigma_ab);
    iof.movrel_hdu(1);


    if(this->sigma_ba_inv_sigma_aa_sigma_ac!=NULL){
      delete [] this->sigma_ba_inv_sigma_aa_sigma_ac;
      this->sigma_ba_inv_sigma_aa_sigma_ac = NULL;
    }

    if(this->sigma_cc!=NULL){
      delete [] this->sigma_cc;
      this->sigma_cc = NULL;
    }

    if(this->target_emitter!=NULL)
      delete this->target_emitter;
    this->target_emitter = emitter::emitter_factory(iof);

    for(int i=0; i<this->guide_star_emitters.size(); i++){
      if(this->guide_star_emitters[i]!=NULL)
	delete this->guide_star_emitters[i];
    }


    this->guide_star_emitters.resize(nguidestars);
    for(int i=0; i<this->guide_star_emitters.size(); i++){
      this->guide_star_emitters[i] = emitter::emitter_factory(iof);
    }

    this->circ_ap.read(iof);

    this->ref_atm_model_at_zenith.read(iof);
    
    if(this->aperture_mask!=NULL)
      delete [] this->aperture_mask;
    this->aperture_mask = NULL;

  }
  
  template<typename precision>
    void pseudo_open_loop_tomographic_reconstructor<precision>::
    write(const char * filename) const {

    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "pseudo_open_loop_tomographic_reconstructor::write - "
	   << "error opening file " << filename << endl;
      throw(string("pseudo_open_loop_tomographic_reconstructor::write"));
    }

    try{this->write(iof);}
    catch(...){
      cerr << "pseudo_open_loop_tomographic_reconstructor::write - "
	   << "error writing "
	   << this->unique_name() << " to file " 
	   << filename << endl;
      throw(string("pseudo_open_loop_tomographic_reconstructor::write"));
    }
  }
  
  template<typename precision>
    void pseudo_open_loop_tomographic_reconstructor<precision>::
    write(Arroyo::iofits & iof) const {

    fits_header_data<precision> fhd;
    fhd.write(iof);
    string type = this->unique_name();
    string comment = "object type";
    iof.write_key("TYPE", type, comment);
    
    comment = "number of guidestars";
    iof.write_key("NGUIDEST", (long)this->guide_star_emitters.size(), comment);

    comment = "number of actuators to control";
    iof.write_key("NACTUATR", this->nactuators, comment);

    comment = "number of steps in numerical integration";
    iof.write_key("NSTEPINT", this->nsteps_in_integration, comment);

    comment = "pupil plane pixel scale (meters)";
    iof.write_key("PUPIXSCL", this->pupil_plane_pixel_scale_meters, comment);

    vector<long> first_axes(2, this->nactuators);
    first_axes[1] *= this->guide_star_emitters.size();
    fhd = fits_header_data<precision>(first_axes);
    fhd.write(iof);
    iof.write_image(0,
		    first_axes[0]*first_axes[1]-1,
		    this->sigma_ba_inv_sigma_aa);

    vector<long> second_axes(2, this->nactuators);
    fhd = fits_header_data<precision>(second_axes);
    fhd.write(iof);
    iof.write_image(0,
		    second_axes[0]*second_axes[1]-1,
		    this->sigma_ba_inv_sigma_aa_sigma_ab);

    this->target_emitter->write(iof);

    for(int i=0; i<this->guide_star_emitters.size(); i++){
      this->guide_star_emitters[i]->write(iof);
    }

    this->circ_ap.write(iof);

    this->ref_atm_model_at_zenith.write(iof);
  }
  
  template<typename precision>
    void pseudo_open_loop_tomographic_reconstructor<precision>::
    print(std::ostream & os, const char * prefix) const {

    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield);
    os << prefix << "TYPE       = " << setw(vlspc) << this->unique_name()
       << "/" << "object type" << endl;
    os << prefix << "NGUIDEST   = " << setw(vlspc) << this->guide_star_emitters.size()
       << "/" << "number of guidestars" << endl;
    os << prefix << "NACTUATR   = " << setw(vlspc) << this->nactuators
       << "/" << "number of actuators" << endl;
    os << prefix << "NSTEPINT   = " << setw(vlspc) << this->nsteps_in_integration
       << "/" << "number of steps in numerical integration" << endl;
    os << prefix << "PUPIXSCL   = " << setw(vlspc) << this->pupil_plane_pixel_scale_meters
       << "/" << "pupil plane pixel scale (meters)" << endl;
    os << prefix << "APDIAM     = " << setw(vlspc) << this->circ_ap.get_diameter()
       << "/" << "aperture diameter (meters)" << endl;
  }

  template<typename precision>
    double pseudo_open_loop_tomographic_reconstructor<precision>::
    aperture_averaged_phase_variance(emitter & emtr,
				     double wavelength_meters) const {
    return(private_aperture_averaged_phase_variance(emtr,
						    wavelength_meters,
						    false));
  }
    
  template<typename precision>
    double pseudo_open_loop_tomographic_reconstructor<precision>::
    aperture_averaged_differential_phase_variance(emitter & emtr,
						  double wavelength_meters) const {
    return(private_aperture_averaged_phase_variance(emtr,
						    wavelength_meters,
						    true));
  }

  template<typename precision>
    double pseudo_open_loop_tomographic_reconstructor<precision>::
    private_aperture_averaged_phase_variance(emitter & emtr,
					     double wavelength_meters,
					     bool differential) const {

    if(wavelength_meters<=0){
      cerr << "pseudo_open_loop_tomographic_reconstructor::private_aperture_averaged_phase_variance error -\n"
	   << "wavelength "
	   << wavelength_meters
	   << " out of range\n";
      throw(string("pseudo_open_loop_tomographic_reconstructor::private_aperture_averaged_phase_variance"));
    }

    this->initialize_emitter(emtr);

    // sum up the diagonal elements of sigma_cc and subtract 
    // the sum of the diagonal elements from sigma_ba_inv_sigma_aa_sigma_ac
    //
    // then divide by the number of points
    int count = 0;
    double variance = 0;
    
    if(differential){
      for(int i=0; i<this->nactuators; i++){
	variance += this->sigma_cc[i*this->nactuators+i];
	variance += this->sigma_ba_inv_sigma_aa_sigma_ab[i*this->nactuators+i];
	variance -= 2*this->sigma_ba_inv_sigma_aa_sigma_ac[i*this->nactuators+i];
	count++;
      }
    } else {
      for(int i=0; i<this->nactuators; i++){
	variance += this->sigma_cc[i*this->nactuators+i];
	count++;
      }
    }      

    cout << "raw variance " << variance << endl;

    /*
    if(variance>-sqrt(three_frame::precision) &&
       variance<0)
      return(0);
    else 
    */
    
    if(variance < 0){
      cerr << "pseudo_open_loop_tomographic_reconstructor::private_aperture_averaged_phase_variance error -\n"
	   << " variance "
	   << variance 
	   << " out of range "
	   << -three_frame::precision
	   << endl;
      return(4*M_PI*M_PI*variance/wavelength_meters/wavelength_meters/(double)count);
      //throw(string("pseudo_open_loop_tomographic_reconstructor::private_aperture_averaged_phase_variance"));
    } else 
      return(4*M_PI*M_PI*variance/wavelength_meters/wavelength_meters/(double)count);
  }
  
  template<typename precision>
    pixel_array<precision> pseudo_open_loop_tomographic_reconstructor<precision>::
    phase_variance(emitter & emtr,
		   double wavelength_meters) const {

    return(this->private_phase_variance(emtr,
					wavelength_meters,
					false));
  }
  
  template<typename precision>
    pixel_array<precision> pseudo_open_loop_tomographic_reconstructor<precision>::
    differential_phase_variance(emitter & emtr,
				double wavelength_meters) const {

    return(this->private_phase_variance(emtr,
					wavelength_meters,
					true));
  }

  template<typename precision>
    pixel_array<precision> pseudo_open_loop_tomographic_reconstructor<precision>::
    private_phase_variance(emitter & emtr,
			   double wavelength_meters,
			   bool differential) const {
    if(wavelength_meters<=0){
      cerr << "pseudo_open_loop_tomographic_reconstructor::private_phase_variance error -\n"
	   << "wavelength "
	   << wavelength_meters
	   << " out of range\n";
      throw(string("pseudo_open_loop_tomographic_reconstructor::private_phase_variance"));
    }

    this->initialize_emitter(emtr);

    vector<long> axes(2,(long)ceil(this->circ_ap.get_diameter()/this->pupil_plane_pixel_scale_meters));
    pixel_array<precision> pixarr(axes);

    int count = 0;
    int index;
    double fac = 4*M_PI*M_PI/wavelength_meters/wavelength_meters;

    if(differential){
      for(int i=0; i<axes[1]; i++){
	for(int j=0; j<axes[0]; j++){
	  index = i*axes[0]+j;
	  if(this->aperture_mask[index]){
	    pixarr.set_data(index, 
			    fac*(sigma_cc[count*this->nactuators+count] +
				 this->sigma_ba_inv_sigma_aa_sigma_ab[count*this->nactuators+count] -
				 2*sigma_ba_inv_sigma_aa_sigma_ac[count*this->nactuators+count]));
	    count++;
	  }
	}
      }
    } else {
      for(int i=0; i<axes[1]; i++){
	for(int j=0; j<axes[0]; j++){
	  index = i*axes[0]+j;
	  if(this->aperture_mask[index]){
	    pixarr.set_data(index, 
			    fac*sigma_cc[count*this->nactuators+count]);
	    count++;
	  }
	}
      }
    }

    return(pixarr);
  }
  
  template<typename precision>
    pixel_array<precision> pseudo_open_loop_tomographic_reconstructor<precision>::
    phase_covariance(emitter & emtr,
		     int xindex,
		     int yindex,
		     double wavelength_meters) const {

    return(private_phase_covariance(emtr,
				    xindex,
				    yindex,
				    wavelength_meters,
				    false));
  }

  template<typename precision>
    pixel_array<precision> pseudo_open_loop_tomographic_reconstructor<precision>::
    differential_phase_covariance(emitter & emtr,
				  int xindex,
				  int yindex,
				  double wavelength_meters) const {
    
    return(private_phase_covariance(emtr,
				    xindex,
				    yindex,
				    wavelength_meters,
				    true));
  }

  template<typename precision>
    pixel_array<precision> pseudo_open_loop_tomographic_reconstructor<precision>::
    private_phase_covariance(emitter & emtr,
			     int xindex,
			     int yindex,
			     double wavelength_meters,
			     bool differential) const {
    
    if(wavelength_meters<=0){
      cerr << "pseudo_open_loop_tomographic_reconstructor::private_phase_covariance error -\n"
	   << "wavelength "
	   << wavelength_meters
	   << " out of range\n";
      throw(string("pseudo_open_loop_tomographic_reconstructor::private_phase_covariance"));
    }

    this->initialize_emitter(emtr);
    vector<long> axes(2,(long)ceil(this->circ_ap.get_diameter()/this->pupil_plane_pixel_scale_meters));
    pixel_array<precision> pixarr(axes);

    int ref_index = xindex*axes[0]+yindex;
    if(ref_index>pixarr.total_space()){
      cerr << "pseudo_open_loop_tomographic_reconstructor::private_phase_covariance error -\n"
	   << "index "
	   << xindex 
	   << ", "
	   << yindex
	   << " lies outside of the array\n";
      throw(string("pseudo_open_loop_tomographic_reconstructor::private_phase_covariance"));
    }
      
    if(this->aperture_mask[ref_index]==0){
      cerr << "pseudo_open_loop_tomographic_reconstructor::private_phase_covariance error -\n"
	   << "index "
	   << xindex 
	   << ", "
	   << yindex
	   << " lies outside of aperture\n";
      throw(string("pseudo_open_loop_tomographic_reconstructor::private_phase_covariance"));
    }

    // Convert xindex,yindex into an array index
    int array_index = 0;
    for(int i=0; i<xindex*axes[0]+yindex; i++)
      if(this->aperture_mask[i])
	array_index++;
    
    int index, count=0;
    double fac = 4*M_PI*M_PI/wavelength_meters/wavelength_meters;

    if(differential){
      for(int i=0; i<axes[1]; i++){
	for(int j=0; j<axes[0]; j++){
	  index = i*axes[0]+j;
	  if(this->aperture_mask[index]){
	    pixarr.set_data(index, 
			    fac*(sigma_cc[array_index*this->nactuators+count] + 
				 sigma_ba_inv_sigma_aa_sigma_ab[array_index*this->nactuators+count] -
				 sigma_ba_inv_sigma_aa_sigma_ac[array_index*this->nactuators+count] -
				 sigma_ba_inv_sigma_aa_sigma_ac[count*this->nactuators+array_index]));
	    count++;
	  }
	}
      }
    } else {
      for(int i=0; i<axes[1]; i++){
	for(int j=0; j<axes[0]; j++){
	  index = i*axes[0]+j;
	  if(this->aperture_mask[index]){
	    pixarr.set_data(index, 
			    fac*sigma_cc[array_index*this->nactuators+count]);
	    count++;
	  }
	}
      }
    }
    return(pixarr);
  }
  
  template<typename precision>
    pixel_array<precision> pseudo_open_loop_tomographic_reconstructor<precision>::
    phase_structure_function(const emitter & emtr,
			     double wavelength_meters,
			     int xindex,
			     int yindex) const {

    if(wavelength_meters<=0){
      cerr << "pseudo_open_loop_tomographic_reconstructor::phase_structure_function error -\n"
	   << "wavelength "
	   << wavelength_meters
	   << " out of range\n";
      throw(string("pseudo_open_loop_tomographic_reconstructor::phase_structure_function"));
    }

    this->initialize_emitter(emtr);


    pixel_array<precision> pixarr;

    return(pixarr); 

 }
  
  template<typename precision>
    basic_otf<precision> pseudo_open_loop_tomographic_reconstructor<precision>::
    optical_transfer_function(const emitter & emtr,
			      double wavelength_meters) const {

    if(wavelength_meters<=0){
      cerr << "pseudo_open_loop_tomographic_reconstructor::optical_transfer_function error -\n"
	   << "wavelength "
	   << wavelength_meters
	   << " out of range\n";
      throw(string("pseudo_open_loop_tomographic_reconstructor::optical_transfer_function"));
    }

    this->initialize_emitter(emtr);

    basic_otf<precision> otf;

    return(otf);

  }
  
  
  template<typename precision>
    basic_observation<precision> pseudo_open_loop_tomographic_reconstructor<precision>::
    point_spread_function(const emitter & emtr,
			  double wavelength_meters,
			  double field_size_arcsecs,
			  double oversampling_factor) const {

    if(wavelength_meters<=0){
      cerr << "pseudo_open_loop_tomographic_reconstructor::point_spread_function error -\n"
	   << "wavelength "
	   << wavelength_meters
	   << " out of range\n";
      throw(string("pseudo_open_loop_tomographic_reconstructor::point_spread_function"));
    }

    this->initialize_emitter(emtr);


    basic_observation<precision> basic_obs;

    return(basic_obs);
  }

  template<typename precision>
    void pseudo_open_loop_tomographic_reconstructor<precision>::
    reconstruct(const vector<Arroyo::pixel_array<double> > & measured_phases,
		Arroyo::pixel_array<double> & commands) const {

  }
  


}

#endif
