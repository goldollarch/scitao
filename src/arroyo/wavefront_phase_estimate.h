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

#ifndef WAVEFRONT_PHASE_ESTIMATE_H
#define WAVEFRONT_PHASE_ESTIMATE_H

#include <string>
#include <iostream>
#include <vector>
#include "AO_sim_base.h"
//#include "phase_estimate.h"
#include "covariance.h"
#include "aperture.h"
#include "observation.h"

namespace Arroyo {

  class iofits;
  using std::string;

  ///    
  /// A class to hold a wavefront phase
  /// estimate formed from the linear 
  /// combination of an arbitrary number of 
  /// tip tilt guide stars and an arbitrary
  /// number of high order guide stars
  ///

  template<typename precision, typename aperture_type>
    class wavefront_phase_estimate :
    //public phase_estimate_base {
    public AO_sim_base {

    private:
    
    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("wavefront phase estimate"));};

    //////////////////////////////////////////////////////////////////
    ///  Returns a pixel array containing the long exposure point
    ///  spread function
    basic_observation<precision> private_point_spread_function(double field_size_arcsecs,
							       double oversampling_factor,
							       const basic_otf<precision> & long_exposure_OTF) const;


    protected:

    refractive_atmospheric_model ref_atm_model;
    aperture_type ap;

    mutable emitter * stored_emtr;

    mutable phase_covariance<precision, aperture_type> phase_covariance_aa;
    mutable tilt_covariance<precision, aperture_type> tilt_covariance_aa;

    mutable vector<phase_covariance<precision, aperture_type> > ho_phase_covariance_ai;
    mutable vector<tilt_covariance<precision, aperture_type> > tt_tilt_covariance_ai;
    mutable vector<tilt_covariance<precision, aperture_type> > ho_tilt_covariance_ai;

    vector<phase_covariance<precision, aperture_type> > ho_phase_covariance_ij;
    vector<tilt_covariance<precision, aperture_type> > tt_tilt_covariance_ij;
    vector<tilt_covariance<precision, aperture_type> > ho_tilt_covariance_ij;

    vector<double> tip_tilt_weights;
    vector<double> high_order_weights;

    // Stored phase variance
    mutable pixel_array<precision> stored_phase_variance_over_k_squared;

    // Stored OTF
    mutable pixel_array<precision> stored_OTF;
    mutable double stored_OTF_wavelength;

    bool emitter_check(const emitter & emtr) const;

    void initialize_emitter(const emitter & emtr) const;

    double aperture_averaged_tilt_variance_calculation(double wavelength_meters,
						       int nsteps_in_integration,
						       bool tip_tilt_stars) const;

    double tilt_variance_calculation(double wavelength_meters,
				     int nsteps_in_integration,
				     const three_point & tp,
				     bool tip_tilt_stars) const;

    pixel_array<precision> tilt_variance_calculation(double pixel_scale_meters,
						     double wavelength_meters,
						     int nsteps_in_integration,
						     bool tip_tilt_stars) const;


    double get_aperture_outer_diameter() const;

    public:

    ///////////////////////////////////////////
    ///  Null constructor
    wavefront_phase_estimate(){
      stored_emtr = NULL;
    };

    ///////////////////////////////////////////
    ///  
    /// Construct a null instance of a wavefront phase estimate.
    /// At present, the aperture ap must be a circular_aperture, or
    /// this constructor throws an error
    wavefront_phase_estimate(const refractive_atmospheric_model & ref_atm_model,
			     const aperture_type & ap_type);

    ///////////////////////////////////////////
    ///  
    /// Construct an instance of a wavefront phase estimate
    /// that uses linear combinations of guide stars to generate
    /// tip tilt and high order estimates.
    /// At present, the aperture ap must be a circular_aperture, or
    /// this constructor throws an error
    wavefront_phase_estimate(const vector<emitter *> & tip_tilt_guide_stars,
			     const vector<double> & tip_tilt_weights,
			     const vector<emitter *> & high_order_guide_stars,
			     const vector<double> & high_order_weights,
			     const refractive_atmospheric_model & ref_atm_model,
			     const aperture_type & ap_type);

    ///////////////////////////////////////////
    ///  Destructor
    ~wavefront_phase_estimate(){
      delete stored_emtr;
    };

    ///////////////////////////////////////////
    ///  Copy constructor
    wavefront_phase_estimate(const wavefront_phase_estimate & wpe){
      this->operator=(wpe);
    };

    ///////////////////////////////////////////
    ///  Construct from file
    wavefront_phase_estimate(const char * filename){
      this->read(filename);
    };

    ///////////////////////////////////////////
    ///  Construct from an iofits object
    wavefront_phase_estimate(const iofits & iof){
      this->read(iof);
    };

    ///////////////////////////////////////////
    ///  Operator =
    wavefront_phase_estimate & operator=(const wavefront_phase_estimate & wpe);

    ///////////////////////////////////////////
    ///  read from a file
    void read(const char * filename);

    ///////////////////////////////////////////
    ///  read from an iofits object
    void read(const Arroyo::iofits & iof);

    ///////////////////////////////////////////
    ///  write to a file
    void write(const char * filename) const;

    ///////////////////////////////////////////
    ///  write to an iofits object
    void write(Arroyo::iofits & iof) const;

    ///////////////////////////////////////////
    ///  Function to print the coefficients
    void print(std::ostream & os, const char * prefix = "") const;

    ///////////////////////////////////////////
    ///  Function to print the coefficients
    refractive_atmospheric_model get_refractive_atmospheric_model() const {
      return(this->ref_atm_model);
    }

    //////////////////////////////////////////////////////////////////
    ///  Return the aperture averaged differential phase variance
    ///  between the emitter emtr and the wavefront phase estimate
    double aperture_averaged_differential_phase_variance(const emitter & emtr,
							 double wavelength_meters,
							 int nsteps_in_integration) const;

    //////////////////////////////////////////////////////////////////
    ///  Return the aperture averaged differential tilt phase variance
    ///  between the emitter emtr and the wavefront phase estimate
    double aperture_averaged_differential_tilt_phase_variance(const emitter & emtr,
							      double wavelength_meters,
							      int nsteps_in_integration) const;

    //////////////////////////////////////////////////////////////////
    ///  Return the differential phase variance between the emitter
    ///  emtr and the wavefront phase estimate
    double differential_phase_variance(const emitter & emtr,
				       const three_point & tp,
				       double wavelength_meters,
				       int nsteps_in_integration) const;

    //////////////////////////////////////////////////////////////////
    ///  Return a pixel array containing the differential phase
    ///  variance between the emitter emtr and the wavefront phase
    ///  estimate at point tp in the pupil plane.  The differential
    ///  phase is computed on a grid that covers the aperture with
    ///  sampling set by the argument pixel_scale_meters
    pixel_array<precision> differential_phase_variance(const emitter & emtr,
						       double pixel_scale_meters,
						       double wavelength_meters,
						       int nsteps_in_integration) const;
    
    //////////////////////////////////////////////////////////////////
    ///  Return the differential tilt phase variance between the
    ///  emitter emtr and the wavefront phase estimate at point tp in
    ///  the pupil plane
    double differential_tilt_phase_variance(const emitter & emtr,
					    const three_point & tp,
					    double wavelength_meters,
					    int nsteps_in_integration) const;

    //////////////////////////////////////////////////////////////////
    ///  Return a pixel array containing the differential tilt phase
    ///  variance between the emitter emtr and the wavefront phase
    ///  estimate at point tp in the pupil plane.  The differential
    ///  phase is computed on a grid that covers the aperture with
    ///  sampling set by the argument pixel_scale_meters
    pixel_array<precision> differential_tilt_phase_variance(const emitter & emtr,
							    double pixel_scale_meters,
							    double wavelength_meters,
							    int nsteps_in_integration) const;


    //////////////////////////////////////////////////////////////////
    ///  Return the value of the structure function computed between
    ///  points tp1 and tp2.  These points must lie within the
    ///  aperture, or this function throws an error.
    double phase_structure_function(const emitter & emtr,
				    double wavelength_meters,
				    int nsteps_in_integration,
				    const three_point & tp1,
				    const three_point & tp2) const;


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
						    double pixel_scale_meters,
						    double wavelength_meters,
						    int nsteps_in_integration,
						    const three_point & tp) const;
    
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
						    double pixel_scale_meters,
						    double wavelength_meters,
						    int nsteps_in_integration,
						    int xindex,
						    int yindex) const;
    
    //////////////////////////////////////////////////////////////////
    ///  Returns a pixel array containing the long exposure optical
    ///  transfer function
    basic_otf<precision> NGS_optical_transfer_function(const emitter & emtr,
						       double wavelength_meters,
						       int nsteps_in_integration,
						       double pupil_plane_pixel_scale_meters,
						       double cutoff = 30) const;

    //////////////////////////////////////////////////////////////////
    ///  Returns a pixel array containing the long exposure optical
    ///  transfer function
    basic_otf<precision> optical_transfer_function(const emitter & emtr,
						   double wavelength_meters,
						   int nsteps_in_integration,
						   double pupil_plane_pixel_scale_meters,
						   double cutoff = 30) const;
    

    //////////////////////////////////////////////////////////////////
    ///  Returns a pixel array containing the long exposure point
    ///  spread function
    basic_observation<precision> point_spread_function(const emitter & emtr,
						       double wavelength_meters,
						       int nsteps_in_integration,
						       double pupil_plane_pixel_scale_meters,
						       double field_size_arcsecs,
						       double oversampling_factor,
						       double cutoff=30) const;


    //////////////////////////////////////////////////////////////////
    ///  Returns a pixel array containing the long exposure point
    ///  spread function
    basic_observation<precision> NGS_point_spread_function(const emitter & emtr,
							   double wavelength_meters,
							   int nsteps_in_integration,
							   double pupil_plane_pixel_scale_meters,
							   double field_size_arcsecs,
							   double oversampling_factor,
							   double cutoff=30) const;
  };


  template<typename precision, typename aperture_type>
    double wavefront_phase_estimate<precision, aperture_type>::get_aperture_outer_diameter() const {
    const circular_aperture * circ_ap = dynamic_cast<const circular_aperture *>(&(this->ap));
    const annular_aperture * ann_ap = dynamic_cast<const annular_aperture *>(&(this->ap));
    if(circ_ap!=NULL)
      return(circ_ap->get_diameter());
    else if(ann_ap!=NULL)
      return(ann_ap->get_outer_diameter());
    else {
      cerr << "wavefront_phase_estimate::get_aperture_outer_diameter error\n";
      throw(string("wavefront_phase_estimate::get_aperture_outer_diameter"));
    }
  }


  template<typename precision, typename aperture_type>
    bool wavefront_phase_estimate<precision, aperture_type>::
    emitter_check(const emitter & emtr) const {

    if(stored_emtr==NULL) return(false);
    
    plane_wave_emitter * stored_pwe=dynamic_cast<plane_wave_emitter *>(this->stored_emtr);
    const plane_wave_emitter * arg_pwe=dynamic_cast<const plane_wave_emitter *>(&emtr);
    spherical_wave_emitter * stored_swe=dynamic_cast<spherical_wave_emitter *>(this->stored_emtr);
    const spherical_wave_emitter * arg_swe=dynamic_cast<const spherical_wave_emitter *>(&emtr);

    if(stored_pwe!=NULL && arg_pwe!=NULL && *stored_pwe==*arg_pwe) return true;
    if(stored_swe!=NULL && arg_swe!=NULL && *stored_swe==*arg_swe) return true;
    return false;
  }

  template<typename precision, typename aperture_type>
    void wavefront_phase_estimate<precision, aperture_type>::
    initialize_emitter(const emitter & emtr) const {
    
    try{

      if(this->emitter_check(emtr)) return;

      delete this->stored_emtr;
      this->stored_emtr = emitter::emitter_factory(&emtr);

      stored_OTF = pixel_array<double>();
      stored_OTF_wavelength = -1;

    
      phase_covariance_aa = 
	phase_covariance<precision, aperture_type>(*(this->stored_emtr),
						   *(this->stored_emtr),
						   this->ref_atm_model,
						   this->ap);
      
      tilt_covariance_aa = 
	tilt_covariance<precision, aperture_type>(*(this->stored_emtr),
						  *(this->stored_emtr),
						  this->ref_atm_model,
						  this->ap);
      
      int nhogs = high_order_weights.size();
      ho_phase_covariance_ai.resize(nhogs);
      ho_tilt_covariance_ai.resize(nhogs);

      for(int i=0; i<nhogs; i++){
	ho_phase_covariance_ai[i] = 
	  phase_covariance<precision, aperture_type>(*(this->stored_emtr), 
						     ho_phase_covariance_ij[i].get_second_emitter(),
						     this->ref_atm_model,
						     this->ap);
	ho_tilt_covariance_ai[i] = 
	  tilt_covariance<precision, aperture_type>(*(this->stored_emtr), 
						    ho_tilt_covariance_ij[i].get_second_emitter(),
						    this->ref_atm_model,
						    this->ap);
      }
      
      int nttgs = tip_tilt_weights.size();
      tt_tilt_covariance_ai.resize(nttgs);
      for(int i=0; i<nttgs; i++){
	tt_tilt_covariance_ai[i] = 
	  tilt_covariance<precision, aperture_type>(*(this->stored_emtr), 
						    tt_tilt_covariance_ij[i].get_second_emitter(),
						    this->ref_atm_model,
						    this->ap);
      }
    } catch(...) {
      cerr << "wavefront_phase_estimate::intitialize_emitter error\n";
      throw(string("wavefront_phase_estimate::initialize_emitter"));
    }
  }

  template<typename precision, typename aperture_type>
    double wavefront_phase_estimate<precision, aperture_type>::
    aperture_averaged_tilt_variance_calculation(double wavelength_meters,
						int nsteps_in_integration,
						bool tip_tilt_stars) const {


    try{
      if(this->stored_emtr==NULL){
	cerr << "wavefront_phase_estimate::aperture_averaged_tilt_variance_calculation error - "
	     << "emitter is not initialized\n";
	throw(string("wavefront_phase_estimate::aperture_averaged_tilt_variance_calculation"));
      }
      
      const vector<tilt_covariance<precision, aperture_type> > & tilt_covariance_ai = 
	tip_tilt_stars ? tt_tilt_covariance_ai : ho_tilt_covariance_ai;
      const vector<tilt_covariance<precision, aperture_type> > & tilt_covariance_ij = 
	tip_tilt_stars ? tt_tilt_covariance_ij : ho_tilt_covariance_ij;
      const vector<double> & weights = 
	tip_tilt_stars ? tip_tilt_weights : high_order_weights;
      
      int ngs = weights.size();
      
      double val = 
	tilt_covariance_aa.aperture_averaged_variance(wavelength_meters,
						      nsteps_in_integration);
      

      for(int i=0; i<ngs; i++)
	val -= 2*weights[i]*
	  tilt_covariance_ai[i].aperture_averaged_variance(wavelength_meters,
							   nsteps_in_integration);

      int count = 0;
      double fac;
      for(int i=0; i<ngs; i++){
	for(int j=i; j<ngs; j++){
	  fac = i==j ? 1. : 2.;
	  val += fac*weights[i]*weights[j]*
	    tilt_covariance_ij[count].aperture_averaged_variance(wavelength_meters,
								 nsteps_in_integration);
	  count++;
	}
      }

      return(val);
    } catch(...){
      cerr << "wavefront_phase_estimate::aperture_averaged_tilt_variance_calculation error\n";
      throw(string("wavefront_phase_estimate::aperture_averaged_tilt_variance_calculation"));
    }
  }

  template<typename precision, typename aperture_type>
    double wavefront_phase_estimate<precision, aperture_type>::
    tilt_variance_calculation(double wavelength_meters,
			      int nsteps_in_integration,
			      const three_point & tp,
			      bool tip_tilt_stars) const {

    try{

      double ap_diameter = this->get_aperture_outer_diameter();
      if((tp-this->ap).length_squared()>.25*ap_diameter*ap_diameter){
	cerr << "wavefront_phase_estimate::tilt_variance_calculation error - "
	     << "three point lies outside of aperture\n";
	cerr << "ap diameter " << ap_diameter << endl << endl;
	tp.print(cerr,"tp ");
	throw(string("wavefront_phase_estimate::tilt_variance_calculation"));
      }

      if(this->stored_emtr==NULL){
	cerr << "wavefront_phase_estimate::tilt_variance_calculation error - "
	     << "emitter is not initialized\n";
	throw(string("wavefront_phase_estimate::tilt_variance_calculation"));
      }

      const vector<tilt_covariance<precision, aperture_type> > & tilt_covariance_ai = 
	tip_tilt_stars ? tt_tilt_covariance_ai : ho_tilt_covariance_ai;
      const vector<tilt_covariance<precision, aperture_type> > & tilt_covariance_ij = 
	tip_tilt_stars ? tt_tilt_covariance_ij : ho_tilt_covariance_ij;
      const vector<double> & weights = 
	tip_tilt_stars ? tip_tilt_weights : high_order_weights;

      int ngs = tilt_covariance_ai.size();

      double val = 
	tilt_covariance_aa.variance(wavelength_meters,
				    nsteps_in_integration,
				    tp);
      for(int i=0; i<ngs; i++)
	val -= 2*weights[i]*
	  tilt_covariance_ai[i].variance(wavelength_meters,
					 nsteps_in_integration,
					 tp);
      int count = 0;
      double fac;
      for(int i=0; i<ngs; i++){
	for(int j=i; j<ngs; j++){
	  fac = i==j ? 1 : 2;
	    val += fac*weights[i]*weights[j]*
	      tilt_covariance_ij[count].variance(wavelength_meters,
						 nsteps_in_integration,
						 tp);
	  count++;
	}
      }

      if(val<-three_frame::precision){
	cerr << "wavefront_phase_estimate::tilt_variance_calculation error - val " 
	     << val 
	     << " less than zero\n";
	cerr << "aa: " 
	     << tilt_covariance_aa.variance(wavelength_meters,
					    nsteps_in_integration,
					    tp)
	     << endl;

	
	for(int i=0; i<ngs; i++)
	  cerr << "a" << i << ": " 
	       << tilt_covariance_ai[i].variance(wavelength_meters,
						 nsteps_in_integration,
						 tp)
	       << endl;
	

	int count = 0;
	double fac;
	for(int i=0; i<ngs; i++){
	  for(int j=i; j<ngs; j++){
	    fac = i==j ? 1 : 2;
	    cerr << i << " " << j << "\t"
		 << tilt_covariance_ij[count].variance(wavelength_meters,
						       nsteps_in_integration,
						       tp)
		 << endl;
	    count++;
	  }
	}
	throw(string("wavefront_phase_estimate::tilt_variance_calculation error\n"));
      }
      /*
      cout << "TILT VARIANCE - tip tilt stars "
	   << tip_tilt_stars << " RETURNING " << val << endl;
      */
      return(val);
    } catch(...){
      cerr << "wavefront_phase_estimate::tilt_variance_calculation error\n";
      throw(string("wavefront_phase_estimate::tilt_variance_calculation"));
    }
  }
  
  template<typename precision, typename aperture_type>
    pixel_array<precision> wavefront_phase_estimate<precision, aperture_type>::
    tilt_variance_calculation(double pixel_scale_meters,
			      double wavelength_meters,
			      int nsteps_in_integration,
			      bool tip_tilt_stars) const {

    try{
      if(this->stored_emtr==NULL){
	cerr << "wavefront_phase_estimate::tilt_variance_calculation error - "
	     << "emitter is not initialized\n";
	throw(string("wavefront_phase_estimate::tilt_variance_calculation"));
      }

      const vector<tilt_covariance<precision, aperture_type> > & tilt_covariance_ai = 
	tip_tilt_stars ? this->tt_tilt_covariance_ai : this->ho_tilt_covariance_ai;
      const vector<tilt_covariance<precision, aperture_type> > & tilt_covariance_ij = 
	tip_tilt_stars ? this->tt_tilt_covariance_ij : this->ho_tilt_covariance_ij;
      const vector<double> & weights = 
	tip_tilt_stars ? this->tip_tilt_weights : this->high_order_weights;

      int ngs = tilt_covariance_ai.size();

      pixel_array<precision> val = 
	tilt_covariance_aa.variance(pixel_scale_meters,
				    wavelength_meters,
				    nsteps_in_integration);

      pixel_array<precision> xxx;
      for(int i=0; i<ngs; i++){
	xxx = 
	  tilt_covariance_ai[i].variance(pixel_scale_meters,
					 wavelength_meters,
					 nsteps_in_integration);
	xxx *= (2*(weights[i]));
	val -= xxx;
      }
      
      int count = 0;
      double fac;
      for(int i=0; i<ngs; i++){
	for(int j=i; j<ngs; j++){
	  fac = i==j ? 1 : 2;
	  xxx = 
	    tilt_covariance_ij[count].variance(pixel_scale_meters,
					       wavelength_meters,
					       nsteps_in_integration);
	  xxx *= (fac*weights[i]*weights[j]);
	  val += xxx;
	  count++;
	}
      }
      return(val);
    } catch(...){
      cerr << "wavefront_phase_estimate::tilt_variance_calculation error\n";
      throw(string("wavefront_phase_estimate::tilt_variance_calculation"));
    }
  }


  template<typename precision, typename aperture_type>
    wavefront_phase_estimate<precision, aperture_type>::
    wavefront_phase_estimate(const refractive_atmospheric_model & ref_atm_model,
			     const aperture_type & ap_type) {

    this->ref_atm_model = ref_atm_model;
    this->stored_emtr = NULL;

    this->stored_OTF_wavelength = -1;

    this->ap = ap_type;

  }

  template<typename precision, typename aperture_type>
    wavefront_phase_estimate<precision, aperture_type>::
    wavefront_phase_estimate(const vector<emitter *> & ttgs,
			     const vector<double> & ttwts,
			     const vector<emitter *> & hogs,
			     const vector<double> & howts,
			     const refractive_atmospheric_model & ref_atm_model,
			     const aperture_type & ap_type) {

    if(ttgs.size()!=ttwts.size()){
      cerr << "wavefront_phase_estimate::wavefront_phase_estimate error -\n"
	   << "unequal number of tip tilt emitters and weights supplied to this function\n";
      throw(string("wavefront_phase_estimate::wavefront_phase_estimate"));
    }
    if(hogs.size()!=howts.size()){
      cerr << "wavefront_phase_estimate::wavefront_phase_estimate error -\n"
	   << "unequal number of high order emitters and weights supplied to this function\n";
      throw(string("wavefront_phase_estimate::wavefront_phase_estimate"));
    }

    tip_tilt_weights = ttwts;
    high_order_weights = howts;

    this->ref_atm_model = ref_atm_model;
    this->stored_emtr = NULL;

    this->stored_OTF_wavelength = -1;

    this->ap = ap_type;

    for(int i=0; i<hogs.size(); i++)
      for(int j=i; j<hogs.size(); j++)
	ho_phase_covariance_ij.push_back(phase_covariance<precision, aperture_type>(*(hogs[i]),
										    *(hogs[j]),
										    ref_atm_model,
										    ap_type));

    for(int i=0; i<ttgs.size(); i++)
      for(int j=i; j<ttgs.size(); j++)
	tt_tilt_covariance_ij.push_back(tilt_covariance<precision, aperture_type>(*(ttgs[i]),
										  *(ttgs[j]),
										  ref_atm_model,
										  ap_type));

    for(int i=0; i<hogs.size(); i++)
      for(int j=i; j<hogs.size(); j++)
	ho_tilt_covariance_ij.push_back(tilt_covariance<precision, aperture_type>(*(hogs[i]),
										  *(hogs[j]),
										  ref_atm_model,
										  ap_type));

  }
  
  template<typename precision, typename aperture_type>
    wavefront_phase_estimate<precision, aperture_type> & 
    wavefront_phase_estimate<precision, aperture_type>::operator=(const wavefront_phase_estimate<precision, aperture_type> & wpe){
    if(this==&wpe)
      return(*this);

    this->ref_atm_model = wpe.ref_atm_model;
    this->ap = wpe.ap;
    this->stored_emtr = emitter::emitter_factory(wpe.stored_emtr);

    tip_tilt_weights = wpe.tip_tilt_weights;
    high_order_weights = wpe.high_order_weights;
    
    phase_covariance_aa = wpe.phase_covariance_aa;
    tilt_covariance_aa = wpe.tilt_covariance_aa;

    ho_phase_covariance_ai = wpe.ho_phase_covariance_ai;
    ho_tilt_covariance_ai = wpe.ho_tilt_covariance_ai;
    tt_tilt_covariance_ai = wpe.tt_tilt_covariance_ai;

    ho_phase_covariance_ij = wpe.ho_phase_covariance_ij;
    ho_tilt_covariance_ij = wpe.ho_tilt_covariance_ij;
    tt_tilt_covariance_ij = wpe.tt_tilt_covariance_ij;

    return(*this);
  }
  
  template<typename precision, typename aperture_type>
    void wavefront_phase_estimate<precision, aperture_type>::read(const char * filename) {
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "wavefront_phase_estimate::read - "
	   << "error opening file " << filename << endl;
      throw(string("wavefront_phase_estimate::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "wavefront_phase_estimate::read - "
	   << "error reading "
	   << this->unique_name() << " from file " 
	   << filename << endl;
      throw(string("wavefront_phase_estimate::read"));
    }
  }

  template<typename precision, typename aperture_type>
    void wavefront_phase_estimate<precision, aperture_type>::read(const iofits & iof) {

    try{
      if(!iof.key_exists("TYPE")){
	cerr << "wavefront_phase_estimate::read error - "
	     << "unrecognized type of file\n";
	throw(string("wavefront_phase_estimate::read"));
      }
      
      string type, comment;
      iof.read_key("TYPE", type, comment);
      if(type!=this->unique_name()){
	cerr << "wavefront_phase_estimate::read error - file of type " 
	     << type << " rather than of type "
	     << this->unique_name() << endl;
	throw(string("wavefront_phase_estimate::read"));
      }

      long nhogs, nttgs;
      iof.read_key("TYPE", type, comment);
      iof.read_key("NTTGS", nttgs, comment);
      tip_tilt_weights.resize(nttgs);

      stringstream ss;
      for(int i=0; i<nttgs; i++){
	ss.str("");
	ss << "TTWT" << i;
	iof.read_key(ss.str().c_str(), tip_tilt_weights[i], comment);
      }

      iof.read_key("NHOGS", nhogs, comment);
      high_order_weights.resize(nhogs);
      for(int i=0; i<nhogs; i++){
	ss.str("");
	ss << "HOWT" << i;
	iof.read_key(ss.str().c_str(), high_order_weights[i], comment);
      }

      tt_tilt_covariance_ij = vector<tilt_covariance<precision, aperture_type> >(nttgs*(nttgs+1)/2);

      iof.movrel_hdu(1);

      this->ap.read(iof);
      this->ref_atm_model.read(iof);

      // Loop over covariances, reading them from subsequent HDUS
      int count = 0;
      for(int i=0; i<nttgs; i++)
	for(int j=i; j<nttgs; j++){
	  tt_tilt_covariance_ij[count].read(iof);
	  count++;
	}

      ho_tilt_covariance_ij = vector<tilt_covariance<precision, aperture_type> >(nhogs*(nhogs+1)/2);
      ho_phase_covariance_ij = vector<phase_covariance<precision, aperture_type> >(nhogs*(nhogs+1)/2);
      count = 0;
      for(int i=0; i<nhogs; i++){
	for(int j=i; j<nhogs; j++){
	  ho_tilt_covariance_ij[count].read(iof);
	  ho_phase_covariance_ij[count].read(iof);
	  count++;
	}
      }

      stored_emtr = NULL;
    } catch(...){
      cerr << "wavefront_phase_estimate::read error - failed on read\n";
      throw(string("wavefront_phase_estimate::read"));
    }
  }

  template<typename precision, typename aperture_type>
    void wavefront_phase_estimate<precision, aperture_type>::write(const char * filename) const {
    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "wavefront_phase_estimate::write - "
	   << "error opening file " << filename << endl;
      throw(string("wavefront_phase_estimate::write"));
    }

    try{this->write(iof);}
    catch(...){
      cerr << "wavefront_phase_estimate::write - "
	   << "error writing "
	   << this->unique_name() << " to file " 
	   << filename << endl;
      throw(string("wavefront_phase_estimate::write"));
    }
  }

  template<typename precision, typename aperture_type>
    void wavefront_phase_estimate<precision, aperture_type>::write(iofits & iof) const {

    try{
      fits_header_data<precision> fhd;
      fhd.write(iof);
      string type = this->unique_name();
      string comment = "object type";
      iof.write_key("TYPE", type, comment);
      
      long nttgs = tip_tilt_weights.size();
      iof.write_key("NTTGS", nttgs, string("number of tip tilt guidestars"));
      
      stringstream ss;
      for(int i=0; i<tip_tilt_weights.size(); i++){
	ss.str("");
	ss << "TTWT" << i;
	iof.write_key(ss.str().c_str(), tip_tilt_weights[i], string("tip tilt guidestar weight"));
      }
      
      long nhogs = high_order_weights.size();
      iof.write_key("NHOGS", nhogs, string("number of high order guidestars"));

      for(int i=0; i<nhogs; i++){
	ss.str("");
	ss << "HOWT" << i;
	iof.write_key(ss.str().c_str(), high_order_weights[i], string("high order guidestar weight"));
      }
      
      this->ap.write(iof);
      ref_atm_model.write(iof);

      // Loop over covariances, writing them to subsequent HDUS
      int count = 0;
      for(int i=0; i<nttgs; i++)
	for(int j=i; j<nttgs; j++){
	  tt_tilt_covariance_ij[count].write(iof);
	  count++;
	}
      
      count = 0;
      for(int i=0; i<nhogs; i++){
	for(int j=i; j<nhogs; j++){
	  ho_tilt_covariance_ij[count].write(iof);
	  ho_phase_covariance_ij[count].write(iof);
	  count++;
	}
      }
    } catch(...){
      cerr << "wavefront_phase_estimate::write error - failed on write\n";
      throw(string("wavefront_phase_estimate::write"));
    }
  }

  template<typename precision, typename aperture_type>
    void wavefront_phase_estimate<precision, aperture_type>::print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield);
    os << prefix << "TYPE       = " << setw(vlspc) << this->unique_name()
       << "/" << "object type" << endl;

    stringstream ss;
    os << prefix << "NTTGS      = " << setw(vlspc) << this->tip_tilt_weights.size()
       << "/" << "object type" << endl;
    for(int i=0; i<tip_tilt_weights.size(); i++){
      ss.str("");
      ss << "TTWT" << i;
      os << prefix << setw(11) << ss.str() << "= " << setw(vlspc) << this->tip_tilt_weights[i]
	 << "/" << "tip tilt guidestar weight" << endl;
    }

    os << prefix << "NHOGS      = " << setw(vlspc) << this->high_order_weights.size()
       << "/" << "object type" << endl;
    for(int i=0; i<high_order_weights.size(); i++){
      ss.str("");
      ss << "HOWT" << i;
      os << prefix << setw(11) << ss.str() << "= " << setw(vlspc) << this->high_order_weights[i]
	 << "/" << "high order guidestar weight" << endl;
    }
    /*    
    this->circ_ap.print(os, prefix);
    this->ref_atm_model.print(os, prefix);
    */

  }

  template<typename precision, typename aperture_type>
    double wavefront_phase_estimate<precision, aperture_type>::
    aperture_averaged_differential_tilt_phase_variance(const emitter & emtr,
						       double wavelength_meters,
						       int nsteps_in_integration) const {

    try{

      this->initialize_emitter(emtr);
      double val = 
	this->aperture_averaged_tilt_variance_calculation(wavelength_meters,
							  nsteps_in_integration,
							  true);
      if(!finite(val) || val<-three_frame::precision){
	cerr << "wavefront_phase_estimate::aperture_averaged_differential_tilt_phase_variance error - val " 
	     << val
	     << " is less than zero\n";
	throw(string("wavefront_phase_estimate::aperture_averaged_differential_tilt_phase_variance"));
      }

      return(val);

    } catch(...) {
      cerr << "wavefront_phase_estimate::aperture_averaged_differential_tilt_phase_variance error\n";
      emtr.print(cerr, "emitter ");
      throw(string("wavefront_phase_estimate::aperture_averaged_differential_tilt_phase_variance"));
    }
  }

  template<typename precision, typename aperture_type>
    double wavefront_phase_estimate<precision, aperture_type>::
    aperture_averaged_differential_phase_variance(const emitter & emtr,
						  double wavelength_meters,
						  int nsteps_in_integration) const {

    try{
      
      this->initialize_emitter(emtr);
      
      // The tip tilt tilt part
      double val = 
	this->aperture_averaged_tilt_variance_calculation(wavelength_meters,
							  nsteps_in_integration,
							  true);
      
      // The high order tilt part
      val -= this->aperture_averaged_tilt_variance_calculation(wavelength_meters,
							       nsteps_in_integration,
							       false);

      int nhogs = ho_tilt_covariance_ai.size();        
      
      // The high order phase variance part
      val +=
	phase_covariance_aa.aperture_averaged_variance(wavelength_meters,
						       nsteps_in_integration);

      for(int i=0; i<nhogs; i++)
	val -= 2*high_order_weights[i]*
	  ho_phase_covariance_ai[i].aperture_averaged_variance(wavelength_meters,
							       nsteps_in_integration);
      
      int count = 0;
      double fac;
      for(int i=0; i<nhogs; i++){
	for(int j=i; j<nhogs; j++){
	  fac = i==j ? 1. : 2.;
	  val += fac*high_order_weights[i]*high_order_weights[j]*
	    ho_phase_covariance_ij[count].aperture_averaged_variance(wavelength_meters,
								     nsteps_in_integration);
	  count++;
	}
      }

      if(!finite(val) || val<-three_frame::precision){
	cerr << "wavefront_phase_estimate::aperture_averaged_differential_phase_variance error - val " 
	     << val
	     << " is less than zero\n";
	throw(string("wavefront_phase_estimate::aperture_averaged_differential_phase_variance"));
      }

      return(val);
    } catch(...) {
      cerr << "wavefront_phase_estimate::aperture_averaged_differential_phase_variance error\n";
      emtr.print(cerr, "emitter ");
      throw(string("wavefront_phase_estimate::aperture_averaged_differential_phase_variance"));
    }
  }

  template<typename precision, typename aperture_type>
    double wavefront_phase_estimate<precision, aperture_type>::
    differential_tilt_phase_variance(const emitter & emtr,
				     const three_point & tp,
				     double wavelength_meters,
				     int nsteps_in_integration) const {
    
    try{
      this->initialize_emitter(emtr);
      return(this->tilt_variance_calculation(wavelength_meters,
					     nsteps_in_integration,
					     tp,
					     true));
    } catch(...) {
      cerr << "wavefront_phase_estimate::differential_tilt_phase_variance error\n";
      throw(string("wavefront_phase_estimate::differential_tilt_phase_variance"));
    }
  }

  template<typename precision, typename aperture_type>
    pixel_array<precision> wavefront_phase_estimate<precision, aperture_type>::
    differential_tilt_phase_variance(const emitter & emtr,
				     double pixel_scale_meters,
				     double wavelength_meters,
				     int nsteps_in_integration) const {
    try{
      this->initialize_emitter(emtr);
      pixel_array<precision> val = 
	this->tilt_variance_calculation(pixel_scale_meters,
					wavelength_meters,
					nsteps_in_integration,
					true);

      double min, max;      
      val.min_and_max(min, max);
      if(!finite(min) || !finite(max) || min<-three_frame::precision){
	cerr << "wavefront_phase_estimate::differential_tilt_phase_variance error - val min " 
	     << min
	     << " is less than zero\n";
	throw(string("wavefront_phase_estimate::differential_tilt_phase_variance"));
      }

      return(val);

    } catch(...) {
      cerr << "wavefront_phase_estimate::differential_tilt_phase_variance error\n";
      emtr.print(cerr, "emitter ");
      throw(string("wavefront_phase_estimate::differential_tilt_phase_variance"));
    }
  }

  template<typename precision, typename aperture_type>
    double wavefront_phase_estimate<precision, aperture_type>::
    differential_phase_variance(const emitter & emtr,
				const three_point & tp,
				double wavelength_meters,
				int nsteps_in_integration) const {
    
    try{
      this->initialize_emitter(emtr);

      // The tip tilt tilt part
      double tt_tilt_variance = 
	this->tilt_variance_calculation(wavelength_meters,
					nsteps_in_integration,
					tp,
					true);
      
      // The high order tilt part
      double ho_tilt_variance = 
	this->tilt_variance_calculation(wavelength_meters,
					nsteps_in_integration,
					tp,
					false);


      double val = tt_tilt_variance - ho_tilt_variance;

      int nhogs = ho_tilt_covariance_ai.size();        

      // The high order phase variance part
      val += 
	phase_covariance_aa.variance(wavelength_meters,
				     nsteps_in_integration,
				     tp);

      for(int i=0; i<nhogs; i++)
	val -= 2*high_order_weights[i]*
	  ho_phase_covariance_ai[i].variance(wavelength_meters,
					     nsteps_in_integration,
					     tp);

      int count = 0;
      double fac;
      for(int i=0; i<nhogs; i++){
	for(int j=i; j<nhogs; j++){
	  fac = i==j ? 1 : 2;
	  val += fac*high_order_weights[i]*high_order_weights[j]*
	    ho_phase_covariance_ij[count].variance(wavelength_meters,
						   nsteps_in_integration,
						   tp);
	  count++;
	}
      }

      if(!finite(val) || val<-three_frame::precision){
	cerr << "wavefront_phase_estimate::differential_phase_variance error 1 - val " 
	     << val
	     << " is less than zero\n";

	cerr << "differential tilt phase variances:  tt " 
	     << tt_tilt_variance
	     << "\tho " << ho_tilt_variance 
	     << endl;

	cerr << "differential high order variance: " 
	     << val - tt_tilt_variance + ho_tilt_variance << endl;
	
	cout << "aa: " 
	     << phase_covariance_aa.variance(wavelength_meters,
					     nsteps_in_integration,
					     tp)
	     << "\t" 
	     << tilt_covariance_aa.variance(wavelength_meters,
					    nsteps_in_integration,
					    tp)
	     << endl;

	for(int i=0; i<nhogs; i++)
	  cout << "a " << i << ": " 
	       << ho_phase_covariance_ai[i].variance(wavelength_meters,
						     nsteps_in_integration,
						     tp)
	       << "\t" 
	       << ho_tilt_covariance_ai[i].variance(wavelength_meters,
						    nsteps_in_integration,
						    tp)
	       << endl;
	
	int count = 0;
	for(int i=0; i<nhogs; i++){
	  for(int j=i; j<nhogs; j++){
	    cout << i << " " << j << ": " 
		 << ho_phase_covariance_ij[count].variance(wavelength_meters,
							   nsteps_in_integration,
							   tp)
		 << "\t" 
		 << ho_tilt_covariance_ij[count].variance(wavelength_meters,
							  nsteps_in_integration,
							  tp)
		 << endl;
	    count++;
	  }
	}

	throw(string("wavefront_phase_estimate::differential_phase_variance"));
      }

      return(val);
    } catch(...) {
      cerr << "wavefront_phase_estimate::differential_phase_variance error\n";
      tp.print(cerr, "three point ");
      emtr.print(cerr, "emitter ");
      throw(string("wavefront_phase_estimate::differential_phase_variance"));
    }
  }


  template<typename precision, typename aperture_type>
    pixel_array<precision> wavefront_phase_estimate<precision, aperture_type>::
    differential_phase_variance(const emitter & emtr,
				double pixel_scale_meters,
				double wavelength_meters,
				int nsteps_in_integration) const {

    try{
      this->initialize_emitter(emtr);

      // The tip tilt tilt part
      pixel_array<precision> val =
	this->tilt_variance_calculation(pixel_scale_meters,
					wavelength_meters,
					nsteps_in_integration,
					true);

      // The high order tilt part
      val -= this->tilt_variance_calculation(pixel_scale_meters,
					     wavelength_meters,
					     nsteps_in_integration,
					     false);

      int nhogs = ho_tilt_covariance_ai.size();        

      // The high order phase variance part
      val += 
	phase_covariance_aa.variance(pixel_scale_meters,
				     wavelength_meters,
				     nsteps_in_integration);

      pixel_array<precision> xxx;

      for(int i=0; i<nhogs; i++){
	xxx = 
	  ho_phase_covariance_ai[i].variance(pixel_scale_meters,
					     wavelength_meters,
					     nsteps_in_integration);
	xxx *= (2*high_order_weights[i]);
	val -= xxx;
      }

      int count = 0;
      double fac;
      for(int i=0; i<nhogs; i++){
	for(int j=i; j<nhogs; j++){
	  fac = i==j ? 1 : 2;
	  xxx = ho_phase_covariance_ij[count].variance(pixel_scale_meters,
						       wavelength_meters,
						       nsteps_in_integration);
	  xxx *= (fac*high_order_weights[i]*high_order_weights[j]);
	  
	  val += xxx;

	  count++;
	}
      }

      double min, max;      
      val.min_and_max(min, max);
      if(!finite(min) || !finite(max) || min<-three_frame::precision){
	cerr << "wavefront_phase_estimate::differential_phase_variance error 2 - val min " 
	     << min
	     << " is less than zero\n";
	throw(string("wavefront_phase_estimate::differential_phase_variance"));
      }

      return(val);
    } catch(...) {
      cerr << "wavefront_phase_estimate::differential_phase_variance error\n";
      emtr.print(cerr, "emitter ");
      throw(string("wavefront_phase_estimate::differential_phase_variance"));
    }
  }

  template<typename precision, typename aperture_type>
    double wavefront_phase_estimate<precision, aperture_type>::
    phase_structure_function(const emitter & emtr,
			     double wavelength_meters,
			     int nsteps_in_integration,
			     const three_point & tp1,
			     const three_point & tp2) const {

    try{
      this->initialize_emitter(emtr);

      
      double val = 
	this->differential_phase_variance(emtr,
					  tp1,
					  wavelength_meters,
					  nsteps_in_integration);

      val += 
	this->differential_phase_variance(emtr,
					  tp2,
					  wavelength_meters,
					  nsteps_in_integration);

      val -= 
	2*phase_covariance_aa.covariance(wavelength_meters,
					 nsteps_in_integration,
					 tp1,
					 tp2);

      int nhogs = ho_tilt_covariance_ai.size();        
      int nttgs = tt_tilt_covariance_ai.size();        

      int count = 0;
      double fac;
      for(int i=0; i<nhogs; i++){
	for(int j=i; j<nhogs; j++){
	  fac = i==j ? 1 : 2;
	  val += 2*fac*high_order_weights[i]*high_order_weights[j]*
	    (ho_tilt_covariance_ij[count].covariance(wavelength_meters,
						     nsteps_in_integration,
						     tp1,
						     tp2) -
	     ho_phase_covariance_ij[count].covariance(wavelength_meters,
						      nsteps_in_integration,
						      tp1,
						      tp2));
	  count++;
	}
      }

      count = 0;
      for(int i=0; i<nttgs; i++){
	for(int j=i; j<nttgs; j++){
	  fac = i==j ? 1 : 2;
	  val -= 2*fac*tip_tilt_weights[i]*tip_tilt_weights[j]*
	    tt_tilt_covariance_ij[count].covariance(wavelength_meters,
						    nsteps_in_integration,
						    tp1,
						    tp2);

	  count++;
	}
      }

      for(int i=0; i<nhogs; i++){
	val += 2*high_order_weights[i]*
	  (ho_phase_covariance_ai[i].covariance(wavelength_meters,
						nsteps_in_integration,
						tp1,
						tp2) -
	   ho_tilt_covariance_ai[i].covariance(wavelength_meters,
					       nsteps_in_integration,
					       tp1,
					       tp2) +
	   ho_phase_covariance_ai[i].covariance(wavelength_meters,
						nsteps_in_integration,
						tp2,
						tp1) -
	   ho_tilt_covariance_ai[i].covariance(wavelength_meters,
					       nsteps_in_integration,
					       tp2,
					       tp1));
      }

      for(int i=0; i<nttgs; i++){
	val += 2*tip_tilt_weights[i]*
	  (tt_tilt_covariance_ai[i].covariance(wavelength_meters,
					       nsteps_in_integration,
					       tp1,
					       tp2) +
	   tt_tilt_covariance_ai[i].covariance(wavelength_meters,
					       nsteps_in_integration,
					       tp2,
					       tp1));
      }

      if(!finite(val) || val<-three_frame::precision){
      //if(!finite(val) || val<-1000*three_frame::precision){
	cerr << "wavefront_phase_estimate::phase_structure_function 1 error - val " 
	     << val
	     << " is less than zero\n";
	throw(string("wavefront_phase_estimate::phase_structure_function 1"));
      }

      return(val);
    } catch(...) {
      cerr << "wavefront_phase_estimate::phase_structure_function 1 error\n";
      emtr.print(cerr, "emitter ");
      throw(string("wavefront_phase_estimate::phase_structure_function"));
    }
  }
 
  template<typename precision, typename aperture_type>
    pixel_array<precision> wavefront_phase_estimate<precision, aperture_type>::
    phase_structure_function(const emitter & emtr,
			     double pixel_scale_meters,
			     double wavelength_meters,
			     int nsteps_in_integration,
			     const three_point & tp) const {
     
    try{
      this->initialize_emitter(emtr);

      pixel_array<precision> val = 
	this->differential_phase_variance(emtr,
					  pixel_scale_meters,
					  wavelength_meters,
					  nsteps_in_integration);



      double tmp = this->differential_phase_variance(emtr,
						     tp,
						     wavelength_meters,
						     nsteps_in_integration);
      add_val_to_pixarr(val,
			tmp,
			2*pixel_scale_meters/this->get_aperture_outer_diameter());


      pixel_array<precision> xxx = 
	phase_covariance_aa.covariance(pixel_scale_meters,
				       wavelength_meters,
				       nsteps_in_integration,
				       tp);

      xxx *= 2.;
      val -= xxx;

      int nhogs = ho_tilt_covariance_ai.size();        
      int nttgs = tt_tilt_covariance_ai.size();        

      int count = 0;
      double fac;
      for(int i=0; i<nhogs; i++){
	for(int j=i; j<nhogs; j++){
	  fac = i==j ? 1 : 2;
	  xxx = 
	    (ho_tilt_covariance_ij[count].covariance(pixel_scale_meters,
						     wavelength_meters,
						     nsteps_in_integration,
						     tp) -
	     ho_phase_covariance_ij[count].covariance(pixel_scale_meters,
						      wavelength_meters,
						      nsteps_in_integration,
						      tp));
	  xxx *= (2*fac*high_order_weights[i]*high_order_weights[j]);
	  val += xxx;
	  count++;
	}
      }

      count = 0;
      for(int i=0; i<nttgs; i++){
	for(int j=i; j<nttgs; j++){
	  fac = i==j ? 1 : 2;
	  xxx = 
	    tt_tilt_covariance_ij[count].covariance(pixel_scale_meters,
						    wavelength_meters,
						    nsteps_in_integration,
						    tp);
	  xxx *= (2*fac*tip_tilt_weights[i]*tip_tilt_weights[j]);
	  val -= xxx;

	  count++;
	}
      }

      for(int i=0; i<nhogs; i++){
	xxx = 
	  (ho_phase_covariance_ai[i].covariance(pixel_scale_meters,
						wavelength_meters,
						nsteps_in_integration,
						tp,
						true) -
	   ho_tilt_covariance_ai[i].covariance(pixel_scale_meters,
					       wavelength_meters,
					       nsteps_in_integration,
					       tp,
					       true) +
	   ho_phase_covariance_ai[i].covariance(pixel_scale_meters,
						wavelength_meters,
						nsteps_in_integration,
						tp, 
						false) -
	   ho_tilt_covariance_ai[i].covariance(pixel_scale_meters,
					       wavelength_meters,
					       nsteps_in_integration,
					       tp,
					       false));
	xxx *= (2*high_order_weights[i]);
	val += xxx;
      }

      for(int i=0; i<nttgs; i++){
	xxx = 
	  (tt_tilt_covariance_ai[i].covariance(pixel_scale_meters,
					       wavelength_meters,
					       nsteps_in_integration,
					       tp,
					       true) +
	   tt_tilt_covariance_ai[i].covariance(pixel_scale_meters,
					       wavelength_meters,
					       nsteps_in_integration,
					       tp,
					       false));
	xxx *= (2*tip_tilt_weights[i]);
	val += xxx;
      }

      double min, max;      
      val.min_and_max(min, max);
      if(!finite(min) || !finite(max) || min<-three_frame::precision){
	//if(!finite(min) || !finite(max) || min<-1000*three_frame::precision){
	cerr << "wavefront_phase_estimate::phase_structure_function 2 error - val min " 
	     << min
	     << " is less than zero\n";

	/*
	structure_function strfn(val, 1e-6);
	strfn.write("strfn_2_error.fits");
	*/

	throw(string("wavefront_phase_estimate::phase_structure_function 2"));
      }

      return(val);
    } catch(...) {
      cerr << "wavefront_phase_estimate::phase_structure_function 2 error\n";
      emtr.print(cerr, "emitter ");
      throw(string("wavefront_phase_estimate::phase_structure_function"));
    }
  }

  template<typename precision, typename aperture_type>
    pixel_array<precision> wavefront_phase_estimate<precision, aperture_type>::
    phase_structure_function(const emitter & emtr,
			     double pixel_scale_meters,
			     double wavelength_meters,
			     int nsteps_in_integration,
			     int xindex,
			     int yindex) const {

    try{
      this->initialize_emitter(emtr);

      pixel_array<precision> val = 
	this->differential_phase_variance(emtr,
					  pixel_scale_meters,
					  wavelength_meters,
					  nsteps_in_integration);


      
      vector<long> axes = val.get_axes();
      double x_halfpix=0, y_halfpix=0;
      if(axes[1]%2==0){
	x_halfpix = .5;
      }
      if(axes[0]%2==0){
	y_halfpix = .5;
      }

      three_point tp((xindex+x_halfpix)*pixel_scale_meters,
		     (yindex+y_halfpix)*pixel_scale_meters,
		     0,
		     this->phase_covariance_aa.get_aperture());

      double tmp = this->differential_phase_variance(emtr,
						     tp,
						     wavelength_meters,
						     nsteps_in_integration);

      add_val_to_pixarr(val,
			tmp,
			2*pixel_scale_meters/this->get_aperture_outer_diameter());
    

      pixel_array<precision> xxx = 
	phase_covariance_aa.covariance(pixel_scale_meters,
				       wavelength_meters,
				       nsteps_in_integration,
				       xindex,
				       yindex);

      xxx *= 2;
      val -= xxx;

      int nhogs = ho_tilt_covariance_ai.size();        
      int nttgs = tt_tilt_covariance_ai.size();        

      int count = 0;
      double fac;
      for(int i=0; i<nhogs; i++){
	for(int j=i; j<nhogs; j++){
	  fac = i==j ? 1 : 2;

	  xxx = 
	    (ho_tilt_covariance_ij[count].covariance(pixel_scale_meters,
						     wavelength_meters,
						     nsteps_in_integration,
						     xindex,
						     yindex) -
	     ho_phase_covariance_ij[count].covariance(pixel_scale_meters,
						      wavelength_meters,
						      nsteps_in_integration,
						      xindex,
						      yindex));

	  xxx *= (2*fac*high_order_weights[i]*high_order_weights[j]);
	  val += xxx;
	  count++;
	}
      }


      count = 0;
      for(int i=0; i<nttgs; i++){
	for(int j=i; j<nttgs; j++){
	  fac = i==j ? 1 : 2;
	  //val -= //XXX2*fac*tip_tilt_weights[i]*tip_tilt_weights[j]*
	  xxx = 
	    tt_tilt_covariance_ij[count].covariance(pixel_scale_meters,
						    wavelength_meters,
						    nsteps_in_integration,
						    xindex,
						    yindex);
	  xxx *= (2*fac*tip_tilt_weights[i]*tip_tilt_weights[j]);
	  val -= xxx;
	  count++;
	}
      }

      for(int i=0; i<nhogs; i++){
	//val += //XXX2*high_order_weights[i]*
	xxx = 
	  (ho_phase_covariance_ai[i].covariance(pixel_scale_meters,
						wavelength_meters,
						nsteps_in_integration,
						xindex,
						yindex,
						true) -
	   ho_tilt_covariance_ai[i].covariance(pixel_scale_meters,
					       wavelength_meters,
					       nsteps_in_integration,
					       xindex,
					       yindex,
					       true) +
	   ho_phase_covariance_ai[i].covariance(pixel_scale_meters,
						wavelength_meters,
						nsteps_in_integration,
						xindex,
						yindex,
						false) -
	   ho_tilt_covariance_ai[i].covariance(pixel_scale_meters,
					       wavelength_meters,
					       nsteps_in_integration,
					       xindex,
					       yindex,
					       false));
	xxx *= (2*high_order_weights[i]);
	val += xxx;
      }

      for(int i=0; i<nttgs; i++){
	xxx = 
	  (tt_tilt_covariance_ai[i].covariance(pixel_scale_meters,
					       wavelength_meters,
					       nsteps_in_integration,
					       xindex,
					       yindex,
					       true) +
	   tt_tilt_covariance_ai[i].covariance(pixel_scale_meters,
					       wavelength_meters,
					       nsteps_in_integration,
					       xindex,
					       yindex,
					       false));
	xxx *= (2*tip_tilt_weights[i]);
	val += xxx;
      }

      double min, max;      
      val.min_and_max(min, max);
      if(!finite(min) || !finite(max) || min<-1e-6){
      //if(!finite(min) || !finite(max) || min<-three_frame::precision){
	//if(!finite(min) || !finite(max) || min*min<-three_frame::precision){
	cerr << "wavefront_phase_estimate::phase_structure_function 3 error - val min " 
	     << min
	     << " is less than zero\n";

	/*
	structure_function strfn(val, 1e-6);
	strfn.write("strfn_3_error.fits");
	*/

	throw(string("wavefront_phase_estimate::phase_structure_function 3"));
      }

      return(val);
    } catch(...) {
      cerr << "wavefront_phase_estimate::phase_structure_function 3 error\n";
      emtr.print(cerr, "emitter ");
      throw(string("wavefront_phase_estimate::phase_structure_function"));
    }
  }

  // specialized version for 2 NGS to take advantage of stationarity
  // of the structure function over the pupil plane.
  template<typename precision, typename aperture_type>
    basic_otf<precision> wavefront_phase_estimate<precision, aperture_type>::
    NGS_optical_transfer_function(const emitter & emtr,
				  double wavelength_meters,
				  int nsteps_in_integration,
				  double pupil_plane_pixel_scale_meters,
				  double cutoff) const {


    if(ho_phase_covariance_ij.size()!=1 || tt_tilt_covariance_ij.size()!=1){
      cerr << "wavefront_phase_estimate::NGS_optical_transfer_function error "
	   << " - this function can only operate on phase estimates generated from a single plane wave emitter\n";
      throw(string("wavefront_phase_estimate::NGS_optical_transfer_function"));
    }      

    try{dynamic_cast<const plane_wave_emitter &>(ho_phase_covariance_ij[0].get_first_emitter());}
    catch(...){
      cerr << "wavefront_phase_estimate::NGS_optical_transfer_function error "
	   << " - this function can only operate on phase estimates generated from a single plane wave emitter\n";
      throw(string("wavefront_phase_estimate::NGS_optical_transfer_function"));      
    }

    try{dynamic_cast<const plane_wave_emitter &>(tt_tilt_covariance_ij[0].get_first_emitter());}
    catch(...){
      cerr << "wavefront_phase_estimate::NGS_optical_transfer_function error "
	   << " - this function can only operate on phase estimates generated from a single plane wave emitter\n";
      throw(string("wavefront_phase_estimate::NGS_optical_transfer_function"));      
    }

    if(dynamic_cast<const plane_wave_emitter &>(ho_phase_covariance_ij[0].get_first_emitter()) !=
       dynamic_cast<const plane_wave_emitter &>(tt_tilt_covariance_ij[0].get_first_emitter())){
      cerr << "wavefront_phase_estimate::NGS_optical_transfer_function error "
	   << " - this function can only operate on phase estimates generated from a single plane wave emitter\n";
      throw(string("wavefront_phase_estimate::NGS_optical_transfer_function"));      
    }

    try{dynamic_cast<const plane_wave_emitter &>(emtr);}
    catch(...){
      cerr << "wavefront_phase_estimate::NGS_optical_transfer_function error "
	   << " - argument not of type plane wave emitter\n";
      emtr.print(cerr, "arg emitter");
      throw(string("wavefront_phase_estimate::NGS_optical_transfer_function"));
    }

    try{

      // Check arguments
      if(wavelength_meters<=0){
	cerr << "wavefront_phase_estimate::optical_transfer_function error - wavelength "
	     << wavelength_meters
	     << " out of range\n";
	throw(string("wavefront_phase_estimate::optical_transfer_function"));
      }

      if(nsteps_in_integration<=0){
	cerr << "wavefront_phase_estimate::optical_transfer_function error - "
	     << "number of steps in numerical integration "
	     << nsteps_in_integration
	     << " out of range\n";
	throw(string("wavefront_phase_estimate::optical_transfer_function"));
      }

      if(pupil_plane_pixel_scale_meters<=0){
	cerr << "wavefront_phase_estimate::optical_transfer_function error - pupil plane pixel scale "
	     << pupil_plane_pixel_scale_meters
	     << " out of range\n";
	throw(string("wavefront_phase_estimate::optical_transfer_function"));
      }

      if(cutoff<=0){
	cerr << "wavefront_phase_estimate::optical_transfer_function error - cutoff "
	     << cutoff
	     << " out of range\n";
	throw(string("wavefront_phase_estimate::optical_transfer_function"));
      }

      // Compute the OTF
      double aperture_diameter_meters = this->get_aperture_outer_diameter();
      vector<long> strfn_axes(2,(long)ceil(aperture_diameter_meters/pupil_plane_pixel_scale_meters));
      pixel_array<precision> phase_structure_function_pixarr(strfn_axes);
      int nsf_elems = phase_structure_function_pixarr.total_space();

      // otf array always has odd dimensions
      int otf_extrapix=1;
      double otf_halfpix = 0;

      vector<long> otf_axes(2,2*strfn_axes[0]+1);
      int notf_elems = otf_axes[0]*otf_axes[1];

      // If we've already initialized the OTF, return it
      if(stored_OTF_wavelength==wavelength_meters &&
	 stored_OTF.total_space()==notf_elems)
	return(basic_otf<precision>(0,
				    0,
				    pupil_plane_pixel_scale_meters,
				    wavelength_meters,
				    stored_OTF,
				    pixel_phase_array<precision>(stored_OTF.get_axes())));

      precision * data;
      try{
	data = new precision[notf_elems];
      } catch(...) {
	cerr << "wavefront_phase_estimate::optical_transfer_function - error allocating memory\n";
	throw(string("wavefront_phase_estimate::optical_transfer_function"));
      }


      three_frame model_tf = this->ref_atm_model.get_three_frame();
      double secant_zenith_angle;
      double max_range_meters;
      double min_range_meters;
      three_vector little_omega;

      this->ref_atm_model.Tyler_get_constants(emtr,
					      ho_phase_covariance_ij[0].get_first_emitter(),
					      model_tf,
					      aperture_diameter_meters,
					      secant_zenith_angle,
					      max_range_meters,
					      min_range_meters,
					      little_omega);

      
      vector<double> layer_heights = ref_atm_model.get_layer_heights();
      vector<power_spectrum *> pspectra = ref_atm_model.get_power_spectra();
      vector<double> cn2_coeffs(pspectra.size());
      int npspectra = pspectra.size();
      for(int i=0; i<npspectra; i++){

	isotropic_power_law_spectrum<power_law, null_inner_scale> * isospec = 
	  dynamic_cast<isotropic_power_law_spectrum<power_law, null_inner_scale> *>(pspectra[i]);
	if(!isospec){
	  cerr << "wavefront_phase_estimate::NGS_optical_transfer_function error - ref atm model power spectrum "
	       << i 
	       << " not an isotropic power_law_spectrum\n";
	  throw(string("wavefront_phase_estimate::NGS_optical_transfer_function"));
	}

	if(fabs(isospec->get_power_law().get_exponent()+11/3.)>three_frame::precision){
	  cerr << "wavefront_phase_estimate::NGS_optical_transfer_function error - ref atm model power spectrum "
	       << i 
	       << " has a non-Komolgorov power law "
	       << isospec->get_power_law().get_exponent()
	       << endl;
	  throw(string("wavefront_phase_estimate::NGS_optical_transfer_function"));
	}

	double fac = 2*M_PI*5*gamma_function(5/6.)/pow(2,4/3.)/pow(M_PI,3/2.)/9./gamma_function(2/3.);
	cn2_coeffs[i] = pspectra[i]->get_coefficient()/fac;
	delete pspectra[i];
      }

      double sum_omega = 0;
      three_vector omega;
      for(int i=0; i<npspectra; i++){
	omega = little_omega*layer_heights[i];
	sum_omega += cn2_coeffs[i]*
	  2*pow(omega.length(),5/3.);
      }

      // m,n label possible values of the x and y components of the vector x
      int index;
      double sum_dr;
      double sum_dr_plus_omega;
      double sum_dr_minus_omega;
      double limit = aperture_diameter_meters/pupil_plane_pixel_scale_meters;
      double val;
      double tmp_omega_x, tmp_omega_y;
      double prefac = 2*get_Xi()*pow(aperture_diameter_meters, 5/3.)*4*M_PI*M_PI/wavelength_meters/wavelength_meters;
      double alpha;
      for(int m=-otf_axes[1]/2; m<otf_axes[1]/2+otf_extrapix; m++){
	for(int n=-otf_axes[0]/2; n<otf_axes[0]/2+otf_extrapix; n++){

	  index = (m+otf_axes[1]/2)*otf_axes[0]+n+otf_axes[0]/2;
	    
	  if((m*m+n*n - limit*limit)>-three_frame::precision){
	    data[index] = 0;
	  } else {

	    sum_dr = sum_dr_plus_omega = sum_dr_minus_omega = 0;
	    for(int i=0; i<npspectra; i++){
	      sum_dr += cn2_coeffs[i]*2*pow((m*m+n*n)*16/(double)(otf_axes[0]*otf_axes[0]),5/6.);
	      
	      tmp_omega_x = little_omega.x(model_tf)*layer_heights[i];
	      tmp_omega_y = little_omega.y(model_tf)*layer_heights[i];
	      sum_dr_plus_omega += 
		cn2_coeffs[i]*pow(((m*4/(double)otf_axes[0]+tmp_omega_x)*(m*4/(double)otf_axes[0]+tmp_omega_x)+
				   (n*4/(double)otf_axes[0]+tmp_omega_y)*(n*4/(double)otf_axes[0]+tmp_omega_y)),5/6.);
	      sum_dr_minus_omega += 
		cn2_coeffs[i]*pow(((m*4/(double)otf_axes[0]-tmp_omega_x)*(m*4/(double)otf_axes[0]-tmp_omega_x)+
				   (n*4/(double)otf_axes[0]-tmp_omega_y)*(n*4/(double)otf_axes[0]-tmp_omega_y)),5/6.);
	      
	    }
	    
	    val = prefac*(sum_omega + sum_dr - sum_dr_plus_omega - sum_dr_minus_omega);

	    alpha = sqrt(m*m+n*n)*pupil_plane_pixel_scale_meters/aperture_diameter_meters;
	    data[index] = exp(-.5*val)*(2/M_PI)*(acos(alpha) - alpha*sqrt(1-alpha*alpha));
	  }
	}
      }

      // normalize array, either using the peak of the OTF or
      // by constraining total power to be unity
      val = data[notf_elems/2];
      for(int i=0; i<notf_elems; i++){
	data[i] /= val;
      }

      stored_OTF = pixel_array<precision>(otf_axes,data);
      stored_OTF_wavelength = wavelength_meters;

      pixel_amp_array<precision> tmp_amp(otf_axes,data);
      pixel_phase_array<precision> tmp_phase(otf_axes);

      delete [] data;

      return(basic_otf<precision>(0,
				  0,
				  pupil_plane_pixel_scale_meters,
				  wavelength_meters,
				  tmp_amp,
				  tmp_phase));
    } catch(...) {
      cerr << "wavefront_phase_estimate::NGS_optical_transfer_function error\n";
      throw(string("wavefront_phase_estimate::NGS_optical_transfer_function"));
    }
  }

  template<typename precision, typename aperture_type>
    basic_otf<precision> 
    wavefront_phase_estimate<precision, aperture_type>::
    optical_transfer_function(const emitter & emtr,
			      double wavelength_meters,
			      int nsteps_in_integration,
			      double pupil_plane_pixel_scale_meters,
			      double cutoff) const {
    

    try{

      // Check arguments
      if(wavelength_meters<=0){
	cerr << "wavefront_phase_estimate::optical_transfer_function error - wavelength "
	     << wavelength_meters
	     << " out of range\n";
	throw(string("wavefront_phase_estimate::optical_transfer_function"));
      }

      if(nsteps_in_integration<=0){
	cerr << "wavefront_phase_estimate::optical_transfer_function error - "
	     << "number of steps in numerical integration "
	     << nsteps_in_integration
	     << " out of range\n";
	throw(string("wavefront_phase_estimate::optical_transfer_function"));
      }

      if(pupil_plane_pixel_scale_meters<=0){
	cerr << "wavefront_phase_estimate::optical_transfer_function error - pupil plane pixel scale "
	     << pupil_plane_pixel_scale_meters
	     << " out of range\n";
	throw(string("wavefront_phase_estimate::optical_transfer_function"));
      }

      if(cutoff<=0){
	cerr << "wavefront_phase_estimate::optical_transfer_function error - cutoff "
	     << cutoff
	     << " out of range\n";
	throw(string("wavefront_phase_estimate::optical_transfer_function"));
      }

      // Initialize emitter
      try{
	this->initialize_emitter(emtr);
      } catch(...) {
	cerr << "wavefront_phase_estimate::optical_transfer_function - "
	     << "error initializing emitter\n";
	throw(string("wavefront_phase_estimate::optical_transfer_function"));
      }


      if(ho_phase_covariance_ij.size()==1 && 
	 tt_tilt_covariance_ij.size()==1 &&
	 dynamic_cast<const plane_wave_emitter *>(&(ho_phase_covariance_ij[0].get_first_emitter())) &&
	 dynamic_cast<const plane_wave_emitter *>(&(tt_tilt_covariance_ij[0].get_first_emitter())) &&
	 (dynamic_cast<const plane_wave_emitter &>(ho_phase_covariance_ij[0].get_first_emitter()) ==
	  dynamic_cast<const plane_wave_emitter &>(tt_tilt_covariance_ij[0].get_first_emitter())) &&
	 dynamic_cast<const plane_wave_emitter *>(&emtr))
	return(this->NGS_optical_transfer_function(emtr,
						   wavelength_meters,
						   nsteps_in_integration,
						   pupil_plane_pixel_scale_meters,
						   cutoff));


      // Compute the OTF
      double aperture_diameter_meters = this->get_aperture_outer_diameter();
      vector<long> strfn_axes(2,(long)ceil(aperture_diameter_meters/pupil_plane_pixel_scale_meters));
      pixel_array<precision> phase_structure_function_pixarr(strfn_axes);
      int nsf_elems = phase_structure_function_pixarr.total_space();

      double x_halfpix=0, y_halfpix=0;
      int x_extrapix=1, y_extrapix=1;
      if(strfn_axes[1]%2==0){
	x_halfpix = .5;
	x_extrapix = 0;
      }
      if(strfn_axes[0]%2==0){
	y_halfpix = .5;
	y_extrapix = 0;
      }
    
      // otf array always has odd dimensions
      int otf_extrapix=1;
      double otf_halfpix = 0;

      int index, index2;
      double val;
      double limit = pow(aperture_diameter_meters/2./pupil_plane_pixel_scale_meters,2.);

      vector<long> otf_axes(2,2*strfn_axes[0]+1);
      int notf_elems = otf_axes[0]*otf_axes[1];

      // If we've already initialized the OTF, return it
      if(stored_OTF_wavelength==wavelength_meters &&
	 stored_OTF.total_space()==notf_elems)
	return(basic_otf<precision>(0,
				    0,
				    pupil_plane_pixel_scale_meters,
				    wavelength_meters,
				    stored_OTF,
				    pixel_phase_array<precision>(stored_OTF.get_axes())));

      precision * data, *strfn_data;
      try{
	data = new precision[notf_elems];
	strfn_data = new precision[nsf_elems];
      } catch(...) {
	cerr << "wavefront_phase_estimate::optical_transfer_function - error allocating memory\n";
	throw(string("wavefront_phase_estimate::optical_transfer_function"));
      }
      for(int i=0; i<notf_elems; i++)
	data[i] = 0;

      int half_strfn_axes_0 = strfn_axes[0]/2;
      int half_strfn_axes_1 = strfn_axes[1]/2;
      int half_otf_axes_0 = otf_axes[0]/2;
      int half_otf_axes_1 = otf_axes[1]/2;

      // i,j label possible locations of the three point r1 in the pupil plane
      for(int i=-half_strfn_axes_1; i<half_strfn_axes_1+x_extrapix; i++){
	for(int j=-half_strfn_axes_0; j<half_strfn_axes_0+y_extrapix; j++){

	  // This is the test on whether r1 lies within the pupil
	  if(((i+x_halfpix)*(i+x_halfpix)+(j+y_halfpix)*(j+y_halfpix))>limit)
	    continue;
	
	  try{
	    phase_structure_function_pixarr = 
	      this->phase_structure_function(emtr,
					     pupil_plane_pixel_scale_meters,
					     wavelength_meters,
					     nsteps_in_integration,
					     i,
					     j);
	  } catch(...) {
	    cerr << "wavefront_phase_estimate::optical_transfer_function error - could not get structure function\n";
	    throw(string("wavefront_phase_estimate::optical_transfer_function"));
	  }

	  // copy out of pixel array for optimization purposes in the loop below
	  for(int x=0; x<nsf_elems; x++)
	    strfn_data[x] = phase_structure_function_pixarr.data(x);

	  // m,n label possible values of the x and y components of the vector x
	  for(int m=-half_otf_axes_1; m<half_otf_axes_1+otf_extrapix; m++){
	    for(int n=-half_otf_axes_0; n<half_otf_axes_0+otf_extrapix; n++){
	    
	      // This is the test on whether r2 lies within the pupil
	      val = (i+x_halfpix+m-otf_halfpix)*(i+x_halfpix+m-otf_halfpix)+
		(j+y_halfpix+n-otf_halfpix)*(j+y_halfpix+n-otf_halfpix);
	      if(val>limit)
		continue;

	      index = (m+half_otf_axes_1)*otf_axes[0]+n+half_otf_axes_0;
	      index2 = (i + half_strfn_axes_1 + m)*strfn_axes[0]+(j + half_strfn_axes_0 + n);

	      val = strfn_data[index2];

	      if(val<cutoff)
		data[index] += exp(-.5*val);
	    }	
	  }
	}
      }

      // normalize array, either using the peak of the OTF or
      // by constraining total power to be unity
      double data_max = -DBL_MAX;
      for(int i=0; i<notf_elems; i++){
	if(data_max<data[i]) data_max=data[i];
      }

      for(int i=0; i<notf_elems; i++)
	data[i] /= data_max;

      stored_OTF = pixel_array<precision>(otf_axes,data);
      stored_OTF_wavelength = wavelength_meters;

      pixel_amp_array<precision> tmp_amp(otf_axes,data);
      pixel_phase_array<precision> tmp_phase(otf_axes);

      delete [] data;
      delete [] strfn_data;

      return(basic_otf<precision>(0,
				  0,
				  pupil_plane_pixel_scale_meters,
				  wavelength_meters,
				  tmp_amp,
				  tmp_phase));
    } catch(...) {
      cerr << "wavefront_phase_estimate::optical_transfer_function error\n";
      throw(string("wavefront_phase_estimate::optical_transfer_function"));
    }
  }

  template<typename precision, typename aperture_type>
    basic_observation<precision> wavefront_phase_estimate<precision, aperture_type>::
    private_point_spread_function(double field_size_arcsecs,
				  double oversampling_factor,
				  const basic_otf<precision> & long_exposure_OTF) const {
    
    try{
      // Check arguments
      if(field_size_arcsecs<=0){
	cerr << "wavefront_phase_estimate::private_point_spread_function error - field size "
	     << field_size_arcsecs
	     << " out of range\n";
	throw(string("wavefront_phase_estimate::private_point_spread_function"));
      }

      if(oversampling_factor<=0){
	cerr << "wavefront_phase_estimate::private_point_spread_function error - oversampling factor "
	     << oversampling_factor
	     << " out of range\n";
	throw(string("wavefront_phase_estimate::private_point_spread_function"));
      }


      // set up the storage
      vector<long> otf_axes = long_exposure_OTF.get_axes();

      // Define parameters for the Goertzel Reinsch propagator
      double rad_to_arcsec = 180*3600/M_PI;
      double distance = 1e10;
      double initial_pixel_scale = long_exposure_OTF.get_pixel_scale();

      // The otf axes are always 2*(pupil_plane_axes)+1.  Here we convert 
      // back to pupil plane axes in computing the nyquist pixel scale
      double nyquist_pixel_scale = 
	long_exposure_OTF.get_wavelength() * distance/ (double) (otf_axes[0]/2) / initial_pixel_scale;

      double final_pixel_scale = nyquist_pixel_scale/oversampling_factor;
      double final_pixel_scale_arcsecs = rad_to_arcsec * final_pixel_scale / distance;
    
      vector<long> psf_axes = vector<long>(2,(int)ceil(field_size_arcsecs/final_pixel_scale_arcsecs));
    
      int notf_elems = otf_axes[0]*otf_axes[1];
      int npsf_elems = psf_axes[0]*psf_axes[1];
      precision *otf_data;
      precision *psf_data;    
      try{
	otf_data = new precision[2*notf_elems];
	psf_data = new precision[2*npsf_elems];
      } catch(...) {
	cerr << "wavefront_phase_estimate::private_point_spread_function - error allocating memory\n";
	throw(string("wavefront_phase_estimate::private_point_spread_function"));
      }


      // Slope correction for even dimension arrays
      double xslope = psf_axes[1]%2==1 ? 0 : M_PI*final_pixel_scale/nyquist_pixel_scale/(double)(otf_axes[1]/2);
      double yslope = psf_axes[0]%2==1 ? 0 : M_PI*final_pixel_scale/nyquist_pixel_scale/(double)(otf_axes[0]/2);


      // Halfpixel information
      // Note - right now the OTF axes always have odd dimensions
      double otf_x_halfpix=0, otf_y_halfpix=0;
      int otf_x_extrapix=1, otf_y_extrapix=1;
      if(otf_axes[1]%2==0){
	otf_x_halfpix = .5;
	otf_x_extrapix = 0;
      }
      if(otf_axes[0]%2==0){
	otf_y_halfpix = .5;
	otf_y_extrapix = 0;
      }
      double psf_x_halfpix=0, psf_y_halfpix=0;
      int psf_x_extrapix=1, psf_y_extrapix=1;
      if(psf_axes[1]%2==0){
	psf_x_halfpix = .5;
	psf_x_extrapix = 0;
      }
      if(psf_axes[0]%2==0){
	psf_y_halfpix = .5;
	psf_y_extrapix = 0;
      }

      // Initialize the complex OTF data array
      int index;
      double twopi = 2*M_PI;
      double amp, phase;
      for(int i=-otf_axes[1]/2; i<otf_axes[1]/2+otf_x_extrapix; i++){
	for(int j=-otf_axes[0]/2; j<otf_axes[0]/2+otf_y_extrapix; j++){
	  index = (i+otf_axes[1]/2)*otf_axes[0]+j+otf_axes[0]/2;
	  amp = abs(long_exposure_OTF.data(index));
	  phase = arg(long_exposure_OTF.data(index)) + 
	    fmod(-xslope*(i+otf_x_halfpix) - yslope*(j+otf_y_halfpix), twopi);
	  otf_data[2*index] = amp*cos(phase);
	  otf_data[2*index+1] = amp*sin(phase);
	}
      }

      double sampling_factor = -initial_pixel_scale*final_pixel_scale / long_exposure_OTF.get_wavelength() / distance;
      goertzel_reinsch_transform(otf_axes, psf_axes, sampling_factor, otf_data, psf_data);

      delete [] otf_data;

      // No need to fix up the slopes here - we just want the amplitude
      for(int i=0; i<npsf_elems; i++)
	psf_data[i] = sqrt(psf_data[2*i]*psf_data[2*i] + psf_data[2*i+1]*psf_data[2*i+1]);

      pixel_amp_array<precision> long_exposure_PSF(psf_axes,psf_data);
      delete [] psf_data;
      return(basic_observation<precision>(0,
					  0,
					  final_pixel_scale/distance,
					  long_exposure_OTF.get_wavelength(),
					  long_exposure_PSF));
    } catch(...) {
      cerr << "wavefront_phase_estimate::private_point_spread_function error\n";
      throw(string("wavefront_phase_estimate::private_point_spread_function"));
    }
  }


  template<typename precision, typename aperture_type>
    basic_observation<precision> wavefront_phase_estimate<precision, aperture_type>::
    point_spread_function(const emitter & emtr,
			  double wavelength_meters,
			  int nsteps_in_integration,
			  double pupil_plane_pixel_scale_meters,
			  double field_size_arcsecs,
			  double oversampling_factor,
			  double cutoff) const {

    if(nsteps_in_integration<=0){
      cerr << "wavefront_phase_estimate::point_spread_function error - number of steps in numerical integration "
	   << nsteps_in_integration
	   << " out of range\n";
      throw(string("wavefront_phase_estimate::point_spread_function"));
    }
    
    basic_otf<precision> long_exposure_OTF = 
      this->optical_transfer_function(emtr,
				      wavelength_meters,
				      nsteps_in_integration,
				      pupil_plane_pixel_scale_meters,
				      cutoff);

    return(this->private_point_spread_function(field_size_arcsecs,
					       oversampling_factor,
					       long_exposure_OTF));
    
  }

  template<typename precision, typename aperture_type>
    basic_observation<precision> wavefront_phase_estimate<precision, aperture_type>::
    NGS_point_spread_function(const emitter & emtr,
			      double wavelength_meters,
			      int nsteps_in_integration,
			      double pupil_plane_pixel_scale_meters,
			      double field_size_arcsecs,
			      double oversampling_factor,
			      double cutoff) const {
    
    if(nsteps_in_integration<=0){
      cerr << "wavefront_phase_estimate::NGS_point_spread_function error - number of steps in numerical integration "
	   << nsteps_in_integration
	   << " out of range\n";
      throw(string("wavefront_phase_estimate::NGS_point_spread_function"));
    }

    basic_otf<precision> long_exposure_OTF = 
      this->NGS_optical_transfer_function(emtr,
					  wavelength_meters,
					  nsteps_in_integration,
					  pupil_plane_pixel_scale_meters,
					  cutoff);

    return(this->private_point_spread_function(field_size_arcsecs,
					       oversampling_factor,
					       long_exposure_OTF));
    
  }
}
#endif
