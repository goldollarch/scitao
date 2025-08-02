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

#ifndef COVARIANCE_H
#define COVARIANCE_H

#include <string>
#include <iostream>
#include <vector>
#include "AO_sim_base.h"
#include "pixel_array.h"
#include "aperture.h"
#include "refractive_atmosphere.h"

namespace Arroyo {
  class iofits;
}

namespace Arroyo {

  using std::string;

  ///    
  ///  A class to represent the phase covariance
  ///  between two emitters
  ///

  template<typename precision, typename aperture_type>
    class phase_covariance :
    public AO_sim_base {

    private:

    phase_covariance(){
      cerr << "phase_covariance::phase_covariance error - unsupported\n";
      exit(-1);
    }
  };


  ///    
  ///  A class to represent the phase covariance between two emitters.
  ///  Template specialization for circular apertures
  ///

  template<typename precision>
    class phase_covariance<precision, circular_aperture> :
    public AO_sim_base {
    
    private:

    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("circular phase covariance"));};

    protected:

    refractive_atmospheric_model ref_atm_model;
    circular_aperture circ_ap;

    mutable double pixel_scale_meters;
    mutable int nsteps_in_integration;

    mutable pixel_array<precision> stored_caliph_A;
    mutable pixel_array<precision> stored_caliph_B;
    mutable precision stored_caliph_D;
    mutable precision stored_caliph_M;

    emitter * emtr_a;
    emitter * emtr_b;

    static double Xi;

    void initialize_integration(int nsteps_in_integration) const;

    void initialize_pixel_scale(double pixel_scale_meters) const;

    public:

    ///////////////////////////////////////////
    ///  Null constructor
    phase_covariance(){
      emtr_a = NULL;
      emtr_b = NULL;
    };

    ///////////////////////////////////////////
    ///  Constructor
    phase_covariance(const emitter & emtr1,
		     const emitter & emtr2,
		     const refractive_atmospheric_model & ref_atm_model,
		     const circular_aperture & circ_ap);
    
    ///////////////////////////////////////////
    ///  Destructor
    ~phase_covariance(){
      delete emtr_a;
      delete emtr_b;
    };

    ///////////////////////////////////////////
    ///  Copy constructor
    phase_covariance(const phase_covariance & pcv){
      this->operator=(pcv);
    };

    ///////////////////////////////////////////
    ///  Construct from file
    phase_covariance(const char * filename){
      this->read(filename);
    };

    ///////////////////////////////////////////
    ///  Construct from an iofits object
    phase_covariance(const iofits & iof){
      this->read(iof);
    };

    ///////////////////////////////////////////
    ///  Operator =
    phase_covariance & operator=(const phase_covariance & pcv);

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
    ///  Get a reference to the refractive atmospheric model
    const refractive_atmospheric_model & get_model() const {
      return(this->ref_atm_model);
    };

    ///////////////////////////////////////////
    ///  Get a reference to the aperture
    circular_aperture get_aperture() const {
      return(this->circ_ap);
    };

    ///////////////////////////////////////////
    ///  Get a reference to the first emitter
    const emitter & get_first_emitter() const {
      return(*(this->emtr_a));
    };

    ///////////////////////////////////////////
    ///  Get a reference to the second emitter
    const emitter & get_second_emitter() const {
      return(*(this->emtr_b));
    };

    ///////////////////////////////////////////
    ///  Function to return the aperture averaged variance
    double aperture_averaged_variance(double wavelength_meters,
				      int nsteps_in_integration) const;

    ///////////////////////////////////////////
    ///  Function to return the variance at the
    ///  point tp
    double variance(double wavelength_meters,
		    int nsteps_in_integration,
		    const three_point & tp) const;

    ///////////////////////////////////////////
    ///  Function to return a pixel array containing the variances
    pixel_array<precision> variance(double pixel_scale_meters,
				    double wavelength_meters,
				    int nsteps_in_integration) const;
    
    ///////////////////////////////////////////
    ///  Function to return the covariance between the
    ///  points tp1 and tp2
    double covariance(double wavelength_meters,
		      int nsteps_in_integration,
		      const three_point & tp1,
		      const three_point & tp2) const;

    ///////////////////////////////////////////
    ///  Function to return a pixel array containing the covariances
    ///  between the point tp and every other point in the aperture
    pixel_array<precision> covariance(double pixel_scale_meters,
				      double wavelength_meters,
				      int nsteps_in_integration,
				      const three_point & tp,
				      bool first_arg = true) const;

    ///////////////////////////////////////////
    ///  Function to return a pixel array containing the covariances
    ///  between the point with index (xindex,yindex) and every other
    ///  point in the aperture
    pixel_array<precision> covariance(double pixel_scale_meters,
				      double wavelength_meters,
				      int nsteps_in_integration,
				      int xindex,
				      int yindex,
				      bool first_arg = true) const;

  };

  /*
  ///    
  ///  A class to represent the phase covariance between two emitters.
  ///  Template specialization for annular apertures
  ///

  template<typename precision>
    class phase_covariance<precision, annular_aperture> :
    public AO_sim_base {
    
    private:

    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("annular phase covariance"));};

    protected:


    public:

    ///////////////////////////////////////////
    ///  Null constructor
    phase_covariance(){};

    ///////////////////////////////////////////
    ///  Constructor
    phase_covariance(const emitter & emtr1,
		     const emitter & emtr2,
		     const refractive_atmospheric_model & ref_atm_model,
		     const annular_aperture & annular_ap);
    
    ///////////////////////////////////////////
    ///  Destructor
    ~phase_covariance(){};

    ///////////////////////////////////////////
    ///  Copy constructor
    phase_covariance(const phase_covariance & pcv){
      this->operator=(pcv);
    };

    ///////////////////////////////////////////
    ///  Construct from file
    phase_covariance(const char * filename){
      this->read(filename);
    };

    ///////////////////////////////////////////
    ///  Construct from an iofits object
    phase_covariance(const iofits & iof){
      this->read(iof);
    };

    ///////////////////////////////////////////
    ///  Operator =
    phase_covariance & operator=(const phase_covariance & pcv);

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
    ///  Get a reference to the refractive atmospheric model
    const refractive_atmospheric_model & get_model() const {
      return(this->outer_aperture_phase_covariance.get_model());
    };

    ///////////////////////////////////////////
    ///  Get a reference to the aperture
    annular_aperture get_aperture() const {
      annular_aperture annular_ap(this->inner_aperture_phase_covariance.get_aperture().get_diameter(),
				  this->outer_aperture_phase_covariance.get_aperture().get_diameter());
      annular_ap.three_frame::operator=(this->outer_aperture_phase_covariance.get_aperture());
      return(annular_ap);
    };

    ///////////////////////////////////////////
    ///  Get a reference to the first emitter
    const emitter & get_first_emitter() const {
      return(this->outer_aperture_phase_covariance.get_first_emitter());
    };

    ///////////////////////////////////////////
    ///  Get a reference to the second emitter
    const emitter & get_second_emitter() const {
      return(this->outer_aperture_phase_covariance.get_second_emitter());
    };

    ///////////////////////////////////////////
    ///  Function to return the aperture averaged variance
    double aperture_averaged_variance(double wavelength_meters,
				      int nsteps_in_integration) const;

    ///////////////////////////////////////////
    ///  Function to return the variance at the
    ///  point tp
    double variance(double wavelength_meters,
		    int nsteps_in_integration,
		    const three_point & tp) const;

    ///////////////////////////////////////////
    ///  Function to return a pixel array containing the variances
    pixel_array<precision> variance(double pixel_scale_meters,
				    double wavelength_meters,
				    int nsteps_in_integration) const;
    
    ///////////////////////////////////////////
    ///  Function to return the covariance between the
    ///  points tp1 and tp2
    double covariance(double wavelength_meters,
		      int nsteps_in_integration,
		      const three_point & tp1,
		      const three_point & tp2) const;

    ///////////////////////////////////////////
    ///  Function to return a pixel array containing the covariances
    ///  between the point tp and every other point in the aperture
    pixel_array<precision> covariance(double pixel_scale_meters,
				      double wavelength_meters,
				      int nsteps_in_integration,
				      const three_point & tp,
				      bool first_arg = true) const;

    ///////////////////////////////////////////
    ///  Function to return a pixel array containing the covariances
    ///  between the point with index (xindex,yindex) and every other
    ///  point in the aperture
    pixel_array<precision> covariance(double pixel_scale_meters,
				      double wavelength_meters,
				      int nsteps_in_integration,
				      int xindex,
				      int yindex,
				      bool first_arg = true) const;

  };
  */


  ///    
  ///  A class to represent the tilt covariance
  ///  between two emitters
  ///

  template<typename precision, typename aperture_type>
    class tilt_covariance :
    public AO_sim_base {

    private:

    tilt_covariance(){
      cerr << "tilt_covariance::tilt_covariance error - unsupported\n";
      exit(-1);
    }
  };


  ///    
  ///  A class to represent the phase covariance between two emitters.
  ///  Template specialization for circular apertures
  ///

  template<typename precision>
    class tilt_covariance<precision, circular_aperture> :
    public AO_sim_base {
    
    private:

    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("circular tilt covariance"));};

    protected:

    refractive_atmospheric_model ref_atm_model;
    circular_aperture circ_ap;

    mutable double pixel_scale_meters;
    mutable int nsteps_in_integration;

    mutable pixel_array<precision> stored_caliph_G;
    mutable pixel_array<precision> stored_caliph_H;
    mutable pixel_array<precision> stored_caliph_I;
    mutable pixel_array<precision> stored_caliph_J;
    mutable precision stored_caliph_E;
    mutable precision stored_caliph_F;
    mutable precision stored_caliph_F_bar;
    mutable precision stored_caliph_K;
    mutable precision stored_caliph_L;

    emitter * emtr_a;
    emitter * emtr_b;

    three_vector omega_hat;

    static double Xi;

    void initialize_integration(int nsteps_in_integration) const;

    void initialize_pixel_scale(double pixel_scale_meters) const;

    public:

    ///////////////////////////////////////////
    ///  Null constructor
    tilt_covariance(){
      emtr_a = NULL;
      emtr_b = NULL;
    };

    ///////////////////////////////////////////
    ///  Constructor
    tilt_covariance(const emitter & emtr1,
		    const emitter & emtr2,
		    const refractive_atmospheric_model & ref_atm_model,
		    const circular_aperture & circ_ap);
    
    ///////////////////////////////////////////
    ///  Destructor
    ~tilt_covariance(){
      delete emtr_a;
      delete emtr_b;
    };

    ///////////////////////////////////////////
    ///  Copy constructor
    tilt_covariance(const tilt_covariance & tcv){
      this->operator=(tcv);
    };

    ///////////////////////////////////////////
    ///  Construct from file
    tilt_covariance(const char * filename){
      this->read(filename);
    };

    ///////////////////////////////////////////
    ///  Construct from an iofits object
    tilt_covariance(const iofits & iof){
      this->read(iof);
    };

    ///////////////////////////////////////////
    ///  Operator =
    tilt_covariance & operator=(const tilt_covariance & tcv);

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
    ///  Get a reference to the refractive atmospheric model
    const refractive_atmospheric_model & get_model() const {
      return(this->ref_atm_model);
    };

    ///////////////////////////////////////////
    ///  Get a reference to the aperture
    circular_aperture get_aperture() const {
      return(this->circ_ap);
    };

    ///////////////////////////////////////////
    ///  Get a reference to the first emitter
    const emitter & get_first_emitter() const {
      return(*(this->emtr_a));
    };

    ///////////////////////////////////////////
    ///  Get a reference to the second emitter
    const emitter & get_second_emitter() const {
      return(*(this->emtr_b));
    };

    ///////////////////////////////////////////
    ///  Function to return the aperture averaged variance
    double aperture_averaged_variance(double wavelength_meters,
				      int nsteps_in_integration) const;

    ///////////////////////////////////////////
    ///  Function to return the variance at the
    ///  point tp
    double variance(double wavelength_meters,
		    int nsteps_in_integration,
		    const three_point & tp) const;

    ///////////////////////////////////////////
    ///  Function to return a pixel array containing the variances
    pixel_array<precision> variance(double pixel_scale_meters,
				    double wavelength_meters,
				    int nsteps_in_integration) const;
    
    ///////////////////////////////////////////
    ///  Function to return the covariance between the
    ///  points tp1 and tp2
    double covariance(double wavelength_meters,
		      int nsteps_in_integration,
		      const three_point & tp1,
		      const three_point & tp2) const;

    ///////////////////////////////////////////
    ///  Function to return a pixel array containing the covariances
    ///  between the point tp and every other point in the aperture
    pixel_array<precision> covariance(double pixel_scale_meters,
				      double wavelength_meters,
				      int nsteps_in_integration,
				      const three_point & tp,
				      bool first_arg = true) const;

    ///////////////////////////////////////////
    ///  Function to return a pixel array containing the covariances
    ///  between the point with index (xindex,yindex) and every other
    ///  point in the aperture
    pixel_array<precision> covariance(double pixel_scale_meters,
				      double wavelength_meters,
				      int nsteps_in_integration,
				      int xindex,
				      int yindex,
				      bool first_arg = true) const;

  };

  /*
  ///    
  ///  A class to represent the phase covariance between two emitters.
  ///  Template specialization for annular apertures
  ///

  template<typename precision>
    class tilt_covariance<precision, annular_aperture> :
    public AO_sim_base {
    
    private:

    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("annular tilt covariance"));};

    protected:
    
    tilt_covariance<precision, circular_aperture> outer_aperture_tilt_covariance;
    tilt_covariance<precision, circular_aperture> inner_aperture_tilt_covariance;

    public:

    ///////////////////////////////////////////
    ///  Null constructor
    tilt_covariance(){};

    ///////////////////////////////////////////
    ///  Constructor
    tilt_covariance(const emitter & emtr1,
		    const emitter & emtr2,
		    const refractive_atmospheric_model & ref_atm_model,
		    const annular_aperture & annular_ap);
    
    ///////////////////////////////////////////
    ///  Destructor
    ~tilt_covariance(){};

    ///////////////////////////////////////////
    ///  Copy constructor
    tilt_covariance(const tilt_covariance & tcv){
      this->operator=(tcv);
    };

    ///////////////////////////////////////////
    ///  Construct from file
    tilt_covariance(const char * filename){
      this->read(filename);
    };

    ///////////////////////////////////////////
    ///  Construct from an iofits object
    tilt_covariance(const iofits & iof){
      this->read(iof);
    };

    ///////////////////////////////////////////
    ///  Operator =
    tilt_covariance & operator=(const tilt_covariance & tcv);

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
    ///  Get a reference to the refractive atmospheric model
    const refractive_atmospheric_model & get_model() const {
      return(this->outer_aperture_tilt_covariance.get_model());
    };

    ///////////////////////////////////////////
    ///  Get a reference to the aperture
    annular_aperture get_aperture() const {
      annular_aperture annular_ap(this->inner_aperture_tilt_covariance.get_aperture().get_diameter(),
				  this->outer_aperture_tilt_covariance.get_aperture().get_diameter());
      annular_ap.three_frame::operator=(this->outer_aperture_tilt_covariance.get_aperture());
      return(annular_ap);
    };

    ///////////////////////////////////////////
    ///  Get a reference to the first emitter
    const emitter & get_first_emitter() const {
      return(this->outer_aperture_tilt_covariance.get_first_emitter());
    };

    ///////////////////////////////////////////
    ///  Get a reference to the second emitter
    const emitter & get_second_emitter() const {
      return(this->outer_aperture_tilt_covariance.get_second_emitter());
    };

    ///////////////////////////////////////////
    ///  Function to return the aperture averaged variance
    double aperture_averaged_variance(double wavelength_meters,
				      int nsteps_in_integration) const;

    ///////////////////////////////////////////
    ///  Function to return the variance at the
    ///  point tp
    double variance(double wavelength_meters,
		    int nsteps_in_integration,
		    const three_point & tp) const;

    ///////////////////////////////////////////
    ///  Function to return a pixel array containing the variances
    pixel_array<precision> variance(double pixel_scale_meters,
				    double wavelength_meters,
				    int nsteps_in_integration) const;
    
    ///////////////////////////////////////////
    ///  Function to return the covariance between the
    ///  points tp1 and tp2
    double covariance(double wavelength_meters,
		      int nsteps_in_integration,
		      const three_point & tp1,
		      const three_point & tp2) const;

    ///////////////////////////////////////////
    ///  Function to return a pixel array containing the covariances
    ///  between the point tp and every other point in the aperture
    pixel_array<precision> covariance(double pixel_scale_meters,
				      double wavelength_meters,
				      int nsteps_in_integration,
				      const three_point & tp,
				      bool first_arg = true) const;

    ///////////////////////////////////////////
    ///  Function to return a pixel array containing the covariances
    ///  between the point with index (xindex,yindex) and every other
    ///  point in the aperture
    pixel_array<precision> covariance(double pixel_scale_meters,
				      double wavelength_meters,
				      int nsteps_in_integration,
				      int xindex,
				      int yindex,
				      bool first_arg = true) const;

  };
  */

  namespace {

    template <typename precision> 
      void add_val_to_pixarr(pixel_array<precision> & pixarr,
			     precision & val,
			     double normalization_factor){

      vector<long> axes = pixarr.get_axes();
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
      three_frame tf;
      double val2;
      int index;
      three_vector pixel_vector;	
      for(int m=-axes[1]/2; m<axes[1]/2+x_extrapix; m++){
	for(int n=-axes[0]/2; n<axes[0]/2+y_extrapix; n++){
	  pixel_vector = three_vector((m+x_halfpix)*normalization_factor,
				      (n+y_halfpix)*normalization_factor,
				      0,
				      tf);
	  index = (m+axes[1]/2)*axes[0]+n+axes[0]/2;
	  if(pixel_vector.length_squared()<=1){
	    val2 = pixarr.data(index);
	    pixarr.set_data(index, val+val2);
	  }
	}
      }
    }

    void get_omega_hat(const emitter & emtr_a,
		       const emitter & emtr_b,
		       const three_frame & tf,
		       three_vector & omega_hat) {

      three_vector emtr_a_unit_vector = 
	emtr_a.get_emission_vector(static_cast<const three_point>(tf));
      three_vector emtr_b_unit_vector = 
	emtr_b.get_emission_vector(static_cast<const three_point>(tf));

      double emtr_a_range_meters = 1e300;
      double emtr_b_range_meters = 1e300;

      const spherical_wave_emitter * swe;
      if((swe=dynamic_cast<const spherical_wave_emitter *>(&emtr_a)))
	emtr_a_range_meters = (*swe - tf).length();
      if((swe=dynamic_cast<const spherical_wave_emitter *>(&emtr_b)))
	emtr_b_range_meters = (*swe - tf).length();

      three_vector max_range_unit_vector, min_range_unit_vector;
      if(emtr_a_range_meters>emtr_b_range_meters){
	max_range_unit_vector = emtr_a_unit_vector;
	min_range_unit_vector = emtr_b_unit_vector;
      } else {
	max_range_unit_vector = emtr_b_unit_vector;
	min_range_unit_vector = emtr_a_unit_vector;
      }

      if((max_range_unit_vector - min_range_unit_vector).length()<three_frame::precision){
	omega_hat = three_vector();
      } else {
	omega_hat = (max_range_unit_vector - min_range_unit_vector);
	omega_hat *= (1/omega_hat.length());
      }
    }
  }

  template<typename precision>
  double phase_covariance<precision, circular_aperture>::Xi = get_Xi();

  template<typename precision>
    void phase_covariance<precision, circular_aperture>::
    initialize_integration(int nsteps_in_integration) const {

    try{

      if(this->emtr_a==NULL || this->emtr_b==NULL){
	cerr << "phase_covariance::initialize_integration error - uninitialized emitters\n";
	throw(string("phase_covariance::initialize_integration"));
      }

      if(this->nsteps_in_integration==nsteps_in_integration)
	return;
      
      this->nsteps_in_integration = nsteps_in_integration;    
      stored_caliph_D = this->ref_atm_model.caliph_D(*(this->emtr_a),
						     *(this->emtr_b),
						     this->circ_ap.get_diameter(),
						     this->nsteps_in_integration);
    } catch(...) {
      cerr << "phase_covariance::initialize_integration error\n";
      this->emtr_a->print(cerr, "\temtra");
      this->emtr_b->print(cerr, "\temtrb");
      cout << "ap diameter " << this->circ_ap.get_diameter()
	   << "\tnsteps in integration " << this->nsteps_in_integration
	   << endl;
      throw(string("phase_covariance::initialize_integration"));
    }
  }
  
  template<typename precision>
    void phase_covariance<precision, circular_aperture>::
    initialize_pixel_scale(double pixel_scale_meters) const {

    try{

      if(this->emtr_a==NULL || this->emtr_b==NULL){
	cerr << "phase_covariance::initialize_pixel_scale error - uninitialized emitters\n";
	throw(string("phase_covariance::initialize_pixel_scale"));
      }

      if(this->pixel_scale_meters==pixel_scale_meters)
	return;
      
      this->pixel_scale_meters = pixel_scale_meters;
      stored_caliph_A = 
	this->ref_atm_model.template caliph_A<precision>(*(this->emtr_a),
							 *(this->emtr_b),
							 this->circ_ap.get_diameter(),
							 this->pixel_scale_meters);
      stored_caliph_B = 
	this->ref_atm_model.template caliph_B<precision>(*(this->emtr_a),
							 *(this->emtr_b),
							 this->circ_ap.get_diameter(),
							 this->pixel_scale_meters);
    } catch(...) {
      cerr << "phase_covariance::initialize_pixel_scale error\n";
      this->emtr_a->print(cerr, "\temtra");
      this->emtr_b->print(cerr, "\temtrb");
      cout << "ap diameter " << this->circ_ap.get_diameter()
	   << "\tnsteps in integration " << this->pixel_scale_meters
	   << endl;
      throw(string("phase_covariance::initialize_pixel_scale"));
    }
  }

  template<typename precision>
    phase_covariance<precision, circular_aperture>::
    phase_covariance(const emitter & emtr_a,
		     const emitter & emtr_b,
		     const refractive_atmospheric_model & ref_atm_model,
		     const circular_aperture & circ_ap) {

    this->emtr_a = emitter::emitter_factory(&emtr_a);
    this->emtr_b = emitter::emitter_factory(&emtr_b);

    this->ref_atm_model = ref_atm_model;

    this->circ_ap = circ_ap;

    stored_caliph_D = -1;
    stored_caliph_M = -1;
    pixel_scale_meters = -1;
    nsteps_in_integration = -1;
  }

  template<typename precision>
    phase_covariance<precision, circular_aperture> &
    phase_covariance<precision, circular_aperture>::operator=(const phase_covariance & pcv) {

    this->emtr_a = emitter::emitter_factory(pcv.emtr_a);
    this->emtr_b = emitter::emitter_factory(pcv.emtr_b);

    this->ref_atm_model = pcv.ref_atm_model;
    this->circ_ap = pcv.circ_ap;
    
    this->pixel_scale_meters = pcv.pixel_scale_meters;
    this->nsteps_in_integration = pcv.nsteps_in_integration;

    this->stored_caliph_A = pcv.stored_caliph_A;
    this->stored_caliph_B = pcv.stored_caliph_B;
    this->stored_caliph_D = pcv.stored_caliph_D;
    this->stored_caliph_M = pcv.stored_caliph_M;
  }

  template<typename precision>
    void phase_covariance<precision, circular_aperture>::read(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "phase_covariance::read - "
	   << "error opening file " << filename << endl;
      throw(string("phase_covariance::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "phase_covariance::read - "
	   << "error reading "
	   << this->unique_name() << " from file " 
	   << filename << endl;
      throw(string("phase_covariance::read"));
    }
  }
  
  template<typename precision>
    void phase_covariance<precision, circular_aperture>::read(const Arroyo::iofits & iof){

    string type, comment;
    if (!iof.key_exists("TYPE"))
      {
	cerr << this->unique_name() << "::read error - unrecognized "
	     << "file type" << endl;
	throw(this->unique_name() + string("::read"));
      }
    iof.read_key("TYPE",type,comment);
    if(type!=this->unique_name()){
      cerr << this->unique_name() << "::read error - file of type: " 
	   << type << " rather than type " << this->unique_name() << endl;
      throw(this->unique_name() + string("::read"));
    }

    // Move to the next HDU
    iof.movrel_hdu(1);

    try{
      this->emtr_a = emitter::emitter_factory(iof);
    } catch(...){
      cerr << "phase_covariance::read error - could not read first emitter\n";
      throw(string("phase_covariance::read"));
    }

    try{
      this->emtr_b = emitter::emitter_factory(iof);
    } catch(...){
      cerr << "phase_covariance::read error - could not read second emitter\n";
      throw(string("phase_covariance::read"));
    }

    try{
      this->circ_ap.read(iof);
    } catch(...){
      cerr << "phase_covariance::read error - could not read aperture\n";
      throw(string("phase_covariance::read"));
    }

    try{
      this->ref_atm_model.read(iof);
    } catch(...){
      cerr << "phase_covariance::read error - could not read refractive atmospheric model\n";
      throw(string("phase_covariance::read"));
    }

  }

  template<typename precision>
    void phase_covariance<precision, circular_aperture>::write(const char * filename) const {

    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "phase_covariance::write - "
	   << "error opening file " << filename << endl;
      throw(string("phase_covariance::write"));
    }

    try{this->write(iof);}
    catch(...){
      cerr << "phase_covariance::write - "
	   << "error writing "
	   << this->unique_name() << " to file " 
	   << filename << endl;
      throw(string("phase_covariance::write"));
    }
  }

  template<typename precision>
    void phase_covariance<precision, circular_aperture>::write(Arroyo::iofits & iof) const {

    fits_header_data<double> tmphdr;
    tmphdr.write(iof);

    string type, comment;

    type = this->unique_name();
    comment = "object type";
    iof.write_key("TYPE", type, comment);

    try{this->emtr_a->write(iof);}
    catch(...){
      cerr << "phase_covariance::write error - could not write first emitter\n";
      throw(string("phase_covariance::write"));
    }
      
    try{this->emtr_b->write(iof);}
    catch(...){
      cerr << "phase_covariance::write error - could not write second emitter\n";
      throw(string("phase_covariance::write"));
    }

    try{this->circ_ap.write(iof);}
    catch(...){
      cerr << "phase_covariance::write error - could not write aperture\n";
      throw(string("phase_covariance::write"));
    }

    try{this->ref_atm_model.write(iof);}
    catch(...){
      cerr << "phase_covariance::write error - could not write refractive atmospheric model\n";
      throw(string("phase_covariance::write"));
    }
  }

  template<typename precision>
    void phase_covariance<precision, circular_aperture>::print(std::ostream & os, 
					    const char * prefix) const {

    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << this->unique_name()
       << "/" << "object type" << endl;
    this->emtr_a->print(os, prefix);
    this->emtr_b->print(os, prefix);
    this->circ_ap.print(os, prefix);
    this->ref_atm_model.print(os, prefix);
  }

  template<typename precision>
    double phase_covariance<precision, circular_aperture>::aperture_averaged_variance(double wavelength_meters,
										      int nsteps_in_integration) const {

    try{

      this->initialize_integration(nsteps_in_integration);
      if(stored_caliph_M==-1)
	stored_caliph_M = this->ref_atm_model.caliph_M(*(this->emtr_a),
							*(this->emtr_b),
							this->circ_ap.get_diameter());
      
      return(Xi*pow(this->circ_ap.get_diameter(),5/3.)*
	     4*M_PI*M_PI/wavelength_meters/wavelength_meters*
	     (stored_caliph_D - stored_caliph_M));
    } catch(...){
      cerr << "phase_covariance::aperture_averaged_variance error - could not form variance\n";
      throw(string("phase_covariance::aperture_averaged_variance"));
    }      

  }

  template<typename precision>
    double phase_covariance<precision, circular_aperture>::variance(double wavelength_meters,
								    int nsteps_in_integration,
								    const three_point & tp) const {

    try{
      double val;
      double ap_diameter_meters = this->circ_ap.get_diameter();
      
      this->initialize_integration(nsteps_in_integration);
      
      three_vector pixel_vector = tp - static_cast<three_point>(this->circ_ap);
      pixel_vector *= (2/ap_diameter_meters);
      if((pixel_vector.length()-1)>-three_frame::precision){
	cerr << "phase_covariance::variance error - "
	     << "three_point lies outside of aperture\n";
	throw(string("phase_covariance::variance"));
      }
      
      
      val = this->ref_atm_model.caliph_A(*(this->emtr_a),
					 *(this->emtr_b), 
					 ap_diameter_meters,
					 pixel_vector);
      
      val += this->ref_atm_model.caliph_B(*(this->emtr_a), 
					  *(this->emtr_b), 
					  ap_diameter_meters,
					  pixel_vector);
      
      val -= this->ref_atm_model.caliph_C(*(this->emtr_a), 
					  *(this->emtr_b), 
					  ap_diameter_meters,
					  pixel_vector,
					  pixel_vector);

      val -= stored_caliph_D;

      return(val*Xi*pow(ap_diameter_meters,5/3.)*4*M_PI*M_PI/wavelength_meters/wavelength_meters);
    } catch(...) {
      cerr << "phase_covariance::variance error - could not form variance\n";
      throw(string("phase_covariance::variance"));
    }      
  }

  template<typename precision>
  pixel_array<precision> 
  phase_covariance<precision, circular_aperture>::variance(double pixel_scale_meters,
							   double wavelength_meters,
							   int nsteps_in_integration) const {
    
    try{
      this->initialize_integration(nsteps_in_integration);
      this->initialize_pixel_scale(pixel_scale_meters);
      
      double ap_diameter_meters = this->circ_ap.get_diameter();
      
      pixel_array<precision> variance_pixarr = this->stored_caliph_A;
      variance_pixarr += this->stored_caliph_B;
      variance_pixarr -= this->ref_atm_model.template caliph_C<precision>(*(emtr_a),
									  *(emtr_b),
									  ap_diameter_meters,
									  pixel_scale_meters);
      
      precision tmp = -1*this->stored_caliph_D;

      add_val_to_pixarr(variance_pixarr,
			tmp,
			2*pixel_scale_meters/ap_diameter_meters);

      variance_pixarr *= Xi*pow(ap_diameter_meters,5/3.)*4*M_PI*M_PI/wavelength_meters/wavelength_meters;

      return(variance_pixarr);

    } catch(...) {
      cerr << "phase_covariance::variance error - could not form variance\n";
      throw(string("phase_covariance::variance"));
    }      
  }

  template<typename precision>
    double phase_covariance<precision, circular_aperture>::covariance(double wavelength_meters,
								      int nsteps_in_integration,
								      const three_point & tp1,
								      const three_point & tp2) const {
    try {

      this->initialize_integration(nsteps_in_integration);
      double val;
      double ap_diameter_meters = this->circ_ap.get_diameter();

      this->initialize_integration(nsteps_in_integration);

      three_vector pixel_vector_1 = tp1 - static_cast<three_point>(this->circ_ap);
      three_vector pixel_vector_2 = tp2 - static_cast<three_point>(this->circ_ap);
      pixel_vector_1 *= (2/ap_diameter_meters);
      pixel_vector_2 *= (2/ap_diameter_meters);

      if((pixel_vector_1.length()-1)>-three_frame::precision){
	cerr << "phase_covariance::covariance error - "
	     << "first three_point lies outside of aperture\n";
	throw(string("phase_covariance::covariance"));
      }
      if((pixel_vector_2.length()-1)>-three_frame::precision){
	cerr << "phase_covariance::covariance error - "
	     << "second three_point lies outside of aperture\n";
	throw(string("phase_covariance::covariance"));
      }

      val = this->ref_atm_model.caliph_A(*(this->emtr_a), 
					  *(this->emtr_b), 
					  ap_diameter_meters,
					  pixel_vector_1);
    
      val += this->ref_atm_model.caliph_B(*(this->emtr_a), 
					   *(this->emtr_b), 
					   ap_diameter_meters,
					   pixel_vector_2);
      
      val -= this->ref_atm_model.caliph_C(*(this->emtr_a), 
					   *(this->emtr_b), 
					   ap_diameter_meters,
					   pixel_vector_1,
					   pixel_vector_2);

      val -= stored_caliph_D;

      return(val*Xi*pow(ap_diameter_meters,5/3.)*4*M_PI*M_PI/wavelength_meters/wavelength_meters);
    } catch(...) {
      cerr << "phase_covariance::covariance error - could not form covariance\n";
      throw(string("phase_covariance::covariance"));
    }      
  }

  template<typename precision>
    pixel_array<precision> phase_covariance<precision, circular_aperture>::covariance(double pixel_scale_meters,
										      double wavelength_meters,
										      int nsteps_in_integration,
										      const three_point & tp,
										      bool first_arg) const {
    
    try {
      pixel_array<precision> covariance_pixarr;

      this->initialize_integration(nsteps_in_integration);
      this->initialize_pixel_scale(pixel_scale_meters);

      if(first_arg){
	covariance_pixarr = this->stored_caliph_B;
      } else {
	covariance_pixarr = this->stored_caliph_A;
      }
    
      vector<long> axes = covariance_pixarr.get_axes();
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

      double ap_diameter_meters = this->circ_ap.get_diameter();
      double normalization_factor = 2*pixel_scale_meters/ap_diameter_meters;
      int index;
      double val;
      double fac = Xi*pow(ap_diameter_meters,5/3.)*4*M_PI*M_PI/wavelength_meters/wavelength_meters;
      three_vector pixel_vector;
      three_vector normalized_tv = (tp - this->circ_ap)*(2/this->circ_ap.get_diameter());
      if((normalized_tv.length()-1)>-three_frame::precision){
	cerr << "phase_covariance::covariance error - "
	     << "three_point lies outside of aperture\n";
	throw(string("phase_covariance::covariance"));
      }

      for(int m=-axes[1]/2; m<axes[1]/2+x_extrapix; m++){
	for(int n=-axes[0]/2; n<axes[0]/2+y_extrapix; n++){
	
	  index = (m+axes[1]/2)*axes[0]+n+axes[0]/2;
	
	  pixel_vector = three_vector((m+x_halfpix)*normalization_factor,
				      (n+y_halfpix)*normalization_factor,
				      0,
				      this->circ_ap);
	
	  if(pixel_vector.length()<=1){
	    if(first_arg){
	      val = this->ref_atm_model.caliph_A(*(this->emtr_a), 
						 *(this->emtr_b), 
						 ap_diameter_meters,
						 normalized_tv);
	      val -= this->ref_atm_model.caliph_C(*(this->emtr_a), 
						  *(this->emtr_b), 
						  ap_diameter_meters,
						  normalized_tv,
						  pixel_vector);
	    } else {
	      val = this->ref_atm_model.caliph_B(*(this->emtr_a), 
						 *(this->emtr_b), 
						 ap_diameter_meters,
						 normalized_tv);
	      val -= this->ref_atm_model.caliph_C(*(this->emtr_a), 
						  *(this->emtr_b), 
						  ap_diameter_meters,
						  pixel_vector,
						  normalized_tv);
	    }
	    val -= stored_caliph_D;
	    covariance_pixarr.set_data(index, fac*(covariance_pixarr.data(index) + val));
	  }
	}
      }
      return(covariance_pixarr);
    } catch(...) {
      cerr << "phase_covariance::covariance error - could not form covariance\n";
      throw(string("phase_covariance::covariance"));
    }      
  }

  template<typename precision>
    pixel_array<precision> phase_covariance<precision, circular_aperture>::covariance(double pixel_scale_meters,
										      double wavelength_meters,
										      int nsteps_in_integration,
										      int xindex,
										      int yindex,
										      bool first_arg) const {

    try{
      pixel_array<precision> covariance_pixarr;

      this->initialize_integration(nsteps_in_integration);
      this->initialize_pixel_scale(pixel_scale_meters);

      double caliph_val;
      vector<long> axes = this->stored_caliph_B.get_axes();
      int index = (xindex+axes[1]/2)*axes[0]+yindex+axes[0]/2;
      if(first_arg){
	covariance_pixarr = this->stored_caliph_B;
	caliph_val = this->stored_caliph_A.data(index) - stored_caliph_D;
      } else {
	covariance_pixarr = this->stored_caliph_A;
	caliph_val = this->stored_caliph_B.data(index) - stored_caliph_D; 
      }

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

      double ap_diameter_meters = this->circ_ap.get_diameter();
      double normalization_factor = 2*pixel_scale_meters/ap_diameter_meters;
      double val;
      double fac = Xi*pow(ap_diameter_meters,5/3.)*4*M_PI*M_PI/wavelength_meters/wavelength_meters;

      three_vector pixel_vector_one((xindex+x_halfpix)*normalization_factor,
				    (yindex+y_halfpix)*normalization_factor,
				    0,
				    this->circ_ap);
      if((pixel_vector_one.length()-1)>-three_frame::precision){
	cerr << "phase_covariance::covariance error - "
	     << "three_point lies outside of aperture\n";
	throw(string("phase_covariance::covariance"));
      }
      three_vector pixel_vector_two;

      int half_axes_1 = axes[1]/2;
      int half_axes_0 = axes[0]/2;

      for(int m=-half_axes_1; m<half_axes_1+x_extrapix; m++){
	for(int n=-half_axes_0; n<half_axes_0+y_extrapix; n++){
	
	  index = (m+half_axes_1)*axes[0]+n+half_axes_0;
	
	  pixel_vector_two = three_vector((m+x_halfpix)*normalization_factor,
					  (n+y_halfpix)*normalization_factor,
					  0,
					  this->circ_ap);
	
	  if(pixel_vector_two.length()<=1){
	    val = caliph_val;
	    if(first_arg){
	      val -= this->ref_atm_model.caliph_C(*(this->emtr_a), 
						  *(this->emtr_b), 
						  ap_diameter_meters,
						  pixel_vector_one,
						  pixel_vector_two);
	    } else {
	      val -= this->ref_atm_model.caliph_C(*(this->emtr_a), 
						  *(this->emtr_b), 
						  ap_diameter_meters,
						  pixel_vector_two,
						  pixel_vector_one);
	    }
	    covariance_pixarr.set_data(index, fac*(covariance_pixarr.data(index) + val));
	  }
	}
      }
      return(covariance_pixarr);
    } catch(...) {
      cerr << "phase_covariance::covariance error - could not form covariance\n";
      throw(string("phase_covariance::covariance"));
    }      
  }

  /*
  template<typename precision>
    phase_covariance<precision, annular_aperture>::
    phase_covariance(const emitter & emtr_a,
		     const emitter & emtr_b,
		     const refractive_atmospheric_model & ref_atm_model,
		     const annular_aperture & annular_ap) {


    circular_aperture circ_ap(annular_ap.get_outer_diameter());
    this->outer_aperture_phase_covariance = 
      phase_covariance<precision, circular_aperture>(emtr_a,
						     emtr_b,
						     ref_atm_model,
						     circ_ap);
    
    circ_ap = circular_aperture(annular_ap.get_inner_diameter());
    this->inner_aperture_phase_covariance = 
      phase_covariance<precision, circular_aperture>(emtr_a,
						     emtr_b,
						     ref_atm_model,
						     circ_ap);
    
  }

  template<typename precision>
    phase_covariance<precision, annular_aperture> &
    phase_covariance<precision, annular_aperture>::operator=(const phase_covariance & pcv) {
    this->outer_aperture_phase_covariance = pcv.outer_aperture_phase_covariance;
    this->inner_aperture_phase_covariance = pcv.inner_aperture_phase_covariance;
  }

  template<typename precision>
    void phase_covariance<precision, annular_aperture>::read(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "phase_covariance::read - "
	   << "error opening file " << filename << endl;
      throw(string("phase_covariance::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "phase_covariance::read - "
	   << "error reading "
	   << this->unique_name() << " from file " 
	   << filename << endl;
      throw(string("phase_covariance::read"));
    }
  }
  
  template<typename precision>
    void phase_covariance<precision, annular_aperture>::read(const Arroyo::iofits & iof){

    string type, comment;
    if (!iof.key_exists("TYPE"))
      {
	cerr << this->unique_name() << "::read error - unrecognized "
	     << "file type" << endl;
	throw(this->unique_name() + string("::read"));
      }
    iof.read_key("TYPE",type,comment);
    if(type!=this->unique_name()){
      cerr << this->unique_name() << "::read error - file of type: " 
	   << type << " rather than type " << this->unique_name() << endl;
      throw(this->unique_name() + string("::read"));
    }

    // Move to the next HDU
    iof.movrel_hdu(1);
    
    try{
      this->outer_aperture_phase_covariance.read(iof);
    } catch(...){
      cerr << "phase_covariance::read error - could not read outer aperture phase covariance\n";
      throw(string("phase_covariance::read"));
    }

    try{
      this->inner_aperture_phase_covariance.read(iof);
    } catch(...){
      cerr << "phase_covariance::read error - could not read inner aperture phase covariance\n";
      throw(string("phase_covariance::read"));
    }

  }

  template<typename precision>
    void phase_covariance<precision, annular_aperture>::write(const char * filename) const {

    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "phase_covariance::write - "
	   << "error opening file " << filename << endl;
      throw(string("phase_covariance::write"));
    }

    try{this->write(iof);}
    catch(...){
      cerr << "phase_covariance::write - "
	   << "error writing "
	   << this->unique_name() << " to file " 
	   << filename << endl;
      throw(string("phase_covariance::write"));
    }
  }

  template<typename precision>
    void phase_covariance<precision, annular_aperture>::write(Arroyo::iofits & iof) const {

    fits_header_data<double> tmphdr;
    tmphdr.write(iof);

    string type, comment;

    type = this->unique_name();
    comment = "object type";
    iof.write_key("TYPE", type, comment);

    try{this->outer_aperture_phase_covariance.write(iof);}
    catch(...){
      cerr << "phase_covariance::write error - could not write outer aperture phase covariance\n";
      throw(string("phase_covariance::write"));
    }

    try{this->inner_aperture_phase_covariance.write(iof);}
    catch(...){
      cerr << "phase_covariance::write error - could not write inner aperture phase covariance\n";
      throw(string("phase_covariance::write"));
    }

  }

  template<typename precision>
    void phase_covariance<precision, annular_aperture>::print(std::ostream & os, 
							      const char * prefix) const {
    
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << this->unique_name()
       << "/" << "object type" << endl;
    this->outer_aperture_phase_covariance.print(os, prefix);
    this->inner_aperture_phase_covariance.print(os, prefix);
  }

  template<typename precision>
    double phase_covariance<precision, annular_aperture>::
    aperture_averaged_variance(double wavelength_meters,
			       int nsteps_in_integration) const {

    try{
      return(this->outer_aperture_phase_covariance.aperture_averaged_variance(wavelength_meters,
									      nsteps_in_integration) -
	     this->inner_aperture_phase_covariance.aperture_averaged_variance(wavelength_meters,
									      nsteps_in_integration));
    } catch(...){
      cerr << "phase_covariance::aperture_averaged_variance error - could not form variance\n";
      throw(string("phase_covariance::aperture_averaged_variance"));
    }      

  }

  template<typename precision>
    double phase_covariance<precision, annular_aperture>::variance(double wavelength_meters,
								   int nsteps_in_integration,
								   const three_point & tp) const {

    try{
      return(this->outer_aperture_phase_covariance.variance(wavelength_meters,
							    nsteps_in_integration,
							    tp) -
	     this->inner_aperture_phase_covariance.variance(wavelength_meters,
							    nsteps_in_integration,
							    tp));
    } catch(...) {
      cerr << "phase_covariance::variance error - could not form variance\n";
      throw(string("phase_covariance::variance"));
    }      
  }

  template<typename precision>
  pixel_array<precision> 
  phase_covariance<precision, annular_aperture>::variance(double pixel_scale_meters,
							  double wavelength_meters,
							  int nsteps_in_integration) const {
    
    try{
      return(this->outer_aperture_phase_covariance.variance(pixel_scale_meters,
							    wavelength_meters,
							    nsteps_in_integration) -
	     this->inner_aperture_phase_covariance.variance(pixel_scale_meters,
							    wavelength_meters,
							    nsteps_in_integration));
    } catch(...) {
      cerr << "phase_covariance::variance error - could not form variance\n";
      throw(string("phase_covariance::variance"));
    }      
  }

  template<typename precision>
    double phase_covariance<precision, annular_aperture>::covariance(double wavelength_meters,
								     int nsteps_in_integration,
								     const three_point & tp1,
								     const three_point & tp2) const {
    try {
      return(this->outer_aperture_phase_covariance.covariance(wavelength_meters,
							      nsteps_in_integration,
							      tp1,
							      tp2) -
	     this->inner_aperture_phase_covariance.covariance(wavelength_meters,
							      nsteps_in_integration,
							      tp1,
							      tp2));
    } catch(...) {
      cerr << "phase_covariance::covariance error - could not form covariance\n";
      throw(string("phase_covariance::covariance"));
    }      
  }

  template<typename precision>
    pixel_array<precision> phase_covariance<precision, annular_aperture>::covariance(double pixel_scale_meters,
										      double wavelength_meters,
										      int nsteps_in_integration,
										      const three_point & tp,
										      bool first_arg) const {
    
    try {
      return(this->outer_aperture_phase_covariance.covariance(pixel_scale_meters,
							      wavelength_meters,
							      nsteps_in_integration,
							      tp,
							      first_arg) -
	     this->inner_aperture_phase_covariance.covariance(pixel_scale_meters,
							      wavelength_meters,
							      nsteps_in_integration,
							      tp,
							      first_arg));
    } catch(...) {
      cerr << "phase_covariance::covariance error - could not form covariance\n";
      throw(string("phase_covariance::covariance"));
    }      
  }

  template<typename precision>
    pixel_array<precision> phase_covariance<precision, annular_aperture>::covariance(double pixel_scale_meters,
										      double wavelength_meters,
										      int nsteps_in_integration,
										      int xindex,
										      int yindex,
										      bool first_arg) const {

    try{
      return(this->outer_aperture_phase_covariance.covariance(pixel_scale_meters,
							      wavelength_meters,
							      nsteps_in_integration,
							      xindex,
							      yindex,
							      first_arg) -
	     this->inner_aperture_phase_covariance.covariance(pixel_scale_meters,
							      wavelength_meters,
							      nsteps_in_integration,
							      xindex,
							      yindex,
							      first_arg));
    } catch(...) {
      cerr << "phase_covariance::covariance error - could not form covariance\n";
      throw(string("phase_covariance::covariance"));
    }      
  }
  */

  template<typename precision>
  double tilt_covariance<precision, circular_aperture>::Xi = get_Xi();
  
  template<typename precision>
    void tilt_covariance<precision, circular_aperture>::
    initialize_integration(int nsteps_in_integration) const {

    try{
      if(this->nsteps_in_integration==nsteps_in_integration)
	return;
    
      if(this->emtr_a==NULL || this->emtr_b==NULL){
	cerr << "tilt_covariance::initialize_integration error - uninitialized emitters\n";
	throw(string("tilt_covariance::initialize_integration"));
      }
  
      this->nsteps_in_integration = nsteps_in_integration;
      
      stored_caliph_E = this->ref_atm_model.caliph_E(*(this->emtr_a), 
						     *(this->emtr_b), 
						     this->circ_ap.get_diameter(),
						     this->nsteps_in_integration);
      
      stored_caliph_F = this->ref_atm_model.caliph_F(*(this->emtr_a), 
						     *(this->emtr_b), 
						     this->circ_ap.get_diameter(),
						     this->nsteps_in_integration);      
      
      stored_caliph_F_bar = this->ref_atm_model.caliph_F_bar(*(this->emtr_a), 
							     *(this->emtr_b), 
							     this->circ_ap.get_diameter(),
							     this->nsteps_in_integration);
      
      stored_caliph_K = this->ref_atm_model.caliph_K(*(this->emtr_a), 
						     *(this->emtr_b), 
						     this->circ_ap.get_diameter(),
						     this->nsteps_in_integration);
      
      stored_caliph_L = this->ref_atm_model.caliph_L(*(this->emtr_a), 
						     *(this->emtr_b), 
						     this->circ_ap.get_diameter(),
						     this->nsteps_in_integration);
    } catch(...) {
      cerr << "tilt_covariance::initialize_integration error\n";
      this->emtr_a->print(cerr, "\temtra");
      this->emtr_b->print(cerr, "\temtrb");
      cout << "ap diameter " << this->circ_ap.get_diameter()
	   << "\tnsteps in integration " << this->nsteps_in_integration
	   << endl;
      throw(string("tilt_covariance::initialize_integration"));
    }

  }

  template<typename precision>
    void tilt_covariance<precision, circular_aperture>::
    initialize_pixel_scale(double pixel_scale_meters) const {
    
    try{
      if(this->pixel_scale_meters==pixel_scale_meters)
	return;
      
      if(this->emtr_a==NULL || this->emtr_b==NULL){
	cerr << "tilt_covariance::initialize_pixel_scale error - uninitialized emitters\n";
	throw(string("tilt_covariance::initialize_pixel_scale"));
      }

      this->pixel_scale_meters = pixel_scale_meters;
      stored_caliph_G = 
	this->ref_atm_model.template caliph_G<precision>(*(this->emtr_a),
							 *(this->emtr_b),
							 this->circ_ap.get_diameter(),
							 this->pixel_scale_meters);
      stored_caliph_H = 
	this->ref_atm_model.template caliph_H<precision>(*(this->emtr_a),
							 *(this->emtr_b),
							 this->circ_ap.get_diameter(),
							 this->pixel_scale_meters);
      stored_caliph_I = 
	this->ref_atm_model.template caliph_I<precision>(*(this->emtr_a),
							 *(this->emtr_b),
							 this->circ_ap.get_diameter(),
							 this->pixel_scale_meters);
      stored_caliph_J = 
	this->ref_atm_model.template caliph_J<precision>(*(this->emtr_a),
							 *(this->emtr_b),
							 this->circ_ap.get_diameter(),
							 this->pixel_scale_meters);
    } catch(...) {
      cerr << "tilt_covariance::initialize_pixel_scale error\n";
      this->emtr_a->print(cerr, "\temtra");
      this->emtr_b->print(cerr, "\temtrb");
      cout << "ap diameter " << this->circ_ap.get_diameter()
	   << "\tnsteps in integration " << this->nsteps_in_integration
	   << endl;
      throw(string("tilt_covariance::initialize_pixel_scale"));
    }
  }

  template<typename precision>
    tilt_covariance<precision, circular_aperture>::
    tilt_covariance(const emitter & emtr_a,
		    const emitter & emtr_b,
		    const refractive_atmospheric_model & ref_atm_model,
		    const circular_aperture & circ_ap) {

    this->emtr_a = emitter::emitter_factory(&emtr_a);
    this->emtr_b = emitter::emitter_factory(&emtr_b);

    this->ref_atm_model = ref_atm_model;
    
    this->circ_ap = circ_ap;

    this->stored_caliph_E = -1;
    this->stored_caliph_F = -1;
    this->stored_caliph_F_bar = -1;
    this->stored_caliph_K = -1;
    this->stored_caliph_L = -1;

    this->pixel_scale_meters = -1;
    this->nsteps_in_integration = -1;

    get_omega_hat(*(this->emtr_a), 
		  *(this->emtr_b),
		  this->circ_ap,
		  omega_hat);

  }

  template<typename precision>
    tilt_covariance<precision, circular_aperture> &
    tilt_covariance<precision, circular_aperture>::operator=(const tilt_covariance & tcv) {

    this->emtr_a = emitter::emitter_factory(tcv.emtr_a);
    this->emtr_b = emitter::emitter_factory(tcv.emtr_b);

    this->ref_atm_model = tcv.ref_atm_model;
    this->circ_ap = tcv.circ_ap;
    
    this->pixel_scale_meters = tcv.pixel_scale_meters;
    this->nsteps_in_integration = tcv.nsteps_in_integration;

    this->stored_caliph_G = tcv.stored_caliph_G;
    this->stored_caliph_H = tcv.stored_caliph_H;
    this->stored_caliph_I = tcv.stored_caliph_I;
    this->stored_caliph_J = tcv.stored_caliph_J;

    this->stored_caliph_E = tcv.stored_caliph_E;
    this->stored_caliph_F = tcv.stored_caliph_F;
    this->stored_caliph_F_bar = tcv.stored_caliph_F_bar;
    this->stored_caliph_K = tcv.stored_caliph_K;
    this->stored_caliph_L = tcv.stored_caliph_L;

    this->omega_hat = tcv.omega_hat;
  }

  template<typename precision>
    void tilt_covariance<precision, circular_aperture>::read(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "tilt_covariance::read - "
	   << "error opening file " << filename << endl;
      throw(string("tilt_covariance::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "tilt_covariance::read - "
	   << "error reading "
	   << this->unique_name() << " from file " 
	   << filename << endl;
      throw(string("tilt_covariance::read"));
    }
  }
  
  template<typename precision>
    void tilt_covariance<precision, circular_aperture>::read(const Arroyo::iofits & iof){
    string type, comment;
    if (!iof.key_exists("TYPE"))
      {
	cerr << this->unique_name() << "::read error - unrecognized "
	     << "file type" << endl;
	throw(this->unique_name() + string("::read"));
      }
    iof.read_key("TYPE",type,comment);
    if(type!=this->unique_name()){
      cerr << this->unique_name() << "::read error - file of type: " 
	   << type << " rather than type " << this->unique_name() << endl;
      throw(this->unique_name() + string("::read"));
    }

    // Move to the next HDU
    iof.movrel_hdu(1);

    try{
      emtr_a = emitter::emitter_factory(iof);
    } catch(...){
      cerr << "tilt_covariance::read error - could not read first emitter\n";
      throw(string("tilt_covariance::read"));
    }

    try{
      emtr_b = emitter::emitter_factory(iof);
    } catch(...){
      cerr << "tilt_covariance::read error - could not read second emitter\n";
      throw(string("tilt_covariance::read"));
    }

    try{
      this->circ_ap.read(iof);
    } catch(...){
      cerr << "tilt_covariance::read error - could not read aperture\n";
      throw(string("tilt_covariance::read"));
    }

    try{
      this->ref_atm_model.read(iof);
    } catch(...){
      cerr << "tilt_covariance::read error - could not read refractive atmospheric model\n";
      throw(string("tilt_covariance::read"));
    }

  }

  template<typename precision>
    void tilt_covariance<precision, circular_aperture>::write(const char * filename) const {

    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "tilt_covariance::write - "
	   << "error opening file " << filename << endl;
      throw(string("tilt_covariance::write"));
    }

    try{this->write(iof);}
    catch(...){
      cerr << "tilt_covariance::write - "
	   << "error writing "
	   << this->unique_name() << " to file " 
	   << filename << endl;
      throw(string("tilt_covariance::write"));
    }
  }

  template<typename precision>
    void tilt_covariance<precision, circular_aperture>::write(Arroyo::iofits & iof) const {

    fits_header_data<double> tmphdr;
    tmphdr.write(iof);

    string type, comment;

    type = this->unique_name();
    comment = "object type";
    iof.write_key("TYPE", type, comment);

    try{this->emtr_a->write(iof);}
    catch(...){
      cerr << "tilt_covariance::write error - could not write first emitter\n";
      throw(string("tilt_covariance::write"));
    }
      
    try{this->emtr_b->write(iof);}
    catch(...){
      cerr << "tilt_covariance::write error - could not write second emitter\n";
      throw(string("tilt_covariance::write"));
    }

    try{this->circ_ap.write(iof);}
    catch(...){
      cerr << "tilt_covariance::write error - could not write aperture\n";
      throw(string("tilt_covariance::write"));
    }

    try{this->ref_atm_model.write(iof);}
    catch(...){
      cerr << "tilt_covariance::write error - could not write refractive atmospheric model\n";
      throw(string("tilt_covariance::write"));
    }
  }

  template<typename precision>
    void tilt_covariance<precision, circular_aperture>::print(std::ostream & os, 
					   const char * prefix) const {

    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << this->unique_name()
       << "/" << "object type" << endl;
    this->emtr_a->print(os, prefix);
    this->emtr_b->print(os, prefix);
    this->circ_ap.print(os, prefix);
    this->ref_atm_model.print(os, prefix);
  }

  template<typename precision>
    double tilt_covariance<precision, circular_aperture>::aperture_averaged_variance(double wavelength_meters,
								   int nsteps_in_integration) const {

    try{
    this->initialize_integration(nsteps_in_integration);
    
    return(Xi*pow(this->circ_ap.get_diameter(),5/3.)*
	   4*M_PI*M_PI/wavelength_meters/wavelength_meters*
	   (.25*stored_caliph_F - .5*stored_caliph_E));
    } catch(...){
      cerr << "tilt_covariance::aperture_averaged_variance error - could not form variance\n";
      throw(string("tilt_covariance::aperture_averaged_variance"));
    }      
  }

  template<typename precision>
    double tilt_covariance<precision, circular_aperture>::variance(double wavelength_meters,
								   int nsteps_in_integration,
								   const three_point & tp) const {

    try{
      double val;
      double ap_diameter_meters = this->circ_ap.get_diameter();

      this->initialize_integration(nsteps_in_integration);

      three_vector pixel_vector = tp - static_cast<three_point>(this->circ_ap);
      pixel_vector *= (2/ap_diameter_meters);
      
      if((pixel_vector.length()-1)>-three_frame::precision){
	cerr << "tilt_covariance::variance error - "
	     << "three_point lies outside of aperture\n";
	cerr << "aperture diameter " << ap_diameter_meters 
	     << " three point radius " << pixel_vector.length()*ap_diameter_meters/2. 
	     << " difference " << ap_diameter_meters*(1 - pixel_vector.length())/2.
	     << endl;
	throw(string("tilt_covariance::variance"));
      }
      
      double dot_prod = dot_product(pixel_vector, omega_hat);

      val = pixel_vector.length_squared()*
	(this->stored_caliph_E -
	 this->ref_atm_model.caliph_G(*(this->emtr_a),
				      *(this->emtr_b),
				      ap_diameter_meters,
				      pixel_vector) -
	 this->ref_atm_model.caliph_I(*(this->emtr_a),
				      *(this->emtr_b),
				      ap_diameter_meters,
				      pixel_vector));


      val -= dot_prod*dot_prod*this->stored_caliph_F;

      val += (dot_prod*dot_prod - cross_product(pixel_vector,omega_hat).length_squared())*
	this->stored_caliph_F_bar;

      val += dot_prod*
	(this->ref_atm_model.caliph_H(*(this->emtr_a), 
				      *(this->emtr_b), 
				      ap_diameter_meters,
				      pixel_vector) -
	 this->stored_caliph_K -
	 this->ref_atm_model.caliph_J(*(this->emtr_a), 
				      *(this->emtr_b), 
				      ap_diameter_meters,
				      pixel_vector) +
	 this->stored_caliph_L);

      return(val*Xi*pow(ap_diameter_meters,5/3.)*4*M_PI*M_PI/wavelength_meters/wavelength_meters);
    } catch(...) {
      cerr << "tilt_covariance::variance error - could not form variance\n";
      throw(string("tilt_covariance::variance"));
    }      

  }

  template<typename precision>
    pixel_array<precision> 
    tilt_covariance<precision, circular_aperture>::variance(double pixel_scale_meters,
							    double wavelength_meters,
							    int nsteps_in_integration) const {

    try{
      this->initialize_integration(nsteps_in_integration);
      this->initialize_pixel_scale(pixel_scale_meters);

      double ap_diameter_meters = this->circ_ap.get_diameter();

      vector<long> axes(2,(long)ceil(ap_diameter_meters/pixel_scale_meters));
      pixel_array<precision> variance_pixarr(axes);

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
      double dot_prod, val;
      double normalization_factor = 2*pixel_scale_meters/ap_diameter_meters;
      double fac = 
	Xi*pow(ap_diameter_meters,5/3.)*4*M_PI*M_PI/wavelength_meters/wavelength_meters;
      three_vector pixel_vector, cross_prod;

      for(int m=-axes[1]/2; m<axes[1]/2+x_extrapix; m++){
	for(int n=-axes[0]/2; n<axes[0]/2+y_extrapix; n++){
	
	  index = (m+axes[1]/2)*axes[0]+n+axes[0]/2;
	
	  pixel_vector = three_vector((m+x_halfpix)*normalization_factor,
				      (n+y_halfpix)*normalization_factor,
				      0,
				      this->circ_ap);
	
	  if(pixel_vector.length()<=1){

	    dot_prod = dot_product(pixel_vector, omega_hat);
	    cross_prod = cross_product(pixel_vector, omega_hat);
	  
	    val = pixel_vector.length_squared()*
	      (this->stored_caliph_E -
	       this->stored_caliph_G.data(index) -
	       this->stored_caliph_I.data(index));
	    val -= dot_prod*dot_prod*this->stored_caliph_F;
	    val += (dot_prod*dot_prod - cross_prod.length_squared())*this->stored_caliph_F_bar;
	    val += dot_prod*(stored_caliph_H.data(index) - 
			     stored_caliph_K - 
			     stored_caliph_J.data(index) + 
			     stored_caliph_L);

	    variance_pixarr.set_data(index,fac*val);
	  }
	}
      }
      return(variance_pixarr);
    } catch(...) {
      cerr << "tilt_covariance::variance error - could not form variance\n";
      throw(string("tilt_covariance::variance"));
    }      
  }
    
  template<typename precision>
    double tilt_covariance<precision, circular_aperture>::covariance(double wavelength_meters,
						  int nsteps_in_integration,
						  const three_point & tp1,
						  const three_point & tp2) const {

    try{
      this->initialize_integration(nsteps_in_integration);
      double val;
      double ap_diameter_meters = circ_ap.get_diameter();

      three_vector pixel_vector_1 = tp1 - static_cast<three_point>(this->circ_ap);
      three_vector pixel_vector_2 = tp2 - static_cast<three_point>(this->circ_ap);
      pixel_vector_1 *= (2/ap_diameter_meters);
      pixel_vector_2 *= (2/ap_diameter_meters);

      if((pixel_vector_1.length()-1)>-three_frame::precision){
	cerr << "tilt_covariance::covariance error - "
	     << "first three_point lies outside of aperture\n";
	throw(string("tilt_covariance::covariance"));
      }
      if((pixel_vector_2.length()-1)>-three_frame::precision){
	cerr << "tilt_covariance::covariance error - "
	     << "second three_point lies outside of aperture\n";
	throw(string("tilt_covariance::covariance"));
      }

      double dot_prod_1 = dot_product(pixel_vector_1, omega_hat);
      double dot_prod_2 = dot_product(pixel_vector_2, omega_hat);
      double cross_prod_fac = 
	dot_product(cross_product(pixel_vector_1, omega_hat),
		    cross_product(pixel_vector_2, omega_hat));

      val = 
	-dot_prod_1*dot_prod_2*this->stored_caliph_F +
	(dot_prod_1*dot_prod_2 - cross_prod_fac)*
	this->stored_caliph_F_bar;
      
      val += dot_product(pixel_vector_1, pixel_vector_2)*
      (this->stored_caliph_E -
	 this->ref_atm_model.caliph_G(*(this->emtr_a), 
				      *(this->emtr_b), 
				      ap_diameter_meters,
				      pixel_vector_2) -
	 this->ref_atm_model.caliph_I(*(this->emtr_a), 
				      *(this->emtr_b), 
				      ap_diameter_meters,
				      pixel_vector_1));
       
      val += dot_prod_1*
	(this->ref_atm_model.caliph_H(*(this->emtr_a), 
				      *(this->emtr_b), 
				      ap_diameter_meters,
				      pixel_vector_2) -
	 this->stored_caliph_K);

      val -= dot_prod_2*
	(this->ref_atm_model.caliph_J(*(this->emtr_a), 
				      *(this->emtr_b), 
				      ap_diameter_meters,
				      pixel_vector_1) -
	 this->stored_caliph_L);

      return(val*Xi*pow(ap_diameter_meters,5/3.)*4*M_PI*M_PI/wavelength_meters/wavelength_meters);
    } catch(...) {
      cerr << "tilt_covariance::covariance error - could not form covariance\n";
      throw(string("tilt_covariance::covariance"));
    }      
  }

  template<typename precision>
    pixel_array<precision> 
    tilt_covariance<precision, circular_aperture>::covariance(double pixel_scale_meters,
							      double wavelength_meters,
							      int nsteps_in_integration,
							      const three_point & tp,
							      bool first_arg) const {
    
    try{
      this->initialize_integration(nsteps_in_integration);
      this->initialize_pixel_scale(pixel_scale_meters);

      double ap_diameter_meters = circ_ap.get_diameter();
      vector<long> axes(2,(long)ceil(ap_diameter_meters/pixel_scale_meters));
      pixel_array<precision> covariance_pixarr(axes);

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

      double normalization_factor = 2*pixel_scale_meters/ap_diameter_meters;
      int loop_index;
      double val;
      double fac = Xi*pow(ap_diameter_meters,5/3.)*4*M_PI*M_PI/wavelength_meters/wavelength_meters;
      three_vector pixel_vector;
      three_vector normalized_tv = 
	(tp - static_cast<three_point>(this->circ_ap))*(2/ap_diameter_meters);

      if((normalized_tv.length()-1)>-three_frame::precision){
	cerr << "tilt_covariance::covariance error - "
	     << "three_point lies outside of aperture\n";
	throw(string("tilt_covariance::covariance"));
      }

      double dot_prod_1 = dot_product(normalized_tv, omega_hat);
      double dot_prod_2;
      three_vector cross_prod = cross_product(normalized_tv, omega_hat);
      double cross_prod_fac;

      for(int m=-axes[1]/2; m<axes[1]/2+x_extrapix; m++){
	for(int n=-axes[0]/2; n<axes[0]/2+y_extrapix; n++){
	
	  loop_index = (m+axes[1]/2)*axes[0]+n+axes[0]/2;
	
	  pixel_vector = three_vector((m+x_halfpix)*normalization_factor,
				      (n+y_halfpix)*normalization_factor,
				      0,
				      this->circ_ap);
	
	  if(pixel_vector.length_squared()<=1){

	    cross_prod_fac = 
	      dot_product(cross_product(pixel_vector, omega_hat),
			  cross_prod);
	    dot_prod_2 = dot_product(pixel_vector, omega_hat);

	    val = 
	      -dot_prod_1*dot_prod_2*this->stored_caliph_F +
	      (dot_prod_1*dot_prod_2 - cross_prod_fac)*
	      this->stored_caliph_F_bar;

	    if(first_arg){

	      val += dot_product(pixel_vector, normalized_tv)*
		(this->stored_caliph_E -
		 this->stored_caliph_G.data(loop_index) -
		 this->ref_atm_model.caliph_I(*(this->emtr_a), 
					      *(this->emtr_b), 
					      ap_diameter_meters,
					      normalized_tv));

	      val += dot_prod_1*
		(this->stored_caliph_H.data(loop_index) -
		 this->stored_caliph_K);
	    
	      val -= dot_prod_2*
		(this->ref_atm_model.caliph_J(*(this->emtr_a), 
					      *(this->emtr_b), 
					      ap_diameter_meters,
					      normalized_tv) -
		 this->stored_caliph_L);
	    } else {
	      val += dot_product(pixel_vector, normalized_tv)*
		(this->stored_caliph_E -
		 this->ref_atm_model.caliph_G(*(this->emtr_a), 
					      *(this->emtr_b), 
					      ap_diameter_meters,
					      normalized_tv) -
		 this->stored_caliph_I.data(loop_index));
	    
	      val += dot_prod_2*
		(this->ref_atm_model.caliph_H(*(this->emtr_a), 
					      *(this->emtr_b), 
					      ap_diameter_meters,
					      normalized_tv) -
		 this->stored_caliph_K);

	      val -= dot_prod_1*
		(this->stored_caliph_J.data(loop_index) -
		 this->stored_caliph_L);
	    }

	    covariance_pixarr.set_data(loop_index, fac*val);
	  }
	}
      }
      return(covariance_pixarr);
    } catch(...) {
      cerr << "tilt_covariance::covariance error - could not form covariance\n";
      throw(string("tilt_covariance::covariance"));
    }      
  }

  template<typename precision>
    pixel_array<precision> 
    tilt_covariance<precision, circular_aperture>::covariance(double pixel_scale_meters,
							      double wavelength_meters,
							      int nsteps_in_integration,
							      int xindex,
							      int yindex,
							      bool first_arg) const {

    try{
      this->initialize_integration(nsteps_in_integration);
      this->initialize_pixel_scale(pixel_scale_meters);

      double ap_diameter_meters = circ_ap.get_diameter();
      vector<long> axes(2,(long)ceil(ap_diameter_meters/pixel_scale_meters));
      pixel_array<precision> covariance_pixarr(axes);

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

      double normalization_factor = 2*pixel_scale_meters/ap_diameter_meters;
      int arg_index = (xindex+axes[1]/2)*axes[0]+yindex+axes[0]/2;
      int loop_index;
      double val;
      double fac = Xi*pow(ap_diameter_meters,5/3.)*4*M_PI*M_PI/wavelength_meters/wavelength_meters;
      three_vector pixel_vector;
      three_vector normalized_tv = (three_point((xindex+x_halfpix)*normalization_factor,
						(yindex+y_halfpix)*normalization_factor,
						0,
						this->circ_ap) - 
				    this->circ_ap);

      if((normalized_tv.length()-1)>-three_frame::precision){
	cerr << "tilt_covariance::covariance error - "
	     << "three_point lies outside of aperture\n";
	throw(string("tilt_covariance::covariance"));
      }

      double dot_prod_1 = dot_product(normalized_tv, omega_hat);
      double dot_prod_2;
      three_vector cross_prod = cross_product(normalized_tv, omega_hat);
      double cross_prod_fac;

      int half_axes_1 = axes[1]/2;
      int half_axes_0 = axes[0]/2;

      double arg_index_caliph_G = this->stored_caliph_G.data(arg_index);
      double arg_index_caliph_H = this->stored_caliph_H.data(arg_index);
      double arg_index_caliph_I = this->stored_caliph_I.data(arg_index);
      double arg_index_caliph_J = this->stored_caliph_J.data(arg_index);


      for(int m=-half_axes_1; m<half_axes_1+x_extrapix; m++){
	for(int n=-half_axes_0; n<half_axes_0+y_extrapix; n++){
	
	  loop_index = (m+half_axes_1)*axes[0]+n+half_axes_0;
	
	  pixel_vector = three_vector((m+x_halfpix)*normalization_factor,
				      (n+y_halfpix)*normalization_factor,
				      0,
				      this->circ_ap);
	
	  if(pixel_vector.length_squared()<=1){

	    cross_prod_fac = 
	      dot_product(cross_product(pixel_vector, omega_hat),
			  cross_prod);
	    dot_prod_2 = dot_product(pixel_vector, omega_hat);

	    val = 
	      -dot_prod_1*dot_prod_2*this->stored_caliph_F +
	      (dot_prod_1*dot_prod_2 - cross_prod_fac)*
	      this->stored_caliph_F_bar;

	    if(first_arg){

	      val += dot_product(pixel_vector, normalized_tv)*
		(this->stored_caliph_E -
		 this->stored_caliph_G.data(loop_index) -
		 arg_index_caliph_I);

	      val += dot_prod_1*
		(this->stored_caliph_H.data(loop_index) -
		 this->stored_caliph_K);

	      val -= dot_prod_2*
		(arg_index_caliph_J -
		 this->stored_caliph_L);
	    } else {
	      val += dot_product(pixel_vector, normalized_tv)*
		(this->stored_caliph_E -
		 arg_index_caliph_G - 
		 this->stored_caliph_I.data(loop_index));
	    
	      val += dot_prod_2*
		(arg_index_caliph_H - 
		 this->stored_caliph_K);

	      val -= dot_prod_1*
		(this->stored_caliph_J.data(loop_index) -
		 this->stored_caliph_L);
	    }

	    covariance_pixarr.set_data(loop_index, fac*val);
	  }
	}
      }
      return(covariance_pixarr);
    } catch(...) {
      cerr << "tilt_covariance::covariance error - could not form covariance\n";
      throw(string("tilt_covariance::covariance"));
    }      
  }

  /*
  template<typename precision>
    tilt_covariance<precision, annular_aperture>::
    tilt_covariance(const emitter & emtr_a,
		    const emitter & emtr_b,
		    const refractive_atmospheric_model & ref_atm_model,
		    const annular_aperture & annular_ap) {

    circular_aperture circ_ap(annular_ap.get_outer_diameter());
    this->outer_aperture_tilt_covariance = 
      tilt_covariance<precision, circular_aperture>(emtr_a,
						    emtr_b,
						    ref_atm_model,
						    circ_ap);
    
    circ_ap = circular_aperture(annular_ap.get_inner_diameter());
    this->inner_aperture_tilt_covariance = 
      tilt_covariance<precision, circular_aperture>(emtr_a,
						    emtr_b,
						    ref_atm_model,
						    circ_ap);

  }

  template<typename precision>
    tilt_covariance<precision, annular_aperture> &
    tilt_covariance<precision, annular_aperture>::operator=(const tilt_covariance & tcv) {

    this->outer_aperture_tilt_covariance = tcv.outer_aperture_tilt_covariance;
    this->inner_aperture_tilt_covariance = tcv.inner_aperture_tilt_covariance;
  }

  template<typename precision>
    void tilt_covariance<precision, annular_aperture>::read(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "tilt_covariance::read - "
	   << "error opening file " << filename << endl;
      throw(string("tilt_covariance::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "tilt_covariance::read - "
	   << "error reading "
	   << this->unique_name() << " from file " 
	   << filename << endl;
      throw(string("tilt_covariance::read"));
    }
  }
  
  template<typename precision>
    void tilt_covariance<precision, annular_aperture>::read(const Arroyo::iofits & iof){
    string type, comment;
    if (!iof.key_exists("TYPE"))
      {
	cerr << this->unique_name() << "::read error - unrecognized "
	     << "file type" << endl;
	throw(this->unique_name() + string("::read"));
      }
    iof.read_key("TYPE",type,comment);
    if(type!=this->unique_name()){
      cerr << this->unique_name() << "::read error - file of type: " 
	   << type << " rather than type " << this->unique_name() << endl;
      throw(this->unique_name() + string("::read"));
    }

    // Move to the next HDU
    iof.movrel_hdu(1);


    try{
      this->outer_aperture_tilt_covariance.read(iof);
    } catch(...){
      cerr << "tilt_covariance::read error - could not read outer aperture tilt covariance\n";
      throw(string("tilt_covariance::read"));
    }

    try{
      this->inner_aperture_tilt_covariance.read(iof);
    } catch(...){
      cerr << "tilt_covariance::read error - could not read inner aperture tilt covariance\n";
      throw(string("tilt_covariance::read"));
    }

  }

  template<typename precision>
    void tilt_covariance<precision, annular_aperture>::write(const char * filename) const {

    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "tilt_covariance::write - "
	   << "error opening file " << filename << endl;
      throw(string("tilt_covariance::write"));
    }

    try{this->write(iof);}
    catch(...){
      cerr << "tilt_covariance::write - "
	   << "error writing "
	   << this->unique_name() << " to file " 
	   << filename << endl;
      throw(string("tilt_covariance::write"));
    }
  }

  template<typename precision>
    void tilt_covariance<precision, annular_aperture>::write(Arroyo::iofits & iof) const {

    fits_header_data<double> tmphdr;
    tmphdr.write(iof);

    string type, comment;

    type = this->unique_name();
    comment = "object type";
    iof.write_key("TYPE", type, comment);

    try{this->outer_aperture_tilt_covariance.write(iof);}
    catch(...){
      cerr << "tilt_covariance::write error - could not write outer aperture tilt covariance\n";
      throw(string("tilt_covariance::write"));
    }
      
    try{this->inner_aperture_tilt_covariance.write(iof);}
    catch(...){
      cerr << "tilt_covariance::write error - could not write inner aperture tilt covariance\n";
      throw(string("tilt_covariance::write"));
    }
      
  }

  template<typename precision>
    void tilt_covariance<precision, annular_aperture>::print(std::ostream & os, 
							     const char * prefix) const {

    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << this->unique_name()
       << "/" << "object type" << endl;
    this->outer_aperture_tilt_covariance.print(os, prefix);
    this->inner_aperture_tilt_covariance.print(os, prefix);
  }

  template<typename precision>
    double tilt_covariance<precision, annular_aperture>::aperture_averaged_variance(double wavelength_meters,
										    int nsteps_in_integration) const {
    
    try{
      return(this->outer_aperture_tilt_covariance.aperture_averaged_variance(wavelength_meters,
									     nsteps_in_integration) - 
	     this->inner_aperture_tilt_covariance.aperture_averaged_variance(wavelength_meters,
									     nsteps_in_integration));
    } catch(...){
      cerr << "tilt_covariance::aperture_averaged_variance error - could not form variance\n";
      throw(string("tilt_covariance::aperture_averaged_variance"));
    }      
  }

  template<typename precision>
    double tilt_covariance<precision, annular_aperture>::variance(double wavelength_meters,
								  int nsteps_in_integration,
								  const three_point & tp) const {
    
    try{
      return(this->outer_aperture_tilt_covariance.variance(wavelength_meters,
							   nsteps_in_integration,
							   tp) - 
	     this->inner_aperture_tilt_covariance.variance(wavelength_meters,
							   nsteps_in_integration,
							   tp));
    } catch(...) {
      cerr << "tilt_covariance::variance error - could not form variance\n";
      throw(string("tilt_covariance::variance"));
    }      

  }

  template<typename precision>
    pixel_array<precision> 
    tilt_covariance<precision, annular_aperture>::variance(double pixel_scale_meters,
							   double wavelength_meters,
							   int nsteps_in_integration) const {

    try{
      return(this->outer_aperture_tilt_covariance.variance(pixel_scale_meters,
							   wavelength_meters,
							   nsteps_in_integration) - 
	     this->inner_aperture_tilt_covariance.variance(pixel_scale_meters,
							   wavelength_meters,
							   nsteps_in_integration));
    } catch(...) {
      cerr << "tilt_covariance::variance error - could not form variance\n";
      throw(string("tilt_covariance::variance"));
    }      
  }
    
  template<typename precision>
    double tilt_covariance<precision, annular_aperture>::covariance(double wavelength_meters,
								    int nsteps_in_integration,
								    const three_point & tp1,
								    const three_point & tp2) const {
    
    try{
      return(this->outer_aperture_tilt_covariance.covariance(wavelength_meters,
							     nsteps_in_integration,
							     tp1,
							     tp2) - 
	     this->inner_aperture_tilt_covariance.covariance(wavelength_meters,
							     nsteps_in_integration,
							     tp1,
							     tp2));
    } catch(...) {
      cerr << "tilt_covariance::covariance error - could not form covariance\n";
      throw(string("tilt_covariance::covariance"));
    }      
  }

  template<typename precision>
    pixel_array<precision> 
    tilt_covariance<precision, annular_aperture>::covariance(double pixel_scale_meters,
							     double wavelength_meters,
							     int nsteps_in_integration,
							     const three_point & tp,
							     bool first_arg) const {
    
    try{
      return(this->outer_aperture_tilt_covariance.covariance(pixel_scale_meters,
							     wavelength_meters,
							     nsteps_in_integration,
							     tp,
							     first_arg) - 
	     this->inner_aperture_tilt_covariance.covariance(pixel_scale_meters,
							     wavelength_meters,
							     nsteps_in_integration,
							     tp,
							     first_arg));
    } catch(...) {
      cerr << "tilt_covariance::covariance error - could not form covariance\n";
      throw(string("tilt_covariance::covariance"));
    }      
  }

  template<typename precision>
    pixel_array<precision> 
    tilt_covariance<precision, annular_aperture>::covariance(double pixel_scale_meters,
							     double wavelength_meters,
							     int nsteps_in_integration,
							     int xindex,
							     int yindex,
							     bool first_arg) const {

    try{
      return(this->outer_aperture_tilt_covariance.covariance(pixel_scale_meters,
							     wavelength_meters,
							     nsteps_in_integration,
							     xindex,
							     yindex,
							     first_arg) - 
	     this->inner_aperture_tilt_covariance.covariance(pixel_scale_meters,
							     wavelength_meters,
							     nsteps_in_integration,
							     xindex,
							     yindex,
							     first_arg));
    } catch(...) {
      cerr << "tilt_covariance::covariance error - could not form covariance\n";
      throw(string("tilt_covariance::covariance"));
    }      
  }
  */
}

#endif
