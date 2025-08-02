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

#ifndef REFRACTIVE_ATMOSPHERE_H
#define REFRACTIVE_ATMOSPHERE_H

#include "three_frame.h"
#include "region_base.h"
#include "computational_geometry.h"
#include "AO_sim_base.h"
#include "fits_header_data.h"
#include "pixel_array.h"
#include "power_spectrum.h"
#include "propagation_plan.h"
#include "aperture.h"
#include "emitter.h"
#include "special_functions.h"

namespace Arroyo {

  using std::vector;
  using std::ostream;

  class subharmonic_method;

  /// A class to represent a refractive atmospheric model
  /// as an array of power spectra at different heights.
  ///  This class serves as a base class for a hierarchy of
  ///  atmospheric models.  This class is expected to contain most
  ///  or all of the data members, and the derived classes serve as
  ///  a useful way for the user to instantiate the different models.

  class refractive_atmospheric_model :
    virtual public AO_sim_base {

    private:
    
    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("refractive atmospheric model"));};

    protected:

    /// The array of power spectra - one for each layer
    vector<power_spectrum *> power_spectra_;

    /// The layer heights in meters
    vector<double> layer_heights_;

    /// The reference frame for the ground
    three_frame ground_ref_frame_;

    ///////////////////////////////////////////
    /// Read the data members from file.  This member
    /// function may be called directly from the
    /// read member function of a derived class,
    /// after the unique type key has been checked.
    void read_common_data(const iofits & iof);

    ///////////////////////////////////////////
    /// Write the data members to file.  This member
    /// function may be called directly from the
    /// write member function of a derived class,
    /// after the unique type key has been written.
    void write_common_data(iofits & iof) const;

    public:
  
    ///////////////////////////////////////////
    ///  Null constructor
    refractive_atmospheric_model(){};

    ///////////////////////////////////////////
    ///  Copy constructor
    refractive_atmospheric_model(const refractive_atmospheric_model & ref_atm_model);

    ///////////////////////////////////////////
    ///  Construct from a file
    refractive_atmospheric_model(const char * filename);

    ///////////////////////////////////////////
    ///  Construct from an iofits object
    refractive_atmospheric_model(const iofits & iof);

    ///////////////////////////////////////////
    ///  Construct from the bits
    ///
    ///  Layer heights are measured in meters.  The ground reference
    ///  frame serves to define the zero points for the heights, which
    ///  are measured along the positive z axis of this frame
    ///
    refractive_atmospheric_model(const vector<power_spectrum *> & power_spectra, 
				 const vector<double> & layer_heights,
				 const three_frame & ground_ref_frame);

    ///////////////////////////////////////////
    ///  Construct a model with Komolgorov power spectrum from the bits
    ///
    ///  Layer heights are measured in meters.  
    ///
    ///  Cn2 dz values are measured in meters^{2/3}
    ///
    ///  The ground reference frame serves to define the zero points
    ///  for the heights, which are measured along the positive z axis
    ///  of this frame
    ///
    refractive_atmospheric_model(const vector<double> & layer_heights,
				 const vector<double> & Cn2_dz,
				 const three_frame & ground_ref_frame);

    ///////////////////////////////////////////
    ///  Destructor
    ~refractive_atmospheric_model();

    ///////////////////////////////////////////
    ///  Operator = 
    refractive_atmospheric_model & 
      operator=(const refractive_atmospheric_model & ref_atm_model);

    ///////////////////////////////////////////
    /// Virtual clone method
    ///
    /// Calling routine is responsible for deleting memory
    virtual refractive_atmospheric_model * clone() const {
      return new refractive_atmospheric_model(*this);
    };

    ///////////////////////////////////////////
    ///  Read from file
    void read(const char * filename);

    ///////////////////////////////////////////
    ///  Read from iofits
    void read(const iofits & iof);

    ///////////////////////////////////////////
    ///  Write to file
    void write(const char * filename) const;

    ///////////////////////////////////////////
    ///  Write to iofits
    void write(iofits & iof) const;

    ///////////////////////////////////////////
    ///  Print
    void print(ostream & os, const char * prefix="") const;

    ///////////////////////////////////////////
    ///  Get the number of layers in the model
    long get_number_of_layers() const {return layer_heights_.size();};

    ///////////////////////////////////////////
    ///  Get the number of layers in the model
    vector<double> get_layer_heights() const {return layer_heights_;};

    ///////////////////////////////////////////
    ///  Get the power spectra in the model
    ///  Memory is dynamically allocated, and must be freed.
    vector<power_spectrum *> get_power_spectra() const {
      vector<power_spectrum *> pspec(this->power_spectra_.size());
      for(uint i=0; i<pspec.size(); i++)
	pspec[i] = power_spectrum::power_spectrum_factory(this->power_spectra_[i]);
      return(pspec);
    };

    ///////////////////////////////////////////
    ///  Get the three frame.
    three_frame get_three_frame() const {return this->ground_ref_frame_;};

    ///////////////////////////////////////////
    ///  Get a turbulence moment
    ///
    ///  The value returned has units of 
    ///  meters^{1/3 + moment}
    ///
    ///  This function throws an error unless 
    ///  all layers in the model have a Komolgorov
    ///  turbulence power spectrum.
    double turbulence_moment(double moment, 
			     double zenith_angle_degrees=0) const;

    ///////////////////////////////////////////
    ///  Get a velocity moment
    ///
    ///  The wind velocities should be in meters/sec,
    ///  and are assumed to be ordered in the same
    ///  way as the heights and power spectra
    ///
    ///  The value returned has units of 
    ///  meters^{1/3 + moment} sec^{-moment}
    ///
    ///  This function throws an error unless 
    ///  all layers in the model have a Komolgorov
    ///  turbulence power spectrum.
    double velocity_moment(const vector<three_vector> & layer_wind_velocities_meters_per_sec,
			   double moment, 
			   double azimuth_angle_degrees=0,
			   double zenith_angle_degrees=0) const;

    ///////////////////////////////////////////
    ///  Get the Fried parameter in an upward-looking
    ///  scenario.
    ///
    ///  The value returned is in meters
    ///
    ///  This function throws an error unless 
    ///  all layers in the model have a Komolgorov
    ///  turbulence power spectrum.
    double fried_parameter(double wavelength_meters, 
			   double zenith_angle_degrees=0, 
			   double guide_star_height_meters=-1) const;


    ///////////////////////////////////////////
    ///  Get the Fried parameter in a downward-looking
    ///  scenario.
    ///
    ///  The value returned is in meters
    ///
    ///  This function throws an error unless 
    ///  all layers in the model have a Komolgorov
    ///  turbulence power spectrum.
    double fried_parameter_downward(double wavelength_meters, 
				    double nadir_angle_degrees,
				    double altitude_meters) const;
    
    ///////////////////////////////////////////
    ///  Get the two axis tilt jitter in an upward-looking
    ///  scenario.
    ///
    ///  The value returned is in radians
    ///
    ///  This function throws an error unless 
    ///  all layers in the model have a Komolgorov
    ///  turbulence power spectrum.
    double two_axis_tilt_jitter(double wavelength_meters, 
				double aperture_diameter_meters,
				double zenith_angle_degrees=0, 
				double guide_star_height_meters=-1) const;


    ///////////////////////////////////////////
    ///  Get the two axis tilt jitter in a downward-looking
    ///  scenario.
    ///
    ///  The value returned is in radians
    ///
    ///  This function throws an error unless 
    ///  all layers in the model have a Komolgorov
    ///  turbulence power spectrum.
    double two_axis_tilt_jitter_downward(double wavelength_meters, 
					 double aperture_diameter_meters,
					 double nadir_angle_degrees,
					 double altitude_meters) const;
    
    ///////////////////////////////////////////
    ///  Get the seeing
    ///
    ///  The value returned is in radians
    ///
    ///  This function throws an error unless 
    ///  all layers in the model have a Komolgorov
    ///  turbulence power spectrum.
    double seeing(double wavelength_meters, 
		  double zenith_angle_degrees=0) const;

    ///////////////////////////////////////////
    ///  Get the isoplanatic angle in an upward-looking
    ///  scenario
    ///
    ///  The value returned is in radians
    ///
    ///  This function throws an error unless 
    ///  all layers in the model have a Komolgorov
    ///  turbulence power spectrum.
    double isoplanatic_angle(double wavelength_meters, 
			     double zenith_angle_degrees=0,
			     double guide_star_height_meters=-1) const;

    ///////////////////////////////////////////
    ///  Get the isoplanatic angle in an downward-looking
    ///  scenario
    ///
    ///  The value returned is in radians
    ///
    ///  This function throws an error unless 
    ///  all layers in the model have a Komolgorov
    ///  turbulence power spectrum.
    double isoplanatic_angle_downward(double wavelength_meters, 
				      double nadir_angle_degrees,
				      double altitude_meters) const;

    ///////////////////////////////////////////
    ///  Get the isokinetic angle in an upward-looking
    ///  scenario
    ///
    ///  The value returned is in radians
    ///
    ///  This function throws an error unless 
    ///  all layers in the model have a Komolgorov
    ///  turbulence power spectrum.
    double isokinetic_angle(double wavelength_meters,
			    double aperture_diameter_meters,
			    double zenith_angle_degrees=0,
			    double guide_star_height_meters=-1) const;

    ///////////////////////////////////////////
    ///  Get the isokinetic angle in a downward-looking
    ///  scenario
    ///
    ///  The value returned is in radians
    ///
    ///  This function throws an error unless 
    ///  all layers in the model have a Komolgorov
    ///  turbulence power spectrum.
    double isokinetic_angle_downward(double wavelength_meters,
				     double aperture_diameter_meters,
				     double nadir_angle_degrees,
				     double altitude_meters) const;
    
    ///////////////////////////////////////////
    ///  Get the Greenwood frequency
    ///
    ///  The wind velocities should be in meters/sec,
    ///  and ordered in the same way as the layer
    ///  heights
    ///
    ///  The value returned is in Hertz
    ///
    ///  This function throws an error unless 
    ///  all layers in the model have a Komolgorov
    ///  turbulence power spectrum.
    double greenwood_frequency(vector<three_vector> & layer_wind_velocities_meters_per_sec,
			       double wavelength_meters, 
			       double azimuth_angle_degrees=0,
			       double zenith_angle_degrees=0) const;

    ///////////////////////////////////////////
    ///  Get the value of d_0 (See Hardy eq 7.36)
    ///
    ///  The value returned is in meters
    ///
    ///  This function throws an error unless 
    ///  all layers in the model have a Komolgorov
    ///  turbulence power spectrum.
    double d_0(double guide_star_height_meters,
               double wavelength_meters,
	       double zenith_angle_degrees=0) const;


    ///////////////////////////////////////////
    ///  Get the value of the log amplitude variance
    ///  in an upward-looking scenario.
    ///
    ///  The value returned is dimensionless
    ///
    ///  This function throws an error unless 
    ///  all layers in the model have a Komolgorov
    ///  turbulence power spectrum.
    double log_amplitude_variance(double wavelength_meters,
				  double zenith_angle_degrees,
				  double altitude_meters=-1) const;

    ///////////////////////////////////////////
    ///  Get the value of the log amplitude variance
    ///  in a downward-looking scenario.
    ///
    ///  The value returned is dimensionless
    ///
    ///  This function throws an error unless 
    ///  all layers in the model have a Komolgorov
    ///  turbulence power spectrum.
    double log_amplitude_variance_downward(double wavelength_meters,
					   double nadir_angle_degrees,
					   double altitude_meters) const;

    ///////////////////////////////////////////
    ///  Get the value of the scintillation quenching diameter
    ///
    ///  The value returned is in meters
    ///
    ///  This function throws an error unless 
    ///  all layers in the model have a Komolgorov
    ///  turbulence power spectrum.
    double scintillation_quenching_diameter(double wavelength_meters,
					    double zenith_angle_degrees=0) const;

    ///////////////////////////////////////////
    ///  Get the value of the scintillation isoplanatic angle
    ///
    ///  The value returned is in radians
    ///
    ///  This function throws an error unless 
    ///  all layers in the model have a Komolgorov
    ///  turbulence power spectrum.
    double scintillation_isoplanatic_angle(double wavelength_meters,
					   double zenith_angle_degrees=0) const;

    ///////////////////////////////////////////
    ///  Tyler F_1
    double Tyler_F_1(double x) const;

    ///////////////////////////////////////////
    ///  Tyler 
    double Tyler_H(double & rho, 
	     double & omega,
	     vector<double> & numerator_args,
	     vector<double> & denominator_args) const;

    ///////////////////////////////////////////
    ///  Tyler 
    double Tyler_K_1(double rho, double q) const;

    ///////////////////////////////////////////
    ///  Tyler 
    double Tyler_F_2(double q, 
	       double omega, 
	       int nsamples_in_integration) const;
      
    ///////////////////////////////////////////
    ///  Tyler 
    double Tyler_G_hat(double rho) const;

    ///////////////////////////////////////////
    ///  Tyler 
    double Tyler_G_1(three_vector r1,
	       three_vector r2) const;

    ///////////////////////////////////////////
    ///  Tyler 
    double Tyler_G_2(three_vector r, 
	       double q, 
	       three_vector omega, 
	       int nsamples_in_integration) const;

    ///////////////////////////////////////////
    ///  Tyler 
    double Tyler_G_3(double q, 
	       double omega,
	       int nsamples_in_integration) const;
    
    ///////////////////////////////////////////
    ///  Tyler 
    double Tyler_G_4(three_vector r1, 
	       three_vector r2, 
	       double q, 
	       three_vector omega, 
	       int nsamples_in_integration) const;
    
    ///////////////////////////////////////////
    ///  Tyler 
    void Tyler_get_constants(const emitter & emtr_a,
		       const emitter & emtr_b,
		       const three_frame & tf,
		       double aperture_diameter_meters,
		       double & secant_zenith_angle,
		       double & max_range_meters,
		       double & min_range_meters,
		       three_vector & little_omega) const;
    
    ///////////////////////////////////////////
    ///  Get the vector of cn2 coefficients
    ///
    ///  These are in units of meters^{1/3}
    vector<double> get_cn2_coefficients() const;

    ///////////////////////////////////////////
    ///  A_{ab} term in phase covariance
    /// 
    ///  The vector rho must have amplitude less than unity
    double caliph_A(const emitter & emtr_a,
		    const emitter & emtr_b,
		    double aperture_diameter_meters,
		    const three_vector & rho) const;

    ///////////////////////////////////////////
    ///  A_{ab} term in phase covariance
    template<class T>
    pixel_array<T> caliph_A(const emitter & emtr_a,
			    const emitter & emtr_b,
			    double aperture_diameter_meters,
			    double pixel_scale_meters) const;

    ///////////////////////////////////////////
    ///  B_{ab} term in phase covariance
    /// 
    ///  The vector rho must have amplitude less than unity
    double caliph_B(const emitter & emtr_a,
		    const emitter & emtr_b,
		    double aperture_diameter_meters,
		    const three_vector & rho) const;

    ///////////////////////////////////////////
    ///  B_{ab} term in phase covariance
    template<class T>
    pixel_array<T> caliph_B(const emitter & emtr_a,
			    const emitter & emtr_b,
			    double aperture_diameter_meters,
			    double pixel_scale_meters) const;

    ///////////////////////////////////////////
    ///  C_{ab} term in phase covariance
    /// 
    ///  The vectors rho_1 and rho_2 must have amplitudes less than
    ///  unity
    double caliph_C(const emitter & emtr_a,
		    const emitter & emtr_b,
		    double aperture_diameter_meters,
		    const three_vector & rho_1,
		    const three_vector & rho_2) const;

    ///////////////////////////////////////////
    ///  C_{ab} term in phase covariance
    template<class T>
    pixel_array<T> caliph_C(const emitter & emtr_a,
			    const emitter & emtr_b,
			    double aperture_diameter_meters,
			    double pixel_scale_meters) const;

    ///////////////////////////////////////////
    ///  D_{ab} term in phase covariance
    double caliph_D(const emitter & emtr_a,
		    const emitter & emtr_b,
		    double aperture_diameter_meters,
		    int nsteps_in_integration) const;

    ///////////////////////////////////////////
    ///  E_{ab} term in phase covariance
    double caliph_E(const emitter & emtr_a,
		    const emitter & emtr_b,
		    double aperture_diameter_meters,
		    int nsteps_in_integration) const;

    ///////////////////////////////////////////
    ///  F_{ab} term in phase covariance
    double caliph_F(const emitter & emtr_a,
		    const emitter & emtr_b,
		    double aperture_diameter_meters,
		    int nsteps_in_integration) const;

    ///////////////////////////////////////////
    ///  \bar{F}_{ab} term in phase covariance
    double caliph_F_bar(const emitter & emtr_a,
			const emitter & emtr_b,
			double aperture_diameter_meters,
			int nsteps_in_integration) const;

    /*
    ///////////////////////////////////////////
    ///  Test of Tyler's G_1
    double tmp_tyler_G_1_plus(const emitter & emtr_a,
			      const emitter & emtr_b,
			      double aperture_diameter_meters,
			      const three_vector & rho) const;

    ///////////////////////////////////////////
    ///  Test of Tyler's G_1
    double tmp_tyler_G_1_minus(const emitter & emtr_a,
			       const emitter & emtr_b,
			       double aperture_diameter_meters,
			       const three_vector & rho) const;

    ///////////////////////////////////////////
    ///  Test of Tyler's G_3
    double caliph_G_3(const emitter & emtr_a,
		      const emitter & emtr_b,
		      double aperture_diameter_meters,
		      int nsteps_in_integration) const;
    */

    ///////////////////////////////////////////
    ///  G_{ab} term in phase covariance
    /// 
    ///  The vector rho must have amplitude less than unity
    double caliph_G(const emitter & emtr_a,
		    const emitter & emtr_b,
		    double aperture_diameter_meters,
		    const three_vector & rho) const;

    ///////////////////////////////////////////
    ///  G_{ab} term in phase covariance
    template<class T>
    pixel_array<T> caliph_G(const emitter & emtr_a,
			    const emitter & emtr_b,
			    double aperture_diameter_meters,
			    double pixel_scale_meters) const;

    ///////////////////////////////////////////
    ///  H_{ab} term in phase covariance
    /// 
    ///  The vector rho must have amplitude less than unity
    double caliph_H(const emitter & emtr_a,
		    const emitter & emtr_b,
		    double aperture_diameter_meters,
		    const three_vector & rho) const;

    ///////////////////////////////////////////
    ///  H_{ab} term in phase covariance
    template<class T>
    pixel_array<T> caliph_H(const emitter & emtr_a,
			    const emitter & emtr_b,
			    double aperture_diameter_meters,
			    double pixel_scale_meters) const;

    ///////////////////////////////////////////
    ///  I_{ab} term in phase covariance
    /// 
    ///  The vector rho must have amplitude less than unity
    double caliph_I(const emitter & emtr_a,
		    const emitter & emtr_b,
		    double aperture_diameter_meters,
		    const three_vector & rho) const;

    ///////////////////////////////////////////
    ///  I_{ab} term in phase covariance
    template<class T>
    pixel_array<T> caliph_I(const emitter & emtr_a,
			    const emitter & emtr_b,
			    double aperture_diameter_meters,
			    double pixel_scale_meters) const;

    ///////////////////////////////////////////
    ///  J_{ab} term in phase covariance
    /// 
    ///  The vector rho must have amplitude less than unity
    double caliph_J(const emitter & emtr_a,
		    const emitter & emtr_b,
		    double aperture_diameter_meters,
		    const three_vector & rho) const;

    ///////////////////////////////////////////
    ///  J_{ab} term in phase covariance
    template<class T>
    pixel_array<T> caliph_J(const emitter & emtr_a,
			    const emitter & emtr_b,
			    double aperture_diameter_meters,
			    double pixel_scale_meters) const;
    
    ///////////////////////////////////////////
    ///  K_{ab} term in phase covariance
    double caliph_K(const emitter & emtr_a,
		    const emitter & emtr_b,
		    double aperture_diameter_meters,
		    int nsteps_in_integration) const;

    ///////////////////////////////////////////
    ///  L_{ab} term in phase covariance
    double caliph_L(const emitter & emtr_a,
		    const emitter & emtr_b,
		    double aperture_diameter_meters,
		    int nsteps_in_integration) const;

    ///////////////////////////////////////////
    ///  M_{ab} term in phase covariance
    double caliph_M(const emitter & emtr_a,
		    const emitter & emtr_b,
		    double aperture_diameter_meters) const;

    ///////////////////////////////////////////
    ///  Get the piston removed phase covariance between
    ///  two beams, averaged over a circular aperture.  
    ///
    ///  This quantity is computed from Tyler, 
    ///  JOSA 11 p 339 1994
    ///  
    ///  The value returned is in radians squared
    ///
    ///  This function throws an error unless 
    ///  all layers in the model have a Komolgorov
    ///  turbulence power spectrum.
    double aperture_averaged_phase_covariance(const emitter & emtr_a,
					      const emitter & emtr_b,
					      const aperture & ap,
					      double wavelength_meters,
					      int nsteps_in_integration=1000) const;

    ///////////////////////////////////////////
    ///  Get the tilt covariance between two beams, 
    ///  averaged over a circular aperture.
    ///
    ///  This quantity is computed from Tyler, 
    ///  JOSA 11 p 339 1994
    ///  
    ///  The value returned is in radians squared
    ///
    ///  This function throws an error unless 
    ///  all layers in the model have a Komolgorov
    ///  turbulence power spectrum.
    double aperture_averaged_tilt_phase_covariance(const emitter & emtr_a,
						   const emitter & emtr_b,
						   const aperture & ap,
						   double wavelength_meters,
						   int nsteps_in_integration=1000) const;
      
    ///////////////////////////////////////////
    ///  Get the perpendicular component of the tilt 
    ///  covariance between two beams, averaged 
    ///  over a circular aperture.
    ///  
    ///  The value returned is in radians squared
    ///
    ///  This function throws an error unless 
    ///  all layers in the model have a Komolgorov
    ///  turbulence power spectrum.
    double aperture_averaged_parallel_tilt_phase_covariance(const emitter & emtr_a,
							    const emitter & emtr_b,
							    const aperture & ap,
							    double wavelength_meters,
							    int nsteps_in_integration=1000) const;
      
    ///////////////////////////////////////////
    ///  Get the perpendicular component of the tilt 
    ///  covariance between two beams, averaged 
    ///  over a circular aperture.
    ///  
    ///  The value returned is in radians squared
    ///
    ///  This function throws an error unless 
    ///  all layers in the model have a Komolgorov
    ///  turbulence power spectrum.
    double aperture_averaged_perpendicular_tilt_phase_covariance(const emitter & emtr_a,
								 const emitter & emtr_b,
								 const aperture & ap,
								 double wavelength_meters,
								 int nsteps_in_integration=1000) const;
    
    ///////////////////////////////////////////
    ///  Get the piston removed phase covariance between two beams
    ///
    ///  The three_points pupil_location_one and pupil_location_two must
    ///  lie within the aperture diameter, whose origin is assumed to 
    ///  that of the three frame in the refractive atmospheric model
    ///
    ///  This quantity is computed from Tyler, 
    ///  JOSA 11 p 409 1994
    ///  
    ///  The value returned is in radians squared
    ///
    ///  This function throws an error unless 
    ///  all layers in the model have a Komolgorov
    ///  turbulence power spectrum.
    double phase_covariance(const emitter & emtr_a,
			    const three_point & pupil_location_one,
			    const emitter & emtr_b,	
			    const three_point & pupil_location_two,
			    const aperture & ap,
			    double wavelength_meters,
			    int nsteps_in_integration=1000) const;
    
    ///////////////////////////////////////////
    ///  Get the tilt covariance between two beams.
    ///
    ///  This quantity is computed from Tyler, 
    ///  JOSA 11 p 409 1994
    ///
    ///  The three_points pupil_location_one and pupil_location_two must
    ///  lie within the aperture diameter, whose origin is assumed to 
    ///  that of the three frame in the refractive atmospheric model
    ///  
    ///  The value returned is in radians squared
    ///
    ///  This function throws an error unless 
    ///  all layers in the model have a Komolgorov
    ///  turbulence power spectrum.
    double tilt_phase_covariance(const emitter & emtr_a,
				 const three_point & pupil_location_one,
				 const emitter & emtr_b,
				 const three_point & pupil_location_two,
				 const aperture & ap,
				 double wavelength_meters,
				 int nsteps_in_integration=1000) const;

    ///////////////////////////////////////////
    ///  Get a diffractive wavefront header 
    ///  suitable for this atmospheric model.
    ///  The three frame of the wavefront header
    ///  has an origin at the topmost layer of the
    ///  model, with z axis directed towards the
    ///  center of the aperture.  The wavefront
    ///  is sized according to the covering region
    ///  of the aperture and the padding required
    ///  by the propagation plan.  This padding itself
    ///  depends on the electromagnetic wavelength and
    ///  the distance from the uppermost layer to the 
    ///  aperture.  Finally, the transverse axes are 
    ///  chosen according to whether layer foreshortening
    ///  has been selected
    ///
    ///  Wavelength and pixel scale must be specified in meters
    template<class T>
      diffractive_wavefront_header<T> get_diffractive_wavefront_header(double wavelength, 
								       double pixscale,
								       const emitter * emtr,
								       const aperture * ap,
								       bool layer_foreshortening, 
								       propagation_plan * pplan) const;


    ///////////////////////////////////////////
    ///  Get a randomly generated set of 
    ///  refractive atmospheric layers
    ///
    ///  This function selects random velocities
    ///  for the layers according to the wind model
    ///  contained in this refractive atmospheric model.
    ///  This function determines the size of the layers
    ///  so that throughout the specified time_interval
    ///  and for each layer, every wavefront (as defined 
    ///  by their headers in dwfhs) - after being propagated 
    ///  along its own z axis to the height of the layer - 
    ///  will be contained by the layer.  
    ///
    ///  The flag layer_axes_wind_vector_aligned determines
    ///  whether one of the layer axes will be aligned with 
    ///  the wind vector of the layer.  This choice will save
    ///  RAM, but requires much more computation.
    ///
    ///  The flag layer_foreshortening determines whether
    ///  foreshortening of the layers will be used.  This is
    ///  another computationally expensive choice.
    ///
    ///  Layer pixel scales should be specified in
    ///  meters
    ///
    template<class T, class U>
      void get_refractive_atmospheric_layers(const vector<double> & layer_pixscales,
					     const subharmonic_method & subm,
					     const vector<diffractive_wavefront_header<T> > dwfhdrs,
					     const vector<three_vector> layer_wind_vectors,
					     double time_interval,
					     bool layer_axes_wind_vector_aligned,
					     bool layer_foreshortening,
					     vector<refractive_atmospheric_layer<U> > & ref_atm_layers) const;

    static int verbose_level;

  };

  // useful anonymous namespace functions for the turbulence
  // computations below
  namespace {

    void check_wavelength(double wavelength_meters){
      if(wavelength_meters<=0){
	cerr << "check_wavelength error\n"
	     << "wavelength " 
	     << wavelength_meters
	     << " out of range\n";
	throw(string("check_wavelength"));
      }
    }

    void check_zenith(double zenith_angle_degrees){
      if(zenith_angle_degrees<0 ||
	 zenith_angle_degrees>=90){
	cerr << "check_zenith error\n"
	     << "zenith angle " 
	     << zenith_angle_degrees
	     << " out of range\n";
	throw(string("check_zenith"));
      }
    }
    
    void check_wavelength_zenith(double wavelength_meters, 
				 double zenith_angle_degrees){
      check_wavelength(wavelength_meters);
      check_zenith(zenith_angle_degrees);
    }

    /*
    double get_cn2_coefficient(const power_spectrum & pspec){
      try{
	const isotropic_power_law_spectrum<power_law, null_inner_scale> & iso_pspec = 
	  dynamic_cast<const isotropic_power_law_spectrum<power_law, null_inner_scale> &>(pspec);
	// Defined in Sasiela eq 2.18
	double fac = 5*gamma_function(5/6.)/pow(2,4/3.)/pow(M_PI,3/2.)/9./gamma_function(2/3.);
	return(iso_pspec.get_power_law().get_coefficient()/2./M_PI/fac);
      } catch(...) {
	cerr << "get_cn2_coefficient error - power spectrum not Komolgorov\n";
	pspec.print(cerr, "pspectrum ");
      }
    }
    */

    double get_cn2_coefficient(double power_law_coefficient){
      static double fac = 2*M_PI*5*gamma_function(5/6.)/pow(2,4/3.)/pow(M_PI,3/2.)/9./gamma_function(2/3.);
      return(power_law_coefficient/fac);
    }
  }



  template<class T>
  pixel_array<T> refractive_atmospheric_model::caliph_A(const emitter & emtr_a,
							const emitter & emtr_b,
							double aperture_diameter_meters,
							double pixel_scale_meters) const {
    
    if(aperture_diameter_meters<=0){
      cerr << "refractive_atmospheric_model::caliph_A error - \n"
	   << "invalid aperture diameter "
	   << aperture_diameter_meters
	   << " passed to this function\n";
      throw(string("refractive_atmospheric_model::caliph_A"));
    }
			       
    if(pixel_scale_meters<=0){
      cerr << "refractive_atmospheric_model::caliph_A error -\n"
	   << "pixel scale "
	   << pixel_scale_meters
	   << " out of range\n";
      throw(string("refractive_atmospheric_model::caliph_A"));
    }

    T * data;
    vector<long> axes(2,(long)ceil(aperture_diameter_meters/pixel_scale_meters));
    try{
      data = new T[axes[0]*axes[1]];
    } catch(...) {
      cerr << "refractive_atmospheric_model::caliph_A - error allocating memory\n";
      throw(string("refractive_atmospheric_model::caliph_A"));
    }

    try{
      double secant_zenith_angle;
      double max_range_meters;
      double min_range_meters;
      three_vector little_omega;
      
      Tyler_get_constants(emtr_a,
		    emtr_b,
		    this->ground_ref_frame_,
		    aperture_diameter_meters,
		    secant_zenith_angle,
		    max_range_meters,
		    min_range_meters,
		    little_omega);

      double val;

      double Q;

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

      three_vector pixel_vector;
      double normalization_factor = 2*pixel_scale_meters/aperture_diameter_meters;
      int index;

      for(int i=-axes[1]/2; i<axes[1]/2+x_extrapix; i++){
	for(int j=-axes[0]/2; j<axes[0]/2+y_extrapix; j++){
	  pixel_vector = 
	    this->ground_ref_frame_.x()*((i+x_halfpix)*normalization_factor) +
	    this->ground_ref_frame_.y()*((j+y_halfpix)*normalization_factor);

	  index = (i+axes[1]/2)*axes[0]+j+axes[0]/2;

	  if(pixel_vector.length_squared()>1)
	    data[index]=0;
	  else {
	    val = 0;
	    for(int k=0; k<this->power_spectra_.size(); k++){
	      if(layer_heights_[k]>min_range_meters) continue;
	      
	      Q = (1 - layer_heights_[k]*secant_zenith_angle/min_range_meters) /
		(1 - layer_heights_[k]*secant_zenith_angle/max_range_meters);
	      
	      val += secant_zenith_angle*get_cn2_coefficient(power_spectra_[k]->get_coefficient()) * 
		pow((1 - layer_heights_[k]*secant_zenith_angle/max_range_meters),5/3.) *
		pow(Q,5/3.) *
		Tyler_F_1((pixel_vector + 
		     (little_omega*(layer_heights_[k]/(1-layer_heights_[k]/max_range_meters)))).length()
		    /Q);
	    }
	    data[index] = val;
	  }
	}
      }
    } catch(...) {
      cerr << "refractive_atmospheric_model::caliph_A error\n";
      throw(string("refractive_atmospheric_model::caliph_A"));
    }

    pixel_array<T> pixarr(axes,data);
    delete [] data;
    return(pixarr);
  }

  template<class T>
  pixel_array<T> refractive_atmospheric_model::caliph_B(const emitter & emtr_a,
							const emitter & emtr_b,
							double aperture_diameter_meters,
							double pixel_scale_meters) const {
    
    if(aperture_diameter_meters==0){
      cerr << "refractive_atmospheric_model::caliph_B error - \n"
	   << "invalid aperture diameter "
	   << aperture_diameter_meters
	   << " passed to this function\n";
      throw(string("refractive_atmospheric_model::caliph_B"));
    }
			       
    if(pixel_scale_meters<=0){
      cerr << "refractive_atmospheric_model::caliph_B error -\n"
	   << "pixel scale "
	   << pixel_scale_meters
	   << " out of range\n";
      throw(string("refractive_atmospheric_model::caliph_B"));
    }

    T * data;
    vector<long> axes(2,(long)ceil(aperture_diameter_meters/pixel_scale_meters));
    try{
      data = new T[axes[0]*axes[1]];
    } catch(...) {
      cerr << "refractive_atmospheric_model::caliph_B - error allocating memory\n";
      throw(string("refractive_atmospheric_model::caliph_B"));
    }

    try{
      double secant_zenith_angle;
      double max_range_meters;
      double min_range_meters;
      three_vector little_omega;
      
      Tyler_get_constants(emtr_a,
		    emtr_b,
		    this->ground_ref_frame_,
		    aperture_diameter_meters,
		    secant_zenith_angle,
		    max_range_meters,
		    min_range_meters,
		    little_omega);

      double val;

      double Q;

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

      three_vector pixel_vector;
      double normalization_factor = 2*pixel_scale_meters/aperture_diameter_meters;
      int index;

      for(int i=-axes[1]/2; i<axes[1]/2+x_extrapix; i++){
	for(int j=-axes[0]/2; j<axes[0]/2+y_extrapix; j++){
	  pixel_vector = 
	    this->ground_ref_frame_.x()*((i+x_halfpix)*normalization_factor) +
	    this->ground_ref_frame_.y()*((j+y_halfpix)*normalization_factor);

	  index = (i+axes[1]/2)*axes[0]+j+axes[0]/2;
	  if(pixel_vector.length_squared()>1)
	    data[index]=0;
	  else {
	    val = 0;
	    for(int k=0; k<this->power_spectra_.size(); k++){
	      if(layer_heights_[k]>min_range_meters) continue;
	      
	      Q = (1 - layer_heights_[k]*secant_zenith_angle/min_range_meters) /
		(1 - layer_heights_[k]*secant_zenith_angle/max_range_meters);
	      
	      val += secant_zenith_angle*get_cn2_coefficient(power_spectra_[k]->get_coefficient()) * 
		pow((1 - layer_heights_[k]*secant_zenith_angle/max_range_meters),5/3.) *
		Tyler_F_1(((Q*pixel_vector) - 
		     (little_omega*(layer_heights_[k]/(1-layer_heights_[k]/max_range_meters)))).length());
	    }
	    data[index] = val;
	  }
	}
      }
    } catch(...) {
      cerr << "refractive_atmospheric_model::caliph_B error\n";
      throw(string("refractive_atmospheric_model::caliph_B"));
    }
    pixel_array<T> pixarr(axes,data);
    delete [] data;
    return(pixarr);
  }

  template<class T>
  pixel_array<T> refractive_atmospheric_model::caliph_C(const emitter & emtr_a,
							const emitter & emtr_b,
							double aperture_diameter_meters,
							double pixel_scale_meters) const {
    
    if(aperture_diameter_meters==0){
      cerr << "refractive_atmospheric_model::caliph_C error - \n"
	   << "invalid aperture diameter "
	   << aperture_diameter_meters
	   << " passed to this function\n";
      throw(string("refractive_atmospheric_model::caliph_C"));
    }
			       
    if(pixel_scale_meters<=0){
      cerr << "refractive_atmospheric_model::caliph_C error -\n"
	   << "pixel scale "
	   << pixel_scale_meters
	   << " out of range\n";
      throw(string("refractive_atmospheric_model::caliph_C"));
    }

    T * data;
    vector<long> axes(2,(long)ceil(aperture_diameter_meters/pixel_scale_meters));
    try{
      data = new T[axes[0]*axes[1]];
    } catch(...) {
      cerr << "refractive_atmospheric_model::caliph_C - error allocating memory\n";
      throw(string("refractive_atmospheric_model::caliph_C"));
    }

    try{
      double secant_zenith_angle;
      double max_range_meters;
      double min_range_meters;
      three_vector little_omega;
      
      Tyler_get_constants(emtr_a,
		    emtr_b,
		    this->ground_ref_frame_,
		    aperture_diameter_meters,
		    secant_zenith_angle,
		    max_range_meters,
		    min_range_meters,
		    little_omega);

      double val;

      double Q;

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

      three_vector pixel_vector;
      double normalization_factor = 2*pixel_scale_meters/aperture_diameter_meters;
      int index;

      for(int i=-axes[1]/2; i<axes[1]/2+x_extrapix; i++){
	for(int j=-axes[0]/2; j<axes[0]/2+y_extrapix; j++){
	  pixel_vector = 
	    this->ground_ref_frame_.x()*((i+x_halfpix)*normalization_factor) +
	    this->ground_ref_frame_.y()*((j+y_halfpix)*normalization_factor);

	  index = (i+axes[1]/2)*axes[0]+j+axes[0]/2;

	  if(pixel_vector.length_squared()>1)
	    data[index]=0;
	  else {
	    val = 0;
	    for(int k=0; k<this->power_spectra_.size(); k++){
	      if(layer_heights_[k]>min_range_meters) continue;
	      
	      Q = (1 - layer_heights_[k]*secant_zenith_angle/min_range_meters) /
		(1 - layer_heights_[k]*secant_zenith_angle/max_range_meters);
	      
	      val += secant_zenith_angle*get_cn2_coefficient(power_spectra_[k]->get_coefficient()) * 
		pow((1 - layer_heights_[k]*secant_zenith_angle/max_range_meters),5/3.) *
		pow((pixel_vector*(1 - Q) +
		     (little_omega*(layer_heights_[k]/(1-layer_heights_[k]/max_range_meters)))).length(),5/3.);

	    }
	    data[index] = val;
	    if(!finite(val)){
	      cerr << "refractive_atmospheric_model::caliph_C error - value " 
		   << val 
		   << " not finite\n";
	      throw(string("refractive_atmospheric_model::caliph_C"));
	    }
	  }
	}
      }
    } catch(...) {
      cerr << "refractive_atmospheric_model::caliph_C error\n";
      throw(string("refractive_atmospheric_model::caliph_C"));
    }
    pixel_array<T> pixarr(axes,data);
    delete [] data;
    return(pixarr);
  }

  template<class T>
  pixel_array<T> refractive_atmospheric_model::caliph_G(const emitter & emtr_a,
							const emitter & emtr_b,
							double aperture_diameter_meters,
							double pixel_scale_meters) const {

    if(aperture_diameter_meters==0){
      cerr << "refractive_atmospheric_model::caliph_G error - \n"
	   << "invalid aperture diameter "
	   << aperture_diameter_meters
	   << " passed to this function\n";
      throw(string("refractive_atmospheric_model::caliph_G"));
    }
			       
    if(pixel_scale_meters<=0){
      cerr << "refractive_atmospheric_model::caliph_G error -\n"
	   << "pixel scale "
	   << pixel_scale_meters
	   << " out of range\n";
      throw(string("refractive_atmospheric_model::caliph_G"));
    }

    T * data;
    vector<long> axes(2,(long)ceil(aperture_diameter_meters/pixel_scale_meters));
    try{
      data = new T[axes[0]*axes[1]];
    } catch(...) {
      cerr << "refractive_atmospheric_model::caliph_G - error allocating memory\n";
      throw(string("refractive_atmospheric_model::caliph_G"));
    }

    try{
      double secant_zenith_angle;
      double max_range_meters;
      double min_range_meters;
      three_vector little_omega;
      
      Tyler_get_constants(emtr_a,
		    emtr_b,
		    this->ground_ref_frame_,
		    aperture_diameter_meters,
		    secant_zenith_angle,
		    max_range_meters,
		    min_range_meters,
		    little_omega);

      double ghat_arg;
      double val;

      double Q;

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

      three_vector pixel_vector;
      double normalization_factor = 2*pixel_scale_meters/aperture_diameter_meters;
      int index;

      for(int i=-axes[1]/2; i<axes[1]/2+x_extrapix; i++){
	for(int j=-axes[0]/2; j<axes[0]/2+y_extrapix; j++){
	  pixel_vector = 
	    this->ground_ref_frame_.x()*((i+x_halfpix)*normalization_factor) +
	    this->ground_ref_frame_.y()*((j+y_halfpix)*normalization_factor);

	  index = (i+axes[1]/2)*axes[0]+j+axes[0]/2;

	  if(pixel_vector.length_squared()>1)
	    data[index]=0;
	  else {
	    val = 0;
	    for(int k=0; k<this->power_spectra_.size(); k++){
	      if(layer_heights_[k]>min_range_meters) continue;
	      
	      Q = (1 - layer_heights_[k]*secant_zenith_angle/min_range_meters) /
		(1 - layer_heights_[k]*secant_zenith_angle/max_range_meters);
	      
	      ghat_arg = 
		((pixel_vector*Q) - (little_omega*(layer_heights_[k]/(1-layer_heights_[k]/max_range_meters)))).length();

	      val += secant_zenith_angle*get_cn2_coefficient(power_spectra_[k]->get_coefficient()) * 
		pow((1 - layer_heights_[k]*secant_zenith_angle/max_range_meters),5/3.) *
		Q * Tyler_G_hat(ghat_arg);
	    }

	    data[index] = 4*val;
	  }
	}
      }
    } catch(...) {
      cerr << "refractive_atmospheric_model::caliph_G error\n";
      throw(string("refractive_atmospheric_model::caliph_G"));
    }
    pixel_array<T> pixarr(axes,data);
    delete [] data;
    return(pixarr);
  }

  template<class T>
  pixel_array<T> refractive_atmospheric_model::caliph_H(const emitter & emtr_a,
							const emitter & emtr_b,
							double aperture_diameter_meters,
							double pixel_scale_meters) const {
    
    if(aperture_diameter_meters==0){
      cerr << "refractive_atmospheric_model::caliph_H error - \n"
	   << "invalid aperture diameter "
	   << aperture_diameter_meters
	   << " passed to this function\n";
      throw(string("refractive_atmospheric_model::caliph_H"));
    }
			       
    if(pixel_scale_meters<=0){
      cerr << "refractive_atmospheric_model::caliph_H error -\n"
	   << "pixel scale "
	   << pixel_scale_meters
	   << " out of range\n";
      throw(string("refractive_atmospheric_model::caliph_H"));
    }

    T * data;
    vector<long> axes(2,(long)ceil(aperture_diameter_meters/pixel_scale_meters));
    try{
      data = new T[axes[0]*axes[1]];
    } catch(...) {
      cerr << "refractive_atmospheric_model::caliph_H - error allocating memory\n";
      throw(string("refractive_atmospheric_model::caliph_H"));
    }

    try{
      double secant_zenith_angle;
      double max_range_meters;
      double min_range_meters;
      three_vector little_omega;
      
      Tyler_get_constants(emtr_a,
		    emtr_b,
		    this->ground_ref_frame_,
		    aperture_diameter_meters,
		    secant_zenith_angle,
		    max_range_meters,
		    min_range_meters,
		    little_omega);

      double ghat_arg;
      double val;

      double Q;

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

      three_vector pixel_vector;
      double normalization_factor = 2*pixel_scale_meters/aperture_diameter_meters;
      int index;

      for(int i=-axes[1]/2; i<axes[1]/2+x_extrapix; i++){
	for(int j=-axes[0]/2; j<axes[0]/2+y_extrapix; j++){
	  pixel_vector = 
	    this->ground_ref_frame_.x()*((i+x_halfpix)*normalization_factor) +
	    this->ground_ref_frame_.y()*((j+y_halfpix)*normalization_factor);

	  index = (i+axes[1]/2)*axes[0]+j+axes[0]/2;

	  if(pixel_vector.length_squared()>1)
	    data[index]=0;
	  else {
	    val = 0;
	    for(int k=0; k<this->power_spectra_.size(); k++){
	      if(layer_heights_[k]>min_range_meters) continue;
	      
	      Q = (1 - layer_heights_[k]*secant_zenith_angle/min_range_meters) /
		(1 - layer_heights_[k]*secant_zenith_angle/max_range_meters);
	      
	      ghat_arg = 
		((pixel_vector*Q) - (little_omega*(layer_heights_[k]/(1-layer_heights_[k]/max_range_meters)))).length();

	      val += secant_zenith_angle*get_cn2_coefficient(power_spectra_[k]->get_coefficient()) * 
		pow((1 - layer_heights_[k]*secant_zenith_angle/max_range_meters),5/3.) *
		little_omega.length()*(layer_heights_[k]/(1-layer_heights_[k]/max_range_meters))
		* Tyler_G_hat(ghat_arg);
	    }
	    data[index] = 4*val;
	  }
	}
      }
    } catch(...) {
      cerr << "refractive_atmospheric_model::caliph_H error\n";
      throw(string("refractive_atmospheric_model::caliph_H"));
    }
    pixel_array<T> pixarr(axes,data);
    delete [] data;
    return(pixarr);
  }

  template<class T>
  pixel_array<T> refractive_atmospheric_model::caliph_I(const emitter & emtr_a,
							const emitter & emtr_b,
							double aperture_diameter_meters,
							double pixel_scale_meters) const {
    
    if(aperture_diameter_meters==0){
      cerr << "refractive_atmospheric_model::caliph_I error - \n"
	   << "invalid aperture diameter "
	   << aperture_diameter_meters 
	   << " passed to this function\n";
      throw(string("refractive_atmospheric_model::caliph_I"));
    }
			       
    if(pixel_scale_meters<=0){
      cerr << "refractive_atmospheric_model::caliph_I error -\n"
	   << "pixel scale "
	   << pixel_scale_meters
	   << " out of range\n";
      throw(string("refractive_atmospheric_model::caliph_I"));
    }

    T * data;
    vector<long> axes(2,(long)ceil(aperture_diameter_meters/pixel_scale_meters));
    try{
      data = new T[axes[0]*axes[1]];
    } catch(...) {
      cerr << "refractive_atmospheric_model::caliph_I - error allocating memory\n";
      throw(string("refractive_atmospheric_model::caliph_I"));
    }

    try{
      double secant_zenith_angle;
      double max_range_meters;
      double min_range_meters;
      three_vector little_omega;
      
      Tyler_get_constants(emtr_a,
		    emtr_b,
		    this->ground_ref_frame_,
		    aperture_diameter_meters,
		    secant_zenith_angle,
		    max_range_meters,
		    min_range_meters,
		    little_omega);

      double ghat_arg;
      double val;

      double Q;
      
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

      three_vector pixel_vector;
      double normalization_factor = 2*pixel_scale_meters/aperture_diameter_meters;
      int index;

      for(int i=-axes[1]/2; i<axes[1]/2+x_extrapix; i++){
	for(int j=-axes[0]/2; j<axes[0]/2+y_extrapix; j++){
	  pixel_vector = 
	    this->ground_ref_frame_.x()*((i+x_halfpix)*normalization_factor) +
	    this->ground_ref_frame_.y()*((j+y_halfpix)*normalization_factor);

	  index = (i+axes[1]/2)*axes[0]+j+axes[0]/2;

	  if(pixel_vector.length_squared()>1)
	    data[index]=0;
	  else {
	    val = 0;
	    for(int k=0; k<this->power_spectra_.size(); k++){
	      if(layer_heights_[k]>min_range_meters) continue;
	      
	      Q = (1 - layer_heights_[k]*secant_zenith_angle/min_range_meters) /
		(1 - layer_heights_[k]*secant_zenith_angle/max_range_meters);
	      
	      ghat_arg = 
		(pixel_vector + (little_omega*(layer_heights_[k]/(1-layer_heights_[k]/max_range_meters)))).length()/Q;

	      val += secant_zenith_angle*get_cn2_coefficient(power_spectra_[k]->get_coefficient()) * 
		pow((1 - layer_heights_[k]*secant_zenith_angle/max_range_meters),5/3.) *
		pow(Q,2/3.) * 
		Tyler_G_hat(ghat_arg);
	    }
	    data[index] = 4*val;
	  }
	}
      }
    } catch(...) {
      cerr << "refractive_atmospheric_model::caliph_I error\n";
      throw(string("refractive_atmospheric_model::caliph_I"));
    }
    pixel_array<T> pixarr(axes,data);
    delete [] data;
    return(pixarr);
  }

  template<class T>
  pixel_array<T> refractive_atmospheric_model::caliph_J(const emitter & emtr_a,
							const emitter & emtr_b,
							double aperture_diameter_meters,
							double pixel_scale_meters) const {
    
    if(aperture_diameter_meters==0){
      cerr << "refractive_atmospheric_model::caliph_J error - \n"
	   << "invalid aperture diameter "
	   << aperture_diameter_meters
	   << " passed to this function\n";
      throw(string("refractive_atmospheric_model::caliph_J"));
    }
			       
    if(pixel_scale_meters<=0){
      cerr << "refractive_atmospheric_model::caliph_J error -\n"
	   << "pixel scale "
	   << pixel_scale_meters
	   << " out of range\n";
      throw(string("refractive_atmospheric_model::caliph_J"));
    }

    T * data;
    vector<long> axes(2,(long)ceil(aperture_diameter_meters/pixel_scale_meters));
    try{
      data = new T[axes[0]*axes[1]];
    } catch(...) {
      cerr << "refractive_atmospheric_model::caliph_J - error allocating memory\n";
      throw(string("refractive_atmospheric_model::caliph_J"));
    }

    try{
      double secant_zenith_angle;
      double max_range_meters;
      double min_range_meters;
      three_vector little_omega;
      
      Tyler_get_constants(emtr_a,
		    emtr_b,
		    this->ground_ref_frame_,
		    aperture_diameter_meters,
		    secant_zenith_angle,
		    max_range_meters,
		    min_range_meters,
		    little_omega);

      double ghat_arg;
      double val;

      double Q;

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

      three_vector pixel_vector;
      double normalization_factor = 2*pixel_scale_meters/aperture_diameter_meters;
      int index;

      for(int i=-axes[1]/2; i<axes[1]/2+x_extrapix; i++){
	for(int j=-axes[0]/2; j<axes[0]/2+y_extrapix; j++){
	  pixel_vector = 
	    this->ground_ref_frame_.x()*((i+x_halfpix)*normalization_factor) +
	    this->ground_ref_frame_.y()*((j+y_halfpix)*normalization_factor);

	  index = (i+axes[1]/2)*axes[0]+j+axes[0]/2;

	  if(pixel_vector.length_squared()>1)
	    data[index]=0;
	  else {
	    val = 0;
	    for(int k=0; k<this->power_spectra_.size(); k++){
	      if(layer_heights_[k]>min_range_meters) continue;
	      
	      Q = (1 - layer_heights_[k]*secant_zenith_angle/min_range_meters) /
		(1 - layer_heights_[k]*secant_zenith_angle/max_range_meters);
	      
	      ghat_arg = 
		(pixel_vector + (little_omega*(layer_heights_[k]/(1-layer_heights_[k]/max_range_meters)))).length()/Q;

	      val += secant_zenith_angle*get_cn2_coefficient(power_spectra_[k]->get_coefficient()) * 
		pow((1 - layer_heights_[k]*secant_zenith_angle/max_range_meters),5/3.) *
		pow(Q, 2/3.) *
		little_omega.length()*(layer_heights_[k]/(1-layer_heights_[k]/max_range_meters)) * 
		Tyler_G_hat(ghat_arg);
	    }
	    data[index] = 4*val;
	  }
	}
      }
    } catch(...) {
      cerr << "refractive_atmospheric_model::caliph_J error\n";
      throw(string("refractive_atmospheric_model::caliph_J"));
    }
    pixel_array<T> pixarr(axes,data);
    delete [] data;
    return(pixarr);
  }

  template<class T>
  diffractive_wavefront_header<T> 
  refractive_atmospheric_model::get_diffractive_wavefront_header(double wavelength, 
								 double pixscale,
								 const emitter * emtr,
								 const aperture * ap,
								 bool layer_foreshortening, 
								 propagation_plan * pplan) const {

    // Define the wavefront reference frame.  The choice of transverse axes
    // depends on whether we choose to foreshorten the layer 
    three_vector emitter_direction_vector = -1*emtr->get_emission_vector(*ap);

    double cos_angle = dot_product(-1*emitter_direction_vector, ap->z());
    three_frame wf_frame(*ap-(layer_heights_[0]/cos_angle)*emitter_direction_vector,
			 -1*ground_ref_frame_.x(), 
			 -1*ground_ref_frame_.y(), 
			 -1*ground_ref_frame_.z());
    if(layer_foreshortening) {
      three_vector axis_in_layer = cross_product(ground_ref_frame_.z(), emitter_direction_vector);
      if(axis_in_layer.length()>three_frame::precision){
	wf_frame = three_frame(*ap-(layer_heights_[0]/cos_angle)*emitter_direction_vector,
			       cross_product(axis_in_layer, emitter_direction_vector),
			       axis_in_layer, 
			       emitter_direction_vector);
      }
    } else {
      three_vector rotation_vector = cross_product(emitter_direction_vector, 
						   wf_frame.z());
      if(rotation_vector.length()!=0){
	three_rotation trot(wf_frame, 
			    rotation_vector, 
			    rotation_vector.length());
	trot.transform(wf_frame);
      }
    }

    // Get the aperture covering region for this frame, and
    // round its dimensions up to an integral number of wavefront
    // pixels
    rectangular_region wf_region = ap->get_covering_region(wf_frame);
    wf_region = rectangular_region(wf_region, pixscale);

    // Find the size of the wavefront array for which we need 
    // valid data.
    vector<three_point> region_corners = wf_region.get_corners();
    long dimen = (long)ceil((region_corners[1]-region_corners[0]).length()/pixscale);
    if((region_corners[3]-region_corners[0]).length() >
       (region_corners[1]-region_corners[0]).length())
      dimen = (long)ceil((region_corners[3]-region_corners[0]).length()/pixscale);

    vector<long> wf_axes(2, dimen);

    // set up the curvature and pixel scale
    double curvature = 0;
    double init_pixscale = pixscale;
    const spherical_wave_emitter * swe;
    if(swe=dynamic_cast<const spherical_wave_emitter *>(emtr)){
      curvature = 1/(wf_frame-*swe).length();
      init_pixscale *= (wf_frame-*swe).length()/(*ap-*swe).length();
    }

    // Declare the wavefront header
    diffractive_wavefront_header<T> dwf(wf_axes, wf_frame, wavelength, init_pixscale);

    dwf.set_curvature(curvature);

    // Declare a wavefront header for which we have added in the 
    // padding required by the propagation plan
    diffractive_wavefront_header<T> padded_dwf = pplan->pad(dwf, layer_heights_[0]);

    return(padded_dwf);
  }

  template<class T, class U>
  void refractive_atmospheric_model::
  get_refractive_atmospheric_layers(const vector<double> & layer_pixscales,
				    const subharmonic_method & subm,
				    const vector<diffractive_wavefront_header<T> > dwfhdrs,
				    const vector<three_vector> layer_wind_vectors,
				    double time_interval,
				    bool layer_axes_wind_vector_aligned,
				    bool layer_foreshortening,
				    vector<refractive_atmospheric_layer<U> > & ref_atm_layers) const {

    if(time_interval < 0){
      cerr << "refractive_atmospheric_model::get_refractive_atmospheric_layers error - "
	   << "time interval " << time_interval << " provided to this function "
	   << "is less than or equal to zero\n";
      throw(string("refractive_atmospheric_model::get_refractive_atmospheric_layers"));
    }

    if(layer_pixscales.size() != power_spectra_.size()){
      cerr << "refractive_atmospheric_model::get_refractive_atmospheric_layers error - "
	   << layer_pixscales.size() << " layer pixel scales were provided, but there are "
	   << power_spectra_.size() << " layers in the model\n";
      throw(string("refractive_atmospheric_model::get_refractive_atmospheric_layers"));
    }

    if(layer_wind_vectors.size()!=this->get_number_of_layers()){
      cerr << "refractive_atmospheric_model::get_refractive_atmospheric_layers error - "
	   << "number of layer wind vectors " << layer_wind_vectors.size()
	   << " does not match the number of layers " << this->get_number_of_layers()
	   << " in the atmospheric model\n";
      throw(string("refractive_atmospheric_model::get_refractive_atmospheric_layers"));
    }

    // Here for each layer we do the following:
    // Get a random wind vector
    // Initialize the three_frame of the layer
    // Initialize the rectangular_region that will serve to define the size of the layer
    int nlayers = power_spectra_.size();
    vector<three_frame> layer_three_frames(nlayers, ground_ref_frame_);
    vector<rectangular_region> layer_rec_regions;
    double windspeed, x_wind_cmpnt, y_wind_cmpnt, angle;
    for(int i=0; i<nlayers; i++){

      if(refractive_atmospheric_model::verbose_level) 
	cout << "refractive_atmospheric_model::get_refractive_atmospheric_layers - initializing layer three frame...";

      // If the layer axes are to be aligned with the wind vector, we
      // redefine the layer three frame to have this alignment.
      if(layer_axes_wind_vector_aligned){
	layer_three_frames[i] = three_frame(ground_ref_frame_,
					    cross_product(layer_wind_vectors[i], ground_ref_frame_.z()), 
					    layer_wind_vectors[i], 
					    ground_ref_frame_.z());
      } 

      // Bump the layer three frame up to the correct altitude
      layer_three_frames[i] += layer_heights_[i]*ground_ref_frame_.z();

      if(refractive_atmospheric_model::verbose_level) cout << " initializing rectangular region...";

      // Initialize the layer's rectangular region by translating the
      // first wavefront down to the layer height, and then finding
      // the rectangular region that will encompass this wavefront
      diffractive_wavefront_header<T> tmp_header = dwfhdrs[0];
      if(refractive_atmospheric_model::verbose_level) cout << "getting distance...";
      tmp_header.three_point::operator=(get_ray_plane_intersection(tmp_header, tmp_header.z(), 
								   layer_three_frames[i], layer_three_frames[i].z()));

      if(refractive_atmospheric_model::verbose_level) cout << "making rec region...";
      layer_rec_regions.push_back(tmp_header.get_covering_region(layer_three_frames[i], layer_foreshortening));
      
      if(refractive_atmospheric_model::verbose_level) cout << "initialization complete\n";
    }


    // Here for each layer we loop over the wavefront headers finding the
    // rectangular region containing the wavefront at the layer height
    // and forming the union with the region from the previous iteration
    vector<three_point> tp;
    ref_atm_layers.clear();
    vector<long> layer_axes(2);
    diffractive_wavefront_header<T> tmp_header;
    for(int i=0; i<nlayers; i++){
      if(refractive_atmospheric_model::verbose_level) cout << "Layer " << i << " constructing rectangular region...";
      for(int j=1; j<dwfhdrs.size(); j++){


	tmp_header = dwfhdrs[j];
	tmp_header.three_point::operator=(get_ray_plane_intersection(dwfhdrs[j], dwfhdrs[j].z(), 
								     layer_three_frames[i], layer_three_frames[i].z()));
 
	layer_rec_regions[i] = 
	  region_union(layer_rec_regions[i], 
		       tmp_header.get_covering_region(layer_three_frames[i], layer_foreshortening));
      }

      if(refractive_atmospheric_model::verbose_level) cout << " enlarging region for simulation...";

      // Next we will enlarge the rectangular regions to account for the
      // duration of the simulation and the wind speed of the layers.
      three_vector lateral_displacement = time_interval*layer_wind_vectors[i];
      tp = layer_rec_regions[i].get_corners();
      for(int j=0; j<4; j++) tp[j] -= lateral_displacement;
      layer_rec_regions[i] = region_union(layer_rec_regions[i], rectangular_region(tp));

      if(refractive_atmospheric_model::verbose_level)
	layer_rec_regions[i].print(cout, "lrr ");

      if(refractive_atmospheric_model::verbose_level) cout << " rounding up region...";

      // Round the rectangular regions up to size evenly 
      // divisible by the layer pixel scale
      layer_rec_regions[i] = rectangular_region(layer_rec_regions[i], layer_pixscales[i]);
      layer_three_frames[i].three_point::operator=(layer_rec_regions[i].get_center());

      if(refractive_atmospheric_model::verbose_level)
	layer_three_frames[i].print(cout, "layer three frame ");

      // Define the layer axes.  If layer_axes_wind_vector_aligned ==
      // true, we choose the convention that axes[0] is along the wind
      // vector.  If layer_axes_wind_vector_aligned == false, we
      // choose the convention that axes[0] is along the x axis defined
      // by the ground_ref_frame_
      //
      // Note - here we add one more pixel to the axes to beat rounding
      // issues.  We also pad by 10 percent to avoid ringing at the edges
      // due to the subharmonic correction

      if(refractive_atmospheric_model::verbose_level) cout << " getting layer axes...";
      double padding_factor = 1.2;
      tp = layer_rec_regions[i].get_corners();
      three_vector tva = tp[1]-tp[0];
      three_vector tvb = tp[3]-tp[0];
      layer_axes[0] = (long)(padding_factor*(1+ceil(tvb.length()/layer_pixscales[i])));
      layer_axes[1] = (long)(padding_factor*(1+ceil(tva.length()/layer_pixscales[i])));
      if(layer_axes_wind_vector_aligned){
	if(cross_product(tva, layer_wind_vectors[i]).length()<three_frame::precision){
	  layer_axes[0] = (long)(padding_factor*(1+ceil(tva.length()/layer_pixscales[i])));
	  layer_axes[1] = (long)(padding_factor*(1+ceil(tvb.length()/layer_pixscales[i])));
	}
      } else {
	if(cross_product(tva, ground_ref_frame_.x()).length()<three_frame::precision){
	  layer_axes[0] = (long)(padding_factor*(1+ceil(tva.length()/layer_pixscales[i])));
	  layer_axes[1] = (long)(padding_factor*(1+ceil(tvb.length()/layer_pixscales[i])));
	}
      }

      // Finally we're ready to make the layers.  
      if(refractive_atmospheric_model::verbose_level)
	cout << "making layer " << ref_atm_layers.size()
		<< " of size " << layer_axes[0] << " x " << layer_axes[1] 
		<< " with pixscale " << layer_pixscales[i] << endl;

      /*
      refractive_atmospheric_layer<T> ral;
      power_spectra_[i]->get_refractive_atmospheric_layer(
      			layer_axes, layer_pixscales[i], subm, ral);
      ref_atm_layers.push_back(ral);
      */
      //ref_atm_layers.push_back(power_spectra_[i]->get_refractive_atmospheric_layer(
      //			layer_axes, layer_pixscales[i], subm));


      ref_atm_layers.push_back(refractive_atmospheric_layer<U>(power_spectra_[i],
							       subm,
							       layer_axes, 
							       layer_pixscales[i]));
      ref_atm_layers[i].set_wind_vector(layer_wind_vectors[i]);
      ref_atm_layers[i].three_frame::operator=(layer_three_frames[i]);
      ref_atm_layers[i].set_foreshortening(false);
    }
  }

  /*

  /// A class to represent a refractive atmosphere     
  /// as an array of refractive_atmospheric_layers.

  class refractive_atmosphere :
    //public atmosphere 
    virtual public AO_sim_base, 
    public fits_header_data {

    private:
    
    static const bool factory_registration;

    protected:
  
    /// The array of refractive_atmospheric_layers
    vector<refractive_atmospheric_layer > ref_atm_layers;

    public:

    ///////////////////////////////////////////
    ///  Null constructor
    refractive_atmosphere(){};

    ///////////////////////////////////////////
    ///  Copy constructor
    refractive_atmosphere(const refractive_atmosphere & ref_atm);

    ///////////////////////////////////////////
    ///  Construct from a file
    refractive_atmosphere(const char * filename);

    ///////////////////////////////////////////
    ///  Construct from an iofits object
    refractive_atmosphere(const iofits & iof);

    ///////////////////////////////////////////
    ///  Construct from an atmospheric_model
    refractive_atmosphere(const refractive_atmospheric_model & ref_atm_model, 
			  double pixel_scale, 
			  const vector<long> & axes,
			  const subharmonic_method & subm);

    ///////////////////////////////////////////
    ///  Construct from atmospheric_layers
    refractive_atmosphere(const vector<refractive_atmospheric_layer> & ref_atm_layers);

    ///////////////////////////////////////////
    ///  Destructor
    ~refractive_atmosphere(){};

    ///////////////////////////////////////////
    ///  Operator = 
    refractive_atmosphere & operator=(const refractive_atmosphere & ref_atm);

    ///////////////////////////////////////////
    ///  Read from file
    void read(const char * filename);

    ///////////////////////////////////////////
    ///  Read from iofits object
    void read(const iofits & iof);

    ///////////////////////////////////////////
    ///  Write to file
    void write(const char * filename) const;

    ///////////////////////////////////////////
    ///  Write to iofits object
    void write(iofits & iof) const;

    ///////////////////////////////////////////
    ///  Print
    void print(ostream & os, const char * prefix="") const;

    ///////////////////////////////////////////
    ///  Get the number of layers
    long size() const {return(ref_atm_layers.size());};
 
    static int verbose_level;
 
  };

  */

}

#endif
