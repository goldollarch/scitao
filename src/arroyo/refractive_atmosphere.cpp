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
#include "iofits.h"
#include "fits_factory.h"
#include "fits_header_data.h"
#include "power_spectrum.h"
#include "structure_function.h"
#include "refractive_atmospheric_layer.h"
#include "refractive_atmosphere.h"
#include "wind_model.h"
#include "emitter.h"
#include "subharmonic_method.h"

using namespace std;

namespace Arroyo {

  namespace factory_register {
    const fits_keyval_set & get_refractive_atmospheric_model_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "refractive atmospheric model"));
      return *fkvs;
    }
    
    AO_sim_base * create_refractive_atmospheric_model(const iofits & iof) {
      return new refractive_atmospheric_model(iof);
    }
  }
  
  const bool refractive_atmospheric_model::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_refractive_atmospheric_model_keyval_set(), 
				      factory_register::create_refractive_atmospheric_model);

  int refractive_atmospheric_model::verbose_level = 0;

  refractive_atmospheric_model::
  refractive_atmospheric_model(const refractive_atmospheric_model & ref_atm_model) {
    this->operator=(ref_atm_model);
  }

  refractive_atmospheric_model::
  refractive_atmospheric_model(const char * filename) {
    this->read(filename);
  } 

  refractive_atmospheric_model::
  refractive_atmospheric_model(const iofits & iof) {
    this->read(iof);
  } 

  namespace {
    struct sort_pred {
      bool operator()(const std::pair<double, power_spectrum *> & left, 
		      const std::pair<double, power_spectrum *> & right){
	return left.first < right.first;
      }
    };
  }

  refractive_atmospheric_model::
  refractive_atmospheric_model(const vector<power_spectrum *> & power_spectra,
			       const vector<double> & layer_heights,
			       const three_frame & ground_ref_frame)  {
    if(power_spectra.size()!=layer_heights.size()){
      cerr << "refractive_atmospheric_model::refractive_atmospheric_model error - "
	   << "mismatched sizes for vectors of power_spectrum and heights\n";
      throw(string("refractive_atmospheric_model::refractive_atmospheric_model"));
    }

    std::vector<std::pair<double, power_spectrum *> > sortable(layer_heights.size());

    std::vector<double>::const_iterator cdit;
    std::vector<power_spectrum *>::const_iterator cpit;
    std::vector<std::pair<double, power_spectrum * > >::iterator it;
    for(cdit = layer_heights.begin(),
	  cpit = power_spectra.begin(),
	  it = sortable.begin();
	cdit!=layer_heights.end();
	cdit++,cpit++,it++){
      it->first = *(cdit);
      it->second = *(cpit);
    }

    std::sort(sortable.begin(),sortable.end(),sort_pred());

    this->layer_heights_.resize(layer_heights.size());
    this->power_spectra_.resize(layer_heights.size());

    std::vector<double>::iterator dit;
    std::vector<power_spectrum *>::iterator pit;
    std::vector<std::pair<double, power_spectrum * > >::const_iterator cit;
    for(dit = this->layer_heights_.begin(),
	  pit = this->power_spectra_.begin(),
	  cit = sortable.begin();
	dit!=this->layer_heights_.end();
	dit++,pit++,cit++){
      *dit = cit->first;
      *pit = power_spectrum::power_spectrum_factory(cit->second);
    }

    /*
    vector<double> tmp_layer_heights(layer_heights);

    // Sort this->layer_heights_ and this->power_spectra_
    // from highest to lowest altitude
    int max_height_index;
    double max_height;
    for(int i=0; i<layer_heights.size(); i++){
      max_height = -.5;
      for(int j=0; j<layer_heights.size(); j++){
	if(tmp_layer_heights[j]>max_height){
	  max_height = tmp_layer_heights[j];
	  max_height_index = j;
	}
      }

      if(tmp_layer_heights[max_height_index]<0){
	cerr << "refractive_atmospheric_model::refractive_atmospheric_model error - "
	     << "layer height " 
	     << tmp_layer_heights[max_height_index]
	     << " out of range\n";
	throw(string("refractive_atmospheric_model::refractive_atmospheric_model"));
      }

      layer_heights_.push_back(tmp_layer_heights[max_height_index]);
      power_spectra_.push_back(power_spectrum::power_spectrum_factory(power_spectra[max_height_index]));
      tmp_layer_heights[max_height_index] = -1;

    }
    */

    ground_ref_frame_ = ground_ref_frame;
  }

  refractive_atmospheric_model::
  refractive_atmospheric_model(const vector<double> & layer_heights_meters,
			       const vector<double> & Cn2_dz_m_2_3,
			       const three_frame & ground_ref_frame){

    if(layer_heights_meters.size()!=Cn2_dz_m_2_3.size()){
      cerr << "refractive_atmospheric_model::refractive_atmospheric_model error - "
	   << "mismatched sizes for vectors of Cn2 dz and heights\n";
      throw(string("refractive_atmospheric_model::refractive_atmospheric_model"));
    }

    // Defined in Sasiela eq 2.18
    double cn2_to_power_law_coefficient_conversion_factor = 
      2*M_PI*5*gamma_function(5/6.)/pow(2,4/3.)/pow(M_PI,3/2.)/9./gamma_function(2/3.);
    
    std::vector<power_spectrum *> pspecs(layer_heights_meters.size());
    double exponent = -11/3.;
    for(int i=0; i<layer_heights_meters.size(); i++){
      pspecs[i] =
	new isotropic_power_law_spectrum<power_law,null_inner_scale>(power_law(exponent, 
									       Cn2_dz_m_2_3[i]*
									       cn2_to_power_law_coefficient_conversion_factor),
								     null_inner_scale());	  	  
    }
    
    this->operator=(refractive_atmospheric_model(pspecs,
						 layer_heights_meters,
						 ground_ref_frame));
      
    for(int i=0; i<pspecs.size(); i++)
      delete pspecs[i];
  }

  refractive_atmospheric_model::~refractive_atmospheric_model(){
    for(int i=0; i<power_spectra_.size(); i++)
      delete power_spectra_[i];
  }

  refractive_atmospheric_model & refractive_atmospheric_model::
  operator=(const refractive_atmospheric_model & ref_atm_model){
    if(this==&ref_atm_model)
      return(*this);

    power_spectra_.resize(ref_atm_model.power_spectra_.size());
    for(int i=0; i<ref_atm_model.power_spectra_.size(); i++)
      power_spectra_[i] = power_spectrum::power_spectrum_factory(ref_atm_model.power_spectra_[i]);

    layer_heights_ = ref_atm_model.layer_heights_;
    ground_ref_frame_ = ref_atm_model.ground_ref_frame_;

    return(*this);
  }

  void refractive_atmospheric_model::read_common_data(const iofits & iof){
    long npspec;
    string comment;
    iof.read_key("NPSPEC", npspec, comment);
    layer_heights_.resize(npspec);
    power_spectra_.resize(npspec);
    
    ground_ref_frame_.read(iof);
    
    stringstream ss;
    for(int i=0; i<npspec; i++){
      ss.str("");
      ss << "HGHT" << i;
      iof.read_key(ss.str().c_str(), layer_heights_[i], comment);
    }
    
    iof.movrel_hdu(1);
    
    for(int i=0; i<npspec; i++)
      power_spectra_[i] = power_spectrum::power_spectrum_factory(iof);
  }

  void refractive_atmospheric_model::read(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "refractive_atmospheric_model::read - "
	   << "error opening file " << filename << endl;
      throw(string("refractive_atmospheric_model::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "refractive_atmospheric_model::read - "
	   << "error reading "
	   << this->unique_name() << " from file "
	   << filename << endl;
      throw(string("refractive_atmospheric_model::read"));
    }
  }

  void refractive_atmospheric_model::read(const iofits & iof){
    if(!iof.key_exists("TYPE")){
      cerr << "refractive_atmospheric_model::read error - "
	   << "unrecognized type of file\n";
      throw(string("refractive_atmospheric_model::read"));
    }
    string type, comment;
    iof.read_key("TYPE", type, comment);
    if(type!=this->unique_name()){
      cerr << "refractive_atmospheric_model::read error - file of type " 
	   << type << " rather than of type "
	   << this->unique_name() << endl;
      throw(string("refractive_atmospheric_model::read"));
    }

    this->read_common_data(iof);
  }

  void refractive_atmospheric_model::write_common_data(iofits & iof) const {
    long npspec = power_spectra_.size();
    string comment;
    iof.write_key("NPSPEC", npspec, "number of power spectra");
    ground_ref_frame_.write(iof);
    
    stringstream ss;
    for(int i=0; i<npspec; i++){
      ss.str("");
      ss << "HGHT" << i;
      iof.write_key(ss.str().c_str(), layer_heights_[i], comment);
    }
    
    for(int i=0; i<npspec; i++)
      power_spectra_[i]->write(iof);
  }

  void refractive_atmospheric_model::write(const char * filename) const {
    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "refractive_atmospheric_model::write - "
	   << "error opening file " << filename << endl;
      throw(string("refractive_atmospheric_model::write"));
    }
    try{this->write(iof);}
    catch(...){
      cerr << "refractive_atmospheric_model::write - "
	   << "error writing "
	   << this->unique_name() << " to file " 
	   << filename << endl;
      throw(string("refractive_atmospheric_model::write"));
    }
  }

  void refractive_atmospheric_model::write(iofits & iof) const {

    fits_header_data<char> tmphdr;
    tmphdr.write(iof);

    string type = this->unique_name();
    string comment = "object type";
    iof.write_key("TYPE", type, comment);

    this->write_common_data(iof);
  }
  
  void refractive_atmospheric_model::
  print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << this->unique_name()
       << "/" << "object type" << endl;
    os << prefix << "NPSPEC     = " << setw(vlspc) << power_spectra_.size()
       << "/" << "number of power spectra" << endl;
    ground_ref_frame_.print(os, prefix);

    stringstream ss;
    for(int i=0; i<power_spectra_.size(); i++){
      ss.str("");
      ss << "HGHT" << i << setw(8) << "= ";
      os << prefix << ss.str().c_str() << setw(vlspc) << layer_heights_[i]
	 << "/" << "height of layer (meters)" << endl;
      power_spectra_[i]->print(os, prefix);
    }
  }



  double refractive_atmospheric_model::
  turbulence_moment(double moment, 
		    double zenith_angle_degrees) const {

    try{check_zenith(zenith_angle_degrees);}
    catch(...){
      cerr << "refractive_atmospheric_model::turbulence_moment error\n";
      throw(string("refractive_atmospheric_model::turbulence_moment"));
    }
    double secant_zenith = 1/cos(zenith_angle_degrees*M_PI/180.);
    double val = 0;
    int npowerspec = this->power_spectra_.size();
    for(int i=0; i<npowerspec; i++){

      // only include the first layer for the mu_0 moment
      if((layer_heights_[i]==0 && moment==0) || layer_heights_[i]!=0)
	val += get_cn2_coefficient(power_spectra_[i]->get_coefficient()) *
	  pow(layer_heights_[i], moment) *
	  pow(secant_zenith, moment+1);
    }
    return(val);
  }

  double refractive_atmospheric_model::
  velocity_moment(const vector<three_vector> & layer_wind_velocities_meters_per_sec,
		  double moment, 
		  double azimuth_angle_degrees,
		  double zenith_angle_degrees) const {

    try{check_zenith(zenith_angle_degrees);}
    catch(...){
      cerr << "refractive_atmospheric_model::velocity_moment error\n";
      throw(string("refractive_atmospheric_model::turbulence_moment"));
    }

    if(layer_wind_velocities_meters_per_sec.size()!=layer_heights_.size()){
      cerr << "refractive_atmospheric_model::velocity_moment error\n";
      cerr << "number of layers in model " 
	   << layer_heights_.size()
	   << " does not match number of wind velocities "
	   << layer_wind_velocities_meters_per_sec.size()
	   << endl;
      throw(string("refractive_atmospheric_model::turbulence_moment"));
    }      

    double degrees_to_radians = M_PI/180;
    double secant_zenith = 1/cos(zenith_angle_degrees*M_PI/180.);

    Arroyo::three_vector downwind_unit_vector = 
      this->ground_ref_frame_.x()*cos(azimuth_angle_degrees*degrees_to_radians) + 
      this->ground_ref_frame_.y()*sin(azimuth_angle_degrees*degrees_to_radians);
    Arroyo::three_vector crosswind_unit_vector = 
      cross_product(downwind_unit_vector, this->ground_ref_frame_.z());

    double val = 0;
    Arroyo::three_vector wind_vector;
    for(int i=0; i<layer_wind_velocities_meters_per_sec.size(); i++){
      if(dot_product(layer_wind_velocities_meters_per_sec[i],this->ground_ref_frame_.z())>three_frame::precision){
	cerr << "refractive_atmospheric_model::velocity_moment error - \n"
	     << "wind velocity vector " << i 
	     << " is not parallel to the ground\n";
	cerr << "dot product " 
	     << dot_product(layer_wind_velocities_meters_per_sec[i],this->ground_ref_frame_.z())
	     << endl;
	layer_wind_velocities_meters_per_sec[i].print(cerr, "layer velocity vector ");
	this->ground_ref_frame_.print(cerr, "ground reference frame ");
	throw(string("refractive_atmospheric_model::velocity_moment"));
      }
      
    
      wind_vector = 
	fabs(dot_product(layer_wind_velocities_meters_per_sec[i],downwind_unit_vector))*
	cos(zenith_angle_degrees*degrees_to_radians)*
	downwind_unit_vector +
	fabs(dot_product(layer_wind_velocities_meters_per_sec[i],crosswind_unit_vector))*
	crosswind_unit_vector;

      val += get_cn2_coefficient(power_spectra_[i]->get_coefficient())*pow(wind_vector.length(), moment)*secant_zenith;

    }

    return(val);
  }

  namespace {
    double isoangle_coeff(){
      double HJ1 = Arroyo::gamma_function(-5/6.0)/pow(2,8/3.0)/Arroyo::gamma_function(11/6.0);
      double C_A = pow(2,1/3.0)*5/36.0/pow(M_PI,5/3.0)/Arroyo::gamma_function(1/3.0);
      double c_isoangle = 2*pow(2*M_PI,8/3.0)*C_A*fabs(HJ1);
      return(c_isoangle);
    }
  }

  double refractive_atmospheric_model::
  fried_parameter(double wavelength_meters, 
		  double zenith_angle_degrees,
		  double guide_star_height_meters) const {

    try{check_wavelength_zenith(wavelength_meters,zenith_angle_degrees);}
    catch(...){
      cerr << "refractive_atmospheric_model::fried_parameter error\n";
      throw(string("refractive_atmospheric_model::fried_parameter"));
    }

    double val;
    double wavenumber = 2*M_PI/wavelength_meters;
    double fac = isoangle_coeff()/(2*pow(24/5.*gamma_function(6/5.),5/6.));

    if(guide_star_height_meters==-1){
      try{
	// Tyler JOSA 11 p 409 equations A18, A27 and A28
	val=pow(fac*wavenumber*wavenumber*turbulence_moment(0,zenith_angle_degrees), -3/5.0);
      } catch(...){
	cerr << "refractive_atmospheric_model::fried_parameter error\n";
	throw(string("refractive_atmospheric_model::fried_parameter"));
      }
    } else {
      try{check_zenith(zenith_angle_degrees);}
      catch(...){
	cerr << "refractive_atmospheric_model::fried_parameter error\n";
	throw(string("refractive_atmospheric_model::fried_parameter"));
      }

      if(guide_star_height_meters<=this->layer_heights_[0]){
	cerr << "refractive_atmospheric_model::fried_parameter error - guide star height invalid\n";
	throw(string("refractive_atmospheric_model::fried_parameter"));
      }


      std::vector<double>::const_iterator lit;
      std::vector<power_spectrum *>::const_iterator pit;
      std::vector<double>::const_iterator endlit = 
	std::lower_bound(this->layer_heights_.begin(), 
			 this->layer_heights_.end(), 
			 guide_star_height_meters);

      double secant_zenith = 1/cos(zenith_angle_degrees*M_PI/180.);
      double modified_turbulence_moment = 0;
      for(lit = this->layer_heights_.begin(),
	    pit = this->power_spectra_.begin();
	  lit!=endlit;
	  lit++,pit++){
	modified_turbulence_moment += 
	  get_cn2_coefficient((*pit)->get_coefficient()) *
	  pow((guide_star_height_meters - (*lit))/guide_star_height_meters, 5./3.) *
	  secant_zenith;
      }	
      val=pow(fac*wavenumber*wavenumber*modified_turbulence_moment, -3./5.0);
    }

    return(val);
  }

  double refractive_atmospheric_model::
  fried_parameter_downward(double wavelength_meters, 
			   double nadir_angle_degrees,
			   double altitude_meters) const {

    try{check_wavelength_zenith(wavelength_meters,nadir_angle_degrees);}
    catch(...){
      cerr << "refractive_atmospheric_model::fried_parameter error\n";
      throw(string("refractive_atmospheric_model::fried_parameter"));
    }

    if(altitude_meters<=this->layer_heights_[0]){
      cerr << "refractive_atmospheric_model::fried_parameter error - altitude invalid\n";
      throw(string("refractive_atmospheric_model::fried_parameter"));
    }

    double val;
    double wavenumber = 2*M_PI/wavelength_meters;
    double fac = isoangle_coeff()/(2*pow(24/5.*gamma_function(6/5.),5/6.));

    std::vector<double>::const_iterator lit;
    std::vector<power_spectrum *>::const_iterator pit;
    std::vector<double>::const_iterator endlit = 
      std::lower_bound(this->layer_heights_.begin(), 
		       this->layer_heights_.end(), 
		       altitude_meters);
    
    double secant_nadir = 1/cos(nadir_angle_degrees*M_PI/180.);
    double modified_turbulence_moment = 0;
    for(lit = this->layer_heights_.begin(),
	  pit = this->power_spectra_.begin();
	lit!=endlit;
	lit++,pit++){

      modified_turbulence_moment += 
	get_cn2_coefficient((*pit)->get_coefficient()) *
	pow(*lit/altitude_meters, 5/3.) * secant_nadir;
    }	
    
    val=pow(fac*wavenumber*wavenumber*modified_turbulence_moment, -3/5.0);
    return(val);
  }

  double refractive_atmospheric_model::
  two_axis_tilt_jitter(double wavelength_meters, 
		       double aperture_diameter_meters,
		       double zenith_angle_degrees,
		       double guide_star_height_meters) const {
    
    try{check_wavelength_zenith(wavelength_meters,zenith_angle_degrees);}
    catch(...){
      cerr << "refractive_atmospheric_model::two_axis_tilt_jitter error\n";
      throw(string("refractive_atmospheric_model::two_axis_tilt_jitter"));
    }

    double val;

    if(guide_star_height_meters==-1){
      try{
	// Sasiela eq 4.15
	double r0 = this->fried_parameter(wavelength_meters, zenith_angle_degrees);
	val=.3641*pow((aperture_diameter_meters) / r0, 5./3.)*pow((wavelength_meters / aperture_diameter_meters),2.);
	val = sqrt(val);
      } catch(...){
	cerr << "refractive_atmospheric_model::two_axis_tilt_jitter error\n";
	throw(string("refractive_atmospheric_model::two_axis_tilt_jitter"));
      }
    } else {
      try{check_zenith(zenith_angle_degrees);}
      catch(...){
	cerr << "refractive_atmospheric_model::two_axis_tilt_jitter error\n";
	throw(string("refractive_atmospheric_model::two_axis_tilt_jitter"));
      }

      if(guide_star_height_meters<=this->layer_heights_[0]){
	cerr << "refractive_atmospheric_model::two_axis_tilt_jitter error - guide star height invalid\n";
	throw(string("refractive_atmospheric_model::two_axis_tilt_jitter"));
      }


      std::vector<double>::const_iterator lit;
      std::vector<power_spectrum *>::const_iterator pit;
      std::vector<double>::const_iterator endlit = 
	std::lower_bound(this->layer_heights_.begin(), 
			 this->layer_heights_.end(), 
			 guide_star_height_meters);

      double secant_zenith = 1/cos(zenith_angle_degrees*M_PI/180.);
      double modified_turbulence_moment = 0;
      for(lit = this->layer_heights_.begin(),
	    pit = this->power_spectra_.begin();
	  lit!=endlit;
	  lit++,pit++){
	modified_turbulence_moment += 
	  get_cn2_coefficient((*pit)->get_coefficient()) *
	  pow((guide_star_height_meters - (*lit))/guide_star_height_meters, 5./3.) *
	  secant_zenith;
      }	

      // Sasiela between eq 4.14 and 4.15
      val=105.1*gamma_function(1./6.)*gamma_function(7./3.)*modified_turbulence_moment / 
	gamma_function(29./6.) / gamma_function(17./6.) / 2 / sqrt(M_PI) / pow(aperture_diameter_meters,1./3.);
      val = sqrt(val);
    }

    return(val);
  }

  double refractive_atmospheric_model::
  two_axis_tilt_jitter_downward(double wavelength_meters, 
				double aperture_diameter_meters,
				double nadir_angle_degrees,
				double altitude_meters) const {

    try{check_wavelength_zenith(wavelength_meters,nadir_angle_degrees);}
    catch(...){
      cerr << "refractive_atmospheric_model::two_axis_tilt_jitter error\n";
      throw(string("refractive_atmospheric_model::two_axis_tilt_jitter"));
    }

    if(altitude_meters<=this->layer_heights_[0]){
      cerr << "refractive_atmospheric_model::two_axis_tilt_jitter error - altitude invalid\n";
      throw(string("refractive_atmospheric_model::two_axis_tilt_jitter"));
    }

    double val;

    std::vector<double>::const_iterator lit;
    std::vector<power_spectrum *>::const_iterator pit;
    std::vector<double>::const_iterator endlit = 
      std::lower_bound(this->layer_heights_.begin(), 
		       this->layer_heights_.end(), 
		       altitude_meters);
    
    double secant_nadir = 1/cos(nadir_angle_degrees*M_PI/180.);
    double modified_turbulence_moment = 0;
    for(lit = this->layer_heights_.begin(),
	  pit = this->power_spectra_.begin();
	lit!=endlit;
	lit++,pit++){

      modified_turbulence_moment += 
	get_cn2_coefficient((*pit)->get_coefficient()) *
	pow(*lit/altitude_meters, 5/3.) * secant_nadir;
    }	
     
    // Sasiela between eq 4.14 and 4.15
    val=105.1*gamma_function(1./6.)*gamma_function(7./3.)*modified_turbulence_moment / 
      gamma_function(29./6.) / gamma_function(17./6.) / 2 / sqrt(M_PI) / pow(aperture_diameter_meters,1./3.);
    val = sqrt(val);

    return(val);
  }
 
  double refractive_atmospheric_model::
  seeing(double wavelength_meters, 
	 double zenith_angle_degrees) const {

    try{check_wavelength_zenith(wavelength_meters,zenith_angle_degrees);}
    catch(...){
      cerr << "refractive_atmospheric_model::seeing error\n";
      throw(string("refractive_atmospheric_model::seeing"));
    }

    double val;
    double wavenumber = 2*M_PI/wavelength_meters;
    try{
      val=wavelength_meters/this->fried_parameter(wavelength_meters, zenith_angle_degrees);
    }
    catch(...){
      cerr << "refractive_atmospheric_model::seeing error\n";
      throw(string("refractive_atmospheric_model::seeing"));
    }
    return(val);
  }

  double refractive_atmospheric_model::
  isoplanatic_angle(double wavelength_meters, 
		    double zenith_angle_degrees,
		    double guide_star_height_meters) const {
    try{check_wavelength_zenith(wavelength_meters,zenith_angle_degrees);}
    catch(...){
      cerr << "refractive_atmospheric_model::isoplanatic_angle error\n";
      throw(string("refractive_atmospheric_model::isoplanatic_angle"));
    }

    double val;
    double wavenumber = 2*M_PI/wavelength_meters;

    if(guide_star_height_meters==-1){

      try{val=pow(isoangle_coeff()*wavenumber*wavenumber*turbulence_moment(5/3.,zenith_angle_degrees), -3/5.0);}
      catch(...){
	cerr << "refractive_atmospheric_model::isoplanatic_angle error\n";
	throw(string("refractive_atmospheric_model::isoplanatic_angle"));
      }

    } else {
      std::vector<double>::const_iterator lit;
      std::vector<power_spectrum *>::const_iterator pit;
      std::vector<double>::const_iterator endlit = 
	std::lower_bound(this->layer_heights_.begin(), 
			 this->layer_heights_.end(), 
			 guide_star_height_meters);

      double range_meters;
      double secant_zenith = 1/cos(zenith_angle_degrees*M_PI/180.);
      double modified_turbulence_moment = 0;
      for(lit = this->layer_heights_.begin(),
	    pit = this->power_spectra_.begin();
	  lit!=endlit;
	  lit++,pit++){
	range_meters = (*lit)*secant_zenith;
	modified_turbulence_moment += 
	  get_cn2_coefficient((*pit)->get_coefficient()) *
	  pow(range_meters*(guide_star_height_meters - (*lit))/guide_star_height_meters, 5/3.) *
	  secant_zenith;
      }	
      val=pow(isoangle_coeff()*wavenumber*wavenumber*modified_turbulence_moment, -3/5.0);
    }

    return(val);
  }

  double refractive_atmospheric_model::
  isoplanatic_angle_downward(double wavelength_meters, 
			     double nadir_angle_degrees,
			     double altitude_meters) const {
    
    try{check_wavelength_zenith(wavelength_meters,nadir_angle_degrees);}
    catch(...){
      cerr << "refractive_atmospheric_model::isoplanatic_angle_downward error\n";
      throw(string("refractive_atmospheric_model::isoplanatic_angle_downward"));
    }

    if(altitude_meters<=this->layer_heights_[0]){
      cerr << "refractive_atmospheric_model::isoplanatic_angle_downward error - altitude invalid\n";
      throw(string("refractive_atmospheric_model::isoplanatic_angle_downward"));
    }

    double val;
    double wavenumber = 2*M_PI/wavelength_meters;

    std::vector<double>::const_iterator lit;
    std::vector<power_spectrum *>::const_iterator pit;
    std::vector<double>::const_iterator endlit = 
      std::lower_bound(this->layer_heights_.begin(), 
		       this->layer_heights_.end(), 
		       altitude_meters);
    
    double range_meters;
    double secant_nadir = 1/cos(nadir_angle_degrees*M_PI/180.);
    double modified_turbulence_moment = 0;
    for(lit = this->layer_heights_.begin(),
	  pit = this->power_spectra_.begin();
	lit!=endlit;
	lit++,pit++){
      range_meters = (*lit)*secant_nadir;
      modified_turbulence_moment += 
	get_cn2_coefficient((*pit)->get_coefficient()) *
	pow(range_meters*(*lit)/altitude_meters, 5/3.) * secant_nadir;
    }	
    
    val=pow(isoangle_coeff()*wavenumber*wavenumber*modified_turbulence_moment, -3/5.0);
    return(val);
  }
   
  double refractive_atmospheric_model::
  isokinetic_angle(double wavelength_meters,
		   double aperture_diameter_meters,
		   double zenith_angle_degrees,
		   double guide_star_height_meters) const {

    try{check_wavelength_zenith(wavelength_meters,zenith_angle_degrees);}
    catch(...){
      cerr << "refractive_atmospheric_model::isokinetic_angle error\n";
      throw(string("refractive_atmospheric_model::isokinetic_angle"));
    }

    if(aperture_diameter_meters<=0){
      cerr << "refractive_atmospheric_model::isokinetic_angle error - aperture diameter "
	   << aperture_diameter_meters
	   << " out of range\n";
      throw(string("refractive_atmospheric_model::isokinetic_angle"));
    }

    double val;
    double wavenumber = 2*M_PI/wavelength_meters;

    if(guide_star_height_meters==-1){
      try{
	
	/*
	// Hardy equation 7.62
	val=pow(.668*wavenumber*wavenumber*
		turbulence_moment(2.,zenith_angle_degrees)*
		pow(aperture_diameter_meters,1/3.), -.5);
	*/

	// Sasiela equation 7.74
	val = .184*wavelength_meters*pow(aperture_diameter_meters,1/6.) / 
	  sqrt(turbulence_moment(2.,zenith_angle_degrees));
      } catch(...){
	cerr << "refractive_atmospheric_model::isokinetic_angle error\n";
	throw(string("refractive_atmospheric_model::isokinetic_angle"));
      }
    } else {
      std::vector<double>::const_iterator lit;
      std::vector<power_spectrum *>::const_iterator pit;
      std::vector<double>::const_iterator endlit = 
	std::lower_bound(this->layer_heights_.begin(), 
			 this->layer_heights_.end(), 
			 guide_star_height_meters);

      double range_meters;
      double secant_zenith = 1/cos(zenith_angle_degrees*M_PI/180.);
      double modified_turbulence_moment = 0;
      for(lit = this->layer_heights_.begin(),
	    pit = this->power_spectra_.begin();
	  lit!=endlit;
	  lit++,pit++){
	range_meters = (*lit)*secant_zenith;
	modified_turbulence_moment += 
	  get_cn2_coefficient((*pit)->get_coefficient()) *
	  pow((guide_star_height_meters - (*lit))/guide_star_height_meters, 5/3.) * range_meters * range_meters *
	  secant_zenith;
      }	
      val = .184*wavelength_meters*pow(aperture_diameter_meters,1/6.) / 
	sqrt(modified_turbulence_moment);
    }

    return(val);
  }

  double refractive_atmospheric_model::
  isokinetic_angle_downward(double wavelength_meters, 
			    double aperture_diameter_meters,
			    double nadir_angle_degrees,
			    double altitude_meters) const {
    
    try{check_wavelength_zenith(wavelength_meters,nadir_angle_degrees);}
    catch(...){
      cerr << "refractive_atmospheric_model::isokinetic_angle_downward error\n";
      throw(string("refractive_atmospheric_model::isokinetic_angle_downward"));
    }

    if(altitude_meters<=this->layer_heights_[0]){
      cerr << "refractive_atmospheric_model::isokinetic_angle_downward error - altitude invalid\n";
      throw(string("refractive_atmospheric_model::isokinetic_angle_downward"));
    }

    double val;
    double wavenumber = 2*M_PI/wavelength_meters;

    std::vector<double>::const_iterator lit;
    std::vector<power_spectrum *>::const_iterator pit;
    std::vector<double>::const_iterator endlit = 
      std::lower_bound(this->layer_heights_.begin(), 
		       this->layer_heights_.end(), 
		       altitude_meters);

    double range_meters;
    double secant_nadir = 1/cos(nadir_angle_degrees*M_PI/180.);
    double modified_turbulence_moment = 0;
    for(lit = this->layer_heights_.begin(),
	  pit = this->power_spectra_.begin();
	lit!=endlit;
	lit++,pit++){
      range_meters = (*lit)*secant_nadir;
      modified_turbulence_moment += 
	get_cn2_coefficient((*pit)->get_coefficient()) *
	pow((*lit)/altitude_meters, 5/3.) * range_meters * range_meters * secant_nadir;
    }	
    
    val = .184*wavelength_meters*pow(aperture_diameter_meters,1/6.) / 
      sqrt(modified_turbulence_moment);

    return(val);
  }
 
  double refractive_atmospheric_model::
  greenwood_frequency(vector<three_vector> & layer_wind_velocities_meters_per_sec,
		      double wavelength_meters, 
		      double azimuth_angle_degrees,
		      double zenith_angle_degrees) const {

    try{check_wavelength_zenith(wavelength_meters,zenith_angle_degrees);}
    catch(...){
      cerr << "refractive_atmospheric_model::greenwood_frequency error\n";
      throw(string("refractive_atmospheric_model::greenwood_frequency"));
    }
    

    double wavenumber = 2*M_PI/wavelength_meters;
    double c_greenwood = pow(2,1/3.0)*gamma_function(1/6.0)/18.0/pow(M_PI,7/6.0);
    double v_five_over_three = this->velocity_moment(layer_wind_velocities_meters_per_sec,
						     5/3.,
						     azimuth_angle_degrees,
						     zenith_angle_degrees);
    double greenwood_frequency = 
      pow(c_greenwood*wavenumber*wavenumber*v_five_over_three,3/5.0);

    return(greenwood_frequency);
  }
  
  double refractive_atmospheric_model::
  d_0(double guide_star_height_meters,
      double wavelength_meters,
      double zenith_angle_degrees) const {
    
    try{check_wavelength_zenith(wavelength_meters,zenith_angle_degrees);}
    catch(...){
      cerr << "refractive_atmospheric_model::d_0 error\n";
      throw(string("refractive_atmospheric_model::d_0"));
    }

    if(guide_star_height_meters<=0){
      cerr << "refractive_atmospheric_model::d_0 error\n"
	   << "guide star height " 
	   << guide_star_height_meters
	   << " out of range\n";
      throw(string("refractive_atmospheric_model::d_0"));
    }
 
    double val = 0;
    double wavenumber = 2*M_PI/wavelength_meters;
    double secant_zenith = 1/cos(zenith_angle_degrees*M_PI/180.);
    int npowerspec = this->power_spectra_.size();
    for(int i=0; i<npowerspec; i++){
      if(layer_heights_[i]>=guide_star_height_meters)
	val += wavenumber*wavenumber*.057*get_cn2_coefficient(power_spectra_[i]->get_coefficient());
      else
	val += wavenumber*wavenumber*
	  get_cn2_coefficient(power_spectra_[i]->get_coefficient())*
	  (.5*pow(layer_heights_[i]/guide_star_height_meters,5/3.)*pow(secant_zenith,8/3.)
	   - .452*pow(layer_heights_[i]/guide_star_height_meters,2.)*pow(secant_zenith,3));
    }

    val = pow(val, -3/5.);

    return(val);
  }

  double refractive_atmospheric_model::
  scintillation_quenching_diameter(double wavelength_meters,
				   double zenith_angle_degrees) const {
    
    try{check_wavelength_zenith(wavelength_meters,zenith_angle_degrees);}
    catch(...){
      cerr << "refractive_atmospheric_model::scintillation_quenching_diameter error\n";
      throw(string("refractive_atmospheric_model::scintillation_quenching_diameter"));
    }

    double val;
    double wavenumber = 2*M_PI/wavelength_meters;
    // Sasiela equation 7.96
    try{
      val=.957*
	pow(turbulence_moment(2.,zenith_angle_degrees)/turbulence_moment(5/6.,zenith_angle_degrees),3/7.)*
	sqrt(wavelength_meters);
    } catch(...){
      cerr << "refractive_atmospheric_model::scintillation_quenching_diameter error\n";
      throw(string("refractive_atmospheric_model::scintillation_quenching_diameter"));
    }

    return(val);
  }

  double refractive_atmospheric_model::
  scintillation_isoplanatic_angle(double wavelength_meters,
				  double zenith_angle_degrees) const {
    
    try{check_wavelength_zenith(wavelength_meters,zenith_angle_degrees);}
    catch(...){
      cerr << "refractive_atmospheric_model::scintillation_isoplanatic_angle error\n";
      throw(string("refractive_atmospheric_model::scintillation_isoplanatic_angle"));
    }

    double val;
    double wavenumber = 2*M_PI/wavelength_meters;
    // Sasiela equation 7.103
    try{
      val=.0957*
	pow(turbulence_moment(-1/3.,zenith_angle_degrees)/turbulence_moment(5/6.,zenith_angle_degrees),3/7.)*
	sqrt(wavelength_meters);
    } catch(...){
      cerr << "refractive_atmospheric_model::scintillation_isoplanatic_angle error\n";
      throw(string("refractive_atmospheric_model::scintillation_isoplanatic_angle"));
    }

    return(val);
  }

  double refractive_atmospheric_model::
  log_amplitude_variance(double wavelength_meters,
			 double zenith_angle_degrees,
			 double altitude_meters) const {


    double fac = .5631 * pow(2*M_PI / wavelength_meters, 7/6.);
    double val = 0;

    std::vector<double>::const_iterator lit;
    std::vector<power_spectrum *>::const_iterator pit;
    std::vector<double>::const_iterator endlit;

    if(altitude_meters == -1){
      endlit = this->layer_heights_.end();
      altitude_meters =  *std::max_element(this->layer_heights_.begin(), 
				  this->layer_heights_.end());
    } else {
      endlit = std::lower_bound(this->layer_heights_.begin(), 
				this->layer_heights_.end(), 
				altitude_meters);
    }    

    double secant_zenith = 1/cos(zenith_angle_degrees*M_PI/180.);
    double modified_turbulence_moment = 0;
    for(lit = this->layer_heights_.begin(),
	  pit = this->power_spectra_.begin();
	lit!=endlit;
	lit++,pit++){
      val += 
	get_cn2_coefficient((*pit)->get_coefficient()) *
	pow((altitude_meters - (*lit))*(*lit)/altitude_meters, 5/6.);
    }	

    return(val*fac * pow(secant_zenith,11/6.));

  };

  double refractive_atmospheric_model::
  log_amplitude_variance_downward(double wavelength_meters,
				  double nadir_angle_degrees,
				  double altitude_meters) const {

    try{check_wavelength_zenith(wavelength_meters,nadir_angle_degrees);}
    catch(...){
      cerr << "refractive_atmospheric_model::log_amplitude_variance_downward error\n";
      throw(string("refractive_atmospheric_model::log_amplitude_variance_downward"));
    }

    double fac = .5631 * pow(2*M_PI / wavelength_meters, 7/6.);
    double val = 0;

    std::vector<double>::const_iterator lit;
    std::vector<power_spectrum *>::const_iterator pit;
    std::vector<double>::const_iterator endlit = 
      std::lower_bound(this->layer_heights_.begin(), 
		       this->layer_heights_.end(), 
		       altitude_meters);
    
    double secant_nadir = 1/cos(nadir_angle_degrees*M_PI/180.);
    double modified_turbulence_moment = 0;
    for(lit = this->layer_heights_.begin(),
	  pit = this->power_spectra_.begin();
	lit!=endlit;
	lit++,pit++){
      val += 
	get_cn2_coefficient((*pit)->get_coefficient()) *
	pow((altitude_meters - (*lit))*(*lit)/altitude_meters, 5/6.);
    }	
    
    return(val*fac * pow(secant_nadir,11/6.));
  };


  // Tyler JOSA A v 11 p 343 1994
  // Equation 34
  double refractive_atmospheric_model::Tyler_F_1(double x) const {
    
    if(x<0){
      cerr << "F_1 error - invalid argument "
	   << x << endl;
      throw(string("F_1"));
    }
      
    double val;
    vector<double> numerator_vals(2,-5/6.);
    vector<double> denominator_vals(1,2);
      
    try{
      if(x<1){
	numerator_vals[0] = -11/6.;
	denominator_vals[0] = 1;
	val = (6/11.)*generalized_hypergeometric_function(x*x,
							  numerator_vals,
							  denominator_vals);
      } else {
	val = pow(x,5/3.0)*generalized_hypergeometric_function(1/x/x,
							       numerator_vals,
							       denominator_vals);
      }
    } catch(...){
      cerr << "F_1 error\n";
      throw(string("F_1"));
    }
    return(val);
  }

  // Tyler JOSA A v 11 p 343 1994
  // Equation 36
  double refractive_atmospheric_model::Tyler_H(double & rho, 
					 double & omega,
					 vector<double> & numerator_args,
					 vector<double> & denominator_args) const {
      
    double h;
    try{
      if(rho>omega)
	h = pow(rho,8/3.)*
	  generalized_hypergeometric_function(omega*omega/rho/rho,
					      numerator_args,
					      denominator_args);
      else 
	h = rho*pow(omega,5/3.)*
	  generalized_hypergeometric_function(rho*rho/omega/omega,
					      numerator_args,
					      denominator_args);
    } catch(...){
      cerr << "H error\n";
      throw(string("H"));
    }
    return(h);
  }

  // Tyler JOSA A v 11 p 343 1994
  // Equations 37-39
  double refractive_atmospheric_model::Tyler_K_1(double rho, double q) const {
      
    if(fabs(rho)<=1-q)
      return(M_PI*q*q);
    else if(fabs(rho)>1-q && 
	    fabs(rho)<1+q){
      double s_o = (-(1 - q*q)/2./rho + rho/2.)/q;
      double t_o = (1 - q*q)/2./rho + rho/2.;
      double val = acos(t_o) - t_o * sqrt(1-t_o*t_o) + q*q*(acos(s_o) - (s_o)*sqrt(1 - s_o*s_o));

      return(val);
    } else 
      return(0);
  }
    
  // Tyler JOSA A v 11 p 343 1994
  // Equation 35    
  double refractive_atmospheric_model::Tyler_F_2(double q, 
					   double omega, 
					   int nsamples_in_integration) const {

    double val = 0;
    double step_size = (1+q) / (double)nsamples_in_integration;
    double rho = .5*step_size;
    vector<double> numerator_args(2,-5/6.);
    vector<double> denominator_args(1,1);

    try{
      while(rho<(1+q)){
	val += Tyler_K_1(rho, q)*Tyler_H(rho, omega, numerator_args, denominator_args)*step_size;
	rho += step_size;
      }    
      val *= 2/M_PI/q/q;	
    } catch(...){
      cerr << "F_2 error - could not compute value\n";
      throw(string("F_2"));
    }
    return(val);
  }

  // Tyler JOSA A v 11 p 409 1994
  // Equation A81 defines this quantity.
  double refractive_atmospheric_model::Tyler_G_hat(double rho) const {

    if(rho<0){
      cerr << "G_hat error - argument "
	   << rho 
	   << " is not positive definite\n";
      throw(string("G_hat"));
    }

    vector<double> hypergeometric_vals(2,1/6.0);
    fabs(rho)<=1 ? 
      hypergeometric_vals[0] = -11/6. :
      hypergeometric_vals[0] = -5/6.;
	

    try{
      double val = rho<=1 ? 
	(-5/11.0)*
	generalized_hypergeometric_function(rho*rho,
					    hypergeometric_vals,
					    vector<double>(1,2))
	:
	(-5/12.0)*pow(rho, -1/3.)*
	generalized_hypergeometric_function(1/rho/rho,	
					    hypergeometric_vals,
					    vector<double>(1,3));
      if(!finite(val)){
	cerr << "G_hat error - val not finite\n";
	cerr << "rho " << rho << endl;
	throw(string("G_hat"));
      }
      return(val);

    } catch(...) {
      cerr << "G_hat failure for argument "
	   << rho 
	   << endl;
      throw(string("G_hat"));
    }
  }
 

  // Tyler JOSA A v 11 p 409 1994
  // Equation A75 
  // (Note:  typo in equation 31)
  double refractive_atmospheric_model::Tyler_G_1(three_vector r1, 
					   three_vector r2) const {
    double val = dot_product(r1,r2)*Tyler_G_hat(r2.length());
    return(val);
  }


  // Tyler JOSA A v 11 p 409 1994
  // Equations 32-35 are inconsistent
  // with equations A81 and A94
  double refractive_atmospheric_model::Tyler_G_2(three_vector r, 
					   double q, 
					   three_vector omega, 
					   int nsamples_in_integration) const {

    double omega_magnitude = omega.length();
    if(omega_magnitude<three_frame::precision) 
      return(0);

    double omega_squared = omega.length_squared();

    double rho = fabs(q - omega_magnitude);
    double val = 0;
    double step_size = (q + omega_magnitude - rho)/(double)nsamples_in_integration;

    do {
      val += step_size *
	rho*rho *
	sqrt(1-(omega_squared + rho*rho - q*q)/2.0/omega_magnitude/rho)*
	Tyler_G_hat(rho);
      rho += step_size;
    } while(rho<(q+omega_magnitude));
    val *= 2*dot_product(omega, r)/M_PI/q/q/omega_magnitude;
      
    if(!finite(val)){
      cerr << "G_2 error - val not finite\n";
      cerr << endl;
      cerr << "q " << q << endl;
      cerr << endl;
      r.print(cerr, "r ");
      cerr << endl;
      omega.print(cerr, "omega ");
      throw(string("G_2"));
    }
    return(val);
  }


  // Tyler JOSA A v 11 p 343 1994
  // Equations 40-41
  double refractive_atmospheric_model::Tyler_G_3(double q, 
					   double omega,
					   int nsamples_in_integration) const {

    double phi_bar = 0;
    if(omega>=(1-q) && omega<=(1+q))
      phi_bar = acos((1+q*q-omega*omega)/(2*q));
    else if(omega>=(1+q))
      phi_bar = M_PI;
      
    double step_size = M_PI/(double)nsamples_in_integration;
    double phi = .5*step_size;
    double val1 = 0;
    double val2 = 0;
      
    vector<double> numerator_vals(2,1/6.0);
    vector<double> denominator_vals(1,3);

    if(phi_bar>phi){
      do {
	try{
	  val1 += step_size * 
	    pow(sin(phi),4.0) *
	    generalized_hypergeometric_function((1+q*q-2*q*cos(phi))/omega/omega,
						numerator_vals,
						denominator_vals);
	} catch(...) {
	  cerr << "G_3 error in first integral\n";
	  throw(string("G_3"));
	}
	  
	phi += step_size;
      } while(phi < phi_bar);
	
      val1 *= pow(omega, -1/3.0);
    }
      
    numerator_vals[0] = 1/6.0;
    numerator_vals[1] = -11/6.0;
    denominator_vals[0] = 1;
      
    if(phi<M_PI){
      do {
	  
	try{
	  val2 += step_size * 
	    pow(sin(phi),4.0) *
	    pow((1+q*q-2*q*cos(phi)),-1/6.0) *
	    generalized_hypergeometric_function(omega*omega/(1+q*q-2*q*cos(phi)),
						numerator_vals,
						denominator_vals);
	} catch(...) {
	  cerr << "G_3 error in second integral\n";
	  throw(string("G_3"));
	}
	  
	phi += step_size;
      } while(phi < M_PI);
	
      val2 *= 72/55.0;
    }
    return(-(val1+val2)*q*25/54.0/M_PI);
  }

  // Tyler JOSA A v 11 p 409 1994
  // Equation A101
  double refractive_atmospheric_model::Tyler_G_4(three_vector r1, 
					   three_vector r2, 
					   double q, 
					   three_vector omega, 
					   int nsamples_in_integration) const {

    double rho;
    double step_size;
    double omega_magnitude = omega.length();
    double omega_squared = omega.length_squared();

    double val1=0, val2=0, val3=0;


    // First term
    if(q-omega_magnitude>0){
      rho = 0;
      step_size = (q-omega_magnitude)/(double)nsamples_in_integration;
      do {
	val1 += step_size *
	  rho*rho*rho* Tyler_G_hat(rho);
	rho += step_size;
      } while(rho<(q-omega_magnitude));
      val1 *= dot_product(r1,r2)/q/q/q;
    }

    if(!finite(val1)){
      cerr << "G_4 error - first term is not finite\n";
      cerr << endl;
      cerr << "q " << q << endl;
      cerr << endl;
      r1.print(cerr, "r1 ");
      cerr << endl;
      r2.print(cerr, "r2 ");
      cerr << endl;
      omega.print(cerr, "omega ");
      throw(string("G_4"));
    }

    if(omega_magnitude<three_frame::precision)
      return(val1);

    // second and third terms
    double fac_a, fac_b;
    rho = fabs(q-omega_magnitude);
    step_size = (q + omega.length() - rho)/(double)nsamples_in_integration;
    do {
      fac_a = step_size*rho*rho*rho*Tyler_G_hat(rho);
      fac_b = (omega_squared + rho*rho-q*q)/2.0/omega_magnitude/rho;
      val2 += fac_a * (1 - fabs(fac_b))*sqrt(fabs(1-fac_b*fac_b));
      val3 += fac_a * (fabs(fac_b) - omega_magnitude/rho)*sqrt(1-fac_b*fac_b);
    } while(rho<(q+omega_magnitude));

    val2 *= dot_product(r1,r2)/M_PI/q/q/q;
    val3 *= 2*dot_product(r1,omega)*dot_product(r2,omega)/M_PI/q/q/q/omega_squared;

    if(!finite(val2) || 
       !finite(val3)){
      cerr << "G_4 error - one of the last two terms is not finite\n";
      cerr << "\tval 2 " << val2
	   << "\tval 3 " << val3
	   << endl;
      cerr << endl;
      cerr << "q " << q << endl;
      cerr << endl;
      r1.print(cerr, "r1 ");
      cerr << endl;
      r2.print(cerr, "r2 ");
      cerr << endl;
      omega.print(cerr, "omega ");
      throw(string("G_4"));
    }
    return(val1+val2+val3);
  }

  void refractive_atmospheric_model::Tyler_get_constants(const emitter & emtr_a,
							 const emitter & emtr_b,
							 const three_frame & tf,
							 double aperture_diameter_meters,
							 double & secant_zenith_angle,
							 double & max_range_meters,
							 double & min_range_meters,
							 three_vector & little_omega) const {

    // Set up beams...
    three_vector beam_a_unit_vector = 
      emtr_a.get_emission_vector(static_cast<const three_point>(tf));
    three_vector beam_b_unit_vector = 
      emtr_b.get_emission_vector(static_cast<const three_point>(tf));
      
    three_vector mean_pointing_unit_vector = 
      beam_a_unit_vector + beam_b_unit_vector;
    mean_pointing_unit_vector *= (1/mean_pointing_unit_vector.length());
      
    // three frame is pointing up, so switch sign on z axis to get angle
    // w.r.t mean pointing vector
    secant_zenith_angle = 
      1/dot_product(-1*tf.z(),mean_pointing_unit_vector);

    if(!finite(secant_zenith_angle)){
      cerr << "get_constants error - secant zenith angle "
	   << secant_zenith_angle
	   << " not finite\n";
      beam_a_unit_vector.print(cerr, "beam a vector ");
      beam_b_unit_vector.print(cerr, "beam b vector ");
      mean_pointing_unit_vector.print(cerr, "pointing vector ");
      throw(string("get_constants"));
    }
    /*
      double beam_a_range_meters = 1e20;
      double beam_b_range_meters = 1e20;
    */

    double beam_a_range_meters = 1e300;
    double beam_b_range_meters = 1e300;

    const spherical_wave_emitter * swe;
    if((swe=dynamic_cast<const spherical_wave_emitter *>(&emtr_a)))
      beam_a_range_meters = (*swe - tf).length();
    if(swe=dynamic_cast<const spherical_wave_emitter *>(&emtr_b))
      beam_b_range_meters = (*swe - tf).length();

    three_vector max_range_unit_vector, min_range_unit_vector;
    if(beam_a_range_meters>beam_b_range_meters){
      max_range_meters = beam_a_range_meters;
      min_range_meters = beam_b_range_meters;
      max_range_unit_vector = beam_a_unit_vector;
      min_range_unit_vector = beam_b_unit_vector;
    } else {
      max_range_meters = beam_b_range_meters;
      min_range_meters = beam_a_range_meters;
      max_range_unit_vector = beam_b_unit_vector;
      min_range_unit_vector = beam_a_unit_vector;
    }

    if((max_range_unit_vector - min_range_unit_vector).length()<three_frame::precision)
      little_omega = three_vector();
    else 
      little_omega = (max_range_unit_vector - min_range_unit_vector)*(2/aperture_diameter_meters);
  }
  
  vector<double> refractive_atmospheric_model::get_cn2_coefficients() const {

    vector<double> cn2_coefficients;
    for(int k=0; k<this->power_spectra_.size(); k++)
      cn2_coefficients.push_back(get_cn2_coefficient(power_spectra_[k]->get_coefficient()));
    return(cn2_coefficients);
  }

  double refractive_atmospheric_model::caliph_A(const emitter & emtr_a,
						const emitter & emtr_b,
						double aperture_diameter_meters,
						const three_vector & rho) const {

    if(rho.length_squared()>1){
      cerr << "refractive_atmospheric_model::caliph_A error - \n"
	   << "three vector passed to this function has an amplitude "
	   << rho.length()
	   << " greater than unity\n";
      rho.print(cerr, "rho ");
      throw(string("refractive_atmospheric_model::caliph_A"));
    }

    if(aperture_diameter_meters<=0){
      cerr << "refractive_atmospheric_model::caliph_A error - \n"
	   << "invalid aperture diameter "
	   << aperture_diameter_meters
	   << " passed to this function\n";
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

      double val = 0;

      double Q;
      int npowerspec = this->power_spectra_.size();
      for(int i=0; i<npowerspec; i++){
	if(layer_heights_[i]>min_range_meters) continue;

	Q = (1 - layer_heights_[i]*secant_zenith_angle/min_range_meters) /
	  (1 - layer_heights_[i]*secant_zenith_angle/max_range_meters);

	val += secant_zenith_angle*get_cn2_coefficient(power_spectra_[i]->get_coefficient()) * 
	  pow((1 - layer_heights_[i]*secant_zenith_angle/max_range_meters),5/3.) *
	  pow(Q,5/3.) *
	  Tyler_F_1((rho + 
	       (little_omega*(layer_heights_[i]/(1-layer_heights_[i]/max_range_meters)))).length()
	      /Q);
      }
      if(!finite(val)){
	cerr << "refractive_atmospheric_model::caliph_A error - value " 
	     << val 
	     << " not finite\n";
	throw(string("refractive_atmospheric_model::caliph_A"));
      }

      return(val);
    } catch(...) {
      cerr << "refractive_atmospheric_model::caliph_A error\n";
      throw(string("refractive_atmospheric_model::caliph_A"));
    }
  }

  double refractive_atmospheric_model::caliph_B(const emitter & emtr_a,
						const emitter & emtr_b,
						double aperture_diameter_meters,
						const three_vector & rho) const {

    if(rho.length_squared()>1){
      cerr << "refractive_atmospheric_model::caliph_B error - \n"
	   << "three vector passed to this function has an amplitude "
	   << rho.length()
	   << " greater than unity\n";
      rho.print(cerr, "rho ");
      throw(string("refractive_atmospheric_model::caliph_B"));
    }

    if(aperture_diameter_meters<=0){
      cerr << "refractive_atmospheric_model::caliph_B error - \n"
	   << "invalid aperture diameter "
	   << aperture_diameter_meters
	   << " passed to this function\n";
      throw(string("refractive_atmospheric_model::caliph_B"));
    }

    try {
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

      double val = 0;

      double Q;
      int npowerspec = this->power_spectra_.size();
      for(int i=0; i<npowerspec; i++){
	if(layer_heights_[i]>min_range_meters) continue;

	Q = (1 - layer_heights_[i]*secant_zenith_angle/min_range_meters) /
	  (1 - layer_heights_[i]*secant_zenith_angle/max_range_meters);

	val += secant_zenith_angle*get_cn2_coefficient(power_spectra_[i]->get_coefficient()) * 
	  pow((1 - layer_heights_[i]*secant_zenith_angle/max_range_meters),5/3.) *
	  Tyler_F_1(((Q*rho) - 
	       (little_omega*(layer_heights_[i]/(1-layer_heights_[i]/max_range_meters)))).length());
      }
      if(!finite(val)){
	cerr << "refractive_atmospheric_model::caliph_B error - value " 
	     << val 
	     << " not finite\n";
	throw(string("refractive_atmospheric_model::caliph_B"));
      }

      return(val);
    } catch(...) {
      cerr << "refractive_atmospheric_model::caliph_B error\n";
      throw(string("refractive_atmospheric_model::caliph_B"));
    }
  }
  
  double refractive_atmospheric_model::caliph_C(const emitter & emtr_a,
						const emitter & emtr_b,
						double aperture_diameter_meters,
						const three_vector & rho_1,
						const three_vector & rho_2) const {

    if(rho_1.length_squared()>1){
      cerr << "refractive_atmospheric_model::caliph_C error - \n"
	   << "three vector passed to this function has an amplitude "
	   << rho_1.length()
	   << " greater than unity\n";
      rho_1.print(cerr, "rho 1 ");
      throw(string("refractive_atmospheric_model::caliph_C"));
    }
    if(rho_2.length_squared()>1){
      cerr << "refractive_atmospheric_model::caliph_C error - \n"
	   << "three vector passed to this function has an amplitude "
	   << rho_2.length()
	   << " greater than unity\n";
      rho_2.print(cerr, "rho 1 ");
      throw(string("refractive_atmospheric_model::caliph_C"));
    }

    if(aperture_diameter_meters<=0){
      cerr << "refractive_atmospheric_model::caliph_C error - \n"
	   << "invalid aperture diameter "
	   << aperture_diameter_meters
	   << " passed to this function\n";
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

      double val = 0;

      double Q;
      double fac;
      int npowerspec = this->power_spectra_.size();
      for(int i=0; i<npowerspec; i++){
	if(layer_heights_[i]>min_range_meters) continue;

	fac = (1 - layer_heights_[i]*secant_zenith_angle/max_range_meters);

	Q = (1 - layer_heights_[i]*secant_zenith_angle/min_range_meters) / fac;

	val += secant_zenith_angle*get_cn2_coefficient(power_spectra_[i]->get_coefficient()) * 
	  pow(fac*fac*
	      (rho_1 - (Q*rho_2) +
	       (little_omega*(layer_heights_[i]/(1-layer_heights_[i]/max_range_meters)))).length_squared(),5/6.);
      }
      if(!finite(val)){
	cerr << "refractive_atmospheric_model::caliph_C error - value " 
	     << val 
	     << " not finite\n";
	throw(string("refractive_atmospheric_model::caliph_C"));
      }

      return(val);
    } catch(...) {
      cerr << "refractive_atmospheric_model::caliph_C error\n";
      throw(string("refractive_atmospheric_model::caliph_C"));
    }
  }


  double refractive_atmospheric_model::caliph_D(const emitter & emtr_a,
						const emitter & emtr_b,
						double aperture_diameter_meters,
						int nsteps_in_integration) const {

    if(nsteps_in_integration<=1){
      cerr << "refractive_atmospheric_model::caliph_D error - \n"
	   << "number of steps in integration "
	   << nsteps_in_integration
	   << " out of range\n";
	throw(string("refractive_atmospheric_model::caliph_D"));
    }

    if(aperture_diameter_meters<=0){
      cerr << "refractive_atmospheric_model::caliph_D error - \n"
	   << "invalid aperture diameter "
	   << aperture_diameter_meters
	   << " passed to this function\n";
      throw(string("refractive_atmospheric_model::caliph_D"));
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

      double val = 0;

      double Q, omega_amplitude;

      vector<double> numerator_args(2,-5/6.);
      vector<double> denominator_args(1,1.);
      int npowerspec = this->power_spectra_.size();

      for(int i=0; i<npowerspec; i++){
	if(layer_heights_[i]>min_range_meters) continue;

	Q = (1 - layer_heights_[i]*secant_zenith_angle/min_range_meters) /
	  (1 - layer_heights_[i]*secant_zenith_angle/max_range_meters);

	omega_amplitude = 
	  little_omega.length()*(layer_heights_[i]/(1-(layer_heights_[i]/max_range_meters)));

	val += secant_zenith_angle*get_cn2_coefficient(power_spectra_[i]->get_coefficient()) * 
	  pow((1 - layer_heights_[i]*secant_zenith_angle/max_range_meters),5/3.) *
	  Tyler_F_2(Q, 
	      omega_amplitude,
	      nsteps_in_integration);
      }
      if(!finite(val)){
	cerr << "refractive_atmospheric_model::caliph_D error - value " 
	     << val 
	     << " not finite\n";
	throw(string("refractive_atmospheric_model::caliph_D"));
      }

      return(val);
    } catch(...) {
      cerr << "refractive_atmospheric_model::caliph_D error\n";
      throw(string("refractive_atmospheric_model::caliph_D"));
    }
  }

  /*
  double refractive_atmospheric_model::tmp_tyler_G_1_minus(const emitter & emtr_a,
							   const emitter & emtr_b,
							   double aperture_diameter_meters,
							   const three_vector & rho) const {

    if(rho.length_squared()>1){
      cerr << "refractive_atmospheric_model::tmp_tyler_G_1_minus error - \n"
	   << "three vector passed to this function has an amplitude "
	   << rho.length()
	   << " greater than unity\n";
      rho.print(cerr, "rho ");
      throw(string("refractive_atmospheric_model::tmp_tyler_G_1_minus"));
    }

    if(aperture_diameter_meters<=0){
      cerr << "refractive_atmospheric_model::tmp_tyler_G_1_minus error - \n"
	   << "invalid aperture diameter "
	   << aperture_diameter_meters
	   << " passed to this function\n";
      throw(string("refractive_atmospheric_model::tmp_tyler_G_1_minus"));
    }

    try{
      double secant_zenith_angle;
      double max_range_meters;
      double min_range_meters;
      double ghat_arg;
      three_vector little_omega;

      get_constants(emtr_a,
		    emtr_b,
		    this->ground_ref_frame_,
		    aperture_diameter_meters,
		    secant_zenith_angle,
		    max_range_meters,
		    min_range_meters,
		    little_omega);

      double val = 0;
      double Q;
      three_vector tmp;
      int npowerspec = this->power_spectra_.size();

      for(int i=0; i<npowerspec; i++){
	if(layer_heights_[i]>min_range_meters) continue;

	Q = (1 - layer_heights_[i]*secant_zenith_angle/min_range_meters) /
	  (1 - layer_heights_[i]*secant_zenith_angle/max_range_meters);

	tmp = rho*Q - little_omega*(layer_heights_[i]/(1-layer_heights_[i]/max_range_meters));

	ghat_arg = tmp.length();

	val += secant_zenith_angle*get_cn2_coefficient(power_spectra_[i]->get_coefficient()) * 
	  pow((1 - layer_heights_[i]*secant_zenith_angle/max_range_meters),5/3.) *
	  dot_product(rho, tmp)*
	  G_hat(ghat_arg);
      }
      val *= 4;

      if(!finite(val)){
	cerr << "refractive_atmospheric_model::tmp_tyler_G_1_minus error - value " 
	     << val 
	     << " not finite\n";
	throw(string("refractive_atmospheric_model::tmp_tyler_G_1_minus"));
      }

      return(val);
    } catch(...) {
      cerr << "refractive_atmospheric_model::tmp_tyler_G_1_minus error\n";
      throw(string("refractive_atmospheric_model::tmp_tyler_G_1_minus"));
    }
  }

  double refractive_atmospheric_model::tmp_tyler_G_1_plus(const emitter & emtr_a,
							  const emitter & emtr_b,
							  double aperture_diameter_meters,
							  const three_vector & rho) const {

    if(rho.length_squared()>1){
      cerr << "refractive_atmospheric_model::tmp_tyler_G_1_plus error - \n"
	   << "three vector passed to this function has an amplitude "
	   << rho.length()
	   << " greater than unity\n";
      rho.print(cerr, "rho ");
      throw(string("refractive_atmospheric_model::tmp_tyler_G_1_plus"));
    }

    if(aperture_diameter_meters<=0){
      cerr << "refractive_atmospheric_model::tmp_tyler_G_1_plus error - \n"
	   << "invalid aperture diameter "
	   << aperture_diameter_meters
	   << " passed to this function\n";
      throw(string("refractive_atmospheric_model::tmp_tyler_G_1_plus"));
    }

    try{
      double secant_zenith_angle;
      double max_range_meters;
      double min_range_meters;
      double ghat_arg;
      three_vector little_omega;

      get_constants(emtr_a,
		    emtr_b,
		    this->ground_ref_frame_,
		    aperture_diameter_meters,
		    secant_zenith_angle,
		    max_range_meters,
		    min_range_meters,
		    little_omega);

      double val = 0;
      double Q;
      three_vector tmp;
      int npowerspec = this->power_spectra_.size();

      for(int i=0; i<npowerspec; i++){
	if(layer_heights_[i]>min_range_meters) continue;

	Q = (1 - layer_heights_[i]*secant_zenith_angle/min_range_meters) /
	  (1 - layer_heights_[i]*secant_zenith_angle/max_range_meters);

	tmp = rho + little_omega*(layer_heights_[i]/(1-layer_heights_[i]/max_range_meters));
	tmp *= (1/Q);

	ghat_arg = tmp.length();

	val += secant_zenith_angle*get_cn2_coefficient(power_spectra_[i]->get_coefficient()) * 
	  pow((1 - layer_heights_[i]*secant_zenith_angle/max_range_meters),5/3.) *
	  dot_product(rho, tmp)*
	  pow(Q,5/3.) * 
	  G_hat(ghat_arg);
      }
      val *= 4;

      if(!finite(val)){
	cerr << "refractive_atmospheric_model::tmp_tyler_G_1_plus error - value " 
	     << val 
	     << " not finite\n";
	throw(string("refractive_atmospheric_model::tmp_tyler_G_1_plus"));
      }

      return(val);
    } catch(...) {
      cerr << "refractive_atmospheric_model::tmp_tyler_G_1_plus error\n";
      throw(string("refractive_atmospheric_model::tmp_tyler_G_1_plus"));
    }
  }

  double refractive_atmospheric_model::caliph_G_3(const emitter & emtr_a,
						const emitter & emtr_b,
						double aperture_diameter_meters,
						int nsteps_in_integration) const {

    if(nsteps_in_integration<=1){
      cerr << "refractive_atmospheric_model::caliph_G_3 error - \n"
	   << "number of steps in integration "
	   << nsteps_in_integration
	   << " out of range\n";
	throw(string("refractive_atmospheric_model::caliph_G_3"));
    }

    if(aperture_diameter_meters<=0){
      cerr << "refractive_atmospheric_model::caliph_G_3 error - \n"
	   << "invalid aperture diameter "
	   << aperture_diameter_meters
	   << " passed to this function\n";
      throw(string("refractive_atmospheric_model::caliph_G_3"));
    }
    try{
      double secant_zenith_angle;
      double max_range_meters;
      double min_range_meters;
      three_vector little_omega;

      get_constants(emtr_a,
		    emtr_b,
		    this->ground_ref_frame_,
		    aperture_diameter_meters,
		    secant_zenith_angle,
		    max_range_meters,
		    min_range_meters,
		    little_omega);

      double val = 0;
      double step_size, integral_value, rho;
      double Q, omega_amplitude, omega_amplitude_squared, fac;
      int npowerspec = this->power_spectra_.size();
      for(int i=0; i<npowerspec; i++){
	if(layer_heights_[i]>min_range_meters) continue;

	Q = (1 - (layer_heights_[i]*secant_zenith_angle/min_range_meters)) /
	  (1 - (layer_heights_[i]*secant_zenith_angle/max_range_meters));

	omega_amplitude = 
	  little_omega.length()*(layer_heights_[i]/(1-(layer_heights_[i]/max_range_meters)));

	val += -4*secant_zenith_angle*get_cn2_coefficient(power_spectra_[i]->get_coefficient()) * 
	  pow((1 - layer_heights_[i]*secant_zenith_angle/max_range_meters),5/3.) * 
	  G_3(Q,omega_amplitude,nsteps_in_integration);
      }

      if(!finite(val)){
	cerr << "refractive_atmospheric_model::caliph_G_3 error - value " 
	     << val 
	     << " not finite\n";
	throw(string("refractive_atmospheric_model::caliph_G_3"));
      }

      return(val);
    } catch(...) {
      cerr << "refractive_atmospheric_model::caliph_G_3 error\n";
      throw(string("refractive_atmospheric_model::caliph_G_3"));
    }
  }
  */

  double refractive_atmospheric_model::caliph_E(const emitter & emtr_a,
						const emitter & emtr_b,
						double aperture_diameter_meters,
						int nsteps_in_integration) const {

    if(nsteps_in_integration<=1){
      cerr << "refractive_atmospheric_model::caliph_E error - \n"
	   << "number of steps in integration "
	   << nsteps_in_integration
	   << " out of range\n";
	throw(string("refractive_atmospheric_model::caliph_E"));
    }

    if(aperture_diameter_meters<=0){
      cerr << "refractive_atmospheric_model::caliph_E error - \n"
	   << "invalid aperture diameter "
	   << aperture_diameter_meters
	   << " passed to this function\n";
      throw(string("refractive_atmospheric_model::caliph_E"));
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

      double val = 0;
      double step_size, integral_value, rho;
      double Q, omega_amplitude, omega_amplitude_squared, fac;
      int npowerspec = this->power_spectra_.size();
      for(int i=0; i<npowerspec; i++){
	if(layer_heights_[i]>min_range_meters) continue;

	Q = (1 - (layer_heights_[i]*secant_zenith_angle/min_range_meters)) /
	  (1 - (layer_heights_[i]*secant_zenith_angle/max_range_meters));

	omega_amplitude = 
	  little_omega.length()*(layer_heights_[i]/(1-(layer_heights_[i]/max_range_meters)));
	omega_amplitude_squared = omega_amplitude*omega_amplitude;

	integral_value = 0;

	if(Q > omega_amplitude){
	  step_size = (Q - omega_amplitude)/(double)nsteps_in_integration;
	  rho = .5*step_size;

	  do {
	    integral_value += (step_size*rho*rho*rho*Tyler_G_hat(rho));
	    rho += step_size;
	  } while(rho<(Q-omega_amplitude));
	}

	if(omega_amplitude>three_frame::precision){
	  step_size = (Q + omega_amplitude - fabs(Q - omega_amplitude))/(double)nsteps_in_integration;
	  rho = fabs(Q - omega_amplitude) + .5*step_size;
	  while(rho<(Q+omega_amplitude)){

	    fac = (omega_amplitude_squared + rho*rho-Q*Q)/2.0/omega_amplitude/rho;
	    if(fac>1){
	      cerr << "caliph E imminent failure\n";
	      exit(-1);
	    }
	    // This one is strongly preferred based on comparison to G_3
	    integral_value += (step_size*rho*rho*rho*Tyler_G_hat(rho)/M_PI) * acos(fac);

	    // This one seems mathematically correct...
	    //integral_value += (step_size*rho*rho*rho*G_hat(rho)/M_PI) * asin(sqrt(1-fac*fac));

	    rho += step_size;
	  }
	}

	val += secant_zenith_angle*get_cn2_coefficient(power_spectra_[i]->get_coefficient()) * 
	  pow((1 - layer_heights_[i]*secant_zenith_angle/max_range_meters),5/3.) * 
	  integral_value / Q / Q / Q ;
      }
      val *= 16;

      if(!finite(val)){
	cerr << "refractive_atmospheric_model::caliph_E error - value " 
	     << val 
	     << " not finite\n";
	throw(string("refractive_atmospheric_model::caliph_E"));
      }

      return(val);
    } catch(...) {
      cerr << "refractive_atmospheric_model::caliph_E error\n";
      throw(string("refractive_atmospheric_model::caliph_E"));
    }
  }

  double refractive_atmospheric_model::caliph_F(const emitter & emtr_a,
						const emitter & emtr_b,
						double aperture_diameter_meters,
						int nsteps_in_integration) const {

    if(nsteps_in_integration<=1){
      cerr << "refractive_atmospheric_model::caliph_F error - \n"
	   << "number of steps in integration "
	   << nsteps_in_integration
	   << " out of range\n";
	throw(string("refractive_atmospheric_model::caliph_F"));
    }

    if(aperture_diameter_meters<=0){
      cerr << "refractive_atmospheric_model::caliph_F error - \n"
	   << "invalid aperture diameter "
	   << aperture_diameter_meters
	   << " passed to this function\n";
      throw(string("refractive_atmospheric_model::caliph_F"));
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

      

      double val = 0;
      double step_size, integral_value, rho;
      double Q, omega_amplitude, omega_amplitude_squared, fac;
      double little_omega_length = little_omega.length();

      if(little_omega_length<three_frame::precision) 
	return(0);

      int npowerspec = this->power_spectra_.size();
      for(int i=0; i<npowerspec; i++){
	if(layer_heights_[i]>min_range_meters) continue;

	Q = (1 - layer_heights_[i]*secant_zenith_angle/min_range_meters) /
	  (1 - layer_heights_[i]*secant_zenith_angle/max_range_meters);

	omega_amplitude = 
	  little_omega_length*(layer_heights_[i]/(1-layer_heights_[i]/max_range_meters));

	omega_amplitude_squared = omega_amplitude*omega_amplitude;

	integral_value = 0;

	step_size = (Q + omega_amplitude - fabs(Q - omega_amplitude))/(double)nsteps_in_integration;
	rho = fabs(Q - omega_amplitude) + .5*step_size;
	while(rho<(Q+omega_amplitude)){
	  fac = fabs((omega_amplitude_squared + rho*rho-Q*Q)/2.0/omega_amplitude/rho);
	  if(fac>1){
	    cerr << "caliph F imminent failure\n";
	    exit(-1);
	  }
	  integral_value += (step_size*rho*rho*Tyler_G_hat(rho)) *sqrt(1-fac*fac);
	  rho += step_size;
	}

	val += secant_zenith_angle*get_cn2_coefficient(power_spectra_[i]->get_coefficient()) * 
	  pow((1 - layer_heights_[i]*secant_zenith_angle/max_range_meters),5/3.) * 
	  integral_value * omega_amplitude / Q / Q / Q;
      }
      val *= 32/M_PI;

      if(!finite(val)){
	cerr << "refractive_atmospheric_model::caliph_F error - value " 
	     << val 
	     << " not finite\n";
	throw(string("refractive_atmospheric_model::caliph_F"));
      }

      return(val);
    } catch(...) {
      cerr << "refractive_atmospheric_model::caliph_F error\n";
      throw(string("refractive_atmospheric_model::caliph_F"));
    }
  }

  double refractive_atmospheric_model::caliph_F_bar(const emitter & emtr_a,
						    const emitter & emtr_b,
						    double aperture_diameter_meters,
						    int nsteps_in_integration) const {

    if(nsteps_in_integration<=1){
      cerr << "refractive_atmospheric_model::caliph_F_bar error - \n"
	   << "number of steps in integration "
	   << nsteps_in_integration
	   << " out of range\n";
	throw(string("refractive_atmospheric_model::caliph_F_bar"));
    }

    if(aperture_diameter_meters<=0){
      cerr << "refractive_atmospheric_model::caliph_F_bar error - \n"
	   << "invalid aperture diameter "
	   << aperture_diameter_meters
	   << " passed to this function\n";
      throw(string("refractive_atmospheric_model::caliph_F_bar"));
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

      double val = 0;
      double step_size, integral_value, rho;
      double Q, omega_amplitude, omega_amplitude_squared, fac;
      int npowerspec = this->power_spectra_.size();
      for(int i=0; i<npowerspec; i++){
	if(layer_heights_[i]>min_range_meters) continue;

	Q = (1 - layer_heights_[i]*secant_zenith_angle/min_range_meters) /
	  (1 - layer_heights_[i]*secant_zenith_angle/max_range_meters);

	omega_amplitude = 
	  little_omega.length()*(layer_heights_[i]/(1-layer_heights_[i]/max_range_meters));

	if(omega_amplitude==0) 
	  continue;

	omega_amplitude_squared = omega_amplitude*omega_amplitude;

	integral_value = 0;

	step_size = (Q + omega_amplitude - fabs(Q - omega_amplitude))/(double)nsteps_in_integration;
	rho = fabs(Q - omega_amplitude) + .5*step_size;
	while(rho<(Q+omega_amplitude)){
	  fac = (omega_amplitude_squared + rho*rho-Q*Q)/2.0/omega_amplitude/rho;
	  if(fac>1){
	    cerr << "caliph_F_bar imminent failure\n";
	    exit(-1);
	  }
	  integral_value += (step_size*rho*rho*rho*Tyler_G_hat(rho)) * fac*sqrt(1-fac*fac);
	  rho += step_size;
	}

	val += secant_zenith_angle*get_cn2_coefficient(power_spectra_[i]->get_coefficient()) * 
	  pow((1 - layer_heights_[i]*secant_zenith_angle/max_range_meters),5/3.) * 
	  integral_value / Q / Q / Q;
      }
      val *= 16/M_PI;

      if(!finite(val)){
	cerr << "refractive_atmospheric_model::caliph_F_bar error - value " 
	     << val 
	     << " not finite\n";
	throw(string("refractive_atmospheric_model::caliph_F_bar"));
      }
      return(val);
    } catch(...) {
      cerr << "refractive_atmospheric_model::caliph_F_bar error\n";
      throw(string("refractive_atmospheric_model::caliph_F_bar"));
    }
  }

  double refractive_atmospheric_model::caliph_G(const emitter & emtr_a,
						const emitter & emtr_b,
						double aperture_diameter_meters,
						const three_vector & rho) const {

    if(rho.length_squared()>1){
      cerr << "refractive_atmospheric_model::caliph_G error - \n"
	   << "three vector passed to this function has an amplitude "
	   << rho.length()
	   << " greater than unity\n";
      rho.print(cerr, "rho ");
      throw(string("refractive_atmospheric_model::caliph_G"));
    }

    if(aperture_diameter_meters<=0){
      cerr << "refractive_atmospheric_model::caliph_G error - \n"
	   << "invalid aperture diameter "
	   << aperture_diameter_meters
	   << " passed to this function\n";
      throw(string("refractive_atmospheric_model::caliph_G"));
    }

    try{
      double secant_zenith_angle;
      double max_range_meters;
      double min_range_meters;
      double ghat_arg;
      three_vector little_omega;

      Tyler_get_constants(emtr_a,
		    emtr_b,
		    this->ground_ref_frame_,
		    aperture_diameter_meters,
		    secant_zenith_angle,
		    max_range_meters,
		    min_range_meters,
		    little_omega);

      double val = 0;

      double Q;
      int npowerspec = this->power_spectra_.size();
      for(int i=0; i<npowerspec; i++){
	if(layer_heights_[i]>min_range_meters) continue;

	Q = (1 - layer_heights_[i]*secant_zenith_angle/min_range_meters) /
	  (1 - layer_heights_[i]*secant_zenith_angle/max_range_meters);

	ghat_arg = (rho*Q - little_omega*(layer_heights_[i]/(1-layer_heights_[i]/max_range_meters))).length();

	val += secant_zenith_angle*get_cn2_coefficient(power_spectra_[i]->get_coefficient()) * 
	  pow((1 - layer_heights_[i]*secant_zenith_angle/max_range_meters),5/3.) *
	  Q * Tyler_G_hat(ghat_arg);
      }
      val *= 4;

      if(!finite(val)){
	cerr << "refractive_atmospheric_model::caliph_G error - value " 
	     << val 
	     << " not finite\n";
	throw(string("refractive_atmospheric_model::caliph_G"));
      }

      return(val);
    } catch(...) {
      cerr << "refractive_atmospheric_model::caliph_G error\n";
      throw(string("refractive_atmospheric_model::caliph_G"));
    }
  }

  double refractive_atmospheric_model::caliph_H(const emitter & emtr_a,
						const emitter & emtr_b,
						double aperture_diameter_meters,
						const three_vector & rho) const {

    if(rho.length_squared()>1){
      cerr << "refractive_atmospheric_model::caliph_H error - \n"
	   << "three vector passed to this function has an amplitude "
	   << rho.length()
	   << " greater than unity\n";
      rho.print(cerr, "rho ");
      throw(string("refractive_atmospheric_model::caliph_H"));
    }

    if(aperture_diameter_meters<=0){
      cerr << "refractive_atmospheric_model::caliph_H error - \n"
	   << "invalid aperture diameter "
	   << aperture_diameter_meters
	   << " passed to this function\n";
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

      double little_omega_length = little_omega.length();
      if(little_omega_length < three_frame::precision)
	return(0);

      double ghat_arg;
      double val = 0;

      double Q;
      int npowerspec = this->power_spectra_.size();
      for(int i=0; i<npowerspec; i++){
	if(layer_heights_[i]>min_range_meters) continue;

	Q = (1 - layer_heights_[i]*secant_zenith_angle/min_range_meters) /
	  (1 - layer_heights_[i]*secant_zenith_angle/max_range_meters);

	ghat_arg = (rho*Q - little_omega*(layer_heights_[i]/(1-layer_heights_[i]/max_range_meters))).length();

	val += secant_zenith_angle*get_cn2_coefficient(power_spectra_[i]->get_coefficient()) * 
	  pow((1 - layer_heights_[i]*secant_zenith_angle/max_range_meters),5/3.) *
	  little_omega_length*(layer_heights_[i]/(1-layer_heights_[i]/max_range_meters))
	  * Tyler_G_hat(ghat_arg);
      }
      val *= 4;

      if(!finite(val)){
	cerr << "refractive_atmospheric_model::caliph_H error - value " 
	     << val 
	     << " not finite\n";
	throw(string("refractive_atmospheric_model::caliph_H"));
      }

      return(val);
    } catch(...) {
      cerr << "refractive_atmospheric_model::caliph_H error\n";
      throw(string("refractive_atmospheric_model::caliph_H"));
    }
  }

  double refractive_atmospheric_model::caliph_I(const emitter & emtr_a,
						const emitter & emtr_b,
						double aperture_diameter_meters,
						const three_vector & rho) const {

    if(rho.length_squared()>1){
      cerr << "refractive_atmospheric_model::caliph_I error - \n"
	   << "three vector passed to this function has an amplitude "
	   << rho.length()
	   << " greater than unity\n";
      rho.print(cerr, "rho ");
      throw(string("refractive_atmospheric_model::caliph_I"));
    }

    if(aperture_diameter_meters<=0){
      cerr << "refractive_atmospheric_model::caliph_I error - \n"
	   << "invalid aperture diameter "
	   << aperture_diameter_meters
	   << " passed to this function\n";
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
      double val = 0;

      double Q;
      int npowerspec = this->power_spectra_.size();
      for(int i=0; i<npowerspec; i++){
	if(layer_heights_[i]>min_range_meters) continue;

	Q = (1 - layer_heights_[i]*secant_zenith_angle/min_range_meters) /
	  (1 - layer_heights_[i]*secant_zenith_angle/max_range_meters);

	ghat_arg = (rho + little_omega*(layer_heights_[i]/(1-layer_heights_[i]/max_range_meters))).length()/Q;

	val += secant_zenith_angle*get_cn2_coefficient(power_spectra_[i]->get_coefficient()) * 
	  pow((1 - layer_heights_[i]*secant_zenith_angle/max_range_meters),5/3.) *
	  pow(Q,2/3.) * 
	  Tyler_G_hat(ghat_arg);
      }
      val *= 4;

      if(!finite(val)){
	cerr << "refractive_atmospheric_model::caliph_I error - value " 
	     << val 
	     << " not finite\n";
	throw(string("refractive_atmospheric_model::caliph_I"));
      }

      return(val);
    } catch(...) {
      cerr << "refractive_atmospheric_model::caliph_I error\n";
      throw(string("refractive_atmospheric_model::caliph_I"));
    }

  }

  double refractive_atmospheric_model::caliph_J(const emitter & emtr_a,
						const emitter & emtr_b,
						double aperture_diameter_meters,
						const three_vector & rho) const {

    if(rho.length_squared()>1){
      cerr << "refractive_atmospheric_model::caliph_J error - \n"
	   << "three vector passed to this function has an amplitude "
	   << rho.length()
	   << " greater than unity\n";
      rho.print(cerr, "rho ");
      throw(string("refractive_atmospheric_model::caliph_J"));
    }

    if(aperture_diameter_meters<=0){
      cerr << "refractive_atmospheric_model::caliph_J error - \n"
	   << "invalid aperture diameter "
	   << aperture_diameter_meters
	   << " passed to this function\n";
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

      double little_omega_length = little_omega.length();
      if(little_omega_length < three_frame::precision)
	return(0);

      double ghat_arg;
      double val = 0;

      double Q;
      int npowerspec = this->power_spectra_.size();
      for(int i=0; i<npowerspec; i++){
	if(layer_heights_[i]>min_range_meters) continue;

	Q = (1 - layer_heights_[i]*secant_zenith_angle/min_range_meters) /
	  (1 - layer_heights_[i]*secant_zenith_angle/max_range_meters);

	ghat_arg = (rho + little_omega*(layer_heights_[i]/(1-layer_heights_[i]/max_range_meters))).length()/Q;

	val += secant_zenith_angle*get_cn2_coefficient(power_spectra_[i]->get_coefficient()) * 
	  pow((1 - layer_heights_[i]*secant_zenith_angle/max_range_meters),5/3.) *
	  pow(Q, 2/3.) *
	  little_omega_length*(layer_heights_[i]/(1-layer_heights_[i]/max_range_meters)) * 
	  Tyler_G_hat(ghat_arg);
      }
      val *= 4;

      if(!finite(val)){
	cerr << "refractive_atmospheric_model::caliph_J error - value " 
	     << val 
	     << " not finite\n";
	throw(string("refractive_atmospheric_model::caliph_J"));
      }

      return(val);
    } catch(...) {
      cerr << "refractive_atmospheric_model::caliph_J error\n";
      throw(string("refractive_atmospheric_model::caliph_J"));
    }
  }

  double refractive_atmospheric_model::caliph_K(const emitter & emtr_a,
						const emitter & emtr_b,
						double aperture_diameter_meters,
						int nsteps_in_integration) const {

    if(nsteps_in_integration<=1){
      cerr << "refractive_atmospheric_model::caliph_K error - \n"
	   << "number of steps in integration "
	   << nsteps_in_integration
	   << " out of range\n";
	throw(string("refractive_atmospheric_model::caliph_K"));
    }

    if(aperture_diameter_meters<=0){
      cerr << "refractive_atmospheric_model::caliph_K error - \n"
	   << "invalid aperture diameter "
	   << aperture_diameter_meters
	   << " passed to this function\n";
      throw(string("refractive_atmospheric_model::caliph_K"));
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

      double val = 0;
      double step_size, integral_value, rho;
      double Q, omega_amplitude, omega_amplitude_squared, fac;
      int npowerspec = this->power_spectra_.size();
      for(int i=0; i<npowerspec; i++){
	if(layer_heights_[i]>min_range_meters) continue;

	Q = (1 - layer_heights_[i]*secant_zenith_angle/min_range_meters) /
	  (1 - layer_heights_[i]*secant_zenith_angle/max_range_meters);

	omega_amplitude = 
	  little_omega.length()*(layer_heights_[i]/(1-layer_heights_[i]/max_range_meters));
	if(omega_amplitude==0) 
	  continue;

	omega_amplitude_squared = omega_amplitude*omega_amplitude;

	step_size = 2*omega_amplitude/(double)nsteps_in_integration;
	integral_value = 0;
	rho = fabs(Q - omega_amplitude) + .5*step_size;
	while(rho<(Q+omega_amplitude)){
	  fac = (omega_amplitude_squared + rho*rho-Q*Q)/2.0/omega_amplitude/rho;
	  if(fac>1){
	    cerr << "caliph K imminent failure\n";
	    exit(-1);
	  }
	  integral_value += (step_size*rho*rho*Tyler_G_hat(rho)) * 
	    sqrt(1-fac*fac);
	  rho += step_size;
	}

	val += secant_zenith_angle*get_cn2_coefficient(power_spectra_[i]->get_coefficient()) * 
	  pow((1 - layer_heights_[i]*secant_zenith_angle/max_range_meters),5/3.) * 
	  integral_value / Q / Q;
      }
      val *= 8/M_PI;

      if(!finite(val)){
	cerr << "refractive_atmospheric_model::caliph_K error - value " 
	     << val 
	     << " not finite\n";
	throw(string("refractive_atmospheric_model::caliph_K"));
      }

      return(val);
    } catch(...) {
      cerr << "refractive_atmospheric_model::caliph_K error\n";
      throw(string("refractive_atmospheric_model::caliph_K"));
    }
  }

  double refractive_atmospheric_model::caliph_L(const emitter & emtr_a,
						const emitter & emtr_b,
						double aperture_diameter_meters,
						int nsteps_in_integration) const {

    if(nsteps_in_integration<=1){
      cerr << "refractive_atmospheric_model::caliph_L error - \n"
	   << "number of steps in integration "
	   << nsteps_in_integration
	   << " out of range\n";
	throw(string("refractive_atmospheric_model::caliph_L"));
    }

    if(aperture_diameter_meters<=0){ 
      cerr << "refractive_atmospheric_model::caliph_L error - \n"
	   << "invalid aperture diameter "
	   << aperture_diameter_meters
	   << " passed to this function\n";
      throw(string("refractive_atmospheric_model::caliph_L"));
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

      double val = 0;
      double step_size, integral_value, rho;
      double Q, omega_amplitude, omega_amplitude_squared, fac;
      int npowerspec = this->power_spectra_.size();
      for(int i=0; i<npowerspec; i++){
	if(layer_heights_[i]>min_range_meters) continue;

	Q = (1 - layer_heights_[i]*secant_zenith_angle/min_range_meters) /
	  (1 - layer_heights_[i]*secant_zenith_angle/max_range_meters);

	omega_amplitude = 
	  little_omega.length()*(layer_heights_[i]/(1-layer_heights_[i]/max_range_meters));
	if(omega_amplitude==0) 
	  continue;

	omega_amplitude_squared = omega_amplitude*omega_amplitude;

	step_size = 2*omega_amplitude/(double)nsteps_in_integration;
	integral_value = 0;
	rho = fabs(1 - omega_amplitude)/Q + .5*step_size;
	while(rho<(1+omega_amplitude)/Q){
	  fac = (omega_amplitude_squared + Q*Q*rho*rho-1)/2.0/Q/omega_amplitude/rho;
	  if(fac>1){
	    cerr << "caliph L imminent failure\n";
	    exit(-1);
	  }
	  integral_value += (step_size*rho*rho*Tyler_G_hat(rho)) * sqrt(1-fac*fac);
	  rho += step_size;
	}

	val += secant_zenith_angle*get_cn2_coefficient(power_spectra_[i]->get_coefficient()) * 
	  pow((1 - layer_heights_[i]*secant_zenith_angle/max_range_meters),5/3.) * 
	  integral_value * pow(Q, 11/3.);
      }
      val *= 8/M_PI;

      if(!finite(val)){
	cerr << "refractive_atmospheric_model::caliph_L error - value " 
	     << val 
	     << " not finite\n";
	throw(string("refractive_atmospheric_model::caliph_L"));
      }

      return(val);
    } catch(...) {
      cerr << "refractive_atmospheric_model::caliph_L error\n";
      throw(string("refractive_atmospheric_model::caliph_L"));
    }
  }

  double refractive_atmospheric_model::caliph_M(const emitter & emtr_a,
						const emitter & emtr_b,
						double aperture_diameter_meters) const {

    if(aperture_diameter_meters<=0){ 
      cerr << "refractive_atmospheric_model::caliph_M error - \n"
	   << "invalid aperture diameter "
	   << aperture_diameter_meters
	   << " passed to this function\n";
      throw(string("refractive_atmospheric_model::caliph_M"));
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

      double val = 0;
      double Q, omega, tmp_f1;
      int npowerspec = this->power_spectra_.size();
      for(int i=0; i<npowerspec; i++){
	if(layer_heights_[i]>min_range_meters) continue;

	Q = (1 - layer_heights_[i]*secant_zenith_angle/min_range_meters) /
	  (1 - layer_heights_[i]*secant_zenith_angle/max_range_meters);

	omega = little_omega.length()*(layer_heights_[i]/(1-layer_heights_[i]/max_range_meters));

	if(fabs(Q-1)<three_frame::precision) 
	  tmp_f1 = pow(omega,5/3.);
	else 
	  tmp_f1 = pow(1-Q, 5/3.) * Tyler_F_1(omega/(1-Q));

	val += secant_zenith_angle*get_cn2_coefficient(power_spectra_[i]->get_coefficient()) * 
	  pow((1 - layer_heights_[i]*secant_zenith_angle/max_range_meters),5/3.) * 
	  tmp_f1;
      }
      if(!finite(val)){
	cerr << "refractive_atmospheric_model::caliph_M error - value " 
	     << val 
	     << " not finite\n";
	throw(string("refractive_atmospheric_model::caliph_M"));
      }

      return(val);
    } catch(...) {
      cerr << "refractive_atmospheric_model::caliph_M error\n";
      throw(string("refractive_atmospheric_model::caliph_M"));
    }
  }


  double refractive_atmospheric_model::
  aperture_averaged_phase_covariance(const emitter & emtr_a,
				     const emitter & emtr_b,
				     const aperture & ap,
				     double wavelength_meters,
				     int nsteps_in_integration) const {
    try{check_wavelength(wavelength_meters);}
    catch(...){
      cerr << "refractive_atmospheric_model::phase_covariance error\n";
      throw(string("refractive_atmospheric_model::aperture_averaged_phase_covariance"));
    }
    if(nsteps_in_integration<=0){
      cerr << "refractive_atmospheric_model::aperture_averaged_phase_covariance error\n"
	   << "number of samples in integration " 
	   << nsteps_in_integration
	   << " out of range\n";
      throw(string("refractive_atmospheric_model::aperture_averaged_phase_covariance"));
    }

    circular_aperture circ_ap;
    try{circ_ap=dynamic_cast<const circular_aperture &>(ap);}
    catch(...){
      cerr << "refractive_atmospheric_model::aperture_averaged_phase_covariance error - \n"
	   << "an aperture was passed to this function that is not a circular_aperture\n";
	throw(string("refractive_atmospheric_model::aperture_averaged_phase_covariance"));
    }
    double aperture_diameter_meters = circ_ap.get_diameter();
    if(aperture_diameter_meters<=0){ 
      cerr << "refractive_atmospheric_model::aperture_averaged_phase_covariance error - \n"
	   << "invalid aperture diameter "
	   << aperture_diameter_meters
	   << " passed to this function\n";
      throw(string("refractive_atmospheric_model::aperture_averaged_phase_covariance"));
    }

    double wavenumber = 2 * M_PI / wavelength_meters;

    double val = 
      get_Xi()*pow(aperture_diameter_meters,5/3.)*wavenumber*wavenumber*
      (this->caliph_D(emtr_a, emtr_b, aperture_diameter_meters, nsteps_in_integration) -
       this->caliph_M(emtr_a, emtr_b, aperture_diameter_meters));

    if(!finite(val)){
      cerr << "refractive_atmospheric_model::aperture_averaged_phase_covariance error - value " 
	   << val 
	   << " not finite\n";
      throw(string("refractive_atmospheric_model::aperture_averaged_phase_covariance"));
    }

    return(val);
	   
  }
  
  double refractive_atmospheric_model::
  aperture_averaged_tilt_phase_covariance(const emitter & emtr_a,
					  const emitter & emtr_b,
					  const aperture & ap,
					  double wavelength_meters,
					  int nsteps_in_integration) const {
    try{check_wavelength(wavelength_meters);}
    catch(...){
      cerr << "refractive_atmospheric_model::aperture_averaged_tilt_phase_covariance error\n";
      throw(string("refractive_atmospheric_model::aperture_averaged_tilt_phase_covariance"));
    }
    if(nsteps_in_integration<=0){
      cerr << "refractive_atmospheric_model::aperture_averaged_tilt_phase_covariance error\n"
	   << "number of samples in integration " 
	   << nsteps_in_integration
	   << " out of range\n";
      throw(string("refractive_atmospheric_model::aperture_averaged_tilt_phase_covariance"));
    }
    circular_aperture circ_ap;
    try{circ_ap=dynamic_cast<const circular_aperture &>(ap);}
    catch(...){
      cerr << "refractive_atmospheric_model::aperture_averaged_tilt_phase_covariance error - \n"
	   << "an aperture was passed to this function that is not a circular_aperture\n";
	throw(string("refractive_atmospheric_model::aperture_averaged_tilt_phase_covariance"));
    }
    double aperture_diameter_meters = circ_ap.get_diameter();
    if(aperture_diameter_meters<=0){
      cerr << "refractive_atmospheric_model::aperture_averaged_tilt_phase_covariance error - \n"
	   << "invalid aperture diameter "
	   << aperture_diameter_meters
	   << " passed to this function\n";
      throw(string("refractive_atmospheric_model::aperture_averaged_tilt_phase_covariance"));
    }


    double wavenumber = 2 * M_PI / wavelength_meters;

    double val = 
      get_Xi()*pow(aperture_diameter_meters,5/3.)*wavenumber*wavenumber*
      (.25*this->caliph_F(emtr_a, emtr_b, aperture_diameter_meters, nsteps_in_integration) -
       .5*this->caliph_E(emtr_a, emtr_b, aperture_diameter_meters, nsteps_in_integration));

    if(!finite(val)){
      cerr << "refractive_atmospheric_model::aperture_averaged_tilt_phase_covariance error - value " 
	   << val 
	   << " not finite\n";
      throw(string("refractive_atmospheric_model::aperture_averaged_tilt_phase_covariance"));
    }
	   
    return(val);

  }
  
  double refractive_atmospheric_model::
  aperture_averaged_parallel_tilt_phase_covariance(const emitter & emtr_a,
					  const emitter & emtr_b,
					  const aperture & ap,
					  double wavelength_meters,
					  int nsteps_in_integration) const {
    try{check_wavelength(wavelength_meters);}
    catch(...){
      cerr << "refractive_atmospheric_model::aperture_averaged_parallel_tilt_phase_covariance error\n";
      throw(string("refractive_atmospheric_model::aperture_averaged_parallel_tilt_phase_covariance"));
    }
    if(nsteps_in_integration<=0){
      cerr << "refractive_atmospheric_model::aperture_averaged_parallel_tilt_phase_covariance error\n"
	   << "number of samples in integration " 
	   << nsteps_in_integration
	   << " out of range\n";
      throw(string("refractive_atmospheric_model::aperture_averaged_parallel_tilt_phase_covariance"));
    }
    circular_aperture circ_ap;
    try{circ_ap=dynamic_cast<const circular_aperture &>(ap);}
    catch(...){
      cerr << "refractive_atmospheric_model::aperture_averaged_parallel_tilt_phase_covariance error - \n"
	   << "an aperture was passed to this function that is not a circular_aperture\n";
	throw(string("refractive_atmospheric_model::aperture_averaged_parallel_tilt_phase_covariance"));
    }
    double aperture_diameter_meters = circ_ap.get_diameter();
    if(aperture_diameter_meters<=0){
      cerr << "refractive_atmospheric_model::aperture_averaged_parallel_tilt_phase_covariance error - \n"
	   << "invalid aperture diameter "
	   << aperture_diameter_meters
	   << " passed to this function\n";
      throw(string("refractive_atmospheric_model::aperture_averaged_parallel_tilt_phase_covariance"));
    }


    double wavenumber = 2 * M_PI / wavelength_meters;

    double val = 
      get_Xi()*pow(aperture_diameter_meters,5/3.)*wavenumber*wavenumber*
      (this->caliph_F(emtr_a, emtr_b, aperture_diameter_meters, nsteps_in_integration)*5/16. -
       this->caliph_F_bar(emtr_a, emtr_b, aperture_diameter_meters, nsteps_in_integration)*9/16. -
       this->caliph_E(emtr_a, emtr_b, aperture_diameter_meters, nsteps_in_integration)/4.);

    if(!finite(val)){
      cerr << "refractive_atmospheric_model::aperture_averaged_parallel_tilt_phase_covariance error - value " 
	   << val 
	   << " not finite\n";
      throw(string("refractive_atmospheric_model::aperture_averaged_parallel_tilt_phase_covariance"));
    }
	   
    return(val);

  }
  
  double refractive_atmospheric_model::
  aperture_averaged_perpendicular_tilt_phase_covariance(const emitter & emtr_a,
							const emitter & emtr_b,
							const aperture & ap,
							double wavelength_meters,
							int nsteps_in_integration) const {
    try{check_wavelength(wavelength_meters);}
    catch(...){
      cerr << "refractive_atmospheric_model::aperture_averaged_perpendicular_tilt_phase_covariance error\n";
      throw(string("refractive_atmospheric_model::aperture_averaged_perpendicular_tilt_phase_covariance"));
    }
    if(nsteps_in_integration<=0){
      cerr << "refractive_atmospheric_model::aperture_averaged_perpendicular_tilt_phase_covariance error\n"
	   << "number of samples in integration " 
	   << nsteps_in_integration
	   << " out of range\n";
      throw(string("refractive_atmospheric_model::aperture_averaged_perpendicular_tilt_phase_covariance"));
    }
    circular_aperture circ_ap;
    try{circ_ap=dynamic_cast<const circular_aperture &>(ap);}
    catch(...){
      cerr << "refractive_atmospheric_model::aperture_averaged_perpendicular_tilt_phase_covariance error - \n"
	   << "an aperture was passed to this function that is not a circular_aperture\n";
	throw(string("refractive_atmospheric_model::aperture_averaged_perpendicular_tilt_phase_covariance"));
    }
    double aperture_diameter_meters = circ_ap.get_diameter();
    if(aperture_diameter_meters<=0){
      cerr << "refractive_atmospheric_model::aperture_averaged_perpendicular_tilt_phase_covariance error - \n"
	   << "invalid aperture diameter "
	   << aperture_diameter_meters
	   << " passed to this function\n";
      throw(string("refractive_atmospheric_model::aperture_averaged_perpendicular_tilt_phase_covariance"));
    }


    double wavenumber = 2 * M_PI / wavelength_meters;

    double val = 
      get_Xi()*pow(aperture_diameter_meters,5/3.)*wavenumber*wavenumber*
      (-this->caliph_F(emtr_a, emtr_b, aperture_diameter_meters, nsteps_in_integration)/16. -
       this->caliph_F_bar(emtr_a, emtr_b, aperture_diameter_meters, nsteps_in_integration)*9/16. -
       this->caliph_E(emtr_a, emtr_b, aperture_diameter_meters, nsteps_in_integration)/4.);

    if(!finite(val)){
      cerr << "refractive_atmospheric_model::aperture_averaged_perpendicular_tilt_phase_covariance error - value " 
	   << val 
	   << " not finite\n";
      throw(string("refractive_atmospheric_model::aperture_averaged_perpendicular_tilt_phase_covariance"));
    }
	   
    return(val);

  }
  
  double refractive_atmospheric_model::
  phase_covariance(const emitter & emtr_a,
		   const three_point & pupil_location_one,
		   const emitter & emtr_b,
		   const three_point & pupil_location_two,
		   const aperture & ap,
		   double wavelength_meters,
		   int nsteps_in_integration) const {

    try{check_wavelength(wavelength_meters);}
    catch(...){
      cerr << "refractive_atmospheric_model::phase_covariance error\n";
      throw(string("refractive_atmospheric_model::phase_covariance"));
    }
    if(nsteps_in_integration<=0){
      cerr << "refractive_atmospheric_model::phase_covariance error\n"
	   << "number of samples in integration " 
	   << nsteps_in_integration
	   << " out of range\n";
      throw(string("refractive_atmospheric_model::phase_covariance"));
    }

    if(fabs(pupil_location_one.z(ground_ref_frame_))>three_frame::precision ||
       fabs(pupil_location_two.z(ground_ref_frame_))>three_frame::precision){
      cerr << "refractive_atmospheric_model::phase_covariance error - vector does not lie in plane of aperture\n";
      pupil_location_one.print(cerr, "pupil location a ");
      pupil_location_two.print(cerr, "pupil location b ");
      throw(string("refractive_atmospheric_model::phase_covariance"));
    }

    circular_aperture circ_ap;
    try{circ_ap=dynamic_cast<const circular_aperture &>(ap);}
    catch(...){
      cerr << "refractive_atmospheric_model::phase_covariance error - \n"
	   << "an aperture was passed to this function that is not a circular_aperture\n";
	throw(string("refractive_atmospheric_model::phase_covariance"));
    }
    double aperture_diameter_meters = circ_ap.get_diameter();
    if(aperture_diameter_meters<=0){
      cerr << "refractive_atmospheric_model::phase_covariance error - \n"
	   << "invalid aperture diameter "
	   << aperture_diameter_meters
	   << " passed to this function\n";
      throw(string("refractive_atmospheric_model::phase_covariance"));
    }

    if((pupil_location_one-this->ground_ref_frame_).length()>.5*aperture_diameter_meters){
      cerr << "refractive_atmospheric_model::phase_covariance error - \n"
	   << " first point outside of the pupil\n";
      pupil_location_one.print(cerr, "pupil location one ");
      throw(string("refractive_atmospheric_model::phase_covariance"));
    }

    if((pupil_location_two-this->ground_ref_frame_).length()>.5*aperture_diameter_meters){
      cerr << "refractive_atmospheric_model::phase_covariance error - \n"
	   << " first point outside of the pupil\n";
      pupil_location_two.print(cerr, "pupil location two ");
      throw(string("refractive_atmospheric_model::phase_covariance"));
    }

    double wavenumber = 2 * M_PI / wavelength_meters;
    double fac = get_Xi()*pow(aperture_diameter_meters,5/3.)*wavenumber*wavenumber;

    three_vector normalized_pupil_location_one = (pupil_location_one - ground_ref_frame_)*(2/aperture_diameter_meters);
    three_vector normalized_pupil_location_two = (pupil_location_two - ground_ref_frame_)*(2/aperture_diameter_meters);

    double val = 
      caliph_A(emtr_a, emtr_b, aperture_diameter_meters, normalized_pupil_location_one) + 
      caliph_B(emtr_a, emtr_b, aperture_diameter_meters, normalized_pupil_location_two) -
      caliph_C(emtr_a, emtr_b, aperture_diameter_meters, normalized_pupil_location_one, normalized_pupil_location_two) - 
      caliph_D(emtr_a, emtr_b, aperture_diameter_meters, nsteps_in_integration);

    val *= fac;

    if(!finite(val)){
      cerr << "refractive_atmospheric_model::phase_covariance error - value " 
	   << val 
	   << " not finite\n";
      throw(string("refractive_atmospheric_model::phase_covariance"));
    }
    return(val);
  } 

  double refractive_atmospheric_model::
  tilt_phase_covariance(const emitter & emtr_a,
			const three_point & pupil_location_one,
			const emitter & emtr_b,
			const three_point & pupil_location_two,
			const aperture & ap,
			double wavelength_meters,
			int nsteps_in_integration) const {


    try{check_wavelength(wavelength_meters);}
    catch(...){
      cerr << "refractive_atmospheric_model::tilt_phase_covariance error\n";
      throw(string("refractive_atmospheric_model::tilt_phase_covariance"));
    }
    if(nsteps_in_integration<=0){
      cerr << "refractive_atmospheric_model::tilt_phase_covariance error\n"
	   << "number of samples in integration " 
	   << nsteps_in_integration
	   << " out of range\n";
      throw(string("refractive_atmospheric_model::tilt_phase_covariance"));
    }

    circular_aperture circ_ap;
    try{circ_ap=dynamic_cast<const circular_aperture &>(ap);}
    catch(...){
      cerr << "refractive_atmospheric_model::tilt_phase_covariance error - \n"
	   << "an aperture was passed to this function that is not a circular_aperture\n";
	throw(string("refractive_atmospheric_model::tilt_phase_covariance"));
    }
    double aperture_diameter_meters = circ_ap.get_diameter();

    if(aperture_diameter_meters<=0){
      cerr << "refractive_atmospheric_model::tilt_phase_covariance error - \n"
	   << "invalid aperture diameter "
	   << aperture_diameter_meters
	   << " passed to this function\n";
      throw(string("refractive_atmospheric_model::tilt_phase_covariance"));
    }

    if((pupil_location_one-this->ground_ref_frame_).length()>.5*aperture_diameter_meters){
      cerr << "refractive_atmospheric_model::phase_covariance error - \n"
	   << " first point outside of the pupil\n";
      pupil_location_one.print(cerr, "pupil location one ");
      throw(string("refractive_atmospheric_model::tilt_phase_covariance"));
    }

    if((pupil_location_two-this->ground_ref_frame_).length()>.5*aperture_diameter_meters){
      cerr << "refractive_atmospheric_model::phase_covariance error - \n"
	   << " first point outside of the pupil\n";
      pupil_location_two.print(cerr, "pupil location two ");
      throw(string("refractive_atmospheric_model::tilt_phase_covariance"));
    }

    double secant_zenith_angle, max_range_meters, min_range_meters;
    three_vector omega;
    Tyler_get_constants(emtr_a,
		  emtr_b,
		  this->ground_ref_frame_,
		  aperture_diameter_meters,
		  secant_zenith_angle,
		  max_range_meters,
		  min_range_meters,
		  omega);

    if(omega.length()>three_frame::precision)
      omega *= (1/omega.length());

    double wavenumber = 2 * M_PI / wavelength_meters;
    double fac = get_Xi()*pow(aperture_diameter_meters,5/3.)*wavenumber*wavenumber;

    three_vector normalized_pupil_location_one = 
      (pupil_location_one - ground_ref_frame_)*(2/aperture_diameter_meters);
    three_vector normalized_pupil_location_two = 
      (pupil_location_two - ground_ref_frame_)*(2/aperture_diameter_meters);

    double val = 0;

    double rho_1_rho_2_dot_product = dot_product(normalized_pupil_location_one,
						 normalized_pupil_location_two);

    double rho_1_omega_dot_product = dot_product(normalized_pupil_location_one,
						 omega);

    double rho_2_omega_dot_product = dot_product(normalized_pupil_location_two,
						 omega);

    double cross_prod_term = dot_product(cross_product(normalized_pupil_location_one,omega),
					 cross_product(normalized_pupil_location_two,omega));

    if(rho_1_rho_2_dot_product!=0){
      val += rho_1_rho_2_dot_product*
	(caliph_E(emtr_a, emtr_b, aperture_diameter_meters, nsteps_in_integration) - 
	 caliph_G(emtr_a, emtr_b, aperture_diameter_meters, normalized_pupil_location_two) - 
	 caliph_I(emtr_a, emtr_b, aperture_diameter_meters, normalized_pupil_location_one));
    }

    if(rho_1_omega_dot_product!=0 && rho_2_omega_dot_product!=0){
      val += -rho_1_omega_dot_product*rho_2_omega_dot_product*
	caliph_F(emtr_a, emtr_b, aperture_diameter_meters, nsteps_in_integration)
      + (rho_1_omega_dot_product*rho_2_omega_dot_product -
	 cross_prod_term)*
	caliph_F_bar(emtr_a, emtr_b, aperture_diameter_meters, nsteps_in_integration);      
    }

    if(rho_1_omega_dot_product!=0){
      val += rho_1_omega_dot_product*
	(caliph_H(emtr_a, emtr_b, aperture_diameter_meters, normalized_pupil_location_two)-
	 caliph_K(emtr_a, emtr_b, aperture_diameter_meters, nsteps_in_integration));
    }

    if(rho_2_omega_dot_product!=0){
      val -= rho_2_omega_dot_product*
	(caliph_J(emtr_a, emtr_b, aperture_diameter_meters, normalized_pupil_location_one)-
	 caliph_L(emtr_a, emtr_b, aperture_diameter_meters, nsteps_in_integration));
    }

    if(!finite(val)){
      cerr << "refractive_atmospheric_model::tilt_phase_covariance error - value " 
	   << val 
	   << " not finite\n";
      throw(string("refractive_atmospheric_model::tilt_phase_covariance"));
    }
	   
    return(fac*val);

  } 


  /*

  THIS CLASS IS REALLY AN OPTICAL SYSTEM, AND REQUIRES ADDITIONAL WORK BEFORE IT WILL BE USEFUL

  namespace {
  const keyval_map & get_refractive_atmosphere_keyval_map(){
  static keyval_map * kvm = new keyval_map;
  kvm->insert(pair<const string, string>("TYPE", "refractive atmosphere"));
  return *kvm;
  }
    
  AO_sim_base * create_refractive_atmosphere(const iofits & iof) {
  return new refractive_atmosphere(iof);
  }

  const bool refractive_atmosphere::factory_registration = 
  fits_factory<AO_sim_base>::Register(get_refractive_atmosphere_keyval_map(), create_refractive_atmosphere);
  } 

  int refractive_atmosphere::verbose_level = 0;

  refractive_atmosphere::
  refractive_atmosphere(const refractive_atmosphere & ref_atm) {
    this->operator=(ref_atm);
  }

  refractive_atmosphere::refractive_atmosphere(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "refractive_atmosphere::refractive_atmosphere - "
	   << "error opening file " << filename << endl;
      throw(string("refractive_atmosphere::refractive_atmosphere"));
    }
    this->read(iof);
  } 

  refractive_atmosphere::refractive_atmosphere(const iofits & iof) {
    this->read(iof);
  } 

  refractive_atmosphere::refractive_atmosphere(const refractive_atmospheric_model & ref_atm_model, 
					       double pixel_scale, 
					       const vector<long> & axes,
					       const subharmonic_method & subm){

    cerr << "refractive_atmosphere::refractive_atmosphere - not yet coded\n";
    throw(string("refractive_atmosphere::refractive_atmosphere"));


    //for(int i=0; i<ref_atm_model.size(); i++)
    //ref_atm_layers.push_back(ref_atm_model.get_power_spectrum(i)->get_refractive_atmospheric_layer(axes, pixel_scale, subm));

    // do whatever you have to do to fix up the three_frame using the heights

  }

  refractive_atmosphere::
  refractive_atmosphere(const vector<refractive_atmospheric_layer> & in_ref_atm_layers) {
    ref_atm_layers = in_ref_atm_layers;
  }

  refractive_atmosphere & refractive_atmosphere::
  operator=(const refractive_atmosphere & ref_atm){
    if(this==&ref_atm)
      return(*this);
    ref_atm_layers = ref_atm.ref_atm_layers;
    return(*this);
  }

  void refractive_atmosphere::read(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "refractive_atmosphere::read - "
	   << "error opening file " << filename << endl;
      throw(string("refractive_atmosphere::read"));
    }
    try{this->read(iof);}
    catch(...){
    cerr << "refractive_atmosphere::read - "
    << "error reading "	
	   << this->unique_name() << " from file "
	   << filename << endl;
	   throw(string("refractive_atmosphere::read"));
    }
  }

  void refractive_atmosphere::read(const iofits & iof){

    if(!iof.key_exists("TYPE")){
      cerr << "refractive_atmosphere::read error - "
	   << "unrecognized type of file\n";
      throw(string("refractive_atmosphere::read"));
    }
    string type, comment;
    iof.read_key("TYPE", type, comment);
    if(type!="refractive atmosphere"){
      cerr << "refractive_atmosphere::read error - file of type " 
	   << type << " rather than of type refractive atmosphere\n";
      throw(string("refractive_atmosphere::read"));
    }

    long nlayers;
    iof.read_key("NLAYERS", nlayers, comment);
    iof.movrel_hdu(1);

    ref_atm_layers.resize(nlayers);
    for(int i=0; i<nlayers; i++)
      ref_atm_layers[i].read(iof);
  }

  void refractive_atmosphere::write(const char * filename) const {
    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "refractive_atmosphere::write - "
	   << "error opening file " << filename << endl;
      throw(string("refractive_atmosphere::write"));
    }
    try{this->write(iof);}
    catch(...){
      cerr << "refractive_atmosphere::write - "
	   << "error writing "
	   << this->unique_name() << " to file " 
	   << filename << endl;
      throw(string("refractive_atmosphere::write"));
    }
  }

  void refractive_atmosphere::write(iofits & iof) const {

    fits_header_data tmphdr(iofits::BYTEIMG, vector<long>());
    tmphdr.write(iof);

    string type = "refractive atmosphere";
    string comment = "object type";
    iof.write_key("TYPE", type, comment);
    comment = "number of atmospheric layers";
    iof.write_key("NLAYERS", (long)ref_atm_layers.size(), comment);
    for(int i=0; i<ref_atm_layers.size(); i++)
      ref_atm_layers[i].write(iof);
  }

  void refractive_atmosphere::
  print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << "refractive atmosphere"
       << "/" << "object type" << endl;
    os << prefix << "NLAYERS    = " << setw(vlspc) << ref_atm_layers.size()
       << "/" << "number of atmospheric layers" << endl << endl;
  
    for(int i=0; i<ref_atm_layers.size(); i++){
      ref_atm_layers[i].print(os, prefix);
      os << endl;
    }
  }
  */

}  
