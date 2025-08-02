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

#include "fits_factory.h"
#include "SLCSAT_model.h"

using namespace std;

namespace Arroyo {

  namespace factory_register {

    const fits_keyval_set & get_SLCSAT_day_model_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "SLCSAT day model"));
      return *fkvs;
    }
    
    AO_sim_base * create_SLCSAT_day_model(const iofits & iof) {
      return new SLCSAT_day_model(iof);
    }
  }

  const bool SLCSAT_day_model::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_SLCSAT_day_model_keyval_set(), 
				      factory_register::create_SLCSAT_day_model);


  SLCSAT_day_model::SLCSAT_day_model(const SLCSAT_day_model & slcday_model){
    this->operator=(slcday_model);
  }

  SLCSAT_day_model::SLCSAT_day_model(const char * filename){
    this->read(filename);
  }

  SLCSAT_day_model::SLCSAT_day_model(const iofits & iof){
    this->read(iof);
  }

  namespace {
    // Returns the integral of the SLCSAT day model from 
    // the ground to height, which is in meters.
    double SLCSAT_day_Cn2_integral(double height){
      if(height<0){
	cerr << "SLCSAT_day_Cn2_integral error - height " << height 
	     << " is invalid\n";
	throw(string("SLCSAT_day_Cn2_integral"));
      }

      vector<double> h(5);
      h[0] = 18.5;
      h[1] = 232;
      h[2] = 880;
      h[3] = 7220;
      h[4] = 20500;
      vector<double> c(4);
      c[0] = 3.96e-13;
      c[1] = 1.3e-15;
      c[2] = 8.87e-7;
      c[3] = 2.0e-16;

      double Cn2_integral_value = 0;
      if(height<h[0]) return(0);

      if(height<h[1]) return((c[0]/.05)*(pow(h[0],-.05) - pow(height,-.05)));
      else Cn2_integral_value += (c[0]/.05)*(pow(h[0],-.05) - pow(h[1],-.05));

      if(height<h[2]) return(Cn2_integral_value + c[1]*(height-h[1]));
      else Cn2_integral_value += c[1]*(h[2]-h[1]);

      if(height<h[3]) return(Cn2_integral_value + .5*c[2]*(1/(h[2]*h[2]) - 1/(height*height)));
      else Cn2_integral_value += .5*c[2]*(1/(h[2]*h[2]) - 1/(height*height));
			    
      if(height<h[4]) return(Cn2_integral_value + 2*c[3]*(sqrt(height) - sqrt(h[3])));
      else return(Cn2_integral_value + 2*c[3]*(sqrt(h[4]) - sqrt(h[3])));
			    
    }
  }
  
  SLCSAT_day_model::SLCSAT_day_model(const three_frame & ground_ref_frame,
				     const std::vector<double> & layer_heights){

    if(layer_heights.size()==0){
      cerr << "SLCSAT_day_model::SLCSAT_day_model error - "
	   << "vector of layer altitudes supplied to this constructor has zero size\n";
      throw(string("SLCSAT_day_model::SLCSAT_day_model"));
    }

    for(int i=0; i<layer_heights.size(); i++){
      if(layer_heights[i]<0){
	cerr << "SLCSAT_day_model::SLCSAT_day_model error - "
	     << "layer altitude " << i << " has negative value " 
	     << layer_heights[i] << endl;
	throw(string("SLCSAT_day_model::SLCSAT_day_model"));
      }
    }

    ground_ref_frame_ = ground_ref_frame;

    // Initialize layer heights
    layer_heights_ = layer_heights;

    // Sort the layers
    sort(layer_heights_.begin(),layer_heights_.end());

    // Define the power law coefficients
    int nlayers = layer_heights_.size();
    double integrated_cn2_profile;
    vector<double> power_law_coefficients(nlayers);
    power_law_coefficients[0] = 
      2*M_PI*.033*(SLCSAT_day_Cn2_integral(.5*(layer_heights_[0]+layer_heights_[1])));

    for(int i=1; i<nlayers-1; i++)
	power_law_coefficients[i] = 	   
	  2*M_PI*.033*(SLCSAT_day_Cn2_integral(.5*(layer_heights_[i]+layer_heights_[i+1]))) -
	  power_law_coefficients[i-1];

    power_law_coefficients.back() = 
      2*M_PI*.033*(SLCSAT_day_Cn2_integral(DBL_MAX)) - 
      power_law_coefficients[layer_heights_.size()-2];


    // Reverse sort the layers and coefficients,
    // as this is the convention set by previous
    // models    
    sort(layer_heights_.rbegin(),layer_heights_.rend());
    vector<double> tmp_coeffs = power_law_coefficients;
    for(int i=0; i<nlayers; i++){
      power_law_coefficients[nlayers-i-1] = tmp_coeffs[i];
    }

    // Construct the power spectra
    power_spectra_.resize(layer_heights_.size());
    double exponent = -11/3.0;
    for(int i=0; i<nlayers; i++){
      power_law plaw(exponent, power_law_coefficients[i]);
      power_spectra_[i] = new isotropic_power_law_spectrum<power_law, null_inner_scale>(plaw, null_inner_scale());
    }
	     
  }

  SLCSAT_day_model & 
  SLCSAT_day_model::operator=(const SLCSAT_day_model & slcday_model){
    if(this==&slcday_model)
      return(*this);
    this->refractive_atmospheric_model::operator=(slcday_model);
    return(*this);
  }

  void SLCSAT_day_model::read(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "SLCSAT_day_model::read - "
	   << "error opening file " << filename << endl;
      throw(string("SLCSAT_day_model::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "SLCSAT_day_model::read - "
	   << "error reading "
	   << this->unique_name() << " from file " 
	   << filename << endl;
      throw(string("SLCSAT_day_model::read"));
    }
  }

  void SLCSAT_day_model::read(const iofits & iof){
    if(!iof.key_exists("TYPE")){
      cerr << "SLCSAT_day_model::read error - "
	   << "unrecognized type of file\n";
      throw(string("SLCSAT_day_model::read"));
    }
    string type, comment;
    iof.read_key("TYPE", type, comment);
    if(type!=this->unique_name()){
      cerr << "SLCSAT_day_model::read error - file of type " 
	   << type << " rather than of type "
	   << this->unique_name() << endl;
      throw(string("SLCSAT_day_model::read"));
    }
    this->read_common_data(iof);
  }

  void SLCSAT_day_model::write(const char * filename) const {
    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "SLCSAT_day_model::write - "
	   << "error opening file " << filename << endl;
      throw(string("SLCSAT_day_model::write"));
    }
    try{this->write(iof);}
    catch(...){
      cerr << "SLCSAT_day_model::write - "
	   << "error writing "
	   << this->unique_name() << " to file " 
	   << filename << endl;
      throw(string("SLCSAT_day_model::write"));
    }
  }

  void SLCSAT_day_model::write(iofits & iof) const {

    fits_header_data<char> tmphdr;
    tmphdr.write(iof);

    string type = this->unique_name();
    string comment = "object type";
    iof.write_key("TYPE", type, comment);

    this->write_common_data(iof);
  }

  void SLCSAT_day_model::
  print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << this->unique_name()
       << "/" << "object type" << endl;
    os << prefix << "NPSPEC     = " << setw(vlspc) << power_spectra_.size()
       << "/" << "number of power spectra" << endl;
    ground_ref_frame_.print(os, prefix);

    for(int i=0; i<power_spectra_.size(); i++){
      power_spectra_[i]->print(os, prefix);
      os << prefix << "HEIGHT     = " << setw(vlspc) << layer_heights_[i]
	 << "/" << "height of layer (meters)" << endl << endl;
    }
  }

  namespace factory_register {

    const fits_keyval_set & get_SLCSAT_night_model_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "SLCSAT night model"));
      return *fkvs;
    }
    
    AO_sim_base * create_SLCSAT_night_model(const iofits & iof) {
      return new SLCSAT_night_model(iof);
    }
  }

  const bool SLCSAT_night_model::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_SLCSAT_night_model_keyval_set(), 
				      factory_register::create_SLCSAT_night_model);


  SLCSAT_night_model::SLCSAT_night_model(const SLCSAT_night_model & slcnight_model){
    this->operator=(slcnight_model);
  }

  SLCSAT_night_model::SLCSAT_night_model(const char * filename){
    this->read(filename);
  }

  SLCSAT_night_model::SLCSAT_night_model(const iofits & iof){
    this->read(iof);
  }
  
  namespace {
    // Returns the integral of the SLCSAT night model from 
    // the ground to height, which is in meters.
    double SLCSAT_night_Cn2_integral(double height){
      if(height<0){
	cerr << "SLCSAT_night_Cn2_integral error - height " << height 
	     << " is invalid\n";
	throw(string("SLCSAT_night_Cn2_integral"));
      }
      
      vector<double> h(5);
      h[0] = 18.5;
      h[1] = 110;
      h[2] = 850;
      h[3] = 7000;
      h[4] = 20500;
      vector<double> c(5);
      c[0] = 5.0e-15;
      c[1] = 2.875e-12;
      c[2] = 2.5e-16;
      c[3] = 8.87e-7;
      c[4] = 2.0e-16;
      
      double Cn2_integral_value = 0;
      if(height<h[0]) return(c[0]*height);
      else Cn2_integral_value += c[0]*h[0];
      
      if(height<h[1]) return(c[1]*(1/h[0] - 1/height));
      else Cn2_integral_value += c[1]*(1/h[0] - 1/h[1]);
      
      if(height<h[2]) return(Cn2_integral_value + c[2]*(height-h[1]));
      else Cn2_integral_value += c[2]*(h[2]-h[1]);
      
      if(height<h[3]) return(Cn2_integral_value + .5*c[3]*(1/(h[2]*h[2]) - 1/(height*height)));
      else Cn2_integral_value += .5*c[3]*(1/(h[2]*h[2]) - 1/(h[3]*h[3]));
      
      if(height<h[4]) return(Cn2_integral_value + 2*c[4]*(sqrt(height) - sqrt(h[3])));
      else return(Cn2_integral_value + 2*c[4]*(sqrt(h[4]) - sqrt(h[3])));
      
    }
  }
  
  SLCSAT_night_model::SLCSAT_night_model(const three_frame & ground_ref_frame,
					 const std::vector<double> & layer_heights){

    if(layer_heights.size()==0){
      cerr << "SLCSAT_night_model::SLCSAT_night_model error - "
	   << "vector of layer altitudes supplied to this constructor has zero size\n";
      throw(string("SLCSAT_night_model::SLCSAT_night_model"));
    }

    for(int i=0; i<layer_heights.size(); i++){
      if(layer_heights[i]<0){
	cerr << "SLCSAT_night_model::SLCSAT_night_model error - "
	     << "layer altitude " << i << " has negative value " 
	     << layer_heights[i] << endl;
	throw(string("SLCSAT_night_model::SLCSAT_night_model"));
      }
    }

    ground_ref_frame_ = ground_ref_frame;

    // Initialize layer heights
    layer_heights_ = layer_heights;

    // sort the layers
    sort(layer_heights_.begin(),layer_heights_.end());

    int nlayers = layer_heights_.size();
    double integrated_cn2_profile;
    vector<double> power_law_coefficients(nlayers);
    power_law_coefficients[0] = 
      2*M_PI*.033*(SLCSAT_night_Cn2_integral(.5*(layer_heights_[0]+layer_heights_[1])));

    for(int i=1; i<nlayers-1; i++)
	power_law_coefficients[i] = 	   
	  2*M_PI*.033*(SLCSAT_night_Cn2_integral(.5*(layer_heights_[i]+layer_heights_[i+1]))) -
	  power_law_coefficients[i-1];

    power_law_coefficients.back() = 
      2*M_PI*.033*(SLCSAT_night_Cn2_integral(DBL_MAX)) - 
      power_law_coefficients[nlayers-2];


    // Reverse sort the layers and coefficients,
    // as this is the convention set by previous
    // models    
    sort(layer_heights_.rbegin(),layer_heights_.rend());
    vector<double> tmp_coeffs = power_law_coefficients;
    for(int i=0; i<nlayers; i++){
      power_law_coefficients[nlayers-i-1] = tmp_coeffs[i];
    }

    power_spectra_.resize(nlayers);
    double exponent = -11/3.0;
    for(int i=0; i<nlayers; i++){
      power_law plaw(exponent, power_law_coefficients[i]);
      power_spectra_[i] = new isotropic_power_law_spectrum<power_law, null_inner_scale>(plaw, null_inner_scale());
    }
	     
  }

  SLCSAT_night_model & 
  SLCSAT_night_model::operator=(const SLCSAT_night_model & slcnight_model){
    if(this==&slcnight_model)
      return(*this);
    this->refractive_atmospheric_model::operator=(slcnight_model);
    return(*this);
  }

  void SLCSAT_night_model::read(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "SLCSAT_night_model::read - "
	   << "error opening file " << filename << endl;
      throw(string("SLCSAT_night_model::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "SLCSAT_night_model::read - "
	   << "error reading "
	   << this->unique_name() << " from file " 
	   << filename << endl;
      throw(string("SLCSAT_night_model::read"));
    }
  }

  void SLCSAT_night_model::read(const iofits & iof){
    if(!iof.key_exists("TYPE")){
      cerr << "SLCSAT_night_model::read error - "
	   << "unrecognized type of file\n";
      throw(string("SLCSAT_night_model::read"));
    }
    string type, comment;
    iof.read_key("TYPE", type, comment);
    if(type!=this->unique_name()){
      cerr << "SLCSAT_night_model::read error - file of type " 
	   << type << " rather than of type "
	   << this->unique_name() << endl;
      throw(string("SLCSAT_night_model::read"));
    }
    this->read_common_data(iof);
  }

  void SLCSAT_night_model::write(const char * filename) const {
    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "SLCSAT_night_model::write - "
	   << "error opening file " << filename << endl;
      throw(string("SLCSAT_night_model::write"));
    }
    try{this->write(iof);}
    catch(...){
      cerr << "SLCSAT_night_model::write - "
	   << "error writing "
	   << this->unique_name() << " to file " 
	   << filename << endl;
      throw(string("SLCSAT_night_model::write"));
    }
  }

  void SLCSAT_night_model::write(iofits & iof) const {

    fits_header_data<char> tmphdr;
    tmphdr.write(iof);

    string type = this->unique_name();
    string comment = "object type";
    iof.write_key("TYPE", type, comment);

    this->write_common_data(iof);
  }

  void SLCSAT_night_model::
  print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << this->unique_name()
       << "/" << "object type" << endl;
    os << prefix << "NPSPEC     = " << setw(vlspc) << power_spectra_.size()
       << "/" << "number of power spectra" << endl;
    ground_ref_frame_.print(os, prefix);

    for(int i=0; i<power_spectra_.size(); i++){
      power_spectra_[i]->print(os, prefix);
      os << prefix << "HEIGHT     = " << setw(vlspc) << layer_heights_[i]
	 << "/" << "height of layer (meters)" << endl << endl;
    }
  }

}
