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
#include "Hufnagel_Valley_model.h"

using namespace std;

namespace Arroyo {

  namespace factory_register {

    const fits_keyval_set & get_Hufnagel_Valley_model_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "Hufnagel Valley model"));
      return *fkvs;
    }
    
    AO_sim_base * create_Hufnagel_Valley_model(const iofits & iof) {
      return new Hufnagel_Valley_model(iof);
    }
  }
  
  const bool Hufnagel_Valley_model::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_Hufnagel_Valley_model_keyval_set(), 
				      factory_register::create_Hufnagel_Valley_model);


  Hufnagel_Valley_model::Hufnagel_Valley_model(const Hufnagel_Valley_model & hv_model){
    this->operator=(hv_model);
  }

  Hufnagel_Valley_model::Hufnagel_Valley_model(const char * filename){
    this->read(filename);
  }

  Hufnagel_Valley_model::Hufnagel_Valley_model(const iofits & iof){
    this->read(iof);
  }
  
  namespace {
    long factorial(int n){
      if(n<0){
	cerr << "factorial error - can't return factorial of " << n << endl;
	throw(string("factorial"));
      }
      long val = 1;
      for(int i=2; i<=n; i++) val *= i;
      return(val);
    }

    // Returns the integral of the Hufnagel Valley model from the
    // ground to height, which is in meters.
    double Hufnagel_Valley_Cn2_integral(double pseudowind, 
					double A, 
					double height){
      if(height<0){
	cerr << "Hufnagel_Valley_Cn2_integral error - height " << height 
	     << " is invalid\n";
	throw(string("Hufnagel_Valley_Cn2_integral"));
      }
      
      vector<double> c(3);
      c[0] = .00594*pseudowind*pseudowind/27./27.;
      c[1] = 2.7e-16;
      c[2] = A;

      vector<double> e(3);
      e[0] = 1e-2;
      e[1] = 1500;
      e[2] = 100;

      double ten_factorial = factorial(10);
      double Cn2_integral_value = c[0]*ten_factorial * pow(e[0],11.) * 1e5;
      double tmp = 0;
      for(int i=0; i<=10; i++)
	tmp += pow(1e-5*height, i)*pow(e[0], 10-i+1)*ten_factorial/(double)factorial(i);
      Cn2_integral_value -= c[0]*tmp*exp(-1e-5*height/e[0])*1e5;
	
      Cn2_integral_value += e[1]*c[1]*(1-exp(-height/e[1]));
      Cn2_integral_value += e[2]*c[2]*(1-exp(-height/e[2]));

      return(Cn2_integral_value);
    }
  }


  Hufnagel_Valley_model::Hufnagel_Valley_model(const three_frame & ground_ref_frame, 
					       const std::vector<double> & layer_heights,
					       double pseudowind, double A){

    if(layer_heights.size()==0){
      cerr << "Hufnagel_Valley_model::Hufnagel_Valley_model error - "
	   << "vector of layer altitudes supplied to this constructor has zero size\n";
      throw(string("Hufnagel_Valley_model::Hufnagel_Valley_model"));
    }

    for(int i=0; i<layer_heights.size(); i++){
      if(layer_heights[i]<0){
	cerr << "Hufnagel_Valley_model::Hufnagel_Valley_model error - "
	     << "layer altitude " << i << " has negative value " 
	     << layer_heights[i] << endl;
	throw(string("Hufnagel_Valley_model::Hufnagel_Valley_model"));
      }
    }

    ground_ref_frame_ = ground_ref_frame;

    // Initialize layer heights
    layer_heights_ = layer_heights;

    sort(layer_heights_.begin(),layer_heights_.end()-1);

    // sort the layers
    sort(layer_heights_.begin(),layer_heights_.end());

    int nlayers = layer_heights_.size();
    double integrated_cn2_profile;
    vector<double> power_law_coefficients(nlayers);
    power_law_coefficients[0] = 
      2*M_PI*.033*(Hufnagel_Valley_Cn2_integral(pseudowind, 
						A, 
						.5*(layer_heights_[0]+layer_heights_[1])));

    for(int i=1; i<nlayers-1; i++)
	power_law_coefficients[i] = 	   
	  2*M_PI*.033*(Hufnagel_Valley_Cn2_integral(pseudowind, 
						    A, 
						    .5*(layer_heights_[i]+layer_heights_[i+1]))) -
	  power_law_coefficients[i-1];

    power_law_coefficients.back() = 
      2*M_PI*.033*(Hufnagel_Valley_Cn2_integral(pseudowind, 
						A, 
						1e6)) - 
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
	     
    for(int i=0; i<nlayers; i++)
      cout << i << "\t" << layer_heights_[i] << "\t" << power_law_coefficients[i] << endl;
  }

  Hufnagel_Valley_model & 
  Hufnagel_Valley_model::operator=(const Hufnagel_Valley_model & hv_model){
    if(this==&hv_model)
      return(*this);
    this->refractive_atmospheric_model::operator=(hv_model);
    return(*this);
  }

  void Hufnagel_Valley_model::read(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "Hufnagel_Valley_model::read - "
	   << "error opening file " << filename << endl;
      throw(string("Hufnagel_Valley_model::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "Hufnagel_Valley_model::read - "
	   << "error reading "
	   << this->unique_name() << " from file " 
	   << filename << endl;
      throw(string("Hufnagel_Valley_model::read"));
    }
  }

  void Hufnagel_Valley_model::read(const iofits & iof){
    if(!iof.key_exists("TYPE")){
      cerr << "Hufnagel_Valley_model::read error - "
	   << "unrecognized type of file\n";
      throw(string("Hufnagel_Valley_model::read"));
    }
    string type, comment;
    iof.read_key("TYPE", type, comment);
    if(type!=this->unique_name()){
      cerr << "Hufnagel_Valley_model::read error - file of type " 
	   << type << " rather than of type "
	   << this->unique_name() << endl;
      throw(string("Hufnagel_Valley_model::read"));
    }
    this->read_common_data(iof);
  }

  void Hufnagel_Valley_model::write(const char * filename) const {
    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "Hufnagel_Valley_model::write - "
	   << "error opening file " << filename << endl;
      throw(string("Hufnagel_Valley_model::write"));
    }
    try{this->write(iof);}
    catch(...){
      cerr << "Hufnagel_Valley_model::write - "
	   << "error writing "
	   << this->unique_name() << " to file " 
	   << filename << endl;
      throw(string("Hufnagel_Valley_model::write"));
    }
  }

  void Hufnagel_Valley_model::write(iofits & iof) const {

    fits_header_data<char> tmphdr;
    tmphdr.write(iof);

    string type = this->unique_name();
    string comment = "object type";
    iof.write_key("TYPE", type, comment);

    this->write_common_data(iof);
  }

  void Hufnagel_Valley_model::
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
