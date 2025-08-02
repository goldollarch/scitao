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
#include "NGAO_system_design_model.h"

using namespace std;

namespace Arroyo {

  namespace factory_register {

    const fits_keyval_set & get_NGAO_system_design_model_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "NGAO system design model"));
      return *fkvs;
    }
    
    AO_sim_base * create_NGAO_system_design_model(const iofits & iof) {
      return new NGAO_system_design_model(iof);
    }
  }

  const bool NGAO_system_design_model::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_NGAO_system_design_model_keyval_set(), 
				      factory_register::create_NGAO_system_design_model);


  NGAO_system_design_model::NGAO_system_design_model(const NGAO_system_design_model & cn2_model){
    this->operator=(cn2_model);
  }

  NGAO_system_design_model::NGAO_system_design_model(const char * filename){
    this->read(filename);
  }

  NGAO_system_design_model::NGAO_system_design_model(const iofits & iof){
    this->read(iof);
  }
  
  NGAO_system_design_model::NGAO_system_design_model(const three_frame & ground_ref_frame, 
						     double r_0_meters, 
						     double r_0_ref_wavelength_meters){

    ground_ref_frame_ = ground_ref_frame;

    // Initialize layer heights
    layer_heights_.resize(7);
    layer_heights_[6] = 0;
    layer_heights_[5] = 2100;
    layer_heights_[4] = 4100;
    layer_heights_[3] = 6500;
    layer_heights_[2] = 9000;
    layer_heights_[1] = 12000;
    layer_heights_[0] = 14800;

    // initialize layer weights
    vector<double> layer_weights(7);
    layer_weights[6] = .47;
    layer_weights[5] = .18;
    layer_weights[4] = .11;
    layer_weights[3] = .09;
    layer_weights[2] = .04;
    layer_weights[1] = .09;
    layer_weights[0] = .02;

    power_spectra_.resize(7);
    double exponent = -11/3.0;
    double layer_r_0_meters;
    for(int i=0; i<7; i++){
      layer_r_0_meters = r_0_meters * pow(layer_weights[i], -3/5.0);
      power_spectra_[i] = new isotropic_power_law_spectrum< power_law, null_inner_scale>(power_law(exponent, 
												   layer_r_0_meters,
												   r_0_ref_wavelength_meters),
											 null_inner_scale());
    }
		     
  }

  NGAO_system_design_model & 
  NGAO_system_design_model::operator=(const NGAO_system_design_model & cn2_model){
    if(this==&cn2_model)
      return(*this);
    this->refractive_atmospheric_model::operator=(cn2_model);
    return(*this);
  }

  void NGAO_system_design_model::read(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "NGAO_system_design_model::read - "
	   << "error opening file " << filename << endl;
      throw(string("NGAO_system_design_model::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "NGAO_system_design_model::read - "
	   << "error reading "
	   << this->unique_name() << " from file " 
	   << filename << endl;
      throw(string("NGAO_system_design_model::read"));
    }
  }

  void NGAO_system_design_model::read(const iofits & iof){
    if(!iof.key_exists("TYPE")){
      cerr << "NGAO_system_design_model::read error - "
	   << "unrecognized type of file\n";
      throw(string("NGAO_system_design_model::read"));
    }
    string type, comment;
    iof.read_key("TYPE", type, comment);
    if(type!=this->unique_name()){
      cerr << "NGAO_system_design_model::read error - file of type " 
	   << type << " rather than of type "
	   << this->unique_name() << endl;
      throw(string("NGAO_system_design_model::read"));
    }

    this->read_common_data(iof);

  }

  void NGAO_system_design_model::write(const char * filename) const {
    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "NGAO_system_design_model::write - "
	   << "error opening file " << filename << endl;
      throw(string("NGAO_system_design_model::write"));
    }
    try{this->write(iof);}
    catch(...){
      cerr << "NGAO_system_design_model::write - "
	   << "error writing "
	   << this->unique_name() << " to file " 
	   << filename << endl;
      throw(string("NGAO_system_design_model::write"));
    }
  }

  void NGAO_system_design_model::write(iofits & iof) const {

    fits_header_data<char> tmphdr;
    tmphdr.write(iof);

    string type = this->unique_name();
    string comment = "object type";
    iof.write_key("TYPE", type, comment);

    this->write_common_data(iof);
  }

  void NGAO_system_design_model::
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
