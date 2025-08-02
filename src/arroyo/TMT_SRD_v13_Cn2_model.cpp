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
#include "TMT_SRD_v13_Cn2_model.h"

using namespace std;

namespace Arroyo {

  namespace factory_register {

    const fits_keyval_set & get_TMT_SRD_v13_Cn2_model_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "TMT SRD v13 Cn2 model"));
      return *fkvs;
    }
    
    AO_sim_base * create_TMT_SRD_v13_Cn2_model(const iofits & iof) {
      return new TMT_SRD_v13_Cn2_model(iof);
    }
  }

  const bool TMT_SRD_v13_Cn2_model::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_TMT_SRD_v13_Cn2_model_keyval_set(), 
				      factory_register::create_TMT_SRD_v13_Cn2_model);


  TMT_SRD_v13_Cn2_model::TMT_SRD_v13_Cn2_model(const TMT_SRD_v13_Cn2_model & cn2_model){
    this->operator=(cn2_model);
  }

  TMT_SRD_v13_Cn2_model::TMT_SRD_v13_Cn2_model(const char * filename){
    this->read(filename);
  }

  TMT_SRD_v13_Cn2_model::TMT_SRD_v13_Cn2_model(const iofits & iof){
    this->read(iof);
  }
  
  TMT_SRD_v13_Cn2_model::TMT_SRD_v13_Cn2_model(const three_frame & ground_ref_frame, 
							       double r_0_meters, 
							       double r_0_ref_wavelength_meters){

    ground_ref_frame_ = ground_ref_frame;

    // Initialize layer heights
    layer_heights_.resize(7);
    layer_heights_[6] = 0;
    layer_heights_[5] = 1800;
    layer_heights_[4] = 3300;
    layer_heights_[3] = 5800;
    layer_heights_[2] = 7400;
    layer_heights_[1] = 13100;
    layer_heights_[0] = 15800;

    // initialize layer weights
    vector<double> layer_weights(7);
    layer_weights[6] = .646;
    layer_weights[5] = .078;
    layer_weights[4] = .119;
    layer_weights[3] = .035;
    layer_weights[2] = .025;
    layer_weights[1] = .082;
    layer_weights[0] = .015;

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

  TMT_SRD_v13_Cn2_model & 
  TMT_SRD_v13_Cn2_model::operator=(const TMT_SRD_v13_Cn2_model & cn2_model){
    if(this==&cn2_model)
      return(*this);
    this->refractive_atmospheric_model::operator=(cn2_model);
    return(*this);
  }

  void TMT_SRD_v13_Cn2_model::read(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "TMT_SRD_v13_Cn2_model::read - "
	   << "error opening file " << filename << endl;
      throw(string("TMT_SRD_v13_Cn2_model::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "TMT_SRD_v13_Cn2_model::read - "
	   << "error reading "
	   << this->unique_name() << " from file " 
	   << filename << endl;
      throw(string("TMT_SRD_v13_Cn2_model::read"));
    }
  }

  void TMT_SRD_v13_Cn2_model::read(const iofits & iof){
    if(!iof.key_exists("TYPE")){
      cerr << "TMT_SRD_v13_Cn2_model::read error - "
	   << "unrecognized type of file\n";
      throw(string("TMT_SRD_v13_Cn2_model::read"));
    }
    string type, comment;
    iof.read_key("TYPE", type, comment);
    if(type!=this->unique_name()){
      cerr << "TMT_SRD_v13_Cn2_model::read error - file of type " 
	   << type << " rather than of type "
	   << this->unique_name() << endl;
      throw(string("TMT_SRD_v13_Cn2_model::read"));
    }

    this->read_common_data(iof);

  }

  void TMT_SRD_v13_Cn2_model::write(const char * filename) const {
    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "TMT_SRD_v13_Cn2_model::write - "
	   << "error opening file " << filename << endl;
      throw(string("TMT_SRD_v13_Cn2_model::write"));
    }
    try{this->write(iof);}
    catch(...){
      cerr << "TMT_SRD_v13_Cn2_model::write - "
	   << "error writing "
	   << this->unique_name() << " to file " 
	   << filename << endl;
      throw(string("TMT_SRD_v13_Cn2_model::write"));
    }
  }

  void TMT_SRD_v13_Cn2_model::write(iofits & iof) const {

    fits_header_data<char> tmphdr;
    tmphdr.write(iof);

    string type = this->unique_name();
    string comment = "object type";
    iof.write_key("TYPE", type, comment);

    this->write_common_data(iof);
  }

  void TMT_SRD_v13_Cn2_model::
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
