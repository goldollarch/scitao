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

#include <cmath>
#include <iomanip>
#include "fits_factory.h"
#include "iofits.h"
#include "fits_header_data.h"
#include "AO_algo.h"
#include "wind_model.h"

using namespace std;

namespace Arroyo {

  int wind_model::verbose_level = 0;

  wind_model * wind_model::wind_model_factory(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...) {
      cerr << "wind_model::wind_model_factory - "
	   << "error opening file " << filename << endl;
      throw(string("wind_model::wind_model_factory"));
    }
    return(wind_model::wind_model_factory(iof));
  }

  wind_model * wind_model::wind_model_factory(const iofits & iof){
    AO_sim_base * aosb = AO_sim_base::AO_sim_base_factory(iof);
    wind_model * wm = dynamic_cast<wind_model *>(aosb);
    if(wm==NULL)
      throw(string("wind_model::wind_model_factory"));
     return(wm);    
  }

  wind_model * wind_model::wind_model_factory(const wind_model & wm){
    return(wm.clone());
  } 

  namespace factory_register {
    const fits_keyval_set & get_Hardy_wind_model_keyval_set(){
      static fits_keyval_set * fvks = new fits_keyval_set;
      fvks->push_back(fits_keyval_entry("TYPE", "Hardy wind model"));
      return *fvks;
    }
    
    AO_sim_base * create_Hardy_wind_model(const iofits & iof) {
      return new Hardy_wind_model(iof);
    }
  }

  const bool Hardy_wind_model::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_Hardy_wind_model_keyval_set(), 
				      factory_register::create_Hardy_wind_model);

  Hardy_wind_model::Hardy_wind_model(){
    ground_layer_wind_velocity = 0;
    tropopause_wind_velocity = 0;
    tropopause_height = 0;
    tropopause_thickness = 0;
  }

  Hardy_wind_model::Hardy_wind_model(const Hardy_wind_model & twm){
    this->operator=(twm);
  }
  
  Hardy_wind_model::Hardy_wind_model(const char * filename){
    this->read(filename);
  }

  Hardy_wind_model::Hardy_wind_model(const iofits & iof){
    this->read(iof);
  }

  Hardy_wind_model & Hardy_wind_model::operator=(const Hardy_wind_model & twm){
    if(this==&twm)
      return(*this);
    ground_layer_wind_velocity = twm.ground_layer_wind_velocity;
    tropopause_wind_velocity = twm.tropopause_wind_velocity;
    tropopause_height = twm.tropopause_height;
    tropopause_thickness = twm.tropopause_thickness;
    return(*this);
  }

  Hardy_wind_model::Hardy_wind_model(double grnd_lyr_wind_velocity,
				     double trpse_height,
				     double trpse_thickness,
				     double trpse_wind_velocity){
    if(grnd_lyr_wind_velocity<0){
      cerr << "Hardy_wind_model::Hardy_wind_model error - "
	   << "ground layer velocity " << grnd_lyr_wind_velocity
	   << " provided to constructor is less than zero\n";
      throw(string("Hardy_wind_model::Hardy_wind_model"));
    }

    if(trpse_wind_velocity<0){
      cerr << "Hardy_wind_model::Hardy_wind_model error - "
	   << "tropopause velocity " << trpse_wind_velocity
	   << " provided to constructor is less than zero\n";
      throw(string("Hardy_wind_model::Hardy_wind_model"));
    }

    if(trpse_wind_velocity>0 && trpse_height<0){
      cerr << "Hardy_wind_model::Hardy_wind_model error - "
	   << "tropopause height " << trpse_height
	   << " provided to constructor is less than zero\n";
      throw(string("Hardy_wind_model::Hardy_wind_model"));
    }

    if(trpse_wind_velocity>0 && trpse_thickness<=0){
      cerr << "Hardy_wind_model::Hardy_wind_model error - "
	   << "tropopause thickness " << trpse_thickness
	   << " provided to constructor is less than or equal to zero\n";
      throw(string("Hardy_wind_model::Hardy_wind_model"));
    }

    ground_layer_wind_velocity = grnd_lyr_wind_velocity;
    tropopause_wind_velocity = trpse_wind_velocity;
    tropopause_height = trpse_height;
    tropopause_thickness = trpse_thickness;    
  }

  void Hardy_wind_model::read(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "Hardy_wind_model::read - "
	   << "error opening file " << filename << endl;
      throw(string("Hardy_wind_model::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "Hardy_wind_model::read - "
	   << "error reading "
	   << this->unique_name() << " from file " 
	   << filename << endl;
      throw(string("Hardy_wind_model::read"));
    }
  }

  void Hardy_wind_model::read(const iofits & iof) {

    if(!iof.key_exists("TYPE")){
      cerr << "Hardy_wind_model::read error - "
	   << "unrecognized type of file\n";
      throw(string("Hardy_wind_model::read"));
    }
    string type, comment;
    iof.read_key("TYPE", type, comment);
    if(type!=this->unique_name()){
      cerr << "Hardy_wind_model::read error - file of type " 
	   << type << " rather than of type "
	   << this->unique_name() << endl;
      throw(string("Hardy_wind_model::read"));
    }
    iof.read_key("GRNDVEL", ground_layer_wind_velocity, comment);
    iof.read_key("TRPSEVEL", tropopause_wind_velocity, comment);
    iof.read_key("TRPSEHGT", tropopause_height, comment);
    iof.read_key("TRPSETHK", tropopause_thickness, comment);
    // leave iof pointing to the next header, if it exists
    if(iof.get_hdu_num()!=iof.get_num_hdus())
      iof.movrel_hdu(1);
  }

  void Hardy_wind_model::write(const char * filename) const {

    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "Hardy_wind_model::write - "
	   << "error opening file " << filename << endl;
      throw(string("Hardy_wind_model::write"));
    }

    try{this->write(iof);}
    catch(...){
      cerr << "Hardy_wind_model::write - "
	   << "error writing "
	   << this->unique_name() << " to file "
	   << filename << endl;
      throw(string("Hardy_wind_model::write"));
    }
  }

  void Hardy_wind_model::write(iofits & iof) const {
    fits_header_data<char> fhd;
    fhd.write(iof);

    string type = this->unique_name();
    string comment = "object type";
    iof.write_key("TYPE", type, comment);
    iof.write_key("GRNDVEL", ground_layer_wind_velocity, "ground layer wind velocity (meters per sec)");
    iof.write_key("TRPSEVEL", tropopause_wind_velocity, "tropopause wind velocity (meters per sec)");
    iof.write_key("TRPSEHGT", tropopause_height, "tropopause height (meters)");
    iof.write_key("TRPSETHK", tropopause_thickness, "tropopause thickness (meters)");
  }

  void Hardy_wind_model::print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << this->unique_name()
       << "/" << "object type" << endl;
    os << prefix << "GRNDVEL    = " << setw(vlspc) << ground_layer_wind_velocity
       << "/" << "ground layer wind velocity in meters per second" << endl;
    os << prefix << "TRPSEVEL   = " << setw(vlspc) << tropopause_wind_velocity
       << "/" << "tropopause wind velocity in meters per second" << endl;
    os << prefix << "TRPSEHGT   = " << setw(vlspc) << tropopause_height
       << "/" << "tropopause height in meters" << endl;
    os << prefix << "TRPSETHK   = " << setw(vlspc) << tropopause_thickness
       << "/" << "tropopause thickness in meters" << endl;
  }

  vector<three_vector> Hardy_wind_model::get_random_wind_vectors(const vector<double> & heights, const three_frame & ref_frame) const {

    vector<three_vector> wind_vectors(heights.size());

    if(heights.size()<=0){
      cerr << "Hardy_wind_model::get_random_wind_vectors error - "
	   << "number of requested heights " << heights.size() << " is not positive\n";
      throw(string("Hardy_wind_model::get_random_wind_vectors"));
    }

    for(int i=0; i<heights.size(); i++){
      if(heights[i]<0){
	cerr << "Hardy_wind_model::get_random_wind_vectors error - "
	     << "requested height " << heights[i] << " is not positive\n";
	throw(string("Hardy_wind_model::get_random_wind_vectors"));
      }
    }

    double ground_vx, ground_vy;
    double tropo_vx, tropo_vy;
    box_mueller(ground_vx, ground_vy);
    box_mueller(tropo_vx, tropo_vy);

    ground_vx *= ground_layer_wind_velocity;
    ground_vy *= ground_layer_wind_velocity;

    tropo_vx *= tropopause_wind_velocity;
    tropo_vy *= tropopause_wind_velocity;

    if(wind_model::verbose_level>=2){
      cerr << "Hardy_wind_model::get_random_wind_vectors - ground wind velocity " << ground_vx << ", " << ground_vy << endl;
      cerr << "Hardy_wind_model::get_random_wind_vectors - tropopause wind velocity " << tropo_vx << ", " << tropo_vy << endl;
    }

    double tmp;
    for(int i=0; i<heights.size(); i++){
      if(tropopause_thickness>0){
	tmp = (heights[i] - tropopause_height)/tropopause_thickness;
	wind_vectors[i] = three_vector(ground_vx + tropo_vx*exp(-tmp*tmp), 
				       ground_vy + tropo_vy*exp(-tmp*tmp), 0, ref_frame);
      } else 
	wind_vectors[i] = three_vector(ground_vx, ground_vy, 0, ref_frame);
    }

    return(wind_vectors);
  }   
}
