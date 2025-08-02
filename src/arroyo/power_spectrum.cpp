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
#include "AO_algo.h"
#include "fft_manager.h"
#include "iofits.h"
#include "fits_factory.h"
#include "pixel_array.h"
#include "special_functions.h"
#include "power_spectrum.h"
#include "subharmonic_method.h"
#include "structure_function.h"
#include "refractive_atmospheric_layer.h"
#include "diffractive_wavefront.h"

using namespace std;

namespace Arroyo {

  inner_scale * 
  inner_scale::inner_scale_factory(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...) {
      cerr << "inner_scale::inner_scale_factory - "
	   << "error opening file " << filename << endl;
      throw(string("inner_scale::inner_scale_factory"));
    }
    return(inner_scale::inner_scale_factory(iof));
  }

  inner_scale * 
  inner_scale::inner_scale_factory(const iofits & iof){
    AO_sim_base * aosb = AO_sim_base::AO_sim_base_factory(iof);
    inner_scale * inscle = dynamic_cast<inner_scale *>(aosb);
    if(inscle==NULL)
      throw(string("inner_scale::inner_scale_factory"));
    return(inscle);
  }

  namespace factory_register {
    const fits_keyval_set & get_null_inner_scale_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "null inner scale"));
      return *fkvs;
    }
    
    AO_sim_base * create_null_inner_scale(const iofits & iof) {
      return new null_inner_scale(iof);
    }
  }

  const bool null_inner_scale::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_null_inner_scale_keyval_set(), 
				      factory_register::create_null_inner_scale);

  null_inner_scale::
  null_inner_scale(const null_inner_scale & null_inscle){
    this->operator=(null_inscle);
  }

  null_inner_scale::null_inner_scale(const char * filename){
    this->read(filename);
  }

  null_inner_scale::null_inner_scale(const iofits & iof){
    this->read(iof);
  }

  null_inner_scale & null_inner_scale::
  operator=(const null_inner_scale & null_inscle){
    if(this==&null_inscle)
      return(*this);

    return(*this);
  }

  void null_inner_scale::read(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "null_inner_scale::read - "
	   << "error opening file " << filename << endl;
      throw(string("null_inner_scale::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "null_inner_scale::read - "
	   << "error reading "
	   << this->unique_name() << " from file "
	   << filename << endl;
      throw(string("null_inner_scale::read"));
    }
  }

  void null_inner_scale::read(const iofits & iof){

    if(!iof.key_exists("TYPE")){
      cerr << "null_inner_scale::read error - "
	   << "unrecognized type\n";
      throw(string("null_inner_scale::read"));
    } else { 
      string type, comment;
      iof.read_key("TYPE", type, comment);
      if(type!=this->unique_name()){
	cerr << "null_inner_scale::read error - file of type " 
	     << type << " rather than type " 
	     << this->unique_name() << endl;
	throw(string("null_inner_scale::read"));
      } 
    }
  }

  void null_inner_scale::write(const char * filename) const {
    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "null_inner_scale::write - "
	   << "error opening file " << filename << endl;
      throw(string("null_inner_scale::write"));
    }

    try{this->write(iof);}
    catch(...){
      cerr << "null_inner_scale::write - "
	   << "error writing "
	   << this->unique_name() << " to file " 
	   << filename << endl;
      throw(string("null_inner_scale::write"));
    }
  }

  void null_inner_scale::write(iofits & iof) const {
    fits_header_data<char> tmphdr;
    tmphdr.write(iof);
    string comment = "type";
    string type = this->unique_name();
    iof.write_key("TYPE", type, comment); 
  }

  void null_inner_scale::print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << this->unique_name()
       << "/" << "object type" << endl;
  }

  bool operator==(const null_inner_scale & eis1, const null_inner_scale & eis2) {
    return true;
  }

  bool operator!=(const null_inner_scale & eis1, const null_inner_scale & eis2) {
    return !operator==(eis1, eis2);    
  }

  namespace factory_register {
    const fits_keyval_set & get_exponential_inner_scale_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "exponential inner scale"));
      return *fkvs;
    }
    
    AO_sim_base * create_exponential_inner_scale(const iofits & iof) {
      return new exponential_inner_scale(iof);
    }
  }

  const bool exponential_inner_scale::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_exponential_inner_scale_keyval_set(), 
				      factory_register::create_exponential_inner_scale);

  exponential_inner_scale::exponential_inner_scale(){
    inner_scale_value = 0;
  }

  exponential_inner_scale::
  exponential_inner_scale(const exponential_inner_scale & exponential_inscle){
    this->operator=(exponential_inscle);
  }

  exponential_inner_scale::exponential_inner_scale(const char * filename){
    this->read(filename);
  }

  exponential_inner_scale::exponential_inner_scale(const iofits & iof){
    this->read(iof);
  }

  exponential_inner_scale::exponential_inner_scale(double inscle){
    if(inscle <= 0){
      cerr << "exponential_inner_scale::exponential_inner_scale error - "
	   << "inner scale may not be less than zero\n";
      throw(string("exponential_inner_scale::exponential_inner_scale"));
    }
    inner_scale_value = inscle;
  }

  exponential_inner_scale & exponential_inner_scale::
  operator=(const exponential_inner_scale & exponential_inscle){
    if(this==&exponential_inscle)
      return(*this);
    inner_scale_value = exponential_inscle.inner_scale_value;
    return(*this);
  }

  void exponential_inner_scale::read(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "exponential_inner_scale::read - "
	   << "error opening file " << filename << endl;
      throw(string("exponential_inner_scale::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "exponential_inner_scale::read - "
	   << "error reading "
	   << this->unique_name() << " from file "
	   << filename << endl;
      throw(string("exponential_inner_scale::read"));
    }
  }

  void exponential_inner_scale::read(const iofits & iof){

    if(!iof.key_exists("TYPE")){
      cerr << "exponential_inner_scale::read error - "
	   << "unrecognized type\n";
      throw(string("exponential_inner_scale::read"));
    } else { 
      string type, comment;
      iof.read_key("TYPE", type, comment);
      if(type!=this->unique_name()){
	cerr << "exponential_inner_scale::read error - file of type " 
	     << type << " rather than type exponential inner scale\n";
	throw(string("exponential_inner_scale::read"));
      } else 
	iof.read_key("INSCALE", inner_scale_value, comment);
    }
  }

  void exponential_inner_scale::write(const char * filename) const {
    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "exponential_inner_scale::write - "
	   << "error opening file " << filename << endl;
      throw(string("exponential_inner_scale::write"));
    }

    try{this->write(iof);}
    catch(...){
      cerr << "exponential_inner_scale::write - "
	   << "error writing "
	   << this->unique_name() << " to file " 
	   << filename << endl;
      throw(string("exponential_inner_scale::write"));
    }
  }

  void exponential_inner_scale::write(iofits & iof) const {
    fits_header_data<char> tmphdr;
    tmphdr.write(iof);
    string comment = "type";
    string type = this->unique_name();
    iof.write_key("TYPE", type, comment); 
    comment = "inner scale value (meters)";
    iof.write_key("INSCALE", inner_scale_value, comment); 
  }

  void exponential_inner_scale::print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << this->unique_name()
       << "/" << "object type" << endl;
    os << prefix << "INSCALE    = " << setw(vlspc) << inner_scale_value
       << "/" << "inner scale value (meters)" << endl;
  }

  bool operator==(const exponential_inner_scale & eis1, const exponential_inner_scale & eis2) {
    if(eis1.inner_scale_value!=eis2.inner_scale_value) return false;
    return true;
  }

  bool operator!=(const exponential_inner_scale & eis1, const exponential_inner_scale & eis2) {
    return !operator==(eis1, eis2);
  }

  namespace factory_register {
    const fits_keyval_set & get_frehlich_inner_scale_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "frehlich inner scale"));
      return *fkvs;
    }
    
    AO_sim_base * create_frehlich_inner_scale(const iofits & iof) {
      return new frehlich_inner_scale(iof);
    }
  }

  const bool frehlich_inner_scale::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_frehlich_inner_scale_keyval_set(), 
				      factory_register::create_frehlich_inner_scale);

  frehlich_inner_scale::frehlich_inner_scale(){
    inner_scale_value = 0;
  }

  frehlich_inner_scale::
  frehlich_inner_scale(const frehlich_inner_scale & frehlich_inscle){
    this->operator=(frehlich_inscle);
  }

  frehlich_inner_scale::frehlich_inner_scale(const char * filename){
    this->read(filename);
  }

  frehlich_inner_scale::frehlich_inner_scale(const iofits & iof){
    this->read(iof);
  }

  frehlich_inner_scale::frehlich_inner_scale(double inscle){
    if(inscle <= 0){
      cerr << "frehlich_inner_scale::frehlich_inner_scale error - "
	   << "inner scale may not be less than zero\n";
      throw(string("frehlich_inner_scale::frehlich_inner_scale"));
    }
    inner_scale_value = inscle;
  }

  frehlich_inner_scale & 
  frehlich_inner_scale::operator=(const frehlich_inner_scale & frehlich_inscle){
    if(this==&frehlich_inscle)
      return(*this);
    inner_scale_value = frehlich_inscle.inner_scale_value;
    return(*this);
  }

  void frehlich_inner_scale::read(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "frehlich_inner_scale::read - "
	   << "error opening file " << filename << endl;
      throw(string("frehlich_inner_scale::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "frehlich_inner_scale::read - "
	   << "error reading "
	   << this->unique_name() << " from file "
	   << filename << endl;
      throw(string("frehlich_inner_scale::read"));
    }
  }

  void frehlich_inner_scale::read(const iofits & iof){
    if(!iof.key_exists("TYPE")){
      cerr << "frehlich_inner_scale::read error - "
	   << "no inner scale type specified\n";
      throw(string("frehlich_inner_scale::read"));
    } else {
      string type, comment;
      iof.read_key("TYPE", type, comment);
      if(type!=this->unique_name()){
	cerr << "frehlich_inner_scale::read error - file of type " 
	     << type << " rather than type "
	     << this->unique_name() << endl;
	throw(string("frehlich_inner_scale::read"));
      } else 
	iof.read_key("INSCALE", inner_scale_value, comment);
    }
  }
  
  void frehlich_inner_scale::write(const char * filename) const {
    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "frehlich_inner_scale::write - "
	   << "error opening file " << filename << endl;
      throw(string("frehlich_inner_scale::write"));
    }

    try{this->write(iof);}
    catch(...){
      cerr << "frehlich_inner_scale::write - "
	   << "error writing "
	   << this->unique_name() << " to file " 
	   << filename << endl;
      throw(string("frehlich_inner_scale::write"));
    }
  }

  void frehlich_inner_scale::write(iofits & iof) const {
    fits_header_data<char> tmphdr;
    tmphdr.write(iof);
    string comment = "type";
    string type = this->unique_name();
    iof.write_key("TYPE", type, comment); 
    comment = "inner scale value (meters)";
    iof.write_key("INSCALE", inner_scale_value, comment); 
  }

  void frehlich_inner_scale::print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << this->unique_name()
       << "/" << "object type" << endl;
    os << prefix << "INSCALE    = " << setw(vlspc) << inner_scale_value
       << "/" << "inner scale value (meters)" << endl;
  }

  bool operator==(const frehlich_inner_scale & fis1, const frehlich_inner_scale & fis2) {
    if(fis1.inner_scale_value!=fis2.inner_scale_value) return false;
    return true;
  }

  bool operator!=(const frehlich_inner_scale & fis1, const frehlich_inner_scale & fis2) {
    return !operator==(fis1, fis2);
  }

  namespace factory_register {
    const fits_keyval_set & get_power_law_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "power law"));
      return *fkvs;
    }
    
    AO_sim_base * create_power_law(const iofits & iof) {
      return new power_law(iof);
    }
  }

  const bool power_law::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_power_law_keyval_set(), 
				      factory_register::create_power_law);

  power_law::power_law(){
    exponent = 0;
    coefficient = 0;
  }

  power_law::power_law(const power_law & pl){
    this->operator=(pl);
  }

  power_law::power_law(const iofits & iof){
    this->read(iof);
  }

  power_law::power_law(double expnt, double coeff){

    if(expnt > 0){
      cerr << "power_law::power_law error - exponent " << expnt
	   << " cannot be less than zero\n";
      throw(string("power_law::power_law"));
    }
    exponent = expnt;
    if(coeff < 0){
      cerr << "power_law::power_law error - coefficient " << coeff
	   << " cannot be less than zero\n";
      throw(string("power_law::power_law"));
    }
    coefficient = coeff;
  }

  power_law::power_law(double expnt, 
		       double r_0_meters,
		       double r_0_ref_wavelength_meters){

    if(expnt > 0){
      cerr << "power_law::power_law error - exponent " 
	   << expnt
	   << " cannot be less than zero\n";
      throw(string("power_law::power_law"));
    }
    exponent = expnt;
    if(r_0_meters <= 0){
      cerr << "power_law::power_law error - r_0 " 
	   << r_0_meters
	   << " cannot be less than zero\n";
      throw(string("power_law::power_law"));
    }
    if(r_0_ref_wavelength_meters <= 0){
      cerr << "power_law::power_law error - r_0 reference wavelength " 
	   << r_0_ref_wavelength_meters
	   << " cannot be less than zero\n";
      throw(string("power_law::power_law"));
    }
    double integrated_cn2_profile = 
      pow(r_0_meters, -5/3.0)*r_0_ref_wavelength_meters*r_0_ref_wavelength_meters/.423/4/M_PI/M_PI; 
    coefficient = 2*M_PI*.033*integrated_cn2_profile;
  }

  power_law & power_law::operator=(const power_law & pl){
    if(this==&pl)
      return(*this);
    exponent = pl.exponent;
    coefficient = pl.coefficient;
    return(*this);
  }

  void power_law::print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << this->unique_name()
       << "/" << "type" << endl;
    os << prefix << "EXPONENT   = " << setw(vlspc) << exponent
       << "/" << "power law exponent" << endl;
    os << prefix << "COEFF      = " << setw(vlspc) << coefficient
       << "/" << "power law coefficient" << endl;
  }

  void power_law::read(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "power_law::read - "
	   << "error opening file " << filename << endl;
      throw(string("power_law::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "power_law::read - "
	   << "error reading "
	   << this->unique_name() << " from file "
	   << filename << endl;
      throw(string("power_law::read"));
    }
  }

  void power_law::read(const iofits & iof){
    if(!iof.key_exists("TYPE")){
      cerr << "power_law::read error - "
	   << "no type specified\n";
      throw(string("power_law::read"));
    } else {
      string type, comment;
      iof.read_key("TYPE", type, comment);
      if(type!=this->unique_name()){
	cerr << "power_law::read error - file of type " 
	     << type << " rather than type " 
	     << this->unique_name() << endl;
	throw(string("power_law::read"));
      } else {
	iof.read_key("EXPONENT", exponent, comment); 
	iof.read_key("COEFF", coefficient, comment); 
      }
    }
  }

  void power_law::write(const char * filename) const {
    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "power_law::write - "
	   << "error opening file " << filename << endl;
      throw(string("power_law::write"));
    }

    try{this->write(iof);}
    catch(...){
      cerr << "power_law::write - "
	   << "error writing "
	   << this->unique_name() << " to file " 
	   << filename << endl;
      throw(string("power_law::write"));
    }
  }

  void power_law::write(iofits & iof) const {
    fits_header_data<char> tmphdr;
    tmphdr.write(iof);
    string comment = "type";
    string type = this->unique_name();
    iof.write_key("TYPE", type, comment); 
    comment = "power law exponent";
    iof.write_key("EXPONENT", exponent, comment); 
    comment = "power law coefficient";
    iof.write_key("COEFF", coefficient, comment); 
  }

  power_law * 
  power_law::power_law_factory(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...) {
      cerr << "power_law::power_law_factory - "
	   << "error opening file " << filename << endl;
      throw(string("power_law::power_law_factory"));
    }
    return(power_law::power_law_factory(iof));
  }

  power_law * power_law::power_law_factory(const iofits & iof){
    AO_sim_base * aosb = AO_sim_base::AO_sim_base_factory(iof);
    power_law * plaw = dynamic_cast<power_law *>(aosb);
    if(plaw==NULL)
      throw(string("power_law::power_law_factory"));
    return(plaw);
  }

  bool operator==(const power_law & plaw1, 
		  const power_law & plaw2) {
    if(plaw1.exponent!=plaw2.exponent) return(false);
    if(plaw1.coefficient!=plaw2.coefficient) return(false);
    return(true);
  }

  bool operator!=(const power_law & plaw1, 
		  const power_law & plaw2) {
    return !operator==(plaw1, plaw2);
  }

  namespace factory_register {
    const fits_keyval_set & get_von_karman_power_law_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "von karman power law"));
      return *fkvs;
    }
    
    AO_sim_base * create_von_karman_power_law(const iofits & iof) {
      return new von_karman_power_law(iof);
    }
  }

  const bool von_karman_power_law::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_von_karman_power_law_keyval_set(), 
				      factory_register::create_von_karman_power_law); 

  von_karman_power_law::von_karman_power_law() {
    outer_scale_value = 0;
  }

  von_karman_power_law::
  von_karman_power_law(const von_karman_power_law & vk_power_law){
    this->operator=(vk_power_law);
  }

  von_karman_power_law::von_karman_power_law(const iofits & iof){
    this->read(iof);
  }

  von_karman_power_law::
  von_karman_power_law(double expnt, 
		       double coeff, 
		       double outscale){
    if(expnt > 0){
      cerr << "von_karman_power_law::von_karman_power_law error - "
	   << "exponent " << expnt
	   << " cannot be less than zero\n";
      throw(string("von_karman_power_law::von_karman_power_law"));
    }
    exponent = expnt;
    if(coeff < 0){
      cerr << "von_karman_power_law::von_karman_power_law error - "
	   << "coefficient " << coeff
	   << " cannot be less than zero\n";
      throw(string("von_karman_power_law::von_karman_power_law"));
    }
    coefficient = coeff;
    if(outscale < 0){
      cerr << "von_karman_power_law::von_karman_power_law error - "
	   << "outer scale value " << outscale
	   << " cannot be less than zero\n";
      throw(string("von_karman_power_law::von_karman_power_law"));
    }
    outer_scale_value = outscale;
  }

  von_karman_power_law::von_karman_power_law(double expnt, 
					     double r_0_meters,
					     double r_0_ref_wavelength_meters,
					     double outscale){

    if(expnt > 0){
      cerr << "von_karman_power_law::von_karman_power_law error - exponent " 
	   << expnt
	   << " cannot be less than zero\n";
      throw(string("von_karman_power_law::von_karman_power_law"));
    }
    exponent = expnt;
    if(r_0_meters <= 0){
      cerr << "von_karman_power_law::von_karman_power_law error - r_0 " 
	   << r_0_meters
	   << " cannot be less than zero\n";
      throw(string("von_karman_power_law::von_karman_power_law"));
    }
    if(r_0_ref_wavelength_meters <= 0){
      cerr << "von_karman_power_law::von_karman_power_law error - r_0 reference wavelength " 
	   << r_0_ref_wavelength_meters
	   << " cannot be less than zero\n";
      throw(string("von_karman_power_law::von_karman_power_law"));
    }
    double integrated_cn2_profile = 
      pow(r_0_meters, -5/3.0)*r_0_ref_wavelength_meters*r_0_ref_wavelength_meters/.423/4/M_PI/M_PI; 
    coefficient = 2*M_PI*.033*integrated_cn2_profile;

    if(outscale < 0){
      cerr << "von_karman_power_law::von_karman_power_law error - "
	   << "outer scale value " << outscale
	   << " cannot be less than zero\n";
      throw(string("von_karman_power_law::von_karman_power_law"));
    }
    outer_scale_value = outscale;
  }

  von_karman_power_law & von_karman_power_law::
  operator=(const von_karman_power_law & vk_power_law){
    if(this==&vk_power_law)
      return(*this);
    this->power_law::operator=(vk_power_law);
    outer_scale_value = vk_power_law.outer_scale_value;
    return(*this);
  }

  void von_karman_power_law::read(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "von_karman_power_law::read - "
	   << "error opening file " << filename << endl;
      throw(string("von_karman_power_law::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "von_karman_power_law::read - "
	   << "error reading "
	   << this->unique_name() << " from file "
	   << filename << endl;
      throw(string("von_karman_power_law::read"));
    }
  }

  void von_karman_power_law::read(const iofits & iof){ 
    if(!iof.key_exists("TYPE")){
      cerr << "von_karman_power_law::read error - "
	   << "no power law type specified\n";
      throw(string("von_karman_power_law::read"));
    } else {
      string type, comment;
      iof.read_key("TYPE", type, comment);
      if(type!=this->unique_name()){
	cerr << "von_karman_power_law::read error - file of type " 
	     << type << " rather than type "
	     << this->unique_name() << endl;
	throw(string("von_karman_power_law::read"));
      } else {
	iof.read_key("EXPONENT", exponent, comment); 
	iof.read_key("COEFF", coefficient, comment); 
	iof.read_key("OSCALE", outer_scale_value, comment);
      }
    }
  }

  void von_karman_power_law::write(const char * filename) const {
    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "von_karman_power_law::write - "
	   << "error opening file " << filename << endl;
      throw(string("von_karman_power_law::write"));
    }

    try{this->write(iof);}
    catch(...){
      cerr << "von_karman_power_law::write - "
	   << "error writing "
	   << this->unique_name() << " to file " 
	   << filename << endl;
      throw(string("von_karman_power_law::write"));
    }
  }

  void von_karman_power_law::write(iofits & iof) const {
    fits_header_data<char> tmphdr;
    tmphdr.write(iof);
    string type = this->unique_name();
    string comment = "type";
    iof.write_key("TYPE", type, comment);
    comment = "power law exponent";
    iof.write_key("EXPONENT", exponent, comment); 
    comment = "power law coefficient";
    iof.write_key("COEFF", coefficient, comment); 
    comment = "outer scale value (meters)";
    iof.write_key("OSCALE", outer_scale_value, comment);
  }

  void von_karman_power_law::print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << this->unique_name()
       << "/" << "type" << endl;
    os << prefix << "OSCALE     = " << setw(vlspc) << outer_scale_value
       << "/" << "outer scale value (meters)" << endl;
  }

  bool operator==(const von_karman_power_law & vkplaw1, 
		  const von_karman_power_law & vkplaw2) {
    if(vkplaw1.outer_scale_value!=vkplaw2.outer_scale_value) return(false);
    return(operator==(static_cast<power_law>(vkplaw1),static_cast<power_law>(vkplaw2)));
  }

  bool operator!=(const von_karman_power_law & vkplaw1, 
		  const von_karman_power_law & vkplaw2) {
    return !operator==(vkplaw1, vkplaw2);
  }

  namespace factory_register {
    const fits_keyval_set & get_greenwood_power_law_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "greenwood power law"));
      return *fkvs;
    }
    
    AO_sim_base * create_greenwood_power_law(const iofits & iof) {
      return new greenwood_power_law(iof);
    }
  }

  const bool greenwood_power_law::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_greenwood_power_law_keyval_set(), 
				      factory_register::create_greenwood_power_law);

  greenwood_power_law::greenwood_power_law() {
    outer_scale_value = 0;
  }

  greenwood_power_law::
  greenwood_power_law(const greenwood_power_law & vk_power_law){
    this->operator=(vk_power_law);
  }

  greenwood_power_law::greenwood_power_law(const iofits & iof){
    this->read(iof);
  }

  greenwood_power_law::
  greenwood_power_law(double expnt, 
		      double coeff, 
		      double outscale){
    if(expnt > 0){
      cerr << "greenwood_power_law::greenwood_power_law error - "
	   << "exponent " << expnt
	   << " cannot be less than zero\n";
      throw(string("greenwood_power_law::greenwood_power_law"));
    }
    exponent = expnt;
    if(coeff < 0){
      cerr << "greenwood_power_law::greenwood_power_law error - "
	   << "coefficient " << coeff
	   << " cannot be less than zero\n";
      throw(string("greenwood_power_law::greenwood_power_law"));
    }
    coefficient = coeff;
    if(outscale < 0){
      cerr << "greenwood_power_law::greenwood_power_law error - "
	   << "outer scale value " << outscale
	   << " cannot be less than zero\n";
      throw(string("greenwood_power_law::greenwood_power_law"));
    }
    outer_scale_value = outscale;
  }

  greenwood_power_law::greenwood_power_law(double expnt, 
					   double r_0_meters,
					   double r_0_ref_wavelength_meters,
					   double outscale){

    if(expnt > 0){
      cerr << "greenwood_power_law::greenwood_power_law error - exponent " 
	   << expnt
	   << " cannot be less than zero\n";
      throw(string("greenwood_power_law::greenwood_power_law"));
    }
    exponent = expnt;
    if(r_0_meters <= 0){
      cerr << "greenwood_power_law::greenwood_power_law error - r_0 " 
	   << r_0_meters
	   << " cannot be less than zero\n";
      throw(string("greenwood_power_law::greenwood_power_law"));
    }
    if(r_0_ref_wavelength_meters <= 0){
      cerr << "greenwood_power_law::greenwood_power_law error - r_0 reference wavelength " 
	   << r_0_ref_wavelength_meters
	   << " cannot be less than zero\n";
      throw(string("greenwood_power_law::greenwood_power_law"));
    }
    double integrated_cn2_profile = 
      pow(r_0_meters, -5/3.0)*r_0_ref_wavelength_meters*r_0_ref_wavelength_meters/.423/4/M_PI/M_PI; 
    coefficient = 2*M_PI*.033*integrated_cn2_profile;

    if(outscale < 0){
      cerr << "greenwood_power_law::greenwood_power_law error - "
	   << "outer scale value " << outscale
	   << " cannot be less than zero\n";
      throw(string("greenwood_power_law::greenwood_power_law"));
    }
    outer_scale_value = outscale;
  }

  greenwood_power_law & greenwood_power_law::
  operator=(const greenwood_power_law & greenwood_power_law){
    if(this==&greenwood_power_law)
      return(*this);
    this->power_law::operator=(greenwood_power_law);
    outer_scale_value = greenwood_power_law.outer_scale_value;
    return(*this);
  }

  void greenwood_power_law::read(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "greenwood_power_law::read - "
	   << "error opening file " << filename << endl;
      throw(string("greenwood_power_law::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "greenwood_power_law::read - "
	   << "error reading "
	   << this->unique_name() << " from file "
	   << filename << endl;
      throw(string("greenwood_power_law::read"));
    }
  }

  void greenwood_power_law::read(const iofits & iof){ 
    if(!iof.key_exists("TYPE")){
      cerr << "greenwood_power_law::read error - "
	   << "no power law type specified\n";
      throw(string("greenwood_power_law::read"));
    } else {
      string type, comment;
      iof.read_key("TYPE", type, comment);
      if(type!=this->unique_name()){
	cerr << "greenwood_power_law::read error - outer scale of type " 
	     << type << " rather than type "
	     << this->unique_name() << endl;
	throw(string("greenwood_power_law::read"));
      } else {
	iof.read_key("EXPONENT", exponent, comment); 
	iof.read_key("COEFF", coefficient, comment); 
	iof.read_key("OSCALE", outer_scale_value, comment);
      }
    }
  }

  void greenwood_power_law::write(const char * filename) const {
    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "greenwood_power_law::write - "
	   << "error opening file " << filename << endl;
      throw(string("greenwood_power_law::write"));
    }

    try{this->write(iof);}
    catch(...){
      cerr << "greenwood_power_law::write - "
	   << "error writing "
	   << this->unique_name() << " to file " 
	   << filename << endl;
      throw(string("greenwood_power_law::write"));
    }
  }

  void greenwood_power_law::write(iofits & iof) const {
    fits_header_data<char> tmphdr;
    tmphdr.write(iof);
    string type = this->unique_name();
    string comment = "type";
    iof.write_key("TYPE", type, comment);
    comment = "power law exponent";
    iof.write_key("EXPONENT", exponent, comment); 
    comment = "power law coefficient";
    iof.write_key("COEFF", coefficient, comment); 
    comment = "outer scale value (meters)";
    iof.write_key("OSCALE", outer_scale_value, comment);
  }

  void greenwood_power_law::print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << this->unique_name()
       << "/" << "outer scale type" << endl;
    os << prefix << "OSCALE     = " << setw(vlspc) << outer_scale_value
       << "/" << "outer scale value (meters)" << endl;
  }

  bool operator==(const greenwood_power_law & gwplaw1, 
		  const greenwood_power_law & gwplaw2) {
    if(gwplaw1.outer_scale_value!=gwplaw2.outer_scale_value) return(false);
    return(operator==(static_cast<power_law>(gwplaw1),static_cast<power_law>(gwplaw2)));
  }

  bool operator!=(const greenwood_power_law & gwplaw1, 
		  const greenwood_power_law & gwplaw2) {
    return !operator==(gwplaw1, gwplaw2);
  }

  int power_spectrum::verbose_level = 0;
  
  power_spectrum * power_spectrum::
  power_spectrum_factory(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "power_spectrum::power_spectrum_factory - error opening file " 
	   << filename << endl;
      throw(string("power_spectrum::power_spectrum_factory"));
    }
    return(power_spectrum::power_spectrum_factory(iof));
  }

  power_spectrum * power_spectrum::power_spectrum_factory(const iofits & iof){
    AO_sim_base * aosb = AO_sim_base::AO_sim_base_factory(iof);
    power_spectrum * pspec = dynamic_cast<power_spectrum *>(aosb);
    if(pspec==NULL)
      throw(string("power_spectrum::power_spectrum_factory"));
    return(pspec);
  } 

  namespace factory_register {

    template<class power_law_type, class inner_scale_type> 
    AO_sim_base* create_isotropic_power_law_spectrum(const iofits& iof){
      return new isotropic_power_law_spectrum<power_law_type, inner_scale_type>(iof);
    }

    const fits_keyval_set & get_isotropic_power_law_spectrum_power_law_null_inscle_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE","isotropic power law spectrum"));
      fkvs->push_back(fits_keyval_entry("TYPE","power law",1));
      fkvs->push_back(fits_keyval_entry("TYPE","null inner scale",2));
      return *fkvs;
    }
    
    const fits_keyval_set & get_isotropic_power_law_spectrum_power_law_exponential_inscle_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE","isotropic power law spectrum"));
      fkvs->push_back(fits_keyval_entry("TYPE","power law",1));
      fkvs->push_back(fits_keyval_entry("TYPE","exponential inner scale",2));
      return *fkvs;
    }
    
    const fits_keyval_set & get_isotropic_power_law_spectrum_power_law_frehlich_inscle_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE","isotropic power law spectrum"));
      fkvs->push_back(fits_keyval_entry("TYPE","power law",1));
      fkvs->push_back(fits_keyval_entry("TYPE","frehlich inner scale",2));
      return *fkvs;
    }
    

    const fits_keyval_set & get_isotropic_power_law_spectrum_von_karman_power_law_null_inscle_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE","isotropic power law spectrum"));
      fkvs->push_back(fits_keyval_entry("TYPE","von karman power law",1));
      fkvs->push_back(fits_keyval_entry("TYPE","null inner scale",2));
      return *fkvs;
    }
    
    const fits_keyval_set & get_isotropic_power_law_spectrum_von_karman_power_law_exponential_inscle_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE","isotropic power law spectrum"));
      fkvs->push_back(fits_keyval_entry("TYPE","von karman power law",1));
      fkvs->push_back(fits_keyval_entry("TYPE","exponential inner scale",2));
      return *fkvs;
    }
    
    const fits_keyval_set & get_isotropic_power_law_spectrum_von_karman_power_law_frehlich_inscle_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE","isotropic power law spectrum"));
      fkvs->push_back(fits_keyval_entry("TYPE","von karman power law",1));
      fkvs->push_back(fits_keyval_entry("TYPE","frehlich inner scale",2));
      return *fkvs;
    }
    

    const fits_keyval_set & get_isotropic_power_law_spectrum_greenwood_power_law_null_inscle_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE","isotropic power law spectrum"));
      fkvs->push_back(fits_keyval_entry("TYPE","greenwood power law",1));
      fkvs->push_back(fits_keyval_entry("TYPE","null inner scale",2));
      return *fkvs;
    }
    
    const fits_keyval_set & get_isotropic_power_law_spectrum_greenwood_power_law_exponential_inscle_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE","isotropic power law spectrum"));
      fkvs->push_back(fits_keyval_entry("TYPE","greenwood power law",1));
      fkvs->push_back(fits_keyval_entry("TYPE","exponential inner scale",2));
      return *fkvs;
    }
    
    const fits_keyval_set & get_isotropic_power_law_spectrum_greenwood_power_law_frehlich_inscle_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE","isotropic power law spectrum"));
      fkvs->push_back(fits_keyval_entry("TYPE","greenwood power law",1));
      fkvs->push_back(fits_keyval_entry("TYPE","frehlich inner scale",2));
      return *fkvs;
    }

  }
    
  template<>
  const bool isotropic_power_law_spectrum<power_law, null_inner_scale>::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_isotropic_power_law_spectrum_power_law_null_inscle_keyval_set(), 
				      factory_register::create_isotropic_power_law_spectrum<power_law, null_inner_scale>);
  
  template<>
  const bool isotropic_power_law_spectrum<power_law, exponential_inner_scale>::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_isotropic_power_law_spectrum_power_law_exponential_inscle_keyval_set(), 
				      factory_register::create_isotropic_power_law_spectrum<power_law, exponential_inner_scale>);
  
  template<>
  const bool isotropic_power_law_spectrum<power_law, frehlich_inner_scale>::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_isotropic_power_law_spectrum_power_law_frehlich_inscle_keyval_set(), 
				      factory_register::create_isotropic_power_law_spectrum<power_law, frehlich_inner_scale>);
  

  template<>
  const bool isotropic_power_law_spectrum<von_karman_power_law, null_inner_scale>::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_isotropic_power_law_spectrum_von_karman_power_law_null_inscle_keyval_set(), 
				      factory_register::create_isotropic_power_law_spectrum<von_karman_power_law, null_inner_scale>);
  
  template<>
  const bool isotropic_power_law_spectrum<von_karman_power_law, exponential_inner_scale>::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_isotropic_power_law_spectrum_von_karman_power_law_exponential_inscle_keyval_set(), 
				      factory_register::create_isotropic_power_law_spectrum<von_karman_power_law, exponential_inner_scale>);
  
  template<>
  const bool isotropic_power_law_spectrum<von_karman_power_law, frehlich_inner_scale>::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_isotropic_power_law_spectrum_von_karman_power_law_frehlich_inscle_keyval_set(), 
				      factory_register::create_isotropic_power_law_spectrum<von_karman_power_law, frehlich_inner_scale>);


  template<>
  const bool isotropic_power_law_spectrum<greenwood_power_law, null_inner_scale>::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_isotropic_power_law_spectrum_greenwood_power_law_null_inscle_keyval_set(), 
				      factory_register::create_isotropic_power_law_spectrum<greenwood_power_law, null_inner_scale>);

  template<>
  const bool isotropic_power_law_spectrum<greenwood_power_law, exponential_inner_scale>::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_isotropic_power_law_spectrum_greenwood_power_law_exponential_inscle_keyval_set(), 
				      factory_register::create_isotropic_power_law_spectrum<greenwood_power_law, exponential_inner_scale>);

  template<>
  const bool isotropic_power_law_spectrum<greenwood_power_law, frehlich_inner_scale>::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_isotropic_power_law_spectrum_greenwood_power_law_frehlich_inscle_keyval_set(), 
				      factory_register::create_isotropic_power_law_spectrum<greenwood_power_law, frehlich_inner_scale>);

}

