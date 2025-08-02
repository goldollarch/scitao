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
#include "fits_factory.h"
#include "emitter.h"

using namespace std;

namespace Arroyo {

  emitter * emitter::emitter_factory(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...) {
      cerr << "emitter::emitter_factory - "
	   << "error opening file " << filename << endl;
      throw(string("emitter::emitter_factory"));
    }
    return(emitter::emitter_factory(iof));
  }

  emitter * emitter::emitter_factory(const iofits & iof){
    AO_sim_base * aosb = AO_sim_base::AO_sim_base_factory(iof);
    emitter * e = dynamic_cast<emitter *>(aosb);
    if(e==NULL){
      cerr << "emitter::emitter_factory error - "
	   << "could not read emitter from iofits object\n";
      throw(string("emitter::emitter_factory"));
    }
    return(e);    
  }

  // This is a poor implementation, but will have to do for the moment
  emitter * emitter::emitter_factory(const emitter * emtr){

    if(emtr==NULL) return(NULL);

    const plane_wave_emitter * pwe = dynamic_cast<const plane_wave_emitter *>(emtr);
    const spherical_wave_emitter * swe = dynamic_cast<const spherical_wave_emitter *>(emtr);
    if(pwe!=NULL) return(new plane_wave_emitter(*pwe));
    else if(swe!=NULL) return(new spherical_wave_emitter(*swe));
    else{
      cerr << "emitter::emitter_factory error - unrecognized emitter\n";
      emtr->print(cerr, "emtr ");
      throw(string("emitter::emitter_factory"));
    }
  }

  void emitter::write(iofits & iof) const {
    fits_header_data<char> tmphdr;
    tmphdr.write(iof);
  }

  namespace factory_register {
    const fits_keyval_set & get_plane_wave_emitter_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "plane wave emitter"));
      return *fkvs;
    }
    
    AO_sim_base * create_plane_wave_emitter(const iofits & iof) {
      return new plane_wave_emitter(iof);
    }
  }

  const bool plane_wave_emitter::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_plane_wave_emitter_keyval_set(), 
				      factory_register::create_plane_wave_emitter);

  plane_wave_emitter::plane_wave_emitter(const three_vector & tv) {
    if(tv.length()<three_frame::precision){
      cerr << "plane_wave_emitter::plane_wave_emitter error - "
	   << "null three_vector supplied to constructor\n";
    }
    this->three_vector::operator=(tv*(1/tv.length()));
  }

  plane_wave_emitter::plane_wave_emitter(const char * filename) {
    this->read(filename);
  }

  plane_wave_emitter::plane_wave_emitter(const plane_wave_emitter & pwe) {
    this->operator=(pwe);
  }

  plane_wave_emitter & plane_wave_emitter::operator=(const plane_wave_emitter & pwe) {
    if(this==&pwe)
      return(*this);
    this->three_vector::operator=(pwe);
    return(*this);
  }

  void plane_wave_emitter::read(const char * filename) {
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "plane_wave_emitter::read - "
	   << "error opening file " << filename << endl;
      throw(string("plane_wave_emitter::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "plane_wave_emitter::read - "
	   << "error reading "
	   << this->unique_name() << " from file "
	   << filename << endl;
      throw(string("plane_wave_emitter::read"));
    }
  }

  void plane_wave_emitter::read(const iofits & iof) {
    if(!iof.key_exists("TYPE")){
      cerr << "plane_wave_emitter::read error - "
	   << "unrecognized type of file\n";
      throw(string("plane_wave_emitter::read"));
    }
    string type, comment;
    iof.read_key("TYPE", type, comment);
    if(type!=this->unique_name()){
      cerr << "plane_wave_emitter::read error - file of type " 
	   << type << " rather than of type "
	   << this->unique_name() << endl;
      throw(string("plane_wave_emitter::read"));
    }
    this->three_vector::read(iof);
    // leave iof pointing to the next header, if it exists
    if(iof.get_hdu_num()!=iof.get_num_hdus())
      iof.movrel_hdu(1);
  }

  void plane_wave_emitter::write(const char * filename) const {
   iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "plane_wave_emitter::write - "
	   << "error opening file " << filename << endl;
      throw(string("plane_wave_emitter::write"));
    }

    try{this->write(iof);}
    catch(...){
      cerr << "plane_wave_emitter::write - "
	   << "error writing "
	   << this->unique_name() << " to file " 
	   << filename << endl;
      throw(string("plane_wave_emitter::write"));
    }
  }

  void plane_wave_emitter::write(iofits & iof) const {
    this->emitter::write(iof);
    string type = this->unique_name();
    string comment = "object type";
    iof.write_key("TYPE", type, comment);
    this->three_vector::write(iof);    
  }

  void plane_wave_emitter::print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << this->unique_name()
       << "/" << "object type" << endl;
    this->three_vector::print(os, prefix);
  }

  bool operator==(const plane_wave_emitter & pwe1, const plane_wave_emitter & pwe2) {
    if(operator==(static_cast<const three_vector>(pwe1),
		  static_cast<const three_vector>(pwe2))) return true;
    return false;
  }

  bool operator!=(const plane_wave_emitter & pwe1, const plane_wave_emitter & pwe2) {
    return !operator==(pwe1, pwe2);
  }

  namespace factory_register {
    const fits_keyval_set & get_spherical_wave_emitter_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "spherical wave emitter"));
      return *fkvs;
    }
    
    AO_sim_base * create_spherical_wave_emitter(const iofits & iof) {
      return new spherical_wave_emitter(iof);
    }
  }

  const bool spherical_wave_emitter::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_spherical_wave_emitter_keyval_set(), 
				      factory_register::create_spherical_wave_emitter);

  spherical_wave_emitter::spherical_wave_emitter(const three_point & tp) {
    this->three_point::operator=(tp);
  }

  spherical_wave_emitter::spherical_wave_emitter(const char * filename) {
    this->read(filename);
  }

  spherical_wave_emitter::spherical_wave_emitter(const spherical_wave_emitter & swe) {
    this->operator=(swe);
  }

  spherical_wave_emitter & spherical_wave_emitter::operator=(const spherical_wave_emitter & swe) {
    if(this==&swe)
      return(*this);
    this->three_point::operator=(swe);
    return(*this);
  }

  void spherical_wave_emitter::read(const char * filename) {
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "spherical_wave_emitter::read - "
	   << "error opening file " << filename << endl;
      throw(string("spherical_wave_emitter::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "spherical_wave_emitter::read - "
	   << "error reading "
	   << this->unique_name() << " from file "
	   << filename << endl;
      throw(string("spherical_wave_emitter::read"));
    }
  }

  void spherical_wave_emitter::read(const iofits & iof) {
    if(!iof.key_exists("TYPE")){
      cerr << "spherical_wave_emitter::read error - "
	   << "unrecognized type of file\n";
      throw(string("spherical_wave_emitter::read"));
    }
    string type, comment;
    iof.read_key("TYPE", type, comment);
    if(type!=this->unique_name()){
      cerr << "spherical_wave_emitter::read error - file of type " 
	   << type << " rather than of type "
	   << this->unique_name() << endl;
      throw(string("spherical_wave_emitter::read"));
    }
    this->three_point::read(iof);
    // leave iof pointing to the next header, if it exists
    if(iof.get_hdu_num()!=iof.get_num_hdus())
      iof.movrel_hdu(1);
  }

  void spherical_wave_emitter::write(const char * filename) const {
   iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "spherical_wave_emitter::write - "
	   << "error opening file " << filename << endl;
      throw(string("spherical_wave_emitter::write"));
    }

    try{this->write(iof);}
    catch(...){
      cerr << "spherical_wave_emitter::write - "
	   << "error writing "
	   << this->unique_name() << " to file " 
	   << filename << endl;
      throw(string("spherical_wave_emitter::write"));
    }
  }

  void spherical_wave_emitter::write(iofits & iof) const {
    this->emitter::write(iof);
    string type = this->unique_name();
    string comment = "object type";
    iof.write_key("TYPE", type, comment);
    this->three_point::write(iof);    
  }

  void spherical_wave_emitter::print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << this->unique_name()
       << "/" << "object type" << endl;
    this->three_point::print(os, prefix);
  }

  bool operator==(const spherical_wave_emitter & swe1, const spherical_wave_emitter & swe2) {
    if(operator==(static_cast<const three_point>(swe1), static_cast<const three_point>(swe2))) return true;
    return false;
  }

  bool operator!=(const spherical_wave_emitter & swe1, const spherical_wave_emitter & swe2) {
    return !operator==(swe1, swe2);
  }
}
