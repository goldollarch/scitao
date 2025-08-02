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
#include <vector>
#include "fits_factory.h"
#include "region_base.h"
#include "diffractive_wavefront.h"
#include "Shack_Hartmann_centroids.h"

using namespace std;

namespace Arroyo {

  namespace factory_register {
     const fits_keyval_set & get_Shack_Hartmann_centroids_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "Shack Hartmann centroids"));
      return *fkvs;
    }
    
    AO_sim_base * create_Shack_Hartmann_centroids(const iofits & iof) {
      return new Shack_Hartmann_centroids(iof);
    }
  }

  const bool Shack_Hartmann_centroids::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_Shack_Hartmann_centroids_keyval_set(), 
				      factory_register::create_Shack_Hartmann_centroids);

  Shack_Hartmann_centroids & Shack_Hartmann_centroids::operator=(const Shack_Hartmann_centroids & shcentroids) {
    if(this==&shcentroids) return(*this);
    this->pixel_array<double>::operator=(shcentroids);
    return(*this);
  }

  int Shack_Hartmann_centroids::verbose_level = 0;

  Shack_Hartmann_centroids::Shack_Hartmann_centroids(const vector<long> & lenslet_axes) {    
    if(lenslet_axes.size()!=2 || lenslet_axes[0]<=0 || lenslet_axes[1]<=0){
      cerr << "Shack_Hartmann_centroids::Shack_Hartmann_centroids error - "
	   << "invalid lenslet dimensions\n";
      throw(string("Shack_Hartmann_centroids::Shack_Hartmann_centroids"));      
    }

    vector<long> shaxes = lenslet_axes;
    shaxes[1] *= 2;
    this->pixel_array<double>::set_axes(shaxes);
  }

  void Shack_Hartmann_centroids::read(const char * filename) {
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "Shack_Hartmann_centroids::read - "
	   << "error opening file " << filename << endl;
      throw(string("Shack_Hartmann_centroids::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "Shack_Hartmann_centroids::read - "
	   << "error reading "
	   << this->unique_name() << " from file " 
	   << filename << endl;
      throw(string("Shack_Hartmann_centroids::read"));
    }
  }

  void Shack_Hartmann_centroids::read(const iofits & iof) {
    if(!iof.key_exists("TYPE")){
      cerr << "Shack_Hartmann_centroids::read error - "
	   << "no type identifier specified\n";
      throw(string("Shack_Hartmann_centroids::read"));
    } else {
      string type, comment;
      iof.read_key("TYPE", type, comment);
      if(type!=this->unique_name()){
	cerr << "Shack_Hartmann_centroids::read error - file of type " 
	     << type << " rather than type "
	     << this->unique_name() << endl;
	throw(string("Shack_Hartmann_centroids::read"));
      } else {
	this->pixel_array<double>::read(iof);
      }
    }

  }
 
  void Shack_Hartmann_centroids::write(const char * filename) const {
    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "Shack_Hartmann_centroids::write - "
	   << "error opening file " << filename << endl;
      throw(string("Shack_Hartmann_centroids::write"));
    }

    try{this->write(iof);}
    catch(...){
      cerr << "Shack_Hartmann_centroids::write - "
	   << "error writing "
	   << this->unique_name() << " to file "
	   << filename << endl;
      throw(string("Shack_Hartmann_centroids::write"));
    }
  }

  void Shack_Hartmann_centroids::write(iofits & iof) const {
    string type = this->unique_name();
    string comment = "object type";
    fits_header_data<double> fhd(this->get_axes());
    fhd.write(iof);
    iof.write_key("TYPE", type, comment);
    this->pixel_array<double>::write(iof);
  }

  void Shack_Hartmann_centroids::print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 

    fits_header_data<double> fhd(this->get_axes());
    os << prefix << "TYPE       = " << setw(vlspc) << this->unique_name()
       << "/" << "object type" << endl;
    fhd.print(os, prefix);
  }

}
