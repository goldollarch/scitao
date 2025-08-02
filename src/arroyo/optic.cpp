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
#include "iofits.h"
#include "computational_geometry.h"
#include "optic.h"
//#include "geometric_wavefront.h"
#include "diffractive_wavefront.h"

using namespace std;

namespace Arroyo {

  int optic::verbose_level = 0;

  optic * optic::optic_factory(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...) {
      cerr << "optic::optic_factory - "
	   << "error opening file " << filename << endl;
      throw(string("optic::optic_factory"));
    }
    return(optic::optic_factory(iof));
  }

  optic * optic::optic_factory(const iofits & iof){
    AO_sim_base * aosb = AO_sim_base::AO_sim_base_factory(iof);
    optic * o = dynamic_cast<optic *>(aosb);
    if(o==NULL)
      throw(string("optic::optic_factory"));
     return(o);    
  }

  optic::optic(const optic & op) {
    this->operator=(op);
  }

  optic & optic::operator=(const optic & op) {
    if(this==&op)
      return(*this);
    foreshortening = op.foreshortening;
    return(*this);
  }

  void optic::read(const iofits & iof){ 
    string comment;
    iof.read_key("FRSHRTEN", foreshortening, comment);
  }
 
  void optic::write(iofits & iof) const {
    iof.write_key("FRSHRTEN", foreshortening, "foreshorten optic");
  }

  void optic::print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "FRSHRTEN   = " << setw(vlspc) << foreshortening
       << "/" << "foreshorten optic" << endl;
  }

  plane_optic * plane_optic::plane_optic_factory(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...) {
      cerr << "plane_optic::plane_optic_factory - "
	   << "error opening file " << filename << endl;
      throw(string("plane_optic::plane_optic_factory"));
    }
    return(plane_optic::plane_optic_factory(iof));
  }

  plane_optic * plane_optic::plane_optic_factory(const iofits & iof){

    optic * o = optic::optic_factory(iof);
    plane_optic * po = dynamic_cast<plane_optic *>(o);
    if(po==NULL)
      throw(string("plane_optic::plane_optic_factory"));
     return(po);
  }

  plane_optic::plane_optic(const plane_optic & plane_op) {
    this->operator=(plane_op);
  }

  plane_optic & plane_optic::operator=(const plane_optic & plane_op) {
    if(this==&plane_op)
      return(*this);
    this->three_frame::operator=(plane_op);
    this->optic::operator=(plane_op);
    return(*this);
  }

  void plane_optic::read(const iofits & iof){ 
    this->three_frame::read(iof);
    this->optic::read(iof);
  }
 
  void plane_optic::write(iofits & iof) const {
    this->three_frame::write(iof);
    this->optic::write(iof);
  }

  void plane_optic::print(ostream & os, const char * prefix) const {
    this->optic::print(os, prefix);
    this->three_frame::print(os, prefix);
  }

  three_point plane_optic::get_point_of_intersection(const three_point & tp, 
						     const three_vector & tv) const {

    three_point intersection_point;
    try{intersection_point = Arroyo::get_ray_plane_intersection(tp, tv, *this, this->z());}
    catch(...){
      cerr << "plane_optic::get_point_of_intersection error - "
	   << "could not get point of intersection\n";
      throw(string("plane_optic::get_point_of_intersection"));
    }
    return(intersection_point);
  } 
 
  //geometric_ray * one_to_one_optic::get_wavefront_data(const geometric_wavefront & gwf) const {
  //  return(gwf.rays);
  //}
}
