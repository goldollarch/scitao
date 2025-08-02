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
#include <iomanip>
#include "fits_factory.h"
#include "geometric_wavefront.h"
#include "sim_utils.h"

using namespace std;

namespace Arroyo {

  namespace factory_register {
    const fits_keyval_set & get_geometric_wavefront_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "geometric wavefront"));
      return *fkvs;
    }
    
    AO_sim_base * create_geometric_wavefront(const iofits & iof) {
      return new geometric_wavefront(iof);
    }
  } 

  const bool geometric_wavefront::factory_registration = 
    fits_factory<AO_sim_base>::Register(factory_register::get_geometric_wavefront_keyval_set(),
					factory_register::create_geometric_wavefront);

  geometric_wavefront::geometric_wavefront() {
    rays = NULL;
  }

  geometric_wavefront::geometric_wavefront(const geometric_wavefront & gwf) {
    rays = NULL;
    this->operator=(gwf);
  }

  geometric_wavefront::geometric_wavefront(const char * filename) {
    rays = NULL;
    this->read(filename);
  }

  geometric_wavefront::geometric_wavefront(const iofits & iof) {
    rays = NULL;
    this->read(iof);
  }

  geometric_wavefront::geometric_wavefront(const geometric_wavefront_header & gwfh) {
    this->geometric_wavefront_header::operator=(gwfh);
    if(nrays>0){
      try{rays = new geometric_ray[nrays];}
      catch(...){
	cerr << "geometric_wavefront::geometric_wavefront error - "
	     << "could not allocate memory for geometric rays\n";
	throw;
      }
    } else rays = NULL;
  }

  geometric_wavefront::~geometric_wavefront() {
    if(rays!=NULL) delete [] rays;
  }

  geometric_wavefront & geometric_wavefront::operator=(const geometric_wavefront & gwf) {
    if(this==&gwf)
      return(*this);
    if(nrays!=gwf.nrays && rays!=NULL){
      delete [] rays;
      nrays = 0;
    }
    if(gwf.nrays!=0) rays = new geometric_ray[gwf.nrays];
    nrays = gwf.nrays;
    return(*this);
  }

  void geometric_wavefront::read(const char * filename) {
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "geometric_wavefront::read - "
	   << "error opening file " << filename << endl;
      throw(string("geometric_wavefront::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "geometric_wavefront::read - "
	   << "error reading "
	   << this->unique_name() << " from file "
	   << filename << endl;
      throw(string("geometric_wavefront::read"));
    }
  }

  void geometric_wavefront::read(const iofits & iof) {

    string comment;
    long tmp_nrays;
    iof.read_key("NRAYS", tmp_nrays, comment);

    if(tmp_nrays==0 || nrays!=tmp_nrays){
      if(rays!=NULL) delete [] rays;
      rays = NULL;
      nrays = 0;
    }

    this->geometric_wavefront_header::read(iof);
    if(nrays == 0) return;

    try{
      rays = new geometric_ray[nrays];
    } catch(...) {
      cerr << "geometric_wavefront::read error - "
	   << "could not allocate memory\n";
      throw;
    }

    // The strategy here is that we're going to read a single,
    // contiguous array of doubles from the iofits object into
    // the array rawdata, then parse out the elements of this 
    // array into the objects that compose an individual ray.
    //
    // There are 8 numbers for each ray, which are stored in 
    // the following order:  
    // 3 coordinates of its point
    // 3 components of its normal 
    // one pathlength
    // one number to determine if the ray has been lost
    //
    // The array rawdata holds the concatenated numbers for each 
    // of the nray geometric rays
    double * rawdata;
    try{rawdata = new double[nrays*8];} 
    catch(...){
      cerr << "geometric_wavefront::read error - "
	   << "could not allocate memory\n";
      throw;
    }

    try{
      if(wavefront::verbose_level) 
	cout << "geometric_wavefront::read - reading data\n";
      iof.read_image(0, nrays*8-1, rawdata);
    } catch(...){
      cerr << "geometric_wavefront::read - error reading image\n";
      throw(string("geometric_wavefront::read"));
    }
  
    // stored values are components relative to the default frame
    three_frame tf;
    for(int i=0; i<nrays; i++){
      rays[i].three_point::operator=(three_point(rawdata[i*8],rawdata[i*8+1],rawdata[i*8+2],tf));
      rays[i].three_vector::operator=(three_vector(rawdata[i*8+3],rawdata[i*8+4],rawdata[i*8+5],tf));
      rays[i].pathlength = rawdata[i*8+6];
      rays[i].lost = rawdata[i*8+7];
    }
    delete [] rawdata;
  }

  void geometric_wavefront::write(const char * filename) const {
    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "geometric_wavefront::write - "
	   << "error opening file " << filename << endl;
      throw(string("geometric_wavefront::write"));
    }
    try{this->write(iof);}
    catch(...){
      cerr << "geometric_wavefront::write - "
	   << "error writing "
	   << this->unique_name() << " to file "
	   << filename << endl;
      iof.print_header(cerr, "error");
      throw(string("geometric_wavefront::write"));
    }
  }

  void geometric_wavefront::write(iofits & iof) const {

    this->geometric_wavefront_header::write(iof);

    // The strategy here is that we're going to write a single,
    // contiguous array of doubles from the iofits object into
    // the array rawdata, then parse out the elements of this 
    // array into the objects that compose an individual ray.
    //
    // There are 8 numbers for each ray, which are stored in 
    // the following order:  
    // 3 coordinates of its point
    // 3 components of its normal 
    // one pathlength
    //
    // The 8 numbers for each ray are concatenated into the 
    // array rawdata
    double * rawdata;
    try{rawdata = new double[nrays*8];} 
    catch(...){
      cerr << "geometric_wavefront::read error - "
	   << "could not allocate memory\n";
      throw;
    }

    // stored values are components relative to the default frame
    three_frame tf;
    for(int i=0; i<nrays; i++){
      rawdata[i*8] = rays[i].three_point::x(tf);
      rawdata[i*8+1] = rays[i].three_point::y(tf);
      rawdata[i*8+2] = rays[i].three_point::z(tf);
      rawdata[i*8+3] = rays[i].three_vector::x(tf);
      rawdata[i*8+4] = rays[i].three_vector::y(tf);
      rawdata[i*8+5] = rays[i].three_vector::z(tf);
      rawdata[i*8+6] = rays[i].get_pathlength();
      rawdata[i*8+7] = rays[i].is_lost();
    }

    if(nrays==0) iof.write_image(0,0,rawdata);
    else iof.write_image(0,nrays*8-1,rawdata);

    delete [] rawdata;
  }

  void geometric_wavefront::print(ostream & os, const char * prefix) const {
    this->geometric_wavefront_header::print(os, prefix);
  }

  void geometric_wavefront::propagate(const three_point & o, const three_vector & n) {
    double distance;
    for(int i=0; i<nrays; i++){
      try{
	rays[i].propagate(rays[i].distance_along_normal(o,n));}
      catch(...){
	rays[i].set_lost(true);
      }
    }
  }

  bool operator==(const geometric_wavefront & gwf1, const geometric_wavefront & gwf2) {
    if(!operator==(static_cast<geometric_wavefront_header>(gwf1),
		   static_cast<geometric_wavefront_header>(gwf2)))
      return(false);
    if(gwf1.nrays!=gwf2.nrays) return(false);
    for(int i=0; i<gwf1.nrays; i++)
      if(gwf1.rays[i]!=gwf2.rays[i]) return(false);
    return(true);
  }

  int geometric_ray::verbose_level = 0;

  geometric_ray::geometric_ray() :
    three_vector(0,0,1,three_frame()) {
    pathlength = 0;
    lost = false;
  }

  geometric_ray::geometric_ray(const geometric_ray & gr) {
    this->operator=(gr);
  }

  geometric_ray::geometric_ray(const three_point & tp, const three_vector & tv, double in_pathlength) {
    this->three_point::operator=(tp);
    this->three_vector::operator=(tv);
    pathlength = in_pathlength;
    lost = false;
  }

  geometric_ray & geometric_ray::operator=(const geometric_ray & gr) {
    if(this==&gr)
      return(*this);
    this->three_point::operator=(gr);
    this->three_vector::operator=(gr);
    pathlength = gr.pathlength;
    lost = gr.lost;
    return(*this);
  }

  void geometric_ray::print(ostream & os, const char * prefix) const {
    int w = 15;
    os.setf(ios::right, ios::adjustfield); 
    this->three_point::print(cout, three_frame());
    this->three_vector::print(cout, three_frame());
    os << setw(w) << this->pathlength << endl;
    os << setw(w) << this->lost << endl;
  }

  void geometric_ray::propagate(double distance) {
    three_vector propagation_vector = *this;
    propagation_vector *= distance;
    this->three_point::operator+=(propagation_vector);
  }

  bool operator==(const geometric_ray & gr1, const geometric_ray & gr2) {
    if(!operator==(static_cast<const three_point>(gr1), static_cast<const three_point>(gr2))) return(false);
    if(!operator==(static_cast<const three_vector>(gr1), static_cast<const three_vector>(gr2))) return(false);
    if(gr1.pathlength != gr2.pathlength) return(false);
    if(gr1.lost != gr2.lost) return(false);
    return(true);
  }

  bool operator!=(const geometric_ray & gr1, const geometric_ray & gr2) {
    if(operator==(gr1,gr2)) return(false);
    return(true);
  }

  double geometric_ray::distance_along_normal(const three_point & p, const three_vector & n) const {
    throw(string("geometric_ray::distance_along_normal"));
    /*
    if(dot_product(*this, n)==0) 
      throw(string("geometric_ray::distance_along_normal"));
    return(get_distance_along_normal(static_cast<const three_point & >(*this), 
				     static_cast<const three_vector & >(*this), p, n));
    */
  }

}
