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

#include <math.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include "fits_factory.h"
#include "fits_header_data.h"
#include "conic_section.h"

using namespace std;

namespace Arroyo {

  int conic_section::verbose_level = 0;

  namespace factory_register {
    const fits_keyval_set & get_conic_section_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "conic section"));
      return *fkvs;
    }
    
    AO_sim_base * create_conic_section(const iofits & iof) {
      return new conic_section(iof);
    }
  }

  const bool conic_section::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_conic_section_keyval_set(), 
				      factory_register::create_conic_section);


  conic_section::conic_section(const conic_section & conic_sctn) {
    this->operator=(conic_sctn);
  }

  conic_section::conic_section(const char * filename){
    this->read(filename);
  }

  conic_section::conic_section(const iofits & iof){
    this->read(iof);
  }

  conic_section::conic_section(const three_point & vtx, 
			       const three_point & focus, 
			       double eccty) {
    if((focus-vtx).length()<three_frame::precision){
      cerr << "conic_section::conic_section error - vertex and focal three_points supplied to constructor are identical\n";
      throw(string("conic_section::conic_section"));
    }

    this->vertex = vtx;
    this->eccentricity = eccty;
    this->focus = focus;
  }
 
  conic_section & conic_section::operator=(const conic_section & cmr) {
    if(this==&cmr) return(*this);
    this->vertex = cmr.vertex;
    this->eccentricity = cmr.eccentricity;
    this->focus = cmr.focus;
    return(*this);
  }

  void conic_section::read(const char * filename) {
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "conic_section::read - "
	   << "error opening file " << filename << endl;
      throw(string("conic_section::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "conic_section::read - "
	   << "error reading "
	   << this->unique_name() << " from file " 
	   << filename << endl;
      throw(string("conic_section::read"));
    }
  }

  void conic_section::read(const iofits & iof) {
    if(!iof.key_exists("TYPE")){
      cerr << "conic_section::read error - "
	   << "no power law type specified\n";
      throw(string("conic_section::read"));
    } else {
      string type, comment;
      iof.read_key("TYPE", type, comment);
      if(type!=this->unique_name()){
	cerr << "conic_section::read error - outer scale of type " 
	     << type << " rather than type "
	     << this->unique_name() << endl;
	throw(string("conic_section::read"));
      } else {
	iof.read_key("ECCTY", this->eccentricity, comment);
	this->vertex.read(iof);
	this->focus.read(iof);
      }
    }

  }
 
  void conic_section::write(const char * filename) const {
    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "conic_section::write - "
	   << "error opening file " << filename << endl;
      throw(string("conic_section::write"));
    }

    try{this->write(iof);}
    catch(...){
      cerr << "conic_section::write - "
	   << "error writing "
	   << this->unique_name() << " to file " 
	   << filename << endl;
      throw(string("conic_section::write"));
    }
  }

  void conic_section::write(iofits & iof) const {
    Arroyo::fits_header_data<char> tmphdr;
    tmphdr.write(iof);

    string type = this->unique_name();
    string comment = "object type";
    iof.write_key("TYPE", type, comment);

    comment = "eccentricity";
    iof.write_key("ECCTY", this->eccentricity, comment);
    this->vertex.write(iof);
    this->focus.write(iof);
  }

  void conic_section::print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(std::ios::left, std::ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << this->unique_name()
       << "/" << "object type" << endl;
    os << prefix << "ECCTY      = " << setw(vlspc) << eccentricity
       << "/" << "conic eccentricity" << endl;
    this->vertex.print(os, prefix);
    this->focus.print(os, prefix);
  }

  double conic_section::get_local_curvature(three_point & tp) const {

    three_vector tv = tp - this->vertex;
    three_vector optical_axis = this->focus - this->vertex;
    optical_axis = (1/optical_axis.length())*optical_axis;

    double optical_axis_offset = 
      cross_product((tp - this->vertex), optical_axis).length();

    double four_p_squared = .25*this->get_latus_rectum()*this->get_latus_rectum();

    /*
    cout << "conic_section::get_local_curvature - 4p^{2} " << four_p_squared
	 << " optical_axis_offset " << optical_axis_offset 
	 << " ecc " << this->eccentricity
	 << endl;
    */

    // Stavroudis expression for curvature, p 90
    double curvature = 
      four_p_squared / pow((four_p_squared + pow(this->eccentricity*optical_axis_offset,2.)),3/2.);

    /*
    cout << "conic_section::get_local_curvature - returning "
	 << curvature << endl;
    */

    return(curvature);
  }

  bool conic_section::on_conic(const three_point & tp) const {
    if((tp-this->vertex).length()<three_frame::precision)
      return(true);

    three_vector axial_vector = this->focus - this->vertex;
    three_vector tv = tp - this->vertex;

    double axial_offset = 
      (tv - dot_product(tv,axial_vector)*(1/axial_vector.length())*axial_vector).length();

    if(fabs(axial_offset - this->get_latus_rectum()*axial_vector.length() + 
       (1 - this->eccentricity*this->eccentricity)*axial_vector.length_squared()) < three_frame::precision)
      return(true);

    return(false);
  }
 
  void conic_section::raytrace(const three_point & tp, 
			       const three_vector & tv,
			       double & distance_to_conic,
			       three_point & point_of_intersection,
			       three_vector & conic_unit_normal,
			       double & R_squared,
			       double & V) const {

    // This function is based on Chapter 6 of Stavroudis,
    // "The Optics of Rays, Wavefronts and Caustics"

    double length = tv.length();
    if(length<three_frame::precision){
      cerr << "conic_section::get_point_of_intersection error - "
	   << "null three vector supplied to this function\n";
      throw(string("conic_section::get_point_of_intersection"));
    }

    // The three vector W in Stavroudis is the vector from the
    // vertex of the conic to the three point
    // With this definition, t = 0 in Stavroudis' equations
    three_vector W = tp - this->vertex;

    // This is N in Stavroudis
    three_vector normal_tv = (1/length)*tv;

    // This is A in Stavroudis, which denotes the optical axis
    three_vector optical_axis = this->focus - this->vertex;
    optical_axis = (1/optical_axis.length())*optical_axis;
    
    // This is the vertex curvature, which Stavroudis calls c
    double vertex_curvature = 2/this->get_latus_rectum();

    // Stavroudis, eq. VI-35 a
    R_squared = 
      (1-this->eccentricity*this->eccentricity*dot_product(normal_tv,optical_axis)*dot_product(normal_tv,optical_axis))*
      (1-(cross_product(optical_axis,normal_tv)+vertex_curvature*cross_product(normal_tv, W)).length_squared()) +
      this->eccentricity*this->eccentricity*
      pow((1-dot_product(cross_product(optical_axis,normal_tv),
			 (cross_product(optical_axis,normal_tv)+
			  vertex_curvature*cross_product(normal_tv,W)))),2.);

    if(conic_section::verbose_level){
      cout << endl << endl;
      tp.print(cout, "conic_section::raytrace - init tp ");
      tv.print(cout, "conic_section::raytrace - init tv ");
      cout << "conic_section::raytrace - eccentricity " << this->eccentricity << endl;
      cout << "conic_section::raytrace - vertex curvature " << vertex_curvature << endl;
      W.print(cout, "conic_section::raytrace - W ");
      normal_tv.print(cout, "conic_section::raytrace - N ");
      this->vertex.print(cout, "conic_section::raytrace - vertex ");
      this->focus.print(cout, "conic_section::raytrace - focus ");
      optical_axis.print(cout, "conic_section::raytrace - A ");
      cout << "conic_section::raytrace - R_squared "
	   << R_squared
	   << endl;
    }

    // Stavroudis, eq. VI-35 b
    /*
    V = 
      dot_product(optical_axis,normal_tv) -
      vertex_curvature*dot_product(normal_tv,W) +
      sqrt(R_squared) +
      vertex_curvature*this->eccentricity*this->eccentricity*
      dot_product(optical_axis,W)*
      dot_product(optical_axis,normal_tv);
    */

    double positive_root_V = 
      dot_product(optical_axis,normal_tv) -
      vertex_curvature*dot_product(normal_tv,W) +
      vertex_curvature*this->eccentricity*this->eccentricity*
      dot_product(optical_axis,W)*
      dot_product(optical_axis,normal_tv);

    // Here there are two roots for R_squared, corresponding to the
    // two points of intersection with the conic.  We would like to
    // take the point of intersection that is closest to tp.  Since V
    // appears only in the denominator of distance_to_conic, this
    // corresponds to maximizing the magnitude of V

    double negative_root_V = positive_root_V - sqrt(R_squared);
    positive_root_V += sqrt(R_squared);

    V = fabs(positive_root_V) > fabs(negative_root_V) ? positive_root_V : negative_root_V;

    if(conic_section::verbose_level)
      cout << "conic_section::raytrace - A*N "
	   << dot_product(optical_axis,normal_tv)
	   << "\tN*W/R_{c} " 
	   << vertex_curvature*dot_product(normal_tv,W)
	   << "\tsqrt(R_squared) " 
	   << sqrt(R_squared)
	   << "\tlast term " 
	   << vertex_curvature*this->eccentricity*this->eccentricity*
	dot_product(optical_axis,W)*
	dot_product(optical_axis,normal_tv)
	   << "\t+root "
	   << positive_root_V
	   << "\t-root " 
	   << negative_root_V
	   << "\tV "
	   << V
	   << endl;


    // Stavroudis, eq. VI-35 c (distance_to_conic is lambda_bar) For a
    // parabola, V vanishes when tp is on the optical axis.  This
    // special case requires one not to divide out V, but to instead
    // use the unnormalized expression for the distance to the conic
    // in Stavroudis eq. VI-33.  And this expression is also singular
    // when the propagation direction lies along the optical axis.
    // So there are two special cases considered below.

    // Well, its worse than this.  The special cases depend on tp,
    // and some cases for eccentricity==1 are valid.  (e.g. plane
    // wave incident on the parabola.  Try to exclude cases for which
    // V==0

    if(this->on_conic(tp)){
      distance_to_conic = 0;

      if(conic_section::verbose_level)
	cout << "conic_secton::raytrace - on conic - distance 0\n";

    } else if(fabs(V)>three_frame::precision){
      distance_to_conic = 
	(1/V)*(vertex_curvature*W.length_squared() -
	       2*dot_product(optical_axis,W) -
	       vertex_curvature*this->eccentricity*this->eccentricity*pow(dot_product(optical_axis,W),2.));

      if(conic_section::verbose_level)
	cout << "conic_section::raytrace - case a "
	     << W.length()
	     << "\t"
	     << vertex_curvature
	     << "\t"
	     << vertex_curvature*W.length_squared()
	     << "\t" 
	     << 2*dot_product(optical_axis,W)
	     << "\t" 
	     << vertex_curvature*this->eccentricity*this->eccentricity*pow(dot_product(optical_axis,W),2.)
	     << "\t" 
	     << distance_to_conic 
	     << endl;

    } else if(cross_product(optical_axis,W).length()<three_frame::precision && 
	      cross_product(optical_axis,normal_tv).length()<three_frame::precision) {

      distance_to_conic = dot_product(W,normal_tv);

      if(conic_section::verbose_level)
	cout << "conic_section::raytrace - case b " << distance_to_conic << endl;

    } else {

      distance_to_conic = 
	(dot_product(optical_axis,normal_tv) - 
	vertex_curvature*dot_product(normal_tv,W) +
	vertex_curvature*this->eccentricity*this->eccentricity*dot_product(optical_axis,W)*dot_product(optical_axis,normal_tv) - 
	 sqrt(R_squared))/
	(vertex_curvature*(1-this->eccentricity*this->eccentricity*pow(dot_product(normal_tv,optical_axis),2.)));

      if(conic_section::verbose_level)
	cout << "conic_section::raytrace - case c "
	     << "\tV "
	     << V
	     << "\tdenom "
	     << vertex_curvature*(1-this->eccentricity*this->eccentricity*pow(dot_product(normal_tv,optical_axis),2.))
	     << "\t"
	     << distance_to_conic 
	     << endl;
    }

    // Stavroudis, eq. VI-35 d
    three_vector R_bar = W + distance_to_conic * normal_tv;

    point_of_intersection = tp + distance_to_conic*normal_tv;

    // Stavroudis, eq. VI-35 e (conic_unit_normal is N_bar)
    // Note discrepancy in denominators of eqs. VI-31 and VI-35 e

    conic_unit_normal = 
      (1/sqrt(1 + pow(vertex_curvature*this->eccentricity*dot_product(R_bar,optical_axis),2.)))*
      (-vertex_curvature*R_bar + 
       optical_axis + 
       vertex_curvature*this->eccentricity*this->eccentricity*
       dot_product(optical_axis,R_bar)*optical_axis);

    /*
    conic_unit_normal = 
      (1/sqrt(1 + pow(vertex_curvature*this->eccentricity*dot_product(R_bar,optical_axis),2.)))*
      (-vertex_curvature*R_bar + 
       optical_axis + 
       vertex_curvature*this->eccentricity*this->eccentricity*
       dot_product(optical_axis,R_bar)*optical_axis);
    */


    // Stavroudis leaves this unnormalized in Equation VI-35 e
    // Here we explicitly normalize it
    conic_unit_normal = (1/conic_unit_normal.length())*conic_unit_normal;

    if(conic_section::verbose_level){
      point_of_intersection.print(cout, "conic_section::raytrace - poi ");
      conic_unit_normal.print(cout, "conic_section::raytrace - cun ");
      cout << "conic_section::raytrace - distance to conic " 
	   << distance_to_conic 
	   << endl;
      cout << "conic_section::raytrace - R_squared " 
	   << R_squared
	   << endl;
      cout << "conic_section::raytrace - V " 
	   << V
	   << endl;
    }
  }

}
