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

#ifndef CONIC_SECTION_H
#define CONIC_SECTION_H

#include "AO_sim_base.h"
#include "three_frame.h"

namespace Arroyo {

  ///
  /// A class to represent a conic section
  ///
  class conic_section :
    virtual public AO_sim_base {

  private:
    
    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("conic section"));};

  protected:

    // Plane wave flag
    bool plane_wave;

    // The vertex of the conic section
    three_point vertex;

    // The focus of the conic section closest to the vertex
    three_point focus;

    // The eccentricity:
    // For a spherical mirror, eccentricity = 0
    // For an elliptical mirror, 0 < eccentricity < 1
    // For a parabolic mirror, eccentricity = 1
    // For a hyperbolic mirror, eccentricity > 1
    double eccentricity;

  public:

    ///////////////////////////////////////////
    ///  Null constructor
    conic_section(){
      this->eccentricity = -1;
    };

    ///////////////////////////////////////////
    ///  Copy constructor
    conic_section(const conic_section & cmr);

    ///////////////////////////////////////////
    ///  Construct from file
    conic_section(const char * filename);

    ///////////////////////////////////////////
    ///  Construct from iofits object
    conic_section(const iofits & iof);

    ///////////////////////////////////////////
    ///  Construct from the bits
    /// 
    ///  The vertex of the conic lies along the 
    ///  axis of symmetry. 
    ///
    ///  The focus of the conic is the focal point
    ///  closest to the vertex.  
    ///
    ///  The eccentricity of the conic is as follows:
    ///
    ///     For a spherical conic, eccentricity = 0
    ///     For an elliptical conic, 0 < eccentricity < 1
    ///     For a parabolic conic, eccentricity = 1
    ///     For a hyperbolic conic, eccentricity > 1
    ///
    conic_section(const three_point & vertex, 
		  const three_point & focus, 
		  double eccty);

    ///////////////////////////////////////////
    ///  Destructor
    ~conic_section(){};

    ///////////////////////////////////////////
    ///  Operator = 
    conic_section & operator=(const conic_section & cmr);

    ///////////////////////////////////////////
    ///  Read from file
    void read(const char * filename);

    ///////////////////////////////////////////
    ///  Read from an iofits object
    void read(const iofits & iof);
 
    ///////////////////////////////////////////
    ///  Write to file
    void write(const char * filename) const;

    ///////////////////////////////////////////
    ///  Write to an iofits object
    void write(iofits & iof) const;

    ///////////////////////////////////////////
    ///  Print
    void print(ostream & os, const char * prefix="") const;

    ///////////////////////////////////////////
    ///  Get eccentricity
    double get_eccentricity() const {
      return(this->eccentricity);
    };

    ///////////////////////////////////////////
    ///  Get vertex
    three_point get_vertex() const {
      return(this->vertex);
    };

    ///////////////////////////////////////////
    ///  Get focal point nearest vertex
    three_point get_near_focus() const {
      return(this->focus);
    };

    ///////////////////////////////////////////
    ///  Get focal point farthest from vertex
    three_point get_far_focus() const {
      if(this->eccentricity==1){
	cerr << "conic_section::get_far_focus error - this instance is a parabola, and the focus is at infinity\n";
	throw(string("conic_section::get_far_focus"));
      }
      double dist = this->get_latus_rectum()/2./(1-this->eccentricity);
      three_vector unit_vector = this->focus - this->vertex;
      unit_vector = (1/unit_vector.length())*unit_vector;
      return(this->vertex + dist*unit_vector);
    };

    ///////////////////////////////////////////
    ///  Get latus rectum
    double get_latus_rectum() const {
      return(2*(this->focus - this->vertex).length()*(1+this->eccentricity));
    };

    ///////////////////////////////////////////
    ///  Get radius of curvature at vertex
    double get_local_curvature(three_point & tp) const;

    ///////////////////////////////////////////
    ///  State whether the three point is on the conic surface
    bool on_conic(const three_point & tp) const;

    ///////////////////////////////////////////
    ///  Raytrace a line extending from three_point tp in the
    ///  direction of the three_vector tv and the conic section.  If
    ///  there is no intersection point, this function throws an error
    ///
    ///  distance_to_conic is the distance from tp to the conic point
    ///  of intersection in whatever units were used to define the
    ///  focus and vertex of the conic section
    ///
    ///  point_of_intersection is the point at which the ray
    ///  intersects the conic surface
    ///
    ///  conic_unit_normal is the unit vector normal to the conic
    ///  surface at the point of intersection
    ///
    ///  The other two arguments are defined in Stavroudis 
    void raytrace(const three_point & tp, 
		  const three_vector & tv,
		  double & distance_to_conic,
		  three_point & point_of_intersection,
		  three_vector & conic_unit_normal,
		  double & R_squared,
		  double & V) const;

    ///////////////////////////////////////////
    ///  Verbose level
    static int verbose_level;

  };

}

#endif
  
