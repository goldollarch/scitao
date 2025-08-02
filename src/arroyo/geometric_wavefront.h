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

#ifndef GEOMETRIC_WAVEFRONT_H
#define GEOMETRIC_WAVEFRONT_H

#include <complex>
#include "three_point.h"
#include "three_vector.h"
#include "wavefront.h"

namespace Arroyo {

  using std::ostream;
  using std::string;

  class optic;
  class geometric_ray;

  ///
  ///     A class to represent a geometric wavefront.
  ///

  class geometric_wavefront :
    virtual public geometric_wavefront_header, 
    virtual public wavefront {

    private:

    static const bool factory_registration;

    /// The optics are going to transform the wavefront,
    /// and any practical, efficient implementation requires access to
    /// the protected wfdata array.  So we make the base class of the optic
    /// inheritance hierarchy a friend of wavefront, and then write 
    /// protected member functions in optic that return protected data
    /// from wavefront.  
    friend class one_to_one_optic;

    /// The optics are going to transform the wavefront,
    /// and any practical, efficient implementation requires access to
    /// the protected wfdata array.  So we make the base class of the optic
    /// inheritance hierarchy a friend of wavefront, and then write 
    /// protected member functions in optic that return protected data
    /// from wavefront.  
    friend class one_to_many_optic;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("geometric wavefront"));};

     protected:

    /// The ray array
    geometric_ray * rays;

    public:

    ///////////////////////////////////////////
    ///  Null constructor
    geometric_wavefront();

    ///////////////////////////////////////////
    ///  Copy constructor
    geometric_wavefront(const geometric_wavefront & wf);

    ///////////////////////////////////////////
    ///  Construct from a file
    geometric_wavefront(const char * filename);

    ///////////////////////////////////////////
    ///  Construct from an iofits object
    geometric_wavefront(const iofits & iof);

    ///////////////////////////////////////////
    ///  Construct from a header
    geometric_wavefront(const geometric_wavefront_header & gwfh);

    ///////////////////////////////////////////
    ///  Null destructor
    ~geometric_wavefront();

    ///////////////////////////////////////////
    ///  Operator = 
    geometric_wavefront & operator=(const geometric_wavefront & wf);

    ///////////////////////////////////////////
    ///  Read geometric_wavefront from a file
    void read(const char * filename);

    ///////////////////////////////////////////
    ///  Read geometric_wavefront from an iofits object
    void read(const iofits & iof);

    ///////////////////////////////////////////
    ///  Write geometric_wavefront to a file
    void write(const char * filename) const;

    ///////////////////////////////////////////
    ///  Write geometric_wavefront to an iofits object
    void write(iofits & iof) const;

    ///////////////////////////////////////////
    ///  Print information about the geometric_wavefront 
    void print(ostream & os, const char * prefix) const;

    ///////////////////////////////////////////
    ///  Propagate geometric wavefront to a plane
    ///  specified by origin o and normal vector n
    void propagate(const three_point & o, 
		   const three_vector & n);

    ///////////////////////////////////////////
    ///  Friend operator==  for geometric_wavefronts
    friend bool operator==(const geometric_wavefront & wf1, const geometric_wavefront & wf2);
  };

  bool operator!=(const geometric_wavefront & wf1, const geometric_wavefront & wf2);

///
/// A class to represent a geometric ray.
///

  class geometric_ray : 
    public three_point, 
    public three_vector {

    protected:
  
    /// The geometric path length traversed by the ray 
    double pathlength;

    /// A flag to indicate if the ray has been lost, via aperture
    /// or propagation orthogonal to some surface
    bool lost;

    public:

    ///////////////////////////////////////////
    ///  Null constructor
    geometric_ray();

    ///////////////////////////////////////////
    ///  Copy constructor
    geometric_ray(const geometric_ray & wf);

    ///////////////////////////////////////////
    ///  Construct from the bits
    geometric_ray(const three_point & tp, 
		  const three_vector & tv,
		  double in_pathlength);

    ///////////////////////////////////////////
    ///  Null destructor
    ~geometric_ray(){};

    ///////////////////////////////////////////
    ///  Operator = 
    geometric_ray & operator=(const geometric_ray & wf);

    ///////////////////////////////////////////
    ///  Print information about the geometric_ray 
    void print(ostream & os, const char * prefix) const;

    ///////////////////////////////////////////
    ///  Get the pathlength
    double get_pathlength() const {return(pathlength);};

    ///////////////////////////////////////////
    ///  Set the pathlength
    void set_pathlength(double plength) {pathlength = plength;};

    ///////////////////////////////////////////
    ///  Check if the ray is lost
    bool is_lost() const {return(lost);};

    ///////////////////////////////////////////
    ///  Set the lost flag 
    ///  lost (b=true) 
    ///  found (b=false)
    void set_lost(bool b) {lost = b;};

    ///////////////////////////////////////////
    ///  Propagate using geometric ray propagation
    void propagate(double distance);

    ///////////////////////////////////////////
    ///  Find the distance along the propagation 
    ///  direction to a plane specified by a point
    ///  p and a normal vector n
    double distance_along_normal(const three_point & p, const three_vector & n) const;

    ///////////////////////////////////////////
    ///  Verbose level
    static int verbose_level;

    ///////////////////////////////////////////
    ///  Friend operator ==  for geometric_rays
    friend bool operator ==(const geometric_ray & wf1, const geometric_ray & wf2);

    ///////////////////////////////////////////
    ///  Friendship for geometric_wavefront::read
    ///  to implement efficient storage of rays
    friend void geometric_wavefront::read(const iofits & iof);

  };

  ///////////////////////////////////////////
  ///  Operator !=  for geometric_rays
  bool operator!=(const geometric_ray & ray1, const geometric_ray & ray2);

}

#endif
