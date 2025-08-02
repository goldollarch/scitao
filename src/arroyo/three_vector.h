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

#ifndef THREE_VECTOR_H
#define THREE_VECTOR_H

#include <math.h>
#include "AO_cpp.h"
#include "iofits.h"
#include "three_point.h"
#include "three_transformation.h"

namespace Arroyo {

  using std::ostream;

  ///
  /// A class to represent a vector in three dimensional space.
  ///

  class three_vector {

  protected:

    /// X components defined relative to global coordinate system
    double x_;

    /// Y components defined relative to global coordinate system
    double y_;

    /// Z components defined relative to global coordinate system
    double z_;

    ///////////////////////////////////////////
    ///  Construct from components wrt global frame 
    three_vector(double x, double y, double z);

  public:
  
    ///////////////////////////////////////////
    ///  Null constructor
    three_vector(){x_=y_=z_=0;};

    ///////////////////////////////////////////
    ///  Copy constructor
    three_vector(const three_vector & tv);

    ///////////////////////////////////////////
    ///  Construct from components wrt frame f
    three_vector(double x, double y, double z, const three_frame & tf);

    ///////////////////////////////////////////
    ///  Destructor
    ~three_vector(){};

    ///////////////////////////////////////////
    ///  Operator = 
    three_vector & operator=(const three_vector & tv);

    ///////////////////////////////////////////
    ///  Read from iofits object
    void read(const iofits & iof);

    ///////////////////////////////////////////
    ///  Write to iofits object
    void write(iofits & iof) const;

    ///////////////////////////////////////////
    ///  Print components wrt default frame
    void print(ostream & os, const char * prefix="", long precision = 6) const;

    ///////////////////////////////////////////
    ///  Print components wrt frame f
    void print(ostream & os, const three_frame & tf,
    			const char * prefix="", long precision = 6) const;

    ///////////////////////////////////////////
    ///  Get x component wrt frame f
    double x(const three_frame & tf) const;

    ///////////////////////////////////////////
    ///  Get y component wrt frame f
    double y(const three_frame & tf) const;

    ///////////////////////////////////////////
    ///  Get z component wrt frame f
    double z(const three_frame & tf) const;

    ///////////////////////////////////////////
    ///  Returns the length of the vector
    double length() const {return(sqrt(x_*x_+y_*y_+z_*z_));};

    ///////////////////////////////////////////
    ///  Returns the squared length of the vector
    ///  This is faster than length() by a sqrt call
    double length_squared() const {return(x_*x_+y_*y_+z_*z_);};

    ///////////////////////////////////////////
    ///  Add a vector to this one 
    three_vector & operator+=(const three_vector & tv);

    ///////////////////////////////////////////
    ///  Subtract a vector from this one 
    three_vector & operator-=(const three_vector & tv);

    ///////////////////////////////////////////
    ///  Multiply vector by a scalar
    three_vector & operator*=(double d);

    ///////////////////////////////////////////
    ///  Subtract two points to get a vector 
    friend three_vector operator-(const three_point & tp1,
    					const three_point & tp2);

    ///////////////////////////////////////////
    ///  Dot product of two vectors
    friend double dot_product(const three_vector & tv1,
    					const three_vector & tv2);

    ///////////////////////////////////////////
    ///  Cross product of two vectors
    ///  Specifically, it returns tv1 x tv2
    friend three_vector cross_product(const three_vector & tv1,
    					const three_vector & tv2);

    ///////////////////////////////////////////
    ///  Add a vector to a point to get another point
    friend three_point & three_point::operator+=(const three_vector & tv);

    ///////////////////////////////////////////
    ///  Subtract a vector from a point to get another point
    friend three_point & three_point::operator-=(const three_vector & tv);

    ///////////////////////////////////////////
    ///  Friend member function of three_transformation
    ///  that will effect that transformation on a three_vector
    friend void three_transformation::transform(three_vector & tv) const;

    ///////////////////////////////////////////
    ///  Operator==
    friend bool operator==(const three_vector & tv1, const three_vector & tv2);

    /// A verbose_level for printing messages
    static int verbose_level;

  };

  ///////////////////////////////////////////
  ///  Operator!=
  bool operator!=(const three_vector & tv1, const three_vector & tv2);

  three_vector operator+(const three_vector & tv1, const three_vector & tv2);
 
  three_vector operator-(const three_vector & tv1, const three_vector & tv2);

  three_vector operator*(double d, const three_vector & tv);
 
  three_vector operator*(const three_vector & tv, double d);
 
  ///////////////////////////////////////////
  ///  Return the dot product of two vectors
  inline double dot_product(const three_vector & tv1, const three_vector & tv2) {
    return(tv1.x_ * tv2.x_ + tv1.y_ * tv2.y_ + tv1.z_ * tv2.z_);
  }

  ///////////////////////////////////////////
  ///  Return the cross product of two vectors
  inline three_vector cross_product(const three_vector & tv1,
  					const three_vector & tv2) {
    return(three_vector(tv1.y_ * tv2.z_ - tv1.z_ * tv2.y_,
			tv1.z_ * tv2.x_ - tv1.x_ * tv2.z_,
			tv1.x_ * tv2.y_ - tv1.y_ * tv2.x_));
  }


  ///////////////////////////////////////////
  ///  Return the projection of the three vector
  ///  tv into a plane orthogonal to the vector
  ///  nrml, along the direction specified
  ///  by the three vector drctn.  If drctn==normal
  ///  then the parallel projection yields the component
  ///  of tv in the plane orthogonal to nrml.
  three_vector parallel_projection(const three_vector & tv, 
				   const three_vector & nrml, 
				   const three_vector & drctn);
}
#endif

