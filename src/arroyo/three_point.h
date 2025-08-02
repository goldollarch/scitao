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

#ifndef THREE_POINT_H
#define THREE_POINT_H

#include "AO_cpp.h"
#include "iofits.h"
#include "three_transformation.h"

namespace Arroyo {

  // forward declarations
  class three_vector;
  class three_frame;
 
  using std::ostream;

  ///
  /// Aclass to represent a point in three dimensional space.
  ///

  class three_point {

  protected:

    /// X coordinate defined relative to global coordinate system
    double x_;

    /// Y coordinate defined relative to global coordinate system
    double y_;

    /// Z coordinate defined relative to global coordinate system
    double z_;

    ///////////////////////////////////////////
    ///  Construct from coordinates wrt underlying frame
    three_point(double x, double y, double z);

  public:
  
    ///////////////////////////////////////////
    ///  Null constructor
    three_point(){x_=y_=z_=0;};

    ///////////////////////////////////////////
    ///  Copy constructor
    three_point(const three_point & tp);

    ///////////////////////////////////////////
    ///  Construct from coordinates wrt frame f
    three_point(double x, double y, double z, const three_frame & tf);

    ///////////////////////////////////////////
    ///  Destructor
    ~three_point(){};

    ///////////////////////////////////////////
    ///  Operator = 
    three_point & operator=(const three_point & tp);

    ///////////////////////////////////////////
    ///  Read from iofits object
    void read(const iofits & iof);

    ///////////////////////////////////////////
    ///  Write to iofits object
    void write(iofits & iof) const;

    ///////////////////////////////////////////
    ///  Print coordinates wrt default frame
    void print(ostream & os, const char * prefix = "",
    					long precision = 6) const;

    ///////////////////////////////////////////
    ///  Print coordinates wrt frame tf
    void print(ostream & os, const three_frame & tf,
    			const char * prefix = "", long precision = 6) const;

    ///////////////////////////////////////////
    ///  Get x coordinate wrt frame f
    double x(const three_frame & tf) const;

    ///////////////////////////////////////////
    ///  Get y coordinate wrt frame f
    double y(const three_frame & tf) const;

    ///////////////////////////////////////////
    ///  Get z coordinate wrt frame f
    double z(const three_frame & tf) const;

    ///////////////////////////////////////////
    ///  Add a vector to a point to get another point
    three_point & operator+=(const three_vector & tv);

    ///////////////////////////////////////////
    ///  Subtract a vector from a point to get another point
    three_point & operator-=(const three_vector & tv);

    ///////////////////////////////////////////
    ///  Subtract two points to get a vector 
    friend three_vector operator-(const three_point & tp1,
    					const three_point & tp2);

    ///////////////////////////////////////////
    ///  Friend member function of three_transformation
    ///  that will effect that transformation on a three_point
    friend void three_transformation::transform(three_point & tp) const;

    ///////////////////////////////////////////
    ///  Operator==
    friend bool operator==(const three_point & tp1,
    					const three_point & tp2);

    /// A verbose_level for printing messages
    static int verbose_level;

  };

  ///////////////////////////////////////////
  ///  Operator!=
  bool operator!=(const three_point & tp1, const three_point & tp2);

  three_point operator+(const three_point & tp, const three_vector & tv);
  three_point operator-(const three_point & tp, const three_vector & tv);

}

#endif
