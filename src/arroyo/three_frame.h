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

#ifndef THREE_FRAME_H
#define THREE_FRAME_H

#include "AO_cpp.h"
#include "iofits.h"
#include "three_point.h"
#include "three_vector.h"

namespace Arroyo {

  using std::ostream;

  ///
  /// A class to represent a frame in three dimensional space.
  ///

  class three_frame :
    public three_point {

    protected:

    /// Nine components defined 
    /// relative to a global coordinate 
    /// system
    double basis[9];

    public:
  
    ///////////////////////////////////////////
    ///  Null constructor
    three_frame();

    ///////////////////////////////////////////
    ///  Copy constructor
    three_frame(const three_frame & tp);

    ///////////////////////////////////////////
    ///  Construct from a point and 3 vectors
    ///  These vectors must be orthogonal, but
    ///  need not be normalized
    three_frame(const three_point & tp, const three_vector & x, 
		const three_vector & y, const three_vector & z);

    ///////////////////////////////////////////
    ///  Destructor
    ~three_frame(){};

    ///////////////////////////////////////////
    ///  Operator = 
    three_frame & operator=(const three_frame & tp);

    ///////////////////////////////////////////
    ///  Read from iofits object
    void read(const iofits & iof);

    ///////////////////////////////////////////
    ///  Write to iofits object
    void write(iofits & iof) const;

    ///////////////////////////////////////////
    ///  Print components wrt default frame
    void print(ostream & os, const char * prefix="",
    			long precision = 6) const;

    ///////////////////////////////////////////
    ///  Print components wrt frame tf
    void print(ostream & os, const three_frame & tf,
    			const char * prefix="", long precision = 6) const;

    ///////////////////////////////////////////
    ///  Get x basis vector wrt frame f
    three_vector x() const;

    ///////////////////////////////////////////
    ///  Get y basis vector wrt frame f
    three_vector y() const;

    ///////////////////////////////////////////
    ///  Get z basis vector wrt frame f
    three_vector z() const;

    ///////////////////////////////////////////
    ///  Check if this frame is a right-handed frame.
    ///  That is, whether X x Y = Z
    bool right_handed_frame() const;

    ///////////////////////////////////////////
    ///  Operator==
    friend bool operator==(const three_frame & tf1, const three_frame & tf2);

    ///////////////////////////////////////////
    ///  Construct a three vector from coordinates
    ///  specified relative to a frame tf
    friend three_vector::three_vector(double x, double y, double z,
    						const three_frame & tf);

    /// The precision to use in comparison calculations
    static double precision;

    /// A verbose_level for printing messages
    static int verbose_level;

  };

  ///////////////////////////////////////////
  ///  Operator!=
  bool operator!=(const three_frame & tf1, const three_frame & tf2);
}

#endif
