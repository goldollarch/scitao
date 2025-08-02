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

#ifndef REGION_BASE_H
#define REGION_BASE_H

#include <iostream>
#include <string>
#include "three_point.h"

namespace Arroyo {

  class three_point;

  using std::string;
  using std::ostream;
  using std::vector;

  ///
  /// virtual base class to specify a region in space
  ///

  class region_base {

  public:
  
    ////////////////////////////
    ///  Null constructor
    region_base(){};

    ////////////////////////////
    ///  Virtual destructor
    virtual ~region_base(){};

    ////////////////////////////
    ///  Virtual print
    virtual void print(ostream & os, const char * prefix = "",
    					long precision = 6) const = 0;

    ////////////////////////////
    ///  Print
    virtual void print(ostream & os, const three_frame & tf,
    		const char * prefix = "", long precision = 6) const = 0;
    
    //////////////////////////// 
    ///  At some point in the future
    ///  we may want to add an alexandrescu-like
    ///  double dispatcher here, which is why
    ///  this base class exists.

    /// Verbose level for printing messages
    static int verbose_level;

  };

  ///
  /// A class to specify a 2d rectangular region in space
  ///

  class rectangular_region :
    public region_base {

    /// the 4 three_points that specify the corners of the region
    vector<three_point> corners;

    ////////////////////////////////////////////////////////////////////
    /// This function will find the overlap status of the two rectangles
    /// There are a number of return values:
    /// no overlap              - no overlap between the rectangles
    /// match                   - rectangles identical
    /// four corner container   - first contains second
    /// four corner containee   - second contains first
    /// two corner container    - two corners of the second are contained by the first
    /// two corner containee    - two corners of the first are contained by the second
    /// one corner              - one corner of each rectangle is contained by the other
    string region_status(const rectangular_region & rec_region) const;

    ////////////////////////////
    ///  Function to sort corners
    ///  so that the index on the 
    ///  corners runs clockwise or
    ///  counterclockwise around the
    ///  rectangle
    void sort_corners();

    public:

    ////////////////////////////
    ///  Null constructor
    rectangular_region(){};
    
    ////////////////////////////
    ///  Copy constructor
    rectangular_region(const rectangular_region & rec_region);
    
    ////////////////////////////
    ///  Construct from four three_points,
    ///  which must be coplanar and form
    ///  a rectangular region
    rectangular_region(const vector<three_point> & in_corners);
    
    ////////////////////////////
    ///  Construct from a three_frame,
    ///  a set of axes, and a pixel scale
    rectangular_region(const three_frame & tf, vector<long> axes, double pix);
    
    ////////////////////////////
    ///  Construct from another rectangular region
    ///  by projecting it onto a plane with normal
    ///  vector nrml.  The normal vector must be
    ///  orthogonal to one of the two transverse
    ///  vectors defined by the sides of the rectangle
    ///  in order for the projected area to be
    ///  itself rectangular.  If this is not the
    ///  case, an error is thrown.
    ///
    ///  The bool along_nrml determines whether the 
    ///  projection is carried out along the three_vector 
    ///  nrml, or along the vector normal to the 
    ///  rectangular region.
    ///
    ///  If the projection is carried out along the three_vector
    ///  nrml, the resulting rectangular region is smaller than 
    ///  the original.  This method is used to find the rectangular
    ///  region that covers a foreshortened aperture.
    ///
    ///  If the projection is carried out along the normal 
    ///  to the rectangular region, the resulting rectangular 
    ///  region is larger than the original.  This method is used
    ///  to find the rectangular region of a wavefront projected
    ///  onto an optic.
    ///
    ///  In both cases, the new rectangular region has 
    ///  the same center as the original.
    ///
    ///  Constructing a rectangular region with along_nrml==false
    ///  and then using that region to construct another with 
    ///  along_nrml==true yields the original region.
    rectangular_region(const rectangular_region & rec_region,
    			const three_vector & nrml, bool along_nrml);
    
    ////////////////////////////
    ///  Construct from another rectangular_region
    ///  so that the resulting rectangular_region
    ///  has dimensions evenly divisible by pix.
    ///  This constructor rounds the dimensions of 
    ///  the rectangular_region up to get the smallest
    ///  possible rectangular region.  The center of
    ///  the original rectangular region is preserved
    rectangular_region(const rectangular_region & rec_region,
		       double pix);
    
    ////////////////////////////
    ///  Destructor
    ~rectangular_region(){};
    
    ////////////////////////////
    ///  Operator = 
    rectangular_region & operator=(const rectangular_region & rec_region);

    ////////////////////////////
    ///  Determine whether this
    ///  rectangular region is aligned
    ///  with rec_region.
    ///  Aligned means that the axes
    ///  defined by the edges of the
    ///  rectangles must match.
    bool aligned(const rectangular_region & rec_region) const;

    ////////////////////////////
    ///  Determine whether this
    ///  region contains rec_region.
    ///  This rectangular_region and rec_region must
    ///  be aligned 
    bool contains(const rectangular_region & rec_region) const;

    ////////////////////////////
    ///  Determine whether this
    ///  region is contained by rec_region.
    ///  This rectangular_region and rec_region must
    ///  be aligned 
    bool is_contained(const rectangular_region & rec_region) const;

    ////////////////////////////
    ///  Determine whether this
    ///  region has no overlap with rec_region.
    ///  This rectangular_region and rec_region must
    ///  be aligned 
    bool is_disjoint(const rectangular_region & rec_region) const;

    ////////////////////////////
    ///  Returns the three_points that 
    ///  identify this region
    vector<three_point> get_corners() const;

    ////////////////////////////
    ///  Returns the center point
    ///  of this region
    three_point get_center() const;    

    ////////////////////////////
    ///  Print
    void print(ostream & os, const char * prefix = "", long precision = 6) const;

    ////////////////////////////
    ///  Print
    void print(ostream & os, const three_frame & tf,
    		const char * prefix = "", long precision = 6) const;

    ////////////////////////////
    ///  Return the rectangular_region
    ///  of smallest size that contains
    ///  both rectangular regions supplied
    ///  as arguments.
    ///  rec_region1 and rec_region2 must
    ///  be aligned 
    friend rectangular_region region_union(const rectangular_region & rec_region1,
					   const rectangular_region & rec_region2);

    ////////////////////////////
    ///  Return the rectangular_region
    ///  of largest size that is contained
    ///  by both of the rectangular regions 
    ///  supplied as arguments.
    ///  rec_region1 and rec_region2 must
    ///  be aligned 
    friend rectangular_region region_intersection(
    			const rectangular_region & rec_region1,
			const rectangular_region & rec_region2);

    ////////////////////////////
    ///  Operator == 
    friend bool operator==(const rectangular_region & rec_region1,
			   const rectangular_region & rec_region2);

  };

  bool operator!=(const rectangular_region & rec_region1,
		  const rectangular_region & rec_region2);

}

#endif
