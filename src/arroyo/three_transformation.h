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

#ifndef THREE_TRANSFORMATION_H
#define THREE_TRANSFORMATION_H

#include "AO_cpp.h"
#include "iofits.h"

namespace Arroyo {

  using std::ostream;

  class three_point;
  class three_vector;
  class three_frame;

  ///
  /// A class to represent a transformation in three dimensional space.
  ///

  class three_transformation {

  protected:

    /// Nine components defined relative to a global coordinate system
    double xx_, xy_, xz_, dx_;
    double yx_, yy_, yz_, dy_;
    double zx_, zy_, zz_, dz_;

    ///////////////////////////////////////////
    ///  Protected constructor
    three_transformation(double xx, double xy, double xz, double dx,
			 double yx, double yy, double yz, double dy,
			 double zx, double zy, double zz, double dz){
      xx_ = xx; xy_ = xy; xz_ = xz; dx_ = dx;
      yx_ = yx; yy_ = yy; yz_ = yz; dy_ = dy;
      zx_ = zx; zy_ = zy; zz_ = zz; dz_ = dz;
    };

    ///////////////////////////////////////////
    ///  Protected modifier
    void set_three_transformation(double xx, double xy, double xz, double dx,
				  double yx, double yy, double yz, double dy,
				  double zx, double zy, double zz, double dz){
      xx_ = xx; xy_ = xy; xz_ = xz; dx_ = dx;
      yx_ = yx; yy_ = yy; yz_ = yz; dy_ = dy;
      zx_ = zx; zy_ = zy; zz_ = zz; dz_ = dz;
    };

  public:
  
    ///////////////////////////////////////////
    ///  Null constructor
    three_transformation();

    ///////////////////////////////////////////
    ///  Copy constructor
    three_transformation(const three_transformation & tt);

    ///////////////////////////////////////////
    ///  Destructor
    ~three_transformation(){};

    ///////////////////////////////////////////
    ///  Operator = 
    three_transformation & operator=(const three_transformation & tt);

    ///////////////////////////////////////////
    ///  Reset transformation to tt
    void set_transformation(const three_transformation & tt) {
      xx_ = tt.xx_; xy_ = tt.xy_; xz_ = tt.xz_; dx_ = tt.dx_;
      yx_ = tt.yx_; yy_ = tt.yy_; yz_ = tt.yz_; dy_ = tt.dy_;
      zx_ = tt.zx_; zy_ = tt.zy_; zz_ = tt.zz_; dz_ = tt.dz_;
    };

    ///////////////////////////////////////////
    ///  Read from iofits object
    void read(const iofits & iof);

    ///////////////////////////////////////////
    ///  Write to iofits object
    void write(iofits & iof) const;

    ///////////////////////////////////////////
    ///  Write to iofits object
    void print(ostream & os, const char * prefix = "") const;

    ///////////////////////////////////////////
    ///  Apply this transformation to a point
    void transform(three_point & tp) const;

    ///////////////////////////////////////////
    ///  Apply this transformation to a vector
    void transform(three_vector & tv) const;

    ///////////////////////////////////////////
    ///  Get the inverse transformation
    three_transformation inverse() const;

    ///////////////////////////////////////////
    ///  Multiply this transformation with another.
    ///  The application of these transformations
    ///  to a point, vector, or frame is defined
    ///  in the sense that first this transformation
    ///  is applied, and then the transformation tt
    ///  is applied.  In matrix form, this means
    ///  that this function performs a left matrix
    ///  multiply of this transformation by the transformation
    ///  tt.
    three_transformation & operator*=(const three_transformation & tt);

    ///////////////////////////////////////////
    ///  Operator==
    friend bool operator==(const three_transformation & tt1,
    				const three_transformation & tt2);

    ///////////////////////////////////////////
    ///  Multiply two transformations to yield a 3rd.
    ///  The application of these transformations
    ///  to a point, vector, or frame is defined
    ///  in the sense that first the transformation tt2
    ///  is applied, and then the transformation tt1
    ///  is applied.  In matrix form, this means
    ///  that this function performs a left matrix
    ///  multiply of the transformation tt2 by the 
    ///  transformation tt1.  That is, the arguments
    ///  are read right to left.
    friend three_transformation operator*(const three_transformation & tt1,
    				const three_transformation & tt2);

    /// A verbose_level for printing messages
    static int verbose_level;

  };

  ///////////////////////////////////////////
  ///  Operator!=
  bool operator!=(const three_transformation & tv1,
  				const three_transformation & tv2);

  three_transformation operator*(const three_transformation & tt1,
  				const three_transformation & tt2);

  ///
  /// A class to represent an orthonormal transformation in three dimensional space.
  ///

  class three_orthonormal_transformation :
    public three_transformation {

    public:

    ///////////////////////////////////////////
    ///  Null constructor
    three_orthonormal_transformation(){};

    ///////////////////////////////////////////
    ///  Copy constructor
    three_orthonormal_transformation(const three_orthonormal_transformation & tt) :
      three_transformation(tt){};

    ///////////////////////////////////////////
    ///  Destructor
    ~three_orthonormal_transformation(){};

    ///////////////////////////////////////////
    ///  Get the inverse transformation
    three_orthonormal_transformation inverse() const;

    ///////////////////////////////////////////
    ///  Apply this transformation to a three point
    ///
    ///  I shouldn't have to write this function, 
    ///  but for some reason the compiler won't skip
    ///  from the classes derived from this one 
    ///  back to the three_transformation class.
    void transform(three_point & tp) const {
    		this->three_transformation::transform(tp);};

    ///////////////////////////////////////////
    ///  Apply this transformation to a three vector
    ///
    ///  I shouldn't have to write this function, 
    ///  but for some reason the compiler won't skip
    ///  from the classes derived from this one 
    ///  back to the three_transformation class.
    void transform(three_vector & tv) const {
    			this->three_transformation::transform(tv);};

    ///////////////////////////////////////////
    ///  Apply this transformation to a three frame
    void transform(three_frame & tf) const;

  };


  /// A Class to represent a translation in three dimensional space.
  ///
  /// Affects a three_point.
  /// Has no effect on a three_vector.
  /// Affects a three_frame.

  class three_translation :
    public three_orthonormal_transformation {

    public:
  
    ///////////////////////////////////////////
    ///  Null constructor
    three_translation(){};

    ///////////////////////////////////////////
    ///  Copy constructor
    three_translation(const three_translation & tt) :
      three_orthonormal_transformation(tt){};

    ///////////////////////////////////////////
    ///  Construct a translation by a vector v
    three_translation(const three_vector & v);

    ///////////////////////////////////////////
    ///  Destructor
    ~three_translation(){};

    ///////////////////////////////////////////
    ///  Get the inverse translation
    three_translation inverse() const;

  };
 

  /// A class to represent a rotation in three dimensional space.
  ///
  /// Affects a three_point.
  /// Affects a three_vector.
  /// Affects a three_frame.

  class three_rotation :
    public three_orthonormal_transformation {

    public:
  
    ///////////////////////////////////////////
    ///  Null constructor
    three_rotation(){};

    ///////////////////////////////////////////
    ///  Copy constructor
    three_rotation(const three_rotation & tr) :
      three_orthonormal_transformation(tr){};

    ///////////////////////////////////////////
    ///  Construct a rotation about an axis 
    ///  specified by a three_vector, by an angle.
    ///  The angle is specified in radians.
    ///
    ///  This function throws an error if tv is
    ///  a null vector
    three_rotation(const three_point & tp, const three_vector & tv, double angle);

    ///////////////////////////////////////////
    ///  Destructor
    ~three_rotation(){};

    ///////////////////////////////////////////
    ///  Get the inverse rotation
    three_rotation inverse() const;

  };

  /// A class to represent a reflection in three dimensional space.
  ///
  /// Affects a three_point.
  /// Affects a three_vector.
  /// Affects a three_frame.

  class three_reflection :
    public three_orthonormal_transformation {

    public:
  
    ///////////////////////////////////////////
    ///  Null constructor
    three_reflection(){};

    ///////////////////////////////////////////
    ///  Copy constructor
    three_reflection(const three_reflection & tr) :
      three_orthonormal_transformation(tr){};

    ///////////////////////////////////////////
    ///  Construct a reflection about a plane
    ///  specified by a three_vector and a point
    ///  
    ///  This function throws an error if tv is a
    ///  null vector
    three_reflection(const three_point & tp, const three_vector & tv);

    ///////////////////////////////////////////
    ///  Destructor
    ~three_reflection(){};

    ///////////////////////////////////////////
    ///  Get the inverse reflection
    three_reflection inverse() const;

  };

  /// A class to represent a scaling in three dimensional space.
  ///
  /// Affects a three_point.
  /// Affects a three_vector.
  /// Cannot be applied to a three_frame.

  class three_scaling :
    public three_transformation {

    public:
  
    ///////////////////////////////////////////
    ///  Null constructor
    three_scaling(){};

    ///////////////////////////////////////////
    ///  Copy constructor
    three_scaling(const three_scaling & ts) :
      three_transformation(ts){};

    ///////////////////////////////////////////
    ///  Construct a uniform scale transformation 
    ///
    ///  This function throws an error if scale <= 0
    three_scaling(const three_point & tp, double scale);

    ///////////////////////////////////////////
    ///  Construct a scale transformation along 
    ///  an axis specified by a three vector
    ///  by a scalar factor scale.  
    ///
    ///  This function throws an error if
    ///  tv is a null vector, or if scale <= 0
    three_scaling(const three_point & tp, const three_vector & tv, double scale);

    ///////////////////////////////////////////
    ///  Destructor
    ~three_scaling(){};

    ///////////////////////////////////////////
    ///  Get the inverse scaling
    three_scaling inverse() const;

  };

}
#endif
