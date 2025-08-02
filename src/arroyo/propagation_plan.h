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

#ifndef PROPAGATION_PLAN_H
#define PROPAGATION_PLAN_H

#include "wavefront_header.h"

namespace Arroyo {


  ///
  /// A virtual base class for different types of propagation plans
  ///

  class propagation_plan {

  public:

    ///////////////////////////////////////////
    ///  Null constructor
    propagation_plan(){};

    ///////////////////////////////////////////
    ///  Virtual destructor
    virtual ~propagation_plan(){};

    ///////////////////////////////////////////
    ///  Virtual member function to return
    ///  a diffractive wavefront header that
    ///  has been padded enough that subsequent
    ///  propagation through a distance dist will 
    ///  not corrupt the data within a subarray 
    ///  equal to the original axes of the diffractive
    ///  wavefront header dwf
    ///
    ///  This member function is defined explicitly for
    ///  the template argument to diffractive_wavefront_header
    ///  because virtual template member functions are not
    ///  supported in C++
    virtual diffractive_wavefront_header<float> 
      pad(const diffractive_wavefront_header<float> & dwf, double dist) const = 0;

    ///////////////////////////////////////////
    ///  Virtual member function to return
    ///  a diffractive wavefront header that
    ///  has been padded enough that subsequent
    ///  propagation through a distance dist will 
    ///  not corrupt the data within a subarray 
    ///  equal to the original axes of the diffractive
    ///  wavefront header dwf
    ///
    ///  This member function is defined explicitly for
    ///  the template argument to diffractive_wavefront_header
    ///  because virtual template member functions are not
    ///  supported in C++
    virtual diffractive_wavefront_header<double> 
      pad(const diffractive_wavefront_header<double> & dwf, double dist) const = 0;

  };

  ///
  /// A class to represent a geometric propagation plan
  ///

  class geometric_propagation_plan :
    public propagation_plan {

  public:
    
    ///////////////////////////////////////////
    ///  Null constructor
    geometric_propagation_plan(){};

    ///////////////////////////////////////////
    ///  Destructor
    ~geometric_propagation_plan(){};

    ///////////////////////////////////////////
    ///  Return an exact copy of the original
    ///  wavefront header, since geometric 
    ///  propagation requires no padding.
    diffractive_wavefront_header<float> 
      pad(const diffractive_wavefront_header<float> & dwf, double dist) const {
      return (dwf);
    };

    ///////////////////////////////////////////
    ///  Return an exact copy of the original
    ///  wavefront header, since geometric 
    ///  propagation requires no padding.
    diffractive_wavefront_header<double> 
      pad(const diffractive_wavefront_header<double> & dwf, double dist) const {
      return (dwf);
    };

  };


  ///
  /// A class to represent a near field fresnel propagation plan
  ///

  class near_field_fresnel_propagation_plan :
    public propagation_plan {

  private:

    ///////////////////////////////////////////
    ///  A private template member function to
    ///  cover the explicit public versions
    template<class T>
      diffractive_wavefront_header<T> private_pad(const diffractive_wavefront_header<T> & dwf, double dist) const;
    
  public:

    ///////////////////////////////////////////
    ///  Null constructor
    near_field_fresnel_propagation_plan(){};

    ///////////////////////////////////////////
    ///  Destructor
    ~near_field_fresnel_propagation_plan(){};

    ///////////////////////////////////////////
    ///  Return a diffractive wavefront header
    ///  padded by an amount equal to 
    ///  wavelength*dist/pixscale/pixscale
    diffractive_wavefront_header<float> 
      pad(const diffractive_wavefront_header<float> & dwf, double dist) const;

    ///////////////////////////////////////////
    ///  Return a diffractive wavefront header
    ///  padded by an amount equal to 
    ///  wavelength*dist/pixscale/pixscale
    diffractive_wavefront_header<double> 
      pad(const diffractive_wavefront_header<double> & dwf, double dist) const;

  };

}

#endif
