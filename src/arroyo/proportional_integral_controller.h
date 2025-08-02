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

#ifndef PROPORTIONAL_INTEGRAL_CONTROLLER_H
#define PROPORTIONAL_INTEGRAL_CONTROLLER_H

#include <vector>
#include <ostream>
#include "iofits.h"

namespace Arroyo {

  using std::vector;
  using std::ostream;

  ///
  /// A template class to represent a proportional integral
  /// controller.
  ///
  /// Note that the units on the proportional gain and integral
  /// gain differ.
  ///
  /// This class requires that the following
  /// operators be defined
  ///
  ///   input::operator*=(gain)
  ///   input::operator-=(input)
  ///   output::operator+=(input)
  ///   input::operator=(input)
  ///   output::operator=(output)
  ///   gain::operator=(gain)
  ///
  /// Illustrative examples:
  ///
  /// Simplest possible controllers
  ///   proportional_integral_controller(double, double, double, double)
  ///   proportional_integral_controller(float, float, float, float)
  ///
  /// SCAO zonal control, single gain
  ///   proportional_integral_controller(pixel_array<T>, pixel_array<U>, double, double)
  ///
  /// SCAO zonal control, independent gains
  ///   proportional_integral_controller(pixel_array<T>, pixel_array<U>, pixel_array<V>, pixel_array<W>)
  ///
  /// SCAO modal control, single gain
  ///   proportional_integral_controller(zernike, zernike, double, double)
  ///
  /// SCAO modal control, independent gains
  ///   proportional_integral_controller(zernike, zernike, pixel_array<T>, pixel_array<T>)
  ///
  /// atmospheric reconstruction control, single gain
  ///   proportional_integral_controller(refractive_atmospheric_layer, refractive_atmospheric_layer, double, double)


  template<class input, class output, class proportional_gain, class integral_gain>
    class proportional_integral_controller {

    private:

    static bool factory_registration;

    protected:

    // proportional gain
    proportional_gain propgain;

    // last input for proportional update
    input last_input;

    // integral gain 
    integral_gain intgain;

    ///////////////////////////////////////////
    ///  Null constructor
    proportional_integral_controller(){};

    public:
    
    ///////////////////////////////////////////
    ///  Copy constructor
    proportional_integral_controller(const proportional_integral_controller & picntrlr);

    ///////////////////////////////////////////
    ///  Construct from file
    proportional_integral_controller(const char * filename);

    ///////////////////////////////////////////
    ///  Construct from iofits object
    proportional_integral_controller(const iofits & iof);

    ///////////////////////////////////////////
    ///  Construct from the bits
    proportional_integral_controller(const input & init,
				     const proportional_gain & pgain, 
				     const integral_gain & igain);

    ///////////////////////////////////////////
    ///  Destructor
    ~proportional_integral_controller(){};

    ///////////////////////////////////////////
    ///  Operator = 
    proportional_integral_controller & operator=(const proportional_integral_controller & picntrlr);

    ///////////////////////////////////////////
    ///  Generate commands from residuals
    void update(const input & i, output & o);
    
  };


  template<class input, class output, class proportional_gain, class integral_gain> 
    proportional_integral_controller<input, output, proportional_gain, integral_gain>::
    proportional_integral_controller(const proportional_integral_controller<input, output, proportional_gain, integral_gain> & picntrlr){
    this->operator=(picntrlr);
  }

  template<class input, class output, class proportional_gain, class integral_gain> 
    proportional_integral_controller<input, output, proportional_gain, integral_gain>::
    proportional_integral_controller(const char * filename){
    this->read(filename);
  }

  template<class input, class output, class proportional_gain, class integral_gain> 
    proportional_integral_controller<input, output, proportional_gain, integral_gain>::
    proportional_integral_controller(const iofits & iof){
    this->read(iof);
  }

  template<class input, class output, class proportional_gain, class integral_gain> 
    proportional_integral_controller<input, output, proportional_gain, integral_gain>::
    proportional_integral_controller(const input & init, 
				     const proportional_gain & propg, 
				     const integral_gain & intg){
    last_input = init;
    propgain = propg;
    intgain = intg;
  }

  template<class input, class output, class proportional_gain, class integral_gain> 
    proportional_integral_controller<input, output, proportional_gain, integral_gain> & 
    proportional_integral_controller<input, output, proportional_gain, integral_gain>::
    operator=(const proportional_integral_controller<input, output, proportional_gain, integral_gain> & picntrlr){
    if(this==&picntrlr)
      return(*this);

    last_input = picntrlr.last_input;
    propgain = picntrlr.propgain;
    intgain = picntrlr.intgain;

    return(*this);
  }

  template<class input, class output, class proportional_gain, class integral_gain> 
    void proportional_integral_controller<input, output, proportional_gain, integral_gain>::
    update(const input & i, output & o){
    o += i*intgain + (i - last_input)*propgain;
    last_input = i;
  }

}

#endif
