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

#include "fits_factory.h"
#include "arroyo_least_squares_reconstructor.h"

using namespace std;

namespace Arroyo {
  namespace {

    template<class T>
    AO_sim_base * create_arroyo_least_squares_reconstructor(const iofits & iof){
      return new arroyo_least_squares_reconstructor<T>(iof);
    }

    const fits_keyval_set & get_arroyo_least_squares_reconstructor_float_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "arroyo least squares reconstructor"));
      fkvs->push_back(fits_keyval_entry("BITPIX", "-32"));
      return *fkvs;
    }
    const fits_keyval_set & get_arroyo_least_squares_reconstructor_double_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "arroyo least squares reconstructor"));
      fkvs->push_back(fits_keyval_entry("BITPIX", "-64"));
      return *fkvs;
    }

    template<>
    const bool arroyo_least_squares_reconstructor<float>::factory_registration = 
    fits_factory<AO_sim_base>::Register(get_arroyo_least_squares_reconstructor_float_keyval_set(), 
					create_arroyo_least_squares_reconstructor<float>);
    template<>
    const bool arroyo_least_squares_reconstructor<double>::factory_registration = 
    fits_factory<AO_sim_base>::Register(get_arroyo_least_squares_reconstructor_double_keyval_set(), 
					create_arroyo_least_squares_reconstructor<double>);
  } 
}
