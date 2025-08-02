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

#include <math.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include "wavefront_phase_estimate.h"
#include "fits_factory.h"

using namespace std;

namespace Arroyo {

  namespace {

    template<typename precision, typename aperture_type>
    AO_sim_base * create_wavefront_phase_estimate(const iofits & iof){
      return new wavefront_phase_estimate<precision, aperture_type>(iof);
    }

    const fits_keyval_set & get_wavefront_phase_estimate_float_circular_aperture_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "wavefront phase estimate"));
      fkvs->push_back(fits_keyval_entry("BITPIX", "-32"));
      fkvs->push_back(fits_keyval_entry("TYPE", "circular aperture", 1));
      return *fkvs;
    }
    
    const fits_keyval_set & get_wavefront_phase_estimate_double_circular_aperture_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "wavefront phase estimate"));
      fkvs->push_back(fits_keyval_entry("BITPIX", "-64"));
      fkvs->push_back(fits_keyval_entry("TYPE", "circular aperture", 1));
      return *fkvs;
    }

    template<>
    const bool wavefront_phase_estimate<float, circular_aperture>::factory_registration = 
    fits_factory<AO_sim_base>::Register(get_wavefront_phase_estimate_float_circular_aperture_keyval_set(), 
					create_wavefront_phase_estimate<float, circular_aperture>);
    template<>
    const bool wavefront_phase_estimate<double, circular_aperture>::factory_registration = 
    fits_factory<AO_sim_base>::Register(get_wavefront_phase_estimate_double_circular_aperture_keyval_set(), 
					create_wavefront_phase_estimate<double, circular_aperture>);
  } 

}
