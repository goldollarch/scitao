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

#include <iostream>
#include <string>
#include "fits_factory.h"
#include "iofits.h"
#include "fits_header_data.h"
#include "pixel_array.h"
#include "region_base.h"
#include "three_vector.h"
#include "AO_algo.h"
#include "refractive_atmospheric_layer.h"
#include "diffractive_wavefront.h"

using namespace std;

namespace Arroyo {
  namespace {

    template<class T>
    AO_sim_base * create_refractive_atmospheric_layer(const iofits & iof) {
      return new refractive_atmospheric_layer<T>(iof);
    }

    const fits_keyval_set & get_refractive_atmospheric_layer_float_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "refractive atmospheric layer"));
      fkvs->push_back(fits_keyval_entry("BITPIX", "-32"));
      return *fkvs;
    }
    
    const fits_keyval_set & get_refractive_atmospheric_layer_double_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "refractive atmospheric layer"));
      fkvs->push_back(fits_keyval_entry("BITPIX", "-64"));
      return *fkvs;
    }
    
    template<>
    const bool refractive_atmospheric_layer<float>::factory_registration = 
    fits_factory<AO_sim_base>::Register(get_refractive_atmospheric_layer_float_keyval_set(), 
					create_refractive_atmospheric_layer<float>);

    template<>
    const bool refractive_atmospheric_layer<double>::factory_registration = 
    fits_factory<AO_sim_base>::Register(get_refractive_atmospheric_layer_double_keyval_set(), 
					create_refractive_atmospheric_layer<double>);

  } 

  /*
  void refractive_atmospheric_layer::transform(geometric_wavefront & gwf) const {
    long nrays = gwf.get_number_of_rays();
    geometric_ray * g_ray = this->get_wavefront_data(gwf);
    int x_index, y_index;
    for(int i=0; i<nrays; i++){
      if(g_ray[i].is_lost()) return;
      x_index = (int)(pixel_array<double>::axes[0]/2.0 + 
		      g_ray[i].three_point::x(*this)/pixel_scale);
      y_index = (int)(pixel_array<double>::axes[1]/2.0 + 
		      g_ray[i].three_point::y(*this)/pixel_scale);
      // we need to put a factor of lambda in here somewhere
      g_ray[i].set_pathlength(pixel_array<double>::pixeldata[x_index*axes[0]+y_index] + 
			      g_ray[i].get_pathlength());
    }
  }
  */

}
