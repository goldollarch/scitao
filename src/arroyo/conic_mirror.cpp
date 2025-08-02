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
#include <vector>
#include "fits_factory.h"
#include "region_base.h"
#include "aperture.h"
#include "conic_mirror.h"

using namespace std;

namespace Arroyo {

  namespace {

    template<class aperture_type>
    AO_sim_base * create_conic_mirror(const iofits & iof) {
      return new conic_mirror<aperture_type>(iof);
    }

    const fits_keyval_set & get_conic_mirror_circular_aperture_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "conic mirror"));
      fkvs->push_back(fits_keyval_entry("TYPE", "circular aperture", 1));
      return *fkvs;
    }
    const fits_keyval_set & get_conic_mirror_annular_aperture_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "conic mirror"));
      fkvs->push_back(fits_keyval_entry("TYPE", "annular aperture", 1));
      return *fkvs;
    }
    const fits_keyval_set & get_conic_mirror_spidered_annular_aperture_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "conic mirror"));
      fkvs->push_back(fits_keyval_entry("TYPE", "spidered annular aperture", 1));
      return *fkvs;
    }
    const fits_keyval_set & get_conic_mirror_rectangular_aperture_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "conic mirror"));
      fkvs->push_back(fits_keyval_entry("TYPE", "rectangular aperture", 1));
      return *fkvs;
    }
    const fits_keyval_set & get_conic_mirror_hexagonal_aperture_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "conic mirror"));
      fkvs->push_back(fits_keyval_entry("TYPE", "hexagonal aperture", 1));
      return *fkvs;
    }

    template<>
    const bool conic_mirror<circular_aperture>::factory_registration = 
    fits_factory<AO_sim_base>::Register(get_conic_mirror_circular_aperture_keyval_set(), 
					create_conic_mirror<circular_aperture>);
    template<>
    const bool conic_mirror<annular_aperture>::factory_registration = 
    fits_factory<AO_sim_base>::Register(get_conic_mirror_annular_aperture_keyval_set(), 
					create_conic_mirror<annular_aperture>);
    template<>
    const bool conic_mirror<spidered_annular_aperture>::factory_registration = 
    fits_factory<AO_sim_base>::Register(get_conic_mirror_spidered_annular_aperture_keyval_set(), 
					create_conic_mirror<spidered_annular_aperture>);
    template<>
    const bool conic_mirror<rectangular_aperture>::factory_registration = 
    fits_factory<AO_sim_base>::Register(get_conic_mirror_rectangular_aperture_keyval_set(), 
					create_conic_mirror<rectangular_aperture>);
    template<>
    const bool conic_mirror<hexagonal_aperture>::factory_registration = 
    fits_factory<AO_sim_base>::Register(get_conic_mirror_hexagonal_aperture_keyval_set(), 
					create_conic_mirror<hexagonal_aperture>);
  }
}
