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
#include "deformable_mirror.h"
#include "aperture.h"

using namespace std;

namespace Arroyo 
{
  
  /////////////////////////////////////////////////////////////////////
  // deformable_mirror_base class function declarations.
  /////////////////////////////////////////////////////////////////////
  deformable_mirror_base* deformable_mirror_base::deformable_mirror_base_factory(const char* filename)
  {
    iofits iof;
    try
    {
      iof.open(filename);
    }
    catch(...)
    {
      cerr << "deformable_mirror_base::deformable_mirror_base_factory - "
           << "error opening file: " << filename << endl;
      throw(string("deformable_mirror_base::deformable_mirror_base_factory"));
    }
    return (deformable_mirror_base::deformable_mirror_base_factory(iof));
  }

  deformable_mirror_base* deformable_mirror_base::deformable_mirror_base_factory(const iofits& iof)
  {
    AO_sim_base* aosb = AO_sim_base::AO_sim_base_factory(iof);
    deformable_mirror_base* dmbase = dynamic_cast<deformable_mirror_base *>(aosb);

    if (dmbase==NULL) throw(string("deformable_mirror_base::deformable_mirror_base_factory"));

    return(dmbase);
  }

  /////////////////////////////////////////////////////////////////////
  // ideal_deformable_mirror class function declarations.
  /////////////////////////////////////////////////////////////////////

  ///////////////////////////////////
  // First do the class registration:
  // Use an anonymous namespace....
  ///////////////////////////////////
  namespace 
  {
    ////////////////////////////////////////////////////////
    // 1. Define the function to be executed at registration time
    //    if many registrations are required, this is the case
    //    when templated on objects ocurred.
    ////////////////////////////////////////////////////////
    template<class aperture_type> 
    AO_sim_base* create_ideal_deformable_mirror(const iofits& iof)
    {
      return new ideal_deformable_mirror<aperture_type>(iof);
    }

    ////////////////////////////////////////////////////////
    // 2. Set the TYPE of the object and templated objects.
    ////////////////////////////////////////////////////////
    const fits_keyval_set & 
    get_ideal_deformable_mirror_circular_aperture_keyval_set()
    {
      static fits_keyval_set *fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE","ideal_deformable_mirror"));
      fkvs->push_back(fits_keyval_entry("TYPE","circular aperture",2));
      return *fkvs;
    }

    const fits_keyval_set & 
    get_ideal_deformable_mirror_annular_aperture_keyval_set()
    {
      static fits_keyval_set *fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE","ideal_deformable_mirror"));
      fkvs->push_back(fits_keyval_entry("TYPE","annular aperture",2));
      return *fkvs;
    }

    const fits_keyval_set & 
    get_ideal_deformable_mirror_spidered_annular_aperture_keyval_set()
    {
      static fits_keyval_set *fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE","ideal_deformable_mirror"));
      fkvs->push_back(fits_keyval_entry("TYPE","spidered annular aperture",2));
      return *fkvs;
    }

    const fits_keyval_set & 
    get_ideal_deformable_mirror_rectangular_aperture_keyval_set()
    {
      static fits_keyval_set *fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE","ideal_deformable_mirror"));
      fkvs->push_back(fits_keyval_entry("TYPE","rectangular aperture",2));
      return *fkvs;
    }

    const fits_keyval_set & 
    get_ideal_deformable_mirror_hexagonal_aperture_keyval_set()
    {
      static fits_keyval_set *fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE","ideal_deformable_mirror"));
      fkvs->push_back(fits_keyval_entry("TYPE","hexagonal aperture",2));
      return *fkvs;
    }


    ////////////////////////////////////////////////////////
    // 3. Finally register all the objects
    ////////////////////////////////////////////////////////
    template<>
    const bool ideal_deformable_mirror<circular_aperture>::factory_registration = 
    fits_factory<AO_sim_base>::Register(get_ideal_deformable_mirror_circular_aperture_keyval_set(),
					create_ideal_deformable_mirror<circular_aperture>);
    
    template<>
    const bool ideal_deformable_mirror<annular_aperture>::factory_registration = 
    fits_factory<AO_sim_base>::Register(get_ideal_deformable_mirror_annular_aperture_keyval_set(),
					create_ideal_deformable_mirror<annular_aperture>);

    template<>
    const bool ideal_deformable_mirror<spidered_annular_aperture>::factory_registration = 
    fits_factory<AO_sim_base>::Register(get_ideal_deformable_mirror_spidered_annular_aperture_keyval_set(),
					create_ideal_deformable_mirror<spidered_annular_aperture>);

    template<>
    const bool ideal_deformable_mirror<rectangular_aperture>::factory_registration = 
    fits_factory<AO_sim_base>::Register(get_ideal_deformable_mirror_rectangular_aperture_keyval_set(),
					create_ideal_deformable_mirror<rectangular_aperture>);
    
    template<>
    const bool ideal_deformable_mirror<hexagonal_aperture>::factory_registration = 
    fits_factory<AO_sim_base>::Register(get_ideal_deformable_mirror_hexagonal_aperture_keyval_set(),
					create_ideal_deformable_mirror<hexagonal_aperture>);
  }
}
