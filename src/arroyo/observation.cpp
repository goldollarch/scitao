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
#include "observation.h"

using namespace std;

namespace Arroyo 
{
  
  /////////////////////////////////////////////////////////////////////
  // observation_base class factories
  /////////////////////////////////////////////////////////////////////
  observation_base* observation_base::observation_base_factory(const char* filename)
  {
    iofits iof;
    try
    {
      iof.open(filename);
    }
    catch(...)
    {
      cerr << "observation_base::observation_base_factory - "
           << "error opening file: " 
	   << filename 
	   << endl;
      throw(string("observation_base::observation_base_factory"));
    }
    return (observation_base::observation_base_factory(iof));
  }

  observation_base* observation_base::observation_base_factory(const iofits & iof)
  {
    AO_sim_base* aosb = AO_sim_base::AO_sim_base_factory(iof);
    observation_base* obase = dynamic_cast<observation_base *>(aosb);

    if (obase==NULL) throw(string("observation_base::observation_base_factory"));

    return(obase);
  }

  namespace {
    template<class type>
    AO_sim_base * create_basic_observation(const iofits & iof){
      return new basic_observation<type>(iof);
    }

    const fits_keyval_set & get_basic_observation_float_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "basic observation"));
      fkvs->push_back(fits_keyval_entry("BITPIX", "-32"));
      return *fkvs;
    }
    const fits_keyval_set & get_basic_observation_double_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "basic observation"));
      fkvs->push_back(fits_keyval_entry("BITPIX", "-64"));
      return *fkvs;
    }

    template <>
    const bool basic_observation<float>::factory_registration = 
    fits_factory<AO_sim_base>::Register(get_basic_observation_float_keyval_set(), create_basic_observation<float>);
    template <>
    const bool basic_observation<double>::factory_registration = 
    fits_factory<AO_sim_base>::Register(get_basic_observation_double_keyval_set(), create_basic_observation<double>);
  } 


  namespace {
    template<class type>
    AO_sim_base * create_basic_otf(const iofits & iof){
      return new basic_otf<type>(iof);
    }

    const fits_keyval_set & get_basic_otf_float_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "basic otf"));
      fkvs->push_back(fits_keyval_entry("BITPIX", "-32"));
      return *fkvs;
    }
    const fits_keyval_set & get_basic_otf_double_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "basic otf"));
      fkvs->push_back(fits_keyval_entry("BITPIX", "-64"));
      return *fkvs;
    }

    template <>
    const bool basic_otf<float>::factory_registration = 
    fits_factory<AO_sim_base>::Register(get_basic_otf_float_keyval_set(), create_basic_otf<float>);
    template <>
    const bool basic_otf<double>::factory_registration = 
    fits_factory<AO_sim_base>::Register(get_basic_otf_double_keyval_set(), create_basic_otf<double>);
  } 
}

