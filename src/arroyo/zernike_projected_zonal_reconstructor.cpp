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
#include "fits_factory.h"
#include "region_base.h"
#include "zernike_projected_zonal_reconstructor.h"

using namespace std;

namespace Arroyo {

  zernike_projected_zonal_reconstructor * 
  zernike_projected_zonal_reconstructor::zernike_projected_zonal_reconstructor_factory(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...) {
      cerr << "zernike_projected_zonal_reconstructor::zernike_projected_zonal_reconstructor_factory - "
	   << "error opening file " << filename << endl;
      throw(string("zernike_projected_zonal_reconstructor::zernike_projected_zonal_reconstructor_factory"));
    }
    return(zernike_projected_zonal_reconstructor::zernike_projected_zonal_reconstructor_factory(iof));
  }

  zernike_projected_zonal_reconstructor * 
  zernike_projected_zonal_reconstructor::zernike_projected_zonal_reconstructor_factory(const iofits & iof){
    AO_sim_base * aosb = AO_sim_base::AO_sim_base_factory(iof);
    zernike_projected_zonal_reconstructor * zpzr = dynamic_cast<zernike_projected_zonal_reconstructor *>(aosb);
    if(zpzr==NULL)
      throw(string("zernike_projected_zonal_reconstructor::zernike_projected_zonal_reconstructor_factory"));
     return(zpzr);    
  }

} 
