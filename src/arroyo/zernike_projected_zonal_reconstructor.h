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

#ifndef ZERNIKE_PROJECTED_ZONAL_RECONSTRUCTOR_H
#define ZERNIKE_PROJECTED_ZONAL_RECONSTRUCTOR_H

#include "AO_sim_base.h"
#include "pixel_array.h"
#include <map>
#include <vector>
#include "computational_geometry.h"
#include "zernike.h"
#include "Shack_Hartmann_centroids.h"
#include "deformable_mirror.h"

namespace Arroyo {

  ///
  /// A virtual base class to hold single conjugate AO reconstructors.
  /// This class aims to represent reconstructors that yield zonal
  /// residuals as well as some number of zernike mode residuals.
  ///

  class zernike_projected_zonal_reconstructor :
    public AO_sim_base {

    public:

    ///////////////////////////////////////////
    ///  Null constructor
    zernike_projected_zonal_reconstructor(){};

    ///////////////////////////////////////////
    ///  Virtual destructor
    virtual ~zernike_projected_zonal_reconstructor(){};
    
    ///////////////////////////////////////////
    ///  Read from file
    virtual void read(const char * filename) = 0;

    ///////////////////////////////////////////
    ///  Read from iofits
    virtual void read(const Arroyo::iofits & iof) = 0;

    ///////////////////////////////////////////
    ///  Write to file
    virtual void write(const char * filename) const = 0;
 
    ///////////////////////////////////////////
    ///  Write to iofits
    virtual void write(Arroyo::iofits & iof) const = 0;
 
    ///////////////////////////////////////////
    ///  Print
    virtual void print(std::ostream & os, const char * prefix="") const = 0;

    ///////////////////////////////////////////
    ///  Get dimensions of centroid measurements passed
    ///  to the reconstructor
    virtual vector<long> get_centroid_axes() const = 0;

    ///////////////////////////////////////////
    ///  Get dimensions of pixel array returned
    ///  by the reconstructor
    virtual vector<long> get_actuator_axes() const = 0;

    ///////////////////////////////////////////
    ///  Get a zernike instance that contains information
    ///  about which modes are returned by the reconstructor.
    /// 
    ///  This instance is minimally sized so as to hold the largest
    ///  mode returned by the reconstructor.  Each element of this
    ///  instance is initialized to unity if the corresponding mode is
    ///  returned by the reconstructor, and to zero if it is not
    virtual Arroyo::zernike get_zernike_modes() const = 0;

    ///////////////////////////////////////////
    ///  Reconstruct the zernike residuals from a 
    ///  Shack Hartmann centroid class instance
    virtual void reconstruct_zernike_residuals(const Arroyo::Shack_Hartmann_centroids & shcentroids, 
					       Arroyo::zernike & znke) const = 0;
    
    ///////////////////////////////////////////
    ///  Reconstruct the zonal residuals from a
    ///  Shack Hartmann centroid class instance
    virtual void reconstruct_zonal_residuals(const Arroyo::Shack_Hartmann_centroids & shcentroids, 
					     Arroyo::pixel_array<double> & pixarr) const = 0;

    ///////////////////////////////////////////
    ///  Reconstruct the zernike and zonal residuals from a
    ///  Shack Hartmann centroid class instance
    virtual void reconstruct_residuals(const Arroyo::Shack_Hartmann_centroids & shcentroids, 
				       Arroyo::zernike & znke, 
				       Arroyo::pixel_array<double> & pixarr) const = 0;

    ///////////////////////////////////////////
    ///  Factory to construct reconstructors from file
    static zernike_projected_zonal_reconstructor * 
      zernike_projected_zonal_reconstructor_factory(const char * filename);

    ///////////////////////////////////////////
    ///  Factory to construct reconstructors from an iofits object
    static zernike_projected_zonal_reconstructor * 
      zernike_projected_zonal_reconstructor_factory(const iofits & iof);

  };

}

#endif
