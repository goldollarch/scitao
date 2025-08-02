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

#ifndef ARCADIA_MCAO_RECONSTRUCTOR_H
#define ARCADIA_MCAO_RECONSTRUCTOR_H

#include "AO_sim_base.h"
#include "pixel_array.h"
#include <vector>
#include "aperture.h"
#include "refractive_atmosphere.h"
#include "lenslet_array.h"
#include "tip_tilt_mirror.h"
#include "deformable_mirror.h"
#include "zernike.h"
#include "Shack_Hartmann_centroids.h"
 
namespace Arroyo {

  ///
  /// A class to hold an Arcadia MCAO reconstructor.
  /// 

  class arcadia_mcao_reconstructor :
    public AO_sim_base {

    private:

    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    std::string unique_name() const {return(std::string("Arcadia MCAO reconstructor"));};

    protected:

    ///////////////////////////////////////////
    ///  Null constructor
    arcadia_mcao_reconstructor(){};

    Arroyo::circular_aperture circ_ap;
    Arroyo::refractive_atmospheric_model * ref_atm_model;
    std::vector<Arroyo::emitter *> guide_stars;
    std::vector<bool> tilt_removal_flags;
    std::vector<Arroyo::square_lenslet_array> lenslet_arrays;
    std::vector<Arroyo::ideal_tip_tilt_mirror<circular_aperture> > ttms;
    std::vector<Arroyo::ideal_deformable_mirror<circular_aperture> > dms;
    std::vector<std::string> slaving_files;

    Arroyo::pixel_array<double> mcao_reconstructor;

    public:

    ///////////////////////////////////////////
    ///  Copy constructor
    arcadia_mcao_reconstructor(const arcadia_mcao_reconstructor & arcadia_mcao_recon);

    ///////////////////////////////////////////
    ///  Construct from fits file
    arcadia_mcao_reconstructor(const char * filename);

    ///////////////////////////////////////////
    ///  Construct from iofits object
    arcadia_mcao_reconstructor(const Arroyo::iofits & iof);

    ///////////////////////////////////////////
    ///  Construct from arcadia input file and reconstructor
    arcadia_mcao_reconstructor(const char * arcadia_input_file,
			       const char * arcadia_mcao_reconstructor_file);

    ///////////////////////////////////////////
    ///  Virtual destructor
    ~arcadia_mcao_reconstructor(){};

    ///////////////////////////////////////////
    ///  Operator = 
    arcadia_mcao_reconstructor & operator=(const arcadia_mcao_reconstructor & arcadia_mcao_recon);

    ///////////////////////////////////////////
    ///  Read from file
    void read(const char * filename);

    ///////////////////////////////////////////
    ///  Read from iofits
    void read(const Arroyo::iofits & iof);

    ///////////////////////////////////////////
    ///  Write to file
    void write(const char * filename) const;
 
    ///////////////////////////////////////////
    ///  Write to iofits
    void write(Arroyo::iofits & iof) const;
 
    ///////////////////////////////////////////
    ///  Print
    void print(std::ostream & os, const char * prefix="") const;

    ///////////////////////////////////////////
    ///  Reconstruct the zonal residuals from 
    ///  Shack Hartmann centroid data
    void reconstruct_zernike_residuals(const std::vector<Arroyo::Shack_Hartmann_centroids> & shcentroids, 
				       std::vector<Arroyo::zernike> & znkes) const;

    ///////////////////////////////////////////
    ///  Reconstruct the residuals from 
    ///  Shack Hartmann centroid data
    void reconstruct_residuals(const std::vector<Arroyo::Shack_Hartmann_centroids> & shcentroids,
			       std::vector<Arroyo::zernike> & znke, 
			       std::vector<Arroyo::pixel_array<double> > & pixarrs) const;

    static int verbose_level;
  };


  // Slaving thresholds correspond to areal overlap of edge supapertures.  This may permit
  // us to utilize tiled primaries by computing areal overlap with hexagons.
  // 
  // We need a variable to pass along that specifies whether slaving should be nearest neighbor
  // or komolgorov extrapolation.  The example arcadia.in file uses the latter.
  void write_arcadia_input_file(const Arroyo::refractive_atmospheric_model & ref_atm_model,
				const Arroyo::circular_aperture & circ_ap,
				const std::vector<Arroyo::emitter *> & guide_stars,
				const std::vector<Arroyo::square_lenslet_array> & lenslet_arrays,
				const std::vector<bool> & tilt_removal_flags,
				const std::vector<Arroyo::ideal_tip_tilt_mirror<circular_aperture> > & ttms,
				const std::vector<Arroyo::ideal_deformable_mirror<circular_aperture> > & dms,
				const std::vector<double> & slaving_thresholds,
				const char * filename);

}

#endif
