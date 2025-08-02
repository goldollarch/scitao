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

#ifndef APLUSPLUS_RECONSTRUCTOR_H
#define APLUSPLUS_RECONSTRUCTOR_H

#include "zernike_projected_zonal_reconstructor.h"

namespace Arroyo {

  ///
  /// A class to hold reconstructors generated using A++
  /// 
  /// This class aims to represent reconstructors that yield zonal
  /// residuals as well as an arbitrary number of zernike mode
  /// residuals.
  ///

  class a_plusplus_reconstructor :
    public zernike_projected_zonal_reconstructor,
    public Arroyo::pixel_array<double> {

    private:

    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    std::string unique_name() const {return(std::string("A++ reconstructor"));};

    protected:

    /// A map from reconstructor array indices to zernike modes.  This
    /// data member maps the reconstructor index corresponding to the
    /// output residual of the modal part of the reconstructor to a
    /// specific zernike mode.  The zernike mode is labelled by a
    /// pair: the first element is a char (either c or s) and the
    /// second element is a vector with two elements (zernike order,
    /// level).  This mapping allows this class to support
    /// reconstructors that return residuals for any number of
    /// possibly noncontiguous zernike modes.
    std::map<int, std::pair<char, std::vector<long> > > reconstructor_to_zernike_mode_map;

    ///////////////////////////////////////////
    ///  Null constructor
    a_plusplus_reconstructor(){};

    public:

    ///////////////////////////////////////////
    ///  Copy constructor
    a_plusplus_reconstructor(const a_plusplus_reconstructor & appr);

    ///////////////////////////////////////////
    ///  Construct from file
    a_plusplus_reconstructor(const char * filename);

    ///////////////////////////////////////////
    ///  Construct from an A++ ascii file
    /// 
    ///  You must specify the file containing the reconstructor, along
    ///  with the number of slopes (twice the number of subaps) and
    ///  the number of actuators.  
    ///
    ///  Note that the A++ reconstructor must be fully populated with
    ///  subapertures and actuators throughout the grid.  This is due
    ///  to the fact that Arroyo assumes an instance of the Shack
    ///  Hartmann centroid class will be passed into this function,
    ///  and this class populates a rectangular grid, rather than
    ///  leaving out the unilluminated subaps.  Similarly, the zonal
    ///  residuals returned by the reconstructor lie on a rectangular
    ///  grid, meaning that residuals for uncontrolled actuators are
    ///  returned by this reconstructor.  One of the specific
    ///  consequences of this assumption is that user-defined wiring
    ///  schemes in A++ are not supported by this class.
    ///
    ///  This constructor attempts to perform all checks that it
    ///  possibly can, as the A++ interface really must be considered
    ///  fragile and largely untested.
    a_plusplus_reconstructor(const char * APP_zonal_reconstructor_file, 
			     int nactuators, int nslopes);

    ///////////////////////////////////////////
    ///  Construct from two A++ ascii files
    /// 
    ///  You must specify the file containing the reconstructor with
    ///  the zernike modes projected out, along with the number of
    ///  slopes (twice the number of subaps) and the number of
    ///  actuators.  You must also specify the file containing the
    ///  reconstructor that converts slopes to zernike modes, along
    ///  with the number of modes in the projection.  This reconstructor
    ///  assumes that the zernike modes projected out of the reconstructor
    ///  run from lowest to highest using the A++ counting scheme, and 
    ///  exclude Bal. SpAb (mode number 9).
    ///
    ///  It is not necessarily the case that one will always want to
    ///  make a reconstructor that generates the zonal residuals with
    ///  the zernike modes projected out.  For example, one could call
    ///  this reconstructor by supplying an APP zonal reconstructor
    ///  file without any zernike mode projection together with an APP
    ///  zernike reconstructor that generates the first 10 zernike
    ///  modes.  This would allow one to monitor statistical
    ///  fluctuations in these modes by examining the time histories
    ///  of their coefficients.
    ///
    ///  Note that the A++ reconstructor must be fully populated with
    ///  subapertures and actuators throughout the grid.  This is due
    ///  to the fact that Arroyo assumes an instance of the Shack
    ///  Hartmann centroid class will be passed into this function,
    ///  and this class populates a rectangular grid, rather than
    ///  leaving out the unilluminated subaps.  Similarly, the zonal
    ///  residuals returned by the reconstructor lie on a rectangular
    ///  grid, meaning that residuals for uncontrolled actuators are
    ///  returned by this reconstructor.  One of the specific
    ///  consequences of this assumption is that user-defined wiring
    ///  schemes in A++ are not supported by this class.
    ///
    ///  This constructor attempts to perform all checks that it
    ///  possibly can, as the A++ interface really must be considered
    ///  fragile and largely untested. 
    a_plusplus_reconstructor(const char * APP_zonal_reconstructor_file, 
					  int nactuators, int nslopes,
					  const char * APP_zernike_reconstructor_file, 
					  int nzernike_modes);

    ///////////////////////////////////////////
    ///  Construct from iofits object
    a_plusplus_reconstructor(const Arroyo::iofits & iof);

    ///////////////////////////////////////////
    ///  Destructor
    ~a_plusplus_reconstructor(){};

    ///////////////////////////////////////////
    ///  Operator = 
    a_plusplus_reconstructor & operator=(const a_plusplus_reconstructor & appr);

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
    ///  Get dimensions of centroid measurements passed
    ///  to the reconstructor
    vector<long> get_centroid_axes() const;

    ///////////////////////////////////////////
    ///  Get dimensions of pixel array returned
    ///  by the reconstructor
    vector<long> get_actuator_axes() const;

    ///////////////////////////////////////////
    ///  Get a zernike instance that contains information
    ///  about which modes are returned by the reconstructor.
    /// 
    ///  This instance is minimally sized so as to hold the largest
    ///  mode returned by the reconstructor.  Each element of this
    ///  instance is initialized to unity if the corresponding mode is
    ///  returned by the reconstructor, and to zero if it is not
    Arroyo::zernike get_zernike_modes() const;

    ///////////////////////////////////////////
    ///  Reconstruct the zonal residuals from 
    ///  Shack Hartmann centroid data
    ///
    /// This reconstructor reconstructs a zernike instance
    /// with a number of modes specified by the reconstructor
    void reconstruct_zernike_residuals(const Arroyo::Shack_Hartmann_centroids & shcentroids, 
				       Arroyo::zernike & znke) const;

    ///////////////////////////////////////////
    ///  Reconstruct the zonal residuals from a
    ///  Shack Hartmann centroid class instance
    void reconstruct_zonal_residuals(const Arroyo::Shack_Hartmann_centroids & shcentroids, 
				     Arroyo::pixel_array<double> & pixarr) const;

    ///////////////////////////////////////////
    ///  Reconstruct the residuals from 
    ///  Shack Hartmann centroid data
    ///
    /// This reconstructor reconstructs a zernike instance with a
    /// number of modes specified by the reconstructor, and a pixel
    /// array containing zonal residuals
    void reconstruct_residuals(const Arroyo::Shack_Hartmann_centroids & shcentroids, 
			       Arroyo::zernike & znke, 
			       Arroyo::pixel_array<double> & pixarr) const;
  };


}

#endif
