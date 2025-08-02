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

#ifndef AO_SIMULATION_H
#define AO_SIMULATION_H

#include "refractive_atmosphere.h"

namespace Arroyo {

  using std::vector;

  ///  A class to hold the information necessary to
  ///  specify a whole simulation.  This involves
  ///  a lot of information.

  class ao_simulation {

  protected:

    /// The refractive atmospheric model
    refractive_atmospheric_model * ref_atm_model;

    /// The sodium atmospheric model
    sodium_atmospheric_model * sodium_atm_model;

    /// The rayleigh-mie atmospheric model
    rayleigh_mie_atmospheric_model * rayliegh_mie_atm_model;

    /// Random realization of the sodium atmosphere
    sodium_atmosphere * sodium_atm;

    /// Random realization of the rayleigh mie atmosphere
    rayleigh_mie_atmosphere * rayleigh_mie_atm;

    /// The simulation step interval, in seconds
    double simulation_step_interval;

    /// Optical systems used in the simulation
    vector<optical_system> op_sys;

    /// Wavefront sensors used in the simulation
    vector<wavefront_sensor> wfsensors;

    /// Readout times for the wavefront sensors,
    /// measured in units of the simulation step
    /// interval
    vector<long> wfs_readout_times;

    /// The science camera detector
    detector science_camera;

    /// Readout times for the science camera,
    /// measured in units of the simulation step
    /// interval
    long science_camera_readout_time;

    /// Emitters in space, the header information for the
    /// wavefronts that they will emit, the optical
    /// trajectories for these wavefronts, and the detectors 
    /// that will receive these wavefronts.
    /// Each vector has a number of elements equal to the
    /// number of extraterrestrial emitters
    /// The pointers to optical_systems address instances
    /// in the op_sys vector above
    /// The pointers to detectors address instances
    /// in the wfsensors vector above or the science_camera
    /// vector<plane_wave_emitter> extraterrestrial_emitters;
    /// vector<wavefront_header> extraterrestrial_wavefront_headers;
    /// vector<vector<optical_system *> > extraterrestrial_optics;
    /// vector<detector *> extraterrestrial_detectors;

    /// Laser emitters, the header information for the wavefronts
    /// that they will emit, the optical trajectories for these
    /// wavefronts, and the heights for which backscattered 
    /// wavefronts will be propagated to the detectors.
    /// Header information for the backscattered wavefronts,
    /// optical trajectories for these wavefronts, and the
    /// detectors for these wavefronts.
    /// Each vector has a number of elements equal to the
    /// number of laser emitters
    /// The pointers to optical_systems address instances
    /// in the op_sys vector above
    /// The pointers to detectors address instances
    /// in the wfsensors vector above or the science_camera
    /// vector<emitter> laser_emitters;
    /// vector<wavefront_header> laser_wavefront_headers;
    /// vector<vector<optical_system *> > laser_optics;
    /// vector<vector<double> > laser_backscatter_heights;
    /// vector<wavefront_header> laser_backscatter_wavefront_headers;
    /// vector<vector<optical_system *> > laser_backscatter_optics;
    /// vector<detector *> laser_backscatter_detectors;

    /// the output plan for the simulation
    output_plan sim_output;

  public:

    ///////////////////////////////////////////
    ///  Null constructor
    ao_simulation();

    ///////////////////////////////////////////
    ///  Copy constructor
    ao_simulation(const ao_simulation & ao_sim);

    ///////////////////////////////////////////
    ///  Construct from a file
    ao_simulation(const char * filename);

    ///////////////////////////////////////////
    ///  Destructor
    ~ao_simulation();

    ///////////////////////////////////////////
    ///  Operator = 
    ao_simulation & operator=(const ao_simulation & ao_sim);

  }

}

#endif
