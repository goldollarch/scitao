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

#ifndef ATMOSPHERIC_MODEL_H
#define ATMOSPHERIC_MODEL_H


namespace Arroyo {

  // forward class declaration
  class spherical_wave_emitter;

  ///
  /// A class to represent a sodium atmospheric model.
  ///

  class sodium_atmospheric_model {

  public:
  
    ///////////////////////////////////////////
    ///  Null constructor
    sodium_atmospheric_model();

    ///////////////////////////////////////////
    ///  Copy constructor
    sodium_atmospheric_model(const sodium_atmospheric_model & ref_atm_model);

    ///////////////////////////////////////////
    ///  Construct from a file
    sodium_atmospheric_model(const char * filename);

    ///////////////////////////////////////////
    ///  Construct from an iofits object
    sodium_atmospheric_model(iofits & iof);

    ///////////////////////////////////////////
    ///  Destructor
    ~sodium_atmospheric_model();

    ///////////////////////////////////////////
    ///  Operator = 
    sodium_atmospheric_model & operator=(const sodium_atmospheric_model & ref_atm_model);

    ///////////////////////////////////////////
    ///  Read from file
    void read(const char * filename);

    ///////////////////////////////////////////
    ///  Read from iofits
    void read(const iofits & iof);

    ///////////////////////////////////////////
    ///  Write to file
    void write(const char * filename) const;

    ///////////////////////////////////////////
    ///  Write to iofits
    void write(const iofits & iof) const;

    ///////////////////////////////////////////
    ///  Return an emitter by backscattering a wavefront
    spherical_wave_emitter backscatter(const wavefront & wf) const;

  };

  ///
  /// A class to represent a rayleigh_mie atmospheric model.
  ///

  class rayleigh_mie_atmospheric_model {

  public:
  
    ///////////////////////////////////////////
    ///  Null constructor
    rayleigh_mie_atmospheric_model();

    ///////////////////////////////////////////
    ///  Copy constructor
    rayleigh_mie_atmospheric_model(const rayleigh_mie_atmospheric_model & ref_atm_model);

    ///////////////////////////////////////////
    ///  Construct from a file
    rayleigh_mie_atmospheric_model(const char * filename);

    ///////////////////////////////////////////
    ///  Construct from an iofits object
    rayleigh_mie_atmospheric_model(iofits & iof);

    ///////////////////////////////////////////
    ///  Destructor
    ~rayleigh_mie_atmospheric_model();

    ///////////////////////////////////////////
    ///  Operator = 
    rayleigh_mie_atmospheric_model & operator=(const rayleigh_mie_atmospheric_model & ref_atm_model);

    ///////////////////////////////////////////
    ///  Read from file
    void read(const char * filename);

    ///////////////////////////////////////////
    ///  Read from iofits
    void read(const iofits & iof);

    ///////////////////////////////////////////
    ///  Write to file
    void write(const char * filename) const;

    ///////////////////////////////////////////
    ///  Write to iofits
    void write(const iofits & iof) const;

    ///////////////////////////////////////////
    ///  Return an emitter by backscattering a wavefront
    spherical_wave_emitter backscatter(const wavefront & wf) const;

  };

}

#endif
