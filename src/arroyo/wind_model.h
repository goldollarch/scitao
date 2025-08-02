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

#ifndef WIND_MODEL_H
#define WIND_MODEL_H

#include <string>
#include <vector>
#include "three_frame.h"
#include "AO_sim_base.h"

namespace Arroyo {

  using std::string;

  ///
  /// A base class to represent atmospheric wind models
  ///

  class wind_model :
    public AO_sim_base {

    public:
  
    ///////////////////////////////////////////
    ///  Null constructor
    wind_model(){};
    
    ///////////////////////////////////////////
    ///  Virtual destructor
    virtual ~wind_model(){};

    ///////////////////////////////////////////
    /// Virtual clone method
    ///
    /// Calling routine is responsible for deleting memory
    virtual wind_model * clone() const = 0;
    
    ///////////////////////////////////////////
    ///  Read from file
    virtual void read(const char * filename) = 0;

    ///////////////////////////////////////////
    ///  Read from iofits
    virtual void read(const iofits & iof) = 0;

    ///////////////////////////////////////////
    ///  Write to file
    virtual void write(const char * filename) const = 0;

    ///////////////////////////////////////////
    ///  Write to iofits
    virtual void write(iofits & iof) const = 0;

    ///////////////////////////////////////////
    ///  Print
    virtual void print(ostream & os, const char * prefix) const = 0;

    ///////////////////////////////////////////
    ///  Retrieve a set of random wind vectors
    ///  from the model.  The wind vectors are 
    ///  potentially composed of contributions from
    ///  a number of sources (e.g. ground layer, 
    ///  tropospheric layer, etc.).  The wind
    ///  model determines at what level the wind
    ///  vectors are correlated in height.  Because
    ///  of this, the wind vectors for a
    ///  particular realization of a layered 
    ///  atmosphere should be retrieved at the 
    ///  same time, so that these correlations may
    ///  be enforced by the wind model.
    ///
    ///  wind vectors are measured in meters per second.
    ///
    ///  Height is in meters
    virtual std::vector<Arroyo::three_vector> get_random_wind_vectors(const std::vector<double> & heights, 
								      const Arroyo::three_frame & ref_frame) const = 0;
    
    ///////////////////////////////////////////
    ///  Factory to construct wind models from file
    static wind_model * wind_model_factory(const char * filename);

    ///////////////////////////////////////////
    ///  Factory to construct wind models from an iofits object
    static wind_model * wind_model_factory(const iofits & iof);

    ///////////////////////////////////////////
    ///  Factory to construct wind models from another instance
    static wind_model * wind_model_factory(const wind_model & wm);

    static int verbose_level;

  };


  ///
  /// A class to represent the atmospheric wind model
  /// described by Hardy eq. 3.20
  ///

  class Hardy_wind_model :
    public wind_model {

    private:
    
    static const bool factory_registration;
    
    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("Hardy wind model"));};
    
    protected:

    ///////////////////////////////////////////
    ///  Null constructor
    Hardy_wind_model();
    
    ///  The wind velocity of the ground layer
    double ground_layer_wind_velocity;

    ///  The wind velocity at the tropopause
    double tropopause_wind_velocity;

    ///  The height of the tropopause
    double tropopause_height;

    ///  The thickness of the tropopause
    double tropopause_thickness;

    public:
  
    ///////////////////////////////////////////
    ///  destructor
    ~Hardy_wind_model(){};

    ///////////////////////////////////////////
    ///  Copy constructor
    Hardy_wind_model(const Hardy_wind_model & twm);
    
    ///////////////////////////////////////////
    ///  Construct from a file
    Hardy_wind_model(const char * filename);
    
    ///////////////////////////////////////////
    ///  Construct from an iofits object
    Hardy_wind_model(const iofits & iof);
    
    ///////////////////////////////////////////
    ///  Operator=
    Hardy_wind_model & operator=(const Hardy_wind_model & twm);
    
    ///////////////////////////////////////////
    ///  Construct from the bits
    ///
    ///  The ground layer velocity is in meters per second
    ///
    ///  The tropopause velocity is in meters per second
    ///
    ///  The tropopause height is in meters
    ///
    ///  The tropopause thickness is in meters
    Hardy_wind_model(double grnd_lyr_velocity,
		     double trpse_height = 0,
		     double trpse_thickness = 0,
		     double trpse_velocity = 0);

    ///////////////////////////////////////////
    /// Virtual clone method
    ///
    /// Calling routine is responsible for deleting memory
    Hardy_wind_model * clone() const {
      return new Hardy_wind_model(*this);
    };
    
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
    void write(iofits & iof) const;

    ///////////////////////////////////////////
    ///  Print
    void print(ostream & os, const char * prefix) const;

    ///////////////////////////////////////////
    ///  Retrieve a set of random wind vectors
    ///  from the Hardy wind model.  These wind 
    ///  vectors contain a random ground layer
    ///  component and a component arising from
    ///  the tropospheric layer, which is modulated
    ///  by a gaussian in height.  Thus, there
    ///  are two random components to the wind vector
    ///  that are mapped to the layer heights in 
    ///  a way that depends on the tropospheric 
    ///  parameters of the model.
    ///
    ///  For a given atmospheric realization, you should
    ///  get all the wind vectors at the same time, since
    ///  these two random components determine the wind
    ///  vectors at all altitudes.
    ///
    ///  wind vectors are measured in meters per second.
    ///
    ///  Height is in meters
    std::vector<Arroyo::three_vector> get_random_wind_vectors(const std::vector<double> & heights, 
							      const Arroyo::three_frame & ref_frame) const;
    
  };

}

#endif
