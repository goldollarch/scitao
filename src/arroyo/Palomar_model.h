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

#ifndef PALOMAR_MODEL_H
#define PALOMAR_MODEL_H

#include "three_vector.h"
#include "refractive_atmosphere.h"

namespace Arroyo {

  ///
  /// A class to hold a model for the turbulence profile
  /// from Palomar, formulated from 2 months of DIMM/MASS
  /// data
  ///

  class Palomar_DIMM_MASS_model :
    public refractive_atmospheric_model {

    private:
    
    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("Gemini GLAO study model"));};
    
    ///////////////////////////////////////////
    ///  Null constructor
    Palomar_DIMM_MASS_model(){};

    public:
  
    ///////////////////////////////////////////
    ///  Copy constructor
    Palomar_DIMM_MASS_model(const Palomar_DIMM_MASS_model & cn2_model);

    ///////////////////////////////////////////
    ///  Construct from a file
    Palomar_DIMM_MASS_model(const char * filename);

    ///////////////////////////////////////////
    ///  Construct from an iofits object
    Palomar_DIMM_MASS_model(const iofits & iof);
    
    ///////////////////////////////////////////
    ///  Construct from a three frame and
    ///  a pair of strings.  The strings take
    ///  one of three values:
    ///  "good", "typical" or "bad"
    Palomar_DIMM_MASS_model(const three_frame & tf,
			    double r_0_meters = .087, 
			    double r_0_ref_wavelength_meters = .5e-6);

    ///////////////////////////////////////////
    ///  Destructor
    ~Palomar_DIMM_MASS_model(){};

    ///////////////////////////////////////////
    ///  Operator = 
    Palomar_DIMM_MASS_model & 
      operator=(const Palomar_DIMM_MASS_model & cn2_model);

    ///////////////////////////////////////////
    /// Clone method
    ///
    /// Calling routine is responsible for deleting memory
    Palomar_DIMM_MASS_model * clone() const {
      return(new Palomar_DIMM_MASS_model(*this));
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
    void print(ostream & os, const char * prefix="") const;

  };

}

#endif
