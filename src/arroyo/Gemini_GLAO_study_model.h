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

#ifndef GEMINI_GLAO_MODELS_H
#define GEMINI_GLAO_MODELS_H

#include "three_vector.h"
#include "refractive_atmosphere.h"

namespace Arroyo {

  ///
  /// A class to instantiate one of the nine refractive atmospheric
  /// models used in the Gemini GLAO study.
  ///

  class Gemini_GLAO_study_model :
    public refractive_atmospheric_model {

    private:
    
    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("Gemini GLAO study model"));};
    
    ///////////////////////////////////////////
    ///  Null constructor
    Gemini_GLAO_study_model(){};

    public:
  
    ///////////////////////////////////////////
    ///  Copy constructor
    Gemini_GLAO_study_model(const Gemini_GLAO_study_model & cn2_model);

    ///////////////////////////////////////////
    ///  Construct from a file
    Gemini_GLAO_study_model(const char * filename);

    ///////////////////////////////////////////
    ///  Construct from an iofits object
    Gemini_GLAO_study_model(const iofits & iof);
    
    ///////////////////////////////////////////
    ///  Construct from a three frame and
    ///  a pair of strings.  The strings take
    ///  one of three values:
    ///  "good", "typical" or "bad"
    Gemini_GLAO_study_model(const three_frame & tf,
			    const string & ground_layer,
			    const string & focal_anisoplanatism,
			    bool extended_profile = false,
			    double outer_scale = -1);

    ///////////////////////////////////////////
    ///  Destructor
    ~Gemini_GLAO_study_model(){};

    ///////////////////////////////////////////
    ///  Operator = 
    Gemini_GLAO_study_model & 
      operator=(const Gemini_GLAO_study_model & cn2_model);

    ///////////////////////////////////////////
    /// Clone method
    ///
    /// Calling routine is responsible for deleting memory
    Gemini_GLAO_study_model * clone() const {
      return(new Gemini_GLAO_study_model(*this));
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
