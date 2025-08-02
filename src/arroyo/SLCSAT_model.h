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

#ifndef SLCSAT_MODEL_H
#define SLCSAT_MODEL_H

#include "three_vector.h"
#include "refractive_atmosphere.h"

namespace Arroyo {

  ///
  /// A class to implement a multilayer approximation to the analytic
  /// SLCSAT daytime atmospheric model.  Sasiela p. 61 states the
  /// analytic model as
  ///
  /// <table border=0 width=500>
  /// <tr class=r1 height=12>
  ///    <td></td>             <td>3.96e-13/h^{1.05}</td> <td>18.5 <= h < 232</td>   </tr><tr>
  ///    <td>Cn^{2}(h) = </td> <td>1.3e-15</td>           <td>232 <= h < 880</td>    </tr><tr>
  ///    <td></td>             <td>8.87e-7/h^{3}</td>     <td>880 <= h < 7220</td>   </tr><tr>
  ///    <td></td>             <td>2.0e-16/h^{.5}</td>    <td>7220 <= h < 20500</td> </tr></table>
  ///
  /// where the height h is in meters.
  ///
  /// To generate a multilayer approximation to this model, the user
  /// must supply a vector of layer heights to the constructor.  The
  /// convention for assigning turbulence weights to the layers is as
  /// follows.  The C_{n}^{2} profile is integrated from the ground to
  /// the midpoint between the first and second layer, and this weight
  /// is assigned to the first layer.  The profile is integrated from
  /// this point to the midpoint between the second and third layers,
  /// and this weight is assigned to the second layer.  For the
  /// highest layer, the profile is integrated from the midpoint
  /// between the second highest layer and the highest layer to
  /// infinity.
  ///

  class SLCSAT_day_model :
    public refractive_atmospheric_model {

    private:
    
    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("SLCSAT day model"));};
    
    ///////////////////////////////////////////
    ///  Null constructor
    SLCSAT_day_model(){};

    public:
  
    ///////////////////////////////////////////
    ///  Copy constructor
    SLCSAT_day_model(const SLCSAT_day_model & slcday_model);

    ///////////////////////////////////////////
    ///  Construct from a file
    SLCSAT_day_model(const char * filename);

    ///////////////////////////////////////////
    ///  Construct from an iofits object
    SLCSAT_day_model(const iofits & iof);
    
    ///////////////////////////////////////////
    ///  Construct from a three frame.
    ///
    ///  The layer heights are in meters
    SLCSAT_day_model(const three_frame & ground_ref_frame, 
		     const std::vector<double> & layer_heights);

    ///////////////////////////////////////////
    ///  Destructor
    ~SLCSAT_day_model(){};

    ///////////////////////////////////////////
    ///  Operator = 
    SLCSAT_day_model & 
      operator=(const SLCSAT_day_model & slcday_model);

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

  ///
  /// A class to implement a multilayer approximation to the analytic
  /// SLCSAT nighttime atmospheric model.  Sasiela p. 61 states the
  /// analytic model as
  ///
  /// <table border=0 width=500>
  /// <tr class=r1 height=12>
  ///    <td></td>             <td>5e-15</td>           <td>h<18.5</td>            </tr><tr>  
  ///    <td></td>             <td>2.875e-12/h^{2}</td> <td>18.5 <= h < 110</td>   </tr><tr>
  ///    <td>Cn^{2}(h) = </td> <td>2.5e-16</td>         <td>110 <= h < 850</td>    </tr><tr>
  ///    <td></td>             <td>8.87e-7/h^{3}</td>   <td>850 <= h < 7000</td>   </tr><tr>
  ///    <td></td>             <td>2.0e-16/h^{.5}</td>  <td>7000 <= h < 20500</td> </tr></table>
  ///
  /// where the height h is in meters.  (Note - there is an inconsistency
  /// in Sasiela's statement of the model.  Check original source).
  ///
  /// To generate a multilayer approximation to this model, the user
  /// must supply a vector of layer heights to the constructor.  The
  /// convention for assigning turbulence weights to the layers is as
  /// follows.  The C_{n}^{2} profile is integrated from the ground to
  /// the midpoint between the first and second layer, and this weight
  /// is assigned to the first layer.  The profile is integrated from
  /// this point to the midpoint between the second and third layers,
  /// and this weight is assigned to the second layer.  For the
  /// highest layer, the profile is integrated from the midpoint
  /// between the second highest layer and the highest layer to
  /// infinity.
  ///

  class SLCSAT_night_model :
    public refractive_atmospheric_model {

    private:
    
    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("SLCSAT night model"));};
    
    ///////////////////////////////////////////
    ///  Null constructor
    SLCSAT_night_model(){};

    public:
  
    ///////////////////////////////////////////
    ///  Copy constructor
    SLCSAT_night_model(const SLCSAT_night_model & slcnight_model);

    ///////////////////////////////////////////
    ///  Construct from a file
    SLCSAT_night_model(const char * filename);

    ///////////////////////////////////////////
    ///  Construct from an iofits object
    SLCSAT_night_model(const iofits & iof);
    
    ///////////////////////////////////////////
    ///  Construct from a three frame.
    ///
    ///  The layer heights are in meters
    SLCSAT_night_model(const three_frame & ground_ref_frame, 
		       const std::vector<double> & layer_heights);

    ///////////////////////////////////////////
    ///  Destructor
    ~SLCSAT_night_model(){};

    ///////////////////////////////////////////
    ///  Operator = 
    SLCSAT_night_model & 
      operator=(const SLCSAT_night_model & slcnight_model);

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
