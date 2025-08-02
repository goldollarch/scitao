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

#ifndef ELLERBROEK_CERRO_PACHON_MODEL_H
#define ELLERBROEK_CERRO_PACHON_MODEL_H

#include "three_vector.h"
#include "refractive_atmosphere.h"

namespace Arroyo {

  ///
  /// A class to implement the refractive atmospheric model
  /// proposed by Ellerbroek, JOSA v 19, p 1803.  This model
  /// consists of a six layer atmosphere, with altitudes and
  /// relative weights given in the table below:
  ///
  /// <table border=0 width=200>
  /// <tr class=r1 height=12>
  ///  <td>layer</td>       <td>altitude (km)</td>   <td>weight</td> </tr><tr>
  ///    <td>1</td>           <td>0.00</td>          <td>.652</td>   </tr><tr>  
  ///    <td>2</td>           <td>2.58</td>          <td>.172</td>   </tr><tr>
  ///    <td>3</td>           <td>5.16</td>          <td>.055</td>   </tr><tr>
  ///    <td>4</td>           <td>7.73</td>          <td>.025</td>   </tr><tr>
  ///    <td>5</td>           <td>12.89</td>         <td>.074</td>   </tr><tr>
  ///    <td>6</td>           <td>15.46</td>         <td>.022</td>   </tr></table>
  ///
  ///  The power spectrum assumed for these layers is not stated
  ///  in the paper, so I've assumed Komolgorov.  To instantiate
  ///  this class you must provide a reference three frame

  class Ellerbroek_Cerro_Pachon_model :
    public refractive_atmospheric_model {

    private:
    
    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("Ellerbroek Cerro Pachon model"));};
    
    ///////////////////////////////////////////
    ///  Null constructor
    Ellerbroek_Cerro_Pachon_model(){};

    public:
  
    ///////////////////////////////////////////
    ///  Copy constructor
    Ellerbroek_Cerro_Pachon_model(const Ellerbroek_Cerro_Pachon_model & e_cp_model);

    ///////////////////////////////////////////
    ///  Construct from a file
    Ellerbroek_Cerro_Pachon_model(const char * filename);

    ///////////////////////////////////////////
    ///  Construct from an iofits object
    Ellerbroek_Cerro_Pachon_model(const iofits & iof);
    
    ///////////////////////////////////////////
    ///  Construct from a three frame
    ///  The strength of the turbulence may be set by 
    ///  specifying r_0 at a particular wavelength.  
    ///  These values default to those proposed by Brent 
    ///  Ellerbroek and Francois Rigaut in SPIE 4007 p 1088
    ///  With these values, the isoplanatic angle is 
    ///  2.29 arcsecs
    Ellerbroek_Cerro_Pachon_model(const three_frame & ground_ref_frame,
				  double r_0_meters = .166, 
				  double r_0_ref_wavelength_meters = .5e-6);

    ///////////////////////////////////////////
    ///  Destructor
    ~Ellerbroek_Cerro_Pachon_model(){};

    ///////////////////////////////////////////
    ///  Operator = 
    Ellerbroek_Cerro_Pachon_model & 
      operator=(const Ellerbroek_Cerro_Pachon_model & e_cp_model);

    ///////////////////////////////////////////
    /// Clone method
    ///
    /// Calling routine is responsible for deleting memory
    Ellerbroek_Cerro_Pachon_model * clone() const {
      return(new Ellerbroek_Cerro_Pachon_model(*this));
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
