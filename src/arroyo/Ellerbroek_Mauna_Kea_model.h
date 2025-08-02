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

#ifndef ELLERBROEK_MAUNA_KEA_MODEL_H
#define ELLERBROEK_MAUNA_KEA_MODEL_H

#include "three_vector.h"
#include "refractive_atmosphere.h"

namespace Arroyo {

  ///
  /// A class to implement the refractive atmospheric model
  /// proposed by Ellerbroek.  This model
  /// consists of an eleven layer atmosphere, with altitudes and
  /// relative weights given in the table below:
  ///
  /// <table border=0 width=200>
  /// <tr class=r1 height=12>
  ///  <td>layer</td>       <td>altitude (km)</td>  <td>weight</td> </tr><tr>
  ///    <td>1</td>           <td>.09</td>          <td>.003</td>   </tr><tr>  
  ///    <td>2</td>           <td>1.826</td>        <td>.136</td>   </tr><tr>
  ///    <td>3</td>           <td>2.72</td>         <td>.163</td>   </tr><tr>
  ///    <td>4</td>           <td>4.256</td>        <td>.161</td>   </tr><tr>
  ///    <td>5</td>           <td>6.269</td>        <td>.167</td>   </tr><tr>
  ///    <td>6</td>           <td>8.34</td>         <td>.235</td>   </tr><tr>
  ///    <td>7</td>           <td>10.546</td>       <td>.068</td>   </tr><tr>
  ///    <td>8</td>           <td>12.375</td>       <td>.032</td>   </tr><tr>
  ///    <td>9</td>           <td>14.61</td>        <td>.023</td>   </tr><tr>
  ///    <td>10</td>          <td>16.471</td>       <td>.006</td>   </tr><tr>
  ///    <td>11</td>          <td>17.028</td>       <td>.007</td>   </tr></table>
  ///
  ///  The power spectrum assumed for these layers is not stated
  ///  in the paper, so I've assumed Komolgorov.  To instantiate
  ///  this class you must provide a reference three frame

  class Ellerbroek_Mauna_Kea_model :
    public refractive_atmospheric_model {

    private:
    
    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("Ellerbroek Mauna Kea model"));};
    
    ///////////////////////////////////////////
    ///  Null constructor
    Ellerbroek_Mauna_Kea_model(){};

    public:
  
    ///////////////////////////////////////////
    ///  Copy constructor
    Ellerbroek_Mauna_Kea_model(const Ellerbroek_Mauna_Kea_model & e_cp_model);

    ///////////////////////////////////////////
    ///  Construct from a file
    Ellerbroek_Mauna_Kea_model(const char * filename);

    ///////////////////////////////////////////
    ///  Construct from an iofits object
    Ellerbroek_Mauna_Kea_model(const iofits & iof);
    
    ///////////////////////////////////////////
    ///  Construct from a three frame.
    ///  The strength of the turbulence may be set by 
    ///  specifying r_0 at a particular wavelength.  
    ///  These values default to those proposed by Brent 
    ///  Ellerbroek and Francois Rigaut in SPIE 4007 p 1088.
    ///  With these values, the isoplanatic angle is 
    ///  2.74 arcsecs
    Ellerbroek_Mauna_Kea_model(const three_frame & ground_ref_frame, 
			       double r_0_meters = .236, 
			       double r_0_ref_wavelength_meters = .5e-6);

    ///////////////////////////////////////////
    ///  Destructor
    ~Ellerbroek_Mauna_Kea_model(){};

    ///////////////////////////////////////////
    ///  Operator = 
    Ellerbroek_Mauna_Kea_model & 
      operator=(const Ellerbroek_Mauna_Kea_model & e_cp_model);

    ///////////////////////////////////////////
    /// Clone method
    ///
    /// Calling routine is responsible for deleting memory
    Ellerbroek_Mauna_Kea_model * clone() const {
      return(new Ellerbroek_Mauna_Kea_model(*this));
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
