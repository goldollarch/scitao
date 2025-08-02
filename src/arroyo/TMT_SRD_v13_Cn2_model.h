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

#ifndef TMT_SRD_v13_CN2_MODEL_H
#define TMT_SRD_v13_CN2_MODEL_H

#include "three_vector.h"
#include "refractive_atmosphere.h"

namespace Arroyo {

  ///
  /// A class to implement the refractive atmospheric model
  /// proposed by the Science Advisory Committe for the Thirty
  /// Meter Telescope project in version 13 of the TMT Science 
  /// Requirements Document.  This model consists of a seven 
  /// layer atmosphere, with altitudes and relative weights 
  /// given in the table below:
  ///
  /// <table border=0 width=200>
  /// <tr class=r1 height=12>
  ///  <td>layer</td>       <td>altitude (km)</td>  <td>weight</td> </tr><tr>
  ///    <td>1</td>           <td>0.00</td>         <td>.646</td>   </tr><tr>  
  ///    <td>2</td>           <td>1.8</td>          <td>.078</td>   </tr><tr>
  ///    <td>3</td>           <td>3.3</td>          <td>.119</td>   </tr><tr>
  ///    <td>4</td>           <td>5.8</td>          <td>.035</td>   </tr><tr>
  ///    <td>5</td>           <td>7.4</td>          <td>.025</td>   </tr><tr>
  ///    <td>6</td>           <td>13.1</td>         <td>.082</td>   </tr><tr>
  ///    <td>7</td>           <td>15.8</td>         <td>.015</td>   </tr></table>
  ///
  ///  The power spectrum assumed for these layers is not stated
  ///  in the paper, so I've assumed Komolgorov.  To instantiate
  ///  this class you must provide a reference three frame

  class TMT_SRD_v13_Cn2_model :
    public refractive_atmospheric_model {

    private:
    
    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("TMT SRD v13 Cn2 model"));};
    
    ///////////////////////////////////////////
    ///  Null constructor
    TMT_SRD_v13_Cn2_model(){};

    public:
  
    ///////////////////////////////////////////
    ///  Copy constructor
    TMT_SRD_v13_Cn2_model(const TMT_SRD_v13_Cn2_model & cn2_model);

    ///////////////////////////////////////////
    ///  Construct from a file
    TMT_SRD_v13_Cn2_model(const char * filename);

    ///////////////////////////////////////////
    ///  Construct from an iofits object
    TMT_SRD_v13_Cn2_model(const iofits & iof);
    
    ///////////////////////////////////////////
    ///  Construct from a three frame
    ///  The strength of the turbulence may be set by 
    ///  specifying r_0 at a particular wavelength.  
    ///  The nominal Cn2 integral of 3.535e-13 m^{1/3}
    ///  specified in version 13 of the TMT SRD
    ///  corresponds to an r0 of 14.99294 cm at .5 microns.
    ///  This is the default.
    TMT_SRD_v13_Cn2_model(const three_frame & ground_ref_frame,
				  double r_0_meters = .1499294, 
				  double r_0_ref_wavelength_meters = .5e-6);

    ///////////////////////////////////////////
    ///  Destructor
    ~TMT_SRD_v13_Cn2_model(){};

    ///////////////////////////////////////////
    ///  Operator = 
    TMT_SRD_v13_Cn2_model & 
      operator=(const TMT_SRD_v13_Cn2_model & cn2_model);

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
