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

#ifndef NGAO_SYSTEM_DESIGN_MODEL_H
#define NGAO_SYSTEM_DESIGN_MODEL_H

#include "three_vector.h"
#include "refractive_atmosphere.h"

namespace Arroyo {

  ///
  /// A class to implement the refractive atmospheric model used in
  /// the NGAO system design study.  This model consists of a seven
  /// layer atmosphere, with altitudes and relative weights given in
  /// the table below:
  ///
  /// <table border=0 width=200>
  /// <tr class=r1 height=12>
  ///  <td>layer</td>       <td>altitude (km)</td>  <td>weight</td> </tr><tr>
  ///    <td>1</td>           <td>0.00</td>         <td>.47</td>   </tr><tr>  
  ///    <td>2</td>           <td>2.1</td>          <td>.18</td>   </tr><tr>
  ///    <td>3</td>           <td>4.1</td>          <td>.11</td>   </tr><tr>
  ///    <td>4</td>           <td>6.5</td>          <td>.09</td>   </tr><tr>
  ///    <td>5</td>           <td>9.0</td>          <td>.04</td>   </tr><tr>
  ///    <td>6</td>           <td>12.0</td>         <td>.09</td>   </tr><tr>
  ///    <td>7</td>           <td>14.8</td>         <td>.02</td>   </tr></table>
  ///
  ///  The power spectrum assumed for these layers is not stated
  ///  in the paper, so I've assumed Komolgorov.  To instantiate
  ///  this class you must provide a reference three frame

  class NGAO_system_design_model :
    public refractive_atmospheric_model {

    private:
    
    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("NGAO system design model"));};
    
    ///////////////////////////////////////////
    ///  Null constructor
    NGAO_system_design_model(){};

    public:
  
    ///////////////////////////////////////////
    ///  Copy constructor
    NGAO_system_design_model(const NGAO_system_design_model & cn2_model);

    ///////////////////////////////////////////
    ///  Construct from a file
    NGAO_system_design_model(const char * filename);

    ///////////////////////////////////////////
    ///  Construct from an iofits object
    NGAO_system_design_model(const iofits & iof);
    
    ///////////////////////////////////////////
    ///  Construct from a three frame
    ///  The strength of the turbulence may be set by 
    ///  specifying r_0 at a particular wavelength.  
    NGAO_system_design_model(const three_frame & ground_ref_frame,
				  double r_0_meters = .18, 
				  double r_0_ref_wavelength_meters = .5e-6);

    ///////////////////////////////////////////
    ///  Destructor
    ~NGAO_system_design_model(){};

    ///////////////////////////////////////////
    ///  Operator = 
    NGAO_system_design_model & 
      operator=(const NGAO_system_design_model & cn2_model);

    ///////////////////////////////////////////
    /// Clone method
    ///
    /// Calling routine is responsible for deleting memory
    NGAO_system_design_model * clone() const {
      return(new NGAO_system_design_model(*this));
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
