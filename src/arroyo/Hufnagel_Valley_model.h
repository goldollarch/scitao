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

#ifndef HUFNAGEL_VALLEY_MODEL_H
#define HUFNAGEL_VALLEY_MODEL_H

#include "three_vector.h"
#include "refractive_atmosphere.h"

namespace Arroyo {

  ///
  /// A class to implement a multilayer approximation to the analytic
  /// atmospheric model proposed by Hufnagel and Valley.  Paraphrasing
  /// Sasiela p. 60, the analytic model is
  ///
  /// <table border=0 width=500>
  /// <tr class=r1 height=12>
  ///    <td></td>             <td>.00594 (W/27)^{2} (1e-5*h)^{10}exp(-h/1000) + </td>    </tr><tr>
  ///    <td>Cn^{2}(h) = </td> <td>2.7e-16 exp(-h/1500) +</td>   </tr><tr>
  ///    <td></td>             <td>A exp(-h/100)</td> </tr></table>
  ///
  /// where W is the pseudowind and A is a parameter that is usually
  /// set equal to 1.7e-14.  The HV-21 model has this value for A, and
  /// W is equal to 21.  This model is sometimes referred to as the HV
  /// 5/7 model.  (Note - h is apparently in meters)
  ///
  /// Sasiela tabulates the turbulence moments for 4 different values of 
  /// the pseudowind W.  They are reproduced in the table below.
  ///
  /// <table border=0 width=500>
  /// <tr class=r1 height=12>
  ///  <td></td>                 <td>HV 21</td>  <td>HV 35</td>  <td>HV 54</td> <td>HV 72</td> </tr><tr>
  ///    <td>r_0 (cm)</td>       <td>4.96</td>   <td>4.67</td>   <td>4.18</td>  <td>3.70</td>  </tr><tr> 
  ///    <td>theta_0 (urad)</td> <td>6.90</td>   <td>3.95</td>   <td>2.40</td>  <td>1.71</td>  </tr><tr> <td></td> <td></td> <td></td> <td></td> <td></td> </tr>
  /// </table>
  ///
  /// where r_0 and theta_0 are quoted at .5 microns.
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
  /// Note that this procedure does not necessarily preserve theta_0
  /// in the table above, but does preserve r_0.
  ///

  class Hufnagel_Valley_model :
    public refractive_atmospheric_model {

    private:
    
    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("Hufnagel Valley model"));};
    
    ///////////////////////////////////////////
    ///  Null constructor
    Hufnagel_Valley_model(){};

    public:
  
    ///////////////////////////////////////////
    ///  Copy constructor
    Hufnagel_Valley_model(const Hufnagel_Valley_model & hv_model);

    ///////////////////////////////////////////
    ///  Construct from a file
    Hufnagel_Valley_model(const char * filename);

    ///////////////////////////////////////////
    ///  Construct from an iofits object
    Hufnagel_Valley_model(const iofits & iof);
    
    ///////////////////////////////////////////
    ///  Construct from a three frame.
    Hufnagel_Valley_model(const three_frame & ground_ref_frame,
			  const std::vector<double> & layer_heights,
			  double pseudowind, double A = 1.7e-14);

    ///////////////////////////////////////////
    ///  Destructor
    ~Hufnagel_Valley_model(){};

    ///////////////////////////////////////////
    ///  Operator = 
    Hufnagel_Valley_model & 
      operator=(const Hufnagel_Valley_model & hv_model);

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
