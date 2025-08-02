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

#ifndef Cn2_DISTRIBUTION_H
#define Cn2_DISTRIBUTION_H

#include <iostream>
#include <vector>

namespace Arroyo {

  using std::vector;
  using std::ostream;

  ///
  ///  A virtual base class to hold Cn2 distributions
  ///

  class Cn2_distribution {

  private:

    ///////////////////////////////////////////
    ///  Make a copy of this Cn2 distribution.
    ///  You are responsible for calling delete
    ///  on the instance returned by this function 
    virtual Cn2_distribution * clone() const = 0;

  public:

    ///////////////////////////////////////////
    ///  Null constructor
    Cn2_distribution(){};

    ///////////////////////////////////////////
    ///  Virtual destructor
    virtual ~Cn2_distribution(){};

    ///////////////////////////////////////////
    ///  Return a moment of this turbulent distribution.
    ///  This moment is calculated as the integral of 
    ///  Cn2 multiplied by the altitude raised to the exponent
    ///  power.  
    ///
    ///  This function throws an error if the exponent is negative.
    virtual double get_moment(double exponent) const = 0;

    ///////////////////////////////////////////
    ///  Get the value of the Cn2 distribution at
    ///  a particular height.  
    ///
    ///  Height is in meters
    virtual double val(double height) const = 0;

    ///////////////////////////////////////////
    ///  Print the turbulence distribution
    virtual void print(ostream & os, const char * prefix = "") const = 0;

    ///////////////////////////////////////////
    ///  Determine if this Cn2 distribution is
    ///  equal to another.  Equality means the
    ///  two instances are of the same derived
    ///  type and have the same data members
    virtual bool equals(Cn2_distribution * cn2_dist) const = 0;

    ///////////////////////////////////////////
    ///  Factory to construct a Cn2_distribution from another
    static Cn2_distribution * Cn2_distribution_factory(const Cn2_distribution * cn2dist){
      return cn2dist->clone(); 
    };

  };

  ///  A class to hold the Cn2 distribution model from 
  ///  Hardy's book "Adaptive Optics for Astronomical 
  ///  Telescopes", eq. 3.19

  class Hardy_Cn2_distribution :
    public Cn2_distribution {

    private:

    ///////////////////////////////////////////
    ///  Make a copy of the Cn2 distribution.
    ///  This copy must be deleted or a memory
    ///  leak will result
    Cn2_distribution * clone() const;

    ///////////////////////////////////////////
    ///  Check if two Cn2 distributions are equal.
    ///  This check is polymorphic
    bool equals(Cn2_distribution * cn2_dist) const;

    protected:

    /// The boundary layer coefficient, in units of m^{-2/3}
    double boundary_coeff_;
    /// The boundary layer scale height, in meters
    double boundary_scale_;

    /// The tropospheric layer coefficient, in units of m^{-2/3}
    double troposphere_coeff_;
    /// The tropospheric layer scale height, in meters
    double troposphere_scale_;

    /// The tropopause coefficient, in units of m^{-2/3}
    double tropopause_coeff_;
    /// The tropopause height, in meters
    double tropopause_height_;

    /// A vector containing information on additional layers (coeff, height, thickness)
    vector<vector<double> > additional_layers_;

    public:

    ///////////////////////////////////////////
    ///  Null constructor
    Hardy_Cn2_distribution();

    ///////////////////////////////////////////
    ///  Copy constructor
    Hardy_Cn2_distribution(const Hardy_Cn2_distribution & tcn2d);

    ///////////////////////////////////////////
    ///  Construct from the bits
    Hardy_Cn2_distribution(double boundary_coeff, double boundary_scale,
			   double troposphere_coeff, double troposphere_scale,
			   double tropopause_coeff, double tropopause_height,
			   vector<vector<double> > additional_layers = vector<vector<double> >(0));

    ///////////////////////////////////////////
    ///  Destructor
    ~Hardy_Cn2_distribution(){};

    ///////////////////////////////////////////
    ///  Add a turbulent layer to the model
    ///
    ///  Coeff is in m^{-2/3}
    ///  Height is in meters
    ///  Thickness is in meters
    void add_layer(double coeff, double height, double thickness);

    ///////////////////////////////////////////
    ///  Return a moment of this turbulent distribution.
    ///  This moment is calculated as the integral of 
    ///  Cn2 multiplied by the altitude raised to the exponent
    ///  power.  
    ///
    ///  This function throws an error if the exponent is negative.
    double get_moment(double exponent) const;

    ///////////////////////////////////////////
    ///  Get the value of the Cn2 distribution
    ///  at a particular height
    double val(double height) const;

    ///////////////////////////////////////////
    ///  Print a description of the Cn2 distribution
    void print(ostream & os, const char * prefix = "") const;

    ///////////////////////////////////////////
    ///  Check if two Cn2 distributions are equal.
    friend bool operator==(const Hardy_Cn2_distribution & t1, const Hardy_Cn2_distribution & t2);
  };


  /// A class to hold the SLCSAT night model
  /// for the Cn2 distribution.  From
  /// Sasiela's book "Electromagnetic Wave 
  /// Propagation in Turbulence" eq. 4.12.

  class SLCSAT_day_Cn2_distribution :
    public Cn2_distribution {

    private:

    ///////////////////////////////////////////
    ///  Check if two Cn2 distributions are equal.
    ///  This check is polymorphic
    Cn2_distribution * clone() const;
  
    ///////////////////////////////////////////
    ///  Check if two Cn2 distributions are equal.
    bool equals(Cn2_distribution * cn2_dist) const;

    public:

    ///////////////////////////////////////////
    ///  Null constructor
    SLCSAT_day_Cn2_distribution(){};

    ///////////////////////////////////////////
    ///  Destructor
    ~SLCSAT_day_Cn2_distribution(){};

    ///////////////////////////////////////////
    ///  Return a moment of this turbulent distribution.
    ///  This moment is calculated as the integral of 
    ///  Cn2 multiplied by the altitude raised to the exponent
    ///  power.  
    ///
    ///  This function throws an error if the exponent is negative.
    double get_moment(double exponent) const;

    ///////////////////////////////////////////
    ///  Get the value of the Cn2 distribution
    ///  at a particular height
    double val(double height) const;

    ///////////////////////////////////////////
    ///  Print a description of the Cn2 distribution
    void print(ostream & os, const char * prefix = "") const;

  };

  /// A Class to hold the SLCSAT night model
  /// for the Cn2 distribution.  From
  /// Sasiela's book "Electromagnetic Wave 
  /// Propagation in Turbulence" eq. 4.13.
  /// Note that this equation has an inconsistency
  /// in the limits 

  class SLCSAT_night_Cn2_distribution :
    public Cn2_distribution {

    private:

    ///////////////////////////////////////////
    ///  Check if two Cn2 distributions are equal.
    ///  This check is polymorphic
    Cn2_distribution * clone() const;

    ///////////////////////////////////////////
    ///  Check if two Cn2 distributions are equal.
    bool equals(Cn2_distribution * cn2_dist) const;

    public:

    ///////////////////////////////////////////
    ///  Null constructor
    SLCSAT_night_Cn2_distribution(){};

    ///////////////////////////////////////////
    ///  Destructor
    ~SLCSAT_night_Cn2_distribution(){};

    ///////////////////////////////////////////
    ///  Return a moment of this turbulent distribution.
    ///  This moment is calculated as the integral of 
    ///  Cn2 multiplied by the altitude raised to the exponent
    ///  power.  
    ///
    ///  This function throws an error if the exponent is negative.
    double get_moment(double exponent) const;

    ///////////////////////////////////////////
    ///  Get the value of the Cn2 distribution
    ///  at a particular height
    double val(double height) const;

    ///////////////////////////////////////////
    ///  Print a description of the Cn2 distribution
    void print(ostream & os, const char * prefix = "") const;

  };

}

#endif

