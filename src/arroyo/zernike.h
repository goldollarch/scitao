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

#ifndef ZERNIKE_H
#define ZERNIKE_H

#include <string>
#include <iostream>
#include <vector>
#include "AO_sim_base.h"
#include "pixel_array.h"
#include "modal_expansion.h"

namespace Arroyo {
  class iofits;
}

namespace Arroyo {

  using std::string;
  class aperture;
  class circular_aperture;

  ///    
  /// A class to hold zernike polynomial coefficients.
  /// The polynomials used for this class are those of
  /// Born & Wolf - section 9.2.1 equations 5 and 10.
  /// This class uses the real polynomial
  /// expansion in which radial polynomials are multiplied
  /// by cos(m*phi) and sin(m*phi).
  /// For each order n of the expansion, there are n/2+1
  /// possible levels m, given by m = n%2, n%2+2, ... n
  /// for both cos and sin terms.  However, for m==0 there
  /// is no sin term.  
  ///
  /// Coefficients are stored in the one dimensional 
  /// pixel_array inherited by this class.
  ///

  class zernike :
    public AO_sim_base,
    public modal_expansion, 
    public Arroyo::pixel_array<double> {

    private:

    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("zernike expansion"));};

    ///////////////////////////////////////////
    ///  Returns a pixel_array with the
    ///  given axes and pixel scale, and defined
    ///  on the circle specified by aperture.
    ///  This function throws an error if the 
    ///  aperture is not a circular aperture.
    Arroyo::pixel_array<double> get_pixel_array(const std::vector<long> & axes, 
						double pixscale,
						const aperture & ap) const;

    ///////////////////////////////////////////
    ///  Function to resize zernike to a particular order.
    void resize(long order);

    public:

    ///////////////////////////////////////////
    ///  Null constructor
    zernike();

    ///////////////////////////////////////////
    ///  Construct a null instance sized to hold
    ///  order orders of zernike polynomials.  
    zernike(long order){
      this->resize(order);
    };

    ///////////////////////////////////////////
    ///  Destructor
    ~zernike(){};

    ///////////////////////////////////////////
    ///  Copy constructor
    zernike(const zernike & iznke){
      this->operator=(iznke);
    };

    ///////////////////////////////////////////
    ///  Construct from file
    zernike(const char * filename){
      this->read(filename);
    };

    ///////////////////////////////////////////
    ///  Construct from an iofits object
    zernike(const iofits & iof){
      this->read(iof);
    };

    ///////////////////////////////////////////
    ///  Construct from a pixel array
    ///  by expanding the pixel array up 
    ///  to the given order.  The pixel array
    ///  is assumed to be bound by the circular
    ///  aperture circ_ap, with pixscale specifying the
    ///  pixel scale.
    ///
    ///  Currently, pixels on the edge whose centers are
    ///  outside the aperture are not included in the 
    ///  expansion, regardless of whether circ_ap has
    ///  its areal_weighting flag set 
    zernike(const Arroyo::pixel_array<double> & pixarr, 
	    double pixscale,
	    long order, 
	    const circular_aperture & circ_ap);

    ///////////////////////////////////////////
    ///  Operator =
    ///
    ///  Forced to do this because compiler cannot
    ///  seem to find this function itself.
    zernike & operator=(const Arroyo::pixel_array<double> & pixarr) {
      this->Arroyo::pixel_array<double>::operator=(pixarr);
      return(*this);
    };

    ///////////////////////////////////////////
    ///  Operator =
    zernike & operator=(const zernike & znke);

    ///////////////////////////////////////////
    ///  read from a file
    void read(const char * filename);

    ///////////////////////////////////////////
    ///  read from an iofits object
    void read(const Arroyo::iofits & iof);

    ///////////////////////////////////////////
    ///  write to a file
    void write(const char * filename) const;

    ///////////////////////////////////////////
    ///  write to an iofits object
    void write(Arroyo::iofits & iof) const;

    ///////////////////////////////////////////
    ///  Function to get a coefficient of particular order and level to cos zernike
    double get_cos_coeff(int order, int level) const;

    ///////////////////////////////////////////
    ///  Function to set a coefficient of particular order and level to cos zernike
    void set_cos_coeff(int order, int level, double coeff);

    ///////////////////////////////////////////
    ///  Function to get a coefficient of particular order and level to sin zernike
    double get_sin_coeff(int order, int level) const;

    ///////////////////////////////////////////
    ///  Function to set a coefficient of particular order and level to sin zernike
    void set_sin_coeff(int order, int level, double coeff);

    ///////////////////////////////////////////
    ///  Function to print the coefficients
    void print(std::ostream & os, const char * prefix = "") const;

    //////////////////////////////////////////
    ///  Function to return the largest radial 
    ///  order in the expansion.  
    ///  
    ///  If there are no terms in the expansion,
    ///  this function returns -1
    long get_order() const;

    ///////////////////////////////////////////
    ///  Returns a pixel_array with the
    ///  given axes and pixel scale, defined
    ///  by the circular aperture circ_ap.
    ///
    ///  If circ_ap has its areal_weighting
    ///  flag set, then this function will take
    ///  into account pixels that lie on the boundary
    ///  but whose centers lie outside the aperture.
    ///  In these cases the phase will be set to the
    ///  value given by the expansion.
    Arroyo::pixel_array<double> get_pixel_array(const std::vector<long> & axes,
						double pixscale,
						const circular_aperture & circ_ap) const;


    ///  Verbose level
    static int verbose_level;

  };
}

#endif
