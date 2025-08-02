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

#ifndef LENSLET_ARRAY_H
#define LENSLET_ARRAY_H

#include "optic.h"

namespace Arroyo {

  using std::string;
  using std::vector;
  using std::ostream;

  ///
  /// A virtual base class to represent a lenslet array.
  ///

  class lenslet_array_base :
    public plane_optic,
    public one_to_one_optic {

    public:

    ///////////////////////////////////////////
    ///  Null constructor
    lenslet_array_base(){};

    ///////////////////////////////////////////
    ///  Copy constructor
    lenslet_array_base(const lenslet_array_base & labase);

    ///////////////////////////////////////////
    ///  Virtual destructor
    virtual ~lenslet_array_base(){};

    ///////////////////////////////////////////
    ///  Write to an iofits object
    virtual void write(iofits & iof) const = 0;

    ///////////////////////////////////////////
    ///  Operator = 
    lenslet_array_base & operator=(const lenslet_array_base & labase);

    ///////////////////////////////////////////
    ///  Factory to construct lenslet_array_bases from file
    static lenslet_array_base * lenslet_array_base_factory(const char * filename);

    ///////////////////////////////////////////
    ///  Factory to construct lenslet_array_bases from file
    static lenslet_array_base * lenslet_array_base_factory(const iofits & iof);

  };

  ///
  /// A class to represent an ideal square lenslet array.
  ///
  /// This class contains a subtle implementation of the
  /// one_to_one_optic::transform member function.  The
  /// usual implementation leaves the wavefront at the same
  /// location in space, but changes the amplitude and/or phase
  /// values to effect the transformation.  In a lenslet
  /// array, this paradigm is not so good because the
  /// wavefront pixels may straddle the lenslet array boundaries.
  /// Thus, straight modification of the wavefront phase would 
  /// result in unwelcome edge effects.
  ///
  /// To avoid this issue, I've included in the transform member 
  /// function the free space propagation to the far field
  /// of the lenslet array.  This propagation is effected using
  /// the Goertzel-Reinsch algorithm, and wavefront pixels that
  /// straddle boundaries are properly weighted.  By default
  /// the propagation leaves the wavefront in the focal plane
  /// of the lenslet array, with Nyquist sampling.  These
  /// parameters may be changed by using the member functions
  ///
  /// set_final_wavefront_propagation_distance(double)
  /// set_final_wavefront_pixels_per_lenslet(long)
  /// set_final_wavefront_pixels_per_transform(long)
  ///

  class square_lenslet_array :
    public lenslet_array_base {

    private:

    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("square lenslet array"));};
    
    /// A template member function to perform
    /// transform on both float and double
    /// instantiations of wavefront.  This 
    /// is necessary because there is no 
    /// mechanism in C++ for virtual template
    /// member functions.
    template<class T>
      void private_transform(diffractive_wavefront<T> & wf) const;

    protected:

    /// Number of lenslets across array
    vector<long> axes;

    /// Focal length of each lenslet, in meters
    double focal_length;
    
    /// Physical dimensions of each lenslet, in meters
    double lenslet_pitch;

    /// Parameter to set the location where
    /// transform member function leaves
    /// the wavefront.  Measured in meters
    /// from the lenslet array.  This parameter
    /// should be chosen to leave the array
    /// in the far field.  This parameter
    /// defaults to the focal length of the 
    /// lenslet array.
    double final_wavefront_propagation_distance;
    
    /// Axes to use for the final wavefront
    /// in the transform member function.
    long final_wavefront_pixels_per_lenslet;

    /// Oversampling factor to use for the final 
    /// wavefrontin the transform member function.  
    /// This parameter is unity for Nyquist sampling,
    /// which is the default.
    long final_wavefront_pixels_per_transform;

    ///////////////////////////////////////////
    ///  Null constructor
    square_lenslet_array(){};

    public:

    ///////////////////////////////////////////
    ///  Copy constructor
    square_lenslet_array(const square_lenslet_array & sq_lns_arr);

    ///////////////////////////////////////////
    ///  Construct from file
    square_lenslet_array(const char * filename);

    ///////////////////////////////////////////
    ///  Construct from iofits object
    square_lenslet_array(const iofits & iof);

    ///////////////////////////////////////////
    ///  Construct from the bits
    ///
    ///  array_axes is a vector containing the
    ///  number of lenslets along each axis
    ///
    ///  flength is the focal lenth of the lenslets,
    ///  in meters.
    ///  
    ///  lnslt_pitch is the size of the square lenslet,
    ///  in meters
    square_lenslet_array(vector<long> array_axes,
			 double flength,
			 double lnslt_pitch,
			 long pix_per_lenslet,
			 long pix_per_xform);

    ///////////////////////////////////////////
    ///  Destructor
    ~square_lenslet_array(){};

    ///////////////////////////////////////////
    ///  Operator = 
    square_lenslet_array & operator=(const square_lenslet_array & sq_lns_arr);

    ///////////////////////////////////////////
    ///  Read from file
    void read(const char * filename);

    ///////////////////////////////////////////
    ///  Read from an iofits object
    void read(const iofits & iof);
 
    ///////////////////////////////////////////
    ///  Write to file
    void write(const char * filename) const;

    ///////////////////////////////////////////
    ///  Write to an iofits object
    void write(iofits & iof) const;

    ///////////////////////////////////////////
    ///  Print
    void print(ostream & os, const char * prefix="") const;

    ///////////////////////////////////////////
    ///  Get the final wavefront propagation distance, in meters
    double get_final_wavefront_propagation_distance() const {return final_wavefront_propagation_distance;};
   
    ///////////////////////////////////////////
    ///  Set the final wavefront propagation distance.
    ///
    ///  wf_prop_dist is in meters
    void set_final_wavefront_propagation_distance(double wf_prop_dist) {final_wavefront_propagation_distance = wf_prop_dist;};
   
    ///////////////////////////////////////////
    ///  Get the final wavefront pixels per lenslet
    long get_final_wavefront_pixels_per_lenslet() const {return final_wavefront_pixels_per_lenslet;};
   
    ///////////////////////////////////////////
    ///  Set the final wavefront pixels per lenslet
    void set_final_wavefront_pixels_per_lenslet(long pix_per_lnslt);
   
    ///////////////////////////////////////////
    ///  Get the final wavefront pixels per transform
    long get_final_wavefront_pixels_per_transform() const {return final_wavefront_pixels_per_transform;};
   
    ///////////////////////////////////////////
    ///  Set the final wavefront pixels per transform
    void set_final_wavefront_pixels_per_transform(long pix_per_xform);
   
    ///////////////////////////////////////////
    ///  Get the lenslet axes
    vector<long> get_axes() const {return(axes);};
   
    ///////////////////////////////////////////
    ///  Get the lenslet pitch
    double get_lenslet_pitch() const {return(lenslet_pitch);};
   
    ///////////////////////////////////////////
    ///  Get the focal length
    double get_focal_length() const {return(focal_length);};
   
    ///////////////////////////////////////////
    ///  Get a rectangular region guaranteed to cover
    ///  the lenslet array.  The resulting region will
    ///  have its edges aligned with the x and y axes
    ///  of the three_frame tf.  
    ///
    ///  If the z axis of the three_frame is orthogonal
    ///  to the z axis of the aperture, this function
    ///  throws an error
    rectangular_region get_covering_region(const three_frame & tf) const;

    ///////////////////////////////////////////
    ///  Apply the lenslet array to the wavefront
    void transform(diffractive_wavefront<float> & wf) const;

    ///////////////////////////////////////////
    ///  Apply the aperture to the wavefront
    void transform(diffractive_wavefront<double> & wf) const;

  };
}

#endif
