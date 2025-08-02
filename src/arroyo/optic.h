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

#ifndef OPTIC_H
#define OPTIC_H

#include "AO_sim_base.h"
#include "fits_header_data.h"
#include "three_frame.h"

namespace Arroyo {

  using std::ostream;
  using std::vector;
  using std::string;

  class iofits;
  class rectangular_region;
  class three_vector;
  class three_scaling;
  class wavefront_header;
  template<typename T> class diffractive_wavefront;

  ///
  /// A virtual base class for an optic.
  ///

  class optic :
    virtual public AO_sim_base {

    protected:
    
    /// \brief
    /// A flag to indicate whether to account for 
    /// foreshortening that results when a wavefront
    /// strikes an optic at an oblique angle.
    /// Default is true
    bool foreshortening;

    public:

    ///////////////////////////////////////////
    ///  Null constructor
    optic(){foreshortening = true;};

    ///////////////////////////////////////////
    ///  Copy constructor
    optic(const optic & op);

    ///////////////////////////////////////////
    ///  Virtual destructor
    virtual ~optic(){};

    ///////////////////////////////////////////
    ///  Operator = 
    optic & operator=(const optic & op);

    //////////////////////////////////////////// 
    /// Query whether the optic appears 
    /// foreshortened to a wavefront incident at
    /// an angle to the optic
    bool get_foreshortening() const {return(foreshortening);};

    //////////////////////////////////////////// 
    /// Choose whether the optic appears 
    /// foreshortened to a wavefront incident at
    /// an angle to the optic
    void set_foreshortening(bool fshrtn) {foreshortening = fshrtn;};

    ///////////////////////////////////////////
    ///  Read from an iofits object
    virtual void read(const iofits & iof);
 
    ///////////////////////////////////////////
    ///  Write to an iofits object
    virtual void write(iofits & iof) const;

    ///////////////////////////////////////////
    ///  Print
    virtual void print(ostream & os, const char * prefix="") const;

    ///////////////////////////////////////////
    ///  Get a rectangular region guaranteed to cover
    ///  the entire optic.  The resulting region will
    ///  have its edges aligned with the x and y axes
    ///  of the three_frame tf, regardless of the
    ///  foreshortening status of the optic.  However,
    ///  if foreshortening is on then the rectangular
    ///  region returned by this function covers the
    ///  optic when projected onto the plane containing
    ///  the optic.  If foreshortening is off, then
    ///  this function returns a region that covers the
    ///  optic when rotated into the plane of the optic.
    ///  That is, when we just assume the z axis defining
    ///  the plane containing the region is the same
    ///  as the z axis defining the optic.)
    ///
    virtual rectangular_region get_covering_region(const three_frame & tf) const = 0;

    ///////////////////////////////////////////
    ///  Get the point of intersection of
    ///  a line extending from three_point tp
    ///  in the direction of the three_vector
    ///  tv and this optic.  If there is no
    ///  intersection point, this function throws
    ///  an error
    virtual three_point get_point_of_intersection(const three_point & tp, 
						  const three_vector & tv) const = 0;

    ///////////////////////////////////////////
    ///  Verbose level
    static int verbose_level;

    ///////////////////////////////////////////
    ///  Factory to construct optics from file
    static optic * optic_factory(const char * filename);

    ///////////////////////////////////////////
    ///  Factory to construct optics from file
    static optic * optic_factory(const iofits & iof);

  };

  ///
  /// A virtual base class for a plane_optic.
  ///

  class plane_optic :
    virtual public optic, 
    public three_frame {

    protected:

    /// A simple function to find the projected wavefront
    /// pixel spacing in the plane containing the plane_optic
    template<class T>
      void get_projected_wavefront_pixel_spacing(const diffractive_wavefront<T> & wf, 
						 three_vector & origin_offset, 
						 three_vector & dx, 
						 three_vector & dy) const;
  
    public:

    ///////////////////////////////////////////
    ///  Null constructor
    plane_optic(){};

    ///////////////////////////////////////////
    ///  Copy constructor
    plane_optic(const plane_optic & plane_op);

    ///////////////////////////////////////////
    ///  Virtual destructor
    virtual ~plane_optic(){};

    ///////////////////////////////////////////
    ///  Operator = 
    plane_optic & operator=(const plane_optic & plane_op);

    ///////////////////////////////////////////
    ///  Read from an iofits object
    virtual void read(const iofits & iof);
 
    ///////////////////////////////////////////
    ///  Write to an iofits object
    virtual void write(iofits & iof) const;

    ///////////////////////////////////////////
    ///  Print
    virtual void print(ostream & os, const char * prefix="") const;

    ///////////////////////////////////////////
    ///  Get the point of intersection of
    ///  a line extending from three_point tp
    ///  in the direction of the three_vector
    ///  tv and this optic.  If there is no
    ///  intersection point, this function throws
    ///  an error
    virtual three_point get_point_of_intersection(const three_point & tp, 
						  const three_vector & tv) const;

    ///////////////////////////////////////////
    ///  Factory to construct optics from file
    static plane_optic * plane_optic_factory(const char * filename);

    ///////////////////////////////////////////
    ///  Factory to construct optics from file
    static plane_optic * plane_optic_factory(const iofits & iof);

  };

  template<class T>
  void plane_optic::get_projected_wavefront_pixel_spacing(const diffractive_wavefront<T> & wf, 
							  three_vector & origin_offset, 
							  three_vector & dx, 
							  three_vector & dy) const {

    // First, ensure that wavefront normal is not orthogonal to the
    // plane_optic normal
    if(fabs(dot_product(wf.three_frame::z(), this->three_frame::z()))<three_frame::precision){
      cerr << "plane_optic::get_projected_wavefront_pixel_spacing error - "
	   << "diffractive_wavefront orthogonal to circular plane_optic normal\n";
      throw(string("plane_optic::get_projected_wavefront_pixel_spacing"));
    }
  
    // Next, ensure that the wavefront origin lies in the transverse
    // frame of the optic.
    origin_offset = static_cast<three_point>(wf) - static_cast<const three_point>(*this);

    if(fabs(dot_product(origin_offset, this->three_frame::z()))>three_frame::precision){
      cerr << "plane_optic::get_projected_wavefront_pixel_spacing error - "
	   << "diffractive_wavefront center not in transverse plane of circular plane_optic\n";
      throw(string("plane_optic::get_projected_wavefront_pixel_spacing"));
    }

    double wf_pixel_scale = wf.get_pixel_scale();
    dx = wf.three_frame::x()*wf_pixel_scale;
    dy = wf.three_frame::y()*wf_pixel_scale;

    try{
      dx = parallel_projection(dx, this->z(), wf.z());
      dy = parallel_projection(dy, this->z(), wf.z());
    } catch(...){
      cerr << "plane_optic::get_projected_wavefront_pixel_spacing error - "
	   << "could not form parallel projection\n";
      throw(string("plane_optic::get_projected_wavefront_pixel_spacing"));
    }
    
    // If there's no foreshortening, we explicitly set the length to be
    // equal to the wf pixel scale
    if(!foreshortening){
      dx = (wf_pixel_scale/dx.length())*dx;
      dy = (wf_pixel_scale/dy.length())*dy;
    }
  }


  ///
  /// A virtual base class for an optic
  /// that maps one wavefront to another.
  ///

  class one_to_one_optic :
    virtual public optic {

    protected:

    ///////////////////////////////////////////
    ///  Return pointer to the raw wavefront data
    ///  array.  
    ///
    ///  This member function is used in the one_to_one_optic inheritance
    ///  hierarchy to implement the virtual member function
    ///
    ///  void optic::one_to_one_transform(diffractive_wavefront<T> &)
    ///
    ///  It is intended that this member function modify the wavefront
    ///  data directly to apply the optical transformation.
    template<class T>
      T * get_wavefront_data(diffractive_wavefront<T> & wf) const {
      return(wf.wfdata);
    }

    ///////////////////////////////////////////
    ///  return storage method - real-imag or
    ///  amp-phase
    ///
    ///  This member function is used in the optic inheritance hierarchy
    ///  to implement the virtual member function 
    ///  void one_to_one_optic::transform(diffractive_wavefront<T>)
    template<class T>
      bool is_real_imag_storage(const diffractive_wavefront<T> & wf) const {
      return(wf.real_imag);
    }

    ///////////////////////////////////////////
    ///  return storage method - interleaved
    ///  or non-interleaved 
    ///
    ///  This member function is used in the optic inheritance hierarchy
    ///  to implement the virtual member function 
    ///  void one_to_one_optic::transform(diffractive_wavefront<T>)
    template<class T>
      bool is_interleaved_storage(const diffractive_wavefront<T> & wf) const {
      return(wf.interleaved);
    }

    ///////////////////////////////////////////
    ///  Convert wavefront to real/imag storage.
    ///  This conversion is idempotent.
    template<class T>
      void wavefont_real_imag_conversion(diffractive_wavefront<T> & wf) const {
      wf.real_imag_conversion();
    }
 
    ///////////////////////////////////////////
    ///  Convert wavefront to amp/phase storage.
    ///  This conversion is idempotent.
    template<class T>
      void wavefront_amp_phase_conversion(diffractive_wavefront<T> & wf) const {
      wf.amp_phase_conversion();
    }

    ///////////////////////////////////////////
    ///  return pointer to geometric rays
    ///
    ///  This member function is used in the optic inheritance hierarchy
    ///  to implement the virtual member function 
    ///  void one_to_one_optic::transform(geometric_wavefront)
    // geometric_ray * get_wavefront_data(const geometric_wavefront & gwf) const;

    public:

    ///////////////////////////////////////////
    ///  Null constructor
    one_to_one_optic(){};

    ///////////////////////////////////////////
    ///  Virtual destructor
    virtual ~one_to_one_optic(){};

    ///////////////////////////////////////////
    ///  Virtual member function to apply this
    ///  optic to a diffractive_wavefront<float>
    virtual void transform(diffractive_wavefront<float> & wf) const = 0;

    ///////////////////////////////////////////
    ///  Virtual member function to apply this
    ///  optic to a diffractive_wavefront<double>
    virtual void transform(diffractive_wavefront<double> & wf) const = 0;

    ///////////////////////////////////////////
    ///  Virtual member function to apply this
    ///  optic to a geometric_wavefront
    //virtual void transform(geometric_wavefront & gwf) const = 0;

  };

  ///
  /// A virtual base class for an optic
  /// that maps one wavefront to many
  /// wavefronts.

  class one_to_many_optic :
    virtual public optic {

    protected:

    ///////////////////////////////////////////
    ///  Return pointer to the raw wavefront data
    ///  array.  
    ///
    ///  This member function is used in the one_to_one_optic inheritance
    ///  hierarchy to implement the virtual member function
    ///
    ///  vector<diffractive_wavefront<T> > one_to_many_optic::transform(const diffractive_wavefront<T> &)
    ///
    ///  It is intended that the above member function leave the
    ///  diffractive_wavefront data unchanged, instead constructing diffractive_wavefronts
    ///  corresponding to the multiple outputs of this optic.  To help
    ///  enforce this, this member function returns a const pointer to a
    ///  const.  You shouldn't need to cast away this constness to
    ///  implement the behavior described above
    template<class T>
      const T * const get_wavefront_data(const diffractive_wavefront<T> & wf) const {
      return(wf.wfdata);
    }

    ///////////////////////////////////////////
    ///  return storage method - real-imag or
    ///  amp-phase
    ///
    ///  This member function is used in the optic inheritance hierarchy
    ///  to implement the virtual member function vector<diffractive_wavefront<T> >
    ///  one_to_many_optic::transform(const diffractive_wavefront<T>)
    template<class T>
      bool real_imag_storage(const diffractive_wavefront<T> & wf) const {
      return(wf.real_imag);
    }

    ///////////////////////////////////////////
    ///  return storage method - interleaved
    ///  or non-interleaved 
    ///
    ///  This member function is used in the optic inheritance hierarchy
    ///  to implement the virtual member function vector<diffractive_wavefront<T> >
    ///  one_to_many_optic::transform(const diffractive_wavefront<T>)
    template<class T>
      bool interleaved_storage(const diffractive_wavefront<T> & wf) const {
      return(wf.interleaved);
    }

    public:
  
    ///////////////////////////////////////////
    ///  Null constructor
    one_to_many_optic(){};

    ///////////////////////////////////////////
    ///  Virtual destructor
    virtual ~one_to_many_optic(){};

    ///////////////////////////////////////////
    ///  Virtual member function to report
    ///  the number of wavefronts returned
    ///  by optic::transform
    virtual long number_of_outputs() const = 0;

    ///////////////////////////////////////////
    ///  Virtual member function to return
    ///  the wavefront header from the nth
    ///  output that will result from transforming
    ///  a wavefront with header wfhdr 
    ///virtual wavefront_header output_header(const wavefront_header & input_wfhdr, 
    ///long output_index) const = 0;

    ///////////////////////////////////////////
    ///  Virtual member function to apply this
    ///  optic to a diffractive_wavefront<double>
    virtual vector<diffractive_wavefront<float> > transform(const diffractive_wavefront<float> & wf) const = 0;

    ///////////////////////////////////////////
    ///  Virtual member function to apply this
    ///  optic to a diffractive_wavefront<double>
    virtual vector<diffractive_wavefront<double> > transform(const diffractive_wavefront<double> & wf) const = 0;

    ///////////////////////////////////////////
    ///  Virtual member function to apply this
    ///  optic to a geometric_ray
    //virtual vector<geometric_wavefront> transform(const geometric_wavefront & gwf) const = 0;

  };

  /// others:  mirror, oap, pyramid, 
  ///          shearing interferometer,
  ///          tiled hexagonal aperture
  ///          deformable mirror

}

#endif

