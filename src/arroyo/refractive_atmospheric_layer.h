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

#ifndef REFRACTIVE_ATMOSPHERIC_LAYER_H
#define REFRACTIVE_ATMOSPHERIC_LAYER_H

#include <iostream>
#include "region_base.h"
#include "power_spectrum.h"
#include "pixel_array.h"
#include "optic.h"
#include "diffractive_wavefront.h"

namespace Arroyo {

  using std::string;
  using std::vector;
  using std::ostream;

  template<typename U> class diffractive_wavefront;
  //class geometric_wavefront;

  ///
  ///  A class to represent a refractive atmospheric layer
  ///

  template<class T>
  class refractive_atmospheric_layer :
    public plane_optic,
    public one_to_one_optic, 
    public pixel_array<T> {

    private:

    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("refractive atmospheric layer"));};

    protected:
    /// The pixel scale in meters
    double pixel_scale;

    /// A vector in the direction of the layer's motion
    /// Through the constructor this vector is always coaligned
    /// with one of the axes of the optical_path_differences array
    three_vector wind_vector;

    /// A template member function to perform
    /// transform on both float and double
    /// instantiations of wavefront.  This 
    /// is necessary because there is no 
    /// mechanism in C++ for virtual template
    /// member functions.
    template<class U>
    void private_transform(diffractive_wavefront<U> & wf) const;

    /// A template member function to perform
    /// an aligned transform on both float and double
    /// instantiations of wavefront.  This 
    /// function throws an error if the transverse
    /// axes of the wavefront projected onto the
    /// transverse plane of the layer are not
    /// aligned with the transverse axes of the 
    /// layer.  This condition is met only if either
    /// the x or the y basis vectors of each frame are
    /// identical, or both are identical.
    template<class U>
    void aligned_private_transform(diffractive_wavefront<U> & wf) const;

    public:

    ///////////////////////////////////////////
    ///  Null constructor
    refractive_atmospheric_layer();

    ///////////////////////////////////////////
    ///  Copy constructor
    refractive_atmospheric_layer(const refractive_atmospheric_layer<T> & ref_atm_layer);

    ///////////////////////////////////////////
    ///  Construct from a subregion of another
    ///  refractive atmospheric layer
    refractive_atmospheric_layer(const refractive_atmospheric_layer<T> & ref_atm_layer, 
				 const rectangular_region & subregion);

    ///////////////////////////////////////////
    ///  Construct from file
    refractive_atmospheric_layer<T>(const char * filename);

    ///////////////////////////////////////////
    ///  Construct from iofits object
    refractive_atmospheric_layer<T>(const iofits & iof);

    ///////////////////////////////////////////
    ///  Construct from the bits
    refractive_atmospheric_layer<T>(const pixel_array<T> & pixarr, 
				    double pixel_scale);

    ///////////////////////////////////////////
    ///  Construct from a power spectrum and a 
    ///  subharmonic method
    refractive_atmospheric_layer(const power_spectrum * pspectrum,
				 const subharmonic_method & subm,
				 const vector<long> & axes, 
				 double pixscale);

    ///////////////////////////////////////////
    ///  Construct a refractive_atmospheric_layer 
    ///  of smaller size from the argument
    ///  The corners argument must be a vector 
    ///  containing 4 three_points, which form a 
    ///  rectangular region contained by ref_atm_layer
    refractive_atmospheric_layer<T>(
    		const refractive_atmospheric_layer<T> & ref_atm_layer,
		vector<three_point> & corners);

    ///////////////////////////////////////////
    ///  Destructor
    ~refractive_atmospheric_layer(){};

    ///////////////////////////////////////////
    ///  Operator=
    refractive_atmospheric_layer<T> & operator=(
    		const refractive_atmospheric_layer<T> & ref_atm_layer);

    ///////////////////////////////////////////
    ///  Read from file
    void read(const char * filename);

    ///////////////////////////////////////////
    ///  Read from iofits object
    void read(const iofits & iof);

    ///////////////////////////////////////////
    ///  Write to file
    void write(const char * filename) const;

    ///////////////////////////////////////////
    ///  Write to iofits object
    void write(iofits & iof) const;

    ///////////////////////////////////////////
    ///  Print
    void print(ostream & os, const char * prefix="") const;

    ///////////////////////////////////////////
    ///  Get a rectangular region guaranteed to cover
    ///  the layer.  The resulting region will
    ///  have its edges aligned with the x and y axes
    ///  of the three_frame tf.  
    ///
    ///  If foreshortening is on, the projected
    ///  region is guaranteed to cover the optic
    ///
    ///  If the z axis of the three_frame is orthogonal
    ///  to the z axis of the layer, this function
    ///  throws an error
    rectangular_region get_covering_region(const three_frame & tf) const;

    ///////////////////////////////////////////
    ///  Apply the aperture to the wavefront
    virtual void transform(diffractive_wavefront<float> & wf) const;

    ///////////////////////////////////////////
    ///  Apply the aperture to the wavefront
    virtual void transform(diffractive_wavefront<double> & wf) const;

    ///////////////////////////////////////////
    ///  Apply the aperture to a geometric_wavefront
    //virtual void transform(geometric_wavefront & gwf) const;

    ///////////////////////////////////////////
    ///  Rotate layer by an angle - in radians
    void rotate_by_fft(double angle, bool window=true);

    ///////////////////////////////////////////
    ///  Get the pixel scale
    double get_pixel_scale() const {return(pixel_scale);};

    ///////////////////////////////////////////
    ///  Set the pixel scale
    void set_pixel_scale(double pscale) {pixel_scale = pscale;};

    ///////////////////////////////////////////
    ///  Get the axes
    vector<long> get_axes() const {return(this->axes);};

    ///////////////////////////////////////////
    ///  Set the axes - this destroys and
    ///  reallocates the underlying pixel array
    ///  if in_axes != this->axes
    void set_axes(const vector<long> & in_axes);

    ///////////////////////////////////////////
    ///  Get the wind vector
    three_vector get_wind_vector() const {return(wind_vector);};

    ///////////////////////////////////////////
    ///  Set the wind vector.
    void set_wind_vector(const three_vector & wvec) {wind_vector = wvec;};

  };

  template<class T>
  refractive_atmospheric_layer<T>::refractive_atmospheric_layer() {
    pixel_scale = 0;
  }

  template<class T>
  refractive_atmospheric_layer<T>::
  refractive_atmospheric_layer(
  	const refractive_atmospheric_layer<T> & ref_atm_layer){
    this->operator=(ref_atm_layer);
  }

  template<class T>
  refractive_atmospheric_layer<T>::
  refractive_atmospheric_layer(
  		const refractive_atmospheric_layer<T> & ref_atm_layer, 
		const rectangular_region & subregion) {

    rectangular_region arg_layer_region(ref_atm_layer, 
					ref_atm_layer.get_axes(), 
					ref_atm_layer.get_pixel_scale());

    if(!subregion.aligned(arg_layer_region)){
      cerr << "refractive_atmospheric_layer::refractive_atmospheric_layer error - "
	   << "regions not aligned\n";
      subregion.print(cerr, "this region ");
      arg_layer_region.print(cerr, "arg region  ");
      throw(string("refractive_atmospheric_layer::refractive_atmospheric_layer"));
    }
    
    if(!arg_layer_region.contains(subregion)){
      cerr << "refractive_atmospheric_layer::refractive_atmospheric_layer error - "
	   << "requested region not contained in layer\n";
      subregion.print(cerr, "subregion ");
      arg_layer_region.print(cerr, "arg region  "); 
      throw(string("refractive_atmospheric_layer::refractive_atmospheric_layer"));
    }
   
    // Enlarge subregion so that it is an even number of pixels
    // in size.
    rectangular_region integral_pixel_subregion(subregion,
    					ref_atm_layer.get_pixel_scale());

    // Find the subpixel shift between the two regions
    three_point arg_layer_center = arg_layer_region.get_center();
    three_point integral_pixel_subregion_center =
    				integral_pixel_subregion.get_center();

    three_vector delta_center =
    			integral_pixel_subregion_center - arg_layer_center;

    double delta_center_length = delta_center.length();
    three_vector subpixel_shift;
    if(delta_center_length>three_frame::precision) {
      subpixel_shift = 
	(fmod(delta_center_length, ref_atm_layer.get_pixel_scale()) / 
	 delta_center_length)*delta_center;
    }

    // zero out the subpixel shift if components are less than the limiting precision
    if(fabs(subpixel_shift.x(ref_atm_layer))<three_frame::precision)
      subpixel_shift = three_vector(0, subpixel_shift.y(ref_atm_layer),
      			0, ref_atm_layer);
    if(fabs(subpixel_shift.y(ref_atm_layer))<three_frame::precision)
      subpixel_shift = three_vector(subpixel_shift.x(ref_atm_layer),
      			0, 0, ref_atm_layer);

    if(optic::verbose_level)
      cout << "refractive_atmospheric_layer::refractive_atmospheric_layer - "
           << "subpixel_shift "
	   << subpixel_shift.x(ref_atm_layer) << "\t"
	   << subpixel_shift.y(ref_atm_layer) << endl;

    // Get the corners of the region that we will extract
    vector<three_point> integral_pixel_subregion_corners =
    			integral_pixel_subregion.get_corners();
    for(int i=0; i<4; i++)
      integral_pixel_subregion_corners[i] -= subpixel_shift;

    // Construct a frame of reference centered at the corner of the arg layer
    double pixscale = ref_atm_layer.get_pixel_scale();
    vector<long> arg_layer_axes = ref_atm_layer.get_axes();
    three_vector corner_vec(pixscale*arg_layer_axes[0]/2.0, 
			    pixscale*arg_layer_axes[1]/2.0, 
			    0, ref_atm_layer);
    three_translation ttrans(corner_vec);
    three_frame local_frame(ref_atm_layer);
    ttrans.transform(local_frame);

    // Find the limiting pixel coordinates of the subregion
    vector<long> coord_limits(4);
    coord_limits[0] = coord_limits[1] = 
      (long)(ceil(integral_pixel_subregion_corners[0].x(local_frame)/pixscale));
    coord_limits[2] = coord_limits[3] = 
      (long)(ceil(integral_pixel_subregion_corners[0].y(local_frame)/pixscale));
    long tx, ty;
    for(int i=1; i<4; i++){
      tx = (long)(ceil(integral_pixel_subregion_corners[i].x(local_frame)/pixscale));
      ty = (long)(ceil(integral_pixel_subregion_corners[i].y(local_frame)/pixscale));
      if(tx<coord_limits[0]) coord_limits[0] = tx;
      if(tx>coord_limits[1]) coord_limits[1] = tx;
      if(ty<coord_limits[2]) coord_limits[2] = ty;
      if(ty>coord_limits[3]) coord_limits[3] = ty;
    }

    if(optic::verbose_level)
      cout << "refractive_atmospheric_layer::refractive_atmospheric_layer - "
           << "extracting range "
	   << coord_limits[0] << " - " << coord_limits[1] << "  "
	   << coord_limits[2] << " - " << coord_limits[3] << endl;

    // extract the appropriate subregion and perform the shift by fft to account
    // for the non integral shift
    this->pixel_array<T>::operator=(pixel_array<T>(ref_atm_layer,coord_limits));
    this->shift_by_fft(subpixel_shift.x(local_frame),
    				subpixel_shift.y(local_frame));

    // fix up the three_frame of the layer so that its origin
    // is in the center of the new subregion
    this->three_point::operator=(subregion.get_center());

    // and finally the pixel scale
    this->pixel_scale = pixscale;    
  }

  template<class T>
  refractive_atmospheric_layer<T>::
  refractive_atmospheric_layer(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "refractive_atmospheric_layer::refractive_atmospheric_layer - "
	   << "error opening file " << filename << endl;
      throw(string("refractive_atmospheric_layer::refractive_atmospheric_layer"));
    }
    this->read(iof);
  } 

  template<class T>
  refractive_atmospheric_layer<T>::refractive_atmospheric_layer(const iofits & iof){
    this->read(iof);
  } 

  template<class T>
  refractive_atmospheric_layer<T>::
  refractive_atmospheric_layer(const refractive_atmospheric_layer & ref_atm_layer, 
			       vector<three_point> & corners) {
    pixel_scale = ref_atm_layer.pixel_scale;
    wind_vector = ref_atm_layer.wind_vector;
    this->three_frame::operator=(ref_atm_layer);
  }

  template<class T>
  refractive_atmospheric_layer<T>::
  refractive_atmospheric_layer(const pixel_array<T> & pixarr, 
			       double pixscale)
    :  pixel_array<T>(pixarr) {

    if(pixarr.get_axes().size()!=2){
      cerr << "refractive_atmospheric_layer::refractive_atmospheric_layer error - "
	   << "cannot allocate layer with dimension "
	   << pixarr.get_axes().size() << endl;
      throw(string("refractive_atmospheric_layer::refractive_atmospheric_layer"));
    }
    pixel_scale = pixscale;
  }

  template<class T>
  refractive_atmospheric_layer<T> & refractive_atmospheric_layer<T>::
  operator=(const refractive_atmospheric_layer<T> & ref_atm_layer){
    if(this==&ref_atm_layer)
      return(*this);

    pixel_scale = ref_atm_layer.pixel_scale;
    wind_vector = ref_atm_layer.wind_vector;
    this->plane_optic::operator=(ref_atm_layer);
    this->pixel_array<T>::operator=(ref_atm_layer);
    return(*this);
  }

  template<class T>
  refractive_atmospheric_layer<T>::
    refractive_atmospheric_layer(const power_spectrum * pspectrum,
				 const subharmonic_method & subm,
				 const vector<long> & axes, 
				 double pixscale){

    if(axes.size()!=2){
      cerr << "refractive_atmospheric_layer<T>::refractive_atmospheric_layer error - "
	   << "cannot construct refractive atmospheric layer with "
	   << axes.size() << "axes\n";
      throw(string("refractive_atmospheric_layer<T>::refractive_atmospheric_layer"));
    }

    for(int i=0; i<axes.size(); i++){
      if(axes[i]<=0){
	cerr << "refractive_atmospheric_layer<T>::"
	     << "refractive_atmospheric_layer error - "
	     << " axis " << i << " has value " << axes[i] 
	     << ", which is less than zero - can't make a "
	     << "refractive atmospheric layer with this dimension\n";
	throw(string("refractive_atmospheric_layer<T>::refractive_atmospheric_layer"));
      }
    }

    if(pixscale<=0){
      cerr << "refractive_atmospheric_layer<T>::refractive_atmospheric_layer error -\n"
	   << "invalid pixel scale " << pixscale << endl;
      throw(string("refractive_atmospheric_layer<T>::refractive_atmospheric_layer"));
    }

    // This is a temporary trick to force the
    // computation to always occur on axes suitable for
    // the subharmonic method
    vector<long> working_axes = axes;
    if(subm.intrinsic_dimensionality()!=-1){
      if(working_axes[0]%2!=subm.intrinsic_dimensionality()) working_axes[0]++;
      if(working_axes[1]%2!=subm.intrinsic_dimensionality()) working_axes[1]++;
    } 

    T * data;
    //try{data = new T[2*working_axes[1]*working_axes[0]];}
    //alloc_size sz(working_axes[0], working_axes[1], 2);
    try{
      //data = new T[sz];
      data = new T[2*working_axes[1]*working_axes[0]];
    }
    catch(...){
      cerr << "refractive_atmospheric_layer<T>::refractive_atmospheric_layer error - "
	   << "unable to allocate memory: " << 2*working_axes[1]*working_axes[0]*sizeof(T) << " bytes\n";
      throw;
    }

    initialize_frequency_array<T>(data, pspectrum, working_axes, pixscale, true, subm);

    fft_manager<T> fft_mgr;
    halfpixel_fft<T>(working_axes, data, fft_mgr);

    // Nominally we're losing half the effort here.
    // It would be better to try to do a real to complex
    // fft, but it is not so easy to get it right.
    for(int i=0; i<working_axes[1]; i++)
      for(int j=0; j<working_axes[0]; j++)
	data[i*working_axes[0]+j] = data[2*(i*working_axes[0]+j)];

    // Sort the data back into the original axes
    if(axes!=working_axes){
      for(int i=0; i<axes[1]; i++)
	for(int j=0; j<axes[0]; j++)
	  data[i*axes[0]+j] = data[i*working_axes[0]+j];
    }

    pixel_array<T> pixarr(axes, data);
    delete [] data; 
    this->operator=(refractive_atmospheric_layer<T>(pixarr, pixscale));
  }

  template<class T>
  void refractive_atmospheric_layer<T>::set_axes(const vector<long> & in_axes){
    if(in_axes.size()!=2){
      cerr << "refractive_atmospheric_layer::set_axes error - "
	   << "cannot set axes to dimension "
	   << in_axes.size() << " - must be 2 dimensional\n";
      throw(string("refractive_atmospheric_layer::set_axes"));
    }
    pixel_array<T>::set_axes(in_axes);
  }

  template<class T>
  void refractive_atmospheric_layer<T>::read(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "refractive_atmospheric_layer::read - "
	   << "error opening file " << filename << endl;
      throw(string("refractive_atmospheric_layer::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "refractive_atmospheric_layer::read - "
	   << "error reading "
	   << this->unique_name() << " from file "
	   << filename << endl;
      throw(string("refractive_atmospheric_layer::read"));
    }
  }

  template<class T>
  void refractive_atmospheric_layer<T>::read(const iofits & iof){

    if(!iof.key_exists("TYPE")){
      cerr << "refractive_atmospheric_layer::read error - "
	   << "unrecognized type of file\n";
      iof.print_header(cerr, "hdr dump");
      cerr << "hdu num " << iof.get_hdu_num() << " of "
	   << iof.get_num_hdus() << endl;
      throw(string("refractive_atmospheric_layer::read"));
    }
    string type, comment;
    iof.read_key("TYPE", type, comment);
    if(type!=this->unique_name()){
      cerr << "refractive_atmospheric_layer::read error - file of type " 
	   << type << " rather than of type "
	   << this->unique_name() << endl;
      throw(string("refractive_atmospheric_layer::read"));
    }

    iof.read_key("PIXSCALE", pixel_scale, comment);
    wind_vector.read(iof);
    this->plane_optic::read(iof);

    this->pixel_array<T>::read(iof);

    if(iof.get_hdu_num()<iof.get_num_hdus())
      iof.movrel_hdu(1);
  }

  template<class T>
  void refractive_atmospheric_layer<T>::write(const char * filename) const {
    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "refractive_atmospheric_layer::write - "
	   << "error opening file " << filename << endl;
      throw(string("refractive_atmospheric_layer::write"));
    }
    try{this->write(iof);}
    catch(...){
      cerr << "refractive_atmospheric_layer::write - "
	   << "error writing " 
	   << this->unique_name() << " to file "
	   << filename << endl;
      throw(string("refractive_atmospheric_layer::write"));
    }
  }

  template<class T>
  void refractive_atmospheric_layer<T>::write(iofits & iof) const {

    fits_header_data<double> fhd(this->get_axes());
    fhd.write(iof);
    string type = this->unique_name();
    string comment = "object type";
    iof.write_key("TYPE", type, comment);
    type = "pixel scale (meters)";
    iof.write_key("PIXSCALE", pixel_scale, type);
    wind_vector.write(iof);
    this->plane_optic::write(iof);
    this->pixel_array<T>::write(iof);
  }

  template<class T>
  void refractive_atmospheric_layer<T>::
  print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << this->unique_name()
       << "/" << "object type" << endl;
    fits_header_data<double> fhd(this->get_axes());
    fhd.print(os, prefix);
    os << prefix << "PIXSCALE   = " << setw(vlspc) << pixel_scale
       << "/" << "pixel scale (meters)" << endl;
    three_frame tf;
    wind_vector.print(os, tf, prefix);
    this->plane_optic::print(cout, prefix);
  }

  template<class T>
  void refractive_atmospheric_layer<T>::rotate_by_fft(double angle, bool window){
    this->pixel_array<T>::rotate_by_fft(angle, window);
    three_rotation trot(*this, this->z(), -angle);
    trot.transform(*this);
  }

  template<class T>
  rectangular_region
  refractive_atmospheric_layer<T>::get_covering_region(
  					const three_frame & tf) const {

    // This is the same code as
    // rectangular_region::get_covering_region - you could consolidate
    // it...

    if(fabs(dot_product(tf.z(), this->z()))<three_frame::precision){
      cerr << "refractive_atmospheric_layer::get_covering_region error - "
	   << "z axes of this aperture and three frame provided "
	   << "to this function are orthogonal\n";
      throw(string("refractive_atmospheric_layer::get_covering_region"));
    }

    vector<double> size(2, pixel_scale*this->axes[0]);
    size[1] = pixel_scale*this->axes[1];

    vector<three_point> tp(4);
    tp[0] = *this + .5*size[0]*this->x() + .5*size[1]*this->y();
    tp[1] = *this - .5*size[0]*this->x() + .5*size[1]*this->y();
    tp[2] = *this - .5*size[0]*this->x() - .5*size[1]*this->y();
    tp[3] = *this + .5*size[0]*this->x() - .5*size[1]*this->y();

    if(foreshortening){
      // project off the component along tf.z() to get
      // a three point in the XY plane of tf
      int aligned = 1;
      if(dot_product(this->z(), tf.z())<1) aligned = -1;
      for(int i=0; i<4; i++)
	tp[i] -= aligned*dot_product(tp[i]-*this, tf.z())*tf.z();
    } else {
      three_vector rotation_axis = cross_product(this->z(), tf.z());
      if(rotation_axis.length() > three_frame::precision){
	three_rotation trot(*this, rotation_axis, rotation_axis.length());
	for(int i=0; i<4; i++)
	  trot.transform(tp[i]);
      }
    }
    
    // Next we find the extremal coordinates of the points identified above
    // in the three_frame tf
    double xmin = tp[0].x(tf);
    double xmax = xmin;
    double ymin = tp[0].y(tf);
    double ymax = ymin;
    for(int i=1; i<4; i++){
      if(tp[i].x(tf)>xmax) xmax = tp[i].x(tf);
      if(tp[i].x(tf)<xmin) xmin = tp[i].x(tf);
      if(tp[i].y(tf)>ymax) ymax = tp[i].y(tf);
      if(tp[i].y(tf)<ymin) ymin = tp[i].y(tf);
    }
    
    // Finally we make a rectangular region from these limiting coordinates
    tp[0] = three_point(xmin, ymin, 0, tf);
    tp[1] = three_point(xmax, ymin, 0, tf);
    tp[2] = three_point(xmax, ymax, 0, tf);
    tp[3] = three_point(xmin, ymax, 0, tf);
    return(rectangular_region(tp));
    
  }

  template<class T>
  void refractive_atmospheric_layer<T>::transform(
  			diffractive_wavefront<float> & wf) const {
    try{this->private_transform(wf);}
    catch(...){
      cerr << "refractive_atmospheric_layer::transform error - "
	   << "error transforming a float instantiation of "
	   << "diffractive_wavefront\n";
      throw(string("refractive_atmospheric_layer::transform"));
    }
  }

  template<class T>
  void refractive_atmospheric_layer<T>::transform(
  			diffractive_wavefront<double> & wf) const {
    try{this->private_transform(wf);}
    catch(...){
      cerr << "refractive_atmospheric_layer::transform error - "
	   << "error transforming a float instantiation "
	   << "of diffractive_wavefront\n";
      throw(string("refractive_atmospheric_layer::transform"));
    }
  }

  template<class T> template<class U>
  void refractive_atmospheric_layer<T>::private_transform(
  			diffractive_wavefront<U> & wf) const {

    // This frame will represent the wavefront frame for the transformation.
    // If foreshortening is off, this frame is different than the true
    // wavefront frame
    three_frame wf_frame(wf);
    if(!foreshortening){
      three_vector rotation_axis = cross_product(this->z(), wf.z());
      if(rotation_axis.length()>three_frame::precision){
	three_rotation trot(wf, rotation_axis, rotation_axis.length());
	trot.transform(wf_frame);
      }
    }

    // If the transverse axes are aligned, we can call
    //   aligned_private_transform directly
    if(fabs(dot_product(wf_frame.x(), this->y()))<three_frame::precision && 
       fabs(dot_product(wf_frame.y(), this->x()))<three_frame::precision){
      if(optic::verbose_level)
	cerr << "refractive_atmospheric_layer::private_transform"
	     << " - calling aligned_private_transform directly...";
      try{this->aligned_private_transform(wf);}
      catch(...){
	cerr << "refractive_atmospheric_layer::private_transform error - "
	     << "could not perform transform\n";
	throw(string("refractive_atmospheric_layer::private_transform"));
      }
      if(optic::verbose_level) cout << "returning\n";
      return;      
    }

    // Otherwise, we have to construct a layer with the above property
    // by extracting the minimal part of the layer and rotating it to 
    // align with the wavefront axes.

    // Ensure that at least one basis vector of the wavefront frame
    // lies in the plane of the layer.  In this way we can rotate the
    // layer about its z axis so that it is aligned with this axis.
    // Then the wf pixels projected onto the plane of the layer will
    // be rectangles rather than parallelograms, and we can call the 
    // aligned_private_transform function.
    double dot_xz = dot_product(wf_frame.x(), this->z());
    double dot_yz = dot_product(wf_frame.y(), this->z());
    if(fabs(dot_xz)>three_frame::precision
    				&& fabs(dot_yz)>three_frame::precision){
      cerr << "refractive_atmospheric_layer::private_transform error - "
	   << "at least one of the transverse axes of the wavefront must lie "
	   << "in the plane of the layer.\n";
      wf_frame.three_frame::print(cerr, "wf three_frame ");
      this->three_frame::print(cerr, "layer three_frame ");
      throw(string("refractive_atmospheric_layer::private_transform"));
    }

    // Here we define the rotation angle as negative because
    // we are in essence doing a passive rotation - rotating
    // the pixel frame rather than the image itself
    double rotation_angle;
    if(fabs(dot_xz)<three_frame::precision)
      rotation_angle = -acos(dot_product(wf_frame.x(), this->x()));
    else
      rotation_angle = -acos(dot_product(wf_frame.y(), this->y()));

    // Compute the size of the wavefront projected onto the plane of the
    // layer

    // If foreshortening is on, we need to account for this.
    // If the z axes of the wf is not aligned with that of the layer,
    // project the wf region into the transverse plane of the layer.
    rectangular_region wf_region(*this, wf.get_axes(), wf.get_pixel_scale());
    if(this->foreshortening){
      wf_region = rectangular_region(wf_frame, wf.get_axes(),
      						wf.get_pixel_scale());
      wf_region = rectangular_region(wf_region, this->z(), false);
    } 

    vector<three_point> wf_projected_corners = wf_region.get_corners();

    // The wavefront projected into the plane of the layer is a
    // rectangular region whose edges are not necessarily aligned with
    // the axes of the layer.  We need to find the minimal rectangular
    // region that contains the projected wavefront and whose edges
    // are aligned with the axes of the layer.  To do so, we just 
    // get the extremal coordinates of the wf_projected_corners in 
    // the frame of the layer.
    double x_min = wf_projected_corners[0].x(*this);
    double x_max = x_min;
    double y_min = wf_projected_corners[0].y(*this);
    double y_max = y_min;
    double tmp;
    for(int i=1; i<4; i++){
      tmp = wf_projected_corners[i].x(*this);
      if(tmp<x_min) x_min = tmp;
      if(tmp>x_max) x_max = tmp;
      tmp = wf_projected_corners[i].y(*this);
      if(tmp<y_min) y_min = tmp;
      if(tmp>y_max) y_max = tmp;
    }

    // Get the corner coordinates of the region we'll extract.
    // Here we have to account for the wind, so we shift the 
    // layer_extraction corners by the wind vector.
    three_vector wind_vector = wf.get_timestamp()*this->get_wind_vector();

    vector<three_point> layer_extraction_corners(4);
    layer_extraction_corners[0] =
    		three_point(x_min, y_min, 0, *this) - wind_vector;
    layer_extraction_corners[1] =
    		three_point(x_max, y_min, 0, *this) - wind_vector;
    layer_extraction_corners[2] =
    		three_point(x_max, y_max, 0, *this) - wind_vector;
    layer_extraction_corners[3] =
    		three_point(x_min, y_max, 0, *this) - wind_vector;

    rectangular_region layer_extraction_region(layer_extraction_corners);

    // Finally, we extract this region, rotate it about its center,
    // and extract a final region that is of the same size as the 
    // projected wavefront.  We may then call aligned_private_transform 
    // using this layer.
    refractive_atmospheric_layer<T> xtrctd_layer(*this, layer_extraction_region);
    xtrctd_layer.rotate_by_fft(rotation_angle, false);
    xtrctd_layer.aligned_private_transform(wf);
  }

  /*
  template<class T> template<class U>
  void refractive_atmospheric_layer<T>::aligned_private_transform(
  				diffractive_wavefront<U> & wf) const {

    // First, ensure that the center of this wavefront lies in the
    // plane of the layer
    if(operator!=(static_cast<three_point>(wf),
    				static_cast<three_point>(*this))){
      three_vector origin_offset = 
	static_cast<three_point>(wf) - static_cast<three_point>(*this);
    
      if(fabs(dot_product(origin_offset,
      			this->three_frame::z()))>three_frame::precision){
	cerr << "refractive_atmospheric_layer::aligned_private_transform error - "
	     << "diffractive_wavefront center not in transverse plane of layer\n";
	wf.three_point::print(cerr, *this, "wavefront center in frame of layer ");
	throw(string("refractive_atmospheric_layer::aligned_private_transform"));
      }
    }
    
    // Verify that the transverse axes are aligned. 
    // This check is valid regardless of the foreshortening
    // status.
    if(fabs(dot_product(wf.x(), this->y()))>three_frame::precision &&
       fabs(dot_product(wf.y(), this->x()))<three_frame::precision){
      cerr << "refractive_atmospheric_layer::aligned_private_transform error - "
	   << "the wavefront transverse axes are not aligned "
	   << "with those of the layer\n";
      throw(string("refractive_atmospheric_layer::aligned_private_transform"));
    }

    three_frame wind_blown_frame(*this);
    wind_blown_frame += wf.get_timestamp()*this->get_wind_vector();

    // Verify that the wavefront is fully contained by the layer.
    rectangular_region *layer_region, *wf_region, *modified_wf_region;

    try{
      layer_region = new rectangular_region(wind_blown_frame,
					    this->get_axes(), this->get_pixel_scale());
    } catch(...) {
      cerr << "refractive_atmospheric_layer::aligned_private_transform error - "
	   << "could not get layer region\n";
      throw(string("refractive_atmospheric_layer::aligned_private_transform"));
    }

    try{
      wf_region = new rectangular_region(wf, wf.get_axes(), wf.get_pixel_scale());
    } catch(...) {
      cerr << "refractive_atmospheric_layer::aligned_private_transform error - "
	   << "could not get layer region\n";
      throw(string("refractive_atmospheric_layer::aligned_private_transform"));
    }

    if(this->foreshortening){
      modified_wf_region = new rectangular_region(*wf_region, this->z(), false);
    } else {
      three_frame tmp_frame;
      if(dot_product(wf.z(), this->z())>0)
	tmp_frame = three_frame(wf, this->x(), this->y(), this->z());
      else 
	tmp_frame = three_frame(wf, this->x(), this->y(), -1*this->z());
      modified_wf_region = new rectangular_region(tmp_frame,
      				wf.get_axes(), wf.get_pixel_scale());
    }

    try{
      if(!layer_region->contains(*modified_wf_region)){
	throw(string(""));
      }
    } catch(...) {
      cerr << "refractive_atmospheric_layer::aligned_private_transform error - "
	   << "layer does not contain wavefront\n";
      layer_region->print(cerr, "layer region ");
      modified_wf_region->print(cerr, "wfrnt region ");
      throw(string("refractive_atmospheric_layer::aligned_private_transform"));
    }

    delete layer_region, wf_region, modified_wf_region;

    // the half pixel definitions - so that
    // when you have even axes the pixel
    // centroids lie on half pixel values.
    vector<long> wf_axes = wf.get_axes();
    double wf_x_halfpix=0, wf_y_halfpix=0;
    double layer_x_halfpix=0, layer_y_halfpix=0;
    int wf_x_extrapix=1, wf_y_extrapix=1;
    int layer_x_extrapix=1, layer_y_extrapix=1;
    if(wf_axes[1]%2==0){
      wf_x_halfpix = .5;
      wf_x_extrapix = 0;
    }
    if(wf_axes[0]%2==0){
      wf_y_halfpix = .5;
      wf_y_extrapix = 0;
    }
    if(this->axes[1]%2==0){
      layer_x_halfpix = .5;
      layer_x_extrapix = 0;
    }
    if(this->axes[0]%2==0){
      layer_y_halfpix = .5;
      layer_y_extrapix = 0;
    }
 
    // Finally, we are ready to run through the arrays adding the layer
    // to the wavefront.  If foreshortening is on the pixel scale of the
    // wavefront may differ in the x and y direction, due to the
    // projection effect.  Here we define some quantities to specify the
    // wf pixel scales
    double wf_x_pixscale = wf.get_pixel_scale();
    double wf_y_pixscale = wf.get_pixel_scale();
    double dp;
    if(this->get_foreshortening()){
      dp = dot_product(wf.x(), this->x());
      wf_x_pixscale /= (dp*dp);
      dp = dot_product(wf.y(), this->y());
      wf_y_pixscale /= (dp*dp);
      if(optic::verbose_level)
	cout << "refractive_atmospheric_layer::aligned_private_transform -\n"
	     << "\twavefront pixel scale " << wf.get_pixel_scale() << endl
	     << "\tprojected x pixel scale " << wf_x_pixscale << endl
	     << "\tprojected y pixel scale " << wf_y_pixscale << endl;
    }

    // Convert the wavefront to amplitude and phase.  In this way
    // we can add in the phase from the layer OPD and it won't wrap
    this->wavefront_amp_phase_conversion(wf);
    U * wfdata = this->get_wavefront_data(wf);

    // Next, for each pixel in the wavefront, we are going to average
    // over the pixels in the layer that are contained by the wavefront
    // pixel.  If there are layer pixels that partially overlap the
    // wavefront pixel, we'll downweight them by their areal overlap.

    double wf_wavelength = wf.get_wavelength();
    double wf_x_pixel_coord, wf_y_pixel_coord;
    double tmp_wf_x_pixel_coord, tmp_wf_y_pixel_coord;
    double max_wf_x_pixel_coord, max_wf_y_pixel_coord;
    double area, refractive_optical_path;
    three_vector origin_offset = wf - wind_blown_frame;
    bool is_real_imag = this->is_real_imag_storage(wf);
    bool is_interleaved = this->is_interleaved_storage(wf);
    double tmp, cp, sp;
    int index, wf_nelem = wf_axes[0]*wf_axes[1],
      layer_nelem = this->axes[0]*this->axes[1];

    for(int i=-wf_axes[1]/2; i<wf_axes[1]/2+wf_x_extrapix; i++){
      for(int j=-wf_axes[0]/2; j<wf_axes[0]/2+wf_y_extrapix; j++){

      
	// These are the coordinates of the lower corner of the wf pixel
	// in the wind blown frame, measured in units of the layer pixel scale.
	wf_x_pixel_coord = ((j+wf_x_halfpix-.5)*wf_x_pixscale + 
			    origin_offset.x(wind_blown_frame))/this->pixel_scale;
	wf_y_pixel_coord = ((i+wf_y_halfpix-.5)*wf_y_pixscale + 
			    origin_offset.y(wind_blown_frame))/this->pixel_scale;


	// These are the coordinates of the upper corner of the wf pixel
	// in the wind blown frame, measured in units of the layer pixel scale.
	max_wf_x_pixel_coord =
			wf_x_pixel_coord + wf_x_pixscale/this->pixel_scale;
	max_wf_y_pixel_coord =
			wf_y_pixel_coord + wf_y_pixscale/this->pixel_scale;
	
	// These are the coordinates of the upper corner of the layer
	// pixel in the wind blown frame, measured in units of the layer
	// pixel scale.  
	//
	// If the lower corner of the wf pixel happens to lie on a layer
	// pixel boundary, then ceil() won't change the values and so we
	// make sure to explicitly increment these.
	tmp_wf_x_pixel_coord = ceil(wf_x_pixel_coord);
	if(fabs(tmp_wf_x_pixel_coord-wf_x_pixel_coord)<three_frame::precision)
	  tmp_wf_x_pixel_coord++;
	tmp_wf_y_pixel_coord = ceil(wf_y_pixel_coord);
	if(fabs(tmp_wf_y_pixel_coord-wf_y_pixel_coord)<three_frame::precision)
	  tmp_wf_y_pixel_coord++;

	refractive_optical_path = 0;
	area = 0;

	while(wf_x_pixel_coord < max_wf_x_pixel_coord){

	  if(tmp_wf_x_pixel_coord > max_wf_x_pixel_coord)
	    tmp_wf_x_pixel_coord = max_wf_x_pixel_coord;
	  
	  while(wf_y_pixel_coord < max_wf_y_pixel_coord){
	  
	    if(tmp_wf_y_pixel_coord > max_wf_y_pixel_coord)
	      tmp_wf_y_pixel_coord = max_wf_y_pixel_coord;
	    
	    index = (int)((ceil(tmp_wf_y_pixel_coord)-1
	    		+layer_y_extrapix+this->axes[1]/2)*this->axes[0] +
			ceil(tmp_wf_x_pixel_coord)-1
			+layer_x_extrapix+this->axes[0]/2);
	    
	    if(index<0 || index>=layer_nelem){
	      cerr << "\nrefractive_atmospheric_layer::aligned_private_transform"
	           << "- impending doom\n";
	      cerr << i << "\t" << j << "\t" << index << "\t" << layer_nelem
	           << endl;
	      cerr << wf_x_pixel_coord << "\t" << tmp_wf_x_pixel_coord 
		   << "\t" << ceil(tmp_wf_x_pixel_coord) << "\t" 
		   << max_wf_x_pixel_coord << endl;
	      cerr << wf_y_pixel_coord << "\t" << tmp_wf_y_pixel_coord 
		   << "\t" << ceil(tmp_wf_y_pixel_coord) << "\t" 
		   << max_wf_y_pixel_coord << endl;
	      cerr << this->axes[0] << "\t" << this->axes[1] << endl;
	      cerr << wf_axes[0] << "\t" << wf_axes[1] << endl;
	      throw(string("xxx"));
	    }

	    refractive_optical_path += 
	      this->pixeldata[index]*
	      (tmp_wf_x_pixel_coord - wf_x_pixel_coord)
	      *(tmp_wf_y_pixel_coord - wf_y_pixel_coord);

	    area += (tmp_wf_x_pixel_coord - wf_x_pixel_coord)
	    	*(tmp_wf_y_pixel_coord - wf_y_pixel_coord);

	    wf_y_pixel_coord = tmp_wf_y_pixel_coord;
	    tmp_wf_y_pixel_coord++;
	  }
	  wf_x_pixel_coord = tmp_wf_x_pixel_coord;
	  tmp_wf_x_pixel_coord++;
	}

	// add in the phase 
	tmp = refractive_optical_path*2*M_PI/wf_wavelength/area;
	if(is_interleaved){
	  index = (i+wf_axes[1]/2)*wf_axes[0] + j + wf_axes[0]/2;
	  wfdata[2*index+1] += tmp;
	} else if(!is_interleaved){
	  index = (i+wf_axes[1]/2)*wf_axes[0] + j + wf_axes[0]/2;
	  wfdata[index+wf_nelem] += tmp;
	}	
      }
    }
  }
  */

  template<class T> template<class U>
  void refractive_atmospheric_layer<T>::
    aligned_private_transform(diffractive_wavefront<U> & wf) const {

    // First, ensure that the center of this wavefront lies in the
    // plane of the layer
    if(operator!=(static_cast<three_point>(wf),
		  static_cast<three_point>(*this))){
      three_vector origin_offset = 
	static_cast<three_point>(wf) - static_cast<three_point>(*this);
    
      if(fabs(dot_product(origin_offset,
			  this->three_frame::z()))>three_frame::precision){
	cerr << "refractive_atmospheric_layer::aligned_private_transform error - "
	     << "diffractive_wavefront center not in transverse plane of layer\n";
	wf.three_point::print(cerr, *this, "wavefront center in frame of layer ");
	throw(string("refractive_atmospheric_layer::aligned_private_transform"));
      }
    }
    
    // Verify that the transverse axes are aligned. 
    // This check is valid regardless of the foreshortening
    // status.
    if(fabs(dot_product(wf.x(), this->y()))>three_frame::precision &&
       fabs(dot_product(wf.y(), this->x()))<three_frame::precision){
      cerr << "refractive_atmospheric_layer::aligned_private_transform error - "
	   << "the wavefront transverse axes are not aligned "
	   << "with those of the layer\n";
      throw(string("refractive_atmospheric_layer::aligned_private_transform"));
    }

    three_frame wind_blown_frame(*this);
    wind_blown_frame += wf.get_timestamp()*this->get_wind_vector();

    //(wf.get_timestamp()*this->get_wind_vector()).print(cerr, "wind vector ");

    // Verify that the wavefront is fully contained by the layer.
    rectangular_region layer_region, wf_region, modified_wf_region;

    try{
      layer_region = rectangular_region(wind_blown_frame,
					this->get_axes(), 
					this->get_pixel_scale());
    } catch(...) {
      cerr << "refractive_atmospheric_layer::aligned_private_transform error - "
	   << "could not get layer region\n";
      throw(string("refractive_atmospheric_layer::aligned_private_transform"));
    }

    try{
      wf_region = rectangular_region(wf, 
				     wf.get_axes(), 
				     wf.get_pixel_scale());
    } catch(...) {
      cerr << "refractive_atmospheric_layer::aligned_private_transform error - "
	   << "could not get layer region\n";
      throw(string("refractive_atmospheric_layer::aligned_private_transform"));
    }

    if(this->foreshortening){
      modified_wf_region = rectangular_region(wf_region, 
					      this->z(), 
					      false);
    } else {
      three_frame tmp_frame;
      if(dot_product(wf.z(), this->z())>0)
	tmp_frame = three_frame(wf, this->x(), this->y(), this->z());
      else 
	tmp_frame = three_frame(wf, this->x(), this->y(), -1*this->z());
      modified_wf_region = rectangular_region(tmp_frame,
					      wf.get_axes(), 
					      wf.get_pixel_scale());
    }

    if(!layer_region.contains(modified_wf_region)){
      cerr << "refractive_atmospheric_layer::aligned_private_transform error - "
	   << "layer does not contain wavefront\n";
      layer_region.print(cerr, "layer region ");
      modified_wf_region.print(cerr, "wfrnt region ");
      throw(string("refractive_atmospheric_layer::aligned_private_transform"));
    }

    // the half pixel definitions - so that
    // when you have even axes the pixel
    // centroids lie on half pixel values.
    vector<long> wf_axes = wf.get_axes();
    double wf_x_halfpix=0, wf_y_halfpix=0;
    double layer_x_halfpix=0, layer_y_halfpix=0;
    int wf_x_extrapix=1, wf_y_extrapix=1;
    int layer_x_extrapix=1, layer_y_extrapix=1;
    if(wf_axes[1]%2==0){
      wf_x_halfpix = .5;
      wf_x_extrapix = 0;
    }
    if(wf_axes[0]%2==0){
      wf_y_halfpix = .5;
      wf_y_extrapix = 0;
    }
    if(this->axes[1]%2==0){
      layer_x_halfpix = .5;
      layer_x_extrapix = 0;
    }
    if(this->axes[0]%2==0){
      layer_y_halfpix = .5;
      layer_y_extrapix = 0;
    }

    double x_sign = 1;
    if(dot_product(wf.x(),this->x())<0)
      x_sign = -1;
    double y_sign = 1;
    if(dot_product(wf.y(),this->y())<0)
      y_sign = -1;

    //cerr << "x sign " << x_sign << " y sign " << y_sign << endl;

    // Finally, we are ready to run through the arrays adding the layer
    // to the wavefront.  If foreshortening is on the pixel scale of the
    // wavefront may differ in the x and y direction, due to the
    // projection effect.  Here we define some quantities to specify the
    // wf pixel scales
    double wf_x_pixscale = wf.get_pixel_scale();
    double wf_y_pixscale = wf.get_pixel_scale();
    double dp;
    if(this->get_foreshortening()){
      dp = dot_product(wf.x(), this->x());
      wf_x_pixscale /= (dp*dp);
      dp = dot_product(wf.y(), this->y());
      wf_y_pixscale /= (dp*dp);
      if(optic::verbose_level)
	cout << "refractive_atmospheric_layer::aligned_private_transform -\n"
	     << "\twavefront pixel scale " << wf.get_pixel_scale() << endl
	     << "\tprojected x pixel scale " << wf_x_pixscale << endl
	     << "\tprojected y pixel scale " << wf_y_pixscale << endl;
    }

    // Convert the wavefront to amplitude and phase.  In this way
    // we can add in the phase from the layer OPD and it won't wrap
    this->wavefront_amp_phase_conversion(wf);
    U * wfdata = this->get_wavefront_data(wf);

    // Next, for each pixel in the wavefront, we are going to average
    // over the pixels in the layer that are contained by the wavefront
    // pixel.  If there are layer pixels that partially overlap the
    // wavefront pixel, we'll downweight them by their areal overlap.

    double wf_wavelength = wf.get_wavelength();
    double lower_wf_x_pixel_coord, lower_wf_y_pixel_coord;
    double upper_wf_x_pixel_coord, upper_wf_y_pixel_coord;
    long index_wf_x_pixel_coord, index_wf_y_pixel_coord;

    double min_wf_x_pixel_coord, min_wf_y_pixel_coord;
    double max_wf_x_pixel_coord, max_wf_y_pixel_coord;

    double area, refractive_optical_path;

    three_vector origin_offset = wf - wind_blown_frame;
    bool is_real_imag = this->is_real_imag_storage(wf);
    bool is_interleaved = this->is_interleaved_storage(wf);
    double tmp;
    int index, wf_nelem = wf_axes[0]*wf_axes[1],
      layer_nelem = this->axes[0]*this->axes[1];

    /*
    origin_offset.print(cerr, "origin offset ");
    cerr << "x offset "
	 << origin_offset.x(wind_blown_frame)
	 << endl;
    cerr << "y offset "
	 << origin_offset.y(wind_blown_frame)
	 << endl;
    cerr << endl << endl;
    */

    for(int i=-wf_axes[1]/2; i<wf_axes[1]/2+wf_x_extrapix; i++){
      for(int j=-wf_axes[0]/2; j<wf_axes[0]/2+wf_y_extrapix; j++){
      
	// These are the coordinates of the lower corner of the wf pixel
	// in the wind blown frame, measured in units of the layer pixel scale.
	/*
	min_wf_x_pixel_coord = ((j+wf_x_halfpix-.5)*wf_x_pixscale + 
				origin_offset.x(wind_blown_frame))/this->pixel_scale;
	min_wf_y_pixel_coord = ((i+wf_y_halfpix-.5)*wf_y_pixscale + 
				origin_offset.y(wind_blown_frame))/this->pixel_scale;
	*/
	min_wf_x_pixel_coord = (((j+wf_x_halfpix)*x_sign - .5)*wf_x_pixscale + 
				origin_offset.x(wind_blown_frame))/this->pixel_scale;
	min_wf_y_pixel_coord = (((i+wf_y_halfpix)*y_sign - .5)*wf_y_pixscale + 
				origin_offset.y(wind_blown_frame))/this->pixel_scale;
	/*
	if(i==0 && j==0){
	  cerr << "xxx - min wf x coord " << min_wf_x_pixel_coord << endl;
	  cerr << "xxx - min wf y coord " << min_wf_y_pixel_coord << endl;
	}
	*/

	// These are the coordinates of the upper corner of the wf pixel
	// in the wind blown frame, measured in units of the layer pixel scale.
	max_wf_x_pixel_coord = min_wf_x_pixel_coord + wf_x_pixscale/this->pixel_scale;
	max_wf_y_pixel_coord = min_wf_y_pixel_coord + wf_y_pixscale/this->pixel_scale;

	/*
	if(i==0 && j==0){
	  cerr << "xxx - max wf x coord " << max_wf_x_pixel_coord << endl;
	  cerr << "xxx - max wf y coord " << max_wf_y_pixel_coord << endl;
	}
	*/

	refractive_optical_path = 0;

	index_wf_x_pixel_coord = (long)floor(min_wf_x_pixel_coord);
	while(index_wf_x_pixel_coord < max_wf_x_pixel_coord){
	  index_wf_y_pixel_coord = (long)floor(min_wf_y_pixel_coord);
	  while(index_wf_y_pixel_coord < max_wf_y_pixel_coord){

	    // This is the index into the layer array
	    index = 
	      (index_wf_y_pixel_coord+layer_y_extrapix+this->axes[1]/2)*this->axes[0] +
	      index_wf_x_pixel_coord+layer_x_extrapix+this->axes[0]/2;
	    
	    lower_wf_x_pixel_coord = 
	      index_wf_x_pixel_coord < min_wf_x_pixel_coord ? 
	      min_wf_x_pixel_coord : index_wf_x_pixel_coord;

	    lower_wf_y_pixel_coord = 
	      index_wf_y_pixel_coord < min_wf_y_pixel_coord ? 
	      min_wf_y_pixel_coord : index_wf_y_pixel_coord;

	    upper_wf_x_pixel_coord = 
	      index_wf_x_pixel_coord + 1 < max_wf_x_pixel_coord ? 
	      index_wf_x_pixel_coord + 1 : max_wf_x_pixel_coord;

	    upper_wf_y_pixel_coord = 
	      index_wf_y_pixel_coord + 1 < max_wf_y_pixel_coord ? 
	      index_wf_y_pixel_coord + 1 : max_wf_y_pixel_coord;

	    /*
	    if(i==0 && j==0){
	      cerr << "xxx - lower wf coords " 
		   << lower_wf_x_pixel_coord 
		   << ","
		   << lower_wf_y_pixel_coord 
		   << endl;
	      cerr << "xxx - index wf coords " 
		   << index_wf_x_pixel_coord 
		   << ","
		   << index_wf_y_pixel_coord 
		   << endl;
	      cerr << "xxx - upper wf coords " 
		   << upper_wf_x_pixel_coord 
		   << ","
		   << upper_wf_y_pixel_coord 
		   << endl;
	    }
	    */
	
	    area = (upper_wf_x_pixel_coord - lower_wf_x_pixel_coord)
	      *(upper_wf_y_pixel_coord - lower_wf_y_pixel_coord);

	    /*
	    if(i==0 && j==0){
	      cerr << "xxx - area " << area << endl;
	    }
	    */

	    if(index<0 || index>=layer_nelem || area<=0){
	      cerr << "\nrefractive_atmospheric_layer::aligned_private_transform"
	           << " - error transforming wavefront\n";
	      cerr << "indices "
		   << i 
		   << "\t" << j 
		   << "\t" << index 
		   << "\t" << layer_nelem
	           << endl;
	      cerr << "x coords "
		   << "\tlower " << lower_wf_x_pixel_coord 
		   << "\tupper " << upper_wf_x_pixel_coord 
		   << "\tmin   " << min_wf_x_pixel_coord 
		   << "\tmax   " << max_wf_x_pixel_coord 
		   << endl;
	      cerr << "y coords "
		   << "\tlower " << lower_wf_y_pixel_coord 
		   << "\tupper " << upper_wf_y_pixel_coord 
		   << "\tmin   " << min_wf_y_pixel_coord 
		   << "\tmax   " << max_wf_y_pixel_coord 
		   << endl;
	      cerr << "layer axes "
		   << this->axes[0] 
		   << "\t" << this->axes[1] 
		   << endl;
	      cerr << "wavefront axes "
		   << wf_axes[0] 
		   << "\t" << wf_axes[1] 
		   << endl;
	      cerr << "area "
		   << area
		   << endl;
	      throw(string("refractive_atmospheric_layer::aligned_private_transform"));
	    }

	    refractive_optical_path += 
	      this->pixeldata[index]*area;

	    /*
	    if(i==0 && j==0){
	      cerr << "xxx - refractive optical path " 
		   << this->pixeldata[index] 
		   << "\t"
		   << this->pixeldata[index]*area
		   << "\t"
		   << refractive_optical_path 
		   << endl;
	    }
	    */

	    index_wf_y_pixel_coord++;
	  }
	  index_wf_x_pixel_coord++;
	}

	// add in the phase 
	tmp = refractive_optical_path*2*M_PI/wf_wavelength;

	/*
	if(i==0 && j==0){
	  cerr << "xxx - layer x index " << index_wf_y_pixel_coord+layer_y_extrapix+this->axes[1]/2 << endl;
	  cerr << "xxx - layer y index " << index_wf_x_pixel_coord+layer_x_extrapix+this->axes[0]/2 << endl;
	  cerr << "xxx - refractive optical path " << refractive_optical_path << endl;
	  cerr << "xxx - phase difference " << tmp << endl << endl;
	}
	*/

	if(is_interleaved){
	  index = (i+wf_axes[1]/2)*wf_axes[0] + j + wf_axes[0]/2;
	  wfdata[2*index+1] += tmp;
	} else if(!is_interleaved){
	  index = (i+wf_axes[1]/2)*wf_axes[0] + j + wf_axes[0]/2;
	  wfdata[index+wf_nelem] += tmp;
	}	
      }
    }
  }
}

#endif
