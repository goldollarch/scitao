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

#include <iostream>
#include "fits_factory.h"
#include "region_base.h"
#include "computational_geometry.h"
#include "diffractive_wavefront.h"
#include "lenslet_array.h"
#include "sim_utils.h"

using namespace std;

namespace Arroyo {

  lenslet_array_base * lenslet_array_base::lenslet_array_base_factory(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...) {
      cerr << "lenslet_array_base::lenslet_array_base_factory - "
	   << "error opening file " << filename << endl;
      throw(string("lenslet_array_base::lenslet_array_base_factory"));
    }
    return(lenslet_array_base::lenslet_array_base_factory(iof));
  }

  lenslet_array_base * lenslet_array_base::lenslet_array_base_factory(const iofits & iof){
    AO_sim_base * aosb = AO_sim_base::AO_sim_base_factory(iof);
    lenslet_array_base * labase = dynamic_cast<lenslet_array_base *>(aosb);
    if(labase==NULL)
      throw(string("lenslet_array_base::lenslet_array_base_factory"));
     return(labase);    
  }

  lenslet_array_base::lenslet_array_base(const lenslet_array_base & labase) {
    this->operator=(labase);
  }

  lenslet_array_base & lenslet_array_base::operator=(const lenslet_array_base & labase) {
    if(this==&labase) 
      return(*this);
    this->plane_optic::operator=(labase);
    return(*this);
  }

  void lenslet_array_base::write(iofits & iof) const {
    fits_header_data<char> tmphdr;
    tmphdr.write(iof);
    this->plane_optic::write(iof);
  }

  namespace factory_register {
     const fits_keyval_set & get_square_lenslet_array_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "square lenslet array"));
      return *fkvs;
    }
    
    AO_sim_base * create_square_lenslet_array(const iofits & iof) {
      return new square_lenslet_array(iof);
    }
  }

  const bool square_lenslet_array::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_square_lenslet_array_keyval_set(), 
				      factory_register::create_square_lenslet_array);

  square_lenslet_array::square_lenslet_array(const square_lenslet_array & sq_lens_arr){
    this->operator=(sq_lens_arr);
  }

  square_lenslet_array::square_lenslet_array(const char * filename){
    axes.resize(2);
    this->read(filename);
  }

  square_lenslet_array::square_lenslet_array(const iofits & iof){
    axes.resize(2);
    this->read(iof);
  }

  square_lenslet_array::square_lenslet_array(vector<long> array_axes,
					     double flength,
					     double lnslt_pitch,
					     long pix_per_lenslet,
					     long pix_per_xform){
    
    if(array_axes.size()!=2){
      cerr << "square_lenslet_array::square_lenslet_array error - "
	   << "cannot construct lenslet array with " 
	   << array_axes.size() << " - must have 2 axes\n";
      throw(string("square_lenslet_array::square_lenslet_array"));
    }

    if(array_axes[0]<=0 || array_axes[1]<=0){
      cerr << "square_lenslet_array::square_lenslet_array error - "
	   << "cannot construct lenslet array with axes " 
	   << array_axes[0] << "x" << array_axes[1] << endl;
      throw(string("square_lenslet_array::square_lenslet_array"));
    }

    if(flength<=0){
      cerr << "square_lenslet_array::square_lenslet_array error - "
	   << "focal length " << flength 
	   << " supplied to constructor is not positive\n";
      throw(string("square_lenslet_array::square_lenslet_array"));
    }
    if(lnslt_pitch<=0){
      cerr << "square_lenslet_array::square_lenslet_array error - "
	   << "lenslet pitch " << lnslt_pitch
	   << " supplied to constructor is not positive\n";
      throw(string("square_lenslet_array::square_lenslet_array"));
    }
    if(pix_per_lenslet<=0){
      cerr << "square_lenslet_array::square_lenslet_array error - "
	   << "pixels per lenslet " << pix_per_lenslet
	   << " supplied to constructor is not positive\n";
      throw(string("square_lenslet_array::square_lenslet_array"));
    }
    if(pix_per_xform<=0){
      cerr << "square_lenslet_array::square_lenslet_array error - "
	   << "pixels per xform " << pix_per_xform
	   << " supplied to constructor is not positive\n";
      throw(string("square_lenslet_array::square_lenslet_array"));
    }

    this->set_foreshortening(false);

    axes = array_axes;
    focal_length = flength;
    lenslet_pitch = lnslt_pitch;
    final_wavefront_propagation_distance = flength;
    final_wavefront_pixels_per_lenslet = pix_per_lenslet;
    final_wavefront_pixels_per_transform = pix_per_xform;
  }

  square_lenslet_array & square_lenslet_array::operator=(const square_lenslet_array & sq_lens_arr){
    if(this==&sq_lens_arr)
      return(*this);

    axes = sq_lens_arr.axes;
    focal_length = sq_lens_arr.focal_length;
    lenslet_pitch = sq_lens_arr.lenslet_pitch;
    final_wavefront_propagation_distance = 
      sq_lens_arr.final_wavefront_propagation_distance;
    final_wavefront_pixels_per_lenslet = sq_lens_arr.final_wavefront_pixels_per_lenslet;
    final_wavefront_pixels_per_transform = sq_lens_arr.final_wavefront_pixels_per_transform;

    this->lenslet_array_base::operator=(sq_lens_arr);
    return(*this);
  }

  void square_lenslet_array::set_final_wavefront_pixels_per_lenslet(long pix_per_lnslt) {
    if(pix_per_lnslt>final_wavefront_pixels_per_transform){
      cerr << "square_lenslet_array::set_final_wavefront_pixels_per_lenslet error - "
	   << "pixels per lenslet " << pix_per_lnslt 
	   << " cannot exceed pixels per transform " << final_wavefront_pixels_per_transform
	   << endl;
      throw(string("square_lenslet_array::set_final_wavefront_pixels_per_lenslet"));
    }
    final_wavefront_pixels_per_lenslet = pix_per_lnslt;
  }

  void square_lenslet_array::set_final_wavefront_pixels_per_transform(long pix_per_xform) {
    if(pix_per_xform<final_wavefront_pixels_per_lenslet){
      cerr << "square_lenslet_array::set_final_wavefront_pixels_per_transform error - "
	   << "pixels per transform " << pix_per_xform 
	   << " cannot be less than pixels per lenslet " << final_wavefront_pixels_per_lenslet
	   << endl;
      throw(string("square_lenslet_array::set_final_wavefront_pixels_per_transform"));
    }
    final_wavefront_pixels_per_transform = pix_per_xform;
  }

  void square_lenslet_array::read(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "square_lenslet_array::read - "
	   << "error opening file " << filename << endl;
      throw(string("square_lenslet_array::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "square_lenslet_array::read - "
	   << "error reading "
	   << this->unique_name() << " from file "
	   << filename << endl;
      throw(string("square_lenslet_array::read"));
    }
  }

  void square_lenslet_array::read(const iofits & iof) {
    if(!iof.key_exists("TYPE")){
      cerr << "square_lenslet_array::read error - "
	   << "unrecognized type of file\n";
      throw(string("square_lenslet_array::read"));
    }
    string type, comment;
    iof.read_key("TYPE", type, comment);
    if(type!=this->unique_name()){
      cerr << "square_lenslet_array::read error - file of type " 
	   << type << " rather than of type "
	   << this->unique_name() << endl;
      throw(string("square_lenslet_array::read"));
    }

    this->plane_optic::read(iof);

    iof.read_key("FOCLNGTH", focal_length, comment);
    iof.read_key("LNSPITCH", lenslet_pitch, comment);
    iof.read_key("LNSDIMX", axes[0], comment);
    iof.read_key("LNSDIMY", axes[1], comment);
    iof.read_key("WFTPROPD", final_wavefront_propagation_distance, comment);
    iof.read_key("WFTPPLNS", final_wavefront_pixels_per_lenslet, comment);
    iof.read_key("WFTPPXFM", final_wavefront_pixels_per_transform, comment);
    
    // leave iof pointing to the next header, if it exists
    if(iof.get_hdu_num()!=iof.get_num_hdus())
      iof.movrel_hdu(1);
  }

  void square_lenslet_array::write(const char * filename) const {

    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "square_lenslet_array::write - "
	   << "error opening file " << filename << endl;
      throw(string("square_lenslet_array::write"));
    }

    try{this->write(iof);}
    catch(...){
      cerr << "square_lenslet_array::write - "
	   << "error writing "
	   << this->unique_name() << " to file " 
	   << filename << endl;
      throw(string("square_lenslet_array::write"));
    }
  }

  void square_lenslet_array::write(iofits & iof) const {
    this->lenslet_array_base::write(iof);
    string type = this->unique_name();
    string comment = "object type";
    iof.write_key("TYPE", type, comment);

    iof.write_key("FOCLNGTH", focal_length, "lenslet focal length (meters)");
    iof.write_key("LNSPITCH", lenslet_pitch, "lenslet pitch (meters)");
    iof.write_key("LNSDIMX", axes[0], "number of lenslets along x axis");
    iof.write_key("LNSDIMY", axes[1], "number of lenslets along y axis");
    iof.write_key("WFTPROPD", final_wavefront_propagation_distance, "final wavefront propagation distance (meters)");
    iof.write_key("WFTPPLNS", final_wavefront_pixels_per_lenslet, "final wavefront pixels per lenslet");
    iof.write_key("WFTPPXFM", final_wavefront_pixels_per_transform, "final wavefront pixels per transform");
  }

  void square_lenslet_array::print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << this->unique_name()
       << "/" << "object type" << endl;
    this->plane_optic::print(os, prefix);
    os << prefix << "FOCLNGTH   = " << setw(vlspc) << focal_length
       << "/" << "lenslet focal length (meters)" << endl;
    os << prefix << "LNSPITCH   = " << setw(vlspc) << lenslet_pitch
       << "/" << "lenslet pitch (meters)" << endl;
    os << prefix << "LNSDIMX    = " << setw(vlspc) << axes[0]
       << "/" << "number of lenslets along x axis " << endl;
    os << prefix << "LNSDIMY    = " << setw(vlspc) << axes[1]
       << "/" << "number of lenslets along y axis " << endl;
    os << prefix << "WFTPROPD   = " << setw(vlspc) << final_wavefront_propagation_distance
       << "/" << "final wavefront propagation distance (meters) " << endl;
    os << prefix << "WFTPPLNS   = " << setw(vlspc) << final_wavefront_pixels_per_lenslet
       << "/" << "final wavefront pixels per lenslet" << endl;
    os << prefix << "WFTPPXFM   = " << setw(vlspc) << final_wavefront_pixels_per_transform
       << "/" << "final wavefront pixels per transform" << endl;
  }

  rectangular_region square_lenslet_array::get_covering_region(const three_frame & tf) const {

    if(fabs(dot_product(tf.z(), this->z()))<three_frame::precision){
      cerr << "square_lenslet_array::get_covering_region error - "
	   << "z axes of this aperture and three frame provided to this function are orthogonal\n";
      throw(string("square_lenslet_array::get_covering_region"));
    }

    vector<three_point> tp(4);
    tp[0] = *this + .5*lenslet_pitch*axes[0]*this->x() + .5*lenslet_pitch*axes[1]*this->y();
    tp[1] = *this - .5*lenslet_pitch*axes[0]*this->x() + .5*lenslet_pitch*axes[1]*this->y();
    tp[2] = *this - .5*lenslet_pitch*axes[0]*this->x() - .5*lenslet_pitch*axes[1]*this->y();
    tp[3] = *this + .5*lenslet_pitch*axes[0]*this->x() - .5*lenslet_pitch*axes[1]*this->y();

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

  void square_lenslet_array::transform(diffractive_wavefront<float> & wf) const {
    try{this->private_transform(wf);}
    catch(...){
      cerr << "square_lenslet_array::transform error - "
	   << "error transforming a float instantiation of diffractive_wavefront\n";
      throw(string("square_lenslet_array::transform"));
    }
  }

  void square_lenslet_array::transform(diffractive_wavefront<double> & wf) const {
    try{this->private_transform(wf);}
    catch(...){
      cerr << "square_lenslet_array::transform error - "
	   << "error transforming a double instantiation of diffractive_wavefront\n";
      throw(string("square_lenslet_array::transform"));
    }
  }

  template<class T>
  void square_lenslet_array::private_transform(diffractive_wavefront<T> & wf) const { 

    three_vector origin_offset, dx, dy;
    try{
      this->plane_optic::get_projected_wavefront_pixel_spacing(wf, origin_offset, dx, dy);
    } catch(...) {
      cerr << "square_lenslet_array::private_transform error - could not get projected wavefront pixel spacing\n";
      throw(string("square_lenslet_array::private_transform"));
    }

    if(wf.get_pixel_scale()>=lenslet_pitch){
      cerr << "square_lenslet_array::private_transform error - "
	   << "lenslet pitch " << lenslet_pitch
	   << " is smaller than the wavefront pixel scale "
	   << wf.get_pixel_scale() << endl;
      throw(string("square_lenslet_array::private_transform"));
    }

    // Check whether there is any overlap between the 
    // aperture and the diffractive_wavefront.  If not, zero the wavefront
    if(origin_offset.x(*this)>lenslet_pitch*axes[0] ||
       origin_offset.y(*this)>lenslet_pitch*axes[1]){
      wf*=complex<double>(0,0);
      return;
    }

    bool real_imag_storage = is_real_imag_storage(wf);
    bool interleaved_storage = is_interleaved_storage(wf);

    // Project lenslet edges into the plane of the wavefront
    three_vector x_lenslet_vector = lenslet_pitch*this->x();
    three_vector y_lenslet_vector = lenslet_pitch*this->y();

    if(foreshortening){
      x_lenslet_vector = parallel_projection(x_lenslet_vector, wf.z(), this->z());
      y_lenslet_vector = parallel_projection(y_lenslet_vector, wf.z(), this->z());
    }

    if(optic::verbose_level){
      x_lenslet_vector.print(cout, "x lenslet vector ");
      y_lenslet_vector.print(cout, "y lenslet vector ");
    }

    // Determine size of array required to hold enough
    // wavefront pixels to completely cover a single lenslet
    vector<long> lenslet_covering_axes(2, (long)(ceil((x_lenslet_vector+y_lenslet_vector).x(*this)/wf.get_pixel_scale())+1));

    // Allocate memory for small array and for an array to 
    // hold the wf data overlapping the lenslet
    T * tmp_wf_data;
    //int tmp_nelem = lenslet_covering_axes[0]*lenslet_covering_axes[1];
    alloc_size sz(lenslet_covering_axes[0], lenslet_covering_axes[1], 2);
    //try{tmp_wf_data = new T[2*tmp_nelem];}
    try{
      tmp_wf_data = new T[sz];
    }
    catch(...){
      cerr << "square_lenslet_array::private_transform error - "
	   << "could not allocate memory for temporary wavefront data array\n";
      throw;
    }

    // Allocate memory for small array and for an array to 
    // hold the propagated wavefront from each lenslet
    T * tmp_wf_data_psf;
    //tmp_nelem = final_wavefront_pixels_per_transform*final_wavefront_pixels_per_transform;
    alloc_size sz2(final_wavefront_pixels_per_transform,
    		final_wavefront_pixels_per_transform, 2);
    //try{tmp_wf_data_psf = new T[2*tmp_nelem];}
    try{
      tmp_wf_data_psf = new T[sz2];
    }
    catch(...){
      cerr << "square_lenslet_array::private_transform error - "
	   << "could not allocate memory for temporary wavefront data array\n";
      throw;
    }

    // Allocate memory for an array to hold the entire
    // propagated wavefront
    T * final_wf_data;
    int pad = final_wavefront_pixels_per_transform > final_wavefront_pixels_per_lenslet ? 
      final_wavefront_pixels_per_transform - final_wavefront_pixels_per_lenslet : 0;
    vector<long> final_wf_axes(2, final_wavefront_pixels_per_lenslet*axes[0]+pad);
    final_wf_axes[1] = final_wavefront_pixels_per_lenslet*axes[1]+pad;
    //int final_nelem = final_wf_axes[0]*final_wf_axes[1];
    alloc_size sz3(final_wf_axes[0], final_wf_axes[1], 2);
    //try{final_wf_data = new T[2*final_nelem];} 
    try{
      final_wf_data = new T[sz3];
    } 
    catch(...){
      cerr << "square_lenslet_array::private_transform error - "
	   << "could not allocate memory for final wavefront data array\n";
      throw;
    }
    unsigned long final_nelem = final_wf_axes[0]*final_wf_axes[1];
    for(int i=0; i<final_nelem; i++) {
      final_wf_data[2*i] = final_wf_data[2*i+1] = 0;
    }

    // The lenslet halfpixel information
    double lenslet_x_halfpix=0, lenslet_y_halfpix=0;
    int lenslet_x_extrapix=1, lenslet_y_extrapix=1;
    if(this->axes[1]%2==0){
      lenslet_x_halfpix = .5;
      lenslet_x_extrapix = 0;
    }
    if(this->axes[0]%2==0){
      lenslet_y_halfpix = .5;
      lenslet_y_extrapix = 0;
    }

    // The wf halfpixel information
    double wf_x_halfpix=0, wf_y_halfpix=0;
    int wf_x_extrapix=1, wf_y_extrapix=1;
    if(wf.get_axes()[1]%2==0){
      wf_x_halfpix = .5;
      wf_x_extrapix = 0;
    }
    if(wf.get_axes()[0]%2==0){
      wf_y_halfpix = .5;
      wf_y_extrapix = 0;
    }

    // The covering halfpixel information
    double covering_x_halfpix=0, covering_y_halfpix=0;
    int covering_x_extrapix=1, covering_y_extrapix=1;
    if(lenslet_covering_axes[1]%2==0){
      covering_x_halfpix = .5;
      covering_x_extrapix = 0;
    }
    if(lenslet_covering_axes[0]%2==0){
      covering_y_halfpix = .5;
      covering_y_extrapix = 0;
    }

    // The final wf halfpixel information
    double final_wf_x_halfpix=0, final_wf_y_halfpix=0;
    int final_wf_x_extrapix=1, final_wf_y_extrapix=1;
    if(final_wf_axes[1]%2==0){
      final_wf_x_halfpix = .5;
      final_wf_x_extrapix = 0;
    }
    if(final_wf_axes[0]%2==0){
      final_wf_y_halfpix = .5;
      final_wf_y_extrapix = 0;
    }

    double final_wf_pixel_scale = lenslet_pitch/(double)final_wavefront_pixels_per_lenslet;

    // definitions to account for halfpixel shifts for even initial array dimensions
    double initial_nyquist_pixel_scale =
      wf.get_wavelength() * final_wavefront_propagation_distance / 
      (double)lenslet_covering_axes[0] / wf.get_pixel_scale();
    double initial_xslope = 
      final_wavefront_pixels_per_transform%2==1 ? 0 : 
      M_PI*final_wf_pixel_scale/initial_nyquist_pixel_scale/(double)lenslet_covering_axes[1];
    double initial_yslope = 
      final_wavefront_pixels_per_transform%2==1 ? 0 : 
      M_PI*final_wf_pixel_scale/initial_nyquist_pixel_scale/(double)lenslet_covering_axes[0];

    // definitions to account for halfpixel shifts for even final array dimensions
    double final_nyquist_pixel_scale = 
      wf.get_wavelength() * final_wavefront_propagation_distance / 
      (double)final_wavefront_pixels_per_transform / final_wf_pixel_scale;
    double final_xslope = lenslet_covering_axes[1]%2==1 ? 0 : 
      M_PI*wf.get_pixel_scale()/final_nyquist_pixel_scale/(double)final_wavefront_pixels_per_transform;
    double final_yslope = lenslet_covering_axes[0]%2==1 ? 0 : 
      M_PI*wf.get_pixel_scale()/final_nyquist_pixel_scale/(double)final_wavefront_pixels_per_transform;

    if(optic::verbose_level){
      cout << "square_lenslet_array::private_transform " 
	   << "lenslet covering axes " << lenslet_covering_axes[0] << " x " << lenslet_covering_axes[1] << endl;
      cout << "square_lenslet_array::private_transform " 
	   << "final wf pix per xform " << final_wavefront_pixels_per_transform << endl;
      cout << "square_lenslet_array::private_transform " 
	   << "final pixel scale " << final_wf_pixel_scale << endl;
      cout << "square_lenslet_array::private_transform " 
	   << "initial nyquist " << initial_nyquist_pixel_scale << endl;
      cout << "square_lenslet_array::private_transform " 
	   << "slope " << initial_xslope << "\t" << initial_yslope << endl;
    }
    
    // Information for the GR transform
    // from far field propagator: double sampling_factor = pixel_scale*final_pixel_scale / wavelength / distance;
    double gr_sampling_factor = wf.get_pixel_scale() * final_wf_pixel_scale / 
      final_wavefront_propagation_distance / wf.get_wavelength();
    vector<long> final_transform_axes(2, final_wavefront_pixels_per_transform);

    // other definitions
    three_vector lenslet_center, wf_pixel_center;
    three_point wf_pixel_coords;
    three_frame lenslet_centered_frame;
    int min_wf_x_index, min_wf_y_index;
    int max_wf_x_index, max_wf_y_index;
    T * orig_wf_data = get_wavefront_data(wf);
    vector<long> orig_wf_axes = wf.get_axes();
    double lenslet_weight, wf_amp, wf_phase;
    long index, wf_offset_nelem = wf.get_axes()[0]*wf.get_axes()[1];

    three_vector tmps = .5*(dy+dx);
    three_vector tmpd = .5*(dy-dx);
    vector<three_point> pixel_vertices(4), lenslet_vertices(4), intersection_vertices;

    for(int i=-this->axes[1]/2; i<this->axes[1]/2+lenslet_x_extrapix; i++){
      for(int j=-this->axes[0]/2; j<this->axes[0]/2+lenslet_y_extrapix; j++){

	// Check to see if any wf pixels lie in this lenslet.  If not,
	// continue.
	lenslet_center = three_vector((i+lenslet_x_halfpix)*lenslet_pitch, 
				      (j+lenslet_y_halfpix)*lenslet_pitch, 
				      0, *this);

	if(foreshortening)
	  lenslet_center = parallel_projection(lenslet_center, wf.z(), this->z());
	
	lenslet_centered_frame = *this;
	lenslet_centered_frame += lenslet_center;

	lenslet_vertices[0] = lenslet_centered_frame + .5*lenslet_centered_frame.x() + .5*lenslet_centered_frame.y();
	lenslet_vertices[1] = lenslet_centered_frame - .5*lenslet_centered_frame.x() + .5*lenslet_centered_frame.y();
	lenslet_vertices[2] = lenslet_centered_frame - .5*lenslet_centered_frame.x() - .5*lenslet_centered_frame.y();
	lenslet_vertices[3] = lenslet_centered_frame + .5*lenslet_centered_frame.x() - .5*lenslet_centered_frame.y();

	// Find the wf pixel coordinates of the corner of the lenslet covering array
	min_wf_x_index = (long)((lenslet_centered_frame-.5*x_lenslet_vector-.5*y_lenslet_vector).x(wf)/wf.get_pixel_scale()+wf.get_axes()[1]/2);
	min_wf_y_index = (long)((lenslet_centered_frame-.5*x_lenslet_vector-.5*y_lenslet_vector).y(wf)/wf.get_pixel_scale()+wf.get_axes()[0]/2);

	// Transfer the appropriate wavefront data into the temporary array,
	// taking into account the effect of the lenslet
	for(int k=0; k<lenslet_covering_axes[1]; k++){
	  for(int l=0; l<lenslet_covering_axes[0]; l++){
	    if((k+min_wf_x_index>=orig_wf_axes[1]) || (k+min_wf_x_index<0) ||
	       (l+min_wf_y_index>=orig_wf_axes[0]) || (l+min_wf_y_index<0)) {
	      tmp_wf_data[2*(k*lenslet_covering_axes[0]+l)] = 
		tmp_wf_data[2*(k*lenslet_covering_axes[0]+l)+1] = 0;
	    } else {
	      // the wavefront pixel coordinates in the plane of the lenslet array
	      // INCORRECT
	      wf_pixel_center = wf.get_pixel_scale()*((min_wf_x_index+k+wf_x_halfpix-wf.get_axes()[1]/2)*wf.x() + 
						      (min_wf_y_index+l+wf_y_halfpix-wf.get_axes()[0]/2)*wf.y());
	      if(foreshortening)
		wf_pixel_center = parallel_projection(wf_pixel_center, this->z(), wf.z());

	      wf_pixel_coords = static_cast<three_point>(wf) + wf_pixel_center;

	      pixel_vertices[0] = wf_pixel_coords + tmps;
	      pixel_vertices[1] = wf_pixel_coords + tmpd;
	      pixel_vertices[2] = wf_pixel_coords - tmps;
	      pixel_vertices[3] = wf_pixel_coords - tmpd;
	      
	      intersection_vertices = get_convex_polygon_intersection(pixel_vertices, lenslet_vertices);

	      lenslet_weight = 0;
	      if(intersection_vertices.size()>2)
		lenslet_weight = get_area_of_polygon(intersection_vertices);

	      index = (min_wf_x_index+k)*wf.get_axes()[0]+min_wf_y_index+l;

	      if(real_imag_storage){
		if(interleaved_storage){
		  wf_amp = sqrt(orig_wf_data[2*index]*orig_wf_data[2*index]+
				orig_wf_data[2*index+1]*orig_wf_data[2*index+1]);
		  wf_phase = atan2(orig_wf_data[2*index]*orig_wf_data[2*index+1],
				   orig_wf_data[2*index]*orig_wf_data[2*index]);
		} else {
		  wf_amp = sqrt(orig_wf_data[index]*orig_wf_data[index]+
				orig_wf_data[index+wf_offset_nelem]*orig_wf_data[index+wf_offset_nelem]);
		  wf_phase = atan2(orig_wf_data[index]*orig_wf_data[index+wf_offset_nelem],
				   orig_wf_data[index]*orig_wf_data[index]);
		}
	      } else {
		if(interleaved_storage){
		  wf_amp = orig_wf_data[2*index];
		  wf_phase = orig_wf_data[2*index+1];
		} else {
		  wf_amp = orig_wf_data[index];
		  wf_phase = orig_wf_data[index+wf_offset_nelem];
		}
	      }

	      // Here we add in the phase slope that will set the transformed
	      // array to be centered on the offsets appropriate for this function,
	      // as was done in diffractive_wavefront::far_field_propagator

	      tmp_wf_data[2*(k*lenslet_covering_axes[0]+l)] = 
		lenslet_weight * wf_amp * cos(wf_phase + 
					      initial_xslope*(k+covering_x_halfpix) +
					      initial_yslope*(l+covering_y_halfpix));
	      tmp_wf_data[2*(k*lenslet_covering_axes[0]+l)+1] = 
		lenslet_weight * wf_amp * sin(wf_phase + 
					      initial_xslope*(k+covering_x_halfpix) +
					      initial_yslope*(l+covering_y_halfpix));

	    }
	  }
	}

	// Perform the Goertzel Reinsch transform and add the data into
	// the array containing the final wavefront
	goertzel_reinsch_transform(lenslet_covering_axes, final_transform_axes, gr_sampling_factor, tmp_wf_data, tmp_wf_data_psf);

	min_wf_x_index = (long)((i+axes[1]/2)*final_wavefront_pixels_per_lenslet);
	min_wf_y_index = (long)((j+axes[0]/2)*final_wavefront_pixels_per_lenslet);
	for(int k=0; k<final_transform_axes[1]; k++){
	  for(int l=0; l<final_transform_axes[0]; l++){
	    // Here we subtract off the quadratic phase and any phase slope 
	    // that arises from even initial array dimensions
	    index = (min_wf_x_index+k)*final_wf_axes[0]+min_wf_y_index+l;
	    final_wf_data[2*index] += tmp_wf_data_psf[2*(k*final_transform_axes[0]+l)];
	    final_wf_data[2*index+1] += tmp_wf_data_psf[2*(k*final_transform_axes[0]+l)+1];
	  }
	}
      }
    }
 
    delete [] tmp_wf_data;
    delete [] tmp_wf_data_psf;

    // Set up the final wavefront header
    diffractive_wavefront_header<T> final_wfh(wf);
    final_wfh.set_axes(final_wf_axes);
    final_wfh.set_pixel_scale(final_wf_pixel_scale);
    final_wfh += final_wavefront_propagation_distance*final_wfh.z();

    wf = diffractive_wavefront<T>(final_wfh, final_wf_data, true, true);
    delete [] final_wf_data;
  }
}
