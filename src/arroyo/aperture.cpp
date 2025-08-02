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
//#include "geometric_wavefront.h"
#include "aperture.h"
#include "sim_utils.h"

using namespace std;

namespace Arroyo {

  aperture * aperture::aperture_factory(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...) {
      cerr << "aperture::aperture_factory - "
	   << "error opening file " << filename << endl;
      throw(string("aperture::aperture_factory"));
    }
    return(aperture::aperture_factory(iof));
  }

  aperture * aperture::aperture_factory(const iofits & iof){
    AO_sim_base * aosb = AO_sim_base::AO_sim_base_factory(iof);
    aperture * ap = dynamic_cast<aperture *>(aosb);
    if(ap==NULL)
      throw(string("aperture::aperture_factory"));
     return(ap);    
  }

  aperture::aperture(const aperture & ap) {
    this->operator=(ap);
  }

  aperture & aperture::operator=(const aperture & ap) {
    if(this==&ap) 
      return(*this);
    this->plane_optic::operator=(ap);
    areal_weighting = ap.areal_weighting;
    return(*this);
  }

  void aperture::read(const iofits & iof){ 
    string comment;
    iof.read_key("AREALWTG", areal_weighting, comment);
    this->plane_optic::read(iof);
  }
 
  void aperture::write(iofits & iof) const {
    fits_header_data<char> tmphdr;
    tmphdr.write(iof);
    iof.write_key("AREALWTG", areal_weighting, "weight edge pixels by area of overlap");
    this->plane_optic::write(iof);
  }

  void aperture::print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    this->plane_optic::print(os, prefix);
    os << prefix << "AREALWTG   = " << setw(vlspc) << areal_weighting
       << "/" << "weight diffractive wavefront edge pixels by area of overlap" << endl;
  }

  namespace factory_register {
    const fits_keyval_set & get_circular_aperture_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "circular aperture"));
      return *fkvs;
    }
    
    AO_sim_base * create_circular_aperture(const iofits & iof) {
      return new circular_aperture(iof);
    }

  } 

  const bool circular_aperture::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_circular_aperture_keyval_set(), 
				      factory_register::create_circular_aperture);


  circular_aperture::circular_aperture(){
    diameter = 0;
  }

  circular_aperture::circular_aperture(const circular_aperture & circ_ap){
    this->operator=(circ_ap);
  }

  circular_aperture::circular_aperture(const char * filename){
    this->read(filename);
  }

  circular_aperture::circular_aperture(const iofits & iof){
    this->read(iof);
  }

  circular_aperture::circular_aperture(double in_diameter){
    if(in_diameter < 0){
      cerr << "circular_aperture::circular_aperture error - "
	   << "diameter " << diameter 
	   << " supplied to constructor is less than zero\n";
      throw(string("circular_aperture::circular_aperture"));
    }

    diameter = in_diameter;
  }

  circular_aperture & circular_aperture::operator=(const circular_aperture & circ_ap){
    if(this==&circ_ap)
      return(*this);
    diameter = circ_ap.diameter;
    this->aperture::operator=(circ_ap);
    return(*this);
  }

  void circular_aperture::read(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "circular_aperture::read - "
	   << "error opening file " << filename << endl;
      throw(string("circular_aperture::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "circular_aperture::read - "
	   << "error reading "
	   << this->unique_name() << " from file "
	   << filename << endl;
      throw(string("circular_aperture::read"));
    }
  }

  void circular_aperture::read(const iofits & iof) {
    if(!iof.key_exists("TYPE")){
      cerr << "circular_aperture::read error - "
	   << "unrecognized type of file\n";
      throw(string("circular_aperture::read"));
    }
    string type, comment;
    iof.read_key("TYPE", type, comment);
    if(type!=this->unique_name()){
      cerr << "circular_aperture::read error - file of type " 
	   << type << " rather than of type "
	   << this->unique_name() << endl;
      throw(string("circular_aperture::read"));
    }
    this->aperture::read(iof);
    iof.read_key("DIAMETER", diameter, comment);
    // leave iof pointing to the next header, if it exists
    if(iof.get_hdu_num()!=iof.get_num_hdus())
      iof.movrel_hdu(1);
  }

  void circular_aperture::write(const char * filename) const {

    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "circular_aperture::write - "
	   << "error opening file " << filename << endl;
      throw(string("circular_aperture::write"));
    }

    try{this->write(iof);}
    catch(...){
      cerr << "circular_aperture::write - "
	   << "error writing "
	   << this->unique_name() << " to file " 
	   << filename << endl;
      throw(string("circular_aperture::write"));
    }
  }

  void circular_aperture::write(iofits & iof) const {
    this->aperture::write(iof);
    string type = this->unique_name();
    string comment = "object type";
    iof.write_key("TYPE", type, comment);

    iof.write_key("DIAMETER", diameter, "aperture diameter (meters)");
  }

  void circular_aperture::print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << this->unique_name()
       << "/" << "object type" << endl;
    this->aperture::print(os, prefix);
    os << prefix << "DIAMETER   = " << setw(vlspc) << diameter
       << "/" << "aperture diameter (meters)" << endl;
  }

  rectangular_region circular_aperture::get_covering_region(const three_frame & tf) const {

    if(fabs(dot_product(tf.z(), this->z()))<three_frame::precision){
      cerr << "circular_aperture::get_covering_region error - "
	   << "z axes of this aperture and three frame provided to this function are orthogonal\n";
      throw(string("circular_aperture::get_covering_region"));
    }

    /*
    if(foreshortening && cross_product(this->z(), tf.z()).length()>three_frame::precision){

      // In this case the circle projects into an ellipse in the tf
      // frame, with major axis orthogonal to the plane containing
      // tf.z() and this->z().  First we find the three points corresponding
      // to the projected semimajor and semiminor axes
      vector<three_point> limiting_points(4);
      three_vector tmp = cross_product(this->z(), tf.z());
      tmp = tmp*(1/tmp.length());
      limiting_points[0] = *this + .5*diameter*tmp;
      limiting_points[2] = *this - .5*diameter*tmp;
      double projected_diameter = fabs(diameter/dot_product(this->z(), tf.z()));
      cerr << " projected diameter " << setprecision(15) << projected_diameter << endl;
      tmp = cross_product(tf.z(), tmp);
      tmp = tmp*(1/tmp.length());      
      limiting_points[1] = *this + .5*projected_diameter*tmp;
      limiting_points[3] = *this - .5*projected_diameter*tmp;

      for(int i=0; i<4; i++)
	limiting_points[i].print(cerr, "limpt ");

      // Next we find the extremal coordinates of the points identified above
      // in the three_frame tf
      double xmin = limiting_points[0].x(tf);
      double xmax = xmin;
      double ymin = limiting_points[0].y(tf);
      double ymax = ymin;
      for(int i=1; i<4; i++){
	if(limiting_points[i].x(tf)>xmax) xmax = limiting_points[i].x(tf);
	if(limiting_points[i].x(tf)<xmin) xmin = limiting_points[i].x(tf);
	if(limiting_points[i].y(tf)>ymax) ymax = limiting_points[i].y(tf);
	if(limiting_points[i].y(tf)<ymin) ymin = limiting_points[i].y(tf);
      }

      // Finally we make a rectangular region from these limiting coordinates
      limiting_points[0] = three_point(xmin, ymin, 0, tf);
      limiting_points[1] = three_point(xmax, ymin, 0, tf);
      limiting_points[2] = three_point(xmax, ymax, 0, tf);
      limiting_points[3] = three_point(xmin, ymax, 0, tf);
      return(rectangular_region(limiting_points));

    } else {

      // Here we just make a rectangular region with sides of length
      // equal to the diameter
      vector<three_point> corners(4);
      corners[0] = *this + .5*diameter*tf.x() + .5*diameter*tf.y();
      corners[1] = *this - .5*diameter*tf.x() + .5*diameter*tf.y();
      corners[2] = *this - .5*diameter*tf.x() - .5*diameter*tf.y();
      corners[3] = *this + .5*diameter*tf.x() - .5*diameter*tf.y();
      return(rectangular_region(corners));
    }
    */
    

    double projected_radius;

    if(foreshortening && cross_product(this->z(), tf.z()).length()>three_frame::precision)
      projected_radius = fabs(.5*diameter/dot_product(this->z(), tf.z()));
    else 
      projected_radius = .5*diameter;

    three_point intersection_tp = this->get_point_of_intersection(tf, tf.z());

    three_vector offset_vector = *this - intersection_tp;
    double xmax = fabs(offset_vector.x(*this)) + projected_radius;
    double ymax = fabs(offset_vector.y(*this)) + projected_radius;

    if(foreshortening && cross_product(this->z(), tf.z()).length()>three_frame::precision){
      xmax /= dot_product(this->z(), tf.z());
      ymax /= dot_product(this->z(), tf.z());
    }

    vector<three_point> limiting_points(4);
    limiting_points[0] = intersection_tp + xmax*tf.x() + ymax*tf.y();
    limiting_points[1] = intersection_tp + xmax*tf.x() - ymax*tf.y();
    limiting_points[2] = intersection_tp - xmax*tf.x() - ymax*tf.y();
    limiting_points[3] = intersection_tp - xmax*tf.x() + ymax*tf.y();
    return(rectangular_region(limiting_points));

  }

  void circular_aperture::transform(diffractive_wavefront<float> & wf) const {
    try{this->private_transform(wf);}
    catch(...){
      cerr << "circular_aperture::transform error - "
	   << "error transforming a float instantiation of diffractive_wavefront\n";
      throw(string("circular_aperture::transform"));
    }
  }

  void circular_aperture::transform(diffractive_wavefront<double> & wf) const {
    try{this->private_transform(wf);}
    catch(...){
      cerr << "circular_aperture::transform error - "
	   << "error transforming a double instantiation of diffractive_wavefront\n";
      throw(string("circular_aperture::transform"));
    }
  }

  /*
  void circular_aperture::transform(geometric_ray & gray) const {
    if(gray.is_lost()) return;
    three_vector odiff = static_cast<three_point>(*this) - static_cast<three_point>(gray);
    if(dot_product(odiff, this->z())!=0){
      cerr << "circular_aperture::transform error - "
	   << "geometric ray does not lie in plane of aperture\n";
      throw(string("circular_aperture::transform"));
    }
    if(odiff.length_squared() > .25*diameter*diameter)
      gray.set_lost(true);
  }

  void circular_aperture::transform(geometric_wavefront & gwf) const {
    long nrays = gwf.get_number_of_rays();
    geometric_ray * gray = this->get_wavefront_data(gwf);
    for(int i=0; i<nrays; i++)
      this->transform(gray[i]);
  }
  */

  template<class T>
  void circular_aperture::private_transform(diffractive_wavefront<T> & wf) const { 
    
    // perform the aperture check
    three_vector origin_offset, dx, dy;
    this->get_projected_wavefront_pixel_spacing(wf, origin_offset, dx, dy);

    // Check whether there is any overlap between the 
    // aperture and the diffractive_wavefront.  If not, zero the wavefront
    if(origin_offset.length()>diameter){
      wf*=complex<double>(0,0);
      return;
    }

    // The halfpixel information
    vector<long> wf_axes = wf.get_axes();
    double x_halfpix=0, y_halfpix=0;
    int x_extrapix=1, y_extrapix=1;
    if(wf_axes[1]%2==0){
      x_halfpix = .5;
      x_extrapix = 0;
    }
    if(wf_axes[0]%2==0){
      y_halfpix = .5;
      y_extrapix = 0;
    }

    bool real_imag_storage = is_real_imag_storage(wf);
    bool interleaved_storage = is_interleaved_storage(wf);
    T * wfdata = get_wavefront_data(wf);

    int index;
    int nelem = wf_axes[0]*wf_axes[1];
    three_point center_pixel_coords;
    double val;
    three_vector tmps = .5*(dy+dx);
    three_vector tmpd = .5*(dy-dx);
    vector<three_point> pixel_vertices(4);
    double norm = cross_product(dx, dy).length();
    double buffer = tmps.length() > tmpd.length() ? 2*tmps.length() : 2*tmpd.length();
    double upper_limit = (.5*this->diameter+buffer)*(.5*this->diameter+buffer);
    double lower_limit = (.5*this->diameter-buffer)*(.5*this->diameter-buffer);
    double center_distance_squared;
    double radius_squared = .25*diameter*diameter;

    for(int i=-wf_axes[1]/2; i<wf_axes[1]/2+x_extrapix; i++){
      for(int j=-wf_axes[0]/2; j<wf_axes[0]/2+y_extrapix; j++){
	index = (i+wf_axes[1]/2)*wf_axes[0]+j+wf_axes[0]/2;

	center_pixel_coords = static_cast<const three_point>(*this) + 
	  (origin_offset + (i+x_halfpix)*dx + (j+y_halfpix)*dy);

	center_distance_squared = (center_pixel_coords - *this).length_squared();

	if(areal_weighting){
	  if(center_distance_squared>upper_limit) 
	    val = 0;
	  else if(center_distance_squared<lower_limit) 
	    val = 1;
	  else {
	    pixel_vertices[0] = center_pixel_coords + tmps;
	    pixel_vertices[1] = center_pixel_coords + tmpd;
	    pixel_vertices[2] = center_pixel_coords - tmps;
	    pixel_vertices[3] = center_pixel_coords - tmpd;
	    
	    val = convex_polygon_overlap(pixel_vertices)/norm;
	  }
	} else {
	  val = 1;
	  if(center_distance_squared > radius_squared)
	    val = 0;
	}

	if(val!=1){
	  if(interleaved_storage){
	    if(val==0){
	      wfdata[2*index] = wfdata[2*index+1] = 0;
	      continue;
	    }
	    wfdata[2*index] *= val;
	    if(real_imag_storage)
	      wfdata[2*index+1] *= val;
	  } else {
	    if(val==0){
	      wfdata[2*index] = wfdata[2*index+1] = 0;
	      continue;
	    }
	    wfdata[index] *= val;
	    if(real_imag_storage)
	      wfdata[index+nelem] *= val;
	  }
	}
      }
    }
  }

  double circular_aperture::convex_polygon_overlap(const vector<three_point> & polygon_vertices) const { 
    try{
      double val;
      vector<three_point> intersection_vertices = 
	get_convex_polygon_circle_intersection(polygon_vertices, *this, .5*diameter);
      if(intersection_vertices.size()<=2) val = 0;
      else val = get_area_of_polygon(intersection_vertices);
      return(val);
    } catch(...){
      cerr << "circular_aperture::convex_polygon_overlap error\n";
      throw(string("circular_aperture::convex_polygon_overlap"));
    }
  }

  namespace factory_register {
    const fits_keyval_set & get_annular_aperture_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "annular aperture"));
      return *fkvs;
    }
    
    AO_sim_base * create_annular_aperture(const iofits & iof) {
      return new annular_aperture(iof);
    }
  }

  const bool annular_aperture::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_annular_aperture_keyval_set(), 
				      factory_register::create_annular_aperture);

  annular_aperture::annular_aperture(){
    inner_diameter = 0;
    outer_diameter = 0;
  }

  annular_aperture::annular_aperture(const annular_aperture & annular_ap){
    this->operator=(annular_ap);
  }

  annular_aperture::annular_aperture(const char * filename){
    this->read(filename);
  }

  annular_aperture::annular_aperture(const iofits & iof){
    this->read(iof);
  }

  annular_aperture::annular_aperture(double in_diameter, double out_diameter){

    if(in_diameter < 0){
      cerr << "annular_aperture::annular_aperture error - "
	   << "inner diameter " << in_diameter 
	   << " supplied to constructor is less than zero\n";
      throw(string("annular_aperture::annular_aperture"));
    }

    if(out_diameter < 0){
      cerr << "annular_aperture::annular_aperture error - "
	   << "outer diameter " << out_diameter 
	   << " supplied to constructor is less than zero\n";
      throw(string("annular_aperture::annular_aperture"));
    }

    if(in_diameter > out_diameter) {
      cerr << "annular_aperture::annular_aperture error - "
	   << "inner diameter " << in_diameter 
	   << " supplied to constructor is greater than or equal to the outer diameter " 
	   << out_diameter << " supplied to constructor\n";
      throw(string("annular_aperture::annular_aperture"));
    }
    inner_diameter = in_diameter;
    outer_diameter = out_diameter;
  }

  annular_aperture & annular_aperture::operator=(const annular_aperture & annular_ap){
    if(this==&annular_ap)
      return(*this);
    inner_diameter = annular_ap.inner_diameter;
    outer_diameter = annular_ap.outer_diameter;
    this->aperture::operator=(annular_ap);
    return(*this);
  }

  void annular_aperture::read(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "annular_aperture::read - "
	   << "error opening file " << filename << endl;
      throw(string("annular_aperture::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "annular_aperture::read - "
	   << "error reading "
	   << this->unique_name() << " from file "
	   << filename << endl;
      throw(string("annular_aperture::read"));
    }
  }

  void annular_aperture::read(const iofits & iof) {
    if(!iof.key_exists("TYPE")){
      cerr << "annular_aperture::read error - "
	   << "unrecognized type of file\n";
      throw(string("annular_aperture::read"));
    }
    string type, comment;
    iof.read_key("TYPE", type, comment);
    if(type!=this->unique_name()){
      cerr << "annular_aperture::read error - file of type " 
	   << type << " rather than of type "
	   << this->unique_name() << endl;
      throw(string("annular_aperture::read"));
    }
    this->aperture::read(iof);
    iof.read_key("INDMTR", inner_diameter, comment);
    iof.read_key("OUTDMTR", outer_diameter, comment);
    // leave iof pointing to the next header, if it exists
    if(iof.get_hdu_num()!=iof.get_num_hdus())
      iof.movrel_hdu(1);
  }

  void annular_aperture::write(const char * filename) const {

    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "annular_aperture::write - "
	   << "error opening file " << filename << endl;
      throw(string("annular_aperture::write"));
    }

    try{this->write(iof);}
    catch(...){
      cerr << "annular_aperture::write - "
	   << "error writing "
	   << this->unique_name() << " to file " 
	   << filename << endl;
      throw(string("annular_aperture::write"));
    }
  }

  void annular_aperture::write(iofits & iof) const {
    this->aperture::write(iof);
    string type = this->unique_name();
    string comment = "object type";
    iof.write_key("TYPE", type, comment);
    iof.write_key("INDMTR", inner_diameter, "aperture inner diameter (meters)");
    iof.write_key("OUTDMTR", outer_diameter, "aperture outer diameter (meters)");
  }

  void annular_aperture::print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << this->unique_name()
       << "/" << "object type" << endl;
    this->aperture::print(os, prefix);
    os << prefix << "INDMTR     = " << setw(vlspc) << inner_diameter
       << "/" << "aperture inner diameter (meters)" << endl;
    os << prefix << "OUTDMTR    = " << setw(vlspc) << outer_diameter
       << "/" << "aperture outer diameter (meters)" << endl;
  }

  rectangular_region annular_aperture::get_covering_region(const three_frame & tf) const {
    if(fabs(dot_product(tf.z(), this->z()))<three_frame::precision){
      cerr << "annular_aperture::get_covering_region error - "
	   << "z axes of this aperture and three frame provided to this function are orthogonal\n";
      throw(string("annular_aperture::get_covering_region"));
    }
    circular_aperture circ_ap(outer_diameter);
    circ_ap.set_foreshortening(this->get_foreshortening());
    return(circ_ap.get_covering_region(tf));
  }

  void annular_aperture::transform(diffractive_wavefront<float> & wf) const {
    try{this->private_transform(wf);}
    catch(...){
      cerr << "annular_aperture::transform error - "
	   << "error transforming a float instantiation of diffractive_wavefront\n";
      throw(string("annular_aperture::transform"));
    }
  }

  void annular_aperture::transform(diffractive_wavefront<double> & wf) const {
    try{this->private_transform(wf);}
    catch(...){
      cerr << "annular_aperture::transform error - "
	   << "error transforming a double instantiation of diffractive_wavefront\n";
      throw(string("annular_aperture::transform"));
    }
  }

  /*
  void annular_aperture::transform(geometric_ray & gray) const {
    if(gray.is_lost()) return;
    three_vector odiff = static_cast<three_point>(*this) - static_cast<three_point>(gray);
    if(dot_product(odiff, this->z())!=0){
      cerr << "annular_aperture::transform error - "
	   << "geometric ray does not lie in plane of aperture\n";
      throw(string("annular_aperture::transform"));
    }
    double tmp = odiff.length_squared();
    if(tmp < .25*inner_diameter*inner_diameter || 
       tmp > .25*outer_diameter*outer_diameter)
      gray.set_lost(true);
  }

  void annular_aperture::transform(geometric_wavefront & gwf) const {
    long nrays = gwf.get_number_of_rays();
    geometric_ray * gray = this->get_wavefront_data(gwf);
    for(int i=0; i<nrays; i++)
      this->transform(gray[i]);
  }
  */

  template<class T>
  void annular_aperture::private_transform(diffractive_wavefront<T> & wf) const { 

    // perform the aperture check
    three_vector origin_offset, dx, dy;
    this->get_projected_wavefront_pixel_spacing(wf, origin_offset, dx, dy);

    // Check whether there is any overlap between the 
    // aperture and the diffractive_wavefront.  If not, zero the wavefront
    if(origin_offset.length()>outer_diameter){
      wf*=complex<double>(0,0);
      return;
    }

    // The halfpixel information
    vector<long> wf_axes = wf.get_axes();
    double x_halfpix=0, y_halfpix=0;
    int x_extrapix=1, y_extrapix=1;
    if(wf_axes[1]%2==0){
      x_halfpix = .5;
      x_extrapix = 0;
    }
    if(wf_axes[0]%2==0){
      y_halfpix = .5;
      y_extrapix = 0;
    }

    bool real_imag_storage = is_real_imag_storage(wf);
    bool interleaved_storage = is_interleaved_storage(wf);
    T * wfdata = get_wavefront_data(wf);

    int index;
    int nelem = wf_axes[0]*wf_axes[1];
    three_point center_pixel_coords;
    double val;
    three_vector tmps = .5*(dy+dx);
    three_vector tmpd = .5*(dy-dx);
    vector<three_point> pixel_vertices(4);
    double norm = cross_product(dx, dy).length();
    double buffer = tmps.length() > tmpd.length() ? 2*tmps.length() : 2*tmpd.length();
    double outer_upper_limit = (.5*this->outer_diameter+buffer)*(.5*this->outer_diameter+buffer);
    double outer_lower_limit = (.5*this->outer_diameter-buffer)*(.5*this->outer_diameter-buffer);
    double inner_upper_limit = (.5*this->inner_diameter+buffer)*(.5*this->inner_diameter+buffer);
    double inner_lower_limit = (.5*this->inner_diameter-buffer)*(.5*this->inner_diameter-buffer);
    double center_distance_squared;
    double outer_radius_squared = .25*outer_diameter*outer_diameter;
    double inner_radius_squared = .25*inner_diameter*inner_diameter;

    for(int i=-wf_axes[1]/2; i<wf_axes[1]/2+x_extrapix; i++){
      for(int j=-wf_axes[0]/2; j<wf_axes[0]/2+y_extrapix; j++){
	index = (i+wf_axes[1]/2)*wf_axes[0]+j+wf_axes[0]/2;

	center_pixel_coords = static_cast<const three_point>(*this) + 
	  (origin_offset + (i+x_halfpix)*dx + (j+y_halfpix)*dy);

	center_distance_squared = (center_pixel_coords - *this).length_squared();

	if(areal_weighting){
	  if(center_distance_squared > outer_upper_limit ||
	     center_distance_squared < inner_lower_limit) 
	    val = 0;
	  else if(center_distance_squared < outer_lower_limit &&
		  center_distance_squared > inner_upper_limit)
	    val = 1;
	  else {
	    pixel_vertices[0] = center_pixel_coords + tmps;
	    pixel_vertices[1] = center_pixel_coords + tmpd;
	    pixel_vertices[2] = center_pixel_coords - tmps;
	    pixel_vertices[3] = center_pixel_coords - tmpd;
	    
	    val = convex_polygon_overlap(pixel_vertices)/norm;
	  }
	} else {
	  val = 1;
	  if(center_distance_squared < inner_radius_squared ||
	     center_distance_squared > outer_radius_squared)
	    val = 0;
	}

	if(val!=1){
	  if(interleaved_storage){
	    if(val==0){
	      wfdata[2*index] = wfdata[2*index+1] = 0;
	      continue;
	    }
	    wfdata[2*index] *= val;
	    if(real_imag_storage)
	      wfdata[2*index+1] *= val;
	  } else {
	    if(val==0){
	      wfdata[2*index] = wfdata[2*index+1] = 0;
	      continue;
	    }
	    wfdata[index] *= val;
	    if(real_imag_storage)
	      wfdata[index+nelem] *= val;
	  }
	}
      }
    }
  }

  double annular_aperture::convex_polygon_overlap(const vector<three_point> & polygon_vertices) const { 
    try{
      double val = 0;
      vector<three_point> intersection_vertices = 
	get_convex_polygon_circle_intersection(polygon_vertices, *this, .5*outer_diameter);
      if(intersection_vertices.size()>2)
	val = get_area_of_polygon(intersection_vertices);
      intersection_vertices = 
	get_convex_polygon_circle_intersection(polygon_vertices, *this, .5*inner_diameter);
      if(intersection_vertices.size()>2)
	val -= get_area_of_polygon(intersection_vertices);
      return(val);
    } catch(...){
      cerr << "annular_aperture::convex_polygon_overlap error\n";
      throw(string("annular_aperture::convex_polygon_overlap"));
    }

  }

  namespace factory_register {
    const fits_keyval_set & get_rectangular_aperture_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "rectangular aperture"));
      return *fkvs;
    }
    
    AO_sim_base * create_rectangular_aperture(const iofits & iof) {
      return new rectangular_aperture(iof);
    }

  }

  const bool rectangular_aperture::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_rectangular_aperture_keyval_set(), 
				      factory_register::create_rectangular_aperture);

  rectangular_aperture::rectangular_aperture(){
    size = vector<double>(2,0);
  }

  rectangular_aperture::rectangular_aperture(const rectangular_aperture & rect_ap){
    this->operator=(rect_ap);
  }

  rectangular_aperture::rectangular_aperture(const char * filename){
    this->read(filename);
  }

  rectangular_aperture::rectangular_aperture(const iofits & iof){
    this->read(iof);
  }

  rectangular_aperture::rectangular_aperture(double x_size, double y_size){
    if(x_size < 0){
      cerr << "rectangular_aperture::rectangular_aperture error - "
	   << "x size " << x_size
	   << " supplied to constructor is less than zero\n";
      throw(string("rectangular_aperture::rectangular_aperture"));
    }
    if(y_size < 0){
      cerr << "rectangular_aperture::rectangular_aperture error - "
	   << "y size " << y_size
	   << " supplied to constructor is less than zero\n";
      throw(string("rectangular_aperture::rectangular_aperture"));
    }

    size.resize(2);
    size[0] = x_size;
    size[1] = y_size;
  }

  rectangular_aperture & rectangular_aperture::operator=(const rectangular_aperture & rect_ap){
    if(this==&rect_ap)
      return(*this);
    size = rect_ap.size;
    this->aperture::operator=(rect_ap);
    return(*this);
  }

  void rectangular_aperture::read(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "rectangular_aperture::read - "
	   << "error opening file " << filename << endl;
      throw(string("rectangular_aperture::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "rectangular_aperture::read - "
	   << "error reading "
	   << this->unique_name() << " from file "
	   << filename << endl;
      throw(string("rectangular_aperture::read"));
    }
  }

  void rectangular_aperture::read(const iofits & iof) {
    if(!iof.key_exists("TYPE")){
      cerr << "rectangular_aperture::read error - "
	   << "unrecognized type of file\n";
      throw(string("rectangular_aperture::read"));
    }
    string type, comment;
    iof.read_key("TYPE", type, comment);
    if(type!=this->unique_name()){
      cerr << "rectangular_aperture::read error - file of type " 
	   << type << " rather than of type "
	   << this->unique_name() << endl;
      throw(string("rectangular_aperture::read"));
    }
    this->aperture::read(iof);
    size.resize(2);
    iof.read_key("XSIZE", size[0], comment);
    iof.read_key("YSIZE", size[1], comment);
    // leave iof pointing to the next header, if it exists
    if(iof.get_hdu_num()!=iof.get_num_hdus())
      iof.movrel_hdu(1);
  }

  void rectangular_aperture::write(const char * filename) const {

    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "rectangular_aperture::write - "
	   << "error opening file " << filename << endl;
      throw(string("rectangular_aperture::write"));
    }

    try{this->write(iof);}
    catch(...){
      cerr << "rectangular_aperture::write - "
	   << "error writing "
	   << this->unique_name() << " to file " 
	   << filename << endl;
      throw(string("rectangular_aperture::write"));
    }
  }

  void rectangular_aperture::write(iofits & iof) const {
    this->aperture::write(iof);
    string type = this->unique_name();
    string comment = "object type";
    iof.write_key("TYPE", type, comment);
    iof.write_key("XSIZE", size[0], "aperture width (meters)");
    iof.write_key("YSIZE", size[1], "aperture height (meters)");
  }

  void rectangular_aperture::print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << this->unique_name()
       << "/" << "object type" << endl;
    this->aperture::print(os, prefix);
    os << prefix << "XSIZE      = " << setw(vlspc) << size[0]
       << "/" << "aperture width (meters)" << endl;
    os << prefix << "YSIZE      = " << setw(vlspc) << size[1]
       << "/" << "aperture height (meters)" << endl;
  }

  rectangular_region rectangular_aperture::get_covering_region(const three_frame & tf) const {

    if(fabs(dot_product(tf.z(), this->z()))<three_frame::precision){
      cerr << "rectangular_aperture::get_covering_region error - "
	   << "z axes of this aperture and three frame provided to this function are orthogonal\n";
      throw(string("rectangular_aperture::get_covering_region"));
    }

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

  void rectangular_aperture::transform(diffractive_wavefront<float> & wf) const {
    try{this->private_transform(wf);}
    catch(...){
      cerr << "rectangular_aperture::transform error - "
	   << "error transforming a float instantiation of diffractive_wavefront\n";
      throw(string("rectangular_aperture::transform"));
    }
  }

  void rectangular_aperture::transform(diffractive_wavefront<double> & wf) const {
    try{this->private_transform(wf);}
    catch(...){
      cerr << "rectangular_aperture::transform error - "
	   << "error transforming a double instantiation of diffractive_wavefront\n";
      throw(string("rectangular_aperture::transform"));
    }
  }

  /*
  void rectangular_aperture::transform(geometric_ray & gray) const {
    if(gray.is_lost()) return;
    three_vector odiff = static_cast<three_point>(*this) - static_cast<three_point>(gray);
    if(dot_product(odiff, this->z())!=0){
      cerr << "rectangular_aperture::transform error - "
	   << "geometric ray does not lie in plane of aperture\n";
      throw(string("rectangular_aperture::transform"));
    }
    if(fabs(odiff.x(*this))>size[0]/2.0 ||
       fabs(odiff.y(*this))>size[1]/2.0)
      gray.set_lost(true);
  }

  void rectangular_aperture::transform(geometric_wavefront & gwf) const {
    long nrays = gwf.get_number_of_rays();
    geometric_ray * gray = this->get_wavefront_data(gwf);
    for(int i=0; i<nrays; i++)
      this->transform(gray[i]);
  }
  */

  template<class T>
  void rectangular_aperture::private_transform(diffractive_wavefront<T> & wf) const { 
    
    // perform the aperture check
    three_vector origin_offset, dx, dy;
    this->get_projected_wavefront_pixel_spacing(wf, origin_offset, dx, dy);

    // Check whether there is any overlap between the 
    // aperture and the diffractive_wavefront.  If not, zero the wavefront
    if(origin_offset.x(*this)>size[0] &&
       origin_offset.y(*this)>size[1]){
      wf*=complex<double>(0,0);
      return;
    }

    // The halfpixel information
    vector<long> wf_axes = wf.get_axes();
    double x_halfpix=0, y_halfpix=0;
    int x_extrapix=1, y_extrapix=1;
    if(wf_axes[1]%2==0){
      x_halfpix = .5;
      x_extrapix = 0;
    }
    if(wf_axes[0]%2==0){
      y_halfpix = .5;
      y_extrapix = 0;
    }

    bool real_imag_storage = is_real_imag_storage(wf);
    bool interleaved_storage = is_interleaved_storage(wf);
    T * wfdata = get_wavefront_data(wf);

    int index;
    int nelem = wf_axes[0]*wf_axes[1];
    three_point center_pixel_coords;
    double val;
    three_vector tmps = .5*(dy+dx);
    three_vector tmpd = .5*(dy-dx);
    vector<three_point> pixel_vertices(4);
    double norm = cross_product(dx, dy).length();
    double buffer = tmps.length() > tmpd.length() ? 2*tmps.length() : 2*tmpd.length();
    double upper_limit_a = (.5*this->size[0]+buffer);
    double lower_limit_a = (.5*this->size[0]-buffer);
    double upper_limit_b = (.5*this->size[1]+buffer);
    double lower_limit_b = (.5*this->size[1]-buffer);

    for(int i=-wf_axes[1]/2; i<wf_axes[1]/2+x_extrapix; i++){
      for(int j=-wf_axes[0]/2; j<wf_axes[0]/2+y_extrapix; j++){
	index = (i+wf_axes[1]/2)*wf_axes[0]+j+wf_axes[0]/2;
	center_pixel_coords = static_cast<const three_point>(*this)
	  + (origin_offset + (i+x_halfpix)*dx + (j+y_halfpix)*dy);

	if(areal_weighting){
	  if(fabs(center_pixel_coords.x(*this))>upper_limit_a ||
	     fabs(center_pixel_coords.y(*this))>upper_limit_b)
	    val = 0;
	  else if(fabs(center_pixel_coords.x(*this))<lower_limit_a &&
		  fabs(center_pixel_coords.y(*this))<lower_limit_b)
	    val = 1;
	  else {
	    pixel_vertices[0] = center_pixel_coords + tmps;
	    pixel_vertices[1] = center_pixel_coords + tmpd;
	    pixel_vertices[2] = center_pixel_coords - tmps;
	    pixel_vertices[3] = center_pixel_coords - tmpd;
	    
	    val = convex_polygon_overlap(pixel_vertices)/norm;
	  }
	} else {
	  val = 1;
	  if(fabs(center_pixel_coords.x(*this)) > .5*this->size[0] ||
	     fabs(center_pixel_coords.y(*this)) > .5*this->size[1])
	    val = 0;
	}

	if(val!=1){
	  if(interleaved_storage){
	    if(val==0){
	      wfdata[2*index] = wfdata[2*index+1] = 0;
	      continue;
	    }
	    wfdata[2*index] *= val;
	    if(real_imag_storage)
	      wfdata[2*index+1] *= val;
	  } else {
	    if(val==0){
	      wfdata[2*index] = wfdata[2*index+1] = 0;
	      continue;
	    }
	    wfdata[index] *= val;
	    if(real_imag_storage)
	      wfdata[index+nelem] *= val;
	  }
	}
      }
    }
  }

  double rectangular_aperture::convex_polygon_overlap(const vector<three_point> & polygon_vertices) const { 
    try{
      double val;
      vector<three_point> rectangle_vertices(4);
      rectangle_vertices[0] = *this + .5*this->size[0]*this->x() + .5*this->size[1]*this->y();
      rectangle_vertices[1] = *this - .5*this->size[0]*this->x() + .5*this->size[1]*this->y();
      rectangle_vertices[2] = *this - .5*this->size[0]*this->x() - .5*this->size[1]*this->y();
      rectangle_vertices[3] = *this + .5*this->size[0]*this->x() - .5*this->size[1]*this->y();
      
      vector<three_point> intersection_vertices = 
	get_convex_polygon_intersection(polygon_vertices, rectangle_vertices);
      if(intersection_vertices.size()<=2) val = 0;
      else val = get_area_of_polygon(intersection_vertices);
      return(val);
    } catch(...){
      cerr << "rectangular_aperture::convex_polygon_overlap error\n";
      throw(string("rectangular_aperture::convex_polygon_overlap"));
    }
  }
  
  namespace factory_register {
    const fits_keyval_set & get_spidered_annular_aperture_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "spidered annular aperture"));
      return *fkvs;
    }
    
    AO_sim_base * create_spidered_annular_aperture(const iofits & iof) {
      return new spidered_annular_aperture(iof);
    }
  }

  const bool spidered_annular_aperture::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_spidered_annular_aperture_keyval_set(), 
				      factory_register::create_spidered_annular_aperture);

  spidered_annular_aperture::spidered_annular_aperture(){
    nspiders = 0;
    spider_width = 0;
  }

  spidered_annular_aperture::spidered_annular_aperture(const spidered_annular_aperture & spidered_annular_ap){
    this->operator=(spidered_annular_ap);
  }

  spidered_annular_aperture::spidered_annular_aperture(const char * filename){
    this->read(filename);
  }

  spidered_annular_aperture::spidered_annular_aperture(const iofits & iof){
    this->read(iof);
  }

  spidered_annular_aperture::spidered_annular_aperture(double in_diameter, 
						       double out_diameter,
						       int in_nspiders, 
						       double in_spider_width) :
    annular_aperture(in_diameter, out_diameter) {

    if(in_nspiders <= 0){
      cerr << "spidered_annular_aperture::spidered_annular_aperture error - "
	   << "number of spiders " << in_nspiders
	   << " supplied to constructor is less than or equal to zero\n";
      throw(string("spidered_annular_aperture::spidered_annular_aperture"));
    }

    if(in_spider_width < 0){
      cerr << "spidered_annular_aperture::spidered_annular_aperture error - "
	   << "spider width " << in_spider_width 
	   << " supplied to constructor is less than zero\n";
      throw(string("spidered_annular_aperture::spidered_annular_aperture"));
    }

    nspiders = in_nspiders;
    spider_width = in_spider_width;
  }

  spidered_annular_aperture & spidered_annular_aperture::operator=(const spidered_annular_aperture & spidered_annular_ap){
    if(this==&spidered_annular_ap)
      return(*this);
    nspiders = spidered_annular_ap.nspiders;
    spider_width = spidered_annular_ap.spider_width;
    this->annular_aperture::operator=(spidered_annular_ap);
    return(*this);
  }

  void spidered_annular_aperture::read(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "spidered_annular_aperture::read - "
	   << "error opening file " << filename << endl;
      throw(string("spidered_annular_aperture::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "spidered_annular_aperture::read - "
	   << "error reading "
	   << this->unique_name() << " from file "
	   << filename << endl;
      throw(string("spidered_annular_aperture::read"));
    }
  }

  void spidered_annular_aperture::read(const iofits & iof) {
    if(!iof.key_exists("TYPE")){
      cerr << "spidered_annular_aperture::read error - "
	   << "unrecognized type of file\n";
      throw(string("spidered_annular_aperture::read"));
    }
    string type, comment;
    iof.read_key("TYPE", type, comment);
    if(type!=this->unique_name()){
      cerr << "spidered_annular_aperture::read error - file of type " 
	   << type << " rather than of type "
	   << this->unique_name() << endl;
      throw(string("spidered_annular_aperture::read"));
    }
    this->aperture::read(iof);
    iof.read_key("INDMTR", inner_diameter, comment);
    iof.read_key("OUTDMTR", outer_diameter, comment);
    iof.read_key("NSPIDERS", nspiders, comment);
    iof.read_key("SPDRWDTH", spider_width, comment);
    // leave iof pointing to the next header, if it exists
    if(iof.get_hdu_num()!=iof.get_num_hdus())
      iof.movrel_hdu(1);
  }

  void spidered_annular_aperture::write(const char * filename) const {

    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "spidered_annular_aperture::write - "
	   << "error opening file " << filename << endl;
      throw(string("spidered_annular_aperture::write"));
    }

    try{this->write(iof);}
    catch(...){
      cerr << "spidered_annular_aperture::write - "
	   << "error writing "
	   << this->unique_name() << " to file " 
	   << filename << endl;
      throw(string("spidered_annular_aperture::write"));
    }
  }

  void spidered_annular_aperture::write(iofits & iof) const {
    this->aperture::write(iof);
    string type = this->unique_name();
    string comment = "object type";
    iof.write_key("TYPE", type, comment);
    iof.write_key("INDMTR", inner_diameter, "aperture inner diameter (meters)");
    iof.write_key("OUTDMTR", outer_diameter, "aperture outer diameter (meters)");
    iof.write_key("NSPIDERS", nspiders, "number of spiders");
    iof.write_key("SPDRWDTH", spider_width, "width of spiders (meters)");
  }

  void spidered_annular_aperture::print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << this->unique_name()
       << "/" << "object type" << endl;
    this->aperture::print(os, prefix);
    os << prefix << "INDMTR     = " << setw(vlspc) << inner_diameter
       << "/" << "aperture inner diameter (meters)" << endl;
    os << prefix << "OUTDMTR    = " << setw(vlspc) << outer_diameter
       << "/" << "aperture outer diameter (meters)" << endl;
    os << prefix << "NSPIDERS   = " << setw(vlspc) << nspiders
       << "/" << "number of spiders" << endl;
    os << prefix << "SPDRWDTH   = " << setw(vlspc) << spider_width
       << "/" << "width of spiders (meters)" << endl;
  }

  void spidered_annular_aperture::transform(diffractive_wavefront<float> & wf) const {
    try{this->private_transform(wf);}
    catch(...){
      cerr << "spidered_annular_aperture::transform error - "
	   << "error transforming a float instantiation of diffractive_wavefront\n";
      throw(string("spidered_annular_aperture::transform"));
    }
  }

  void spidered_annular_aperture::transform(diffractive_wavefront<double> & wf) const {
    try{this->private_transform(wf);}
    catch(...){
      cerr << "spidered_annular_aperture::transform error - "
	   << "error transforming a double instantiation of diffractive_wavefront\n";
      throw(string("spidered_annular_aperture::transform"));
    }
  }

  /*
  void spidered_annular_aperture::transform(geometric_ray & gray) const {
    if(gray.is_lost()) return;
    three_vector odiff = static_cast<three_point>(*this) - static_cast<three_point>(gray);
    if(dot_product(odiff, this->z())!=0){
      cerr << "spidered_annular_aperture::transform error - "
	   << "geometric ray does not lie in plane of aperture\n";
      throw(string("spidered_annular_aperture::transform"));
    }
    double tmp = odiff.length_squared();
    if(tmp < .25*inner_diameter*inner_diameter || 
       tmp > .25*outer_diameter*outer_diameter){
      gray.set_lost(true);
      return;
    }
  
    three_rotation spider_rot(*this, this->z(), 2*M_PI/(double)nspiders);
    three_vector spider_vec(0,1,0,*this);
    for(int i=0; i<nspiders; i++){
      if((odiff - dot_product(odiff, spider_vec)*spider_vec).length()>spider_width/2.0){
	gray.set_lost(true);
	return;
      }
      spider_rot.transform(spider_vec);
    }
  }

  void spidered_annular_aperture::transform(geometric_wavefront & gwf) const {
    long nrays = gwf.get_number_of_rays();
    geometric_ray * gray = this->get_wavefront_data(gwf);
    for(int i=0; i<nrays; i++)
      this->transform(gray[i]);
  }
  */

  template<class T>
  void spidered_annular_aperture::private_transform(diffractive_wavefront<T> & wf) const { 
    
    // perform the aperture check
    three_vector origin_offset, dx, dy;
    this->get_projected_wavefront_pixel_spacing(wf, origin_offset, dx, dy);

    // Check whether there is any overlap between the 
    // aperture and the diffractive_wavefront.  If not, zero the wavefront
    if(origin_offset.length()>outer_diameter){
      wf*=complex<double>(0,0);
      return;
    }

    // The halfpixel information
    vector<long> wf_axes = wf.get_axes();
    double x_halfpix=0, y_halfpix=0;
    int x_extrapix=1, y_extrapix=1;
    if(wf_axes[1]%2==0){
      x_halfpix = .5;
      x_extrapix = 0;
    }
    if(wf_axes[0]%2==0){
      y_halfpix = .5;
      y_extrapix = 0;
    }

    bool real_imag_storage = is_real_imag_storage(wf);
    bool interleaved_storage = is_interleaved_storage(wf);
    T * wfdata = get_wavefront_data(wf);

    int index;
    int nelem = wf_axes[0]*wf_axes[1];
    three_point center_pixel_coords;
    double val;
    three_vector tmps = .5*(dy+dx);
    three_vector tmpd = .5*(dy-dx);
    vector<three_point> pixel_vertices(4);
    double norm = cross_product(dx, dy).length();
    double buffer = tmps.length() > tmpd.length() ? 2*tmps.length() : 2*tmpd.length();
    double outer_upper_limit = (.5*this->outer_diameter+buffer)*(.5*this->outer_diameter+buffer);
    double outer_lower_limit = (.5*this->outer_diameter-buffer)*(.5*this->outer_diameter-buffer);
    double inner_upper_limit = (.5*this->inner_diameter+buffer)*(.5*this->inner_diameter+buffer);
    double inner_lower_limit = (.5*this->inner_diameter-buffer)*(.5*this->inner_diameter-buffer);
    double center_distance_squared;
    double outer_radius_squared = .25*outer_diameter*outer_diameter;
    double inner_radius_squared = .25*inner_diameter*inner_diameter;

    // make unit vectors that lie along the directions of the spiders
    double spider_angle = 2*M_PI/(double)nspiders;
    three_rotation spider_rot(*this, this->z(), spider_angle);
    vector<three_vector> spider_unit_vectors(nspiders);
    spider_unit_vectors[0] = three_vector(0,1,0,*this);
    for(int i=1; i<nspiders; i++) {
      spider_unit_vectors[i] = spider_unit_vectors[i-1];
      spider_rot.transform(spider_unit_vectors[i]);
    }
    
    three_vector center_pixel_vector;
    double spider_limit = .25*spider_width*spider_width;
    double weighted_spider_upper_limit = (.5*spider_width + buffer)*(.5*spider_width + buffer);
    double weighted_spider_lower_limit = (.5*spider_width - buffer)*(.5*spider_width - buffer);
    double spider_distance_squared, closest_spider_distance_squared;

    for(int i=-wf_axes[1]/2; i<wf_axes[1]/2+x_extrapix; i++){
      for(int j=-wf_axes[0]/2; j<wf_axes[0]/2+y_extrapix; j++){
	index = (i+wf_axes[1]/2)*wf_axes[0]+j+wf_axes[0]/2;
	center_pixel_coords = static_cast<const three_point>(*this) +
	  (origin_offset + (i+x_halfpix)*dx + (j+y_halfpix)*dy);

	center_distance_squared = (center_pixel_coords - *this).length_squared();

	if(areal_weighting){
	  if(center_distance_squared > outer_upper_limit ||
	     center_distance_squared < inner_lower_limit) 
	    val = 0;
	  else {
	    
	    center_pixel_vector = center_pixel_coords - *this;

	    if(dot_product(spider_unit_vectors[0], center_pixel_vector)<0)
	      closest_spider_distance_squared = center_pixel_vector.length();
	    else 
	      closest_spider_distance_squared = cross_product(spider_unit_vectors[0],
							      center_pixel_vector).length_squared();
	    for(int k=1; k<nspiders; k++){
	      if(dot_product(spider_unit_vectors[k], center_pixel_vector)<0)
		spider_distance_squared = center_pixel_vector.length();
	      else
		spider_distance_squared = cross_product(spider_unit_vectors[k],
							center_pixel_vector).length_squared();
	      if(closest_spider_distance_squared>spider_distance_squared)
		closest_spider_distance_squared=spider_distance_squared;
	    }

	    if(closest_spider_distance_squared < weighted_spider_lower_limit)
	      val = 0;
	    else if(center_distance_squared < outer_lower_limit &&
	       center_distance_squared > inner_upper_limit &&
	       closest_spider_distance_squared > weighted_spider_upper_limit)
	      val = 1;
	    else {
	      pixel_vertices[0] = center_pixel_coords + tmps;
	      pixel_vertices[1] = center_pixel_coords + tmpd;
	      pixel_vertices[2] = center_pixel_coords - tmps;
	      pixel_vertices[3] = center_pixel_coords - tmpd;
	      
	      try{val = convex_polygon_overlap(pixel_vertices)/norm;}
	      catch(...){
		cerr << "spidered_annular_aperture::private_transform error\n";
		stringstream ss;
		for(int k=0; k<4; k++){
		  ss.str("");
		  ss << "\tpixel vertex " << k << " ";
		  pixel_vertices[k].print(cerr, ss.str().c_str());
		}
		throw(string("spidered_annular_aperture::private_transform"));
	      }
	    }
	  }
	} else {
	  val = 1;
	  if(center_distance_squared < inner_radius_squared ||
	     center_distance_squared > outer_radius_squared)
	    val = 0;
	  else {
	    center_pixel_vector = center_pixel_coords - *this;
	    for(int k=0; k<nspiders; k++){
	      if(cross_product(center_pixel_vector,
			       spider_unit_vectors[k]).length_squared()
		 < spider_limit){
		val = 0;
		break;
	      }
	    }
	  }
	}

	if(val!=1){
	  if(interleaved_storage){
	    if(val==0){
	      wfdata[2*index] = wfdata[2*index+1] = 0;
	      continue;
	    }
	    wfdata[2*index] *= val;
	    if(real_imag_storage)
	      wfdata[2*index+1] *= val;
	  } else {
	    if(val==0){
	      wfdata[2*index] = wfdata[2*index+1] = 0;
	      continue;
	    }
	    wfdata[index] *= val;
	    if(real_imag_storage)	
	      wfdata[index+nelem] *= val;
	  }
	}
      }
    }
  }

  double spidered_annular_aperture::convex_polygon_overlap(const vector<three_point> & polygon_vertices) const { 

    try{
      double polygon_overlap_with_outer_diameter;
      double polygon_overlap_with_inner_diameter = 0;
      double polygon_overlap_with_nearest_spider = 0;
      double inner_diameter_overlap_with_nearest_spider = 0;
      
      vector<three_point> outer_intersection_vertices = 
	get_convex_polygon_circle_intersection(polygon_vertices, *this, .5*outer_diameter);
      
      if(outer_intersection_vertices.size()<=2) return(0);
      polygon_overlap_with_outer_diameter = get_area_of_polygon(outer_intersection_vertices);

      vector<three_point> inner_intersection_vertices = 
	get_convex_polygon_circle_intersection(outer_intersection_vertices, *this, .5*inner_diameter);
      if(inner_intersection_vertices.size()>2) 
	polygon_overlap_with_inner_diameter = get_area_of_polygon(inner_intersection_vertices);
      
      double spider_angle = 2*M_PI/(double)nspiders;
      three_rotation spider_rot(*this, this->z(), spider_angle);
      three_vector spider_unit_vector(0,1,0,*this);
      int closest_spider_index = 0;
      double spider_dot_product;
      double closest_spider_dot_product = dot_product(polygon_vertices[0]-*this,spider_unit_vector);
      for(int i=1; i<nspiders; i++) {
	spider_rot.transform(spider_unit_vector);
	spider_dot_product = dot_product(polygon_vertices[0]-*this, spider_unit_vector);
	if(closest_spider_dot_product<spider_dot_product){
	  closest_spider_index = i;
	  closest_spider_dot_product = spider_dot_product;
	}
      }

      vector<three_point> spider_vertices(4);
      spider_vertices[0] = *this - .5*spider_width*this->x();
      spider_vertices[1] = *this + .5*spider_width*this->x();
      spider_vertices[2] = spider_vertices[1] + 1.1*outer_diameter*this->y();
      spider_vertices[3] = spider_vertices[0] + 1.1*outer_diameter*this->y();
    
      spider_rot = three_rotation(*this, this->z(), closest_spider_index*spider_angle);
      for(int i=0; i<4; i++)
	spider_rot.transform(spider_vertices[i]);
      
      vector<three_point> spider_intersection_vertices = 
	get_convex_polygon_intersection(outer_intersection_vertices, spider_vertices);
      if(spider_intersection_vertices.size()>2) 
	polygon_overlap_with_nearest_spider = get_area_of_polygon(spider_intersection_vertices);
      
      if(polygon_overlap_with_nearest_spider>0){
	if(polygon_overlap_with_inner_diameter>0){
	  vector<three_point> intersection_vertices = 
	    get_convex_polygon_intersection(inner_intersection_vertices, spider_intersection_vertices);
	  if(intersection_vertices.size()>2) 
	  inner_diameter_overlap_with_nearest_spider = get_area_of_polygon(intersection_vertices);
	}
      } 

      return(polygon_overlap_with_outer_diameter - 
	     polygon_overlap_with_inner_diameter -
	     polygon_overlap_with_nearest_spider +
	     inner_diameter_overlap_with_nearest_spider);
    } catch(...){
      cerr << "spidered_annular_aperture::convex_polygon_overlap error\n";
      throw(string("spidered_annular_aperture::convex_polygon_overlap"));
    }
  }


  namespace factory_register {
    const fits_keyval_set & get_hexagonal_aperture_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "hexagonal aperture"));
      return *fkvs;
    }
    
    AO_sim_base * create_hexagonal_aperture(const iofits & iof) {
      return new hexagonal_aperture(iof);
    }
  }

  const bool hexagonal_aperture::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_hexagonal_aperture_keyval_set(), 
				      factory_register::create_hexagonal_aperture);

  hexagonal_aperture::hexagonal_aperture(){
    edge_length = 0;
  }

  hexagonal_aperture::hexagonal_aperture(const hexagonal_aperture & hex_ap){
    this->operator=(hex_ap);
  }

  hexagonal_aperture::hexagonal_aperture(const char * filename){
    this->read(filename);
  }

  hexagonal_aperture::hexagonal_aperture(const iofits & iof){
    this->read(iof);
  }

  hexagonal_aperture::hexagonal_aperture(double in_edge_length){
    if(in_edge_length < 0){
      cerr << "hexagonal_aperture::hexagonal_aperture error - "
	   << "edge length " << in_edge_length
	   << " supplied to constructor is less than zero\n";
      throw(string("hexagonal_aperture::hexagonal_aperture"));
    }
    edge_length = in_edge_length;
  }

  hexagonal_aperture & hexagonal_aperture::operator=(const hexagonal_aperture & hex_ap){
    if(this==&hex_ap)
      return(*this);
    edge_length = hex_ap.edge_length;  
    this->aperture::operator=(hex_ap);
    return(*this);
  }

  void hexagonal_aperture::read(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "hexagonal_aperture::read - "
	   << "error opening file " << filename << endl;
      throw(string("hexagonal_aperture::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "hexagonal_aperture::read - "
	   << "error reading "
	   << this->unique_name() << " from file "
	   << filename << endl;
      throw(string("hexagonal_aperture::read"));
    }
  }

  void hexagonal_aperture::read(const iofits & iof) {
    if(!iof.key_exists("TYPE")){
      cerr << "hexagonal_aperture::read error - "
	   << "unrecognized type of file\n";
      throw(string("hexagonal_aperture::read"));
    }
    string type, comment;
    iof.read_key("TYPE", type, comment);
    if(type!=this->unique_name()){
      cerr << "hexagonal_aperture::read error - file of type " 
	   << type << " rather than of type "
	   << this->unique_name() << endl;
      throw(string("hexagonal_aperture::read"));
    }
    this->aperture::read(iof);
    iof.read_key("EDGELNTH", edge_length, comment);
    // leave iof pointing to the next header, if it exists
    if(iof.get_hdu_num()!=iof.get_num_hdus())
      iof.movrel_hdu(1);
  }

  void hexagonal_aperture::write(const char * filename) const {

    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "hexagonal_aperture::write - "
	   << "error opening file " << filename << endl;
      throw(string("hexagonal_aperture::write"));
    }

    try{this->write(iof);}
    catch(...){
      cerr << "hexagonal_aperture::write - "
	   << "error writing "
	   << this->unique_name() << " to file " 
	   << filename << endl;
      throw(string("hexagonal_aperture::write"));
    }
  }

  void hexagonal_aperture::write(iofits & iof) const {
    this->aperture::write(iof);
    string type = this->unique_name();
    string comment = "object type";
    iof.write_key("TYPE", type, comment);
    iof.write_key("EDGELNTH", edge_length, "hexagon edge length (meters)");
  }

  void hexagonal_aperture::print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << this->unique_name()
       << "/" << "object type" << endl;
    this->aperture::print(os, prefix);
    os << prefix << "EDGELNTH   = " << setw(vlspc) << edge_length
       << "/" << "hexagon edge length (meters)" << endl;
  }

  void hexagonal_aperture::transform(diffractive_wavefront<float> & wf) const {
    try{this->private_transform(wf);}
    catch(...){
      cerr << "hexagonal_aperture::transform error - "
	   << "error transforming a float instantiation of diffractive_wavefront\n";
      throw(string("hexagonal_aperture::transform"));
    }
  }

  rectangular_region hexagonal_aperture::get_covering_region(const three_frame & tf) const {

    if(fabs(dot_product(tf.z(), this->z()))<three_frame::precision){
      cerr << "hexagonal_aperture::get_covering_region error - "
	   << "z axes of this aperture and three frame provided to this function are orthogonal\n";
      throw(string("hexagonal_aperture::get_covering_region"));
    }

    vector<three_point> hex_corners(6);
    hex_corners[0] = *this + three_vector(0, edge_length, 0, *this);
    hex_corners[1] = *this + three_vector(-edge_length*sqrt(3.0)/2.0, edge_length/2.0, 0, *this);
    hex_corners[2] = *this + three_vector(-edge_length*sqrt(3.0)/2.0, -edge_length/2.0, 0, *this);
    hex_corners[3] = *this + three_vector(0, -edge_length, 0, *this);
    hex_corners[4] = *this + three_vector(edge_length*sqrt(3.0)/2.0, -edge_length/2.0, 0, *this);
    hex_corners[5] = *this + three_vector(edge_length*sqrt(3.0)/2.0, edge_length/2.0, 0, *this);

    if(foreshortening){
      // project off the component along tf.z() to get
      // a three point in the XY plane of tf
      int aligned = 1;
      if(dot_product(this->z(), tf.z())<1) aligned = -1;
      for(int i=0; i<6; i++)
	hex_corners[i] -= aligned*dot_product(hex_corners[i]-*this, tf.z())*tf.z();
    } else {
      three_vector rotation_axis = cross_product(this->z(), tf.z());
      if(rotation_axis.length()>three_frame::precision){
	three_rotation trot(*this, rotation_axis, rotation_axis.length());
	for(int i=0; i<6; i++)
	  trot.transform(hex_corners[i]);
      }
    }
    
    // Next we find the extremal coordinates of the points identified above
    // in the three_frame tf
    double xmin = hex_corners[0].x(tf);
    double xmax = xmin;
    double ymin = hex_corners[0].y(tf);
    double ymax = ymin;
    for(int i=1; i<6; i++){
      if(hex_corners[i].x(tf)>xmax) xmax = hex_corners[i].x(tf);
      if(hex_corners[i].x(tf)<xmin) xmin = hex_corners[i].x(tf);
      if(hex_corners[i].y(tf)>ymax) ymax = hex_corners[i].y(tf);
      if(hex_corners[i].y(tf)<ymin) ymin = hex_corners[i].y(tf);
    }
    
    // Finally we make a rectangular region from these limiting coordinates
    vector<three_point> tp(4);
    tp[0] = three_point(xmin, ymin, 0, tf);
    tp[1] = three_point(xmax, ymin, 0, tf);
    tp[2] = three_point(xmax, ymax, 0, tf);
    tp[3] = three_point(xmin, ymax, 0, tf);
    return(rectangular_region(tp));

  }

  void hexagonal_aperture::transform(diffractive_wavefront<double> & wf) const {
    try{this->private_transform(wf);}
    catch(...){
      cerr << "hexagonal_aperture::transform error - "
	   << "error transforming a double instantiation of diffractive_wavefront\n";
      throw(string("hexagonal_aperture::transform"));
    }
  }

  /*
  void hexagonal_aperture::transform(geometric_wavefront & gwf) const {
    long nrays = gwf.get_number_of_rays();
    geometric_ray * gray = this->get_wavefront_data(gwf);
    for(int i=0; i<nrays; i++)
      this->transform(gray[i]);
  }

  void hexagonal_aperture::transform(geometric_ray & gray) const {
    if(gray.is_lost()) return;
    three_vector odiff = static_cast<three_point>(*this) - static_cast<three_point>(gray);
    if(dot_product(odiff, this->z())!=0){
      cerr << "hexagonal_aperture::transform error - "
	   << "geometric ray does not lie in plane of aperture\n";
      throw(string("hexagonal_aperture::transform"));
    }
  
    if(odiff.length_squared()>edge_length*edge_length ||
       odiff.length_squared()<3*edge_length*edge_length/4.0){
      gray.set_lost(true);
      return;
    }

    int max_overlap_index = 0;
    double overlap;
    three_rotation hex_rot(*this, this->z(), 2*M_PI/6.0);
    three_vector hex_vec[6];
    hex_vec[0] = three_vector(0,1,0,*this);
    double max_overlap = dot_product(hex_vec[0], odiff);
    for(int k=1; k<6; k++){
      hex_vec[k] = hex_vec[k-1];
      hex_rot.transform(hex_vec[k]);
      if((overlap=dot_product(hex_vec[k],odiff))>max_overlap){
	max_overlap_index = k;
	max_overlap = overlap;
      }
    }	  
    // max_overlap_index now labels the index of the triangle
    // in which pixel_coord lies
    if(dot_product(hex_vec[(max_overlap_index+1)%6],odiff) <
       dot_product(hex_vec[(max_overlap_index+5)%6],odiff))
      max_overlap_index = (max_overlap_index+5)%6;

    three_vector edge_vector = hex_vec[(max_overlap_index+1)%6] - hex_vec[max_overlap_index];
    double cos_beta = dot_product(edge_vector, odiff);
    double sin_beta = sqrt(1-cos_beta*cos_beta);
    double distance_to_edge = sin(M_PI/3.0)*edge_length/sin_beta;
    if(distance_to_edge < odiff.length())
      gray.set_lost(true);
  }
  */

  template<class T>
  void hexagonal_aperture::private_transform(diffractive_wavefront<T> & wf) const { 
    
    // perform the aperture check
    three_vector origin_offset, dx, dy;
    this->get_projected_wavefront_pixel_spacing(wf, origin_offset, dx, dy);

    // Check whether there is any overlap between the 
    // aperture and the diffractive_wavefront.  If not, zero the wavefront
    if(origin_offset.length()>2*edge_length){
      wf*=complex<double>(0,0);
      return;
    }

    // The halfpixel information
    vector<long> wf_axes = wf.get_axes();
    double x_halfpix=0, y_halfpix=0;
    int x_extrapix=1, y_extrapix=1;
    if(wf_axes[1]%2==0){
      x_halfpix = .5;
      x_extrapix = 0;
    }
    if(wf_axes[0]%2==0){
      y_halfpix = .5;
      y_extrapix = 0;
    }

    bool real_imag_storage = is_real_imag_storage(wf);
    bool interleaved_storage = is_interleaved_storage(wf);
    T * wfdata = get_wavefront_data(wf);

    int index;
    int nelem = wf_axes[0]*wf_axes[1];
    three_point center_pixel_coords;
    double val;
    three_vector tmps = .5*(dy+dx);
    three_vector tmpd = .5*(dy-dx);
    vector<three_point> pixel_vertices(4);
    double norm = cross_product(dx, dy).length();
    double buffer = tmps.length() > tmpd.length() ? 2*tmps.length() : 2*tmpd.length();
    double upper_limit = (edge_length+buffer)*(edge_length+buffer);
    double lower_limit = (.5*sqrt(3.0)*edge_length-buffer)*(.5*sqrt(3.0)*edge_length-buffer);
    double center_distance_squared;

    // definitions for non-areal weighting case
    double sqrt_three = sqrt(3.0);
    double coord_lim = sqrt_three*edge_length/2.0;
    double abs_xcoord, abs_ycoord, tmp;
    three_vector tv(.5,.5*sqrt_three,0,*this); 

    for(int i=-wf_axes[1]/2; i<wf_axes[1]/2+x_extrapix; i++){
      for(int j=-wf_axes[0]/2; j<wf_axes[0]/2+y_extrapix; j++){
	index = (i+wf_axes[1]/2)*wf_axes[0]+j+wf_axes[0]/2;

	center_pixel_coords = static_cast<const three_point>(*this) + 
	  (origin_offset + (i+x_halfpix)*dx + (j+y_halfpix)*dy);

	center_distance_squared = (center_pixel_coords - *this).length_squared();

	if(areal_weighting){
	  if(center_distance_squared>upper_limit) 
	    val = 0;
	  else if(center_distance_squared<lower_limit) 
	    val = 1;
	  else {
	    pixel_vertices[0] = center_pixel_coords + tmps;
	    pixel_vertices[1] = center_pixel_coords + tmpd;
	    pixel_vertices[2] = center_pixel_coords - tmps;
	    pixel_vertices[3] = center_pixel_coords - tmpd;
	    
	    val = convex_polygon_overlap(pixel_vertices)/norm;
	  }
	} else {
	  val = 1;
	  abs_xcoord = fabs(center_pixel_coords.x(*this));
	  abs_ycoord = fabs(center_pixel_coords.y(*this));
	  if(abs_xcoord>sqrt_three*abs_ycoord){
	    if(abs_xcoord>coord_lim) val = 0;
	  } else {
	    if(dot_product(tv, three_vector(abs_xcoord, abs_ycoord, 0, *this))>coord_lim)
	      val = 0;
	  }
	}

	if(val!=1){
	  if(interleaved_storage){
	    if(val==0){
	      wfdata[2*index] = wfdata[2*index+1] = 0;
	      continue;
	    }
	    wfdata[2*index] *= val;
	    if(real_imag_storage)
	      wfdata[2*index+1] *= val;
	  } else {
	    if(val==0){
	      wfdata[2*index] = wfdata[2*index+1] = 0;
	      continue;
	    }
	    wfdata[index] *= val;
	    if(real_imag_storage)
	      wfdata[index+nelem] *= val;
	  }
	}
      }
    }
  }

  double hexagonal_aperture::convex_polygon_overlap(const vector<three_point> & polygon_vertices) const { 
    try{
      double val;
      vector<three_point> hexagon_vertices(6);
      hexagon_vertices[0] = three_point(0, edge_length, 0, *this);
      hexagon_vertices[1] = three_point(-edge_length*sqrt(3.0)/2.0, edge_length/2.0, 0, *this);
      hexagon_vertices[2] = three_point(-edge_length*sqrt(3.0)/2.0, -edge_length/2.0, 0, *this);
      hexagon_vertices[3] = three_point(0, -edge_length, 0, *this);
      hexagon_vertices[4] = three_point(edge_length*sqrt(3.0)/2.0, -edge_length/2.0, 0, *this);
      hexagon_vertices[5] = three_point(edge_length*sqrt(3.0)/2.0, edge_length/2.0, 0, *this);
      
      vector<three_point> intersection_vertices = 
	get_convex_polygon_intersection(polygon_vertices, hexagon_vertices);
      if(intersection_vertices.size()<=2) val = 0;
      else val = get_area_of_polygon(intersection_vertices);
      return(val);
    } catch(...){
      cerr << "hexagonal_aperture::convex_polygon_overlap error\n";
      throw(string("hexagonal_aperture::convex_polygon_overlap"));
    }
  }

  namespace factory_register {
    const fits_keyval_set & get_tiled_hexagonal_aperture_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "tiled hexagonal aperture"));
      return *fkvs;
    }
    
    AO_sim_base * create_tiled_hexagonal_aperture(const iofits & iof) {
      return new tiled_hexagonal_aperture(iof);
    }
  }

  const bool tiled_hexagonal_aperture::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_tiled_hexagonal_aperture_keyval_set(), 
				      factory_register::create_tiled_hexagonal_aperture);

  namespace {

    void hex_to_rect(long jagged_index, 
		     long straight_index, 
		     double edge_length, 
		     double & jagged_coord, 
		     double & straight_coord){
      jagged_coord = 1.5*jagged_index*edge_length;
      straight_coord = jagged_index%2 ? 
	sqrt(3.0)*edge_length*(straight_index-.5) :
	sqrt(3.0)*edge_length*straight_index;

    }


    void hex_coords(double jagged_coord, 
		    double straight_coord,
		    double edge_length,
		    long & h, long & h_base,
		    long & e, long & e_base,
		    long & x, long & x_base){

      edge_length *= sqrt(3.0)/2.0;

      h = (long)floor(straight_coord/edge_length);
      e = (long)floor(.5*(sqrt(3.0)*jagged_coord + straight_coord)/edge_length);
      x = (long)floor(.5*(sqrt(3.0)*jagged_coord - straight_coord)/edge_length);

      h_base = -1;
      e_base = 0;
      long s1 = e-h-x;
      long s2 = (long)fmod((fmod(e-2*h,3.0)+3),3.0);
      if((s1==0 && s2!=0) || (s1==1 && s2==1)) e_base = -1;
      if(s2==0 || (s2==2 && s1==0)) h_base = 0;
      x_base = e_base-h_base-s1;
    }

    void rect_to_hex(double jagged_coord, 
		     double straight_coord,
		     double edge_length, 
		     long & jagged_index, 
		     long & straight_index,
		     long & nearest_jagged_index,
		     long & nearest_straight_index,
		     long & second_nearest_jagged_index,
		     long & second_nearest_straight_index){

      long h, e, x;
      long h_base, e_base, x_base;
      hex_coords(jagged_coord, straight_coord, edge_length,
		 h, h_base,
		 e, e_base,
		 x, x_base);

      int n = (e-2*h - e_base + 2*h_base)/3;
      int k = h - h_base;

      jagged_index = 2*n+k;
      straight_index = fmod(k,2.0)==0 ? k/2 : (k+1)/2;

      double hex_center_jagged_coord, hex_center_straight_coord;
      hex_to_rect(jagged_index, straight_index, edge_length,
		  hex_center_jagged_coord, hex_center_straight_coord);

      jagged_coord -= hex_center_jagged_coord;
      straight_coord -= hex_center_straight_coord;
      

      if(h_base==0 && e_base==0 && x_base==0){
	nearest_jagged_index = jagged_index+1;
	if(fmod(jagged_index,2.0)==0) nearest_straight_index = straight_index+1;
	else nearest_straight_index = straight_index;
	if(jagged_coord>straight_coord){
	  second_nearest_jagged_index = nearest_jagged_index;
	  second_nearest_straight_index = nearest_straight_index-1;
	} else {
	  second_nearest_jagged_index = jagged_index;
	  second_nearest_straight_index = straight_index+1;
	}
      } else if(h_base==0 && e_base==0 && x_base==-1){
	nearest_jagged_index = jagged_index;
	nearest_straight_index = straight_index+1;
	if(jagged_coord>0) second_nearest_jagged_index = jagged_index+1;
	else second_nearest_jagged_index = jagged_index-1;
	if(jagged_index%2==0) second_nearest_straight_index = nearest_straight_index;
	else second_nearest_straight_index = straight_index;
      } else if(h_base==0 && e_base==-1 && x_base==-1){
	nearest_jagged_index = jagged_index-1;
	if(fmod(jagged_index,2.0)==0) nearest_straight_index = straight_index+1;
	else nearest_straight_index = straight_index;	
	if(-jagged_coord>straight_coord){
	  second_nearest_jagged_index = nearest_jagged_index;
	  second_nearest_straight_index = nearest_straight_index-1;
	} else {
	  second_nearest_jagged_index = jagged_index;
	  second_nearest_straight_index = straight_index+1;
	}
      } else if(h_base==-1 && e_base==-1 && x_base==-1){
	nearest_jagged_index = jagged_index-1;
	if(fmod(jagged_index,2.0)==0) nearest_straight_index = straight_index;
	else nearest_straight_index = straight_index-1;
	if(jagged_coord<straight_coord){
	  second_nearest_jagged_index = nearest_jagged_index;
	  second_nearest_straight_index = nearest_straight_index+1;
	} else {
	  second_nearest_jagged_index = jagged_index;
	  second_nearest_straight_index = straight_index-1;
	}
      } else if(h_base==-1 && e_base==-1 && x_base==0){
	nearest_jagged_index = jagged_index;
	nearest_straight_index = straight_index-1;
	if(jagged_coord>0) second_nearest_jagged_index = jagged_index+1;
	else second_nearest_jagged_index = jagged_index-1;
	if(jagged_index%2==0) second_nearest_straight_index = straight_index;
	else second_nearest_straight_index = nearest_straight_index;
     } else if(h_base==-1 && e_base==0 && x_base==0){
	nearest_jagged_index = jagged_index+1;
	if(fmod(jagged_index,2.0)==0) nearest_straight_index = straight_index;
	else nearest_straight_index = straight_index-1;	
	if(jagged_coord>-straight_coord){
	  second_nearest_jagged_index = nearest_jagged_index;
	  second_nearest_straight_index = nearest_straight_index+1;
	} else {
	  second_nearest_jagged_index = jagged_index;
	  second_nearest_straight_index = straight_index-1;
	}
      } else {
	cerr << "rect_to_hex error - unexpected case\n";
	cerr << "h_base " << h_base << "\te_base " << e_base << "\tx_base " << x_base << endl;
	throw(string("rect_to_hex"));
      }
    }
  }

  tiled_hexagonal_aperture::tiled_hexagonal_aperture(const tiled_hexagonal_aperture & tiled_hex_ap){
    this->operator=(tiled_hex_ap);
  }

  tiled_hexagonal_aperture::tiled_hexagonal_aperture(const char * filename){
    this->read(filename);
  }

  tiled_hexagonal_aperture::tiled_hexagonal_aperture(const iofits & iof){
    this->read(iof);
  }

  tiled_hexagonal_aperture::tiled_hexagonal_aperture(double inner_diameter, double outer_diameter,
						     double in_edge_length, double in_gap_size) {

    if(inner_diameter<0 || outer_diameter<=0 || outer_diameter < inner_diameter){
      cerr << "tiled_hexagonal_aperture::tiled_hexagonal_aperture error - "
	   << "outer diameter must be positive, inner diameter must be nonnegative, "
	   << "and the former must be larger than the latter\n";
      cerr << "outer diameter " << outer_diameter << endl;
      cerr << "inner diameter " << inner_diameter << endl;
      throw(string("tiled_hexagonal_aperture::tiled_hexagonal_aperture"));
    }

    if(in_gap_size<0){
      cerr << "tiled_hexagonal_aperture::tiled_hexagonal_aperture error - "
	   << " gap size " << in_gap_size << " less than zero\n";
      throw(string("tiled_hexagonal_aperture::tiled_hexagonal_aperture"));
    }

    if(in_edge_length<0){
      cerr << "tiled_hexagonal_aperture::tiled_hexagonal_aperture error - "
	   << " edge length " << in_edge_length << " less than zero\n";
      throw(string("tiled_hexagonal_aperture::tiled_hexagonal_aperture"));
    }

    gap_size = in_gap_size;
    edge_length = in_edge_length;

    //long dimen = 2*(long)floor(.5*outer_diameter/(sqrt(3.0)*(edge_length+gap_size)))+1;
    long dimen = 2*(long)ceil(.5*outer_diameter/(sqrt(3.0)*(edge_length+2*gap_size/sqrt(3.0))))+3;
    tilemap = pixel_array<long>(vector<long>(2,dimen));
    
    three_point nearest_hex_center_coords;
    double edge_length_with_gaps = edge_length + 2*gap_size/sqrt(3.0);
    double jagged_coord, straight_coord, dist_from_center;
    for(int i=-dimen/2; i<=dimen/2; i++){
      for(int j=-dimen/2; j<=dimen/2; j++){
	hex_to_rect(i, j, edge_length_with_gaps, jagged_coord, straight_coord);
	dist_from_center = sqrt(jagged_coord*jagged_coord+straight_coord*straight_coord);
	if(dist_from_center<.5*outer_diameter && dist_from_center>=.5*inner_diameter)
	  tilemap.set_data((i+dimen/2)*dimen+j+dimen/2,1);
      }
    } 
  } 

  tiled_hexagonal_aperture & tiled_hexagonal_aperture::operator=(const tiled_hexagonal_aperture & tiled_hex_ap){
    if(this==&tiled_hex_ap)
      return(*this);

    tilemap = tiled_hex_ap.tilemap;
    gap_size = tiled_hex_ap.gap_size;
    edge_length = tiled_hex_ap.edge_length;

    this->aperture::operator=(tiled_hex_ap);
    return(*this);
  }

  void tiled_hexagonal_aperture::read(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "tiled_hexagonal_aperture::read - "
	   << "error opening file " << filename << endl;
      throw(string("tiled_hexagonal_aperture::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "tiled_hexagonal_aperture::read - "
	   << "error reading "
	   << this->unique_name() << " from file "
	   << filename << endl;
      throw(string("tiled_hexagonal_aperture::read"));
    }
  }

  void tiled_hexagonal_aperture::read(const iofits & iof) {
    if(!iof.key_exists("TYPE")){
      cerr << "tiled_hexagonal_aperture::read error - "
	   << "unrecognized type of file\n";
      throw(string("tiled_hexagonal_aperture::read"));
    }
    string type, comment;
    iof.read_key("TYPE", type, comment);
    if(type!=this->unique_name()){
      cerr << "tiled_hexagonal_aperture::read error - file of type " 
	   << type << " rather than of type "
	   << this->unique_name() << endl;
      throw(string("tiled_hexagonal_aperture::read"));
    }

    iof.read_key("AREALWTG", areal_weighting, comment);
    this->plane_optic::read(iof);
    iof.read_key("GAPSIZE", gap_size, comment);
    iof.read_key("EDGELNTH", edge_length, comment);
    tilemap.read(iof);

    // leave iof pointing to the next header, if it exists
    if(iof.get_hdu_num()!=iof.get_num_hdus())
      iof.movrel_hdu(1);
  }

  void tiled_hexagonal_aperture::write(const char * filename) const {

    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "tiled_hexagonal_aperture::write - "
	   << "error opening file " << filename << endl;
      throw(string("tiled_hexagonal_aperture::write"));
    }

    try{this->write(iof);}
    catch(...){
      cerr << "tiled_hexagonal_aperture::write - "
	   << "error writing "
	   << this->unique_name() << " to file " 
	   << filename << endl;
      throw(string("tiled_hexagonal_aperture::write"));
    }
  }

  void tiled_hexagonal_aperture::write(iofits & iof) const {
    fits_header_data<long> tmphdr(tilemap.get_axes());
    tmphdr.write(iof);

    string type = this->unique_name();
    string comment = "object type";
    iof.write_key("TYPE", type, comment);
    iof.write_key("GAPSIZE", gap_size, "gap size in meters");
    iof.write_key("EDGELNTH", edge_length, "edge length in meters");
    iof.write_key("AREALWTG", areal_weighting, "weight edge pixels by area of overlap");
    this->plane_optic::write(iof);
    tilemap.write(iof);
  }

  void tiled_hexagonal_aperture::print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << this->unique_name()
       << "/" << "object type" << endl;
    this->aperture::print(os, prefix);
    os << prefix << "GAPSIZE    = " << setw(vlspc) << gap_size
       << "/" << "gap size in meters" << endl;
    os << prefix << "EDGELNTH   = " << setw(vlspc) << edge_length
       << "/" << "hexagon edge length (meters)" << endl;
  }

  rectangular_region tiled_hexagonal_aperture::get_covering_region(const three_frame & tf) const {
    vector<long> axes = tilemap.get_axes();
    int dimen = axes[0];
    hexagonal_aperture hex_ap(edge_length);
    hex_ap.three_frame::operator=(*this);

    rectangular_region covering_region(hex_ap.get_covering_region(*this));
    double edge_length_with_gaps = edge_length + 2*gap_size/sqrt(3.0);
    double jagged_coord, straight_coord;

    for(int i=-dimen/2; i<=dimen/2; i++){
      for(int j=-dimen/2; j<=dimen/2; j++){
	if(tilemap.data((i+dimen/2)*dimen+j+dimen/2)==1){
	  hex_to_rect(i, j, edge_length_with_gaps, jagged_coord, straight_coord);
	  hex_ap.three_point::operator=(three_point(straight_coord, jagged_coord, 0, *this));
	  covering_region = region_union(covering_region, hex_ap.get_covering_region(*this));
	}
      }
    }
    return covering_region;
  }

  void tiled_hexagonal_aperture::transform(diffractive_wavefront<float> & wf) const {
    try{this->private_transform(wf);}
    catch(...){
      cerr << "tiled_hexagonal_aperture::transform error - "
	   << "error transforming a float instantiation of diffractive_wavefront\n";
      throw(string("tiled_hexagonal_aperture::transform"));
    }
  }

  void tiled_hexagonal_aperture::transform(diffractive_wavefront<double> & wf) const {
    try{this->private_transform(wf);}
    catch(...){
      cerr << "tiled_hexagonal_aperture::transform error - "
	   << "error transforming a double instantiation of diffractive_wavefront\n";
      throw(string("tiled_hexagonal_aperture::transform"));
    }
  }

  template<class T>
  void tiled_hexagonal_aperture::private_transform(diffractive_wavefront<T> & wf) const { 

    // get projected wf pixel spacing
    three_vector origin_offset, dx, dy;
    this->get_projected_wavefront_pixel_spacing(wf, origin_offset, dx, dy);

    // The halfpixel information
    vector<long> wf_axes = wf.get_axes();
    double x_halfpix=0, y_halfpix=0;
    int x_extrapix=1, y_extrapix=1;
    if(wf_axes[1]%2==0){
      x_halfpix = .5;
      x_extrapix = 0;
    }
    if(wf_axes[0]%2==0){
      y_halfpix = .5;
      y_extrapix = 0;
    }

    bool real_imag_storage = is_real_imag_storage(wf);
    bool interleaved_storage = is_interleaved_storage(wf);
    T * wfdata = get_wavefront_data(wf);

    long jagged_index, straight_index;
    long nearest_jagged_index, nearest_straight_index;
    long second_nearest_jagged_index, second_nearest_straight_index;
    double jagged_coord, straight_coord;
    double edge_length_with_gaps = edge_length + 2*gap_size/sqrt(3.0);
    vector<long> tilemap_axes = tilemap.get_axes();

    long dimen = tilemap.get_axes()[0];
    long dimen_lim = (dimen-1)/2;


    int index;
    int nelem = wf_axes[0]*wf_axes[1];
    three_point center_pixel_coords;
    three_vector tmps = .5*(dy+dx);
    three_vector tmpd = .5*(dy-dx);
    vector<three_point> pixel_vertices(4), intersection_vertices;
    double norm = cross_product(dx, dy).length();
    double buffer = tmps.length() > tmpd.length() ? 2*tmps.length() : 2*tmpd.length();
    double upper_hex_limit = (edge_length+gap_size+buffer)*(edge_length+gap_size+buffer);
    double lower_hex_limit = (.5*sqrt(3.0)*edge_length-buffer)*(.5*sqrt(3.0)*edge_length-buffer);
    double val;
    double hex_distance_squared;
    bool pixel_vertices_initialized;

    three_point hexagon_center_point;
    three_vector hexagon_offset_vector;
    vector<three_point> hexagon_vertices(6), base_hexagon_vertices(6);
    base_hexagon_vertices[0] = three_point(0, edge_length, 0, *this);
    base_hexagon_vertices[1] = three_point(-edge_length*.5*sqrt(3.0), .5*edge_length, 0, *this);
    base_hexagon_vertices[2] = three_point(-edge_length*.5*sqrt(3.0), -.5*edge_length, 0, *this);
    base_hexagon_vertices[3] = three_point(0, -edge_length, 0, *this);
    base_hexagon_vertices[4] = three_point(edge_length*.5*sqrt(3.0), -.5*edge_length, 0, *this);
    base_hexagon_vertices[5] = three_point(edge_length*.5*sqrt(3.0), .5*edge_length, 0, *this);

    vector<long> jagged_indices(3), straight_indices(3);

    for(int i=-wf_axes[1]/2; i<wf_axes[1]/2+x_extrapix; i++){
      for(int j=-wf_axes[0]/2; j<wf_axes[0]/2+y_extrapix; j++){

	index = (i+wf_axes[1]/2)*wf_axes[0]+j+wf_axes[0]/2;
	center_pixel_coords = *this +
	  (origin_offset + (i+x_halfpix)*dx + (j+y_halfpix)*dy);

	// Find the indices and center coords of the hexagonal aperture
	// containing this pixel and the neighboring hexagonal aperture.
	rect_to_hex(center_pixel_coords.y(*this), center_pixel_coords.x(*this),
		    edge_length_with_gaps, 
		    jagged_indices[0], straight_indices[0],
		    jagged_indices[1], straight_indices[1],
		    jagged_indices[2], straight_indices[2]);

	val = 0;

	pixel_vertices_initialized = false;

	for(int k=0; k<3; k++){

	  if(fabs(1-val)<three_frame::precision) break;

	  if(jagged_indices[k]>=-dimen_lim && jagged_indices[k]<=dimen_lim &&
	     straight_indices[k]>=-dimen_lim && straight_indices[k]<=dimen_lim &&
	     tilemap.data((jagged_indices[k]+dimen_lim)*dimen+straight_indices[k]+dimen_lim)!=0) {
	    
	    hex_to_rect(jagged_indices[k], straight_indices[k], 
			edge_length_with_gaps, 
			jagged_coord, straight_coord);

	    hexagon_center_point = three_point(straight_coord, jagged_coord, 0, *this);
	    
	    hex_distance_squared = (center_pixel_coords-hexagon_center_point).length_squared();
	    
	    if(hex_distance_squared<lower_hex_limit){
	      val += 1;
	    } else if(hex_distance_squared<upper_hex_limit){
	      
	      if(!pixel_vertices_initialized){
		pixel_vertices[0] = center_pixel_coords + tmps;
		pixel_vertices[1] = center_pixel_coords + tmpd;
		pixel_vertices[2] = center_pixel_coords - tmps;
		pixel_vertices[3] = center_pixel_coords - tmpd;
	      }

	      hexagon_offset_vector = hexagon_center_point - *this;
	      for(int k=0; k<6; k++)
		hexagon_vertices[k] = base_hexagon_vertices[k] + hexagon_offset_vector;
	      
	      intersection_vertices = get_convex_polygon_intersection(pixel_vertices, hexagon_vertices);
	      if(intersection_vertices.size()>2)
		val += get_area_of_polygon(intersection_vertices)/norm;
	    }
	  }
	}
	
	if(interleaved_storage){
	  if(val==0){
	    wfdata[2*index] = wfdata[2*index+1] = 0;
	    continue;
	  }
	  wfdata[2*index] *= val;
	  if(real_imag_storage)
	    wfdata[2*index+1] *= val;
	} else {
	  if(val==0){
	    wfdata[2*index] = wfdata[2*index+1] = 0;
	    continue;
	  }
	  wfdata[index] *= val;
	  if(real_imag_storage)
	    wfdata[index+nelem] *= val;
	}
      }
    }
  }

  double tiled_hexagonal_aperture::convex_polygon_overlap(const vector<three_point> & polygon_vertices) const {

    try{
      vector<three_point> hexagon_vertices(6), tmp_vertices(6), intersection_vertices;
      hexagon_vertices[0] = three_point(0, edge_length, 0, *this);
      hexagon_vertices[1] = three_point(-edge_length*sqrt(3.0)/2.0, edge_length/2.0, 0, *this);
      hexagon_vertices[2] = three_point(-edge_length*sqrt(3.0)/2.0, -edge_length/2.0, 0, *this);
      hexagon_vertices[3] = three_point(0, -edge_length, 0, *this);
      hexagon_vertices[4] = three_point(edge_length*sqrt(3.0)/2.0, -edge_length/2.0, 0, *this);
      hexagon_vertices[5] = three_point(edge_length*sqrt(3.0)/2.0, edge_length/2.0, 0, *this);
      
      vector<long> axes = tilemap.get_axes();
      int dimen = axes[0];
      hexagonal_aperture hex_ap(edge_length);
      hex_ap.three_frame::operator=(*this);
      
      rectangular_region covering_region(hex_ap.get_covering_region(*this));
      double edge_length_with_gaps = edge_length + 2*gap_size/sqrt(3.0);
      double jagged_coord, straight_coord;
      three_vector offset_vector;
      double val = 0;
      
      for(int i=-dimen/2; i<=dimen/2; i++){
	for(int j=-dimen/2; j<=dimen/2; j++){
	  if(this->tilemap.data((i+dimen/2)*dimen+j+dimen/2)==1){
	    hex_to_rect(i, j, edge_length_with_gaps, jagged_coord, straight_coord);
	    offset_vector = three_vector(straight_coord, jagged_coord, 0, *this);
	    for(int k=0; k<6; k++)
	      tmp_vertices[k] = hexagon_vertices[k] + offset_vector;
	    
	    intersection_vertices = 
	      get_convex_polygon_intersection(polygon_vertices, tmp_vertices);
	    if(intersection_vertices.size()>2)
	      val += get_area_of_polygon(intersection_vertices);
	  }
	}
      }
      return(val);
    } catch(...){
      cerr << "tiled_hexagonal_aperture::convex_polygon_overlap error\n";
      throw(string("tiled_hexagonal_aperture::convex_polygon_overlap"));
    }
  }
}
