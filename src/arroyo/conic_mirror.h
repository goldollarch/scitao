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

#ifndef CONIC_MIRROR_H
#define CONIC_MIRROR_H

#include <iomanip>
#include "optic.h"
#include "diffractive_wavefront.h"
#include "conic_section.h"

namespace Arroyo {


  class conic_mirror_base :
    public one_to_one_optic,
    public conic_section  {

    public:

    ///////////////////////////////////////////
    ///  Null constructor
    conic_mirror_base(){};

    ///////////////////////////////////////////
    ///  Virtual destructor
    virtual ~conic_mirror_base(){};

    ///////////////////////////////////////////
    ///  Read from file
    virtual void read(const char * filename) = 0;

    ///////////////////////////////////////////
    ///  Read from an iofits object
    virtual void read(const iofits & iof) = 0;
 
    ///////////////////////////////////////////
    ///  Write to file
    virtual void write(const char * filename) const = 0;

    ///////////////////////////////////////////
    ///  Write to an iofits object
    virtual void write(iofits & iof) const = 0;

    ///////////////////////////////////////////
    ///  Print
    virtual void print(ostream & os, const char * prefix="") const = 0;

    ///////////////////////////////////////////
    ///  Get the point of intersection of
    ///  a line extending from three_point tp
    ///  in the direction of the three_vector
    ///  tv and this optic.  If there is no
    ///  intersection point, this function throws
    ///  an error
    virtual three_point get_point_of_intersection(const three_point & tp, 
						  const three_vector & tv) const = 0;

  };


  ///
  /// A class to represent a mirror formed from a conic shape.
  /// If the mirror is constructed using an aperture with a
  /// three frame of non-common origin, this will model
  /// an off-axis mirror (e.g. an off-axis parabola)
  ///

  template<class aperture_type>
    class conic_mirror :
    public conic_mirror_base {

    private:

    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("conic mirror"));};

    protected:

    // The aperture
    aperture_type ap;

    // The transmission spectrum of this mirror
    //mirror_reflectivity_spectrum * mirror_reflec_spec;
    
    // Here we can either choose to enforce
    // the law of reflection, or we can 
    // require that the raytrace occurs along
    // a ray perpendicular to the wavefront 
    // surface.
    //
    // raytrace_policy=true => law of reflection
    // raytrace_policy=false => orthogonal final ray
    bool raytrace_policy;

    ///////////////////////////////////////////
    ///  Null constructor
    ///
    ///  Protected because a class instance requires
    ///  initialization
    conic_mirror(){};

    ///////////////////////////////////////////
    /// A template member function to perform
    /// transform on both float and double
    /// instantiations of wavefront.  This 
    /// is necessary because there is no 
    /// mechanism in C++ for virtual template
    /// member functions.
    template<class T>
      void private_transform(diffractive_wavefront<T> & wf) const;

    public:

    ///////////////////////////////////////////
    ///  Copy constructor
    conic_mirror(const conic_mirror & cmr);

    ///////////////////////////////////////////
    ///  Construct from file
    conic_mirror(const char * filename);

    ///////////////////////////////////////////
    ///  Construct from iofits object
    conic_mirror(const iofits & iof);

    ///////////////////////////////////////////
    ///  Construct from the bits
    /// 
    ///  The vertex of the conic lies along the 
    ///  axis of symmetry. 
    ///
    ///  The focus of the conic is the focal point
    ///  closest to the vertex.  
    ///
    ///  The eccentricity of the conic is as follows:
    ///
    ///     For a spherical mirror, eccentricity = 0
    ///     For an elliptical mirror, 0 < eccentricity < 1
    ///     For a parabolic mirror, eccentricity = 1
    ///     For a hyperbolic mirror, eccentricity > 1
    ///
    /// The concave flag denotes the reflective surface of the conic.
    /// If the flag is true, the conic reflects waves incident from
    /// the general direction of the focal point.
    conic_mirror(const three_point & vertex, 
		 const three_point & focus, 
		 double eccty, 
		 const aperture_type & in_ap);

    ///////////////////////////////////////////
    ///  Destructor
    ~conic_mirror(){};

    ///////////////////////////////////////////
    ///  Operator = 
    conic_mirror & operator=(const conic_mirror & cmr);

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
    ///  Get the raytrace policy
    ///
    ///  Returns true if the policy is to use the law of reflection
    ///
    ///  Returns false if the policy is to force the final ray to be 
    ///  orthogonal to the final spherical wavefront
    ///
    ///  These differ due to numerical precision
    bool get_raytrace_policy() const {
      return(this->raytrace_policy);
    };

    ///////////////////////////////////////////
    ///  Set the raytrace policy
    void set_raytrace_policy(bool raytrace_policy){
      this->raytrace_policy = raytrace_policy;
    };

    ///////////////////////////////////////////
    ///  Get a rectangular region guaranteed to cover
    ///  the aperture.  The resulting region will
    ///  have its edges aligned with the x and y axes
    ///  of the three_frame tf.  
    ///
    ///  If foreshortening is on, the projected
    ///  region is guaranteed to cover the optic
    ///
    ///  If the z axis of the three_frame is orthogonal
    ///  to the z axis of the aperture, this function
    ///  throws an error
    rectangular_region get_covering_region(const three_frame & tf) const;

    ///////////////////////////////////////////
    ///  Get the point of intersection of
    ///  a line extending from three_point tp
    ///  in the direction of the three_vector
    ///  tv and the conic section.  If there is no
    ///  intersection point, this function throws
    ///  an error
    three_point get_point_of_intersection(const three_point & tp, 
					  const three_vector & tv) const;
    
    ///////////////////////////////////////////
    ///  Apply the aperture to the wavefront
    void transform(diffractive_wavefront<float> & wf) const;

    ///////////////////////////////////////////
    ///  Apply the aperture to the wavefront
    void transform(diffractive_wavefront<double> & wf) const;

  };


  template<class aperture_type>
    conic_mirror<aperture_type>::conic_mirror(const conic_mirror & cmr) {
    this->operator=(cmr);
  }

  template<class aperture_type>
    conic_mirror<aperture_type>::conic_mirror(const char * filename){
    this->read(filename);
  }

  template<class aperture_type>
    conic_mirror<aperture_type>::conic_mirror(const iofits & iof){
    this->read(iof);
  }

  template<class aperture_type>
    conic_mirror<aperture_type>::conic_mirror(const three_point & vtx, 
					      const three_point & focus, 
					      double eccty, 
					      const aperture_type & in_ap) : ap(in_ap){

    this->conic_section::operator=(conic_section(vtx,focus,eccty));

    this->raytrace_policy = true;

    if((this->get_point_of_intersection(ap, ap.z()) - ap).length()>three_frame::precision){
      cerr << "conic_mirror::conic_mirror error - "
	   << "aperture center does not lie on the surface of the conic mirror\n";
      throw(string("conic_mirror::conic_mirror"));
    }
  }
 
  template<class aperture_type>
    conic_mirror<aperture_type> & conic_mirror<aperture_type>::operator=(const conic_mirror & cmr) {
    if(this==&cmr) return(*this);
    this->conic_section::operator=(cmr);
    this->ap = cmr.ap;
    this->raytrace_policy=cmr.raytrace_policy;
    return(*this);
  }

  template<class aperture_type>
    void conic_mirror<aperture_type>::read(const char * filename) {
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "conic_mirror::read - "
	   << "error opening file " << filename << endl;
      throw(string("conic_mirror::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "conic_mirror::read - "
	   << "error reading "
	   << this->unique_name() << " from file " 
	   << filename << endl;
      throw(string("conic_mirror::read"));
    }
  }

  template<class aperture_type>
    void conic_mirror<aperture_type>::read(const iofits & iof) {
    if(!iof.key_exists("TYPE")){
      cerr << "conic_mirror::read error - "
	   << "no power law type specified\n";
      throw(string("conic_mirror::read"));
    } else {
      string type, comment;
      iof.read_key("TYPE", type, comment);
      if(type!=this->unique_name()){
	cerr << "conic_mirror::read error - outer scale of type " 
	     << type << " rather than type "
	     << this->unique_name() << endl;
	throw(string("conic_mirror::read"));
      } else {

	iof.read_key("RAYPLCY", this->raytrace_policy, comment);

	// move to the aperture header
	try{
	  iof.movrel_hdu(1);
	} catch(...){
	  cerr << "conic_mirror::read - error moving to aperture hdu\n";
	  throw(string("conic_mirror::read"));
	}	  

	try{
	  this->ap.read(iof);
	} catch(...){
	  cerr << "conic_mirror::read - error reading aperture\n";
	  throw(string("conic_mirror::read"));
	}	  

	try{
	  this->conic_section::read(iof);
	} catch(...){
	  cerr << "conic_mirror::read - error reading conic section\n";
	  throw(string("conic_mirror::read"));
	}	
      }
    }

  }
 
  template<class aperture_type>
    void conic_mirror<aperture_type>::write(const char * filename) const {
    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "conic_mirror::write - "
	   << "error opening file " << filename << endl;
      throw(string("conic_mirror::write"));
    }

    try{this->write(iof);}
    catch(...){
      cerr << "conic_mirror::write - "
	   << "error writing "
	   << this->unique_name() << " to file " 
	   << filename << endl;
      throw(string("conic_mirror::write"));
    }
  }

  template<class aperture_type>
    void conic_mirror<aperture_type>::write(iofits & iof) const {
    Arroyo::fits_header_data<char> tmphdr;
    tmphdr.write(iof);

    string type = this->unique_name();
    string comment = "object type";
    iof.write_key("TYPE", type, comment);
    iof.write_key("RAYPLCY", this->raytrace_policy, string("raytrace policy"));

    this->ap.write(iof);
    this->conic_section::write(iof);
  }

  template<class aperture_type>
    void conic_mirror<aperture_type>::print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(std::ios::left, std::ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << this->unique_name()
       << "/" << "object type" << endl;
    os << prefix << "RAYPLCY    = " << setw(vlspc) << this->raytrace_policy
       << "/" << "raytrace policy" << endl;
    this->conic_section::print(os, prefix);
    this->ap.print(os, prefix);
  }

  template<class aperture_type>
    rectangular_region conic_mirror<aperture_type>::get_covering_region(const three_frame & tf) const {
    return(this->ap.get_covering_region(tf));
  }

  template<class aperture_type>
    three_point conic_mirror<aperture_type>::get_point_of_intersection(const three_point & tp,
								       const three_vector & tv) const {

    double distance_to_conic, R_squared, V;
    three_vector conic_unit_normal;
    three_point point_of_intersection;

    this->conic_section::raytrace(tp,
				  tv,
				  distance_to_conic,
				  point_of_intersection,
				  conic_unit_normal,
				  R_squared,
				  V);

    return(point_of_intersection);
  }

  template<class aperture_type>
    void conic_mirror<aperture_type>::transform(diffractive_wavefront<float> & wf) const {
    try{this->private_transform(wf);}
    catch(...){
      cerr << "conic_mirror::transform error - "
	   << "error transforming a float instantiation of diffractive_wavefront\n";
      throw(string("conic_mirror::transform"));
    }
  }

  template<class aperture_type>
    void conic_mirror<aperture_type>::transform(diffractive_wavefront<double> & wf) const {
    try{this->private_transform(wf);}
    catch(...){
      cerr << "conic_mirror::transform error - "
	   << "error transforming a double instantiation of diffractive_wavefront\n";
      throw(string("conic_mirror::transform"));
    }
  }
  
  template<class aperture_type>
    template<class T>
    void conic_mirror<aperture_type>::private_transform(diffractive_wavefront<T> & wf) const {

    // Here we do this off the bat to avoid sign convention issues in
    // the raytrace that I can't figure out right now.
    //wf.set_wavefront_curvature(0);

    //  First, test whether the center of the wavefront 
    //  lies on the conic surface
    double init_wf_to_conic_OPD, conic_to_final_wf_OPD, conic_R_squared, conic_V;
    three_vector conic_unit_normal;
    three_point conic_point_of_intersection;    

    try{
      this->conic_section::raytrace(wf,
				    wf.z(),
				    init_wf_to_conic_OPD,
				    conic_point_of_intersection,
				    conic_unit_normal,
				    conic_R_squared,
				    conic_V);
    } catch(...) {
      cerr << "conic_mirror::private_transform - error tracing central ray\n";
      throw(string("conic_mirror::private_transform"));
    }

    if(fabs(init_wf_to_conic_OPD)>three_frame::precision){
      cerr << "conic_mirror::private_transform error - wavefront does not lie on conic surface\n";
      cerr << "distance to conic " << init_wf_to_conic_OPD << endl;

      throw(string("conic_mirror::private_transform"));
    }
    
    // Define the three_frame of the final wavefront using the law of
    // reflection: rotate the incident wavefront three_frame about the
    // vector orthogonal to the plane formed from wf.z() and the
    // conic_unit_normal, and then flip the sign on wf.z()
    three_frame final_wf_tf(wf);
    three_vector perpendicular_vector;

    try{
      perpendicular_vector = cross_product(wf.z(),
					   conic_unit_normal);
      
      if(perpendicular_vector.length()>three_frame::precision){
	three_rotation trot(wf, 
			  perpendicular_vector, 
			    -2*perpendicular_vector.length());
	trot.transform(final_wf_tf);
      }
      final_wf_tf = three_frame(final_wf_tf, 
				final_wf_tf.x(),
				final_wf_tf.y(),
				-1*final_wf_tf.z());
    } catch(...) {
      cerr << "conic_mirror::private_transform - error forming final three frame\n";
      throw(string("conic_mirror::private_transform"));
    }

    // Define the final wavefront curvature based on the paraxial approximation
    // and the local curvature.
    //
    // Here the wavefront curvature is positive for diverging wavefronts,
    // and the conic curvature is positive for concave conics.
    //
    // Thus, a diverging wavefront meeting a concave conic with half
    // the curvature will yield a flat wavefront, and a flat wavefront
    // reflecting off a concave conic will yield a converging wavefront
    double initial_wavefront_curvature = wf.get_curvature();
    double final_wavefront_curvature = 
      initial_wavefront_curvature - 
      2*this->conic_section::get_local_curvature(wf);

    if(optic::verbose_level)
      cout << "conic_mirror::private_transform: "
	   << "initial curvature " << initial_wavefront_curvature
	   << " local curvature " << this->conic_section::get_local_curvature(wf)
	   << " final curvature " << final_wavefront_curvature
	   << endl;

    // Define a conic section to represent the final spherical wavefront
    three_point final_focus;
    conic_section final_wavefront_conic;
    if(fabs(final_wavefront_curvature)<three_frame::precision){
      final_focus = final_wf_tf + 1e15*wf.z();
      final_wavefront_conic = conic_section(wf,
					    final_focus,
					    0);
    } else {
      final_focus = final_wf_tf - (1/final_wavefront_curvature)*final_wf_tf.z();
      final_wavefront_conic = conic_section(wf,
					    final_focus,
					    0);
    }

    try{

      // Make the transformation to amp, phase storage in the wf, and
      // get the pointer to the raw wavefront data
      bool interleaved_storage = is_interleaved_storage(wf);
      this->wavefront_amp_phase_conversion(wf);
      T * wfdata = get_wavefront_data(wf);

      // The halfpixel information
      vector<long> wf_axes = wf.get_axes();
      long nelem = wf_axes[0]*wf_axes[1];
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
      
      int index;
      double residual_OPD;
      double wavelength_meters = wf.get_wavelength();
      double arc_length_meters;
      double twopi = 2*M_PI;
      double final_wf_R_squared, final_wf_V;
      
      three_vector rotation_vector;
      three_point initial_pixel_tp, final_pixel_tp;
      three_vector initial_pixel_tv, final_pixel_tv;
      three_vector perpendicular_vector;
      three_point final_wf_point_of_intersection;
      three_vector final_wf_unit_normal;

      for(int i=-wf_axes[1]/2; i<wf_axes[1]/2+x_extrapix; i++){
	for(int j=-wf_axes[0]/2; j<wf_axes[0]/2+y_extrapix; j++){

	  try{
	    if(fabs(initial_wavefront_curvature)<three_frame::precision){
	      
	      initial_pixel_tv = wf.z();
	      initial_pixel_tp = 
		wf + 
		(i+x_halfpix)*wf.get_pixel_scale()*wf.x()+
		(j+y_halfpix)*wf.get_pixel_scale()*wf.y();
	      
	    } else {
	      
	      initial_pixel_tv = wf.z();

	      three_point wavefront_center_of_curvature = 
		wf - fabs(1/initial_wavefront_curvature)*wf.z();
	      
	      // the length of the arc
	      arc_length_meters = 
		sqrt((i+x_halfpix)*(i+x_halfpix)+
		     (j+y_halfpix)*(j+y_halfpix))*
		wf.get_pixel_scale();
	      
	      // the vector about which to rotate wf.z() to
	      // get the ray propagation direction

	      if(fabs(arc_length_meters)>three_frame::precision){
		rotation_vector = cross_product(wf.z(),
						(i+x_halfpix)*wf.x()+
						(j+y_halfpix)*wf.y());
		
		three_rotation trot(wavefront_center_of_curvature,
				    rotation_vector,
				    -1*initial_wavefront_curvature*arc_length_meters);
		
		trot.transform(initial_pixel_tv);
	      }

	      initial_pixel_tp = 
		wavefront_center_of_curvature +
		(1/fabs(initial_wavefront_curvature))*initial_pixel_tv;
	    }
	  } catch(...) {
	    cerr << "conic_mirror::private_transform - error finding initial point on wavefront\n";
	    throw(string("conic_mirror::private_transform"));
	  }


	  // If the initial point is on the conic, we can skip the
	  // analysis
	  if(this->conic_section::on_conic(initial_pixel_tp))
	    continue;

	  residual_OPD = 0;

	  try{
	    this->conic_section::raytrace(initial_pixel_tp,
					  initial_pixel_tv,
					  init_wf_to_conic_OPD,
					  conic_point_of_intersection,
					  conic_unit_normal,
					  conic_R_squared,
					  conic_V);

	    residual_OPD += init_wf_to_conic_OPD;
	    final_pixel_tp = conic_point_of_intersection;

	  } catch(...) {
	    cerr << "conic_mirror::private_transform - error tracing initial rays\n";
	    throw(string("conic_mirror::private_transform"));
	  }


	  try{
	    if(this->raytrace_policy){
	      
	      final_pixel_tv = initial_pixel_tv;

	      three_reflection tref(conic_point_of_intersection,
				    conic_unit_normal);
	      
	      tref.transform(final_pixel_tv);

	      // Make sure we've reflected properly
	      if(fabs(dot_product(initial_pixel_tv,
				  conic_unit_normal) +
		      dot_product(final_pixel_tv,
				  conic_unit_normal))>three_frame::precision){
		cerr << "conic_mirror::private_transform - error forming final ray vector\n";
		
		three_vector tmp = -1*initial_pixel_tv;

		tmp.print(cerr, "\treflected initial tv ");
		conic_unit_normal.print(cerr, "\tconic normal tv ");
		perpendicular_vector.print(cerr, "\tperp tv ");
		final_pixel_tv.print(cerr, "\tfinal tv ");
		
		cout << "\tlengths: init " 
		     << initial_pixel_tv.length()
		     << "\treflected "
		     << tmp.length()
		     << "\tconic unit normal " 
		     << conic_unit_normal.length()
		     << "\tfinal "
		     << final_pixel_tv.length()
		     << endl;

		cout << "\tdot products: " 
		     << dot_product(tmp,conic_unit_normal) 
		     << "\t" 
		     << dot_product(conic_unit_normal, final_pixel_tv) 
		     << "\t" 
		     << dot_product(tmp,conic_unit_normal) - dot_product(conic_unit_normal, final_pixel_tv) 
		     << endl;

		throw(string("conic_mirror::private_transform"));
	      }

	    } else {
	      // Requirement that the ray is perpendicular to the
	      // spherical wavefront surface
	      final_pixel_tv = final_focus - conic_point_of_intersection;
	      final_pixel_tv = (1/final_pixel_tv.length())*final_pixel_tv;
	    }
	  } catch (...){ 
	    cerr << "conic_mirror::private_transform - error finding final pixel vector\n";
	    throw(string("conic_mirror::private_transform"));
	  }

	  try{
	    final_wavefront_conic.raytrace(final_pixel_tp,
					   final_pixel_tv,
					   conic_to_final_wf_OPD,
					   final_wf_point_of_intersection,
					   final_wf_unit_normal,
					   final_wf_R_squared,
					   final_wf_V);

	    residual_OPD += conic_to_final_wf_OPD;
	  } catch(...) {
	    cerr << "conic_mirror::private_transform - error tracing final rays\n";
	    final_pixel_tv.print(cerr, "fpv ");
	    final_focus.print(cerr, "ffocus ");
	    final_wf_point_of_intersection.print(cerr, "cnc pt int");
	    throw(string("conic_mirror::private_transform"));
	  }


	  if(optic::verbose_level){
	    cout << "conic_mirror::private_transform - "
		 << (i+x_halfpix)*wf.get_pixel_scale()
		 << "\t" 
		 << (j+y_halfpix)*wf.get_pixel_scale()
		 << "\tinitial OPD " 
		 << init_wf_to_conic_OPD
		 << "\tfinal OPD " 
		 << conic_to_final_wf_OPD
		 << " net OPD "
		 << residual_OPD
		 << " phase " 
		 << twopi*fmod(residual_OPD/wavelength_meters, 1.)
		 << endl;
	  }

	  // Some tests:
	  if(cross_product(.5*(final_pixel_tv - initial_pixel_tv),conic_unit_normal).length()/1e6>
	     three_frame::precision){
	    cout << "conic_mirror::private_transform error - cross product "
		 << cross_product(.5*(final_pixel_tv - initial_pixel_tv),conic_unit_normal).length()
		 << " indicates initial, final, and conic normal vectors are not coplanar\n";
	    
	    initial_pixel_tp.print(cout, "\tinitial pixel tp ");
	    initial_pixel_tv.print(cout, "\tinitial pixel tv ");
	    final_pixel_tp.print(cout, "\tfinal pixel tp ");
	    final_pixel_tv.print(cout, "\tfinal pixel tv ");
	    conic_point_of_intersection.print(cout, "\tconic pt of intersection ");
	    conic_unit_normal.print(cout, "\tconic unit normal ");
	    final_wf_point_of_intersection.print(cout, "\tfinal wf pt of intersection ");
	    final_wf_unit_normal.print(cout, "\tfinal wf unit normal ");
	    cout << "\tconic params:  R_squared " << conic_R_squared << "\tV " << conic_V << endl;
	    cout << "\tfinal wf params:  R_squared " << final_wf_R_squared << "\tV " << final_wf_V << endl;
	    cout << endl << endl;
	    three_vector tmp = .5*(final_pixel_tv - initial_pixel_tv);

	    tmp.print(cout, "\t.5*(final - initial tv) ");
	    tmp = (1/tmp.length())*tmp;
	    tmp.print(cout, "\tnormalized .5*(final - initial tv) ");
	    conic_unit_normal.print(cout, "\tconic unit normal ");
	    (tmp-conic_unit_normal).print(cout, "\tdiff vector ");

	    cross_product(.5*(final_pixel_tv - initial_pixel_tv),conic_unit_normal).print(cout, "\txproduct " );
	    cout << "\tInit length " << initial_pixel_tv.length() << endl;
	    cout << "\tFinal length " << final_pixel_tv.length() << endl;
	    cout << "\tDiff length " << .5*(final_pixel_tv - initial_pixel_tv).length() << endl;
	    cout << "\tConic normal length " << conic_unit_normal.length() << endl;
	    
	    throw(string("conic_mirror::private_transform"));
	  }


	  // Set the new phase correctly
	  //
	  // A positive residual OPD indicates that the ray is advanced
	  // wrt the wavefront, so that OPD/wavelength should be
	  // added to the wavefront phase
	  try{
	    index = (i+wf_axes[1]/2)*wf_axes[0]+j+wf_axes[0]/2;
	    if(interleaved_storage){
	      wfdata[2*index+1] = 
		fmod(wfdata[2*index+1]+twopi*residual_OPD/wavelength_meters, twopi);
	      if(wfdata[2*index+1]<0) wfdata[2*index+1]+=twopi;
	    } else {
	      wfdata[index+nelem] = 
		fmod(wfdata[index+nelem]+twopi*residual_OPD/wavelength_meters, twopi);
	      if(wfdata[index+nelem]<0) wfdata[index+nelem]+=twopi;
	    }
	  } catch(...) {
	    cerr << "conic_mirror::private_transform - error updating phases\n";
	    throw(string("conic_mirror::private_transform"));
	  }

	    
	  // set the curvature
	  wf.set_curvature(final_wavefront_curvature);

	}
      }
    } catch(...) {
      cerr << "conic_mirror::private_transform - error raytracing wavefront\n";
      throw(string("conic_mirror::private_transform"));
    }

    // Do the aperture transformation.  Translate the aperture
    // along ap->z() to the transverse plane of the wavefront.  Apply
    // the aperture, and translate it back
    try{
      three_point wf_ap_intsctn_pt = 
	ap.get_point_of_intersection(wf,
				     ap.z());

      three_point orig_wf_pt(wf);

      wf.three_point::operator=(wf_ap_intsctn_pt);
      ap.transform(wf);

      wf.three_point::operator=(orig_wf_pt);

    } catch(...) {
      cerr << "conic_mirror::private_transform - error aperturing wavefront\n";
      throw(string("conic_mirror::private_transform"));
    }

    // Assign the final wavefront three frame
    wf.three_frame::operator=(final_wf_tf);


  }
}
#endif
