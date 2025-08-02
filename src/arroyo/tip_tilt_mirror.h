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

#ifndef TIP_TILT_MIRROR_H
#define TIP_TILT_MIRROR_H

#include <iomanip>
#include "region_base.h"
#include "three_transformation.h"
#include "optic.h"
#include "aperture.h"
#include "diffractive_wavefront.h"

namespace Arroyo {

  using std::vector;
  using std::ostream;

  template<typename U> class diffractive_wavefront;

  ///
  /// A virtual base class to represent tip tilt mirrors.
  ///

  class tip_tilt_mirror_base :
    virtual public plane_optic,
    virtual public one_to_one_optic {

    public:
    
    ///////////////////////////////////////////
    ///  Null constructor
    tip_tilt_mirror_base(){};

    ///////////////////////////////////////////
    ///  Virtual destructor
    virtual ~tip_tilt_mirror_base(){};

     ///////////////////////////////////////////
    ///  Read from file
    virtual void read(const char * filename) = 0;
  
    ///////////////////////////////////////////
    ///  Read from iofits
    virtual void read(const iofits & iof) = 0;
 
    ///////////////////////////////////////////
    ///  Write to file
    virtual void write(const char * filename) const = 0;
 
    ///////////////////////////////////////////
    ///  Write to iofits
    virtual void write(iofits & iof) const = 0;
 
    ///////////////////////////////////////////
    ///  Print
    virtual void print(ostream & os, const char * prefix="") const = 0;

    ///////////////////////////////////////////
    ///  Factory to construct tip_tilt_mirror_bases from file
    static tip_tilt_mirror_base * tip_tilt_mirror_base_factory(const char * filename);

    ///////////////////////////////////////////
    ///  Factory to construct tip_tilt_mirror_bases from an iofits instance
    static tip_tilt_mirror_base * tip_tilt_mirror_base_factory(const iofits & iof);

  };

  /////////////////////////////////////////////////////////////
  ///
  ///  A template class to represent an idealized tip tilt
  ///  mirror.  One may instantiate this class using an 
  ///  aperture type. 
  ///
  ///  The dynamic model does not include effects
  ///  of hysteresis, and uses a constant angular 
  ///  velocity approximation to the mirror orientation.
  ///  (i.e. the mirror moves at constant angular velocity 
  ///  towards its commanded position until it reaches 
  ///  that position, and then stops).
  ///
  ///  This class inherits three_frame through planar_optic,
  ///  and the orientation of this three frame represents
  ///  the mirror orientation at a particular time.  This
  ///  time is held as a data member within the class.  The
  ///  class also holds a three vector that represents the
  ///  latest mirror orientation command.  
  ///
  ///  The member function 
  ///
  ///  ideal_tip_tilt_mirror::update(three_vector & tv, double & timestamp)
  ///
  ///  may be used to issue a new command to the tip tilt mirror.
  ///  To effect an update, this member function computes the current
  ///  orientation of the mirror at the time given by the timestamp
  ///  passed to this function.  It then moves the aperture to this
  ///  orientation, sets the internal timestamp to that passed into
  ///  the update function, and sets the latest mirror command to be
  ///  equal to the three_vector passed into the update function. 
  ///
  ///  Attempting to issue an update command with a timestamp that
  ///  predates that stored within the class is an error.
  ///
  ///
  ///  To perform the transformation on a wavefront using the
  ///  transform member function, the class uses the mirror
  ///  angular velocity to interpolate to the timestamp of 
  ///  the wavefront in the transform member function.  
  ///  The transformation is performed by applying the aperture
  ///  transformation, and then adjusting the wavefront three_frame
  ///  to account for the reflection.
  ///
  ///  Supplying a wavefront with a timestamp that preceeds 
  ///  the internal one is an error.
  ///
  ///

  template<class aperture_type>
  class ideal_tip_tilt_mirror :
    public tip_tilt_mirror_base,
    public aperture_type {

    private:

    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("ideal tip tilt mirror"));};

    protected:

    /// The mirror orientation timestamp.
    /// At this time, the mirror orientation
    /// is defined by the aperture three frame.
    double mirror_timestamp;

    /// The current orientation command, specified
    /// using a vector normal to the surface of
    /// the mirror.  The mirror is heading towards
    /// this orientation with constant angular
    /// velocity.
    three_vector mirror_orientation_command;

    /// The angular velocity, in rad/sec
    double mirror_angular_velocity;

    /// A template member function to perform
    /// transform on both float and double
    /// instantiations of wavefront.  This 
    /// is necessary because there is no 
    /// mechanism in C++ for virtual template
    /// member functions.
    template<class T>
      void private_transform(diffractive_wavefront<T> & wf) const;

    ///////////////////////////////////////////
    ///  Null constructor.
    ///
    ///  Protected because instances must be initialized
    ideal_tip_tilt_mirror(){};

    public:

    ///////////////////////////////////////////
    ///  Copy constructor
    ideal_tip_tilt_mirror(const ideal_tip_tilt_mirror & ideal_ttm);

    ///////////////////////////////////////////
    ///  Construct from file
    ideal_tip_tilt_mirror(const char * filename);

    ///////////////////////////////////////////
    ///  Construct from iofits object
    ideal_tip_tilt_mirror(const iofits & iof);

    ///////////////////////////////////////////
    ///  Construct from the bits
    ideal_tip_tilt_mirror(const aperture_type & ap, 
			  double angular_velocity,
			  double timestamp = 0);

    ///////////////////////////////////////////
    ///  Destructor
    ~ideal_tip_tilt_mirror(){};

    ///////////////////////////////////////////
    ///  Operator=
    ideal_tip_tilt_mirror & operator=(const ideal_tip_tilt_mirror & ideal_ttm);

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
    ///  Apply the aperture to the wavefront
    virtual void transform(diffractive_wavefront<float> & wf) const;

    ///////////////////////////////////////////
    ///  Apply the aperture to the wavefront
    virtual void transform(diffractive_wavefront<double> & wf) const;

    ///////////////////////////////////////////
    ///  Update the mirror orientation.  Starting
    ///  at the timestamp provided, the mirror will
    ///  start to move towards a new orientation at
    ///  a constant angular velocity.  The new orientation
    ///  is defined by orientation_command, which
    ///  is a three_vector orthogonal to the plane of the
    ///  mirror. 
    ///  
    void update(const three_vector & orientation_command, double timestamp);

    ///////////////////////////////////////////
    ///  Get the last mirror orientation command
    three_vector get_orientation_command() const {
      return(mirror_orientation_command);
    }

    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    void set_mirror_command_vector(three_vector ttm_command_vector)	{
      this->mirror_orientation_command = ttm_command_vector;
    }
    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////
    ///  Get the mirror orientation at the
    ///  timestamp provided.
    ///
    ///  If the timestamp provided preceeds the mirror
    ///  timestamp, an error is thrown
    three_frame get_orientation(double timestamp) const;

    ///////////////////////////////////////////
    ///  Get the mirror angular velocity
    ///
    ///  Angular velocity is in rad/sec
    double get_angular_velocity() const {return(mirror_angular_velocity);};

    ///////////////////////////////////////////
    ///  Set the mirror angular velocity
    ///
    ///  Angular velocity is in rad/sec
    void set_angular_velocity(double angular_velocity);

    ///////////////////////////////////////////
    ///  Get the mirror timestamp
    ///
    ///  Timestamp is in seconds
    double get_timestamp() const {return(mirror_timestamp);};

    ///////////////////////////////////////////
    ///  Get the point of intersection 
    /// 
    ///  There's a problem here.  Normally one would
    ///  call optic::get_point_of_intersection() to 
    ///  determine the point of intersection of a ray
    ///  with the optic, and then propagate a wavefront
    ///  to this point.  By default, this function is 
    ///  covered by plane_optic::get_point_of_intersection
    ///  However, this will return the point of intersection 
    ///  at the mirror timestamp.  What we probably want is
    ///  the point of intersection at the time the wavefront
    ///  arrives at the optic.  However, we don't know 
    ///  what time it will arrive because we don't know the 
    ///  propagation distance...so we need a new function
    ///  that understands the speed of light and can solve for 
    ///  the point of intersection at a particular timestamp.
    ///  Something like
    //three_point get_point_of_intersection(const three_point & tp, 
    //const three_vector & tv,
    //double timestamp) const;

  };

  template<class aperture_type>
  ideal_tip_tilt_mirror<aperture_type>::
    ideal_tip_tilt_mirror(const ideal_tip_tilt_mirror<aperture_type> & ttm){
    this->operator=(ttm);
  }

  template<class aperture_type>
  ideal_tip_tilt_mirror<aperture_type>::
    ideal_tip_tilt_mirror(const char * filename){
    this->read(filename);
  }

  template<class aperture_type>
  ideal_tip_tilt_mirror<aperture_type>::
    ideal_tip_tilt_mirror(const iofits & iof){
    this->read(iof);
  }

  template<class aperture_type>
  ideal_tip_tilt_mirror<aperture_type>::
    ideal_tip_tilt_mirror(const aperture_type & ap, 
			  double angular_velocity,
			  double timestamp){

    /*
    const aperture_type * aptype;
    try{aptype = dynamic_cast<const aperture_type *>(&ap);}
    catch(...){
      cerr << "ideal_tip_tilt_mirror::ideal_tip_tilt_mirror error - "
	   << "clash between instantiated type of mirror and aperture provided to constructor\n";
      ap.print(cerr, "ap provided to constructor ");
      throw(string("ideal_tip_tilt_mirror::ideal_tip_tilt_mirror"));
    }
    this->aperture_type::operator=(*aptype);
    */

    this->aperture_type::operator=(ap);
    mirror_timestamp = timestamp;
    mirror_angular_velocity = angular_velocity;
    mirror_orientation_command = ap.z();
  }
  
  template<class aperture_type>
  ideal_tip_tilt_mirror<aperture_type> & 
    ideal_tip_tilt_mirror<aperture_type>::
    operator=(const ideal_tip_tilt_mirror<aperture_type> & ttm){
    if(this==&ttm) return(*this);
    mirror_timestamp = ttm.mirror_timestamp;
    mirror_orientation_command = ttm.mirror_orientation_command;
    mirror_angular_velocity = ttm.mirror_angular_velocity;
    this->aperture_type::operator=(ttm);
    return(*this);
  }

  template<class aperture_type>
  void ideal_tip_tilt_mirror<aperture_type>::
    read(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "ideal_tip_tilt_mirror::read - "
	   << "error opening file " << filename << endl;
      throw(string("ideal_tip_tilt_mirror::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "ideal_tip_tilt_mirror::read - "
	   << "error reading "
	   << this->unique_name() << " from file "
	   << filename << endl;
      throw(string("ideal_tip_tilt_mirror::read"));
    }
  }

  template<class aperture_type>
  void ideal_tip_tilt_mirror<aperture_type>::
    read(const iofits & iof) {
    if(!iof.key_exists("TYPE")){
      cerr << "ideal_tip_tilt_mirror::read error - "
	   << "unrecognized type of file\n";
      throw(string("ideal_tip_tilt_mirror::read"));
    }
    string type, comment;
    iof.read_key("TYPE", type, comment);
    if(type!="ideal tip tilt mirror"){
      cerr << "ideal_tip_tilt_mirror::read error - file of type " 
	   << type << " rather than of type circular aperture\n";
      throw(string("ideal_tip_tilt_mirror::read"));
    }

    iof.read_key("TIMESTMP", mirror_timestamp, comment);
    iof.read_key("ANGVEL", mirror_angular_velocity, comment);
    mirror_orientation_command.read(iof);

    // Here we advance to the next HDU, where the aperture is stored.
    // This is stored here to avoid fits header key nameclashes with
    // the three frame that holds the last orientation
    iof.movrel_hdu(1);
    this->aperture_type::read(iof);

    // leave iof pointing to the next header, if it exists
    if(iof.get_hdu_num()!=iof.get_num_hdus())
      iof.movrel_hdu(1);
  }

  template<class aperture_type>
  void ideal_tip_tilt_mirror<aperture_type>::
    write(const char * filename) const {

    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "ideal_tip_tilt_mirror::write - "
	   << "error opening file " << filename << endl;
      throw(string("ideal_tip_tilt_mirror::write"));
    }

    try{this->write(iof);}
    catch(...){
      cerr << "ideal_tip_tilt_mirror::write - "
	   << "error writing "
	   << this->unique_name() << " to file "
	   << filename << endl;
      throw(string("ideal_tip_tilt_mirror::write"));
    }
  }

  template<class aperture_type>
  void ideal_tip_tilt_mirror<aperture_type>::
    write(iofits & iof) const {

    Arroyo::fits_header_data<char> fhd;
    fhd.write(iof);
    string type = "ideal tip tilt mirror";
    string comment = "object type";
    iof.write_key("TYPE", type, comment);
    iof.write_key("TIMESTMP", mirror_timestamp, "timestamp (secs)");
    iof.write_key("ANGVEL", mirror_angular_velocity, "mirror angular velocity (rad/sec)");
    mirror_orientation_command.write(iof);

    // Write the aperture to the next HDU
    this->aperture_type::write(iof);
  }

  template<class aperture_type>
  void ideal_tip_tilt_mirror<aperture_type>::
    print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(std::ios::left, std::ios::adjustfield); 
    os << prefix << "TYPE       = " << std::setw(vlspc) << this->unique_name()
       << "/" << "object type" << endl;
    os << prefix << "TIMESTMP   = " << std::setw(vlspc) << mirror_timestamp
       << "/" << "timestamp for last orientation (secs)" << endl;
    os << prefix << "ANGVEL     = " << std::setw(vlspc) << mirror_angular_velocity
       << "/" << "mirror angular velocity (rad/sec)" << endl;
    mirror_orientation_command.print(os, prefix);
    this->aperture_type::print(os, prefix);
  }

  template<class aperture_type>
    three_frame ideal_tip_tilt_mirror<aperture_type>::
    get_orientation(double timestamp) const {

    if(timestamp<this->mirror_timestamp){
      cerr << "ideal_tip_tilt_mirror::get_orientation error - "
	   << "last mirror orientation timestamp " 
	   << mirror_timestamp 
	   << " postdates timestamp "
	   << timestamp 
	   << " provided to this function\n";
      throw(string("ideal_tip_tilt_mirror::get_orientation"));
    }

    three_vector axis_of_rotation = 
      cross_product(this->aperture_type::z(), mirror_orientation_command);

    // Special case - mirror is already at the mirror orientation command.
    if(axis_of_rotation.length() < three_frame::precision){
      return(three_frame(*this));
    }

    double rotation_angle = asin(axis_of_rotation.length());
    if(rotation_angle > this->mirror_angular_velocity*(timestamp-mirror_timestamp))
      rotation_angle = this->mirror_angular_velocity*(timestamp-mirror_timestamp);

    // Make a copy of the aperture and rotate it so that it is
    // oriented along the direction it reached at the timestamp
    // in the wavefront
    three_rotation trot(*this, axis_of_rotation, rotation_angle);

    aperture_type aptype(*this);
    trot.transform(aptype);
    return(aptype);
  }

  template<class aperture_type>
    void ideal_tip_tilt_mirror<aperture_type>::
    set_angular_velocity(double angular_velocity){
    if(angular_velocity<=0){
      cerr << "ideal_tip_tilt_mirror::set_angular_velocity error - "
	   << "cannot set mirror angular velocity to " 
	   << angular_velocity << endl;
      throw(string("ideal_tip_tilt_mirror::set_angular_velocity"));
    }
    this->mirror_angular_velocity = angular_velocity;
  } 

  template<class aperture_type>
  void ideal_tip_tilt_mirror<aperture_type>::
    update(const three_vector & orientation_command, 
	   double timestamp){

    if(orientation_command.length()<three_frame::precision){
      cerr << "ideal_tip_tilt_mirror::update error - "
	   << "orientation vector provided to this function is null\n";
      throw(string("ideal_tip_tilt_mirror::update"));
    } 

    if(timestamp<this->mirror_timestamp){
      cerr << "ideal_tip_tilt_mirror::update error - "
	   << "last mirror orientation timestamp " 
	   << mirror_timestamp 
	   << " postdates timestamp "
	   << timestamp 
	   << " provided to this function\n";
      throw(string("ideal_tip_tilt_mirror::update"));
    }

    three_vector normalized_orientation_command = 
      orientation_command*(1/orientation_command.length());

    this->three_frame::operator=(this->get_orientation(timestamp));
    mirror_timestamp = timestamp;
    mirror_orientation_command = normalized_orientation_command;
  }

  template<class aperture_type>
  void ideal_tip_tilt_mirror<aperture_type>::
    transform(diffractive_wavefront<float> & wf) const {
    try{this->private_transform(wf);}
    catch(...){
      cerr << "ideal_tip_tilt_mirror::transform error - "
	   << "error transforming a float instantiation of diffractive_wavefront\n";
      throw(string("ideal_tip_tilt_mirror::transform"));
    }
  }

  template<class aperture_type>
  void ideal_tip_tilt_mirror<aperture_type>::
    transform(diffractive_wavefront<double> & wf) const {
    try{this->private_transform(wf);}
    catch(...){
      cerr << "ideal_tip_tilt_mirror::transform error - "
	   << "error transforming a double instantiation of diffractive_wavefront\n";
      throw(string("ideal_tip_tilt_mirror::transform"));
    }
  }

  template<class aperture_type>
  template<class T>
  void ideal_tip_tilt_mirror<aperture_type>::
    private_transform(diffractive_wavefront<T> & wf) const { 

    // Check to ensure that the wavefront incident on the mirror
    // arrives at a time which postdates the last mirror update.
    if(wf.get_timestamp()<this->mirror_timestamp){
      cerr << "ideal_tip_tilt_mirror::private_transform error - "
	   << "mirror orientation timestamp " 
	   << mirror_timestamp 
	   << " postdates timestamp "
	   << wf.get_timestamp()
	   << " provided to this function\n";
      throw(string("ideal_tip_tilt_mirror::private_transform"));
    }

    // Check that the wavefront is incident on the reflective side
    // of the mirror, which is defined by the z basis vector of the
    // aperture
    if(dot_product(wf.z(), this->z())>=0){
      cerr << "ideal_tip_tilt_mirror::private_transform error - " << endl
	   << "wavefront supplied to this function is not incident "
	   << "on the reflective side of the mirror,\n";

      wf.z().print(cerr, "wf z ");
      this->z().print(cerr, "tt z ");
      throw(string("ideal_tip_tilt_mirror::private_transform"));
    }

    // Make a copy of the aperture and set its three_frame equal to
    // the mirror orientation at the time of the wavefront timestamp
    aperture_type aptype(*this);
    aptype.three_frame::operator=(this->get_orientation(wf.get_timestamp()));

    // Use the aperture copy to transform the wavefront 
    try{aptype.transform(wf);}
    catch(...){
      cerr << "ideal_tip_tilt_mirror::private_transform error - "
	   << "failed to transform wavefront using aperture\n";
	wf.print(cerr, "wf ");
	aptype.print(cerr, "aperture ");
	throw(string("ideal_tip_tilt_mirror::private_transform"));
    }

    // Reflect the wavefront
    three_reflection tref(wf, aptype.z());
    tref.transform(wf);
  }
}

#endif
