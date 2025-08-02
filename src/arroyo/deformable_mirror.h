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

#ifndef DEFORMABLE_MIRROR_H
#define DEFORMABLE_MIRROR_H

#include <iostream>
#include <string>
#include <iomanip>
#include <complex>
#include "region_base.h"
#include "three_transformation.h"
#include "aperture.h"
#include "diffractive_wavefront.h"

namespace Arroyo 
{
  using std::string;
  using std::vector;
  using std::ostream;

  ///
  /// A virtual base class to represent a deformable_mirror
  ///

  class deformable_mirror_base : 
    virtual public plane_optic, 
    virtual public one_to_one_optic {
    
    public:

    ///////////////////////////////////////////
    ///  Null constructor
    deformable_mirror_base(){};

    ///////////////////////////////////////////
    ///  Virtual destructor
    virtual ~deformable_mirror_base(){};

    ///////////////////////////////////////////
    // Pure virtual functions, need to be 
    // defined in derived classes.
    ///////////////////////////////////////////

    ///////////////////////////////////////////
    /// Read from file
    virtual void read(const char* filename) = 0;

    ///////////////////////////////////////////
    /// Read from iofits object
    virtual void read(const iofits& iof) = 0;

    ///////////////////////////////////////////
    /// Write to file
    virtual void write(const char* filename) const = 0;

    ///////////////////////////////////////////
    /// Write to iofits object
    virtual void write(iofits& iof) const = 0;

    ///////////////////////////////////////////
    /// Define how to output to stdout
    virtual void print(ostream& os, const char* prefix="") const = 0;

    ///////////////////////////////////////////
    ///  Factory to construct deformable_mirror_bases from file
    static deformable_mirror_base* deformable_mirror_base_factory(const char* filename);

    ///////////////////////////////////////////
    ///  Factory to construct deformable_mirror_bases from an iofits instance
    static deformable_mirror_base* deformable_mirror_base_factory(const iofits& iof);
  };

  //////////////////////////////////////////////
  ///
  /// A template class to represent and ideal 
  /// deformable mirror.
  ///
  /// The dynamic model for the motion of the
  /// actuators does not contain any fancy
  /// motion, such as:  hysteresis, and uses
  /// a constant velocity approximation to
  /// move the actuators.
  ///
  /// As with the tip_tilt_mirror class, this class
  /// inherits a three_frame through the planar_optic
  /// and the orientation of this three_frame represents 
  /// that of the deformable mirror
  ///
  /// The member function
  /// 
  /// ideal_deformable_mirror::
  ///   update(pixel_array<double>& pixarr, double timestamp)
  ///
  /// may be used to issue a new command to the DM
  /// To effect an update, this member function 
  /// computes the current actuator positions at the time
  /// given by the timestamp passed to this function.  It then
  /// updates the actuator positions via a linear function
  /// of time.
  ///
  /// Attempting to issue an update command with a timestamp
  /// that predates the one store in this class is flagged
  /// as an error.
  ///
  /// To perform the transformation on a wavefront using the
  /// trasform member function, the class uses the actuator
  /// velocities to interpolate to the timestamp of the wavefront.
  /// The transformation is effected by applying the moving
  /// the actuator positions and substracting it from
  /// the wavefront phase.
  ///
  /// Supplying a wavefront with a timestamp that precedes
  /// the internal one is an error.
  ///

  template<class aperture_type> 
    class ideal_deformable_mirror : 
    public deformable_mirror_base,
    public aperture_type {

    private:

    static const bool factory_registration;
    string unique_name() const
      {
        return(string("ideal deformable mirror"));
      };

    protected:

    /// Actuator pitch, in meters
    double mirror_actuator_pitch;

    /// The mirror timestamp.
    /// At this time, the mirror surface
    /// is defined by the mirror actuator
    /// positions
    double mirror_timestamp;

    /// Actuator velocity, in meters per second
    double mirror_actuator_velocity;

    /// The current actuator commands, specified
    /// using a pixel_array.  The mirror surface 
    /// is heading towards this orientation with 
    /// constant linear velocity.  These commands
    /// are in meters.
    pixel_array<double> mirror_actuator_commands;

    /// Actuator positions at the time of the
    /// mirror timestamp
    pixel_array<double> mirror_actuator_positions;
      
    /// A template member function to perform
    /// transform on both float and double
    /// instantiations of wavefront.  This 
    /// is necessary because there is no 
    /// mechanism in C++ for virtual template
    /// member functions.
    template<class T> 
      void private_transform(diffractive_wavefront<T>& wf) const;
    
    ///////////////////////////////////////////
    ///  Null constructor
    ideal_deformable_mirror();

    public:

    ///////////////////////////////////////////
    ///  Copy constructor
    ideal_deformable_mirror(const ideal_deformable_mirror& ideal_dm);

    ///////////////////////////////////////////
    ///  Construct from pixel_array
    ideal_deformable_mirror(const aperture& ap, 
			    const vector<long> & nact,
			    const double pitch, 
			    const double vel, 
			    const double ts=0);
      
    ///////////////////////////////////////////
    ///  Construct from file
    ideal_deformable_mirror(const char *filename);

    ///////////////////////////////////////////
    ///  Construct from iofits
    ideal_deformable_mirror(const iofits &iof);

    ///////////////////////////////////////////
    ///  Virtual destructor
    virtual ~ideal_deformable_mirror(){};

    ///////////////////////////////////////////
    ///  Operator = 
    ideal_deformable_mirror& operator=(const ideal_deformable_mirror& ideal_dm);

    ///////////////////////////////////////////
    ///  Read from file
    virtual void read(const char* filename);

    ///////////////////////////////////////////
    ///  Read from an iofits object
    virtual void read(const iofits& iof);
 
    ///////////////////////////////////////////
    ///  Write to file
    void write(const char* filename) const;

    ///////////////////////////////////////////
    ///  Write to an iofits object
    void write(iofits& iof) const;

    ///////////////////////////////////////////
    ///  Print
    void print(ostream& os, const char* prefix="") const;

    ///////////////////////////////////////////
    /// Get the actuator pitch
    ///
    /// Returns the pitch in meters
    vector<long> get_axes(void) const 
      {
	return mirror_actuator_commands.get_axes();
      }

    ///////////////////////////////////////////
    /// Get the actuator pitch
    ///
    /// Returns the pitch in meters
    double get_actuator_pitch(void) const 
      {
	return mirror_actuator_pitch;
      }

    ///////////////////////////////////////////
    /// Set the actuator pitch
    ///
    /// act_pitch is in meters 
    void set_actuator_pitch(double act_pitch);

    ///////////////////////////////////////////
    /// Get the actuator velocity
    ///
    /// Returns the velocity in meters per second
    double get_actuator_velocity(void) const
      {
	return mirror_actuator_velocity;
      }

    ///////////////////////////////////////////
    /// Set the actuator velocity
    ///
    /// act_velocity is in meters per second
    void set_actuator_velocity(double act_velocity);

    ///////////////////////////////////////////
    ///  Get the last mirror actuator commands
    pixel_array<double> get_actuator_commands() const {
      return(mirror_actuator_commands);
    }

    ///////////////////////////////////////////
    ///  Get the mirror actuator positions at the
    ///  timestamp provided.
    ///
    ///  If the timestamp provided preceeds the mirror
    ///  timestamp, an error is thrown
    pixel_array<double> get_actuator_positions(double timestamp) const;

    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    void set_mirror_actuator_positions(const pixel_array<double>& pixarr) {
      this->mirror_actuator_positions=pixarr;
    }

    void set_mirror_actuator_commands(const pixel_array<double>& pixarr) {
      this->mirror_actuator_commands=pixarr;
    }
    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////
    /// Update the actuators by the output of the
    /// reconstructor at a given timestamp
    void update(const pixel_array<double>& pixarr, double timestamp);

    ///////////////////////////////////////////
    ///  Apply the aperture to the wavefront
    void transform(diffractive_wavefront<float> & wf) const;

    ///////////////////////////////////////////
    ///  Apply the aperture to the wavefront
    void transform(diffractive_wavefront<double> & wf) const;

    ///////////////////////////////////////////
    ///  Get the mirror timestamp
    ///
    ///  Timestamp is in seconds
    double get_timestamp() const {return(mirror_timestamp);};

  };

  template<class aperture_type> 
    ideal_deformable_mirror<aperture_type>::
    ideal_deformable_mirror()
    {
      mirror_actuator_pitch     = 0.;
      mirror_actuator_velocity  = 0.;
      mirror_timestamp          = 0.;
    }

  template<class aperture_type> 
    ideal_deformable_mirror<aperture_type>::
    ideal_deformable_mirror(const ideal_deformable_mirror<aperture_type> &ideal_dm)
    {
      this->operator=(ideal_dm);
    }

  template<class aperture_type> 
    ideal_deformable_mirror<aperture_type>::
    ideal_deformable_mirror(const char *filename)
    {
      this->read(filename);
    }

  template<class aperture_type> 
    ideal_deformable_mirror<aperture_type>::
    ideal_deformable_mirror(const iofits &iof)
    {
      this->read(iof);
    }

  template<class aperture_type> 
    ideal_deformable_mirror<aperture_type>::
    ideal_deformable_mirror(const aperture& ap, 
			    const vector<long> & actuator_dimensions,
			    const double actuator_pitch, 
			    const double actuator_velocity, 
			    const double actuator_timestamp)
    {
      const aperture_type* aptype;
      try{aptype = dynamic_cast<const aperture_type*>(&ap);}
      catch(...) {
	cerr << "ideal_deformable_mirror::ideal_deformable_mirror error -\n"
	     << "clash between instantiated type of "
	     << "mirror and aperture " << endl;
	cerr << "provided to constructor" << endl;
	ap.print(cerr,"ap provided to constructor");
	throw(this->unique_name() + string("::ideal_deformable_mirror constructor"));
      }
      this->aperture_type::operator=(*aptype);


      if(actuator_dimensions.size()!=2){
	cerr << "ideal_deformable_mirror::ideal_deformable_mirror error -\n"
	     << " cannot instantiate deformable mirror with " 
	     << actuator_dimensions.size() << " actuator dimensions\n";
	throw(string("ideal_deformable_mirror::ideal_deformable_mirror"));
      }
      if(actuator_dimensions[0]<=0 || actuator_dimensions[1]<=0){
	cerr << "ideal_deformable_mirror::ideal_deformable_mirror error -\n"
	     << " cannot instantiate deformable mirror with actuator dimensions "
	     << actuator_dimensions[0] << "x" << actuator_dimensions[1] << endl;
	throw(string("ideal_deformable_mirror::ideal_deformable_mirror"));
      }

      mirror_actuator_pitch     = actuator_pitch;
      mirror_actuator_velocity  = actuator_velocity;
      mirror_timestamp          = actuator_timestamp;

      mirror_actuator_positions = pixel_array<double>(actuator_dimensions);
      mirror_actuator_commands = pixel_array<double>(actuator_dimensions);
    }
  
  template<class aperture_type> 
    ideal_deformable_mirror<aperture_type> & 
    ideal_deformable_mirror<aperture_type>::
    operator=(const ideal_deformable_mirror<aperture_type>& ideal_dm)
    {
      if(this == &ideal_dm) return(*this);
      this->aperture_type::operator=(*this);
      mirror_actuator_pitch     = ideal_dm.mirror_actuator_pitch;
      mirror_actuator_velocity  = ideal_dm.mirror_actuator_velocity;
      mirror_timestamp          = ideal_dm.mirror_timestamp;
      mirror_actuator_positions = ideal_dm.mirror_actuator_positions;
      mirror_actuator_commands  = ideal_dm.mirror_actuator_commands;
      return(*this);
    }

  template<class aperture_type> 
    void ideal_deformable_mirror<aperture_type>::
    read(const char* filename)
    { 
      iofits iof;
      try
	{
	  iof.open(filename);
	}
      catch(...)
	{
	  cerr << this->unique_name() << "::read - error opening file " 
	       << filename << endl;
	  throw(this->unique_name() + string("::read"));
	}
    
      try
	{
	  this->read(iof);
	}
      catch(...)
	{
	  cerr << this->unique_name() << "::read - error reading " 
	       << this->unique_name() << " from file " << filename << endl;
	  throw(this->unique_name() + string("::read"));
	}
    }

  template<class aperture_type>
    void ideal_deformable_mirror<aperture_type>::
    read(const iofits& iof)
    { 
      string type, comment;
      if (!iof.key_exists("TYPE"))
	{
	  cerr << this->unique_name() << "::read error - unrecognized "
	       << "file type" << endl;
	  throw(this->unique_name() + string("::read"));
	}
      iof.read_key("TYPE",type,comment);
      if (type != "ideal_deformable_mirror")
	{
	  cerr << this->unique_name() << "::read error - file of type: " 
	       << type << " rather than type " << this->unique_name() << endl;
	  throw(this->unique_name() + string("::read"));
	}
    
      iof.read_key("ACTPITCH",mirror_actuator_pitch,comment);
      iof.read_key("ACT_VEL",mirror_actuator_velocity,comment);
      iof.read_key("TIMESTMP",mirror_timestamp,comment);
    
      mirror_actuator_positions.read(iof);

      // Here we advance to the next HDU, where the mirror actuator 
      // commands are stored.
      mirror_actuator_commands.read(iof);

      // Here we advance to the next HDU, where the aperture is stored.
      // This is stored here to avoid fits header key nameclashes with
      // the pixel_array that holds the last actuator positions
      this->aperture_type::read(iof);
    
      // leave iof pointing to the next header, if it exists
      if(iof.get_hdu_num()!=iof.get_num_hdus()) iof.movrel_hdu(1);
    }
  
  template<class aperture_type>
    void ideal_deformable_mirror<aperture_type>::
    write(const char* filename) const 
    {
      iofits iof;
      try
	{
	  iof.create(filename);
	}
      catch(...)
	{
	  cerr << this->unique_name() << "::write - error creating file " 
	       << filename << endl;
	  throw(this->unique_name() + string("::write"));
	}
    
      try
	{
	  this->write(iof);
	}
      catch(...)
	{
	  cerr << this->unique_name() << "::write - error writing " 
	       << this->unique_name() << " to file " << filename << endl;
	  throw(this->unique_name() + string("::write"));
	}
    }

  template<class aperture_type>
    void ideal_deformable_mirror<aperture_type>::
    write(iofits & iof) const 
    {
      string type, comment;
      vector<long> axes = mirror_actuator_positions.get_axes();
    
      fits_header_data<double> tmphdr(axes);
      tmphdr.write(iof);
    
      type    = "ideal_deformable_mirror";
      comment = "object type";
      iof.write_key("TYPE", type, comment);
      comment = "actuator pitch (meters)";
      iof.write_key("ACTPITCH",mirror_actuator_pitch,comment);
      comment = "actuator velocity (meters/secs)";
      iof.write_key("ACT_VEL",mirror_actuator_velocity,comment);
      comment = "timestamp (secs)";
      iof.write_key("TIMESTMP",mirror_timestamp,comment);
    
      mirror_actuator_positions.write(iof);
	
      // Write the mirror actuator commands to the next HDU
      tmphdr.write(iof);
      mirror_actuator_commands.write(iof);
    
      // Write the aperture to the next HDU
      this->aperture_type::write(iof);
 
    }
  
  template<class aperture_type>
    void ideal_deformable_mirror<aperture_type>::
    print(ostream & os, const char * prefix) const 
    {
      int vlspc = 30;
      vector<long> axes = mirror_actuator_positions.get_axes();

      os.setf(std::ios::left, std::ios::adjustfield); 
      os << prefix << "TYPE       = " << setw(vlspc) 
	 << this->unique_name() << "/" << "object type" << endl;
      fits_header_data<double> tmphdr(axes);
      tmphdr.print(os, prefix);
      os << prefix << "ACTPITCH   = " << setw(vlspc) << mirror_actuator_pitch
	 << "   actuator pitch (meters)" << endl;
      os << prefix << "ACT_VEL    = " << setw(vlspc) << mirror_actuator_velocity
	 << "   actuator velocity (meters/secs)" << endl;
      os << prefix << "TIMESTMP   = " << setw(vlspc)
	 << mirror_timestamp
	 << "   actuators time_stamp (secs)" << endl;
      this->aperture_type::print(os, prefix);
    }

  template<class aperture_type>
    void ideal_deformable_mirror<aperture_type>::
    set_actuator_pitch(double actuator_pitch)
    {
      if(actuator_pitch<0){
	cerr << "ideal_deformable_mirror::set_actuator_pitch error - "
	     << " cannot set actuator pitch to value " << actuator_pitch << endl;
	throw(string("ideal_deformable_mirror::set_actuator_pitch"));
      }
      mirror_actuator_pitch = actuator_pitch;
    }

  template<class aperture_type>
    void ideal_deformable_mirror<aperture_type>::
    set_actuator_velocity(double actuator_velocity)
    {
      if(actuator_velocity<0){
	cerr << "ideal_deformable_mirror::set_actuator_velocity error - "
	     << " cannot set actuator velocity to value " << actuator_velocity << endl;
	throw(string("ideal_deformable_mirror::set_actuator_velocity"));
      }
      mirror_actuator_velocity = actuator_velocity;
    }

  template<class aperture_type> 
    pixel_array<double> ideal_deformable_mirror<aperture_type>::
    get_actuator_positions(double timestamp) const
    {
      if (timestamp < this->mirror_timestamp)
	{
	  cerr << this->unique_name() << "::get_actuator_positions error - timestamp: " 
	       << timestamp << " predates the internal ideal_deformable_mirror timestamp "
	       << mirror_timestamp << endl;
	  throw(this->unique_name() + string("::get_actuator_positions"));
	}

      vector<long> axes = mirror_actuator_positions.get_axes();
      pixel_array<double> new_actuator_positions(axes);

      double displacement_magnitude = mirror_actuator_velocity*(timestamp-mirror_timestamp);
      double last_actuator_position, actuator_command_position, new_actuator_position;
      int nelem = mirror_actuator_positions.total_space();
      for (int i=0; i<nelem; i++)
	{
	  last_actuator_position = mirror_actuator_positions.data(i);
	  actuator_command_position = mirror_actuator_commands.data(i);
	  if(last_actuator_position<actuator_command_position)
	    new_actuator_position = 
	      last_actuator_position+displacement_magnitude<actuator_command_position ?
	      last_actuator_position+displacement_magnitude : actuator_command_position;
	  else 
	    new_actuator_position = 
	      last_actuator_position-displacement_magnitude>actuator_command_position ?
	      last_actuator_position-displacement_magnitude : actuator_command_position;
	  
	  new_actuator_positions.set_data(i,new_actuator_position);
	}
      return(new_actuator_positions);
    }

  template<class aperture_type> 
    void ideal_deformable_mirror<aperture_type>::
    update(const pixel_array<double>& act_cmds, double timestamp)
    {
      mirror_actuator_positions = this->get_actuator_positions(timestamp);
      mirror_actuator_commands = act_cmds;
      mirror_timestamp = timestamp;
    }

  template<class aperture_type>
    void ideal_deformable_mirror<aperture_type>::
    transform(diffractive_wavefront<float> & wf) const 
    {
      try
	{
	  this->private_transform(wf);
	}
      catch(...)
	{
	  cerr << this->unique_name() << "::transform error - "
	       << "error transforming a float instantiation " << endl;
	  cerr << "of diffractive_wavefront" << endl;
	  throw(this->unique_name() + string("::transform"));
	}
    }

  template<class aperture_type>
    void ideal_deformable_mirror<aperture_type>::
    transform(diffractive_wavefront<double> & wf) const 
    {
      try
	{
	  this->private_transform(wf);
	}
      catch(...)
	{
	  cerr << this->unique_name() << "::transform error - "
	       << "error transforming a double instantiation " << endl;
	  cerr << "of diffractive_wavefront" << endl;
	  throw(this->unique_name() + string("::transform"));
	}
    }


  template<class aperture_type>
    template<class T> void 
    ideal_deformable_mirror<aperture_type>::
    private_transform(diffractive_wavefront<T> & wf) const {
    three_vector origin_offset, dx, dy;

    // For the time being assume that the wavefront is normally
    // incident to the deformable_mirror.  The next code snipet tries to
    // insure this case for now.
    if (cross_product(this->z(),wf.z()).length() > 
	three_frame::precision && foreshortening) {
	cerr << this->unique_name() << "::private_transform error - "
	     << "The wavefront is not normally " << endl;
	cerr << "incident onto the deformable_mirror.  Aborting..." << endl;
	wf.z().print(cerr, "wf z ");
	this->z().print(cerr, "tt z ");
	throw(this->unique_name()+string("::private_transform"));
    }

    if (dot_product(wf.z(),this->z()) >= 0.) {
	cerr << this->unique_name() << "::private_transform error - "
	     << "The wavefront provided to this function is " << endl;
	cerr << "travelling in the wrong direction, i.e., is "
	     << "not incident on the reflective side of "
	     << "the mirror." << endl;
	wf.z().print(cerr, "wf z ");
	this->z().print(cerr, "tt z ");
	throw(this->unique_name()+string("::private_transform"));
      }

    if (wf.get_timestamp() < mirror_timestamp) {
	cerr << this->unique_name() << "::private_transform error - "
	     << "deformable_mirror last timestamp postdates wavefront " << endl;
	cerr << "timestamp: " << wf.get_timestamp() << " provided "
	     << "to this function." << endl;
	throw(string(this->unique_name()+"::private_transform"));
      }

    if (wf.get_pixel_scale() > mirror_actuator_pitch) {
	cerr << this->unique_name() << "::private_transform error - "
	     << "actuator pitch: " << mirror_actuator_pitch
	     << " is larger than " << endl;
	cerr << "the wavefront pixel scale: " << wf.get_pixel_scale()
	     << endl;
	throw(this->unique_name()+string("::private_transform"));
      }

    try{get_projected_wavefront_pixel_spacing(wf, origin_offset, dx, dy);} 
    catch(...){
      cerr << this->unique_name()
	   << "::get_projected_wavefront_pixel_spacing error - "
	   << endl;
      cerr << "could not get a projected pixel scale" << endl;
      throw(this->unique_name()+string("::private_transform"));
    }

    // Do the transformation now:

    // a) First apply the aperture:
    //
    // Do the aperture transformation which it would zero out the
    // wavefront if it is incident outside the optic.

    this->aperture_type::transform(wf);

    // b) apply the actuators.
    //
    // Move the actuators from the store positions to the 
    // position corresponding to the wf.timestamp.

    pixel_array<double> current_actuator_positions(this->get_actuator_positions(wf.get_timestamp()));

    // The halfpixel information
    vector<long> wf_axes  = wf.get_axes();
    vector<long> dm_axes  = mirror_actuator_positions.get_axes();

    int wf_X_min     = -wf_axes[1]/2;
    int wf_X_max     =  wf_axes[1]/2 + wf_axes[1]%2;
    int wf_Y_min     = -wf_axes[0]/2;
    int wf_Y_max     =  wf_axes[0]/2 + wf_axes[0]%2;
    double wf_X_hpix = ((wf_axes[1]%2==0) ? 0.5 : 0.);
    double wf_Y_hpix = ((wf_axes[0]%2==0) ? 0.5 : 0.);

    wavefront_amp_phase_conversion(wf);
    T *wfdata = get_wavefront_data(wf);

    three_point center_pixel; 
    double act_X_hpix = ((dm_axes[1]%2==0) ? 0.5 : 0.);
    int act_X_extrapix = ((dm_axes[1]%2==0) ? 0 : 1);
    double act_Y_hpix = ((dm_axes[0]%2==0) ? 0.5 : 0.);
    int act_Y_extrapix = ((dm_axes[0]%2==0) ? 0 : 1);

    int iact, jact;
    int stride = is_interleaved_storage(wf) ? 2 : 1;
    int wf_amp_index = -stride;
    int wf_phase_index  = is_interleaved_storage(wf) ? 1-stride : wf_axes[0]*wf_axes[1]-stride;
    double a, b, x, y;
    vector<int> actuator_index(4);
    vector<double> actuator_weight(4);
    double kvector = 2*M_PI/wf.get_wavelength();

    for(int i=wf_X_min; i<wf_X_max; i++){
      for(int j=wf_Y_min; j<wf_Y_max; j++){
	  
	wf_amp_index += stride;
	wf_phase_index += stride;
	  
	// Don't bother with the phase if the amplitude is zero
	if(wfdata[wf_amp_index]==0) continue;
	  
	center_pixel  = *this+(origin_offset + +(i+wf_X_hpix)*dx + (j+wf_Y_hpix)*dy);
	// Coordinates of wf pixel center in dm frame, measured
	// in units of the mirror actuator pitch
	x = center_pixel.x(*this)/mirror_actuator_pitch;
	y = center_pixel.y(*this)/mirror_actuator_pitch;
	  
	// Condition for the wf pixel to lie beyond the
	// rectangular boundary defined by the limits of the
	// pyramidal influence functions of the dm edge
	// actuators.  In this case we may skip this pixel.
	if(fabs(x)>dm_axes[1]/2+act_X_extrapix+act_X_hpix ||
	   fabs(y)>dm_axes[0]/2+act_Y_extrapix+act_Y_hpix){
	  continue;
	}
	  
	// Condition for the wf pixel to lie within the
	// rectangular boundary defined by the centers of the dm
	// edge actuators.  This is the case in which 4
	// actuators contribute to the wf phase.
	if(fabs(x)<dm_axes[1]/2-act_X_hpix &&
	   fabs(y)<dm_axes[0]/2-act_Y_hpix){
	    
	  iact = (int)(x+dm_axes[1]/2-act_X_hpix);
	  jact = (int)(y+dm_axes[0]/2-act_Y_hpix);
	    
	  a = fabs(x+dm_axes[1]/2-act_X_hpix-iact);
	  b = fabs(y+dm_axes[0]/2-act_Y_hpix-jact);
	    
	  actuator_weight[0] = (1-a)*(1-b);
	  actuator_weight[1] = a*(1-b);
	  actuator_weight[2] = (1-a)*b;
	  actuator_weight[3] = a*b;
	  actuator_index[0] = jact*dm_axes[0]+iact;
	  actuator_index[1] = actuator_index[0]+1;
	  actuator_index[2] = (jact+1)*dm_axes[0]+iact;
	  actuator_index[3] = actuator_index[2]+1;

	  try{
	    for(int k=0; k<4; k++){
	      wfdata[wf_phase_index] -= 
		2*kvector*current_actuator_positions.data(actuator_index[k])*actuator_weight[k];
	    }
	  } catch(...) {
	    cerr << "ideal_deformable_mirror::private_transform error - first error\n";
	    for(int k=0; k<4; k++)
	      cerr << "actuator " << k
		   << "\tindex " << actuator_index[k]
		   << "\tactuator position axes " << current_actuator_positions.get_axes()[0]
		   << "x" << current_actuator_positions.get_axes()[1]
		   << endl;
	    throw(string("ideal_deformable_mirror::private_transform"));
	  }
	} else {
	  // In the remaining cases, the wf pixel lies outside
	  // the boundary defined by the centers of the dm edge
	  // actuators, but inside the boundary defined by the
	  // pyramidal influence functions of the dm edge
	  // actuators.  In these cases, one or two actuators
	  // contribute to the wf phase.
	    
	  // In these 4 cases, only the corner actuator contributes to the wf phase
	  try{
	    if(x<-dm_axes[1]/2+act_X_hpix && y<-dm_axes[0]/2+act_Y_hpix){
	      wfdata[wf_phase_index] -= 
		2*kvector*current_actuator_positions.data(0)*
		(1+dm_axes[1]/2-act_X_hpix+x)*(1+dm_axes[0]/2-act_Y_hpix+y);
	    } else if(x<-dm_axes[1]/2+act_X_hpix && y>dm_axes[0]/2-act_Y_hpix){
	      wfdata[wf_phase_index] -= 
		2*kvector*current_actuator_positions.data(dm_axes[1]*(dm_axes[0]-1))*
		(1+dm_axes[1]/2-act_X_hpix+x)*(1+dm_axes[0]/2-act_Y_hpix-y);
	    } else if(x>dm_axes[1]/2-act_X_hpix && y<-dm_axes[0]/2+act_Y_hpix){
	      wfdata[wf_phase_index] -= 
		2*kvector*current_actuator_positions.data(dm_axes[0]-1)*
		(1+dm_axes[1]/2-act_X_hpix-x)*(1+dm_axes[0]/2-act_Y_hpix+y);
	    } else if(x>dm_axes[1]/2-act_X_hpix && y>dm_axes[0]/2-act_Y_hpix){
	      wfdata[wf_phase_index] -= 
		2*kvector*current_actuator_positions.data(dm_axes[0]*dm_axes[1]-1)*
		(1+dm_axes[1]/2-act_X_hpix-x)*(1+dm_axes[0]/2-act_Y_hpix-y);
	    }
	  
	    // In these 4 cases, two edge actuators contribute to the wf phase
	    else {
	      if(x<-dm_axes[0]/2+act_X_hpix){
		jact = (int)(y+dm_axes[0]/2-act_Y_hpix);
		a = -dm_axes[1]/2+act_X_hpix-x;
		b = fabs(y+dm_axes[0]/2-act_Y_hpix-jact);
		actuator_weight[0] = (1-a)*(1-b);
		actuator_weight[1] = (1-a)*b;
		actuator_index[0] = jact*dm_axes[0];
		actuator_index[1] = (jact+1)*dm_axes[0];
	      } else if(x>dm_axes[1]/2-act_X_hpix){
		jact = (int)(y+dm_axes[0]/2-act_Y_hpix);
		a = -dm_axes[1]/2+act_X_hpix+x;
		b = fabs(y+dm_axes[0]/2-act_Y_hpix-jact);
		actuator_weight[0] = (1-a)*(1-b);
		actuator_weight[1] = (1-a)*b;
		actuator_index[0] = (jact+1)*dm_axes[0]-1;
		actuator_index[1] = (jact+2)*dm_axes[0]-1;
	      } else if(y<-dm_axes[0]/2+act_Y_hpix){
		iact = (int)(x+dm_axes[1]/2-act_X_hpix);
		a = fabs(x+dm_axes[1]/2-act_X_hpix-iact);
		b = -dm_axes[0]/2+act_Y_hpix-y;
		actuator_weight[0] = (1-a)*(1-b);
		actuator_weight[1] = a*(1-b);
		actuator_index[0] = iact;
		actuator_index[1] = iact+1;
	      } else if(y>dm_axes[0]/2-act_Y_hpix){
		iact = (int)(x+dm_axes[1]/2-act_X_hpix);
		a = fabs(x+dm_axes[1]/2-act_X_hpix-iact);
		b = -dm_axes[0]/2+act_Y_hpix+y;
		actuator_weight[0] = (1-a)*(1-b);
		actuator_weight[1] = a*(1-b);
		actuator_index[0] = (dm_axes[1]-1)*dm_axes[0]+iact;
		actuator_index[1] = (dm_axes[1]-1)*dm_axes[0]+iact+1;
	      }

	      for(int k=0; k<2; k++)
		wfdata[wf_phase_index] -= 
		  2*kvector*current_actuator_positions.data(actuator_index[k])*actuator_weight[k];		  
	    }
	  } catch(...) {
	    cerr << "ideal_deformable_mirror::private_transform error - third error\n";
	    throw(string("ideal_deformable_mirror::private_transform"));
	  }  
	}
      }
    }
    // Reflect the wavefront
    three_reflection tref(wf, this->z());
    tref.transform(wf);
  }
}
#endif
