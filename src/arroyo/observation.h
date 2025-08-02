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

#ifndef OBSERVATION_H
#define OBSERVATION_H

#include <complex>
#include <fstream>
#include "colormap.h"
#include "fft_manager.h"
#include "pixel_amp_array.h"
#include "pixel_phase_array.h"
#include "diffractive_wavefront.h"

namespace Arroyo {


  template <typename precision> class basic_otf;

  ///
  /// A base class for the observation class hierarchy
  ///

  class observation_base {

    public:

    ///////////////////////////////////////////
    ///  Null constructor
    observation_base(){};

    ///////////////////////////////////////////
    ///  Virtual destructor
    virtual ~observation_base(){};

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
    virtual void write(iofits & iof) const = 0;

    ///////////////////////////////////////////
    /// Define how to output to stdout
    virtual void print(ostream & os, const char* prefix="") const = 0;

    ///////////////////////////////////////////
    ///  Factory to construct observation_bases from file
    static observation_base* observation_base_factory(const char* filename);

    ///////////////////////////////////////////
    ///  Factory to construct observation_bases from an iofits instance
    static observation_base* observation_base_factory(const iofits& iof);
  };


  template<class precision>
    class basic_observation :
    public AO_sim_base,
    public observation_base, 
    public Arroyo::pixel_amp_array<precision> {

    private:

    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("basic observation"));};

    protected:

    double right_ascension_radians;
    
    double declination_radians;
    
    double pixel_scale_radians_per_pixel;
    
    double wavelength_meters;

    public:

    ///////////////////////////////////////////
    ///  Null constructor
    basic_observation();

    ///////////////////////////////////////////
    ///  Destructor
    ~basic_observation(){};

    ///////////////////////////////////////////
    ///  Copy constructor
    basic_observation(const basic_observation & basic_obs){
      this->operator=(basic_obs);
    };

    ///////////////////////////////////////////
    ///  Construct from file
    basic_observation(const char * filename){
      this->read(filename);
    };

    ///////////////////////////////////////////
    ///  Construct from an iofits object
    basic_observation(const iofits & iof){
      this->read(iof);
    };

    ///////////////////////////////////////////
    ///  Construct from the bits
    basic_observation(double ra_radians,
		      double dec_radians,
		      double pixel_scale_radians_per_pixel,
		      double wavelength_meters,
		      const pixel_amp_array<precision> & pixamparr);


    ///////////////////////////////////////////
    ///  Construct from an OTF
    basic_observation(double aperture_diameter_meters, 
		      double focal_plane_image_size_arcsecs,
		      double oversampling_factor,
		      const basic_otf<precision> & otf);

    /*
    ///////////////////////////////////////////
    ///  Construct from an OTF
    basic_observation(int psf_dimen,
		      double wavelength_meters,
		      double pixel_scale_radians_per_pixel,
		      double aperture_diameter_meters,
		      pixel_amp_array<precision> OTF_amps,
		      pixel_phase_array<precision> OTF_phases);
    */

    ///////////////////////////////////////////
    ///  Operator =
    basic_observation & operator=(const basic_observation & basic_obs);

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
    ///  Function to print the coefficients
    void print(std::ostream & os, const char * prefix = "") const;

    ///////////////////////////////////////////
    ///  Get right ascension in radians
    double get_ra() const {return(right_ascension_radians);};
   
    ///////////////////////////////////////////
    ///  Get declination in radians
    double get_dec() const {return(declination_radians);};
   
    ///////////////////////////////////////////
    ///  Get pixel scale in radians per pixel
    double get_pixel_scale() const {return(pixel_scale_radians_per_pixel);};
   
    ///////////////////////////////////////////
    ///  Get wavelength in meters
    double get_wavelength() const {return(wavelength_meters);};
   
    ///////////////////////////////////////////
    ///  Function to return ensquared energy as a function
    ///  of aperture width.
    ///
    ///  First element of the pair is the slit width
    ///  in arcseconds
    ///  Second element is the ensquared energy
    vector<pair<double, double> > ensquared_energy() const;


    /*
    ///////////////////////////////////////////
    ///  Function to compute the optical transfer function
    void optical_transfer_function(int otf_dimen,
				   double aperture_diameter_meters,
				   pixel_amp_array<precision> & OTF_amps,
				   pixel_phase_array<precision> & OTF_phases) const;
    */
  };

  template <class precision> 
    basic_otf<precision> operator * (const basic_otf<precision> & p1,
				     double & fac);
  
  template <class precision> 
    basic_otf<precision> operator * (const basic_otf<precision> & p1,
				     const basic_otf<precision> & p2);
  
  template <class precision> 
    basic_otf<precision> operator / (const basic_otf<precision> & p1,
				     const basic_otf<precision> & p2);
  
  template <class precision> 
    basic_otf<precision> operator + (const basic_otf<precision> & p1,
				     const basic_otf<precision> & p2);
  

  template<class precision>
    class basic_otf :
    public AO_sim_base {

    private:

    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("optical transfer function"));};

    protected:

    double right_ascension_radians;
    
    double declination_radians;
    
    double pupil_plane_pixel_scale_meters;
    
    double wavelength_meters;

    pixel_amp_array<precision> OTF_amps;
    pixel_phase_array<precision> OTF_phases;

    public:

    ///////////////////////////////////////////
    ///  Null constructor
    basic_otf();

    ///////////////////////////////////////////
    ///  Destructor
    ~basic_otf(){};

    ///////////////////////////////////////////
    ///  Copy constructor
    basic_otf(const basic_otf & otf){
      this->operator=(otf);
    };

    ///////////////////////////////////////////
    ///  Construct from file
    basic_otf(const char * filename){
      this->read(filename);
    };

    ///////////////////////////////////////////
    ///  Construct from an iofits object
    basic_otf(const iofits & iof){
      this->read(iof);
    };

    ///////////////////////////////////////////
    ///  Construct from the bits
    basic_otf(double ra_radians,
	      double dec_radians,
	      double pupil_plane_pixel_scale_meters,
	      double wavelength_meters,
	      const pixel_amp_array<precision> & otf_amps,
	      const pixel_phase_array<precision> & otf_phases);
    
    ///////////////////////////////////////////
    ///  Construct from a basic observation
    basic_otf(double aperture_diameter_meters, 
	      double pupil_plane_pixel_scale_meters, 
	      const basic_observation<precision> & basic_obs);

    ///////////////////////////////////////////
    ///  Operator =
    basic_otf & operator=(const basic_otf & otf);

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
    ///  Function to print the otf
    void print(std::ostream & os, const char * prefix = "") const;

    ///////////////////////////////////////////
    ///  Get right ascension in radians
    double get_ra() const {return(right_ascension_radians);};
   
    ///////////////////////////////////////////
    ///  Get declination in radians
    double get_dec() const {return(declination_radians);};
   
    ///////////////////////////////////////////
    ///  Get pixel scale in radians per pixel
    double get_pixel_scale() const {return(pupil_plane_pixel_scale_meters);};
   
    ///////////////////////////////////////////
    ///  Get wavelength in meters
    double get_wavelength() const {return(wavelength_meters);};
   
    ///////////////////////////////////////////
    ///  Get axes
    vector<long> get_axes() const {return(this->OTF_amps.get_axes());};
   
    ///////////////////////////////////////////
    ///  Get data
    complex<precision> data(int index) const {
      return(polar(this->OTF_amps.data(index), 
		   this->OTF_phases.data(index)));
    };

    ///////////////////////////////////////////
    ///  Get OTF amps
    pixel_amp_array<precision> get_amps() const {
      return(OTF_amps);
    };

    ///////////////////////////////////////////
    ///  Get OTF phases
    pixel_phase_array<precision> get_phases() const {
      return(OTF_phases);
    };

    ///////////////////////////////////////////
    ///  Set data
    void set_data(int index, complex<precision> cmplx) const {
      this->OTF_amps.set_data(index, abs(cmplx)); 
      this->OTF_phases.set_data(index, arg(cmplx)); 
    };

    ///////////////////////////////////////////
    ///  pad the arrays
    void pad(int padlength) {
      this->OTF_amps.pad_array(padlength);
      this->OTF_phases.pad_array(padlength);
    }

    ///////////////////////////////////////////
    ///  clip the arrays
    void clip(int cliplength) {
      this->OTF_amps.clip_array(cliplength);
      this->OTF_phases.clip_array(cliplength);
    }

    ///////////////////////////////////////////
    ///  Operator *=  for otfs
    basic_otf<precision> & operator*=(const basic_otf<precision> & rhs);

    ///////////////////////////////////////////
    ///  Operator /=  for otfs
    basic_otf<precision> & operator/=(const basic_otf<precision> & rhs);

    ///////////////////////////////////////////
    ///  Operator *=  for otfs and doubles
    basic_otf<precision> & operator*=(double fac);

    ///////////////////////////////////////////
    ///  Operator +=  for otfs
    basic_otf<precision> & operator+=(const basic_otf<precision> & rhs);

    ///////////////////////////////////////////
    ///  Operator -=  for otfs
    basic_otf<precision> & operator-=(const basic_otf<precision> & rhs);

  };

  template<class precision>
    basic_observation<precision>::
    basic_observation(){
    
    this->right_ascension_radians = 0;
    this->declination_radians = 0;
    this->pixel_scale_radians_per_pixel = -1;
    this->wavelength_meters = -1;

  }

  template<class precision>
    basic_observation<precision>::
    basic_observation(double ra_radians,
		      double dec_radians,
		      double pixscale_radians_per_pixel,
		      double wavelength_meters,
		      const pixel_amp_array<precision> & pixamparr){

    this->right_ascension_radians = ra_radians;
    this->declination_radians = dec_radians;
    this->pixel_scale_radians_per_pixel = pixscale_radians_per_pixel;
    this->wavelength_meters = wavelength_meters;

    this->pixel_amp_array<precision>::operator=(pixamparr);

  }

  /*
  template<typename precision>
    basic_observation<precision>::basic_observation(int psf_dimen,
						    double wavelength_meters,
						    double pixel_scale_radians_per_pixel,
						    double aperture_diameter_meters,
						    pixel_amp_array<precision> OTF_amps,
						    pixel_phase_array<precision> OTF_phases) {


    try{
      
      if(OTF_amps.get_axes()!=OTF_phases.get_axes()){
	cerr << "basic_observation::basic_observation error - axes mismatch\n";
	throw(string("basic_observation::basic_observation"));
      }

      if(aperture_diameter_meters<=0){
	cerr << "basic_observation::basic_observation error - aperture diameter "
	     << aperture_diameter_meters
	     << " out of range\n";
	throw(string("basic_observation::basic_observation"));
      }

      this->wavelength_meters = wavelength_meters;
      this->pixel_scale_radians_per_pixel = pixel_scale_radians_per_pixel;
      this->set_axes(vector<long>(2,psf_dimen));
      
      // set up the storage
      vector<long> otf_axes = OTF_amps.get_axes();
      vector<long> psf_axes(2, psf_dimen);
      
      int notf_elems = OTF_amps.total_space();
      int npsf_elems = this->axes[0]*this->axes[1];
      precision *otf_data;
      precision *psf_data;    
      try{
	otf_data = new precision[2*notf_elems];
	psf_data = new precision[2*npsf_elems];
      } catch(...) {
	cerr << "basic_observation::basic_observation - error allocating memory\n";
	throw(string("basic_observation::basic_observation"));
      }
      
      // Slope correction for even dimension arrays
      double xslope = this->axes[1]%2==1 ? 0 : M_PI/(double)(psf_axes[1]/2);
      double yslope = this->axes[0]%2==1 ? 0 : M_PI/(double)(psf_axes[0]/2);

      // Halfpixel information
      double otf_x_halfpix=0, otf_y_halfpix=0;
      int otf_x_extrapix=1, otf_y_extrapix=1;
      if(otf_axes[1]%2==0){
	otf_x_halfpix = .5;
	otf_x_extrapix = 0;
      }
      if(otf_axes[0]%2==0){
	otf_y_halfpix = .5;
	otf_y_extrapix = 0;
      }
      double psf_x_halfpix=0, psf_y_halfpix=0;
      int psf_x_extrapix=1, psf_y_extrapix=1;
      if(this->axes[1]%2==0){
	psf_x_halfpix = .5;
	psf_x_extrapix = 0;
      }
      if(this->axes[0]%2==0){
	psf_y_halfpix = .5;
	psf_y_extrapix = 0;
      }

      // Initialize the complex PSF data array
      int index;
      double twopi = 2*M_PI;
      double amp, phase;
      for(int i=-otf_axes[1]/2; i<otf_axes[1]/2+otf_x_extrapix; i++){
	for(int j=-otf_axes[0]/2; j<otf_axes[0]/2+otf_y_extrapix; j++){
	  index = (i+otf_axes[1]/2)*otf_axes[0]+j+otf_axes[0]/2;
	  amp = OTF_amps.data(index);
	  phase = OTF_phases.data(index) + fmod(-xslope*(i+otf_x_halfpix) - yslope*(j+otf_y_halfpix), twopi);
	  otf_data[2*index] = amp*cos(phase);
	  otf_data[2*index+1] = amp*sin(phase);
	}
      }

      double oversampling_factor = 
	wavelength_meters / pixel_scale_radians_per_pixel / aperture_diameter_meters;

      //cerr << "oversampling factor " << oversampling_factor << endl;

      double extra_fac = 2*(double)this->axes[0]/(double)otf_axes[0];

      //double sampling_factor = -2/(double)otf_axes[0]/oversampling_factor;
      double sampling_factor = -extra_fac/(double)this->axes[0]/oversampling_factor;

      //double normalization_factor = fabs(sampling_factor);

      goertzel_reinsch_transform(otf_axes, this->axes, sampling_factor, otf_data, psf_data);

      delete [] otf_data;

      // No need to fix up the slopes here - we just want the amplitude
      for(int i=0; i<npsf_elems; i++)
	this->pixeldata[i] = 
	  sqrt(psf_data[2*i]*psf_data[2*i] + psf_data[2*i+1]*psf_data[2*i+1]);

      delete [] psf_data;

      } catch(...) {
      cerr << "basic_observation::basic_observation error\n";
      throw(string("basic_observation::basic_observation"));
    }
  }
  */


  template<typename precision>
    basic_observation<precision>::basic_observation(double aperture_diameter_meters,
						    double focal_plane_image_size_arcsecs,
						    double oversampling_factor,
						    const basic_otf<precision> & otf) {


    try{
      
      if(aperture_diameter_meters<=0){
	cerr << "basic_observation::basic_observation error - aperture diameter "
	     << aperture_diameter_meters
	     << " out of range\n";
	throw(string("basic_observation::basic_observation"));
      }

      if(focal_plane_image_size_arcsecs<=0){
	cerr << "basic_observation::basic_observation error - focal plane image size "
	     << focal_plane_image_size_arcsecs
	     << " out of range\n";
	throw(string("basic_observation::basic_observation"));
      }

      if(oversampling_factor<=0){
	cerr << "basic_observation::basic_observation error - oversampling factor "
	     << oversampling_factor
	     << " out of range\n";
	throw(string("basic_observation::basic_observation"));
      }

      double arcsecs_to_radians = M_PI/180./3600.;

      vector<long> otf_axes = otf.get_axes();
      double nyquist_pixel_scale_radians_per_pixel = 
	otf.get_wavelength() / 2. / aperture_diameter_meters;
      this->wavelength_meters = otf.get_wavelength();
      this->pixel_scale_radians_per_pixel = nyquist_pixel_scale_radians_per_pixel/oversampling_factor;
      long psf_dimen = (long)(ceil(focal_plane_image_size_arcsecs*arcsecs_to_radians/this->pixel_scale_radians_per_pixel));
      vector<long> psf_axes(2,psf_dimen);
      this->set_axes(psf_axes);
      
      // set up the storage      
      int notf_elems = otf_axes[0]*otf_axes[1];
      int npsf_elems = psf_axes[0]*psf_axes[1];
      precision *otf_data;
      precision *psf_data;    
      try{
	otf_data = new precision[2*notf_elems];
	psf_data = new precision[2*npsf_elems];
      } catch(...) {
	cerr << "basic_observation::basic_observation - error allocating memory\n";
	throw(string("basic_observation::basic_observation"));
      }
      
      // Slope correction for even dimension arrays
      double fac = 2*this->pixel_scale_radians_per_pixel/nyquist_pixel_scale_radians_per_pixel;
      double xslope = psf_axes[1]%2==1 ? 0 : M_PI*fac/(double)(otf_axes[1]);
      double yslope = psf_axes[0]%2==1 ? 0 : M_PI*fac/(double)(otf_axes[0]);

      // Halfpixel information
      double otf_x_halfpix=0, otf_y_halfpix=0;
      int otf_x_extrapix=1, otf_y_extrapix=1;
      if(otf_axes[1]%2==0){
	otf_x_halfpix = .5;
	otf_x_extrapix = 0;
      }
      if(otf_axes[0]%2==0){
	otf_y_halfpix = .5;
	otf_y_extrapix = 0;
      }
      double psf_x_halfpix=0, psf_y_halfpix=0;
      int psf_x_extrapix=1, psf_y_extrapix=1;
      if(this->axes[1]%2==0){
	psf_x_halfpix = .5;
	psf_x_extrapix = 0;
      }
      if(this->axes[0]%2==0){
	psf_y_halfpix = .5;
	psf_y_extrapix = 0;
      }

      // Initialize the complex PSF data array
      int index;
      double twopi = 2*M_PI;
      double amp, phase;
      complex<precision> cmp;
      for(int i=-otf_axes[1]/2; i<otf_axes[1]/2+otf_x_extrapix; i++){
	for(int j=-otf_axes[0]/2; j<otf_axes[0]/2+otf_y_extrapix; j++){
	  index = (i+otf_axes[1]/2)*otf_axes[0]+j+otf_axes[0]/2;
	  cmp = otf.data(index);
	  amp = abs(cmp);
	  phase = fmod(arg(cmp) - xslope*(i+otf_x_halfpix) - yslope*(j+otf_y_halfpix), twopi);
	  otf_data[2*index] = amp*cos(phase);
	  otf_data[2*index+1] = amp*sin(phase);
	}
      }

      //double sampling_factor = -2/(double)otf_axes[0]/oversampling_factor;
      double sampling_factor = -1/(double)otf_axes[0]/oversampling_factor;

      goertzel_reinsch_transform(otf_axes, psf_axes, sampling_factor, otf_data, psf_data);

      delete [] otf_data;

      // No need to fix up the slopes here - we just want the amplitude
      for(int i=0; i<npsf_elems; i++)
	this->pixeldata[i] = 
	  sqrt(psf_data[2*i]*psf_data[2*i] + psf_data[2*i+1]*psf_data[2*i+1]);

      delete [] psf_data;

      } catch(...) {
      cerr << "basic_observation::basic_observation error\n";
      throw(string("basic_observation::basic_observation"));
    }
  }


  template<class precision>
  basic_observation<precision> & basic_observation<precision>::
  operator=(const basic_observation<precision> & basic_obs){
    if(this==&basic_obs)
      return(*this);

    this->right_ascension_radians = basic_obs.right_ascension_radians;
    this->declination_radians = basic_obs.declination_radians;
    this->pixel_scale_radians_per_pixel = basic_obs.pixel_scale_radians_per_pixel;
    this->wavelength_meters = basic_obs.wavelength_meters;

    this->pixel_amp_array<precision>::operator=(basic_obs);

    return(*this);
  }

  template<class precision>
  void basic_observation<precision>::read(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "basic_observation::read - "
	   << "error opening file " << filename << endl;
      throw(string("basic_observation::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "basic_observation::read - "
	   << "error reading "
	   << this->unique_name() << " from file "
	   << filename << endl;
      throw(string("basic_observation::read"));
    }
  }

  template<class precision>
  void basic_observation<precision>::read(const iofits & iof){

    if(!iof.key_exists("TYPE")){
      cerr << "basic_observation::read error - "
	   << "unrecognized type of file\n";
      iof.print_header(cerr, "hdr dump");
      cerr << "hdu num " << iof.get_hdu_num() << " of "
	   << iof.get_num_hdus() << endl;
      throw(string("basic_observation::read"));
    }

    string type, comment;
    iof.read_key("TYPE", type, comment);
    if(type!=this->unique_name()){
      cerr << "basic_observation::read error - file of type " 
	   << type << " rather than of type "
	   << this->unique_name() << endl;
      throw(string("basic_observation::read"));
    }

    iof.read_key("CDELT1", this->pixel_scale_radians_per_pixel, comment);
    iof.read_key("CRVAL1", this->right_ascension_radians, comment);
    iof.read_key("CRVAL2", this->declination_radians, comment);

    double radians_to_degrees = 180./M_PI;
    double degrees_to_radians = M_PI/180.;


    this->pixel_scale_radians_per_pixel *= -1*degrees_to_radians;
    this->right_ascension_radians *= degrees_to_radians;
    this->declination_radians *= degrees_to_radians;

    comment = "observation wavelength (meters)";
    iof.read_key("WVLNGTH", this->wavelength_meters, comment);

    this->pixel_amp_array<precision>::read(iof);

  }

  template<class precision>
  void basic_observation<precision>::write(const char * filename) const {
    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "basic_observation::write - "
	   << "error opening file " << filename << endl;
      throw(string("basic_observation::write"));
    }
    try{this->write(iof);}
    catch(...){
      cerr << "basic_observation::write - "
	   << "error writing " 
	   << this->unique_name() << " to file "
	   << filename << endl;
      throw(string("basic_observation::write"));
    }
  }

  template<class precision>
  void basic_observation<precision>::write(iofits & iof) const {

    double radians_to_degrees = 180./M_PI;

    fits_header_data<double> fhd(this->get_axes());
    fhd.write(iof);
    string type = this->unique_name();
    string comment = "object type";
    iof.write_key("TYPE", type, comment);

    iof.write_key("RADECSYS", string("FK4"), string("SYSTEM OF REF. COORD"));
    iof.write_key("EQUINOX", (long)2000, string("EQUINOX OF REF. COORD"));
    iof.write_key("CTYPE1", string("RA---TAN"), string("AXIS TYPE"));
    iof.write_key("CTYPE2", string("DEC--TAN"), string("AXIS TYPE"));
    iof.write_key("CD001001", (long)1, string("UNITLESS"));
    iof.write_key("CD002001", (long)0, string("UNITLESS"));
    iof.write_key("CD001002", (long)0, string("UNITLESS"));
    iof.write_key("CD002002", (long)1, string("UNITLESS"));
    iof.write_key("CDELT1", -1*this->pixel_scale_radians_per_pixel*radians_to_degrees, string("DEGREES/PIXEL"));
    iof.write_key("CDELT2", this->pixel_scale_radians_per_pixel*radians_to_degrees, string("DEGREES/PIXEL"));
    iof.write_key("CRPIX1", (long)(this->get_axes()[1]/2), string("REFERENCE PIXEL IN X"));
    iof.write_key("CRPIX2", (long)(this->get_axes()[0]/2), string("REFERENCE PIXEL IN Y"));
    iof.write_key("CRVAL1", this->right_ascension_radians*radians_to_degrees, string("R.A. (DEGREES)"));
    iof.write_key("CRVAL2", this->declination_radians*radians_to_degrees, string("DEC. (DEGREES)"));

    comment = "observation wavelength (meters)";
    iof.write_key("WVLNGTH", this->wavelength_meters, comment);

    this->pixel_amp_array<precision>::write(iof);
  }

  template<class precision>
  void basic_observation<precision>::
  print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << this->unique_name()
       << "/" << "object type" << endl;
    fits_header_data<double> fhd(this->get_axes());
    fhd.print(os, prefix);
    os << prefix << "RA         = " << setw(vlspc) << this->right_ascension_radians
       << "/" << "right ascension (radians)" << endl;

    os << prefix << "DEC        = " << setw(vlspc) << this->declination_radians
       << "/" << "declination (radians)" << endl;

    os << prefix << "PIXSCALE   = " << setw(vlspc) << this->pixel_scale_radians_per_pixel
       << "/" << "pixel scale (radians)" << endl;

    os << prefix << "WVLNGTH    = " << setw(vlspc) << this->wavelength_meters
       << "/" << "wavelength (meters)" << endl;

  }


  template<class precision>
    vector<pair<double, double> > 
    basic_observation<precision>::ensquared_energy() const {

    vector<long> axes = this->get_axes();
    int extrapix = 1;
    double halfpix = 0;
    if(axes[0]%2==0){
      extrapix = 0;
      halfpix = .5;
    }

    int min_slit_width_pixels = extrapix==0 ? 2 : 1;
    int max_slit_width_pixels = axes[0];

    int upper_left_index, upper_right_index, lower_left_index, lower_right_index;
    double ensquared_energy = 0;

    vector<pair<double, double> > slit_widths_and_EEs;

    for(int k=min_slit_width_pixels; k<max_slit_width_pixels; k+=2){
      
      if(k==1) 
	ensquared_energy += this->pixel_amp_array<double>::pixeldata[(axes[0]+1)*(axes[0]/2)];
      else {
	upper_left_index = ((axes[0]-k)/2)*(axes[0]+1);
	upper_right_index = upper_left_index + (k-1);
	lower_left_index = upper_left_index + axes[0]*(k-1);
	lower_right_index = lower_left_index + (k-1);
	    
	for(int m=upper_left_index; m<=upper_right_index; m++){
	  ensquared_energy += this->pixel_amp_array<double>::pixeldata[m];
	}
	
	for(int m=1; m<=(k-2); m++){
	  ensquared_energy += this->pixel_amp_array<double>::pixeldata[upper_left_index+axes[0]*m];
	  ensquared_energy += this->pixel_amp_array<double>::pixeldata[upper_right_index+axes[0]*m];
	}
	
	for(int m=lower_left_index; m<=lower_right_index; m++){
	  ensquared_energy += this->pixel_amp_array<double>::pixeldata[m];
	}
      }
      slit_widths_and_EEs.push_back(pair<double,double>(k*this->pixel_scale_radians_per_pixel, 
							ensquared_energy));
      
    }
    return(slit_widths_and_EEs);
  }

  /*
  template<typename precision>
    void basic_observation<precision>::
    (int otf_dimen,
	      double aperture_diameter_meters,
	      pixel_amp_array<precision> & OTF_amps,
	      pixel_phase_array<precision> & OTF_phases) const {

    try{
     if(otf_dimen<=0){
	cerr << "basic_observation::basic_otf error - OTF dimension "
	     << otf_dimen
	     << " out of range\n";
	throw(string("basic_observation::basic_otf"));
      }
     if(aperture_diameter_meters<=0){
	cerr << "basic_observation::basic_otf error - aperture diameter "
	     << aperture_diameter_meters
	     << " out of range\n";
	throw(string("basic_observation::basic_otf"));
      }

     // set up the storage
     vector<long> otf_axes(2, otf_dimen);
     
     int notf_elems = otf_axes[0]*otf_axes[1];
     int npsf_elems = this->axes[0]*this->axes[1];
     precision *otf_data;
     precision *psf_data;    
     try{
       otf_data = new precision[2*notf_elems];
       psf_data = new precision[2*npsf_elems];
     } catch(...) {
       cerr << "wavefront_phase_estimate::point_spread_function - error allocating memory\n";
	throw(string("wavefront_phase_estimate::point_spread_function"));
      }

      // Slope correction for even dimension arrays
      double xslope = this->axes[1]%2==1 ? 0 : M_PI/(double)(otf_axes[1]/2);
      double yslope = this->axes[0]%2==1 ? 0 : M_PI/(double)(otf_axes[0]/2);

      // Halfpixel information
      double psf_x_halfpix=0, psf_y_halfpix=0;
      int psf_x_extrapix=1, psf_y_extrapix=1;
      if(this->axes[1]%2==0){
	psf_x_halfpix = .5;
	psf_x_extrapix = 0;
      }
      if(this->axes[0]%2==0){
	psf_y_halfpix = .5;
	psf_y_extrapix = 0;
      }
      double otf_x_halfpix=0, otf_y_halfpix=0;
      int otf_x_extrapix=1, otf_y_extrapix=1;
      if(otf_axes[1]%2==0){
	otf_x_halfpix = .5;
	otf_x_extrapix = 0;
      }
      if(otf_axes[0]%2==0){
	otf_y_halfpix = .5;
	otf_y_extrapix = 0;
      }

      // Initialize the complex PSF data array
      int index;
      double twopi = 2*M_PI;
      double amp, phase;
      for(int i=-this->axes[1]/2; i<this->axes[1]/2+psf_x_extrapix; i++){
	for(int j=-this->axes[0]/2; j<this->axes[0]/2+psf_y_extrapix; j++){
	  index = (i+this->axes[1]/2)*this->axes[0]+j+this->axes[0]/2;
	  amp = this->pixeldata[index];
	  phase = fmod(-xslope*(i+psf_x_halfpix) - yslope*(j+psf_y_halfpix), twopi);
	  psf_data[2*index] = amp*cos(phase);
	  psf_data[2*index+1] = amp*sin(phase);
	}
      }

      double oversampling_factor = 
	wavelength_meters / pixel_scale_radians_per_pixel / aperture_diameter_meters;

      //cerr << "oversampling factor " << oversampling_factor << endl;

      //double sampling_factor = 2/(double)this->axes[0]/oversampling_factor;
      double sampling_factor = 2/(double)otf_axes[0]/oversampling_factor;

      goertzel_reinsch_transform(this->axes, otf_axes, sampling_factor, psf_data, otf_data);

      delete [] psf_data;

      OTF_amps = pixel_amp_array<precision>(otf_axes);
      OTF_phases = pixel_phase_array<precision>(otf_axes);

      for(int i=0; i<notf_elems; i++){
	amp = sqrt(otf_data[2*i]*otf_data[2*i]+otf_data[2*i+1]*otf_data[2*i+1]);
	phase = atan2(otf_data[2*i+1],otf_data[2*i]);
	OTF_amps.set_data(i, amp);
	OTF_phases.set_data(i, phase);
      }
    } catch(...) {
      cerr << "basic_observation::basic_otf error\n";
      throw(string("basic_observation::basic_otf"));
    }
  }
  */















  template<class precision>
    basic_otf<precision>::
    basic_otf(){

    this->OTF_amps = pixel_amp_array<precision>(vector<long>(2,0));
    this->OTF_phases = pixel_phase_array<precision>(vector<long>(2,0));
    
    this->right_ascension_radians = 0;
    this->declination_radians = 0;
    this->pupil_plane_pixel_scale_meters = -1;
    this->wavelength_meters = -1;
  }

  template<class precision>
    basic_otf<precision>::
    basic_otf(double ra_radians,
	      double dec_radians,
	      double pupil_plane_pixscale_meters,
	      double wavelength_meters,
	      const pixel_amp_array<precision> & otf_amps,
	      const pixel_phase_array<precision> & otf_phases){
    
    this->right_ascension_radians = ra_radians;
    this->declination_radians = dec_radians;
    this->pupil_plane_pixel_scale_meters = pupil_plane_pixscale_meters;
    this->wavelength_meters = wavelength_meters;

    this->OTF_amps = otf_amps;
    this->OTF_phases = otf_phases;

  }


  template<typename precision>
    basic_otf<precision>::basic_otf(double aperture_diameter_meters,
				    double pupil_plane_pixscale_meters,
				    const basic_observation<precision> & basic_obs) {
    

    try{
      
      if(aperture_diameter_meters<=0){
	cerr << "basic_otf::basic_otf error - pupil plane pixel scale "
	     << aperture_diameter_meters
	     << " out of range\n";
	throw(string("basic_otf::basic_otf"));
      }
      if(pupil_plane_pixscale_meters<=0){
	cerr << "basic_otf::basic_otf error - pupil plane pixel scale "
	     << pupil_plane_pixscale_meters
	     << " out of range\n";
	throw(string("basic_otf::basic_otf"));
      }

      this->wavelength_meters = basic_obs.get_wavelength();
      this->pupil_plane_pixel_scale_meters = pupil_plane_pixscale_meters;

      // set up the storage
      int otf_dimen = (int)ceil(2*aperture_diameter_meters/pupil_plane_pixscale_meters);
      if(otf_dimen%2==0) otf_dimen++;
      vector<long> otf_axes(2,otf_dimen);
      vector<long> psf_axes = basic_obs.get_axes();
      
      int notf_elems = otf_dimen*otf_dimen;
      int npsf_elems = psf_axes[0]*psf_axes[1];
      precision *otf_data;
      precision *psf_data;    
      try{
	otf_data = new precision[2*notf_elems];
	psf_data = new precision[2*npsf_elems];
      } catch(...) {
	cerr << "basic_otf::basic_otf - error allocating memory\n";
	throw(string("basic_otf::basic_otf"));
      }
      
      // Slope correction for even dimension arrays
      double xslope = otf_axes[1]%2==1 ? 0 : M_PI/(double)(psf_axes[1]);
      double yslope = otf_axes[0]%2==1 ? 0 : M_PI/(double)(psf_axes[0]);

      // Halfpixel information
      double otf_x_halfpix=0, otf_y_halfpix=0;
      int otf_x_extrapix=1, otf_y_extrapix=1;
      if(otf_axes[1]%2==0){
	otf_x_halfpix = .5;
	otf_x_extrapix = 0;
      }
      if(otf_axes[0]%2==0){
	otf_y_halfpix = .5;
	otf_y_extrapix = 0;
      }
      double psf_x_halfpix=0, psf_y_halfpix=0;
      int psf_x_extrapix=1, psf_y_extrapix=1;
      if(psf_axes[1]%2==0){
	psf_x_halfpix = .5;
	psf_x_extrapix = 0;
      }
      if(psf_axes[0]%2==0){
	psf_y_halfpix = .5;
	psf_y_extrapix = 0;
      }

      // Initialize the complex PSF data array
      int index;
      double twopi = 2*M_PI;
      double amp, phase;
      for(int i=-psf_axes[1]/2; i<psf_axes[1]/2+psf_x_extrapix; i++){
	for(int j=-psf_axes[0]/2; j<psf_axes[0]/2+psf_y_extrapix; j++){
	  index = (i+psf_axes[1]/2)*psf_axes[0]+j+psf_axes[0]/2;
	  amp = basic_obs.data(index);
	  phase = (i+psf_x_halfpix)*xslope + (j+psf_y_halfpix)*yslope;
	  psf_data[2*index] = amp*cos(phase);
	  psf_data[2*index+1] = amp*sin(phase);
	}
      }

      /*
      double oversampling_factor = 
	wavelength_meters / pixel_scale_radians_per_pixel / aperture_diameter_meters;
      double sampling_factor = 2/(double)otf_axes[0]/oversampling_factor;
      */

      double oversampling_factor = 
	wavelength_meters / basic_obs.get_pixel_scale() / aperture_diameter_meters;
      double sampling_factor = 2/(double)otf_axes[0]/oversampling_factor;
      //double sampling_factor = 1/(double)psf_axes[0]/oversampling_factor;

      goertzel_reinsch_transform(psf_axes, otf_axes, sampling_factor, psf_data, otf_data);

      delete [] psf_data;

      this->OTF_amps = pixel_amp_array<precision>(vector<long>(2,otf_dimen));
      this->OTF_phases = pixel_phase_array<precision>(vector<long>(2,otf_dimen));

      double nyquist_pixel_scale_radians_per_pixel = 
	basic_obs.get_wavelength() / (double) (otf_axes[0]/2) / this->pupil_plane_pixel_scale_meters;
      double fac = basic_obs.get_pixel_scale()/nyquist_pixel_scale_radians_per_pixel;

      xslope = psf_axes[1]%2==1 ? 0 : -2*M_PI*fac/(double)otf_axes[1];
      yslope = psf_axes[0]%2==1 ? 0 : -2*M_PI*fac/(double)otf_axes[0];

      for(int i=-otf_axes[1]/2; i<otf_axes[1]/2+otf_x_extrapix; i++){
	for(int j=-otf_axes[0]/2; j<otf_axes[0]/2+otf_y_extrapix; j++){
	  index = (i+otf_axes[1]/2)*otf_axes[0]+j+otf_axes[0]/2;
	  amp = sqrt(otf_data[2*index]*otf_data[2*index] + otf_data[2*index+1]*otf_data[2*index+1]);
	  phase = fmod(atan2(otf_data[2*index+1],otf_data[2*index]) -
		       xslope*(i+otf_x_halfpix) - 
		       yslope*(j+otf_y_halfpix), twopi);
	  this->OTF_amps.set_data(index, amp);
	  this->OTF_phases.set_data(index, phase);
	}
      }

      delete [] otf_data;

      } catch(...) {
      cerr << "basic_otf::basic_otf error\n";
      throw(string("basic_otf::basic_otf"));
    }
  }


  template<class precision>
  basic_otf<precision> & basic_otf<precision>::
  operator=(const basic_otf<precision> & otf){
    if(this==&otf)
      return(*this);

    this->right_ascension_radians = otf.right_ascension_radians;
    this->declination_radians = otf.declination_radians;
    this->pupil_plane_pixel_scale_meters = otf.pupil_plane_pixel_scale_meters;
    this->wavelength_meters = otf.wavelength_meters;

    this->OTF_amps = otf.OTF_amps;
    this->OTF_phases = otf.OTF_phases;

    return(*this);
  }

  template<class precision>
  void basic_otf<precision>::read(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "basic_otf::read - "
	   << "error opening file " << filename << endl;
      throw(string("basic_otf::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "basic_otf::read - "
	   << "error reading "
	   << this->unique_name() << " from file "
	   << filename << endl;
      throw(string("basic_otf::read"));
    }
  }

  template<class precision>
  void basic_otf<precision>::read(const iofits & iof){

    if(!iof.key_exists("TYPE")){
      cerr << "basic_otf::read error - "
	   << "unrecognized type of file\n";
      iof.print_header(cerr, "hdr dump");
      cerr << "hdu num " << iof.get_hdu_num() << " of "
	   << iof.get_num_hdus() << endl;
      throw(string("basic_otf::read"));
    }

    string type, comment;
    iof.read_key("TYPE", type, comment);
    if(type!=this->unique_name()){
      cerr << "basic_otf::read error - file of type " 
	   << type << " rather than of type "
	   << this->unique_name() << endl;
      throw(string("basic_otf::read"));
    }

    iof.read_key("PIXSCALE", this->pupil_plane_pixel_scale_meters, comment);
    iof.read_key("CRVAL1", this->right_ascension_radians, comment);
    iof.read_key("CRVAL2", this->declination_radians, comment);

    double radians_to_degrees = 180./M_PI;
    double degrees_to_radians = M_PI/180.;

    comment = "observation wavelength (meters)";
    iof.read_key("WVLNGTH", this->wavelength_meters, comment);

    this->OTF_amps.read(iof);

    this->OTF_phases.read(iof);

  }

  template<class precision>
  void basic_otf<precision>::write(const char * filename) const {
    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "basic_otf::write - "
	   << "error opening file " << filename << endl;
      throw(string("basic_otf::write"));
    }
    try{this->write(iof);}
    catch(...){
      cerr << "basic_otf::write - "
	   << "error writing " 
	   << this->unique_name() << " to file "
	   << filename << endl;
      throw(string("basic_otf::write"));
    }
  }

  template<class precision>
  void basic_otf<precision>::write(iofits & iof) const {

    double radians_to_degrees = 180./M_PI;

    fits_header_data<double> fhd(this->OTF_amps.get_axes());
    fhd.write(iof);
    string type = this->unique_name();
    string comment = "object type";
    iof.write_key("TYPE", type, comment);

    comment = "pupil plane pixel scale (meters)";
    iof.write_key("PIXSCALE", this->pupil_plane_pixel_scale_meters, comment);
    iof.write_key("CRVAL1", this->right_ascension_radians*radians_to_degrees, string("R.A. (DEGREES)"));
    iof.write_key("CRVAL2", this->declination_radians*radians_to_degrees, string("DEC. (DEGREES)"));

    comment = "observation wavelength (meters)";
    iof.write_key("WVLNGTH", this->wavelength_meters, comment);

    this->OTF_amps.write(iof);

    fhd.write(iof);
    this->OTF_phases.write(iof);
  }

  template<class precision>
  void basic_otf<precision>::
  print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << this->unique_name()
       << "/" << "object type" << endl;
    fits_header_data<double> fhd(this->OTF_amps.get_axes());
    fhd.print(os, prefix);
    os << prefix << "RA         = " << setw(vlspc) << this->right_ascension_radians
       << "/" << "right ascension (radians)" << endl;

    os << prefix << "DEC        = " << setw(vlspc) << this->declination_radians
       << "/" << "declination (radians)" << endl;

    os << prefix << "PIXSCALE   = " << setw(vlspc) << this->pupil_plane_pixel_scale_meters
       << "/" << "pupil plane pixel scale (meters)" << endl;

    os << prefix << "WVLNGTH    = " << setw(vlspc) << this->wavelength_meters
       << "/" << "wavelength (meters)" << endl;

  }

  template<class precision>
  basic_otf<precision> & 
    basic_otf<precision>::operator*=(const basic_otf<precision> & rhs){

    if(this->get_axes()!=rhs.get_axes()){
      cerr << "basic_otf<precision>::operator*= error - mismatched axes\n";
      this->print(cerr, "lhs ");
      rhs.print(cerr, "rhs ");
      throw(string("basic_otf<precision>::operator*="));
    }
    vector<long> axes = rhs.get_axes();
    int nelems = axes[0]*axes[1];
    complex<precision> cmp;
    for(int i=0; i<nelems; i++){
      if(abs(this->data(i))!=0  && abs(rhs.data(i))!=0)
	this->set_data(i, this->data(i)*rhs.data(i));
      else
	this->set_data(i,0);
    }
  }

  template<class precision>
  basic_otf<precision> & 
    basic_otf<precision>::operator/=(const basic_otf<precision> & rhs){

    if(this->get_axes()!=rhs.get_axes()){
      cerr << "basic_otf<precision>::operator/= error - mismatched axes\n";
      this->print(cerr, "lhs ");
      rhs.print(cerr, "rhs ");
      throw(string("basic_otf<precision>::operator/="));
    }
    vector<long> axes = rhs.get_axes();
    int nelems = axes[0]*axes[1];
    for(int i=0; i<nelems; i++){
      if(abs(this->data(i))!=0)
	this->set_data(i, this->data(i)/rhs.data(i));
      else
	this->set_data(i,0);
    }
  }

  template<class precision>
  basic_otf<precision> & 
    basic_otf<precision>::operator*=(double fac){
    this->OTF_amps *= fac;
  }

  template<class precision>
  basic_otf<precision> & 
    basic_otf<precision>::operator+=(const basic_otf<precision> & rhs){

    if(this->get_axes()!=rhs.get_axes()){
      cerr << "basic_otf<precision>::operator+= error - mismatched axes\n";
      this->print(cerr, "lhs ");
      rhs.print(cerr, "rhs ");
      throw(string("basic_otf<precision>::operator+="));
    }
    vector<long> axes = rhs.get_axes();
    int nelems = axes[0]*axes[1];
    for(int i=0; i<nelems; i++)
      this->set_data(i, this->data(i)+rhs.data(i));
  }

  template<class precision>
  basic_otf<precision> & 
    basic_otf<precision>::operator-=(const basic_otf<precision> & rhs){

    if(this->get_axes()!=rhs.get_axes()){
      cerr << "basic_otf<precision>::operator+= error - mismatched axes\n";
      this->print(cerr, "lhs ");
      rhs.print(cerr, "rhs ");
      throw(string("basic_otf<precision>::operator+="));
    }
    vector<long> axes = rhs.get_axes();
    int nelems = axes[0]*axes[1];
    for(int i=0; i<nelems; i++)
      this->set_data(i, this->data(i)-rhs.data(i));
  }

  template <class precision> 
    basic_otf<precision> operator * (const basic_otf<precision> &p1, 
				     double & fac){
    basic_otf<precision> tmp(p1);
    tmp.operator*=(fac);
    return(tmp);
  }

  template <class precision> 
    basic_otf<precision> operator * (const basic_otf<precision> &p1, 
				     const basic_otf<precision> &p2){
    basic_otf<precision> tmp(p1);
    tmp.operator*=(p2);
    return(tmp);
  }

  template <class precision> 
    basic_otf<precision> operator / (const basic_otf<precision> &p1, 
				     const basic_otf<precision> &p2){
    basic_otf<precision> tmp(p1);
    tmp.operator/=(p2);
    return(tmp);
  }

  template <class precision> 
    basic_otf<precision> operator + (const basic_otf<precision> &p1, 
				     const basic_otf<precision> &p2){
    basic_otf<precision> tmp(p1);
    tmp.operator+=(p2);
    return(tmp);
  }

  template <class precision> 
    basic_otf<precision> operator - (const basic_otf<precision> &p1, 
				     const basic_otf<precision> &p2){
    basic_otf<precision> tmp(p1);
    tmp.operator-=(p2);
    return(tmp);
  }

}
#endif
