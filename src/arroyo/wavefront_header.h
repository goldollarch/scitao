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

#ifndef WAVEFRONT_HEADER_H
#define WAVEFRONT_HEADER_H

#include <iostream>
#include <iomanip>
#include "three_frame.h"
#include "region_base.h"
#include "fits_header_data.h"

namespace Arroyo {

  using std::string;
  using std::vector;
  using std::ostream;

  ///
  ///  A virtual base class for diffractive and geometric wavefront
  ///  headers

  class wavefront_header {

  protected:
  
    /// The timestamp, in seconds
    double timestamp;

  public:
  
    ///////////////////////////////////////////
    ///  Null constructor
    wavefront_header(){timestamp = 0;};

    ///////////////////////////////////////////
    ///  Virtual destructor
    virtual ~wavefront_header(){};

    ///////////////////////////////////////////
    ///  Operator = 
    wavefront_header & operator=(const wavefront_header & wfh);

    ///////////////////////////////////////////
    ///  Virtual read from iofits object
    virtual void read(const iofits & iof);
  
    ///////////////////////////////////////////
    ///  Virtual write to iofits object
    virtual void write(iofits & iof) const;

    ///////////////////////////////////////////
    ///  Virtual print
    virtual void print(ostream & os, const char * prefix="") const;

    ///////////////////////////////////////////
    ///  Virtual function to retrieve the wavelength of the radiation
    virtual double get_wavelength() const = 0;  

    ///////////////////////////////////////////
    ///  Virtual function to retrieve the wavelength of the radiation
    virtual void set_wavelength(double wvlngth) = 0;  

    ///////////////////////////////////////////
    ///  Function to return the timestamp on the wavefront header
    double get_timestamp() const {return(timestamp);};   

    ///////////////////////////////////////////
    ///  Function to set the timestamp on the wavefront header
    void set_timestamp(double tstamp) {timestamp = tstamp;};

    ///////////////////////////////////////////
    ///  Factory constructor from file
    static wavefront_header * wavefront_header_factory(const char * filename);

    ///////////////////////////////////////////
    ///  Factory constructor from iofits
    static wavefront_header * wavefront_header_factory(const iofits & iof);

    ///////////////////////////////////////////
    /// Function to compare wavefront_headers
    /// based on their timestamp
    friend int operator<(const wavefront_header & wfh1, const wavefront_header & wfh2);

    ///////////////////////////////////////////
    /// Function to compare wavefront_headers
    /// based on their timestamp
    friend int operator>(const wavefront_header & wfh1, const wavefront_header & wfh2);

    ///////////////////////////////////////////
    /// Function to compare wavefront_headers
    /// for equal timestamps
    friend bool operator==(const wavefront_header & wfh1, const wavefront_header & wfh2);

    ///////////////////////////////////////////
    ///  verbose level
    static int verbose_level;

  };

  ///////////////////////////////////////////
  /// Function to compare wavefront_headers
  /// for unequal timestamps
  bool operator!=(const wavefront_header & wfh1, const wavefront_header & wfh2);

  ///
  /// A class to hold the geometric_wavefront parameters.
  ///

  class geometric_wavefront_header :
    virtual public wavefront_header {

    protected:

    /// The number of rays
    long nrays;

    /// The wavelength of the radiation in meters
    double wavelength;

    public:

    ///////////////////////////////////////////
    ///  Null constructor
    geometric_wavefront_header();

    ///////////////////////////////////////////
    ///  Construct from an iofits object
    geometric_wavefront_header(const iofits & iof);

    ///////////////////////////////////////////
    ///  Copy constructor
    geometric_wavefront_header(const geometric_wavefront_header & gwfh);

    ///////////////////////////////////////////
    ///  Construct from the bits
    ///  wavelength in meters
    ///  pixel_scale in meters
    geometric_wavefront_header(long in_nrays, double in_wavelength);

    ///////////////////////////////////////////
    ///  Destructor
    ~geometric_wavefront_header(){};

    ///////////////////////////////////////////
    ///  Operator = 
    geometric_wavefront_header & operator=(const geometric_wavefront_header & gwfh);

    ///////////////////////////////////////////
    ///  Read geometric_wavefront_header from iofits
    void read(const iofits & iof);

    ///////////////////////////////////////////
    ///  Write geometric_wavefront_header to iofits
    void write(iofits & iof) const;

    ///////////////////////////////////////////
    ///  Print
    void print(ostream & os, const char * prefix="") const;

    ///////////////////////////////////////////
    ///  Function to retrieve the wavelength of the radiation
    double get_wavelength() const {return(wavelength);};  

    ///////////////////////////////////////////
    ///  Function to set the wavelength of the radiation.
    ///  The function argument wvlngth is in meters
    void set_wavelength(double wvlngth) {wavelength = wvlngth;};  

    ///////////////////////////////////////////
    ///  Function to retrieve the wavelength of the radiation
    long get_number_of_rays() const {return(nrays);};  

    ///////////////////////////////////////////
    ///  Friend operator ==  for geometric_wavefronts
    friend bool operator ==(const geometric_wavefront_header & gwfh1, const geometric_wavefront_header & gwfh2);

  };

  ///////////////////////////////////////////
  ///  operator !=  for geometric_wavefronts
  bool operator !=(const geometric_wavefront_header & gwfh1, const geometric_wavefront_header & gwfh2);

///
/// class to hold the diffractive_wavefront parameters.
///

  template<class T>
    class diffractive_wavefront_header :
    public three_frame,
    virtual public wavefront_header {

    protected:

    /// The wavelength of the radiation in meters
    double wavelength;

    /// The pixel scale of the array
    double pixel_scale;

    /// parameter to specify vergence of the beam.
    /// Positive curvature implies diverging beam.
    /// Zero curvature implies a plane wave.
    double curvature;

    /// To be specified in the future.
    /// A parameter to specify the conversion
    /// to units of power, or some such unit
    double radiometric_conversion;

    /// The dimensions of the array
    vector<long> axes;

    public:

    ///////////////////////////////////////////
    ///  Null constructor
    diffractive_wavefront_header();

    ///////////////////////////////////////////
    ///  Construct from an iofits object
    diffractive_wavefront_header(const iofits & iof);

    ///////////////////////////////////////////
    ///  Copy constructor
    diffractive_wavefront_header(const diffractive_wavefront_header<T> & dwfh);

    ///////////////////////////////////////////
    ///  Construct from the bits
    ///  wavelength in meters
    ///  pixel_scale in meters
    diffractive_wavefront_header(const vector<long> & in_axes, 
				 const three_frame & tf,
				 double in_wavelength, 
				 double in_pixel_scale, 
				 double in_curvature = 0);

    ///////////////////////////////////////////
    ///  Destructor
    ~diffractive_wavefront_header(){};

    ///////////////////////////////////////////
    ///  Operator = 
    diffractive_wavefront_header & operator=(const diffractive_wavefront_header<T> & dwfh);

    ///////////////////////////////////////////
    ///  Read diffractive_wavefront_header from iofits
    void read(const iofits & iof);

    ///////////////////////////////////////////
    ///  Write diffractive_wavefront_header to iofits
    void write(iofits & iof) const;

    ///////////////////////////////////////////
    ///  Print the diffractive wavefront header
    void print(ostream & os, const char * prefix="") const;

    ///////////////////////////////////////////
    ///  Function to retrieve the wavelength of the radiation
    double get_wavelength() const {return(wavelength);};  

    ///////////////////////////////////////////
    ///  Function to set the wavelength of the radiation.
    ///  The function argument wvlngth is in meters
    void set_wavelength(double wvlngth) {wavelength = wvlngth;};  

    ///////////////////////////////////////////
    ///  
    rectangular_region get_covering_region(const three_frame & tf, bool foreshortening) const;  

    ///////////////////////////////////////////
    ///  Returns the number of elements in the data array
    long total_space() const;

    ///////////////////////////////////////////
    ///  Returns the axes 
    vector<long> get_axes() const {return(axes);};

    ///////////////////////////////////////////
    ///  Sets the axes 
    void set_axes(const vector<long> & in_axes);

    ///////////////////////////////////////////
    ///  Returns the pixel scale of the wavefront
    double get_pixel_scale() const {return pixel_scale;};

    ///////////////////////////////////////////
    ///  Returns the pixel scale of the wavefront
    void set_pixel_scale(double pixscale) {
      if(pixscale<=0){
	cerr << "diffractive_wavefront::set_pixel_scale error - invalid pixel scale " << pixscale << endl;
	throw(string("diffractive_wavefront::set_pixel_scale"));
      }
      pixel_scale=pixscale;
    };

    ///////////////////////////////////////////
    ///  Returns the curvature of the wavefront, in 1/meters
    double get_curvature() const {return curvature;};

    ///////////////////////////////////////////
    ///  Sets the curvature of the wavefront
    ///
    ///  Units of in_curvature are 1/meters
    void set_curvature(double in_curvature) {curvature = in_curvature;};

    ///////////////////////////////////////////
    ///  Friend operator==  for diffractive wavefront headers
    friend bool operator==(const diffractive_wavefront_header<T> & dwfh1, const diffractive_wavefront_header<T> & dwfh2){

      if(dwfh1.timestamp!=dwfh2.timestamp) return(false);
      if(operator!=(static_cast<const three_frame>(dwfh1), 
		    static_cast<const three_frame>(dwfh2))) 
	return(false);
      if(dwfh1.wavelength!=dwfh2.wavelength) return(false);
      if(dwfh1.pixel_scale!=dwfh2.pixel_scale) return(false);
      if(dwfh1.curvature!=dwfh2.curvature) return(false);
      if(dwfh1.axes!=dwfh2.axes) return(false);
      return(true);
    };

  };

  ///////////////////////////////////////////
  ///  Operator!=  for diffractive wavefront headers
  template<class T>
    bool operator !=(const diffractive_wavefront_header<T> & dwfh1, const diffractive_wavefront_header<T> & dwfh2){
    return(!(dwfh1==dwfh2));
  }

  template<class T>
    diffractive_wavefront_header<T>::diffractive_wavefront_header() {
    wavelength = 0;
    pixel_scale = 0;
    curvature = 0;
    // ensure that the fits_header_data data
    // member axes has the correct dimensionality
    axes = vector<long>(2,0);
  }

  template<class T>
    diffractive_wavefront_header<T>::diffractive_wavefront_header(const iofits & iof) {
    // ensure that the fits_header_data data
    // member axes has the correct dimensionality
    axes = vector<long>(2,0);
    this->read(iof);
  }

  template<class T>
    diffractive_wavefront_header<T>::diffractive_wavefront_header(const diffractive_wavefront_header & wfh) {
    // ensure that the fits_header_data data
    // member axes has the correct dimensionality
    axes = vector<long>(2,0);
    this->operator=(wfh);
  }

  template<class T>
    diffractive_wavefront_header<T>::diffractive_wavefront_header(const vector<long> & in_axes, 
								  const three_frame & tf,
								  double in_wavelength, 
								  double in_pixel_scale,
								  double in_curvature){
    if(in_axes.size()!=2){
      std::cerr << "diffractive_wavefront_header::diffractive_wavefront_header error - "
		<< "axes supplied to this constructor have dimensionality " << in_axes.size()
		<< " rather than dimensionality 2\n";
      throw(string("diffractive_wavefront_header::diffractive_wavefront_header"));
    }
    if(in_axes[0]<0 || in_axes[1]<0){
      std::cerr << "diffractive_wavefront_header::diffractive_wavefront_header error - "
		<< "axes do not have positive definite values\n";
      throw(string("diffractive_wavefront_header::diffractive_wavefront_header"));
    }
  
    axes = in_axes;
    this->three_frame::operator=(tf);
    wavelength = in_wavelength;  
    pixel_scale = in_pixel_scale;
    curvature = in_curvature;
  }

  template<class T>
    diffractive_wavefront_header<T> & diffractive_wavefront_header<T>::operator=(const diffractive_wavefront_header<T> & wfh) {
    if(this==&wfh)
      return(*this);
    wavelength = wfh.wavelength;
    pixel_scale = wfh.pixel_scale;
    curvature = wfh.curvature;
    this->three_frame::operator=(wfh);  
    this->wavefront_header::operator=(wfh);
    axes = wfh.axes;
    return(*this);
  }

  template<class T>
    void diffractive_wavefront_header<T>::read(const iofits & iof) {
    if(!iof.key_exists("TYPE")){
      std::cerr << "diffractive_wavefront_header::read error - "
	   << "unrecognized type of file\n";
      throw(std::string("diffractive_wavefront::read"));
    }
    std::string type, comment;
    iof.read_key("TYPE", type, comment);
    if(type!="diffractive wavefront"){
      std::cerr << "diffractive_wavefront_header::read error - file of type " 
	   << type << " rather than of type diffractive_wavefront\n";
      throw(std::string("diffractive_wavefront::read"));
    }

    Arroyo::fits_header_data<T> tmphdr(iof);
    axes = tmphdr.get_axes();

    this->wavefront_header::read(iof);
    this->three_frame::read(iof);
    iof.read_key("WVLNGTH", wavelength, comment);
    iof.read_key("PIXSCALE", pixel_scale, comment);
    iof.read_key("CURVATUR", curvature, comment);
  }

  template<class T>
    void diffractive_wavefront_header<T>::write(iofits & iof) const {
  
    fits_header_data<T> fhd(this->get_axes()); 

    fhd.write(iof);
    string type = "diffractive wavefront";
    string comment = "object type";
    iof.write_key("TYPE", type, comment);
    this->wavefront_header::write(iof);
    this->three_frame::write(iof);
    iof.write_key("WVLNGTH", wavelength, "wavelength (microns)");
    iof.write_key("PIXSCALE", pixel_scale, "pixel scale (meters)");
    iof.write_key("CURVATUR", curvature, "curvature (1/meters)");
  }

  template<class T>
    void diffractive_wavefront_header<T>::print(ostream & os, const char * prefix) const {

    int vlspc = 30;
    os.setf(std::ios::left, std::ios::adjustfield); 
    os << prefix << "TYPE       = " << std::setw(vlspc) << "diffractive wavefront"
       << "/" << "object type" << std::endl;
    fits_header_data<T> fhd(axes); 
    fhd.print(os, prefix);
    this->wavefront_header::print(os, prefix);
    three_frame::print(os, prefix);
    os << prefix << "WVLNGTH    = " << std::setw(vlspc) << wavelength*1e6
       << "/" << "wavelength (microns)" << std::endl;
    os << prefix << "PIXSCALE   = " << std::setw(vlspc) << pixel_scale
       << "/" << "pixel scale (meters)" << std::endl;
    os << prefix << "CURVATUR   = " << std::setw(vlspc) << curvature
       << "/" << "curvature (1/meters)" << std::endl;
  }

  template<class T>
    long diffractive_wavefront_header<T>::total_space() const {
    if(axes.size()==0) return(0);
    int nelem = 1;
    for(uint i=0; i<axes.size(); i++)
      nelem *= axes[i];
    return(nelem);
  }

  template<class T>
    void diffractive_wavefront_header<T>::set_axes(const vector<long> & in_axes) {

    if(in_axes.size()!=2){
      cerr << "diffractive_wavefront_header::set_axes error - "
	   << "axes have dimension " << axes.size() << " rather than dimension 2\n";
      throw(string("diffractive_wavefront_header::set_axes")); 
    }
    if(in_axes[0]<=0 || in_axes[1]<=0){
      cerr << "diffractive_wavefront_header::set_axes error - "
	   << "dimensions " << axes[0] << "x" << axes[1] 
	   << "provided to this function are not positive\n";
      throw(string("diffractive_wavefront_header::set_axes")); 
    } 
    axes = in_axes;
  } 

//#################################### gdc ####################################################### 
  three_point get_ray_plane_intersection(const three_point & o_a, const three_vector & n_a,
    const three_point & o_b, const three_vector & n_b);
//#################################### gdc ####################################################### 

  template<class T>
    rectangular_region diffractive_wavefront_header<T>::get_covering_region(const three_frame & tf, 
									    bool foreshortening) const {
    
    if(axes[0]==0 || axes[1]==0){
      cerr << "diffractive_wavefront_header::get_covering_region error - "
	   << "wavefront axes " << axes[0] << "x" << axes[1] 
	   << " cannot be used to determine a covering region\n";
      throw(string("diffractive_wavefront_header::get_covering_region"));
    }

    if(fabs(dot_product(tf.z(), this->z()))<three_frame::precision){
      cerr << "diffractive_wavefront_header::get_covering_region error - "
	   << "z axes of this aperture and three frame provided to this function are orthogonal\n";
      throw(string("diffractive_wavefront_header::get_covering_region"));
    }

    vector<double> size(2, axes[0]*pixel_scale);
    size[1] = axes[1]*pixel_scale;
    vector<three_point> tp(4);
    tp[0] = *this + .5*size[0]*this->x() + .5*size[1]*this->y();
    tp[1] = *this - .5*size[0]*this->x() + .5*size[1]*this->y();
    tp[2] = *this - .5*size[0]*this->x() - .5*size[1]*this->y();
    tp[3] = *this + .5*size[0]*this->x() - .5*size[1]*this->y();

    if(foreshortening){
      // add a component along this->z() to get
      // a three point in the XY plane of tf
      double distance;
      for(int i=0; i<4; i++)
	      tp[i] = get_ray_plane_intersection(tp[i], this->z(), tf, tf.z());

    } else {
      Arroyo::three_vector rotation_axis = cross_product(this->z(), tf.z());
      if(rotation_axis.length() > three_frame::precision) {
	Arroyo::three_rotation trot(*this, rotation_axis, rotation_axis.length());
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
}

#endif
