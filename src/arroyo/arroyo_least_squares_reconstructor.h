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

#ifndef ARROYO_LEAST_SQUARES_RECONSTRUCTOR_H
#define ARROYO_LEAST_SQUARES_RECONSTRUCTOR_H

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <iostream>
#include "region_base.h"
#include "linear_algebra.h"
#include "zernike_projected_zonal_reconstructor.h"
#include "aperture.h"
#include "lenslet_array.h"

namespace Arroyo {
  ///
  /// A class to hold least squares reconstructors generated using Arroyo
  ///

  template<class T>
  class arroyo_least_squares_reconstructor :
    public zernike_projected_zonal_reconstructor,
    public Arroyo::pixel_array<T> {

    private:

    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    std::string unique_name() const {return(std::string("arroyo least squares reconstructor"));};

    protected:

    /// Actuator axes
    vector<long> actuator_axes;

    /// Lenslet axes
    vector<long> centroid_axes;

    /// Projected modes
    zernike projected_modes;

    // A flag to indicate whether the eigenmodes were preserved
    // when the reconstructor was generated
    bool eigenmodes_stored_in_instance;

    // The eigenvalue threshold used in forming the reconstructor
    double eigenvalue_threshold;

    // A flag to indicate whether the geometry matrix was preserved
    // when the reconstructor was generated
    bool geometry_matrix_stored_in_instance;

    // A flag to indicate whether the A++ sign convention is used
    bool app_sign_convention;

    // The subap illumination threshold used in forming the reconstructor
    double subaperture_illumination_threshold;

    /// pixel array to contain the geometry matrix G
    pixel_array<T> geometry_matrix;

    /// pixel array to contain the eigenmodes of the matrix GG^{T}
    pixel_array<T> g_gtranspose_eigenmodes;
    
    /// pixel array to contain the matrix G^{T}G
    pixel_array<T> gtranspose_g_eigenmodes;

    /// pixel array to contain the square roots of the eigenvalues
    pixel_array<T> sqrt_eigenvalues;


    ///////////////////////////////////////////
    ///  Null constructor
    arroyo_least_squares_reconstructor(){};

    public:

    ///////////////////////////////////////////
    ///  Copy constructor
    arroyo_least_squares_reconstructor(const arroyo_least_squares_reconstructor & arroyo_lsq_recon);

    ///////////////////////////////////////////
    ///  Construct from file
    arroyo_least_squares_reconstructor(const char * filename);

    ///////////////////////////////////////////
    ///  Construct from iofits object
    arroyo_least_squares_reconstructor(const Arroyo::iofits & iof);

    ///////////////////////////////////////////
    ///  Construct a least squares reconstructor from an aperture, an
    ///  ideal deformable mirror, and a square lenslet array.  The
    ///  influence function of the ideal deformable mirror is assumed
    ///  to be pyramidal.
    ///
    ///  The aperture dimensions, deformable mirror actuator pitch,
    ///  and the lenslet pitch are used together with the three frames
    ///  from each of these objects to define the geometrical
    ///  registration of the pupil, actuators and subapertures.
    ///  Arbitrary misregistration, magnification, tip and tilt, and
    ///  rotation are handled by this constructor.  A geometry matrix
    ///  for the configuration is generated and then inverted using the
    ///  singular value decomposition algortihm from LAPACK.
    ///
    ///  NOTE - currently only Fried geometry is supported
    ///
    ///  The subaperture illumination threshold argument to this
    ///  function allows one to zero geometry matrix elements
    ///  corresponding to actuator influence functions and
    ///  subapertures that overlap by less than the threshold value.
    ///  This value should lie between zero and one.
    ///
    ///  The eigenvalue threshold argument to this function allows one
    ///  to remove modes from the reconstructor with eigenvalues lower
    ///  than the threshold.  This argument is interpreted as the
    ///  ratio between the threshold eigenvalue and the largest
    ///  eigenvalue, and its value should lie between zero and one.
    template<class aperture_type, class dm_ap>
    arroyo_least_squares_reconstructor(const aperture_type & ap, 
				       const ideal_deformable_mirror<dm_ap> & idm,
				       const square_lenslet_array & sq_lnslt_arr,
				       const zernike & projection_modes,
				       bool modal_projection_in_actuator_space,
				       double subaperture_illumination_threshold,
				       double eigenvalue_threshold,
				       bool aplusplus_sign_convention = false,
				       bool store_eigenmodes = false,
				       bool store_geometry_matrix = false);

    ///////////////////////////////////////////
    ///  Destructor
    ~arroyo_least_squares_reconstructor(){};

    ///////////////////////////////////////////
    ///  Operator = 
    arroyo_least_squares_reconstructor & 
      operator=(const arroyo_least_squares_reconstructor & arroyo_lsq_recon);

     ///////////////////////////////////////////
    ///  Read from file
    void read(const char * filename);

    ///////////////////////////////////////////
    ///  Read from iofits
    void read(const Arroyo::iofits & iof);

    ///////////////////////////////////////////
    ///  Write to file
    void write(const char * filename) const;
 
    ///////////////////////////////////////////
    ///  Write to iofits
    void write(Arroyo::iofits & iof) const;
 
    ///////////////////////////////////////////
    ///  Print
    void print(std::ostream & os, const char * prefix="") const;

    ///////////////////////////////////////////
    ///  Get the eigenvalue threshold used 
    ///  in constructing the reconstructor
    double get_eigenvalue_threshold() const;

    ///////////////////////////////////////////
    ///  Get the subaperture illumination threshold used in
    ///  constructing the reconstructor
    double get_subaperture_illumination_threshold() const;

    ///////////////////////////////////////////
    ///  Get dimensions of centroid measurements passed
    ///  to the reconstructor
    vector<long> get_centroid_axes() const;

    ///////////////////////////////////////////
    ///  Get dimensions of pixel array returned
    ///  by the reconstructor
    vector<long> get_actuator_axes() const;

    ///////////////////////////////////////////
    ///  Get a zernike instance that contains information
    ///  about which modes are returned by the reconstructor.
    /// 
    ///  This instance is minimally sized so as to hold the largest
    ///  mode returned by the reconstructor.  Each element of this
    ///  instance is initialized to unity if the corresponding mode is
    ///  returned by the reconstructor, and to zero if it is not
    Arroyo::zernike get_zernike_modes() const;

    ///////////////////////////////////////////
    ///  Indicate whether the geometry matrix was stored
    bool geometry_matrix_stored() const;

    ///////////////////////////////////////////
    ///  Indicate whether the geometry matrix was stored
    bool eigenmodes_stored() const;

    ///////////////////////////////////////////
    ///  Retrieve the pixel array containing the geometry matrix G
    ///
    ///  If the reconstructor was constructed without saving the
    ///  geometry matrix, this function throws an error.
    ///  
    pixel_array<T> get_geometry_matrix() const;

    ///////////////////////////////////////////
    ///  Get the number of actuator eigenmodes
    ///
    ///  If the reconstructor was constructed without saving the
    ///  eigenmodes, this function throws an error.
    int get_nactuator_eigenmodes() const;

    ///////////////////////////////////////////
    ///  Get the number of centroid eigenmodes
    ///
    ///  If the reconstructor was constructed without saving the
    ///  eigenmodes, this function throws an error.
    int get_ncentroid_eigenmodes() const;

    ///////////////////////////////////////////
    ///  Retrieve the eigenvalue for a particular mode
    ///
    ///  If the reconstructor was constructed without saving the
    ///  eigenmodes, this function throws an error.
    ///  
    ///  If the mode number is out of range, this function throws an
    ///  error.  Note - number of eigenvalues is equal to the number
    ///  of actuator eigenmodes
    double get_eigenvalue(int mode_number) const;

    ///////////////////////////////////////////
    ///  Retrieve an actuator eigenmode.  This function returns the
    ///  eigenvalue of the mode.
    ///
    ///  If the reconstructor was constructed without saving the
    ///  eigenmodes, this function throws an error.
    ///  
    ///  If the mode number is out of range, this function throws an
    ///  error.
    template<class U>
    double get_actuator_eigenmode(int mode_number, 
				  pixel_array<U> & actuator_mode) const;

    ///////////////////////////////////////////
    ///  Retrieve a centroid eigenmode.  This function returns the
    ///  eigenvalue of the mode.
    ///
    ///  If the reconstructor was constructed without saving the
    ///  eigenmodes, this function throws an error.
    ///  
    ///  If the mode number is out of range, this function throws an
    ///  error.
    double get_centroid_eigenmode(int mode_number, 
				  Shack_Hartmann_centroids & centroid_mode) const;

    ///////////////////////////////////////////
    ///  Reconstruct the zernike residuals from 
    ///  Shack Hartmann centroid data
    ///
    /// This reconstructor reconstructs a zernike instance
    /// with tip and tilt modes
    void reconstruct_zernike_residuals(const Arroyo::Shack_Hartmann_centroids & shcentroids, 
				       Arroyo::zernike & znke) const;

    ///////////////////////////////////////////
    ///  Reconstruct the zonal residuals from a
    ///  Shack Hartmann centroid class instance
    void reconstruct_zonal_residuals(const Arroyo::Shack_Hartmann_centroids & shcentroids, 
				     Arroyo::pixel_array<double> & pixarr) const;

    ///////////////////////////////////////////
    ///  Reconstruct the residuals from 
    ///  Shack Hartmann centroid data
    ///
    /// This reconstructor reconstructs a zernike instance
    /// with tip and tilt modes, and a pixel array containing
    /// zonal residuals
    void reconstruct_residuals(const Arroyo::Shack_Hartmann_centroids & shcentroids, 
			       Arroyo::zernike & znke, 
			       Arroyo::pixel_array<double> & pixarr) const;

    ///////////////////////////////////////////
    //  Compute the set of centroids consistent with 
    //  the phase map provided in the argument
    //  to this function.  The slope discrepancy
    //  may be computed by subtracting these
    //  centroids from the set of centroids used 
    //  to reconstruct the phase map.
    //
    //  This operation utilizes the geometry matrix,
    //  and the function throws an error if this 
    //  matrix was not preserved during construction
    //  of this instance.
    Arroyo::Shack_Hartmann_centroids get_slope_consistency(Arroyo::pixel_array<double> & pixarr) const;

  };

  template<class T> 
  arroyo_least_squares_reconstructor<T>::
  arroyo_least_squares_reconstructor(const arroyo_least_squares_reconstructor & arroyo_lsq_recon){
    this->operator=(arroyo_lsq_recon);
  }

  template<class T> 
  arroyo_least_squares_reconstructor<T>::
  arroyo_least_squares_reconstructor(const char * filename){
    this->axes.resize(2);
    this->read(filename);
  }

  template<class T> 
  arroyo_least_squares_reconstructor<T>::
  arroyo_least_squares_reconstructor(const iofits & iof){
    this->axes.resize(2);
    this->read(iof);
  }

  template<class T> 
  arroyo_least_squares_reconstructor<T> & 
  arroyo_least_squares_reconstructor<T>::
  operator=(const arroyo_least_squares_reconstructor & arroyo_lsq_recon){
    if(this==&arroyo_lsq_recon)
      return(*this);
    actuator_axes = arroyo_lsq_recon.actuator_axes;
    centroid_axes = arroyo_lsq_recon.centroid_axes;
    projected_modes = arroyo_lsq_recon.projected_modes;
    eigenmodes_stored_in_instance = arroyo_lsq_recon.eigenmodes_stored_in_instance;
    geometry_matrix_stored_in_instance = arroyo_lsq_recon.geometry_matrix_stored_in_instance;
    geometry_matrix = arroyo_lsq_recon.geometry_matrix;
    app_sign_convention = arroyo_lsq_recon.app_sign_convention;
    g_gtranspose_eigenmodes = arroyo_lsq_recon.g_gtranspose_eigenmodes; 
    gtranspose_g_eigenmodes = arroyo_lsq_recon.gtranspose_g_eigenmodes; 
    sqrt_eigenvalues = arroyo_lsq_recon.sqrt_eigenvalues;

    this->pixel_array<T>::operator=(arroyo_lsq_recon);
    return(*this);
  }

  template<class T> 
  void arroyo_least_squares_reconstructor<T>::read(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "arroyo_least_squares_reconstructor::read - "
	   << "error opening file " << filename << endl;
      throw(string("arroyo_least_squares_reconstructor::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "arroyo_least_squares_reconstructor::read - "
	   << "error reading "
	   << this->unique_name() << " from file " 
	   << filename << endl;
      throw(string("arroyo_least_squares_reconstructor::read"));
    }
  }

  template<class T> 
  void arroyo_least_squares_reconstructor<T>::read(const iofits & iof) {
    if(!iof.key_exists("TYPE")){
      cerr << "arroyo_least_squares_reconstructor::read error - "
	   << "unrecognized type of file\n";
      throw(string("arroyo_least_squares_reconstructor::read"));
    }
    string type, comment;
    iof.read_key("TYPE", type, comment);
    if(type!=this->unique_name()){
      cerr << "arroyo_least_squares_reconstructor::read error - file of type " 
	   << type << " rather than of type "
	   << this->unique_name() << endl;
      throw(string("arroyo_least_squares_reconstructor::read"));
    }
    
    actuator_axes.resize(2);
    centroid_axes.resize(2);

    iof.read_key("ACTUATR1", actuator_axes[0], comment);
    iof.read_key("ACTUATR2", actuator_axes[1], comment);

    iof.read_key("CNTROID1", centroid_axes[0], comment);
    iof.read_key("CNTROID2", centroid_axes[1], comment);

    iof.read_key("EGNTHRSH", eigenvalue_threshold, comment);
    iof.read_key("SBPTHRSH", subaperture_illumination_threshold, comment);

    if(iof.key_exists("EIGENMDS"))
       iof.read_key("EIGENMDS", eigenmodes_stored_in_instance, comment);
    else eigenmodes_stored_in_instance = false;

    if(iof.key_exists("GEOMETRY"))
      iof.read_key("GEOMETRY", geometry_matrix_stored_in_instance, comment);
    else geometry_matrix_stored_in_instance = false;

    if(iof.key_exists("APPSIGN"))
      iof.read_key("APPSIGN", app_sign_convention, comment);
    else app_sign_convention = false;

    this->pixel_array<T>::read(iof);

    // Work out the number of projected zernike modes, and initialize
    // the projected_modes data member
    int nelem = this->pixel_array<T>::total_space();
    nelem -= actuator_axes[0]*actuator_axes[1]*centroid_axes[0]*centroid_axes[1];
    if(nelem%(centroid_axes[0]*centroid_axes[1])!=0){
      cerr << "arroyo_least_squares_reconstructor::read error - "
	   << "unexpected number of elements " << nelem << " in the matrix\n";
      throw(string("arroyo_least_squares_reconstructor::read"));
    }
      
    int nmodes = nelem/(centroid_axes[0]*centroid_axes[1]);
    int zernike_order = -1;
    while(nmodes>0){
      zernike_order++;
      if(zernike_order%2)
	nmodes -= 2*(zernike_order/2+1);
      else 
	nmodes -= (zernike_order+1);
    }

    if(zernike_order>=0)
      projected_modes = zernike(zernike_order);

    if(geometry_matrix_stored_in_instance){
      try{
	geometry_matrix.read(iof);
      } catch(...) {
	cerr << "arroyo_least_squares_reconstructor::read - error reading geometry matrix\n";
	throw(string("arroyo_least_squares_reconstructor::read"));
      }
    }

    if(eigenmodes_stored_in_instance){

      try{
	sqrt_eigenvalues.read(iof);
      } catch(...) {
	cerr << "arroyo_least_squares_reconstructor::read - error reading eigenvalues\n";
	throw(string("arroyo_least_squares_reconstructor::read"));
      }

      try{
	g_gtranspose_eigenmodes.read(iof);
	gtranspose_g_eigenmodes.read(iof);
      } catch(...) {
	cerr << "arroyo_least_squares_reconstructor::read - error reading GGtranspose matrices\n";
	throw(string("arroyo_least_squares_reconstructor::read"));
      }
    }

    // leave iof pointing to the next header, if it exists
    if(iof.get_hdu_num()!=iof.get_num_hdus())
      iof.movrel_hdu(1);
  }

  template<class T> 
  void arroyo_least_squares_reconstructor<T>::write(const char * filename) const {

    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "arroyo_least_squares_reconstructor::write - "
	   << "error opening file " << filename << endl;
      throw(string("arroyo_least_squares_reconstructor::write"));
    }

    try{this->write(iof);}
    catch(...){
      cerr << "arroyo_least_squares_reconstructor::write - "
	   << "error writing "
	   << this->unique_name() << " to file "
	   << filename << endl;
      throw(string("arroyo_least_squares_reconstructor::write"));
    }
  }

  template<class T> 
  void arroyo_least_squares_reconstructor<T>::write(iofits & iof) const {
    fits_header_data<T> fhd(this->axes);
    fhd.write(iof);
    string type = this->unique_name();
    string comment = "object type";
    iof.write_key("TYPE", type, comment);

    iof.write_key("ACTUATR1", actuator_axes[0], "length of actuator axis 1");
    iof.write_key("ACTUATR2", actuator_axes[1], "length of actuator axis 2");

    iof.write_key("CNTROID1", centroid_axes[0], "length of centroid axis 1");
    iof.write_key("CNTROID2", centroid_axes[1], "length of centroid axis 2");

    iof.write_key("EGNTHRSH", eigenvalue_threshold, "threshold for keeping eigenvalues");
    iof.write_key("SBPTHRSH", subaperture_illumination_threshold, "threshold for keeping partially illuminated subaps");

    iof.write_key("APPSIGN", app_sign_convention, "A++ sign convention");

    if(sqrt_eigenvalues.get_axes().size())
      iof.write_key("EIGENMDS", true, "eigenmodes stored");
    else 
      iof.write_key("EIGENMDS", false, "eigenmodes stored");

    if(geometry_matrix.get_axes().size())
      iof.write_key("GEOMETRY", true, "geometry matrix stored");
    else
      iof.write_key("GEOMETRY", false, "geometry matrix stored");

    this->pixel_array<T>::write(iof);

    if(geometry_matrix.get_axes().size()){
      fits_header_data<T> fhd(geometry_matrix.get_axes());
      fhd.write(iof);
      geometry_matrix.write(iof);
    }
    
    if(sqrt_eigenvalues.get_axes().size()){
      fits_header_data<T> fhd1(sqrt_eigenvalues.get_axes());
      fhd1.write(iof);
      sqrt_eigenvalues.write(iof);

      fits_header_data<T> fhd2(g_gtranspose_eigenmodes.get_axes());
      fhd2.write(iof);
      g_gtranspose_eigenmodes.write(iof);

      fits_header_data<T> fhd3(gtranspose_g_eigenmodes.get_axes());
      fhd3.write(iof);
      gtranspose_g_eigenmodes.write(iof);
    }
  }

  template<class T> 
  void arroyo_least_squares_reconstructor<T>::print(ostream & os, const 
						 char * prefix) const {
    int vlspc = 30;
    os.setf(std::ios::left, std::ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << this->unique_name()
       << "/" << "object type" << endl;
    fits_header_data<T> fhd(this->axes);
    fhd.print(os, prefix);
    os << prefix << "ACTUATR1   = " << setw(vlspc) << this->actuator_axes[0]
       << "/" << "length of actuator axis 1" << endl;
    os << prefix << "ACTUATR2   = " << setw(vlspc) << this->actuator_axes[1]
       << "/" << "length of actuator axis 2" << endl;
    os << prefix << "CNTROID1   = " << setw(vlspc) << this->centroid_axes[0]
       << "/" << "length of centroid axis 1" << endl;
    os << prefix << "CNTROID2   = " << setw(vlspc) << this->centroid_axes[1]
       << "/" << "length of centroid axis 2" << endl;
    os << prefix << "EGNTHRSH   = " << setw(vlspc) << this->eigenvalue_threshold
       << "/" << "threshold for keeping eigenvalues" << endl;
    os << prefix << "APPSIGN    = " << setw(vlspc) << this->app_sign_convention
       << "/" << "A++ sign convention" << endl;
    os << prefix << "SBPTHRSH   = " << setw(vlspc) << this->subaperture_illumination_threshold
       << "/" << "threshold for keeping partially illuminated subaps" << endl;
    os << prefix << "GEOMETRY   = " << setw(vlspc) << this->geometry_matrix_stored_in_instance
       << "/" << "geometry matrix storage flag" << endl;
    os << prefix << "EIGENMDS   = " << setw(vlspc) << this->eigenmodes_stored_in_instance
       << "/" << "eigenmodes storage flag" << endl;
  }

  template<class T> 
  vector<long> arroyo_least_squares_reconstructor<T>::get_centroid_axes() const {
    return(centroid_axes);
  }

  template<class T> 
  vector<long> arroyo_least_squares_reconstructor<T>::get_actuator_axes() const {
    return(actuator_axes);
  }

  template<class T> 
  zernike arroyo_least_squares_reconstructor<T>::get_zernike_modes() const {
    return(projected_modes);
  }

  template<class T> 
  double arroyo_least_squares_reconstructor<T>::get_eigenvalue_threshold() const {
    return(eigenvalue_threshold);
  }

  template<class T> 
  double arroyo_least_squares_reconstructor<T>::get_subaperture_illumination_threshold() const {
    return(subaperture_illumination_threshold);
  }

  template<class T> 
  int arroyo_least_squares_reconstructor<T>::get_nactuator_eigenmodes() const {
    return(actuator_axes[0]*actuator_axes[1]);
  }

  template<class T> 
  int arroyo_least_squares_reconstructor<T>::get_ncentroid_eigenmodes() const {
    return(centroid_axes[0]*centroid_axes[1]);
  }

  template<class T>
    bool arroyo_least_squares_reconstructor<T>::
    geometry_matrix_stored() const {
    return(this->geometry_matrix_stored_in_instance);
  }

  template<class T>
    bool arroyo_least_squares_reconstructor<T>::
    eigenmodes_stored() const {
    return(this->eigenmodes_stored_in_instance);
  }

  template<class T>
    pixel_array<T> arroyo_least_squares_reconstructor<T>::
    get_geometry_matrix() const {
    if(this->geometry_matrix_stored())
      return(this->geometry_matrix);
    else {
      cerr << "arroyo_least_squares_reconstructor::get_geometry_matrix error - " 
	   << "this class instance generated without saving geometry matrix\n";
      throw(string("arroyo_least_squares_reconstructor::get_geometry_matrix"));
    }
  }

  template<class T> 
  double arroyo_least_squares_reconstructor<T>::
    get_eigenvalue(int mode_number) const {

    if(!this->eigenmodes_stored_in_instance){
      cerr << "arroyo_least_squares_reconstructor::get_eigenvalue error - "
	   << " this class instance was constructed without saving the eigenmodes\n";
      throw(string("arroyo_least_squares_reconstructor::get_eigenvalue"));
    }

    if(mode_number<0 || mode_number>this->get_nactuator_eigenmodes()){
      cerr << "arroyo_least_squares_reconstructor::get_eigenvalue error - "
	   << " mode number " << mode_number << " out of range\n";
      throw(string("arroyo_least_squares_reconstructor::get_eigenvalue"));
    }

    return(sqrt_eigenvalues.data(mode_number)*sqrt_eigenvalues.data(mode_number));
  }

  template<class T> 
    template<class U>
  double arroyo_least_squares_reconstructor<T>::
    get_actuator_eigenmode(int mode_number, pixel_array<U> & actuators) const {

    if(!this->eigenmodes_stored_in_instance){
      cerr << "arroyo_least_squares_reconstructor::get_actuator_eigenmode error - "
	   << " this class instance was constructed without saving the eigenmodes\n";
      throw(string("arroyo_least_squares_reconstructor::get_actuator_eigenmode"));
    }

    if(mode_number<0 || mode_number>this->get_nactuator_eigenmodes()){
      cerr << "arroyo_least_squares_reconstructor::get_actuator_eigenmode error - "
	   << " mode number " << mode_number << " out of range\n";
      throw(string("arroyo_least_squares_reconstructor::get_actuator_eigenmode"));
    }

    if(actuators.get_axes()!=actuator_axes){
      cerr << "arroyo_least_squares_reconstructor::get_actuator_eigenmode error - "
	   << "pixel array instance passed to this function "
	   << "has a different dimensionality than expected\n";
      cerr << "\texpected " << actuator_axes[0] << "x" << actuator_axes[1] << endl;
      cerr << "\tgot " << actuators.get_axes()[0] << "x" << actuators.get_axes()[1] << endl;
      throw(string("arroyo_least_squares_reconstructor::get_actuator_eigenmode"));
    }

    int nactuators = actuator_axes[0]*actuator_axes[1];
    for(int i=0; i<nactuators; i++)
      actuators.set_data(i, gtranspose_g_eigenmodes.data(mode_number+i*nactuators));

    return(sqrt_eigenvalues.data(mode_number)*sqrt_eigenvalues.data(mode_number));
  }

  template<class T> 
  double arroyo_least_squares_reconstructor<T>::
    get_centroid_eigenmode(int mode_number, Shack_Hartmann_centroids & shcentroids) const {

    if(!this->eigenmodes_stored_in_instance){
      cerr << "arroyo_least_squares_reconstructor::get_centroid_eigenmode error - "
	   << " this class instance was constructed without saving the eigenmodes\n";
      throw(string("arroyo_least_squares_reconstructor::get_centroid_eigenmode"));
    }

    if(mode_number<0 || mode_number>this->get_ncentroid_eigenmodes()){
      cerr << "arroyo_least_squares_reconstructor::get_centroid_eigenmode error - "
	   << " mode number " << mode_number << " out of range\n";
      throw(string("arroyo_least_squares_reconstructor::get_centroid_eigenmode"));
    }

    if(shcentroids.get_axes()!=centroid_axes){
      cerr << "arroyo_least_squares_reconstructor::get_centroid_eigenmode error - "
	   << "Shack Hartmann centroid instance passed to this function "
	   << "has a different dimensionality than expected\n";
      cerr << "\texpected " << centroid_axes[0] << "x" << centroid_axes[1] << endl;
      cerr << "\tgot " << shcentroids.get_axes()[0] << "x" << shcentroids.get_axes()[1] << endl;
      throw(string("arroyo_least_squares_reconstructor::get_centroid_eigenmode"));
    }

    int ncentroids = centroid_axes[0]*centroid_axes[1];
    for(int i=0; i<ncentroids; i++)
      shcentroids.set_data(i, g_gtranspose_eigenmodes.data(mode_number*ncentroids+i));

    if(mode_number<sqrt_eigenvalues.get_axes()[0])
      return(sqrt_eigenvalues.data(mode_number)*sqrt_eigenvalues.data(mode_number));
    else 
      return(0);
  }

  template<class T> 
  void arroyo_least_squares_reconstructor<T>::
  reconstruct_zernike_residuals(const Arroyo::Shack_Hartmann_centroids & shcentroids, 
				zernike & znke) const {

    if(shcentroids.total_space()!=this->axes[0]){
      cerr << "arroyo_least_squares_reconstructor::reconstruct_zernike_residuals error - "
	   << " Shack Hartmann centroids instance supplied to this function has " 
	   << shcentroids.total_space() << " elements, whereas " 
	   << this->axes[0] << " were expected by the reconstructor\n";
      throw(string("arroyo_least_squares_reconstructor::reconstruct_zernike_residuals"));
    }

    if(znke.get_order()!=projected_modes.get_order()){
      cerr << "arroyo_least_squares_reconstructor::reconstruct_zernike_residuals error - "
	   << " zernike order supplied to this function " 
	   << znke.get_order() << " does not match that expected by the reconstructor "
	   << projected_modes.get_order() << endl;
      throw(string("arroyo_least_squares_reconstructor::reconstruct_zernike_residuals"));
    }

    // For now, this is implemented for tip and tilt only, as the constructor only handles
    // this case.

    int xindex = centroid_axes[0]*centroid_axes[1]*(actuator_axes[0]*actuator_axes[1]+1);
    int yindex = xindex + centroid_axes[0]*centroid_axes[1];
    double xresid, yresid;
    xresid = yresid = 0;
    
    for(int j=0; j<this->axes[0]; j++){
      xresid -= this->pixeldata[xindex+j]*shcentroids.data(j);
      yresid -= this->pixeldata[yindex+j]*shcentroids.data(j);
    }

    znke.set_cos_coeff(1,1,xresid);
    znke.set_sin_coeff(1,1,yresid);
  }

  template<class T> 
  void arroyo_least_squares_reconstructor<T>::
  reconstruct_zonal_residuals(const Arroyo::Shack_Hartmann_centroids & shcentroids, 
			      pixel_array<double> & pixarr) const {

    vector<long> pixarr_axes = pixarr.get_axes();
    if(pixarr_axes[0]!=actuator_axes[0] || 
       pixarr_axes[1]!=actuator_axes[1]){
      cerr << "arroyo_least_squares_reconstructor::reconstruct_zonal_residuals error - " << endl
	   << "pixel array instance passed to reconstructor has axes "
	   << pixarr_axes[0] << "x" << pixarr_axes[1] 
	   << " elements, rather than axes " 
	   << actuator_axes[0] << "x" << actuator_axes[1]
	   << " expected by this reconstructor\n";
      throw(string("arroyo_least_squares_reconstructor::reconstruct_zonal_residuals"));
    }

    long index = 0;
    double resid;
    long nsubaps = centroid_axes[0]*centroid_axes[1]/2;
    int nactuators = actuator_axes[0]*actuator_axes[1];
    int app_sign = app_sign_convention ? -1 : 1;

    for(int i=0; i<actuator_axes[1]; i++){
      for(int j=0; j<actuator_axes[0]; j++){
	resid = 0;
	for(int k=0; k<nsubaps; k++){

	  resid -= app_sign*this->pixeldata[index+k]*shcentroids.data(k);
	  resid -= this->pixeldata[index+nsubaps+k]*shcentroids.data(k+nsubaps);
	}
	pixarr.set_data(j*actuator_axes[1]+i, resid);
	index += this->axes[0];
      }
    }
  }

  template<class T> 
  void arroyo_least_squares_reconstructor<T>::
  reconstruct_residuals(const Arroyo::Shack_Hartmann_centroids & shcentroids, 
			zernike & znke, 
			pixel_array<double> & pixarr) const {
    this->reconstruct_zernike_residuals(shcentroids, znke);
    this->reconstruct_zonal_residuals(shcentroids, pixarr);
  }

  template<class T>
  Arroyo::Shack_Hartmann_centroids 
    arroyo_least_squares_reconstructor<T>::get_slope_consistency(Arroyo::pixel_array<double> & phase_map) const {

    vector<long> lnslt_axes(2,this->centroid_axes[0]);
    Arroyo::Shack_Hartmann_centroids consistency_centroids(lnslt_axes);

    if(!this->geometry_matrix_stored_in_instance){
      cerr << "arroyo_least_squares_reconstructor::get_slope_consistency error - geometry matrix is required "
	   << "to compute slope consistency, but this matrix was not preserved during construction\n";
      throw(std::string("arroyo_least_squares_reconstructor::get_slope_consistency"));
    }

    vector<long> phase_map_axes = phase_map.get_axes();
    if(phase_map_axes[0]!=actuator_axes[0] || 
       phase_map_axes[1]!=actuator_axes[1]){
      cerr << "arroyo_least_squares_reconstructor::get_slope_consistency error - " << endl
	   << "pixel array instance passed to reconstructor has axes "
	   << phase_map_axes[0] << "x" << phase_map_axes[1] 
	   << " elements, rather than axes " 
	   << actuator_axes[0] << "x" << actuator_axes[1]
	   << " expected by this reconstructor\n";
      throw(string("arroyo_least_squares_reconstructor::get_slope_consistency"));
    }

    long nsubaps = centroid_axes[0]*centroid_axes[1]/2;
    long nslopes = 2*nsubaps;
    int nactuators = actuator_axes[0]*actuator_axes[1];
    
    // First version - straight loop, but mixes X and Y and switches sign
    /*
    double tmp;
    for(int k=0; k<nslopes; k++){	
      tmp = 0;
      for(int i=0; i<nactuators; i++){
	tmp += geometry_matrix.data(i*nslopes+k)*phase_map.data(i);
      }
      consistency_centroids.set_data(k, tmp);
    }
    */

    // Switches X and Y, but these are sorted row major instead of column major, and switches sign
    /*
    double tmpx, tmpy;
    for(int k=0; k<nsubaps; k++){	
      tmpx = tmpy = 0;
      for(int i=0; i<nactuators; i++){

	tmpx += geometry_matrix.data(i*nslopes+k+nsubaps)*phase_map.data(i);
	tmpy += geometry_matrix.data(i*nslopes+k)*phase_map.data(i);
      }
      consistency_centroids.set_data(k, tmpx);
      consistency_centroids.set_data(k+nsubaps, tmpy);
    }
    */

    // Version with correct sorting.  With these conventions, one can
    // reconstruct a set of phases from a set of non-circulant
    // centroids using this->reconstruct_residuals() and then recover
    // the input centroids using this->get_slope_consistency().
    // Here, by non-circulant I mean there's no slope discrepancy component
    // to the centroids.
    int app_sign = app_sign_convention ? -1 : 1;
    double tmpx, tmpy;
    for(int k=0; k<centroid_axes[0]; k++){	
      for(int l=0; l<centroid_axes[0]; l++){	
	tmpx = tmpy = 0;
	for(int i=0; i<nactuators; i++){
	  tmpx -= geometry_matrix.data(i*nslopes+k*centroid_axes[0]+l+nsubaps)*phase_map.data(i);
	  tmpy -= app_sign*geometry_matrix.data(i*nslopes+k*centroid_axes[0]+l)*phase_map.data(i);
	}
	consistency_centroids.set_data(l*centroid_axes[0]+k, tmpx);
	consistency_centroids.set_data(l*centroid_axes[0]+k+nsubaps, tmpy);
      }
    }

    return(consistency_centroids);
  };

  template<class T>
  template<class pupil_aperture, class dm_aperture>
    arroyo_least_squares_reconstructor<T>::
    arroyo_least_squares_reconstructor(const pupil_aperture & ap, 
				       const ideal_deformable_mirror<dm_aperture> & ideal_dm,
				       const square_lenslet_array & sq_lnslt_arr,
				       const zernike & projected_modes,
				       bool modal_projection_in_actuator_space,
				       double subaperture_illumination_threshold,
				       double eigenvalue_threshold,
				       bool aplusplus_sign_convention,
				       bool store_eigenmodes,
				       bool store_geometry_matrix){

    this->projected_modes = projected_modes;
    this->actuator_axes = ideal_dm.get_axes();
    this->centroid_axes = sq_lnslt_arr.get_axes();
    this->centroid_axes[1]*=2;
    this->subaperture_illumination_threshold = subaperture_illumination_threshold;
    this->eigenvalue_threshold = eigenvalue_threshold;
    this->eigenmodes_stored_in_instance = store_eigenmodes;
    this->geometry_matrix_stored_in_instance = store_geometry_matrix;
    this->app_sign_convention = aplusplus_sign_convention;

    // enforce requirement that lenslet array is in the plane of the aperture
    if(fabs(sq_lnslt_arr.three_point::z(ap))>three_frame::precision){
      cerr << "arroyo_least_squares_reconstructor::arroyo_least_squares_reconstructor error - \n"
	   << "lenslet array does not lie in the plane of the aperture\n";
      throw(string("arroyo_least_squares_reconstructor::arroyo_least_squares_reconstructor"));      
    }

    // enforce requirement that deformable mirror is in the plane of the aperture
    if(fabs(ideal_dm.three_point::z(ap))>three_frame::precision){
      cerr << "arroyo_least_squares_reconstructor::arroyo_least_squares_reconstructor error - \n"
	   << "deformable mirror does not lie in the plane of the aperture\n";
      throw(string("arroyo_least_squares_reconstructor::arroyo_least_squares_reconstructor"));      
    }

    // enforce requirement that aperture and lenslet array are oriented in the same direction
    if(fabs(1-dot_product(ap.z(), sq_lnslt_arr.z()))>three_frame::precision){
      cerr << "arroyo_least_squares_reconstructor::arroyo_least_squares_reconstructor error - \n"
	   << "cannot construct reconstructor with aperture and lenslet axes misaligned\n";
      throw(string("arroyo_least_squares_reconstructor::arroyo_least_squares_reconstructor"));      
    }

    // enforce (temporary) requirement that aperture and deformable mirror are oriented in the same direction
    if(fabs(1-dot_product(ap.z(), ideal_dm.z()))>three_frame::precision){
      cerr << "arroyo_least_squares_reconstructor::arroyo_least_squares_reconstructor error - \n"
	   << "cannot construct reconstructor with aperture and deformable mirror misaligned\n";
      throw(string("arroyo_least_squares_reconstructor::arroyo_least_squares_reconstructor"));      
    }

    if(projected_modes.get_order()>1){
      cerr << "arroyo_least_squares_reconstructor::arroyo_least_squares_reconstructor error - "
	   << "projections for modes other than tip and tilt are not yet supported\n";
      throw(string("arroyo_least_squares_reconstructor::arroyo_least_squares_reconstructor"));
    }

    // Check for valid geometry threshold
    if(subaperture_illumination_threshold<0 || subaperture_illumination_threshold>1){
      cerr << "arroyo_least_squares_reconstructor::arroyo_least_squares_reconstructor error - "
	   << "invalid subaperture illumination threshold " << subaperture_illumination_threshold << endl;
      throw(string("arroyo_least_squares_reconstructor::arroyo_least_squares_reconstructor"));
    }

    // Check for valid eigenvalue threshold
    if(eigenvalue_threshold<0 || eigenvalue_threshold>1){
      cerr << "arroyo_least_squares_reconstructor::arroyo_least_squares_reconstructor error - "
	   << "invalid eigenvalue threshold " << eigenvalue_threshold << endl;
      throw(string("arroyo_least_squares_reconstructor::arroyo_least_squares_reconstructor"));
    }

    // get lenslet geometry information
    three_vector lenslet_array_origin_offset = sq_lnslt_arr - ap;
    three_vector lenslet_dx = sq_lnslt_arr.get_lenslet_pitch()*sq_lnslt_arr.x();
    three_vector lenslet_dy = sq_lnslt_arr.get_lenslet_pitch()*sq_lnslt_arr.y();
    three_vector lenslet_dx_unit_vector = lenslet_dx*(1/lenslet_dx.length());
    three_vector lenslet_dy_unit_vector = lenslet_dy*(1/lenslet_dy.length());

    double lenslet_buffer = (lenslet_dy+lenslet_dx).length() > (lenslet_dy-lenslet_dx).length() ? 
      .5*(lenslet_dy+lenslet_dx).length() : .5*(lenslet_dy-lenslet_dx).length();

    // get projected actuator spacing
    // SIMPLIFYING ASSUMPTION - FIX THIS LATER
    three_vector dm_origin_offset = ideal_dm - ap;
    three_vector dm_dx = ideal_dm.get_actuator_pitch()*ideal_dm.x();
    three_vector dm_dy = ideal_dm.get_actuator_pitch()*ideal_dm.y();

    double dm_buffer = (dm_dy+dm_dx).length() > (dm_dy-dm_dx).length() ? 
      .5*(dm_dy+dm_dx).length() : .5*(dm_dy-dm_dx).length();

    double geometry_matrix_buffer = 2*(lenslet_buffer+dm_buffer)*(lenslet_buffer+dm_buffer);

    // The lenslet array halfpixel information
    vector<long> lenslet_axes = sq_lnslt_arr.get_axes();
    double lenslet_x_halfpix = lenslet_axes[1]%2==1 ? 0 : .5;
    double lenslet_y_halfpix = lenslet_axes[0]%2==1 ? 0 : .5;
    int lenslet_x_extrapix = lenslet_axes[1]%2==1 ? 1 : 0;
    int lenslet_y_extrapix = lenslet_axes[0]%2==1 ? 1 : 0;

    // The deformable mirror halfpixel information
    vector<long> dm_axes = ideal_dm.get_axes();
    double dm_x_halfpix = dm_axes[1]%2==1 ? 0 : .5;
    double dm_y_halfpix = dm_axes[0]%2==1 ? 0 : .5;
    int dm_x_extrapix = dm_axes[1]%2==1 ? 1 : 0;
    int dm_y_extrapix = dm_axes[0]%2==1 ? 1 : 0;

    long nlenslets = lenslet_axes[0]*lenslet_axes[1];
    long nactuators = dm_axes[0]*dm_axes[1];


    ///////////////////////////////////
    // Construct the geometry matrix //
    ///////////////////////////////////

    vector<long> geometry_matrix_axes(2,2*nlenslets);
    geometry_matrix_axes[1] = nactuators;

    T * geometry_matrix_data;

    try{
      int nelem = geometry_matrix_axes[0]*geometry_matrix_axes[1];
      geometry_matrix_data = new T[nelem];
      for(int i=0; i<nelem; i++)
	geometry_matrix_data[i] = 0;
    } catch(...){
      cerr << "arroyo_least_squares_reconstructor::arroyo_least_squares_reconstructor error - "
	   << "unable to allocate memory for a geometry matrix with " 
	   << geometry_matrix_axes[0]*geometry_matrix_axes[1] << " elements\n";
      throw(string("arroyo_least_squares_reconstructor::arroyo_least_squares_reconstructor"));
    }

    three_point subaperture_center_point;
    three_vector subap_corner_offset_a = .5*(lenslet_dy+lenslet_dx);
    three_vector subap_corner_offset_b = .5*(lenslet_dy-lenslet_dx);
    three_vector actuator_offset_a = dm_dy+dm_dx;
    three_vector actuator_offset_b = dm_dy-dm_dx;

    vector<three_point> subaperture_vertices(4);
    vector<three_point> actuator_influence_function_vertices(4);
    vector<three_point> intersection_vertices;


    long index;
    double overlap;
    double subaperture_area = cross_product(lenslet_dx, lenslet_dy).length();
    double sign_convention = aplusplus_sign_convention ? -1 : 1;

    three_frame local_actuator_frame, local_subap_frame;
    vector<long> unilluminated_subaps;

    double lenslet_pitch = sq_lnslt_arr.get_lenslet_pitch();
    double norm = 1/(lenslet_pitch*lenslet_pitch*lenslet_pitch);

    // These are the geometric relationships between the dm and lenslet coordinates
    double a_0 = dot_product(ideal_dm.x(), sq_lnslt_arr.x());
    double a_1 = dot_product(ideal_dm.x(), sq_lnslt_arr.y());
    double b_0 = dot_product(ideal_dm.y(), sq_lnslt_arr.x());
    double b_1 = dot_product(ideal_dm.y(), sq_lnslt_arr.y());
    double a_2, b_2;

    // These hold the integrals over the polygon
    double S_0, S_x, S_y;

    // Signs for the 4 pyramid faces
    double s_u, s_v;

    try{
      // loop over subaps
      for(int l=-lenslet_axes[1]/2; l<lenslet_axes[1]/2+lenslet_x_extrapix; l++){
	for(int m=-lenslet_axes[0]/2; m<lenslet_axes[0]/2+lenslet_y_extrapix; m++){
	  
	  subaperture_center_point = ap + 
	    (lenslet_array_origin_offset + (l+lenslet_x_halfpix)*lenslet_dx + (m+lenslet_y_halfpix)*lenslet_dy);

	  subaperture_vertices[0] = subaperture_center_point + subap_corner_offset_a;
	  subaperture_vertices[1] = subaperture_center_point + subap_corner_offset_b;
	  subaperture_vertices[2] = subaperture_center_point - subap_corner_offset_a;
	  subaperture_vertices[3] = subaperture_center_point - subap_corner_offset_b;

	  /*
	  subaperture_vertices[0] = subaperture_center_point + subap_corner_offset_a;
	  subaperture_vertices[1] = subaperture_center_point - subap_corner_offset_b;
	  subaperture_vertices[2] = subaperture_center_point - subap_corner_offset_a;
	  subaperture_vertices[3] = subaperture_center_point + subap_corner_offset_b;
	  */

	  overlap = ap.convex_polygon_overlap(subaperture_vertices);

	  if(overlap<subaperture_area*subaperture_illumination_threshold){
	    unilluminated_subaps.push_back((l+lenslet_axes[1]/2)*lenslet_axes[0]+m+lenslet_axes[0]/2);
	    continue;
	  }

	  local_subap_frame = sq_lnslt_arr;
	  local_subap_frame.three_point::operator=(subaperture_vertices[2]);

	  // loop over actuators
	  for(int i=-dm_axes[1]/2; i<dm_axes[1]/2+dm_x_extrapix; i++){
	    for(int j=-dm_axes[0]/2; j<dm_axes[0]/2+dm_y_extrapix; j++){
	      
	      // one pyramid vertex is always coincident with the actuator location
	      actuator_influence_function_vertices[0] = ap + 
		(dm_origin_offset + (i+dm_x_halfpix)*dm_dx + (j+dm_y_halfpix)*dm_dy);
	
	      if((actuator_influence_function_vertices[0]-subaperture_center_point).length_squared()>geometry_matrix_buffer)
		continue;
	      
	      // loop over facets of the pyramid
	      for(int k=0; k<4; k++){
		  
		if(k==0){
		  actuator_influence_function_vertices[1] = actuator_influence_function_vertices[0] + dm_dx;
		  actuator_influence_function_vertices[2] = actuator_influence_function_vertices[1] + dm_dy;
		  actuator_influence_function_vertices[3] = actuator_influence_function_vertices[2] - dm_dx;
		  s_u = s_v = -1;
		} else if(k==1){
		  actuator_influence_function_vertices[1] = actuator_influence_function_vertices[0] + dm_dy;
		  actuator_influence_function_vertices[2] = actuator_influence_function_vertices[1] - dm_dx;
		  actuator_influence_function_vertices[3] = actuator_influence_function_vertices[2] - dm_dy;
		  s_u = 1; s_v = -1;
		} else if(k==2){
		  actuator_influence_function_vertices[1] = actuator_influence_function_vertices[0] - dm_dx;
		  actuator_influence_function_vertices[2] = actuator_influence_function_vertices[1] - dm_dy;
		  actuator_influence_function_vertices[3] = actuator_influence_function_vertices[2] + dm_dx;
		  s_u = s_v = 1;
		} else if(k==3){
		  actuator_influence_function_vertices[1] = actuator_influence_function_vertices[0] - dm_dy;
		  actuator_influence_function_vertices[2] = actuator_influence_function_vertices[1] + dm_dx;
		  actuator_influence_function_vertices[3] = actuator_influence_function_vertices[2] + dm_dy;
		  s_u = -1; s_v = 1;
		}

		intersection_vertices = 
		  get_convex_polygon_intersection(actuator_influence_function_vertices, subaperture_vertices);

		overlap = 0;
		if(intersection_vertices.size()>2){
		  try{overlap = ap.convex_polygon_overlap(intersection_vertices);}
		  catch(...) {
		    cerr << "facet " << k << endl;
		    for(uint i=0; i<actuator_influence_function_vertices.size(); i++)
		      actuator_influence_function_vertices[i].print(cerr, "aif ");
		    cerr << endl;
		    for(uint i=0; i<subaperture_vertices.size(); i++)
		      subaperture_vertices[i].print(cerr, "subap ");
		    cerr << endl;
		    for(uint i=0; i<intersection_vertices.size(); i++)
		      intersection_vertices[i].print(cerr, "int ");
		    cerr << endl;
		    throw(string("overlap"));
		  }
		}

		if(overlap>0){
		  index = ((i+dm_axes[1]/2)*dm_axes[0]+j+dm_axes[0]/2)*2*nlenslets + 
		    (l+lenslet_axes[1]/2)*lenslet_axes[0]+m+lenslet_axes[0]/2;


		  local_actuator_frame = three_frame(ideal_dm);
		  local_actuator_frame.three_point::operator=(actuator_influence_function_vertices[0]);
		  a_2 = local_subap_frame.three_point::x(local_actuator_frame);
		  b_2 = local_subap_frame.three_point::y(local_actuator_frame);

		  S_0 = get_area_of_polygon(intersection_vertices);

		  /*
		  if(l==0 && m==0)
		    convex_polygon_integration(local_subap_frame,
					       intersection_vertices,
					       S_x, 
					       S_y,1);
		  else
		  */
		  convex_polygon_integration(local_subap_frame,
					     intersection_vertices,
					     S_x, 
					     S_y);

		  geometry_matrix_data[index] += 
		    -1*sign_convention*norm*((s_u*a_0*(lenslet_pitch+s_v*b_2)+s_v*b_0*(lenslet_pitch+s_u*a_2))*S_0 +
					     2*s_u*s_v*a_0*b_0*S_x + s_u*s_v*(a_0*b_1+a_1*b_0)*S_y);
		  
		  geometry_matrix_data[index+nlenslets] += 
		    -1*norm*((s_u*a_1*(lenslet_pitch+s_v*b_2)+s_v*b_1*(lenslet_pitch+s_u*a_2))*S_0 +
			     s_u*s_v*(a_0*b_1+a_1*b_0)*S_x + 2*s_u*s_v*a_1*b_1*S_y);
		  /*
		  if(l==0 && m==0){
		    cerr << i << "\t" << j << "\t" << k << endl;
		    
		    cerr << "a0 " << a_0 << "\tb0 " << b_0 << endl;
		    cerr << "a1 " << a_1 << "\tb1 " << b_1 << endl;
		    cerr << "a2 " << a_2/lenslet_pitch << "\tb2 " << b_2/lenslet_pitch << endl;
		    cerr << "su " << s_u << " sv " << s_v << endl;
		    
		    cerr << "S0 " << S_0/lenslet_pitch/lenslet_pitch
			 << " Sx " << S_x/lenslet_pitch/lenslet_pitch/lenslet_pitch
			 << " Sy " << S_y/lenslet_pitch/lenslet_pitch/lenslet_pitch
			 << endl;
		    
		    cerr << "facet "
			 << norm*((s_u*a_0*(lenslet_pitch+s_v*b_2)+s_v*b_0*(lenslet_pitch+s_u*a_2))*S_0 +
				  2*s_u*s_v*a_0*b_0*S_x + s_u*s_v*(a_0*b_1+a_1*b_0)*S_y)
			 << "\t" 
			 << norm*((s_u*a_1*(lenslet_pitch+s_v*b_2)+s_v*b_1*(lenslet_pitch+s_u*a_2))*S_0 +
				  s_u*s_v*(a_0*b_1+a_1*b_0)*S_x + 2*s_u*s_v*a_1*b_1*S_y)
			 << endl;
		    
		    cerr << "geo " << geometry_matrix_data[index] << "\t" 
		    << geometry_matrix_data[index+nlenslets] << endl << endl << endl;
		  }
		  */
		}
	      }
	    }
	  }
	}
      }
    } catch(...){
      cerr << "arroyo_least_squares_reconstructor::arroyo_least_squares_reconstructor error - "
	   << "could not form geometry matrix\n";
      throw(string("arroyo_least_squares_reconstructor::arroyo_least_squares_reconstructor"));
    }

    // Warning - the svd appears to overwrite the input array.  Here we make the copy
    // to preserve the geometry matrix, if this has been requested
    if(store_geometry_matrix)
      geometry_matrix = pixel_array<T>(geometry_matrix_axes, geometry_matrix_data);




    /////////////////////
    // Perform the SVD //
    /////////////////////
    
    char jobu[1];
    char jobvt[1];
    if(store_eigenmodes)
      jobu[0] = jobvt[0] = 'A';
    else 
      jobu[0] = jobvt[0] = 'S';
    int lwork=-1, info;
    T optimal_workspace_size;
    T *singular_values, *g_gtranspose_eigenmode_data, *gtranspose_g_eigenmode_data, *work;

    int geo_mtx_axes_0 = geometry_matrix_axes[0];
    int geo_mtx_axes_1 = geometry_matrix_axes[1];

    singular_value_decomposition<T>(jobu, 
				    jobvt, 
				    geo_mtx_axes_0,
				    geo_mtx_axes_1, 
				    geometry_matrix_data, 
				    geo_mtx_axes_0, 
				    singular_values, 
				    g_gtranspose_eigenmode_data, 
				    geo_mtx_axes_0, 
				    gtranspose_g_eigenmode_data, 
				    geo_mtx_axes_1,
				    &optimal_workspace_size, 
				    lwork, 
				    info);
    
    lwork = (int)optimal_workspace_size;

    try{
      singular_values = new T[geometry_matrix_axes[1]];
      if(store_eigenmodes)
	g_gtranspose_eigenmode_data = new T[geometry_matrix_axes[0]*geometry_matrix_axes[0]];
      else
	g_gtranspose_eigenmode_data = new T[geometry_matrix_axes[0]*geometry_matrix_axes[1]];
      gtranspose_g_eigenmode_data = new T[geometry_matrix_axes[1]*geometry_matrix_axes[1]];
      work = new T[(int)(optimal_workspace_size)];
    } catch(...) {
      cerr << "arroyo_least_squares_reconstructor::arroyo_least_squares_reconstructor error - "
	   << "unable to allocate memory to perform the singular value decomposition\n";
      throw(string("arroyo_least_squares_reconstructor::arroyo_least_squares_reconstructor"));
    }

    singular_value_decomposition<T>(jobu, 
				    jobvt, 
				    geo_mtx_axes_0, 
				    geo_mtx_axes_1, 
				    geometry_matrix_data, 
				    geo_mtx_axes_0, 
				    singular_values, 
				    g_gtranspose_eigenmode_data, 
				    geo_mtx_axes_0, 
				    gtranspose_g_eigenmode_data, 
				    geo_mtx_axes_1,
				    work, 
				    lwork, 
				    info);
    
    delete [] work;
    delete [] geometry_matrix_data;


    // Here we add the number of projected modes to the short axis to
    // allow room in the reconstructor matrix for modal reconstruction
    // from the centroids
    vector<long> reconstructor_axes(geometry_matrix_axes);
    reconstructor_axes[1] += projected_modes.get_axes()[0];
    this->pixel_array<T>::set_axes(reconstructor_axes); 

    // The eigenvalues are sorted in descending order.  We only process
    // entries for which the eigenvalues are above the threshold
    int limit = 0;
    double scaled_singular_value_threshold = 
      sqrt(eigenvalue_threshold*singular_values[0]*singular_values[0]);
    while(limit<geometry_matrix_axes[1] && 
	  singular_values[limit]>scaled_singular_value_threshold) 
      limit++;

    int indexa;
    double tmp;
    // Form the pseudoinverse from V S^{+} U^{T}
    for(int i=0; i<geometry_matrix_axes[1]; i++){
      for(int j=0; j<geometry_matrix_axes[0]; j++){
	tmp = 0;
	indexa = i*geometry_matrix_axes[1];
	for(int k=0; k<limit; k++)
	  tmp += gtranspose_g_eigenmode_data[indexa+k]*
	    g_gtranspose_eigenmode_data[j+k*geometry_matrix_axes[0]]/singular_values[k];

	this->pixeldata[i*geometry_matrix_axes[0]+j]=tmp;
      }
    }

    if(store_eigenmodes)
      sqrt_eigenvalues = pixel_array<T>(vector<long>(1,geometry_matrix_axes[1]), singular_values);
    delete [] singular_values;

    if(store_eigenmodes)
      gtranspose_g_eigenmodes = pixel_array<T>(vector<long>(2,geometry_matrix_axes[1]), gtranspose_g_eigenmode_data);
    delete [] gtranspose_g_eigenmode_data;

    if(store_eigenmodes)
      g_gtranspose_eigenmodes = pixel_array<T>(vector<long>(2,geometry_matrix_axes[0]), g_gtranspose_eigenmode_data);
    delete [] g_gtranspose_eigenmode_data;


    ///////////////////////////////////////////
    // Zero columns in the reconstructor     //
    // corresponding to unilluminated subaps //
    ///////////////////////////////////////////
    /*
    for(int i=0; i<unilluminated_subaps.size(); i++){
      for(int j=0; j<reconstructor_axes[1]; j++){
	this->pixeldata[geometry_matrix_axes[0]*j + unilluminated_subaps[i]] = 0;
	//this->pixeldata[geometry_matrix_axes[1]*unilluminated_subaps[i]+nlenslets] = 0;
      }
    }
    */




    //////////////////////////////////
    // Perform the modal projection //
    //////////////////////////////////

    if(projected_modes.get_order()>0){

      int array_index = 2*nactuators*nlenslets;	
      // Skip piston
      array_index+=2*nlenslets;

      int nunilluminated_subaps = unilluminated_subaps.size();
      int subap_index;

      if(modal_projection_in_actuator_space){
	
	for(int i=0; i<nunilluminated_subaps; i++)
	  unilluminated_subaps.push_back(unilluminated_subaps[i]+nlenslets);

	// Fill in rows of the reconstructor to perform x and y tilt computation
	subap_index = 0;
	for(int i=0; i<geometry_matrix_axes[0]; i++){
	  while(i>unilluminated_subaps[subap_index] && 
		subap_index<nunilluminated_subaps) 
	    subap_index++;
	  if(i!=unilluminated_subaps[subap_index]){
	    for(int j=0; j<nactuators; j++){
	      this->pixeldata[array_index+i] += 
		(2*(j%dm_axes[0])/(double)dm_axes[0] - 1)*this->pixeldata[i+j*geometry_matrix_axes[0]];
	      this->pixeldata[array_index+2*nlenslets+i] += 
		(2*(j/dm_axes[0])/(double)dm_axes[0] - 1)*this->pixeldata[i+j*geometry_matrix_axes[0]]; 
	    }
	  }
	}

	// Project off x and y tilt from the zonal part of the reconstructor
	

	
      } else {
	
	double tmp, mean_fac = 1/(double)(nlenslets - unilluminated_subaps.size());
	
	if(projected_modes.get_cos_coeff(1,1)!=0){	  

	  // Project the x tilt off the zonal elements
	  for(int i=0; i<geometry_matrix_axes[1]; i++){
	    subap_index = 0;
	    tmp = 0;
	    for(int j=0; j<nlenslets; j++)
	      tmp += this->pixeldata[i*geometry_matrix_axes[0]+j];

	    tmp *= mean_fac;

	    for(int j=0; j<nlenslets; j++){
	      while(j>unilluminated_subaps[subap_index] && 
		    subap_index<nunilluminated_subaps) 
		subap_index++;
	      if(j!=unilluminated_subaps[subap_index])
		this->pixeldata[i*geometry_matrix_axes[0]+j] -= tmp;

	    }
	  }

	  // Fill in a row of the reconstructor to perform x tilt computation
	  subap_index = 0;
	  for(int j=0; j<nlenslets; j++){
	    while(j>unilluminated_subaps[subap_index] && 
		  subap_index<nunilluminated_subaps) 
	      subap_index++;
	    if(j!=unilluminated_subaps[subap_index])
	      this->pixeldata[array_index+j] = 2*mean_fac;
	  }
	  array_index += 2*nlenslets;
	}

	if(projected_modes.get_sin_coeff(1,1)!=0){	  

	  // Project the y tilt off the zonal elements
	  for(int i=0; i<geometry_matrix_axes[1]; i++){
	    subap_index = 0;
	    tmp = 0;
	    for(int j=nlenslets; j<geometry_matrix_axes[0]; j++)
	      tmp += this->pixeldata[i*geometry_matrix_axes[0]+j];

	    tmp *= mean_fac;

	    for(int j=nlenslets; j<geometry_matrix_axes[0]; j++){
	      while(j>(unilluminated_subaps[subap_index]+nlenslets) && 
		    subap_index<nunilluminated_subaps) 
		subap_index++;
	      if(j!=(unilluminated_subaps[subap_index]+nlenslets))
		this->pixeldata[i*geometry_matrix_axes[0]+j] -= tmp;

	    }
	  }

	  // Fill in a row of the reconstructor to perform y tilt computation
	  subap_index = 0;
	  for(int j=nlenslets; j<geometry_matrix_axes[0]; j++){
	    while(j>(unilluminated_subaps[subap_index]+nlenslets) && 
		  subap_index<nunilluminated_subaps) 
	      subap_index++;
	    if(j!=(unilluminated_subaps[subap_index]+nlenslets))
	      this->pixeldata[array_index+j] = 2*mean_fac;
	  }
	  array_index += 2*nlenslets;
	}
      }
    }
  }
}

#endif
