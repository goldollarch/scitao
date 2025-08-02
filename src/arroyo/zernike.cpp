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

#include <math.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include "sim_utils.h"
#include "zernike.h"
#include "aperture.h"
#include "pixel_array.h"
#include "fits_factory.h"

using namespace std;

namespace Arroyo {

  int zernike::verbose_level = 0;

  namespace factory_register {
    const fits_keyval_set & get_zernike_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "zernike expansion"));
      return *fkvs;
    }
    
    AO_sim_base * create_zernike(const iofits & iof) {
      return new zernike(iof);
    }

  } 

  const bool zernike::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_zernike_keyval_set(), 
				      factory_register::create_zernike);
  

  namespace {

    double factorial_combination(long order, long level, long s){
      double val = 1;
      
      for(long i=2; i<=order; i++){
	if(i<=order-s) val *= i;
	if(i<=s) val /= (double)i;
	if(i<=((order+level)/2-s)) val /= (double)i;
	if(i<=((order-level)/2-s)) val /= (double)i;
      }
      return(val);
    }

    // This function constructs a radial zernike polynomial on a 
    // grid according to the generating function 
    // R^level_order = SUM^(order-level)/2_0 
    //           (-1)^i (order - i)!/(((order+level)/2-i)! ((order-level)/2-i)! rho^(order-2i)
    // See Born and Wolf pp 464-466
    // For convenience, (order+level)/2-i)! is defined as order_plus_level_over_two_minus_index
    //                  (order-level)/2-i)! is defined as order_minus_level_over_two_minus_index
    double radial_zernike_polynomial(long order, long level, double rho){
    
      if(order<0 || level<0 || order < level || order%2 != level%2 || (order-level)%2!=0){
	cerr << "zernike::radial_zernike_polynomial error - invalid zernike indices: "
	     << "order " << order << " level " << level << endl;
	throw(string("zernike::radial_zernike_polynomial"));
      }
      
      if(rho<0 || rho>1){
	cerr << "zernike::radial_zernike_polynomial error - invalid value "
	     << rho << " for the radial coordinate, which must lie between zero and one\n";
	throw(string("zernike::radial_zernike_polynomial"));
      }
      
      // Factors that appear in the zernike generating function
      int minus_one = -1;
      double val = 0;
      for(int s=0; s<=(order-level)/2; s++){
	minus_one *= -1;
	val += minus_one*factorial_combination(order, level, s)*pow(rho, (double)(order-2*s));
      }
      // Normalization factor - with this, the radial zernike
      // polynomials are orthonormal See Born & Wolf section
      // 9.2.1 eqn 2
      val *= sqrt((order+1)/M_PI);
      return(val);
    }

    /*
    void inner_product(const double * const fn_array, 
		       long dimen, 
		       long aperture_diameter_pixels, 
		       long order, 
		       long level,
		       double & cos_coeff, 
		       double & sin_coeff){

      if(zernike::verbose_level)
	cout << "zernike::inner_product:  dimen " << dimen
	     << " aperture diameter in pixels " << aperture_diameter_pixels 
	     << " order " << order
	     << " level " << level
	     << endl;

      if(dimen <=0){
	cerr << "inner_product error - dimension " << dimen << " out of range\n";
	throw(string("inner_product"));
      }

      if(dimen < aperture_diameter_pixels){
	cerr << "inner_product error - dimension " << dimen 
	     << " less than aperture diameter " << aperture_diameter_pixels << endl;
	throw(string("inner_product"));
      }

      double halfpix=0;
      int extrapix=1;
      if(dimen%2==0){
	halfpix = .5;
	extrapix = 0;
      }

      cos_coeff = 0;
      sin_coeff = 0;
      double xcoord, ycoord, rho, phi;
      double delta = 2.0/aperture_diameter_pixels;
      double delta_sq = delta*delta;
      long index;

      // This factor normalizes the angular part of the polynomial
      // There's a factor of 1/2 if level!=0 for the averaging
      // of the trig terms (cos^{2} or sin^{2})
      double fac = 1;
      if(level!=0) fac = sqrt(2.0);
      double radial_val;

      for(int k=-dimen/2; k<dimen/2 + extrapix; k++){
	for(int l=-dimen/2; l<dimen/2 + extrapix; l++){
	  xcoord = (k+halfpix)*delta;
	  ycoord = (l+halfpix)*delta;
	  rho = xcoord*xcoord + ycoord*ycoord;
	  if(rho<1){
	    rho = sqrt(rho);
	    radial_val = radial_zernike_polynomial(order, level, rho);
	    index = (k+dimen/2)*dimen+l+dimen/2;
	    phi = atan2(ycoord,xcoord);
	    cos_coeff += 
	      fac*radial_val*cos(level*phi)*fn_array[index]*delta_sq;
	    sin_coeff += 
	      fac*radial_val*sin(level*phi)*fn_array[index]*delta_sq;

	    if(zernike::verbose_level)
	      cout << k << "\t" << l 
		   << "\t" << fac 
		   << "\t" << radial_val 
		   << "\t" << cos(level*phi) 
		   << "\t" << sin(level*phi) 
		   << "\t" << delta_sq 
		   << "\t" << fn_array[index] 
		   << endl;

	  }
	}
      }

      if(fabs(cos_coeff)<1e-15) cos_coeff = 0;
      if(fabs(sin_coeff)<1e-15) sin_coeff = 0;
    }
  }
    */

    void inner_product(pixel_array<double> & pixarr,
		       long aperture_diameter_pixels, 
		       long order, 
		       long level,
		       double & cos_coeff, 
		       double & sin_coeff){


      vector<long> axes = pixarr.get_axes();
      if(axes.size()!=2 || axes[0]!=axes[1]){
	cerr << "zernike::inner_product error - axes mismatch\n";
	cerr << "\tsize " << axes.size() << endl;
	for(int i=0; i<axes.size(); i++)
	  cerr << "\tvalue " << axes[i] << endl;
      }

      int dimen = axes[0];

      if(zernike::verbose_level)
	cout << "zernike::inner_product:  dimen " << dimen
	     << " aperture diameter in pixels " << aperture_diameter_pixels 
	     << " order " << order
	     << " level " << level
	     << endl;

      if(dimen <=0){
	cerr << "inner_product error - dimension " << dimen << " out of range\n";
	throw(string("inner_product"));
      }

      if(dimen < aperture_diameter_pixels){
	cerr << "inner_product error - dimension " << dimen 
	     << " less than aperture diameter " << aperture_diameter_pixels << endl;
	throw(string("inner_product"));
      }

      double halfpix=0;
      int extrapix=1;
      if(dimen%2==0){
	halfpix = .5;
	extrapix = 0;
      }

      cos_coeff = 0;
      sin_coeff = 0;
      double xcoord, ycoord, rho, phi;
      double delta = 2.0/aperture_diameter_pixels;
      double delta_sq = delta*delta;
      long index;

      // This factor normalizes the angular part of the polynomial
      // There's a factor of 1/2 if level!=0 for the averaging
      // of the trig terms (cos^{2} or sin^{2})
      double fac = 1;
      if(level!=0) fac = sqrt(2.0);
      double radial_val;

      for(int k=-dimen/2; k<dimen/2 + extrapix; k++){
	for(int l=-dimen/2; l<dimen/2 + extrapix; l++){
	  xcoord = (k+halfpix)*delta;
	  ycoord = (l+halfpix)*delta;
	  rho = xcoord*xcoord + ycoord*ycoord;
	  if(rho<1){
	    rho = sqrt(rho);
	    radial_val = radial_zernike_polynomial(order, level, rho);
	    index = (k+dimen/2)*dimen+l+dimen/2;
	    phi = atan2(ycoord,xcoord);
	    cos_coeff += 
	      fac*radial_val*cos(level*phi)*pixarr.data(index)*delta_sq;
	    sin_coeff += 
	      fac*radial_val*sin(level*phi)*pixarr.data(index)*delta_sq;

	    if(zernike::verbose_level)
	      cout << k << "\t" << l 
		   << "\t" << fac 
		   << "\t" << radial_val 
		   << "\t" << cos(level*phi) 
		   << "\t" << sin(level*phi) 
		   << "\t" << delta_sq 
		   << "\t" << pixarr.data(index)
		   << endl;
	  }
	}
      }
      if(fabs(cos_coeff)<1e-15) cos_coeff = 0;
      if(fabs(sin_coeff)<1e-15) sin_coeff = 0;
    }
  }


  zernike::zernike(){
    this->set_axes(vector<long>(1,0));
  }

  zernike::zernike(const Arroyo::pixel_array<double> & pixarr, 
		   double pixscale,
		   long order,
		   const circular_aperture & circ_ap){

    try{
      if(order<0){
	cerr << "zernike::zernike error - order " 
	     << order << " less than zero\n";
	throw(string("zernike::zernike"));
      }

      vector<long> pixarr_axes = pixarr.get_axes();

      if(pixarr_axes.size()!=2){
	cerr << "zernike::zernike error - "
	     << "pixel_array does not have 2 axes\n";
	throw(string("zernike::zernike"));
      }

      if(pixarr_axes[0]!=pixarr_axes[1]){
	cerr << "zernike::zernike error - "
	     << "for the moment this function is restriced to square arrays\n" 
	     << " but the pixel array supplied to this function has dimensions "
	     << pixarr_axes[0] << " x " << pixarr_axes[1] << endl;
	throw(string("zernike::zernike"));
      }

      if(circ_ap.get_diameter()/pixscale > pixarr_axes[0]){
	cerr << "zernike::zernike error - "
	     << "pixel array is smaller than the circular aperture\n";
	cerr << "aperture diameter " << circ_ap.get_diameter() << endl;
	cerr << "pixel scale " << pixscale << endl;
	cerr << "number of pixels across aperture " << circ_ap.get_diameter()/pixscale << endl;
	cerr << "pixel array axes " << pixarr_axes[0] << endl;
	throw(string("zernike::zernike"));
      }

      this->resize(order);

      Arroyo::pixel_array<double> tmp_pixarr(pixarr), tmp_modal_pixarr;

      double cos_coeff, sin_coeff;
      for(int i=0; i<=order; i++){
	for(int j=i%2; j<=i; j+=2){
	  inner_product(tmp_pixarr,
			(long)(ceil(circ_ap.get_diameter()/pixscale)),
			i, j, 
			cos_coeff, sin_coeff);
	  this->set_cos_coeff(i,j,cos_coeff);
	  if(j!=0)
	    this->set_sin_coeff(i,j,sin_coeff);
	}

	// If I don't subtract out the phase surfaces I have fitted
	// order by order, there is a lot of cross-contamination.
	// The big problem here is that atmospheric turbulence is 
	// dominated by low spatial frequency modes, so that even a 
	// little cross-contamination can destroy the fit for high
	// spatial frequency modes.

	zernike tmp_zernike(i);
	for(int j=i%2; j<=i; j+=2){
	  tmp_zernike.set_cos_coeff(i, j, this->get_cos_coeff(i,j));
	  if(j!=0)
	    tmp_zernike.set_sin_coeff(i, j, this->get_sin_coeff(i,j));
	}

	if(zernike::verbose_level)
	  tmp_zernike.print(cout, "tmp znke ");

	tmp_modal_pixarr = tmp_zernike.get_pixel_array(pixarr_axes, 
						       pixscale, 
						       circ_ap);
	tmp_pixarr -= tmp_modal_pixarr;
      }
    } catch(...) {
      cerr << "zernike::zernike error - could not expand pixel phase array\n";
      throw(string("zernike::zernike"));
    }
  }

  pixel_array<double> zernike::get_pixel_array(const vector<long> & axes, 
					       double pixscale,
					       const aperture & aper) const {
    try{
      dynamic_cast<const circular_aperture &>(aper);
    } catch(...) {
      cerr << "zernike::get_pixel_array error - \n"
	   << "cannot construct pixel array on an aperture that is not circular\n";
      aper.print(cerr, "invalid aperture ");
      throw(string("zernike::get_pixel_array"));
    }
    return(this->get_pixel_array(axes, pixscale, 
				 dynamic_cast<const circular_aperture &>(aper)));
  }
    

  pixel_array<double> zernike::get_pixel_array(const vector<long> & axes,
					       double pixscale,
					       const circular_aperture & circ_ap) const {


    if(axes.size()!=2){
      cerr << "zernike::get_pixel_array error - "
	   << "axes have dimension " << axes.size() 
	   << " rather than dimension 2\n";
      throw(string("zernike::get_pixel_array"));
    }

    if(axes[0]!=axes[1]){
      cerr << "zernike::get_pixel_array error - "
	   << "for the moment this function is restriced to square arrays,\n" 
	   << "but the axes supplied to this function have dimensions " 
	   << axes[0] << " x " << axes[1] << endl;
      throw(string("zernike::get_pixel_array"));
    }

    if(ceil(circ_ap.get_diameter()/pixscale) > axes[0]){
      cerr << "zernike::get_pixel_array error - "
	   << "pixel phase array is smaller than the circular aperture\n";
      cerr << "aperture diameter " << circ_ap.get_diameter() << endl;
      cerr << "pixel scale " << pixscale << endl;
      cerr << "number of pixels across aperture " << circ_ap.get_diameter()/pixscale << endl;
      cerr << "axes " << axes[0] << endl;
      throw(string("zernike::get_pixel_array"));
    }

    pixel_array<double> pixarr(axes);
    
    double xcoord, ycoord, rho, phi;
    double delta = 2.0/ceil(circ_ap.get_diameter()/pixscale);

    double halfpix=0;
    int extrapix=1;
    if(axes[0]%2==0){
      halfpix = .5;
      extrapix = 0;
    }

    long index;
    long order = this->get_order();
    double radial_val;
    bool areal_weighting = circ_ap.get_areal_weighting();
    double val, valsum;
    try{ 
      for(int i=0; i<=order; i++){  
	for(int j=i%2; j<=i; j+=2){
	  if(this->get_cos_coeff(i,j)!=0 || (j!=0 && this->get_sin_coeff(i,j)!=0)){
	    // This factor normalizes the angular part of the polynomial
	    // There's a factor of 1/2 if level!=0 for the averaging
	    // of the trig terms (cos^{2} or sin^{2})
	    double fac = 1;
	    if(j!=0) fac = sqrt(2.0);
	    
	    valsum=0;

	    for(int k=-axes[1]/2; k<axes[1]/2+extrapix; k++){
	      for(int l=-axes[0]/2; l<axes[0]/2+extrapix; l++){
		xcoord = (k+halfpix)*delta;
		ycoord = (l+halfpix)*delta;
		rho = sqrt(xcoord*xcoord+ycoord*ycoord);
		index = (k+axes[1]/2)*axes[0] + l + axes[0]/2;
		if(rho<1){
		  radial_val = radial_zernike_polynomial(i, j, rho);
		  phi = atan2(ycoord,xcoord);
		  if(j!=0)
		    val = fac*radial_val*
		      (cos(j*phi)*this->get_cos_coeff(i,j) + 
		       sin(j*phi)*this->get_sin_coeff(i,j));
		  else
		    val = fac*radial_val*
		      cos(j*phi)*this->get_cos_coeff(i,j);
		    
		} else val = 0;

		pixarr.set_data(index, pixarr.data(index)+ val);

		/*		
		else if(areal_weighting && 
			  (((fabs(xcoord)-.5)*(fabs(xcoord)-.5)+
			  (fabs(ycoord)-.5)*(fabs(ycoord)-.5)) < 1)){
		  radial_val = radial_zernike_polynomial(i, j, 1);
		  index = (k+axes[1]/2)*axes[0] + l + axes[0]/2;
		  phi = atan2(ycoord,xcoord);
		  tmp_array[index] += fac*radial_val*
		    (cos(j*phi)*cos_zernike_coeffs[i][j/2] + 
		     sin(j*phi)*sin_zernike_coeffs[i][j/2]);
		}
		*/
	      }
	    }
	  } 
	}
      } 
    } catch(...) {
      cerr << "zernike::get_pixel_array - error constructing array\n";
      throw(string("zernike::get_pixel_array"));
    }
    return(pixarr);
  }

  zernike & zernike::operator=(const zernike & znke){
    if(this==&znke)
      return(*this);

    //this->pixel_array<double>::verbose_level = 1;
    this->pixel_array<double>::operator=(znke);

    return(*this);
  }
  
  void zernike::resize(long order){

    if(order<0){
      cerr << "zernike::resize error - cannot resize to order "
	   << order << endl;
      throw(string("zernike::resize"));
    }
    
    if(this->get_order()==order) return;
    int old_nelem = this->total_space();
    int new_nelem = 0;
    for(int i=0; i<=order; i++){
      if(i%2) new_nelem += i/2+i/2+2; 
      else new_nelem += i/2+i/2+1;
    }

    zernike tmp_znke(*this);

    vector<long> new_axes(1,new_nelem);
    this->pixel_array<double>::set_axes(new_axes);
    for(int i=0; i<new_nelem && i<old_nelem; i++)
      this->pixeldata[i] = tmp_znke.pixel_array<double>::data(i);
  }

  void zernike::read(const char * filename) {
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "zernike::read - "
	   << "error opening file " << filename << endl;
      throw(string("zernike::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "zernike::read - "
	   << "error reading "
	   << this->unique_name() << " from file " 
	   << filename << endl;
      throw(string("zernike::read"));
    }
  }

  void zernike::read(const iofits & iof) {
    if(!iof.key_exists("TYPE")){
      cerr << "zernike::read error - "
	   << "unrecognized type of file\n";
      throw(string("zernike::read"));
    }
    string type, comment;
    iof.read_key("TYPE", type, comment);
    if(type!=this->unique_name()){
      cerr << "zernike::read error - file of type " 
	   << type << " rather than of type "
	   << this->unique_name() << endl;
      throw(string("zernike::read"));
    }

    this->pixel_array<double>::read(iof);

    // leave iof pointing to the next header, if it exists
    if(iof.get_hdu_num()!=iof.get_num_hdus())
      iof.movrel_hdu(1);
  }

  void zernike::write(const char * filename) const {
   iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "zernike::write - "
	   << "error opening file " << filename << endl;
      throw(string("zernike::write"));
    }

    try{this->write(iof);}
    catch(...){
      cerr << "zernike::write - "
	   << "error writing "
	   << this->unique_name() << " to file " 
	   << filename << endl;
      throw(string("zernike::write"));
    }
  }

  void zernike::write(iofits & iof) const {
    fits_header_data<double> fhd(axes);
    fhd.write(iof);
    string type = this->unique_name();
    string comment = "object type";
    iof.write_key("TYPE", type, comment);
    iof.write_key("MAXORDER", this->get_order(), string("maximum order in zernike expansion"));
    this->pixel_array<double>::write(iof);
  }

  namespace {

    bool check_cos_order_and_level(int order, int level){
      if(order<0) return(false);
      if(level<0 || level>order || (level%2!=order%2)) return(false);
      return(true);
    }
    
    bool check_sin_order_and_level(int order, int level){
      if(level==0) return(false);
      if(!check_cos_order_and_level(order, level)) return(false);
      return(true);
    }

    int get_cos_index(int order, int level){
      if(!check_cos_order_and_level(order, level)){
	cerr << "get_cos_index - " 
	     << " error for order " << order << " level " << level << endl;
	throw(string("get_cos_index"));
      }
      int index = 0;
      for(int i=0; i<order; i++)
	index += i+1;
      
      int tmp_level = order;
      while(tmp_level > level){
	if(tmp_level==0) index+=1;
	else index+=2;
	tmp_level-=2;
      }
      return(index);
    }

    int get_sin_index(int order, int level){
      if(!check_cos_order_and_level(order, level)){
	cerr << "get_sin_index - " 
	     << " error for order " << order << " level " << level << endl;
	throw(string("get_sin_index"));
      }
      int index = get_cos_index(order, level)+1;
      return(index);
    }

  }
  
  double zernike::get_cos_coeff(int order, 
				int level) const {
    return(pixeldata[get_cos_index(order, level)]);
  }

  void zernike::set_cos_coeff(int order, 
			      int level, 
			      double coeff){
    if(this->get_order()<order)
      this->resize(order);
    pixeldata[get_cos_index(order, level)] = coeff;
  }
  
  double zernike::get_sin_coeff(int order, 
				int level) const {
    return(pixeldata[get_sin_index(order, level)]);
  }

  void zernike::set_sin_coeff(int order, 
			      int level, 
			      double coeff){

    if(this->get_order()<order)
      this->resize(order);
    pixeldata[get_sin_index(order, level)] = coeff;
  }

  void zernike::print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield);
    fits_header_data<double> fhd(axes);
    os << prefix << "TYPE       = " << setw(vlspc) << this->unique_name()
       << "/" << "object type" << endl;
    fhd.print(os, prefix);
    os << prefix << "MAXORDER   = " << setw(vlspc) << this->get_order()
       << "/" << "maximum order in expansion" << endl;
    for(int i=0; i<=this->get_order(); i++){
      for(int j=i; j>=i%2; j-=2){
	os << prefix << "COEFF " << i << ", " << j << "\t" << this->get_cos_coeff(i,j);
	if(j==0) os << endl;
	else os << "\t" << this->get_sin_coeff(i,j) << endl;
      }
    }
  }

  long zernike::get_order() const {
    int nelem = this->total_space();
    if(nelem==0) return(-1);

    int order = -1;
    while(nelem>0){
      order++;
      if(order%2) nelem -= order/2+order/2+2; 
      else nelem -= order/2+order/2+1;
    } 
    if(nelem!=0){
      cerr << "zernike::get_order error - could not compute order\n";
      throw(string("zernike::get_order"));
    }

    return(order);
  }
}



