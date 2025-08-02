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

#include <cmath>
#include <iomanip>
#include "Cn2_distribution.h"

using namespace std;

namespace Arroyo {

  Hardy_Cn2_distribution::Hardy_Cn2_distribution() {
    boundary_coeff_ = boundary_scale_ = 0;
    troposphere_coeff_ = troposphere_scale_ = 0;
    tropopause_coeff_ = tropopause_height_ = 0;
  }

  Hardy_Cn2_distribution::Hardy_Cn2_distribution(const Hardy_Cn2_distribution & tcn2d) {
    boundary_coeff_ = tcn2d.boundary_coeff_;
    boundary_scale_ = tcn2d.boundary_scale_;   
    troposphere_coeff_ = tcn2d.troposphere_coeff_;
    troposphere_scale_ = tcn2d.troposphere_scale_;
    tropopause_coeff_ = tcn2d.tropopause_coeff_;
    tropopause_height_ = tcn2d.tropopause_height_;
    additional_layers_ = tcn2d.additional_layers_;
  }

  Hardy_Cn2_distribution::Hardy_Cn2_distribution(double boundary_coeff, double boundary_scale,
						 double troposphere_coeff, double troposphere_scale,
						 double tropopause_coeff, double tropopause_height, 
						 vector<vector<double> > additional_layers) {

    if(boundary_coeff_<0 || troposphere_coeff<0 || tropopause_coeff<0){
      cerr << "Hardy_Cn2_distribution::Hardy_Cn2_distribution error - "
	   << "coefficients supplied to constructor are not all positive\n";
      throw(string("Hardy_Cn2_distribution::Hardy_Cn2_distribution"));
    }
    if((boundary_scale<=0 && boundary_coeff>0) ||
       (troposphere_scale<=0 && troposphere_coeff>0) ||
       (tropopause_height<=0 && tropopause_coeff>0)){
      cerr << "Hardy_Cn2_distribution::Hardy_Cn2_distribution error - "
	   << "bad coeff/scale combination\n";
      throw(string("Hardy_Cn2_distribution::Hardy_Cn2_distribution"));
    }

    boundary_coeff_ = boundary_coeff;
    boundary_scale_ = boundary_scale;
    troposphere_coeff_ = troposphere_coeff;
    troposphere_scale_ = troposphere_scale;
    tropopause_coeff_ = tropopause_coeff;
    tropopause_height_ = tropopause_height;

    for(int i=0; i<additional_layers.size(); i++){
      if(additional_layers[i].size()!=3){
	cerr << "Hardy_Cn2_distribution::Hardy_Cn2_distribution error - "
	     << "each additional layer must have 3 entries, but " 
	     << additional_layers[i].size() << " were provided\n";
	throw(string("Hardy_Cn2_distribution::Hardy_Cn2_distribution"));
      }
      for(int j=0; j<3; j++){
	if(additional_layers[i][j]<0){
	  cerr << "Hardy_Cn2_distribution::Hardy_Cn2_distribution error - "
	       << "each additional layer must have positive entries, but the value " 
	       << additional_layers[i][j] << " was provided for the ";
	  if(j==0) cerr << "coefficient of layer " << i << endl;
	  if(j==1) cerr << "height of layer " << i << endl;
	  if(j==2) cerr << "thickness of layer " << i << endl;
	  throw(string("Hardy_Cn2_distribution::Hardy_Cn2_distribution"));
	}
      }
    }
    additional_layers_ = additional_layers;
  }

  void Hardy_Cn2_distribution::add_layer(double coeff, double height, double thickness){
    if(coeff < 0 || height < 0 || thickness < 0){
      cerr << "Hardy_Cn2_distribution::add_layer error - " << endl;
      if(coeff<0) cerr << "\tcoefficient " << coeff << " less than zero\n";
      if(height<0) cerr << "\theight " << height << " less than zero\n";
      if(thickness<0) cerr << "\tthickness " << thickness << " less than zero\n";
      throw(string("Hardy_Cn2_distribution::add_layer"));
    }
    vector<double> v(3,coeff);
    v[1] = height;
    v[2] = thickness;
    additional_layers_.push_back(v);
  }

  double Hardy_Cn2_distribution::val(double height) const {

    if(height<0){
      cerr << "Hardy_Cn2_distribution::val error - "
	   << "height " << height << " specified in argument "
	   << "is not positive\n";
      throw(string("Hardy_Cn2_distribution::val"));
    }

    double val = 0;

    if(boundary_coeff_>0)
      val += boundary_coeff_*exp(-height/boundary_scale_);
    if(troposphere_coeff_>0)
      val += troposphere_coeff_*exp(-height/troposphere_scale_);
    if(tropopause_coeff_>0)
      val += tropopause_coeff_*pow(height,10.0)*exp(-height/tropopause_height_);

    double expntl;
    for(int i=0; i<additional_layers_.size() && additional_layers_[i][0]>0; i++){
      expntl = (height-additional_layers_[i][1])*(height-additional_layers_[i][1])/
	(2*additional_layers_[i][2]*additional_layers_[i][2]);
      val += additional_layers_[i][0]*exp(-expntl);
    }
    return(val);
  }

  double Hardy_Cn2_distribution::get_moment(double exponent) const {

    // Here we make sure we integrate up to at least to the highest layer
    double max_height = boundary_scale_;
    if(max_height<troposphere_scale_) max_height = troposphere_scale_;
    if(max_height<tropopause_height_) max_height = tropopause_height_;
    for(int i=0; i<additional_layers_.size(); i++)
      if(max_height<additional_layers_[i][1]) max_height = additional_layers_[i][1];

    double height=0, step=2; // in meters
    double moment=0, step_value;

    do {
      step_value = this->val(height)*step;
      if(exponent>0)
	step_value *= pow(height,exponent);
      moment += step_value;
      height += step;
    } while(step_value>moment*1e-8 || height<max_height);
    // above, terminate when additional steps contribute less than the
    // precision of a float to the value of the moment.
    return(moment);
  }

  void Hardy_Cn2_distribution::print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "       TYPE = " << setw(vlspc) << "Hardy Cn2 distribution"
       << "/" << "object type" << endl;
    os << prefix << "   BDRYCOEF = " << setw(vlspc) << boundary_coeff_
       << "/" << "boundary coefficient (m^{2/3})" << endl;
    os << prefix << "   BDRYSCLE = " << setw(vlspc) << boundary_scale_
       << "/" << "boundary scale height (m)" << endl;
    os << prefix << "   TPSPCOEF = " << setw(vlspc) << troposphere_coeff_
       << "/" << "troposphere coefficient (m^{2/3})" << endl;
    os << prefix << "   TPSPSCLE = " << setw(vlspc) << troposphere_scale_
       << "/" << "troposphere scale height (m)" << endl;
    os << prefix << "   TPPSCOEF = " << setw(vlspc) << tropopause_coeff_
       << "/" << "tropopause coefficient (m^{2/3})" << endl;
    os << prefix << "   TPPSSCLE = " << setw(vlspc) << tropopause_height_
       << "/" << "tropopause height (m)" << endl;
    for(int i=0; i<additional_layers_.size(); i++){
      os << prefix << "      COEFF = " << setw(vlspc) << additional_layers_[i][0]
	 << "/" << "layer " << i+1 << " coefficient (m^{2/3})" << endl;
      os << prefix << "     HEIGHT = " << setw(vlspc) << additional_layers_[i][1]
	 << "/" <<  "layer " << i+1 << " height (m)" << endl;
      os << prefix << "   THCKNESS = " << setw(vlspc) << additional_layers_[i][2]
	 << "/" << "layer " << i+1 << " thickness (m)" << endl;
    }       
  }

  Cn2_distribution * Hardy_Cn2_distribution::clone() const {
    return(new Hardy_Cn2_distribution(*this));
  }

  bool Hardy_Cn2_distribution::equals(Cn2_distribution * cn2dist) const {
    Hardy_Cn2_distribution * tmp;
    if((tmp=dynamic_cast<Hardy_Cn2_distribution *>(cn2dist)) && operator==(*tmp, *this))
      return(true);
    return(false);
  }

  bool operator==(const Hardy_Cn2_distribution & t1, const Hardy_Cn2_distribution & t2){
    if(t1.boundary_coeff_ == t2.boundary_coeff_ && 
       t1.boundary_scale_ == t2.boundary_scale_  && 
       t1.troposphere_coeff_ == t2.troposphere_coeff_ &&
       t1.troposphere_scale_ == t2.troposphere_scale_ &&
       t1.tropopause_coeff_ == t2.tropopause_coeff_ &&
       t1.tropopause_height_ == t2.tropopause_height_ &&
       t1.additional_layers_ == t2.additional_layers_)
      return(true);
    return(false);
  }

  double SLCSAT_day_Cn2_distribution::val(double height) const {
    if(height<0){
      cerr << "SLCSAT_day_Cn2_distribution::val error - "
	   << "height " << height << " supplied as argument is negative\n";
      throw(string("SLCSAT_day_Cn2_distribution::val"));
    }

    if(height < 18.5) return(0);
    else if(height < 232) return(3.96e-13/pow(height, 1.05));
    else if(height < 880) return(1.3e-15);
    else if(height < 7220) return(8.87e-7/height/height/height);
    else if(height < 20500) return(2e-16/sqrt(height));
    else return(0);
  }

  double SLCSAT_day_Cn2_distribution::get_moment(double exponent) const {
    double step=1; // in meters
    double moment=0, step_value;
    for(double height=18.5; height<20500; height+=step)
      moment += this->val(height)*pow(height,exponent)*step;
    return(moment);
  }

  void SLCSAT_day_Cn2_distribution::print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "       TYPE = " << setw(vlspc) << "SLCSAT day Cn2 distribution"
       << "/" << "object type" << endl;
  }

  Cn2_distribution * SLCSAT_day_Cn2_distribution::clone() const {
    return(new SLCSAT_day_Cn2_distribution());
  }

  bool SLCSAT_day_Cn2_distribution::equals(Cn2_distribution * cn2dist) const {
    if(dynamic_cast<SLCSAT_day_Cn2_distribution *>(cn2dist))
      return(true);
    return(false);
  }

  double SLCSAT_night_Cn2_distribution::val(double height) const {
    if(height<0){
      cerr << "SLCSAT_night_Cn2_distribution::val error - "
	   << "height " << height << " supplied as argument is negative\n";
      throw(string("SLCSAT_night_Cn2_distribution::val"));
    }

    if(height<18.5) return(5e-15);
    else if(height<110) return(2.875e-12/height/height);
    else if(height<1500) return(2.5e-16);
    else if(height<7000) return(8.87e-7/height/height/height);
    else if(height<20500) return(2e-16/sqrt(height));
    else return(0);
  }

  double SLCSAT_night_Cn2_distribution::get_moment(double exponent) const {
    double step = 1; // in meters
    double moment=0, step_value;
    for(double height=18.5; height<20500; height+=step)
      moment += this->val(height)*pow(height,exponent)*step;
    return(moment);
  }

  void SLCSAT_night_Cn2_distribution::print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "       TYPE = " << setw(vlspc) << "SLCSAT night Cn2 distribution"
       << "/" << "object type" << endl;
  }

  Cn2_distribution * SLCSAT_night_Cn2_distribution::clone() const {
    return(new SLCSAT_night_Cn2_distribution());
  }

  bool SLCSAT_night_Cn2_distribution::equals(Cn2_distribution * cn2dist) const {
    if(dynamic_cast<SLCSAT_night_Cn2_distribution *>(cn2dist))
      return(true);
    return(false);
  }

}
