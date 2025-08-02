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
#include <iomanip>
#include "three_frame.h"
#include "three_transformation.h"

using namespace std;

namespace Arroyo {

  int three_transformation::verbose_level = 0;

  three_transformation::three_transformation(){
    xx_ = yy_ = zz_ = 1;
    xy_ = xz_ = yx_ = yz_ = zx_ = zy_ = dx_ = dy_ = dz_ = 0;
  }

  three_transformation::three_transformation(const three_transformation & tt){
    this->operator=(tt);
  }

  three_transformation & three_transformation::operator=(const three_transformation & tt) {
    if(this==&tt) 
      return(*this);
    xx_ = tt.xx_; xy_ = tt.xy_; xz_ = tt.xz_;
    yx_ = tt.yx_; yy_ = tt.yy_; yz_ = tt.yz_;
    zx_ = tt.zx_; zy_ = tt.zy_; zz_ = tt.zz_;
    dx_ = tt.dx_; dy_ = tt.dy_; dz_ = tt.dz_;
    return(*this);
  }

  void three_transformation::read(const iofits & iof) {
    string comment;
    iof.read_key("XX", xx_, comment);
    iof.read_key("XY", xy_, comment);
    iof.read_key("XZ", xz_, comment);
    iof.read_key("YX", yx_, comment);
    iof.read_key("YY", yy_, comment);
    iof.read_key("YZ", yz_, comment);
    iof.read_key("ZX", zx_, comment);
    iof.read_key("ZY", zy_, comment);
    iof.read_key("ZZ", zz_, comment);
    iof.read_key("DX", dx_, comment);
    iof.read_key("DY", dy_, comment);
    iof.read_key("DZ", dz_, comment);
  }

  void three_transformation::write(iofits & iof) const {
    string comment = "element of transformation matrix";
    iof.write_key("XX", xx_, comment);
    iof.write_key("XY", xy_, comment);
    iof.write_key("XZ", xz_, comment);
    iof.write_key("YX", yx_, comment);
    iof.write_key("YY", yy_, comment);
    iof.write_key("YZ", yz_, comment);
    iof.write_key("ZX", zx_, comment);
    iof.write_key("ZY", zy_, comment);
    iof.write_key("ZZ", zz_, comment);
    iof.write_key("DX", dx_, comment);
    iof.write_key("DY", dy_, comment);
    iof.write_key("DZ", dz_, comment);
  }

  void three_transformation::print(ostream & os, const char * prefix) const {
    int w = 20;
    os.setf(ios::right, ios::adjustfield); 
    cout << prefix << setw(w) << xx_ << setw(w) << yx_ << setw(w) << zx_ << setw(w) << dx_ << endl;
    cout << prefix << setw(w) << xy_ << setw(w) << yy_ << setw(w) << zy_ << setw(w) << dy_ << endl;
    cout << prefix << setw(w) << xz_ << setw(w) << yz_ << setw(w) << zz_ << setw(w) << dz_ << endl;
  }

  three_transformation three_transformation::inverse() const {

    double det = xx_*(yy_*zz_-yz_*zy_) - yx_*(xy_*zz_-xz_*zy_) + zx_*(xy_*yz_-xz_*yy_);
 
    if(det==0){
      cerr << "three_transformation::inverse error - "
	   << "transformation is not invertible\n";
      throw(string("three_transformation::inverse"));
    }

    // this function relies on the fact that one can decompose an
    // arbitrary transformation into a translation times a
    // transformation that has no translation component (i.e. the final
    // column has all zeroes).  
    //
    // Here we invert each of these components separately, and return the
    // product of the inverses.  
    // Note that in this scheme where we include the translation component,
    // operations are not commutative

    double idet = 1/det;

    // To get the inverse, for each entry in the 3x3 matrix eliminate
    // the row and column containing the entry and form the
    // determinant of the leftover 2x2 matrix.  Count the steps to the
    // entry from the upper left corner by moving horizontally or
    // vertically along columns or rows, and multiply it by -1 raised
    // to the power of the number of steps.
    three_transformation tr(idet*(yy_*zz_-yz_*zy_), idet*(-zz_*xy_+xz_*zy_), idet*(xy_*yz_-xz_*yy_),0,
			    idet*(-zz_*yx_+zx_*yz_), idet*(xx_*zz_-xz_*zx_), idet*(-xx_*yz_+yx_*xz_),0,
			    idet*(yx_*zy_-zx_*yy_), idet*(-xx_*zy_+zx_*xy_), idet*(xx_*yy_-xy_*yx_),0);

    three_transformation tt(1,0,0,-dx_,
			    0,1,0,-dy_,
			    0,0,1,-dz_);  

    return(tr*tt);
  }

  three_transformation & three_transformation::operator*=(const three_transformation & tt){
    *this = tt * *this;
    return(*this);
  }

  three_transformation operator*(const three_transformation & tt1, const three_transformation & tt2){

    return(three_transformation(tt1.xx_*tt2.xx_+tt1.xy_*tt2.yx_+tt1.xz_*tt2.zx_,
				tt1.xx_*tt2.xy_+tt1.xy_*tt2.yy_+tt1.xz_*tt2.zy_,
				tt1.xx_*tt2.xz_+tt1.xy_*tt2.yz_+tt1.xz_*tt2.zz_,
				tt1.xx_*tt2.dx_+tt1.xy_*tt2.dy_+tt1.xz_*tt2.dz_ + tt1.dx_,
				tt1.yx_*tt2.xx_+tt1.yy_*tt2.yx_+tt1.yz_*tt2.zx_,
				tt1.yx_*tt2.xy_+tt1.yy_*tt2.yy_+tt1.yz_*tt2.zy_,
				tt1.yx_*tt2.xz_+tt1.yy_*tt2.yz_+tt1.yz_*tt2.zz_,
				tt1.yx_*tt2.dx_+tt1.yy_*tt2.dy_+tt1.yz_*tt2.dz_ + tt1.dy_,
				tt1.zx_*tt2.xx_+tt1.zy_*tt2.yx_+tt1.zz_*tt2.zx_,
				tt1.zx_*tt2.xy_+tt1.zy_*tt2.yy_+tt1.zz_*tt2.zy_,
				tt1.zx_*tt2.xz_+tt1.zy_*tt2.yz_+tt1.zz_*tt2.zz_,
				tt1.zx_*tt2.dx_+tt1.zy_*tt2.dy_+tt1.zz_*tt2.dz_ + tt1.dz_));
  }

  bool operator==(const three_transformation & tt1, const three_transformation & tt2) {
    if(tt1.xx_ == tt2.xx_ && tt1.xy_ == tt2.xy_ && tt1.xz_ == tt2.xz_ && 
       tt1.yx_ == tt2.yx_ && tt1.yy_ == tt2.yy_ && tt1.yz_ == tt2.yz_ && 
       tt1.zx_ == tt2.zx_ && tt1.zy_ == tt2.zy_ && tt1.zz_ == tt2.zz_ &&
       tt1.dx_ == tt2.dx_ && tt1.dy_ == tt2.dy_ && tt1.dz_ == tt2.dz_)
      return(true);
    return(false);
  }

  bool operator!=(const three_transformation & tt1, const three_transformation & tt2) {
    return(!operator==(tt1, tt2));
  }

  void three_transformation::transform(three_point & tp) const {
    double tmpx = tp.x_, tmpy = tp.y_, tmpz=tp.z_;
    tp.x_ = xx_*tmpx + xy_*tmpy + xz_*tmpz + dx_; 
    tp.y_ = yx_*tmpx + yy_*tmpy + yz_*tmpz + dy_;
    tp.z_ = zx_*tmpx + zy_*tmpy + zz_*tmpz + dz_;

    /*
    tp.x_ = xx_*tmpx + yx_*tmpy + zx_*tmpz + dx_;
    tp.y_ = xy_*tmpx + yy_*tmpy + zy_*tmpz + dy_; 
    tp.z_ = xz_*tmpx + yz_*tmpy + zz_*tmpz + dz_;
    */
  } 

  void three_transformation::transform(three_vector & tv) const {
    double tmpx = tv.x_, tmpy = tv.y_, tmpz=tv.z_;
    tv.x_ = xx_*tmpx + xy_*tmpy + xz_*tmpz;
    tv.y_ = yx_*tmpx + yy_*tmpy + yz_*tmpz; 
    tv.z_ = zx_*tmpx + zy_*tmpy + zz_*tmpz;
 
    /*
    tv.x_ = xx_*tmpx + yx_*tmpy + zx_*tmpz;
    tv.y_ = xy_*tmpx + yy_*tmpy + zy_*tmpz; 
    tv.z_ = xz_*tmpx + yz_*tmpy + zz_*tmpz;
    */
  }

  void three_orthonormal_transformation::transform(three_frame & tf) const {
  
    three_orthonormal_transformation inv_tf = this->inverse();

    three_point tp(static_cast<three_point>(tf));
    three_vector tmpx = tf.x(), tmpy = tf.y(), tmpz = tf.z();
    inv_tf.transform(tp);  
    inv_tf.transform(tmpx);  
    inv_tf.transform(tmpy);  
    inv_tf.transform(tmpz);  
    tf = three_frame(tp, tmpx, tmpy, tmpz);
  }

  three_orthonormal_transformation three_orthonormal_transformation::inverse() const {
    three_transformation tt = this->three_transformation::inverse();
    three_orthonormal_transformation tot;
    tot.set_transformation(tt);
    return(tot);
  }

  three_translation::three_translation(const three_vector & v){
    three_frame tf;
    set_three_transformation(1,0,0, v.x(tf),
			     0,1,0, v.y(tf),
			     0,0,1, v.z(tf));
  }

  three_translation three_translation::inverse() const {
    three_transformation tt = this->inverse();
    three_translation ttrans;
    ttrans.set_transformation(tt);
    return(ttrans);
  }

  three_rotation::three_rotation(const three_point & tp, const three_vector & tv, double angle){

    if(tv.length_squared()==0){
      cerr << "three_rotation::three_rotation error - "
	   << "three_vector supplied to this function is a null vector\n";
      throw(string("three_rotation::three_rotation"));
    }

    three_vector tmpv(tv);
    tmpv *= 1/tmpv.length();
    three_frame tf;
    double px = tp.x(tf), py = tp.y(tf), pz = tp.z(tf);
    double vx = tmpv.x(tf), vy = tmpv.y(tf), vz = tmpv.z(tf); 
    double ca = cos(angle), sa = sin(angle);

    set_three_transformation(   ca + (1-ca)*vx*vx, (1-ca)*vx*vy + sa*vz, (1-ca)*vx*vz - sa*vy, 
				px*(1-ca - (1-ca)*vx*vx) - py*((1-ca)*vx*vy + sa*vz) - pz*((1-ca)*vx*vz - sa*vy), 
				(1-ca)*vy*vx - sa*vz,    ca + (1-ca)*vy*vy, (1-ca)*vy*vz + sa*vx,  
				-px*((1-ca)*vy*vx - sa*vz) + py*(1-ca - (1-ca)*vy*vy) - pz*((1-ca)*vy*vz + sa*vx),
				(1-ca)*vz*vx + sa*vy, (1-ca)*vz*vy - sa*vx,    ca + (1-ca)*vz*vz,
				-px*((1-ca)*vz*vx + sa*vy) - py*((1-ca)*vz*vy - sa*vx) + pz*(1-ca - (1-ca)*vz*vz));
  }

  three_rotation three_rotation::inverse() const {
    three_transformation tt = this->inverse();
    three_rotation tr;
    tr.set_transformation(tt);
    return(tr);
  }

  three_reflection::three_reflection(const three_point & tp, const three_vector & tv){

    if(tv.length_squared()==0){
      cerr << "three_reflection::three_reflection error - "
	   << "three_vector supplied to this function is a null vector\n";
      throw(string("three_reflection::three_reflection"));
    }

    three_vector tmpv(tv);
    tmpv *= 1/tmpv.length();
    three_frame tf;
    double px = tp.x(tf), py = tp.y(tf), pz = tp.z(tf);
    double vx = tmpv.x(tf), vy = tmpv.y(tf), vz = tmpv.z(tf); 
    double dot = px*vx+py*vy+pz*vz;

    set_three_transformation(1-2*vx*vx,  -2*vx*vy,  -2*vx*vz, 2*dot*vx,
			     -2*vy*vx, 1-2*vy*vy,  -2*vy*vz, 2*dot*vy,
			     -2*vz*vx,  -2*vz*vy, 1-2*vz*vz, 2*dot*vz);
  }

  three_reflection three_reflection::inverse() const {
    three_transformation tt = this->inverse();
    three_reflection tr;
    tr.set_transformation(tt);
    return(tr);
  }

  three_scaling::three_scaling(const three_point & tp, double scale){

    if(scale<=0){
      cerr << "three_scaling::three_scaling error - "
	   << "scale factor " << scale << " supplied to this function "
	   << "is not a positive number\n";
      throw(string("three_scaling::three_scaling"));
    }

    three_frame tf;
    set_three_transformation(scale,0,0,tp.x(tf)*(1-scale),
			     0,scale,0,tp.y(tf)*(1-scale),
			     0,0,scale,tp.z(tf)*(1-scale));
  }

  three_scaling::three_scaling(const three_point & tp, const three_vector & tv, double scale){

    if(scale<=0){
      cerr << "three_scaling::three_scaling error - "
	   << "scale factor " << scale << " supplied to this function "
	   << "is not a positive number\n";
      throw(string("three_scaling::three_scaling"));
    }

    if(tv.length_squared()==0){
      cerr << "three_scaling::three_scaling error - "
	   << "three_vector supplied to this function is a null vector\n";
      throw(string("three_scaling::three_scaling"));
    }

    three_vector tmpv(tv);
    tmpv *= 1/tmpv.length();
    three_frame tf;
    double px = tp.x(tf), py = tp.y(tf), pz = tp.z(tf);
    double vx = tmpv.x(tf), vy = tmpv.y(tf), vz = tmpv.z(tf); 
    double dot = px*vx+py*vy+pz*vz;
    double fac = 1-scale;

    set_three_transformation(1 - fac*vx*vx,    -fac*vx*vy,    -fac*vx*vz, fac*dot*vx,
			     -fac*vx*vy, 1 - fac*vy*vy,    -fac*vy*vz, fac*dot*vy,
			     -fac*vx*vz,    -fac*vy*vz, 1 - fac*vz*vz, fac*dot*vz);
  }

  three_scaling three_scaling::inverse() const {
    three_transformation tt = this->inverse();
    three_scaling ts;
    ts.set_transformation(tt);
    return(ts);
  }

}
