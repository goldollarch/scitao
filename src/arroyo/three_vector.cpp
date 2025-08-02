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
#include "three_point.h"
#include "three_vector.h"
#include "three_frame.h"

using namespace std;
namespace Arroyo {

  int three_vector::verbose_level = 0;

  three_vector::three_vector(const three_vector & tv){
    this->operator=(tv);
  }

  three_vector::three_vector(double x, double y, double z){
    x_ = x;
    y_ = y;
    z_ = z;
  }

  three_vector::three_vector(double x, double y, double z, const three_frame & tf){
    this->operator=(three_vector(x*tf.basis[0]+y*tf.basis[3]+z*tf.basis[6], 
				 x*tf.basis[1]+y*tf.basis[4]+z*tf.basis[7], 
				 x*tf.basis[2]+y*tf.basis[5]+z*tf.basis[8]));
  }

  three_vector & three_vector::operator=(const three_vector & tv) {
    if(this==&tv) 
      return(*this);
    x_ = tv.x_;
    y_ = tv.y_;
    z_ = tv.z_;
    return(*this);
  }

  void three_vector::read(const iofits & iof) {
    string comment;
    iof.read_key("XCMPNT", x_, comment);
    iof.read_key("YCMPNT", y_, comment);
    iof.read_key("ZCMPNT", z_, comment);
  }

  void three_vector::write(iofits & iof) const {
    string comment = "X component of vector";
    iof.write_key("XCMPNT", x_, comment);
    comment = "Y component of vector";
    iof.write_key("YCMPNT", y_, comment); 
    comment = "Z component of vector";
    iof.write_key("ZCMPNT", z_, comment);
  }

  void three_vector::print(ostream & os, const char * prefix, long precision) const {
    three_frame tf;
    this->print(os, tf, prefix, precision);
  }

  void three_vector::print(ostream & os, const three_frame & tf, const char * prefix, long precision) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "XCMPNT     = " << setw(vlspc) << setprecision(precision) << this->x(tf)
       << "/" << "X component of vector" << endl;
    os << prefix << "YCMPNT     = " << setw(vlspc) << setprecision(precision) << this->y(tf)
       << "/" << "Y component of vector" << endl;
    os << prefix << "ZCMPNT     = " << setw(vlspc) << setprecision(precision) << this->z(tf)
       << "/" << "Z component of vector" << endl;
  }

  double three_vector::x(const three_frame & tf) const {
    return(dot_product(*this, tf.x()));
  }

  double three_vector::y(const three_frame & tf) const {
    return(dot_product(*this, tf.y()));
  }

  double three_vector::z(const three_frame & tf) const {
    return(dot_product(*this, tf.z()));
  }

  three_vector & three_vector::operator+=(const three_vector & tv) {
    x_ += tv.x_;
    y_ += tv.y_;
    z_ += tv.z_;
    return(*this);
  }

  three_vector & three_vector::operator-=(const three_vector & tv) {
    x_ -= tv.x_;
    y_ -= tv.y_;
    z_ -= tv.z_;
    return(*this);
  }
 
  three_vector & three_vector::operator*=(double d){
    x_ *= d;
    y_ *= d;
    z_ *= d;
    return(*this);
  }

  bool operator==(const three_vector & tv1, const three_vector & tv2) {
    if(fabs(tv1.x_-tv2.x_)<three_frame::precision &&
       fabs(tv1.y_-tv2.y_)<three_frame::precision &&
       fabs(tv1.z_-tv2.z_)<three_frame::precision) 
      return(true);
    return(false);
  }

  bool operator!=(const three_vector & tv1, const three_vector & tv2) {
    return(!operator==(tv1, tv2));
  }

  three_vector operator+(const three_vector & tv1, const three_vector & tv2){
    three_vector tv(tv1);
    tv += tv2;
    return(tv);
  }
 
  three_vector operator-(const three_vector & tv1, const three_vector & tv2){
    three_vector tv(tv1);
    tv -= tv2;
    return(tv);
  }

  three_vector operator*(double d, const three_vector & tv){
    three_vector v(tv);
    return(v*=d);
  }
 
  three_vector operator*(const three_vector & tv, double d){
    three_vector v(tv);
    return(v*=d);
  }
 
  three_vector parallel_projection(const three_vector & tv, 
				   const three_vector & nrml, 
				   const three_vector & drctn){
    
    if(fabs(dot_product(nrml, drctn))<three_frame::precision){
      cerr << "parallel_projection error - cannot project a vector along a vector orthogonal to it\n";
      nrml.print(cerr, "normal ");
      tv.print(cerr, "vector to project ");
      drctn.print(cerr, "projection direction vector ");
      throw(string("parallel_projection"));
    }
    return(tv - (dot_product(tv,nrml)/dot_product(drctn,nrml))*drctn);
  }

}
