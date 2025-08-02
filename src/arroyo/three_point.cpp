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

  int three_point::verbose_level = 0;

  three_point::three_point(const three_point & tp){
    this->operator=(tp);
  }

  three_point::three_point(double x, double y, double z){
    x_ = x;
    y_ = y;
    z_ = z;
  }

  three_point::three_point(double x, double y, double z, const three_frame & tf){
    three_vector v(x,y,z,tf);
    v += (tf-three_point());
    this->operator=(three_point());
    this->operator+=(v);
  }

  three_point & three_point::operator=(const three_point & tp) {
    if(this==&tp) 
      return(*this);
    x_ = tp.x_;
    y_ = tp.y_;
    z_ = tp.z_;
    return(*this);
  }

  void three_point::read(const iofits & iof) {
    string comment;
    iof.read_key("XCOORD", x_, comment);
    iof.read_key("YCOORD", y_, comment);
    iof.read_key("ZCOORD", z_, comment);
  }

  void three_point::write(iofits & iof) const {
    string comment = "X coordinate of point";
    iof.write_key("XCOORD", x_, comment);
    comment = "Y coordinate of point";
    iof.write_key("YCOORD", y_, comment); 
    comment = "Z coordinate of point";
    iof.write_key("ZCOORD", z_, comment);
  }

  void three_point::print(ostream & os, const char * prefix, long precision) const {
    three_frame tf;
    this->print(os, tf, prefix, precision);
  }

  void three_point::print(ostream & os, const three_frame & tf, const char * prefix, long precision) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "XCOORD     = " << setw(vlspc) << setprecision(precision) << this->x(tf)
       << "/" << "X coordinate of point" << endl;
    os << prefix << "YCOORD     = " << setw(vlspc) << setprecision(precision) << this->y(tf)
       << "/" << "Y coordinate of point" << endl;
    os << prefix << "ZCOORD     = " << setw(vlspc) << setprecision(precision) << this->z(tf)
       << "/" << "Z coordinate of point" << endl;
  }

  double three_point::x(const three_frame & tf) const {
    // *this-tf is three_point subtraction that results in
    // a vector defining this point wrt the tf frame
    // The dot product then gets the component along the x axis
    // in the tf frame
    return(dot_product(*this-tf,tf.x()));
  }

  double three_point::y(const three_frame & tf) const {
    return(dot_product(*this-tf,tf.y()));
  }

  double three_point::z(const three_frame & tf) const {
    return(dot_product(*this-tf,tf.z()));
  }

  three_point & three_point::operator+=(const three_vector & tv){
    x_ += tv.x_;
    y_ += tv.y_;
    z_ += tv.z_;
    return(*this);
  }

  three_point & three_point::operator-=(const three_vector & tv){
    x_ -= tv.x_;
    y_ -= tv.y_;
    z_ -= tv.z_;
    return(*this);
  }

  three_vector operator-(const three_point & tp1, const three_point & tp2){
    return(three_vector(tp1.x_-tp2.x_,tp1.y_-tp2.y_,tp1.z_-tp2.z_));
  }

  bool operator==(const three_point & tp1, const three_point & tp2) {
    if(fabs(tp1.x_-tp2.x_) < three_frame::precision &&
       fabs(tp1.y_-tp2.y_) < three_frame::precision &&
       fabs(tp1.z_-tp2.z_) < three_frame::precision)
      return(true);
    return(false);
  }

  bool operator!=(const three_point & tp1, const three_point & tp2) {
    return(!operator==(tp1, tp2));
  }

  three_point operator+(const three_point & tp, const three_vector & tv) {
    three_point tmp(tp);
    tmp += tv;
    return(tmp);
  }

  three_point operator-(const three_point & tp, const three_vector & tv) {
    three_point tmp(tp);
    tmp -= tv;
    return(tmp);
  }

}
