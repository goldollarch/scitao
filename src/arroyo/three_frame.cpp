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

  int three_frame::verbose_level = 0;

  // Choosing 1e-9.  Intrinsic double precision is accurate 
  // to 1e-15.  However, repeated operations on doubles cause
  // numerical accuracy to degrade.  If we assume it degrades
  // as the sqrt of the number of operations, 1e-14 permits
  // about 100 operations before precision is reached.
  // However, this number really is used for geometric accuracy,
  // where all numbers are typically in meters.  So this is .1 nm
  double three_frame::precision = 1e-10;

  three_frame::three_frame() {
    for(int i=0; i<9; i++) basis[i] = 0;

    basis[0] = basis[4] = basis[8] = 1;
  }

  three_frame::three_frame(const three_frame & tf){
    this->operator=(tf);
  }

  three_frame::three_frame(const three_point & tp, const three_vector & x, 
			   const three_vector & y, const three_vector & z) :
    three_point(tp) {

    if(x.length()==0 || y.length()==0 || z.length()==0){
      cerr << "three_frame::three_frame error - "
	   << "null vectors were supplied to this function\n";
      cerr << "\tx length " << x.length() << endl;
      cerr << "\ty length " << y.length() << endl;
      cerr << "\tz length " << z.length() << endl;
      throw(string("three_frame::three_frame"));
    }
    
    double r = x.length();
    three_vector tx(x); tx *= 1/r;
    r = y.length();
    three_vector ty(y); ty *= 1/r;
    r = z.length();
    three_vector tz(z); tz *= 1/r;
  
    if(fabs(dot_product(tx,ty))>precision ||
       fabs(dot_product(ty,tz))>precision ||
       fabs(dot_product(tz,tx))>precision) {
      cerr << "three_frame::three_frame error - "
	   << "the vectors supplied to this function are not mutually orthogonal\n";
      three_frame tf;
      tx.print(cerr, tf, "X hat ");
      ty.print(cerr, tf, "Y hat ");
      tz.print(cerr, tf, "Z hat ");
      cerr << "X dot Y " << dot_product(tx,ty) 
	   << "\tY dot Z " << dot_product(ty,tz) 
	   << "\tZ dot X " << dot_product(tz,tx) 
	   << endl;
      throw(string("three_frame::three_frame"));
    }

    three_frame dflt_frame;

    basis[0] = tx.x(dflt_frame);
    basis[1] = tx.y(dflt_frame);
    basis[2] = tx.z(dflt_frame);

    basis[3] = ty.x(dflt_frame);
    basis[4] = ty.y(dflt_frame);
    basis[5] = ty.z(dflt_frame);

    basis[6] = tz.x(dflt_frame);
    basis[7] = tz.y(dflt_frame);
    basis[8] = tz.z(dflt_frame);

  }

  three_frame & three_frame::operator=(const three_frame & tf) {
    if(this==&tf) 
      return(*this);
    this->three_point::operator=(tf);
    for(int i=0; i<9; i++)
      basis[i] = tf.basis[i];
    return(*this);
  }

  void three_frame::read(const iofits & iof) {
    string comment;
    this->three_point::read(iof);
    iof.read_key("XX", basis[0], comment);
    iof.read_key("XY", basis[1], comment);
    iof.read_key("XZ", basis[2], comment);
    iof.read_key("YX", basis[3], comment);
    iof.read_key("YY", basis[4], comment);
    iof.read_key("YZ", basis[5], comment);
    iof.read_key("ZX", basis[6], comment);
    iof.read_key("ZY", basis[7], comment);
    iof.read_key("ZZ", basis[8], comment);
  }

  void three_frame::write(iofits & iof) const {
    string comment;
    this->three_point::write(iof);
    comment = "X coordinate of X basis vector";
    iof.write_key("XX", basis[0], comment);
    comment = "Y coordinate of X basis vector";
    iof.write_key("XY", basis[1], comment);
    comment = "Z coordinate of X basis vector";
    iof.write_key("XZ", basis[2], comment);
    comment = "X coordinate of Y basis vector";
    iof.write_key("YX", basis[3], comment);
    comment = "Y coordinate of Y basis vector";
    iof.write_key("YY", basis[4], comment);
    comment = "Z coordinate of Y basis vector";
    iof.write_key("YZ", basis[5], comment);
    comment = "X coordinate of Z basis vector";
    iof.write_key("ZX", basis[6], comment);
    comment = "Y coordinate of Z basis vector";
    iof.write_key("ZY", basis[7], comment);
    comment = "Z coordinate of Z basis vector";
    iof.write_key("ZZ", basis[8], comment);

  }

  void three_frame::print(ostream & os, const char * prefix, long precision) const {
    three_frame tf;
    this->print(os, tf, prefix, precision);
  }

  void three_frame::print(ostream & os, const three_frame & tf, const char * prefix, long precision) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    three_point t(static_cast<const three_point>(*this));

    os << prefix << "X          = " << setw(vlspc) << setprecision(precision) << t.x(tf)
       << "/" << "X coordinate of origin" << endl;
    os << prefix << "Y          = " << setw(vlspc) << setprecision(precision) << t.y(tf)
       << "/" << "Y coordinate of origin" << endl;
    os << prefix << "Z          = " << setw(vlspc) << setprecision(precision) << t.z(tf)
       << "/" << "Z coordinate of origin" << endl;
  
    three_vector tmp = this->x();
    os << prefix << "XX         = " << setw(vlspc) << setprecision(precision) << tmp.x(tf)
       << "/" << "X coordinate of X basis vector" << endl;
    os << prefix << "YX         = " << setw(vlspc) << setprecision(precision) << tmp.y(tf)
       << "/" << "Y coordinate of X basis vector" << endl;
    os << prefix << "ZX         = " << setw(vlspc) << setprecision(precision) << tmp.z(tf)
       << "/" << "Z coordinate of X basis vector" << endl;

    tmp = this->y();
    os << prefix << "XY         = " << setw(vlspc) << setprecision(precision) << tmp.x(tf)
       << "/" << "X coordinate of Y basis vector" << endl;
    os << prefix << "YY         = " << setw(vlspc) << setprecision(precision) << tmp.y(tf)
       << "/" << "Y coordinate of Y basis vector" << endl;
    os << prefix << "ZY         = " << setw(vlspc) << setprecision(precision) << tmp.z(tf)
       << "/" << "Z coordinate of Y basis vector" << endl;

    tmp = this->z();
    os << prefix << "XZ         = " << setw(vlspc) << setprecision(precision) << tmp.x(tf)
       << "/" << "X coordinate of Z basis vector" << endl;
    os << prefix << "YZ         = " << setw(vlspc) << setprecision(precision) << tmp.y(tf)
       << "/" << "Y coordinate of Z basis vector" << endl;
    os << prefix << "ZZ         = " << setw(vlspc) << setprecision(precision) << tmp.z(tf)
       << "/" << "Z coordinate of Z basis vector" << endl;

  }

  three_vector three_frame::x() const {
    return(three_vector(basis[0],basis[1],basis[2], three_frame()));
  }

  three_vector three_frame::y() const {
    return(three_vector(basis[3],basis[4],basis[5], three_frame()));
  }

  three_vector three_frame::z() const {
    return(three_vector(basis[6],basis[7],basis[8], three_frame()));
  }

  bool three_frame::right_handed_frame() const {
    three_vector tmp = cross_product(this->x(), this->y());
    if(dot_product(tmp, this->z())>0) return(true);
    return(false);
  }

  bool operator==(const three_frame & tf1, const three_frame & tf2) {
    const three_point * tp1 = dynamic_cast<const three_point *>(&tf1);
    const three_point * tp2 = dynamic_cast<const three_point *>(&tf2);
    if(operator!=(*tp1,*tp2)) return(false);
    for(int i=0; i<9; i++)
      if(tf1.basis[i]!=tf2.basis[i]) return(false);
    return(true);
  }

  bool operator!=(const three_frame & tf1, const three_frame & tf2) {
    return(!operator==(tf1, tf2));
  }

}
