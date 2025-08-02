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

#ifndef DEREFERENCED_H
#define DEREFERENCED_H

namespace Arroyo {

  ///////////////////////////////////////////
  ///  Template class to test equality of objects
  ///  referred to by pointers
  template<class T>
    class DereferencedEqual {
    public:
    DereferencedEqual(const T* p) : p_(p) { }
    bool operator() (const T* p2) const { return *p_ == *p2; }
    private:
    const T* p_;
  }; 

  ///////////////////////////////////////////
  ///  Template class to test operator> for objects
  ///  referred to by pointers
  template<class T>
    class DereferencedGreater {
    public:
    DereferencedGreater() { p_ = NULL;};
    DereferencedGreater(const T* p) : p_(p) { }
    bool operator() (const T* p1, const T* p2) { return *p1 > *p2; }
    private:
    const T* p_;
  }; 

  ///////////////////////////////////////////
  ///  Template class to test operator< for objects
  ///  referred to by pointers
  template<class T>
    class DereferencedLess {
    public:
    DereferencedLess() { p_ = NULL;};
    DereferencedLess(const T* p) : p_(p) { }
    bool operator() (const T* p1, const T* p2) const { return *p1 < *p2; }
    private:
    const T* p_;
  }; 

}

#endif
