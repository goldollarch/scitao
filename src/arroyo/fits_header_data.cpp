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
#include <string>
#include "fits_header_data.h"
#include "iofits.h"

using namespace std;

namespace Arroyo {

  fits_scale_factor::fits_scale_factor(){
    bscale = 0;     
    bzero = 0;
  }

  fits_scale_factor::fits_scale_factor(const iofits & iof){
    iof.read_key("BSCALE",   bscale,    bscale_comment);
    iof.read_key("BZERO",    bzero,     bzero_comment);
  }

  fits_scale_factor & fits_scale_factor::operator= (const fits_scale_factor & fhd){

    if(this == &fhd) 
      return(*this);
    bscale    = fhd.bscale;         bscale_comment    = fhd.bscale_comment;   
    bzero     = fhd.bzero;          bzero_comment     = fhd.bzero_comment;    
    return(*this);
  }

  void fits_scale_factor::write(iofits & iof) const {
    iof.write_key("BSCALE",   bscale,    bscale_comment);
    iof.write_key("BZERO",    bzero,     bzero_comment);
  }

  void fits_scale_factor::print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    int i;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "BSCALE     = " << setw(vlspc) << bscale      
       << "/" << bscale_comment      << endl;
    os << prefix << "BZERO      = " << setw(vlspc) << bzero 
       << "/" << bzero_comment      << endl;
  }

  template<>
  iofits::imagetype fits_header_data<short>::get_image_type() const {
    return(iofits::SHORTIMG);
  }

  template<>
  iofits::imagetype fits_header_data<unsigned short>::get_image_type() const {
    return(iofits::USHORTIMG);
  }

  template<>
  iofits::imagetype fits_header_data<long>::get_image_type() const {
    return(iofits::LONGIMG);
  }

  template<>
  iofits::imagetype fits_header_data<char>::get_image_type() const {
    return(iofits::BYTEIMG);
  }

  template<>
  iofits::imagetype fits_header_data<float>::get_image_type() const {
    return(iofits::FLOATIMG);
  }

  template<>
  iofits::imagetype fits_header_data<double>::get_image_type() const {
    return(iofits::DOUBLEIMG);
  }
}
