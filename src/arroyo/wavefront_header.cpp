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
#include "wavefront_header.h"

using namespace std;

namespace Arroyo {

  int wavefront_header::verbose_level = 0;

  void wavefront_header::read(const iofits & iof) {
    string comment;
    iof.read_key("TIMESTMP", timestamp, comment);
  }

  void wavefront_header::write(iofits & iof) const {
    string comment = "timestamp (secs)";
    iof.write_key("TIMESTMP", timestamp, comment);
  }

  void wavefront_header::print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "TIMESTMP   = " << setw(vlspc) << timestamp
       << "/" << "timestamp (secs)" << endl;
  }

  wavefront_header & wavefront_header::operator=(const wavefront_header & wfh){
    if(this==&wfh)
      return(*this);
    timestamp = wfh.timestamp;
    return(*this);
  }

  wavefront_header * wavefront_header::wavefront_header_factory(const char * filename) {
    iofits iof;
    try{iof.open(filename);}
    catch(...) {
      cerr << "wavefront_header::wavefront_header_factory - "
	   << "error opening file " << filename << endl;
      throw(string("wavefront_header::wavefront_header_factory"));
    }
    return(wavefront_header::wavefront_header_factory(iof));
  }

  wavefront_header * wavefront_header::wavefront_header_factory(const iofits & iof) {
    string type, comment;
    iof.read_key("TYPE", type, comment);
    //if(type=="geometric wavefront")
    //  return(new geometric_wavefront_header(iof));
    long bitpix;
    iof.read_key("BITPIX", bitpix, comment);
    if(bitpix==iofits::DOUBLEIMG)
      return(new diffractive_wavefront_header<double>(iof));
    else if(bitpix==iofits::FLOATIMG)
      return(new diffractive_wavefront_header<float>(iof));
  
    cerr << "wavefront_header::wavefront_header_factory error - "
	 << "could not retrieve wavefront_header from iofits object\n";
    throw(string("wavefront_header::wavefront_header_factory"));
  }

  int operator<(const wavefront_header & wfh1, const wavefront_header & wfh2) {
    if(wfh1.timestamp < wfh2.timestamp) return(true);
    return(false);
  }

  int operator>(const wavefront_header & wfh1, const wavefront_header & wfh2) {
    if(wfh1.timestamp > wfh2.timestamp) return(true);
    return(false);
  }

  bool operator==(const wavefront_header & wfh1, const wavefront_header & wfh2) {
    if(wfh1.timestamp==wfh2.timestamp) return(true);
    return(false);
  }

  bool operator!=(const wavefront_header & wfh1, const wavefront_header & wfh2) {
    if(operator==(wfh1,wfh2)) return(false);
    return(true);
  }

  /*
  geometric_wavefront_header::geometric_wavefront_header() {
    nrays = 0;
    wavelength = 0;
  }

  geometric_wavefront_header::geometric_wavefront_header(const iofits & iof) {
    this->read(iof);
  }

  geometric_wavefront_header::geometric_wavefront_header(const geometric_wavefront_header & wfh) {
    nrays = 0;
    this->operator=(wfh);
  }

  geometric_wavefront_header::geometric_wavefront_header(long in_nrays, double in_wavelength){
    nrays = in_nrays;
    wavelength = in_wavelength;  
  }

  geometric_wavefront_header & geometric_wavefront_header::operator=(const geometric_wavefront_header & wfh) {
    if(this==&wfh)
      return(*this);
    this->wavefront_header::operator=(wfh);
    nrays = wfh.nrays;
    wavelength = wfh.wavelength;
    return(*this);
  }

  void geometric_wavefront_header::read(const iofits & iof) {
    if(!iof.key_exists("TYPE")){
      cerr << "geometric_wavefront_header::read error - "
	   << "unrecognized type of file\n";
      throw(string("geometric_wavefront::read"));
    }
    string type, comment;
    iof.read_key("TYPE", type, comment);
    if(type!="geometric wavefront"){
      cerr << "geometric_wavefront_header::read error - file of type " 
	   << type << " rather than of type geometric_wavefront\n";
      throw(string("geometric_wavefront::read"));
    }
    this->wavefront_header::read(iof);
    iof.read_key("NRAYS", nrays, comment);
    iof.read_key("WVLNGTH", wavelength, comment);
  }

  void geometric_wavefront_header::write(iofits & iof) const {

    fits_header_data tmphdr(iofits::DOUBLEIMG, vector<long>(1,nrays*8));
    tmphdr.write(iof);

    string type = "geometric wavefront";
    string comment = "object type";
    iof.write_key("TYPE", type, comment);
    this->wavefront_header::write(iof);
    iof.write_key("NRAYS", nrays, comment);
    iof.write_key("WVLNGTH", wavelength, "wavelength (microns)");
  }

  void geometric_wavefront_header::print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << "geometric wavefront"
       << "/" << "object type" << endl;
    this->wavefront_header::print(os, prefix);
    os << prefix << "NRAYS      = " << setw(vlspc) << nrays
       << "/" << "wavelength (microns)" << endl;
    os << prefix << "WVLNGTH    = " << setw(vlspc) << wavelength*1e6
       << "/" << "wavelength (microns)" << endl;
  }

  bool operator ==(const geometric_wavefront_header & wfh1, const geometric_wavefront_header & wfh2){
    if(wfh1.timestamp!=wfh2.timestamp) return(false);
    if(wfh1.nrays!=wfh2.nrays) return(false);
    if(wfh1.wavelength!=wfh2.wavelength) return(false);
    return(true);
  }

  bool operator !=(const geometric_wavefront_header & wfh1, const geometric_wavefront_header & wfh2){
    return(!(wfh1==wfh2));
  }
  */

}
