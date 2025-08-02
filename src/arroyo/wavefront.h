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

#ifndef WAVEFRONT_H
#define WAVEFRONT_H

#include "AO_sim_base.h"
#include "wavefront_header.h"

namespace Arroyo {

  using std::ostream;

  class optic;

  ///
  /// A virtual base class for diffractive and geometric wavefronts
  ///

  class wavefront :
    virtual public AO_sim_base, 
    virtual public wavefront_header {
  
    public:
  
    ///////////////////////////////////////////
    ///  Null constructor
    wavefront(){};
  
    ///////////////////////////////////////////
    ///  Virtual destructor
    virtual ~wavefront(){};

    ///////////////////////////////////////////
    ///  Virtual read from file
    virtual void read(const char * filename) = 0;

    ///////////////////////////////////////////
    ///  Virtual read from iofits object
    virtual void read(const iofits & iof) = 0;

    ///////////////////////////////////////////
    ///  Virtual write to file
    virtual void write(const char * filename) const = 0;

    ///////////////////////////////////////////
    ///  Virtual write to iofits object
    virtual void write(iofits & iof) const = 0;

    ///////////////////////////////////////////
    ///  Virtual print
    virtual void print(ostream & os, const char * prefix) const = 0;

    ///virtual propagate(const optic & op) = 0;

    ////virtual propagate(const detector & dtr) = 0;

    ///////////////////////////////////////////
    ///  Factory constructor from file
    static wavefront * wavefront_factory(const char * filename);
  
    ///////////////////////////////////////////
    ///  Factory constructor from iofits
    static wavefront * wavefront_factory(const iofits & iof);

    ///////////////////////////////////////////
    ///  Verbose level
    static int verbose_level;

  };

}

#endif

