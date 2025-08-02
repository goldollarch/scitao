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

#ifndef WAVEFRONT_SENSOR_H
#define WAVEFRONT_SENSOR_H

namespace Arroyo {

  ///
  /// A virtual base class to hold a wavefront sensor.
  ///

  class wavefront_sensor {

  public:
    ///////////////////////////////////////////
    ///  Null constructor
    wavefront_sensor(){};

    ///////////////////////////////////////////
    ///  Destructor
    virtual ~wavefront_sensor() = 0;

    ///////////////////////////////////////////
    ///  Read from iofits object
    virtual void read(iofits & iof) = 0;

    ///////////////////////////////////////////
    ///  Write to iofits object
    virtual void write(iofits & iof) const = 0;

    ///////////////////////////////////////////
    ///  Detect a wavefront
    virtual void detect(const wavefront & wf) = 0;

    ///////////////////////////////////////////
    ///  Convolve with an emitter
    virtual void convolve(const extended_emitter & xtnd_emtr) = 0;

  };

}

#endif
