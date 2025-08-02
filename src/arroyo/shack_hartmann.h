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

#ifndef SHACK_HARTMANN_H
#define SHACK_HARTMANN_H

namespace Arroyo {

  ///
  /// A class to hold a shack-hartmann wavefront sensor.
  ///

  class shack_hartmann_wavefront_sensor :
    public wavefront_sensor {

    protected:
    /// The rectangular lenslet array
    rectangular_lenslet_array rec_lens_array;
    /// The detector - could be optical or IR
    detector sh_dtctr;

    public:

    ///////////////////////////////////////////
    ///  Null constructor
    shack_hartmann_wavefront_sensor();

    ///////////////////////////////////////////
    ///  Copy constructor
    shack_hartmann_wavefront_sensor(const shack_hartmann_wavefront_sensor & shwfs);

    ///////////////////////////////////////////
    ///  Construct from iofits object
    shack_hartmann_wavefront_sensor(iofits & iof);

    ///////////////////////////////////////////
    ///  Construct from iofits object
    shack_hartmann_wavefront_sensor(const char * filename);

    ///////////////////////////////////////////
    ///  Destructor
    ~shack_hartmann_wavefront_sensor();

    ///////////////////////////////////////////
    ///  Operator = 
    shack_hartmann_wavefront_sensor & operator=(const shack_hartmann_wavefront_sensor & shwfs);

    ///////////////////////////////////////////
    ///  Read from iofits object
    void read(iofits & iof);

    ///////////////////////////////////////////
    ///  Write to iofits object
    void write(iofits & iof) const;

    ///////////////////////////////////////////
    ///  Detect a wavefront
    ///  performs the propagation through the 
    ///  lenslet array and the detection 
    ///  by the optical_detector
    void detect(const wavefront & wf);

    ///////////////////////////////////////////
    ///  Convolve with an emitter
    void convolve(const extended_emitter & xtnd_emtr);

    ///////////////////////////////////////////
    ///  Detect a wavefront
    pixel_amp_array readout() const;

    ///////////////////////////////////////////
    ///  Clear detector
    void clear();

  };

}

#endif
