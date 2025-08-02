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

#ifndef OPTICAL_SYSTEM_H
#define OPTICAL_SYSTEM_H

#include "optic.h"

namespace Arroyo {

  ///
  /// A virtual base class to hold a system of optics.
  ///

  class optical_system :
    public vector<optic> {

    protected:
  
    /// The near field propagation method
    int near_field_propagator;

    /// The far field propagation method
    int far_field_propagator;

    /// Flag to specify whether the
    /// current optic has been applied to
    /// the wavefront
    mutable bool current_optic_applied;

    public:

    ///////////////////////////////////////////
    ///  The near field propagation method
    enum near_field_propagation_method {GEOMETRIC, FRESNEL, VERGING_FRESNEL, EXACT};

    ///////////////////////////////////////////
    ///  The far field propagation method
    enum far_field_propagation_method {GEOMETRIC, FRESNEL, EXACT};

    ///////////////////////////////////////////
    ///  Null constructor
    optical_system();

    ///////////////////////////////////////////
    ///  Copy constructor
    optical_system(const optical_system & op_sys);

    ///////////////////////////////////////////
    ///  Construct from a file
    optical_system(const char * filename);

    ///////////////////////////////////////////
    ///  Construct from an iofits object
    optical_system(iofits & iof);

    ///////////////////////////////////////////
    ///  Destructor
    ~optical_system();

    ///////////////////////////////////////////
    ///  Operator = 
    optical_system & operator=(const optical_system & opt_sys);

    ///////////////////////////////////////////
    ///  Prepend an optical_system to the optical system
    void prepend(optical_system & pre_opt_sys);

    ///////////////////////////////////////////
    /// Insert an optical_system into the optical 
    /// system after the first pos optics.
    ///    pos < 0                 an error
    ///    pos == 0                equivalent to prepend 
    ///    pos > number of optics  equivalent to append
    void insert(optical_system & in_opt_sys, long pos);

    ///////////////////////////////////////////
    ///  Append an optical_system to the optical system
    void append(optical_system & post_opt_sys);

    ///////////////////////////////////////////
    ///  Choose near field propagation preference
    void set_near_field_propagation(near_field_propagation_method & nrfld_prop_method);
  
    ///////////////////////////////////////////
    ///  Choose far field propagation preference
    void set_far_field_propagation(far_field_propagation_method & frfld_prop_method);
  
    ///////////////////////////////////////////
    ///  Propagate the wavefront forward through the 
    ///  entire optical_system
    void propagate_forward(wavefront & wf) const;
  
    ///////////////////////////////////////////
    ///  Propagate the wavefront backward through the 
    ///  entire optical_system 
    void propagate_backward(wavefront & wf) const;
  
    ///////////////////////////////////////////
    ///  Propagate the wavefront forward through the 
    ///  optical_system from the current optic
    ///  to the end.
    void propagate_to_end(wavefront & wf) const;
  
    ///////////////////////////////////////////
    ///  Propagate the wavefront backward through the 
    ///  optical_system from the current optic
    ///  to the beginning.
    void propagate_to_beginning(wavefront & wf) const;
  
    ///////////////////////////////////////////
    ///  Propagate the wavefront forward
    ///  to the next optic in the optical_system.
    ///  This optic is not applied to the wavefront.
    ///  If we're at the end of the optical_system,
    ///  this function returns without doing anything.
    void step_forward(wavefront & wf) const;
  
    ///////////////////////////////////////////
    ///  Propagate the wavefront backward through the 
    ///  optical_system to the previous optic.
    ///  The effect of this optic is removed from the wavefront
    ///  If we're at the beginning of the optical_system,
    ///  this function returns without doing anything.
    void step_backward(wavefront & wf) const;
  
    ///////////////////////////////////////////
    ///  Transforms the wavefront using the current 
    ///  optic.  If the wavefront has already been
    ///  transformed by the optic, this function 
    ///  returns without doing anything.
    void apply_current_optical_transformation(wavefront & wf) const;

    ///////////////////////////////////////////
    ///  Remove the transformation affected by
    ///  the current optic from the wavefront.
    ///  If this transformation has already been
    ///  removed, this function returns without 
    ///  doing anything.
    void remove_current_optical_transformation(wavefront & wf) const;

    ///////////////////////////////////////////
    ///  Return a copy of the current optic
    optic current_optic() const;
  

  };

}

#endif
