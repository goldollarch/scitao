/*
scitao - A toolbox based on Scilab/Scicos environment for the simulation  
of wave optics, especially for the simulation of adaptive optics .

Copyright (c) 2005-2006 IAPCM, Beijing, China.  Written by
Chen jingyuan.  For comments or questions about this software,
please contact the author at jingyuan_chen@yahoo.com.cn.

This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as  published by the
Free Software Foundation; either version 2 of the License, or (at your
option) any later version.

This program is provided "as is" and distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, 
Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
*/

#ifndef WRAP_DEFORMABLE_MIRROR_H
#define WRAP_DEFORMABLE_MIRROR_H

#include "arroyo_wrap.h"

typedef struct{ 
	APERTURE AP;
	int actuator_x_dim, actuator_y_dim;
	double timestamp, actuator_pitch, actuator_velocity;
	double *actuator_commands, *actuator_positions;
} DeformableMirror;

#ifdef __cplusplus

#include "../arroyo/arroyo.h"

using namespace Arroyo;

extern "C" {
#endif

	DeformableMirror create_DeformableMirror( APERTURE AP, int actuator_x_dim, 
		int actuator_y_dim, double actuator_pitch, double actuator_velocity );

	void set_dm_timestamp( DeformableMirror* dm,double timestamp );
	void set_actuator_positions( DeformableMirror* dm,double *actuator_positions);
	void set_actuator_commands( DeformableMirror* dm, double *actuator_commands );

	WaveFront DeformableMirror_transform( DeformableMirror* dm, WaveFront WF );
	void DeformableMirror_updata( DeformableMirror* dm, double *actuator_commands, double timestamp );

	DeformableMirror array2DeformableMirror( double *inptr );
	double *DeformableMirror2array( DeformableMirror dm );

#ifdef __cplusplus
}
#endif

#endif
