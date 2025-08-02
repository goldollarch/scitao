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

#ifndef WRAP_WIND_MODEL_H
#define WRAP_WIND_MODEL_H

#include "arroyo_wrap.h"

typedef struct{ 
	double  ground_layer_wind_velocity;
	double  tropopause_wind_velocity;
	double  tropopause_height;
	double  tropopause_thickness;
} HardyWind;

#ifdef __cplusplus

#include "../arroyo/arroyo.h"
using namespace Arroyo;

Hardy_wind_model c2cpp_HardyWind( HardyWind HWM);
// HardyWind cpp2c_HardyWind( Hardy_wind_model hwm);

extern "C" {
#endif

	HardyWind construct_HardyWind(double ground_vel,double tropopause_vel,
		double tropopause_height, double tropopause_thickness);

	ThreeVector *get_HardyWind_velocities( HardyWind HWM,
		       int nlayers, double *height, ThreeFrame TF );

#ifdef __cplusplus
}
#endif

#endif
