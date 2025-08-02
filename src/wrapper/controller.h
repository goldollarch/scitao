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

#ifndef WRAP_CONTROLLER_H
#define WRAP_CONTROLLER_H

#include "arroyo_wrap.h"  

typedef struct{ 
	ZERNIKE last_input;
	double propgain, intgain;
} TTM_PI_controller;

typedef struct{ 
	PixelArray last_input;
	double propgain, intgain;
} DM_PI_controller;

#ifdef __cplusplus

#include "../arroyo/arroyo.h"

using namespace Arroyo;

proportional_integral_controller<zernike,zernike, double, double> 
c2cpp_TTM_PI_controller(TTM_PI_controller ttm_pi_controller);

proportional_integral_controller<pixel_array<double>, pixel_array<double>, double, double> 
c2cpp_DM_PI_controller(DM_PI_controller dm_pi_controller);

extern "C" {
#endif

	TTM_PI_controller create_TTM_PI_controller(	double pgain, double igain );
	DM_PI_controller create_DM_PI_controller( DeformableMirror dm, double pgain, double igain );

	TTM_PI_controller array2TTM_PI_controller (double *inptr);
	double *TTM_PI_controller2array ( TTM_PI_controller ttm_controller );

	DM_PI_controller array2DM_PI_controller (double *inptr);
	double *DM_PI_controller2array ( DM_PI_controller dm_controller );

	ZERNIKE update_TTM_PI_controller( TTM_PI_controller *ttm_controller, ZERNIKE input,ZERNIKE TTM_commands );
	PixelArray update_DM_PI_controller( DM_PI_controller *dm_controller, PixelArray input,PixelArray DM_commands );

#ifdef __cplusplus
}
#endif

#endif
