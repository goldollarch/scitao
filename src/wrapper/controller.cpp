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

#include "arroyo_wrap.h"
using namespace Arroyo;

TTM_PI_controller create_TTM_PI_controller( double pgain, double igain )
{
	TTM_PI_controller ttm_controller;

	ttm_controller.propgain=pgain;
	ttm_controller.intgain=igain;

	ZERNIKE ZNK=create_ZERNIKE(1);
	set_ZNK_cos_coef(&ZNK,1,1,0);
	set_ZNK_sin_coef(&ZNK,1,1,0);

	ttm_controller.last_input=ZNK;

	return ttm_controller;
}

DM_PI_controller create_DM_PI_controller
( DeformableMirror dm, double pgain, double igain )
{
	DM_PI_controller dm_controller;

	dm_controller.propgain=pgain;	
	dm_controller.intgain=igain;

	dm_controller.last_input=create_PixelArray
		(dm.actuator_x_dim,dm.actuator_y_dim,0);

	return dm_controller;
}

proportional_integral_controller<zernike,zernike, double, double> 
c2cpp_TTM_PI_controller(TTM_PI_controller ttm_pi)
{
	zernike znke=c2cpp_zernike(ttm_pi.last_input);
	proportional_integral_controller<zernike, zernike, double, double> 
		TT_controller( znke, ttm_pi.propgain,ttm_pi.intgain);
	
	return TT_controller;
}

proportional_integral_controller<pixel_array<double>, pixel_array<double>, double, double> 
c2cpp_DM_PI_controller(DM_PI_controller dm_pi)
{
	pixel_array<double> parray
		=c2cpp_PixelArray(dm_pi.last_input);

	proportional_integral_controller<pixel_array<double>, pixel_array<double>, double, double> 
		DM_controller( parray, dm_pi.propgain, dm_pi.intgain );
	
	return DM_controller;
}

ZERNIKE update_TTM_PI_controller
( TTM_PI_controller *ttm_controller, ZERNIKE input,ZERNIKE TTM_commands )
{
	zernike tt_residuals=c2cpp_zernike(input);
	proportional_integral_controller<zernike, zernike, double, double> 
		TT_controller=c2cpp_TTM_PI_controller(*ttm_controller);

	zernike tt_commands=c2cpp_zernike(TTM_commands);
	TT_controller.update(tt_residuals,tt_commands);
	ttm_controller->last_input=input;

	return (cpp2c_zernike(tt_commands));
}

PixelArray update_DM_PI_controller
( DM_PI_controller *dm_controller, PixelArray input, PixelArray DM_commands )
{
	pixel_array<double> dm_residuals=c2cpp_PixelArray(input);
	proportional_integral_controller<pixel_array<double>, pixel_array<double>, double, double> 
		DM_controller=c2cpp_DM_PI_controller(*dm_controller);

	pixel_array<double> dm_commands=c2cpp_PixelArray(DM_commands);
	DM_controller.update(dm_residuals,dm_commands);
	dm_controller->last_input=input;

	return(cpp2c_PixelArray(dm_commands));
}

double *TTM_PI_controller2array ( TTM_PI_controller ttm_controller )
{
	int i, nelem;
	double *tmp=NULL;

	nelem=total_znk_space( ttm_controller.last_input.order );
	tmp= (double *)	calloc( nelem+3, sizeof(double) );

	tmp[0]=ttm_controller.propgain; tmp[1]=ttm_controller.intgain;
	tmp[2]=ttm_controller.last_input.order;
	for(i=0;i<nelem;i++) tmp[i+3]=ttm_controller.last_input.coefficients[i];

	return tmp;
}

TTM_PI_controller array2TTM_PI_controller (double *inptr)
{
	TTM_PI_controller ttm_controller;

	ttm_controller.propgain=inptr[0]; ttm_controller.intgain=inptr[1];
	ttm_controller.last_input=array2ZERNIKE(inptr+2);

	return ttm_controller;
}

double *DM_PI_controller2array ( DM_PI_controller dm_controller )
{
	int i, nelem;
	double *tmp=NULL;

	nelem=dm_controller.last_input.x_axes*dm_controller.last_input.y_axes;
	tmp= (double *)	calloc( nelem+4, sizeof(double) );

	tmp[0]=dm_controller.last_input.x_axes; 
	tmp[1]=dm_controller.last_input.y_axes;
	tmp[2]=dm_controller.propgain; tmp[3]=dm_controller.intgain;
	for(i=0;i<nelem;i++) 
		tmp[i+4]=dm_controller.last_input.pixeldata[i];

	return tmp;
}

DM_PI_controller array2DM_PI_controller (double *inptr)
{
	int x_axes, y_axes;
	DM_PI_controller dm_controller;

	x_axes=(int)inptr[0]; y_axes=(int)inptr[1];
	dm_controller.propgain=inptr[2]; 
	dm_controller.intgain=inptr[3];

	dm_controller.last_input=create_PixelArray(x_axes,y_axes,0);
	dm_controller.last_input.pixeldata=inptr+4;

	return dm_controller;
}
