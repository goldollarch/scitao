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

#include "optics_scicos.h"

void controller(scicos_block *block,int flag) 
{
	int *ipar; 
	double *rpar;

	int type,nx_actuator,ny_actuator;
	double proportional_gain ;
	double integral_gain;
	static int tt_pi_created=0;
	static int dm_pi_created=0;

	static TTM_PI_controller ttm_controller;
	static ZERNIKE tt_commands;
	ZERNIKE tt_residuals;

	static DM_PI_controller dm_controller;
	static PixelArray dm_commands;
	PixelArray dm_residuals;

	rpar=block->rpar; 
	--rpar;
	ipar=block->ipar; 
	--ipar;

	type=ipar[1];

	if (flag==1) {

		if(type==0) {

			tt_residuals = block_inptr_ZERNIKE( block, 0 );

			if (tt_residuals.order) {
				tt_commands=update_TTM_PI_controller(&ttm_controller,
					tt_residuals,tt_commands);
                ZERNIKE_block_outptr(tt_commands,block,0);
			}

		}

		else if(type==1) {

			dm_residuals =block_inptr_PixelArray( block, 0 );

			if(dm_residuals.x_axes && dm_residuals.y_axes) {
				dm_commands=update_DM_PI_controller(&dm_controller,
					dm_residuals,dm_commands);
				PixelArray_block_outptr(dm_commands,block,0);
			}

		}

	}

	else if ( flag==6 ) {
		
		proportional_gain=rpar[1];
		integral_gain=rpar[2];

		if( type==0 && (!tt_pi_created) ) {

			tt_pi_created=1;
			tt_commands=create_ZERNIKE(1);
			ttm_controller=	create_TTM_PI_controller
				(proportional_gain,integral_gain);
		
		}

		else if( type==1 && (!dm_pi_created) && (scao.global_DM_created==1)) {

			dm_pi_created=1;

			nx_actuator=scao.global_DM.actuator_x_dim;
			ny_actuator=scao.global_DM.actuator_y_dim;
			dm_commands=create_PixelArray(nx_actuator,ny_actuator,0);

			dm_controller=create_DM_PI_controller
				(scao.global_DM, proportional_gain,integral_gain);
		
		}

	}
}
