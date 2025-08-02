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

void tt_mirror(scicos_block *block,int flag) 
{
	int *ipar,ap_type; 
	double *rpar,data1,data2,data3,data4;

	double angular_velocity;
	double x_orient=0,y_orient=0;

	static int created=0;
	double arcsecs_to_radians = 3.1415927/180./3600.;

	APERTURE AP; 
	static TipTiltMirror TTM;

	ThreeVector mirror_command_vector;
	ZERNIKE tt_commands;

	FIELD field;
	WaveFront WF_in, WF_out;

	if (flag==1) 
	{
		field=block_inptr_field(block,0);
		if(field.number) {
			WF_in = FIELD_WaveFront(field);
			WF_out =TipTiltMirror_transform(&TTM, WF_in);
			WaveFront_FIELD(WF_out,&field);
		}
		field_block_outptr(field,block,0);	
        
		field=block_inptr_field(block,1);
		if(field.number) {
			WF_in = FIELD_WaveFront(field);
			WF_out =TipTiltMirror_transform(&TTM, WF_in);
			WaveFront_FIELD(WF_out,&field);
		}
		field_block_outptr(field,block,1);

		tt_commands = block_inptr_ZERNIKE(block,2);
		if(tt_commands.order) {
			x_orient = -get_ZNK_cos_coeff(tt_commands,1,1)*arcsecs_to_radians;
			y_orient = -get_ZNK_sin_coeff(tt_commands,1,1)*arcsecs_to_radians;
		}
		mirror_command_vector=construct_ThreeVector(x_orient,y_orient,1);

		TipTiltMirror_updata(&TTM,mirror_command_vector,field.timestamp);
	}

	else if ( flag==6 && (!created) ) 
	{
		rpar=block->rpar; 
		--rpar;
		ipar=block->ipar; 
		--ipar;
		
		ap_type=ipar[1];
        data1=rpar[1];data2=rpar[2];
		data3=rpar[3];data4=rpar[4];
		
		angular_velocity=rpar[5];
	
		if(scao.global_AP_created==1) AP=scao.global_AP;
		else AP=create_APERTURE(ap_type,data1,data2,data3,data4);

		TTM=create_TipTiltMirror(AP,angular_velocity);

		created=1;
	}
}
