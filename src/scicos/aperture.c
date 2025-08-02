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

void aperture(scicos_block *block,int flag)
{
	int *ipar,type; 
	double *rpar,data1,data2,data3,data4;

	FIELD field;
	APERTURE AP; 

	WaveFront WF_in,WF_out;

	rpar=block->rpar;
	--rpar;
	ipar=block->ipar;
	--ipar;

	type=ipar[1];
	data1=rpar[1]; data2=rpar[2];
	data3=rpar[3]; data4=rpar[4];
	
	AP=create_APERTURE(type,data1,data2,data3,data4);

	if(flag==1)
	{
		field=block_inptr_field(block,0);

		if(field.number) {
			WF_in = FIELD_WaveFront(field);
			WF_out =APERTURE_transform(AP, WF_in);
			WaveFront_FIELD(WF_out,&field);
		}
		
		field_block_outptr(field,block,0);

		if(block->nin==2) {

			field=block_inptr_field(block,1);

			if(field.number) {
				WF_in = FIELD_WaveFront(field);
				WF_out =APERTURE_transform(AP, WF_in);
				WaveFront_FIELD(WF_out,&field);
			}
			
			field_block_outptr(field,block,1);

		}
	}

	else if (flag==4)
	{
		if(!scao.global_AP_created) {
			scao.global_AP = AP;
			scao.global_AP_created=1;
		}
	}

	else if (flag==5) {
		scao.global_AP_created=0;
	}

	else if (flag==6) {
		if(scao.global_AP.ap_data1<data1)
			scao.global_AP = AP;
	}

}
