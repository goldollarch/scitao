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

void df_mirror(scicos_block *block,int flag) 
{
	int i,j,index,*ipar,ap_type; 
	double *rpar,data1,data2,data3,data4;

	int actuator_x_dim, actuator_y_dim;
	double actuator_pitch, actuator_velocity;

	APERTURE dm_AP; 
	static DeformableMirror DM;
	PixelArray dm_commands;

	ThreeFrame tmp_TF;
	ThreeVector tmp_TVx,tmp_TVy,tmp_TVz;
	static ThreeVector AP_TVz;

	FIELD field;
	WaveFront  WF_in, WF_out;

	if (flag==1) 
	{
		field=block_inptr_field(block,0);
		if(field.number) {
			WF_in = FIELD_WaveFront(field);
			WF_in=set_WaveFront_propagation_direction(WF_in,AP_TVz);
			WF_out =DeformableMirror_transform(&DM, WF_in);
			WF_out=set_WaveFront_propagation_direction(WF_out,WF_in.WfH.TF.TVz);
			WaveFront_FIELD(WF_out,&field);
		}
		field_block_outptr(field,block,0);

		field=block_inptr_field(block,1);
		if(field.number) {
			WF_in = FIELD_WaveFront(field);
			WF_in=set_WaveFront_propagation_direction(WF_in,AP_TVz);
			WF_out =DeformableMirror_transform(&DM, WF_in);
			WaveFront_FIELD(WF_out,&field);
		}
		field_block_outptr(field,block,1);

		dm_commands = block_inptr_PixelArray(block,2);
		for(i=0;i<dm_commands.x_axes;i++) {
			for(j=0;j<dm_commands.y_axes;j++) {
				index=i*dm_commands.y_axes+j;
				dm_commands.pixeldata[index] *= -1;
			}
		}

		DeformableMirror_updata(&DM,dm_commands.pixeldata,field.timestamp);
	}

	else if (flag==6 && (scao.global_DM_created!=1) ) 
	{
		rpar=block->rpar; 
		--rpar;
		ipar=block->ipar; 
		--ipar;

		ap_type=ipar[1];
		actuator_x_dim=ipar[2];
		actuator_y_dim=ipar[3];

		data1=rpar[1];data2=rpar[2];
		data3=rpar[3];data4=rpar[4];

		actuator_pitch=rpar[5];
		actuator_velocity=rpar[6];

		if ( scao.reconstructor_nsubaps ) {
			actuator_x_dim=scao.reconstructor_nsubaps+1; 
			actuator_y_dim=scao.reconstructor_nsubaps+1;
		}
	
		if(scao.global_AP_created) {
			dm_AP=scao.global_AP;
			data1=scao.global_AP.ap_data1;
			if(scao.global_AP.aperture_type==2)
				data1=sqrt(scao.global_AP.ap_data1*scao.global_AP.ap_data1+\
					scao.global_AP.ap_data2*scao.global_AP.ap_data2);
			else if(scao.global_AP.aperture_type==3) 
				data1=2*scao.global_AP.ap_data1;
		}
		else dm_AP=create_APERTURE(ap_type,data1,data2,data3,data4);
        
		AP_TVz=dm_AP.TF.TVz;

		tmp_TVx=construct_ThreeVector
			(-dm_AP.TF.TVx.x,-dm_AP.TF.TVx.y,-dm_AP.TF.TVx.z);
		tmp_TVy=construct_ThreeVector
			(-dm_AP.TF.TVy.x,-dm_AP.TF.TVy.y,-dm_AP.TF.TVy.z);
		tmp_TVz=construct_ThreeVector
			(-dm_AP.TF.TVz.x,-dm_AP.TF.TVz.y,-dm_AP.TF.TVz.z);
		tmp_TF=construct_ThreeFrame
			(dm_AP.TF.TP,tmp_TVx,tmp_TVy,tmp_TVz);
		dm_AP.TF=tmp_TF;

		actuator_pitch = data1/(double)(actuator_x_dim-1);

		DM = create_DeformableMirror( dm_AP, actuator_x_dim,
			actuator_y_dim, actuator_pitch, actuator_velocity);

		scao.global_DM=DM;

		scao.global_DM_created=1;
	}

	else if (flag == 5) scao.global_DM_created=0;

}
