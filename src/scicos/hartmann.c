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

void hartmann(scicos_block *block,int flag) 
{
	char str[50];
	int *ipar,job = 1;
	double *rpar;

	double data1;
	int lenslet_x_axes, lenslet_y_axes;
    int final_wf_pix_per_lenslet,final_wf_pix_per_xform;
	double lenslet_pitch_meters,lenslet_focal_length_meters;
	
	FIELD field;
	WaveFront  WF_in, WF_out;

	static LensletArray lnslt_array;
	SHartmannCentroids shcentroid;

	rpar=block->rpar; 
	--rpar;
	ipar=block->ipar; 
	--ipar;

	lenslet_x_axes=ipar[1];
	lenslet_y_axes=ipar[2];

	final_wf_pix_per_lenslet=ipar[3];
	final_wf_pix_per_xform=ipar[4]; 
    
	F2C(cvstr)(&(ipar[5]),&(ipar[6]), str,&job,strlen(str));
	str[ipar[5]] = '\0';

	lenslet_pitch_meters=rpar[1]; 
	lenslet_focal_length_meters=rpar[2];
	
	if (flag==1) {

		field=block_inptr_field(block,0);

		if(field.number) {

			WF_in = FIELD_WaveFront(field);

			if (scao.global_AP_created)
				WF_in.WfH.TF=scao.global_AP.TF;
			if (scao.global_DM_created==1) WF_in.WfH.pixscale *= 
				lenslet_pitch_meters/scao.global_DM.actuator_pitch;

			WF_out = LensletArray_transform( lnslt_array, WF_in );
			if( (strstr(str,"void")) == NULL )
				write_WaveFront_file( WF_out,str,field.timestamp);

			WF_out = WaveFront_clip_array( WF_out, 
				(final_wf_pix_per_xform-final_wf_pix_per_lenslet)/2);
			shcentroid = create_SHartmannCentroids 
				( lenslet_x_axes, lenslet_y_axes, WF_out );
			SHartmannCentroids_block_outptr(shcentroid,block,0);
		}
		
	}

	else if ( flag==6 && (scao.global_LA_created!=1) ) {

		if ( scao.reconstructor_nsubaps ) {
			lenslet_x_axes=scao.reconstructor_nsubaps; 
			lenslet_y_axes=scao.reconstructor_nsubaps;
		}

		lnslt_array=create_LensletArray
			( lenslet_x_axes, lenslet_y_axes,
			lenslet_focal_length_meters, lenslet_pitch_meters,
			final_wf_pix_per_lenslet, final_wf_pix_per_xform);

		if (scao.global_AP_created) {
			data1=scao.global_AP.ap_data1;
			if(scao.global_AP.aperture_type==2)
				data1=sqrt(scao.global_AP.ap_data1*scao.global_AP.ap_data1+
					scao.global_AP.ap_data2*scao.global_AP.ap_data2);
			else if(scao.global_AP.aperture_type==3)
				data1=2*scao.global_AP.ap_data1;
			lenslet_pitch_meters = data1/(double)lenslet_x_axes;
		}
		
		scao.global_LenArray = create_LensletArray
			( lenslet_x_axes, lenslet_y_axes,
			lenslet_focal_length_meters,lenslet_pitch_meters,
			final_wf_pix_per_lenslet, final_wf_pix_per_xform);

		scao.global_LA_created=1;

	}

	else if (flag == 5)  
		scao.global_LA_created=0;

}
