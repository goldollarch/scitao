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

void reconstructor(scicos_block *block,int flag) 
{
	char str[50];
	int *ipar,job = 1; 
	double *rpar;

	static int steps=0;
	int nsubaps,ttm_delay=1,dm_delay=1;

	double data1,eigenvalue_threshold;
	double subaperture_illumination_threshold;

	APERTURE tmp_AP; 
	DeformableMirror tmp_DM;
	LensletArray tmp_lnslt_array;
	SHartmannCentroids shcentroid;

	PixelArray dm_residuals;
	ZERNIKE projection_modes,tip_tilt_residuals;

	rpar=block->rpar; 
	--rpar;
	ipar=block->ipar; 
	--ipar;

	nsubaps=ipar[1];
	
	F2C(cvstr)(&(ipar[2]),&(ipar[3]), str,&job,strlen(str));
	str[ipar[2]] = '\0';

	if (flag==1) {

		shcentroid=block_inptr_SHartmannCentroids(block,0);

		if(shcentroid.pixarr_x_axes && shcentroid.pixarr_y_axes ) 
		{
			steps++;
			if(steps>dm_delay) {
				dm_residuals=arroyo_reconstruct_zonal_residuals
					( scao.arroyo_Reconstructor, shcentroid );
				PixelArray_block_outptr(dm_residuals,block,0);
			}
			if(steps>ttm_delay) {
				tip_tilt_residuals = arroyo_reconstruct_zernike_residuals
					( scao.arroyo_Reconstructor, shcentroid );
				ZERNIKE_block_outptr(tip_tilt_residuals,block,1);
			}
		}
	}

	else if (flag == 4)  
		if(nsubaps) 
			scao.reconstructor_nsubaps=nsubaps;
		else {
			scao.arroyo_Reconstructor =
				RECONSTRUCTOR_file_constructor(simple_filename(str,".fits"));

			scao.reconstructor_nsubaps=
				scao.arroyo_Reconstructor.LenArray.lenslet_x_axes;
			scao.arroyo_Reconstructor_created=1;
		}

	else if (flag == 5)  
		scao.arroyo_Reconstructor_created=0;

	else if ( flag==6 && (!scao.arroyo_Reconstructor_created) ) 
	{
		eigenvalue_threshold=rpar[1];
		subaperture_illumination_threshold=rpar[2];

		if( (scao.global_AP_created==1) 
			&& (scao.global_LA_created==1) 
			&& (scao.global_DM_created==1) ) 
		{
			projection_modes=create_ZERNIKE(1);
			set_ZNK_cos_coef(&projection_modes,1,1,1);
			set_ZNK_sin_coef(&projection_modes,1,1,1);

			tmp_AP=scao.global_AP;
			data1=scao.global_AP.ap_data1;
			if(scao.global_AP.aperture_type==2) {
				data1=sqrt(scao.global_AP.ap_data1*scao.global_AP.ap_data1+
					scao.global_AP.ap_data2*scao.global_AP.ap_data2);
				tmp_AP=create_APERTURE(0,data1,0,0,0);
			}
			else if(scao.global_AP.aperture_type==3) {
				data1=2*scao.global_AP.ap_data1;
				tmp_AP=create_APERTURE(0,data1,0,0,0);
			}

			tmp_DM=scao.global_DM; tmp_DM.AP=tmp_AP;
			tmp_lnslt_array=scao.global_LenArray;

			scao.arroyo_Reconstructor=create_RECONSTRUCTOR
				( tmp_AP,tmp_DM,tmp_lnslt_array,projection_modes,
				subaperture_illumination_threshold,eigenvalue_threshold );

			if( (strstr(str, "void")) == NULL )
				write_RECONSTRUCTOR_file(scao.arroyo_Reconstructor,str);
            
			scao.arroyo_Reconstructor_created=1;
		}

	}

}
