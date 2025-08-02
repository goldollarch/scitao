/*
scitao - A toolbox based on Scilab/Scicos environment for the simulation of 
wave optics, especially for the simulation of adaptive optics .

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

void file_fits(scicos_block *block,int flag)
{
	double t;
	int *ipar;

	char str[50];
	int type,job = 1;

	FIELD field;
	WaveFront WF;
	SHartmannCentroids shcentroids;

	ZERNIKE ZNK;
	PixelArray PixArr;

	ipar=block->ipar;
	--ipar;

	type=ipar[1];

	F2C(cvstr)(&(ipar[2]),&(ipar[3]), str,&job,strlen(str));
	str[ipar[2]] = '\0';

	t=get_scicos_time();

	if (flag==1) 
	{
		if(type==0) {
			field=block_inptr_field(block,0);
			if(field.number) {
				WF = FIELD_WaveFront(field);
				write_WaveFront_file(WF,str,field.timestamp);
				field_block_outptr(field,block,0);
			}
		}

		else if(type==1) {
			shcentroids=block_inptr_SHartmannCentroids(block,0);
			if(shcentroids.pixarr_x_axes && shcentroids.pixarr_y_axes) 
			{
				write_SHartmannCentroids(shcentroids,str,t);
				SHartmannCentroids_block_outptr(shcentroids,block,0);
			}
		}

		else if(type==2) {
			ZNK=block_inptr_ZERNIKE(block,0);
			if(ZNK.order) {
				zernike_write(ZNK,str,t);
				ZERNIKE_block_outptr(ZNK,block,0);
			}
		}

		else if(type==3) {
			PixArr=block_inptr_PixelArray(block,0);
			if(PixArr.x_axes && PixArr.y_axes) {
				write_PixelArray_file(PixArr,str,t);
				PixelArray_block_outptr(PixArr,block,0);
			}
		}

	} 

}
