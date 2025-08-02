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

extern int C2F(sciwin)(), C2F(dr1)(), C2F(xgray)();

void display(scicos_block *block,int flag)
{
	double t;
	double *z__;
	int nipar, *ipar;
 
	/* Initialized data */

	static int cur = 0;
	static int verb = 0;
	
	int i,type;

	static int v;
	static int na;
	static double dv;
	static int wid, iwd;
	static int iwp;

	FIELD field;
	PixelArray PixArr;
	SHartmannCentroids shcentroids;

	/*     ipar(1) = win_num */
	/*     ipar(2:3) : window position */
	/*     ipar(4:5) : window dimension */
	
	double *z;
	int n_x, n_y;

	char strflag[3]="021";
	int aaint[4]={10,10,10,10};
	static double brect[4];
	long int l1=0;

	ipar=block->ipar;
	nipar=block->nipar;

	t=get_scicos_time();
  
	--ipar;

	iwp = nipar - 3;
	iwd = nipar - 1;

	type=ipar[1];
	wid = ipar[2];

	if(wid==-1)  
		wid=10+get_block_number();
 
	if (flag == 2) {

		C2F(dr1)("xget\000", "window\000", &verb, &cur, &na, &v, &v, &v, &dv, &dv, &dv, &dv);
		if (cur != wid) C2F(dr1)("xset\000", "window\000", &wid, &v, &v, &v, &v, &v, &dv, &dv, &dv, &dv);

		C2F(dr1)("xclear\000", "v\000", &v, &v, &v, &v, &v, &v, &dv, &dv, &dv, &dv);

		z__=*block->work; --z__;
		for(i=1;i<=scao.max_number;i++) z__[i]=i;
		++z__;

		if(type==4) {
			PixArr=block_inptr_PixelArray(block,0);
			n_x = PixArr.x_axes; n_y = PixArr.y_axes;
			z = PixArr.pixeldata;
		}
		else if(type==5) {
			shcentroids=block_inptr_SHartmannCentroids(block,0);
			n_x=shcentroids.pixarr_x_axes; n_y=shcentroids.pixarr_y_axes;
			z=shcentroids.x_centroids;
		}

		else {
			field=block_inptr_field(block,0);
			n_x=field.number; n_y=field.number2;

			if(type==1) z = field.real;
			else if(type==2) z=field.imaginary;
			else {
				amp_phase_field( &field );
				if(type==3) z=field.imaginary;
				else z=field.real;
			}
		}

		C2F(xgray)(  z__, z__, z, &n_x, &n_y, strflag, &brect, &aaint, 0, l1);

	} 

	else if (flag == 4) {
		if ((*block->work = scicos_malloc(sizeof(double)*(scao.max_number)))== NULL ) {
			set_block_error(-16);
			return;
		}

		C2F(sciwin)();
		C2F(dr1)("xset\000", "window\000", &wid, &v, &v, &v, &v, &v, &dv, &dv, &dv, &dv);

		if (ipar[iwp] >= 0) C2F(dr1)("xset\000", "wpos\000", &ipar[iwp], &ipar[iwp + 1], &v, &v, &v, &v, &dv, &dv, &dv, &dv);
 		if (ipar[iwd] >= 0) {
			C2F(dr1)("xset\000", "wdim\000", &ipar[iwd], &ipar[iwd + 1], &v, &v, &v, &v, &dv, &dv, &dv, &dv);
			C2F(dr1)("xset\000", "window\000", &wid, &v, &v, &v, &v, &v, &dv, &dv, &dv, &dv);
		}

	} 

	else if (flag == 5) {
		scicos_free(*block->work);
	}
}
