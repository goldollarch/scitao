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

#ifndef WRAP_SOME_UTILS_H
#define WRAP_SOME_UTILS_H

#include "arroyo_wrap.h"

typedef struct{ 
	ThreeFrame TF;
	int foreshortening;
} Optic;

typedef struct{ 
	int wt_alloc; 
	int x_axes, y_axes;
	double *pixeldata;
	float *pixelwts;
} PixelArray;
//now wt_alloc must be 0,ie. false;

#ifdef __cplusplus

#include "../arroyo/arroyo.h"
using namespace Arroyo;

PixelArray cpp2c_PixelArray(pixel_array<double> pixarr);
pixel_array<double> c2cpp_PixelArray(PixelArray PixArr); 

char *Current_filename(char *fname,double timestamp,char *type);

extern "C" {
#endif

	Optic create_Optic(int foreshortening,ThreeFrame TF);

	PixelArray create_PixelArray(int x_axes, int y_axes, int wt_alloc);
	void set_PixelArray_data(PixelArray *PixArr,double *data);

	PixelArray array2PixelArray( double *inptr);
	double *PixelArray2array ( PixelArray PixArr );
	void PixelArray_array(PixelArray PixArr,double *outptr);

	char *simple_filename(char *fname,char *type);
	char *number_filename(char *fname,int n, char *type);
	char *current_filename(char *fname,double timestamp,char *type);

#ifdef __cplusplus
}
#endif

#endif
