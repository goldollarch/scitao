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

#ifndef WRAP_SHACK_HARTMANN_H
#define WRAP_SHACK_HARTMANN_H

#include "arroyo_wrap.h"

typedef struct{ 
	int  lenslet_x_axes, lenslet_y_axes;
	double flength, lnslt_pitch;
	int  pix_per_lenslet,  pix_per_xform;
} LensletArray;

typedef struct{ 
	int  pixarr_x_axes, pixarr_y_axes;
	double *x_centroids, *y_centroids;
} SHartmannCentroids;

#ifdef __cplusplus

#include "../arroyo/arroyo.h"

using namespace Arroyo;

square_lenslet_array
c2cpp_LensletArray(const LensletArray lnslt_array);

Shack_Hartmann_centroids 
c2cpp_SHartmannCentroids(SHartmannCentroids shcentroids);

SHartmannCentroids 
cpp2c_SHartmannCentroids(Shack_Hartmann_centroids shack_hartmann_centroids);

extern "C" {
#endif

	LensletArray create_LensletArray
		( int lenslet_x_axes,int lenslet_y_axes,double flength, 
		double lnslt_pitch, int pix_per_lenslet, int pix_per_xform);

	LensletArray array2LensletArray( double *inptr);
	double *LensletArray2array ( LensletArray lnslt_array );

	WaveFront LensletArray_transform( LensletArray lnslt_array, WaveFront WF );

	SHartmannCentroids create_SHartmannCentroids
		( int pixarr_x_axes, int pixarr_y_axes, WaveFront WF );

	void SHartmannCentroids_array(SHartmannCentroids shcentroids,double *outptr);
	double *SHartmannCentroids2array( SHartmannCentroids shcentroids );
	SHartmannCentroids array2SHartmannCentroids(double *inptr);

	void write_SHartmannCentroids( SHartmannCentroids shcentroids, char *fname, double timestamp );

#ifdef __cplusplus
}
#endif

#endif
