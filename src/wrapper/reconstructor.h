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

#ifndef WRAP_RECONSTRUCTOR_H
#define WRAP_RECONSTRUCTOR_H

#include "arroyo_wrap.h"

typedef struct{ 
	APERTURE AP;
	LensletArray LenArray;
	ZERNIKE projected_modes;
	DeformableMirror DM;
	PixelArray reconstructor_matrix;
	double subaperture_illumination_threshold;
	double eigenvalue_threshold;
} RECONSTRUCTOR;

#ifdef __cplusplus

#include "../arroyo/arroyo.h"
using namespace Arroyo;

zernike_projected_zonal_reconstructor
*c2cpp_RECONSTRUCTOR( RECONSTRUCTOR Reconstructor);

arroyo_least_squares_reconstructor<float> 
arroyo_least_squares_reconstructor_file_constructor( char *fname );

extern "C" {
#endif

	RECONSTRUCTOR create_RECONSTRUCTOR 
		( APERTURE AP, DeformableMirror DM, 
		LensletArray LenArray, ZERNIKE projected_modes, 
		double subaperture_illumination_threshold, 
		double eigenvalue_threshold );

	long total_Reconstructor_array_nelem(double *inptr);
	long total_Reconstructor_array_space( RECONSTRUCTOR Reconstructor );

	RECONSTRUCTOR array2RECONSTRUCTOR( double *inptr);
	double *RECONSTRUCTOR2array ( RECONSTRUCTOR Reconstructor );

	ZERNIKE arroyo_reconstruct_zernike_residuals
		( RECONSTRUCTOR Reconstructor,SHartmannCentroids SHcentroids);

	PixelArray arroyo_reconstruct_zonal_residuals
		( RECONSTRUCTOR Reconstructor,SHartmannCentroids SHcentroids);

	void write_RECONSTRUCTOR_file(RECONSTRUCTOR Reconstructor, char *fname);
	RECONSTRUCTOR RECONSTRUCTOR_file_constructor(char *fname);

#ifdef __cplusplus
}
#endif

#endif
