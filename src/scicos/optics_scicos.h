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

#ifndef OPTICS_SCICOS_H
#define OPTICS_SCICOS_H

#include "machine.h"
#include "scicos_block.h"

#include "../lightPipes/lightPipes.h"
#include "../wrapper/arroyo_wrap.h"

extern int C2F(cvstr)();
extern void sciprint __PARAMS((char *fmt,...));

typedef struct{ 
	int max_number,num_Rlayers,	
		reconstructor_nsubaps,global_AP_created,
		sensing_emt_created, detected_emt_created,
		atm_mod_ref_lay_created,global_DM_created,
		global_LA_created,arroyo_Reconstructor_created;
	double sensing_wavelength,detected_wavelengths,
		sensing_pixel_scale,detected_pixel_scale;
	APERTURE global_AP;
	RefractiveLayer *Rlayers;
	DeformableMirror global_DM;
	LensletArray global_LenArray;
	RECONSTRUCTOR arroyo_Reconstructor;
	Emitter sensing_emt,detected_emt;
	WavefrontHeader sensing_WfH,detected_WfH;
} SCAO_MODEL;

SCAO_MODEL scao;

FIELD block_inptr_field( scicos_block *block, int port );
void field_block_outptr( FIELD field, scicos_block *block, int port );

SHartmannCentroids block_inptr_SHartmannCentroids( scicos_block *block, int port );
void SHartmannCentroids_block_outptr( SHartmannCentroids shcentroids, scicos_block *block, int port );

PixelArray block_inptr_PixelArray( scicos_block *block, int port );
void PixelArray_block_outptr( PixelArray PixArr, scicos_block *block, int port );

ZERNIKE block_inptr_ZERNIKE( scicos_block *block, int port );
void ZERNIKE_block_outptr( ZERNIKE ZNK, scicos_block *block, int port );

#endif
