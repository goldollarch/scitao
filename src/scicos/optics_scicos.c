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
with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
*/

#include "optics_scicos.h"

FIELD block_inptr_field( scicos_block *block, int port )
{
	return array2field(block->inptr[port]);
}

void field_block_outptr( FIELD field, scicos_block *block, int port )
{
	field_array (field, block->outptr[port]);
}

SHartmannCentroids 
block_inptr_SHartmannCentroids( scicos_block *block, int port )
{
	return array2SHartmannCentroids(block->inptr[port]);
}

void SHartmannCentroids_block_outptr
( SHartmannCentroids shcentroids, scicos_block *block, int port )
{
	SHartmannCentroids_array(shcentroids, block->outptr[port]);
}

PixelArray block_inptr_PixelArray( scicos_block *block, int port )
{
	return array2PixelArray(block->inptr[port]);
}

void PixelArray_block_outptr( PixelArray PixArr, scicos_block *block, int port )
{
	PixelArray_array(PixArr, block->outptr[port]);
}

ZERNIKE block_inptr_ZERNIKE( scicos_block *block, int port )
{
	return array2ZERNIKE(block->inptr[port]);
}

void ZERNIKE_block_outptr( ZERNIKE ZNK, scicos_block *block, int port )
{
	ZERNIKE_array(ZNK, block->outptr[port]);
}
