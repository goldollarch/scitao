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

#ifndef WRAP_ZERNIKE_H
#define WRAP_ZERNIKE_H

#include "arroyo_wrap.h"

typedef struct{ 
	int order;
	double *coefficients;
} ZERNIKE;

#ifdef __cplusplus

#include "../arroyo/arroyo.h"

using namespace Arroyo;

ZERNIKE cpp2c_zernike(zernike znk);
zernike c2cpp_zernike(ZERNIKE ZNK);

extern "C" {
#endif

	int total_znk_space( int order );
	ZERNIKE create_ZERNIKE(int order);

	void set_ZNK_cos_coef(ZERNIKE *ZNK,int order, int level, double coeff);
	void set_ZNK_sin_coef(ZERNIKE *ZNK,int order, int level, double coeff);

	double get_ZNK_cos_coeff (ZERNIKE ZNK, int order, int level) ;
	double get_ZNK_sin_coeff (ZERNIKE ZNK, int order, int level) ;

	void zernike_write( ZERNIKE ZNK, char *fname, double timestamp );

	void ZERNIKE_array(ZERNIKE ZNK,double *outptr);
	double *ZERNIKE2array( ZERNIKE ZNK );
	ZERNIKE array2ZERNIKE( double *inptr );

#ifdef __cplusplus
}
#endif

#endif
