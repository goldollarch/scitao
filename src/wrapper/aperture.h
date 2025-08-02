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

#ifndef WRAP_APERTURE_H
#define WRAP_APERTURE_H

#include "arroyo_wrap.h"

typedef struct{ 
	ThreeFrame TF;
	int aperture_type;
	double ap_data1,ap_data2,ap_data3,ap_data4;
} APERTURE;

#ifdef __cplusplus

#include "../arroyo/arroyo.h"

using namespace Arroyo;

aperture *c2cpp_APERTURE(APERTURE AP);

extern "C" {
#endif

	APERTURE create_APERTURE( int aperture_type, double ap_data1,
		          double ap_data2, double ap_data3, double ap_data4 );

	WaveFront APERTURE_transform(APERTURE AP,WaveFront WF);

	double *APERTURE2array ( APERTURE AP );
	APERTURE array2APERTURE (double *inptr);

#ifdef __cplusplus
}
#endif

#endif
