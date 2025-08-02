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

#ifndef WRAP_WAVEFRONT_H
#define WRAP_WAVEFRONT_H

#include "arroyo_wrap.h"

typedef struct{ 
	ThreeFrame TF;
	int axes_x,axes_y;
	double lambda,pixscale,curvature,timestamp;
} WavefrontHeader;

typedef struct{ 
	WavefrontHeader WfH;
    int real_imag, interleaved;
    double *real, *imaginary; 
} WaveFront;

#ifdef __cplusplus

#include "../arroyo/arroyo.h"

using namespace Arroyo;

diffractive_wavefront_header<double> c2cpp_WavefrontHeader(WavefrontHeader WfH); 
WavefrontHeader cpp2c_WavefrontHeader(diffractive_wavefront_header<double> dwfh);

WaveFront cpp2c_WaveFront(diffractive_wavefront<double> dwf);
diffractive_wavefront<double> c2cpp_WaveFront( WaveFront WF);

extern "C" {
#endif

	WavefrontHeader construct_WavefrontHeader( ThreeFrame TF, 
		int axes_x, int axes_y, double lambda, double pixscale );

	void set_FIELD_WavefrontHeader( unsigned int number, 
		unsigned int number2,double lambda, double pixscale, 
		double curvature, double timestamp, FIELD *field );

	WaveFront shift_WaveFront ( WaveFront WF, double dx, double dy );
	WaveFront rotate_WaveFront ( WaveFront WF, double angle);
	WaveFront reflect_WaveFront ( WaveFront WF, double x, double y, double z );

	WaveFront WaveFront_clip_array( WaveFront WF, int nclip );
	WaveFront set_dwf_pixel_scale( WaveFront WF, double pxlscl);
	WaveFront set_WaveFront_propagation_direction( WaveFront,ThreeVector);

	WaveFront FIELD_WaveFront( FIELD field );
	WavefrontHeader FIELD_WavefrontHeader( FIELD field);
	void WaveFront_FIELD( WaveFront WF, FIELD *field );

	void write_WaveFront_file(WaveFront WF,char *fname,double timestamp);

	double *WavefrontHeader2array( WavefrontHeader WFH );
	WavefrontHeader array2WavefrontHeader( double *inptr );

#ifdef __cplusplus
}
#endif

#endif
