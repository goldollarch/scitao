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

#ifndef WRAP_GEOMETRY_H
#define WRAP_GEOMETRY_H

#include "arroyo_wrap.h"

typedef struct{ 
	double x, y, z;
} ThreePoint;

typedef struct{ 
	double x, y, z;
} ThreeVector;

typedef struct{ 
	ThreePoint TP;
	ThreeVector TVx, TVy, TVz;
} ThreeFrame;

#ifdef __cplusplus

#include "../arroyo/arroyo.h"

using namespace Arroyo;

three_point c2cpp_ThreePoint(ThreePoint TP);
three_vector c2cpp_ThreeVector(ThreeVector TV);
three_frame c2cpp_ThreeFrame(ThreeFrame TF);

ThreePoint cpp2c_ThreePoint(three_point tp);
ThreeVector cpp2c_ThreeVector(three_vector tv);
ThreeFrame cpp2c_ThreeFrame(three_frame tf);

three_frame FIELD_three_frame( FIELD field ); 
void three_frame_FIELD( three_frame tf, FIELD *field ); 

extern "C" {
#endif

	ThreeFrame default_ThreeFrame();

	ThreePoint construct_ThreePoint( double x,double y,double z );
	ThreeVector construct_ThreeVector( double x,double y,double z );
	ThreeFrame construct_ThreeFrame ( ThreePoint TP,
		ThreeVector TVx, ThreeVector TVy, ThreeVector TVz );

	void set_FIELD_ThreeFrame( ThreeFrame TF, FIELD *field );

	ThreeFrame FIELD_ThreeFrame( FIELD field );

	ThreeFrame ThreeTranslaton ( ThreeFrame TF, ThreeVector TV );
	ThreeFrame ThreeReflection ( ThreeFrame TF, ThreePoint TP, ThreeVector TV );
	ThreeFrame ThreeRotation ( ThreeFrame TF, ThreePoint TP, ThreeVector TV, double angle );

	double *ThreeFrame2array( ThreeFrame TF );
	ThreeFrame array2ThreeFrame( double *inptr );

#ifdef __cplusplus
}
#endif

#endif

 
