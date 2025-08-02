/*
scitao - A toolbox based on Scilab/Scicos environment for the simulation of 
wave optics, especially for the simulation of adaptive optics .

Copyright (c) 2000-2006 IAPCM, Beijing, China.  Written by
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

/* 
    Based on the software LightPipes, there are some modifications to accord 
with the Scilab/Scicos environment. the author thanks Dr. Gleb Vdovin 
for developing the excellent software package LightPipes.
*/

#ifndef LIGHTPIPES_PIPES_H
#define LIGHTPIPES_PIPES_H
    
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h> 

#define Pi 3.1415926535897932384626
#define Pi2 2*Pi

/* The structure FIELD contains the characteristics
of the light beam: number of points along the side 
of a  square grid, wavelength and side length of the
square grid, then two huge arrays of Re and Im data 
*/

typedef struct{ 

    int number;
    double size,lambda;
    int int1,int2,int3;
    double double1,double2,double3;

    int number2,real_imag,interleaved;
	double x,y,z,xx,xy,xz,yx,yy,yz,zx,zy,zz;
	double pixscale,curvature,timestamp;

    double *real, *imaginary; 

} FIELD;

double phase( double y,double x );  
void amp_phase_field( FIELD * ); 
void real_imaginary_field( FIELD * ); 

FIELD create_field( int number, int number2 );
void field_header_array( FIELD field, double *outptr );
void field_array( FIELD field, double *outptr );

FIELD array2field( double *inptr );
void header_array_field( FIELD *field, double *outptr );
double *field2array( FIELD field );

#endif
