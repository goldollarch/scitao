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

#include "fftn.h"
#include "pipes.h"

double phase( double y,double x )
{
	double pp=0.; 

    if(x==0.){
       if (y>0) pp=0.5*Pi;
       if (y==0) pp=0.;
       if (y<0) pp=-0.5*Pi;
     }	else {
		if(y!=0) pp=atan2(y,x);
		else pp=0.;
     }
	return pp;
}

void fft3( int ind, FIELD field )
{
	int dims[2];
	dims[0]=field.number; dims[1]=field.number2;
    fftn( 2, dims, field.real, field.imaginary, ind, (double) field.number );
}

void amp_phase_field( FIELD *field )
{
	int i;
	long  size;
	double tmp_real,tmp_imag;

	if(field->real_imag) {
		size=(field->number)*(field->number2);
		for(i=0;i<size;i++)  {
			tmp_real=field->real[i];
			tmp_imag=field->imaginary[i];
			field->real[i]=sqrt(tmp_real*tmp_real + tmp_imag*tmp_imag);
			field->imaginary[i]=phase(tmp_imag,tmp_real);
		}
		field->real_imag=0;
	}
}

void real_imaginary_field( FIELD *field )
{
	int i;
	long  size;
	double tmp_real,tmp_imag;

	if(!field->real_imag) {
		size=(field->number)*(field->number2);
		for(i=0;i<size;i++)  {
			tmp_real=field->real[i];
			tmp_imag=field->imaginary[i];
			field->real[i]=tmp_real*cos(tmp_imag);
			field->imaginary[i]=tmp_real*sin(tmp_imag);
		}
		field->real_imag=1;
	}
}

FIELD create_field( int number, int number2 )
{
 	FIELD field;

    field.number=number;
	field.number2=number2;

	field.real = (double *) calloc( (number)*(number2), sizeof(double) );
	field.imaginary = (double *) calloc( (number)*(number2), sizeof(double) );

	return field;
}

void header_array_field( FIELD *field, double *inptr )
{
	int number,number2;

	number=(int)inptr[0];
	number2=(int)inptr[10];

	field->number=number;
	field->size=inptr[1];field->lambda=inptr[2];
	field->int1=(int)inptr[3];field->int2=(int)inptr[4];field->int3=(int)inptr[5];
    field->double1=inptr[6];field->double2=inptr[7];field->double3=inptr[8];

	// inptr[9] not used 

	field->number2=number2;	
	field->real_imag=(int)inptr[11];field->interleaved=(int)inptr[12];

    field->x=inptr[13];field->y=inptr[14];field->z=inptr[15];
    field->xx=inptr[16];field->xy=inptr[17];field->xz=inptr[18];
    field->yx=inptr[19];field->yy=inptr[20];field->yz=inptr[21];
    field->zx=inptr[22];field->zy=inptr[23];field->zz=inptr[24];

    field->pixscale=inptr[25];
	field->curvature=inptr[26];	
	field->timestamp=inptr[27];

	// inptr[28] to inptr[29] not used 
}

FIELD array2field( double *inptr )
{
	FIELD field;
    unsigned long  size;
	int number,number2;

	number=(int)inptr[0];
	number2=(int)inptr[10];
	size=number*number2;

	field = create_field( number, number2 );

	header_array_field(&field,inptr);

    field.real=inptr + 30;
	field.imaginary=field.real + size;

	return field;
}

void field_header_array( FIELD field, double *outptr )
{
	outptr[0]=field.number;
	outptr[1]=field.size;outptr[2]=field.lambda;

    outptr[3]=field.int1; outptr[4]=field.int2;	outptr[5]=field.int3;
    outptr[6]=field.double1;outptr[7]=field.double2;outptr[8]=field.double3;

	outptr[10]=field.number2;
	outptr[11]=field.real_imag;outptr[12]=field.interleaved;

    outptr[13]=field.x;outptr[14]=field.y;outptr[15]=field.z;
    outptr[16]=field.xx;outptr[17]=field.xy;outptr[18]=field.xz;
    outptr[19]=field.yx;outptr[20]=field.yy;outptr[21]=field.yz;
    outptr[22]=field.zx;outptr[23]=field.zy;outptr[24]=field.zz;

    outptr[25]=field.pixscale;outptr[26]=field.curvature;
	outptr[27]=field.timestamp;
}

void field_array( FIELD field, double *outptr )
{
	unsigned long  size; 
    int i,j,number,number2;

	number=field.number; number2=field.number2;
	size=number*number2;

	field_header_array(field,outptr);

	for (i=0;i<field.number;i++) {
		for (j=0;j<field.number2 ;j++) {
			outptr[i*field.number2+j+30] =
				field.real[i*field.number2+j];
			outptr[i*field.number2+j+30+size] =
				field.imaginary[i*field.number2+j];
		}
	}
}

double *field2array( FIELD field )
{
	double *outptr;
    int number,number2;
	unsigned long  size; 

	number=field.number; 
	number2=field.number2;
	size=number*number2;

	outptr = (double*)calloc( 2*size+30, sizeof(double));

	field_array(field,outptr);

	return outptr;
}
