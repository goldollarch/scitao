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

/*--------------------------------------------------------------*/
/*      (C) Gleb Vdovin 1993-1999                               */
/*      This file is a part of LightPipes package               */
/*      Send bug reports to gleb@okotech.com                    */
/*                                                              */
/*--------------------------------------------------------------*/
    
#include "lightPipes.h"

FIELD *b_split( FIELD field, double split )
{
    int i,j;  
    double asplit, bsplit, csplit;
	long ik1;

	FIELD *outptr;

	outptr=(FIELD*) calloc( 2, sizeof(FIELD) );

    asplit=sqrt(split);
    bsplit=sqrt(1-split);
    csplit=bsplit/asplit;

	for (i=1;i<=field.number ;i++) 
		for (j=1;j<=field.number ;j++){
			ik1=(i-1)*field.number+j-1;
			field.real[ik1] *= asplit;
			field.imaginary[ik1] *=asplit;
		}

    outptr[0]= field;

    for (i=1;i<=field.number ;i++)
		for (j=1;j<=field.number ;j++){
			ik1=(i-1)*field.number+j-1;
			field.real[ik1] *= csplit;
			field.imaginary[ik1] *= csplit;
		}

    outptr[1]= field;

	return outptr;
}

double* lp_b_split( FIELD field, double split )
{
    int i,j;  
    double asplit, bsplit, csplit;
	long ik1, array_space;

	double *outptr;

	array_space = 2*(field.number)*(field.number)+30;

	outptr=(double*) calloc( 2*array_space, sizeof(double) );

    asplit=sqrt(split);
    bsplit=sqrt(1-split);
    csplit=bsplit/asplit;

	for (i=1;i<=field.number ;i++) 
		for (j=1;j<=field.number ;j++){
			ik1=(i-1)*field.number+j-1;
			field.real[ik1] *= asplit;
			field.imaginary[ik1] *=asplit;
		}

    field_array(field,outptr);

    for (i=1;i<=field.number ;i++)
		for (j=1;j<=field.number ;j++){
			ik1=(i-1)*field.number+j-1;
			field.real[ik1] *= csplit;
			field.imaginary[ik1] *= csplit;
		}

    field_array(field,outptr+array_space);

	return outptr;

}
