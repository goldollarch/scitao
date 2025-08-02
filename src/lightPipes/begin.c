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

FIELD lp_begin(  double size, double lambda, int number, int number2 )
{
    int i,j;
    long ik;

	FIELD field;

    field=create_field( number, number2 );

    field.size=size; field.lambda=lambda;
    field.int1=0; field.int2=0; field.int3=0;
    field.double1=0; field.double2=0; field.double3=0;

	field.real_imag=1; field.interleaved=1;

	/*  Here the initial field is formed   */

    ik=0;
    for (i=1;i<=field.number ;i++)
	{
		for (j=1;j<=field.number2 ;j++)
		{
			field.real[ik]=1.;
			field.imaginary[ik]=0.;
			ik++;
		}
	}

	return ( field );

}
