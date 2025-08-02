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

FIELD lp_b_mix(FIELD field, FIELD input)
{
    long ik1;
    int i,j, number;

	if(field.number != input.number){
	    printf("b_mix: You can not mix grids with \
					   different dimensions!\n");
	    exit(1);
	}

	number = field.number;

	if(field.size != input.size){
	    printf("b_mix: You can not mix grids with \
					   different sizes!\n");
	    exit(1);
	}

	if(field.lambda != input.lambda){
	    printf("b_mix: You can not mix grids with \
					   different wavelengths!\n");
	    exit(1);
	}

    if(input.int2 >field.int2) field.int2 = input.int2;

	if(field.double1 != input.double1){
	    printf("b_mix: You can not mix grids originated \
					   from different coordinate systems!\n");
	    exit(1);
	}

	for (i=1;i<=number ;i++) {
	    for (j=1;j<=number ;j++){
			ik1=(i-1)*number+j-1;
			field.real[ik1] += input.real[ik1];
			field.imaginary[ik1] += input.imaginary[ik1];
	    }
	}

	return ( field );

}

