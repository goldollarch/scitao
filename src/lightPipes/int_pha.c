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
    Based on the software LightPipes, there are some changes to accord 
with the Scilab/Scicos environment. the author thanks Dr. Gleb Vdovin 
for developing the excellent software package LightPipes.
*/

#include "lightPipes.h"

double *lp_field_int( double *inptr, int amp)
{
	int i;
	long  size;

	FIELD field;     
    field=array2field( inptr ); 

	size=field.number*field.number2;

	amp_phase_field( &field );

	if(!amp) for(i=0;i<size;i++)
		field.real[i] *= field.real[i];

	return(field.real);

}

double *lp_field_pha( double *inptr)
{
	FIELD field;     
    field=array2field( inptr ); 
	amp_phase_field( &field );
	return(field.imaginary);
}
