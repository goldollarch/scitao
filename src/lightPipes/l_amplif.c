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

#include "lightPipes.h"

FIELD lp_l_amplif( FIELD field, double gain, double base, double i_sat )
{
    int i,j;
    double ss, intensity, dx,dx2;
    long ik1;

    dx =field.size/(field.number);
    dx2 = dx*dx;
       
    ik1=0;
    for (i=1;i<=field.number ;i++){
        for (j=1;j<=field.number ;j++){
			intensity=(field.real[ik1]*field.real[ik1]+
				field.imaginary[ik1]*field.imaginary[ik1]);
			ss=exp(base*(gain/(1.+(intensity/i_sat))));

			field.real[ik1] *= ss;
			field.imaginary[ik1] *= ss;
			ik1++;
		}
	}

	return ( field );

}
