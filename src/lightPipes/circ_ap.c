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

FIELD lp_circ_ap(FIELD field,double R,double x_shift,double y_shift )
{
    int i,j,i2;
    double dx,x,y,rr;
    long ik1;
	
    rr=R*R;

    dx =field.size/(field.number);
    i2=field.number/2;

    /* Cuttitng the aperture      */

    for (i=1;i<=field.number; i++)
	{
        x=(i-i2-1)*dx-x_shift;
        for (j=1;j<=field.number; j++)
		{
            y=(j-i2-1)*dx-y_shift;
            ik1=(i-1)*field.number+j-1; 
            if(x*x+y*y > rr) 
			{
                field.real[ik1]=0.;
                field.imaginary[ik1]=0.; 
            }
        }
    }

	return ( field );

}
