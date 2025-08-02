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

extern void fresnel();

FIELD lp_lens_fresn( FIELD field, double f, double z )
{
    int i,j;
    long ik;
    double  z1, f1, ampl_scale;
	
	if(f==0)
		printf("lens_forvard: Lens with \
			   ZERO focal length does not exist.\n");

    f1=0.;

    if (field.double1 !=0. ) 
		f1=1./field.double1;
    else 
		f1=10000000.* field.size*field.size/field.lambda;
    if( (f+f1) != 0.)
		f=(f*f1)/(f+f1);
    else
		f=10000000.* field.size*field.size/field.lambda;
	
	z1=-z*f/(z-f);
  
	if(z1 < 0. ) {
		printf("Sorry, lens_fresn can not propagate behind\n\
			   the focal point, use lens_forvard instead, exiting.\n\n");
		exit (1);
    }

    fresnel(z1,field); 

    ampl_scale=(f-z)/f;
    field.size *= ampl_scale;
    field.double1= -1./(z-f);

    ik=0;
    for (i=1;i<=field.number; i++){
        for (j=1;j<=field.number; j++){
			field.real[ik] = field.real[ik]/ampl_scale;
            field.imaginary[ik] =field.imaginary[ik]/ampl_scale; 
            ik++;
        }
    }

	return ( field );

}
