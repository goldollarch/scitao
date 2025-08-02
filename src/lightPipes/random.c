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
 
#include <time.h>
#include "lightPipes.h"

FIELD lp_random( FIELD field, char *job, double ampl, unsigned int seed )
{
    int i,j;
    
    long ik;

	if(seed==0)   
		seed=(unsigned int) time((unsigned int) 0);
	srand(seed);
    
    if ((strstr(job, "ph"))!= NULL){
		ik=0;
		for (i=1;i<= field.number;i += 1){
			for (j=1;j <= field.number;j += 1){ 
				double fi, cab,sab, my_rnd, cc;
				my_rnd= ((double) rand()) / ((double) RAND_MAX)-0.5;
				fi=my_rnd*ampl;
				cab=cos(fi);
				sab=sin(fi);
				cc=field.real[ik]*cab-field.imaginary[ik]*sab;
				field.imaginary[ik]=field.real[ik]*sab+field.imaginary[ik]*cab;
				field.real[ik]=cc;
				ik++;
			}
		}
	}

	else if ((strstr(job, "in"))!= NULL){
		double maxint =0.;
		ik=0;
		for (i=1;i<= field.number;i += 1){
			for (j=1;j <= field.number;j += 1){ 
				double   cc;
				cc=field.imaginary[ik]*field.imaginary[ik] + 
					field.real[ik]*field.real[ik];
				if (cc > maxint) maxint=cc; 
				ik++;
			}
		}
		
		ik=0;
		for (i=1;i<= field.number;i += 1){
			for (j=1;j <= field.number;j += 1){ 
				double  cab, sab, my_rnd, cc, phi;
				my_rnd= ((double) rand()) / ((double) RAND_MAX);
				phi=phase(field.imaginary[ik],field.real[ik]);
				cc=(field.imaginary[ik]*field.imaginary[ik] + 
					field.real[ik]*field.real[ik])+ maxint*ampl*my_rnd;
				cc=sqrt(cc);
				cab=cos(phi);
				sab=sin(phi);
				field.imaginary[ik]=cc*cab;
				field.real[ik]=cc*sab;
				ik++;
			}
		}
	}

	return ( field );

}
