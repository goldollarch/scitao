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

extern char* simple_filename();
extern char* current_filename();

FIELD lp_cros_out( FIELD field,char *fnam, double time )
{
    FILE *fr;char *fname;
    double *int1, int2, *phase1,phase2, dx;
    int i,j, jj;  long ik1;

    dx=field.size/(field.number-1.);

    /* allocating memory for output arrays  */

    int1=(double *) calloc(field.number, sizeof(double));
	if(int1 == NULL) { 
		printf( "Unable to alloc memory for the int1\n" );
		exit(1);
	}

    phase1=(double *) calloc(field.number, sizeof(double));
	if(phase1 == NULL) { 
		printf( "Unable to alloc memory for the phase1\n" );
		exit(1);
	}

    /* writing the intensity into arrays  */

    i=field.number/2+1;
    for (j=1;j<=field.number;j += 1){ 
		ik1=(i-1)*field.number+j-1;
		jj=j-1; 
		int1[jj]=field.real[ik1]*field.real[ik1]+field.imaginary[ik1]*field.imaginary[ik1]; 
		phase1[jj]=phase(field.imaginary[ik1],field.real[ik1]);
	}

	if( time<0 ) 
		fname = simple_filename(fnam,".txt");
	else 
		fname = current_filename(fnam,time,".txt");

	if((fr=fopen(fname,"w"))==NULL) {
		printf( "error opening file %s\n",fname );
		exit(1);
	} 
 
    j=field.number/2+1;
    for (i=1;i<=field.number;i += 1){
		double cc;
		ik1=(i-1)*field.number+j-1;
		jj=i-1;
		cc=dx*(i-field.number/2-1);

		int2=field.real[ik1]*field.real[ik1]+field.imaginary[ik1] *field.imaginary[ik1]; 
		phase2=phase(field.imaginary[ik1], field.real[ik1]);
		fprintf(fr," %e %e %e %e %e\n", cc, int1[jj], int2, phase1[jj], phase2);

    }

	free(int1);	free(phase1);
    fclose(fr);

	return ( field );

}
