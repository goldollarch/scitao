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

FIELD lp_file_int( FIELD field, char *fnam, int imax, double time, int gnuplot )
{
    FILE *fr;char *fname;
    int i,j, istep;
	double dx,dx2;
    long ik1;

	dx=field.size/(field.number-1.);
	dx2=dx*dx;

    istep=1;

    if ( imax==0 ) // equals to grid sampling if you pass the zero to imax
		imax=field.number;

	if(imax>field.number)
		imax=field.number;

    if(field.number/imax >1) {
		istep= field.number/imax;
		imax=field.number/istep;
	}

	if( time<0 ) 
		fname = simple_filename(fnam,".txt");
	else 
		fname = current_filename(fnam,time,".txt");

	if((fr=fopen(fname,"w"))==NULL){
		printf("error opening file %s \n", fname);
		exit(1);
	} 
 
    /* writing the intensity     */
    for (i=1;i<= field.number;i += istep) {
		for (j=1;j <= field.number;j += istep) { 
			double sum;
			sum=0;
		    ik1=(i-1)*field.number +j- 1; 
		    sum += field.real[ik1] *field.real[ik1]+
				field.imaginary[ik1] *field.imaginary[ik1];
			if(gnuplot) fprintf(fr,"%e\n", sum );
			else fprintf(fr,"%e\t", sum );;
		}
		
		fprintf(fr,"\n");
	}

    fclose(fr);

	return ( field );

}
