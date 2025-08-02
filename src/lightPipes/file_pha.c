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

FIELD lp_file_pha( FIELD  field, char *fnam, int imax, double time, int gnuplot  )
{
    FILE *fr;char *fname;
    int i, j, istep, ii, jj;
    double  im, re;
    long ik1;
	
    istep=1;

    if ( imax==0 ) // equals to grid sampling if you pass the zero to imax
		imax=field.number;

    if (imax>field.number) 
		imax=field.number;

    if(field.number/imax >1) {
		istep= field.number/imax;
		imax=field.number/istep;
	}

	if( time<0 ) 
		fname = simple_filename(fnam,".txt");
	else 
		fname = current_filename(fnam,time,".txt");

	if ( (fr=fopen(fname,"w"))==NULL )
	{
		printf("Can't create file %s.\n", fname);
		exit(0);
	}
 
    if(istep != 1)
	{
		for (i=1;i < field.number;i+=istep)
		{
			for (j=1;j< field.number;j+=istep)
			{
				re=im=0.;
				for(ii=i; ii<i+istep; ii++)
				{
					for(jj=j; jj<j+istep; jj++)
					{
						ik1=(ii-1)*field.number+jj-1;
						re += field.real[ik1];
						im += field.imaginary[ik1];
					}
				}
				if(gnuplot) fprintf(fr,"%e\n",phase(im,re));
				else fprintf(fr,"%e\t",phase(im,re));;
			}
			fprintf(fr,"\n");
		}
	}
	
	else 
	{ 
		for (i=1;i<= field.number;i+=istep)
		{
			for (j=1;j<= field.number;j+=istep)
			{
				re=im=0.;
				ik1=(i-1)*field.number+j-1;
				re = field.real[ik1];
				im = field.imaginary[ik1];
				if(gnuplot) fprintf(fr,"%e\n",phase(im,re));
				else fprintf(fr,"%e\t",phase(im,re));;
			}
			fprintf(fr,"\n");
		}
	}
	
	fclose(fr);

	return ( field );

}

