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

/*         file_ps writes intensity into *.ps file F
         parameter: F N G. 
             where F is the output filename,
                   N is the grid size, N=128 if defaulted  
                   G is the $\gamma$ parameter, [0.1...10],
				      higher G gives better contrast in low intensities, 
					  default G=2.0
         Output file F can be processed with any postscript device
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

FIELD lp_file_ps( FIELD field, char *fnam, double gamma, int imax, double time)
{
	FILE *fr;char *fname;
	int i,j,ii,jj,istep,i0;
	long ik1;float max_int;
	
	max_int=0;
   
    istep=1;

    if ( imax==0 ) // equals to grid sampling if you pass the zero to imax
		imax=field.number;

	if(imax>field.number)
		imax=field.number;

    if(field.number/imax >1) 
	{
		istep= field.number/imax;
		imax=(int)ceil ((float)field.number/(float) istep);
	}

	if( time<0 ) 
		fname = simple_filename( fnam,".ps" );
	else 
		fname = current_filename( fnam,time,".ps" );

	if ( (fr=fopen(fname,"w"))==NULL )
	{
		printf("Can't create file %s.\n", fname);
		exit(0);
	}

    fprintf(fr,"%%!PS-Adobe-2.0 EPSF-2.0\n");
    fprintf(fr,"%%%%BoundingBox: 0 0 256 256  \n");
    fprintf(fr,"%%%%Creator: scitao  2006 Chen jingyuan, beta version 2006 \n");
    fprintf(fr,"%%%%EndComments\n");
    fprintf(fr,"%%%%EndProlog\n\n");

    fprintf(fr,"/origstate save def\n20 dict begin\n");
    fprintf(fr,"/picstr %d string def\n", (imax-1));
    fprintf(fr,"256  256   scale\n");
    fprintf(fr,"%d %d  8\n", (imax-1), (imax-1));
    fprintf(fr,"[ %d 0 0 %d 0 %d]\n",(imax-1), -(imax-1), (imax-1));
    fprintf(fr,"{currentfile\npicstr readhexstring pop}\nimage\n ");
    
    for (i=1; i<= field.number-istep; i +=  istep){
		for (j=1; j <= field.number-istep; j += istep) {
			double sum;
			sum=0;
			for (ii=i; ii<=i+istep; ii++)
				for(jj=j; jj<=j+istep; jj++){
					ik1=(ii-1)*field.number +jj- 1;
					sum += field.real[ik1] *field.real[ik1]+ 
						field.imaginary[ik1] *field.imaginary[ik1];
				}
				
			sum=sum/(istep*istep);
			if(sum>max_int) max_int=sum;
		}
	}
	
	for (i = 1; i<= field.number-istep; i +=  istep) {
		for (j=1;j <= field.number-istep;j += istep) { 
			double sum;
			sum=0;
			for (ii=i; ii<=i+istep; ii++)
				for(jj=j; jj<=j+istep; jj++){
					ik1=(ii-1)*field.number +jj- 1;
					sum +=field.real[ik1] *field.real[ik1]+ 
						field.imaginary[ik1] *field.imaginary[ik1];
				}
			sum=sum/(istep*istep);
			i0=(int) floor(pow((sum/max_int),1./(gamma+0.00001))*255);
			fprintf(fr,"%02x", i0);
		}
	}

    fprintf(fr,"\n\nshowpage\n%%%%Trailer\n");

    fclose(fr);

	return ( field );

}
