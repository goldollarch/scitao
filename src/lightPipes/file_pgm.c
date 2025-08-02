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

FIELD lp_file_pgm
(FIELD field, char *fnam, double gamma, int imax, int max_val, double time )
{
    int i,j,ii,jj,istep,i0, i_i;
	long ik1; float max_int=0;
    FILE *fr; char *fname;
	
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
		fname=simple_filename( fnam,".pgm" );
	else 
		fname=current_filename( fnam,time,".pgm" );

	if ( (fr=fopen(fname,"w"))==NULL )
	{
		printf("Can't create file %s.\n", fname);
		exit(0);
	}
  
    if(istep != 1){
		fprintf(fr,"P2\n");
		fprintf(fr,"#Creator: scitao (C) 2006, Chen jingyuan\n");
		fprintf(fr,"%d %d\n", imax-1, imax-1);
		fprintf(fr,"%d\n", max_val);
	}
    else{ 
		fprintf(fr,"P2\n");
		fprintf(fr,"#Creator: scitao (C) 2006, Chen jingyuan\n");
		fprintf(fr,"%d %d\n", imax, imax);
		fprintf(fr,"%d\n", max_val);
	}

    if( istep != 1){
		for (i=1 ; i<= field.number-istep; i +=  istep){
			for (j=1;j <= field.number-istep;j += istep) {
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

		i_i=1;
		for (i=1; i<= field.number-istep; i += istep){
			for (j=1;j <= field.number-istep;j += istep){ 
				double sum;
				sum=0;
				for (ii=i; ii<=i+istep; ii++)
					for(jj=j; jj<=j+istep; jj++){
						ik1=(ii-1)*field.number +jj- 1;
						sum +=field.real[ik1] *field.real[ik1]+ 
							field.imaginary[ik1] *field.imaginary[ik1];
					}
					
				sum=sum/(istep*istep);
				i0=(int) floor(pow((sum/max_int),1./(gamma+0.0001))*max_val);
				fprintf(fr,"%d ", i0);
				i_i++;
				if (i_i == 40){
					fprintf(fr,"\n");
					i_i=1;
				}
			}
		}
	}
	
	else{
        
		for (i=1; i<= field.number; i++){
			for (j=1 ;j <= field.number; j++ ){
				double sum;
				ik1=(i-1)*field.number +j- 1;
				sum =field.real[ik1] *field.real[ik1]+ 
					field.imaginary[ik1] *field.imaginary[ik1];
				if(sum>max_int) max_int=sum;
			}
		}
		
		i_i=1;
		
		for (i=1; i<= field.number; i++ ){
			for (j=1; j <= field.number; j++ ){
				double sum;
				ik1=(i-1)*field.number +j- 1;
				sum =field.real[ik1] *field.real[ik1]+ 
					field.imaginary[ik1] *field.imaginary[ik1];
				i0=(int) floor(pow((sum/max_int),1./(gamma+0.00001))*max_val);
                
				fprintf(fr,"%d ", i0);
				i_i++;
				if (i_i == 40){
					fprintf(fr,"\n");
					i_i=1;
				}
			}
		}
	}

    fclose(fr);

	return ( field );

}
