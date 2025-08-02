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

/*
   allows for operations on the phase/intensity
 
   parameter: C1 C2 F N
        C1 is character constant with valid values int and pha
        C2 is character constant with valid values mult and subst
        F is the name of a file containing intensity/phase mask
          It may be portable graymap or anymap (*.pgm, *.pnm) 
		  with grayscale data
          The number of values in F must be the same as the grid dimension
        N is argument indicating that the intensity mask should
          be normalized before applying

    Examples:
      int mult aa: filter THROUGH  intensity mask from file aa
      pha subst aa: substitute the phase with phase taken from file aa
      pha mult aa: filter the field through phase filter aa
      int subst aa: substitute intensity with one taken from file aa
      int subst aa haha: substitute intensity with a normalized
         one taken from file aa
*/

/*--------------------------------------------------------------*/
/*      (C) Gleb Vdovin 1993-1999                               */
/*      This file is a part of LightPipes package               */
/*      Send bug reports to gleb@okotech.com                    */
/*                                                              */
/*--------------------------------------------------------------*/

#include "lightPipes.h"

FIELD lp_fil_ter(FIELD field, char *c1,char *c2,char *fname,int norm )
{
	FILE *fr; long ik;
	double *sum, cc, smax, ssum, fi, cab,sab,phi,ccc;
	int i,j, bin_ind, ind;
	char *first_string, c_int;
	unsigned char b_in;

	ind=0;

	first_string = (char *) calloc (110,1);

	if((fr=fopen(fname,"r"))==NULL){
		printf("error opening file %s \n", fname);
		exit(1);
	} 
	
	bin_ind =0;

	i=0;
	while((c_int=getc(fr)) != '\n') {
		first_string[i]=c_int;
		i++;
		if(i>100) i=100;
	}
	
	if((strstr(first_string, "P2")) != NULL ) {
	//	printf("Portable ASCII graymap detected \n"); 
	}
	else if((strstr(first_string, "P5")) != NULL ){
	//	printf(stderr,"Portable binary graymap detected \n");
		bin_ind=1;
	}    
	else goto l1;
 
	do{
		i=0;
		while((c_int=getc(fr)) != '\n') {
			first_string[i]=c_int;
			i++;
			if(i>100) i=100;
		}
		/*  printf(stderr, "%s\n", first_string); */
	} while (first_string[0] == '#');
 
	do{ 
		i=0;
		while((c_int=getc(fr)) != '\n') {
			first_string[i]=c_int;
			i++;
			if(i>100) i=100;
		}
		/*  printf(stderr, "%s\n", first_string); */
	} while (first_string[0] == '#');
	goto l2;

l1:
	rewind(fr);

l2:
	if ((strstr(c1, "in"))!= NULL) {
		if ((strstr(c2, "su"))!=NULL){
			/*  int subst here */
			sum  = (double *) calloc( (field.number)*(field.number), sizeof(double) );
			if(sum == NULL){
				printf("fil_ter int subst: Allocation error, exiting\n");
				exit(1);
			}

			smax=0;   
			
			ik=0;
			for (i=1;i<= field.number;i += 1){
				for (j=1;j <= field.number;j += 1){
					if(bin_ind) {
						if(fread (&b_in, sizeof(unsigned char), 1, fr) != 1){
							fprintf(stderr,"Error reading portable bitmap\n");
							exit(1);
						}
						sum[ik]= (double) b_in;
					}
					else{
						if ((fscanf(fr,"%le",&fi))==EOF){
							printf("fil_ter int subst: end of input file reached, exiting\n");
							exit(1);
						}
						sum[ik]=fi;
					}

					if(sum[ik] < 0 && ind == 0){ 
						printf("fil_ter int subs  warning: the\
									   intensity is negative\n"); 
						ind=1;
					}
					
					if(smax < sum[ik]) smax=sum[ik];
					ik++;
				}
			}
			
			if (norm) smax=1.;

			ik=0;
			for (i=1;i<= field.number;i += 1){
				for (j=1;j <= field.number;j += 1){
					ssum=sqrt(sum[ik]/smax);
					phi=phase(field.imaginary[ik], field.real[ik]); 
					field.real[ik]=ssum*cos(phi);
					field.imaginary[ik] = ssum*sin(phi);
					ik++;
				}
			}
			
			free(sum);
		}
		
		else { 
			/* int mult here */
			sum  = (double *) calloc( (field.number)*(field.number), sizeof(double) );
			if(sum == NULL){
				printf("fil_ter int subst: Allocation error, exiting\n");
				exit(1);
			}
			smax=0;
			
			ik=0;
			for (i=1;i<= field.number;i += 1){
				for (j=1;j <= field.number;j += 1){
					if(bin_ind) {
						if(fread (&b_in, sizeof(unsigned char), 1, fr) != 1){
							printf("Error reading portable bitmap\n");
							exit(1);
						}
						
						sum[ik]= (double) b_in;
					}
					else{
						if ((fscanf(fr,"%le",&fi))==EOF){
							printf("fil_ter int subst: end of input file reached, exiting\n");
							exit(1);
						}
						sum[ik]=fi;
					} 
             
					if(sum[ik] < 0 && ind == 0){
						printf("fil_ter int mult  warning: the\
									   intensity is negative\n"); 
						ind = 1;
					}
					
					if(smax < sum[ik]) smax=sum[ik];
					ik++;
				}
			}

			if (norm) smax=1.;
			
			ik=0;
			for (i=1;i<= field.number;i += 1){
				for (j=1;j <= field.number;j += 1){
					ssum=sqrt(sum[ik]/smax);
                    field.real[ik] *= ssum;
					field.imaginary[ik] *= ssum;
					ik++;
				}
			}
			free(sum);
		}
    }

    else { 
		if (strstr(c2, "su")!=NULL){
			/*  pha  subst here */
			ik=0;
			for (i=1;i<= field.number;i += 1){
				for (j=1;j <= field.number;j += 1){
					if(bin_ind) {
						if(fread (&b_in, sizeof(unsigned char), 1, fr) != 1){
							printf("Error reading portable bitmap\n");
							exit(1);
						}
						fi = (double) b_in;
					}
					else{
						if ((fscanf(fr,"%le",&fi))==EOF){
							printf("fil_ter pha subst: end of input file reached, exiting\n");
							exit(1);
						}
					} 

					cab=cos(fi);
					sab=sin(fi);
					ccc=sqrt(field.real[ik]*field.real[ik]+field.imaginary[ik]*field.imaginary[ik]);
					field.real[ik]= ccc*cab;
					field.imaginary[ik] = ccc*sab;
					ik ++;
				}
			}
		}
		
		else { 
			/* pha  mult here */
			ik=0;
			for (i=1;i<= field.number;i += 1){
				for (j=1;j <= field.number;j += 1){
					if(bin_ind) {
						if(fread (&b_in, sizeof(unsigned char), 1, fr) != 1){
							printf("Error reading portable bitmap\n");
							exit(1);
						}
						fi = (double) b_in;
					}
					else{
						if ((fscanf(fr,"%le",&fi))==EOF){
							printf("fil_ter int subst: end of input file reached, exiting\n");
							exit(1);
						}
					} 

					cab=cos(fi);
					sab=sin(fi);
					cc=field.real[ik]*cab-field.imaginary[ik]*sab;
					field.imaginary[ik]=field.real[ik]*sab+field.imaginary[ik]*cab;
					field.real[ik]=cc;
					ik++;
				}
			}
		}
	}

	fclose(fr);

	free(first_string);

	return ( field );

}
