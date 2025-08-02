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

double d3();int unf3_i_min(),ubbf2();
double d5();int unf4_i_min(),ubbf3();

double* unf3(double koeff, double rad, int nn, double *x)
{
	int  i, j, ii, ii1, ii2,  ilong=2, icount;
	long ik ;

	for(i=1; i<=nn-1; i++){
		for (j=1; j<=nn-1; j++){
			int ii, ij;
			ii = i-nn/2+1;
			ij=j-nn/2+1; 
			ii1=(i-1)*nn+j-1;
			if (ii*ii+ij*ij > rad*rad*nn*nn/4) x[ii1]=0.;
		}
	}

	for (icount=1; icount <= ilong; icount ++){
		for(i=nn/2; i<=nn-1; i++){
			for (j=nn/2; j<=nn-1; j++){
				ii1=(i-1)*nn+j-1;
				ii2=(i-1)*nn +j-1+1;
				ii=ubbf2(&x[ii1],&x[ii2]);
			}
		} 

		for(i=nn/2; i>=2; i--){
			for (j=nn/2; j<=nn-1; j++){
				ii1=(i-1)*nn+j-1;
				ii2=(i-1-1)*nn +j-1;
				ii=ubbf2(&x[ii1],&x[ii2]);
			}
		} 
		
		for(i=nn/2; i<=nn-1; i++){
			for (j=nn/2; j>=2; j--){
				ii1=(i-1)*nn+j-1;
				ii2=(i-1+1)*nn +j-1;
				ii=ubbf2(&x[ii1],&x[ii2]);
			}
		} 
		
		for(i=nn/2; i>=2; i--){
			for (j=nn/2; j>=2; j--){
				ii1=(i-1)*nn+j-1;
				ii2=(i-1)*nn +j-1-1;
                
				ii=ubbf2(&x[ii1],&x[ii2]);
			}
		} 

	}

	for(i=1; i<=nn; i++){
		for (j=1; j<=nn; j++){
			ik=(i-1)*nn+j-1;
			x[ik] *= koeff;
		}
	}

	return x;

}

int ubbf2( a1, a2)
double *a1, *a2; { 
	double s[5], aa2; 
	int ic,  ik;
	
l1:
	ik=1;

	for(aa2= *a2 - Pi2; aa2<= *a2 + Pi2+0.1; aa2 += Pi2){
		s[ik]=d3( *a1, aa2 );
		ik++;
	}
 
	ic=unf3_i_min(s[1],s[2],s[3]);
	if (ic==1){
		*a2 -= Pi2;
		goto l1;
	}

	if (ic==2){
		goto l2;
	}
	
	if (ic==3){
		*a2 += Pi2;
		goto l1;
	}

l2: 
	return ic;
}

double d3(a1,a2)
double a1,a2;
{
	return fabs(a1-a2);
}

int unf3_i_min(a1,a2,a3)
double a1,a2,a3;
{ 
	int ii;
	double dmin;
	if (a2 < a1) {
		dmin=a2; 
		ii=2;
	}
	
	else {
		dmin=a1;  
		ii=1;
	}
	
	if (a3<dmin) {
		dmin=a3;
		ii=3;
	}
	
	return ii;
}

////////////////////////////////////////////////////////////////////

double* unf4(double koeff, int nn, double *x)
{
	int  i, j, ii, ii1, ii2, ii3;
	long ik;
	
	for(i=1; i<=nn-1; i++){
		for (j=1; j<=nn-1; j++){
			ii1=(i-1)*nn+j-1;
			ii2=(i-1+1)*nn +j-1;
			ii3=(i-1)*nn +j-1+1;
			ii=ubbf3(&x[ii1],&x[ii2],&x[ii3]);
		}
	} 
	
	for(i=nn; i>=2; i--){
		for (j=1; j<=nn-1; j++){
			ii1=(i-1)*nn+j-1;
			ii2=(i-1-1)*nn +j-1;
			ii3=(i-1)*nn +j-1+1;
			ii=ubbf3(&x[ii1],&x[ii2],&x[ii3]);
		}
	} 
	
	for(i=1; i<=nn-1; i++){
		for (j=nn; j>=2; j--){
			ii1=(i-1)*nn+j-1;
			ii2=(i-1+1)*nn +j-1;
			ii3=(i-1)*nn +j-1-1;
			ii=ubbf3(&x[ii1],&x[ii2],&x[ii3]);
		}
	} 
	
	for(i=1; i<=nn; i++){
		for (j=1; j<=nn; j++){
			ik=(i-1)*nn+j-1;
			x[ik] *= koeff;
		}
	}
	
	return x;
    
}

int ubbf3( a1, a2, a3)
double *a1, *a2, *a3; { 
	double s[10], aa2, aa3;
	int ic,  ik;
 
l1:
    ik=1;
	for(aa2= *a2 - Pi2; aa2<= *a2 + Pi2+0.1; aa2 += Pi2)
	{
		for(aa3= *a3 - Pi2; aa3<= *a3 + Pi2+0.1; aa3 += Pi2)
		{
			s[ik]=d5( *a1, aa2, aa3 );
			ik++;
		}
	}
 
	ic=unf4_i_min(s[1],s[2],s[3],s[4],s[5],s[6],s[7],s[8],s[9]);

	if (ic==1){	*a2 -= Pi2;	*a3 -= Pi2;	goto l1;}
	if (ic==2){	*a2 -= Pi2;	goto l1;}
	if (ic==3){	*a2 -= Pi2;	*a3 += Pi2;	goto l1;}
	if (ic==4){	*a3 -= Pi2;	goto l1;}
	if (ic==5){	goto l2;}
	if (ic==6){	*a3 += Pi2;	goto l1;}
	if (ic==7){	*a2 += Pi2;	*a3 -= Pi2;	goto l1;}
	if (ic==8){	*a2 += Pi2;	goto l1;}
	if (ic==9){	*a2 += Pi2;	*a3 += Pi2;	goto l1;}

l2: 
	return ic;
}

double d5(a1,a2,a3)
double a1,a2,a3;
{
	return fabs(a1-a2)+fabs(a1-a3);
}

int unf4_i_min(a1,a2,a3, a4, a5, a6, a7, a8, a9)
double a1,a2,a3,a4,a5,a6,a7,a8,a9;
{ 
	int ii;
	double dmin;

	if (a2 < a1) {	dmin=a2; ii=2;}
	else { dmin=a1; ii=1; }

	if (a3<dmin) {	dmin=a3; ii=3;	}
	if (a4<dmin) {	dmin=a4; ii=4;	}
	if (a5<dmin) {	dmin=a5; ii=5;	}
	if (a6<dmin) {	dmin=a6; ii=6;	}
	if (a7<dmin) {	dmin=a7; ii=7;	}
	if (a8<dmin) {	dmin=a8; ii=8;	}
	if (a9<dmin) {	dmin=a9; ii=9;	}  
	
	return ii;
}
