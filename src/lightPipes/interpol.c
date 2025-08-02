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
   transfers the field into a grid with different size and dimension

    parameters: B [N, X_shift, Y_shift, A, M]  
       where B is the new size of the grid and  N is new dimension
         you can transfer 0 to B and N to preserve the existing values 
         (if you only want to shift, rotate or scale)
       X_shift and Y_shift are the shifts of the field in a new grid
       A is the angle of rotation (after shifts) and
       M (M><0) is the magnification (the last applied)
*/

/*--------------------------------------------------------------*/
/*      (C) Gleb Vdovin 1993-1999                               */
/*      This file is a part of LightPipes package               */
/*      Send bug reports to gleb@okotech.com                    */
/*                                                              */
/*--------------------------------------------------------------*/

#include "lightPipes.h"

double *lp_interpol( FIELD field, double size_new, int new_number, 
		double x_shift,	double  y_shift, double  angle, double magnif )
{
	double dx_new, dx_old, x_new, x_old, y_new, y_old, 
		old_number, newfloat, inv_squares(), lower,
		upper, on21, nn21, ss, cc, x0, y0, change, int16();
	long  n_old_max, kk;

	double z00, z01, z02 ,z03,z10, z11 ,z12, z13, 
		z20 ,z21, z22, z23,z30, z31, z32, z33;
	int i, j, i_old, j_old;

	double *outptr;

	if (size_new==0) size_new=field.size;
	if (new_number==0) new_number=field.number;

	outptr=(double*) calloc( 2*(new_number)*(new_number)+30, sizeof(double) );

    old_number=field.number;
    n_old_max=old_number*(old_number-1)-1;

	dx_new=size_new/(new_number-1.); 
    dx_old=field.size/(old_number-1.);

	change= 1.;
    angle *= Pi2/360.;

    field.number=new_number;
    field.number2=new_number;
    field.size=size_new;

	field_header_array(field,outptr);

    on21= (int) old_number/2+1;
    nn21= (int) new_number/2+1;
    lower= (1-on21)*dx_old;
    upper= (old_number-on21)*dx_old;

    cc=cos(angle);
    ss=sin(angle);

	kk=30;
    for (i=1; i<=new_number; i++)
	{
        for (j=1; j<=new_number; j++)
		{
			long n_tmp;
			int i_small, i_local;

            x0=(i-nn21)*dx_new-x_shift;
            y0=(j-nn21)*dx_new-y_shift;
            x_new=(x0*cc+y0*ss)/magnif;
            y_new=(-x0*ss+y0*cc)/magnif; 

            i_old=(int) floor(x_new/dx_old)+on21;
            x_old=(i_old-on21)*dx_old; 
            j_old= (int)floor( y_new/dx_old)+on21;
            y_old=(j_old-on21)*dx_old;
			
			i_small=0;
			i_local = 0;

			/*          first row  */	
			n_tmp=(i_old-2)*old_number+j_old-2;

			if( n_tmp < 0 || n_tmp > n_old_max)	{ i_small=1; i_local=1;	} 
			
			if (i_local != 1) z00=field.real[n_tmp];
			else z00=0;
			
			i_local=0;

            n_tmp += old_number;
			
			if( n_tmp < 0 || n_tmp > n_old_max) { i_small=1; i_local=1;	}
			if (i_local != 1) 
				z10=field.real[n_tmp];
			else z10=0;i_local=0;

            n_tmp += old_number;
			if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
			if (i_local != 1) z20=field.real[n_tmp];
			else z20=0;i_local=0;

			n_tmp += old_number;
			if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
			if (i_local != 1) z30=field.real[n_tmp];
			else z30=0;i_local=0;

			/*          second row */
			n_tmp=(i_old-2)*old_number+j_old-1;
			if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
			if (i_local != 1) z01=field.real[n_tmp];
			else z01=0;i_local=0;

            n_tmp += old_number;
			if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
			if (i_local != 1) z11=field.real[n_tmp];
			else z11=0;i_local=0;

            n_tmp += old_number;
			if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
			if (i_local != 1) z21=field.real[n_tmp];
			else z21=0;i_local=0;

			n_tmp += old_number;
			if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
			if (i_local != 1) z31=field.real[n_tmp];
			else z31=0;i_local=0;

			/*          third row */
			n_tmp=(i_old-2)*old_number+j_old;
			if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
			if (i_local != 1) z02=field.real[n_tmp];
			else z02=0;i_local=0;

            n_tmp += old_number;
			if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
			if (i_local != 1) z12=field.real[n_tmp];
			else z12=0;i_local=0;

            n_tmp += old_number;
			if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
			if (i_local != 1) z22=field.real[n_tmp];
			else z22=0;i_local=0;

			n_tmp += old_number;
			if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
			if (i_local != 1) z32=field.real[n_tmp];
			else z32=0;i_local=0;

			/*          fourth     row */
			n_tmp=(i_old-2)*old_number+j_old+1;
			if( n_tmp < 0 || n_tmp > n_old_max){i_small=1; i_local=1;}
			if (i_local != 1) z03=field.real[n_tmp];
			else z03=0;i_local=0;

            n_tmp += old_number;
			if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
			if (i_local != 1) z13=field.real[n_tmp];
			else z13=0;i_local=0;

            n_tmp += old_number;
			if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
			if (i_local != 1) z23=field.real[n_tmp];
			else z23=0;i_local=0;

			n_tmp += old_number;
			if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
			if (i_local != 1) z33=field.real[n_tmp];
			else z33=0;i_local=0;

			if (i_small == 1)
			{
				if(x_new > lower && x_new < upper && y_new >lower && y_new < upper)
				{
					newfloat = inv_squares(x_old,y_old,dx_old,z11,z21,z12,z22,x_new,y_new)/magnif;
				}
				else newfloat=0.;
			}
			else{ 
				if(x_new > lower && x_new < upper && y_new >lower && y_new < upper)
				{
					newfloat=int16(x_old, y_old, x_new, y_new, dx_old,z00, z01, z02 ,z03,
						z10, z11 ,z12, z13,z20 ,z21, z22, z23,z30, z31, z32, z33)/magnif;
				}
				else newfloat=0.;                
			}
			
			outptr[kk] =newfloat;
			kk ++;
        }
	}

	/* phase is here */
	for (i=1; i<=new_number; i++)
	{
        for (j=1; j<=new_number; j++)
		{
			long n_tmp;
			int i_small, i_local;

            x0=(i-nn21)*dx_new-x_shift;
            y0=(j-nn21)*dx_new-y_shift;
            x_new=(x0*cc+y0*ss)/magnif;
            y_new=(-x0*ss+y0*cc)/magnif; 

            i_old=(int) floor(x_new/dx_old)+on21;
            x_old=(i_old-on21)*dx_old; 
            j_old= (int)floor( y_new/dx_old)+on21;
            y_old=(j_old-on21)*dx_old;

			i_small=0;
			i_local = 0;

			/*          first row  */	
			n_tmp=(i_old-2)*old_number+j_old-2;

			if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;} 
			if (i_local != 1) z00=field.imaginary[n_tmp];
			else z00=0;
			i_local=0;

            n_tmp += old_number;
			if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
			if (i_local != 1) z10=field.imaginary[n_tmp];
			else z10=0;i_local=0;

            n_tmp += old_number;
			if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
			if (i_local != 1) z20=field.imaginary[n_tmp];
			else z20=0;i_local=0;

			n_tmp += old_number;
			if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
			if (i_local != 1) z30=field.imaginary[n_tmp];
			else z30=0;i_local=0;
			
			/*          second row */
			n_tmp=(i_old-2)*old_number+j_old-1;
			if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
			if (i_local != 1) z01=field.imaginary[n_tmp];
			else z01=0;i_local=0;

            n_tmp += old_number;
			if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
			if (i_local != 1) z11=field.imaginary[n_tmp];
			else z11=0;i_local=0;

            n_tmp += old_number;
			if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
			if (i_local != 1) z21=field.imaginary[n_tmp];
			else z21=0;i_local=0;

			n_tmp += old_number;
			if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
			if (i_local != 1) z31=field.imaginary[n_tmp];
			else z31=0;i_local=0;

			/*          third row */
			n_tmp=(i_old-2)*old_number+j_old;
			if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
			if (i_local != 1) z02=field.imaginary[n_tmp];
			else z02=0;i_local=0;

            n_tmp += old_number;
			if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
			if (i_local != 1) z12=field.imaginary[n_tmp];
			else z12=0;i_local=0;

            n_tmp += old_number;
			if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
			if (i_local != 1) z22=field.imaginary[n_tmp];
			else z22=0;i_local=0;

			n_tmp += old_number;
			if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
			if (i_local != 1) z32=field.imaginary[n_tmp];
			else z32=0;i_local=0;

			/*          fourth     row */
			n_tmp=(i_old-2)*old_number+j_old+1;
			if( n_tmp < 0 || n_tmp > n_old_max){i_small=1; i_local=1;}
			if (i_local != 1) z03=field.imaginary[n_tmp];
			else z03=0;i_local=0;

            n_tmp += old_number;
			if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
			if (i_local != 1) z13=field.imaginary[n_tmp];
			else z13=0;i_local=0;

            n_tmp += old_number;
			if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
			if (i_local != 1) z23=field.imaginary[n_tmp];
			else z23=0;i_local=0;

			n_tmp += old_number;
			if( n_tmp < 0 || n_tmp > n_old_max) {i_small=1; i_local=1;}
			if (i_local != 1) z33=field.imaginary[n_tmp];
			else z33=0;i_local=0;

			if (i_small == 1){
				if(x_new > lower && x_new < upper && y_new >lower && y_new < upper)
				{
					newfloat=inv_squares(x_old,y_old,dx_old,z11,z21,z12,z22,x_new,y_new)/magnif;
				}
				else newfloat=0.;
			}
			else
			{ 
				if(x_new > lower && x_new < upper && y_new >lower && y_new < upper)
				{
					newfloat=int16(x_old, y_old, x_new, y_new, dx_old,
						z00, z01, z02 ,z03,z10, z11 ,z12, z13,z20 ,z21, z22, z23,
						z30, z31, z32, z33)/magnif;
				}
				else newfloat=0.;	
			}
			
			outptr[kk] =newfloat;
			kk ++;
	    }
	}

	return outptr;

}


////////////////////////////////////////////////////////////////////////

double int16( xc, yc, xp, yp, dx, z00, z01, z02 ,z03,
	 z10, z11 ,z12, z13, z20 ,z21, z22, z23, z30, z31, z32, z33)
double  xc, yc, xp, yp, dx,z00, z01, z02 ,z03, 
z10, z11 ,z12, z13, z20 ,z21, z22, z23, z30, z31, z32, z33;
{ 
	double z1, z2, z3, z4, zz1, zz2, zz3, zz4, int4();

	z1=z00; z2=z01; z3=z02, z4=z03;
	zz1=int4(yc, dx, z1, z2, z3, z4, yp);

	z1=z10; z2=z11; z3=z12, z4=z13;
	zz2=int4(yc, dx, z1, z2, z3, z4, yp);
	
	z1=z20; z2=z21; z3=z22, z4=z23;
	zz3=int4(yc, dx, z1, z2, z3, z4, yp);

	z1=z30; z2=z31; z3=z32, z4=z33;
	zz4=int4(yc, dx, z1, z2, z3, z4, yp);

	return int4(xc,dx,zz1,zz2,zz2,zz3,xp);

}

/* cubic four-point interpolation 

a+b*x1+c*x1^2+d*x1^3=y1
......................
a+b*x4+c*x4^2+d*x4^3=y4

where 
x1=x2-dx; 
x3=x2+dx; 
x4=x2+2*dx; 
the grid is uniform

*/

double int4(x2,dx,y1,y2,y3,y4,xz)
double dx,x2,y1,y2,y3,y4,xz;
{
	double  t0,a,b,c,d;
	t0=1./(dx*dx*dx); 

	a = t0*(3.0*y2*x2*dx*dx+x2*x2*x2*y1+3.0*x2*x2*y1*dx+2.0*x2*y1*\
		dx*dx-x2*x2*x2*y4+3.0*x2*x2*x2*y3+3.0*x2*x2*y3*dx-6.0*x2*y3*dx*dx+x2*y4*dx*dx+\
		6.0*y2*dx*dx*dx-3.0*y2*x2*x2*x2-6.0*y2*x2*x2*dx)/6;

	b = -t0*(-6.0*y3*dx*dx+y4*dx*dx+3.0*y2*dx*dx+9.0*x2*x2*y3-3.0*\
		x2*x2*y4-9.0*y2*x2*x2-12.0*y2*x2*dx+6.0*y3*x2*dx+3.0*x2*x2*\
		y1+6.0*y1*x2*dx+2.0*y1*dx*dx)/6;

	d = -t0*(-3.0*y2-y4+y1+3.0*y3)/6;
	c = t0*(-3.0*y2*x2+y3*dx+y1*dx-2.0*y2*dx-x2*y4+3.0*x2*y3+x2*y1)/2;
 
	return a+xz*(b+xz*(c+xz*d)); 

}
