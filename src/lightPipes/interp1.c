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

double *lp_interp1(FIELD field, double size_new, int new_number, 
		double x_shift,	double  y_shift,double  angle, double magnif )
{
    double dx_new, dx_old, x_new, x_old, y_new, y_old, 
		old_number, newfloat, inv_squares(), lower,
		upper, on21, nn21, ss, cc, x0, y0, change;
    int i,j,i_old, j_old;
    long n_old, n_oldx, n_oldy, n_oldxy, kk;

	double *outptr;

	if (size_new==0) size_new=field.size;
	if (new_number==0) new_number=field.number;

	outptr=(double*) calloc( 2*(new_number)*(new_number)+30, sizeof(double) );

	old_number=field.number;

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
            x0=(i-nn21)*dx_new-x_shift;
            y0=(j-nn21)*dx_new-y_shift;
            x_new=(x0*cc+y0*ss)/magnif;
            y_new=(-x0*ss+y0*cc)/magnif; 

            i_old=(int) floor(x_new/dx_old)+on21;
            x_old=(i_old-on21)*dx_old; 
            j_old= (int)floor( y_new/dx_old)+on21;
            y_old=(j_old-on21)*dx_old;
            n_old=(i_old-1)*old_number+j_old-1;
            n_oldx=n_old+old_number;
            n_oldy=n_old+1;
            n_oldxy=n_oldx+1;

			if(x_new > lower && x_new < upper && y_new >lower && y_new < upper)
			{
				newfloat=inv_squares(x_old,y_old,dx_old,field.real[n_old],
					field.real[n_oldx],field.real[n_oldy],field.real[n_oldxy],
					x_new,y_new)/magnif;
				newfloat *= change;
            }
            else newfloat=0.;

			outptr[kk] =newfloat;
			kk ++;
		}
	}

    for (i=1; i<=new_number; i++)
	{
        for (j=1; j<=new_number; j++)
		{
            x0=(i-nn21)*dx_new-x_shift;
            y0=(j-nn21)*dx_new-y_shift;
            x_new=(x0*cc+y0*ss)/magnif;
            y_new=(-x0*ss+y0*cc)/magnif; 

            i_old=(int) floor(x_new/dx_old)+on21;
            x_old=(i_old-on21)*dx_old; 
            j_old= (int) floor( y_new/dx_old)+on21;
            y_old=(j_old-on21)*dx_old;
            n_old=(i_old-1)*old_number+j_old-1;
            n_oldx=n_old+old_number;
            n_oldy=n_old+1;
            n_oldxy=n_oldx+1;

            if(x_new > lower && x_new < upper && y_new >lower && y_new < upper)
			{
                newfloat=inv_squares(x_old,y_old,dx_old,field.imaginary[n_old],
					field.imaginary[n_oldx],field.imaginary[n_oldy],
					field.imaginary[n_oldxy],x_new,y_new)/magnif;
				newfloat *= change;
            }
            else newfloat=0.;

			outptr[kk] =newfloat;
			kk ++;
        }
    }

	return outptr;

}

/*
*         Inverse square interpolation : 
*         given square (x,y) (x+dx,y) (x,y+dx) (x+dx,y+dx) 
*         with values   z     zx       zy        zxy 
*         the program returns value for Z for arbitrary 
*         x1 and y1 inside the rectangle.
*/

double inv_squares(x,y,dx,z,zx,zy,zxy,x1,y1)
double x,y,dx,z,zx,zy,zxy,x1,y1;
{
    double tol;
    double s1,s2,s3,s4,xlow,xhigh,ylow,yhigh, sum;

    tol=1e-6*dx;
    if(x1< x-tol || x1>x+dx+tol || y1<y-tol || y1>y+dx+tol) 
		exit(1);

    xlow=x1-x;
    xhigh=x+dx-x1;
    ylow=y1-y;
    yhigh=y+dx-y1;

    if(xlow< -tol || xhigh< -tol || ylow < -tol || yhigh < -tol ) 
		exit(1);

    if (fabs(xlow) < tol) 
		return z+ylow*(zy-z)/dx;
    if (fabs(ylow) < tol) 
		return z+xlow*(zx-z)/dx;
    if (fabs(xhigh) < tol) 
		return zx+ylow*(zxy-zx)/dx;
    if (fabs(yhigh) < tol) 
		return zy+xlow*(zxy-zx)/dx;

    s1=1./(xlow*ylow);
    s2=1./(xhigh*ylow);
    s3=1./(xlow*yhigh);
    s4=1./(xhigh*yhigh);

    sum=s1+s2+s3+s4;
    s1=s1/sum;
    s2=s2/sum;
    s3=s3/sum;
    s4=s4/sum;

    return z*s1+zx*s2+zy*s3+zxy*s4;

}
