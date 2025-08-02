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

extern int fresnl();

double *lp_forward(  FIELD field, double z, int new_number, double new_size )
{
    double old_size, dx_new, dx_old, x_new, y_new,old_number, 
		on21, nn21,P1, P2, P3, P4, R22; 

    double fc1, fs1, fc2, fs2, fc3, fs3, fc4, fs4, fr, fi;

    double c4c1, c2c3, c4s1, s4c1, s2c3, c2s1, s4c3, s2c1, c4s3, s2s3,
		s2s1, c2s3, s4s1, c4c3, s4s3, c2c1;

	double *fieldimaginary, fieldreal;

    int i,j,i_old, j_old, io, jo, dum;
    long  kk, kk1,array_space;
	
	double *outptr;

    fs1=0.;	fc1=0.;	fs2=0.;	fc2=0.;
	fs3=0.;	fc3=0.;	fs4=0.;	fc4=0.;

	if (new_size==0) new_size=field.size;
	if (new_number==0) new_number=field.number;

	array_space=new_number*new_number;
	outptr = (double *) calloc( 2*array_space+30, sizeof(double) );

	field_header_array(field,outptr);

    R22=sqrt(1./(2.*field.lambda*z));

    old_size=field.size;
    old_number=field.number;
        
    dx_new=new_size/(new_number-1.); 
    dx_old=old_size/(old_number-1.);
   
    field.number = new_number;
    field.size = new_size;

    on21=old_number/2+1;
    nn21=new_number/2+1;

    fieldimaginary = (double *) calloc( array_space, sizeof(double) );

    kk=0;
    for (i=1; i<=field.number; i++) 
	{
		x_new=(i-nn21)*dx_new;     
        for (j=1; j<=field.number; j++) 
		{
			y_new=(j-nn21)*dx_new;
			fieldimaginary[kk]=0.;
			fieldreal = 0.;

			kk1=0;
			for (i_old=1; i_old<=old_number; i_old++)
			{
				io=i_old-on21;
				for (j_old=1; j_old<=old_number; j_old++)
				{
					jo=j_old-on21;

					P1=R22*(2*(dx_old*io-x_new)+dx_old);
					P2=R22*(2*(dx_old*jo-y_new)-dx_old);
					P3=R22*(2*(dx_old*io-x_new)-dx_old);
					P4=R22*(2*(dx_old*jo-y_new)+dx_old);

					dum=fresnl(P1,&fs1, &fc1);
					dum=fresnl(P2,&fs2, &fc2);
					dum=fresnl(P3,&fs3, &fc3);
					dum=fresnl(P4,&fs4, &fc4);

					fr=0.5*field.real[kk1];
					fi=0.5*field.imaginary[kk1];

					c4c1=fc4*fc1;  c2s3=fc2*fs3;	c4s1=fc4*fs1;
					s4c1=fs4*fc1;  s2c3=fs2*fc3;	c2s1=fc2*fs1;
					s4c3=fs4*fc3;  s2c1=fs2*fc1;	c4s3=fc4*fs3;

					fieldreal += fr*(c2s3+c4s1+s4c1+s2c3-c2s1-s4c3-s2c1-c4s3);

					s2s3=fs2*fs3;	s2s1=fs2*fs1;
					c2c3=fc2*fc3;	s4s1=fs4*fs1;
					c4c3=fc4*fc3;	c4c1=fc4*fc1;
					s4s3=fs4*fs3;	c2c1=fc2*fc1;

					fieldreal += fi*(-s2s3+s2s1+c2c3-s4s1-c4c3+c4c1+s4s3-c2c1);

					fieldimaginary[kk] +=
						fr*(-c4c1+s2s3+c4c3-s4s3+c2c1-s2s1+s4s1-c2c3);
					fieldimaginary[kk] +=
						fi*(c2s3+s2c3+c4s1+s4c1-c4s3-s4c3-c2s1-s2c1); 

					kk1 ++;
				}

			}

			outptr[kk+30]=fieldreal;
			kk ++;
		}
	}

	kk=0;
    for (i=1; i<=new_number; i++)
	{
		for (j=1; j<=new_number; j++)
		{
			outptr[kk+array_space+30] = fieldimaginary[kk];
			kk ++;
		}
	}

	free(fieldimaginary);

	return outptr;
}
