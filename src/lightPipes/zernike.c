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

/*  lightPipes's Factorial function  */
double lp_factorial(int n)
{
    double product;

    if (n<0) {
        printf("lp_factorial: argument is negative, exiting \n");
        exit(1);
    }
    if (n==0) return 1.;
    else{ 
        product =1;
        while(n>=1){
            product *= n;
            --n;
        }
        return product;
    }
}

/***************************************************************/
/* Zernike polynomial 

   +-m
  R    as in Born and Volf, p. 465, sixth edition, Pergamon
     n

The implementation have not been optimized for speed.
 
*/

double Zernike(int n,int m,double rho,double phi)
{
    int s, int_sign, mm, ncheck, ind;
    double sum, product, lp_factorial();

    if(n<0){
        printf("Zernike: n must be >0; \
				 |m| must be less or equal than n\n\
				 if n is odd then m must be odd,\n\
				 if n is even then m must be even\n");

        exit(1);
    }

    ind=0;
    for(ncheck=n; ncheck>=-n; ncheck -= 2){
        if (ncheck == m) ind=1;
    }

    if(ind == 0){
        printf("Zernike: n must be >0; \
				 |m| must be less or equal than n\n\
				 if n is odd then m must be odd,\n\
				 if n is even then m must be even\n");

        exit(1);
    }

    mm=(int) fabs(m);
    sum=0;
    int_sign=1;
    for (s=0; s<= (int)((n-mm)/2); s++){
        if(n-2*s != 0) product=pow( rho, (double)(n-2*s) );
        else product =1;
        product *= lp_factorial(n-s)*int_sign;
        product /= lp_factorial(s)*lp_factorial(((n+mm)/2)-s)*lp_factorial(((n-mm)/2)-s);
        sum += product;
        int_sign = -int_sign;
    }

    if(m<=0) return sum*cos(m*phi);
    else return sum*sin(m*phi); 
}

void Zer(int n,int m,double R,double A,FIELD field)
{ 
    int i,j,n2; long ik;
    double x, y, dx, fi, cab, sab, cc, cc1,rho,phi;

    n2=field.number/2;

    dx=field.size/field.number;

    ik=0;

    for (i=1;i<=field.number; i++){
        x=(i-n2-1)*dx;
        for (j=1;j<=field.number; j++){
            y=(j-n2-1)*dx;
            rho=sqrt((x*x+y*y)/(R*R));

            phi=phase(y,x);

            fi= A*Zernike(n,m,rho,phi);
            cab=cos(fi); sab=sin(fi);

            cc=field.real[ik]*cab-field.imaginary[ik]*sab;
            cc1=field.real[ik]*sab+field.imaginary[ik]*cab;

            field.real[ik]=cc;
			field.imaginary[ik]=cc1;
            ik++;
        }

    }
}

FIELD lp_zernike(FIELD field, int n, int m, double R, double AA )
{
    Zer( n,m,R,AA,field ); 
	return ( field );
}

