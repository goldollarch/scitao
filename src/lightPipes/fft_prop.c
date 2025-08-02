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

extern void fft3();

void forvard(double zz, FIELD field)
{ 
    int i,j, ii, ij, n12;
    long ik, ir;
    double z,z1,cc;
    double sw, sw1, bus, abus, pi2, cab, sab;
    pi2=2.*3.141592654;
    z=fabs(zz);
 
    ik=0;
    ii=ij=1;
    for (i=1;i<=field.number; i++){
        for (j=1;j<=field.number; j++){
            field.real[ik] = field.real[ik]*ii*ij;
            field.imaginary[ik] =field.imaginary[ik]*ii*ij; 
            ik++;
            ij=-ij;
        }
        ii=-ii;
    }

    if(zz>=0.) fft3(1,field);
    else fft3(-1,field);

    if(zz >= 0.){
		z1=z*field.lambda/2.;
		n12=field.number/2;
		ik=0;
		for (i=1;i<=field.number; i++){
			for (j=1;j<=field.number; j++){ 
				sw=((i-n12-1)/field.size);
				sw *= sw;
				sw1=((j-n12-1)/field.size);
				sw1 *= sw1;
				sw += sw1; 
				bus=z1*sw;
				ir = (long) bus;
				abus=pi2*(ir- bus);

				cab=cos(abus);
				sab=sin(abus);
				cc=field.real[ik]*cab-field.imaginary[ik]*sab;
				field.imaginary[ik]=field.real[ik]*sab+field.imaginary[ik]*cab;
				field.real[ik]=cc;
				ik++;
			}
		}
	}

    else { 
      z1=z*field.lambda/2.;
      n12=field.number/2;
      ik=0;
      for (i=1;i<=field.number; i++){
		  for (j=1;j<=field.number; j++){ 
            sw=((i-n12-1)/field.size);
            sw *= sw;
            sw1=((j-n12-1)/field.size);
            sw1 *= sw1;
            sw += sw1; 
            bus=z1*sw;
            ir = (long) bus;
            abus=pi2*(ir- bus);

            cab=cos(abus);
            sab=sin(abus);
            cc=field.real[ik]*cab + field.imaginary[ik]*sab;
            field.imaginary[ik]= field.imaginary[ik]*cab-field.real[ik]*sab;
            field.real[ik]=cc;
            ik++;
		  }
	  }
      
	}
  
    if(zz >= 0.) fft3(-1,field);
    else fft3(1,field);   
    
	ik=0;
    ii=ij=1;
    for (i=1;i<=field.number; i++){
        for (j=1;j<=field.number; j++){
            field.real[ik] = field.real[ik]*ii*ij;
            field.imaginary[ik] =field.imaginary[ik]*ii*ij;
     
            ik++;
            ij=-ij;
        }
        ii=-ii;
    }
    
}
