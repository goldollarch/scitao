/*
scitao - A toolbox based on Scilab/Scicos environment for the simulation 
of wave optics, especially for the simulation of adaptive optics .

Copyright (c) 2005-2006 IAPCM, Beijing, China.  Written by
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
with this program; if not, write to the Free Software Foundation,
Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
*/

#include "int_optics.h"

//////////////////////////////////////////////////////////////////////////////////

int create_pixel_array_int(char *fname) 
{
	int  m1, n1, l1, m2, n2, l2, m3, n3, 
		minlhs=1, maxlhs=1, minrhs=1, maxrhs=2;
	PixelArray tmp_PixArr;
	int x_axes, y_axes;
	double *outptr=NULL;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "i", &m1, &n1, &l1); 
	if ( *istk(l1) < 0 )  {
		sciprint("Error: x_axes arguments must be positive\r\n");
		return 0;
	}
	x_axes=(*istk(l1)); y_axes=x_axes;

	if (Rhs>1) {
		GetRhsVar(2, "i", &m2, &n2, &l2); 
		if ( *istk(l2) < 0 )  {
			sciprint("Error: y_axes arguments must be positive\r\n");
			return 0;
		}
		y_axes=(*istk(l2));
	}

	tmp_PixArr=create_PixelArray(x_axes,y_axes,0);
	outptr = PixelArray2array(tmp_PixArr);

	m3=1; n3=x_axes*y_axes+2;
	CreateVarFromPtr(3, "d", &m3, &n3, &outptr);   

	LhsVar(1) = 3;
	return 0;
}

int set_pixel_array_data_int(char *fname) 
{
	int  m1, n1, l1, m2, n2, l2,
        minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
	PixelArray tmp_PixArr;
	int x_axes, y_axes;
	double *outptr=NULL;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);  

	x_axes=(int)(*stk(l1));
	y_axes=(int)(*stk(l1+1));

	if (  (m1 !=1)||(n1 !=(x_axes*y_axes+2)) )  {
		sciprint("Error: first arguments must be pixel_array\r\n");
		return 0;
	}
	tmp_PixArr=array2PixelArray(stk(l1));

	GetRhsVar(2, "d", &m2, &n2, &l2); 
	if (  ( m2 != 1)||( n2 != x_axes*y_axes ) )  {
		sciprint("Error: second arguments must be %d elements\r\n",x_axes*y_axes);
		return 0;
	}
	set_PixelArray_data( &tmp_PixArr, stk(l2) );

	outptr = PixelArray2array(tmp_PixArr);
	CreateVarFromPtr(3, "d", &m1, &n1, &outptr);  
	LhsVar(1) = 3;
	return 0;
}

int multi_pixarr_int(char *fname) 
{
	int  m1, n1, l1, m2, n2, l2, x_axes, y_axes, nelem, i,
        minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
	double gain, *outptr=NULL;
	PixelArray tmp_PixArr;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);  

	x_axes=(int)(*stk(l1));	y_axes=(int)(*stk(l1+1));
	nelem=x_axes*y_axes;
	if (  (m1 !=1)||(n1 !=(nelem+2)) )  {
		sciprint("Error: first arguments must be pixel_array\r\n");
		return 0;
	}
	tmp_PixArr=array2PixelArray(stk(l1));

	GetRhsVar(2, "d", &m2, &n2, &l2);
	gain = *stk(l2);

	for(i=0;i<nelem;i++) tmp_PixArr.pixeldata[i] *= gain;

	outptr = PixelArray2array(tmp_PixArr);
	CreateVarFromPtr(3, "d", &m1, &n1, &outptr);  
	LhsVar(1) = 3;
	return 0;
}

int pixarr_data_int(char *fname) 
{
	int  m1, n1, l1, m2, n2, l2,
        minlhs=1, maxlhs=1, minrhs=1, maxrhs=1;
	PixelArray tmp_PixArr;
	int x_axes, y_axes;
	double *outptr=NULL;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);  

	x_axes=(int)(*stk(l1));
	y_axes=(int)(*stk(l1+1));

	if (  (m1 !=1)||(n1 !=(x_axes*y_axes+2)) )  {
		sciprint("Error: first arguments must be pixel_array\r\n");
		return 0;
	}
	tmp_PixArr=array2PixelArray(stk(l1));

	m2=x_axes;n2=x_axes;
	outptr = tmp_PixArr.pixeldata;
	CreateVarFromPtr(2, "d", &m2, &n2, &outptr);  
	LhsVar(1) = 2;
	return 0;
}
