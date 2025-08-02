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

int lenslet_array_int(char *fname) 
{
	int  m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4,  m5, n5, l5, 
		m6, n6, l6, m7, n7, minlhs=1, maxlhs=1, minrhs=6, maxrhs=6;
	int lenslet_x_axes,lenslet_y_axes, pix_per_lenslet, pix_per_xform;
	double flength, lnslt_pitch, *outptr=NULL;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "i", &m1, &n1, &l1);  
	if ( *istk(l1) < 0 )  {
		sciprint("Error: arguments lenslet_x_axes must be positive integer\r\n");
		return 0;
	}
	lenslet_x_axes= *istk(l1);

	GetRhsVar(2, "i", &m2, &n2, &l2); 
	if ( *istk(l2) < 0 )  {
		sciprint("Error: arguments lenslet_y_axes must be positive integer\r\n");
		return 0;
	}
	lenslet_y_axes= *istk(l2);
	
	GetRhsVar(3, "d", &m3, &n3, &l3); 
	if ( *stk(l3) < 0 )  {
		sciprint("Error: arguments flength must be positive\r\n");
		return 0;
	}
	flength= *stk(l3);
	
	GetRhsVar(4, "d", &m4, &n4, &l4); 
	if ( *stk(l4) < 0 )  {
		sciprint("Error: arguments lnslt_pitch must be position interger\r\n");
		return 0;
	}
	lnslt_pitch= *stk(l4);

	GetRhsVar(5, "i", &m5, &n5, &l5); 
	if ( *istk(l5) < 0 )  {
		sciprint("Error: arguments pix_per_lenslet must be positive integer\r\n");
		return 0;
	}
	pix_per_lenslet= *istk(l5);

	GetRhsVar(6, "i", &m6, &n6, &l6); 
	if ( *istk(l6) < 0 )  {
		sciprint("Error: arguments pix_per_xform must be positive integer\r\n");
		return 0;
	}
	pix_per_xform= *istk(l6);

	outptr = (double*)calloc( 6, sizeof(double));
	outptr[0]=lenslet_x_axes; outptr[1]=lenslet_y_axes;
	outptr[2]=flength; outptr[3]=lnslt_pitch;
	outptr[4]=pix_per_lenslet; outptr[5]=pix_per_xform;

	m7=1; n7=6;
	CreateVarFromPtr(7, "d", &m7, &n7, &outptr);  
	LhsVar(1) = 7;
	return 0;
}

int lenslet_array_transform_int(char *fname) 
{
	FIELD field;
	LensletArray lnslt_array;
	WaveFront  WF_in, WF_out;

	int  m1, n1, l1, m2, n2, l2, m3, n3,
		minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
	double *outptr=NULL;
	int number,number2;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);  
    number=(int)(*stk(l1)); number2=(int)(*stk(l1+10));
	if ( (m1!=1)||(n1 != (2*number*number2+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}
	field = array2field( stk(l1) );
	WF_in = FIELD_WaveFront(field);
	
	GetRhsVar(2, "d", &m2, &n2, &l2); 
	if ( m2!=1|| n2!=6 )  {
		sciprint("Error: second arguments must be lenslet array\r\n");
		return 0;
	}
	lnslt_array=create_LensletArray( (int)(*stk(l2)), (int)(*stk(l2+1)),
		*stk(l2+2), *stk(l2+3), (int)(*stk(l2+4)), (int)(*stk(l2+5)));

	WF_out =LensletArray_transform(lnslt_array,WF_in);
	WaveFront_FIELD(WF_out,&field);

	outptr = field2array( field );

	number=(int)outptr[0];  number2=(int)outptr[10];
	m3=1; n3=2*number*number2+30;
	
	CreateVarFromPtr(3, "d", &m3, &n3, &outptr);   
	LhsVar(1) = 3;
	return 0;
}

int create_shcentroids_int(char *fname) 
{
	FIELD field;
	WaveFront WF_in;
	SHartmannCentroids shcentroid;
	int lnslt_x, lnslt_y, number, number2, m1, n1, l1, m2, n2, l2,
		m3, n3, minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
	double *outptr=NULL;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1); 
    number=(int)(*stk(l1)); number2=(int)(*stk(l1+10));
	if ( (m1!=1)||(n1 != (2*number*number2+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}
	field = array2field( stk(l1) );
	WF_in = FIELD_WaveFront(field);

	GetRhsVar(2, "d", &m2, &n2, &l2); 
	if ( m2!=1|| n2!=6 )  {
		sciprint("Error: second arguments must be lenslet array\r\n");
		return 0;
	}
	lnslt_x=(int)(*stk(l2)); 
	lnslt_y=(int)(*stk(l2+1));

	shcentroid = create_SHartmannCentroids ( lnslt_x, lnslt_y, WF_in );
	outptr = SHartmannCentroids2array( shcentroid );

	number=shcentroid.pixarr_x_axes; 
	number2=shcentroid.pixarr_y_axes;
	m3=1; n3=number*number2+2;
	CreateVarFromPtr(3, "d", &m3, &n3, &outptr);   

	LhsVar(1) = 3;
	return 0;
}

int shc_fits_int(char *fname) 
{
	int  m1, n1, l1, m2, n2, l2, m3, n3, l3, pixarr_x, pixarr_y,
		minlhs=1, maxlhs=1, minrhs=2, maxrhs=3;
	SHartmannCentroids shcentroids;
	double *outptr=NULL;
	double time=-1;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);
	pixarr_x=(int)(*stk(l1)); pixarr_y=(int)(*stk(l1+1));
	if ( (m1!=1)||(n1 != pixarr_x*pixarr_y+2) )  {
		sciprint("Error: first arguments must be SHartmann Centroids \r\n");
		return 0;
	}
	shcentroids = array2SHartmannCentroids( stk(l1) );

	GetRhsVar(2, "c", &m2, &n2, &l2); 

	if (Rhs>2) {
		GetRhsVar(3, "d", &m3, &n3, &l3); 
		time= *stk(l3);
	}

	write_SHartmannCentroids(shcentroids,cstk(l2),time);

	LhsVar(1) = 1;
	return 0;
}
