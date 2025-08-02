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

int diffractive_wavefront_header_int(char *fname) 
{
	int m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4, m5, n5, l5,	
		m6=1, n6=18, minlhs=1, maxlhs=1, minrhs=2, maxrhs=5,
		number=256, number2=256;
	ThreeFrame TF=default_ThreeFrame();
	WavefrontHeader WfH;
	double *outptr;

	CheckRhs(minrhs,maxrhs);
	CheckLhs(minlhs,maxlhs);

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	if ( *stk(l1) <= 0 )  {
		sciprint("Error: wavelength must be positive\r\n");
		return 0;
	}

	GetRhsVar(2, "d", &m2, &n2, &l2);
	if ( *stk(l2) <= 0 )  {
		sciprint("Error: pixscale must be positive\r\n");
		return 0;
	}

	if (Rhs>2) {
		GetRhsVar(3, "i", &m3, &n3, &l3); 
		number= *istk(l3);
		if ( number <= 0 )  {
			sciprint("Error: the dimension of axes must be positive\r\n");
			return 0;
		}
		number2=number;
	}

	if (Rhs>3) {
		GetRhsVar(4, "i", &m4, &n4, &l4);
		number2= *istk(l4);
		if ( number2 <= 0 )  {
			sciprint("Error: the dimension of  y_axes must be positive\r\n");
			return 0;
		}
	}
	
	if (Rhs>4) {
		GetRhsVar(5, "d", &m5, &n5, &l5);  
		if ( (m5!=1)||(n5 != 12 ) )  {
			sciprint("Error: the last of arguments must be three frame\r\n");
			return 0;
		}
		TF=array2ThreeFrame(stk(l5));
	}

	WfH=construct_WavefrontHeader( TF, number, number2, *stk(l1), *stk(l2) );	

	outptr = WavefrontHeader2array(WfH);

	CreateVarFromPtr(6, "d", &m6, &n6, &outptr);   
    LhsVar(1) = 6;
	return 0;
}

int set_dwf_timestamp_int(char *fname) 
{
	int  m1, n1, l1, m2, n2, l2, number, number2,
		minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
	WaveFront WF_in; FIELD field;
	double *outptr=NULL;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	number=(int)(*stk(l1)); number2=(int)(*stk(l1+10));
	if ( (m1!=1)||(n1 != (2*number*number2+30) ) )  {
		sciprint("Error: first arguments must be diffractive_wavefront\r\n");
		return 0;
	}
	field = array2field( stk(l1) );
	WF_in = FIELD_WaveFront(field);

	GetRhsVar(2, "d", &m2, &n2, &l2); 
	if ( *stk(l2) < 0 )  {
		sciprint("Error: arguments timestamp must be positive\r\n");
		return 0;
	}
	WF_in.WfH.timestamp = *stk(l2);

	WaveFront_FIELD(WF_in,&field);
	outptr = field2array( field );
	CreateVarFromPtr(3, "d", &m1, &n1, &outptr);  
	LhsVar(1) = 3;
	return 0;
}

int set_dwf_frame_int(char *fname) 
{
	int  m1, n1, l1, m2, n2, l2, number, number2,
		minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
	WaveFront WF_in;  FIELD field;
	ThreeFrame TF; double *outptr=NULL;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	number=(int)(*stk(l1)); number2=(int)(*stk(l1+10));
	if ( (m1!=1)||(n1 != (2*number*number2+30) ) )  {
		sciprint("Error: first arguments must be diffractive_wavefront\r\n");
		return 0;
	}
	field = array2field( stk(l1) );
	WF_in = FIELD_WaveFront(field);

	GetRhsVar(2, "d", &m2, &n2, &l2); 
	if ( (m2!=1)||(n2 != 12 ) )  {
		sciprint("Error: second arguments must be 3D frame\r\n");
		return 0;
	}
	TF = array2ThreeFrame( stk(l2) );
	WF_in.WfH.TF=TF;

	WaveFront_FIELD(WF_in,&field);
	outptr = field2array( field );
	CreateVarFromPtr(3, "d", &m1, &n1, &outptr);  
	LhsVar(1) = 3;
	return 0;
}

int set_dwf_direction_int(char *fname) 
{
	int  m1, n1, l1, m2, n2, l2, number, number2,
		minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
	ThreeVector TV;	FIELD field; 
	WaveFront WF_in,WF_out;
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
	if ( (m2!=1)||(n2 != 3 ) )  {
		sciprint("Error: second arguments must be 3D vector\r\n");
		return 0;
	}
	TV=construct_ThreeVector(*stk(l2),*stk(l2+1),*stk(l2+2));

	WF_out=set_WaveFront_propagation_direction(WF_in,TV);

	WaveFront_FIELD(WF_out,&field);
	outptr = field2array( field );
	CreateVarFromPtr(3, "d", &m1, &n1, &outptr);   
	LhsVar(1) = 3;
	return 0;
}

int dwf_clip_array_int(char *fname) 
{
	int  m1, n1, l1, m2, n2, l2, m3=1, n3, number, number2, nclip,
		minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
	FIELD field; WaveFront WF_in,WF_out;
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

	GetRhsVar(2, "i", &m2, &n2, &l2); 
	nclip = *istk(l2);

	WF_out=WaveFront_clip_array(WF_in,nclip);
	WaveFront_FIELD(WF_out,&field);

	outptr = field2array( field );
	number -= 2*nclip;number2 -= 2*nclip;n3 = 2*number*number2+30;
	CreateVarFromPtr(3, "d", &m3, &n3, &outptr);   
	LhsVar(1) = 3;
	return 0;
}

int set_dwf_pixel_scale_int(char *fname) 
{
	int  m1, n1, l1, m2, n2, l2, number, number2,
		minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
	FIELD field; WaveFront WF_in,WF_out;
	double pxlscl,*outptr=NULL;

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
	pxlscl = *stk(l2);

	WF_out=set_dwf_pixel_scale(WF_in,pxlscl);
	WaveFront_FIELD(WF_out,&field);

	outptr = field2array( field );
	CreateVarFromPtr(3, "d", &m1, &n1, &outptr);   
	LhsVar(1) = 3;
	return 0;
}

int dwf_fits_int(char *fname) 
{
	int  m1, n1, l1, m2, n2, l2, m3, n3, l3, number, number2,
		minlhs=1, maxlhs=1, minrhs=2, maxrhs=3;
	FIELD field; WaveFront WF_in;
	double *outptr=NULL;
	double time=-1;

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

	GetRhsVar(2, "c", &m2, &n2, &l2);  

	if (Rhs>2) {
		GetRhsVar(3, "d", &m3, &n3, &l3); 
		time= *stk(l3);
	}

	write_WaveFront_file( WF_in, cstk(l2), time );

	LhsVar(1) = 1;
	return 0;
}
