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

int geom_prop_int(char *fname) 
{
	FIELD field;
	WaveFront  WF_in, WF_out;

	int   m1, n1, l1, m2, n2, l2, number, number2,
		minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
	double z, *outptr=NULL;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	number=(int)(*stk(l1)); number2=(int)(*stk(l1+10));
	if ( (m1!=1)||(n1 != (2*number*number2+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}

	field=array2field( stk(l1) );
	WF_in=FIELD_WaveFront(field);

	GetRhsVar(2, "d", &m2, &n2, &l2); 
	z= *stk(l2);

	WF_out=geometric_propagate_transform(z,WF_in);
	WaveFront_FIELD( WF_out, &field );
	outptr = field2array( field );
    
	CreateVarFromPtr(3, "d", &m1, &n1, &outptr);   
	LhsVar(1) = 3;
	return 0;
}

int near_ang_int(char *fname) 
{
	FIELD field;
	WaveFront  WF_in, WF_out;

	int   m1, n1, l1, m2, n2, l2, number, number2,
		minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
	double z, *outptr=NULL;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	number=(int)(*stk(l1)); number2=(int)(*stk(l1+10));
	if ( (m1!=1)||(n1 != (2*number*number2+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}

	field=array2field( stk(l1) );
	WF_in=FIELD_WaveFront(field);

	GetRhsVar(2, "d", &m2, &n2, &l2); 
	z= *stk(l2);

	WF_out=near_field_angular_transform(z,WF_in);
	WaveFront_FIELD( WF_out, &field );
	outptr = field2array( field );
    
	CreateVarFromPtr(3, "d", &m1, &n1, &outptr);   
	LhsVar(1) = 3;
	return 0;
}

int near_fresnel_int(char *fname) 
{
	FIELD field;
	WaveFront  WF_in, WF_out;

	int   m1, n1, l1, m2, n2, l2, number, number2,
		minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
	double z, *outptr=NULL;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	number=(int)(*stk(l1)); number2=(int)(*stk(l1+10));
	if ( (m1!=1)||(n1 != (2*number*number2+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}

	field=array2field( stk(l1) );
	WF_in=FIELD_WaveFront(field);

	GetRhsVar(2, "d", &m2, &n2, &l2); 
	z= *stk(l2);

	WF_out=near_field_fresnel_transform(z,WF_in);
	WaveFront_FIELD( WF_out, &field );
	outptr = field2array( field );
    
	CreateVarFromPtr(3, "d", &m1, &n1, &outptr);   
	LhsVar(1) = 3;
	return 0;
}

int far_fresnel_int(char *fname) 
{
	FIELD field;
	WaveFront  WF_in, WF_out;

	int   m1, n1, l1, m2, n2, l2, number, number2,
		minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
	double z, *outptr=NULL;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	number=(int)(*stk(l1)); number2=(int)(*stk(l1+10));
	if ( (m1!=1)||(n1 != (2*number*number2+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}

	field=array2field( stk(l1) );
	WF_in=FIELD_WaveFront(field);

	GetRhsVar(2, "d", &m2, &n2, &l2); 
	z= *stk(l2);

	WF_out=far_field_fresnel_transform(z,WF_in);
	WaveFront_FIELD( WF_out, &field );
	outptr = field2array( field );
    
	CreateVarFromPtr(3, "d", &m1, &n1, &outptr);   
	LhsVar(1) = 3;
	return 0;
}

int far_fraunhoffer_int(char *fname) 
{
	FIELD field;
	WaveFront  WF_in, WF_out;

	int   m1, n1, l1, m2, n2, l2, number, number2,
		minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
	double z, *outptr=NULL;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	number=(int)(*stk(l1)); number2=(int)(*stk(l1+10));
	if ( (m1!=1)||(n1 != (2*number*number2+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}

	field=array2field( stk(l1) );
	WF_in=FIELD_WaveFront(field);

	GetRhsVar(2, "d", &m2, &n2, &l2); 
	z= *stk(l2);

	WF_out=far_field_fraunhoffer_transform(z,WF_in);
	WaveFront_FIELD( WF_out, &field );
	outptr = field2array( field );
    
	CreateVarFromPtr(3, "d", &m1, &n1, &outptr);   
	LhsVar(1) = 3;
	return 0;
}

int far_fresnel_GRT_int(char *fname) 
{
	FIELD field;
	WaveFront  WF_in, WF_out;
	
	int   m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4, m5, n5, l5,
		m6, n6, number, number2, final_x_axes,final_y_axes,
		minlhs=1, maxlhs=1, minrhs=4, maxrhs=5;
	double distance,final_pixel_scale, *outptr=NULL;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	number=(int)(*stk(l1)); number2=(int)(*stk(l1+10));
	if ( (m1!=1)||(n1 != (2*number*number2+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}
	field=array2field( stk(l1) );
	WF_in=FIELD_WaveFront(field);

	GetRhsVar(2, "d", &m2, &n2, &l2); 
	distance= *stk(l2);

	GetRhsVar(3, "d", &m3, &n3, &l3); 
	if ( *stk(l3) < 0 )  {
		sciprint("Error: arguments final_pixel_scale must be position interger\r\n");
		return 0;
	}
	final_pixel_scale= *stk(l3);

	GetRhsVar(4, "i", &m4, &n4, &l4);
	if ( *istk(l4) < 0 )  {
		sciprint("Error: arguments final_x_axes must be position interger\r\n");
		return 0;
	}
	final_x_axes= *istk(l4);
	final_y_axes= final_x_axes;

	if (Rhs>4) {
		GetRhsVar(5, "i", &m5, &n5, &l5); 
		if ( *istk(l5) < 0 )  {
			sciprint("Error: arguments final_y_axes must be position interger\r\n");
			return 0;
		}
		final_y_axes= *istk(l5);
	}

	WF_out=far_fresnel_goertzel_reinsch_transform( distance,
		   final_pixel_scale,final_x_axes,final_y_axes,WF_in);

	field=create_field(final_x_axes, final_y_axes);
	WaveFront_FIELD( WF_out, &field );
	outptr = field2array( field );

	m6=1; n6=2*final_x_axes*final_y_axes+30;
	CreateVarFromPtr(6, "d", &m6, &n6, &outptr);   

	LhsVar(1) = 6;
	return 0;
}

int far_fraunhoffer_GRT_int(char *fname) 
{
	FIELD field;
	WaveFront  WF_in, WF_out;
	
	int   m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4, m5, n5, l5,
		m6, n6, number, number2, final_x_axes,final_y_axes,
		minlhs=1, maxlhs=1, minrhs=4, maxrhs=5;
	double distance,final_pixel_scale,*outptr=NULL;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	number=(int)(*stk(l1)); number2=(int)(*stk(l1+10));
	if ( (m1!=1)||(n1 != (2*number*number2+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}
	field=array2field( stk(l1) );
	WF_in=FIELD_WaveFront(field);

	GetRhsVar(2, "d", &m2, &n2, &l2); 
	distance= *stk(l2);

	GetRhsVar(3, "d", &m3, &n3, &l3); 
	if ( *stk(l3) < 0 )  {
		sciprint("Error: arguments final_pixel_scale must be position interger\r\n");
		return 0;
	}
	final_pixel_scale= *stk(l3);

	GetRhsVar(4, "i", &m4, &n4, &l4);
	if ( *istk(l4) < 0 )  {
		sciprint("Error: arguments final_x_axes must be position interger\r\n");
		return 0;
	}
	final_x_axes= *istk(l4);
	final_y_axes= final_x_axes;

	if (Rhs>4) {
		GetRhsVar(5, "i", &m5, &n5, &l5); 
		if ( *istk(l5) < 0 )  {
			sciprint("Error: arguments final_y_axes must be position interger\r\n");
			return 0;
		}
		final_y_axes= *istk(l5);
	}

	WF_out=far_fraunhoffer_goertzel_reinsch_transform( distance,
		   final_pixel_scale,final_x_axes,final_y_axes,WF_in);

	field=create_field(final_x_axes, final_y_axes);
	WaveFront_FIELD( WF_out, &field );
	outptr = field2array( field );

	m6=1; n6=2*final_x_axes*final_y_axes+30;
	CreateVarFromPtr(6, "d", &m6, &n6, &outptr);   

	LhsVar(1) = 6;
	return 0;
}
