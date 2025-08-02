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

int construct_Emitter_int(char *fname) 
{
	int   m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4, 
		m5, n5, minlhs=1, maxlhs=1, minrhs=4, maxrhs=4;
	double *outptr=NULL; Emitter emt;
	
	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "i", &m1, &n1, &l1);  
	GetRhsVar(2, "d", &m2, &n2, &l2);  
	GetRhsVar(3, "d", &m3, &n3, &l3);  
	GetRhsVar(4, "d", &m4, &n4, &l4);  

	emt = construct_Emitter( *istk(l1), *stk(l2), *stk(l3), *stk(l4) );
	outptr = Emitter2array(emt);
	
	m5=1; n5=4;
	CreateVarFromPtr(5, "d", &m5, &n5, &outptr);  
	LhsVar(1) = 5;
	return 0;
}

int plane_wave_emit_int(char *fname) 
{
	int   m1, n1, l1, m2, n2, l2, m3, n3, number, number2,
		minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
	double *outptr=NULL; Emitter Emt;
	FIELD field; WaveFront WF;
	WavefrontHeader WfH; 

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	if ( (m1!=1)||( n1 != 18 ) )  {
		sciprint("Error: first arguments must be diffractive_wavefront_header\r\n");
		return 0;
	}
	WfH=array2WavefrontHeader(stk(l1));

	GetRhsVar(2, "d", &m2, &n2, &l2); 
	if ( (m2!=1)||(n2 != 4 ) )  {
		sciprint("Error: second arguments must be plane wave emitter\r\n");
		return 0;
	}
	Emt=construct_Emitter(0,*stk(l2+1),*stk(l2+2),*stk(l2+3));

	WF=wave_emit(Emt,WfH);

	number=WF.WfH.axes_x;
	number2=WF.WfH.axes_y;

	field=create_field(number, number2);
	WaveFront_FIELD(WF,&field);

	outptr = field2array( field );
	
	m3=1; n3= 2*number*number2+30;
	CreateVarFromPtr(3, "d", &m3, &n3, &outptr);   
	LhsVar(1) = 3;
	return 0;
}

int spherical_wave_emit_int(char *fname) 
{
	int   m1, n1, l1, m2, n2, l2, m3, n3, number, number2,
		minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
	double *outptr=NULL; Emitter Emt;
	FIELD field; WaveFront WF;
	WavefrontHeader WfH; 

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	if ( (m1!=1)||( n1 != 18 ) )  {
		sciprint("Error: first arguments must be diffractive_wavefront_header\r\n");
		return 0;
	}
	WfH=array2WavefrontHeader(stk(l1));

	GetRhsVar(2, "d", &m2, &n2, &l2); 
	if ( (m2!=1)||(n2 != 4 ) )  {
		sciprint("Error: second arguments must be spherical wave emitter\r\n");
		return 0;
	}
	Emt=construct_Emitter(1,*stk(l2+1),*stk(l2+2),*stk(l2+3));

	WF=wave_emit(Emt,WfH);

	number=WF.WfH.axes_x;
	number2=WF.WfH.axes_y;

	field=create_field(number, number2);
	WaveFront_FIELD(WF,&field);

	outptr = field2array( field );
	
	m3=1; n3= 2*number*number2+30;
	CreateVarFromPtr(3, "d", &m3, &n3, &outptr);   
	LhsVar(1) = 3;
	return 0;
}
