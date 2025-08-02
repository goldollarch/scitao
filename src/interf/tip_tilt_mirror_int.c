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

int create_ideal_ttm_int(char *fname) 
{
	int  m1, n1, l1, m2, n2, l2, m3, n3,
		minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
	double *outptr=NULL; APERTURE AP;
	TipTiltMirror ttm;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	if ( m1!=1|| n1!=17 )  {
		sciprint("Error: first arguments must be aperture \r\n");
		return 0;
	}
	AP=array2APERTURE( stk(l1) );

	GetRhsVar(2, "d", &m2, &n2, &l2); 
	if ( *stk(l2)<0 )  {
		sciprint("Error: arguments angular_velocity must be positive\r\n");
		return 0;
	}

	ttm = create_TipTiltMirror( AP,*stk(l2) );
	outptr = TipTiltMirror2array(ttm);

	m3=1; n3=22;
	CreateVarFromPtr(3, "d", &m3, &n3, &outptr);  
    LhsVar(1) = 3;
	return 0;
}

int set_ttm_timestamp_int(char *fname) 
{
	int  m1, n1, l1, m2, n2, l2,
		minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
	double *outptr=NULL;
	TipTiltMirror ttm;
	
	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1); 
	if ( ( m1!=1 ) || ( n1 != 22 ) )  {
		sciprint("Error: first arguments must be tip_tilt mirror\r\n");
		return 0;
	}
	ttm=array2TipTiltMirror(stk(l1));

	GetRhsVar(2, "d", &m2, &n2, &l2); 

	set_ttm_timestamp(&ttm,*stk(l2));
	
	outptr = TipTiltMirror2array(ttm);
	
	CreateVarFromPtr(3, "d", &m1, &n1, &outptr);   
	LhsVar(1) = 3;
	return 0;
}

int set_ttm_commands_vector_int(char *fname) 
{
	int  m1, n1, l1, m2, n2, l2,
		minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
	double *outptr=NULL;
	TipTiltMirror ttm;
	
	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1); 
	if ( ( m1!=1 ) || ( n1 != 22 ) )  {
		sciprint("Error: first arguments must be tip_tilt mirror\r\n");
		return 0;
	}
	ttm=array2TipTiltMirror(stk(l1));

	GetRhsVar(2, "d", &m2, &n2, &l2);
	if ( ( m2!=1 )||( n2 != 3 ) )  {
		sciprint("Error: first arguments must be commands vector \r\n");
		return 0;
	}

	set_command_vector( &ttm, *stk(l2), *stk(l2+1), *stk(l2+2) );
	
	outptr = TipTiltMirror2array(ttm);
	CreateVarFromPtr(3, "d", &m1, &n1, &outptr);   
	LhsVar(1) = 3;
	return 0;

}

int TipTiltMirror_update_int(char *fname) 
{
	int  m1, n1, l1, m2, n2, l2, m3, n3, l3,
		minlhs=1, maxlhs=1, minrhs=3, maxrhs=3;
	ThreeVector ttm_command;
	TipTiltMirror ttm;
	double *outptr=NULL;
	
	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1); 
	if ( ( m1!=1 ) || ( n1 != 22 ) )  {
		sciprint("Error: first arguments must be tip_tilt mirror\r\n");
		return 0;
	}
	ttm=array2TipTiltMirror(stk(l1));

	GetRhsVar(2, "d", &m2, &n2, &l2);
	if ( ( m2!=1 )||( n2 != 3 ) )  {
		sciprint("Error: first arguments must be commands vector \r\n");
		return 0;
	}
	ttm_command = construct_ThreeVector(*stk(l2), *stk(l2+1), *stk(l2+2));

	GetRhsVar(3, "d", &m3, &n3, &l3); 
	
	TipTiltMirror_updata(&ttm,ttm_command,*stk(l3));

	outptr = TipTiltMirror2array(ttm);
	CreateVarFromPtr(4, "d", &m1, &n1, &outptr);   
	LhsVar(1) = 4;

	return 0;
}

int ideal_ttm_transform_int(char *fname) 
{
	FIELD field;
	TipTiltMirror TTM;
	WaveFront WF_in, WF_out;
	
	int  m1, n1, l1, m2, n2, l2, number, number2,
		minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
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
	if ( (m2!=1)||( n2 != 22 ) )  {
		sciprint("Error: second arguments must be tip_tilt mirror\r\n");
		return 0;
	}
	TTM=array2TipTiltMirror(stk(l2));

	WF_out =TipTiltMirror_transform(&TTM, WF_in);
	WaveFront_FIELD(WF_out,&field);

	outptr = field2array( field );
	CreateVarFromPtr(3, "d", &m1, &n1, &outptr); 
	LhsVar(1) = 3;
	return 0;
}

