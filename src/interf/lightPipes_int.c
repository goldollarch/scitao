
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

int begin_int(char *fname) 
{
	int  m_size, n_size, l_size, m_lambda, n_lambda, l_lambda, m3, n3, l3,
		m4, n4, l4, minlhs=1, maxlhs=1, minrhs=2, maxrhs=4;
	double size, lambda, pixscale,curvature,timestamp;
	int  number, number2, m_outptr, n_outptr;
	ThreeFrame TF; FIELD field;
	double *outptr=NULL;

	CheckRhs(minrhs,maxrhs);
	CheckLhs(minlhs,maxlhs);

	GetRhsVar( 1, "d", &m_size, &n_size, &l_size ); 
	size= *stk(l_size);

	GetRhsVar( 2, "d", &m_lambda, &n_lambda, &l_lambda );  
	lambda= *stk(l_lambda);

	if ( *stk(l_size) <= 0 || *stk(l_lambda) <= 0 )  {
		sciprint("Error: arguments must be positive\r\n");
		return 0;
	}
	
	number=256; number2=256;

	if (Rhs>2) {
		GetRhsVar( 3, "i", &m3, &n3, &l3 ); 
		number=*istk( l3 ); number2=number;
		if ( number%2 != 0 ) {
			sciprint("Error: The number of grid must be even\r\n");
			return 0;
		}
	}

	if( Rhs>3 ) {
		GetRhsVar( 4, "i", &m4,  &n4, &l4 ); 
		number2=*istk( l4 );
		if ( number2%2 != 0 ) {
			sciprint("Error: The number2 of grid must be even\r\n");
			return 0;
		}
		if ( number2 != number ) {
			sciprint("Notice: the dimension of  y axes are not equal to x's,\
					 \n some lightPipes functions may not work properly.\
					 \n you 'd better let them equate !!");
		}
	}
  
	field=lp_begin( size, lambda, number, number2 );

	timestamp=0;curvature=0;pixscale=size/number;	
	TF=default_ThreeFrame(); set_FIELD_ThreeFrame( TF, &field );
	set_FIELD_WavefrontHeader( number,number2,lambda,pixscale,curvature,timestamp,&field);

	outptr = field2array( field );

	m_outptr=1; n_outptr = 2*number*number2+30;
	CreateVarFromPtr( 5, "d", &m_outptr, &n_outptr, &outptr );   

	LhsVar(1) = 5;
	return 0;
}

int circ_ap_int(char *fname) 
{
	int  m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4,
		minlhs=1, maxlhs=1, minrhs=2, maxrhs=4, number, number2;
	FIELD field_in, field_out; double *outptr=NULL;
	double x_shift=0, y_shift=0;
  
	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1); 
	number=(int)(*stk(l1)); number2=(int)(*stk(l1+10));
	if ( (m1!=1)||(n1 != (2*number*number2+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}
	field_in=array2field( stk(l1) );

	GetRhsVar(2, "d", &m2, &n2, &l2); 
	if ( *stk(l2) < 0 )  {
		sciprint("Error: arguments R must be positive\r\n");
		return 0;
	}

	if (Rhs>2)
	{
		GetRhsVar(3, "d", &m3, &n3, &l3); 
		x_shift=*stk(l3);
	}
	
	if (Rhs>3)
	{
		GetRhsVar(4, "d", &m4, &n4, &l4);
		y_shift=*stk(l4);
	}

	field_out=lp_circ_ap(field_in,*stk(l2),x_shift,y_shift);
	outptr = field2array( field_out );
	
	CreateVarFromPtr(5, "d", &m1, &n1, &outptr);   
	LhsVar(1) = 5;
	return 0;
}

int circ_screen_int(char *fname) 
{
	int  m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4,
		minlhs=1, maxlhs=1, minrhs=2, maxrhs=4, number, number2;
	FIELD field_in, field_out; double *outptr=NULL;
	double x_shift=0, y_shift=0;
  
	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1); 
	number=(int)(*stk(l1)); number2=(int)(*stk(l1+10));
	if ( (m1!=1)||(n1 != (2*number*number2+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}
	field_in=array2field( stk(l1) );

	GetRhsVar(2, "d", &m2, &n2, &l2); 
	if ( *stk(l2) < 0 )  {
		sciprint("Error: arguments R must be positive\r\n");
		return 0;
	}

	if (Rhs>2)
	{
		GetRhsVar(3, "d", &m3, &n3, &l3); 
		x_shift=*stk(l3);
	}
	
	if (Rhs>3)
	{
		GetRhsVar(4, "d", &m4, &n4, &l4);
		y_shift=*stk(l4);
	}

	field_out=lp_circ_screen(field_in,*stk(l2),x_shift,y_shift);
	outptr = field2array( field_out );
	
	CreateVarFromPtr(5, "d", &m1, &n1, &outptr);   
	LhsVar(1) = 5;
	return 0;
}

int absorber_int(char *fname) 
{
	int   m1, n1, l1, m2, n2, l2, 
		minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
	FIELD field_in, field_out;
	double *outptr=NULL;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	if ( (m1!=1)||(n1 != (2*((int)(*stk(l1)))*((int)(*stk(l1+10)))+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}
	field_in=array2field( stk(l1) );

	GetRhsVar(2, "d", &m2, &n2, &l2); 
	if ( *stk(l2) < 0 )  {
		sciprint("Error: arguments absorber coef must be positive\r\n");
		return 0;
	}

	field_out =lp_absorber( field_in,*stk(l2) );
	outptr = field2array( field_out );

	CreateVarFromPtr(3, "d", &m1, &n1, &outptr);   
	
	LhsVar(1) = 3;
	return 0;
}

int cros_out_int(char *fname) 
{
	int   m1, n1, l1, m2, n2, l2, m3, n3, l3,
		minlhs=1, maxlhs=1, minrhs=2, maxrhs=3;
	FIELD field_in, field_out;
	double time=-1;
	
	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	if ( (m1!=1)||(n1 != (2*((int)(*stk(l1)))*((int)(*stk(l1+10)))+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}
	field_in=array2field( stk(l1) );

	GetRhsVar(2, "c", &m2, &n2, &l2); 

	if(Rhs>2) {
		GetRhsVar(3, "d", &m3, &n3, &l3); 
		time = *stk(l3);
	}

	field_out=lp_cros_out( field_in, cstk(l2), time );

	LhsVar(1) = 1;
	return 0;
}

int file_pgm_int(char *fname) 
{
	int   m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4, m5, n5, l5, 
		m6, n6, l6, minlhs=1, maxlhs=1, minrhs=2, maxrhs=6;
	double gamma=2.0; int imax=128; int max_val=255;
	FIELD field_in, field_out;
	double time=-1;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	if ( (m1!=1)||(n1 != (2*((int)(*stk(l1)))*((int)(*stk(l1+10)))+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}
	field_in=array2field( stk(l1) );

	GetRhsVar(2, "c", &m2, &n2, &l2);  

	if(Rhs>2) {
		GetRhsVar(3, "d", &m3, &n3, &l3); 
		time = *stk(l3);
	}

	if (Rhs>3) {
		GetRhsVar(4, "d", &m4, &n4, &l4); 
		gamma=*stk(l4);
	}
	
	if (Rhs>4) {
		GetRhsVar(5, "i", &m5, &n5, &l5); 
		imax=*istk(l5); 
		  // equals to grid sampling if you pass the zero to imax
	}
	
	if (Rhs>5) {
		GetRhsVar(6, "i", &m6, &n6, &l6); 
		max_val = *istk(l6);
	}

	field_out=lp_file_pgm( field_in, cstk(l2), gamma, imax, max_val, time );

	LhsVar(1) = 1;
	return 0;
}

int file_ps_int(char *fname) 
{
	int   m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4,
		m5, n5, l5, minlhs=1, maxlhs=1, minrhs=2, maxrhs=5;
	double gamma=2.0; int imax=128; 
	FIELD field_in, field_out;
	double time=-1;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	if ( (m1!=1)||(n1 != (2*((int)(*stk(l1)))*((int)(*stk(l1+10)))+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}
	field_in=array2field( stk(l1) );

	GetRhsVar(2, "c", &m2, &n2, &l2);  

	if(Rhs>2) {
		GetRhsVar(3, "d", &m3, &n3, &l3); 
		time = *stk(l3);
	}

	if (Rhs>3) {
		GetRhsVar(4, "d", &m4, &n4, &l4); 
		gamma=*stk(l4);
	}
	
	if (Rhs>4) {
		GetRhsVar(5, "i", &m5, &n5, &l5);
		imax=*istk(l5); 
		  // equals to grid sampling if you pass the zero to imax
	}
	
	field_out=lp_file_ps( field_in, cstk(l2), gamma, imax, time );

	LhsVar(1) = 1;
	return 0;
}

int file_int_int(char *fname) 
{
	int  m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4, m5, n5, l5,
		gnuplot=0, imax=64, minlhs=1, maxlhs=1, minrhs=2, maxrhs=5;
	FIELD field_in, field_out;
	double time=-1;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);
	if ( (m1!=1)||(n1 != (2*((int)(*stk(l1)))*((int)(*stk(l1+10)))+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}
	field_in=array2field( stk(l1) );

	GetRhsVar(2, "c", &m2, &n2, &l2);  

	if (Rhs>2) {
		GetRhsVar(3, "d", &m3, &n3, &l3);
		time = *stk(l3);
	}
	
	if (Rhs>3) {
		GetRhsVar(4, "i", &m4, &n4, &l4); 
		imax=*istk(l4);  // equals to grid sampling if you pass the zero to imax
	}
	
	if (Rhs>4) {
		GetRhsVar(5, "i", &m5, &n5, &l5); 
		gnuplot = *istk(l5);  // equals to grid sampling if you pass the zero to imax
	}
	
	field_out=lp_file_int( field_in, cstk(l2), imax, time, gnuplot );
	
	LhsVar(1) = 1;
	return 0;
}

int file_pha_int(char *fname) 
{
	int  m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4, m5, n5, l5,
		gnuplot=0, imax=64, minlhs=1, maxlhs=1, minrhs=2, maxrhs=5;
	FIELD field_in, field_out;
	double time=-1;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);
	if ( (m1!=1)||(n1 != (2*((int)(*stk(l1)))*((int)(*stk(l1+10)))+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}
	field_in=array2field( stk(l1) );

	GetRhsVar(2, "c", &m2, &n2, &l2);  

	if (Rhs>2) {
		GetRhsVar(3, "d", &m3, &n3, &l3);
		time = *stk(l3);
	}
	
	if (Rhs>3) {
		GetRhsVar(4, "i", &m4, &n4, &l4); 
		imax=*istk(l4);  // equals to grid sampling if you pass the zero to imax
	}
	
	if (Rhs>4) {
		GetRhsVar(5, "i", &m5, &n5, &l5); 
		gnuplot = *istk(l5);  // equals to grid sampling if you pass the zero to imax
	}

	field_out=lp_file_pha( field_in, cstk(l2), imax, time, gnuplot );
	
	LhsVar(1) = 1;
	return 0;
}

int forvard_int(char *fname) 
{
	int   m1, n1, l1, m2, n2, l2,
		minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
	FIELD field_in, field_out;
	double *outptr=NULL;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);
	if ( (m1!=1)||(n1 != (2*((int)(*stk(l1)))*((int)(*stk(l1+10)))+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}
	field_in=array2field( stk(l1) );
    
	GetRhsVar(2, "d", &m2, &n2, &l2); 

	field_out=lp_forvard( field_in, *stk(l2) );
	outptr = field2array( field_out );

	CreateVarFromPtr(3, "d", &m1, &n1, &outptr);  
	
	LhsVar(1) = 3;
	return 0;
}

int fresnel_int(char *fname) 
{
	int   m1, n1, l1, m2, n2, l2,
		minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
	FIELD field_in, field_out;
	double *outptr=NULL;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);
	if ( (m1!=1)||(n1 != (2*((int)(*stk(l1)))*((int)(*stk(l1+10)))+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}
	field_in=array2field( stk(l1) );
    
	GetRhsVar(2, "d", &m2, &n2, &l2); 

	field_out=lp_fresnel( field_in, *stk(l2) );
	outptr = field2array( field_out );

	CreateVarFromPtr(3, "d", &m1, &n1, &outptr);  
	
	LhsVar(1) = 3;
	return 0;
}

int forward_int(char *fname) 
{
	int   m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4,
		m5, n5, minlhs=1, maxlhs=1, minrhs=2, maxrhs=4;
	int new_number=0; double new_size=0;
	FIELD field_in;	double *outptr=NULL;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);
	if ( (m1!=1)||(n1 != (2*((int)(*stk(l1)))*((int)(*stk(l1+10)))+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}
	field_in=array2field( stk(l1) );

	GetRhsVar(2, "d", &m2, &n2, &l2); 

	if (Rhs>2) {
		GetRhsVar(3, "d", &m3, &n3, &l3); 
		new_size=*stk(l3);
        
	}
	
	if (Rhs>3) {
		GetRhsVar(4, "i", &m4, &n4, &l4);
		new_number=*istk(l4);
		if ( new_number%2 != 0 ) {
			sciprint("Error: The number of grid must be even\r\n");
			return 0;
		}
		  // equals to grid sampling if you pass the zero to imax
	}

	outptr = lp_forward( field_in, *stk(l2), new_number, new_size );

	if(new_number==0)
		new_number=(int)outptr[0];
	m5=1; n5=2*new_number*new_number+30;

	CreateVarFromPtr(5, "d", &m5, &n5, &outptr);  

	LhsVar(1) = 5;
	return 0;
}

int gauss_int(char *fname) 
{
	int  m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4,
		m5, n5, l5, minlhs=1, maxlhs=1, minrhs=2, maxrhs=5;
	double x_shift=0, y_shift=0, AA=1;
	FIELD field_in, field_out;
	double *outptr=NULL;

	CheckRhs(minrhs,maxrhs);
	CheckLhs(minlhs,maxlhs);

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	if ( (m1!=1)||(n1 != (2*((int)(*stk(l1)))*((int)(*stk(l1+10)))+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}
	field_in=array2field( stk(l1) );

	GetRhsVar(2, "d", &m2, &n2, &l2); 
	if ( *stk(l2) < 0 )  {
		sciprint("Error: arguments R must be positive\r\n");
		return 0;
	}

	if (Rhs>2) {
		GetRhsVar(3, "d", &m3, &n3, &l3);
		x_shift= *stk(l3);
	}
	
	if (Rhs>3) {
		GetRhsVar(4, "d", &m4, &n4, &l4);
		y_shift= *stk(l4);
	}
	
	if (Rhs>4) {
		GetRhsVar(5, "d", &m5, &n5, &l5); 
		AA=*stk(l5);
		if ( AA < 0 )  {
			sciprint("Error: maximum intensity transmission \
					 arguments must be positive\r\n");
			return 0;
		}
	}

	field_out=lp_gauss( field_in, *stk(l2), x_shift, y_shift, AA );
	outptr = field2array( field_out );
	
	CreateVarFromPtr(6, "d", &m1, &n1, &outptr);  

	LhsVar(1) = 6;
	return 0;
}

int gauss_screen_int(char *fname) 
{
	int  m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4,
		m5, n5, l5, minlhs=1, maxlhs=1, minrhs=2, maxrhs=5;
	double x_shift=0, y_shift=0, AA=1;
	FIELD field_in, field_out;
	double *outptr=NULL;

	CheckRhs(minrhs,maxrhs);
	CheckLhs(minlhs,maxlhs);

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	if ( (m1!=1)||(n1 != (2*((int)(*stk(l1)))*((int)(*stk(l1+10)))+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}
	field_in=array2field( stk(l1) );

	GetRhsVar(2, "d", &m2, &n2, &l2); 
	if ( *stk(l2) < 0 )  {
		sciprint("Error: arguments R must be positive\r\n");
		return 0;
	}

	if (Rhs>2) {
		GetRhsVar(3, "d", &m3, &n3, &l3);
		x_shift= *stk(l3);
	}
	
	if (Rhs>3) {
		GetRhsVar(4, "d", &m4, &n4, &l4);
		y_shift= *stk(l4);
	}
	
	if (Rhs>4) {
		GetRhsVar(5, "d", &m5, &n5, &l5); 
		AA=*stk(l5);
		if ( AA < 0 )  {
			sciprint("Error: maximum intensity transmission \
					 arguments must be positive\r\n");
			return 0;
		}
	}

	field_out=lp_gauss_screen( field_in, *stk(l2), x_shift, y_shift, AA );
	outptr = field2array( field_out );
	
	CreateVarFromPtr(6, "d", &m1, &n1, &outptr);  

	LhsVar(1) = 6;
	return 0;
}

int rect_ap_int(char *fname) 
{
	int   m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4, m5, n5, l5,
		m6, n6, l6, minlhs=1, maxlhs=1, minrhs=2, maxrhs=6;
	double sy, x_shift=0, y_shift=0, angle=0;
	FIELD field_in, field_out;
	double *outptr;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	if ( (m1!=1)||(n1 != (2*((int)(*stk(l1)))*((int)(*stk(l1+10)))+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}
	field_in=array2field( stk(l1) );
	
	GetRhsVar(2, "d", &m2, &n2, &l2);
	if ( *stk(l2) < 0 )  {
		sciprint("Error: sx arguments must be positive\r\n");
		return 0;
	}

	sy=*stk(l2);

	if (Rhs>2) {
		GetRhsVar(3, "d", &m3, &n3, &l3);
		if ( *stk(l3) <= 0 )  {
			sciprint("Error: sy arguments must be positive\r\n");
			return 0;
		}
		sy=*stk(l3);        
	} 
	
	if (Rhs>3) {
		GetRhsVar(4, "d", &m4, &n4, &l4);
		x_shift=*stk(l4);
	}
	
	if (Rhs>4) {
		GetRhsVar(5, "d", &m5, &n5, &l5);
		y_shift=*stk(l5);
	}
	
	if (Rhs>5) {
		GetRhsVar(6, "d", &m6, &n6, &l6);
		angle=*stk(l6);
	}
    
	field_out=lp_rect_ap( field_in, *stk(l2), sy, x_shift, y_shift, angle );

	outptr = field2array( field_out );

	CreateVarFromPtr(7, "d", &m1, &n1, &outptr);   
	
	LhsVar(1) = 7;
	return 0;
}

int rect_screen_int(char *fname) 
{
	int   m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4, m5, n5, l5,
		m6, n6, l6, minlhs=1, maxlhs=1, minrhs=2, maxrhs=6;
	double sy, x_shift=0, y_shift=0, angle=0;
	FIELD field_in, field_out;
	double *outptr;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	if ( (m1!=1)||(n1 != (2*((int)(*stk(l1)))*((int)(*stk(l1+10)))+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}
	field_in=array2field( stk(l1) );
	
	GetRhsVar(2, "d", &m2, &n2, &l2);
	if ( *stk(l2) < 0 )  {
		sciprint("Error: sx arguments must be positive\r\n");
		return 0;
	}

	sy=*stk(l2);

	if (Rhs>2) {
		GetRhsVar(3, "d", &m3, &n3, &l3);
		if ( *stk(l3) <= 0 )  {
			sciprint("Error: sy arguments must be positive\r\n");
			return 0;
		}
		sy=*stk(l3);        
	} 
	
	if (Rhs>3) {
		GetRhsVar(4, "d", &m4, &n4, &l4);
		x_shift=*stk(l4);
	}
	
	if (Rhs>4) {
		GetRhsVar(5, "d", &m5, &n5, &l5);
		y_shift=*stk(l5);
	}
	
	if (Rhs>5) {
		GetRhsVar(6, "d", &m6, &n6, &l6);
		angle=*stk(l6);
	}
    
	field_out=lp_rect_screen( field_in, *stk(l2), sy, x_shift, y_shift, angle );

	outptr = field2array( field_out );

	CreateVarFromPtr(7, "d", &m1, &n1, &outptr);   
	
	LhsVar(1) = 7;
	return 0;
}

int random_int(char *fname) 
{
	int  m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4,
		minlhs=1, maxlhs=1, minrhs=3, maxrhs=4;
	FIELD field_in, field_out;
	double *outptr;
	int seed=0;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	if ( (m1!=1)||(n1 != (2*((int)(*stk(l1)))*((int)(*stk(l1+10)))+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}
	field_in=array2field( stk(l1) );
	
	GetRhsVar(2, "c", &m2, &n2, &l2); 
	GetRhsVar(3, "d", &m3, &n3, &l3);  
	
	if (Rhs>3) {
		GetRhsVar(4, "i", &m4, &n4, &l4); 
		seed=*istk(l4);
	}

	field_out=lp_random( field_in, cstk(l2), *stk(l3), seed );
	outptr = field2array( field_out );

	CreateVarFromPtr(5, "d", &m1, &n1, &outptr);  
	
	LhsVar(1) = 5;
	return 0;
}

int lp_zernike_int(char *fname) 
{
	int  m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4, n, m,
		m5, n5, l5, minlhs=1, maxlhs=1, minrhs=5, maxrhs=5;
	FIELD field_in, field_out;
	double R, AA, *outptr;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	if ( (m1!=1)||(n1 != (2*((int)(*stk(l1)))*((int)(*stk(l1+10)))+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}
	field_in=array2field( stk(l1) );
    
	GetRhsVar(2, "i", &m2, &n2, &l2); 
	n = *istk(l2);
    if(n<0){
        sciprint("Error: n must be positive\n");
        return 0;
    }

	GetRhsVar(3, "i", &m3, &n3, &l3);
	m = *istk(l3); 
    if(abs(m)>n){
        sciprint("Error: |m| must be less or equal than n\n");
        return 0;
    }
    if((m%2)!=(n%2)){
        sciprint("Error: n and m must all odd or even\n");
        return 0;
    }

	GetRhsVar(4, "d", &m4, &n4, &l4); R = *stk(l4); 
	GetRhsVar(5, "d", &m5, &n5, &l5); AA=*stk(l5);

	field_out=lp_zernike( field_in, n, m, R, AA );
	outptr = field2array( field_out );

	CreateVarFromPtr(6, "d", &m1, &n1, &outptr);   
	LhsVar(1) = 6;
	return 0;
}

int l_amplif_int(char *fname) 
{
	int  m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4,
		minlhs=1, maxlhs=1, minrhs=4, maxrhs=4;
	FIELD field_in, field_out;
	double *outptr;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);
	if ( (m1!=1)||(n1 != (2*((int)(*stk(l1)))*((int)(*stk(l1+10)))+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}
	field_in=array2field( stk(l1) );

	GetRhsVar(2, "d", &m2, &n2, &l2); 
	GetRhsVar(3, "d", &m3, &n3, &l3);  
	if ( *stk(l3) < 0 )  {
		sciprint("Error: arguments base  must be positive\r\n");
		return 0;
	}

	GetRhsVar(4, "d", &m4, &n4, &l4);  
	if ( *stk(l4) < 0 )  {
		sciprint("Error: arguments i_sat  must be positive\r\n");
		return 0;
	}

	field_out=lp_l_amplif( field_in, *stk(l2), *stk(l3), *stk(l4) );
	outptr = field2array( field_out );

	CreateVarFromPtr(5, "d", &m1, &n1, &outptr);   
	
	LhsVar(1) = 5;
	return 0;
}

int lens_int(char *fname) 
{
	int   m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4,
		minlhs=1, maxlhs=1, minrhs=2, maxrhs=4;
	double f, x_shift=0, y_shift=0;
	FIELD field_in, field_out;
	double *outptr;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	if ( (m1!=1)||(n1 != (2*((int)(*stk(l1)))*((int)(*stk(l1+10)))+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}
	field_in=array2field( stk(l1) );

	GetRhsVar(2, "d", &m2, &n2, &l2); 
	f = *stk(l2);

	if (Rhs>2) {
		GetRhsVar(3, "d", &m3, &n3, &l3);
		x_shift=*stk(l3);
	}
	
	if (Rhs>3) {
		GetRhsVar(4, "d", &m4, &n4, &l4); 
		y_shift=*stk(l4);
	}

	field_out=lp_lens( field_in, f, x_shift, y_shift );
	outptr = field2array( field_out );

	CreateVarFromPtr(5, "d", &m1, &n1, &outptr);  
	
	LhsVar(1) = 5;
	return 0;
}

int lens_forvard_int(char *fname) 
{
	int   m1, n1, l1, m2, n2, l2, m3, n3, l3,
		minlhs=1, maxlhs=1, minrhs=3, maxrhs=3;
	FIELD field_in, field_out;
	double f, z, *outptr;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	if ( (m1!=1)||(n1 != (2*((int)(*stk(l1)))*((int)(*stk(l1+10)))+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}
	field_in=array2field( stk(l1) );

	GetRhsVar(2, "d", &m2, &n2, &l2); 
	f = *stk(l2);

	GetRhsVar(3, "d", &m3, &n3, &l3);  
	z = *stk(l3);

	field_out=lp_lens_forvard( field_in, f, z );
	outptr = field2array( field_out );
	
	CreateVarFromPtr(4, "d", &m1, &n1, &outptr);

	LhsVar(1) = 4;
	return 0;
}

int lens_fresn_int(char *fname) 
{
	int   m1, n1, l1, m2, n2, l2, m3, n3, l3,
		minlhs=1, maxlhs=1, minrhs=3, maxrhs=3;
	FIELD field_in, field_out;
	double f, z, *outptr;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	if ( (m1!=1)||(n1 != (2*((int)(*stk(l1)))*((int)(*stk(l1+10)))+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}
	field_in=array2field( stk(l1) );

	GetRhsVar(2, "d", &m2, &n2, &l2); 
	f = *stk(l2);

	GetRhsVar(3, "d", &m3, &n3, &l3);  
	z = *stk(l3);

	field_out=lp_lens_fresn( field_in, f, z );
	outptr = field2array( field_out );
	
	CreateVarFromPtr(4, "d", &m1, &n1, &outptr);

	LhsVar(1) = 4;
	return 0;
}

int convert_int(char *fname) 
{
	int   m1, n1, l1,  minlhs=1, 
		maxlhs=1, minrhs=1, maxrhs=1;
	FIELD field_in, field_out;
	double *outptr=NULL;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	if ( (m1!=1)||(n1 != (2*((int)(*stk(l1)))*((int)(*stk(l1+10)))+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}
	field_in=array2field( stk(l1) );

	field_out=lp_convert( field_in );
	outptr = field2array( field_out );

	CreateVarFromPtr(2, "d", &m1, &n1, &outptr);  
	LhsVar(1) = 2;
	return 0;
}

int normal_int(char *fname) 
{
	int   m1, n1, l1,  minlhs=1, 
		maxlhs=1, minrhs=1, maxrhs=1;
	FIELD field_in, field_out;
	double *outptr=NULL;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	if ( (m1!=1)||(n1 != (2*((int)(*stk(l1)))*((int)(*stk(l1+10)))+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}
	field_in=array2field( stk(l1) );

	field_out=lp_normal( field_in );
	outptr = field2array( field_out );

	CreateVarFromPtr(2, "d", &m1, &n1, &outptr);  
	LhsVar(1) = 2;
	return 0;
}

int tilt_int(char *fname) 
{
	int  m1, n1, l1, m2, n2, l2, m3, n3, l3,
	  minlhs=1, maxlhs=1, minrhs=3, maxrhs=3;
	FIELD field_in, field_out;
	double tx,ty,*outptr;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	if ( (m1!=1)||(n1 != (2*((int)(*stk(l1)))*((int)(*stk(l1+10)))+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}
	field_in=array2field( stk(l1) );

	GetRhsVar(2, "d", &m2, &n2, &l2);
	tx= *stk(l2);

	GetRhsVar(3, "d", &m3, &n3, &l3);  
	ty= *stk(l3);

	field_out=lp_tilt( field_in, tx, ty );
	outptr = field2array( field_out );

	CreateVarFromPtr(4, "d", &m1, &n1, &outptr);   
	LhsVar(1) = 4;
	return 0;
}

int tor_lens_int(char *fname) 
{
	int   m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4,
		m5, n5, l5, minlhs=1, maxlhs=1, minrhs=3, maxrhs=5;
	double R, f, x_shift=0, y_shift=0;
	FIELD field_in, field_out;
	double *outptr;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;
	
	GetRhsVar(1, "d", &m1, &n1, &l1);  
	if ( (m1!=1)||(n1 != (2*((int)(*stk(l1)))*((int)(*stk(l1+10)))+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}
	field_in=array2field( stk(l1) );

	GetRhsVar(2, "d", &m2, &n2, &l2);
	R= *stk(l2);
	if ( R < 0 )  {
		sciprint("Error: arguments R  must be positive\r\n");
		return 0;
	}

	GetRhsVar(3, "d", &m3, &n3, &l3);
	f= *stk(l3);

	if (Rhs>3) 
	{
		GetRhsVar(4, "d", &m4, &n4, &l4);
		x_shift= *stk(l4);
	}
	
	if (Rhs>4) 
	{
		GetRhsVar(5, "d", &m5, &n5, &l5);
		y_shift= *stk(l5);
	}
	
	field_out=lp_tor_lens( field_in, R, f, x_shift, y_shift );
	outptr = field2array( field_out );
	
	CreateVarFromPtr(6, "d", &m1, &n1, &outptr);  
	LhsVar(1) = 6;
	return 0;
}

int pip_fft_int(char *fname) 
{
	int m1, n1, l1, m2, n2, l2, ind,
		minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
	FIELD field_in, field_out;
	double *outptr;
	
	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;
	
	GetRhsVar(1, "d", &m1, &n1, &l1);
	if ( (m1!=1)||(n1 != (2*((int)(*stk(l1)))*((int)(*stk(l1+10)))+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}
	field_in=array2field( stk(l1) );

	GetRhsVar(2, "i", &m2, &n2, &l2); 
	ind = *istk(l2); 

	field_out=lp_pip_fft( field_in, ind );
	outptr = field2array( field_out );
	
	CreateVarFromPtr(3, "d", &m1, &n1, &outptr);  
	LhsVar(1) = 3;
	return 0;
}

int interp1_int(char *fname) 
{
	int m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4,m5, n5,
		l5, m6, n6, l6, m7, n7, l7, m8, n8,  new_number=0,
		minlhs=1, maxlhs=1, minrhs=2, maxrhs=7;
	double size_new, x_shift=0, y_shift=0, angle=0, magnif=1;
	double *outptr=NULL;
	FIELD field_in;
    
	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	if ( (m1!=1)||(n1 != (2*((int)(*stk(l1)))*((int)(*stk(l1+10)))+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}
	field_in=array2field( stk(l1) );

	GetRhsVar(2, "d", &m2, &n2, &l2); 
	size_new= *stk(l2);

	if (Rhs>2) {
		GetRhsVar(3, "i", &m3, &n3, &l3); 
		new_number=*istk(l3);
		if ( new_number < 0 )  {
			sciprint("Error: arguments new_number  must be positive\r\n");
			return 0;
		}
	}
	
	if (Rhs>3) {
		GetRhsVar(4, "d", &m4, &n4, &l4); 
		x_shift=*stk(l4);
	}
	
	if (Rhs>4) {
		GetRhsVar(5, "d", &m5, &n5, &l5);
		y_shift=*stk(l5);
	}
	
	if (Rhs>5) {
		GetRhsVar(6, "d", &m6, &n6, &l6); 
		angle=*stk(l6);
	}
	
	if (Rhs>6) {
		GetRhsVar(7, "d", &m7, &n7, &l7);  
		magnif=*stk(l7);
	}

	outptr = lp_interp1( field_in, size_new, 
		new_number, x_shift, y_shift, angle, magnif);

	if (new_number==0)
		new_number=(int) outptr[0];
	m8=1; n8=2*new_number*new_number+30;
	CreateVarFromPtr(8, "d", &m8, &n8, &outptr);   

	LhsVar(1) = 8;
	return 0;
}

int interpol_int(char *fname) 
{
	int m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4,m5, n5,
		l5, m6, n6, l6, m7, n7, l7, m8, n8,  new_number=0,
		minlhs=1, maxlhs=1, minrhs=2, maxrhs=7;
	double size_new, x_shift=0, y_shift=0, angle=0, magnif=1;
	double *outptr=NULL;
	FIELD field_in;
    
	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	if ( (m1!=1)||(n1 != (2*((int)(*stk(l1)))*((int)(*stk(l1+10)))+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}
	field_in=array2field( stk(l1) );

	GetRhsVar(2, "d", &m2, &n2, &l2); 
	size_new= *stk(l2);

	if (Rhs>2) {
		GetRhsVar(3, "i", &m3, &n3, &l3); 
		new_number=*istk(l3);
		if ( new_number < 0 )  {
			sciprint("Error: arguments new_number  must be positive\r\n");
			return 0;
		}
	}
	
	if (Rhs>3) {
		GetRhsVar(4, "d", &m4, &n4, &l4); 
		x_shift=*stk(l4);
	}
	
	if (Rhs>4) {
		GetRhsVar(5, "d", &m5, &n5, &l5);
		y_shift=*stk(l5);
	}
	
	if (Rhs>5) {
		GetRhsVar(6, "d", &m6, &n6, &l6); 
		angle=*stk(l6);
	}
	
	if (Rhs>6) {
		GetRhsVar(7, "d", &m7, &n7, &l7);  
		magnif=*stk(l7);
	}

	outptr = lp_interpol( field_in, size_new, 
		new_number, x_shift, y_shift, angle, magnif);

	if (new_number==0) 
		new_number=(int) outptr[0];
	m8=1; n8=2*new_number*new_number+30;
	CreateVarFromPtr(8, "d", &m8, &n8, &outptr);   

	LhsVar(1) = 8;
	return 0;
}

int b_mix_int(char *fname) 
{
	int   m1, n1, l1, m2, n2, l2,
		minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
	FIELD field_in_1, field_in_2, field_out;
	double *outptr=NULL;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	if ( (m1!=1)||(n1 != (2*((int)(*stk(l1)))*((int)(*stk(l1+10)))+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}
	field_in_1=array2field( stk(l1) );

	GetRhsVar(2, "d", &m2, &n2, &l2); 
	if ( (m2!=1)||(n2 != (2*((int)(*stk(l2)))*((int)(*stk(l2+10)))+30) ) )  {
		sciprint("Error: second arguments must be diffractive wavefront \r\n");
		return 0;
	}
	field_in_2=array2field( stk(l2) );
	
	field_out=lp_b_mix(field_in_1,field_in_2);
	outptr = field2array( field_out );
	
	CreateVarFromPtr(3, "d", &m1, &n1, &outptr);  
	LhsVar(1) = 3;
	return 0;
}

int b_split_int(char *fname) 
{
	int m_split, n_split, l_split, m1, n1, l1,
        minlhs=1, maxlhs=2, minrhs=1, maxrhs=2;
	double *outptr, *outptr1, *outptr2;
	FIELD field_in;  double split=0.5;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1); 
	if ( (m1!=1)||(n1 != (2*((int)(*stk(l1)))*((int)(*stk(l1+10)))+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}
	field_in=array2field( stk(l1) );

	if (Rhs==2) {
		GetRhsVar(2, "d", &m_split, &n_split, &l_split);
		if ( *stk(l_split) < 0 )  {
			sciprint("Error: arguments must be positive\r\n");
			return 0;
		}
		split=*stk(l_split);
		if( split== 0.) split=0.000001;
	}

	outptr =lp_b_split(field_in,split);

	outptr1=outptr; 
	CreateVarFromPtr(3, "d", &m1, &n1, &outptr1);   

	LhsVar(1) = 3;
	if (Lhs==2) {
		outptr2=outptr+n1;
		CreateVarFromPtr(4, "d", &m1, &n1, &outptr2); 
		LhsVar(2) = 4;
	}

	return 0;
}

int fil_ter_int(char *fname) 
{
	int   m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4,
		m5, n5, l5, minlhs=1, maxlhs=1, minrhs=4, maxrhs=5;
	FIELD field_in, field_out;
	double *outptr=NULL;
	int norm=0;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	if ( (m1!=1)||(n1 != (2*((int)(*stk(l1)))*((int)(*stk(l1+10)))+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}
	field_in=array2field( stk(l1) );

	GetRhsVar(2, "c", &m2, &n2, &l2); 
	GetRhsVar(3, "c", &m3, &n3, &l3);  
	GetRhsVar(4, "c", &m4, &n4, &l4);  

	if (Rhs>4) {
		GetRhsVar(5, "d", &m5, &n5, &l5); 
		norm=1;
	}

	field_out=lp_fil_ter( field_in, cstk(l2), cstk(l3), cstk(l4), norm );
	outptr = field2array( field_out );
	
	CreateVarFromPtr(6, "d", &m1, &n1, &outptr);   
	LhsVar(1) = 6;
	return 0;
}

int steps_int(char *fname) 
{
	int   m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4,
		m5, n5, l5, m6, n6, l6, m7, n7, l7, nstep=1, n_st=1, 
		minlhs=1, maxlhs=1, minrhs=2, maxrhs=7;
	char *rfname="void", *afname="void", *ofname=NULL;
	FIELD field_in, field_out;
	double z, *outptr=NULL;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	if ( (m1!=1)||(n1 != (2*((int)(*stk(l1)))*((int)(*stk(l1+10)))+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}
	field_in=array2field( stk(l1) );
	
	GetRhsVar(2, "d", &m2, &n2, &l2); 
	z= *stk(l2);
    
	if (Rhs>2) {
		GetRhsVar(3, "i", &m3, &n3, &l3); 
		nstep=*istk(l3);
		if ( nstep < 0 )  {
			sciprint("Error: arguments  nstep  must be positive\r\n");
			return 0;
		}
	}
	if (Rhs>3) {
		GetRhsVar(4, "c", &m4, &n4, &l4);
		rfname=cstk(l4);
	}
	  
	if (Rhs>4) {
		GetRhsVar(5, "c", &m5, &n5, &l5);
		afname=cstk(l5);
	}
	if (Rhs>5) {
		GetRhsVar(6, "c", &m6, &n6, &l6);
		ofname=cstk(l6);
	}
	if (Rhs>6) {
		GetRhsVar(7, "i", &m7, &n7, &l7);
		n_st=*istk(l7);
		if ( n_st < 0 )  {
			sciprint("Error: arguments n_st  must be positive\r\n");
			return 0;
		}
	}

	field_out=lp_steps( field_in, z, nstep, rfname, afname, ofname, n_st );
	outptr = field2array( field_out );
	
	CreateVarFromPtr(8, "d", &m1, &n1, &outptr);  
	LhsVar(1) = 8;
	return 0;
}

int strehl_int(char *fname) 
{
	int m1, n1, l1, m2, n2,  
		minlhs=1, maxlhs=2, minrhs=1, maxrhs=1;
	FIELD field_in;
	double *outptr;

	CheckRhs(minrhs,maxrhs);
	CheckLhs(minlhs,maxlhs);

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	if ( (m1!=1)||(n1 != (2*((int)(*stk(l1)))*((int)(*stk(l1+10)))+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}
	field_in=array2field( stk(l1) );

	outptr = lp_strehl( field_in );
	
	m2=1; n2=7;
	CreateVarFromPtr(2, "d", &m2, &n2, &outptr);   

	LhsVar(1) = 2;
	if(Lhs==2) LhsVar(2) = 1;
	return 0;
}

int field_int_int(char *fname) 
{
	int m1, n1, l1, m2, n2, l2, m3, n3, amp=0,
		minlhs=1, maxlhs=1, minrhs=1, maxrhs=2;
	double *outptr;

	CheckRhs(minrhs,maxrhs);
	CheckLhs(minlhs,maxlhs);

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	if ( (m1!=1)||(n1 != (2*((int)(*stk(l1)))*((int)(*stk(l1+10)))+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}

	if (Rhs>1) {
		GetRhsVar(2, "i", &m2, &n2, &l2); 
		amp = *istk(l2);
	}

	outptr=lp_field_int( stk(l1), amp );

	m3=(int)(*stk(l1)); n3=(int)(*stk(l1+10));
	CreateVarFromPtr(3, "d", &m3, &n3, &outptr);   

	LhsVar(1) = 3;
	return 0;
}

int field_pha_int(char *fname) 
{
	int m1, n1, l1, m2, n2,
		minlhs=1, maxlhs=1, minrhs=1, maxrhs=1;
	double *outptr;

	CheckRhs(minrhs,maxrhs);
	CheckLhs(minlhs,maxlhs);

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	if ( (m1!=1)||(n1 != (2*((int)(*stk(l1)))*((int)(*stk(l1+10)))+30) ) )  {
		sciprint("Error: first arguments must be diffractive wavefront \r\n");
		return 0;
	}

	outptr=lp_field_pha( stk(l1) );

	m2=(int)(*stk(l1)); n2=(int)(*stk(l1+10));
	CreateVarFromPtr(2, "d", &m2, &n2, &outptr);   

	LhsVar(1) = 2;
	return 0;
}

int create_field_int(char *fname) 
{
	int m1, n1, l1, m2, n2,l2, n3,m3,l3, n4, m4, l4, m5, n5, nelem,
		number,number2,real_imag=0, minlhs=1, maxlhs=1, minrhs=3, maxrhs=4;
	double *outptr=NULL;
	FIELD field;

	CheckRhs(minrhs,maxrhs);
	CheckLhs(minlhs,maxlhs);

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	if ( (m1!=1)||(n1 != 30 ) )  {
		sciprint("Error: first arguments must be diffractive wavefront header\r\n");
		return 0;
	}
	number=(int)(*stk(l1));number2=(int)(*stk(l1+10));
	nelem=number*number2;

	GetRhsVar(2, "d", &m2, &n2, &l2);  
	if ( (m2!=number)||(n2 != number2 ) )  {
		sciprint("Error: first arguments must be field amplitude matrix\r\n");
		return 0;
	}

	GetRhsVar(3, "d", &m3, &n3, &l3);  
	if ( (m3!=number)||(n3 != number2 ) )  {
		sciprint("Error: first arguments must be field phase matrix\r\n");
		return 0;
	}

	if (Rhs>3) {
		GetRhsVar(4, "i", &m4, &n4, &l4);
		real_imag = *istk(l4);
	}

	field=create_field( number,number2);
	header_array_field(&field,stk(l1));

	field.real_imag=real_imag;

    field.real=stk(l2);
	field.imaginary=stk(l3);
	real_imaginary_field(&field);

	outptr = field2array( field );

	m5=1; n5 = 2*nelem+30;
	CreateVarFromPtr( 5, "d", &m5, &n5, &outptr );   

	LhsVar(1) = 5;
	return 0;
}

int field_contents_int(char *fname) 
{
	int m1, n1, l1, m2, n2, m3, n3, number,	number2, 
		nelem, minlhs=1, maxlhs=3, minrhs=1, maxrhs=1;
	double *outptr1,*outptr2,*outptr3;
	FIELD field_in;

	CheckRhs(minrhs,maxrhs);
	CheckLhs(minlhs,maxlhs);

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	number=(int)(*stk(l1));number2=(int)(*stk(l1+10));
	if ( (m1!=1)||(n1 != 2*number*number2+30 ) )  {
		sciprint("Error: argument must be diffractive wavefront \r\n");
		return 0;
	}
	field_in=array2field( stk(l1) );
	nelem=number*number2;

	outptr1=(double*) calloc( 30, sizeof(double) );
	field_header_array(field_in,outptr1);

	outptr2=field_in.real;
	outptr3=field_in.imaginary;

	m2=1;n2=30;
	CreateVarFromPtr(2, "d", &m2, &n2, &outptr1);   

	LhsVar(1) = 2;
	
	if(Lhs>1) {
		m3=number; n3=number2;
		CreateVarFromPtr(3, "d", &m3, &n3, &outptr2);   
		LhsVar(2) = 3;
	}
	if(Lhs>2) {
		CreateVarFromPtr(4, "d", &m3, &n3, &outptr3);   
		LhsVar(3) = 4;
	}

	return 0;
}

int unfold_phase_int(char *fname) 
{
	int m1, n1, l1, m2, n2,l2, m3, n3, l3, m4, n4, l4, 
		type=0, minlhs=1, maxlhs=1, minrhs=4, maxrhs=5;
	double *outptr;

	CheckRhs(minrhs,maxrhs);
	CheckLhs(minlhs,maxlhs);

	GetRhsVar(1, "i", &m1, &n1, &l1);  

	GetRhsVar(2, "d", &m2, &n2, &l2);  
	if ( (m2 != *istk(l1) )||(n2 != *istk(l1) ) )  {
		sciprint("Error: second arguments must be phase matrix \r\n");
		return 0;
	}
	GetRhsVar(3, "d", &m3, &n3, &l3);
	GetRhsVar(4, "d", &m4, &n4, &l4);  

	if (Rhs>4) type= 1;

	if(type) outptr=unf4( *stk(l3),*istk(l1),stk(l2) );
	else outptr=unf3( *stk(l3),*stk(l4),*istk(l1),stk(l2) );

	CreateVarFromPtr(6, "d", &m2, &n2, &outptr);   

	LhsVar(1) = 6;
	return 0;
}

int c_scilab_int(char *fname) 
{
	int  m1, n1, l1, m2, n2, l2, 
		minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "c", &m1, &n1, &l1);  
	GetRhsVar(2, "c", &m2, &n2, &l2);  

	c_scilab( cstk(l1), cstk(l2) );

	LhsVar(1) = 2;
	return 0;
}
