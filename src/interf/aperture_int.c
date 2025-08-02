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

double *create_aperture( int aperture_type, double ap_data1,
		 double ap_data2, double ap_data3, double ap_data4 )
{
	APERTURE AP=create_APERTURE(aperture_type,
		    ap_data1, ap_data2, ap_data3, ap_data4);
	return (APERTURE2array(AP));
}

int circular_aperture_int(char *fname) 
{
  int   m1, n1, l1, m2, n2, 
	  minlhs=1, maxlhs=1, minrhs=1, maxrhs=1;
  double *outptr=NULL;

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  GetRhsVar(1, "d", &m1, &n1, &l1);  
  if ( *stk(l1) < 0 )  {
	  sciprint("Error: arguments diameter must be positive\r\n");
	  return 0;
  }

  outptr = create_aperture( 0, *stk(l1), 0, 0, 0 );

  m2=1; n2=17;
  CreateVarFromPtr(2, "d", &m2, &n2, &outptr);   

  LhsVar(1) = 2;
  return 0;
}

int annular_aperture_int(char *fname) 
{
  int   m1, n1, l1, m2, n2, l2, m3, n3, 
	  minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
  double *outptr=NULL;

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  GetRhsVar(1, "d", &m1, &n1, &l1);  
  if ( *stk(l1) < 0 )  {
	  sciprint("Error: arguments out_diameter must be positive\r\n");
	  return 0;
  }

  GetRhsVar(2, "d", &m2, &n2, &l2); 
  if ( *stk(l2) < 0 )  {
	  sciprint("Error: arguments in_diameter must be positive\r\n");
	  return 0;
  }

  if ( *stk(l1) < *stk(l2) )  {
	  sciprint("Error: out_diameter must greater than in_diameter \r\n");
	  return 0;
  }

  outptr = create_aperture( 1, *stk(l1), *stk(l2), 0, 0 );

  m3=1; n3=17;
  CreateVarFromPtr(3, "d", &m3, &n3, &outptr);   

  LhsVar(1) = 3;
  return 0;
}

int rectangular_aperture_int(char *fname) 
{
  int   m1, n1, l1, m2, n2, l2, m3, n3, 
	  minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
  double *outptr=NULL;

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  GetRhsVar(1, "d", &m1, &n1, &l1);  
  if ( *stk(l1) < 0 )  {
	  sciprint("Error: arguments x_size must be positive\r\n");
	  return 0;
  }

  GetRhsVar(2, "d", &m2, &n2, &l2); 
  if ( *stk(l2) < 0 )  {
	  sciprint("Error: arguments y_size must be positive\r\n");
	  return 0;
  }

  outptr = create_aperture( 2, *stk(l1), *stk(l2), 0, 0 );

  m3=1; n3=17;
  CreateVarFromPtr(3, "d", &m3, &n3, &outptr);   

  LhsVar(1) = 3;
  return 0;
}

int hexagonal_aperture_int(char *fname) 
{
  int   m1, n1, l1, m2, n2, 
	  minlhs=1, maxlhs=1, minrhs=1, maxrhs=1;
  double *outptr=NULL;

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  GetRhsVar(1, "d", &m1, &n1, &l1);  
  if ( *stk(l1) < 0 )  {
	  sciprint("Error: arguments in_edge_length must be positive\r\n");
	  return 0;
  }

  outptr = create_aperture( 3, *stk(l1), 0, 0, 0 );

  m2=1; n2=17;
  CreateVarFromPtr(2, "d", &m2, &n2, &outptr);   

  LhsVar(1) = 2;
  return 0;
}

int spidered_annular_aperture_int(char *fname) 
{
  int   m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4,
	  m5, n5, minlhs=1, maxlhs=1, minrhs=4, maxrhs=4;
  double *outptr=NULL;

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  GetRhsVar(1, "d", &m1, &n1, &l1);  

  if ( *stk(l1) < 0 )  {
	  sciprint("Error: arguments out_diameter must be positive\r\n");
	  return 0;
  }

  GetRhsVar(2, "d", &m2, &n2, &l2); 
  if ( *stk(l2) < 0 )  {
	  sciprint("Error: arguments in_diameter must be positive\r\n");
	  return 0;
  }

  if ( *stk(l1) < *stk(l2) )  {
	  sciprint("Error: out_diameter must greater than in_diameter \r\n");
	  return 0;
  }

  GetRhsVar(3, "i", &m3, &n3, &l3); 
  if ( *istk(l3) < 0 )  {
	  sciprint("Error: arguments nspiders must be position interger\r\n");
	  return 0;
  }

  GetRhsVar(4, "d", &m4, &n4, &l4); 
  if ( *stk(l4) < 0 )  {
	  sciprint("Error: arguments spider_width must be positive\r\n");
	  return 0;
  }

  outptr = create_aperture( 4, *stk(l1), *stk(l2), *istk(l3), *stk(l4) );

  m5=1; n5=17;
  CreateVarFromPtr(5, "d", &m5, &n5, &outptr);   

  LhsVar(1) = 5;
  return 0;
}

int tiled_hexagonal_aperture_int(char *fname) 
{
  int   m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4,
	  m5, n5, minlhs=1, maxlhs=1, minrhs=4, maxrhs=4;
  double *outptr=NULL;

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  GetRhsVar(1, "d", &m1, &n1, &l1);  

  if ( *stk(l1) < 0 )  {
	  sciprint("Error: arguments out_diameter must be positive\r\n");
	  return 0;
  }

  GetRhsVar(2, "d", &m2, &n2, &l2); 
  if ( *stk(l2) < 0 )  {
	  sciprint("Error: arguments in_diameter must be positive\r\n");
	  return 0;
  }

  if ( *stk(l1) < *stk(l2) )  {
	  sciprint("Error: out_diameter must greater than in_diameter \r\n");
	  return 0;
  }

  GetRhsVar(3, "d", &m3, &n3, &l3); 
  if ( *stk(l3) < 0 )  {
	  sciprint("Error: arguments in_edge_length must be position \r\n");
	  return 0;
  }

  GetRhsVar(4, "d", &m4, &n4, &l4); 
  if ( *stk(l4) < 0 )  {
	  sciprint("Error: arguments in_gap_size must be positive\r\n");
	  return 0;
  }

  outptr = create_aperture( 5, *stk(l1), *stk(l2), *stk(l3), *stk(l4) );

  m5=1; n5=17;
  CreateVarFromPtr(5, "d", &m5, &n5, &outptr);   

  LhsVar(1) = 5;
  return 0;
}

int aperture_transform_int(char *fname) 
{
	int  m1, n1, l1, m2, n2, l2, number, number2,
		minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
	FIELD field; APERTURE AP;
	WaveFront  WF_in, WF_out;
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
	if ( (m2!=1)||(n2 != 17 ) )  {
		sciprint("Error: second arguments must be aperture \r\n");
		return 0;
	}
	AP=array2APERTURE(stk(l2));

	WF_out =APERTURE_transform(AP, WF_in);
	WaveFront_FIELD(WF_out,&field);

	outptr = field2array( field );
	CreateVarFromPtr(3, "d", &m1, &n1, &outptr);   
	LhsVar(1) = 3;
	return 0;
}

int ap_frame_int(char *fname) 
{
	int  m1, n1, l1, m2, n2, l2, minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
	APERTURE AP; ThreeFrame TF;
	double *outptr=NULL;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1); 
	if ( (m1!=1)||( n1 != 17) )  {
		sciprint("Error: first arguments must be aperture \r\n");
		return 0;
	}
	AP=array2APERTURE(stk(l1));

	GetRhsVar(2, "d", &m2, &n2, &l2); 
	if ( (m2!=1)||(n2 != 12 ) )  {
		sciprint("Error: second arguments must be 3D frame \r\n");
		return 0;
	}
	TF = array2ThreeFrame( stk(l2) );

	AP.TF=TF;
	outptr = APERTURE2array( AP );
	CreateVarFromPtr(3, "d", &m1, &n1, &outptr);   
	LhsVar(1) = 3;
	return 0;
}

