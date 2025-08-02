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

int create_actuator_arr_int(char *fname) 
{
  int  m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4,  m5, n5, 
	   minlhs=1, maxlhs=1, minrhs=4, maxrhs=4;
  double *outptr=NULL;

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  GetRhsVar(1, "i", &m1, &n1, &l1);  
  if ( *istk(l1) < 0 )  {
	  sciprint("Error: arguments actuator_x_dim must be positive integer\r\n");
	  return 0;
  }

  GetRhsVar(2, "i", &m2, &n2, &l2); 
  if ( *istk(l2) < 0 )  {
	  sciprint("Error: arguments actuator_y_dim must be positive integer\r\n");
	  return 0;
  }

  GetRhsVar(3, "d", &m3, &n3, &l3); 
  if ( *stk(l3) < 0 )  {
	  sciprint("Error: arguments actuator_pitch must be positive\r\n");
	  return 0;
  }

  GetRhsVar(4, "d", &m4, &n4, &l4); 
  if ( *stk(l4) < 0 )  {
	  sciprint("Error: arguments actuator_velocity must be position interger\r\n");
	  return 0;
  }

  outptr = (double*)calloc( 4, sizeof(double) );
  outptr[0]= *istk(l1); outptr[1]= *istk(l2);
  outptr[2]= *stk(l3); outptr[3]= *stk(l4);

  m5=1; n5=4;
  CreateVarFromPtr(5, "d", &m5, &n5, &outptr);   

  LhsVar(1) = 5;  
  return 0;
}


int create_ideal_dm_int(char *fname) 
{
  int  m1, n1, l1, m2, n2, l2, m3, n3,
	  minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
  APERTURE AP;  DeformableMirror dm;
  double *outptr=NULL; 

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  GetRhsVar(1, "d", &m1, &n1, &l1);  
  if ( m1!=1|| n1!=17 )  {
	  sciprint("Error: first arguments must be aperture \r\n");
	  return 0;
  }
  AP=array2APERTURE(stk(l1));

  GetRhsVar(2, "d", &m2, &n2, &l2); 
  if ( m2!=1|| n2!=4 )  {
	  sciprint("Error: second arguments must be actuator array\r\n");
	  return 0;
  }

  dm=create_DeformableMirror( AP,(int)(*stk(l2)),
		  (int)(*stk(l2+1)), *stk(l2+2), *stk(l2+3));

  outptr = DeformableMirror2array(dm);

  m3=1; n3=2*(dm.actuator_x_dim)*(dm.actuator_y_dim)+22;
  CreateVarFromPtr(3, "d", &m3, &n3, &outptr);   

  LhsVar(1) = 3;
  return 0;
}

int set_dm_timestamp_int(char *fname) 
{
  int  m1, n1, l1, m2, n2, l2, actuator_x_dim,actuator_y_dim,
	  nactuator, minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
  DeformableMirror dm;  double *outptr=NULL;

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  GetRhsVar(1, "d", &m1, &n1, &l1);  
  actuator_x_dim=(int)(*stk(l1)); actuator_y_dim=(int)(*(stk(l1)+1));
  nactuator=actuator_x_dim*actuator_y_dim;

  if ( (m1!=1)||(n1 != (2*nactuator+22)) )  {
	  sciprint("Error: first arguments must be deformable mirror\r\n");
	  return 0;
  }
  dm=array2DeformableMirror(stk(l1));

  GetRhsVar(2, "d", &m2, &n2, &l2); 
  
  set_dm_timestamp( &dm,*stk(l2) );

  outptr=DeformableMirror2array(dm);

  CreateVarFromPtr(3, "d", &m1, &n1, &outptr);   

  LhsVar(1) = 3;
  return 0;
}

int set_dm_actuator_positions_int(char *fname) 
{
  int  m1, n1, l1, m2, n2, l2, actuator_x_dim,actuator_y_dim,
	  nactuator, minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
  double *outptr=NULL; DeformableMirror dm;

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  GetRhsVar(1, "d", &m1, &n1, &l1);  
  actuator_x_dim=(int)(*stk(l1)); actuator_y_dim=(int)(*(stk(l1)+1));
  nactuator=actuator_x_dim*actuator_y_dim;

  if ( (m1!=1)||(n1 != (2*nactuator+22)) )  {
	  sciprint("Error: first arguments must be deformable mirror\r\n");
	  return 0;
  }

  dm=array2DeformableMirror(stk(l1));

  GetRhsVar(2, "d", &m2, &n2, &l2); 
  if ( (m2!=1)||(n2 != nactuator+2) )  {
	  sciprint("Error: first arguments must be actuator positions array\r\n");
	  return 0;
  }

  set_actuator_positions( &dm,stk(l2+2) );

  outptr=DeformableMirror2array(dm);

  CreateVarFromPtr(3, "d", &m1, &n1, &outptr);   

  LhsVar(1) = 3;
  return 0;
}


int set_dm_actuator_commands_int(char *fname) 
{
  int  m1, n1, l1, m2, n2, l2, actuator_x_dim,actuator_y_dim,
	  nactuator, minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
  double *outptr=NULL; DeformableMirror dm;

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  GetRhsVar(1, "d", &m1, &n1, &l1);  
  actuator_x_dim=(int)(*stk(l1)); actuator_y_dim=(int)(*(stk(l1)+1));
  nactuator=actuator_x_dim*actuator_y_dim;

  if ( (m1!=1)||(n1 != (2*nactuator+22)) )  {
	  sciprint("Error: first arguments must be deformable mirror\r\n");
	  return 0;
  }

  dm=array2DeformableMirror(stk(l1));

  GetRhsVar(2, "d", &m2, &n2, &l2); 
  if ( (m2!=1)||(n2 != nactuator+2) )  {
	  sciprint("Error: second arguments must be actuator commands array\r\n");
	  return 0;
  }

  set_actuator_commands(&dm,stk(l2+2));

  outptr=DeformableMirror2array(dm);

  CreateVarFromPtr(3, "d", &m1, &n1, &outptr);   

  LhsVar(1) = 3;
  return 0;
}

int DeformableMirror_update_int(char *fname) 
{
  int  m1, n1, l1, m2, n2, l2, m3, n3, l3, actuator_x_dim,actuator_y_dim,
	  nactuator, minlhs=1, maxlhs=1, minrhs=3, maxrhs=3;
  DeformableMirror dm;  double *outptr=NULL;

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  GetRhsVar(1, "d", &m1, &n1, &l1);  
  actuator_x_dim=(int)(*stk(l1)); actuator_y_dim=(int)(*(stk(l1)+1));
  nactuator=actuator_x_dim*actuator_y_dim;

  if ( (m1!=1)||(n1 != (2*nactuator+22)) )  {
	  sciprint("Error: first arguments must be deformable mirror\r\n");
	  return 0;
  }
  dm=array2DeformableMirror(stk(l1));

  GetRhsVar(2, "d", &m2, &n2, &l2); 
  if ( (m2!=1)||(n2 != nactuator+2) )  {
	  sciprint("Error: second arguments must be actuator commands array\r\n");
	  return 0;
  }

  GetRhsVar(3, "d", &m3, &n3, &l3); 
  
  DeformableMirror_updata( &dm,stk(l2+2),*stk(l3) );

  outptr=DeformableMirror2array(dm);
  CreateVarFromPtr(3, "d", &m1, &n1, &outptr);   
  LhsVar(1) = 3;

  return 0;
}

int ideal_dm_transform_int(char *fname) 
{
	int  m1, n1, l1, m2, n2, l2, actuator_x_dim,actuator_y_dim, nactuator,
		number, number2, minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
	double *outptr=NULL; FIELD field;
	WaveFront  WF_in, WF_out;
	DeformableMirror DM;
	
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
	actuator_x_dim=(int)(*stk(l2)); actuator_y_dim=(int)(*(stk(l2)+1));
	nactuator=actuator_x_dim*actuator_y_dim;
	if ( (m2!=1)||(n2 != (2*nactuator+22)) )  {
		sciprint("Error: second arguments must be deformable mirror\r\n");
		return 0;
	}
	DM=array2DeformableMirror(stk(l2));

	WF_out =DeformableMirror_transform(&DM, WF_in);
	WaveFront_FIELD(WF_out,&field);

	outptr = field2array( field );

	CreateVarFromPtr(3, "d", &m1, &n1, &outptr);  
	LhsVar(1) = 3;
	return 0;
}

