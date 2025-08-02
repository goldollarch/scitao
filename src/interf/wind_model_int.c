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

int Hardy_wind_model_int(char *fname) 
{
  int   m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4,
	  m5, n5, minlhs=1, maxlhs=1, minrhs=1, maxrhs=4;
  double tropopause_vel=0, tropopause_height=0, tropopause_thickness=0;
  double *outptr=NULL;

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  GetRhsVar(1, "d", &m1, &n1, &l1);  
  if ( *stk(l1) < 0 )  {
	  sciprint("Error: arguments ground_layer_wind_velocity must be position \r\n");
	  return 0;
  }

  if (Rhs>1) {
	  GetRhsVar(2, "d", &m2, &n2, &l2); 
	  if ( *stk(l2) < 0 )  {
		  sciprint("Error: arguments tropopause_wind_velocity must be position \r\n");
		  return 0;
	  }
	  tropopause_vel = *stk(l2);
  }

  if (Rhs>2) {
	  GetRhsVar(3, "d", &m3, &n3, &l3); 
	  if ( *stk(l3) < 0 )  {
		  sciprint("Error: arguments tropopause_height must be position \r\n");
		  return 0;
	  }
	  tropopause_height = *stk(l3);
  }

  if (Rhs>3) {
	  GetRhsVar(4, "d", &m4, &n4, &l4); 
	  if ( *stk(l4) < 0 )  {
		  sciprint("Error: arguments tropopause_thickness must be position \r\n");
		  return 0;
	  }
	  tropopause_thickness = *stk(l4);
  }

  outptr = (double*) calloc( 4, sizeof(double) );
  outptr[0] = *stk(l1);
  outptr[1]=tropopause_vel;
  outptr[2]=tropopause_height;
  outptr[3]=tropopause_thickness;

  m5=1; n5=4;
  CreateVarFromPtr(5, "d", &m5, &n5, &outptr);   

  LhsVar(1) = 5;
  return 0;
}


int get_HardyWind_velocities_int(char *fname) 
{
  int   m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4,
	  m5, n5, minlhs=1, maxlhs=1, minrhs=3, maxrhs=4;
  HardyWind HWM;  ThreeVector *wind_TVs;
  int i, nlayers; double *outptr=NULL;
  ThreeFrame TF=default_ThreeFrame();

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  GetRhsVar(1, "d", &m1, &n1, &l1);  
  if ( m1 !=1 || n1 !=4 ) {
	  sciprint("Error: first arguments must be hardy wind model\r\n");
	  return 0;
  }

  GetRhsVar(2, "i", &m2, &n2, &l2);
  nlayers = *istk(l2);
  if ( nlayers < 0 ) {
	  sciprint("Error: arguments nlayers must be positive\r\n");
	  return 0;
  }

  GetRhsVar(3, "d", &m3, &n3, &l3); 
  if ( m3 != 1 || n3 != nlayers )  {
	  sciprint("Error: arguments heights must be %d position \r\n",nlayers);
	  return 0;
  }

  if (Rhs>3) {
	  GetRhsVar(4, "d", &m4, &n4, &l4); 
	  if ( m4 != 1 || n4 != 12 )  {
		  sciprint("Error: fouth arguments must be three_frame\r\n");
		  return 0;
	  }
	  TF = array2ThreeFrame( stk(l4) );
  }

  HWM = construct_HardyWind( *stk(l1),*stk(l1+1),*stk(l1+2),*stk(l1+3) );
  wind_TVs = get_HardyWind_velocities( HWM, nlayers, stk(l3), TF );

  outptr = (double*) calloc( 3*nlayers, sizeof(double) );
  for( i=0;i<nlayers;i++) {
	  outptr[3*i]=wind_TVs[i].x;
	  outptr[3*i+1]=wind_TVs[i].y;
	  outptr[3*i+2]=wind_TVs[i].z;
  }

  m5=3; n5=nlayers;
  CreateVarFromPtr(5, "d", &m5, &n5, &outptr);   

  LhsVar(1) = 5;
  return 0;
}

