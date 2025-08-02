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

int PI_controller_int(char *fname) 
{
  int  m1, n1, l1, m2, n2, l2, m3=1, n3=2,
	  minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
  double pgain,igain,*outptr=NULL;

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  GetRhsVar(1, "d", &m1, &n1, &l1);  
  if ( *stk(l1) < 0 )  {
	  sciprint("Error: arguments proportional_gain must be positive\r\n");
	  return 0;
  }
  pgain= *stk(l1);

  GetRhsVar(2, "d", &m2, &n2, &l2);  
  if ( *stk(l2) < 0 )  {
	  sciprint("Error: arguments integral_gain must be positive\r\n");
	  return 0;
  }
  igain= *stk(l2);

  outptr = (double*) calloc( 2, sizeof(double) );
  outptr[0]=pgain; outptr[1]=igain;

  CreateVarFromPtr( 3, "d", &m3, &n3, &outptr );   

  LhsVar(1) = 3;
  return 0;
}

int ttm_PI_controller_int(char *fname) 
{
	int  m1, n1, l1, m2, n2, minlhs=1, maxlhs=1, minrhs=1, maxrhs=1;
	TTM_PI_controller ttm_controller;
	double *outptr=NULL;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	if ( m1 !=1 || n1!=2 )  {
		sciprint("Error: arguments  must be proportional_integral_controller\r\n");
		return 0;
	}

	ttm_controller=create_TTM_PI_controller(*stk(l1),*stk(l1+1));
	outptr = TTM_PI_controller2array(ttm_controller);

	m2=1; n2=total_znk_space(1)+3;
	CreateVarFromPtr( 2, "d", &m2, &n2, &outptr );   

	LhsVar(1) = 2;
	return 0;
}

int dm_PI_controller_int(char *fname) 
{
	int  m1, n1, l1, m2, n2, l2, m3, n3,
		minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
	int actuator_x_dim,actuator_y_dim,nactuator;
	double *outptr=NULL; DeformableMirror dm;
	DM_PI_controller dm_controller;

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
	if ( m2 !=1 || n2!=2 )  {
		sciprint("Error: second argument must be proportional integral controller\r\n");
		return 0;
	}

	dm_controller=create_DM_PI_controller(dm,*stk(l2),*stk(l2+1));
	outptr = DM_PI_controller2array(dm_controller);

	m3=1; n3=nactuator+4;
	CreateVarFromPtr( 3, "d", &m3, &n3, &outptr );   

	LhsVar(1) = 3;
	return 0;
}

int update_ttm_PI_controller_int(char *fname) 
{
	int  m1, n1, l1, m2, n2, l2, m3, n3, l3, nelem,
		minlhs=2, maxlhs=2, minrhs=3, maxrhs=3;
	double *ttm_controller_outptr,*command_outptr; 
	ZERNIKE input, TTM_commands;
	TTM_PI_controller ttm_controller;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	nelem=total_znk_space(1);
	GetRhsVar(1, "d", &m1, &n1, &l1);  
	if ( (m1!=1)||(n1 != nelem+3 ) )  {
		sciprint("Error: first arguments must be ttm PI controller\r\n");
		return 0;
	}
	ttm_controller=array2TTM_PI_controller(stk(l1));
    
	GetRhsVar(2, "d", &m2, &n2, &l2);  
	if ( m2 !=1 || n2!=nelem+1 )  {
		sciprint("Error: second argument must be tip_tilt (ie.zernike(1)) residuals\r\n");
		return 0;
	}
	input=array2ZERNIKE(stk(l2));

	GetRhsVar(3, "d", &m3, &n3, &l3);  
	if ( m3 !=1 || n3!=nelem+1 )  {
		sciprint("Error: third argument must be tip_tilt (ie.zernike(1)) commands\r\n");
		return 0;
	}
	TTM_commands=array2ZERNIKE(stk(l3));

	TTM_commands=update_TTM_PI_controller(&ttm_controller,input,TTM_commands);
	ttm_controller_outptr = TTM_PI_controller2array(ttm_controller);
	command_outptr = ZERNIKE2array(TTM_commands);

	CreateVarFromPtr( 4, "d", &m1, &n1, &ttm_controller_outptr );  
	CreateVarFromPtr( 5, "d", &m3, &n3, &command_outptr );  

	LhsVar(1) = 4;
	LhsVar(2) = 5;

	return 0;
}

int update_dm_PI_controller_int(char *fname) 
{
  int  m1, n1, l1, m2, n2, l2, m3, n3, l3,
	  minlhs=2, maxlhs=2, minrhs=3, maxrhs=3;
  int actuator_x_dim,actuator_y_dim,nactuator;
  DM_PI_controller dm_controller;
  PixelArray input,DM_commands;

  double *dm_controller_outptr=NULL;
  double *commands_outptr=NULL;

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  GetRhsVar(1, "d", &m1, &n1, &l1);  

  actuator_x_dim=(int)(*stk(l1)); actuator_y_dim=(int)(*(stk(l1)+1));
  nactuator=actuator_x_dim*actuator_y_dim;
  if ( (m1!=1)||(n1 != (nactuator+4)) )  {
	  sciprint("Error: first arguments must be dm PI controller\r\n");
	  return 0;
  }
  dm_controller=array2DM_PI_controller(stk(l1));

  GetRhsVar(2, "d", &m2, &n2, &l2);  
  if ( m2 !=1 || n2!=nactuator+2 )  {
	  sciprint("Error: second argument must be dm residuals (pixel array)\r\n");
	  return 0;
  }
  input=array2PixelArray(stk(l2));

  GetRhsVar(3, "d", &m3, &n3, &l3);  
  if ( m3 !=1 || n3 !=nactuator+2 )  {
	  sciprint("Error: third argument must be dm command (pixel array)\r\n");
	  return 0;
  }
  DM_commands=array2PixelArray(stk(l3));

  DM_commands=update_DM_PI_controller(&dm_controller,input,DM_commands);
  dm_controller_outptr = DM_PI_controller2array(dm_controller);
  commands_outptr = PixelArray2array(DM_commands);

  CreateVarFromPtr( 4, "d", &m1, &n1, &dm_controller_outptr );   
  CreateVarFromPtr( 5, "d", &m3, &n3, &commands_outptr );   

  LhsVar(1) = 4;
  LhsVar(2) = 5;

  return 0;
}
