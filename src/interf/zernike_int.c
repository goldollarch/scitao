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

int create_zernike_int(char *fname) 
{
  int   m1, n1, l1, m2, n2, order,
	  minlhs=1, maxlhs=1, minrhs=1, maxrhs=1;
  double *outptr=NULL; ZERNIKE ZNK;

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  GetRhsVar(1, "i", &m1, &n1, &l1);  
  order= *istk(l1);
  if ( order <  0  )  {
      sciprint("Error: order must be positive\r\n");
      return 0;
  }
  ZNK=create_ZERNIKE( order );

  outptr = ZERNIKE2array( ZNK );
  m2=1; n2=total_znk_space( order )+1;
  CreateVarFromPtr(2, "d", &m2, &n2, &outptr);   
  LhsVar(1) = 2;
  return 0;
}

int get_cos_coef_int(char *fname) 
{
  int  m1, n1, l1, m2, n2, l2, m3, n3, l3, m4=1, n4=1, 
	  minlhs=1, maxlhs=1, minrhs=3, maxrhs=3;
  int nelem; ZERNIKE ZNK;
  double *outptr,coef;

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  GetRhsVar(1, "d", &m1, &n1, &l1); 
  nelem=total_znk_space((int)(*stk(l1)))+1;
  if ( m1!=1|| n1!=nelem)  {
	  sciprint("Error: first arguments must be zernike polynomial \r\n");
	  return 0;
  }
  ZNK = array2ZERNIKE( stk(l1) );

  GetRhsVar(2, "i", &m2, &n2, &l2); 
  GetRhsVar(3, "i", &m3, &n3, &l3);  
  
  coef=get_ZNK_cos_coeff ( ZNK, *istk(l2), *istk(l3) );
  outptr=(double*)malloc(sizeof(double));outptr[0]=coef;

  CreateVarFromPtr(4, "d", &m4, &n4, &outptr);   
  LhsVar(1) = 4;
  return 0;
}

int set_cos_coef_int(char *fname) 
{
  int  m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4, 
	  minlhs=1, maxlhs=1, minrhs=4, maxrhs=4;
  int nelem; ZERNIKE ZNK;
  double *outptr;

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  GetRhsVar(1, "d", &m1, &n1, &l1); 
  nelem=total_znk_space((int)(*stk(l1)))+1;
  if ( m1!=1|| n1!=nelem)  {
	  sciprint("Error: first arguments must be zernike polynomial \r\n");
	  return 0;
  }
  ZNK = array2ZERNIKE( stk(l1) );

  GetRhsVar(2, "i", &m2, &n2, &l2); 
  GetRhsVar(3, "i", &m3, &n3, &l3);  
  GetRhsVar(4, "d", &m4, &n4, &l4);  
  
  set_ZNK_cos_coef ( &ZNK, *istk(l2), *istk(l3), *stk(l4) );
  outptr = ZERNIKE2array( ZNK );
  CreateVarFromPtr(5, "d", &m1, &n1, &outptr);   
  LhsVar(1) = 5;
  return 0;
}

int get_sin_coef_int(char *fname) 
{
  int  m1, n1, l1, m2, n2, l2, m3, n3, l3, m4=1, n4=1, 
	  minlhs=1, maxlhs=1, minrhs=3, maxrhs=3;
  int nelem; ZERNIKE ZNK;
  double *outptr,coef;

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  GetRhsVar(1, "d", &m1, &n1, &l1); 
  nelem=total_znk_space((int)(*stk(l1)))+1;
  if ( m1!=1|| n1!=nelem)  {
	  sciprint("Error: first arguments must be zernike polynomial \r\n");
	  return 0;
  }
  ZNK = array2ZERNIKE( stk(l1) );

  GetRhsVar(2, "i", &m2, &n2, &l2); 
  GetRhsVar(3, "i", &m3, &n3, &l3);  
  
  coef=get_ZNK_sin_coeff ( ZNK, *istk(l2), *istk(l3) );
  outptr=(double*)malloc(sizeof(double));outptr[0]=coef;

  CreateVarFromPtr(4, "d", &m4, &n4, &outptr);   
  LhsVar(1) = 4;
  return 0;
}

int set_sin_coef_int(char *fname) 
{
  int  m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4, 
	  minlhs=1, maxlhs=1, minrhs=4, maxrhs=4;
  int nelem; ZERNIKE ZNK;
  double *outptr;

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  GetRhsVar(1, "d", &m1, &n1, &l1); 
  nelem=total_znk_space((int)(*stk(l1)))+1;
  if ( m1!=1|| n1!=nelem)  {
	  sciprint("Error: first arguments must be zernike polynomial \r\n");
	  return 0;
  }
  ZNK = array2ZERNIKE( stk(l1) );

  GetRhsVar(2, "i", &m2, &n2, &l2); 
  GetRhsVar(3, "i", &m3, &n3, &l3);  
  GetRhsVar(4, "d", &m4, &n4, &l4);  
  
  set_ZNK_sin_coef ( &ZNK, *istk(l2), *istk(l3), *stk(l4) );
  outptr = ZERNIKE2array( ZNK );
  CreateVarFromPtr(5, "d", &m1, &n1, &outptr);   
  LhsVar(1) = 5;
  return 0;
}

int zernike_fits_int(char *fname) 
{
  int  m1, n1, l1, m2, n2, l2, m3, n3, l3, 
	  minlhs=1, maxlhs=1, minrhs=2, maxrhs=3;
  int nelem; ZERNIKE ZNK;
  double time=-1;

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  GetRhsVar(1, "d", &m1, &n1, &l1); 
  nelem=total_znk_space((int)(*stk(l1)))+1;
  if ( m1!=1|| n1!=nelem)  {
	  sciprint("Error: first arguments must be zernike polynomial \r\n");
	  return 0;
  }
  ZNK = array2ZERNIKE( stk(l1) );

  GetRhsVar(2, "c", &m2, &n2, &l2); 

	if (Rhs>2) {
		GetRhsVar(3, "d", &m3, &n3, &l3); 
		time= *stk(l3);
	}

  zernike_write ( ZNK, cstk(l2), time );
 
  LhsVar(1) = 1;
  return 0;
}

