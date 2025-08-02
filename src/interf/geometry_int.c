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

int three_point_int(char *fname) 
{
  int m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, 
	  minlhs=1, maxlhs=1, minrhs=3, maxrhs=3;
  double *outptr;

  CheckRhs(minrhs,maxrhs);
  CheckLhs(minlhs,maxlhs);

  GetRhsVar(1, "d", &m1, &n1, &l1);  
  GetRhsVar(2, "d", &m2, &n2, &l2);  
  GetRhsVar(3, "d", &m3, &n3, &l3);  

  outptr = (double*) calloc( 3, sizeof(double) );
  outptr[0]=*stk(l1); outptr[1]=*stk(l2); outptr[2]=*stk(l3);

  m4=1; n4=3;
  CreateVarFromPtr(4, "d", &m4, &n4, &outptr);   

  LhsVar(1) = 4;
  return 0;
}

int three_vector_int(char *fname) 
{
  int m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, 
	  minlhs=1, maxlhs=1, minrhs=3, maxrhs=3;
  double *outptr;

  CheckRhs(minrhs,maxrhs);
  CheckLhs(minlhs,maxlhs);

  GetRhsVar(1, "d", &m1, &n1, &l1);  
  GetRhsVar(2, "d", &m2, &n2, &l2);  
  GetRhsVar(3, "d", &m3, &n3, &l3);  

  outptr = (double*) calloc( 3, sizeof(double) );
  outptr[0]=*stk(l1); outptr[1]=*stk(l2); outptr[2]=*stk(l3);

  m4=1; n4=3;
  CreateVarFromPtr(4, "d", &m4, &n4, &outptr);   

  LhsVar(1) = 4;
  return 0;
}

int three_frame_int(char *fname) 
{
	int m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4, m5, n5,
		minlhs=1, maxlhs=1, minrhs=0, maxrhs=4;
	double tp[3]={0,0,0}, tvx[3]={1,0,0}, tvy[3]={0,1,0}, 
		tvz[3]={0,0,1}, *outptr;
	ThreePoint TP; ThreeVector TVx, TVy, TVz; ThreeFrame TF;

	CheckRhs(minrhs,maxrhs);
	CheckLhs(minlhs,maxlhs);

	if (Rhs>0) {
		GetRhsVar(1, "d", &m1, &n1, &l1);
		if ( (m1!=1)||(n1 != 3 ) )  {
			sciprint("Error: first arguments must be three point\r\n");
			return 0;
		}
		tp[0]=*stk(l1);tp[1]=*stk(l1+1);tp[2]=*stk(l1+2);
	}
	TP=construct_ThreePoint(tp[0],tp[1],tp[2]);

	if (Rhs>1) {
		GetRhsVar(2, "d", &m2, &n2, &l2);
		if ( (m2!=1)||(n2 != 3 ) )  {
			sciprint("Error: second arguments must be three vector\r\n");
			return 0;
		}
		tvx[0]=*stk(l2);tvx[1]=*stk(l2+1);tvx[2]=*stk(l2+2);
	}
	TVx=construct_ThreeVector(tvx[0],tvx[1],tvx[2]);
    
	if (Rhs>2) {
		GetRhsVar(3, "d", &m3, &n3, &l3);
		if ( (m3!=1)||(n3 != 3 ) )  {
			sciprint("Error: third arguments must be three vector\r\n");
			return 0;
		}
		tvy[0]=*stk(l3);tvy[1]=*stk(l3+1);tvy[2]=*stk(l3+2);
	}
	TVy=construct_ThreeVector(tvy[0],tvy[1],tvy[2]);

	if (Rhs>3) {
		GetRhsVar(4, "d", &m4, &n4, &l4);
		if ( (m4!=1)||(n4 != 3 ) )  {
			sciprint("Error: fouth arguments must be three vector\r\n");
			return 0;
		}
		tvz[0]=*stk(l4);tvz[1]=*stk(l4+1);tvz[2]=*stk(l4+2);
	}
	TVz=construct_ThreeVector(tvz[0],tvz[1],tvz[2]);

	TF=construct_ThreeFrame(TP,TVx,TVy,TVz);
	outptr = ThreeFrame2array( TF );

	m5=1; n5=12;
	CreateVarFromPtr(5, "d", &m5, &n5, &outptr);   

	LhsVar(1) = 5;
	return 0;
}

int three_reflection_int(char *fname)
{
	int m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, 
		minlhs=1, maxlhs=1, minrhs=3, maxrhs=3;
	ThreePoint TP, tf_tp; ThreeFrame TF, TF_out;
	ThreeVector TV, tf_tvx,tf_tvy,tf_tvz;
	double *outptr;
	
	CheckRhs(minrhs,maxrhs);
	CheckLhs(minlhs,maxlhs);

	GetRhsVar(1, "d", &m1, &n1, &l1); 
	if ( (m1!=1)||(n1 != 12 ) )  {
		sciprint("Error: first arguments must be three frame\r\n");
		return 0;
	}

	tf_tp=construct_ThreePoint(*stk(l1),*stk(l1+1),*stk(l1+2));

	tf_tvx=construct_ThreeVector(*stk(l1+3),*stk(l1+4),*stk(l1+5));
	tf_tvy=construct_ThreeVector(*stk(l1+6),*stk(l1+7),*stk(l1+8));
	tf_tvz=construct_ThreeVector(*stk(l1+9),*stk(l1+10),*stk(l1+11));

	TF=construct_ThreeFrame ( tf_tp,tf_tvx,tf_tvy,tf_tvz );

	GetRhsVar(2, "d", &m2, &n2, &l2);  
	if ( (m2!=1)||(n2 != 3 ) ) {
		sciprint("Error: first arguments must be three point\r\n");
		return 0;
	}

	TP=construct_ThreePoint(*stk(l2),*stk(l2+1),*stk(l2+2));

	GetRhsVar(3, "d", &m3, &n3, &l3);  
	if ( (m3!=1)||(n3 != 3 ) ) {
		sciprint("Error: first arguments must be three vector\r\n");
		return 0;
	}

	TV=construct_ThreeVector(*stk(l3),*stk(l3+1),*stk(l3+2));

	TF_out=ThreeReflection(TF,TP,TV);
	outptr = ThreeFrame2array(TF_out);

	m4=1; n4=12;
	CreateVarFromPtr(4, "d", &m4, &n4, &outptr); 
	LhsVar(1) = 4;
	return 0;
}

int three_translation_int(char *fname)
{
	int m1, n1, l1, m2, n2, l2, m3, n3, 
		minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
	ThreePoint tf_tp; ThreeFrame TF, TF_out;
	ThreeVector TV, tf_tvx,tf_tvy,tf_tvz;
	double *outptr;
	
	CheckRhs(minrhs,maxrhs);
	CheckLhs(minlhs,maxlhs);

	GetRhsVar(1, "d", &m1, &n1, &l1); 
	if ( ( m1!=1)||( n1 != 12 ) )  {
		sciprint("Error: first arguments must be three frame\r\n");
		return 0;
	}

	tf_tp=construct_ThreePoint(*stk(l1),*stk(l1+1),*stk(l1+2));

	tf_tvx=construct_ThreeVector(*stk(l1+3),*stk(l1+4),*stk(l1+5));
	tf_tvy=construct_ThreeVector(*stk(l1+6),*stk(l1+7),*stk(l1+8));
	tf_tvz=construct_ThreeVector(*stk(l1+9),*stk(l1+10),*stk(l1+11));

	TF=construct_ThreeFrame ( tf_tp,tf_tvx,tf_tvy,tf_tvz );

	GetRhsVar(2, "d", &m2, &n2, &l2);  
	if ( ( m2!=1)||( n2 != 3 ) ) {
		sciprint("Error: first arguments must be three vector\r\n");
		return 0;
	}

	TV=construct_ThreeVector(*stk(l2),*stk(l2+1),*stk(l2+2));

	TF_out=ThreeTranslaton(TF, TV);
	outptr = ThreeFrame2array(TF_out);

	m3=1; n3=12;
	CreateVarFromPtr(3, "d", &m3, &n3, &outptr); 
	LhsVar(1) = 3;
	return 0;
}

int three_rotation_int(char *fname)
{
	int m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4, 
		m5, n5, minlhs=1, maxlhs=1, minrhs=4, maxrhs=4;
	ThreePoint TP, tf_tp; ThreeFrame TF, TF_out;
	ThreeVector TV, tf_tvx,tf_tvy,tf_tvz;
	double angle, *outptr;
	
	CheckRhs(minrhs,maxrhs);
	CheckLhs(minlhs,maxlhs);

	GetRhsVar(1, "d", &m1, &n1, &l1); 
	if ( ( m1!=1)||( n1 != 12 ) )  {
		sciprint("Error: first arguments must be three frame\r\n");
		return 0;
	}

	tf_tp=construct_ThreePoint(*stk(l1),*stk(l1+1),*stk(l1+2));

	tf_tvx=construct_ThreeVector(*stk(l1+3),*stk(l1+4),*stk(l1+5));
	tf_tvy=construct_ThreeVector(*stk(l1+6),*stk(l1+7),*stk(l1+8));
	tf_tvz=construct_ThreeVector(*stk(l1+9),*stk(l1+10),*stk(l1+11));

	TF=construct_ThreeFrame ( tf_tp,tf_tvx,tf_tvy,tf_tvz );

	GetRhsVar(2, "d", &m2, &n2, &l2);  
	if ( (m2!=1)||(n2 != 3 ) ) {
		sciprint("Error: first arguments must be three point\r\n");
		return 0;
	}

	TP=construct_ThreePoint(*stk(l2),*stk(l2+1),*stk(l2+2));

	GetRhsVar(3, "d", &m3, &n3, &l3);  
	if ( (m3!=1)||(n3 != 3 ) ) {
		sciprint("Error: first arguments must be three vector\r\n");
		return 0;
	}

	TV=construct_ThreeVector(*stk(l3),*stk(l3+1),*stk(l3+2));

	GetRhsVar(4, "d", &m4, &n4, &l4);  
	angle = *stk(l4);

	TF_out=ThreeRotation(TF,TP,TV,angle);
	outptr = ThreeFrame2array(TF_out);

	m5=1; n5=12;
	CreateVarFromPtr(5, "d", &m5, &n5, &outptr); 
	LhsVar(1) = 5;
	return 0;
}
