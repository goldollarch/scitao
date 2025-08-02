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

int create_RECONSTRUCTOR_int(char *fname) 
{
  int AP_nelem,DM_nelem,LenArray_nelem,ZNK_nelem;
  int m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4, m5, n5, l5,
	  m6, n6, l6, m7, n7, minlhs=1, maxlhs=1, minrhs=6, maxrhs=6;

  double subaperture_illumination_threshold, eigenvalue_threshold;
  APERTURE AP; DeformableMirror DM; LensletArray LenArray; 
  ZERNIKE projected_modes; RECONSTRUCTOR Reconstructor; 

  long total_space;double *outptr;

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  GetRhsVar(1, "d", &m1, &n1, &l1); 
  AP_nelem=17;
  if ( ( m1 != 1 )||( n1 != AP_nelem ) )  {
	  sciprint("Error: first arguments must be aperture \r\n");
	  return 0;
  }
  AP = array2APERTURE( stk(l1) );

  GetRhsVar(2, "d", &m2, &n2, &l2);
  DM_nelem=2*((int)(*stk(l2))*(int)(*(stk(l2)+1)))+22;
  if ( (m2!=1)||(n2 != DM_nelem ) ) {
	  sciprint("Error: second arguments must be deformable mirror\r\n");
	  return 0;
  }
  DM = array2DeformableMirror( stk(l2) );

  GetRhsVar(3, "d", &m3, &n3, &l3);
  LenArray_nelem=6;
  if ( m3 != 1|| n3 != LenArray_nelem )  {
	  sciprint("Error: third arguments must be lenslet array\r\n");
	  return 0;
  }
  LenArray=create_LensletArray( (int)(*stk(l3)), (int)(*stk(l3+1)),
		*stk(l3+2), *stk(l3+3), (int)(*stk(l3+4)), (int)(*stk(l3+5)) );

  GetRhsVar(4, "d", &m4, &n4, &l4);
  ZNK_nelem=total_znk_space((int)(*stk(l4)))+1;
  if ( m4!=1|| n4!=ZNK_nelem)  {
	  sciprint("Error: fouth arguments must be zernike polynomial \r\n");
	  return 0;
  }
  projected_modes = array2ZERNIKE( stk(l4) );

  GetRhsVar(5, "d", &m5, &n5, &l5);
  if (  *stk(l5) < 0 ) {
	  sciprint("Error: fifth arguments subaperture_illumination_threshold must be positive\r\n");
	  return 0;
  }
  subaperture_illumination_threshold = *stk(l5);

  GetRhsVar(6, "d", &m6, &n6, &l6);
  if (  *stk(l6) < 0 ) {
	  sciprint("Error: sixth arguments eigenvalue_threshold must be positive\r\n");
	  return 0;
  }

  eigenvalue_threshold = *stk(l6);

  Reconstructor = create_RECONSTRUCTOR ( AP, DM, LenArray,	projected_modes, 
	  subaperture_illumination_threshold, eigenvalue_threshold );

  outptr = RECONSTRUCTOR2array( Reconstructor );
  total_space=total_Reconstructor_array_space(Reconstructor);

  m7=1; n7=total_space;
  CreateVarFromPtr( 7, "d", &m7, &n7, &outptr );   

  LhsVar(1) = 7;
  return 0;
}

int arroyo_reconstructor_fits_int(char *fname) 
{
  int  m2, n2, l2, minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
  RECONSTRUCTOR Reconstructor; double *outptr=NULL; 
  long m1, n1, l1, nelem;

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  GetRhsVar(1, "d", &m1, &n1, &l1);  
  nelem=total_Reconstructor_array_nelem(stk(l1));
  if ( ( m1 != 1 )||( n1 != nelem ) )  {
	  sciprint("Error: first arguments must be arroyo reconstructor \r\n");
	  return 0;
  }
  Reconstructor=array2RECONSTRUCTOR(stk(l1));

  GetRhsVar(2, "c", &m2, &n2, &l2);  

  write_RECONSTRUCTOR_file( Reconstructor, cstk(l2) );

  LhsVar(1) = 1;
  return 0;
}

int arroyo_reconstruct_zernike_residuals_int(char *fname) 
{
  int  m2, n2, l2,m3, n3, minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
  RECONSTRUCTOR Reconstructor; ZERNIKE ZNK;
  SHartmannCentroids SHcentroids;
  long  m1, n1, l1, nelem;
  int pixarr_x, pixarr_y;
  double *outptr=NULL;

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  GetRhsVar(1, "d", &m1, &n1, &l1); 
  nelem=total_Reconstructor_array_nelem(stk(l1));
  if ( ( m1 != 1 )||( n1 != nelem ) )  {
	  sciprint("Error: first arguments must be arroyo reconstructor \r\n");
	  return 0;
  }
  Reconstructor=array2RECONSTRUCTOR(stk(l1));

  GetRhsVar(2, "d", &m2, &n2, &l2); 
  pixarr_x=(int)(*stk(l2)); pixarr_y=(int)(*stk(l2+1));
  if ( ( m2!=1 )||( n2 != (pixarr_x*pixarr_y+2 ) ) )  {
	  sciprint("Error: second arguments must be Shack Hartmann Centroids \r\n");
	  return 0;
  }
  SHcentroids=array2SHartmannCentroids( stk(l2) );

  ZNK=arroyo_reconstruct_zernike_residuals( Reconstructor,SHcentroids );

  m3=1; n3=total_znk_space(ZNK.order)+1;
  outptr = ZERNIKE2array( ZNK );

  CreateVarFromPtr(3, "d", &m3, &n3, &outptr);   

  LhsVar(1) = 3;
  return 0;
}

int arroyo_reconstruct_zonal_residuals_int(char *fname) 
{
  int  m2, n2, l2,m3, n3, minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
  RECONSTRUCTOR Reconstructor; PixelArray dm_residuals;
  SHartmannCentroids SHcentroids;
  long  m1, n1, l1, nelem;
  int pixarr_x, pixarr_y;
  double *outptr=NULL;

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  GetRhsVar(1, "d", &m1, &n1, &l1); 
  nelem=total_Reconstructor_array_nelem(stk(l1));
  if ( ( m1 != 1 )||( n1 != nelem ) )  {
	  sciprint("Error: first arguments must be arroyo reconstructor \r\n");
	  return 0;
  }
  Reconstructor=array2RECONSTRUCTOR(stk(l1));

  GetRhsVar(2, "d", &m2, &n2, &l2); 
  pixarr_x=(int)(*stk(l2)); pixarr_y=(int)(*stk(l2+1));
  if ( ( m2!=1 )||( n2 != (pixarr_x*pixarr_y+2 ) ) )  {
	  sciprint("Error: second arguments must be Shack Hartmann Centroids \r\n");
	  return 0;
  }
  SHcentroids=array2SHartmannCentroids( stk(l2) );

  dm_residuals=arroyo_reconstruct_zonal_residuals( Reconstructor,SHcentroids );

  m3=1; n3=(dm_residuals.x_axes*dm_residuals.y_axes)+2;
  outptr = PixelArray2array( dm_residuals );

  CreateVarFromPtr(3, "d", &m3, &n3, &outptr);   

  LhsVar(1) = 3;
  return 0;
}

