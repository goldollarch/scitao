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

int AtmosphericModel_int(char *fname) 
{
  int i, nlayers;  
  double *outptr=NULL;  ThreeFrame TF;
  int  m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4, m5, n5,
	  minlhs=1, maxlhs=1, minrhs=3, maxrhs=4;
  AtmosphericModel AtmModel;  PowerSpectrum *power_spectra;

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  GetRhsVar(1, "i", &m1, &n1, &l1);  
  nlayers = *istk(l1);
  if ( nlayers < 0 ) {
	  sciprint("Error: arguments nlayers must be positive\r\n");
	  return 0;
  }

  GetRhsVar(2, "d", &m2, &n2, &l2);
  if ( m2 !=1 || n2 != nlayers ) {
	  sciprint("Error: arguments heights must be %d position \r\n",nlayers);
	  return 0;
  }

  GetRhsVar(3, "d", &m3, &n3, &l3);  
  if ( m3 != 1 || n3 != 7*nlayers )  {
	  sciprint("Error: second arguments must be %d power spectra\r\n",nlayers);
	  return 0;
  }

  power_spectra=(PowerSpectrum *)calloc( nlayers, sizeof( PowerSpectrum) );
  for(i=0;i<nlayers;i++) power_spectra[i]=array2PowerSpectrum( stk(l3+7*i) );

  TF=default_ThreeFrame();
  if (Rhs>3) {
	  GetRhsVar(4, "d", &m4, &n4, &l4); 
	  if ( m4 != 1 || n4 != 12 )  {
		  sciprint("Error: fouth arguments must be three_frame\r\n");
		  return 0;
	  }
	  TF = array2ThreeFrame( stk(l4) );
  }

  AtmModel = construct_AtmosphericModel( *istk(l1),power_spectra,stk(l2),TF );
  outptr = AtmosphericModel2array( AtmModel );

  m5=1; n5=8*nlayers+21;
  CreateVarFromPtr(5, "d", &m5, &n5, &outptr);   

  LhsVar(1) = 5;
  return 0;
}

int Ellerbroek_Cerro_Pachon_model_int(char *fname) 
{
  int  m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4,
	  minlhs=1, maxlhs=1, minrhs=2, maxrhs=3;
  AtmosphericModel AtmModel; 
  int nlayers=0;  
  double *outptr=NULL;  ThreeFrame TF;

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  GetRhsVar(1, "d", &m1, &n1, &l1);  
  if ( *stk(l1)<0 ) {
	  sciprint("Error: arguments r_0_meters must be positive\r\n");
	  return 0;
  }

  GetRhsVar(2, "d", &m2, &n2, &l2);
  if ( *stk(l2)<0 ) {
	  sciprint("Error: arguments r_0_ref_wavelength_meters must be position \r\n",nlayers);
	  return 0;
  }

  TF=default_ThreeFrame();
  if (Rhs>2) {
	  GetRhsVar(3, "d", &m3, &n3, &l3); 
	  if ( m3 != 1 || n3 != 12 )  {
		  sciprint("Error: third arguments must be three_frame\r\n");
		  return 0;
	  }
	  TF = array2ThreeFrame( stk(l3) );
  }

  AtmModel = construct_Ellerbroek_Cerro_Pachon_model ( TF, *stk(l1), *stk(l2) );
  outptr = AtmosphericModel2array( AtmModel );

  m4=1; n4=8*nlayers+21;
  CreateVarFromPtr(4, "d", &m4, &n4, &outptr);   

  LhsVar(1) = 4;
  return 0;
}

int Ellerbroek_Mauna_Kea_model_int(char *fname) 
{
  int  m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4,
	  minlhs=1, maxlhs=1, minrhs=2, maxrhs=3;
  AtmosphericModel AtmModel; 
  int nlayers=0;  
  double *outptr=NULL;  ThreeFrame TF;

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  GetRhsVar(1, "d", &m1, &n1, &l1);  
  if ( *stk(l1)<0 ) {
	  sciprint("Error: arguments r_0_meters must be positive\r\n");
	  return 0;
  }

  GetRhsVar(2, "d", &m2, &n2, &l2);
  if ( *stk(l2)<0 ) {
	  sciprint("Error: arguments r_0_ref_wavelength_meters must be position \r\n",nlayers);
	  return 0;
  }

  TF=default_ThreeFrame();
  if (Rhs>2) {
	  GetRhsVar(3, "d", &m3, &n3, &l3); 
	  if ( m3 != 1 || n3 != 12 )  {
		  sciprint("Error: third arguments must be three_frame\r\n");
		  return 0;
	  }
	  TF = array2ThreeFrame( stk(l3) );
  }

  AtmModel = construct_Ellerbroek_Mauna_Kea_model ( TF, *stk(l1), *stk(l2) );
  outptr = AtmosphericModel2array( AtmModel );

  m4=1; n4=8*nlayers+21;
  CreateVarFromPtr(4, "d", &m4, &n4, &outptr);   

  LhsVar(1) = 4;
  return 0;
}

int Palomar_DIMM_MASS_model_int(char *fname) 
{
  int  m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4,
	  minlhs=1, maxlhs=1, minrhs=2, maxrhs=3;
  AtmosphericModel AtmModel; 
  int nlayers=0;  
  double *outptr=NULL;  ThreeFrame TF;

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  GetRhsVar(1, "d", &m1, &n1, &l1);  
  if ( *stk(l1)<0 ) {
	  sciprint("Error: arguments r_0_meters must be positive\r\n");
	  return 0;
  }

  GetRhsVar(2, "d", &m2, &n2, &l2);
  if ( *stk(l2)<0 ) {
	  sciprint("Error: arguments r_0_ref_wavelength_meters must be position \r\n",nlayers);
	  return 0;
  }

  TF=default_ThreeFrame();
  if (Rhs>2) {
	  GetRhsVar(3, "d", &m3, &n3, &l3); 
	  if ( m3 != 1 || n3 != 12 )  {
		  sciprint("Error: third arguments must be three_frame\r\n");
		  return 0;
	  }
	  TF = array2ThreeFrame( stk(l3) );
  }

  AtmModel = construct_Palomar_DIMM_MASS_model ( TF, *stk(l1), *stk(l2) );
  outptr = AtmosphericModel2array( AtmModel );

  m4=1; n4=8*nlayers+21;
  CreateVarFromPtr(4, "d", &m4, &n4, &outptr);   

  LhsVar(1) = 4;
  return 0;
}

int Hufnagel_Valley_model_int(char *fname) 
{
  int  m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4,
	  m5, n5, minlhs=1, maxlhs=1, minrhs=3, maxrhs=4;
  AtmosphericModel AtmModel;  int nlayers;  
  double *outptr=NULL;  ThreeFrame TF;

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  GetRhsVar(1, "i", &m1, &n1, &l1);  
  nlayers = *istk(l1);
  if ( nlayers<0 ) {
	  sciprint("Error: arguments nlayers must be positive\r\n");
	  return 0;
  }

  GetRhsVar(2, "d", &m2, &n2, &l2);
  if ( m2 !=1 || n2 != nlayers ) {
	  sciprint("Error: arguments heights must be %d position \r\n",nlayers);
	  return 0;
  }

  GetRhsVar(3, "d", &m3, &n3, &l3);  

  TF=default_ThreeFrame();
  if (Rhs>3) {
	  GetRhsVar(4, "d", &m4, &n4, &l4); 
	  if ( m4 != 1 || n4 != 12 )  {
		  sciprint("Error: third arguments must be three_frame\r\n");
		  return 0;
	  }
	  TF = array2ThreeFrame( stk(l4) );
  }

  AtmModel = construct_Hufnagel_Valley_model ( TF, nlayers, stk(l2), *stk(l3) );
  outptr = AtmosphericModel2array( AtmModel );

  m5=1; n5=8*nlayers+21;
  CreateVarFromPtr(5, "d", &m5, &n5, &outptr);   

  LhsVar(1) = 5;
  return 0;
}

int SLCSAT_day_model_int(char *fname) 
{
  int  m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4,
	  minlhs=1, maxlhs=1, minrhs=2, maxrhs=3;
  AtmosphericModel AtmModel;  int nlayers;  
  double *outptr=NULL;  ThreeFrame TF;

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  GetRhsVar(1, "i", &m1, &n1, &l1);  
  nlayers = *istk(l1);
  if ( nlayers<0 ) {
	  sciprint("Error: arguments nlayers must be positive\r\n");
	  return 0;
  }

  GetRhsVar(2, "d", &m2, &n2, &l2);
  if ( m2 !=1 || n2 != nlayers ) {
	  sciprint("Error: arguments heights must be %d position \r\n",nlayers);
	  return 0;
  }

  TF=default_ThreeFrame();
  if (Rhs>2) {
	  GetRhsVar(3, "d", &m3, &n3, &l3); 
	  if ( m3 != 1 || n3 != 12 )  {
		  sciprint("Error: third arguments must be three_frame\r\n");
		  return 0;
	  }
	  TF = array2ThreeFrame( stk(l3) );
  }

  AtmModel = construct_SLCSAT_day_model ( TF, nlayers, stk(l2) );
  outptr = AtmosphericModel2array( AtmModel );

  m4=1; n4=8*nlayers+21;
  CreateVarFromPtr(4, "d", &m4, &n4, &outptr);   

  LhsVar(1) = 4;
  return 0;
}

int SLCSAT_night_model_int(char *fname) 
{
  int  m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4,
	  minlhs=1, maxlhs=1, minrhs=2, maxrhs=3;
  AtmosphericModel AtmModel;  int nlayers;  
  double *outptr=NULL;  ThreeFrame TF;

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  GetRhsVar(1, "i", &m1, &n1, &l1);  
  nlayers = *istk(l1);
  if ( nlayers<0 ) {
	  sciprint("Error: arguments nlayers must be positive\r\n");
	  return 0;
  }

  GetRhsVar(2, "d", &m2, &n2, &l2);
  if ( m2 !=1 || n2 != nlayers ) {
	  sciprint("Error: arguments heights must be %d position \r\n",nlayers);
	  return 0;
  }

  TF=default_ThreeFrame();
  if (Rhs>2) {
	  GetRhsVar(3, "d", &m3, &n3, &l3); 
	  if ( m3 != 1 || n3 != 12 )  {
		  sciprint("Error: third arguments must be three_frame\r\n");
		  return 0;
	  }
	  TF = array2ThreeFrame( stk(l3) );
  }

  AtmModel = construct_SLCSAT_night_model ( TF, nlayers, stk(l2) );
  outptr = AtmosphericModel2array( AtmModel );

  m4=1; n4=8*nlayers+21;
  CreateVarFromPtr(4, "d", &m4, &n4, &outptr);   

  LhsVar(1) = 4;
  return 0;
}

int TMT_SRD_v13_Cn2_model_int(char *fname) 
{
  int  m1, n1, l1, m2, n2, 
	  minlhs=1, maxlhs=1, minrhs=0, maxrhs=1;
  AtmosphericModel AtmModel;  
  int nlayers=0;  
  double *outptr=NULL;  ThreeFrame TF;

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  TF=default_ThreeFrame();
  if (Rhs==1) {
	  GetRhsVar(1, "d", &m1, &n1, &l1); 
	  if ( m1 != 1 || n1 != 12 )  {
		  sciprint("Error: arguments must be three_frame\r\n");
		  return 0;
	  }
	  TF = array2ThreeFrame( stk(l1) );
  }

  AtmModel = construct_TMT_SRD_v13_Cn2_model ( TF );
  outptr = AtmosphericModel2array( AtmModel );

  m2=1; n2=8*nlayers+21;
  CreateVarFromPtr(2, "d", &m2, &n2, &outptr);   

  LhsVar(1) = 2;
  return 0;
}

int Gemini_GLAO_study_model_int(char *fname) 
{
  int  m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4, m5, n5, l5,
	  m6, n6, minlhs=1, maxlhs=1, minrhs=4, maxrhs=5;
  AtmosphericModel AtmModel;  ThreeFrame TF;
  double *outptr=NULL; 
  int nlayers=0;  

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  GetRhsVar(1, "i", &m1, &n1, &l1);  
  GetRhsVar(2, "i", &m2, &n2, &l2);
  GetRhsVar(3, "i", &m3, &n3, &l3);
  GetRhsVar(4, "d", &m4, &n4, &l4);
  if ( *stk(l4) < 0 ) {
	  sciprint("Error: arguments Gemini_outer_scale_meters must be positive\r\n");
	  return 0;
  }

  TF=default_ThreeFrame();
  if (Rhs>4) {
	  GetRhsVar(5, "d", &m5, &n5, &l5); 
	  if ( m5 != 1 || n5 != 12 )  {
		  sciprint("Error: third arguments must be three_frame\r\n");
		  return 0;
	  }
	  TF = array2ThreeFrame( stk(l5) );
  }

  AtmModel = construct_Gemini_GLAO_study_model 
	  ( TF, *istk(l1), *istk(l2), *istk(l3), *stk(l4) );
  outptr = AtmosphericModel2array( AtmModel );

  m6=1; n6=8*nlayers+21;
  CreateVarFromPtr(6, "d", &m6, &n6, &outptr);   

  LhsVar(1) = 6;
  return 0;
}

int AtmosphericModel_layer_number_int(char *fname) 
{
	int nlayers; int *outptr;
	AtmosphericModel AtmModel;

	int m1, n1, l1, m2, n2, minlhs=1, 
		maxlhs=1, minrhs=1, maxrhs=1;
	
	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1); 
	nlayers=(int)(*stk(l1));
	if (  m1 != 1 || n1 != 8*nlayers+21 ) {
		sciprint("Error: first arguments must be Atmospheric Model\r\n");
		return 0;
	}
	AtmModel = array2AtmosphericModel( stk(l1) );

	nlayers = AtmosphericModel_layer_number(AtmModel);
	outptr = (int*) malloc (sizeof(int));
	outptr[0] = nlayers;
	
	m2=1; n2=1;
	CreateVarFromPtr(2, "i", &m2, &n2, &outptr);   

	LhsVar(1) = 2;
	return 0;
}

int AtmosphericModel_layer_heights_int(char *fname) 
{
	int nlayers; double *outptr;
	AtmosphericModel AtmModel;

	int m1, n1, l1, m2, n2, minlhs=1, 
		maxlhs=1, minrhs=1, maxrhs=1;
	
	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1); 
	nlayers=(int)(*stk(l1));
	if (  m1 != 1 || n1 != 8*nlayers+21 ) {
		sciprint("Error: first arguments must be Atmospheric Model\r\n");
		return 0;
	}
	AtmModel = array2AtmosphericModel( stk(l1) );

	nlayers=AtmosphericModel_layer_number(AtmModel);
	outptr=AtmosphericModel_layer_heights(AtmModel);
	
	m2=1; n2=nlayers;
	CreateVarFromPtr(2, "d", &m2, &n2, &outptr);   

	LhsVar(1) = 2;
	return 0;
}

int get_AtmosphericModel_WavefrontHeader_int(char *fname) 
{
  int  m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4, m5, n5, l5,  m6, n6, l6, 
	  m7, n7, l7, m8, n8, minlhs=1, maxlhs=1, minrhs=7, maxrhs=7;
  AtmosphericModel AtmModel;  Emitter Emt;  APERTURE AP;
  double *outptr=NULL;   int nlayers;  
  WavefrontHeader WfH;

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  GetRhsVar(1, "d", &m1, &n1, &l1); 
  nlayers=(int)(*stk(l1));
  if (  m1 != 1 || n1 != 8*nlayers+21 ) {
	  sciprint("Error: first arguments must be Atmospheric Model\r\n");
	  return 0;
  }
  AtmModel = array2AtmosphericModel( stk(l1) );

  GetRhsVar(2, "d", &m2, &n2, &l2);
  if ( *stk(l2) < 0 ) {
	  sciprint("Error: second arguments wavelength must be positive\r\n");
	  return 0;
  }

  GetRhsVar(3, "d", &m3, &n3, &l3);
  if ( *stk(l3) < 0 ) {
	  sciprint("Error: third arguments pixscale must be positive\r\n");
	  return 0;
  }

  GetRhsVar(4, "d", &m4, &n4, &l4);
  if ( m4 != 1 || n4 != 4 ) {
	  sciprint("Error: fouth arguments must be Emitter\r\n");
	  return 0;
  }
  Emt = construct_Emitter( (int)(*stk(l4)), *stk(l4+1), *stk(l4+2), *stk(l4+3) );

  GetRhsVar(5, "d", &m5, &n5, &l5);
  if ( m5 != 1 || n5 != 17 )  {
	  sciprint("Error: fifth arguments must be aperture\r\n");
	  return 0;
  }
  AP = array2APERTURE( stk(l5) );

  GetRhsVar(6, "i", &m6, &n6, &l6);
  GetRhsVar(7, "i", &m7, &n7, &l7);

  WfH = get_AtmosphericModel_WavefrontHeader 
	  ( AtmModel, *stk(l2), *stk(l3), Emt, AP, *istk(l6), *istk(l7) );

  outptr = WavefrontHeader2array( WfH );

  m8=1; n8=18;
  CreateVarFromPtr(8, "d", &m8, &n8, &outptr);   

  LhsVar(1) = 8;
  return 0;
}

int get_AtmosphericModel_RefractiveLayer_int(char *fname) 
{
  int m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4, m5, n5, l5, m6, n6, l6, 
	  m7, n7, l7, m8, n8, l8, m9, n9, l9, i, nlayers, num_dwfhs,
	  minlhs=1, maxlhs=1, minrhs=9, maxrhs=9 ;
  double *tmp, *outptr, *layer_heights, *layer_pixscales;
  AtmosphericModel AtmModel;  SubharmonicMethod SubM; 
  WavefrontHeader *dwfhdrs; ThreeVector *layer_wind_vectors;
  RefractiveLayer *Rlayer;  HardyWind HWM;

  long  m10, n10, j, k, kk, total_Rlayer_space;

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  GetRhsVar(1, "d", &m1, &n1, &l1); 
  nlayers=(int)(*stk(l1));
  if (  m1 != 1 || n1 != 8*nlayers+21 ) {
	  sciprint("Error: first arguments must be Atmospheric Model\r\n");
	  return 0;
  }
  AtmModel = array2AtmosphericModel( stk(l1) );
  nlayers = AtmosphericModel_layer_number(AtmModel);
  layer_heights=AtmosphericModel_layer_heights(AtmModel);

  GetRhsVar(2, "d", &m2, &n2, &l2);
  if ( (m2 != 1)||(n2 != 4 ) ) {
	  sciprint("Error: second arguments wavelength must be Subharmonic method\r\n");
	  return 0;
  }
  SubM = array2SubharmonicMethod( stk(l2) );

  GetRhsVar(3, "i", &m3, &n3, &l3);
  num_dwfhs = *istk(l3);
  if (  num_dwfhs < 0 ) {
	  sciprint("Error: third arguments num_dwfhs must be positive\r\n");
	  return 0;
  }

  GetRhsVar(4, "d", &m4, &n4, &l4);
  if ( m4 != 1 || n4 != 18*num_dwfhs ) {
	  sciprint("Error: fouth arguments must be %d diffractive wavefront header\r\n",num_dwfhs);
	  return 0;
  }
  dwfhdrs=(WavefrontHeader*) calloc( num_dwfhs, sizeof( WavefrontHeader ) );
  for (i=0; i<num_dwfhs;i++) dwfhdrs[i]=array2WavefrontHeader( stk(l4+18*i) );

  GetRhsVar(5, "d", &m5, &n5, &l5);
  if ( m5 < 0 )  {
	  sciprint("Error: fifth arguments must be layer pixscale\r\n");
	  return 0;
  }
  layer_pixscales = (double*) malloc( nlayers*sizeof(double) );
  for(i=0;i<nlayers;i++) layer_pixscales[i]= *stk(l5);

  GetRhsVar(6, "d", &m6, &n6, &l6);
  if ( m6 != 1 || n6 != 4 )  {
	  sciprint("Error: sixth arguments must be Hardy wind model\r\n");
	  return 0;
  }

  HWM = construct_HardyWind( *stk(l6),*stk(l6+1),*stk(l6+2),*stk(l6+3) );
  layer_wind_vectors = get_HardyWind_velocities
	  ( HWM,nlayers,layer_heights, default_ThreeFrame( ) );

  GetRhsVar(7, "d", &m7, &n7, &l7);
  if (  *stk(l7) < 0 ) {
	  sciprint("Error: seventh arguments time_interval must be positive\r\n");
	  return 0;
  }

  GetRhsVar(8, "i", &m8, &n8, &l8);
  GetRhsVar(9, "i", &m9, &n9, &l9);

  Rlayer = get_AtmosphericModel_RefractiveLayer( AtmModel, SubM, num_dwfhs, 
	  dwfhdrs, layer_pixscales, layer_wind_vectors, *stk(l7), *istk(l8), *istk(l9) );

  total_Rlayer_space=1;
  for(i=0;i<nlayers;i++)  {
	  kk = ( Rlayer[i].PixArr.x_axes )*( Rlayer[i].PixArr.y_axes ) + 19;
	  total_Rlayer_space  += kk;
  }

  outptr = (double*) calloc( total_Rlayer_space, sizeof(double) );

  k=1;
  outptr[0]=nlayers;
  for(i=0;i<nlayers;i++) {
	  tmp=RefractiveLayer2array(Rlayer[i]);
	  kk = ( Rlayer[i].PixArr.x_axes )*( Rlayer[i].PixArr.y_axes )+19;
	  for(j=0;j<kk;j++) {
		  outptr[k]=tmp[j];  
		  k++;
	  }
  } 

  m10=1; n10=total_Rlayer_space;
  CreateVarFromPtr( 10, "d", &m10, &n10, &outptr );   

  LhsVar(1) = 10;
  return 0;
}

long total_atm_mod_Rlayer_space( double *inptr )
{
	int i, nlayers, x_axes, y_axes;
	long kk, total_Rlayer_space;

	nlayers=(int)(inptr[0]);

	total_Rlayer_space=1;
	for( i=0; i<nlayers; i++ ) {
        x_axes=(int)(inptr[total_Rlayer_space]);
		y_axes=(int)(inptr[total_Rlayer_space+1]);

		kk = x_axes*y_axes + 19;

		total_Rlayer_space += kk;
	}

	return total_Rlayer_space; 
}

int atm_mod_layers_fits_int(char *fname) 
{
  int  m2, n2, l2, minlhs=1, 
	  maxlhs=1, minrhs=2, maxrhs=2;
  long m1, n1, l1, k, total_Rlayer_space;
  RefractiveLayer *Rlayer;
  int i,nlayers;

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  GetRhsVar(1, "d", &m1, &n1, &l1); 
  total_Rlayer_space=total_atm_mod_Rlayer_space(stk(l1));
  if (  m1 != 1 || n1 != total_Rlayer_space ) {
	  sciprint("Error: first arguments must be Atmospheric Model refractive layers\r\n");
	  sciprint("n1 = %ld,  total_Rlayer_space=%ld\r\n",n1,total_Rlayer_space);
	  return 0;
  }

  nlayers=(int)(*stk(l1));
  Rlayer=(RefractiveLayer*) calloc( nlayers,sizeof(RefractiveLayer) );

  GetRhsVar(2, "c", &m2, &n2, &l2);  

  k=1;
  for( i=0; i<nlayers; i++ ) {
	  Rlayer[i]= array2RefractiveLayer( stk(l1+k) ) ;
	  k += (int)(*stk(l1+k)) *(int)(*stk(l1+k+1))+19;
	  write_atm_mod_ref_layer_file( Rlayer[i],i,cstk(l2) );
  }

  LhsVar(1) = 1;
  return 0;
}

int AtmosphericModel_to_APERTURE_transform_int(char *fname) 
{
	int  m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4, 
		minlhs=1, maxlhs=1, minrhs=4, maxrhs=4;
	int diffractive_near_field_propagator;
	int i, nlayers, number, number2;
	long k,total_Rlayer_space;
	double *outptr=NULL;

	FIELD field;
	RefractiveLayer *Rlayers;
	WaveFront  WF_in, WF_out;
	APERTURE AP;

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
	total_Rlayer_space=total_atm_mod_Rlayer_space(stk(l2));
	if (  m2 != 1 || n2 != total_Rlayer_space ) {
	  sciprint("Error: second arguments must be Atmospheric Model refractive layers\r\n");
	  sciprint("n2 = %ld,  total_Rlayer_space=%ld\r\n",n2,total_Rlayer_space);
	  return 0;
	}
	nlayers=(int)(*stk(l2));
	Rlayers=(RefractiveLayer*) calloc( nlayers,sizeof(RefractiveLayer) );
	
	k=1;
	for( i=0; i<nlayers; i++ ) {
	  Rlayers[i] = array2RefractiveLayer( stk(l2+k) ) ;
	  k += (int)(*stk(l2+k)) *(int)(*stk(l2+k+1))+19 ;
	}

	GetRhsVar(3, "d", &m3, &n3, &l3); 
	if ( m3 != 1 || n3 != 17 )  {
		sciprint("Error: third arguments must be aperture \r\n");
		return 0;
	}
	AP=array2APERTURE(stk(l3));

	GetRhsVar(4, "i", &m4, &n4, &l4); 
	diffractive_near_field_propagator = *istk(l4);

	WF_out =AtmosphericModel_to_APERTURE_transform
		( nlayers, Rlayers, AP, diffractive_near_field_propagator, WF_in );
	WaveFront_FIELD( WF_out, &field );
	
	outptr = field2array( field );
	CreateVarFromPtr( 5, "d", &m1, &n1, &outptr ); 
	LhsVar(1) = 5;

	return 0;
}
