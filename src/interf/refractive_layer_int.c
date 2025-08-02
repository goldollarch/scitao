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

double *null_inner_scale()
{
	double *tmp;
	int inner_scale_type=0;

	tmp=(double *) calloc( 2, sizeof(double) );

	tmp[0]=inner_scale_type;
	tmp[1]=0;	

	return tmp;
}

double *null_subharmonic_method()
{
	double *tmp;
	int subharmonic_method_type=0;

	tmp=(double *) calloc( 4, sizeof(double) );

	tmp[0]=subharmonic_method_type;
	tmp[1]=0; tmp[2]=0; tmp[3]=0;	

	return tmp;
}

double *power_law( double exponent, double r_0_meters,
				  double r_0_ref_wavelength_meters )
{
	double *tmp;
	int power_law_type=0;

	tmp=(double *) calloc( 5, sizeof(double) );

	tmp[0]=power_law_type;
	tmp[1]=exponent;	
	tmp[2]=r_0_meters;
	tmp[3]=r_0_ref_wavelength_meters;	
	tmp[4]=0;	

	return tmp;
}

double *create_power_spectrum( double *power_law_data, 
							  double *inner_scale_data)
{
	double *tmp;
	int power_law_type,inner_scale_type;

	tmp=(double *) calloc( 7, sizeof(double) );

	power_law_type=(int)power_law_data[0];
	tmp[0]=power_law_type;
	tmp[1]=power_law_data[1]; tmp[2]=power_law_data[2];
	tmp[3]=power_law_data[3]; tmp[4]=power_law_data[4];

	inner_scale_type=(int)inner_scale_data[0];
	tmp[5]=inner_scale_type;	
	tmp[6]=inner_scale_data[1];

	return tmp;
}

int power_law_int(char *fname) 
{
  int   m1, n1, l1, m2, n2, l2, m3, n3, l3, m4=1, n4=5,
	  minlhs=1, maxlhs=1, minrhs=1, maxrhs=3;
  double exponent = -11/3.0;
  double *outptr=NULL;

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  GetRhsVar(1, "d", &m1, &n1, &l1);  
  if ( *stk(l1) < 0 )  {
	  sciprint("Error: arguments r_0_meters must be positive\r\n");
	  return 0;
  }

  GetRhsVar(2, "d", &m2, &n2, &l2);  
  if ( *stk(l2) < 0 )  {
	  sciprint("Error: arguments r_0_ref_wavelength_meters must be positive\r\n");
	  return 0;
  }

  if (Rhs>2) {
	  GetRhsVar(3, "d", &m3, &n3, &l3);
	  exponent=(*stk(l3));
  }

  outptr = power_law( exponent, *stk(l1),  *stk(l2) );

  CreateVarFromPtr(4, "d", &m4, &n4, &outptr);   

  LhsVar(1) = 4;
  return 0;
}

int von_karman_power_law_int(char *fname) 
{
	int   m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4, power_law_type=1,
		m5=1, n5=5,  minlhs=1, maxlhs=1, minrhs=3, maxrhs=4;
	double r_0_meters,r_0_ref_wavelength_meters,outer_scale;
	double exponent = -11/3.0, *outptr=NULL;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	if ( *stk(l1) < 0 )  {
		sciprint("Error: arguments r_0_meters must be positive\r\n");
		return 0;
	}
	r_0_meters= *stk(l1);

	GetRhsVar(2, "d", &m2, &n2, &l2);  
	if ( *stk(l2) < 0 )  {
		sciprint("Error: arguments r_0_ref_wavelength_meters must be positive\r\n");
		return 0;
	}
	r_0_ref_wavelength_meters= *stk(l2);

	GetRhsVar(3, "d", &m3, &n3, &l3);  
	if ( *stk(l3) < 0 )  {
		sciprint("Error: arguments outer_scale must be positive\r\n");
		return 0;
	}
	outer_scale= *stk(l3);

	if (Rhs>3) {
		GetRhsVar(4, "d", &m4, &n4, &l4);
		exponent=(*stk(l4));
	}

	outptr = (double *) calloc( 5, sizeof(double) );
	outptr[0]=power_law_type;	outptr[1]=exponent;	
	outptr[2]=r_0_meters;	outptr[3]=r_0_ref_wavelength_meters;	
	outptr[4]=outer_scale;	
    
	CreateVarFromPtr(5, "d", &m5, &n5, &outptr);   
    LhsVar(1) = 5;
	return 0;
}

int greenwood_power_law_int(char *fname) 
{
	int   m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4, power_law_type=2,
		m5=1, n5=5,  minlhs=1, maxlhs=1, minrhs=3, maxrhs=4;
	double r_0_meters,r_0_ref_wavelength_meters,outer_scale;
	double exponent = -11/3.0, *outptr=NULL;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	if ( *stk(l1) < 0 )  {
		sciprint("Error: arguments r_0_meters must be positive\r\n");
		return 0;
	}
	r_0_meters= *stk(l1);

	GetRhsVar(2, "d", &m2, &n2, &l2);  
	if ( *stk(l2) < 0 )  {
		sciprint("Error: arguments r_0_ref_wavelength_meters must be positive\r\n");
		return 0;
	}
	r_0_ref_wavelength_meters= *stk(l2);

	GetRhsVar(3, "d", &m3, &n3, &l3);  
	if ( *stk(l3) < 0 )  {
		sciprint("Error: arguments outer_scale must be positive\r\n");
		return 0;
	}
	outer_scale= *stk(l3);

	if (Rhs>3) {
		GetRhsVar(4, "d", &m4, &n4, &l4);
		exponent=(*stk(l4));
	}

	outptr = (double *) calloc( 5, sizeof(double) );
	outptr[0]=power_law_type;	outptr[1]=exponent;	
	outptr[2]=r_0_meters;	outptr[3]=r_0_ref_wavelength_meters;	
	outptr[4]=outer_scale;	
    
	CreateVarFromPtr(5, "d", &m5, &n5, &outptr);   
    LhsVar(1) = 5;
	return 0;
}

int null_inner_scale_int(char *fname) 
{
  int   m1=1, n1=2, minlhs=1, maxlhs=1, minrhs=0, maxrhs=0;
  double *outptr=NULL;

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  outptr = null_inner_scale();

  CreateVarFromPtr(1, "d", &m1, &n1, &outptr);   

  LhsVar(1) = 1;
  return 0;
}

int exponential_inner_scale_int(char *fname) 
{
	int   m1, n1, l1, m2=1, n2=2, inner_scale_type=1,
		minlhs=1, maxlhs=1, minrhs=1, maxrhs=1;
	double *outptr=NULL;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	if ( *stk(l1) < 0 )  {
		sciprint("Error: arguments inner_scale must be positive\r\n");
		return 0;
	}

	outptr = (double *) calloc( 2, sizeof(double) );
	outptr[0]=inner_scale_type;
	outptr[1]=*stk(l1);	
	
	CreateVarFromPtr(2, "d", &m2, &n2, &outptr);
	LhsVar(1) = 2;
	return 0;
}

int frehlich_inner_scale_int(char *fname) 
{
	int   m1, n1, l1, m2=1, n2=2, inner_scale_type=2,
		minlhs=1, maxlhs=1, minrhs=1, maxrhs=1;
	double *outptr=NULL;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "d", &m1, &n1, &l1);  
	if ( *stk(l1) < 0 )  {
		sciprint("Error: arguments inner_scale must be positive\r\n");
		return 0;
	}

	outptr = (double *) calloc( 2, sizeof(double) );
	outptr[0]=inner_scale_type;
	outptr[1]=*stk(l1);	
	
	CreateVarFromPtr(2, "d", &m2, &n2, &outptr);
	LhsVar(1) = 2;
	return 0;
}

int create_power_spectrum_int(char *fname) 
{
  int  m1, n1, l1, m2, n2, l2, m3=1, n3=7,
	  minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
  double *outptr=NULL;

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  GetRhsVar(1, "d", &m1, &n1, &l1);  
  if ( (m1!=1)||(n1 != 5 ) )  {
	  sciprint("Error: first arguments must be power law \r\n");
	  return 0;
  }

  GetRhsVar(2, "d", &m2, &n2, &l2);  
  if ( (m2!=1)||(n2 != 2 ) )  {
	  sciprint("Error: second arguments must be inner scale \r\n");
	  return 0;
  }

  outptr = create_power_spectrum( stk(l1), stk(l2) );

  CreateVarFromPtr( 3, "d", &m3, &n3, &outptr );   

  LhsVar(1) = 3;
  return 0;
}

int null_subharmonic_method_int(char *fname) 
{
  int   m1=1, n1=4, minlhs=1, maxlhs=1, minrhs=0, maxrhs=0;
  double *outptr=NULL;

  CheckRhs(minrhs,maxrhs) ;
  CheckLhs(minlhs,maxlhs) ;

  outptr = null_subharmonic_method( );

  CreateVarFromPtr(1, "d", &m1, &n1, &outptr);   

  LhsVar(1) = 1;
  return 0;
}

int quad_pixel_subharmonic_method_int(char *fname) 
{
	int   m1, n1, l1, m2=1, n2=4, subharmonic_method_type=1,
		minlhs=1, maxlhs=1, minrhs=1, maxrhs=1;
	double *outptr=NULL;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "i", &m1, &n1, &l1);  
	if ( *istk(l1) < 0 )  {
		sciprint("Error: arguments subharmonic_depth must be positive\r\n");
		return 0;
	}

	outptr = (double *) calloc( 4, sizeof(double) );
	outptr[0]=subharmonic_method_type;
	outptr[1]= *istk(l1); 
	outptr[2]=0; outptr[3]=0;	

	CreateVarFromPtr(2, "d", &m2, &n2, &outptr);   
	LhsVar(1) = 2;
	return 0;
}

int Lane_subharmonic_method_int(char *fname) 
{
	int   m1, n1, l1, m2=1, n2=4, subharmonic_method_type=2,
		minlhs=1, maxlhs=1, minrhs=1, maxrhs=1;
	double *outptr=NULL;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "i", &m1, &n1, &l1);  
	if ( *istk(l1) < 0 )  {
		sciprint("Error: arguments subharmonic_depth must be positive\r\n");
		return 0;
	}

	outptr = (double *) calloc( 4, sizeof(double) );
	outptr[0]=subharmonic_method_type;
	outptr[1]= *istk(l1); 
	outptr[2]=0; outptr[3]=0;	

	CreateVarFromPtr(2, "d", &m2, &n2, &outptr);   
	LhsVar(1) = 2;
	return 0;
}

int Johansson_Gavel_subharmonic_method_int(char *fname) 
{
	int   m1, n1, l1, m2=1, n2=4, subharmonic_method_type=3,
		minlhs=1, maxlhs=1, minrhs=1, maxrhs=1;
	double *outptr=NULL;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "i", &m1, &n1, &l1);  
	if ( *istk(l1) < 0 )  {
		sciprint("Error: arguments subharmonic_depth must be positive\r\n");
		return 0;
	}

	outptr = (double *) calloc( 4, sizeof(double) );
	outptr[0]=subharmonic_method_type;
	outptr[1]= *istk(l1); 
	outptr[2]=0; outptr[3]=0;	

	CreateVarFromPtr(2, "d", &m2, &n2, &outptr);   
	LhsVar(1) = 2;
	return 0;
}

int generalized_subharmonic_method_int(char *fname) 
{
	int   m1, n1, l1, m2, n2, l2, m3, n3, l3, m4=1, n4=4, 
		minlhs=1, maxlhs=1, minrhs=3, maxrhs=3;
	int subharmonic_depth,subpixels_per_level,
		subpixels_per_pixel,subharmonic_method_type=3;
	double *outptr=NULL;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "i", &m1, &n1, &l1);  
	if ( *istk(l1) < 0 )  {
		sciprint("Error: arguments subharmonic_depth must be positive\r\n");
		return 0;
	}
	subharmonic_depth= *istk(l1);

	GetRhsVar(2, "i", &m2, &n2, &l2);  
	if ( *istk(l2) < 3 || *istk(l2)%2==0 )  {
		sciprint("Error: subpixels_per_level must be  an odd number and >= 3\r\n");
		return 0;
	}
	subpixels_per_level= *istk(l2);

	GetRhsVar(3, "i", &m3, &n3, &l3);  
	if ( *istk(l3) < 3 || *istk(l3)%2==0 )  {
		sciprint("Error: subpixels_per_pixel must be  an odd number of and >= 3\r\n");
		return 0;
	}
	subpixels_per_pixel= *istk(l3);

	outptr = (double *) calloc( 4, sizeof(double) );
	outptr[0]=subharmonic_method_type;
	outptr[1]=subharmonic_depth; 
	outptr[2]=subpixels_per_level; 
	outptr[3]=subpixels_per_pixel;	

	CreateVarFromPtr(4, "d", &m4, &n4, &outptr);   
	LhsVar(1) = 4;
	return 0;
}

int ref_atm_lay_int(char *fname) 
{
	int m1, n1, l1, m2, n2, l2, m3, n3, l3, m4, n4, l4, m5, n5, l5, m6, n6, l6, 
		m7, n7, minlhs=1, maxlhs=1, minrhs=1, maxrhs=6,	
		x_axes=256, y_axes=256,  seed=0;
	double *outptr, *plaw, *inscl, *PS_data, *SubM_data, pixscale;
	PowerSpectrum PS; SubharmonicMethod SubM;
	RefractiveLayer Rlayer;	

	CheckRhs(minrhs,maxrhs);
	CheckLhs(minlhs,maxlhs);

	plaw = power_law(-11/3.0,0.15,0.5e-6);
	inscl = null_inner_scale();
    
	PS_data = create_power_spectrum( plaw, inscl );
	SubM_data = null_subharmonic_method();

	GetRhsVar(1, "d", &m1, &n1, &l1);
	if ( *stk(l1)<0 )  {
		sciprint("Error: arguments pixscale, must be positive\r\n");
		return 0;
	}
	pixscale = (*stk(l1));

	if (Rhs>1) {
		GetRhsVar(2, "i", &m2, &n2, &l2);
		if ( *istk(l2)<0  )  {
			sciprint("Error: the dimenstion of x_axes must be positive\r\n");
			return 0;
		}
		x_axes=(*istk(l2)); 
		y_axes=x_axes;  
	}

	if (Rhs>2) {
		GetRhsVar(3, "i", &m3, &n3, &l3);
		if ( *istk(l3)<0  )  {
			sciprint("Error: arguments random seed must be positive\r\n");
			return 0;
		}
		seed=(*istk(l3)); 
	}

	if (Rhs>3) {
		GetRhsVar(4, "i", &m4, &n4, &l4);
		if ( *istk(l4) < 0  )  {
			sciprint("Error: the dimenstion of y_axes  must be positive\r\n");
			return 0;
		}
		y_axes=(*istk(l4)); 
	}

	if (Rhs>4) {
		GetRhsVar(5, "d", &m5, &n5, &l5);
		if ( (m5!=1)||(n5 != 7 ) )  {
			sciprint("Error: fifth arguments must be power spectrum\r\n");
			return 0;
		}
		PS_data = stk(l5);
	}
	PS=array2PowerSpectrum( PS_data );

	if (Rhs>5) {
		GetRhsVar(6, "d", &m6, &n6, &l6);
		if ( (m6!=1)||(n6 != 4 ) )  {
			sciprint("Error: sixth arguments must be sunhormonic method\r\n");
			return 0;
		}
		SubM_data = stk(l6);
	}
	SubM=array2SubharmonicMethod( SubM_data );

	Rlayer=construct_RefractiveLayer( seed, PS, SubM, x_axes, y_axes, pixscale);
	outptr = RefractiveLayer2array( Rlayer );

	m7=1; n7=x_axes*y_axes+19;
	CreateVarFromPtr(7,"d", &m7, &n7, &outptr);   
	LhsVar(1) = 7;
	return 0;
}

int set_lay_wind_int(char *fname) 
{
	int  m1, n1, l1, m2, n2, l2, minlhs=1,
		 maxlhs=1, minrhs=2, maxrhs=2;
	RefractiveLayer layer; ThreeVector tv;
	int x_axes, y_axes;	double *outptr;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;
	
	GetRhsVar(1, "d", &m1, &n1, &l1);  
	x_axes=(int)(*stk(l1)); y_axes=(int)(*stk(l1+1));
	if ( (m1!=1)||(n1 != (x_axes*y_axes+19) ) )  {
		sciprint("Error: first arguments must be refractive atmospheric layer \r\n");
		return 0;
	}
	layer=array2RefractiveLayer( stk(l1) );

	GetRhsVar(2, "d", &m2, &n2, &l2); 
	if ( (m2 !=1)||(n2 != 3) )  {
		sciprint("Error: second arguments must be 3D wind vector\r\n");
		return 0;
	}
	tv=construct_ThreeVector( *stk(l2),*stk(l2+1),*stk(l2+2) );

	layer.wind_TV=tv;
	outptr = RefractiveLayer2array(layer);
	CreateVarFromPtr(3,"d", &m1, &n1,&outptr);   
	LhsVar(1) = 3;
	return 0;
}


int set_lay_frame_int(char *fname) 
{
	int  m1, n1, l1, m2, n2, l2, minlhs=1,
		 maxlhs=1, minrhs=2, maxrhs=2;
	RefractiveLayer layer; ThreeFrame tf;
	int x_axes, y_axes;	double *outptr;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;
	
	GetRhsVar(1, "d", &m1, &n1, &l1);  
	x_axes=(int)(*stk(l1)); y_axes=(int)(*stk(l1+1));
	if ( (m1!=1)||(n1 != (x_axes*y_axes+19) ) )  {
		sciprint("Error: first arguments must be refractive atmospheric layer \r\n");
		return 0;
	}
	layer=array2RefractiveLayer( stk(l1) );

	GetRhsVar(2, "d", &m2, &n2, &l2); 
	if ( (m2 !=1)||(n2 != 12) )  {
		sciprint("Error: second arguments must be 3D frame\r\n");
		return 0;
	}
	tf=array2ThreeFrame( stk(l2) );

	layer.Op.TF=tf;
	outptr = RefractiveLayer2array(layer);
	CreateVarFromPtr(3,"d", &m1, &n1,&outptr);   
	LhsVar(1) = 3;
	return 0;
}

int lay_pixarr_int(char *fname) 
{
	int  m1, n1, l1, m2, n2, minlhs=1,
		maxlhs=1, minrhs=1, maxrhs=1;
	RefractiveLayer layer;
	int x_axes, y_axes;
	double *outptr;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;
	
	GetRhsVar(1, "d", &m1, &n1, &l1);  
	x_axes=(int)(*stk(l1)); y_axes=(int)(*stk(l1+1));
	if ( (m1!=1)||(n1 != (x_axes*y_axes+19) ) )  {
		sciprint("Error: first arguments must be refractive atmospheric layer \r\n");
		return 0;
	}

	layer=array2RefractiveLayer( stk(l1) );
	outptr = PixelArray2array(layer.PixArr);

	m2=1; n2=x_axes*y_axes+2;
	CreateVarFromPtr(2,"d", &m2, &n2,&outptr);   
	LhsVar(1) = 2;
	return 0;
}

int ref_atm_lay_fits_int(char *fname) 
{
	int  m1, n1, l1, m2, n2, l2, 
		minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
	RefractiveLayer layer;
	int x_axes, y_axes;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;
	
	GetRhsVar(1, "d", &m1, &n1, &l1);  
	x_axes=(int)(*stk(l1)); y_axes=(int)(*stk(l1+1));
	if ( (m1!=1)||(n1 != (x_axes*y_axes+19) ) )  {
		sciprint("Error: first arguments must be refractive atmospheric layer \r\n");
		return 0;
	}

	GetRhsVar(2, "c", &m2, &n2, &l2);  
	layer=array2RefractiveLayer( stk(l1) );
	write_refractive_layer_file( layer, cstk(l2) );
	LhsVar(1) = 1;
	return 0;
}

int ref_atm_lay_transform_int(char *fname) 
{
	FIELD field;
	RefractiveLayer Rlayer;
	WaveFront  WF_in, WF_out;

	int  m1, n1, l1, m2, n2, l2,
		minlhs=1, maxlhs=1, minrhs=2, maxrhs=2;
	int x_axes, y_axes, number, number2;
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
	x_axes=(int)(*stk(l2)); y_axes=(int)(*stk(l2+1));
	if ( (m2 !=1)||(n2 != (x_axes*y_axes+19) ) )  {
		sciprint("Error: second arguments must be refractive atmospheric layer \r\n");
		return 0;
	}
	Rlayer=array2RefractiveLayer( stk(l2) );
	WF_out =RefractiveLayer_transform( Rlayer, WF_in );
	WaveFront_FIELD( WF_out, &field );
	
	outptr = field2array( field );
	CreateVarFromPtr( 3, "d", &m1, &n1, &outptr ); 
	LhsVar(1) = 3;
	return 0;
}
