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

#ifndef WRAP_REFRACTIVE_LAYER_H
#define WRAP_REFRACTIVE_LAYER_H

#include "arroyo_wrap.h"

typedef struct{ 
	int power_law_type;
	int inner_scale_type;
	double outer_scale,inner_scale;
	double exponent,r_0_meters,r_0_ref_wavelength_meters;
} PowerSpectrum;

typedef struct{ 
	int subharmonic_method_type;
	int subharmonic_depth,subpixels_per_level,subpixels_per_pixel;
} SubharmonicMethod;

typedef struct{ 
	Optic Op;
	PixelArray PixArr;
	double pixel_scale;
	ThreeVector wind_TV;
} RefractiveLayer;

#ifdef __cplusplus
#include "../arroyo/arroyo.h"

using namespace Arroyo;

power_spectrum 
*c2cpp_PowerSpectrum( PowerSpectrum PS );

subharmonic_method 
*c2cpp_SubharmonicMethod( SubharmonicMethod SubM );

refractive_atmospheric_layer<double> 
c2cpp_RefractiveLayer( RefractiveLayer Rlayer );

RefractiveLayer 
cpp2c_RefractiveLayer( refractive_atmospheric_layer<double> ref_atm_layer );

refractive_atmospheric_layer<double> 
refractive_layer_file_constructor(char *fname);

extern "C" {
#endif

	PowerSpectrum create_PowerSpectrum(	int power_law_type,	double exponent,
		double r_0_meters,double r_0_ref_wavelength_meters, double outer_scale,
		int inner_scale_type, double inner_scale );

	double *PowerSpectrum2array ( PowerSpectrum PS );
	PowerSpectrum array2PowerSpectrum (double *inptr);

	SubharmonicMethod create_SubharmonicMethod(int subharmonic_method_type,
		int subharmonic_depth,int subpixels_per_level,int subpixels_per_pixel);

	double *SubharmonicMethod2array ( SubharmonicMethod SubM );
	SubharmonicMethod array2SubharmonicMethod (double *inptr);

	RefractiveLayer construct_RefractiveLayer( int seed, PowerSpectrum PS,
		SubharmonicMethod SubM,	int x_axes, int y_axes, double pixscale );
	RefractiveLayer PixArr_RefractiveLayer (PixelArray PixArr, double pixscale );

	double *RefractiveLayer2array ( RefractiveLayer Rlayer );
	RefractiveLayer array2RefractiveLayer (double *inptr);

	WaveFront RefractiveLayer_transform(RefractiveLayer Rlayer, WaveFront WF);
	void write_refractive_layer_file( RefractiveLayer layer, char *fname );
	
	void write_PixelArray_file( PixelArray PixArr,char *fname,double timestamp );

#ifdef __cplusplus
}
#endif

#endif 
