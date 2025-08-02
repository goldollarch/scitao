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

#include <ctime>
#include "arroyo_wrap.h"

using namespace Arroyo;

PowerSpectrum create_PowerSpectrum(	int power_law_type,	double exponent,
		double r_0_meters,double r_0_ref_wavelength_meters, double outer_scale,
		int inner_scale_type, double inner_scale )
{
	PowerSpectrum PS;

	PS.power_law_type=power_law_type;
	PS.inner_scale_type=inner_scale_type;
	PS.exponent=exponent;
	PS.inner_scale=inner_scale;
	PS.outer_scale=outer_scale;
	PS.r_0_meters=r_0_meters;
	PS.r_0_ref_wavelength_meters=r_0_ref_wavelength_meters;

	return PS;
}

double *PowerSpectrum2array ( PowerSpectrum PS )
{
	double *tmp;
	tmp=(double*) calloc( 7, sizeof(double) );

	tmp[0]=PS.power_law_type;
	tmp[1]=PS.exponent;
	tmp[2]=PS.r_0_meters;
	tmp[3]=PS.r_0_ref_wavelength_meters;
	tmp[4]=PS.outer_scale;
	tmp[5]=PS.inner_scale_type;
	tmp[6]=PS.inner_scale;

	return tmp;
}

PowerSpectrum array2PowerSpectrum (double *inptr)
{
	PowerSpectrum PS = create_PowerSpectrum	( (int)inptr[0],
		inptr[1],inptr[2],inptr[3],inptr[4],(int)inptr[5],inptr[6]);
	return PS;
}

power_spectrum *c2cpp_PowerSpectrum( PowerSpectrum PS )
{
    power_spectrum *pspectrum=NULL;

	if (PS.power_law_type==1) 
	{
		von_karman_power_law vkplaw( PS.exponent, 
			PS.r_0_meters, PS.r_0_ref_wavelength_meters, PS.outer_scale );

		if(PS.inner_scale_type==1) {
			exponential_inner_scale expinscle( PS.inner_scale );
			pspectrum=new isotropic_power_law_spectrum
			<von_karman_power_law, exponential_inner_scale> (vkplaw,expinscle);
		}
		else if (PS.inner_scale_type==2) {
			frehlich_inner_scale frinscle( PS.inner_scale );
			pspectrum=new isotropic_power_law_spectrum
			<von_karman_power_law, frehlich_inner_scale> (vkplaw,frinscle);
		}
		else 
			pspectrum=new isotropic_power_law_spectrum
			<von_karman_power_law, null_inner_scale> (vkplaw,null_inner_scale());
	}

	else if (PS.power_law_type==2) 
	{
		greenwood_power_law gwplaw( PS.exponent, 
			PS.r_0_meters, PS.r_0_ref_wavelength_meters, PS.outer_scale );

		if(PS.inner_scale_type==1) {
			exponential_inner_scale expinscle( PS.inner_scale );
			pspectrum=new isotropic_power_law_spectrum
			<greenwood_power_law, exponential_inner_scale> (gwplaw,expinscle);
		}
		else if (PS.inner_scale_type==2) {
			frehlich_inner_scale frinscle( PS.inner_scale );
			pspectrum=new isotropic_power_law_spectrum
			<greenwood_power_law, frehlich_inner_scale>(gwplaw,frinscle);
		}
		else
			pspectrum=new isotropic_power_law_spectrum
			<greenwood_power_law, null_inner_scale>(gwplaw,null_inner_scale());
	}

	else {
		power_law plaw = power_law( PS.exponent, 
			PS.r_0_meters, PS.r_0_ref_wavelength_meters );

		if(PS.inner_scale_type==1) {
			exponential_inner_scale expinscle( PS.inner_scale );
			pspectrum=new isotropic_power_law_spectrum
			<power_law, exponential_inner_scale> (plaw,expinscle);
		}
		else if (PS.inner_scale_type==2) {
			frehlich_inner_scale frinscle( PS.inner_scale );
			pspectrum=new isotropic_power_law_spectrum
			<power_law, frehlich_inner_scale>(plaw,frinscle);
		}
		else 
			pspectrum=new isotropic_power_law_spectrum
			<power_law, null_inner_scale> (plaw,null_inner_scale());
	}

	return(pspectrum);
}

SubharmonicMethod create_SubharmonicMethod(int subharmonic_method_type,
	int subharmonic_depth,int subpixels_per_level,int subpixels_per_pixel)
{
	SubharmonicMethod SubM;

	SubM.subharmonic_method_type=subharmonic_method_type;
	SubM.subharmonic_depth=subharmonic_depth;
	SubM.subpixels_per_level=subpixels_per_level;
	SubM.subpixels_per_pixel=subpixels_per_pixel;

	return SubM;
}

double *SubharmonicMethod2array ( SubharmonicMethod SubM )
{
	double *tmp;
	tmp=(double*) calloc(4, sizeof(double) );
	tmp[0]=SubM.subharmonic_method_type;
	tmp[1]=SubM.subharmonic_depth;
	tmp[2]=SubM.subpixels_per_level;
	tmp[3]=SubM.subpixels_per_pixel;
	return tmp;
}

SubharmonicMethod array2SubharmonicMethod (double *inptr)
{
	SubharmonicMethod SubM = create_SubharmonicMethod	
		( (int)inptr[0], (int)inptr[1], (int)inptr[2], (int)inptr[3] );
	return SubM;
}

subharmonic_method *c2cpp_SubharmonicMethod( SubharmonicMethod SubM)
{
	
	subharmonic_method * subm = NULL;
	switch(SubM.subharmonic_method_type) {
		case 1:
			subm = new 
				quad_pixel_subharmonic_method(SubM.subharmonic_depth);
			break;
		case 2:
			subm = new
				Lane_subharmonic_method(SubM.subharmonic_depth);
			break;
		case 3:
			subm = new 
				Johansson_Gavel_subharmonic_method(SubM.subharmonic_depth);
			break;
		case 4:
			subm = new generalized_subharmonic_method(
				SubM.subharmonic_depth,
				SubM.subpixels_per_level,
				SubM.subpixels_per_pixel);
			break;
		default:
			subm = new null_subharmonic_method();
			break;
	}

	return(subm);

}

RefractiveLayer PixArr_RefractiveLayer (PixelArray PixArr, double pixscale )
{
	pixel_array<double> pixarr=c2cpp_PixelArray( PixArr );
	refractive_atmospheric_layer<double> ref_layer( pixarr, pixscale );
	return cpp2c_RefractiveLayer(ref_layer);
}

double *RefractiveLayer2array ( RefractiveLayer Rlayer )
{
	double *tmp,*outptr;
	int i, n_pixarr;

	n_pixarr=(Rlayer.PixArr.x_axes)*(Rlayer.PixArr.y_axes);

	outptr=(double*) calloc( n_pixarr+19, sizeof(double) );

	outptr[0]=Rlayer.PixArr.x_axes;outptr[1]=Rlayer.PixArr.y_axes;
	outptr[2]=Rlayer.pixel_scale;
	outptr[3]=Rlayer.wind_TV.x;outptr[4]=Rlayer.wind_TV.y;outptr[5]=Rlayer.wind_TV.z;

	outptr[6]=Rlayer.Op.foreshortening;
	tmp=ThreeFrame2array(Rlayer.Op.TF);
	for(i=0;i<12;i++) outptr[i+7]=tmp[i];

	for(i=0;i<n_pixarr;i++) 
		outptr[i+19]=Rlayer.PixArr.pixeldata[i];

	return outptr;
}

RefractiveLayer array2RefractiveLayer (double *inptr)
{
	PixelArray PixArr;
	int x_axes, y_axes, forsh;
	ThreeVector wind_TV;
	RefractiveLayer Rlayer; 	

	x_axes=(int)inptr[0]; y_axes=(int)inptr[1];
	Rlayer.pixel_scale=inptr[2];
	wind_TV=construct_ThreeVector( inptr[3], inptr[4], inptr[5] );

	forsh=(int)inptr[6];
	ThreeFrame TF=array2ThreeFrame(inptr+7);
	Rlayer.Op = create_Optic(forsh,TF);

	PixArr=create_PixelArray( x_axes, y_axes, 0 );
	set_PixelArray_data( &PixArr, inptr+19 );

	Rlayer.PixArr=PixArr;
	Rlayer.wind_TV=wind_TV;

	return Rlayer;
}

RefractiveLayer construct_RefractiveLayer( int seed, PowerSpectrum PS,
	SubharmonicMethod SubM, int x_axes, int y_axes, double pixscale)
{
	if(seed==0)
		seed = ( int) time(NULL);

# if (defined( __unix__ ) || (defined(linux)) )
	srandom(seed);
# else
	srand(seed);
# endif

	vector<long> layer_axes(2);
	layer_axes[0] = x_axes;  layer_axes[1] = y_axes;

	power_spectrum *ipls = c2cpp_PowerSpectrum( PS );
	
	subharmonic_method * subm = c2cpp_SubharmonicMethod( SubM );

	refractive_atmospheric_layer<double> 
		ref_layer( ipls, *subm, layer_axes,	pixscale );

	RefractiveLayer Rlayer=cpp2c_RefractiveLayer( ref_layer );
	return Rlayer;
}

refractive_atmospheric_layer<double> 
c2cpp_RefractiveLayer( RefractiveLayer Rlayer )
{
	pixel_array<double> pixarr=c2cpp_PixelArray( Rlayer.PixArr );
	refractive_atmospheric_layer<double> ref_layer( pixarr, Rlayer.pixel_scale );

	three_vector wind_vector = c2cpp_ThreeVector( Rlayer.wind_TV );
	ref_layer.set_wind_vector( wind_vector );

	three_frame tf=c2cpp_ThreeFrame(Rlayer.Op.TF);
	ref_layer.three_frame::operator=(tf);
	ref_layer.set_foreshortening(Rlayer.Op.foreshortening);

	return ref_layer;
}

RefractiveLayer 
cpp2c_RefractiveLayer( refractive_atmospheric_layer<double> ref_atm_layer )
{
	int forsh;
	RefractiveLayer Rlayer;

	forsh=ref_atm_layer.get_foreshortening();
	ThreeFrame TF=cpp2c_ThreeFrame(ref_atm_layer);
	Rlayer.Op = create_Optic(forsh,TF);

	Rlayer.PixArr = cpp2c_PixelArray( ref_atm_layer );
	Rlayer.pixel_scale = ref_atm_layer.get_pixel_scale();
	Rlayer.wind_TV = cpp2c_ThreeVector(ref_atm_layer.get_wind_vector());

	return Rlayer;
}

refractive_atmospheric_layer<double> 
refractive_layer_file_constructor(char *fname)
{
	stringstream filename_stream;
	filename_stream.str(""); filename_stream << fname << ".fits";
	refractive_atmospheric_layer<double> 
		ref_atm_layer(filename_stream.str().c_str());
	//  it can not run under debug mode. I don't know why! 

	return ref_atm_layer;
}

WaveFront RefractiveLayer_transform( RefractiveLayer Rlayer, WaveFront WF )
{
	refractive_atmospheric_layer<double>  
		ref_layer = c2cpp_RefractiveLayer( Rlayer );
	diffractive_wavefront<double> dwf=c2cpp_WaveFront( WF );

	ref_layer.transform(dwf);

	WaveFront tmp_WF=cpp2c_WaveFront(dwf);
	return tmp_WF;
}

void write_refractive_layer_file( RefractiveLayer layer, char *fname )
{
	refractive_atmospheric_layer<double> 
		ref_layer = c2cpp_RefractiveLayer(layer);

	ref_layer.write(simple_filename(fname,".fits"));
}

void write_PixelArray_file( PixelArray PixArr, char *fname, double timestamp )
{
	RefractiveLayer layer=PixArr_RefractiveLayer(PixArr,1);

	refractive_atmospheric_layer<double> 
		ref_layer = c2cpp_RefractiveLayer(layer);

    if( timestamp <0 )
		ref_layer.write(simple_filename(fname,".fits"));
	else
		ref_layer.write(Current_filename(fname,timestamp,".fits"));
}
