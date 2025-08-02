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

#ifndef WRAP_ATMOSPHERIC_MODEL_H
#define WRAP_ATMOSPHERIC_MODEL_H

#include "arroyo_wrap.h"

typedef struct{ 
	int type;
	int nlayers;
	ThreeFrame TF;
	double *layer_heights;
	PowerSpectrum *power_spectra;
	double r_0_meters, r_0_ref_wavelength_meters;
	double Hufnagel_Valley_pseudowind, Gemini_outer_scale_meters;
	int Gemini_ground_layer_quality,Gemini_focal_anisoplanatism_quality, 
		Gemini_extended_profile; 
} AtmosphericModel;

#ifdef __cplusplus 

#include "../arroyo/arroyo.h"
using namespace Arroyo;

refractive_atmospheric_model 
*c2cpp_AtmosphericModel( AtmosphericModel AtmModel );

extern "C" {
#endif

	AtmosphericModel construct_AtmosphericModel( int nlayers,
		PowerSpectrum *power_spectra,double *layer_heights, ThreeFrame TF );

	AtmosphericModel construct_Ellerbroek_Cerro_Pachon_model
		( ThreeFrame TF, double r_0_meters,double r_0_ref_wavelength_meters );

	AtmosphericModel construct_Ellerbroek_Mauna_Kea_model
		( ThreeFrame TF, double r_0_meters,double r_0_ref_wavelength_meters );

	AtmosphericModel construct_Palomar_DIMM_MASS_model
		( ThreeFrame TF, double r_0_meters,double r_0_ref_wavelength_meters );

	AtmosphericModel construct_Hufnagel_Valley_model( ThreeFrame TF,  
		int nlayers, double *layer_heights, double Hufnagel_Valley_pseudowind );

	AtmosphericModel construct_SLCSAT_day_model
		( ThreeFrame TF, int nlayers, double *layer_heights );

	AtmosphericModel construct_SLCSAT_night_model
		( ThreeFrame TF, int nlayers, double *layer_heights );

	AtmosphericModel construct_TMT_SRD_v13_Cn2_model( ThreeFrame TF );

	AtmosphericModel construct_Gemini_GLAO_study_model
		( ThreeFrame TF, int Gemini_ground_layer_quality, 
		int Gemini_focal_anisoplanatism_quality, 
		int Gemini_extended_profile,
		double Gemini_outer_scale_meters );

	int AtmosphericModel_layer_number( AtmosphericModel AtmModel );
	double *AtmosphericModel_layer_heights( AtmosphericModel AtmModel );

	WavefrontHeader get_AtmosphericModel_WavefrontHeader( 
	    AtmosphericModel AtmModel, double wavelength, double pixscale, 
		Emitter Emt, APERTURE ap, int layer_foreshortening, int pplan_type );

	RefractiveLayer *get_AtmosphericModel_RefractiveLayer( AtmosphericModel AtmModel, 
		SubharmonicMethod SubM, int num_dwfhs, WavefrontHeader *dwfhdrs,
		double *layer_pixscales, ThreeVector *layer_wind_vectors, double time_interval,
		int layer_axes_wind_vector_aligned, int layer_foreshortening );
	void write_atm_mod_ref_layer_file( RefractiveLayer layer, int n, char *fname );

	WaveFront AtmosphericModel_to_APERTURE_transform(
		int nlayers, RefractiveLayer *Rlayers, APERTURE AP,
		int diffractive_near_field_propagator,WaveFront WF );

	double *AtmosphericModel2array ( AtmosphericModel AtmModel );
	AtmosphericModel array2AtmosphericModel (double *inptr);

#ifdef __cplusplus
}
#endif

#endif
