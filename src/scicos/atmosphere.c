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

#include "optics_scicos.h"

void atmosphere(scicos_block *block,int flag) 
{
	char str[50];

	int i,*ipar,job = 1,mtype, power_law_type,inner_scale_type,
		Gemini_ground_layer_quality, Gemini_focal_anisoplanatism_quality, 
		Gemini_extended_profile, subharmonic_method_type, subharmonic_depth,
		subpixels_per_level, subpixels_per_pixel,layer_foreshortening,
		propagation_plan_type,layer_axes_wind_vector_aligned;

	double *rpar,*layer_heights,*layer_pixscales,tmp_pixsl,
		Gemini_outer_scale_meters,Hufnagel_Valley_pseudowind,outer_scale,
		inner_scale,exponent,r_0_meters, r_0_ref_wavelength_meters,
		tropopause_height,tropopause_thickness,ground_layer_wind_velocity,
		tropopause_wind_velocity,total_simu_time_interval;

	FIELD field;

	ThreeFrame TF;
	WavefrontHeader WfH[2];
	WaveFront WF_in, WF_out;

	HardyWind HWM;
	PowerSpectrum *PS,tmp_ps;
	SubharmonicMethod SubM;
 
	AtmosphericModel AtmModel;
	ThreeVector *layer_wind_vectors;

	rpar=block->rpar; 
	--rpar;
	ipar=block->ipar; 
	--ipar;

	if (flag==1) {

		propagation_plan_type=ipar[14];

		field=block_inptr_field(block,0);

		if(field.number) {
			WF_in = FIELD_WaveFront(field);
			WF_out=AtmosphericModel_to_APERTURE_transform(scao.num_Rlayers,
				scao.Rlayers,scao.global_AP,propagation_plan_type,WF_in);
			WaveFront_FIELD( WF_out,&field );
		}
		
		field_block_outptr(field,block,0);

		if(block->nin==2){

			field=block_inptr_field(block,1);

			if(field.number) {
				WF_in = FIELD_WaveFront(field);
				WF_out=AtmosphericModel_to_APERTURE_transform(scao.num_Rlayers,
					scao.Rlayers,scao.global_AP,propagation_plan_type,WF_in);
				WaveFront_FIELD( WF_out,&field );
			}
			
			field_block_outptr(field,block,1);

		}

	}

	else if (flag==5)  {
		scao.atm_mod_ref_lay_created=0;
	}

	else if (flag==6 && !scao.atm_mod_ref_lay_created) {

		if(( block->nin==2 && scao.detected_emt_created && scao.sensing_emt_created )
			|| 	( block->nin==1 && scao.detected_emt_created ) ){
			
			mtype = ipar[1]; 
			scao.num_Rlayers = ipar[2];

			Gemini_ground_layer_quality = ipar[3];
			Gemini_focal_anisoplanatism_quality = ipar[4];
			Gemini_extended_profile = ipar[5]; 
			
			power_law_type=ipar[6];
			inner_scale_type=ipar[7];
	
			subharmonic_method_type=ipar[8];
			subharmonic_depth=ipar[9];
			subpixels_per_level=ipar[10];
			subpixels_per_pixel=ipar[11];

			layer_foreshortening=ipar[12];
			layer_axes_wind_vector_aligned=ipar[13];

			propagation_plan_type=ipar[14];

			F2C(cvstr)(&(ipar[15]),&(ipar[16]), str,&job,strlen(str));
			str[ipar[15]] = '\0';

			tmp_pixsl=rpar[1];
			Hufnagel_Valley_pseudowind=rpar[2];
			
			r_0_meters=rpar[3];
			r_0_ref_wavelength_meters=rpar[4];
			
			Gemini_outer_scale_meters=rpar[5];
			
			exponent=rpar[6];
			inner_scale=rpar[7];
			outer_scale=rpar[8];

			ground_layer_wind_velocity=rpar[9];
			tropopause_wind_velocity=rpar[10];
			tropopause_height=rpar[11];
			tropopause_thickness=rpar[12];

			total_simu_time_interval=rpar[13];

            tmp_ps = create_PowerSpectrum ( 
				power_law_type,exponent,
				r_0_meters,	r_0_ref_wavelength_meters,
				outer_scale,inner_scale_type,inner_scale);
			
			layer_heights=(double*)malloc(scao.num_Rlayers*sizeof(double));
			PS=(PowerSpectrum*)malloc(scao.num_Rlayers*sizeof(PowerSpectrum));

			TF = default_ThreeFrame( );

			for(i=0;i<scao.num_Rlayers;i++) {
				PS[i]=tmp_ps;
				layer_heights[i]=rpar[14+i];
			}
			
			switch(mtype)
			{
			case 0:
				AtmModel=construct_AtmosphericModel
					(scao.num_Rlayers,PS,layer_heights,TF);
				break;
			case 1:
				AtmModel=construct_Ellerbroek_Cerro_Pachon_model
					(TF,r_0_meters,r_0_ref_wavelength_meters) ;
				break;		
			case 2:
				AtmModel=construct_Ellerbroek_Mauna_Kea_model
					(TF,r_0_meters,r_0_ref_wavelength_meters) ;
				break;		
			case 3:
				AtmModel=construct_Palomar_DIMM_MASS_model
					(TF,r_0_meters,r_0_ref_wavelength_meters) ;
				break;		
			case 4:
				AtmModel=construct_Hufnagel_Valley_model
					(TF,scao.num_Rlayers,layer_heights,Hufnagel_Valley_pseudowind) ;
				break;		
			case 5:
				AtmModel=construct_SLCSAT_day_model
					(TF,scao.num_Rlayers,layer_heights) ;
				break;		
			case 6:
				AtmModel=construct_SLCSAT_night_model
					(TF,scao.num_Rlayers,layer_heights) ;
				break;		
			case 7:
				AtmModel=construct_TMT_SRD_v13_Cn2_model(TF) ;
				break;		
			case 8:
				AtmModel=construct_Gemini_GLAO_study_model
					( TF, Gemini_ground_layer_quality, 
					Gemini_focal_anisoplanatism_quality,
					Gemini_extended_profile,
					Gemini_outer_scale_meters );
				break;		
			}

			scao.detected_WfH = get_AtmosphericModel_WavefrontHeader
				( AtmModel, scao.detected_wavelengths, 
				scao.detected_pixel_scale, scao.detected_emt,scao.global_AP,
				layer_foreshortening,propagation_plan_type );
			WfH[0]=scao.detected_WfH;

			if(scao.max_number<scao.detected_WfH.axes_x)
				scao.max_number=scao.detected_WfH.axes_x;
			if(scao.max_number<scao.detected_WfH.axes_y) 
				scao.max_number=scao.detected_WfH.axes_y;

			if(block->nin==2) {
				scao.sensing_WfH = get_AtmosphericModel_WavefrontHeader
					( AtmModel, scao.sensing_wavelength, 
					scao.sensing_pixel_scale,scao.sensing_emt,scao.global_AP,
					layer_foreshortening,propagation_plan_type );
				WfH[1]=scao.sensing_WfH;

				if(scao.max_number<scao.sensing_WfH.axes_x) 
					scao.max_number=scao.sensing_WfH.axes_x;
				if(scao.max_number<scao.sensing_WfH.axes_y) 
					scao.max_number=scao.sensing_WfH.axes_y;
			}

			HWM = construct_HardyWind ( ground_layer_wind_velocity,
				tropopause_wind_velocity,tropopause_height,
				tropopause_thickness );

			SubM = create_SubharmonicMethod ( subharmonic_method_type,
				subharmonic_depth, subpixels_per_level, subpixels_per_pixel);

			scao.num_Rlayers=AtmosphericModel_layer_number(AtmModel);
			layer_heights=AtmosphericModel_layer_heights(AtmModel);

			layer_pixscales=(double*)malloc(scao.num_Rlayers*sizeof(double));
			for(i=0;i<scao.num_Rlayers;i++) layer_pixscales[i]=tmp_pixsl;

			layer_wind_vectors = get_HardyWind_velocities
				( HWM,scao.num_Rlayers,layer_heights,TF );

			if(block->nin==2)
				scao.Rlayers=get_AtmosphericModel_RefractiveLayer 
				( AtmModel, SubM,  2, WfH, layer_pixscales,
				layer_wind_vectors,	total_simu_time_interval,
				layer_axes_wind_vector_aligned, layer_foreshortening );
			else
				scao.Rlayers=get_AtmosphericModel_RefractiveLayer 
				( AtmModel, SubM,  1, WfH, layer_pixscales,
				layer_wind_vectors,	total_simu_time_interval,
				layer_axes_wind_vector_aligned, layer_foreshortening );
			
			if( (strstr(str, "void")) == NULL )
				for(i=0;i<scao.num_Rlayers;i++) 
					write_atm_mod_ref_layer_file(scao.Rlayers[i],i,str);

			scao.atm_mod_ref_lay_created=1;

		}
	}

}
