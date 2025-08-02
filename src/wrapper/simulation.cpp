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

#include "simulation.h"

 AtmosphericModel SimuPara_AtmModel(SIMULATION_PARAMETER para) 
  {
	  AtmosphericModel AtmModel;

	  three_frame global_tf;
	  ThreeFrame TF=cpp2c_ThreeFrame(global_tf);
	  AtmModel.TF=TF;

	  AtmModel.type=para.atm_type;
	  AtmModel.nlayers=para.nlayers;
	  AtmModel.r_0_meters=para.r_0_meters;
	  AtmModel.r_0_ref_wavelength_meters
		  =para.r_0_ref_wavelength_meters;
	  AtmModel.Hufnagel_Valley_pseudowind
		  =para.Hufnagel_Valley_pseudowind;
	  AtmModel.Gemini_outer_scale_meters
		  =para.Gemini_outer_scale_meters;
	  AtmModel.Gemini_ground_layer_quality
		  =para.Gemini_ground_layer_quality;
	  AtmModel.Gemini_focal_anisoplanatism_quality
		  =para.Gemini_focal_anisoplanatism_quality;
	  AtmModel.Gemini_extended_profile
		  =para.Gemini_extended_profile; 
	  AtmModel.layer_heights=para.layer_heights;

	  return AtmModel;
}

double* SimuPara2array(SIMULATION_PARAMETER para)
{
	int i,total;
	double *tmp;

	total=54+para.nlayers+3*para.detected_wavelengths_number;
	tmp=(double*)calloc(total,sizeof(double));

	tmp[0]=para.nsteps_in_simulation;
	tmp[1]=para.nsteps_to_delay_closing_dm_loop;
	tmp[2]=para.interval_for_propagating_science_wavefronts;
	tmp[3]=para.nsteps_to_delay_propagating_science_wavefronts;
	tmp[4]=para.timestep_seconds;

	tmp[5]=para.seed;

	tmp[6]=para.tip_tilt_angular_velocity_rad_per_sec;
	tmp[7]=para.tip_tilt_proportional_gain;
	tmp[8]=para.tip_tilt_integral_gain;
	tmp[9]=para.dm_actuator_linear_velocity_meters_per_sec;
	tmp[10]=para.dm_proportional_gain;
	tmp[11]=para.dm_integral_gain;

	tmp[12]=para.sensing_wavelength_meters;
	tmp[13]=para.detected_wavelengths_number;

	tmp[14]=para.subharmonic_method_code;
	tmp[15]=para.subharmonic_depth;
	tmp[16]=para.generalized_subharmonic_subpixels_per_level;
	tmp[17]=para.generalized_subharmonic_subpixels_per_pixel;

	tmp[18]=para.tropospheric_height_meters;
	tmp[19]=para.tropospheric_thickness_meters;
	tmp[20]=para.rms_ground_wind_speed_meters_per_sec;
	tmp[21]=para.rms_tropospheric_wind_speed_meters_per_sec;

	tmp[22]=para.wavefront_pixel_scale_meters;
	tmp[23]=para.layer_pixel_scale_meters;

	tmp[24]=para.n_x_emitters;
	tmp[25]=para.n_y_emitters;
	tmp[26]=para.emitter_spacing_arcsecs;

	tmp[27]=para.verbose;
	tmp[28]=para.vverbose;

	tmp[29]=para.aperture_code;
	tmp[30]=para.inner_diameter_meters;
	tmp[31]=para.outer_diameter_meters;
	tmp[32]=para.nspiders;
	tmp[33]=para.spider_width;
	tmp[34]=para.edge_length;
	tmp[35]=para.gap_size;

	tmp[36]=para.lenslet_focal_length_meters;
	tmp[37]=para.lenslet_pitch_meters;
	tmp[38]=para.final_wf_pixels_per_lenslet;
	tmp[39]=para.final_wf_pixels_per_xform;

	tmp[40]=para.diffractive_near_field_propagator;
	tmp[41]=para.geometric_near_field_propagator;

	tmp[42]=para.reconstructor_nsubaps;
	tmp[43]=para.reconstructor_areal_threshold;
	tmp[44]=para.reconstructor_eigenvalue_threshold;

	tmp[45]=para.atm_type;
	tmp[46]=para.nlayers;
	tmp[47]=para.r_0_meters;
	tmp[48]=para.r_0_ref_wavelength_meters;
	tmp[49]=para.Hufnagel_Valley_pseudowind;
	tmp[50]=para.Gemini_outer_scale_meters;
	tmp[51]=para.Gemini_ground_layer_quality;
	tmp[52]=para.Gemini_focal_anisoplanatism_quality;
	tmp[53]=para.Gemini_extended_profile;

	for(i=0;i<para.nlayers;i++) 
		tmp[54+i]=para.layer_heights[i];

	for(i=0;i<para.detected_wavelengths_number;i++) {
		tmp[54+para.nlayers+i]=para.detected_wavelengths_meters[i];
		tmp[54+para.nlayers+para.detected_wavelengths_number+i]=
			para.focal_plane_image_sizes_arcsecs[i];
		tmp[54+para.nlayers+2*para.detected_wavelengths_number+i]=
			para.oversampling_factors[i];
	}

	return tmp;
}

SIMULATION_PARAMETER array2SimuPara(double* tmp)
{
	int i;
	SIMULATION_PARAMETER para;

	para.nsteps_in_simulation=(int)tmp[0];
	para.nsteps_to_delay_closing_dm_loop=(int)tmp[1];
	para.interval_for_propagating_science_wavefronts=(int)tmp[2];
	para.nsteps_to_delay_propagating_science_wavefronts=(int)tmp[3];
	para.timestep_seconds=tmp[4];

	para.seed=(int)tmp[5];

	para.tip_tilt_angular_velocity_rad_per_sec=tmp[6];
	para.tip_tilt_proportional_gain=tmp[7];
	para.tip_tilt_integral_gain=tmp[8];
	para.dm_actuator_linear_velocity_meters_per_sec=tmp[9];
	para.dm_proportional_gain=tmp[10];
	para.dm_integral_gain=tmp[11];

	para.sensing_wavelength_meters=tmp[12];
	para.detected_wavelengths_number=(int)tmp[13];

	para.subharmonic_method_code=(int)tmp[14];
	para.subharmonic_depth=(int)tmp[15];
	para.generalized_subharmonic_subpixels_per_level=(int)tmp[16];
	para.generalized_subharmonic_subpixels_per_pixel=(int)tmp[17];

	para.tropospheric_height_meters=tmp[18];
	para.tropospheric_thickness_meters=tmp[19];
	para.rms_ground_wind_speed_meters_per_sec=tmp[20];
	para.rms_tropospheric_wind_speed_meters_per_sec=tmp[21];

	para.wavefront_pixel_scale_meters=tmp[22];
	para.layer_pixel_scale_meters=tmp[23];

	para.n_x_emitters=(int)tmp[24];
	para.n_y_emitters=(int)tmp[25];
	para.emitter_spacing_arcsecs=tmp[26];

	para.verbose=(int)tmp[27];
	para.vverbose=(int)tmp[28];

	para.aperture_code=(int)tmp[29];
	para.inner_diameter_meters=tmp[30];
	para.outer_diameter_meters=tmp[31];
	para.nspiders=(int)tmp[32];
	para.spider_width=tmp[33];
	para.edge_length=tmp[34];
	para.gap_size=tmp[35];

	para.lenslet_focal_length_meters=tmp[36];
	para.lenslet_pitch_meters=tmp[37];
	para.final_wf_pixels_per_lenslet=(int)tmp[38];
	para.final_wf_pixels_per_xform=(int)tmp[39];

	para.diffractive_near_field_propagator=(int)tmp[40];
	para.geometric_near_field_propagator=(int)tmp[41];

	para.reconstructor_nsubaps=(int)tmp[42];
	para.reconstructor_areal_threshold=tmp[43];
	para.reconstructor_eigenvalue_threshold=tmp[44];

	para.atm_type=(int)tmp[45];
	para.nlayers=(int)tmp[46];
	para.r_0_meters=tmp[47];
	para.r_0_ref_wavelength_meters=tmp[48];
	para.Hufnagel_Valley_pseudowind=tmp[49];
	para.Gemini_outer_scale_meters=tmp[50];
	para.Gemini_ground_layer_quality=(int)tmp[51];
	para.Gemini_focal_anisoplanatism_quality=(int)tmp[52];
	para.Gemini_extended_profile=(int)tmp[53];

	para.layer_heights=(double*)malloc(para.nlayers*sizeof(double));
	for(i=0;i<para.nlayers;i++) 
		para.layer_heights[i]=tmp[54+i];

	para.detected_wavelengths_meters=(double*)malloc(para.detected_wavelengths_number*sizeof(double));
	para.focal_plane_image_sizes_arcsecs=(double*)malloc(para.detected_wavelengths_number*sizeof(double));
	para.oversampling_factors=(double*)malloc(para.detected_wavelengths_number*sizeof(double));
	for(i=0;i<para.detected_wavelengths_number;i++) {
		para.detected_wavelengths_meters[i]=tmp[54+para.nlayers+i];
		para.focal_plane_image_sizes_arcsecs[i]=
			tmp[54+para.nlayers+para.detected_wavelengths_number+i];
		para.oversampling_factors[i]=
			tmp[54+para.nlayers+2*para.detected_wavelengths_number+i];
	}

	return para;
}
  