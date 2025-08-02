
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

You should have received a copy of the GNU General Public License aint
with this program; if not, write to the Free Software Foundation, 
Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
*/

#ifndef SIMULATION_PARAMETER_H
#define SIMULATION_PARAMETER_H

#include "arroyo_wrap.h"

typedef struct{ 

	int nsteps_in_simulation;
	int nsteps_to_delay_closing_dm_loop;
	int interval_for_propagating_science_wavefronts;
	int nsteps_to_delay_propagating_science_wavefronts;
	double timestep_seconds;

	int seed;

	double tip_tilt_angular_velocity_rad_per_sec;
	double tip_tilt_proportional_gain;
	double tip_tilt_integral_gain;
	double dm_actuator_linear_velocity_meters_per_sec;
	double dm_proportional_gain;
	double dm_integral_gain;

	double sensing_wavelength_meters;
	int detected_wavelengths_number;

	int subharmonic_method_code;
	int subharmonic_depth;
	int generalized_subharmonic_subpixels_per_level; 
    int generalized_subharmonic_subpixels_per_pixel;

	double tropospheric_height_meters;
	double tropospheric_thickness_meters;
	double rms_ground_wind_speed_meters_per_sec;
	double rms_tropospheric_wind_speed_meters_per_sec;

	double wavefront_pixel_scale_meters;
	double layer_pixel_scale_meters;

	int n_x_emitters;
	int n_y_emitters;
	double emitter_spacing_arcsecs;

	int verbose;
	int vverbose;

	int aperture_code;
	double inner_diameter_meters;
	double outer_diameter_meters;
	int nspiders; 
	double spider_width;
	double edge_length;
	double gap_size;

	double lenslet_focal_length_meters;
	double lenslet_pitch_meters;
	int final_wf_pixels_per_lenslet;
	int final_wf_pixels_per_xform;

	int diffractive_near_field_propagator;
	int geometric_near_field_propagator;

	int reconstructor_nsubaps;
	double reconstructor_areal_threshold;
	double reconstructor_eigenvalue_threshold;

	int atm_type;
	int nlayers;
	double r_0_meters;
	double r_0_ref_wavelength_meters;
	double Hufnagel_Valley_pseudowind; 
	double Gemini_outer_scale_meters;
	int Gemini_ground_layer_quality;
	int Gemini_focal_anisoplanatism_quality;
	int Gemini_extended_profile; 

	double *layer_heights;

	double* detected_wavelengths_meters;
	double* focal_plane_image_sizes_arcsecs;
	double* oversampling_factors;

} SIMULATION_PARAMETER;

#ifdef __cplusplus

#include "../arroyo/arroyo.h"
using namespace Arroyo;

extern "C" {
#endif

	double* SimuPara2array(SIMULATION_PARAMETER);
	SIMULATION_PARAMETER array2SimuPara(double*);

	AtmosphericModel SimuPara_AtmModel(SIMULATION_PARAMETER);

	void scao_aperture(SIMULATION_PARAMETER);
	void scao_reconstructor(SIMULATION_PARAMETER,char*);
	void scao_emitter(SIMULATION_PARAMETER);
	void scao_wavefront_header(SIMULATION_PARAMETER);
	void scao_ref_atm_model(SIMULATION_PARAMETER);
	void scao_ref_atm_layer(SIMULATION_PARAMETER);
	void scao_ttm_and_pi(SIMULATION_PARAMETER);
	void scao_dm_and_pi(SIMULATION_PARAMETER);
	void scao_lenslet_array(SIMULATION_PARAMETER);

	void scao_unaberrated_image(SIMULATION_PARAMETER);
	void scao_detected_aperture(SIMULATION_PARAMETER,int,int,int);
	void scao_uncorrected_detected_far_field(SIMULATION_PARAMETER,int,int,int);
	void scao_correction_detected(SIMULATION_PARAMETER,int,int,int);
	void scao_corrected_detected_far_field(SIMULATION_PARAMETER,int,int,int);
	void scao_sensing_aperture(SIMULATION_PARAMETER,int);
	void scao_correction_sensing(SIMULATION_PARAMETER,int);
	void scao_sensing_lenslet_array(SIMULATION_PARAMETER,int);
	void scao_construct_centroids(SIMULATION_PARAMETER,int);
	void scao_reconstruct_residuals(SIMULATION_PARAMETER,int);
	void scao_update_corrector_mirrors(SIMULATION_PARAMETER,int);

	void scao_write_wavefront(int,char*,double,int);

#ifdef __cplusplus
}
#endif

#endif
