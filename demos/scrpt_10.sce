
mode(-1);
//
// -------------------------------------------------------------------------
// scitao - Scilab/Scicos Adaptive Optics tooolbox
//
// Copyright (C) 2006  IAPCM , Beijing, China.  Written by
// Chen jingyuan.  For comments or questions about this software,
// please contact the author at jingyuan_chen@yahoo.com.cn.
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software Foundation, 
// Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// -------------------------------------------------------------------------
//

nsteps_in_simulation=10;
nsteps_to_delay_closing_dm_loop=5;
timestep_seconds=0.002;

tip_tilt_angular_velocity_rad_per_sec=100;
tip_tilt_proportional_gain=0;tip_tilt_integral_gain=5e-2;
dm_actuator_linear_velocity_meters_per_sec=2.67e-4;
dm_proportional_gain=0;dm_integral_gain=1e-7;

sensing_wavelength_meters=0.6e-6;

tropospheric_height_meters=10000;
tropospheric_thickness_meters=5000;
rms_ground_wind_speed_meters_per_sec=5;
rms_tropospheric_wind_speed_meters_per_sec=30;

subharmonic_method_code=2;subharmonic_depth=3;

wavefront_pixel_scale_meters=0.02;
layer_pixel_scale_meters=0.02;

inner_diameter_meters=0.5;
outer_diameter_meters=2;

lenslet_focal_length_meters=0.012246;
lenslet_pitch_meters=0.000252;
final_wf_pixels_per_lenslet=32;
final_wf_pixels_per_xform=32;

reconstructor_nsubaps=16;
reconstructor_areal_threshold=0.5;
reconstructor_eigenvalue_threshold=1e-6;

r_0_meters=0.05;
r_0_ref_wavelength_meters=0.5e-6;

detected_wavelengths_meters=2.2e-6;
focal_plane_image_sizes_arcsecs=8;
oversampling_factors=8;

/////////////////////////////////////////////////

ap=annul_ap(outer_diameter_meters,inner_diameter_meters);
tmp_ap=circle(2*outer_diameter_meters);
ap_tf=three_frame();propagation_vector=three_vector(0,0,1);

atm_mod=ell_cer_pac_mod(r_0_meters,r_0_ref_wavelength_meters);
nlayers=atm_mod_lay_number(atm_mod);hghts=atm_mod_lay_heights(atm_mod);

act_dim=reconstructor_nsubaps+1;
lenslet_axes=reconstructor_nsubaps;
pitch=outer_diameter_meters/reconstructor_nsubaps;

num_dwfhs=2;
sensing_emt=emitter(0,0,0,-1);detected_emt=sensing_emt;
act_arr=actuator_array(act_dim,act_dim,pitch,dm_actuator_linear_velocity_meters_per_sec);

tmp_dm=ideal_dm(tmp_ap,act_arr);
tmp_lnslt_arr=lnslt_array(lenslet_axes,lenslet_axes,lenslet_focal_length_meters,pitch,final_wf_pixels_per_lenslet,final_wf_pixels_per_xform);
znke=znk_mod(1);znke=set_znk_cos_coef(znke,1,1,1);znke=set_znk_sin_coef(znke,1,1,1);
zpz_reconstrucor=reconstructor(ap,tmp_dm,tmp_lnslt_arr,znke, reconstructor_areal_threshold, reconstructor_eigenvalue_threshold);

tt_residuals=znk_mod(1);tt_commands=znk_mod(1);
ttm=ideal_ttm(tmp_ap,tip_tilt_angular_velocity_rad_per_sec);
tmp_pi=pi_controller(tip_tilt_proportional_gain,tip_tilt_integral_gain);
ttm_controller=ttm_pi_controller(tmp_pi);

dm_residuals=pixel_array(act_dim);dm_commands=pixel_array(act_dim);
tmp_tp=three_point(0,0,0);tmp_tv_x=three_vector(-1,0,0);tmp_tv_y=three_vector(0,-1,0);tmp_tv_z=three_vector(0,0,-1);
tmp_tf=three_frame(tmp_tp,tmp_tv_x,tmp_tv_y,tmp_tv_z);dm_ap=set_ap_frame(tmp_ap,tmp_tf);
dm=ideal_dm(dm_ap,act_arr);
dm_pi=pi_controller(dm_proportional_gain,dm_integral_gain);
dm_controller=dm_pi_controller(dm,dm_pi);

rescaled_wf_pixel_scale_meters = wavefront_pixel_scale_meters*lenslet_pitch_meters/pitch;
lnslt_arr=lnslt_array(lenslet_axes,lenslet_axes,lenslet_focal_length_meters,lenslet_pitch_meters,final_wf_pixels_per_lenslet,final_wf_pixels_per_xform);

sensing_dwfh=get_atm_mod_dwfh( atm_mod, sensing_wavelength_meters, wavefront_pixel_scale_meters, sensing_emt , ap , 0, 0 );
detected_dwfh=get_atm_mod_dwfh( atm_mod, detected_wavelengths_meters, wavefront_pixel_scale_meters, detected_emt , ap , 0, 0 );
dwfhs=[sensing_dwfh,detected_dwfh];

subm=lane_subm(subharmonic_depth);
hwm=hardy_wind(rms_ground_wind_speed_meters_per_sec,rms_tropospheric_wind_speed_meters_per_sec,tropospheric_height_meters,tropospheric_thickness_meters);
layer_wind_velocities=get_hardy_wind_vectors(hwm,nlayers,hghts);
time_interval=nsteps_in_simulation*timestep_seconds;
layers=get_atm_mod_layers(atm_mod,subm,num_dwfhs,dwfhs,layer_pixel_scale_meters,hwm,time_interval,0,0);

//////////////////////////////////////////////////////

M_PI=3.1415926;
arcsecs_to_radians = M_PI/180./3600.;
rad_to_arcsec = 180*3600/M_PI;
propagation_distance_meters = 1e10;

final_pixel_scales_meters = detected_wavelengths_meters*propagation_distance_meters /detected_dwfh(1)/detected_dwfh(4)/oversampling_factors;
final_pixel_scale_arcsecs = rad_to_arcsec*final_pixel_scales_meters/propagation_distance_meters;
final_array_dimensions = ceil(focal_plane_image_sizes_arcsecs/final_pixel_scale_arcsecs);

init_detected_dwf=pl_wave(detected_dwfh,detected_emt);
init_sensing_dwf=pl_wave(sensing_dwfh,sensing_emt);

////////////////////////////////////////////////////////

for i=0:nsteps_in_simulation-1

	timestamp=i*timestep_seconds;

	dwf=set_dwf_timestamp(init_detected_dwf,timestamp);
	dwf=atm_mod_lays_transform(dwf,layers,ap,0);
	dwf=aperture_transform(dwf,ap);
//	dwf_fits(dwf,"uncorrected_pupil_wf_",timestamp);

	tmp_dwf=fresnel_grt(dwf,propagation_distance_meters,final_pixel_scales_meters,final_array_dimensions);
	dwf_fits(tmp_dwf,"uncorrected_psf",timestamp);

	tt_dwf=ideal_ttm_transform(dwf,ttm);
	dwf=set_dwf_direction(tt_dwf,propagation_vector);
	dm_dwf=ideal_dm_transform(dwf,dm);
//	dwf_fits(dm_dwf,"corrected_pupil_wf_",timestamp);

	dwf=set_dwf_direction(dm_dwf,propagation_vector);
	tmp_dwf=fresnel_grt(dwf,propagation_distance_meters,final_pixel_scales_meters,final_array_dimensions);
	dwf_fits(tmp_dwf,"corrected_psf",timestamp);

	///////////////////////////////////////////////////////

	dwf=set_dwf_timestamp(init_sensing_dwf,timestamp);
	dwf=atm_mod_lays_transform(dwf,layers,ap,0);
	dwf=aperture_transform(dwf,ap);

	tt_dwf=ideal_ttm_transform(dwf,ttm);
	dwf=set_dwf_direction(tt_dwf,propagation_vector);
	dm_dwf=ideal_dm_transform(dwf,dm);

	tmp_dwf=set_dwf_frame(dm_dwf,ap_tf);
	tmp_dwf=set_dwf_pixel_scale(tmp_dwf,rescaled_wf_pixel_scale_meters);

	dwf=lnslt_arr_transform(tmp_dwf,lnslt_arr);

	nclip=(final_wf_pixels_per_xform-final_wf_pixels_per_lenslet)/2;
	dwf=dwf_clip_array(dwf,nclip);
	shc=create_shcentroids(dwf,lnslt_arr);

	tt_residuals = zernike_residuals( zpz_reconstrucor, shc );
	dm_residuals=zonal_residuals( zpz_reconstrucor, shc );

	[ttm_controller,tt_commands]=ttm_pi_update(ttm_controller,tt_residuals,tt_commands);
	x_orient=-get_znk_cos_coef(tt_commands,1,1)*arcsecs_to_radians;
	y_orient=-get_znk_sin_coef(tt_commands,1,1)*arcsecs_to_radians;
	ttm_command_vector=three_vector(x_orient,y_orient,1);
	ttm=ttm_update(ttm,ttm_command_vector,timestamp);

	if i>=nsteps_to_delay_closing_dm_loop then
		[dm_controller,dm_commands]=dm_pi_update(dm_controller,dm_residuals,dm_commands);
		tmp_commands=multi_pixarr(dm_commands,-1);
		dm=dm_update(dm,tmp_commands,timestamp);
	end
		
end
