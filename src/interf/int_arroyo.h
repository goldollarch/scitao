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

#ifndef INT_ARROYO_H
#define INT_ARROYO_H

extern Gatefunc dwf_fits_int;
extern Gatefunc circular_aperture_int;
extern Gatefunc rectangular_aperture_int;
extern Gatefunc annular_aperture_int;
extern Gatefunc hexagonal_aperture_int;
extern Gatefunc spidered_annular_aperture_int;
extern Gatefunc tiled_hexagonal_aperture_int;
extern Gatefunc plane_wave_emit_int;
extern Gatefunc spherical_wave_emit_int;
extern Gatefunc geom_prop_int;
extern Gatefunc near_ang_int;
extern Gatefunc near_fresnel_int;
extern Gatefunc far_fresnel_int;
extern Gatefunc far_fraunhoffer_int;
extern Gatefunc far_fresnel_GRT_int;
extern Gatefunc far_fraunhoffer_GRT_int;
extern Gatefunc lenslet_array_int;
extern Gatefunc lenslet_array_transform_int;
extern Gatefunc create_shcentroids_int;
extern Gatefunc shc_fits_int;
extern Gatefunc aperture_transform_int;
extern Gatefunc create_actuator_arr_int;
extern Gatefunc create_ideal_dm_int;
extern Gatefunc ideal_dm_transform_int;
extern Gatefunc set_dm_timestamp_int;
extern Gatefunc set_dm_actuator_positions_int;
extern Gatefunc set_dm_actuator_commands_int;
extern Gatefunc three_point_int;
extern Gatefunc three_vector_int;
extern Gatefunc three_frame_int;
extern Gatefunc three_reflection_int;
extern Gatefunc diffractive_wavefront_header_int;
extern Gatefunc set_dwf_timestamp_int;
extern Gatefunc set_dwf_frame_int;
extern Gatefunc create_ideal_ttm_int;
extern Gatefunc set_ttm_timestamp_int;
extern Gatefunc set_ttm_commands_vector_int;
extern Gatefunc ideal_ttm_transform_int;
extern Gatefunc create_zernike_int;
extern Gatefunc set_cos_coef_int;
extern Gatefunc set_sin_coef_int;
extern Gatefunc zernike_fits_int;
extern Gatefunc PI_controller_int;
extern Gatefunc ttm_PI_controller_int;
extern Gatefunc dm_PI_controller_int;
extern Gatefunc update_ttm_PI_controller_int;
extern Gatefunc update_dm_PI_controller_int;
extern Gatefunc TipTiltMirror_update_int;
extern Gatefunc DeformableMirror_update_int;
extern Gatefunc power_law_int;
extern Gatefunc von_karman_power_law_int;
extern Gatefunc greenwood_power_law_int;
extern Gatefunc null_inner_scale_int;
extern Gatefunc exponential_inner_scale_int;
extern Gatefunc frehlich_inner_scale_int;
extern Gatefunc create_power_spectrum_int;
extern Gatefunc null_subharmonic_method_int;
extern Gatefunc quad_pixel_subharmonic_method_int;
extern Gatefunc Lane_subharmonic_method_int;
extern Gatefunc Johansson_Gavel_subharmonic_method_int;
extern Gatefunc generalized_subharmonic_method_int;
extern Gatefunc create_pixel_array_int;
extern Gatefunc set_pixel_array_data_int;
extern Gatefunc ref_atm_lay_int;
extern Gatefunc ref_atm_lay_fits_int;
extern Gatefunc ref_atm_lay_transform_int;
extern Gatefunc Hardy_wind_model_int;
extern Gatefunc get_HardyWind_velocities_int;
extern Gatefunc AtmosphericModel_int;
extern Gatefunc Ellerbroek_Cerro_Pachon_model_int;
extern Gatefunc Ellerbroek_Mauna_Kea_model_int;
extern Gatefunc Palomar_DIMM_MASS_model_int;
extern Gatefunc Hufnagel_Valley_model_int;
extern Gatefunc SLCSAT_day_model_int;
extern Gatefunc SLCSAT_night_model_int;
extern Gatefunc TMT_SRD_v13_Cn2_model_int;
extern Gatefunc Gemini_GLAO_study_model_int;
extern Gatefunc construct_Emitter_int;
extern Gatefunc get_AtmosphericModel_WavefrontHeader_int;
extern Gatefunc get_AtmosphericModel_RefractiveLayer_int;
extern Gatefunc atm_mod_layers_fits_int;
extern Gatefunc create_RECONSTRUCTOR_int;
extern Gatefunc arroyo_reconstructor_fits_int;
extern Gatefunc arroyo_reconstruct_zernike_residuals_int;
extern Gatefunc arroyo_reconstruct_zonal_residuals_int;
extern Gatefunc three_translation_int;
extern Gatefunc three_rotation_int;
extern Gatefunc AtmosphericModel_layer_number_int;
extern Gatefunc AtmosphericModel_layer_heights_int;
extern Gatefunc AtmosphericModel_to_APERTURE_transform_int;

static GenericTable arroyo_Tab[]={
  {(Myinterfun)sci_gateway,dwf_fits_int,"dwf_fits"},
  {(Myinterfun)sci_gateway,circular_aperture_int,"circle"},
  {(Myinterfun)sci_gateway,rectangular_aperture_int,"rectang"},
  {(Myinterfun)sci_gateway,annular_aperture_int,"annul_ap"},
  {(Myinterfun)sci_gateway,hexagonal_aperture_int,"hex_ap"},
  {(Myinterfun)sci_gateway,spidered_annular_aperture_int,"spid_annul_ap"},
  {(Myinterfun)sci_gateway,tiled_hexagonal_aperture_int,"til_hex_ap"},
  {(Myinterfun)sci_gateway,plane_wave_emit_int,"pl_wave"},
  {(Myinterfun)sci_gateway,spherical_wave_emit_int,"sp_wave"},
  {(Myinterfun)sci_gateway,geom_prop_int,"geom_propagation"},
  {(Myinterfun)sci_gateway,near_ang_int,"near_angular"},
  {(Myinterfun)sci_gateway,near_fresnel_int,"near_fresnel"},
  {(Myinterfun)sci_gateway,far_fresnel_int,"far_fresnel"},
  {(Myinterfun)sci_gateway,far_fraunhoffer_int,"far_fraunhoffer"},
  {(Myinterfun)sci_gateway,far_fresnel_GRT_int,"fresnel_grt"},
  {(Myinterfun)sci_gateway,far_fraunhoffer_GRT_int,"fraunhoffer_grt"},
  {(Myinterfun)sci_gateway,lenslet_array_int,"lnslt_array"},
  {(Myinterfun)sci_gateway,lenslet_array_transform_int,"lnslt_arr_transform"},
  {(Myinterfun)sci_gateway,create_shcentroids_int,"create_shcentroids"},
  {(Myinterfun)sci_gateway,shc_fits_int,"shc_fits"},
  {(Myinterfun)sci_gateway,aperture_transform_int,"aperture_transform"},
  {(Myinterfun)sci_gateway,create_actuator_arr_int,"actuator_array"},
  {(Myinterfun)sci_gateway,create_ideal_dm_int,"ideal_dm"},
  {(Myinterfun)sci_gateway,ideal_dm_transform_int,"ideal_dm_transform"},
  {(Myinterfun)sci_gateway,set_dm_timestamp_int,"set_dm_timestamp"},
  {(Myinterfun)sci_gateway,set_dm_actuator_positions_int,"set_dm_actuator_positions"},
  {(Myinterfun)sci_gateway,set_dm_actuator_commands_int,"set_dm_actuator_commands"},
  {(Myinterfun)sci_gateway,three_point_int,"three_point"},
  {(Myinterfun)sci_gateway,three_vector_int,"three_vector"},
  {(Myinterfun)sci_gateway,three_frame_int,"three_frame"},
  {(Myinterfun)sci_gateway,three_reflection_int,"three_reflection"},
  {(Myinterfun)sci_gateway,diffractive_wavefront_header_int,"wavefront_header"},
  {(Myinterfun)sci_gateway,set_dwf_timestamp_int,"set_dwf_timestamp"},
  {(Myinterfun)sci_gateway,set_dwf_frame_int,"set_dwf_frame"},
  {(Myinterfun)sci_gateway,create_ideal_ttm_int,"ideal_ttm"},
  {(Myinterfun)sci_gateway,set_ttm_timestamp_int,"set_ttm_timestamp"},
  {(Myinterfun)sci_gateway,set_ttm_commands_vector_int,"set_ttm_commands_vector"},
  {(Myinterfun)sci_gateway,ideal_ttm_transform_int,"ideal_ttm_transform"},
  {(Myinterfun)sci_gateway,create_zernike_int,"znk_mod"},
  {(Myinterfun)sci_gateway,set_cos_coef_int,"set_znk_cos_coef"},
  {(Myinterfun)sci_gateway,set_sin_coef_int,"set_znk_sin_coef"},
  {(Myinterfun)sci_gateway,zernike_fits_int,"znk_fits"},
  {(Myinterfun)sci_gateway,PI_controller_int,"pi_controller"},
  {(Myinterfun)sci_gateway,ttm_PI_controller_int,"ttm_pi_controller"},
  {(Myinterfun)sci_gateway,dm_PI_controller_int,"dm_pi_controller"},
  {(Myinterfun)sci_gateway,update_ttm_PI_controller_int,"ttm_pi_update"},
  {(Myinterfun)sci_gateway,update_dm_PI_controller_int,"dm_pi_update"},
  {(Myinterfun)sci_gateway,TipTiltMirror_update_int,"ttm_update"},
  {(Myinterfun)sci_gateway,DeformableMirror_update_int,"dm_update"},
  {(Myinterfun)sci_gateway,power_law_int,"power_law"},
  {(Myinterfun)sci_gateway,von_karman_power_law_int,"von_karman_power"},
  {(Myinterfun)sci_gateway,greenwood_power_law_int,"greenwood_power"},
  {(Myinterfun)sci_gateway,null_inner_scale_int,"null_inner"},
  {(Myinterfun)sci_gateway,exponential_inner_scale_int,"exponent_inner"},
  {(Myinterfun)sci_gateway,frehlich_inner_scale_int,"frehlich_inner"},
  {(Myinterfun)sci_gateway,create_power_spectrum_int,"power_spectrum"},
  {(Myinterfun)sci_gateway,null_subharmonic_method_int,"null_subm"},
  {(Myinterfun)sci_gateway,quad_pixel_subharmonic_method_int,"quad_pix_subm"},
  {(Myinterfun)sci_gateway,Lane_subharmonic_method_int,"lane_subm"},
  {(Myinterfun)sci_gateway,Johansson_Gavel_subharmonic_method_int,"johansson_gavel_subm"},
  {(Myinterfun)sci_gateway,generalized_subharmonic_method_int,"general_subm"},
  {(Myinterfun)sci_gateway,create_pixel_array_int,"pixel_array"},
  {(Myinterfun)sci_gateway,set_pixel_array_data_int,"set_pix_arr_data"},
  {(Myinterfun)sci_gateway,ref_atm_lay_int,"ref_atm_lay"},
  {(Myinterfun)sci_gateway,ref_atm_lay_transform_int,"ref_atm_lay_transform"},
  {(Myinterfun)sci_gateway,ref_atm_lay_fits_int,"ref_atm_lay_fits"},
  {(Myinterfun)sci_gateway,Hardy_wind_model_int,"hardy_wind"},
  {(Myinterfun)sci_gateway,get_HardyWind_velocities_int,"get_hardy_wind_vectors"},
  {(Myinterfun)sci_gateway,AtmosphericModel_int,"ref_atm_model"},
  {(Myinterfun)sci_gateway,Ellerbroek_Cerro_Pachon_model_int,"ell_cer_pac_mod"},
  {(Myinterfun)sci_gateway,Ellerbroek_Mauna_Kea_model_int,"ell_mau_kea_mod"},
  {(Myinterfun)sci_gateway,Palomar_DIMM_MASS_model_int,"pal_dimm_mass_mod"},
  {(Myinterfun)sci_gateway,Hufnagel_Valley_model_int,"huf_val_mod"},
  {(Myinterfun)sci_gateway,SLCSAT_day_model_int,"slcsat_day_mod"},
  {(Myinterfun)sci_gateway,SLCSAT_night_model_int,"slcsat_night_mod"},
  {(Myinterfun)sci_gateway,TMT_SRD_v13_Cn2_model_int,"tmt_srd_v13_cn2_mod"},
  {(Myinterfun)sci_gateway,Gemini_GLAO_study_model_int,"gemini_glao_study_mod"},
  {(Myinterfun)sci_gateway,construct_Emitter_int,"emitter"},
  {(Myinterfun)sci_gateway,get_AtmosphericModel_WavefrontHeader_int,"get_atm_mod_dwfh"},
  {(Myinterfun)sci_gateway,get_AtmosphericModel_RefractiveLayer_int,"get_atm_mod_layers"},
  {(Myinterfun)sci_gateway,atm_mod_layers_fits_int,"atm_mod_layers_fits"},
  {(Myinterfun)sci_gateway,create_RECONSTRUCTOR_int,"reconstructor"},
  {(Myinterfun)sci_gateway,arroyo_reconstructor_fits_int,"reconstructor_fits"},
  {(Myinterfun)sci_gateway,arroyo_reconstruct_zernike_residuals_int,"zernike_residuals"},
  {(Myinterfun)sci_gateway,arroyo_reconstruct_zonal_residuals_int,"zonal_residuals"},
  {(Myinterfun)sci_gateway,three_translation_int,"three_translation"},
  {(Myinterfun)sci_gateway,three_rotation_int,"three_rotation"},
  {(Myinterfun)sci_gateway,AtmosphericModel_layer_number_int,"atm_mod_lay_number"},
  {(Myinterfun)sci_gateway,AtmosphericModel_layer_heights_int,"atm_mod_lay_heights"},
  {(Myinterfun)sci_gateway,AtmosphericModel_to_APERTURE_transform_int,"atm_mod_lays_transform"},
};

#endif
