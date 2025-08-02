##
## -------------------------------------------------------------------------
## scitao- Scilab/Scicos Adaptive Optics toolbox
## Copyright (C) 2006  IAPCM, Beijing, China.
## Some inspiration came from SIP's contribution 
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
## See the GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software Foundation,
## Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
## -------------------------------------------------------------------------
##

images  = coordinates.gif ell_cer_pac.gif ell_mau_kea.gif forvard_1.gif forvard_2.gif forvard_3.gif fresnel_1.gif fresnel_2.gif fresnel_3.gif fresnel_4.gif huf_val.gif l_amplify.gif normal.gif scao_ui.gif slcsat_day.gif slcsat_night.gif steps_1.gif steps_2.gif steps_3.gif steps_4.gif strehl.gif tmt_srd_v13_cn2.gif

manpages_xml  = absorber.xml actuator_array.xml annul_ap.xml aperture_transform.xml atm_mod_layers_fits.xml atm_mod_lay_heights.xml atm_mod_lay_number.xml atm_mod_lays_transform.xml begin.xml b_mix.xml b_split.xml circ_ap.xml circle.xml circ_screen.xml convert.xml create_field.xml create_shcentroids.xml cros_out.xml c_scilab.xml dm_pi_controller.xml dm_pi_update.xml dm_update.xml dwf_clip_array.xml dwf_fits.xml ell_cer_pac_mod.xml ell_mau_kea_mod.xml emitter.xml exponent_inner.xml far_fraunhoffer.xml far_fresnel.xml field_contents.xml field_int.xml field_pha.xml field_plot.xml file_int.xml file_pgm.xml file_pha.xml file_ps.xml fil_ter.xml fits_head.xml fits_image_info.xml fitsinfo.xml fitsread.xml fits_structure.xml forvard.xml forward.xml fraunhoffer_grt.xml frehlich_inner.xml fresnel_grt.xml fresnel.xml gauss_screen.xml gauss.xml gemini_glao_study_mod.xml general_subm.xml geom_propagation.xml get_atm_mod_dwfh.xml get_atm_mod_layers.xml get_fits_image.xml get_hardy_wind_vectors.xml get_znk_cos_coef.xml get_znk_sin_coef.xml greenwood_power.xml hardy_wind.xml hex_ap.xml huf_val_mod.xml ideal_dm_transform.xml ideal_dm.xml ideal_ttm_transform.xml ideal_ttm.xml interpol.xml johansson_gavel_subm.xml l_amplify.xml lane_subm.xml lay_pixarr.xml lens_forvard.xml lens_fresnel.xml lens.xml lnslt_array.xml lnslt_arr_transform.xml lp_interp.xml multi_pixarr.xml near_angular.xml near_fresnel.xml normal.xml null_inner.xml null_subm.xml pal_dimm_mass_mod.xml pi_controller.xml pip_fft.xml pixarr_data.xml pixel_array.xml pl_wave.xml power_law.xml power_spectrum.xml quad_pix_subm.xml random.xml reconstructor_fits.xml reconstructor.xml rectang.xml rect_ap.xml rect_screen.xml ref_atm_lay_fits.xml ref_atm_lay_transform.xml ref_atm_lay.xml ref_atm_model.xml scao_simulation.xml set_ap_frame.xml set_dm_actuator_commands.xml set_dm_actuator_positions.xml set_dm_timestamp.xml set_dwf_direction.xml set_dwf_frame.xml set_dwf_pixel_scale.xml set_dwf_timestamp.xml set_lay_frame.xml set_lay_wind.xml set_pix_arr_data.xml set_ttm_commands_vector.xml set_ttm_timestamp.xml set_znk_cos_coef.xml set_znk_sin_coef.xml shc_fits.xml slcsat_day_mod.xml slcsat_night_mod.xml spid_annul_ap.xml sp_wave.xml steps.xml strehl.xml three_frame.xml three_point.xml three_reflection.xml three_rotation.xml three_translation.xml three_vector.xml til_hex_ap.xml tilt.xml tmt_srd_v13_cn2_mod.xml tor_lens.xml ttm_pi_controller.xml ttm_pi_update.xml ttm_update.xml unfold_phase.xml von_karman_power.xml wavefront_header.xml zernike_residuals.xml zernike.xml znk_fits.xml znk_mod.xml zonal_residuals.xml

manpages_html = $(manpages_xml:.xml=.htm)

helpdir=@TOOLBOXDIR@/man/eng/optics
help_DATA = $(manpages_html) $(images) whatis.htm
EXTRA_DIST = $(manpages_xml) $(images) scilab_rev.xsl

.xml.htm: $(manpages_xml)
	@echo "Creating html help files"
	@$(SCIDIR)/bin/scilab -nwni -nb -e "exec('../../sao_xmltohtml.sci');%helps=[PWD,'Optics';%helps];sao_xmltohtml(PWD,'Optics','scilab_rev.xsl');quit"

clean-local:
	rm -f *.htm
