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

// New functions
lightPipes_interf=[..
           'begin';..
           'circ_ap';..
           'circ_screen'; ..
           'absorber';..
           'cros_out'; ..
           'file_pgm';..
           'file_ps'; ..
           'file_int'; ..
           'file_pha'; ..
           'forvard'; ..
           'fresnel';  ..
           'forward';..
           'gauss'; ..
           'gauss_screen';..
           'rect_ap';..
           'rect_screen';..
           'random';..
           'zernike';..
           'l_amplify';..
           'lens';..
           'lens_forvard';..
           'lens_fresnel';..
           'convert'; ..
           'normal';..
           'tilt'; ..
           'tor_lens';..
           'pip_fft';..
           'interpol';
           'lp_interp';..
           'b_mix'; ..
           'b_split'; ..
           'fil_ter';..
           'steps'; ..
           'strehl'; ..
           'field_int'; ..
           'field_pha'; ..
           'create_field'; ..
           'field_contents'; ..
           'unfold_phase'; ..
           'c_scilab'];

arroyo_interf=[..
           'dwf_fits';..
           'circle';..
           'rectang'; ..
           'annul_ap';..
           'hex_ap'; ..
           'spid_annul_ap';..
           'til_hex_ap'; ..
           'pl_wave'; ..
           'sp_wave';..
           'geom_propagation'; ..
           'near_angular'; ..
           'near_fresnel'; ..
           'far_fresnel';..
           'far_fraunhoffer';..
           'fresnel_grt';..
           'fraunhoffer_grt'; ..
           'lnslt_array';..
           'lnslt_arr_transform'; ..
           'create_shcentroids'; ..
           'shc_fits'; ..
           'aperture_transform';..
           'actuator_array'; ..
           'ideal_dm'; ..
           'ideal_dm_transform'; ..
           'set_dm_timestamp';..
           'set_dm_actuator_positions'; ..
           'set_dm_actuator_commands'; ..
           'three_point';..
           'three_vector'; ..
           'three_frame'; ..
           'three_reflection';..
           'wavefront_header';..
           'set_dwf_timestamp'; ..
           'set_dwf_frame';..
           'ideal_ttm'; ..
           'set_ttm_timestamp';..
           'set_ttm_commands_vector'; ..
           'ideal_ttm_transform';..
           'znk_mod'; ..
           'set_znk_cos_coef';..
           'set_znk_sin_coef'; ..
           'znk_fits'; ..
           'pi_controller'; ..
           'ttm_pi_controller';..
           'dm_pi_controller'; ..
           'ttm_pi_update'; ..
           'dm_pi_update'; ..
           'ttm_update';..
           'dm_update'; ..
           'power_law'; ..
           'von_karman_power'; ..
           'greenwood_power';..
           'null_inner'; ..
           'exponent_inner'; ..
           'frehlich_inner'; ..
           'power_spectrum'; ..
           'null_subm'; ..
           'quad_pix_subm'; ..
           'lane_subm'; ..
           'johansson_gavel_subm';..
           'general_subm'; ..
           'pixel_array'; ..
           'set_pix_arr_data'; ..
           'ref_atm_lay';..
           'ref_atm_lay_transform'; ..
           'ref_atm_lay_fits'; ..
           'hardy_wind';..
           'get_hardy_wind_vectors'; ..
           'ref_atm_model';..
           'ell_cer_pac_mod';..
           'ell_mau_kea_mod'; ..
           'pal_dimm_mass_mod'; ..
           'huf_val_mod';..
           'slcsat_day_mod'; ..
           'slcsat_night_mod'; ..
           'tmt_srd_v13_cn2_mod'; ..
           'gemini_glao_study_mod';..
           'emitter'; ..
           'get_atm_mod_dwfh'; ..
           'get_atm_mod_layers'; ..
           'atm_mod_layers_fits';..
           'reconstructor'; ..
           'reconstructor_fits'; ..
           'zernike_residuals';..
           'zonal_residuals'; ..
           'three_translation'; ..
           'three_rotation'; ..
           'atm_mod_lay_number';..
           'atm_mod_lay_heights'; ..
           'atm_mod_lays_transform'];

other_interf=[..
           'fits_head';..
           'fits_structure';..
           'fits_image_info';..
           'get_fits_image';..
           'scao_simulation';..
           'pixarr_data';..
           'lay_pixarr';..
           'set_ap_frame';..
           'set_dwf_direction';..
           'dwf_clip_array';..
           'set_dwf_pixel_scale';..
           'get_znk_cos_coef';..
           'get_znk_sin_coef';..
           'multi_pixarr';..
           'set_lay_wind';..
           'set_lay_frame'];

scicos_functions=['start';..
           'absorber'; ..
           'mixer';..
           'normal'; ..
           'controller'; ..
           'Lens'; ..
           'splitter';..
           'convert'; ..
           'diagnose'; ..
           'display'; ..
           'zernike';..
           'pip_fft'; ..
           'hartmann'; ..
           'emitter'; ..
           'file_fits';..
           'atmosphere'; ..
           'geometry'; ..
           'propagate';..
           'df_mirror';..
           'aperture'; ..
           'interpolate';..
           'tt_mirror';..
           'gauss'; ..
           'reconstructor'];

[units,typs,nams]=file();
clear units typs
for k=size(nams,'*'):-1:1
	l=strindex(nams(k),'loader.sce');
	if l<>[] then
		DIR=part(nams(k),1:l($)-1);
		break
	end
end

if ~MSDOS then
	if part(DIR,1)<>'/' then
		DIR=getcwd()+'/'+DIR
	end

	SRC=DIR+'src/'
	SCICOS=SRC+'scicos/'
	INTERF=SRC+'interf/'
	
	MACROS=DIR+'macros/'
	
	DEMOS=DIR+'demos/'
	EXAMPLES=DIR+'examples/'
	MAN=DIR+'man/'

else

 	if part(DIR,2)<>':' then
		DIR=getcwd()+'\'+DIR
	end

	SRC=DIR+'src\'
	SCICOS=SRC+'scicos\'
	INTERF=SRC+'interf\'
	
	MACROS=DIR+'macros\'
	
	DEMOS=DIR+'demos\' 
	EXAMPLES=DIR+'examples\'
	MAN=DIR+'man\'

end

disp('load scitao toolbox');
//---------------------------------------------------------------------
// Under Unix and Linux systems
//---------------------------------------------------------------------
if ~MSDOS then

	link(+DIR+'lib/libarroyo.so');
	link(+DIR+'lib/libarroyo_wrap.so');
	link(+DIR+'lib/liblightpipes.so');

	link(+DIR+'lib/libscicos.so',scicos_functions,'C');
	addinter(+DIR+'lib/libinterf.so','int_others',other_interf);
	addinter(+DIR+'lib/libinterf.so','int_arroyo',arroyo_interf);
	addinter(+DIR+'lib/libinterf.so','int_lightPipes',lightPipes_interf);

//---------------------------------------------------------------------
// Under Windows systems
//---------------------------------------------------------------------
else

	link(+SCICOS+'scicos.dll',scicos_functions,'C');
	addinter(+INTERF+'interf.dll','int_others',other_interf);
	addinter(+INTERF+'interf.dll','int_arroyo',arroyo_interf);
	addinter(+INTERF+'interf.dll','int_lightPipes',lightPipes_interf);
end

// load the macros library
genlib('optics',MACROS);
load(+MACROS+'lib');

// Add palette

chdir(MACROS);
errcatch(-1,'continue','nomessage');
load('.scicos_pal')
if iserror(-1) then errclear(-1);end
scicos_pal=[scicos_pal; ['Optics',+MACROS+'Optics.cosf'] ]	
[x,k]=gsort(scicos_pal(:,2))
keq=find(x(2:$)==x(1:$-1))
k(keq)=[];k=-sort(-k);
scicos_pal=scicos_pal(k,:);
errcatch(-1,'continue')
save('.scicos_pal',scicos_pal)
if iserror(-1) then
  errclear(-1)
  message('Cannot save .scicos_pal in current directory')
end

chdir(MAN)
exec('loader.sce');

demotitle="Wave/Adaptive Optics"
execstr("global demolist","errcatch") 

if ~or(demolist(:,1)==demotitle) then
  add_demo(demotitle, DIR+pathconvert("/demos/optics.dem",%f))
end

disp('scitao is loaded !!!');
disp('Enjoy it !!!');
disp('');

chdir(EXAMPLES);
stacksize(5e7);

clear lightPipes_interf arroyo_interf other_interf scicos_functions demotitle;
clear LIBNAME DIR SRC SCICOS INTERF MACROS MAN DEMOS EXAMPLES;  

