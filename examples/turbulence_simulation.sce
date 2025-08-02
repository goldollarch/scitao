
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

outer_diameter = 2;
wavelength = 2.2e-6;
focal_plane_image_size_arcsecs = 4;
oversampling_factor = 4;
nsteps = 30;timestep_secs = .01;
wf_pixel_scale_meters = .02;
layer_pixel_scale_meters = .02;
exponent = -11/3.0;r_0_meters = .15;
r_0_ref_wavelength_meters = .5e-6;
seed = 1;subharmonic_depth = 3;
wind_velocity_meters_per_sec = 10;
   
tf=three_frame();
ap = circle(outer_diameter);

pl=power_law(r_0_meters,r_0_ref_wavelength_meters,exponent);
nin=null_inner();ipls=power_spectrum(pl,nin);

subm = lane_subm(subharmonic_depth);   
dimen_x = ceil(outer_diameter/layer_pixel_scale_meters);
dimen_y = ceil((outer_diameter+abs(wind_velocity_meters_per_sec)*nsteps*timestep_secs)/layer_pixel_scale_meters);
layer_axes_x=ceil(1.1*dimen_x);layer_axes_y=ceil(1.1*dimen_y);

ref_layer=ref_atm_lay(layer_pixel_scale_meters,layer_axes_x,seed,layer_axes_y,ipls,subm);
tv_wind=three_vector(0,wind_velocity_meters_per_sec,0);
ref_layer=set_lay_wind(ref_layer,tv_wind);

translation_vector = three_vector(0,0.5*nsteps*timestep_secs*wind_velocity_meters_per_sec,0);
tf=three_translation(tf,translation_vector);
ref_layer=set_lay_frame(ref_layer,tf);

//ref_atm_lay_fits(ref_layer,"ref_layer");

wf_axes=ceil(outer_diameter/wf_pixel_scale_meters);
dwfh=wavefront_header(wavelength,wf_pixel_scale_meters,wf_axes);

initial_dwf = begin(outer_diameter,wavelength,wf_axes);
//dwf_fits(initial_dwf,"initial_dwf");

propagation_distance = 1e10;
final_axes=ceil(oversampling_factor*outer_diameter/wf_pixel_scale_meters);
final_pixel_scale = wavelength*propagation_distance/final_axes/wf_pixel_scale_meters/oversampling_factor;

dwf=aperture_transform(initial_dwf,ap);
unaberrated_dwf=fresnel_grt(dwf,propagation_distance,final_pixel_scale,final_axes);
//dwf_fits(unaberrated_dwf,"unaberrated_focal_wf");

x=1:wf_axes;y=x;
grid=[1,1,wf_axes,wf_axes];

for i=1:nsteps
	time=i*timestep_secs;
	dwf=set_dwf_timestamp(initial_dwf,time);

	dwf=ref_atm_lay_transform(dwf,ref_layer);
	dwf=aperture_transform(dwf,ap);
//	dwf_fits(dwf,"pupil_wf",time);

	pha=field_pha(dwf);
	grayplot(x,y,pha,strf="030",rect=grid);
	
//	tmp_dwf=fresnel_grt(dwf,propagation_distance,final_pixel_scale,final_axes);
//	dwf_fits(tmp_dwf,"focal_wf",time);

end