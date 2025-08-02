
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

nsteps = 16;
step_distance = 1000;

wavelength = 2.2e-6;
wf_pixel_scale = .01;
layer_pixel_scale = .01;

seed=1;
outer_diameter = 1;

exponent = -11/3.0;
r_0_meters = .15;
r_0_ref_wavelength_meters = .5e-6;
    
subharmonic_depth = 3;

pl=power_law(r_0_meters,r_0_ref_wavelength_meters,exponent);
nin=null_inner();ipls=power_spectrum(pl,nin);
subm = lane_subm(subharmonic_depth);     

wf_padding = ceil(wavelength*nsteps*step_distance/wf_pixel_scale/wf_pixel_scale);
layer_padding = ceil(wavelength*nsteps*step_distance/layer_pixel_scale/layer_pixel_scale);
  
dimen = ceil(outer_diameter/layer_pixel_scale + 2*layer_padding);
layer_axes=ceil(1.05*dimen);

ref_layer=ref_atm_lay(layer_pixel_scale,layer_axes,seed,layer_axes,ipls,subm);
//ref_atm_lay_fits(ref_layer,"ref_layer");

tp=three_point(0,0,0);
tv_x=three_vector(1,0,0);tv_y=three_vector(0,1,0);tv_z=three_vector(0,0,-1);
wf_frame=three_frame(tp,tv_x,tv_y,tv_z);
N=ceil(outer_diameter/wf_pixel_scale);
wf_axes = N+2*wf_padding;

dwfh=wavefront_header(wavelength,wf_pixel_scale,wf_axes,wf_axes,wf_frame);

pwemtr=emitter(0,0,0,-1);
dwf = pl_wave(dwfh,pwemtr);

initial_dwf=ref_atm_lay_transform(dwf,ref_layer);
//dwf_fits(initial_dwf,"initial_dwf");

xbasc();
drawlater();
ncolor=60;
h=hotcolormap(ncolor);
xset('colormap',h);

for i=1:nsteps
	dwf=near_fresnel(initial_dwf,i*step_distance); 
	dwf=dwf_clip_array(dwf,wf_padding);
//	dwf_fits(dwf,"dwf",i);

	I=field_int(dwf);
	subplot(4,4,i);
	grayplot(1:N,1:N,I,strf="030",rect=[1,1,N,N]);
	Str=sprintf('z=%e m',i*step_distance)
	xtitle(Str);
end

drawnow();
