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

#ifndef INT_LIGHTPIPES_H
#define INT_LIGHTPIPES_H

extern Gatefunc begin_int;
extern Gatefunc circ_ap_int;
extern Gatefunc circ_screen_int;
extern Gatefunc absorber_int;
extern Gatefunc cros_out_int;
extern Gatefunc file_pgm_int;
extern Gatefunc file_ps_int;
extern Gatefunc file_int_int;
extern Gatefunc file_pha_int;
extern Gatefunc forvard_int;
extern Gatefunc fresnel_int;
extern Gatefunc forward_int;
extern Gatefunc gauss_int;
extern Gatefunc gauss_screen_int;
extern Gatefunc rect_ap_int;
extern Gatefunc rect_screen_int;
extern Gatefunc random_int;
extern Gatefunc lp_zernike_int;
extern Gatefunc l_amplif_int;
extern Gatefunc lens_int;
extern Gatefunc lens_forvard_int;
extern Gatefunc lens_fresn_int;
extern Gatefunc convert_int;
extern Gatefunc normal_int;
extern Gatefunc tilt_int;
extern Gatefunc tor_lens_int;
extern Gatefunc pip_fft_int;
extern Gatefunc interpol_int;
extern Gatefunc interp1_int;
extern Gatefunc b_mix_int;
extern Gatefunc b_split_int;
extern Gatefunc fil_ter_int;
extern Gatefunc steps_int;
extern Gatefunc strehl_int;
extern Gatefunc field_int_int;
extern Gatefunc field_pha_int;
extern Gatefunc create_field_int;
extern Gatefunc field_contents_int;
extern Gatefunc unfold_phase_int;
extern Gatefunc c_scilab_int;

static GenericTable lightPipes_Tab[]={
  {(Myinterfun)sci_gateway,begin_int,"begin"},
  {(Myinterfun)sci_gateway,circ_ap_int,"circ_ap"},
  {(Myinterfun)sci_gateway,circ_screen_int,"circ_screen"},
  {(Myinterfun)sci_gateway,absorber_int,"absorber"},
  {(Myinterfun)sci_gateway,cros_out_int,"cros_out"},
  {(Myinterfun)sci_gateway,file_pgm_int,"file_pgm"},
  {(Myinterfun)sci_gateway,file_ps_int,"file_ps"},
  {(Myinterfun)sci_gateway,file_int_int,"file_int"},
  {(Myinterfun)sci_gateway,file_pha_int,"file_pha"},
  {(Myinterfun)sci_gateway,forvard_int,"forvard"},
  {(Myinterfun)sci_gateway,fresnel_int,"fresnel"},
  {(Myinterfun)sci_gateway,forward_int,"forward"},
  {(Myinterfun)sci_gateway,gauss_int,"gauss"},
  {(Myinterfun)sci_gateway,gauss_screen_int,"gauss_screen"},
  {(Myinterfun)sci_gateway,rect_ap_int,"rect_ap"},
  {(Myinterfun)sci_gateway,rect_screen_int,"rect_screen"},
  {(Myinterfun)sci_gateway,random_int,"random"},
  {(Myinterfun)sci_gateway,lp_zernike_int,"zernike"},
  {(Myinterfun)sci_gateway,l_amplif_int,"l_amplify"},
  {(Myinterfun)sci_gateway,lens_int,"lens"},
  {(Myinterfun)sci_gateway,lens_forvard_int,"lens_forvard"},
  {(Myinterfun)sci_gateway,lens_fresn_int,"lens_fresnel"},
  {(Myinterfun)sci_gateway,convert_int,"convert"},
  {(Myinterfun)sci_gateway,normal_int,"normal"},
  {(Myinterfun)sci_gateway,tilt_int,"tilt"},
  {(Myinterfun)sci_gateway,tor_lens_int,"tor_lens"},
  {(Myinterfun)sci_gateway,pip_fft_int,"pip_fft"},
  {(Myinterfun)sci_gateway,interpol_int,"interpol"},
  {(Myinterfun)sci_gateway,interp1_int,"interp"},
  {(Myinterfun)sci_gateway,b_mix_int,"b_mix"},
  {(Myinterfun)sci_gateway,b_split_int,"b_split"},
  {(Myinterfun)sci_gateway,fil_ter_int,"fil_ter"},
  {(Myinterfun)sci_gateway,steps_int,"steps"},
  {(Myinterfun)sci_gateway,strehl_int,"strehl"},
  {(Myinterfun)sci_gateway,field_int_int,"field_int"},
  {(Myinterfun)sci_gateway,field_pha_int,"field_pha"},
  {(Myinterfun)sci_gateway,create_field_int,"create_field"},
  {(Myinterfun)sci_gateway,field_contents_int,"field_contents"},
  {(Myinterfun)sci_gateway,unfold_phase_int,"unfold_phase"},
  {(Myinterfun)sci_gateway,c_scilab_int,"c_scilab"},
};

#endif
