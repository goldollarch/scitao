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

#ifndef INT_OTHERS_H
#define INT_OTHERS_H

extern Gatefunc list_fits_head_int;
extern Gatefunc list_fits_structure_int;
extern Gatefunc fits_image_info_int;
extern Gatefunc get_fits_image_int;

extern Gatefunc scao_simulation_int;
extern Gatefunc pixarr_data_int;
extern Gatefunc lay_pixarr_int;
extern Gatefunc ap_frame_int;
extern Gatefunc set_dwf_direction_int;
extern Gatefunc dwf_clip_array_int;
extern Gatefunc set_dwf_pixel_scale_int;
extern Gatefunc get_cos_coef_int;
extern Gatefunc get_sin_coef_int;
extern Gatefunc multi_pixarr_int;
extern Gatefunc set_lay_wind_int;
extern Gatefunc set_lay_frame_int;

static GenericTable others_Tab[]={
  {(Myinterfun)sci_gateway,list_fits_head_int,"fits_head"},
  {(Myinterfun)sci_gateway,list_fits_structure_int,"fits_structure"},
  {(Myinterfun)sci_gateway,fits_image_info_int,"fits_image_info"},
  {(Myinterfun)sci_gateway,get_fits_image_int,"get_fits_image"},
  {(Myinterfun)sci_gateway,scao_simulation_int,"scao_simulation"},
  {(Myinterfun)sci_gateway,pixarr_data_int,"pixarr_data"},
  {(Myinterfun)sci_gateway,lay_pixarr_int,"lay_pixarr"},
  {(Myinterfun)sci_gateway,ap_frame_int,"set_ap_frame"},
  {(Myinterfun)sci_gateway,set_dwf_direction_int,"set_dwf_direction"},
  {(Myinterfun)sci_gateway,dwf_clip_array_int,"dwf_clip_array"},
  {(Myinterfun)sci_gateway,set_dwf_pixel_scale_int,"set_dwf_pixel_scale"},
  {(Myinterfun)sci_gateway,get_cos_coef_int,"get_znk_cos_coef"},
  {(Myinterfun)sci_gateway,get_sin_coef_int,"get_znk_sin_coef"},
  {(Myinterfun)sci_gateway,multi_pixarr_int,"multi_pixarr"},
  {(Myinterfun)sci_gateway,set_lay_wind_int,"set_lay_wind"},
  {(Myinterfun)sci_gateway,set_lay_frame_int,"set_lay_frame"},
};

#endif
