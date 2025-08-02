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

#include "optics_scicos.h"

void propagate(scicos_block *block,int flag)
{
	int *ipar; 
	double *rpar;

	double distance,focal_pixel;
	int type,final_x_axes,final_y_axes;

    double propagation_distance,final_pixel_scale;
	int final_array_dimensions;
	
	FIELD field;
	WaveFront WF_in, WF_out;

   /* reading the data from interface */

	rpar=block->rpar;
	--rpar;
	ipar=block->ipar;
	--ipar;

    type = ipar[1];
	final_x_axes = ipar[2];
	final_y_axes = ipar[3];

    distance = rpar[1];
	focal_pixel = rpar[2];

	if(flag==1)
	{
		field=block_inptr_field(block,0);

		if(field.number) 
		{
			WF_in = FIELD_WaveFront(field);

			switch(type) 
			{
			case 1:
				field = lp_forvard( field,distance );
				break;
		    case 2:
				field = lp_fresnel( field,distance );
				break;

		    case 3:
				WF_out = near_field_angular_transform( distance, WF_in );
				WaveFront_FIELD(WF_out,&field);
				break;
		    case 4:
				WF_out = near_field_fresnel_transform( distance, WF_in );
				WaveFront_FIELD(WF_out,&field);
				break;
		    case 5:
				WF_out = far_field_fresnel_transform( distance, WF_in );
				WaveFront_FIELD(WF_out,&field);
				break;
		    case 6:
				WF_out = far_field_fraunhoffer_transform( distance, WF_in );
				WaveFront_FIELD(WF_out,&field);
				break;
		    case 7:
				WF_out = far_fresnel_goertzel_reinsch_transform( distance,
					focal_pixel, final_x_axes, final_y_axes, WF_in );
				WaveFront_FIELD(WF_out,&field);
				break;
		    case 8:
				WF_out = far_fraunhoffer_goertzel_reinsch_transform( distance,
					focal_pixel, final_x_axes, final_y_axes, WF_in );
				WaveFront_FIELD(WF_out,&field);
				break;

			case 9:
				field = lp_lens_forvard( field,focal_pixel,distance );
				break;
		    case 10:
				field = lp_lens_fresn( field,focal_pixel,distance );
				break;

		    case 11:

				propagation_distance=1e10;
				final_array_dimensions=(int)(scao.global_AP.ap_data1/WF_in.WfH.pixscale);
				final_pixel_scale = WF_in.WfH.lambda * propagation_distance / 
					(double) final_array_dimensions / WF_in.WfH.pixscale;

				WF_out = far_fresnel_goertzel_reinsch_transform( 
					propagation_distance,final_pixel_scale, 
					final_array_dimensions, final_array_dimensions,
					WF_in );

				WaveFront_FIELD(WF_out,&field);

				break;

		    default:
				WF_out = geometric_propagate_transform( distance,WF_in );
				WaveFront_FIELD(WF_out,&field);
				break;
			}
			
			field_block_outptr(field,block,0);

		}

	}

}
