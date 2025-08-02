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

#include "arroyo_wrap.h"

using namespace Arroyo;

//////////////////////////////////////////////////////////////////////////////////
propagation_plan *get_propagation_plan(int type)
{
  propagation_plan * pplan;

  if(type==1)
    pplan = new near_field_fresnel_propagation_plan();
  else 
    pplan = new geometric_propagation_plan();

  return pplan;

}

WaveFront exact_propagate_transform	( double distance,
	double final_pixel_scale, int final_x_axes,	int final_y_axes, WaveFront WF)
{
	std::vector<long> final_axes(2);
	final_axes[0]=final_x_axes; final_axes[1]=final_y_axes;

	WaveFront tmp_WF;
	diffractive_wavefront<double> dwf=c2cpp_WaveFront( WF );

	dwf.exact_propagator (distance,final_pixel_scale,final_axes );

	tmp_WF=cpp2c_WaveFront(dwf);
	return tmp_WF;
}

WaveFront geometric_propagate_transform(double distance, WaveFront WF)
{
	WaveFront tmp_WF;
	diffractive_wavefront<double> dwf=c2cpp_WaveFront( WF );
	dwf.geometric_propagator(distance);
	tmp_WF=cpp2c_WaveFront(dwf);
	return tmp_WF;
}

WaveFront near_field_angular_transform(double distance, WaveFront WF)
{
	WaveFront tmp_WF;
	diffractive_wavefront<double> dwf=c2cpp_WaveFront( WF );
	dwf.near_field_angular_propagator(distance);
	tmp_WF=cpp2c_WaveFront(dwf);
	return tmp_WF;
}

WaveFront near_field_fresnel_transform(double distance, WaveFront WF)
{
	WaveFront tmp_WF;
	diffractive_wavefront<double> dwf=c2cpp_WaveFront( WF );
	dwf.near_field_fresnel_propagator(distance);
	tmp_WF=cpp2c_WaveFront(dwf);
	return tmp_WF;
}

WaveFront far_field_fresnel_transform(double distance, WaveFront WF)
{
	WaveFront tmp_WF;
	diffractive_wavefront<double> dwf=c2cpp_WaveFront( WF );
	dwf.far_field_fresnel_propagator(distance);
	tmp_WF=cpp2c_WaveFront(dwf);
	return tmp_WF;
}

WaveFront far_field_fraunhoffer_transform(double distance, WaveFront WF)
{
	WaveFront tmp_WF;
	diffractive_wavefront<double> dwf=c2cpp_WaveFront( WF );
	dwf.far_field_fraunhoffer_propagator(distance);
	tmp_WF=cpp2c_WaveFront(dwf);
	return tmp_WF;
}

WaveFront far_fresnel_goertzel_reinsch_transform(double distance, 
		double final_pixel_scale, int final_x_axes,int final_y_axes, WaveFront WF )
{
	std::vector<long> final_axes(2);
	final_axes[0]=final_x_axes; final_axes[1]=final_y_axes;

	WaveFront tmp_WF;
	diffractive_wavefront<double> dwf=c2cpp_WaveFront( WF );
	dwf.far_field_fresnel_goertzel_reinsch_propagator(distance,final_pixel_scale,final_axes );
	tmp_WF=cpp2c_WaveFront(dwf);
	return tmp_WF;
}

WaveFront  far_fraunhoffer_goertzel_reinsch_transform(double distance, 
		double final_pixel_scale, int final_x_axes,int final_y_axes, WaveFront WF )
{
	std::vector<long> final_axes(2);
	final_axes[0]=final_x_axes; final_axes[1]=final_y_axes;

	WaveFront tmp_WF;
	diffractive_wavefront<double> dwf=c2cpp_WaveFront( WF );
	dwf.far_field_fraunhoffer_goertzel_reinsch_propagator (distance,final_pixel_scale,final_axes );
	tmp_WF=cpp2c_WaveFront(dwf);
	return tmp_WF;
}
