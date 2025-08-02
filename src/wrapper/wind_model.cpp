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

HardyWind construct_HardyWind( double ground_vel,double tropopause_vel,
		double tropopause_height, double tropopause_thickness )
{
	HardyWind HWM;
	HWM.ground_layer_wind_velocity = ground_vel;
	HWM.tropopause_height = tropopause_height;
	HWM.tropopause_thickness = tropopause_thickness;
	HWM.tropopause_wind_velocity = tropopause_vel;
	return HWM;
}

Hardy_wind_model c2cpp_HardyWind( HardyWind HWM )
{
	return Hardy_wind_model(HWM.ground_layer_wind_velocity,
		HWM.tropopause_height,HWM.tropopause_thickness,
		HWM.tropopause_wind_velocity);
}

ThreeVector *get_HardyWind_velocities( HardyWind HWM,
		    int nlayers, double *height, ThreeFrame TF )
{
	ThreeVector *TVs;
	Hardy_wind_model hwm=c2cpp_HardyWind( HWM );

	vector< double > heights( nlayers );
	for(int i=0;i<nlayers;i++) 	heights[i]=height[i];

	three_frame ref_frame = c2cpp_ThreeFrame( TF );
	vector<three_vector> wind_velocities
		= hwm.get_random_wind_vectors ( heights, ref_frame );
	TVs=(ThreeVector *)calloc( nlayers,sizeof(ThreeVector) );

	for(int i=0;i<nlayers;i++) 
		TVs[i] = cpp2c_ThreeVector( wind_velocities[i] );

	return TVs;
}
