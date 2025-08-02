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

DeformableMirror create_DeformableMirror( APERTURE AP, int actuator_x_dim, 
		int actuator_y_dim, double actuator_pitch, double actuator_velocity )
{
	DeformableMirror dm;

	dm.timestamp=0;

	dm.actuator_x_dim=actuator_x_dim; 
	dm.actuator_y_dim=actuator_y_dim;

	dm.AP = AP;
	dm.actuator_pitch=actuator_pitch;
	dm.actuator_velocity=actuator_velocity;

	dm.actuator_commands=(double *) calloc( actuator_x_dim*actuator_y_dim, sizeof(double) );
	dm.actuator_positions=(double *) calloc( actuator_x_dim*actuator_y_dim, sizeof(double) );

	return dm;
}

void set_dm_timestamp( DeformableMirror* dm, double timestamp)
{
	dm->timestamp=timestamp;
}
void set_actuator_positions( DeformableMirror* dm, double *actuator_positions)
{ 
	dm->actuator_positions=actuator_positions;
}
void set_actuator_commands( DeformableMirror* dm, double *actuator_commands )
{
	dm->actuator_commands=actuator_commands;
}

WaveFront DeformableMirror_transform( DeformableMirror* dm, WaveFront WF )
{
	int nspiders;

	double diameter,in_diameter,out_diameter,
		x_size,y_size,in_edge_length,spider_width,in_gap_size;

	diffractive_wavefront<double> dwf=c2cpp_WaveFront( WF );

	Arroyo::three_frame oriented_frame=c2cpp_ThreeFrame(dm->AP.TF);

	vector<long> actuator_dimensions(2);
	actuator_dimensions[0]=dm->actuator_x_dim; 
	actuator_dimensions[1]=dm->actuator_y_dim;

	long index;
	pixel_array<double> actuator_commands(actuator_dimensions);
	pixel_array<double> actuator_positions(actuator_dimensions);

	for(int i=0; i<dm->actuator_x_dim; i++) {
		for(int j=0; j<dm->actuator_y_dim; j++) {
			index=i*dm->actuator_x_dim+j;
			actuator_commands.set_data(index,dm->actuator_commands[index]);
			actuator_positions.set_data(index,dm->actuator_positions[index]);
		}
	}
	
	switch(dm->AP.aperture_type)  {
		case 1: 
			{
				in_diameter=dm->AP.ap_data2;
				out_diameter=dm->AP.ap_data1;
				ideal_deformable_mirror<annular_aperture>
					deformable_mirror (annular_aperture(in_diameter,out_diameter),
					actuator_dimensions, dm->actuator_pitch, dm->actuator_velocity,
					dm->timestamp);

				deformable_mirror.Arroyo::three_frame::operator=(oriented_frame);
				deformable_mirror.set_mirror_actuator_positions(actuator_positions);
				deformable_mirror.set_mirror_actuator_commands(actuator_commands);

				deformable_mirror.transform(dwf);
				break;
			}
		case 2:
			{
				x_size=dm->AP.ap_data1;
				y_size=dm->AP.ap_data2;
				ideal_deformable_mirror<rectangular_aperture> 
					deformable_mirror (rectangular_aperture(x_size,y_size),
					actuator_dimensions, dm->actuator_pitch, 
					dm->actuator_velocity,dm->timestamp);

				deformable_mirror.Arroyo::three_frame::operator=(oriented_frame);
				deformable_mirror.set_mirror_actuator_positions(actuator_positions);
				deformable_mirror.set_mirror_actuator_commands(actuator_commands);

				deformable_mirror.transform(dwf);
				break;
			}
		case 3:
			{
				in_edge_length=dm->AP.ap_data1;
				ideal_deformable_mirror<hexagonal_aperture> 
					deformable_mirror (hexagonal_aperture(in_edge_length),
					actuator_dimensions, dm->actuator_pitch, 
					dm->actuator_velocity,dm->timestamp);

				deformable_mirror.Arroyo::three_frame::operator=(oriented_frame);
				deformable_mirror.set_mirror_actuator_positions(actuator_positions);
				deformable_mirror.set_mirror_actuator_commands(actuator_commands);

				deformable_mirror.transform(dwf);
				break;
			}
		case 4:
			{
				in_diameter=dm->AP.ap_data2;
				out_diameter=dm->AP.ap_data1;
				nspiders=(int)(dm->AP.ap_data3); 
				spider_width=dm->AP.ap_data4;
				ideal_deformable_mirror<spidered_annular_aperture> 
					deformable_mirror (spidered_annular_aperture
					(in_diameter,out_diameter,nspiders,spider_width),
					actuator_dimensions, dm->actuator_pitch, 
					dm->actuator_velocity,dm->timestamp);

				deformable_mirror.Arroyo::three_frame::operator=(oriented_frame);
				deformable_mirror.set_mirror_actuator_positions(actuator_positions);
				deformable_mirror.set_mirror_actuator_commands(actuator_commands);

				deformable_mirror.transform(dwf);
				break;
			}

		default:
			{
				diameter=dm->AP.ap_data1;
				ideal_deformable_mirror<circular_aperture> 
					deformable_mirror ( circular_aperture(diameter),
					actuator_dimensions, dm->actuator_pitch, 
					dm->actuator_velocity,dm->timestamp);

				deformable_mirror.Arroyo::three_frame::operator=(oriented_frame);
				deformable_mirror.set_mirror_actuator_positions(actuator_positions);
				deformable_mirror.set_mirror_actuator_commands(actuator_commands);

				deformable_mirror.transform(dwf);
				break;
			}
	}
    
	WaveFront tmp_WF=cpp2c_WaveFront(dwf);
	return tmp_WF;
}

void DeformableMirror_updata( DeformableMirror* dm, 
		double *actuator_commands_ptr, double timestamp )
{
	int nspiders;

	double diameter,in_diameter,out_diameter,
		x_size,y_size,in_edge_length,spider_width,in_gap_size;

	Arroyo::three_frame oriented_frame=c2cpp_ThreeFrame(dm->AP.TF);

	vector<long> actuator_dimensions(2);
	actuator_dimensions[0]=dm->actuator_x_dim; 
	actuator_dimensions[1]=dm->actuator_y_dim;

	long index;
	pixel_array<double> actuator_commands(actuator_dimensions);
	pixel_array<double> actuator_positions(actuator_dimensions);
	PixelArray tmp_commands_PixArr,tmp_positions_PixArr;

	for(int i=0; i<dm->actuator_x_dim; i++) {
		for(int j=0; j<dm->actuator_y_dim; j++) {
			index=i*dm->actuator_x_dim+j;
			actuator_commands.set_data(index,dm->actuator_commands[index]);
			actuator_positions.set_data(index,dm->actuator_positions[index]);
		}
	}

	PixelArray actuator_commands_PixArr=create_PixelArray
		( actuator_dimensions[0], actuator_dimensions[1],0);
	set_PixelArray_data(&actuator_commands_PixArr,actuator_commands_ptr);
	pixel_array<double> new_commands=c2cpp_PixelArray(actuator_commands_PixArr);
	
	switch(dm->AP.aperture_type)  {
		case 1: 
			{
				in_diameter=dm->AP.ap_data2;
				out_diameter=dm->AP.ap_data1;
				ideal_deformable_mirror<annular_aperture>
					deformable_mirror (annular_aperture(in_diameter,out_diameter),
					actuator_dimensions, dm->actuator_pitch, dm->actuator_velocity,
					dm->timestamp);

				deformable_mirror.Arroyo::three_frame::operator=(oriented_frame);
				deformable_mirror.set_mirror_actuator_positions(actuator_positions);
				deformable_mirror.set_mirror_actuator_commands(actuator_commands);

				deformable_mirror.update(new_commands,timestamp);
				tmp_commands_PixArr=cpp2c_PixelArray(deformable_mirror.get_actuator_commands());
				tmp_positions_PixArr=cpp2c_PixelArray(deformable_mirror.get_actuator_positions(timestamp));

				break;
			}
		case 2:
			{
				x_size=dm->AP.ap_data1;
				y_size=dm->AP.ap_data2;
				ideal_deformable_mirror<rectangular_aperture> 
					deformable_mirror (rectangular_aperture(x_size,y_size),
					actuator_dimensions, dm->actuator_pitch, 
					dm->actuator_velocity,dm->timestamp);

				deformable_mirror.Arroyo::three_frame::operator=(oriented_frame);
				deformable_mirror.set_mirror_actuator_positions(actuator_positions);
				deformable_mirror.set_mirror_actuator_commands(actuator_commands);

				deformable_mirror.update(new_commands,timestamp);
				tmp_commands_PixArr=cpp2c_PixelArray(deformable_mirror.get_actuator_commands());
				tmp_positions_PixArr=cpp2c_PixelArray(deformable_mirror.get_actuator_positions(timestamp));

				break;
			}
		case 3:
			{
				in_edge_length=dm->AP.ap_data1;
				ideal_deformable_mirror<hexagonal_aperture> 
					deformable_mirror (hexagonal_aperture(in_edge_length),
					actuator_dimensions, dm->actuator_pitch, 
					dm->actuator_velocity,dm->timestamp);

				deformable_mirror.Arroyo::three_frame::operator=(oriented_frame);
				deformable_mirror.set_mirror_actuator_positions(actuator_positions);
				deformable_mirror.set_mirror_actuator_commands(actuator_commands);

				deformable_mirror.update(new_commands,timestamp);
				tmp_commands_PixArr=cpp2c_PixelArray(deformable_mirror.get_actuator_commands());
				tmp_positions_PixArr=cpp2c_PixelArray(deformable_mirror.get_actuator_positions(timestamp));

				break;
			}
		case 4:
			{
				in_diameter=dm->AP.ap_data2;
				out_diameter=dm->AP.ap_data1;
				nspiders=(int)(dm->AP.ap_data3); 
				spider_width=dm->AP.ap_data4;
				ideal_deformable_mirror<spidered_annular_aperture> 
					deformable_mirror (spidered_annular_aperture
					(in_diameter,out_diameter,nspiders,spider_width),
					actuator_dimensions, dm->actuator_pitch, 
					dm->actuator_velocity,dm->timestamp);

				deformable_mirror.Arroyo::three_frame::operator=(oriented_frame);
				deformable_mirror.set_mirror_actuator_positions(actuator_positions);
				deformable_mirror.set_mirror_actuator_commands(actuator_commands);

				deformable_mirror.update(new_commands,timestamp);
				tmp_commands_PixArr=cpp2c_PixelArray(deformable_mirror.get_actuator_commands());
				tmp_positions_PixArr=cpp2c_PixelArray(deformable_mirror.get_actuator_positions(timestamp));

				break;
			}

		default:
			{
				diameter=dm->AP.ap_data1;
				ideal_deformable_mirror<circular_aperture> 
					deformable_mirror ( circular_aperture(diameter),
					actuator_dimensions, dm->actuator_pitch, 
					dm->actuator_velocity,dm->timestamp);

				deformable_mirror.Arroyo::three_frame::operator=(oriented_frame);
				deformable_mirror.set_mirror_actuator_positions(actuator_positions);
				deformable_mirror.set_mirror_actuator_commands(actuator_commands);

				deformable_mirror.update(new_commands,timestamp);
				tmp_commands_PixArr=cpp2c_PixelArray(deformable_mirror.get_actuator_commands());
				tmp_positions_PixArr=cpp2c_PixelArray(deformable_mirror.get_actuator_positions(timestamp));

				break;
			}
	}
	
	set_dm_timestamp(dm,timestamp);
	set_actuator_commands(dm,tmp_commands_PixArr.pixeldata);
	set_actuator_positions(dm,tmp_positions_PixArr.pixeldata);

}

DeformableMirror array2DeformableMirror( double *inptr )
{
	int nactuator;
	DeformableMirror dm;

	dm.actuator_x_dim=(int)inptr[0];
	dm.actuator_y_dim=(int)inptr[1];
	nactuator=(dm.actuator_x_dim)*(dm.actuator_y_dim);

	dm.AP=array2APERTURE(inptr+2);

	dm.actuator_pitch=inptr[19];
	dm.actuator_velocity=inptr[20];
	dm.timestamp=inptr[21];

	dm.actuator_commands=inptr+22;
	dm.actuator_positions=dm.actuator_commands+nactuator;

	return dm;
}

double *DeformableMirror2array( DeformableMirror dm )
{
	int i, nactuator;
	double *tmp,*dm_ap; 

	nactuator=(dm.actuator_x_dim)*(dm.actuator_y_dim);

	tmp=(double *) calloc( 2*nactuator+22, sizeof(double) );

	tmp[0]=dm.actuator_x_dim; 
	tmp[1]=dm.actuator_y_dim;

	dm_ap=APERTURE2array(dm.AP);
	for(i=0;i<17;i++) tmp[2+i]=dm_ap[i];
	free(dm_ap);

	tmp[19]=dm.actuator_pitch;
	tmp[20]=dm.actuator_velocity;
	tmp[21]=dm.timestamp;

	for( i=0;i<nactuator;i++) {
		tmp[22+i]=dm.actuator_commands[i];
		tmp[22+nactuator+i]=dm.actuator_positions[i];
	}

	return tmp;
}

