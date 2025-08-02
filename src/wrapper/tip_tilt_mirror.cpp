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

TipTiltMirror create_TipTiltMirror( APERTURE AP, double angular_velocity )
{
	TipTiltMirror TTM;

	TTM.AP = AP;
	TTM.angular_velocity=angular_velocity;

	TTM.timestamp=0;
	
	ThreeVector TV=AP.TF.TVz;
    set_command_vector(&TTM,TV.x,TV.y,TV.z);

	return TTM;
}

void set_ttm_timestamp( TipTiltMirror* TTM, double timestamp )
{
	TTM->timestamp=timestamp;
}

void set_command_vector( TipTiltMirror* TTM,	 
	double x_orient, double y_orient, double z_orient)
{
	TTM->x_orient=x_orient;
	TTM->y_orient=y_orient;
	TTM->z_orient=z_orient;
}

WaveFront TipTiltMirror_transform( TipTiltMirror* TTM, WaveFront WF )
{
	int nspiders;

	double diameter,in_diameter,out_diameter,
		x_size,y_size,in_edge_length,spider_width,in_gap_size;

	diffractive_wavefront<double> dwf=c2cpp_WaveFront( WF );

	ThreeVector TV=construct_ThreeVector(TTM->x_orient,TTM->y_orient,TTM->z_orient);
	three_vector ttm_command_vector=c2cpp_ThreeVector(TV);

	Arroyo::three_frame oriented_frame=c2cpp_ThreeFrame(TTM->AP.TF);

	switch(TTM->AP.aperture_type)  {
		case 1: 
			{
				in_diameter=TTM->AP.ap_data2;
				out_diameter=TTM->AP.ap_data1;
				ideal_tip_tilt_mirror<annular_aperture> tip_tilt_mirror
					(annular_aperture(in_diameter,out_diameter),
					TTM->angular_velocity,TTM->timestamp);
				tip_tilt_mirror.Arroyo::three_frame::operator=(oriented_frame);
				tip_tilt_mirror.set_mirror_command_vector(ttm_command_vector);
				tip_tilt_mirror.transform(dwf);
				break;
			}
		case 2:
			{
				x_size=TTM->AP.ap_data1;
				y_size=TTM->AP.ap_data2;
				ideal_tip_tilt_mirror<rectangular_aperture> tip_tilt_mirror
					(rectangular_aperture(x_size,y_size),
					TTM->angular_velocity,TTM->timestamp);
				tip_tilt_mirror.Arroyo::three_frame::operator=(oriented_frame);
				tip_tilt_mirror.set_mirror_command_vector(ttm_command_vector);
				tip_tilt_mirror.transform(dwf);
				break;
			}
		case 3:
			{
				in_edge_length=TTM->AP.ap_data1;
				ideal_tip_tilt_mirror<hexagonal_aperture> tip_tilt_mirror
					(hexagonal_aperture(in_edge_length),
					TTM->angular_velocity,TTM->timestamp);
				tip_tilt_mirror.Arroyo::three_frame::operator=(oriented_frame);
				tip_tilt_mirror.set_mirror_command_vector(ttm_command_vector);
				tip_tilt_mirror.transform(dwf);
				break;
			}
		case 4:
			{
				in_diameter=TTM->AP.ap_data2;
				out_diameter=TTM->AP.ap_data1;
				nspiders=(int)(TTM->AP.ap_data3); 
				spider_width=TTM->AP.ap_data4;
				ideal_tip_tilt_mirror<spidered_annular_aperture> 
					tip_tilt_mirror(spidered_annular_aperture
					(in_diameter,out_diameter,nspiders,spider_width),
					TTM->angular_velocity,TTM->timestamp);
				tip_tilt_mirror.Arroyo::three_frame::operator=(oriented_frame);
				tip_tilt_mirror.set_mirror_command_vector(ttm_command_vector);
				tip_tilt_mirror.transform(dwf);
				break;
			}

		default:
			{
				diameter=TTM->AP.ap_data1;
				ideal_tip_tilt_mirror<circular_aperture> tip_tilt_mirror
					(circular_aperture(diameter), 
					TTM->angular_velocity,TTM->timestamp);
				tip_tilt_mirror.Arroyo::three_frame::operator=(oriented_frame);
				tip_tilt_mirror.set_mirror_command_vector(ttm_command_vector);
				tip_tilt_mirror.transform(dwf);
				break;
			}
	}

	WaveFront tmp_WF=cpp2c_WaveFront(dwf);

	return tmp_WF;
}

void TipTiltMirror_updata( TipTiltMirror* TTM, 
		ThreeVector ttm_command_TV, double timestamp )
{
	int nspiders;

	ThreeFrame AP_TF;

	double diameter,in_diameter,out_diameter,
		x_size,y_size,in_edge_length,spider_width,in_gap_size;

	ThreeVector tmp_TV=construct_ThreeVector(TTM->x_orient,TTM->y_orient,TTM->z_orient);
	three_vector ttm_orientation=c2cpp_ThreeVector(tmp_TV);

	three_vector ttm_command_vector=c2cpp_ThreeVector(ttm_command_TV);
	Arroyo::three_frame aperture_frame=c2cpp_ThreeFrame(TTM->AP.TF);

	switch(TTM->AP.aperture_type)  {
		case 1: 
			{
				in_diameter=TTM->AP.ap_data2;
				out_diameter=TTM->AP.ap_data1;
				ideal_tip_tilt_mirror<annular_aperture> tip_tilt_mirror
					(annular_aperture(in_diameter,out_diameter),
					TTM->angular_velocity,TTM->timestamp);
				tip_tilt_mirror.Arroyo::three_frame::operator=(aperture_frame);
				tip_tilt_mirror.set_mirror_command_vector(ttm_orientation);

				tip_tilt_mirror.update(ttm_command_vector,timestamp);

				AP_TF=cpp2c_ThreeFrame(tip_tilt_mirror);
				ttm_command_vector=tip_tilt_mirror.get_orientation_command();
				tmp_TV=cpp2c_ThreeVector(ttm_command_vector);

				break;
			}
		case 2:
			{
				x_size=TTM->AP.ap_data1;
				y_size=TTM->AP.ap_data2;
				ideal_tip_tilt_mirror<rectangular_aperture> tip_tilt_mirror
					(rectangular_aperture(x_size,y_size),
					TTM->angular_velocity,TTM->timestamp);
				tip_tilt_mirror.Arroyo::three_frame::operator=(aperture_frame);
				tip_tilt_mirror.set_mirror_command_vector(ttm_orientation);

				tip_tilt_mirror.update(ttm_command_vector,timestamp);

				AP_TF=cpp2c_ThreeFrame(tip_tilt_mirror);
				ttm_command_vector=tip_tilt_mirror.get_orientation_command();
				tmp_TV=cpp2c_ThreeVector(ttm_command_vector);

				break;
			}
		case 3:
			{
				in_edge_length=TTM->AP.ap_data1;
				ideal_tip_tilt_mirror<hexagonal_aperture> tip_tilt_mirror
					(hexagonal_aperture(in_edge_length),
					TTM->angular_velocity,TTM->timestamp);
				tip_tilt_mirror.Arroyo::three_frame::operator=(aperture_frame);
				tip_tilt_mirror.set_mirror_command_vector(ttm_orientation);

				tip_tilt_mirror.update(ttm_command_vector,timestamp);

				AP_TF=cpp2c_ThreeFrame(tip_tilt_mirror);
				ttm_command_vector=tip_tilt_mirror.get_orientation_command();
				tmp_TV=cpp2c_ThreeVector(ttm_command_vector);

				break;
			}
		case 4:
			{
				in_diameter=TTM->AP.ap_data2;
				out_diameter=TTM->AP.ap_data1;
				nspiders=(int)(TTM->AP.ap_data3); 
				spider_width=TTM->AP.ap_data4;
				ideal_tip_tilt_mirror<spidered_annular_aperture> 
					tip_tilt_mirror(spidered_annular_aperture
					(in_diameter,out_diameter,nspiders,spider_width),
					TTM->angular_velocity,TTM->timestamp);
				tip_tilt_mirror.Arroyo::three_frame::operator=(aperture_frame);
				tip_tilt_mirror.set_mirror_command_vector(ttm_orientation);

				tip_tilt_mirror.update(ttm_command_vector,timestamp);

				AP_TF=cpp2c_ThreeFrame(tip_tilt_mirror);
				ttm_command_vector=tip_tilt_mirror.get_orientation_command();
				tmp_TV=cpp2c_ThreeVector(ttm_command_vector);

				break;
			}

		default:
			{
				diameter=TTM->AP.ap_data1;
				ideal_tip_tilt_mirror<circular_aperture> tip_tilt_mirror
					(circular_aperture(diameter), 
					TTM->angular_velocity,TTM->timestamp);
				tip_tilt_mirror.Arroyo::three_frame::operator=(aperture_frame);
				tip_tilt_mirror.set_mirror_command_vector(ttm_orientation);

				tip_tilt_mirror.update(ttm_command_vector,timestamp);

				AP_TF=cpp2c_ThreeFrame(tip_tilt_mirror);
				ttm_command_vector=tip_tilt_mirror.get_orientation_command();
				tmp_TV=cpp2c_ThreeVector(ttm_command_vector);

				break;
			}
	}

	TTM->AP.TF=AP_TF;
	TTM->timestamp=timestamp;
	TTM->x_orient=tmp_TV.x;	TTM->y_orient=tmp_TV.y;	TTM->z_orient=tmp_TV.z;

}

TipTiltMirror array2TipTiltMirror( double *inptr )
{
	TipTiltMirror ttm;

	ttm.AP = array2APERTURE(inptr);
	ttm.angular_velocity =inptr[17]; ttm.timestamp=inptr[18];
	ttm.x_orient=inptr[19];	ttm.y_orient=inptr[20];	ttm.z_orient=inptr[21];

	return ttm;
}

double *TipTiltMirror2array( TipTiltMirror ttm )
{
	int i; double *tmp,*ttm_ap; 

	tmp=(double *) calloc( 22, sizeof(double) );
	ttm_ap=APERTURE2array(ttm.AP);
	for(i=0;i<17;i++) tmp[i]=ttm_ap[i];
	free(ttm_ap);

	tmp[17]=ttm.angular_velocity; tmp[18]=ttm.timestamp;
	tmp[19]=ttm.x_orient;tmp[20]=ttm.y_orient;tmp[21]=ttm.z_orient;

	return tmp;
}