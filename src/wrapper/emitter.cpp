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

emitter *c2cpp_Emitter(Emitter Emt)
{
	emitter *emt=NULL;

	if(Emt.type==1) {
		three_point swemitter_tp
			=c2cpp_ThreePoint(Emt.spherical_wave_emission_point);
		emt=new spherical_wave_emitter(swemitter_tp);
	}
	else {
		three_vector emission_direction
			=c2cpp_ThreeVector(Emt.plane_wave_emission_direction);
		emt=new plane_wave_emitter(emission_direction);
	}

	return emt;
}

Emitter construct_Emitter(int type,double x,double y,double z)
{
	Emitter Emt;
	Emt.type=type;
	if(type==1) 
		Emt.spherical_wave_emission_point=construct_ThreePoint(x,y,z);
	else 
		Emt.plane_wave_emission_direction=construct_ThreeVector(x,y,z);
	return Emt;
}

WaveFront wave_emit( Emitter Emt, WavefrontHeader WfH )
{
	WaveFront WF;
	diffractive_wavefront<double> dwf;

	if(Emt.type!=1) WfH.curvature=0;
	diffractive_wavefront_header<double> dwfh=c2cpp_WavefrontHeader(WfH);

	emitter *emt=c2cpp_Emitter(Emt);
	dwf=emt->emit(dwfh);

	WF=cpp2c_WaveFront(dwf);

	return WF;
}

double *Emitter2array(Emitter emt)
{
	double *tmp;

	tmp = (double*) calloc ( 4,sizeof(double) );

	tmp[0]=emt.type;
	if (emt.type==1) {
		tmp[1]=emt.spherical_wave_emission_point.x;
		tmp[2]=emt.spherical_wave_emission_point.y;
		tmp[3]=emt.spherical_wave_emission_point.z;
	} else {
		tmp[1]=emt.plane_wave_emission_direction.x;
		tmp[2]=emt.plane_wave_emission_direction.y;
		tmp[3]=emt.plane_wave_emission_direction.z;
	}

	return tmp;
}