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

void emitter(scicos_block *block,int flag) 
{
    int type,nx_grid,ny_grid,beacon;
	double x,y,z,lambda,pixscale;
	
	double *rpar;
	int *ipar;

	FIELD field;
	ThreeFrame TF;

	WavefrontHeader WfH;
	WaveFront WF;

	static Emitter Emt;
	ThreeVector TV;

	rpar=block->rpar; 
	--rpar;
	ipar=block->ipar; 
	--ipar;

	type=ipar[1];
	beacon=ipar[2];
	nx_grid=ipar[3]; 
	ny_grid=ipar[4];

	x=rpar[1]; y=rpar[2]; z=rpar[3];	
	lambda=rpar[4]; pixscale=rpar[5];

	if (flag==1) {

		field = create_field( scao.max_number, scao.max_number );

		field.size=nx_grid*pixscale;
		field.int1=0; field.int2=0; field.int3=0;
		field.double1=0; field.double2=0; field.double3=0;

		if(scao.atm_mod_ref_lay_created!=1) {
			TF=default_ThreeFrame();
			WfH = construct_WavefrontHeader( TF, 
			   nx_grid, ny_grid, lambda, pixscale );
		}
		else if(beacon) 
			WfH = scao.sensing_WfH;
		else 
			WfH = scao.detected_WfH;

		WfH.timestamp=get_scicos_time();

		WF = wave_emit( Emt,WfH );

		WaveFront_FIELD( WF,&field );

		field_block_outptr(field,block,0);

	}

	else if (flag==4) {
		scao.max_number=nx_grid;   
		if(scao.max_number<ny_grid) scao.max_number=ny_grid;
	}

	else if (flag==6 ) {	
		
		Emt = construct_Emitter( type,x,y,z );	
		
		if(scao.global_AP_created && type==0) {
			TV = construct_ThreeVector(-scao.global_AP.TF.TVz.x,
				-scao.global_AP.TF.TVz.y,-scao.global_AP.TF.TVz.z);
			Emt.plane_wave_emission_direction = TV;
		}

		if(beacon && !scao.sensing_emt_created) {
			scao.sensing_emt = Emt;
			scao.sensing_wavelength = lambda;
			scao.sensing_pixel_scale = pixscale;
			scao.sensing_emt_created=1;
		}
		else if (!scao.detected_emt_created) {
			scao.detected_emt = Emt;
			scao.detected_wavelengths = lambda;
			scao.detected_pixel_scale = pixscale;
			scao.detected_emt_created=1;
		}
	} 
}
