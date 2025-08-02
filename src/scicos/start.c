/*
scitao - A toolbox based on Scilab/Scicos environment for the simulation of 
wave optics, especially for the simulation of adaptive optics .

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

void start(scicos_block *block,int flag) 
{
	FIELD field;
	ThreeFrame TF;

    int nx_grid, ny_grid;
	double size_grid, lambda;
	double pixscale, curvature, timestamp;
	
	int *ipar;	
	double *rpar;

	rpar=block->rpar;
	--rpar;
	ipar=block->ipar;
	--ipar;

	nx_grid=ipar[1]; 
	ny_grid=ipar[2];

	size_grid=rpar[1]; 
	lambda=rpar[2];	

	pixscale=size_grid/nx_grid;
	curvature=0;

	timestamp=get_scicos_time();

	if (flag==1) {

		TF = default_ThreeFrame( );

		field = lp_begin( size_grid, lambda, nx_grid, ny_grid );

		set_FIELD_ThreeFrame( TF, &field );

		set_FIELD_WavefrontHeader( nx_grid, nx_grid,
			lambda,	pixscale,curvature,timestamp, &field );

		field_block_outptr( field,block,0 ); 

	}

	else if (flag==4) {
		scao.detected_emt = construct_Emitter( 0,0,0,1 );
		scao.detected_wavelengths = lambda;
		scao.detected_pixel_scale = pixscale;

		scao.max_number=nx_grid;   
		if(scao.max_number<ny_grid)scao.max_number=ny_grid;
	} 
	
}
