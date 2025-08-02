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

#ifndef WRAP_TIP_TILT_MIRROR_H
#define WRAP_TIP_TILT_MIRROR_H

#include "arroyo_wrap.h"

typedef struct{
	APERTURE AP;
	double timestamp,angular_velocity;
	double x_orient,y_orient,z_orient;
} TipTiltMirror;

#ifdef __cplusplus

#include "../arroyo/arroyo.h"

using namespace Arroyo;

extern "C" {
#endif

	TipTiltMirror create_TipTiltMirror( APERTURE AP, double angular_velocity );

	void set_ttm_timestamp( TipTiltMirror* ttm, double timestamp );
	void set_command_vector( TipTiltMirror* ttm,
		double x_orient, double y_orient, double z_orient );

	WaveFront TipTiltMirror_transform( TipTiltMirror* TTM, WaveFront WF );
	void TipTiltMirror_updata( TipTiltMirror* TTM, ThreeVector ttm_command_TV, double timestamp );

	TipTiltMirror array2TipTiltMirror( double *inptr );
	double *TipTiltMirror2array( TipTiltMirror TTM );

#ifdef __cplusplus
}
#endif

#endif
