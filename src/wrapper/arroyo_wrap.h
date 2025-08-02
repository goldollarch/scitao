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

#ifndef ARROYO_WRAP_H
#define ARROYO_WRAP_H

#include <string.h>

#include "../lightPipes/pipes.h"

#include "zernike.h"
#include "geometry.h"
#include "some_utils.h"

#include "wavefront.h"
#include "emitter.h"
#include "aperture.h"
#include "propagation.h"

#include "shack_hartmann.h"
#include "refractive_layer.h"
#include "atmospheric_model.h"
#include "tip_tilt_mirror.h"
#include "deformable_mirror.h"
#include "controller.h"
#include "reconstructor.h"
#include "wind_model.h"

#include "simulation.h"

#endif
