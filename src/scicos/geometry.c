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

void geometry(scicos_block *block,int flag) 
{
	double *rpar;
	double x,y,z;
	double dx,dy,angle;

	FIELD field;
	WaveFront WF_in,WF1,WF2,WF_out;

	rpar=block->rpar;
	--rpar;

	dx = rpar[1]; dy = rpar[2]; 
	angle = rpar[3];
	x = rpar[4]; y = rpar[5]; z = rpar[6]; 

	if (flag==1) {

		field = block_inptr_field( block,0 );
 
		if(field.number) {

			WF_in = FIELD_WaveFront(field);

			WF1 = shift_WaveFront(WF_in,dx,dy);

			WF2 = rotate_WaveFront(WF1,angle);

			WF_out = reflect_WaveFront(WF2,x,y,z);

			WaveFront_FIELD( WF_out,&field );
			
			field_block_outptr( field,block,0 ); 

		}

	}

}
