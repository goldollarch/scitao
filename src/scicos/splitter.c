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

void splitter(scicos_block *block,int flag)
{
	double *rpar;
	double split;
	
	FIELD field,*tmp;

	rpar=block->rpar; 
	--rpar;

	/*   rpar[1]:  split       */

	split=rpar[1];
	
	if(flag==1)
	{
		field=block_inptr_field( block,0 );

		tmp=b_split( field, split );

		field_block_outptr( tmp[0],block,0 );
		field_block_outptr( tmp[1],block,1 );

	}

}
