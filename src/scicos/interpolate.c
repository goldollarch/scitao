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

void interpolate(scicos_block *block,int flag)
{
	int *ipar;
	double *rpar, *tmp;

    double new_size, magnif;
    int new_number;

	FIELD field;

	rpar=block->rpar; 
	--rpar;
	ipar=block->ipar; 
	--ipar;
  
	new_number=ipar[1];  

	new_size=rpar[1];
	magnif=rpar[2];

	if(flag==1)
	{
		field=block_inptr_field( block,0 );
		tmp=lp_interp1(field,new_size,new_number,0,0,0,magnif);

		field=array2field(tmp);
		field_block_outptr( field,block,0 );
	}

	else if (flag==4) {
		if( new_number > scao.max_number ) scao.max_number=new_number;
	}

}

