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

void mixer(scicos_block *block,int flag)
{
	int i,number,number2,nport; 
	long k,nelem;
	int *ipar;

	ipar=block->ipar; 
	--ipar;

	nport=ipar[1];
	
	if(flag==1) 
	{
		number = (int)block->inptr[0][0];
		number2 = (int)block->inptr[0][10];

		if( number && number2 ) 
		{
			nelem=30+2*number*number2;
			if (block->nin == 1)
				for(k=0;k<nelem;k++)
					for(i=0;i<nport;++i)	
						block->outptr[i][k] = block->inptr[0][k]; 
			else 
			{
				for(k=0;k<nelem;k++) 
					block->outptr[0][k] = block->inptr[0][k];
				for(k=30;k<nelem;k++) 
					for(i=1;i<nport;i++)
						block->outptr[0][k] += block->inptr[i][k]; 
			}
		}
	}

}
