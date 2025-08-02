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

/*
  introduces arbitrary Zernike aberration into the field distribution
  parameter:  n m R A
    where  n and m  are the integer orders 
      - see Born and Volf p. 465, sixth edition, Pergamon, 1993
    R is the radius at which the phase amplitude reaches A (in radians);
*/

#include "optics_scicos.h"

void zernike(scicos_block *block,int flag) 
{
	double *rpar;
	int *ipar;

    int n,m;
    double R,A;

 	FIELD field;

  /* reading the data from interface */
  /*   ipar[1-2]:  n,m       */
  /*   rpar[1-2]:  R,A       */

	rpar=block->rpar; 
	--rpar;
	ipar=block->ipar;
	--ipar;

	n=ipar[1];	
	m=ipar[2];

	R=rpar[1];	
	A=rpar[2];

	if(flag==1)
	{
		field = block_inptr_field( block,0 );

		field = lp_zernike( field,n,m,R,A );

		field_block_outptr( field,block,0 );
	}

}
