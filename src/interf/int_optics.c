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

#include "mex.h" 

#include "int_arroyo.h"
#include "int_lightPipes.h"
#include "int_others.h"

static int direct_gateway(char *fname,void F(void)) 
{ 
	F();
	return 0; 
};
 
int C2F(int_lightPipes)()
{
  Rhs = Max(0, Rhs);
  (*(lightPipes_Tab[Fin-1].f))(lightPipes_Tab[Fin-1].name,lightPipes_Tab[Fin-1].F);
  return 0;
}
 
int C2F(int_arroyo)()
{
  Rhs = Max(0, Rhs);
  (*(arroyo_Tab[Fin-1].f))(arroyo_Tab[Fin-1].name,arroyo_Tab[Fin-1].F);
  return 0;
}

int C2F(int_others)()
{
  Rhs = Max(0, Rhs);
  (*(others_Tab[Fin-1].f))(others_Tab[Fin-1].name,others_Tab[Fin-1].F);
  return 0;
}
