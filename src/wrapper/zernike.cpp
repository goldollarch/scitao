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

#include "arroyo_wrap.h"
using namespace Arroyo;

int get_cos_index(int order, int level)
{
	int index = 0;
	for(int i=0; i<order; i++)
		index += i+1;
	int tmp_level = order;
	while(tmp_level > level)
	{
		if(tmp_level==0) index+=1;
		else index+=2;
		tmp_level-=2;
	}
	return(index);
}

int get_sin_index(int order, int level)
{
	int index = get_cos_index(order, level)+1;
	return(index);
}

int total_znk_space( int order )
{
	int nelem = 0;
	for(int i=0; i<=order; i++)
	{
		if(i%2) nelem += i/2+i/2+2; 
		else nelem += i/2+i/2+1;
	}
	return nelem;
}

ZERNIKE create_ZERNIKE(int order)
{
	int nelem;
	ZERNIKE ZNK;

	ZNK.order=order;
	nelem = total_znk_space( order );
	ZNK.coefficients=(double *)	calloc( nelem, sizeof(double) );

	return(ZNK);
}

ZERNIKE cpp2c_zernike(zernike znke)
{
	int order=znke.get_order();

	ZERNIKE ZNK=create_ZERNIKE(order);

	int index;
	for(int i=0; i<=order; i++){
		for(int j=i%2; j<=i; j+=2){
			index=get_cos_index(i,j);
			ZNK.coefficients[index]=znke.get_cos_coeff(i,j);
			if(j!=0){
				index=get_sin_index(i,j);
				ZNK.coefficients[index]=znke.get_sin_coeff(i,j);
			}
		}
	}
	return(ZNK);
}

zernike c2cpp_zernike(ZERNIKE ZNK)
{
	zernike znke(ZNK.order);

	int index;
	for(int i=0; i<=ZNK.order; i++){
		for(int j=i%2; j<=i; j+=2){
			index=get_cos_index(i,j);
			znke.set_cos_coeff(i,j,ZNK.coefficients[index]);
			if(j!=0){
				index=get_sin_index(i,j);
				znke.set_sin_coeff(i,j,ZNK.coefficients[index]);
			}
		}
	}

	return(znke);
}

void set_ZNK_cos_coef(ZERNIKE *ZNK, int order, int level, double coeff)
{
	int index;
	index=get_cos_index(order, level);
	ZNK->coefficients[index]=coeff;
}

void set_ZNK_sin_coef(ZERNIKE *ZNK, int order, int level, double coeff)
{
	int index;
	index=get_sin_index(order, level);
	ZNK->coefficients[index]=coeff;
}

double get_ZNK_cos_coeff (ZERNIKE ZNK,int order, int level) 
{
	int index;
	index=get_cos_index(order, level);
	return ZNK.coefficients[index];
}

double get_ZNK_sin_coeff (ZERNIKE ZNK,int order, int level)
{
	int index;
	index=get_sin_index(order, level);
	return ZNK.coefficients[index];
}

void zernike_write(ZERNIKE ZNK,char *fname,double timestamp )
{
	zernike znke=c2cpp_zernike(ZNK);
   if( timestamp <0 )
	   znke.write(simple_filename(fname,".fits"));
   else
	   znke.write(Current_filename(fname,timestamp,".fits"));
}

void ZERNIKE_array(ZERNIKE ZNK,double *outptr)
{
	int i, nelem;
	nelem=total_znk_space( ZNK.order );
	outptr[0]=ZNK.order;
	for(i=0;i<nelem;i++) outptr[i+1]=ZNK.coefficients[i];
}

double *ZERNIKE2array( ZERNIKE ZNK )
{
	int nelem;
	double *outptr;

	nelem=total_znk_space( ZNK.order );
	outptr= (double *)	calloc( nelem+1, sizeof(double) );
	ZERNIKE_array(ZNK,outptr);

	return outptr;
}

ZERNIKE array2ZERNIKE( double *inptr )
{
	ZERNIKE ZNK=create_ZERNIKE( (int)inptr[0] );
	ZNK.coefficients=inptr+1;
	return ZNK;
}

