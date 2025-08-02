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

LensletArray create_LensletArray(int lenslet_x_axes,int lenslet_y_axes,
		double flength, double lnslt_pitch, int pix_per_lenslet, int pix_per_xform)
{
	LensletArray lnslt_array;
	lnslt_array.lenslet_x_axes=lenslet_x_axes;
	lnslt_array.lenslet_y_axes=lenslet_y_axes;
	lnslt_array.flength=flength;
	lnslt_array.lnslt_pitch=lnslt_pitch;
	lnslt_array.pix_per_lenslet=pix_per_lenslet;
	lnslt_array.pix_per_xform=pix_per_xform;
	return(lnslt_array);
}

LensletArray array2LensletArray( double *inptr)
{
	return create_LensletArray
		((int)inptr[0],(int)inptr[1],
		inptr[2], inptr[3],
		(int)inptr[4],(int)inptr[5]);
}

double *LensletArray2array ( LensletArray lnslt_array )
{
	double *tmp=(double*)calloc(6,sizeof(double));

	tmp[0]=lnslt_array.lenslet_x_axes;
	tmp[1]=lnslt_array.lenslet_y_axes;
	tmp[2]=lnslt_array.flength;
	tmp[3]=lnslt_array.lnslt_pitch;
	tmp[4]=lnslt_array.pix_per_lenslet;
	tmp[5]=lnslt_array.pix_per_xform;

	return tmp;
}

square_lenslet_array c2cpp_LensletArray(const LensletArray lnslt_array)
{
	std::vector<long> lenslet_axes(2);
	lenslet_axes[0]=lnslt_array.lenslet_x_axes; 
	lenslet_axes[1]=lnslt_array.lenslet_y_axes;

	square_lenslet_array sq_lnslt_array(lenslet_axes,
					     lnslt_array.flength,
					     lnslt_array.lnslt_pitch,
					     lnslt_array.pix_per_lenslet,
					     lnslt_array.pix_per_xform);

	return sq_lnslt_array;
}

WaveFront LensletArray_transform( LensletArray lnslt_array, WaveFront WF )
{
	square_lenslet_array sq_lnslt_array=c2cpp_LensletArray(lnslt_array);
	diffractive_wavefront<double> dwf=c2cpp_WaveFront( WF );
	sq_lnslt_array.transform(dwf);
	WaveFront tmp_WF=cpp2c_WaveFront(dwf);
	return tmp_WF;
}

Shack_Hartmann_centroids c2cpp_SHartmannCentroids(SHartmannCentroids shcentroids)
{
	std::vector<long> lenslet_axes(2);
	lenslet_axes[0]=shcentroids.pixarr_x_axes; 
	lenslet_axes[1]=(shcentroids.pixarr_y_axes)/2;

    long nlenslets = lenslet_axes[0]*lenslet_axes[1];

	Shack_Hartmann_centroids shack_hartmann_centroids(lenslet_axes);
	long index;
	for(int i=0; i<lenslet_axes[1]; i++) {
		for(int j=0; j<lenslet_axes[0]; j++) {
			index=i*lenslet_axes[0]+j;
			shack_hartmann_centroids.set_data(index,shcentroids.x_centroids[index]);
			shack_hartmann_centroids.set_data(index+nlenslets,shcentroids.y_centroids[index]);
		}
	}

	return(shack_hartmann_centroids);
}

SHartmannCentroids cpp2c_SHartmannCentroids(Shack_Hartmann_centroids shack_hartmann_centroids)
{ 
	SHartmannCentroids shcentroids;

	shcentroids.pixarr_x_axes=shack_hartmann_centroids.get_axes()[0];
	shcentroids.pixarr_y_axes=shack_hartmann_centroids.get_axes()[1];

	std::vector<long> lenslet_axes(2);
	lenslet_axes[0]=shack_hartmann_centroids.get_axes()[0]; 
	lenslet_axes[1]=(shack_hartmann_centroids.get_axes()[1])/2;

    long nlenslets = lenslet_axes[0]*lenslet_axes[1];

	shcentroids.x_centroids=(double *) calloc( nlenslets, sizeof(double) );
	shcentroids.y_centroids=(double *) calloc( nlenslets, sizeof(double) );

	long index;
	for(int i=0; i<lenslet_axes[1]; i++) {
		for(int j=0; j<lenslet_axes[0]; j++) {
			index=i*lenslet_axes[0]+j;
			shcentroids.x_centroids[index]=shack_hartmann_centroids.data(index);
			shcentroids.y_centroids[index]=shack_hartmann_centroids.data(index+nlenslets);
		}
	}
	
	return(shcentroids);
}

SHartmannCentroids create_SHartmannCentroids
(int lenslet_x_axes, int lenslet_y_axes, WaveFront WF )
{
	std::vector<long> lenslet_axes(2);
	lenslet_axes[0]=lenslet_x_axes; lenslet_axes[1]=lenslet_y_axes;

	diffractive_wavefront<double> dwf=c2cpp_WaveFront( WF );
	Shack_Hartmann_centroids shack_hartmann_centroids(lenslet_axes,dwf);
	return (cpp2c_SHartmannCentroids(shack_hartmann_centroids));
}

void SHartmannCentroids_array(SHartmannCentroids shcentroids,double *outptr)
{
	std::vector<long> lenslet_axes(2);
	lenslet_axes[0]=shcentroids.pixarr_x_axes; 
	lenslet_axes[1]=(shcentroids.pixarr_y_axes)/2;

    long nlenslets = lenslet_axes[0]*lenslet_axes[1];

	outptr[0]=shcentroids.pixarr_x_axes; 
	outptr[1]=shcentroids.pixarr_y_axes;

	long index;
	for(int i=0; i<lenslet_axes[1]; i++) {
		for(int j=0; j<lenslet_axes[0]; j++) {
			index=i*lenslet_axes[0]+j;
			outptr[index+2]=shcentroids.x_centroids[index];
			outptr[index+2+nlenslets]=shcentroids.y_centroids[index];
		}
	}
}

double *SHartmannCentroids2array( SHartmannCentroids shcentroids )
{
	double *outptr;

	std::vector<long> lenslet_axes(2);
	lenslet_axes[0]=shcentroids.pixarr_x_axes; 
	lenslet_axes[1]=(shcentroids.pixarr_y_axes)/2;

    long nlenslets = lenslet_axes[0]*lenslet_axes[1];
	outptr=(double *) calloc( 2*nlenslets+2, sizeof(double) );

	SHartmannCentroids_array(shcentroids,outptr);

	return outptr;
}

SHartmannCentroids array2SHartmannCentroids(double *inptr)
{
	SHartmannCentroids shcentroids;

	shcentroids.pixarr_x_axes=(int)inptr[0];
	shcentroids.pixarr_y_axes=(int)inptr[1];

	std::vector<long> lenslet_axes(2);
	lenslet_axes[0]=shcentroids.pixarr_x_axes; 
	lenslet_axes[1]=(shcentroids.pixarr_y_axes)/2;

    long nlenslets = lenslet_axes[0]*lenslet_axes[1];

	shcentroids.x_centroids=inptr+2;
	shcentroids.y_centroids=shcentroids.x_centroids+nlenslets;

	return(shcentroids);
}

void write_SHartmannCentroids(SHartmannCentroids shcentroids,char *fname,double timestamp)
{
	Shack_Hartmann_centroids shack_hartmann_centroids=c2cpp_SHartmannCentroids(shcentroids);

	if( timestamp <0 )
		shack_hartmann_centroids.write(simple_filename(fname,".fits"));
	else
		shack_hartmann_centroids.write(Current_filename(fname,timestamp,".fits"));
}

