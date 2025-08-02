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

Optic create_Optic(int foreshortening,ThreeFrame TF)
{
	Optic Op;
	Op.foreshortening=foreshortening;
	Op.TF=TF;
	return Op;
}

PixelArray create_PixelArray(int x_axes, int y_axes,int wt_alloc)
{
	PixelArray tmp_PixArr;
	tmp_PixArr.x_axes=x_axes; tmp_PixArr.y_axes=y_axes;
	tmp_PixArr.wt_alloc=wt_alloc;
	
	long nelem=tmp_PixArr.x_axes*tmp_PixArr.y_axes;
	tmp_PixArr.pixeldata=(double*)malloc(nelem*sizeof(double));
	for(long i=0;i<nelem;i++) tmp_PixArr.pixeldata[i]=0;

	if( wt_alloc ) tmp_PixArr.pixelwts=(float*)malloc(nelem*sizeof(float));
	return tmp_PixArr;
}

void set_PixelArray_data(PixelArray *PixArr,double *data)
{
	int index;
	for(int i=0;i<PixArr->x_axes;i++)
		for(int j=0;j<PixArr->y_axes;j++)
		{
			index=i*PixArr->y_axes+j;
			PixArr->pixeldata[index]=data[index];
		}
}

PixelArray array2PixelArray( double *inptr)
{
	int wt_alloc=0, x_axes, y_axes;
	PixelArray tmp_PixArr;

	x_axes=(int)inptr[0]; y_axes=(int)inptr[1];

	tmp_PixArr=create_PixelArray(x_axes,y_axes,wt_alloc);
	tmp_PixArr.pixeldata=inptr+2;

	return tmp_PixArr;
}

void PixelArray_array(PixelArray PixArr,double *outptr)
{
	int i;
	long nelem=PixArr.x_axes*PixArr.y_axes;

	outptr[0]=PixArr.x_axes; outptr[1]=PixArr.y_axes;
	for(i=0;i<nelem;i++)  
		outptr[i+2]=PixArr.pixeldata[i];
}

double *PixelArray2array ( PixelArray PixArr )
{
	double *outptr;

	long nelem=PixArr.x_axes*PixArr.y_axes;
	outptr=(double*)calloc( nelem+2, sizeof(double) );

	PixelArray_array(PixArr,outptr);
	return outptr;
}

PixelArray cpp2c_PixelArray(pixel_array<double> pixarr)
{
	int axes_x,axes_y;
	axes_x=pixarr.get_axes()[0]; axes_y=pixarr.get_axes()[1];
	long nelem=axes_x*axes_y;

	PixelArray tmp_PixArr;
    tmp_PixArr=create_PixelArray(axes_x,axes_y,pixarr.weights_allocated());
	
	for(int i=0;i<nelem;i++) {
		tmp_PixArr.pixeldata[i]=pixarr.data(i);
		if( pixarr.weights_allocated())	
			tmp_PixArr.pixelwts[i]=pixarr.wt(i);
	}
	return tmp_PixArr;
}

pixel_array<double> c2cpp_PixelArray(PixelArray PixArr)
{
	vector<long> pixarr_axes(2);
	pixarr_axes[0] = PixArr.x_axes;  pixarr_axes[1] = PixArr.y_axes;
	pixel_array<double> tmp_pixarr(pixarr_axes,PixArr.pixeldata,NULL);
	//now wt_alloc must be 0 ie. false;

	return tmp_pixarr;
}

char *simple_filename(char *fname,char *type)
{
   char *tmp;
   tmp=(char*)malloc(100*sizeof(char));

   stringstream filename_stream;
   filename_stream.str("");
   filename_stream << fname<< type<<"\0";
   strcpy(tmp, filename_stream.str().c_str());

   return(tmp);
}

char *number_filename(char *fname,int n, char *type)
{
   char *tmp;
   tmp=(char*)malloc(100*sizeof(char));

   stringstream filename_stream;
   filename_stream.str("");
   filename_stream << fname<<"_"<<n<< type<<"\0";
   strcpy(tmp, filename_stream.str().c_str());

   return(tmp);
}

char* Current_filename(char *fname,double timestamp,char *type)
{
   char *tmp;
   tmp=(char*)malloc(100*sizeof(char));

   stringstream filename_stream;
   filename_stream.setf(ios::fixed, 
		       ios::floatfield);
   filename_stream.str("");
   filename_stream << fname<<"_"
	        << setfill('0')
			<< setprecision(4) 
		    << setw(4) 
		    << timestamp
		    << type<<"\0";
   strcpy(tmp, filename_stream.str().c_str());

   return(tmp);
}

char* current_filename(char *fname,double timestamp,char *type)
{
   return(Current_filename(fname,timestamp,type));
}
