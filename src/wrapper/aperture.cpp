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

//////////////////////////////////////////////////////////////////////////////////

APERTURE create_APERTURE(int aperture_type,
		double ap_data1,double ap_data2,
		double ap_data3,double ap_data4 )
{
	APERTURE aperture;

	ThreeFrame TF;
	TF=default_ThreeFrame();
	aperture.TF=TF;

	aperture.aperture_type=aperture_type;
	aperture.ap_data1=ap_data1;	aperture.ap_data2=ap_data2;
	aperture.ap_data3=ap_data3;	aperture.ap_data4=ap_data4;

	return(aperture);
}

double *APERTURE2array ( APERTURE AP )
{
	int i;
	double *tmp,*tf;

	tmp=(double *) calloc( 17, sizeof(double) );

	tmp[0]=AP.aperture_type;
	tmp[1]=AP.ap_data1;	tmp[2]=AP.ap_data2;
	tmp[3]=AP.ap_data3;	tmp[4]=AP.ap_data4;

	tf=ThreeFrame2array(AP.TF);
	for(i=0;i<12;i++) tmp[5+i]=tf[i];
	free(tf);

	return tmp;
}

APERTURE array2APERTURE (double *inptr)
{
	APERTURE aperture;

	aperture.aperture_type=(int)inptr[0];

	aperture.ap_data1=inptr[1];	aperture.ap_data2=inptr[2];
	aperture.ap_data3=inptr[3];	aperture.ap_data4=inptr[4];

	aperture.TF = array2ThreeFrame(inptr+5);

	return(aperture);
}

aperture *c2cpp_APERTURE(APERTURE AP)
{
	int nspiders;

	double diameter,in_diameter,out_diameter,
		x_size,y_size,in_edge_length,spider_width,in_gap_size;

	aperture *ap = NULL;
	
	switch(AP.aperture_type) 
	{
		case 1:
			out_diameter=AP.ap_data1;
			in_diameter=AP.ap_data2;
			ap=new annular_aperture(in_diameter,out_diameter);
			break;
		case 2:
			x_size=AP.ap_data1;	
			y_size=AP.ap_data2;
			ap=new rectangular_aperture(x_size,y_size);
			break;
		case 3:
			in_edge_length=AP.ap_data1;
			ap=new hexagonal_aperture (in_edge_length);
			break;
		case 4:
			out_diameter=AP.ap_data1;
			in_diameter=AP.ap_data2;
			nspiders=(int)(AP.ap_data3); 
			spider_width=AP.ap_data4;
			ap=new spidered_annular_aperture 
				(in_diameter,out_diameter,nspiders,spider_width);
			break;
	    ///////////////////////////////////////////////////////////////
		case 5:
			in_diameter=AP.ap_data2;
			out_diameter=AP.ap_data1;
			in_edge_length=AP.ap_data3;	
			in_gap_size=AP.ap_data4;
			ap = new tiled_hexagonal_aperture 
				(in_diameter,out_diameter,in_edge_length,in_gap_size);
			break;	
		///////////////////////////////////////////////////////////////
		default:
			diameter=AP.ap_data1;
			ap=new circular_aperture(diameter);
			break;
	}

	return ap;

}

WaveFront APERTURE_transform(APERTURE AP,WaveFront WF)
{
	aperture *ap = c2cpp_APERTURE(AP);
	diffractive_wavefront<double> dwf=c2cpp_WaveFront( WF );

	ap->transform(dwf);
	WaveFront tmp_WF=cpp2c_WaveFront(dwf);
	return tmp_WF;
}
