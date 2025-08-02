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

diffractive_wavefront_header<double> c2cpp_WavefrontHeader(WavefrontHeader WfH)
{
   double wavelength, wf_pixel_scale;
   wavelength = WfH.lambda;wf_pixel_scale =WfH.pixscale;

   std::vector<long> wf_axes(2);
   wf_axes[0]=WfH.axes_x; wf_axes[1]=WfH.axes_y;

   three_frame wf_frame;
   wf_frame=c2cpp_ThreeFrame(WfH.TF);

   diffractive_wavefront_header<double> 
	   dwfh(wf_axes,wf_frame,wavelength,wf_pixel_scale);

   dwfh.set_curvature(WfH.curvature);
   dwfh.set_timestamp(WfH.timestamp);

   return(dwfh);
} 

WavefrontHeader
cpp2c_WavefrontHeader(diffractive_wavefront_header<double> dwfh)
{
	WavefrontHeader WfH;

	WfH.TF=cpp2c_ThreeFrame(dwfh);

	WfH.axes_x=dwfh.get_axes()[0];
	WfH.axes_y=dwfh.get_axes()[1];

	WfH.lambda=dwfh.get_wavelength();
	WfH.pixscale=dwfh.get_pixel_scale();

	WfH.curvature=dwfh.get_curvature();
	WfH.timestamp=dwfh.get_timestamp();

	return WfH;
}

WaveFront cpp2c_WaveFront(diffractive_wavefront<double> dwf)
{
	WaveFront tmp_WF;
	tmp_WF.WfH=cpp2c_WavefrontHeader(dwf);
	tmp_WF.real_imag = 0; tmp_WF.interleaved = 1;

	pixel_amp_array<double> pixamparr(dwf.extract_amps());
	pixel_phase_array<double> pixpharr(dwf.extract_phases());

	long nelem = dwf.get_axes()[0]*dwf.get_axes()[1];
	tmp_WF.real = (double *) calloc( nelem, sizeof(double) );
	tmp_WF.imaginary = (double *) calloc( nelem, sizeof(double) );

	for(int i=0;i<nelem;i++) 
	{
		tmp_WF.real[i]=pixamparr.data(i);
		tmp_WF.imaginary[i]=pixpharr.data(i);
	}

	return tmp_WF;
}

diffractive_wavefront<double> c2cpp_WaveFront( WaveFront WF)
{
   diffractive_wavefront_header<double> dwfh=c2cpp_WavefrontHeader(WF.WfH);
   diffractive_wavefront<double> dwf(dwfh,NULL,WF.real_imag,WF.interleaved);

   if( WF.real_imag ){  
	   std::complex<double> c;
	   for(int i=0;i<WF.WfH.axes_x;i++) {
		   for(int j=0;j<WF.WfH.axes_y;j++) {
			   c = std::complex<double>
				   (WF.real[i*WF.WfH.axes_y+j],WF.imaginary[i*WF.WfH.axes_y+j]);
			   dwf.set_data(i*WF.WfH.axes_y+j,c);
		   }
	   }
   }
   else {
	   vector<long> pixarr_axes(dwfh.get_axes());
	   pixel_amp_array<double> pixamparr(pixarr_axes,WF.real,NULL);
	   pixel_phase_array<double> pixpharr(pixarr_axes,WF.imaginary,NULL);
	   dwf.install(pixamparr);dwf.install(pixpharr);
   }

   return(dwf);
}

//////////////////////////////////////////////////////////////////////////////////

WavefrontHeader construct_WavefrontHeader( ThreeFrame TF,
		int axes_x,int axes_y,double lambda, double pixscale)
{
	WavefrontHeader WfH;
    
	WfH.TF=TF;WfH.curvature=0;
	WfH.axes_x=axes_x;WfH.axes_y=axes_y;
	WfH.lambda=lambda; WfH.pixscale=pixscale; 

	return WfH;
}

void set_FIELD_WavefrontHeader(unsigned int number, unsigned int number2,
	  double lambda, double pixscale, double curvature,
	  double timestamp, FIELD *field)
{
    field->number=number;
    field->number2=number2;

    field->lambda=lambda;
    field->pixscale=pixscale;
    field->curvature=curvature;
    field->timestamp=timestamp;
}

void write_WaveFront_file(WaveFront WF,char *fname,double timestamp)
{
   diffractive_wavefront<double> dwf=c2cpp_WaveFront(WF);
//   dwf.set_timestamp(timestamp);

   if( timestamp <0 )
	   dwf.write(simple_filename(fname,".fits"));
   else
	   dwf.write(Current_filename(fname,timestamp,".fits"));
}

WavefrontHeader FIELD_WavefrontHeader( FIELD field)
{
	WavefrontHeader WfH;

	WfH.TF=FIELD_ThreeFrame( field );

	WfH.axes_x=field.number;
	WfH.axes_y=field.number2;

	WfH.lambda=field.lambda;
	WfH.pixscale=field.pixscale;
	WfH.curvature=field.curvature;
	WfH.timestamp=field.timestamp;

	return WfH;
}

WaveFront FIELD_WaveFront( FIELD field)
{
	WaveFront WF;

	WF.WfH = FIELD_WavefrontHeader( field );
	WF.interleaved = field.interleaved;
	WF.real_imag = field.real_imag;
	WF.real = field.real;
	WF.imaginary = field.imaginary;

	return WF;
}

void WaveFront_FIELD( WaveFront WF, FIELD *field )
{
	set_FIELD_ThreeFrame( WF.WfH.TF,field );

	field->number = WF.WfH.axes_x;
	field->number2 = WF.WfH.axes_y;

	field->lambda = WF.WfH.lambda;
	field->pixscale = WF.WfH.pixscale;
	field->curvature = WF.WfH.curvature;
	field->timestamp = WF.WfH.timestamp;

	field->interleaved=WF.interleaved;
	field->real_imag=WF.real_imag;

	field->size=WF.WfH.axes_x*WF.WfH.pixscale;
	field->int1=0;field->int2=0;field->int3=0;
	field->double1=0;field->double2=0;field->double3=0;

	field->real = WF.real;
	field->imaginary = WF.imaginary;
}

double *WavefrontHeader2array( WavefrontHeader WFH )
{
	double *wfh, *tmp;
	wfh=(double*) calloc(18,sizeof(double));

	wfh[0]=WFH.axes_x;wfh[1]=WFH.axes_y;
	wfh[2]=WFH.lambda;wfh[3]=WFH.pixscale;
	wfh[4]=WFH.curvature;wfh[5]=WFH.timestamp;

	tmp=ThreeFrame2array(WFH.TF);
	for(int i=0;i<12;i++) wfh[i+6]=tmp[i];

	return wfh;
}

WavefrontHeader array2WavefrontHeader( double *inptr )
{
	WavefrontHeader WfH;
	WfH.axes_x=(int)inptr[0];WfH.axes_y=(int)inptr[1];
	WfH.lambda=inptr[2]; WfH.pixscale=inptr[3]; 
	WfH.curvature=inptr[4]; WfH.timestamp=inptr[5];
	WfH.TF=array2ThreeFrame(inptr+6);
	return WfH;
}

WaveFront shift_WaveFront ( WaveFront WF, double dx, double dy )
{
	diffractive_wavefront<double> dwf=c2cpp_WaveFront(WF);

	three_vector tv(dx,dy,0,dwf);
	three_translation ttrans(tv);
	ttrans.transform(dwf);

	return cpp2c_WaveFront(dwf);
}

WaveFront rotate_WaveFront ( WaveFront WF, double angle)
{
	diffractive_wavefront<double> dwf=c2cpp_WaveFront(WF);

	three_rotation trot(dwf,dwf.z(),angle);
	trot.transform(dwf);

	return cpp2c_WaveFront(dwf);
}

WaveFront reflect_WaveFront ( WaveFront WF, double x, double y, double z )
{
	diffractive_wavefront<double> dwf=c2cpp_WaveFront(WF);

	if(x*x+y*y+z*z) {
		three_vector tv(x,y,z,dwf);
		three_reflection tref(dwf,tv);
		tref.transform(dwf);
	}

	return cpp2c_WaveFront(dwf);
}

WaveFront set_WaveFront_propagation_direction( WaveFront WF, ThreeVector TV )
{
	diffractive_wavefront<double> dwf=c2cpp_WaveFront(WF);
	three_vector tv=c2cpp_ThreeVector(TV);
	dwf.set_propagation_direction(tv);
	return cpp2c_WaveFront(dwf);
}

WaveFront WaveFront_clip_array( WaveFront WF, int nclip )
{
	diffractive_wavefront<double> dwf=c2cpp_WaveFront(WF);
	dwf.clip_array(nclip);
	return cpp2c_WaveFront(dwf);
}

WaveFront set_dwf_pixel_scale( WaveFront WF, double pxlscl)
{
	diffractive_wavefront<double> dwf=c2cpp_WaveFront(WF);
	dwf.set_pixel_scale(pxlscl);
	return cpp2c_WaveFront(dwf);
}
