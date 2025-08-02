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

three_point c2cpp_ThreePoint(ThreePoint TP)
{
   three_frame tf;
   three_point tp(TP.x,TP.y,TP.z,tf);
   return tp;
}

three_vector c2cpp_ThreeVector(ThreeVector TV)
{
   three_frame tf; 
   three_vector tv(TV.x,TV.y,TV.z,tf); 
   return tv;
}

three_frame c2cpp_ThreeFrame(ThreeFrame TF)
{
	three_point tp=c2cpp_ThreePoint(TF.TP);
	three_vector tv_x=c2cpp_ThreeVector(TF.TVx); 
	three_vector tv_y=c2cpp_ThreeVector(TF.TVy); 
	three_vector tv_z=c2cpp_ThreeVector(TF.TVz); 
	three_frame tf(tp,tv_x,tv_y,tv_z); 
	return(tf);
}

ThreePoint cpp2c_ThreePoint(three_point tp)
{
	ThreePoint TP;
	three_frame tf;
	TP.x=tp.x(tf); TP.y=tp.y(tf); TP.z=tp.z(tf);
	return TP;
}

ThreeVector cpp2c_ThreeVector(three_vector tv)
{
	ThreeVector TV;
	three_frame tf;
	TV.x=tv.x(tf);	TV.y=tv.y(tf);	TV.z=tv.z(tf);
	return TV;
}

ThreeFrame cpp2c_ThreeFrame(three_frame tf)
{
	ThreeFrame TF;
	TF.TP=cpp2c_ThreePoint(tf);
	TF.TVx=cpp2c_ThreeVector(tf.x());
	TF.TVy=cpp2c_ThreeVector(tf.y());
	TF.TVz=cpp2c_ThreeVector(tf.z());
	return TF;
}

three_point field_three_point(FIELD field) 
{
   three_frame tf;
   three_point tp(field.x,field.y,field.z,tf);
   return tp;
}

void three_point_field(three_point tp, FIELD *field)
{
   three_frame tf;
   field->x=tp.x(tf);
   field->y=tp.y(tf);
   field->z=tp.z(tf);
}

three_vector field_x_vector(FIELD field) {
   three_frame tf; 
   three_vector tv(field.xx,field.xy,field.xz,tf); 
   return tv;
}

void vector_x_field(three_vector tv, FIELD *field)
{
   three_frame tf;
   field->xx=tv.x(tf);
   field->xy=tv.y(tf);
   field->xz=tv.z(tf);
}

three_vector field_y_vector(FIELD field) {
   three_frame tf; 
   three_vector tv(field.yx,field.yy,field.yz,tf); 
   return tv;
}

void vector_y_field(three_vector tv, FIELD *field)
{
   three_frame tf;
   field->yx=tv.x(tf);
   field->yy=tv.y(tf);
   field->yz=tv.z(tf);
}

three_vector field_z_vector(FIELD field) {
   three_frame tf; 
   three_vector tv(field.zx,field.zy,field.zz,tf); 
   return tv;
}

void vector_z_field(three_vector tv, FIELD *field)
{
   three_frame tf;
   field->zx=tv.x(tf);
   field->zy=tv.y(tf);
   field->zz=tv.z(tf);
}

three_frame  FIELD_three_frame(FIELD field)
{
   three_point  tp=field_three_point(field);
   three_vector tv_x=field_x_vector(field); 
   three_vector tv_y=field_y_vector(field); 
   three_vector tv_z=field_z_vector(field); 

   three_frame tf(tp,tv_x,tv_y,tv_z); 
   return(tf);
}

void three_frame_FIELD(three_frame tf, FIELD *field)
{
	three_point_field(tf,field);
	vector_x_field(tf.x(),field);
	vector_y_field(tf.y(),field);
	vector_z_field(tf.z(),field);
}

/////////////////////////////////////////////////////////////////////////

ThreePoint construct_ThreePoint(double x,double y,double z)
{
	ThreePoint TP;
	TP.x=x;TP.y=y;TP.z=z;
	return TP;
}

ThreeVector construct_ThreeVector(double x,double y,double z)
{
	ThreeVector TV;
	TV.x=x;TV.y=y;TV.z=z;
	return TV;
}

ThreeFrame default_ThreeFrame()
{
	ThreeFrame TF;

	TF.TP=construct_ThreePoint(0,0,0);
	TF.TVx=construct_ThreeVector(1,0,0);	
	TF.TVy=construct_ThreeVector(0,1,0);	
	TF.TVz=construct_ThreeVector(0,0,1);

	return TF;
}

ThreeFrame construct_ThreeFrame( ThreePoint TP,
		ThreeVector TVx, ThreeVector TVy, ThreeVector TVz)
{
	ThreeFrame TF;
	TF.TP=TP;
	TF.TVx=TVx;	TF.TVy=TVy;	TF.TVz=TVz;
	return TF;
}

void set_FIELD_ThreeFrame( ThreeFrame TF, FIELD *field )
{
	field->x=TF.TP.x;field->y=TF.TP.y;field->z=TF.TP.z;

	field->xx=TF.TVx.x;field->xy=TF.TVx.y;field->xz=TF.TVx.z;
	field->yx=TF.TVy.x;field->yy=TF.TVy.y;field->yz=TF.TVy.z;
	field->zx=TF.TVz.x;field->zy=TF.TVz.y;field->zz=TF.TVz.z;
}

double *ThreeFrame2array( ThreeFrame TF )
{
	double *tmp=NULL;

	tmp=(double*) calloc( 12,sizeof(double));

	tmp[0]=TF.TP.x;tmp[1]=TF.TP.y;tmp[2]=TF.TP.z;
	tmp[3]=TF.TVx.x;tmp[4]=TF.TVx.y;tmp[5]=TF.TVx.z;
	tmp[6]=TF.TVy.x;tmp[7]=TF.TVy.y;tmp[8]=TF.TVy.z;
	tmp[9]=TF.TVz.x;tmp[10]=TF.TVz.y;tmp[11]=TF.TVz.z;

	return tmp;
}

ThreeFrame array2ThreeFrame( double *inptr )
{
	ThreePoint TP=construct_ThreePoint(inptr[0],inptr[1],inptr[2]);
	ThreeVector TVx=construct_ThreeVector(inptr[3],inptr[4],inptr[5]);
	ThreeVector TVy=construct_ThreeVector(inptr[6],inptr[7],inptr[8]);
	ThreeVector TVz=construct_ThreeVector(inptr[9],inptr[10],inptr[11]);

	return construct_ThreeFrame( TP,TVx,TVy,TVz);
}

ThreeFrame FIELD_ThreeFrame( FIELD field )
{
	ThreePoint TP=construct_ThreePoint(field.x,field.y,field.z);
	ThreeVector TVx=construct_ThreeVector(field.xx,field.xy,field.xz);
	ThreeVector TVy=construct_ThreeVector(field.yx,field.yy,field.yz);
	ThreeVector TVz=construct_ThreeVector(field.zx,field.zy,field.zz);

	return construct_ThreeFrame( TP,TVx,TVy,TVz);
}

ThreeFrame ThreeTranslaton ( ThreeFrame TF, ThreeVector TV )
{
	three_vector tv=c2cpp_ThreeVector(TV);
	three_frame tf=c2cpp_ThreeFrame(TF);

	three_translation ttrans(tv);
	ttrans.transform(tf);

	return(cpp2c_ThreeFrame(tf));
}

ThreeFrame ThreeReflection ( ThreeFrame TF, ThreePoint TP, ThreeVector TV )
{
	three_point tp=c2cpp_ThreePoint(TP);
	three_vector tv=c2cpp_ThreeVector(TV);
	three_frame tf=c2cpp_ThreeFrame(TF);

	three_reflection tref(tp, tv);
	tref.transform(tf);

	return(cpp2c_ThreeFrame(tf));
}

ThreeFrame ThreeRotation ( ThreeFrame TF, ThreePoint TP, ThreeVector TV, double angle )
{
	three_point tp=c2cpp_ThreePoint(TP);
	three_vector tv=c2cpp_ThreeVector(TV);
	three_frame tf=c2cpp_ThreeFrame(TF);

	three_rotation trot(tp, tv, angle);
	trot.transform(tf);

	return(cpp2c_ThreeFrame(tf));
}
