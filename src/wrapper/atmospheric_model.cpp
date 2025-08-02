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

AtmosphericModel construct_AtmosphericModel( int nlayers,
		PowerSpectrum *power_spectra,double *layer_heights, ThreeFrame TF)
{
	AtmosphericModel AtmModel;

	AtmModel.type=0;

	AtmModel.nlayers=nlayers;

	AtmModel.Gemini_extended_profile=0;
	AtmModel.Gemini_focal_anisoplanatism_quality=0;
	AtmModel.Gemini_ground_layer_quality=0;
	AtmModel.Gemini_outer_scale_meters=0;
	AtmModel.Hufnagel_Valley_pseudowind=0;

	AtmModel.r_0_meters=0;
	AtmModel.r_0_ref_wavelength_meters=0;

	AtmModel.layer_heights=layer_heights;
	AtmModel.power_spectra=power_spectra;
	AtmModel.TF=TF;

	return AtmModel;
}

AtmosphericModel construct_Ellerbroek_Cerro_Pachon_model
		( ThreeFrame TF, double r_0_meters,double r_0_ref_wavelength_meters)
{
	AtmosphericModel PAtmModel;

	PAtmModel.type=1;

	PAtmModel.nlayers=0;
	PAtmModel.TF=TF;

	PAtmModel.Gemini_extended_profile=0;
	PAtmModel.Gemini_focal_anisoplanatism_quality=0;
	PAtmModel.Gemini_ground_layer_quality=0;
	PAtmModel.Gemini_outer_scale_meters=0;
	PAtmModel.Hufnagel_Valley_pseudowind=0;

	PAtmModel.r_0_meters=r_0_meters;
	PAtmModel.r_0_ref_wavelength_meters=r_0_ref_wavelength_meters;

	return PAtmModel;
}

AtmosphericModel construct_Ellerbroek_Mauna_Kea_model
		( ThreeFrame TF, double r_0_meters,double r_0_ref_wavelength_meters)
{
	AtmosphericModel PAtmModel;

	PAtmModel.type=2;

	PAtmModel.nlayers=0;
	PAtmModel.TF=TF;

	PAtmModel.Gemini_extended_profile=0;
	PAtmModel.Gemini_focal_anisoplanatism_quality=0;
	PAtmModel.Gemini_ground_layer_quality=0;
	PAtmModel.Gemini_outer_scale_meters=0;
	PAtmModel.Hufnagel_Valley_pseudowind=0;

	PAtmModel.r_0_meters=r_0_meters;
	PAtmModel.r_0_ref_wavelength_meters=r_0_ref_wavelength_meters;

	return PAtmModel;
}

AtmosphericModel construct_Palomar_DIMM_MASS_model
		( ThreeFrame TF, double r_0_meters,double r_0_ref_wavelength_meters)
{
	AtmosphericModel PAtmModel;

	PAtmModel.type=3;

	PAtmModel.nlayers=0;
	PAtmModel.TF=TF;

	PAtmModel.Gemini_extended_profile=0;
	PAtmModel.Gemini_focal_anisoplanatism_quality=0;
	PAtmModel.Gemini_ground_layer_quality=0;
	PAtmModel.Gemini_outer_scale_meters=0;
	PAtmModel.Hufnagel_Valley_pseudowind=0;

	PAtmModel.r_0_meters=r_0_meters;
	PAtmModel.r_0_ref_wavelength_meters=r_0_ref_wavelength_meters;

	return PAtmModel;
}

AtmosphericModel construct_Hufnagel_Valley_model( ThreeFrame TF, int nlayers,
					double *layer_heights, double Hufnagel_Valley_pseudowind )
{
	AtmosphericModel PAtmModel;

	PAtmModel.type=4;

	PAtmModel.TF=TF;
	PAtmModel.nlayers=nlayers;

	PAtmModel.Gemini_extended_profile=0;
	PAtmModel.Gemini_focal_anisoplanatism_quality=0;
	PAtmModel.Gemini_ground_layer_quality=0;
	PAtmModel.Gemini_outer_scale_meters=0;

	PAtmModel.r_0_meters=0;
	PAtmModel.r_0_ref_wavelength_meters=0;

	PAtmModel.layer_heights=(double*)malloc(nlayers*sizeof(double));
	for(int i=0;i<nlayers;i++)
		PAtmModel.layer_heights[i]=layer_heights[i];

	PAtmModel.Hufnagel_Valley_pseudowind=Hufnagel_Valley_pseudowind;

	return PAtmModel;
}

AtmosphericModel construct_SLCSAT_day_model
		( ThreeFrame TF, int nlayers, double *layer_heights )
{
	AtmosphericModel PAtmModel;
	PAtmModel.type=5;

	PAtmModel.TF=TF;
	PAtmModel.nlayers=nlayers;

	PAtmModel.Gemini_extended_profile=0;
	PAtmModel.Gemini_focal_anisoplanatism_quality=0;
	PAtmModel.Gemini_ground_layer_quality=0;
	PAtmModel.Gemini_outer_scale_meters=0;
	PAtmModel.Hufnagel_Valley_pseudowind=0;

	PAtmModel.r_0_meters=0;
	PAtmModel.r_0_ref_wavelength_meters=0;

	PAtmModel.layer_heights=(double*)malloc(nlayers*sizeof(double));
	for(int i=0;i<nlayers;i++)
		PAtmModel.layer_heights[i]=layer_heights[i];

	return PAtmModel;
}

AtmosphericModel construct_SLCSAT_night_model
		( ThreeFrame TF, int nlayers, double *layer_heights )
{
	AtmosphericModel PAtmModel;
	PAtmModel.type=6;

	PAtmModel.TF=TF;
	PAtmModel.nlayers=nlayers;

	PAtmModel.Gemini_extended_profile=0;
	PAtmModel.Gemini_focal_anisoplanatism_quality=0;
	PAtmModel.Gemini_ground_layer_quality=0;
	PAtmModel.Gemini_outer_scale_meters=0;
	PAtmModel.Hufnagel_Valley_pseudowind=0;

	PAtmModel.r_0_meters=0;
	PAtmModel.r_0_ref_wavelength_meters=0;

	PAtmModel.layer_heights=(double*)malloc(nlayers*sizeof(double));
	for(int i=0;i<nlayers;i++)
		PAtmModel.layer_heights[i]=layer_heights[i];

	return PAtmModel;
}

AtmosphericModel construct_TMT_SRD_v13_Cn2_model( ThreeFrame TF )
{
	AtmosphericModel PAtmModel;

	PAtmModel.type=7;
	PAtmModel.nlayers=0;

	PAtmModel.Gemini_extended_profile=0;
	PAtmModel.Gemini_focal_anisoplanatism_quality=0;
	PAtmModel.Gemini_ground_layer_quality=0;
	PAtmModel.Gemini_outer_scale_meters=0;
	PAtmModel.Hufnagel_Valley_pseudowind=0;

	PAtmModel.r_0_meters=0;
	PAtmModel.r_0_ref_wavelength_meters=0;

	PAtmModel.TF=TF;

	return PAtmModel;
}

AtmosphericModel construct_Gemini_GLAO_study_model
		(ThreeFrame TF, int Gemini_ground_layer_quality, 
		int Gemini_focal_anisoplanatism_quality, 
		int Gemini_extended_profile,
		double Gemini_outer_scale_meters )
{
	AtmosphericModel PAtmModel;
	PAtmModel.type=8;
	PAtmModel.nlayers=0;

	PAtmModel.r_0_meters=0;
	PAtmModel.r_0_ref_wavelength_meters=0;
	PAtmModel.Hufnagel_Valley_pseudowind=0;

	PAtmModel.TF=TF;
	PAtmModel.Gemini_extended_profile
		=Gemini_extended_profile;
	PAtmModel.Gemini_ground_layer_quality
		=Gemini_ground_layer_quality;
	PAtmModel.Gemini_focal_anisoplanatism_quality
		=Gemini_focal_anisoplanatism_quality;
	PAtmModel.Gemini_outer_scale_meters
		=Gemini_outer_scale_meters;

	return PAtmModel;
}

refractive_atmospheric_model
*c2cpp_AtmosphericModel( AtmosphericModel AtmModel )
{

	three_frame tf=c2cpp_ThreeFrame(AtmModel.TF);

	double fried_parameter_meters=AtmModel.r_0_meters;
	double fried_parameter_reference_wavelength_meters=AtmModel.r_0_ref_wavelength_meters;

	double pseudowind=AtmModel.Hufnagel_Valley_pseudowind;
	double Gemini_outer_scale_meters=AtmModel.Gemini_outer_scale_meters;

	bool Gemini_extended_profile = false;
	if ( AtmModel.Gemini_extended_profile)	Gemini_extended_profile = true;

	string Gemini_ground_layer_quality;
	switch(AtmModel.Gemini_ground_layer_quality) 
	{
		case 0:
			Gemini_ground_layer_quality="good";
			break;
		case 1:
			Gemini_ground_layer_quality="typical";
			break;
		case 2:
			Gemini_ground_layer_quality="bad";
			break;
	}

	string Gemini_focal_anisoplanatism_quality;
	switch(AtmModel.Gemini_focal_anisoplanatism_quality) 
	{
		case 0:
			Gemini_focal_anisoplanatism_quality="good";
			break;
		case 1:
			Gemini_focal_anisoplanatism_quality="typical";
			break;
		case 2:
			Gemini_focal_anisoplanatism_quality="bad";
			break;
	}

	refractive_atmospheric_model *ratmdl=NULL; 

	switch(AtmModel.type) 
	{
		case 0:
			{
				vector<double> heights(AtmModel.nlayers);
				vector<power_spectrum *> pspectra( AtmModel.nlayers);
				for(int i=0; i<AtmModel.nlayers; i++)
				{
					heights[i] = AtmModel.layer_heights[i];
					pspectra[i]= c2cpp_PowerSpectrum( AtmModel.power_spectra[i] );
				}

				ratmdl=new refractive_atmospheric_model( pspectra,heights,tf );
				break;
			}
		case 1:
			ratmdl=new Ellerbroek_Cerro_Pachon_model(tf, 
				fried_parameter_meters,
				fried_parameter_reference_wavelength_meters);
			break;
		case 2:
			ratmdl=new Ellerbroek_Mauna_Kea_model(tf, 
				fried_parameter_meters,
				fried_parameter_reference_wavelength_meters);
			break;
		case 3:
			ratmdl=new Palomar_DIMM_MASS_model(tf, 
				fried_parameter_meters,
				fried_parameter_reference_wavelength_meters);
			break;
		case 4:
			{
				vector<double> heights(AtmModel.nlayers);
				for(int i=0; i<AtmModel.nlayers; i++)
					heights[i] = AtmModel.layer_heights[i];
				ratmdl=new Hufnagel_Valley_model(tf, heights, pseudowind);
				break;
			}
		case 5:
			{
				vector<double> heights(AtmModel.nlayers);
				for(int i=0; i<AtmModel.nlayers; i++)
					heights[i] = AtmModel.layer_heights[i];
				ratmdl=new SLCSAT_day_model(tf,	heights );
				break;
			}
		case 6:
			{
				vector<double> heights(AtmModel.nlayers);
				for(int i=0; i<AtmModel.nlayers; i++)
					heights[i] = AtmModel.layer_heights[i];
				ratmdl=new SLCSAT_night_model(tf, heights );
				break;
			}
		case 7:
			ratmdl=new TMT_SRD_v13_Cn2_model(tf);
			break;
		case 8:
			ratmdl=new Gemini_GLAO_study_model(tf,
				Gemini_ground_layer_quality,
				Gemini_focal_anisoplanatism_quality,
				Gemini_extended_profile,
				Gemini_outer_scale_meters);
			break;
	}

	return(ratmdl);
}

int AtmosphericModel_layer_number( AtmosphericModel AtmModel )
{
	refractive_atmospheric_model * ref_atm_model 
		= c2cpp_AtmosphericModel(AtmModel);
	return ref_atm_model->get_number_of_layers();
}

double *AtmosphericModel_layer_heights( AtmosphericModel AtmModel )
{
	int i,nlayers; double *tmp;

	refractive_atmospheric_model * ref_atm_model 
		= c2cpp_AtmosphericModel(AtmModel);

	nlayers=ref_atm_model->get_number_of_layers();
	vector<double> heights( ref_atm_model->get_layer_heights() );

	tmp=(double*)malloc(nlayers*sizeof(double));
	for(i=0;i<nlayers;i++) tmp[i]=heights[i];

	return tmp;
}

WavefrontHeader get_AtmosphericModel_WavefrontHeader( 
	    AtmosphericModel AtmModel, double wavelength, double pixscale, 
		Emitter Emt, APERTURE AP, int layer_foreshortening, int pplan_type )
{
	refractive_atmospheric_model * ref_atm_model = c2cpp_AtmosphericModel(AtmModel);
	emitter *emtr=c2cpp_Emitter(Emt); aperture *ap=c2cpp_APERTURE(AP);
	propagation_plan * pplan=get_propagation_plan(pplan_type);

	diffractive_wavefront_header<double> dwfh = 
		ref_atm_model->get_diffractive_wavefront_header<double>
		( wavelength, pixscale, emtr, ap, layer_foreshortening, pplan );

	WavefrontHeader WfH=cpp2c_WavefrontHeader(dwfh);
	return WfH;
}

RefractiveLayer *get_AtmosphericModel_RefractiveLayer( AtmosphericModel AtmModel, 
		SubharmonicMethod SubM, int num_dwfhs, WavefrontHeader *dwfhdrs,
		double *layer_pixscales, ThreeVector *layer_wind_vectors, double time_interval,
		int layer_axes_wind_vector_aligned, int layer_foreshortening )
{
	refractive_atmospheric_model * ref_atm_model 
		= c2cpp_AtmosphericModel(AtmModel);

	subharmonic_method * subm = c2cpp_SubharmonicMethod(SubM);

	vector<diffractive_wavefront_header<double> > dwfhs( num_dwfhs );
	for(int i=0;i<num_dwfhs;i++) 
		dwfhs[i]=c2cpp_WavefrontHeader(dwfhdrs[i]);

	int nlayers = ref_atm_model->get_number_of_layers();
	vector< double > pixscales(nlayers);
	vector<three_vector> wind_velocities(nlayers);
	for(int i=0;i<nlayers;i++) {
		pixscales[i]=layer_pixscales[i];
		wind_velocities[i]=c2cpp_ThreeVector(layer_wind_vectors[i]);
	}

	vector<refractive_atmospheric_layer<double> > ref_atm_layers;
	ref_atm_model->get_refractive_atmospheric_layers
		( pixscales, *subm,dwfhs, wind_velocities, time_interval,
		layer_axes_wind_vector_aligned, layer_foreshortening, ref_atm_layers);

	RefractiveLayer *Rlayers = (RefractiveLayer *)calloc(nlayers,sizeof(RefractiveLayer));
	for(int i=0;i<ref_atm_layers.size();i++)
		Rlayers[i]=cpp2c_RefractiveLayer(ref_atm_layers[i]);

	return Rlayers;
}

template<class T, class U>
void propagate_wavefront_to_aperture(
	const vector<refractive_atmospheric_layer<T> > & ref_atm_layers, 
	const aperture * ap, int diffractive_near_field_propagator,
	diffractive_wavefront<U> & dwf)
{
	long nlayers = ref_atm_layers.size();

	for(int k=0; k<nlayers; k++){
		try{
			ref_atm_layers[k].transform(dwf);
		} catch(...) {
			cerr << "Error in applying transformation from layer " << k 
				<< " to the wavefront\n";
			exit(-1);
		}
 
		Arroyo::three_point intersection_tp;
		if(k==nlayers-1)
			intersection_tp = ap->get_point_of_intersection(dwf, dwf.z());
		else 
			intersection_tp = ref_atm_layers[k+1].get_point_of_intersection(dwf, dwf.z());

		if(diffractive_near_field_propagator)
			dwf.near_field_angular_propagator((intersection_tp - dwf).length());
		else 
			dwf.geometric_propagator((intersection_tp - dwf).length());
		
		dwf.Arroyo::three_point::operator=(intersection_tp);
	}

}

WaveFront AtmosphericModel_to_APERTURE_transform( 
	int nlayers, RefractiveLayer *Rlayers, APERTURE AP,
	int diffractive_near_field_propagator, WaveFront WF )
{
	int i;

	aperture *ap=c2cpp_APERTURE(AP);
	diffractive_wavefront<double> dwf=c2cpp_WaveFront( WF );
	vector<Arroyo::refractive_atmospheric_layer<double> > ref_atm_layers(nlayers);
	for(i=0;i<nlayers;i++)	ref_atm_layers[i]=c2cpp_RefractiveLayer( Rlayers[i] );

	propagate_wavefront_to_aperture( ref_atm_layers,
		ap, diffractive_near_field_propagator, dwf);

	WaveFront tmp_WF=cpp2c_WaveFront(dwf);
	return tmp_WF;
}

void write_atm_mod_ref_layer_file( RefractiveLayer layer, int n, char *fname )
{
	refractive_atmospheric_layer<double> 
		ref_layer = c2cpp_RefractiveLayer(layer);

	ref_layer.write(number_filename(fname,n,".fits"));
}

double *AtmosphericModel2array ( AtmosphericModel AtmModel )
{
	double *outptr, *tmp;
	int  i, j, type, nlayers;

	type = AtmModel.type; 
	nlayers = AtmModel.nlayers;

	outptr = (double*)calloc( 8*nlayers+21, sizeof(double) );
	
	outptr[0]=nlayers; 
	outptr[1]=type;
	outptr[2]=AtmModel.r_0_meters;
	outptr[3]=AtmModel.r_0_ref_wavelength_meters;
	outptr[4]=AtmModel.Hufnagel_Valley_pseudowind;
	outptr[5]=AtmModel.Gemini_outer_scale_meters;
	outptr[6]=AtmModel.Gemini_ground_layer_quality;
	outptr[7]=AtmModel.Gemini_focal_anisoplanatism_quality;
	outptr[8]=AtmModel.Gemini_extended_profile;

	tmp=ThreeFrame2array( AtmModel.TF );
	for(i=0;i<12;i++) outptr[i+9] = tmp[i];

	if( nlayers != 0 )	{
		for(i=0;i<nlayers;i++)	{
			outptr[i+21] = AtmModel.layer_heights[i];
			if( type==0 ) {
				tmp = PowerSpectrum2array( AtmModel.power_spectra[i] );
				for(j=0;j<7;j++) outptr[7*i+j+nlayers+21] = tmp[j];
			}
		}
	}

	return outptr;
}

AtmosphericModel array2AtmosphericModel (double *inptr)
{
	AtmosphericModel AtmModel;
	int i,type,nlayers;

	nlayers=(int)inptr[0];
	type=(int)inptr[1];

	AtmModel.type=type;
	AtmModel.nlayers=nlayers;
	AtmModel.r_0_meters=inptr[2];
	AtmModel.r_0_ref_wavelength_meters=inptr[3];
	AtmModel.Hufnagel_Valley_pseudowind=inptr[4];
	AtmModel.Gemini_outer_scale_meters=inptr[5];
	AtmModel.Gemini_ground_layer_quality=(int)inptr[6];
	AtmModel.Gemini_focal_anisoplanatism_quality=(int)inptr[7];
	AtmModel.Gemini_extended_profile=(int)inptr[8];

	AtmModel.TF=array2ThreeFrame(inptr+9);
	if( nlayers != 0 )	{
		AtmModel.layer_heights=inptr+21;
		if( type==0 ) {
			AtmModel.power_spectra
				=(PowerSpectrum*) calloc ( nlayers, sizeof( PowerSpectrum ) );
			for(i=0;i<nlayers;i++) 
				AtmModel.power_spectra[i]
				=array2PowerSpectrum(inptr+7*i+nlayers+21);
		}
	}

	return AtmModel;
}
