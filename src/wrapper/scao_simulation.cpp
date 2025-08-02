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

#include <sys/types.h>
#include <sys/stat.h>
#include <ctime>

#include <iostream>
#include <sstream>
#include <vector>

# if (defined( __unix__ ) || (defined(linux)) )
#include <unistd.h>
# else
#include "../arroyo/win-util/some_defines.h"
# endif

#include "simulation.h"
#include "../arroyo/arroyo.h"

using namespace std;

APERTURE AP; aperture *ap;
int reconstructor_nsubaps;
RECONSTRUCTOR arroyo_Reconstructor;
Arroyo::zernike_projected_zonal_reconstructor * zpz_reconstructor;

Emitter sensing_emt,*detected_emt;
emitter *guide_star_emitter;

WavefrontHeader *WfH;
RefractiveLayer *Rlayers;
ZERNIKE TTM_residuals,TTM_commands; 
TipTiltMirror TTM;
TTM_PI_controller ttm_controller;
PixelArray DM_residuals,DM_commands;
DeformableMirror C_DM;
DM_PI_controller dm_controller;
LensletArray lnslt_array;
SHartmannCentroids shcentroid;
Shack_Hartmann_centroids shcentroids;

clock_t ci;
stringstream ss;

char *tmp_char=(char*)malloc(100*sizeof(char));
double arcsecs_to_radians = M_PI/180./3600.;
double rad_to_arcsec = 180*3600/M_PI;

double propagation_distance_meters = 1e10;
double final_pixel_scale_arcsecs;

AtmosphericModel AtmModel;
refractive_atmospheric_model *ref_atm_model;
int num_dwfhs,nlayers;

three_vector mirror_command_vector;
zernike modal_residuals, modal_commands;

double dm_actuator_pitch_meters;
double rescaled_wf_pixel_scale_meters;

diffractive_wavefront<double> dwf, tmp_dwf;
WaveFront WF, tmp_WF;

Arroyo::diffractive_wavefront<double> tt_flat_wf;
WaveFront  tt_flat_WF;

Arroyo::diffractive_wavefront<double> dm_flat_wf;
WaveFront  dm_flat_WF;

Arroyo::pixel_array<double> dm_residuals, dm_commands;

void scao_aperture( SIMULATION_PARAMETER scao_para )
{
	switch(scao_para.aperture_code) 
	{
	case 1:
		AP.ap_data2=scao_para.inner_diameter_meters;
		AP.ap_data1=scao_para.outer_diameter_meters;
		break;
	case 4:
		AP.ap_data2=scao_para.inner_diameter_meters;
		AP.ap_data1=scao_para.outer_diameter_meters;
		AP.ap_data3=scao_para.nspiders;
		AP.ap_data4=scao_para.spider_width;
		break;
	case 5:
		AP.ap_data2=scao_para.inner_diameter_meters;
		AP.ap_data1=scao_para.outer_diameter_meters;
		AP.ap_data3=scao_para.edge_length;
		AP.ap_data4=scao_para.gap_size;
		break;
	default:
		AP.ap_data1=scao_para.outer_diameter_meters;
		break;
	}
    
	ap = c2cpp_APERTURE(AP);
	ap->set_foreshortening(false);

}

void scao_reconstructor( SIMULATION_PARAMETER scao_para, char *rec_file_name )
{

	if(rec_file_name != NULL ){
		zpz_reconstructor = 
			Arroyo::zernike_projected_zonal_reconstructor::
			zernike_projected_zonal_reconstructor_factory(rec_file_name);

		int nmodes = zpz_reconstructor->get_zernike_modes().get_order();

		if(nmodes!=1){
			cerr << "Error = reconstructor " << rec_file_name
				<< " projects out " << nmodes 
				<< " zernike modes, but this program supports projection of tip and tilt only\n";
			exit(-1);
		}
		reconstructor_nsubaps = zpz_reconstructor->get_centroid_axes()[0];
	} else {
		reconstructor_nsubaps = scao_para.reconstructor_nsubaps;
		double pitch = scao_para.outer_diameter_meters/(double)reconstructor_nsubaps;

		Arroyo::ideal_deformable_mirror<Arroyo::circular_aperture>
            temp_DM(Arroyo::circular_aperture(2*scao_para.outer_diameter_meters),
			vector<long>(2,reconstructor_nsubaps+1),  pitch,
			scao_para.dm_actuator_linear_velocity_meters_per_sec);

		Arroyo::square_lenslet_array temp_lnslt_arr(
			vector<long>(2,reconstructor_nsubaps),
			scao_para.lenslet_focal_length_meters,
			pitch,
			scao_para.final_wf_pixels_per_lenslet,
			scao_para.final_wf_pixels_per_xform);

		Arroyo::zernike znke(1);
		znke.set_cos_coeff(1,1,1);
		znke.set_sin_coeff(1,1,1);
		
		zpz_reconstructor = 
			new Arroyo::arroyo_least_squares_reconstructor<float>(*ap,
			temp_DM, temp_lnslt_arr,   znke,  false,
			scao_para.reconstructor_areal_threshold, 
			scao_para.reconstructor_eigenvalue_threshold);

		ss.str("");
		ss << "arroyo_least_squares_reconstructor_" 
			<< reconstructor_nsubaps  << "x" 
			<< reconstructor_nsubaps  << ".fits";
		zpz_reconstructor->write(ss.str().c_str());

  }

	strcpy(tmp_char, ss.str().c_str());
	arroyo_Reconstructor=RECONSTRUCTOR_file_constructor(tmp_char);
  
}

void scao_emitter( SIMULATION_PARAMETER scao_para )
{
	sensing_emt.type=0;
	sensing_emt.plane_wave_emission_direction=cpp2c_ThreeVector(-1*ap->z());

	guide_star_emitter=c2cpp_Emitter(sensing_emt);
	Emitter tmp_emt;
	detected_emt=(Emitter*)malloc(scao_para.n_x_emitters*scao_para.n_y_emitters*sizeof(Emitter));

	int num_science_emitters=scao_para.n_x_emitters*scao_para.n_y_emitters;

	Arroyo::three_rotation trot;
	Arroyo::three_vector emission_direction;

	double x_halfpix = scao_para.n_x_emitters%2==0 ? .5 : 0;
	double x_extrapix = scao_para.n_x_emitters%2==0 ? 0 : 1;
	double y_halfpix = scao_para.n_y_emitters%2==0 ? .5 : 0;
	double y_extrapix = scao_para.n_y_emitters%2==0 ? 0 : 1;
	
	for(int i=-scao_para.n_y_emitters/2; i<scao_para.n_y_emitters/2+y_extrapix; i++){
		for(int j=-scao_para.n_x_emitters/2; j<scao_para.n_x_emitters/2+x_extrapix; j++){ 
            emission_direction = -1*ap->z();
			Arroyo::three_vector rotation_vector = 
				cross_product((j+x_halfpix)*scao_para.emitter_spacing_arcsecs*arcsecs_to_radians*ap->x()-
				(i+y_halfpix)*scao_para.emitter_spacing_arcsecs*arcsecs_to_radians*ap->y(), ap->z());
            if(rotation_vector.length()>Arroyo::three_frame::precision){
				trot = Arroyo::three_rotation(*ap, rotation_vector, -rotation_vector.length());
				trot.transform(emission_direction);
			}

			tmp_emt.type=0;
			tmp_emt.plane_wave_emission_direction= cpp2c_ThreeVector(emission_direction);
			detected_emt[(i+scao_para.n_y_emitters/2)*scao_para.n_x_emitters+j+scao_para.n_x_emitters/2]
				= tmp_emt;
		}
	}
}

void scao_wavefront_header( SIMULATION_PARAMETER scao_para )
{
	int layer_foreshortening=0;
	int pplan_type;
	if(scao_para.diffractive_near_field_propagator)  pplan_type = 1;
	else  pplan_type = 0;

	int num_science_emitters=scao_para.n_x_emitters*scao_para.n_y_emitters;
	num_dwfhs=num_science_emitters*scao_para.detected_wavelengths_number+1;
	WfH=(WavefrontHeader *)malloc(num_dwfhs*sizeof(WavefrontHeader));

	WfH[0] = get_AtmosphericModel_WavefrontHeader( AtmModel,
		scao_para.sensing_wavelength_meters, scao_para.wavefront_pixel_scale_meters,
		sensing_emt,  AP, layer_foreshortening, pplan_type);

	for(int i=0; i<scao_para.detected_wavelengths_number; i++){
		for(int j=0; j<num_science_emitters; j++){
			WfH[i*num_science_emitters+j+1] = 
				get_AtmosphericModel_WavefrontHeader(AtmModel,
				scao_para.detected_wavelengths_meters[i], 
				scao_para.wavefront_pixel_scale_meters,
				detected_emt[j], AP, layer_foreshortening,
				pplan_type);
		} 
	}
}

void scao_ref_atm_layer( SIMULATION_PARAMETER scao_para )
{
	ofstream ofs;
	stringstream ss;

	double azimuth_angle_degrees = 0;
	double zenith_angle_degrees = 0;
	double guide_star_height_meters = 90000;

	double wavelength_meters;

	SubharmonicMethod SubM=create_SubharmonicMethod(
		scao_para.subharmonic_method_code, scao_para.subharmonic_depth, 
		scao_para.generalized_subharmonic_subpixels_per_level,
		scao_para.generalized_subharmonic_subpixels_per_pixel);

	Arroyo::subharmonic_method * subm = c2cpp_SubharmonicMethod(SubM);
	
	Arroyo::wind_model * wm  = 
		new Hardy_wind_model(scao_para.rms_ground_wind_speed_meters_per_sec,
			scao_para.tropospheric_height_meters, 
			scao_para.tropospheric_thickness_meters, 
			scao_para.rms_tropospheric_wind_speed_meters_per_sec);

  //////////////////////////////////////////////////
  // Construct the refractive atmospheric layers  //
  //////////////////////////////////////////////////

# if (defined( __unix__ ) || (defined(linux)) )
	srandom(scao_para.seed);
# else
	srand(scao_para.seed);
# endif

	nlayers = ref_atm_model->get_number_of_layers();
	vector<Arroyo::three_vector> layer_wind_velocities = 
		wm->get_random_wind_vectors(ref_atm_model->get_layer_heights(), *ap);

	double *layer_pixscales;
	ThreeVector *layer_wind_vectors;
	layer_pixscales=(double*)malloc(nlayers*sizeof(double));
	layer_wind_vectors=(ThreeVector *)malloc(nlayers*sizeof(ThreeVector));
	for(int i=0;i<nlayers;i++) {
		layer_pixscales[i]=scao_para.layer_pixel_scale_meters;
		layer_wind_vectors[i]=cpp2c_ThreeVector(layer_wind_velocities[i]);
	}

	double time_interval=scao_para.nsteps_in_simulation*scao_para.timestep_seconds;
	vector<Arroyo::refractive_atmospheric_layer<double> > ref_atm_layers(nlayers);
	
	Rlayers=get_AtmosphericModel_RefractiveLayer
		( AtmModel, SubM,  num_dwfhs, WfH, layer_pixscales,
		layer_wind_vectors,	time_interval,0,0);

	for(int i=0;i<nlayers;i++) 
		ref_atm_layers[i]=c2cpp_RefractiveLayer(Rlayers[i]);

	if(scao_para.vverbose) {
		for(int i=0; i<ref_atm_layers.size(); i++){
			ss.str("");
			ss << "refractive_atmospheric_layer_"  << i << ".fits";			
			ref_atm_layers[i].write(ss.str().c_str());
		}
	}

	for(int i=0; i<=scao_para.detected_wavelengths_number; i++){
		
		wavelength_meters=scao_para.sensing_wavelength_meters;

		if(i<scao_para.detected_wavelengths_number)
			wavelength_meters = scao_para.detected_wavelengths_meters[i];

		ss.str("");
		ss << "atmospheric_params_" << scientific << wavelength_meters  << ".txt";
		ofs.open(ss.str().c_str(), ios_base::out);
		ofs << endl << endl;
		ofs << setw(20)  << "Fried parameter "
			<< ref_atm_model->fried_parameter(wavelength_meters,zenith_angle_degrees)
			<< " meters\n";

		ofs << setw(20) << "Seeing "
			<< ref_atm_model->seeing(wavelength_meters,zenith_angle_degrees)*rad_to_arcsec 
			<< " arcsecs\n";
		ofs << setw(20) << "Isoplanatic angle "
			<< ref_atm_model->isoplanatic_angle(wavelength_meters,zenith_angle_degrees)*rad_to_arcsec 
			<< " arcsecs\n";
		ofs << setw(20)<< "Isokinetic angle " 
			<< ref_atm_model->isokinetic_angle(wavelength_meters,
				scao_para.outer_diameter_meters,
				zenith_angle_degrees)*rad_to_arcsec 
            << " arcsecs\n";
		ofs << setw(20)	<< "Greenwood frequency " 
			<< ref_atm_model->greenwood_frequency(layer_wind_velocities,
					      wavelength_meters, 
					      azimuth_angle_degrees,
					      zenith_angle_degrees)
			<< " Hz\n";
		ofs << setw(20) << "D0 for "<< guide_star_height_meters/1000.<< " km guidestar " 
			<< ref_atm_model->d_0( guide_star_height_meters,
				wavelength_meters, zenith_angle_degrees)
			<< " meters\n";
		ofs.close();
  }

}

void scao_ttm_and_pi( SIMULATION_PARAMETER scao_para )
{
	modal_residuals = Arroyo::zernike(1);
	modal_commands = Arroyo::zernike(1);

	TTM_residuals=create_ZERNIKE(1); TTM_commands=create_ZERNIKE(1);
    APERTURE tmp_AP=create_APERTURE(0,2*scao_para.outer_diameter_meters,0,0,0);
	TTM=create_TipTiltMirror(tmp_AP,scao_para.tip_tilt_angular_velocity_rad_per_sec);

	ttm_controller=create_TTM_PI_controller(
		scao_para.tip_tilt_proportional_gain, 
		scao_para.tip_tilt_integral_gain);
}

void scao_dm_and_pi( SIMULATION_PARAMETER scao_para )
{
	vector<long> dm_actuator_axes(2,reconstructor_nsubaps+1);
	if(dm_actuator_axes[0]!=dm_actuator_axes[1]){
		cout << "Error - reconstructor assumes unexpected deformable mirror actuator axes: "
			<< dm_actuator_axes[0] << "x"  << dm_actuator_axes[1] << endl;
		exit(-1);
	}

	dm_actuator_pitch_meters = scao_para.outer_diameter_meters/(double)(reconstructor_nsubaps);
	
	DM_residuals=create_PixelArray(dm_actuator_axes[0],dm_actuator_axes[1],0);
	DM_commands=create_PixelArray(dm_actuator_axes[0],dm_actuator_axes[1],0);

	APERTURE tmp_AP=create_APERTURE(0,2*scao_para.outer_diameter_meters,0,0,0);
    Arroyo::three_frame downward_oriented_frame(*ap, -1*ap->x(), -1*ap->y(), -1*ap->z());
	ThreeFrame DM_TF=cpp2c_ThreeFrame(downward_oriented_frame);
	tmp_AP.TF=DM_TF;

	C_DM = create_DeformableMirror( tmp_AP, dm_actuator_axes[0],dm_actuator_axes[1], 
		dm_actuator_pitch_meters, scao_para.dm_actuator_linear_velocity_meters_per_sec);

	dm_controller=create_DM_PI_controller(C_DM,scao_para.dm_proportional_gain,scao_para.dm_integral_gain);
}

void scao_lenslet_array( SIMULATION_PARAMETER scao_para )
{
	vector<long> lenslet_axes(2,reconstructor_nsubaps);
	rescaled_wf_pixel_scale_meters = 
		scao_para.wavefront_pixel_scale_meters*scao_para.lenslet_pitch_meters/dm_actuator_pitch_meters;

	lnslt_array=create_LensletArray( lenslet_axes[0],lenslet_axes[1],
		scao_para.lenslet_focal_length_meters,scao_para.lenslet_pitch_meters,
		scao_para.final_wf_pixels_per_lenslet, scao_para.final_wf_pixels_per_xform);
}

void scao_ref_atm_model( SIMULATION_PARAMETER scao_para)
{
	AtmModel=SimuPara_AtmModel(scao_para);
	ref_atm_model=c2cpp_AtmosphericModel(AtmModel);
}

void scao_unaberrated_image( SIMULATION_PARAMETER scao_para)
{
	vector<double> final_pixel_scales_meters(scao_para.detected_wavelengths_number);
	vector<vector<long> > final_array_dimensions(scao_para.detected_wavelengths_number);

	Arroyo::diffractive_wavefront_header<double> tmp_dwfh;
	for(int i=0; i<scao_para.detected_wavelengths_number; i++){
		tmp_dwfh=c2cpp_WavefrontHeader(WfH[1]);
		final_pixel_scales_meters[i] = 
			scao_para.detected_wavelengths_meters[i] * propagation_distance_meters /
			(double) tmp_dwfh.get_axes()[0] / tmp_dwfh.get_pixel_scale() / scao_para.oversampling_factors[i];
        final_pixel_scale_arcsecs = 
			rad_to_arcsec*final_pixel_scales_meters[i]/propagation_distance_meters;
        final_array_dimensions[i] = 
			vector<long>(2,(long)ceil(scao_para.focal_plane_image_sizes_arcsecs[i]/final_pixel_scale_arcsecs));
	}

	int num_science_emitters=scao_para.n_x_emitters*scao_para.n_y_emitters;
	Arroyo::pixel_amp_array<double> palomar_ap_pixamparr;

	for(int i=0; i<scao_para.detected_wavelengths_number; i++){
		try {
			tmp_dwfh=c2cpp_WavefrontHeader(WfH[i*num_science_emitters+1]);
			tmp_dwf = Arroyo::diffractive_wavefront<double>(tmp_dwfh);
			tmp_dwf += complex<double>(1,0);     

			tmp_dwf.three_frame::operator=(*ap);
			ap->transform(tmp_dwf);
            
			tmp_dwf.far_field_fresnel_goertzel_reinsch_propagator(
				propagation_distance_meters, 
				final_pixel_scales_meters[i],
				final_array_dimensions[i]);
			
			ss.str("");
			ss << "unaberrated_focal_wf_" << scientific << setprecision(2) 
				<< scao_para.detected_wavelengths_meters[i] << ".fits";
			tmp_dwf.write(ss.str().c_str());

		} catch(...) {
			cerr << "Error forming unaberrated image " 	<< i   << endl;
			exit(-1);
		} 
	}
}

void scao_detected_aperture
( SIMULATION_PARAMETER scao_para,int i, int j, int k)
{
	int num_science_emitters=scao_para.n_x_emitters*scao_para.n_y_emitters;
	try{
		WF=wave_emit(detected_emt[k],WfH[j*num_science_emitters+k+1]);
		dwf=c2cpp_WaveFront(WF); 
	} catch(...) {
		cerr << "Error emitting wavefront from science emitter " << j << endl;
		exit(-1);
	}

	dwf.set_timestamp(i*scao_para.timestep_seconds);
	WF.WfH.timestamp=i*scao_para.timestep_seconds;
	try{
		WF=AtmosphericModel_to_APERTURE_transform(nlayers,Rlayers,AP,
			scao_para.diffractive_near_field_propagator,WF);
		WF=APERTURE_transform(AP,WF);
		dwf=c2cpp_WaveFront(WF);
	  } catch(...) {
	    cerr << "Error propagating wavefront to aperture\n";
	    exit(-1);
	  }
}

void scao_uncorrected_detected_far_field
( SIMULATION_PARAMETER scao_para,int i, int j, int k)
{

	vector<double> final_pixel_scales_meters(scao_para.detected_wavelengths_number);
	vector<vector<long> > final_array_dimensions(scao_para.detected_wavelengths_number);

	Arroyo::diffractive_wavefront_header<double> tmp_dwfh;
	for(int i=0; i<scao_para.detected_wavelengths_number; i++){
		tmp_dwfh=c2cpp_WavefrontHeader(WfH[1]);
		final_pixel_scales_meters[i] = 
			scao_para.detected_wavelengths_meters[i] * propagation_distance_meters / 
			(double) tmp_dwfh.get_axes()[0] / tmp_dwfh.get_pixel_scale() / scao_para.oversampling_factors[i];
		final_pixel_scale_arcsecs = 
			rad_to_arcsec*final_pixel_scales_meters[i]/propagation_distance_meters;
		final_array_dimensions[i] = 
			vector<long>(2,(long)ceil(scao_para.focal_plane_image_sizes_arcsecs[i]/final_pixel_scale_arcsecs));
	}

	if(scao_para.verbose){
		tmp_dwf = dwf;
	    try{
			tmp_dwf.far_field_fresnel_goertzel_reinsch_propagator(
				propagation_distance_meters, 
				final_pixel_scales_meters[j],
				final_array_dimensions[j]);
		} catch(...) {
			cerr << "Error propagating wavefront to far field\n";
			exit(-1);
		}
	}
}

void scao_correction_detected
( SIMULATION_PARAMETER scao_para,int i, int j, int k)
{
	if(k==0){
		tt_flat_wf = dwf;
		tt_flat_wf.install(Arroyo::pixel_phase_array<double>(tt_flat_wf.get_axes()));
		tt_flat_WF = cpp2c_WaveFront(tt_flat_wf);
	}

	Arroyo::three_vector propagation_vector = dwf.z();
	ThreeVector propagation_TV=cpp2c_ThreeVector(propagation_vector);

	try{
		WF =TipTiltMirror_transform(&TTM, WF);
		dwf=c2cpp_WaveFront(WF);
        if(k==0) {
			tt_flat_WF=TipTiltMirror_transform(&TTM, tt_flat_WF);
			tt_flat_wf=c2cpp_WaveFront(tt_flat_WF);
		}
	} catch(...) {
		cerr << "Error transforming wavefront with tip tilt mirror\n";
	    exit(-1);
	  }

	try{
		dwf.set_propagation_direction(ap->z());
		WF=set_WaveFront_propagation_direction(WF,cpp2c_ThreeVector(ap->z()));
		if(k==0) {
			tt_flat_wf.set_propagation_direction(ap->z());
			tt_flat_WF=set_WaveFront_propagation_direction(tt_flat_WF,cpp2c_ThreeVector(ap->z()));
		}
	} catch(...) {
		cerr << "Error setting propagation direction\n";
	    exit(-1);
	}

	if(k==0){
		dm_flat_wf = tt_flat_wf;
		dm_flat_wf.install(Arroyo::pixel_phase_array<double>(dm_flat_wf.get_axes()));
		dm_flat_WF = cpp2c_WaveFront(dm_flat_wf);
	}

	try{
		WF =DeformableMirror_transform(&C_DM, WF);
		dwf=c2cpp_WaveFront(WF);
		if(k==0){
			dm_flat_WF =DeformableMirror_transform(&C_DM, dm_flat_WF);
			dm_flat_wf=c2cpp_WaveFront(dm_flat_WF);
		}
	} catch(...) {
		cerr << "Error transforming wavefront with deformable mirror\n";
	    exit(-1);
	}

    try{
		dwf.set_propagation_direction(propagation_vector);
		WF=set_WaveFront_propagation_direction(WF,propagation_TV);
	  } catch(...) {
	    cerr << "Error setting propagation direction\n";
	    exit(-1);
	  }
}

void scao_corrected_detected_far_field
( SIMULATION_PARAMETER scao_para,int i, int j, int k)
{
	vector<double> final_pixel_scales_meters(scao_para.detected_wavelengths_number);
	vector<vector<long> > final_array_dimensions(scao_para.detected_wavelengths_number);

	Arroyo::diffractive_wavefront_header<double> tmp_dwfh;
	for(int i=0; i<scao_para.detected_wavelengths_number; i++){
		tmp_dwfh=c2cpp_WavefrontHeader(WfH[1]);
		final_pixel_scales_meters[i] = 
			scao_para.detected_wavelengths_meters[i] * propagation_distance_meters /
			(double) tmp_dwfh.get_axes()[0] / tmp_dwfh.get_pixel_scale() / scao_para.oversampling_factors[i];
 
		final_pixel_scale_arcsecs = 
			rad_to_arcsec*final_pixel_scales_meters[i]/propagation_distance_meters;
        final_array_dimensions[i] = 
			vector<long>(2,(long)ceil(scao_para.focal_plane_image_sizes_arcsecs[i]/final_pixel_scale_arcsecs));
	} 
	
	try{
		dwf.far_field_fresnel_goertzel_reinsch_propagator(propagation_distance_meters,
             final_pixel_scales_meters[j], final_array_dimensions[j]);
		WF=cpp2c_WaveFront(dwf);
	  } catch(...) {
	    cerr << "Error propagating wavefront to far field\n";
	    exit(-1);
	  }

}

void scao_sensing_aperture( SIMULATION_PARAMETER scao_para,int i)
{
	Arroyo::diffractive_wavefront_header<double> tmp_dwfh;
	tmp_dwfh=c2cpp_WavefrontHeader(WfH[0]);

    dwf = guide_star_emitter->emit(tmp_dwfh);
    dwf.set_timestamp(i*scao_para.timestep_seconds);
	WF=cpp2c_WaveFront(dwf);

    try{
		WF=AtmosphericModel_to_APERTURE_transform(nlayers,Rlayers,
			AP,scao_para.diffractive_near_field_propagator,WF);
		WF=APERTURE_transform(AP,WF);

		dwf=c2cpp_WaveFront(WF);

    } catch(...) {
      cerr << "Error propagating wavefront to aperture\n";
      exit(-1);
    }

}

void scao_correction_sensing( SIMULATION_PARAMETER scao_para,int i)
{

	WF=cpp2c_WaveFront(dwf);
    tt_flat_wf = dwf;
    tt_flat_wf.install(Arroyo::pixel_phase_array<double>(tt_flat_wf.get_axes()));
	tt_flat_WF = cpp2c_WaveFront(tt_flat_wf);
	  
    try{
		WF =TipTiltMirror_transform(&TTM, WF);
		dwf=c2cpp_WaveFront(WF);
		tt_flat_WF=TipTiltMirror_transform(&TTM, tt_flat_WF);
		tt_flat_wf=c2cpp_WaveFront(tt_flat_WF);
	} catch(...) {
		cerr << "Error transforming wavefront with tip tilt mirror\n";
		exit(-1);
	}
      
    try{
		dwf.set_propagation_direction(ap->z());
		tt_flat_wf.set_propagation_direction(ap->z());
		ThreeVector tmp_direction=cpp2c_ThreeVector(ap->z());
		WF=set_WaveFront_propagation_direction(WF,tmp_direction);
		tt_flat_WF=set_WaveFront_propagation_direction(tt_flat_WF,tmp_direction);
    } catch(...) {
      cerr << "Error setting propagation direction\n";
      exit(-1);
    }

    dm_flat_wf = tt_flat_wf;
    dm_flat_wf.install(Arroyo::pixel_phase_array<double>(dm_flat_wf.get_axes()));
	dm_flat_WF = cpp2c_WaveFront(dm_flat_wf);
	  
    try{
		WF =DeformableMirror_transform(&C_DM, WF);
		dwf=c2cpp_WaveFront(WF);
		dm_flat_WF =DeformableMirror_transform(&C_DM, dm_flat_WF);
		tt_flat_wf=c2cpp_WaveFront(tt_flat_WF);
	} catch(...) {
		cerr << "Error transforming wavefront with deformable mirror\n";
		exit(-1);
	}

}

void scao_sensing_lenslet_array( SIMULATION_PARAMETER scao_para,int i)
{
    dwf.three_frame::operator=(*ap);
    dwf.set_pixel_scale(rescaled_wf_pixel_scale_meters);
	WF=cpp2c_WaveFront(dwf);

    try{
		WF=LensletArray_transform( lnslt_array, WF );
		dwf=c2cpp_WaveFront(WF);
    } catch(...) {
      cerr << "Error transforming wavefront with lenslet array\n";
      exit(-1);
    }

}

void scao_construct_centroids( SIMULATION_PARAMETER scao_para,int i)
{
	vector<long> centroid_axes = zpz_reconstructor->get_centroid_axes();
	vector<long> lenslet_axes(2,centroid_axes[0]);

    dwf.clip_array((scao_para.final_wf_pixels_per_xform-scao_para.final_wf_pixels_per_lenslet)/2);
	WF=cpp2c_WaveFront(dwf);

    try{
		shcentroid = create_SHartmannCentroids( lenslet_axes[0], lenslet_axes[1], WF );
		shcentroids=c2cpp_SHartmannCentroids(shcentroid);
    } catch(...) {
      cerr << "Error constructing Shack Hartmann residuals\n";
      exit(-1);
    }

}

void scao_reconstruct_residuals( SIMULATION_PARAMETER scao_para,int i)
{
	DM_residuals=arroyo_reconstruct_zonal_residuals
		( arroyo_Reconstructor, shcentroid );
	TTM_residuals = arroyo_reconstruct_zernike_residuals
		( arroyo_Reconstructor, shcentroid );

	dm_residuals=c2cpp_PixelArray(DM_residuals);
	modal_residuals = c2cpp_zernike(TTM_residuals);
}

void scao_update_corrector_mirrors( SIMULATION_PARAMETER scao_para,int i)
{
    try{
	  TTM_commands=update_TTM_PI_controller(&ttm_controller,TTM_residuals,TTM_commands);
	  modal_commands=c2cpp_zernike(TTM_commands);
    } catch(...) {
      cerr << "Error passing residuals to tip tilt controller\n";
      exit(-1);
    }
    
    ThreeVector tmp_command_vector=construct_ThreeVector(-modal_commands.get_cos_coeff(1,1)*arcsecs_to_radians, 
			   -modal_commands.get_sin_coeff(1,1)*arcsecs_to_radians, 1);
    mirror_command_vector =c2cpp_ThreeVector(tmp_command_vector); 
    
    try{
		TipTiltMirror_updata(&TTM,tmp_command_vector,i*scao_para.timestep_seconds);
		mirror_command_vector=c2cpp_ThreeVector(tmp_command_vector);
	} catch(...) {
      cerr << "Error updating tip tilt mirror\n";
      exit(-1);
    }
    
    if(i>=scao_para.nsteps_to_delay_closing_dm_loop){
		try{
			DM_commands=update_DM_PI_controller(&dm_controller,DM_residuals,DM_commands);
			dm_commands=c2cpp_PixelArray(DM_commands);
		} catch(...) {
			cerr << "Error passing residuals to tip tilt controller\n";
			exit(-1);
		}

		dm_commands*=-1;
		DM_commands=cpp2c_PixelArray(dm_commands);

		try{
			DeformableMirror_updata(&C_DM,DM_commands.pixeldata,i*scao_para.timestep_seconds);
			dm_commands=dm_commands=c2cpp_PixelArray(DM_commands);
		} catch(...) {
			cerr << "Error updating deformable mirror\n";
			exit(-1);
		}
      
		dm_commands*=-1;
		DM_commands=cpp2c_PixelArray(dm_commands);
	} 
}

void scao_write_wavefront( int type,char *fname,double time,int k)
{
	ss.str("");

	if(!k) ss << fname;
	else ss << fname<< k << "_" ;

	ss << scientific << setprecision(2) 
		<< dwf.get_wavelength();
	if(time<0) ss << ".fits";
	else ss<< "_" << fixed << setprecision(3) 
		<< setfill('0') << time << ".fits";
	
	if(type==0) 
		dwf.write(ss.str().c_str());
	else if(type==1) 
		tmp_dwf.write(ss.str().c_str());
	else if(type==2) 
		tt_flat_wf.write(ss.str().c_str());
	else if(type==3) 
		dm_flat_wf.write(ss.str().c_str());

}
