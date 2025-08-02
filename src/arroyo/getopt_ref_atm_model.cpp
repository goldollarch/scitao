/*
Arroyo - software for the simulation of electromagnetic wave propagation
through turbulence and optics.

Copyright (c) 2000-2004 California Institute of Technology.  Written by
Dr. Matthew Britton.  For comments or questions about this software,
please contact the author at mbritton@astro.caltech.edu.

This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as  published by the
Free Software Foundation; either version 2 of the License, or (at your
option) any later version.

This program is provided "as is" and distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  In no
event shall California Institute of Technology be liable to any party
for direct, indirect, special, incidental or consequential damages,
including lost profits, arising out of the use of this software and its
documentation, even if the California Institute of Technology has been
advised of the possibility of such damage.   The California Institute of
Technology has no obligation to provide maintenance, support, updates,
enhancements or modifications.  See the GNU General Public License for
more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
*/

#include "getopt_ref_atm_model.h"

namespace Arroyo {

  // A convenience function to encapsulate the help for all refractive atmospheric
  // models supported by Arroyo
  void usage_ref_atm_model(){
    cout << "    -M <code><options>     Select a refractive atmospheric model\n";
    cout << "                           The available models are\n";
    cout << "                             c - Ellerbroek Cerro Pachon\n";
    cout << "                             d - SLC day model\n";
    cout << "                             f - refractive atmospheric model stored in a fits file\n";
    cout << "                             g - Scidar-based Cerro Pachon model used for GLAO modeling\n";
    cout << "                             h - Hufnagel-Valley model\n";
    cout << "                             k - NGAO system design model\n";
    cout << "                             m - Ellerbroek Mauna Kea model\n";
    cout << "                             n - SLC night model\n";
    cout << "                             p - Palomar turbulence model\n";
    cout << "                             s - Single layer model\n";
    cout << "                             t - TMT turbulence model\n";
    cout << "                             u - user defined\n";
    cout << "                           The options for each of the codes are as follows:\n";
    cout << "                             c: <r0> <r0 ref wvlngth>\n";
    cout << "                             d: <nlayers> <height 1> ... <height N>\n";
    cout << "                             f: <filename>\n";
    cout << "                             g: <extended><ground layer><upper layers> O <outer scale>\n";
    cout << "                             h: <pseudowind> <nlayers> <height 1> ... <height N>\n";
    cout << "                             k: <r0> <r0 ref wvlngth>\n";
    cout << "                             m: <r0> <r0 ref wvlngth>\n"; 
    cout << "                             n: <nlayers> <height 1> ... <height N>\n";
    cout << "                             p: <r0> <r0 ref wvlngth>\n";
    cout << "                             s: <height (km)><Cn2 dz>\n";
    cout << "                             t: <none>\n";
    cout << "                             u: <nlayers><height 0 (km)><Cn2 0>...<height n (km)><Cn2 n>\n";
    cout << "                           where <r0> and <r0 ref wvlngth> are in meters and \n";
    cout << "                           the layer heights are in kilometers.  For the Gemini\n";
    cout << "                           GLAO models, the extended keyword is optional, the options\n";
    cout << "                           <ground layer> and <upper layers> take the value good, typical\n";
    cout << "                           or bad, the O is optional and indicates that an outer scale\n"
	 << "                           will follow, and <outer scale> is in meters\n";
    cout << "                           For TMT profiles, the options are\n";
    cout << "                             SRD\n";
    cout << "                             f <filename><index>\n";
    cout << "                           where <filename> is the name of a file containing turbulence\n"
	 << "                           profiles distributed by the TMT site group\n";
  }


  // A convenience function to encapsulate command line parsing for all refractive atmospheric
  // models supported by Arroyo
  refractive_atmospheric_model * parse_ref_atm_model(int argc, 
						     char * const argv[], 
						     three_frame tf){

    char atmospheric_model_code = '\0';
    stringstream ss;

    double fried_parameter_meters;
    double fried_parameter_reference_wavelength_meters;
    long nlayers_in_atmospheric_model;
    double tmp_layer_height_kilometers;
    double pseudowind;
    vector<double> atmospheric_model_layer_heights_meters;
    string Gemini_ground_layer_quality, Gemini_focal_anisoplanatism_quality;
    bool Gemini_extended_profile = false;
    double Gemini_outer_scale_meters = -1;

    string TMT_turbulence_profile_filename;
    string TMT_turbulence_profile_file_entry;
    string junk;
    int TMT_turbulence_profile_index;
    ifstream TMT_ifs;
    int TMT_profile_count;
    double TMT_exponent = -11/3.;
    vector<double> TMT_Cn2_coeffs_m_2_3;
    vector<double> TMT_layer_heights_meters;
    vector<power_spectrum *> TMT_power_spectra;
    // Defined in Sasiela eq 2.18
    double cn2_to_power_law_coefficient_conversion_factor = 
      2*M_PI*5*gamma_function(5/6.)/pow(2,4/3.)/pow(M_PI,3/2.)/9./gamma_function(2/3.);

    double single_layer_height_km;
    double single_layer_cn2_dz_m1_over_3;

    refractive_atmospheric_model * ref_atm_model;

    if(sscanf(optarg, "%c", &atmospheric_model_code)!=1 || 
       (atmospheric_model_code!='c' && 
	atmospheric_model_code!='d' &&
	atmospheric_model_code!='f' &&
	atmospheric_model_code!='g' &&
	atmospheric_model_code!='h' && 
	atmospheric_model_code!='k' && 
	atmospheric_model_code!='m' &&
	atmospheric_model_code!='n' && 
	atmospheric_model_code!='p' && 
	atmospheric_model_code!='s' && 
	atmospheric_model_code!='t' &&
	atmospheric_model_code!='u')){
      cerr << "parse_ref_atm_model - error scanning " 
	   << optarg 
	   << " as an atmospheric model code\n";
      throw(string("parse_ref_atm_model"));
    }
    if(atmospheric_model_code=='c' ||
       atmospheric_model_code=='k' ||
       atmospheric_model_code=='m' ||
       atmospheric_model_code=='p'){
      if(sscanf(argv[optind], "%lf", &fried_parameter_meters)!=1 || 
	 fried_parameter_meters<=0){
	cerr << "parse_ref_atm_model - error scanning " 
	     << argv[optind] 
	     << " as a Fried parameter for the Cerro Pachon atmospheric model\n";
	throw(string("parse_ref_atm_model"));
      }
      optind++;
      if(sscanf(argv[optind], "%lf", &fried_parameter_reference_wavelength_meters)!=1 || 
	 fried_parameter_reference_wavelength_meters<=0){
	cerr << "parse_ref_atm_model - error scanning " 
	     << argv[optind] 
	     << " as a Fried parameter reference wavelength for the Cerro Pachon atmospheric model\n";
	throw(string("parse_ref_atm_model"));
      }
      optind++;

      if(atmospheric_model_code=='c')
	ref_atm_model = new Arroyo::Ellerbroek_Cerro_Pachon_model(tf, 
								  fried_parameter_meters,
								  fried_parameter_reference_wavelength_meters); 
      else if(atmospheric_model_code=='m')
	ref_atm_model = new Arroyo::Ellerbroek_Mauna_Kea_model(tf, 
							       fried_parameter_meters,
							       fried_parameter_reference_wavelength_meters); 
      else if(atmospheric_model_code=='p')
	ref_atm_model = new Arroyo::Palomar_DIMM_MASS_model(tf, 
							    fried_parameter_meters,
							    fried_parameter_reference_wavelength_meters); 
      else if(atmospheric_model_code=='k')
	ref_atm_model = new Arroyo::NGAO_system_design_model(tf,
							     fried_parameter_meters,
							     fried_parameter_reference_wavelength_meters); 
    } else if(atmospheric_model_code=='t'){
      if(string(argv[optind])==string("SRD")){
	ref_atm_model = new Arroyo::TMT_SRD_v13_Cn2_model(tf); 
	optind++;
      } else if(string(argv[optind])==string("f")){
	optind++;
	TMT_turbulence_profile_filename = argv[optind];
	optind++;
	if(sscanf(argv[optind], "%d", &TMT_turbulence_profile_index)!=1 || 
	   TMT_turbulence_profile_index<0){
	  cerr << "parse_ref_atm_model - error scanning " 
	       << argv[optind] 
	       << " as a TMT turbulence profile index\n";
	  throw(string("parse_ref_atm_model"));
	}
	optind++;
	
	TMT_ifs.open(TMT_turbulence_profile_filename.c_str(), ifstream::in);
	if(!TMT_ifs.is_open()){
	  cerr << "parse_ref_atm_model - error opening file "
	       << TMT_turbulence_profile_filename
	       << endl;
	  throw(string("parse_ref_atm_model"));
	}
	
	TMT_profile_count=-1;
	while(!TMT_ifs.eof()){
	  getline(TMT_ifs,TMT_turbulence_profile_file_entry);
	  if(!TMT_turbulence_profile_file_entry.empty() &&
	     TMT_turbulence_profile_file_entry[0]!='#' &&
	     TMT_profile_count<TMT_turbulence_profile_index)
	    TMT_profile_count++;
	  if(TMT_profile_count==TMT_turbulence_profile_index)
	    break;
	}
	
	if(TMT_profile_count!=TMT_turbulence_profile_index){
	  cerr << "parse_ref_atm_model - error finding profile "
	       << TMT_turbulence_profile_index
	       << " in file "
	       << TMT_turbulence_profile_filename
	       << endl;
	  throw(string("parse_ref_atm_model"));
	}
	
	ss.str(TMT_turbulence_profile_file_entry);
	ss >> junk >> junk >> junk >> junk >> junk >> junk;
	if(ss.bad()){
	  cerr << "parse_ref_atm_model - error parsing profile from line\n"
	       << TMT_turbulence_profile_file_entry
	       << endl;
	  throw(string("parse_ref_atm_model"));
	}	  
	    

	TMT_layer_heights_meters.resize(7);
	TMT_layer_heights_meters[0] = 16000;
	TMT_layer_heights_meters[1] = 8000;
	TMT_layer_heights_meters[2] = 4000;
	TMT_layer_heights_meters[3] = 2000;
	TMT_layer_heights_meters[4] = 1000;
	TMT_layer_heights_meters[5] = 500;
	TMT_layer_heights_meters[6] = 0;
	TMT_power_spectra.resize(TMT_layer_heights_meters.size());
	TMT_Cn2_coeffs_m_2_3.resize(7);

	for(int i=TMT_layer_heights_meters.size()-1; i>=0; i--){
	  ss >> TMT_Cn2_coeffs_m_2_3[i];
	  if(ss.bad()){
	    cerr << "parse_ref_atm_model - error parsing Cn2 coefficient " 
		 << TMT_layer_heights_meters.size()-i
		 << " from line\n"
		 << TMT_turbulence_profile_file_entry
		 << endl;
	    throw(string("parse_ref_atm_model"));
	  }	  
	}

	ref_atm_model = new refractive_atmospheric_model(TMT_layer_heights_meters,
							 TMT_Cn2_coeffs_m_2_3,
							 tf);
      }
    } else if(atmospheric_model_code=='f'){
      try{
	ref_atm_model = new refractive_atmospheric_model(argv[optind]);
      } catch(...) {
	cerr << "parse_ref_atm_model - error loading model from file "
	     << argv[optind]
	     << endl;
	throw(string("parse_ref_atm_model"));
      }	
      optind++;
    } else if(atmospheric_model_code=='g'){
      
      if(string(argv[optind])==string("extended")){
	Gemini_extended_profile = true;
	optind++;
      }

      Gemini_ground_layer_quality = argv[optind];
      optind++;
      if(Gemini_ground_layer_quality!="good" &&
	 Gemini_ground_layer_quality!="typical" &&
	 Gemini_ground_layer_quality!="bad"){
	cerr << "parse_ref_atm_model - error scanning Gemini ground layer quality "
	     << Gemini_ground_layer_quality
	     << " as one of good, typical or bad\n";
	exit(-1);
      }

      Gemini_focal_anisoplanatism_quality = argv[optind];
      optind++;
      if(Gemini_focal_anisoplanatism_quality!="good" &&
	 Gemini_focal_anisoplanatism_quality!="typical" &&
	 Gemini_focal_anisoplanatism_quality!="bad"){
	cerr << "parse_ref_atm_model - error scanning Gemini focal anisoplanatism quality "
	     << Gemini_focal_anisoplanatism_quality
	     << " as one of good, typical or bad\n";
	exit(-1);
      }

      if(optind<argc && argv[optind][0]=='O'){
	optind++;
	if(sscanf(argv[optind], "%lf", &Gemini_outer_scale_meters)!=1 ||
	   Gemini_outer_scale_meters<=0){
	  cerr << "parse_ref_atm_model - error scanning " 
	       << argv[optind] 
	       << " as an outer scale for the Gemini profile\n";
	  throw(string("parse_ref_atm_model"));
	}
	optind++;
      }
	  
      ref_atm_model = 
	new Arroyo::Gemini_GLAO_study_model(tf, 
					    Gemini_ground_layer_quality,
					    Gemini_focal_anisoplanatism_quality,
					    Gemini_extended_profile,
					    Gemini_outer_scale_meters);
       
    } else if(atmospheric_model_code=='s'){
      if(sscanf(argv[optind], "%lf", &single_layer_height_km)!=1 || 
	 single_layer_height_km<0){
	cerr << "parse_ref_atm_model - error scanning " 
	     << argv[optind] 
	     << " as a layer height\n";
	throw(string("parse_ref_atm_model"));
      }
      optind++;

      if(sscanf(argv[optind], "%lf", &single_layer_cn2_dz_m1_over_3)!=1 || 
	 single_layer_cn2_dz_m1_over_3<0){
	cerr << "parse_ref_atm_model - error scanning " 
	     << argv[optind] 
	     << " as a layer Cn2 value\n";
	throw(string("parse_ref_atm_model"));
      }
      optind++;

      double exponent = -11/3.;
      power_spectrum * pspec = 
	new isotropic_power_law_spectrum<power_law,null_inner_scale>(power_law(exponent, 
									       single_layer_cn2_dz_m1_over_3*
									       cn2_to_power_law_coefficient_conversion_factor),
								     null_inner_scale());	  	  

      vector<power_spectrum *> pspecs(1,pspec);
      vector<double> layer_heights_meters(1,single_layer_height_km*1000);
      ref_atm_model = new refractive_atmospheric_model(pspecs,
						       layer_heights_meters,
						       tf);
      delete pspec;

    } else if(atmospheric_model_code=='u'){

      if(sscanf(argv[optind], "%ld", &nlayers_in_atmospheric_model)!=1 || 
	 nlayers_in_atmospheric_model<=0){
	cerr << "parse_ref_atm_model - error scanning " 
	     << argv[optind] 
	     << " as a number of layers in the atmospheric model\n";
	throw(string("parse_ref_atm_model"));
      }
      optind++;

      vector<double> layer_heights_meters(nlayers_in_atmospheric_model);
      vector<power_spectrum *> pspecs(nlayers_in_atmospheric_model);

      double exponent = -11/3.;
      double cn2_dz;
      for(int i=0; i<nlayers_in_atmospheric_model; i++){
	
	if(sscanf(argv[optind], "%lf", &layer_heights_meters[i])!=1 || 
	   layer_heights_meters[i]<0){
	  cerr << "parse_ref_atm_model - error scanning " 
	       << argv[optind] 
	       << " as a layer height in meters\n";
	  throw(string("parse_ref_atm_model"));
	}
	layer_heights_meters[i]*=1000.;

	optind++;
	
	
	if(sscanf(argv[optind], "%lf", &cn2_dz)!=1 || 
	   cn2_dz<=0){
	  cerr << "parse_ref_atm_model - error scanning " 
	       << argv[optind] 
	       << " as a Cn2 coefficient\n";
	  throw(string("parse_ref_atm_model"));
	}
	optind++;

	pspecs[i] =
	  new isotropic_power_law_spectrum<power_law,null_inner_scale>(power_law(exponent, 
										 cn2_dz*
										 cn2_to_power_law_coefficient_conversion_factor),
								       null_inner_scale());	  	  
      }


      ref_atm_model = new refractive_atmospheric_model(pspecs,
						       layer_heights_meters,
						       tf);

      for(int i=0; i<pspecs.size(); i++)
	delete pspecs[i];

    } else {
      if(atmospheric_model_code=='h'){
	if(sscanf(argv[optind], "%lf", &pseudowind)!=1 || 
	   pseudowind<=0){
	  cerr << "parse_ref_atm_model - error scanning " 
	       << argv[optind] 
	       << " as a pseudowind for the Hufnagel Valley atmospheric model\n";
	  throw(string("parse_ref_atm_model"));
	}
	optind++;
      }

      if(sscanf(argv[optind], "%ld", &nlayers_in_atmospheric_model)!=1 || 
	 nlayers_in_atmospheric_model<=0){
	cerr << "parse_ref_atm_model - error scanning " 
	     << argv[optind] 
	     << " as a number of atmospheric layers\n";
	throw(string("parse_ref_atm_model"));
      }
      optind++;
      atmospheric_model_layer_heights_meters.resize(nlayers_in_atmospheric_model);
      for(int i=0; i<nlayers_in_atmospheric_model; i++){
	if(sscanf(argv[optind], "%lf", &tmp_layer_height_kilometers)!=1 || 
	   tmp_layer_height_kilometers<0){
	  cerr << "parse_ref_atm_model - error scanning " 
	       << argv[optind] 
	       << " as a height for atmospheric layer " 
	       << i+1 
	       << endl;
	  throw(string("parse_ref_atm_model"));
	}
	optind++;
	atmospheric_model_layer_heights_meters[i] = tmp_layer_height_kilometers*1000;
      }

      if(atmospheric_model_code=='h')
	ref_atm_model = new Arroyo::Hufnagel_Valley_model(tf, 
							  atmospheric_model_layer_heights_meters,
							  pseudowind);
      else if(atmospheric_model_code=='d')
	ref_atm_model = new Arroyo::SLCSAT_day_model(tf,
						     atmospheric_model_layer_heights_meters);
      else if(atmospheric_model_code=='n')
	ref_atm_model = new Arroyo::SLCSAT_night_model(tf,
						       atmospheric_model_layer_heights_meters);
	
    }
    return(ref_atm_model);
  }
} 
