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

#include <iostream>
#include <fstream>
#include <map>
#include "fits_factory.h"
#include "region_base.h"
#include "diffractive_wavefront.h"
#include "arcadia_mcao_reconstructor.h"
#include "sim_utils.h"
#include "Shack_Hartmann_centroids.h"

using namespace std;

namespace Arroyo {

  int arcadia_mcao_reconstructor::verbose_level = 0;

  namespace factory_register {
    const fits_keyval_set & get_arcadia_mcao_reconstructor_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "arcadia mcao reconstructor"));
      return *fkvs;
    }
    
    AO_sim_base * create_arcadia_mcao_reconstructor(const iofits & iof) {
      return new arcadia_mcao_reconstructor(iof);
    }
  }

  const bool arcadia_mcao_reconstructor::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_arcadia_mcao_reconstructor_keyval_set(), 
				      factory_register::create_arcadia_mcao_reconstructor);


  // Some random functions to help scan the arcadia.in file
  namespace {

    bool parse_boolean(char c){
      if(c=='T' || c=='t') return(true);
      else if(c=='F' || c=='f') return(false);
      else {
	cerr << "parse_boolean error - could not parse character " 
	     << c 
	     << " as true or false\n";
	throw(string("parse_boolean"));
      }
    }

    char output_boolean(bool b){
      if(b) return('t');
      else return('f');
    }
    
    double convert_fortran_exponential_notation(string s){
      double coeff;
      long exp;
      char c;
      stringstream ss(s);
      ss >> coeff;
      ss >> c;
      ss >> exp;
      
      if(c!='d'){
	cerr << "convert_fortran_exponential_notation error - unexpected char " << c << endl;
	throw(string("convert_fortran_exponential_notation"));
      }
      return(coeff*pow(10,(double)exp));
    }

    void parse_arcadia_input_file_error(const char * msg){
      cerr << "arcadia_mcao_reconstructor::arcadia_mcao_reconstructor error - "
	   << msg << endl;
      throw(string("arcadia_mcao_reconstructor::arcadia_mcao_reconstructor"));
    }
    

    class arcadia_wfs_line {

    protected:

      vector<long> subap_indices;
      double subap_area;
      vector<double> snrs;

    public:
     
      arcadia_wfs_line(){
	subap_indices.resize(2);
	snrs.resize(4);
      };

      arcadia_wfs_line(const arcadia_wfs_line & awfsl){
	this->operator=(awfsl);
      }
 
      arcadia_wfs_line(istream & is){
	subap_indices.resize(2);
	is >> subap_indices[0];
	is >> subap_indices[1];
	is >> subap_area;
	snrs.resize(4);
	for(int i=0; i<4; i++)
	  is >> snrs[i];
      }

      arcadia_wfs_line & operator=(const arcadia_wfs_line & awfsl){
	if(this==&awfsl) return(*this);
	subap_indices = awfsl.subap_indices;
	subap_area = awfsl.subap_area;
	snrs = awfsl.snrs;
	return(*this);
      }

      void print(ostream & os) const {
	os << subap_indices[0] << " " 
	   << subap_indices[1] << " " 
	   << subap_area << " ";
	for(int i=0; i<4; i++)
	  os << snrs[i] << " ";
	os << endl;
      }

    };
  }


  arcadia_mcao_reconstructor::
  arcadia_mcao_reconstructor(const arcadia_mcao_reconstructor & arcadia_mcao_recon){
    this->operator=(arcadia_mcao_recon);
  }

  arcadia_mcao_reconstructor::arcadia_mcao_reconstructor(const char * filename){
  }

  arcadia_mcao_reconstructor::arcadia_mcao_reconstructor(const iofits & iof){
    this->read(iof);
  }

  arcadia_mcao_reconstructor::arcadia_mcao_reconstructor(const char * arcadia_input_file,
							 const char * arcadia_mcao_reconstructor_file){
    ////////////////////////////////////
    ///                              ///
    ///  READ THE ARCADIA INPUT FILE ///
    ///                              ///
    ////////////////////////////////////

    /////////////////
    /////////////////
    /////////////////  NOTE - YOU MUST CHECK VARIABLES TO GUARANTEE THE RESTRICTED ASSUMPTIONS OF TAOS
    /////////////////
    /////////////////


    ifstream ifs(arcadia_input_file);

    // Local variables for arcadia.in file
    char tmp_char;
    string tmp_string;
    vector<long> tmp_long_vec(2);

    ////////////////
    // First line //
    ////////////////
    double aperture_diameter_meters;   	      
    bool round_aperture_flag;
    double zenith_angle_degrees;               
    double slew_rate_millirad_per_sec;        
    double reference_wavelength_meters;        
    double nominal_fried_parameter_meters;     
    double nominal_greenwood_frequency_hertz; 
    double outer_scale_meters;                 
    long wind_season;
    bool random_wind_direction_flag; 
    long turbulence_model;  
    vector<double> altitude_limits_meters(2);  
    long number_of_layers;
    bool debug_flag;

    // parse first line
    ifs >> aperture_diameter_meters;
    ifs >> tmp_char;
    round_aperture_flag = parse_boolean(tmp_char);
    ifs >> zenith_angle_degrees;
    ifs >> slew_rate_millirad_per_sec;
    ifs >> tmp_string;
    reference_wavelength_meters = convert_fortran_exponential_notation(tmp_string);
    ifs >> nominal_fried_parameter_meters;
    ifs >> nominal_greenwood_frequency_hertz;
    ifs >> outer_scale_meters;
    ifs >> wind_season;
    ifs >> tmp_char;
    random_wind_direction_flag = parse_boolean(tmp_char);
    ifs >> turbulence_model;
    ifs >> altitude_limits_meters[0];
    ifs >> altitude_limits_meters[1];
    ifs >> number_of_layers;
    ifs >> tmp_char;
    debug_flag = parse_boolean(tmp_char);

    // check that parsing succeeded
    if(!ifs.good()){
      cerr << "arcadia_mcao_reconstructor::arcadia_mcao_reconstructor error - "
	   << "failed to parse first line of file " << arcadia_input_file << endl;
      throw(string("arcadia_mcao_reconstructor::arcadia_mcao_reconstructor"));
    }

    // print out contents if requested
    if(arcadia_mcao_reconstructor::verbose_level)
      cout << aperture_diameter_meters << " "
	   << output_boolean(round_aperture_flag) << " "
	   << zenith_angle_degrees << " "
	   << slew_rate_millirad_per_sec << " "
	   << reference_wavelength_meters << " "
	   << nominal_fried_parameter_meters << " "
	   << nominal_greenwood_frequency_hertz << " "
	   << outer_scale_meters << " "
	   << wind_season << " "
	   << output_boolean(random_wind_direction_flag) << " "
	   << turbulence_model << " "
	   << altitude_limits_meters[0] << " "
	   << altitude_limits_meters[1] << " "
	   << number_of_layers << " "
	   << debug_flag << endl;
      
    // test variables to ensure assumptions are enforced
    if(round_aperture_flag!=true)
      parse_arcadia_input_file_error("round aperture flag set false");
    if(zenith_angle_degrees!=0)
      parse_arcadia_input_file_error("zenith angle differs from zero degrees");
    if(slew_rate_millirad_per_sec!=0)
      parse_arcadia_input_file_error("slew rate differs from zero milliradians per second");
    if(nominal_greenwood_frequency_hertz!=0)
      parse_arcadia_input_file_error("nominal greenwood frequency differs from zero Hz");
    if(wind_season!=7999)
      parse_arcadia_input_file_error("nominal wind season differs from 7999");
    if(random_wind_direction_flag!=true)
      parse_arcadia_input_file_error("random wind direction flag set false");
    if(turbulence_model!=7999)
      parse_arcadia_input_file_error("nominal turbulence model differs from 7999");


    /////////////////
    // Second line //
    /////////////////
    long number_of_evaluation_directions;
    long number_of_evaluation_wavelengths;
    long number_of_noise_cases;
    long nsysjit;
    bool calculate_otf;
    bool calculate_otf_squared;
    bool saevl;
    bool leotf;
    bool seotf;
    bool xyotf;
    bool otfplt;
    bool saotfplt;
    bool eib;
    string otffile;
    string simfile;
    string bwcvrfil;
    double deleib;
    long iotf;
    long myaps;
    long xyaps;

    // parse second line
    ifs >> number_of_evaluation_directions;
    ifs >> number_of_evaluation_wavelengths;
    ifs >> number_of_noise_cases;
    ifs >> nsysjit;
    ifs >> tmp_char;
    calculate_otf = parse_boolean(tmp_char);
    ifs >> tmp_char;
    calculate_otf_squared = parse_boolean(tmp_char);
    ifs >> tmp_char;
    saevl = parse_boolean(tmp_char);
    ifs >> tmp_char;
    leotf = parse_boolean(tmp_char);
    ifs >> tmp_char;
    seotf = parse_boolean(tmp_char);
    ifs >> tmp_char;
    xyotf = parse_boolean(tmp_char);
    ifs >> tmp_char;
    otfplt = parse_boolean(tmp_char);
    ifs >> tmp_char;
    saotfplt = parse_boolean(tmp_char);
    ifs >> tmp_char;
    eib = parse_boolean(tmp_char);
    ifs >> otffile;
    ifs >> simfile;
    ifs >> bwcvrfil;
    ifs >> deleib;
    ifs >> iotf;
    ifs >> myaps;
    ifs >> xyaps;

    // check that parsing succeeded
    if(!ifs.good()){
      cerr << "arcadia_mcao_reconstructor::arcadia_mcao_reconstructor error - "
	   << "failed to parse second line of file " << arcadia_input_file << endl;
      throw(string("arcadia_mcao_reconstructor::arcadia_mcao_reconstructor"));
    }

    // print out contents if requested
    if(arcadia_mcao_reconstructor::verbose_level)
      cout << number_of_evaluation_directions << " "
	   << number_of_evaluation_wavelengths << " "
	   << number_of_noise_cases << " "
	   << nsysjit << " "
	   << output_boolean(calculate_otf) << " "
	   << output_boolean(calculate_otf_squared) << " "
	   << output_boolean(saevl) << " "
	   << output_boolean(leotf) << " "
	   << output_boolean(seotf) << " "
	   << output_boolean(xyotf) << " "
	   << output_boolean(otfplt) << " "
	   << output_boolean(saotfplt) << " "
	   << output_boolean(eib) << " "
	   << otffile << " "
	   << simfile << " "
	   << bwcvrfil << " "
	   << deleib << " "
	   << iotf << " "
	   << myaps << " "
	   << xyaps << endl;

    // test variables to ensure assumptions are enforced
    if(number_of_noise_cases!=1)
      parse_arcadia_input_file_error("number of noise cases does not equal one");
    if(nsysjit!=1)
      parse_arcadia_input_file_error("nsysjit does not equal one");
    if(calculate_otf!=false)
      parse_arcadia_input_file_error("calculate_otf is true");
    if(calculate_otf_squared!=false)
      parse_arcadia_input_file_error("calculate_otf_squared is true");
    if(saevl!=false)
      parse_arcadia_input_file_error("saevl is true");
    if(leotf!=false)
      parse_arcadia_input_file_error("leotf is true");
    if(seotf!=false)
      parse_arcadia_input_file_error("seotf is true");
    if(xyotf!=false)
      parse_arcadia_input_file_error("xyotf is true");
    if(otfplt!=false)
      parse_arcadia_input_file_error("otfplt is true");
    if(saotfplt!=false)
      parse_arcadia_input_file_error("saotfplt is true");
    if(eib!=false)
      parse_arcadia_input_file_error("eib is true");
    if(otffile!="'none'")
      parse_arcadia_input_file_error("otffile differs from 'none'");
    if(simfile!="'none'")
      parse_arcadia_input_file_error("simfile differs from 'none'");
    if(bwcvrfil!="'none'")
      parse_arcadia_input_file_error("bwcvrfil differs from 'none'");
    

      
    // Third line
    vector<double> evaluation_wavelengths_meters(number_of_evaluation_wavelengths);

    for(int i=0; i<number_of_evaluation_wavelengths; i++){
      ifs >> evaluation_wavelengths_meters[i];
    }

    if(arcadia_mcao_reconstructor::verbose_level){
      for(int i=0; i<number_of_evaluation_wavelengths; i++)
	cout << evaluation_wavelengths_meters[i] << " ";
      cout << endl;
    }
      
    // Fourth line
    vector<double> telescope_jitter_levels(3);
    for(int i=0; i<telescope_jitter_levels.size(); i++)
      ifs >> telescope_jitter_levels[i];

    if(arcadia_mcao_reconstructor::verbose_level){
      for(int i=0; i<telescope_jitter_levels.size(); i++)
	cout << telescope_jitter_levels[i] << " ";
      cout << endl;
    }

    // Fifth line
    double aperture_radius_grid_points;       
    double central_obscuration_grid_points;   
    bool finite_altitude_aperture_stop_flag; 
    double aperture_stop_altitude_meters;      
    bool dm_is_system_pupil_flag;
    double minimum_subaperture_area_meters;    
    vector<long> aperture_limits_grid_points(2); 
      
    ifs >> aperture_radius_grid_points;
    ifs >> central_obscuration_grid_points;
    ifs >> tmp_char;
    finite_altitude_aperture_stop_flag = parse_boolean(tmp_char);
    ifs >> aperture_stop_altitude_meters;
    ifs >> tmp_char;
    dm_is_system_pupil_flag = parse_boolean(tmp_char);
    ifs >> minimum_subaperture_area_meters;
    ifs >> aperture_limits_grid_points[0];
    ifs >> aperture_limits_grid_points[1];

    if(arcadia_mcao_reconstructor::verbose_level)  
      cout << aperture_radius_grid_points << " " 
	   << central_obscuration_grid_points << " " 
	   << output_boolean(finite_altitude_aperture_stop_flag) << " " 
	   << aperture_stop_altitude_meters << " " 
	   << output_boolean(dm_is_system_pupil_flag) << " " 
	   << minimum_subaperture_area_meters << " " 
	   << aperture_limits_grid_points[0] << " " 
	   << aperture_limits_grid_points[1] << endl; 

    // Sixth line
    long number_of_dms;

    ifs >> number_of_dms;

    if(arcadia_mcao_reconstructor::verbose_level)  
      cout << number_of_dms << endl;

    // DM variables
    vector<double> dm_conjugate_altitudes_meters(number_of_dms); 
    vector<double> mirror_radii_meters(number_of_dms);           
    vector<bool> zonal_dm_flags(number_of_dms);
    vector<string> yaps_dm_description_files(number_of_dms);
    vector<string> slaving_files(number_of_dms);
    vector<bool> defdm(number_of_dms);
    vector<long> actuator_spacing_grid_points(number_of_dms);
    vector<long> nactuators(number_of_dms);
    vector<vector<pair<long, long> > > actuator_maps(number_of_dms);
    vector<long> nslaved_actuators(number_of_dms);
    vector<vector<pair<long, double> > > actuator_slaving_maps(number_of_dms);

    // LOOP over DM's
    for(int i=0; i<number_of_dms; i++){
      ifs >> dm_conjugate_altitudes_meters[i];
      ifs >> mirror_radii_meters[i];
      ifs >> tmp_char;
      zonal_dm_flags[i] = parse_boolean(tmp_char);
      ifs >> yaps_dm_description_files[i];
      ifs >> slaving_files[i];
      ifs >> tmp_char;
      defdm[i] = parse_boolean(tmp_char);
      ifs >> actuator_spacing_grid_points[i];
      ifs >> nactuators[i];

      if(arcadia_mcao_reconstructor::verbose_level){
	cout << dm_conjugate_altitudes_meters[i] << " "
	     << mirror_radii_meters[i] << " "
	     << output_boolean(zonal_dm_flags[i]) << " "
	     << yaps_dm_description_files[i] << " "
	     << slaving_files[i] << " "
	     << output_boolean(defdm[i]) << endl;
	cout << actuator_spacing_grid_points[i] << " "
	     << nactuators[i] << endl;
      }

      long first, second;
      for(int j=0; j<nactuators[i]; j++){
	ifs >> first;
	ifs >> second;
	actuator_maps[i].push_back(pair<long,long>(first, second));
      }

      ifs >> nslaved_actuators[i];

      double dsecond;
      for(int j=0; j<nslaved_actuators[i]; j++){
	ifs >> first;
	ifs >> second;
	for(int k=0; k<second; k++){
	  ifs >> first;
	  ifs >> dsecond;
	  actuator_slaving_maps[i].push_back(pair<long,double>(first, dsecond));
	}
      }
    }

    string apifile, apofile;
    ifs >> apifile;
    ifs >> apofile;

    if(arcadia_mcao_reconstructor::verbose_level)  
      cout << apifile << " " << apofile << endl;

    long number_of_beacons;
    ifs >> number_of_beacons;
      
    if(arcadia_mcao_reconstructor::verbose_level)  
      cout << number_of_beacons << endl;

    // Loop over evaluation directions
    vector<double> evaluation_altitudes_meters(number_of_evaluation_directions);
    vector<vector<double> > evaluation_directions_radians(number_of_evaluation_directions);
    vector<double> evaluation_outer_obscurations(number_of_evaluation_directions);
    vector<double> evaluation_weights(number_of_evaluation_directions);
    vector<bool> evaluation_flags(number_of_evaluation_directions);
    vector<string> noncommon_path_files(number_of_evaluation_directions);      
    for(int i=0; i<number_of_evaluation_directions; i++){
      ifs >> tmp_string;
      evaluation_altitudes_meters[i] = convert_fortran_exponential_notation(tmp_string);
      evaluation_directions_radians[i].resize(2);
      ifs >> evaluation_directions_radians[i][0];
      ifs >> evaluation_directions_radians[i][1];
      ifs >> evaluation_outer_obscurations[i];
      ifs >> evaluation_weights[i];
      ifs >> tmp_char;
      evaluation_flags[i] = parse_boolean(tmp_char);
      ifs >> noncommon_path_files[i];
	
      if(arcadia_mcao_reconstructor::verbose_level)  
	cout << evaluation_altitudes_meters[i] << " "
	     << evaluation_directions_radians[i][0] << " "
	     << evaluation_directions_radians[i][1] << " "
	     << evaluation_outer_obscurations[i] << " "
	     << evaluation_weights[i] << " "
	     << output_boolean(evaluation_flags[i]) << " "
	     << noncommon_path_files[i] << endl;
    }

    /// Loop over guide stars
    long number_of_guide_stars = number_of_beacons - number_of_evaluation_directions;
    vector<double> guide_star_altitudes_meters(number_of_guide_stars);
    vector<vector<double> > guide_star_directions_radians(number_of_guide_stars);
    vector<double> guide_star_outer_obscurations(number_of_guide_stars);
    vector<double> beacon_subtenses_radians(number_of_guide_stars);
    vector<long> tilt_removal_flags(number_of_guide_stars);
    vector<bool> hartmann_sensor_flags(number_of_guide_stars);
    vector<bool> wfs_stop_at_dm_flags(number_of_guide_stars);
    vector<double> shear_widths(number_of_guide_stars);
    vector<string> yaps_wavefront_sensor_description_files(number_of_guide_stars);
    vector<bool> radiometric_calculation_flags(number_of_guide_stars);
    vector<long> servo_filters(number_of_guide_stars);
    vector<bool> define_wavefront_sensor_flags(number_of_guide_stars);
    vector<long> sampling_rates_hertz(number_of_guide_stars);
    vector<double> pulse_offsets_seconds(number_of_guide_stars);
    vector<double> intentional_time_delays_seconds(number_of_guide_stars);
    vector<vector<double> > beacon_NEA_values(number_of_guide_stars);
    vector<long> isgb(number_of_guide_stars);
    vector<long> nsns(number_of_guide_stars);
    vector<vector<arcadia_wfs_line> > awfs_lines(number_of_guide_stars);
    vector<string> snifiles(number_of_guide_stars);
    vector<string> snofiles(number_of_guide_stars);
      
    for(int i=0; i<number_of_guide_stars; i++){
      ifs >> guide_star_altitudes_meters[i];
      guide_star_directions_radians[i].resize(2);
      ifs >> guide_star_directions_radians[i][0];
      ifs >> guide_star_directions_radians[i][1];
      ifs >> guide_star_outer_obscurations[i];
      ifs >> beacon_subtenses_radians[i];
      ifs >> tilt_removal_flags[i];
      ifs >> tmp_char;
      hartmann_sensor_flags[i] = parse_boolean(tmp_char);
      ifs >> tmp_char;
      wfs_stop_at_dm_flags[i] = parse_boolean(tmp_char);
      ifs >> shear_widths[i];
      ifs >> yaps_wavefront_sensor_description_files[i];
      ifs >> tmp_char;
      radiometric_calculation_flags[i] = parse_boolean(tmp_char);
      ifs >> servo_filters[i];
      ifs >> tmp_char;
      define_wavefront_sensor_flags[i] = parse_boolean(tmp_char);

      if(arcadia_mcao_reconstructor::verbose_level)  
	cout << guide_star_altitudes_meters[i] << " " 
	     << guide_star_directions_radians[i][0] << " "
	     << guide_star_directions_radians[i][1] << " "
	     << guide_star_outer_obscurations[i] << " "
	     << beacon_subtenses_radians[i] << " "
	     << tilt_removal_flags[i] << " "
	     << output_boolean(hartmann_sensor_flags[i]) << " "
	     << output_boolean(wfs_stop_at_dm_flags[i]) << " "
	     << shear_widths[i] << " "
	     << yaps_wavefront_sensor_description_files[i] << " "
	     << output_boolean(radiometric_calculation_flags[i]) << " "
	     << servo_filters[i] << " "
	     << output_boolean(define_wavefront_sensor_flags[i]) << endl;

      ifs >> sampling_rates_hertz[i];
      ifs >> pulse_offsets_seconds[i];
      ifs >> intentional_time_delays_seconds[i];

      if(arcadia_mcao_reconstructor::verbose_level)  
	cout << sampling_rates_hertz[i] << " "
	     << pulse_offsets_seconds[i] << " "
	     << intentional_time_delays_seconds[i] << endl;
	
      beacon_NEA_values[i].resize(2);
      ifs >> beacon_NEA_values[i][0];
      ifs >> beacon_NEA_values[i][1];

      if(arcadia_mcao_reconstructor::verbose_level)  
	cout << beacon_NEA_values[i][0] << " " 
	     << beacon_NEA_values[i][1] << endl;
	
      ifs >> isgb[i];
      ifs >> nsns[i];
	
      if(arcadia_mcao_reconstructor::verbose_level)  
	cout << isgb[i] << " " 
	     << nsns[i] << endl;
	
      awfs_lines[i].resize(nsns[i]);
      for(int j=0; j<nsns[i]; j++){
	awfs_lines[i][j] = arcadia_wfs_line(ifs);
	if(arcadia_mcao_reconstructor::verbose_level)  
	  awfs_lines[i][j].print(cout);
      }

      ifs >> snifiles[i];
      ifs >> snofiles[i];

      if(arcadia_mcao_reconstructor::verbose_level)  
	cout << snifiles[i] << " "
	     << snofiles[i] << endl;
    }


    bool goto_adaptive_optics_flag;
    bool zxform;
    bool input_filter_functions;
    long navg;
    long number_of_servo_filters;
    ifs >> tmp_char;
    goto_adaptive_optics_flag = parse_boolean(tmp_char);
    ifs >> tmp_char;
    zxform = parse_boolean(tmp_char);
    ifs >> tmp_char;
    input_filter_functions = parse_boolean(tmp_char);
    ifs >> navg;
    ifs >> number_of_servo_filters;

    if(arcadia_mcao_reconstructor::verbose_level)  
      cout << output_boolean(goto_adaptive_optics_flag) << " "
	   << output_boolean(zxform) << " "
	   << output_boolean(input_filter_functions) << " " 
	   << navg << " "
	   << number_of_servo_filters << endl;
	
    long nps;
    long nbeacons;
    ifs >> nps;
    ifs >> nbeacons;
      
    if(arcadia_mcao_reconstructor::verbose_level)  
      cout << nps << " "
	   << nbeacons << endl;

    vector<double> beacon_servo_values(nbeacons);
    for(int i=0; i<beacon_servo_values.size(); i++)
      ifs >> beacon_servo_values[i];

    bool least_squares_estimate;
    double least_squares_threshold;
    bool clcon;
    bool supmodes;
    long minslv;
    bool tslv;

    ifs >> tmp_char;
    least_squares_estimate = parse_boolean(tmp_char);
    ifs >> least_squares_threshold;
    ifs >> tmp_char;
    clcon = parse_boolean(tmp_char);
    ifs >> tmp_char;
    supmodes = parse_boolean(tmp_char);
    ifs >> minslv;
    ifs >> tmp_char;
    tslv = parse_boolean(tmp_char);
	
    if(arcadia_mcao_reconstructor::verbose_level)  
      cout << output_boolean(least_squares_estimate) << " "
	   << least_squares_threshold << " "
	   << output_boolean(clcon) << " "
	   << output_boolean(supmodes) << " "
	   << minslv << " "
	   << output_boolean(tslv) << endl;

    string reifile;
    string reofile;
    ifs >> reifile;
    ifs >> reofile;

    if(arcadia_mcao_reconstructor::verbose_level)  
      cout << reifile << " "
	   << reofile << endl;
      
    string gifile;
    string gofile;
    bool snsp;
    bool rcsp;
      
    ifs >> gifile;
    ifs >> gofile;
    ifs >> tmp_char;
    snsp = parse_boolean(tmp_char);
    ifs >> tmp_char;
    rcsp = parse_boolean(tmp_char);
      
    if(arcadia_mcao_reconstructor::verbose_level)  
      cout << gifile << " "
	   << gofile << " "
	   << output_boolean(snsp) << " " 
	   << output_boolean(rcsp) << endl;
      
    string acvrfile;
    string bcvrfile;
    string ccvrfile;

    ifs >> acvrfile;
    ifs >> bcvrfile;
    ifs >> ccvrfile;

    if(arcadia_mcao_reconstructor::verbose_level)  
      cout << acvrfile << " "
	   << bcvrfile << " "
	   << ccvrfile << endl;

    vector<double> layer_heights(number_of_layers);
    vector<double> layer_weights(number_of_layers);
    vector<double> layer_wind_speeds(number_of_layers);
    for(int i=0; i<number_of_layers; i++){
      ifs >> layer_heights[i];
      ifs >> layer_weights[i];
      ifs >> layer_wind_speeds[i];
      if(arcadia_mcao_reconstructor::verbose_level)  
	cout << layer_heights[i] << " "
	     << layer_weights[i] << " "
	     << layer_wind_speeds[i] << endl;
    }
      

    // Other file contents read by subroutines within arcadia (?)

    ifs.close();


    ////////////////////////////////////
    ///                              ///
    ///  READ IN THE RECONSTRUCTOR   ///
    ///                              ///
    ////////////////////////////////////
    
    ifs.open(arcadia_mcao_reconstructor_file);
    ifs.close();


    ////////////////////////////////////
    ///                              ///
    ///  CONSTRUCT THE CLASS MEMBERS ///
    ///                              ///
    ////////////////////////////////////


    // refractive atmospheric model
    // aperture
    // guide stars
    // lenslet arrays
    // tilt removal flags
    // tip tilt mirrors
    // deformable mirrors
    // slaving thresholds
    // grid point spacing
      
  }

  arcadia_mcao_reconstructor & 
  arcadia_mcao_reconstructor::operator=(const arcadia_mcao_reconstructor & arcadia_mcao_recon){
    if(this==&arcadia_mcao_recon)
      return(*this);

    return(*this);
  }

  void arcadia_mcao_reconstructor::read(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "arcadia_mcao_reconstructor::read - "
	   << "error opening file " << filename << endl;
      throw(string("arcadia_mcao_reconstructor::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "arcadia_mcao_reconstructor::read - "
	   << "error reading "
	   << this->unique_name() << " from file " 
	   << filename << endl;
      throw(string("arcadia_mcao_reconstructor::read"));
    }
  }

  void arcadia_mcao_reconstructor::read(const iofits & iof) {
    if(!iof.key_exists("TYPE")){
      cerr << "arcadia_mcao_reconstructor::read error - "
	   << "unrecognized type of file\n";
      throw(string("arcadia_mcao_reconstructor::read"));
    }
    string type, comment;
    iof.read_key("TYPE", type, comment);
    if(type!=this->unique_name()){
      cerr << "arcadia_mcao_reconstructor::read error - file of type " 
	   << type << " rather than of type "
	   << this->unique_name() << endl;
      throw(string("arcadia_mcao_reconstructor::read"));
    }


    // leave iof pointing to the next header, if it exists
    if(iof.get_hdu_num()!=iof.get_num_hdus())
      iof.movrel_hdu(1);
  }

  void arcadia_mcao_reconstructor::write(const char * filename) const {

    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "arcadia_mcao_reconstructor::write - "
	   << "error opening file " << filename << endl;
      throw(string("arcadia_mcao_reconstructor::write"));
    }

    try{this->write(iof);}
    catch(...){
      cerr << "arcadia_mcao_reconstructor::write - "
	   << "error writing "
	   << this->unique_name() << " to file "
	   << filename << endl;
      throw(string("arcadia_mcao_reconstructor::write"));
    }
  }

  void arcadia_mcao_reconstructor::write(iofits & iof) const {

  }

  void arcadia_mcao_reconstructor::print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << this->unique_name()
       << "/" << "object type" << endl;
  }

  void arcadia_mcao_reconstructor::
  reconstruct_zernike_residuals(const vector<Arroyo::Shack_Hartmann_centroids> & shcentroids, 
				vector<zernike> & znke) const {

  }

  void arcadia_mcao_reconstructor::
  reconstruct_residuals(const vector<Arroyo::Shack_Hartmann_centroids> & shcentroids, 
			vector<zernike> & znke, 
			vector<pixel_array<double> > & pixarrs) const {

  }


  void write_arcadia_input_file(const refractive_atmospheric_model & ref_atm_model,
				const circular_aperture & circ_ap,
				const vector<emitter *> & guide_stars,
				const vector<square_lenslet_array> & lenslet_arrays,
				const vector<bool> & tilt_removal_flags,
				const vector<ideal_tip_tilt_mirror<circular_aperture> > & ttms,
				const vector<ideal_deformable_mirror<circular_aperture> > & dms,
				const vector<double> & slaving_thresholds,
				const char * filename){
    


    ofstream ofs(filename);
    double min_height=0, max_height=15460;  // FIX THIS
    ofs << circ_ap.get_diameter() << " "               	  	// aperture_diameter_meters
	<< "T "                                        	  	// round_aperture_flag
	<< "0 "                                        	  	// zenith_angle_degrees
	<< "0 "                                        	  	// slew_rate_millirad_per_sec
	<< "0.5d-6 "                                   	  	// reference_wavelength_meters
	<< "0.160 "                                    	  	// nominal_fried_parameter_meters
	<< "0 "                                        	  	// nominal_greenwood_frequency_hertz
	<< "1e13 "                                     	  	// outer_scale_meters
	<< "7999 "                                     	  	// wind_season
	<< "t "                                        	  	// random_wind_direction_flag
	<< "7999 "                                     	  	// turbulence_model
	<< min_height << " "                           	  	// altitude_limits_meters[0]
	<< max_height << " "                           	  	// altitude_limits_meters[1]
	<< ref_atm_model.get_number_of_layers() << " " 	  	// number_of_layers
	<< "f "                                        	  	// debug_flag
	<< endl;

    ofs << "0 "                                        	  	// number_of_evaluation_directions
	<< "0 "					       	  	// number_of_evaluation_wavelengths
	<< "1 "					       	  	// number_of_noise_cases
	<< "1 "					       	  	// nsysjit
	<< "f "					       	  	// calculate_otf
	<< "f "					       	  	// calculate_otf_squared
	<< "f "					       	  	// saevl
	<< "f "					       	  	// leotf
	<< "f "					       	  	// seotf
	<< "f "					       	  	// xyotf
	<< "f "					       	  	// otfplt
	<< "f "					       	  	// saotfplt
	<< "f "					       	  	// eib
	<< "'none' "				       	  	// otffile
	<< "'none' "				       	  	// simfile
	<< "'none' "                                   	  	// bwcvrfil
 	<< "0.5 "				       	  	// deleib
	<< "2 "					       	  	// iotf
	<< "0 "					       	  	// myaps
	<< "0 "					       	  	// xyaps
	<< endl;


    ofs << "0 0 0\n";                                  	  	// telescope_jitter_levels
      
    int arcadia_grid_size = 101;
    double arcadia_grid_point_spacing = .125;
    double aperture_radius_grid_points = 
      .5*circ_ap.get_diameter()/arcadia_grid_point_spacing;
    int mingrid = (int)(arcadia_grid_size/2+1 - aperture_radius_grid_points);
    int maxgrid = (int)(arcadia_grid_size/2+1 + aperture_radius_grid_points);
    ofs << aperture_radius_grid_points << " "                   // aperture_radius_grid_points
	<< "0 "                                                 // central_obscuration_grid_points
	<< "f "                                                 // finite_altitude_aperture_stop_flag
	<< "0 "                                                 // aperture_stop_altitude_meters
	<< "f "                                                 // dm_is_system_pupil_flag
	<< mingrid << " "                                       // aperture_limits_grid_points[0]
	<< maxgrid                                              // aperture_limits_grid_points[1]
	<< endl;

    ofs << dms.size()                                           // number_of_dms
	<< endl;

    for(int i=0; i<dms.size(); i++){
      ofs << dms[i].three_point::z(circ_ap) << " "              // dm_conjugate_altitudes_meters
	  << dms[i].get_diameter() << " "                       // mirror_radii_meters
	  << "t "                                               // zonal_dm_flags
	  << "'none' "                                          // yaps_dm_description_files
	  << "'inline' "                                        // slaving_files
	  << "t "                                               // defdm
	  << endl;
      
    }




    ofs.close();
  }
}
