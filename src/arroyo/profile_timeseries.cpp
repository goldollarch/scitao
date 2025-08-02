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

#include "profile_timeseries.h"

using namespace std;

namespace Arroyo {

  namespace {
    
    time_t parse_time(stringstream & ss){
 
      time_t unix_time;
      try{
	struct tm utc_tm;
	int tmp;
      
	ss >> utc_tm.tm_year;
	ss >> utc_tm.tm_mon;
	ss >> utc_tm.tm_mday;
	ss >> utc_tm.tm_hour;
	ss >> utc_tm.tm_min;
	ss >> utc_tm.tm_sec;
	if(ss.bad()){
	  cerr << "parse_time error - stringstream in a bad state\n";
	  throw(string("parse_time"));
	}
      
	utc_tm.tm_year -= 1900;
	utc_tm.tm_mon -=1;

	if((unix_time = my_timegm(&utc_tm))==-1){
	  cerr << "parse_time - error converting struct tm to time_t\n";
	  cerr << ss.str() << endl;
	  cerr << "year " << utc_tm.tm_year << endl;
	  cerr << "month " << utc_tm.tm_mon << endl;
	  cerr << "mday " << utc_tm.tm_mday << endl;
	  cerr << "wday " << utc_tm.tm_wday << endl;
	  cerr << "yday " << utc_tm.tm_yday << endl;
	  cerr << "hour " << utc_tm.tm_hour << endl;
	  cerr << "min " << utc_tm.tm_min << endl;
	  cerr << "sec " << utc_tm.tm_sec << endl;
	  cerr << "unix time " << unix_time << endl;
	
	  throw(string("parse_time"));
	}
      } catch(...){
	cerr << "parse_time - error parsing time\n";
	throw(string("parse_time"));
      }
      return(unix_time);
    }
  }


  ///////////////////////////////////////////
  ///  Construct from the bits
  profile_timeseries::profile_timeseries(const vector<time_t> & profile_timestamps,
					 const vector<vector<double> > & layer_heights_meters,
					 const vector<vector<double> > & layer_Cn2_coeffs){
    
    if(profile_timestamps.size()!=layer_heights_meters.size()){
      cerr << "profile_timeseries::profile_timeseries error - "
		<< "mismatch between the number of timestamps " 
		<< profile_timestamps.size()
		<< " and the number of layer height entries "
		<< layer_heights_meters.size()
		<< endl;
      throw(string("profile_timeseries::profile_timeseries"));
    }

    if(profile_timestamps.size()!=layer_Cn2_coeffs.size()){
      cerr << "profile_timeseries::profile_timeseries error - "
		<< "mismatch between the number of timestamps " 
		<< profile_timestamps.size()
		<< " and the number of layer Cn2 coefficient entries "
		<< layer_Cn2_coeffs.size()
		<< endl;
      throw(string("profile_timeseries::profile_timeseries"));
    }

    for(int i=0; i<layer_heights_meters.size(); i++){
      if(layer_heights_meters[i].size()!=layer_Cn2_coeffs[i].size()){
	cerr << "profile_timeseries::profile_timeseries error - "
		  << "mismatch between the number of layers " 
		  << layer_heights_meters[i].size()
		  << " and the number of layer Cn2 coefficients "
		  << layer_Cn2_coeffs[i].size()
		  << " at timestep " 
		  << i
		  << endl;
	throw(string("profile_timeseries::profile_timeseries"));
      }	
    }

    this->profile_timestamps = profile_timestamps;
    this->layer_heights_meters = layer_heights_meters;
    this->layer_Cn2_coeffs = layer_Cn2_coeffs;

  }

  ///////////////////////////////////////////
  ///  Construct from dimm and mass files
  profile_timeseries::profile_timeseries(const char * dimm_filename,
					 const char * mass_filename){


    // The strategy here is to load the mass and dimm data, interpolate the
    // latter onto the timestamps of the former, and save this information
    // in the class data members

    /////////////////////////
    // Parse the DIMM file //
    /////////////////////////

    ifstream dimm_ifs(dimm_filename, ios::in); 
    if(!dimm_ifs.is_open()){
      cerr << "profile_timeseries::profile_timeseries - error opening file " 
		<< dimm_filename
		<< endl;
      throw(string("profile_timeseries::profile_timeseries"));
    }

    // This dodge seems to avoid daylight savings time issues
    // when timestamps are parsed from files and converted
    // with mktime and localtime
    time_t time_tee;
    struct tm tmp_timestamp;    
    time(&time_tee);
    tmp_timestamp = *(gmtime(&time_tee));

    vector<time_t> dimm_timestamps;
    vector<double> dimm_6pt_seeing;
    string dimm_entry, tmp_str;
    double tmp_dbl;
    while(!dimm_ifs.eof()){
      getline(dimm_ifs,dimm_entry);
      if(!dimm_entry.empty() && 
	 dimm_entry[0]!='#'){

	stringstream ss(dimm_entry);
	time_tee = parse_time(ss);
      
	dimm_timestamps.push_back(time_tee);

	ss >> tmp_str >> tmp_str;
	ss >> tmp_dbl;
	ss >> tmp_dbl;

	dimm_6pt_seeing.push_back(tmp_dbl);
      }
    }
    dimm_ifs.close();

    /*
      cerr << "Found "
      << dimm_timestamps.size()
      << " entries in the dimm file "
      << dimm_filename
      << endl;
    */

    /////////////////////////
    // Parse the mass file //
    /////////////////////////

    ifstream mass_ifs(mass_filename, ios::in); 
    if(!mass_ifs.is_open()){
      cerr << "profile_timeseries::profile_timeseries - error opening file " 
		<< mass_filename
		<< endl;
      throw(string("profile_timeseries::profile_timeseries"));
    }

    string mass_entry;
    while(!mass_ifs.eof()){
      getline(mass_ifs,mass_entry);
      if(!mass_entry.empty() && 
	 mass_entry[0]!='#'){

	stringstream ss(mass_entry);
	time_tee = parse_time(ss);
	
	ss >> tmp_str;
	if(tmp_str=="X"){	
	  this->profile_timestamps.push_back(time_tee);
	  this->layer_Cn2_coeffs.push_back(vector<double>(7));
	  ss >> tmp_dbl >> tmp_dbl >> tmp_dbl;
	  for(int i=1; i<7; i++){
	    ss >> tmp_dbl;
	    ss >> this->layer_Cn2_coeffs.back()[i];
	  }
	}
      }
    }
    mass_ifs.close();

    /*
      cerr << "Found "
      << this->profile_timestamps.size()
      << " entries in the mass file "
      << mass_filename
      << endl;
    */

    // Some range checking
    if(dimm_timestamps.size()<1){
      cerr << "profile_timeseries::profile_timeseries - error - no valid dimm lines in file "
		<< dimm_filename
		<< endl;
      throw(string("profile_timeseries::profile_timeseries"));
    }

    if(this->profile_timestamps.size()<1){
      cerr << "profile_timeseries::profile_timeseries - error - no valid mass lines in file "
		<< mass_filename
		<< endl;
      throw(string("profile_timeseries::profile_timeseries"));
    }

    if(dimm_timestamps.back()<this->profile_timestamps[0] ||
       this->profile_timestamps.back()<dimm_timestamps[0]){
      cerr << "profile_timeseries::profile_timeseries - error - no overlapping range of times between files "
		<< dimm_filename
		<< " and "
		<< mass_filename
		<< endl;

      cerr << dimm_timestamps[0] << endl;
      cerr << dimm_timestamps.back() << endl;
    
      cerr << this->profile_timestamps[0] << endl;
      cerr << this->profile_timestamps.back() << endl;
    
      throw(string("profile_timeseries::profile_timeseries"));
    }

  
    int first_mass_index = 0;
    int last_mass_index = profile_timestamps.size();
    while(dimm_timestamps[0]>this->profile_timestamps[first_mass_index])
      first_mass_index++;

    while(dimm_timestamps.back()<this->profile_timestamps[last_mass_index])
      last_mass_index--;


    int post_dimm_index;
    double slope, intercept;
    double fried_parameter;
    double arcsecs_to_radians = M_PI/180./3600.;
    double radians_to_arcsecs = 1/arcsecs_to_radians; 
    double ref_wavelength_meters = .5e-6;
    double wavenumber = 2*M_PI/ref_wavelength_meters;
    double mu_0, ground_layer_coeff;

    double interpolated_dimm_seeing;
    double mass_seeing;
    for(int i=first_mass_index; i<last_mass_index; i++){

      post_dimm_index=0;
      while(dimm_timestamps[post_dimm_index]<this->profile_timestamps[i])
	post_dimm_index++;

      if(dimm_timestamps[post_dimm_index]==this->profile_timestamps[i]){
	interpolated_dimm_seeing = dimm_6pt_seeing[post_dimm_index];
      } else {
	slope = 
	  (dimm_6pt_seeing[post_dimm_index] - dimm_6pt_seeing[post_dimm_index-1])/
	  (dimm_timestamps[post_dimm_index] - dimm_timestamps[post_dimm_index-1]);
	intercept = 
	  dimm_6pt_seeing[post_dimm_index-1] - slope*dimm_timestamps[post_dimm_index-1];

	interpolated_dimm_seeing = slope*profile_timestamps[i] + intercept;
      }

      fried_parameter = ref_wavelength_meters/(interpolated_dimm_seeing*arcsecs_to_radians);

      mu_0 = pow(fried_parameter,-5/3.)/.423 / wavenumber / wavenumber;

      ground_layer_coeff = mu_0;

      double mass_mu_0 = 0;
      for(int j=1; j<7; j++)
	mass_mu_0 += this->layer_Cn2_coeffs[i][j];

      mass_seeing = 
	radians_to_arcsecs * ref_wavelength_meters / 
	pow((.423 * mass_mu_0 * wavenumber * wavenumber),-3/5.);

      ground_layer_coeff = mu_0 - mass_mu_0;

      if(ground_layer_coeff<0) 
	ground_layer_coeff=0;

      this->layer_Cn2_coeffs[i][0] = ground_layer_coeff;
	
    }

    vector<double> tmp_layer_heights_meters(7);
    tmp_layer_heights_meters[0] = 0;
    tmp_layer_heights_meters[1] = 500;
    tmp_layer_heights_meters[2] = 1000;
    tmp_layer_heights_meters[3] = 2000;
    tmp_layer_heights_meters[4] = 4000;
    tmp_layer_heights_meters[5] = 8000;
    tmp_layer_heights_meters[6] = 16000;

    this->layer_heights_meters = vector<vector<double> >(this->profile_timestamps.size(), tmp_layer_heights_meters);
  };

  ///////////////////////////////////////////
  ///  Operator =
  profile_timeseries & profile_timeseries::operator=(const profile_timeseries & pts) {
    this->profile_timestamps = pts.profile_timestamps;
    this->layer_heights_meters = pts.layer_heights_meters;
    this->layer_Cn2_coeffs = pts.layer_Cn2_coeffs;
  };

  ///////////////////////////////////////////
  ///  Print
  void profile_timeseries::print(ostream & os, const char * prefix) const {

    for(int i=0; i<profile_timestamps.size(); i++){
      os << prefix 
	 << setw(20) 
	 << gmtime(&(this->profile_timestamps[i]))
	 << endl;
      for(int j=0; j<this->layer_heights_meters[i].size(); j++)
	os << setw(8) << this->layer_heights_meters[i][j]
	   << setw(12) << this->layer_Cn2_coeffs[i][j]
	   << endl;
    }    
  };

  ///////////////////////////////////////////
  ///  Get timestamp of a particular profile
  time_t profile_timeseries::get_timestamp(int index) const {
    if(index<0 || index>profile_timestamps.size()-1){
      cerr << "profile_timeseries::get_timestamp error - "
		<< "index "
		<< index
		<< " out of range\n";
      throw(string("profile_timeseries::get_timestamp"));
    }
    return(this->profile_timestamps[index]);
  };

  namespace {

    double get_fac(){
      double HJ1 = Arroyo::gamma_function(-5/6.0)/pow(2,8/3.0)/Arroyo::gamma_function(11/6.0);
      double C_A = pow(2,1/3.0)*5/36.0/pow(M_PI,5/3.0)/Arroyo::gamma_function(1/3.0);
      double isoangle_coeff = 2*pow(2*M_PI,8/3.0)*C_A*fabs(HJ1);
      double fac = isoangle_coeff/(2*pow(24/5.*gamma_function(6/5.),5/6.));
      return(fac);
    }
    
    refractive_atmospheric_model make_ref_atm_model(const vector<double> & layer_heights_meters,
						    const vector<double> & layer_Cn2_coeffs){

      double fac = get_fac();
      double ref_wavelength_meters = .5e-6;
      double wavenumber = 2*M_PI/ref_wavelength_meters;
      
      vector<Arroyo::power_spectrum*> pspectra(layer_heights_meters.size());
      for(int j=0; j<pspectra.size(); j++){
	double layer_r_0_meters = pow(fac*wavenumber*wavenumber*layer_Cn2_coeffs[j], -3/5.0);
	pspectra[j] = new isotropic_power_law_spectrum< power_law, null_inner_scale>(power_law(-11/3.,
											       layer_r_0_meters,
											       ref_wavelength_meters),
										     null_inner_scale());
	
      }
      
      refractive_atmospheric_model ref_atm_model(pspectra,
						 layer_heights_meters,
						 three_frame());

      
      for(int i=0; i<pspectra.size(); i++)
	delete pspectra[i];
    
      return(ref_atm_model);
    }
  }


  ///////////////////////////////////////////
  ///  Return a refractive atmospheric model at a particular time
  refractive_atmospheric_model 
  profile_timeseries::get_refractive_atmospheric_model(const time_t & timestamp) const {

    if(timestamp < this->profile_timestamps[0] ||
       timestamp > this->profile_timestamps.back()){
      cerr << "profile_timeseries::get_refractive_atmospheric_model error - timestamp "
		<< timestamp
		<< " out of range "
		<< this->profile_timestamps[0] 
		<< " - " 
		<< this->profile_timestamps.back()
		<< endl;
      throw(string("profile_timeseries::get_refractive_atmospheric_model"));
    }

    vector<double> interpolated_layer_heights_meters;
    vector<double> interpolated_layer_Cn2_coeffs;

    double slope, intercept;
    int post_timestamp_index = 0;
    while(this->profile_timestamps[post_timestamp_index]<timestamp)
      post_timestamp_index++;
    
    if(this->profile_timestamps[post_timestamp_index]==timestamp){
      interpolated_layer_heights_meters = this->layer_heights_meters[post_timestamp_index];
      interpolated_layer_Cn2_coeffs = this->layer_Cn2_coeffs[post_timestamp_index];
    } else {

      if(this->layer_heights_meters[post_timestamp_index-1].size()!=
	 this->layer_heights_meters[post_timestamp_index].size()){
	cerr << "profile_timeseries::get_refractive_atmospheric_model error - "
		  << " mismatched dimensionality of layer heights\n";
	throw(string("profile_timeseries::get_refractive_atmospheric_model"));
      }

      if(this->layer_Cn2_coeffs[post_timestamp_index-1].size()!=
	 this->layer_Cn2_coeffs[post_timestamp_index].size()){
	cerr << "profile_timeseries::get_refractive_atmospheric_model error - "
		  << " mismatched dimensionality of layer Cn2 coefficients\n";
	throw(string("profile_timeseries::get_refractive_atmospheric_model"));
      }

      interpolated_layer_heights_meters.resize(this->layer_Cn2_coeffs[post_timestamp_index].size());
      interpolated_layer_Cn2_coeffs.resize(this->layer_Cn2_coeffs[post_timestamp_index].size());

      for(int i=0; i<this->layer_heights_meters[post_timestamp_index].size(); i++){
	
	// interpolate layer heights
	slope = 
	  (this->layer_heights_meters[post_timestamp_index][i] - this->layer_heights_meters[post_timestamp_index-1][i])/
	  (this->profile_timestamps[post_timestamp_index] - this->profile_timestamps[post_timestamp_index-1]);
	intercept = 
	  this->layer_heights_meters[post_timestamp_index-1][i] - slope*this->profile_timestamps[post_timestamp_index-1];
	
	interpolated_layer_heights_meters[i] = slope*timestamp + intercept;
	
	// interpolate layer Cn2 coeffs
	slope = 
	  (this->layer_Cn2_coeffs[post_timestamp_index][i] - this->layer_Cn2_coeffs[post_timestamp_index-1][i])/
	  (this->profile_timestamps[post_timestamp_index] - this->profile_timestamps[post_timestamp_index-1]);
	intercept = 
	  this->layer_Cn2_coeffs[post_timestamp_index-1][i] - slope*this->profile_timestamps[post_timestamp_index-1];
	
	interpolated_layer_Cn2_coeffs[i] = slope*timestamp + intercept;
      }
    }

    refractive_atmospheric_model ref_atm_model = 
      make_ref_atm_model(interpolated_layer_heights_meters,
			 interpolated_layer_Cn2_coeffs);

    return(ref_atm_model);
  };


  refractive_atmospheric_model 
  profile_timeseries::get_refractive_atmospheric_model(int index) const {
    refractive_atmospheric_model ref_atm_model;
    try{
      int nprofiles = this->get_nprofiles();    
      ref_atm_model = this->get_refractive_atmospheric_model(this->get_timestamp(index));
    } catch(...) {
      cerr << "profile_timeseries::get_refractive_atmospheric_model error\n";
      throw(string("profile_timeseries::get_refractive_atmospheric_model"));
    }
    return(ref_atm_model);
  };
    
  refractive_atmospheric_model 
  profile_timeseries::get_refractive_atmospheric_model(const time_t & start,
						       const time_t & end,
						       int & nprofiles_in_average) const {
    
    try{
      int start_index=0;
      while(this->profile_timestamps[start_index]<start)
	start_index++;

      if(start_index==this->profile_timestamps.size()){
	cerr << "profile_timeseries::get_refractive_atmospheric_model error - "
	     << "start time "
	     << start
	     << "is greater than the timestamp "
	     << end
	     << " on the last profile in the timeseries\n";
	throw(string("profile_timeseries::get_refractive_atmospheric_model"));
      }
      
      int end_index=start_index;
      while(this->profile_timestamps[end_index]<=end)
	end_index++;
      end_index--;

      if(end_index<start_index){
	cerr << "profile_timeseries::get_refractive_atmospheric_model error - "
	     << "no profiles in timeseries lie between times "
	     << start
	     << " and " 
	     << end
	     << endl;
	cerr << "Nearest profile at " << this->profile_timestamps[start_index] << endl;
	throw(string("profile_timeseries::get_refractive_atmospheric_model"));
      }


      int nlayers = layer_Cn2_coeffs[start_index].size();
      vector<double> ave_Cn2_coeffs(nlayers,0);

      for(int i=start_index; i<=end_index; i++){
	if(this->layer_Cn2_coeffs[i].size()!=nlayers){
	  cerr << "profile_timeseries::get_refractive_atmospheric_model error - "
	       << "unequal number of layers between profiles at times "
	       << this->profile_timestamps[start_index] 
	       << " and "
	       << this->profile_timestamps[i] 
	       << endl;
	  throw(string("profile_timeseries::get_refractive_atmospheric_model"));
	}	

	for(int j=0; j<nlayers; j++)
	  ave_Cn2_coeffs[j] += this->layer_Cn2_coeffs[i][j];
      }

      nprofiles_in_average = end_index - start_index + 1;
      for(int j=0; j<nlayers; j++)
	ave_Cn2_coeffs[j] /= (double)(nprofiles_in_average);
      
      refractive_atmospheric_model ref_atm_model = 
	make_ref_atm_model(this->layer_heights_meters[start_index],
			   ave_Cn2_coeffs);

      return(ref_atm_model);

    } catch(...) {
      cerr << "profile_timeseries::get_refractive_atmospheric_model error\n";
      throw(string("profile_timeseries::get_refractive_atmospheric_model"));
    }
  };
}

