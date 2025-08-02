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
#include "fits_factory.h"
#include "region_base.h"
#include "aplusplus_reconstructor.h"

using namespace std;

namespace Arroyo {

  namespace factory_register {
    const fits_keyval_set & get_a_plusplus_reconstructor_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "A++ reconstructor"));
      return *fkvs;
    }
    
    AO_sim_base * create_a_plusplus_reconstructor(const iofits & iof) {
      return new a_plusplus_reconstructor(iof);
    }
  }

  const bool a_plusplus_reconstructor::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_a_plusplus_reconstructor_keyval_set(), 
				      factory_register::create_a_plusplus_reconstructor);


  a_plusplus_reconstructor::
  a_plusplus_reconstructor(const a_plusplus_reconstructor & appr){
    this->operator=(appr);
  }

  a_plusplus_reconstructor::a_plusplus_reconstructor(const char * filename){
    axes.resize(2);
    this->read(filename);
  }

  a_plusplus_reconstructor::a_plusplus_reconstructor(const iofits & iof){
    axes.resize(2);
    this->read(iof);
  }

  namespace {

    void parse_a_plusplus_ascii_file(const char * filename, int nrows, int ncolumns, double * data){
      fstream fs;

      fs.open(filename);
      if(fs.fail()){
	cerr << "parse_a_plusplus_ascii_file - error opening file " << filename << endl;
	throw(string("parse_a_plusplus_ascii_file"));
      }

      int nelem = nrows*ncolumns;
      for(int i=0; i<nelem; i++){
	fs >> data[i];
      }

      if(fs.fail()){
	cerr << "parse_a_plusplus_ascii_file error - "
	     << "could not parse A++ reconstructor file " << filename << endl;
	throw(string("parse_a_plusplus_ascii_file"));
      } 
      // this must fail - otherwise the file is too large
      double tst;
      fs >> tst;
      if(!fs.fail()){
	cerr << "parse_a_plusplus_ascii_file error - "
	     << "A++ reconstructor file " << filename 
	     << " contains too many entries.  An array of size " 
	     << nrows << "x" << ncolumns << " was expected\n";
	throw(string("parse_a_plusplus_ascii_file")); 
      }
    }

    // when we figure out A++'s binary format, we can write an analogous function here...

  }

  a_plusplus_reconstructor::
  a_plusplus_reconstructor(const char * a_plusplus_reconstructor_file, 
			   int nactuators, int nslopes){

    // ensure that this reconstructor is formed for a square array of
    // actuators and subapertures
    int actuator_grid_dimension = (int)(sqrt((double)nactuators));
    int subaperture_grid_dimension = (int)(sqrt(.5*nslopes));
    
    if(actuator_grid_dimension*actuator_grid_dimension != nactuators){
      cerr << "a_plusplus_reconstructor::a_plusplus_reconstructor error - number of actuators "
	   << nactuators << " passed to this function is not the square of an integer\n";
      throw(string("a_plusplus_reconstructor::a_plusplus_reconstructor"));
    }

    if(2*subaperture_grid_dimension*subaperture_grid_dimension != nslopes){
      cerr << "a_plusplus_reconstructor::a_plusplus_reconstructor error - number of slopes "
	   << nslopes << " passed to this function is not twice the square of an integer\n";
      throw(string("a_plusplus_reconstructor::a_plusplus_reconstructor"));
    }

    vector<long> reconstructor_axes(2,nslopes);
    reconstructor_axes[1] = nactuators;

    this->pixel_array<double>::set_axes(reconstructor_axes);

    try{parse_a_plusplus_ascii_file(a_plusplus_reconstructor_file, 
				    nactuators, nslopes, 
				    this->pixeldata);}
    catch(...){
      cerr << "a_plusplus_reconstructor::a_plusplus_reconstructor error - "
	   << "could not parse A++ reconstructor file\n";
      throw(string("a_plusplus_reconstructor::a_plusplus_reconstructor"));
    } 
  }
  
  a_plusplus_reconstructor::
  a_plusplus_reconstructor(const char * a_plusplus_zonal_reconstructor_file, 
			   int nactuators, int nslopes,
			   const char * a_plusplus_zernike_reconstructor_file, 
			   int nzernike_modes){
    
    // ensure that this reconstructor is formed for a square array of
    // actuators and subapertures
    int actuator_grid_dimension = (int)(sqrt((double)nactuators));
    int subaperture_grid_dimension = (int)(sqrt(.5*nslopes));
    
    if(actuator_grid_dimension*actuator_grid_dimension != nactuators){
      cerr << "a_plusplus_reconstructor::a_plusplus_reconstructor error - number of actuators "
	   << nactuators << " passed to this function is not the square of an integer\n";
      throw(string("a_plusplus_reconstructor::a_plusplus_reconstructor"));
    }

    if(2*subaperture_grid_dimension*subaperture_grid_dimension != nslopes){
      cerr << "a_plusplus_reconstructor::a_plusplus_reconstructor error - number of slopes "
	   << nslopes << " passed to this function is not twice the square of an integer\n";
      throw(string("a_plusplus_reconstructor::a_plusplus_reconstructor"));
    }


    // The mode ordering scheme in A++ appears to be to number
    // by order, and within each order to number from lowest
    // to highest azimuthal mode
    vector<long> zmodes(2,-1);
    for(int i=0; i<nzernike_modes;){
      zmodes[0]++;
      for(zmodes[1]=zmodes[0]%2; zmodes[1]<=zmodes[0] && i<nzernike_modes; zmodes[1]+=2){
	reconstructor_to_zernike_mode_map.
	  insert(pair<int, pair<char, vector<long> > >(i, pair<char, vector<long> >('c', zmodes)));
	i++;
	if(zmodes[1]!=0 && i<nzernike_modes){
	  reconstructor_to_zernike_mode_map.
	    insert(pair<int, pair<char, vector<long> > >(i, pair<char, vector<long> >('s', zmodes)));
	  i++;
	}	  
      }      
    }

    vector<long> reconstructor_axes(2,nslopes);
    reconstructor_axes[1] = nactuators + nzernike_modes;

    this->pixel_array<double>::set_axes(reconstructor_axes);
    
    try{parse_a_plusplus_ascii_file(a_plusplus_zonal_reconstructor_file, 
				    nactuators, nslopes, 
				    this->pixeldata);}
    catch(...){
      cerr << "a_plusplus_reconstructor::a_plusplus_reconstructor error - "
	   << "could not parse A++ reconstructor file\n";
      throw(string("a_plusplus_reconstructor::a_plusplus_reconstructor"));
    } 

    if(nzernike_modes>0){
      try{parse_a_plusplus_ascii_file(a_plusplus_zernike_reconstructor_file, 
				      nzernike_modes, nslopes, 
				      &(this->pixeldata[nactuators*nslopes]));}
      catch(...){
	cerr << "a_plusplus_reconstructor::a_plusplus_reconstructor error - "
	     << "could not parse A++ reconstructor file\n";
	throw(string("a_plusplus_reconstructor::a_plusplus_reconstructor"));
      } 
    }
  }
  
  a_plusplus_reconstructor & 
  a_plusplus_reconstructor::operator=(const a_plusplus_reconstructor & appr){
    if(this==&appr)
      return(*this);
    reconstructor_to_zernike_mode_map = appr.reconstructor_to_zernike_mode_map;
    this->pixel_array<double>::operator=(appr);
    return(*this);
  }

  void a_plusplus_reconstructor::read(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "a_plusplus_reconstructor::read - "
	   << "error opening file " << filename << endl;
      throw(string("a_plusplus_reconstructor::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "a_plusplus_reconstructor::read - "
	   << "error reading "
	   << this->unique_name() << " from file " 
	   << filename << endl;
      throw(string("a_plusplus_reconstructor::read"));
    }
  }

  void a_plusplus_reconstructor::read(const iofits & iof) {
    if(!iof.key_exists("TYPE")){
      cerr << "a_plusplus_reconstructor::read error - "
	   << "unrecognized type of file\n";
      throw(string("a_plusplus_reconstructor::read"));
    }
    string type, comment;
    iof.read_key("TYPE", type, comment);
    if(type!=this->unique_name()){
      cerr << "a_plusplus_reconstructor::read error - file of type " 
	   << type << " rather than of type "
	   << this->unique_name() << endl;
      throw(string("a_plusplus_reconstructor::read"));
    }

    long nprojected_zernikes;
    iof.read_key("NPRJZNKE", nprojected_zernikes, comment);

    int count=0;
    string val, zmode;
    char zchar;
    vector<long> mode_indices(2);
    stringstream ss_key;
    for(int i=0; i<nprojected_zernikes; i++){
      ss_key.str("");
      ss_key << "M" << count;
      iof.read_key(ss_key.str().c_str(), val, comment);

      stringstream ss_val(val);

      zmode.clear();
      ss_val >> zmode >> mode_indices[0] >> mode_indices[1];
      if(zmode.find("c")!=string::npos) zchar = 'c';
      else if(zmode.find("s")!=string::npos) zchar = 's';
      else {
	cerr << "a_plusplus_reconstructor::read error - "
	     << "could not parse reconstructor to zernike mode mapping character code as 'c' or 's'\n";
	throw(string("a_plusplus_reconstructor::read"));
      }

      pair<char,vector<long> > pcv(zchar, mode_indices);
      pair<int,pair<char,vector<long> > > pip(count, pcv);

      reconstructor_to_zernike_mode_map.insert(pip);
      count++;
    }

    this->pixel_array<double>::read(iof);

    // leave iof pointing to the next header, if it exists
    if(iof.get_hdu_num()!=iof.get_num_hdus())
      iof.movrel_hdu(1);
  }

  void a_plusplus_reconstructor::write(const char * filename) const {

    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "a_plusplus_reconstructor::write - "
	   << "error opening file " << filename << endl;
      throw(string("a_plusplus_reconstructor::write"));
    }

    try{this->write(iof);}
    catch(...){
      cerr << "a_plusplus_reconstructor::write - "
	   << "error writing "
	   << this->unique_name() << " to file "
	   << filename << endl;
      throw(string("a_plusplus_reconstructor::write"));
    }
  }

  void a_plusplus_reconstructor::write(iofits & iof) const {
    fits_header_data<double> fhd(this->axes);
    fhd.write(iof);
    string type = this->unique_name();
    string comment = "object type";
    iof.write_key("TYPE", type, comment);

    long nprojected_zernikes = reconstructor_to_zernike_mode_map.size();
    iof.write_key("NPRJZNKE", nprojected_zernikes, "number of projected zernike modes");

    typedef map<int, pair<char, vector<long> > >::const_iterator map_iterator;
    comment = string("reconstructor to zernike mapping");
    stringstream ss_key, ss_val;
    int count=0;
    for(map_iterator p=reconstructor_to_zernike_mode_map.begin(); p!=reconstructor_to_zernike_mode_map.end(); ++p){
      ss_key.str("");
      ss_key << "M" << count;
      ss_val.str("");
      ss_val << p->second.first << " " << p->second.second[0] << " " << p->second.second[1];
      iof.write_key(ss_key.str().c_str(), ss_val.str(), comment);
      count++;
    }

    this->pixel_array<double>::write(iof);
  }

  void a_plusplus_reconstructor::print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << this->unique_name()
       << "/" << "object type" << endl;
    fits_header_data<double> fhd(this->axes);
    fhd.print(os, prefix);

    os << prefix << "NPRJZNKE   = " << setw(vlspc) << reconstructor_to_zernike_mode_map.size()
       << "/" << "number of projected zernikes" << endl;

    typedef map<int, pair<char, vector<long> > >::const_iterator map_iterator;
    stringstream ss_val;
    int count=0;
    for(map_iterator p=reconstructor_to_zernike_mode_map.begin(); p!=reconstructor_to_zernike_mode_map.end(); ++p){
      ss_val.str("");
      ss_val << p->second.first << " " << p->second.second[0] << " " << p->second.second[1];
      os << prefix << "M" << setw(10) << count << "= " << setw(vlspc) << ss_val.str()
	 << "/" << "reconstructor to zernike mapping" << endl;
      count++;
    }
  }

  vector<long> a_plusplus_reconstructor::get_centroid_axes() const {
    vector<long> reconstructor_axes = this->pixel_array<double>::get_axes();
    long subap_dimen = (long)sqrt(.5*reconstructor_axes[0]);
    if(2*subap_dimen*subap_dimen!=reconstructor_axes[0]){
      cerr << "a_plusplus_reconstructor::get_centroid_axes error - "
	   << "inconsistent reconstructor dimensions\n";
      throw(string("a_plusplus_reconstructor::get_centroid_axes"));
    }
    vector<long> centroid_axes(2,subap_dimen);
    centroid_axes[1] *= 2;
    return(centroid_axes);
  }

  vector<long> a_plusplus_reconstructor::get_actuator_axes() const {
    vector<long> reconstructor_axes = this->pixel_array<double>::get_axes();
    long actuator_dimen = (long)sqrt((double)(reconstructor_axes[1]-
				     reconstructor_to_zernike_mode_map.size()));
    if(actuator_dimen*actuator_dimen!=reconstructor_axes[1]-reconstructor_to_zernike_mode_map.size()){
      cerr << "a_plusplus_reconstructor::get_actuator_axes error - "
	   << "inconsistent reconstructor dimensions\n";
      throw(string("a_plusplus_reconstructor::get_actuator_axes"));
    }
    return(vector<long>(2,actuator_dimen));
  }

  zernike a_plusplus_reconstructor::get_zernike_modes() const {

    int max_order = -1;
    typedef map<int, pair<char, vector<long> > >::const_iterator map_iterator;
    for(map_iterator p=reconstructor_to_zernike_mode_map.begin(); p!=reconstructor_to_zernike_mode_map.end(); ++p)
      if(p->second.second[0]>max_order) max_order = p->second.second[0];

    zernike znke(max_order);

    for(map_iterator p=reconstructor_to_zernike_mode_map.begin(); p!=reconstructor_to_zernike_mode_map.end(); ++p){
      if(p->second.first=='c') 
	znke.set_cos_coeff(p->second.second[0],p->second.second[1],1);
      else if(p->second.first=='s'){
	znke.set_sin_coeff(p->second.second[0],p->second.second[1],1);
      } else {
	cerr << "a_plusplus_reconstructor::get_zernike_modes error - unanticipated letter " 
	     << p->second.first << endl;
	throw(string("a_plusplus_reconstructor::get_zernike_modes"));
      }
    }
    return(znke);
  }

  void a_plusplus_reconstructor::
  reconstruct_zernike_residuals(const Arroyo::Shack_Hartmann_centroids & shcentroids, 
				zernike & znke) const {

    if(shcentroids.total_space()!=this->axes[0]){
      cerr << "a_plusplus_reconstructor::reconstruct_zernike_residuals error - "
	   << " Shack Hartmann centroids instance supplied to this function has " 
	   << shcentroids.total_space() << " elements, whereas " 
	   << this->axes[0] << " were expected by the reconstructor\n";
      throw(string("a_plusplus_reconstructor::reconstruct_zernike_residuals"));
    }

    int max_order = -1;
    typedef map<int, pair<char, vector<long> > >::const_iterator map_iterator;
    for(map_iterator p=reconstructor_to_zernike_mode_map.begin(); p!=reconstructor_to_zernike_mode_map.end(); ++p)
      if(p->second.second[0]>max_order) max_order = p->second.second[0];

    if(znke.get_order()!=max_order){
      cerr << "a_plusplus_reconstructor::reconstruct_zernike_residuals error - "
	   << " zernike order supplied to this function " 
	   << znke.get_order() << " does not match that expected by the reconstructor "
	   << max_order << endl;
      throw(string("a_plusplus_reconstructor::reconstruct_zernike_residuals"));
    }

    double resid;
    int level_index = this->axes[0]*(this->axes[1]-reconstructor_to_zernike_mode_map.size());
    long halfaxes = this->axes[0]/2;
    for(map_iterator p=reconstructor_to_zernike_mode_map.begin(); p!=reconstructor_to_zernike_mode_map.end(); ++p){
      resid = 0;

      for(int j=0; j<this->axes[0]; j++)
	resid -= this->pixeldata[level_index+j]*shcentroids.data(j);

      if(p->second.first=='c') 
	znke.set_cos_coeff(p->second.second[0],p->second.second[1],resid);
      else if(p->second.first=='s'){
	znke.set_sin_coeff(p->second.second[0],p->second.second[1],resid);
      } else {
	cerr << "a_plusplus_reconstructor::reconstruct_zernike_residuals error - unanticipated letter " 
	     << p->second.first << endl;
	throw(string("a_plusplus_reconstructor::reconstruct_zernike_residuals"));
      }
      
      level_index += this->axes[0];
    }
  }

  void a_plusplus_reconstructor::
  reconstruct_zonal_residuals(const Arroyo::Shack_Hartmann_centroids & shcentroids, 
			      pixel_array<double> & pixarr) const {

    vector<long> pixarr_axes = pixarr.get_axes();
    if(pixarr_axes[0]*pixarr_axes[1]!=(this->axes[1]-reconstructor_to_zernike_mode_map.size())){
      cerr << "a_plusplus_reconstructor::reconstruct_zonal_residuals error - " << endl
	   << "pixel array instance passed to reconstructor has "
	   << pixarr_axes[0]*pixarr_axes[1] 
	   << " elements, rather than the " 
	   << (this->axes[0]-reconstructor_to_zernike_mode_map.size())
	   << " elements expected by this reconstructor\n";
      throw(string("a_plusplus_reconstructor::reconstruct_zonal_residuals"));
    }

    long index = 0;
    double total=0, tmp;
    long halfaxes = this->axes[0]/2;
    int nzonal_elem = this->axes[1]-reconstructor_to_zernike_mode_map.size();
    int actuator_grid_dimension = (long)(sqrt((double)nzonal_elem));
    
    for(int i=0; i<actuator_grid_dimension; i++){
      for(int j=0; j<actuator_grid_dimension; j++){
	tmp = 0;
	for(int k=0; k<halfaxes; k++){
	  tmp -= this->pixeldata[index+halfaxes+k]*shcentroids.data(k);
	  tmp += this->pixeldata[index+k]*shcentroids.data(k+halfaxes);
	}
	pixarr.set_data(j*actuator_grid_dimension+i, tmp);
	index += this->axes[0];
      }
    }
  }

  void a_plusplus_reconstructor::
  reconstruct_residuals(const Arroyo::Shack_Hartmann_centroids & shcentroids, 
			zernike & znke, 
			pixel_array<double> & pixarr) const {
    this->reconstruct_zernike_residuals(shcentroids, znke);
    this->reconstruct_zonal_residuals(shcentroids, pixarr);
  }

} 
