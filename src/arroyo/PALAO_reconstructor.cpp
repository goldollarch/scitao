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
#include "PALAO_reconstructor.h"

using namespace std;

namespace Arroyo {
  namespace factory_register {
    const fits_keyval_set & get_PALAO_reconstructor_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "PALAO reconstructor"));
      return *fkvs;
    }
    
    AO_sim_base * create_PALAO_reconstructor(const iofits & iof) {
      return new PALAO_reconstructor(iof);
    }
  }

  const bool PALAO_reconstructor::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_PALAO_reconstructor_keyval_set(), 
				      factory_register::create_PALAO_reconstructor);


  PALAO_reconstructor::PALAO_reconstructor(const PALAO_reconstructor & palao_recon){
    this->operator=(palao_recon);
  }

  PALAO_reconstructor::PALAO_reconstructor(const char * filename){
    axes.resize(2);
    this->read(filename);
  }

  PALAO_reconstructor::PALAO_reconstructor(const iofits & iof){
    axes.resize(2);
    this->read(iof);
  }

  PALAO_reconstructor & PALAO_reconstructor::operator=(const PALAO_reconstructor & palao_recon){
    if(this==&palao_recon)
      return(*this);
    this->pixel_array<double>::operator=(palao_recon);
    return(*this);
  }

  void PALAO_reconstructor::read(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "PALAO_reconstructor::read - "
	   << "error opening file " << filename << endl;
      throw(string("PALAO_reconstructor::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "PALAO_reconstructor::read - "
	   << "error reading "
	   << this->unique_name() << " from file " 
	   << filename << endl;
      throw(string("PALAO_reconstructor::read"));
    }
  }

  void PALAO_reconstructor::read(const iofits & iof) {
    if(!iof.key_exists("TYPE")){
      cerr << "PALAO_reconstructor::read error - "
	   << "unrecognized type of file\n";
      throw(string("PALAO_reconstructor::read"));
    }
    string type, comment;
    iof.read_key("TYPE", type, comment);
    if(type!=this->unique_name()){
      cerr << "PALAO_reconstructor::read error - file of type " 
	   << type << " rather than of type "
	   << this->unique_name() << endl;
      throw(string("PALAO_reconstructor::read"));
    }

    string format;
    iof.read_key("FORMAT", format, comment);
    if(format!=string("May 1999") && format!=string("March 2003")){
      cerr << "PALAO_reconstructor::read error - unrecognized format " 
	   << format << endl;
      throw(string("PALAO_reconstructor"));
    }

    this->pixel_array<double>::read(iof);

    if(format=="May 1999"){
      int nelem = 512*241;
      for(int i=1; i<nelem; i+=2)
	this->pixeldata[i] *= -1;
      
      for(int i=nelem; i<nelem+512; i+=2)
	this->pixeldata[i] *= -1;

    }
   
    // leave iof pointing to the next header, if it exists
    if(iof.get_hdu_num()!=iof.get_num_hdus())
      iof.movrel_hdu(1);
  }

  void PALAO_reconstructor::write(const char * filename) const {

    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "PALAO_reconstructor::write - "
	   << "error opening file " << filename << endl;
      throw(string("PALAO_reconstructor::write"));
    }

    try{this->write(iof);}
    catch(...){
      cerr << "PALAO_reconstructor::write - "
	   << "error writing "
	   << this->unique_name() << " to file "
	   << filename << endl;
      throw(string("PALAO_reconstructor::write"));
    }
  }

  void PALAO_reconstructor::write(iofits & iof) const {
    fits_header_data<double> fhd(this->axes);
    fhd.write(iof);
    string type = this->unique_name();
    string comment = "object type";
    iof.write_key("TYPE", type, comment);
    string format("March 2003");
    comment = "reconstructor_format";
    iof.write_key("FORMAT", format, comment);
    this->pixel_array<double>::write(iof);
  }

  void PALAO_reconstructor::print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << this->unique_name()
       << "/" << "object type" << endl;
    fits_header_data<double> fhd(this->axes);
    fhd.print(os, prefix);
  }

  vector<long> PALAO_reconstructor::get_centroid_axes() const {
    vector<long> v(2,16);
    v[1]=32;
    return(v);
  }

  vector<long> PALAO_reconstructor::get_actuator_axes() const {
    return(vector<long>(2,17));
  }

  zernike PALAO_reconstructor::get_zernike_modes() const {
    zernike znke(1);
    znke.set_cos_coeff(0,0,0);
    znke.set_cos_coeff(1,1,1);
    znke.set_sin_coeff(1,1,1);
    return(znke);
  }

  void PALAO_reconstructor::reconstruct_zernike_residuals(const Arroyo::Shack_Hartmann_centroids & shcentroids, 
							  zernike & znke) const {

    if(znke.get_order()!=1){
      cerr << "PALAO_reconstructor::reconstruct_zernike_residuals error - "
	   << "zernike instance passed to this function has maximum order "
	   << znke.get_order() << " rather than order 1 expected by this function\n";
      throw(string("PALAO_reconstructor::reconstruct_zernike_residuals"));
    }

    if(centroid_weighting_lookup_table.empty())
      this->create_centroid_weighting_lookup_table();

    // reconstruct the tip tilt residuals
    int index;
    double a_tip_tilt_resid = 0;
    double b_tip_tilt_resid = 0;
    int tta_index=241*this->axes[0];
    int ttb_index=242*this->axes[0]+1;
    for(int i=0; i<this->axes[0]/2; i++){
      a_tip_tilt_resid -= this->pixeldata[tta_index+2*i]*shcentroids.data(i)*centroid_weighting_lookup_table[i];
      b_tip_tilt_resid += this->pixeldata[ttb_index+2*i]*shcentroids.data(i+256)*centroid_weighting_lookup_table[i];
    }
    znke.set_cos_coeff(1,1,a_tip_tilt_resid);
    znke.set_sin_coeff(1,1,b_tip_tilt_resid);
  }

  void PALAO_reconstructor::reconstruct_zonal_residuals(const Arroyo::Shack_Hartmann_centroids & shcentroids, 
							pixel_array<double> & pixarr) const {

    if(dm_actuator_lookup_table.empty())
      this->create_dm_actuator_lookup_table();

    int nactuators = 17;
    int ncentroids = 16;
    vector<long> centroid_axes = shcentroids.get_axes();
    if(centroid_axes[0]!=ncentroids || centroid_axes[1]!=2*ncentroids){
      cerr << "PALAO_reconstructor::reconstruct_zonal_residuals error - " << endl
	   << "Shack Hartmann centroids instance passed to reconstructor has dimensions "
	   << centroid_axes[0] << "x" << centroid_axes[1] 
	   << "rather than dimensions " 
	   << ncentroids << "x" << 2*ncentroids
	   << " expected by this function\n";
      throw(string("PALAO_reconstructor::reconstruct_zonal_residuals"));
    }

    vector<long> pixarr_axes = pixarr.get_axes();
    if(pixarr_axes[0]!=nactuators || pixarr_axes[1]!=nactuators){
      cerr << "PALAO_reconstructor::reconstruct_zonal_residuals error - " << endl
	   << "pixel array instance passed to reconstructor has dimensions "
	   << pixarr_axes[0] << "x" << pixarr_axes[1] 
	   << "rather than dimensions " 
	   << nactuators << "x" << nactuators 
	   << " expected by this function\n";
      throw(string("PALAO_reconstructor::reconstruct_zonal_residuals"));
    }

    long index, nelem = ncentroids*ncentroids;
    long nactive_actuators = 241;
    double total=0, tmp;
    for(int i=0; i<nactive_actuators; i++){
      tmp = 0;
      for(int j=0; j<nelem; j++){
	tmp += shcentroids.data(j+nelem)*this->pixeldata[2*(i*nelem+j)];
	tmp -= shcentroids.data(j)*this->pixeldata[2*(i*nelem+j)+1];
      }
      total += tmp;
      pixarr.set_data(nactuators*(dm_actuator_lookup_table[i]%nactuators)+dm_actuator_lookup_table[i]/nactuators,tmp);
    }
    // remove the piston component
    //total /= (double)nactive_actuators;
    //pixarr -= total;
  }

  void PALAO_reconstructor::reconstruct_residuals(const Arroyo::Shack_Hartmann_centroids & shcentroids, 
						  zernike & znke, 
						  pixel_array<double> & pixarr) const {
    this->reconstruct_zernike_residuals(shcentroids, znke);
    this->reconstruct_zonal_residuals(shcentroids, pixarr);
  }

  void PALAO_reconstructor::create_dm_actuator_lookup_table() const {
    dm_actuator_lookup_table.resize(241);

    dm_actuator_lookup_table[0] = 5;
    dm_actuator_lookup_table[1] = 6;
    dm_actuator_lookup_table[2] = 7;
    dm_actuator_lookup_table[3] = 8;
    dm_actuator_lookup_table[4] = 9;
    dm_actuator_lookup_table[5] = 10;
    dm_actuator_lookup_table[6] = 11;
    dm_actuator_lookup_table[7] = 20;
    dm_actuator_lookup_table[8] = 21;
    dm_actuator_lookup_table[9] = 22;
    dm_actuator_lookup_table[10] = 23;
    dm_actuator_lookup_table[11] = 24;
    dm_actuator_lookup_table[12] = 25;
    dm_actuator_lookup_table[13] = 26;
    dm_actuator_lookup_table[14] = 27;
    dm_actuator_lookup_table[15] = 28;
    dm_actuator_lookup_table[16] = 29;
    dm_actuator_lookup_table[17] = 30;
    dm_actuator_lookup_table[18] = 36;
    dm_actuator_lookup_table[19] = 37;
    dm_actuator_lookup_table[20] = 38;
    dm_actuator_lookup_table[21] = 39;
    dm_actuator_lookup_table[22] = 40;
    dm_actuator_lookup_table[23] = 41;
    dm_actuator_lookup_table[24] = 42;
    dm_actuator_lookup_table[25] = 43;
    dm_actuator_lookup_table[26] = 44;
    dm_actuator_lookup_table[27] = 45;
    dm_actuator_lookup_table[28] = 46;
    dm_actuator_lookup_table[29] = 47;
    dm_actuator_lookup_table[30] = 48;
    dm_actuator_lookup_table[31] = 52;
    dm_actuator_lookup_table[32] = 53;
    dm_actuator_lookup_table[33] = 54;
    dm_actuator_lookup_table[34] = 55;
    dm_actuator_lookup_table[35] = 56;
    dm_actuator_lookup_table[36] = 57;
    dm_actuator_lookup_table[37] = 58;
    dm_actuator_lookup_table[38] = 59;
    dm_actuator_lookup_table[39] = 60;
    dm_actuator_lookup_table[40] = 61;
    dm_actuator_lookup_table[41] = 62;
    dm_actuator_lookup_table[42] = 63;
    dm_actuator_lookup_table[43] = 64;
    dm_actuator_lookup_table[44] = 65;
    dm_actuator_lookup_table[45] = 66;
    dm_actuator_lookup_table[46] = 69;
    dm_actuator_lookup_table[47] = 70;
    dm_actuator_lookup_table[48] = 71;
    dm_actuator_lookup_table[49] = 72;
    dm_actuator_lookup_table[50] = 73;
    dm_actuator_lookup_table[51] = 74;
    dm_actuator_lookup_table[52] = 75;
    dm_actuator_lookup_table[53] = 76;
    dm_actuator_lookup_table[54] = 77;
    dm_actuator_lookup_table[55] = 78;
    dm_actuator_lookup_table[56] = 79;
    dm_actuator_lookup_table[57] = 80;
    dm_actuator_lookup_table[58] = 81;
    dm_actuator_lookup_table[59] = 82;
    dm_actuator_lookup_table[60] = 83;
    dm_actuator_lookup_table[61] = 85;
    dm_actuator_lookup_table[62] = 86;
    dm_actuator_lookup_table[63] = 87;
    dm_actuator_lookup_table[64] = 88;
    dm_actuator_lookup_table[65] = 89;
    dm_actuator_lookup_table[66] = 90;
    dm_actuator_lookup_table[67] = 91;
    dm_actuator_lookup_table[68] = 92;
    dm_actuator_lookup_table[69] = 93;
    dm_actuator_lookup_table[70] = 94;
    dm_actuator_lookup_table[71] = 95;
    dm_actuator_lookup_table[72] = 96;
    dm_actuator_lookup_table[73] = 97;
    dm_actuator_lookup_table[74] = 98;
    dm_actuator_lookup_table[75] = 99;
    dm_actuator_lookup_table[76] = 100;
    dm_actuator_lookup_table[77] = 101;
    dm_actuator_lookup_table[78] = 102;
    dm_actuator_lookup_table[79] = 103;
    dm_actuator_lookup_table[80] = 104;
    dm_actuator_lookup_table[81] = 105;
    dm_actuator_lookup_table[82] = 106;
    dm_actuator_lookup_table[83] = 107;
    dm_actuator_lookup_table[84] = 108;
    dm_actuator_lookup_table[85] = 109;
    dm_actuator_lookup_table[86] = 110;
    dm_actuator_lookup_table[87] = 111;
    dm_actuator_lookup_table[88] = 112;
    dm_actuator_lookup_table[89] = 113;
    dm_actuator_lookup_table[90] = 114;
    dm_actuator_lookup_table[91] = 115;
    dm_actuator_lookup_table[92] = 116;
    dm_actuator_lookup_table[93] = 117;
    dm_actuator_lookup_table[94] = 118;
    dm_actuator_lookup_table[95] = 119;
    dm_actuator_lookup_table[96] = 120;
    dm_actuator_lookup_table[97] = 121;
    dm_actuator_lookup_table[98] = 122;
    dm_actuator_lookup_table[99] = 123;
    dm_actuator_lookup_table[100] = 124;
    dm_actuator_lookup_table[101] = 125;
    dm_actuator_lookup_table[102] = 126;
    dm_actuator_lookup_table[103] = 127;
    dm_actuator_lookup_table[104] = 128;
    dm_actuator_lookup_table[105] = 129;
    dm_actuator_lookup_table[106] = 130;
    dm_actuator_lookup_table[107] = 131;
    dm_actuator_lookup_table[108] = 132;
    dm_actuator_lookup_table[109] = 133;
    dm_actuator_lookup_table[110] = 134;
    dm_actuator_lookup_table[111] = 135;
    dm_actuator_lookup_table[112] = 136;
    dm_actuator_lookup_table[113] = 137;
    dm_actuator_lookup_table[114] = 138;
    dm_actuator_lookup_table[115] = 139;
    dm_actuator_lookup_table[116] = 140;
    dm_actuator_lookup_table[117] = 141;
    dm_actuator_lookup_table[118] = 142;
    dm_actuator_lookup_table[119] = 143;
    dm_actuator_lookup_table[120] = 144;
    dm_actuator_lookup_table[121] = 145;
    dm_actuator_lookup_table[122] = 146;
    dm_actuator_lookup_table[123] = 147;
    dm_actuator_lookup_table[124] = 148;
    dm_actuator_lookup_table[125] = 149;
    dm_actuator_lookup_table[126] = 150;
    dm_actuator_lookup_table[127] = 151;
    dm_actuator_lookup_table[128] = 152;
    dm_actuator_lookup_table[129] = 153;
    dm_actuator_lookup_table[130] = 154;
    dm_actuator_lookup_table[131] = 155;
    dm_actuator_lookup_table[132] = 156;
    dm_actuator_lookup_table[133] = 157;
    dm_actuator_lookup_table[134] = 158;
    dm_actuator_lookup_table[135] = 159;
    dm_actuator_lookup_table[136] = 160;
    dm_actuator_lookup_table[137] = 161;
    dm_actuator_lookup_table[138] = 162;
    dm_actuator_lookup_table[139] = 163;
    dm_actuator_lookup_table[140] = 164;
    dm_actuator_lookup_table[141] = 165;
    dm_actuator_lookup_table[142] = 166;
    dm_actuator_lookup_table[143] = 167;
    dm_actuator_lookup_table[144] = 168;
    dm_actuator_lookup_table[145] = 169;
    dm_actuator_lookup_table[146] = 170;
    dm_actuator_lookup_table[147] = 171;
    dm_actuator_lookup_table[148] = 172;
    dm_actuator_lookup_table[149] = 173;
    dm_actuator_lookup_table[150] = 174;
    dm_actuator_lookup_table[151] = 175;
    dm_actuator_lookup_table[152] = 176;
    dm_actuator_lookup_table[153] = 177;
    dm_actuator_lookup_table[154] = 178;
    dm_actuator_lookup_table[155] = 179;
    dm_actuator_lookup_table[156] = 180;
    dm_actuator_lookup_table[157] = 181;
    dm_actuator_lookup_table[158] = 182;
    dm_actuator_lookup_table[159] = 183;
    dm_actuator_lookup_table[160] = 184;
    dm_actuator_lookup_table[161] = 185;
    dm_actuator_lookup_table[162] = 186;
    dm_actuator_lookup_table[163] = 187;
    dm_actuator_lookup_table[164] = 188;
    dm_actuator_lookup_table[165] = 189;
    dm_actuator_lookup_table[166] = 190;
    dm_actuator_lookup_table[167] = 191;
    dm_actuator_lookup_table[168] = 192;
    dm_actuator_lookup_table[169] = 193;
    dm_actuator_lookup_table[170] = 194;
    dm_actuator_lookup_table[171] = 195;
    dm_actuator_lookup_table[172] = 196;
    dm_actuator_lookup_table[173] = 197;
    dm_actuator_lookup_table[174] = 198;
    dm_actuator_lookup_table[175] = 199;
    dm_actuator_lookup_table[176] = 200;
    dm_actuator_lookup_table[177] = 201;
    dm_actuator_lookup_table[178] = 202;
    dm_actuator_lookup_table[179] = 203;
    dm_actuator_lookup_table[180] = 205;
    dm_actuator_lookup_table[181] = 206;
    dm_actuator_lookup_table[182] = 207;
    dm_actuator_lookup_table[183] = 208;
    dm_actuator_lookup_table[184] = 209;
    dm_actuator_lookup_table[185] = 210;
    dm_actuator_lookup_table[186] = 211;
    dm_actuator_lookup_table[187] = 212;
    dm_actuator_lookup_table[188] = 213;
    dm_actuator_lookup_table[189] = 214;
    dm_actuator_lookup_table[190] = 215;
    dm_actuator_lookup_table[191] = 216;
    dm_actuator_lookup_table[192] = 217;
    dm_actuator_lookup_table[193] = 218;
    dm_actuator_lookup_table[194] = 219;
    dm_actuator_lookup_table[195] = 222;
    dm_actuator_lookup_table[196] = 223;
    dm_actuator_lookup_table[197] = 224;
    dm_actuator_lookup_table[198] = 225;
    dm_actuator_lookup_table[199] = 226;
    dm_actuator_lookup_table[200] = 227;
    dm_actuator_lookup_table[201] = 228;
    dm_actuator_lookup_table[202] = 229;
    dm_actuator_lookup_table[203] = 230;
    dm_actuator_lookup_table[204] = 231;
    dm_actuator_lookup_table[205] = 232;
    dm_actuator_lookup_table[206] = 233;
    dm_actuator_lookup_table[207] = 234;
    dm_actuator_lookup_table[208] = 235;
    dm_actuator_lookup_table[209] = 236;
    dm_actuator_lookup_table[210] = 240;
    dm_actuator_lookup_table[211] = 241;
    dm_actuator_lookup_table[212] = 242;
    dm_actuator_lookup_table[213] = 243;
    dm_actuator_lookup_table[214] = 244;
    dm_actuator_lookup_table[215] = 245;
    dm_actuator_lookup_table[216] = 246;
    dm_actuator_lookup_table[217] = 247;
    dm_actuator_lookup_table[218] = 248;
    dm_actuator_lookup_table[219] = 249;
    dm_actuator_lookup_table[220] = 250;
    dm_actuator_lookup_table[221] = 251;
    dm_actuator_lookup_table[222] = 252;
    dm_actuator_lookup_table[223] = 258;
    dm_actuator_lookup_table[224] = 259;
    dm_actuator_lookup_table[225] = 260;
    dm_actuator_lookup_table[226] = 261;
    dm_actuator_lookup_table[227] = 262;
    dm_actuator_lookup_table[228] = 263;
    dm_actuator_lookup_table[229] = 264;
    dm_actuator_lookup_table[230] = 265;
    dm_actuator_lookup_table[231] = 266;
    dm_actuator_lookup_table[232] = 267;
    dm_actuator_lookup_table[233] = 268;
    dm_actuator_lookup_table[234] = 277;
    dm_actuator_lookup_table[235] = 278;
    dm_actuator_lookup_table[236] = 279;
    dm_actuator_lookup_table[237] = 280;
    dm_actuator_lookup_table[238] = 281;
    dm_actuator_lookup_table[239] = 282;
    dm_actuator_lookup_table[240] = 283;    
  }

  void PALAO_reconstructor::create_centroid_weighting_lookup_table() const {
    centroid_weighting_lookup_table.resize(256);

    bool zero_inner_annular_centroids = true;

    centroid_weighting_lookup_table[0] = 0;
    centroid_weighting_lookup_table[1] = 0;
    centroid_weighting_lookup_table[2] = 0;
    centroid_weighting_lookup_table[3] = 0;
    centroid_weighting_lookup_table[4] = 0;
    centroid_weighting_lookup_table[5] = 1;
    centroid_weighting_lookup_table[6] = 1;
    centroid_weighting_lookup_table[7] = 1;
    centroid_weighting_lookup_table[8] = 1;
    centroid_weighting_lookup_table[9] = 1;
    centroid_weighting_lookup_table[10] = 1;
    centroid_weighting_lookup_table[11] = 0;
    centroid_weighting_lookup_table[12] = 0;
    centroid_weighting_lookup_table[13] = 0;
    centroid_weighting_lookup_table[14] = 0;
    centroid_weighting_lookup_table[15] = 0;
    centroid_weighting_lookup_table[16] = 0;
    centroid_weighting_lookup_table[17] = 0;
    centroid_weighting_lookup_table[18] = 0;
    centroid_weighting_lookup_table[19] = 1;
    centroid_weighting_lookup_table[20] = 1;
    centroid_weighting_lookup_table[21] = 1;
    centroid_weighting_lookup_table[22] = 1;
    centroid_weighting_lookup_table[23] = 1;
    centroid_weighting_lookup_table[24] = 1;
    centroid_weighting_lookup_table[25] = 1;
    centroid_weighting_lookup_table[26] = 1;
    centroid_weighting_lookup_table[27] = 1;
    centroid_weighting_lookup_table[28] = 1;
    centroid_weighting_lookup_table[29] = 0;
    centroid_weighting_lookup_table[30] = 0;
    centroid_weighting_lookup_table[31] = 0;
    centroid_weighting_lookup_table[32] = 0;
    centroid_weighting_lookup_table[33] = 0;
    centroid_weighting_lookup_table[34] = 1;
    centroid_weighting_lookup_table[35] = 1;
    centroid_weighting_lookup_table[36] = 1;
    centroid_weighting_lookup_table[37] = 1;
    centroid_weighting_lookup_table[38] = 1;
    centroid_weighting_lookup_table[39] = 1;
    centroid_weighting_lookup_table[40] = 1;
    centroid_weighting_lookup_table[41] = 1;
    centroid_weighting_lookup_table[42] = 1;
    centroid_weighting_lookup_table[43] = 1;
    centroid_weighting_lookup_table[44] = 1;
    centroid_weighting_lookup_table[45] = 1;
    centroid_weighting_lookup_table[46] = 0;
    centroid_weighting_lookup_table[47] = 0;
    centroid_weighting_lookup_table[48] = 0;
    centroid_weighting_lookup_table[49] = 1;
    centroid_weighting_lookup_table[50] = 1;
    centroid_weighting_lookup_table[51] = 1;
    centroid_weighting_lookup_table[52] = 1;
    centroid_weighting_lookup_table[53] = 1;
    centroid_weighting_lookup_table[54] = 1;
    centroid_weighting_lookup_table[55] = 1;
    centroid_weighting_lookup_table[56] = 1;
    centroid_weighting_lookup_table[57] = 1;
    centroid_weighting_lookup_table[58] = 1;
    centroid_weighting_lookup_table[59] = 1;
    centroid_weighting_lookup_table[60] = 1;
    centroid_weighting_lookup_table[61] = 1;
    centroid_weighting_lookup_table[62] = 1;
    centroid_weighting_lookup_table[63] = 0;
    centroid_weighting_lookup_table[64] = 0;
    centroid_weighting_lookup_table[65] = 1;
    centroid_weighting_lookup_table[66] = 1;
    centroid_weighting_lookup_table[67] = 1;
    centroid_weighting_lookup_table[68] = 1;
    centroid_weighting_lookup_table[69] = 1;
    centroid_weighting_lookup_table[70] = 1;
    centroid_weighting_lookup_table[71] = 1;
    centroid_weighting_lookup_table[72] = 1;
    centroid_weighting_lookup_table[73] = 1;
    centroid_weighting_lookup_table[74] = 1;
    centroid_weighting_lookup_table[75] = 1;
    centroid_weighting_lookup_table[76] = 1;
    centroid_weighting_lookup_table[77] = 1;
    centroid_weighting_lookup_table[78] = 1;
    centroid_weighting_lookup_table[79] = 0;
    centroid_weighting_lookup_table[80] = 1;
    centroid_weighting_lookup_table[81] = 1;
    centroid_weighting_lookup_table[82] = 1;
    centroid_weighting_lookup_table[83] = 1;
    centroid_weighting_lookup_table[84] = 1;
    centroid_weighting_lookup_table[85] = 1;
    centroid_weighting_lookup_table[86] = zero_inner_annular_centroids ? 0 : 1;
    centroid_weighting_lookup_table[87] = zero_inner_annular_centroids ? 0 : 1;
    centroid_weighting_lookup_table[88] = zero_inner_annular_centroids ? 0 : 1;
    centroid_weighting_lookup_table[89] = zero_inner_annular_centroids ? 0 : 1;
    centroid_weighting_lookup_table[90] = 1;
    centroid_weighting_lookup_table[91] = 1;
    centroid_weighting_lookup_table[92] = 1;
    centroid_weighting_lookup_table[93] = 1;
    centroid_weighting_lookup_table[94] = 1;
    centroid_weighting_lookup_table[95] = 1;
    centroid_weighting_lookup_table[96] = 1;
    centroid_weighting_lookup_table[97] = 1;
    centroid_weighting_lookup_table[98] = 1;
    centroid_weighting_lookup_table[99] = 1;
    centroid_weighting_lookup_table[100] = 1;
    centroid_weighting_lookup_table[101] = zero_inner_annular_centroids ? 0 : 1;
    centroid_weighting_lookup_table[102] = zero_inner_annular_centroids ? 0 : 1;
    centroid_weighting_lookup_table[103] = zero_inner_annular_centroids ? 0 : 1;
    centroid_weighting_lookup_table[104] = zero_inner_annular_centroids ? 0 : 1;
    centroid_weighting_lookup_table[105] = zero_inner_annular_centroids ? 0 : 1;
    centroid_weighting_lookup_table[106] = zero_inner_annular_centroids ? 0 : 1;
    centroid_weighting_lookup_table[107] = 1;
    centroid_weighting_lookup_table[108] = 1;
    centroid_weighting_lookup_table[109] = 1;
    centroid_weighting_lookup_table[110] = 1;
    centroid_weighting_lookup_table[111] = 1;
    centroid_weighting_lookup_table[112] = 1;
    centroid_weighting_lookup_table[113] = 1;
    centroid_weighting_lookup_table[114] = 1;
    centroid_weighting_lookup_table[115] = 1;
    centroid_weighting_lookup_table[116] = 1;
    centroid_weighting_lookup_table[117] = zero_inner_annular_centroids ? 0 : 1;
    centroid_weighting_lookup_table[118] = zero_inner_annular_centroids ? 0 : 1;
    centroid_weighting_lookup_table[119] = zero_inner_annular_centroids ? 0 : 1;
    centroid_weighting_lookup_table[120] = zero_inner_annular_centroids ? 0 : 1;
    centroid_weighting_lookup_table[121] = zero_inner_annular_centroids ? 0 : 1;
    centroid_weighting_lookup_table[122] = zero_inner_annular_centroids ? 0 : 1;
    centroid_weighting_lookup_table[123] = 1;
    centroid_weighting_lookup_table[124] = 1;
    centroid_weighting_lookup_table[125] = 1;
    centroid_weighting_lookup_table[126] = 1;
    centroid_weighting_lookup_table[127] = 1;
    centroid_weighting_lookup_table[128] = 1;
    centroid_weighting_lookup_table[129] = 1;
    centroid_weighting_lookup_table[130] = 1;
    centroid_weighting_lookup_table[131] = 1;
    centroid_weighting_lookup_table[132] = 1;
    centroid_weighting_lookup_table[133] = zero_inner_annular_centroids ? 0 : 1;
    centroid_weighting_lookup_table[134] = zero_inner_annular_centroids ? 0 : 1;
    centroid_weighting_lookup_table[135] = zero_inner_annular_centroids ? 0 : 1;
    centroid_weighting_lookup_table[136] = zero_inner_annular_centroids ? 0 : 1;
    centroid_weighting_lookup_table[137] = zero_inner_annular_centroids ? 0 : 1;
    centroid_weighting_lookup_table[138] = zero_inner_annular_centroids ? 0 : 1;
    centroid_weighting_lookup_table[139] = 1;
    centroid_weighting_lookup_table[140] = 1;
    centroid_weighting_lookup_table[141] = 1;
    centroid_weighting_lookup_table[142] = 1;
    centroid_weighting_lookup_table[143] = 1;
    centroid_weighting_lookup_table[144] = 1;
    centroid_weighting_lookup_table[145] = 1;
    centroid_weighting_lookup_table[146] = 1;
    centroid_weighting_lookup_table[147] = 1;
    centroid_weighting_lookup_table[148] = 1;
    centroid_weighting_lookup_table[149] = zero_inner_annular_centroids ? 0 : 1;
    centroid_weighting_lookup_table[150] = zero_inner_annular_centroids ? 0 : 1;
    centroid_weighting_lookup_table[151] = zero_inner_annular_centroids ? 0 : 1;
    centroid_weighting_lookup_table[152] = zero_inner_annular_centroids ? 0 : 1;
    centroid_weighting_lookup_table[153] = zero_inner_annular_centroids ? 0 : 1;
    centroid_weighting_lookup_table[154] = zero_inner_annular_centroids ? 0 : 1;
    centroid_weighting_lookup_table[155] = 1;
    centroid_weighting_lookup_table[156] = 1;
    centroid_weighting_lookup_table[157] = 1;
    centroid_weighting_lookup_table[158] = 1;
    centroid_weighting_lookup_table[159] = 1;
    centroid_weighting_lookup_table[160] = 1;
    centroid_weighting_lookup_table[161] = 1;
    centroid_weighting_lookup_table[162] = 1;
    centroid_weighting_lookup_table[163] = 1;
    centroid_weighting_lookup_table[164] = 1;
    centroid_weighting_lookup_table[165] = 1;
    centroid_weighting_lookup_table[166] = zero_inner_annular_centroids ? 0 : 1;
    centroid_weighting_lookup_table[167] = zero_inner_annular_centroids ? 0 : 1;
    centroid_weighting_lookup_table[168] = zero_inner_annular_centroids ? 0 : 1;
    centroid_weighting_lookup_table[169] = zero_inner_annular_centroids ? 0 : 1;
    centroid_weighting_lookup_table[170] = 1;
    centroid_weighting_lookup_table[171] = 1;
    centroid_weighting_lookup_table[172] = 1;
    centroid_weighting_lookup_table[173] = 1;
    centroid_weighting_lookup_table[174] = 1;
    centroid_weighting_lookup_table[175] = 1;
    centroid_weighting_lookup_table[176] = 0;
    centroid_weighting_lookup_table[177] = 1;
    centroid_weighting_lookup_table[178] = 1;
    centroid_weighting_lookup_table[179] = 1;
    centroid_weighting_lookup_table[180] = 1;
    centroid_weighting_lookup_table[181] = 1;
    centroid_weighting_lookup_table[182] = 1;
    centroid_weighting_lookup_table[183] = 1;
    centroid_weighting_lookup_table[184] = 1;
    centroid_weighting_lookup_table[185] = 1;
    centroid_weighting_lookup_table[186] = 1;
    centroid_weighting_lookup_table[187] = 1;
    centroid_weighting_lookup_table[188] = 1;
    centroid_weighting_lookup_table[189] = 1;
    centroid_weighting_lookup_table[190] = 1;
    centroid_weighting_lookup_table[191] = 0;
    centroid_weighting_lookup_table[192] = 0;
    centroid_weighting_lookup_table[193] = 1;
    centroid_weighting_lookup_table[194] = 1;
    centroid_weighting_lookup_table[195] = 1;
    centroid_weighting_lookup_table[196] = 1;
    centroid_weighting_lookup_table[197] = 1;
    centroid_weighting_lookup_table[198] = 1;
    centroid_weighting_lookup_table[199] = 1;
    centroid_weighting_lookup_table[200] = 1;
    centroid_weighting_lookup_table[201] = 1;
    centroid_weighting_lookup_table[202] = 1;
    centroid_weighting_lookup_table[203] = 1;
    centroid_weighting_lookup_table[204] = 1;
    centroid_weighting_lookup_table[205] = 1;
    centroid_weighting_lookup_table[206] = 1;
    centroid_weighting_lookup_table[207] = 0;
    centroid_weighting_lookup_table[208] = 0;
    centroid_weighting_lookup_table[209] = 0;
    centroid_weighting_lookup_table[210] = 1;
    centroid_weighting_lookup_table[211] = 1;
    centroid_weighting_lookup_table[212] = 1;
    centroid_weighting_lookup_table[213] = 1;
    centroid_weighting_lookup_table[214] = 1;
    centroid_weighting_lookup_table[215] = 1;
    centroid_weighting_lookup_table[216] = 1;
    centroid_weighting_lookup_table[217] = 1;
    centroid_weighting_lookup_table[218] = 1;
    centroid_weighting_lookup_table[219] = 1;
    centroid_weighting_lookup_table[220] = 1;
    centroid_weighting_lookup_table[221] = 1;
    centroid_weighting_lookup_table[222] = 0;
    centroid_weighting_lookup_table[223] = 0;
    centroid_weighting_lookup_table[224] = 0;
    centroid_weighting_lookup_table[225] = 0;
    centroid_weighting_lookup_table[226] = 0;
    centroid_weighting_lookup_table[227] = 1;
    centroid_weighting_lookup_table[228] = 1;
    centroid_weighting_lookup_table[229] = 1;
    centroid_weighting_lookup_table[230] = 1;
    centroid_weighting_lookup_table[231] = 1;
    centroid_weighting_lookup_table[232] = 1;
    centroid_weighting_lookup_table[233] = 1;
    centroid_weighting_lookup_table[234] = 1;
    centroid_weighting_lookup_table[235] = 1;
    centroid_weighting_lookup_table[236] = 1;
    centroid_weighting_lookup_table[237] = 0;
    centroid_weighting_lookup_table[238] = 0;
    centroid_weighting_lookup_table[239] = 0;
    centroid_weighting_lookup_table[240] = 0;
    centroid_weighting_lookup_table[241] = 0;
    centroid_weighting_lookup_table[242] = 0;
    centroid_weighting_lookup_table[243] = 0;
    centroid_weighting_lookup_table[244] = 0;
    centroid_weighting_lookup_table[245] = 1;
    centroid_weighting_lookup_table[246] = 1;
    centroid_weighting_lookup_table[247] = 1;
    centroid_weighting_lookup_table[248] = 1;
    centroid_weighting_lookup_table[249] = 1;
    centroid_weighting_lookup_table[250] = 1;
    centroid_weighting_lookup_table[251] = 0;
    centroid_weighting_lookup_table[252] = 0;
    centroid_weighting_lookup_table[253] = 0;
    centroid_weighting_lookup_table[254] = 0;
    centroid_weighting_lookup_table[255] = 0;

  }

}
