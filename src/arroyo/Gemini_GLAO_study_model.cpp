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

#include "fits_factory.h"
#include "Gemini_GLAO_study_model.h"

  /*
GL_bad_FA_bad 
0	1.38e-13	0.26	8.2	243
25	1.08e-13	0.31	7.0	235
50	1.53e-13	0.25	6.4	225
100	1.58e-13	0.24	8.2	216
200	1.03e-13	0.31	9.4	210
400	6.46e-14	0.42	10.5	230
800	7.29e-14	0.39	10.4	233
1600	6.77e-14	0.40	11.8	238
3600	3.2e-13		0.16	19.3	287

0       1.38e-13       0.26    8.2     243
25      1.08e-13       0.31    7.0     235
50      1.53e-13       0.25    6.4     225
100     1.58e-13       0.24    8.2     216
200     1.03e-13       0.31    9.4     210
400     6.46e-14       0.42    10.5    230
800     7.29e-14       0.39    10.4    233
1600    6.77e-14       0.40    11.8    238
3400    2.00e-13        0.21    7.2     270
6000    7.95e-14        0.37    16.5    269
7600    2.84e-14        0.68    23.2    259
13300   7.70e-15        1.49    32.7    259
16000   4.81e-15        1.98    5.7     320

GL_bad_FA_typical
0	1.38e-13	0.26	8.2	243
25	1.08e-13	0.31	7.0	235
50	1.53e-13	0.25	6.4	225
100	1.58e-13	0.24	8.2	216
200	1.03e-13	0.31	9.4	210
400	6.46e-14	0.42	10.5	230
800	7.29e-14	0.39	10.4	233
1600	6.77e-14	0.40	11.8	238
5500	1.7e-13		0.23	27.0	283

0       1.38e-13       0.26    8.2     243
25      1.08e-13       0.31    7.0     235
50      1.53e-13       0.25    6.4     225
100     1.58e-13       0.24    8.2     216
200     1.03e-13       0.31    9.4     210
400     6.46e-14       0.42    10.5    230
800     7.29e-14       0.39    10.4    233
1600    6.77e-14       0.40    11.8    238
3400    6.97e-14        0.40    7.2     270
6000    5.47e-14        0.50    16.5    269
7600    1.95e-14        0.85    23.2    259
13300   1.30e-14        1.09    32.7    259
16000   1.31e-14        1.08    5.7     320

GL_bad_FA_good 
0	1.38e-13	0.26	8.2	243
25	1.08e-13	0.31	7.0	235
50	1.53e-13	0.25	6.4	225
100	1.58e-13	0.24	8.2	216
200	1.03e-13	0.31	9.4	210
400	6.46e-14	0.42	10.5	230
800	7.29e-14	0.39	10.4	233
1600	6.77e-14	0.40	11.8	238
8400	9.0e-14		0.34	39.3	270

0       1.38e-13       0.26    8.2     243
25      1.08e-13       0.31    7.0     235
50      1.53e-13       0.25    6.4     225
100     1.58e-13       0.24    8.2     216
200     1.03e-13       0.31    9.4     210
400     6.46e-14       0.42    10.5    230
800     7.29e-14       0.39    10.4    233
1600    6.77e-14       0.40    11.8    238
3400    2.95e-14       0.67    7.2     270
6000    1.44e-14       1.02    16.5    269
7600    1.01e-14       1.27    23.2    259
13300   2.98e-14       0.66    32.7    259
16000   6.15e-15       1.71    5.7     320

GL_typical_FA_bad 
0	7.04e-14	0.40	6.9	284
25	2.25e-14	0.78	7.5	267
50	1.35e-14	1.07	7.8	244
100	1.24e-14	1.12	8.3	267
200	1.99e-14	0.84	9.6	237
400	2.87e-14	0.68	9.9	232
800	3.02e-14	0.66	9.6	286
1600	1.75e-14	0.91	10.1	293
3600	3.2e-13		0.16	19.4	275

0       7.04e-14       0.40    6.9     284
25      2.25e-14       0.78    7.5     267
50      1.35e-14       1.07    7.8     244
100     1.24e-14       1.12    8.3     267
200     1.99e-14       0.84    9.6     237
400     2.87e-14       0.68    9.9     232
800     3.02e-14       0.66    9.6     286
1600    1.75e-14       0.91    10.1    293
3400    2.00e-13        0.21    7.2     270
6000    7.95e-14        0.37    16.5    269
7600    2.84e-14        0.68    23.2    259
13300   7.70e-15        1.49    32.7    259
16000   4.81e-15        1.98    5.7     320
 
GL_typical_FA_typical
0	7.04e-14	0.40	6.9	284
25	2.25e-14	0.78	7.5	267
50	1.35e-14	1.07	7.8	244
100	1.24e-14	1.12	8.3	267
200	1.99e-14	0.84	9.6	237
400	2.87e-14	0.68	9.9	232
800	3.02e-14	0.66	9.6	286
1600	1.75e-14	0.91	10.1	293
5500	1.7e-13		0.23	31.5	274

0       7.04e-14       0.40    6.9     284
25      2.25e-14       0.78    7.5     267
50      1.35e-14       1.07    7.8     244
100     1.24e-14       1.12    8.3     267
200     1.99e-14       0.84    9.6     237
400     2.87e-14       0.68    9.9     232
800     3.02e-14       0.66    9.6     286
1600    1.75e-14       0.91    10.1    293
3400    6.97e-14        0.40    7.2     270
6000    5.47e-14        0.50    16.5    269
7600    1.95e-14        0.85    23.2    259
13300   1.30e-14        1.09    32.7    259
16000   1.31e-14        1.08    5.7     320
 
GL_typical_FA_good
0	7.04e-14	0.40	6.9	284
25	2.25e-14	0.78	7.5	267
50	1.35e-14	1.07	7.8	244
100	1.24e-14	1.12	8.3	267
200	1.99e-14	0.84	9.6	237
400	2.87e-14	0.68	9.9	232
800	3.02e-14	0.66	9.6	286
1600	1.75e-14	0.91	10.1	293
8400	9.0e-14		0.34	41.9	274

0       7.04e-14       0.40    6.9     284
25      2.25e-14       0.78    7.5     267
50      1.35e-14       1.07    7.8     244
100     1.24e-14       1.12    8.3     267
200     1.99e-14       0.84    9.6     237
400     2.87e-14       0.68    9.9     232
800     3.02e-14       0.66    9.6     286
1600    1.75e-14       0.91    10.1    293
3400    2.95e-14       0.67    7.2     270
6000    1.44e-14       1.02    16.5    269
7600    1.01e-14       1.27    23.2    259
13300   2.98e-14       0.66    32.7    259
16000   6.15e-15       1.71    5.7     320
 

GL_good_FA_bad
0	9.26e-14	0.34	1.7	230
25	1.83e-14	0.89	2.8	242
50	5.74e-15	1.78	4.0	256
100	3.62e-15	2.34	5.5	170
200	6.14e-15	1.71	5.6	162
400	9.60e-15	1.31	5.1	198
800	1.18e-14	1.15	5.4	320
1600	9.13e-15	1.35	7.1	310
3600	3.2e-13		0.16	12.8	278

0       9.26e-14       0.34    1.7     230
25      1.83e-14       0.89    2.8     242
50      5.74e-15       1.78    4.0     256
100     3.62e-15       2.34    5.5     170
200     6.14e-15       1.71    5.6     162
400     9.60e-15       1.31    5.1     198
800     1.18e-14       1.15    5.4     320
1600    9.13e-15       1.35    7.1     310
3400    2.00e-13        0.21    7.2     270
6000    7.95e-14        0.37    16.5    269
7600    2.84e-14        0.68    23.2    259
13300   7.70e-15        1.49    32.7    259
16000   4.81e-15        1.98    5.7     320

GL_good_FA_typical
0	9.26e-14	0.34	1.7	230
25	1.83e-14	0.89	2.8	242
50	5.74e-15	1.78	4.0	256
100	3.62e-15	2.34	5.5	170
200	6.14e-15	1.71	5.6	162
400	9.60e-15	1.31	5.1	198
800	1.18e-14	1.15	5.4	320
1600	9.13e-15	1.35	7.1	310
5500	1.7e-13		0.23	19.7	262

0       9.26e-14       0.34    1.7     230
25      1.83e-14       0.89    2.8     242
50      5.74e-15       1.78    4.0     256
100     3.62e-15       2.34    5.5     170
200     6.14e-15       1.71    5.6     162
400     9.60e-15       1.31    5.1     198
800     1.18e-14       1.15    5.4     320
1600    9.13e-15       1.35    7.1     310
3400    6.97e-14        0.40    7.2     270
6000    5.47e-14        0.50    16.5    269
7600    1.95e-14        0.85    23.2    259
13300   1.30e-14        1.09    32.7    259
16000   1.31e-14        1.08    5.7     320
 
GL_good_FA_good
0	9.26e-14	0.34	1.7	230
25	1.83e-14	0.89	2.8	242
50	5.74e-15	1.78	4.0	256
100	3.62e-15	2.34	5.5	170
200	6.14e-15	1.71	5.6	162
400	9.60e-15	1.31	5.1	198
800	1.18e-14	1.15	5.4	320
1600	9.13e-15	1.35	7.1	310
8400	9.0e-14		0.34	30.6	259

0       9.26e-14       0.34    1.7     230
25      1.83e-14       0.89    2.8     242
50      5.74e-15       1.78    4.0     256
100     3.62e-15       2.34    5.5     170
200     6.14e-15       1.71    5.6     162
400     9.60e-15       1.31    5.1     198
800     1.18e-14       1.15    5.4     320
1600    9.13e-15       1.35    7.1     310
3400    2.95e-14       0.67    7.2     270
6000    1.44e-14       1.02    16.5    269
7600    1.01e-14       1.27    23.2    259
13300   2.98e-14       0.66    32.7    259
16000   6.15e-15       1.71    5.7     320
 
*/



using namespace std;

namespace Arroyo {

  namespace factory_register {

    const fits_keyval_set & get_Gemini_GLAO_study_model_keyval_set(){
      static fits_keyval_set * fkvs = new fits_keyval_set;
      fkvs->push_back(fits_keyval_entry("TYPE", "Gemini GLAO study model"));
      return *fkvs;
    }
    
    AO_sim_base * create_Gemini_GLAO_study_model(const iofits & iof) {
      return new Gemini_GLAO_study_model(iof);
    }
  }

  const bool Gemini_GLAO_study_model::factory_registration = 
  fits_factory<AO_sim_base>::Register(factory_register::get_Gemini_GLAO_study_model_keyval_set(), 
				      factory_register::create_Gemini_GLAO_study_model);


  Gemini_GLAO_study_model::Gemini_GLAO_study_model(const Gemini_GLAO_study_model & cn2_model){
    this->operator=(cn2_model);
  }

  Gemini_GLAO_study_model::Gemini_GLAO_study_model(const char * filename){
    this->read(filename);
  }

  Gemini_GLAO_study_model::Gemini_GLAO_study_model(const iofits & iof){
    this->read(iof);
  }
  
  Gemini_GLAO_study_model::Gemini_GLAO_study_model(const three_frame & ground_ref_frame, 
						   const string & ground_layer,
						   const string & upper_layers,
						   bool extended_profile,
						   double outer_scale_meters){

    if(ground_layer!="good" &&
       ground_layer!="typical" &&
       ground_layer!="bad"){
      cerr << "Gemini_GLAO_study_model::Gemini_GLAO_study_model error - "
	   << "cannot parse ground layer string " 
	   << ground_layer
	   << " as one of good, typical or bad\n";
      throw(string("Gemini_GLAO_study_model::Gemini_GLAO_study_model"));
    }

    if(upper_layers!="good" &&
       upper_layers!="typical" &&
       upper_layers!="bad"){
      cerr << "Gemini_GLAO_study_model::Gemini_GLAO_study_model error - "
	   << "cannot parse focal anisoplanatism string " 
	   << upper_layers
	   << " as one of good, typical or bad\n";
      throw(string("Gemini_GLAO_study_model::Gemini_GLAO_study_model"));
    }

    if(outer_scale_meters<=0 && outer_scale_meters!=-1){
      cerr << "Gemini_GLAO_study_model::Gemini_GLAO_study_model error - "
	   << "invalid outer scale "
	   << outer_scale_meters
	   << endl;
      throw(string("Gemini_GLAO_study_model::Gemini_GLAO_study_model"));
    }
    
    ground_ref_frame_ = ground_ref_frame;

    int index;

    if(extended_profile){
      layer_heights_.resize(13);
      index = 5;
      layer_heights_[4] = 3400;
      layer_heights_[3] = 6000;
      layer_heights_[2] = 7600;
      layer_heights_[1] = 13300;
      layer_heights_[0] = 16000;
    } else {
      layer_heights_.resize(9);
      index = 1;
      if(upper_layers=="bad")
	layer_heights_[0] = 3600;
      else if(upper_layers=="typical")
	layer_heights_[0] = 5500;
      else if(upper_layers=="good")
	layer_heights_[0] = 8400;
    }
    layer_heights_[index+7] = 0;
    layer_heights_[index+6] = 25;
    layer_heights_[index+5] = 50;
    layer_heights_[index+4] = 100;
    layer_heights_[index+3] = 200;
    layer_heights_[index+2] = 400;
    layer_heights_[index+1] = 800;
    layer_heights_[index  ] = 1600;


    vector<double> layer_r0_meters(layer_heights_.size());

    if(extended_profile){
      index = 5;
      if(ground_layer=="bad" && 
	 upper_layers=="bad"){	
	layer_r0_meters[4] = 0.21;
	layer_r0_meters[3] = 0.37;
	layer_r0_meters[2] = 0.68;
	layer_r0_meters[1] = 1.49;
	layer_r0_meters[0] = 1.98;
      } else if(ground_layer=="bad" && 
		upper_layers=="typical"){
	layer_r0_meters[4] = 0.40;
	layer_r0_meters[3] = 0.50;
	layer_r0_meters[2] = 0.85;
	layer_r0_meters[1] = 1.09;
	layer_r0_meters[0] = 1.08;
      } else if(ground_layer=="bad" && 
		upper_layers=="good"){
	layer_r0_meters[4] = 0.67;
	layer_r0_meters[3] = 1.02;
	layer_r0_meters[2] = 1.27;
	layer_r0_meters[1] = 0.66;
	layer_r0_meters[0] = 1.71;
      } else if(ground_layer=="typical" && 
	 upper_layers=="bad"){
	layer_r0_meters[4] = 0.21;
	layer_r0_meters[3] = 0.37;
	layer_r0_meters[2] = 0.68;
	layer_r0_meters[1] = 1.49;
	layer_r0_meters[0] = 1.98;
      } else if(ground_layer=="typical" && 
		upper_layers=="typical"){
	layer_r0_meters[4] = 0.40;
	layer_r0_meters[3] = 0.50;
	layer_r0_meters[2] = 0.85;
	layer_r0_meters[1] = 1.09;
	layer_r0_meters[0] = 1.08;
      } else if(ground_layer=="typical" && 
		upper_layers=="good"){
	layer_r0_meters[4] = 0.67;
	layer_r0_meters[3] = 1.02;
	layer_r0_meters[2] = 1.27;
	layer_r0_meters[1] = 0.66;
	layer_r0_meters[0] = 1.71;
      } else if(ground_layer=="good" && 
	 upper_layers=="bad"){
	layer_r0_meters[4] = 0.21;
	layer_r0_meters[3] = 0.37;
	layer_r0_meters[2] = 0.68;
	layer_r0_meters[1] = 1.49;
	layer_r0_meters[0] = 1.98;
      } else if(ground_layer=="good" && 
		upper_layers=="typical"){
	layer_r0_meters[4] = 0.40;
	layer_r0_meters[3] = 0.50;
	layer_r0_meters[2] = 0.85;
	layer_r0_meters[1] = 1.09;
	layer_r0_meters[0] = 1.08;
      } else if(ground_layer=="good" && 
		upper_layers=="good"){
	layer_r0_meters[4] = 0.67;
	layer_r0_meters[3] = 1.02;
	layer_r0_meters[2] = 1.27;
	layer_r0_meters[1] = 0.66;
	layer_r0_meters[0] = 1.71;
      }
    } else {
      index = 1;
      if(upper_layers=="good")
	layer_r0_meters[0] = 0.34;
      else if(upper_layers=="typical")
	layer_r0_meters[0] = 0.23;
      else if(upper_layers=="bad")
	layer_r0_meters[0] = 0.16;
    }

    if(ground_layer=="bad"){
      layer_r0_meters[index+7] = 0.26;
      layer_r0_meters[index+6] = 0.31;
      layer_r0_meters[index+5] = 0.25;
      layer_r0_meters[index+4] = 0.24;
      layer_r0_meters[index+3] = 0.31;
      layer_r0_meters[index+2] = 0.42;
      layer_r0_meters[index+1] = 0.39;
      layer_r0_meters[index  ] = 0.40;
    } else if(ground_layer=="typical"){
      layer_r0_meters[index+7] = 0.40;
      layer_r0_meters[index+6] = 0.78;
      layer_r0_meters[index+5] = 1.07;
      layer_r0_meters[index+4] = 1.12;
      layer_r0_meters[index+3] = 0.84;
      layer_r0_meters[index+2] = 0.68;
      layer_r0_meters[index+1] = 0.66;
      layer_r0_meters[index  ] = 0.91;
    } else if(ground_layer=="good"){
      layer_r0_meters[index+7] = 0.34;
      layer_r0_meters[index+6] = 0.89;
      layer_r0_meters[index+5] = 1.78;
      layer_r0_meters[index+4] = 2.34;
      layer_r0_meters[index+3] = 1.71;
      layer_r0_meters[index+2] = 1.31;
      layer_r0_meters[index+1] = 1.15;
      layer_r0_meters[index  ] = 1.35;
    }

    /*
    // Initialize layer heights
    layer_heights_.resize(9);
    layer_heights_[8] = 0;
    layer_heights_[7] = 25;
    layer_heights_[6] = 50;
    layer_heights_[5] = 100;
    layer_heights_[4] = 200;
    layer_heights_[3] = 400;
    layer_heights_[2] = 800;
    layer_heights_[1] = 1600;

    if(upper_layers=="good")
      layer_heights_[0] = 3600;
    else if(upper_layers=="typical"){
      if(extended_profile){
	layer_heights_.resize(13);
	layer_heights_[12] = 0;
	layer_heights_[11] = 25;
	layer_heights_[10] = 50;
	layer_heights_[9] = 100;
	layer_heights_[8] = 200;
	layer_heights_[7] = 400;
	layer_heights_[6] = 800;
	layer_heights_[5] = 1600;
	layer_heights_[4] = 3400;
	layer_heights_[3] = 6000;
	layer_heights_[2] = 7600;
	layer_heights_[1] = 13300;
	layer_heights_[0] = 16000;
      } else {
	layer_heights_[0] = 5500;
      }
    } else if(upper_layers=="bad")
      layer_heights_[0] = 8400;

    vector<double> layer_r0_meters(layer_heights_.size());

    if(ground_layer=="good"){
      layer_r0_meters[8] = 0.34;
      layer_r0_meters[7] = 0.89;
      layer_r0_meters[6] = 1.78;
      layer_r0_meters[5] = 2.34;
      layer_r0_meters[4] = 1.71;
      layer_r0_meters[3] = 1.31;
      layer_r0_meters[2] = 1.15;
      layer_r0_meters[1] = 1.35;
    } else if(ground_layer=="typical"){
      layer_r0_meters[8] = 0.40;
      layer_r0_meters[7] = 0.78;
      layer_r0_meters[6] = 1.07;
      layer_r0_meters[5] = 1.12;
      layer_r0_meters[4] = 0.84;
      layer_r0_meters[3] = 0.68;
      layer_r0_meters[2] = 0.66;
      layer_r0_meters[1] = 0.91;
    } else if(ground_layer=="bad"){
      layer_r0_meters[8] = 0.26;
      layer_r0_meters[7] = 0.31;
      layer_r0_meters[6] = 0.25;
      layer_r0_meters[5] = 0.24;
      layer_r0_meters[4] = 0.31;
      layer_r0_meters[3] = 0.42;
      layer_r0_meters[2] = 0.39;
      layer_r0_meters[1] = 0.40;
    }			    

    if(upper_layers=="good")
      layer_r0_meters[0] = 0.34;
    else if(upper_layers=="typical"){
      if(extended_profile){
	layer_r0_meters.resize(13);
	layer_r0_meters[12] = 0.40;
	layer_r0_meters[11] = 0.78;
	layer_r0_meters[10] = 1.07;
	layer_r0_meters[9]  = 1.12;
	layer_r0_meters[8]  = 0.84;
	layer_r0_meters[7]  = 0.68;
	layer_r0_meters[6]  = 0.66;
	layer_r0_meters[5]  = 0.91;
	layer_r0_meters[4]  = 0.40;
	layer_r0_meters[3]  = 0.46;
	layer_r0_meters[2]  = 0.85;
	layer_r0_meters[1]  = 1.09;
	layer_r0_meters[0]  = 1.08;
      } else {
	layer_r0_meters[0] = 0.23;
      }
    } else if(upper_layers=="bad")
      layer_r0_meters[0] = 0.16;
    */

  
    power_spectra_.resize(layer_heights_.size());
    double exponent = -11/3.0;
    double r0_ref_wavelength_meters = .5e-6;
    for(int i=0; i<layer_r0_meters.size(); i++){
      if(outer_scale_meters==-1){
	try {
	  power_spectra_[i] = 
	    new isotropic_power_law_spectrum< power_law, null_inner_scale>(power_law(exponent, 
										     layer_r0_meters[i],
										     r0_ref_wavelength_meters),
									   null_inner_scale());
	} catch(...) {
	  cerr << "Gemini_GLAO_study_model::Gemini_GLAO_study_model error - "
	       << "could not instantiate power law";
	  throw(string("Gemini_GLAO_study_model::Gemini_GLAO_study_model"));
	}
      } else {
	try {
	  power_spectra_[i] = 
	    new isotropic_power_law_spectrum<von_karman_power_law, null_inner_scale>(von_karman_power_law(exponent, 
													  layer_r0_meters[i], 
													  r0_ref_wavelength_meters, 
													  outer_scale_meters),
										     null_inner_scale());
	} catch(...) {
	  cerr << "Gemini_GLAO_study_model::Gemini_GLAO_study_model error - "
	       << "could not instantiate von Karman power law";
	  throw(string("Gemini_GLAO_study_model::Gemini_GLAO_study_model"));
	}
      }
    }

    /*
    for(int i=layer_heights_.size()-1; i>=0; i--)
      cerr << i
	   << "\t" 
	   << setw(20)
	   << layer_heights_[i] 
	   << setw(20)
	   << layer_r0_meters[i]
	   << setw(20)
	   << outer_scale_meters
	   << endl;
    */
  }

  Gemini_GLAO_study_model & 
  Gemini_GLAO_study_model::operator=(const Gemini_GLAO_study_model & cn2_model){
    if(this==&cn2_model)
      return(*this);
    this->refractive_atmospheric_model::operator=(cn2_model);
    return(*this);
  }

  void Gemini_GLAO_study_model::read(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "Gemini_GLAO_study_model::read - "
	   << "error opening file " << filename << endl;
      throw(string("Gemini_GLAO_study_model::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "Gemini_GLAO_study_model::read - "
	   << "error reading "
	   << this->unique_name() << " from file " 
	   << filename << endl;
      throw(string("Gemini_GLAO_study_model::read"));
    }
  }

  void Gemini_GLAO_study_model::read(const iofits & iof){
    if(!iof.key_exists("TYPE")){
      cerr << "Gemini_GLAO_study_model::read error - "
	   << "unrecognized type of file\n";
      throw(string("Gemini_GLAO_study_model::read"));
    }
    string type, comment;
    iof.read_key("TYPE", type, comment);
    if(type!=this->unique_name()){
      cerr << "Gemini_GLAO_study_model::read error - file of type " 
	   << type << " rather than of type "
	   << this->unique_name() << endl;
      throw(string("Gemini_GLAO_study_model::read"));
    }

    this->read_common_data(iof);

  }

  void Gemini_GLAO_study_model::write(const char * filename) const {
    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "Gemini_GLAO_study_model::write - "
	   << "error opening file " << filename << endl;
      throw(string("Gemini_GLAO_study_model::write"));
    }
    try{this->write(iof);}
    catch(...){
      cerr << "Gemini_GLAO_study_model::write - "
	   << "error writing "
	   << this->unique_name() << " to file " 
	   << filename << endl;
      throw(string("Gemini_GLAO_study_model::write"));
    }
  }

  void Gemini_GLAO_study_model::write(iofits & iof) const {

    fits_header_data<char> tmphdr;
    tmphdr.write(iof);

    string type = this->unique_name();
    string comment = "object type";
    iof.write_key("TYPE", type, comment);

    this->write_common_data(iof);
  }

  void Gemini_GLAO_study_model::
  print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << this->unique_name()
       << "/" << "object type" << endl;
    os << prefix << "NPSPEC     = " << setw(vlspc) << power_spectra_.size()
       << "/" << "number of power spectra" << endl;
    ground_ref_frame_.print(os, prefix);

    for(int i=0; i<power_spectra_.size(); i++){
      power_spectra_[i]->print(os, prefix);
      os << prefix << "HEIGHT     = " << setw(vlspc) << layer_heights_[i]
	 << "/" << "height of layer (meters)" << endl << endl;
    }
  }

}
