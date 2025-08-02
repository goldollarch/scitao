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
#include <iomanip>
#include <algorithm>
#include "AO_algo.h"

using namespace std;

namespace Arroyo {

  double chisq_2gf(const vector<float> & data, vector<double> & centers, 
		   vector<double> & widths, vector<double> & amps, 
		   double & dcconst, int verbose){
    double tmp;
    double chisq = 0;
    for(int i=0; i<data.size(); i++){
      tmp = data[i] -
	amps[0] * exp(-(i-centers[0])*(i-centers[0])/(2*widths[0]*widths[0])) -
	amps[1] * exp(-(i-centers[1])*(i-centers[1])/(2*widths[1]*widths[1])) - dcconst;
      if(verbose==2)
	cout << data[i] << "\t" 
	     << amps[0] * exp(-(i-centers[0])*(i-centers[0])/(2*widths[0]*widths[0])) +
	  amps[1] * exp(-(i-centers[1])*(i-centers[1])/(2*widths[1]*widths[1])) + dcconst << endl;
      chisq += tmp*tmp;
    }
    if(verbose==2)
      cout << chisq << "\t" << amps[0] << "\t" << amps[1] << "\t" << centers[0] << "\t"
	   <<  centers[1] << "\t" << widths[0] << "\t" << widths[1] << "\t" << dcconst << endl;
    return(chisq);
  }

  void twogauss_fit(const vector<float> & data, vector<double> & centers, 
		    vector<double> & widths, vector<double> & amps, double & dcconst){

    double chisq, schisq;
    vector<double> deriv_amps(2), deriv_widths(2), deriv_centers(2);
    double deriv_dcconst;

    double delta_amp = 10.0, delta_center = .1, delta_width = .1, delta_dcconst = 1;
    double grad, frac = .1;
    int count = 0;
    while(1){
      count ++;
      if(count%20==0)
	chisq = chisq_2gf(data, centers, widths, amps, dcconst, 0);
      else
	chisq = chisq_2gf(data, centers, widths, amps, dcconst, 0);
      for(int i=0; i<2; i++){
	amps[i] += delta_amp*frac;
	schisq = chisq_2gf(data, centers, widths, amps, dcconst, 0);
	deriv_amps[i] = schisq - chisq;
	amps[i] -= delta_amp*frac;      

	widths[i] += delta_width*frac;
	schisq = chisq_2gf(data, centers, widths, amps, dcconst, 0);
	deriv_widths[i] = schisq - chisq;
	widths[i] -= delta_width*frac;

	centers[i] += delta_center*frac;
	schisq = chisq_2gf(data, centers, widths, amps, dcconst, 0);
	deriv_centers[i] = schisq - chisq;
	centers[i] -= delta_center*frac;

	if(i==0){
	  dcconst += delta_dcconst*frac;
	  schisq = chisq_2gf(data, centers, widths, amps, dcconst, 0);
	  deriv_dcconst = schisq - chisq;
	  dcconst -= delta_dcconst*frac;
	}	
      }
      grad = 0;
      for(int i=0; i<2; i++){
	grad += deriv_amps[i]*deriv_amps[i]*delta_amp*delta_amp;
	grad += deriv_widths[i]*deriv_widths[i]*delta_width*delta_width;
	grad += deriv_centers[i]*deriv_centers[i]*delta_center*delta_center;
	if(i==0) grad += deriv_dcconst*deriv_dcconst*delta_dcconst*delta_dcconst;
      }
      grad = sqrt(grad);

      for(int i=0; i<2; i++){
	amps[i] -= deriv_amps[i] * delta_amp / grad;
	widths[i] -= deriv_widths[i] * delta_width / grad;
	centers[i] -= deriv_centers[i] * delta_center / grad;
	if(i==0)
	  dcconst -= deriv_dcconst * delta_dcconst / grad;
      }    
      schisq = chisq_2gf(data, centers, widths, amps, dcconst, 0);
      if(chisq < schisq || fabs(schisq - chisq) < 1e-2){
	//      cout << "twogaussfit - chisq " << chisq << "\t";
	//      chisq_2gf(data, centers, widths, amps, dcconst, 2);
	return;
      }
    }
  }

}
