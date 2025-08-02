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

#ifndef AO_ALGO_H
#define AO_ALGO_H

#include <cmath>
#include <vector>
#include <string> 
#include "sim_utils.h"
#include "AO_cpp.h"
#include "three_frame.h"


namespace Arroyo {

  ///////////////////////////////////////////////
  /// function used in twogauss_fit
  double chisq_2gf(const vector<float> & data, vector<double> & centers, 
		   vector<double> & widths, vector<double> & amps, 
		   double & dcconst, int verbose = 0);

  ///////////////////////////////////////////////
  /// fit gaussians
  void twogauss_fit(const vector<float> & data, vector<double> & centers, 
	vector<double> & widths, vector<double> & amps, double & dcconst);

  ///////////////////////////////////////////////
  /// return a random variable with 2d gaussian distribution
  inline void box_mueller(double & r1, double & r2) {
    double r, fac;

    do {
#if defined(__sun)
      r1 = 2*(random()/(double)LONG_MAX) - 1;
      r2 = 2*(random()/(double)LONG_MAX) - 1;
#else
      r1 = 2*(random()/(double)RAND_MAX) - 1;
      r2 = 2*(random()/(double)RAND_MAX) - 1;
#endif
      r = r1*r1+r2*r2;
    } while(r>=1 || r==0);

    fac = sqrt(-2*log(r)/r);
    r1 *= fac;
    r2 *= fac;
  }

  ///////////////////////////////////////////////
  /// Modification of numerical recipes algorithm
  /// for fitting a straight line
  template<class U, class V>
    void fit_straight_line(const long & nelems, const double * x,
			const U * y, const V * wts, 
			double &slope, double & sigsq_slope,
			double & yintercept, double & sigsq_yintercept, 
			double & chisq){

    slope = sigsq_slope = yintercept = sigsq_yintercept = chisq = 0;

    if(nelems <2){
      cerr << "fit_straight_line error - only " << nelems << " points\n";
      throw(string("fit_straight_line"));
    }

    double s, sx, sy, s_tt, ti;

    int nnonzerowts = 0;
    for(int i=0; i<nelems; i++)
      if(wts[i] != 0) nnonzerowts++;

    // If two or more weights are nonzero, proceed with weights
    if(nnonzerowts>=3){
      s = sx = sy = s_tt = 0;
      for(int i=0; i<nelems; i++){
	s += wts[i]*wts[i];
	sx += x[i]*wts[i]*wts[i];
	sy += y[i]*wts[i]*wts[i];
      }
      for(int i=0; i<nelems; i++){
	ti = wts[i]*(x[i] - sx/s);
	s_tt += ti*ti;
	slope += wts[i]*ti*y[i];
      }
      /*      cout << endl << "weights: s " << s << "\t sx " << sx << "\tsy " << sy */
      /*  	 << "\ts_tt " << s_tt; */
    } 
    // otherwise proceed w/o weights 
    else if(nnonzerowts<=2){  
      sx = sy = s_tt = 0;
      s = (double)nelems;
      for(int i=0; i<nelems; i++){
	sx += x[i];
	sy += y[i];
      }
      for(int i=0; i<nelems; i++){
	ti = x[i] - sx/s;
	s_tt += ti*ti;
	slope += ti*y[i];
      }
      /*      cout << "no weights: s " << s << "\t sx " << sx << "\tsy " << sy */
      /*  	 << "\ts_tt " << s_tt; */
    }

    slope /= s_tt;
    yintercept = (sy - sx*slope)/s;


    double chi = 0;
    if(nnonzerowts>=3){
      sigsq_yintercept = (1+sx*sx/(s*s_tt))/s;
      sigsq_slope = 1/s_tt; 
      for(int i=0; i<nelems; i++){
	chi = (y[i] - yintercept - slope*x[i])*wts[i];
	chisq += chi*chi;
      }
    } else {
      for(int i=0; i<nelems; i++){
	chi = (y[i] - yintercept - slope*x[i]);
	chisq += chi*chi;
      }
    } 
  }


  ///////////////////////////////////////////////
  /// This isn't a particularly safe function...
  /// It cyclically permutes data by xshift and yshift.
  /// The zeroth element of the axes vector corresponds
  /// to the y direction, and the first element of the 
  /// axes vector corresponds to the x direction
  template<class T>
    void cyclic_permutation(vector<long> axes, long xshift,
    					long yshift, T * data){

    while(xshift<0) xshift += axes[1];
    while(yshift<0) yshift += axes[0];

    if(xshift==0 && yshift==0) return;

    long tmpelem = axes[0];
    if(axes[1]>axes[0]) tmpelem=axes[1];

    T * tmparr;
    try{tmparr = new T[tmpelem];}
    catch(...){
      cerr << "cyclic_permutation error - could not allocate memory: "
      			<< tmpelem << " bytes\n";
      throw;
    }

    if(yshift!=0){
      for(int i=0; i<axes[1]; i++){
	for(int j=0; j<axes[0]; j++)
	  tmparr[j] = data[i*axes[0]+j];
	for(int j=0; j<axes[0]; j++)
	  data[i*axes[0]+j] = tmparr[(j+yshift)%axes[0]];
      }
    }
    if(xshift!=0){
      for(int j=0; j<axes[0]; j++){
	for(int i=0; i<axes[1]; i++)
	  tmparr[i] = data[i*axes[0]+j];
	for(int i=0; i<axes[1]; i++)
	  data[i*axes[0]+j] = tmparr[(i+xshift)%axes[1]];
      }
    }
    delete [] tmparr;
  }


  ///////////////////////////////////////////////
  /// This isn't a particularly safe function...
  /// It cyclically permutes data by xshift and yshift.
  /// The zeroth element of the axes vector corresponds
  /// to the y direction, and the first element of the 
  /// axes vector corresponds to the x direction
  template<class T>
    void complex_cyclic_permutation(vector<long> axes,
    			long xshift, long yshift, T * data){

    while(xshift<0) xshift += axes[1];
    while(yshift<0) yshift += axes[0];

    if(xshift==0 && yshift==0) return;

    long tmpelem = 2*axes[0];
    if(axes[1]>axes[0]) tmpelem=2*axes[1];

    T * tmparr;
    //try{tmparr = new T[2*tmpelem];}
    //alloc_size sz(tmpelem, 2);
    try{
      //tmparr = new T[sz];
      tmparr = new T[2*tmpelem];
    }
    catch(...){
      cerr << "complex_cyclic_permutation error - could not allocate memory: "
	   << 2*tmpelem*sizeof(T) << " bytes\n";
      throw;
    }

    if(yshift!=0){
      for(int i=0; i<axes[1]; i++){
	for(int j=0; j<axes[0]; j++){
	  tmparr[2*j] = data[2*(i*axes[0]+j)];
	  tmparr[2*j+1] = data[2*(i*axes[0]+j)+1];
	}
	for(int j=0; j<axes[0]; j++){
	  data[2*(i*axes[0]+j)] = tmparr[2*((j+yshift)%axes[0])];
	  data[2*(i*axes[0]+j)+1] = tmparr[2*((j+yshift)%axes[0])+1];
	}
      }
    }

    if(xshift!=0){
      for(int j=0; j<axes[0]; j++){
	for(int i=0; i<axes[1]; i++){
	  tmparr[2*i] = data[2*(i*axes[0]+j)];
	  tmparr[2*i+1] = data[2*(i*axes[0]+j)+1];
	}
	for(int i=0; i<axes[1]; i++){
	  data[2*(i*axes[0]+j)] = tmparr[2*((i+xshift)%axes[1])];
	  data[2*(i*axes[0]+j)+1] = tmparr[2*((i+xshift)%axes[1])+1];
	}
      }
    }

    delete [] tmparr;
  }


}

#endif
