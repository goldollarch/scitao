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
#include "special_functions.h"

using namespace std;

namespace Arroyo {

  namespace {
    void simple_modulo(double arg, int intgr, double frac){
      while(arg>1.5){
	arg -= 1;
	intgr++;
      }
      while(arg<.5){
	intgr--;
	arg += 1;
      }
      frac = arg;
    }
  }

  double gamma_function(double arg) {

    if(fmod(arg, 1.0)==0){
      if(arg<=0){
	cerr << "gamma_function error - argument " << arg 
	     << " is an integer less than or equal to zero, where the gamma function has a pole\n";
	throw(string("gamma_function"));
      }
      double val = 1;
      for(int i=2; i<arg; i++) val *= i;
      return(val);
    }

    // Series expansion from Abramowitz and Stegun,
    // section 6.1.34
    double data[26] = {1,
		        0.5772156649015329,
		       -0.6558780715202538,		    
		       -0.0420026350340952,
		        0.1665386113822915,
		       -0.0421977345555443,
		       -0.0096219715278770,
		        0.0072189432466630,
		       -0.0011651675918591,
		       -0.0002152416741149,
		        0.0001280502823882,
		       -0.0000201348547807,
		       -0.0000012504934821,
		        0.0000011330272320,
		       -0.0000002056338417,
		        0.0000000061160950,
		        0.0000000050020075,
		       -0.0000000011812746,
		        0.0000000001043427,
		        0.0000000000077823,
		       -0.0000000000036968,
		        0.0000000000005100,
		       -0.0000000000000206,
		       -0.0000000000000054,
		        0.0000000000000014,
		        0.0000000000000001};

    double fac = 1;
    while(arg>1.5){
      arg--;
      fac*=arg;
    }
    while(arg<.5){
      fac/=arg;
      arg++;
    }

    double sum = 0;
    for(int i=25; i>=0; i--)
      sum = (sum+data[i])*arg;

    return(fac/sum);
  } 

  double gamma_ratio(const vector<double> & numerator_vals, 
		     const vector<double> & denominator_vals){

    int p = numerator_vals.size();
    int q = denominator_vals.size();
    double fac = 1;

    /*
    for(int i=0; i<p; i++)
      cerr << "\t\tnumerator " << numerator_vals[i] << endl;
    for(int i=0; i<q; i++)
      cerr << "\t\tdenominator " << denominator_vals[i] << endl;
    */

    for(int i=0; i<p && i<q; i++)
      fac *= gamma_function(numerator_vals[i])/gamma_function(denominator_vals[i]);

    if(p>q)
      for(int i=q; i<p; i++) 
	fac *= gamma_function(numerator_vals[i]);
    else 
      for(int i=p; i<q; i++) 
	fac /= gamma_function(denominator_vals[i]);
    return(fac);
  }

  double generalized_hypergeometric_function(double arg, 
					     const vector<double> & numerator_vals,
					     const vector<double> & denominator_vals, 
					     double error_tolerance) {

    if(error_tolerance<=0){
      cerr << "generalized_hypergeometric_function error - error tolerance of "
	   << error_tolerance << " is not valid\n";
      throw(string("generalized_hypergeometric_function"));
    }

    int p = numerator_vals.size();
    int q = denominator_vals.size();

    if((p>q+1) || 
       (p==q+1 && fabs(arg)>1)){
      cerr << "generalized_hypergeometric_function error - the function "
	   << p << "F" << q;
      if(p>q+1)
	cerr << " will not converge for any argument\n";
      else 
	cerr << " will not converge for argument " << arg << endl;
      throw(string("generalized_hypergeometric_function"));
    }

    double sum = 0;         // this is the answer
    double term = 1;
    long term_index = 1;    // this is k
    while(fabs(term)>.1*error_tolerance){
      sum += term;
      for(int i=0; i<numerator_vals.size(); i++)
	term *= numerator_vals[i]+term_index-1;
      for(int i=0; i<denominator_vals.size(); i++)
	term /= denominator_vals[i]+term_index-1;
      term *= arg/(double)term_index;
      term_index++;
    }
    return(sum);
  }


  double bessel_Jnu(double order, double arg){
    return(pow(.5*arg, order)*
	   generalized_hypergeometric_function(-.25*arg*arg, vector<double>(), vector<double>(1,1+order))/
	   gamma_function(1+order));
  }

  double bessel_Knu(double order, double arg){
    if(arg==0){
      cerr << "bessel_Knu error - argument supplied to this function has value zero, where this function has a pole\n";
      throw(string("bessel_Knu"));
    }

    if(fmod(order,1)==0){
      cerr <<  "bessel_Knu error - this function not yet able to handle bessel functions of integral order\n";
      throw(string("bessel_Knu"));

      // Abramowitz and Stegun eq 13.6.21
      //vector<double> pre_index(2,.5+order);
      //pre_index[1] = .5-order;
      //return(pow(2*arg, -order)*generalized_hypergeometric_function(-1/arg, pre_index, vector<double>())/sqrt(M_PI)/exp(arg));
    } else {
      // Sasiela eq. B19
      return(.5*pow(.5*arg, order)*gamma_function(-order)*
	     generalized_hypergeometric_function(.25*arg*arg, vector<double>(), vector<double>(1,1+order)) +
	     .5*pow(.5*arg, -order)*gamma_function(order)*
	     generalized_hypergeometric_function(.25*arg*arg, vector<double>(), vector<double>(1,1-order)));
    }
  }

  double get_Xi(){
    return(pow(gamma_function(1/6.),2) / gamma_function(1/3.)/5/pow(2,7/3.));
  }
}  
