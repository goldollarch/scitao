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

//#include <g2c.h>
#include <iostream>
#include "linear_algebra.h"

extern "C" {void dgesvd_(char*, char*, int*, int*, double *, int*, double *, double *, int*, double *, int*, double *, int*, int*);}
extern "C" {void sgesvd_(char*, char*, int*, int*, float *, int*, float *, float *, int*, float *, int*, float *, int*, int*);}


using namespace std;

namespace Arroyo {

  namespace {
    
    void check_svd_args(char * jobu, char * jobvt, int m, int n, int lda, int ldu, int ldvt, int lwork, int info){
      
      if(jobu[0]!='A' && jobu[0]!='S' && jobu[0]!='O' && jobu[0]!='N'){
	cerr << "check_svd_args error - invalid jobu code " << jobu << endl;
	throw(string("check_svd_args"));
      }
      if(jobvt[0]!='A' && jobvt[0]!='S' && jobvt[0]!='O' && jobvt[0]!='N'){
	cerr << "check_svd_args error - invalid jobvt code " << jobvt << endl;
	throw(string("check_svd_args"));
      }
      if(m<0){
	cerr << "check_svd_args error - invalid M dimension " << m << endl;
	throw(string("check_svd_args"));
      }
      if(n<0){
	cerr << "check_svd_args error - invalid N dimension " << n << endl;
	throw(string("check_svd_args"));
      }
    }
  }
  
  template<>
    void singular_value_decomposition<float>(char * jobu,
					     char * jobvt,
					     int & m,
					     int & n,
					     float * a,
					     int & lda,
					     float * s,
					     float * u,
					     int & ldu,
					     float * vt,
					     int & ldvt,
					     float * work,
					     int & lwork,
					     int & info) {

    try{
      check_svd_args(jobu, jobvt, m, n, lda, ldu, ldvt, lwork, info);
    } catch(...){
      cerr << "singular_value_decomposition error - invalid arguments passed to this function\n";
      throw(string("singular_value_decomposition"));
    }

    sgesvd_(jobu, jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, &info);
    if(info<0){
      cerr << "singular_value_decomposition<double> - error calling LAPACK dgesvd\n";
      throw(string("singular_value_decomposition"));
    }
  }

  template<>
  void singular_value_decomposition<double>(char * jobu,
					    char * jobvt,
					    int & m,
					    int & n,
					    double * a,
					    int & lda,
					    double * s,
					    double * u,
					    int & ldu,
					    double * vt,
					    int & ldvt,
					    double * work,
					    int & lwork,
					    int & info) {

    try{
      check_svd_args(jobu, jobvt, m, n, lda, ldu, ldvt, lwork, info);
    } catch(...){
      cerr << "singular_value_decomposition error - invalid arguments passed to this function\n";
      throw(string("singular_value_decomposition"));
    }
    
    dgesvd_(jobu, jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, &info);
    if(info<0){
      cerr << "singular_value_decomposition<double> - error calling LAPACK dgesvd\n";
      throw(string("singular_value_decomposition"));
    }
  }
}

