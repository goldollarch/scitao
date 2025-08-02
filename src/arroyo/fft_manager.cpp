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

#include <assert.h>
#include <string.h>
#include "fft_manager.h"

using namespace std;

namespace Arroyo {

  // These template specializations need to be in a 
  // .C file rather than in the .h file with the 
  // generic template definitions.  Otherwise they
  // violate the one definition rule


  template<>
  fft_manager<float>::fft_manager(){
  }

  template<>
  fft_manager<double>::fft_manager(){
  }

  template<>
  void fft_manager<float>::forward_fft(const vector<long> & in_array_dimens,
			bool estimate, bool in_place,
  			int howmany, float * in, int istride, int idist, 
			float * out, int ostride, int odist) {

    if(verbose_level) cout << "fft_manager::forward_fft - performing fft\n";

    int rank = in_array_dimens.size();
    int * intarray_dimens = new int[rank];
    assert(intarray_dimens != NULL);
    int nelem = 1;
    for (int i = 0; i < rank; i++) {
      intarray_dimens[i] = in_array_dimens[i];
      nelem *= in_array_dimens[i];
    }
    int flags = 0;
    if (estimate) flags |= FFTW_ESTIMATE;
    else flags |= FFTW_MEASURE;
    fftwf_plan plan = NULL;
    if (in_place) {
      float * buf = new float[nelem*2];
      memcpy(buf, in, sizeof(float)*nelem*2);
      plan = fftwf_plan_many_dft(rank, intarray_dimens, howmany,
      		(fftwf_complex*)in, NULL, istride, idist,
		(fftwf_complex*)in, NULL, istride, idist,
		FFTW_FORWARD, flags);
      memcpy(in, buf, sizeof(float)*nelem*2);
      delete [] buf;
    } else {
      plan = fftwf_plan_many_dft(rank, intarray_dimens, howmany,
      		(fftwf_complex*)in, NULL, istride, idist,
		(fftwf_complex*)out, NULL, ostride, odist,
		FFTW_FORWARD, flags);
      plan = fftwf_plan_many_dft(rank, intarray_dimens, howmany,
      		(fftwf_complex*)in, NULL, istride, idist,
		(fftwf_complex*)out, NULL, ostride, odist,
		FFTW_FORWARD, flags);

    }
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
    delete [] intarray_dimens;
  } 

  template<>
  void fft_manager<double>::forward_fft(const vector<long> & in_array_dimens,
			bool estimate, bool in_place,
  			int howmany, double * in, int istride, int idist, 
			double * out, int ostride, int odist) {

    if(verbose_level) cout << "fft_manager::forward_fft - performing fft\n";

    int rank = in_array_dimens.size();
    int * intarray_dimens = new int[rank];
    assert(intarray_dimens != NULL);
    int nelem = 1;
    for (int i = 0; i < rank; i++) {
      intarray_dimens[i] = in_array_dimens[i];
      nelem *= in_array_dimens[i];
    }
    int flags = 0;
    if (estimate) flags |= FFTW_ESTIMATE;
    else flags |= FFTW_MEASURE;
    fftw_plan plan = NULL;
    if (in_place) {
      double * buf = new double[nelem*2];
      memcpy(buf, in, sizeof(double)*nelem*2);
      plan = fftw_plan_many_dft(rank, intarray_dimens, howmany,
      		(fftw_complex*)in, NULL, istride, idist,
		(fftw_complex*)in, NULL, istride, idist,
		FFTW_FORWARD, flags);
      memcpy(in, buf, sizeof(double)*nelem*2);
      delete [] buf;
    } else {
      plan = fftw_plan_many_dft(rank, intarray_dimens, howmany,
      		(fftw_complex*)in, NULL, istride, idist,
		(fftw_complex*)out, NULL, ostride, odist,
		FFTW_FORWARD, flags);
    }
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    delete [] intarray_dimens;
  } 

  template<>
  void fft_manager<float>::backward_fft(const vector<long> & in_array_dimens,
			bool estimate, bool in_place,
  			int howmany, float * in,  int istride, int idist, 
			float * out, int ostride, int odist) {

    if(verbose_level) cout << "fft_manager::forward_fft - performing fft\n";

    int rank = in_array_dimens.size();
    int * intarray_dimens = new int[rank];
    assert(intarray_dimens != NULL);
    int nelem = 1;
    for (int i = 0; i < rank; i++) {
      intarray_dimens[i] = in_array_dimens[i];
      nelem *= in_array_dimens[i];
    }
    int flags = 0;
    if (estimate) flags |= FFTW_ESTIMATE;
    else flags |= FFTW_MEASURE;
    fftwf_plan plan = NULL;
    if (in_place) {
      float * buf = new float[nelem*2];
      memcpy(buf, in, sizeof(float)*nelem*2);
      plan = fftwf_plan_many_dft(rank, intarray_dimens, howmany,
      		(fftwf_complex*)in, NULL, istride, idist,
		(fftwf_complex*)in, NULL, istride, idist,
		FFTW_FORWARD, flags);
      memcpy(in, buf, sizeof(float)*nelem*2);
      delete [] buf;
    } else {
      plan = fftwf_plan_many_dft(rank, intarray_dimens, howmany,
      		(fftwf_complex*)in, NULL, istride, idist,
		(fftwf_complex*)out, NULL, ostride, odist,
		FFTW_FORWARD, flags);
    }
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
    delete [] intarray_dimens;
  }

  template<>
  void fft_manager<double>::backward_fft(const vector<long> & in_array_dimens,
			bool estimate, bool in_place,
  			int howmany, double * in, int istride, int idist, 
			double * out, int ostride, int odist) {

    if(verbose_level) cout << "fft_manager::backward_fft - performing fft\n";

    int rank = in_array_dimens.size();
    int * intarray_dimens = new int[rank];
    assert(intarray_dimens != NULL);
    int nelem = 1;
    for (int i = 0; i < rank; i++) {
      intarray_dimens[i] = in_array_dimens[i];
      nelem *= in_array_dimens[i];
    }
    int flags = 0;
    if (estimate) flags |= FFTW_ESTIMATE;
    else flags |= FFTW_MEASURE;
    fftw_plan plan = NULL;
    if (in_place) {
      double * buf = new double[nelem*2];
      memcpy(buf, in, sizeof(double)*nelem*2);
      plan = fftw_plan_many_dft(rank, intarray_dimens, howmany,
      		(fftw_complex*)in, NULL, istride, idist,
		(fftw_complex*)in, NULL, istride, idist,
		FFTW_BACKWARD, flags);
      memcpy(in, buf, sizeof(double)*nelem*2);
      delete [] buf;
    } else {
      plan = fftw_plan_many_dft(rank, intarray_dimens, howmany,
      		(fftw_complex*)in, NULL, istride, idist,
		(fftw_complex*)out, NULL, ostride, odist,
		FFTW_BACKWARD, flags);
    }
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    delete [] intarray_dimens;
  }

  // specialize the float and double instantiations,
  // so as to disallow other possibile instantiations
  template<>
  rfft_manager<float>::rfft_manager(){
  }

  template<>
  rfft_manager<double>::rfft_manager(){
  }

  template<>
  void rfft_manager<float>::real_to_complex_fft(const vector<long> & in_array_dimens,
			bool estimate, bool in_place,
  			int howmany, float * in, int istride, int idist,
			float * out, int ostride, int odist) {

    int rank = in_array_dimens.size();
    int * intarray_dimens = new int[rank];
    assert(intarray_dimens != NULL);
    int nelem = 1;
    for (int i = 0; i < rank; i++) {
      intarray_dimens[i] = in_array_dimens[i];
      nelem *= in_array_dimens[i];
    }
    int flags = 0;
    if (estimate) flags |= FFTW_ESTIMATE;
    else flags |= FFTW_MEASURE;
    fftwf_plan plan = NULL;
    if (in_place) {
      float * buf = new float[nelem];
      memcpy(buf, in, sizeof(float)*nelem);
      plan = fftwf_plan_many_dft_r2c(rank, intarray_dimens, howmany,
      		in, NULL, istride, idist,
		(fftwf_complex*)in, NULL, ostride, odist,
		flags);
      memcpy(in, buf, sizeof(float)*nelem);
      delete [] buf;
    } else {
      plan = fftwf_plan_many_dft_r2c(rank, intarray_dimens, howmany,
      		in, NULL, istride, idist,
		(fftwf_complex*)out, NULL, ostride, odist,
		flags);
    }
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
    delete [] intarray_dimens;
  }

  template<>
  void rfft_manager<double>::real_to_complex_fft(const vector<long> & in_array_dimens,
			bool estimate, bool in_place,
  			int howmany, double * in, int istride, int idist, 
			double * out, int ostride, int odist) {

    int rank = in_array_dimens.size();
    int * intarray_dimens = new int[rank];
    assert(intarray_dimens != NULL);
    int nelem = 1;
    for (int i = 0; i < rank; i++) {
      intarray_dimens[i] = in_array_dimens[i];
      nelem *= in_array_dimens[i];
    }
    int flags = 0;
    if (estimate) flags |= FFTW_ESTIMATE;
    else flags |= FFTW_MEASURE;
    fftw_plan plan = NULL;
    if (in_place) {
      double * buf = new double[nelem];
      memcpy(buf, in, sizeof(double)*nelem);
      plan = fftw_plan_many_dft_r2c(rank, intarray_dimens, howmany,
      		in, NULL, istride, idist,
		(fftw_complex*)in, NULL, ostride, odist,
		flags);
      memcpy(in, buf, sizeof(double)*nelem);
      delete [] buf;
    } else {
      plan = fftw_plan_many_dft_r2c(rank, intarray_dimens, howmany,
      		in, NULL, istride, idist,
		(fftw_complex*)out, NULL, ostride, odist,
		flags);
    }
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    delete [] intarray_dimens;
  }

  template<>
  void rfft_manager<float>::complex_to_real_fft(const vector<long> & in_array_dimens,
			bool estimate, bool in_place,
  			int howmany,  float * in, int istride, int idist, 
			float * out, int ostride, int odist) {

    int rank = in_array_dimens.size();
    int * intarray_dimens = new int[rank];
    assert(intarray_dimens != NULL);
    int nelem = 1;
    for (int i = 0; i < rank; i++) {
      intarray_dimens[i] = in_array_dimens[i];
      nelem *= in_array_dimens[i];
    }
    int flags = 0;
    if (estimate) flags |= FFTW_ESTIMATE;
    else flags |= FFTW_MEASURE;
    fftwf_plan plan = NULL;
    if (in_place) {
      float * buf = new float[nelem];
      memcpy(buf, in, sizeof(float)*nelem);
      plan = fftwf_plan_many_dft_c2r(rank, intarray_dimens, howmany,
      		(fftwf_complex*)in, NULL, istride, idist,
		in, NULL, ostride, odist,
		flags);
      memcpy(in, buf, sizeof(float)*nelem);
      delete [] buf;
    } else {
      plan = fftwf_plan_many_dft_c2r(rank, intarray_dimens, howmany,
      		(fftwf_complex*)in, NULL, istride, idist,
		out, NULL, ostride, odist,
		flags);
    }
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
    delete [] intarray_dimens;
  }

  template<>
  void rfft_manager<double>::complex_to_real_fft(const vector<long> & in_array_dimens,
			bool estimate, bool in_place,
  			int howmany, double * in,  int istride, int idist, 
			double * out, int ostride, int odist) {

    int rank = in_array_dimens.size();
    int * intarray_dimens = new int[rank];
    assert(intarray_dimens != NULL);
    int nelem = 1;
    for (int i = 0; i < rank; i++) {
      intarray_dimens[i] = in_array_dimens[i];
      nelem *= in_array_dimens[i];
    }
    int flags = 0;
    if (estimate) flags |= FFTW_ESTIMATE;
    else flags |= FFTW_MEASURE;
    fftw_plan plan = NULL;
    if (in_place) {
      double * buf = new double[nelem];
      memcpy(buf, in, sizeof(double)*nelem);
      plan = fftw_plan_many_dft_c2r(rank, intarray_dimens, howmany,
      		(fftw_complex*)in, NULL, istride, idist,
		in, NULL, ostride, odist,
		flags);
      memcpy(in, buf, sizeof(double)*nelem);
      delete [] buf;
    } else {
      plan = fftw_plan_many_dft_c2r(rank, intarray_dimens, howmany,
      		(fftw_complex*)in, NULL, istride, idist,
		out, NULL, ostride, odist,
		flags);
    }
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    delete [] intarray_dimens;
  }
}
