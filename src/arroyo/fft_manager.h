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

#ifndef FFT_MANAGER_H
#define FFT_MANAGER_H

#ifndef FFTW_REAL_TO_COMPLEX
#define FFTW_REAL_TO_COMPLEX FFTW_FORWARD
#define FFTW_COMPLEX_TO_REAL FFTW_BACKWARD
#endif

#include <fftw3.h>
#include <iostream>
#include <vector>
#include <string>
#include "AO_cpp.h"


namespace Arroyo {

  using std::string;
  using std::vector;
  using std::ostream;

  /// A template class to facilitate complex to complex forward and backward fft's.
  /// Currently this class supports fftw, but support for additional fft libraries
  /// may be added at a later time.  For example, conditional compilation of native
  /// fft libraries should be possible.
  ///
  /// This class aims to make arbitrary-dimensional fast fourier transforms
  /// relatively straightforward, and supports all the FFTW options 

  template<class T>
    class fft_manager {

    private:
  

    protected:


    public:

    ///////////////////////////////////////////
    ///  Null constructor
    fft_manager();

    ///////////////////////////////////////////
    ///  Copy constructor
    fft_manager(const fft_manager<T> & fft_mgr);

    ///////////////////////////////////////////
    ///  Destructor
    ~fft_manager();

    ///////////////////////////////////////////
    ///  Operator =
    fft_manager<T> & operator=(const fft_manager<T> & fft_mgr);

    ///////////////////////////////////////////
    ///  Perform a complex to complex forward transform
    ///  input data should be sorted real, imag, real, imag...
    ///  The "out" argument is optional - if you include it an
    ///  out of place transform will be performed and the results
    ///  will be stored there.  
    void forward_fft(const vector<long> & in_array_dimens,
			bool estimate, bool in_place,
    			T * in, T * out = NULL);
  
    ///////////////////////////////////////////
    ///  Perform a complex to complex forward transform
    ///  input data should be sorted real, imag, real, imag...
    ///  The "out" argument is optional - if you include it an
    ///  out of place transform will be performed and the results
    ///  will be stored there.  
    ///  See the FFTW documentation for a description of istride, idist,
    ///  ostride and odist
    void forward_fft(const vector<long> & in_array_dimens,
			bool estimate, bool in_place,
    			int howmany, T * in, int istride, int idist, 
			T * out=NULL, int ostride=1, int odist=1);
  
    ///////////////////////////////////////////
    ///  Perform a complex to complex backward transform
    ///  input data should be sorted real, imag, real, imag...
    ///  The "out" argument is optional - if you include it an
    ///  out of place transform will be performed and the results
    ///  will be stored there.  
    void backward_fft(const vector<long> & in_array_dimens,
			bool estimate, bool in_place,
    			T * in, T * out = NULL);

    ///////////////////////////////////////////
    ///  Perform a complex to complex forward transform
    ///  input data should be sorted real, imag, real, imag...
    ///  The "out" argument is optional - if you include it an
    ///  out of place transform will be performed and the results
    ///  will be stored there.  
    ///  See the FFTW documentation for a description of istride, idist,
    ///  ostride and odist
    void backward_fft(const vector<long> & in_array_dimens,
			bool estimate, bool in_place,
    			int howmany, T * in, int istride, int idist, 
			T * out=NULL, int ostride=1, int odist=1);

    /// A verbose_level for printing messages
    static int verbose_level;

  };

  template<class T>
    int fft_manager<T>::verbose_level = 0;

  template<class T>
    fft_manager<T>::fft_manager(const fft_manager<T> & fft_mgr){

    this->operator=(fft_mgr);
  }

  template<class T>
    fft_manager<T>::~fft_manager(){
  }

  template<class T>
    fft_manager<T> & fft_manager<T>::operator=(const fft_manager<T> & fft_mgr){
    if(this==&fft_mgr) 
      return(*this);

    return(*this);
  }

  template<class T>
    void fft_manager<T>::forward_fft(const vector<long> & in_array_dimens,
			bool estimate, bool in_place,
    			T * in, T * out) {
      this->forward_fft(in_array_dimens, estimate, in_place, 1, in, 1, 1, out, 1, 1);
  }

  template<class T>
    void fft_manager<T>::backward_fft(const vector<long> & in_array_dimens,
			bool estimate, bool in_place,
    			T * in, T * out) {
    this->backward_fft(in_array_dimens, estimate, in_place, 1, in, 1, 1, out, 1, 1);
  }


  /// A template class to facilitate real to complex forward and backward fft's
  /// using fftw.  This class aims to make arbitrary-dimensional fast
  /// fourier transforms relatively straightforward, and supports all
  /// the FFTW options 

  template<class T>
    class rfft_manager {

    private:
  
 
    protected:


    public:

    ///////////////////////////////////////////
    ///  Null constructor
    rfft_manager();

    ///////////////////////////////////////////
    ///  Copy constructor
    rfft_manager(const rfft_manager<T> & rfft_mgr);

    ///////////////////////////////////////////
    ///  Destructor
    ~rfft_manager();

    ///////////////////////////////////////////
    ///  Operator =
    rfft_manager<T> & operator=(const rfft_manager<T> & rfft_mgr);

    ///////////////////////////////////////////
    ///  Perform a real to complex transform
    ///  input data should be sorted real, real, real...
    ///  Some warnings are in order here:
    ///    The transformed array requires slightly
    ///    more storage than the original if the final 
    ///    axis is even, because of nyquist.
    ///    See the FFTW home page for a discussion of this
    ///    issue and for the output array format.
    void real_to_complex_fft(const vector<long> & in_array_dimens,
			bool estimate, bool in_place,
    			T * in, T * out = NULL) ;

    ///////////////////////////////////////////
    ///  Perform a real to complex transform
    ///  as above, with more options
    ///  The "out" argument is optional - if you include it an
    ///  out of place transform will be performed and the results
    ///  will be stored there.  
    ///  See the FFTW documentation for a description of istride, idist,
    ///  ostride and odist
    void real_to_complex_fft(const vector<long> & in_array_dimens,
			bool estimate, bool in_place,
    			int howmany, T * in, int istride, int idist, 
			T * out=NULL, int ostride=1, int odist=1) ;
 
    ///////////////////////////////////////////
    ///  Perform a complex to real transform
    ///  input data should be sorted real, imag, real, imag...
    ///  The "out" argument is optional - if you include it an
    ///  out of place transform will be performed and the results
    ///  will be stored there.  
    void complex_to_real_fft(const vector<long> & in_array_dimens,
			bool estimate, bool in_place,
    			T * in, T * out = NULL) ;

    ///////////////////////////////////////////
    ///  Perform a complex to real transform
    ///  as above, with more options
    ///  The "out" argument is optional - if you include it an
    ///  out of place transform will be performed and the results
    ///  will be stored there.  
    ///  See the FFTW documentation for a description of istride, idist,
    ///  ostride and odist
    void complex_to_real_fft(const vector<long> & in_array_dimens,
			bool estimate, bool in_place,
    			int howmany, T * in, int istride, int idist, 
			T * out=NULL, int ostride=1, int odist=1) ;
 
    /// A verbose_level for printing messages
    static int verbose_level;

  };

  template<class T>
    int rfft_manager<T>::verbose_level = 0;


  template<class T>
    rfft_manager<T>::rfft_manager(const rfft_manager<T> & rfft_mgr){

    this->operator=(rfft_mgr);
  }

  template<class T>
    rfft_manager<T>::~rfft_manager(){

  }

  template<class T>
    rfft_manager<T> & rfft_manager<T>::operator=(const rfft_manager<T> & rfft_mgr){
    if(this==&rfft_mgr) 
      return(*this);

    return(*this);
  }

  template<class T>
    void rfft_manager<T>::real_to_complex_fft(
    		const vector<long> & in_array_dimens,
		bool estimate, bool in_place,
		T * in, T * out) {
      this->real_to_complex_fft(in_array_dimens, estimate, in_place,
    				1, in, 1, 1, out, 1, 1);
  }

  template<class T>
    void rfft_manager<T>::complex_to_real_fft(
    		const vector<long> & in_array_dimens,
		bool estimate, bool in_place,
		T * in, T * out) {
      this->complex_to_real_fft(in_array_dimens, estimate, in_place,
      				1, in, 1, 1, out, 1, 1);
  }

}

#endif
