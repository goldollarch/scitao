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

#ifndef STRUCTURE_FUNCTION_H
#define STRUCTURE_FUNCTION_H

#include <ostream>
#include "AO_sim_base.h"
#include "iofits.h"
#include "fits_header_data.h"
#include "pixel_array.h"
#include "power_spectrum.h"

namespace Arroyo {

  using std::string;
  using std::ostream;
  using std::vector;

  // class forward declaration
  template<class T>
  class refractive_atmospheric_layer;

  /// A class to represent structure functions.  The structure function can
  /// be one or two dimensional, or zero dimensional when created by the
  /// null constructor.
  ///
  /// Also, the add_statistics member function can leave the structure
  /// function in a state in which weights for the underlying pixel_array
  /// have been allocated and the data itself is unnormalized by these
  /// weights.  Alternatively, the data may have been normalized by a call
  /// to pixel_array<double>::normalize_by_weights.  Finally, no weights
  /// may exist.  Thus there are 3 possible states for an instance of this
  /// class.  However, none of the member functions need to check the
  /// state - they all initialize the pixel array to whatever state they
  /// require.  In the future, I forsee only the plotting member function
  /// being affected by the state.  This function will normalize by
  /// weights if the weights exist, and then will plot the data.

  class structure_function :
    virtual public AO_sim_base, 
    public pixel_array<double> {
    
    private:

    static const bool factory_registration;
  
    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("structure function"));};

    protected:

    /// The pixel scale in meters
    double pixel_scale;

    public:
  
    ///////////////////////////////////////////
    ///  Null constructor
    structure_function();
 
    ///////////////////////////////////////////
    ///  Copy constructor
    structure_function(const structure_function & structfn);
  
    ///////////////////////////////////////////
    ///  Construct from file
    structure_function(const char * filename);
  
    ///////////////////////////////////////////
    ///  Construct from an iofits object
    structure_function(const iofits & iof);
  
    ///////////////////////////////////////////
    ///  Construct from the pieces
    structure_function(const pixel_array<double> & pixarr, double pixscale);
  
    ///////////////////////////////////////////
    ///  Construct the theoretical structure function from a power
    ///  spectrum.
    ///  
    ///  This constructor throws an error unless the power spectrum
    ///  is isotropic with a von karman outer scale and no inner scale.
    structure_function(const power_spectrum * pspectrum, 
		       vector<long> axes, 
		       double pixscale);
  
    ///////////////////////////////////////////
    ///  Construct the expected structure function from a power
    ///  spectrum.
    structure_function(const power_spectrum * pspectrum, 
		       const subharmonic_method & subm,
		       vector<long> axes, 
		       double pixscale);
  
    ///////////////////////////////////////////
    ///  Destructor
    ~structure_function(){};
  
    ///////////////////////////////////////////
    ///  Operator = 
    structure_function &
    operator=(const structure_function & structfn);
  
    ///////////////////////////////////////////
    ///  Set the axes - this destroys and
    ///  reallocates the underlying pixel array
    ///  if in_axes != this->axes
    void set_axes(const vector<long> & in_axes);

    ///////////////////////////////////////////
    ///  Return the pixel scale
    double get_pixel_scale() const {return(pixel_scale);};

    ///////////////////////////////////////////
    ///  Return the axes
    vector<long> get_axes() const {return(pixel_array<double>::axes);};

    ///////////////////////////////////////////
    ///  Read from file
    void read(const char * filename);
  
    ///////////////////////////////////////////
    ///  Read from iofits
    void read(const iofits & iof);
 
    ///////////////////////////////////////////
    ///  Write to file
    void write(const char * filename) const;
 
    ///////////////////////////////////////////
    ///  Write to iofits
    void write(iofits & iof) const;
 
    ///////////////////////////////////////////
    ///  Print
    void print(ostream & os, const char * prefix="") const;

    ///////////////////////////////////////////
    ///  Incorporate statistics of a
    ///  refractive atmospheric layer into
    ///  the structure function.  
    ///  Within this function, the axes of this 
    ///  structure function are reset to a size
    ///  corresponding to the maximum allowed 
    ///  by the dimensions of the 
    ///  refractive_atmospheric_layer
    ///  The pixel_scale of the structure function
    ///  is also set to that of the 
    ///  refractive_atmospheric_layer
    template<class T>
    void add_statistics(const refractive_atmospheric_layer<T> & ref_atm_layer){

      vector<long> ref_layer_axes = 
        ref_atm_layer.get_axes();    
      double ref_layer_pixel_scale = ref_atm_layer.get_pixel_scale();

      if(this->get_axes().size()==1){

        // If this is the first time, and weights haven't been allocated,
        // we reinitialize the pixel array and set the weights
        // Also, if the weights have been allocated but the ref_atm_layer
        // is of a different size, we reinitialize as above
        vector<long> sorted_axes = ref_layer_axes;
        sort(sorted_axes.begin(), sorted_axes.end());

        // if the long axis is even, subtract one element.
        long truncated_axis = sorted_axes[0];
        if(sorted_axes[0]%2==0) truncated_axis--;
        long half_truncated_axis = truncated_axis/2;

        if(!this->weights_allocated() || 
	   this->get_axes()[0]!=half_truncated_axis+1 ||
	   ref_layer_pixel_scale!=this->pixel_scale){
	  this->set_axes(vector<long>(1,half_truncated_axis+1));
	  this->allocate_weights(0);
	  this->pixel_scale = ref_layer_pixel_scale;
        }

        double tmp, centerval;
        if(ref_layer_axes[0]>=ref_layer_axes[1]){
	  for(int i=0; i<ref_layer_axes[1]; i++){
	    centerval = ref_atm_layer.data(i*ref_layer_axes[0]+half_truncated_axis);
	    for(int j=0; j<truncated_axis; j++){
	      tmp = ref_atm_layer.data(i*ref_layer_axes[0] + j) - centerval;
	      this->pixeldata[abs(j-half_truncated_axis)] += tmp*tmp;
	      this->pixelwts[abs(j-half_truncated_axis)]++;
	    }
	  }
        } else {
	  for(int j=0; j<ref_layer_axes[0]; j++){
	    centerval = ref_atm_layer.data((half_truncated_axis)*ref_layer_axes[0]+j);
	    for(int i=0; i<truncated_axis; i++){	
	      tmp = ref_atm_layer.data(i*ref_layer_axes[0] + j) - centerval;
	      this->pixeldata[abs(i-half_truncated_axis)] += tmp*tmp;
	      this->pixelwts[abs(i-half_truncated_axis)]++;
	    }
	  }
        }

        /*
	  double tmp;
	  for(int i=0; i<ref_layer_axes[1] && i<sorted_axes[0]; i++){

	  for(int j=i; j<ref_layer_axes[0]; j++){
	  tmp = 
	    ref_atm_layer.pixel_array<double>::data(i*ref_layer_axes[0] + j) -
	    ref_atm_layer.pixel_array<double>::data(i*ref_layer_axes[0] + i);
	    this->pixel_array<double>::pixeldata[j-i] += tmp*tmp;
	    this->pixel_array<double>::pixelwts[j-i]++;
	  }

	  for(int j=i; j<ref_layer_axes[1]; j++){
	    tmp = 
	      ref_atm_layer.pixel_array<double>::data(i*ref_layer_axes[0] + j) -
	      ref_atm_layer.pixel_array<double>::data(i*ref_layer_axes[0] + i);
	      this->pixel_array<double>::pixeldata[j-i] += tmp*tmp;
	      this->pixel_array<double>::pixelwts[j-i]++;
	  }
	}
      */

      } else {

        // As above, check for dimensional consistency with the ref_atm_layer
        vector<long> midpoint(2);
        if(!this->weights_allocated() || 
	   this->get_axes()!=ref_layer_axes){
	  for(int i=0; i<ref_layer_axes.size(); i++){
	    if(ref_layer_axes[i]%2==0) midpoint[i] = ref_layer_axes[i]/2-1;
	    else midpoint[i] = ref_layer_axes[i]/2;
	  }
	  this->set_axes(ref_layer_axes);
	  this->pixel_array<double>::allocate_weights(0);
        }

        double tmp;
        double midpoint_val = 
	  ref_atm_layer.pixel_array<double>::data(midpoint[1]*ref_layer_axes[0] + 
						midpoint[0]);
        for(int i=0; i<ref_layer_axes[1]; i++){
	  for(int j=0; j<ref_layer_axes[0]; j++){
	    tmp = 
	      ref_atm_layer.pixel_array<double>::data(i*ref_layer_axes[0] + j) -
	      midpoint_val;
	    this->pixel_array<double>::pixeldata[i*ref_layer_axes[0] + j] +=
	      tmp*tmp;
	    this->pixel_array<double>::pixelwts[i*ref_layer_axes[0] + j]++;
	  }
        }
      }
    }
  };
}

#endif
