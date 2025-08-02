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

#ifndef SHACK_HARTMANN_CENTROIDS_H
#define SHACK_HARTMANN_CENTROIDS_H

#include <iostream>
#include <vector>
#include "AO_sim_base.h"
#include "diffractive_wavefront.h"

namespace Arroyo {


  ///    
  /// A class to hold Shack Hartmann centroids
  ///

  class Shack_Hartmann_centroids :
    public AO_sim_base, 
    public pixel_array<double> {

    private:

    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("Shack Hartmann centroids"));};

    public:

    ///////////////////////////////////////////
    ///  Null constructor
    Shack_Hartmann_centroids(){};

    ///////////////////////////////////////////
    ///  Destructor
    ~Shack_Hartmann_centroids(){};

    ///////////////////////////////////////////
    ///  Copy constructor
    Shack_Hartmann_centroids(const Shack_Hartmann_centroids & shcentroids){
      this->operator=(shcentroids);
    };

    ///////////////////////////////////////////
    ///  Construct from file
    Shack_Hartmann_centroids(const char * filename){
      this->read(filename);
    };

    ///////////////////////////////////////////
    ///  Construct from an iofits object
    Shack_Hartmann_centroids(const iofits & iof){
      this->read(iof);
    };

    ///////////////////////////////////////////
    ///  Construct a null-initialized instance
    Shack_Hartmann_centroids(const vector<long> & lenslet_axes);

    ///////////////////////////////////////////
    ///  Construct from a wavefront
    template<class T>
    Shack_Hartmann_centroids(const vector<long> & lenslet_axes, const diffractive_wavefront<T> & dwf);

    ///////////////////////////////////////////
    ///  Construct from an as-yet-to-be-written
    /// simulated observation
    //Shack_Hartmann_centroids(const vector<long> & lenslet_axes, const simulated_observation & simobs);

    ///////////////////////////////////////////
    ///  Operator =
    Shack_Hartmann_centroids & operator=(const Shack_Hartmann_centroids & iznke);

    ///////////////////////////////////////////
    ///  read from a file
    void read(const char * filename);

    ///////////////////////////////////////////
    ///  read from an iofits object
    void read(const Arroyo::iofits & iof);

    ///////////////////////////////////////////
    ///  write to a file
    void write(const char * filename) const;

    ///////////////////////////////////////////
    ///  write to an iofits object
    void write(Arroyo::iofits & iof) const;

    ///////////////////////////////////////////
    ///  Function to print the coefficients
    void print(std::ostream & os, const char * prefix = "") const;

    static int verbose_level;

  };

  template<class T>
    Shack_Hartmann_centroids::Shack_Hartmann_centroids(const vector<long> & lenslet_axes,
						       const diffractive_wavefront<T> & dwf){


    if(lenslet_axes.size()!=2 || lenslet_axes[0]<=0 || lenslet_axes[1]<=0){
      cerr << "Shack_Hartmann_centroids::Shack_Hartmann_centroids error - "
	   << "invalid lenslet dimensions\n";
      throw(string("Shack_Hartmann_centroids::Shack_Hartmann_centroids"));      
    }

    vector<long> wf_axes = dwf.get_axes();
    if(lenslet_axes[0]*(wf_axes[0]/lenslet_axes[0])!=wf_axes[0] ||
       lenslet_axes[1]*(wf_axes[1]/lenslet_axes[1])!=wf_axes[1]){
      cerr << "Shack_Hartmann_centroids::Shack_Hartmann_centroids error - "
	   << "wavefront axes not evenly divisible by lenslet axes\n";
      cerr << "wavefront dimensions " << wf_axes[0] << "x" << wf_axes[1] << endl;
      cerr << "lenslet dimensions " << lenslet_axes[0] << "x" << lenslet_axes[1] << endl;
      throw(string("Shack_Hartmann_centroids::Shack_Hartmann_centroids"));
    }

    vector<long> pixarr_axes(lenslet_axes);
    pixarr_axes[1]*=2;

    this->set_axes(pixarr_axes);

    long xpix_per_lenslet = wf_axes[1]/lenslet_axes[1];
    long ypix_per_lenslet = wf_axes[0]/lenslet_axes[0];

    // The halfpixel information
    double x_halfpix=0, y_halfpix=0;
    int x_extrapix=1, y_extrapix=1;
    if(xpix_per_lenslet%2==0){
      x_halfpix = .5;
      x_extrapix = 0;
    }
    if(ypix_per_lenslet%2==0){
      y_halfpix = .5;
      y_extrapix = 0;
    }

    long nlenslets = lenslet_axes[0]*lenslet_axes[1];
    double xcentroid, ycentroid, total, intensity;
    for(int i=0; i<lenslet_axes[1]; i++){
      for(int j=0; j<lenslet_axes[0]; j++){
	xcentroid = ycentroid = total = 0;
	for(int k=0; k<xpix_per_lenslet; k++){
	  for(int l=0; l<ypix_per_lenslet; l++){
	    intensity = norm(dwf.data((i*xpix_per_lenslet+k)*wf_axes[0]+j*ypix_per_lenslet+l));
	    total+=intensity;
	    xcentroid += intensity*(k-xpix_per_lenslet/2+x_halfpix);
	    ycentroid += intensity*(l-ypix_per_lenslet/2+y_halfpix);
	  }
	}
	total==0 ? pixeldata[i*lenslet_axes[0]+j] = 0 :
	  pixeldata[i*lenslet_axes[0]+j] = xcentroid/total;
	total==0 ? pixeldata[i*lenslet_axes[0]+j+nlenslets] = 0 :
	  pixeldata[i*lenslet_axes[0]+j+nlenslets] = ycentroid/total;
	if(Shack_Hartmann_centroids::verbose_level)
	  cout << i << "\t" << j 
	       << "\t" << pixeldata[i*lenslet_axes[0]+j] 
	       << "\t" << pixeldata[i*lenslet_axes[0]+j+nlenslets] 
	       << endl;

      }
    }
  }
}

#endif
