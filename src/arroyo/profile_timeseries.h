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

#ifndef PROFILE_TIMESERIES_H
#define PROFILE_TIMESERIES_H

#include <iostream>
#include <fstream>
#include <ctime>
#include "refractive_atmosphere.h"
#include "time_val.h"

namespace Arroyo {

  class profile_timeseries {

  protected:

    std::vector<time_t> profile_timestamps;
    std::vector<std::vector<double> > layer_heights_meters;
    std::vector<std::vector<double> > layer_Cn2_coeffs;
  
  public:
  
    ///////////////////////////////////////////
    ///  NULL constructor
    profile_timeseries(){};

    ///////////////////////////////////////////
    ///  Copy constructor
    profile_timeseries(const profile_timeseries & pts) {
      this->operator=(pts);
    };

    ///////////////////////////////////////////
    ///  Construct from the bits
    profile_timeseries(const std::vector<time_t> & profile_timestamps,
		       const std::vector<std::vector<double> > & layer_heights_meters,
		       const std::vector<std::vector<double> > & layer_Cn2_coeffs);

    ///////////////////////////////////////////
    ///  Construct from dimm and mass files
    profile_timeseries(const char * dimm_filename,
		       const char * mass_filename);

    ///////////////////////////////////////////
    ///  Destructor
    ~profile_timeseries(){};

    ///////////////////////////////////////////
    ///  Operator =
    profile_timeseries & operator=(const profile_timeseries & pts);

    ///////////////////////////////////////////
    ///  Print
    void print(std::ostream & os, const char * prefix = "") const;

    ///////////////////////////////////////////
    ///  Get number of profiles
    int get_nprofiles() const {
      return(this->profile_timestamps.size());
    };

    ///////////////////////////////////////////
    ///  Get timestamp of a particular profile
    time_t get_timestamp(int index) const;

    ///////////////////////////////////////////
    ///  Return a refractive atmospheric model at a particular time
    refractive_atmospheric_model get_refractive_atmospheric_model(const time_t & timestamp) const;
    
    ///////////////////////////////////////////
    ///  Return a refractive atmospheric model for a particular index
    refractive_atmospheric_model get_refractive_atmospheric_model(int index) const;
    
    ///////////////////////////////////////////
    ///  Return a refractive atmospheric model averaged
    ///  over a period of time.
    ///
    ///  The function computes the average of all profiles in the
    ///  timeseries with timestamps >= start and <= end, and returns
    ///  this profile in ave_ref_atm_model.  The function returns the
    ///  number of profiles in the average.
    ///
    ///  If there are no profiles in this time period, the function
    ///  throws an error.
    refractive_atmospheric_model get_refractive_atmospheric_model(const time_t & start,
								  const time_t & end,
								  int & nprofiles_in_average) const;
    
  };
}
#endif
