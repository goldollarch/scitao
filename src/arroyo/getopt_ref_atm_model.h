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

#ifndef GETOPT_REF_ATM_MODEL_H
#define GETOPT_REF_ATM_MODEL_H
#include <unistd.h>
#include "Ellerbroek_Cerro_Pachon_model.h"
#include "Ellerbroek_Mauna_Kea_model.h"
#include "Hufnagel_Valley_model.h"
#include "SLCSAT_model.h"
#include "TMT_SRD_v13_Cn2_model.h"
#include "NGAO_system_design_model.h"
#include "Gemini_GLAO_study_model.h"
#include "Palomar_model.h"

namespace Arroyo {

  // A convenience function to encapsulate the help for all refractive atmospheric
  // models supported by Arroyo
  void usage_ref_atm_model();


  // A convenience function to encapsulate command line parsing for all refractive atmospheric
  // models supported by Arroyo
  refractive_atmospheric_model * parse_ref_atm_model(int argc, 
						     char * const argv[], 
						     three_frame tf);
} 
#endif
