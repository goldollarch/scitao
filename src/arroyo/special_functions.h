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

#ifndef SPECIAL_FUNCTIONS_H
#define SPECIAL_FUNCTIONS_H

#include <cmath>
#include <vector>

namespace Arroyo {

  ///////////////////////////////////////////////
  /// Returns the gamma function evaluated at arg
  /// 
  /// This function uses the series expansion from
  /// Abramowitz and Stegun, section 6.1.34
  ///
  /// This function throws an error if arg is
  /// a negative integer or is equal to zero
  double gamma_function(double arg);

  ///////////////////////////////////////////////
  /// Returns the generalized hypergeometric function
  /// p_F_q evaluated at arg.  Here p and q are the
  /// sizes of the vectors numerator_vals and denominator_vals,
  /// respectively.  Specifically, this function computes
  ///
  /// p_F_q(arg) = sum over k=0 to infinity of 
  ///
  /// (a1,k*a2,k*...*ap,k)/(b1,k*b2,k*...*bq,k)*(arg^{k}/k!)
  ///
  /// where the factors an,k = gamma_function(n+k)/gamma_function(n)
  /// are the elements of the vector numerator_vals.  Similarly,
  /// the factors bn,k are the elements of the vector denominator_vals
  /// and are likewise defined.
  ///
  /// Convergence criteria:  according to Sasiela section 1.3,
  ///
  /// p < q + 1  converges for all arg values
  /// p = q + 1  converges for abs(arg) <= 1 
  /// p > q + 1  does not converge
  ///
  /// See Sasiela eqn 1.3.1 for a full definition
  ///
  /// The series will be evaluated until successive contributions
  /// fall below the error tolerance.
  ///
  /// This function throws an error if abs(arg) >= 1
  double generalized_hypergeometric_function(double arg,
  			const std::vector<double> & numerator_vals,
			const std::vector<double> & denominator_vals,
			double error_tolerance = 1e-10);

  ///////////////////////////////////////////////
  /// Returns the Bessel function of the first 
  /// kind for arbitrary order evaluated at arg
  ///
  /// This function uses the representation in
  /// terms of generalized hypergeometric functions
  /// from Sasiela Table 1.2
  double bessel_Jnu(double order, double arg);

  ///////////////////////////////////////////////
  /// Returns the Bessel function of the third 
  /// kind for arbitrary order evaluated at arg
  ///
  /// This function uses the representation in
  /// terms of generalized hypergeometric functions.
  double bessel_Knu(double order, double arg);

  ///////////////////////////////////////////////
  /// Returns the numerical constant
  ///
  /// (1/5)*(1/2)^{7/3}*(gamma_function(1/6))^{2}/gamma_function(1/3)
  ///
  /// This constant appears frequently in phase covariance
  /// calculations
  double get_Xi();
}

#endif
