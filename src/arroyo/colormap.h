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

#ifndef COLORMAP_H
#define COLORMAP_H

#include <vector>
#include <ostream>

namespace Arroyo {

  ///
  /// An abstract base class to hold a colormap
  ///

  class colormap {

  public:
    ///////////////////////////////////////////
    ///  Null constructor
    colormap(){};

    ///////////////////////////////////////////
    ///  Trivial destructor
    virtual ~colormap(){};

    ///////////////////////////////////////////
    ///  Operator=
    colormap & operator= (const colormap & fhd);
 
    ///////////////////////////////////////////
    ///  Read from file
    virtual void read(const char * filename) = 0;

    ///////////////////////////////////////////
    ///  Print
    virtual void print(std::ostream & os, const char * prefix = "") const = 0;

    ///////////////////////////////////////////
    ///  Invert the colormap
    virtual void invert() = 0;

    ///////////////////////////////////////////
    ///  Get R value - 0 to 255
    virtual char get_R(double val, double min, double max,
				bool logscale = false) const = 0;

    ///////////////////////////////////////////
    ///  Get G value - 0 to 255
    virtual char get_G(double val, double min, double max,
				bool logscale = false) const = 0;

    ///////////////////////////////////////////
    ///  Get B value - 0 to 255
    virtual char get_B(double val, double min, double max,
				bool logscale = false) const = 0;

  };

  ///
  /// A class to hold an sao colormap.  This type of 
  /// colormap is defined at
  /// http://tdc-www.harvard.edu/software/saoimage/saoimage.color.html#cmap
  ///  
  /// This class does not yet support the GAMMA entry, and will throw
  /// an error if this entry appears in the colormap
  ///

  class sao_colormap :
    public colormap {

  private:

  protected:
    
    std::vector<std::vector<double> > R_points;
    std::vector<std::vector<double> > G_points;
    std::vector<std::vector<double> > B_points;

  public:

    ///////////////////////////////////////////
    ///  Null constructor - initialize to greyscale
    sao_colormap(){
      R_points = std::vector<std::vector<double> >(3);
      R_points[0] = std::vector<double>(2,0);
      R_points[1] = std::vector<double>(2,1);
      R_points[1][0] = .5;
      R_points[2] = std::vector<double>(2,1);

      G_points = R_points;
      B_points = R_points;

    };

    ///////////////////////////////////////////
    ///  Construct from the bits
    sao_colormap(std::vector<std::vector<double> > & in_R_points, 
		 std::vector<std::vector<double> > & in_G_points, 
		 std::vector<std::vector<double> > & in_B_points);

    ///////////////////////////////////////////
    ///  Construct from file
    sao_colormap(const char * filename);

    ///////////////////////////////////////////
    ///  Copy constructor
    sao_colormap(const sao_colormap & scm){
      sao_colormap::operator=(scm);
    };

    ///////////////////////////////////////////
    ///  Trivial destructor
    ~sao_colormap(){};

    ///////////////////////////////////////////
    ///  Operator=
    sao_colormap & operator= (const sao_colormap & scm);
 
    ///////////////////////////////////////////
    ///  Read from file
    void read(const char * filename);

    ///////////////////////////////////////////
    ///  Print
    void print(std::ostream & os, const char * prefix = "") const;

    ///////////////////////////////////////////
    ///  Invert the colormap
    void invert();

    ///////////////////////////////////////////
    ///  Get R value - 0 to 255
    char get_R(double val, double min, double max, bool logscale) const;

    ///////////////////////////////////////////
    ///  Get G value - 0 to 255
    char get_G(double val, double min, double max, bool logscale) const;

    ///////////////////////////////////////////
    ///  Get B value - 0 to 255
    char get_B(double val, double min, double max, bool logscale) const;

  };

}

#endif
