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

#ifndef APERTURE_H
#define APERTURE_H

#include "pixel_array.h"
#include "optic.h"

namespace Arroyo {

  using std::string;
  using std::vector;
  using std::ostream;

  ///
  /// A virtual base class to represent an aperture.
  ///

  class aperture :
    virtual public plane_optic,
    virtual public one_to_one_optic {

    protected:

    /// \brief
    /// A flag to specify whether partial illumination of edge pixels
    /// should be calculated. Default is true.
    bool areal_weighting;

    public:

    ///////////////////////////////////////////
    ///  Null constructor
    aperture(){areal_weighting = true;};

    ///////////////////////////////////////////
    ///  Copy constructor
    aperture(const aperture & ap);

    ///////////////////////////////////////////
    ///  Virtual destructor
    virtual ~aperture(){};

    ///////////////////////////////////////////
    ///  Operator = 
    aperture & operator=(const aperture & ap);

    ///////////////////////////////////////////
    /// Virtual clone method
    ///
    /// Calling routine is responsible for deleting memory
    virtual aperture * clone() const = 0;

    //////////////////////////////////////////// 
    /// Query whether wavefront pixels incident
    /// on the edge of the aperture will be weighted
    /// according to the area of overlap.  The 
    /// alternative is to zero them or leave them
    /// unchanged according to whether the center
    /// of the pixel lies outside or inside the
    /// aperture
    bool get_areal_weighting() const {return(areal_weighting);};

    //////////////////////////////////////////// 
    /// Choose whether wavefront pixels incident
    /// on the edge of the aperture will be weighted
    /// according to the area of overlap.  The 
    /// alternative is to zero them or leave them
    /// unchanged according to whether the center
    /// of the pixel lies outside or inside the
    /// aperture
    void set_areal_weighting(bool aewtg) {areal_weighting = aewtg;};

    ///////////////////////////////////////////
    ///  Read from file
    virtual void read(const char * filename) = 0;

    ///////////////////////////////////////////
    ///  Read from an iofits object
    virtual void read(const iofits & iof) = 0;
 
    ///////////////////////////////////////////
    ///  Write to file
    virtual void write(const char * filename) const = 0;

    ///////////////////////////////////////////
    ///  Write to an iofits object
    virtual void write(iofits & iof) const = 0;

    ///////////////////////////////////////////
    ///  Print
    virtual void print(ostream & os, const char * prefix="") const;

    ///////////////////////////////////////////
    /// Virtual function to return the overlapping area between a
    /// convex polygon and the circular aperture.
    ///
    /// The vertices of the polygon must lie in the plane of the aperture,
    /// or this function throws an error
    virtual double convex_polygon_overlap(const vector<three_point> & polygon_vertices) const = 0;

    ///////////////////////////////////////////
    ///  Factory to construct apertures from file
    static aperture * aperture_factory(const char * filename);

    ///////////////////////////////////////////
    ///  Factory to construct apertures from file
    static aperture * aperture_factory(const iofits & iof);

  };

  ///
  /// A class to represent a circular aperture.
  ///

  class circular_aperture :
    public aperture {

    private:

    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("circular aperture"));};

    /// The diameter of the aperture in meters
    double diameter;

    /// A template member function to perform
    /// transform on both float and double
    /// instantiations of wavefront.  This 
    /// is necessary because there is no 
    /// mechanism in C++ for virtual template
    /// member functions.
    template<class T>
      void private_transform(diffractive_wavefront<T> & wf) const;

    public:

    ///////////////////////////////////////////
    ///  Null constructor
    circular_aperture();

    ///////////////////////////////////////////
    ///  Copy constructor
    circular_aperture(const circular_aperture & circ_ap);

    ///////////////////////////////////////////
    ///  Construct from file
    circular_aperture(const char * filename);

    ///////////////////////////////////////////
    ///  Construct from iofits object
    circular_aperture(const iofits & iof);

    ///////////////////////////////////////////r
    ///  Construct from the bits
    ///
    ///  The diameter should be specified in
    ///  meters
    circular_aperture(double in_diameter);

    ///////////////////////////////////////////
    ///  Destructor
    ~circular_aperture(){};

    ///////////////////////////////////////////
    ///  Operator = 
    circular_aperture & operator=(const circular_aperture & circ_ap);

    ///////////////////////////////////////////
    /// Clone method
    ///
    /// Calling routine is responsible for deleting memory
    circular_aperture * clone() const {
      return(new circular_aperture(*this));
    };

    ///////////////////////////////////////////
    ///  Read from file
    virtual void read(const char * filename);

    ///////////////////////////////////////////
    ///  Read from an iofits object
    virtual void read(const iofits & iof);
 
    ///////////////////////////////////////////
    ///  Write to file
    virtual void write(const char * filename) const;

    ///////////////////////////////////////////
    ///  Write to an iofits object
    virtual void write(iofits & iof) const;

    ///////////////////////////////////////////
    ///  Print
    void print(ostream & os, const char * prefix="") const;

    ///////////////////////////////////////////
    ///  Return the diameter of the aperture
    ///
    ///  Diameter returned in meters
    double get_diameter() const {return(diameter);};

    ///////////////////////////////////////////
    ///  Get a rectangular region guaranteed to cover
    ///  the aperture.  The resulting region will
    ///  have its edges aligned with the x and y axes
    ///  of the three_frame tf.  
    ///
    ///  If foreshortening is on, the projected
    ///  region is guaranteed to cover the optic
    ///
    ///  If the z axis of the three_frame is orthogonal
    ///  to the z axis of the aperture, this function
    ///  throws an error
    rectangular_region get_covering_region(const three_frame & tf) const;

    ///////////////////////////////////////////
    ///  Apply the aperture to the wavefront
    void transform(diffractive_wavefront<float> & wf) const;

    ///////////////////////////////////////////
    ///  Apply the aperture to the wavefront
    void transform(diffractive_wavefront<double> & wf) const;

    ///////////////////////////////////////////
    ///  Apply the aperture to a geometric_ray
    //void transform(geometric_ray & gray) const;

    ///////////////////////////////////////////
    ///  Apply the aperture to a geometric_wavefront
    //void transform(geometric_wavefront & gwf) const;

    ///////////////////////////////////////////
    /// Return the overlapping area between a convex polygon and the
    /// circular aperture.
    ///
    /// The vertices of the polygon must lie in the plane of the aperture,
    /// or this function throws an error
    double convex_polygon_overlap(const vector<three_point> & polygon_vertices) const;

  };

  ///
  /// A class to represent an annular aperture.
  ///

  class annular_aperture :
    public aperture {

    private:

    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("annular aperture"));};

    protected:

    /// The inner_diameter of the aperture in meters
    double inner_diameter;

    /// The outer_diameter of the aperture in meters
    double outer_diameter;

    /// A template member function to perform
    /// transform on both float and double
    /// instantiations of wavefront.  This 
    /// is necessary because there is no 
    /// mechanism in C++ for virtual template
    /// member functions.
    template<class T>
      void private_transform(diffractive_wavefront<T> & wf) const;

    public:

    ///////////////////////////////////////////
    ///  Null constructor
    annular_aperture();

    ///////////////////////////////////////////
    ///  Copy constructor
    annular_aperture(const annular_aperture & annular_ap);

    ///////////////////////////////////////////
    ///  Construct from file
    annular_aperture(const char * filename);

    ///////////////////////////////////////////
    ///  Construct from iofits object
    annular_aperture(const iofits & iof);

    ///////////////////////////////////////////
    ///  Construct from the bits
    ///
    ///  The inner and outer diameters should be
    ///  specified in meters
    annular_aperture(double in_diameter, double out_diameter);

    ///////////////////////////////////////////
    ///  Destructor
    ~annular_aperture(){};

    ///////////////////////////////////////////
    ///  Operator = 
    annular_aperture & operator=(const annular_aperture & annular_ap);

    ///////////////////////////////////////////
    /// Clone method
    ///
    /// Calling routine is responsible for deleting memory
    annular_aperture * clone() const {
      return(new annular_aperture(*this));
    };

    ///////////////////////////////////////////
    ///  Read from file
    virtual void read(const char * filename);

    ///////////////////////////////////////////
    ///  Read from an iofits object
    virtual void read(const iofits & iof);
 
    ///////////////////////////////////////////
    ///  Write to file
    virtual void write(const char * filename) const;

    ///////////////////////////////////////////
    ///  Write to an iofits object
    virtual void write(iofits & iof) const;

    ///////////////////////////////////////////
    ///  Print
    void print(ostream & os, const char * prefix="") const;

    ///////////////////////////////////////////
    ///  Return the inner diameter of the aperture
    ///
    ///  Diameter returned in meters
    double get_inner_diameter() const {return(inner_diameter);};

    ///////////////////////////////////////////
    ///  Return the outer diameter of the aperture
    ///
    ///  Diameter returned in meters
    double get_outer_diameter() const {return(outer_diameter);};

    ///////////////////////////////////////////
    ///  Get a rectangular region guaranteed to cover
    ///  the aperture.  The resulting region will
    ///  have its edges aligned with the x and y axes
    ///  of the three_frame tf.  
    ///
    ///  If foreshortening is on, the projected
    ///  region is guaranteed to cover the optic
    ///
    ///  If the z axis of the three_frame is orthogonal
    ///  to the z axis of the aperture, this function
    ///  throws an error
    rectangular_region get_covering_region(const three_frame & tf) const;

    ///////////////////////////////////////////
    ///  Apply the aperture to the wavefront
    virtual void transform(diffractive_wavefront<float> & wf) const;

    ///////////////////////////////////////////
    ///  Apply the aperture to the wavefront
    virtual void transform(diffractive_wavefront<double> & wf) const;

    ///////////////////////////////////////////
    ///  Apply the aperture to a geometric_ray
    //virtual void transform(geometric_ray & gray) const;

    ///////////////////////////////////////////
    ///  Apply the aperture to a geometric_wavefront
    //virtual void transform(geometric_wavefront & gwf) const;

    ///////////////////////////////////////////
    /// Return the overlapping area between a convex polygon and the
    /// annular aperture
    ///
    /// The vertices of the polygon must lie in the plane of the aperture,
    /// or this function throws an error
    double convex_polygon_overlap(const vector<three_point> & polygon_vertices) const;

  };

  ///
  /// A class to represent a rectangular aperture.
  ///

  class rectangular_aperture :
    public aperture {

    private:

    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("rectangular aperture"));};

    protected:

    /// The physical dimensions of the aperture
    vector<double> size;

    /// A template member function to perform
    /// transform on both float and double
    /// instantiations of wavefront.  This 
    /// is necessary because there is no 
    /// mechanism in C++ for virtual template
    /// member functions.
    template<class T>
      void private_transform(diffractive_wavefront<T> & wf) const;

    public:

    ///////////////////////////////////////////
    ///  Null constructor
    rectangular_aperture();

    ///////////////////////////////////////////
    ///  Copy constructor
    rectangular_aperture(const rectangular_aperture & rectangular_ap);

    ///////////////////////////////////////////
    ///  Construct from file
    rectangular_aperture(const char * filename);

    ///////////////////////////////////////////
    ///  Construct from iofits object
    rectangular_aperture(const iofits & iof);

    ///////////////////////////////////////////
    ///  Construct from the bits
    ///
    ///  The dimensions of the aperture should be
    ///  specified in meters
    rectangular_aperture(double x_size, double y_size);

    ///////////////////////////////////////////
    ///  Destructor
    ~rectangular_aperture(){};

    ///////////////////////////////////////////
    ///  Operator = 
    rectangular_aperture & operator=(const rectangular_aperture & rectangular_ap);

    ///////////////////////////////////////////
    /// Clone method
    ///
    /// Calling routine is responsible for deleting memory
    rectangular_aperture * clone() const {
      return(new rectangular_aperture(*this));
    };

    ///////////////////////////////////////////
    ///  Read from file
    virtual void read(const char * filename);

    ///////////////////////////////////////////
    ///  Read from an iofits object
    virtual void read(const iofits & iof);
 
    ///////////////////////////////////////////
    ///  Write to file
    virtual void write(const char * filename) const;

    ///////////////////////////////////////////
    ///  Write to an iofits object
    virtual void write(iofits & iof) const;

    ///////////////////////////////////////////
    ///  Print
    void print(ostream & os, const char * prefix="") const;

    ///////////////////////////////////////////
    ///  Return the physical dimensions of the aperture
    ///
    ///  Size returned in meters
    vector<double> get_size() const {return(size);};

    ///////////////////////////////////////////
    ///  Get a rectangular region guaranteed to cover
    ///  the aperture.  The resulting region will
    ///  have its edges aligned with the x and y axes
    ///  of the three_frame tf.  
    ///
    ///  If foreshortening is on, the projected
    ///  region is guaranteed to cover the optic
    ///
    ///  If the z axis of the three_frame is orthogonal
    ///  to the z axis of the aperture, this function
    ///  throws an error
    rectangular_region get_covering_region(const three_frame & tf) const;

    ///////////////////////////////////////////
    ///  Apply the aperture to the wavefront
    virtual void transform(diffractive_wavefront<float> & wf) const;

    ///////////////////////////////////////////
    ///  Apply the aperture to the wavefront
    virtual void transform(diffractive_wavefront<double> & wf) const;

    ///////////////////////////////////////////
    ///  Apply the aperture to a geometric_ray
    //virtual void transform(geometric_ray & gray) const;

    ///////////////////////////////////////////
    ///  Apply the aperture to a geometric_wavefront
    //virtual void transform(geometric_wavefront & gwf) const;

    ///////////////////////////////////////////
    /// Return the overlapping area between a convex polygon and the
    /// rectangular aperture
    ///
    /// The vertices of the polygon must lie in the plane of the aperture,
    /// or this function throws an error
    double convex_polygon_overlap(const vector<three_point> & polygon_vertices) const;

  };

  ///
  /// A class to represent an annular aperture with spiders
  /// of arbitrary width.
  ///

  class spidered_annular_aperture :
    public annular_aperture {

    private:

    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("spidered annular aperture"));};

    protected:

    /// The number of spiders. 
    long nspiders;

    /// The width of the spiders. 
    double spider_width;

    /// A template member function to perform
    /// transform on both float and double
    /// instantiations of wavefront.  This 
    /// is necessary because there is no 
    /// mechanism in C++ for virtual template
    /// member functions.
    template<class T>
      void private_transform(diffractive_wavefront<T> & wf) const;

    public:

    ///////////////////////////////////////////
    ///  Null constructor
    spidered_annular_aperture();

    ///////////////////////////////////////////
    ///  Copy constructor
    spidered_annular_aperture(const spidered_annular_aperture & spdrd_annular_ap);

    ///////////////////////////////////////////
    ///  Construct from file
    spidered_annular_aperture(const char * filename);

    ///////////////////////////////////////////
    ///  Construct from iofits object
    spidered_annular_aperture(const iofits & iof);

    ///////////////////////////////////////////
    ///  Construct from the bits
    ///
    ///  The inner and outer diameters and the width
    ///  of the spiders should be specified in meters
    spidered_annular_aperture(double in_diameter, double out_diameter, 
			      int nspiders, double spider_width);

    ///////////////////////////////////////////
    ///  Destructor
    ~spidered_annular_aperture(){};

    ///////////////////////////////////////////
    ///  Operator = 
    spidered_annular_aperture & operator=(const spidered_annular_aperture & spdrd_annular_ap);

    ///////////////////////////////////////////
    /// Clone method
    ///
    /// Calling routine is responsible for deleting memory
    spidered_annular_aperture * clone() const {
      return(new spidered_annular_aperture(*this));
    };

    ///////////////////////////////////////////
    ///  Read from file
    void read(const char * filename);

    ///////////////////////////////////////////
    ///  Read from an iofits object
    void read(const iofits & iof);
 
    ///////////////////////////////////////////
    ///  Write to file
    void write(const char * filename) const;

    ///////////////////////////////////////////
    ///  Write to an iofits object
    void write(iofits & iof) const;

    ///////////////////////////////////////////
    ///  Print
    void print(ostream & os, const char * prefix="") const;

    ///////////////////////////////////////////
    ///  Get a rectangular region guaranteed to cover
    ///  the aperture.  The resulting region will
    ///  have its edges aligned with the x and y axes
    ///  of the three_frame tf.  
    ///
    ///  If foreshortening is on, the projected
    ///  region is guaranteed to cover the optic
    ///
    ///  If the z axis of the three_frame is orthogonal
    ///  to the z axis of the aperture, this function
    ///  throws an error
    //rectangular_region get_covering_region(const three_frame & tf) const;
    
    ///////////////////////////////////////////
    ///  Apply the aperture to the wavefront
    void transform(diffractive_wavefront<float> & wf) const;

    ///////////////////////////////////////////
    ///  Apply the aperture to the wavefront
    void transform(diffractive_wavefront<double> & wf) const;

    ///////////////////////////////////////////
    ///  Apply the aperture to a geometric_ray
    //virtual void transform(geometric_ray & gray) const;

    ///////////////////////////////////////////
    ///  Apply the aperture to a geometric_wavefront
    //virtual void transform(geometric_wavefront & gwf) const;

    ///////////////////////////////////////////
    /// Return the overlapping area between a convex polygon and the
    /// spidered annular aperture
    ///
    /// The vertices of the polygon must lie in the plane of the aperture,
    /// or this function throws an error
    double convex_polygon_overlap(const vector<three_point> & polygon_vertices) const;

  };

  ///
  /// A class to represent a hexagonal aperture.
  ///     

  class hexagonal_aperture :
    public aperture {

    private:

    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("hexagonal aperture"));};

    /// A template member function to perform
    /// transform on both float and double
    /// instantiations of wavefront.  This 
    /// is necessary because there is no 
    /// mechanism in C++ for virtual template
    /// member functions.
    template<class T>
      void private_transform(diffractive_wavefront<T> & wf) const;

    protected:

    /// The length of an edge
    double edge_length;

    public:

    ///////////////////////////////////////////
    ///  Null constructor
    hexagonal_aperture();

    ///////////////////////////////////////////
    ///  Copy constructor
    hexagonal_aperture(const hexagonal_aperture & hexagonal_ap);

    ///////////////////////////////////////////
    ///  Construct from file
    hexagonal_aperture(const char * filename);

    ///////////////////////////////////////////
    ///  Construct from iofits object
    hexagonal_aperture(const iofits & iof);

    ///////////////////////////////////////////
    ///  Construct from the bits
    ///
    ///  The edge length should be specified in 
    ///  meters
    hexagonal_aperture(double in_edge_length);

    ///////////////////////////////////////////
    ///  Destructor
    ~hexagonal_aperture(){};

    ///////////////////////////////////////////
    ///  Operator = 
    hexagonal_aperture & operator=(const hexagonal_aperture & hexagonal_ap);

    ///////////////////////////////////////////
    /// Clone method
    ///
    /// Calling routine is responsible for deleting memory
    hexagonal_aperture * clone() const {
      return(new hexagonal_aperture(*this));
    };

    ///////////////////////////////////////////
    ///  Read from file
    virtual void read(const char * filename);

    ///////////////////////////////////////////
    ///  Read from an iofits object
    virtual void read(const iofits & iof);
 
    ///////////////////////////////////////////
    ///  Write to file
    virtual void write(const char * filename) const;

    ///////////////////////////////////////////
    ///  Write to an iofits object
    virtual void write(iofits & iof) const;

    ///////////////////////////////////////////
    ///  Print
    void print(ostream & os, const char * prefix="") const;

    ///////////////////////////////////////////
    ///  Return the edge length
    ///
    ///  Length returned in meters
    double get_edge_length() const {return(edge_length);};

    ///////////////////////////////////////////
    ///  Get a rectangular region guaranteed to cover
    ///  the aperture.  The resulting region will
    ///  have its edges aligned with the x and y axes
    ///  of the three_frame tf.  
    ///
    ///  If foreshortening is on, the projected
    ///  region is guaranteed to cover the optic
    ///
    ///  If the z axis of the three_frame is orthogonal
    ///  to the z axis of the aperture, this function
    ///  throws an error
    rectangular_region get_covering_region(const three_frame & tf) const;

    ///////////////////////////////////////////
    ///  Apply the aperture to the wavefront
    virtual void transform(diffractive_wavefront<float> & wf) const;

    ///////////////////////////////////////////
    ///  Apply the aperture to the wavefront
    virtual void transform(diffractive_wavefront<double> & wf) const;

    ///////////////////////////////////////////
    ///  Apply the aperture to a geometric_ray
    //virtual void transform(geometric_ray & gray) const;

    ///////////////////////////////////////////
    ///  Apply the aperture to a geometric_wavefront
    //virtual void transform(geometric_wavefront & gwf) const;

    ///////////////////////////////////////////
    /// Return the overlapping area between a convex polygon and the
    /// hexagonal aperture
    ///
    /// The vertices of the polygon must lie in the plane of the aperture,
    /// or this function throws an error
    double convex_polygon_overlap(const vector<three_point> & polygon_vertices) const;

  };


  /// 
  /// A class to represent a tiled hexagonalagonal aperture
  /// 

  class tiled_hexagonal_aperture :
    public aperture {

    private:
 
    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("tiled hexagonal aperture"));};

    ///////////////////////////////////////////
    ///  Null constructor - disabled because we require initialization
    ///  using a prototype hexagonal aperture
    tiled_hexagonal_aperture(){};

    protected:

    /// A 2d map of zeroes or ones to indicate whether tiles exist
    Arroyo::pixel_array<long> tilemap;

    /// The edge length of a tile, in meters
    double edge_length;

    /// The gap between tiles, in meters
    double gap_size;

    /// template member function to perform
    /// transform on both float and double
    /// instantiations of wavefront.  This 
    /// is necessary because there is no 
    /// mechanism in C++ for virtual template
    /// member functions.
    template<class T>
      void private_transform(diffractive_wavefront<T> & dwf) const;

    public:

    ///////////////////////////////////////////
    ///  Copy constructor
    tiled_hexagonal_aperture(const tiled_hexagonal_aperture & tiled_hexagonal_ap);

    ///////////////////////////////////////////
    ///  Construct from file
    tiled_hexagonal_aperture(const char * filename);

    ///////////////////////////////////////////
    ///  Construct from iofits object
    tiled_hexagonal_aperture(const iofits & iof);

    ///////////////////////////////////////////
    ///  Construct an aperture like at Keck,
    ///  with rings of hexagonal apertures.
    ///  The central hexagon is labelled as ring
    ///  zero.  So for the Keck primary,
    ///  inner_ring = 1 and outer_ring = 3
    /// 
    ///  According to Mitch, the actual glass 
    ///  hexagon has a 90 cm edge.  But the outer
    ///  1 mm is beveled, so is optically opaque.
    ///  Then there's a 5 mm gap between hexes.
    ///  So its like having a 7 mm gap, with a hexagon
    ///  that has an edge length of 
    ///  90 cm - 2*(.1 cm) /sqrt(3) = 89.88453 cm
    ///
    ///  Note - edge length and gap size must be 
    ///  positive.  Both are in meters.
    //tiled_hexagonal_aperture(int inner_ring, int outer_ring, 
    //double in_edge_length, double in_gap_size);

    ///////////////////////////////////////////
    ///  Construct an aperture like that proposed
    ///  for Celt, where hexagons tile an annulus.
    ///  All hexagons with centers that lie at radii
    ///  between inner_diameter and outer_diameter
    ///  are included in the aperture.
    ///  According to Mitch, one proposal for celt:
    ///  4mm gaps, hex side .5 meters, maybe 20 rings.
    ///
    ///  Note - edge length and gap size must be 
    ///  positive.  
    ///
    ///  All arguments are in meters
    tiled_hexagonal_aperture(double inner_diameter, double outer_diameter,
			     double in_edge_length, double in_gap_size);

    ///////////////////////////////////////////
    ///  Destructor
    ~tiled_hexagonal_aperture(){};

    ///////////////////////////////////////////
    ///  Operator = 
    tiled_hexagonal_aperture & operator=(const tiled_hexagonal_aperture & tiled_hexagonal_ap);
 
    ///////////////////////////////////////////
    /// Clone method
    ///
    /// Calling routine is responsible for deleting memory
    tiled_hexagonal_aperture * clone() const {
      return(new tiled_hexagonal_aperture(*this));
    };

    ///////////////////////////////////////////
    ///  Read from file
    virtual void read(const char * filename);

    ///////////////////////////////////////////
    ///  Read from an iofits object
    virtual void read(const iofits & iof);
 
    ///////////////////////////////////////////
    ///  Write to file
    virtual void write(const char * filename) const;

    ///////////////////////////////////////////
    ///  Write to an iofits object
    virtual void write(iofits & iof) const;

    ///////////////////////////////////////////
    ///  Print
    void print(ostream & os, const char * prefix="") const;

    ///////////////////////////////////////////
    ///  Return the edge length in meters
    double get_edge_length() const {return(edge_length);};

    ///////////////////////////////////////////
    ///  Return the gap between tiles, in meters
    double get_gap_size() const {return(gap_size);};

    ///////////////////////////////////////////
    ///  Get a rectangular region guaranteed to cover
    ///  the aperture.  The resulting region will
    ///  have its edges aligned with the x and y axes
    ///  of the three_frame tf.  
    ///
    ///  If foreshortening is on, the projected
    ///  region is guaranteed to cover the optic
    ///
    ///  If the z axis of the three_frame is orthogonal
    ///  to the z axis of the aperture, this function
    ///  throws an error
    rectangular_region get_covering_region(const three_frame & tf) const;

    ///////////////////////////////////////////
    ///  Apply the aperture to the wavefront
    virtual void transform(diffractive_wavefront<float> & dwf) const;

    ///////////////////////////////////////////
    ///  Apply the aperture to the wavefront
    virtual void transform(diffractive_wavefront<double> & dwf) const;

    ///////////////////////////////////////////
    /// Return the overlapping area between a convex polygon and the
    /// tiled hexagonal aperture
    ///
    /// The vertices of the polygon must lie in the plane of the aperture,
    /// or this function throws an error
    double convex_polygon_overlap(const vector<three_point> & polygon_vertices) const;
  };
}

#endif
