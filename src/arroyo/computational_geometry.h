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

#ifndef COMPUTATIONAL_GEOMETRY_H
#define COMPUTATIONAL_GEOMETRY_H

#include <cmath>
#include <vector>
#include <string> 
#include "AO_cpp.h"
#include "three_frame.h"

namespace Arroyo {

  ///////////////////////////////////////////////
  // Find the point of intersection between two line segments defined
  // by their endpoints a1, a2 and b1, b2.  The result is contained in
  // the vector returned by this function, which may have size 0, 1, or
  // 2.
  //
  //  0             No intersection
  //  1             Intersection is a single point
  //  2             Intersection is a line segment
  //
  // The last result may only occur if the line segments are collinear.
  //
  // WARNING:  this function can return inaccurate results due to 
  // numerical roundoff.  Quantifying the limits of validity of this
  // function is a work in progress.
  std::vector<Arroyo::three_point>
  get_line_segment_intersection(const three_point & a1, 
				const three_point & a2, 
				const three_point & b1,
				const three_point & b2);


  ///////////////////////////////////////////////
  /// Return the point of intersection between a ray extending from
  /// o_a in the direction n_a, and a ray extending from o_b in the
  /// direction n_b.
  //
  // WARNING:  this function can return inaccurate results due to 
  // numerical roundoff.  Quantifying the limits of validity of this
  // function is a work in progress.
  Arroyo::three_point get_ray_ray_intersection(
  				const Arroyo::three_point & o_a, 
				const Arroyo::three_vector & n_a,
				const Arroyo::three_point & o_b, 
				const Arroyo::three_vector & n_b);

  ///////////////////////////////////////////////
  /// Return the point of intersection between a ray extending from
  /// o_a in the direction n_a and a plane defined by normal n_b and
  /// origin o_b.
  //
  // WARNING:  this function can return inaccurate results due to 
  // numerical roundoff.  Quantifying the limits of validity of this
  // function is a work in progress.
  Arroyo::three_point get_ray_plane_intersection(
				const Arroyo::three_point & o_a, 
				const Arroyo::three_vector & n_a,
				const Arroyo::three_point & o_b, 
				const Arroyo::three_vector & n_b);

  ///////////////////////////////////////////////
  // Find the directed area of a polygon specified by its vertices.
  // The vector polygon_vertices must have size greater than or equal
  // to three, and no three vertices may be collinear.  The vertices
  // must be sorted cyclically, with the function returning a positive
  // value for counterclockwise sorting and a negative value for
  // clockwise sorting.
  //
  // WARNING:  this function can return inaccurate results due to 
  // numerical roundoff.  Quantifying the limits of validity of this
  // function is a work in progress.
  double get_area_of_polygon(
		const std::vector<Arroyo::three_point> & polygon_vertices);

  ///////////////////////////////////////////////
  /// Modified from "Ray Tracing News Volume 5, Number 3.  Some explanation
  /// also appears on page 117 of "Object oriented ray tracing in C++" by
  /// Nicholas Wilt
  /// 
  /// Shoot a test ray along +X axis and count the crossings Note: the
  /// points in vtp are interpreted as ordered vertices of the polygon.
  /// They may be ordered clockwise or counterclockwise
  ///
  /// If the point lies at a distance less than three_frame::precision from
  /// a vertex or a segment edge of the polygon, it is defined as lying within
  /// the polygon.  These two cases are indicated by the value of the booleans 
  /// passed to the function
  ///
  /// WARNING:  this function can return inaccurate results due to 
  /// numerical roundoff.  Quantifying the limits of validity of this
  /// function is a work in progress.
  bool point_within_polygon(const three_point & tp, 
			    const vector<three_point> & vertices,
			    bool & point_on_edge,
			    bool & point_on_vertex);
    
  ///////////////////////////////////////////////
  /// Same as above, but no degeneracy information
  bool point_within_polygon(
		const three_point & tp, const vector<three_point> & vtp);

  ///////////////////////////////////////////////
  /// Partially based on "Computational Geometry in C" by O'Rourke
  /// Section 7.4.
  ///
  /// This function returns the vertices of a convex polygon formed
  /// from the intersection of two convex polygons, each specified by
  /// their vertices.  The first_polygon_vertices and
  /// second_polygon_vertices must each contain at least three unique
  /// three_points, no three of which may be collinear.  The vertices
  /// of each polygon must be ordered in the same sense of rotation.
  /// If either of these conditions are violated, this function throws
  /// an error.
  ///
  /// WARNING:  this function can return inaccurate results due to 
  /// numerical roundoff.  Quantifying the limits of validity of this
  /// function is a work in progress.
  vector<three_point> get_convex_polygon_intersection(
  			const vector<three_point> & first_polygon_vertices,
			const vector<three_point> & second_polygon_vertices,
			int verbose = 0);

  ///////////////////////////////////////////////
  // Find the point of intersection between a line segment defined by
  // its endpoints a1 and a2 and a circle with origin circle_origin
  // and the given radius.  The result is contained in the vector of
  // three_points returned by this function, which may have size 0, 1,
  // or 2.
  //
  //  0             No intersection
  //  1             Intersection is a single point
  //  2             Intersection contains two points
  //
  // WARNING:  this function can return inaccurate results due to 
  // numerical roundoff.  Quantifying the limits of validity of this
  // function is a work in progress.
  vector<three_point> get_line_segment_circle_intersection(
				const three_point & a1, 
				const three_point & a2, 
				const three_point & circle_origin,
				const double & radius, 
				vector<bool> & intersection_point_on_boundary,
				int verbose=0);


  ///////////////////////////////////////////////
  /// Find the overlap between a convex polygon with vertices
  /// polygon_vertices and a circular region with the given radius.
  /// The origin of the reference frame is taken to be the center of
  /// the circle, which is assumed to be in the x-y plane of this
  /// reference frame.  The polygon vertices must lie in this plane or
  /// this function throws an error.  If the number of vertices in the
  /// polygon is less than three, an error is thrown.
  /// 
  /// This function is meant to apply to the case when the polygon is
  /// much smaller than the circle, so that we can assume that
  /// polygons which overlap the circular border are divided by a
  /// chord of the circle rather than an arc.  The function returns
  /// the vertices of the resulting convex polygon.
  vector<three_point> get_convex_polygon_circle_intersection(
			const vector<three_point> & polygon_vertices,
			const three_frame & circle_tf, 
			double radius);


  ///////////////////////////////////////////////
  ///  This function analytically integrates x and y over a region of
  ///  support defined by the convex polygon.  The three frame tf is
  ///  used to define x and y, and the origin of this frame defines
  ///  the zero point.  All vertices of the polygon must lie in the
  ///  x-y plane of this three frame.  The results of the integration
  ///  are returned in x_intgrl and y_intgrl, respectively.
  ///
  ///  This function is used to compute elements of the geometry matrix for
  ///  pyramidal actuator influence functions
  template<class T>
    void convex_polygon_integration(const three_frame & tf,
			    const vector<three_point> & polygon_vertices,
			    T & x_intgrl,
			    T & y_intgrl,
			    int verbose = 0){
    try{
      x_intgrl = y_intgrl = 0;
      int nvertices = polygon_vertices.size();

      // If the polygon is degenerate, return zero for the integrand values
      if(nvertices<=2) return;

      double slope_b, slope_c;
      double intercept_b, intercept_c;
      double x_a, y_a, x_b, y_b, x_c, y_c;
      double last_x_b, last_y_b, last_x_c, last_y_c;
      double xmin, ymin, xmax, ymax;
      double sign;
      int index;
      double tmp;
      
      /////////////////////////////
      //  Perform the x integral //
      /////////////////////////////
      // Find the point with smallest x coordinate
      index = 0;
      x_a = polygon_vertices[0].x(tf);
      for(int i=1; i<nvertices; i++){
	tmp = polygon_vertices[i].x(tf);
	if(tmp<x_a){
	  index = i;
	  x_a = tmp;
	}
      }
      int indexb=index, indexc=index;

      do {
	indexb = (indexb+1)%nvertices;
	x_b = polygon_vertices[indexb].x(tf);
      } while((x_b-x_a)<three_frame::precision);
    
      do {
	indexc = (indexc-1+nvertices)%nvertices;
	x_c = polygon_vertices[indexc].x(tf);
      } while((x_c-x_a)<three_frame::precision);

      sign = 1;
      if(polygon_vertices[indexc].y(tf)>polygon_vertices[indexb].y(tf))
	sign = -1;

      if(verbose) cerr << "x input indices " << indexb
		       << "\t" << indexc << endl;

      while(1){
	if(indexb==indexc) break;

	y_b = polygon_vertices[indexb].y(tf);
	last_x_b = polygon_vertices[(indexb-1+nvertices)%nvertices].x(tf);
	last_y_b = polygon_vertices[(indexb-1+nvertices)%nvertices].y(tf);

	slope_b = (y_b-last_y_b)/(x_b-last_x_b);
	intercept_b = y_b - slope_b*x_b;

	y_c = polygon_vertices[indexc].y(tf);
	last_x_c = polygon_vertices[(indexc+1)%nvertices].x(tf);
	last_y_c = polygon_vertices[(indexc+1)%nvertices].y(tf);

	slope_c = (y_c-last_y_c)/(x_c-last_x_c);
	intercept_c = y_c - slope_c*x_c;

	xmin = last_x_b > last_x_c ? last_x_b : last_x_c;
	xmax = x_b < x_c ? x_b : x_c;
	ymin = y_b;
	ymax = y_c;
	if(y_b>y_c){
	  ymin = y_c;
	  ymax = y_b;
	}
      
	x_intgrl += sign*((slope_b-slope_c)*(xmax*xmax*xmax - xmin*xmin*xmin)/3.0 +
			  (intercept_b - intercept_c)*(xmax*xmax - xmin*xmin)/2.0);

	if(verbose){
	  cerr << "\tindices " << indexb << "\t" << indexc << endl; 
	  cerr << "\tlast indices " << (indexb-1+nvertices)%nvertices
	       << "\t" << (indexc+1)%nvertices << endl;
	  cerr << "\tb coords " << x_b << ", " << y_b << "\t"
	       << " c coords " << x_c << ", " << y_c << endl;
	  cerr << "\tlast b coords " << last_x_b << ", " << last_y_b << "\t" 
	       << "  last c coords " << last_x_c << ", " << last_y_c << endl;
	  cerr << "\tx min " << xmin << "\tx max " << xmax << endl;
	  cerr << "\tsign " << sign << "\tslopes " << slope_b << "\t"
	       << slope_c << endl;
	  cerr << "\tintercepts " << intercept_b << "\t" << intercept_c << endl;
	  cerr << "\tintgrl " << x_intgrl << endl;
	}

	if(fabs(x_b-x_c)<three_frame::precision){
	  indexb = (indexb+1)%nvertices;
	  if(indexb==indexc) break;
	  indexc = (indexc-1+nvertices)%nvertices;
	}	else if(x_b<x_c){
	  indexb = (indexb+1)%nvertices;
	} else {
	  indexc = (indexc-1+nvertices)%nvertices;
	}
	if(verbose)
	  cerr << "\tupdated indices " << indexb << "\t"
	       << indexc << "\t" << y_intgrl << endl;
      }


      /////////////////////////////
      //  Perform the y integral //
      /////////////////////////////
      // Find the point with smallest y coordinate
      index = 0;
      y_a = polygon_vertices[0].y(tf);
      for(int i=1; i<nvertices; i++){
	tmp = polygon_vertices[i].y(tf);
	if(tmp<y_a){
	  index = i;
	  y_a = tmp;
	}
      }
      indexb=index;
      indexc=index;

      do {
	indexb = (indexb+1)%nvertices;
	y_b = polygon_vertices[indexb].y(tf);
      } while((y_b-y_a)<three_frame::precision);
    
      do {
	indexc = (indexc-1+nvertices)%nvertices;
	y_c = polygon_vertices[indexc].y(tf);
      } while((y_c-y_a)<three_frame::precision);

      sign = 1;
      if(polygon_vertices[indexc].x(tf)>polygon_vertices[indexb].x(tf))
	sign = -1;

      if(verbose) cerr << "y input indices " << indexb << "\t" << indexc << endl;

      while(1){
	if(indexb==indexc) break;

	x_b = polygon_vertices[indexb].x(tf);
	last_x_b = polygon_vertices[(indexb-1+nvertices)%nvertices].x(tf);
	last_y_b = polygon_vertices[(indexb-1+nvertices)%nvertices].y(tf);

	slope_b = (x_b-last_x_b)/(y_b-last_y_b);
	intercept_b = x_b - slope_b*y_b;

	x_c = polygon_vertices[indexc].x(tf);
	last_x_c = polygon_vertices[(indexc+1)%nvertices].x(tf);
	last_y_c = polygon_vertices[(indexc+1)%nvertices].y(tf);

	slope_c = (x_c-last_x_c)/(y_c-last_y_c);
	intercept_c = x_c - slope_c*y_c;

	ymin = last_y_b > last_y_c ? last_y_b : last_y_c;
	ymax = y_b < y_c ? y_b : y_c;
	xmin = x_b;
	xmax = x_c;
	if(x_b>x_c){
	  xmin = x_c;
	  xmax = x_b;
	}
      
	y_intgrl += sign*((slope_b-slope_c)*(ymax*ymax*ymax - ymin*ymin*ymin)/3.0 +
			  (intercept_b - intercept_c)*(ymax*ymax - ymin*ymin)/2.0);
            
	if(verbose){
	  cerr << "\tindices " << indexb << "\t" << indexc << endl; 
	  cerr << "\tlast indices " << (indexb-1+nvertices)%nvertices
	       << "\t" << (indexc+1)%nvertices << endl;
	  cerr << "\tb coords " << x_b << ", " << y_b << "\t"
	       << " c coords " << x_c << ", " << y_c << endl;
	  cerr << "\tlast b coords " << last_x_b << ", " << last_y_b << "\t" 
	       << "  last c coords " << last_x_c << ", " << last_y_c << endl;
	  cerr << "\ty min " << ymin << "\ty max " << ymax << endl;
	  cerr << "\tsign " << sign << "\tslopes " << slope_b << "\t"
	       << slope_c << endl;
	  cerr << "\tintercepts " << intercept_b << "\t" << intercept_c << endl;
	  cerr << "\tintgrl " << y_intgrl << endl;
	}

	if(fabs(y_b-y_c)<three_frame::precision){
	  indexb = (indexb+1)%nvertices;
	  if(indexb==indexc) break;
	  indexc = (indexc-1+nvertices)%nvertices;
	}	else if(y_b<y_c){
	  indexb = (indexb+1)%nvertices;
	} else {
	  indexc = (indexc-1+nvertices)%nvertices;
	}      
	if(verbose)
	  cerr << "\tupdated indices " << indexb << "\t"
	       << indexc << "\t" << y_intgrl << endl;
      }

      if(verbose)
	cerr << "returning " << x_intgrl << "\t" << y_intgrl << endl << endl;
    } catch(...) {
      cerr << "convex_polygon_integration error\n";
      throw(string("convex_polygon_integration"));
    }
  }
}

#endif
