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

#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include "computational_geometry.h"

using namespace std;

namespace Arroyo {

  std::vector<Arroyo::three_point> get_line_segment_intersection(const three_point & a1, 
								 const three_point & a2, 
								 const three_point & b1,
								 const three_point & b2){
    // The most common return case
    std::vector<Arroyo::three_point> intersection_points(1);

    // The degenerate case where the first segment is a single point
    if(a1==a2){
      cerr << "get_line_segment_intersection error - first line segment is degenerate\n";
      a1.print(cerr, "first endpoint ");
      a2.print(cerr, "second endpoint ");
      throw(string("get_line_segment_intersection"));
    }

    // The degenerate case where the second segment is a single point
    if(b1==b2){
      cerr << "get_line_segment_intersection error - second line segment is degenerate\n";
      b1.print(cerr, "first endpoint ");
      b2.print(cerr, "second endpoint ");
      throw(string("get_line_segment_intersection"));
    }

    // Check coplanarity.  If not coplanar, the segments do not intersect
    if(fabs(dot_product(cross_product(a2-a1, b2-b1), a1-b1)/(a2-a1).length()/(b2-b1).length()/(a1-b1).length())>three_frame::precision){
      intersection_points.resize(0);
      return(intersection_points);
    }

    double la, lb, lc;
    la = (b1-a1).length();
    lb = (b2-a1).length();
    lc = (a2-b1).length();
    
    double a, b, c, d;
    // Here we normalize by the lengths of the vectors (a2-a1) and
    // (b2-b1).  This normalizes the comparisons so that the fixed
    // value of three_frame::precision doesn't cause numerical
    // overflow for large vectors.
    a = cross_product(b1-a1,a2-a1).length()/(a2-a1).length(); 
    b = cross_product(b2-a1,a2-a1).length()/(a2-a1).length();
    c = cross_product(a1-b1,b2-b1).length()/(b2-b1).length();
    d = cross_product(a2-b1,b2-b1).length()/(b2-b1).length();

    // Here we additionally normalize out the lengths of the
    // vectors like (b1-a1), but only if they are non null.
    if(la>three_frame::precision){
      a/=la;
      c/=la;
    }
    if(lb>three_frame::precision) b/=lb;      
    if(lc>three_frame::precision) b/=lc;
      

    // The degenerate case where the two segments are collinear
    if(a < Arroyo::three_frame::precision && 
       b < Arroyo::three_frame::precision &&
       c < Arroyo::three_frame::precision &&
       d < Arroyo::three_frame::precision){

      if((a1==b1 && a2==b2) || (a1==b2 && a2==b1)){
	intersection_points[0] = a1;
	intersection_points.push_back(a2);
      } else if(a1==b1) {
	intersection_points[0] = a1;
	if(dot_product(a2-a1,b2-b1)>0){
	  if((a2-a1).length()<(b2-b1).length()) intersection_points.push_back(a2);
	  else intersection_points.push_back(b2);
	}
      } else if(a1==b2) {
	intersection_points[0] = a1;
	if(dot_product(a2-a1,b2-b1)<0){
	  if((a2-a1).length()<(b2-b1).length()) intersection_points.push_back(a2);
	  else intersection_points.push_back(b1);
	}
      } else if(a2==b1) {
	intersection_points[0] = a2;
	if(dot_product(a2-a1,b2-b1)<0){
	  if((a2-a1).length()<(b2-b1).length()) intersection_points.push_back(a1);
	  else intersection_points.push_back(b2);
	}
      } else if(a2==b2) {
	intersection_points[0] = a2;
	if(dot_product(a2-a1,b2-b1)>0){
	  if((a2-a1).length()<(b2-b1).length()) intersection_points.push_back(a1);
	  else intersection_points.push_back(b1);
	} 
      } else {

	double e, f, g, h;
	e = dot_product(b1-a1,a2-a1);
	f = dot_product(b2-a1,a2-a1);
	g = dot_product(b1-a2,a2-a1);
	h = dot_product(b2-a2,a2-a1);
	
	if((e<0 && f<0 && g<0 && h<0) || (e>0 && f>0 && g>0 && h>0)){
	  intersection_points.resize(0);
	} else if(e>0 && f>0 && g<0 && h<0){
	  intersection_points[0] = b1;
	  intersection_points.push_back(b2);
	} else if((e<0 && f>0 && g<0 && h>0) ||
		  (e>0 && f<0 && g>0 && h<0)){
	  intersection_points[0] = a1;
	  intersection_points.push_back(a2);
	} else if(e<0 && f>0 && g<0 && h<0){
	  intersection_points[0] = a1;
	  intersection_points.push_back(b2);
	} else if(e>0 && f>0 && g<0 && h>0){
	  intersection_points[0] = b1;
	  intersection_points.push_back(a2);
	} else if(e>0 && f<0 && g<0 && h<0){
	  intersection_points[0] = b1;
	  intersection_points.push_back(a1);
	} else if(e>0 && f>0 && g>0 && h<0){
	  intersection_points[0] = a2;
	  intersection_points.push_back(b2);
	}
      }
    } else {
      if((a<Arroyo::three_frame::precision && c<Arroyo::three_frame::precision) ||
	 (b<Arroyo::three_frame::precision && c<Arroyo::three_frame::precision))
	intersection_points[0] = a1;
      else if((a<Arroyo::three_frame::precision && d<Arroyo::three_frame::precision) ||
	      (b<Arroyo::three_frame::precision && d<Arroyo::three_frame::precision))
	intersection_points[0] = a2;
      else if(a<Arroyo::three_frame::precision) {
	if(dot_product(b1-a1,b1-a2)>0) intersection_points.resize(0);
	else intersection_points[0] = b1;
      } else if(b<Arroyo::three_frame::precision) {
	if(dot_product(b2-a1,b2-a2)>0) intersection_points.resize(0);
	else intersection_points[0] = b2;
      } else if(c<Arroyo::three_frame::precision) {
	if(dot_product(b1-a1,b2-a1)>0) intersection_points.resize(0);
	else intersection_points[0] = a1;
      } else if(d<Arroyo::three_frame::precision) {
	if(dot_product(b1-a2,b2-a2)>0) intersection_points.resize(0);
	else intersection_points[0] = a2;
      } else if(dot_product(cross_product(a1-b1,b2-b1), cross_product(a2-b1,b2-b1)) > 0 ||
		dot_product(cross_product(b1-a1,a2-a1), cross_product(b2-a1,a2-a1)) > 0)
	intersection_points.resize(0);
      else 
	intersection_points[0] = b1 + 
	  (cross_product(b1-a1,a2-a1).length()/cross_product(a2-a1,b2-b1).length())*(b2-b1);
    }
    return(intersection_points);
  }

  three_point get_ray_plane_intersection(const three_point & o_a, const three_vector & n_a,
					 const three_point & o_b, const three_vector & n_b) {

    if(n_a.length_squared()==0 || n_b.length_squared()==0){
      cerr << "get_ray_plane_intersection error - "
	   << "vector supplied to this function has zero length\n";
      throw(string("get_ray_plane_intersection"));
    }

    // Special case - distance along normal is zero
    // because o_a is in the plane defined by o_b and n_b
    three_vector odiff = o_b - o_a;
    if(fabs(dot_product(odiff, n_b))<three_frame::precision){
      return o_a;
    }

    // normalize the vectors n_a and n_b
    three_vector na(n_a);
    na *= 1/n_a.length();
    three_vector nb(n_b);
    nb *= 1/n_b.length();

    // Ensure nb is directed towards oa
    // This cuts down on the number of cases one must consider
    if(dot_product(nb, odiff)<0) nb *=-1;

  
    // Special case - na and nb are perpendicular
    // and the distance along na isn't zero
    if(fabs(dot_product(na, nb))<three_frame::precision){
      cerr << "get_ray_plane_intersection error - "
	   << "vectors supplied to this function are orthogonal\n";
      throw(string("get_ray_plane_intersection"));
    }

    // Special case - na, nb and odiff are all parallel
    if(cross_product(na, nb).length()<three_frame::precision &&
       cross_product(na, odiff).length()<three_frame::precision){
      return(o_b);
    }

    // The three point that results from dropping a perpendicular
    // to the plane from oa
    three_point perpendicular_point = o_a - dot_product(odiff, nb)*nb;

    
    double distance_along_normal = (o_a - perpendicular_point).length()/dot_product(na,nb);

    return(o_a + distance_along_normal*na);

  } 

  three_point get_ray_ray_intersection(const three_point & tp1, 
				       const three_vector & tv1, 
				       const three_point & tp2, 
				       const three_vector & tv2){

    // Ensure these aren't null three vectors
    if(tv1.length()<three_frame::precision || tv2.length()<three_frame::precision){
      cerr << "get_ray_ray_intersection error - a three_vector supplied to this function is null\n";
      tv1.print(cerr, "tv1 ", 15);
      tv2.print(cerr, "tv2 ", 15);
      cerr << "\ttv 1 length " << setprecision(15) << tv1.length() 
	   << "\ttv 2 length " << setprecision(15) << tv2.length() << endl;
      throw(string("get_ray_ray_intersection"));
    }

    // Ensure three vectors aren't parallel
    if(cross_product(tv1,tv2).length()*(1/tv1.length()/tv2.length())<three_frame::precision){
      cerr << "get_ray_ray_intersection error - three_vectors supplied to this function are parallel\n";
      tv1.print(cerr, "tv1 ", 15);
      tv2.print(cerr, "tv2 ", 15);
      cerr << "\tcross product " << setprecision(15) << cross_product(tv1,tv2).length() << endl;
      throw(string("get_ray_ray_intersection"));
    }

    three_vector tv3 = tp1 - tp2;

    // Ensure three points are unique
    if(tv3.length()<three_frame::precision){
      cerr << "get_ray_ray_intersection error - three_points supplied to this function are identical\n";
      tv3.print(cerr, "tv3 ", 15);
      cerr << "\ttv 3 length " << setprecision(15) << tv3.length() << endl;
      throw(string("get_ray_ray_intersection"));
    }

    // Ensure vectors are coplanar
    if(fabs(dot_product(cross_product(tv1, tv2), tv3)/tv1.length()/tv2.length()/tv3.length())>three_frame::precision){
      cerr << "get_ray_ray_intersection error - vectors are not coplanar\n";
      tv1.print(cerr, "tv1 ", 15);
      tv2.print(cerr, "tv2 ", 15);
      tv3.print(cerr, "tv3 ", 15);
      cerr << "\t" << setprecision(15) << dot_product(cross_product(tv1, tv2), tv3) << endl;
      throw(string("get_ray_ray_intersection"));
    }

    // law of sines, using cross products, to get the vector connecting tp1 to the
    // point of intersection.  The direction of the vector is also checked to ensure
    // that we get the right one.
    three_vector mag = (cross_product(tv2,tv3).length()/cross_product(tv1,tv2).length())*tv1;
    if(cross_product(tv2, ((tp1+mag) - tp2)).length()>
       cross_product(tv2, ((tp1-mag) - tp2)).length())
      return(tp1 - mag);
    return(tp1 + mag);
  }

  double get_area_of_polygon(const std::vector<Arroyo::three_point> & polygon_vertices){
    
    int nvertices = polygon_vertices.size();
    if(nvertices<3){
      cerr << "get_area_of_polygon error - number of vertices in polygon "
	   << polygon_vertices.size() << " less than three\n";
      throw(string("get_area_of_polygon"));
    }
 
    int index=1;
    three_vector v1,v2;
    while((v1=(polygon_vertices[index]-polygon_vertices[0])).length()<three_frame::precision)
      index++;
    int index2 = index+1;
    while((v2=(polygon_vertices[index2]-polygon_vertices[index])).length()<three_frame::precision)
      index2++;

    three_vector out_of_plane_normal = cross_product(v1,v2);

    out_of_plane_normal*=(1/out_of_plane_normal.length());

    double triangle_area, polygon_area=0;
    for(int i=2; i<nvertices; i++){
      triangle_area = .5*dot_product(cross_product(polygon_vertices[i-1]-polygon_vertices[0], 
						   polygon_vertices[i]-polygon_vertices[0]),
				     out_of_plane_normal);

      // This check ensures that 3 vertices are not collinear,
      // and that the vertices are consistently sorted either
      // clockwise or counterclockwise.
      //
      // commented out due to normalization issue
      /*
      if(triangle_area<three_frame::precision){
	cerr << "get_area_of_polygon error - "
	     << "area of triangle " 
	     << triangle_area << " is less than or equal to zero.\n";
	stringstream ss;
	for(int i=0; i<polygon_vertices.size(); i++){
	  ss.str("");
	  ss << "vertex " << i+1 << " ";
	  polygon_vertices[i].print(cerr, ss.str().c_str());
	}
	throw(string("get_area_of_polygon"));
      }
      */

      if(!finite(triangle_area)){
	cerr << "get_area_of_polygon error \n";
	(polygon_vertices[i-1]-polygon_vertices[0]).print(cerr, "v1 ");
	(polygon_vertices[i]-polygon_vertices[0]).print(cerr, "v2 ");
	out_of_plane_normal.print(cerr, "oopn");
	throw(string("get_area_of_polygon"));
      }

      polygon_area += triangle_area;
      
    }
    return(polygon_area);
  }

  bool point_within_polygon(const three_point & tp, const vector<three_point> & vertices) {
    bool a, b;
    return(point_within_polygon(tp, vertices, a, b));
  }
  
  bool point_within_polygon(const three_point & tp, 
			    const vector<three_point> & vertices,
			    bool & point_on_edge,
			    bool & point_on_vertex) {
    
    double x, dy;
    int crossings;
    bool xflag0, yflag0, yflag1;
    const three_point *vertex_one, *vertex_two;
    int nvertices = vertices.size();
    three_vector v1, v2;

    if(nvertices < 3){
      cerr << "point_within_polygon error - " << nvertices 
	   << " vertices do not form polygon\n";
      throw(string("point_within_polygon"));
    }

    // check for common points
    for(int i=0; i<nvertices; i++){
      for(int j=i+1; j<nvertices; j++){
	if(vertices[i]==vertices[j]){
	  cerr << "point_within_polygon error - duplicate vertices passed to this function\n";
	  throw(string("point_within_polygon"));
	}
      }
    }

    // Ensure that all vertices are coplanar
    three_vector first = vertices[1] - vertices[0];
    three_vector second = vertices[2] - vertices[0];
    three_vector cross = cross_product(first, second);
    cross *= (1/cross.length());
    three_vector tmp_unit_vector;
    for(int i=3; i<nvertices; i++){
      tmp_unit_vector = vertices[i]-vertices[0];
      tmp_unit_vector *= (1/tmp_unit_vector.length());
      if(fabs(dot_product(cross, tmp_unit_vector))>three_frame::precision){
	cerr << "point_within_polygon error - " 
	     << " vertices are not all coplanar\n";
	for(int j=0; j<nvertices; j++)
	  vertices[j].print(cerr);
	throw(string("point_within_polygon"));
      }
    }

    // If the three_point is equal to one of the vertices, return true
    for(int i=0; i<nvertices; i++)
      if((tp-vertices[i]).length()<three_frame::precision){
	point_on_edge = point_on_vertex = true;
	return true;
      }

    // Ensure that the three_point is coplanar
    tmp_unit_vector = tp-vertices[0];
    tmp_unit_vector *= (1/tmp_unit_vector.length());    
    if(fabs(dot_product(cross, tmp_unit_vector))>three_frame::precision){
      cerr << "point within polygon - out of plane component " 
	   << fabs(dot_product(cross, tmp_unit_vector)) 
	   << " three frame precision " << three_frame::precision << endl;
      cerr << "point_within_polygon error - " 
	   << " point in question doesn't lie in the plane of the polygon\n";
      throw(string("point_within_polygon"));
    }      

    three_frame tf(tp, first, cross_product(cross, first), cross);

    vertex_one = &(vertices[nvertices-1]);
    // get test bit for above/below Y axis 
    yflag0 = ((dy=vertex_one->y(tf) - tp.y(tf)) >= 0.0);
  
    crossings = 0;
    for (int i=0; i<nvertices; i++) {
      // cleverness:  bobble between filling endpoints of edges, so
      // that the previous edge's shared endpoint is maintained.
      if(i%2==1){
	vertex_one = &(vertices[i]);
	yflag0 = ((dy=(vertex_one->y(tf) - tp.y(tf))) >= 0.0);
      } else { 
	vertex_two = &(vertices[i]);
	yflag1 = (vertex_two->y(tf) >= tp.y(tf));
      }
      
      // test to see if the point is on the segment.
      // If this is the case, we return 1 crossing
      v1 = (*vertex_two - *vertex_one);
      v2 = (tp - *vertex_one);
      if(cross_product(v1, v2).length()/v1.length()/v2.length() < three_frame::precision){
	if(dot_product(v1, v2) > 0 && (v1.length() > v2.length())){
	  point_on_vertex = false;
	  point_on_edge = true;
	  return(1);
	}
      }

      // First test - if y coordinates of the two vertices have the 
      // same sign, the ray cannot hit the edge
      if ( yflag0 != yflag1 ) {
	// Second test - if x coordinates are both positive, the ray
	// hits the edge
	if ( (xflag0=(vertex_one->x(tf)>=tp.x(tf))) == (vertex_two->x(tf)>=tp.x(tf)) ) {
	  if ( xflag0 ) crossings++;
	} else {
	  // Final test - if one x coordinate is positive and the other
	  // negative, check to see whether intersection of ray and line
	  // segment has a positive value of x.
	  crossings += ((vertex_one->x(tf) - dy*(vertex_two->x(tf)-vertex_one->x(tf))/(vertex_two->y(tf)-vertex_one->y(tf))) >= tp.x(tf));
	}
      }
    }
    point_on_vertex = point_on_edge = false;
    return(crossings%2) ;
  }
  
  namespace {
    vector<three_point> gcp_return(int first_counter, int second_counter,
				   const vector<three_point> intersection_vertices){
      cout << "\tFinal first counter " << first_counter 
	   << " final second counter " << second_counter << endl;

      stringstream ss;
      for(int i=0; i<intersection_vertices.size(); i++){
	ss.str("");
	ss << "\tintersection vertex " << i << " ";
	intersection_vertices[i].print(cout, ss.str().c_str());
      }
      return(intersection_vertices);
    }

    vector<three_point> reorder(vector<three_point> & intersection_vertices,
				three_vector & sense_of_rotation){
      
      three_point tmp;
      for(int i=1; i<intersection_vertices.size(); i++){
	for(int j=i+1; j<intersection_vertices.size(); j++){
	  if(dot_product(cross_product(intersection_vertices[i]-intersection_vertices[0],
				     intersection_vertices[j]-intersection_vertices[0]),
		       sense_of_rotation)<0){
	    tmp = intersection_vertices[i];
	    intersection_vertices[i] = intersection_vertices[j];
	    intersection_vertices[j] = tmp;
	  }
	}
      }
      return(intersection_vertices);
    }
  }

  vector<three_point> get_convex_polygon_intersection(const vector<three_point> & first_vertices,
						      const vector<three_point> & second_vertices,
						      int verbose){
    
    vector<three_point> intersection_vertices;
    int nfirst_vertices = first_vertices.size();
    int nsecond_vertices = second_vertices.size();
    bool same_orientation = true;
    
    if(verbose>1) cout << endl << "\tget_convex_polygon_intersection\n";

    // Ensure the vectors contain the minimum number of vertices
    if(nfirst_vertices<3 || nsecond_vertices<3){
      cerr << "get_convex_polygon_intersection error - "
	   << "polygons contain too few vertices\n";
      cerr << "First polygon: " << nfirst_vertices << " vertices " << endl;
      cerr << "Second polygon: " << nsecond_vertices << " vertices " << endl;
      throw(string("get_convex_polygon_intersection"));
    }

    // A vector to use in assessing the polygon sense of rotation
    three_vector sense_of_rotation = 
      cross_product(first_vertices[1]-first_vertices[0], first_vertices[2]-first_vertices[0]);

    // Ensure that the polygons are oriented with the same sense of rotation
    if(dot_product(sense_of_rotation,
		   cross_product(second_vertices[1]-second_vertices[0],
				 second_vertices[2]-second_vertices[0]))
       <=0){
      if(sense_of_rotation.length()<three_frame::precision){
	cerr << "get_convex_polygon_intersection error - "
	     << "first three vertices of first polygon are collinear\n";
	first_vertices[0].print(cerr, "first vertex ");
	first_vertices[1].print(cerr, "second vertex ");
	first_vertices[2].print(cerr, "third vertex ");
      }
      else if(cross_product(second_vertices[1]-second_vertices[0],
			    second_vertices[2]-second_vertices[0]).length()
	      <three_frame::precision){
	cerr << "get_convex_polygon_intersection error - "
	     << "first three vertices of second polygon are collinear\n";
	second_vertices[0].print(cerr, "first vertex ");
	second_vertices[1].print(cerr, "second vertex ");
	second_vertices[2].print(cerr, "third vertex ");
      } else {
	cerr << "get_convex_polygon_intersection error - "
	     << "polygons are oriented with the opposite sense of rotation\n";
	stringstream ss;
	for(int i=0; i<first_vertices.size(); i++){
	  ss.str("");
	  ss << "first polygon vertex " << i << " ";
	  first_vertices[i].print(cerr, ss.str().c_str());
	}
	cerr << endl << endl;
	for(int i=0; i<second_vertices.size(); i++){
	  ss.str("");
	  ss << "second polygon vertex " << i << " ";
	  second_vertices[i].print(cerr, ss.str().c_str());
	}
      }
      throw(string("get_convex_polygon_intersection")); 
    }     

    // Normalize to form a unit vector
    sense_of_rotation *= (1/sense_of_rotation.length());
    
    vector<three_point> segment_intersection_points;
    vector<three_point>::const_iterator ci;
    bool vertex_on_vertex, vertex_on_edge;
    bool skip_vertex;
    int first_counter = 0;
    int second_counter = 0;
    bool first_vertex_in_left_halfplane_of_second_segment;
    bool second_vertex_in_left_halfplane_of_first_segment;
    double edge_cross_product;
    while((first_counter<nfirst_vertices || second_counter<nsecond_vertices) &&
	  first_counter < 2*nfirst_vertices && 
	  second_counter<2*nsecond_vertices){

      edge_cross_product = dot_product(cross_product(first_vertices[(first_counter+1)%nfirst_vertices] - 
						     first_vertices[first_counter%nfirst_vertices], 
						     second_vertices[(second_counter+1)%nsecond_vertices] - 
						     second_vertices[second_counter%nsecond_vertices]),
				       sense_of_rotation); 
	
      if(dot_product(cross_product(first_vertices[(first_counter+1)%nfirst_vertices] - 
				   second_vertices[second_counter%nsecond_vertices], 
				   second_vertices[(second_counter+1)%nsecond_vertices] - 
				   second_vertices[second_counter%nsecond_vertices]),
		     sense_of_rotation)<three_frame::precision)	
	first_vertex_in_left_halfplane_of_second_segment = true;
      else
	first_vertex_in_left_halfplane_of_second_segment = false;

      if(dot_product(cross_product(second_vertices[(second_counter+1)%nsecond_vertices] - 
				   first_vertices[first_counter%nfirst_vertices], 
				   first_vertices[(first_counter+1)%nfirst_vertices] - 
				   first_vertices[first_counter%nfirst_vertices]),
		     sense_of_rotation)<three_frame::precision)
	second_vertex_in_left_halfplane_of_first_segment = true;
      else
	second_vertex_in_left_halfplane_of_first_segment = false;

      if(verbose>1) 
	cout << "\tfirst counter " << first_counter 
	     << "\tsecond counter " << second_counter 
	     << "\tedge cross product " << edge_cross_product 
	     << "\tfirsthsecond " << first_vertex_in_left_halfplane_of_second_segment
	     << "\tsecondhfirst" << second_vertex_in_left_halfplane_of_first_segment
	     << endl;
      
      if(fabs(edge_cross_product)>three_frame::precision){

	segment_intersection_points = 
	  get_line_segment_intersection(first_vertices[first_counter%nfirst_vertices],
					first_vertices[(first_counter+1)%nfirst_vertices],
					second_vertices[second_counter%nsecond_vertices],
					second_vertices[(second_counter+1)%nsecond_vertices]);
	
	if(verbose>1) 	
	  cout << "\tsegment intersection size " << segment_intersection_points.size() << endl;
	
	if(segment_intersection_points.size()==1){
	  // We never count the head or tail end of a line segment as an
	  // intersection vertex, as they will be identified by the point
	  // within polygon tests below
	  if(segment_intersection_points[0]!=first_vertices[first_counter%nfirst_vertices] &&
	     segment_intersection_points[0]!=first_vertices[(first_counter+1)%nfirst_vertices] &&
	     segment_intersection_points[0]!=second_vertices[second_counter%nsecond_vertices] &&
	     segment_intersection_points[0]!=second_vertices[(second_counter+1)%nsecond_vertices]){	    

	    // Here we return if we've come full circle
	    if(intersection_vertices.size()>0 && 
	       segment_intersection_points[0]==intersection_vertices[0]){
	      if(verbose){
		cout << "Returning from a\n";
		return(gcp_return(first_counter, second_counter,intersection_vertices));
	      } else return(reorder(intersection_vertices, sense_of_rotation));
	    }
	    
	    // Here we add the segment intersection point to the list of intersection vertices
	    intersection_vertices.push_back(segment_intersection_points[0]);
	    if(verbose>1){
	      intersection_vertices.back().print(cout, "\ta - pushed back vertex ");
	      cout << "\tintersection vertices size " << intersection_vertices.size() << endl;
	    }
	  }
	}
      }

      // The advance rules 
      if(edge_cross_product >= 0){
	if(second_vertex_in_left_halfplane_of_first_segment){
	  if(point_within_polygon(first_vertices[(first_counter+1)%nfirst_vertices], 
				  second_vertices, 
				  vertex_on_vertex,
				  vertex_on_edge)){

	    skip_vertex = false;
	    if(vertex_on_vertex){
	      for(ci=intersection_vertices.begin(); ci!=intersection_vertices.end(); ci++)
		if((*ci)==first_vertices[(first_counter+1)%nfirst_vertices])
		  skip_vertex = true;
	    }

	    if(!skip_vertex){
	      // Here we return if we've come full circle
	      if(intersection_vertices.size()>0 && 
		 intersection_vertices[0]==first_vertices[(first_counter+1)%nfirst_vertices]){
		if(verbose){
		  cout << "Returning from b\n";
		  return(gcp_return(first_counter, second_counter,intersection_vertices));
		} else return(reorder(intersection_vertices, sense_of_rotation));
	      }
	      
	      // Here we add this vertex to the list
	      intersection_vertices.push_back(first_vertices[(first_counter+1)%nfirst_vertices]);
	      if(verbose>1){
		intersection_vertices.back().print(cout, "\tb - pushed back vertex ");
		cout << "\tintersection vertices size " << intersection_vertices.size() << endl;
	      }
	    }
	  }
	  if(verbose>1) cout << "\t\tb - advancing first vertex\n\n";
	  first_counter++;
	} else {
	  if(point_within_polygon(second_vertices[(second_counter+1)%nsecond_vertices], 
				  first_vertices, 
				  vertex_on_vertex,
				  vertex_on_edge)){

	    skip_vertex = false;
	    if(vertex_on_vertex){
	      for(ci=intersection_vertices.begin(); ci!=intersection_vertices.end(); ci++)
		if((*ci)==second_vertices[(second_counter+1)%nsecond_vertices])
		  skip_vertex = true;
	    }

	    if(!skip_vertex){
	      // Here we return if we've come full circle
	      if(intersection_vertices.size()>0 && 
		 intersection_vertices[0]==second_vertices[(second_counter+1)%nsecond_vertices]){
		if(verbose){
		  cout << "Returning from c\n";
		  return(gcp_return(first_counter, second_counter,intersection_vertices));
		} else return(reorder(intersection_vertices, sense_of_rotation));
	      }
	      
	      // Here we add this vertex to the list
	      intersection_vertices.push_back(second_vertices[(second_counter+1)%nsecond_vertices]);
	      if(verbose>1){
		intersection_vertices.back().print(cout, "\tc - pushed back vertex ");
		cout << "\tintersection vertices size " << intersection_vertices.size() << endl;
	      }
	    }
	  }
	  if(verbose>1) cout << "\t\tc - advancing second vertex\n\n";
	  second_counter++;
	}
      } else {
	if(first_vertex_in_left_halfplane_of_second_segment){
	  if(point_within_polygon(second_vertices[(second_counter+1)%nsecond_vertices], 
				  first_vertices, 
				  vertex_on_vertex,
				  vertex_on_edge)){

	    skip_vertex = false;
	    if(vertex_on_vertex){
	      for(ci=intersection_vertices.begin(); ci!=intersection_vertices.end(); ci++)
		if((*ci)==second_vertices[(second_counter+1)%nsecond_vertices])
		  skip_vertex = true;
	    }

	    if(!skip_vertex){
	      // Here we return if we've come full circle
	      if(intersection_vertices.size()>0 && 
		 intersection_vertices[0]==second_vertices[(second_counter+1)%nsecond_vertices]){
		if(verbose){
		  cout << "Returning from c\n";
		  return(gcp_return(first_counter, second_counter,intersection_vertices));
		} else return(reorder(intersection_vertices, sense_of_rotation));
	      }
	      
	      // Here we add this vertex to the list
	      intersection_vertices.push_back(second_vertices[(second_counter+1)%nsecond_vertices]);
	      if(verbose>1){
		intersection_vertices.back().print(cout, "\td - pushed back vertex ");
		cout << "\tintersection vertices size " << intersection_vertices.size() << endl;
	      }
	    }
	  }
	  if(verbose>1) cout << "\t\td - advancing second vertex\n\n";
	  second_counter++;
	} else {
	  if(point_within_polygon(first_vertices[(first_counter+1)%nfirst_vertices], 
				  second_vertices, 
				  vertex_on_vertex,
				  vertex_on_edge)){

	    skip_vertex = false;
	    if(vertex_on_vertex){
	      for(ci=intersection_vertices.begin(); ci!=intersection_vertices.end(); ci++)
		if((*ci)==first_vertices[(first_counter+1)%nfirst_vertices])
		  skip_vertex = true;
	    }

	    if(!skip_vertex){
	      // Here we return if we've come full circle
	      if(intersection_vertices.size()>0 && 
		 intersection_vertices[0]==first_vertices[(first_counter+1)%nfirst_vertices]){
		if(verbose){
		  cout << "Returning from e\n";
		  return(gcp_return(first_counter, second_counter,intersection_vertices));
		} else return(reorder(intersection_vertices, sense_of_rotation));
	      }
	      
	      // Here we add this vertex to the list
	      intersection_vertices.push_back(first_vertices[(first_counter+1)%nfirst_vertices]);
	      if(verbose>1){
		intersection_vertices.back().print(cout, "\te - pushed back vertex ");
		cout << "\tintersection vertices size " << intersection_vertices.size() << endl;
	      }
	    }
	  }
	  if(verbose>1) cout << "\t\te - advancing first vertex\n\n";
	  first_counter++;
	}
      }
    }
    if(verbose){
      cout << "Returning from b\n";
      return(gcp_return(first_counter, second_counter,intersection_vertices));
    } else return(reorder(intersection_vertices, sense_of_rotation));

  } 


  vector<three_point> get_line_segment_circle_intersection(const three_point & a1, 
							   const three_point & a2, 
							   const three_point & circle_origin,
							   const double & radius,
							   int verbose){
    vector<bool> b;
    get_line_segment_circle_intersection(a1, a2, circle_origin, radius, b, verbose);
  }
    
  vector<three_point> get_line_segment_circle_intersection(const three_point & a1, 
							   const three_point & a2, 
							   const three_point & circle_origin,
							   const double & radius,
							   vector<bool> & intersection_point_on_boundary,
							   int verbose){
    
    vector<three_point> intersection_points;
    intersection_point_on_boundary.resize(0);

    // check for invalid arguments
    if(radius<=three_frame::precision){
      cerr << "get_line_segment_circle_intersection error - radius of circle "
	   << radius << " is less than the precision limit\n";
      throw(string("get_line_segment_circle_intersection"));
    }
    
    if((a1-a2).length()<three_frame::precision){
      cerr << "get_line_segment_circle_intersection error - line segment is a single point\n";
      a1.print(cerr, "a1 ");
      a2.print(cerr, "a2 ");
      throw(string("get_line_segment_circle_intersection"));
    }
    
    double a1_origin_distance = (a1 - circle_origin).length();
    double a2_origin_distance = (a2 - circle_origin).length();
    
    /////////////////////////////////////////////////////
    // Case 1 
    // Both points interior to circle
    if((radius - a1_origin_distance)>three_frame::precision &&
       (radius - a2_origin_distance)>three_frame::precision){
      if(verbose) cout << "get_line_segment_circle_intersection: returning from case 1\n";
      return(intersection_points);
    }
    
    /////////////////////////////////////////////////////
    // Case 2
    // Both points lie on the edge of the circle
    if(fabs(radius - a1_origin_distance)<three_frame::precision &&
       fabs(radius - a2_origin_distance)<three_frame::precision){
      intersection_points.push_back(a1);
      intersection_points.push_back(a2);
      intersection_point_on_boundary.push_back(true);
      intersection_point_on_boundary.push_back(true);
      if(verbose) cout << "get_line_segment_circle_intersection: returning from case 2\n";
      return(intersection_points);
    }
    
    double tmp;
    double a1_a2_distance = (a2 - a1).length();
    if(fabs(radius-a1_origin_distance)<three_frame::precision){
      intersection_points.push_back(a1);
      intersection_point_on_boundary.push_back(true);
      if((radius - a2_origin_distance)>three_frame::precision ||
	 dot_product(a2-a1, a1-circle_origin)>-three_frame::precision){
	/////////////////////////////////////////////////////
	// Case 3 
	// The point a1 lies on the edge of the circle, and the point
	// a2 lies within the circle, or the point a2 lies far enough
	// outside the circle so that the segment does not intersect
	// the circle at a second point
	if(verbose) cout << "get_line_segment_circle_intersection: returning from case 3(a)\n";
	return(intersection_points);
      } else {
	/////////////////////////////////////////////////////
	// Case 4 
	// The point a1 lies on the edge of the circle, and the point
	// a2 lies outside the circle but close enough so that the
	// segment intersects the circle at a second point
	tmp = cross_product(circle_origin-a1, a2-a1).length()*(1/a1_a2_distance);
	tmp = sqrt(fabs(radius*radius-tmp*tmp));
	if(!finite(tmp)){
	  cerr << "get_line_segment_circle_intersection internal error in case 4(a)\n";
	  a1.print(cerr, "a1 ", 15);
	  a2.print(cerr, "a2 ", 15);
	  circle_origin.print(cerr, "origin ", 15);
	  cerr << "radius " << setprecision(15) << radius << endl;
	  exit(-1);
	}
	intersection_points.push_back(a1 + (a2-a1)*(2*tmp/a1_a2_distance));
	intersection_point_on_boundary.push_back(false);
	if(verbose) cout << "get_line_segment_circle_intersection: returning from case 4(a)\n";
	return(intersection_points);
      }
    }
       
    if(fabs(radius-a2_origin_distance)<three_frame::precision){
      intersection_points.push_back(a2);
      intersection_point_on_boundary.push_back(true);
      if((radius - a1_origin_distance)>three_frame::precision ||
	 dot_product(a1-a2, a2-circle_origin)>-three_frame::precision){
	/////////////////////////////////////////////////////
	// Case 3 
	// The point a2 lies on the edge of the circle, and the point
	// a1 lies within the circle, or the point a1 lies far enough
	// outside the circle so that the segment does not intersect
	// the circle at a second point
	if(verbose) cout << "get_line_segment_circle_intersection: returning from case 3(b)\n";
	return(intersection_points);
      } else {
	/////////////////////////////////////////////////////
	// Case 4 
	// The point a2 lies on the edge of the circle, and the point
	// a1 lies outside the circle but close enough so that the
	// segment intersects the circle at a second point
	tmp = cross_product(circle_origin-a2, a1-a2).length()*(1/a1_a2_distance);
	tmp = sqrt(fabs(radius*radius-tmp*tmp));
	if(!finite(tmp)){
	  cerr << "get_line_segment_circle_intersection internal error in case 4(b)\n";
	  a1.print(cerr, "a1 ", 15);
	  a2.print(cerr, "a2 ", 15);
	  circle_origin.print(cerr, "origin ", 15);
	  cerr << "radius " << setprecision(15) << radius << endl;
	  exit(-1);
	}
	intersection_points.push_back(a2 + (a1-a2)*(2*tmp/a1_a2_distance));
	intersection_point_on_boundary.push_back(false);
	if(verbose) cout << "get_line_segment_circle_intersection: returning from case 4(b)\n";
	return(intersection_points);
      }
    }
       
    // Both points lie outside the circle
    if(radius-a1_origin_distance<-three_frame::precision &&
       radius-a2_origin_distance<-three_frame::precision){

      // distance of perpendicular from circle origin to line defined
      // by the line segment
      tmp = cross_product(circle_origin-a1, a2-a1).length()*(1/a1_a2_distance);

      /////////////////////////////////////////////////////
      // Case 5 
      // Segment doesn't intersect circle because distance of closest
      // approach is greater than the circle radius
      if((tmp-radius) > three_frame::precision){
	if(verbose) cout << "get_line_segment_circle_intersection: returning from case 5\n";
	return(intersection_points);
      }

      /////////////////////////////////////////////////////
      // Case 6 
      // Segment doesn't intersect circle because both points lie on
      // the same side of the line defining distance of closest
      // approach
      if(dot_product(a2-a1, circle_origin-a1)*dot_product(a2-a1,circle_origin-a2) > 0){
	if(verbose) cout << "get_line_segment_circle_intersection: returning from case 6\n";
	return(intersection_points);
      }
		     
      /////////////////////////////////////////////////////
      // Case 7
      // Segment intersects circle at tangent point
      if(fabs(radius-tmp)<three_frame::precision){
	tmp = sqrt(a1_origin_distance*a1_origin_distance-tmp*tmp);
	intersection_points.push_back(a1 + (a2-a1)*(tmp/a1_a2_distance));
	intersection_point_on_boundary.push_back(false);
	if(verbose) cout << "get_line_segment_circle_intersection: returning from case 7\n";
	return(intersection_points);
      }
      
      /////////////////////////////////////////////////////
      // Case 8
      // Segment intersects circle twice
      tmp = sqrt(a1_origin_distance*a1_origin_distance-tmp*tmp);
      three_point point_of_closest_approach = a1 + (a2-a1)*(tmp/a1_a2_distance);
      
      tmp = (point_of_closest_approach - circle_origin).length();
      tmp = sqrt(radius*radius - tmp*tmp)/(a1 - point_of_closest_approach).length();

      intersection_points.push_back(point_of_closest_approach + (a1 - point_of_closest_approach)*tmp);
      intersection_points.push_back(point_of_closest_approach + (point_of_closest_approach-a1)*tmp);
      intersection_point_on_boundary.push_back(false);
      intersection_point_on_boundary.push_back(false);
      if(verbose) cout << "get_line_segment_circle_intersection: returning from case 8\n";
      return(intersection_points);
    }


    /////////////////////////////////////////////////////
    // Case 9a
    // The point a1 lies within the circle and the point
    // a2 lies outside the circle
    if((radius-a1_origin_distance)>three_frame::precision && 
       (radius-a2_origin_distance)<-three_frame::precision){

      if(a1_origin_distance<three_frame::precision){
	intersection_points.push_back(circle_origin + (a2-circle_origin)*(radius/a2_origin_distance));
	intersection_point_on_boundary.push_back(false);
      } else {
	
	double cos_gamma = dot_product(a2-a1, circle_origin-a1)/a1_a2_distance/a1_origin_distance;
	double sin_gamma = cross_product(circle_origin-a1, a2-a1).length()/a1_a2_distance/a1_origin_distance;
	double sin_eta = a1_origin_distance*sin_gamma/radius;
	double cos_eta = sqrt(1-sin_eta*sin_eta);
	
	tmp = (radius*radius - a1_origin_distance*a1_origin_distance) / (radius*cos_eta - a1_origin_distance*cos_gamma);
	
	intersection_points.push_back(a1 + (a2-a1)*(tmp/a1_a2_distance));
	intersection_point_on_boundary.push_back(false);
      }

      if(verbose) cout << "get_line_segment_circle_intersection: returning from case 9a\n";
      return(intersection_points);
    }

    /////////////////////////////////////////////////////
    // Case 9b
    // The point a2 lies within the circle and the point
    // a1 lies outside the circle
    if((radius-a2_origin_distance)>three_frame::precision && 
       (radius-a1_origin_distance)<-three_frame::precision){

      if(a2_origin_distance<three_frame::precision){
	intersection_points.push_back(circle_origin + (a1-circle_origin)*(radius/a1_origin_distance));
	intersection_point_on_boundary.push_back(false);
      } else {
	
	double cos_gamma = dot_product(a1-a2, circle_origin-a2)/a1_a2_distance/a2_origin_distance;
	double sin_gamma = cross_product(circle_origin-a2, a1-a2).length()/a1_a2_distance/a2_origin_distance;
	double sin_eta = a2_origin_distance*sin_gamma/radius;
	double cos_eta = sqrt(1-sin_eta*sin_eta);
	
	tmp = (radius*radius - a2_origin_distance*a2_origin_distance) / (radius*cos_eta - a2_origin_distance*cos_gamma);
	
	intersection_points.push_back(a2 + (a1-a2)*(tmp/a1_a2_distance));
	intersection_point_on_boundary.push_back(false);
      }

      if(verbose) cout << "get_line_segment_circle_intersection: returning from case 9b\n";
      return(intersection_points);
    }

    cerr << "get_line_segment_circle_intersection error - unexpected case\n";
    throw(string("get_line_segment_circle_intersection"));
  }


  vector<three_point> get_convex_polygon_circle_intersection(const vector<three_point> & polygon_vertices,
							     const three_frame & circle_tf, 
							     double radius){

    try{
      vector<three_point> polygon_intersection_vertices;
      vector<three_point> segment_intersection_vertices;
      vector<bool> intersection_point_on_boundary;

      if(radius<=0){
	cerr << "get_convex_polygon_circle_intersection error - radius " << radius 
	     << " of circle is not positive valued\n";
	throw(string("get_convex_polygon_circle_intersection"));
      }

      int nvertices = polygon_vertices.size();

      // Ensure the vector contains the minimum number of vertices
      if(nvertices<3){
	cerr << "get_convex_polygon_circle_intersection error - "
	     << "polygon contains too few vertices\n";
	throw(string("get_convex_polygon_circle_intersection"));
      }


      bool first_vertex_on_circle, second_vertex_on_circle;
      three_vector orientation_vector = cross_product(polygon_vertices[1]-polygon_vertices[0],
						      polygon_vertices[2]-polygon_vertices[0]);
      orientation_vector *= (1/orientation_vector.length());
    
      for(int i=0; i<nvertices; i++){
      
	if(fabs(dot_product(polygon_vertices[i] - circle_tf, circle_tf.z()))>three_frame::precision){
	  cerr << "get_convex_polygon_circle_intersection error - "
	       << " vertex " << i << " is not in the plane of the three frame\n";
	  polygon_vertices[i].print(cerr, "polygon vertex ");
	  circle_tf.print(cerr, "three frame ");
	  cerr << " dot product " << dot_product(polygon_vertices[i] - circle_tf, circle_tf.z()) << endl;
	  throw(string("get_convex_polygon_circle_intersection"));
	}
    
	if(i>1 && 
	   dot_product(orientation_vector,cross_product(polygon_vertices[i-1]-polygon_vertices[0],
							polygon_vertices[i]-polygon_vertices[0]))
	   <= 0){
	  cerr << "get_convex_polygon_circle_intersection error - "
	       << "vertices are collinear or inconsistently ordered\n";
	  cerr << "\tfailed on vertex " << i 
	       << " with dot product " 
	       << dot_product(orientation_vector,cross_product(polygon_vertices[i-1]-polygon_vertices[0],
							       polygon_vertices[i]-polygon_vertices[0]))
	       << endl;
	  stringstream ss;
	  for(int j=0; j<nvertices; j++){
	    ss.str("");
	    ss << "vertex " << j << " ";
	    polygon_vertices[j].print(cerr, ss.str().c_str());
	  }
	  throw(string("get_convex_polygon_circle_intersection"));
	}

	segment_intersection_vertices = get_line_segment_circle_intersection(polygon_vertices[i],
									     polygon_vertices[(i+1)%nvertices],
									     circle_tf,
									     radius,
									     intersection_point_on_boundary);

	// This check is supposed to add circle/segment intersection points to the
	// list of polygon_intersection_vertices, but only if the intersection
	// point is not on the boundary.  However, there seems to be some cases
	// where the same point gets added in the first and second checks.  The
	// duplicate points don't matter in the function get_area_of_polygon, but
	// this behavior is not really acceptable.  Try to debug this later...
	for(int j=0; j<segment_intersection_vertices.size(); j++){	
	  if(!intersection_point_on_boundary[j]){
	    polygon_intersection_vertices.push_back(segment_intersection_vertices[j]);
	  }
	}
      
	if((radius-(polygon_vertices[(i+1)%nvertices]-circle_tf).length())>-three_frame::precision){
	  polygon_intersection_vertices.push_back(polygon_vertices[(i+1)%nvertices]);
	}
      }
      return(polygon_intersection_vertices);
    } catch(...){
      cerr << "get_convex_polygon_circle_intersection error\n";
      throw(string("get_convex_polygon_circle_intersection"));
    }
  }
}
