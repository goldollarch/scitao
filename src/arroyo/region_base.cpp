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

#include <sstream>
#include <cfloat>
#include "region_base.h"
#include "three_frame.h"
#include "computational_geometry.h"

using namespace std;

namespace Arroyo {

  int region_base::verbose_level = 0;

  void rectangular_region::sort_corners(){
    vector<three_point> tmp_corners = corners;
    int max_index;
    double length, max_length = 0;
    for(int i=1; i<4; i++){
      length = (tmp_corners[0]-tmp_corners[i]).length();
      if(length>max_length){
	max_length = length;
	max_index = i;
      }
    }

    corners[2] = tmp_corners[max_index];

    bool one = true;
    for(int i=1; i<4; i++){
      if(i==max_index) continue;
      if(one){
	corners[1] = tmp_corners[i];
	one = false;
      } else corners[3] = tmp_corners[i];
    }
  }

  string rectangular_region::region_status(const rectangular_region & rec_region) const {

    if(!this->aligned(rec_region)){
      cerr << "rectangular_region::region_status error - "
	   << "rectangular region supplied to this function is not aligned "
	   << "with this rectangular region\n";
      throw(string("rectangular_region::region_status"));
    }

    int contains = 0, contained = 0;
    for(int i=0; i<4; i++){
      try{
	if(point_within_polygon(this->corners[i], rec_region.corners)) contained++;
      } catch(...){
	cerr << "rectangular_region::region_status error - "
	     << "could not get status for this three point " << i << endl;
	this->corners[i].print(cerr, "three_point ");
	rec_region.print(cerr, "rec region ");
	throw(string("rectangular_region::region_status"));
      }
      try{
	if(point_within_polygon(rec_region.corners[i], this->corners)) contains++;
      } catch(...){
	cerr << "rectangular_region::region_status error - "
	     << "could not get status for rec_region three point " << i << endl;
	rec_region.corners[i].print(cerr, "three_point ");
	this->print(cerr, "this region ");
	throw(string("rectangular_region::region_status"));
      }
    }

    if(region_base::verbose_level)
      cout << "rectangular_region::region_status - \n\tThis region contains " << contains 
	   << " corners from the arg region.\n\tThe arg region contains " << contained
	   << " corners from this region\n";

    if(contains == 4 && contained == 4) return string("match");
    else if (contains == 0 && contained == 0) return string("no overlap");
    else if (contains == 4) return string("four corner container");
    else if (contained == 4)return string("four corner containee");
    else if (contains == 2 && contained<=2) return string("two corner container");
    else if (contains <= 1 && contained == 2) return string("two corner containee");
    else if (contains == 1 && contained == 1) return string("one corner");
    else {
      cerr << "rectangular_region::region_status error - "
	   << "inconsistent corner container count.\n"
	   << "this region contains " << contains << " corners, and "
	   << contained << " corners are contained by the argument region\n\n";
      this->print(cerr, "this region ");
      rec_region.print(cerr, " arg region ");
      for(int i=0; i<4; i++)
	if(point_within_polygon(this->corners[i], rec_region.corners))
	  cerr << "arg region contains corner " << i << " of this region\n";

      for(int i=0; i<4; i++)
	if(point_within_polygon(rec_region.corners[i], this->corners))
	  cerr << "this region contains corner " << i << " of arg region\n";

      throw(string("rectangular_region::region_status"));
    }
  }  

  rectangular_region::rectangular_region(const rectangular_region & rec_region){
    corners = rec_region.corners;
  }

  rectangular_region::rectangular_region(const vector<three_point> & in_corners){
    if(in_corners.size()!=4){
      cerr << "rectangular_region::rectangular_region error - "
	   << "cannot construct rectangular region with " 
	   << in_corners.size() << " corners\n";
      throw(string("rectangular_region::rectangular_region"));
    }

    // Sort the corners so that the index runs over the corners
    // in a clockwise or counterclockwise direction
    corners = in_corners;
    this->sort_corners();

    // Verify that none of the corners match
    for(int i=0; i<4; i++){
      for(int j=i+1; j<4; j++){
	if(corners[i]==corners[j]){
	  cerr << "rectangular_region::rectangular_region error - "
	       << "two identical corners supplied to this constructor\n";
	  stringstream ss;
	  ss << "corner " << i << "  ";
	  corners[i].print(cerr, ss.str().c_str());
	  ss.str("");
	  ss << "corner " << j << "  ";
	  corners[j].print(cerr, ss.str().c_str());
	  throw(string("rectangular_region::rectangular_region"));
	}
      }
    }

    three_vector tmp_1_unit_vector = (corners[1]-corners[0]);
    tmp_1_unit_vector *= (1/tmp_1_unit_vector.length());
    three_vector tmp_2_unit_vector = (corners[2]-corners[1]);
    tmp_2_unit_vector *= (1/tmp_2_unit_vector.length());
    three_vector tmp_3_unit_vector = (corners[3]-corners[2]);
    tmp_3_unit_vector *= (1/tmp_3_unit_vector.length());
    three_vector orthogonal_vector = cross_product(tmp_1_unit_vector, tmp_2_unit_vector);

    // Verify that these corners are coplanar
    if(fabs(dot_product(orthogonal_vector, tmp_3_unit_vector))>three_frame::precision){
      cerr << "rectangular_region::rectangular_region error - "
	   << "the four three_points supplied to this constructor are not coplanar\n";
      this->print(cerr);
      throw(string("rectangular_region::rectangular_region"));
    }

    // Verify that the edges are orthogonal, indicating that the 
    // corners form a rectangle
    if(fabs(dot_product(tmp_1_unit_vector, tmp_2_unit_vector))>three_frame::precision ||
       fabs(dot_product(tmp_2_unit_vector, tmp_3_unit_vector))>three_frame::precision){
      cerr << "rectangular_region::rectangular_region error - "
	   << "three_points supplied to this constructor do not form a rectangular region\n";
      throw(string("rectangular_region::rectangular_region"));
    }
  }
  
  rectangular_region::rectangular_region(const three_frame & tf, vector<long> axes, double pix){
    if(axes.size()!=2){
      cerr << "rectangular_region::rectangular_region error - "
	   << "cannot construct a rectangular region with " 
	   << axes.size() << " dimensions\n";
      throw(string("rectangular_region::rectangular_region"));
    }

    if(axes[0]<=0 || axes[1]<=0){
      cerr << "rectangular_region::rectangular_region error - "
	   << "axes " << axes[0] << ", " << axes[1]
	   << " are not positive definite\n";
      throw(string("rectangular_region::rectangular_region"));
    }

    if(pix<=0){
      cerr << "rectangular_region::rectangular_region error - "
	   << "cannot construct a rectangular region with pixel scale " 
	   << pix << endl;
      throw(string("rectangular_region::rectangular_region"));
    }

    vector<three_point> corners(4);
    double x = axes[0]*pix/2.0;
    double y = axes[1]*pix/2.0;
    corners[0] = three_point(-x, -y, 0, tf);
    corners[1] = three_point(-x,  y, 0, tf);
    corners[2] = three_point( x,  y, 0, tf);
    corners[3] = three_point( x, -y, 0, tf);

    this->operator=(rectangular_region(corners));
  }

  rectangular_region::rectangular_region(const rectangular_region & rec_region, 
					 const three_vector & nrml,
					 bool along_nrml){

    if(nrml.length()==0){
      cerr << "rectangular_region::rectangular_region error - "
	   << "a null vector has been provided as the normal vector\n";
      throw(string("rectangular_region::rectangular_region"));
    }

    three_vector nrml_hat = nrml*(1/nrml.length());

    vector<three_point> rec_region_corners = rec_region.get_corners();
    three_vector xhat = rec_region_corners[1] - rec_region_corners[0];
    xhat *= (1/xhat.length());
    three_vector yhat = rec_region_corners[3] - rec_region_corners[0]; 
    yhat *= (1/yhat.length());
    three_vector zhat = cross_product(xhat, yhat);
    three_point rec_region_center = rec_region.get_center();

    // if the rectangular region is already in the plane defined
    // by nrml_hat, return a copy
    if(cross_product(zhat, nrml_hat).length() < three_frame::precision){
      this->operator=(rec_region);
      return;
    }
       
    if(fabs(dot_product(zhat, nrml_hat))<three_frame::precision){
      cerr << "rectangular_region::rectangular_region error - "
	   << "cannot project region onto plane orthogonal to the region\n";
      throw(string("rectangular_region::rectangular_region"));
    }

    if(fabs(dot_product(nrml_hat, xhat))>three_frame::precision &&
       fabs(dot_product(nrml_hat, yhat))>three_frame::precision){
      cerr << "rectangular_region::rectangular_region error - "
	   << "projection of this region onto the plane specified by the "
	   << "normal vector provided to this constructor would "
	   << "result in a non-rectangular region\n";
      throw(string("rectangular_region::rectangular_region"));
    }

    /*
    three_vector projection_axis;
    if(along_nrml)
      projection_axis = cross_product(cross_product(zhat, nrml_hat), zhat);
    else
      projection_axis = cross_product(cross_product(zhat, nrml_hat), nrml_hat);
    projection_axis = (1/projection_axis.length())*projection_axis;

    double cos_angle_squared = dot_product(zhat, nrml_hat);
    cos_angle_squared *= cos_angle_squared;

    // The projection vector, when added to the corners of the region,
    // results in a three_point that lies in the plane defined by nrml_hat
    // and is located at the corner of the projected wavefront.  In
    // words, we dot the projection axis into one of the corners and
    // divide by cos_angle_squared to get the length of the component
    // along the projection axis, projected into the plane defined by
    // the normal.  We then multiply this by the dot product of zhat with
    // the projection axis to find the length of the vector that may
    // be added to the corner to get into the plane defined by the
    // normal.  We fabs this to resolve its direction, and make it a
    // vector by multiplying this by zhat
    three_vector projection_vector = 
      fabs((dot_product(rec_region_corners[0]-rec_region_center, projection_axis) / cos_angle_squared) *
	   (dot_product(zhat, projection_axis))) * zhat;
    */

    // New method - find the vector connecting two corners such that 
    // this vector is orthogonal to the cross product of zhat and nrml_hat.
    three_vector half_edge = .5*(rec_region_corners[1] - rec_region_corners[0]);
    if(fabs(dot_product(cross_product(zhat, nrml_hat), half_edge))>three_frame::precision)
      half_edge = .5*(rec_region_corners[2] - rec_region_corners[1]);

    double cos_projection_angle = dot_product(zhat, nrml_hat);
    double sin_projection_angle = sqrt(1-cos_projection_angle*cos_projection_angle);
    three_vector projection_vector;
    if(along_nrml)
      projection_vector = (half_edge.length()*sin_projection_angle)*nrml_hat;
    else
      projection_vector = (half_edge.length()*sin_projection_angle/cos_projection_angle)*zhat;

    for(int i=0; i<4; i++){
      if(dot_product(nrml_hat, rec_region_corners[i]-rec_region_center)>0)
	rec_region_corners[i]-=projection_vector;
      else 
	rec_region_corners[i]+=projection_vector;
    }
    this->operator=(rectangular_region(rec_region_corners));
  }

  rectangular_region::rectangular_region(const rectangular_region & rec_region,
					 double pix){
    
    if(pix<=0){
      cerr << "rectangular_region::rectangular_region error - "
	   << "cannot construct region with dimensions evenly divisible by "
	   << pix << endl;
      throw(string("rectangular_region::rectangular_region"));
    }

    vector<three_point> in_corners = rec_region.get_corners();

    double length = (in_corners[1] - in_corners[0]).length();
    if(fmod(length, pix) > three_frame::precision && 
       fmod(length, pix) < (pix-three_frame::precision)){
      three_vector delta = .5*(in_corners[1]-in_corners[0])*((pix-fmod(length, pix))/length);
      in_corners[0] -= delta;
      in_corners[1] += delta;
      in_corners[2] += delta;
      in_corners[3] -= delta;
    }

    // check...
    length = (in_corners[1] - in_corners[0]).length();
    if(fmod(length, pix) > three_frame::precision && 
       fmod(length, pix) < (pix-three_frame::precision)){
      cerr << "rectangular_region::rectangular_region - error rounding corners\n";
      throw(string("rectangular_region::rectangular_region"));
    }

    length = (in_corners[3] - in_corners[0]).length();
    if(fmod(length, pix) > three_frame::precision && 
       fmod(length, pix) < (pix-three_frame::precision)){
      three_vector delta = .5*(in_corners[3]-in_corners[0])*((pix-fmod(length, pix))/length);
      in_corners[0] -= delta;
      in_corners[1] -= delta;
      in_corners[2] += delta;
      in_corners[3] += delta;
    }

    // check...
    length = (in_corners[3] - in_corners[0]).length();
    if(fmod(length, pix) > three_frame::precision && 
       fmod(length, pix) < (pix-three_frame::precision)){
      cerr << "rectangular_region::rectangular_region - error rounding corners\n";
      throw(string("rectangular_region::rectangular_region"));
    }    
       
    this->operator=(rectangular_region(in_corners));
  }


  rectangular_region & rectangular_region::operator=(const rectangular_region & rec_region){
    if(this==&rec_region) 
      return(*this);
    corners = rec_region.corners;
    return(*this);
  }

  bool rectangular_region::aligned(const rectangular_region & rec_region) const {

    three_vector this_first_edge = this->corners[1] - this->corners[0];
    three_vector this_second_edge = this->corners[2] - this->corners[1];
    three_vector rec_region_first_edge = rec_region.corners[1] - rec_region.corners[0];
    three_vector rec_region_second_edge = rec_region.corners[2] - rec_region.corners[1];
  
    // Make these unit vectors so we don't bias the check
    // when these are vectors with large amplitudes
    this_first_edge *= (1/this_first_edge.length());
    this_second_edge *= (1/this_second_edge.length());
    rec_region_first_edge *= (1/rec_region_first_edge.length());
    rec_region_second_edge *= (1/rec_region_second_edge.length());

    // If the first two edges are parallel, then we check to see if the
    // next two edges are also parallel
    // Likewise, if the first two edges are perpendicular, then we check
    // to see if the next two edges are perpendicular
    if(cross_product(this_first_edge, rec_region_first_edge).length()<three_frame::precision){
      if(cross_product(this_second_edge, rec_region_second_edge).length()<three_frame::precision)
	return true;
      else return false;
    } else if(fabs(dot_product(this_second_edge, rec_region_second_edge))<three_frame::precision){
      if(fabs(dot_product(this_second_edge, rec_region_second_edge))<three_frame::precision)
	return true;
      else return false;
    } else 
      return false;
  }

  bool rectangular_region::contains(const rectangular_region & rec_region) const {
    try{
      string status = this->region_status(rec_region);
      if(status == "match" || status == "four corner container") return true;
      return false;
    } catch(...){
      cerr << "rectangular_region::contains error - "
	   << "could not get region status\n";
      throw(string("rectangular_region::contains"));
    }
  }

  bool rectangular_region::is_contained(const rectangular_region & rec_region) const {
    try{
      string status = this->region_status(rec_region);
      if(status == "match" || status == "four corner containee") return true;
      return false;  
    } catch(...){
      cerr << "rectangular_region::is_contained error - "
	   << "could not get region status\n";
      throw(string("rectangular_region::is_contained"));
    }
  }

  bool rectangular_region::is_disjoint(const rectangular_region & rec_region) const {
    try{
      string status = this->region_status(rec_region);
      if(status == "no overlap") return true;
      return false;   
    } catch(...){
      cerr << "rectangular_region::is_disjoint error - "
	   << "could not get region status\n";
      throw(string("rectangular_region::is_disjoint"));
    }
  }

  vector<three_point> rectangular_region::get_corners() const {
    return(corners);
  }

  three_point rectangular_region::get_center() const {
    return corners[0] + .5*(corners[2]-corners[0]);
  }    

  void rectangular_region::print(ostream & os, const char * prefix, long precision) const {
    for(int i=0; i<4; i++){
      corners[i].print(os, prefix, precision);
      os << endl;
    }
  }

  void rectangular_region::print(ostream & os, const three_frame & tf, const char * prefix, long precision) const {
    for(int i=0; i<4; i++){
      corners[i].print(os, tf, prefix, precision);
      os << endl;
    }
  }

  rectangular_region region_union(const rectangular_region & rec_region1,
				  const rectangular_region & rec_region2){

    // Here we form a three frame in which the regions lie in the
    // x-y plane of the frame.  We then simply find the min and 
    // max x and y values of the three_points in the two regions in 
    // this frame and use them to define the union.  

    three_vector x = rec_region1.corners[1]-rec_region1.corners[0];
    three_vector y = rec_region1.corners[3]-rec_region1.corners[0];
    three_frame tf(rec_region1.get_center(), x, y, cross_product(x,y));
    
    double xmin = rec_region1.corners[0].x(tf);
    double xmax = rec_region1.corners[0].x(tf);
    double ymin = rec_region1.corners[0].y(tf);
    double ymax = rec_region1.corners[0].y(tf);
    
    for(int i=0; i<4; i++){
      if(xmin>rec_region1.corners[i].x(tf)) xmin=rec_region1.corners[i].x(tf);
      if(xmax<rec_region1.corners[i].x(tf)) xmax=rec_region1.corners[i].x(tf);

      if(xmin>rec_region2.corners[i].x(tf)) xmin=rec_region2.corners[i].x(tf);
      if(xmax<rec_region2.corners[i].x(tf)) xmax=rec_region2.corners[i].x(tf);
      
      if(ymin>rec_region1.corners[i].y(tf)) ymin=rec_region1.corners[i].y(tf);
      if(ymax<rec_region1.corners[i].y(tf)) ymax=rec_region1.corners[i].y(tf);

      if(ymin>rec_region2.corners[i].y(tf)) ymin=rec_region2.corners[i].y(tf);
      if(ymax<rec_region2.corners[i].y(tf)) ymax=rec_region2.corners[i].y(tf);
    }

    vector<three_point> union_corners(4);
    union_corners[0] = three_point(xmin, ymin, 0, tf);
    union_corners[1] = three_point(xmax, ymin, 0, tf);
    union_corners[2] = three_point(xmax, ymax, 0, tf);
    union_corners[3] = three_point(xmin, ymax, 0, tf);

    rectangular_region final_rec_region(union_corners);

    if(region_base::verbose_level){
      cout << "region_union - determined union of regions\n";
      rec_region1.print(cout, "region 1     ");
      rec_region2.print(cout, "region 2     ");
      final_rec_region.print(cout, "union region  ");
     }

    return final_rec_region;
  }

  rectangular_region region_intersection(const rectangular_region & rec_region1,
					 const rectangular_region & rec_region2){

    string status;

    try{
      status = rec_region1.region_status(rec_region2);
    } catch(...) {
      cerr << "region_intersection error - could not get status of rectangular regions\n";
      throw(string("region_intersection"));
    }
    if(status == "match"){
      if(region_base::verbose_level) cout << "region_intersection - identical coordinate ranges\n";
      return(rec_region1);
    } else if(status == "four corner containee"){
      if(region_base::verbose_level) 
	cout << "region_intersection - the second region completely contains the first\n"; 
      return(rec_region1);
    } else if(status == "four corner container"){
      if(region_base::verbose_level) 
	cout << "region_intersection - the first region completely contains the second\n"; 
      return(rec_region2);
    } else if(status == "no overlap"){
      cout << "region_intersection error - no overlap between regions\n";
      throw(string("region_intersection"));
    }

    // So much for the easy cases.  Now we do the tough ones.
    
    vector<three_point> intersection_corners;
    
    if(status == "two corner container" || status == "two corner containee"){
  
      rectangular_region container_region, containee_region;

      if(status == "two corner containee"){
	if(region_base::verbose_level)
	  cout << "region_intersection - second region contains 2 corners from the first\n";
	container_region = rec_region2;
	containee_region = rec_region1;
      } else {
	if(region_base::verbose_level)
	  cout << "region_intersection - first region contains 2 corners from the second\n";
	container_region = rec_region1;
	containee_region = rec_region2;
      }

      vector<three_point> exterior_corners;

      for(int i=0; i<4; i++){
	if(point_within_polygon(containee_region.corners[i], container_region.corners))
	  intersection_corners.push_back(containee_region.corners[i]);
	else 
	  exterior_corners.push_back(containee_region.corners[i]);
      }
      if(intersection_corners.size()!=2){
	cerr << "region_intersection error - found " << intersection_corners.size()
	     << " interior corners when expecting two\n";
	throw(string("region_intersection"));
      }

      // This is the vector that forms the edge of the containee
      // rectangle that crosses the boundary of the container
      // region.  The interior and exterior corners are not sorted
      // in the above arrays, so we need to choose the combination
      // with the smallest length.
      three_vector containee_edge_vector = 
	(exterior_corners[0]-intersection_corners[0]).length() >	
	(exterior_corners[0]-intersection_corners[1]).length() ? 
	 exterior_corners[0]-intersection_corners[1] :
	 exterior_corners[0]-intersection_corners[0];

      // Find the two points of intersection that define the rest
      // of the intersecting region
      three_vector unit_containee_edge_vector = containee_edge_vector*(1/containee_edge_vector.length());
      double length, max_length = dot_product(unit_containee_edge_vector, (container_region.corners[0] - intersection_corners[0]));
      for(int i=1; i<4; i++){
	length = dot_product(unit_containee_edge_vector, (container_region.corners[i] - intersection_corners[0]));
	if(length>max_length) max_length = length;
      }

      intersection_corners.push_back(intersection_corners[0] + max_length*unit_containee_edge_vector);
      intersection_corners.push_back(intersection_corners[1] + max_length*unit_containee_edge_vector);

    } else if(status == "one corner") {

      intersection_corners.resize(4);

      int index;
      for(int i=0; i<4; i++){
	if(point_within_polygon(rec_region2.corners[i], rec_region1.corners))
	  intersection_corners[0] = rec_region2.corners[i];
	
	if(point_within_polygon(rec_region1.corners[i], rec_region2.corners))
	  intersection_corners[2] = rec_region1.corners[i];
      }

      double length, min_length = DBL_MAX;
      for(int i=1; i<3; i++){
	length = (rec_region1.corners[i] - rec_region1.corners[0]).length();
	if(length < min_length){
	  index = i;
	  min_length = length;
	}
      }

      three_vector intersection_unit_vector = (rec_region1.corners[index] - rec_region1.corners[0])*(1/min_length);

      three_vector intersection_edge = 
	dot_product(intersection_unit_vector, (intersection_corners[2] - intersection_corners[0]))*intersection_unit_vector;
      intersection_corners[1] = intersection_corners[0] + intersection_edge;
      intersection_corners[3] = intersection_corners[0] + (intersection_corners[2] - intersection_corners[1]);

    }

    rectangular_region final_rec_region(intersection_corners);

    if(region_base::verbose_level){
      cout << "region_intersection - determined intersection of regions\n";
      rec_region1.print(cout, "region 1             ");
      rec_region2.print(cout, "region 2             ");
      final_rec_region.print(cout, "intersection region  ");
    }
  
    if(rec_region1.region_status(final_rec_region)!="four corner container" || 
       rec_region2.region_status(final_rec_region)!="four corner container"){
      cerr << "region_intersection error - "
	   << "in attempting to find the intersection between two sets of " << endl
	   << "corner coords, the resulting corner coords do not contain " << endl
	   << "the first set\n";
      rec_region1.print(cerr, "region 1             ");
      rec_region2.print(cerr, "region 2             ");
      final_rec_region.print(cerr, "intersection region  ");
      throw(string("region_intersection"));
    }
    return rectangular_region(intersection_corners);
  }

  bool operator==(const rectangular_region & rec_region1,
		  const rectangular_region & rec_region2){
    bool match;
    for(int i=0; i<4; i++){
      match = false;
      for(int j=0; j<4; j++){
	if(rec_region1.corners[i]==rec_region2.corners[j]){
	  match = true;
	  break;
	}
      }
      if(match==false) return(false);
    }
    return(true);
  }

  bool operator!=(const rectangular_region & rec_region1,
		  const rectangular_region & rec_region2) {
    return(!operator==(rec_region1, rec_region2));
  }

}
