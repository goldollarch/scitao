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

#include <cmath>
#include <fstream>
#include <string>
#include <iostream>
#include "colormap.h"

using namespace std;

namespace Arroyo {

  sao_colormap::sao_colormap(vector<vector<double> > & in_R_points, 
			     vector<vector<double> > & in_G_points, 
			     vector<vector<double> > & in_B_points) {

    R_points = in_R_points;
    G_points = in_G_points;
    B_points = in_B_points;
  }

  sao_colormap::sao_colormap(const char * filename) {
    this->read(filename);
  }

  sao_colormap & sao_colormap::operator=(const sao_colormap & scm) {
    if(this==&scm)
      return(*this);
    R_points = scm.R_points;
    G_points = scm.G_points;
    B_points = scm.B_points;
    return(*this);
  }
 
  namespace {

    void parse_sao_colormap_entry_line(ifstream & fs, vector<vector<double> > & points){
      char c;
      fs.get(c);
      double d;
      points.clear();
      while(c=='('){
	points.push_back(vector<double>());
	fs >> d;
	points.back().push_back(d);
	fs.get(c);
	if(c!=','){
	  cerr << "parse_sao_colormap_entry_line error - couldn't parse comma\n";
	  throw(string("parse_sao_colormap_entry_line"));
	}
	fs >> d;
	points.back().push_back(d);
	fs.get(c);
	if(c!=')'){
	  cerr << "parse_sao_colormap_entry_line error - couldn't parse closing parenthesis\n";
	  throw(string("parse_sao_colormap_entry_line"));
	}
	fs.get(c);
      }
      if(points.size()<2){
	cerr << "parse_sao_colormap_entry_line error - " << points.size() 
	     << " entries in the colormap, when there should be at least 2\n";
	throw(string("parse_sao_colormap_entry_line"));
      }
      double level = points[0][0];
      for(int i=1; i<points.size(); i++){
	if(points[i][0]<points[i-1][0]){
	  cerr << "parse_sao_colormap_entry_line error - " 
	       << "levels not in ascending order\n";
	  throw(string("parse_sao_colormap_entry_line"));
	}
      }
    }
  }

  void sao_colormap::read(const char * filename) {

    string fname(filename);

    /*
    string suffix = string(fname.begin()+fname.rfind('.'), fname.end());
    if(suffix!=".sao"){
      cerr << "sao_colormap::read error - "
	   << "cannot load colormap from file with suffix "
	   << suffix << endl
	   << "as this class enforces the .sao suffix convention\n";
      throw(string("sao_colormap::read"));
    }
    */

    ifstream fs(filename);
    if(!fs){
      cerr << "sao_colormap::read error - "
	   << "could not open file " << filename << endl;
      throw(string("sao_colormap::read"));
    }

    int max=100;
    char * space = new char[max];

    do {
      fs.getline(space, max);
    } while(space[0]=='#');
      
    if(space!=string("PSEUDOCOLOR")){
      cerr << "sao_colormap::read error - "
	   << "couldn't parse PSEUDOCOLOR\n";
      throw(string("sao_colormap::read"));
    }

    fs.getline(space, max);
    if(space!=string("RED:")){
      cerr << "sao_colormap::read error - couldn't parse RED:\n";
      throw(string("sao_colormap::read"));
    }
    parse_sao_colormap_entry_line(fs, R_points);

    fs.getline(space, max);
    if(space!=string("GREEN:")){
      cerr << "sao_colormap::read error - couldn't parse GREEN:\n";
      throw(string("sao_colormap::read"));
    }
    parse_sao_colormap_entry_line(fs, G_points);

    fs.getline(space, max);
    if(space!=string("BLUE:")){
      cerr << "sao_colormap::read error - couldn't parse GREEN:\n";
      throw(string("sao_colormap::read"));
    }
    parse_sao_colormap_entry_line(fs, B_points);
    delete [] space;
  }

  void sao_colormap::print(ostream & os, const char * prefix) const {
    cout << "# SAOimage color table" << endl;
    cout << "PSEUDOCOLOR" << endl;
    cout << "RED:" << endl;
    for(int i=0; i<R_points.size(); i++)
      cout << "(" << R_points[i][0] << "," << R_points[i][1] << ")";
    cout << endl;
    cout << "GREEN:" << endl;
    for(int i=0; i<G_points.size(); i++)
      cout << "(" << G_points[i][0] << "," << G_points[i][1] << ")";
    cout << endl;
    cout << "BLUE:" << endl;
    for(int i=0; i<B_points.size(); i++)
      cout << "(" << B_points[i][0] << "," << B_points[i][1] << ")";
    cout << endl;
  }

  void sao_colormap::invert() {
    /*
    for(int i=0; i<R_points.size(); i++)
      R_points[i][1] = 1-R_points[i][1];
    for(int i=0; i<G_points.size(); i++)
      G_points[i][1] = 1-G_points[i][1];
    for(int i=0; i<B_points.size(); i++)
      B_points[i][1] = 1-B_points[i][1];
    */
    vector<vector<double> > tmp;

    tmp = R_points;
    for(int i=0; i<R_points.size(); i++){
      R_points[i] = tmp[R_points.size()-i-1];
      R_points[i][0] = 1-R_points[i][0];
    }
    tmp = G_points;
    for(int i=0; i<G_points.size(); i++){
      G_points[i] = tmp[G_points.size()-i-1];
      G_points[i][0] = 1-G_points[i][0];
    }
    tmp = B_points;
    for(int i=0; i<B_points.size(); i++){
      B_points[i] = tmp[B_points.size()-i-1];
      B_points[i][0] = 1-B_points[i][0];
    }
  }

  char sao_colormap::get_R(double val, double min, double max, bool logscale) const {

    if(min>=max){
      cerr << "sao_colormap::get_R error - min " << min << " and max " << max 
	   << " supplied to this function are inconsistent\n";
      throw(string("sao_colormap::get_R"));
    }

    if(val<min) val = min;
    if(val>max) val = max;

    double frac = (val - min)/(max-min);
    if(logscale)
      frac = .25*log10(frac)+1;

    if(frac<1/(double)255) frac = 1/(double)255;

    int index=0;
    while(R_points[index][0]<frac) index++;

    // perform the linear interpolation
    int intensity;
    return (int)(255*(R_points[index-1][1] + 
		      ((R_points[index][1]-R_points[index-1][1]) *
		       (frac - R_points[index-1][0])/(R_points[index][0] - R_points[index-1][0]))));      
  }

  char sao_colormap::get_G(double val, double min, double max, bool logscale) const {

    if(min>=max){
      cerr << "sao_colormap::get_G error - min " << min << " and max " << max 
	   << " supplied to this function are inconsistent\n";
      throw(string("sao_colormap::get_G"));
    }

    if(val<min) val = min;
    if(val>max) val = max;

    double frac = (val - min)/(max-min);
    if(logscale)
      frac = .25*log10(frac)+1;

    if(frac<1/(double)255) frac = 1/(double)255;

    int index=0;
    while(G_points[index][0]<frac) index++;

    // perform the linear interpolation
    int intensity;
    return (int)(255*(G_points[index-1][1] + 
		      ((G_points[index][1]-G_points[index-1][1]) *
		       (frac - G_points[index-1][0])/(G_points[index][0] - G_points[index-1][0]))));
  }

  char sao_colormap::get_B(double val, double min, double max, bool logscale) const {

    if(min>=max){
      cerr << "sao_colormap::get_B error - min " << min << " and max " << max 
	   << " supplied to this function are inconsistent\n";
      throw(string("sao_colormap::get_B"));
    }

    if(val<min) val = min;
    if(val>max) val = max;

    double frac = (val - min)/(max-min);
    if(logscale)
      frac = .25*log10(frac)+1;

    if(frac<1/(double)255) frac = 1/(double)255;

    int index=0;
    while(B_points[index][0]<frac) index++;

    // perform the linear interpolation
    return (int)(255*(B_points[index-1][1] + 
		      ((B_points[index][1]-B_points[index-1][1]) *
		       (frac - B_points[index-1][0])/(B_points[index][0] - B_points[index-1][0]))));
      
  }

}
  
