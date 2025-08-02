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

#ifndef FITS_FACTORY_H
#define FITS_FACTORY_H

#include <map>
#include <string>
#include <vector>
#include <iterator>
#include <iostream>
#include "iofits.h"

namespace Arroyo {

  using std::string;
  using std::vector;
  using std::cerr;
  using std::cout;
  using std::endl;

  class iofits;


  ////////////////////////////////////////
  ///  A class to hold a fits keyval pair
  ///  and an index that holds a relative
  ///  HDU offset

  class fits_keyval_entry {

  private:

    ///////////////////////////
    ///  Private null constructor
    fits_keyval_entry(){};

  protected:

    /// The HDU offset
    int header_offset;

    string key;
    string val;

  public:

    ///////////////////////////
    ///  Construct from the bits
    fits_keyval_entry(const string & k, const string & v, int hdr_offset = 0) {
      if(hdr_offset<0){
	cerr << "fits_keyval_entry::fits_keyval_entry error - "
		  << "could not initialize fits_keyval_entry with header offset "
		  << hdr_offset << endl;
	throw(string("fits_keyval_entry::fits_keyval_entry"));
      }
      header_offset = hdr_offset;
      key = k;
      val = v;
    }

    ///////////////////////////
    ///  Copy constructor
    fits_keyval_entry(const fits_keyval_entry & fkve) {
      this->operator=(fkve);
    }

    ///////////////////////////
    ///  Destrutor
    ~fits_keyval_entry(){};

    ///////////////////////////
    ///  Operator=
    fits_keyval_entry & operator=(const fits_keyval_entry & fkve){
      if(this==&fkve) return(*this);
      header_offset = fkve.header_offset;
      key = fkve.key;
      val = fkve.val;
      return(*this);
    }

    ///////////////////////////
    ///  Returns HDU offset
    int get_offset() const {return(header_offset);};
    
    ///////////////////////////
    ///  Returns key
    string get_key() const {return(key);};

    ///////////////////////////
    ///  Returns val
    string get_val() const {return(val);};

    ///////////////////////////
    ///  Print fits keyval header entry
    void print(std::ostream & os) const {
      os << "\t" << key << "\t" << val
	 << "\tHDU offset " << header_offset << "\t" 
	 << endl;
    }

    ///////////////////////////
    ///  Comparison function to 
    ///  permit sorting
    friend bool operator<(const fits_keyval_entry & kve1, 
			  const fits_keyval_entry & kve2){
      if(kve1.header_offset<kve2.header_offset) return(true);
      return(false);
    }

    ///////////////////////////////////////////
    ///  Friend operator ==
    friend bool operator==(const fits_keyval_entry &f1, const fits_keyval_entry &f2){
      if(f1.header_offset!=f2.header_offset) return(false);
      if(f1.key!=f2.key) return(false);
      if(f1.val!=f2.val) return(false);
      return(true);
    }
  };


  ////////////////////////////////////////
  ///  A class to hold a set of fits keyval pairs

  class fits_keyval_set :
    public vector<fits_keyval_entry> {

  public:

    ///////////////////////////
    ///  Construct from the bits
    fits_keyval_set(){}

    ///////////////////////////
    ///  Copy constructor
    fits_keyval_set(const fits_keyval_set & fkvs) {
      this->operator=(fkvs);
    }

    ///////////////////////////
    ///  Destrutor
    ~fits_keyval_set(){};

    ///////////////////////////
    ///  Operator=
    fits_keyval_set & operator=(const fits_keyval_set & fkvs){
      if(this==&fkvs) return(*this);
      this->vector<fits_keyval_entry>::operator=(fkvs);
      return(*this);
    }

    ///////////////////////////
    ///  Print fits keyval header entry
    void print(std::ostream & os) const {      
      for(uint i=0; i<this->size(); i++)
	(*this)[i].print(os);
    };

    ///////////////////////////
    ///  Friend operator <
    friend bool operator<(const fits_keyval_set & kvs1, 
			  const fits_keyval_set & kvs2){
      return(true);
      /*
      cout << "op<\n";
      if(kvs1.size()<kvs2.size()) return(true);
      for(int i=0; i<kvs1.size(); i++)
	if(kvs1[i]<kvs2[i]) return(true);
      return(false);
      */
    }

    ///////////////////////////////////////////
    ///  Friend operator ==
    friend bool operator==(const fits_keyval_set &s1, const fits_keyval_set &s2){
      s1.print(cout);
      cout << endl;
      s2.print(cout);
      cout << endl << endl;

      if(s1.size()!=s2.size()) return(false);
      for(uint i=0; i<s1.size(); i++)
	if(!(s1[i]==s2[i])) return(false);
      return(true);
    }
  };

  ///
  ///  A template factory class used for loading classes from a fits
  ///  file.  This factory is a singleton.  You instantiate it for the
  ///  particular type of base class that you want to return from this
  ///  factory.  You then register your derived classes using the
  ///  register function, by supplying a map from a
  ///  fits_keyval_header_map to a function pointer.  The
  ///  fits_keyval_header_map is itself a map between pairs of strings
  ///  and integers.  The pairs of strings are interpreted as fits
  ///  keyval pairs, and the int is interpreted as an HDU offset in
  ///  the fits file where the keyval pair is found.  This function
  ///  must take as its argument an iofits reference.  This factory
  ///  then runs through its registry searching the fits file for each
  ///  set of keyval pairs.  If it finds a match, it calls the
  ///  function supplied to the registry.
  ///
  ///  The fits_keyval_header_map between pairs of strings and integers
  ///  is intended to facilitate polymorphic reading and writing of 
  ///  aggregate and multiply inherited objects.  The idea is that for
  ///  one of these compound objects, the first fits header would hold
  ///  information necessary to identify the compound object, and subsequent
  ///  fits HDU's would hold the atomic objects.  Thus the factory must
  ///  know which keyval pairs to look for in each header, and the int
  ///  contains this information.  Since every compound object may itself 
  ///  become incorporated into a larger object through inheritance or aggregation,
  ///  we interpret the integer identifying the fits HDU as a relative offset
  ///  from the current HDU.  In order for this scheme to work, the 
  ///  fits_keyval_header_map must be sorted from smallest header offset 
  ///  (zero offset) to largest header offset.  The Register function ensures that
  ///  this is the case by sorting the vector of fits
  ///
  ///  The first entry in the pair of strings must match the fits header key exactly.
  ///  The second entry in the map must either match the fits header value exactly,
  ///  or must be empty.  If the second entry in the map is empty then
  ///  this class will consider the entry a match if the first entry
  ///  matches the fits header key.
  ///
  ///  This class was modelled after the Object Factory class
  ///  described by Alexandrescu in "Modern C++ design".  While
  ///  this class could almost have been generated from an instantiation
  ///  of the more general class he proposed, this instantiation
  ///  would have relied on his functor, typelist, and singleton 
  ///  classes.  Since no version of his library is available as
  ///  yet, I wrote this simpler version.
  ///

  template<class abstract_product>
  class fits_factory {

  private:

    ///////////////////////////
    ///  Private null constructor
    fits_factory(){};

    ///////////////////////////
    ///  Private destructor
    ~fits_factory(){};

    static fits_factory &instance()
      {
	static fits_factory f;
	return f;
      };
  
    typedef std::map<fits_keyval_set, 
      abstract_product *(*)(const iofits &) > assoc_map;

    /// The map from a vector of fits keyval header entries to a 
    /// pointer to a function that returns a pointer to the abstract
    /// product.
    assoc_map creators_;

  public:

    static bool Register(const fits_keyval_set &fkvs, 
			 abstract_product *(*creator)(const iofits &) ) {
	
      //typedef std::pair<const std::map<const std::string, std::string>, abstract_product *(*)(const char * )> value_type;
      
      /*
	std::cout << "called register creator\n";
	assoc_map::const_iterator e = instance().creators_.end();
        assoc_map::const_iterator i = instance().creators_.find(fits_keyval_map);
	if(i!=e){
	cout << "keyval map already registered\n";
	  throw(string("fits_factory::register"));
	  }
      */
	
      //fits_keyval_entries entries(fkves);
      //std::sort(entries.begin(), entries.end());

      //cout << "REGISTERING\n";
      //fkvs.print(cout);

      if(!(instance().creators_.insert(typename assoc_map::value_type(fkvs, creator))).second){
	cerr << "fits_factory::Register error - could not register class\n";
	fkvs.print(cerr);
	throw(string("fits_factory::Register"));
      }
      /*
      cout << "REGISTERED\n";
      typename assoc_map::const_iterator i = instance().creators_.begin();
      while(i!=instance().creators_.end()){
	i->first.print(cout);
	cout << endl << endl;
	i++;
      } 
      */     
      return true;
    }
    
    static abstract_product *create(const iofits & iof) {
      // iterate through the entries in the assoc_map, 
      // checking for a match
      string val, comment;
      typename assoc_map::const_iterator i = instance().creators_.begin();
      bool match;

      //cout << "checking " << instance().creators_.size() << " entries\n";

      while(i!=instance().creators_.end()){
	match = true;

	//cout << "checking " << i->first.size() << " fits_keyval_entries\n";

	fits_keyval_set::const_iterator j = (i->first).begin();
	while(j!=(i->first).end()){
	  iof.movrel_hdu(j->get_offset());
	  /*
	  cout << "checking key x" << j->get_key() << "x"
	       << " val x" << j->get_val() << "x"
	       << " hdu offset " << j->get_offset()
	       << endl;
	  */
	  if(!iof.key_exists(j->get_key().c_str())){
	    //cout << "\tcould not find key " << j->get_key() << endl;
	    match = false;
	    iof.movrel_hdu(-j->get_offset());
	    break;
	  }
	  // Here if the second entry in the map is empty
	  // we assume that the existence of the first
	  // key in the header is what the user wants to test.
	  // So we skip checking the value
	  if(j->get_val().length()){
	    iof.read_key(j->get_key().c_str(), val, comment);
	    if(j->get_val()!=val){
	      //cout << "\tfound val " << val << " rather than val " << j->get_val() << endl;
	      match = false;
	      iof.movrel_hdu(-j->get_offset());
	      break;
	    }
	  }
	  iof.movrel_hdu(-j->get_offset());
	  j++;
	}
	if(match) return (i->second)(iof);
	i++;
      }
      cout << "fits_factory_create error - "
           << "could not find matching map in fits header - registered entries:\n\n";
      i = instance().creators_.begin();
      while(i!=instance().creators_.end()){
	i->first.print(cout);
	cout << endl << endl;
	i++;
      }
      throw string("fits_factory::create");
    }    
  };

}
#endif 

