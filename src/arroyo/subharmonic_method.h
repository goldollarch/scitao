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

#ifndef SUBHARMONIC_METHOD_H
#define SUBHARMONIC_METHOD_H

#include <iostream>
#include <iomanip>
#include <vector>
#include "AO_sim_base.h"

namespace Arroyo {

  using std::ostream;
  using std::cerr;
  using std::endl;
  using std::string;
  using std::ios;
  using std::setw;
  using std::vector;

  class power_spectrum;

  ///
  /// A virtual base class to represent methods for performing 
  /// subharmonic corrections when constructing refractive_atmospheric_layers
  /// from power spectra
  ///

  class subharmonic_method :
    public AO_sim_base {

  private:

    ///////////////////////////////////////////
    ///  Returns false if the subharmonic method is
    ///  not of the same derived type as this.
    /// 
    ///  If this and subm are of the same derived
    ///  type, returns true if they are equal
    virtual bool equal(const subharmonic_method & subm) const = 0;

  public:

    ///////////////////////////////////////////
    ///  Virtual function to apply the subharmonic
    ///  correction to a raw array of data.
    ///  
    ///  This array is assumed to be halfcomplex
    ///  The axes contain the two dimensions of the
    ///  halfcomplex array, with the second dimension
    ///  corresponding to the shortened axis.
    ///  These axes count complex elements
    template<class T>
    void apply_subharmonic_correction_impl(const power_spectrum & pspec, 
				      const vector<long> & axes, 
				      double pixscale,
				      bool random,
				      T * data) const;
    virtual void apply_subharmonic_correction(const power_spectrum & pspec, 
					      const vector<long> & axes, 
					      double pixscale, 
					      bool random,
					      float * data) const = 0;
    virtual void apply_subharmonic_correction(const power_spectrum & pspec, 
					      const vector<long> & axes, 
					      double pixscale, 
					      bool random,
					      double * data) const = 0;

    ///////////////////////////////////////////
    ///  Null constructor
    subharmonic_method(){};

    ///////////////////////////////////////////
    ///  Virtual destructor
    virtual ~subharmonic_method(){};

    ///////////////////////////////////////////
    ///  Virtual clone method
    ///
    /// Calling routine is responsible for deleting memory
    virtual subharmonic_method * clone() const = 0;

    ///////////////////////////////////////////
    ///  Read from file
    virtual void read(const char * filename) = 0;

    ///////////////////////////////////////////
    ///  Read from iofits
    virtual void read(const iofits & iof) = 0;

    ///////////////////////////////////////////
    ///  Write to file
    virtual void write(const char * filename) const = 0;

    ///////////////////////////////////////////
    ///  Write to iofits
    virtual void write(iofits & iof) const = 0;

    ///////////////////////////////////////////
    ///  Print
    virtual void print(ostream & os, const char * prefix="") const = 0;

    ///////////////////////////////////////////
    ///  Subharmonic methods are intended to be
    ///  used on arrays with either even or odd 
    ///  dimensional axes.  This function returns
    ///  zero if the method should be used on even
    ///  axes, and 1 if the method should be used
    ///  on odd axes.  If either is acceptable, this
    ///  function returns -1
    virtual int intrinsic_dimensionality() const = 0;

    ///////////////////////////////////////////
    ///  Get the number of subharmonics
    virtual long get_subharmonic_depth() const = 0;

    ///////////////////////////////////////////
    ///  Set the number of subharmonics
    virtual void set_subharmonic_depth(long nsubm) = 0;

    ///////////////////////////////////////////
    ///  Factory to construct subharmonic methods from file
    static subharmonic_method * subharmonic_method_factory(const char * filename);

    ///////////////////////////////////////////
    ///  Factory to construct subharmonic methods from an iofits object
    static subharmonic_method * subharmonic_method_factory(const iofits & iof);

    friend bool operator==(const subharmonic_method & subm1, 
			   const subharmonic_method & subm2);


    static int verbose_level;

  };

  bool operator!=(const subharmonic_method & subm1, 
		  const subharmonic_method & subm2);

  ///
  /// A class to indicate that no subharmonic correction will be applied
  ///

  class null_subharmonic_method :
    public subharmonic_method {

  private:

    private:
    
    static const bool factory_registration;    

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("null subharmonic method"));};

    ///////////////////////////////////////////
    ///  Function to apply the subharmonic
    ///  correction to a raw array of data.
    ///  
    ///  This function returns without doing anything
    template<class T>
    void apply_subharmonic_correction_impl(const power_spectrum & pspec, 
				      const vector<long> & axes, 
				      double pixscale,
				      bool random,
				      T * data) const;
    void apply_subharmonic_correction(const power_spectrum & pspec, 
				      const vector<long> & axes, 
				      double pixscale, 
				      bool random,
				      float * data) const {return;};
    void apply_subharmonic_correction(const power_spectrum & pspec, 
				      const vector<long> & axes, 
				      double pixscale, 
				      bool random,
				      double * data) const {return;};

    ///////////////////////////////////////////
    ///  Returns false if the inner scale iscle is
    ///  not of type null_subharmonic_method.  
    ///  Otherwise, returns true.
    bool equal(const subharmonic_method & subm) const {
      const null_subharmonic_method * nsubm = dynamic_cast<const null_subharmonic_method *>(&subm);
      if(nsubm) return(true);
      return(false);    
    };
    
  public:

    ///////////////////////////////////////////
    ///  Null constructor
    null_subharmonic_method(){};

    ///////////////////////////////////////////
    ///  Destructor
    ~null_subharmonic_method(){};

    ///////////////////////////////////////////
    ///  Copy constructor
    null_subharmonic_method(const null_subharmonic_method & nullsubm){
      this->operator=(nullsubm);
    };

    ///////////////////////////////////////////
    ///  Construct from file
    null_subharmonic_method(const char * filename){
      this->read(filename);
    };

    ///////////////////////////////////////////
    ///  Construct from iofits 
    null_subharmonic_method(const iofits & iof){
      this->read(iof);
    };

    ///////////////////////////////////////////
    ///  operator =
    null_subharmonic_method & operator=(const null_subharmonic_method & nullsubm){
      return(*this);
    }

    ///////////////////////////////////////////
    ///  Clone method
    ///
    /// Calling routine is responsible for deleting memory
    null_subharmonic_method * clone() const {
      return new null_subharmonic_method(*this);
    };

    ///////////////////////////////////////////
    ///  Read from file
    void read(const char * filename){
      return;
    }

    ///////////////////////////////////////////
    ///  Read from iofits
    void read(const iofits & iof);

    ///////////////////////////////////////////
    ///  Write to file
    void write(const char * filename) const;

    ///////////////////////////////////////////
    ///  Write to iofits
    void write(iofits & iof) const;

    ///////////////////////////////////////////
    ///  Print
    void print(ostream & os, const char * prefix="") const;

    ///////////////////////////////////////////
    ///  The null subharmonic correction may be 
    ///  used on arrays with arbitrary dimensionality
    int intrinsic_dimensionality() const {return -1;};

    ///////////////////////////////////////////
    ///  Get the number of subharmonics
    long get_subharmonic_depth() const {return 0;};

    ///////////////////////////////////////////
    ///  Set the number of subharmonics
    void set_subharmonic_depth(long nsubm) {
      if(nsubm!=0){
	cerr << "null_subharmonic_method::set_subharmonic_depth error - "
	  << "cannot set subharmonic depth to "
	  << nsubm << " - only zero subharmonics permitted for this class\n";
	throw(string("null_subharmonic_method::set_subharmonic_depth"));
      }
    };
  };

  ///
  /// A class to represent a generalization of the methods 
  /// of subharmonic correction described by Johansson
  /// and Gavel, 1994, and Lane 1992
  ///

  class generalized_subharmonic_method :
    public subharmonic_method {

    private:

    static const bool factory_registration;
    
    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("generalized subharmonic method"));};

    protected:
    
    /// the number of subharmonic levels to correct
    long subharmonic_depth;

    /// the number of subpixels per pixel
    long nsubpixels_per_pixel;

    /// the number of subpixels for each subharmonic
    /// level
    long nsubpixels_per_level;

    ///////////////////////////////////////////
    ///  Function to apply the subharmonic
    ///  correction to a raw array of data.
    ///  
    ///  This array is assumed to be halfcomplex
    template<class T>
    void apply_subharmonic_correction_impl(const power_spectrum & pspec, 
				      const vector<long> & axes, 
				      double pixscale,
				      bool random,
				      T * data) const;
    void apply_subharmonic_correction(const power_spectrum & pspec, 
				      const vector<long> & axes, 
				      double pixscale,
				      bool random,
				      float * data) const;
    void apply_subharmonic_correction(const power_spectrum & pspec, 
				      const vector<long> & axes, 
				      double pixscale,
				      bool random,
				      double * data) const;

    ///////////////////////////////////////////
    ///  Returns false if the subharmonic method
    ///  subm is not of type 
    ///  generalized_subharmonic_method
    /// 
    ///  If subm is of type 
    ///  generalized_subharmonic_method,
    ///  returns true if the instances are equal
    bool equal(const subharmonic_method & subm) const{
      const generalized_subharmonic_method * lsubm = 
	dynamic_cast<const generalized_subharmonic_method *>(&subm);
      if(lsubm && lsubm->get_subharmonic_depth()==this->get_subharmonic_depth()) return(true);
      return(false);
    };

    public:
  
    ///////////////////////////////////////////
    ///  Null constructor
    generalized_subharmonic_method(){
      subharmonic_depth = 0;
      nsubpixels_per_pixel = 3;
      nsubpixels_per_level = 3;
    };

    ///////////////////////////////////////////
    ///  Construct from the bits
    ///
    ///  sbdpth is the number of subharmonic levels 
    ///  to use in the correction
    ///
    ///  nsubpix_per_pix is the number of subpixels 
    ///  in each of the original pixels.  This number
    ///  must be odd
    ///
    ///  nsubpix_per_lvl is the number of subpixels 
    ///  to use in each subharmonic level.  This number
    ///  must be odd and must be greater than 
    ///  nsubpix_per_pix and less than 3*nsubpix_per_pix
    generalized_subharmonic_method(long sbdpth, 
				   long nsubpix_per_lvl,
				   long nsubpix_per_pix);

    ///////////////////////////////////////////
    ///  destructor
    ~generalized_subharmonic_method(){};

    ///////////////////////////////////////////
    ///  Copy constructor
    generalized_subharmonic_method(const generalized_subharmonic_method & 
				   generalized_subm){
      this->operator=(generalized_subm);
    };

    ///////////////////////////////////////////
    ///  Construct from file
    generalized_subharmonic_method(const char * filename){
      this->read(filename);
    };

    ///////////////////////////////////////////
    ///  Construct from iofits 
    generalized_subharmonic_method(const iofits & iof){
      this->read(iof);
    };

    ///////////////////////////////////////////
    ///  operator =
    generalized_subharmonic_method & 
      operator=(const generalized_subharmonic_method & generalized_subm){
      if(this==&generalized_subm)
	return(*this);
      subharmonic_depth = generalized_subm.subharmonic_depth;
      nsubpixels_per_pixel = generalized_subm.nsubpixels_per_pixel;
      nsubpixels_per_level = generalized_subm.nsubpixels_per_level;
      return(*this);
    }

    ///////////////////////////////////////////
    ///  Clone method
    ///
    /// Calling routine is responsible for deleting memory
    generalized_subharmonic_method * clone() const {
      return new generalized_subharmonic_method(*this);
    };

    ///////////////////////////////////////////
    ///  Read from file
    void read(const char * filename);

    ///////////////////////////////////////////
    ///  Read from iofits
    void read(const iofits & iof);

    ///////////////////////////////////////////
    ///  Write to file
    void write(const char * filename) const;

    ///////////////////////////////////////////
    ///  Write to iofits
    void write(iofits & iof) const;

    ///////////////////////////////////////////
    ///  Print
    void print(ostream & os, const char * prefix="") const;

    ///////////////////////////////////////////
    ///  The generalized Johansson Gavel subharmonic correction may be 
    ///  used on arrays with odd dimensionality
    int intrinsic_dimensionality() const {return 1;};

    ///////////////////////////////////////////
    ///  Get the number of subharmonics
    long get_subharmonic_depth() const {return subharmonic_depth;};

    ///////////////////////////////////////////
    ///  Set the number of subharmonics
    void set_subharmonic_depth(long nsubm) {
      if(nsubm<0){
	cerr << "generalized_subharmonic_method::set_subharmonic_depth error - "
	  << "cannot set subharmonic depth to "
	  << nsubm << " - only positive subharmonics permitted for this class\n";
	throw(string("generalized_subharmonic_method::set_subharmonic_depth"));
      }
      subharmonic_depth = nsubm;
    };
 };

  ///
  /// A class to represent the method of subharmonic correction 
  /// described by Lane et al. 1992, waves in random media
  ///

  class Lane_subharmonic_method :
    public generalized_subharmonic_method {

    private:

    static const bool factory_registration;
    
    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("Lane subharmonic method"));};

    public:
  
    ///////////////////////////////////////////
    ///  Null constructor
    Lane_subharmonic_method(){
      subharmonic_depth = 0;
      nsubpixels_per_pixel = 3;
      nsubpixels_per_level = 3;
    };

    ///////////////////////////////////////////
    ///  Construct from the bits
    Lane_subharmonic_method(long sbdpth){
      if(sbdpth<0){
	cerr << "Lane_subharmonic_method::Lane_subharmonic_method error - "
	     << "can't construct with subharmonic_level " << sbdpth 
	     << " less than zero\n";
	throw(string("Lane_subharmonic_method::Lane_subharmonic_method"));
      }
      subharmonic_depth = sbdpth;
      nsubpixels_per_pixel = 3;
      nsubpixels_per_level = 3;
    };

    ///////////////////////////////////////////
    ///  destructor
    ~Lane_subharmonic_method(){};

    ///////////////////////////////////////////
    ///  Copy constructor
    Lane_subharmonic_method(const Lane_subharmonic_method & Lanesubm){
      this->operator=(Lanesubm);
    };

    ///////////////////////////////////////////
    ///  Construct from file
    Lane_subharmonic_method(const char * filename){
      this->read(filename);
    };

    ///////////////////////////////////////////
    ///  Construct from iofits 
    Lane_subharmonic_method(const iofits & iof){
      this->read(iof);
    };

    ///////////////////////////////////////////
    ///  operator =
    Lane_subharmonic_method & operator=(const Lane_subharmonic_method & Lanesubm){
      if(this==&Lanesubm)
	return(*this);
      subharmonic_depth = Lanesubm.subharmonic_depth;
      return(*this);
    }

    ///////////////////////////////////////////
    ///  Clone method
    ///
    /// Calling routine is responsible for deleting memory
    Lane_subharmonic_method * clone() const {
      return new Lane_subharmonic_method(*this);
    };

    ///////////////////////////////////////////
    ///  Read from file
    void read(const char * filename);

    ///////////////////////////////////////////
    ///  Read from iofits
    void read(const iofits & iof);

    ///////////////////////////////////////////
    ///  Write to file
    void write(const char * filename) const;

    ///////////////////////////////////////////
    ///  Write to iofits
    void write(iofits & iof) const;

    ///////////////////////////////////////////
    ///  Print
    void print(ostream & os, const char * prefix="") const;

 };


  ///
  /// A class to represent the method of subharmonic correction 
  /// described by Johansson and Gavel 1994
  ///

  class Johansson_Gavel_subharmonic_method :
    public generalized_subharmonic_method {

    private:

    static const bool factory_registration;
    
    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("Johansson and Gavel subharmonic method"));};

    public:
  
    ///////////////////////////////////////////
    ///  Null constructor
    Johansson_Gavel_subharmonic_method(){
      subharmonic_depth = 0;
      nsubpixels_per_pixel = 3;
      nsubpixels_per_level = 5;
    };

    ///////////////////////////////////////////
    ///  Construct from the bits
    Johansson_Gavel_subharmonic_method(long sbdpth){
      if(sbdpth<0){
	cerr << "Johansson_Gavel_subharmonic_method::Johansson_Gavel_subharmonic_method error - "
	     << "can't construct with subharmonic_level " << sbdpth 
	     << " less than zero\n";
	throw(string("Johansson_Gavel_subharmonic_method::Johansson_Gavel_subharmonic_method"));
      }
      subharmonic_depth = sbdpth;
      nsubpixels_per_pixel = 3;
      nsubpixels_per_level = 5;
    };

    ///////////////////////////////////////////
    ///  destructor
    ~Johansson_Gavel_subharmonic_method(){};

    ///////////////////////////////////////////
    ///  Copy constructor
    Johansson_Gavel_subharmonic_method(const Johansson_Gavel_subharmonic_method & JG_subm){
      this->operator=(JG_subm);
    };

    ///////////////////////////////////////////
    ///  Construct from file
    Johansson_Gavel_subharmonic_method(const char * filename){
      this->read(filename);
    };

    ///////////////////////////////////////////
    ///  Construct from iofits 
    Johansson_Gavel_subharmonic_method(const iofits & iof){
      this->read(iof);
    };

    ///////////////////////////////////////////
    ///  operator =
    Johansson_Gavel_subharmonic_method & operator=(const Johansson_Gavel_subharmonic_method & JG_subm){
      if(this==&JG_subm)
	return(*this);
      subharmonic_depth = JG_subm.subharmonic_depth;
      return(*this);
    }

    ///////////////////////////////////////////
    ///  Clone method
    ///
    /// Calling routine is responsible for deleting memory
    Johansson_Gavel_subharmonic_method * clone() const {
      return new Johansson_Gavel_subharmonic_method(*this);
    };

    ///////////////////////////////////////////
    ///  Read from file
    void read(const char * filename);

    ///////////////////////////////////////////
    ///  Read from iofits
    void read(const iofits & iof);

    ///////////////////////////////////////////
    ///  Write to file
    void write(const char * filename) const;

    ///////////////////////////////////////////
    ///  Write to iofits
    void write(iofits & iof) const;

    ///////////////////////////////////////////
    ///  Print
    void print(ostream & os, const char * prefix="") const;

 };

  ///
  /// A class to represent an experimental method for
  /// subharmonic correction for arrays with even axes
  ///

  class quad_pixel_subharmonic_method :
    public subharmonic_method {

    private:

    static const bool factory_registration;
    
    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("quad pixel subharmonic method"));};

    /// the number of subharmonic levels to correct
    long subharmonic_depth;

    ///////////////////////////////////////////
    ///  Function to apply the subharmonic
    ///  correction to a raw array of data.
    ///  
    ///  This array is assumed to be halfcomplex
    template<class T>
    void apply_subharmonic_correction_impl(const power_spectrum & pspec, 
				      const vector<long> & axes, 
				      double pixscale,
				      bool random,
				      T * data) const;
    void apply_subharmonic_correction(const power_spectrum & pspec, 
				      const vector<long> & axes, 
				      double pixscale,
				      bool random,
				      float * data) const;
    void apply_subharmonic_correction(const power_spectrum & pspec, 
				      const vector<long> & axes, 
				      double pixscale,
				      bool random,
				      double * data) const;

    ///////////////////////////////////////////
    ///  Returns false if the subharmonic method is
    ///  not of type quad_pixel_subharmonic_method
    /// 
    ///  If subm is of type quad_pixel_subharmonic_method,
    ///  returns true if the number of subharmonics are
    ///  equal
    bool equal(const subharmonic_method & subm) const{
      const quad_pixel_subharmonic_method * qpsubm = 
	dynamic_cast<const quad_pixel_subharmonic_method *>(&subm);
      if(qpsubm && qpsubm->get_subharmonic_depth()==this->get_subharmonic_depth()) return(true);
      return(false);
    };

    public:
  
    ///////////////////////////////////////////
    ///  Null constructor
    quad_pixel_subharmonic_method(){subharmonic_depth = 0;};

    ///////////////////////////////////////////
    ///  Construct from the bits
    quad_pixel_subharmonic_method(long sbdpth){
      if(sbdpth<0){
	cerr << "quad_pixel_subharmonic_method::quad_pixel_subharmonic_method error - "
	     << "can't construct with subharmonic_level " << sbdpth 
	     << " less than zero\n";
	throw(string("quad_pixel_subharmonic_method::quad_pixel_subharmonic_method"));
      }
      subharmonic_depth = sbdpth;
    };

    ///////////////////////////////////////////
    ///  destructor
    ~quad_pixel_subharmonic_method(){};

    ///////////////////////////////////////////
    ///  Copy constructor
    quad_pixel_subharmonic_method(const quad_pixel_subharmonic_method & qpsubm){
      this->operator=(qpsubm);
    };

    ///////////////////////////////////////////
    ///  Construct from file
    quad_pixel_subharmonic_method(const char * filename){
      this->read(filename);
    };

    ///////////////////////////////////////////
    ///  Construct from iofits 
    quad_pixel_subharmonic_method(const iofits & iof){
      this->read(iof);
    };

    ///////////////////////////////////////////
    ///  operator =
    quad_pixel_subharmonic_method & operator=(const quad_pixel_subharmonic_method & qpsubm){
      if(this==&qpsubm)
	return(*this);
      subharmonic_depth = qpsubm.subharmonic_depth;
      return(*this);      
    }

    ///////////////////////////////////////////
    ///  Clone method
    ///
    /// Calling routine is responsible for deleting memory
    quad_pixel_subharmonic_method * clone() const {
      return new quad_pixel_subharmonic_method(*this);
    };

    ///////////////////////////////////////////
    ///  Read from file
    void read(const char * filename);

    ///////////////////////////////////////////
    ///  Read from iofits
    void read(const iofits & iof);

    ///////////////////////////////////////////
    ///  Write to file
    void write(const char * filename) const;

    ///////////////////////////////////////////
    ///  Write to iofits
    void write(iofits & iof) const;

    ///////////////////////////////////////////
    ///  Print
    void print(ostream & os, const char * prefix="") const;

    ///////////////////////////////////////////
    ///  The quad pixel subharmonic correction may be 
    ///  used on arrays with even dimensionality
    int intrinsic_dimensionality() const {return 0;};

    ///////////////////////////////////////////
    ///  Get the number of subharmonics
    long get_subharmonic_depth() const {return subharmonic_depth;};

    ///////////////////////////////////////////
    ///  Set the number of subharmonics
    void set_subharmonic_depth(long nsubm) {
      if(nsubm<0){
	cerr << "quad_pixel_subharmonic_method::set_subharmonic_depth error - "
	  << "cannot set subharmonic depth to "
	  << nsubm << " - only positive subharmonics permitted for this class\n";
	throw(string("quad_pixel_subharmonic_method::set_subharmonic_depth"));
      }
      subharmonic_depth = nsubm;
    };
 };

  // To do - Herman & Strugala subharmonic method class - but not a subharmonic method
  // Use Goertzel Reinsch recursor to perform
  // interpolation of coarse sampled grid onto fine
  // sampled grid
}

#endif
