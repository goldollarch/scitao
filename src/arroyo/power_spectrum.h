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

#ifndef POWER_SPECTRUM_H
#define POWER_SPECTRUM_H

#include <cmath>
#include <iostream>
#include "fft_manager.h"
#include "AO_algo.h"
#include "AO_sim_base.h"
#include "fits_header_data.h"
#include "AO_algo.h"
#include "subharmonic_method.h"

namespace Arroyo {

  using std::ostream;
  using std::cerr;
  using std::endl;
  using std::string;
  using std::vector;

  // class forward declarations
  class structure_function;
  
///
/// A virtual base class to represent types of inner scales.
///

  class inner_scale :
    virtual public AO_sim_base {

  public:
  
    ///////////////////////////////////////////
    ///  Null constructor
    inner_scale(){};

    ///////////////////////////////////////////
    ///  Virtual destructor
    virtual ~inner_scale(){};

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
    ///  Return the value at a given spatial frequency
    virtual double value(double spatial_frequency) const = 0;

    ///////////////////////////////////////////
    ///  Factory to construct an inner_scale from a file
    static inner_scale * inner_scale_factory(const char * filename); 

    ///////////////////////////////////////////
    ///  Factory to construct an inner_scale from an iofits object
    static inner_scale * inner_scale_factory(const iofits & iof); 

  };

  ///
  /// A class to represent null inner scales.
  ///

  class null_inner_scale :
    public inner_scale {

    private:
  
    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("null inner scale"));};

    public:

    ///////////////////////////////////////////
    ///  Null constructor
    null_inner_scale(){};

    ///////////////////////////////////////////
    ///  Copy constructor
    null_inner_scale(const null_inner_scale & nullinscle);

    ///////////////////////////////////////////
    ///  Construct from file
    null_inner_scale(const char * filename);

    ///////////////////////////////////////////
    ///  Construct from iofits object
    null_inner_scale(const iofits & iof);

    ///////////////////////////////////////////
    ///  Virtual destructor
    ~null_inner_scale(){};

    ///////////////////////////////////////////
    ///  Operator = 
    null_inner_scale &
      operator=(const null_inner_scale & nullinscle);

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
    ///  Return the value at a given spatial frequency
    double value(double spatial_frequency) const {
      return(1);
    };

    ///////////////////////////////////////////
    ///  operator == - always returns true for null inner scales
    friend bool operator==(const null_inner_scale & nis1, const null_inner_scale & nis2);
  };

  ///////////////////////////////////////////
  ///  operator != - always returns false for null inner scales
  bool operator!=(const null_inner_scale & eis1, const null_inner_scale & eis2);


  ///
  /// A class to represent exponential inner scales.
  ///

  class exponential_inner_scale :
    public inner_scale {

    private:
  
    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("exponential inner scale"));};

    protected:

    /// The value of the inner scale in meters
    double inner_scale_value;
  
    public:

    ///////////////////////////////////////////
    ///  Null constructor
    exponential_inner_scale();

    ///////////////////////////////////////////
    ///  Copy constructor
    exponential_inner_scale(const exponential_inner_scale & expinscle);

    ///////////////////////////////////////////
    ///  Construct from file
    exponential_inner_scale(const char * filename);

    ///////////////////////////////////////////
    ///  Construct from iofits object
    exponential_inner_scale(const iofits & iof);

    ///////////////////////////////////////////
    ///  Construct from a value
    exponential_inner_scale(double inscle);

    ///////////////////////////////////////////
    ///  Virtual destructor
    ~exponential_inner_scale(){};

    ///////////////////////////////////////////
    ///  Operator = 
    exponential_inner_scale &
      operator=(const exponential_inner_scale & expinscle);

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
    ///  Return the value at a given spatial frequency
    double value(double spatial_frequency) const {
      if(spatial_frequency<0) {
	cerr << "exponential_inner_scale::value error - cannot return value for negative spatial frequency " 
		  << spatial_frequency << endl;
	throw(string("exponential_inner_scale::value"));
      }
      // see Sasiela eq 2.25, Hardy eq 3.96 for the factor 5.91
      return(exp(-spatial_frequency*spatial_frequency*inner_scale_value*inner_scale_value/(5.91*5.91)));
    };

    ///////////////////////////////////////////
    ///  operator == - returns true if inner scales are equal
    friend bool operator==(const exponential_inner_scale & eis1, const exponential_inner_scale & eis2);
  };

  ///////////////////////////////////////////
  ///  operator != - returns true if inner scales are not equal
  bool operator!=(const exponential_inner_scale & eis1, const exponential_inner_scale & eis2);

  ///
  /// A class to represent frehlich inner scales.
  ///

  class frehlich_inner_scale : 
    public inner_scale {

    private:

    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("frehlich inner scale"));};

    protected:

    /// The value of the inner scale in meters
    double inner_scale_value;

    public:

    ///////////////////////////////////////////
    ///  Null constructor
    frehlich_inner_scale();

    ///////////////////////////////////////////
    ///  Copy constructor
    frehlich_inner_scale(const frehlich_inner_scale & frehlich_inscle);

    ///////////////////////////////////////////
    ///  Construct from file
    frehlich_inner_scale(const char * filename);

    ///////////////////////////////////////////
    ///  Construct from iofits object
    frehlich_inner_scale(const iofits & iof);

    ///////////////////////////////////////////
    ///  Construct from a value
    frehlich_inner_scale(double inscle);

    ///////////////////////////////////////////
    ///  Virtual destructor
    ~frehlich_inner_scale(){};

    ///////////////////////////////////////////
    ///  Operator = 
    frehlich_inner_scale &
      operator=(const frehlich_inner_scale & frehlich_inscle);

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
    ///  Return the value at a given spatial frequency
    double value(double spatial_frequency) const {
      if(spatial_frequency<0) {
	cerr << "frehlich_inner_scale::value error - cannot return value for negative spatial frequency " 
		  << spatial_frequency << endl;
	throw(string("frehlich_inner_scale::value"));
      }
    
      double delta = 1.109;
      double a[4] = {.70937, 2.8235, -.28086, -.08277};
      double d = exp(-delta * spatial_frequency * inner_scale_value);
      double s = 0;
      for(int i=0; i<4; i++)
	s += a[i]*pow(spatial_frequency*inner_scale_value, i);
    
      return(d*s);
    };
    
    ///////////////////////////////////////////
    ///  operator == - returns true if inner scales are equal
    friend bool operator==(const frehlich_inner_scale & fis1, const frehlich_inner_scale & fis2);

  };

  ///////////////////////////////////////////
  ///  operator != - returns true if inner scales are not equal
  bool operator!=(const frehlich_inner_scale & fis1, const frehlich_inner_scale & fis2);

  ///
  /// A class to represent a power law
  ///

  class power_law :
    virtual public AO_sim_base {

  private:

    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("power law"));};

  protected:
    /// The exponent on the power law
    /// Komolgorov spectra have exponent = 11/3rds
    double exponent;

    /// The coefficient on the power law
    /// Komolgorov spectra have exponent C_n^2
    double coefficient;

  public:

    ///////////////////////////////////////////
    ///  Null constructor
    power_law();

    ///////////////////////////////////////////
    ///  Copy constructor
    power_law(const power_law & plaw);

    ///////////////////////////////////////////
    ///  Construct from file
    power_law(const char * filename);

    ///////////////////////////////////////////
    ///  Construct from iofits object
    power_law(const iofits & iof);

    ///////////////////////////////////////////
    ///  Construct from values
    power_law(double expnt, double coeff);

    ///////////////////////////////////////////
    ///  Construct from values
    power_law(double expnt, 
	      double r_0_meters,
	      double r_0_ref_wavelength_meters);

    ///////////////////////////////////////////
    ///  Virtual destructor
    ~power_law(){};

    ///////////////////////////////////////////
    ///  Copy constructor
    power_law & operator=(const power_law & plaw);

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
    ///  Return the value at a given spatial frequency
    virtual double value(double spatial_frequency) const {
      if(spatial_frequency<0) {
	cerr << "power_law::value error - cannot return value for negative spatial frequency " 
		  << spatial_frequency << endl;
	throw(string("power_law::value"));
      }
      return(coefficient*pow(spatial_frequency, exponent));
    };

    ///////////////////////////////////////////
    ///  Print
    virtual void print(ostream & os, const char * prefix="") const;

    ///////////////////////////////////////////
    ///  Indicates whether there's a pole at zero
    ///  spatial frequency
    virtual bool pole_at_zero_spatial_frequency() const {return true;};

    ///////////////////////////////////////////
    ///  Get the coefficient on the power law
    double get_coefficient() const {return coefficient;};
 
    ///////////////////////////////////////////
    ///  Get the exponent on the power law
    double get_exponent() const {return exponent;};

    ///////////////////////////////////////////
    ///  Factory to construct a power_law from a file
    static power_law * power_law_factory(const char * filename); 

    ///////////////////////////////////////////
    ///  Factory to construct a power_law from an iofits object
    static power_law * power_law_factory(const iofits & iof); 

    ///////////////////////////////////////////
    ///  operator == - returns true if power laws are equal
    friend bool operator==(const power_law & plaw1, 
			   const power_law & plaw2);

  };

  ///////////////////////////////////////////
  ///  operator != - returns true if power laws are not equal
  bool operator!=(const power_law & plaw1, 
		  const power_law & plaw2);

  ///  
  /// A class to represent a power law with a von karman outer scale
  ///

  class von_karman_power_law :
    public power_law {

    private:

    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("von karman power law"));};

    protected:

    /// The value of the outer scale in meters
    double outer_scale_value;

    public:

    ///////////////////////////////////////////
    ///  Null constructor
    von_karman_power_law();

    ///////////////////////////////////////////
    ///  Copy constructor
    von_karman_power_law(const von_karman_power_law & vk_power_law);

    ///////////////////////////////////////////
    ///  Construct from file
    von_karman_power_law(const char * filename);

    ///////////////////////////////////////////
    ///  Construct from iofits object
    von_karman_power_law(const iofits & iof);

    ///////////////////////////////////////////
    ///  Construct from values
    von_karman_power_law(double expnt, 
			 double coeff,
			 double outscle);

    ///////////////////////////////////////////
    ///  Construct from values
    von_karman_power_law(double expnt, 
			 double r_0_meters,
			 double r_0_ref_wavelength_meters,
			 double outscle);

    ///////////////////////////////////////////
    ///  Virtual destructor
    ~von_karman_power_law(){};

    ///////////////////////////////////////////
    ///  Operator = 
    von_karman_power_law &
      operator=(const von_karman_power_law & vk_power_law);
 
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
    ///  Indicates whether there's a pole at zero
    ///  spatial frequency
    bool pole_at_zero_spatial_frequency() const {return false;};

    ///////////////////////////////////////////
    ///  Return the value at a given spatial_frequency
    double value(double spatial_frequency) const {
      if(spatial_frequency<0) {
	cerr << "von_karman_power_law::value error - cannot return value for negative spatial frequency " 
		  << spatial_frequency << endl;
	throw(string("von_karman_power_law::value"));
      }  
      return(coefficient * pow((spatial_frequency*spatial_frequency + 
				4*M_PI*M_PI/(outer_scale_value*outer_scale_value)),
			       exponent/2)); 
    };

    ///////////////////////////////////////////
    ///  Return the outer scale
    double get_outer_scale() const {return outer_scale_value;};
  
    ///////////////////////////////////////////
    ///  operator == - returns true if power laws are equal
    friend bool operator==(const von_karman_power_law & vkplaw1, const von_karman_power_law & vkplaw2);
  };

  ///////////////////////////////////////////
  ///  operator != - returns true if power laws are not equal
  bool operator!=(const von_karman_power_law & vkplaw1, const von_karman_power_law & vkplaw2);

  ///
  /// A class to represent greenwood outer scales.
  ///

  class greenwood_power_law :
    public power_law {

    private:

    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("greenwood power law"));};

    protected:

    /// The value of the outer scale in meters
    double outer_scale_value;

    public:

    ///////////////////////////////////////////
    ///  Null constructor
    greenwood_power_law();

    ///////////////////////////////////////////
    ///  Copy constructor
    greenwood_power_law(const greenwood_power_law & gw_power_law);

    ///////////////////////////////////////////
    ///  Construct from file
    greenwood_power_law(const char * filename);

    ///////////////////////////////////////////
    ///  Construct from iofits object
    greenwood_power_law(const iofits & iof);

    ///////////////////////////////////////////
    ///  Construct from values
    greenwood_power_law(double expnt, 
			double coeff, 
			double outscle);

    ///////////////////////////////////////////
    ///  Construct from values
    greenwood_power_law(double expnt, 
			double r_0_meters,
			double r_0_ref_wavelength_meters,
			double outscle);

    ///////////////////////////////////////////
    ///  Virtual destructor
    ~greenwood_power_law(){};

    ///////////////////////////////////////////
    ///  Operator = 
    greenwood_power_law &
      operator=(const greenwood_power_law & gw_power_law);

    ///////////////////////////////////////////
    ///  Read from file
    void read(const char * filename);

    ///////////////////////////////////////////
    ///  Read from iofits
    void read(const iofits & iof);

    ///////////////////////////////////////////
    ///  Write to iofits
    void write(const char * filename) const;

    ///////////////////////////////////////////
    ///  Write to iofits
    void write(iofits & iof) const;

    ///////////////////////////////////////////
    ///  Print
    void print(ostream & os, const char * prefix="") const;

    ///////////////////////////////////////////
    ///  Indicates whether there's a pole at zero
    ///  spatial frequency
    bool pole_at_zero_spatial_frequency() const {return false;};

    ///////////////////////////////////////////
    ///  Return the value at a given spatial frequency
    double value(double spatial_frequency) const {
      if(spatial_frequency<0) {
	cerr << "greenwood_power_law::value error - cannot return value for negative spatial frequency " 
		  << spatial_frequency << endl;
	throw(string("greenwood_power_law::value"));
      }  
      return(coefficient * pow((spatial_frequency*(spatial_frequency + 2*M_PI/outer_scale_value)),exponent/2)); 
    };

    ///////////////////////////////////////////
    ///  Return the outer scale
    double get_outer_scale() const {return outer_scale_value;};

    ///////////////////////////////////////////
    ///  operator == - returns true if power laws are equal
    friend bool operator==(const greenwood_power_law & gwplaw1, const greenwood_power_law & gwplaw2);
  };

  ///////////////////////////////////////////
  ///  operator != - returns true if power laws are not equal
  bool operator!=(const greenwood_power_law & gwplaw1, const greenwood_power_law & gwplaw2);


  /// class forward declaration
  template<class T>
  class refractive_atmospheric_layer;

///
/// A base class to represent power spectra.
///  

  class power_spectrum :
    virtual public AO_sim_base {

    private:

    ///////////////////////////////////////////
    ///  Virtual member function to create a clone
    virtual power_spectrum * clone() const = 0;

    public:

    ///////////////////////////////////////////
    ///  Null constructor
    power_spectrum(){};

    ///////////////////////////////////////////
    ///  Virtual destructor
    virtual ~power_spectrum(){};

    ///////////////////////////////////////////
    ///  Virtual read from file
    virtual void read(const char * filename) = 0;

    ///////////////////////////////////////////
    ///  Virtual read from iofits
    virtual void read(const iofits & iof) = 0;

    ///////////////////////////////////////////
    ///  Virtual write to file
    virtual void write(const char * filename) const = 0;

    ///////////////////////////////////////////
    ///  Virtual write to iofits
    virtual void write(iofits & iof) const = 0;

    ///////////////////////////////////////////
    ///  Print
    virtual void print(ostream & os, const char * prefix="") const = 0;

    ///////////////////////////////////////////
    ///  Return the value at a given spatial frequency
    virtual double value(double spatial_frequency) const = 0;

    ///////////////////////////////////////////
    ///  Return the power law coefficient
    virtual double get_coefficient() const = 0;

    ///////////////////////////////////////////
    ///  Indicates whether there's a pole at zero
    ///  spatial frequency
    virtual bool pole_at_zero_spatial_frequency() const = 0;

    ///////////////////////////////////////////
    ///  Factory to construct a power spectrum from file
    static power_spectrum * power_spectrum_factory(const char * filename); 

    ///////////////////////////////////////////
    ///  Factory to construct a power spectrum from an iofits object
    static power_spectrum * power_spectrum_factory(const iofits & iof); 

    ///////////////////////////////////////////
    ///  Factory to construct a power spectrum from an iofits object
    static power_spectrum * power_spectrum_factory(const power_spectrum * pspec){
      return pspec->clone(); 
    }; 

    /// verbose level for messages
    static int verbose_level;

  };

  ///
  /// A class to represent an isotropic power law spectrum,
  /// optionally including inner and outer scales.
  ///

  template <class power_law_type, class inner_scale_type>
    class isotropic_power_law_spectrum :
    public power_spectrum {
    
    private:
    
    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("isotropic power law spectrum"));};

    ///////////////////////////////////////////
    ///  Null constructor
    isotropic_power_law_spectrum(){};

    ///////////////////////////////////////////
    ///  Virtual member function to create a clone
    isotropic_power_law_spectrum * clone() const {
      return new isotropic_power_law_spectrum(*this);
    }

    protected:

    power_law_type plaw;
    
    inner_scale_type inscle;

    public:

    ///////////////////////////////////////////
    ///  Copy constructor
    isotropic_power_law_spectrum(const isotropic_power_law_spectrum & ipls);

    ///////////////////////////////////////////
    ///  Construct from iofits object
    isotropic_power_law_spectrum(const iofits & iof);

    ///////////////////////////////////////////
    ///  Construct from a file
    isotropic_power_law_spectrum(const char * filename);

    ///////////////////////////////////////////
    ///  Construct from the pieces
    isotropic_power_law_spectrum(const power_law_type & plaw_type, 
				 const inner_scale_type & inscle_type);

    ///////////////////////////////////////////
    ///  Destructor
    ~isotropic_power_law_spectrum(){};

    ///////////////////////////////////////////
    ///  Operator = 
    isotropic_power_law_spectrum &
      operator=(const isotropic_power_law_spectrum & ipls);

    ///////////////////////////////////////////
    ///  Read from a file
    void read(const char * filename);

    ///////////////////////////////////////////
    ///  Read from an iofits object
    void read(const iofits & iof);

    ///////////////////////////////////////////
    ///  Write to a file
    void write(const char * filename) const;

    ///////////////////////////////////////////
    ///  Write to an iofits object
    void write(iofits & iof) const;

    ///////////////////////////////////////////
    ///  Print
    void print(ostream & os, const char * prefix="") const;

    ///////////////////////////////////////////
    ///  Return the power law coefficient
    double get_coefficient() const {return plaw.get_coefficient();};

    ///////////////////////////////////////////
    ///  Return the value at a given spatial frequency
    double value(double spatial_frequency) const {
      return(inscle.value(spatial_frequency)*plaw.value(spatial_frequency));
    };

    ///////////////////////////////////////////
    ///  Return the power law
    power_law_type get_power_law() const {return plaw;};
    
    ///////////////////////////////////////////
    ///  Return the power law
    inner_scale_type get_inner_scale() const {return inscle;};
    
    ///////////////////////////////////////////
    ///  Indicates whether there's a pole at zero
    ///  spatial frequency
    bool pole_at_zero_spatial_frequency() const {
    	return plaw.pole_at_zero_spatial_frequency();
    }

  };

  ///////////////////////////////////////////
  ///  Check for equality of two isotropic power law spectra
  template<class power_law_type, class inner_scale_type>
    bool operator==(const isotropic_power_law_spectrum<power_law_type, inner_scale_type> & iso1,
		    const isotropic_power_law_spectrum<power_law_type, inner_scale_type> & iso2){
    if(iso1.get_power_law()!=iso2.get_power_law()) return(false);
    if(iso1.get_inner_scale()!=iso2.get_inner_scale()) return(false);
    return(true);
  }

  ///////////////////////////////////////////
  ///  Check for inequality of two isotropic power law spectra
  template<class power_law_type, class inner_scale_type>
    bool operator!=(const isotropic_power_law_spectrum<power_law_type, inner_scale_type> & iso1,
		    const isotropic_power_law_spectrum<power_law_type, inner_scale_type> & iso2){
    return(!operator==(iso1, iso2));
  }

  template<class power_law_type, class inner_scale_type> 
    isotropic_power_law_spectrum<power_law_type, inner_scale_type>::
    isotropic_power_law_spectrum(const isotropic_power_law_spectrum<power_law_type, inner_scale_type> & ipls) {
    this->operator=(ipls);
  }

  template<class power_law_type, class inner_scale_type> 
  isotropic_power_law_spectrum<power_law_type, inner_scale_type>::
  isotropic_power_law_spectrum(const iofits & iof) {
    this->read(iof);
  }

  template<class power_law_type, class inner_scale_type> 
  isotropic_power_law_spectrum<power_law_type, inner_scale_type>::
  isotropic_power_law_spectrum(const char * filename) {
    this->read(filename);
  }

  template<class power_law_type, class inner_scale_type> 
  isotropic_power_law_spectrum<power_law_type, inner_scale_type>::
  isotropic_power_law_spectrum(const power_law_type & plaw, 
			       const inner_scale_type & inscle) {

    // Check the args
    try{dynamic_cast<const power_law &>(plaw);}
    catch(...){
      cerr << "isotropic_power_law_spectrum::isotropic_power_law_spectrum error - argument not of type power law\n";
      throw(string("isotropic_power_law_spectrum::isotropic_power_law_spectrum"));
    }
    try{dynamic_cast<const inner_scale &>(inscle);}
    catch(...){
      cerr << "isotropic_power_law_spectrum::isotropic_power_law_spectrum error - argument not of type inner scale\n";
      throw(string("isotropic_power_law_spectrum::isotropic_power_law_spectrum"));
    }
  
    this->plaw = plaw;
    this->inscle = inscle;
  }

  template<class power_law_type, class inner_scale_type> 
  isotropic_power_law_spectrum<power_law_type, inner_scale_type> & 
    isotropic_power_law_spectrum<power_law_type, inner_scale_type>::
    operator=(const isotropic_power_law_spectrum<power_law_type, inner_scale_type> & ipls){

    if(this==&ipls)
      return(*this);
    plaw = ipls.plaw;
    inscle = ipls.inscle;
    return(*this);
  }

  template<class power_law_type, class inner_scale_type> 
  void isotropic_power_law_spectrum<power_law_type, inner_scale_type>::
    read(const char * filename){
    iofits iof;
    try{iof.open(filename);}
    catch(...){
      cerr << "isotropic_power_law_spectrum::read - "
	   << "error opening file " << filename << endl;
      throw(string("isotropic_power_law_spectrum::read"));
    }
    try{this->read(iof);}
    catch(...){
      cerr << "isotropic_power_law_spectrum::read - "
	   << "error reading "
	   << this->unique_name() << " from file "
	   << filename << endl;
      throw(string("isotropic_power_law_spectrum::read"));
    }
  }
  
  template<class power_law_type, class inner_scale_type> 
  void isotropic_power_law_spectrum<power_law_type, inner_scale_type>::
    read(const iofits & iof){

    if(!iof.key_exists("TYPE")){
      cerr << "isotropic_power_law_spectrum::read error - "
	   << "unrecognized type\n";
      throw(string("isotropic_power_law_spectrum::read"));
    }

    string type, comment;
    iof.read_key("TYPE", type, comment);
    if(type!=this->unique_name()){
      cerr << "isotropic_power_law_spectrum::read error - file of type " 
	   << type << " rather than of type "
	   << this->unique_name() << endl;
      throw(string("isotropic_power_law_spectrum::read"));
    }

    // Move to the power law header
    iof.movrel_hdu(1);
    try{plaw.read(iof);}
    catch(...){
      cerr << "isotropic_power_law_spectrum::read error - "
	   << "could not read power law\n";
      throw(string("isotropic_power_law_spectrum::read"));
    }

    // Move to the inner scale header
    iof.movrel_hdu(1);
    try{inscle.read(iof);}
    catch(...){
      cerr << "isotropic_power_law_spectrum::read error - "
	   << "could not read inner scale\n";
      throw(string("isotropic_power_law_spectrum::read"));
    }

    // leave iof pointing to the next header, if it exists
    if(iof.get_hdu_num()!=iof.get_num_hdus())
      iof.movrel_hdu(1);
  }

  template<class power_law_type, class inner_scale_type> 
  void isotropic_power_law_spectrum<power_law_type, inner_scale_type>::
    write(const char * filename) const {

    iofits iof;
    try{iof.create(filename);}
    catch(...){
      cerr << "isotropic_power_law_spectrum::write - "
	   << "error opening file " << filename << endl;
      throw(string("isotropic_power_law_spectrum::write"));
    }

    try{this->write(iof);}
    catch(...){
      cerr << "isotropic_power_law_spectrum::write - "
	   << "error writing "
	   << this->unique_name() << " to file " 
	   << filename << endl;
      throw(string("isotropic_power_law_spectrum::write"));
    }
  }

  template<class power_law_type, class inner_scale_type>  
    void isotropic_power_law_spectrum<power_law_type, inner_scale_type>::
    write(iofits & iof) const {
    
    fits_header_data<char> tmphdr;
    tmphdr.write(iof);

    string type = this->unique_name();
    string comment = "object type";
    iof.write_key("TYPE", type, comment);

    try{plaw.write(iof);}
    catch(...){
      cerr << "isotropic_power_law_spectrum::write error - "
	   << "could not write power law\n";
      throw(string("isotropic_power_law_spectrum::write"));
    }
    try{inscle.write(iof);}
    catch(...){
      cerr << "isotropic_power_law_spectrum::write error - "
	   << "could not write inner scale\n";
      throw(string("isotropic_power_law_spectrum::write"));
    } 
  }

  template<class power_law_type, class inner_scale_type>  
    void isotropic_power_law_spectrum<power_law_type, inner_scale_type>::
    print(ostream & os, const char * prefix) const {
    int vlspc = 30;
    os.setf(ios::left, ios::adjustfield); 
    os << prefix << "TYPE       = " << setw(vlspc) << this->unique_name()
       << "/" << "object type" << endl;
    plaw.print(os, prefix);
    inscle.print(os, prefix);
  }


  template<class T>
  void initialize_frequency_array(T * data, 
				  const power_spectrum * pspectrum,
				  const vector<long> & axes, 
				  double pixscale, bool random, 
				  const subharmonic_method & subm) {
    
    if(axes.size()!=2){
      cerr << "isotropic_power_law_spectrum::initialize_frequency_array error -\n"
	   << "axes of dimension " << axes.size()
	   << " but must be of dimension 2" << endl;
      throw(string("isotropic_power_law_spectrum::initialize_frequency_array"));
    }

    if(axes[0]<=0 || axes[1]<=0){
      cerr << "isotropic_power_law_spectrum::initialize_frequency_array error -\n"
	   << "invalid axes " << axes[0] << "\t" << axes[1] << endl;
      throw(string("isotropic_power_law_spectrum::initialize_frequency_array"));
    }

    if(pixscale<=0){
      cerr << "isotropic_power_law_spectrum::initialize_frequency_array error -\n"
	   << "invalid pixel scale " << pixscale << endl;
      throw(string("isotropic_power_law_spectrum::initialize_frequency_array"));
    }

    // the coordinate locations of pixels and subpixels
    double pixel_location[2];

    double x_halfpix=0, y_halfpix=0;
    int x_extrapix=1, y_extrapix=1;
    if(axes[1]%2==0){
      x_halfpix = .5;
      x_extrapix = 0;
    }
    if(axes[0]%2==0){
      y_halfpix = .5;
      y_extrapix = 0;
    }
  
    int index;
    double xfreq_interval = 2*M_PI/(pixscale*axes[1]);
    double yfreq_interval = 2*M_PI/(pixscale*axes[0]);
    double val, r1, r2;

    // The factor to fix the units on the delta function
    // contribution to make it look like a smooth probability
    // distribution.
    // See Johansson & Gavel
    double norm = sqrt(xfreq_interval*yfreq_interval);

    // Here if the zero frequency value of the power spectrum
    // is infinite, we will not fill it in
    bool skip_zero_frequency = false;
    if(pspectrum->pole_at_zero_spatial_frequency()) skip_zero_frequency = true;

    for(int i=-axes[1]/2; i<axes[1]/2+x_extrapix; i++){
      for(int j=-axes[0]/2; j<axes[0]/2+y_extrapix; j++){

	pixel_location[0] = (i+x_halfpix)*xfreq_interval;
	pixel_location[1] = (j+y_halfpix)*yfreq_interval;
      
 	// Here, if both axes are odd we have to skip the middle pixel,
	// since its nominal value is infinite.
	if(pixel_location[0]==0 && pixel_location[1]==0 && skip_zero_frequency){
	  index = (axes[1]/2)*axes[0] + axes[0]/2;
	  data[2*index] = data[2*index+1] = 0;
	  continue;
	} 

	val = norm*sqrt(pspectrum->value(sqrt(pixel_location[0]*pixel_location[0]+
					 pixel_location[1]*pixel_location[1])));

	index = (i+axes[1]/2)*axes[0] + j + axes[0]/2;
	if(random){
	  box_mueller(r1,r2);
	  data[2*index] = r1*val;
	  data[2*index+1] = r2*val;
	} else {
	  data[2*index] = val*val;
	  data[2*index+1] = 0;
	} 
      }
    }

    if(random)
      subm.apply_subharmonic_correction(*pspectrum, axes, pixscale, true, data);
    else {							  
      subm.apply_subharmonic_correction(*pspectrum, axes, pixscale, false, data);
    }	

  }

  // Here we compensate for the halfpixel shift that arises if
  // the axes are even.  This shift was discussed in the comments
  // to the function diffractive_wavefront<T>::far_field_propagator
  template<class T>
    void halfpixel_fft(const vector<long> & axes, 
		       T * data,
		       fft_manager<T> & fft_mgr){

    double x_halfpix=0, y_halfpix=0;
    int x_extrapix=1, y_extrapix=1;
    if(axes[1]%2==0){
      x_halfpix = .5;
      x_extrapix = 0;
    }
    if(axes[0]%2==0){
      y_halfpix = .5;
      y_extrapix = 0;
    }
    int index;
    double xslope = M_PI/(double)axes[1];
    double yslope = M_PI/(double)axes[0];
    double cp, sp, tmp;
    if(axes[1]%2==1) xslope = 0;
    if(axes[0]%2==1) yslope = 0;

    if(xslope!=0 || yslope!=0){
      for(int i=-axes[1]/2; i<axes[1]/2+x_extrapix; i++){
	for(int j=-axes[0]/2; j<axes[0]/2+y_extrapix; j++){
	  cp = cos(xslope*(i+x_halfpix) + yslope*(j+y_halfpix));
	  sp = sin(xslope*(i+x_halfpix) + yslope*(j+y_halfpix));
	  index = (i+axes[1]/2)*axes[0]+j+axes[0]/2;
	  tmp = cp*data[2*index] + sp*data[2*index+1];
	  data[2*index+1] = -sp*data[2*index] + cp*data[2*index+1];
	  data[2*index] = tmp;
	}
      }
    }

    complex_cyclic_permutation(axes, axes[1]/2, axes[0]/2, data);
      
    vector<long> flipped_axes(2, axes[1]);
    flipped_axes[1]=axes[0];
      
    fft_mgr.forward_fft(flipped_axes, true, true, data);
      
    complex_cyclic_permutation(axes, -axes[1]/2, -axes[0]/2, data);
      
    if(xslope!=0 || yslope!=0){
      for(int i=-axes[1]/2; i<axes[1]/2+x_extrapix; i++){
	for(int j=-axes[0]/2; j<axes[0]/2+y_extrapix; j++){
	  cp = cos(xslope*(i+x_halfpix) + yslope*(j+y_halfpix));
	  sp = sin(xslope*(i+x_halfpix) + yslope*(j+y_halfpix));
	  index = (i+axes[1]/2)*axes[0]+j+axes[0]/2;
	  tmp = cp*data[2*index] + sp*data[2*index+1];
	  data[2*index+1] = -sp*data[2*index] + cp*data[2*index+1];
	  data[2*index] = tmp;
	}
      }
    }
  }


  /// COMMENTED OUT FOR THE FOLLOWING REASON:
  /// The thing that looks like it would be difficult to resolve:
  /// how to represent a generic power spectrum that may
  /// have an infinity at zero spatial frequency?  It may be 
  /// that this is actually unnecessary - as pictured in Lane et. al,
  /// they recover a komolgorov power spectrum from simulation 
  /// in which the value at zero frequency is finite.
  /// We can leave this for the moment...

  /*
    CLASS 
    generic_power_spectrum
  
    class to represent generic power spectra

    class generic_power_spectrum {

    protected:

    /// the number of points
    long npts;

    /// the data
    float * pspec_values;

    /// the corresponding spatial frequencies
    float * pspec_freqs;

    public:

    ///////////////////////////////////////////
    ///  null constructor
    generic_power_spectrum();

    ///////////////////////////////////////////
    ///  copy constructor
    generic_power_spectrum(const generic_power_spectrum & gps);

    ///////////////////////////////////////////
    ///  destructor
    ~generic_power_spectrum();

    ///////////////////////////////////////////
    ///  operator = 
    generic_power_spectrum &
    operator=(const generic_power_spectrum & gps);

    ///////////////////////////////////////////
    ///  read from file
    void read(const char * filename);

    ///////////////////////////////////////////
    ///  read from iofits
    void read(const iofits & iof);

    ///////////////////////////////////////////
    ///  write to file
    void write(const char * filename);

    ///////////////////////////////////////////
    ///  write to iofits
    void write(iofits & iof);

    ///////////////////////////////////////////
    ///  function to get the structure function
    one_d_structure_function 
    get_one_d_structure_function(long npts, double spacing) const;

    ///////////////////////////////////////////
    ///  function to get the structure function
    two_d_structure_function 
    get_two_d_structure_function(vector<long> & npts, 
    vector<double> & spacing) const;

    ///////////////////////////////////////////
    ///  function to get a
    ///  refractive_atmospheric_layer
    refractive_atmospheric_layer 
    get_random_layer(double pixel_scale,
    vector<long> pixel_dimensions) const;
  

    }
  */

}

#endif
