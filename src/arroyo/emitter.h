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

#ifndef EMITTER_H
#define EMITTER_H

#include "three_vector.h"
#include "diffractive_wavefront.h"
#include "AO_sim_base.h"

namespace Arroyo {
  class iofits;
  class three_point;
}

namespace Arroyo {


  ///
  /// A virtual base class to represent an emitter of wavefronts.
  ///

  class emitter :
    virtual public AO_sim_base {

  public:
    
    ///////////////////////////////////////////
    ///  Null constructor
    emitter(){};

    ///////////////////////////////////////////
    ///  Destructor
    virtual ~emitter(){};

    ///////////////////////////////////////////
    ///  Read from Arroyo::iofits object
    virtual void read(const char * filename) = 0;

    ///////////////////////////////////////////
    ///  Read from iofits object
    virtual void read(const Arroyo::iofits & iof) = 0;

    ///////////////////////////////////////////
    ///  Write to file
    virtual void write(const char * filename) const = 0;

    ///////////////////////////////////////////
    ///  Write to iofits object
    virtual void write(Arroyo::iofits & iof) const = 0;

    ///////////////////////////////////////////
    ///  Print
    virtual void print(ostream & os, const char * prefix="") const = 0;

    ///////////////////////////////////////////
    ///  Emit a wavefront
    virtual diffractive_wavefront<float> emit(const diffractive_wavefront_header<float> & dwfh) const = 0;

    ///////////////////////////////////////////
    ///  Emit a wavefront
    virtual diffractive_wavefront<double> emit(const diffractive_wavefront_header<double> & dwfh) const = 0;

    ///////////////////////////////////////////
    ///  Get the direction vector from the emitter
    ///  towards a three_point tp
    virtual three_vector get_emission_vector(const three_point & tp) const = 0;

    ///////////////////////////////////////////
    ///  Factory to construct emitters from file
    static emitter * emitter_factory(const char * filename);

    ///////////////////////////////////////////
    ///  Factory to construct emitters from file
    static emitter * emitter_factory(const iofits & iof);

    ///////////////////////////////////////////
    ///  Factory to construct emitters from file
    static emitter * emitter_factory(const emitter * emtr);

  };

  ///
  /// A class to represent plane wave emitters
  ///

  class plane_wave_emitter :
    public emitter, 
    public Arroyo::three_vector {

    private:

    static const bool factory_registration;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("plane wave emitter"));};

    ///////////////////////////////////////////
    ///  template member function to dodge lack
    ///  of virtual template member functions
    template<class T>
      diffractive_wavefront<T> private_emit(const diffractive_wavefront_header<T> & dwfh) const;

    protected:

    ///////////////////////////////////////////
    ///  null constructor 
    plane_wave_emitter(){};

    public:

    ///////////////////////////////////////////
    ///  construct from the bits
    plane_wave_emitter(const three_vector & tv);

    ///////////////////////////////////////////
    ///  construct from file
    plane_wave_emitter(const char * filename);

    ///////////////////////////////////////////
    ///  copy constructor
    plane_wave_emitter(const plane_wave_emitter & pwe);

    ///////////////////////////////////////////
    ///  Construct from iofits object
    plane_wave_emitter(const iofits & iof) {
      this->read(iof);
    };

    ///////////////////////////////////////////
    ///  destructor
    ~plane_wave_emitter(){};

    ///////////////////////////////////////////
    ///  operator=
    plane_wave_emitter & operator=(const plane_wave_emitter & pwe);

    ///////////////////////////////////////////
    ///  read from file
    void read(const char * filename);

    ///////////////////////////////////////////
    ///  read from iofits object
    void read(const Arroyo::iofits & iof);

    ///////////////////////////////////////////
    ///  write to file
    void write(const char * filename) const;

    ///////////////////////////////////////////
    ///  write to iofits object
    void write(Arroyo::iofits & iof) const;

    ///////////////////////////////////////////
    ///  Print
    void print(ostream & os, const char * prefix="") const;

    ///////////////////////////////////////////
    ///  emit a wavefront
    diffractive_wavefront<float> emit(const diffractive_wavefront_header<float> & dwfh) const {
      return this->private_emit(dwfh);
    };

    ///////////////////////////////////////////
    ///  emit a wavefront
    diffractive_wavefront<double> emit(const diffractive_wavefront_header<double> & dwfh) const {
      return this->private_emit(dwfh);
    };

    ///////////////////////////////////////////
    ///  Get the direction vector from the emitter
    ///  towards a three_point tp
    three_vector get_emission_vector(const three_point & tp) const {
      return(*this);
    }

    ///////////////////////////////////////////
    ///  operator==
    friend bool operator==(const plane_wave_emitter & pwe1, const plane_wave_emitter & pwe2);

  };

  bool operator!=(const plane_wave_emitter & pwe1, const plane_wave_emitter & pwe2);

  ///
  /// A class to represent spherical wave emitters
  ///

  class spherical_wave_emitter :
    public emitter, 
    public Arroyo::three_point {

    private:

    static const bool factory_registration;

    Arroyo::three_vector emission_direction;

    ////////////////////////////
    ///  Return a name unique to the class
    string unique_name() const {return(string("spherical wave emitter"));};

    ///////////////////////////////////////////
    ///  template member function to dodge lack
    ///  of virtual template member functions
    template<class T>
      diffractive_wavefront<T> private_emit(const diffractive_wavefront_header<T> & dwfh) const;

    protected:

    ///////////////////////////////////////////
    ///  null constructor 
    spherical_wave_emitter(){};

    public:

    ///////////////////////////////////////////
    ///  construct from the bits
    spherical_wave_emitter(const three_point & tp);

    ///////////////////////////////////////////
    ///  construct from file
    spherical_wave_emitter(const char * filename);

    ///////////////////////////////////////////
    ///  copy constructor
    spherical_wave_emitter(const spherical_wave_emitter & swe);

    ///////////////////////////////////////////
    ///  Construct from iofits object
    spherical_wave_emitter(const iofits & iof) {
      this->read(iof);
    };

    ///////////////////////////////////////////
    ///  destructor
    ~spherical_wave_emitter(){};

    ///////////////////////////////////////////
    ///  operator=
    spherical_wave_emitter & operator=(const spherical_wave_emitter & swe);

    ///////////////////////////////////////////
    ///  read from file
    void read(const char * filename);

    ///////////////////////////////////////////
    ///  read from iofits object
    void read(const Arroyo::iofits & iof);

    ///////////////////////////////////////////
    ///  write to file
    void write(const char * filename) const;

    ///////////////////////////////////////////
    ///  write to iofits object
    void write(Arroyo::iofits & iof) const;

    ///////////////////////////////////////////
    ///  Print
    void print(ostream & os, const char * prefix="") const;

    ///////////////////////////////////////////
    ///  Get the direction vector from the emitter
    ///  towards a three_point tp
    three_vector get_emission_vector(const three_point & tp) const {
      if(operator==(tp, static_cast<const three_point>(*this))){
	std::cerr << "spherical_wave_emitter::get_emission_vector error - "
		  << "cannot get an emission vector at the location of the "
		  << "spherical wave emitter\n";
	throw(string("spherical_wave_emitter::get_emission_vector"));
      }
      three_vector tv = tp - *this;
      tv = tv*(1/tv.length());
      return(tv);
    };

    ///////////////////////////////////////////
    ///  emit a wavefront
    diffractive_wavefront<float> emit(const diffractive_wavefront_header<float> & dwfh) const {
      return this->private_emit(dwfh);
    };

    ///////////////////////////////////////////
    ///  emit a wavefront
    diffractive_wavefront<double> emit(const diffractive_wavefront_header<double> & dwfh) const {
      return this->private_emit(dwfh);
    };

    ///////////////////////////////////////////
    ///  operator==
    friend bool operator==(const spherical_wave_emitter & swe1, const spherical_wave_emitter & swe2);

  };

  bool operator!=(const spherical_wave_emitter & swe1, const spherical_wave_emitter & swe2);

  template<class T>
  diffractive_wavefront<T> plane_wave_emitter::private_emit(const diffractive_wavefront_header<T> & dwfh) const {
    
    //if(cross_product(*this, dwfh.z()).length()>three_frame::precision){
    if(cross_product(*this, dwfh.z()).length()>1e-6){
      cerr << "plane_wave_emitter::private_emit error - "
	   << "three_frame of diffractive wavefront header misaligned with the emission direction of this emitter\n";
      this->three_vector::print(cerr, "emitter vector ");
      dwfh.z().print(cerr, "header vector ");
      cerr << "cross product " << cross_product(*this, dwfh.z()).length() << endl;
      throw(string("plane_wave_emitter::emit"));
    }
    diffractive_wavefront<T> dwf(dwfh);

    // In the future, when the emitter class accounts for radiosity, this will be more elaborate
    complex<T> c(1,0);
    dwf += c;
    return(dwf);
  }

  template<class T>
  diffractive_wavefront<T> spherical_wave_emitter::private_emit(const diffractive_wavefront_header<T> & dwfh) const {
    
    diffractive_wavefront<T> dwf(dwfh);

    // Set the curvature of the wavefront
    if((*this - dwfh).length()<three_frame::precision){
      cerr << "spherical_wave_emitter::private_emit error -"
	   << "cannot emit wavefront at the location of the spherical wave emitter since curvature would be infinite\n";
      throw(string("spherical_wave_emitter::private_emit"));
    }
    dwf.set_curvature(1/(*this - dwfh).length());

    // In the future, when the emitter class accounts for radiosity, this will be more elaborate
    complex<T> c(1,0);
    dwf += c;
    return(dwf);
  }
}

#endif
