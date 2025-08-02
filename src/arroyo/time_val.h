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

#ifndef TIME_VAL_H
#define TIME_VAL_H

#include <time.h>
#include <iomanip>
#include "iofits.h"

namespace Arroyo {

  time_t my_timegm (struct tm *tm);

  class time_val {

  protected:

    time_t sec;
    suseconds_t usec;

  public:

    time_val(){
      this->sec = -1;
      this->usec = -1;
    };

    time_val(const time_val & tv){
      this->operator=(tv);
    };

    time_val(const time_t & tt){
      this->sec = tt;
      this->usec = 0;
    };

    time_val(struct timeval tval){
      this->sec = tval.tv_sec;
      this->usec = tval.tv_usec;
    };

    time_val(struct tm & utc_tm){
      this->sec = my_timegm(&utc_tm);
      this->usec = 0;
    };

    time_val(time_t sec,
	     suseconds_t usec){
      this->sec = sec;
      this->usec = usec;
    };

    time_val & operator=(const time_val & tv){
      if(this==&tv)
	return(*this);
      this->sec = tv.sec;
      this->usec = tv.usec;
    };

    struct timeval get_timeval() const {
      struct timeval tv;
      tv.tv_sec = this->sec;
      tv.tv_usec = this->usec;
      return(tv);
    }

    void get_tm(struct tm & time_tm) const {
      gmtime_r(&(this->sec), &time_tm);
    }

    double get_time() const {
      return(this->sec + 1e-6*this->usec);
    }

    time_t get_secs() const {
      return(this->sec);
    }

    time_t get_usecs() const {
      return(this->usec);
    }

    string get_timestr() const {
      struct tm frame_time;
      this->get_tm(frame_time);
      string st(asctime(&frame_time));
      st.erase(st.length()-1);
      return(st);
    }

    void read(const iofits & iof) {
      long tmp_sec;
      long tmp_usec;
      string comment;
      iof.read_key("TIMESEC", 
		    tmp_sec, 
		    comment);
      iof.read_key("TIMEUSEC", 
		    tmp_usec, 
		    comment);
      this->operator=(time_val((time_t)tmp_sec, (suseconds_t)tmp_usec));
    }

    void write(iofits & iof) const {
      timeval tv = this->get_timeval();
      iof.write_key("TIMESEC", 
		    tv.tv_sec, 
		    string("timestamp (secs)"));
      iof.write_key("TIMEUSEC", 
		    tv.tv_usec, 
		    string("timestamp (microsecs)"));
    }

    // operator overloads
    friend bool operator==(const time_val & tv1, 
			   const time_val & tv2);

    friend time_val operator+(const time_val & tv1, 
			      const time_val & tv2);

    friend time_val operator-(const time_val & tv1, 
			      const time_val & tv2);

    friend bool operator>(const time_val & tv1, 
			  const time_val & tv2);

    friend bool operator>=(const time_val & tv1, 
			   const time_val & tv2);

    friend bool operator<(const time_val & tv1, 
			  const time_val & tv2);

    friend bool operator<=(const time_val & tv1, 
			   const time_val & tv2);

    friend std::ostream & operator<<(std::ostream & os,
				     const time_val & tv);

  };

}
#endif
