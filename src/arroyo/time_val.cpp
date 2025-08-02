#include "time_val.h"

namespace Arroyo {

  time_t my_timegm (struct tm *tm) {
    time_t ret;
    char *tz;

    tz = getenv("TZ");
    setenv("TZ", "", 1);
    tzset();
    ret = mktime(tm);
    if (tz)
      setenv("TZ", tz, 1);
    else
      unsetenv("TZ");
    tzset();
    return ret;
  }


  time_val operator+(const time_val & tv1,
		     const time_val & tv2){
    time_t sum_sec = tv1.sec+tv2.sec;
    suseconds_t sum_usec = tv1.usec+tv2.usec;
    while(sum_usec>1e6){
      sum_sec++;
      sum_usec = sum_usec - (suseconds_t)1e6;
    }
    return(time_val(sum_sec,
		    sum_usec));
  }

  time_val operator-(const time_val & tv1,
		     const time_val & tv2){
    double diff_sec = difftime(tv1.sec,tv2.sec);
    suseconds_t diff_usec = tv1.usec-tv2.usec;
    while(diff_usec<0){
      diff_sec--;
      diff_usec += (suseconds_t)1e6;
    }
    return(time_val((time_t)diff_sec,
		    diff_usec));
  }

  bool operator==(const time_val & tv1,
		  const time_val & tv2){
    if(tv1.sec==tv2.sec &&
       tv1.usec==tv2.usec)
      return(true);
    return(false);
  }

  bool operator!=(const time_val & tv1,
		  const time_val & tv2){
    return(!operator==(tv1,tv2));
  }

  bool operator>(const time_val & tv1, 
		 const time_val & tv2){

    if((tv1.sec>tv2.sec) ||
       (tv1.sec==tv2.sec && 
	tv1.usec>tv2.usec))
      return(true);
    return(false);
  }

  bool operator>=(const time_val & tv1, 
		  const time_val & tv2){
    if(operator>(tv1,tv2) ||
       operator==(tv1,tv2))
      return(true);
    return(false);
  }

  bool operator<(const time_val & tv1, 
		 const time_val & tv2){
    return(!operator>=(tv1,tv2));
  }

  bool operator<=(const time_val & tv1, 
		  const time_val & tv2){
    return(!operator>(tv1,tv2));
  }

  std::ostream & operator<<(std::ostream & os,
			    const time_val & tv){
    return(os << std::setw(12) << tv.sec << std::setw(12) << tv.usec);
  }
}
