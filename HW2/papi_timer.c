/***************************************************************
   File:     papi_timer.c

   Function: get_cur_time()

   Return:   the time since some arbitrary starting point in seconds
             as a double value.
             PAPI timers use the most accurate timers available 
             on the platform in use.

   Note:     must be used on the platform where PAPI is installed 
             (e.g. torc1), and remember to call PAPI_library_init() 
             first.

****************************************************************/

#include <papi.h>

#include <sys/time.h>

double get_cur_time_papi() {
  return PAPI_get_real_usec() / 1000000.0;
} 

double get_cur_time() {
      struct timeval   tv;
        struct timezone  tz;
          double cur_time;
            
            gettimeofday(&tv, &tz);
              cur_time = tv.tv_sec + tv.tv_usec / 1000000.0;
                
                return cur_time;
} 
