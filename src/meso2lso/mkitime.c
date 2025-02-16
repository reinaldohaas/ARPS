/*  This c subroutine is an interface between
    FORTRAN and the C library routine mktime
    Keith Brewster March 1995 */
#include <time.h>
#ifdef UNDERSCORE
#define mkitime mkitime_
#endif
void mkitime(int *iyr,int *imon,int *iday,int *ihour,
                int *imin, int *isec, int *itime)
{
  struct tm time_str;
  extern time_t timezone;

  time_str.tm_year       = *iyr - 1900;
  time_str.tm_mon        = *imon - 1;
  time_str.tm_mday       = *iday;
  time_str.tm_hour       = *ihour;
  time_str.tm_min        = *imin;
  time_str.tm_sec        = *isec;
  time_str.tm_isdst      = 0;

  *itime=mktime(&time_str);
  *itime=*itime-timezone;
}
