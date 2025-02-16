/*  This c subroutine is an interface between
    FORTRAN and the C library routine mktime
    Keith Brewster March 1995 */
#include <stdio.h>
#include <time.h>
#ifdef UNDERSCORE
#define cvtitoi cvtitoi_
#endif
void cvtitoi(int *itime, int *iyr,int *imon,int *iday,int *ihour,
             int *imin, int *isec)
{
  struct tm time_str, *tmptr;
  time_t timet;

  tmptr=&time_str;
  timet=(time_t)(*itime);
/*  printf(" timet = %ld, input time = %ld.\n",timet,*itime); */
  tmptr=gmtime(&timet);

  *iyr = (tmptr->tm_year + 1900);
  *imon = (tmptr->tm_mon + 1);
  *iday = (tmptr->tm_mday);
  *ihour = (tmptr->tm_hour);
  *imin = (tmptr->tm_min);
  *isec = (tmptr->tm_sec);

/*  printf(" iyr = %ld, iday = %ld, imin = %ld.\n",*iyr,*iday,*imin); */
}
