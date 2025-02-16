/*
*
*     ##################################################################
*     ##################################################################
*     ######                                                      ######
*     ######                      88D2ARPS                        ######
*     ######                                                      ######
*     ######                     Developed by                     ######
*     ######    Center for Analysis and Prediction of Storms      ######
*     ######                University of Oklahoma                ######
*     ######                                                      ######
*     ##################################################################
*     ##################################################################
*
*     fakeio.c
*
*#######################################################################
*
*     PURPOSE:
*
*     Emulates the behavior of the NSSL WSR-88D realtime and 
*     Archive II tape reading software.  Arrays filled with fake
*     data are passed to the calling program as if they were
*     being read from the realtime buffer or tape stream.
*
*     To use: Replace links to NSSL 88D libraries with 
*             object code from this file.
*
*     Some shortcomings of version 1.0:
*        A few calls are not supported.
*        Set_radar_name doesn't read the radar table and change location.
*        Storm mode only is emulated.
*        Runs faster than realtime (sleeps can be inserted in read_radial).
*
*     Things for the future:
*        Address shortcomings, above.
*        More sophisticated fake data, such as read from a file,
*        analytic field or even a cloud model.
*
*#######################################################################
*
*     AUTHOR: Keith Brewster
*     March, 1995   version 1.0
*     kbrewster@ou.edu
*
*     MODIFICATION HISTORY:
*
*#######################################################################
*
*/

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include <time.h>

struct tm rtime, *radtime;
struct count_struct {
     int iazi;
     int iangle;
     int itilt;
     int iscan;
} radcount;

struct loc_struct {
     char rsite[80];
     char radid[4];
     int latitude;
     int longitude;
     int altitude;
} radinfo;

/* ****  End of structure definitions **** */

/*     set_radar_name    */

/*     Assigns radar name in location structure
       but does not look up the lat,lon,elev like the
       real radar reader.                               */

void set_radar_name(site)
char site[80];
{
      extern struct loc_struct radinfo;
      strcpy(radinfo.rsite,site);
}

/*     get_field_number  */

/*    Assigns field number based on input string.
      Must agree with action of get_data_field.  */

int get_field_num(fchar)
char fchar[3];
{
      if(strcmp(fchar,"DBZ") == 0)
        return(0);
      else if(strcmp(fchar,"VEL") == 0)
        return(1);
      else if(strcmp(fchar,"SPW") == 0)
        return(2);
      else if(strcmp(fchar,"SNR") == 0)
        return(3);
      else
        return(-1);
}

/*     radar_init  */

/*    Initializes the radar counters.
      prints a warning message that fake data is being used.  */

void radar_init(drive)
char drive[80];
{
      extern struct count_struct radcount;
      extern struct tm rtime, *radtime;

      printf("\n\n *************FAKE DATA***********************\n\n");
      printf("  Initializing a2io_fake for file testing only.\n");
      printf("\n\n *************FAKE DATA***********************\n\n");

      radtime =&rtime;

      radcount.iazi=0;
      radcount.iangle=100;
      radcount.itilt=1;
      radcount.iscan=1;
}

/*     read_radial    */

/*     Emulates the act of reading a radial.
       Gets the time from the system and increments
       counters that pretend that the radar is spinning
       and going through a normal sequence of elev angles.

       Unlike the true read_radial, data are not buffered
       for use by get_data_field, rather get_data_field makes
       things up as it goes along                       */

/*     These are the tilts for the WSR-88D storm mode   */
/*     The precise behavior of the 0.5 degree scans     */
/*     is not replicated                                */

#define MAXTILT 16
static int kelev[MAXTILT] =
         {50,50,50,150,240,340,430,530,620,750,
          870,1000,1200,1400,1670,1950};

int read_radial()
{
      extern struct tm rtime, *radtime;
      extern struct count_struct radcount;
      extern int kelev[];
      int istat;
      long curtime[2], *tp;
      long timez[2], *tzp;

      tp = curtime;
      tzp = timez;

      radcount.iazi=radcount.iazi+100;
      if ( radcount.iazi > 35999 ) {
        radcount.iazi=radcount.iazi-36000;
        radcount.itilt++;
        if ( radcount.itilt > MAXTILT) {
          radcount.itilt=1;
          radcount.iscan++;
        }
      }
      radcount.iangle=kelev[radcount.itilt - 1];
      istat=gettimeofday(tp,tzp);
      radtime=gmtime(tp);
      return(0);
}
int get_first_gate(ifield)
int ifield;
{
    if( ifield == 0 )
      return(1000);
    else if ( ifield == 1)
      return(1000);
    else if ( ifield == 2)
      return(1000);
    else if ( ifield == 3)
      return(1000);
    else
      return(-999);
}
int get_gate_spacing(ifield)
int ifield;
{
    if( ifield == 0 )
      return(1000);
    else if ( ifield == 1)
      return(250);
    else if ( ifield == 2)
      return(250);
    else if ( ifield == 3)
      return(1000);
    else
      return(-999);
}

/*     get_data_field     */

/*     returns fake data

       Some day this could read from a storm model or other
       more realistic fake data generator.

       Reflectivity data are assigned an increasing function of range
       Velocity data are assigned a value like a constant horiz wind
       Spectrum Width data are constant = 0.5 ms-1
       SNR data are constant = 40.                                */

int get_data_field(ifield,field,nsize)
int ifield;
float field[];
int nsize;
{
    extern struct count_struct radcount;
    float twopi,dtr;
    float sinaz,cosaz,dxdr,dydr,dsdr;
    float mpgate,mgate0,x0,y0,xscale,yscale,x,y,u,v;
    int i;
    float *ptr;

    ptr=field;
 
    if( ifield == 0 )
      for (i = 0; i < nsize; i++){
        *ptr = (0.1 * (float) i);
        ptr++;
      }
    else if( ifield == 1 )
      for (i = 0; i < nsize; i++){
/*      *ptr =-abs (20.0 * sin (dtr * 0.01 * (float) radcount.iazi) );*/
/*      *ptr = 20.0 * sin (dtr * 0.01 * (float) radcount.iazi);*/
        dtr=0.017453292;
        twopi = 6.283185307;
        mpgate=250.;
        mgate0=1000.;
        x0 = 15000.;
        y0 = 15000.;
        xscale = twopi / 50000.;
        yscale = twopi / 50000.;
        cosaz = cos (dtr * 0.01 * (float) radcount.iazi);
        sinaz = sin (dtr * 0.01 * (float) radcount.iazi);
        dxdr = mpgate * cosaz ;
        dydr = mpgate * sinaz ;
        dsdr = cos (dtr * 0.01 * (float) radcount.iangle);
        for (i = 0; i < nsize; i++){
          x = x0 + dxdr * (mgate0 + (float) i);
          y = y0 + dydr * (mgate0 + (float) i);
          u = 10.0 * sin (yscale * y);
          v =-10.0 * sin (xscale * x);
          *ptr = 10.0 * dsdr * (u * sinaz + v * cosaz);
          ptr++;
        }
      }
    else if( ifield == 2 )
      for (i = 0; i < nsize; i++){
        *ptr = 0.5;
        ptr++;
      }
    else if( ifield == 3 )
      for (i = 0; i < nsize; i++){
        *ptr = 40.;
        ptr++;
      }
    else
      return(-1);

    return(0);
}

/*     get_year    */

/*     returns year 
       time is controlled by read_radial */

int get_year()
{
      extern struct tm *radtime;
      return(radtime->tm_year);
}

/*     get_month   */

/*     returns number of month (1-12)
       time is controlled by read_radial */

int get_month()
{
      extern struct tm *radtime;
      int imon;
      imon=radtime->tm_mon + 1;
      return(imon);
}

/*     get_day   */

/*     returns day
       time is controlled by read_radial */

int get_day()
{
      extern struct tm *radtime;
      return(radtime->tm_mday);
}

/*     get_hour   */

/*     returns hour
       time is controlled by read_radial */

int get_hour()
{
      extern struct tm *radtime;
      return(radtime->tm_hour);
}

/*     get_min   */

/*     returns minute 
       time is controlled by read_radial */

int get_min()
{
      extern struct tm *radtime;
      return(radtime->tm_min);
}

/*     get_sec   */

/*     returns seconds 
       time is controlled by read_radial */

int get_sec()
{
      extern struct tm *radtime;
      return(radtime->tm_sec);
}

/*     get_azi   */

/*     returns radar azimuth angle
       angle is controlled by read_radial  */

int get_azi()
{
      extern struct count_struct radcount;
      return(radcount.iazi);
}

/*     get_fixed_angle   */

/*     returns radar elevation angle
       angle is controlled by read_radial  */

int get_fixed_angle()
{
      extern struct count_struct radcount;
      return(radcount.iangle);
}

/*     get_scan   */

/*     returns radar scan number
       number is controlled by read_radial  */

int get_scan()
{
      extern struct count_struct radcount;
      return(radcount.iscan);
}

/*     get_tilt   */

/*     returns radar tilt number
       number is controlled by read_radial  */

int get_tilt()
{
      extern struct count_struct radcount;
      return(radcount.itilt);
}

/*     get_status   */

/*     status returned is fixed at 1,
       status is always good in radar nirvana */

int get_status(index)
int index;
{
      return(1);
}

/*     get_rt_mode   */

/*     rt_mode returned is fixed at 1, real-time /

int get_rt_mode()
{
      return(1);
}

/*     get_nyquist   */

/*     nyquist velocity returned is fixed at 35.0 ms-1  */

 
int get_nyquist()
{
      return(3500);
}

/*     get_vcp   */

/*     vcp returned is fixed at 41  */

int get_vcp()
{
      return(41);
}

/*     get_latitude   */

/*     latitude returned is fixed at 35.23N  */

int get_latitude()
{
/*    return(3523000); */
/*    return(3674000); */
      return(3523000);
}

/*     get_longitude  */

/*     longitude returned is fixed at 97.46W  */

int get_longitude()
{
/*    return(9746000); */
/*    return(9812750); */
      return(9746000);   
}

/*     get_altitude  */

/*     altitude returned is fixed at 381 m */

int get_altitude()
{
      return(381);
}

/*  get_number_of_gates
 *  
 *  Need to be checked for correctness -- Y. Wang.
 *  Just added for compiling 88d2arps_fake.
 *
*/

int get_number_of_gates(int n) {
  return 0;
}

/*  Always return 0. */

int get_restart_flag() {
  return 0;
}

/* Dummy subroutine for SOLOIO */
int init_solo(char *radar_name,int i_lat, int i_lon, int i_alt,int i_vcp,
	        int iyear,int imon, int iday,int ihr, int imin, double second,
	        int i_scan, int i_tilt,int i_angle,
	        int n_ref_gates,int n_vel_gates, int n_max_rays) {
	return 0;
}

int add_radial_to_sweep(int i_tilt,float *eleva_ptr,float *azi_ptr,int n_azim,
		               int ref_ok, int vel_ok,int n_ref_gates,int n_vel_gates,int n_spw_gates,
	                  float *v_nyq,
	                 int iyear,int imon, int iday,int ihr, int imin, double second,
	                int rfrst_ref,int gsp_ref,int rfrst_vel, int gsp_vel,
	               float *ref_ptr,float *vel_ptr,float *spw_ptr) {
	return 0;
}

int write_sweep() {

	return 0;
}

void clean_up(void) {
}
