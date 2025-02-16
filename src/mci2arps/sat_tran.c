/*********************************************************************

   ROUTINE: sat_tran.c
   PROGRAMMER: Dan Vietor, D. Farmer, Eric A. Smith & D. R. Phillips
   PROGRAM TYPE: WXP_LIBRARY
   VERSION: 1.0
   WXP VERSION: 4.8
   DATE: 920515

   DESCRIPTION: General routines to convert satellite navigation to 
      earth coordinates.

   modifications by D. Farmer
   Unidata Program Center in Boulder, CO

   Adapted C code from original Fortran code written by
   Eric A. Smith & D. R. Phillips
   Dept. of Atmospheric Science
   CSU/Foothills Campus
   Fort Collins, CO 80523

   COMPUTER:
      IBM RS/6000 running AIX 3.1 C Compiler
      IBM RT/PC running AIX 2.2.1 C Compiler
      IBM RT/PC running BSD 4.3 C Compiler
      IBM AT compatible running MSDOS 3.x using the
         Microsoft C 6.0 Compiler.
      SUN 3 and 4 running SUNOS 4.1
      DEC VaxStation running VMS 5.3 and GKS
      Interprocess communications based on either ATT UNIX System V
         Message queues or BSD Sockets.
      Asynchronous data stream function interface based on either
         System V ioctl, BSD 4.3 ioctl or Blaise Asynch Manager
         for MSDOS 3.x.
      Graphics calls based on ANSI C language Graphical Kernel
         System binding (revision dated July 1988)

   REVISIONS:                                    DATE:       INIT:
Version 1                                        920515      DEV

*********************************************************************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mc_area.h"
#include "wxp.h"

#define GOES 1
#define GVAR 2
#define PSTEREO 10

struct pt {
   double x,y,z;
   };

struct sp_pt {
   double lat,lon,hgt;
   };

struct angle {
   double ang;
   double sin,cos;
   };

struct datim {
   long date;        /* Date YYDDD */
   double tim;       /* Time HH.HHH */
   };
/*
********** GOES Parameters *************
*/
#define EL_LL 1
#define LL_EL 2
#define RADE_EQUA    6378.388
#define RADE_POLE    6356.912
#define RADE_MEAN    6371.221
#define SOLAR_YR     365.24219879
#define SIDER_YR     366.24219879
#define SOLSID       (SIDER_YR/SOLAR_YR)
#define REF_DATE     90001L
#define REF_TIME     171537L
#define PREC_VER_EQ  25781.
#define OBLIQ_ECLIP  23.45
#define NOMORB 42164.365      /* Nominal distance of datellite (km) */
#define REQ     6378.137
#define RPL     6356.7533
#define REQRPL2 ((REQ*REQ)/(RPL*RPL))
#define FER     (RPL/REQ)
#define FER2    (1-RPL/REQ)
#define REQRPL3 (REQRPL2-1)
#define REQRPL4 (FER*FER*FER*FER-1)

int image_type = 0;

struct goes_orbit {
   int sat_type;    /* Satellite type - iosat
                       set positive for initializing new satellite type
                       set negative for retaining old satellite type with new 
                          orbit parms, type is then set positive */
   int anom_type;   /* Anomaly type - imort
                       0 for orbit anomaly given as mean anomaly (e.g. NASA)
                       1 for orbit anomaly given as true anomaly (e.g. ESA) */
   int sec_ord;     /* Secular theory order - iosec
                       0 for zero order secular perturbation theory
                       1 for first order secular perturbation theory
                       2 for second order secular perturbation theory */
   long epoch_date; /* Epoch date - iedate (yymmdd in calendar form)
                       date for which following orbital parameters are valid */
   long epoch_time; /* Epoch time - ietime (hhmmss in GMT)
                       date for which following orbital parameters are valid */
   double semi_maj; /* Semi-major axis (km) - semima
                       half the distance between two apses of apo-focus and 
                       peri-focus */
   double eccent;   /* Eccentricity of earth orbit (unitless) - oeccen
                       degree of ellipticity of orbit */
   double inclin;   /* Orbit inclination (degrees) - orbinc
                       angle between the orbit and equatorial planes */
   double anomaly;  /* Orbit anomaly at epoch time (degrees) - oanoml
                       angle in orbital plane between peri-focus and 
                       satellite position given as either a mean anomaly 
                       or a true anomaly */
   double perigee;  /* Argument of perigee at epoch time (degrees) - perige
                       angle in orbit plane from ascending node to peri-focus */
   double asc_node; /* Right ascension of ascending node at epoch 
                       time (degrees) - asnode
                       angle in equatorial plane between vernal equinox 
                       (principle axis) and northward equator crossing */
   double period;   /* Period (minutes) - period
                       statement of Kepler's Third Law
                       this parameter calculated in satpos */
   double anom_per; /* Anomalistic period (minutes) - aperod
                       time between the passage from one peri-focus to the next
                       this parameter calculated in satpos */
   double nodal_per;/* Nodal period (minutes) - eperod
                       time between the passage from one equator passing 
                       to the next this parameter calculated in eqcros */
   } g_orb;

struct goes_navigate {
   long num_lines;  /* Total number of scan lines per frame and 
                       sensors (lines) - lintot */
   double deg_line; /* Angle of north-south frame scan (degrees) - deglin */
   long num_elem;   /* Total number of elements per sscan line 
                       (elements) - num_elem */
   double deg_elem; /* Angle of east-west frame scan (degrees) - degele */
   double spin_rate;/* Spin rate (milliseconds/line) - spinra */
   double declin;   /* Spin axis declination (degrees) - declin */
   double right_asc;/* Spin axis right ascension in celestial coor 
                       system (degrees) - rascen */
   long center_line;/* Picture center line (lines) - piclin */
   double prec_rate;/* Precession rate (degrees/day) - prerat */
   double prec_dir; /* Precession direction (degrees) - predir */
   double pitch;    /* North-south misalignment in principle axis of 
                       camera (degrees) - pitch */
   double yaw;      /* Skew misalignment in principle axis of 
                       camera (degrees) - yaw */
   double roll;     /* East-west misalignment in principle axis of 
                       camera (degrees) - roll */
   } g_nav;

static struct datim sat_time;
static double gamma_init;  /* Initial gamma value */
static double gamma_dot;   /* Derivative of gamma */
static double atms_hght;   /* Atmospheric height */
static double earth_rot;   /* Earth rotation speed */
static double prec_veq;
static struct datim ref;
static struct mc_area *sat_image;
static struct pt sat;
static struct sp_pt sp_sat;
static float lres,eres;
static float lcor,ecor;
/*
************ GVAR Parameters ************
*/
struct gvar_navigate {
  char sttype[TYPELEN];     /*   1     STTYPE = Satellite type */
  long idntfr;     /*   2     IDNTFR = */
  long imcact;     /*   3     IMCACT = IMC active flag */
  float reflon;    /*   6    +REFLON = Reference longitude */
  float refdis;    /*   7    +REFDIS = Reference distance from nominal */
  float reflat;    /*   8    +REFLAT = Reference latitude */
  float refyaw;    /*   9    +REFYAW = Reference yaw */
  float ratrol;    /*  10    +RATROL = Reference attitude roll */
  float ratptc;    /*  11    +RATPTC = Reference attitude pitch */
  float ratyaw;    /*  12    +RATYAW = Reference attitude yaw */
  double epoch;    /*  13-14  ETIME  = Epoch time */
  float edtime;    /*  15    +EDTIME = Delta from epoch time */
  float imcrol;    /*  16    +IMCROL = Image motion compensation roll */
  float imcptc;    /*  17    +IMCPTC = Image motion compensation pitch */
  float imcyaw;    /*  18    +IMCYAW = Image motion compensation yaw */
  float ldr[13];   /*  19-31 +LDR    = Longitude delta from ref parameters */
  float rddr[11];  /*  32-42 +RDDR   = Radial distance delta from ref params */
  float dgl[9];    /*  43-51 +DGL    = Geocentric latitude delta parameters */
  float doy[9];    /*  52-60 +DOY    = Orbit yaw delta parameters */
  float solinc;    /*  62    +EXPTIM = Exponential start time from epoch */
  float exptim;    /*  62    +EXPTIM = Exponential start time from epoch */
  struct ang_wds {
     float expn[3];    /*   1- 3 */
     long nsinus;      /*   4    */
     float sinus[30];  /*   5-34 */
     long nmsin;       /*  35    */
     struct mono_sin {
        long ord;      /*  36    */
        long mord;     /*  37    */
        float sinus[3];/*  38-40 */
        } msin[4];
     } raawds,     /*  63-129 RAAWDS = Roll attitude angle words */
       paawds,     /* 130-184 PAAWDS = Pitch attitude angle words */
       yaawds,     /* 185-257 YAAWDS = Yaw attitude angle words */
       rmawds,     /* 258-312 RMAWDS = Roll misalignment angle words */
       pmawds;     /* 313-367 PMAWDS = Pitch misalignment angle words */
  double imgtim;   /* 368     IMGDAY = Image day value (YYDDD) */
                   /* 369     IMGTM  = Image time value (HHMMSS) */
  long imgsnd;     /* 370     IMGSND = Imager/sounder instrument flag */
  } gv_nav;

/*  ELCOMM include variables */
double bt[3][3];             /* Instrument to ECEF coordinates transformation */
double q3;                   /* Used in function lpoint */
double pitch,roll,yaw;       /* Pitch,roll,yaw angles of instrument (rad) */
float pma,rma;               /* Pitch,roll misalignments of instrument (rad) */

/*  INSTCO include variables */
long inc_max[2];             /* Number of increments per cycle */
float elev_bnds[2];          /* Bounds in elevation (radians) */
float scan_bnds[2];          /* Bounds in scan angle (radians) */
float elev_inc[2];           /* Change in elevation angle per increment (rad) */
float scan_inc[2];           /* Change in scan angle per increment (radians) */
float elev_dln[2];           /* Elevation angle per detector line (radians) */
float scan_pix[2];           /* Scan angle per pixel (radians) */

/*  GVRCOM common variables */
int itype;
int instr;
double sublat;
double sublon;

/*  SAVCOM common variables */
struct pt xs;               /* Normalized S/C position in ECEF coordinates */
double b[3][3];             /* Spacecraft to ECEF coordinates transformation */
double dr;
struct angle phi;
struct angle psi;

/*********************************************************
   NUMYR - Number of days in a year
*********************************************************/
long numyr( year )

long year;

{
if( year < 50 )
   year += 2000;
else
   year += 1900;
if(( year % 4 == 0 && year % 100 != 0 ) || year % 400 == 0 )
   return 366;
else
   return 365;
}

/*********************************************************
   NUMDY - Number of days between 2 dates
*********************************************************/
long numdy(date1,date2)

long date1,date2;

{
long year1,day1,year2,day2,temp,days;
   
year1 = date1/1000;
day1 = date1 % 1000;
year2 = date2/1000;
day2 = date2 % 1000;
temp = 1;
/*
   Switch values if date1 > date2
*/
if(( year1 > year2 ) ||
   ( year1 == year2 && day1 > day2)){
   temp = year2;
   year2 = year1;
   year1 = temp;
   temp = day2;
   day2 = day1;
   day1 = temp;
   temp = -1;
   }
for( days = 0; year1 < year2; year1++ ){
   days += numyr(year1) - day1 + 1;
   day1 = 1;
   }
days += day2 - day1;
days *= temp;
return days;
}

/*********************************************************
   TIMDIF - Number of minutes between 2 dates
*********************************************************/
double timdif( date1,time1,date2,time2 )

long date1,date2;
double time1,time2;

{
return 1440.0*numdy(date1,date2) + 60.0*(time2 - time1);
}

/*********************************************************
   GEOLAT - Geoedtic-geocentric latitude conersion
*********************************************************/
double geolat( type,lat )

int type;
double lat;

{
double fac;
   
fac = (RADE_EQUA - RADE_POLE) / RADE_EQUA;
fac = (1.0 - fac)*(1.0 - fac);
lat = DRC * lat;
/*
   Geodetic to geocentric conversion 
*/
if( type == 1 )
   return atan(tan(lat)*fac)*RDC;
/*
   Geocentric to geodetic conversion 
*/
else
   return atan(tan(lat)/fac)*RDC;
}

/*********************************************************
   MDCON - Conversion between YYMMDD to YYDDD for
   Julian date conversions
*********************************************************/
long mdcon( type,date )

int type;
long date;

{
long year,month,day,leap,julian,tyear;
static int num[]   = { 0,31,59,90,120,151,181,212,243,273,304,334,365 };
static int num_l[] = { 0,31,60,91,121,152,182,213,244,274,305,335,366 };
/*
   YYMMDD to YYDDD conversion
*/
if( type == 1 ){
   year = date/10000;
   date %= 10000;
   month = date/100;
   date %= 100;
   day = date;
   if( month < 1 )
      month = 1;
   if( month > 12 )
      month = 12;
   if( year < 50 )
      tyear = 2000 + year;
   else
      tyear = 1900 + year;
   if(( tyear % 4 == 0 && tyear % 100 ) || tyear % 400 == 0 )
      julian = year * 1000 + num_l[month-1];
   else
      julian = year * 1000 + num[month-1];
   julian += day;
   return julian;
   }
/*
   YYDDD to YYMMDD conversion
*/
else { 
   year = date/1000;
   julian = date % 1000;
   if( year < 50 )
      tyear = 2000 + year;
   else
      tyear = 1900 + year;
   if(( tyear % 4 == 0 && tyear % 100 ) || tyear % 400 == 0 )
      leap = 1;
   else
      leap = 0;
   if( julian < 1 )
      julian = 1;
   if( leap && julian > 366 )
      julian = 366;
   else if( !leap && julian > 365 )
      julian = 365;
   if( leap ){
      for( month = 0; month < 12 && num_l[month] < julian; month++ );
      day = julian - num_l[month-1];
      }
   else {
      for( month = 0; month < 12 && num[month] < julian; month++ );
      day = julian - num[month-1];
      }
   return year*10000 + month*100 + day;
   }
}

/*********************************************************
   INT2DEG - Converts HH.MMSS to decimal degrees
*********************************************************/
double int2deg( ang )

long ang;

{
double deg;
int sign;

if( ang < 0 ){
   ang = -ang;
   sign = -1;
   }
else
   sign = 1;

deg = (double)(ang % 100) / 3600.;
ang /= 100;
deg += (double)(ang % 100) / 60.;
ang /= 100;
deg += ang;
return sign*deg;
}

/*********************************************************
   DEG2INT - Converts decimal degrees to HH.MMSS
*********************************************************/
long deg2int( ang )

double ang;

{
long lang, deg;
int sign;

if( ang < 0 ){
   ang = -ang;
   sign = -1;
   }
else
   sign = 1;

ang += 1.38888E-4;

lang = (long)(ang*3600);
deg = (long)ang;
lang -= deg*3600;
deg = deg * 100 + (lang/60);
lang -= (lang/60)*60;
deg = deg * 100 + lang;
return sign*deg;
}

/*********************************************************
   ROUND - Rounds to the nearest whole number
*********************************************************/
long round_long( x )

double x;

{
if( x < 0 )
   return (long)( x - 0.5 );
else if( x == 0 )
   return 0L;
else
   return (long)( x + 0.5 );
}

/************************************************************
   elev_scan2line_pix - Converts GVAR elevation and scan line
   information to line and element.
************************************************************/
void elev_scan2line_pix(elev,scan,line,elem)

float elev;       /* Elevation angle in radians */
float scan;       /* Scan angle in radians */
float *line;      /* Line number (OUTPUT) */
float *elem;      /* Element number (OUTPUT) */
{
/*
  Compute fractional line number
*/
*line = (elev_bnds[instr] - elev) / elev_dln[instr];
if( instr == 0 )
  *line = *line + 4.5;
else
  *line = *line + 2.5;
/*
  Compute fractional element number
*/
*elem = (scan_bnds[instr] + scan) / scan_pix[instr] + 1;
}

/*********************************************************
   sat_line_elem - Converts image line,element to satellite
   line,element.
*********************************************************/
void sat_line_elem( type, ilin,iele, olin,oele )

int type;
float ilin,iele, *olin,*oele;

{
if( type == 1 ){
   *oele = ecor + iele*eres;
   *olin = lcor + ilin*lres;
   }
else {
   *oele = (iele - ecor)/eres;
   *olin = (ilin - lcor)/lres;
   }
}

/*********************************************************
   SET_SAT_RES - Sets up the satellite resolution for
   line element conversion.
*********************************************************/
void set_sat_res( ilcor,iecor,ilres,ieres )

float ilcor,iecor,ilres,ieres;

{
lres = ilres;
eres = ieres;
lcor = ilcor;
ecor = iecor;
}
/*********************************************************
   SPHERE_CVT - Converts from cartesian to spherical 
   coordinates.
*********************************************************/
void sphere_cvt( type,pnt,sph_pnt )

int type;
struct pt *pnt;
struct sp_pt *sph_pnt;

{
double temp;
/*
   Convert to spherical coord's
*/
if( type == 1 ){
   temp = pnt->x * pnt->x + pnt->y * pnt->y;
   sph_pnt->lat = atan2(pnt->z, sqrt(temp)) * RDC;
   sph_pnt->lon = atan2(pnt->y, pnt->x) * RDC;
   sph_pnt->hgt = sqrt(temp + pnt->z * pnt->z);
   }
}

/************************************************************
   gatt_init - Initializes attitude and misalignment attributes
************************************************************/
void gatt_init( words, data )

struct ang_wds words;
long *data;

{
int i,j,k;

words.expn[0] = data[0]/10000000.;
words.expn[1] = data[1]/100.;
words.expn[2] = data[2]/10000000.;

words.nsinus = data[3];
for( i = 0; i < 30; i++ )
   words.sinus[i] = data[i+4]/10000000.;

words.nmsin = data[34];
for( i = 0, k = 35; i < 4; i++ ){
   words.msin[i].ord = data[k++];
   words.msin[i].mord = data[k++];
   for( j = 0; j < 3; j++ )
      words.msin[i].sinus[j] = data[k++]/10000000.;
   }
}

/************************************************************
   gatt - Calculates attitude and misalignment attributes
************************************************************/
double gatt(parms,wa,te)

struct ang_wds parms;
double wa;            /* Input solar orbit angle in radians */
double te;            /* Input exponential time delay from epoch (minutes) */

{
int i;
double att;

att = parms.expn[2];
/*
  Computes the exponential term
*/
if( te >= 0 ) 
   att += parms.expn[0] * exp(-te / parms.expn[1]);
/*
  Calculation of sinusoids
*/
for( i = 0; i < parms.nsinus; i++ )
   att += parms.sinus[2*i] * cos(wa*i+parms.sinus[2*i+1]);
/*
  Computes monomial sinusoids
*/
for( i = 0; i < parms.nmsin; i++ ){
   att += parms.msin[i].sinus[0] * 
          pow((wa - parms.msin[i].sinus[2]),(double)parms.msin[i].mord) * 
          cos((double)parms.msin[i].ord*wa+parms.msin[i].sinus[1]);
   }
return att;
}

/************************************************************
   inst2e - Calculates instrument to earth coordinate transformation
   matrix
************************************************************/
inst2e(roll,pitch,yaw,a,at)

double roll;      /* Roll angle in radians */
double pitch;     /* Pitch angle in radians */
double yaw;       /* Yaw angle in radians */
double a[3][3];   /* Spacecraft to ECEF coordinates transformation matrix */
double at[3][3];  /* Instrument to ECEF coordinates transformation matrix */

{
double rpy[3][3];
int i,j;
/*
  We compute instrument to body coordinates transformation
  matrix by using a small angle approximation of trigonometric
  functions of the roll, pitch and yaw.
*/
rpy[0][0] = 1 - .5 * (pitch*pitch + yaw*yaw);
rpy[0][1] = -yaw;
rpy[0][2] = pitch;
rpy[1][0] = yaw + pitch * roll;
rpy[1][1] = 1. - 0.5 * (yaw * yaw + roll * roll);
rpy[1][2] = -roll;
rpy[2][0] = -pitch + roll * yaw;
rpy[2][1] = roll + pitch * yaw;
rpy[2][2] = 1. - 0.5 * (pitch * pitch + roll * roll);
/*
  Multiplication of matrices a and rpy
*/
for( i = 0; i < 3; i++ )
   for( j = 0; j < 3; j++ )
       at[i][j] = a[i][0]*rpy[0][j] + a[i][1]*rpy[1][j] + a[i][2]*rpy[2][j];

return(0);
}

/************************************************************
   set_gvar_con - Initializes GVAR constants
************************************************************/
void set_gvar_con(instr, nadnsc, nadnsi, nadewc, nadewi)

long instr, nadnsc, nadnsi, nadewc, nadewi;

{
inc_max[0]   = 6136;
inc_max[1]   = 2805;
elev_inc[0]  = 8.E-6;
elev_inc[1]  = 17.5E-6;
scan_inc[0]  = 16.E-6;
scan_inc[1]  = 35.E-6;
elev_dln[0]  = 28.E-6;
elev_dln[1]  = 280.E-6;
scan_pix[0]  = 16.E-6;
scan_pix[1]  = 280.E-6;
elev_bnds[0] = 0.220896;
elev_bnds[1] = 0.22089375;
scan_bnds[0] = 0.24544;
scan_bnds[1] = 0.2454375;
/*
  New code because of change to elug; nadir position is available
  in the signal, so should compute 4 values above using them
  instead of having them hardwired

  New code from Kathy Kelly for sounder nav - 10/27
*/
if( nadnsc != 0 && nadnsi != 0 && nadewc != 0 && nadewi != 0 ){
   if( gv_nav.imgsnd == 0 )
       elev_bnds[instr] = (nadnsc*inc_max[instr]+nadnsi)*elev_inc[instr];
   else
       elev_bnds[instr] = ((9-nadnsc)*inc_max[instr]-nadnsi)*elev_inc[instr];

   scan_bnds[instr] = (nadewc*inc_max[instr]+nadewi)*scan_inc[instr];
   }
}

/************************************************************
   timex - General time conversion utility to convert year, 
   day, hour, minute and second to minutes from Jan 1, 1950
************************************************************/
double timex( year, day, hour, min, sec )

int year, day, hour, min;
float sec;

{
unsigned long temp;
/*
  Here we convert year and day of year to number of
  days from 0 hour UT, 1950 Jan. 1.0
  This conversion is based on an algorithm by Fliegel and Van
  Flandern, Comm. of ACM, Vol.11, No. 10, Oct. 1968 (P.657)
*/
temp = day + 1461 * (year + 4799) / 4 - 3 * ((year + 4899) / 100) / 4 - 2465022;
/*
  Compute time in minutes from January 1.0, 1950
*/
return (double)temp * 1440 + hour * 60 + min + sec / 60;
}

/************************************************************
   time50 - Converts file based date to epoch time in minutes
   from Jan 1, 1950
************************************************************/
double time50( epoch )

long epoch[2];

{
unsigned long temp;
int year, day, hour, min;
float sec;

temp = epoch[0];

hour = (temp % 16) * 10; temp >>= 4;
day = (temp % 16); temp >>= 4;
day += (temp % 16) * 10; temp >>= 4;
day += (temp % 16) * 100; temp >>= 4;
year = (temp % 16); temp >>= 4;
year += (temp % 16) * 10; temp >>= 4;
year += (temp % 16) * 100; temp >>= 4;
year += (temp % 16) * 1000; temp >>= 4;

temp = epoch[1];

sec = (temp % 16) * .001; temp >>= 4;
sec += (temp % 16) * .01; temp >>= 4;
sec += (temp % 16) * .1; temp >>= 4;
sec += (temp % 16); temp >>= 4;
sec += (temp % 16) * 10; temp >>= 4;
min = (temp % 16); temp >>= 4;
min += (temp % 16) * 10; temp >>= 4;
hour += (temp % 16); temp >>= 4;
/*
  Here we convert year and day of year to number of
  days from 0 hour UT, 1950 Jan. 1.0
  This conversion is based on an algorithm by Fliegel and Van
  Flandern, Comm. of ACM, Vol.11, No. 10, Oct. 1968 (P.657)
*/
temp = day + 1461 * (year + 4799) / 4 - 3 * ((year + 4899) / 100) / 4 - 2465022;
/*
  Compute time in minutes from January 1.0, 1950
*/
return (double)temp * 1440 + hour * 60 + min + sec / 60;
}

/************************************************************
   ELEV_LINE
************************************************************/
double elev_line(instr,line)

int instr;         /* Instrument code (0-imager, 1-sounder) */
double line;       /* Fractional line number */

{
if (instr == 1)
  return elev_bnds[instr] - (line - 4.5) * elev_dln[instr];
else
  return elev_bnds[instr] - (line - 2.5) * elev_dln[instr];
}

/************************************************************
   SCAN_ANGLE
************************************************************/
double scan_angle(instr,pixel)

int instr;         /* Instrument code (0-imager, 1-sounder) */
double pixel;      /* Fractional pixel number */

{
return (pixel - 1.) * scan_pix[instr] - scan_bnds[instr];
}

/************************************************************
   lpoint - Converts GVAR elevation and scan line to lat and lon
************************************************************/
int lpoint(elev,scan,lat,lon)

float elev;       /* Elevation angle (rad) */
float scan;       /* Scan angle (rad) */
float *lat;       /* Latitude in radians (OUTPUT) */
float *lon;       /* Longitude in radians (OUTPUT) */

{
/*
  Output status; 0 - point on the earth found,
                 1 - instrument points off earth
*/
struct angle alpha;
struct angle zeta;
struct pt g,g1,u;
double q1,q2,d,d1,h;
/*
  Computes trigonometric functions of the scan and elevation
  angles corrected for the roll and pitch misalignments
*/
alpha.ang = elev;
alpha.sin = sin( alpha.ang );
alpha.cos = cos( alpha.ang );

zeta.ang = scan;
zeta.sin = sin( zeta.ang );
zeta.cos = cos( zeta.ang );

alpha.ang = alpha.ang - pma*alpha.sin*(1./zeta.cos + tan(zeta.ang)) -
            rma*(1 - alpha.cos/zeta.cos);
zeta.ang = zeta.ang + rma*alpha.sin;
/*
  Corrected scan angle
*/
alpha.sin = sin( alpha.ang );
alpha.cos = cos( alpha.ang );
zeta.sin = sin( zeta.ang );
zeta.cos = cos( zeta.ang );
/*
  Computes pointing vector in instrument coordinates
*/
g.x =  zeta.sin;
g.y = -zeta.cos * alpha.sin;
g.z =  zeta.cos * alpha.cos;
/*
  Transforms the pointing vector to earth fixed coordinates
*/
g1.x = bt[0][0]*g.x + bt[0][1]*g.y + bt[0][2]*g.z;
g1.y = bt[1][0]*g.x + bt[1][1]*g.y + bt[1][2]*g.z;
g1.z = bt[2][0]*g.x + bt[2][1]*g.y + bt[2][2]*g.z;
/*
  Computes coefficients and solves a quadratic equation to
  find the intersect of the pointing vector with the earth
  surface
*/
q1 = g1.x*g1.x + g1.y*g1.y + REQRPL2*g1.z*g1.z;
q2 = xs.x*g1.x + xs.y*g1.y + REQRPL2*xs.z*g1.z;
d  = q2 * q2 - q1 * q3;

if( d < 1e-9 && d > -1e-9 )
   d = 0;
/*
  If the discriminant of the equation, d, is negative, the
  instrument points off the earth
*/
if( d < 0 ){
   *lat = MISS;
   *lon = MISS;
   return 1;
   }
d = sqrt(d);
/*
  Slant distance from the satellite to the earth point
*/
h = -(q2 + d) / q1;
/*
  Cartesian coordinates of the earth point
*/
u.x = xs.x + h*g1.x;
u.y = xs.y + h*g1.y;
u.z = xs.z + h*g1.z;
/*
  Sinus of geocentric latitude
*/
d1 = u.z / sqrt(u.x*u.x + u.y*u.y + u.z*u.z);
/*
  Geographic (geodetic) coordinates of the point
*/
*lat = atan( REQRPL2 * d1 / sqrt(1. - d1 * d1));
*lon = atan2( u.y, u.x );
return 0;
}

/************************************************************
   sat2earth (NVXSAE) - Converts line, elem on GVAR image to lat,lon
************************************************************/
int sat2earth( line, elem, lat, lon)

float line, elem;
float *lat, *lon;

{
float elev, scan;
/*
  If doing sounder nav, have to trick routines into thinking image is
  at res 1, because nav routines take sounder res into account
*/
if( instr == 1 ){
   line = (line+9)/10;
   line = (line+9)/10;
   }
/*
  Compute elevation and scan angles related to input
  line and pixel numbers
*/
elev = elev_line(instr,line);
scan = scan_angle(instr,elem);
/*
  Transform elevation and scan angles to geographic coordinates
*/
if( lpoint(elev,scan,lat,lon) == 0 ){
   *lat = *lat * RDC;
   *lon = *lon * RDC;
   return 0;
   }
return -1;
}

/************************************************************
   gpoint - Converts GVAR lat and lon to elevation and scan angle
************************************************************/
int gpoint(lat,lon,elev,scan)

float lat;     /* Geographic latitude in radians (input) */
float lon;     /* Geographic longitude in radians (input) */
float *elev;   /* Elevation angle in radians (output) */
float *scan;   /* Scan angle in radians (output) */

{
/*
  Output status; 0 - successful completion,
                 1 - point with given lat/lon is invisible
*/
struct pt u,f,ft;
double sing;
double slat;
double w1,w2;
/*
  Computes sinus of geographic (geodetic) latitude
*/
sing = sin(lat);
w1 = REQRPL4 * sing * sing;
/*
  Sinus of the geocentric latitude
*/
slat = ((0.375 * w1 - 0.5) * w1 + 1.) * sing / REQRPL2;
/*
  Computes local earth radius at specified point
*/
w2 = slat * slat;
w1 = REQRPL3 * w2;
w1 = (0.375 * w1 - 0.5) * w1 + 1;
/*
  Computes cartesian coordinates of the point
*/
u.z = slat * w1;
w2  = w1 * sqrt(1 - w2);
u.x = w2 * cos(lon);
u.y = w2 * sin(lon);
/*
  Pointing vector from satellite to the earth point
*/
f.x = u.x - xs.x;
f.y = u.y - xs.y;
f.z = u.z - xs.z;
w2  = u.x*f.x + u.y*f.y + u.z*f.z*REQRPL2;
/*
  Verifies visibility of the point
*/
if( w2 > 0 ){
   *elev  = MISS;
   *scan  = MISS;
   return 1;
   }
/*
  Converts pointing vector to instrument coordinates
*/
ft.x = bt[0][0]*f.x + bt[1][0]*f.y + bt[2][0]*f.z;
ft.y = bt[0][1]*f.x + bt[1][1]*f.y + bt[2][1]*f.z;
ft.z = bt[0][2]*f.x + bt[1][2]*f.y + bt[2][2]*f.z;
/*
  Converts pointing vector to scan and elevation angles and
  corrects for the roll and pitch misalignments
*/
*scan  =  atan( ft.x / sqrt(ft.y*ft.y + ft.z*ft.z));
*elev  = -atan( ft.y / ft.z );
w1    = sin(*elev);
w2    = cos(*scan);

*elev  = *elev + rma * (1 - cos(*elev) / w2) + pma * w1 * (1/w2 + tan(*scan));
*scan  = *scan - rma * w1;
return 0;
}

/************************************************************
  earth2sat (NVXEAS) - Converts lat, lon to line, elem on GVAR image
************************************************************/
void earth2sat( lat,lon,line,elem )

float lat,lon,*line,*elem;

{
float elev, scan;

if( lat > 90 || lat < -90 ) {
   *line=MISS;
   *elem=MISS;
   return;
  }
/*
  Transform lat/lon to elevation and scan angles
*/
if( gpoint(lat*DRC,lon*DRC,&elev,&scan)) {
   *line=MISS;
   *elem=MISS;
   return;
   }
/*
  Convert elevation and scan angles to line/pixel coordinates
*/
elev_scan2line_pix(elev,scan,line,elem);
/*
  If doing sounder nav, change lin & elem return to res 10 values
*/
if( instr == 1 ){
   *line = *line*10-9;
   *elem = *elem*10-9;
   }
}

#define GRACON 0.07436574  /* Gravitation constant */
#define J2 1082.28E-6      /* 2nd harmonic coefficient of earths 
                              aspherical gravitational potential */
#define J4 -2.12E-6        /* 4th harmonic coefficient of earths 
                              aspherical gravitational potential */
#define EPSILON 1.0E-8     /* Convergence criterion used for calculating 
                              eccentric anomaly */
#undef NUM_ITER
#define NUM_ITER 20

/************************************************************************
   SATPOS - Determines the location of the satellite based on orbit
   navigation information. (GOES only)
************************************************************************/
void satpos( sattim, coord, sat )

struct datim *sattim;
int coord;
struct pt *sat;

{
int i;
double temp;
long itemp;
static int initial = 0;
static long key = 0;
int year;
int day;
int iday;
int sign;
/*
   GEOS declarations
*/
static struct datim epoch;      /* Epoch time */
static struct datim peri_focus; /* Peri-focus time */
static struct angle inc;        /* Inclination angle */
static struct angle earth;      /* Earth position angle */
static struct angle per;        /* Orbit perigee position */
static struct angle asn;        /* Orbit ascending mode */
static struct pt p,q;
static struct pt tpt;           /* Temporary point */
static double ecc_anom;
static double ecc_fac;          /* eccentricity factor */
static double orb_sparam;       /* orbital semi-parameter */
static double mean_motion;      /* Mean motion constant */
static double an_mean_motion;   /* anamalistic mean motion constant */
static double mean_anom;        /* Mean anomaly */
static double del_per;
static double del_asn;
static double time_diff;
static double adj_perigee;
static double adj_ascnode;

/*
  Test to see if day or satellite has changed necessitating param update
*/
if( !initial ){
   if( g_orb.sat_type < 0 )
      g_orb.sat_type = -g_orb.sat_type;
/*
   Convert epoch to julian day - time 
*/
   epoch.date = mdcon(1, g_orb.epoch_date);
   epoch.tim = int2deg(g_orb.epoch_time);
/*
   Define mean anomaly

   Explicit relationships between v,e, and m are given by the following:

   cos(v)=(cos(e)-i)/(1-i*cos(e)
   sin(v)=sqrt(1-i*i)*sin(e)/(1-i*cos(e))
   cos(e)=(cos(v)+i)/(1+i*cos(v))
   sin(e)=(sqrt(1-i*i)*sin(v)/(1+i*cos(v))
   m=e-i*sin(e)   
*/
   if( g_orb.anom_type == 0 )
      mean_anom = g_orb.anomaly;
   else {
      temp = cos(DRC * g_orb.anomaly);
      ecc_anom = acos((temp + g_orb.eccent) / (1.0 + g_orb.eccent * temp));
      mean_anom = (ecc_anom - g_orb.eccent * sin(ecc_anom)) * RDC;
      }
   mean_anom = fmod(mean_anom,360.); 
   if( mean_anom < 0. ) mean_anom += 360.0;
/*
   Define eccentricity factor and orbital semi-parameter
*/
   ecc_fac = sqrt(1.0 - (g_orb.eccent * g_orb.eccent));
   orb_sparam = (g_orb.semi_maj / RADE_EQUA) * ecc_fac * ecc_fac;
/*
   Calculate inclination sin and cos terms
*/
   inc.ang = g_orb.inclin*DRC;
   inc.sin = sin(inc.ang);
   inc.cos = cos(inc.ang);
/*
   Mean motion constant - speed of tranversal around earth in rad/min
*/
   mean_motion = GRACON * pow((RADE_EQUA / g_orb.semi_maj), 1.5);
/*
   Calculate orbital period
*/
   g_orb.period = TWO_PI / mean_motion;
/*
   Calculate anomalistic mean motion constant and derivitives based
   on selected order of secular perturbation theory

   Zero order
*/
   if( g_orb.sec_ord == 0 ){
      an_mean_motion = mean_motion;
      del_per = 0.0;
      del_asn = 0.0;
      }
/*
   First order 
*/
   else if( g_orb.sec_ord == 1 ){
      an_mean_motion = mean_motion*(1.0 + (1.5*J2*ecc_fac/
         (orb_sparam*orb_sparam))*(1.0 - 1.5*inc.sin*inc.sin));
      del_per = (1.5*J2*(2.0 - 2.5*inc.sin*inc.sin)/
         (orb_sparam*orb_sparam))*an_mean_motion*RDC;
      del_asn = -(1.5*J2*inc.cos/(orb_sparam*orb_sparam))*an_mean_motion*RDC;
      }
/*
   Second order 

   else if( g_orb.sec_ord == 2 ){
      an_mean_motion = mean_motion*(1.0 + (1.5*J2*ecc_fac/
             (orb_sparam*orb_sparam))*(1.0 - 1.5*inc.sin*inc.sin) +
             (0.0234375*J2*J2*ecc_fac/pow(orb_sparam,4.))*(16.*ecc_fac +
             25.*ecc_fac*ecc_fac - 15. + (30. - 96.*ecc_fac - 
             90.*ecc_fac*ecc_fac)*inc.cos*inc.cos + (105. + 144.*ecc_fac + 
             25.*ecc_fac*ecc_fac)*pow(inc.cos,4.)) -
             (.3515625*J4*ecc_fac*g_orb.eccent*g_orb.eccent/pow(orb_sparam,4.))*
             (3. - 30.*inc.cos*inc.cos + 35.*pow(inc.cos,4.)));
      del_per = ((1.5*J2*an_mean_motion/(orb_sparam*orb_sparam))*(2. - 
             2.5*inc.sin*inc.sin)*(1. + (1.5*J2/(orb_sparam*orb_sparam))*
             (2. + g_orb.eccent*g_orb.eccent/2. - 2.*ecc_fac - (1.791666667 - 
             g_orb.eccent*g_orb.eccent/48. - 3.*ecc_fac)*inc.sin*inc.sin)) - 
             1.25*J2*J2*g_orb.eccent*g_orb.eccent*mean_motion*pow(inc.cos,4.)/
             pow(orb_sparam,4.) - (4.375*J4*mean_motion/pow(orb_sparam,4.))*
             (1.714285714 - 6.642857143*inc.sin*inc.sin + 
             5.25*pow(inc.sin,4.) + g_orb.eccent*g_orb.eccent*(1.928571429 - 
             6.75*inc.sin*inc.sin + 5.0625*pow(inc.sin,4.))))*RDC;
      del_asn = -((1.25*J2*an_mean_motion*inc.cos/(orb_sparam*orb_sparam))*
             (1. + (1.5*J2/(orb_sparam*orb_sparam))*(1.5 + 
             g_orb.eccent*g_orb.eccent/6. - 2.*ecc_fac - (1.666666667 - 
             .208333333*g_orb.eccent*g_orb.eccent - 3.*ecc_fac)*
             inc.sin*inc.sin)) + (4.375*J4*mean_motion/pow(orb_sparam,4.))*
             (1. + 1.5*g_orb.eccent*g_orb.eccent)*(.857142857 - 
             1.5*inc.sin*inc.sin)*inc.cos)*RDC;
      }
*/
/*
   Calculate anomalistic period 
*/
   g_orb.anom_per = TWO_PI / an_mean_motion;
/*
   Determine time of peri-focal passage 
*/
   year = epoch.date / 1000;
   day = epoch.date % 1000;
   peri_focus.tim = epoch.tim - DRC * mean_anom / (60.0 * an_mean_motion);
   if( peri_focus.tim >= 0.0 )
      sign = 1;
   else
      sign = -1;
   iday = sign * (fabs(peri_focus.tim) / 24.0 + 1.0);
   if( iday > 0 ) iday--;
   peri_focus.tim -= iday * 24.0;
   if( iday != 0 ){
      day += iday;
      if( day < 1 ){
         year--;
         day += numyr((long)year);
         }
      else if( day > numyr((long)year)){
         day -= numyr((long)year);
         year++;
         }
      }
   itemp = deg2int(peri_focus.tim);
   peri_focus.tim = int2deg(itemp);
   peri_focus.date = 1000L * year + day;
/*
   Adjust perigee and ascending node to time of peri-focal passage
*/
   time_diff = timdif(epoch.date, epoch.tim, peri_focus.date, peri_focus.tim);
   adj_perigee = g_orb.perigee + del_per * time_diff;
   adj_perigee = fmod(adj_perigee,360.);
   if( adj_perigee < 0. ) adj_perigee += 360.;
   adj_ascnode = g_orb.asc_node + del_asn * time_diff;
   adj_ascnode = fmod(adj_ascnode,360.);
   if( adj_ascnode < 0. ) adj_ascnode += 360.;
   key = 1;
   }
/*
   Calculate delta-time ( from time of peri-focus to specified time )
*/
time_diff = timdif(peri_focus.date, peri_focus.tim, sattim->date, sattim->tim);
if( g_orb.sec_ord != 0 || key != 0 ){
   key = 0;
/*
   Calculate time dependent values of perigee and ascending node
*/
   per.ang = DRC * (adj_perigee + del_per * time_diff);
   asn.ang = DRC * (adj_ascnode + del_asn * time_diff);
/*
   Calculate perigee and ascending node sin and cos terms
*/
   per.sin = sin(per.ang);
   per.cos = cos(per.ang);
   asn.sin = sin(asn.ang);
   asn.cos = cos(asn.ang);
/*
   Calculate the (P,Q,W) orthogonal orientation vectors P
   points toward peri-focus Q is in the orbit plane advanced
   from P by a right angle in the direction of increasing
   true anomaly W completes a right handed coordinate system
*/
   p.x =  per.cos*asn.cos - per.sin*asn.sin*inc.cos;
   p.y =  per.cos*asn.sin + per.sin*asn.cos*inc.cos;
   p.z =  per.sin*inc.sin;
   q.x = -per.sin*asn.cos - per.cos*asn.sin*inc.cos;
   q.y = -per.sin*asn.sin + per.cos*asn.cos*inc.cos;
   q.z =  per.cos*inc.sin;
/*
   The W terms were commented out in the original
   wx =  asn.sin*inc.sin;
   wy = -asn.cos*inc.sin;
   wz =  inc.cos;
*/
   initial = 1;
   }
/*
   define mean anomaly (m) at specified time 
*/
mean_anom = fmod((an_mean_motion*time_diff),TWO_PI);
/*
   Calculate eccentric anomaly (e) at specified time

   The solution is given by a simplified numerical (Newtons) method. An
   explicit relationship involves a bessel function of the first
   kind.
 
   e = m+2*sum(n=1,infinity)(j(n)(n*i)*sin(n*m))
*/
temp = mean_anom;
for( i = 0; i < NUM_ITER; i++ ){
   ecc_anom = mean_anom + g_orb.eccent * sin(temp);
   if( fabs(ecc_anom - temp) < EPSILON ) break;
   temp = ecc_anom;
   }
/*
   Expression for magnitude of satellite radius vector (r)

   r = re*orb_sparam/(1.0+g_orb.eccent*cos(ecc_anom))

   Generate a position vector with respect to the focus and in the
   orbital plane. Note that the Z coordinate is by definition zero.
*/
tpt.x = g_orb.semi_maj * (cos(ecc_anom) - g_orb.eccent);
tpt.y = g_orb.semi_maj * (sin(ecc_anom) * ecc_fac);
/*
   Transformation to a celestial pointing vector by utilization of
   the transpose of the (P,Q,W) orthoganal transfromation matrix.
   Note that the third row containing W is not required because
   zomega is zero.
*/
sat->x = tpt.x * p.x + tpt.y * q.x;
sat->y = tpt.x * p.y + tpt.y * q.y;
sat->z = tpt.x * p.z + tpt.y * q.z;
if( coord == 0 ){
/*
   Determine transformation matrix for rotation to terrestrial coord's
*/
   time_diff = timdif( ref.date, ref.tim, sattim->date, sattim->tim);
   earth.ang = time_diff * (earth_rot - prec_veq);
   earth.ang = fmod(earth.ang,TWO_PI);
   if( earth.ang < 0 ) earth.ang += TWO_PI;
   earth.sin = sin(earth.ang);
   earth.cos = cos(earth.ang);
   tpt = *sat;
/*
   Rotation to terrestrial pointing vector 
*/
   sat->x = earth.cos * tpt.x + earth.sin * tpt.y;
   sat->y = -earth.sin * tpt.x + earth.cos * tpt.y;
   }
return;
}

/************************************************************************
   SAT_TRAN - Performs coordinate transformations between earth and
   satellite coordinate systems based on orbit navigation information
   for the satellite.  This is based on the Wisconsin SSEC/McIDAS area
   navigation information.

   type = transformation type
        = 1 for satellite to earth transformation (specify xlin,xele)
        = 2 for earth to satellite transformation (specify xlat,xlon)
        = 3 for subpoint calculation (specify nothing)
************************************************************************/
#undef NUM_ITER
#define NUM_ITER 5

void sat_tran( type,lat,lon,lin,ele )

int type;
float *lat,*lon,*lin,*ele;

{
double temp;
static int initial = 0;
static double center_elem;       /* Picture center element */
static double rad_per_line;      /* Radians per scan line */
static double rad_per_elem;      /* Radians per element */
static double time_diff;         /* Time difference in minutes */
static double rot_ax[3][2];      /* Rotation axis transformation matrix */
static double rot_sat[3][3];     /* Satellite axis transformation matrix */
static int num_sensor;           /* Number of sensors on satellite */
static struct datim pix;         /* Pixel sample point time */
static struct angle line;        /* */
static struct angle elem;
static struct angle lats;
static struct angle lons;
static struct angle ras;
static struct angle dec;
static struct angle roll;
static struct angle pitch;
static struct angle yaw;
static struct angle prec;        /* Precession direction */
static struct angle ptot;        /* Total precession to sample point time */
static struct angle gamma;       /* Image east-west drift angle */
static struct angle earth;       /* Earth rotation */
static struct pt spin_ax;        /* Spin axis in earth ref frame at t(0) */
static struct pt u_sat;
static struct pt v_sat;
static struct pt pts[3];
static struct pt pt_vec,tmp_vec;
static long iter;
static double angx,angy;
static double x,y;
static double par_line,par_elem; /* Partial line,element counts */
static double prev_scan,new_scan;

if( image_type == GOES ){
   if( !initial ){
/*
   Expanded set of transformation parameters 
*/
      num_sensor = (g_nav.num_lines/100000L)%100;
      if( num_sensor < 1 ) num_sensor = 1;
      g_nav.num_lines = num_sensor*(g_nav.num_lines%100000L); 
      center_elem = (1 + g_nav.num_elem)/2.;
      rad_per_line = g_nav.deg_line*DRC/(g_nav.num_lines - 1.);
      rad_per_elem = g_nav.deg_elem*DRC/(g_nav.num_elem - 1.);

      dec.ang = g_nav.declin*DRC;
      dec.sin = sin(dec.ang);
      dec.cos = cos(dec.ang);

      time_diff = timdif(ref.date,ref.tim,sat_time.date,0.);
      ras.ang = g_nav.right_asc*DRC - time_diff*(earth_rot - prec_veq);
      ras.ang = fmod(ras.ang,TWO_PI);
      if( ras.ang < 0 ) ras.ang += TWO_PI;
      ras.sin = sin(ras.ang);
      ras.cos = cos(ras.ang);

      spin_ax.x = dec.cos*ras.cos;
      spin_ax.y = dec.cos*ras.sin;
      spin_ax.z = dec.sin;

      g_nav.prec_rate = g_nav.prec_rate*DRC/24.;
      prec.ang = g_nav.prec_dir*DRC;
      prec.sin = sin(prec.ang);
      prec.cos = cos(prec.ang);

      pitch.ang = g_nav.pitch*DRC;
      pitch.sin = sin(pitch.ang); 
      pitch.cos = cos(pitch.ang);

      yaw.ang = g_nav.yaw*DRC;
      yaw.sin = sin(yaw.ang);
      yaw.cos = cos(yaw.ang);

      roll.ang = g_nav.roll*DRC;
      roll.sin = sin(roll.ang);
      roll.cos = cos(roll.ang);
/*
   Compute rotational matrix for pinciple axis misalignment 
*/
      rot_ax[0][0] = roll.cos*pitch.cos;
      rot_ax[0][1] = yaw.sin*roll.sin*pitch.cos + yaw.cos*pitch.sin;
      rot_ax[1][0] = -roll.sin;
      rot_ax[1][1] = yaw.sin*roll.cos;
      rot_ax[2][0] = -roll.cos*pitch.sin;
      rot_ax[2][1] = yaw.cos*pitch.cos - yaw.sin*roll.sin*pitch.sin;
      initial = 1;
      }
/*
   Normalize satellite coordinates and convert earth coordinates to radians
   transform latitude to geocentric coordinates check for transform direction 
*/
   switch( type ){
      case 1:
         line.ang = rad_per_line*(*lin - g_nav.center_line);
         elem.ang = rad_per_elem*(*ele - center_elem);
         line.sin = sin(line.ang);
         line.cos = cos(line.ang);
         elem.sin = sin(elem.ang);
         elem.cos = cos(elem.ang);
         iter = NUM_ITER;
         break;
      case 2:
         *lin = g_nav.center_line;
         *ele = center_elem;
         lats.ang = DRC*geolat(1,*lat);
         lons.ang = DRC*(*lon);
         lats.sin = sin(lats.ang);
         lats.cos = cos(lats.ang);
         lons.sin = sin(lons.ang);
         lons.cos = cos(lons.ang);
         iter = 0; 
         break;
      case 3: 
         *lin = g_nav.center_line; 
         *ele = center_elem;
         line.ang = rad_per_line*(*lin - g_nav.center_line);
         elem.ang = rad_per_elem*(*ele - center_elem);
         line.sin = sin(line.ang);
         line.cos = cos(line.ang);
         elem.sin = sin(elem.ang);
         elem.cos = cos(elem.ang);
         iter = 0; 
         break;
      }
   do {
/*
   Compute time dependent parameters 
*/
      par_line = (double)((round_long(*lin) - 1)/num_sensor); 
      par_elem = (*ele - 1.)/(g_nav.num_elem - 1.)*(g_nav.deg_elem/360.);
      pix.date = sat_time.date;
      pix.tim = sat_time.tim + g_nav.spin_rate*(par_line + par_elem);
      ptot.ang = pix.tim*g_nav.prec_rate;
      ptot.sin = sin(ptot.ang);
      ptot.cos = cos(ptot.ang);
      gamma.ang = rad_per_elem*(gamma_init + gamma_dot*pix.tim);
      gamma.sin = sin(gamma.ang);
      gamma.cos = cos(gamma.ang);
      earth.ang = PI*SOLSID/12.*pix.tim;
      earth.sin = sin(earth.ang);
      earth.cos = cos(earth.ang);
/*
   Determine satellite position vector (keplerian model) 
*/
      satpos(&pix, 0, &sat);
      sphere_cvt( 1, &sat, &sp_sat );
/*
   Compute unit pointing vector 
*/
      u_sat.x = sat.x/sp_sat.hgt;
      u_sat.y = sat.y/sp_sat.hgt;
      u_sat.z = sat.z/sp_sat.hgt;
/*
   Repoint spin axis as function of precession 
*/
      pts[0].x = sqrt(1./(1. + (spin_ax.x/spin_ax.z)*(spin_ax.x/spin_ax.z)));
      pts[0].y = 0.;
      pts[0].z = -(spin_ax.x*pts[0].x/spin_ax.z);
      pts[1].x = spin_ax.y*pts[0].z;
      pts[1].y = spin_ax.z*pts[0].x - spin_ax.x*pts[0].z;
      pts[1].z = -spin_ax.y*pts[0].x;
      pts[0].x = prec.cos*pts[0].x + prec.sin*pts[1].x;
      pts[0].y = prec.cos*pts[0].y + prec.sin*pts[1].y;
      pts[0].z = prec.cos*pts[0].z + prec.sin*pts[1].z;
      spin_ax.x = ptot.cos*spin_ax.x + ptot.sin*pts[0].x;
      spin_ax.y = ptot.cos*spin_ax.y + ptot.sin*pts[0].y;
      spin_ax.z = ptot.cos*spin_ax.z + ptot.sin*pts[0].z;
/*
   Compute nominal satellite position rotational matrix 
*/
      rot_sat[2][0] = earth.cos*spin_ax.x + earth.sin*spin_ax.y;
      rot_sat[2][1] = -earth.sin*spin_ax.x + earth.cos*spin_ax.y;
      rot_sat[2][2] = spin_ax.z;
      temp = u_sat.x*rot_sat[2][0] +
             u_sat.y*rot_sat[2][1] +
             u_sat.z*rot_sat[2][2];
      v_sat.x = u_sat.x - temp*rot_sat[2][0];
      v_sat.y = u_sat.y - temp*rot_sat[2][1];
      v_sat.z = u_sat.z - temp*rot_sat[2][2];
/*
   Normalize v_sat
*/
      temp = -1./sqrt(v_sat.x*v_sat.x + v_sat.y*v_sat.y + v_sat.z*v_sat.z);
      rot_sat[0][0] = v_sat.x*temp;
      rot_sat[0][1] = v_sat.y*temp;
      rot_sat[0][2] = v_sat.z*temp;
      rot_sat[1][0] = rot_sat[2][1]*rot_sat[0][2] -
                      rot_sat[2][2]*rot_sat[0][1];
      rot_sat[1][1] = rot_sat[2][2]*rot_sat[0][0] -
                      rot_sat[2][0]*rot_sat[0][2];
      rot_sat[1][2] = rot_sat[2][0]*rot_sat[0][1] -
                      rot_sat[2][1]*rot_sat[0][0];
      rot_sat[0][0] = rot_sat[0][0]*gamma.cos - gamma.sin*rot_sat[1][0];
      rot_sat[1][0] = rot_sat[0][0]*gamma.sin + gamma.cos*rot_sat[1][0];
      rot_sat[0][1] = rot_sat[0][1]*gamma.cos - gamma.sin*rot_sat[1][1];
      rot_sat[1][1] = rot_sat[0][1]*gamma.sin + gamma.cos*rot_sat[1][1];
      rot_sat[0][2] = rot_sat[0][2]*gamma.cos - gamma.sin*rot_sat[1][2];
      rot_sat[1][2] = rot_sat[0][2]*gamma.sin + gamma.cos*rot_sat[1][2];
/*
   Switch on transformation type 
*/
      switch( type ){
/*
   Case 1 - Transform from satellite coordinates to earth coordinates 
*/
         case 1: 
            {
            double a,b,c,rad;
            double equa_sq, pole_sq, ratio, in_ratio;
/*
   Compute pointing vector in satellite coordinate system at element 0 
*/
            pt_vec.x = rot_ax[0][0]*line.cos - rot_ax[0][1]*line.sin;
            pt_vec.y = rot_ax[1][0]*line.cos - rot_ax[1][1]*line.sin;
            pt_vec.z = rot_ax[2][0]*line.cos - rot_ax[2][1]*line.sin;
/*
   Adjust pointing vector for element count 
*/
            tmp_vec = pt_vec;
            pt_vec.x =  elem.cos*tmp_vec.x + elem.sin*tmp_vec.y;
            pt_vec.y = -elem.sin*tmp_vec.x + elem.cos*tmp_vec.y;
            pt_vec.z =  tmp_vec.z;
/*
   Compute pointing vector in earth coordinate system 
*/
            tmp_vec = pt_vec;
            pt_vec.x = rot_sat[0][0]*tmp_vec.x +
                       rot_sat[1][0]*tmp_vec.y +
                       rot_sat[2][0]*tmp_vec.z;
            pt_vec.y = rot_sat[0][1]*tmp_vec.x +
                       rot_sat[1][1]*tmp_vec.y +
                       rot_sat[2][1]*tmp_vec.z;
            pt_vec.z = rot_sat[0][2]*tmp_vec.x +
                       rot_sat[1][2]*tmp_vec.y +
                       rot_sat[2][2]*tmp_vec.z;
/*
   Adjust for oblateness of earth sphere and atmospheric height 
*/
            equa_sq = (RADE_EQUA + atms_hght)*(RADE_EQUA + atms_hght);
            pole_sq = (RADE_POLE + atms_hght)*(RADE_POLE + atms_hght);
            ratio = pole_sq/equa_sq;
            in_ratio = 1. - ratio;
            a = ratio + in_ratio*pt_vec.z*pt_vec.z;
            b = 2.*((pt_vec.x*sat.x + pt_vec.y*sat.y)*ratio + pt_vec.z*sat.z);
            c = (sat.x*sat.x + sat.y*sat.y)*ratio + sat.z*sat.z - pole_sq;
/*
   Compute radical 
*/
            rad = b*b - 4.*a*c;
/*
   Check if point is off earth and if so set rejection values 
*/
            if( rad >= 1. ){
/*
   Find point along pointing vector intersecting earth surface 
*/
               rad = -(b + sqrt(rad))/(2.*a);
               sat.x = sat.x + pt_vec.x*rad;
               sat.y = sat.y + pt_vec.y*rad;
               sat.z = sat.z + pt_vec.z*rad;
/*
   Convert to earth coordinates 
*/
               lats.ang = atan2(sat.z,sqrt(sat.x*sat.x + sat.y*sat.y));
               lons.ang = atan2(sat.y,sat.x);
/*
   Convert from radians to degrees and convert latitude to geodedtic coor 
*/
               *lat = geolat(2,lats.ang*RDC);
               *lon = lons.ang*RDC;
               return;
               }
            *lat = MISS;
            *lon = MISS;
            return;
            }
         case 2:
/*
   Transform from earth coordinates to satellite coordinates 
*/
            iter++;
            if( iter == 1 ){
/*
   Compute earth coordinate vector 
*/
               pt_vec.x = lats.cos*lons.cos;
               pt_vec.y = lats.cos*lons.sin;
               pt_vec.z = lats.sin; 
/*
   Check if point is out of saatellite view and 
   if so set rejection value and return 
*/
               temp = RADE_MEAN/sp_sat.hgt;
               if( pt_vec.x*u_sat.x + pt_vec.y*u_sat.y + pt_vec.z*u_sat.z < temp ){
                  *lin = MISS;
                  *ele = MISS;
                  return;
                  }
/*
   Adjust for oblateness of earth sphere and atmospheric height 
*/
               temp = (lats.sin/lats.cos)*(lats.sin/lats.cos);
               temp = sqrt((1. + temp)/
                      (RADE_POLE*RADE_POLE + RADE_EQUA*RADE_EQUA*temp))*
                      RADE_EQUA*RADE_POLE + atms_hght;
               pt_vec.x = temp*pt_vec.x;
               pt_vec.y = temp*pt_vec.y;
               pt_vec.z = temp*pt_vec.z;
               }
            break;
         case 3:
            iter++;
            if( iter == 1 ){
               pt_vec.x = 0.;
               pt_vec.y = 0.;
               pt_vec.z = 0.;
               }
            break;
         }
/*
   Compute vector from satellite to earth in earth coordinate system 
*/
      pts[0].y = pt_vec.x - sat.x;
      pts[1].y = pt_vec.y - sat.y;
      pts[2].y = pt_vec.z - sat.z;
/*
   Compute normalized pointing vector in earth coordinate system 
*/
      temp = 1./sqrt(pts[0].y*pts[0].y + pts[1].y*pts[1].y + pts[2].y*pts[2].y);
      pts[0].y = pts[0].y*temp;
      pts[1].y = pts[1].y*temp;
      pts[2].y = pts[2].y*temp;
/*
   Compute vector from satellite to earth in earth coordinate system 
*/
      pts[0].x = rot_sat[0][0]*pts[0].y +
                 rot_sat[0][1]*pts[1].y +
                 rot_sat[0][2]*pts[2].y;
      pts[1].x = rot_sat[1][0]*pts[0].y +
                 rot_sat[1][1]*pts[1].y +
                 rot_sat[1][2]*pts[2].y;
      pts[2].x = rot_sat[2][0]*pts[0].y +
                 rot_sat[2][1]*pts[1].y +
                 rot_sat[2][2]*pts[2].y;
/*
   Convert pointing vector to line-element coordinate 
*/
      temp = sqrt(rot_ax[2][0]*rot_ax[2][0] + rot_ax[2][1]*rot_ax[2][1] -
                  pts[2].x*pts[2].x);
      angx = atan2(pts[2].x,temp) - atan2(rot_ax[2][0],rot_ax[2][1]);
      x = rot_ax[0][0]*cos(angx) + rot_ax[0][1]*sin(angx);
      y = rot_ax[1][0]*cos(angx) + rot_ax[1][1]*sin(angx);
      angy = atan2(pts[1].x,pts[0].x) - atan2(y,x);
      prev_scan = (round_long(*lin) - 1)/num_sensor;
      *lin = g_nav.center_line - angx/rad_per_line;
      *ele = center_elem - angy/rad_per_elem;
      new_scan = (round_long(*lin) - 1)/num_sensor;
      }
   while( iter < NUM_ITER && prev_scan != new_scan );
   switch( type ){
      case 2:
/*
   Check if point is off frame and if so set rejection values 
*/
         if( *lin < 1 || *lin > g_nav.num_lines ||
             *ele < 1 || *ele > g_nav.num_elem ){
            *lin = MISS;
            *ele = MISS;
            }
         return;
      case 3:
/*
   Compute sub-satellite point from satellite position vector 
*/
         *lat = sp_sat.lat;
         *lon = sp_sat.lon;
         return;
      }
   }
/************************************************************
   GVAR Code
************************************************************/
else if( image_type == GVAR ){
/*
   Switch on transformation type 
*/
   switch( type ){
/*
   Case 1 - Transform from satellite coordinates to earth coordinates 
*/
      case 1: 
         sat2earth( *lin, *ele, lat, lon);
         break;
/*
   Transform from earth coordinates to satellite coordinates
*/
      case 2:
         earth2sat( *lat, *lon, lin, ele );
         break;
/*
   Compute sub-satellite point from satellite position vector 
*/
      case 3:
         *lat = sublat*RDC;
         *lon = sublon*RDC;
         return;
      }
   }
}
/************************************************************************
   SAT_INIT - Initializes the satellite parameters.
************************************************************************/
int sat_init( area )

struct mc_area *area;

{
int i;

if( !strncmp( area->type,"GOES",4 )){
/*
   Initialize the gv_nav structure from the data read in from image file
*/
   printf(" Initializing GOES data navigation\n");
   image_type = GOES;
   sat_image       = area;
   sat_time.date   = area->navg->iddate%100000;
   sat_time.tim    = int2deg((long)area->navg->itime );
   g_orb.sat_type    = 1;
   g_orb.sec_ord     = 0;
   g_orb.anom_type   = 0;
   g_orb.epoch_date  = area->navg->edate;
   g_orb.epoch_time  = area->navg->etime;
   g_orb.semi_maj    = area->navg->semima / 100.;
   g_orb.eccent      = area->navg->eccen / 1E6;
   g_orb.inclin      = area->navg->orbinc / 1000.;
   g_orb.anomaly     = area->navg->meana / 1000.;
   g_orb.perigee     = area->navg->perigee / 1000.;
   g_orb.asc_node    = area->navg->asnode / 1000.;
   g_nav.declin      = int2deg((long)area->navg->declin);
   g_nav.right_asc   = int2deg((long)area->navg->rascen);
   g_nav.center_line = area->navg->piclin;
   g_nav.spin_rate   = area->navg->spinp / 3.6E9;
   g_nav.deg_line    = int2deg((long)area->navg->deglin);
   g_nav.num_lines   = area->navg->lintot;
   g_nav.deg_elem    = int2deg((long)area->navg->degele);
   g_nav.num_elem    = area->navg->eletot;
   g_nav.pitch       = int2deg((long)area->navg->pitch);
   g_nav.yaw         = int2deg((long)area->navg->yaw);
   g_nav.roll        = int2deg((long)area->navg->roll);
   g_nav.prec_rate   = 0;
   g_nav.prec_dir    = 0;
   ecor = sat_image->dir->ecor;
   lcor = sat_image->dir->lcor;
   eres = sat_image->dir->eres;
   lres = sat_image->dir->lres;
   gamma_init      = area->navg->gamma / 100.;
   gamma_dot       = area->navg->gamdot / 100.;
   atms_hght       = 0;
   ref.date = REF_DATE;
   ref.tim = int2deg((long)REF_TIME);
/*
   Rotation rate of the vernal equinox in terms of sidereal time 
*/
   prec_veq = TWO_PI*SOLSID/(PREC_VER_EQ*SOLAR_YR*1440.);
   prec_veq = 0;
/*
   Terrestrial rotation rate in terms of sidereal time 
*/ 
   earth_rot = TWO_PI*SOLSID/1440.;

   printf( "date            = %ld\n", sat_time.date );
   printf( "pict_time       = %f\n", sat_time.tim );
   printf( "orb.sat_type    = %d\n", g_orb.sat_type );
   printf( "orb.sec_ord     = %d\n", g_orb.sec_ord );
   printf( "orb.anom_type   = %d\n", g_orb.anom_type );
   printf( "orb.epoch_date  = %ld\n", g_orb.epoch_date );
   printf( "orb.epoch_time  = %ld\n", g_orb.epoch_time );
   printf( "orb.semi_maj    = %f\n", g_orb.semi_maj );
   printf( "orb.eccent      = %f\n", g_orb.eccent );
   printf( "orb.inclin      = %f\n", g_orb.inclin );
   printf( "orb.anomaly     = %f\n", g_orb.anomaly );
   printf( "orb.perigee     = %f\n", g_orb.perigee );
   printf( "orb.asc_node    = %f\n", g_orb.asc_node );
   printf( "nav.declin      = %f\n", g_nav.declin );
   printf( "nav.right_asc   = %f\n", g_nav.right_asc );
   printf( "nav.center_line = %ld\n", g_nav.center_line );
   printf( "nav.spin_rate   = %f\n", g_nav.spin_rate );
   printf( "nav.deg_line    = %f\n", g_nav.deg_line );
   printf( "nav.num_lines   = %ld\n", g_nav.num_lines );
   printf( "nav.deg_elem    = %f\n", g_nav.deg_elem );
   printf( "nav.num_elem    = %ld\n", g_nav.num_elem );
   printf( "nav.pitch       = %f\n", g_nav.pitch );
   printf( "nav.yaw         = %f\n", g_nav.yaw );
   printf( "nav.roll        = %f\n", g_nav.roll );
   printf( "nav.prec_rate   = %f\n", g_nav.prec_rate );
   printf( "nav.prec_dir    = %f\n", g_nav.prec_dir );
   printf( "gamma_init      = %f\n", gamma_init );
   printf( "gamma_dot       = %f\n", gamma_dot );
   printf( "atms_hght       = %f\n", atms_hght );
   return 0;
   }
/************************************************************
   GVAR Code
************************************************************/
else if( !strncmp( area->type,"GVAR",4 )){
   int imc;
   double lam, dlat, dyaw;
   double r;
   double te, ts, wa;
   struct angle oi, u, asc;
/*
   Initialize the gv_nav structure from the data read in from image file
*/
   printf(" Initializing GVAR data navigation\n");
   image_type = GVAR;
   ecor = area->dir->ecor;
   lcor = area->dir->lcor;
   eres = area->dir->eres;
   lres = area->dir->lres;
   strcpy( area->navgv->sttype, gv_nav.sttype );
   gv_nav.idntfr = area->navgv->idntfr;
   gv_nav.imcact = area->navgv->imcact;
   gv_nav.reflon = area->navgv->reflon/10000000.;
   gv_nav.refdis = area->navgv->refdis/10000000.;
   gv_nav.reflat = area->navgv->reflat/10000000.;
   gv_nav.refyaw = area->navgv->refyaw/10000000.;
   gv_nav.ratrol = area->navgv->ratrol/10000000.;
   gv_nav.ratptc = area->navgv->ratptc/10000000.;
   gv_nav.ratyaw = area->navgv->ratyaw/10000000.;
   gv_nav.epoch = time50( area->navgv->etime );
   gv_nav.edtime = area->navgv->edtime/100.;
   gv_nav.imcrol = area->navgv->imcrol/10000000.;
   gv_nav.imcptc = area->navgv->imcptc/10000000.;
   gv_nav.imcyaw = area->navgv->imcyaw/10000000.;
   for( i = 0; i < 13; i++ )
      gv_nav.ldr[i] = area->navgv->ldr[i]/10000000.;
   for( i = 0; i < 11; i++ )
      gv_nav.rddr[i] = area->navgv->rddr[i]/10000000.;
   for( i = 0; i < 9; i++ )
      gv_nav.dgl[i] = area->navgv->dgl[i]/10000000.;
   for( i = 0; i < 9; i++ )
      gv_nav.doy[i] = area->navgv->doy[i]/10000000.;
   gv_nav.exptim = area->navgv->exptim/100.;
   gv_nav.solinc = area->navgv->solinc/10000000.;

   gatt_init( gv_nav.raawds, area->navgv->raawds );
   gatt_init( gv_nav.paawds, area->navgv->paawds );
   gatt_init( gv_nav.yaawds, area->navgv->yaawds );
   gatt_init( gv_nav.rmawds, area->navgv->rmawds );
   gatt_init( gv_nav.pmawds, area->navgv->pmawds );

   {
   int year, day, hour, min;
   float sec;
   
   year = area->navgv->imgday / 1000 + 1900;
   day = area->navgv->imgday % 1000;

   hour = area->navgv->imgtm / 10000000;
   min = area->navgv->imgtm % 10000000 / 100000;
   sec = area->navgv->imgtm % 100000 / 1000.;
   gv_nav.imgtim = timex( year, day, hour, min, sec );
   }

   instr = gv_nav.imgsnd = area->navgv->imgsnd-1;
/*
   Initialize constants 
*/
   set_gvar_con( gv_nav.imgsnd, area->navgv->iofnc, 
      area->navgv->iofni, area->navgv->iofec, area->navgv->iofei );
   if ( area->navgv->imcact & (1<<7))
      imc = 0;
   else
      imc = 1;
/*
  Assign reference values to the subsatellite longitude and
  latitude, the radial distance and the orbit yaw.
*/
   ts = 0;
   lam = gv_nav.reflon;
   dr  = gv_nav.refdis;
   phi.ang = gv_nav.reflat;
   psi.ang = gv_nav.refyaw;
/*
  Assign reference values to the attitudes and misalignments
*/
   roll  = gv_nav.ratrol;
   pitch = gv_nav.ratptc;
   yaw   = gv_nav.ratyaw;
   rma   = 0;
   pma   = 0;
/*
  If imc is off, compute changes in the satellite orbit
*/
   if( imc ){
      struct angle orb,orb1,orb2,orb3;
/*
  Set reference radial distance, latitude and orbit yaw to zero
*/
      dr  = 0;
      phi.ang = 0;
      psi.ang = 0;
/*
  Compute time since epoch (in minutes)
*/
      ts = gv_nav.imgtim - gv_nav.epoch;
/*
  Computes orbit angle and the related trigonometric functions.
  earth rotational rate=.729115e-4 (rad/s)
*/
      orb.ang   = (double)0.729115 - 240 * ts;
      orb.sin = sin( orb.ang );
      orb.cos = cos( orb.ang );
      orb1.ang = .927*orb.ang;
      orb1.sin = sin( orb1.ang );
      orb1.cos = cos( orb1.ang );
      orb2.ang = 2*orb.ang;
      orb2.sin = sin( orb2.ang );
      orb2.cos = cos( orb2.ang );
      orb3.ang = 1.9268*orb.ang;
      orb3.sin = sin( orb3.ang );
      orb3.cos = cos( orb3.ang );
/*
  Computes change in the imc longitude from the reference
*/
      lam = lam + gv_nav.ldr[0] + 
            (gv_nav.ldr[1] + gv_nav.ldr[2]*orb.ang)*orb.ang + 
            (gv_nav.ldr[9]*orb1.sin + gv_nav.ldr[10]*orb1.cos + 
            gv_nav.ldr[3]*orb.sin + gv_nav.ldr[4]*orb.cos + 
            gv_nav.ldr[5]*orb2.sin + gv_nav.ldr[6]*orb2.cos + 
            gv_nav.ldr[7]*orb3.sin + gv_nav.ldr[8]*orb3.cos + 
            orb.ang*(gv_nav.ldr[11]*orb.sin + gv_nav.ldr[12]*orb.cos))*2;
/*
  Computes change in radial distance from the reference (km)
*/
      dr = dr + gv_nav.rddr[0] + gv_nav.rddr[1]*orb.cos + 
           gv_nav.rddr[2]*orb.sin + gv_nav.rddr[3]*orb2.cos + 
           gv_nav.rddr[4]*orb2.sin + gv_nav.rddr[5]*orb3.cos +
           gv_nav.rddr[6]*orb3.sin + gv_nav.rddr[7]*orb1.cos + 
           gv_nav.rddr[8]*orb1.sin + 
           orb.ang*(gv_nav.rddr[9]*orb.cos + gv_nav.rddr[10]*orb.sin);
/*
  Computes the sine of the change in the geocentric latitude
*/
      dlat = gv_nav.dgl[0] + gv_nav.dgl[1]*orb.cos + gv_nav.dgl[2]*orb.sin + 
             gv_nav.dgl[3]*orb2.cos + gv_nav.dgl[4]*orb2.sin + 
             orb.ang*(gv_nav.dgl[5]*orb.cos + gv_nav.dgl[6]*orb.sin) + 
             gv_nav.dgl[7]*orb1.cos + gv_nav.dgl[8]*orb1.sin;
/*
  Computes geocentric latitude by using an expansion for arcsine
*/
      phi.ang = phi.ang + dlat * (1. + dlat * dlat / 6.);
/*
  Computes sine of the change in the orbit yaw
*/
      dyaw = gv_nav.doy[0] + gv_nav.doy[1]*orb.sin + gv_nav.doy[2]*orb.cos + 
             gv_nav.doy[3]*orb2.sin + gv_nav.doy[4]*orb2.cos + 
             orb.ang*(gv_nav.doy[5]*orb.sin + gv_nav.doy[6]*orb.cos) + 
             gv_nav.doy[7]*orb1.sin + gv_nav.doy[8]*orb1.cos;
/*
  Computes the orbit yaw by using an expansion for arcsine
*/
      psi.ang = psi.ang + dyaw * (1. + dyaw * dyaw / 6.);
      }
/*
  Conversion of the imc longitude and orbit yaw to the subsatellite
  longitude and the orbit inclination (Ref: GOES-PCC-TM-2473, inputs
  required for earth location and gridding by SPS, June 6, 1988)
*/
   phi.sin  = sin(phi.ang);
   psi.sin  = sin(psi.ang);
   oi.ang = phi.sin*phi.sin + psi.sin*psi.sin;
   oi.cos = sqrt(1-oi.ang);
   oi.sin = sqrt(oi.ang);

   if( phi.sin == 0 && psi.sin == 0 )
      u.ang = 0;
   else
      u.ang = atan2(phi.sin,psi.sin);

   u.sin  = sin(u.ang);
   u.cos  = cos(u.ang);
/*
  Computes longitude of the ascending node
*/
   asc.ang = lam-u.ang;
   asc.sin = sin(asc.ang);
   asc.cos = cos(asc.ang);
/*
  Computes the subsatellite geographic latitude
*/
   sublat = atan(REQRPL2*tan(phi.ang));
/*
  Computes the subsatellite longitude
*/
   sublon = asc.ang+atan2(oi.cos*u.sin,u.cos);
/*
  Computes the spacecraft to earth fixed coordinates transformation
  matrix:
      (vector in ecef coordinates) = b * (vector in s/c coordinates)
*/
   b[0][1] = -asc.sin*oi.sin;
   b[1][1] =  asc.cos*oi.sin;
   b[2][1] = -oi.cos;
   b[0][2] = -asc.cos*u.cos+asc.sin*u.sin*oi.cos;
   b[1][2] = -asc.sin*u.cos-asc.cos*u.sin*oi.cos;
   b[2][2] = -phi.sin;
   b[0][0] = -asc.cos*u.sin-asc.sin*u.cos*oi.cos;
   b[1][0] = -asc.sin*u.sin+asc.cos*u.cos*oi.cos;
   b[2][0] =  u.cos*oi.sin;
/*
  Computes the normalized spacecraft position vector in earth fixed
  coordinates - xs.
*/
   r     = (NOMORB+dr)/RADE_EQUA;
   xs.x = -b[0][2]*r;
   xs.y = -b[1][2]*r;
   xs.z = -b[2][2]*r;
/*
  Precomputes q3 (used in lpoint)
*/
   q3 = xs.x*xs.x + xs.y*xs.y + REQRPL2 * xs.z*xs.z - 1.0;
/*
  Computes the attitudes and misalignments if imc is off
*/
   if( imc ){
/*
  Computes the solar orbit angle
*/
      wa = gv_nav.solinc*ts;
/*
  Computes the difference between current time, ts, and the
  exponential time, iparms(62). Note that both times are since epoch.
*/
      te = ts - gv_nav.exptim;
/*
  Computes roll + roll misalignment
*/
      roll = roll + gatt(gv_nav.raawds,wa,te);
/*
  Computes pitch + pitch misalignment
*/
      pitch = pitch + gatt(gv_nav.paawds,wa,te);
/*
  Computes yaw
*/
      yaw = yaw + gatt(gv_nav.yaawds,wa,te);
/*
  Computes roll misalignment
*/
      rma = gatt(gv_nav.rmawds,wa,te);
/*
  Computes pitch misalignment
*/
      pma = gatt(gv_nav.pmawds,wa,te);
/*
  Apply the earth sensor compensation if needed
*/
      roll   = roll + gv_nav.imcrol;
      pitch  = pitch + gv_nav.imcptc;
      yaw    = yaw + gv_nav.imcyaw;
      }
/*
  Computes the instrument to earth fixed coordinates transformation
  matrix - bt
*/
   inst2e(roll,pitch,yaw,b,bt);
   }
return 0;
}
