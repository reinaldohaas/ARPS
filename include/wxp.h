/*********************************************************************
**********************************************************************
*                                                                    *
*                        C O P Y R I G H T                           *
*                                                                    *
*                    Copyright 1988,1989,1990 (C) by                 *
*                    Purdue Research Foundation                      *
*                    Purdue University                               *
*                    West Lafayette, Indiana 47907                   *
*                                                                    *
* This software,  in whole or in part,  may be used  and copied only *
* with the written permission of the  Dept. of Earth and Atmospheric *
* Sciences  via  E.  M. Agee,  Purdue  University,  West  Lafayette, *
* Indiana,  and with the  inclusion  of the above  copyright notice. *
* This  software  or  any  copies  thereof  may  not be  provided or *
* otherwise made  available  to any  other person.  No title  to and *
* ownership of the software is hereby transferred.                   *
*                                                                    *
**********************************************************************

   HEADER FILE: wxp.h
   PROGRAMMER: Dan Vietor
   VERSION: 4.8
   WXP VERSION: 4.8
   DATE: 931120

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

   DESCRIPTION: This header file defines and declares structures
      needed by WXP.

******************************************************************/
#ifndef _h_STDIO
#include <stdio.h>
#endif
#ifndef TYPES_INCLUDED
#include <sys/types.h>
#define TYPES_INCLUDED
#endif

#ifdef vms
#ifndef _h_TIME
#include <time.h>
#endif
typedef char * caddr_t;
#endif

#ifdef AIX_GSL
#define AIX
#define SYSV
#define NETCDF
#define MQUEUE
#define GSLGKS
#define GKS88
#endif

#ifdef AIX_GSS
#define AIX
#define SYSV
#define NETCDF
#define MQUEUE
#define GSSGKS
#define GKS88
#endif

#ifdef AIX_GDT
#define AIX
#define SYSV
#define NETCDF
#define MQUEUE
#define GDTGKS
#define GKS88
#endif

#ifdef AIX_XGKS
#define AIX
#define SYSV
#define NETCDF
#define MQUEUE
#define XGKS
#define GKS86
#define XDUMP
#endif

#ifdef AIX_X11GKS
#define AIX
#define SYSV
#define NETCDF
#define MQUEUE
#define X11GKS
#define GKS88
#define RECALL
#endif

#ifdef ULTRIX_X11GKS
#define ULTRIX
#define SYSV
#define MQUEUE
#define X11GKS
#define NETCDF
#define GKS88
#define NO_USLEEP
#define RECALL
#endif

#ifdef HPUX_X11GKS
#define HPUX
#define SDM
#define SYSV
#define NETCDF
#define MQUEUE
#define X11GKS
#define GKS88
#define NO_USLEEP
#define USE_GETITIMER
#define RECALL
#endif

#ifdef ATT_X11GKS
#define ATT
#define SYSV
#define NETCDF
#define MQUEUE
#define X11GKS
#define GKS88
#define NO_USLEEP
#define NO_STRSTR
#define RECALL
#endif

#ifdef ATT_XGKS
#define ATT
#define SYSV
#define NETCDF
#define MQUEUE
#define XGKS
#define GKS86
#define XDUMP
#define NO_USLEEP
#define NO_STRSTR
#endif

#ifdef PS2AIX
#define GRIB_DUP_WRT
#endif

#ifdef SGI_XGKS
#define SGI
#define SYSV
#define NETCDF
#define MQUEUE
#define XGKS
#define GKS86
#define XDUMP
#define FILE_LOCK
#define GRIB_DUP_WRT
#endif

#ifdef SGI_GL
#define SGI
#define SYSV
#define NETCDF
#define MQUEUE
#define GLGKS
#define GKS88
#define FILE_LOCK
#define GRIB_DUP_WRT
#endif

#ifdef SGI_X11GKS
#define SGI
#define SYSV
#define NETCDF
#define MQUEUE
#define X11GKS
#define GKS88
#define FILE_LOCK
#define GRIB_DUP_WRT
#define NO_USLEEP
#define RECALL
#endif

#ifdef BSD_XGKS
#define BSD
#define NETCDF
#define SOCKET
#define XGKS
#define XDUMP
#define GKS86
#endif

#ifdef SUNOS_XGKS
#define BSD
#define ANSI_C
#define NETCDF
#define MQUEUE
#define XGKS
#define XDUMP
#define GKS86
#define SDM
#endif

#ifdef SOLARIS_X11GKS
#define SOLARIS
#define SYSV
#define X11GKS
#define GKS88
#define NETCDF
#define MQUEUE
#define RECALL
#endif

#ifdef SUNOS_X11GKS
#define SUNOS
#define SYSV
#define X11GKS
#define GKS88
#define NETCDF
#define HAVE_STRICMP
#define MQUEUE
#define RECALL
#endif

#ifdef MSDOS_MSC
#define ANSI_C
#define GKS88
#define MSCGKS
#define HAVE_STRICMP
#define GRIB_DUP_WRT
#endif

#ifdef MSDOS_GSS
#define ANSI_C
#define GSSGKS
#define GKS88
#define HAVE_STRICMP
#define GRIB_DUP_WRT
#endif

#ifdef vms
#define SYSV
#define ANSI_C
#define STDARGS
#define NETCDF
#define VMSGKS
#define GKS86
#define SDM
#define GUISE_RESRC
#endif

#ifdef SDM
#define GUISE_RESRC
#define MCIDAS
#endif

#ifdef NOSDM
#undef SDM
#undef NETCDF
#endif
#ifdef NONETCDF
#undef NETCDF
#endif

#ifdef GUISE_RESRC
#include "udres.h"           /* GUISE definitions */
#endif

#ifdef XGKS
#define WINDOW_GKS
#endif
#ifdef X11GKS
#define WINDOW_GKS
#endif
#ifdef GLGKS
#define WINDOW_GKS
#endif
#ifdef MSCGKS
#define WINDOW_GKS
#endif
/*
   Structure for domain location
*/
typedef struct {
   int proj;                 /* Grid projection */
   float lat[2],lon[2];      /* Latitude and longitude information of grid */
   float dx,dy;              /* Grid spacing */
   int nx,ny;                /* Number of grid points in each direction */
   float param[5];           /* Extra projection parameters */
   } DOM_params;
/*
   Surface data structure
*/
#define COMMENT_LENGTH 35
typedef struct {
   char id[11];              /* Station identifier */
   char area[11];            /* Regional identifier */
   int type;                 /* Report type (e.g. SA, RS, SP, COR) */
   long wmo;                 /* WMO station number */
   float lat;                /* Station latitude */
   float lon;                /* Station longitude */
   float elev;               /* Station elevation */
   int prior;                /* Station plotting priority */
   int time;                 /* Observation time */
   float cld_height[5];      /* Cloud height in feet */
   char cld_cover[5];        /* Cloud cover C,S,s,B,b,O,o,X,x */
   int num_cld;              /* Number of cloud layers */
   float vis;                /* Visibility in miles */
   char weather[11];         /* Type of weather */
   float pmb;                /* Sea level pressure */
   float t,td;               /* Temperature and dewpoint */
   float w_dir,w_spd,w_gust; /* Wind direction, speed, and gust */
   float alt;                /* Altimeter setting */
   int ptend_t;              /* Type of pressure tendency */
   float ptend;              /* Pressure tendency */
   float rain6;              /* 3 or 6 hour precipitation */
   float rain12;             /* 12 hour precipitation */
   float rain24;             /* 24 hour precipitation */
   char cld_type[3];         /* Low, medium and high level cloud types */
   float max_t;              /* Maximum temperature */
   float min_t;              /* Minimum temperature */
   float max6_t;             /* 6 hour Maximum temperature */
   float min6_t;             /* 6 hour Minimum temperature */
   float scover;             /* Snow cover information in inches */
   float solar;              /* Total minutes of sunshine */
   float sst;                /* Sea surface temperature */
   float wav_per;            /* Sea wave period */
   float wav_hgt;            /* Sea wave height */
   char comment[COMMENT_LENGTH];/* The comment field */
   } SFC_report;
/*
   Set up the number of levels in the UA_report structure
*/
#ifdef MSDOS
#define NUMMAN 20            /* Number of mandatory levels */
#define NUMSIG 28            /* Number of significant levels */
#define NUMWND 28            /* Number of wind levels */
#define NUMMRG 70            /* Number of merged levels */
#else
#define NUMMAN 20
#define NUMSIG 50
#define NUMWND 50
#define NUMMRG 120
#endif
/*
   Upper air data structure
*/
typedef struct {
   char id[11];              /* Station identifier */
   char area[11];            /* Regional identifier */
   long wmo;                 /* WMO station number */
   float lat;                /* Station latitude */
   float lon;                /* Station longitude */
   float elev;               /* Station elevation */
   int time;                 /* Observation time */
   int type;                 /* Observation type */
   struct {                  /* Mandatory level information */
      float p;               /*   Pressure */
      float h;               /*   Height */
      float t;               /*   Temperature */
      float td;              /*   Dewpoint */
      float dir;             /*   Wind direction */
      float spd;             /*   Wind speed */
      } man[NUMMAN];
   int numman;               /* Number of mandatory levels */
   struct {                  /* Significant level information */
      float p;               /*   Pressure */
      float t;               /*   Temperature */
      float td;              /*   Dewpoint */
      } sig[NUMSIG];
   int numsig;               /* Number of significant levels */
   struct {                  /* Wind level information */
      float h;               /*   Height */
      float dir;             /*   Wind direction */
      float spd;             /*   Wind speed */
      } wnd[NUMWND]; 
   int numwnd;               /* Number of wind levels */
   } UA_report;
/*
   Upper air merged data structure
*/
typedef struct {
   float p[NUMMRG];          /* Pressure */
   float h[NUMMRG];          /* Height */
   float t[NUMMRG];          /* Temperature */
   float td[NUMMRG];         /* Dewpoint */
   float dir[NUMMRG];        /* Wind direction */
   float spd[NUMMRG];        /* Wind speed */
   int num;                  /* Number of merged levels */
   } UA_merge;
/*
   Upper air level data structure
*/
typedef struct {
   float p;                  /* Pressure */
   float h;                  /* Height */
   float t;                  /* Temperature */
   float td;                 /* Dewpoint */
   float dir;                /* Wind direction */
   float spd;                /* Wind speed */
   } UA_level;
/*
   Structure for surface information
*/
typedef struct {
   char datim[20];           /* Date and time information */
   char file[20];            /* Data file name */
   char *cmdline;            /* Command line history information */
   int year;                 /* Year of grid file */
   int month;                /* Month of grid file */
   int day;                  /* Day of grid file */
   int hour;                 /* Hour of grid file */
   int num_stat;             /* Number of stations in file */
   } SFC_info;
/*
   Structure for upper air information
*/
typedef struct {
   char datim[20];           /* Date and time information */
   char file[20];            /* Data file name */
   char *cmdline;            /* Command line history information */
   int year;                 /* Year of grid file */
   int month;                /* Month of grid file */
   int day;                  /* Day of grid file */
   int hour;                 /* Hour of grid file */
   } UA_info;
/*
   Structure for grid information
*/
typedef struct {
   char datim[40];           /* Reference date and time information */
   int max_x,max_y;          /* Maximum x and y dimensions */
   int grid_type;            /* Grid type for file */
   int year;                 /* Year of grid file */
   int month;                /* Month of grid file */
   int day;                  /* Day of grid file */
   int hour;                 /* Hour of grid file */
   } GRD_info;
/*
   Structure for grid data
*/
typedef struct {
   char datim[50];           /* Date and time information */
   char info[50];            /* Label for grid information */
   int source;               /* Source of grid */
   int version;              /* Version of grid */
   int model;                /* Model used to produce grid */
   int type;                 /* Grid type used in grid */
   int ltype;                /* Type of vertical level */
   long level;               /* Level of grid */
   int time;                 /* Time of grid */
   int var;                  /* Grid variable */
   DOM_params domain;        /* The grid location */
   float *grid;              /* A pointer to the grid */
   } GRD_report;
#define PROJ_LATLON  2
#define PROJ_LL      2
#define PROJ_PSTEREO 1
#define PROJ_PS      1
#define PROJ_MERC    3
#define PROJ_LAMB    4
#define PROJ_ORTHO   5
#define PROJ_SAT     9
#define PROJ_is_map(a) ((a) >= PROJ_PS && (a) <= PROJ_SAT )
#define PROJ_PIXEL  99
#define PROJ_XY     20
#define PROJ_XZ     21
#define PROJ_XLOGP  22
#define PROJ_POLAR  23
#define PROJ_CAT    24
#define PROJ_SKEWT  40
#define PROJ_STUVE  41
/*
   Structure for radar information
*/
typedef struct {
   char datim[40];           /* Reference date and time information */
   char *cmdline;            /* Command line history information */
   int year;                 /* Year of grid file */
   int month;                /* Month of grid file */
   int day;                  /* Day of grid file */
   int hour;                 /* Hour of grid file */
   short summary[120][90];   /* Radar summary information */
   } RAD_info;
/*
   Structure for manual digitized radar reports
*/
typedef struct {
   char id[11];          /* Station identifier */
   char region[11];      /* Region to which the station is located */
   float x,y;            /* Station location */
   int type;             /* Station type */
   int time;             /* Report time */
   char config[10];      /* Configuration of echoes */
   struct {
      int area_cov;      /* Area coverage */
      char type[5];      /* Type of precipitation */
      int intens;        /* Intensity of precipitation */
      } prec[3];
   int numprec;          /* Number of precip areas */
   char trend[4];        /* Trend of precip coverage */
   struct {              /* Location of precipitation */
      int dir,rad;       /* Direction and radius of precip area */
      } loc[5];
   struct {              /* Movement of precipitation */
      char type;         /* Type of movement: Area,Cell,Line */
      int dir,spd;       /* Direction and speed */
      } move[3];
   struct {              /* Cloud tops */
      int h;             /* Height */
      int dir,rad;       /* Location */
      } tops;
   char rem[20];         /* Remarks */
   } RAD_report;
/*
   Raw data structure
*/
typedef struct {
   char id[11];              /* Station identifier */
   int type;                 /* Report type (e.g. SA, RS, SP, COR) */
   long wmo;                 /* WMO station number */
   float lat;                /* Station latitude */
   float lon;                /* Station longitude */
   float elev;               /* Station elevation */
   int prior;                /* Station plotting priority */
   float value;              /* Value at station */
   char string[30];          /* String value for a station */
   } RAW_report;
/*
   Structure for raw information
*/
typedef struct {
   char datim[50];           /* Date and time information */
   char info[50];            /* Parameter information */
   char file[20];            /* Data file name */
   char *cmdline;            /* Command line history information */
   int year;                 /* Year of grid file */
   int month;                /* Month of grid file */
   int day;                  /* Day of grid file */
   int hour;                 /* Hour of grid file */
   int num_stat;             /* Number of stations in file */
   } RAW_info;
/*
   Structure for profiler data
*/
typedef struct {
   char id[10];
   float lat,lon,elev;
   struct {
      float h;
      float u,v,w;
      } lev[50];
   } PRF_report;
/*
   Structure for image data
*/
typedef struct {
   int source;               /* Source of image */
   int version;              /* Version of image */
   char datim[50];           /* Date and time information */
   char info[50];            /* Label for image information */
   int type;                 /* Image type used in image */
   DOM_params domain;        /* The image domain */
   unsigned char *image;     /* A pointer to the image */
   } IMG_report;
/*
   Structure for city database information
*/
typedef struct {
   int dtype;                /* Type of city database file */
   } CTY_info;
/*
   Structure for map database information
*/
typedef struct {
   int num_blocks;           /* Number of map blocks */
   int dtype;                /* Type of map database file */
   } MAP_info;
/*
   Structure for MOS data
*/
typedef struct {
   char id[4];               /* The station identifier */
   int t[21];                /* The temperature */
   int td[21];               /* The dewpoint */
   int w_dir[21];            /* The wind direction */
   int w_spd[21];            /* The wind speed */
   int pop6[21];             /* The 6 hour probability of precipitation */
   int pop12[21];            /* The 12 hour probability of precipitation */
   int pot6[21];             /* The 6 hour probability of thunderstorms */
   int pot12[21];            /* The 12 hour probability of thunderstorms */
   int post6[21];            /* The 6 hour prob of severe thunderstorms */
   int post12[21];           /* The 12 hour prob of severe thunderstorms */
   int ptype[21];            /* The precipitation type */
   int ext_t[21];            /* The extreme temperature */
   int cover[21];            /* The cloud cover */
   int hght[21];             /* The cloud ceiling */
   int vis[21];              /* The visibility */
   int obvis[21];            /* The obscuration to visibility */
   int qpf6[21];             /* The 6 hour quantitative precipitation */
   int qpf12[21];            /* The 12 hour quantitative precipitation */
   int qsf6[21];             /* The 6 hour quantitative snow */
   int qsf12[21];            /* The 12 hour quantitative snow */
   } FO_report;

typedef struct {
   char id[4];               /* The station identifier */
   int t[9][4];              /* The temperature */
   int rh[9][4];             /* The relative humdity at 4 levels */
   int vert[9];              /* The 700 mb vertical velocity */
   int li[9];                /* The lifted index */
   int hght[9];              /* The 1000-500 mb thickness */
   int dir[9];               /* The boudary layer wind direction */
   int spd[9];               /* The boudary layer wind speed */
   int psl[9];               /* The sea level pressure */
   int prec[9];              /* The quantitative precipitation */
   } FO_ua_report;

typedef struct {
   char id[4];               /* Station identifier */
   int max[9];               /* Maximum temperature */
   int nmax[9];              /* Normal maximum temperature */
   int min[9];               /* Minimum temperature */
   int nmin[9];              /* Normal minimum temperature */
   int cover[9][2];          /* The cloud cover */
   int w_spd[9][2];          /* The wind speed */
   int pop12[9][2];          /* The 12 hour probability of precipitation */
   int pos12[9][2];          /* The 12 hour probability of snow */
   int pop24[9];             /* The 24 hour probability of precipitation */
   } FO_mrf_report;
/*
   Structure for grid parameters
*/
typedef struct {
   float max,min;            /* The maximum and minimum values */
   float intrvl;             /* The contour interval */
   } GRD_params;
/*
   Objective analysis parameters structure
*/
typedef struct {
   float filter;             /* Filter parameter for smoothing */
   float converg;            /* Convergence parameter for secondary smoothing */
   float rad_inf;            /* Radius of influence */
   int num_pass;             /* Number of passes through data array */
   int num_field;            /* Number of data fields to grid */
   int min_sta;              /* Minimum number of stations affecting each
                                gridpoint */
   } OA_params;
/*
   Structure to save city database info
*/
typedef struct {
   char name[17];            /* Name of city */ 
   char state[3];            /* State city is located */
   char region[3];           /* Region/country city is located */
   char id[6];               /* FAA identifier */
   int prior;                /* Station priority */
   char mos;                 /* Specifies station is a MOS station */
   char ua;                  /* Specifies station is upper air station */
   float lat,lon;            /* Latitude and longitude of the station */
   float el;                 /* Station elevation above sea level */
   long wmo;                 /* 5 digit WMO identifier */
   int num;                  /* Number of occurrences */
   } CTY_data;
/*
   Structure to save city database info
*/
typedef struct {
   char id[6];      /* The station identifier in the city database */
   long wmo;        /* The WMO station number */
   float x;         /* The x p.s. coordinate from the city database */
   float y;         /* The y p.s. coordinate from the city database */
   float el;        /* The elevation from the city database */
   } CTY_locate;
/*
   Structure for program arguments
*/
typedef struct {
   char opt;                 /* Option number */
   int next;                 /* Location of next option */
   int argc;                 /* number of sub-arguments */
   char *argv[20];           /* List of sub-arguments */
   } argument;
/*
   Structure for scanp data
*/
typedef struct {
   int start;            /* The starting point of the match */
   int end;              /* The ending point of the match */
   } parse_param;
/*
   Menu items
*/
typedef struct {
   char *name;
   char *abbr;
   int index;
   int disp[10];
   } menu_item;
typedef menu_item menu_items[];

typedef enum 
   { CV_none,CV_bool,CV_int,CV_short,CV_long,CV_float,CV_string,
     CV_strcpy,CV_path,CV_mess,CV_current,CV_region,CV_name,
     CV_ingestor,CV_device,CV_color,CV_domain,CV_strings } DEF_convert;

#ifdef MSDOS
typedef char * caddr_t;
#endif
#define POS1 (caddr_t)1
#define POS2 (caddr_t)2
#define POS3 (caddr_t)3
#define POS4 (caddr_t)4
#define POS5 (caddr_t)5
#define POS6 (caddr_t)6
#define POS7 (caddr_t)7
#define POS8 (caddr_t)8
#define POS9 (caddr_t)9
#define POS (caddr_t)10
#ifndef VOIDP
#define VOIDP char *
#endif
/*
   WXP resource parameters
*/
#include "wxpresrc.h"
/*
   Structure for defaults
*/
typedef struct {
   char *def;
   char *env;
   char *cmd;
   DEF_convert cvt;
   char *out;
   int  *num;
   } WXP_defobj;
typedef WXP_defobj WXP_deflist[];
/*
   Setup function return types
*/
/* SYSTEM library */
char *WXP_version();
char *WXP_getres();
char *WXP_promptl();
char *WXP_prompts();
float WXP_promptf();
FILE *WXP_fopen();
char *WXP_tab_search();
char *WXP_tab_next();
char *WXP_get_menu_abbr();
char *WXP_proc_name();
long WXP_fsize();
/* FILE library */
char *get_file();
char *get_curfile();
char *get_oldfile();
char *set_filename();
char *set_filedate();
char *set_num();
/* CALC library */
float vapor_pr();
float wind_chill();
float heat_index();
float sfc_pr_alt();
float mix_ratio();
float spec_hum();
float pot_temp();
float virt_temp();
float virt_pot_temp();
float temp_on_satad();
float temp_on_dryad();
float thetae_on_satad();
float temp_vapor_pr();
float temp_on_mixrat();
float rel_hum();
float mean_rel_hum();
float conv_temp();
int lcl_info();
float eq_pot_temp();
float interp_temp();
float pr_ccl();
void layer_mean();
float precip_water();
float freezing_level();
float lifted_index();
float show_index();
float vertical_totals();
float cross_totals();
float total_totals();
float k_index();
float sweat_index();
float energy_index();
/* TOOLS library */
long num_sec();
long num_sec_date();
long gettime();
float interp();
float linear_interp();
char *create_string();
char *create_strcat();
char *get_progname();
char *get_color_name();
struct tm *getgmtime();
float dir_decode();
float map_dir_adj();
float map_uvdir_adj();
float map_dir();
float MOS_cld_hght();
float MOS_vis();
float MOS_quant_prec6();
float MOS_quant_prec12();
long GRD_get_pos();
char *RAW_get_field_info();
#ifdef MSDOS
char *getenv();
#endif
#ifdef NO_STRSTR
char *strstr();
#endif
#ifdef vms
double atof();
#endif
/* Define the missing data value */
#define MISS -9999
/*
   Define some simple unit conversion macros
*/
#define CtoF(t) ((t)!=MISS?((t)*9/5+32):MISS)
#define KtoF(t) ((t)!=MISS?(((t)-273.15)*9/5+32):MISS)
#define CtoK(t) ((t)!=MISS?((t)+273.15):MISS)
#define FtoK(t) ((t)!=MISS?(((t)-32.)*5/9+273.15):MISS)
#define FtoC(t) ((t)!=MISS?(((t)-32.)*5/9):MISS)
#define KtoC(t) ((t)!=MISS?((t)-273.15):MISS)
#define INtoMB(p) ((p)!=MISS?((p)*33.85399):MISS)
#define MBtoIN(p) ((p)!=MISS?((p)*.0295386):MISS)
#define KTtoMPS(w) ((w)!=MISS?((w)*.5144):MISS)
#define MPStoKT(w) ((w)!=MISS?((w)*1.944):MISS)
#define KTtoKMH(l) ((l)!=MISS?((l)*1.852):MISS)
#define KMHtoKT(l) ((l)!=MISS?((l)*0.53996):MISS)
#define FTtoM(l) ((l)!=MISS?((l)*.3048):MISS)
#define MtoFT(l) ((l)!=MISS?((l)*3.2808):MISS)
#define MtoMI(l) ((l)!=MISS?((l)*6.2136E-4):MISS)
#define INtoMM(l) ((l)!=MISS?((l)*25.4):MISS)
#define MMtoIN(l) ((l)!=MISS?((l)*.03937):MISS)
#define INtoCM(l) ((l)!=MISS?((l)*2.54):MISS)
#define CMtoIN(l) ((l)!=MISS?((l)*.3937):MISS)
#define SMtoKM(l) ((l)!=MISS?((l)*1.60934):MISS)
#define KMtoSM(l) ((l)!=MISS?((l)*.62137):MISS)
#define SMtoNM(l) ((l)!=MISS?((l)*.86897):MISS)
#define NMtoSM(l) ((l)!=MISS?((l)*1.1508):MISS)
#define NMtoKM(l) ((l)!=MISS?((l)*1.852):MISS)
#define KMtoNM(l) ((l)!=MISS?((l)*.53996):MISS)
#ifndef MIN
#define MIN(x,y) ((x)<(y)? (x) : (y))
#endif
#ifndef MAX
#define MAX(x,y) ((x)>(y)? (x) : (y))
#endif

/* Define the typical constants */
#define PI      3.1415926535898
#define TWO_PI  (2.*PI)
#define DRC     (PI/180.)
#define RDC     (180./PI)
#define R 287.05              /* Gas constant */
#define G 9.8                 /* Acceleration due to gravity */
#define CP 1004               /* Specific heat at constant pressure */
#define CV 717                /* Specific heat at constant volume */
#define NO 0
#define YES 1
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif
#define WXP_EOF(x) ((x) == -1)
#define OFFSET_LL 360
#define OFFSET_MERC 360*63.71*DRC
/*
   Define data file types
*/
#define FT_unknown 0
#define FT_WXP_1 1
#define FT_WXP_2 2
#define FT_WXP_3 3
#define FT_WXP_4 4
#define FT_netCDF_1 11
#define FT_netCDF_2 12
#define FT_GRIB_1 21
#define FT_GRIB_2 22
#define FT_none 99
#define is_FT_WXP(a) ( a >= FT_WXP_1 && a <= FT_WXP_4 )
#define is_FT_netCDF(a) ( a >= FT_netCDF_1 && a <= FT_netCDF_2 )
#define is_FT_GRIB(a) ( a >= FT_GRIB_1 && a <= FT_GRIB_2 )
/*
  Extension types
*/
#define EXT_WXP 0
#define EXT_CDF 1
#define EXT_WMO 2
#define EXT_NPS 3
#define EXT_MCA 4
/*
   Upper air data types
*/
#define UA_SFC 0
#define UA_1000 1
#define UA_925 2
#define UA_850 3
#define UA_700 4
#define UA_500 5
#define UA_400 6
#define UA_300 7
#define UA_250 8
#define UA_200 9
#define UA_150 10
#define UA_100 11
#define UA_70 12
#define UA_50 13
#define UA_30 14
#define UA_20 15
#define UA_10 16
#define UA_TROP 17
#define UA_WIND 18
#define UA_ALL 30
#define UA_ID 31
#define UA_STN 32
#define UA_SND 33
#define UA_VAR 34
/*
   Surface data types
*/
#define SFC_ALL 0
#define SFC_ID 1
#define SFC_STN 2
/* Define surface report types */
#define SFC_SA 1
#define SFC_SP 2
#define SFC_COR 4
#define SFC_AMOS 8
/*
   Radar data types
*/
#define RAD_ALL 0
#define RAD_ID 1
#define RAD_STN 2
/*
   Profiler data types
*/
#define PRF_ALL 0
#define PRF_ID 1
#define PRF_STN 2
/*
   Raw data types
*/
#define RAW_NONE (char *)0
#define RAW_ALL (char *)1
#define RAW_ID (char *)2
#define RAW_STN (char *)3
#define RAW_DATA (char *)4
#define RAW_LAST (char *)5
#define RAW_VAL(x) (char *)10+x
/*
   Define report writing types
*/
#define REP_NEW 1
#define REP_COR 2
/*
   Grid data types
*/
#define GRD_ALL 0
#define GRD_INQ 1
/*
   Grid spacing types
*/
#define GRD_REG 0
#define GRD_LATLON 1
/*
   Define MOS Forecast times
*/
#define FO_ALL 22
#define FO_ID 21
#define FO_OB 0
#define FO_03 1
#define FO_06 2
#define FO_09 3
#define FO_12 4
#define FO_15 5
#define FO_18 6
#define FO_21 7
#define FO_24 8
#define FO_27 9
#define FO_30 10
#define FO_33 11
#define FO_36 12
#define FO_39 13
#define FO_42 14
#define FO_45 15
#define FO_48 16
#define FO_51 17
#define FO_54 18
#define FO_57 19
#define FO_60 20
/*
   Define MRF MOS Forecast times
*/
#define FOM_OB 0
#define FOM_D1 1
#define FOM_D2 2
#define FOM_D3 3
#define FOM_D4 4
#define FOM_D5 5
#define FOM_D6 6
#define FOM_D7 7
#define FOM_NM 8

#ifdef MSDOS
#define path_delim '\\'
#define path_str "\\"
#define geteuid() 0
#define alarm(a) 0
#else
#define path_delim '/'
#define path_str "/"
#endif

#ifdef BSD
#ifdef IBMRT
double drem();
#define fmod(x,y) (drem(x,y) < 0 ? drem(x,y)+y : drem(x,y))
#endif
#undef tolower
#define tolower(x) (x >= 'A' && x <= 'Z' ? x + 32 : x)
#endif

#define PATH_LENGTH 50
#define FILE_LENGTH 70
/*
   define subprogram flags
*/
#define R_CHDIR      1
#define R_BACKG      2
/*
   define file access flags
*/
#define F_READ       0
#define F_WRITE      1
#define F_APPEND     2
/*
   define wxp flags
*/
#define F_CUR        0
#define F_ING        1
#define F_DISP       2
#define F_NAME       3
#define F_MODEL      4
#define F_REG        5
#define F_LEVEL      6
#define F_TIME       7
#define F_VAR        8
#define F_PRIOR     10
#define F_FILL      11
#define F_HOUR      12
#define F_MIN       13
#define F_FILE      15
#define F_CONV      16
#define F_TYPE      17
#define F_VARNUM    18
/*
   define data sources
*/
#define FAA604       1
#define DDPLUS       2
#define INTERN       3
#define NPS          4
#define UNIDATA      5
#define WXP          6
#define HDS          7
#define DDS          8
#define PPS          9
/*
   define possible ingestors
*/
#define ING_NONE    99
#define ING_RAW     10
#define ING_WXP      1
#define ING_SPLIT    2
#define ING_LDM      5
/*
   define possible file naming conventions
*/
#define FNC_HMDY     1
#define FNC_YMDH     2
#define FNC_HDMY     3
#define FNC_YMDH_T  12
/*
   define file type
*/
#define RAWUNK      -1
#define is_RAWFILE(a) ( a >= 0 && a < 40 )
#define is_WMOFILE(a) ( a >= RAWDAT && a < RAWSUM )
#define RAWDAT       0
#define RAWSFC       1
#define RAWUPA       2
#define RAWUA        2
#define RAWRAD       3
#define RAWSAT       4
#define RAWMOD       5
#define RAWSYN       6
#define RAWFOR       7
#define RAWSEV       8
#define RAWCLI       9
#define RAWSUM      10
#define is_AREAFILE(a) ( a >= RAWAREA && a <= RAWMCRE )
#define RAWAREA     11
#define RAWMCVI     12
#define RAWMCIR     13
#define RAWMCWV     14
#define RAWMCFL     15
#define RAWMCRE     16
#define is_CVTFILE(a) ( a >= 40 && a < 50 )
#define CVTRAW      40
#define CVTSFC      41
#define CVTUPA      42
#define CVTRAD      43
#define CVTSYN      44
#define is_PARFILE(a) ( a >= 50 && a < 60 )
#define PARFRT      51
#define PARWWS      52
#define is_GRBFILE(a) ( a >= 60 && a < 70 )
#define GRID        61
#define GRIB        62
#define RAWGRIB     63
#define is_GRDFILE(a) ( a >= 70 && a < 80 )
#define GRDSFC      71
#define GRDUPA      72
#define GRDMOD      73
#define GRDGRB      74
#define RAWGOES     80
#define TMPUPA      90

#define LATEST    9999
#define LATEST_from(a) (9999+a)
/*
   File parameter
*/
#define FILE_OVER    1
#define FILE_APPEND  2
#define FILE_USE     3
/*
   Model types
*/
#define MOD_RAW     99
#define MOD_LFM      1
#define MOD_NGM      2
#define MOD_ETA      3
#define MOD_RUC      4
#define MOD_INIT     5
#define MOD_MRF     10
#define MOD_MRF1    11
#define MOD_MRF2    12
#define MOD_MRF3    13
#define MOD_MRF4    14
#define MOD_ERF     15
#define MOD_ERF1    16
#define MOD_ERF2    17
#define MOD_ECMWF   20
#define MOD_SWW     31
#define MOD_SST     32
#define MOD_SCA     33
#define MOD_NGM105  42
#define MOD_ETA105  43
/*
   MOS Model types
*/
#define MOS_NONE     0
#define MOS_LFM      1
#define MOS_NGM      2
#define MOS_LFM_UA   3
#define MOS_NGM_UA   4
#define MOS_MRF      5
#define MOS_MRFX     6
#define MOS_ETA      7
#define MOS_ETA_UA   8
#define MOS_ALL     10
/*
   Variable types
*/
#define VAR_PR       1
#define VAR_HH       2
#define VAR_ST       3
#define VAR_RH       4
#define VAR_UW       5
#define VAR_VW       6
#define VAR_VV       7
#define VAR_LI       8
#define VAR_QP       9
#define VAR_SC      10
#define VAR_WH      11
#define VAR_WP      12
#define VAR_WD      13
#define VAR_SS      14
/*
   Pressure level types
*/
#define LEV_1000     1
#define LEV_950      2
#define LEV_900      3
#define LEV_850      4
#define LEV_800      5
#define LEV_750      6
#define LEV_700      7
#define LEV_650      8
#define LEV_600      9
#define LEV_550     10
#define LEV_500     11
#define LEV_450     12
#define LEV_400     13
#define LEV_350     14
#define LEV_300     15
#define LEV_250     16
#define LEV_200     17
#define LEV_150     18
#define LEV_100     19
#define LEV_TROP    20
#define LEV_WIND    21
#define LEV_SND     22
#define LEV_VAR     23
/*
   Forecast time types
*/
#define FTIME_00     1
#define FTIME_06     2
#define FTIME_12     3
#define FTIME_18     4
#define FTIME_24     5
#define FTIME_30     6
#define FTIME_36     7
#define FTIME_42     8
#define FTIME_48     9
#define FTIME_60    10
#define FTIME_D3    11
#define FTIME_E3    12
#define FTIME_D4    13
#define FTIME_E4    14
#define FTIME_D5    15
#define FTIME_E5    16
#define FTIME_D6    17
#define FTIME_E6    18
#define FTIME_D7    19
#define FTIME_E7    20
#define FTIME_D8    21
#define FTIME_E8    22
#define FTIME_D9    23
#define FTIME_E9    24
#define FTIME_D0    25
/*
   grid function types
*/
#define M_ADD        1
#define M_SUB        2
#define M_MULT       3
#define M_DIV        4
#define M_MOD        5
/*
   message types
*/
#define M_NONE       -1
#define M_PRINT      0
#define M_ERROR      1
#define M_WARN       2
#define M_MESS       3
#define M_OUT1       4
#define M_OUT2       5   /* DEFAULT */
#define M_OUT3       6
#define M_OUT4       7
#define M_DEBUG      8
/*
   menu types
*/
#define MENU_EXIT    1
#define MENU_MAIN    2
#define MENU_HELP    4
/*
   exit() error return values
*/
#define WNOERR       0   /* No error */
#define WEMISC       1   /* Miscellaneous error */
#define WNORUN       2   /* Unable to run program */
#define WNOMEM       3   /* Not enough memory */
#define WNOGKS       4   /* Unable to use GKS */
#define WNOLINE      5   /* Dataline not responding */
#define WINGRUN      6   /* Ingest already running */
#define WNOING       7   /* Ingest not running */
#define WNOBATCH     8   /* Cannot run in batch mode */
#define WNOINIT      9   /* Cannot initialize program parameters */
#define WNOINPF     10   /* Input data file not found */
#define WNODATF     11   /* Database file not found */
#define WNOOUTF     12   /* Unable to create output file */
#define WNOTMPF     13   /* Unable to create temporary file */
#define WNOIPC      14   /* Unable to create IPC */
#define WNOINP      15   /* Missing input data */
#define WLPATH      20   /* Pathname too long */
#define WLFILE      21   /* Filename too long */
#define WNODATA     30   /* No data found */
#define WNOSTA      31   /* Station identifier not found */
#define WIFILE      40   /* Invalid filename */
#define WISTA       41   /* Invalid station identifier */
#define WIBUL       42   /* Invalid bulletin identifier */
#define WIREG       43   /* Invalid region */
#define WIVAR       44   /* Invalid variable */
#define WILEV       45   /* Invalid level */
#define WIVAL       46   /* Invalid value */
#define WIPSIZE     47   /* Invalid plot size */
#define WIPRIOR     48   /* Invalid priority */
#define WISCALE     49   /* Invalid scale factor */
#define WIOBJAN     50   /* Invalid objective analysis parameter */
#define WIINT       51   /* Invalid contour interval */
#define WIGRID      52   /* Incompatible grid files */
#define WIDSRC      53   /* Invalid data source */
#define WITIM       45   /* Invalid forecast time */
#define WITEMP      60   /* Invalid temperature */
#define WIDPT       61   /* Invalid dewpoint */
#define WIPRES      62   /* Invalid pressure */
#define WIHGHT      63   /* Invalid height */
#define WISPD       64   /* Invalid wind speed */


