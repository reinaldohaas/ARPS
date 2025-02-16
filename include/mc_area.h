#include <inttypes.h>

/*********************************************************************

   DESCRIPTION: This is a include file defining data structures used
      in McIDAS area files.


   Copyright 1988, University Corporation for Atmospheric Research
      Not for Resale. All copies to include this notice.

   original code by Glenn Davis UPC
   revisions by Doug Farmer UPC
   more revisions by Keith Brewster, CAPS Sep, 1997

******************************************************************/
#define TYPELEN  4            /* Short strings used as identifiers */
#define COMMENTLEN 32         /* longer strings */

/*#define long int    */   
                              /* On some 64bit platform */

/*
   McIdas AREA DIRECTORY, based on documentation dated 5/87 by R. Dengal
*/
struct area_dir {
   uint32_t      status;      /*  1 */
   uint32_t      type;        /*  2 */
   uint32_t      satid;       /*  3 */
   uint32_t      ndate;       /*  4 - YYDDD */
   uint32_t      ntime;       /*  5 - HHMMSS */
         int32_t lcor;        /*  6 */
         int32_t ecor;        /*  7 */
         int32_t zcor;        /*  8 */
   uint32_t      lsiz;        /*  9 */
   uint32_t      esiz;        /* 10 */
   uint32_t      zsiz;        /* 11 */
   uint32_t      lres;        /* 12 */
   uint32_t      eres;        /* 13 */
   uint32_t      bands;       /* 14 */
   uint32_t      yzprefix;    /* 15 */
   uint32_t      projnum;     /* 16 */
   uint32_t      cdate;       /* 17 */
   uint32_t      ctime;       /* 18 */
   uint32_t      filtmap;     /* 19 */
   uint32_t      imageid;     /* 20 */
   uint32_t      resvid[4];   /* 21-24 */
   char comments[COMMENTLEN]; /* 25-32 */
   uint32_t      calkey;      /* 33 */
   uint32_t      datoffst;    /* 34 */
   uint32_t      navoffst;    /* 35 */
   uint32_t      valcode;     /* 36 */
   uint32_t      pdl[8];      /* 37-44 */
   uint32_t      band8;       /* 45 */
   uint32_t      idate;       /* 46 */
   uint32_t      itime;       /* 47 */
   uint32_t      startscan;   /* 48 */
   uint32_t      doclen;      /* 49 */
   uint32_t      callen;      /* 50 */
   uint32_t      levlen;      /* 51 */
   char stype[TYPELEN];       /* 52 */
   char ctype[TYPELEN];       /* 53 */
   uint32_t      reserved[6]; /* 54-59 */
   uint32_t      supoffst;    /* 60 */
   uint32_t      lensup;      /* 61 */
   uint32_t      resvcal;     /* 62 */
   uint32_t      caloffst;    /* 63 */
   uint32_t      ncomment;    /* 64 */
   };

/*
   McIdas NAVIGATION CODICIL, based on documentation dated 5/87 by D. Santek
   type "GOES" structure
*/
struct nav_goes {
   char type[TYPELEN];        /*   1 */
   uint32_t      iddate;      /*   2 */
   uint32_t      itime;       /*   3 */
   uint32_t      otype;       /*   4 */
   uint32_t      edate;       /*   5 */
   uint32_t      etime;       /*   6 */
   uint32_t      semima;      /*   7 */
   uint32_t      eccen;       /*   8 */
   uint32_t      orbinc;      /*   9 */
   uint32_t      meana;       /*  10 */
   uint32_t      perigee;     /*  11 */
   uint32_t      asnode;      /*  12 */
   uint32_t      declin;      /*  13 */
   uint32_t      rascen;      /*  14 */
   uint32_t      piclin;      /*  15 */
   uint32_t      spinp;       /*  16 */
   uint32_t      deglin;      /*  17 */
   uint32_t      lintot;      /*  18 */
   uint32_t      degele;      /*  19 */
   uint32_t      eletot;      /*  20 */
   uint32_t      pitch;       /*  21 */
   uint32_t      yaw;         /*  22 */
   uint32_t      roll;        /*  23 */
   uint32_t      res1;        /*  24 */
   uint32_t      iajust;      /*  25 */
   uint32_t      iajtim;      /*  26 */
   uint32_t      res2;        /*  27 */
   uint32_t      iseang;      /*  28 */
   uint32_t      resskew;     /*  29 */
   uint32_t      res3;        /*  30 */
   uint32_t      bline1;      /*  31 */
   uint32_t      b1timbeg;    /*  32 */
   uint32_t      b1timend;    /*  33 */
   uint32_t      b1count;     /*  34 */
   uint32_t      bline2;      /*  35 */
   uint32_t      b2timbeg;    /*  36 */
   uint32_t      b2timend;    /*  37 */
   uint32_t      b2count;     /*  38 */
   uint32_t      gamma;       /*  39 */
   uint32_t      gamdot;      /*  40 */
   uint32_t      reserved[80]; /*  41 */
   char memo[COMMENTLEN];    /* 121 */
   };
/*
   Definition of GOES 8, type "GVAR", navigation parameters
*/
struct nav_gvar {
   char sttype[TYPELEN]; /*   1     STTYPE = Satellite type */
   int32_t idntfr;       /*   2     IDNTFR = */
   int32_t imcact;       /*   3     IMCACT = IMC active flag */
   int32_t res1[2];      /*   4-5   Reserved */
   int32_t reflon;       /*   6     REFLON = Reference longitude */
   int32_t refdis;       /*   7     REFDIS = Reference distance from nominal */
   int32_t reflat;       /*   8     REFLAT = Reference latitude */
   int32_t refyaw;       /*   9     REFYAW = Reference yaw */
   int32_t ratrol;       /*  10     RATROL = Reference attitude roll */
   int32_t ratptc;       /*  11     RATPTC = Reference attitude pitch */
   int32_t ratyaw;       /*  12     RATYAW = Reference attitude yaw */
   int32_t etime[2];     /*  13-14  ETIME  = Epoch time */
   int32_t edtime;       /*  15     EDTIME = Delta from epoch time */
   int32_t imcrol;       /*  16     IMCROL = Image motion compensation roll */
   int32_t imcptc;       /*  17     IMCPTC = Image motion compensation pitch */
   int32_t imcyaw;       /*  18     IMCYAW = Image motion compensation yaw */
   int32_t ldr[13];      /*  19-31  LDR    = Long delta from ref parameters */
   int32_t rddr[11];     /*  32-42  RDDR   = Radial dist delta from ref params */
   int32_t dgl[9];       /*  43-51  DGL    = Geocentric lat delta parameters */
   int32_t doy[9];       /*  52-60  DOY    = Orbit yaw delta parameters */
   int32_t solinc;       /*  61     Reserved */
   int32_t exptim;       /*  62     EXPTIM = Exp start time from epoch */
   int32_t raawds[55];   /*  63-117 RAAWDS = Roll attitude angle words */
   int32_t spc1[12];     /* 118-129 Reserved */
   int32_t paawds[55];   /* 130-184 PAAWDS = Pitch attitude angle words */
   int32_t yaawds[55];   /* 185-239 YAAWDS = Yaw attitude angle words */
   int32_t spc2[18];     /* 240-257 Reserved */
   int32_t rmawds[55];   /* 258-312 RMAWDS = Roll misalignment angle words */
   int32_t pmawds[55];   /* 313-367 PMAWDS = Pitch misalignment angle words */
   int32_t imgday;       /* 368     IMGDAY = Image day value (YYDDD) */
   int32_t imgtm;        /* 369     IMGTM  = Image time value (HHMMSS) */
   int32_t imgsnd;       /* 370     IMGSND = Imager/sounder instrument flag */
   int32_t res4[9];      /* 371-379 Reserved */
/*
  These four words were added 5-26-94 to comply w/ the new elug
  numbering started at 380 because these same parameters are used
  in the nav message sent from the ingestor to EVX, and we had to
  start somewhere after the 378 nav parameters
*/
   int32_t iofnc;        /* 380     IOFNC = */
   int32_t iofec;        /* 381     IOFEC = */
   int32_t iofni;        /* 382     IOFNI = */
   int32_t iofei;        /* 383     IOFEI = */
   int32_t res5[257];    /* 384-640 Reserved */
  };
/*
   Navigation structure for polar stereographic plots
*/
struct nav_ps {
   char type[TYPELEN];        /*   1 */
   int32_t npline;            /*   2 */
   int32_t npelem;            /*   3 */
   int32_t stdlat;            /*   4 */
   uint32_t      spacing;     /*   5 */
   int32_t nrmlon;            /*   6 */
   uint32_t      equrad;      /*   7 */
   uint32_t      eccent;      /*   8 */
   uint32_t      reserved[112]; /*  9 */
   char memo[COMMENTLEN];    /* 121 */
   };
/*
   Calibration structure for GVAR data
*/
struct cal_gvar {
   float visbias[8];     /*    1-8 */
   float vis1gain[8];    /*   9-16 */
   float vis2gain[8];    /*  17-24 */
   float albedcon;       /*     25 */
   float ir1bias[4];     /*  26-29 */
   float ir2bias[4];     /*  30-33 */
   float ir1gain[4];     /*  34-37 */
   float ir2gain[4];     /*  38-41 */
   float calresvd[87];   /* 42-128 */
};
/*
   A mcidas area file looks like this:
   Area Directory, followed by
   Navigation Codicil, followed by
   the Image.
  
   N.B.
   Use of this a template is compiler dependent.
   Will work most of the time (32 bit architectures) since the
   arms of the struct are divisible by sizeof(word)
*/

struct mc_area {
   struct area_dir *dir;
   char type[TYPELEN+1];
   struct nav_goes *navg;
   struct nav_gvar *navgv;
   struct nav_ps *navps;
   struct cal_gvar *calgv;
   unsigned char *image;     /* image[imagelen] really */
   void *private;
};
