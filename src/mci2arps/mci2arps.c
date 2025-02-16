/*  
*
*   remapsat.c
*
*   remapsat
*
*   Program to read McIDAS AREA files and remap to
*   ARPS grid
*
*   The program also writes an hdf image file for diagnostics.
*
*   McIDAS navigation borrows heavily from WXP xsat software.
*
*   Keith Brewster
*   CAPS/Univ of Oklahoma
*   September, 1997
*
*   Updated for hdf file writing for greater portability.
*   Leilei Wang and Keith Brewster
*   April, 2002
*
*   Updated to run on 64-bit systems.  Fortran code uses type INTEGER which
*   is 4-bytes for both 32-bit and 64-bit compiles.  This code used to use
*   type LONG for the same variables.  Type long is 4-bytes for 32-bit
*   compiles, however, it is 8-bytes for 64-bit compilers.  This leads to
*   core dumps.  File "mc_area.h" has the same problem.
*   Kevin W.Thomas
*   August, 2004
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <ctype.h>
#include <stdlib.h>

#include "hdf.h"
#include "wxp.h"
#include "mc_area.h"

static char datim[30];           /* Date/time string */
static char *mnth[] = {          /* Month labels */
  "JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC" };

static char *out_string = NULL; /* Output data type */

static int jday[] = {0,31,59,90,120,151,181,212,243,273,304,334 };
static int jdayl[] ={0,31,60,91,121,152,182,213,244,274,305,335 };

 
/* Satellite ID numbers from Appendix B of McIDAS Programmers' Manual */
/* http://www.ssec.wisc.edu/mcidas/doc/users_guide/current/app_c-1.html */

#define G06VIS 30   /* GOES-06 Visible */
#define G06IR  31   /* GOES-06 IR      */
#define G07VIS 32   /* GOES-07 Visible */
#define G07IR  33   /* GOES-07 IR      */
#define G08IM  70   /* GOES-08 Imager  */
#define G08SND 71   /* GOES-08 sounder */
#define G09IM  72   /* GOES-09 Imager  */
#define G09SND 73   /* GOES-09 Sounder */
#define G10IM  74   /* GOES-10 Imager  */
#define G10SND 75   /* GOES-10 Sounder */
#define G11IM  76   /* GOES-11 Imager  */
#define G11SND 77   /* GOES-11 Sounder */
#define G12IM  78   /* GOES-12 Imager  */
#define G12SND 79   /* GOES-12 Sounder */
#define G13IM  180  /* GOES-13 Imager  */
#define G13SND 181  /* GOES-13 Sounder */
#define G14IM  182  /* GOES-14 Imager */
#define G14SND 183  /* GOES-14 Sounder */
#define G15IM  184  /* GOES-15 Imager */
#define G15SND 185  /* GOES-15 Sounder */

#define GEORNG 42.164E06  /* Geostationary radius from earth center */

DOM_params domain;    /* Default domain parameters */
DOM_params sat_domain; /* Domain of the raw satellite image */

#ifdef UNDERSCORE
  #define inisatarps inisatarps_
  #define ctim2abss ctim2abss_
  #define solr1r2 solr1r2_
  #define solrsc1 solrsc1_
  #define solrsc2 solrsc2_
  #define solcorset solcorset_
  #define coldfilt coldfilt_
  #define wtsatfld wtsatfld_
  #define lltoxy lltoxy_
#endif
void inisatarps();
void ctim2abss();
void solr1r2();
void solrsc1();
void solrsc2();
void solcorset();
void coldfilt();
void wtsatfld();
void lltoxy();

int AREA_read();
void sat_init();
void sat_tran();
void sat_line_elem();
void solsrc2();
void mapcoord();

int yflip(nx,ny,infld,outfld)

int nx;
int ny;
float *infld;
float *outfld;

{
   int i,j,jin,jout;

   for ( j = 0; j < ny ; j++ ) {
      jin = nx*j;
      jout =  nx*((ny-j)-1);
      for ( i = 0; i < nx ; i++ ) {
          *(outfld+jout+i)=*(infld+jin+i);
          }
      }
   return(0);
}

int rescale(image,len)

uint8 *image;
long len;

{
   int k;
   uint8 brit,immax,immin;
   float scale;

   immax = 0;
   immin = 255;
   for ( k=1 ; k < len ; k++ ) {
      if ( (brit = *(image+k)) != 0 ) {
         if ( brit < immin )
            immin = brit;
         if ( brit > immax )
            immax = brit;
         }
      }
   
   printf( " Before rescale, immin = %d  immax = %d\n",immin,immax);

   scale = 254. / ((float)(immax - immin));

   for ( k=1 ; k < len ; k++ ) {
      if ( (brit = *(image+k)) != 0 ) {
         brit = *(image+k) = 
              1 + (uint8) ( 0.5 + (scale * (float)(*(image+k) - immin)));
         if ( brit > 255 )
           *(image+k) = 255;
         if ( brit < 1 )
           *(image+k) = 1; 
         }
      }

   return (0);
}

int revvideo(image,len)

uint8 *image;
long len;

{
   int k;
   uint8 brit;

   for ( k=1 ; k < len ; k++ ) {
      brit = *(image+k) = 256 - *(image+k);
      if ( brit > 255 )
           *(image+k) = 255;
      if ( brit < 1 )
           *(image+k) = 1; 
      }

   return (0);
}

/************************************************************

  Main program

************************************************************/

main(argc,argv)
 
int argc;
char *argv[];
 
{

/* Input satellite data variables */

  FILE *ifile;
/*  Next line *must* be INT not LONG, as Fortran routines expect INT */
  int julian,iyear,year,month,day,hour,min,sec,i4time;
  long mci_height,mci_width,mci_depth;
  long mci_max;
  struct mc_area *area;

/* Remapping variables */

  float georange,lats,lons,line,elem;
  float lat,lon,latp,lonp;
  float fnx,fny;
  float xpix,ypix;

/*  Next line *must* be INT not LONG, as Fortran routines expect INT */
  int one=1;

/* ARPS grid variables */

  long *grid_knt;
  long *grid_accum;

  int nx_arps;
  int ny_arps;
  long grid_size;
/* INT for Fortran */
  int nfield;
  float dx_arps;
  float dy_arps;
  float xnw,ynw;

/* ARPS output variables */

  long *arps_tem;
  float *arps_grd;
  float *arps_out;
  char arunnam[81] = {81*'\0'};
  char varname[14] = {14*'\0'};
  char filvarn[7] = {7*'\0'};
  char outfname[257] = {257*'\0'};
  char satname[7];

  uint8 colors[256*3];
  uint8 *grid_val;
  uint16 status;
  char hdffname[100] = {100*'\0'};
  int dmpfmt;
  int hdf4cmpr;

/* Misc local variables */

  int cornok, ireturn;
  long i,ii,j,jj,k,kk,koff,kmci,irad,iser,jser;
  long grid_i,grid_j,grid_k;
  long sum,knt;
/*  Must be INT not LONG, as Fortran routine wants INT */
  int brite,brtmin,brtmax;
  int isource;
  long value;
  int irchan,ibird;
  uint8 abyte,bbyte;
  uint16 a2byte,b2byte;
  float filtlen = 10000.;
  float bscale;
  float fi,fj;

  long istep,jstep,iistep,jjstep,ibl,jbl,iright,jbottm;
  float istepinv,jstepinv;
  float di,dj,wnw,wne,wsw,wse;
  float xpnw,xpne,xpsw,xpse;
  float ypnw,ypne,ypsw,ypse;
  float r1nw,r1ne,r1sw,r1se;
  float r2nw,r2ne,r2sw,r2se;
  float r1pix,r2pix;


/* Calibration function prototypes */

  int vis2albedo( struct cal_gvar *, long *, float *, int);
  int ir2bright ( struct cal_gvar *, long *, float *, int, int, int);
  int cnt2bright( long *, float *, int);
  int cnt2albedo( long *, float *, int);

/************************************************************

  Begin executable code

************************************************************/


/* First check argument count and existence of satellite file 
   Default is binary */
  dmpfmt=1;
  hdf4cmpr=0;
  if( argc > 2 ) {
      if( strcmp(argv[2],"-hdf") == 0 ) {  
         dmpfmt=2;
         hdf4cmpr=0;
         if( argc > 3 ) {
             hdf4cmpr=atoi(argv[3]);
             if( hdf4cmpr < 0 || hdf4cmpr > 7 ) hdf4cmpr=0;
         }
      }
      else if ( strcmp(argv[2],"-binary") == 0 ) {
         dmpfmt=1;
         hdf4cmpr=0;
      }
      else { 
      printf(" Usage: mci2arps AREAfile_name <-hdf cmpr_lvl> <-binary> < sat.input\n");
         exit(1);
      }
  }
  if( argc < 2 ) {
      printf(" Usage: mci2arps AREAfile_name <-hdf cmpr_lvl> <-binary> < sat.input\n");
      exit(1);
  }

  if( dmpfmt == 1 ) {
    printf(" Creating binary remapped file\n");
  }
  else {
    printf(" Creating hdf remapped file with compression level: %i \n",hdf4cmpr);
  }

  if( (ifile = fopen(argv[1],"r")) == NULL ) {
      printf("Unable to open McIDAS AREA file %s\n",argv[1] );
      exit(1);
      }

/* Read-in ARPS-grid input parameters */

  inisatarps(arunnam,&nx_arps,&ny_arps,
              &dx_arps,&dy_arps,&xnw,&ynw);
  printf(" Back from inisatarps runname:%s\n",arunnam);
  printf("   nx:%i, ny:%i, dx:%f, dy:%f\n",
           nx_arps,ny_arps,dx_arps,dy_arps);
  printf("   xnw:%f, ynw:%f\n",xnw,ynw);

/* Read-in satellite data from McIDAS AREA file */

  if( AREA_read( ifile,&area ) == -1 ){
      printf("Unable to read McIDAS AREA file %s\n",argv[1] );
      exit(2);
  }

  printf( " Calling sat_init area->type is %s\n",area->type );
  sat_init( area );

/*
   Decode and display header information
*/

  julian = area->dir->ndate % 1000;
 
  if(( area->dir->ndate / 1000 ) % 4 == 0 ){
      for( month = 0; month < 12 && jdayl[month] < julian; month++ );
      day = julian - jdayl[month-1];
      }
  else {
      for( month = 0; month < 12 && jday[month] < julian; month++ );
      day = julian - jday[month-1];
      }
  year=area->dir->ndate/1000;
  hour=area->dir->ntime/10000;
  min=(area->dir->ntime/100) % 100;
  sec=0;
  iyear=year+1900;
  if(iyear < 1960)
     iyear=year+2000;
  ctim2abss(&iyear,&month,&day,&hour,&min,&sec,&i4time);
  sprintf( datim,"%4.4dZ %2d %3s 19%2.2d", (int)(area->dir->ntime/100),
      (int)day,mnth[month-1],(int)(area->dir->ndate/1000) );
  if( out_string && !strcmp( out_string,"date" )){
      printf( "%s\n", datim );
      exit( WNOERR );
      }
  else if( out_string && !strcmp( out_string,"file" )){
         printf( "%2.2i%2.2i%2.2i%2.2i\n", (int)(area->dir->ndate/1000),
            (int)month,(int)day,(int)area->dir->ntime/10000 );
      }
 
  printf( "Date: %i %s %i\n", (int)area->dir->ndate, mnth[month-1], (int)day);
  printf( "%2.2i%2.2i%2.2i%04.4i\n", (int)area->dir->ndate/1000,
      (int)month,(int)day,(int)area->dir->ntime/100 );
  printf( "Size: %ix%i\n", (int)area->dir->esiz, (int)area->dir->lsiz );
  printf( "Data element size: %i bytes\n",(int)area->dir->zsiz );
  printf( "Resolution: %ix%i\n", (int)area->dir->eres, (int)area->dir->lres );
  printf( "Correction: %ix%i\n", (int)area->dir->ecor, (int)area->dir->lcor );
  printf( "Bands: %d\n", (int)area->dir->bands );
  printf( "Band indicator: %i\n", (int)area->dir->filtmap );
  printf( "Prefix length %i\n", (int)area->dir->yzprefix );
  printf( "Nav  offset %i\n", (int)area->dir->navoffst );
  printf( "Cal  offset %i\n", (int)area->dir->caloffst );
  printf( "Data offset %i\n", (int)area->dir->datoffst );
  printf( "Source type: %4.4s\n", area->dir->stype );
  printf( "Calibration type: %4.4s\n", area->dir->ctype );
  printf( "Navigation type: %4.4s\n", area->type );
  printf( "Comments: %32.32s\n", area->dir->comments );
  printf( "\n  startscan: %i\n", (int)area->dir->startscan );
  printf( "  doclen   : %i\n", (int)area->dir->doclen );
  printf( "  callen   : %i\n", (int)area->dir->callen );
  printf( "  levlen   : %i\n", (int)area->dir->levlen );
  printf( "  satid    : %d\n", (int)area->dir->satid );
 
/* Determine bird number from satellite ID number */
 
 
 if (area->dir->satid == G06IR || area->dir->satid == G06VIS)
     ibird=6;
  else if (area->dir->satid == G07IR || area->dir->satid == G07VIS)
     ibird=7;
  else if (area->dir->satid == G08IM || area->dir->satid == G08SND)
     ibird=8;
  else if (area->dir->satid == G09IM || area->dir->satid == G09SND)
     ibird=9;
  else if (area->dir->satid == G10IM || area->dir->satid == G10SND)
     ibird=10;
  else if (area->dir->satid == G11IM || area->dir->satid == G11SND)
     ibird=11;
  else if (area->dir->satid == G12IM || area->dir->satid == G12SND)
     ibird=12;
  else if (area->dir->satid == G13IM || area->dir->satid == G13SND)
     ibird=13;
  else if (area->dir->satid == G14IM || area->dir->satid == G14SND)
     ibird=14;
  else if (area->dir->satid == G15IM || area->dir->satid == G15SND)
     ibird=15;
  printf(" Satellite is GOES%2d\n",ibird);
  sprintf(satname,"goes%2.2d",ibird);

 
/* Check bit set in filtmap and report results 
   Assuming here 1 band per file 
   irchan assigned here is an index in the ir calibration arrays */

  if (( 1 & area->dir->filtmap) != 0 ) {
      irchan=0;
      printf( "\n File contains visible data\n\n");
      }
  else if (( 2 & area->dir->filtmap) != 0 ) {
      irchan=2;
      printf( "\n File contains Channel 2 (3.9 micron) data\n\n");
      }
  else if (( 4 & area->dir->filtmap) != 0 ) {
      irchan=3;
      printf( "\n File contains Channel 3 (6.8 micron) data\n\n");
      }
  else if (( 8 & area->dir->filtmap) != 0 ) {
      irchan=4;
      printf( "\n File contains Channel 4 (10.7 micron) data\n\n");
      }
  else if ((16 & area->dir->filtmap) != 0 ) {
      irchan=5;
      printf( "\n File contains Channel 5 (12.0 micron) data\n\n");
      }

  mci_height = area->dir->lsiz;
  mci_width  = area->dir->esiz;
  mci_depth = area->dir->zsiz;
  // ENatenberg 07-15-13 NOAA GOES CLASS DATA SUPPORT 
  // mci_max = mci_height * mci_width * mci_depth;
  mci_max = mci_height * mci_width * mci_depth + ((long)(area->dir->yzprefix) *  mci_height);

  printf(" Remapping data\n");

/* 
  Determine the satellite sub-point
*/

  georange = GEORNG;
  if( !strncmp( area->type,"GOES",4 ))
      sat_domain.proj = PROJ_SAT;
  else if( !strncmp( area->type,"GVAR",4 ))
      sat_domain.proj = PROJ_SAT;
  else if( !strncmp( area->type,"PS",2 )){
      sat_domain.proj = PROJ_PS;
      lats = area->navps->stdlat/10000.;
      lons = -area->navps->nrmlon/10000.;
      }

  if( sat_domain.proj == PROJ_SAT )
      sat_tran( 3,&lats,&lons,&line,&elem );

  sat_domain.lon[1] = lons;
  sat_domain.lat[1] = lats;
  mapcoord( lats,lons );

  fny = (float)mci_height;
  fnx = (float)mci_width; 
  sat_tran( 2,&lats,&lons,&line,&elem );
  printf(" Sub-satellite line: %10.2f ,elem: %10.2f\n",line,elem);

  sat_line_elem( 1, 1.0, 1.0,&line,&elem);
  sat_tran( 1,&lat,&lon,&line,&elem);
  printf("         1, 1:  lat: %10.2f   lon: %10.2f\n",lat,lon);

  sat_line_elem( 1, fny, 1.0,&line,&elem);
  sat_tran( 1,&lat,&lon,&line,&elem);
  printf("     height,1:  lat: %10.2f   lon: %10.2f\n",lat,lon);
 
  sat_line_elem( 1, 1.0, fnx,&line,&elem);
  sat_tran( 1,&lat,&lon,&line,&elem);
  printf("      1,width:  lat: %10.2f   lon: %10.2f\n",lat,lon);
 
  sat_line_elem( 1, fny, fnx,&line,&elem);
  sat_tran( 1,&lat,&lon,&line,&elem);
  printf(" height,width:  lat: %10.2f   lon: %10.2f\n",lat,lon);

/* Set up variables in solar angle correction routine */

  solcorset(&i4time,&lats,&lons,&georange);

/* Allocate memory for ARPS output */

  grid_size=nx_arps*ny_arps;
  arps_tem = (long *)malloc(sizeof(long)*grid_size );
  arps_grd = (float *)malloc(sizeof(float)*grid_size ); 

/* Allocate memory for grid computations and image */

  grid_knt   = (long *)malloc(sizeof(long)*grid_size );
  grid_accum = (long *)malloc(sizeof(long)*grid_size );
  grid_val = (uint8 *)malloc(sizeof(uint8)*grid_size );

  if ( grid_val != NULL ) {

     for ( k = 0; k < grid_size ; k++ ){
        *(grid_knt+k)   = 0;
        *(grid_accum+k) = 0;
       }

 
     istep=32/area->dir->eres;
     if ( istep < 5 ) 
       istep=5;
     jstep=32/area->dir->lres;
     if ( jstep < 5 ) 
       jstep=5;
     printf("\n Block size for linear pixel x,y calc: %i X %i pixels\n\n",
             (int)istep,(int)jstep);

     printf(" Primary remapping\n");
     brtmin = 999;
     brtmax = -999; 
     if ( mci_depth == 1) {
        bscale=1.0;
        isource = 2;
        //koff  = 0;
        // ENatenberg 07-15-13 NOAA GOES CLASS DATA SUPPORT
        koff = area->dir->yzprefix;
        printf( "  1-byte size, koff = %i\n", (int)koff );

        for (jbl=0; jbl < mci_height; jbl+=jstep) {
           for (ibl=0; ibl < mci_width; ibl+=istep) {
              cornok=0;
              iright=(ibl+istep-1);
              if (iright > mci_width)
                iright = mci_width;

              jbottm=(jbl+jstep-1);
              if (jbottm > mci_height)
                jbottm = mci_height;
/*
Find lat,lon,x,y,r1,r2 at each corner point
Northwest corner:
 */
              fi=(float)ibl;
              fj=(float)jbl;
              sat_line_elem( 1, fj, fi,&line,&elem);
              sat_tran( 1,&latp,&lonp,&line,&elem);
              if ( ( latp > -90. ) &&  ( lonp > -900. ) ) {
                 lltoxy(&one,&one,&latp,&lonp,&xpnw,&ypnw);
                 solr1r2(&i4time,&lons,&latp,&lonp,&r1nw,&r2nw);
              }
              else
                cornok=-1;
/*
Find lat,lon,x,y,r1,r2 at each corner point
Northeast corner:
 */
              fi=(float)iright;
              fj=(float)jbl;
              sat_line_elem( 1, fj, fi,&line,&elem);
              sat_tran( 1,&latp,&lonp,&line,&elem);
              if ( ( latp > -90. ) &&  ( lonp > -900. ) ) {
                 lltoxy(&one,&one,&latp,&lonp,&xpne,&ypne);
                 solr1r2(&i4time,&lons,&latp,&lonp,&r1ne,&r2ne);
              }
              else
                cornok=-1;
/*
Find lat,lon,x,y,r1,r2 at each corner point
Southwest corner:
 */
              fi=(float)ibl;
              fj=(float)jbottm;
              sat_line_elem( 1, fj, fi,&line,&elem);
              sat_tran( 1,&latp,&lonp,&line,&elem);
              if ( ( latp > -90. ) &&  ( lonp > -900. ) ) {
                 lltoxy(&one,&one,&latp,&lonp,&xpsw,&ypsw);
                 solr1r2(&i4time,&lons,&latp,&lonp,&r1sw,&r2sw);
              }
              else
                cornok=-1;
/*
Find lat,lon,x,y,r1,r2 at each corner point
Southeast corner:
 */
              fi=(float)iright;
              fj=(float)jbottm;
              sat_line_elem( 1, fj, fi,&line,&elem);
              sat_tran( 1,&latp,&lonp,&line,&elem);
              if ( ( latp > -90. ) &&  ( lonp > -900. ) ) {
                 lltoxy(&one,&one,&latp,&lonp,&xpse,&ypse);
                 solr1r2(&i4time,&lons,&latp,&lonp,&r1se,&r2se);
              }
              else
                cornok=-1;
 
/*
Loop through all pixels in block
 */
              if ( cornok == 0 ) {

              if ((iistep = iright - ibl) > 0)
                istepinv=1./(float)iistep;
              else
                istepinv=1.;
              
              if ((jjstep = jbottm - jbl) > 0)
                jstepinv=1./(float)jjstep;
              else
                jstepinv=1.;
 
              for (jj=0; jj <= jjstep; jj++) {
                 for (ii=0; ii <= iistep; ii++) {
                    i=ibl+ii;
                    j=jbl+jj;
                    di=ii*istepinv;
                    dj=jj*jstepinv;
                    wnw=(1.-di)*(1.-dj);
                    wne=di*(1.-dj);
                    wsw=(1.-di)*dj;
                    wse=di*dj;
                    xpix=wnw*xpnw + wne*xpne + wsw*xpsw + wse*xpse;
                    ypix=wnw*ypnw + wne*ypne + wsw*ypsw + wse*ypse;
                    r1pix=wnw*r1nw + wne*r1ne + wsw*r1sw + wse*r1se;
                    r2pix=wnw*r2nw + wne*r2ne + wsw*r2sw + wse*r2se;
 
                    k = i + (j*mci_width);
                    //kmci = k + koff;
                    // ENatenberg 07-15-13 NOAA GOES CLASS DATA SUPPORT
                    kmci = 1*k + (j*area->dir->yzprefix) + koff;

		    if ( kmci < mci_max )
			value = *(area->image)+kmci;
		    else
			value = 0;

		    if ( value != 0 ) {
                       brite=(int)*((area->image)+kmci);
                       if ( brite < brtmin ) brtmin = brite;
                       if ( brite > brtmax ) brtmax = brite;

                       if ( ( latp > -90. ) &&  ( lonp > -900. ) ) {
                          grid_i=(long)(((xpix - xnw)/dx_arps)+0.5);
                          grid_j=(long)(((ynw - ypix)/dy_arps)+0.5);
/*    printf( " grid_i: %d  grid_j: %d\n",grid_i,grid_j); */
                          if ( ( grid_i > -1 ) && ( grid_i < nx_arps ) &&
                               ( grid_j > -1 ) && ( grid_j < ny_arps) ){
                             grid_k = grid_i + (grid_j*nx_arps);
/*
                             if ( k%10000 == 0 ) {
    printf(" k: %u image value: %u\n",(k/10000),*((area->image)+kmci));
    printf(" xpix:%f   ypix:%f   grid_i:%i   grid_j:%i\n\n",
           (xpix*0.001),(ypix*0.001),grid_i,grid_j);
                             }
*/
                             *(grid_knt+grid_k) = *(grid_knt+grid_k) + 1;
                              if ( irchan == 0 )
                              solrsc1(&brite,&r1pix,&r2pix,&brite);
                             *(grid_accum+grid_k) = *(grid_accum+grid_k)+brite;
                             }
                          }
                       }
                    }
                 }
                 }
              }
           }
        }
     else {
        isource = 1;
        bscale=255./1023.;
        koff = area->dir->yzprefix;
        printf( "  2-byte size, koff = %i\n", (int)koff );
        for (jbl=0; jbl < mci_height; jbl+=jstep) {
           for (ibl=0; ibl < mci_width; ibl+=istep) {
              cornok=0;

              iright=(ibl+istep-1);
              if (iright > mci_width)
                iright = mci_width;

              jbottm=(jbl+jstep-1);
              if (jbottm > mci_height)
                jbottm = mci_height;
/*
Find lat,lon,x,y,r1,r2 at each corner point
Northwest corner:
 */
              fi=(float)ibl;
              fj=(float)jbl;
              sat_line_elem( 1, fj, fi,&line,&elem);
              sat_tran( 1,&latp,&lonp,&line,&elem);
              if ( ( latp > -90. ) &&  ( lonp > -900. ) ) {
                 lltoxy(&one,&one,&latp,&lonp,&xpnw,&ypnw);
                 solr1r2(&i4time,&lons,&latp,&lonp,&r1nw,&r2nw);
              }
              else
                 cornok=-1;
/*
Find lat,lon,x,y,r1,r2 at each corner point
Northeast corner:
 */
              fi=(float)(iright);
              fj=(float)jbl;
              sat_line_elem( 1, fj, fi,&line,&elem);
              sat_tran( 1,&latp,&lonp,&line,&elem);
              if ( ( latp > -90. ) &&  ( lonp > -900. ) ) {
                 lltoxy(&one,&one,&latp,&lonp,&xpne,&ypne);
                 solr1r2(&i4time,&lons,&latp,&lonp,&r1ne,&r2ne);
              }
              else
                 cornok=-1;
/*
Find lat,lon,x,y,r1,r2 at each corner point
Southwest corner:
 */
              fi=(float)ibl;
              fj=(float)jbottm;
              sat_line_elem( 1, fj, fi,&line,&elem);
              sat_tran( 1,&latp,&lonp,&line,&elem);
              if ( ( latp > -90. ) &&  ( lonp > -900. ) ) {
                 lltoxy(&one,&one,&latp,&lonp,&xpsw,&ypsw);
                 solr1r2(&i4time,&lons,&latp,&lonp,&r1sw,&r2sw);
              }
              else
                 cornok=-1;
/*
Find lat,lon,x,y,r1,r2 at each corner point
Southeast corner:
 */
              fi=(float)iright;
              fj=(float)jbottm;
              sat_line_elem( 1, fj, fi,&line,&elem);
              sat_tran( 1,&latp,&lonp,&line,&elem);
              if ( ( latp > -90. ) &&  ( lonp > -900. ) ) {
                 lltoxy(&one,&one,&latp,&lonp,&xpse,&ypse);
                 solr1r2(&i4time,&lons,&latp,&lonp,&r1se,&r2se);
              }
              else
                 cornok=-1;
 
/*
Loop through all pixels in block
 */
              if( cornok == 0 ) {

              if ((iistep = iright - ibl) > 0)
                istepinv=1./(float)iistep;
              else
                istepinv=1.;
              
              if ((jjstep = jbottm - jbl) > 0)
                jstepinv=1./(float)jjstep;
              else
                jstepinv=1.;

              for (jj=0; jj <= jjstep; jj++) {
                 for (ii=0; ii <= iistep; ii++) {
                    i=ibl+ii;
                    j=jbl+jj;
                    di=ii*istepinv;
                    dj=jj*jstepinv;
                    wnw=(1.-di)*(1.-dj);
                    wne=di*(1.-dj);
                    wsw=(1.-di)*dj;
                    wse=di*dj;

                    k = i + (j*mci_width);
                    kmci = 2*k + (j*area->dir->yzprefix) + koff;

/*
 *  Make sure we don't bust an array.
 */
		    if ( kmci + 1 >= mci_max )
			continue;

/*                  printf(" k: %u image value: %d\n",k,*((area->image)+kmci)); */
                    abyte = *((area->image)+kmci);
                    bbyte = *((area->image)+kmci+1);
                    a2byte = (uint16)abyte;
                    b2byte = (uint16)bbyte;

/*  Calculate total value normalize to 256 and store as long */

                    brite = (int) ((a2byte<<3) + (b2byte>>5));
/*
                    if( ( kmci % 1000) == 0 ) {
                        printf( " a byte = %u   b byte = %u\n",abyte,bbyte); 
                        printf( " brite = %d\n",brite);
                        }
*/
                    if( brite != 0 ) {
                       if ( brite < brtmin ) brtmin = brite;
                       if ( brite > brtmax ) brtmax = brite;
                       xpix=wnw*xpnw + wne*xpne + wsw*xpsw + wse*xpse;
                       ypix=wnw*ypnw + wne*ypne + wsw*ypsw + wse*ypse;
                       grid_i=(long)(((xpix - xnw)/dx_arps)+0.5);
                       grid_j=(long)(((ynw - ypix)/dy_arps)+0.5);
                       if ( ( grid_i > -1 ) && ( grid_i < nx_arps ) &&
                            ( grid_j > -1 ) && ( grid_j < ny_arps) ){

/*                           if ( k%10000 == 0 ) {
    printf(" k: %u image value: %u\n",(k/10000),*((area->image)+kmci));
    printf(" xpix:%f   ypix:%f   grid_i:%i   grid_j:%i\n\n",
           (xpix*0.001),(ypix*0.001),grid_i,grid_j);
                             }
 */
                          r1pix=wnw*r1nw + wne*r1ne + wsw*r1sw + wse*r1se;
                          r2pix=wnw*r2nw + wne*r2ne + wsw*r2sw + wse*r2se;
                          solrsc2(&brite,&r1pix,&r2pix,&brite);
                          grid_k = grid_i + (grid_j*nx_arps);
                          *(grid_knt+grid_k) = *(grid_knt+grid_k) + 1;
                          *(grid_accum+grid_k) = *(grid_accum+grid_k)+brite;

                          }
                       }
                    }
                 }
                 }
              }
           }
        }

/* Print brightness stats */
     printf(" Brightness count range for entire file:\n");
     printf("    min=%ld     max=%ld\n",brtmin,brtmax);

/* Find average in each cell */

     printf(" Average calc\n");
     for ( k = 0; k < grid_size ; k++ ){
/*
           if ( k%500 == 0 )
    printf(" k: %i accum: %i  knt: %i\n",k,*(grid_accum+k),*(grid_knt+k));
*/
        if ( *(grid_accum+k) > 0 ) {

          *(grid_val+k) = (uint8) ( 0.5 + 
             bscale*( ((float)*(grid_accum+k)) /
                      ((float)*(grid_knt+k)) ) );
          *(arps_tem+k) = (long) ( 0.5 + ( ((float)*(grid_accum+k)) /
                                           ((float)*(grid_knt+k)) ) );
           }
        else {
           *(grid_val+k) = 0;
           *(arps_tem+k) = 0;
           }
        }

/* Fill any missing using average of neighbors */

     printf(" Hole Filling\n");
     for ( j = 0; j < ny_arps ; j++ ) {
        for ( i = 0; i < nx_arps ; i++ ) {
           k = i + (j*nx_arps);
           if ( *(grid_knt+k) == 0 ) {
              for ( irad = 1; irad < 5 ; irad++) {
                 sum = 0;
                 knt = 0;
                 for ( ii = -irad; ii < (irad+1) ; ii++ ) {
                    iser = i+ii;
                    if ( iser > -1 && iser < nx_arps ) {
                       for ( jj = -irad; jj < (irad+1) ; jj++) {
                          jser = j+jj;
                          if( jser > -1 && jser < ny_arps ) {
                             kk = iser + (jser*nx_arps);
                             sum = sum + *(grid_accum+kk);
                             knt = knt + *(grid_knt+kk);
                             }
                          }
                       }
                    }
                 if ( knt > 0 ) {
                    *(grid_val+k) = (uint8)(0.5+
                             bscale*(((float)sum)/((float)knt)));
                    *(arps_tem+k) = (long)(0.5+(((float)sum)/((float)knt)));
                    break;
                    }
                 }
              }
           }
        }
  
     free(grid_accum);
     free(grid_knt);

     printf(" Calibrating data for ARPS input file\n");
     if ( mci_depth == 1) {
        if ( irchan == 0 ) {
            printf(" Converting visible counts to albedo\n");
            ireturn = cnt2albedo( arps_tem, arps_grd, grid_size);
            strcpy(varname,"satalb");
            strcpy(filvarn,"satalb");
            }
        else {
            printf(" Converting IR counts to brightness temp\n");
            ireturn = cnt2bright( arps_tem, arps_grd, grid_size);
            strcpy(varname,"satctt");
            strcpy(filvarn,"satctt");
            }
        }
     else {
        if ( irchan == 0 ) {
            printf(" Converting raw visible image to albedo\n");
	    if ( area->calgv == NULL ) {
    		area->calgv =
          	(struct cal_gvar *)malloc(sizeof(struct cal_gvar));
	    }
            ireturn = vis2albedo( area->calgv, arps_tem, arps_grd,
                              grid_size);
            strcpy(varname,"satalb");
            strcpy(filvarn,"satalb");
            }
        else {
            printf(" Converting raw IR image to brightness temp\n");
	    if ( area->calgv == NULL ) {
    		area->calgv =
          	(struct cal_gvar *)malloc(sizeof(struct cal_gvar));
	    }
            ireturn = ir2bright( area->calgv, arps_tem, arps_grd, 
                             ibird, irchan, grid_size);
            strcpy(varname,"satctt");
            strcpy(filvarn,"satctt");
            }
        }

     if( irchan == 4) {
       nfield=2;
       strcpy(varname+6,"satcft");
       arps_out = (float *)malloc( 2*sizeof(float)*grid_size ); 
       ireturn = yflip(nx_arps,ny_arps,arps_grd,arps_out);
       coldfilt(&nx_arps,&ny_arps,&filtlen,
                arps_out,(arps_out+grid_size));
     } else {
       nfield=1;
       arps_out = (float *)malloc(sizeof(float)*grid_size ); 
       ireturn = yflip(nx_arps,ny_arps,arps_grd,arps_out);
     }

     if(dmpfmt==1)
/*
     sprintf( outfname, "%s.%02i%02i%02i%02i%02i.goes%02i.%s",
              arunnam,(int)year,(int)month,(int)day,
              (int)hour,(int)min,ibird,varname );
 */
     sprintf( outfname, "%s_%s.%s", arunnam, satname, filvarn );
     else
/*
     sprintf( outfname, "%s.%02i%02i%02i%02i%02i.goes%02i.%s.hdf4",
              arunnam,(int)year,(int)month,(int)day,
              (int)hour,(int)min,ibird,varname );
 */
     sprintf( outfname, "%s_%s.%s.hdf4", arunnam, satname, filvarn );
     printf( " Output file name: %s\n",outfname);
/*
     sprintf( hdffname, "%s.%02i%02i%02i%02i%02i.goes%02i.%s.hdf",
              arunnam,(int)year,(int)month,(int)day,
              (int)hour,(int)min,ibird,varname );
 */
     sprintf( hdffname, "%s_%s.%s.hdf", arunnam, satname, filvarn );
     printf( " hdf image file name: %s\n",hdffname);

/*   Call FORTRAN routine to write ARPS gridded data */

     printf(" Calling FORTRAN write routine 1\n");
     wtsatfld(&nx_arps,&ny_arps,&nfield,
               outfname,satname,&lats,&lons,
               &year,&month,&day,&hour,&min,&sec,&isource,
               &dmpfmt,&hdf4cmpr,
               varname,arps_out);
     printf(" Back from FORTRAN write routine\n");

/*   Rescale image data for a better greyscale depiction

     printf(" Rescaling image\n");
     status = rescale( grid_val, grid_size);
*/

/*   Get the palette data from file. */

      status = DFPgetpal("palgrey.hdf",colors);
      printf("DFPgetpal return status: %d\n",status);

/*   Set the current palette. */

      status = DFR8setpalette(colors);
      printf("DFR8setpalette return status: %d\n",status);
 
/*   Write the image data to the file. */

      status=DFR8putimage(hdffname,grid_val,nx_arps,ny_arps,0);
      printf("DFR8putimage return status: %d\n",status);

      }
   else {
      printf(" Error allocating memory");
      exit(1);
      }
   
  exit(0);

}
