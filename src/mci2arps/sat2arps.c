/*  wdt Copyright (c) 2001 Weather Decision Technologies, Inc. 2001-03-23 GMB
*
*   unisat2arps.c
*
*   Program to read NOAAPORT satellite images (in raw, 8-bit format) 
*   and remap to ARPS grid.
*
*   The program also writes an hdf image file for diagnostics.
*
*   This program was adapted from mci2arps written by:
*     Keith Brewster
*     CAPS/Univ of Oklahoma
*     September, 1997
*
*   WARNING - This code is untested.
*   Updated to handle 64-bit issues.  Changes are same as in "mci2arps.c".
*   Kevin W. Thomas
*   August, 2004.
*
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <malloc.h>
#include <unistd.h>

#include "hdf.h"

/* Default satellite image information */

/*int sat_nx = 5120, sat_ny = 5120;*/ /* Vis, GOES East */
/*int sat_nx = 4400, sat_ny = 5120;*/ /* Vis, GOES West */
int sat_nx = 1280, sat_ny = 1280; /* IR, GOES East */
/*int sat_nx = 1100, sat_ny = 1280;*/ /* IR, GOES West */

/*float sat_res_m = 1015.9;*/ /* Vis */
float sat_res_m = 4063.5; /* IR */

int sat_mapproj = 2;
float sat_scale = 1.;
float sat_ctrlat = 40.538, sat_ctrlon = -87.597; /* GOES East */
/*float sat_ctrlat = 39.256, sat_ctrlon = -117.485;*/ /* GOES West */
float sat_trulat[2] = {25,25};
float sat_trulon = -95.;

int fill_flg = 0; /* no hole filling - 0, fill missing - 1 */

float sat_xnw, sat_ynw;

#define NUM_TBLS 6

int sat_indx = 3;
     /* sat_indx=0  NESDIS IR -> cloud top temperature (equation)
      * sat_indx=1  Unisys IR -> cloud top temperature 
      * sat_indx=2  NESDIS Vis -> albedo (equation) 
      */

int map_type[NUM_TBLS] = { 0, 0, 1}; /* map_type 0-IR, 1-Vis */
          /* for sat_indx: 0  1  2 */

int mxindx[NUM_TBLS]   = { 254, 200, 254, };
          /* for sat_indx:   0    1    2  */

int sat_table[NUM_TBLS][256] = {
      /* Table 0 */
      /* NESDIS IR -> cloud top temperature 
       * from Howard Carney at NESDIS 
       *    (2 Oct 2001 email "Format of NOAAPORT satellite data")
       * Count range 176 to 255, T = 418 - C 
       * Count range 0 to 176,   T = (660 - C)/2   */
     -999.00, 329.50, 329.00, 328.50, 328.00, 327.50, 327.00, 326.50,
      326.00, 325.50, 325.00, 324.50, 324.00, 323.50, 323.00, 322.50,
      322.00, 321.50, 321.00, 320.50, 320.00, 319.50, 319.00, 318.50,
      318.00, 317.50, 317.00, 316.50, 316.00, 315.50, 315.00, 314.50,
      314.00, 313.50, 313.00, 312.50, 312.00, 311.50, 311.00, 310.50,
      310.00, 309.50, 309.00, 308.50, 308.00, 307.50, 307.00, 306.50,
      306.00, 305.50, 305.00, 304.50, 304.00, 303.50, 303.00, 302.50,
      302.00, 301.50, 301.00, 300.50, 300.00, 299.50, 299.00, 298.50,
      298.00, 297.50, 297.00, 296.50, 296.00, 295.50, 295.00, 294.50,
      294.00, 293.50, 293.00, 292.50, 292.00, 291.50, 291.00, 290.50,
      290.00, 289.50, 289.00, 288.50, 288.00, 287.50, 287.00, 286.50,
      286.00, 285.50, 285.00, 284.50, 284.00, 283.50, 283.00, 282.50,
      282.00, 281.50, 281.00, 280.50, 280.00, 279.50, 279.00, 278.50,
      278.00, 277.50, 277.00, 276.50, 276.00, 275.50, 275.00, 274.50,
      274.00, 273.50, 273.00, 272.50, 272.00, 271.50, 271.00, 270.50,
      270.00, 269.50, 269.00, 268.50, 268.00, 267.50, 267.00, 266.50,
      266.00, 265.50, 265.00, 264.50, 264.00, 263.50, 263.00, 262.50,
      262.00, 261.50, 261.00, 260.50, 260.00, 259.50, 259.00, 258.50,
      258.00, 257.50, 257.00, 256.50, 256.00, 255.50, 255.00, 254.50,
      254.00, 253.50, 253.00, 252.50, 252.00, 251.50, 251.00, 250.50,
      250.00, 249.50, 249.00, 248.50, 248.00, 247.50, 247.00, 246.50,
      246.00, 245.50, 245.00, 244.50, 244.00, 243.50, 243.00, 242.50,
      242.00, 241.00, 240.00, 239.00, 238.00, 237.00, 236.00, 235.00,
      234.00, 233.00, 232.00, 231.00, 230.00, 229.00, 228.00, 227.00,
      226.00, 225.00, 224.00, 223.00, 222.00, 221.00, 220.00, 219.00,
      218.00, 217.00, 216.00, 215.00, 214.00, 213.00, 212.00, 211.00,
      210.00, 209.00, 208.00, 207.00, 206.00, 205.00, 204.00, 203.00,
      202.00, 201.00, 200.00, 199.00, 198.00, 197.00, 196.00, 195.00,
      194.00, 193.00, 192.00, 191.00, 190.00, 189.00, 188.00, 187.00,
      186.00, 185.00, 184.00, 183.00, 182.00, 181.00, 180.00, 179.00,
      178.00, 177.00, 176.00, 175.00, 174.00, 173.00, 172.00, 171.00,
      170.00, 169.00, 168.00, 167.00, 166.00, 165.00, 164.00,-999.00,

      /* Table 1 */
      /* Unisys IR image to Brightness Temperature table (from unisat2arps.c).
       * Derived by comparing mci2arps output with USIRLAMB images. */
      -999.00, 325.68, 325.13, 324.57, 324.01,
       323.45, 322.90, 322.34, 321.78, 321.23,
       320.67, 320.11, 319.56, 319.00, 318.44,
       317.88, 317.33, 316.77, 316.21, 315.66,
       315.10, 314.54, 313.99, 313.43, 312.87,
       312.31, 311.76, 311.20, 310.64, 310.09,
       309.53, 308.97, 308.42, 307.86, 307.30,
       306.75, 306.19, 305.63, 305.07, 304.52,
       303.96, 303.40, 302.85, 302.29, 301.73,
       301.18, 300.62, 300.06, 299.50, 298.95,
       298.39, 297.83, 297.28, 296.72, 296.16,
       295.61, 295.05, 294.49, 293.93, 293.38,
       292.82, 292.32, 291.79, 291.29, 290.74,
       290.08, 289.54, 289.03, 288.44, 287.91,
       287.35, 286.79, 286.22, 285.65, 285.04,
       284.50, 283.87, 283.25, 282.77, 282.29,
       281.74, 281.18, 280.66, 280.18, 279.76,
       279.19, 278.68, 278.09, 277.52, 277.07,
       276.57, 276.14, 275.60, 275.11, 274.56,
       274.10, 273.54, 273.01, 272.45, 271.98,
       271.43, 270.90, 270.28, 269.74, 269.04,
       268.74, 268.29, 267.69, 267.24, 266.74,
       266.18, 265.65, 265.27, 264.69, 264.16,
       263.75, 263.28, 262.76, 262.29, 261.94,
       261.34, 260.88, 260.43, 259.92, 259.44,
       258.96, 258.53, 258.13, 257.64, 257.23,
       256.75, 256.25, 255.81, 255.28, 254.83,
       254.34, 253.91, 253.40, 253.06, 252.64,
       252.07, 251.57, 251.07, 250.73, 250.19,
       249.77, 249.27, 248.67, 248.46, 247.91,
       247.41, 247.01, 246.32, 245.88, 245.53,
       245.01, 244.49, 243.76, 243.41, 242.71,
       242.38, 241.76, 241.24, 240.63, 240.17,
       239.58, 239.06, 238.45, 237.89, 237.33,
       236.86, 236.24, 235.63, 235.14, 234.36,
       233.87, 232.94, 232.07, 231.00, 230.00,
       229.02, 227.89, 226.88, 225.89, 224.88,
       223.90, 222.85, 222.04, 220.92, 219.96,
       218.86, 217.82, 216.78, 215.74, 214.55,
       213.53, 213.04, 212.36, 211.39, 210.44,
       209.48,-999.00,-999.00,-999.00,-999.00,
      -999.00,-999.00,-999.00,-999.00,-999.00,
      -999.00,-999.00,-999.00,-999.00,-999.00,
      -999.00,-999.00,-999.00,-999.00,-999.00,
      -999.00,-999.00,-999.00,-999.00,-999.00,
      -999.00,-999.00,-999.00,-999.00,-999.00,
      -999.00,-999.00,-999.00,-999.00,-999.00,
      -999.00,-999.00,-999.00,-999.00,-999.00,
      -999.00,-999.00,-999.00,-999.00,-999.00,
      -999.00,-999.00,-999.00,-999.00,-999.00,
      -999.00,-999.00,-999.00,-999.00,-999.00,
      -999.00,

      /* Table 2 */
      /* NESDIS Vis 8-bit Brightness Albedo table: 
       * from Howard Carney at NESDIS 
       *    (2 Oct 2001 email "Format of NOAAPORT satellite data")
       * Albedo = (C/255)**2   */
      -999.000, 0.000015, 0.000062, 0.000138, 0.000246, 0.000384, 0.000554,
      0.000754, 0.000984, 0.001246, 0.001538, 0.001861, 0.002215, 0.002599,
      0.003014, 0.003460, 0.003937, 0.004444, 0.004983, 0.005552, 0.006151,
      0.006782, 0.007443, 0.008135, 0.008858, 0.009612, 0.010396, 0.011211,
      0.012057, 0.012933, 0.013841, 0.014779, 0.015748, 0.016747, 0.017778,
      0.018839, 0.019931, 0.021053, 0.022207, 0.023391, 0.024606, 0.025852,
      0.027128, 0.028435, 0.029773, 0.031142, 0.032541, 0.033972, 0.035433,
      0.036924, 0.038447, 0.040000, 0.041584, 0.043199, 0.044844, 0.046521,
      0.048228, 0.049965, 0.051734, 0.053533, 0.055363, 0.057224, 0.059116,
      0.061038, 0.062991, 0.064975, 0.066990, 0.069035, 0.071111, 0.073218,
      0.075356, 0.077524, 0.079723, 0.081953, 0.084214, 0.086505, 0.088827,
      0.091180, 0.093564, 0.095978, 0.098424, 0.100900, 0.103406, 0.105944,
      0.108512, 0.111111, 0.113741, 0.116401, 0.119093, 0.121815, 0.124567,
      0.127351, 0.130165, 0.133010, 0.135886, 0.138793, 0.141730, 0.144698,
      0.147697, 0.150727, 0.153787, 0.156878, 0.160000, 0.163153, 0.166336,
      0.169550, 0.172795, 0.176071, 0.179377, 0.182714, 0.186082, 0.189481,
      0.192910, 0.196371, 0.199862, 0.203383, 0.206936, 0.210519, 0.214133,
      0.217778, 0.221453, 0.225160, 0.228897, 0.232664, 0.236463, 0.240292,
      0.244152, 0.248043, 0.251965, 0.255917, 0.259900, 0.263914, 0.267958,
      0.272034, 0.276140, 0.280277, 0.284444, 0.288643, 0.292872, 0.297132,
      0.301423, 0.305744, 0.310096, 0.314479, 0.318893, 0.323337, 0.327812,
      0.332318, 0.336855, 0.341423, 0.346021, 0.350650, 0.355309, 0.360000,
      0.364721, 0.369473, 0.374256, 0.379070, 0.383914, 0.388789, 0.393695,
      0.398631, 0.403599, 0.408597, 0.413626, 0.418685, 0.423775, 0.428897,
      0.434048, 0.439231, 0.444444, 0.449689, 0.454963, 0.460269, 0.465606,
      0.470973, 0.476371, 0.481799, 0.487259, 0.492749, 0.498270, 0.503822,
      0.509404, 0.515017, 0.520661, 0.526336, 0.532042, 0.537778, 0.543545,
      0.549343, 0.555171, 0.561030, 0.566920, 0.572841, 0.578793, 0.584775,
      0.590788, 0.596832, 0.602907, 0.609012, 0.615148, 0.621315, 0.627512,
      0.633741, 0.640000, 0.646290, 0.652611, 0.658962, 0.665344, 0.671757,
      0.678201, 0.684675, 0.691180, 0.697716, 0.704283, 0.710880, 0.717509,
      0.724168, 0.730857, 0.737578, 0.744329, 0.751111, 0.757924, 0.764767,
      0.771642, 0.778547, 0.785483, 0.792449, 0.799446, 0.806474, 0.813533,
      0.820623, 0.827743, 0.834894, 0.842076, 0.849289, 0.856532, 0.863806,
      0.871111, 0.878447, 0.885813, 0.893210, 0.900638, 0.908097, 0.915586,
      0.923106, 0.930657, 0.938239, 0.945852, 0.953495, 0.961169, 0.968874,
      0.976609, 0.984375, 0.992172, -999.000};

#ifdef UNDERSCORE
  #define inisatarps inisatarps_
  #define ctim2abss ctim2abss_
  #define solr1r2 solr1r2_
  #define solrsc1 solrsc1_
  #define solrsc2 solrsc2_
  #define solcorset solcorset_
  #define wtsatfld wtsatfld_
  #define lltoxy lltoxy_
  #define xytoll xytoll_
  #define setmapr setmapr_
  #define setorig setorig_
#endif
void inisatarps();
void ctim2abss();
void solr1r2();
void solrsc1();
void solrsc2();
void solcorset();
void wtsatfld();
void lltoxy();
void xytoll();
void setmapr();
void setorig();

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

  int fd_sat;
  uint8 *sat_image;
  float *sat_lat;
  float *sat_lon;

/* Remapping variables */

  int one=1;
  int two=2;
  float xp,yp;

/* ARPS grid variables */

  long  *grid_knt;
  long  *grid_accum;

  int nx_arps;
  int ny_arps;
  long grid_size;
  long nfield;
  float dx_arps;
  float dy_arps;
  float xnw,ynw;

/* ARPS output variables */

  float *arps_tem;
  float *arps_grd;
  float *arps_out;
  char arunnam[81] = {81*'\0'};
  char varname[7] = {7*'\0'};
  char outfname[101] = {101*'\0'};
  char satname[7];

  uint8 colors[256*3];
  uint8 *grid_val;
  uint16 status;
  char hdffname[100] = {100*'\0'};

  char rawfname[100] = {100*'\0'};
  int fd_raw;

/* Misc local variables */

  int ireturn;
  long i,ii,j,jj,k,kk,koff,kmci,irad,iser,jser;
  long grid_i,grid_j,grid_k;
  long sum,knt,brite;
  uint8 abyte,bbyte;
  uint16 a2byte,b2byte;
  float fi,fj;
  /* Don't worry about source or time info in adas sat data file */
  long isource=0;
  long year=0,month=0,day=0,hour=0,min=0,sec=0;

  long n_read;
  float f;
  int itbl;
  float xpix,ypix;

  int c;
  extern char *optarg;
  extern int optind;
  int errflg=0;

  int dmpfmt;
  int hdf4cmpr;

/************************************************************

  Begin executable code

************************************************************/
  dmpfmt=1;
  hdf4cmpr=0;

/* Read in arguments */

    while ((c = getopt(argc, argv, "n:r:p:c:t:fbh:")) != EOF) {
      switch (c) {
      case 'n':  /* sat_nx, sat_ny */
	sscanf(optarg,"%d",&sat_nx);
	sscanf(argv[optind],"%d",&sat_ny); optind = optind + 1;
        printf("sat_nx %d sat_ny %d\n",sat_nx,sat_ny);
        break;
      case 'r':  /* sat_res_m */
	sscanf(optarg,"%f",&sat_res_m);
        printf("sat_res_m %f\n",sat_res_m);
        break;
      case 'p':  /* sat_mapproj, sat_trulat[2], sat_trulon */
	sscanf(optarg,"%d",&sat_mapproj);
	sscanf(argv[optind],"%f",sat_trulat); optind = optind + 1;
	sscanf(argv[optind],"%f",sat_trulat+1); optind = optind + 1;
	sscanf(argv[optind],"%f",&sat_trulon); optind = optind + 1;
        printf("sat_mapproj %d sat_trulat %f %f sat_trulon\n",
           sat_mapproj,sat_trulat[0],sat_trulat[1],sat_trulon);
        break;
      case 'c':  /* sat_ctrlat, sat_ctrlon */
	sscanf(optarg,"%f",&sat_ctrlat);
	sscanf(argv[optind],"%f",&sat_ctrlon); optind = optind + 1;
        printf("sat_ctrlat %f sat_ctrlon %f\n",sat_ctrlat,sat_ctrlon);
        break;
      case 't':  /* sat_indx */
	sscanf(optarg,"%d",&sat_indx);
        if (sat_indx >= NUM_TBLS) { errflg++; }
        printf("sat_indx %d\n",sat_indx);
        break;
      case 'f':  /* fill_flg */
	fill_flg = 1;
        printf ("fill_flg set to true\n");
        break;
      case 'h':
	dmpfmt = 2;
	sscanf(optarg,"%d",&hdf4cmpr);
	printf("Write HDF4 format data with compression option %d.\n",hdf4cmpr);
	break;
      case 'b':
	dmpfmt = 1;
	printf("Write binary format data.\n");
	break;
      case '?':
        errflg++;
        break;
      }
    }

/* Check argument count and existance of satellite file */

  if( errflg || argc-optind != 2 ) {
      fprintf(stderr, " Usage: %s [options] ", argv[0]);
      fprintf(stderr, "<raw_image_file> <output_header> < <sat.input>\n");
      fprintf(stderr,"Options:\n  -n <sat_points_wide> <sat_points_high>\n");
         /* sat_nx, sat_ny */
      fprintf(stderr,"  -r <sat_res_m>\n");
         /* sat_res_m */
      fprintf(stderr,"  -p <sat_map_proj> <trulat1> <trulat2> <trulon>\n");
         /* sat_mapproj, sat_trulat[2], sat_trulon */
      fprintf(stderr,"  -c <sat_ctrlat> <sat_ctrlon>\n");
         /* sat_ctrlat, sat_ctrlon */
      fprintf(stderr,"  -f                 [hole filling]\n");
         /* fill_flg */
      fprintf(stderr,"  -t <sat_index>\n");
      fprintf(stderr,
         "      0-NESDIS IR, 1-Unisys IR, 2-NESDIS Vis\n");
      fprintf(stderr,"  -b                 [binary format]\n");
      fprintf(stderr,"  -h <hdf4compr>\n");
      fprintf(stderr,
         "      0-- 6, see hdmpfmt in arps.input\n");
         /* sat_table */
      exit(1);

   }

  if( (fd_sat = open(argv[optind],0)) == -1 ) {
      fprintf(stderr,"Unable to open satellite image %s\n",argv[optind] );
      exit(1);
   }

  sat_image = (uint8 *)malloc(sizeof(uint8)*sat_nx*sat_ny);
  sat_lat = (float *)malloc(sizeof(float)*sat_nx*sat_ny);
  sat_lon = (float *)malloc(sizeof(float)*sat_nx*sat_ny);

/* Read-in satellite data */

  printf(" Reading in satellite data\n");

  n_read = read(fd_sat,sat_image,sat_nx*sat_ny);
  if (n_read != sat_nx*sat_ny) {
     fprintf(stderr,
        "Error reading satellite image data.  Read %d of %d bytes\n",
        n_read,sat_nx*sat_ny);
     exit(1);
  }

/* Remap data from image projection to 
 * ARPS grid.  Note that origin is set to NW corner point. */

  setmapr(&sat_mapproj,&sat_scale,sat_trulat,&sat_trulon);
  lltoxy(&one,&one,&sat_ctrlat,&sat_ctrlon,&sat_xnw,&sat_ynw);
  sat_xnw = sat_xnw - 0.5 * sat_nx * sat_res_m;
  sat_ynw = sat_ynw + 0.5 * sat_ny * sat_res_m;

  for (j=0; j<sat_ny; j++) {
     for (i=0; i<sat_nx; i++) {
        xpix = sat_xnw + sat_res_m * (i + 0.5);
        /*ypix = sat_ynw - sat_res_m * (sat_ny - j - 1);*/
        ypix = sat_ynw - sat_res_m * (j + 0.5);
        ii = i + sat_nx*j;
        xytoll(&one,&one,&xpix,&ypix,sat_lat+ii,sat_lon+ii);
/*if (i==sat_nx/2) {
printf("XXX j %d lat %f lon %f xnw %f ynw %f\n",j,sat_latnw[ii],sat_lonnw[ii],xpix,ypix);
}*/
     }
  }

/* Read-in ARPS-grid input parameters */

  inisatarps(arunnam,&nx_arps,&ny_arps,
              &dx_arps,&dy_arps,&xnw,&ynw);

  printf(" Back from inisatarps runname:%s\n",arunnam);
  printf("   nx:%i, ny:%i, dx:%f, dy:%f\n",
           nx_arps,ny_arps,dx_arps,dy_arps);
  printf("   xnw:%f, ynw:%f\n",xnw,ynw);

  printf(" Remapping data\n");

/* Allocate memory for ARPS output */

  grid_size=nx_arps*ny_arps;
  arps_tem = (float *)malloc(sizeof(float)*grid_size );
  arps_grd = (float *)malloc(sizeof(float)*grid_size ); 

/* Allocate memory for grid computations and image */

  grid_knt   = (long *)malloc(sizeof(long)*grid_size );
  grid_accum = (long *)malloc(sizeof(long)*grid_size );
  grid_val = (uint8 *)malloc(sizeof(uint8)*grid_size );

/*for (j=0; j<sat_ny; j++) {
ii = sat_nx/2 + sat_nx*j;
printf ("XXX sat_image %d\n",sat_image[ii]);
} */
  if ( grid_val != NULL ) {

     for ( k = 0; k < grid_size ; k++ ){
        *(grid_knt+k)   = 0;
        *(grid_accum+k) = 0;
     }

     for (j=0; j<sat_ny; j++) {
        for (i=0; i<sat_nx; i++) {
           ii = i + sat_nx*j;
           if ( sat_image[ii] != 0 ) {
              lltoxy(&one,&one,sat_lat+ii,sat_lon+ii,&xp,&yp);
              if (xp >= xnw) {
                /* avoid problem with int of -0.999 to 0.9999 
                 * all mapping to 0: */
                grid_i=(long)(((xp - xnw)/dx_arps));
              } else {
                grid_i=-1;
              }
              if (yp <= ynw) {
                grid_j=(long)(((ynw - yp)/dy_arps));
              } else {
                grid_j=-1;
              }
/*if (grid_i+1 != i || sat_ny-grid_j != j) {printf ("XXX mapping error: grid_i+1 %d i %d sat_ny-grid_j %d j %d i %d j %d grid_i %d grid_j %d dx %f dy %f\n",grid_i+1,i,sat_ny-grid_j,j,i,j,grid_i,grid_j,xp-xnw-i*dx_arps,yp-ynw-j*dy_arps);}*/
/*if (i==sat_nx/2) {
printf ( "XXX i %d j %d xpnw %f ypnw %f gi %d gj %d\n",i,j,xpnw,ypnw,grid_i,grid_j);
}*/
              if ( ( grid_i > -1 ) && ( grid_i < nx_arps ) &&
                   ( grid_j > -1 ) && ( grid_j < ny_arps) ) {

                 grid_k = grid_i + (grid_j*nx_arps);
                 grid_knt[grid_k] = grid_knt[grid_k] + 1;
                 grid_accum[grid_k] = grid_accum[grid_k]+sat_image[ii];
              }
           }
        }
     }

/* Find average in each cell */

     printf(" Average calc\n");
     for ( k = 0; k < grid_size ; k++ ){
/*
           if ( k%500 == 0 )
    printf(" k: %i accum: %i  knt: %i\n",k,*(grid_accum+k),*(grid_knt+k));
*/
        if ( *(grid_accum+k) > 0 ) {

          *(grid_val+k) = (uint8) ( 0.5 + 
             ( ((float)*(grid_accum+k)) /
               ((float)*(grid_knt+k)) ) );
          *(arps_tem+k) = ((float)*(grid_accum+k)) /
                                           ((float)*(grid_knt+k));
        } else {
           *(grid_val+k) = 0;
           *(arps_tem+k) = 0.;
        }
     }

/*for (j=0; j<ny_arps; j++) {
ii = nx_arps/2+nx_arps*j;
printf ("XXX1 grid_accum %d grid_knt %d arps_tem %f\n",grid_accum[ii],grid_knt[ii],arps_tem[ii]);
}*/
/* Fill any missing using average of neighbors */

     if (fill_flg == 1) {
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
                               if( *(grid_knt+kk) > 0) {
                                  sum = sum + *(grid_accum+kk);
                                  knt = knt + *(grid_knt+kk);
                               }
                            }
                         }
                      }
                   }
                   if ( knt > 0 ) {
                      *(grid_val+k) = (uint8)(0.5+
                                      (((float)sum)/((float)knt)));
                      *(arps_tem+k) = ((float)sum)/((float)knt);
                      break;
                   }
                }
             }
          }
       }
     }
  
/*for (j=0; j<ny_arps; j++) {
ii = nx_arps/2+nx_arps*j;
printf ("XXX grid_accum %d grid_knt %d arps_tem %f\n",grid_accum[ii],grid_knt[ii],arps_tem[ii]);
}*/
     free(grid_accum);
     free(grid_knt);


     for (ii=0; ii<nx_arps*ny_arps; ii++) {

        if (arps_tem[ii] >= 1. && arps_tem[ii] <= mxindx[sat_indx]) {
           itbl = (int)arps_tem[ii];
           f = (arps_tem[ii] - itbl);
           arps_grd[ii] = (1.-f)*sat_table[sat_indx][itbl] 
                          + f*sat_table[sat_indx][itbl+1];
        } else {
           arps_grd[ii] = sat_table[sat_indx][0];
        }

     }

/* Write out results */

     arps_out = (float *)malloc(sizeof(float)*grid_size ); 
     ireturn = yflip(nx_arps,ny_arps,arps_grd,arps_out);

     if (map_type[sat_indx] == 0) {   /* map_type 0-IR, 1-Vis */
       strcpy(varname, "satctt\0");
     } else {
       strcpy(varname, "satalb\0");
     }

     sprintf( outfname, "%s.%s", argv[optind+1],varname );
     printf( " Output file name: %s\n",outfname);
     sprintf( hdffname, "%s.%s.hdf", argv[optind+1],varname );
     printf( " hdf image file name: %s\n",hdffname);
     sprintf( rawfname, "%s.%s.raw", argv[optind+1],varname );
     printf( " raw image file name: %s\n",rawfname);

/*   Call FORTRAN routine to write ARPS gridded data */

     printf(" Calling FORTRAN write routine\n");
     nfield=1;
     wtsatfld(&nx_arps,&ny_arps,&nfield,
               outfname,satname,&sat_ctrlat,&sat_ctrlon,
               &year,&month,&day,&hour,&min,&sec,&isource,
	       &dmpfmt,&hdf4cmpr,varname,arps_out);
     printf(" Back from FORTRAN write routine\n");

/*   Rescale image data for a better greyscale depiction

     printf(" Rescaling image\n");
     status = rescale( grid_val, grid_size);
*/

/*   Reverse video IR pictures so warm is black */
/*FIXME
     if ( irchan != 0 ) {
       printf(" Applying reverse video\n");
       status = revvideo( grid_val, grid_size);
     }
*/

/*   Write out raw image file */

     fd_raw = creat(rawfname,0664);
     write(fd_raw,grid_val,grid_size);
     close(fd_raw);

/*   Get the palette data from file. */

      status = DFPgetpal("palgrey.hdf",colors);
      printf("DFPgetpal return status: %d\n",status);

/*   Set the current palette. */

      status = DFR8setpalette(colors);
      printf("DFR8setpalette return status: %d\n",status);
 
/*   Write the image data to the file. */

      status=DFR8putimage(hdffname,grid_val,nx_arps,ny_arps,0);
      printf("DFR8putimage return status: %d\n",status);

   } else {
      fprintf(stderr," Error allocating memory");
      exit(1);
   }
   
  exit(0);

}
