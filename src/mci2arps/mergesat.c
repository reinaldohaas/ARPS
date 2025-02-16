/*  wdt Copyright (c) 2001 Weather Decision Technologies, Inc. 2001-03-23 GMB
*
*   mergesat.c
*
*   Program to merge multiple ARPS satellite format files which cover the
*   same domain.  Non-missing values are averaged.
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
#include <unistd.h>

#include "hdf.h"

#ifdef UNDERSCORE
  #define wtsatfld wtsatfld_
  #define ld2dgrid ld2dgrid_
  #define checkgrid2d checkgrid2d_
#endif
void wtsatfld();
void ld2dgrid();
void checkgrid2d();

int fill_flg = 1; /* no hole filling - 0, fill missing - 1 */

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
  uint8 pix_val;
  int pix_int;

/* Remapping variables */

  int one=1;
  int two=2;
  float xp,yp;

/* ARPS grid variables */

  int nx_arps;
  int ny_arps;
  long grid_size;
  long nfield;
  float dx_arps;
  float dy_arps;
  float xnw,ynw;

/* ARPS output variables */

  float *satfld;
  float *satfld_in;
  float *sattot;
  int *satcnt;
  char arunnam[81] = {81*'\0'};
  char varname[7] = {7*'\0'};
  char outfname[101] = {101*'\0'};
  char satname[7];

  uint8 colors[256*3];
  uint16 status;
  char hdffname[100] = {100*'\0'};

  char rawfname[100] = {100*'\0'};
  int fd_raw;

  float miss_val=-999;

/* Misc local variables */

  int ireturn;
  long i,ii,ii2,j,jj,k,kk,koff,kmci,irad,iser,jser;
  long grid_i,grid_j,grid_k;
  long sum,knt,brite;
  uint8 abyte,bbyte;
  uint16 a2byte,b2byte;
  float fi,fj;
  /* Don't worry about source or time info in adas sat data file */
  long isource=0;
  long year=0,month=0,day=0,hour=0,min=0,sec=0;
  float sat_ctrlat=0,sat_ctrlon=0; /* variables meaningless in merged file */

  int mapproj;
  float dx,dy,ctrlat,ctrlon,trulat1,trulat2,trulon,sclfct;
  int nx_in,ny_in,nfield_in,mapproj_in;
  float dx_in,dy_in,ctrlat_in,ctrlon_in,trulat1_in,trulat2_in,
        trulon_in,sclfct_in;
  int grid_stat;

  long n_read;
  float f;
  int itbl;
  float xpix,ypix;

  int c;
  extern char *optarg;
  extern int optind;
  int errflg=0;
  int firstind;

  int dmpfmt;
  int hdf4cmpr;

/************************************************************

  Begin executable code

************************************************************/

  dmpfmt = 1;
  hdf4cmpr = 0;

  if( errflg || argc <= 2 ) {
      fprintf(stderr, " Usage: %s ", argv[0]);
      fprintf(stderr, 
         "[-b|-h 0-6] <file1> [<file2> ...] <Merged_arps_sat_file>\n");
      exit(1);

   }

  while ((c = getopt(argc, argv, "bh:")) != EOF) {
    switch (c) {
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

  firstind = optind;
  for (; optind < argc-1; optind++) { /* process each data file */

    printf("Reading %s\n",argv[optind]);
    if( (fd_sat = open(argv[optind],0)) == -1 ) {
        fprintf(stderr,"Unable to open arps format satellite data file %s\n",
           argv[optind] );
        exit(1);
     }

     /* read in header information */

  /* From SUBROUTINE rdsatfld in rdsatfld.f90:

  CHARACTER (LEN=6) :: satname
  CHARACTER (LEN=6) :: fldname(nfield)
  REAL :: satfld(nx,ny,nfield)
  CHARACTER (LEN=80  ) :: runname

    READ (iunit,ERR=200) satname
    READ (iunit,ERR=200) nxin,nyin,nfieldin,itime,idummy,               &
                         idummy,idummy,idummy,idummy,idummy
    READ (iunit,ERR=200) dummy_name
    READ (iunit,ERR=200) idummy,idummy,mprojin,idummy,idummy,          &
                 idummy,idummy,idummy,idummy,idummy
    READ (iunit,ERR=200) dxin,dyin,rdummy,rdummy,ctrlatin,  &
                 ctrlonin,trlat1in,trlat2in,trlonin,sclfctin,  &
                 latsat,lonsat,rdummy,rdummy,rdummy
    READ(iunit,ERR=200) fldname
    READ(iunit,ERR=200) satfld
  */

     /* satname - skip */
     lseek(fd_sat,14L,1);

     /* nx,ny,nfield */

     lseek(fd_sat,4L,1);
     read(fd_sat,&nx_in,4); 
     read(fd_sat,&ny_in,4); 
     read(fd_sat,&nfield_in,4); 
     lseek(fd_sat,32L,1);

     if (nfield_in != 1) {
       fprintf(stderr,"ERROR, nfield != 1 (%d)\n",nfield_in);
     }

     if (optind == firstind) {
       nx_arps = nx_in;
       ny_arps = ny_in;
       nfield = 1;
     }

     /* runname - skip */

     lseek(fd_sat,88L,1);

     /*  mapproj */

     lseek(fd_sat,12L,1);
     read(fd_sat,&mapproj_in,4);
     lseek(fd_sat,32L,1);
     if (optind == firstind) {
       mapproj = mapproj_in;
     }

     /* dx,dy,ctrlat,ctrlon,trulat1,trulat2,trulon,sclfct */

     lseek(fd_sat,4L,1);
     read(fd_sat,&dx_in,4);
     read(fd_sat,&dy_in,4);
     lseek(fd_sat,8L,1);
     read(fd_sat,&ctrlat_in,4);
     read(fd_sat,&ctrlon_in,4);
     read(fd_sat,&trulat1_in,4);
     read(fd_sat,&trulat2_in,4);
     read(fd_sat,&trulon_in,4);
     read(fd_sat,&sclfct_in,4);
     lseek(fd_sat,24L,1);
     if (optind == firstind) {
       dx = dx_in;
       dy = dy_in;
       ctrlat = ctrlat_in;
       ctrlon = ctrlon_in;
       trulat1 = trulat1_in;
       trulat2 = trulat2_in;
       trulon = trulon_in;
       sclfct = sclfct_in;
     }

     /* check grid information */
     
     if (optind == firstind) {
       ld2dgrid(&dx,&dy,&ctrlat,&ctrlon,
                &mapproj,&trulat1,&trulat2,&trulon,&sclfct);
     } else {
       checkgrid2d(&nx_arps,&ny_arps,&nx_in,&ny_in,
                   &dx,&dy,&ctrlat,&ctrlon,
                   &mapproj,&trulat1,&trulat2,&trulon,&sclfct,
                   &dx_in,&dy_in,&ctrlat_in,&ctrlon_in,
                   &mapproj_in,&trulat1_in,&trulat2_in,&trulon_in,&sclfct_in,
                   &grid_stat);
       if (grid_stat != 0) {
         fprintf(stderr,"ERROR: grids don't match between files (%s & %s)\n",
                 argv[firstind],argv[optind]);
       }
     }

     /* fldname */

     lseek(fd_sat,4L,1);
     if (optind == firstind) {
       n_read = read(fd_sat,varname,6);
       if (n_read != 6) {
         fprintf(stderr,"ERROR reading varname (only %d bytes read)\n",n_read);
         exit(1);
       }
     } else {
       lseek(fd_sat,6L,1);
     }
     lseek(fd_sat,4L,1);

     /* satfld */

     if (optind == firstind) { /* Allocate memory for ARPS output */
       grid_size=nx_arps*ny_arps;
       satfld_in = (float *)malloc(sizeof(float)*grid_size );
       satfld = (float *)malloc(sizeof(float)*grid_size ); 
       sattot = (float *)malloc(sizeof(float)*grid_size ); 
       satcnt = (int *)malloc(sizeof(int)*grid_size ); 
       for (ii=0; ii<grid_size; ii++) {
         satfld[ii] = miss_val;
         sattot[ii] = 0.;
         satcnt[ii] = 0;
       }
     }

     lseek(fd_sat,4L,1);

     n_read = read(fd_sat,satfld_in,nx_arps*ny_arps*nfield*4);
     if (n_read != nx_arps*ny_arps*4) {
       fprintf(stderr,
          "ERROR reading satellite data (read only %d of %d bytes)\n",
          n_read,nx_arps*ny_arps*4);
       exit(1);
     }
     lseek(fd_sat,4L,1);

     for (ii=0; ii<grid_size; ii++) {
       if (satfld_in[ii] != miss_val) {
/*if (satcnt[ii] == 1) {
 *printf("XXX %5.1f %5.1f\n",sattot[ii],satfld_in[ii]); } */
         sattot[ii] = sattot[ii] + satfld_in[ii];
         satcnt[ii] = satcnt[ii] + 1;
       }
     }

  } /* end loop for processing each data file */

  for (ii=0; ii<grid_size; ii++) {
    if (satcnt[ii] > 0) {
      satfld[ii] = sattot[ii]/satcnt[ii];
    }
  }

  strcpy(satname, "merged\0");

/* Fill any missing using average of neighbors */

  if (fill_flg == 1) {
    printf(" Hole Filling\n");
    for ( j = 0; j < ny_arps ; j++ ) {
       for ( i = 0; i < nx_arps ; i++ ) {
          k = i + (j*nx_arps);
          if ( *(satcnt+k) == 0 ) {
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
                            if( *(satcnt+kk) > 0) {
                               sum = sum + *(sattot+kk);
                               knt = knt + *(satcnt+kk);
                            }
                         }
                      }
                   }
                }
                if ( knt > 0 ) {
                   *(satfld+k) = ((float)sum)/((float)knt);
                   break;
                }
             }
          }
       }
    }
  }
 
/* Write out results */


  sat_image = (uint8 *)malloc(sizeof(uint8)*grid_size );

  if (strncmp(varname,"satctt",6) == 0) {
    for (ii=0; ii<grid_size; ii++) {
      if (satfld[ii] == miss_val) {
        sat_image[ii] = 0;
      } else {
        if (satfld[ii] >= 240.) {
          pix_int = (420. - satfld[ii] + 0.5);
          if (pix_int < 255) {
            sat_image[ii] = pix_int;
          } else {
            sat_image[ii] = 1;
          }
        } else {
          pix_int = (660. - 2.*satfld[ii] + 0.5);
          if (pix_int > 1) {
            sat_image[ii] = pix_int;
          } else {
            sat_image[ii] = 255;
          }
        }
      }
    }
  } else if (strncmp(varname,"satalb",6) == 0) {
    for (ii=0; ii<grid_size; ii++) {
      if (satfld[ii] == miss_val) {
        sat_image[ii] = 0;
      } else {
        pix_int = (satfld[ii] * 254. + 0.5);
        if (pix_int > 254) {
          pix_int = 1;
        } else if (pix_int < 2) {
          pix_int = 255;
        }
        sat_image[ii] = pix_int;
      }
    }
  } else {
    fprintf(stderr, 
       "WARNING: unknown varname %s, blank 8-bit image made.\n", varname);
    for (ii=0; ii<grid_size; ii++) {
      sat_image[ii] = 0;
    }
  }

  /* flip the image */
  for (j=0; j<ny_arps/2; j++) {
    ii = j*nx_arps;
    ii2 = (ny_arps-j-1)*nx_arps;
    for (i=0; i<nx_arps; i++) {
      pix_val = sat_image[ii+i];
      sat_image[ii+i] = sat_image[ii2+i];
      sat_image[ii2+i] = pix_val;
    }
  }

/*
 *  sprintf( outfname, "%s.%s", argv[argc-1],varname );
 *  printf( " Output file name: %s\n",outfname);
 *  sprintf( hdffname, "%s.%s.hdf", argv[argc-1],varname );
 *  printf( " hdf image file name: %s\n",hdffname);
 *  sprintf( rawfname, "%s.%s.raw", argv[argc-1],varname );
 *  printf( " raw image file name: %s\n",rawfname);
 */
  sprintf( outfname, "%s", argv[argc-1] );
  printf( " Output file name: %s\n",outfname);
  sprintf( hdffname, "%s.hdf", argv[argc-1] );
  printf( " hdf image file name: %s\n",hdffname);
  sprintf( rawfname, "%s.raw", argv[argc-1] );
  printf( " raw image file name: %s\n",rawfname);

/*   Call FORTRAN routine to write ARPS gridded data */

  printf(" Calling FORTRAN write routine\n");
  nfield=1;
  wtsatfld(&nx_arps,&ny_arps,&nfield,
            outfname,satname,&sat_ctrlat,&sat_ctrlon,
            &year,&month,&day,&hour,&min,&sec,&isource,
            &dmpfmt,&hdf4cmpr,varname,satfld);
  printf(" Back from FORTRAN write routine\n");

/*   Write out raw image file */

  fd_raw = creat(rawfname,0664);
  write(fd_raw,sat_image,grid_size);
  close(fd_raw);

/*   Get the palette data from file. */

  status = DFPgetpal("palgrey.hdf",colors);
  printf("DFPgetpal return status: %d\n",status);

/*   Set the current palette. */

  status = DFR8setpalette(colors);
  printf("DFR8setpalette return status: %d\n",status);

/*   Write the image data to the file. */

  status=DFR8putimage(hdffname,sat_image,nx_arps,ny_arps,0);
  printf("DFR8putimage return status: %d\n",status);

  exit(0);

}
