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

   ROUTINE: sat_read.c
   PROGRAMMER: Dan Vietor
   PROGRAM TYPE: WXP_LIBRARY
   VERSION: 1.3
   WXP VERSION: 4.8
   DATE: 930915

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

   DESCRIPTION: This subroutine reads McIDAS area files and passes
      a data structure back with the image and navigation 
      information.

   REVISIONS:                                    DATE:       INIT:
Version 1                                        920115      DEV
Added navigation info                            920515      DEV  1.1
Debugged read in routines                        921015      DEV  1.2
Debugged read in routines,
Added SAT_read interface for generic images      930915      DEV  1.3

***********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include "wxp.h"
#include "mc_area.h"

static FILE *ifile;        /* Input file */
static unsigned char data[4000]; /* Input data from file */
static unsigned char scan_line[1600]; /* Current scan line */
static int dat_ind;        /* Index into data array */
static int dat_buf;        /* Current number of bytes in data array */
static int scan_ind;       /* Index to current scan line */
static int nib_bytes;      /* Current number of 4 bit nibbles read */
static int low = 0;        /* Specifies which 4 bit nibble to read */

/****************************************************************
   GETWORD - Gets a 4 byte word from file
****************************************************************/
unsigned int getword()

{
unsigned int output;

if( dat_ind >= dat_buf ){
   dat_buf = fread( data,1,4000,ifile );
   dat_ind = 0;
   }
output = *(int *)(data+dat_ind);
dat_ind += 4;
return output;
}

/****************************************************************
   GETBYTE - Gets a single byte from file
****************************************************************/
int getbyte()

{
if( dat_ind >= dat_buf ){
   dat_buf = fread( data,1,4000,ifile );
   if( dat_buf == -1 ) return -1;
   dat_ind = 0;
   }
return (int)data[dat_ind++];
}

/****************************************************************
   GETNIB - Gets a 4 bit nibble from file
****************************************************************/
int getnib()

{
if( dat_ind >= dat_buf ){
   dat_buf = fread( data,1,4000,ifile );
   if( dat_buf == -1 ) return -1;
   dat_ind = 0;
   }
if( !low ){
   low = 1;
   return data[dat_ind] >> 4;
   }
else {
   low = 0;
   nib_bytes++;
   return data[dat_ind++] & 0xF;
   }
}

/****************************************************************
   GETSYNC - Reads in the synch code from file
****************************************************************/
getsync()

{
int number;
int err;

number = 0;
while( number < 1024 ){
  number++;
  if(( err=getbyte()) == -1 ) return -1;
  if( err != 0x00 ) continue;
  if( getbyte() != 0xFF ) continue;
  if( getbyte() != 0x00 ) continue;
  if( getbyte() != 0xFF ) continue;
  if( getbyte() != 0x00 ) continue;
  if( getbyte() != 0xFF ) continue;
  if( getbyte() != 0x00 ) continue;
  if( getbyte() != 0xFF ) continue;
  return 1;
  }
return 0;
}

/****************************************************************
   SWAPWORD - swaps the bytes in a word
****************************************************************/
void swapword( ptr )

unsigned char *ptr;

{
unsigned char ch;

ch = *ptr;
*ptr = *(ptr+3);
*(ptr+3) = ch;
ch = *(ptr+1);
*(ptr+1) = *(ptr+2);
*(ptr+2) = ch;
}

/****************************************************************
   AREA_READ - Read in the area from file, determines whether
   packed or unpacked and performs unpacking.
****************************************************************/
int AREA_read( file, area )

FILE *file;
struct mc_area **area;    /* Area data structure */

{
int i,j,k;                 /* Loop control variables */
int mode;                /* Mode of image in bits */
int size;                /* Size of packed line in nibbles */
int pre_scan;            /* Previous scan line */
int scan;                /* Current scan line in image */
int width;               /* Width of image */
int memwidth;            /* Width of image memory*/
int height;              /* Height of image */
int depth;               /* data element size */
int file_type;           /* File type */
int cal;                 /* calibration block location */
int datloc;              /* data block location */
int status;              /* Function status */
unsigned int temp;       /* Temporary integer variable */
unsigned int ch;         /* Bytes read in from file */
unsigned int ref;        /* Previous pixel value */
unsigned char *ptr;      /* Pointer into structures */
unsigned char *glddat;   /* Pointer to bytes in Gould format */
unsigned long swaptest = 1; /* Used for determination of byte order */

double GouldToNative (unsigned char *);

swaptest = *(char *)&swaptest;
/*
   Set the file pointer
*/
ifile = file;
if( ifile == NULL )
   return -1;
rewind( file );
/*
   Allocate the area data structures
*/
*area = (struct mc_area *)malloc(sizeof(struct mc_area));

dat_ind = dat_buf = 0;
/*
   Determine file type
   if 4 or 0x4000000 => unpacked
   else              => packed
*/
temp = getword();
temp = getword();
if( temp == 4 )
   file_type = 1;
else if( temp == 0x4000000 )
   file_type = 2;
else
   file_type = 0;
dat_ind = 0;
/*
   If packed, search for sync bytes
*/
if( !file_type && getsync() != 1 ) return -1;
/*
  Read in area directory
*/
(*area)->dir = (struct area_dir *)malloc(sizeof(struct area_dir));
memcpy( (*area)->dir,data+dat_ind,sizeof(struct area_dir));
if(( file_type == 0 && swaptest ) || file_type == 2 ){
   ptr = (unsigned char *)(*area)->dir;
   for( i = 0; i < 24; i++,ptr+=4 )
      swapword( ptr );
   ptr = (unsigned char *)(*area)->dir+32*4;
   for( i = 32; i < 64; i++,ptr+=4 )
      swapword( ptr );
   }
if( !file_type && getsync() != 1 ) return -1;
/*
   Read in area navigation
*/
dat_ind = (*area)->dir->navoffst;
if( !strncmp( (const char*)data+dat_ind,"GOES",4 )){
   printf(" Nav data is GOES\n");
   strcpy( (*area)->type,"GOES" );
   (*area)->navg = (struct nav_goes *)malloc(sizeof(struct nav_goes));
   memcpy( (*area)->navg,data+dat_ind,sizeof(struct nav_goes));
   if(( file_type == 0 && swaptest ) || file_type == 2 ){
      ptr = (unsigned char *)(*area)->navg + 4;
      for( i = 1; i < 128; i++,ptr+=4 )
         swapword( ptr );
      }
   dat_ind+=512;
   }
else if( !strncmp( (const char *)data+dat_ind,"GVAR",4 )){
   printf(" Nav data is GVAR\n");
   strcpy( (*area)->type,"GVAR" );
   (*area)->navgv = (struct nav_gvar *)malloc(sizeof(struct nav_gvar));
   memcpy( (*area)->navgv,data+dat_ind,sizeof(struct nav_gvar));
   if(( file_type == 0 && swaptest ) || file_type == 2 ){
      ptr = (unsigned char *)(*area)->navgv + 4;
      for( i = 1; i < 640; i++,ptr+=4 )
         swapword( ptr );
      }
   dat_ind+=2560;
   }

/* Read calibration data, if any */

if( (*area)->dir->caloffst != 0 ) {

/*
   Allocate calibration data
*/
    (*area)->calgv =
          (struct cal_gvar *)malloc(sizeof(struct cal_gvar));
/* 
   Position file pointer to the beginning of the calibration block
*/
    cal = (*area)->dir->caloffst;
    status = fseek( ifile, cal, 0);
/*
   Allocate dummy byte data array for reading Gould-formatted reals
   as a sequence of bytes.
*/

    glddat = (unsigned char *)(malloc(sizeof(struct cal_gvar)));

/*
   Read calibration data
*/

    printf("Reading and converting calibration data");
    fread( glddat, 1, sizeof(struct cal_gvar), ifile );

/*
   Translate Gould floating point numbers to IEEE
   A Gould floating-point number is 4-bytes long, hence the steps
   by four.

*/
    for ( i = 0; i < 8 ; i++ ) {
        k=i*4;
        ((*area)->calgv->visbias[i])=(float)GouldToNative(glddat+k);
        printf("visbias[%d]=%10.4f\n",i,(*area)->calgv->visbias[i]);
        }

    for ( i = 0; i < 8 ; i++ ) {
        k=32+i*4;
        (*area)->calgv->vis1gain[i]=(float)GouldToNative(glddat+k);
        printf("vis1gain[%d]=%10.4f\n",i,(*area)->calgv->vis1gain[i]);
        }

    for ( i = 0; i < 8 ; i++ ) {
        k=64+i*4;
        (*area)->calgv->vis2gain[i]=(float)GouldToNative(glddat+k);
        printf("vis2gain[%d]=%10.4f\n",i,(*area)->calgv->vis2gain[i]);
         }

    (*area)->calgv->albedcon=(float)GouldToNative(glddat+96);
    printf("albedo constant=%10.4f\n",(*area)->calgv->albedcon);

    for ( i = 0; i < 4 ; i++ ) {
        k=100+i*4;
        (*area)->calgv->ir1bias[i]=(float)GouldToNative(glddat+k);
        printf("ir1bias[%d]=%10.4f\n",i,(*area)->calgv->ir1bias[i]);
        }

    for ( i = 0; i < 4 ; i++ ) {
        k=116+i*4;
        (*area)->calgv->ir2bias[i]=(float)GouldToNative(glddat+k);
        printf("ir2bias[%d]=%10.4f\n",i,(*area)->calgv->ir2bias[i]);
        }

    for ( i = 0; i < 4 ; i++ ) {
        k=132+i*4;
        (*area)->calgv->ir1gain[i]=(float)GouldToNative(glddat+k);
        printf("ir1gain[%d]=%10.4f\n",i,(*area)->calgv->ir1gain[i]);
        }

    for ( i = 0; i < 4 ; i++ ) {
        k=148+i*4;
        (*area)->calgv->ir2gain[i]=(float)GouldToNative(glddat+k);
        printf("ir2gain[%d]=%10.4f\n",i,(*area)->calgv->ir2gain[i]);
        }

    for ( i = 0; i < 87; i++ ) {
        (*area)->calgv->calresvd[i]=0.;
        }

/*
    Free-up the temporary memory used for reading.
*/

    free (glddat);
    
    }

width = (*area)->dir->esiz;
height = (*area)->dir->lsiz;
depth = (*area)->dir->zsiz;
memwidth = (*area)->dir->yzprefix+(width*depth);

if( width <= 0 || width > 100000 ) return -1;
if( height <= 0 || height > 100000 ) return -1;
/*
   Allocate image data
*/

ptr = (*area)->image =
   (unsigned char *)malloc( memwidth * height );

/* 
   Position file pointer to the beginning of the data block
*/

datloc = (*area)->dir->datoffst;
status = fseek( ifile, datloc, 0);

/*
   If unpacked, just read in
*/

if( file_type != 0 ){
   fread( ptr, 1, (memwidth*height), ifile );
   }
/*
   Unpack data
*/
else {
   scan = 0;
   for( j = 0; j < height; j++ ){
      for( i = 0; i < width; i++ )
         scan_line[i] = 0;
      pre_scan = scan;
/*
   Get the number of the next scan line in the image
*/
      if( getsync() == 1 ){
         scan = (unsigned long)getbyte() << 24;
         scan |= (unsigned long)getbyte() << 16;
         scan |= (unsigned long)getbyte() << 8;
         scan |= (unsigned long)getbyte();
         if(( pre_scan > 0 && scan > pre_scan + 10 ) ||
              scan > height || scan < pre_scan ) scan = pre_scan;
         }
      else 
         scan = height;
/*
   Fill image array for missing scan lines
*/
      for( i = pre_scan; i < scan-1 && j < height; i++ ){
         memcpy( ptr, scan_line, width );
         ptr += width;
         j++;
         }
     if( scan == height ) continue;
/*
   Find out how many nibbles in line and mode of data
*/
      mode = (unsigned long)getbyte() << 8;
      mode |= getbyte();
      size = (unsigned long)getbyte() << 8;
      size |= getbyte();
      low = 0;
      scan_ind = 0;
      ref = 0;
/*
   Read through each nibble
*/
      for( nib_bytes = 0; nib_bytes < size; nib_bytes++) {
         ch = getnib();
         if( ch == 0xFF || scan_ind >= width ) break;
         switch( ch ){
/*
   New value case, read in next 8 bits
*/
            case 0x0:
               ref = getnib() << 4;
               ref |= getnib();
               scan_line[scan_ind++] = ref; break;
/*
   Other cases, differences from previous pixel value
*/
            case 0x1: ref = scan_line[scan_ind++] = ref - 6; break;
            case 0x2: ref = scan_line[scan_ind++] = ref + 6; break;
            case 0x3: ref = scan_line[scan_ind++] = ref - 5; break;
            case 0x4: ref = scan_line[scan_ind++] = ref + 5; break;
            case 0x5: ref = scan_line[scan_ind++] = ref - 4; break;
            case 0x6: ref = scan_line[scan_ind++] = ref + 4; break;
            case 0x7: ref = scan_line[scan_ind++] = ref - 3; break;
            case 0x8: ref = scan_line[scan_ind++] = ref + 3; break;
            case 0x9: ref = scan_line[scan_ind++] = ref - 2; break;
            case 0xA: ref = scan_line[scan_ind++] = ref + 2; break;
            case 0xB: ref = scan_line[scan_ind++] = ref - 1; break;
            case 0xC: ref = scan_line[scan_ind++] = ref + 1; break;
            case 0xD: scan_line[scan_ind++] = ref; break;
            case 0xE:
               scan_line[scan_ind++] = ref;
               scan_line[scan_ind++] = ref; break;
            case 0xF:
               scan_line[scan_ind++] = ref;
               scan_line[scan_ind++] = ref;
               scan_line[scan_ind++] = ref; break;
            }
         }
/*
   If 6 bit mode, shift to 8 bit mode
*/
      if( mode == 6 )
         for( i = 0; i < width; i++ )
            scan_line[i] = scan_line[i] << 2;
      memcpy( ptr, scan_line, width );
      ptr += width;
      }
   }
return 0;
}
