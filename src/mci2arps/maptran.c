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

   ROUTINE: maptran.c
   PROGRAMMER: Dan Vietor
   PROGRAM TYPE: WXP_LIBRARY
   VERSION: 1.3
   WXP VERSION: 4.8
   DATE: 930925

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

   DESCRIPTION: This subroutine calculates map projection variables.

   REVISIONS:                                    DATE:       INIT:
See map.c for other revision information
Version 1                                        910701      DEV
Debugged satellite transformation                920520      DEV 1.1
Added new projections                            930815      DEV 1.2
Added new projection parameters                  930925      DEV 1.3

****************************************************************************/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "wxp.h"

static int map_projection = PROJ_PSTEREO; /* The map projection */
static float proj_lat;         /* The projection base latitude */
static float proj_lon;         /* The projection base longitude */
static float sat_rad = 421.61;
static float sat_rad2 = 177754.9921;
static float sat_rad4 = 3.159683722e10;
static float sat2_earth2 = 173696.028;
static float earth_rad = 63.71;
static float sat_ratio = MISS;
static float earth_ratio = MISS;
static float tlat1 = 30;
static float tlat2 = 60;
static float cfactor = MISS;
static float scale;

/*******************************************************************
   MAPPROJ - This function sets the type of projection used in plots
********************************************************************/
int mapproj_str( proj )

char *proj;

{
int type;
if( !strncmp( proj, "ll", 2 ))
   type = PROJ_LATLON;
else if( !strncmp( proj, "ps", 2 ))
   type = PROJ_PSTEREO;
else if( !strncmp( proj, "sat", 3 ))
   type = PROJ_SAT;
else if( !strncmp( proj, "lc", 2 ))
   type = PROJ_LAMB;
else if( !strncmp( proj, "me", 2 ))
   type = PROJ_MERC;
else if( !strncmp( proj, "pix", 3 ))
   type = PROJ_PIXEL;
else if( !strncmp( proj, "xy", 2 ))
   type = PROJ_XY;
else if( !strcmp( proj, "xlogp" ))
   type = PROJ_XLOGP;
else if( !strcmp( proj, "cat" ))
   type = PROJ_CAT;
else if( !strcmp( proj, "therm" ))
   type = PROJ_STUVE;
else if( !strcmp( proj, "hodo" ))
   type = PROJ_POLAR;
else if( !strcmp( proj, "vert" ))
   type = PROJ_XLOGP;
else if( !strcmp( proj, "miss" ))
   type = MISS;
else
   sscanf( proj,"%d",&type );
return type;
}

/*******************************************************************
   MAPPROJ - This function sets the type of projection used in plots
********************************************************************/
void mapproj( proj )

int proj;

{
map_projection = proj;
cfactor = MISS;
sat_ratio = MISS;
}

/*******************************************************************
   MAPGETPROJ - This function gets the type of projection
********************************************************************/
int mapgetproj( )

{
return map_projection;
}

/*******************************************************************
   MAPPROJCOORD - This function sets extra map projection coordinates.
********************************************************************/
void mapprojcoord( value )

float value[];

{
if( map_projection == PROJ_LAMB ){
   tlat1 = value[0];
   tlat2 = value[1];
   }
}

/*******************************************************************
   MAPCOORD - This function sets the map coordinates.
********************************************************************/
void mapcoord( rlat,rlon )

float rlat,rlon;

{
proj_lat = rlat;
proj_lon = rlon;
cfactor = MISS;
}

/*******************************************************************
   MAPUNTRAN - This function converts projection coordinates back into 
   latitude and longitude.
********************************************************************/
int mapuntran( x,y,lat,lon )

float x,y;
float *lat,*lon;

{
float rad,ang;

/*
   Polar stereographic projection
*/
if( map_projection == PROJ_PSTEREO ){
   rad = sqrt( x*x + y*y );
   if( proj_lat >= 0 )
      scale = 1;
   else
      scale = -1;
   ang = atan2( x,-scale*y )*RDC;
   *lat = scale*90 - 2*scale*atan(rad/(1.866*earth_rad))*RDC;
   *lon = ang+proj_lon;
   }
/*
   Lat lon projection
*/
else if( map_projection == PROJ_LATLON ){
   *lon = x + proj_lon;
   if( *lon > 180 ) *lon -= 360;
   else if( *lon < -180 ) *lon += 360;
   *lat = y;
   }
/*
   Mercator
*/
else if( map_projection == PROJ_MERC ){
   *lon = (x/(earth_rad*DRC)) + proj_lon;
   if( *lon > 180 ) *lon -= 360;
   else if( *lon < -180 ) *lon += 360;
   *lat = 90 - 2.*RDC*atan(1./exp(y/earth_rad));
   }
/*
   Lambert conformal  
*/
else if( map_projection == PROJ_LAMB ){
   if( cfactor == MISS ){
      if( proj_lat >= 0 )
         scale = 1;
      else
         scale = -1;     
      if( tlat1 == tlat2 )
         cfactor = cos(( 90 - scale*tlat1 )*DRC );
      else
         cfactor = log( cos( tlat1*DRC )) - log( cos( tlat2*DRC ))/
            (log( tan((45-tlat1/2)*DRC)) -
             log( tan((45-tlat2/2)*DRC)));
      }
   rad = sqrt( x*x + y*y );
   ang = atan2( x,-scale*y )*RDC;
   *lat = scale*90 - scale*2*RDC*atan( pow( rad / (2*earth_rad), 1/cfactor ));
   *lon = ang/cfactor+proj_lon;
   }
/*
   Satellite projection
*/
else if( map_projection == PROJ_SAT ){
   if( sat_ratio == MISS ){
      sat_ratio = sat_rad / earth_rad;
      earth_ratio = earth_rad / sat_rad;
      }
   rad = x*x + y*y + sat_rad*sat_rad;
   scale = sat_rad4 - rad*sat2_earth2;
   if( scale < 0 ){
      *lat = MISS;
      *lon = MISS;
      return(0);
      }
   scale = (sat_rad2 - sqrt( scale )) / rad;
   x *= scale;
   y *= scale;
   *lat = asin( y/earth_rad );
   *lon = asin( x/(cos(*lat)*earth_rad))*RDC + proj_lon;
   *lat *= RDC;
   }
   return(1);
}

/*******************************************************************
   MAPTRAN - This function converts latitude and longitude for a 
   particular point into projection coordinates.
********************************************************************/
int maptran( lat,lon,x,y )

float lat,lon;
float *x,*y;

{
float rad,ang;
float factor;

/*
   Polar stereographic projection
*/
if( map_projection == PROJ_PSTEREO ){
   if( proj_lat >= 0 ){
      if( lat < -45 ) lat = -45;
      rad = 1.866*earth_rad*tan( ( 90-lat )*DRC/2 );
      }
   else {
      if( lat > 45 ) lat = 45;
      rad = 1.866*earth_rad*tan( ( 90+lat )*DRC/2 );
      }
   ang = ( lon-proj_lon )*DRC;
   *x = rad*sin(ang);
   *y = rad*cos(ang);
   if( proj_lat >= 0 )
      *y = -*y;
   }
/*
   Lat lon projection
*/
else if( map_projection == PROJ_LATLON ){
#define OFFSET 180
   *x = lon - proj_lon;
   if( *x > OFFSET ) *x -= 2*OFFSET;
   if( *x < -OFFSET ) *x += 2*OFFSET;
   *y = lat;
#undef OFFSET
   }
/*
   Mercator
*/
else if( map_projection == PROJ_MERC ){
#define OFFSET 180*earth_rad*DRC
   *x = earth_rad*( lon-proj_lon ) * DRC;
   if( *x > OFFSET ) *x -= 2*OFFSET;
   else if( *x < -OFFSET ) *x += 2*OFFSET;
   *y = earth_rad*log( 1./tan((45-lat/2)*DRC));
#undef OFFSET
   }
/*
   Lambert conformal  
*/
else if( map_projection == PROJ_LAMB ){
   if( cfactor == MISS ){
      if( proj_lat >= 0 )
         scale = 1;
      else
         scale = -1;     
      if( tlat1 == tlat2 )
         cfactor = cos(( 90 - scale*tlat1 )*DRC );
      else
         cfactor = log( cos( tlat1*DRC )) - log( cos( tlat2*DRC ))/
            (log( tan((45-tlat1/2)*DRC)) -
             log( tan((45-tlat2/2)*DRC)));
      }
   factor = 2*earth_rad*pow(tan((45-scale*lat/2)*DRC),cfactor );
   *x = factor*sin( cfactor*( lon-proj_lon )*DRC );
   *y = -scale*factor*cos( cfactor*( lon-proj_lon )*DRC );
   }
/*
   Satellite projection
*/
else if( map_projection == PROJ_SAT ){
   if( sat_ratio == MISS ){
      sat_ratio = sat_rad / earth_rad;
      earth_ratio = earth_rad / sat_rad;
      }
   factor = cos(lat*DRC)*cos((lon-proj_lon)*DRC);
   if( factor < earth_ratio ){
      *x = MISS; *y = MISS;
      return(0);
      }
   factor = sat_rad - earth_rad*factor;
   scale = sat_rad;

   *x = scale * earth_rad * sin((lon-proj_lon)*DRC) * cos(lat*DRC)/factor;
   *y = scale * earth_rad * sin(lat*DRC)/factor;
   }
/*
   PIXEL projection
*/
else if( map_projection == PROJ_PIXEL ){
   *x = lat;
   *y = -lon;
   }
/*
   XY projection
*/
else if( map_projection == PROJ_XY || map_projection == PROJ_CAT ){
   *x = lat;
   *y = lon;
   }
/*
   XLOGP projection
*/
else if( map_projection == PROJ_XLOGP ){
   *x = lat;
   *y = -log(lon);
   }
return(1);
}

/*******************************************************************
   MAP_OFFSET - This routine offsets a point in a specific direction
   and distance on the earth and returns the point in projection space.
********************************************************************/
void map_offset( xi,yi,dir,rad,xf,yf )

float xi;               /* The x location based on projection */
float yi;               /* The y location based on projection */
float dir;              /* The actual direction based on true north */
float rad;              /* The actual distance in km */
float *xf;              /* The x location based on projection */
float *yf;              /* The y location based on projection */

{
float lat,lon;
float xoff,yoff;

mapuntran( xi,yi,&lat,&lon );
xoff = rad*sin( dir*DRC );
yoff = rad*cos( dir*DRC );
lon += xoff/(111.195*cos(lat*DRC));
lat += yoff/111.195;
maptran( lat,lon,xf,yf );
}

/*******************************************************************
   MAP_DIR_ADJ - This routine adjusts the wind direction based on
   location and map projection.
********************************************************************/
float map_dir_adj( dir,x,y )

float dir;              /* The actual direction based on true north */
float x;                /* The x location based on projection */
float y;                /* The y location based on projection */

{
float lat,lon;
float xoff,yoff;

if( map_projection == PROJ_PSTEREO ){
   if( x == 0 && y == 0 ){
      if( proj_lat >= 0 )
         dir = dir - proj_lon;
      else
         dir = dir + proj_lon;
      }
   else if( proj_lat >= 0 )
      dir = dir - atan2( x,-y ) * RDC;
   else
      dir = dir + atan2( x,y ) * RDC;
   }
else if( map_projection == PROJ_LAMB ){
   if( cfactor == MISS ){
      if( proj_lat >= 0 )
         scale = 1;
      else
         scale = -1;     
      if( tlat1 == tlat2 )
         cfactor = cos(( 90 - scale*tlat1 )*DRC );
      else
         cfactor = log( cos( tlat1*DRC )) - log( cos( tlat2*DRC ))/
            (log( tan((45-scale*tlat1/2)*DRC)) -
             log( tan((45-scale*tlat2/2)*DRC)));
      }
   dir = dir - atan2( x,-scale*y ) * RDC;
   }
else if( map_projection == PROJ_SAT ){
   mapuntran( x,y,&lat,&lon );
   xoff = sin( dir*DRC );
   yoff = cos( dir*DRC );
   lon += xoff/(111.195*cos(lat*DRC));
   lat += yoff/111.195;
   maptran( lat,lon,&xoff,&yoff );
   if( xoff == MISS || yoff == MISS )
      dir = MISS;
   else
      dir = atan2(xoff-x,yoff-y)*RDC;
   }
return dir;
}

/*******************************************************************
   MAP_UVDIR_ADJ - This routine adjusts the wind direction based on
   location and map projection.
********************************************************************/
float map_uvdir_adj( u,v,x,y )

float u,v;              /* The actual direction based on true north */
float x;                /* The x location based on projection */
float y;                /* The y location based on projection */

{
float dir;

dir = 0;
if( map_projection == PROJ_PSTEREO ){
   if( x == 0 && y == 0 ){
      if( proj_lat >= 0 )
         dir = dir - proj_lon;
      else
         dir = dir + proj_lon;
      }
   else if( proj_lat >= 0 )
      dir = dir - atan2( x,-y ) * RDC;
   else
      dir = dir + atan2( x,y ) * RDC;
   }
return dir;
}

/*******************************************************************
   MAP_DIR - This routine adjusts the wind direction based on
   location and map projection.
********************************************************************/
float map_dir( x1,y1,x2,y2 )

float x1,y1;             /* The first point */
float x2,y2;             /* The second point */

{
float dir;

if( map_projection == PROJ_PSTEREO ){
   dir = atan2(x2-x1,y2-y1)*RDC;
   if( proj_lat >= 0 )
      dir = dir + atan2( x1,-y1 ) * RDC;
   else
      dir = dir - atan2( x1,y1 ) * RDC;
   }
else
   dir = atan2(x2-x1,y2-y1)*RDC;
return dir;
}

/*******************************************************************
   DOMAIN_INIT - This routine initializes the plot domain parameters.
********************************************************************/
void domain_init( domain )

DOM_params *domain;      /* The domain */

{
domain->proj = MISS;
domain->lat[0] = MISS;
domain->lat[1] = MISS;
domain->lon[0] = MISS;
domain->lon[1] = MISS;
domain->dx = MISS;
domain->dy = MISS;
domain->nx = MISS;
domain->ny = MISS;
domain->param[0] = MISS;
domain->param[1] = MISS;
domain->param[2] = MISS;
domain->param[3] = MISS;
domain->param[4] = MISS;
}
