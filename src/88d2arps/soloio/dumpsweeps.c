
#include <ctype.h>
#include <float.h>

#include "geography.h"
#include "dorade.h"

#define NO_DATA -FLT_MAX
#define MISSING_DATA -999.0

#define VERBOSE 0

/**********************************************************************
 *
 * inital NEXRAD sweep structure (ref and vel)
 *
 * ********************************************************************/

void init_swp_nexrad(struct Dorade_Sweep *swp,
          char *proj_name,char *radar_name,
          GeoPt radarLoc, float radar_alt_m,int i_vcp,
          struct GeoTime_Jul time,
          int n_parms, char *parmNames[], char *parmDesc[], char *parmUnits[],
          short i_scan, Angle fixed_angle, int i_tilt,
          int n_max_rays, int n_max_gates) {

  size_t sz;
  int n,i,j;

  Dorade_InitSweep(swp);

  strcpy(swp->radar_name,radar_name);
  strcpy(swp->proj_name,proj_name);
  //sprintf(swp->proj_name,"%s 88D2ARPS vcp=%d",proj_name,i_vcp);
  swp->n_parms = n_parms;                /* for Nexrads only */
  swp->n_use_parms = n_parms;
  swp->vol_num = i_scan;

  swp->radar_type = 0;             /* it should be ground for NEXRAD */

  swp->radar_location = radarLoc;
  swp->radar_altitude = radar_alt_m;

  swp->time = time;

  swp->scan_mode = DORADE_PPI;      /* A communication with Kevin */
  swp->fixed_angle = fixed_angle;
  swp->sweep_num = i_tilt;

  swp->n_rays  = n_max_rays;         /* n_good_rays to be set later */
  swp->n_cells = n_max_gates;
  swp->n_good_rays = 0;

  /* Allocate distP */
  swp->distP = (float *)MALLOC(swp->n_cells * sizeof(float));

  /* Initialize Dorade_ParmDesc structures */
  sz = swp->n_parms * sizeof(struct Dorade_ParmDesc);
  swp->parmP = (struct Dorade_ParmDesc *)MALLOC(sz);

  struct Dorade_ParmDesc *parmP;
  for (n = 0; n < swp->n_parms; n++) {
    parmP = swp->parmP+n;
    Dorade_InitParm(parmP);
    parmP->use = 1;
    strcpy(parmP->name,parmNames[n]);
    strcpy(parmP->description,parmDesc[n]);
    strcpy(parmP->units,parmUnits[n]);
    parmP->binary_format = DD_16_BITS;
    parmP->bad_data = -999;
  }

  /* Initialize structures rayHdrP */
  struct Dorade_RayHdr *rhp, *rhp1;

  sz = swp->n_rays * sizeof(struct Dorade_RayHdr);
  swp->rayHdrP = (struct Dorade_RayHdr *)MALLOC(sz);
  for (rhp = swp->rayHdrP, rhp1 = rhp+swp->n_rays; rhp < rhp1; rhp++) {
    Dorade_InitRayHdr(rhp);
  }

  /* Allocate data array */
  swp->dat = Dorade_AllocDat(swp->n_parms,swp->n_rays,swp->n_cells);

  for(n=0;n<swp->n_parms;n++) {
    for(i=0;i<swp->n_rays;i++) {
      for(j=0;j<swp->n_cells;j++) {
        swp->dat[n][i][j] = NO_DATA;
      }
    }
  }
  //Dorade_HashParms(swp);
  return;
}

/*
 * This structure stores the sweep that the application is currently working on.
 */

static int initialized;

static struct Dorade_Sweep swp_ref, swp_vel;
static int have_sweep_ref, have_sweep_vel;

/**********************************************************************
 *
 * inital sweep structure
 *
 * ********************************************************************/

int init_solo(char *proj_name, char *radar_name,
        int i_lat, int i_lon, int i_alt,int i_vcp,
        int iyear,int imon, int iday,int ihr, int imin, double second,
        int i_scan, int i_tilt,int i_angle,
        int n_ref_gates,int n_vel_gates, int n_max_rays) {

  char *parmNames_ref[]={"DBZ"};
  char *parmDesc_ref []={"Reflectivity"};
  char *parmUnits_ref[]={"DBZ"};

  char *parmNames_vel[]={"VEL","SPW"};
  char *parmDesc_vel []={"Radial Velocity","Spectral Width"};
  char *parmUnits_vel[]={"m/s"," "};

  struct GeoTime_Cal timeCal;
  struct GeoTime_Jul timeJul;

  GeoPt radarLoc;
  float radarAlt_m;
  Angle fixed_angle;

  short vol_num;

  //printf("== SOLO == Initializing the sweep strutures.\n");

  radarLoc.lat = i_lat*10;        /* Latitude in microdegreee */
  radarLoc.lon = i_lon*-10;
  radarAlt_m   = i_alt*0.001;

  timeCal.year   = iyear ;
  timeCal.month  = imon ;
  timeCal.day    = iday ;
  timeCal.hour   = ihr ;
  timeCal.minute = imin ;
  timeCal.second = second ;
  timeJul = GeoTime_CalToJul(timeCal);

  fixed_angle = AngleFmDeg(i_angle*0.01);

  vol_num = i_scan;

  init_swp_nexrad(&swp_ref,proj_name,radar_name,radarLoc,radarAlt_m,i_vcp,
      timeJul, 1, parmNames_ref,parmDesc_ref,parmUnits_ref,
      vol_num,fixed_angle,i_tilt,
      n_max_rays,n_ref_gates);
  have_sweep_ref = 0;

  init_swp_nexrad(&swp_vel,proj_name,radar_name,radarLoc,radarAlt_m,i_vcp,
      timeJul, 2, parmNames_vel,parmDesc_vel,parmUnits_vel,
      vol_num,fixed_angle,i_tilt,
      n_max_rays,n_vel_gates);
  have_sweep_vel = 0;

  initialized = 1;

  return 0;
}

/**********************************************************************
 *
 * Main Function to fill sweep structure
 *
 * ********************************************************************/

int add_radial_to_sweep(int i_tilt, double elevin, double azimin,
               int ref_ok, int vel_ok,int n_ref_gates,int n_vel_gates,int n_spw_gates,
               float *v_nyq,
               int iyear,int imon, int iday,int ihr, int imin, double second,
               int rfrst_ref,int gsp_ref,int rfrst_vel, int gsp_vel,
               float *ref_ptr,float *vel_ptr,float *spw_ptr) {

  int n;
  Angle elev, azim;
  float *dataP, *dataP1;

  if (!initialized) {   /* Allocate space for sweep */
    ErrMsg_Append("Sweep structure is still not initialized when adding radial.\n");
    return -1;
  }

  /*
  printf("radar lat = %d\n",swp_ref.radar_location.lat);
  printf("radar lon = %d\n",swp_ref.radar_location.lon);
  printf("vol_num = %d\n",swp_ref.vol_num);
  exit(0);
  */

  /*======================== BEGIN ==================================*/

  if (ref_ok || vel_ok) {
    /* First set ray header */
    struct GeoTime_Cal timeCal;

    timeCal.year = iyear ;
    timeCal.month = imon ;
    timeCal.day = iday ;
    timeCal.hour = ihr ;
    timeCal.minute = imin ;
    timeCal.second = second ;

    azim = AngleFmDeg( azimin );
    elev = AngleFmDeg( elevin );

    //swp_ref.eff_unamb_vel = v_nyq;
    swp_vel.eff_unamb_vel = *v_nyq;

    struct Dorade_RayHdr *rhp;
    int rayno;

    /* Reflectivity ray header */
    rayno = swp_ref.n_good_rays;

    rhp = swp_ref.rayHdrP + rayno;
    rhp->good = 1;

    rhp->time =  GeoTime_CalToJul(timeCal);

    rhp->azimuth   = azim;
    rhp->elevation = elev;

    rhp->ray_status = 0;

    /* Velocity ray header */
    rhp = swp_vel.rayHdrP + rayno;
    rhp->good = 1;

    rhp->time =  GeoTime_CalToJul(timeCal);

    rhp->azimuth   = azim;
    rhp->elevation = elev;

    rhp->ray_status = 0;


    /* Add REF ray data */

    if (ref_ok) {

      /* Initialize cell range vector */
      if (! have_sweep_ref) {
        float *fP;
        fP = swp_ref.distP;
        for (n = 0; n < swp_ref.n_cells; n++) {
          *(fP+n) = rfrst_ref + n * gsp_ref;
        }
      }

      /*
      printf("== SOLO == Adding a REF radial to tilt %d, elevation = %d, azi = %d.\n",
              i_tilt,elev,azim);
      */

      /*
      printf("== SOLO == Sample data (%d): %d %d %d\n%f %f %f %f.\n",
        n_ref_gates,swp_ref.n_parms,swp_ref.n_rays,swp_ref.n_cells,
        *(ref_ptr+20),*(ref_ptr+40),
        *(ref_ptr+60),*(ref_ptr+80));
      */

      dataP = swp_ref.dat[0][rayno];
      for ( n = 0; n < n_ref_gates; n++ ){
        if (*(ref_ptr+n) != MISSING_DATA && n < swp_ref.n_cells) {
          *(dataP+n) = *(ref_ptr+n);
        }
      }

      /*
      printf("== SOLO == Sample data (%d): %d %d %d\n%f %f %f %f.\n",
        n_ref_gates,swp_ref.n_parms,swp_ref.n_rays,swp_ref.n_cells,
        *(swp_ref.dat[0][rayno]+20),*(swp_ref.dat[0][rayno]+40),
        *(swp_ref.dat[0][rayno]+60),*(swp_ref.dat[0][rayno]+80));
      if (swp_ref.n_good_rays == 26) exit(0);
      */

      have_sweep_ref++;     /* ray counter */
    }

    if (vel_ok) {

      /* Initialize cell range vector */
      if (! have_sweep_vel) {
        float *fP;
        fP = swp_vel.distP;
        for (n = 0; n < swp_vel.n_cells; n++) {
          *(fP+n) = rfrst_vel + n * gsp_vel;
        }
      }

      /*
      printf("== SOLO == Adding a VEL/SPW radial to tilt %d, elevation = %d, azi = %d.\n",
              i_tilt,elev,azim);
     */

      dataP  = swp_vel.dat[0][rayno];
      dataP1 = swp_vel.dat[1][rayno];
      for ( n=0; n < n_vel_gates; n++){
        if (n < swp_vel.n_cells) {
          if (*(vel_ptr+n) != MISSING_DATA ) *(dataP+n)  = *(vel_ptr+n);
          if (*(spw_ptr+n) != MISSING_DATA ) *(dataP1+n) = *(spw_ptr+n);
        }
      }

      /*
      printf("== SOLO == Sample data (%d): %d %d %d\n%f %f %f %f.\n",
        n_vel_gates,swp_vel.n_parms,swp_vel.n_rays,swp_vel.n_cells,
        *(swp_vel.dat[0][rayno]+20),*(swp_vel.dat[0][rayno]+40),
        *(swp_vel.dat[0][rayno]+60),*(swp_vel.dat[0][rayno]+80));
     */

      /* at last, add counter of good ray of this sweep */
      have_sweep_vel++;
    }

    swp_ref.n_good_rays++;
    swp_vel.n_good_rays++;
  }

  /*======================== END ====================================*/

  return 0;
}

/**********************************************************************
 *
 * Write a NEXRAD sweep
 *
 * ********************************************************************/

int write_sweep_nexrad(struct Dorade_Sweep *swp) {

  char  *drd_fnm = NULL;       /* Name of DORADE file to write */
  size_t drd_fsz = 0;          /* Allocation at drd_fnm */

  FILE *out;

  drd_fnm = Dorade_SweepFName(drd_fnm, &drd_fsz, swp);

  //printf("== SOLO == Writing to sweep file %s.\n",drd_fnm);

  if ( !(out = fopen(drd_fnm, "w")) ) {
     ErrMsg_Append("Could not open sweep file ");
     ErrMsg_Append(drd_fnm);
     ErrMsg_Append(".\n");
     return -1;
  }

  //printf("== SOLO == Writing to sweep file %s.\n",drd_fnm);

  if ( !Dorade_WriteSweep(out, swp) ) {
     ErrMsg_Append("Could not write sweep to ");
     ErrMsg_Append(drd_fnm);
     ErrMsg_Append(".\n");
     return -2;
  }

  FREE(drd_fnm);
  fclose(out);
  Dorade_FreeSweep(swp);

  return 0;
}

/**********************************************************************
 *
 * Funciton to write a sweep file
 *
 * ********************************************************************/

int write_sweep() {

  int istatus;

  if (!initialized) {
    ErrMsg_Append("Sweep structure is still not initialized when calling write_sweep.\n");
    return -3;
  }

  istatus = 0;

  if (have_sweep_ref) {
    if (VERBOSE) Dorade_ToASCII(stdout, &swp_ref);
    istatus = write_sweep_nexrad( &swp_ref );
  }
  if (istatus != 0) return istatus;

  if (have_sweep_vel) {
    if (VERBOSE) Dorade_ToASCII(stdout, &swp_vel);
    istatus = write_sweep_nexrad( &swp_vel );
  }
  if (istatus != 0) return istatus;

  initialized = 0;

  return istatus;
}

/**********************************************************************
 *
 * clean up before exit
 *
 * ********************************************************************/

void clean_up(void) {
  //printf("== SOLO == Cleaning up SOLO IO.\n");
  ErrMsg_Destroy();
  Dorade_DestroyLib();
}
