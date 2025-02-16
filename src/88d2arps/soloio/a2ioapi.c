#include <ctype.h>
#include <float.h>

#include "geography.h"
#include "dorade.h"

#define N_INT(f) ((f) > 0.0 ? (f) + 0.5 : (f) - 0.5)

#define NO_DATA -FLT_MAX

#define MISSING_DATA -999
#define REF 1
#define VEL 2
#define SPW 3

#define REFIPARM 0    // Parameter order in the sweep for reflectivity
#define VELIPARM 1    // Parameter order in the sweep for radial velocity
#define SPWIPARM 2    // Parameter order in the sweep for SPW

#define GOOD 1
#define BAD  0

#define TRUE  1
#define FALSE 0

#define REF_ONLY  1
#define VEL_ONLY  2
#define BOTH      3

#define VERBOSE 0

char *trim(char *);
char *strhead(char *,int);
char *strtail(char *,int);

char *radar_name;

void set_radar98( int n )
{
	return;
}

/**********************************************************************
 *
 * Split a string into two file names
 *
 *********************************************************************/

int splitFNstr(char *str, char *FNstr1, char *FNstr2) {
  char *temstr;
  int slen, iret;
  int i;

  iret = 1;

  temstr = trim(str);
  slen = strlen(temstr);

  if (slen < 44) return 0;    /* At least contains 44 characters */

  i = 0;
  while (temstr[i] != ' ' && temstr[i] != '\t'){
    FNstr1[i] = temstr[i];
    i++;
  }
  FNstr1[i] = '\0';

  while (temstr[i] == ' ' || temstr[i] == '\t'){
    i++;
  }
  slen = i;

  while (temstr[i] != ' ' && temstr[i] != '\t' && temstr[i] != '\0' && temstr[i] != '\n'){
    FNstr2[i-slen] = temstr[i];
    i++;
  }
  FNstr2[i-slen] = '\0';

  temstr = strhead(FNstr1,4);
  //if (strcmp(temstr,"swp.") && strcmp(temstr,"null") ) {
  if (!strcmp(temstr,"null")) {
    strcpy(FNstr1,"NULL");
  } else if (!strstr(FNstr1,"swp.") ) {
    iret = 0;
  } else {
    if ( !strstr(FNstr1,radar_name) ) return 0;   /* File name must contain radar name */
    temstr = strtail(FNstr1,4);
    if ( strcmp(temstr,"_dbz") ) iret = 0;
  }

  temstr = strhead(FNstr2,4);
  //if (strcmp(temstr,"swp.") && strcmp(temstr,"null") ) {
  if (!strcmp(temstr,"null")) {
    strcpy(FNstr2,"NULL");
  } else if ( !strstr(FNstr2,"swp.") ) {
    iret = 0;
  } else {
    if ( !strstr(FNstr2,radar_name) ) return 0;   /* File name must contain radar name */
    temstr = strtail(FNstr2,4);
    if ( strcmp(temstr,"_vel") ) iret = 0;
  }

  free(temstr);

  return iret;
}

/*
 * This structure stores the sweep that the application is currently working on.
 */

static struct Dorade_Sweep swp_ref, swp_vel;
static struct Dorade_Sweep *swp;
static int have_sweep;

/* = 0, Above sweep structures are still not initialized;
 * = 1, Only swp_ref is initialized;
 * = 2, Only swp_vel is initialized;
 * = 3, Both swp_ref, swp_vel are initialized;
 */

static int solo;
/* = 0, Reading 88d2arps file in DORADE format; ref and vel are in separate files.
 * = 1, Reading original SOLO files.
 */

static char *refFiles[80], *velFiles[80];      /* Sweep file names */
static int  swpCT, swpFP;                      /* Sweep file pointer */
static int  rayNo = -1;

/**********************************************************************
 *
 * Set radar name
 *
 *********************************************************************/

int set_radar_name(char *rname_ptr, int slen) {

  //slen = strlen(rname_ptr);
  radar_name = (char *)malloc(slen+1);

  strncpy(radar_name,rname_ptr, slen);
  radar_name[slen] = '\0';

  return 0;
}

/**********************************************************************
 *
 * Open the file that contain file names for sweeps
 *
 *********************************************************************/

int radar_init(char *filelist_ptr) {

  FILE *fP;
  char tempstr[256], *sP1, *sP2;
  char *retstr;

  if ( !(fP = fopen(filelist_ptr, "r")) ) {
    printf("Could not open file: %s.\n",filelist_ptr);
    exit(0);
    //return -1;
  }

  /* Check the radar file type. must be the first 4 characters of the file */
  fgets(tempstr,256,fP);
  sP1 = strhead(tempstr,4);
  if (!strcmp(sP1,"solo")) solo = 1;
  free(sP1);

  /* Reading file names into the string arrays */
  sP1 = (char *)malloc(129);
  sP2 = (char *)malloc(129);
  swpCT = 0;

  retstr = fgets(tempstr,256,fP);
  //printf("reading in %s. %s\n",tempstr,retstr);
  //if (retstr == NULL) break;
  while (!feof(fP)) {
    //fscanf(flP,"%s",tempstr);

    if (solo) {   // read one file only
      int i;
      int j = strlen(tempstr)-1;
      while (tempstr[j] == '\n' || tempstr[j] == '\r' || tempstr[j] == '\t') j--;

      for(i=0;i<=j;i++) sP1[i] = tempstr[i];
      sP1[i] = '\0';

      if ( !strstr(sP1,radar_name) ) goto readagain;
                                    /* File name must contain radar name */
      refFiles[swpCT] = sP1;
      sP1 = (char *)malloc(129);
      swpCT++;
    } else {
      //printf("tempstr = %s\n",tempstr);
      if (splitFNstr(tempstr,sP1,sP2)) {
        refFiles[swpCT] = sP1;
        velFiles[swpCT] = sP2;
        swpCT++;
        sP1 = (char *)malloc(129);
        sP2 = (char *)malloc(129);
      }
    }
    //fgetc(flP);
readagain:
    fgets(tempstr,256,fP);
  }
  free(sP1);
  free(sP2);

  if (swpCT <= 0) {
    printf("== SOLOIO == No valid sweep file was found in file %s.\n",filelist_ptr);
    exit(1);
  }

  printf("== SOLOIO == %d sweep files were found in file \"%s\" and they are:\n",
          swpCT, filelist_ptr);
  int n;
  for (n=0; n<swpCT; n++) {
    printf("   %d: %s",n+1,refFiles[n]);
    if (solo) printf("\n");
    else printf("\t%s\n",velFiles[n]);
  }

  fclose(fP);

  swpFP = 0;
  have_sweep = 0;

  return 0;
}

/**********************************************************************
 *
 * Read in current sweep file
 *
 *********************************************************************/

int read_sweep() {

  FILE *in;
  int ret;

  if ( swpFP >= swpCT ) {
    have_sweep = FALSE;
    return BAD;
  }

  if (have_sweep) {
    if (have_sweep == 1 || have_sweep == 3) Dorade_FreeSweep(&swp_ref);
    if (have_sweep == 2 || have_sweep == 3) Dorade_FreeSweep(&swp_vel);
  }
  have_sweep = FALSE;

  Dorade_InitSweep(&swp_ref);
  if (!solo) Dorade_InitSweep(&swp_vel);

  if ( strcmp(refFiles[swpFP],"NULL") ) {

    if ( !(in = fopen(refFiles[swpFP], "r")) )  {
        printf("Could not open file: %s.\n",refFiles[swpFP]);
        //ErrMsg_Append("Could not open file: ");
        //ErrMsg_Append(filelist_ptr);
        //ErrMsg_Append(".\n");
	      return 0;
    }

    if ( !Dorade_ReadSweep(in, &swp_ref) ) {
      ErrMsg_Append("Could not read sweep from ");
      ErrMsg_Append(refFiles[swpFP]);
      ErrMsg_Append(".\n");
      return 0;
    }
    have_sweep = REF_ONLY;
    swp = &swp_ref;

    fclose(in);

    if (VERBOSE) Dorade_ToASCII(stdout, &swp_ref);
  }

  if (!solo &&  strcmp(velFiles[swpFP],"NULL") ) {

    if ( !(in = fopen(velFiles[swpFP], "r")) ) {
        printf("Could not open file: %s.\n",velFiles[swpFP]);
        //ErrMsg_Append("Could not open file: ");
        //ErrMsg_Append(filelist_ptr);
        //ErrMsg_Append(".\n");
        return 0;
    }

    if ( !Dorade_ReadSweep(in, &swp_vel) ) {
      ErrMsg_Append("Could not read sweep from ");
      ErrMsg_Append(velFiles[swpFP]);
      ErrMsg_Append(".\n");
      return 0;
    }

    if (have_sweep == REF_ONLY) {
      have_sweep = BOTH;
    } else {
      have_sweep = VEL_ONLY;
      swp = &swp_vel;
    }

    fclose(in);
    if (VERBOSE) Dorade_ToASCII(stdout, &swp_vel);
  }

  swpFP++;

  return GOOD;
}

/**********************************************************************
 *
 * Get radar altitude
 *
 *********************************************************************/

int get_altitude() {

  int istat = GOOD;

  if ( !have_sweep ) istat = read_sweep();

  if (istat == GOOD) {
    return N_INT(swp->radar_altitude*1000);   /* from km to meter */
  } else {
    return MISSING_DATA;
  }
}

/**********************************************************************
 *
 * Get radar latitude
 *
 *********************************************************************/

int get_latitude() {

  int istat = GOOD;

  if (!have_sweep) istat = read_sweep();

  GeoPt radarLoc;
  double lat;

  if (istat == GOOD) {
    radarLoc = swp->radar_location;
  } else {
    return MISSING_DATA;
  }

  lat = AngleToDeg(radarLoc.lat);
  return N_INT(lat*1.E5);
}

/**********************************************************************
 *
 * Get radar longitude
 *
 *********************************************************************/

int get_longitude() {

  int istat = GOOD;

  if (!have_sweep) istat = read_sweep();

  GeoPt radarLoc;
  double lon;

  if (istat == GOOD) {
    radarLoc = swp->radar_location;
  } else {
    return MISSING_DATA;
  }

  lon = AngleToDeg(radarLoc.lon);

  //lon = -1*lon;

  return N_INT(lon*1.E5);
}

/**********************************************************************
 *
 * Get Fixed angle
 *
 *********************************************************************/

int get_fixed_angle() {

  int iret;
  double angle;

  if (have_sweep) {
    if (rayNo >= 0) {    /* ray fixed angle */
      struct Dorade_RayHdr *rayHdrP;
      rayHdrP = swp->rayHdrP+rayNo;
      angle = AngleToDeg(rayHdrP->elevation);
    } else {
      angle = AngleToDeg(swp->fixed_angle);
    }
    iret = N_INT(angle*100);
  } else {
    return MISSING_DATA;
  }

  return iret;
}

/**********************************************************************
 *
 * Get scan
 *
 *********************************************************************/

int get_scan() {

  int iret;

  if (have_sweep) {
    //iret = swp->vol_num;
    iret = 1;
  } else {
    iret = MISSING_DATA;
  }

  return iret;
}

/**********************************************************************
 *
 * Get year
 *
 *********************************************************************/

int get_year() {

  struct GeoTime_Jul rtime;
  struct GeoTime_Cal caltime;
  struct Dorade_RayHdr *rhd;
  int iret = GOOD;

  if (!have_sweep || rayNo <0 || rayNo >= swp->n_rays)
    iret = read_sweep();

  if (iret == GOOD) {
    rhd = swp->rayHdrP+rayNo;
    rtime = rhd->time;
    caltime = GeoTime_JulToCal(rtime);
    iret = caltime.year-1900;
    if (iret >= 100) iret = iret-100;
  }
  return iret;
}

/**********************************************************************
 *
 * Get Month
 *
 *********************************************************************/

int get_month() {

  struct GeoTime_Jul rtime;
  struct GeoTime_Cal caltime;
  struct Dorade_RayHdr *rhd;
  int iret = GOOD;

  if (!have_sweep || rayNo <0 || rayNo >= swp->n_rays)
    iret = read_sweep();

  if (iret == GOOD) {
    rhd = swp->rayHdrP+rayNo;
    rtime = rhd->time;
    caltime = GeoTime_JulToCal(rtime);
    iret = caltime.month;
  }

  return iret;
}

/**********************************************************************
 *
 * Get day
 *
 *********************************************************************/

int get_day() {

  struct GeoTime_Jul rtime;
  struct GeoTime_Cal caltime;
  struct Dorade_RayHdr *rhd;
  int iret = GOOD;

  if (!have_sweep || rayNo <0 || rayNo >= swp->n_rays)
    iret = read_sweep();

  if (iret == GOOD) {
    rhd = swp->rayHdrP+rayNo;
    rtime = rhd->time;
    caltime = GeoTime_JulToCal(rtime);
    iret = caltime.day;
  }
  return iret;
}

/**********************************************************************
 *
 * Get hour
 *
 *********************************************************************/

int get_hour() {

  struct GeoTime_Jul rtime;
  struct GeoTime_Cal caltime;
  struct Dorade_RayHdr *rhd;
  int iret = GOOD;

  if (!have_sweep || rayNo <0 || rayNo >= swp->n_rays)
    iret = read_sweep();

  if (iret == GOOD) {
    rhd = swp->rayHdrP+rayNo;
    rtime = rhd->time;
    caltime = GeoTime_JulToCal(rtime);
    iret = caltime.hour;
  }
  return iret;
}

/**********************************************************************
 *
 * Get minute
 *
 *********************************************************************/

int get_min() {

  struct GeoTime_Jul rtime;
  struct GeoTime_Cal caltime;
  struct Dorade_RayHdr *rhd;
  int iret = GOOD;

  if (!have_sweep || rayNo <0 || rayNo >= swp->n_rays)
    iret = read_sweep();

  if (iret == GOOD) {
    rhd = swp->rayHdrP+rayNo;
    rtime = rhd->time;
    caltime = GeoTime_JulToCal(rtime);
    iret = caltime.minute;
  }
  return iret;
}

/**********************************************************************
 *
 * Get second
 *
 *********************************************************************/

int get_sec() {

  struct GeoTime_Jul rtime;
  struct GeoTime_Cal caltime;
  struct Dorade_RayHdr *rhd;
  int iret = GOOD;

  if (!have_sweep || rayNo <0 || rayNo >= swp->n_rays)
    iret = read_sweep();

  if (iret == GOOD) {
    rhd = swp->rayHdrP+rayNo;
    rtime = rhd->time;
    caltime = GeoTime_JulToCal(rtime);
    iret = N_INT(caltime.second);
  }
  return iret;
}

/**********************************************************************
 *
 * Get vcp
 *
 *********************************************************************/

int get_vcp() {
  int iret = 0;
  char *result,*end;

  if (!solo && have_sweep) {
    result = strstr(swp->proj_name,"vcp=");
    if (result == NULL) iret = 0;
    else iret = atoi(result+4);
  }

  return iret;
}

/**********************************************************************
 *
 * Get nyquist
 *
 *********************************************************************/

int get_nyquist() {

  int iret = MISSING_DATA;
  float v_nyq;

  if (solo) {
    v_nyq = swp_ref.eff_unamb_vel;
    iret = N_INT(v_nyq*100);
  } else if (have_sweep == VEL_ONLY || have_sweep == BOTH) {
    v_nyq = swp_vel.eff_unamb_vel;
    iret = N_INT(v_nyq*100);
  }

  if (iret == 0) iret = MISSING_DATA;

  return iret;
}

/**********************************************************************
 *
 * Get tilt
 *
 *********************************************************************/

int get_tilt() {

  int iret;

  if (have_sweep) {
    iret = swp->sweep_num;
  } else {
    iret = MISSING_DATA;
  }

  return iret;
}

/**********************************************************************
 *
 * Get azi
 *
 *********************************************************************/

int get_azi() {

  int iret;
  double angleD;

  if ( !have_sweep || (rayNo < 0 || rayNo >= swp->n_good_rays) ) {
    iret = read_sweep();
    if (iret != GOOD) return MISSING_DATA;
  }

  struct Dorade_RayHdr *rayHdrP;
  rayHdrP = swp->rayHdrP+rayNo;
  angleD = AngleToDeg(rayHdrP->azimuth);

  iret = N_INT(angleD*100);

  /*
    if (rayNo == 60) {
      printf("ray %d azim = %d,%d,%f,%d\n",60,
          (swp_ref.rayHdrP+rayNo)->azimuth,(swp_ref.rayHdrP+rayNo)->azimuth,angleD,iret);
      exit(0);
    }
  */


  return iret;
}


/**********************************************************************
 *
 * read radial
 *
 *********************************************************************/

#define RR_DATA         0
#define RR_EOD          1

int read_radial() {

  int iret = RR_DATA;

  if (have_sweep) rayNo++;

  if (!have_sweep || rayNo >= swp->n_good_rays) {
    iret = read_sweep();
    if (iret == GOOD) {
      rayNo = 0;
      iret = RR_DATA;
    } else {
      rayNo = -1;
      iret = RR_EOD;
    }
  }

  //printf("  Reading rayNo = %d.\n",rayNo);
  return iret;
}

/**********************************************************************
 *
 * Get number of gates
 *
 *********************************************************************/

int get_number_of_gates(int indx) {

  int iret = 0;

  switch (indx) {
    case REF:
      if (have_sweep == REF_ONLY || have_sweep == BOTH) iret = swp_ref.n_cells;
      break;
    case VEL:
    case SPW:
      if (solo) {
        iret = swp_ref.n_cells;
      } else {
        if (have_sweep == VEL_ONLY || have_sweep == BOTH )
          iret = swp_vel.n_cells;
      }
      break;
    default:
      iret = 0;
      break;
  }

  return iret;
}

/**********************************************************************
 *
 * Get fist gate
 *
 *********************************************************************/

int get_first_gate(int indx) {

  int iret;
  float tmp = MISSING_DATA;

  switch (indx) {
    case REF:
      if (have_sweep == REF_ONLY || have_sweep == BOTH) tmp = *swp_ref.distP;
      break;
    case VEL:
    case SPW:
      if (solo) {
        tmp = *swp_ref.distP;
      } else {
        if (have_sweep == VEL_ONLY || have_sweep == BOTH )
          tmp = *swp_vel.distP;
      }
      break;
    default:
      tmp = MISSING_DATA;
      break;
  }

  iret = N_INT(tmp);

  return iret;
}

/**********************************************************************
 *
 * Get gate spacing
 *
 *********************************************************************/

int get_gate_spacing(int indx) {

  int iret;
  float tmp;

  switch (indx) {
    case REF:
      if (have_sweep == REF_ONLY || have_sweep == BOTH )
        tmp =( *(swp_ref.distP+1)-*(swp_ref.distP) );
      break;
    case VEL:
    case SPW:
      if (solo) {
        tmp =( *(swp_ref.distP+1)-*(swp_ref.distP) );
      } else {
        if (have_sweep == VEL_ONLY || have_sweep == BOTH )
          tmp =( *(swp_vel.distP+1)-*(swp_vel.distP) );
      }
      break;
    default:
      tmp = MISSING_DATA;
      break;
  }

  iret = N_INT(tmp);

  return iret;
}

/**********************************************************************
 *
 * Get field number
 *
 *********************************************************************/

int get_field_num(char *fName) {

  int iret = -1;

  if (!strcmp(fName,"DBZ")) {
    iret = REF;
  } else if (!strcmp(fName,"VEL")) {
    iret = VEL;
  } else if (!strcmp(fName,"SPW")) {
    iret = SPW;
  }

  return iret;
}

/**********************************************************************
 *
 * Get status
 *
 *********************************************************************/

int get_status(int indx) {

  int istat = BAD;

  if (swpFP <= swpCT) {
    if (indx == REF) {         /* Reflectivity */
      if (have_sweep == 1 || have_sweep == 3) {
        istat = GOOD;
      }
    } else if (indx == VEL || indx == SPW) {  /* velocity */
      if (have_sweep == 1) {
        if (solo) istat = GOOD;
      } else if (have_sweep == 2 || have_sweep == 3 ) {
        istat = GOOD;
      }
    }
  }
  return istat;

//kzhao temporayly modified
/*
  if(solo)
  {
    if (swpFP == 1 || swpFP == 3)
    {
      if(indx==VEL || indx == SPW) istat=FALSE;
    }

    if (swpFP ==2 || swpFP ==4)
    {
      if(indx ==REF) istat=FALSE;
    }
  }
  return istat;
*/

}

/**********************************************************************
 *
 * Get data field
 *
 *********************************************************************/

int get_data_field(int indx, float *out_ptr, int n_gates) {

  int iret = BAD;
  int nmax_gates;
  int i, iparm;
  struct Dorade_Sweep *swp;
  float tmp;

  switch (indx) {
    case REF:
      if (have_sweep == REF_ONLY || have_sweep == BOTH) {
        swp = &swp_ref;
        iparm = REFIPARM;
        iret = GOOD;
      }
      break;
    case VEL:
      if (solo) {
        swp = &swp_ref;
        iparm = VELIPARM;
        iret = GOOD;
      } else {
        if (have_sweep == VEL_ONLY || have_sweep == BOTH ) {
          swp = &swp_vel;
          iparm = VELIPARM - 1;
          iret = GOOD;
        }
      }
      break;
    case SPW:
      if (solo) {
        swp = &swp_ref;
        iparm = SPWIPARM;
        iret = GOOD;
      } else {
        if (have_sweep == VEL_ONLY || have_sweep == BOTH ) {
          swp = &swp_vel;
          iparm = SPWIPARM - 1;
          iret = GOOD;
        }
      }
      break;
    default:
      iret = BAD;
      break;
  }

  if (iret == GOOD) {
    nmax_gates = n_gates>swp->n_cells ? swp->n_cells : n_gates;

    for (i=0;i<nmax_gates;i++) {
      tmp = *(swp->dat[iparm][rayNo]+i);
      /*
      if (rayNo == 26 && i==40) {
        printf("read in %f.\n",tmp);
        exit(0);
      }
      */
      if ( tmp == NO_DATA) {
         *(out_ptr+i) = MISSING_DATA;
      } else {
         *(out_ptr+i) = tmp;
      }
    }
  } else {
    for (i=0; i< n_gates; i++) {
      *(out_ptr+i) = MISSING_DATA;
    }
  }

  return iret;
}

int get_restart_flag()
{
  return( 0 );
}

/**********************************************************************
 *
 * Clean up at last
 *
 *********************************************************************/

void solo_clean_up(void) {
    //printf("== SOLO == Cleaning up SOLO IO.\n");
    //
    if (have_sweep >= REF_ONLY) Dorade_FreeSweep(&swp_ref);
    if (have_sweep >= VEL_ONLY) Dorade_FreeSweep(&swp_vel);

    ErrMsg_Destroy();
    Dorade_DestroyLib();
}
