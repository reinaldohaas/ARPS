/*
 *  dumpsweeps_f.c
 *
 *  This is a Fortran interface version.
 *
 */

#ifdef UNDERSCORE
#define INITSOLO_F init_solo_f_
#define ADDSWEEP_F add_radial_to_sweep_f_
#define WRITSOLO_F write_sweep_f_
#else
#define INITSOLO_F init_solo_f
#define ADDSWEEP_F add_radial_to_sweep_f
#define WRITSOLO_F write_sweep_f
#endif

void INITSOLO_F(char *proj_name,int *plen, char *radar_name,int* rlen,
        float* lat, float* lon, float* alt,int* vcp,
        int* iyear,int* imon, int* iday,int* ihr, int* imin, int* isecond,
        int* i_scan, int* i_tilt,int* i_angle,
        int* n_ref_gates,int* n_vel_gates, int* n_max_rays, int* istatus)
{
  int i_lat, i_lon, i_alt;
  int i_vcp;
  double second;

  i_lat = (int)(*lat*100000);
  i_lon = (int)(*lon*100000);
  i_alt = (int)(*alt);
  i_vcp = *vcp;
  second = (double)(*isecond);
  printf("==solo== plen=%d\n",*plen);
  *(proj_name+(*plen))  = '\0';
  *(radar_name+(*rlen)) = '\0';

  *istatus = init_solo(proj_name,radar_name,i_lat, i_lon, i_alt,i_vcp,
        *iyear, *imon, *iday, *ihr, *imin, second,
        *i_scan, *i_tilt, *i_angle,
        *n_ref_gates,*n_vel_gates, *n_max_rays);
}

void ADDSWEEP_F(int* i_tilt, float *eleva_ptr, float *azi_ptr,
               int* ref_ok, int* vel_ok,int* n_ref_gates,int* n_vel_gates,int* n_spw_gates,
               float *v_nyq,
               int* iyear,int* imon, int* iday,int* ihr, int* imin, int* isecond,
               int* rfrst_ref,int* gsp_ref,int* rfrst_vel, int* gsp_vel,
               float *ref_ptr,float *vel_ptr,float *spw_ptr, int* istatus)
{
  double second;
  double elevin, azimin;

  second = (double)(*isecond);
  elevin = (double)(*eleva_ptr);
  azimin = (double)(*azi_ptr);

  *istatus = add_radial_to_sweep(*i_tilt, elevin, azimin,
               *ref_ok, *vel_ok, *n_ref_gates, *n_vel_gates, *n_spw_gates,
               v_nyq,
               *iyear,*imon,*iday,*ihr,*imin,second,
               *rfrst_ref,*gsp_ref,*rfrst_vel, *gsp_vel,
               ref_ptr,vel_ptr,spw_ptr);
}

void WRITSOLO_F (int* istatus)
{
  *istatus = write_sweep();
}
