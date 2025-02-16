!
!-----------------------------------------------------------------------
!
!  ADAS Sizing parameters
!
!  Note in the following, nvar_max is the max of
!  nvar_sng, nvar_ua and nvar_anx, and is used for the size of the
!  analysis work arrays.
!
!  ntime is the number of time levels of obs read-in.
!  ntime must be at least 2 to accomodate time consistency check.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nsrc_sng,nsrc_ua,nsrc_rad,nsrc_ret
  PARAMETER (nsrc_sng=9,    & ! number of sources of single-lvl data
         nsrc_ua=5,         & ! number of sources for multi-lvl data
         nsrc_rad=5,        & ! number of sources of radar data
         nsrc_ret=1)          ! number of sources of retrieval data
!
  INTEGER :: ntime,nz_ua,nz_tab,nvar_sng,nvar_ua,                       &
          nvar_radin,nvar_rad_adas,nvar_retin,nvar_ret,                 &
          nvar_anx_adas
  PARAMETER (ntime=2,       & ! number of time levels of data
         nz_ua=200,         & ! max number of vertical levels in sounding
         nz_tab=40,         & ! max number of levels in obs err table
         nvar_sng=14,       & ! number of variables for single-lvl data
         nvar_ua=5,         & ! number of variables for multi-lvl data
         nvar_radin=4,      & ! number of variables for radar data input
         nvar_rad_adas=3,   & ! number of variables for radar data anx
         nvar_retin=5,      & ! number of variables for retrieval data input
         nvar_ret=5,        & ! number of variables for retrieval data anx
         nvar_anx_adas=5)          ! number of analysis variables

  INTEGER, PARAMETER :: nvar_rad_3dvar = 4
  INTEGER, PARAMETER :: nvar_anx_3dvar = 6
  INTEGER, PARAMETER :: mx_nvar_anx = MAX(nvar_anx_adas,nvar_anx_3dvar)

  INTEGER :: mx_sng_file,mx_ua_file,mx_rad,mx_ret
  PARAMETER( mx_sng_file=6,      & ! max number of single-level data files
             mx_ua_file=5,       & ! max number of uppper-level data files
             mx_rad=160,         & ! max number of radar sites
             mx_ret=1)             ! max number of retrieval radar sites

  INTEGER :: mx_sng,mx_ua,mx_colrad,nz_rdr,nz_ret,mx_colret

  PARAMETER (mx_sng=8000,        & ! max number of single-lvl data points
             mx_ua=200,          & ! max number of upper level sites
             mx_colrad=80000,   & ! max number of radar columns
             nz_rdr=53 ,         & ! max number of radar vertical levels
             nz_ret=1,           & ! max number of retrieval vertical levels
             mx_colret=1)          ! max number of retrieval columns

  INTEGER :: mx_pass
  PARAMETER (mx_pass = 10)
!
!-----------------------------------------------------------------------
!
!  ADAS control parameters
!  Input via namelist in INITADAS
!
!-----------------------------------------------------------------------
!
  INTEGER :: npass,nsngfil,nuafil,nradfil,nretfil
  INTEGER :: ccatopt,spradopt
  REAL :: sprdist,wlim,zwlim,thwlim,sfcqcrng

  COMMON /adas_parm/ npass,nsngfil,nuafil,nradfil,nretfil,              &
                     sprdist,wlim,zwlim,thwlim,sfcqcrng,                &
                     ccatopt,spradopt

  INTEGER :: raduvobs,radrhobs
  REAL :: refrh,rhradobs
  INTEGER :: radistride,radkstride
  INTEGER :: radcldopt,radqvopt,radqcopt,radqropt,radptopt
  REAL :: refsat,rhrad
  REAL :: refcld,cldrad
  INTEGER :: ceilopt
  REAL :: ceilmin,dzfill,refrain,radsetrat,radreflim,radptgain

  COMMON /adas_radopt/ raduvobs,radrhobs,refrh,rhradobs,                &
                 radistride,radkstride,                                 &
                 radcldopt,radqvopt,radqcopt,radqropt,radptopt,         &
                 refsat,rhrad,                                          &
                 refcld,cldrad,ceilopt,ceilmin,dzfill,                  &
                 refrain,radsetrat,radreflim,radptgain


  INTEGER :: cloudopt
  INTEGER :: clddiag
  REAL :: range_cld
  REAL :: refthr1,refthr2
  REAL :: hgtrefthr
  REAL :: wmhr_cu
  REAL :: wmhr_sc
  REAL :: wc_st
  INTEGER :: bgqcopt
  INTEGER :: cldqvopt
  INTEGER :: cldqcopt
  INTEGER :: cldqropt
  INTEGER :: cldwopt
  INTEGER :: cldptopt
  INTEGER :: smth_opt
  REAL :: thresh_cvr
  REAL :: rh_thr1
  REAL :: cvr2rh_thr1
  REAL :: rh_thr2
  REAL :: cvr2rh_thr2
  REAL :: frac_qw_2_pt
  REAL :: frac_qc_2_lh
  REAL :: max_lh_2_pt
  REAL :: qvslimit_2_qc
  REAL :: qrlimit
  REAL :: frac_qr_2_qc
  INTEGER :: cld_files
  COMMON /cmpcld_vars/ cloudopt,clddiag,range_cld,                      &
                        refthr1,refthr2,hgtrefthr,                      &
                        wmhr_cu,wmhr_sc,wc_st,bgqcopt,                  &
                        cldqvopt,cldqcopt,                              &
                        cldqropt,cldwopt,cldptopt,thresh_cvr,           &
                        rh_thr1,cvr2rh_thr1,rh_thr2,cvr2rh_thr2,        &
                        qvslimit_2_qc,qrlimit,frac_qr_2_qc,             &
                        frac_qw_2_pt,frac_qc_2_lh,max_lh_2_pt,          &
                        smth_opt,cld_files
!
!-----------------------------------------------------------------------
!
!  ADAS cloud analysis parameters
!
!-----------------------------------------------------------------------
!
!c Bad or no data flag.

  REAL :: r_missing     ! bad or no-data flag.
  PARAMETER (r_missing=-9999.0)
!
!c The maximum possible number of stations with cloud coverage reports.
!
  INTEGER :: max_cld_snd

  PARAMETER (max_cld_snd = 3000)
!
!c The maximum possible number of cloudy grid points in 3D domain
!c (should be about max_cld_snd*nz)
!
  INTEGER :: max_obs
  PARAMETER (max_obs = 100000)
!
!c The size for extended ADAS domain (extension for searching the
!c observations just out side the ADAS domain) ->"ADAS+ domain".

  INTEGER i_perimeter
  PARAMETER (i_perimeter = 10) ! (or 0) the extended ARPS domain
!
!c The radius of influence for the radar data interpolations.

  REAL :: ri_h, ri_v
  !wdt keep or move to ri_h = 12000.0?
  PARAMETER(ri_h = 20000.0, ri_v = 2500.0)
  !PARAMETER(ri_h = 12000.0, ri_v = 3000.0)
!
!c The refrence level (M AGL) for computing the lifting condensation
!c level

  REAL :: z_ref_lcl
  PARAMETER(z_ref_lcl = 180.0)
!
  CHARACTER (LEN=256) :: sngfname(mx_sng_file)
  CHARACTER (LEN=256) :: sngtmchk(mx_sng_file)
  CHARACTER (LEN=256) :: uafname(mx_ua_file)
  CHARACTER (LEN=256) :: radfname(mx_rad)
  CHARACTER (LEN=256) :: retfname(mx_ret)

  CHARACTER (LEN=256) :: backerrfil
  CHARACTER (LEN=256) :: blackfil
  CHARACTER (LEN=256) :: sngerrfil(nsrc_sng)
  CHARACTER (LEN=256) :: uaerrfil(nsrc_ua)
  CHARACTER (LEN=256) :: reterrfil(nsrc_ret)
  CHARACTER (LEN=256) :: raderrfil(nsrc_rad)

  COMMON /adas_fname/ sngfname,sngtmchk,uafname,radfname,retfname,      &
         backerrfil,blackfil,sngerrfil,uaerrfil,reterrfil,raderrfil

  INTEGER :: iusesng(0:nsrc_sng,mx_pass)
  INTEGER :: iuseua (0:nsrc_ua,mx_pass)
  INTEGER :: iuserad(0:nsrc_rad,mx_pass)
  INTEGER :: iuseret(0:nsrc_ret,mx_pass)

  INTEGER :: ianxtyp(mx_pass)
  REAL :: xyrange(mx_pass)
  REAL :: kpvar(mx_nvar_anx)
  REAL :: kpvrsq(mx_nvar_anx)
  REAL :: zrange(mx_pass)
  REAL :: thrng(mx_pass)
  INTEGER :: trnropt(mx_pass)
  REAL :: trnrcst(mx_pass)
  REAL :: trnrng(mx_pass)

  COMMON /adas_vars/ ianxtyp,iusesng,iuseua,iuserad,iuseret,            &
                     xyrange,kpvar,kpvrsq,zrange,thrng,                 &
                     trnropt,trnrcst,trnrng

  CHARACTER (LEN=8) :: srcsng(nsrc_sng)
  CHARACTER (LEN=8) :: srcua(nsrc_ua)
  CHARACTER (LEN=8) :: srcrad(nsrc_rad)
  CHARACTER (LEN=8) :: srcret(nsrc_ret)

  COMMON /adas_srcs/ srcsng,srcua,srcrad,srcret

  CHARACTER (LEN=256) :: incdmpf

  COMMON /adas_incfil/ incdmpf


  INTEGER :: incrdmp         ! Option to write ADAS analysis increments
                             ! to a file
                             ! = 0, Don't create increment file;
                             ! = 1, write increment file
                             !      (unformatted IEEE binary);
                             ! = 3, write HDF format increment file;

  INTEGER :: incrhdfcompr    ! HDF4 compression option for incrdmp = 3.

  INTEGER :: uincdmp,vincdmp,wincdmp,                                   &
          pincdmp,ptincdmp,qvincdmp,                                    &
          qcincdmp,qrincdmp,qiincdmp,qsincdmp,qhincdmp

  COMMON /adas_incr/ incrdmp,incrhdfcompr,                              &
          uincdmp,vincdmp,wincdmp,pincdmp,ptincdmp,                     &
          qvincdmp,qcincdmp,qrincdmp,qiincdmp,qsincdmp,qhincdmp
