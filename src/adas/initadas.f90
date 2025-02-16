!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE INITADAS                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE initadas(namelist_filename)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Read in analysis variables in namelist format from standard input.
!
!  AUTHOR:
!  Keith Brewster, CAPS, January, 1996
!
!  MODIFICATION HISTORY:
!
!  Added read of variables previously read by calling initpara
!  Added new src and iuse variables.
!  Keith Brewster, CAPS, October, 1997
!
!  Added lines to initialize mgrid and nstgrid which are used
!  by the GRIB writer.
!  Keith Brewster, CAPS, March, 1998
!
!  Keith Brewster, CAPS, April 27,1998
!  Added new parameters to output namelist to match ARPS
!  input file options.
!
!  2000/05/10 (Gene Bassett)
!     Merged the ADAS input file back into the ARPS input file.
!
!  2000-05-18 (Gene Bassett)
!     Moved hydrostatic and wind adjustment parameters to an
!     adjust namelist block.
!
!  2003-09-03 (Steve Leyton, CAPS)
!     Added MPI capabilities
!
!  2004/02/10 (Dan Weber, CAPS)
!     Added code for the root processor to read in the namelists
!     (following ARPS convention) and send the data to the other
!     processors via the mpupdate funtion.  Also wrapped each write to
!     standard out with myproc=0 if statements, thus only the root
!     processor will write to the standard output file.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INCLUDE 'adas.inc'
  INCLUDE 'adassat.inc'
  INCLUDE 'adjust.inc'
  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!  Misc local variables
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,isat
  INTEGER :: lenstr
  LOGICAL :: iexist
  CHARACTER (LEN=19) :: initime  ! Real time in form of 'year-mo-dy:hr:mn:ss'

  INTEGER :: err_no
  DATA err_no /0/

  CHARACTER(LEN=*), INTENT(IN) :: namelist_filename
!
!-----------------------------------------------------------------------
!
!  ADAS namelists
!
!-----------------------------------------------------------------------
!
  NAMELIST /adas_const/ npass,sprdist,wlim,zwlim,thwlim,                &
                 ccatopt,spradopt
  NAMELIST /adjust/ hydradj,wndadj,obropt,obrzero
  NAMELIST /adas_radaropt/ raduvobs,radrhobs,refrh,rhradobs,            &
                 radistride,radkstride,                                 &
                 radcldopt,radqvopt,radqcopt,radqropt,radptopt,         &
                 refsat,rhrad,                                          &
                 refcld,cldrad,ceilopt,ceilmin,dzfill,                  &
                 refrain,radsetrat,radreflim,radptgain
  NAMELIST /adas_typ/   ianxtyp
  NAMELIST /adas_range/  sfcqcrng,xyrange
  NAMELIST /adas_kpvar/ kpvar
  NAMELIST /adas_zrange/ zrange
  NAMELIST /adas_thrng/  thrng
  NAMELIST /adas_trnrng/  trnropt,trnrcst,trnrng
  NAMELIST /adas_backerf/ backerrfil
  NAMELIST /adas_sng/ nsngfil,sngfname,sngtmchk,blackfil,               &
                      srcsng,sngerrfil,iusesng
  NAMELIST /adas_ua/ nuafil,uafname,srcua,uaerrfil,iuseua
  NAMELIST /adas_radar/ nradfil,radfname,srcrad,raderrfil,iuserad
  NAMELIST /adas_retrieval/                                             &
                        nretfil,retfname,srcret,reterrfil,iuseret
  NAMELIST /adas_cloud/ cloudopt,clddiag,range_cld,                     &
                        refthr1,refthr2,hgtrefthr,                      &
                        wmhr_cu,wmhr_sc,wc_st,bgqcopt,                  &
                        cldqvopt,cldqcopt,                              &
                        cldqropt,cldwopt,cldptopt,thresh_cvr,           &
                        rh_thr1,cvr2rh_thr1,rh_thr2,cvr2rh_thr2,        &
                        qvslimit_2_qc,qrlimit,frac_qr_2_qc,             &
                        frac_qw_2_pt,frac_qc_2_lh,max_lh_2_pt,          &
                        smth_opt,                                       &
                        nvisfiles,vis_fname,nirfiles,ir_fname,cld_files, &
                        viscalname,ircalname
  NAMELIST /incr_out/ incdmpf,incrdmp,incrhdfcompr,                     &
                      uincdmp,vincdmp,wincdmp,                          &
                      pincdmp,ptincdmp,qvincdmp,                        &
                      qcincdmp,qrincdmp,qiincdmp,qsincdmp,qhincdmp

  INTEGER :: unum
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  unum = 0
  IF (myproc == 0) THEN
    IF (LEN_TRIM(namelist_filename) <= 0 .OR. namelist_filename == ' ') THEN
      unum = 5
      WRITE(6,'(2(1x,a,/))') 'Waiting ADAS namelist from standard input ... ', &
                             '========================================'
    ELSE
      CALL getunit( unum )
      OPEN(unum,FILE=TRIM(namelist_filename),STATUS='OLD',FORM='FORMATTED')
      WRITE(6,'(1x,3a,/,1x,a,/)') 'Reading ADAS namelist from file - ',   &
              TRIM(namelist_filename),' ... ','========================================'
    END IF
  END IF

!-----------------------------------------------------------------------
!  Assign default values to the ADAS input variables
!-----------------------------------------------------------------------

  npass=4
  sprdist = 15000.
  wlim = 1.e-04
  zwlim = 1.e-03
  thwlim = 1.e-03
  ccatopt = 1
  spradopt = 2

  hydradj = 0
  wndadj = 2
  obropt = 2
  obrzero = 12000.

  raduvobs = 0
  radrhobs = 0
  radistride = 2
  radkstride = 2
  refrh = 20.
  rhrad = 0.90

  radcldopt = 0
  radqvopt = 0
  radqcopt = 0
  radqropt = 0
  radptopt = 0
  refcld = 20.
  cldrad = 1.0E-03
  ceilopt = 1
  ceilmin = 1500.
  dzfill = 3000.
  refrain = 30.
  radsetrat = 0.50
  radreflim = 50.
  radptgain = 1.0

  cloudopt = 0
  clddiag = 1
  range_cld = 100.0E03
  refthr1 = 25.0
  refthr2 = 10.0
  hgtrefthr = 2000.0
  wmhr_cu = 0.0005
  wmhr_sc = 0.00005
  wc_st = 0.05
  bgqcopt  = 1
  cldqvopt = 1
  cldwopt  = 0
  cldqcopt = 1
  cldqropt = 1
  cldptopt = 1

  smth_opt = 1
  thresh_cvr = 0.65
  rh_thr1 = 0.5
  cvr2rh_thr1 = 0.2
  rh_thr2 = 0.95
  cvr2rh_thr2 = 0.65
  frac_qw_2_pt = 1.0
  frac_qc_2_lh = 1.0
  max_lh_2_pt = 4.0
  qvslimit_2_qc = 1.0
  qrlimit = 0.010
  frac_qr_2_qc = 0.5


  DO i=1,mx_nvar_anx
    kpvar(i)=1.0
  END DO
!
  sfcqcrng = 100.e03
  DO i=1,mx_pass
    ianxtyp(i)=21
    xyrange(i)=100.e03
    zrange(i)=500.
    thrng(i)=3.0
    trnropt(i)=1
    trnrcst(i)=500.
    trnrng(i)=1.2
  END DO
!
  backerrfil='ruc3herr.adastab'
!
  blackfil='blacklist.sfc'
  nsngfil=0
  DO i=1,mx_sng_file
    sngfname(i)='NULL'
    sngtmchk(i)='NULL'
  END DO
  DO i=1,nsrc_sng
    srcsng(i)='NULL'
    sngerrfil(i)='NULL'
  END DO
  DO i=0,nsrc_sng
    DO j=1,npass
      iusesng(i,j)=0
    END DO
  END DO
!
  nsngfil=0
  DO i=1,mx_ua_file
    uafname(i)='NULL'
  END DO
  DO i=1,nsrc_ua
    srcua(i)='NULL'
    uaerrfil(i)='NULL'
  END DO
  DO i=0,nsrc_ua
    DO j=1,npass
      iuseua(i,j)=0
    END DO
  END DO
!
  nradfil=0
  DO i=1,mx_rad
    radfname(i)='NULL'
  END DO
  DO i=1,nsrc_rad
    srcrad(i)='NULL'
    raderrfil(i)='NULL'
  END DO
  DO i=0,nsrc_rad
    DO j=1,npass
      iuserad(i,j)=0
    END DO
  END DO
!
  nretfil=0
  DO i=1,mx_ret
    retfname(i)='NULL'
  END DO
  DO i=1,nsrc_ret
    srcret(i)='NULL'
    reterrfil(i)='NULL'
  END DO
  DO i=0,nsrc_ret
    DO j=1,npass
      iuseret(i,j)=0
    END DO
  END DO
!
  cld_files=0
  nirfiles=0
  nvisfiles=0
  DO isat=1,mx_sat
    ir_fname(isat)='NULL'
    vis_fname(isat)='NULL'
    viscalname(isat)='NULL'
    ircalname(isat)='NULL'
  END DO
!
  incdmpf='NULL'
  incrdmp = 0
  incrhdfcompr = 0
  uincdmp = 1
  vincdmp = 1
  wincdmp = 0
  pincdmp = 1
  ptincdmp= 1
  qvincdmp= 1
  qcincdmp= 1
  qrincdmp= 1
  qiincdmp= 1
  qsincdmp= 0
  qhincdmp= 0

!
!-----------------------------------------------------------------------
!  read in ADAS namelists
!-----------------------------------------------------------------------
!
  IF(myproc == 0) THEN
    READ(unum, incr_out,  END=201)
    WRITE(6,*) 'Namelist block incr_out sucessfully read.'
    201 CONTINUE
  END IF

  CALL mpupdatei(incrdmp,1)
  CALL mpupdatei(incrhdfcompr,1)
  CALL mpupdatec(incdmpf,256)
  CALL mpupdatei(uincdmp,1)
  CALL mpupdatei(vincdmp,1)
  CALL mpupdatei(wincdmp,1)
  CALL mpupdatei(pincdmp,1)
  CALL mpupdatei(ptincdmp,1)
  CALL mpupdatei(qvincdmp,1)
  CALL mpupdatei(qcincdmp,1)
  CALL mpupdatei(qrincdmp,1)
  CALL mpupdatei(qiincdmp,1)
  CALL mpupdatei(qsincdmp,1)
  CALL mpupdatei(qhincdmp,1)

  IF(myproc == 0) THEN
    READ(unum, adas_const,  END=200)
    WRITE(6,*) 'Namelist block adas_const sucessfully read.'
    200 CONTINUE
  END IF

  CALL mpupdatei(npass,1)
  CALL mpupdater(sprdist,1)
  CALL mpupdater(wlim,1)
  CALL mpupdater(zwlim,1)
  CALL mpupdater(thwlim,1)
  CALL mpupdatei(spradopt,1)
  CALL mpupdatei(ccatopt,1)

  IF(myproc == 0) THEN
    READ(unum, adjust,  END=205)
    WRITE(6,*) 'Namelist block adjust sucessfully read.'
    205 CONTINUE
  END IF

  CALL mpupdatei(hydradj,1)
  CALL mpupdatei(wndadj,1)
  CALL mpupdatei(obropt,1)
  CALL mpupdater(obrzero,1)

  IF(myproc == 0) THEN
    READ(unum, adas_radaropt,  END=210)
    WRITE(6,*) 'Namelist block adas_radaropt sucessfully read.'
    210 CONTINUE
  END IF

  CALL mpupdatei(raduvobs,1)
  CALL mpupdatei(radrhobs,1)
  CALL mpupdatei(radistride,1)
  CALL mpupdatei(radkstride,1)
  CALL mpupdater(refrh,1)
  CALL mpupdater(rhradobs,1)
  CALL mpupdatei(radcldopt,1)
  CALL mpupdatei(radqvopt,1)
  CALL mpupdatei(radqcopt,1)
  CALL mpupdatei(radqropt,1)
  CALL mpupdatei(radptopt,1)
  CALL mpupdater(refsat,1)
  CALL mpupdater(rhrad,1)
  CALL mpupdater(refcld,1)
  CALL mpupdater(cldrad,1)
  CALL mpupdater(ceilopt,1)
  CALL mpupdater(ceilmin,1)
  CALL mpupdater(dzfill,1)
  CALL mpupdater(refrain,1)
  CALL mpupdater(radsetrat,1)
  CALL mpupdater(radreflim,1)
  CALL mpupdater(radptgain,1)

  IF(myproc == 0) THEN
    READ(unum, adas_cloud, END=220)
    WRITE(6,*) 'Namelist block adas_cloud sucessfully read.'
    220 CONTINUE
  END IF

  CALL mpupdatei(cloudopt,1)
  CALL mpupdatei(clddiag,1)
  CALL mpupdatei(cld_files,1)
  CALL mpupdater(range_cld,1)
  CALL mpupdater(refthr1,1)
  CALL mpupdater(refthr2,1)
  CALL mpupdater(hgtrefthr,1)
  CALL mpupdater(thresh_cvr,1)
  CALL mpupdatei(bgqcopt,1)
  CALL mpupdatei(cldqvopt,1)
  CALL mpupdater(rh_thr1,1)
  CALL mpupdater(cvr2rh_thr1,1)
  CALL mpupdater(rh_thr2,1)
  CALL mpupdater(cvr2rh_thr2,1)
  CALL mpupdatei(cldqcopt,1)
  CALL mpupdater(qvslimit_2_qc,1)
  CALL mpupdatei(cldqropt,1)
  CALL mpupdater(qrlimit,1)
  CALL mpupdater(frac_qr_2_qc,1)
  CALL mpupdatei(cldwopt,1)
  CALL mpupdater(wmhr_Cu,1)
  CALL mpupdater(wmhr_Sc,1)
  CALL mpupdater(wc_St,1)
  CALL mpupdatei(cldptopt,1)
  CALL mpupdater(frac_qw_2_pt,1)
  CALL mpupdater(frac_qc_2_lh,1)
  CALL mpupdater(max_lh_2_pt,1)
  CALL mpupdatei(smth_opt,1)
  CALL mpupdatei(nirfiles,1)

  CALL mpupdatec(ir_fname,256*mx_sat)
  CALL mpupdatec(ircalname,256*mx_sat)
  CALL mpupdatei(nvisfiles,1)
  CALL mpupdatec(vis_fname,256*mx_sat)
  CALL mpupdatec(viscalname,256*mx_sat)

  IF(myproc == 0) THEN
    READ(unum, adas_typ,    END=240)
    WRITE(6,*) 'Namelist block adas_typ sucessfully read.'
    240 CONTINUE
  END IF

  CALL mpupdatei(ianxtyp,mx_pass)

  IF(myproc == 0) THEN
    READ(unum, adas_range,  END=250)
    WRITE(6,*) 'Namelist block adas_range sucessfully read.'
    250 CONTINUE
  END IF

  CALL mpupdater(sfcqcrng,1)
  CALL mpupdater(xyrange,mx_pass)

  IF(myproc == 0) THEN
    READ(unum, adas_kpvar,  END=255)
    WRITE(6,*) 'Namelist block adas_kpvar sucessfully read.'
    255 CONTINUE
  END IF
  CALL mpupdater(kpvar, mx_nvar_anx)

  IF(myproc == 0) THEN
    READ(unum, adas_zrange, END=260)
    WRITE(6,*) 'Namelist block adas_zrange sucessfully read.'
    260 CONTINUE
  END IF
  CALL mpupdater(zrange,mx_pass)

  IF(myproc == 0) THEN
    READ(unum, adas_thrng,  END=270)
    WRITE(6,*) 'Namelist block adas_thrng sucessfully read.'
    270 CONTINUE
  END IF
  CALL mpupdater(thrng,mx_pass)

  IF(myproc == 0) THEN
    READ(unum, adas_trnrng,  END=280)
    WRITE(6,*) 'Namelist block adas_trnrng sucessfully read.'
    280 CONTINUE
  END IF
  CALL mpupdatei(trnropt,mx_pass)
  CALL mpupdater(trnrcst,mx_pass)
  CALL mpupdater(trnrng,mx_pass)


  IF(myproc == 0) THEN
    READ(unum, adas_backerf,  END=300)
    WRITE(6,*) 'Namelist block adas_backerf sucessfully read.'
    300 CONTINUE
  END IF
  CALL mpupdatec(backerrfil,256)

  IF(myproc == 0) THEN
    READ(unum, adas_sng,  END=310)
    WRITE(6,*) 'Namelist block adas_sng sucessfully read.'
    310 CONTINUE
  END IF
  CALL mpupdatec(srcsng,8*nsrc_sng)
  CALL mpupdatec(sngerrfil,256*nsrc_sng)
  CALL mpupdatei(iusesng,(nsrc_sng+1)*mx_pass)

  IF(myproc == 0) THEN
    READ(unum, adas_ua,  END=320)
    WRITE(6,*) 'Namelist block adas_ua sucessfully read.'
    320 CONTINUE
  END IF
  CALL mpupdatec(srcua,8*nsrc_sng)
  CALL mpupdatec(uaerrfil,256*mx_ua_file)
  CALL mpupdatei(iuseua,(nsrc_ua+1)*mx_pass)

  IF(myproc == 0) THEN
    READ(unum, adas_radar,  END=330)
    WRITE(6,*) 'Namelist block adas_radar sucessfully read.'
    330 CONTINUE
  END IF
  CALL mpupdatei(nradfil,1)
!
!    Turn off "raduvobs" if we're not processing any radar files, in case
!    the user accidentally left it enable.
!
!    This will prevent possible unpleasant side effects, and in the
!    very least wasted computations.
!
  IF (nradfil == 0 ) raduvobs = 0

  CALL mpupdatec(radfname,256*mx_rad)
  CALL mpupdatec(srcrad,8*nsrc_rad)
  CALL mpupdatec(raderrfil,256*nsrc_rad)
  CALL mpupdatei(iuserad,(nsrc_rad+1)*mx_pass)

  IF(myproc == 0) THEN
     READ(unum, adas_retrieval,  END=340)
     WRITE(6,*) 'Namelist block adas_retrieval sucessfully read.'
     340 CONTINUE
  END IF

  CALL mpupdatei(nretfil,1)
  CALL mpupdatec(retfname,256*mx_ret)
  CALL mpupdatec(srcret,8*nsrc_ret)
  CALL mpupdatec(reterrfil,256*mx_ret)
  CALL mpupdatei(iuseret,(nsrc_ret+1)*mx_pass)

!-----------------------------------------------------------------------
!  Compute squared input variables
!-----------------------------------------------------------------------

  DO i=1,mx_nvar_anx
    kpvrsq(i)=kpvar(i)*kpvar(i)
  END DO

  IF (unum /= 5) THEN
    CLOSE(unum)
    CALL retunit(unum)
  END IF

  RETURN
END SUBROUTINE initadas

!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE RADINF                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE radinf (nvar_anx, rlim )

!-----------------------------------------------------------------------
!
!  PURPOSE:
!  For MPI message passing efficiency, we need to know the maximum
!  possible radii of influence that any any grid point will have so we
!  know processors need to communicate.
!
!  For the non-MPI case, this subroutine has no meaning, and isn't used.
!
!  AUTHOR:
!  Steve Leyton, CAPS, August, 2005.
!
!  MODIFICATION HISTORY:
!  10/06/2006   Add cloud information into the computations.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid parameters
  INCLUDE 'phycst.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
  INCLUDE 'adas.inc'
  INCLUDE 'adassat.inc'
  INCLUDE 'adjust.inc'
!
  REAL :: rlim                        ! 1st pass radius of influence
  REAL :: rpass, rsq, rlimsq
  REAL :: rngsqi
  REAL :: max_xyrange
  INTEGER :: i, nvar_anx
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

! Identify the maximum xyrange for all passes, then use to compute
! the maximum radius of influence for the given data type

  max_xyrange=xyrange(1)
  DO i = 2,npass
    IF (xyrange(i) .gt. max_xyrange) THEN
      max_xyrange=xyrange(i)
    END IF
  END DO
  rpass=max_xyrange*max_xyrange
  rlimsq=0.
  DO i=1,nvar_anx
    rsq=(kpvrsq(i)*rpass)
    rngsqi =1./rsq
    rlimsq=AMAX1(rlimsq,rsq)
  END DO
  rlimsq=-rlimsq*ALOG(wlim)
  rlim=SQRT(rlimsq)

  IF (cloudopt > 0) THEN
!
!   "dx" and "dy" should always be the same, but just to be safe, make sure
!   we have the larger of the two.
!
    rsq = dx
    IF (dy > rsq) rsq = dy
    rsq = rsq * i_perimeter

    IF (rsq > rlim) rlim = rsq
  END IF

  IF (myproc == 0)                                                         &
  WRITE(6,'(/1x,a,f10.2,a)') 'Influence cutoff radius (all data): ',rlim,' m.'

  RETURN

END SUBROUTINE radinf

!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE RADINFV                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE radinfv (nvar_anx,rlim)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!  For MPI message passing efficiency, we need to know the maximum
!  possible radii of influence that any any grid point will have so we
!  know processors need to communicate.
!
!  For the non-MPI case, this subroutine has no meaning, and isn't used.
!
!  This is the radial velocity version, which doesn't worry about things
!  such as cloud soundings.  Although this could be merged into radinf,
!  I've chosen not to go this way.
!
!  AUTHOR:
!  Steve Leyton, CAPS, August, 2005.
!
!  MODIFICATION HISTORY:
!  10/06/2006   Add cloud information into the computations.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid parameters
  INCLUDE 'phycst.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
  INCLUDE 'adas.inc'
  INCLUDE 'adassat.inc'
  INCLUDE 'adjust.inc'
!
  REAL :: rlim                        ! 1st pass radius of influence
  REAL :: rpass, rsq, rlimsq
  REAL :: rngsqi
  REAL :: max_xyrange
  INTEGER :: i, k
  INTEGER :: nvar_anx
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  rlim = -999.9
  IF (raduvobs == 0) RETURN               ! Skip if we're not using the data

! Identify the maximum xyrange for all radial velocity, then use to compute
! the maximum radius of influence

  max_xyrange = rlim
  DO k = 1,npass
    DO i=1,nsrc_rad
      IF (iuserad(i,k) == 0) CYCLE
      IF (xyrange(k) .gt. max_xyrange) THEN
        max_xyrange=xyrange(k)
      END IF
    END DO
  END DO

  IF (max_xyrange < 0.0) THEN
    IF (myproc == 0) WRITE(6,'(/1x,a)')                                 &
      'Raduvobs enabled, all values of iuserad are zero.  Disabling raduvobs'
  END IF
  rpass=max_xyrange*max_xyrange
  rlimsq=0.
  DO i=1,nvar_anx
    rsq=(kpvrsq(i)*rpass)
    rngsqi =1./rsq
    rlimsq=AMAX1(rlimsq,rsq)
  END DO
  rlimsq=-rlimsq*ALOG(wlim)
  rlim=SQRT(rlimsq)

  IF (myproc == 0)                                                      &
    WRITE(6,'(1x,a,f10.2,a)') 'Velocity influence cutoff radius: ',rlim,' m.'

  RETURN

END SUBROUTINE radinfv
