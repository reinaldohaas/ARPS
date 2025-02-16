!
!                               Richard Carpenter
!                               Univ. of Oklahoma
!                                  August 1992
!
!
!                   ########################################
!                   ########################################
!                   ########################################
!                   ########                        ########
!                   ########          SKEWT         ########
!                   ########                        ########
!                   ########################################
!                   ########################################
!                   ########################################

PROGRAM skewt

!-----------------------------------------------------------------------
!
!  HISTORY:
! 09/13/93  Modified somewhat for parallel ARPSTools version.
! 12/08/93  Added option for plotting Tv with or without water ldng.
! 1995/01/24  -tvnowl now also plots Tv of environment.
! 1995/12/07  In Derive added code not to bomb when Td = 0.
! 1997/03/06  Ported to throttle (Sun) from rossby (Sun).
!      Change istnm (INTEGER) to stnm (CHAR).
! 1997/03/21  Hodograph now in upper left corner of main plot.
!      Using color table input file (hardwired).
! 1997/04/22  Heights are now MSL, not AGL.
! 1997/04/23  Add a 2nd color scheme - all lines of each snd are 1 color.
! 1997/05/22  Improve handling of multiple soundings.
! 1999/10/11  Merged in COMET-Tinker version ("skewti"). Last COMET mod
!             was 1999/01/31.
!             Changed default tv_wl to False. [RLC]
!             Converted to Fortran-90 using TO_F90 by Alan Miller
! 1999/10/14  Numerous changes, merged version checks out.
!              Color file may be specified by -colorfile or
!              setenv SKEWT_COLORFILE. [RLC
!
!  Skewed temperatures are indicated by _sk.
!
!  Theta-E lines calculated using Bolton (1980), eq. 43. The slightly less
!  accurate method, eq. 35, dIFfers by no more that 0.2 K. Overall, they seem
!  to by within 0.5 K of the DOD Skew-T chart.
!
!  UNITS: All quantities are SI units unless labeled otherwise (e.g.,
!  presmb, tempc). Mixing ratios are kg/kg.
!
!  m/s = mph / 2.237 = kts / 1.944
!  mph = kts / 1.151
!
!----------------------------------------------------------------------=

  IMPLICIT NONE
  INCLUDE 'thermo.consts'

  LOGICAL, PARAMETER :: thetalines=.true.,  thetaelines=.true.,         &
                        mixinglines=.true., do_legend=.true.
  INTEGER, PARAMETER :: nmax=1000, npmx=200, ntmx=40, nticmx=80, nsmax=8
  REAL,    PARAMETER :: km2hft=32.8084

!----------------------------------------------------------------------=

  LOGICAL ::  &
      plot_sfc_parcel, plot_cb_parcel, plot_frame_2,   plot_tv,         &
      has_wind       , orig          , gempak_fmt,                      &
      weightq        , cape_use_t    , cape_use_irrev, write_modsnd  ,  &
      has_pres       , has_zz        , arpstools = .false.,             &
      do_hodo = .false., hodo_denmwind = .false., helcontrs = .false.,  &
      verbose   = .false., modify_snd = .false., tv_wl      = .false.,  &
      wrt_indcs = .true.,  logfiles   = .false., do_indices = .false.,  &
      plot_wind = .false., print_info = .false.

  INTEGER :: i , j        , k        , n        ,                       &
      ifcb     , nxtic    , nptic    , npres    ,                       &
      isnd     , nztic    , nc=0     , kc       ,                       &
      np       , klcl     , jgray    ,lcolors(nsmax), jhunits  ,        &
      kp       , irev     , ltypes(nsmax), kk   , nargs    , jdefcolor, &
      ksc      , n_td_sk  , kmod     , nsnd=0   , kmix     , nitems   , &
      iyr      , imon     , iday     , ihr      , imin     , ifont    , &
      jtcolor  , jtdcolor , jparcolor, jtvcolor , jhodoclr ,            &
      jcolor(nsmax), jscheme,ii,time(6,4)

  INTEGER :: preciptype_me, preciptype_laps(nmax),meprectypeskewt
  REAL    :: tw

  REAL ::  &
    tmin=0.0 , tmax=0.0 , tinc=10.0, pmin=0.0 , pref=1050.0         ,  &
    skewfac  , tskew    , dtskew   , t1       , t2       , p1       ,  &
    p2       , thetae   , pinc=100.0,zcb      ,  &
    dz       , stnelev  , px       , py       , pres     , es       ,  &
    plogmin  , plogmax  , theta    , dp       , tc       ,qqsgkg(14),  &
    pmbcb    , tccb     , esmb     , preslcl  , templcl  , tlcl     ,  &
    tvlcl    , zlcl     , deg2rad  , pi       , xx       , yy       ,  &
    psfc     , px1      , px2      , py1      , py2      , temp     ,  &
    siz0     , qs       ,  &
    zmin     , zmax     , qmin     , qmax     , siz08    , pyoffset ,  &
    yykm1    , dlnp     , sfctf    , sfctdf   ,rad2deg   , pmax     ,  &
    pspacing , pydelta  , tvavg    , y1       , y2       ,  &
    sfctf0   , sfctdf0  , pmb_mod  , pxoffset , zlnb     , preslnb  ,  &
    capetvwl , cintvwl  ,capetvnowl, cintvnowl, capet    , cint     ,  &
    capetvwlc, cintvwlc,capetvnowlc,cintvnowlc, capetc   , cintc    ,  &
    pxhodo1  , pxhodo2  , pyhodo1  , pyhodo2  , x1       , x2

  REAL :: presmb_s(nmax), tempc_s(nmax), tdewc_s(nmax)  ,            &
    temp_c(nmax)  , preslog(nmax) , pp(npmx)     ,th_sk(npmx,ntmx),  &
    xtic(nticmx)  , ptic(nticmx)  , pticlab(nticmx),pplog(npmx)   ,  &
    the_sk(npmx,ntmx),tq(npmx,ntmx) , alwcpr(nmax)  , alwcpa(nmax),  &
    pres_c(nmax)  , zzkm(nmax)    ,  &
    plog_c(nmax)  , ztic(nticmx)  , zticlab(nticmx),zz_s(nmax)    ,  &
    rh_s(nmax)    ,  &
    theta_s(nmax) , alwc_cr(nmax) , tempp(nmax)   , plogp(nmax)   ,  &
    tvp(nmax)     , tvwlp(nmax)   , thete(nmax)   , thetq(nmax)   ,  &
    tv_c(nmax)    , tvwl_c(nmax)  ,  &
    tempc_sk(nmax), tdewc_sk(nmax), tvc_sk(nmax)  , tvcp_sk(nmax) ,  &
    tcp_sk(nmax)  ,tvcwlp_sk(nmax), tcc_sk(nmax)  ,tvcwlc_sk(nmax),  &
    tvcnowlc_sk(nmax),  &
    pres_s(nmax)  , temp_s(nmax)  , tdew_s(nmax)  , qv_s(nmax)    ,  &
    qs_s(nmax)    , tv_s(nmax)    , qvgkg_s(nmax) ,  &
    uu_s(nmax)    , vv_s(nmax)    , dir_s(nmax)   , spd_s(nmax)   ,  &
    alwcgkgpa(nmax),alwcgkgpr(nmax),  &
    tempc_sk_orig(nmax), tdewc_sk_orig(nmax), presmb_s_orig(nmax),  &
    preslog_orig(nmax), wbt(nmax)
!    Array wbt(nmax) added for precip type algorithms.  1/18/98 EMK

  REAL :: pres500,ht500,temp500,tdewp500,spd500,dir500,pres700,ht700,   &
    temp700,tdewp700,spd700,dir700,pres850,ht850,temp850,tdewp850,      &
    spd850,dir850,showalter,kindex,liftedindex,totaltotals,sweat,       &
    brn,brnshear,wetbulbzero,convtemp,zlfc,preslfc,mixrat,              &
    mixrattot,prestopbl,xmxrat,precipwat,pccl,slat,slon,                &
    lidstrength,dirmean,spmean,helicity

  REAL :: indx(6,27),lat(6),lon(6)

  CHARACTER (LEN=255) :: sndfile(nsmax), strings(nsmax),file1,file2, indoutfile
  CHARACTER (LEN=3)   :: stid
  CHARACTER (LEN=6)   :: stnm
  CHARACTER (LEN=2)   :: chr,cmin,cmoi,cdayi,chri,cmini
  CHARACTER (LEN=3)   :: cstni
  CHARACTER (LEN=132) :: string, string2
  CHARACTER (LEN=100) :: capestring, color_file
  CHARACTER (LEN=19)  :: nm(33)
  DATA nm /                                                             &
    'ht LCL (m AGL)     ', 'ht LFC (m AGL)     ', 'ht LNB (m AGL)     ',&
    'pres LCL (mb)      ', 'pres LFC (mb)      ', 'pres LNB (mb)      ',&
    'CAPE-Tv-w/loading  ', 'CAPE-Tv-w/o loading', 'CAPE-using T only  ',&
    'CIN-Tv-w/loading   ', 'CIN-Tv-w/o loading ', 'CIN-using T only   ',&
    'Showalter Index    ', 'K-Index            ', 'Lifted Index       ',&
    'Total-totals Index ', 'SWEAT Index        ', 'Bulk Rich Shear    ',&
    'Precip. Water (cm) ', 'Wet-bulb zero (m)  ', 'Convective Temp    ',&
    'BL mr (sfc- 50mb)  ', 'pres CCL (mb)      ', 'Lidstrength Index  ',&
    'Storm Direction    ', 'Storm Speed (m/s)  ', 'SR Helicity        ',&
    'Place-Time:        ', 'Date:              ', 'Latitude:          ',&
    'Longitude:         ', 'MEta Sfc Prcp Type ', 'LAPS Sfc Precip    '/

  CHARACTER (LEN=14) :: precindex_me(6),precindex_laps(6)
  CHARACTER (LEN=3) :: place(6)

  INTEGER  :: tstart, tstop, tstep
  INTEGER  :: istatus

  ! jtdcolor was 12, jt was 20, jpar was 22
  DATA  jdefcolor/01/, jgray/04/, jtcolor/39/, jtdcolor/33/,            &
       jparcolor/22/, jtvcolor/22/, jhodoclr/21/
  DATA jcolor /30, 01, 20, 33, 22, 21, 03, 05/ !best for jscheme=2
  !DATA jcolor /01, 30, 20, 33, 22, 21, 03, 05/ !best for jscheme=2
  !DATA jcolor /01, 02, 30, 33, 22, 21, 03, 05/ !best for jscheme=1
  ! jscheme: 1 = first snd is multi-color
  ! jscheme: 2 = all snd are single color. Do not plot sfc parcel for snd>2
  ! jscheme: 3 = all snd are single color. Plot sfc parcel etc for snd>2

  DATA jscheme /2/  !color scheme -- 2 for single color per snd
  ! Update the following for each version
  DATA color_file /'../data/arpsplt/skewt.pltcbar'/
  DATA plot_sfc_parcel, plot_cb_parcel, plot_frame_2, plot_tv           &
      /.false., .false., .false., .false./
  DATA weightq, cape_use_t, cape_use_irrev /.false., .false., .TRUE./
  DATA nargs /0/, siz0/0.02/, ifont /2/
  DATA pspacing/50.0/ ! normalized vert spacing of lvls;
                      ! larger number -> closer spacing
  DATA qqsgkg /0.05,  0.1,  0.2,  0.5,  1.0,  2.0, 5.0,                 &
               10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0/ !Mixing ratio ref lines

  INCLUDE 'thermo.stfunc'

  CHARACTER(LEN=256) :: outfilename
  INTEGER :: lenstr

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!----------------------------------------------------------------------=
!   Set up
!----------------------------------------------------------------------=

  pi = 4. * ATAN (1.)
  deg2rad = pi/180.
  rad2deg = 1. / deg2rad

!  Get input
  outfilename = 'gmeta'
!  PRINT *, 'Call ParseInput'
  CALL parseinput (                                                     &
          nsmax         , nsnd          , sndfile       , has_wind   ,  &
          print_info    , verbose       , logfiles      , do_indices ,  &
          color_file    , outfilename   ,                               &
          do_hodo       , hodo_denmwind , jhunits       , helcontrs  ,  &
          plot_cb_parcel,plot_sfc_parcel, plot_tv       ,               &
          modify_snd    , weightq       , write_modsnd  ,               &
          tv_wl         , cape_use_t    , cape_use_irrev,               &
          tmin          , tmax          , pmin          , pmax          )

!  Confirm options

  PRINT '(/,a)', ' === Confirm user options ==='
  PRINT *, 'Done ParseInput, color file is ',TRIM(color_file)
  PRINT *, 'output file is ',TRIM(outfilename)

  IF(plot_tv)         PRINT *, 'Plotting Tv'
  IF(plot_frame_2)    PRINT *, 'Plotting frame 2'
  IF(plot_cb_parcel)  PRINT *, 'Plotting parcel lifted from cloud base'
  IF(plot_sfc_parcel) PRINT *, 'Plotting parcel lifted from surface'
  IF(weightq)         PRINT *, 'Weighting mixing ratio when modifying sounding'
  IF(cape_use_t)      PRINT *, 'Using T instead of Tv when computing CAPE'
  IF(cape_use_irrev)  PRINT *, 'Using pseudoadiabat'

! string describing CAPE

  capestring = 'Parcel CAPE computed by integ. '
  capestring = 'CAPE: '
  string = capestring

  IF (cape_use_t) THEN
!    string2 = ' T (no virt. correction)'
    string2 = ' T'
  ELSE IF (tv_wl) THEN
!    string2 = ' Tv (incl. condensate loading)'
    string2 = ' Tv*'
  ELSE
!    string2 = ' Tv (no condensate loading)'
    string2 = ' Tv'
  END IF

  capestring = trim(string) // trim(string2)

  string = capestring
  IF (cape_use_irrev) THEN
!    string2 = ' psuedoadiabat.'
    string2 = ' psuedo'
  ELSE
!    string2 = ' reversible adiabat.'
    string2 = ' rev'
  END IF

  capestring = trim(string) // trim(string2)

  PRINT *, TRIM(capestring)
  PRINT *, '=== End of user options ==='
  PRINT *, ' '

  IF (nsnd == 1) jscheme = 1

  DO i=1,6,1
    place(i)='---'
  END DO

!----------------------------------------------------------------------=
!   LOOP THROUGH SOUNDINGS
!----------------------------------------------------------------------=

  DO isnd=1,nsnd

  !...Read the sounding

  !  ReadSound returns P, Z, Theta, Qv, U, V.
  !  P or Z is computed hydrostatically, IF needed.

    PRINT *, 'Calling ReadSound, file = ', trim(sndfile(isnd))
    CALL readsound (pres_s, zz_s, theta_s, qv_s, uu_s, vv_s,  &
        has_pres, has_zz, pmbcb, tccb, dz, stnelev, psfc,  &
        sndfile(isnd), nmax, n, ifcb, plot_wind, gempak_fmt, has_wind,  &
        stid,stnm,slat,slon, iyr,imon,iday,ihr,imin,istatus)

    IF (istatus /= 0) CYCLE

  !...Compute other quantities from P,Z,Th,Qv,U,V.

    CALL derive (n, pres_s, zz_s, theta_s, qv_s, uu_s, vv_s,  &
        temp_s, tempc_s, tdew_s, tdewc_s, presmb_s, spd_s, dir_s, qs_s, rh_s, tv_s)

  ! Graphics stuff -- first sounding only

    IF (isnd == 1) THEN

  !  Figure out coords based on pres at top of sounding

      IF (tmin == 0. .AND. tmax == 0.) THEN
        tmin = -20.
        tmax = 40.
        pmin = 300.
        pmax = pref
        IF (presmb_s(n) < pmin) THEN
          pmin = pmin - pinc
        END IF
        IF (presmb_s(n) < pmin) THEN
          pmin = pmin - pinc
          tmin = tmin - tinc
        END IF
        IF (pmin == 300.) tmin = -10.
        IF (presmb_s(1) > pmax) pmax = pref + 50.
      END IF

      plogmin = LOG(pmin)
      plogmax = LOG(pmax)
      dtskew = tmax - tmin
      IF (plot_wind) dtskew = dtskew + tinc
      IF (pmin >= 200.) dtskew = dtskew - tinc
      IF (pmin >= 300.) dtskew = dtskew - tinc
      skewfac   = dtskew/(plogmax - plogmin)

    END IF !isnd


  !----------------------------------------------------------------------=
  !    ANALYSIS
  !----------------------------------------------------------------------=

    PRINT *, 'Writing sounding info to 1.skewt, 2.skewt'

  !  Save original values; truncate sounding above 100 mb

    DO k=1,n
      presmb_s_orig(k) = presmb_s(k)
      tempc_sk_orig(k) = tempc_s(k)
      tdewc_sk_orig(k) = tdewc_s(k)


  !     print *, k,pmin,presmb_S(k)
      IF (presmb_s(k) <= pmin) THEN
        n = k
        PRINT *, 'Truncating sounding to n,p(n): ', n, presmb_s(n)
        GO TO 8100
      END IF
    END DO
    8100 CONTINUE

  !  ::::::::::::::::::::::::::
  !  :::::::: SOUNDING ::::::::
  !  ::::::::::::::::::::::::::

  !  Note: both alwcPA and alwcPR are no longer obtained (01/28/93)

    irev = 1
    IF (cape_use_irrev) irev = 0
    PRINT *, 'CALLing AnalyzSnd2, irev = ', irev
    CALL analyzsnd2 (n,pres_s,temp_s,tdew_s, irev,  &
        zz_s,qv_s,qs_s,rh_s,tv_s,theta_s, preslcl,templcl,zlcl,preslnb,zlnb,  &
        alwcpr,tempp,tvp,tvwlp,  &
        cintvwl,capetvwl, cintvnowl,capetvnowl, cint,capet, 0)

    DO k=1,n
      qvgkg_s(k)= qv_s(k) * 1000.
    END DO

  !----------------------------------------------------------------------=
  !  Modify sounding (first sounding only)
  !----------------------------------------------------------------------=

    IF (isnd == 1 .AND. modify_snd) THEN
      PRINT *
      CALL sfcupd (presmb_s,zz_s,tempc_s,tdewc_s,qvgkg_s,theta_s,  &
          uu_s,vv_s, sfctf,sfctdf,orig,n,nmax, kmod,pmb_mod,sfctf0,sfctdf0,weightq)
      PRINT *, 'Done modifying sounding. n,kmod: ', n, kmod
      PRINT '(1x,a,4f8.2)', 'sfctf,sfctdf,sfctf0,sfctdf0: ',  &
          sfctf,sfctdf,sfctf0,sfctdf0
      PRINT *

  !  ...debrief the modified values (obtain p,Z,Th,Qv,U,V)

      DO k=1,n
        pres_s(k) = presmb_s(k) * 100.
        qv_s(k) = qvgkg_s(k) * 1.e-3
      END DO

      CALL derive (n, pres_s, zz_s, theta_s, qv_s, uu_s, vv_s,  &
          temp_s, tempc_s, tdew_s, tdewc_s, presmb_s, spd_s, dir_s,  &
          qs_s, rh_s, tv_s)

  ! Compute heights hydrostatiCALLy
  !
      zz_s(1) = 0.

      DO k=2,kmod
        tvavg = .5 * (tv_s(k-1) + tv_s(k))
        dz = rd/grav*tvavg * LOG(pres_s(k-1)/pres_s(k))
        zz_s(k) = zz_s(k-1) + dz
      END DO

      PRINT *, 'Calling AnalyzSnd2 after modifying sounding. irev=',irev
      CALL analyzsnd2(n,pres_s,temp_s,tdew_s, irev,                      &
                      zz_s,qv_s,qs_s,rh_s,tv_s,theta_s,                  &
                      preslcl,templcl,zlcl,preslnb,zlnb,                 &
                      alwcpr,tempp,tvp,tvwlp,                            &
                      cintvwl,capetvwl,cintvnowl,capetvnowl,cint,capet,0)

    END IF  ! modify_snd

  !----------------------------------------------------------------------=


    tvlcl = ftvirtnowl(templcl, fmixrat(preslcl,fsvpres(templcl-tfrz)))


    DO k=1,n
      zzkm(k) = zz_s(k) * 1.e-3
    END DO

  !----------------------------------------------------------------------=
  ! WRITE SOUNDING INFO TO FILE
  !----------------------------------------------------------------------=

    PRINT *
    PRINT *, 'SOUNDING PARAMETERS'
    PRINT *

    IF (logfiles) THEN
      WRITE (file1,'(A,I1,A)') 'a-', isnd, '.skewt'
      PRINT *, 'Writing to ', trim(file1)
      OPEN (UNIT=1, FILE=file1)
      WRITE (1,'(2a)') '# sounding file: ', sndfile(isnd)
      WRITE (1,9903)
      WRITE (1,9909)
    END IF
    PRINT *
    9904 FORMAT (4(a,f8.2))
    PRINT 9904, 'CAPE: ', capetvwl, ' CIN: ', cintvwl,  &
        '   Max updraft:', SQRT(2.*capetvwl)
    IF (cape_use_t) PRINT *, 'Values using T, not Tv*:'
    PRINT 9904, 'CAPE: ', capet, ' CIN: ', cint,  &
        '   Max updraft:', SQRT(2.*capet)
    IF (verbose) THEN
      PRINT *
      PRINT 9903
      PRINT 9909
    END IF

    9903 FORMAT ('# Pres   Z    T[C]   Td[C]   Qv[g/kg]  Qs  ',  &
        ' RH   Tv[C]  Theta   U    V')
    9909 FORMAT ('#',77('-'))
    9901 FORMAT (f7.2,f7.0,12F7.2)
    9902 FORMAT (4F8.2,2X,8F8.2)

    DO k=1,n
      IF (logfiles) WRITE (1,9901)  &
          pres_s(k)*1.e-2,zz_s(k),temp_s(k)-tfrz,tdew_s(k)-tfrz,  &
          qv_s(k)*1000.,qs_s(k)*1000.,rh_s(k),tv_s(k)-tfrz,theta_s(k),  &
          uu_s(k), vv_s(k)
      IF (verbose) THEN
        PRINT *
        PRINT 9901, pres_s(k)*1.e-2,zz_s(k),temp_s(k)-tfrz, tdew_s(k)-tfrz,  &
            qv_s(k)*1000.,qs_s(k)*1000.,rh_s(k),tv_s(k)-tfrz,theta_s(k),  &
            uu_s(k), vv_s(k)
      END IF
    END DO
    IF (logfiles) CLOSE (1)

  ! WRITE to ARPS input format file

    IF (logfiles) THEN
      WRITE (file1,'(A,I1,A)') 'a-', isnd, '.skewt'
      CALL writearpssnd (file1, stnelev, psfc, n,  &
          pres_s, zz_s, theta_s, temp_s, qv_s, rh_s, tdew_s, uu_s, vv_s)
    END IF

    9950 FORMAT (i3,2X,f8.2,f7.0,2X,2F8.2,2X,2F8.2)
    IF (plot_wind .AND. verbose) THEN
      PRINT *
      DO k=1,n
        PRINT 9950, k,pres_s(k)*1.e-2,zz_s(k), dir_s(k),spd_s(k),uu_s(k),vv_s(k)
      END DO
      PRINT *
    END IF


  !  ::::::::::::::::::::::::::::::::::::::::::::
  !  :::::::: PARCEL LIFTED FROM SURFACE ::::::::
  !  ::::::::::::::::::::::::::::::::::::::::::::

    PRINT *
    PRINT *, 'PARAMETERS FOR PARCEL LIFTED FROM SFC '
    PRINT *
    IF (logfiles) THEN
      WRITE (file1,'(A,I1,A)') 'b-', isnd, '.skewt'
      PRINT *, 'Writing ', trim(file1)
      OPEN (UNIT=1, FILE=file1)
      WRITE (1,'(2a)') '# sounding file: ', sndfile(isnd)
      WRITE (1,9910) preslcl*1.e-2,templcl-tfrz, zlcl
      WRITE (1,9911)
      WRITE (1,9909)
    END IF
    PRINT *
    PRINT 9910, preslcl*1.e-2,templcl-tfrz, zlcl
    IF (verbose) THEN
      PRINT *
      PRINT 9911
      PRINT 9909
    END IF

    9910 FORMAT ('# Parcel lifted from surface: LCL P,t,z: ',2F7.2,f7.0)
    9911 FORMAT ('# Pres   Z     ALWC[g/kg]  TP[C]  dTvP[C]   dTv*P[C]',  &
        '   dT')

    DO k=1,n
      IF (pres_s(k) <= preslcl) THEN
        IF (verbose) THEN
          PRINT *
          PRINT 9902, presmb_s(k),zz_s(k),alwcpr(k)*1000., tempp(k)-tfrz,  &
              tvp(k)-tv_s(k),tvwlp(k)-tv_s(k), tempp(k)-temp_s(k)
        END IF
        IF (logfiles) WRITE (1,9902)  &
            presmb_s(k),zz_s(k),alwcpr(k)*1000.,tempp(k)-tfrz,  &
            tvp(k)-tv_s(k),tvwlp(k)-tv_s(k), tempp(k)-temp_s(k)
      END IF
    END DO
    IF (logfiles) CLOSE (1)


  !  :::::::::::::::::::::::::::::::::::::::::::::::::::::
  !  :::::::: PARCEL LIFTED FROM KNOWN CLOUD BASE ::::::::
  !  :::::::::::::::::::::::::::::::::::::::::::::::::::::

    IF (ifcb == 1) THEN
      PRINT *
      PRINT *, 'PARAMETERS FOR PARCEL LIFTED FROM CLOUD BASE '
      PRINT *

      IF (logfiles) THEN
        WRITE ( file1,'(A,I1,A)') 'c-', isnd, '.skewt'
        WRITE (file2,'(A,I1,A)') 'd-', isnd, '.skewt'
        PRINT *, 'Writing ', trim(file1), ', ', trim(file2)
        OPEN (UNIT=1, FILE='3.skewt')
        OPEN (UNIT=2, FILE='4.skewt')
        WRITE (1,'(2a)') '# sounding file: ', sndfile(isnd)
        WRITE (2,'(2a)') '# sounding file: ', sndfile(isnd)
        WRITE (1,9911)
        DO k=1,n
          WRITE (2,9902) pres_s(k)/100., zz_s(k),1.00
        END DO
      END IF

  !  point kC = 1 (cloud base)

      kc = 1
      pres_c(kc) = pmbcb * 100.
      temp_c(kc) = tccb + tfrz
      qs = fmixrat(pres_c(kc),fsvpres(tccb))
      tv_c(kc) = ftvirt(temp_c(kc),qs,0.)
      tvwl_c(kc)= ftvirt(temp_c(kc),qs,0.)
      capetvwlc = 0.
      capetvnowlc = 0.
      capetc = 0.


  !  points kC = 2:nC

  !  point kSC is the first sounding point below the cloud base.
  !  note that k=kSC will not be used.

      DO k=1,n
        IF (presmb_s(k) < pmbcb) THEN
          kc = kc + 1
          pres_c(kc)= pres_s(k)
        END IF
      END DO
      nc = kc
      ksc = n - nc + 1
      PRINT *,  'nC,kSC: ', nc,ksc

      PRINT 9902, pres_c(1)/100., zcb,alwc_cr(1)*1000.,temp_c(1)-tfrz
      PRINT *, 'CALLing AnalyzSnd2 (cloud base), irev = ', irev
      CALL analyzsnd2 (nc,pres_s(ksc),temp_s(ksc),tdew_s(ksc), irev,  &
          zz_s(ksc),qv_s(ksc),qs_s(ksc),rh_s(ksc),tv_s(ksc),theta_s(ksc),  &
          pres_c,temp_c,zcb,xx,xx, alwc_cr,temp_c,tv_c,tvwl_c,  &
          cintvwlc,capetvwlc, cintvnowlc,capetvnowlc, cint,capet, 1)


      IF (logfiles) THEN
        WRITE (1,9922) pmbcb, tccb, zcb
        WRITE (1,9923)
        WRITE (1,9909)
      END IF
      PRINT *
      PRINT 9922, pmbcb, tccb, zcb
      IF (verbose) THEN
        PRINT *
        PRINT 9923
        PRINT 9909
      END IF

      9922 FORMAT ('# Cloud base p,T,z: ',2F7.2,f9.0)
      9923 FORMAT ('# Pres   Z     ALWC[g/kg]  TP[C]')

      PRINT 9904, 'CAPE: ', capetvwlc, ' CIN: ', cintvwlc

      kc = 1
      IF (verbose) THEN
        PRINT *
        PRINT 9902, pres_c(kc)/100., -99.,alwc_cr(kc)*1000.,temp_c(kc)-tfrz,  &
            tv_c(kc),tvwl_c(kc), temp_c(kc)
      END IF
      IF (logfiles) WRITE (1,9902) pres_c(kc)/100.,  &
          -99.,alwc_cr(kc)*1000.,temp_c(kc)-tfrz, tv_c(kc),tvwl_c(kc), temp_c(kc)

      DO kc=2,nc
        k  = kc + ksc - 1 !k=kC will be bogus.
        IF (verbose) THEN
          PRINT 9902, pres_c(kc)/100.,  &
              zz_s(k),alwc_cr(kc)*1000.,temp_c(kc)-tfrz,  &
              tv_c(kc)-tv_s(k),tvwl_c(kc)-tv_s(k), temp_c(kc)-temp_s(k)
        END IF
        IF (logfiles) WRITE (1,9902) pres_c(kc)/100.,  &
            zz_s(k),alwc_cr(kc)*1000.,temp_c(kc)-tfrz,  &
            tv_c(kc)-tv_s(k),tvwl_c(kc)-tv_s(k), temp_c(kc)-temp_s(k)
        IF (logfiles) WRITE (2,9902) pres_c(kc)/100.,  &
            zz_s(k),alwc_cr(kc)*1000.
      END DO
      IF (logfiles) THEN
        CLOSE (1)
        CLOSE (2)
      END IF
    END IF

  !  Write the modified sounding in model format.
  !  Remove point k=kmod, which had been added to the sounding.

    IF (isnd == 1 .AND. write_modsnd .AND. modify_snd) THEN
      file1 = trim(sndfile(isnd)) // '.mod'
      PRINT *, 'Writing modified sounding to: ', file1
      OPEN (UNIT=2, FILE=file1)
      WRITE (2,*) 'Modified sounding:', trim(file1)
      WRITE (2,8801) (n-1)/2, 2.*dz, stnelev, psfc, ifcb, pmbcb, tccb
      WRITE (2,8802) (theta_s(k),k=1,kmod-1), (theta_s(k),k=kmod+1,n)
      WRITE (2,8802) (qv_s(k),k=1,kmod-1), (qv_s(k),k=kmod+1,n)
      WRITE (2,8802) (uu_s(k),k=1,kmod-1), (uu_s(k),k=kmod+1,n)
      WRITE (2,8802) (vv_s(k),k=1,kmod-1), (vv_s(k),k=kmod+1,n)
      8801   FORMAT (//,i4, 3F9.2, i8,2F7.2)
      8802   FORMAT (/(1P,5G15.6))
      CLOSE (UNIT=2)
    END IF

  !----------------------------------------------------------------------=

    IF (isnd == 1) THEN

  !  Height axis

      kk = 1
      yykm1 = 10000.
  !dlnp = (LOG(pref)-LOG(pmin))/70. !spacing between labels
  !dlnp = (LOG(pref)-LOG(pmin))/30. !RLC 1995/12/07
      dlnp = (LOG(pref)-LOG(pmin))/pspacing

      DO k=1,n
        yy = LOG (presmb_s(k))
        IF (yykm1-yy > dlnp) THEN
          ztic(kk)= yy
          zticlab(kk)= zz_s(k) * 1.0E-3
          IF (jhunits == 1) zticlab(kk)= zticlab(kk) * km2hft
          yykm1 = yy
          kk = kk + 1
        END IF
      END DO
      nztic = kk - 1

    END IF !isnd

  !----------------------------------------------------------------------=
  !   SKEW THE TEMPERATURES
  !----------------------------------------------------------------------=

  !  Sounding

    n_td_sk = n

    DO k=1,n
      preslog(k)= LOG (presmb_s(k))
      preslog_orig(k) = LOG (presmb_s_orig(k))
      tskew = skewfac * (plogmax - preslog(k))
      tvc_sk(k) = ftvirtnowl(temp_s(k),qv_s(k)) + tskew - tfrz
      tempc_sk(k) = tempc_s(k) + tskew
      tdewc_sk(k) = tdewc_s(k) + tskew
      IF (n_td_sk == n .AND. tdewc_s(k) < -90.) n_td_sk = k - 1

      tskew = skewfac * (plogmax - preslog_orig(k))
      tempc_sk_orig(k) = tempc_sk_orig(k) + tskew
      tdewc_sk_orig(k) = tdewc_sk_orig(k) + tskew
    END DO

  !  Parcel
  !  Insert an extra point at the LCL

    np = n + 1

  !  ...below LCL

    klcl = 0
    DO kp=1,np
      k  = kp
      IF (klcl == 0 .AND. preslcl > pres_s(k)) THEN
        klcl = k
        GO TO 8800
      END IF
      IF (klcl > 0) k = kp - 1
      plogp(kp) = preslog(k)
      tskew = skewfac * (plogmax - plogp(kp))
      tvcwlp_sk(kp) = tvwlp(k) - tfrz + tskew
      tvcp_sk(kp) = tvp(k) - tfrz + tskew
      tcp_sk(kp)= tempp(k) - tfrz + tskew
      8800 CONTINUE
    END DO
    IF (klcl == 0) klcl = np ! EMK (1/31/98)

  !  ...at LCL

    kp = klcl
    PRINT *, 'kLCL', klcl
    plogp(kp) = LOG (preslcl*1.e-2)
    tskew = skewfac * (plogmax - plogp(kp))
    tvcwlp_sk(kp)= tvlcl - tfrz + tskew
    tvcp_sk(kp)= tvlcl - tfrz + tskew
    tcp_sk(kp)= templcl - tfrz + tskew

  !  Known Cloud Base
  !  the option tvnowl, CB, is not available  RLC 1994/03/04

    DO kc=1,nc
      plog_c(kc) = LOG (pres_c(kc)*1.e-2)
      tskew = skewfac * (plogmax - plog_c(kc))
      tvcwlc_sk(kc) = tvwl_c(kc) - tfrz + tskew
  !tvcnowlC_sk(kC) = tvnowl_C(kC) - tfrz + tskew
      tvcnowlc_sk(kc) = tvwl_c(kc) - tfrz + tskew
      tcc_sk(kc) = temp_c(kc) - tfrz + tskew
    END DO

  !-----------------------------------------------------------------------
  ! EXTENDED CALCULATION OF STABILITY PARAMETERS
  !-----------------------------------------------------------------------

    IF (do_indices) THEN

  !----------------------------------------------------------------------=
  !     INTERPOLATE THE VARIABLES TO 850, 700 and 500 mbs. JJM
  !----------------------------------------------------------------------=

      CALL interpress (n,pres_s,tempc_s,tdewc_s,spd_s,zz_s,  &
          dir_s,pres850,ht850,temp850,tdewp850,spd850,dir850,  &
          pres700,ht700,temp700,tdewp700,spd700,dir700,pres500,  &
          ht500,temp500,tdewp500,spd500,dir500,nmax)

  !----------------------------------------------------------------------=
  !     COMPUTE AN AVERAGE MIXING RATIO IN THE LOWEST 50 mb FOR USE AS
  !     A MIXED SFC PARCEL
  !----------------------------------------------------------------------=

      mixrattot=0.
      prestopbl=0.
      ii=0

      DO i=1,n,1
        IF (pres_s(i) > (pres_s(1)-5000.0)) THEN
          mixrat=xmxrat((pres_s(i)/100.),(tdewc_s(i)+273.16))
          mixrattot=mixrattot+mixrat
          ii=ii+1
        END IF
      END DO
      IF (ii > 1) THEN
        mixrat=mixrattot/(ii)
        prestopbl=pres_s(ii)
      ELSE
        mixrat=mixrattot
        prestopbl = pres_s(1) ! EMK 1/31/98
  !        prestopBL=pres_S(1)+1000. ! Original.  Why set this to be below
  !                                  ! the lowest sounding level?
      END IF


  !----------------------------------------------------------------------=
  !     CALCULATE THE PERTINENT INDICES (JJM)
  !----------------------------------------------------------------------=

      CALL indices(pres850,ht850,temp850,tdewp850,spd850,dir850,  &
          pres700,ht700,temp700,tdewp700,spd700,dir700,pres500,ht500,  &
          temp500,tdewp500,spd500,dir500,showalter,kindex,liftedindex,  &
          totaltotals,sweat,preslcl,zlcl,templcl,mixrat,prestopbl,  &
          pres_s,tempc_s,zz_s,tdewc_s,n,nmax,convtemp,precipwat,  &
          wetbulbzero,pccl,zlfc,preslfc,zlnb,uu_s,vv_s,brnshear,  &
          lidstrength,dirmean,spmean,helicity)

      IF ( pres_s(1) < 85000. ) THEN ! cms

        showalter   = 1000000000000000.
        kindex      = 1000000000000000.
        totaltotals = 1000000000000000.
        sweat       = 1000000000000000.

      END IF

  !     Determine the wet bulb temperature profile (1/18/98) EMK

      k = 1
      DO WHILE (k <= n)
        wbt(k) = tw((tempc_s(k)),(tdewc_s(k)),(pres_s(k)/100.))
        wbt(k) = wbt(k) + 273.15
  !       print*,'k = ',k
  !       print*,'Wet Bulb Temperature (K) = ',wbt(k)
        k = k + 1
      END DO

  !     Call precip. type algorithms (10/24/97, 1/18/98) EMK

      PRINT*, 'Calling meprectypeskewt...'

      preciptype_me = meprectypeskewt(nmax,n,pres_s,zz_s,tempc_s, tdewc_s,wbt)

      PRINT*,'Finished with meprectypeskewt'
  !      print*,'MesoEta Precip. Type = ',preciptype_me

      PRINT*, 'Calling prectype_laps_skewt...'

      CALL prectype_laps_skewt(nmax,n,tempc_s,pres_s,wbt, preciptype_laps)

      PRINT*,'Finished with prectype_laps_skewt'
  !      print*,'LAPS surface Precip. Type = ',preciptype_laps(1)


  !----------------------------------------------------------------------=
  !     FORMATTED PRINTOUTS (USED BY POST-PROCESSING PROGRAMS)
  !----------------------------------------------------------------------=

      indx(isnd,1)=zlcl
      2001 FORMAT(a19,4X,f6.0,4X,f6.0,4X,f6.0,4X,f6.0,4X,f6.0,4X,f6.0)
      indx(isnd,2)=zlfc
      2002 FORMAT(a19,4X,f6.1,4X,f6.1,4X,f6.1,4X,f6.1,4X,f6.1,4X,f6.1)
      indx(isnd,3)=zlnb
      indx(isnd,4)=preslcl/100.0
      indx(isnd,5)=preslfc/100.0
      indx(isnd,6)=preslnb/100.0
      indx(isnd,7)=capetvwl
      2003 FORMAT(a19,5X,f5.0,5X,f5.0,5X,f5.0,5X,f5.0,5X,f5.0,5X,f5.0)
      indx(isnd,8)=capetvnowl
      indx(isnd,9)=capet
      indx(isnd,10)=cintvwl
      indx(isnd,11)=cintvnowl
      indx(isnd,12)=cint
      indx(isnd,13)=showalter
      2004 FORMAT(a19,5X,f5.1,5X,f5.1,5X,f5.1,5X,f5.1,5X,f5.1,5X,f5.1)
      indx(isnd,14)=kindex
      indx(isnd,15)=liftedindex
      indx(isnd,16)=totaltotals
      indx(isnd,17)=sweat
      indx(isnd,18)=brnshear
      indx(isnd,19)=precipwat
      indx(isnd,20)=wetbulbzero
      indx(isnd,21)=convtemp
      indx(isnd,22)=mixrat*1000
      indx(isnd,23)=pccl
      indx(isnd,24)=lidstrength
      indx(isnd,25)=dirmean
      indx(isnd,26)=spmean
      indx(isnd,27)=helicity

  !     Adding precindex_me and precindex_laps for precip. type
  !     algorithms. (EMK 10/23/1997, 1/18/1998)

      IF (preciptype_me == 1) THEN
        precindex_me(isnd) = 'SN'
      ELSE IF (preciptype_me == 2) THEN
        precindex_me(isnd) = 'IP'
      ELSE IF (preciptype_me == 3) THEN
        precindex_me(isnd) = 'FZ RN'
      ELSE
        precindex_me(isnd) = 'RN'
      END IF

      IF (preciptype_laps(1) == 1) THEN
        precindex_laps(isnd) = 'SN'
      ELSE IF (preciptype_laps(1) == 2) THEN
        precindex_laps(isnd) = 'IP'
      ELSE IF (preciptype_laps(1) == 3) THEN
        precindex_laps(isnd) = 'FZ RN'
      ELSE
        precindex_laps(isnd) = 'RN'
      END IF

      place(isnd)=stid
      time(isnd,1)=iday
      time(isnd,2)=ihr
      time(isnd,3)=imin
      time(isnd,4)=imon
      lat(isnd)=slat
      lon(isnd)=slon

      2005 FORMAT(a19,a4,1X,i2.2,':',i2.2,1X,a3,1X,i2.2,':',i2.2,1X,a3,1X,  &
          i2.2,':',i2.2,1X,a3,1X,i2.2,':',i2.2,1X,a3,1X,i2.2,':',i2.2,  &
          1X,a3,1X,i2.2,':',i2.2)
      2006 FORMAT(a19)
      2007 FORMAT(a19,2X,i2.2,'/',i2.2,5X,i2.2,'/',i2.2,5X,i2.2,'/',i2.2,  &
          5X,i2.2,'/',i2.2,5X,i2.2,'/',i2.2,5X,i2.2,'/',i2.2)
      2010 FORMAT(a19,3X,f7.2,3X,f7.2,3X,f7.2,3X,f7.2,3X,f7.2,3X,f7.2)
      2011 FORMAT(a,8X,a5,5X,a5,5X,a5,5X,a5,5X,a5,5X,a5)

      IF (isnd == nsnd) THEN
        PRINT *
        PRINT *
        PRINT 2007, nm(29), (time(i,4),time(i,1), i=1,6)
        PRINT 2005, nm(28), (place(i),time(i,2),time(i,3), i=1,6)
        PRINT 2010, nm(30), (lat(i), i=1,6)
        PRINT 2010, nm(31), (lon(i), i=1,6)
        PRINT *

        DO j=1,3
          PRINT 2001, nm(j), (indx(i,j), i=1,6)
        END DO
        DO j=4,6
          PRINT 2002, nm(j), (indx(i,j), i=1,6)
        END DO
        DO j=7,9
          PRINT 2003, nm(j), (indx(i,j), i=1,6)
        END DO
        DO j=10,12
          PRINT 2002, nm(j), (indx(i,j), i=1,6)
        END DO
        DO j=13,16
          PRINT 2004, nm(j), (indx(i,j), i=1,6)
        END DO

        PRINT 2002, nm(17), (indx(i,17), i=1,6)
        PRINT 2010, nm(18), (indx(i,18), i=1,6)
        PRINT 2004, nm(19), (indx(i,19), i=1,6)
        PRINT 2001, nm(20), (indx(i,20), i=1,6)
        PRINT 2004, nm(21), (indx(i,21), i=1,6)
        PRINT 2010, nm(22), (indx(i,22), i=1,6)
        PRINT 2002, nm(23), (indx(i,23), i=1,6)
        PRINT 2010, nm(24), (indx(i,24), i=1,6)
        PRINT 2003, nm(25), (indx(i,25), i=1,6)
        PRINT 2010, nm(26), (indx(i,26), i=1,6)
        PRINT 2010, nm(27), (indx(i,27), i=1,6)
        PRINT 2011, nm(32), (precindex_me(i), i=1,6)
        PRINT 2011, nm(33), (precindex_laps(i), i=1,6)
      END IF

      IF (wrt_indcs) THEN

        indoutfile=stid//'_ind'//'.txt'
        OPEN (UNIT=7,FILE=TRIM(indoutfile),STATUS='unknown')

        IF (nsnd == 1) THEN
          WRITE (7,2057) nm(29), time(1,4), time(1,1), iyr
          WRITE (7,2055) nm(28), place(1), time(1,2), time(1,3)
          WRITE (7,2058) nm(30), lat(1)
          WRITE (7,2060) nm(31), lon(1)
          WRITE (7,*)

          DO j=1,27
            WRITE (7,2051) nm(j), indx(1,j)
          END DO
          WRITE (7,2061) nm(32), precindex_me(1)
          WRITE (7,2061) nm(33), precindex_laps(1)

          2051   FORMAT(a19,4X,f6.0)
          2052   FORMAT(a19,4X,f6.1)
          2053   FORMAT(a19,5X,f5.0)
          2054   FORMAT(a19,5X,f5.1)
          2055   FORMAT(a19,a4,1X,i2.2,':',i2.2,'Z')
          2057   FORMAT(a19,2X,i2.2,'/',i2.2,'/',i2.2)
          2058   FORMAT(a19,5X,f5.2)
          2059   FORMAT(a19,4X,f6.2)
          2060   FORMAT(a19,3X,f7.2)
          2061   FORMAT(a19,8X,a5)

        ELSE IF (isnd == nsnd) THEN

          WRITE (7,2007) nm(29), (time(i,4),time(i,1), i=1,6)
          WRITE (7,2005) nm(28), (place(i),time(i,2),time(i,3), i=1,6)
          WRITE (7,2010) nm(30), (lat(i), i=1,6)
          WRITE (7,2010) nm(31), (lon(i), i=1,6)
          WRITE (7,*)

          DO j=1,3
            WRITE (7,2001) nm(j), (indx(i,j), i=1,6)
          END DO
          DO j=4,6
            WRITE (7,2002) nm(j), (indx(i,j), i=1,6)
          END DO
          DO j=7,9
            WRITE (7,2003) nm(j), (indx(i,j), i=1,6)
          END DO
          DO j=10,12
            WRITE (7,2002) nm(j), (indx(i,j), i=1,6)
          END DO
          DO j=13,16
            WRITE (7,2004) nm(j), (indx(i,j), i=1,6)
          END DO

          WRITE (7,2002) nm(17), (indx(i,17), i=1,6)
          WRITE (7,2010) nm(18), (indx(i,18), i=1,6)
          WRITE (7,2004) nm(19), (indx(i,19), i=1,6)
          WRITE (7,2001) nm(20), (indx(i,20), i=1,6)
          WRITE (7,2004) nm(21), (indx(i,21), i=1,6)
          WRITE (7,2010) nm(22), (indx(i,22), i=1,6)
          WRITE (7,2002) nm(23), (indx(i,23), i=1,6)
          WRITE (7,2010) nm(24), (indx(i,24), i=1,6)
          WRITE (7,2003) nm(25), (indx(i,25), i=1,6)
          WRITE (7,2010) nm(26), (indx(i,26), i=1,6)
          WRITE (7,2010) nm(27), (indx(i,27), i=1,6)
          WRITE (7,2011) nm(32), (precindex_me(i), i=1,6)
          WRITE (7,2011) nm(33), (precindex_laps(i), i=1,6)

  !----------------------------------------------------------------------=
  !                    END FORMATTED PRINTOUTS (JJM)
  !----------------------------------------------------------------------=

          CLOSE(7)
        END IF
      END IF
    END IF ! do_indices

  !----------------------------------------------------------------------=
  !    PLOTTING
  !----------------------------------------------------------------------=

    PRINT *, 'Begin plotting'

    IF (isnd == 1) THEN

  !  X-tics

      nxtic = nint ((tmax-tmin)/tinc) + 1
      DO i=1,nxtic
        xtic(i) = (i-1)*tinc + tmin
      END DO


  !  P-tics

      nptic = INT((pmax-pmin)/pinc + 1.0)
      DO i=1,nptic
        pticlab(i)= pmin + (i-1)*pinc
        ptic(i) = LOG(pticlab(i))
      END DO


  !  Z-tics

      px1 = 0.1
      IF (do_indices) px1 = 0.065 ! EMK (1/18/98)
      px2 = 0.94
      py1 = 0.10
      py2 = 0.97
      IF (plot_wind) THEN
        px2 = px2 - 0.12
        IF (nsnd > 1) px2 = px2 - 0.05
      END IF
      PRINT *, px1, px2, py1, py2

      IF (jscheme > 1) THEN
        jtcolor = jcolor(isnd)
        jtvcolor = jcolor(isnd)
        jtdcolor = jcolor(isnd)
        jparcolor = jcolor(isnd)
      END IF

      lenstr = LEN_TRIM(outfilename)
      CALL xdevic_new(1,outfilename,lenstr,11)
      CALL xstctfn (color_file)
      CALL setcolors (-1)  ! Use -1 for file, or >=1 for predefined
  ! reverse black/white
!      CALL gscr (1, 0, 1.0, 1.0, 1.0)   ! uncomment them for NCARG
!      CALL gscr (1, 1, 0.0, 0.0, 0.0)
      IF (arpstools) ifont = 2
      CALL xcfont (ifont)
      CALL color (jdefcolor)
      CALL xpspac (px1,px2,py1,py2)
      CALL xmap   (tmin,tmax,plogmax,plogmin)

  ! Protected region for hodograph

      pxhodo1 = px1
      pxhodo2 = pxhodo1 + 0.3
      pyhodo2 = py2
      pyhodo1 = pyhodo2 - 0.3
      IF (do_hodo) THEN
        PRINT *, 'Protected region for hodograph'
        CALL plt2wrld (pxhodo1,pyhodo1,x1,y1)
        CALL plt2wrld (pxhodo2,pyhodo2,x2,y2)
  !print *, tmin,tmax,plogmax,plogmin
  !print *, x1,x2,y1,y2
        CALL xmask (x1,x2,y2,y1) !must reverse ys; see CALL XMAP
      END IF

      CALL xchmag (0.9*siz0)
      CALL xaxnmg (0.9*siz0)
      CALL xaxtik (0,0)
      CALL xaxfmt ('(I3)')
      CALL xxaxis (xtic, xtic, nxtic, plogmax)
  !     CALL XAXISX (tmin, plogmax, tinc)
      CALL xaxfmt ('(I4)')
      CALL xyaxis (tmin, ptic, pticlab, nptic) !pres

  ! Height axis along right edge

      CALL xaxnmg (0.65*siz0)
      CALL xaxnmg (0.90*siz0) !RLC 1995/12/07
      CALL xaxtik (1,-1)
      CALL xaxant (1,1)
      CALL xaxfmt ('(F4.1)')
      string = 'km'
      IF (jhunits == 1) THEN
        CALL xaxfmt ('(I3.3)')
        string = 'hft'
      END IF
      CALL xyaxis (tmax, ztic, zticlab, nztic) !hgt
      CALL wrld2plt (tmax,ztic(1),px,py)
      CALL plt2wrld (px+0.035,py-0.03,xx,yy) ! Original
      CALL xcharc (xx,yy,trim(string))
      CALL xbordr

      pyoffset = 0.03
      pxoffset = 0.01
      pydelta = 0.02

  !  Clip lines outside window

      CALL xwindw (tmin,tmax,plogmin,plogmax)

  !  P-lines

      CALL color (jgray)
      DO i=2,nptic
        CALL myline (tmin,ptic(i),tmax,ptic(i))
      END DO

  !  T-lines

#ifdef INTLOOP
      tstart = NINT(tmin-4*tinc)
      tstop  = NINT(tmax-tinc)
      tstep  = NINT(tinc)

      DO i=tstart,tstop,tstep
#else
      DO i=tmin-4.*tinc,tmax-tinc,tinc
#endif

        t1 = FLOAT(i)
        t2 = FLOAT(i) + dtskew
        p1 = plogmax
        p2 = plogmin
  !     IF (t1.LT.tmin) THEN
  !       p1 = p1 + (tmin-t1)/(t2-t1) * (p2 - p1)
  !       t1 = tmin
  !     END IF
  !     IF (t2.GT.tmax) THEN
  !       p2 = p1 + (tmax-t1)/(t2-t1) * (p2 - p1)
  !       t2 = tmax
  !     END IF
        CALL myline (t1,p1,t2,p2)
      END DO

  !  Theta-lines
  !  Note that th,the are the (skewed) values of T[C] corresp to theta,thetae.

      PRINT *, '... Plotting theta lines'
  !     CALL XBROKN (1,6,5,6)
      dp = 25.
      npres = nint ((pref-pmin)/dp) + 1
  !     PRINT *, 'npres ', npres
      DO k=1,npres
        pp(k) = pmin + (k-1) * dp
        pplog(k) = LOG (pp(k))
      END DO

      IF (thetalines) THEN
        DO i=1,nxtic+4
          theta = tfrz + tmin + i*tinc
          DO k=npres,1,-1
            tskew = skewfac * (plogmax - pplog(k))
            tc = theta * (pp(k)/p00mb) ** rcp - tfrz
            th_sk(k,i) = tc + tskew
          END DO
          CALL xcurve (th_sk(1,i), pplog, npres, 0)
        END DO
      END IF


  !  ThetaE-lines
  !  Procedure: given theta(1000mb), set a reference Theta-E. THEN loop through
  !  other pres levels, changing T,Theta, until we find the T at which Th-E = the
  !  reference value.
  !  See Bolton (1980).


      IF (thetaelines) THEN
        CALL xbrokn (1,6,1,6)
        PRINT *, '... Plotting theta-e lines'
  !     DO i=1,2*nxtic-3
        DO i=1,nxtic-1

  !  lowest point (reference)

          k = npres
          pres = pp(k) * 100.
  !       theta = tfrz + tmin + i*tinc*.5
          theta = tfrz + tmin + i*tinc
          temp = theta * (pres/p00) ** rcp
          tc = temp - tfrz
          es = fsvpres(tc)
          qs = fmixrat(pres,es)
          thetae = fthetae(pres,temp,temp,qs)
  !       PRINT *, 'Theta-E line:', i, theta-tfrz, thetae-tfrz
          tskew = skewfac * (plogmax - pplog(k))
          the_sk(k,i)= tc + tskew

  !  other points (iterate)

          DO k=npres-1,1,-1
            pres = pp(k) * 100.
            CALL the2t (temp, thetae, pres)
            tskew = skewfac * (plogmax - pplog(k))
            the_sk(k,i) = temp - tfrz + tskew
          END DO
          CALL xcurve (the_sk(1,i), pplog, npres, 0)
        END DO
        !CALL Color (jdefcolor)
      END IF

!-----------------------------------------------------------------------
!
!  Mixing ratio lines
!
!-----------------------------------------------------------------------

      IF (mixinglines) THEN
        CALL xchmag (0.6*siz0)
        !CALL Color (3)
        PRINT *, '... Plotting mixing ratio lines'
        !CALL XDASH
        CALL xbrokn (6,16,6,16)

        DO i=1,14

          DO k=1,npres
            tskew = skewfac * (plogmax - pplog(k))
            esmb = pp(k) * qqsgkg(i) / (1000.*cp622 + qqsgkg(i))
            tc = 243.5 / (17.67/LOG(esmb/6.112) - 1.)
            tq(k,i) = tc + tskew
            IF ( pp(k) <= 490.0 ) kmix = k
          END DO

          !  kmix is the uppermost point to plot
          !kmix = MIN (kmix,npres/2)

          CALL xcurve (tq(kmix,i), pplog(kmix), npres-kmix+1, 0)
          IF (qqsgkg(i) >= 1.0) THEN
            WRITE (string,'(I2)') INT(qqsgkg(i))
          ELSE IF (qqsgkg(i) >= 0.1) THEN
            WRITE (string,'(F3.1)') qqsgkg(i)
          ELSE
            WRITE (string,'(F4.2)') qqsgkg(i)
          END IF
          CALL xcharc( tq(kmix,i),pplog(kmix),trim(string) )
        END DO
      END IF

      CALL xwdwof  !Turn off window clipping
      IF (do_hodo) CALL xunmsk (1) ! Un-protect hodograph region
    END IF !isnd=1

  ! Label the plot

    CALL color (jcolor(isnd))
    CALL plot_strings (sndfile(isnd), stid, do_indices,                 &
        cape_use_t, gempak_fmt, modify_snd, plot_cb_parcel,             &
        plot_sfc_parcel, print_info, tv_wl,                             &
        isnd , nsnd , iyr      , imon     , iday     , ihr      ,       &
        imin     , px1      , px2      , py1      , py2      ,          &
        capet    , capetc   ,capetvnowl,capetvnowlc,capetvwl , capetvwlc,  &
        cint     , cintc    , cintvnowl,cintvnowlc, cintvwl  , cintvwlc ,  &
        pmb_mod  , pmbcb    , preslcl  , preslnb  ,                     &
        sfctf    , sfctf0   , sfctdf   , sfctdf0  , zcb, zlcl, zlnb ,   &
        preslfc , zlfc , convtemp , precipwat , liftedindex ,           &
        totaltotals , kindex , sweat , wetbulbzero , brnshear ,         &
        lidstrength , dirmean , spmean , helicity , preciptype_me,      &
        preciptype_laps,nmax,n,jcolor(isnd))


  !  Linetypes: 0,1 full; 2 dash; 3 dot; 4 my dot; 5 my dash; 6 my dash-dot;
  !  negative bold.

!----------------------------------------------------------------------=
!
!   FINALLY, PLOT THE CURVES!
!
!----------------------------------------------------------------------=

    CALL color (jcolor(isnd))

  ! Hodograph

    IF (do_hodo) THEN
      IF (isnd <= 2) THEN
        CALL hodograph (pres_s, zz_s, theta_s, qv_s, uu_s, vv_s,        &
                        n,isnd, sndfile(isnd), helcontrs,hodo_denmwind, &
                        jdefcolor, jgray, jcolor(isnd),                 &
                        pxhodo1,pxhodo2,pyhodo1,pyhodo2)
      END IF

      CALL xpspac (px1,px2,py1,py2)
      CALL xmap   (tmin,tmax,plogmax,plogmin)

    END IF

!----------------------------------------------------------------------=

    CALL xwindw (tmin,tmax,plogmin,plogmax)

!-----------------------------------------------------------------------
!
!  Plot T, Td, Tparcel
!
!-----------------------------------------------------------------------

    CALL color (jcolor(isnd))
    CALL xthick(2) !thickness of lines, 1 or 2
    IF (isnd > 2) CALL xthick(1)
    ! T
    CALL xfull
    IF (isnd == 1 .AND. jscheme == 1) CALL color (jtcolor)
    CALL xcurve (tempc_sk, preslog, n, 0)
    ! Td
    CALL xdash
    IF (isnd == 1 .AND. jscheme == 1) CALL color (jtdcolor)
    CALL xcurve (tdewc_sk, preslog, n_td_sk, 0)
    ! Tv
    IF (plot_tv .AND. (isnd == 1 .OR. jscheme == 3) ) THEN
      IF (isnd == 1 .AND. jscheme == 1) CALL color (jtvcolor)
      CALL xbrokn (10,5,1,5)
      CALL xcurve (tvc_sk, preslog, n, 0)
    END IF

    !  Original part of modified sounding

    IF (modify_snd .AND. (isnd == 1 .OR. jscheme == 3) ) THEN
      CALL xthick (1)
      CALL xfull
      IF (isnd == 1 .AND. jscheme == 1) CALL color (jtcolor)
      CALL xcurve (tempc_sk_orig, preslog_orig, kmod, 0)

      CALL xdash
      IF (isnd == 1 .AND. jscheme == 1) CALL color (jtdcolor)
      CALL xcurve (tdewc_sk_orig, preslog_orig, kmod, 0)
      IF (isnd <= 2) CALL xthick (2)
    END IF


  !  Parcel lifted from surface

    IF (plot_sfc_parcel .AND. (isnd == 1 .OR. jscheme == 3) ) THEN
      IF (isnd == 1 .AND. jscheme == 1) CALL color (jparcolor)
      CALL lintyp (03)
      IF (plot_tv) THEN
        IF (tv_wl) THEN
          CALL xcurve (tvcwlp_sk, plogp, np, 0)
          CALL xscatter (tvcwlp_sk(klcl),plogp(klcl),1,4, 0.005)
        ELSE
          CALL xcurve (tvcp_sk, plogp, np, 0)
          CALL xscatter (tvcp_sk(klcl),plogp(klcl),1,4, 0.005)
        END IF
      ELSE
        CALL xcurve (tcp_sk, plogp, np, 0)
        CALL xscatter (tcp_sk(klcl),plogp(klcl),1,4, 0.005)
      END IF
    END IF


  !  Parcel lifted from cloud base

    IF (ifcb == 1 .AND. plot_cb_parcel .AND.  &
          (isnd == 1 .OR. jscheme == 3) ) THEN
      CALL xdot
      IF (isnd == 1 .AND. jscheme == 1) CALL color (jparcolor)
      IF (plot_tv) THEN
        IF (tv_wl) THEN
          CALL xcurve (tvcwlc_sk, plog_c, nc, 0)
          CALL xfull
          CALL xscatter (tvcwlc_sk(1),plog_c(1),1,4, 0.005)
        ELSE
          CALL xcurve (tvcnowlc_sk, plog_c, nc, 0)
          CALL xfull
          CALL xscatter (tvcnowlc_sk(1),plog_c(1),1,4, 0.005)
        END IF
      ELSE
        CALL xcurve (tcc_sk, plog_c, nc, 0)
        CALL xfull
        CALL xscatter (tcc_sk(1),plog_c(1),1,4, 0.005)
      END IF
    END IF

    CALL xwdwof !Turn off window clipping

  !...Plot wind

    CALL color (jcolor(isnd))
    CALL xfull
    CALL xthick(1)
    IF (isnd <= 2 .AND. plot_wind) THEN
      PRINT *, '...Plotting wind barbs'
      string = ' '
      IF (nsnd > 1) WRITE (string,'(I1)') isnd
      CALL plotwind (presmb_s, dir_s,spd_s, n, pspacing,  &
          jscheme, string, do_indices)
    END IF


  !...Plot LAPS Precip (EMK, 1/18/98)

    IF (do_indices) THEN
      CALL color (jcolor(isnd))
      CALL xfull
      CALL xthick(1)
      IF (isnd <= 2) THEN
        PRINT *, '...Plotting LAPS Precip. Type'
        string = ' '
        IF (nsnd > 1) WRITE (string,'(I1)') isnd
        CALL plotprecip (presmb_s, n, pspacing, jscheme, string,preciptype_laps)
      END IF
    END IF

  !...Plot legend

    IF (isnd == 1 .AND. do_legend) THEN

      nitems = 2
      strings(1)= 'T'
      strings(2)= 'Td'
      ltypes(1) = -1
      ltypes(2) = -2
      lcolors(1)= jtcolor
      lcolors(2)= jtdcolor

      IF (modify_snd) THEN
        nitems  = 4
        strings(1) = 'T (modified)'
        strings(2) = 'Td (modified)'
        strings(3) = 'T (original)'
        strings(4) = 'Td (original)'
        ltypes(1) = -1
        ltypes(2) = -2
        ltypes(3) = 1
        ltypes(4) = 2
        lcolors(1) = jtcolor
        lcolors(2) = jtdcolor
        lcolors(3) = jtcolor
        lcolors(4) = jtdcolor
      END IF

      IF (plot_tv) THEN
        nitems  = nitems + 1
        strings(nitems) = 'Tv'
        ltypes(nitems) = -6
        lcolors(nitems) = jtvcolor
      END IF

      IF (plot_sfc_parcel .OR. plot_cb_parcel) THEN
        nitems  = nitems + 1
        strings(nitems) = 'T (parcel)'
        IF (plot_tv) strings(nitems) = 'Tv* (parcel)'
        ltypes(nitems) = 03
        lcolors(nitems) = jparcolor
      END IF

      CALL legend (nitems, strings, ltypes,lcolors,  &
          px1+0.01, 0.5*(py1+py2), 0, 0.9*siz0)
    END IF

!    CALL sflush     ! uncomment this for NCARG

  END DO !isnd

!----------------------------------------------------------------------=
!   FRAME 2
!----------------------------------------------------------------------=

  IF (.NOT. plot_frame_2) GO TO 9000
  CALL xframe

  PRINT *, '... Plotting frame 2'

  siz08 = 0.8 * siz0
  CALL xchmag (0.9*siz0)
  CALL xaxnmg (0.9*siz0)

  !  plot ALWC

  DO k=1,n
    alwcgkgpa(k) = alwcpa(k) * 1000.
    alwcgkgpr(k) = alwcpr(k) * 1000.
    qvgkg_s(k) = qv_s(k) * 1000.
  END DO

  px1 = 0.10
  px2 = 0.95
  py2 = 0.90
  py1 = 0.50
  zmin = 0.
  zmax = 4.
  DO i=1,10
    IF (zzkm(n) > zmax) zmax = zmax + 1.
  END DO
  qmin = 0.
  qmax = 10.
  CALL pcurves (px1,px2,py1,py2, zmin,zmax,qmin,qmax, 0,  &
      2,'Height [km]','*',0.,   1,'Mixing ratio [g/kg]','(I2)',0.,  &
      n,3,zzkm, qvgkg_s,alwcgkgpa,alwcgkgpr,xx, 2,1,4,0, 2,1,4,0)

  DO k=1,n
    WRITE (2,*) zz_s(k), alwcpr(k)*1000.
  END DO

  strings(1)= 'Water vapor'
  strings(2)= 'LWC (pseudoadiab)'
  strings(3)= 'LWC (reversible)'
  ltypes(1) = 2
  ltypes(2) = 1
  ltypes(3) = 4
  CALL legend (3,strings,ltypes,ltypes, px1,py2-0.02,0,siz08)


  !  compute potential temperatures

  py2 = 0.40
  py1 = 0.10
  DO k=1,n
    tlcl = ftlcl(temp_s(k),tdew_s(k))
    thete(k) = fthetae(pres_s(k),temp_s(k),tlcl,qv_s(k))
    thetq(k) = fthetaq(pres_s(k),temp_s(k),qv_s(k), fthq_cwcptrm(qv_s(k)))
  END DO

  qmin = 320.
  qmax = 320.
  DO i=1,10
    DO k=1,n
      IF (thete(k) > qmax) qmax = qmax + 10.
      IF (thetq(k) < qmin) qmin = qmin - 10.
    END DO
  END DO
  CALL pcurves (px1,px2,py1,py2, zmin,zmax,qmin,qmax, 0,  &
      1,'Height [km]','*',0., 1,'Potential Temperature [K]','(I3)',0.,  &
      n,2,zzkm, thete,thetq,xx,xx, 1,4,0,0, 1,4,0,0)

  strings(1)= 'Theta-E'
  strings(2)= 'Theta-Q'
  ltypes(1) = 1
  ltypes(2) = 4
  CALL legend (2,strings,ltypes,ltypes, px1,py1+.10,0,siz08)

  string = 'File: ' // sndfile(isnd)
  CALL xchmag (0.6*siz0)
  CALL plt2wrld (px1,0.03,xx,yy)
  CALL xcharl (xx,yy, trim(string))
!  CALL sflush                ! uncommented this for NCARG

  9000 CONTINUE

!  CALL sflush                ! uncommented this for NCARG
  CALL xgrend
  PRINT *,  'SKEWT: Normal completion'

  STOP
END PROGRAM skewt
!
!                   ########################################
!                   ########################################
!                   ########################################
!                   ########                        ########
!                   ########        PLOTWIND        ########
!                   ########                        ########
!                   ########################################
!                   ########################################
!                   ########################################

SUBROUTINE plotwind (presmb, dir, spd, n, pspacing, jscheme,            &
                     label, do_indices)

  IMPLICIT NONE
  INTEGER, INTENT(IN)                      :: n
  REAL, INTENT(IN OUT)                     :: presmb(n)
  REAL, INTENT(IN OUT)                     :: dir(n)
  REAL, INTENT(IN OUT)                     :: spd(n)
  REAL, INTENT(IN)                         :: pspacing
  INTEGER, INTENT(IN OUT)                  :: jscheme
  CHARACTER (LEN=*), INTENT(IN)            :: label
  LOGICAL, INTENT(IN OUT)                  :: do_indices

  INTEGER :: k

  REAL :: a,xx,yy,yykm1,xmin,xmax,ymin,ymax, dlnp, px1,px2,py1,py2

  REAL, SAVE :: px=0.0

  !  spacing between printed levels

  CALL xqmap (xmin,xmax,ymin,ymax)
  dlnp = ABS(ymax-ymin)/pspacing
  yykm1 = 10000.

  !  location of wind barbs

  IF (px == 0.0) THEN
    CALL xqpspc (px1,px2,py1,py2)
    a = 0.11
    IF (do_indices) a = 0.12
    px = px2 + a
    CALL plt2wrld (px,0.0,xx,yy)
  ELSE
    a = 0.05
    IF (do_indices) a = 0.08
    px = px + a
    CALL plt2wrld (px,0.0,xx,yy)
  END IF

  DO k=1,n
    IF (ABS(dir(k)) > 360.0 .OR. ABS(spd(k)) > 250.0) GO TO 99
    yy = LOG (presmb(k))
    IF ( yykm1-yy > dlnp ) THEN
      yykm1 = yy
  ! IF (jscheme.NE.2) THEN
  !        IF (spd(k).GT.117.5/2.) THEN
  !   CALL Color (04)
  !        Else IF (spd(k).GT.77.5/2.) THEN
  !   CALL Color (24)
  !        Else IF (spd(k).GT.37.5/2.) THEN
  !   CALL Color (23)
  !        Else
  !   CALL Color (01)
  ! END IF
  ! END IF
      CALL barb (xx,yy,dir(k),spd(k), 0,0.035)
    END IF
    99   CONTINUE
  END DO

  CALL plt2wrld (px,py1-0.02,xx,yy)
  CALL xcharc (xx,yy,trim(label))


  IF (jscheme /= 2) CALL color (01)
END SUBROUTINE plotwind




!                   ########################################
!                   ########################################
!                   ########################################
!                   ########                        ########
!                   ########         DERIVE         ########
!                   ########                        ########
!                   ########################################
!                   ########################################
!                   ########################################


SUBROUTINE derive (n, pres, zz, theta, qv, uu, vv,  &
                   temp, tempc, tdew, tdewc, presmb, spd, dir, qs, rh, tv)

  IMPLICIT NONE
  INTEGER, INTENT(IN)                      :: n
  REAL, INTENT(IN), DIMENSION(n) :: zz, pres, theta, qv, uu, vv
  REAL, INTENT(OUT), DIMENSION(n) :: &
      temp, tempc, tdew, tdewc, presmb, spd, dir, qs, rh, tv
  INCLUDE 'thermo.consts'
  INTEGER :: k


  REAL, PARAMETER :: rad2deg=180/3.14159265
  REAL :: vpres
  INCLUDE 'thermo.stfunc'

  !  Input: n, pres, zz, theta, qv, uu, vv
  !  Output: temp, tempc, tdew, tdewc, presmb, spd, dir, qs, rh, tv

  DO k=1,n

    temp(k) = theta(k) * (pres(k) * p00inv) ** rcp
    tempc(k)= temp(k) - tfrz
    IF (qv(k) <= 0.) THEN
      tdewc(k) = -243.5
    ELSE
      vpres = fvpres(pres(k),qv(k))
      IF (vpres == 611.2) vpres = 611.3 !rlc 1995/12/07
      tdewc(k) = ftdewc(vpres)
    END IF
    tdew(k) = tdewc(k) + tfrz
    presmb(k) = pres(k) * 1.e-2

    spd(k)  = SQRT (uu(k)**2 + vv(k)**2 )
    dir(k) = 0.
    IF (spd(k) /= 0.) dir(k) = 270.- ATAN2(vv(k),uu(k)) * rad2deg
    dir(k)  = MOD (dir(k)+360.,360.)

    qs(k) = fmixrat(pres(k),fsvpres(temp(k)-tfrz))
    rh(k) = frh_q(qv(k),qs(k))
    tv(k) = ftvirtnowl(temp(k),qv(k))

  END DO


END SUBROUTINE derive
!
!
!
!
!                   ########################################
!                   ########################################
!                   ########################################
!                   ########                        ########
!                   ########      PLOT_STRINGS      ########
!                   ########                        ########
!                   ########################################
!                   ########################################
!                   ########################################
!
!

SUBROUTINE plot_strings (sndfile, stid, do_indices,                     &
      cape_use_t, gempak_fmt, modify_snd, plot_cb_parcel,               &
      plot_sfc_parcel, print_info, tv_wl,                               &
      isnd , nsnd , iyr      , imon     , iday     , ihr      ,         &
      imin     , px1      , px2      , py1      , py2      ,            &
      capet    , capetc   ,capetvnowl,capetvnowlc,capetvwl , capetvwlc, &
      cint     , cintc    , cintvnowl,cintvnowlc, cintvwl  , cintvwlc , &
      pmb_mod  , pmbcb    , preslcl  , preslnb  ,                       &
      sfctf    , sfctf0   , sfctdf   , sfctdf0  , zcb, zlcl, zlnb ,     &
      preslfc , zlfc , convtemp , precipwat , liftedindex ,             &
      totaltotals , kindex , sweat , wetbulbzero , brnshear ,           &
      lidstrength , dirmean , spmean , helicity , preciptype_me,        &
      preciptype_laps,nmax,n,sndcolor)
!
!======================================================================-
!      &&&&    G E N E R A L    I N F O R M A T I O N    &&&&
!======================================================================-
!
!  PURPOSE: This subroutine does x.
!
!  AUTHOR:  Richard Carpenter, Univ. of Oklahoma (rcarpenter@ou.edu)
!
!  HISTORY:
! 1997/04/24  First written
!
!  INPUT ARGUMENTS:
! x()   : Array
!
!  INPUT/OUTPUT ARGUMENTS:
! <none>
!
!  OUTPUT ARGUMENTS:
! <none>
!
!  I/O:
! <none>
!
!  SPECIAL REQUIREMENTS:
! <none>
!
!  OTHER INFORMATION:
! <none>
!
!======================================================================-
!      %%%%    D E C L A R A T I O N S    %%%%
!======================================================================-
!

  IMPLICIT NONE
  CHARACTER (LEN=*), INTENT(IN) :: sndfile, stid
  LOGICAL, INTENT(IN) :: &
      do_indices, cape_use_t, gempak_fmt, modify_snd, plot_cb_parcel, &
      plot_sfc_parcel, print_info, tv_wl
  INTEGER, INTENT(IN) :: isnd, nsnd, iyr, imon, iday, ihr, imin
  REAL, INTENT(IN) :: &
      px1, px2, py1, py2, capet, capetc, capetvnowl, capetvnowlc, capetvwl, &
      capetvwlc, cint, cintc, cintvnowl, cintvnowlc, cintvwl, cintvwlc, &
      pmb_mod, pmbcb, preslcl, preslnb, sfctf, sfctf0, sfctdf, sfctdf0, &
      zcb, zlcl, zlnb, preslfc, zlfc, convtemp, precipwat, liftedindex, &
      totaltotals, kindex, sweat, wetbulbzero, brnshear, lidstrength, &
      dirmean, spmean, helicity
  INTEGER, INTENT(IN) :: preciptype_me, nmax, preciptype_laps(nmax), n, sndcolor

  ! - - - - - - - - -  INCLUDE files - - - - - - - - - - - - - - - - - - -

  ! - - - - - - - - -  Constant declarations - - - - - - - - - - - - - - -

  ! - - - - - - - - -  Argument declarations - - - - - - - - - - - - - - -

  CHARACTER (LEN=85) :: key_laps
  REAL ::  spmkts, a, b

  ! - - - - - - - - -  Global/External declarations  - - - - - - - - - - -

  ! - - - - - - - - -  Local declarations  - - - - - - - - - - - - - - - -

  INTEGER :: i
  REAL :: xx, yy, pxoffset, pyoffset, pydelta, siz0
  CHARACTER :: runtype*20, file0*40
  CHARACTER (LEN=132) :: string, string2
  CHARACTER (LEN=2) :: prcp_string(4)
  DATA prcp_string / 'SN', 'IP', 'FZ RN', 'RN' /

  ! - - - - - - - - -  DATA statements - - - - - - - - - - - - - - - - - -

  DATA pyoffset/0.03/, pxoffset/0.01/, pydelta/0.02/, siz0/0.02/
  SAVE pyoffset, pxoffset

  ! - - - - - - - - -  Statement functions - - - - - - - - - - - - - - - -

  !6   & 12345678 , 12345678 , 12345678 , 12345678 , 12345678 , 12345678 ,
  !======================================================================-
  !   @@@@    E X E C U T A B L E    C O D E    @@@@
  !======================================================================-
  !
  PRINT *, 'PLOT_STRINGS'

  file0 = sndfile(INDEX(sndfile,'/',back=.true.)+1:) ! filename minus path

  ! Determine run type string (GEMPAK)

  IF (gempak_fmt) THEN
    IF (file0(:2) == 'ar') THEN
      i = INDEX(file0,'+')
      IF (i == 0) THEN
        string = '00'
      ELSE
        string = file0(INDEX(file0,'+')+1:INDEX(file0,'_',back=.true.)-1)
      END IF
      runtype = trim(string) // 'h FCST'
    ELSE IF (file0(:2) == 'ad') THEN
      runtype = 'ADAS'
    ELSE IF (file0(:2) == 'rc') THEN
      runtype = 'RUC'
    ELSE
      runtype = ' '
    END IF
  END IF


  ! Label for first file

  IF (isnd == 1) THEN
    CALL color (01)
    string = 'File: ' // file0
    IF (gempak_fmt) THEN
      WRITE (string,9988) stid, iyr,imon,iday, ihr,imin, trim(runtype)
      9988     FORMAT (a,1X, i4.2,2('/',i2.2), 1X,2I2.2, 'Z', 1X,a)
    END IF
    CALL plt2wrld (0.5,0.02,xx,yy)
    IF (do_indices) THEN
      a = 1.0
      b = 0.03
    ELSE
      a = 1.2
      b = 0.0
    END IF
    CALL xchmag (a*siz0)
    CALL xcharc (xx,yy+b,trim(string))
    PRINT *, 'Comment string = ',string
  END IF

  !  Key for LAPS Precip. Type Symbols (EMK, 1/18/98)

  IF (isnd==1 .AND. do_indices) THEN
    CALL xchmag (.75*siz0) ! (EMK)
    key_laps= 'LAPS Potnl Prcp Type:  R=.  ZR=~  IP=p  S=*'
    CALL xcharc(xx,yy-0.075,trim(key_laps))
  END IF

  ! other labels

  CALL color (sndcolor)
  CALL xchmag (0.8*siz0)

  IF (gempak_fmt) THEN
    WRITE (string,9987) isnd, stid, iyr,imon,iday, ihr,imin, trim(runtype)
    9987 FORMAT ('(',i1,') ',a,1X,i4.2,2('/',i2.2),1X,2I2.2,'Z',1X,a,1X,a)
    !CALL Plt2Wrld (px1+pxoffset,py1+0.088-0.02*isnd,xx,yy)
    CALL plt2wrld (px2-pxoffset,py2-pyoffset,xx,yy)
    yy = yy - 0.035 ! EMK
    CALL xcharr (xx,yy,trim(string))
    PRINT *, string
    pyoffset = pyoffset + pydelta
  END IF

  string = 'File: ' // trim(file0)
  CALL plt2wrld (px2-pxoffset,py2-pyoffset,xx,yy)
  yy = yy - 0.035 ! EMK
  CALL xcharr (xx,yy,trim(string))
  pyoffset = pyoffset + pydelta

  IF (plot_sfc_parcel .OR. print_info) THEN
    WRITE (string,'(a,i4,a,i5,a)')  &
        ' LCL:', nint(preslcl*1.e-2),' mb,', nint(zlcl),' m'
    CALL plt2wrld (px2-pxoffset,py2-pyoffset,xx,yy)
    yy = yy - 0.035 ! EMK
    CALL xcharr (xx,yy,trim(string))
    pyoffset= pyoffset + pydelta

    WRITE (string,'(a,i4,a,i5,a)')  &
        'LFC:', nint(preslfc*1.e-2),' mb,', nint(zlfc),' m'
    CALL plt2wrld (px2-pxoffset,py2-pyoffset,xx,yy)
    yy = yy - 0.035 ! EMK
    CALL xcharr (xx,yy,trim(string))
    pyoffset= pyoffset + pydelta

    WRITE (string,'(a,i4,a,i5,a)')  &
        'LNB:', nint(preslnb*1.e-2),' mb,', nint(zlnb),' m'
    CALL plt2wrld (px2-pxoffset,py2-pyoffset,xx,yy)
    yy = yy - 0.035 ! EMK
    CALL xcharr (xx,yy,trim(string))
    pyoffset    = pyoffset + pydelta

    IF (do_indices) THEN

      WRITE (string,'(a,f4.1,a,f4.1,a)')  &
          'Conv.T:',convtemp,' C, Prec.Wat.:',precipwat,' cm'
      CALL plt2wrld (px2-pxoffset,py2-pyoffset,xx,yy)
      yy = yy - 0.035 ! EMK
      CALL xcharr (xx,yy,trim(string))
      pyoffset = pyoffset + pydelta

      WRITE (string,'(a,f5.1,a,i2,a,i2)')  &
          'LI: ',liftedindex,'C, TT: ',nint(totaltotals),', KI: ' ,nint(kindex)
      CALL plt2wrld (px2-pxoffset,py2-pyoffset,xx,yy)
      yy = yy - 0.035 ! EMK
      CALL xcharr (xx,yy,trim(string))
      pyoffset    = pyoffset + pydelta

      WRITE (string,'(a,i4,a,i4,a)')  &
          'SWEAT: ',nint(sweat),', W/B zero: ',nint(wetbulbzero), ' m'
      CALL plt2wrld (px2-pxoffset,py2-pyoffset,xx,yy)
      yy = yy - 0.035 ! EMK
      CALL xcharr (xx,yy,trim(string))
      pyoffset    = pyoffset + pydelta

      WRITE (string,'(a,f4.1,a,f6.1,a)')  &
          'BRNshear: ',brnshear,' m/s, LSI:',lidstrength,'K'
      CALL plt2wrld (px2-pxoffset,py2-pyoffset,xx,yy)
      yy = yy - 0.035 ! EMK
      CALL xcharr (xx,yy,trim(string))
      pyoffset    = pyoffset + pydelta

      spmkts = spmean * 1.944
      WRITE (string,'(a,i3,a,i2,a)')  &
          'Mean Storm Motion: ',nint(dirmean),' at ',nint(spmean), ' m/s'
      CALL plt2wrld (px2-pxoffset,py2-pyoffset,xx,yy)
      yy = yy - 0.035 ! EMK
      CALL xcharr (xx,yy,trim(string))
      pyoffset    = pyoffset + pydelta

      WRITE (string,'(a,i3,a)') 'Storm Rel. Hel.:',nint(helicity),' m2/s2'
      CALL plt2wrld (px2-pxoffset,py2-pyoffset,xx,yy)
      yy = yy - 0.035 ! EMK
      CALL xcharr (xx,yy,trim(string))
      pyoffset    = pyoffset + pydelta

      ! Precipitation Types Added to Skew-T.  EMK 10/24/1997, 1/18/98

      WRITE (string,'(A)') 'MEta/LAPS Sfc Cond Prcp: ' // &
          TRIM(prcp_string(preciptype_me)) // ' / ' // &
          TRIM(prcp_string(preciptype_laps(1)))
      CALL plt2wrld (px2-pxoffset,py2-pyoffset,xx,yy)
      yy = yy - 0.035 ! EMK
      CALL xcharr (xx,yy,trim(string))
      pyoffset= pyoffset + pydelta

    END IF


  ! ORIGINAL PLOT LABELS  CMS

  !      WRITE (string,'(a,i4,a,F5.1,a)')
  !     >  'LCL:', NINT(presLCL*1.e-2),' mb,', zLCL/1000.0,' km'
  !      CALL Plt2Wrld (px2-pxoffset,py2-pyoffset,xx,yy)
  !      CALL XCHARR (xx,yy,TRIM(string))
  !      pyoffset = pyoffset + pydelta
  !c
  !      WRITE (string,'(a,i4,a,F5.1,a)')
  !     >  'LNB:', NINT(presLNB*1.e-2),' mb,', zLNB/1000.0,' km'
  !      CALL Plt2Wrld (px2-pxoffset,py2-pyoffset,xx,yy)
  !      CALL XCHARR (xx,yy,TRIM(string))
  !      pyoffset = pyoffset + pydelta
  !c
  ! CAPE string

    string = ' '
    string2 = ' '
    IF (cape_use_t) THEN
      WRITE (string,9951) nint(capet), nint(cint)
      WRITE (string2,9951) nint(capetc), nint(cintc)
    ELSE IF (tv_wl) THEN
      WRITE (string,9951) nint(capetvwl), nint(cintvwl)
      WRITE (string2,9951) nint(capetvwlc), nint(cintvwlc)
    ELSE
      WRITE (string,9951) nint(capetvnowl), nint(cintvnowl)
      WRITE (string2,9951) nint(capetvnowlc), nint(cintvnowlc)
    END IF
    9951 FORMAT ('CAPE:', i5, ' J/kg   CIN:', i4, ' J/kg')
    CALL plt2wrld (px2-pxoffset,py2-pyoffset,xx,yy)
    yy = yy - 0.035 ! EMK
    CALL xcharr (xx,yy,trim(string))
    pyoffset = pyoffset + pydelta

  END IF

  ! CAPE of CB parcel

  IF (plot_cb_parcel) THEN
    WRITE (string,'(a,i4,a,i5,a)')  &
        'Cloud base:', nint(pmbcb), ' mb,', nint(zcb), ' m'
    CALL plt2wrld (px2-pxoffset,py2-pyoffset,xx,yy)
    yy = yy - 0.035 ! EMK
    CALL xcharr (xx,yy,trim(string))
    pyoffset = pyoffset + pydelta

    IF (cape_use_t) THEN
      WRITE (string2,9951) nint(capetc), nint(cintc)
    ELSE IF (tv_wl) THEN
      WRITE (string2,9951) nint(capetvwlc), nint(cintvwlc)
    ELSE
      WRITE (string2,9951) nint(capetvnowlc), nint(cintvnowlc)
    END IF

    CALL plt2wrld (px2-pxoffset,py2-pyoffset,xx,yy)
    yy = yy - 0.035 ! EMK
    CALL xcharr (xx,yy,trim(string2))
    pyoffset = pyoffset + pydelta
  END IF

  ! Modified snd

  IF (modify_snd) THEN
    9011 FORMAT (a,f5.1,'/',f4.1,'F',1X, f5.1,'/',f4.1,'C')
    WRITE (string,9011) 'Init sfc T/Td:',  &
        sfctf0, sfctdf0, (sfctf0-32.)*5./9., (sfctdf0-32.)*5./9.
    CALL plt2wrld (px2-pxoffset,py2-pyoffset,xx,yy)
    yy = yy - 0.035 ! EMK
    CALL xcharr (xx,yy,trim(string))
    pyoffset = pyoffset + pydelta

    WRITE (string,9011) 'Mod sfc T/Td:',  &
        sfctf, sfctdf, (sfctf-32.)*5./9., (sfctdf-32.)*5./9.
    CALL plt2wrld (px2-pxoffset,py2-pyoffset,xx,yy)
    yy = yy - 0.035 ! EMK
    CALL xcharr (xx,yy,trim(string))
    pyoffset = pyoffset + pydelta

    WRITE (string,'(a,i4,a)') 'Top of modified layer:', nint(pmb_mod), ' mb'
    CALL plt2wrld (px2-pxoffset,py2-pyoffset,xx,yy)
    yy = yy - 0.035 ! EMK
    CALL xcharr (xx,yy,trim(string))
    pyoffset = pyoffset + pydelta
  END IF

  pyoffset = pyoffset + pydelta ! ready for next file

END SUBROUTINE plot_strings
