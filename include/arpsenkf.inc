!
!##################################################################
!##################################################################
!######                                                      ######
!######                    ARPSINKF.INC                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!---------------------------------------------------------------------------
!
!  Include file 'ARPSENKF.INC' for ensemble Kalman Filter data assimilation
!
!--------------------------------------------------------------------------
!
!
!  AUTHOR: Mingjing Tong
!  Recent update
!  09/20/2006.
!
!  MODIFIED: Youngsun Jung (07/20/2007)
!            Added polarimetric variables and new error model parameters
!            Ting Lei( 02/2009)
!            added a block for sfc enkf
!
!            Youngsun Jung (03/20/2009)
!            Unified radar processing block
!
!            Youngsun Jung (08/07/2009)
!            Added new parameters: covflgopt, inf_h_width, inf_v_width
!
!--------------------------------------------------------------------------
!
  INTEGER :: nen                ! Number of ensemble members
  CHARACTER (LEN=30) :: casenam ! case name
  REAL :: assim_time        ! The model time of the current analysis
  INTEGER :: anaopt         ! Option of analysis method
  REAL :: fatime            ! The model time of the first analysis

  INTEGER :: modopt         ! Option of forecast model

  INTEGER :: timeext
  INTEGER :: nr

  INTEGER :: mult_inflat    ! Option to do multiplicative covariance inflation
  INTEGER :: adapt_inflat   ! Option to do adaptive covariance inflation
  INTEGER :: add_inflat     ! Option to do additive covariance inflation
  INTEGER :: relax_inflat   ! Option to use relaxation covariance inflation method
  REAL :: zcovinf        ! Ref thredshold for covariance inflation
  REAL :: cinf           ! Covariance inflation factor
  REAL :: inf_h_width    ! Scale factor of transition zone width
  REAL :: inf_v_width    ! Scale factor of transition zone depth
  REAL :: rlxf           ! Covariance relaxation factor
  REAL :: taperf         ! Covariance inflation factor for tapering function
  REAL :: lnsigcovinfcutoff   !length for vertical taper in inflation in ln(sigma)
  INTEGER :: infp        ! Option to apply covariance inflation to pressure
  CHARACTER (LEN=100) :: pertdir  ! Directory saving random
                                  !  perturbations (add_inflat = 3 or 4)
  REAL :: mdeu           ! std of random perturbation added to u
  REAL :: mdev           ! std of random perturbation added to v
  REAL :: mdew           ! std of random perturbation added to w
  REAL :: mdeptprt       ! std of random perturbation added to ptprt
  REAL :: mdepprt        ! std of random perturbation added to pprt
  REAL :: mdeqv          ! std of random perturbation added to qv
  REAL :: mdeqc          ! std of random perturbation added to qc
  REAL :: mdeqr          ! std of random perturbation added to qr
  REAL :: mdeqi          ! std of random perturbation added to qi
  REAL :: mdeqs          ! std of random perturbation added to qs
  REAL :: mdeqg          ! std of random perturbation added to qg
  REAL :: mdeqh          ! std of random perturbation added to qh

  INTEGER :: time_filter  ! option to use filter for AEnSRF
  REAL :: r0time          ! radius of influence in time dimension

  INTEGER :: hrchcl_filter   !Option to use hierarchical filter (Anderson, 2007)
  INTEGER :: grpnum          ! number of groups used for hierarchical filter
  REAL :: ratio2dw           ! ratio between  hrchcl_filter and distance depend filter

  INTEGER :: filter  ! option to use filter when computing K
  REAL :: r0hvr      ! radius of influence in horizontal direction for Vr
  REAL :: r0vvr      ! radius of influence in vertical direction for Vr
  REAL :: r0hz       ! radius of influence in horizontal direction for Z
  REAL :: r0vz       ! radius of influence in vertical direction for Z

  INTEGER :: covflgopt ! Option to read in covariance flag file

  REAL :: Vr_updtime_uv     ! Update time for u, v when assimilating Vr
  REAL :: Vr_updtime_w      ! Update time for w when assimilating Vr
  REAL :: Vr_updtime_ptprt  ! Update time for ptprt when assimilating Vr
  REAL :: Vr_updtime_pprt   ! Update time for pprt when assimilating Vr
  REAL :: Vr_updtime_qv     ! Update time for qv when assimilating Vr
  REAL :: Vr_updtime_qcrish ! Update time for qc, qr, qi, qs, qh when assimilating Vr
  REAL :: Vr_updtime_ncrish ! Update time for Nt(c,r,i,s,h) when assimilating Vr

  REAL :: Z_updtime_uv     ! Update time for u, v when assimilating Z
  REAL :: Z_updtime_w      ! Update time for w when assimilating Z
  REAL :: Z_updtime_ptprt  ! Update time for ptprt when assimilating Z
  REAL :: Z_updtime_pprt   ! Update time for pprt when assimilating Z
  REAL :: Z_updtime_qv     ! Update time for qv when assimilating Z
  REAL :: Z_updtime_qc     ! Update time for qc when assimilating Z
  REAL :: Z_updtime_qi     ! Update time for qi when assimilating Z
  REAL :: Z_updtime_qrsh   ! Update time for qr, qs, qh when assimilating Z
  REAL :: Z_updtime_ncrish ! Update time for Nt(c,r,i,s,h) when assimilating Z

  REAL :: DP_updtime_uv ! Update time for u, v when assimilating dual-pol
  REAL :: DP_updtime_w ! Update time for w when assimilating dual-pol
  REAL :: DP_updtime_ptprt ! Update time for ptprt when assimilating dual-pol
  REAL :: DP_updtime_pprt ! Update time for pprt when assimilating dual-pol
  REAL :: DP_updtime_qv ! Update time for qv when assimilating dual-pol
  REAL :: DP_updtime_qc ! Update time for qc when assimilating dual-pol
  REAL :: DP_updtime_qi ! Update time for qi when assimilating dual-pol
  REAL :: DP_updtime_qrsh ! Update time for qr, qs, qh when assimilating dual-pol
  REAL :: DP_updtime_ncrish ! Update time for Nt(c,r,i,s,h) when assimilating dual-pol

  INTEGER :: DP_delay_opt ! Option to delay the assimilaton of polarimetric data
  REAL :: DP_update_delay ! Time delayed for polarimetric data

  INTEGER :: vrqcopt       ! Option to apply simple QC to Vr obs
  REAL    :: Vr_QC_delay   ! Time delayed for Vr QC
  REAL    :: Vr_factor     ! Multiplicative factor for Vr observation thresholding
  CHARACTER (LEN = 256) :: Vr_err_fn   ! File name that contains Vr error threshold

  REAL :: hgt_thr(5), vr_thr(5)  ! Arrays to store height and Vr error thresholds

  REAL :: updt_hgt       ! The upper bound of height, below which model
                         ! state variables are updated

  INTEGER :: enforcebc   ! Option to enforce BC

  REAL :: epsQ           ! Minimun mixing ratio to compensate spurious
                         ! error growth
  REAL :: epsN           ! Minimun number concentration to compensate
                         ! spurious error growth

  CHARACTER(LEN=256) :: inidirname
  COMMON /enkf001/ inidirname, nen, anaopt, mult_inflat, add_inflat, relax_inflat, &
                   adapt_inflat, infp, filter, enforcebc, timeext, nr,  &
                   epsQ, epsN,                                          &
                   modopt,time_filter,r0time,hrchcl_filter,grpnum,      &
                   ratio2dw

  COMMON /enkf002/ casenam, pertdir

  COMMON /enkf003/ assim_time,                                          &
                   fatime, zcovinf, cinf, inf_h_width, inf_v_width,     &
                   rlxf, taperf, lnsigcovinfcutoff, mdeu, mdev, mdew,   &
                   mdeptprt, mdepprt, mdeqv, mdeqc, mdeqr, mdeqi,       &
                   mdeqs, mdeqg, mdeqh

  COMMON /enkf004/ r0hvr, r0vvr, r0hz, r0vz, covflgopt

  COMMON /enkf005/ Vr_updtime_uv, Vr_updtime_w, Vr_updtime_ptprt,       &
                   Vr_updtime_pprt, Vr_updtime_qv, Vr_updtime_qcrish,   &
                   Vr_updtime_ncrish,                                   &
                   Z_updtime_uv, Z_updtime_w, Z_updtime_ptprt,          &
                   Z_updtime_pprt, Z_updtime_qv, Z_updtime_qc,          &
                   Z_updtime_qi, Z_updtime_qrsh, Z_updtime_ncrish,      &
                   DP_updtime_uv, DP_updtime_w, DP_updtime_ptprt,       &
                   DP_updtime_pprt, DP_updtime_qv, DP_updtime_qc,       &
                   DP_updtime_qi, DP_updtime_qrsh, DP_updtime_ncrish,   &
                   DP_delay_opt, DP_update_delay, vrqcopt,              &
                   Vr_QC_delay, Vr_factor, Vr_err_fn, hgt_thr, vr_thr,  &
                   updt_hgt

! Routine data related parameters

  CHARACTER (LEN = 256) :: cvndatadir, blackfile
  INTEGER :: sfcassim, sndassim, proassim, innovcvn
  INTEGER :: sfc_intopt

  INTEGER :: use_sfc_uv, use_sfc_t, use_sfc_td, use_sfc_p
  REAL    :: err_sfc_uv, err_sfc_t, err_sfc_td, err_sfc_p, err_sfc_t_meso
  REAL    :: dif_sfc_uv, dif_sfc_t, dif_sfc_td, dif_sfc_p

  INTEGER :: use_snd_uv, use_snd_t, use_snd_td, use_snd_p
  REAL    :: err_snd_uv, err_snd_t, err_snd_td, err_snd_p
  REAL    :: dif_snd_uv, dif_snd_t, dif_snd_td, dif_snd_p

  INTEGER :: use_pro_uv
  REAL    :: err_pro_uv
  REAL    :: dif_pro_uv

  INTEGER :: sfcweight, sndweight, proweight
  REAL    :: sfcr0h, sfcr0h_meso, sfcr0v, sndr0h, sndr0v, pror0h, pror0v

  INTEGER :: sfc_updt_u, sfc_updt_v, sfc_updt_w, sfc_updt_pprt, sfc_updt_ptprt,    &
             sfc_updt_qv, sfc_updt_qc, sfc_updt_qr, sfc_updt_qi,                   &
             sfc_updt_qs, sfc_updt_qh, sfc_updt_tsoil, sfc_updt_qsoil
  INTEGER :: snd_updt_u, snd_updt_v, snd_updt_w, snd_updt_pprt, snd_updt_ptprt,    &
             snd_updt_qv, snd_updt_qc, snd_updt_qr, snd_updt_qi,                   &
             snd_updt_qs, snd_updt_qh
  INTEGER :: pro_updt_u, pro_updt_v, pro_updt_w, pro_updt_pprt, pro_updt_ptprt,    &
             pro_updt_qv, pro_updt_qc, pro_updt_qr, pro_updt_qi,                   &
             pro_updt_qs, pro_updt_qh

  REAL    :: err_sfc_pt, err_sfc_qv, err_sfc_rh, dif_sfc_pt, dif_sfc_qv, dif_sfc_rh
  REAL    :: err_snd_pt, err_snd_qv, err_snd_rh, dif_snd_pt, dif_snd_qv, dif_snd_rh

  INTEGER :: cvn_iskip, cvn_jskip
  REAL    :: cvn_var

  COMMON /enkf_cvn01/ cvndatadir, blackfile,                                 &
                      sfcassim, sndassim, proassim, innovcvn

  COMMON /enkf_cvn02/ use_sfc_uv, use_sfc_t, use_sfc_td, use_sfc_p,          &
                      sfc_intopt,                                            &
                      err_sfc_uv, err_sfc_t, err_sfc_td, err_sfc_p, err_sfc_t_meso,      &
                      dif_sfc_uv, dif_sfc_t, dif_sfc_td, dif_sfc_p,          &
                      sfcr0h, sfcr0h_meso, sfcr0v, sfcweight,                &
                      sfc_updt_u, sfc_updt_v, sfc_updt_w, sfc_updt_pprt, sfc_updt_ptprt, &
                      sfc_updt_qv, sfc_updt_qc, sfc_updt_qr, sfc_updt_qi,    &
                      sfc_updt_qs, sfc_updt_qh, sfc_updt_tsoil, sfc_updt_qsoil
  COMMON /enkf_cvn03/ use_snd_uv, use_snd_t, use_snd_td, use_snd_p,          &
                      err_snd_uv, err_snd_t, err_snd_td, err_snd_p,          &
                      dif_snd_uv, dif_snd_t, dif_snd_td, dif_snd_p,          &
                      sndr0h, sndr0v, sndweight,                             &
                      snd_updt_u, snd_updt_v, snd_updt_w, snd_updt_pprt, snd_updt_ptprt, &
                      snd_updt_qv, snd_updt_qc, snd_updt_qr, snd_updt_qi,    &
                      snd_updt_qs, snd_updt_qh
  COMMON /enkf_cvn05/ use_pro_uv, err_pro_uv, dif_pro_uv,                    &
                      pror0h, pror0v, proweight,                             &
                      pro_updt_u, pro_updt_v, pro_updt_w, pro_updt_pprt, pro_updt_ptprt, &
                      pro_updt_qv, pro_updt_qc, pro_updt_qr, pro_updt_qi,    &
                      pro_updt_qs, pro_updt_qh

  COMMON /enkf_cvn06/ err_sfc_pt, err_sfc_qv, err_sfc_rh, dif_sfc_pt, dif_sfc_qv, dif_sfc_rh, &
                      err_snd_pt, err_snd_qv, err_snd_rh, dif_snd_pt, dif_snd_qv, dif_snd_rh

  COMMON /enkf_cvn07/ cvn_iskip, cvn_jskip, cvn_var
!

  INTEGER :: rdrdasrc
  INTEGER :: rdrfmver
  INTEGER :: ntiltmax
  REAL :: vrupbd
  REAL :: reflwbd
  CHARACTER (LEN = 256) :: rdrobsdir
  CHARACTER (LEN = 256) :: trusimdir
  CHARACTER (LEN = 100) :: trusimnam
  CHARACTER (LEN = 256) :: hdmpdir
  CHARACTER (LEN = 256) :: hdmpfheader
  INTEGER :: state2obsopt
  INTEGER :: crsbmvopt
  INTEGER :: rmsdopt
  INTEGER :: histsctr
  INTEGER :: intrpopt

  COMMON /radar01/ rdrdasrc, rdrfmver, ntiltmax, state2obsopt, crsbmvopt, &
                   rmsdopt, histsctr,intrpopt

  COMMON /radar02/ rdrobsdir, trusimdir, trusimnam, hdmpdir, hdmpfheader

  COMMON /radar03/ vrupbd, reflwbd

  INTEGER :: vtopt   ! Option to include terminal velocity
  INTEGER :: zwopt   ! Option to apply Ref weight when calculating Vr
  INTEGER :: rfopt   ! Option for radar scattering calculation method
  CHARACTER (LEN = 256) :: rsadir
  REAL    :: wavelen
  INTEGER :: dsdparaopt
  INTEGER :: randomerropt ! Option to add random errors to Vr and Z obs
  INTEGER :: seedobs ! seed for generating obs perturbation
  INTEGER :: nongaussopt  ! Option to the Non-Gaussian error in reflectivity
  INTEGER, PARAMETER :: ndualvarmax = 7   ! maximum number of dual-pol variables
  INTEGER, PARAMETER :: rdrstdmax = 200
  REAL    :: rdrstd_vr (rdrstdmax) ! std of Gaussian noise of radar data
  REAL    :: rdrstd_zhh(rdrstdmax) ! std of Gaussian noise of radar data
  REAL    :: rdrstd_zvv(rdrstdmax) ! std of Gaussian noise of radar data
  REAL    :: rdrstd_zdr(rdrstdmax) ! std of Gaussian noise of radar data
  REAL    :: rdrstd_zdp(rdrstdmax) ! std of Gaussian noise of radar data
  REAL    :: rdrstd_kdp(rdrstdmax) ! std of Gaussian noise of radar data
  REAL    :: rdrstd_rhv(rdrstdmax) ! std of Gaussian noise of radar data
  REAL    :: sigmaCO ! std of correlated error of reflectivity
  REAL    :: sigmaUC ! std of uncorrelated error of reflectivity
  INTEGER :: attenopt    ! attenuation option
  REAL :: drinteg        ! dr lenght of integration element (m)
  INTEGER :: numvar      ! Number of radar observations available
  CHARACTER (LEN=12) :: dualvarname(ndualvarmax)   ! name of dual-pol variables

  COMMON /radar04/ vtopt, zwopt, rfopt, rsadir, wavelen, dsdparaopt,  &
                   randomerropt, seedobs, nongaussopt, attenopt

  COMMON /radar05/ numvar, dualvarname, rdrstd_vr, rdrstd_zhh, rdrstd_zvv, &
                   rdrstd_zdr, rdrstd_zdp, rdrstd_kdp, rdrstd_rhv,         &
                   sigmaCO, sigmaUC, drinteg

  INTEGER :: paraestopt  ! Option to estimate parameters using EnKF
  INTEGER :: paranum     ! Number of parameters to be estimated
  INTEGER :: n0west
  INTEGER :: n0sest
  INTEGER :: n0gest
  INTEGER :: roqsest
  INTEGER :: roqgest

  INTEGER :: lgn0w
  INTEGER :: lgn0s
  INTEGER :: lgn0g
  INTEGER :: lgroqs
  INTEGER :: lgroqg

  INTEGER :: paestdata(ndualvarmax)
  REAL    :: paestthd(ndualvarmax)  ! The lower bound of each radar parameter
                                    ! for the data to be considered for
                                    ! parameter estimation
  REAL    :: updtpara_time  ! Update time for parameters
  INTEGER :: estsqnc     ! Option of estimation sequence

  INTEGER :: painoprtr   ! Optiona to use member or mean in observation operator
  REAL :: n0wt           ! Assumed true n0rain when paraestopt = 1
  REAL :: n0st           ! Assumed true n0snow when paraestopt = 1
  REAL :: n0gt           ! Assumed true n0hail when paraestopt = 1
  REAL :: roqst          ! Assumed true rhosnow when paraestopt = 1
  REAL :: roqgt          ! Assumed true rhohail when paraestopt = 1

  INTEGER :: parainflat  ! Option to do variance inflation for parameters

  REAL :: minstd_n0w
  REAL :: minstd_n0s
  REAL :: minstd_n0g
  REAL :: minstd_roqs
  REAL :: minstd_roqg

  REAL :: infn0w
  REAL :: infn0s
  REAL :: infn0g
  REAL :: infroqs
  REAL :: infroqg

  INTEGER :: corparaopt ! Option to calculate correlation coefficient
                        ! between parameter and Vr/Z
  INTEGER :: corflgopt  ! corcoe flag
  REAL :: cor0          ! correlation coefficient threshold
  INTEGER :: cornum     !

  COMMON /enkf006/ paraestopt,paranum,                                &
                   n0west,n0sest,n0gest,roqsest,                      &
                   roqgest,lgn0w,lgn0s,lgn0g,lgroqs,lgroqg,paestdata, &
                   estsqnc,painoprtr,parainflat,corparaopt,corflgopt, &
                   cornum

  COMMON /enkf007/ paestthd,updtpara_time,n0wt,n0st,n0gt,roqst,roqgt, &
                   minstd_n0w,minstd_n0s,minstd_n0g,minstd_roqs,      &
                   minstd_roqg,infn0w,infn0s,infn0g,infroqs,infroqg,  &
                   cor0

  INTEGER, PARAMETER :: rdrmax = 200
  INTEGER :: radardaopt
  INTEGER :: exptype
  INTEGER :: rdrcomopt, rdrfmt, numrdra
  CHARACTER (LEN=256) :: loc_rdrinfo
  INTEGER :: nrdrused
  CHARACTER (LEN=12) :: radarname(rdrmax)
  INTEGER :: ntwtype(rdrmax)
  INTEGER :: vcpmode(rdrmax)

  INTEGER,PARAMETER :: nsweepmax = 22
  INTEGER :: U_nsweep(rdrmax)
  REAL    :: U_elvswp(rdrmax,nsweepmax)
  INTEGER :: U_ngate(rdrmax)
  REAL    :: U_gatesp(rdrmax)
  REAL    :: U_beamwid(rdrmax)
  REAL    :: U_delaz(rdrmax)
  REAL    :: U_sradmul(rdrmax)

  INTEGER :: rdrlocopt(rdrmax)
  REAL :: xpos(rdrmax)
  REAL :: ypos(rdrmax)
  REAL :: hpos(rdrmax)
  INTEGER :: ngridx(rdrmax)
  INTEGER :: ngridy(rdrmax)
  INTEGER :: nzlwx(rdrmax)
  INTEGER :: nzlwy(rdrmax)
  INTEGER, DIMENSION(rdrmax) :: rdrdaopt1,rdrdaopt2,rdrdaopt3,rdrdaopt4,     &
                                rdrdaopt5,rdrdaopt6,rdrdaopt7
  REAL, DIMENSION(rdrmax) :: rdrthd1,rdrthd2,rdrthd3,rdrthd4,rdrthd5,        &
                             rdrthd6,rdrthd7
  REAL :: z04vr
  REAL :: Zupprlim
  INTEGER :: paramest(rdrmax)

  COMMON /rdarpara01/ radardaopt,exptype,loc_rdrinfo,nrdrused,radarname,     &
                     ntwtype,vcpmode,U_nsweep,U_elvswp,U_ngate,U_gatesp,     &
                     U_beamwid,U_delaz,U_sradmul,                            &
                     rdrlocopt,xpos,ypos,hpos,rdrcomopt,rdrfmt,numrdra,      &
                     ngridx,ngridy,nzlwx,nzlwy,                              &
                     rdrdaopt1,rdrdaopt2,rdrdaopt3,rdrdaopt4,                &
                     rdrdaopt5,rdrdaopt6,rdrdaopt7,                          &
                     rdrthd1,rdrthd2,rdrthd3,rdrthd4,rdrthd5,rdrthd6,rdrthd7,&
                     z04vr,Zupprlim,paramest

  !REAL :: datauplev   ! The upper boundary of CASA radar data

  INTEGER :: innovsmth
  INTEGER :: smnum
  INTEGER :: innovenv
  CHARACTER (LEN = 100) :: sndsimnam

  INTEGER :: sprd3Ddmp ! Option to dump 3D ensemble spread
  CHARACTER (LEN = 256) :: sprd3Ddir
  INTEGER :: rmsesprddmp   ! Option to dump rms errors

  REAL :: innov_thd_Z

  COMMON /enkf008/ innovsmth, smnum, innovenv, sprd3Ddmp, rmsesprddmp,  &
                   innov_thd_Z

  COMMON /enkf009/ sndsimnam, sprd3Ddir

  REAL :: refthrsld
  INTEGER :: rmsfcst
  INTEGER :: rmsxbgn
  INTEGER :: rmsxend
  INTEGER :: rmsybgn
  INTEGER :: rmsyend
  INTEGER :: rmszbgn
  INTEGER :: rmszend

  INTEGER :: wrtmeanopt

  COMMON /rms001/ refthrsld

  COMMON /rms002/ rmsfcst, rmsxbgn, rmsxend,  &
                  rmsybgn, rmsyend, rmszbgn, rmszend

  COMMON /rms003/ wrtmeanopt

  INTEGER, PARAMETER :: maxrdr = 10, maxenda = 20
  INTEGER :: endtardopt
  CHARACTER (LEN = 10) :: enpostfix
  INTEGER :: meandmp
  INTEGER :: nenfile
  CHARACTER (LEN = 256) :: engrdbasfn
  CHARACTER (LEN = 256) :: endafile(maxenda)
  INTEGER :: obsdtaopt
  INTEGER :: rdrnum
  INTEGER :: rdrdanam
  CHARACTER (LEN = 256) :: rdrdaheader(maxrdr)
  CHARACTER (LEN = 256) :: rdrdafile(maxrdr)
  INTEGER :: rdrdtaopt
  INTEGER, PARAMETER :: maxobs = 10
  INTEGER :: covcoropt
  INTEGER :: covvrnum, covznum
  CHARACTER (LEN = 256) :: dirnamevrc, dirnamezc
  INTEGER, DIMENSION(maxobs) :: ivr, jvr, kvr, iz, jz, kz
  INTEGER, DIMENSION(maxobs) :: corvru, corvrv, corvrw, corvrptprt, corvrpprt
  INTEGER, DIMENSION(maxobs) :: corvrqv, corvrqc, corvrqr, corvrqi, corvrqs
  INTEGER, DIMENSION(maxobs) :: corvrqg, corvrqh, corvrvag, corvrvae
  INTEGER, DIMENSION(maxobs) :: corzu, corzv, corzw, corzptprt, corzpprt
  INTEGER, DIMENSION(maxobs) :: corzqv, corzqc, corzqr, corzqi, corzqs
  INTEGER, DIMENSION(maxobs) :: corzqg, corzqh
  INTEGER :: spreadopt
  INTEGER :: sprdavg
  INTEGER :: mdivopt
  INTEGER :: ni1s, ni2s, nj1s, nj2s, nk1s, nk2s
  INTEGER :: adjwopt

  COMMON /enkf010/ endtardopt,meandmp,nenfile,obsdtaopt,rdrnum,         &
                   rdrdanam,rdrdtaopt,covcoropt,covvrnum,covznum,       &
                   ivr,jvr,kvr,iz,jz,kz,corvru,corvrv,corvrw,           &
                   corvrptprt,corvrpprt,corvrqv,corvrqc,corvrqr,        &
                   corvrqi,corvrqs,corvrqg,corvrqh,                     &
                   corvrvag,corvrvae,corzu,                             &
                   corzv,corzw,corzptprt,corzpprt,corzqv,corzqc,        &
                   corzqr,corzqi,corzqs,corzqg,corzqh,                  &
                   spreadopt,sprdavg,                                   &
                   mdivopt,ni1s,ni2s,nj1s,nj2s,nk1s,nk2s,adjwopt

  CHARACTER (LEN = 256) :: enkfdtadir

  COMMON /enkf011/ enpostfix,engrdbasfn,endafile,rdrdaheader,           &
                   rdrdafile,dirnamevrc, dirnamezc, enkfdtadir

  REAL :: tbgn_dmpin
  REAL :: tend_dmpin
  REAL :: tintv_dmpin

  COMMON /cycle2/ tbgn_dmpin, tend_dmpin, tintv_dmpin

  INTEGER, PARAMETER :: iniana_acct=1,inputobs_acct=2,                  &
       inputbg_acct=3, enkfsub_acct=4,  calbgvr_acct=5,                 &
       calbgref_acct=6,enkfmain_acct=7, miscen_acct=8,                  &
       fornxcyc_acct=9,outputem_acct=10,paraest_acct=11

  INTEGER, PARAMETER  :: maxcasenum = 30
  INTEGER :: sstvtopt
  INTEGER :: expn
  REAL    :: acntl
  REAL    :: apert(maxcasenum)
  CHARACTER(LEN = 20)  :: sstvtid(maxcasenum)
  CHARACTER(LEN = 256) :: cntlfilnam
  CHARACTER(LEN = 256) :: pertfilnam(maxcasenum)
  INTEGER :: vrzwrtopt

  COMMON /enkf012/ sstvtopt, expn, acntl, apert, sstvtid,               &
                   cntlfilnam, pertfilnam, vrzwrtopt
!-----------------------------------------------------------------------

