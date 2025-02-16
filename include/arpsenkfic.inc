!
!##################################################################
!##################################################################
!######                                                      ######
!######                    ARPSENKFIC.INC                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!--------------------------------------------------------------------------
!
!  Include file 'ARPSENKFIC.INC' for arpsenkfic.f90
!
!--------------------------------------------------------------------------
!
!
!  AUTHOR: Mingjing Tong
!  09/20/2006.
!
!--------------------------------------------------------------------------
!

  INTEGER :: modelopt         ! Option of forecast model

  INTEGER :: iniprtopt, iniprt_ptprt, iniprt_qv, flip_qcrish, smoothopt
  REAL :: lhor, lver

  INTEGER :: prtrgn_uv, prtrgn_w, prtrgn_ptprt, prtrgn_pprt
  INTEGER :: prtrgn_qv, prtrgn_qli
  INTEGER :: prtibgn, prtiend, prtjbgn, prtjend
  INTEGER :: prtkbgn, prtkend, prtibgnu, prtiendu
  INTEGER :: prtjbgnv, prtjendv, prtkbgnw, prtkendw

  CHARACTER (LEN=128) :: rdrdtadir ! Directory saving radar data
  INTEGER :: rdrdatype       ! Spatial distribution of radar data
  INTEGER :: radarnum        ! number of radar data
  INTEGER, PARAMETER :: maxnum = 40
  CHARACTER (LEN=128) :: rdrdtanam(maxnum) ! Radar data file name
  INTEGER :: flgwrtopt       ! Option to dump grid flag file
  REAL :: z0_uv, z0_ptprt, z0_pprt, z0_qv, z0_w, z0_qli ! reflectivity
                             ! threshold used to confine initial 
                             ! perturbations for different variables
  REAL :: r0h_uv, r0v_uv, r0h_w, r0v_w, r0h_ptprt, r0v_ptprt
  REAL :: r0h_pprt, r0v_pprt, r0h_qv, r0v_qv, r0h_qli, r0v_qli

  REAL :: uvzbgn             ! The lowest level the perturbations will be
                             ! added for u and v (in meters)
  REAL :: uvzend             ! The highest level the perturbations will be
                             ! added for u and v (in meters)
  REAL :: wzbgn              ! The lowest level the perturbations will be
                             ! added for w (in meters)
  REAL :: wzend              ! The highest level the perturbations will be
                             ! added for w (in meters)
  REAL :: ptprtzbgn          ! The lowest level the perturbations will be
                             ! added for ptprt (in meters)
  REAL :: ptprtzend          ! The highest level the perturbations will be
                             ! added for ptprt (in meters)
  REAL :: pprtzbgn           ! The lowest level the perturbations will be
                             ! added for pprt (in meters)
  REAL :: pprtzend           ! The highest level the perturbations will be
                             ! added for  (in meters)
  REAL :: qvzbgn             ! The lowest level the perturbations will be
                             ! added for qv (in meters)
  REAL :: qvzend             ! The highest level the perturbations will be
                             ! added for qv (in meters)
  REAL :: qczbgn             ! The lowest level the perturbations will be
                             ! added for qc (in meters)
  REAL :: qczend             ! The highest level the perturbations will be
                             ! added for qc (in meters)
  REAL :: qrzbgn             ! The lowest level the perturbations will be
                             ! added for qr (in meters)
  REAL :: qrzend             ! The highest level the perturbations will be
                             ! added for qr (in meters)
  REAL :: qizbgn             ! The lowest level the perturbations will be
                             ! added for qi (in meters)
  REAL :: qizend             ! The highest level the perturbations will be
                             ! added for qi (in meters)
  REAL :: qszbgn             ! The lowest level the perturbations will be
                             ! added for qs (in meters)
  REAL :: qszend             ! The highest level the perturbations will be
                             ! added for qs (in meters)
  REAL :: qgzbgn             ! The lowest level the perturbations will be
                             ! added for qg (in meters)
  REAL :: qgzend             ! The highest level the perturbations will be
                             ! added for qg (in meters)
  REAL :: qhzbgn             ! The lowest level the perturbations will be
                             ! added for qh (in meters)
  REAL :: qhzend             ! The highest level the perturbations will be
                             ! added for qh (in meters)

  INTEGER :: seeds           ! seed for generating random errors
                             ! which is an arbitrary negative integer

  REAL :: stdu               ! User specified std of Gaussion noise for u
  REAL :: stdv               ! User specified std of G noise for v
  REAL :: stdw               ! User specified std of G noise for w
  REAL :: stdptprt           ! User specified std of G noise for ptprt
  REAL :: stdpprt            ! User specified std of G noise for pprt
  REAL :: stdqv              ! User specified std of G noise for qv
  REAL :: stdqc              ! User specified std of G noise for qc
  REAL :: stdqr              ! User specified std of G noise for qr
  REAL :: stdqi              ! User specified std of G noise for qi
  REAL :: stdqs              ! User specified std of G noise for qs
  REAL :: stdqg              ! User specified std of G noise for qg
  REAL :: stdqh              ! User specified std of G noise for qh

  INTEGER :: soilprtopt
  REAL :: stdtsoilsfc        ! User specified std of G noise for tsoilsfc
  REAL :: stdtsoildp         ! User specified std of G noise for tsoildp
  REAL :: stdqsoilsfc        ! User specified std of G noise for qsoilsfc
  REAL :: stdqsoildp         ! User specified std of G noise for qsoildp

  REAL :: stdqrelative       ! Used to rescale the random perturbations of
                             ! qv/qc/qr/qi/qs/qg/qh/qsoil

  INTEGER :: hydroadj        ! Option to apply hydrostatic relationship
                             ! Same as hydradj in arps.input

  INTEGER :: rndprtdmp

  COMMON /arpsenkfic001/ rdrdtadir,rdrdtanam

  COMMON /arpsenkfic002/ iniprtopt,iniprt_ptprt,iniprt_qv,          &
                         flip_qcrish,smoothopt,prtrgn_uv,prtrgn_w,  &
                         prtrgn_ptprt,prtrgn_pprt,prtrgn_qv,        &
                         prtrgn_qli,prtibgn,                        &
                         prtiend,prtjbgn,prtjend,prtkbgn,prtkend,   &
                         prtibgnu,prtiendu,prtjbgnv,prtjendv,       &
                         prtkbgnw,prtkendw,rdrdatype,               &
                         radarnum,seeds,stdqrelative,hydroadj,      &
                         rndprtdmp,modelopt 

  COMMON /arpsenkfic003/ lhor,lver,flgwrtopt,                       &
                         z0_uv,z0_w,z0_ptprt,z0_pprt,               &
                         z0_qv,z0_qli,r0h_uv,r0v_uv,r0h_w,r0v_w,    &
                         r0h_ptprt,r0v_ptprt,r0h_pprt,r0v_pprt,     &
                         r0h_qv,r0v_qv,r0h_qli,r0v_qli,uvzbgn,      &
                         uvzend,wzbgn,wzend,ptprtzbgn,ptprtzend,    &
                         pprtzbgn,pprtzend,qvzbgn,qvzend,qczbgn,    &
                         qczend,qrzbgn,qrzend,qizbgn,qizend,qszbgn, &
                         qszend,qgzbgn,qgzend,                      &
                         qhzbgn,qhzend,stdu,stdv,stdw,              &
                         stdptprt,stdpprt,stdqv,stdqc,stdqr,stdqi,  &
                         stdqs,stdqg,stdqh,soilprtopt,stdtsoilsfc,  &
                         stdtsoildp,stdqsoilsfc,stdqsoildp

  INTEGER :: paraperturb
  INTEGER :: seedp                  ! The seed for parameter perturbation generation
  CHARACTER (LEN=256) :: dirnamp
  CHARACTER (LEN=80) :: microparafn ! The name of the file containing
                                    ! micro parameters added by mtong
  REAL :: initn0w            ! initial value of tnw
  REAL :: initn0s            ! initial value of tns
  REAL :: initn0g            ! initial value of tng
  REAL :: initroqs           ! initial value of roqs
  REAL :: initroqg           ! initial value of roqg

  REAL :: stdn0w
  REAL :: stdn0s
  REAL :: stdn0g
  REAL :: stdroqs
  REAL :: stdroqg

  COMMON /arpsenkfic004/ paraperturb, seedp

  COMMON /arpsenkfic005/ dirnamp, microparafn 

  COMMON /arpsenkfic006/ initn0w,initn0s,initn0g,initroqs,initroqg, &
                         stdn0w,stdn0s,stdn0g,stdroqs,stdroqg

  INTEGER :: envprtopt
  REAL :: prttop
  REAL :: lev1
  INTEGER :: seede
  CHARACTER (LEN=256) :: direprt
  CHARACTER (LEN=80) :: envprtnam
  INTEGER :: temmstopt
  REAL :: stdtem, stdmst

  CHARACTER (LEN=100 ) :: outdumpdir
  CHARACTER (LEN=80  ) :: outname
  REAL :: tfgs

  COMMON /arpsenkfic007/ envprtopt, seede, temmstopt

  COMMON /arpsenkfic008/ prttop, lev1, stdtem, stdmst, tfgs

  COMMON /arpsenkfic009/ direprt, envprtnam, outdumpdir, outname
