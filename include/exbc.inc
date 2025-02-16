!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                       EXBC.INC                       ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Include file 'exbc.inc' for ARPS
!
!  This file defines the parameters and common blocks for the
!  external boundary conditions.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  08/22/1994
!
!  MODIFICATION HISTORY:
!
!  05/09/1995 (Yuhe Liu)
!  Added a new namelist input parameter, brlxhw.
!  Changed definition of ngbrz from PARAMETER to NAMELIST input.
!  (See initpara3d.f)
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  When EXBC option is turned on, i.e., lbcopt=2, set
!
!    nxebc=nx, nyebc=ny, nzebc=nz.
!
!  They are defined in dims.inc.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Pointers to EXBC varaiables in the EXBC buffer array
!
!-----------------------------------------------------------------------
!
  INTEGER :: nu0exb,nv0exb,nw0exb,npt0exb,npr0exb,nqv0exb,              &
             nqscalar0exb(0:20)

  COMMON /exbcom0/nu0exb,nv0exb,nw0exb,npt0exb,                         &
                  npr0exb,nqv0exb,nqscalar0exb
!
!-----------------------------------------------------------------------
!
!  Pointers to EXBC time tendencies in the EXBC buffer array
!
!-----------------------------------------------------------------------
!
  INTEGER :: nudtexb,nvdtexb,nwdtexb,nptdtexb,nprdtexb,nqvdtexb,        &
             nqscalardtexb(0:20)

  COMMON /exbcom1/nudtexb,nvdtexb,nwdtexb,nptdtexb,                     &
                  nprdtexb,nqvdtexb,nqscalardtexb

  INTEGER :: abstinit              ! Initial time in absolute seconds
  INTEGER :: abststop              ! Stop-searching time in absolute
                                   ! seconds
  INTEGER :: abstfcst0             ! Last EXBC forecast time in
                                   ! absolute seconds
  INTEGER :: abstfcst              ! Current EXBC forecasttime in
                                   ! absolute seconds

  COMMON /exbcom2/ abstinit, abststop, abstfcst0, abstfcst

  INTEGER :: ubcrd                 ! Flag for u  in EXBC file
  INTEGER :: vbcrd                 ! Flag for v  in EXBC file
  INTEGER :: wbcrd                 ! Flag for w  in EXBC file
  INTEGER :: ptbcrd                ! Flag for pt in EXBC file
  INTEGER :: prbcrd                ! Flag for pr in EXBC file
  INTEGER :: qvbcrd                ! Flag for qv in EXBC file
  INTEGER :: qscalarbcrd(20)       ! Flag for qc in EXBC file

  COMMON /exbcom3/ubcrd,vbcrd,wbcrd,ptbcrd,prbcrd,                      &
                  qvbcrd,qscalarbcrd

  CHARACTER (LEN=256) :: exbcname         ! Prefix string of EXBC file name
  CHARACTER (LEN=19)  :: tinitebd         ! Initial time of EXBC files

  COMMON /exbcom4/ exbcname,tinitebd

  INTEGER :: tintvebd              ! Time interval of EXBC files
                                   ! in seconds
  COMMON /exbcom5/ tintvebd

  REAL :: cbcdmp                   ! Coefficient of EXBC relaxation
                                   ! damping function
  REAL :: cbcmix                   ! Non-dimensional coefficient of
                                   ! computational mixing in EXBC zone
  REAL :: cbcmixh                  ! Dimensional coefficient of
                                   ! computational mixing in EXBC zone
  REAL :: brlxhw                   ! Helf-width of the boundary
                                   ! relaxation function in term
                                   ! of the number of grid zones
                                   ! (a real number).

  COMMON /exbcom6/ cbcdmp,cbcmix,cbcmixh, brlxhw

  INTEGER :: ngbrz                 ! Number of grid points in the boundary
                                   ! relaxation zone

  COMMON /exbcom8/ ngbrz
