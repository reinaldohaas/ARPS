!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                     AGRICST.H                        ######
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
!  Include file 'agricst.h' for ARPS AGRI
!
!  This file contains miscellaneous information used by arpsagri.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  03/28/1997 Combined original files, rundat.h, verbose.h, and some
!             additional parameters into one file.
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Define parameters for namelist &arpsagri.
!
!  Notes: The following parameters were defined in nodal.inc:
!
!         intrat
!         intratt
!         kcheck
!         mxnest
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=80) :: runold ! Name of the old run if restart
  INTEGER :: nmlntho           ! Length of runold

  CHARACTER (LEN=80) :: runnew ! Name of this run if restart
  INTEGER :: nmlnthn           ! Length of runnew

  CHARACTER (LEN=80) :: machst  ! ??

  REAL :: rstime          ! Time at which to restart,
                          ! if rstart=.true.

  INTEGER :: nxc          ! Dimension size in X-direction for base
                          ! coarse grid
  INTEGER :: nyc          ! Dimension size in Y-direction for base
                          ! coarse grid
  INTEGER :: nzc          ! Dimension size in Z-direction for base
                          ! coarse grid

  INTEGER :: lexbc        ! Flag for external boundary condition for
                          ! base grid

  INTEGER :: levfix       ! Level above which to regrid at start
  INTEGER :: intrpodr     ! Order of spatial interpolation

  LOGICAL :: rstart       ! Restart from previous AGRI restart data?
  LOGICAL :: rstdump      ! Dump data for later restart at end?
  LOGICAL :: grdsrt       ! Regrid at start?

  LOGICAL :: verbose1     ! .true. -> output grid location information
  LOGICAL :: verbose2     ! .true. -> output grid descriptor information
  LOGICAL :: verbose3     ! .true. -> output detailed io information
  LOGICAL :: verbose4     ! .true. -> output detailed storage information
  LOGICAL :: verbose5     ! .true. -> output solver diagnostic info
  LOGICAL :: verbose6     ! .true. -> output random information we've
                          !           found helpful in the past

  INTEGER :: nfinelv
  INTEGER :: nfinelv_max
  PARAMETER ( nfinelv_max = 10 )

  INTEGER :: ngrdnew(nfinelv_max)
  INTEGER :: ngrdnew_max
  PARAMETER ( ngrdnew_max = 10 )

  REAL :: ixc(ngrdnew_max,nfinelv_max), ixln(ngrdnew_max,nfinelv_max)
  REAL :: jyc(ngrdnew_max,nfinelv_max), jyln(ngrdnew_max,nfinelv_max)
  REAL :: gangle(ngrdnew_max,nfinelv_max)

  COMMON / agri000 / runold,runnew,machst
  COMMON / agri001 / nmlntho,nmlnthn
  COMMON / agri002 / rstime
  COMMON / agri003 / rstart,rstdump,grdsrt

  COMMON / agri100 / nxc,nyc,nzc, lexbc,levfix,intrpodr,                &
                     nfinelv,ngrdnew
  COMMON / agri101 / ixc,ixln,jyc,jyln,gangle

  COMMON / agri200 / verbose1, verbose2, verbose3,                      &
                     verbose4, verbose5, verbose6
