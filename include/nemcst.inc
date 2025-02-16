!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                      NEMCST.INC                      ######
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
!  Include file 'nemcst.inc' for ARPS
!
!  This file initializes values of some physical constants
!  that remain unchanged throughout a model run
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Jason J. Levit
!  11/15/1969
!
!  MODIFICATION HISTORY:
!
!  01/27/97 (J. Levit)
!  Corrected the constant 'qpmin' from .0005 to .0001.  The former
!  was erroneous, and the correction to the latter is now consistent
!  with the MWR paper which describes the scheme.
!
!  01/15/1998 (Yuhe Liu)
!  Removed parameter cp. It was defined also in phycst.inc
!
!-----------------------------------------------------------------------
!
  REAL :: pmas
  PARAMETER (pmas=1.0E-11)

  REAL :: qcmin
  PARAMETER (qcmin = .001)

  REAL :: qpmin
  PARAMETER (qpmin = .0001)

  REAL :: v2p
  PARAMETER (v2p = 25.0)

  REAL :: c2p
  PARAMETER (c2p = .0000167)

  REAL :: c2r
  PARAMETER (c2r = 16.7)

  REAL :: c2s
  PARAMETER (c2s = 33.3)

  REAL :: c2i
  PARAMETER (c2i = 16.7)

  REAL :: p2s
  PARAMETER (p2s = 5.0)

  REAL :: r2i
  PARAMETER (r2i = .00000833)

  REAL :: i2r
  PARAMETER (i2r = .00000833)

  REAL :: r2v
  PARAMETER (r2v = .001)

  REAL :: p2v
  PARAMETER (p2v = .004)

  REAL :: s2v
  PARAMETER (s2v = .002)

  REAL :: i2v
  PARAMETER (i2v = .001)

  REAL :: lvl
  PARAMETER (lvl = 2.500E6)      ! 2.350E6 in RAMS

  REAL :: lvi
  PARAMETER (lvi = 2.834E6)

  REAL :: lli
  PARAMETER (lli = 0.334E6)      ! 0.484E6 in RAMS
