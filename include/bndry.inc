!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                     BNDRY.INC                        ######
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
!  Include file 'bndry.inc' for ARPS
!
!  Boundary condition parameters for the east, west, north,
!  south, top and bottom boundaries.
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/01/1991
!
!  MODIFICATION HISTORY:
!
!  7/25/92 (MX)
!  Added parameter rbcopt
!
!  7/29/97.  Dan Weber
!  Added fftopt for distingushing the type of FFT used in the upper
!  radiation condition.
!
!-----------------------------------------------------------------------
!
  INTEGER :: lbcopt ! Lateral boundary condition option;
                    ! = 1, for internal determined LBC
                    ! = 2, for externally forced LBC

  INTEGER :: ebc    ! east   boundary condition option.
  INTEGER :: wbc    ! west   boundary condition option.
  INTEGER :: nbc    ! north  boundary condition option.
  INTEGER :: sbc    ! south  boundary condition option.
  INTEGER :: tbc    ! top    boundary condition option.
  INTEGER :: bbc    ! bottom boundary condition option.
  INTEGER :: fftopt ! FFT option for tbc=4 configuration.

  INTEGER :: ebc_global    ! east   boundary condition option.
  INTEGER :: wbc_global    ! west   boundary condition option.
  INTEGER :: nbc_global    ! north  boundary condition option.
  INTEGER :: sbc_global    ! south  boundary condition option.
!
!-----------------------------------------------------------------------
!
!  The boundary conditions are defined independent of each
!  boundary, but should be consistent.
!
!  When the value of the above parameters is:
!  = 1, rigid wall boundary condition
!  = 2, periodic boundary condition
!  = 3, zero gradient boundary condition (specified for all variables)
!  = 4, radiation (open) boundary condition
!  = 5, user specified boundary conditon (externally)
!
!-----------------------------------------------------------------------
!
  INTEGER :: rbcopt ! Radiation lateral boundary condition option:
! = 1, Klemp & Wihelmson type with constant phase
                     !      speed c = c_phase
! = 2, Orlanski type, with estimated phase speed c

  COMMON /bccntl/ lbcopt, ebc,wbc,nbc,sbc,tbc,bbc,rbcopt,fftopt
  COMMON /bcglobal/ ebc_global,wbc_global,nbc_global,sbc_global

  REAL :: c_phase   ! Fixed phase speed of waves propagating out of
                    ! lateral boundary. Used by option rbcopt=1.
  COMMON /c_gwsped/ c_phase

  REAL :: rlxlbc    ! Inflow boundary relaxation coefficient.
  COMMON /relaxbc/ rlxlbc

  INTEGER :: rbc_plbc ! Option for pressure LBC 
  COMMON /rbcplbc/ rbc_plbc

  INTEGER :: pdetrnd   ! Option for detrending pressure
  COMMON /pdetrend/ pdetrnd
