!
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  DIMS.INC                            ######
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
!  Include file 'dims.inc' for ARPS utility programs.
!
!  This file defines the array size related parameters that need to
!  set before compilation.
!
!-----------------------------------------------------------------------
!
!  Declaration of parameters:
!
!  nx, ny, nz: Dimensions of computatioal grid.
!
!
!  IMPORTANT NOTE: dims.inc is no longer used by ARPS, ADAS, or other
!  main ARPS programs.  It is currently used only for a few utility
!  programs which have not been converted to using array allocation.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx       ! Number of grid points in the x-direction
  INTEGER :: ny       ! Number of grid points in the y-direction
  INTEGER :: nz       ! Number of grid points in the z-direction
!
!-----------------------------------------------------------------------
!
!  Set the following parameters.
!
!-----------------------------------------------------------------------
!
  PARAMETER( nx = 67, ny = 67, nz = 35 )  ! See comment below:
