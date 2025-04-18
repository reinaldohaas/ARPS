!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                    RADAREMUL.INC                     ######
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
!  Include file 'radaremul.inc' for ARPSENKF package.
!
!  This file contains some constants used for radar measurements.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Youngsun Jung
!  06/48/2010
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
  REAL, PARAMETER :: deg2rad = (3.14592654/180.)

  REAL :: cx(6)          ! (PI/6)*rho_qx
  REAL*8 :: alpha(6)     ! Shape parameter
