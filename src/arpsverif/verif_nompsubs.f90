!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE VERIF_COLLECT             ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE VERIF_COLLECT(model_data,obsrv_data,nhisfile)
               
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  No action is required for non-MPI mode.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Kevin W. Thomas
!
!  Original Coding: 08/02/05
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'vericst.inc'

  INTEGER :: nhisfile

  REAL :: model_data(sfcmax,nhisfile,5)
  REAL :: obsrv_data(sfcmax,nhisfile,5)

  RETURN
END SUBROUTINE verif_collect
