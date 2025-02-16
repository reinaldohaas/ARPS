!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                    INDTFLG.F90                       ######
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
!  Common block holding the input data flags (control parameters).
!
!  Modifications:
!
!  3/12/96 (Ming Xue)
!  Added parameter tkein.
!
!-----------------------------------------------------------------------
!

MODULE indtflg

  IMPLICIT NONE
  SAVE

  INTEGER :: grdin,basin,varin,mstin,rainin,prcin,icein,tkein,trbin,    &
      sfcin,landin,totin,radin,flxin,snowin


END MODULE indtflg
