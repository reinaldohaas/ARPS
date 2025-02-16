!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                    INDTFLG.INC                       ######
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
  INTEGER :: grdin,basin,varin,mstin,rainin,prcin,icein,tkein,trbin,    &
             sfcin,landin,totin,radin,flxin,snowcin,snowin,             &
             nscalarin,p_qcin,p_qrin,p_qiin,p_qsin,p_qhin,p_qgin,       &
             p_ncin,p_nrin,p_niin,p_nsin,p_nhin,p_ngin,                 &
             p_zrin,p_ziin,p_zsin,p_zhin,p_zgin,                        &
             p_ccin

  COMMON /inflags/ grdin,basin,varin,mstin,rainin,prcin,icein,          &
          tkein,trbin,sfcin,landin,totin,radin,flxin,snowcin,snowin,    &
          nscalarin,p_qcin,p_qrin,p_qiin,p_qsin,p_qhin,p_qgin,          &
          p_ncin,p_nrin,p_niin,p_nsin,p_nhin,p_ngin,                    &
          p_zrin,p_ziin,p_zsin,p_zhin,p_zgin,                           &
          p_ccin

