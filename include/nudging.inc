!
!##################################################################
!##################################################################
!######                                                      ######
!######                      NUDGING.INC                     ######
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
!  Include file 'nudging.inc' for ARPS
!
!  This file contains variables related to the nudging data
!  assimilation.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  March, 1998
!
!  MODIFICATION HISTORY:
!
!
!-----------------------------------------------------------------------
!
  INTEGER :: nudgopt         ! Nudging adjustment
                             ! = 0 none
                             ! = 1 nudging adjustment
  INTEGER :: nudgstp         ! big time step interval from ndintvl

  COMMON /nudge01/ nudgopt,nudgstp

  REAL :: ndstart            ! Time (secs) to begin adjustment
  REAL :: ndstop             ! Time (secs) to stop adjustment
  REAL :: ndtime             ! Duration time (secs) of adjustment
  REAL :: ndintvl            ! Interval of nudging adjustment
  REAL :: ndgain             ! Gain multiplier to apply to increments
  REAL :: ndscale            ! Scale factor to apply to correction
                             ! at each application

  COMMON /nudge02/ ndstart, ndstop, ndtime, ndintvl, ndgain, ndscale

  CHARACTER (LEN=256)  :: incrfnam  ! Name of file containing analysis increments

  COMMON /nudge03/ incrfnam

  INTEGER :: nudgu,nudgv,nudgw,nudgp,nudgpt,nudgqv,                     &
          nudgqc,nudgqr,nudgqi,nudgqs,nudgqh

  COMMON /nudge04/ nudgu,nudgv,nudgw,nudgp,nudgpt,nudgqv,               &
          nudgqc,nudgqr,nudgqi,nudgqs,nudgqh
