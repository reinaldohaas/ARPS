!
! reference time
!
  INTEGER :: ireftim
  COMMON /ref_time/ ireftim
!
! Thresholding parameters
!
  REAL :: gcrej,wgtthr,refpct,velpct,rvarthr,vvarthr
  REAL :: refmis,velmis,timmis

  PARAMETER (gcrej=1.6,        & ! ground clutter vel threshold
             wgtthr=1.0E-04,   & ! tot weight threshold for averaging
             refpct=0.25,      & ! coverage ratio limit
             velpct=0.25,      & ! coverage ratio limit
             rvarthr=225.,     & ! dBZ**2
             vvarthr=49.,      & ! m/s**2
             refmis=-999999.,  &
             velmis=-999999.,  &
             timmis=-999999.)

