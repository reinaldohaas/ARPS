
!
! adassat.inc
!
  INTEGER, PARAMETER :: mx_sat = 2         ! Maximum number of satellites
! 
  CHARACTER (LEN=256) :: ir_fname(mx_sat)
  CHARACTER (LEN=256) :: vis_fname(mx_sat)
  CHARACTER (LEN=256) :: ircalname(mx_sat)
  CHARACTER (LEN=256) :: viscalname(mx_sat)
  CHARACTER (LEN=6)   :: satnamir(mx_sat)
  CHARACTER (LEN=6)   :: satnamvis(mx_sat)
  COMMON /satell_vars1/  ir_fname,vis_fname,ircalname,viscalname,       &
                         satnamir,satnamvis
!
  INTEGER :: nirfiles,nvisfiles
  COMMON /satell_vars2/  nirfiles,nvisfiles
!
  REAL :: latsatir(mx_sat),lonsatir(mx_sat)
  REAL :: latsatvis(mx_sat),lonsatvis(mx_sat)
  COMMON /satell_vars3/  latsatir,lonsatir,latsatvis,lonsatvis
!
! Satellite Calibration Constants
!
  REAL :: centw(mx_sat),fk1(mx_sat),fk2(mx_sat)
  REAL :: bc1(mx_sat),bc2(mx_sat),bc2inv(mx_sat)
  COMMON /adas_ircal/ centw,fk1,fk2,bc1,bc2,bc2inv
