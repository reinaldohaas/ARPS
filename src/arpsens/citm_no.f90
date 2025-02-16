SUBROUTINE cira_rtm(isatnum,iemiss,chnl,nx,ny,nz,xlat,xlon,zps,         &
                    tsoil,qv,pres,temp,CRTM_OD,qs,qc,qg,qi,qr,          &
                    emiss,CIRA_BT,CIRA_RD)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Dummy function for cira_rtm.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Fanyou Kong
!  04/14/2012
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INCLUDE 'mp.inc'

  INTEGER :: isatnum,iemiss
  INTEGER :: nx ! Number of grid points in the x-dir
  INTEGER :: ny ! Number of grid points in the y-dir
  INTEGER :: nz ! Number of grid points in the z-dir
  REAL :: chnl  ! channel wavelength

  REAL    :: zps(nx,ny,nz)
  REAL    :: pres(nx,ny,nz),temp(nx,ny,nz)
  REAL    :: qv(nx,ny,nz)
  REAL    :: qc(nx,ny,nz),qr(nx,ny,nz),qi(nx,ny,nz)
  REAL    :: qs(nx,ny,nz),qg(nx,ny,nz)
  REAL    :: tsoil(nx,ny),emiss(nx,ny)
  REAL    :: xlat(nx,ny),xlon(nx,ny)
  REAL    :: CRTM_OD(nx,ny,nz)

  REAL    :: CIRA_BT(nx,ny)
  REAL    :: CIRA_RD(nx,ny)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  WRITE(*,'(1x,a,/,8x,a,/,8x,a,/)')                                     &
     'ERROR: You should NOT have been here because CITM library was not linked.',           &
     'You can either link CITM library by recompiling the program with option "-io citm".', &
     'Or run your program with namelist variable "icitm = 0".'

  RETURN
END SUBROUTINE cira_rtm
