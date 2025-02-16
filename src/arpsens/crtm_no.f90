!
!#######################################################################
!#######################################################################
!######                                                           ######
!######                     SUBROUTINE crtmpost                   ######
!######                                                           ######
!######                        Developed by                       ######
!######        Center for Analysis and Prediction of Storms       ######
!######                   University of Oklahoma                  ######
!######                                                           ######
!#######################################################################
!#######################################################################
!
SUBROUTINE crtmpost(nx,ny,nz,zps,temp,qv,rhoair,qscalar,                &
                    pres,tsoil,soiltyp,vegtyp,snowdpth,user_emis,emiss, &
           icrtm,isatid,nchannel,chbgn,chend,CRTM_BT,CRTM_OD,icloud)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Dummy function for crtmpost.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  02/15/2012
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INCLUDE 'mp.inc'

  INTEGER :: nx ! Number of grid points in the x-dir
  INTEGER :: ny ! Number of grid points in the y-dir
  INTEGER :: nz ! Number of grid points in the z-dir

  REAL :: zps(nx,ny,nz)

  REAL    :: pres(nx,ny,nz),temp(nx,ny,nz)
  REAL    :: qv(nx,ny,nz), rhoair(nx,ny,nz)
  REAL    :: qscalar(nx,ny,nz,1)
  REAL    :: tsoil(nx,ny),snowdpth(nx,ny),emiss(nx,ny)

  INTEGER :: soiltyp(nx,ny),vegtyp(nx,ny)

  INTEGER :: icrtm,isatid,nchannel,chbgn,chend
  INTEGER :: icloud
  INTEGER :: user_emis

  REAL    :: CRTM_BT(nx,ny,nchannel)
  REAL    :: CRTM_OD(nx,ny,nchannel)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  WRITE(*,'(1x,a,/,8x,a,/,8x,a,/)')                                     &
     'ERROR: You should NOT have been here because CRTM library was not linked.',           &
     'You can either link CRTM library by recompiling the program with option "-io crtm".', &
     'Or run your program with namelist variable "icrtm = 0".'

  RETURN
END SUBROUTINE crtmpost
