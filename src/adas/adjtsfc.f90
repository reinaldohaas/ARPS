!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE ADJTSFC                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE adjtsfc(nx,ny,nz,nscalar,rbufsz,                             &
!           ptprt,pprt,qv,qc,qr,qi,qs,qh,                                &
           ptprt,pprt,qv,qscalar,                                       &
           ptbar,pbar,ppi,rhostr,                                       &
           x,y,z,zp,j3inv,                                              &
           soiltyp, tsfc, wetsfc,snowdpth,                              &
           radfrc, radsw, rnflx, radswnet, radlwin,                     &
           rsirbm,rsirdf,rsuvbm,rsuvdf,                                 &
           cosz, cosss,                                                 &
           fdirir,fdifir,fdirpar,fdifpar,                               &
           tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,                     &
           tem9,tem10,tem11,tem12,tem13,tem14,tem15,tem16,              &
           radbuf,sh,tem17)

!-----------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER :: nx,ny,nz,nscalar
  INTEGER :: rbufsz
!
!-----------------------------------------------------------------------
!
!  Define ARPS variables
!
!-----------------------------------------------------------------------
!
  REAL, INTENT(IN) :: ptprt (nx,ny,nz)
  REAL, INTENT(IN) :: pprt  (nx,ny,nz)
  REAL, INTENT(IN) :: qv    (nx,ny,nz)
!  REAL, INTENT(IN) :: qc    (nx,ny,nz)
!  REAL, INTENT(IN) :: qr    (nx,ny,nz)
!  REAL, INTENT(IN) :: qi    (nx,ny,nz)
!  REAL, INTENT(IN) :: qs    (nx,ny,nz)
!  REAL, INTENT(IN) :: qh    (nx,ny,nz)
  REAL, INTENT(IN) :: qscalar(nx,ny,nz,nscalar)

  REAL, INTENT(IN) :: ptbar (nx,ny,nz)
  REAL, INTENT(IN) :: pbar  (nx,ny,nz)
  REAL, INTENT(IN) :: ppi   (nx,ny,nz)
  REAL, INTENT(IN) :: rhostr(nx,ny,nz)

  REAL, INTENT(IN) :: x      (nx)
  REAL, INTENT(IN) :: y      (ny)
  REAL, INTENT(IN) :: z      (nz)

  REAL, INTENT(IN) :: zp    (nx,ny,nz)  ! The physical height coordinate defined at
                                        ! w-point of staggered grid.
  REAL, INTENT(IN) :: j3inv (nx,ny,nz)

  INTEGER :: soiltyp(nx,ny)             ! Soil type at each point
  REAL, INTENT(INOUT) :: tsfc   (nx,ny)
  REAL, INTENT(IN) :: wetsfc (nx,ny)    ! Surface soil moisture in the top 1 cm layer
  REAL, INTENT(IN) :: snowdpth(nx,ny)   ! Snow depth (m)

  REAL, INTENT(OUT) :: radfrc(nx,ny,nz) ! Radiation forcing (K/s)
  REAL, INTENT(OUT) :: radsw  (nx,ny)   ! Solar radiation down to the surface
  REAL, INTENT(OUT) :: rnflx  (nx,ny)   ! Net radiation flux absorbed by surface
  REAL, INTENT(OUT) :: radswnet  (nx,ny)! Net solar radiation
  REAL, INTENT(OUT) :: radlwin  (nx,ny) ! Incoming longwave radiation

  REAL, INTENT(OUT) :: rsirbm(nx,ny)    ! Solar IR surface albedo for beam radiation
  REAL, INTENT(OUT) :: rsirdf(nx,ny)    ! Solar IR surface albedo for diffuse radiation
  REAL, INTENT(OUT) :: rsuvbm(nx,ny)    ! Solar UV surface albedo for beam radiation
  REAL, INTENT(OUT) :: rsuvdf(nx,ny)    ! Solar UV surface albedo for diffuse radiation

  REAL, INTENT(OUT) :: cosz  (nx,ny)    ! Cosine of zenith
  REAL, INTENT(OUT) :: cosss (nx,ny)    ! Cosine of angle between sun light and
                                        ! surface terrain slope

  REAL, INTENT(OUT) :: fdirir (nx,ny)   ! all-sky direct downward IR flux
                                        ! (0.7-10 micron) at the surface
  REAL, INTENT(OUT) :: fdifir (nx,ny)   ! all-sky diffuse downward IR flux
                                        ! at the surface
  REAL, INTENT(OUT) :: fdirpar(nx,ny)   ! all-sky direct downward par flux
                                        ! (0.4-0.7 micron) at the surface
  REAL, INTENT(OUT) :: fdifpar(nx,ny)   ! all-sky diffuse downward par flux
                                        ! at the surface
  REAL, INTENT(INOUT) :: tem1(nx,ny,nz)
  REAL, INTENT(INOUT) :: tem2(nx,ny,nz)
  REAL, INTENT(INOUT) :: tem3(nx,ny,nz)
  REAL, INTENT(INOUT) :: tem4(nx,ny,nz)
  REAL, INTENT(INOUT) :: tem5(nx,ny,nz)
  REAL, INTENT(INOUT) :: tem6(nx,ny,nz)
  REAL, INTENT(INOUT) :: tem7(nx,ny,nz)
  REAL, INTENT(INOUT) :: tem8(nx,ny,nz)
  REAL, INTENT(INOUT) :: tem9(nx,ny,nz)
  REAL, INTENT(INOUT) :: tem10(nx,ny,nz)
  REAL, INTENT(INOUT) :: tem11(nx,ny,nz)
  REAL, INTENT(INOUT) :: tem12(nx,ny,nz)
  REAL, INTENT(INOUT) :: tem13(nx,ny,nz)
  REAL, INTENT(INOUT) :: tem14(nx,ny,nz)
  REAL, INTENT(INOUT) :: tem15(nx,ny,nz)
  REAL, INTENT(INOUT) :: tem16(nx,ny,nz)
  REAL, INTENT(INOUT) :: radbuf(rbufsz) ! temporary arrays used for radiation
                                        ! transfer computing
  REAL, INTENT(INOUT) :: tem17(nx,ny,nz) ! added by DTD
  REAL, INTENT(INOUT) :: sh(nx,ny)       ! added by DTD

!
! Misc Local Variables
!
  INTEGER :: i,j
  REAL :: zs2,zs3,tk2,tk3,dtdz,tdiff,t2m,p0inv
!
! Temperature tuning parameters based on Mesonet data studies
! These are from Marena, 2000 OASIS data.
! A relatively conservative slope for dT/d(rnflx)
! 
  REAL, PARAMETER :: ptd0 = -0.66     ! degrees K
  REAL, PARAMETER :: ptd1 =  0.011    ! degrees K / (W/m2)
  REAL, PARAMETER :: ptdifmax = 10.
  REAL, PARAMETER :: ptdifmin = -5.
!
! Include files
!
  INCLUDE 'phycst.inc'
!
! Misc initializations
!
  p0inv=1./p0
!
! Call radiation package to get radiation at the surface
!
  CALL radiation(nx,ny,nz,rbufsz,                                     &
!                 ptprt,pprt,qv,qc,qr,qi,qs,qh,                        &
                 ptprt,pprt,qv,qscalar,                               &
                 ptbar,pbar,ppi, rhostr,                              &
                 x,y,z,zp, j3inv,                                     &
                 soiltyp,tsfc,wetsfc,snowdpth,                        &
                 radfrc,radsw,rnflx,radswnet,radlwin,                 &
                 rsirbm,rsirdf,rsuvbm,rsuvdf,cosz,cosss,              &
                 fdirir,fdifir,fdirpar,fdifpar,                       &
                 tem1,tem2,tem3,tem4,tem5,                            &
                 tem6,tem7,tem8,tem9,tem10,                           &
                 tem11,tem12,tem13,tem14,tem15,tem16,                 &
                 radbuf,sh,tem17)
  DO j=1,ny-1
    DO i=1,nx-1
!
! Estimate temperature offset from net radiation.
! Extrapolate k=2 and k=3 temperatures to shelter height (2m)
! Then apply tdiff from regression of 2m temps vs. skin temp
!
      IF(soiltyp(i,j) < 12 .OR. soiltyp(i,j) > 13) THEN
        zs3=0.5*(zp(i,j,3)+zp(i,j,4))
        zs2=0.5*(zp(i,j,2)+zp(i,j,3))
        tk3=(ptbar(i,j,3)+ptprt(i,j,3))*                           &
            ((pbar(i,j,3)+pprt(i,j,3))*p0inv)**rddcp
        tk2=(ptbar(i,j,2)+ptprt(i,j,2))*                           &
            ((pbar(i,j,2)+pprt(i,j,2))*p0inv)**rddcp
        dtdz=(tk3-tk2)/(zs3-zs2)
        t2m=tk2+dtdz*((zp(i,j,2)+2.)-zs2)
        tdiff=ptd0+ptd1*rnflx(i,j)
        tsfc(i,j)=t2m+tdiff
      END IF
    END DO
  END DO
  RETURN
END SUBROUTINE adjtsfc
