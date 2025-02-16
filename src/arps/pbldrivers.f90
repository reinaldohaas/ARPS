!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE PBL_DRIVER                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE pbl_driver(tintopt,pblopt,nx,ny,nz,                          &
                  u,v,w,ptprt,pprt,ppi,qv,qscalar,                      &
                  ptbar,pbar,qvbar,rhostr,rhostru,rhostrv,deltat,       &
                  !rublten,rvblten,rthblten,rqvblten,rqcblten,rqiblten,  &
                  zp3d,roufns,hfx,qfx,xland,ptsfc,                      &
                  pbldpth,kmv,rprntl,                                   &
                  u3d,v3d,th3d,p3d,qv3d,qc3d,qi3d,                      &
                  ust,br,psim,psih,wspd,istatus)

!-----------------------------------------------------------------------
!
! PURPOSE this is the driver layer for all PBL parameterization schemes
!
!-----------------------------------------------------------------------
!
! Author: Y. Wang and Xiaoming Hu (2012/11/29)
!   Xiaoming Hu provided the formula for fm and fh (psim, psih).
!
! MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------

  USE module_bl_ysu

  IMPLICIT NONE

  INCLUDE 'globcst.inc'     ! Global constants that control model
  INCLUDE 'phycst.inc'
  INCLUDE 'timelvls.inc'
  INCLUDE 'mp.inc'

  INTEGER, INTENT(IN)    :: tintopt, pblopt
  INTEGER, INTENT(IN)    :: nx, ny, nz           ! Number of grid points in 3 directions
  REAL,    INTENT(IN)    :: deltat

  REAL,    INTENT(IN)    :: u     (nx,ny,nz,nt)  ! Total u-velocity (m/s)
  REAL,    INTENT(IN)    :: v     (nx,ny,nz,nt)  ! Total v-velocity (m/s)
  REAL,    INTENT(IN)    :: w     (nx,ny,nz,nt)  ! Total v-velocity (m/s)

  REAL,    INTENT(IN)    :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature (K)
  REAL,    INTENT(IN)    :: pprt  (nx,ny,nz,nt)  ! Perturbation pressure (Pascal)

  REAL,    INTENT(IN)    :: qv    (nx,ny,nz,nt)  ! Water vapor specific humidity (kg/kg)
  REAL,    INTENT(IN)    :: qscalar(nx,ny,nz,nt,nscalar)

  REAL,    INTENT(IN)    :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL,    INTENT(IN)    :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal).
  REAL,    INTENT(IN)    :: qvbar (nx,ny,nz)     ! Base state water vapor specific humidity
                                                 ! (kg/kg)
  REAL,    INTENT(IN)    :: ppi   (nx,ny,nz)     ! Exner function

  REAL,    INTENT(IN)    :: zp3d  (nx,ny,nz)     ! The physical height coordinate defined at
                                                 ! w-point of the staggered grid.
  REAL,    INTENT(IN)    :: rhostr (nx,ny,nz)     ! Base state density rhobar times j3.
  REAL,    INTENT(IN)    :: rhostru(nx,ny,nz)     ! Base state density rhobar times j3 at U stagger.
  REAL,    INTENT(IN)    :: rhostrv(nx,ny,nz)     ! Base state density rhobar times j3 at V stagger.

  REAL,    INTENT(IN)    :: roufns(nx,ny)
  REAL,    INTENT(IN)    :: xland(nx,ny)
  REAL,    INTENT(IN)    :: ptsfc(nx,ny)
  REAL,    INTENT(IN)    :: hfx(nx,ny), qfx(nx,ny)

  REAL,    INTENT(INOUT) :: pbldpth(nx,ny,nt)    ! Planetary boundary layer depth (m)
  REAL,    INTENT(OUT)   :: kmv    (nx,ny,nz)    ! Vertical turb. mixing coef. for
                                                 ! momentum. ( m**2/s )
  REAL,    INTENT(OUT)   :: rprntl (nx,ny,nz)    ! Vertical turb. mixing coef. for
                                                 ! scalars. ( m**2/s )

  !REAL,    INTENT(OUT)   :: rublten(nx,ny,nz), rvblten(nx,ny,nz),       &
  !                          rthblten(nx,ny,nz),rqvblten(nx,ny,nz),      &
  !                          rqcblten(nx,ny,nz),rqiblten(nx,ny,nz)
  !
  REAL,    INTENT(INOUT) :: u3d(nx,ny,nz), v3d(nx,ny,nz),th3d(nx,ny,nz), &
                            p3d(nx,ny,nz), qv3d(nx,ny,nz),qc3d(nx,ny,nz),&
                            qi3d(nx,ny,nz)

  REAL,    INTENT(INOUT) :: br(nx,ny), ust(nx,ny), psim(nx,ny), psih(nx,ny), &
                            wspd(nx,ny)

  INTEGER, INTENT(OUT)   :: istatus

!-----------------------------------------------------------------------

  INTEGER :: i, j, k
  INTEGER :: tlevel
  INTEGER :: nxlg, nylg

  INTEGER :: ips, ipe, jps, jpe

  LOGICAL :: flag_qi

  REAL    :: tema, temb, psfc

  REAL, PARAMETER :: karman     = 0.4
  !REAL, PARAMETER :: g          = 9.81

  !REAL, PARAMETER :: rd         = 287.
  !REAL, PARAMETER :: cp         = 7.*rd/2.
  !REAL, PARAMETER :: rovcp      = rd/cp

  !REAL, PARAMETER :: p1000mb    = 1.0E5

  REAL :: ZA, GZ1OZ0
  REAL :: THX, GOVRTH
  REAL :: ZOL, RZOL
  INTEGER :: NZOL

  INTEGER :: n
  REAL    :: zoln, x, y
  REAL    :: PSIMTB(0:1000), PSIHTB(0:1000)

  REAL    :: DTG, PSIT

  REAL, ALLOCATABLE, SAVE :: mol(:,:)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  IF (pblopt == 0) RETURN

  IF (.NOT. ALLOCATED(mol)) THEN
    ALLOCATE(mol(nx,ny), STAT = istatus)
    mol(:,:) = 0.0
  END IF

!-----------------------------------------------------------------------
!
! indices
!
!-----------------------------------------------------------------------

  tlevel = tpast
  IF(tintegopt == 2 .or. tintegopt == 3) tlevel = tpresent

  nxlg = (nx-3)*nproc_x + 3
  nylg = (ny-3)*nproc_y + 3

  ips = 1; ipe = nx-1
  !IF (loc_x == 1)       ips = 1
  !IF (loc_x == nproc_x) ipe = nx

  jps = 1; jpe = ny-1
  !IF (loc_y == 1)       jps = 1
  !IF (loc_y == nproc_y) jpe = ny

!-----------------------------------------------------------------------
!
! 3D Arrays
!
!-----------------------------------------------------------------------

  DO k = 1,nz
    DO j = 1, ny
      DO i = 1, nx
        u3d(i,j,k)  = u(i,j,k,tlevel)
        v3d(i,j,k)  = v(i,j,k,tlevel)
        th3d(i,j,k) = ptprt(i,j,k,tlevel)+ptbar(i,j,k)
        p3d(i,j,k)  = pprt(i,j,k,tlevel)+pbar(i,j,k)
        qv3d(i,j,k) = qv(i,j,k,tlevel) ! +qvbar(i,j,k)
      END DO
    END DO
  END DO

  IF (P_QC > 0) THEN
    qc3d(:,:,:) = qscalar(:,:,:,tlevel,P_QC)
  ELSE
    qc3d(:,:,:) = 0.0
  END IF

  IF (P_QI > 0) THEN
    qi3d(:,:,:) = qscalar(:,:,:,tlevel,P_QI)
    flag_qi = .TRUE.
  ELSE
    qi3d(:,:,:) = 0.0
    flag_qi = .FALSE.
  END IF

!-----------------------------------------------------------------------
!
! Compute friction velocity and bulk Richardson number since they were
! were not passed out from sfcphysics.
!
!-----------------------------------------------------------------------

  DO j = 1, ny-1
    DO i = 1, nx-1
      tema = 0.25*(u3d(i+1,j,1)+u3d(i,j,1)+u3d(i+1,j,2)+u3d(i,j,2))
      temb = 0.25*(v3d(i,j+1,1)+v3d(i,j,1)+v3d(i,j+1,2)+v3d(i,j,2))
      wspd(i,j)=MAX(0.1,                                                &
                    SQRT(tema*tema+temb*temb+w(i,j,2,tlevel)*w(i,j,2,tlevel)) )
    END DO
  END DO

  CALL cuc_bulkri(nx,ny,nz,ips,ipe,jps,jpe,zp3d,                        &
                  roufns,wspd,ptsfc,th3d(:,:,2), br, ust )

!-----------------------------------------------------------------------
!
! Compute similarity stability function based on BR (Recommonded by Xiaoming Hu).
!
!under 4 types of conditions (module_sf_sfclay.F):
!
!!-----CLASS 1; STABLE (NIGHTTIME) CONDITIONS:   (BR .GE. 0.2)
!        PSIM(I)=-10.*GZ1OZ0(I)
!!    LOWER LIMIT ON PSI IN STABLE CONDITIONS
!        PSIM(I)=AMAX1(PSIM(I),-10.)
!        PSIH(I)=PSIM(I)
!
!!-----CLASS 2; DAMPED MECHANICAL TURBULENCE:  (BR .LT. 0.2 .AND. BR .GT. 0.0);
!        PSIM(I)=-5.0*BR(I)*GZ1OZ0(I)/(1.1-5.0*BR(I))
!!    LOWER LIMIT ON PSI IN STABLE CONDITIONS
!        PSIM(I)=AMAX1(PSIM(I),-10.)
!!.....AKB(1976), EQ(16).
!        PSIH(I)=PSIM(I)
!
!!-----CLASS 3; FORCED CONVECTION: BR .EQ. 0.0
!        PSIM(I)=0.0
!        PSIH(I)=PSIM(I)
!
!!-----CLASS 4; FREE CONVECTION:       BR .LT. 0.0
!        PSIM(I)=PSIMTB(NZOL)+RZOL*(PSIMTB(NZOL+1)-PSIMTB(NZOL))
!        PSIH(I)=PSIHTB(NZOL)+RZOL*(PSIHTB(NZOL+1)-PSIHTB(NZOL))
!
!-----------------------------------------------------------------------

!sfclayinit
  DO N=0,1000
     ZOLN=-FLOAT(N)*0.01
     X=(1-16.*ZOLN)**0.25
     PSIMTB(N)=2*ALOG(0.5*(1+X))+ALOG(0.5*(1+X*X))-2.*ATAN(X)+2.*ATAN(1.)
     Y=(1-16*ZOLN)**0.5
     PSIHTB(N)=2*ALOG(0.5*(1+Y))
  END DO


  DO j = 1, ny-1
    DO i = 1, nx-1

!-----COMPUTE THE HEIGHT OF FULL- AND HALF-SIGMA LEVELS ABOVE GROUND LEVEL
      ZA =  0.5*(zp3d(i,j,2)+zp3d(i,j,3)) - zp3d(i,j,2)

      GZ1OZ0 = ALOG( ZA/roufns(i,j) )
      THX=TH3D(I,J,2)

      GOVRTH = G / THX

!-----CLASS 1; STABLE (NIGHTTIME) CONDITIONS:   (BR .GE. 0.2)
      IF (br(i,j) >= 0.2) THEN
        PSIM(i,j)=-10.*GZ1OZ0
        ! LOWER LIMIT ON PSI IN STABLE CONDITIONS
        PSIM(i,j)=AMAX1(PSIM(i,j),-10.)
        PSIH(i,j)=PSIM(i,j)
      ELSE IF ( br(i,j) < 0.2 .AND. br(i,j) > 0.0 ) THEN
!-----CLASS 2; DAMPED MECHANICAL TURBULENCE:  (BR .LT. 0.2 .AND. BR .GT. 0.0);
        PSIM(i,j)=-5.0*BR(i,j)*GZ1OZ0/(1.1-5.0*BR(i,j))
        ! LOWER LIMIT ON PSI IN STABLE CONDITIONS
        PSIM(i,j)=AMAX1(PSIM(i,j),-10.)
        !.....AKB(1976), EQ(16).
        PSIH(i,j)=PSIM(i,j)
      ELSE IF ( ABS(br(i,j)) <= 1.0E-6) THEN
!-----CLASS 3; FORCED CONVECTION: BR .EQ. 0.0
        PSIM(i,j)=0.0
        PSIH(i,j)=PSIM(i,j)
      ELSE IF ( br(i,j) < 0.0) THEN
!-----CLASS 4; FREE CONVECTION:       BR .LT. 0.0

        IF(UST(i,j) .LT. 0.01) THEN
          ZOL=BR(I,J)*GZ1OZ0
        ELSE
          ZOL=KARMAN*GOVRTH*ZA*MOL(I,J)/(UST(I,J)*UST(I,J))
        ENDIF
        ZOL=AMIN1(ZOL,0.)
        ZOL=AMAX1(ZOL,-9.9999)
        NZOL=INT(-ZOL*100.)
        RZOL=-ZOL*100.-NZOL
        PSIM(I,J)=PSIMTB(NZOL)+RZOL*(PSIMTB(NZOL+1)-PSIMTB(NZOL))
        PSIH(I,J)=PSIHTB(NZOL)+RZOL*(PSIHTB(NZOL+1)-PSIHTB(NZOL))

!---LIMIT PSIH AND PSIM IN THE CASE OF THIN LAYERS AND HIGH ROUGHNESS
!---  THIS PREVENTS DENOMINATOR IN FLUXES FROM GETTING TOO SMALL
!       PSIH(I)=AMIN1(PSIH(I),0.9*GZ1OZ0(I))
!       PSIM(I)=AMIN1(PSIM(I),0.9*GZ1OZ0(I))
!write(0,*) i,j,ust(i,j),br(i,j),psih(i,j),psihtb(nzol),nzol
        PSIH(I,J)=AMIN1(PSIH(I,J),0.9*GZ1OZ0)
        PSIM(I,J)=AMIN1(PSIM(I,J),0.9*GZ1OZ0)

      END IF

      DTG = THX-ptsfc(i,j)
      PSIT=AMAX1(GZ1OZ0-PSIH(I,J),2.)
      MOL(I,J)=KARMAN*DTG/PSIT

      psih(i,j) = gz1oz0 - psih(i,j)  ! fh  What is actually used in WRFV3.4.1
      psim(i,j) = gz1oz0 - psim(i,j)  ! fm
    END DO
  END DO

!-----------------------------------------------------------------------
!
! Call PBL scheme
!
!-----------------------------------------------------------------------

  SELECT CASE (pblopt)

  CASE (1)

    CALL pbl_ysu(u3d,v3d,th3d,qv3d,qc3d,qi3d,p3d,ppi,                   &
                 !rublten,rvblten,rthblten,rqvblten,rqcblten,rqiblten,   &
                 zp3d,roufns,ust,pbldpth(:,:,2),psim,psih,              &
                 xland,hfx,qfx,br, wspd,                                &
                 deltat,flag_qi,kmv,rprntl,                             &
                 1,nxlg-1, 1,nylg-1, 1,nz-1,                            &
                 1,nx,     1,ny,     1,nz,                              &
                 ips,ipe,  jps,jpe,  1,nz-1 )

  !CASE (2)

  CASE DEFAULT

    WRITE(*,'(1x,a,I0,a)') 'Unknown pbl option : ', pblopt,'. Program aborting...'
    CALL arpsstop('Unknown pbloption',1)

  END SELECT

  !DO k = 1, nz-1
  !  DO j = 1, ny-1
  !    DO i = 1, nx-1
  !      rublten (i,j,k) = rhostru(i,j,k)*rublten (i,j,k)
  !      rvblten (i,j,k) = rhostrv(i,j,k)*rvblten (i,j,k)
  !      rthblten(i,j,k) = rhostr (i,j,k)*rthblten(i,j,k)
  !      rqvblten(i,j,k) = rhostr (i,j,k)*rqvblten(i,j,k)
  !      rqcblten(i,j,k) = rhostr (i,j,k)*rqcblten(i,j,k)
  !      rqiblten(i,j,k) = rhostr (i,j,k)*rqiblten(i,j,k)
  !    END DO
  !  END DO
  !END DO

!-----------------------------------------------------------------------
!
! Handle MPI boundaries
!
!-----------------------------------------------------------------------

  CALL edgfill(kmv,    1,nx,ips,ipe,1,ny,jps,jpe,1,nz,1,nz-1)
  CALL edgfill(rprntl, 1,nx,ips,ipe,1,ny,jps,jpe,1,nz,1,nz-1)

  !CALL edgfill(rthblten,1,nx,ips,ipe,1,ny,jps,jpe,1,nz,1,nz-1)
  !CALL edgfill(rqvblten,1,nx,ips,ipe,1,ny,jps,jpe,1,nz,1,nz-1)
  !CALL edgfill(rqcblten,1,nx,ips,ipe,1,ny,jps,jpe,1,nz,1,nz-1)
  !CALL edgfill(rqiblten,1,nx,ips,ipe,1,ny,jps,jpe,1,nz,1,nz-1)

  RETURN
END SUBROUTINE pbl_driver

!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE CUC_BULKRI                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE cuc_bulkri(nx,ny,nz,ibgn,iend,jbgn,jend,zp,                  &
                      roufns,wspd,ptsfc, pt1, bulkri, c_u)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Compute C_u (friction velocity / wind speed). The quantity C_U is used by
!  the subroutine SFCFLX to obtain surface fluxes for both unstable
!  and stable cases.
!
!-----------------------------------------------------------------------
!
!  AUTHORS: V. Wong, X. Song and N. Lin
!  9/10/1993
!
!  For the unstable case (Bulk Richardson number bulkri < 0 ), the
!  formulation is based on the paper "On the Analytical Solutions of
!  Flux-Profile relationships for the Atmospheric Surface Layer" by
!  D.W. Byun in J. of Applied Meteorlogy, July 1990, pp. 652-657.
!
!  For the stable case, the formulation is based on the
!  paper "A Short History of the Operational PBL - Parameterization
!  at ECMWF" by J.F.Louis, M. Tiedtke and J.F. Geleyn in " Workshop
!  on Planetary Boundary Layer Parameterization", a publication
!  by the European Centre for Medium Range Weather Forecasts,
!  25-27 Nov. 1981.
!
!  MODIFICATION HISTORY:
!
!  9/04/1994 (K. Brewster, V. Wong and X. Song)
!  Facelift.
!
!  2/27/95 (V. Wong and X. Song)
!
!  02/07/96 (V.Wong and X.Song)
!  Set an upper limiter, z1limit, for depth of the surface layer z1.
!
!  05/01/97 (V. Wong and X. Tan)
!  Changed the computation of stabp, use the formulation similar to
!  one used by Bynn(1990), the momentum and thermal roughness
!  lengths are different.
!
!  05/29/97 (V. Wong and X. Tan)
!  Modified the formulation considering the height of the surface
!  layer z1 may equal zero.
!
!  10/22/2012 (Y. Wang)
!  A copy of cuc from sfcphy3d.f90, but passes out bulk Richardson number
!  explicitly.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    ibgn     i-index where evaluation begins.
!    iend     i-index where evaluation ends.
!    jbgn     j-index where evaluation begins.
!    jend     j-index where evaluation ends.
!
!    zp       The physical height coordinate defined at w-point of
!             staggered grid.
!    roufns   Surface roughness
!
!    wspd     Surface wind speed (m/s), defined as
!             sqrt(usuf*usuf + vsuf*vsuf + wsuf*wsuf)
!    ptsfc    Potential temperature at the ground level (K)
!    pt1      Potential temperature at the 1st scalar point above
!             ground level, (K)
!
!    c_uneu   Friction velocity at neutral state
!
!  OUTPUT:
!
!    ustar    Friction velocity
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: nx,ny,nz          ! The number grid points in 3 directions

  INTEGER, INTENT(IN)  :: ibgn              ! i-index where evaluation begins.
  INTEGER, INTENT(IN)  :: iend              ! i-index where evaluation ends.
  INTEGER, INTENT(IN)  :: jbgn              ! j-index where evaluation begins.
  INTEGER, INTENT(IN)  :: jend              ! j-index where evaluation ends.

  REAL,    INTENT(IN)  :: zp    (nx,ny,nz)     ! The physical height coordinate
                                               ! defined at w-point of staggered grid.
  REAL,    INTENT(IN)  :: roufns(nx,ny)        ! Surface roughness length

  REAL,    INTENT(IN)  :: wspd  (nx,ny)        ! Surface wind speed (m/s)

  REAL,    INTENT(IN)  :: ptsfc(nx,ny)         ! Potential temperature at the ground
                                               ! level (K)
  REAL,    INTENT(IN)  :: pt1   (nx,ny)        ! Potential temperature at the
                                               ! 1st scalar
                                               ! point above ground level, (K)

  REAL,    INTENT(OUT) :: bulkri(nx,ny)        ! bulk Richardson number

  REAL,    INTENT(OUT) :: c_u(nx,ny)           ! Frictional velocity (m/s)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j

  REAL :: z1                   ! The height of 1st scalar point above
                               ! ground level (m)
  REAL :: dz0                  ! z1-roufns, where roufns is defined as
                               ! surface roughness length
  REAL :: stabp                ! Monin-Obukhov STABility Parameter
                               ! (zeta)
  REAL :: x1,x0,psim           ! Intermediate parameters needed
  REAL :: z1drou,qb3pb2
  REAL :: c7,c8
  REAL :: sb,qb,pb,thetab,tb   ! During  computations
  REAL :: a,b,c,d
  REAL :: tempan,sqrtqb

  REAL :: zt                   ! Thermal roughness length
  REAL :: dzt
  REAL :: ztdrou

  REAL :: z1droup, ztdroup

  REAL, PARAMETER :: epsilon = 1.0E-6

!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'sfcphycst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
!-----------------------------------------------------------------------
!
!  Following Byun (1990).
!
!-----------------------------------------------------------------------
!
  DO j=jbgn,jend
    DO i=ibgn,iend

      z1  = 0.5*(zp(i,j,3)-zp(i,j,2))
      z1  = MIN(z1,z1limit)

      zt = ztz0*roufns(i,j)           ! Original Zt formulation

!---------------------------------------------------------------------
!     Test of new thermal roughness equation (Chen/Dudhia 2001) (JAB)
!---------------------------------------------------------------------
! The following was commented out temporarily because variable psim
! is used before it is initilized. It need to be taken care of later.
!
!      dz0 = z1-roufns(i,j)
!      z1droup = (z1+roufns(i,j))/roufns(i,j)
!
!      bulkri = (g/ptsfc(i,j))*(pt1(i,j)-ptsfc(i,j))*dz0/                &
!               (wspd(i,j)*wspd(i,j))
!
!      IF (bulkri <= 0.0) THEN
!        c_u(i,j) =kv/(LOG(z1droup)-psim)
!      ELSE
!        a=kv/LOG(z1droup)
!        b=5.0
!        d=5.0
!        c=SQRT(1.0+d*bulkri)
!        c_u(i,j) = a/SQRT(1.0+2.0*b*bulkri/c)
!      END IF
!
!       zt = (0.4*c_u(i,j)/0.000024) + 100.0
!
!       IF (abs(zt) > epsilon ) THEN
!         zt = 1.0/zt
!       ELSE IF (abs(zt) <= epsilon ) THEN
!         zt = ztz0*roufns(i,j)
!       END IF
!
!----------------------------------------------------------------------

      dzt = z1-zt
      ztdrou = z1/zt
      ztdroup = (z1+zt)/zt

      dz0 = z1-roufns(i,j)
      z1drou = z1/roufns(i,j)
      z1droup = (z1+roufns(i,j))/roufns(i,j)

      bulkri(i,j) = (g/ptsfc(i,j))*(pt1(i,j)-ptsfc(i,j))*dz0/           &
                    (wspd(i,j)*wspd(i,j))

      IF (bulkri(i,j) <= 0.0) THEN
!
!-----------------------------------------------------------------------
!
!  Unstable case: See equations (28)-(34) in Byun (1990).
!
!-----------------------------------------------------------------------
!
        bulkri(i,j) = MAX (bulkri(i,j),-10.0)

        sb =bulkri(i,j)/prantl0l

        qb=oned9*(c1l+c2l*sb*sb)
        pb=oned54*(c3l+c4l*sb*sb)

        qb3pb2=qb**3-pb*pb
        c7 = (z1*dzt*LOG(z1droup)*LOG(z1droup))/(dz0*dz0*LOG(ztdroup))

        IF( qb3pb2 >= 0.0 ) THEN

          sqrtqb = SQRT(qb)
          tempan = MAX( -1.0, MIN( 1.0, pb/(sqrtqb**3) ) )

          thetab=ACOS(tempan)
          stabp =c7*(-2.0*sqrtqb*COS(thetab/3.0)+c5l)

        ELSE

          tb    =(SQRT(-qb3pb2)+ABS(pb))**oned3
          stabp =c7*(-(tb+qb/tb)+c5l)

        END IF
!
!-----------------------------------------------------------------------
!
!  According to equation (14) in Byun (1990).
!
!-----------------------------------------------------------------------
!
        c8=gammaml*stabp
        x1=(1. - c8)**0.25
        x0=(1. - c8/z1drou)**0.25

        psim=2.0*LOG((1.0+x1)/(1.0+x0))+LOG((1.+x1*x1)/(1.+x0*x0))-     &
             2.0*ATAN(x1)+2.0*ATAN(x0)

!
!-----------------------------------------------------------------------
!
!  Compute C_u via equation (10) in Byun (1981).
!
!-----------------------------------------------------------------------
!
        c_u(i,j) =kv/(LOG(z1droup)-psim)

      ELSE
!
!-----------------------------------------------------------------------
!
!  Stable case:
!
!-----------------------------------------------------------------------
!
        a=kv/LOG(z1droup)
        b=5.0
        d=5.0
        c=SQRT(1.0+d*bulkri(i,j))

        c_u(i,j) = a/SQRT(1.0+2.0*b*bulkri(i,j)/c)

      END IF

    END DO
  END DO

  RETURN
END SUBROUTINE cuc_bulkri
