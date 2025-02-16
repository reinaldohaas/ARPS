  module module_arpsArray
!
!
! 1.1 background state vector
!     -----------------------
!
  IMPLICIT NONE
  save
!
!
!-----------------------------------------------------------------------
!
!  Arrays defining model grid
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: x     (:)      ! The x-coord. of the physical and
                                      ! computational grid. Defined at u-point.
  REAL, ALLOCATABLE :: y     (:)      ! The y-coord. of the physical and
                                      ! computational grid. Defined at v-point.
  REAL, ALLOCATABLE :: z     (:)      ! The z-coord. of the computational grid.
                                      ! Defined at w-point on the staggered grid.
  REAL, ALLOCATABLE :: zp    (:,:,:)  ! The physical height coordinate defined at
                                      ! w-point of the staggered grid.
  REAL, ALLOCATABLE :: zpsoil(:,:,:)  ! Soil level depth.

  REAL, ALLOCATABLE :: hterain(:,:)   ! The height of the terrain.
  REAL, ALLOCATABLE :: mapfct(:,:,:)  ! Map factors at scalar, u and v points

  REAL, ALLOCATABLE :: j1    (:,:,:)  ! Coordinate transformation Jacobian defined
                                      ! as - d( zp )/d( x ).
  REAL, ALLOCATABLE :: j2    (:,:,:)  ! Coordinate transformation Jacobian defined
                                      ! as - d( zp )/d( y ).
  REAL, ALLOCATABLE :: j3    (:,:,:)  ! Coordinate transformation Jacobian defined
                                      ! as d( zp )/d( z ).
  REAL, ALLOCATABLE :: j3soil(:,:,:)  ! Coordinate transformation Jacobian defined
                                      ! as d( zpsoil )/d( z ).
  REAL, ALLOCATABLE :: aj3x  (:,:,:)  ! Coordinate transformation Jacobian defined
                                      ! as d( zp )/d( z ) AVERAGED IN THE X-DIR.
  REAL, ALLOCATABLE :: aj3y  (:,:,:)  ! Coordinate transformation Jacobian defined
                                      ! as d( zp )/d( z ) AVERAGED IN THE Y-DIR.
  REAL, ALLOCATABLE :: aj3z  (:,:,:)  ! Coordinate transformation Jacobian defined
                                      ! as d( zp )/d( z ) AVERAGED IN THE Z-DIR.
  REAL, ALLOCATABLE :: j3inv (:,:,:)  ! Inverse of j3
  REAL, ALLOCATABLE :: j3soilinv (:,:,:)! Inverse of j3soil
!
!-----------------------------------------------------------------------
!
!  ARPS Time-dependent variables
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: u     (:,:,:)  ! Total u-velocity (m/s)
  REAL, ALLOCATABLE :: v     (:,:,:)  ! Total v-velocity (m/s)
  REAL, ALLOCATABLE :: w     (:,:,:)  ! Total w-velocity (m/s)
  REAL, ALLOCATABLE :: wcont (:,:,:)  ! Contravariant vertical velocity (m/s)
  REAL, ALLOCATABLE :: ptprt (:,:,:)  ! Perturbation potential temperature (K)
  REAL, ALLOCATABLE :: pprt  (:,:,:)  ! Perturbation pressure (Pascal)

  REAL, ALLOCATABLE :: qv    (:,:,:)  ! Water vapor specific humidity (kg/kg)
  REAL, ALLOCATABLE :: qscalar(:,:,:,:)
  REAL, ALLOCATABLE :: tke   (:,:,:)  ! Turbulent Kinetic Energy ((m/s)**2)

  REAL, ALLOCATABLE :: ubar  (:,:,:)  ! Base state u-velocity (m/s)
  REAL, ALLOCATABLE :: vbar  (:,:,:)  ! Base state v-velocity (m/s)
  REAL, ALLOCATABLE :: ptbar (:,:,:)  ! Base state potential temperature (K)
  REAL, ALLOCATABLE :: pbar  (:,:,:)  ! Base state pressure (Pascal).
  REAL, ALLOCATABLE :: rhostr(:,:,:)  ! Base state density rhobar times j3.
  REAL, ALLOCATABLE :: qvbar (:,:,:)  ! Base state water vapor specific
                                      ! humidity(kg/kg)

  REAL, ALLOCATABLE :: udteb (:,:)    ! T-tendency of u at e-boundary (m/s**2)
  REAL, ALLOCATABLE :: udtwb (:,:)    ! T-tendency of u at w-boundary (m/s**2)
  REAL, ALLOCATABLE :: vdtnb (:,:)    ! T-tendency of v at n-boundary (m/s**2)
  REAL, ALLOCATABLE :: vdtsb (:,:)    ! T-tendency of v at s-boundary (m/s**2)

  REAL, ALLOCATABLE :: pdteb (:,:)    ! T-tendency of pprt at e-boundary (PASCAL/s)
  REAL, ALLOCATABLE :: pdtwb (:,:)    ! T-tendency of pprt at w-boundary (PASCAL/s)
  REAL, ALLOCATABLE :: pdtnb (:,:)    ! T-tendency of pprt at n-boundary (PASCAL/s)
  REAL, ALLOCATABLE :: pdtsb (:,:)    ! T-tendency of pprt at s-boundary (PASCAL/s)
  REAL, ALLOCATABLE :: trigs1(:)      ! Array containing pre-computed trig
                                      ! function for use in fft.
  REAL, ALLOCATABLE :: trigs2(:)      ! Array containing pre-computed trig
                                      ! function for use in fft.
  INTEGER, ALLOCATABLE :: ifax1(:)    ! Array containing the factors of nx.
  INTEGER, ALLOCATABLE :: ifax2(:)    ! Array containing the factors of ny.
  REAL, ALLOCATABLE :: vwork1 (:,:)   ! 2-D work array for fftopt=2.
  REAL, ALLOCATABLE :: vwork2 (:,:)   ! 2-D work array for fftopt=2.
  REAL, ALLOCATABLE :: wsave1 (:)     ! Work array for fftopt=2.
  REAL, ALLOCATABLE :: wsave2 (:)     ! Work array for fftopt=2.

  REAL, ALLOCATABLE :: ppi(:,:,:)     ! Exner function.
  REAL, ALLOCATABLE :: csndsq(:,:,:)  ! Speed of sound squared

  REAL, ALLOCATABLE :: radfrc(:,:,:)  ! Radiation forcing (K/s)
  REAL, ALLOCATABLE :: radsw(:,:)     ! Solar radiation reaching the surface
  REAL, ALLOCATABLE :: rnflx(:,:)     ! Net absorbed radiation by the surface
  REAL, ALLOCATABLE :: radswnet(:,:)  ! Net shortwave radiation
  REAL, ALLOCATABLE :: radlwin(:,:)   ! Incominging longwave radiation
!

  CONTAINS

  SUBROUTINE allocatearpsArray(nx,ny,nz,nzsoil)

    INTEGER:: nx,ny,nz,nzsoil,istatus

    INCLUDE 'globcst.inc'
!
!-----------------------------------------------------------------------
!
!  Allocate adas arrays
!
!-----------------------------------------------------------------------
!
    ALLOCATE(x(nx),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:x")
    ALLOCATE(y(ny),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:y")
    ALLOCATE(z(nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:z")

    ALLOCATE(hterain(nx,ny),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:hterain")
    ALLOCATE(mapfct (nx,ny,8),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:mapfct")

    ALLOCATE(zp  (nx,ny,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:zp")
    ALLOCATE(zpsoil(nx,ny,nzsoil),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:zpsoil")
    ALLOCATE(j1  (nx,ny,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:j1")
    ALLOCATE(j2  (nx,ny,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:j2")
    ALLOCATE(j3  (nx,ny,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:j3")
    ALLOCATE(j3soil(nx,ny,nzsoil),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:j3soil")
    ALLOCATE(aj3x(nx,ny,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:aj3x")
    ALLOCATE(aj3y(nx,ny,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:aj3y")
    ALLOCATE(aj3z(nx,ny,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:aj3z")
    ALLOCATE(j3inv(nx,ny,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:j3inv")
    ALLOCATE(j3soilinv(nx,ny,nzsoil),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:j3soilinv")

    ALLOCATE(u    (nx,ny,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:u")
    ALLOCATE(v    (nx,ny,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:v")
    ALLOCATE(w    (nx,ny,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:w")
    ALLOCATE(wcont(nx,ny,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:wcont")
    ALLOCATE(ptprt(nx,ny,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:ptprt")
    ALLOCATE(pprt (nx,ny,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:pprt")

    ALLOCATE(qv   (nx,ny,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:qv")
    ALLOCATE(qscalar(nx,ny,nz,nscalar),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:qscalar")
    ALLOCATE(tke  (nx,ny,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:tke")

    ALLOCATE(ubar (nx,ny,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:ubar")
    ALLOCATE(vbar (nx,ny,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:vbar")
    ALLOCATE(ptbar(nx,ny,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:ptbar")
    ALLOCATE(pbar (nx,ny,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:pbar")
    ALLOCATE(rhostr(nx,ny,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:rhostr")
    ALLOCATE(qvbar(nx,ny,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:qvbar")

    ALLOCATE(udteb(ny,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:udteb")
    ALLOCATE(udtwb(ny,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:udtbw")
    ALLOCATE(vdtnb(nx,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:vdtnb")
    ALLOCATE(vdtsb(nx,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:vdtsb")
    ALLOCATE(pdteb(ny,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:pdteb")
    ALLOCATE(pdtwb(ny,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:pdtwb")
    ALLOCATE(pdtnb(nx,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:pdtnb")
    ALLOCATE(pdtsb(nx,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:pdtsb")

    ALLOCATE(trigs1(3*(nx-1)/2+1),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:trigs1")
    ALLOCATE(trigs2(3*(ny-1)/2+1),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:trigs2")
    ALLOCATE(ifax1(13),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:ifax1")
    ALLOCATE(ifax2(13),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:ifax2")

    ALLOCATE(vwork1(nx+1,ny+1),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:vwork1")
    ALLOCATE(vwork2(ny,nx+1),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:vwork2")

    ALLOCATE(wsave1(3*(ny-1)+15),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:wsave1")
    ALLOCATE(wsave2(3*(nx-1)+15),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:wsave2")

    ALLOCATE(ppi   (nx,ny,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:ppi")
    ALLOCATE(csndsq(nx,ny,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:csndsq")

    ALLOCATE(radfrc(nx,ny,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:radfrc")
    ALLOCATE(radsw (nx,ny),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:radsw")
    ALLOCATE(rnflx (nx,ny),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:rnflx")
    ALLOCATE(radswnet (nx,ny),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:radswnet")
    ALLOCATE(radlwin (nx,ny),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:radlwin")

  END  SUBROUTINE allocateArpsArray

  SUBROUTINE deallocateArpsArray
!
    DEALLOCATE(x )
    DEALLOCATE(y )
    DEALLOCATE(z )

    DEALLOCATE(hterain )
    DEALLOCATE(mapfct )

    DEALLOCATE(zp )
    DEALLOCATE(zpsoil )
    DEALLOCATE(j1 )
    DEALLOCATE(j2 )
    DEALLOCATE(j3 )
    DEALLOCATE(j3soil )
    DEALLOCATE(aj3x )
    DEALLOCATE(aj3y )
    DEALLOCATE(aj3z )
    DEALLOCATE(j3inv )
    DEALLOCATE(j3soilinv )

    DEALLOCATE(u )
    DEALLOCATE(v )
    DEALLOCATE(w )
    DEALLOCATE(wcont )
    DEALLOCATE(ptprt )
    DEALLOCATE(pprt  )

    DEALLOCATE(qv    )
    DEALLOCATE(qscalar )
    DEALLOCATE(tke   )

    DEALLOCATE(ubar  )
    DEALLOCATE(vbar  )
    DEALLOCATE(ptbar )
    DEALLOCATE(pbar  )
    DEALLOCATE(rhostr )
    DEALLOCATE(qvbar )

    DEALLOCATE(udteb )
    DEALLOCATE(udtwb )
    DEALLOCATE(vdtnb )
    DEALLOCATE(vdtsb )
    DEALLOCATE(pdteb )
    DEALLOCATE(pdtwb )
    DEALLOCATE(pdtnb )
    DEALLOCATE(pdtsb )

    DEALLOCATE(trigs1 )
    DEALLOCATE(trigs2 )
    DEALLOCATE(ifax1 )
    DEALLOCATE(ifax2 )

    DEALLOCATE(vwork1 )
    DEALLOCATE(vwork2 )

    DEALLOCATE(wsave1 )
    DEALLOCATE(wsave2 )

    DEALLOCATE(ppi    )
    DEALLOCATE(csndsq )

    DEALLOCATE(radfrc )
    DEALLOCATE(radsw  )
    DEALLOCATE(rnflx  )
    DEALLOCATE(radswnet )
    DEALLOCATE(radlwin )

  END SUBROUTINE deallocateArpsArray

END module module_arpsArray
