!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE ASSIMDRIV                   ######
!######                                                      ######
!######               Copyright (c) 1996                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE assimdriv(nx,ny,nz,nzsoil,x,y,z,zp,zpsoil,                   &
           u,v,w,wcont,ptprt,pprt,qv,qc,qr,qi,qs,qh,                    &
           tke,kmh,kmv,mapfct,                                          &
           ubar,vbar,ptbar,pbar,rhostr,qvbar,j1,j2,j3,j3soil,           &
           tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,tem9,                &
           tem10,tem11,tem12,tem13,tem14,tem15,tem16)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  This subroutine serves as a driver for data blending and
!  the variational velocity adjustment.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Limin Zhao and Alan Shapiro
!  16/01/96
!
!  MODIFICATION HISTORY:
!
!  01/17/96 (Steven Lazarus)
!  Added 4 new temporary arrays (tem10, tem11, tem12, and tem13).
!  These arrays are passed in from FRCUVW (as uforce, vforce,
!  wforce, and defsq).
!
!  02/15/96 (Limin Zhao)
!  Merged in the code for blending ADAS/model data with single
!  Doppler vector retrievals.
!
!  02/17/96 (Steve Lazarus)
!  Added the option to fill wcont after time interpolation
!  when ivar=0.
!
!  02/20/96 (Limin Zhao)
!  Added the codes for computing the optimal weights.
!
!  02/21/96 (Limin Zhao)
!  Added a seperate code for reading ADAS data.
!
!  02/24/96 (Limin Zhao)
!  Merged in the Vr hole-filler for hole filling the observed Vr.
!
!  04/20/96 (Limin Zhao)
!  Set an option to choose a 3-D hole-filler/2-D hole-filler.
!
!  04/23/96 (Limin Zhao)
!  Made an option to use ADAS/model values as B.C. in hole-filler,
!  which allows the background information be smoothly blended into
!  the small scale retrieval fields through Dirichlet BC.
!
! 1 June 2002 Eric Kemp
! Soil variable updates.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the z-direction (vertical)
!    nzsoil   Number of soil levels
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space (m)
!    zpsoil   Height of soil levels (m)
!    mapfct   Map factors at scalar, u and v points
!
!    u        x component of velocity (m/s)
!    v        y component of velocity (m/s)
!    w        Vertical component of velocity in Cartesian
!             coordinates (m/s)
!    wcont    Contravariant vertical velocity (m/s)
!
!    ptprt    Perturbation potential temperature (K)
!    pprt     Perturbation pressure (Pascal)
!    qv       Water vapor specific humidity (kg/kg)
!    qc       Cloud water mixing ratio (kg/kg)
!    qr       Rainwater mixing ratio (kg/kg)
!    qi       Cloud ice mixing ratio (kg/kg)
!    qs       Snow mixing ratio (kg/kg)
!    qh       Hail mixing ratio (kg/kg)
!
!    kmh      Horizontal turb. mixing coef. for momentum. ( m**2/s )
!    kmv      Vertical turb. mixing coef. for momentum. ( m**2/s )
!
!    tke      Turbulent Kinetic Energy ((m/s)**2)
!
!    ubar     Base state x velocity component (m/s)
!    vbar     Base state y velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    rhostr   Base state air density (kg/m**3) times j3
!    qvbar    Base state water vapor specific humidity (kg/kg)
!
!    j1       Coordinate transformation Jacobian -d(zp)/dx
!    j2       Coordinate transformation Jacobian -d(zp)/dy
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!    j3soil   Coordinate transformation Jacobian  d(zpsoil)/dz
!
!    tem1     Observed radial velocity
!    tem5     Background radial velocity
!
!---------------------------------------------------------------------
!
!  OUTPUT:
!
!    The adjusted velocity fields.
!
!---------------------------------------------------------------------
!
!  WORK ARRAYS:
!
!    tem2      Temporary work array.
!    tem3      Temporary work array.
!    tem4      Temporary work array.
!    tem6      Temporary work array.
!    tem8      Temporary work array.
!    tem9      Temporary work array.
!    tem10     Temporary work array.
!    tem11     Temporary work array.
!    tem12     Temporary work array.
!    tem13     Temporary work array.
!
!    (These arrays are defined and used locally (i.e. inside this
!    subroutine), they may also be passed into routines called by
!    this one. Exiting the call to this subroutine, these temporary
!    work arrays may be used for other purposes, and therefore their
!    contents may be overwritten. Please examine the usage of work
!    arrays before you alter the code.)
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
  IMPLICIT NONE     ! Force explicit declarations
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'      ! Model physical constants
  INCLUDE 'globcst.inc'     ! Model control constants
  INCLUDE 'assim.inc'       ! Assim/Retr control parameters
  INCLUDE 'bndry.inc'       ! Boundary condition parameters
!
!-----------------------------------------------------------------------
!
  INTEGER :: nt                ! The no. of t-levels of t-dependent arrays
  INTEGER :: tpast             ! Index of time level for the past time.
  INTEGER :: tpresent          ! Index of time level for the present time.
  INTEGER :: tfuture           ! Index of time level for the future time.

  REAL :: ptol

  PARAMETER (nt=3, tpast=1, tpresent=2, tfuture=3)

  PARAMETER (ptol=0.1)      ! changed from 0.1 for AMS99 SSW 12-11-98

  INTEGER :: nx, ny, nz        ! Number of grid points in x, y and z directions
  INTEGER :: nzsoil
  INTEGER :: nstyps

  INTEGER :: grdbas

  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.
  REAL :: zpsoil(nx,ny,nzsoil) ! Height of soil levels (m)
  
  REAL :: mapfct(nx,ny,3)      ! Map factors at scalar, u and v points

  REAL :: u     (nx,ny,nz,nt)  ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nt)  ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz,nt)  ! Total w-velocity (m/s)
  REAL :: wcont (nx,ny,nz)     ! Contravariant vertical velocity (m/s)

  REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz,nt)  ! Perturbation pressure (Pascal)
  REAL :: qv    (nx,ny,nz,nt)  ! Water vapor specific humidity (kg/kg)
  REAL :: qc    (nx,ny,nz,nt)  ! Cloud water mixing ratio (kg/kg)
  REAL :: qr    (nx,ny,nz,nt)  ! Rain water mixing ratio (kg/kg)
  REAL :: qi    (nx,ny,nz,nt)  ! Cloud ice mixing ratio (kg/kg)
  REAL :: qs    (nx,ny,nz,nt)  ! Snow mixing ratio (kg/kg)
  REAL :: qh    (nx,ny,nz,nt)  ! Hail mixing ratio (kg/kg)

  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal).
  REAL :: rhostr(nx,ny,nz)     ! Base state air density (kg/m**3) times j3
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific humidity
                               ! (kg/kg)
  REAL :: tke   (nx,ny,nz,nt)  ! Turbulent Kinetic Energy ((m/s)**2)

  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian -d(zp)/d(x)
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian -d(zp)/d(y)
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian  d(zp)/d(z)
  REAL :: j3soil(nx,ny,nzsoil) ! Coordinate transformation Jacobian  d(zpsoil)/d(z)

  REAL :: tem1(nx,ny,nz)
  REAL :: tem2(nx,ny,nz)
  REAL :: tem3(nx,ny,nz)
  REAL :: tem4(nx,ny,nz)
  REAL :: tem5(nx,ny,nz)
  REAL :: tem6(nx,ny,nz)
  REAL :: tem7(nx,ny,nz)
  REAL :: tem8(nx,ny,nz)
  REAL :: tem9(nx,ny,nz)
  REAL :: tem10(nx,ny,nz)
  REAL :: tem11(nx,ny,nz)
  REAL :: tem12(nx,ny,nz)
  REAL :: tem13(nx,ny,nz)
  REAL :: tem14(nx,ny,nz)
  REAL :: tem15(nx,ny,nz)
  REAL :: tem16(nx,ny,nz)
!   real tem17(nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: icount ! Parameter input flag
  DATA icount /1/
  SAVE icount

  INTEGER :: i, j, k, ir, jr
  INTEGER :: istor,jstor,kstor
  INTEGER :: nchout    ! I/O Channel of unformatted binary restart data
  INTEGER :: it,ix
  INTEGER :: itmax     ! Maximum number of iterations
  INTEGER :: istat
  INTEGER :: isrc
  INTEGER :: tim
  INTEGER :: idummy
  INTEGER :: itag      ! Counter indicating input file date/time
  INTEGER :: count
  INTEGER :: mdpt
  INTEGER :: iassim
  INTEGER :: lenstr

  REAL :: t1           ! Time of first data file
  REAL :: t2           ! Time of second data file
  REAL :: t3           ! Time of third data file

  REAL :: assimtim(100)! Time of (all) input data files
  REAL :: adastim

  REAL :: xs           ! Scalar coordinate x-component
  REAL :: ys           ! Scalar coordinate y-component
  REAL :: zs           ! Scalar coordinate z-component
  REAL :: xs1          ! Scalar coordinate x-component
  REAL :: ys1          ! Scalar coordinate y-component
  REAL :: zs1          ! Scalar coordinate z-component

  REAL :: xmove        ! (x,y) position of radar relative to
  REAL :: ymove        ! moving grid (grid origin at (0,0,0))

  REAL :: znz2         ! = z(nz-2) - zshift
  REAL :: znz1         ! = z(nz-1) - zshift
  REAL :: z2           ! = z(2)    - zshift

  REAL :: rad          ! Radius of sphere in Cartesian coord
  REAL :: radinv       ! Inverse radius magnitude

  REAL :: xor,yor

  INTEGER :: newebc    ! East boundary condition parameter.
  INTEGER :: newwbc    ! West boundary condition parameter.
  INTEGER :: newnbc    ! North boundary condition parameter.
  INTEGER :: newsbc    ! South boundary condition parameter.

  REAL, ALLOCATABLE :: tem4dsoil(:,:,:,:)
!
!-----------------------------------------------------------------------
!
!  Routines called:
!
!-----------------------------------------------------------------------
!
  EXTERNAL avgy
  EXTERNAL avgz
  EXTERNAL dtadump
  EXTERNAL gtdmpfn
  EXTERNAL bcu
  EXTERNAL bcv
  EXTERNAL rhouvw
  EXTERNAL lbcw
  EXTERNAL vbcw
  EXTERNAL assimrd
  EXTERNAL adasread
  EXTERNAL assimvel
  EXTERNAL assimfil
  EXTERNAL assimblnd
  EXTERNAL retrint
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!
!----------------------------------------------------------------------
!
!
!  The following is a brief summary of the assimilation flags
!  available within this routine.  For a more detailed review
!  of the options see ASSIMCON.F.
!
!  ivar   Variational adjustment flag:
!
!         = 0, NO variational adjustment at this model time step
!         = 1, Perform variational adjustment at this model time step
!
!  insrt  Insertion flag:
!
!         = 0, Do NOT insert velocities at this time step
!         = 1, Insert velocities at this time step
!
!  A matrix describing the assimilation options herein:
!
!       ivar  insrt  Comments
!       --------------------------------------------------------
!        0         0         Nothing - Exit assimilation mode
!        0         1         Direct insertion only
!        1         0         Variational adjustment only
!        1         1         Variational adjustment + direct insertion
!
!  NOTE: If irecov=1, then the ingest and time interpolation of
!        3 input data files occur herein.
!
!-----------------------------------------------------------------------
!
  DO k=1,nz
    DO j=1,ny
      DO i=1,nx
        tem1(i,j,k) = 0.0
        tem2(i,j,k) = 0.0
        tem3(i,j,k) = 0.0
        tem4(i,j,k) = 0.0
        tem5(i,j,k) = 0.0
        tem6(i,j,k) = 0.0
        tem7(i,j,k) = 0.0
        tem8(i,j,k) = 0.0
        tem9(i,j,k) = 0.0
        tem10(i,j,k) = 0.0
        tem11(i,j,k) = 0.0
        tem12(i,j,k) = 0.0
        tem13(i,j,k) = 0.0
        tem14(i,j,k) = 0.0
        tem15(i,j,k) = 0.0
        tem16(i,j,k) = 0.0
      END DO
    END DO
  END DO

  WRITE(6,*) 'the code enter assimdriv'
  WRITE(6,*) 'insrt,ivar,irecov = ',                                    &
              insrt,ivar,irecov

  IF ( insrt == 0.AND.ivar == 0.AND.irecov == 0) RETURN

  IF ( insrt == 0.AND.ivar == 0.AND.irecov == 1) GO TO 350

  WRITE(6,*) 'data assimilation is turned on'
!
!-----------------------------------------------------------------------
!
!  Read in model background fields or adas outputs.
!
!  If this is the first time to enter this subroutine and the flag
!  for ADAS data is on, read into ADAS data as the first guess fields
!  for the velocity variational adjustment.
!
!  Outputs are saved in tem2,tem3,tem4, which are total u,v and w
!  respectively.
!
!  Note: rhobar in tem10
!
!-----------------------------------------------------------------------
!
  tim = tpresent

  isrc = 3                ! Flags to read 1 input data file

  CALL adasread(nx,ny,nz,nzsoil,x,y,z,zp,zpsoil,                        &
                isrc,istat,adastim,                                     &
                ubar,vbar,pbar,ptbar,rhostr,qvbar,                      &
                u,v,w,pprt,ptprt,qv,qc,qr,qi,qs,qh,j1,j2,j3,            &
                tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,                &
                tem9,tem10)

  WRITE(6,*)'the code enter adasread successfully'

  PRINT*,'zp(10,10,1),zp(10,10,2),zp(10,10,3),zp(10,10,4):'
  PRINT*,zp(10,10,1),zp(10,10,2),zp(10,10,3),zp(10,10,4)

!     open(30,file='data.background',status='unknown',
!  :                     form='unformatted')
!     write(30) tem1
!     write(30) tem2
!     write(30) tem3
!     write(30) tem4
!     close(30)

!
!-------------------------------------------------------------------------
!
!  Read in observations or radar retrieval.
!
!  Note that with the exception of u,v, and w all prognostic
!  variables are directly overwritten. uprt, vprt, and wprt are
!  returned in arrays tem11,tem12, and tem13. u,v and w are
!  reconstructed from these and the ubar and vbar arrays. The
!  ingested model data is not necessarily the same as the current
!  model values.
!
!  Outputs are saved in tem1, tem11, tem12 and tem13, which are
!  Vr, (u,v,w) respectively.
!
!-------------------------------------------------------------------------
!
  isrc = 3                ! Flags to read 1 input data file

  CALL assimrd(nx,ny,nz,nzsoil,x,y,z,zp,zpsoil,                         &
               isrc,idummy,istat,assimtim,                              &
               ubar,vbar,pbar,ptbar,rhostr,qvbar,tem10,                 &
               u,v,w,pprt,ptprt,qv,qc,qr,qi,qs,qh,j1,j2,j3,             &
               tem1,tem5,tem6,tem7,tem8,tem9,tem14,                     &
               tem11,tem12,tem13)

  WRITE(6,*) 'the code enter assimrd successfully'

!   open(31,file='data.obs',status='unknown',form='unformatted')
!   write(31) tem1
!   write(31) tem11
!   write(31) tem12
!   write(31) tem13
!   close(31)
!
!-----------------------------------------------------------------------
!
!  Hole-fill u,v and w outside the rainwater field.  Hole-filling
!  is performed on the perturbation velocity fields only. tem7 serves
!  as a 'template' - delineating the region to be filled.  The
!  Dirichlet boundary conditions are set here, outside of POIS3D.
!
!-----------------------------------------------------------------------
!

  CALL assimfil(nx,ny,nz,x,y,z,zp,ptol,                                 &
                u,v,w,wcont,ptprt,pprt,qv,qc,qr,qi,qs,qh,               &
                ubar,vbar,ptbar,pbar,rhostr,qvbar,j1,j2,j3,             &
                tem5,tem6,tem7,tem8,tem9,tem10,tem11,tem12,             &
                tem13,tem14,tem15,tem16)

  WRITE(6,*) 'the code enter assimfil successfully'
!
!-----------------------------------------------------------------------
!
!  Hole-fill vr outside the rainwater field.  Hole-filling
!  is performed on the perturbation velocity fields only.
!  Vr is decomposited into Cartesian components (ur,vr,wr).
!  The Dirichlet boundary conditions are set here, outside
!  of POIS3D.
!
!  The hole-filled flag is saved in tem5.
!
!-----------------------------------------------------------------------
!
  CALL vrhole(nx,ny,nz,x,y,z,zp,                                        &
              tem14,tem15,tem16,                                        &
              tem1,ubar,vbar,ptol,                                      &
              tem5,tem6,tem7,tem8,tem9,tem10)

  WRITE(6,*)'the code enter vrhole successfully'

!   open(30,file='data.holefill',status='unknown',
!  :                     form='unformatted')
!   write(30) tem1
!   write(30) tem11
!   write(30) tem12
!   write(30) tem13
!   close(30)
!
!-------------------------------------------------------------------------
!
!    Blend the single Doppler retrieval velocity fields with ADAS
!    data. ADAS data is treated as the first guess field, and the OI
!    method is used to blend the two data sets together. The blended
!    velocity fields, which stored in tem7,tem8,tem9, will be sent
!    to the variational adjustment as background fields.
!
!    Note: The retrievals need be put on staggered grids before
!       blended with model background.
!
!-------------------------------------------------------------------------
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=2,nx-1
        tem7(i,j,k)=0.5*(tem11(i,j,k)+tem11(i-1,j,k))
      END DO
    END DO
  END DO

  DO k=1,nz-1
    DO j=1,ny-1
      tem7(1,j,k)=tem11(1,j,k)
      tem7(nx,j,k)=tem11(nx-1,j,k)
    END DO
  END DO

  DO k=1,nz-1
    DO j=2,ny-1
      DO i=1,nx-1
        tem8(i,j,k)=0.5*(tem12(i,j,k)+tem12(i,j-1,k))
      END DO
    END DO
  END DO

  DO k=1,nz-1
    DO i=1,nx-1
      tem8(i,1,k)=tem12(i,1,k)
      tem8(i,ny,k)=tem12(i,ny-1,k)
    END DO
  END DO

  DO k=2,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem9(i,j,k)=0.5*(tem13(i,j,k)+tem13(i,j,k-1))
      END DO
    END DO
  END DO

  DO i=1,nx-1
    DO j=1,ny-1
      tem9(i,j,1)=tem13(i,j,1)
      tem9(i,j,nz)=tem13(i,j,nz-1)
    END DO
  END DO

!   open(30,file='data.ObsStaggerred',status='unknown',
!  :                     form='unformatted')
!   write(30) tem1
!   write(30) tem7
!   write(30) tem8
!   write(30) tem9
!   close(30)

  CALL assimblnd(nx,ny,nz,                                              &
                 tem7,tem2,tem5,                                        &
                 tem6,tem10,tem11,tem12)

  CALL assimblnd(nx,ny,nz,                                              &
                 tem8,tem3,tem5,                                        &
                 tem6,tem10,tem11,tem13)

  WRITE(6,*) 'the code enter assimblnd successfully'

  DO k=1,nz
    DO j=1,ny
      DO i=1,nx
        u(i,j,k,tim) = tem12(i,j,k)
        v(i,j,k,tim) = tem13(i,j,k)
!
!       w(i,j,k,tim) = 0.7*tem9(i,j,k) + 0.3*tem4(i,j,k)
        w(i,j,k,tim) = 0.9*tem9(i,j,k) + 0.1*tem4(i,j,k)
!       temlz(i,j,k) = w(i,j,k,tim)
!       tem17(i,j,k) = w(i,j,k,tim)       ! only for output
!ccxxx          tem2(i,j,k) = u(i,j,k,tim)
!ccxxx          tem3(i,j,k) = v(i,j,k,tim)
!ccxxx          tem4(i,j,k) = w(i,j,k,tim)
      END DO
    END DO
  END DO

  PRINT*,'base state',ubar(4,4,4),vbar(4,4,4)

!   open(30,file='data.blended',status='unknown',
!  :                     form='unformatted')
!   write(30) tem1
!   write(30) tem12
!   write(30) tem13
!   write(30) temlz
!   close(30)
!
!----------------------------------------------------------------------
!
!  Perform the forward velocity assimilation. Two possible choices
!  are provided here. One is the modified Liou/Gal-Chen procedure,
!  which incorporates the radial wind component from a Doppler radar
!  into a numerical model and the velocity fields are adjusted to
!  satisfy the anelastic mass conservation equation while minimizing
!  the difference between the 'first guess' or 'background' winds and
!  the adjusted winds. Another alternative way is to insert the radial
!  velocity component directly into the model while preserving the
!  model's polar and tangential components.
!
!----------------------------------------------------------------------
!
  CALL assimvel(nx,ny,nz,x,y,z,zp,                                      &
                u,v,w,wcont,ptprt,pprt,qv,qc,qr,qi,qs,qh,               &
                ubar,vbar,ptbar,pbar,rhostr,qvbar,j1,j2,j3,             &
                tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,tem9,           &
                tem10,tem11,tem12,tem13)

  WRITE(6,*) 'code enetr assimvel correctly'

!   open(30,file='data.VarAdj',status='unknown',form='unformatted')
!   write(30) tem1
!   write(30) tem2
!   write(30) tem3
!   write(30) tem4
!   close(30)
!
!----------------------------------------------------------------------
!
!  Apply boundary conditions to u,v,w velocities. If ivar=1, these
!  will be the variationally adjusted velocities.
!
!  NOTE: there is no provision for radiation lateral boundary
!       conditions; if radiation conditions are chosen, zero-gradient
!       conditions are applied instead.
!
!-----------------------------------------------------------------------
!
  IF (ebc == 4) THEN
    newebc = 3
  ELSE
    newebc = ebc
  END IF

  IF (wbc == 4) THEN
    newwbc = 3
  ELSE
    newwbc = wbc
  END IF

  IF (nbc == 4) THEN
    newnbc = 3
  ELSE
    newnbc = nbc
  END IF

  IF (sbc == 4) THEN
    newsbc = 3
  ELSE
    newsbc = sbc
  END IF

  CALL bcu(nx,ny,nz,0., tem2,                                           &
           0.,0.,0.,0.,newebc,newwbc,newnbc,newsbc,tbc,bbc,             &
           newebc, newwbc,newnbc,newsbc)

  CALL bcv(nx,ny,nz,0., tem3,                                           &
           0.,0.,0.,0.,newebc,newwbc,newnbc,newsbc,tbc,bbc,             &
           newebc, newwbc,newnbc,newsbc)

  CALL lbcw(nx,ny,nz,0., tem4, tem5,                                    &
           0.,0.,0.,0.,newebc,newwbc,newnbc,newsbc,                     &
           newebc, newwbc,newnbc,newsbc)

  CALL rhouvw(nx,ny,nz,rhostr,tem5,tem6,tem7)

  CALL wcontra(nx,ny,nz,tem2,tem3,tem4,                                 &
              mapfct,                                                   &
              j1,j2,j3,                                                 &
              rhostr,tem5,tem6,tem7,wcont,tem8,tem9)

  CALL vbcw(nx,ny,nz,tem4,wcont,tbc,bbc,tem2,tem3,                      &
             rhostr,tem5,tem6,tem7,                                     &
             j1,j2,j3)

!   open(30,file='data.AfterBC',status='unknown',form='unformatted')
!   write(30) tem1
!   write(30) tem2
!   write(30) tem3
!   write(30) tem4
!   close(30)

!
!----------------------------------------------------------------------
!
!  If insrt = 1, insert the velocities in place of the model
!  velocities. Due to the leapfrog scheme, overwrite the velocities
!  at tpast with those at tpresent.
!
!----------------------------------------------------------------------
!
  IF (insrt == 1.) THEN

    DO k = 1,nz
      DO j = 1,ny
        DO i = 1,nx
          u(i,j,k,tpresent) = tem2(i,j,k)
          v(i,j,k,tpresent) = tem3(i,j,k)
          w(i,j,k,tpresent) = tem4(i,j,k)
!ccxxx          pprt(i,j,k,tpresent) = 0.0
!ccxxx          ptprt(i,j,k,tpresent)= 0.0
!ccxxx          qv(i,j,k,tpresent) = qvbar(i,j,k)
!ccxxx           qc(i,j,k,tpresent) = 0.0
!ccxxx          qr(i,j,k,tpresent) = 0.0
          u(i,j,k,tpast)    = tem2(i,j,k)
          v(i,j,k,tpast)    = tem3(i,j,k)
          w(i,j,k,tpast)    = tem4(i,j,k)
          pprt(i,j,k,tpast) = pprt(i,j,k,tpresent)
          ptprt(i,j,k,tpast)= ptprt(i,j,k,tpresent)
          qv(i,j,k,tpast)   = qv(i,j,k,tpresent)
          qc(i,j,k,tpast)   = qc(i,j,k,tpresent)
          qr(i,j,k,tpast)   = qr(i,j,k,tpresent)
        END DO
      END DO
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Find a unique name hdmpfn(1:ldmpf) for the history dump data
!  set at time 'curtim'. This file name is stored in 'adjdat' which
!  is used by the thermodynamic recovery.
!
!-----------------------------------------------------------------------
!

  lenstr = LEN(assimnm)
  CALL strlnth (assimnm, lenstr)
  CALL gtdmpfn(assimnm(1:lenstr),dirnam,ldirnm,curtim,                  &
               hdmpfmt,mgrid,nestgrd, hdmpfn, ldmpf)

  WRITE(6,'(1x,a,a)')                                                   &
       'Adjusted data dump in file ',hdmpfn(1:ldmpf)

  adjdat(icount) = hdmpfn(1:ldmpf)

  icount = icount + 1    !!! initial value=1

  grdbas = 0      !!! Note:No base state or grid array is dumped.
!
!-----------------------------------------------------------------------
!
!  Write-out the model-adjusted velocity fields (tem2,tem3,and tem4),
!  pressure,temperature and moisture variables.  This is done so that
!  when a recovery is performed, we can obtain an estimate of the
!  velocity time derivatives. The data dump occurs in conjunction
!  with the variational adjustment, namely at t1,t2,t3......tn
!
!        t1             t2             t3             tn
!       --x--------------x--------------x--------------x--
!       dump           dump           dump           dump
!
!  Compute rhobar (tem9) from rhostr and j3.
!
!-----------------------------------------------------------------------
!
  tim = tpresent

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem9(i,j,k)  = rhostr(i,j,k)/j3(i,j,k)
        tem8(i,j,k)  = 0.0
      END DO
    END DO
  END DO

  ALLOCATE(tem4dsoil(nx,ny,nzsoil,0:nstyps))
  tem4dsoil = 0.0

  CALL dtadump(nx,ny,nz,nzsoil,nstyps,                                  &
       hdmpfmt,nchout,hdmpfn(1:ldmpf),grdbas,filcmprs,                  &
       tem2,tem3,tem4,ptprt(1,1,1,tim),                                 &
       pprt(1,1,1,tim),qv(1,1,1,tim),qc(1,1,1,tim),qr(1,1,1,tim),       &
       qi(1,1,1,tim),qs(1,1,1,tim),qh(1,1,1,tim),tke,kmh,kmv,           &
       ubar,vbar,tem8,ptbar,pbar,tem9,qvbar,                            &
       x,y,z,zp,zpsoil,                                                 &
       tem1(1,1,1),tem1(1,1,1),tem1(1,1,1),                             &
       tem1(1,1,1),tem1(1,1,1),tem1(1,1,1),                             &
       tem4dsoil,tem4dsoil,tem1(1,1,1),tem1(1,1,1),                     &
       tem1(1,1,1),tem1(1,1,1),tem1(1,1,1),                             &
       tem1(1,1,1),tem1(1,1,1),tem1(1,1,1),                             &
       tem1(1,1,1),tem1(1,1,1),                                         &
       tem1(1,1,1),tem1(1,1,1),tem1(1,1,1),tem1(1,1,1),                 &
       tem5,tem6,tem7)

  DEALLOCATE(tem4dsoil)
!
!-----------------------------------------------------------------------
!
!  If recovery flag is turned on (irecov=1), read in the adjusted
!  velocity fields and/or moisture fields at three time levels for
!  thermodynamic retrieval.
!
!-----------------------------------------------------------------------
!

  350   CONTINUE

  IF (irecov == 1) THEN

    isrc = 4
    itag = ii  ! Counter indicating input file name and time

    CALL assimrd(nx,ny,nz,nzsoil,x,y,z,zp,zpsoil,                       &
                 isrc,itag,istat,assimtim,                              &
                 ubar,vbar,pbar,ptbar,rhostr,qvbar,tem11,               &
                 u,v,w,pprt,ptprt,qv,qc,qr,qi,qs,qh,j1,j2,j3,           &
                 tem1,tem2,tem3,tem4,tem5,tem6,                         &
                 tem7,tem8,tem9,tem10)

!
!-----------------------------------------------------------------------
!
!  Set the model time back at the middle time level (t2), and perform
!  a time interpolation on the input data to get data ready at three
!  model time steps t2-dtbig, t2, t2+dtbig.
!
!  NOTE:  All prognostic variables with the exception of qi,qs, and
!         qh are overwritten at all three time levels.
!
!-----------------------------------------------------------------------
!
    t1 = assimtim(ii-3)
    t2 = assimtim(ii-2)
    t3 = assimtim(ii-1)

    curtim = INT((t2+.01)/dtbig)*dtbig
    nstep  = INT( ((curtim-tstart)+0.1*dtbig)/dtbig ) + 1

    IF (ABS(t2-curtim) > ABS(t2-(curtim+dtbig))) THEN
      curtim = curtim + dtbig
    END IF

    CALL retrint(nx,ny,nz,                                              &
                 u,v,w,ptprt,pprt,qv,qc,qr,                             &
                 t1,t2,t3,                                              &
                 tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,tem9)

  END IF

  IF(ivar == 0) THEN
    PRINT *,'Fill wcont after time interpolation when ivar=0'
!
!----------------------------------------------------------------------
!
!  WARNING:  Although u,v,w,qc,qv,qr,pprt,ptprt, are over-written
!           wcont has not been. Upon exiting this routine, FORCE3D
!           calls ADVUVW which uses wcont in the calculation of the
!           vertical advection terms.
!
!  Calculate wcont at time tpresent. wcont will be used in FRCUVW
!  and FRCP. Average rhostr to u, v, w points. They are stored
!  in tem5, tem6 and tem7 respectively.
!
!
!-----------------------------------------------------------------------
!
    tim = tpresent

    CALL rhouvw(nx,ny,nz,rhostr,tem5,tem6,tem7)

    CALL wcontra(nx,ny,nz,u(1,1,1,tim),v(1,1,1,tim),w(1,1,1,tim),       &
                 mapfct,                                                &
                 j1,j2,j3,rhostr,tem5,tem6,tem7,wcont,tem8,tem9)

  END IF

!
!-----------------------------------------------------------------------
!
!  These temporary arrays are returned to FORCE3D via the arrays
!  uforce, vforce, wforce and defsq.  In order to avoid contaminating
!  any of these terms we set them to zero prior to returning.
!
!-----------------------------------------------------------------------
!
  PRINT *, 'In ASSIMDRIV fill tem10,tem11,tem12,tem13 with zeros'
  DO k=1,nz
    DO j=1,ny
      DO i=1,nx
        tem10(i,j,k) =  0.0
        tem11(i,j,k) =  0.0
        tem12(i,j,k) =  0.0
        tem13(i,j,k) =  0.0
        tem14(i,j,k) =  0.0
        tem15(i,j,k) =  0.0
        tem16(i,j,k) =  0.0
      END DO
    END DO
  END DO

  RETURN

END SUBROUTINE assimdriv
