PROGRAM radmosaic
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  PROGRAM RADMOSAIC                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!  PURPOSE:
!
!  Create a 3-d radar mosaic from remapped radar files
!  (output from 88d2arps and/or nids2arps).
!
!  AUTHOR:
!
!  Keith Brewster, CAPS, August, 2003
!
!  MODIFICATION HISTORY:
!  2012/05/03 Keith Brewster, CAPS
!  Added grid NAMELIST to obtain complete grid specifications.
!  Requires corresponding change to input file.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: nx       ! Number of grid points in the x-direction
  INTEGER :: ny       ! Number of grid points in the y-direction
  INTEGER :: nz       ! Number of grid points in the z-direction
  INTEGER :: nzsoil   ! Number of grid points in the -z-direction
!
!-----------------------------------------------------------------------
!
! Maximum number of radars in the mosaic
!
!-----------------------------------------------------------------------
!
  INTEGER, PARAMETER :: mx_rad=120
!
  REAL, ALLOCATABLE :: x(:)
  REAL, ALLOCATABLE :: y(:)
  REAL, ALLOCATABLE :: z(:)
  REAL, ALLOCATABLE :: zp(:,:,:)
  REAL, ALLOCATABLE :: zpsoil(:,:,:)
  REAL, ALLOCATABLE :: hterain(:,:)
  REAL, ALLOCATABLE :: mapfct(:,:,:)
  REAL, ALLOCATABLE :: j1(:,:,:)
  REAL, ALLOCATABLE :: j2(:,:,:)
  REAL, ALLOCATABLE :: j3(:,:,:)
  REAL, ALLOCATABLE :: j3soil(:,:,:)
  REAL, ALLOCATABLE :: j3inv(:,:,:)
  REAL, ALLOCATABLE :: j3soilinv(:,:,:)
!
  REAL, ALLOCATABLE :: xs(:)
  REAL, ALLOCATABLE :: ys(:)
  REAL, ALLOCATABLE :: zs(:,:,:)
  REAL, ALLOCATABLE :: refmos(:,:,:)
  REAL, ALLOCATABLE :: refl(:,:,:)
  REAL, ALLOCATABLE :: tem1(:,:,:)
  REAL, ALLOCATABLE :: tem2(:,:,:)
  REAL, ALLOCATABLE :: tem3(:,:,:)
  CHARACTER (LEN=256) :: radfname(mx_rad)
  CHARACTER (LEN=6)   :: varid
  CHARACTER (LEN=20)  :: varname
  CHARACTER (LEN=20)  :: varunits
  INTEGER :: nradfil
  INTEGER :: istatus
  INTEGER :: i,j,k,irad
  REAL :: alatru(2)
  REAL :: ctrx,ctry,swx,swy,dxscl,dyscl
  REAL :: rhinf
  REAL :: rvinf
  REAL :: refmax

  INTEGER :: mosaicfmt
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'mp.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'globcst.inc'
!
  NAMELIST /grid_dims/ nx,ny,nz
  NAMELIST /terrain/ ternopt,mntopt,hmount,mntwidx,mntwidy,             &
            mntctrx,mntctry,terndta,ternfmt
  NAMELIST /grid/ dx,dy,dz,strhopt,dzmin,zrefsfc,dlayer1,dlayer2,       &
            strhtune,zflat,ctrlat,ctrlon, crdorgnopt
  NAMELIST /filelist/ nradfil,radfname
  NAMELIST /output/ runname,dirname,mosaicfmt,lvldbg
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
  CALL mpinit_var

!-----------------------------------------------------------------------
!
! Initializations
!
!-----------------------------------------------------------------------
!
  nx=5
  ny=5
  nz=5
  ternopt=1
  mntopt=0
  hmount=0.
  mntwidx=0.
  mntwidy=0.
  mntctrx=0.
  mntctry=0.
  terndta='NULL'
  ternfmt=1
  nradfil=1
  DO irad=1,mx_rad
    radfname(irad)='NULL'
  END DO
  runname='NULL'
  dirname='NULL'
  mosaicfmt = 3
  lvldbg=0
!
  READ(5,grid_dims)
  READ(5,terrain)
  READ(5,grid)
  READ(5,filelist)
  READ(5,output)
!
  CALL get_grid_from_rad(radfname(1),dx,dy,dz,dzmin,strhopt,           &
                  mapproj,ctrlat,ctrlon,trulat1,trulat2,trulon,sclfct, &
                  istatus)

  print *, ' dz,dzmin= ',dz,dzmin
  print *, ' strhiopt: ',strhopt
  
  alatru(1)=trulat1
  alatru(2)=trulat2
  nzsoil=2
  dzsoil=1.0
!
  ALLOCATE(x(nx),STAT=istatus)
  CALL check_alloc_status(istatus, "radmosaic:x")
  ALLOCATE(y(ny),STAT=istatus)
  CALL check_alloc_status(istatus, "radmosaic:y")
  ALLOCATE(z(nz),STAT=istatus)
  CALL check_alloc_status(istatus, "radmosaic:z")
  ALLOCATE(zp(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "radmosaic:zp")
  ALLOCATE(zpsoil(nx,ny,nzsoil),STAT=istatus)
  CALL check_alloc_status(istatus, "radmosaic:zpsoil")
  ALLOCATE(hterain(nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "radmosaic:hterain")
  ALLOCATE(mapfct(nx,ny,8),STAT=istatus)
  CALL check_alloc_status(istatus, "radmosaic:mapfct")
  ALLOCATE(j1(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "radmosaic:j1")
  ALLOCATE(j2(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "radmosaic:j2")
  ALLOCATE(j3(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "radmosaic:j3")
  ALLOCATE(j3soil(nx,ny,nzsoil),STAT=istatus)
  CALL check_alloc_status(istatus, "radmosaic:j3soil")
  ALLOCATE(j3inv(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "radmosaic:j3inv")
  ALLOCATE(j3soilinv(nx,ny,nzsoil),STAT=istatus)
  CALL check_alloc_status(istatus, "radmosaic:j3invsoil")
!
  ALLOCATE(xs(nx),STAT=istatus)
  CALL check_alloc_status(istatus, "radmosaic:xs")
  ALLOCATE(ys(ny),STAT=istatus)
  CALL check_alloc_status(istatus, "radmosaic:ys")
  ALLOCATE(zs(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "radmosaic:zs")
  ALLOCATE(refmos(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "radmosaic:refmos")
  ALLOCATE(refl(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "radmosaic:refl")
  ALLOCATE(tem1(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "radmosaic:tem1")
  ALLOCATE(tem2(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "radmosaic:tem2")
  ALLOCATE(tem3(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "radmosaic:tem3")
!
!-----------------------------------------------------------------------
!
!  Initialization of model grid definition arrays.
!
!-----------------------------------------------------------------------
!

  CALL inigrd(nx,ny,nz,nzsoil,x,y,z,zp,zpsoil,                           &
                hterain,mapfct,j1,j2,j3,j3soil,                          &
                j3soilinv,tem1,tem2,tem3)

  CALL setmapr( mapproj,sclfct,alatru,trulon )

  CALL lltoxy( 1,1, ctrlat,ctrlon, ctrx, ctry )

  dxscl = dx*sclfct
  dyscl = dy*sclfct

  swx = ctrx - (0.5*FLOAT(nx-3)) * dxscl
  swy = ctry - (0.5*FLOAT(ny-3)) * dyscl

  CALL setorig( 1, swx, swy)
!
!-----------------------------------------------------------------------
!
! Fill xs, ys and zs arrays.
!
!-----------------------------------------------------------------------
!
  DO i = 1, nx-1
    xs(i)=0.5*x(i)+x(i+1)
  END DO
  xs(nx)=2.*x(nx-1)-x(nx-2)
  DO j = 1, ny-1
    ys(j)=0.5*y(j)+y(j+1)
  END DO
  ys(ny)=2.*y(ny-1)-y(ny-2)
!
  DO k=1, nz-1
    DO j=1, ny
      DO i=1, nx
        zs(i,j,k)=0.5*(zp(i,j,k)+zp(i,j,k+1))
      END DO
    END DO
  END DO
  DO j=1, ny
    DO i=1, nx
      zs(i,j,nz)=2.*zp(i,j,nz-1)-zp(i,j,nz-2)
    END DO
  END DO
  rhinf=2.*dz
  rvinf=2.*sqrt(0.5*(dx*dx + dy*dy))
!
  CALL refmosaic(nradfil,nx,ny,nz,mx_rad,                               &
           xs,ys,zs,radfname,lvldbg,refmos,rhinf,rvinf,                 &
           refl,tem1,tem2,istatus)
!
  DO j=1,ny
    DO i=1,nx
      refmax=-99.
      DO k=2,nz-1
        refmax=max(refmax,refmos(i,j,k))
      END DO
      refmos(i,j,1)=refmax
    END DO
  END DO
!
  varid='refmos'
  varname='Reflectivity Mosaic'
  varunits='dBZ'
  curtim=0.
  CALL wrtvar2(nx,ny,nz,refmos,varid,varname,varunits,curtim,runname,   &
               dirname,mosaicfmt,0,0,istatus)
  STOP
END PROGRAM radmosaic
