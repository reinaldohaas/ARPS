!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE INIGRD                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE inigrd(nx,ny,nz,nzsoil,x,y,z,zp,zpsoil,                      &
                  hterain,mapfct,j1,j2,j3,j3soil,j3soilinv,             &
                  zp1d,dzp1d,tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Initialize the model grid variables.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  3/17/1991.
!
!  MODIFICATION HISTORY:
!
!  5/02/92 (M. Xue)
!  Added full documentation.
!
!  5/03/92 (M. Xue)
!  Further documentation.
!
!  4/10/93 (M. Xue & Hao Jin)
!  Add the terrain.
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!  12/22/94 (Yuhe Liu)
!  Changed variable tuning, which was hard wired inside this
!  subroutine, to strhtune, which is input from namelist input file.
!
!  1/22/96 (Donghai Wang & Yuhe Liu)
!  Added the map projection factor to ARPS governing equations.
!
!  11/06/97 (D. Weber)
!  Added three additional levels to the mapfct array.  The three
!  levels (4,5,6) represent the inverse of the first three in order.
!  The inverse map factors are computed to improve efficiency.
!
!  2000/04/11 (Gene Bassett)
!  Only update the terrain data with a call to BCS2D for ternopt=1.
!
!  05/02/2002 (Dan Weber and Jerry Brotzge)
!  Added zpsoil, j3soil, variables for OUSOIL.
!
!  12/7/2010 (Bryan Putnam)
!  Added option to remove hail from initial data file and add it to the
!  graupel for option to turn off hail.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nzsoil   Number of grid points in the soil model in the -z-direction
!
!  OUTPUT:
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space
!             (m)
!    zpsoil   Vertical coordinate of grid points in the soil model
!             in physical space (m).
!    hterain  Terrain height (m)
!    mapfct   Map factors at scalar, u and v points
!
!    j1       Coordinate transformation Jacobian -d(zp)/dx
!    j2       Coordinate transformation Jacobian -d(zp)/dy
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!    j3soil   Soil coordinate transformation Jacobian  d(zpsoil)/dz
!    j3soilinv Inverse of the soil coordinate transformation j3soil
!
!  WORK ARRAY:
!
!    zp1d     Temporary work array.
!    dzp1d    Temporary work array.
!    tem1     Temporary work array.
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

  INTEGER :: nx,ny,nz     ! The number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of grid points in the -z-direction

  REAL :: x(nx)           ! x-coord. of the physical and
                          ! computational grid. Defined at u-point.
  REAL :: y(ny)           ! y-coord. of the physical and
                          ! computational grid. Defined at v-point.
  REAL :: z(nz)           ! z-coord. of the computational grid.
                          ! Defined at w-point on the staggered grid.
  REAL :: zsoil(nzsoil)   ! z-coord. of the soil model computational grid.
                          ! Defined at the center of a soil layer.
  REAL :: zp(nx,ny,nz)    ! Physical height coordinate defined at
                          ! w-point of the staggered grid.
  REAL :: zpsoil (nx,ny,nzsoil) ! The physical height coordinate defined
                          ! at the center of a soil layer(m).

  REAL :: j1(nx,ny,nz)    ! Coordinate transformation Jacobian
                          ! -d(zp)/dx.
  REAL :: j2(nx,ny,nz)    ! Coordinate transformation Jacobian
                          ! -d(zp)/dy.
  REAL :: j3(nx,ny,nz)    ! Coordinate transformation Jacobian
                          ! d(zp)/dz.
  REAL :: j3soil (nx,ny,nzsoil) ! Coordinate transformation Jacobian
                          ! defined as d( zpsoil )/d( zsoil ).
  REAL :: j3soilinv(nx,ny,nzsoil) ! Inverse of J3soil.


  REAL :: hterain(nx,ny)  ! Terrain height.

  REAL :: mapfct(nx,ny,8)      ! Map factors at scalar, u and v points

  REAL :: zp1d (nz)       ! Temporary array
  REAL :: dzp1d(nz)       ! Temporary array
  REAL :: tem1(nx,ny,nz)  ! Temporary array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  REAL :: xs, ys, zs, radnd, pi2
  REAL :: zflat1,z1,z2

  REAL :: zpmax
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'phycst.inc'
  INCLUDE 'bndry.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
!-----------------------------------------------------------------------
!
!  Define a uniform model grid.
!
!-----------------------------------------------------------------------
!
  CALL setgrd ( nx,ny, x,y )
!
!-----------------------------------------------------------------------
!
!  Set map factor at scalar, u, and v points
!
!-----------------------------------------------------------------------
!
  IF ( mpfctopt /= 0 ) THEN

    DO j=1,ny-1
      DO i=1,nx-1
        xs = 0.5*(x(i)+x(i+1))
        ys = 0.5*(y(j)+y(j+1))
        CALL xytomf( 1,1,xs,ys,mapfct(i,j,1) )
        IF(maptest == 1)THEN   ! set up a symmetric field...
          mapfct(i,j,1) = 1.0 + ABS((xs-0.5*(x(1)+x(nx)))               &
                                   /(x(nx)-x(1))                        &
                                   +(ys-0.5*(y(1)+y(ny)))               &
                                   /(y(ny)-y(1)))
        END IF
        mapfct(i,j,4) = 1.0/mapfct(i,j,1)
        mapfct(i,j,7) = mapfct(i,j,1)*mapfct(i,j,1)
        mapfct(i,j,8) = 0.25*mapfct(i,j,1)  ! for use in sovlwpim
                                            ! and wcontra...
      END DO
    END DO

    DO j=1,ny-1
      DO i=1,nx
        ys = 0.5*(y(j)+y(j+1))
        CALL xytomf( 1,1,x(i),ys,mapfct(i,j,2) )
        IF(maptest == 1)THEN   ! set up a symmetric field...
          mapfct(i,j,2) = 1.0 + ABS((x(i)-0.5*(x(1)+x(nx)))             &
                                   /(x(nx)-x(1))                        &
                                   +(ys-0.5*(y(1)+y(ny)))               &
                                   /(y(ny)-y(1)))
        END IF
        mapfct(i,j,5) = 1.0/mapfct(i,j,2)
      END DO
    END DO

    DO j=1,ny
      DO i=1,nx-1
        xs = 0.5*(x(i)+x(i+1))
        CALL xytomf( 1,1,xs,y(j),mapfct(i,j,3) )
        IF(maptest == 1)THEN   ! set up a symmetric field...
          mapfct(i,j,3) = 1.0 + ABS((xs-0.5*(x(1)+x(nx)))               &
                                   /(x(nx)-x(1))                        &
                                   +(y(j)-0.5*(y(1)+y(ny)))             &
                                   /(y(ny)-y(1)))
        END IF
        mapfct(i,j,6) = 1.0/mapfct(i,j,3)
      END DO
    END DO

  ELSE

    DO k=1,7
      DO j=1,ny
        DO i=1,nx
          mapfct(i,j,k) = 1.0
          mapfct(i,j,8) = 0.25
        END DO
      END DO
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Print map factor at scalar, u, and v points for the purpose of
!  debug
!
!-----------------------------------------------------------------------
!
!   CALL getunit( mpunit )
!   CALL asnctl ('NEWLOCAL', 1, ierr)
!   CALL asnfile('mapfactor.data', '-F f77 -N ieee', ierr)
!
!   OPEN (mpunit,file='mapfactor.data',form='unformatted')
!
!   CALL edgfill(mapfct(1,1,1),1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
!   write(mpunit) ((mapfct(i,j,1),i=1,nx),j=1,ny)
!
!   CALL edgfill(mapfct(1,1,2),1,nx,1,nx,   1,ny,1,ny-1, 1,1,1,1)
!   write(mpunit) ((mapfct(i,j,2),i=1,nx),j=1,ny)
!
!   CALL edgfill(mapfct(1,1,3),1,nx,1,nx-1, 1,ny,1,ny,   1,1,1,1)
!   write(mpunit) ((mapfct(i,j,3),i=1,nx),j=1,ny)
!
!   CLOSE( mpunit )
!   CALL retunit( mpunit )
!
!-----------------------------------------------------------------------
!
!  End of debug print.
!
!-----------------------------------------------------------------------
!
  DO k=1,nz
    z(k) = zrefsfc + (k-2) * dz
  END DO
!
!-----------------------------------------------------------------------
!
!  Specify the terrain
!
!-----------------------------------------------------------------------
!
  IF( ternopt == 0 ) THEN     ! No terrain, the ground is flat

    DO j=1,ny-1
      DO i=1,nx-1
        hterain(i,j) = zrefsfc
      END DO
    END DO

  ELSE IF( ternopt == 1 ) THEN  ! Bell-shaped mountain

    IF (mntopt == 1) THEN
!
!-----------------------------------------------------------------------
!
!  Define the bell-shaped mountain
!
!-----------------------------------------------------------------------
!
    pi2 = 1.5707963267949

    DO j=1,ny-1
      DO i=1,nx-1
        xs = (x(i)+x(i+1))*0.5
        ys = (y(j)+y(j+1))*0.5
        zs = z(2)

        IF( mntwidy < 0.0 .OR. runmod == 2 ) THEN ! 2-d terrain in x-z plane.

          radnd = 1.0+((xs-mntctrx)/mntwidx)**2

        ELSE IF( mntwidx < 0.0 .OR. runmod == 3 ) THEN ! 2-d terrain in y-z plane.

          radnd = 1.0+((ys-mntctry)/mntwidy)**2

        ELSE                             ! 3-d terrain

          radnd = 1.0+((xs-mntctrx)/mntwidx)**2                         &
                 +((ys-mntctry)/mntwidy)**2
        END IF

        hterain(i,j) = zrefsfc + hmount/radnd

      END DO
    END DO

    ELSE  ! mntopt == 2 or others

      WRITE(6,'(1x,a,/)') ' Please supply your own terrain function in subroutine inigrd'
      CALL arpsstop('',1)
    END IF
!
!-----------------------------------------------------------------------
!
!  Make sure that the terrain satisfies the specified boundary
!  conditions.
!
!-----------------------------------------------------------------------
!
    IF (mp_opt > 0) THEN
      CALL acct_interrupt(mp_acct)
      CALL mpsendrecv1dew(hterain,nx,ny,ebc,wbc,0,tem1)
      CALL mpsendrecv1dns(hterain,nx,ny,nbc,sbc,0,tem1)
    END IF
    CALL acct_interrupt(bc_acct)
    CALL bcs2d(nx,ny,hterain, ebc,wbc,nbc,sbc)
    CALL acct_stop_inter

  ELSE IF( ternopt == 2 ) THEN          ! Read from terrain data base

!
!-----------------------------------------------------------------------
!
!    Read in the terrain data.
!
!-----------------------------------------------------------------------
!
!  blocking inserted for ordering i/o for message passing
    DO i=0,nprocs-1,readstride
      IF(myproc >= i.AND.myproc <= i+readstride-1)THEN

        IF (mp_opt > 0 .AND. readsplit(FINDX_T) > 0) THEN

        CALL readsplittrn( nx,ny,dx,dy, terndta,                        &
               mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,     &
               hterain )

        ELSE

        CALL readtrn( nx,ny,dx,dy, terndta,                             &
               mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,     &
               hterain )

        END IF

      END IF

      IF (mp_opt > 0) CALL mpbarrier
    END DO

  END IF

!
!-----------------------------------------------------------------------
!
!  Set up a stretched vertical grid.
!
!  For strhopt=1, function y = a+b*x**3 is used to specify dz as a
!                              function of k.
!  For strhopt=2, function y = c + a*tanh(b*x) is used to specify dz
!                              as a function of k.
!
!-----------------------------------------------------------------------
!
  IF ( strhopt == 0 ) THEN

    DO k=1,nz
      zp1d(k) = z(k)
    END DO

  ELSE IF ( strhopt == 1 .OR.strhopt == 2 ) THEN

    z1 = zrefsfc + MAX(0.0, MIN(dlayer1, z(nz-2)-zrefsfc ))
    z2 = z1      + MAX(0.0, MIN(dlayer2, z(nz-1)-z1      ))

    IF( dlayer1 >= (nz-3)*dzmin ) THEN
      WRITE(6,'(/1x,a,f13.3,/a,f13.3,a,a)')                             &
          'Can not setup a vertical grid with uniform dz=',dzmin,       &
          ' over the depth of ',dlayer1,' please specify a smaller ',   &
          'value of dlayer1. Program stopped INIGRD.'
      CALL arpsstop('arpsstop called from INIGRD with ther vertical grid ' &
            ,1)
    END IF

    CALL strhgrd(nz,strhopt,zrefsfc,z1,z2,z(nz-1),                      &
                 dzmin,strhtune, zp1d,dzp1d)

  ELSE

    WRITE(6,'(1x,a,i3,a/)')                                             &
        'Invalid vertical grid stretching option, strhopt was ',strhopt, &
        '. Program stopped in INIGRD.'
      CALL arpsstop('arpsstop called from INIGRD with stretching ' ,1)

  END IF
!
!-----------------------------------------------------------------------
!
!  Physical height of computational grid defined as
!
!  Zp=(z-zrefsfc)*(Zm-hterain)/(Zm-zrefsfc)+hterain for z=<Zm.
!  ZP=z for z>Zm
!
!  where Zm the height at which the grid levels becomes flat.
!  Hm < Zm =< Ztop, hm is the height of mountain and Ztop the height
!  of model top.
!
!-----------------------------------------------------------------------
!
  DO k=nz-1,2,-1
    IF(zp1d(k) <= zflat) THEN
      zflat1 = zp1d(k)
      EXIT
    END IF
  END DO
  zflat1=MAX(MIN(z(nz-1),zflat1),zrefsfc)

  DO k=2,nz-1

    IF(zp1d(k) > zflat1) THEN
      DO j=1,ny-1
        DO i=1,nx-1
          zp(i,j,k)=zp1d(k)
        END DO
      END DO
    ELSE
      DO j=1,ny-1
        DO i=1,nx-1
          zp(i,j,k)=(zp1d(k)-zrefsfc)*(zflat1-hterain(i,j))             &
                     /(zflat1-zrefsfc)+hterain(i,j)
        END DO
      END DO
    END IF

  END DO

  DO j=1,ny-1
    DO i=1,nx-1
      zp(i,j,2)=hterain(i,j)
      zp(i,j,1)=2.0*zp(i,j,2)-zp(i,j,3)
      zp(i,j,nz)=2.0*zp(i,j,nz-1)-zp(i,j,nz-2)
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Calculate transformation Jacobians J1, J2 and J3.
!
!-----------------------------------------------------------------------
!
  CALL jacob(nx,ny,nz,x,y,z,zp,j1,j2,j3,tem1)

  CALL inisoilgrd(nx,ny,nzsoil,hterain,zpsoil,j3soil,j3soilinv)


  RETURN
END SUBROUTINE inigrd

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE STRHGRD                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE strhgrd(nz,strhopt,z0,z1,z2,ztop,dzmin,strhtune, z,dzk)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  To construct a vertically stretched grid.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/17/94.
!
!  MODIFICATION HISTORY:
!
!  05/11/95 (Jinxing Zong and MX)
!
!  A bug fix for the case of nonzero zrefsfc, the reference height of
!  the surface. Results not affected for zrefsfc=0 (default value).
!
!-----------------------------------------------------------------------
!
!
!  INPUT:
!
!    nz       The vertical dimension of ARPS grid.
!
!
!  OUTPUT:
!
!    z        Array containing the height of veritical grid levels.
!    dzk      Array containing the grid spacing between vertical levels
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: nz
  INTEGER :: strhopt
  REAL :: z0
  REAL :: z1
  REAL :: z2
  REAL :: ztop
  REAL :: dzmin
  REAL :: strhtune

  REAL :: z   (nz)
  REAL :: dzk (nz)

  REAL :: rnzh,dzm
  REAL :: a,b,c,hnew,zkmid,dzu
  INTEGER :: nzh,nzl,k
  REAL :: dz

  IF( (z1-z0) == (nz-3)*dzmin.AND.(z1-z0) == (ztop-z0) ) THEN

    dz = (ztop-z0)/(nz-3)
    DO k=1,nz-1
      dzk(k)= dz
    END DO
    DO k=1,nz
      z(k)=z0 + (k-2) * dz
    END DO

    WRITE(6,'(/1x,a,f13.3,/a,f13.3)')                                   &
        'Layer 1 depth was as deep as the entire domain. i',            &
        'A uniform vertical grid is assumed with dz=',dz,               &
        ' over the model depth of ',ztop-z0

    RETURN

  END IF

  IF(z1 < z0) z1 = z0
  IF(z2 > ztop) z2 = ztop

  nzl = (z1-z0)/dzmin

  IF( (z1-z0) >= (nz-3)*dzmin ) THEN
    WRITE(6,'(/1x,a,f13.3,/a,f13.3,a,a)')                               &
        'Can not setup a vertical grid with uniform dz=',dzmin,         &
        ' over the depth of ',z1-z0,' please specify a smaller',        &
        ' value of dlayer1 '
      CALL arpsstop('arpsstop called from STRHGRD with stretching ' ,1)
  END IF

  IF( z2 >= ztop ) THEN
    dzm = (ztop-z0-nzl*dzmin)/(nz-3-nzl)
!    print*, nzl*dzmin + (nz-3-nzl)*dzm
    nzh = 0
    dzu = 2*dzm - dzmin
  ELSE
    a = 2*(nz-3-nzl)
    b = 2*z0-ztop-z2-(nz-3-3*nzl)*dzmin
    c = dzmin*(z2-z0-nzl*dzmin)
    dzm = (-b + SQRT(b*b-4*a*c) )/(2*a)

    rnzh = (ztop-z2)/(2*dzm-dzmin)
    nzh = INT(rnzh)

    hnew = nzl*dzmin + nzh*(2*dzm-dzmin) +                              &
          (nz-3-nzl-nzh)*dzm + z0

    IF( nzh /= 0 ) THEN
      dzu = (2*dzm-dzmin) + (ztop-hnew)/nzh
    ELSE IF( nz-3-nzl-nzh /= 0 ) THEN
      dzm = dzm + (ztop-hnew)/(nz-3-nzl-nzh)
      dzu = (2*dzm-dzmin)
    END IF

  END IF

  DO k=1,nzl+1
    dzk(k)=dzmin
  END DO


  IF( strhopt == 1 ) THEN

    a   = (dzm-dzmin)
    DO k=nzl+2,nz-2-nzh
      dzk(k)= dzm+a*                                                    &
          ((2.0*REAL(k-nzl-2)/REAL(nz-4-nzh-nzl)-1.0) )**3
    END DO


  ELSE

    zkmid=0.5*REAL( nz-nzh+nzl)

    IF( nzl+2-zkmid == 0.0 ) THEN
      b = 0.0
    ELSE
      b= strhtune* 2.0/(nzl+2-zkmid)
    END IF

    a=(dzmin-dzm)/TANH( strhtune* 2.0)

    DO k=nzl+2,nz-2-nzh
      dzk(k)=dzm + a*TANH(b*(FLOAT(k)-zkmid))
    END DO

  END IF

  DO k=nz-2-nzh+1, nz-2
    dzk(k)= dzu
  END DO

  dzk(nz-1) = dzk(nz-2)
  dzk(nz  ) = dzk(nz-1)


  z(2) = z0
  DO k=3,nz-1
    z(k) = z(k-1)+dzk(k-1)
  END DO

  z(1) =z(2)-dzk(1)
  z(nz-1)=ztop
  z(nz)=z(nz-1)+dzk(nz-1)

  RETURN
END SUBROUTINE strhgrd

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE INISOILGRD                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE inisoilgrd(nx,ny,nzsoil,hterain,zpsoil,j3soil,j3soilinv)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  To construct the soil model grid and all associated variables.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Weber
!  05/24/2002.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nzsoil   Number of grid points in the soil
!
!    zp       Vertical coordinate of grid points in physical space (m)
!
!  OUTPUT:
!
!    zpsoil   Vertical coordinate of grid points in the soil (m)
!    j3soil   Coordinate transformation Jacobian  d(zpsoil)/dz
!    j3soilinv Inverse of the coordinate transformation Jacobian  d(zpsoil)/dz
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!
!  INPUT
!

  INTEGER :: nx       ! Number of grid points in the x-direction
  INTEGER :: ny       ! Number of grid points in the y-direction
  INTEGER :: nzsoil   ! Number of grid points in the -z-direction

  REAL :: hterain(nx,ny) ! The terrain height in meters.

!
!  OUTPUT
!

  REAL :: zpsoil (nx,ny,nzsoil) ! The physical height coordinate defined at
                                ! w-point of the soil.
  REAL :: j3soil(nx,ny,nzsoil)  ! Coordinate transformation Jacobian d(zp)/d(z)
  REAL :: j3soilinv (nx,ny,nzsoil) ! Inverse of Soil coord. trans (1.0/j3soil).

!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!

  include 'grid.inc'
  INCLUDE 'mp.inc'

!
!-----------------------------------------------------------------------
!
!  Miscellaneous local variables:
!
!-----------------------------------------------------------------------
!

  INTEGER :: i,j,k

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
!  Calculate the soil model grid variables:
!
!                  zsoil,zpsoil,j3soil,j3soilinv.
!
!-----------------------------------------------------------------------
!
!  Following ARPS convention,
!  Set zsoil at velocity points or at the top of each soil layer.....
!

   IF(soilstrhopt == 0)THEN   !   set up soil model without a stretched grid..

     DO k = 1,nzsoil
       DO j = 1,ny
         DO i = 1,nx
           zpsoil(i,j,k) = hterain(i,j) - (k-1)*dzsoil
           j3soil(i,j,k)    = 1.0
           j3soilinv(i,j,k) = 1.0
         END DO
       END DO
     END DO

   ELSE IF(soilstrhopt == 1) THEN

!  DO k = 2,nzsoil
!    zsoil(k) = zrefsfc - (k-1)*dzsoil
!  END DO

!-----------------------------------------------------------------------
!
!  Physical height of soil model computational grid defined as
!
!  Zpsoil=(zsoil-zrefsfc)*(Zm-hterain)/(Zm-zrefsfc)+hterain for zsoil=<Zm.
!  ZPsoil=zsoil for zsoil>Zm
!
!  where Zm the height at which the grid levels becomes flat.
!  Hm < Zm =< Ztop, hm is the height of mountain and Ztop the height
!  of model top.
!
!-----------------------------------------------------------------------
!
! DO k=nzsoil-1,2,-1
!   IF(zp1dsoil(k) <= zflatsoil) THEN
!     zflat1soil = zp1dsoil(k)
!     EXIT
!   END IF
! END DO
! zflat1soil=MAX(MIN(zsoil(nz-1),zflat1soil),zrefsfc)

! DO k=2,nzsoil-1
!   IF(zp1dsoil(k) > zflat1soil) THEN
!     DO j=1,ny-1
!       DO i=1,nx-1
!         zpsoil(i,j,k)=zp1dsoil(k)
!       END DO
!     END DO
!   ELSE
!     DO j=1,ny-1
!       DO i=1,nx-1
!         zpsoil(i,j,k)=(zp1dsoil(k)-zrefsfc)*(zflat1soil-hterain(i,j))   &
!                    /(zflat1soil-zrefsfc)+hterain(i,j)
!       END DO
!     END DO
!   END IF

! END DO

! old code..
!  DO k = 2,nzsoil
!    zsoil(k) = zrefsfc - (k-2)*dzsoil
!    zpsoil(i,j,k) = hterain(i,j) - (k-1)*dzsoil
!  END DO

!  jerry's code   note zpsoil is set to the top of each layer.
!  question why are we using positive values?  need to use zrefsfc... and
!  the existing coordinate system...

! zpsoil (1)     = 0.0
! DO k=2,nzsoil                         ! set zpsoil
!!        zpsoil(k) = zpsoil(k-1) - dzsoil * k
!        zpsoil(k) = zpsoil(k-1) + dzsoil    !(positive values)
! END DO                                ! done setting zpsoil

! DO k=1,nzsoil
!   j3soil(k) = 1.0
! END DO
! arps code starts here.

! IF( (z1-z0) == (nz-3)*dzmin.AND.(z1-z0) == (ztop-z0) ) THEN

 !! dz = (ztop-z0)/(nz-3)
!   DO k=1,nz-1
!     dzk(k)= dz
!   END DO
!   DO k=1,nz
!     z(k)=z0 + (k-2) * dz
!   END DO

!   WRITE(6,'(/1x,a,f13.3,/a,f13.3)')                                   &
!       'Layer 1 depth was as deep as the entire domain. i',            &
!       'A uniform vertical grid is assumed with dz=',dz,               &
!       ' over the model depth of ',ztop-z0

!   RETURN

! END IF

! IF(z1 < z0) z1 = z0
! IF(z2 > ztop) z2 = ztop

! nzl = (z1-z0)/dzmin

! IF( (z1-z0) >= (nz-3)*dzmin ) THEN
!   WRITE(6,'(/1x,a,f13.3,/a,f13.3,a,a)')                               &
!       'Can not setup a vertical grid with uniform dz=',dzmin,         &
!       ' over the depth of ',z1-z0,' please specify a smaller',        &
!       ' value of dlayer1 '
!     CALL arpsstop('arpsstop called from STRHGRD with stretching ' ,1)
! END IF

! IF( z2 >= ztop ) THEN
!   dzm = (ztop-z0-nzl*dzmin)/(nz-3-nzl)
!    print*, nzl*dzmin + (nz-3-nzl)*dzm
!   nzh = 0
!   dzu = 2*dzm - dzmin
! ELSE
!   a = 2*(nz-3-nzl)
!   b = 2*z0-ztop-z2-(nz-3-3*nzl)*dzmin
!   c = dzmin*(z2-z0-nzl*dzmin)
!   dzm = (-b + SQRT(b*b-4*a*c) )/(2*a)

!   rnzh = (ztop-z2)/(2*dzm-dzmin)
!   nzh = INT(rnzh)

!   hnew = nzl*dzmin + nzh*(2*dzm-dzmin) +                              &
!         (nz-3-nzl-nzh)*dzm + z0

!   IF( nzh /= 0 ) THEN
!     dzu = (2*dzm-dzmin) + (ztop-hnew)/nzh
!   ELSE IF( nz-3-nzl-nzh /= 0 ) THEN
!     dzm = dzm + (ztop-hnew)/(nz-3-nzl-nzh)
!     dzu = (2*dzm-dzmin)
!   END IF

! END IF

! DO k=1,nzl+1
!   dzk(k)=dzmin
! END DO


! IF( strhopt == 1 ) THEN

!   a   = (dzm-dzmin)
!   DO k=nzl+2,nz-2-nzh
!     dzk(k)= dzm+a*                                                    &
!         ((2.0*FLOAT(k-nzl-2)/FLOAT(nz-4-nzh-nzl)-1.0) )**3
!   END DO


! ELSE

!   zkmid=0.5*FLOAT( nz-nzh+nzl)

!   IF( nzl+2-zkmid == 0.0 ) THEN
!     b = 0.0
!   ELSE
!     b= strhtune* 2.0/(nzl+2-zkmid)
!   END IF

!   a=(dzmin-dzm)/TANH( strhtune* 2.0)

!   DO k=nzl+2,nz-2-nzh
!     dzk(k)=dzm + a*TANH(b*(FLOAT(k)-zkmid))
!   END DO

! END IF

! DO k=nz-2-nzh+1, nz-2
!   dzk(k)= dzu
! END DO

! dzk(nz-1) = dzk(nz-2)
! dzk(nz  ) = dzk(nz-1)

! z(2) = z0
! DO k=3,nz-1
!   z(k) = z(k-1)+dzk(k-1)
! END DO

! z(1) =z(2)-dzk(1)
! z(nz-1)=ztop
! z(nz)=z(nz-1)+dzk(nz-1)

   ELSE IF(soilstrhopt == 9)THEN ! KFY Jan 2008: special for NARR, SREF or NAM
     ! to WRF-ARW.
     ! This option is not meant to run ARPS itself, it's only
     ! for producing initial tsoil and qsoil interpolation on
     ! multiple soil layers (set nzsoil = 6 for NARR and NAM data, and
     ! assign zpsoil the same values as in WPS 2.2) and used for
     ! ARPS2WRF program.
     nzsoil = 6
     DO j = 1,ny
       DO i = 1,nx
         zpsoil(i,j,1) = hterain(i,j)
         zpsoil(i,j,2) = hterain(i,j) - 0.05
         zpsoil(i,j,3) = hterain(i,j) - 0.25
         zpsoil(i,j,4) = hterain(i,j) - 0.70
         zpsoil(i,j,5) = hterain(i,j) - 1.50
         zpsoil(i,j,6) = hterain(i,j) - 3.00
       END DO
     END DO

   ELSE IF(soilstrhopt == 10)THEN ! Apr. 2011: special for ECMWF datasets
     ! to WRF-ARW.
     ! This option is not meant to run ARPS itself, it's only
     ! for producing initial tsoil and qsoil interpolation on
     ! multiple soil layers (set nzsoil = 6 for NARR and NAM data, and
     ! assign zpsoil the same values as in WPS 2.2) and used for
     ! ARPS2WRF program.
     nzsoil = 6
     DO j = 1,ny
       DO i = 1,nx
         zpsoil(i,j,1) = hterain(i,j)
         zpsoil(i,j,2) = hterain(i,j) - 0.035
         zpsoil(i,j,3) = hterain(i,j) - 0.175
         zpsoil(i,j,4) = hterain(i,j) - 0.64
         zpsoil(i,j,5) = hterain(i,j) - 1.775
         zpsoil(i,j,6) = hterain(i,j) - 3.00
       END DO
     END DO

   END IF    !  end of soilstrhopt if block

!  soil grid variables are set.

  IF (soilstrhopt /= 0) THEN  ! Added by Y. Wang on Jan 08, 2010
    IF (soilstrhopt /= 9 .AND. soilstrhopt /= 10) THEN
      WRITE(6,'(1x,a,/,7x,a,/,7x,a,/,31x,a)')                           &
        'INFO: soilstrhopt should be 0 only because I do not know ',    &
        'who has commented out the code for soilstrhopt /= 0.',         &
        'There must be a reason though strange?',                       &
        ' - WYH.'
    END IF
    DO k = 1,nzsoil               ! The following is to avoid float exception
      DO j = 1,ny
        DO i = 1,nx
          j3soil(i,j,k)    = 1.0
          j3soilinv(i,j,k) = 1.0
        END DO
      END DO
    END DO
  END IF                      ! Added by Y. Wang on Jan 08, 2010

  RETURN
END SUBROUTINE inisoilgrd
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE JACOB                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE jacob(nx,ny,nz,x,y,z,zp,j1,j2,j3,mp_tem)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate transformation Jacobians J1, J2 and J3.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  3/17/1991.
!
!  MODIFICATION HISTORY:
!
!  4/10/93 (M. Xue & Hao Jin)
!  Add the terrain.
!
!  9/2/94 (M. Xue)
!  Loop 710 that resets j2 on north and south boundary to
!  zero gradient deleted. It shouldn't have been there.
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!  OUTPUT:
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space
!             (m)
!
!    j1       Coordinate transformation Jacobian -d(zp)/dx
!    j2       Coordinate transformation Jacobian -d(zp)/dy
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!
!  WORK ARRAY:
!
!    mp_tem   Used for message passing, NOTE: the shape of this array
!             may different from that in real paramenter.
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

  INTEGER :: nx,ny,nz     ! The number of grid points in 3 directions

  REAL :: x(nx)           ! x-coord. of the physical and
                          ! computational grid. Defined at u-point.
  REAL :: y(ny)           ! y-coord. of the physical and
                          ! computational grid. Defined at v-point.
  REAL :: z(nz)           ! z-coord. of the computational grid.
                          ! Defined at w-point on the staggered grid.
  REAL :: zp(nx,ny,nz)    ! the physical height coordinate defined at
                          ! w-point of the staggered grid.

  REAL :: j1(nx,ny,nz)    ! Coordinate transformation Jacobian
                          ! -d(zp)/dx.
  REAL :: j2(nx,ny,nz)    ! Coordinate transformation Jacobian
                          ! -d(zp)/dy.
  REAL :: j3(nx,ny,nz)    ! Coordinate transformation Jacobian
                          ! d(zp)/dz.

  REAL :: mp_tem(nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k

  INTEGER :: istat
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'bndry.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  IF (mp_opt > 0) THEN
!    ALLOCATE (mp_tem(MAX(nx,ny)*nz),stat=istat)
!    IF (istat /= 0) THEN
!      CALL arpsstop ("arpsstop called from JACOB: ERROR allocating mp_tem" &
!           ,1)
!    END IF
!  END IF
!
!-----------------------------------------------------------------------
!
!  Calculate transformation Jacobian J1 defined as -del(zp)/del(x)
!
!-----------------------------------------------------------------------
!
  DO k=1,nz
    DO j=1,ny-1
      DO i=2,nx-1
        j1(i,j,k) = -2 * (zp(i,j,k)-zp(i-1,j,k)) / (x(i+1)-x(i-1))
      END DO
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  X - boundary conditions of j1
!
!-----------------------------------------------------------------------
!

  IF (mp_opt > 0) THEN
    CALL acct_interrupt(mp_acct)
    CALL mpsendrecv2dew(j1,nx,ny,nz,ebc,wbc,1,mp_tem)
  END IF

  CALL acct_interrupt(bc_acct)

  IF(wbc == 1) THEN            ! Rigid wall boundary condition

    DO k=1,nz
      DO j=1,ny-1
        j1( 1,j,k)=j1( 3 ,j,k)
      END DO
    END DO

  ELSE IF(wbc == 2) THEN        ! Periodic boundary condition.

    IF(mp_opt == 0) THEN
      DO k=1,nz
      DO j=1,ny-1
        j1( 1,j,k)=j1(nx-2,j,k)
      END DO
      END DO
    END IF

  ELSE IF(wbc /= 0) THEN

    DO k=1,nz
      DO j=1,ny-1
        j1( 1,j,k)=j1( 2 ,j,k)
      END DO
    END DO

  END IF

  IF(ebc == 1) THEN             ! Rigid wall boundary condition

    DO k=1,nz
      DO j=1,ny-1
        j1(nx,j,k)=j1(nx-2,j,k)
      END DO
    END DO

  ELSE IF(ebc == 2) THEN         ! Periodic boundary condition.

    IF(mp_opt == 0) THEN
      DO k=1,nz
      DO j=1,ny-1
        j1(nx,j,k)=j1( 3 ,j,k)
      END DO
      END DO
    END IF

  ELSE IF(ebc /= 0) THEN

    DO k=1,nz
      DO j=1,ny-1
        j1(nx,j,k)=j1(nx-1,j,k)
      END DO
    END DO

  END IF

  DO k=1,nz
    DO i=1,nx
      j1(i,ny,k) = j1(i,ny-1,k)
    END DO
  END DO

  CALL acct_stop_inter
!
!-----------------------------------------------------------------------
!
!  Calculate transformation Jacobian J2 defined as -del(zp)/del(y)
!
!-----------------------------------------------------------------------
!
  DO k=1,nz
    DO i=1,nx-1
      DO j=2,ny-1
        j2(i,j,k) = -2 * (zp(i,j,k)-zp(i,j-1,k)) / (y(j+1)-y(j-1))
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Y - boundary conditions of j2
!
!-----------------------------------------------------------------------
!

  IF (mp_opt > 0) THEN
    CALL acct_interrupt(mp_acct)
    CALL mpsendrecv2dns(j2,nx,ny,nz,nbc,sbc,2,mp_tem)
  END IF

  CALL acct_interrupt(bc_acct)

  IF(sbc == 1) THEN            ! Rigid wall boundary condition

    DO k=1,nz
      DO i=1,nx-1
        j2(i, 1,k)=j2(i, 3 ,k)
      END DO
    END DO

  ELSE IF(sbc == 2) THEN        ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO k=1,nz
      DO i=1,nx-1
        j2(i, 1,k)=j2(i,ny-2,k)
      END DO
      END DO
    END IF

  ELSE IF(sbc /= 0) THEN

    DO k=1,nz
      DO i=1,nx-1
        j2(i, 1,k)=j2(i, 2 ,k)
      END DO
    END DO

  END IF

  IF(nbc == 1) THEN             ! Rigid wall boundary condition

    DO k=1,nz
      DO i=1,nx-1
        j2(i,ny,k)=j2(i,ny-2,k)
      END DO
    END DO

  ELSE IF(nbc == 2) THEN         ! Periodic boundary condition.

    IF (mp_opt == 0) THEN
      DO k=1,nz
      DO i=1,nx-1
        j2(i,ny,k)=j2(i, 3 ,k)
      END DO
      END DO
    END IF

  ELSE IF(nbc /= 0) THEN

    DO k=1,nz
      DO i=1,nx-1
        j2(i,ny,k)=j2(i,ny-1,k)
      END DO
    END DO

  END IF

  DO k=1,nz
    DO j=1,ny
      j2(nx,j,k) = j2(nx-1,j,k)
    END DO
  END DO

  CALL acct_stop_inter
!
!-----------------------------------------------------------------------
!
!  Calculate transformation Jacobian J3 defined as del(zp)/del(z)
!
!-----------------------------------------------------------------------
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        j3(i,j,k) = (zp(i,j,k+1)-zp(i,j,k))/(z(k+1)-z(k))
      END DO
    END DO
  END DO

  DO k=1,nz-1
    DO j=1,ny-1
      j3(nx,j,k) = j3(nx-1,j,k)
    END DO
  END DO

  DO k=1,nz-1
    DO i=1,nx
      j3(i,ny,k) = j3(i,ny-1,k)
    END DO
  END DO

  DO j=1,ny
    DO i=1,nx
      j3(i,j,nz) = j3(i,j,nz-1)
    END DO
  END DO

!  IF (mp_opt > 0) THEN
!    DEALLOCATE (mp_tem,stat=istat)
!  END IF

  RETURN
END SUBROUTINE jacob
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE INITGRDVAR                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE initgrdvar(nx,ny,nz,nzsoil,nts,nstyps,exbcbufsz,             &
           x,y,z,zp,zpsoil,hterain,mapfct,                              &
           j1,j2,j3,j3soil,aj3x,aj3y,aj3z,j3inv,j3soilinv,              &
           u,v,w,wcont,ptprt,pprt,qv,qscalar,tke,                       &
           udteb, udtwb, vdtnb, vdtsb,                                  &
           pdteb,pdtwb ,pdtnb ,pdtsb,                                   &
           trigs1,trigs2,ifax1,ifax2,                                   &
           wsave1,wsave2,vwork1,vwork2,                                 &
           ubar,vbar,wbar,ptbar,pbar,ptbari,pbari,                      &
           rhostr,rhostri,qvbar,ppi,csndsq,                             &
           soiltyp,stypfrct,vegtyp,lai,roufns,veg,                      &
           tsoil,qsoil,wetcanp,snowdpth,qvsfc,                          &
           ptcumsrc,qcumsrc,w0avg,nca,kfraincv,                         &
           cldefi,xland,bmjraincv,                                      &
           raing,rainc,prcrate,exbcbuf,                                 &
           radfrc,radsw,rnflx,radswnet,radlwin,                         &
           usflx,vsflx,ptsflx,qvsflx,                                   &
           tem1soil,tem2soil,tem3soil,tem4soil,tem5soil,                &
           tem1,tem2,tem3,tem4,tem5,                                    &
           tem6,tem7,tem8,tem9)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Initialize the model array variables.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  3/17/1992.
!
!  MODIFICATION HISTORY:
!
!  5/02/92 (M. Xue)
!  Added full documentation.
!
!  5/03/92 (M. Xue)
!  Further documentation.
!
!  10/7/1992 (M. Xue)
!  Added call to subroutine extinit, the option three of
!  initialization.
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!  02/07/1995 (Yuhe Liu)
!  Added a new 2-D permanent array, veg(nx,ny), to the argument list
!
!  10/31/95  (D. Weber)
!  Added trigs1,trigs2,ifax1,ifax2 for use in the fft code for the
!  upper radiation condition.
!
!  1/22/96 (Donghai Wang & Yuhe Liu)
!  Added the map projection factor to ARPS governing equations.
!
!  07/22/97 (D. Weber)
!  Added wsave1,wsave2,vwork1,vwork2 for use in the even fft version
!  of the upper w-p radiation condition (fftopt=2).
!
!  08/01/97 (Zonghui Huo)
!  Added Kain-fritsch cumulus parameterization scheme.
!
!  11/06/97 (D. Weber)
!  Added three additional levels to the mapfct array.  The three
!  levels (4,5,6) represent the inverse of the first three in order.
!  The inverse map factors are computed to improve efficiency.
!
!  12/05/97 (K. Brewster)
!  Added argument, nt, so that routines that do not require more
!  than one time level can be initialized using less memory.
!
!  4/15/1998 (Donghai Wang)
!  Added the source terms to the right hand terms of the qc,qr,qi,qs
!  equations due to the K-F cumulus parameterization.
!
!  4/15/1998 (Donghai Wang)
!  Added the running average vertical velocity (array w0avg)
!  for the K-F cumulus parameterization scheme.
!
!  11/18/98 (Keith Brewster)
!  Changed pibar to ppi (full pi).
!
!  12/09/1998 (Donghai Wang)
!  Added the snow cover.
!
!  11/06/2001 (Yunheng Wang)
!  Added mpupdatei calling for ict/icb to solve radiation forcing
!  difference for MPI run.
!
!  2002/02/28 (Gene Bassett)
!  Replaced mpupdatei for ict/icb to a call to mpmax0.
!
!  13 March 2002 (Eric Kemp)
!  Added arrays for WRF BMJ cumulus scheme.
!
!  04/10/2002 (Yunheng Wang)
!  Subsituted mpupdatei calls for mpmax0 calls again
!  because mpmax0 calls were not correct.
!
!  05/18/2002 (Dan Weber and Jerry Brotzge)
!  Added new soil model variables.
!
!-----------------------------------------------------------------------
!
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nts      Number of time levels to be initialized.
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space
!             (m)
!    zpsoil   Vertical coordinate of grid points in the soil model
!             in physical space (m).
!    hterain  The height of terrain (m)
!
!    mapfct   Map factors at scalar, u and v points
!
!    j1       Coordinate transformation Jacobian -d(zp)/dx
!    j2       Coordinate transformation Jacobian -d(zp)/dy
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!    j3soil   Soil coordinate transformation Jacobian  d(zpsoil)/dz
!    aj3x     Avgx of the coordinate transformation Jacobian  d(zp)/dz
!    aj3y     Avgy of the coordinate transformation Jacobian  d(zp)/dz
!    aj3z     Avgz of the coordinate transformation Jacobian  d(zp)/dz
!    j3soilinv Inverse of the soil coordinate transformation j3soil
!
!    trigs1   Array containing pre-computed trig function for fftopt=1.
!    trigs1   Array containing pre-computed trig function for fftopt=1.
!    ifax1    Array containing the factors of nx for fftopt=1.
!    ifax2    Array containing the factors of ny for fftopt=1.
!
!    vwork1   2-D work array for fftopt=2.
!    vwork2   2-D work array for fftopt=2.
!    wsave1   Work array for fftopt=2.
!    wsave2   Work array for fftopt=2.
!
!  OUTPUT:
!
!    u        x component of velocity at all time levels (m/s)
!    v        y component of velocity at all time levels (m/s)
!    w        Vertical component of Cartesian velocity
!             at all time levels (m/s)
!    wcont    Contravariant vertical velocity (m/s)
!    ptprt    Perturbation potential temperature at all time levels
!             (K)
!    pprt     Perturbation pressure at all time levels (Pascal)
!    qv       Water vapor specific humidity at all time levels
!             (kg/kg)
!    qc       Cloud water mixing ratio at all time levels (kg/kg)
!    qr       Rainwater mixing ratio at all time levels (kg/kg)
!    qi       Cloud ice mixing ratio at all time levels (kg/kg)
!    qs       Snow mixing ratio at all time levels (kg/kg)
!    qh       Hail mixing ratio at all time levels (kg/kg)
!    tke      Turbulent Kinetic Energy ((m/s)**2)
!
!    udteb    Time tendency of u field at east boundary (m/s**2)
!    udtwb    Time tendency of u field at west boundary (m/s**2)
!
!    vdtnb    Time tendency of v field at north boundary (m/s**2)
!    vdtsb    Time tendency of v field at south boundary (m/s**2)
!
!    pdteb    Time tendency of pprt field at east boundary (PASCAL/s)
!    pdtwb    Time tendency of pprt field at west boundary (PASCAL/s)
!    pdtnb    Time tendency of pprt field at north boundary
!             (PASCAL/s)
!    pdtsb    Time tendency of pprt field at south boundary
!             (PASCAL/s)
!
!    ubar     Base state zonal velocity component (m/s)
!    vbar     Base state meridional velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    ptbari   Inverse Base state potential temperature (K)
!    pbari    Inverse Base state pressure (Pascal)
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!    rhostri  Inverse base state density rhobar times j3 (kg/m**3)
!    qvbar    Base state water vapor specific humidity (kg/kg)
!    ppi      Exner function.
!    csndsq   Sound wave speed squared.
!
!    soiltyp  Soil type
!    vegtyp   Vegetation type
!    lai      Leaf Area Index
!    roufns   Surface roughness
!    veg      Vegetation fraction
!
!    tsoil    Soil temperature (K)
!    qsoil    Soil moisture (m**3/m**3)
!    wetcanp  Canopy water amount
!    snowdpth Snow depth (m)
!    qvsfc    Effective S.H. at sfc.
!
!    ptcumsrc Source term in pt-equation due to cumulus parameterization
!    qcumsrc Source term in water equations due to cumulus parameterization
!    kfraincv K-F convective rainfall (cm)
!    nca      K-F counter for CAPE release
!    cldefi   BMJ cloud efficiency
!    xland    BMJ land/sea mask
!    bmjraincv BMJ convective rainfall (cm)
!
!    radfrc   Radiation forcing (K)
!    radsw    Solar radiation reaching the surface
!    rnflx    Net absorbed radiation by the surface
!    radswnet Net shortwave radiation
!    radlwin  Incoming longwave radiation
!
!    raing    Grid supersaturation rain
!    rainc    Cumulus convective rain
!
!  WORK ARRAYS:
!
!    tem1soil Soil model temporary work array.
!    tem2soil Soil model temporary work array.
!    tem3soil Soil model temporary work array.
!    tem4soil Soil model temporary work array.
!    tem5soil Soil model temporary work array.
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
!    tem4     Temporary work array.
!    tem5     Temporary work array.
!    tem6     Temporary work array.
!    tem7     Temporary work array.
!    tem8     Temporary work array.
!    tem9     Temporary work array.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'bndry.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nts               ! Number of time levels to be initialized.
  INTEGER :: tpast             ! Index of time level for the past time.
  INTEGER :: tpresent          ! Index of time level for the present
                               ! time.
  INTEGER :: tfuture           ! Index of time level for the future
                               ! time.
  INTEGER :: nx,ny,nz          ! The number of grid points in 3
                               ! directions
  INTEGER :: nzsoil            ! Number of grid points in the -z direction

  REAL :: x     (nx)           ! x-coord. of the physical and compu-
                               ! tational grid. Defined at u-point.
  REAL :: y     (ny)           ! y-coord. of the physical and compu-
                               ! tational grid. Defined at v-point.
  REAL :: z     (nz)           ! z-coord. of the computational grid.
                               ! Defined at w-point on the staggered
                               ! grid.
  REAL :: zp    (nx,ny,nz)     ! Physical height coordinate defined at
                               ! w-point of the staggered grid.
  REAL :: zpsoil(nx,ny,nzsoil) ! The physical height coordinate defined
                               ! at the center of a soil layer(m).
  REAL :: hterain(nx,ny)       ! The height of the terrain.

  REAL :: mapfct(nx,ny,8)      ! Map factors at scalar, u and v points

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! defined as - d( zp )/d( x ).
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! defined as - d( zp )/d( y ).
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! defined as d( zp )/d( z ).
  REAL :: j3soil(nx,ny,nzsoil) ! Coordinate transformation Jacobian
                               ! defined as d( zpsoil )/d( zsoil ).

  REAL :: aj3x  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE X-DIR.
  REAL :: aj3y  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Y-DIR.
  REAL :: aj3z  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Z-DIR.
  REAL :: j3inv (nx,ny,nz)     ! Inverse of j3
  REAL :: j3soilinv(nx,ny,nzsoil) ! Inverse of J3soil.

  REAL :: trigs1(3*(nx-1)/2+1) ! Array containing pre-computed trig
                               ! function for fftopt=1.
  REAL :: trigs2(3*(ny-1)/2+1) ! Array containing pre-computed trig
                               ! function for fftopt=1.
  INTEGER :: ifax1(13)         ! Array containing the factors of nx for
                               ! fftopt=1.
  INTEGER :: ifax2(13)         ! Array containing the factors of ny for
                               ! fftopt=1.

  REAL :: vwork1 (nx+1,ny+1)   ! 2-D work array for fftopt=2.
  REAL :: vwork2 (ny,nx+1)     ! 2-D work array for fftopt=2.
  REAL :: wsave1 (3*(ny-1)+15) ! Work array for fftopt=2.
  REAL :: wsave2 (3*(nx-1)+15) ! Work array for fftopt=2.

  REAL :: u     (nx,ny,nz,nts) ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nts) ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz,nts) ! Total w-velocity (m/s)
  REAL :: wcont (nx,ny,nz)     ! Contravariant vertical velocity (m/s)
  REAL :: ptprt (nx,ny,nz,nts) ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz,nts) ! Perturbation pressure (Pascal)

  REAL :: qv    (nx,ny,nz,nts) ! Water vapor specific humidity (kg/kg)
  REAL :: qscalar(nx,ny,nz,nts,nscalar)
  REAL :: tke   (nx,ny,nz,nts) ! Turbulent Kinetic Energy ((m/s)**2)

  REAL :: udteb (ny,nz)        ! T-tendency of u at e-boundary (m/s**2)
  REAL :: udtwb (ny,nz)        ! T-tendency of u at w-boundary (m/s**2)

  REAL :: vdtnb (nx,nz)        ! T-tendency of v at n-boundary (m/s**2)
  REAL :: vdtsb (nx,nz)        ! T-tendency of v at s-boundary (m/s**2)

  REAL :: pdteb (ny,nz)        ! T-tendency of pprt at e-boundary
                               ! (PASCAL/s)
  REAL :: pdtwb (ny,nz)        ! T-tendency of pprt at w-boundary
                               ! (PASCAL/s)
  REAL :: pdtnb (nx,nz)        ! T-tendency of pprt at n-boundary
                               ! (PASCAL/s)
  REAL :: pdtsb (nx,nz)        ! T-tendency of pprt at s-boundary
                               ! (PASCAL/s)

  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: wbar  (nx,ny,nz)     ! Base state w-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal).
  REAL :: ptbari(nx,ny,nz)     ! Inverse Base state pot. temperature (K)
  REAL :: pbari (nx,ny,nz)     ! Inverse Base state pressure (Pascal).
  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.
  REAL :: rhostri(nx,ny,nz)    ! Inv. base state density rhobar times j3.
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific
                               ! humidity(kg/kg)
  REAL :: ppi   (nx,ny,nz)     ! Exner function.
  REAL :: csndsq(nx,ny,nz)     ! Sound wave speed squared.

  INTEGER :: nstyps                  ! Number of soil types
  INTEGER :: soiltyp(nx,ny,nstyps)   ! Soil types at grid points
  REAL :: stypfrct(nx,ny,nstyps)  ! Fraction of soil types
  INTEGER :: vegtyp (nx,ny)          ! Vegetation type
  REAL :: lai    (nx,ny)          ! Leaf Area Index
  REAL :: roufns (nx,ny)          ! Surface roughness
  REAL :: veg    (nx,ny)          ! Vegetation fraction

  REAL :: qvsfc(nx,ny,0:nstyps)    ! Effective qv at sfc.
  REAL :: tsoil(nx,ny,nzsoil,0:nstyps)  ! Deep soil temperature (K)
                                     ! (in deep 1 m layer)
  REAL :: qsoil(nx,ny,nzsoil,0:nstyps)    ! Surface soil moisture in the top
                                     ! 1 cm layer
  REAL :: wetcanp(nx,ny,0:nstyps)    ! Canopy water amount
  REAL :: snowdpth(nx,ny)            ! Snow depth (m)

  REAL :: ptcumsrc(nx,ny,nz)   ! Source term in pt-equation due
                               ! to cumulus parameterization
  REAL :: qcumsrc(nx,ny,nz,5)  ! Source term in water equations due
                               ! to cumulus parameterization:
                               ! qcumsrc(1,1,1,1) for qv equation
                               ! qcumsrc(1,1,1,2) for qc equation
                               ! qcumsrc(1,1,1,3) for qr equation
                               ! qcumsrc(1,1,1,4) for qi equation
                               ! qcumsrc(1,1,1,5) for qs equation
  REAL :: w0avg(nx,ny,nz)      ! a closing running average vertical
                               ! velocity in 10min for K-F scheme
  REAL :: kfraincv(nx,ny)      ! K-F convective rainfall (cm)
  INTEGER :: nca(nx,ny)        ! K-F counter for CAPE release

!EMK BMJ
  REAL,INTENT(INOUT) :: cldefi(nx,ny)      ! BMJ cloud efficiency
  REAL,INTENT(INOUT) :: xland(nx,ny)       ! BMJ land mask
                                           !   (1.0 = land, 2.0 = sea)
  REAL,INTENT(INOUT) :: bmjraincv(nx,ny)   ! BMJ convective rainfall (cm)
!EMK END

  REAL :: radfrc(nx,ny,nz)     ! Radiation forcing (K/s)
  REAL :: radsw (nx,ny)        ! Solar radiation reaching the surface
  REAL :: rnflx (nx,ny)        ! Net radiation flux absorbed by surface
  REAL :: radswnet(nx,ny)      ! Net shortwave radiation
  REAL :: radlwin(nx,ny)       ! Incoming longwave radiation

  REAL :: usflx (nx,ny)        ! Surface flux of u-momentum (kg/(m*s**2))
  REAL :: vsflx (nx,ny)        ! Surface flux of v-momentum (kg/(m*s**2))
  REAL :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m*s**2))
  REAL :: qvsflx(nx,ny)        ! Surface moisture flux (kg/(m**2*s))

  REAL :: raing(nx,ny)         ! Grid supersaturation rain
  REAL :: rainc(nx,ny)         ! Cumulus convective rain
  REAL :: prcrate(nx,ny,4)     ! precipitation rate (kg/(m**2*s))
                               ! prcrate(1,1,1) = total precipitation rate
                               ! prcrate(1,1,2) = grid scale precip. rate
                               ! prcrate(1,1,3) = cumulus precip. rate
                               ! prcrate(1,1,4) = microphysics precip. rate

  INTEGER :: exbcbufsz         ! EXBC buffer size
  REAL :: exbcbuf( exbcbufsz ) ! EXBC buffer array

  REAL :: tem1soil(nx,ny,nzsoil)  ! Temporary soil model work array.
  REAL :: tem2soil(nx,ny,nzsoil)  ! Temporary soil model work array.
  REAL :: tem3soil(nx,ny,nzsoil)  ! Temporary soil model work array.
  REAL :: tem4soil(nx,ny,nzsoil)  ! Temporary soil model work array.
  REAL :: tem5soil(nx,ny,nzsoil)  ! Temporary soil model work array.

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem4  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem5  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem6  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem7  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem8  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem9  (nx,ny,nz)     ! Temporary work array.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  INTEGER :: iskip
  REAL :: temp

  REAL :: alatpro(2)
  REAL :: sclf
  REAL :: dxscl             ! Model x-direction grid spacing
                            ! normalized by the map scale
                            ! dxscl=dx/sclf
  REAL :: dyscl             ! Model y-direction grid spacing
                            ! normalized by the map scale
                            ! dyscl=dy/sclf
  REAL :: xs,ys, swx,swy, ctrx, ctry
  REAL zpmax

  REAL :: rmin,rmax
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
!  If initopt = 1, initialize the model fields using intialization
!                  routines. Typically from an analytical
!                  definition of initial perturbations.
!
!  If initopt = 2, initialize the model fields from a restart file
!                  produced by a previous model run.
!
!  If initopt = 3, initialize the model fields by reading in an
!                  external data file.
!
!  If initopt = 4, initialize the model fields by reading in restart
!                  file first then an external data file.
!
!-----------------------------------------------------------------------
!
  IF( nts > 1 ) THEN
    tpast=1
    tpresent=2
    tfuture=3
  ELSE
    tpast=1
    tpresent=1
    tfuture=1
  END IF

  IF( initopt == 1 ) THEN
!
!-----------------------------------------------------------------------
!
!  Initialization of model GRiD definition arrays.
!
!-----------------------------------------------------------------------
!

    CALL inigrd(nx,ny,nz,nzsoil,x,y,z,zp,zpsoil,                       &
                hterain,mapfct,j1,j2,j3,j3soil,                        &
                j3soilinv,tem1,tem2,tem3)
!
!-----------------------------------------------------------------------
!
!  Define the base state atmospheric variables.
!
!-----------------------------------------------------------------------
!
!    IF ((inibasopt == 1) .AND. (max_fopen < nprocs)) THEN
!      CALL wrtcomment("ERROR: for message passing version with "//      &
!                    "inibasopt=1, max_fopen (in arps.input)",1)
!      CALL arpsstop ("arpsstop called from initgrdvar due to nproc_x-y  &
!                    & nproc_x*nproc_y (in arps.input).",1)
!    END IF

     IF(inibasopt == 1) THEN
       iskip = nproc_x * nproc_y
     ELSE
       iskip = max_fopen
     END IF

!  blocking inserted for ordering i/o for message passing
    DO i=0,nprocs-1,iskip
      IF(myproc >= i.AND.myproc <= i+iskip-1)THEN
        CALL inibase(nx,ny,nz, ubar,vbar,ptbar,pbar,ptbari,pbari,       &
                   rhostr,rhostri,qvbar,                                &
                   x,y,z,zp,j3, tem1,tem2,tem3,tem4,tem5,tem6)
      END IF
      IF (mp_opt > 0) CALL mpbarrier
    END DO
!
!-----------------------------------------------------------------------
!
!  Initialize time dependent model variables.
!
!-----------------------------------------------------------------------
!
    CALL initdvr(nx,ny,nz,nts,                                          &
                 ubar,vbar,ptbar,pbar,rhostr,qvbar,                     &
                 x,y,z,zp,hterain, j1,j2,j3,                            &
                 u,v,w,ptprt,pprt,qv,qscalar,                           &
                 ptcumsrc,qcumsrc,raing,rainc,prcrate,tem1)
!
!-----------------------------------------------------------------------
!
!  Set the time tendencies of u, v and pprt on the lateral boundaries
!  to zero for the initial time step
!
!  These tendencies will be used by lateral boundary condition options
!  4 and 5 only.
!
!-----------------------------------------------------------------------
!
    DO k=1,nz
      DO j=1,ny
        udteb(j,k) = 0.0
        udtwb(j,k) = 0.0
        pdteb(j,k) = 0.0
        pdtwb(j,k) = 0.0
      END DO
    END DO

    DO k=1,nz
      DO i=1,nx
        vdtnb(i,k) = 0.0
        vdtsb(i,k) = 0.0
        pdtnb(i,k) = 0.0
        pdtsb(i,k) = 0.0
      END DO
    END DO
    wbar(:,:,:) = 0.0

  ELSE IF( initopt == 2 .or. initopt == 4) THEN   ! restart run
!
!-----------------------------------------------------------------------
!
!  Read in restart data from a restart file to initialize fields
!  u,v,w,ptprt,pprt,qv,qc,qr,qi,qs,qh at time tpast and tpresent,
!  the base state variables ubar,vbar,ptbar,pbar,rhostr,qvbar,
!  and the time tendencies of variables at the lateral boundries.
!
!  Fields at tfuture are set to equal to those at tpresent.
!
!  This subroutine also sets the value of tstart to the restart
!  data time. The value from input file in over-written.
!
!-----------------------------------------------------------------------
!
    IF(mp_opt >0 .AND. readsplit(FINDX_R) > 0) THEN

      CALL rstinsplit(nx,ny,nz,nzsoil,nts, nstyps, exbcbufsz,           &
                 u,v,w,ptprt,pprt,qv,qscalar,tke,                       &
                 udteb, udtwb, vdtnb, vdtsb,                            &
                 pdteb, pdtwb, pdtnb, pdtsb,                            &
                 ubar,vbar,ptbar,pbar,rhostr,qvbar,                     &
                 x,y,z,zp,zpsoil,hterain,mapfct,j1,j2,j3,j3soil,        &
                 soiltyp,stypfrct,vegtyp,lai,roufns,veg,                &
                 tsoil,qsoil,wetcanp,snowdpth,qvsfc,                    &
                 ptcumsrc,qcumsrc,w0avg,nca,kfraincv,                   &
                 cldefi,xland,bmjraincv,                                &
                 radfrc,radsw,rnflx,radswnet,radlwin,                   &
                 raing,rainc,prcrate,exbcbuf,tem1,tem2)

    ELSE

      CALL rstin(nx,ny,nz,nzsoil,nts, nstyps, exbcbufsz,                &
                 u,v,w,ptprt,pprt,qv,qscalar,tke,                       &
                 udteb, udtwb, vdtnb, vdtsb,                            &
                 pdteb, pdtwb, pdtnb, pdtsb,                            &
                 ubar,vbar,ptbar,pbar,rhostr,qvbar,                     &
                 x,y,z,zp,zpsoil,hterain,mapfct,j1,j2,j3,j3soil,        &
                 soiltyp,stypfrct,vegtyp,lai,roufns,veg,                &
                 tsoil,qsoil,wetcanp,snowdpth,qvsfc,                    &
                 ptcumsrc,qcumsrc,w0avg,nca,kfraincv,                   &
                 cldefi,xland,bmjraincv,                                &
                 radfrc,radsw,rnflx,radswnet,radlwin,                   &
                 raing,rainc,prcrate,exbcbuf,tem1,tem2)
    END IF
!
!-----------------------------------------------------------------------
!
!  Set map projection parameters which were not stored in restart
!  data file.
!
!-----------------------------------------------------------------------
!
    alatpro(1) = trulat1
    alatpro(2) = trulat2

    IF( sclfct /= 1.0) THEN
      sclf  = 1.0/sclfct
      dxscl = dx*sclf
      dyscl = dy*sclf
    ELSE
      sclf  = 1.0
      dxscl = dx
      dyscl = dy
    END IF

    CALL setmapr( mapproj,sclf,alatpro,trulon )

    CALL lltoxy( 1,1, ctrlat,ctrlon, ctrx, ctry )

!       swx = ctrx - (float(nx-3)/2.) * dxscl
!       swy = ctry - (float(ny-3)/2.) * dyscl
    swx = ctrx - (FLOAT(nproc_x*(nx-3))/2.) * dxscl
    swy = ctry - (FLOAT(nproc_y*(ny-3))/2.) * dyscl

    CALL setorig( 1, swx, swy)
                               ! set up the model origin to the coord.

    CALL setcornerll(nx,ny, x, y)

    DO k=1,nz-1
      DO j=1,ny
        DO i=1,nx
          tem1(i,j,k) = rhostr(i,j,k)/j3(i,j,k)
          tem2(i,j,k) = (zp(i,j,k)+zp(i,j,k+1))*0.5
        END DO
      END DO
    END DO
    CALL writesnd(nx,ny,nz,ubar,vbar,ptbar,pbar,qvbar,zp, tem1, tem2)

  ENDIF

  IF( initopt == 3 .or. initopt == 4) THEN  ! External data input.
!
!-----------------------------------------------------------------------
!
!  Read in externally supplied initial fields.  These fields include
!  u, v, w, ptprt, pprt, qv, qc, qr, qi, qs, and qh at time level
!  tpresent, and the base state variables ubar, vbar, ptbar, pbar,
!  rhostr,qvbar.
!
!  Fields at tpast and tfuture are set to their values at tpresent.
!
!-----------------------------------------------------------------------
!
    CALL extinit(nx,ny,nz,nzsoil,nts,nstyps,                            &
                 u,v,w,ptprt,pprt,qv,qscalar,tke,                       &
                 ubar,vbar,ptbar,pbar,rhostr,qvbar,                     &
                 x,y,z,zp,zpsoil,hterain,j1,j2,j3,j3soil,               &
                 soiltyp,stypfrct,vegtyp,lai,roufns,veg,                &
                 tsoil,qsoil,wetcanp,snowdpth,qvsfc,                    &
                 ptcumsrc,qcumsrc,raing,rainc,prcrate,                  &
                 radfrc,radsw,rnflx,radswnet,radlwin,                   &
                 usflx,vsflx,ptsflx,qvsflx,                             &
                 tem1,tem2,tem3,tem4,tem5,wbar,tem7,tem8,tem9)
!
!-----------------------------------------------------------------------
!
!  Set map projection parameters which were not stored in history
!  data file.
!
!-----------------------------------------------------------------------
!
    alatpro(1) = trulat1
    alatpro(2) = trulat2

    IF( sclfct /= 1.0) THEN
      sclf  = 1.0/sclfct
      dxscl = dx*sclf
      dyscl = dy*sclf
    ELSE
      sclf  = 1.0
      dxscl = dx
      dyscl = dy
    END IF

    CALL setmapr( mapproj,sclf,alatpro,trulon )

    CALL lltoxy( 1,1, ctrlat,ctrlon, ctrx, ctry )

!       swx = ctrx - (float(nx-3)/2.) * dxscl
!       swy = ctry - (float(ny-3)/2.) * dyscl
    swx = ctrx - (FLOAT(nproc_x*(nx-3))/2.) * dxscl
    swy = ctry - (FLOAT(nproc_y*(ny-3))/2.) * dyscl

    CALL setorig( 1, swx, swy)
                               ! set up the model origin to the coord.

    IF ( mpfctopt /= 0 ) THEN

      DO j=1,ny-1
        DO i=1,nx-1
          xs = 0.5*(x(i)+x(i+1))
          ys = 0.5*(y(j)+y(j+1))
          CALL xytomf( 1,1,xs,ys,mapfct(i,j,1) )
          mapfct(i,j,4) = 1.0/mapfct(i,j,1)
          mapfct(i,j,7) = mapfct(i,j,1)*mapfct(i,j,1)
          mapfct(i,j,8) = 0.25*mapfct(i,j,1)
        END DO
      END DO

      DO j=1,ny-1
        DO i=1,nx
          ys = 0.5*(y(j)+y(j+1))
          CALL xytomf( 1,1,x(i),ys,mapfct(i,j,2) )
          mapfct(i,j,5) = 1.0/mapfct(i,j,2)
        END DO
      END DO

      DO j=1,ny
        DO i=1,nx-1
          xs = 0.5*(x(i)+x(i+1))
          CALL xytomf( 1,1,xs,y(j),mapfct(i,j,3) )
          mapfct(i,j,6) = 1.0/mapfct(i,j,3)
        END DO
      END DO

    ELSE

      DO k=1,7
        DO j=1,ny
          DO i=1,nx
            mapfct(i,j,k) = 1.0
            mapfct(i,j,8) = 0.25
          END DO
        END DO
      END DO

    END IF

    CALL setcornerll(nx,ny, x, y)

    IF( initopt == 3 ) then
!
!-----------------------------------------------------------------------
!
!  Set the time tendencies of u, v and pprt on the lateral boundaries
!  to zero for the initial time step
!
!  These tendencies will be used by lateral boundary condition options
!  4 and 5 only.
!
!-----------------------------------------------------------------------
!
      DO k=1,nz
        DO j=1,ny
          udteb(j,k) = 0.0
          udtwb(j,k) = 0.0
          pdteb(j,k) = 0.0
          pdtwb(j,k) = 0.0
        END DO
      END DO

      DO k=1,nz
        DO i=1,nx
          vdtnb(i,k) = 0.0
          vdtsb(i,k) = 0.0
          pdtnb(i,k) = 0.0
          pdtsb(i,k) = 0.0
        END DO
      END DO

    ENDIF

    DO k=1,nz-1
      DO j=1,ny
        DO i=1,nx
          tem1(i,j,k) = rhostr(i,j,k)/j3(i,j,k)
          tem2(i,j,k) = (zp(i,j,k)+zp(i,j,k+1))*0.5
        END DO
      END DO
    END DO
    CALL writesnd(nx,ny,nz,ubar,vbar,ptbar,pbar,qvbar,zp, tem1, tem2)

  END IF
!
!-----------------------------------------------------------------------
!
!  Define the reversed vertical indeces of height cldh2m and cldm2l
!  which represent the levels that separate high, middle, and low
!  clouds in the computation of the solar radiation transfer
!
!-----------------------------------------------------------------------
!
  ict = nz-2
  icb = nz-2

  DO k=nz-2,2,-1
    tem1(1,1,k) = (zp(1,1,k)-zp(1,1,2))*(zflat-zrefsfc)                 &
                 /(zflat-zp(1,1,2))+zrefsfc
  END DO

  DO k=nz-2,2,-1
    IF ( tem1(1,1,k) <= cldh2m) THEN
      ict = k
      EXIT
    END IF
  END DO

! for bit-for-bit MP accuracy:
  CALL mpupdatei(ict, 1)

  DO k=nz-2,2,-1
    IF ( tem1(1,1,k) <= cldm2l) THEN
      icb = k
      EXIT
    END IF
  END DO

! for bit-for-bit MP accuracy:
  CALL mpupdatei(icb, 1)

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        j3inv(i,j,k)=1.0/j3(i,j,k)
      END DO
    END DO
  END DO

  DO k=1,nzsoil
    DO j=1,ny-1
      DO i=1,nx-1
        j3soilinv(i,j,k)=1.0/j3soil(i,j,k)
      END DO
    END DO
  END DO

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=2,nx-1
        aj3x(i,j,k)=0.5*(j3(i,j,k)+j3(i-1,j,k))
      END DO
    END DO
  END DO

  IF (mp_opt > 0) THEN
    CALL acct_interrupt(mp_acct)
    CALL mpsendrecv2dew(aj3x,nx,ny,nz,ebc,wbc,1,tem2)
    !CALL mpsend2dns(aj3x,nx,ny,nz,1,mptag,tem2)
    !CALL mprecv2dns(aj3x,nx,ny,nz,1,mptag,tem2)
  END IF
  CALL acct_interrupt(bc_acct)
  CALL bcsu(nx,ny,nz,1,ny-1,1,nz-1,ebc,wbc,aj3x)
  CALL acct_stop_inter

  DO k=1,nz-1
    DO j=2,ny-1
      DO i=1,nx-1
        aj3y(i,j,k)=0.5*(j3(i,j,k)+j3(i,j-1,k))
      END DO
    END DO
  END DO

  IF (mp_opt > 0) THEN
    CALL acct_interrupt(mp_acct)
    !CALL mpsend2dew(aj3y,nx,ny,nz,2,mptag,tem2)
    !CALL mprecv2dew(aj3y,nx,ny,nz,2,mptag,tem2)
    CALL mpsendrecv2dns(aj3y,nx,ny,nz,nbc,sbc,2,tem2)
  END IF
  CALL acct_interrupt(bc_acct)
  CALL bcsv(nx,ny,nz,1,nx-1,1,nz-1,nbc,sbc,aj3y)
  CALL acct_stop_inter

  DO k=2,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        aj3z(i,j,k)=0.5*(j3(i,j,k)+j3(i,j,k-1))
      END DO
    END DO
  END DO

  CALL acct_interrupt(bc_acct)
  CALL bcsw(nx,ny,nz,1,nx-1,1,ny-1,tbc,bbc,aj3z)
  CALL acct_stop_inter
!
!-----------------------------------------------------------------------
!
!  Calculate the trigs1,trigs2,ifax1,ifax2 arrays by calling set99
!
!                              OR
!
!  calculate wsave1,wsave2 by calling vcosi.
!
!-----------------------------------------------------------------------
!
  DO i=1,13
    ifax1(i) = 0
    ifax2(i) = 0
  END DO

  DO i=1,3*(nx-1)/2+1
    trigs1(i) = 0
  END DO

  DO j=1,3*(ny-1)/2+1
    trigs2(j) = 0
  END DO

  DO j=1,ny+1
    DO i=1,nx+1
      vwork1(i,j) = 0.0
    END DO
  END DO

  DO i=1,nx+1
    DO j=1,ny
      vwork2(j,i) = 0.0
    END DO
  END DO

  DO i=1,3*(nx-1)+15
    wsave2(i) = 0.0
  END DO

  DO i=1,3*(ny-1)+15
    wsave1(i) = 0.0
  END DO

  IF ( tbc == 4 ) THEN  ! set up the fft work arrays for use in the
                        ! upper radiation boundary condition.
    IF ( fftopt == 1 ) THEN     ! set up periodic work arrays...
      IF ( ny == 4 ) THEN       ! set up trigs in x direction only
        CALL set99(trigs1,ifax1,nx-1)    ! NOTE: nx must be ODD!!!!
                                         ! and of special character...
                                         ! see fft99f.f for details....
      ELSE IF ( nx == 4 ) THEN     ! set up trigs in y direction only

        CALL set99(trigs2,ifax2,ny-1)    ! NOTE: ny must be ODD!!!!
                                         ! and of special character...
                                         ! see fft99f.f for details....
      ELSE    ! set up for 2-d transform...

        CALL set99(trigs1,ifax1,nx-1)    ! NOTE: nx must be ODD!!!!
                                         ! and of special character...
                                         ! see fft99f.f for details....
        CALL set99(trigs2,ifax2,ny-1)    ! NOTE: ny must be ODD!!!!
                                         ! and of special character...
                                         ! see fft99f.f for details....
      END IF

    ELSE IF ( fftopt == 2 ) THEN   ! set up the cos fft arrays...

      IF(ny == 4)THEN   ! set up function in x direction only...

        CALL vcosti(nx-1,wsave2)         ! nx should be even.

      ELSE IF(nx == 4)THEN ! set up function in y direction only...

        CALL vcosti(ny-1,wsave1)         ! ny should be even.

      ELSE   ! set up functions for 2-d transform...

        CALL vcosti(ny-1,wsave1)         ! ny should be even.
        CALL vcosti(nx-1,wsave2)         ! nx should be even.

      END IF  ! end of run type if block...

    END IF   ! end of fftopt if block.....

  END IF
!
!-----------------------------------------------------------------------
!
!  Find the lowest model layer (index rayklow) that is entirely or
!  partially contained in the specified Rayleigh damping (sponge)
!  layers.
!
!  The Rayleigh damping is then applied only to layers with
!  k greater than or equal to rayklow.
!
!-----------------------------------------------------------------------
!
  rayklow = nz-1

  DO k=nz-1,2,-1

    zpmax = zp(1,1,k)

    DO j=1,ny-1
      DO i=1,nx-1
        zpmax = MAX( zp(i,j,k), zpmax )
      END DO
    END DO

! for bit-for-bit accuracy with MP version:
    rmin = zpmax
    call mpmax0(zpmax,rmin)

    IF( zpmax < zbrdmp ) THEN
      rayklow = MAX(2, k+1)
      EXIT
    END IF

  END DO

!
!-----------------------------------------------------------------------
!
!  Calculate Exner function and store in ppi
!
!-----------------------------------------------------------------------
!
  CALL setppi(nx,ny,nz,nts,tpresent,pprt,pbar,ppi)
!
!-----------------------------------------------------------------------
!
!  Calculate and store the sound wave speed squared in csndsq.
!
!-----------------------------------------------------------------------
!
  IF(csopt == 1) THEN       ! Original definition of sound speed.
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          csndsq(i,j,k)= cpdcv*pbar(i,j,k)*j3(i,j,k)/rhostr(i,j,k)
        END DO
      END DO
    END DO

  ELSE IF(csopt == 2) THEN   ! Original sound speed multiplied
                             ! by a factor
    temp = csfactr**2*cpdcv
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          csndsq(i,j,k)= temp * pbar(i,j,k)*j3(i,j,k)/rhostr(i,j,k)
        END DO
      END DO
    END DO
  ELSE                      ! Specified constant sound speed.
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          csndsq(i,j,k)= csound * csound
        END DO
      END DO
    END DO
  END IF
!
!-----------------------------------------------------------------------
!
!  Fill the edges of base state arrays that are otherwise undefined.
!  This is done for safty reason only.
!
!-----------------------------------------------------------------------
!
  DO k=1,nz-1
    DO j=1,ny
      vbar  (nx,j,k) = vbar  (nx-1,j,k)
    END DO
  END DO

  DO k=1,nz-1
    DO i=1,nx
      ubar  (i,ny,k) = ubar  (i,ny-1,k)
    END DO
  END DO

  DO i=1,nx
    DO j=1,ny
      ubar  (i,j,nz) = ubar  (i,j,nz-1)
      vbar  (i,j,nz) = vbar  (i,j,nz-1)
    END DO
  END DO

  DO k=1,nz-1
    DO j=1,ny-1
      ptbar (nx,j,k) = ptbar (nx-1,j,k)
      pbar  (nx,j,k) = pbar  (nx-1,j,k)
      ppi   (nx,j,k) = ppi   (nx-1,j,k)
      qvbar (nx,j,k) = qvbar (nx-1,j,k)
      csndsq(nx,j,k) = csndsq(nx-1,j,k)
      rhostr(nx,j,k) = rhostr(nx-1,j,k)
    END DO
  END DO

  DO k=1,nz-1
    DO i=1,nx
      ptbar (i,ny,k) = ptbar (i,ny-1,k)
      pbar  (i,ny,k) = pbar  (i,ny-1,k)
      ppi   (i,ny,k) = ppi   (i,ny-1,k)
      qvbar (i,ny,k) = qvbar (i,ny-1,k)
      csndsq(i,ny,k) = csndsq(i,ny-1,k)
      rhostr(i,ny,k) = rhostr(i,ny-1,k)
    END DO
  END DO

  DO i=1,nx
    DO j=1,ny
      ptbar (i,j,nz) = ptbar (i,j,nz-1)
      pbar  (i,j,nz) = pbar  (i,j,nz-1)
      ppi   (i,j,nz) = ppi   (i,j,nz-1)
      qvbar (i,j,nz) = qvbar (i,j,nz-1)
      csndsq(i,j,nz) = csndsq(i,j,nz-1)
      rhostr(i,j,nz) = rhostr(i,j,nz-1)
    END DO
  END DO

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        ptbari(i,j,k)   = 1.0/ptbar(i,j,k)
        pbari(i,j,k)    = 1.0/pbar(i,j,k)
        rhostri(i,j,k)  = 1.0/rhostr(i,j,k)
      END DO
    END DO
  END DO

  IF ( sfcphy > 0 ) THEN

    CALL initsfc(nx,ny,nz,nzsoil,nstyps,                                &
                 zpsoil,                                                &
                 pbar,pprt(1,1,1,1),                                    &
                 ptbar,ptprt(1,1,1,1),                                  &
                 qvbar,qv(1,1,1,1),                                     &
                 soiltyp,stypfrct, vegtyp, lai,roufns,veg,tem1,         &
                 tsoil,qsoil,wetcanp,snowdpth,qvsfc,tem1soil)
!
! 07/05/2002 Zuwen He
!
! The following code is added for version compatibility.
! Please refer to readsoil subroutine for more detail.
!
! Before fmtver500 there is no zpsoil specified. Therefore
! error may occur when read in soilvar file written prior
! to version500.
!
! We have hard-coded zpsoil in "readsoil", so that dzsoil
! is always 1m for all the old soilvar files. However,
! there is no terrain data passed to readsoil, and zpsoil
! is not correctly initialized.
!
! The following code fixes the problem.
!
    DO k=1, nzsoil
      DO j=1, ny
        DO i=1, nx
          zpsoil(i,j,k)=hterain(i,j)-(zpsoil(i,j,1)-zpsoil(i,j,k))
        END DO
      END DO
    END DO

  END IF

  RETURN
END SUBROUTINE initgrdvar
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE INITDVR                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE initdvr(nx,ny,nz,nts,                                        &
           ubar,vbar,ptbar,pbar,rhostr,qvbar,                           &
           x,y,z,zp,hterain, j1,j2,j3,                                  &
           u,v,w,ptprt,pprt,qv,qscalar,                                 &
           ptcumsrc,qcumsrc,raing,rainc,prcrate,tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Initialize the model time dependent variables.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/1991.
!
!  MODIFICATION HISTORY:
!
!  5/25/92 (M. Xue)
!  Added full documentation.
!
!  5/03/92 (M. Xue)
!  Further documentation.
!
!  4/10/93 (M. Xue & Hao Jin)
!  Add the terrain.
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!  01/28/95 (G. Bassett)
!  Added pt0opt=5 (soup can shaped perturbation to ptprt).
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nts      Number of time levels to be initialized.
!
!    ubar     Base state x-velocity component (m/s)
!    vbar     Base state y-velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!    qvbar    Base state water vapor specific humidity (kg/kg).
!
!    x        x-coordinate of grid points in computational space (m)
!    y        y-coordinate of grid points in computational space (m)
!    z        z-coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space(m)
!    hterain  Terrain height (m)
!
!    j1       Coordinate transformation Jacobian -d(zp)/dx
!    j2       Coordinate transformation Jacobian -d(zp)/dy
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!
!  OUTPUT:
!
!    u        x-component of velocity at all time levels (m/s).
!    v        y-component of velocity at all time levels (m/s).
!    w        z-component of velocity at all time levels (m/s).
!    ptprt    Perturbation potential temperature at all time levels
!             (K)
!    pprt     Perturbation pressure at all time levels (Pascal)
!    qv       Water vapor specific humidity at all time levels
!             (kg/kg)
!    qc       Cloud water mixing ratio at all time levels (kg/kg)
!    qr       Rainwater mixing ratio at all time levels (kg/kg)
!    qi       Cloud ice mixing ratio at all time levels (kg/kg)
!    qs       Snow mixing ratio at all time levels (kg/kg)
!    qh       Hail mixing ratio at all time levels (kg/kg)
!
!    ptcumsrc Source term in pt-equation due to cumulus parameterization
!    qcumsrc Source term in water equations due to cumulus parameterization
!    raing    Grid supersaturation rain
!    rainc    Cumulus convective rain
!    prcrate  Precipitation ratesrain
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'phycst.inc'
  INCLUDE 'bndry.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz          ! The number of grid points in 3
                               ! directions
  INTEGER :: nts               ! Number of time levels to be initialized.
  INTEGER :: tpast             ! Index of time level for the past time.
  INTEGER :: tpresent          ! Index of time level for the present
                               ! time.
  INTEGER :: tfuture           ! Index of time level for the future
                               ! time.

  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal).
  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific
                               ! humidity(kg/kg)

  REAL :: x     (nx)           ! x-coord. of the physical and compu-
                               ! tational grid. Defined at u-point.
  REAL :: y     (ny)           ! y-coord. of the physical and compu-
                               ! tational grid. Defined at v-point.
  REAL :: z     (nz)           ! z-coord. of the computational grid.
                               ! Defined at w-point on the staggered
                               ! grid.
  REAL :: zp    (nx,ny,nz)     ! Physical height coordinate defined at
                               ! w-point of the staggered grid.
  REAL :: hterain(nx,ny)       ! Terrain height (m).

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! -d(zp)/dx.
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! -d(zp)/dy.
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! d(zp)/dz.

  REAL :: u     (nx,ny,nz,nts) ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nts) ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz,nts) ! Total w-velocity (m/s)
  REAL :: ptprt (nx,ny,nz,nts) ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz,nts) ! Perturbation pressure (Pascal)

  REAL :: qv    (nx,ny,nz,nts) ! Water vapor specific humidity (kg/kg)
  REAL :: qscalar(nx,ny,nz,nts,nscalar)

  REAL :: ptcumsrc(nx,ny,nz)   ! Source term in pt-equation due
                               ! to cumulus parameterization
  REAL :: qcumsrc(nx,ny,nz,5)  ! Source term in water equations due
                               ! to cumulus parameterization:
                               ! qcumsrc(1,1,1,1) for qv equation
                               ! qcumsrc(1,1,1,2) for qc equation
                               ! qcumsrc(1,1,1,3) for qr equation
                               ! qcumsrc(1,1,1,4) for qi equation
                               ! qcumsrc(1,1,1,5) for qs equation
  REAL :: raing (nx,ny)        ! Grid supersaturation rain
  REAL :: rainc (nx,ny)        ! Cumulus convective rain
  REAL :: prcrate(nx,ny,4)     ! precipitation rate (kg/(m**2*s))
                               ! prcrate(1,1,1) = total precip. rate
                               ! prcrate(1,1,2) = grid scale precip. rate
                               ! prcrate(1,1,3) = cumulus precip. rate
                               ! prcrate(1,1,4) = microphysics precip. rate

  REAL :: tem1(nx,ny,nz)       ! Temporary work array

!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  REAL :: xs, ys, zs
  REAL :: us, vs, ws, rhobar
  REAL :: radnd , pi,pi2,pi4
  INTEGER :: i,j,k, n, ip
  INTEGER :: iseed,ibgn,iend,jbgn,jend,kbgn,kend
  INTEGER :: ebc1,wbc1,nbc1,sbc1

  REAL :: amplitud
  REAL :: knumx,lnumy,mnumz
  REAL :: lnthx,lnthy,lnthz
  REAL :: lambda,lambdah,lambda2

  INTEGER           :: nxlg, nylg
  REAL, ALLOCATABLE :: temlg1(:,:,:)
  REAL, ALLOCATABLE :: temlg2(:,:,:)
  INTEGER           :: istatus
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF( nts > 1 ) THEN
    tpast=1
    tpresent=2
    tfuture=3
  ELSE
    tpast=1
    tpresent=1
    tfuture=1
  END IF
!
!-----------------------------------------------------------------------
!
!  Specify the initial potential temperature perturbation.
!
!-----------------------------------------------------------------------
!

  IF( pt0opt == 1 .OR. pt0opt == 6 ) THEN  ! Bubble shaped perturbation

!
!-----------------------------------------------------------------------
!
!  Define a potential temperature perturbation for a bubble-shaped
!  disturbance.
!
!  DTD: updated to allow for multiple thermal bubbles
!
!-----------------------------------------------------------------------
!

    pi2 = 1.5707963267949

    ip = 1
    DO WHILE(ptpert0(ip) /= 0.0)

      DO k= 1,nz-1
        DO j= 1,ny-1
          DO i= 1,nx-1

            xs = (x(i)+x(i+1))*0.5
            ys = (y(j)+y(j+1))*0.5
!            xs = (x(i)+x(i+1))*0.5-x(1)
!            ys = (y(j)+y(j+1))*0.5-y(1)
                            !wdt multiple bubbles for timing tests
            zs = (zp(i,j,k)+zp(i,j,k+1))*0.5

            IF( pt0rady(ip) < 0.0 .OR. runmod == 2 ) THEN
                                         ! 2-d bubble in x-z plane.

              radnd = SQRT( ((xs-pt0ctrx(ip))/pt0radx(ip))**2                   &
                         + ((zs-pt0ctrz(ip))/pt0radz(ip))**2 )

            ELSE IF( pt0radx(ip) < 0.0 .OR. runmod == 3 ) THEN
                                         ! 2-d bubble in y-z plane.

              radnd = SQRT( ((ys-pt0ctry(ip))/pt0rady(ip))**2                   &
                         + ((zs-pt0ctrz(ip))/pt0radz(ip))**2 )

            ELSE                         ! 3-d bubble

              radnd = SQRT( ((xs-pt0ctrx(ip))/pt0radx(ip))**2 +                 &
                           ((ys-pt0ctry(ip))/pt0rady(ip))**2 +                  &
                           ((zs-pt0ctrz(ip))/pt0radz(ip))**2 )
            END IF

            IF (ip == 1) THEN ! For first bubble, zero out rest of ptprt field
              ptprt(i,j,k,1) = 0.0
            END IF
            IF(radnd < 1.0) THEN
              ! Note: in case bubbles overlap, the perturbations are simply summed
              ptprt(i,j,k,1) = ptprt(i,j,k,1)+ptpert0(ip)*(COS(pi2*radnd )**2)
              ! DTD: temporarily changed bubble function to be consistent with Straka et al. 1993
!              ptprt(i,j,k,1) = ptprt(i,j,k,1)+ptpert0(ip)*(COS(pi2*2.0*radnd )+1.0)/2.0
            END IF

          END DO
        END DO
      END DO

      IF(pt0opt == 6) THEN ! Perturbation speficied in T'.

        DO k= 1,nz-1
          DO j= 1,ny-1
            DO i= 1,nx-1
              ptprt(i,j,k,1) = ptprt(i,j,k,1)*(p0/pbar(i,j,k))**rddcp
            END DO
          END DO
        END DO

      END IF

      ip = ip + 1
    END DO ! While there are still thermal bubbles


  ELSE IF( pt0opt == 2 .OR. pt0opt == 3 ) THEN
                                          ! Random field
!
!-----------------------------------------------------------------------
!
!  Define a potential temperature perturbation by a random function.
!  This ensures that the horizontal average of the perturbation in a
!  specified domain is zero.
!
!  Fill an array with zeros.
!
!-----------------------------------------------------------------------
!
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          ptprt(i,j,k,1) = 0.0
        END DO
      END DO
    END DO

    nxlg = ((nx-3)*nproc_x+3)   ! will be nx if nproc_x = 1
                   ! which is the case for serial run, see initpara
    nylg = ((ny-3)*nproc_y+3)   ! will be ny if nproc_y = 1
                   ! which is the case for serial run, see initpara
!
!-----------------------------------------------------------------------
!
!    The following parameters define the portion of domain to be
!    initialized with a random potential temperture perturbation.
!    Users can modify them to fit their needs.
!
!    NOTE: if this is an MPI run, ibgn,iend and jbgn,jend should
!    be defined over the entire domain, and not just for one
!    processor!
!
!-----------------------------------------------------------------------
!
    ibgn = 1
    iend = nxlg - 1      ! will be nx-1 if nproc_x = 1
    jbgn = 1
    jend = nylg - 1      ! will be ny-1 if nproc_y = 1
    kbgn = 1
    kend = nz-1

    iseed = -100

    DO k= kbgn,kend

      CALL ranary(nx,ny,ibgn,iend,jbgn,jend,                        &
                  iseed,ptpert0(1),ptprt(1,1,k,1))

    END DO

    IF( pt0opt == 3 ) THEN ! Symmetric random perturbation
!
!-----------------------------------------------------------------------
!
!  Create a random perturbation field symmetric about both the x and y
!  axes.
!
!-----------------------------------------------------------------------

      ALLOCATE ( temlg1(nxlg,nylg,nz), STAT = istatus )
      CALL check_alloc_status(istatus, "initdvr:temlg1")
      ALLOCATE ( temlg2(nxlg,nylg,nz), STAT = istatus )
      CALL check_alloc_status(istatus, "initdvr:temlg2")

      CALL mpimerge3d(ptprt(:,:,:,1), nx, ny, nz, temlg1)

      IF (myproc ==0) THEN
        DO k=1,nz-1
          DO j=1,nylg-1
            DO i=1,nxlg/2
              temlg1(i,j,k) = temlg1(nxlg-i,j,k)
            END DO
          END DO
        END DO

        DO k=1,nz-1
          DO i=1,nxlg-1
            DO j=1,nylg/2
              temlg1(i,j,k) = temlg1(i,nylg-j,k)
            END DO
          END DO
        END DO

        IF( nxlg == nylg ) THEN
          DO k=1,nz-1
            DO j=1,nylg-1
              DO i=1,nxlg-1
                temlg2(i,j,k) = (temlg1(i,j,k)+temlg1(j,i,k))*0.5
              END DO
            END DO
          END DO

          DO k=1,nz-1
            DO i=1,nxlg-1
              DO j=1,nylg-1
                temlg1(i,j,k) = temlg2(i,j,k)
              END DO
            END DO
          END DO
        END IF

      END IF  ! myproc == 0

      CALL mpisplit3d(temlg1,nx,ny,nz,ptprt(:,:,:,1))

      DEALLOCATE(temlg1, temlg2)

    END IF

  ELSE IF( pt0opt == 4 ) THEN       ! Bubble given in Skamarock and
                                    ! Klemp (1994)
    pi2 = 1.5707963267949

    DO k=1,nz-1
      DO i=1,nx-1
        DO j=1,ny-1
          xs = (x(i)+x(i+1))*0.5
          zs = (zp(i,j,k)+zp(i,j,k+1))*0.5
          ptprt(i,j,k,1) = ptpert0(1)                                      &
              *SIN(pi2*2*zs/((nz-3)*dz))/(1+((xs-pt0ctrx(1))/pt0radx(1))**2)
        END DO
      END DO
    END DO

  ELSE IF( pt0opt == 5 ) THEN       ! Soup can shaped perturbation

    DO k= 1,nz-1
      DO j= 1,ny-1
        DO i= 1,nx-1

          xs = (x(i)+x(i+1))*0.5
          ys = (y(j)+y(j+1))*0.5
          zs = (zp(i,j,k)+zp(i,j,k+1))*0.5
          IF ( ABS(zs-pt0ctrz(1))/pt0radz(1) < 1 ) THEN

            IF( pt0rady(1) < 0.0 .OR. runmod == 2 ) THEN
                                         ! 2-d bubble in x-z plane.

              radnd = ABS(xs-pt0ctrx(1))/pt0radx(1)

            ELSE IF( pt0radx(1) < 0.0 .OR. runmod == 3 ) THEN
                                         ! 2-d bubble in y-z plane.

              radnd = ABS(ys-pt0ctry(1))/pt0rady(1)

            ELSE                         ! 3-d bubble

              radnd = SQRT( ((xs-pt0ctrx(1))/pt0radx(1))**2 +                 &
                            ((ys-pt0ctry(1))/pt0rady(1))**2 )

            END IF

            IF(radnd >= 1.0) THEN
              ptprt(i,j,k,1) = 0.0
            ELSE
              ptprt(i,j,k,1) = ptpert0(1)
            END IF

          ELSE

            ptprt(i,j,k,1) = 0.0

          END IF

        END DO
      END DO
    END DO

  END IF

  ebc1=0
  wbc1=0
  sbc1=0
  nbc1=0

  IF( ebc == 1 .OR.ebc == 2 .OR. ebc == 3 )  ebc1=ebc
  IF( wbc == 1 .OR.wbc == 2 .OR. wbc == 3 )  wbc1=wbc
  IF( sbc == 1 .OR.sbc == 2 .OR. sbc == 3 )  sbc1=sbc
  IF( nbc == 1 .OR.nbc == 2 .OR. nbc == 3 )  nbc1=nbc

  IF (mp_opt > 0) THEN
    CALL acct_interrupt(mp_acct)
    CALL mpsendrecv2dew(ptprt,nx,ny,nz,ebc1,wbc1,0,tem1)
    CALL mpsendrecv2dns(ptprt,nx,ny,nz,nbc1,sbc1,0,tem1)
  END IF
  CALL acct_interrupt(bc_acct)
  CALL bcsclr(nx,ny,nz,dtbig,                                           &
              ptprt(1,1,1,1),ptprt(1,1,1,1),                            &
              ptprt(1,1,1,1),tem1,tem1,tem1,tem1,                       &
              ebc1,wbc1,nbc1,sbc1,tbc,bbc,                              &
              ebc_global,wbc_global,nbc_global,sbc_global)
  CALL acct_stop_inter

  DO n=1,nts
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          ptprt(i,j,k,n)=ptprt(i,j,k,1)
          pprt(i,j,k,n)=0.0
        END DO
      END DO
    END DO
  END DO

  DO n=1,nts
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx
          u(i,j,k,n)=ubar(i,j,k)
        END DO
      END DO
    END DO
  END DO

  DO n=1,nts
    DO k=1,nz-1
      DO j=1,ny
        DO i=1,nx-1
          v(i,j,k,n)=vbar(i,j,k)
        END DO
      END DO
    END DO
  END DO

  DO n=1,nts
    DO k=1,nz
      DO j=1,ny-1
        DO i=1,nx-1
          w(i,j,k,n)=0.0
        END DO
      END DO
    END DO
  END DO

  DO n=1,nts
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          qv(i,j,k,n)= qvbar(i,j,k)
        END DO
      END DO
    END DO
  END DO

  qscalar(:,:,:,:,:) = 0.0

!
!-----------------------------------------------------------------------
! TINA
!  Specify the initial concentration field. (granvold)
!
!-----------------------------------------------------------------------
  IF (ccin == 1 .AND. cpoint < 0 .AND. P_CC > 0) THEN
    CALL ccinit(nx,ny,nz,3,x,y,z,zp,qscalar(:,:,:,:,P_CC),tem1)
    IF (myproc == 0) print *,'Initializing cc as bubble'
  END IF
!
!TINA -end granvold modifications
!
!michi
!TINA modified by granvold - cc initial is set above
  IF (P_CC > 0) THEN
    DO n=1,nts
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            qscalar(i,j,k,n,P_CC)= qscalar(i,j,k,1,P_CC)
          END DO
        END DO
      END DO
    END DO
  END IF
!michi

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        ptcumsrc(i,j,k)= 0.0
        qcumsrc(i,j,k,1)= 0.0
        qcumsrc(i,j,k,2)= 0.0
        qcumsrc(i,j,k,3)= 0.0
        qcumsrc(i,j,k,4)= 0.0
        qcumsrc(i,j,k,5)= 0.0
      END DO
    END DO
  END DO

  DO j=1,ny-1
    DO i=1,nx-1
      raing(i,j)= 0.0
      rainc(i,j)= 0.0
      prcrate(i,j,1)= 0.0
      prcrate(i,j,2)= 0.0
      prcrate(i,j,3)= 0.0
      prcrate(i,j,4)= 0.0
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  The following setup is to overwrite the total u, v, w, pprt,
!  ptprt, and qv for Beltrami flow initial conditions.
!
!-----------------------------------------------------------------------
!
  IF ( pt0opt == 7 ) THEN

    pi = 3.1415926535898
    pi2 = 2*pi
    pi4 = 4*pi

    amplitud = 2.0                      ! amplitude A=2 m/s
    tmixopt  = 1                        ! constant viscosity option
    tmixcst  = 1.0                      ! viscosity=1 m**2/s
    tmixvert = 1.0                      ! compute all mixing terms

    wbc = 2                             ! reset boundary conditions
    ebc = 2                             ! to periodical condition
    nbc = 2
    sbc = 2
    tbc = 2
    bbc = 2

    lnthx = dx*(nx-3)                   ! length in x
    lnthy = dy*(ny-3)                   ! length in y
    lnthz = dz*(nz-3)                   ! length in z

    knumx = pi4/lnthx                   ! wave number in x-dir
    lnumy = pi2/lnthy                   ! wave number in y-dir
    mnumz = pi2/lnthz                   ! wave number in z-dir

    lambdah = knumx*knumx + lnumy*lnumy
    lambda2 = lambdah + mnumz*mnumz
    lambda  = SQRT( lambda2 )

!    print *, ' amplitude = ',amplitud
    PRINT *, ' lnthx   = ',lnthx,                                       &
             ' lnthy   = ',lnthy,                                       &
             ' lnthz   = ',lnthz
    PRINT *, ' knumx   = ',knumx,                                       &
             ' lnumy   = ',lnumy,                                       &
             ' mnumz   = ',mnumz
    PRINT *, ' lambda1 = ',lambdah,                                     &
             ' lambda2 = ',lambda2,                                     &
             ' lambda  = ',lambda

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx
          ys = 0.5*(y(j+1)+y(j))
          zs = 0.5*(z(k+1)+z(k))
          u(i,j,k,1) = -amplitud/lambdah                                &
                     *(lambda*lnumy                                     &
                      *COS(knumx*x(i))*SIN(lnumy*ys)*SIN(mnumz*zs)      &
                      +mnumz*knumx                                      &
                      *SIN(knumx*x(i))*COS(lnumy*ys)*COS(mnumz*zs))
        END DO
      END DO
    END DO

    DO k=1,nz-1
      DO j=1,ny
        DO i=1,nx-1
          xs = 0.5*(x(i+1)+x(i))
          zs = 0.5*(z(k+1)+z(k))
          v(i,j,k,1)= amplitud/lambdah                                  &
                    * (lambda*knumx                                     &
                      *SIN(knumx*xs)*COS(lnumy*y(j))*SIN(mnumz*zs)      &
                      -mnumz*lnumy                                      &
                      *COS(knumx*xs)*SIN(lnumy*y(j))*COS(mnumz*zs))
        END DO
      END DO
    END DO

    DO k=1,nz
      DO j=1,ny-1
        DO i=1,nx-1
          xs = 0.5*(x(i+1)+x(i))
          ys = 0.5*(y(j+1)+y(j))
          w(i,j,k,1)=amplitud                                           &
                    *COS(knumx*xs)*COS(lnumy*ys)*SIN(mnumz*z(k))
        END DO
      END DO
    END DO

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          zs = 0.5*(z(k+1)+z(k))

          us = 0.5*(u(i+1,j,k,1)+u(i,j,k,1))
          vs = 0.5*(v(i,j+1,k,1)+v(i,j,k,1))
          ws = 0.5*(w(i,j,k+1,1)+w(i,j,k,1))

          rhobar = rhostr(i,j,k)/j3(i,j,k)
          pprt(i,j,k,1) = p0-rhobar*(0.5*(us*us+vs*vs+ws*ws)+g*zs)      &
                        - pbar(i,j,k)
        END DO
      END DO
    END DO

    DO n=1,nts
      DO k=1,nz
        DO j=1,ny
          DO i=1,nx
            u    (i,j,k,n) = u(i,j,k,1)
            v    (i,j,k,n) = v(i,j,k,1)
            w    (i,j,k,n) = w(i,j,k,1)
            pprt (i,j,k,n) = pprt(i,j,k,1)
            ptprt(i,j,k,n) = 0.0
            qv   (i,j,k,n) = 0.0
            IF (P_CC > 0) qscalar(:,:,:,n,P_CC) = 0.0
          END DO
        END DO
      END DO
    END DO

  END IF

  RETURN
END SUBROUTINE initdvr

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE BUBBLE                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE bubble(nx,ny,nz,nts,                                        &
           ubar,vbar,ptbar,pbar,rhostr,qvbar,                           &
           x,y,z,zp,hterain, j1,j2,j3,                                  &
           u,v,w,ptprt,pprt,qv,qscalar,                                 &
           tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Apply a thermal bubble at any time step at a given location,
!  while leaving the rest of the fields alone.
!  Based on subroutine INITDVR, and only works for pt0opt = 1,6 for now
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Dawson
!  9/5/2007.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nts      Number of time levels to be initialized.
!
!    ubar     Base state x-velocity component (m/s)
!    vbar     Base state y-velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!    qvbar    Base state water vapor specific humidity (kg/kg).
!
!    x        x-coordinate of grid points in computational space (m)
!    y        y-coordinate of grid points in computational space (m)
!    z        z-coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space(m)
!    hterain  Terrain height (m)
!
!    j1       Coordinate transformation Jacobian -d(zp)/dx
!    j2       Coordinate transformation Jacobian -d(zp)/dy
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!
!  OUTPUT:
!
!    u        x-component of velocity at all time levels (m/s).
!    v        y-component of velocity at all time levels (m/s).
!    w        z-component of velocity at all time levels (m/s).
!    ptprt    Perturbation potential temperature at all time levels
!             (K)
!    pprt     Perturbation pressure at all time levels (Pascal)
!    qv       Water vapor specific humidity at all time levels
!             (kg/kg)
!    qc       Cloud water mixing ratio at all time levels (kg/kg)
!    qr       Rainwater mixing ratio at all time levels (kg/kg)
!    qi       Cloud ice mixing ratio at all time levels (kg/kg)
!    qs       Snow mixing ratio at all time levels (kg/kg)
!    qh       Hail mixing ratio at all time levels (kg/kg)
!
!    ptcumsrc Source term in pt-equation due to cumulus parameterization
!    qcumsrc Source term in water equations due to cumulus parameterization
!    raing    Grid supersaturation rain
!    rainc    Cumulus convective rain
!    prcrate  Precipitation ratesrain
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'phycst.inc'
  INCLUDE 'bndry.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz          ! The number of grid points in 3
                               ! directions
  INTEGER :: nts               ! Number of time levels to be initialized.
  INTEGER :: tpast             ! Index of time level for the past time.
  INTEGER :: tpresent          ! Index of time level for the present
                               ! time.
  INTEGER :: tfuture           ! Index of time level for the future
                               ! time.

  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal).
  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific
                               ! humidity(kg/kg)

  REAL :: x     (nx)           ! x-coord. of the physical and compu-
                               ! tational grid. Defined at u-point.
  REAL :: y     (ny)           ! y-coord. of the physical and compu-
                               ! tational grid. Defined at v-point.
  REAL :: z     (nz)           ! z-coord. of the computational grid.
                               ! Defined at w-point on the staggered
                               ! grid.
  REAL :: zp    (nx,ny,nz)     ! Physical height coordinate defined at
                               ! w-point of the staggered grid.
  REAL :: hterain(nx,ny)       ! Terrain height (m).

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! -d(zp)/dx.
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! -d(zp)/dy.
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! d(zp)/dz.

  REAL :: u     (nx,ny,nz,nts) ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nts) ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz,nts) ! Total w-velocity (m/s)
  REAL :: ptprt (nx,ny,nz,nts) ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz,nts) ! Perturbation pressure (Pascal)

  REAL :: qv    (nx,ny,nz,nts) ! Water vapor specific humidity (kg/kg)
  REAL :: qscalar(nx,ny,nz,nts,nscalar)

  REAL :: ptcumsrc(nx,ny,nz)   ! Source term in pt-equation due
                               ! to cumulus parameterization
  REAL :: qcumsrc(nx,ny,nz,5)  ! Source term in water equations due
                               ! to cumulus parameterization:
                               ! qcumsrc(1,1,1,1) for qv equation
                               ! qcumsrc(1,1,1,2) for qc equation
                               ! qcumsrc(1,1,1,3) for qr equation
                               ! qcumsrc(1,1,1,4) for qi equation
                               ! qcumsrc(1,1,1,5) for qs equation
  REAL :: raing (nx,ny)        ! Grid supersaturation rain
  REAL :: rainc (nx,ny)        ! Cumulus convective rain
  REAL :: prcrate(nx,ny,4)     ! precipitation rate (kg/(m**2*s))
                               ! prcrate(1,1,1) = total precip. rate
                               ! prcrate(1,1,2) = grid scale precip. rate
                               ! prcrate(1,1,3) = cumulus precip. rate
                               ! prcrate(1,1,4) = microphysics precip. rate

  REAL :: tem1(nx,ny,nz)       ! Temporary work array

!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  REAL :: xs, ys, zs
  REAL :: us, vs, ws, rhobar
  REAL :: radnd , pi,pi2,pi4
  INTEGER :: i,j,k, n, ip
  INTEGER :: iseed,ibgn,iend,jbgn,jend,kbgn,kend
  INTEGER :: ebc1,wbc1,nbc1,sbc1

  REAL :: amplitud
  REAL :: knumx,lnumy,mnumz
  REAL :: lnthx,lnthy,lnthz
  REAL :: lambda,lambdah,lambda2

  INTEGER           :: nxlg, nylg
  REAL, ALLOCATABLE :: temlg1(:,:,:)
  REAL, ALLOCATABLE :: temlg2(:,:,:)
  INTEGER           :: istatus
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF( nts > 1 ) THEN
    tpast=1
    tpresent=2
    tfuture=3
  ELSE
    tpast=1
    tpresent=1
    tfuture=1
  END IF

!  print*,'BUBBLE: myproc,nx,ny',myproc,nx,ny
!  print*,'BUBBLE: nts =',nts
!  print*,'BUBBLE: pt0opt,ptpert0(1)',pt0opt,ptpert0(1)
!
!-----------------------------------------------------------------------
!
!  Specify the initial potential temperature perturbation.
!
!-----------------------------------------------------------------------
!

  IF( pt0opt == 1 .OR. pt0opt == 6 ) THEN  ! Bubble shaped perturbation

!
!-----------------------------------------------------------------------
!
!  Define a potential temperature perturbation for a bubble-shaped
!  disturbance.
!
!-----------------------------------------------------------------------
!

    pi2 = 1.5707963267949
    ip = 1

    DO WHILE(ptpert0(ip) /= 0.0)

!      print*,'BUBBLE: ip,ptpert0(ip)',ip,ptpert0(ip)

      DO k= 1,nz-1
        DO j= 1,ny-1
          DO i= 1,nx-1

            xs = (x(i)+x(i+1))*0.5
            ys = (y(j)+y(j+1))*0.5
!            xs = (x(i)+x(i+1))*0.5-x(1)
!            ys = (y(j)+y(j+1))*0.5-y(1)
                            !wdt multiple bubbles for timing tests
            zs = (zp(i,j,k)+zp(i,j,k+1))*0.5

            IF( pt0rady(ip) < 0.0 .OR. runmod == 2 ) THEN
                                         ! 2-d bubble in x-z plane.

              radnd = SQRT( ((xs-pt0ctrx(ip))/pt0radx(ip))**2                   &
                         + ((zs-pt0ctrz(ip))/pt0radz(ip))**2 )

            ELSE IF( pt0radx(ip) < 0.0 .OR. runmod == 3 ) THEN
                                         ! 2-d bubble in y-z plane.

              radnd = SQRT( ((ys-pt0ctry(ip))/pt0rady(ip))**2                   &
                         + ((zs-pt0ctrz(ip))/pt0radz(ip))**2 )

            ELSE                         ! 3-d bubble

              radnd = SQRT( ((xs-pt0ctrx(ip))/pt0radx(ip))**2 +                 &
                           ((ys-pt0ctry(ip))/pt0rady(ip))**2 +                  &
                           ((zs-pt0ctrz(ip))/pt0radz(ip))**2 )
            END IF

            IF(radnd < 1.0) THEN
              ptprt(i,j,k,1) = MAX(ptprt(i,k,k,1),ptpert0(ip)*(COS(pi2*radnd )**2))
              IF(nts > 1) THEN
                ptprt(i,j,k,2) = MAX(ptprt(i,j,k,2),ptpert0(ip)*(COS(pi2*radnd )**2))
                ptprt(i,j,k,3) = MAX(ptprt(i,k,k,3),ptpert0(ip)*(COS(pi2*radnd )**2))
              END IF
            END IF

          END DO
        END DO
      END DO

      IF(pt0opt == 6) THEN ! Perturbation speficied in T'.

        DO k= 1,nz-1
          DO j= 1,ny-1
            DO i= 1,nx-1
              ptprt(i,j,k,1) = ptprt(i,j,k,1)*(p0/pbar(i,j,k))**rddcp
              IF(nts > 1) THEN
                ptprt(i,j,k,2) = ptprt(i,j,k,1)
                ptprt(i,j,k,3) = ptprt(i,k,k,1)
              END IF
            END DO
          END DO
        END DO

      END IF

      ip = ip + 1
    END DO

  END IF

  ebc1=0
  wbc1=0
  sbc1=0
  nbc1=0

  IF( ebc == 1 .OR.ebc == 2 .OR. ebc == 3 )  ebc1=ebc
  IF( wbc == 1 .OR.wbc == 2 .OR. wbc == 3 )  wbc1=wbc
  IF( sbc == 1 .OR.sbc == 2 .OR. sbc == 3 )  sbc1=sbc
  IF( nbc == 1 .OR.nbc == 2 .OR. nbc == 3 )  nbc1=nbc

  IF (mp_opt > 0) THEN
    CALL acct_interrupt(mp_acct)
    CALL mpsendrecv2dew(ptprt,nx,ny,nz,ebc1,wbc1,0,tem1)
    CALL mpsendrecv2dns(ptprt,nx,ny,nz,nbc1,sbc1,0,tem1)
  END IF
  CALL acct_interrupt(bc_acct)
  CALL bcsclr(nx,ny,nz,dtbig,                                           &
              ptprt(1,1,1,1),ptprt(1,1,1,1),                            &
              ptprt(1,1,1,1),tem1,tem1,tem1,tem1,                       &
              ebc1,wbc1,nbc1,sbc1,tbc,bbc,                              &
              ebc_global,wbc_global,nbc_global,sbc_global)
  CALL acct_stop_inter

  RETURN
END SUBROUTINE bubble

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE EXTINIT                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE extinit(nx,ny,nz,nzsoil,nts,nstyps,                          &
           u,v,w,ptprt,pprt,qv,qscalar,tke,                             &
           ubar,vbar,ptbar,pbar,rhostr,qvbar,                           &
           x,y,z,zp,zpsoil,hterain,j1,j2,j3,j3soil,                     &
           soiltyp,stypfrct,vegtyp,lai,roufns,veg,                      &
           tsoil,qsoil,wetcanp,snowdpth,qvsfc,                          &
           ptcumsrc,qcumsrc,raing,rainc,prcrate,                        &
           radfrc,radsw,rnflx,radswnet,radlwin,                         &
           usflx,vsflx,ptsflx,qvsflx,                                   &
           uprt,vprt,qvprt,kmh,kmv,wbar,                                &
           tem6,tem7,tem8)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read in the initial fields from an externally supplied initial
!  data file.
!
!  These fields include u,v,w,ptprt,pprt,qv,qc,qr,qi,qs,qh
!  at time level tpresent, and the base state variables
!  ubar,vbar,ptbar,pbar,rhostr,qvbar.
!
!  Fields at tpast and tfuture are set to their values at tpresent.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/7/1992
!
!  MODIFICATION HISTORY:
!
!  3/25/94 (G. Bassett, M. Xue)
!  Modified to read in regular binary history dumps.
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!  01/14/1995 (M. Xue)
!  Corrected where jacob was called. It should be called
!  before do loop 90.
!
!  03/29/1997 (M. Xue)
!  Modification made so that when values of mapproj,trulat1,trulat2,
!  trulon,ctrlat,ctrlon in the input data do not match those in
!  input file, the values in the data are used instead of those
!  in input, as was done before. Warning messages will be printed.
!
!  12/09/1998 (Donghai Wang)
!  Added the snow cover.
!
!  03/19/2002 (Keith Brewster)
!  Corrected print statement related to mis-match in data and user times.
!  05/18/20020 (Dan Weber)
!  Added new soil model variables.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nzsoil   Number of grid points in the soil model in the -z-direction

!    nts      Number of time levels to be initialized.
!
!  OUTPUT:
!
!    u        x component of velocity at times tpast and tpresent
!             (m/s)
!    v        y component of velocity at times tpast and tpresent
!             (m/s)
!    w        Vertical component of Cartesian velocity at times
!             tpast and tpresent (m/s)
!    ptprt    Perturbation potential temperature at times tpast and
!             tpresent (K)
!    pprt     Perturbation pressure at times tpast and tpresent
!             (Pascal)
!
!    qv       Water vapor specific humidity at times tpast and
!             tpresent (kg/kg)
!    qc       Cloud water mixing ratio at times tpast and tpresent
!             (kg/kg)
!    qr       Rainwater mixing ratio at times tpast and tpresent
!             (kg/kg)
!    qi       Cloud ice mixing ratio at times tpast and tpresent
!             (kg/kg)
!    qs       Snow mixing ratio at times tpast and tpresent (kg/kg)
!    qh       Hail mixing ratio at times tpast and tpresent (kg/kg)
!    tke      Turbulence kinetic energy (m**2/s)
!
!    ubar     Base state zonal velocity component (m/s)
!    vbar     Base state meridional velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!    qvbar    Base state water vapor specific humidity (kg/kg)
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space(m)
!    zpsoil   Vertical coordinate of grid points in the soil model
!             in physical space (m).
!    hterain  The height of terrain (m)
!
!    j1       Coordinate transformation Jacobian -d(zp)/dx
!    j2       Coordinate transformation Jacobian -d(zp)/dy
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!    j3soil   Soil coordinate transformation Jacobian  d(zpsoil)/dz
!
!    soiltyp  Soil type
!    vegtyp   Vegetation type
!    lai      Leaf Area Index
!    roufns   Surface roughness
!    veg      Vegetation fraction
!
!    tsoil    Deep soil temperature (K) (in deep 1 m layer)
!    qsoil    Soil moisture (m**3/m**3)
!    wetcanp  Canopy water amount
!    snowdpth Snow depth (m)
!    qvsfc    Effective S.H. at sfc.
!
!    ptcumsrc Source term in pt-equation due to cumulus parameterization
!    qcumsrc Source term in water equations due to cumulus parameterization
!    raing    Grid supersaturation rain
!    rainc    Cumulus convective rain
!    prcrate  Precipitation rates
!
!    radfrc   Radiation forcing (K/s)
!    radsw    Solar radiation reaching the surface
!    rnflx    Net radiation flux absorbed by surface
!    radswnet Net shortwave radiation
!    radlwin  Incoming longwave radiation
!
!    usflx    Surface flux of u-momentum (kg/(m*s**2))
!    vsflx    Surface flux of v-momentum (kg/(m*s**2))
!    ptsflx   Surface heat flux (K*kg/(m**2 * s ))
!    qvsflx   Surface moisture flux of (kg/(m**2 * s))
!
!    tstart   The time when the time integration starts, which is set
!             to the time of the restart data
!
!  WORK ARRAYS:
!
!    tem6     Temporary work array.
!    tem7     Temporary work array.
!    tem8     Temporary work array.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'phycst.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.
  INCLUDE 'indtflg.inc'
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of grid points in the -z-direction

  INTEGER :: nts               ! Number of time levels to be initialized.
  INTEGER :: tpast             ! Index of time level for the past time.
  INTEGER :: tpresent          ! Index of time level for the present
                               ! time.
  INTEGER :: tfuture           ! Index of time level for the future
                               ! time.

  REAL :: u     (nx,ny,nz,nts)  ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nts) ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz,nts) ! Total w-velocity (m/s)
  REAL :: ptprt (nx,ny,nz,nts) ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz,nts) ! Perturbation pressure (Pascal)

  REAL :: qv    (nx,ny,nz,nts) ! Water vapor specific humidity (kg/kg)
  REAL :: qscalar(nx,ny,nz,nts,nscalar)
  REAL :: tke   (nx,ny,nz,nts) ! Turbulence kinetic energy (m**2/s)

  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal)
  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific
                               ! humidity(kg/kg)

  REAL :: x(nx)                ! x-coord. of the physical and compu-
                               ! tational grid. Defined at u-point.
  REAL :: y(ny)                ! y-coord. of the physical and compu-
                               ! tational grid. Defined at v-point.
  REAL :: z(nz)                ! z-coord. of the computational grid.
                               ! Defined at w-point on the staggered
                               ! grid.
  REAL :: zp(nx,ny,nz)         ! Physical height coordinate defined at
                               ! w-point of the staggered grid.
  REAL :: zpsoil(nx,ny,nzsoil) ! The physical height coordinate defined
                               ! at the center of a soil layer(m).

  REAL :: hterain(nx,ny)       ! Terrain height (m).

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! -d(zp)/dx.
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! -d(zp)/dy.
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! d(zp)/dz.
  REAL :: j3soil(nx,ny,nzsoil) ! Coordinate transformation Jacobian
                               ! defined as d( zpsoil )/d( zsoil ).
  REAL :: j3soilinv(nx,ny,nzsoil) ! Inverse of J3soil.

  INTEGER :: nstyps                  ! Number of soil types
  INTEGER :: soiltyp(nx,ny,nstyps)   ! Soil types at grid points
  REAL    :: stypfrct(nx,ny,nstyps)  ! Fraction of soil types
  INTEGER :: vegtyp (nx,ny)          ! Vegetation type
  REAL    :: lai    (nx,ny)          ! Leaf Area Index
  REAL    :: roufns (nx,ny)          ! Surface roughness
  REAL    :: veg    (nx,ny)          ! Vegetation fraction

  REAL :: qvsfc   (nx,ny,       0:nstyps) ! Effective qv at sfc.
  REAL :: tsoil   (nx,ny,nzsoil,0:nstyps) ! Deep soil temperature (K)
                                          ! (in deep 1 m layer)
  REAL :: qsoil   (nx,ny,nzsoil,0:nstyps) ! Soil layer moisture(m**3/m**3)
  REAL :: wetcanp (nx,ny,0:nstyps)      ! Canopy water amount
  REAL :: snowdpth(nx,ny)               ! Snow depth (m)


  REAL :: ptcumsrc(nx,ny,nz)   ! Source term in pt-equation due
                               ! to cumulus parameterization
  REAL :: qcumsrc(nx,ny,nz,5)  ! Source term in water equations due
                               ! to cumulus parameterization:
                               ! qcumsrc(1,1,1,1) for qv equation
                               ! qcumsrc(1,1,1,2) for qc equation
                               ! qcumsrc(1,1,1,3) for qr equation
                               ! qcumsrc(1,1,1,4) for qi equation
                               ! qcumsrc(1,1,1,5) for qs equation
  REAL :: raing(nx,ny)         ! Grid supersaturation rain
  REAL :: rainc(nx,ny)         ! Cumulus convective rain
  REAL :: prcrate(nx,ny,4)     ! precipitation rate (kg/(m**2*s))
                               ! prcrate(1,1,1) = total precip. rate
                               ! prcrate(1,1,2) = grid scale precip. rate
                               ! prcrate(1,1,3) = cumulus precip. rate
                               ! prcrate(1,1,4) = microphysics precip. rate

  REAL :: radfrc(nx,ny,nz)     ! Radiation forcing (K/s)
  REAL :: radsw (nx,ny)        ! Solar radiation reaching the surface
  REAL :: rnflx (nx,ny)        ! Net radiation flux absorbed by surface
  REAL :: radswnet(nx,ny)      ! Net shortwave radiation
  REAL :: radlwin(nx,ny)       ! Incoming longwave radiation

  REAL :: usflx (nx,ny)        ! Surface flux of u-momentum (kg/(m*s**2))
  REAL :: vsflx (nx,ny)        ! Surface flux of v-momentum (kg/(m*s**2))
  REAL :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m*s**2))
  REAL :: qvsflx(nx,ny)        ! Surface moisture flux (kg/(m**2*s))

  REAL :: uprt (nx,ny,nz)      ! Temporary array
  REAL :: vprt (nx,ny,nz)      ! Temporary array
  REAL :: qvprt(nx,ny,nz)      ! Temporary array
  REAL :: kmh  (nx,ny,nz)      ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv  (nx,ny,nz)      ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: wbar (nx,ny,nz)      ! Temporary array
  REAL :: tem6 (nx,ny,nz)      ! Temporary array
  REAL :: tem7 (nx,ny,nz)      ! Temporary array
  REAL :: tem8 (nx,ny,nz)      ! Temporary array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: lengbf, lenfil
  INTEGER :: i, j, k, ireturn, tim, nq
  REAL    :: amax,amin
  CHARACTER (LEN=256) :: runnamesv
  CHARACTER (LEN=80 ) :: cmnt_old(50)  ! String of comments on this model run
  INTEGER             :: nocmnt_old    ! Number of comment lines
  INTEGER :: nchin
  REAL    :: time_tmp
  REAL    :: tstopsv,thisdmpsv,mapprojsv,latitudsv,ctrlatsv,ctrlonsv
  INTEGER :: monthsv,daysv,yearsv,hoursv,minutesv,secondsv
  REAL    :: trulat1sv,trulat2sv,trulonsv
  REAL    :: dxsv,dysv,strhoptsv,dzsv,dzminsv,zrefsfcsv,dlayer1sv,dlayer2sv,  &
             strhtunesv,zflatsv
  REAL    :: p0inv,eps,tvbar

  REAL    :: n0rainsv, n0snowsv, n0hailsv, rhosnowsv, rhohailsv
  REAL    :: ntcloudsv, n0grplsv,rhoicesv,rhogrplsv, alpharainsv
  REAL    :: alphaicesv, alphasnowsv, alphagrplsv, alphahailsv
  INTEGER :: abstsec0,abstsec1
  CHARACTER(LEN=256) :: tmpstr

  REAL, ALLOCATABLE :: temscalar(:,:,:,:)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  p0inv=1./p0

  IF( nts > 1 ) THEN
    tpast=1
    tpresent=2
    tfuture=3
  ELSE
    tpast=1
    tpresent=1
    tfuture=1
  END IF
!
!-----------------------------------------------------------------------
!
!  Read in the initial fields from the input data file.
!
!-----------------------------------------------------------------------
!
  lengbf = 256
  CALL strlnth( inigbf, lengbf )

  lenfil = 256
  CALL strlnth( inifile, lenfil )

  IF (mp_opt > 0 .AND. readsplit(FINDX_H) == 0) THEN
    tmpstr = inigbf
    CALL gtsplitfn(tmpstr,1,1,loc_x,loc_y,1,1,0,0,1,lvldbg,inigbf,ireturn)
    lengbf = LEN_TRIM(inigbf)

    tmpstr = inifile
    CALL gtsplitfn(tmpstr,1,1,loc_x,loc_y,1,1,0,0,1,lvldbg,inifile,ireturn)
    lenfil = LEN_TRIM(inifile)
  END IF

  tim = tpresent

  runnamesv = runname
  tstopsv = tstop
  thisdmpsv = thisdmp
  mapprojsv = mapproj
  trulat1sv = trulat1
  trulat2sv = trulat2
  trulonsv  = trulon
  latitudsv = latitud
  ctrlatsv = ctrlat
  ctrlonsv = ctrlon
  monthsv = month
  daysv = day
  yearsv = year
  hoursv = hour
  minutesv = minute
  secondsv = second
  dxsv = dx
  dysv = dy
  strhoptsv = strhopt
  dzsv = dz
  dzminsv = dzmin
  zrefsfcsv = zrefsfc
  dlayer1sv = dlayer1
  dlayer2sv = dlayer2
  strhtunesv = strhtune
  zflatsv = zflat

  ntcloudsv = ntcloud
  n0rainsv = n0rain
  n0snowsv = n0snow
  n0grplsv = n0grpl
  n0hailsv = n0hail
  rhoicesv = rhoice
  rhosnowsv = rhosnow
  rhogrplsv = rhogrpl
  rhohailsv = rhohail
  alpharainsv = alpharain
  alphaicesv = alphaice
  alphasnowsv = alphasnow
  alphagrplsv = alphagrpl
  alphahailsv = alphahail

  nocmnt_old = nocmnt
  DO i=1,nocmnt_old
    cmnt_old(i)=cmnt(i)
  END DO

  ALLOCATE(temscalar(nx,ny,nz,nscalar), STAT = ireturn)
  temscalar = 0.0

!  blocking inserted for ordering i/o for message passing
  DO i=0,nprocs-1,readstride
    IF(myproc >= i.AND.myproc <= i+readstride-1)THEN

      CALL dtaread(nx,ny,nz,nzsoil,nstyps, inifmt, nchin ,inigbf,lengbf,&
                 inifile,lenfil,time_tmp,                               &
                 x,y,z,zp,zpsoil,uprt,vprt,w(1,1,1,tim),ptprt(1,1,1,tim),&
                 pprt(1,1,1,tim),qvprt,temscalar(:,:,:,:),              &
                 tke(1,1,1,tim),kmh,kmv,                                &
                 ubar,vbar,wbar,ptbar,pbar,rhostr,qvbar,                &
                 soiltyp,stypfrct,vegtyp,lai,roufns,veg,                &
                 tsoil,qsoil,wetcanp,snowdpth,                          &
                 raing,rainc,prcrate,                                   &
                 radfrc,radsw,rnflx,radswnet,radlwin,                   &
                 usflx,vsflx,ptsflx,qvsflx,                             &
                 ireturn, tem6,tem7,tem8)

    END IF
    IF (mp_opt > 0) CALL mpbarrier
  END DO

  IF( ireturn /= 0) THEN

    WRITE(6,'(3(/1x,a))')                                               &
        'Error occured when reading history data ',                     &
        inifile(1:lenfil)//' or '//inigbf(1:lengbf),                    &
        'Program stopped in EXTINIT.'
    CALL arpsstop ("arpsstop called from extinit in reading file",1)

  END IF

! Added 09/08/2006 by DTD to handle case where MY multimoment scheme is used
! and the other two moments are not initialized.
! The mixing ratio is used to diagnose initial values of these moments
! based on the single-moment version of the scheme

  IF(mphyopt >= 9 .and. nscalarin <= 6) THEN
    WRITE(6,*) 'Diagnosing additional moments for multi-moment scheme'
    CALL init_MM(nx,ny,nz,nts,tim,ptprt,ptbar,pprt,pbar,temscalar)
  END IF

  DO nq = 1, nscalar
    DO k = 1, nz
      DO j = 1, ny
        DO i = 1, nx
          qscalar(i,j,k,tim,nq) = temscalar(i,j,k,nq)
        END DO
      END DO
    END DO
  END DO
  DEALLOCATE ( temscalar )

!-----------------------------------------------------------------------
!
!  To restore the original runname and comments, etc.
!  from ARPS input file.
!
!-----------------------------------------------------------------------

  eps = 0.0001

  IF(mapproj /= mapprojsv .OR.ABS(trulat1sv-trulat1) > eps              &
        .OR.ABS(trulat2sv-trulat2) > eps                                &
        .OR.ABS(trulonsv -trulon ) > eps) THEN
    WRITE(6,'(/2(/5x,a),2(/5x,3(a,f13.3))/)')                           &
        'WARNING: Map projection parameters in the input data do not ', &
        'match those specified in input file. Values in data used.',    &
        'In data,  trulat1=',trulat1  ,', trulat2=',trulat2,            &
        ', trulon=',trulon,                                             &
        'In input, trulat1=',trulat1sv,', trulat2=',trulat2sv,          &
        ', trulon=',trulonsv
  END IF

  IF(ABS(ctrlatsv-ctrlat) > eps .OR. ABS(ctrlonsv-ctrlon) > eps ) THEN
    WRITE(6,'(/3(/5x,a),2(/5x,2(a,f13.3))/)')                           &
        'WARNING: Central latitude and/or longitude of the grid ',      &
        'in the input data do not match those specified in input file.', &
        'Values in data used.',                                         &
        'In data , ctrlat=',ctrlat  ,', ctrlon=',ctrlon,                &
        'In input, ctrlat=',ctrlatsv,', ctrlon=',ctrlonsv
  END IF

  CALL ctim2abss(year,month,day,hour,minute,second, abstsec0)
  CALL ctim2abss(yearsv,monthsv,daysv,hoursv,minutesv,secondsv,         &
                 abstsec1)

  abstsec0 = abstsec0 + nint(time_tmp)
  abstsec1 = abstsec1 + nint(tstart)

  IF ( abstsec0 /= abstsec1 ) THEN
    WRITE(6,'(a,2(/a,1x,i4.4,5(a,i2.2),a,f10.3))')                      &
        'WARNING: Data time is inconsistent with user input time.',     &
        '         Data initime:', year,'-',month,'-',day,'.',          &
                                 hour,':',minute,':',second,            &
                                 '     model data time: ', tstart,      &
        '         User initime:', yearsv,'-',monthsv,'-',daysv,'.',  &
                                 hoursv,':',minutesv,':',secondsv,      &
                                 ' model starting time: ', time_tmp

    IF ( timeopt == 1 ) THEN
      WRITE(6,'(a)') 'Program continues using user specified time'
      year = yearsv
      month = monthsv
      day = daysv
      hour = hoursv
      minute = minutesv
      second = secondsv
    ELSE IF ( timeopt == 2 ) THEN
      WRITE(6,'(a)') 'Program continues using data time'
      tstart = time_tmp
    ELSE
      WRITE(6,'(a)') 'Program stopped in subroutine EXTINIT'
      CALL arpsstop ("arpsstop called from extinit due to timeopt",1)
    END IF
  ELSE
    year = yearsv
    month = monthsv
    day = daysv
    hour = hoursv
    minute = minutesv
    second = secondsv
    IF (myproc == 0) &
    WRITE(6,'(1x,a,i4.4,5(a,i2.2),a,f10.3,a)')                          &
         'Use specified initial time: ',                                &
          year,'-',month,'-',day,'.',hour,':',minute,':',second,        &
         ' and model starting time: ', tstart, ' seconds'
  END IF

  runname = runnamesv(1:80)
  tstop = tstopsv
  thisdmp = thisdmpsv
  dx = dxsv
  dy = dysv
  strhopt = strhoptsv
  dz = dzsv
  dzmin = dzminsv
  zrefsfc = zrefsfcsv
  dlayer1 = dlayer1sv
  dlayer2 = dlayer2sv
  strhtune = strhtunesv
  zflat = zflatsv

  nocmnt = nocmnt_old
  DO i=1,nocmnt
    cmnt(i)=cmnt_old(i)
  END DO

  IF ( ABS(ntcloudsv-ntcloud) > eps .OR. ABS(n0rainsv-n0rain) > eps .OR. &
       ABS(n0snowsv-n0snow) > eps .OR.  ABS(n0grplsv-n0grpl) > eps .OR. &
       ABS(n0hailsv-n0hail) > eps ) THEN
    WRITE(6,'(/,1x,a,/,1x,a,/,10x,a,2(/10x,5(a,f13.3))/)')                     &
        'WARNING: Cloud Number concentration and/or intercept parameters', &
        'for rainwater, snow, graupel, hail ',     &
        'in the input data do not match those specified in input file.',&
        'In data , ntcloud=',ntcloud,', n0rain=',n0rain,  ', n0snow=',n0snow, &
        ', n0grpl=',n0grpl,', n0hail=',n0hail, &
        'In input, ntcloud=',ntcloudsv,', n0rain=',n0rainsv,', n0snow=',n0snowsv, &
        ', n0grpl=',n0grplsv,', n0hail=',n0hailsv
    IF (dsdpref == 1) THEN
      WRITE(6,'(10x,a,/)') 'Values from data file will be used.'
    ELSE
      WRITE(6,'(10x,a,/)') 'Values from namelist will be used.'
      ntcloud = ntcloudsv
      n0rain = n0rainsv
      n0snow = n0snowsv
      n0grpl = n0grplsv
      n0hail = n0hailsv
    END IF
  END IF

  IF ( ABS(rhoicesv-rhoice) > eps .OR. ABS(rhosnowsv-rhosnow) > eps .OR. &
       ABS(rhogrplsv-rhogrpl) > eps .OR. ABS(rhohailsv-rhohail) > eps) THEN
    WRITE(6,'(/,1x,a,/,10x,a,2(/10x,4(a,f13.3))/)')                     &
        'WARNING: Ice, snow graupel and/or hail density ',              &
        'in the input data do not match those specified in input file.',&
        'In data , rhoice=',rhoice,', rhosnow=',rhosnow,      &
        ', rhogrpl=',rhogrpl,', rhohail=',rhohail,            &
        'In input, rhoice=',rhoicesv,', rhosnow=',rhosnowsv,  &
        ', rhogrpl=',rhogrplsv,', rhohail=',rhohailsv
    IF (dsdpref == 1) THEN
      WRITE(6,'(10x,a,/)') 'Values from data file will be used.'
    ELSE
      WRITE(6,'(10x,a,/)') 'Values from namelist will be used.'
      rhoice = rhoicesv
      rhosnow = rhosnowsv
      rhogrpl = rhogrplsv
      rhohail = rhohailsv
    END IF
  END IF

  IF ( ABS(alpharainsv-alpharain) > eps .OR. ABS(alphaicesv-alphaice) > eps .OR. &
       ABS(alphasnowsv-alphasnow) > eps .OR. ABS(alphagrplsv-alphagrpl) > eps .OR. &
       ABS(alphahailsv-alphahail) > eps) THEN
    WRITE(6,'(/,1x,a,/,10x,a,2(/10x,4(a,f13.3))/)')                      &
        'WARNING: rain,ice,snow,graupel, and/or hail shape parameters ', &
        'in the input data do not match those specified in input file.', &
        'In data, alpharain=',alpharain,', alphaice=',alphaice,          &
        ', alphasnow=',alphasnow,', alphagrpl=',alphagrpl,               &
        ', alphahail=',alphahail
     IF (dsdpref == 1) THEN
       WRITE(6,'(10x,a,/)') 'Values from data file will be used.'
     ELSE
       WRITE(6,'(10x,a,/)') 'Values from namelist will be used.'
       alpharain=alpharainsv
       alphaice=alphaicesv
       alphasnow=alphasnowsv
       alphagrpl=alphagrplsv
       alphahail=alphahailsv
     END IF
  END IF
!
!
!-----------------------------------------------------------------------
!
!  Calculate rhostr & jacobians, set hterain.
!
!-----------------------------------------------------------------------
!

  CALL jacob(nx,ny,nz,x,y,z,zp,j1,j2,j3,tem6)

  DO k= 1,nz-1
    DO j= 1,ny-1
      DO i= 1,nx-1
        tvbar=(ptbar(i,j,k)*((pbar(i,j,k)*p0inv)**rddcp))*              &
            ((1.0+rvdrd*qvbar(i,j,k))/(1.0+qvbar(i,j,k)))
        rhostr(i,j,k)= ABS(j3(i,j,k))*pbar(i,j,k)/(rd*tvbar)

!The following is used to define rhostr when initializing the base state
!from a sounding while the above defines rhobar as the virtual
!density. rhostr is there somewhat different for initopt=3 and 4 which calls
!extinit.
!For the code to be consistent with the initialization of
!rhostr/rhobar in INIBASE, use the following version.
!
!       rhostr(i,j,k)= abs(j3(i,j,k))*pbar(i,j,k)/                      &
!           ( Rd * ptbar(i,j,k)*(pbar(i,j,k)*p0inv)**rddcp )

      END DO
    END DO
  END DO

  DO j=1,ny-1
    DO i=1,nx-1
      rhostr(i,j,   1)=rhostr(i,j,2)
      rhostr(i,j,nz-1)=rhostr(i,j,nz-2)
    END DO
  END DO

  DO i=1,nx
    DO j=1,ny
      hterain(i,j) = zp(i,j,2)
    END DO
  END DO

!
!  call the soil model grid initialization
!
! call jacobsoil...stopped here...add the code.

 CALL inisoilgrd(nx,ny,nzsoil,hterain,zpsoil,j3soil,j3soilinv)

!
!-----------------------------------------------------------------------
!
!  Put the 3-d variables into their 4-d arrays.
!
!-----------------------------------------------------------------------
!
  tim = tpresent

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx
        u(i,j,k,tim)=uprt(i,j,k)+ubar(i,j,k)
      END DO
    END DO
  END DO


  DO k=1,nz-1
    DO j=1,ny
      DO i=1,nx-1
        v(i,j,k,tim)=vprt(i,j,k)+vbar(i,j,k)
      END DO
    END DO
  END DO

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        qv(i,j,k,tim)=qvprt(i,j,k)+qvbar(i,j,k)
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!TINA initialize cc bubble here if using external data initialization
! but ccin = 1
!This initializes tpast,present,future
!-----------------------------------------------------------------------
  IF (ccin == 1 .AND. cpoint < 0 .AND. P_CC > 0) then
    CALL ccinit(nx,ny,nz,3,x,y,z,zp,qscalar(:,:,:,:,P_CC),tem8)
    IF (myproc == 0) print *,'Initializing cc as bubble 2'
    IF (myproc == 0) print *,'max(cc) =', maxval(qscalar(:,:,:,:,P_CC))
  END IF
!TINA
!
!-----------------------------------------------------------------------
!
!  Set the future values of variables to their current values.
!  This is done primarily for safety reasons.
!
!-----------------------------------------------------------------------
!
  IF( initopt < 0 .or. initopt > 4 ) then
     write(6,'(a,i10)') 'Value of initopt incorrect. It was ', initopt
     CALL arpsstop ("arpsstop called from EXTINIT ",1)
  ENDIF

  IF(initopt == 3) THEN

    IF(nts > 1 ) THEN
      DO k=1,nz
        DO j=1,ny
          DO i=1,nx
            u    (i,j,k,tfuture) = u    (i,j,k,tpresent)
            v    (i,j,k,tfuture) = v    (i,j,k,tpresent)
            w    (i,j,k,tfuture) = w    (i,j,k,tpresent)
            ptprt(i,j,k,tfuture) = ptprt(i,j,k,tpresent)
            pprt (i,j,k,tfuture) = pprt (i,j,k,tpresent)
            qv   (i,j,k,tfuture) = qv   (i,j,k,tpresent)
            tke  (i,j,k,tfuture) = tke  (i,j,k,tpresent)

            u    (i,j,k,tpast  ) = u    (i,j,k,tpresent)
            v    (i,j,k,tpast  ) = v    (i,j,k,tpresent)
            w    (i,j,k,tpast  ) = w    (i,j,k,tpresent)
            ptprt(i,j,k,tpast  ) = ptprt(i,j,k,tpresent)
            pprt (i,j,k,tpast  ) = pprt (i,j,k,tpresent)
            qv   (i,j,k,tpast  ) = qv   (i,j,k,tpresent)
            tke  (i,j,k,tpast  ) = tke  (i,j,k,tpresent)
          END DO
        END DO
      END DO

      DO nq = 1,nscalar
        DO k = 1,nz
          DO j = 1,ny
            DO i = 1,nx
              qscalar(i,j,k,tfuture,nq) = qscalar(i,j,k,tpresent,nq)
              qscalar(i,j,k,tpast,  nq) = qscalar(i,j,k,tpresent,nq)
            END DO
          END DO
        END DO
      END DO
    END IF
!
    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          ptcumsrc(i,j,k)=0.0
          qcumsrc(i,j,k,1)=0.0
          qcumsrc(i,j,k,2)=0.0
          qcumsrc(i,j,k,3)=0.0
          qcumsrc(i,j,k,4)=0.0
          qcumsrc(i,j,k,5)=0.0
        END DO
      END DO
    END DO

  ENDIF

!-----------------------------------------------------------------------
!
!  Print out the max/min of initial arrays read in.
!
!-----------------------------------------------------------------------
!
  tim = tpresent

  IF (myproc ==0) &
    WRITE(6,'(1x,a/)') 'Min. and max. of the data arrays read in:'

  CALL a3dmax0(x,1,nx,1,nx,1,1,1,1, 1,1,1,1, amax,amin)
  IF (myproc ==0) &
    WRITE(6,'(1x,2(a,e13.6))') 'xmin    = ', amin,',  xmax    =',amax

  CALL a3dmax0(y,1,ny,1,ny,1,1,1,1, 1,1,1,1, amax,amin)
  IF (myproc ==0) &
    WRITE(6,'(1x,2(a,e13.6))') 'ymin    = ', amin,',  ymax    =',amax

  CALL a3dmax0(z,1,nz,1,nz,1,1,1,1, 1,1,1,1, amax,amin)
  IF (myproc ==0) &
    WRITE(6,'(1x,2(a,e13.6))') 'zmin    = ', amin,',  zmax    =',amax

  CALL a3dmax0(zp,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz,                    &
               amax,amin)
  IF (myproc ==0) &
    WRITE(6,'(1x,2(a,e13.6))') 'zpmin   = ', amin,',  zpmax   =',amax

  CALL a3dmax0(j3,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,                  &
               amax,amin)
  IF (myproc ==0) &
    WRITE(6,'(1x,2(a,e13.6))') 'j3min   = ', amin,',  j3max   =',amax

  CALL a3dmax0(ubar,1,nx,1,nx,1,ny,1,ny-1,1,nz,1,nz-1,                  &
               amax,amin)
  IF (myproc ==0) &
    WRITE(6,'(1x,2(a,e13.6))') 'ubarmin = ', amin,',  ubarmax =',amax

  CALL a3dmax0(vbar,1,nx,1,nx-1,1,ny,1,ny,1,nz,1,nz-1,                  &
               amax,amin)
  IF (myproc ==0) &
    WRITE(6,'(1x,2(a,e13.6))') 'vbarmin = ', amin,',  vbarmax =',amax

  CALL a3dmax0(ptbar,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,               &
               amax,amin)
  IF (myproc ==0) &
    WRITE(6,'(1x,2(a,e13.6))') 'ptbarmin= ', amin,',  ptbarmax=',amax

  CALL a3dmax0(pbar,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,                &
               amax,amin)
  IF (myproc ==0) &
    WRITE(6,'(1x,2(a,e13.6))') 'pbarmin = ', amin,',  pbarmax =',amax

  CALL a3dmax0(rhostr,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,              &
               amax,amin)
  IF (myproc ==0) &
    WRITE(6,'(1x,2(a,e13.6))') 'rhostrmin=', amin,', rhostrmax=',amax

  CALL a3dmax0(qvbar,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,               &
               amax,amin)
  IF (myproc ==0) &
    WRITE(6,'(1x,2(a,e13.6))') 'qvbarmin= ', amin,',  qvbarmax=',amax

  CALL a3dmax0(u(1,1,1,tim),1,nx,1,nx,1,ny,1,ny-1,1,nz,1,nz-1,          &
               amax,amin)
  IF (myproc ==0) &
    WRITE(6,'(1x,2(a,e13.6))') 'umin    = ', amin,',  umax    =',amax

  CALL a3dmax0(v(1,1,1,tim),1,nx,1,nx-1,1,ny,1,ny,1,nz,1,nz-1,          &
               amax,amin)
  IF (myproc ==0) &
    WRITE(6,'(1x,2(a,e13.6))') 'vmin    = ', amin,',  vmax    =',amax

  CALL a3dmax0(w(1,1,1,tim),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz,          &
               amax,amin)
  IF (myproc ==0) &
    WRITE(6,'(1x,2(a,e13.6))') 'wmin    = ', amin,',  wmax    =',amax

  CALL a3dmax0(ptprt(1,1,1,tim),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,    &
               amax,amin)
  IF (myproc ==0) &
    WRITE(6,'(1x,2(a,e13.6))') 'ptprtmin= ', amin,',  ptprtmax=',amax

  CALL a3dmax0(pprt(1,1,1,tim),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,     &
               amax,amin)
  IF (myproc ==0) &
    WRITE(6,'(1x,2(a,e13.6))') 'pprtmin = ', amin,',  pprtmax =',amax

  CALL a3dmax0(qv(1,1,1,tim),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,       &
               amax,amin)
  IF (myproc ==0) &
    WRITE(6,'(1x,2(a,e13.6))') 'qvmin   = ', amin,',  qvmax   =',amax

  DO nq = 1,nscalar
    CALL a3dmax0(qscalar(1,1,1,tim,nq),1,nx,1,nx-1,1,ny,1,ny-1,         &
                 1,nz,1,nz-1,amax,amin)
    IF (myproc ==0) WRITE(6,'(1x,2(a,e13.6))')     &
                                 TRIM(qnames(nq))//'min   = ', amin,    &
                          ',  '//TRIM(qnames(nq))//'max   =',amax
  END DO

  CALL a3dmax0(tke(1,1,1,tim),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,      &
               amax,amin)
  IF (myproc ==0) &
    WRITE(6,'(1x,2(a,e13.6))') 'tkemin  = ', amin,',  tkemax  =',amax

!  IF ( sfcphy.gt.0 ) THEN

!    CALL a3dmax0(tsoil,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1, amax,amin)
!    write(6,'(1x,2(a,e13.6))') 'tsoilmin= ',amin,', tsoilmax =',amax

!    CALL a3dmax0(qsoil,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1, amax,amin)
!    write(6,'(1x,2(a,e13.6))') 'qsoilmin = ',amin,', qsoilmax  =',amax

!    CALL a3dmax0(wetcanp,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1, amax,amin)
!    write(6,'(1x,2(a,e13.6))') 'wetcmin = ',amin,', wetcmax  =',amax
!  ENDIF

  CALL a3dmax0(ptcumsrc,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,            &
               amax,amin)
  IF (myproc ==0) &
    WRITE(6,'(1x,2(a,e13.6))') 'ptcummin= ', amin,',  ptcummax=',amax

  CALL a3dmax0(qcumsrc(1,1,1,1),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,    &
               amax,amin)
  IF (myproc ==0) &
    WRITE(6,'(1x,2(a,e13.6))') 'qvcummin= ', amin,',  qvcummax=',amax
  CALL a3dmax0(qcumsrc(1,1,1,2),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,    &
               amax,amin)
  IF (myproc ==0) &
    WRITE(6,'(1x,2(a,e13.6))') 'qccummin= ', amin,',  qccummax=',amax
  CALL a3dmax0(qcumsrc(1,1,1,3),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,    &
               amax,amin)
  IF (myproc ==0) &
    WRITE(6,'(1x,2(a,e13.6))') 'qrcummin= ', amin,',  qrcummax=',amax
  CALL a3dmax0(qcumsrc(1,1,1,4),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,    &
               amax,amin)
  IF (myproc ==0) &
    WRITE(6,'(1x,2(a,e13.6))') 'qicummin= ', amin,',  qicummax=',amax
  CALL a3dmax0(qcumsrc(1,1,1,5),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,    &
               amax,amin)
  IF (myproc ==0) &
    WRITE(6,'(1x,2(a,e13.6))') 'qscummin= ', amin,',  qscummax=',amax

  RETURN

END SUBROUTINE extinit
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE INITSFC                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE initsfc(nx,ny,nz,nzsoil,nstyps,                              &
           zpsoil,                                                      &
           pbar,pprt,ptbar,ptprt,qvbar,qv,                              &
           soiltyp,stypfrct,vegtyp,lai,roufns,veg,                      &
           ndvi,                                                        &
           tsoil,qsoil,wetcanp,snowdpth,qvsfc,tem1soil)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Initialize the surface characteristics data and soil model
!  variables according to option parameters sfcdat and soilinit.
!
!  The surface and soil data files are sequential binary files.
!
!  Their structures should be:
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  02/17/94
!
!  MODIFICATION:
!
!  02/07/1995 (Yuhe Liu)
!  Added a new 2-D permanent array, veg(nx,ny), to the soil model and
!  at the same time delete the table data array veg(14).
!
!  01/29/1995 (V. Wong and X. Song)
!  Add a flag wtrexist in "initsfc".
!
!  12/04/1997 (Yuhe Liu)
!  Set the soil variables to their default value read in from input
!  file if they are not included in soilinit file.
!
!  12/09/1998 (Donghai Wang)
!  Added the snow cover.
!
!  2000/01/07 (Gene Bassett)
!  Modified the behavior for the case where soil data read in from file
!  and also present in history file.  If soilinit=2, initopt=3 and
!  sfcin=1, for soil variables that are not in soildata file the values
!  in the history file are used.
!
!  05/17/2002 (Dan Weber)
!  Added new soil model variables.
!
!  15 June 2002 (Eric Kemp)
!  Bug fixes.
!
!-----------------------------------------------------------------------
!
!  INPUT and OUTPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the z-direction (sfc/top)
!    nzsoil   Number of grid points in the soil model in the -z-direction
!
!    pbar     Base state pressure
!    pprt     Preturbation pressure
!
!    ptbar    Base state potential temperature
!    ptprt    Preturbation potential temperature
!
!    qvbar    Base state water vapor mixing ratio
!    qv       Water vapor mixing ratio
!
!    soiltyp  Soil type at the horizontal grid points
!    stypfrct  Soil type fraction
!    vegtyp   Vegetation type at the horizontal grid points
!    lai      Leaf Area Index
!    roufns   Surface roughness
!    veg      Vegetation fraction
!
!    tsoil    Soil temperature (K)
!    qsoil    Soil moisture (m**3/m**3)
!    wetcanp  Canopy water amount
!    qvsfc    Effective S.H. at sfc.
!
!  TEMPORARY WORKING ARRAY
!
!    tem1soil Soil model temporary work array.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations. (Local Variables)
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  REAL :: pi
  PARAMETER ( pi = 3.141592654)

  INTEGER :: nx, ny, nz
  INTEGER :: nzsoil            ! Number of grid points in the -z-direction

  REAL :: zpsoil (nx,ny,nzsoil) ! soil model points elevations. (m)
  REAL :: pbar   (nx,ny,nz) ! Base state pressure (Pascal)
  REAL :: pprt   (nx,ny,nz) ! Perturbation pressure (Pascal)
  REAL :: ptbar  (nx,ny,nx) ! Base state potential temperature (K)
  REAL :: ptprt  (nx,ny,nz) ! Perturbation potential temperature (K)
  REAL :: qvbar  (nx,ny,nz) ! Base state specific humidity (kg/kg)
  REAL :: qv     (nx,ny,nz) ! Total specific humidity (kg/kg)

  INTEGER :: nstyps                  ! Number of soil types
  INTEGER :: soiltyp(nx,ny,nstyps)   ! Soil type at each point
  REAL :: stypfrct(nx,ny,nstyps)  ! Fraction of soil types
  INTEGER :: vegtyp (nx,ny)          ! Vegetation type at each point
  REAL :: lai    (nx,ny)          ! Leaf Area Index
  REAL :: roufns (nx,ny)          ! Surface roughness
  REAL :: veg    (nx,ny)          ! Vegetation fraction
  REAL :: ndvi   (nx,ny)          ! NDVI

  REAL :: tsoil (nx,ny,nzsoil,0:nstyps) ! Soil layer temperature (K)
  REAL :: qsoil (nx,ny,nzsoil,0:nstyps) ! Soil layer moisture(m**3/m**3)

  REAL :: wetcanp(nx,ny,0:nstyps) ! Canopy moisture
  REAL :: snowdpth(nx,ny)         ! Snow depth (m)
  REAL :: qvsfc  (nx,ny,0:nstyps) ! Effective surface specific
                                  ! humidity (kg/kg)

!
!-----------------------------------------------------------------------
!
!  Define local variables
!
!-----------------------------------------------------------------------
!
  REAL :: psfc        ! Surface pressure
  REAL :: rhgs        ! The relative humidity at ground surface
  REAL :: qvsatts     ! Saturated specific humidity at surface

  REAL :: wrmax       ! Maximum value for canopy moisture, wetcanp

  REAL :: pterm       ! Temporary variable, real positive term flag

  REAL :: p0inv       ! Inverse of p0 (1000 mb)

  REAL :: tema        ! Temporary variable
  REAL :: temb        ! Temporary variable
!
!-----------------------------------------------------------------------
!
!  Global constants and parameters, most of them specify the
!  model run options.
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'indtflg.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'soilcst.inc'
  INCLUDE 'bndry.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,is,imid,jmid,ii
  INTEGER :: zpsoilin
  INTEGER :: tsoilin
  INTEGER :: qsoilin
  INTEGER :: wcanpin
  INTEGER :: snowdin

  INTEGER :: astat
  INTEGER, ALLOCATABLE :: mpitmp(:)
  REAL,    ALLOCATABLE :: mprtmp(:)

  REAL :: tem1soil(nx,ny,nzsoil)    ! Temporary soil model work array.
!
!-----------------------------------------------------------------------
!
!  Function f_qvsat and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
  REAL :: f_qvsat

!fpp$ expand (f_qvsat)
!!dir$ inline always f_qvsat
!*$*  inline routine (f_qvsat)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (mp_opt > 0 .AND. (sfcdat == 1 .OR. soilinit == 1)) THEN
    ALLOCATE(mpitmp(MAX(nx,ny)*2),stat=astat)
    CALL check_alloc_status(astat, "initsfc:mpitmp")

    ALLOCATE(mprtmp(MAX(nx,ny)*2),stat=astat)
    CALL check_alloc_status(astat, "initsfc:mprtmp")
  END IF

  p0inv = 1.0/p0

  IF ( initopt == 2 ) THEN
    IF(myproc ==0)   WRITE (6, '(a,a/a/a/a/a)')                         &
        'Model initialization option was set to restart. ',             &
        'All the surface and soil arrays should have been read',        &
        'in from the restart file. No more initialization is ',         &
        'done in subroutine INITSFC.'

    wtrexist=0
    DO j = 1, ny
      DO i = 1, nx
        DO is = 1,nstyp
          IF (soiltyp(i,j,is) == 13) THEN
            wtrexist=1
            RETURN
          END IF
        END DO
      END DO
    END DO

    RETURN

  END IF

  IF ( sfcdat == 1 ) THEN

    nstyp = 1

    DO j=1, ny-1
      DO i=1, nx-1
        soiltyp(i,j,1) = styp
        vegtyp (i,j) = vtyp
        lai    (i,j) = lai0
        roufns (i,j) = roufns0
        veg    (i,j) = veg0
      END DO
    END DO

  ELSE IF (sfcdat == 2 .OR. (sfcdat == 3.AND.landin /= 1) ) THEN
!
!-----------------------------------------------------------------------
!
!  Read the surface property data from file sfcdtfl (passed
!  in globcst.inc).
!
!-----------------------------------------------------------------------
!
!  blocking inserted for ordering i/o for message passing
    DO i=0,nprocs-1,readstride
      IF(myproc >= i.AND.myproc <= i+readstride-1)THEN

        IF (mp_opt > 0 .AND. readsplit(FINDX_T) > 0) THEN

        CALL readsplitsfcdt( nx,ny,nstyps,sfcdtfl,dx,dy,                &
               mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,     &
                      soiltyp,stypfrct,vegtyp,lai,roufns,veg,ndvi )
        ELSE

        CALL readsfcdt( nx,ny,nstyps,sfcdtfl,dx,dy,                     &
               mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,     &
                      soiltyp,stypfrct,vegtyp,lai,roufns,veg,ndvi )

        END IF

      END IF
      IF (mp_opt > 0) CALL mpbarrier
    END DO

  ELSE IF(sfcdat == 3 .AND. landin == 1) THEN

    IF (myproc == 0)                                                    &
      WRITE(6,'(1x,a/,a/)')                                             &
        'Surface property data in the initialization file was used.',   &
        'No more initialization performed.'

  ELSE

    WRITE(6,'(1x,a,i3,a/)')                                             &
        'Invalid surface data input option. sfcdat =',sfcdat,           &
        '. Program stopped in INITSFC.'
      CALL arpsstop ("arpsstop called from initsfc with sfcdat option",1)

  END IF

!  print *, ' nstyps = ',nstyps
!  imid=(nx/2)+1
!  jmid=(ny/2)+1
!  DO is=1,nstyps
!    print *, '  Sample soil ( ',imid,jmid,') = ',                    &
!             soiltyp(imid,jmid,is),stypfrct(imid,jmid,is)
!  END DO

  IF ( nstyp == 1 ) THEN
    DO j=1,ny
      DO i=1,nx
        stypfrct(i,j,1) = 1.0
      END DO
    END DO
  END IF

  IF ( soilinit == 1 ) THEN
!
!-----------------------------------------------------------------------
!
!  All surface variables are specified uniformly in horizontal from
!  the input namelist.
!
!-----------------------------------------------------------------------
    DO k=1,nzsoil
      DO j=1, ny-1
        DO i=1, nx-1
          psfc = .5 * ( (pbar(i,j,1)+pprt(i,j,1))                         &
                      + (pbar(i,j,2)+pprt(i,j,2)) )

          tsoil (i,j,k,0) = tsoilint(k)
          qsoil (i,j,k,0) = qsoilint(k)

          IF( soiltyp(i,j,1) == 13) THEN  ! For water
            tsoil(i,j,1,0) = ptswtr0*(psfc*p0inv)**rddcp
          ELSE                          ! For land
            tsoil(i,j,1,0) = ptslnd0*(psfc*p0inv)**rddcp
          END IF

          wetcanp(i,j,0) = wetcanp0
          snowdpth(i,j)  = snowdpth0
        END DO
      END DO
    END DO

    DO is=1,nstyp
      DO k=1,nzsoil
        DO j=1,ny-1
          DO i=1,nx-1
            tsoil  (i,j,k,is) = tsoil (i,j,k,0)
            qsoil  (i,j,k,is) = qsoil (i,j,k,0)
          END DO
        END DO
      END DO
      DO j=1,ny-1
        DO i=1,nx-1
          wetcanp(i,j,is) = wetcanp(i,j,0)
        END DO
      END DO
    END DO

  ELSE IF ( (soilinit == 2) .OR. (soilinit == 5)                        &
            .OR. (soilinit == 3 .AND. sfcin /= 1) ) THEN
!
!-----------------------------------------------------------------------
!
!  Read the soil variables from file soilinfl. soilinfl is
!  passed in globcst.inc.
!
!-----------------------------------------------------------------------
!
!  blocking inserted for ordering i/o for message passing
    DO ii=0,nprocs-1,readstride
      IF(myproc >= ii.AND.myproc <= ii+readstride-1) THEN

        IF (mp_opt >0 .AND. readsplit(FINDX_S) > 0) THEN

        CALL readsplitsoil(nx,ny,nzsoil,nstyps,soilinfl,dx,dy,zpsoil,   &
               mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,     &
               zpsoilin,tsoilin,qsoilin,wcanpin,snowdin,                &
               tsoil,qsoil,wetcanp,snowdpth,soiltyp)
        ELSE

        CALL readsoil(nx,ny,nzsoil,nstyps,soilinfl,dx,dy,zpsoil,        &
               mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,     &
               zpsoilin,tsoilin,qsoilin,wcanpin,snowdin,                &
               tsoil,qsoil,wetcanp,snowdpth,soiltyp)
        END IF

      END IF
      IF(soilinit == 5)THEN
        DO is=1,nstyp
          DO k=1,nzsoil
            DO j=1,ny-1
              DO i=1,nx-1
                tsoil  (i,j,k,is) = tsoil (i,j,k,0)
                qsoil  (i,j,k,is) = qsoil (i,j,k,0)
              END DO
            END DO
          END DO
        END DO
      END IF

      IF (mp_opt > 0) CALL mpbarrier
    END DO

    IF (sfcin /= 1) THEN
      IF ( tsoilin /= 1 ) THEN
        DO is=1,nstyp
          DO j=1, ny-1
            DO i=1, nx-1
              psfc = .5 * ( (pbar(i,j,1)+pprt(i,j,1))                   &
                          + (pbar(i,j,2)+pprt(i,j,2)) )

              IF( soiltyp(i,j,1) == 13) THEN  ! For water
                tsoil(i,j,1,is) = ptswtr0*(psfc*p0inv)**rddcp
              ELSE                          ! For land
                tsoil(i,j,1,is) = ptslnd0*(psfc*p0inv)**rddcp
              END IF
            END DO
          END DO
        END DO

        DO is=1,nstyp
          DO k=2, nzsoil
            DO j=1, ny-1
              DO i=1, nx-1
                tsoil(i,j,k,is) = tsoilint(k)
              END DO
            END DO
          END DO
        END DO
      END IF

      IF ( qsoilin /= 1 ) THEN
        DO is=1,nstyp
          DO j=1, ny-1
            DO i=1, nx-1
              qsoil(i,j,1,is) = qsoilint(1)
            END DO
          END DO
        END DO
        DO is=1,nstyp
          DO k=2, nzsoil
            DO j=1, ny-1
              DO i=1, nx-1
                qsoil(i,j,k,is) = qsoilint(2)
              END DO
            END DO
          END DO
        END DO
      END IF

      IF ( wcanpin /= 1 ) THEN
        DO is=1,nstyp
          DO j=1, ny-1
            DO i=1, nx-1
              wetcanp(i,j,is) = wetcanp0
            END DO
          END DO
        END DO
      END IF

      IF ( snowdin /= 1 ) THEN
        DO j=1, ny-1
          DO i=1, nx-1
            snowdpth(i,j) = snowdpth0
          END DO
        END DO
      END IF

    END IF

  ELSE IF(soilinit == 3 .AND. sfcin == 1) THEN

    IF (myproc == 0)                                                    &
      WRITE(6,'(1x,a/,a/)')                                             &
        'Surface variable in the initialization file was used.',        &
        'No more initialization performed.'

    IF ( nstyp == 1 ) THEN
      DO k=1,nzsoil
        DO j=1,ny-1
          DO i=1,nx-1
            tsoil  (i,j,k,1) = tsoil  (i,j,k,0)
            qsoil  (i,j,k,1) = qsoil  (i,j,k,0)
          END DO
        END DO
      END DO
      DO j=1,ny-1
        DO i=1,nx-1
          wetcanp(i,j,1) = wetcanp(i,j,0)
        END DO
      END DO
    END IF
  ELSE IF(soilinit == 4) THEN
    p0inv=1./p0
    DO j=1,ny-1
      DO i=1,nx-1
        tema=0.5*((ptbar(i,j,2)+ptprt(i,j,2))                           &
                 +(ptbar(i,j,1)+ptprt(i,j,1)))
        psfc=0.5*((pbar(i,j,2)+pprt(i,j,2))                             &
                 +(pbar(i,j,1)+pprt(i,j,1)))

        temb = tema*(psfc*p0inv)**rddcp

        tsoil  (i,j,1,0) = temb + ttprt
        tsoil  (i,j,nzsoil,0) = temb + tbprt

        qsoil  (i,j,1,0) = 0.0
        qsoil  (i,j,2,0) = 0.0
        wetcanp(i,j,0) = wetcanp0
        snowdpth(i,j)  = snowdpth0
      END DO
    END DO

    DO is=1,nstyp
      DO j=1,ny-1
        DO i=1,nx-1
          tsoil  (i,j,1,is) = tsoil (i,j,1,0)
          IF (soiltyp(i,j,is) == 13) THEN
                             !Andreas 30-10-03: Make moisture of water 1
            qsoil(i,j,1,is) = 1.
          ELSE
            qsoil(i,j,1,is) = wgrat*wsat(soiltyp(i,j,is))
          ENDIF
          qsoil(i,j,1,0) = qsoil(i,j,1,0)+stypfrct(i,j,is)*qsoil(i,j,1,is)
          wetcanp(i,j,is) = wetcanp(i,j,0)
        END DO
      END DO
      DO k=2,nzsoil
        DO j=1,ny-1
          DO i=1,nx-1
            tsoil(i,j,k,is) = tsoil(i,j,1,0) - ((k - 1)/(nzsoil - 1))*  &
                               (tsoil(i,j,1,0)-tsoil(i,j,nzsoil,0))
            IF (soiltyp(i,j,is) == 13) THEN
                                 !Andreas 30-10-03: Make moisture of water 1
               qsoil(i,j,k,is) = 1.
            ELSE
               qsoil(i,j,k,is) = w2rat*wsat(soiltyp(i,j,is))
            ENDIF
            qsoil(i,j,k,0)  = qsoil(i,j,k,0)+stypfrct(i,j,is)*qsoil(i,j,k,is)
          END DO
        END DO
      END DO
    END DO

  ELSE

    WRITE(6,'(1x,a,i3,a/)')                                             &
        'Invalid surface variable input option. soilinit =',soilinit,   &
        '. Program stopped in INITSFC.'
    CALL arpsstop ("arpsstop called from initsfc with soilinit option",1)

  END IF

  IF ( moist /= 0 ) THEN
    DO is=1,nstyp
      DO k = 1, nzsoil
        DO j = 1, ny-1
          DO i = 1, nx-1
            qsoil(i,j,k,is) = MAX( qsoil(i,j,k,is), 0.0 )
            qsoil(i,j,k,is) = MIN( qsoil(i,j,k,is), wsat(soiltyp(i,j,is)) )
          END DO
        END DO
      END DO
      DO j = 1, ny-1
        DO i = 1, nx-1
          wrmax = 0.2*1.0E-3*veg(i,j)*lai(i,j)
          wetcanp(i,j,is) = MAX( wetcanp(i,j,is), 0.0 )
          wetcanp(i,j,is) = MIN( wetcanp(i,j,is), wrmax )
        END DO
      END DO
    END DO
    DO k = 1, nzsoil
      DO j = 1, ny-1
        DO i = 1, nx-1
          qsoil (i,j,k,0) = MAX( qsoil(i,j,k,0), 0.0 )
        END DO
      END DO
    END DO
    DO j = 1, ny-1
      DO i = 1, nx-1
        wetcanp(i,j,0) = MAX( wetcanp(i,j,0), 0.0 )
      END DO
    END DO
  ELSE
    DO is=0,nstyp
      DO k = 1, nzsoil
        DO j = 1, ny-1
          DO i = 1, nx-1
            qsoil(i,j,k,is) = 0.0
          END DO
        END DO
      END DO
      wetcanp(:,:,is) = 0.0
    END DO
  END IF

  DO is=1,nstyp
    DO j = 1, ny-1
      DO i = 1, nx-1

        psfc = .5 * ( (pbar(i,j,1)+pprt(i,j,1))                         &
                    + (pbar(i,j,2)+pprt(i,j,2)) )

        qvsatts = f_qvsat( psfc, tsoil(i,j,1,is) )

        IF ( soiltyp(i,j,is) == 12 .OR. soiltyp(i,j,is) == 13) THEN
                                                     ! Ice and water
          rhgs = 1.0

        ELSE IF ( soiltyp(i,j,is) > 0 ) THEN

          pterm = .5+SIGN(.5,qsoil(i,j,1,is)-1.1*wfc(soiltyp(i,j,is)))

          rhgs = pterm + (1.-pterm)                                     &
                       * 0.25 * ( 1.-COS(qsoil(i,j,1,is)*pi              &
                                    /(1.1*wfc(soiltyp(i,j,is)))) )**2
        ELSE

          rhgs = 0.0

        END IF

        qvsfc(i,j,is) = rhgs * qvsatts

      END DO
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Ensure soil values are consistent with each other.
!
!-----------------------------------------------------------------------
!
  tsoil  (:,:,:,0) =  0.0
  qsoil  (:,:,:,0) =  0.0
  wetcanp(:,:,0) =  0.0

  DO is = 1,nstyp
    DO k=1,nzsoil
      DO j=1,ny-1
        DO i=1,nx-1
          tsoil(i,j,k,0) = tsoil(i,j,k,0)                               &
                         + tsoil(i,j,k,is) * stypfrct(i,j,is)
          qsoil(i,j,k,0) = qsoil(i,j,k,0)                               &
                         + qsoil(i,j,k,is) * stypfrct(i,j,is)
        END DO
      END DO
    END DO
    DO j=1,ny-1
      DO i=1,nx-1
        wetcanp(i,j,0) = wetcanp(i,j,0)                                 &
                       + wetcanp(i,j,is) * stypfrct(i,j,is)
      END DO
    END DO
  END DO

  DO j=1,ny-1
    DO i=1,nx-1
      IF ((vegtyp(i,j) == 14) .OR. (tsoil(i,j,1,0) > 273.16)) THEN
        snowdpth(i,j) = 0.  ! Make sure that snow cover is consistent
      END IF                 ! with vegetation type and soil temperature.
    END DO
  END DO

  DO j = 1, ny-1
    DO i = 1, nx-1
      qvsfc(i,j,0) = 0.0
    END DO
  END DO

  DO is=1,nstyp
    DO j = 1, ny-1
      DO i = 1, nx-1
        qvsfc(i,j,0) = qvsfc(i,j,0) + qvsfc(i,j,is)*stypfrct(i,j,is)
      END DO
    END DO
  END DO

  IF ( sfcdat == 1 ) THEN
    DO is=1,nstyp
      IF (mp_opt > 0) THEN
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv1diew(soiltyp(1,1,is),nx,ny,ebc,wbc,0,mpitmp)
        CALL mpsendrecv1dins(soiltyp(1,1,is),nx,ny,nbc,sbc,0,mpitmp)
        CALL mpsendrecv1dew(stypfrct(1,1,is),nx,ny,ebc,wbc,0,mprtmp)
        CALL mpsendrecv1dns(stypfrct(1,1,is),nx,ny,nbc,sbc,0,mprtmp)
      END IF
      CALL acct_interrupt(bc_acct)
      CALL bcis2d(nx,ny, soiltyp (1,1,is), ebc,wbc,nbc,sbc)
      CALL bcs2d (nx,ny, stypfrct(1,1,is), ebc,wbc,nbc,sbc)
      CALL acct_stop_inter
    END DO

    IF (mp_opt > 0) THEN
      CALL acct_interrupt(mp_acct)
      CALL mpsendrecv1diew(vegtyp,nx,ny,ebc,wbc,0,mpitmp)
      CALL mpsendrecv1dins(vegtyp,nx,ny,nbc,sbc,0,mpitmp)
      CALL mpsendrecv1dew(lai,nx,ny,ebc,wbc,0,mprtmp)
      CALL mpsendrecv1dns(lai,nx,ny,nbc,sbc,0,mprtmp)
      CALL mpsendrecv1dew(roufns,nx,ny,ebc,wbc,0,mprtmp)
      CALL mpsendrecv1dns(roufns,nx,ny,nbc,sbc,0,mprtmp)
      CALL mpsendrecv1dew(veg,nx,ny,ebc,wbc,0,mprtmp)
      CALL mpsendrecv1dns(veg,nx,ny,nbc,sbc,0,mprtmp)
      CALL mpsendrecv1dew(snowdpth,nx,ny,ebc,wbc,0,mprtmp)
      CALL mpsendrecv1dns(snowdpth,nx,ny,nbc,sbc,0,mprtmp)
    END IF
    CALL acct_interrupt(bc_acct)
    CALL bcis2d(nx,ny, vegtyp, ebc,wbc,nbc,sbc)
    CALL bcs2d (nx,ny, lai,    ebc,wbc,nbc,sbc)
    CALL bcs2d (nx,ny, roufns, ebc,wbc,nbc,sbc)
    CALL bcs2d (nx,ny, veg,    ebc,wbc,nbc,sbc)
    CALL bcs2d (nx,ny, snowdpth,ebc,wbc,nbc,sbc)
    CALL acct_stop_inter
  END IF

  IF ( soilinit == 1 ) THEN
    DO is=0,nstyp
      IF (mp_opt > 0) THEN
        CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(tsoil(1,1,1,is),nx,ny,nzsoil,ebc,wbc,0,tem1soil)
        CALL mpsendrecv2dns(tsoil(1,1,1,is),nx,ny,nzsoil,nbc,sbc,0,tem1soil)
        CALL mpsendrecv2dew(qsoil(1,1,1,is),nx,ny,nzsoil,ebc,wbc,0,tem1soil)
        CALL mpsendrecv2dns(qsoil(1,1,1,is),nx,ny,nzsoil,nbc,sbc,0,tem1soil)

        CALL mpsendrecv1dew(wetcanp(1,1,is),nx,ny,ebc,wbc,0,mprtmp)
        CALL mpsendrecv1dns(wetcanp(1,1,is),nx,ny,nbc,sbc,0,mprtmp)
        CALL mpsendrecv1dew(qvsfc(1,1,is),  nx,ny,ebc,wbc,0,mprtmp)
        CALL mpsendrecv1dns(qvsfc(1,1,is),  nx,ny,nbc,sbc,0,mprtmp)
      END IF
      CALL acct_interrupt(bc_acct)

      DO k=1,nzsoil
        CALL bcs2d (nx,ny, tsoil  (1,1,k,is), ebc,wbc,nbc,sbc)
        CALL bcs2d (nx,ny, qsoil  (1,1,k,is), ebc,wbc,nbc,sbc)
      END DO

      CALL bcs2d (nx,ny, wetcanp(1,1,is), ebc,wbc,nbc,sbc)
      CALL bcs2d (nx,ny, qvsfc  (1,1,is), ebc,wbc,nbc,sbc)
      CALL acct_stop_inter
    END DO
  END IF

  IF (mp_opt > 0 .AND. (sfcdat == 1 .OR. soilinit == 1)) THEN
    DEALLOCATE (mpitmp,stat=astat)
    DEALLOCATE (mprtmp,stat=astat)
  END IF

  wtrexist=0
  DO j = 1, ny
    DO i = 1, nx
      DO is = 1,nstyp
        IF (soiltyp(i,j,is) == 13) THEN
          wtrexist=1
          RETURN
        END IF
      END DO
    END DO
  END DO

  RETURN

  WRITE (6,'(/a,i2/a)')                                                 &
         '     Read error in surface data file '                        &
         //soilinfl//' with the I/O unit ',sfcunit,                     &
         'The model will STOP in subroutine INITSFC.'

  CALL arpsstop ("arpsstop called from initsfc reading sfc data file",1)

  WRITE (6,'(/a,i2/a)')                                                 &
         '     Read error in surface initial data file '                &
         //soilinfl//' with the I/O unit ',sfcunit,                     &
         'The model will STOP in subroutine INITSFC.'

  CALL arpsstop ("arpsstop called from initsfc reading sfc data file-2",1)

END SUBROUTINE initsfc
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE FLZERO                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE flzero(a,n)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Fill in an entire array with zeros.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/1991.
!
!  MODIFICATION HISTORY:
!
!  5/02/92 (M. Xue)
!  Added full documentation.
!
!  5/03/92 (M. Xue)
!  Further documentation.
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    n        The dimension of 1-D array a.
!
!  OUTPUT:
!
!    a        1-D array filled with zeros.
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
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, n

  REAL :: a(n)                 ! 1-D array filled with zeros
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  DO i=1,n
    a(i)=0.0
  END DO

  RETURN
END SUBROUTINE flzero
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE INITLKTB                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE initlktb
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Initializes arrays used for lookup table functions.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Adwait Sathye
!  06/02/94
!
!  MODIFICATION HISTORY:
!
!  2/17/97 (J. Zong)
!  Corrected definitions of the power functions in loop 100.
!
!  2/24/97 (Jinxing Zong, Ming Xue and Yuhe Liu)
!  Defined five pwr arrays for lookup tables to replace fractional
!  power calculations in Tao ice microphysics.
!
!  8 November 2002 (Eric Kemp)
!  Added pwr lookup tables for Ferrier fall speed options, and further
!  optimized previous tables.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  OUTPUT:
!     pwr arrays filled with evenly spaced data points
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  INCLUDE FILES
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'

  REAL :: tnw,tns,tng,roqr,roqs,roqg
  COMMON /size/ tnw,tns,tng,roqr,roqs,roqg ! Set in subroutine STCSTICE

!
!-----------------------------------------------------------------------
!
!  LOCAL VARIABLES
!
!-----------------------------------------------------------------------
!
  INTEGER :: i
  REAL :: temper

  REAL :: cpi, rhoqx, lambda ! EMK

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  DO i = 0, 1001

    pwr2046(i) =  ( 0.05 * i / 1000.0 ) ** 0.2046
    pwr525(i)  =  ( 0.05 * i / 1000.0 ) ** 0.525
    pwr875(i)  =  ( 0.05 * i / 1000.0 ) ** 0.875
    pwr1364(i) =  ( 0.05 * i / 1000000.0 ) ** 0.1364
  END DO

!-----------------------------------------------------------------------
!
!  Build a look up table for the Latent heat of vaporation, calculated
!  in the temperature range of (-100,50) degree Celcius
!  according to formula
!
!  Lv = lathv * (273.15/tbar)**(0.167+3.67E-4*tbar)
!
!  Lv for t(K) is given by
!
!  index = max( min(151, NINT(t-172.15)), 1)
!  Lv = latheatv(index)
!
!  where lathv = 2.500780e6 is set in phycst.inc, and t is the
!
!-----------------------------------------------------------------------

  DO i=-100,50
    temper =  i+273.15
    latheatv(i+101)=lathv * (273.15/temper)**(0.167+3.67E-4*temper)
  END DO

!-----------------------------------------------------------------------
!
!  Build lookup tables for T**0.81 and (rhobar/mu/psi**2)**one6th,
!  where mu and psi stand for dynamic viscosity of air and diffusivity
!  of water vapor in air, respectively. The temerature ranges from
!  -100 to 50 degree Celcius, pressure from 0 to 1500 hPa, and air
!  density from 0 to 1.5 kg/m**3 in the calculations of the tables.
!  With the above ranges of T, p and air density, (rhobar/mu/psi**2)
!  is between 0.0 and 3000.0.
!
!-----------------------------------------------------------------------

  DO i = 0, 151
    pwr81(i) = ( 173.15 + i ) ** 0.81
  END DO

  DO i = 0, 10001
    pwr1666(i) = ( 3000.0 * i / 10000.0 ) ** 0.1666667
  END DO


!-----------------------------------------------------------------------
!
!  Lookup tables for (rhobar*q)**a, where q can be qr, qs or qh.
!  rhobar is in g/cm***3. rhobar*q ranges from 0 to 50.e-6. The
!  increment of the tables is 50.e-9.
!
!-----------------------------------------------------------------------

  cpi = 4.*ATAN(1.) ! EMK
  CALL stcstice()   ! EMK Initialize Lin-Tao scheme parameters

  DO i = 0, 10001

!    pwr2    (i) =  ( 0.05 * i / 10000000.0 ) ** 0.2
!    pwr0625 (i) =  ( 0.05 * i / 10000000.0 ) ** 0.0625
!    pwr15625(i) =  ( 0.05 * i / 10000000.0 ) ** 0.15625
!    pwr105(i) = ( 0.05 * i / 10000000.0 ) ** 0.105   ! EMK
!    pwr1596(i) = ( 0.05 * i / 10000000.0 ) ** 0.1596 ! EMK

!    rhoqx = 0.05 * i / 10000000.0 ! EMK cgs units
    rhoqx = 0.05 * i * 1.0E-7 ! EMK cgs units

    pwr2    (i) =  ( rhoqx ) ** 0.2
    pwr0625 (i) =  ( rhoqx ) ** 0.0625
    pwr15625(i) =  ( rhoqx ) ** 0.15625
    pwr105(i) = ( rhoqx ) ** 0.105   ! EMK
    pwr1596(i) = ( rhoqx ) ** 0.1596 ! EMK

! EMK 18 November 2002
    IF(i == 0)THEN
     pwrlam195ratio(i) = 0.0
    ELSE
     lambda = SQRT(SQRT(cpi*roqr*tnw/rhoqx)) ! EMK cgs units
     pwrlam195ratio(i) = ( lambda + 1.95 ) ** 5 ! EMK
     pwrlam195ratio(i) = (lambda**4)/pwrlam195ratio(i)
!    pwrlam195ratio(i) = SQRT(SQRT(rhoqx/(cpi*roqr*tnw))) ! EMK cgs units
!    pwrlam195ratio(i) = ( 1. + (1.95*pwrlam195ratio(i))) ** 5 ! EMK
!    pwrlam195ratio(i) = rhoqx/pwrlam195ratio(i)
    END IF

  END DO

  RETURN
END SUBROUTINE initlktb
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GTSINLAT                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE gtsinlat(nx,ny, x,y, sinlat, xs, ys, temxy)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the sin of the lattitude of each grid point, to be used
!  in the calculation of latitude-dependent Coriolis parameters.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Ming Xue
!  3/21/95.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    x        x-coordinate of grid points in computational space (m)
!    y        y-coordinate of grid points in computational space (m)
!
!  OUTPUT:
!
!    sinlat   Sin of latitude at each grid point
!
!  WORK ARRAYS:
!
!    xs       x-coordinate at scalar points.
!    ys       y-coordinate at scalar points.
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

  INTEGER :: nx,ny          ! Dimensions of model grid.
  REAL :: x     (nx)        ! The x-coord. of the u points.
  REAL :: y     (ny)        ! The y-coord. of the v points.

  REAL :: sinlat(nx,ny)     ! Sin of latitude at each grid point

  REAL :: xs    (nx)        ! x-coord. of scalar points.
  REAL :: ys    (ny)        ! y-coord. of scalar points.

  REAL :: temxy (nx,ny)     ! Work array.

  REAL :: r2d
  INTEGER :: i,j
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO i=1,nx-1
    xs(i) = (x(i)+x(i+1))*0.5
  END DO
  xs(nx) = -xs(nx-2)+2.0*xs(nx-1)

  DO j=1,ny-1
    ys(j) = (y(j)+y(j+1))*0.5
  END DO
  ys(ny) = -ys(ny-2)+2.0*ys(ny-1)
!
  CALL xytoll(nx,ny,xs,ys,sinlat,temxy)

  r2d = ATAN(1.0)*4.0/180.0

  DO j=1,ny
    DO i=1,nx
      sinlat(i,j) = SIN( r2d * sinlat(i,j) )
    END DO
  END DO

  RETURN
END SUBROUTINE gtsinlat
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE SETPPI                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE setppi(nx,ny,nz,nt,tlvl,pprt,pbar,ppi)

  IMPLICIT NONE
  INTEGER :: nx,ny,nz,nt,tlvl
  REAL :: pprt(nx,ny,nz,nt)
  REAL :: pbar(nx,ny,nz)
  REAL :: ppi (nx,ny,nz)

  INCLUDE 'phycst.inc'

  INTEGER :: i,j,k
  REAL :: p0inv

  p0inv=1./p0
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        ppi(i,j,k)=((pprt(i,j,k,tlvl)+pbar(i,j,k))*p0inv)**rddcp
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE setppi

!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE INIT_MM                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE init_MM(nx,ny,nz,nt,tim,ptprt,ptbar,pprt,pbar,qscalar)

  USE my3mom_fncs_mod, ONLY: gammaDP

  IMPLICIT NONE

  INCLUDE  'phycst.inc'
  INCLUDE  'globcst.inc'

  INTEGER :: nx,ny,nz,nt,tim
  REAL :: pprt(nx,ny,nz,nt)
  REAL :: pbar(nx,ny,nz)
  REAL :: ppi(nx,ny,nz)     ! Exner function
  REAL :: ptprt(nx,ny,nz,nt)
  REAL :: ptbar(nx,ny,nz)
  REAL :: qscalar(nx,ny,nz,nscalar)

  REAL :: rho(nx,ny,nz)     ! Air density
  REAL :: t(nx,ny,nz)       ! Air temperature


  REAL :: temNtx            ! temporary number concentration for species x
  REAL :: temZx             ! temporary reflectivity for species x
  REAL :: Gr                ! constant in equation for radar reflectivity
  REAL :: Gi
  REAL :: Gs
  REAL :: Gg
  REAL :: Gh

  REAL, PARAMETER :: pi = 3.14159265

  ! Fixed intercept parameters for rain,snow,graupel,hail

  REAL :: N0r
  REAL :: N0s
  REAL :: N0g
  REAL :: N0h

  ! Fixed cloud number concentration

  REAL :: Ntc

  ! Fixed densities for rain,snow,graupel,hail

  REAL :: rhor
  REAL :: rhoi
  REAL :: rhos
  REAL :: rhog
  REAL :: rhoh

  REAL, PARAMETER :: epsQ = 1.0e-14

  INTEGER :: i,j,k,nq

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  N0r = n0rain
  N0s = n0snow
  N0g = n0grpl
  N0h = n0hail

  Ntc = ntcloud
  rhor = 1000.0
  rhoi = rhoice
  rhos = rhosnow
  rhog = rhogrpl
  rhoh = rhohail

  temNtx = 0.0
  temZx = 0.0

  ! Calculate exner function

  CALL setppi(nx,ny,nz,nt,tim,pprt,pbar,ppi)

  ! Calculate air temperature and density

  DO k = 1,nz-1
    DO j = 1,ny-1
      DO i = 1,nx-1
        t(i,j,k) = (ptprt(i,j,k,tim) + ptbar(i,j,k))*ppi(i,j,k)
        rho(i,j,k) = (pbar(i,j,k)+pprt(i,j,k,tim))/(rd*t(i,j,k))
      END DO
    END DO
  END DO

  ! For 2 or 3 moment scheme, diagnose the additional moments using the initial values of mixing ratios
  ! Assume alpha=0 (the shape parameter) EDIT: 08/08/08, now uses value of alpha specified in namelist

  ! Precalculate constant in radar equation

!  Gx = 20.

  Gr = ((6.+alpharain)*(5+alpharain)*(4+alpharain))/((3.+alpharain)*(2+alpharain)*(1+alpharain))
  Gi = ((6.+alphaice)*(5+alphaice)*(4+alphaice))/((3.+alphaice)*(2+alphaice)*(1+alphaice))
  Gs = ((6.+alphasnow)*(5+alphasnow)*(4+alphasnow))/((3.+alphasnow)*(2+alphasnow)*(1+alphasnow))
  Gg = ((6.+alphagrpl)*(5+alphagrpl)*(4+alphagrpl))/((3.+alphagrpl)*(2+alphagrpl)*(1+alphagrpl))
  Gh = ((6.+alphahail)*(5+alphahail)*(4+alphahail))/((3.+alphahail)*(2+alphahail)*(1+alphahail))

  IF(mphyopt >= 9) THEN

    IF((hail_ON == 0) .and. (graupel_ON == 1)) THEN
      DO k = 1,nz-1
        DO j = 1,ny-1
          DO i = 1,nx-1
             qscalar(i,j,k,P_QG) = qscalar(i,j,k,P_QG) + qscalar(i,j,k,p_QH)
             qscalar(i,j,k,P_QH) = 0.0
          END DO
        END DO
      END DO
    ELSE IF((graupel_ON == 0) .and. (hail_ON == 1) ) THEN
      DO k = 1,nz-1
        DO j = 1,ny-1
          DO i = 1,nx-1
             qscalar(i,j,k,P_QH) = qscalar(i,j,k,P_QG) + qscalar(i,j,k,p_QH)
             qscalar(i,j,k,P_QG) = 0.0
          END DO
        END DO
      END DO
    END IF

    DO nq = 1,nscalarq
      DO k = 1,nz-1
        DO j = 1,ny-1
          DO i = 1,nx-1

            temNtx = 0.0
            temZx = 0.0
            IF(qscalar(i,j,k,nq) >= epsQ) THEN
            IF(nq == P_QC) THEN   ! Set Nc to fixed cloud number concentration
              temNtx = Ntc
            ELSE IF(nq == P_QR) THEN
              temNtx = sngl(gammaDP(1.d0+dble(alpharain)))*(N0r**(3./(4.+alpharain)))* &
                       (rho(i,j,k)*qscalar(i,j,k,nq)/((pi/6.)*rhor* &
                       sngl(gammaDP(4.d0+dble(alpharain)))))**((1.+alpharain)/(4.+alpharain))
              temZx = Gr/(((pi/6.)*rhor)**2.)*((rho(i,j,k)*qscalar(i,j,k,nq))**2.)/temNtx
            ELSE IF(nq == P_QI) THEN
              temNtx = 5.*exp(0.304*(273.15-max(233.,t(i,j,k))))     ! Ice nucleation from Cooper's eqn.
              temZx = Gi/((440.0)**2.)*((rho(i,j,k)*qscalar(i,j,k,nq))**2.)/temNtx
            ELSE IF(nq == P_QS) THEN
              temNtx = sngl(gammaDP(1.d0+dble(alphasnow)))*(N0s**(3./(4.+alphasnow)))* &
                       (rho(i,j,k)*qscalar(i,j,k,nq)/((pi/6.)*rhos* &
                       sngl(gammaDP(4.d0+dble(alphasnow)))))**((1.+alphasnow)/(4.+alphasnow))
              temZx = Gs/(((pi/6.)*rhos)**2.)*((rho(i,j,k)*qscalar(i,j,k,nq))**2.)/temNtx
            ELSE IF(nq == P_QG) THEN
              temNtx = sngl(gammaDP(1.d0+dble(alphagrpl)))*(N0g**(3./(4.+alphagrpl)))* &
                       (rho(i,j,k)*qscalar(i,j,k,nq)/((pi/6.)*rhog* &
                       sngl(gammaDP(4.d0+dble(alphagrpl)))))**((1.+alphagrpl)/(4.+alphagrpl))
              temZx = Gg/(((pi/6.)*rhog)**2.)*((rho(i,j,k)*qscalar(i,j,k,nq))**2.)/temNtx
            ELSE IF(nq == P_QH) THEN
              temNtx = sngl(gammaDP(1.d0+dble(alphahail)))*(N0h**(3./(4.+alphahail)))* &
                       (rho(i,j,k)*qscalar(i,j,k,nq)/((pi/6.)*rhoh* &
                       sngl(gammaDP(4.d0+dble(alphahail)))))**((1.+alphahail)/(4.+alphahail))
              temZx = Gh/(((pi/6.)*rhoh)**2.)*((rho(i,j,k)*qscalar(i,j,k,nq))**2.)/temNtx
            END IF
            END IF
            qscalar(i,j,k,nq+nscalarq) = temNtx
            IF(mphyopt == 11) THEN
              qscalar(i,j,k,nq+nscalarq+(nscalarq-1)) = temZx
            END IF
          END DO
        END DO
      END DO
    END DO

  END IF

  RETURN
END SUBROUTINE init_MM

!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE STCSTICE                  ######
!######                                                      ######
!######                     Developed by                     ######
!######                                                      ######
!######    Goddard Cumulus Ensemble Modeling Group, NASA     ######
!######                                                      ######
!######     Center for Analysis and Prediction of Storms     ######
!######               University of Oklahoma                 ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE stcstice
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set constants used by the ice microphyscs parameterization routine
!  ICECVT
!
!  Lin et.al.  J. Clim. Appl. Meteor.  22, 1065-1092
!  Modified and coded by tao and simpson (JAS, 1989; Tao, 1993)
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: W.G. Tao, Goddard Cumulus Ensemble Modeling Group, NASA
!  01/01/1990
!
!  Modification history:
!
!  9/8/1994 (M. Xue)
!  Defined inline function CVMGP to replace external function CVMGP.
!
!  2/24/1997 (J. Zong, M. Xue and Yuhe Liu)
!  Constants rn5, rn6, rn22, rn17a, rn19a, rn20b and rn101 are
!  redefined to facilitate use of lookup tables to replace power
!  calculations for microphysical conversions and terminal velocity.
!
!  08/27/2002 (D. Weber)
!  Added fallopt option to provide a choice of fall velocity
!  computations.
!
!-----------------------------------------------------------------------
!
!  IMPLICIT NONE

  INCLUDE "globcst.inc"
  INCLUDE "phycst.inc"

  COMMON/size/ tnw,tns,tng,roqr,roqs,roqg
  COMMON/cont/ c76,c358,c172,c409,c218,c580,c610,c149,c879,c141
  COMMON/b3cs/ ag,bg,as,bs,aww,bww,bgh,bgq,bsh,bsq,bwh,bwq
  COMMON/bterv/ zrc,zgc,zsc,vrc,vgc,vsc
  COMMON/bsnw/ alv,alf,als,t0,t00,avc,afc,asc,rn1,bnd1,rn2,bnd2,        &
      rn3,rn4,rn5,rn6,rn7,rn8,rn9,rn10,rn101,rn10a,rn11,rn11a,          &
      rn12,rn12a(31),rn12b(31),rn13(31),rn14,rn15,rn15a,rn16,rn17,      &
      rn17a,rn17b,rn17c,rn18,rn18a,rn19,rn19a,rn19b,rn20,rn20a,rn20b,   &
      bnd3,rn21,rn22,rn23,rn23a,rn23b,rn25,rn25a(31),rn30a,rn30b,       &
      rn30c,rn31,beta,rn32
!
  REAL ::  a1(31),a2(31)
  DATA a1/.7939E-7,.7841E-6,.3369E-5,.4336E-5,.5285E-5,.3728E-5,        &
      .1852E-5,.2991E-6,.4248E-6,.7434E-6,.1812E-5,.4394E-5,.9145E-5,   &
      .1725E-4,.3348E-4,.1725E-4,.9175E-5,.4412E-5,.2252E-5,.9115E-6,   &
      .4876E-6,.3473E-6,.4758E-6,.6306E-6,.8573E-6,.7868E-6,.7192E-6,   &
      .6513E-6,.5956E-6,.5333E-6,.4834E-6/
  DATA a2/.4006,.4831,.5320,.5307,.5319,.5249,.4888,.3894,.4047,        &
      .4318,.4771,.5183,.5463,.5651,.5813,.5655,.5478,.5203,.4906,      &
      .4447,.4126,.3960,.4149,.4320,.4506,.4483,.4460,.4433,.4413,      &
      .4382,.4361/

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
!  Define the density and size distribution of precipitation
!
!-----------------------------------------------------------------------
!
  roqr = 1.

!  tnw = .08
!  roqs = .1
!  tns = .03
!  roqg = .913
!  tng = .0004
!
! Assign intercept and density parameters passed in through common
! block in phycst.inc and convert to cgs unit.

  tnw = n0rain*1.0e-8
  tns = n0snow*1.0e-8
  tng = n0hail*1.0e-8

  roqs = rhosnow * 0.001
  roqg = rhohail * 0.001

!  WRITE(6,*) n0rain,n0snow, n0hail
!
!  WRITE(0,*) '***Inside stcstice,tnw,roqs,tns,roqg,tng=',tnw,roqs,tns,roqg,tng


  cpi = 4.*ATAN(1.)
  cpi2 = cpi*cpi
  grvt = 980.
  tca = 2.43E3
  dwv = .226
  dva = 1.718E-4
  amw = 18.016
  ars = 8.314E7
  scv = 2.2904487
  t0 = 273.16
  t00 = 238.16
  alv = 2.5E10
  alf = 3.336E9
  als = 2.8336E10

!%%%
  cplocal = 1.003E7  ! Missing in the original Tao's code
!%%%

  avc = alv/cplocal
  afc = alf/cplocal
  asc = als/cplocal
  rw = 4.615E6
  cw = 4.187E7
  ci = 2.093E7
  c76 = 7.66
  c358 = 35.86
  c172 = 17.26939
  c409 = 4098.026
  c218 = 21.87456
  c580 = 5807.695
  c610 = 6.1078E3
  c149 = 1.496286E-5
  c879 = 8.794142
  c141 = 1.4144354E7

!
!-----------------------------------------------------------------------
!
!  Define the coefficients used in terminal velocity
!
!-----------------------------------------------------------------------
!

  ag = 1400.
  bg = .5
  as = 152.93
  bs = .25
  aww= 2115.
  bww= .8

!     print *,fallopt
      IF(fallopt.eq.2)THEN   ! compute Ferrier version of the fall vel.
!       ferrier numbers for rain.
        aww= 4854.
        bww= 1.0
!       ferrier numbers for snow.
        as = 129.6
        bs = 0.42
!       ferrier numbers for graupel/hail.
        ag = 1094.3
        bg = 0.6384
      END IF


  bgh = .5*bg
  bsh = .5*bs
  bwh = .5*bww
  bgq = .25*bg
  bsq = .25*bs
  bwq = .25*bww


  ga3 = 2.
  ga4 = 6.
  ga5 = 24.
  ga6 = 120.
  ga7 = 720.
  ga8 = 5040.
  ga9 = 40320.
  ga3b = 4.6941552

  IF(fallopt.eq.2)THEN   ! Ferrier formulation...
    ga4b = 4.0
  ELSE                   ! original Lin configuration
    ga4b = 17.83779
  END IF



  ga4b = 17.83779
  ga6b = 496.6041
  ga5bh = 1.827363
  ga4g = 11.63177
  ga3g = 3.3233625
  ga5gh = 1.608355

  IF (bg == 0.37) THEN
    ga4g = 9.730877
    ga3g = 2.8875
    ga5gh = 1.526425
  END IF

  ga3d = 2.54925
  ga4d = 8.285063
  ga5dh = 1.456943

  IF (bs == 0.57) THEN
    ga3d = 3.59304
    ga4d = 12.82715
    ga5dh = 1.655588
  END IF

  IF(bs == 0.11) THEN
    ga3d = 2.218906
    ga4d = 6.900796
    ga5dh = 1.382792
  END IF
!
!-----------------------------------------------------------------------
!
!  Lin et al., 1983
!
!-----------------------------------------------------------------------
!
  ac1 = aww
  bc1 = bww
  cc1 = as
  dc1 = bs
  cd1 = 6.e-1
  cd2 = 4.*grvt/(3.*cd1)
  zrc = (cpi*roqr*tnw)**0.25
  zsc = (cpi*roqs*tns)**0.25
  zgc = (cpi*roqg*tng)**0.25
  vrc = ac1*ga4b/(6.*zrc**bww)
  vsc = cc1*ga4d/(6.*zsc**bs)
  vgc = ga4g*SQRT(cd2*roqg/zgc)/6.


!  Added by D. Weber  8/27/02
!  for rain we have a more complicated modification
!  vtr =  4854*sqrt(rho0/rho)* [gamma(5)*lambda**(4)/
!                               gamma(4)*(lambda+f)**(5)
!  gamma (5) = 5-1! = 24, gamma (4) = 4-1! = 4

   IF (fallopt == 2) THEN
!    ferrier equation for rain partial calc for lambda....
     vrc = 19416.0/SQRT(SQRT(cpi*roqr*tnw)) ! EMK 19 November 2002
     vsc = 225.12/((cpi*roqs*tns)**0.105)
!    ferrier equation for hail
     vgc = 2577.6419/((cpi*roqg*tng)**.1596)
   END IF

!  print *,ag,ga4d,zgc,bg
!  print *,'vrc,vsc,vgc= ',vrc,vsc,vgc


  rn1 = 1.e-3
  bnd1 = 6.e-4
  rn2 = 1.e-3
  bnd2 = 1.e-3
  rn3 = .25*cpi*tns*cc1*ga3d
  esw = 1.
  rn4 = .25*cpi*esw*tns*cc1*ga3d
  eri = 1.
  rn5 = .25*cpi*eri*tnw*ac1*ga3b / zrc**bww
  ami = 1./(24.*4.19E-10)
  rn6 = cpi2*eri*tnw*ac1*roqr*ga6b*ami / zrc**bww
  esr = 1.
  rn7 = cpi2*esr*tnw*tns*roqs
  rn8 = cpi2*esr*tnw*tns*roqr
  rn9 = cpi2*tns*tng*roqs
  rn10 = 2.*cpi*tns
  rn101 = .31*ga5dh*SQRT(cc1) / zsc**0.625
  rn10a = als*als/rw
  rn11 = 2.*cpi*tns/alf
  rn11a = cw/alf
  ami50 = 4.8E-7
  ami40 = 3.84E-9
  eiw = 1.
  ui50 = 100.
  ri50 = 5.e-3
  cmn = 1.05E-15
  rn12 = cpi*eiw*ui50*ri50**2

  DO k = 1,31

    y1 = 1.-a2(k)
    rn13(k) = a1(k)*y1/(ami50**y1-ami40**y1)
    rn12a(k) = rn13(k)/ami50
    rn12b(k) = a1(k)*ami50**a2(k)
    rn25a(k) = a1(k)*cmn**a2(k)

  END DO

  egw = 1.
  rn14 = .25*cpi*egw*tng*ga3g*SQRT(cd2*roqg)
  egi = .1
  rn15 = .25*cpi*egi*tng*ga3g*SQRT(cd2*roqg)
  egi = 1.
  rn15a = .25*cpi*egi*tng*ga3g*SQRT(cd2*roqg)
  egr = 1.
  rn16 = cpi2*egr*tng*tnw*roqr
  rn17 = 2.*cpi*tng
  rn17a = .31*ga5gh*(cd2*roqg)**.25 * zgc**0.25
  rn17b = cw-ci
  rn17c = cw
  apri = .66
  bpri = 1.e-4
  rn18 = 20.*cpi2*bpri*tnw*roqr
  rn18a = apri
  rn19 = 2.*cpi*tng/alf
  rn19a = .31*ga5gh*(cd2*roqg)**.25 * zgc**0.25
  rn19b = cw/alf
  rn20 = 2.*cpi*tng
  rn20a = als*als/rw
  rn20b = .31*ga5gh*(cd2*roqg)**.25 * zgc**0.25
  bnd3 = 2.e-3
  rn21 = 1.e3*1.569E-12/0.15
  erw = 1.
  rn22 = .25*cpi*erw*ac1*tnw*ga3b / zrc**bww
  rn23 = 2.*cpi*tnw
  rn23a = .31*ga5bh*SQRT(ac1)
  rn23b = alv*alv/rw
!  cn0 = 1.e-8
  cn0 = 1.e-5
  rn25 = cn0/1000.
  rn30a = alv*als*amw/(tca*ars)
  rn30b = alv/tca
  rn30c = ars/(dwv*amw)
  rn31 = 1.e-17
  beta = -.6
  rn32 = 4.*51.545E-4

  RETURN
END SUBROUTINE stcstice
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CCINIT (TINA)              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE ccinit(nx,ny,nz,nts,x,y,z,zp,cc,tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Initialize the scalar concentration for bubble initialization in the
!  middle of a restart run.
!
!  Note: This is the exact copy from Tina except for formatting.
!        It may need checking for MPI validality. (by Y. Wang on 4/25/2012)
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Tina Chow
!  DATE: 2/22/10
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nts      Number of time levels to be initialized.
!
!    cc       Scalar concentration (-)
!
!    x        x-coordinate of grid points in computational space (m)
!    y        y-coordinate of grid points in computational space (m)
!    z        z-coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space(m)
!
!
!  OUTPUT:
!
!    u        x-component of velocity at all time levels (m/s).
!    v        y-component of velocity at all time levels (m/s).
!    w        z-component of velocity at all time levels (m/s).
!    ptprt    Perturbation potential temperature at all time levels
!             (K)
!    pprt     Perturbation pressure at all time levels (Pascal)
!    cc       Pollutant concentration (-) !michi
!    qv       Water vapor specific humidity at all time levels
!             (kg/kg)
!    qc       Cloud water mixing ratio at all time levels (kg/kg)
!    qr       Rainwater mixing ratio at all time levels (kg/kg)
!    qi       Cloud ice mixing ratio at all time levels (kg/kg)
!    qs       Snow mixing ratio at all time levels (kg/kg)
!    qh       Hail mixing ratio at all time levels (kg/kg)
!
!    ptcumsrc Source term in pt-equation due to cumulus parameterization
!    qcumsrc Source term in water equations due to cumulus parameterization
!    raing    Grid supersaturation rain
!    rainc    Cumulus convective rain
!    prcrate  Precipitation ratesrain
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!
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


  INTEGER :: nx,ny,nz          ! The number of grid points in 3
                               ! directions
  INTEGER :: nts               ! Number of time levels to be initialized.
  INTEGER :: tpast             ! Index of time level for the past time.
  INTEGER :: tpresent          ! Index of time level for the present
                               ! time.
  INTEGER :: tfuture           ! Index of time level for the future
                               ! time.

  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal).
  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific
                               ! humidity(kg/kg)

  REAL :: x     (nx)           ! x-coord. of the physical and compu-
                               ! tational grid. Defined at u-point.
  REAL :: y     (ny)           ! y-coord. of the physical and compu-
                               ! tational grid. Defined at v-point.
  REAL :: z     (nz)           ! z-coord. of the computational grid.
                               ! Defined at w-point on the staggered
                               ! grid.
  REAL :: zp    (nx,ny,nz)     ! Physical height coordinate defined at
                               ! w-point of the staggered grid.

  REAL :: cc    (nx,ny,nz,nts) ! Pollutant concentration (-) !michi

  REAL :: tem1(nx,ny,nz)       ! Temporary work array

!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  REAL :: xs, ys, zs
  REAL :: us, vs, ws, rhobar
  REAL :: radnd , pi,pi2,pi4
  INTEGER :: i,j,k, n, ip
  INTEGER :: iseed,ibgn,iend,jbgn,jend,kbgn,kend
  INTEGER :: ebc1,wbc1,nbc1,sbc1

  REAL :: amplitud
  REAL :: knumx,lnumy,mnumz
  REAL :: lnthx,lnthy,lnthz
  REAL :: lambda,lambdah,lambda2

!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'phycst.inc'
  INCLUDE 'bndry.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.

  INTEGER :: nxlg, nylg, ilg,jlg
  INTEGER :: istatus
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF( nts > 1 ) THEN
    tpast=1
    tpresent=2
    tfuture=3
  ELSE
    tpast=1
    tpresent=1
    tfuture=1
  END IF
!
!
!-----------------------------------------------------------------------
! TINA
!  Specify the initial concentration field. (granvold)
!
!-----------------------------------------------------------------------
!TINA test
!          IF (myproc == 0) print *, 'Init cc here'
!TINA test


!      IF (mp_opt > 0) CALL mpbarrier !is this needed?

  IF (ccin == 1 .AND. (cpoint == -2 .OR. cpoint == -3)) THEN
     IF (cpoint == -2) THEN
       IF (myproc == 0) print *, 'Initializing uniform cc bubble in init3d'
       !TINA use ccemit to hold max concentration - note this only
       !works for one initial bubble so use ccemit(1)
       IF (myproc == 0) print *, '  Concentration = ', ccemit(1)
     ELSE
       IF (myproc == 0) print *, 'Initializing non-uniform cc bubble in init3d'
       IF (myproc == 0) print *, '  Maximum concentration = ', ccemit(1)
     END IF

     ip = 1
     ! similar to pt bubble prescribed above but for cc here
     DO k= 1,nz-1
       DO j= 1,ny-1
         DO i= 1,nx-1

           xs = (x(i)+x(i+1))*0.5
           ys = (y(j)+y(j+1))*0.5
           zs = (zp(i,j,k)+zp(i,j,k+1))*0.5
           IF ( ABS(zs-pt0ctrz(ip))/pt0radz(ip) < 1 ) THEN

             IF( pt0rady(ip) < 0.0 .OR. runmod == 2 ) THEN
                                  ! 2-d bubble in x-z plane.

               radnd = SQRT( ((xs-pt0ctrx(ip))/pt0radx(ip))**2          &
                           + ((zs-pt0ctrz(ip))/pt0radz(ip))**2 )

             ELSE IF( pt0radx(ip) < 0.0 .OR. runmod == 3 ) THEN
                                  ! 2-d bubble in y-z plane.

               radnd = SQRT( ((ys-pt0ctry(ip))/pt0rady(ip))**2          &
                           + ((zs-pt0ctrz(ip))/pt0radz(ip))**2 )

             ELSE                         ! 3-d bubble

               radnd = SQRT( ((xs-pt0ctrx(ip))/pt0radx(ip))**2 +        &
                             ((ys-pt0ctry(ip))/pt0rady(ip))**2 +        &
                             ((zs-pt0ctrz(ip))/pt0radz(ip))**2 )
             END IF

             IF(radnd >= 1.0) THEN
               cc(i,j,k,1) = 0.0
             ELSE

               IF (cpoint == -2) THEN
                  cc(i,j,k,1) = ccemit(1)
               ELSE
                  cc(i,j,k,1) = ccemit(1)*(COS(pi2*radnd )**2)

               END IF

             END IF

           END IF

         END DO
       END DO
     END DO

     IF (myproc == 0) print *, '  Maximum value initialized is ', maxval(cc)

  ELSE IF (ccin == 1 .AND. cpoint == -4) THEN

!TINA test
     IF (myproc == 0) print *, 'Init cc here 2'
!TINA test

!            IF (myproc == 0) print *, 'Initializing uniform cc rectangle/cuboid in init3d'
!            IF (myproc == 0) print *, '  Concentration = ', ccemit(1)
!            IF (myproc == 0) print *, '  pt0radx = ', pt0radx
!            IF (myproc == 0) print *, '  pt0rady = ', pt0rady
!            IF (myproc == 0) print *, '  pt0radz = ', pt0radz
!            IF (myproc == 0) print *, '  pt0ctrx = ', pt0ctrx
!            IF (myproc == 0) print *, '  pt0ctry = ', pt0ctry
!            IF (myproc == 0) print *, '  pt0ctrz = ', pt0ctrz

     ip = 1
     DO k= 1,nz-1
       DO j= 1,ny-1
         DO i= 1,nx-1

           xs = (x(i)+x(i+1))*0.5
           ys = (y(j)+y(j+1))*0.5
           zs = (zp(i,j,k)+zp(i,j,k+1))*0.5
           IF ( ABS(zs-pt0ctrz(ip))/pt0radz(ip) < 1 ) THEN

             IF( pt0rady(ip) < 0.0 .OR. runmod == 2 ) THEN
                                          ! 2-d rectangle in x-z plane.

               IF (      abs(xs-pt0ctrx(ip)) <= pt0radx(ip)             &
                   .AND. abs(zs-pt0ctrz(ip)) <= pt0radz(ip)) THEN
                   radnd = 0.0
               ELSE
                   radnd = 1.0
               END IF

             ELSE IF( pt0radx(ip) < 0.0 .OR. runmod == 3 ) THEN
                                          ! 2-d rectangle in y-z plane.

               IF (      abs(ys-pt0ctry(ip)) <= pt0rady(ip)             &
                   .AND. abs(zs-pt0ctrz(ip)) <= pt0radz(ip)) THEN
                   radnd = 0.0
               ELSE
                   radnd = 1.0
               END IF

             ELSE                         ! 3-d cuboid

               IF (      abs(xs-pt0ctrx(ip)) <= pt0radx(ip)             &
                   .AND. abs(ys-pt0ctry(ip)) <= pt0rady(ip)             &
                   .AND. abs(zs-pt0ctrz(ip)) <= pt0radz(ip)) THEN
                   radnd = 0.0
               ELSE
                   radnd = 1.0
               END IF

             END IF

             IF(radnd >= 1.0) THEN
               cc(i,j,k,1) = 0.0
             ELSE
               cc(i,j,k,1) = ccemit(1)
             END IF

           END IF

         END DO
       END DO
     END DO

! print *,'TINA testing after loop'
!      IF (mp_opt > 0) CALL mpbarrier !is this needed?

!if (myproc == 1) write(990,*) 'myproc = ',myproc,', max(cc) after loop in ccinit =',maxval(cc)
!if (myproc == 0) write(991,*) 'myproc = ',myproc,', max(cc) after loop in ccinit =',maxval(cc)

     ELSE IF (ccin == 1 .AND. cpoint == -5) THEN  ! user may specify field here
           ! specify a cylinder for BM cases
       IF (myproc == 0) print *, 'Initializing uniform cc cylinder in init3d'
       IF (myproc == 0) print *, '  Concentration = ', ccemit(1)

       ip = 1
       DO k= 1,nz-1
         DO j= 1,ny-1
           DO i= 1,nx-1

             xs = (x(i)+x(i+1))*0.5
             ys = (y(j)+y(j+1))*0.5
             zs = (zp(i,j,k)+zp(i,j,k+1))*0.5


             IF( pt0rady(ip) < 0.0 .OR. runmod == 2 ) THEN
                                          ! 2-d bubble in x-z plane.

               radnd = SQRT( ((xs-pt0ctrx(ip))/pt0radx(ip))**2  )

             ELSE IF( pt0radx(ip) < 0.0 .OR. runmod == 3 ) THEN
                                          ! 3-d bubble in y-z plane.

               radnd = SQRT( ((ys-pt0ctry(ip))/pt0rady(ip))**2  )

             ELSE                         ! 3-d bubble

               radnd = SQRT( ((xs-pt0ctrx(ip))/pt0radx(ip))**2 +        &
                             ((ys-pt0ctry(ip))/pt0rady(ip))**2 )
             END IF

             IF (zs <= pt0radz(ip) .AND. radnd < 1.0) THEN
                cc(i,j,k,1) = ccemit(1)
                IF (myproc == 0) print *, '--> ', xs, ys, radnd, zs
             ELSE
                cc(i,j,k,1) = 0.0
             END IF

           END DO
         END DO
       END DO

     END IF

!if (myproc == 0) write(990,*) 'myproc = ',myproc,', max(cc) midway in ccinit =',maxval(cc)
!if (myproc == 1) write(991,*) 'myproc = ',myproc,', max(cc) midway in ccinit =',maxval(cc)

  ebc1=0
  wbc1=0
  sbc1=0
  nbc1=0

  IF( ebc == 1 .OR.ebc == 2 .OR. ebc == 3 )  ebc1=ebc
  IF( wbc == 1 .OR.wbc == 2 .OR. wbc == 3 )  wbc1=wbc
  IF( sbc == 1 .OR.sbc == 2 .OR. sbc == 3 )  sbc1=sbc
  IF( nbc == 1 .OR.nbc == 2 .OR. nbc == 3 )  nbc1=nbc

  IF (mp_opt > 0) THEN
    CALL acct_interrupt(mp_acct)
    CALL mpsendrecv2dew(cc,nx,ny,nz,ebc1,wbc1,0,tem1)
    CALL mpsendrecv2dns(cc,nx,ny,nz,nbc1,sbc1,0,tem1)
  END IF
  CALL acct_interrupt(bc_acct)
  CALL bcsclr(nx,ny,nz,dtbig,                                           &
              cc(1,1,1,1),cc(1,1,1,1),                                  &
              cc(1,1,1,1),tem1,tem1,tem1,tem1,                          &
              ebc1,wbc1,nbc1,sbc1,tbc,bbc,                              &
              ebc_global,wbc_global,nbc_global,sbc_global)
  CALL acct_stop_inter

!TINA copy to other times
  DO n=1,nts
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          cc(i,j,k,n)= cc(i,j,k,1)
        END DO
      END DO
    END DO
  END DO
!michi

! IF (mp_opt > 0) CALL mpbarrier !is this needed?
! print *, 'myproc = ',myproc,', max(cc) inside ccinit =',maxval(cc)

  RETURN
END SUBROUTINE ccinit
