!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE AAMULT                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE aamult(a,b,nx,ny,nz,                                         &
           ibgn,iend,jbgn,jend,kbgn,kend, ab)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the element-wise product of arrays a and b.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/91
!
!  MODIFICATION HISTORY:
!
!  5/05/92 (M. Xue)
!  Added full documentation.
!
!  6/01/2 (K. Brewster)
!  Further facelift.
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    a        First multiplier array
!    b        Second multiplier array
!
!    nx       First dimension of arrays a, b and ab
!    ny       Second dimension of arrays a, b and ab
!    nz       Third dimension of arrays a, b and ab
!
!    ibgn     i-index where multiplication begins.
!    iend     i-index where multiplication ends.
!    jbgn     j-index where multiplication begins.
!    jend     j-index where multiplication ends.
!    kbgn     k-index where multiplication begins.
!    kend     k-index where multiplication ends.
!
!  OUTPUT:
!
!    ab       Element-wise product (a*b) over range specified by
!             the starting and ending indices given above.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions

  REAL :: a (nx,ny,nz)         ! Input array 1
  REAL :: b (nx,ny,nz)         ! Input array 2

  INTEGER :: ibgn,iend,jbgn,jend,kbgn,kend

  REAL :: ab(nx,ny,nz)         ! Product array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
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

  DO k=kbgn,kend
    DO j=jbgn,jend
      DO i=ibgn,iend

        ab(i,j,k)=a(i,j,k)*b(i,j,k)

      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE aamult

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE AVGX                       ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE avgx(a, onvf, nx,ny,nz,                                      &
           ibgn,iend,jbgn,jend,kbgn,kend,                               &
           aavg)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Perform a spatial average operation on array (a) in the x direction.
!  The average for the two input values is defined at the midpoint
!  between the two values.
!  If the averaged variable aavg is on the grid volume face, onvf =1
!  otherwise, onvf = 0
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/91
!
!  MODIFICATION HISTORY:
!
!  5/05/92 (M. Xue)
!  Added full documentation.
!
!  6/01/2 (K. Brewster)
!  Further facelift.
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    a        Input array
!
!    onvf     Integer grid point indicator
!             If the averaged variable aavg is on the grid volume face,
!             onvf =1; otherwise, onvf = 0
!
!    nx       First dimension of arrays a and aavg
!    ny       Second dimension of arrays a and aavg
!    nz       Third dimension of arrays a and aavg
!
!    ibgn     i-index where multiplication begins.
!    iend     i-index where multiplication ends.
!    jbgn     j-index where multiplication begins.
!    jend     j-index where multiplication ends.
!    kbgn     k-index where multiplication begins.
!    kend     k-index where multiplication ends.
!
!  OUTPUT:
!
!    aavg     Result of average operation on array
!             (a) in the x direction
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions

  REAL :: a   (nx,ny,nz)       ! Input array
  INTEGER :: onvf              ! Integer grid point indicator

  INTEGER :: ibgn,iend,jbgn,jend,kbgn,kend
                                ! Integer indicator of multiplication
  REAL :: aavg(nx,ny,nz)       ! Averaged array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,iright,ileft
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  IF ( onvf == 1) THEN

    iright = 0
    ileft = -1

  ELSE

    iright = 1
    ileft = 0

  END IF

  DO k=kbgn,kend
    DO j=jbgn,jend
      DO i=ibgn,iend

        aavg(i,j,k)=(a(i+iright,j,k) + a(i+ileft ,j,k))*0.5

      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE avgx

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE AVG2X                      ######
!######                                                      ######
!######                     Developed by                     ######
!######  Center for the Analysis and Prediction of Storms    ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE avg2x(a, nx,ny,nz,                                           &
           ibgn,iend,jbgn,jend,kbgn,kend, aavg)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Perform a spatial average operation on array a in the x direction.
!  Averaging is over the interval 2*deltax.
!  The averaged variable is defined at the midpoint between the two
!  points for the two input values.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:  David E. Jahn
!  9/30/93
!
!  MODIFICATION HISTORY:
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    a        Input array
!
!    nx       First dimension of arrays a and aavg
!    ny       Second dimension of arrays a and aavg
!    nz       Third dimension of arrays a and aavg
!
!    ibgn     Index in first dimension to begin multiplication
!    iend     Index in first dimension to end multiplication
!    jbgn     Index in second dimension to begin multiplication
!    jend     Index in second dimension to end multiplication
!    kbgn     Index in third dimension to begin multiplication
!    kend     Index in third dimension to end multiplication
!
!  OUTPUT:
!
!    aavg     Result of average operation on array a in x direction
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions

  REAL :: a   (nx,ny,nz)       ! Input array

  INTEGER :: ibgn,iend,jbgn,jend,kbgn,kend
                                ! Integer indicator of multiplication

  REAL :: aavg(nx,ny,nz)       ! Averaged array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
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

  DO k=kbgn,kend
    DO j=jbgn,jend
      DO i=ibgn,iend

        aavg(i,j,k)=(a(i-1,j,k)                                         &
                    +a(i+1,j,k))*0.5

      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE avg2x

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE AVGY                       ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE avgy(a, onvf, nx,ny,nz,                                      &
           ibgn,iend,jbgn,jend,kbgn,kend, aavg)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Perform a spatial average operation on array a in the y direction.
!  The average for the two input values is defined at the midpoint
!  between the two values.
!  If the averaged variable aavg is on the grid volume face, onvf =1
!  otherwise, onvf = 0
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/91
!
!  MODIFICATION HISTORY:
!
!  5/05/92 (M. Xue)
!  Added full documentation.
!
!  6/01/2 (K. Brewster)
!  Further facelift.
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    a        Input array
!
!    nx       First dimension of arrays a and aavg
!    ny       Second dimension of arrays a and aavg
!    nz       Third dimension of arrays a and aavg
!
!    onvf     Integer grid point indicator
!             If the averaged variable aavg is on the grid volume face,
!             onvf =1; otherwise, onvf = 0
!
!    ibgn     i-index where multiplication begins.
!    iend     i-index where multiplication ends.
!    jbgn     j-index where multiplication begins.
!    jend     j-index where multiplication ends.
!    kbgn     k-index where multiplication begins.
!    kend     k-index where multiplication ends.
!
!  OUTPUT:
!
!    aavg     Result of average operation on array (a) in the
!             y direction
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions

  REAL :: a   (nx,ny,nz)       ! Input array
  INTEGER :: onvf              ! Integer grid point indicator

  INTEGER :: ibgn,iend,jbgn,jend,kbgn,kend
                                ! Integer indicator of multiplication
  REAL :: aavg(nx,ny,nz)       ! Average array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,jright,jleft
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF ( onvf == 1) THEN

    jright = 0
    jleft = -1

  ELSE

    jright = 1
    jleft = 0

  END IF

  DO k=kbgn,kend
    DO j=jbgn,jend
      DO i=ibgn,iend

        aavg(i,j,k)=(a(i,j+jright,k)+a(i,j+jleft,k))*0.5

      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE avgy

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE AVG2Y                      ######
!######                                                      ######
!######                     Developed by                     ######
!######  Center for the Analysis and Prediction of Storms    ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE avg2y(a, nx,ny,nz,                                           &
           ibgn,iend,jbgn,jend,kbgn,kend, aavg)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Perform a spatial average operation on array a in the y direction.
!  Averaging is over a distance 2y*delta-y.
!  The averaged variable is defined at the midpoint between the two
!  points for the two input values.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:  David E. Jahn
!  9/30/93
!
!  MODIFICATION HISTORY:
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    a        Input array
!
!    nx       First dimension of arrays a and aavg
!    ny       Second dimension of arrays a and aavg
!    nz       Third dimension of arrays a and aavg
!
!    ibgn     Index in first dimension to begin operation
!    iend     Index in first dimension to end operation
!    jbgn     Index in second dimension to begin operation
!    jend     Index in second dimension to end operation
!    kbgn     Index in third dimension to begin operation
!    kend     Index in third dimension to end operation
!
!  OUTPUT:
!
!    aavg     Result of average operation on array a in y direction
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions

  REAL :: a   (nx,ny,nz)       ! Input array

  INTEGER :: ibgn,iend,jbgn,jend,kbgn,kend
                                ! Integer indicator of multiplication
  REAL :: aavg(nx,ny,nz)       ! Average array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
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
  DO k=kbgn,kend
    DO j=jbgn,jend
      DO i=ibgn,iend

        aavg(i,j,k)=(a(i,j+1,k)+a(i,j-1,k))*0.5

      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE avg2y


!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE AVGZ                       ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE avgz(a, onvf, nx,ny,nz,                                      &
           ibgn,iend,jbgn,jend,kbgn,kend, aavg)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Perform a spatial average operation on array (a) in the z direction.
!  The average of the two input values is defined at the midpoint
!  between the two values.
!  If the averaged variable aavg is on the grid volume face, onvf =1
!  otherwise, onvf = 0
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/91
!
!  MODIFICATION HISTORY:
!
!  5/05/92 (M. Xue)
!  Added full documentation.
!
!  6/01/2 (K. Brewster)
!  Further facelift.
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    a        Input array
!
!    nx       First dimension of arrays a and aavg
!    ny       Second dimension of arrays a and aavg
!    nz       Third dimension of arrays a and aavg
!
!    onvf     Integer grid point indicator
!             If the averaged variable aavg is on the grid volume
!             face,onvf =1; otherwise, onvf = 0
!
!    ibgn     i-index where multiplication begins.
!    iend     i-index where multiplication ends.
!    jbgn     j-index where multiplication begins.
!    jend     j-index where multiplication ends.
!    kbgn     k-index where multiplication begins.
!    kend     k-index where multiplication ends.
!
!  OUTPUT:
!
!    aavg     Result of average operation on array (a) in the z direction
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions

  REAL :: a   (nx,ny,nz)       ! Input array
  INTEGER :: onvf              ! Integer grid point indicator

  INTEGER :: ibgn,iend,jbgn,jend,kbgn,kend
                                ! Integer indicator of multiplication
  REAL :: aavg(nx,ny,nz)       ! Average array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,kup,kdown
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF ( onvf == 1) THEN

    kup = 0
    kdown = -1

  ELSE

    kup = 1
    kdown = 0

  END IF

  DO k=kbgn,kend
    DO j=jbgn,jend
      DO i=ibgn,iend

        aavg(i,j,k)=(a(i,j,k+kup)+a(i,j,k+kdown ))*0.5

      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE avgz

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE AVG2Z                      ######
!######                                                      ######
!######                     Developed by                     ######
!######  Center for the Analysis and Prediction of Storms    ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE avg2z(a, nx,ny,nz,                                           &
           ibgn,iend,jbgn,jend,kbgn,kend, aavg)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Perform a spatial average operation on array a in the z direction.
!  Averaging is over interval 2*deltaz.
!  The averaged variable is defined at the midpoint between the two
!  points for the two input values.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:  David E. Jahn
!  9/30/93
!
!  MODIFICATION HISTORY:
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    a        Input array
!
!    nx       First dimension of arrays a and aavg
!    ny       Second dimension of arrays a and aavg
!    nz       Third dimension of arrays a and aavg
!
!    ibgn     Index in first dimension to begin operation
!    iend     Index in first dimension to end operation
!    jbgn     Index in second dimension to begin operation
!    jend     Index in second dimension to end operation
!    kbgn     Index in third dimension to begin operation
!    kend     Index in third dimension to end operation
!
!  OUTPUT:
!
!    aavg     Result of average operation on array a in z direction
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions

  REAL :: a   (nx,ny,nz)       ! Input array

  INTEGER :: ibgn,iend,jbgn,jend,kbgn,kend
                                ! Integer indicator of multiplication
  REAL :: aavg(nx,ny,nz)       ! Average array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
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
  DO k=kbgn,kend
    DO j=jbgn,jend
      DO i=ibgn,iend

        aavg(i,j,k)=(a(i,j,k+1)+a(i,j,k-1))*0.5

      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE avg2z



!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE DIFX                       ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE difx(a, onvf, nx,ny,nz,                                      &
           ibgn,iend,jbgn,jend,kbgn,kend, dx, adifx)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Perform a finite difference operation on an array (a) in the x
!  direction: adifx = d( a ) / dx. The output variable is defined at
!  the midpoint between the two points whose values are differenced
!  (i.e., the output and input arrays are staggered). The output
!  variable adifx is defined on the grid volume face when onvf =1;
!  otherwise, onvf = 0.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/91
!
!  MODIFICATION HISTORY:
!
!  5/05/92 (M. Xue)
!  Added full documentation.
!
!  6/01/2 (K. Brewster)
!  Further facelift.
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    a        Input array
!
!    onvf     Integer grid point indicator
!             If the averaged variable aavg is on the grid volume
!             face, onvf =1; otherwise, onvf = 0
!
!    nx       First dimension of arrays a and adifx
!    ny       Second dimension of arrays a and adifx
!    nz       Third dimension of arrays a and adfix
!
!    ibgn     i-index where multiplication begins.
!    iend     i-index where multiplication ends.
!    jbgn     j-index where multiplication begins.
!    jend     j-index where multiplication ends.
!    kbgn     k-index where multiplication begins.
!    kend     k-index where multiplication ends.
!
!    dx       Grid spacing in x direction (m)
!
!  OUTPUT:
!
!    adifx    Differenced array del(a)/delx
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions

  REAL :: a    (nx,ny,nz)      ! Input array
  INTEGER :: onvf              ! Integer grid point indicator

  INTEGER :: ibgn,iend,jbgn,jend,kbgn,kend
                                ! Integer indicator of multiplication
  REAL :: dx                   ! Grid spacing in x direction (m)
  REAL :: adifx(nx,ny,nz)      ! Differenced array del(a)/delx
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,iright,ileft
  REAL :: dxinv
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF ( onvf == 1) THEN

    iright = 0
    ileft = -1

  ELSE

    iright = 1
    ileft = 0

  END IF

  dxinv = 1.0/dx

  DO k=kbgn,kend
    DO j=jbgn,jend
      DO i=ibgn,iend

        adifx(i,j,k)=(a(i+iright,j,k)-a(i+ileft ,j,k))*dxinv

      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE difx

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE DIF2X                      ######
!######                                                      ######
!######                     Developed by                     ######
!######  Center for the Analysis and Prediction of Storms    ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE dif2x(a, nx,ny,nz,                                           &
           ibgn,iend,jbgn,jend,kbgn,kend, dx, adifx)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Perform a centered finite difference operation on an array a over
!  2 grid distance in the x direction: adifx = d( a ) / dx.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: David E. Jahn
!  4/19/93
!
!  MODIFICATION HISTORY:
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    a        Input array
!
!    nx       First dimension of arrays a and adifx
!    ny       Second dimension of arrays a and adifx
!    nz       Third dimension of arrays a and adfix
!
!    ibgn     Index in first dimension to begin operation
!    iend     Index in first dimension to end operation
!    jbgn     Index in second dimension to begin operation
!    jend     Index in second dimension to end operation
!    kbgn     Index in third dimension to begin operation
!    kend     Index in third dimension to end operation
!
!    dx       Grid spacing in x direction (m)
!
!  OUTPUT:
!
!    adifx    Differenced array del(a)/delx
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions

  REAL :: a    (nx,ny,nz)      ! Input array

  INTEGER :: ibgn,iend,jbgn,jend,kbgn,kend
                                ! Integer indicator of multiplication
  REAL :: dx                   ! Grid spacing in x direction (m)
  REAL :: adifx(nx,ny,nz)      ! Differenced array del(a)/delx

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  REAL :: dxinv
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  dxinv = 0.5/dx

  DO k=kbgn,kend
    DO j=jbgn,jend
      DO i=ibgn,iend

        adifx(i,j,k)=(a(i+1,j,k)-a(i-1,j,k))*dxinv

      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE dif2x



!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE DIFY                       ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE dify(a, onvf, nx,ny,nz,                                      &
           ibgn,iend,jbgn,jend,kbgn,kend, dy, adify)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Perform a finite difference operation on an array (a) in the y
!  direction: adify = d( a ) / dy. The output variable is defined at
!  the midpoint between the two points whose values are differenced
!  (i.e., the output and input arrays are staggered). The output
!  variable adifx is defined on the grid volume face when onvf =1;
!  otherwise, onvf = 0.

!  direction. adify = d( a ) / dy. The output variable is defined at
!  the midpoint between the two points whose values are differenced.
!  The output variable adify is defined on the grid volume face when
!  onvf =1; otherwise, onvf = 0.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/91
!
!  MODIFICATION HISTORY:
!
!  5/05/92 (M. Xue)
!  Added full documentation.
!
!  6/01/2 (K. Brewster)
!  Further facelift.
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    a        Input array
!
!    onvf     Integer grid point indicator
!             If the averaged variable aavg is on the grid volume
!             face, onvf =1; otherwise, onvf = 0
!
!    nx       First dimension of arrays a and adify
!    ny       Second dimension of arrays a and adify
!    nz       Third dimension of arrays a and adify
!
!    ibgn     i-index where multiplication begins.
!    iend     i-index where multiplication ends.
!    jbgn     j-index where multiplication begins.
!    jend     j-index where multiplication ends.
!    kbgn     k-index where multiplication begins.
!    kend     k-index where multiplication ends.
!
!    dy       Grid spacing in y direction (m)
!
!  OUTPUT:
!
!    adify    Differenced array del(a)/dely
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions

  REAL :: a   (nx,ny,nz)       ! Input array
  INTEGER :: onvf              ! Integer grid point indicator

  INTEGER :: ibgn,iend,jbgn,jend,kbgn,kend
                                ! Integer indicator of multiplication
  REAL :: dy                   ! Grid spacing in y direction (m)
  REAL :: adify(nx,ny,nz)      ! Differenced array del(a)/dely
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,jright,jleft
  REAL :: dyinv
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF ( onvf == 1) THEN

    jright = 0
    jleft = -1

  ELSE

    jright = 1
    jleft = 0

  END IF

  dyinv = 1.0/dy

  DO k=kbgn,kend
    DO j=jbgn,jend
      DO i=ibgn,iend

        adify(i,j,k)=(a(i,j+jright,k)-a(i,j+jleft,k))*dyinv

      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE dify

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE DIF2Y                      ######
!######                                                      ######
!######                     Developed by                     ######
!######  Center for the Analysis and Prediction of Storms    ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE dif2y(a, nx,ny,nz,                                           &
           ibgn,iend,jbgn,jend,kbgn,kend, dy, adify)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Perform a centered finite difference operation on an array a over
!  2 grid distance in the y direction: adify = d( a ) / dy.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: David E. Jahn
!  4\19\93
!
!  MODIFICATION HISTORY:
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    a        Input array
!
!    nx       First dimension of arrays a and adify
!    ny       Second dimension of arrays a and adify
!    nz       Third dimension of arrays a and adify
!
!    ibgn     Index in first dimension to begin operation
!    iend     Index in first dimension to end operation
!    jbgn     Index in second dimension to begin operation
!    jend     Index in second dimension to end operation
!    kbgn     Index in third dimension to begin operation
!    kend     Index in third dimension to end operation
!
!    dy       Grid spacing in y direction (m)
!
!  OUTPUT:
!
!    adify    Differenced array del(a)/dely
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions

  REAL :: a   (nx,ny,nz)       ! Input array

  INTEGER :: ibgn,iend,jbgn,jend,kbgn,kend
                                ! Integer indicator of multiplication
  REAL :: dy                   ! Grid spacing in y direction (m)
  REAL :: adify(nx,ny,nz)      ! Differenced array del(a)/dely
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  REAL :: dyinv
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  dyinv = 0.5/dy

  DO k=kbgn,kend
    DO j=jbgn,jend
      DO i=ibgn,iend

        adify(i,j,k)=(a(i,j+1,k)-a(i,j-1,k))*dyinv

      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE dif2y


!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE DIFZ                       ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE difz(a, onvf, nx,ny,nz,                                      &
           ibgn,iend,jbgn,jend,kbgn,kend, dz, adifz)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Perform a finite difference operation on an array (a) in the z
!  direction: adifz = d( a ) / dz. The output variable is defined at
!  the midpoint between the two points whose values are differenced
!  (i.e., the output and input arrays are staggered). The output
!  variable adifx is defined on the grid volume face when onvf =1;
!  otherwise, onvf = 0.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/91
!
!  MODIFICATION HISTORY:
!
!  5/05/92 (M. Xue)
!  Added full documentation.
!
!  6/01/2 (K. Brewster)
!  Further facelift.
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    a        Input array
!
!    nx       First dimension of arrays a and adifz
!    ny       Second dimension of arrays a and adifz
!    nz       Third dimension of arrays a and adifz
!
!    onvf     Integer grid point indicator
!             If the averaged variable aavg is on the grid volume
!             face, onvf =1; otherwise, onvf = 0
!
!    ibgn     i-index where multiplication begins.
!    iend     i-index where multiplication ends.
!    jbgn     j-index where multiplication begins.
!    jend     j-index where multiplication ends.
!    kbgn     k-index where multiplication begins.
!    kend     k-index where multiplication ends.
!
!    dz       Grid spacing in z direction (m)
!
!  OUTPUT:
!
!    adifz    Differenced array del(a)/delz
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions

  REAL :: a    (nx,ny,nz)      ! Input array
  INTEGER :: onvf              ! Integer grid point indicator

  INTEGER :: ibgn,iend,jbgn,jend,kbgn,kend
                                ! Integer indicator of multiplication
  REAL :: dz                   ! Grid spacing in z direction (m)
  REAL :: adifz(nx,ny,nz)      ! Differenced array del(a)/delz
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,kup,kdown
  REAL :: dzinv
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF ( onvf == 1) THEN

    kup = 0
    kdown = -1

  ELSE

    kup = 1
    kdown = 0

  END IF

  dzinv = 1.0/dz

  DO k=kbgn,kend
    DO j=jbgn,jend
      DO i=ibgn,iend

        adifz(i,j,k)=(a(i ,j,k+kup)-a(i ,j,k+kdown ))*dzinv

      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE difz
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE DIF2Z                      ######
!######                                                      ######
!######                     Developed by                     ######
!######  Center for the Analysis and Prediction of Storms    ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE dif2z(a, nx,ny,nz,                                           &
           ibgn,iend,jbgn,jend,kbgn,kend, dz, adifz)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Perform a centered finite difference operation on an array a over
!  2 grid distance in the z direction: adifz = d( a ) / dz.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  11/6/95
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    a        Input array
!
!    nx       First dimension of arrays a and adifz
!    ny       Second dimension of arrays a and adifz
!    nz       Third dimension of arrays a and adfiz
!
!    ibgn     Index in first dimension to begin operation
!    iend     Index in first dimension to end operation
!    jbgn     Index in second dimension to begin operation
!    jend     Index in second dimension to end operation
!    kbgn     Index in third dimension to begin operation
!    kend     Index in third dimension to end operation
!
!    dz       Grid spacing in z direction (m)
!
!  OUTPUT:
!
!    adifz    Differenced array del(a)/delz
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions

  REAL :: a    (nx,ny,nz)      ! Input array

  INTEGER :: ibgn,iend,jbgn,jend,kbgn,kend
                                ! Integer indicator of multiplication
  REAL :: dz                   ! Grid spacing in z direction (m)
  REAL :: adifz(nx,ny,nz)      ! Differenced array del(a)/delz

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  REAL :: dzinv
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  dzinv = 0.5/dz

  DO k=kbgn,kend
    DO j=jbgn,jend
      DO i=ibgn,iend

        adifz(i,j,k)=(a(i,j,k+1)-a(i,j,k-1))*dzinv

      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE dif2z

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE DIFXX                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE difxx(a, nx,ny,nz,                                           &
           ibgn,iend,jbgn,jend,kbgn,kend, dx, adifxx)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the second order difference of an array (a) in the
!  x direction.  The operator is defined as del**2(a)/delx**2.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/91
!
!  MODIFICATION HISTORY:
!
!  5/05/92 (M. Xue)
!  Added full documentation.
!
!  6/01/2 (K. Brewster)
!  Further facelift.
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    a        Input array
!
!    nx       First dimension of arrays a and adifxx
!    ny       Second dimension of arrays a and adifxx
!    nz       Third dimension of arrays a and adifxx
!
!    ibgn     i-index where multiplication begins.
!    iend     i-index where multiplication ends.
!    jbgn     j-index where multiplication begins.
!    jend     j-index where multiplication ends.
!    kbgn     k-index where multiplication begins.
!    kend     k-index where multiplication ends.
!
!    dx       Grid spacing in x direction (m)
!
!  OUTPUT:
!
!    adifxx   Differenced array del**2(a)/delx**2
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions

  REAL :: a    (nx,ny,nz)      ! Input array

  INTEGER :: ibgn,iend,jbgn,jend,kbgn,kend
                                ! Integer indicator of multiplication
  REAL :: dx                   ! Grid spacing in x direction (m)
  REAL :: adifxx(nx,ny,nz)     ! Differenced array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  REAL :: dxinv2
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
!  NOTE:
!
!  The order of calculation in the following formula should not
!  be changed. The order affects the machine trunction error
!  in the result.
!
!-----------------------------------------------------------------------
!
  dxinv2 = 1.0/(dx*dx)

  DO k=kbgn,kend
    DO j=jbgn,jend
      DO i=ibgn,iend
        adifxx(i,j,k)=((a(i+1,j,k)-a(i,j,k))-(a(i,j,k)-a(i-1,j,k)))     &
                     *dxinv2

      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE difxx

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE DIFYY                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE difyy(a, nx,ny,nz,                                           &
           ibgn,iend,jbgn,jend,kbgn,kend, dy, adifyy)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the second order finite difference of an array (a) in the
!  y direction.  The operator is defined as del**2(a)/dely**2.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/91
!
!  MODIFICATION HISTORY:
!
!  5/05/92 (M. Xue)
!  Added full documentation.
!
!  6/01/2 (K. Brewster)
!  Further facelift.
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    a        input array
!
!    nx       First dimension of arrays a and adifyy
!    ny       Second dimension of arrays a and adifyy
!    nz       Third dimension of arrays a and adifyy
!
!    ibgn     i-index where multiplication begins.
!    iend     i-index where multiplication ends.
!    jbgn     j-index where multiplication begins.
!    jend     j-index where multiplication ends.
!    kbgn     k-index where multiplication begins.
!    kend     k-index where multiplication ends.
!
!    dy       grid spacing in y direction (m)
!
!  OUTPUT:
!
!    adifyy   difference array del**2(a)/dely**2
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx, ny, nz         ! Number of grid points in 3 directions

  REAL :: a     (nx,ny,nz)      ! Input array

  INTEGER :: ibgn,iend,jbgn,jend,kbgn,kend
                                 ! Integer indicator of multiplication
  REAL :: dy                    ! Grid spacing in y direction (m)
  REAL :: adifyy(nx,ny,nz)      ! Differenced array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  REAL :: dyinv2
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  dyinv2 = 1.0/(dy*dy)

  DO k=kbgn,kend
    DO j=jbgn,jend
      DO i=ibgn,iend
!
!-----------------------------------------------------------------------
!
!  NOTE:
!
!  The order of calculation in the following formula should not
!  be changed. The order affects the machine trunction error
!  in the result.
!
!-----------------------------------------------------------------------
!
        adifyy(i,j,k)=((a(i,j+1,k)-a(i,j,k))-(a(i,j,k)-a(i,j-1,k)))     &
                     *dyinv2

      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE difyy

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE DIFZZ                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE difzz(a, nx,ny,nz,                                           &
           ibgn,iend,jbgn,jend,kbgn,kend, dz, adifzz)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the second order finite difference of an array (a) in the
!  z direction.  The operator is defined as del**2(a)/delz**2.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/91
!
!  MODIFICATION HISTORY:
!
!  5/05/92 (M. Xue)
!  Added full documentation.
!
!  6/01/2 (K. Brewster)
!  Further facelift.
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    a        Input array
!
!    nx       First dimension of arrays a and adifzz
!    ny       Second dimension of arrays a and adifzz
!    nz       Third dimension of arrays a and adifzz
!
!    ibgn     i-index where multiplication begins.
!    iend     i-index where multiplication ends.
!    jbgn     j-index where multiplication begins.
!    jend     j-index where multiplication ends.
!    kbgn     k-index where multiplication begins.
!    kend     k-index where multiplication ends.
!
!    dz       Grid spacing in the vertical direction (m)
!
!  OUTPUT:
!
!    adifzz   Differenced array del**2(a)/delz**2
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx, ny, nz         ! Number of grid points in 3 directions

  REAL :: a     (nx,ny,nz)      ! Input array

  INTEGER :: ibgn,iend,jbgn,jend,kbgn,kend
                                 ! Integer indicator of multiplication
  REAL :: dz                    ! Grid spacing in the vertical
                                ! direction (m)
  REAL :: adifzz(nx,ny,nz)      ! Output difference array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  REAL :: dzinv2
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  NOTE:
!
!  The order of calculation in the following formula should not
!  be changed. The order affects the machine trunction error
!  in the result.
!
!-----------------------------------------------------------------------
!
  dzinv2 = 1.0/(dz*dz)

  DO k=kbgn,kend
    DO j=jbgn,jend
      DO i=ibgn,iend

        adifzz(i,j,k)=((a(i,j,k+1)-a(i,j,k))-(a(i,j,k)-a(i,j,k-1)))     &
                     *dzinv2

      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE difzz

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE AVGSU                      ######
!######                                                      ######
!######                     Developed by                     ######
!######  Center for the Analysis and Prediction of Storms    ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE avgsu(s,nx,ny,nz,jbgn,jend,kbgn,kend,su,tem)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Average scalar array s to u points, up to the x boundary.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR:
!  2/15/93 (M. Xue and H. Jin)
!
!  MODIFICATION HISTORY:
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!  2000/09/15 (Gene Bassett)
!  Added extra temporary array to argument list for MP option.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    s        An array defined at the scalar point.
!
!    nx       First dimension of arrays a and aavg
!    ny       Second dimension of arrays a and aavg
!    nz       Third dimension of arrays a and aavg
!
!    jbgn     Starting point for j computation
!    jend     Ending point for j computation
!    kbgn     Starting point for k computation
!    kend     Ending point for k computation
!
!  OUTPUT:
!
!    su       An array averaged from array 's' defined at the
!             scalar point to the u-point.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions

  INTEGER :: jbgn,jend         ! Domain of j computations
  INTEGER :: kbgn,kend         ! Domain of k computations

  REAL :: s (nx,ny,nz)         ! Input array
  REAL :: su(nx,ny,nz)         ! Averaged array

  REAL :: tem(nx,ny,nz)        ! Temporary work array for MP option
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: onvf

  INTEGER :: astat
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
  onvf = 1
  CALL avgx(s, onvf, nx,ny,nz, 2,nx-1, jbgn,jend, kbgn,kend, su)

  IF (mp_opt > 0) THEN
    CALL acct_interrupt(mp_acct)
    CALL mpsendrecv2dew(su,nx,ny,nz,ebc,wbc,1,tem)
  END IF
  CALL acct_interrupt(bc_acct)
  CALL bcsu(nx,ny,nz,jbgn,jend,kbgn,kend,ebc,wbc,su)
  CALL acct_stop_inter

  RETURN
END SUBROUTINE avgsu
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE AVGSV                      ######
!######                                                      ######
!######                     Developed by                     ######
!######  Center for the Analysis and Prediction of Storms    ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE avgsv(s,nx,ny,nz,ibgn,iend,kbgn,kend,sv,tem)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Average scalar array s to v points, up to the y boundary.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR:
!  2/15/93 (M. Xue and H. Jin)
!
!  MODIFICATION HISTORY:
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!  2000/09/15 (Gene Bassett)
!  Added extra temporary array to argument list for MP option.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    s        An array defined at the scalar point.
!
!    nx       First dimension of arrays a and aavg
!    ny       Second dimension of arrays a and aavg
!    nz       Third dimension of arrays a and aavg
!
!    ibgn     Starting point for i computation
!    iend     Ending point for i computation
!    kbgn     Starting point for k computation
!    kend     Ending point for k computation
!
!  OUTPUT:
!
!    su       An array averaged from array 's' defined at the
!             scalar point to the v-point.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions

  INTEGER :: ibgn,iend         ! Domain of i computations
  INTEGER :: kbgn,kend         ! Domain of k computations

  REAL :: s (nx,ny,nz)         ! Input array
  REAL :: sv(nx,ny,nz)         ! Averaged array
  REAL :: tem(nx,ny,nz)        ! Temporary work array for MP option
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: onvf

  INTEGER :: astat
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
  onvf = 1
  CALL avgy(s, onvf,                                                    &
            nx,ny,nz, ibgn,iend, 2,ny-1, kbgn,kend, sv)

  IF (mp_opt > 0) THEN
    CALL acct_interrupt(mp_acct)
    !CALL mpsend2dew(sv,nx,ny,nz,2,mptag,tem)
    !CALL mprecv2dew(sv,nx,ny,nz,2,mptag,tem)
    CALL mpsendrecv2dns(sv,nx,ny,nz,nbc,sbc,2,tem)
  END IF

  CALL acct_interrupt(bc_acct)
  CALL bcsv(nx,ny,nz,ibgn,iend,kbgn,kend,nbc,sbc,sv)
  CALL acct_stop_inter

  RETURN
END SUBROUTINE avgsv
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE AVGSW                      ######
!######                                                      ######
!######                     Developed by                     ######
!######  Center for the Analysis and Prediction of Storms    ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE avgsw(s,nx,ny,nz,ibgn,iend,jbgn,jend,sw)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Average scalar array s to w points, up to k=1, and nz.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR:
!  2/15/93 (M. Xue and H. Jin)
!
!  MODIFICATION HISTORY:
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    s        An array defined at the scalar point.
!
!    nx       First dimension of arrays a and aavg
!    ny       Second dimension of arrays a and aavg
!    nz       Third dimension of arrays a and aavg
!
!    ibgn     Starting point for i computation
!    iend     Ending point for i computation
!    jbgn     Starting point for j computation
!    jend     Ending point for j computation
!
!  OUTPUT:
!
!    sw       An array averaged from array 's' defined at the
!             scalar point to the w-point.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions

  INTEGER :: ibgn,iend         ! Domain of i computation
  INTEGER :: jbgn,jend         ! Domain of j computation

  REAL :: s (nx,ny,nz)         ! Input array
  REAL :: sw(nx,ny,nz)         ! Averaged array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: onvf
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
  onvf = 1
  CALL avgz(s, onvf,                                                    &
            nx,ny,nz, ibgn,iend, jbgn,jend, 2,nz-1, sw)

  CALL acct_interrupt(bc_acct)
  CALL bcsw(nx,ny,nz,ibgn,iend,jbgn,jend,tbc,bbc,sw)
  CALL acct_stop_inter

  RETURN
END SUBROUTINE avgsw
