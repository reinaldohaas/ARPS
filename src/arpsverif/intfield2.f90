!File created 27 March 2000.  "Official" routines are in intfield.f...EMK
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE SETDRVY2                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE setdrvy2(nx,ny,nz,                                           &
           ibeg,iend,jbeg,jend,kbeg,kend,                               &
           dyfld,rdyfld,var,                                            &
           slopey,alphay,betay)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Calculate the coefficients of interpolating polynomials
!    in the y-direction.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, CAPS, November, 1996
!
!  MODIFICATION HISTORY:
!  Eric Kemp, November 1999
!  Added checks for missing data (-9999.)
!
!
!-----------------------------------------------------------------------
!
!  INPUT:
!    nx       Number of model grid points in the x-direction (east/west)
!    ny       Number of model grid points in the y-direction (north/south)
!    nz       Number of model grid points in the vertical
!
!    ibeg,iend   Range of x index to do interpolation
!    jbeg,jend   Range of y index to do interpolation
!    kbeg,kend   Range of z index to do interpolation
!
!    dyfld    Vector of delta-y (m) of field to be interpolated
!    rdyfld   Vector of 1./delta-y (1/m) of field to be interpolated
!
!    var      variable to be interpolated
!
!    slopey   Piecewise linear df/dy
!    alphay   Coefficient of y-squared term in y quadratic interpolator
!    betay    Coefficient of y term in y quadratic interpolator
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  INTEGER :: ibeg,iend,jbeg,jend,kbeg,kend
  REAL :: dyfld(ny)
  REAL :: rdyfld(ny)
  REAL :: var(nx,ny,nz)
  REAL :: slopey(nx,ny,nz)
  REAL :: alphay(nx,ny,nz)
  REAL :: betay(nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  INTEGER :: jstart,jstop
  REAL :: rtwody
  REAL,PARAMETER :: miss = -9999.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  jstart=MAX(jbeg,2)
  jstop=MIN((jend-1),(ny-2))
  DO k=kbeg,kend
    DO j=jstart,jstop
      DO i=ibeg,iend
        slopey(i,j,k)=(var(i,j+1,k)-var(i,j,k))*rdyfld(j)
        IF (var(i,j+1,k) == miss .OR. var(i,j,k) == miss) THEN
          slopey(i,j,k) = miss
        END IF
        rtwody=1./(dyfld(j-1)+dyfld(j))
        alphay(i,j,k)=((var(i,j+1,k)-var(i,j,k))*rdyfld(j) +            &
                 (var(i,j-1,k)-var(i,j,k))*rdyfld(j-1))*rtwody
        IF (var(i,j+1,k) == miss .OR. var(i,j,k) == miss .OR.           &
            var(i,j-1,k) == miss) THEN
          alphay(i,j,k) = miss
        END IF
        betay(i,j,k)=(var(i,j+1,k)-var(i,j,k))*rdyfld(j) -              &
                   dyfld(j)*alphay(i,j,k)
        IF (var(i,j+1,k) == miss .OR. var(i,j,k) == miss .OR.           &
            alphay(i,j,k) == miss) THEN
          betay(i,j,k) = miss
        END IF
      END DO
    END DO
  END DO
  RETURN
END SUBROUTINE setdrvy2

!
!##################################################################
!##################################################################
!######                                                      ######
!######                 FUNCTION PNTINT2D2                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

  FUNCTION pntint2d2(vnx,vny,                                           &
             ivbeg,ivend,jvbeg,jvend,                                   &
             iorder,vx,vy,xpnt,ypnt,iloc,jloc,var,                      &
             dxfld,dyfld,rdxfld,rdyfld,                                 &
             slopey,alphay,betay)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    Interpolate a 2-d field for a single point on that plane.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, CAPS, November, 1996
!
!  MODIFICATION HISTORY:
!
!  Eric Kemp, November 1999
!  Added checks for missing data.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!    vnx       Number of model grid points in the x-direction (east/west)
!    vny       Number of model grid points in the y-direction (north/south)
!
!    ivbeg,ivend   Range of x index to use in verification array
!    jvbeg,jvend   Range of y index to use in verification array
!
!    iorder   Interpolation parameter.
!             iorder specifies the order of interpolation
!             1 = bi-linear
!             2 = bi-quadratic
!
!    vx       x coordinate of verif scalar grid points in physical space (m)
!    vy       y coordinate of verif scalar grid points in physical space (m)
!
!    xpnt     x coordinate (m) of interpolation point
!    ypnt     y coordinate (m) of interpolation point
!
!    iloc     I-index of interpolation point in field to be interpolated
!    jloc     J-index of interpolation point in field to be interpolated
!    dxfld    Vector of delta-x (m) of field to be interpolated
!    dyfld    Vector of delta-y (m) of field to be interpolated
!    rdxfld   Vector of 1./delta-x (1/m) of field to be interpolated
!    rdyfld   Vector of 1./delta-y (1/m) of field to be interpolated
!
!    slopey   Piecewise linear df/dy
!    alphay   Coefficient of y-squared term in y quadratic interpolator
!    betay    Coefficient of y term in y quadratic interpolator
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  REAL :: pntint2d2
  INTEGER :: vnx,vny
  INTEGER :: ivbeg,ivend,jvbeg,jvend
  INTEGER :: iorder
  REAL :: vx(vnx)
  REAL :: vy(vny)
  REAL :: xpnt
  REAL :: ypnt
  INTEGER :: iloc
  INTEGER :: jloc
  REAL :: var(vnx,vny)
  REAL :: dxfld(vnx)
  REAL :: dyfld(vny)
  REAL :: rdxfld(vnx)
  REAL :: rdyfld(vny)
  REAL :: slopey(vnx,vny)
  REAL :: alphay(vnx,vny)
  REAL :: betay(vnx,vny)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: ii,jj
  REAL :: delx,dely
  REAL :: alpha,beta,rtwodx
  REAL :: varm1,var00,varp1
  REAL :: varint
  REAL,PARAMETER :: miss = -9999.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  Compute bilinear interpolated value
!
!-----------------------------------------------------------------------
!
  IF(iorder == 1) THEN
    ii=MIN(MAX(iloc,ivbeg),(ivend-1))
    jj=MIN(MAX(jloc,jvbeg),(jvend-1))
    delx=(xpnt-vx(ii))
    dely=(ypnt-vy(jj))
    varint=(1.-delx*rdxfld(ii))*                                        &
               (var(ii  ,jj)+slopey(ii  ,jj)*dely)+                     &
               (delx*rdxfld(ii))*                                       &
               (var(ii+1,jj)+slopey(ii+1,jj)*dely)

    IF (var(ii,jj) == miss .OR. slopey(ii,jj) ==  miss .OR.             &
        var(ii+1,jj) == miss .OR. slopey(ii+1,jj) == miss) THEN
      varint = miss
    END IF

!    WRITE(6,*)
!    WRITE(6,*)'EMK: ii = ',ii,' jj = ',jj
!    WRITE(6,*)'EMK: delx = ',delx,' dely = ',dely
!    WRITE(6,*)'EMK: rdxfld(ii) = ',rdxfld(ii)
!    WRITE(6,*)'EMK: var(ii,jj) = ',var(ii,jj),' var(ii+1,jj) = ',       &
!                    var(ii+1,jj)
!    WRITE(6,*)'EMK: slopey(ii,jj) = ',slopey(ii,jj),                    &
!                    ' slopey(ii+1,jj) = ',slopey(ii+1,jj)
!    WRITE(6,*)'EMK: varint = ',varint
!
!-----------------------------------------------------------------------
!
!  Compute biquadratic
!
!-----------------------------------------------------------------------
!
  ELSE
    ii=MIN(MAX(iloc,(ivbeg+1)),(ivend-1))
    jj=MIN(MAX(jloc,(jvbeg+1)),(jvend-1))
    delx=(xpnt-vx(ii))
    dely=(ypnt-vy(jj))
!
!-----------------------------------------------------------------------
!
!    Stencil is ii-1 to ii+1 and jj-1 to jj + 1
!
!    Interpolate in y.
!
!-----------------------------------------------------------------------
!
    varm1=(alphay(ii-1,jj)*dely+betay(ii-1,jj))*dely+var(ii-1,jj)
    IF (alphay(ii-1,jj) == miss .OR. betay(ii-1,jj) == miss .OR.        &
        var(ii-1,jj) == miss) varm1 = miss

    var00=(alphay(ii  ,jj)*dely+betay(ii  ,jj))*dely+var(ii  ,jj)
    IF (alphay(ii,jj) == miss .OR. betay(ii,jj) == miss .OR.            &
        var(ii,jj) == miss) var00 = miss

    varp1=(alphay(ii+1,jj)*dely+betay(ii+1,jj))*dely+var(ii+1,jj)
    IF (alphay(ii+1,jj) == miss .OR. betay(ii+1,jj) == miss .OR.        &
        var(ii,jj) == miss) varp1 = miss
!
!-----------------------------------------------------------------------
!
!    Interpolate intermediate results in x.
!
!-----------------------------------------------------------------------
!
    rtwodx=1./(dxfld(ii-1)+dxfld(ii))
    alpha=((varp1-var00)*rdxfld(ii  ) +                                 &
           (varm1-var00)*rdxfld(ii-1))*rtwodx
    IF (varp1 == miss .OR. var00 == miss .OR. varm1 == miss) alpha = miss
    beta=(varp1-var00)*rdxfld(ii) -                                     &
             dxfld(ii)*alpha
    IF (varp1 == miss .OR. var00 == miss .OR. alpha == miss) beta = miss
    varint=(alpha*delx+beta)*delx+var00
    IF (alpha == miss .OR. beta == miss .OR. var00 == miss) varint = miss
!        WRITE(6,*)
!        WRITE(6,*)'EMK: varint = ',varint
!        WRITE(6,*)'EMK: var00 = ',var00,' delx = ',delx
!        WRITE(6,*)'EMK: alpha = ',alpha,' beta = ',beta
!        WRITE(6,*)'EMK: ii = ',ii,' jj = ',jj
!        WRITE(6,*)'EMK: dely = ',dely
!        WRITE(6,*)'EMK: xpnt = ',xpnt,' ypnt = ',ypnt
!        WRITE(6,*)'EMK: dxfld(ii) = ',dxfld(ii),' dxfld(ii-1) = ',      &
!                   dxfld(ii-1)
  END IF
  pntint2d2=varint
  RETURN
  END FUNCTION pntint2d2

