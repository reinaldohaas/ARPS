!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CHECKSX                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE checksx(s, nx,ny,nz, ch, tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Check the west-east symmetry of a field defined at a scalar point
!  on the y-z plane containing the maximum magnitude of that field.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue & Hao Jin
!  4/20/93.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    s        An array field defined at the scalar point.
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!  OUTPUT:
!
!    None.
!
!  WORK ARRAY:
!
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

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: s   (nx,ny,nz)       ! An array field defined at the
                               ! scalar point
  CHARACTER (LEN=*) :: ch

  REAL :: tem1(nx,ny,nz)
  INTEGER :: i,j,k
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,(nx-1)/2+1
        tem1(i,j,k)=s( nx-i,j,k) - s(i,j,k)
      END DO
    END DO
  END DO

  j = (ny-1)/2+1

  CALL wrigar(tem1,1,nx,1,ny,1,nz,1,(nx-1)/2,j,j,1,nz-1,ch,0.0,2)

  RETURN
END SUBROUTINE checksx

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CHECKUX                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE checkux(u, nx,ny,nz, ch, tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Check the west-east symmetry of a field defined at a u-point
!  on the y-z plane containing the maximum magnitude of that field.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  11/20/91
!
!  MODIFICATION HISTORY:
!
!  5/06/92 (M. Xue)
!  Added full documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    u        An array field defined at the u grid point.
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!  OUTPUT:
!
!    None.
!
!  WORK ARRAY:
!
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

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: u   (nx,ny,nz)       ! An array field defined at the u
                               ! grid point
  CHARACTER (LEN=*) :: ch

  REAL :: tem1(nx,ny,nz)
  INTEGER :: i,j,k
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,(nx-1)/2+1
        tem1(i,j,k)=u( nx -i+1, j,k) + u(i, j,k)
      END DO
    END DO
  END DO

  j = (ny-1)/2+1

  CALL wrigar(tem1,1,nx,1,ny,1,nz,1,(nx-1)/2+1,j,j,1,nz-1,ch,.0,2)


  RETURN
END SUBROUTINE checkux

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CHECKVX                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE checkvx(v, nx,ny,nz, ch, tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Check the west-east symmetry of a field defined at a v-point
!  on the y-z plane containing the maximum magnitude of that field.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  11/20/91.
!
!  MODIFICATION HISTORY:
!
!  5/06/92 (M. Xue)
!  Added full documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    v        An array field defined at the v grid point.
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!  OUTPUT:
!
!    None.
!
!  WORK ARRAY:
!
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

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: v   (nx,ny,nz)       ! An array field defined at the v
                               ! grid point
  CHARACTER (LEN=*) :: ch

  REAL :: tem1(nx,ny,nz)
  INTEGER :: i,j,k
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO k=1,nz-1
    DO j=1,ny
      DO i=1,(nx-1)/2+1
        tem1(i,j,k)=v( nx-i,j,k) - v(i,j,k)
      END DO
    END DO
  END DO

  j = (ny-1)/2+1

  CALL wrigar(tem1,1,nx,1,ny,1,nz,1,(nx-1)/2,j,j,1,nz-1,ch,0.0,2)

  RETURN
END SUBROUTINE checkvx

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CHECKWX                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE checkwx(w, nx,ny,nz, ch, tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Check the west-east symmetry of a field defined at a w-point
!  on the y-z plane containing the maximum magnitude of that field.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  11/20/91.
!
!  MODIFICATION HISTORY:
!
!  5/06/92 (M. Xue)
!  Added full documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    w        An array field defined at the w grid point.
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!  OUTPUT:
!
!    None.
!
!  WORK ARRAY:
!
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

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: w   (nx,ny,nz)       ! An array field defined at the w
                               ! grid point
  CHARACTER (LEN=*) :: ch

  REAL :: tem1(nx,ny,nz)
  INTEGER :: i,j,k
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  DO k=1,nz
    DO j=1,ny-1
      DO i=1,(nx-1)/2+1
        tem1(i,j,k)=w( nx-i,j,k) - w(i,j,k)
      END DO
    END DO
  END DO

  j = (ny-1)/2+1

  CALL wrigar(tem1,1,nx,1,ny,1,nz,1,(nx-1)/2,j,j,1,nz,ch,.0,2)

  RETURN
END SUBROUTINE checkwx
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CHECKSY                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE checksy(s, nx,ny,nz, ch, tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Check the south-north symmetry of a field defined at a scalar point
!  on the x-z plane containing the maximum magnitude of that field.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue & Hao Jin
!  4/20/93.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    s        An array field defined at the scalar point.
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!  OUTPUT:
!
!    None.
!
!  WORK ARRAY:
!
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

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: s   (nx,ny,nz)       ! An array field defined at the
                               ! scalar point
  CHARACTER (LEN=*) :: ch

  REAL :: tem1(nx,ny,nz)
  INTEGER :: i,j,k
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO k=1,nz-1
    DO i=1,nx-1
      DO j=1,(ny-1)/2+1
        tem1(i,j,k)=s( i,ny-j,k) - s(i,j,k)
      END DO
    END DO
  END DO

  i = (nx-1)/2+1

  CALL wrigar(tem1,1,nx,1,ny,1,nz,i,i,1,(ny-1)/2,1,nz-1,ch,0.0,2)

  RETURN
END SUBROUTINE checksy

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CHECKUY                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE checkuy(u, nx,ny,nz, ch, tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Check the south-north symmetry of a field defined at a u-point
!  on the x-z plane containing the maximum magnitude of that field.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue & Hao Jin
!  4/20/93.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    u        An array field defined at the u grid point.
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!  OUTPUT:
!
!    None.
!
!  WORK ARRAY:
!
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

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: u   (nx,ny,nz)       ! An array field defined at the u
                               ! grid point
  CHARACTER (LEN=*) :: ch

  REAL :: tem1(nx,ny,nz)
  INTEGER :: i,j,k
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  DO k=1,nz-1
    DO i=1,nx-1
      DO j=1,(ny-1)/2+1
        tem1(i,j,k)=u( i,ny-j,k) - u(i, j,k)
      END DO
    END DO
  END DO

  i = (nx-1)/2+1

  CALL wrigar(tem1,1,nx,1,ny,1,nz,i,i,1,(ny-1)/2,1,nz-1,ch,.0,2)


  RETURN
END SUBROUTINE checkuy

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CHECKVY                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE checkvy(v, nx,ny,nz, ch, tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Check the south-north symmetry of a field defined at a v-point
!  on the x-z plane containing the maximum magnitude of that field.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue & Hao Jin
!  4/20/93.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    v        An array field defined at the v grid point.
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!  OUTPUT:
!
!    None.
!
!  WORK ARRAY:
!
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

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: v   (nx,ny,nz)       ! An array field defined at the v
                               ! grid point
  CHARACTER (LEN=*) :: ch

  REAL :: tem1(nx,ny,nz)
  INTEGER :: i,j,k
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO k=1,nz-1
    DO i=1,nx-1
      DO j=1,(ny-1)/2+1
        tem1(i,j,k)=v( i,ny-j+1,k) + v(i,j,k)
      END DO
    END DO
  END DO

  i = (nx-1)/2+1

  CALL wrigar(tem1,1,nx,1,ny,1,nz,i,i,1,(ny-1)/2+1,1,nz-1,ch,0.0,2)

  RETURN
END SUBROUTINE checkvy

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CHECKWY                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE checkwy(w, nx,ny,nz, ch, tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Check the south-north symmetry of a field defined at a w-point
!  on the x-z plane containing the maximum magnitude of that field.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue & Hao Jin
!  4/20/93.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    w        An array field defined at the w grid point.
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!  OUTPUT:
!
!    None.
!
!  WORK ARRAY:
!
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

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: w   (nx,ny,nz)       ! An array field defined at the w
                               ! grid point
  CHARACTER (LEN=*) :: ch

  REAL :: tem1(nx,ny,nz)
  INTEGER :: i,j,k
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  DO k=1,nz
    DO i=1,nx-1
      DO j=1,(ny-1)/2+1
        tem1(i,j,k)=w( i,ny-j,k) - w(i,j,k)
      END DO
    END DO
  END DO

  i = (nx-1)/2+1

  CALL wrigar(tem1,1,nx,1,ny,1,nz,i,i,1,(ny-1)/2,1,nz,ch,.0,2)

  RETURN
END SUBROUTINE checkwy
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CHECKSHX                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE checkshx(s, nx,ny,nz,                                        &
           ibgn,iend,jbgn,jend,kbgn,kend,                               &
           ch, tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Check the west-east symmetry of a field defined at a scalar point
!  on the level containing the maximum magnitude of that field.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue & Hao Jin
!  4/20/93.
!
!  MODIFICATION HISTORY:
!
!  02/12/1996 (Yuhe Liu)
!  Added the beginning and ending index for each direction
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    s        An array field defined at the scalar point.
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    ibgn     Beginning index in x-direction
!    iend     Ending    index in x-direction
!    jbgn     Beginning index in y-direction
!    jend     Ending    index in y-direction
!    kbgn     Beginning index in z-direction
!    kend     Ending    index in z-direction
!
!  OUTPUT:
!
!    None.
!
!  WORK ARRAY:
!
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

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  INTEGER :: ibgn,iend,jbgn,jend,kbgn,kend

  REAL :: s   (nx,ny,nz)       ! An array field defined at
                               ! the scalar point
  CHARACTER (LEN=*) :: ch

  REAL :: tem1(nx,ny,nz)
  INTEGER :: i,j,k
  REAL :: amin, amax
  INTEGER :: imax,jmax,kmax,imin,jmin,kmin
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO k=kbgn,kend
    DO j=jbgn,jend
      DO i=ibgn,(nx-1)/2+1
        tem1(i,j,k)=s( nx-i,j,k) - s(i,j,k)
      END DO
    END DO
  END DO

  CALL a3dmax(tem1,1,nx,ibgn,(nx-1)/2,1,ny,jbgn,jend,1,nz,kbgn,kend,    &
              amax,amin, imax,jmax,kmax, imin,jmin,kmin)

  WRITE(6,'(/2(1x,a,f13.7,3(a,i3)))')                                   &
       'sdmin =',amin,' at i=',imin,', j=',jmin,', k=',kmin,            &
       'sdmax =',amax,' at i=',imax,', j=',jmax,', k=',kmax

  IF ( ABS(amax) > ABS(amin) ) THEN
    k=kmax
  ELSE
    k=kmin
  END IF

  CALL wrigar(tem1,1,nx,1,ny,1,nz,ibgn,(nx-1)/2,jbgn,jend,k,k,          &
              ch,0.0,2)

  RETURN
END SUBROUTINE checkshx

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CHECKUHX                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE checkuhx(u, nx,ny,nz,                                        &
           ibgn,iend,jbgn,jend,kbgn,kend,                               &
           ch, tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Check the west-east symmetry of a field defined at a u-point
!  on the level containing the maximum magnitude of that field.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  11/20/91
!
!  MODIFICATION HISTORY:
!
!  5/06/92 (M. Xue)
!  Added full documentation.
!
!  02/12/1996 (Yuhe Liu)
!  Added the beginning and ending index for each direction
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    u        An array field defined at the u grid point.
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    ibgn     Beginning index in x-direction
!    iend     Ending    index in x-direction
!    jbgn     Beginning index in y-direction
!    jend     Ending    index in y-direction
!    kbgn     Beginning index in z-direction
!    kend     Ending    index in z-direction
!
!  OUTPUT:
!
!    None.
!
!  WORK ARRAY:
!
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

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  INTEGER :: ibgn,iend,jbgn,jend,kbgn,kend

  REAL :: u   (nx,ny,nz)       ! An array field defined at the u
                               ! grid point
  CHARACTER (LEN=*) :: ch

  REAL :: tem1(nx,ny,nz)
  INTEGER :: i,j,k
  REAL :: amin, amax
  INTEGER :: imax,jmax,kmax,imin,jmin,kmin
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  DO k=kbgn,kend
    DO j=jbgn,jend
      DO i=ibgn,(nx-1)/2+1
        tem1(i,j,k)=u( nx -i+1, j,k) + u(i, j,k)
      END DO
    END DO
  END DO

  CALL a3dmax(tem1,1,nx,ibgn,(nx-1)/2,1,ny,jbgn,jend,1,nz,kbgn,kend,    &
              amax,amin, imax,jmax,kmax, imin,jmin,kmin)

  WRITE(6,'(/2(1x,a,f13.7,3(a,i3)))')                                   &
       'udmin =',amin,' at i=',imin,', j=',jmin,', k=',kmin,            &
       'udmax =',amax,' at i=',imax,', j=',jmax,', k=',kmax

  IF ( ABS(amax) > ABS(amin) ) THEN
    k=kmax
  ELSE
    k=kmin
  END IF

  CALL wrigar(tem1,1,nx,1,ny,1,nz,ibgn,(nx-1)/2,jbgn,jend,k,k,          &
              ch,0.0,2)

  RETURN
END SUBROUTINE checkuhx

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CHECKVHX                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE checkvhx(v, nx,ny,nz,                                        &
           ibgn,iend,jbgn,jend,kbgn,kend,                               &
           ch, tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Check the west-east symmetry of a field defined at a v-point
!  on the level containing the maximum magnitude of that field.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  11/20/91.
!
!  MODIFICATION HISTORY:
!
!  5/06/92 (M. Xue)
!  Added full documentation.
!
!  02/12/1996 (Yuhe Liu)
!  Added the beginning and ending index for each direction
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    v        An array field defined at the v grid point.
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    ibgn     Beginning index in x-direction
!    iend     Ending    index in x-direction
!    jbgn     Beginning index in y-direction
!    jend     Ending    index in y-direction
!    kbgn     Beginning index in z-direction
!    kend     Ending    index in z-direction
!
!  OUTPUT:
!
!    None.
!
!  WORK ARRAY:
!
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

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  INTEGER :: ibgn,iend,jbgn,jend,kbgn,kend

  REAL :: v   (nx,ny,nz)
  CHARACTER (LEN=*) :: ch

  REAL :: tem1(nx,ny,nz)
  INTEGER :: i,j,k
  REAL :: amin, amax
  INTEGER :: imax,jmax,kmax,imin,jmin,kmin
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO k=kbgn,kend
    DO j=jbgn,jend
      DO i=ibgn,(nx-1)/2+1
        tem1(i,j,k)=v( nx-i,j,k) - v(i,j,k)
      END DO
    END DO
  END DO

  CALL a3dmax(tem1,1,nx,ibgn,(nx-1)/2,1,ny,jbgn,jend,1,nz,kbgn,kend,    &
              amax,amin, imax,jmax,kmax, imin,jmin,kmin)

  WRITE(6,'(/2(1x,a,f13.7,3(a,i3)))')                                   &
       'vdmin =',amin,' at i=',imin,', j=',jmin,', k=',kmin,            &
       'vdmax =',amax,' at i=',imax,', j=',jmax,', k=',kmax

  IF ( ABS(amax) > ABS(amin) ) THEN
    k=kmax
  ELSE
    k=kmin
  END IF

  CALL wrigar(tem1,1,nx,1,ny,1,nz,ibgn,(nx-1)/2,jbgn,jend,k,k,          &
              ch,0.0,2)

  RETURN
END SUBROUTINE checkvhx

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CHECKWHX                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE checkwhx(w, nx,ny,nz,                                        &
           ibgn,iend,jbgn,jend,kbgn,kend,                               &
           ch, tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Check the west-east symmetry of a field defined at a w-point
!  on the level containing the maximum magnitude of that field.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  11/20/91.
!
!  MODIFICATION HISTORY:
!
!  5/06/92 (M. Xue)
!  Added full documentation.
!
!  02/12/1996 (Yuhe Liu)
!  Added the beginning and ending index for each direction
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    w        An array field defined at the w grid point.
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    ibgn     Beginning index in x-direction
!    iend     Ending    index in x-direction
!    jbgn     Beginning index in y-direction
!    jend     Ending    index in y-direction
!    kbgn     Beginning index in z-direction
!    kend     Ending    index in z-direction
!
!  OUTPUT:
!
!    None.
!
!  WORK ARRAY:
!
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

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  INTEGER :: ibgn,iend,jbgn,jend,kbgn,kend

  REAL :: w   (nx,ny,nz)       ! An array field defined at the w
                               ! grid point
  CHARACTER (LEN=*) :: ch

  REAL :: tem1(nx,ny,nz)
  INTEGER :: i,j,k
  REAL :: amin, amax
  INTEGER :: imax,jmax,kmax,imin,jmin,kmin
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  DO k=kbgn,kend
    DO j=jbgn,jend
      DO i=ibgn,(nx-1)/2+1
        tem1(i,j,k)=w( nx-i,j,k) - w(i,j,k)
      END DO
    END DO
  END DO


  CALL a3dmax(tem1,1,nx,ibgn,(nx-1)/2,1,ny,jbgn,jend,1,nz,kbgn,kend,    &
              amax,amin, imax,jmax,kmax, imin,jmin,kmin)

  WRITE(6,'(/2(1x,a,f13.7,3(a,i3)))')                                   &
       'wdmin =',amin,' at i=',imin,', j=',jmin,', k=',kmin,            &
       'wdmax =',amax,' at i=',imax,', j=',jmax,', k=',kmax

  IF ( ABS(amax) > ABS(amin) ) THEN
    k=kmax
  ELSE
    k=kmin
  END IF

  CALL wrigar(tem1,1,nx,1,ny,1,nz,ibgn,(nx-1)/2,jbgn,jend,k,k,          &
              ch,0.0,2)

  RETURN
END SUBROUTINE checkwhx
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CHECKSHY                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE checkshy(s, nx,ny,nz,                                        &
           ibgn,iend,jbgn,jend,kbgn,kend,                               &
           ch, tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Check the south-north symmetry of a field defined at a scalar point
!  on the level containing the maximum magnitude of that field.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue & Hao Jin
!  4/20/93.
!
!  MODIFICATION HISTORY:
!
!  02/12/1996 (Yuhe Liu)
!  Added the beginning and ending index for each direction
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    s        An array field defined at the scalar point.
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    ibgn     Beginning index in x-direction
!    iend     Ending    index in x-direction
!    jbgn     Beginning index in y-direction
!    jend     Ending    index in y-direction
!    kbgn     Beginning index in z-direction
!    kend     Ending    index in z-direction
!
!  OUTPUT:
!
!    None.
!
!  WORK ARRAY:
!
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

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  INTEGER :: ibgn,iend,jbgn,jend,kbgn,kend

  REAL :: s   (nx,ny,nz)       ! An array field defined at the
                               ! scalar point
  CHARACTER (LEN=*) :: ch

  REAL :: tem1(nx,ny,nz)
  INTEGER :: i,j,k
  REAL :: amin, amax
  INTEGER :: imax,jmax,kmax,imin,jmin,kmin
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO k=kbgn,kend
    DO i=ibgn,iend
      DO j=jbgn,(ny-1)/2+1
        tem1(i,j,k)=s( i,ny-j,k) - s(i,j,k)
      END DO
    END DO
  END DO

  CALL a3dmax(tem1,1,nx,ibgn,iend,1,ny,jbgn,(ny-1)/2,1,nz,kbgn,kend,    &
              amax,amin, imax,jmax,kmax, imin,jmin,kmin)

  WRITE(6,'(/2(1x,a,f13.7,3(a,i3)))')                                   &
       'sdmin =',amin,' at i=',imin,', j=',jmin,', k=',kmin,            &
       'sdmax =',amax,' at i=',imax,', j=',jmax,', k=',kmax

  IF ( ABS(amax) > ABS(amin) ) THEN
    k=kmax
  ELSE
    k=kmin
  END IF

  CALL wrigar(tem1,1,nx,1,ny,1,nz,ibgn,iend,jbgn,(ny-1)/2,k,k,          &
              ch,0.0,2)

  RETURN
END SUBROUTINE checkshy

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CHECKUHY                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE checkuhy(u, nx,ny,nz,                                        &
           ibgn,iend,jbgn,jend,kbgn,kend,                               &
           ch, tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Check the south-north symmetry of a field defined at a u-point
!  on the level containing the maximum magnitude of that field.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue & Hao Jin
!  4/20/93.
!
!  MODIFICATION HISTORY:
!
!  02/12/1996 (Yuhe Liu)
!  Added the beginning and ending index for each direction
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    u        An array field defined at the u grid point.
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    ibgn     Beginning index in x-direction
!    iend     Ending    index in x-direction
!    jbgn     Beginning index in y-direction
!    jend     Ending    index in y-direction
!    kbgn     Beginning index in z-direction
!    kend     Ending    index in z-direction
!
!  OUTPUT:
!
!    None.
!
!  WORK ARRAY:
!
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

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  INTEGER :: ibgn,iend,jbgn,jend,kbgn,kend

  REAL :: u   (nx,ny,nz)       ! An array field defined at the u
                               ! grid point
  CHARACTER (LEN=*) :: ch

  REAL :: tem1(nx,ny,nz)
  INTEGER :: i,j,k
  REAL :: amin, amax
  INTEGER :: imax,jmax,kmax,imin,jmin,kmin
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  DO k=kbgn,kend
    DO i=ibgn,iend
      DO j=jbgn,(ny-1)/2+1
        tem1(i,j,k)=u( i,ny-j,k) - u(i, j,k)
      END DO
    END DO
  END DO

  CALL a3dmax(tem1,1,nx,ibgn,iend,1,ny,ibgn,(ny-1)/2,1,nz,kbgn,kend,    &
              amax,amin, imax,jmax,kmax, imin,jmin,kmin)

  WRITE(6,'(/2(1x,a,f13.7,3(a,i3)))')                                   &
       'udmin =',amin,' at i=',imin,', j=',jmin,', k=',kmin,            &
       'udmax =',amax,' at i=',imax,', j=',jmax,', k=',kmax

  IF ( ABS(amax) > ABS(amin) ) THEN
    k=kmax
  ELSE
    k=kmin
  END IF

  CALL wrigar(tem1,1,nx,1,ny,1,nz,ibgn,iend,ibgn,(ny-1)/2,k,k,          &
              ch,0.0,2)

  RETURN
END SUBROUTINE checkuhy

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CHECKVHY                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE checkvhy(v, nx,ny,nz,                                        &
           ibgn,iend,jbgn,jend,kbgn,kend,                               &
           ch, tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Check the south-north symmetry of a field defined at a v-point
!  on the level containing the maximum magnitude of that field.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue & Hao Jin
!  4/20/93.
!
!  MODIFICATION HISTORY:
!
!  02/12/1996 (Yuhe Liu)
!  Added the beginning and ending index for each direction
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    v        An array field defined at the v grid point.
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    ibgn     Beginning index in x-direction
!    iend     Ending    index in x-direction
!    jbgn     Beginning index in y-direction
!    jend     Ending    index in y-direction
!    kbgn     Beginning index in z-direction
!    kend     Ending    index in z-direction
!
!  OUTPUT:
!
!    None.
!
!  WORK ARRAY:
!
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

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  INTEGER :: ibgn,iend,jbgn,jend,kbgn,kend

  REAL :: v   (nx,ny,nz)       ! An array field defined at the v
                               ! grid point
  CHARACTER (LEN=*) :: ch

  REAL :: tem1(nx,ny,nz)
  INTEGER :: i,j,k
  REAL :: amin, amax
  INTEGER :: imax,jmax,kmax,imin,jmin,kmin
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO k=kbgn,kend
    DO i=ibgn,iend
      DO j=jbgn,(ny-1)/2+1
        tem1(i,j,k)=v( i,ny-j+1,k) + v(i,j,k)
      END DO
    END DO
  END DO

  CALL a3dmax(tem1,1,nx,ibgn,iend,1,ny,jbgn,(ny-1)/2,1,nz,kbgn,kend,    &
              amax,amin, imax,jmax,kmax, imin,jmin,kmin)

  WRITE(6,'(/2(1x,a,f13.7,3(a,i3)))')                                   &
       'vdmin =',amin,' at i=',imin,', j=',jmin,', k=',kmin,            &
       'vdmax =',amax,' at i=',imax,', j=',jmax,', k=',kmax

  IF ( ABS(amax) > ABS(amin) ) THEN
    k=kmax
  ELSE
    k=kmin
  END IF

  CALL wrigar(tem1,1,nx,1,ny,1,nz,ibgn,iend,jbgn,(ny-1)/2,k,k,          &
              ch,0.0,2)

  RETURN
END SUBROUTINE checkvhy

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CHECKWHY                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE checkwhy(w, nx,ny,nz,                                        &
           ibgn,iend,jbgn,jend,kbgn,kend,                               &
           ch, tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Check the south-north symmetry of a field defined at a w-point
!  on the level containing the maximum magnitude of that field.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue & Hao Jin
!  4/20/93.
!
!  MODIFICATION HISTORY:
!
!  02/12/1996 (Yuhe Liu)
!  Added the beginning and ending index for each direction
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    w        An array field defined at the w grid point.
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    ibgn     Beginning index in x-direction
!    iend     Ending    index in x-direction
!    jbgn     Beginning index in y-direction
!    jend     Ending    index in y-direction
!    kbgn     Beginning index in z-direction
!    kend     Ending    index in z-direction
!
!  OUTPUT:
!
!    None.
!
!  WORK ARRAY:
!
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

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  INTEGER :: ibgn,iend,jbgn,jend,kbgn,kend

  REAL :: w   (nx,ny,nz)       ! An array field defined at the w
                               ! grid point
  CHARACTER (LEN=*) :: ch

  REAL :: tem1(nx,ny,nz)
  INTEGER :: i,j,k
  REAL :: amin, amax
  INTEGER :: imax,jmax,kmax,imin,jmin,kmin
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  DO k=kbgn,kend
    DO i=ibgn,iend
      DO j=jbgn,(ny-1)/2+1
        tem1(i,j,k)=w( i,ny-j,k) - w(i,j,k)
      END DO
    END DO
  END DO

  CALL a3dmax(tem1,1,nx,ibgn,iend,jbgn,ny,1,(ny-1)/2,1,nz,kbgn,kend,    &
              amax,amin, imax,jmax,kmax, imin,jmin,kmin)

  WRITE(6,'(/2(1x,a,f13.7,3(a,i3)))')                                   &
       'wdmin =',amin,' at i=',imin,', j=',jmin,', k=',kmin,            &
       'wdmax =',amax,' at i=',imax,', j=',jmax,', k=',kmax

  IF ( ABS(amax) > ABS(amin) ) THEN
    k=kmax
  ELSE
    k=kmin
  END IF

  CALL wrigar(tem1,1,nx,1,ny,1,nz,ibgn,iend,jbgn,(ny-1)/2,k,k,          &
              ch,0.0,2)

  RETURN
END SUBROUTINE checkwhy
