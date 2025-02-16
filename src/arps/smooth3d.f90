!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SMOOTH9P_NOBC              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE smooth9p_nobc( arr, nx,ny,ibgn,iend,jbgn,jend, tem1 )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!                                        1 2 1
!  Smooth a 2-D array by the filter of { 2 4 2 }. This subroutine
!                                        1 2 1
!  does not deal with the boundaries. Users are supposed to deal with
!  the boundaries themselves based on how they want the BC to be.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:       Yuhe Liu
!
!  4/13/99
!
!  Modification History
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  nx       Number of grid points in the x-direction
!  ny       Number of grid points in the y-direction
!  ibgn     First index in x-direction in the soomthing region.
!  iend     Last  index in x-direction in the soomthing region.
!  jbgn     First index in j-direction in the soomthing region.
!  jend     Last  index in j-direction in the soomthing region.
!
!  arr    2-D array
!
!  OUTPUT:
!
!  arr    2-D array
!
!  TEMPORARY:
!
!  tem1     Temporary 2-D array
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: nx         ! Number of grid points in the x-direction
  INTEGER :: ny         ! Number of grid points in the y-direction
  INTEGER :: ibgn
  INTEGER :: iend
  INTEGER :: jbgn
  INTEGER :: jend
!
  REAL :: arr (nx,ny)   ! 2-D array
!
  REAL :: tem1(nx,ny)   ! Temporary array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j
  REAL :: wtf,wtfb,wtfc
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  wtf    = 1.0/16.0
  wtfb = 1.0/12.0
  wtfc = 1.0/9.0

  DO j = jbgn+1, jend-1
    DO i = ibgn+1, iend-1

      tem1(i,j) = wtf                                                   &
          * (    arr(i-1,j-1) + 2.*arr(i,j-1) +    arr(i+1,j-1)         &
          + 2.*arr(i-1,j  ) + 4.*arr(i,j  ) + 2.*arr(i+1,j  )           &
          +    arr(i-1,j+1) + 2.*arr(i,j+1) +    arr(i+1,j+1) )

    END DO
  END DO

  DO j = 1, ny
    DO i = 1, nx
      arr(i,j) = tem1(i,j)
    END DO
  END DO

  RETURN
END SUBROUTINE smooth9p_nobc
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SMOOTH9P                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE smooth9p( arr, nx,ny,ibgn,iend,jbgn,jend, stagdim, tem1 )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!                                        1 2 1
!  Smooth a 2-D array by the filter of { 2 4 2 }
!                                        1 2 1
!
!-----------------------------------------------------------------------
!
!  AUTHOR:       Yuhe Liu
!
!  5/3/94
!
!  Modification History
!  8/20/1995 (M. Xue)
!  Fixed errors in the index bound of loops 100 and 200.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  nx       Number of grid points in the x-direction
!  ny       Number of grid points in the y-direction
!  ibgn     First index in x-direction in the soomthing region.
!  iend     Last  index in x-direction in the soomthing region.
!  jbgn     First index in j-direction in the soomthing region.
!  jend     Last  index in j-direction in the soomthing region.
!
!  arr    2-D array
!
!  OUTPUT:
!
!  arr    2-D array
!
!  TEMPORARY:
!
!  tem1     Temporary 2-D array
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: nx         ! Number of grid points in the x-direction
  INTEGER :: ny         ! Number of grid points in the y-direction
  INTEGER :: ibgn
  INTEGER :: iend
  INTEGER :: jbgn
  INTEGER :: jend
  INTEGER, INTENT(IN) :: stagdim

  REAL :: arr (nx,ny)   ! 2-D array
!
  REAL :: tem1(nx,ny)   ! Temporary array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j
  REAL :: wtf,wtfb,wtfc
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'mp.inc'
  INCLUDE 'bndry.inc'
  INCLUDE 'globcst.inc'   ! just for mp_acct
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  wtf  = 1.0/16.0
  wtfb = 1.0/12.0
  wtfc = 1.0/9.0

  DO j = jbgn+1, jend-1
    DO i = ibgn+1, iend-1

      tem1(i,j) = wtf                                                   &
          * (    arr(i-1,j-1) + 2.*arr(i,j-1) +    arr(i+1,j-1)         &
            + 2.*arr(i-1,j  ) + 4.*arr(i,j  ) + 2.*arr(i+1,j  )         &
            +    arr(i-1,j+1) + 2.*arr(i,j+1) +    arr(i+1,j+1) )

    END DO
  END DO

  DO j = jbgn+1, jend-1
    tem1(ibgn,j) = wtfb                                                 &
        * ( 2.*arr(ibgn,j-1) +    arr(ibgn+1,j-1)                       &
        + 4.*arr(ibgn,j  ) + 2.*arr(ibgn+1,j  )                         &
        + 2.*arr(ibgn,j+1) +    arr(ibgn+1,j+1) )

    tem1(iend,j) = wtfb                                                 &
        * (    arr(iend-1,j-1)  + 2.*arr(iend,j-1)                      &
        + 2.*arr(iend-1,j  )  + 4.*arr(iend,j  )                        &
        +    arr(iend-1,j+1)  + 2.*arr(iend,j+1) )

  END DO

  DO i = ibgn+1, iend-1
    tem1(i,jbgn) = wtfb                                                 &
        * ( 2.*arr(i-1,jbgn  ) + 4.*arr(i,jbgn  ) + 2.*arr(i+1,jbgn  )  &
        +    arr(i-1,jbgn+1) + 2.*arr(i,jbgn+1) +    arr(i+1,jbgn+1) )

    tem1(i,jend) = wtfb                                                 &
        * (    arr(i-1,jend-1) + 2.*arr(i,jend-1) +    arr(i+1,jend-1)  &
        + 2.*arr(i-1,jend  ) + 4.*arr(i,jend  ) + 2.*arr(i+1,jend  ) )

  END DO

  tem1(ibgn,jbgn) = wtfc                                                &
      * ( 2.*arr(ibgn,jbgn+1) +    arr(ibgn+1,jbgn+1)                   &
      + 4.*arr(ibgn,jbgn  ) + 2.*arr(ibgn+1,jbgn  ) )

  tem1(ibgn,jend) = wtfc                                                &
      * ( 4.*arr(ibgn,jend  ) + 2.*arr(ibgn+1,jend  )                   &
      + 2.*arr(ibgn,jend-1) +    arr(ibgn+1,jend-1) )

  tem1(iend,jbgn) = wtfc                                                &
      * (    arr(iend-1,jbgn+1) + 2.*arr(iend,jbgn+1)                   &
      + 2.*arr(iend-1,jbgn  ) + 4.*arr(iend,jbgn  ) )

  tem1(iend,jend) = wtfc                                                &
      * ( 2.*arr(iend-1,jend  ) + 4.*arr(iend,  jend)                   &
       +    arr(iend-1,jend-1) + 2.*arr(iend-1,jend) )

  DO j = 1, ny
    DO i = 1, nx
      arr(i,j) = tem1(i,j)
    END DO
  END DO

!
! Ensure the out array is mpi valid, the input MUST be already MPI valid.
!

  IF (mp_opt > 0) THEN    
    CALL acct_interrupt(mp_acct)
    CALL mpsendrecv1dew(arr,nx,ny,ebc,wbc,stagdim,tem1)
    CALL mpsendrecv1dns(arr,nx,ny,nbc,sbc,stagdim,tem1)
  END IF

  RETURN
END SUBROUTINE smooth9p
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SMOOTH25P                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE smooth25p( arr, nx,ny,ibgn,iend,jbgn,jend, tem1 )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  To perform a 25 point smoothing on a 2-D array.
!
!
!                                        1   4   6  4   1
!                                        4  16  24 16   4
!  Smooth a 2-D array by the filter of { 6  24  36 24   6 }
!                                        4  16  24 16   4
!                                        1   4   6  4   1
!
!
!                                                       1 2 1
!  At the inner boundary 9-point filter will be used: { 2 4 2 }
!                                                       1 2 1
!
!                                                       2 1
!  At the outer boundary 6-point filter will be used: { 4 2 }
!                                                       2 1
!
!                                           2 1
!  At the corners 4-point filter is used: {     }
!                                           4 2
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  5/3/94
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  nx       Number of grid points in the x-direction
!  ny       Number of grid points in the y-direction
!  ibgn     First index in x-direction in the soomthing region.
!  iend     Last  index in x-direction in the soomthing region.
!  jbgn     First index in j-direction in the soomthing region.
!  jend     Last  index in j-direction in the soomthing region.
!
!  arr    2-D array
!
!  OUTPUT:
!
!  arr    2-D array smoothed
!
!  TEMPORARY:
!
!  tem1     Temporary 2-D array
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: nx         ! Number of grid points in the x-direction
  INTEGER :: ny         ! Number of grid points in the y-direction

  INTEGER :: ibgn,iend,jbgn,jend
!
  REAL :: arr(nx,ny)  ! 2-D array
!
  REAL :: tem1(nx,ny)  ! Temporary array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j
  REAL :: wtf
  REAL :: wtfb1
  REAL :: wtfb2
  REAL :: wtfc

  REAL :: wt(5,5)
  DATA wt/ 1.,  4.,  6.,  4.,  1.,                                      &
           4., 16., 24., 16.,  4.,                                      &
           6., 24., 36., 24.,  6.,                                      &
           4., 16., 24., 16.,  4.,                                      &
           1.,  4.,  6.,  4.,  1./

  REAL :: wtb(3,3)
  DATA wtb/ 1., 2., 1.,                                                 &
            2., 4., 2.,                                                 &
            1., 2., 1./
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  wtf   = 1./256.

  DO j = jbgn+2, jend-2
    DO i = ibgn+2, iend-2

      tem1(i,j) = wtf                                                   &
          * ( wt(1,1)*arr(i-2,j-2) + wt(2,1)*arr(i-1,j-2)               &
          + wt(3,1)*arr(i,  j-2) + wt(4,1)*arr(i+1,j-2)                 &
          + wt(5,1)*arr(i+2,j-2)                                        &
          + wt(1,2)*arr(i-2,j-1) + wt(2,2)*arr(i-1,j-1)                 &
          + wt(3,2)*arr(i,  j-1) + wt(4,2)*arr(i+1,j-1)                 &
          + wt(5,2)*arr(i+2,j-1)                                        &
          + wt(1,3)*arr(i-2,j  ) + wt(2,3)*arr(i-1,j  )                 &
          + wt(3,3)*arr(i,  j  ) + wt(4,3)*arr(i+1,j  )                 &
          + wt(5,3)*arr(i+2,j  )                                        &
          + wt(1,4)*arr(i-2,j+1) + wt(2,4)*arr(i-1,j+1)                 &
          + wt(3,4)*arr(i,  j+1) + wt(4,4)*arr(i+1,j+1)                 &
          + wt(5,4)*arr(i+2,j+1)                                        &
          + wt(1,5)*arr(i-2,j+2) + wt(2,5)*arr(i-1,j+2)                 &
          + wt(3,5)*arr(i,  j+2) + wt(4,5)*arr(i+1,j+2)                 &
          + wt(5,5)*arr(i+2,j+2) )

    END DO
  END DO

  wtfb1 = 1./16.
  wtfb2 = 1./12.

  DO j = jbgn+1, jend-1
    tem1(ibgn+1,j) = wtfb1                                              &
        * ( wtb(1,1)*arr(ibgn,  j-1) + wtb(2,1)*arr(ibgn+1,j-1)         &
        + wtb(3,1)*arr(ibgn+2,j-1)                                      &
        + wtb(1,2)*arr(ibgn,  j  ) + wtb(2,2)*arr(ibgn+1,j  )           &
        + wtb(3,2)*arr(ibgn+2,j  )                                      &
        + wtb(1,3)*arr(ibgn  ,j+1) + wtb(2,3)*arr(ibgn+1,j+1)           &
        + wtb(3,3)*arr(ibgn+2,j+1) )

    tem1(ibgn,j) = wtfb2                                                &
        * ( wtb(2,1)*arr(ibgn,j-1) + wtb(3,1)*arr(ibgn+1,j-1)           &
        + wtb(2,2)*arr(ibgn,j  ) + wtb(3,2)*arr(ibgn+1,j  )             &
        + wtb(2,3)*arr(ibgn,j+1) + wtb(3,3)*arr(ibgn+1,j+1) )

    tem1(iend-1,j) = wtfb1                                              &
        * ( wtb(1,1)*arr(iend-2,j-1) + wtb(2,1)*arr(iend-1,j-1)         &
        + wtb(3,1)*arr(iend,  j-1)                                      &
        + wtb(1,2)*arr(iend-2,j  ) + wtb(2,2)*arr(iend-1,j  )           &
        + wtb(3,2)*arr(iend,  j  )                                      &
        + wtb(1,3)*arr(iend-2,j+1) + wtb(2,3)*arr(iend-1,j+1)           &
        + wtb(3,3)*arr(iend,  j+1) )

    tem1(iend,j) = wtfb2                                                &
        * ( wtb(1,1)*arr(iend-1,j-1) + wtb(2,1)*arr(iend,j-1)           &
        + wtb(1,2)*arr(iend-1,j  ) + wtb(2,2)*arr(iend,j  )             &
        + wtb(1,3)*arr(iend-1,j+1) + wtb(2,3)*arr(iend,j+1) )

  END DO

  DO i = ibgn+1, iend-1
    tem1(i,jbgn+1) = wtfb1                                              &
        * ( wtb(1,1)*arr(i-1,jbgn)   + wtb(2,1)*arr(i,jbgn  )           &
        + wtb(3,1)*arr(i+1,jbgn)                                        &
        + wtb(1,2)*arr(i-1,jbgn+1) + wtb(2,2)*arr(i,jbgn+1)             &
        + wtb(3,2)*arr(i+1,jbgn+1)                                      &
        + wtb(1,3)*arr(i-1,jbgn+2) + wtb(2,3)*arr(i,jbgn+2)             &
        + wtb(3,3)*arr(i+1,jbgn+2) )

    tem1(i,jbgn) = wtfb2                                                &
        * ( wtb(1,2)*arr(i-1,jbgn)   + wtb(2,2)*arr(i,jbgn  )           &
        + wtb(3,2)*arr(i+1,jbgn)                                        &
        + wtb(1,3)*arr(i-1,jbgn+1) + wtb(2,3)*arr(i,jbgn+1)             &
        + wtb(3,3)*arr(i+1,jbgn+1) )

    tem1(i,jend-1) = wtfb1                                              &
        * ( wtb(1,1)*arr(i-1,jend-2) + wtb(2,1)*arr(i,jend-2)           &
        + wtb(3,1)*arr(i+1,jend-2)                                      &
        + wtb(1,2)*arr(i-1,jend-1) + wtb(2,2)*arr(i,jend-1)             &
        + wtb(3,2)*arr(i+1,jend-1)                                      &
        + wtb(1,3)*arr(i-1,jend  ) + wtb(2,3)*arr(i,jend  )             &
        + wtb(3,3)*arr(i+1,jend  ) )

    tem1(i,jend) = wtfb2                                                &
        * ( wtb(1,1)*arr(i-1,jend-1) + wtb(2,1)*arr(i,jend-1)           &
        + wtb(3,1)*arr(i+1,jend-1)                                      &
        + wtb(1,2)*arr(i-1,jend  ) + wtb(2,2)*arr(i,jend  )             &
        + wtb(3,2)*arr(i+1,jend  ) )

  END DO

  wtfc  = 1./9.

  tem1(ibgn,jbgn) = wtfc                                                &
      * ( wtb(1,2)*arr(ibgn,jbgn+1) + wtb(1,3)*arr(ibgn+1,jbgn+1)       &
      + wtb(2,2)*arr(ibgn,jbgn  ) + wtb(2,3)*arr(ibgn+1,jbgn  ) )


  tem1(iend,jbgn) = wtfc                                                &
      * ( wtb(1,1)*arr(iend-1,jbgn+1) + wtb(1,2)*arr(iend,jbgn+1)       &
      + wtb(2,1)*arr(iend-1,jbgn  ) + wtb(2,2)*arr(iend,jbgn  ) )

  tem1(ibgn,jend) = wtfc                                                &
      * ( wtb(2,2)*arr(ibgn,jend  ) + wtb(3,2)*arr(ibgn+1,jend  )       &
      + wtb(2,3)*arr(ibgn,jend-1) + wtb(3,3)*arr(ibgn+1,jend-1) )

  tem1(iend,jend) = wtfc                                                &
      * ( wtb(2,1)*arr(iend-1,jend  ) + wtb(2,2)*arr(iend,jend  )       &
      + wtb(3,1)*arr(iend-1,jend-1) + wtb(3,2)*arr(iend,jend-1) )

  DO j = jbgn, jend
    DO i = ibgn, iend
      arr(i,j) = tem1(i,j)
    END DO
  END DO

  RETURN
END SUBROUTINE smooth25p
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SMOOTH3D                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE smooth3d(nx,ny,nz,ibgn,iend,jbgn,jend,kbgn,kend,stagdim,     &
                    s,ht,zin,zwork,zout)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Performs symmetrical 3-dimensional smoothing of input field
!  zin which is output as zout, the smoothed field.  A work array
!  zwork, dimension (nx,ny) is required.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:  K. Brewster
!  3/1989
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  nx       Number of grid points in the x-direction
!  ny       Number of grid points in the y-direction
!  nz       Number of grid points in the z-direction
!  ibgn     First index in i-direction in the soomthing region.
!  iend     Last  index in i-direction in the soomthing region.
!  jbgn     First index in j-direction in the soomthing region.
!  jend     Last  index in j-direction in the soomthing region.
!  kbgn     First index in k-direction in the soomthing region.
!  kend     Last  index in k-direction in the soomthing region.
!
!  s
!  ht
!
!  zin      3-d array to be smoothed.
!
!  OUTPUT:
!
!  zout     3-D smoothed array
!
!  WORK ARRAY:
!
!  zwork    Temporary work array
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz
  INTEGER :: ibgn,iend
  INTEGER :: jbgn,jend
  INTEGER :: kbgn,kend
  INTEGER :: stagdim

  REAL :: s
  REAL :: ht(nx,ny,nz)
  REAL :: zin(nx,ny,nz),zwork(nx,ny,nz),zout(nx,ny,nz)

!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------

  INTEGER :: i,j,k
  REAL    :: wcen,wsid,whigh

  INCLUDE 'mp.inc'
  INCLUDE 'bndry.inc'     ! used wbc, ebc, nbc etc.
  INCLUDE 'globcst.inc'   ! just for mp_acct

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code ... ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!
!  write (6,'(a/a,7x,a,12x,a,12x,a,12x,a,12x,a)')
!    :'zin(i,j,k), k=5', ' j\\i', '39', '40', '41', '42', '43'
!  write (6,'(i3,5e14.6)') (j,(zin(i,j,5),i=nx-4,nx),j=ny-4,ny)

  wcen=1.-s
  wsid=s*0.5
!
!-----------------------------------------------------------------------
!
!  I direction 3pt smoothing
!
!-----------------------------------------------------------------------
!
  DO k=kbgn, kend
    DO j=jbgn, jend
      DO i=ibgn+1, iend-1
        zwork(i,j,k)=zin(i  ,j,k)*wcen +                                &
                     zin(i+1,j,k)*wsid +                                &
                     zin(i-1,j,k)*wsid
      END DO
    END DO
    DO j=jbgn,jend
      zwork(ibgn,j,k)=zin(ibgn,j,k)
      zwork(iend,j,k)=zin(iend,j,k)
    END DO
!
!-----------------------------------------------------------------------
!
!  J direction 3 point smoothing.
!
!-----------------------------------------------------------------------
!
    DO j=jbgn+1, jend-1
      DO i=ibgn, iend
        zout(i,j,k)=zwork(i  ,j,k)*wcen +                               &
                    zwork(i,j+1,k)*wsid +                               &
                    zwork(i,j-1,k)*wsid
      END DO
    END DO
    DO i=ibgn, iend
      zout(i,jbgn,k)=zwork(i,jbgn,k)
      zout(i,jend,k)=zwork(i,jend,k)
    END DO
  END DO
!
!  Transfer output back to work array.
!
  DO k=kbgn,kend
    DO j=jbgn,jend
      DO i=ibgn,iend
        zwork(i,j,k)=zout(i,j,k)
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Smooth in vertical
!
!-----------------------------------------------------------------------
!
  DO j=jbgn,jend
    DO i=ibgn,iend
      DO k=kbgn+1,kend-1
        whigh=(ht(i,j,  k)-ht(i,j,k-1))/                                &
              (ht(i,j,k+1)-ht(i,j,k-1))
        zout(i,j,k)=zwork(i,j,k)*wcen +                                 &
                    s*( whigh*zwork(i,j,k+1) +                          &
                      (1.-whigh)*zwork(i,j,k-1) )
      END DO
      zout(i,j,kbgn)=zwork(i,j,kbgn)
      zout(i,j,kend)=zwork(i,j,kend)
    END DO
  END DO

!  write (6,'(a/a,7x,a,12x,a,12x,a,12x,a,12x,a)')                 &
!    'zout(i,j,k), k=5', ' j\\i', '65', '66', '67', '68', '69'
!  write (6,'(i3,5e14.6)') (j,(zout(i,j,5),i=nx-4,nx),j=ny-4,ny)
!
!  write (6,'(a/a,7x,a,12x,a,12x,a,12x,a,12x,a)')                 &
!    'zwork(i,j,k), k=5', ' j\\i', '65', '66', '67', '68', '69'
!  write (6,'(i3,5e14.6)') (j,(zwork(i,j,5),i=nx-4,nx),j=ny-4,ny)

!
! Ensure the out  array is mpi valid
!
  IF (mp_opt > 0) THEN    ! Ensure output is good for followed smoothing
    CALL acct_interrupt(mp_acct)
    CALL mpsendrecv2dew(zout,nx,ny,nz,ebc,wbc,stagdim,zwork)
    CALL mpsendrecv2dns(zout,nx,ny,nz,nbc,sbc,stagdim,zwork)
  END IF

  RETURN
END SUBROUTINE smooth3d
