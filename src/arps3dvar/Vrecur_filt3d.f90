
SUBROUTINE recurfilt3d(nx,ny,nz,pgrd,ipass_filt,ipass_loop,             &
                       hradius,radius_z,alpha,beta,tem1,tem2)
!-----------------------------------------------------------------------
!
! PURPOSE:
!    Parallel version of the 3D recursive filter.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!
!  Yunheng Wang, CAPS, OU, 11/01/2007
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx,ny,nz
  INTEGER, INTENT(IN) :: ipass_filt, ipass_loop
  REAL,    INTENT(IN) :: hradius
  REAL,    INTENT(IN) :: radius_z(nx,ny,nz)


  REAL,    INTENT(INOUT) :: pgrd(nx,ny,nz)

  ! working arrays
  REAL,    INTENT(INOUT) :: alpha(nx,ny,nz)
  REAL,    INTENT(INOUT) :: beta (nx,ny,nz)

  REAL,    INTENT(INOUT) :: tem1(nx,ny,nz)   ! Temporary array for MPI
                                             ! Can be two dimensions as (max(nx,ny),nz)
  REAL,    INTENT(INOUT) :: tem2(nx,ny,nz)   ! Temporary array for MPI
                                             ! Can be two dimensions as (max(nx,ny),nz)

!-----------------------------------------------------------------------
!
! Misc local variables
!
!-----------------------------------------------------------------------

  REAL :: ee
  REAL :: temp1,temp2

  INTEGER :: i,j,k,n

  INCLUDE 'mp.inc'

  INTEGER :: ic, jc

!  integer :: nxlg, nylg
!  real, allocatable :: outlg(:,:,:)
!  logical :: outp

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF( hradius == 0 ) return

!  outp = .false.
!  if (ipass_filt >= 100) then
!    ipass_filt = ipass_filt-100
!    outp = .true.
!
!    nxlg = (nx-3)*nproc_x+3
!    nylg = (ny-3)*nproc_y+3
!    Allocate(outlg(nxlg,nylg,nz),stat=ic)
!    CALL mpimerge3d(pgrd,nx,ny,nz,outlg)
!
!  if (mp_opt >0 .and. myproc == 0) ThEN
!    DO k = 2,2
!     do j = 1,1
!     do i = 1,nxlg
!       write(100,*) i, outlg(i,j,k)
!     end do
!       write(100,*)
!    end do
!    end do
!  else if (mp_opt == 0) then
!    DO k = 2,2
!     do j = 1,1
!     do i = 1,nxlg
!    write(100,*) i, outlg(i,j,k)
!     end do
!    write(100,*)
!    end do
!    end do
!  end if
!  end if
!
!write(0,*) ipass_filt

  DO n = 1, ipass_loop ! ipass_filt/2

!write(0,*) 'ipass = ',n, loc_x, loc_y
!-----------------------------------------------------------------------
!
! X direction - forward
!
!-----------------------------------------------------------------------

    DO k = 1, nz-1
      DO j = 1, ny-1
        DO i = 1, nx-1
          ee = REAL(ipass_filt) / (hradius*hradius)
          alpha(i,j,k) = 1+ee-SQRT( ee*(ee+2.) )
          beta(i,j,k) = 1.-alpha(i,j,k)
        END DO
      END DO
    END DO

    CALL inctag

    DO ic = 1, nproc_x

      IF (loc_x == ic) THEN
        IF (loc_x == 1) THEN
          CALL set_bdyxs(pgrd,nx,ny,nz,n,alpha,beta)
        ELSE
          CALL receive_bdyxs(pgrd,nx,ny,nz,tem1)
        END IF

        DO k = 1,nz-1
          DO j = 1, ny-1
            DO i = 2, nx-1
             pgrd(i,j,k)=alpha(i,j,k)*pgrd(i-1,j,k)+beta(i,j,k)*pgrd(i,j,k)
            END DO
          END DO
        END DO

        IF (loc_x < nproc_x) THEN
          CALL send_next_bdyxs(pgrd,nx,ny,nz,tem1)
        END IF
      END IF

      CALL mpbarrier
    END DO

!-----------------------------------------------------------------------
!
! X direction - backward
!
!-----------------------------------------------------------------------

    CALL inctag

    DO ic = nproc_x, 1, -1

      IF (loc_x == ic) THEN
        IF (loc_x == nproc_x) THEN
          CALL set_bdyxe(pgrd,nx,ny,nz,n,alpha,beta)
        ELSE
          CALL receive_bdyxe(pgrd,nx,ny,nz,tem1)
        END IF

        DO k = 1,nz-1
          DO j = 1,ny-1
            DO i = nx-2, 1, -1
             pgrd(i,j,k)=alpha(i,j,k)*pgrd(i+1,j,k)+beta(i,j,k)*pgrd(i,j,k)
            END DO
          END DO
        END DO

        IF (loc_x > 1) THEN
          CALL send_previous_bdyxe(pgrd,nx,ny,nz,tem1)
        END IF
      END IF

!write(0,*) '=== after backward:', nx/2,ny/2,nz/2,pgrd(nx/2,ny/2,nz/2)
      CALL mpbarrier
    END DO

!-----------------------------------------------------------------------
!
! Y direction - forward
!
!-----------------------------------------------------------------------

    CALL inctag

    DO jc = 1, nproc_y

      IF (jc == loc_y) THEN
        IF (loc_y == 1) THEN
          CALL set_bdyys(pgrd,nx,ny,nz,n,alpha,beta)
        ELSE
          CALL receive_bdyys(pgrd,nx,ny,nz,tem1)
        END IF

        DO k = 1, nz-1
          DO j = 2, ny-1
            DO i = 1, nx-1
              pgrd(i,j,k)=alpha(i,j,k)*pgrd(i,j-1,k)+beta(i,j,k)*pgrd(i,j,k)
            END DO
          END DO
        END DO

        IF (loc_y < nproc_y) THEN
           CALL send_up_bdyys(pgrd,nx,ny,nz,tem1)
        END IF

      END IF

      CALL mpbarrier
    END DO

!write(0,*) '*--- after forward:', n,' --- ',nx/2,ny/2,nz/2,pgrd(nx/2,ny/2,nz/2)

!-----------------------------------------------------------------------
!
! Y direction - backword
!
!-----------------------------------------------------------------------

    CALL inctag

    DO jc = nproc_y, 1, -1

      IF (jc == loc_y) THEN

        IF (loc_y == nproc_y) THEN
          CALL set_bdyye(pgrd,nx,ny,nz,n,alpha,beta)
        ELSE
          CALL receive_bdyye(pgrd,nx,ny,nz,tem1)
        END IF

        DO k = 1,nz-1
          DO j = ny-2,1,-1
            DO i = 1,nx-1
              pgrd(i,j,k)=alpha(i,j,k)*pgrd(i,j+1,k)+beta(i,j,k)*pgrd(i,j,k)
            END DO
          END DO
        END DO

        IF (loc_y > 1) THEN
          CALL send_down_bdyye(pgrd,nx,ny,nz,tem1)
        END IF
      END IF

      CALL mpbarrier
    END DO

!write(0,*) '*=== after backward:',n,' --- ', nx/2,ny/2,nz/2,pgrd(nx/2,ny/2,nz/2)
!-----------------------------------------------------------------------
!
! Z direction - forward
!
!-----------------------------------------------------------------------

!   IF( radius_z /= 0 ) THEN

      DO k = 1, nz-1
        DO j = 1, ny-1
          DO i = 1, nx-1
            ee = REAL(ipass_filt)/(radius_z(i,j,k)*radius_z(i,j,k))
            alpha(i,j,k) = 1+ee-SQRT( ee*(ee+2.) )
             beta(i,j,k) = 1.-alpha(i,j,k)
          END DO
        END DO
      END DO

      SELECT CASE (n)
      CASE (1)
        DO j = 1, ny-1
          DO i = 1, nx-1
            pgrd(i,j,1) = beta(i,j,1) * pgrd(i,j,1)
          END DO
        END DO
      CASE (2)
        DO j = 1, ny-1
          DO i = 1, nx-1
            temp1 = 1.-alpha(i,j,1)*alpha(i,j,1)
            pgrd(i,j,1) = beta(i,j,1)/temp1 * pgrd(i,j,1)
          END DO
        END DO
      CASE (3)
        DO j = 1, ny-1
          DO i = 1, nx-1
            temp1 = (1-alpha(i,j,1))/                                            &
                ((1-alpha(i,j,1)*alpha(i,j,1))*(1-alpha(i,j,1)*alpha(i,j,1)))
            temp2 =alpha(i,j,1)*alpha(i,j,1)*alpha(i,j,1)
            pgrd(i,j,1) = temp1 * (pgrd(i,j,1)-temp2*pgrd(i,j,2))
          END DO
        END DO
      CASE DEFAULT
        DO j = 1, ny-1
          DO i = 1, nx-1
            temp2 =alpha(i,j,1)*alpha(i,j,1)*alpha(i,j,1)
            temp1 = (1-alpha(i,j,1))/                                            &
               (1-3*alpha(i,j,1)*alpha(i,j,1)+3*temp2*alpha(i,j,1)-temp2*temp2)
            pgrd(i,j,1) = temp1 * (pgrd(i,j,1)-3*temp2*pgrd(i,j,2)+              &
               temp2*alpha(i,j,1)*alpha(i,j,1)*pgrd(i,j,2)                       &
               +temp2*alpha(i,j,1)*pgrd(i,j,3))
          END DO
        END DO
      END SELECT

      DO k = 2, nz, 1
        DO j = 1, ny-1
          DO i = 1, nx-1
           pgrd(i,j,k) = alpha(i,j,k)*pgrd(i,j,k-1)+beta(i,j,k)*pgrd(i,j,k)
          END DO
        END DO
      END DO

!write(0,*) '##--- after forward:', nx/2,ny/2,nz/2,pgrd(nx/2,ny/2,nz/2)

!-----------------------------------------------------------------------
!
! Z direction - backward
!
!-----------------------------------------------------------------------

      SELECT CASE (n)
      CASE (0)
        DO j = 1, ny-1
          DO i = 1, nx-1
            pgrd(i,j,nz-1) = beta(i,j,nz-1) * pgrd(i,j,nz-1)
          END DO
        END DO
      CASE (1)
        DO j = 1, ny-1
          DO i = 1, nx-1
            temp1 = (1.-alpha(i,j,nz-1)*alpha(i,j,nz-1))
            pgrd(i,j,nz-1) = beta(i,j,nz-1)/temp1 * pgrd(i,j,nz-1)
          END DO
        END DO
      CASE (2)
        DO j = 1, ny-1
          DO i = 1, nx-1
            temp1 = (1-alpha(i,j,nz-1))/                                               &
             ((1-alpha(i,j,nz-1)*alpha(i,j,nz-1))*(1-alpha(i,j,nz-1)*alpha(i,j,nz-1)))
            temp2 = alpha(i,j,nz-1)*alpha(i,j,nz-1)*alpha(i,j,nz-1)
            pgrd(i,j,nz-1) = temp1 * (pgrd(i,j,nz-1)-temp2*pgrd(i,j,nz-2))
          END DO
        END DO
      CASE DEFAULT
        DO j = 1, ny-1
          DO i = 1, nx-1
            temp2 = alpha(i,j,nz-1)*alpha(i,j,nz-1)*alpha(i,j,nz-1)
            temp1 = (1-alpha(i,j,nz-1))/(1-3*alpha(i,j,nz-1)*alpha(i,j,nz-1)+          &
                              3*temp2*alpha(i,j,nz-1)-temp2*temp2)
            pgrd(i,j,nz-1) = temp1 * (pgrd(i,j,nz-1)-3*temp2*pgrd(i,j,nz-2)+           &
                   temp2*alpha(i,j,nz-1)*alpha(i,j,nz-1)*pgrd(i,j,nz-2)+               &
                   temp2*alpha(i,j,nz-1)*pgrd(i,j,nz-3))
          END DO
        END DO
      END SELECT

      DO k = nz-2, 1, -1
        DO j = 1, ny-1
          DO i = 1, nx-1
            pgrd(i,j,k)=alpha(i,j,k)*pgrd(i,j,k+1)+beta(i,j,k)*pgrd(i,j,k)
          END DO
        END DO
      END DO

!   ENDIF

  END DO

!
!if (myproc == 0) then
!write(0,*) outlg(2,2,30)
!endif
!
!write(0,*) outp
!  outp = .false.
!  if (outp) then
!  CALL mpimerge3d(pgrd,nx,ny,nz,outlg)
!  if (mp_opt >0 .and. myproc == 0) ThEN
!    DO k = 1,nz
!     do j = 1,nylg
!    write(200,*) k,j
!    write(200,*) outlg(:,j,k)
!     end do
!    write(200,*)
!   end do
!  else if (mp_opt == 0) then
!
!    DO k = 1,nz
!     do j = 1,nylg
!       write(200,*) k,j
!       write(200,*) outlg(:,j,k)
!     end do
!       write(200,*)
!     end do
!   end if
!  call arpsstop(' ',0)
!   end if

  RETURN
END SUBROUTINE recurfilt3d

SUBROUTINE set_bdyxs(pgrd,nx,ny,nz,n,alpha,beta)

! Set XS boundary
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: nx, ny, nz, n
  REAL,    INTENT(IN)  :: alpha(nx,ny,nz)
  REAL,    INTENT(IN)  ::  beta(nx,ny,nz)
  REAL,    INTENT(INOUT) :: pgrd(nx,ny,nz)

  ! local variables

  REAL    :: temp1, temp2
  INTEGER :: k,j

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!XS

      SELECT CASE (n)
        CASE (1)
          DO k = 1,nz-1
            DO j = 1,ny-1
              pgrd(1,j,k) = beta(1,j,k)*pgrd(1,j,k)
            END DO
          END DO
        CASE (2)
          DO k = 1,nz-1
            DO j = 1,ny-1
              temp1 = (1.-alpha(1,j,k)*alpha(1,j,k))
              pgrd(1,j,k) = beta(1,j,k)/temp1 * pgrd(1,j,k)
            END DO
          END DO
        CASE (3)
          DO k = 1,nz-1
            DO j = 1,ny-1
              temp1 = (1-alpha(1,j,k))/                                            &
                  ((1-alpha(1,j,k)*alpha(1,j,k))*(1-alpha(1,j,k)*alpha(1,j,k)))
              temp2 = alpha(1,j,k)*alpha(1,j,k)*alpha(1,j,k)
              pgrd(1,j,k) = temp1 * (pgrd(1,j,k)-temp2*pgrd(2,j,k))
            END DO
          END DO
        CASE DEFAULT
          DO k = 1,nz-1
            DO j = 1,ny-1
              temp2 = alpha(1,j,k)*alpha(1,j,k)*alpha(1,j,k)
              temp1 = (1-alpha(1,j,k))/                                            &
                  (1-3*alpha(1,j,k)*alpha(1,j,k)+3*temp2*alpha(1,j,k)-temp2*temp2)
              pgrd(1,j,k) = temp1 * (pgrd(1,j,k)-3*temp2*pgrd(2,j,k)+              &
                            temp2*alpha(1,j,k)*alpha(1,j,k)*pgrd(2,j,k)+           &
                            temp2*alpha(1,j,k)*pgrd(3,j,k))
            END DO
          END DO
      END SELECT

   RETURN
END SUBROUTINE set_bdyxs

SUBROUTINE set_bdyxe(pgrd,nx,ny,nz,n,alpha,beta)

  INTEGER, INTENT(IN)  :: nx, ny, nz, n
  REAL,    INTENT(IN)  :: alpha(nx,ny,nz)
  REAL,    INTENT(IN)  ::  beta(nx,ny,nz)

  REAL,    INTENT(INOUT)  :: pgrd(nx,ny,nz)

  ! local variables

  REAL    :: temp1, temp2
  INTEGER :: k,j

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!XE
  SELECT CASE (n)
    CASE (0)
      DO k = 1,nz-1
        DO j = 1,ny-1
          pgrd(nx-1,j,k) = beta(nx-1,j,k) * pgrd(nx-1,j,k)
        END DO
      END DO
    CASE (1)
      DO k = 1,nz-1
        DO j = 1,ny-1
          temp1 = 1-alpha(nx-1,j,k)*alpha(nx-1,j,k)
          pgrd(nx-1,j,k) = beta(nx-1,j,k)/temp1 * pgrd(nx-1,j,k)
        END DO
      END DO
    CASE (2)
      DO k = 1,nz-1
        DO j = 1,ny-1
          temp1 = (1-alpha(nx-1,j,k))/                                          &
     ((1-alpha(nx-1,j,k)*alpha(nx-1,j,k))*(1-alpha(nx-1,j,k)*alpha(nx-1,j,k)))
          temp2 = alpha(nx-1,j,k)*alpha(nx-1,j,k)*alpha(nx-1,j,k)
          pgrd(nx-1,j,k) = temp1 * (pgrd(nx-1,j,k)-temp2*pgrd(nx-2,j,k))
        END DO
      END DO

  CASE DEFAULT
    DO k = 1,nz-1
      DO j = 1,ny-1
        temp2 = alpha(nx-1,j,k)*alpha(nx-1,j,k)*alpha(nx-1,j,k)
        temp1 = (1-alpha(nx-1,j,k))/                                            &
      (1-3*alpha(nx-1,j,k)*alpha(nx-1,j,k)+3*temp2*alpha(nx-1,j,k)-temp2*temp2)
        pgrd(nx-1,j,k) = temp1 * (pgrd(nx-1,j,k)-3*temp2*pgrd(nx-2,j,k)+        &
             temp2*alpha(nx-1,j,k)*alpha(nx-1,j,k)*pgrd(nx-2,j,k)+              &
             temp2*alpha(nx-1,j,k)*pgrd(nx-3,j,k))
      END DO
    END DO
  END SELECT

  RETURN
END SUBROUTINE set_bdyxe

SUBROUTINE set_bdyys(pgrd,nx,ny,nz,n,alpha,beta)

! Set XS boundary
  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: nx, ny, nz, n
  REAL,    INTENT(IN)    :: alpha(nx,ny,nz)
  REAL,    INTENT(IN)    ::  beta(nx,ny,nz)
  REAL,    INTENT(INOUT) :: pgrd(nx,ny,nz)

  ! local variables

  REAL    :: temp1, temp2
  INTEGER :: k,i

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!YS

  SELECT CASE (n)
    CASE (1)
      DO k = 1,nz-1
        DO i = 1,nx-1
          pgrd(i,1,k) = beta(i,1,k) * pgrd(i,1,k)
        END DO
      END DO

    CASE (2)
      DO k = 1,nz-1
        DO i = 1,nx-1
          temp1 = (1.-alpha(i,1,k)*alpha(i,1,k))
          pgrd(i,1,k) = beta(i,1,k)/temp1 * pgrd(i,1,k)
        END DO
      END DO

    CASE (3)
      DO k = 1,nz-1
        DO i = 1,nx-1
          temp1 = (1-alpha(i,1,k))/                                            &
               ((1-alpha(i,1,k)*alpha(i,1,k))*(1-alpha(i,1,k)*alpha(i,1,k)))
          temp2 = alpha(i,1,k)*alpha(i,1,k)*alpha(i,1,k)
          pgrd(i,1,k) = temp1 * (pgrd(i,1,k)-temp2*pgrd(i,2,k))
        END DO
      END DO

    CASE DEFAULT
      DO k = 1,nz-1
        DO i = 1,nx-1
          temp2 = alpha(i,1,k)*alpha(i,1,k)*alpha(i,1,k)
          temp1 = (1-alpha(i,1,k))/                                            &
             (1-3*alpha(i,1,k)*alpha(i,1,k)+3*temp2*alpha(i,1,k)-temp2*temp2)
          pgrd(i,1,k) = temp1 * (pgrd(i,1,k)-3*temp2*pgrd(i,2,k)+              &
              temp2*alpha(i,1,k)*alpha(i,1,k)*pgrd(i,2,k)+                     &
              temp2*alpha(i,1,k)*pgrd(i,3,k))
        END DO
      END DO

  END SELECT

  RETURN
END SUBROUTINE set_bdyys

SUBROUTINE set_bdyye(pgrd,nx,ny,nz,n,alpha,beta)

  INTEGER, INTENT(IN)  :: nx, ny, nz, n
  REAL,    INTENT(IN)  ::  alpha(nx,ny,nz)
  REAL,    INTENT(IN)  ::   beta(nx,ny,nz)
  REAL,    INTENT(INOUT) :: pgrd(nx,ny,nz)

  ! local variables

  REAL    :: temp1, temp2
  INTEGER :: k,i

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!YE
  SELECT CASE ( n )

    CASE (0)
      DO k = 1,nz-1
        DO i = 1,nx-1
          pgrd(i,ny-1,k) = beta(i,ny-1,k) * pgrd(i,ny-1,k)
        END DO
      END DO

    CASE ( 1)
      DO k = 1,nz-1
        DO i = 1,nx-1
          temp1 = (1.-alpha(i,ny-1,k)*alpha(i,ny-1,k))
          pgrd(i,ny-1,k) = beta(i,ny-1,k)/temp1 * pgrd(i,ny-1,k)
        END DO
      END DO

    CASE ( 2)
      DO k = 1,nz-1
        DO i = 1,nx-1
         temp1 = (1-alpha(i,ny-1,k))/                                                 &
           ((1-alpha(i,ny-1,k)*alpha(i,ny-1,k))*(1-alpha(i,ny-1,k)*alpha(i,ny-1,k)))
         temp2 = alpha(i,ny-1,k)*alpha(i,ny-1,k)*alpha(i,ny-1,k)
         pgrd(i,ny-1,k) = temp1 * (pgrd(i,ny-1,k)-temp2*pgrd(i,ny-2,k))
        END DO
      END DO

    CASE DEFAULT
      DO k = 1,nz-1
        DO i = 1,nx-1
          temp2 = alpha(i,ny-1,k)*alpha(i,ny-1,k)*alpha(i,ny-1,k)
          temp1 = (1-alpha(i,ny-1,k))/(1-3*alpha(i,ny-1,k)*alpha(i,ny-1,k)+           &
                   3*temp2*alpha(i,ny-1,k)-temp2*temp2)
          pgrd(i,ny-1,k) = temp1 * (pgrd(i,ny-1,k)-3*temp2*pgrd(i,ny-2,k)+            &
               temp2*alpha(i,ny-1,k)*alpha(i,ny-1,k)*pgrd(i,ny-2,k)+                  &
               temp2*alpha(i,ny-1,k)*pgrd(i,ny-3,k))
        END DO
      END DO

  END SELECT

  RETURN
END SUBROUTINE  set_bdyye
