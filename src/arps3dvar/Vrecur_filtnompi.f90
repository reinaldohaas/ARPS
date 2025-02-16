SUBROUTINE send_next_bdyxs(tem1,nx,ny,nz,bdyxs)
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: nx, ny, nz
  REAL,    INTENT(IN)  :: tem1(nx,ny,nz)
  REAL,    INTENT(OUT) :: bdyxs(ny,nz)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  RETURN
END SUBROUTINE send_next_bdyxs

SUBROUTINE receive_bdyxs(tem1,nx,ny,nz,bdyxs)
  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: nx, ny, nz
  REAL,    INTENT(INOUT) :: tem1(nx,ny,nz)
  REAL,    INTENT(OUT)   :: bdyxs(ny,nz)


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  RETURN
END SUBROUTINE receive_bdyxs

! XE
SUBROUTINE send_previous_bdyxe(tem1,nx,ny,nz,bdyxe)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx, ny, nz
  REAL,    INTENT(IN) :: tem1(nx,ny,nz)
  REAL,    INTENT(OUT) :: bdyxe(ny,nz)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  RETURN
END SUBROUTINE send_previous_bdyxe

SUBROUTINE receive_bdyxe(tem1,nx,ny,nz,bdyxe)
  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: nx, ny, nz
  REAL,    INTENT(INOUT) :: tem1(nx,ny,nz)
  REAL,    INTENT(OUT)   :: bdyxe(ny,nz)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  RETURN
END SUBROUTINE receive_bdyxe

SUBROUTINE send_up_bdyys(tem1,nx,ny,nz,bdyys)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx, ny, nz
  REAL,    INTENT(IN) :: tem1(nx,ny,nz)
  REAL,    INTENT(OUT) :: bdyys(nx,nz)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  RETURN
END SUBROUTINE send_up_bdyys

SUBROUTINE receive_bdyys(tem1,nx,ny,nz,bdyys)
  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: nx, ny, nz
  REAL,    INTENT(INOUT) :: tem1(nx,ny,nz)
  REAL,    INTENT(OUT)   :: bdyys(nx,nz)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  RETURN
END SUBROUTINE receive_bdyys

! yE
SUBROUTINE send_down_bdyye(tem1,nx,ny,nz,bdyye)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx, ny, nz
  REAL,    INTENT(IN) :: tem1(nx,ny,nz)
  REAL,    INTENT(OUT) :: bdyye(nx,nz)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  RETURN
END SUBROUTINE send_down_bdyye

SUBROUTINE receive_bdyye(tem1,nx,ny,nz,bdyye)
  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: nx, ny, nz
  REAL,    INTENT(INOUT) :: tem1(nx,ny,nz)
  REAL,    INTENT(OUT)   :: bdyye(nx,nz)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  RETURN
END SUBROUTINE receive_bdyye
