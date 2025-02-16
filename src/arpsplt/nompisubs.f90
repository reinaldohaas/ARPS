SUBROUTINE mpmax0i
  RETURN
END SUBROUTINE mpmax0i

SUBROUTINE globalpbar(pbarmax,ini,inj,klvl,zpc,nx,ny,nz,zpcmax)

  IMPLICIT NONE

  REAL,    INTENT(IN)  :: pbarmax
  INTEGER, INTENT(IN)  :: ini
  INTEGER, INTENT(IN)  :: inj
  INTEGER, INTENT(IN)  :: klvl
  INTEGER, INTENT(IN)  :: nx,ny,nz
  REAL,    INTENT(IN)  :: zpc(nx,ny,nz)
  REAL,    INTENT(OUT) :: zpcmax

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code below ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF (ini /= 0 .AND. inj /= 0) THEN
    zpcmax = zpc(ini,inj,klvl)
  ELSE
    zpcmax = -9999.0   ! missing value, will not be used
  END IF

  RETURN
END SUBROUTINE globalpbar
