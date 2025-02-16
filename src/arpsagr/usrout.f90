
SUBROUTINE usrout1( levels,levele, initdata )
!
! To do data output for all grids by calling grid output routine
! arpsout for ARPS. (Ming Xue, 11/6/1992).
!
  INCLUDE 'nodal.inc'
  INCLUDE 'agrialloc.inc'
  INCLUDE 'agrigrid.inc'
  INCLUDE 'agricpu.inc'

  cpu0 = f_cputime()

!-----------------------------------------------------------------------------
! Search for grids from the starting level (levels) to ending level (levele)
! Then search for each grid on a particular level
!-----------------------------------------------------------------------------
  DO level = levels, levele

    IF(level > lfine ) EXIT

    mptr = lstart(level)
    5       IF(mptr == 0) CYCLE

!-----------------------------------------------------------------------------
! Call data output routine for a grid mptr
!-----------------------------------------------------------------------------
    IF(.true.) PRINT*,'calling arpsout for grid ',                      &
               mptr,' at time step  ',node(13,mptr)

    nestgrd = 1
    CALL arpsout(mptr, initdata , nestgrd )
!
    mptr = node(10,mptr)
    GO TO 5
  END DO

  999   CONTINUE

  cpu_usrout = cpu_usrout + f_cputime() - cpu0

  RETURN
END SUBROUTINE usrout1

SUBROUTINE usrout2( levels,levele )
!
! Equivalent to usrout1 for ARPS
!
! Do nothing for the moment
!
  RETURN
END SUBROUTINE usrout2
