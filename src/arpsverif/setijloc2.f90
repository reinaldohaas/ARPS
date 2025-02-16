!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE SETIJLOC2                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE setijloc2(fnx,fny,anx,any,iloc,jloc,fx2d,fy2d,ax2d,ay2d)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Find i,j indicies in forecast grid of each analysis grid point.
! No indicies are returned if the analysis grid point is outside
! of the forecast grid.
!
! NOTE:  The forecast and analysis grid points are assumed to have
! the same origin and map projection.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Eric Kemp, November 1999.
!
! Based on subroutine setijloc
!
!-----------------------------------------------------------------------
!
! Variable declarations
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER,INTENT(IN) :: fnx,fny,anx,any
  REAL,INTENT(IN) :: fx2d(fnx,fny),fy2d(fnx,fny)
  REAL,INTENT(IN) :: ax2d(anx,any),ay2d(anx,any)
  INTEGER,INTENT(INOUT) :: iloc(anx,any),jloc(anx,any)
                       
!-----------------------------------------------------------------------
!
! Miscellaneous variables
!
!-----------------------------------------------------------------------

  INTEGER :: i,j,m,n
  INTEGER :: fimid,fjmid
  REAL :: fxmid,fymid

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
! Initialize iloc and jloc to -9999 (i.e., outside the forecast grid)
!
!-----------------------------------------------------------------------

  DO n = 1,any
    DO m = 1,anx
      iloc(m,n) = -9999
      jloc(m,n) = -9999
    END DO
  END DO

!-----------------------------------------------------------------------
!
! Determine the scalar coordinates of the middle of the forecast field.
!
!-----------------------------------------------------------------------

  fimid = fnx/2
  fjmid = fny/2

  fxmid = fx2d(fimid,fjmid)
  fymid = fy2d(fimid,fjmid)

!-----------------------------------------------------------------------
!
! Find the forecast i,j of each analysis grid point.
!
!-----------------------------------------------------------------------

  DO n = 1,any-1
    analysism: DO m = 1,anx-1
      IF (ax2d(m,n) < fxmid) THEN
        IF (ay2d(m,n) < fymid) THEN
          DO j = fjmid,2,-1
            DO i = fimid,2,-1
              IF (fx2d(i,j) <= ax2d(m,n) .AND.                       &
                  fx2d(i+1,j) > ax2d(m,n) .AND.                    &
                  fy2d(i,j) <= ay2d(m,n) .AND.                       &
                  fy2d(i,j+1) > ay2d(m,n)) THEN
                iloc(m,n) = i
                jloc(m,n) = j
                CYCLE analysism
              END IF
            END DO
          END DO
        ELSE
          DO j = fjmid,fny-1
            DO i = fimid,2,-1
              IF (fx2d(i,j) <= ax2d(m,n) .AND.                       &
                  fx2d(i+1,j) > ax2d(m,n) .AND.                    &
                  fy2d(i,j-1) < ay2d(m,n) .AND.                       &
                  fy2d(i,j) >= ay2d(m,n)) THEN
                iloc(m,n) = i
                jloc(m,n) = j-1
                CYCLE analysism
              END IF
            END DO
          END DO
        END IF
      ELSE
        IF (ay2d(m,n) < fymid) THEN
          DO j = fjmid,2,-1
            DO i = fimid,fnx-1
              IF (fx2d(i-1,j) < ax2d(m,n) .AND.                       &
                  fx2d(i,j) >= ax2d(m,n) .AND.                        &
                  fy2d(i,j) <= ay2d(m,n) .AND.                        &
                  fy2d(i,j+1) > ay2d(m,n)) THEN
                iloc(m,n) = i-1
                jloc(m,n) = j
                CYCLE analysism
              END IF
            END DO
          END DO    
        ELSE
          DO j = fjmid,fny-1
            DO i = fimid,fnx-1 
              IF (fx2d(i-1,j) < ax2d(m,n) .AND.                       &
                  fx2d(i,j) >= ax2d(m,n) .AND.                        &
                  fy2d(i,j-1) < ay2d(m,n) .AND.                       &
                  fy2d(i,j) >= ay2d(m,n)) THEN
                iloc(m,n) = i-1
                jloc(m,n) = j-1
                CYCLE analysism
              END IF
            END DO
          END DO      
        END IF
      END IF
    END DO analysism
  END DO
  
END SUBROUTINE setijloc2
