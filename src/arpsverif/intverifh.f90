!######################################################################## 
!########################################################################
!######                                                            ######
!######                   SUBROUTINE INTVERIF_H                    ######
!######                                                            ######       
!######                      Developed by                          ######
!######         Center for Analysis and Prediction of Storms       ######
!######                   University of Oklahoma                   ######
!######                                                            ######
!########################################################################
!########################################################################


SUBROUTINE intverif_H(vnx,vny,fnx,fny,fntime,fn_datasets,               &
                      vibeg,viend,vjbeg,vjend,fibeg,fiend,fjbeg,fjend,  &
                      iorder,vx2d,vy2d,fxs,fys,ftem4d,vtem4d)


!-----------------------------------------------------------------------
!   
! PURPOSE:
!
! Reads in forecast fields, and horizontally interpolates them to
! a second verification grid.
!
! AUTHOR:  Eric Kemp, November 1999
!
!----------------------------------------------------------------------- 
!
!----------------------------------------------------------------------- 
!
! Use modules
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 
! Variable declarations
! 
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: vnx,vny,fntime,fnx,fny,fn_datasets

  INTEGER :: vibeg,viend,vjbeg,vjend
  INTEGER :: fibeg,fiend,fjbeg,fjend

  REAL,DIMENSION(fnx,fny,fntime,fn_datasets):: ftem4d
  REAL,DIMENSION(vnx,vny,fntime,fn_datasets):: vtem4d


  REAL :: dxfld(fnx),dyfld(fny),rdxfld(fnx),rdyfld(fny)
  REAL,DIMENSION(fnx,fny) :: slopey,alphay,betay
  INTEGER :: iorder

  REAL,DIMENSION(vnx,vny) :: vx2d,vy2d
  REAL :: fxs(fnx),fys(fny)


!-----------------------------------------------------------------------
! 
! Miscellaneous variables
! 
!-----------------------------------------------------------------------

  INTEGER,DIMENSION(fnx,fny) :: iloc,jloc

  REAL :: ftem2d(fnx,fny),vtem2d(vnx,vny)
  INTEGER :: i,j,k,l

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
! 
! Find i,j indices in forecast grid of each verification point.
! 
!-----------------------------------------------------------------------

  CALL setijloc2(vnx,vny,fnx,fny,vx2d,vy2d,fxs,fys,iloc,jloc)
  CALL setdxdy(fnx,fny,1,fnx-1,1,fny-1,fxs,fys,dxfld,dyfld,             &
               rdxfld,rdyfld)

!-----------------------------------------------------------------------
! 
! Interpolate 2-D fields
! 
!-----------------------------------------------------------------------

  DO l = 1,fn_datasets
    DO k = 1,fntime
      DO j = 1,fny
        DO i = 1,fnx
          ftem2d(i,j) = ftem4d(i,j,k,l)
        END DO
      END DO

      WRITE(6,*)'Calling fldint2d...'
      CALL fldint2d(vnx,vny,fnx,fny,                                   &
                    vibeg,viend,vjbeg,vjend,                           &
                    fibeg,fiend,fjbeg,fjend,                           &
                    iorder,vx2d,vy2d,ftem2d,fxs,fys,iloc,jloc,         &
                    dxfld,dyfld,rdxfld,rdyfld,                         &
                    slopey,alphay,betay,                               &
                    vtem2d)                                            

      DO j = 1,vny
        DO i = 1,vnx
          vtem4d(i,j,k,l) = vtem2d(i,j)
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE intverif_H
