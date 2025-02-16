!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE VERIF_SORT                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE VERIF_SORT(model_data,obsrv_data,tem1,tem2,sfcmax,nhisfile,  &
  nfields,sfcstid,sfcstn,sfcstid_lcl,sfcstn_lcl)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Sort the observations to guarantee that both the non-MPI runs and MPI
!  runs produce identical output files.  Station order really doesn't
!  matter, however, users will complain if they think they have different
!  results.
!
!  Note:  this routine is only called in the MPI case, as by definition,
!  the non-MPI case is already sorted.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Kevin W. Thomas
!
!  Original Coding: 08/03/05
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: sfcmax                     ! number of stations possible
  INTEGER :: nhisfile                   ! number of history files
  INTEGER :: nfields                    ! number of data fields, currently 5
  INTEGER :: sfcstn                     ! number of stations in master list
  INTEGER :: sfcstn_lcl                 ! number of stations in local list

  REAL :: model_data(sfcmax,nhisfile,nfields) ! output model_data array
  REAL :: obsrv_data(sfcmax,nhisfile,nfields) ! output obsrv_data array
  REAL :: tem1(sfcmax,nhisfile,nfields)       ! input model_data array
  REAL :: tem2(sfcmax,nhisfile,nfields)       ! input obsrv_data array

  CHARACTER(LEN=4) :: sfcstid(sfcstn)         ! master output station list
  CHARACTER(LEN=4) :: sfcstid_lcl(sfcstn_lcl) ! local station list

!
! Note:  "tem1" and "model_data" had better not be the say array!
!        The same holds true for "tem2" and "obsrv_data"!
!

  INTEGER :: nstart                     ! starting point in output arrays
  INTEGER :: kount                      ! number of stations to process
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!  None.
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Misc. internal variables
!
!-----------------------------------------------------------------------
!

  INTEGER :: i,j,k,l,istn

  DO l=1,sfcstn_lcl
    DO i=1,sfcstn
      IF (sfcstid(i) == sfcstid_lcl(l)) EXIT
    END DO

    IF (i > sfcstn) THEN
      write(6,*) 'verif_sort:  can''t find ',sfcstid_lcl(l),' in master list???'
      exit
    END IF

    DO k=1,nfields
      DO j=1,nhisfile

        model_data(i,j,k) = tem1(l,j,k)
        obsrv_data(i,j,k) = tem2(l,j,k)

      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE verif_sort
