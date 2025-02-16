!
!-----------------------------------------------------------------------
!
!  Time level declarations
!
!-----------------------------------------------------------------------
!
  INTEGER :: nt                ! The no. of time levels of time-dependent
                               ! arrays.
  INTEGER :: tpast             ! Index of time level for the past time.
  INTEGER :: tpresent          ! Index of time level for the present time.
  INTEGER :: tfuture           ! Index of time level for the future time.

  PARAMETER (nt=3, tpast=1, tpresent=2, tfuture=3)

