!########################################################################
!########################################################################
!#########                                                      #########
!#########                SUBROUTINE MPIPROCESS                 #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE mpiprocess(nobmpi,indexmpi,useobs,np,kitem,kitemmax,         &
                      isrc,item1,nx,ny,xmpi,ympi,xs,ys,iflag)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Determine which processor "owns" an ob.  Later, when we need info,
! we know who to "contact".
!
!-----------------------------------------------------------------------
!
! AUTHOR: Kevin W. Thomas, CAPS
! December 2, 2005
!
! MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
! Force explicit declarations
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
! Subroutine arguments
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nobmpi         ! Total number of single-point
                                        ! observations stored in arrays
  INTEGER, INTENT(INOUT) :: indexmpi(nobmpi) ! Processor Owner of each observation
                                        ! All processor contains the same values
  LOGICAL, INTENT(OUT)   :: useobs(nobmpi)   ! Valid observation for each
                                        ! processor. All processor has its own copy.
  INTEGER, INTENT(IN)    :: np          ! Number of processors (nprocs)
  INTEGER, INTENT(INOUT) :: kitem(np)   ! Number obs handled by each processor
  INTEGER, INTENT(INOUT) :: kitemmax    ! Largest "kitem" value
  INTEGER, INTENT(INOUT) :: isrc(nobmpi)   ! Data source number
  INTEGER, INTENT(INOUT) :: item1(nobmpi)  ! Work array
  REAL,    INTENT(INOUT) :: xmpi(nobmpi)   ! Observation x grid coordinate (m)
  REAL,    INTENT(INOUT) :: ympi(nobmpi)   ! Observation y grid coordinate (m)
  INTEGER, INTENT(IN)    :: nx,ny       ! Grid dimensions.
  REAL,    INTENT(IN)    :: xs(nx)      ! x-coordinates of grid scalar points (m)  
  REAL,    INTENT(IN)    :: ys(ny)      ! y-coordinates of grid scalar points (m)
  INTEGER, INTENT(IN)    :: iflag       ! if radar processing, skip some code

!-----------------------------------------------------------------------
!
! Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER :: ksta
  INTEGER :: ipt,jpt
  INTEGER :: indom
  INTEGER :: ierror

  INTEGER :: k

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
! Single/multi-level data need to find if a point is inside of a domain
! so a processor can claim the data.  For radar, this has already been
! done, so skip this part and proceed to the merge.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  IF (iflag == 0) THEN
    DO ksta=1,nobmpi

      indexmpi(ksta) = -1
      useobs(ksta)   = .FALSE.
      
      !
      ! If we've been declared "bad" or joined via "suprob", then stop here.
      !
      IF (iflag == 0 ) THEN
        IF ( isrc(ksta) == 0 ) CYCLE
      END IF

      CALL findlc(nx,ny,xs,ys,xmpi(ksta),ympi(ksta),ipt,jpt,indom)

!      IF (indom >= 0 )  useobs(ksta) = .TRUE.    ! use more obs (should be put
                                                 ! inside CASE statement below
      SELECT CASE (indom)    ! deterministic ownership
      CASE (0) 
        indexmpi(ksta) = myproc    ! own it and
        useobs(ksta)   = .TRUE.    ! also use it

      ! will be used without conditions, but only be owned by specific processor
      CASE (1)       ! west boundary
        IF (loc_x == 1)       indexmpi(ksta) = myproc
        useobs(ksta)   = .TRUE.
      CASE (10)      ! south boundary
        IF (loc_y == 1)       indexmpi(ksta) = myproc
        useobs(ksta)   = .TRUE.
      CASE (11)      ! southwest corner
        IF (loc_x == 1  .AND. loc_y == 1 )  indexmpi(ksta) = myproc
        useobs(ksta)  = .TRUE.

      ! do not use and be owned except for specifc processor
!      CASE (12)      ! only last processor in X direction use it and also own it
!        IF (loc_x == nproc_x) THEN
!          indexmpi(ksta) = myproc   
!          useobs(ksta)  = .TRUE.
!        END IF
!      CASE (20)     ! only last processor in Y direction use it and also own it
!        IF (loc_y == nproc_y) THEN
!          indexmpi(ksta) = myproc
!          useobs(ksta)  = .TRUE.
!        END IF
!      CASE (21)     ! 
!        IF (loc_x == 1       .AND. loc_y == nproc_y) THEN
!          indexmpi(ksta) = myproc
!          useobs(ksta)  = .TRUE.
!        END IF
!      CASE (22)
!        IF (loc_x == nproc_x .AND. loc_y == 1) THEN
!          indexmpi(ksta) = myproc
!          useobs(ksta)  = .TRUE.
!        END IF
!      CASE (32)
!        IF (loc_x == nproc_x .AND. loc_y == nproc_y) THEN
!          indexmpi(ksta) = myproc
!          useobs(ksta)  = .TRUE.
!        END IF
      CASE DEFAULT
        ! do nothing
      END SELECT

    END DO
!  END IF

  IF (mp_opt > 0) THEN
    DO k=1,nprocs-1
    !
    ! Collect and merge the data.
    !
      IF (myproc > 0 ) THEN
        IF (myproc == k ) CALL mpsendi(indexmpi,nobmpi,0,1000+myproc,ierror)
      ELSE
        CALL mprecvi(item1,nobmpi,k,1000+k,ierror)
        DO ksta=1,nobmpi
          IF ( item1(ksta) == -1 ) CYCLE

!-------------------------------------------------------------------------
!
! Since there are overlapping grids in MPI, it is possible for an ob to
! be available to more than one processor.  We select the first processor.
! There is no need for more than one processor to make identical
! computations.
!
! The WARNING message is commented out, as it is useful for debugging,
! however, it will likely confuse anyone else.
!
!-------------------------------------------------------------------------

!       IF ( indexmpi(ksta) .ne. -1 ) THEN
!         WRITE(6,*) 'WARNING:  station ',ksta,' found in ',             &
!           indexmpi(ksta),' and ',item1(ksta)
!       END IF
          indexmpi(ksta) = item1(ksta)
        END DO
      END IF
     CALL mpbarrier
   END DO

! Dump the station to processor mapping.  Useful only for code debugging.

! if ( myproc == 0 ) then
! write(6,*) 'MAPPING:  '
! do ksta=1,nobmpi
!   write(6,*) ksta,indexmpi(ksta)
! end do
! endif

    CALL mpupdatei(indexmpi,nobmpi)

!
! Everybody computes the same table of how many obs each processor owns.
!

    kitem(:) = 0
    kitemmax = 0

    DO k=1,nobmpi
      IF(indexmpi(k) >= 0) kitem(indexmpi(k)+1) = kitem(indexmpi(k)+1) + 1
      IF( useobs(k) ) THEN
        IF (indexmpi(k) < 0) THEN
          useobs(k) = .FALSE.  ! to avoid using an observation that doesn't 
                               ! own by anyone.
        ELSE
          kitemmax =  kitemmax + 1
        END IF
      END IF
    END DO

    CALL mpmaxi(kitemmax)

!  DO k=1,np
!    IF(kitem(k) > kitemmax) kitemmax = kitem(k)
!  END DO
  ELSE
    DO k=1,nobmpi
      IF( useobs(k) .AND. (indexmpi(k) < 0)) useobs(k) = .FALSE.  
      ! to avoid using an observation that doesn't own by anyone.
    END DO
  END IF

  RETURN
END SUBROUTINE mpiprocess

SUBROUTINE mpiprocess_complex(nobmaxmpi,nobmpi,indexmpi,oindexmpi,useobs, &
                      np,kitem,kitemmax,isrc,item1,nx,ny,xmpi,ympi,       &
                      xs,ys,iflag,istatus)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Determine which processor "owns" an ob.  Later, when we need info,
! we know who to "contact".
!
! This version is more complex than "mpiprocess".  That routine assumes
! that there is a one to one correspondence between all arrays.  This
! version doesn't work this way, as unnecssary are in this case no longer
! fill the arrays.  Instead, there is an obs number mapper, oindexmpi,
! that controls the insanity.  This is motivated by the excessive number
! of radial velocity data points.
!
!-----------------------------------------------------------------------
!
! AUTHOR: Kevin W. Thomas, CAPS
! January 9, 2008
!
! MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
! Force explicit declarations
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
! Subroutine arguments
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nobmaxmpi      ! Maximum number of observations in 
           ! the global domain. Same as "nobmpi" in subroutine mpiprocess.

  INTEGER, INTENT(IN) :: nobmpi         ! Total number of single-point
           ! observations stored in arrays. This is local in this subroutine. It
           ! is global in subroutine mpiprocess.

  INTEGER, INTENT(INOUT) :: indexmpi(nobmpi)  ! Owner
  INTEGER, INTENT(IN)    :: oindexmpi(nobmpi) ! Owner index mapper
  INTEGER, INTENT(IN)    :: np          ! Number of processors (nprocs)
  INTEGER, INTENT(INOUT) :: kitem(np)   ! Number obs handled by each processor
  INTEGER, INTENT(INOUT) :: kitemmax    ! Largest "kitem" value
  INTEGER, INTENT(INOUT) :: isrc(nobmpi) ! Data source number
  INTEGER, INTENT(INOUT) :: item1(nobmpi)    ! Work array
  REAL,    INTENT(INOUT) :: xmpi(nobmpi)   ! Observation x grid coordinate (m)
  REAL,    INTENT(INOUT) :: ympi(nobmpi)   ! Observation y grid coordinate (m)
  INTEGER, INTENT(IN) :: nx,ny          ! Grid dimensions.
  REAL,    INTENT(IN) :: xs(nx) ! x-coordinates of grid scalar points (m)
  REAL,    INTENT(IN) :: ys(ny) ! y-coordinates of grid scalar points (m)
  INTEGER, INTENT(IN) :: iflag          ! if radar processing, skip some code

  LOGICAL, INTENT(OUT) :: useobs(nobmpi)
  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Processors 0 has to have a lot of bookkeeping.
!
!-----------------------------------------------------------------------

  INTEGER, ALLOCATABLE :: nobmpi_global(:)
  INTEGER, ALLOCATABLE :: itmp(:)
  INTEGER, ALLOCATABLE :: indexmpi_global(:)
  INTEGER, ALLOCATABLE :: indexmpi_memory(:,:)
  INTEGER, ALLOCATABLE :: oindexmpi_memory(:,:)

!-----------------------------------------------------------------------
!
! Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER :: ksta
  INTEGER :: ipt,jpt
  INTEGER :: indom
  INTEGER :: ierror

  INTEGER :: i
  INTEGER :: isub
  INTEGER :: k
  INTEGER :: kitemlocal

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
! Single/multi-level data need to find if a point is inside of a domain
! so a processor can claim the data.  For radar, this has already been
! done, so skip this part and proceed to the merge.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  useobs(:) = .FALSE.

!  IF (iflag == 0) THEN    ! we want to recalculate the useobs for 3dVAR
    DO ksta=1,nobmpi

      indexmpi(ksta) = -1

!
!     If we've been declared "bad" or joined via "suprob", then stop here.
!

      IF (iflag == 0) THEN
        IF (isrc(ksta) == 0 ) CYCLE
      END IF

      CALL findlc(nx,ny,xs,ys,xmpi(ksta),ympi(ksta),ipt,jpt,indom)

!      IF (indom >= 0 )  useobs(ksta) = .TRUE.
!      IF (indom >= 0 ) indexmpi(ksta) = myproc

      SELECT CASE (indom)    ! deterministic ownership
      CASE (0) 
        indexmpi(ksta) = myproc    ! own it and
        useobs(ksta)   = .TRUE.    ! also use it

      ! will be used without conditions, but only be owned by specific processor
      CASE (1)       ! west boundary
        IF (loc_x == 1)       indexmpi(ksta) = myproc
        useobs(ksta)   = .TRUE.
      CASE (10)      ! south boundary
        IF (loc_y == 1)       indexmpi(ksta) = myproc
        useobs(ksta)   = .TRUE.
      CASE (11)      ! southwest corner
        IF (loc_x == 1  .AND. loc_y == 1 )  indexmpi(ksta) = myproc
        useobs(ksta)  = .TRUE.

      ! will not be used or be owned except for specifc processor
!      CASE (12)      ! only last processor in X direction use it and also own it
!        IF (loc_x == nproc_x) THEN
!          indexmpi(ksta) = myproc   
!          useobs(ksta)  = .TRUE.
!        END IF
!      CASE (20)     ! only last processor in Y direction use it and also own it
!        IF (loc_y == nproc_y) THEN
!          indexmpi(ksta) = myproc
!          useobs(ksta)  = .TRUE.
!        END IF
!      CASE (21)     ! 
!        IF (loc_x == 1       .AND. loc_y == nproc_y) THEN
!          indexmpi(ksta) = myproc
!          useobs(ksta)  = .TRUE.
!        END IF
!      CASE (22)
!        IF (loc_x == nproc_x .AND. loc_y == 1) THEN
!          indexmpi(ksta) = myproc
!          useobs(ksta)  = .TRUE.
!        END IF
!      CASE (32)
!        IF (loc_x == nproc_x .AND. loc_y == nproc_y) THEN
!          indexmpi(ksta) = myproc
!          useobs(ksta)  = .TRUE.
!        END IF
      CASE DEFAULT
        ! do nothing
      END SELECT

    END DO
!  END IF

!
! Collect the number of obs for each processor.
!

  IF (myproc == 0) THEN
      ALLOCATE(nobmpi_global(0:nprocs-1),STAT=istatus)
      CALL check_alloc_status(istatus,"mpiprocess_complex:nobmpi_global")
      nobmpi_global(0) = nobmpi
      ksta = nobmpi_global(0)
  END IF

  DO k=1,nprocs-1
    IF (myproc > 0) THEN
      IF (myproc == k) CALL mpsendi(nobmpi,1,0,2000+myproc,ierror)
    ELSE
!
! "ksta" is used to store the max number of obs, so the arrays we allocate
! are only as large as necessary.
!
      CALL mprecvi(nobmpi_global(k),1,k,2000+k,ierror)
      IF (nobmpi_global(k) > ksta) ksta = nobmpi_global(k)
    END IF
    CALL mpbarrier
  END DO

  IF (myproc == 0) THEN

!
! This will send/receive data.
!

    ALLOCATE(itmp(ksta),STAT=istatus)
    CALL check_alloc_status(istatus,"mpiprocess_complex:itmp")

!
! Arrays for indexmpi and friends.
!

    ALLOCATE(indexmpi_global(nobmaxmpi),STAT=istatus)
    CALL check_alloc_status(istatus,"mpiprocess_complex:indexmpi_global")

    ALLOCATE(indexmpi_memory(ksta,0:nprocs-1),STAT=istatus)
    CALL check_alloc_status(istatus,"mpiprocess_complex:indexmpi_memory")

!
! Array for oindexmpi.
!

    ALLOCATE(oindexmpi_memory(ksta,0:nprocs-1),STAT=istatus)
    CALL check_alloc_status(istatus,"mpiprocess_complex:oindexmpi_memory")

  END IF

!
!  Collect the the indexmpi/oindexmpi values.  Processors that don't have
!  any values are skipped.
!

  DO k=1,nprocs-1
    IF (myproc > 0) THEN
      IF (nobmpi > 0 .AND. myproc == k) THEN
        CALL mpsendi( indexmpi,nobmpi,0,3000+myproc,ierror)
        CALL mpsendi(oindexmpi,nobmpi,0,4000+myproc,ierror)
      END IF
    ELSE
      IF (nobmpi_global(k) > 0) THEN
        CALL mprecvi(itmp,nobmpi_global(k),k,3000+k,ierror)
        DO i=1,nobmpi_global(k)
          indexmpi_memory(i,k) = itmp(i)
        END DO

        CALL mprecvi(itmp,nobmpi_global(k),k,4000+k,ierror)
        DO i=1,nobmpi_global(k)
          oindexmpi_memory(i,k) = itmp(i)
        END DO
      END IF
    ENDIF
    CALL mpbarrier
  END DO

!
! Fill in our own information.
!

  IF (myproc == 0 .AND. nobmpi > 0) THEN
    DO i=1,nobmpi
      indexmpi_memory (i,0) =  indexmpi(i)
      oindexmpi_memory(i,0) = oindexmpi(i)
    END DO
  END IF

  IF (myproc == 0) THEN
    DO i=1,nobmaxmpi
      indexmpi_global(i) = -1
    END DO

!
!  Search and merge the info.  Unfortunately, all processors but zero are
!  idle, though this should be fast.
!

    DO k=0,nprocs-1
      IF (nobmpi_global(k) == 0) CYCLE
      DO i=1,nobmpi_global(k)
        isub = oindexmpi_memory(i,k)
        IF (indexmpi_memory(i,k) .NE. k) CYCLE
!
!   No complaints (not even for debug right now) for duplicates which will
!   always happen with overlap points.
!
        IF (indexmpi_global(isub) == -1)                                 &
          indexmpi_global(isub) = indexmpi_memory(i,k)
      END DO
    END DO
  END IF

!
! Rewrite the mapping tables, and send them out.
!

  DO k=1,nprocs-1
    IF (myproc > 0) THEN
       IF (nobmpi > 0 .AND. myproc == k) THEN
         CALL mprecvi(indexmpi,nobmpi,0,5000+myproc,ierror)
       END IF
    ELSE
       IF (nobmpi_global(k) > 0) THEN

         DO i=1,nobmpi_global(k)
           isub = oindexmpi_memory(i,k)
           itmp(i) = indexmpi_global(isub)
         END DO
         CALL mpsendi(itmp,nobmpi_global(k),k,5000+k,ierror)
       END IF
    ENDIF
    CALL mpbarrier
  END DO

!
! Fill in our own information.
!

  IF (myproc == 0 .AND. nobmpi > 0) THEN
    DO i=1,nobmpi
      isub = oindexmpi_memory(i,0)
      indexmpi(i) = indexmpi_global(isub)
    END DO
  END IF

!
!   Collect the "ncolitem" values so we can compute the kitem/kitemmax info.
!

  kitemlocal = 0
  DO k=1,nobmpi
    IF(indexmpi(k) == myproc) kitemlocal = kitemlocal + 1
  END DO

  IF (myproc == 0) THEN
      kitem(1) = kitemlocal
      kitemmax = kitem(1)
  END IF

  DO k=1,nprocs-1
    IF (myproc > 0) THEN
      IF (myproc == k) CALL mpsendi(kitemlocal,1,0,6000+myproc,ierror)
    ELSE
      CALL mprecvi(kitem(k+1),1,k,6000+k,ierror)
      IF (kitem(k+1) > kitemmax) kitemmax = kitem(k+1)
    END IF
    CALL mpbarrier
  END DO

  CALL mpupdatei(kitem,nprocs)
  CALL mpupdatei(kitemmax,1)

  IF (myproc == 0) THEN
    DEALLOCATE(nobmpi_global)
    DEALLOCATE(indexmpi_global)
    DEALLOCATE(indexmpi_memory)
  END IF

  RETURN
END SUBROUTINE mpiprocess_complex

!########################################################################
!########################################################################
!#########                                                      #########
!#########              SUBROUTINE MPIPROCESS_UPDATE            #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE mpiprocess_update(nobmpi,indexmpi)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Update the communications index table if it has been altered.  This
! is a hook for cloud soundings.
!
!-----------------------------------------------------------------------
!
! AUTHOR: Kevin W. Thomas, CAPS
! December 2, 2005
!
! MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
! Force explicit declarations
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
! Subroutine arguments
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nobmpi         ! Total number of single-point
                                        ! observations stored in arrays
  INTEGER, INTENT(INOUT) :: indexmpi(nobmpi) ! Owner
  INTEGER :: item1(nobmpi)

!-----------------------------------------------------------------------
!
! Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER :: ksta
  INTEGER :: ipt,jpt
  INTEGER :: indom
  INTEGER :: ierror

  INTEGER :: k

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!
! Collect and merge the data.
!

  DO k=1,nprocs-1
    IF (myproc > 0) THEN
      IF (myproc == k) CALL mpsendi(indexmpi,nobmpi,0,1000+myproc,ierror)
    ELSE
      CALL mprecvi(item1,nobmpi,k,1000+k,ierror)
      DO ksta=1,nobmpi

        IF ( item1(ksta) == -1 ) CYCLE

!-------------------------------------------------------------------------
!
! Since there are overlapping grids in MPI, it is possible for an ob to
! be available to more than one processor.  We select the first processor.
! There is no need for more than one processor to make identical
! computations.
!
! The WARNING message is commented out, as it is useful for debugging,
! however, it will likely confuse anyone else.
!
! There is a potential issue affecting cloud soundings that doesn't exist
! in "mpiprocess" above.  Stations *outside* of the domain are permitted.
! These are assigned a processor inside of "insert_sao1".  It *is* possible
! for a point to be assigned to more than one processors, however, only one
! of them is correct, if the non-MPI and MPI solutions are to be the same.
!
! For now, we'll use the rule that processor 0 will always override any
! other ownership claims.  This is based on one case.  If other problems
! are seen, something else will have to be done.
!
!-------------------------------------------------------------------------

        IF ( indexmpi(ksta) .eq. 0 ) CYCLE

!       IF ( indexmpi(ksta) .ne. -1 ) THEN
!         WRITE(6,*) 'WARNING:  station ',ksta,' found in ',             &
!           indexmpi(ksta),' and ',item1(ksta)
!       END IF
        indexmpi(ksta) = item1(ksta)
      END DO
    END IF
    CALL mpbarrier
  END DO

! Dump the station to processor mapping.  Useful only for code debugging.

! if ( myproc == 0 ) then
! write(6,*) 'MAPPING:  '
! do ksta=1,nobmpi
!   write(6,*) ksta,indexmpi(ksta)
! end do
! endif

  CALL mpupdatei(indexmpi,nobmpi)

  RETURN
END SUBROUTINE mpiprocess_update

SUBROUTINE make_mpi_map(mpi_map,nmap,iproc,jproc,nx,ny)

!-----------------------------------------------------------------------
!
! PURPOSE:
! Build a map of who needs to communicate with other processors.
!
!-----------------------------------------------------------------------
!
! AUTHOR: Kevin W. Thomas, CAPS
! October 6, 2006
!
! MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
! Force explicit declarations
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
! Subroutine arguments
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nmap           ! Number of map entries
  INTEGER, INTENT(INOUT) :: mpi_map(nmap,2) ! The map
  INTEGER, INTENT(IN) :: iproc          ! Number of x-direction offsets
  INTEGER, INTENT(IN) :: jproc          ! Number of y-direction offsets
  INTEGER, INTENT(IN) :: nx             ! Number of x-direction grid pts
  INTEGER, INTENT(IN) :: ny             ! Number of y-direction grid pts

!-----------------------------------------------------------------------
!
! Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER :: i,j,k,n
  INTEGER :: lx, ly

  IF (mp_opt == 0) RETURN

  mpi_map = -1

!
! Receive map.
!

  k = 0
  DO j=-jproc,jproc
    ly = loc_y + j
    IF (ly < 1 .OR. ly > nproc_y) THEN
      k = k + (2 * iproc + 1 )
      CYCLE
    END IF
    DO i=-iproc,iproc
      k = k + 1
      IF (i == 0 .and. j == 0) CYCLE
      lx = loc_x + i
      IF (lx < 1 .OR. lx > nproc_x) CYCLE
      n = myproc + j * nproc_x + i
      IF (n >= nprocs) CYCLE
      mpi_map(k,2) = n
    END DO
  END DO

!
! Send map, just flip flop the receive map.
!

   DO k=1,nmap
     mpi_map(k,1) = mpi_map(nmap-k+1,2)
   END DO
 
END SUBROUTINE make_mpi_map

SUBROUTINE data_filt(nobmpi,isrc,indexmpi,nmap,mpi_map)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! In MPI mode, if a processor will never need an ob, mark it not to be
! used to reduce computation time.  Note that the cloud analysis
! (cloudopt=1) may correctly still use the ob.
!
!-----------------------------------------------------------------------
!
! AUTHOR: Kevin W. Thomas, CAPS
! November 8, 2007
!
! MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
! Force explicit declarations
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
! Subroutine arguments
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN)    :: nobmpi           ! Actual number of stations
  INTEGER, INTENT(INOUT) :: isrc(nobmpi)     ! Data source number
  INTEGER, INTENT(IN)    :: indexmpi(nobmpi) ! Array saying which ob is on which
                                        ! processor.
  INTEGER, INTENT(IN) :: nmap           ! Number of entries in "mpi_map"
  INTEGER, INTENT(IN) :: mpi_map(nmap,2)! Mapping scheme

!-----------------------------------------------------------------------
!
! Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER, ALLOCATABLE :: item1(:)      ! Temporary work array.
  INTEGER :: k
  INTEGER :: istatus

  ALLOCATE(item1(0:nprocs-1),STAT=istatus)
  CALL check_alloc_status(istatus,"adas:data_filt")
  
  item1 = 0

  DO k=1,nmap
    IF (mpi_map(k,1) .NE. -1) item1(mpi_map(k,1)) = 1
  END DO
  item1(myproc) = 1

  DO k=1,nobmpi
    IF (isrc(k) < 1) CYCLE              ! Already marked not to use
    IF (indexmpi(k) == -1) CYCLE        ! No owner obs also already marked
    IF (item1(indexmpi(k)) == 0) isrc(k) = 0
  END DO
  DEALLOCATE (item1)

  RETURN
END SUBROUTINE data_filt

SUBROUTINE mpi_1di_collect(m, nobmpi, indexmpi, np,                     &
                           kdata, kdatamax, mpi_map, nmap, tmps, tmpr )

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Collect all the processor dependent calculations and merge them into a
! single array which we will broadcast.
!
!-----------------------------------------------------------------------
!
! AUTHOR: Kevin W. Thomas, CAPS
! January 10, 2006
!
! MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
! Force explicit declarations
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
! Subroutine arguments
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nobmpi         ! Actual number of stations
  INTEGER, INTENT(IN) :: indexmpi(nobmpi) ! Array saying which ob is on which
                                        ! processor.
  INTEGER, INTENT(INOUT) :: m(nobmpi)   ! Data to be collected and updated.
  INTEGER, INTENT(IN) :: np             ! Just "nprocs"
  INTEGER, INTENT(IN) :: kdata(np)      ! Number of obs owned by each processor
  INTEGER, INTENT(IN) :: kdatamax       ! Largest "kdata" value
  INTEGER, INTENT(IN) :: nmap           ! Number of entries in "mpi_map"
  INTEGER, INTENT(IN) :: mpi_map(nmap,2)! Mapping scheme

  INTEGER, INTENT(INOUT) :: tmps(kdatamax)
  INTEGER, INTENT(INOUT) :: tmpr(kdatamax)

!-----------------------------------------------------------------------
!
! Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER :: istat
  INTEGER :: isum
  INTEGER :: itag

  INTEGER :: ierror

  INTEGER :: i, j, k, l, ksta

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  isum = kdatamax

!
!  Save our data.
!

  k = 0
  DO j=1,nobmpi
    IF (indexmpi(j) == myproc) THEN
      k=k+1
      tmps(k) = m(j)
    END IF
  END DO

!
! Sanity.
!

  IF ( k .ne. kdata(myproc+1) ) THEN
    WRITE(6,*) 'mpi_1di_collect inconsistency:  ',k,kdata(myproc+1)
    CALL arpsstop("mpi_1di_collect",1)
  END IF

  DO k=1,nmap

    CALL inctag
    itag = gentag

    !
    !  Are we a sender?
    !
    IF (mpi_map(k,1) .NE. -1 ) THEN
      CALL mpsendi(tmps,isum,mpi_map(k,1),itag,ierror)
    ENDIF

    !
    !  Are we a receiver?
    !
    IF (mpi_map(k,2) .NE. -1) THEN
      CALL mprecvi(tmpr,isum,mpi_map(k,2),itag,ierror)
    ELSE
      CYCLE
    END IF
 
    l = 0
    DO ksta=1,nobmpi
      !
      !  Make sure we are the right processor.
      !
      IF (indexmpi(ksta) .NE. mpi_map(k,2)) CYCLE
      l = l + 1
      m(ksta) = tmpr(l)
    END DO
  END DO

  RETURN
END SUBROUTINE mpi_1di_collect

SUBROUTINE mpi_1dr_collect(m, nobmpi, indexmpi,                          &
  np, kdata, kdatamax, mpi_map, nmap, tmps, tmpr )

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Collect all the processor dependent calculations and merge them into a
! single array which we will broadcast.
!
!-----------------------------------------------------------------------
!
! AUTHOR: Kevin W. Thomas, CAPS
! January 10, 2006
!
! MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
! Force explicit declarations
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
! Subroutine arguments
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nobmpi         ! Actual number of stations
  INTEGER, INTENT(IN) :: indexmpi(nobmpi) ! Array saying which ob is on which
                                        ! processor.
  REAL, INTENT(INOUT) :: m(nobmpi)      ! Data to be collected and updated.
  INTEGER, INTENT(IN) :: np             ! Just "nprocs"
  INTEGER, INTENT(IN) :: kdata(np)      ! Number of obs owned by each processor
  INTEGER, INTENT(IN) :: kdatamax       ! Largest "kdata" value
  INTEGER, INTENT(IN) :: nmap           ! Number of entries in "mpi_map"
  INTEGER, INTENT(IN) :: mpi_map(nmap,2)! Mapping scheme

  REAL, INTENT(INOUT) :: tmps(kdatamax)
  REAL, INTENT(INOUT) :: tmpr(kdatamax)

!-----------------------------------------------------------------------
!
! Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER :: istat
  INTEGER :: isum
  INTEGER :: itag

  INTEGER :: ierror

  INTEGER :: i, j, k, l, ksta

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  isum = kdatamax

!
!  Save our data.
!

  k = 0
  DO j=1,nobmpi
    IF (indexmpi(j) == myproc) THEN
      k=k+1
      tmps(k) = m(j)
    END IF
  END DO

!
! Sanity.
!

  IF ( k .ne. kdata(myproc+1) ) THEN
    WRITE(6,*) 'mpi_1dr_collect inconsistency:  ',k,kdata(myproc+1)
    CALL arpsstop("mpi_1dr_collect",1)
  END IF

  DO k=1,nmap

    CALL inctag
    itag = gentag

!
!  Are we a sender?
!

    IF (mpi_map(k,1) .NE. -1 ) THEN
      CALL mpsendr(tmps,isum,mpi_map(k,1),itag,ierror)
    ENDIF

!
!  Are we a receiver?
!

    IF (mpi_map(k,2) .NE. -1) THEN
      CALL mprecvr(tmpr,isum,mpi_map(k,2),itag,ierror)
    ELSE
      CYCLE
    END IF
 
    l = 0
    DO ksta=1,nobmpi
!
!  Make sure we are the right processor.
!
      IF (indexmpi(ksta) .NE. mpi_map(k,2)) CYCLE
      l = l + 1
      m(ksta) = tmpr(l)
    END DO
  END DO

  RETURN
END SUBROUTINE mpi_1dr_collect

SUBROUTINE mpi_2dr_collect(q, nvar, mxmpi, nobmpi, indexmpi,             &
  np, kdata, kdatamax, mpi_map, nmap, tmps, tmpr )

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Collect all the processor dependent calculations and merge them into a
! single array which we will broadcast.
!
!-----------------------------------------------------------------------
!
! AUTHOR: Kevin W. Thomas, CAPS
! January 10, 2006
!
! MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
! Force explicit declarations
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
! Subroutine arguments
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nvar           ! Number of variables
  INTEGER, INTENT(IN) :: mxmpi          ! Max possible stations
  INTEGER, INTENT(IN) :: nobmpi         ! Actual number of stations
  INTEGER, INTENT(IN) :: indexmpi(nobmpi) ! Array saying which ob is on which
                                        ! processor.
  REAL, INTENT(INOUT) :: q(nvar,mxmpi)  ! Data to be collected and updated.
  INTEGER, INTENT(IN) :: np             ! Just "nprocs"
  INTEGER, INTENT(IN) :: kdata(np)      ! Number of obs owned by each processor
  INTEGER, INTENT(IN) :: kdatamax       ! Largest "kdata" value
  INTEGER, INTENT(IN) :: nmap           ! Number of entries in "mpi_map"
  INTEGER, INTENT(IN) :: mpi_map(nmap,2)! Mapping scheme

  REAL, INTENT(INOUT) :: tmps(nvar,kdatamax)
  REAL, INTENT(INOUT) :: tmpr(nvar,kdatamax)

!-----------------------------------------------------------------------
!
! Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER :: istat
  INTEGER :: isum
  INTEGER :: itag

  INTEGER :: ierror

  INTEGER :: i, j, k, l, ksta

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  isum = nvar * kdatamax

!
!  Save our data.
!

  k = 0
  DO j=1,nobmpi
    IF (indexmpi(j) == myproc) THEN
      k=k+1
      DO i=1,nvar
        tmps(i,k) = q(i,j)
      END DO
    END IF
  END DO

!
! Sanity.
!

  IF ( k .ne. kdata(myproc+1) ) THEN
    WRITE(6,*) 'mpi_2dr_collect inconsistency:  ',k,kdata(myproc+1)
    CALL arpsstop("mpi_2dr_collect",1)
  END IF

  DO k=1,nmap

    CALL inctag
    itag = gentag

!
!  Are we a sender?
!

    IF (mpi_map(k,1) .NE. -1 ) THEN
      IF (kdata(myproc+1) > 0) THEN
        CALL mpsendr(tmps,isum,mpi_map(k,1),itag,ierror)
      END IF
    ENDIF

!
!  Are we a receiver?
!

    IF (mpi_map(k,2) .NE. -1) THEN
      IF (kdata(mpi_map(k,2)+1) > 0)  THEN
        CALL mprecvr(tmpr,isum,mpi_map(k,2),itag,ierror)
      ELSE
        CYCLE
      END IF
    ELSE
      CYCLE
    END IF
 
    l = 0
    DO ksta=1,nobmpi
!
!  Make sure we are the right processor.
!
      IF (indexmpi(ksta) .NE. mpi_map(k,2)) CYCLE
      l = l + 1
      DO i=1,nvar
        q(i,ksta) = tmpr(i,l)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE mpi_2dr_collect

SUBROUTINE mpi_2dcr_collect(q, mxmpi, nvar, nobmpi, indexmpi )

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Collect all the processor dependent calculations and merge them into a
! single array which we will broadcast.
!
! Similar to "mpi_2dr_collect" except the subscripts are reversed.  Most of
! the arrays have the station index as the last subscript.  For "cloud
! soundings, the station index is the first subscript.
!
!-----------------------------------------------------------------------
!
! AUTHOR: Kevin W. Thomas, CAPS
! September 15, 2006
!
! MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
! Force explicit declarations
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
! Subroutine arguments
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nvar           ! Number of variables
  INTEGER, INTENT(IN) :: mxmpi          ! Max possible stations
  INTEGER, INTENT(IN) :: nobmpi         ! Actual number of stations
  INTEGER, INTENT(IN) :: indexmpi(nobmpi) ! Array saying which ob is on which
                                        ! processor.
  REAL, INTENT(INOUT) :: q(mxmpi,nvar)  ! Data to be collected and updated.

!-----------------------------------------------------------------------
!
! Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------

  REAL, ALLOCATABLE :: tmp(:,:)
  INTEGER :: istat
  INTEGER :: isum

  INTEGER :: ierror

  INTEGER :: i, k, ksta

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ALLOCATE(tmp(mxmpi,nvar),STAT=istat)
  CALL check_alloc_status(istat, "mpi_2dcr_collect:tmp")

  isum = nvar * mxmpi

!
! Collect and merge the data.
!

  IF (myproc > 0 ) THEN
    CALL mpsendr(q,isum,0,4000+myproc,ierror)
  ELSE
    DO k=1,nprocs-1
      CALL mprecvr(tmp,isum,k,4000+k,ierror)
!
!  We only have to process "nobmpi" of the obs, even the array has space for
!  "mxmpi" obs.
!
      DO ksta=1,nobmpi
!
!  Make sure we are the right processor.
!
        IF (indexmpi(ksta) .NE. k) CYCLE
        DO i=1,nvar
          q(ksta,i) = tmp(ksta,i)
        END DO
      END DO
    END DO
  END IF

  CALL mpupdater(q,isum)

  DEALLOCATE(tmp)

  RETURN
END SUBROUTINE mpi_2dcr_collect

SUBROUTINE mpi_3di_collect(q, nvar, nzmpi, mxmpi, nobmpi, indexmpi,    &
  np, kdata, kdatamax, mpi_map, nmap, tmps, tmpr)
!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Collect all the processor dependent calculations and merge them into a
! single array which we will broadcast.
!
!-----------------------------------------------------------------------
!
! AUTHOR: Kevin W. Thomas, CAPS
! January 10, 2006
!
! MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
! Force explicit declarations
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
! Subroutine arguments
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nvar           ! Number of variables
  INTEGER, INTENT(IN) :: nzmpi          ! Number of vertical levels
  INTEGER, INTENT(IN) :: mxmpi          ! Max possible stations
  INTEGER, INTENT(IN) :: nobmpi         ! Actual number of stations
  INTEGER, INTENT(IN) :: indexmpi(nobmpi) ! Array saying which ob is on which
                                        ! processor.
  INTEGER, INTENT(INOUT) :: q(nvar,nzmpi,mxmpi)! Data to be collected & updated.
  INTEGER, INTENT(IN) :: np             ! Just "nprocs"
  INTEGER, INTENT(IN) :: kdata(np)      ! Number of obs owned by each processor
  INTEGER, INTENT(IN) :: kdatamax       ! Largest "kdata" value
  INTEGER, INTENT(IN) :: nmap           ! Number of entries in "mpi_map"
  INTEGER, INTENT(IN) :: mpi_map(nmap,2)! Mapping scheme
  INTEGER, INTENT(INOUT) :: tmps(nvar,nzmpi,kdatamax)
  INTEGER, INTENT(INOUT) :: tmpr(nvar,nzmpi,kdatamax)

!-----------------------------------------------------------------------
!
! Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER :: istat
  INTEGER :: isum
  INTEGER :: itag

  INTEGER :: ierror

  INTEGER :: i, j, k, l, m, ksta

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  isum = nvar * nzmpi *kdatamax

!
!  Save our data.
!

  m = 0
  DO k=1,nobmpi
    IF (indexmpi(k) == myproc) THEN
      m=m+1
      DO j=1,nzmpi
        DO i=1,nvar
          tmps(i,j,m) = q(i,j,k)
        END DO
      END DO
    END IF
  END DO

!
!  Sanity.
!

  IF ( m .ne. kdata(myproc+1) ) THEN
    WRITE(6,*) 'mpi_3di_collect inconsistency:  ',m,kdata(myproc+1)
    CALL arpsstop("mpi_3di_collect",1)
  END IF

  DO l=1,nmap

    CALL inctag
    itag = gentag

!
!  Are we a sender?
!

    IF (mpi_map(l,1) .NE. -1 ) THEN
      IF (kdata(myproc+1) > 0) THEN
        CALL mpsendr(tmps,isum,mpi_map(l,1),itag,ierror)
      END IF
    ENDIF

!
!  Are we a receiver?
!

    IF (mpi_map(l,2) .NE. -1) THEN
      IF (kdata(mpi_map(l,2)+1) > 0) THEN
        CALL mprecvr(tmpr,isum,mpi_map(l,2),itag,ierror)
      ELSE
        CYCLE
      END IF
    ELSE
      CYCLE
    END IF
 
    m = 0
    DO ksta=1,nobmpi
!
!  Make sure we are the right processor.
!
      IF (indexmpi(ksta) .NE. mpi_map(l,2)) CYCLE
      m = m + 1
      DO j=1,nzmpi
        DO i=1,nvar
          q(i,j,ksta) = tmpr(i,j,m)
        END DO
      END DO
    END DO
  END DO

  RETURN

END SUBROUTINE mpi_3di_collect

SUBROUTINE mpi_3dr_collect(q, nvar, nzmpi, mxmpi, nobmpi, indexmpi,    &
  np, kdata, kdatamax, mpi_map, nmap, tmps, tmpr)
!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Collect all the processor dependent calculations and merge them into a
! single array which we will broadcast.
!
!-----------------------------------------------------------------------
!
! AUTHOR: Kevin W. Thomas, CAPS
! January 10, 2006
!
! MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
! Force explicit declarations
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
! Subroutine arguments
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nvar           ! Number of variables
  INTEGER, INTENT(IN) :: nzmpi          ! Number of vertical levels
  INTEGER, INTENT(IN) :: mxmpi          ! Max possible stations
  INTEGER, INTENT(IN) :: nobmpi         ! Actual number of stations
  INTEGER, INTENT(IN) :: indexmpi(nobmpi) ! Array saying which ob is on which
                                        ! processor.
  REAL, INTENT(INOUT) :: q(nvar,nzmpi,mxmpi)! Data to be collected and updated.
  INTEGER, INTENT(IN) :: np             ! Just "nprocs"
  INTEGER, INTENT(IN) :: kdata(np)      ! Number of obs owned by each processor
  INTEGER, INTENT(IN) :: kdatamax       ! Largest "kdata" value
  INTEGER, INTENT(IN) :: nmap           ! Number of entries in "mpi_map"
  INTEGER, INTENT(IN) :: mpi_map(nmap,2)! Mapping scheme
  REAL, INTENT(INOUT) :: tmps(nvar,nzmpi,kdatamax)
  REAL, INTENT(INOUT) :: tmpr(nvar,nzmpi,kdatamax)

!-----------------------------------------------------------------------
!
! Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER :: istat
  INTEGER :: isum
  INTEGER :: itag

  INTEGER :: ierror

  INTEGER :: i, j, k, l, m, ksta

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  isum = nvar * nzmpi * kdatamax

!
!  Save our data.
!

  m = 0
  DO k=1,nobmpi
    IF (indexmpi(k) == myproc) THEN
      m=m+1
      DO j=1,nzmpi
        DO i=1,nvar
          tmps(i,j,m) = q(i,j,k)
        END DO
      END DO
    END IF
  END DO

!
!  Sanity.
!

  IF ( m .ne. kdata(myproc+1) ) THEN
    WRITE(6,*) 'mpi_3dr_collect inconsistency:  ',m,kdata(myproc+1)
    CALL arpsstop("mpi_3dr_collect",1)
  END IF

  DO l=1,nmap

    CALL inctag
    itag = gentag

!
!  Are we a sender?
!

    IF (mpi_map(l,1) .NE. -1 ) THEN
      IF (kdata(myproc+1) > 0) THEN
        CALL mpsendr(tmps,isum,mpi_map(l,1),itag,ierror)
      END IF
    ENDIF

!
!  Are we a receiver?
!

    IF (mpi_map(l,2) .NE. -1) THEN
      IF (kdata(mpi_map(l,2)+1) > 0) THEN
        CALL mprecvr(tmpr,isum,mpi_map(l,2),itag,ierror)
      ELSE
        CYCLE
      END IF
    ELSE
      CYCLE
    END IF
 
    m = 0
    DO ksta=1,nobmpi
!
!  Make sure we are the right processor.
!
      IF (indexmpi(ksta) .NE. mpi_map(l,2)) CYCLE
      m = m + 1
      DO j=1,nzmpi
        DO i=1,nvar
          q(i,j,ksta) = tmpr(i,j,m)
        END DO
      END DO
    END DO
  END DO

  RETURN

END SUBROUTINE mpi_3dr_collect
