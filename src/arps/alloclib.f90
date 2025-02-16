!
!##################################################################
!##################################################################
!######                                                      ######
!######      SUBROUTINEs for array allocation/deallocation   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!-----------------------------------------------------------------------
!
! PURPOSE:
! Suboroutines to allocate and deallocate arrays using pointer.
! They also keep track of total memory usage. 
!-----------------------------------------------------------------------
!
! AUTHOR: Ming Xue
! 07/14/99
!
! MODIFICATION HISTORY:
!
! 2000/04/17 (Gene Bassett)
!   Converted to F90 fixed format.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Allocate real arrays
!-----------------------------------------------------------------------

SUBROUTINE alloc_real_1d_array(array,variable,nx)

  IMPLICIT NONE
  REAL, POINTER :: array(:)
  CHARACTER(len=*), INTENT(IN) :: variable
  INTEGER, INTENT(IN) :: nx
  INTEGER :: istatus

  ALLOCATE(array(nx),STAT=istatus)
  IF(istatus /= 0) CALL alloc_failed(istatus,variable)
  array = 0.

  CALL f_memory_use(nx)

  RETURN

END SUBROUTINE alloc_real_1d_array

SUBROUTINE alloc_real_2d_array(array,variable,nx,ny)

  IMPLICIT NONE
  REAL, POINTER :: array(:,:)
  CHARACTER(len=*), INTENT(IN) :: variable
  INTEGER, INTENT(IN) :: nx,ny
  INTEGER :: istatus

  ALLOCATE(array(nx,ny),STAT=istatus)
  IF(istatus /= 0) CALL alloc_failed(istatus,variable)
  array = 0.

  CALL f_memory_use(nx*ny)

  RETURN

END SUBROUTINE alloc_real_2d_array

SUBROUTINE alloc_real_3d_array(array,variable,nx,ny,nz)

  REAL, POINTER :: array(:,:,:)
  CHARACTER(len=*), INTENT(IN) :: variable
  INTEGER, INTENT(IN) :: nx,ny,nz
  INTEGER :: istatus

  ALLOCATE(array(nx,ny,nz),STAT=istatus)
  IF(istatus /= 0) CALL alloc_failed(istatus,variable)
  array = 0.

  CALL f_memory_use(nx*ny*nz)

  RETURN

END SUBROUTINE alloc_real_3d_array

SUBROUTINE alloc_real_4d_array(array,variable,nx,ny,nz,nt)

  IMPLICIT NONE
  REAL, POINTER :: array(:,:,:,:)
  CHARACTER(len=*), INTENT(IN) :: variable
  INTEGER, INTENT(IN) :: nx,ny,nz,nt
  INTEGER :: istatus

  ALLOCATE(array(nx,ny,nz,nt),STAT=istatus)
  IF(istatus /= 0) CALL alloc_failed(istatus,variable)
  array = 0.

  CALL f_memory_use(nx*ny*nz*nt)

  RETURN

END SUBROUTINE alloc_real_4d_array

!-----------------------------------------------------------------------
! Allocate integer arrays
!-----------------------------------------------------------------------

SUBROUTINE alloc_int_1d_array(array,variable,nx)

  IMPLICIT NONE
  INTEGER, POINTER :: array(:)
  CHARACTER(len=*), INTENT(IN) :: variable
  INTEGER, INTENT(IN) :: nx
  INTEGER :: istatus

  ALLOCATE(array(nx),STAT=istatus)
  IF(istatus /= 0) CALL alloc_failed(istatus,variable)
  array = 0

  CALL f_memory_use(nx)

  RETURN

END SUBROUTINE alloc_int_1d_array

SUBROUTINE alloc_int_2d_array(array,variable,nx,ny)

  IMPLICIT NONE
  INTEGER, POINTER :: array(:,:)
  CHARACTER(len=*), INTENT(IN) :: variable
  INTEGER, INTENT(IN) :: nx,ny
  INTEGER :: istatus

  ALLOCATE(array(nx,ny),STAT=istatus)
  IF(istatus /= 0) CALL alloc_failed(istatus,variable)
  array = 0

  CALL f_memory_use(nx*ny)

  RETURN

END SUBROUTINE alloc_int_2d_array

SUBROUTINE alloc_int_3d_array(array,variable,nx,ny,nz)

  INTEGER, POINTER :: array(:,:,:)
  CHARACTER(len=*), INTENT(IN) :: variable
  INTEGER, INTENT(IN) :: nx,ny,nz
  INTEGER :: istatus

  ALLOCATE(array(nx,ny,nz),STAT=istatus)
  IF(istatus /= 0) CALL alloc_failed(istatus,variable)
  array = 0

  CALL f_memory_use(nx*ny*nz)

  RETURN

END SUBROUTINE alloc_int_3d_array

SUBROUTINE alloc_int_4d_array(array,variable,nx,ny,nz,nt)

  IMPLICIT NONE
  INTEGER, POINTER :: array(:,:,:,:)
  CHARACTER(len=*), INTENT(IN) :: variable
  INTEGER, INTENT(IN) :: nx,ny,nz,nt
  INTEGER :: istatus

  IF (ASSOCIATED(array)) THEN
    WRITE(6,'(/1x,a,a,a,/1x,a/)')                                      &
    'Warning: Array ',variable,' had already been allocated.',         &
    'Allocating it again will distroy its existing content.'            
  ENDIF

  ALLOCATE(array(nx,ny,nz,nt),STAT=istatus)
  IF(istatus /= 0) CALL alloc_failed(istatus,variable)
  array = 0

  CALL f_memory_use(nx*ny*nz*nt)

  RETURN

END SUBROUTINE alloc_int_4d_array

SUBROUTINE alloc_logic_2d_array(array,variable,nx,ny)

  IMPLICIT NONE
  LOGICAL, POINTER :: array(:,:)
  CHARACTER(len=*), INTENT(IN) :: variable
  INTEGER, INTENT(IN) :: nx,ny
  INTEGER :: istatus

  ALLOCATE(array(nx,ny),STAT=istatus)
  IF(istatus /= 0) CALL alloc_failed(istatus,variable)
  array = .false.

  CALL f_memory_use(nx*ny)

  RETURN

END SUBROUTINE alloc_logic_2d_array

!-----------------------------------------------------------------------
! Allocate misc arrays
!-----------------------------------------------------------------------

SUBROUTINE alloc_int2_1d_array(array,variable,nx)

  INTEGER(kind=2), POINTER :: array(:)
  CHARACTER(len=*), INTENT(IN) :: variable
  INTEGER, INTENT(IN) :: nx
  INTEGER :: istatus

  ALLOCATE(array(nx),STAT=istatus)
  IF(istatus /= 0) CALL alloc_failed(istatus,variable)
  array = 0_2

  CALL f_memory_use(nx/2)

  RETURN

END SUBROUTINE alloc_int2_1d_array

SUBROUTINE alloc_int2_2d_array(array,variable,nx,ny)

  INTEGER(kind=2), POINTER :: array(:,:)
  CHARACTER(len=*), INTENT(IN) :: variable
  INTEGER, INTENT(IN) :: nx,ny
  INTEGER :: istatus

  ALLOCATE(array(nx,ny),STAT=istatus)
  IF(istatus /= 0) CALL alloc_failed(istatus,variable)
  array = 0_2

  CALL f_memory_use((nx*ny)/2)

  RETURN

END SUBROUTINE alloc_int2_2d_array

SUBROUTINE alloc_int2_3d_array(array,variable,nx,ny,nz)

  INTEGER(kind=2), POINTER :: array(:,:,:)
  CHARACTER(len=*), INTENT(IN) :: variable
  INTEGER, INTENT(IN) :: nx,ny,nz
  INTEGER :: istatus

  ALLOCATE(array(nx,ny,nz),STAT=istatus)
  IF(istatus /= 0) CALL alloc_failed(istatus,variable)
  array = 0_2

  CALL f_memory_use((nx*ny*nz)/2)

  RETURN

END SUBROUTINE alloc_int2_3d_array

SUBROUTINE alloc_char_1d_array(array,variable,nx)

  CHARACTER, POINTER :: array(:)
  CHARACTER(len=*), INTENT(IN) :: variable
  INTEGER, INTENT(IN) :: nx
  INTEGER :: istatus

  ALLOCATE(array(nx),STAT=istatus)
  IF(istatus /= 0) CALL alloc_failed(istatus,variable)

  CALL f_memory_use(nx)

  RETURN

END SUBROUTINE alloc_char_1d_array

!-----------------------------------------------------------------------
! Deallocate real arrays
!-----------------------------------------------------------------------

SUBROUTINE dealloc_real_1d_array(array,variable)

  IMPLICIT NONE
  REAL, POINTER :: array(:)
  CHARACTER(len=*), INTENT(IN) :: variable
  INTEGER :: istatus, a_size

  IF (ASSOCIATED(array)) THEN
    a_size = size(array)
    DEALLOCATE(array,STAT=istatus)
    IF(istatus /= 0) THEN
      CALL dealloc_failed(istatus,variable)
    ELSE
      CALL f_memory_use(-a_size)
    ENDIF
  ELSE
    WRITE(6,'(/1x,a,a,a,/1x,a/)')                                      &
    'Pointer array ', variable,                                        &
    ' is not currently allocated, therefore',                          &
    'it could not be deallocated. Program will continue.'
  ENDIF

  RETURN

END SUBROUTINE dealloc_real_1d_array

SUBROUTINE dealloc_real_2d_array(array,variable)

  IMPLICIT NONE
  REAL, POINTER :: array(:,:)
  CHARACTER(len=*), INTENT(IN) :: variable
  INTEGER :: istatus, a_size

  IF (ASSOCIATED(array)) THEN
    a_size = size(array)
    DEALLOCATE(array,STAT=istatus)
    IF(istatus /= 0) THEN
      CALL dealloc_failed(istatus,variable)
    ELSE
      CALL f_memory_use(-a_size)
    ENDIF
  ELSE
    WRITE(6,'(/1x,a,a,a,/1x,a/)')                                      &
    'Pointer array ', variable,                                        &
    ' is not currently allocated, therefore',                          &
    'it could not be deallocated. Program will continue.'
  ENDIF

  RETURN

END SUBROUTINE dealloc_real_2d_array

SUBROUTINE dealloc_real_3d_array(array,variable)

  REAL, POINTER :: array(:,:,:)
  CHARACTER(len=*), INTENT(IN) :: variable
  INTEGER :: istatus, a_size

  IF (ASSOCIATED(array)) THEN
    a_size = size(array)
    DEALLOCATE(array,STAT=istatus)
    IF(istatus /= 0) THEN
      CALL dealloc_failed(istatus,variable)
    ELSE
      CALL f_memory_use(-a_size)
    ENDIF
  ELSE
    WRITE(6,'(/1x,a,a,a,/1x,a/)')                                       &
    'Pointer array ', variable,                                         &
    ' is not currently allocated, therefore',                           &
    'it could not be deallocated. Program will continue.'
  ENDIF

  RETURN

END SUBROUTINE dealloc_real_3d_array

SUBROUTINE dealloc_real_4d_array(array,variable)

  IMPLICIT NONE
  REAL, POINTER :: array(:,:,:,:)
  CHARACTER(len=*), INTENT(IN) :: variable
  INTEGER :: istatus, a_size

  IF (ASSOCIATED(array)) THEN
    a_size = size(array)
    DEALLOCATE(array,STAT=istatus)
    IF(istatus /= 0) THEN
      CALL dealloc_failed(istatus,variable)
    ELSE
      CALL f_memory_use(-a_size)
    ENDIF
  ELSE
    WRITE(6,'(/1x,a,a,a,/1x,a/)')                                      &
    'Pointer array ', variable,                                        &
    ' is not currently allocated, therefore',                          &
    'it could not be deallocated. Program will continue.'
  ENDIF

  RETURN

END SUBROUTINE dealloc_real_4d_array

!-----------------------------------------------------------------------
! Deallocate integer arrays
!-----------------------------------------------------------------------

SUBROUTINE dealloc_int_1d_array(array,variable)

  IMPLICIT NONE
  INTEGER, POINTER :: array(:)
  CHARACTER(len=*), INTENT(IN) :: variable
  INTEGER :: istatus, a_size

  IF (ASSOCIATED(array)) THEN
    a_size = size(array)
    DEALLOCATE(array,STAT=istatus)
    IF(istatus /= 0) THEN
      CALL dealloc_failed(istatus,variable)
    ELSE
      CALL f_memory_use(-a_size)
    ENDIF
  ELSE
    WRITE(6,'(/1x,a,a,a,/1x,a/)')                                      &
    'Pointer array ', variable,                                        &
    ' is not currently allocated, therefore',                          &
    'it could not be deallocated. Program will continue.'
  ENDIF

  RETURN

END SUBROUTINE dealloc_int_1d_array

SUBROUTINE dealloc_int_2d_array(array,variable)

  IMPLICIT NONE
  INTEGER, POINTER :: array(:,:)
  CHARACTER(len=*), INTENT(IN) :: variable
  INTEGER :: istatus, a_size

  IF (ASSOCIATED(array)) THEN
    a_size = size(array)
    DEALLOCATE(array,STAT=istatus)
    IF(istatus /= 0) THEN
      CALL dealloc_failed(istatus,variable)
    ELSE
      CALL f_memory_use(-a_size)
    ENDIF
  ELSE
    WRITE(6,'(/1x,a,a,a,/1x,a/)')                                      &
    'Pointer array ', variable,                                        &
    ' is not currently allocated, therefore',                          &
    'it could not be deallocated. Program will continue.'
  ENDIF

  RETURN

END SUBROUTINE dealloc_int_2d_array

SUBROUTINE dealloc_int_3d_array(array,variable)

  INTEGER, POINTER :: array(:,:,:)
  CHARACTER(len=*), INTENT(IN) :: variable
  INTEGER :: istatus, a_size

  IF (ASSOCIATED(array)) THEN
    a_size = size(array)
    DEALLOCATE(array,STAT=istatus)
    IF(istatus /= 0) THEN
      CALL dealloc_failed(istatus,variable)
    ELSE
      CALL f_memory_use(-a_size)
    ENDIF
  ELSE
    WRITE(6,'(/1x,a,a,a,/1x,a/)')                                      &
    'Pointer array ', variable,                                        &
    ' is not currently allocated, therefore',                          &
    'it could not be deallocated. Program will continue.'
  ENDIF

  RETURN

END SUBROUTINE dealloc_int_3d_array

SUBROUTINE dealloc_int_4d_array(array,variable)

  IMPLICIT NONE
  INTEGER, POINTER :: array(:,:,:,:)
  CHARACTER(len=*), INTENT(IN) :: variable
  INTEGER :: istatus, a_size

  IF (ASSOCIATED(array)) THEN
    a_size = size(array)
    DEALLOCATE(array,STAT=istatus)
    IF(istatus /= 0) THEN
      CALL dealloc_failed(istatus,variable)
    ELSE
      CALL f_memory_use(-a_size)
    ENDIF
  ELSE
    WRITE(6,'(/1x,a,a,a,/1x,a/)')                                      &
    'Pointer array ', variable,                                        &
    ' is not currently allocated, therefore',                          &
    'it could not be deallocated. Program will continue.'
  ENDIF

  RETURN

END SUBROUTINE dealloc_int_4d_array

!-----------------------------------------------------------------------
! Deallocate misc arrays
!-----------------------------------------------------------------------

SUBROUTINE dealloc_logic_2d_array(array,variable)

  IMPLICIT NONE
  LOGICAL, POINTER :: array(:,:)
  CHARACTER(len=*), INTENT(IN) :: variable
  INTEGER :: istatus, a_size

  IF (ASSOCIATED(array)) THEN
    a_size = size(array)
    DEALLOCATE(array,STAT=istatus)
    IF(istatus /= 0) THEN
      CALL dealloc_failed(istatus,variable)
    ELSE
      CALL f_memory_use(-a_size)
    ENDIF
  ELSE
    WRITE(6,'(/1x,a,a,a,/1x,a/)')                                      &
    'Pointer array ', variable,                                        &
    ' is not currently allocated, therefore',                          &
    'it could not be deallocated. Program will continue.'
  ENDIF

  RETURN

END SUBROUTINE dealloc_logic_2d_array

SUBROUTINE dealloc_int2_1d_array(array,variable)

  INTEGER(kind=2), POINTER :: array(:)
  CHARACTER(len=*), INTENT(IN) :: variable
  INTEGER :: istatus, a_size

  IF (ASSOCIATED(array)) THEN
    a_size = size(array)
    DEALLOCATE(array,STAT=istatus)
    IF(istatus /= 0) THEN
      CALL dealloc_failed(istatus,variable)
    ELSE
      CALL f_memory_use(-a_size/2)
    ENDIF
  ELSE
    WRITE(6,'(/1x,a,a,a,/1x,a/)')                                      &
    'Pointer array ', variable,                                        &
    ' is not currently allocated, therefore',                          &
    'it could not be deallocated. Program will continue.'
  ENDIF

  RETURN

END SUBROUTINE dealloc_int2_1d_array

SUBROUTINE dealloc_int2_2d_array(array,variable)

  INTEGER(kind=2), POINTER :: array(:,:)
  CHARACTER(len=*), INTENT(IN) :: variable
  INTEGER :: istatus, a_size

  IF (ASSOCIATED(array)) THEN
    a_size = size(array)
    DEALLOCATE(array,STAT=istatus)
    IF(istatus /= 0) THEN
      CALL dealloc_failed(istatus,variable)
    ELSE
      CALL f_memory_use(-a_size/2)
    ENDIF
  ELSE
    WRITE(6,'(/1x,a,a,a,/1x,a/)')                                      &
    'Pointer array ', variable,                                        &
    ' is not currently allocated, therefore',                          &
    'it could not be deallocated. Program will continue.'
  ENDIF

  RETURN

END SUBROUTINE dealloc_int2_2d_array

SUBROUTINE dealloc_int2_3d_array(array,variable)

  INTEGER(kind=2), POINTER :: array(:,:,:)
  CHARACTER(len=*), INTENT(IN) :: variable
  INTEGER :: istatus, a_size

  IF (ASSOCIATED(array)) THEN
    a_size = size(array)
    DEALLOCATE(array,STAT=istatus)
    IF(istatus /= 0) THEN
      CALL dealloc_failed(istatus,variable)
    ELSE
      CALL f_memory_use(-a_size/2)
    ENDIF
  ELSE
    WRITE(6,'(/1x,a,a,a,/1x,a/)')                                      &
    'Pointer array ', variable,                                        &
    ' is not currently allocated, therefore',                          &
    'it could not be deallocated. Program will continue.'
  ENDIF

  RETURN

END SUBROUTINE dealloc_int2_3d_array

SUBROUTINE dealloc_char_1d_array(array,variable)

  CHARACTER, POINTER :: array(:)
  CHARACTER(len=*), INTENT(IN) :: variable
  INTEGER :: istatus, a_size

  IF (ASSOCIATED(array)) THEN
    a_size = size(array)
    DEALLOCATE(array,STAT=istatus)
    IF(istatus /= 0) THEN
      CALL dealloc_failed(istatus,variable)
    ELSE
      CALL f_memory_use(-a_size)
    ENDIF
  ELSE
    WRITE(6,'(/1x,a,a,a,/1x,a/)')                                      &
    'Pointer array ', variable,                                        &
    ' is not currently allocated, therefore',                          &
    'it could not be deallocated. Program will continue.'
  ENDIF

  RETURN

END SUBROUTINE dealloc_char_1d_array

!-----------------------------------------------------------------------
! Array allocation status reporting
!-----------------------------------------------------------------------

SUBROUTINE alloc_failed(istatus,variable)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: istatus
  CHARACTER(len=*), INTENT(IN) :: variable

  WRITE(6,'(/1x,a,a,/1x,a,i2,a/)')                                     &
  'Program failed when allocating memory space for array ',            &
  variable,                                                            &
  'Program stopped. Status of allocation is ', istatus,'.'
  CALL arpsstop('arpsstop called from alloc_failed',1)

END SUBROUTINE alloc_failed

SUBROUTINE dealloc_failed(istatus,variable)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: istatus
  CHARACTER(len=*), INTENT(IN) :: variable

  WRITE(6,'(/1x,a,a,/1x,a,i2,a/)')                                     &                 
  'Program failed when deallocating memory space for array ',          &
  variable,                                                            &
  'Program will continue. Status of allocation is ', istatus,'.'

  RETURN

END SUBROUTINE dealloc_failed

SUBROUTINE Alloc_status_accounting(istatus,a_size,variable)

  CHARACTER(len=*), INTENT(IN) :: variable
  INTEGER :: istatus, a_size

  IF(istatus /= 0) THEN
    IF(a_size.ge.0) then
      CALL alloc_failed(istatus,variable)
    ELSE
      CALL dealloc_failed(istatus,variable)
    ENDIF
  ELSE
    CALL f_memory_use(a_size)
  ENDIF

  RETURN

END SUBROUTINE Alloc_status_accounting 
!
!##################################################################
!##################################################################
!######                                                      ######
!######                MODULE memory_accounting              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
MODULE memory_accounting

!-----------------------------------------------------------------------
!
! PURPOSE:
! Module for saving and passing memory usage information.
!
!-----------------------------------------------------------------------
!
! AUTHOR: Ming Xue
! 07/14/99
!
! MODIFICATION HISTORY:
!
! 2000/04/17 (Gene Bassett)
!   Converted to F90 fixed format.
!-----------------------------------------------------------------------

  IMPLICIT NONE
  SAVE 
  REAL :: current_memory_use = 0.0 
  REAL :: max_memory_use = 0.0 

END MODULE memory_accounting

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE f_memory_use               ######
!######                                                      ######
!######                     Developed by                     ######
!######	    Center for Analysis and Prediction of Storms     ######
!######	             University of Oklahoma                  ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE f_memory_use( n_word )

!-----------------------------------------------------------------------
!
! PURPOSE:
! Subroutine for performing memory accounting 
!
!-----------------------------------------------------------------------
!
! AUTHOR: Ming Xue
! 7/14/99
!
! MODIFICATION HISTORY:
!
! 2000/04/17 (Gene Bassett)
!   Converted to F90 fixed format.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INCLUDE 'alloc.inc'

  INTEGER, INTENT(IN) :: n_word ! Number of words to be allocated (positive)
                                ! or deallocated (nagative).

  current_memory_use = current_memory_use + n_word
  max_memory_use   = max(max_memory_use,current_memory_use)

END SUBROUTINE f_memory_use

!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE check_alloc_status            ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE check_alloc_status(istatus, message )

!-----------------------------------------------------------------------
!
! PURPOSE: Check status of array allocation.
!
!-----------------------------------------------------------------------
!
! AUTHOR: Richard Carpenter, 2001/12/07
!
! MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER,       INTENT(IN) :: istatus
  CHARACTER*(*), INTENT(IN) :: message

  IF (istatus /= 0) THEN
    WRITE(*,'(1x,/,A,I2,3A)') 'ERROR: check_alloc_status, status= ',    &
                              istatus, ', [', TRIM(message), '].'
    CALL arpsstop('FATAL: Unable to allocate array. Program ends',1)
  END IF

  RETURN
END SUBROUTINE check_alloc_status
