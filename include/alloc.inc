!##################################################################
!##################################################################
!######                                                      ######
!######                      ALLOC.INC                       ######
!######                                                      ######
!######                     Developed by                     ######
!######         Center for Analysis and Prediction of Storms ######
!######               University of Oklahoma                 ######
!######                                                      ######
!##################################################################
!##################################################################

!-----------------------------------------------------------------------
!
! PURPOSE:
! Provides interface definitions for array allocation and deallocation
! subroutines. This include file is required whenever those subroutines
! are called.
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

  interface

    SUBROUTINE alloc_real_1d_array(array,variable,nx)
      REAL, pointer :: array(:)
      CHARACTER(LEN=*), intent(in) :: variable
      INTEGER, intent(in) :: nx
    END SUBROUTINE alloc_real_1d_array
    
    SUBROUTINE alloc_real_2d_array(array,variable,nx,ny)
      REAL, pointer :: array(:,:)
      CHARACTER(LEN=*), intent(in) :: variable
      INTEGER, intent(in) :: nx,ny
    END SUBROUTINE alloc_real_2d_array
    
    SUBROUTINE alloc_real_3d_array(array,variable,nx,ny,nz)
      REAL, pointer :: array(:,:,:)
      CHARACTER(LEN=*), intent(in) :: variable
      INTEGER, intent(in) :: nx,ny,nz
    END SUBROUTINE alloc_real_3d_array
    
    SUBROUTINE alloc_real_4d_array(array,variable,nx,ny,nz,nt)
      REAL, pointer :: array(:,:,:,:)
      CHARACTER(LEN=*), intent(in) :: variable
      INTEGER, intent(in) :: nx,ny,nz,nt
    END SUBROUTINE alloc_real_4d_array
    
    SUBROUTINE alloc_int_1d_array(array,variable,nx)
      INTEGER, pointer :: array(:)
      CHARACTER(LEN=*), intent(in) :: variable
      INTEGER, intent(in) :: nx
    END SUBROUTINE alloc_int_1d_array
    
    SUBROUTINE alloc_int_2d_array(array,variable,nx,ny)
      INTEGER, pointer :: array(:,:)
      CHARACTER(LEN=*), intent(in) :: variable
      INTEGER, intent(in) :: nx,ny
    END SUBROUTINE alloc_int_2d_array
    
    SUBROUTINE alloc_int_3d_array(array,variable,nx,ny,nz)
      INTEGER, pointer :: array(:,:,:)
      CHARACTER(LEN=*), intent(in) :: variable
      INTEGER, intent(in) :: nx,ny,nz
    END SUBROUTINE alloc_int_3d_array
    
    SUBROUTINE alloc_int_4d_array(array,variable,nx,ny,nz,nt)
      INTEGER, pointer :: array(:,:,:,:)
      CHARACTER(LEN=*), intent(in) :: variable
      INTEGER, intent(in) :: nx,ny,nz,nt
    END SUBROUTINE alloc_int_4d_array
    
    SUBROUTINE alloc_logic_2d_array(array,variable,nx,ny)
      LOGICAL, pointer :: array(:,:)
      CHARACTER(LEN=*), intent(in) :: variable
      INTEGER, intent(in) :: nx,ny
    END SUBROUTINE alloc_logic_2d_array
    
    SUBROUTINE alloc_int2_1d_array(array,variable,nx)
      INTEGER(kind=2), pointer :: array(:)
      CHARACTER(LEN=*), intent(in) :: variable
      INTEGER, intent(in) :: nx
    END SUBROUTINE alloc_int2_1d_array
    
    SUBROUTINE alloc_int2_2d_array(array,variable,nx,ny)
      INTEGER(kind=2), pointer :: array(:,:)
      CHARACTER(LEN=*), intent(in) :: variable
      INTEGER, intent(in) :: nx,ny
    END SUBROUTINE alloc_int2_2d_array
    
    SUBROUTINE alloc_int2_3d_array(array,variable,nx,ny,nz)
      INTEGER(kind=2), pointer :: array(:,:,:)
      CHARACTER(LEN=*), intent(in) :: variable
      INTEGER, intent(in) :: nx,ny,nz
    END SUBROUTINE alloc_int2_3d_array
    
    SUBROUTINE alloc_char_1d_array(array,variable,nx)
      CHARACTER, pointer :: array(:)
      CHARACTER(LEN=*), intent(in) :: variable
      INTEGER, intent(in) :: nx
    END SUBROUTINE alloc_char_1d_array
    
    SUBROUTINE dealloc_real_1d_array(array,variable)
      REAL, pointer :: array(:)
      CHARACTER(LEN=*), intent(in) :: variable
    END SUBROUTINE dealloc_real_1d_array
    
    SUBROUTINE dealloc_real_2d_array(array,variable)
      REAL, pointer :: array(:,:)
      CHARACTER(LEN=*), intent(in) :: variable
    END SUBROUTINE dealloc_real_2d_array
    
    SUBROUTINE dealloc_real_3d_array(array,variable)
      REAL, pointer :: array(:,:,:)
      CHARACTER(LEN=*), intent(in) :: variable
    END SUBROUTINE dealloc_real_3d_array
    
    SUBROUTINE dealloc_real_4d_array(array,variable)
      REAL, pointer :: array(:,:,:,:)
      CHARACTER(LEN=*), intent(in) :: variable
    END SUBROUTINE dealloc_real_4d_array
    
    SUBROUTINE dealloc_int_1d_array(array,variable)
      INTEGER, pointer :: array(:)
      CHARACTER(LEN=*), intent(in) :: variable
    END SUBROUTINE dealloc_int_1d_array
    
    SUBROUTINE dealloc_int_2d_array(array,variable)
      INTEGER, pointer :: array(:,:)
      CHARACTER(LEN=*), intent(in) :: variable
    END SUBROUTINE dealloc_int_2d_array
    
    SUBROUTINE dealloc_int_3d_array(array,variable)
      INTEGER, pointer :: array(:,:,:)
      CHARACTER(LEN=*), intent(in) :: variable
    END SUBROUTINE dealloc_int_3d_array
    
    SUBROUTINE dealloc_int_4d_array(array,variable)
      INTEGER, pointer :: array(:,:,:,:)
      CHARACTER(LEN=*), intent(in) :: variable
    END SUBROUTINE dealloc_int_4d_array
    
    SUBROUTINE dealloc_logic_2d_array(array,variable)
      LOGICAL, pointer :: array(:,:)
      CHARACTER(LEN=*), intent(in) :: variable
    END SUBROUTINE dealloc_logic_2d_array
    
    SUBROUTINE dealloc_int2_1d_array(array,variable)
      INTEGER(kind=2), pointer :: array(:)
      CHARACTER(LEN=*), intent(in) :: variable
    END SUBROUTINE dealloc_int2_1d_array
    
    SUBROUTINE dealloc_int2_2d_array(array,variable)
      INTEGER(kind=2), pointer :: array(:,:)
      CHARACTER(LEN=*), intent(in) :: variable
    END SUBROUTINE dealloc_int2_2d_array
    
    SUBROUTINE dealloc_int2_3d_array(array,variable)
      INTEGER(kind=2), pointer :: array(:,:,:)
      CHARACTER(LEN=*), intent(in) :: variable
    END SUBROUTINE dealloc_int2_3d_array
    
    SUBROUTINE dealloc_char_1d_array(array,variable)
      CHARACTER, pointer :: array(:)
      CHARACTER(LEN=*), intent(in) :: variable
    END SUBROUTINE dealloc_char_1d_array

  END interface

  REAL :: current_memory_use
  REAL :: max_memory_use
  COMMON /memory_use/ current_memory_use, max_memory_use

