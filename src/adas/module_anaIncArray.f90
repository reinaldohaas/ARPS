  module module_anaIncArray
!
!
! 1.1 background state vector
!     -----------------------
!
  IMPLICIT NONE
  save
!
!-----------------------------------------------------------------------
!
!  Arrays for analysis increment updating (a type of nudging) output
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: uincr(:,:,:)
  REAL, ALLOCATABLE :: vincr(:,:,:)
  REAL, ALLOCATABLE :: wincr(:,:,:)
  REAL, ALLOCATABLE :: pincr(:,:,:)
  REAL, ALLOCATABLE :: ptincr(:,:,:)
  REAL, ALLOCATABLE :: qvincr(:,:,:)
  REAL, ALLOCATABLE :: qcincr(:,:,:)
  REAL, ALLOCATABLE :: qrincr(:,:,:)
  REAL, ALLOCATABLE :: qiincr(:,:,:)
  REAL, ALLOCATABLE :: qsincr(:,:,:)
  REAL, ALLOCATABLE :: qhincr(:,:,:)
!

  CONTAINS

  SUBROUTINE allocateanaIncArray(nxndg,nyndg,nzndg)

    INTEGER:: nxndg,nyndg,nzndg,istatus
!
!-----------------------------------------------------------------------
!
!  Allocate array for analysis increment updating
!
!-----------------------------------------------------------------------
!
    ALLOCATE(uincr(nxndg,nyndg,nzndg),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:uincr")
    ALLOCATE(vincr(nxndg,nyndg,nzndg),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:vincr")
    ALLOCATE(wincr(nxndg,nyndg,nzndg),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:wincr")
    ALLOCATE(pincr(nxndg,nyndg,nzndg),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:pincr")
    ALLOCATE(ptincr(nxndg,nyndg,nzndg),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:ptincr")
    ALLOCATE(qvincr(nxndg,nyndg,nzndg),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:qvincr")
    ALLOCATE(qcincr(nxndg,nyndg,nzndg),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:qcincr")
    ALLOCATE(qrincr(nxndg,nyndg,nzndg),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:qrincr")
    ALLOCATE(qiincr(nxndg,nyndg,nzndg),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:qiincr")
    ALLOCATE(qsincr(nxndg,nyndg,nzndg),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:qsincr")
    ALLOCATE(qhincr(nxndg,nyndg,nzndg),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:qhincr")

  END  SUBROUTINE allocateanaIncArray

  SUBROUTINE deallocateanaIncArray

    DEALLOCATE(uincr )
    DEALLOCATE(vincr )
    DEALLOCATE(wincr )
    DEALLOCATE(pincr )
    DEALLOCATE(ptincr )
    DEALLOCATE(qvincr )
    DEALLOCATE(qcincr )
    DEALLOCATE(qrincr )
    DEALLOCATE(qiincr )
    DEALLOCATE(qsincr )
    DEALLOCATE(qhincr )

  END SUBROUTINE deallocateanaIncArray

END module module_anaIncArray
