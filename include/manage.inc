!
  INTEGER :: idimf  ! the maximum number of free blocks in lfree
!   parameter( idimf = 50 )
!
!  lfree(n,1) is the starting index of free block number n
!  lfree(n,2) is the length free block number n
!
  INTEGER :: lfree(50,2)

  INTEGER :: lenf   ! The acutal number of free blocks in lfree

  COMMON / manage / lfree,lenf,idimf
  LOGICAL :: sprint
  DATA    sprint/.false./
