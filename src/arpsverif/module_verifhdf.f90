

!                   ########################################
!                   ########################################
!                   ####                                ####
!                   ####          hdf_constants         ####
!                   ####                                ####
!                   ########################################
!                   ########################################

MODULE hdf_constants

  IMPLICIT NONE
  INCLUDE 'hdf.f90'		! HDF constants
  !INCLUDE 'dffunc.f90'		! HDF function declarations

  INTEGER, EXTERNAL :: sfstart, sfcreate, sfwdata, sfendacc, &
      sfend, sfselect, sfsnatt, sfscatt, sfscompress, sffinfo, sffattr, &
      sfgainfo, sfginfo, sfrattr, sfrcatt, sfrdata, sfwcdata, sfrcdata, &
      sfn2index

  ! specify GZIP (DEFLATE) compression. comp_prm(1) is deflate_level.
  INTEGER, PARAMETER :: comp_type=COMP_CODE_DEFLATE, comp_prm(1)=6

END MODULE hdf_constants


