!########################################################################
!########################################################################
!#########                                                      #########
!#########                 MODULE hdf_constants                 #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################


MODULE hdf_constants

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Declares functions and settings used with HDF I/O routines.
!
!----------------------------------------------------------------------- 
!
! AUTHOR:  Rich Carpenter, Time Unknown.
!
! MODIFICATION HISTORY:
!
!   Eric Kemp, January 2002.  Documented code, and moved it to 
!   arps5.0.0Beta5 for use by arpscvtobs.
!
!----------------------------------------------------------------------- 
!
! Force explicit declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!----------------------------------------------------------------------- 
!
! Include statements.
!
!-----------------------------------------------------------------------

  INCLUDE 'hdf.f90'		! HDF constants
  !INCLUDE 'dffunc.f90'		! HDF function declarations

!----------------------------------------------------------------------- 
!
! Declare HDF functions.
!
!-----------------------------------------------------------------------

  INTEGER, EXTERNAL :: sfstart, sfcreate, sfwdata, sfendacc, &
      sfend, sfselect, sfsnatt, sfscatt, sfscompress, sffinfo, sffattr, &
      sfgainfo, sfginfo, sfrattr, sfrcatt, sfrdata, sfwcdata, sfrcdata, &
      sfn2index, sfdimid, sfsdmname

!-----------------------------------------------------------------------
!
! Specify GZIP (DEFLATE) compression. comp_prm(1) is deflate_level.
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: comp_type=COMP_CODE_DEFLATE, comp_prm(1)=6

END MODULE hdf_constants
