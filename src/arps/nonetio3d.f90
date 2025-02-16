!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE NETREAD                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE netread
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  A dummy routine to be called instead of the truly functional routine
!  netread when word packing routines are not available on a given system.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  WRITE (6,'(1x,/,a,/,10x,2a,/,/,2a,/,10x,a,/)')                        &
      ' WARNING: NetCDF library was not linked.','Please add ',         &
      '"-io net" option for makearps to link NetCDF 3.0 libraries.',    &
      ' NOTE:    NetCDF 3.0 library has a conflict with HDF 4 ',        &
      'libraries.','So NetCDF file and HDF 4 file cannot work mixed.'
  CALL arpsstop('arpsstop called from NETREAD not supported.',1)

END SUBROUTINE netread
!
!
SUBROUTINE netreadsplit
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  A dummy routine to be called instead of the truly functional routine
!  netread when word packing routines are not available on a given system.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  WRITE (6,'(1x,/,a,/,10x,2a,/,/,2a,/,10x,a,/)')                        &
      ' WARNING: NetCDF library was not linked.','Please add ',         &
      '"-io net" option for makearps to link NetCDF 3.0 libraries.',    &
      ' NOTE:    NetCDF 3.0 library has a conflict with HDF 4 ',        &
      'libraries.','So NetCDF file and HDF 4 file cannot work mixed.'
  CALL arpsstop('arpsstop called from NETREAD not supported.',1)

END SUBROUTINE netreadsplit
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE NETDUMP                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE netdump
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    A dummy routine to be called instead of the functional routine
!    netdump when packing is not avaiable.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  WRITE (6,'(1x,/,a,/,10x,2a,/,/,2a,/,10x,a,/)')                        &
      ' WARNING: NetCDF library was not linked.','Please add ',         &
      '"-io net" option for makearps to link NetCDF 3.0 libraries.',    &
      ' NOTE:    NetCDF 3.0 library has a conflict with HDF 4 ',        &
      'libraries.','So NetCDF file and HDF 4 file cannot work mixed.'

  CALL arpsstop('arpsstop called from NETDUMP not supported.',1)

END SUBROUTINE netdump
!
!
!
SUBROUTINE netjoindump
  IMPLICIT NONE

  WRITE (6,'(1x,/,a,/,10x,2a,/,/,2a,/,10x,a,/)')                        &
      ' WARNING: NetCDF library was not linked.','Please add ',         &
      '"-io net" option for makearps to link NetCDF 3.0 libraries.',    &
      ' NOTE:    NetCDF 3.0 library has a conflict with HDF 4 ',        &
      'libraries.','So NetCDF file and HDF 4 file cannot work mixed.'

  CALL arpsstop('arpsstop called from NETDUMP not supported.',1)

END SUBROUTINE netjoindump

SUBROUTINE netsplitdump
  IMPLICIT NONE

  WRITE (6,'(1x,/,a,/,10x,2a,/,/,2a,/,10x,a,/)')                        &
      ' WARNING: NetCDF library was not linked.','Please add ',         &
      '"-io net" option for makearps to link NetCDF 3.0 libraries.',    &
      ' NOTE:    NetCDF 3.0 library has a conflict with HDF 4 ',        &
      'libraries.','So NetCDF file and HDF 4 file cannot work mixed.'

  CALL arpsstop('arpsstop called from NETDUMP not supported.',1)

END SUBROUTINE netsplitdump

!
!
!
SUBROUTINE netreaddims

  WRITE (6,'(1x,/,a,/,10x,2a,/,/,2a,/,10x,a,/)')                        &
      ' WARNING: NetCDF library was not linked.','Please add ',         &
      '"-io net" option for makearps to link NetCDF 3.0 libraries.',    &
      ' NOTE:    NetCDF 3.0 library has a conflict with HDF 4 ',        &
      'libraries.','So NetCDF file and HDF 4 file cannot work mixed.'
 
  CALL arpsstop('arpsstop called from netreaddims, NetCDF not supported.',1)

END SUBROUTINE netreaddims

SUBROUTINE netreadatti
  RETURN
END SUBROUTINE netreadatti

SUBROUTINE netreadattr
  RETURN
END SUBROUTINE netreadattr

SUBROUTINE netreadattstr
  RETURN
END SUBROUTINE netreadattstr

SUBROUTINE net_define_trn

  WRITE (6,'(1x,/,a,/,10x,2a,/,/,2a,/,10x,a,/)')                        &
      ' WARNING: NetCDF library was not linked.','Please add ',         &
      '"-io net" option for makearps to link NetCDF 3.0 libraries.',    &
      ' NOTE:    NetCDF 3.0 library has a conflict with HDF 4 ',        &
      'libraries.','So NetCDF file and HDF 4 file cannot work mixed.'

  CALL arpsstop('arpsstop called from net_define_trn not supported.',1)

  RETURN
END SUBROUTINE net_define_trn

SUBROUTINE net_get_trn

  WRITE (6,'(1x,/,a,/,10x,2a,/,/,2a,/,10x,a,/)')                        &
      ' WARNING: NetCDF library was not linked.','Please add ',         &
      '"-io net" option for makearps to link NetCDF 3.0 libraries.',    &
      ' NOTE:    NetCDF 3.0 library has a conflict with HDF 4 ',        &
      'libraries.','So NetCDF file and HDF 4 file cannot work mixed.'

  CALL arpsstop('arpsstop called from net_get_trn not supported.',1)

  RETURN
END SUBROUTINE

SUBROUTINE net_define_sfc

  WRITE (6,'(1x,/,a,/,10x,2a,/,/,2a,/,10x,a,/)')                        &
      ' WARNING: NetCDF library was not linked.','Please add ',         &
      '"-io net" option for makearps to link NetCDF 3.0 libraries.',    &
      ' NOTE:    NetCDF 3.0 library has a conflict with HDF 4 ',        &
      'libraries.','So NetCDF file and HDF 4 file cannot work mixed.'

  CALL arpsstop('arpsstop called from net_define_sfc not supported.',1)

  RETURN
END SUBROUTINE net_define_sfc

SUBROUTINE net_get_sfc

  WRITE (6,'(1x,/,a,/,10x,2a,/,/,2a,/,10x,a,/)')                        &
      ' WARNING: NetCDF library was not linked.','Please add ',         &
      '"-io net" option for makearps to link NetCDF 3.0 libraries.',    &
      ' NOTE:    NetCDF 3.0 library has a conflict with HDF 4 ',        &
      'libraries.','So NetCDF file and HDF 4 file cannot work mixed.'

  CALL arpsstop('arpsstop called from net_get_sfc not supported.',1)

  RETURN
END SUBROUTINE net_get_sfc

SUBROUTINE net_define_soil

  WRITE (6,'(1x,/,a,/,10x,2a,/,/,2a,/,10x,a,/)')                        &
      ' WARNING: NetCDF library was not linked.','Please add ',         &
      '"-io net" option for makearps to link NetCDF 3.0 libraries.',    &
      ' NOTE:    NetCDF 3.0 library has a conflict with HDF 4 ',        &
      'libraries.','So NetCDF file and HDF 4 file cannot work mixed.'

  CALL arpsstop('arpsstop called from net_define_soil not supported.',1)

  RETURN
END SUBROUTINE net_define_soil

SUBROUTINE net_get_soil

  WRITE (6,'(1x,/,a,/,10x,2a,/,/,2a,/,10x,a,/)')                        &
      ' WARNING: NetCDF library was not linked.','Please add ',         &
      '"-io net" option for makearps to link NetCDF 3.0 libraries.',    &
      ' NOTE:    NetCDF 3.0 library has a conflict with HDF 4 ',        &
      'libraries.','So NetCDF file and HDF 4 file cannot work mixed.'

  CALL arpsstop('arpsstop called from net_get_soil not supported.',1)

  RETURN
END SUBROUTINE net_get_soil

SUBROUTINE net_define_exbc

  WRITE (6,'(1x,/,a,/,10x,2a,/,/,2a,/,10x,a,/)')                        &
      ' WARNING: NetCDF library was not linked.','Please add ',         &
      '"-io net" option for makearps to link NetCDF 3.0 libraries.',    &
      ' NOTE:    NetCDF 3.0 library has a conflict with HDF 4 ',        &
      'libraries.','So NetCDF file and HDF 4 file cannot work mixed.'

  CALL arpsstop('arpsstop called from net_define_soil not supported.',1)

  RETURN
END SUBROUTINE net_define_exbc

SUBROUTINE net_get_exbc

  WRITE (6,'(1x,/,a,/,10x,2a,/,/,2a,/,10x,a,/)')                        &
      ' WARNING: NetCDF library was not linked.','Please add ',         &
      '"-io net" option for makearps to link NetCDF 3.0 libraries.',    &
      ' NOTE:    NetCDF 3.0 library has a conflict with HDF 4 ',        &
      'libraries.','So NetCDF file and HDF 4 file cannot work mixed.'

  CALL arpsstop('arpsstop called from net_get_soil not supported.',1)

  RETURN
END SUBROUTINE net_get_exbc

SUBROUTINE netopen

  WRITE (6,'(1x,/,a,/,10x,2a,/,/,2a,/,10x,a,/)')                        &
      ' WARNING: NetCDF library was not linked.','Please add ',         &
      '"-io net" option for makearps to link NetCDF 3.0 libraries.',    &
      ' NOTE:    NetCDF 3.0 library has a conflict with HDF 4 ',        &
      'libraries.','So NetCDF file and HDF 4 file cannot work mixed.'
  CALL arpsstop('arpsstop called from NETOPEN not supported.',1)

  RETURN
END SUBROUTINE netopen

SUBROUTINE netclose
  RETURN
END SUBROUTINE netclose

SUBROUTINE netwrt2d
  RETURN
END SUBROUTINE netwrt2d

SUBROUTINE netread2d
  RETURN
END SUBROUTINE netread2d

SUBROUTINE netwrt2di
  RETURN
END SUBROUTINE netwrt2di

SUBROUTINE netread2di
  RETURN
END SUBROUTINE netread2di

SUBROUTINE netwrt3d
  RETURN
END SUBROUTINE netwrt3d

SUBROUTINE netread3d
  RETURN
END SUBROUTINE netread3d

SUBROUTINE netwrt3di
  RETURN
END SUBROUTINE netwrt3di

SUBROUTINE netread3di
  RETURN
END SUBROUTINE netread3di

SUBROUTINE netwrt4d
  RETURN
END SUBROUTINE netwrt4d

SUBROUTINE netread4d
  RETURN
END SUBROUTINE netread4d

SUBROUTINE net_define_onevar
  RETURN
END SUBROUTINE net_define_onevar

SUBROUTINE net_get_onevar
  RETURN
END SUBROUTINE net_get_onevar

SUBROUTINE net_get_unlimit_size( filename, no_times, istatus)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Get the size of unlimitted dimension in the file
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!  Yunheng Wang (04/20/2010)
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN)  :: filename
  INTEGER,          INTENT(OUT) :: no_times
  INTEGER,          INTENT(OUT) :: istatus

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus   = 0
  no_times  = 0

  RETURN
END SUBROUTINE net_get_unlimit_size
