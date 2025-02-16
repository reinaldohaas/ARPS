!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE open_ncd_file               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE open_ncd_file(filename,ncidout)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!    Open a WRF file and return NetCDF file handler.
!    NOTE: it is required to call close_wrf_file explicitly to close
!          the opened file in your calling program.
!
!------------------------------------------------------------------

  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN)  :: filename
  INTEGER,          INTENT(OUT) :: ncidout

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Begining of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


  WRITE(6,'(1x,a)') 'ERROR: NetCDF library is not linked.'
  CALL arpsstop('Please link netCDF libary and try again.',1)

  RETURN

END SUBROUTINE open_ncd_file

SUBROUTINE close_ncd_file
  RETURN
END SUBROUTINE close_ncd_file

SUBROUTINE get_ncd_dom_ti_integer
  RETURN
END SUBROUTINE get_ncd_dom_ti_integer

SUBROUTINE get_ncd_dom_ti_real
  RETURN
END SUBROUTINE get_ncd_dom_ti_real

SUBROUTINE get_ncd_dom_ti_char
  RETURN
END SUBROUTINE get_ncd_dom_ti_char

SUBROUTINE get_ncd_1d
  RETURN
END SUBROUTINE get_ncd_1d

SUBROUTINE get_ncd_scalar
  RETURN
END SUBROUTINE get_ncd_scalar

SUBROUTINE get_ncd_2d
  RETURN
END SUBROUTINE get_ncd_2d

SUBROUTINE get_ncd_2di
  RETURN
END SUBROUTINE get_ncd_2di

SUBROUTINE get_ncd_3d
  RETURN
END SUBROUTINE get_ncd_3d

SUBROUTINE get_ncd_dimension
  RETURN
END SUBROUTINE get_ncd_dimension

SUBROUTINE nf_handle_error
  RETURN
END SUBROUTINE nf_handle_error
