SUBROUTINE hdfread
  WRITE (6,'(/a/a/)')                                                   &
      'HDF library was not linked. No HDF data read created.',          &
      'Make sure you include "-io hdf" or "-io hdfnet" option when doing makearps.'
  CALL arpsstop('arpsstop called from dummy HDFREAD, library not linked',1)
  RETURN
END SUBROUTINE hdfread

SUBROUTINE hdfreadsplit
  WRITE (6,'(/a/a/)')                                                   &
      'HDF library was not linked. No HDF data read created.',          &
      'Make sure you include "-io hdf" or "-io hdfnet" option when doing makearps.'
  CALL arpsstop('arpsstop called from dummy HDFREADSPLIT, library not linked',1)
  RETURN
END SUBROUTINE hdfreadsplit

SUBROUTINE hdfdump
  WRITE (6,'(/a/a/)')                                                   &
      'HDF library was not linked. No HDF data dump created.',          &
      'Make sure you include "-io hdf" or "-io hdfnet" option when doing makearps.'
  CALL arpsstop('arpsstop called from dummy HDFDUMP, library not linked',1)
  RETURN
END SUBROUTINE hdfdump

SUBROUTINE hdfjoindump
  WRITE (6,'(/a/a/)')                                                   &
      'HDF library was not linked. No HDF data dump created.',          &
      'Make sure you include "-io hdf" or "-io hdfnet" option when doing makearps.'
  CALL arpsstop('arpsstop called from dummy HDFJOINDUMP, library not linked',1)
  RETURN
END SUBROUTINE hdfjoindump

SUBROUTINE hdfwrt4d
  RETURN
END SUBROUTINE hdfwrt4d

SUBROUTINE hdfrd4d
  RETURN
END SUBROUTINE hdfrd4d

SUBROUTINE hdfwrt3d
  RETURN
END SUBROUTINE hdfwrt3d

SUBROUTINE hdfwrt3di
  RETURN
END SUBROUTINE hdfwrt3di

SUBROUTINE hdfwrt2d
  RETURN
END SUBROUTINE hdfwrt2d

SUBROUTINE hdfwrt2di
  RETURN
END SUBROUTINE hdfwrt2di

SUBROUTINE hdfwrt1d
  RETURN
END SUBROUTINE hdfwrt1d

SUBROUTINE hdfwrt1di
  RETURN
END SUBROUTINE hdfwrt1di

SUBROUTINE hdfwrtr
  RETURN
END SUBROUTINE hdfwrtr

SUBROUTINE hdfwrti
  RETURN
END SUBROUTINE hdfwrti

SUBROUTINE hdfwrtc
  RETURN
END SUBROUTINE hdfwrtc

SUBROUTINE hdfrd3d
  RETURN
END SUBROUTINE hdfrd3d

SUBROUTINE hdfrd3di
  RETURN
END SUBROUTINE hdfrd3di

SUBROUTINE hdfrd2d
  RETURN
END SUBROUTINE hdfrd2d

SUBROUTINE hdfrd2di
  RETURN
END SUBROUTINE hdfrd2di

SUBROUTINE hdfrd1d
  RETURN
END SUBROUTINE hdfrd1d

SUBROUTINE hdfrd1di
  RETURN
END SUBROUTINE hdfrd1di

SUBROUTINE hdfrdr
  RETURN
END SUBROUTINE hdfrdr

SUBROUTINE hdfrdi
  RETURN
END SUBROUTINE hdfrdi

SUBROUTINE hdfrdc
  RETURN
END SUBROUTINE hdfrdc

SUBROUTINE hdfopen
  WRITE (6,'(/a/a/)')                                                   &
      'HDF library was not linked. No HDF data read created.',          &
      'Make sure you include "-io hdf" or "-io hdfnet" option when doing makearps.'
  CALL arpsstop('arpsstop called from dummy HDFOPEN library not linked',1)
  RETURN
END SUBROUTINE hdfopen

SUBROUTINE hdfclose
  RETURN
END SUBROUTINE hdfclose

INTEGER FUNCTION d8pimg (imgfn , image, idata,jdata, icompres)
  INTEGER :: icompres, idata,jdata
  CHARACTER (LEN=*) :: image
  CHARACTER (LEN=*) :: imgfn

  d8pimg = 1

  RETURN
END FUNCTION d8pimg

SUBROUTINE get_dims_from_hdf
  RETURN
END SUBROUTINE get_dims_from_hdf

SUBROUTINE get_var_attr_from_hdf
  RETURN
END SUBROUTINE get_var_attr_from_hdf

SUBROUTINE hdf_read_head
  RETURN
END SUBROUTINE hdf_read_head

SUBROUTINE hdf_read_grid
  RETURN
END SUBROUTINE hdf_read_grid

SUBROUTINE hdf_read_base
  RETURN
END SUBROUTINE hdf_read_base

SUBROUTINE hdf_read_3d
  RETURN
END SUBROUTINE hdf_read_3d

SUBROUTINE hdf_read_static
  RETURN
END SUBROUTINE hdf_read_static

SUBROUTINE hdf_read_optional
  RETURN
END SUBROUTINE hdf_read_optional

SUBROUTINE hdf_read_dsd
  RETURN
END SUBROUTINE hdf_read_dsd

SUBROUTINE hdf_read_var
  RETURN
END SUBROUTINE hdf_read_var
