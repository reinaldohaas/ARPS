SUBROUTINE xdevic_new(IWTYPE,filename,len_filename,iounitin)
  IMPLICIT NONE

  INTEGER,            INTENT(IN)    :: IWTYPE
  CHARACTER(LEN=256), INTENT(IN)    :: filename             ! output graphical file name (base)
  INTEGER,            INTENT(INOUT) :: len_filename         ! valid character length in filename
                                                            ! = 0 to use the default "gmeta" for backward compatibility.
  INTEGER,            INTENT(IN)    :: iounitin             ! A file unit for output, used only by PS interface.
                                                            ! It is ignored here.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  len_filename = -1      ! to let us know it is void.

  RETURN
END SUBROUTINE xdevic_new

SUBROUTINE xpspac
  RETURN
END SUBROUTINE xpspac

SUBROUTINE xmap
  RETURN
END SUBROUTINE xmap

SUBROUTINE xaxsca
  RETURN
END SUBROUTINE xaxsca

SUBROUTINE xnctrs
  RETURN
END SUBROUTINE xnctrs

SUBROUTINE xconta
  RETURN
END SUBROUTINE xconta

SUBROUTINE xclimt
  RETURN
END SUBROUTINE xclimt

SUBROUTINE xframe
  RETURN
END SUBROUTINE xframe

SUBROUTINE xgrend
  RETURN
END SUBROUTINE xgrend

SUBROUTINE xctrbadv
  RETURN
END SUBROUTINE xctrbadv

SUBROUTINE xvtrbadv
  RETURN
END SUBROUTINE xvtrbadv

SUBROUTINE xbadval
  RETURN
END SUBROUTINE xbadval

SUBROUTINE xcharc
  RETURN
END SUBROUTINE xcharc

SUBROUTINE xaxfmt
  RETURN
END SUBROUTINE xaxfmt

SUBROUTINE xclfmt
  RETURN
END SUBROUTINE xclfmt

SUBROUTINE xcolfil
  RETURN
END SUBROUTINE xcolfil

SUBROUTINE xstpjgrd
  RETURN
END SUBROUTINE xstpjgrd

SUBROUTINE xctrlim
  RETURN
END SUBROUTINE xctrlim

SUBROUTINE xsetclrs_new
  RETURN
END SUBROUTINE xsetclrs_new

SUBROUTINE xctrclr
  RETURN
END SUBROUTINE xctrclr

SUBROUTINE xcolor
  RETURN
END SUBROUTINE xcolor

SUBROUTINE xdrawmap
  RETURN
END SUBROUTINE xdrawmap

SUBROUTINE xcpalet
  RETURN
END SUBROUTINE xcpalet

SUBROUTINE xchsiz
  RETURN
END SUBROUTINE xchsiz

SUBROUTINE xaxsor
  RETURN
END SUBROUTINE xaxsor

SUBROUTINE xcontcopt
  RETURN
END SUBROUTINE xcontcopt
