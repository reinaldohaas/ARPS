
SUBROUTINE open_phdf5_for_write(filename,sysdepinfo,Hndl,initialized,ireturn)
  IMPLICIT NONE

  CHARACTER(*), INTENT(IN)  :: filename
  CHARACTER(*), INTENT(IN)  :: sysdepinfo
  LOGICAL,      INTENT(IN)  :: initialized
  INTEGER,      INTENT(OUT) :: Hndl
  INTEGER,      INTENT(OUT) :: ireturn

  WRITE (6,'(1x,/,a,/,10x,2a,/)')                                       &
      ' WARNING: PHDF5 library was not linked.','Please add ',          &
      '"-io phdf5" option for makearps to link PHDF5 libraries.'
  CALL arpsstop('arpsstop called from open_phdf5_for_write.',1)

  Hndl    = -1
  ireturn = -1
  RETURN
END SUBROUTINE open_phdf5_for_write

SUBROUTINE  close_phdf5_for_write(Hndl,ireturn)

  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: Hndl
  INTEGER, INTENT(OUT) :: ireturn

  ireturn = -1
  RETURN
END SUBROUTINE close_phdf5_for_write

SUBROUTINE shutdown_phdf5_io(ireturn)
  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: ireturn

  ireturn = -1
  RETURN
END SUBROUTINE

SUBROUTINE put_phdf5_dom_ti_char(Hndl,attname,attstr,ireturn)

  IMPLICIT NONE
  INTEGER,      INTENT(IN)  :: Hndl
  CHARACTER(*), INTENT(IN)  :: attname
  CHARACTER(*), INTENT(IN)  :: attstr
  INTEGER,      INTENT(OUT) :: ireturn

  ireturn = -1
  RETURN
END SUBROUTINE put_phdf5_dom_ti_char

SUBROUTINE put_phdf5_dom_ti_integer(Hndl,attname,attval,ireturn)

  IMPLICIT NONE
  INTEGER,      INTENT(IN)  :: Hndl
  CHARACTER(*), INTENT(IN)  :: attname
  INTEGER,      INTENT(IN)  :: attval
  INTEGER,      INTENT(OUT) :: ireturn

  ireturn = -1
  RETURN
END SUBROUTINE put_phdf5_dom_ti_integer

SUBROUTINE put_phdf5_dom_ti_real(Hndl,attname,attval,attsiz,ireturn)

  IMPLICIT NONE
  INTEGER,      INTENT(IN)  :: Hndl
  CHARACTER(*), INTENT(IN)  :: attname
  REAL,         INTENT(IN)  :: attval
  INTEGER,      INTENT(IN)  :: attsiz
  INTEGER,      INTENT(OUT) :: ireturn

  ireturn = -1
  RETURN
END SUBROUTINE put_phdf5_dom_ti_real

SUBROUTINE put_phdf5_dom_td_char(Hndl,varname,CurrDate,DateStr,ireturn)

  IMPLICIT NONE

  INTEGER,      INTENT(IN)  :: Hndl
  CHARACTER(*), INTENT(IN)  :: varname
  CHARACTER(*), INTENT(IN)  :: CurrDate
  CHARACTER(*), INTENT(IN)  :: DateStr
  INTEGER,      INTENT(OUT) :: ireturn

  ireturn = -1
  RETURN
END SUBROUTINE put_phdf5_dom_td_char

SUBROUTINE write_phdf5_field(Hndl,DateStr,VarName,Desc,Units,Stagger,   &
                             Field,FieldType, DomainDesc,               &
                             MemoryOrder,DimNames,                      &
                             DomainStart,DomainEnd,                     &
                             MemoryStart,MemoryEnd,                     &
                             PatchStart,PatchEnd, ireturn)
  IMPLICIT NONE

  INTEGER,       INTENT(IN)  :: Hndl
  CHARACTER*(*), INTENT(IN)  :: DateStr
  CHARACTER*(*), INTENT(IN)  :: VarName
  CHARACTER*(*), INTENT(IN)  :: Desc
  CHARACTER*(*), INTENT(IN)  :: Units
  CHARACTER*(*), INTENT(IN)  :: Stagger

  INTEGER,       INTENT(IN)  :: Field
  INTEGER,       INTENT(IN)  :: FieldType
  INTEGER,       INTENT(IN)  :: DomainDesc
  CHARACTER*(*), INTENT(IN)  :: MemoryOrder
  CHARACTER*(*), INTENT(IN)  :: DimNames
  INTEGER,       INTENT(IN)  :: DomainStart(3), DomainEnd(3)
  INTEGER,       INTENT(IN)  :: MemoryStart(3), MemoryEnd(3)
  INTEGER,       INTENT(IN)  :: PatchStart(3),  PatchEnd(3)
  INTEGER,       INTENT(OUT) :: ireturn

  ireturn = -1

  RETURN
END SUBROUTINE write_phdf5_field

SUBROUTINE write_phdf5_bdy(Hndl,DateStr,VarName,Desc,Units,Stagger,   &
                             Field,FieldType, DomainDesc,               &
                             MemoryOrder,DimNames,                      &
                             DomainStart,DomainEnd,                     &
                             MemoryStart,MemoryEnd,                     &
                             PatchStart,PatchEnd, IOFlag,ireturn)
  IMPLICIT NONE

  INTEGER,       INTENT(IN)  :: Hndl
  CHARACTER*(*), INTENT(IN)  :: DateStr
  CHARACTER*(*), INTENT(IN)  :: VarName
  CHARACTER*(*), INTENT(IN)  :: Desc
  CHARACTER*(*), INTENT(IN)  :: Units
  CHARACTER*(*), INTENT(IN)  :: Stagger

  INTEGER,       INTENT(IN)  :: Field
  INTEGER,       INTENT(IN)  :: FieldType
  INTEGER,       INTENT(IN)  :: DomainDesc
  CHARACTER*(*), INTENT(IN)  :: MemoryOrder
  CHARACTER*(*), INTENT(IN)  :: DimNames
  INTEGER,       INTENT(IN)  :: DomainStart(3), DomainEnd(3)
  INTEGER,       INTENT(IN)  :: MemoryStart(3), MemoryEnd(3)
  INTEGER,       INTENT(IN)  :: PatchStart(3),  PatchEnd(3)
  LOGICAL,       INTENT(IN)  :: IOFlag
  INTEGER,       INTENT(OUT) :: ireturn

  ireturn = -1

  RETURN
END SUBROUTINE write_phdf5_bdy
