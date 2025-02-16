SUBROUTINE open_phdf5_for_read (FileName, SysDepInfo, initialized, Hndl, iStatus)
  IMPLICIT NONE

  CHARACTER(*), INTENT(IN)  :: filename
  CHARACTER(*), INTENT(IN)  :: SysDepInfo
  LOGICAL,      INTENT(IN)  :: initialized 
  INTEGER,      INTENT(OUT) :: Hndl
  INTEGER,      INTENT(OUT) :: iStatus
 

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code ... ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Hndl = -1
  iStatus = 0
  WRITE(6,*) 'ERROR: PHDF5 library was not linked. Use the following option to link:'
  WRITE(6,*) '       makearps -io phdf5 wrf2arps_mpi'
  CALL arpsstop('phdf5 option was not specified.',1) 
  
  RETURN
END SUBROUTINE open_phdf5_for_read

SUBROUTINE  close_phdf5_for_read(Hndl,ireturn)
 
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: Hndl
  INTEGER, INTENT(OUT) :: ireturn
 
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code .....
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
  ireturn = 0
  RETURN
END SUBROUTINE close_phdf5_for_read

SUBROUTINE shutdown_phdf5_io(ireturn)
  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: ireturn
 
  ireturn = 0
  RETURN
END SUBROUTINE

SUBROUTINE get_phdf5_next_time(Hndl,DateStr,istatus)
  IMPLICIT NONE
  INTEGER,      INTENT(IN)  :: Hndl
  CHARACTER(*), INTENT(OUT) :: DateStr
  INTEGER,      INTENT(OUT) :: istatus

  DateStr = ' '
  istatus = 0
  RETURN
END SUBROUTINE

SUBROUTINE get_phdf5_dom_ti_integer(Hndl, Element, Data, Status)

  IMPLICIT NONE
  INTEGER,      INTENT(IN) :: Hndl
  CHARACTER(*), INTENT(IN) :: Element
  INTEGER,      INTENT(OUT):: Data
  INTEGER,      INTENT(OUT):: Status

  Data = 0
  Status = 0
  RETURN
END SUBROUTINE get_phdf5_dom_ti_integer

SUBROUTINE get_phdf5_dom_ti_real(Hndl, Element, Data, Status)

  IMPLICIT NONE
  INTEGER,      INTENT(IN) :: Hndl
  CHARACTER(*), INTENT(IN) :: Element
  REAL,         INTENT(OUT):: Data
  INTEGER,      INTENT(OUT):: Status

  Data = 0
  Status = 0
  RETURN
END SUBROUTINE get_phdf5_dom_ti_real

SUBROUTINE get_phdf5_dom_ti_char(Hndl, Element, Data, Status)

  IMPLICIT NONE
  INTEGER,      INTENT(IN) :: Hndl
  CHARACTER(*), INTENT(IN) :: Element
  CHARACTER(*), INTENT(OUT):: Data
  INTEGER,      INTENT(OUT):: Status

  Data = ' '
  Status = 0
  RETURN
END SUBROUTINE get_phdf5_dom_ti_char

SUBROUTINE get_phdf5_field(Hndl, DateStr, VarName, Field, FieldType,    &
                           DomainDesc,MemoryOrder,Stagger,DimNames,     &
                           DomainStart,DomainEnd,MemoryStart,MemoryEnd, &
                           PatchStart,PatchEnd,                         &
                           iStatus)

  IMPLICIT NONE
  INTEGER,      INTENT(IN)  :: Hndl
  CHARACTER(*), INTENT(IN)  :: DateStr
  CHARACTER(*), INTENT(IN)  :: VarName
  INTEGER,      INTENT(IN)  :: FieldType
  INTEGER,      INTENT(IN)  :: DomainDesc
  CHARACTER(*), INTENT(IN)  :: MemoryOrder
  CHARACTER(*), INTENT(IN)  :: Stagger

  CHARACTER(*), INTENT(IN)  :: DimNames(3)
  INTEGER,      INTENT(IN)  :: DomainStart(3)
  INTEGER,      INTENT(IN)  :: DomainEnd(3)
  INTEGER,      INTENT(IN)  :: MemoryStart(3)
  INTEGER,      INTENT(IN)  :: MemoryEnd(3)
  INTEGER,      INTENT(IN)  :: PatchStart(3)
  INTEGER,      INTENT(IN)  :: PatchEnd(3)

  INTEGER,      INTENT(OUT) :: Field(*)
  INTEGER,      INTENT(OUT) :: iStatus


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Field(1) = 0
  istatus = -1

  RETURN
END SUBROUTINE get_phdf5_field
