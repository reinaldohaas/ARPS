SUBROUTINE open_phdf5_for_read (FileName, SysDepInfo, initialized, Hndl, iStatus)
  IMPLICIT NONE

  CHARACTER(*), INTENT(IN)  :: filename
  CHARACTER(*), INTENT(IN)  :: SysDepInfo
  LOGICAL,      INTENT(IN)  :: initialized 
  INTEGER,      INTENT(OUT) :: Hndl
  INTEGER,      INTENT(OUT) :: iStatus
 
  INCLUDE  'mpif.h'
 
  INTEGER :: Comm_compute , Comm_io

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code ... ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Comm_compute = MPI_COMM_WORLD
  Comm_io      = MPI_COMM_WORLD

  IF(.NOT. initialized) CALL ext_phdf5_ioinit(sysdepinfo,iStatus)

  CALL ext_phdf5_open_for_read ( FileName , Comm_compute, Comm_io, SysDepInfo, &
                               Hndl , iStatus )
  
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
 
  CALL ext_phdf5_ioclose( Hndl, ireturn )
 
  RETURN
END SUBROUTINE close_phdf5_for_read

SUBROUTINE shutdown_phdf5_io(ireturn)
  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: ireturn
 
  CALL ext_phdf5_ioexit(ireturn)
  RETURN
END SUBROUTINE

SUBROUTINE get_phdf5_next_time(Hndl,DateStr,istatus)
  IMPLICIT NONE
  INTEGER,      INTENT(IN)  :: Hndl
  CHARACTER(*), INTENT(OUT) :: DateStr
  INTEGER,      INTENT(OUT) :: istatus

  CALL ext_phdf5_get_next_time( Hndl, DateStr, iStatus )

  RETURN
END SUBROUTINE

SUBROUTINE get_phdf5_dom_ti_integer(Hndl, Element, Data, Status)

  IMPLICIT NONE
  INTEGER,      INTENT(IN) :: Hndl
  CHARACTER(*), INTENT(IN) :: Element
  INTEGER,      INTENT(OUT):: Data
  INTEGER,      INTENT(OUT):: Status

  INTEGER :: Outcount

  CALL ext_phdf5_get_dom_ti_integer ( Hndl, Element,   Data, &
                              1, Outcount, Status )

  RETURN
END SUBROUTINE get_phdf5_dom_ti_integer

SUBROUTINE get_phdf5_dom_ti_real(Hndl, Element, Data, Status)

  IMPLICIT NONE
  INTEGER,      INTENT(IN) :: Hndl
  CHARACTER(*), INTENT(IN) :: Element
  REAL,         INTENT(OUT):: Data
  INTEGER,      INTENT(OUT):: Status

  INTEGER :: Outcount

  CALL ext_phdf5_get_dom_ti_real ( Hndl, Element,   Data, &
                              1, Outcount, Status )

  RETURN
END SUBROUTINE get_phdf5_dom_ti_real

SUBROUTINE get_phdf5_dom_ti_char(Hndl, Element, Data, Status)

  IMPLICIT NONE
  INTEGER,      INTENT(IN) :: Hndl
  CHARACTER(*), INTENT(IN) :: Element
  CHARACTER(*), INTENT(OUT):: Data
  INTEGER,      INTENT(OUT):: Status

  CALL ext_phdf5_get_dom_ti_char ( Hndl, Element,   Data, &
                                   Status )

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

  INTEGER :: Comm, IOComm

  INCLUDE  'mpif.h'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Comm   = MPI_COMM_WORLD
  IOComm = MPI_COMM_WORLD

!  CALL ext_phdf5_get_var_ti_char ( Hndl,'MemoryOrder' ,  Varname, memoryorder, &
!                               iStatus )
!
!  CALL ext_phdf5_get_var_ti_char ( Hndl,'stagger' ,  Varname, stagger, &
!                               iStatus )

  CALL ext_phdf5_read_field   ( Hndl , DateStr , VarName ,              &
                                Field , FieldType ,                     &
                                Comm , IOComm , DomainDesc ,            &
                                MemoryOrder , Stagger , DimNames ,      &
                                DomainStart , DomainEnd ,               &
                                MemoryStart , MemoryEnd ,               &
                                PatchStart , PatchEnd ,                 &
                                iStatus )
  RETURN
END SUBROUTINE get_phdf5_field
