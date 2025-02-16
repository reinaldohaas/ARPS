
SUBROUTINE open_phdf5_for_write(filename,sysdepinfo,Hndl,initialized,ireturn)
  IMPLICIT NONE

  CHARACTER(*), INTENT(IN)  :: filename
  CHARACTER(*), INTENT(IN)  :: sysdepinfo
  LOGICAL,      INTENT(IN)  :: initialized
  INTEGER,      INTENT(OUT) :: Hndl
  INTEGER,      INTENT(OUT) :: ireturn

  INCLUDE  'mpif.h'
  
  INTEGER :: Comm_compute , Comm_io

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code .....
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Comm_compute = MPI_COMM_WORLD
  Comm_io      = MPI_COMM_WORLD

  IF(.NOT. initialized) CALL ext_phdf5_ioinit(sysdepinfo,ireturn)
  CALL ext_phdf5_open_for_write_begin(FileName, Comm_compute, Comm_io,  &
                     SysDepInfo, Hndl, ireturn)

  CALL ext_phdf5_open_for_write_commit ( Hndl , ireturn )

 RETURN
END SUBROUTINE open_phdf5_for_write

SUBROUTINE  close_phdf5_for_write(Hndl,ireturn)

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
END SUBROUTINE close_phdf5_for_write

SUBROUTINE shutdown_phdf5_io(ireturn)
  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: ireturn

  CALL ext_phdf5_ioexit(ireturn)
  RETURN
END SUBROUTINE

SUBROUTINE put_phdf5_dom_ti_char(Hndl,attname,attstr,ireturn)

  IMPLICIT NONE
  INTEGER,      INTENT(IN)  :: Hndl
  CHARACTER(*), INTENT(IN)  :: attname
  CHARACTER(*), INTENT(IN)  :: attstr
  INTEGER,      INTENT(OUT) :: ireturn

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code .....
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL ext_phdf5_put_dom_ti_char(Hndl,TRIM(attname),TRIM(attstr),ireturn)

  RETURN
END SUBROUTINE put_phdf5_dom_ti_char

SUBROUTINE put_phdf5_dom_ti_integer(Hndl,attname,attval,ireturn)

  IMPLICIT NONE
  INTEGER,      INTENT(IN)  :: Hndl
  CHARACTER(*), INTENT(IN)  :: attname
  INTEGER,      INTENT(IN)  :: attval
  INTEGER,      INTENT(OUT) :: ireturn

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code .....
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL ext_phdf5_put_dom_ti_integer(Hndl,TRIM(attname),attval,1,ireturn)

  RETURN
END SUBROUTINE put_phdf5_dom_ti_integer

SUBROUTINE put_phdf5_dom_ti_real(Hndl,attname,attval,attsiz,ireturn)

  IMPLICIT NONE
  INTEGER,      INTENT(IN)  :: Hndl
  CHARACTER(*), INTENT(IN)  :: attname
  REAL,         INTENT(IN)  :: attval
  INTEGER,      INTENT(IN)  :: attsiz
  INTEGER,      INTENT(OUT) :: ireturn

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code .....
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL ext_phdf5_put_dom_ti_real(Hndl,TRIM(attname),attval,attsiz,ireturn)

  RETURN
END SUBROUTINE put_phdf5_dom_ti_real

SUBROUTINE put_phdf5_dom_td_char(Hndl,varname,CurrDate,DateStr,ireturn)

  IMPLICIT NONE

  INTEGER,      INTENT(IN)  :: Hndl
  CHARACTER(*), INTENT(IN)  :: varname
  CHARACTER(*), INTENT(IN)  :: CurrDate
  CHARACTER(*), INTENT(IN)  :: DateStr
  INTEGER,      INTENT(OUT) :: ireturn

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code .....
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL ext_phdf5_put_dom_td_char(Hndl,varname,CurrDate(1:19),           &
                                 DateStr(1:19),ireturn)

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

  INCLUDE  'mpif.h'

  INTEGER  :: Comm, IOComm

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code ... ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Comm   = MPI_COMM_WORLD
  IOComm = MPI_COMM_WORLD

  CALL ext_phdf5_write_field( Hndl , DateStr , VarName ,                &
                              Field , FieldType , Comm , IOComm ,       &
                              DomainDesc , MemoryOrder , Stagger , DimNames , &
                              DomainStart, DomainEnd,                   &
                              MemoryStart, MemoryEnd,                   &
                              PatchStart,  PatchEnd,                    &
                              ireturn )
  IF (ireturn /= 0) WRITE(0,*) 'ERROR in ext_phdf5_write_field for ',VarName

  CALL ext_phdf5_put_var_ti_char(Hndl,'description',Varname,Desc,ireturn)
  IF (ireturn /= 0) WRITE(0,*) 'ERROR in ext_phdf5_put_var_ti_char,     &
                              & description for ',VarName

  CALL ext_phdf5_put_var_ti_char(Hndl,'units',      Varname,Units, ireturn )
  IF (ireturn /= 0) WRITE(0,*) 'ERROR in ext_phdf5_put_var_ti_char,     &
                              & units for ',VarName

  CALL ext_phdf5_put_var_ti_char(Hndl,'stagger',    Varname,Stagger,ireturn)
  IF (ireturn /= 0) WRITE(0,*) 'ERROR in ext_phdf5_put_var_ti_char,     &
                              & stagger for ',VarName

  RETURN
END SUBROUTINE write_phdf5_field

SUBROUTINE write_phdf5_bdy(Hndl,DateStr,VarName,Desc,Units,Stagger,   &
                             Field,FieldType, DomainDesc,               &
                             MemoryOrder,DimNames,                      &
                             DomainStart,DomainEnd,                     &
                             MemoryStart,MemoryEnd,                     &
                             PatchStart,PatchEnd,                       &
                             IOFLAG, ireturn)
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
  LOGICAL,       INTENT(IN)  :: IOFLAG
  INTEGER,       INTENT(OUT) :: ireturn

  INCLUDE  'mpif.h'

  INTEGER  :: Comm, IOComm

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code ... ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Comm   = MPI_COMM_WORLD
  IOComm = MPI_COMM_WORLD

  CALL ext_phdf5_write_bdy( Hndl , DateStr , VarName ,                &
                            Field , FieldType , Comm , IOComm ,       &
                            DomainDesc , MemoryOrder , Stagger , DimNames , &
                            DomainStart, DomainEnd,                   &
                            MemoryStart, MemoryEnd,                   &
                            PatchStart,  PatchEnd,                    &
                            IOFLAG,      ireturn )
  IF (ireturn /= 0) WRITE(0,*) 'ERROR in ext_phdf5_write_field for ',VarName

  CALL ext_phdf5_put_var_ti_char(Hndl,'description',Varname,Desc,ireturn)
  IF (ireturn /= 0) WRITE(0,*) 'ERROR in ext_phdf5_put_var_ti_char,     &
                              & description for ',VarName

  CALL ext_phdf5_put_var_ti_char(Hndl,'units',      Varname,Units, ireturn )
  IF (ireturn /= 0) WRITE(0,*) 'ERROR in ext_phdf5_put_var_ti_char,     &
                              & units for ',VarName

  CALL ext_phdf5_put_var_ti_char(Hndl,'stagger',    Varname,Stagger,ireturn)
  IF (ireturn /= 0) WRITE(0,*) 'ERROR in ext_phdf5_put_var_ti_char,     &
                              & stagger for ',VarName

  RETURN
END SUBROUTINE write_phdf5_bdy
