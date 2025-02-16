SUBROUTINE wtrftiltcdf(nazim,ngate,fname,outdir,varname,radname,        &
                       radlat,radlon,radelv,ivcp,elv,                   &
                       rmisval,rngfval,itimcdf,frtime,initime,          &
                       vnyquist,rfrgate,                                &
                       azim,beamw,gtspc,refl)
  IMPLICIT NONE
  INTEGER :: nazim
  INTEGER :: ngate
  CHARACTER (LEN=256) :: fname
  CHARACTER (LEN=80)  :: outdir
  CHARACTER (LEN=40)  :: varname
  CHARACTER (LEN=4)   :: radname
  REAL :: radlat
  REAL :: radlon
  REAL :: radelv
  INTEGER :: ivcp
  REAL :: elv
  REAL :: rmisval
  REAL :: rngfval
  INTEGER :: itimcdf
  REAL :: frtime
  INTEGER :: initime
  REAL :: vnyquist
  REAL :: rfrgate
!
  REAL :: azim(nazim)
  REAL :: beamw(nazim)
  REAL :: gtspc(nazim)
  REAL :: refl(ngate,nazim)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  WRITE(6,'(//a)') ' You have called a dummy subroutine wrtrftiltcdf'
  WRITE(6,'(a)') ' Recompile the program using -io net'
  STOP
END SUBROUTINE wtrftiltcdf

SUBROUTINE wtvrtiltcdf(nazim,ngate,fname,outdir,varname,radname,        &
                       radlat,radlon,radelv,ivcp,elv,                   &
                       rmisval,rngfval,itimcdf,frtime,initime,          &
                       vnyquist,rfrgate,                                &
                       azim,beamw,gtspc,vnyq,radv)
  IMPLICIT NONE
  INTEGER :: nazim
  INTEGER :: ngate
  CHARACTER (LEN=256) :: fname
  CHARACTER (LEN=80)  :: outdir
  CHARACTER (LEN=40)  :: varname
  CHARACTER (LEN=4)   :: radname
  REAL :: radlat
  REAL :: radlon
  REAL :: radelv
  INTEGER :: ivcp
  REAL :: elv
  REAL :: rmisval
  REAL :: rngfval
  INTEGER :: itimcdf
  REAL :: frtime
  INTEGER :: initime
  REAL :: vnyquist
  REAL :: rfrgate

  REAL :: azim(nazim)
  REAL :: beamw(nazim)
  REAL :: gtspc(nazim)
  REAL :: vnyq(nazim)
  REAL :: radv(ngate,nazim)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  WRITE(6,'(//a)') ' You have called a dummy subroutine wrtrftiltcdf'
  WRITE(6,'(a)') ' Recompile the program using -io net'
  STOP
END SUBROUTINE wtvrtiltcdf

SUBROUTINE wtvvtiltcdf(nazim,ngate,fname,outdir,varname,radname,        &
                       radlat,radlon,radelv,ivcp,elv,                   &
                       rmisval,rngfval,itimcdf,frtime,initime,          &
                       vnyquist,rfrgate,                                &
                       azim,beamw,gtspc,vvort)
  IMPLICIT NONE
  INTEGER :: nazim
  INTEGER :: ngate
  CHARACTER (LEN=256) :: fname
  CHARACTER (LEN=80)  :: outdir
  CHARACTER (LEN=40)  :: varname
  CHARACTER (LEN=4)   :: radname
  REAL :: radlat
  REAL :: radlon
  REAL :: radelv
  INTEGER :: ivcp
  REAL :: elv
  REAL :: rmisval
  REAL :: rngfval
  INTEGER :: itimcdf
  REAL :: frtime
  INTEGER :: initime
  REAL :: vnyquist
  REAL :: rfrgate
!
  REAL :: azim(nazim)
  REAL :: beamw(nazim)
  REAL :: gtspc(nazim)
  REAL :: vvort(ngate,nazim)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  WRITE(6,'(//a)') ' You have called a dummy subroutine wrtrftiltcdf'
  WRITE(6,'(a)') ' Recompile the program using -io net'
  STOP
END SUBROUTINE wtvvtiltcdf
