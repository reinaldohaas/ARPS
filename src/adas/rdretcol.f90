!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE RDRETCOL                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE rdretcol(nx,ny,nz,nvar_ret,                                  &
           mx_ret,nz_ret,mx_colret,nretfil,fretname,                    &
           srcret,isrcret,stnret,latret,lonret,elvret,                  &
           latretc,lonretc,iret,nlevret,hgtretc,obsret,qrret,           &
           ncolret,istatus,tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Reads retrieval data stored as columns, i.e. psuedo-soundings.
!  This allows the retrieval to occur on a different grid
!  than the analysis.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  August, 1995
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nvar_ret number of variables in the obsret array
!    mx_ret     maximum number of retrieval radars
!    nz_ret     maximum number of levels in a retreival columns
!    mx_colret  maximum number of retrieval columns
!
!    nretfil    number of retrieval files
!    fretname   file name for retrieval datasets
!    srcret     name of retrieval sources
!
!  OUTPUT:
!
!    isrcret  index of retrieval source
!    stnret   retrieval radar site name    character*4
!    latret   latitude of retrieval radar  (degrees N)
!    lonret   longitude of retrieval radar (degrees E)
!    elvret   elevation of feed horn of retrieval radar (m MSL)
!    latretc  latitude of retrieval column   (degrees N)
!    lonretc  longitude of retrieval column  (degrees E)
!    iret     retrieval radar number of each column
!    nlevret  number of levels of retrieval data in each column
!    hgtretc  height (m MSL) of retrieval observations
!    obsret   retrieval observations
!    qrret    retrieval qr
!    ncolret  number of retrieval columns read-in
!    istatus  status indicator
!
!    tem1     Temporary work array.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nx,ny,nz,nvar_ret,mx_ret,nz_ret,mx_colret

  INTEGER :: nretfil
  CHARACTER (LEN=*) :: fretname(mx_ret)
!
!-----------------------------------------------------------------------
!
!  Radar site variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=8) :: srcret(mx_ret)
  INTEGER           :: isrcret(mx_ret)
  REAL              :: latret(mx_ret),lonret(mx_ret)
  REAL              :: elvret(mx_ret)
  CHARACTER (LEN=5) :: stnret(mx_ret)
!
!-----------------------------------------------------------------------
!
!  Retrieval observation variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iret(mx_colret)
  INTEGER :: nlevret(mx_colret)
  REAL :: latretc(mx_colret),lonretc(mx_colret)
  REAL :: hgtretc(nz_ret,mx_colret)
  REAL :: obsret(nvar_ret,nz_ret,mx_colret)
  REAL :: qrret(nz_ret,mx_colret)
  INTEGER :: ncolret
  INTEGER :: istatus
!
!-----------------------------------------------------------------------
!
!  Temporary work array
!
!-----------------------------------------------------------------------
!
  REAL :: tem1(nx*ny*nz)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=4)   :: stn
  CHARACTER (LEN=80)  :: runname
  CHARACTER (LEN=256) :: fname
  INTEGER :: ireftim,itime,vcpnum,idummy
  INTEGER :: hdmpfmt,strhopt,mapprin
  INTEGER :: nchanl
  INTEGER :: ierr
  INTEGER :: iyr, imon, idy, ihr, imin, isec
  INTEGER :: i,j,icstrt,icol,klev,kk,kret,nfile,maxk

  REAL :: dxin,dyin,dzin,dzminin,ctrlatin
  REAL :: ctrlonin,tlat1in,tlat2in,tlonin,scalin,rdummy
  REAL :: xrd,yrd,elev
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  maxk=0
  icol=0
  istatus=0
  icstrt=1
  nfile=nretfil
  IF(nretfil > mx_ret) THEN
    WRITE(6,'(a,i3,a,i3/a,i3,a)')                                       &
        ' WARNING nretfil ',nretfil,' exceeds mx_ret dimension ',       &
        mx_ret,' only ',mx_ret,' files will be read.'
    nfile=mx_ret
  END IF
!
!-----------------------------------------------------------------------
!
!  Loop through all retrieval files
!
!-----------------------------------------------------------------------
!
  DO kret=1,nfile
    fname=fretname(kret)
    CALL asnctl ('NEWLOCAL', 1, ierr)
    CALL asnfile(fname, '-F f77 -N ieee', ierr)

    CALL getunit( nchanl )
    OPEN(UNIT=nchanl,FILE=trim(fname),ERR=390,                          &
         FORM='unformatted',STATUS='old')

    istatus=1
    isrcret(kret)=1
!
    READ(nchanl) stn
    stnret(kret)=stn
    READ(nchanl) ireftim,itime,vcpnum,idummy,idummy,                    &
             idummy,idummy,idummy,idummy,idummy
!
    READ(nchanl) runname
    READ(nchanl) hdmpfmt,strhopt,mapprin,idummy,idummy,                 &
             idummy,idummy,idummy,idummy,idummy
    READ(nchanl) dxin,dyin,dzin,dzminin,ctrlatin,                       &
             ctrlonin,tlat1in,tlat2in,tlonin,scalin,                    &
             latret(kret),lonret(kret),elvret(kret),                    &
             rdummy,rdummy
!
    CALL abss2ctim(itime, iyr, imon, idy, ihr, imin, isec )
    iyr=MOD(iyr,100)
    WRITE(6,815) 'Reading retrieval data for: ',                        &
                  imon,idy,iyr,ihr,imin
    815   FORMAT(/a,i2.2,'/',i2.2,'/',i2.2,1X,i2.2,':',i2.2,' UTC')
!
!-----------------------------------------------------------------------
!
!  Note here the retrieval data indices:
!
!       1 u wind component
!       2 u wind component
!       3 pressure
!       4 potential temperature
!       5 specific humidity
!
!-----------------------------------------------------------------------
!
    DO icol=icstrt,mx_colret
      READ(nchanl,END=51) i,j,xrd,yrd,                                  &
          latretc(icol),lonretc(icol),elev,klev
      maxk=MAX(maxk,klev)
      klev=MIN(klev,nz_ret)
      nlevret(icol)=klev
      iret(icol)=kret
!      print *, ' i,j,xrd,yrd: ',i,j,xrd,yrd
!      print *, ' lat,lon,klev: ',
!    :               latretc(icol),lonretc(icol),klev
      READ(nchanl,END=52) (tem1(kk),kk=1,klev)
      READ(nchanl,END=52) (hgtretc(kk,icol),kk=1,klev)
      READ(nchanl,END=52) (obsret(1,kk,icol),kk=1,klev)
      READ(nchanl,END=52) (obsret(2,kk,icol),kk=1,klev)
      READ(nchanl,END=52) (obsret(3,kk,icol),kk=1,klev)
      READ(nchanl,END=52) (obsret(4,kk,icol),kk=1,klev)
      READ(nchanl,END=52) (obsret(5,kk,icol),kk=1,klev)
      READ(nchanl,END=52) (qrret(kk,icol),kk=1,klev)
      READ(nchanl,END=52) (tem1(kk),kk=1,klev)
!      print *, ' icol,hgt(5)=',icol,hgtretc(5,icol)
!      print *, ' u(5),v(5) = ',obsret(1,5,icol),obsret(2,5,icol)
!      print *, ' p(5),pt(5)= ',obsret(3,5,icol),obsret(4,5,icol)
!      print *, ' qv(5)     = ',obsret(5,5,icol)
!
!      Temporary code to set the thermo variables to missing
!      to test, individually, the contributions of each variable.
!
!      DO 45 kk=1,klev
!        obsret(3,kk,icol)=-999.
!        obsret(4,kk,icol)=-999.
!        obsret(5,kk,icol)=-999.
!  45     CONTINUE
!
    END DO
    icol=icol-1
    WRITE(6,'(a,i6,a)')                                                 &
              ' WARNING ran out space, increase mx_colret ',            &
                       icol,' columns'
    GO TO 55
    51   CONTINUE
    icol=icol-1
    WRITE(6,'(a,i6,a)') ' End of file reached after reading',           &
                       icol,' columns'
    GO TO 55
    52   CONTINUE
    WRITE(6,'(a,i6,a)') ' End of file reached while reading',           &
                       icol,' column'
    55   CONTINUE
    CLOSE(nchanl)
    CALL retunit( nchanl )
    icstrt=icol+1
    CYCLE
    390   CONTINUE
    WRITE(6,'(a,a)') ' Error opening file: ',fname
  END DO
  ncolret=icol
  IF (maxk > 0) WRITE(6,'(a,i5)') ' Maximum number of vert levels read:',maxk
  IF(maxk > nz_ret) THEN
    WRITE(6,'(a,i5)')                                                   &
            '   EXCEEDS nz_ret, increase nz_ret:',nz_ret
    STOP
  END IF
  RETURN
END SUBROUTINE rdretcol
