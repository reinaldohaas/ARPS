PROGRAM pltradcol
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 PROGRAM PLTRADCOL                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!  PURPOSE:
!
!  Plots data written by wtradcol.  Utilizes zxplot, can be
!  linked to NCAR Graphics or use postscript format.
!
!  makearps -zxncar pltradcol  OR
!  makearps -zxpost pltradcol (default)
!
!  In addition to this program, the results of the radar remapping
!  can be viewed by using the -reffile -velfile -ref2d and -vel2d
!  options in the remappers and the arbvar features of arpsplt
!  (see the arpsplt.input file for details).
!
!  AUTHOR:
!
!  Keith Brewster, CAPS, April, 1996
!
!  MODIFICATIONS:
!
!  01/14/2003 Keith Brewster
!  Updated to use zxplot calls instead of NCARGraphics routines.
!  Added color plotting options and increased user control of plotting
!  variables and levels without modifying source code.  Update included
!  additional NAMELISTS in the pltradcol.input file.
!
!  12/01/2010 Youngsun Jung
!  Added capability to plot data in EnKF grid tilt-column format.
!  A new variable, gridtilt_opt, is added.
!
!  03/29/2012 Keith Brewster
!  Added check of max range to only plot grid within max range.  Helpful
!  for seeing data remapped to very large grid domains.  Also added radar
!  name and filename as command line variables (overrides input file variables)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
!
! mxpltlvl : Maximum number of levels plotted
!
  INTEGER, PARAMETER :: mxpltlvl = 100
  
  INTEGER :: pltref,pltvel,pltnyq,plttim
  INTEGER :: pltzdr,pltkdp,pltrhv
  INTEGER :: pltlvl(mxpltlvl)
!
  CHARACTER (LEN=4) :: radid
  CHARACTER (LEN=5) :: lvlid
  INTEGER :: ireftim,itime
  INTEGER :: vcpnum
  INTEGER :: nradvr
  INTEGER :: irngmin,irngmax,isource
  REAL    :: refelvmin,refelvmax
  REAL    :: elvrad

!-----------------------------------------------------------------------
!
!  Dimension in x, y, and z direction
!
!-----------------------------------------------------------------------

  INTEGER :: nx, ny, nz

!-----------------------------------------------------------------------
!
!  Read-in variables
!
!-----------------------------------------------------------------------
!
  INTEGER, ALLOCATABLE :: readk(:)
  REAL,    ALLOCATABLE :: readhgt(:)
  REAL,    ALLOCATABLE :: readref(:)
  REAL,    ALLOCATABLE :: readvel(:)
  REAL,    ALLOCATABLE :: readnyq(:)
  REAL,    ALLOCATABLE :: readtim(:)
  REAL,    ALLOCATABLE :: readrng(:)
  REAL,    ALLOCATABLE :: readpol(:)
  REAL,    ALLOCATABLE :: readzdr(:)
  REAL,    ALLOCATABLE :: readkdp(:)
  REAL,    ALLOCATABLE :: readrhv(:)
!
!-----------------------------------------------------------------------
!
!  Plotting variables
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: gridref(:,:,:)
  REAL, ALLOCATABLE :: gridvel(:,:,:)
  REAL, ALLOCATABLE :: gridnyq(:,:,:)
  REAL, ALLOCATABLE :: gridtim(:,:,:)
  REAL, ALLOCATABLE :: gridzdr(:,:,:)
  REAL, ALLOCATABLE :: gridkdp(:,:,:)
  REAL, ALLOCATABLE :: gridrhv(:,:,:)
!
!-----------------------------------------------------------------------
!
!  hdf variables and temporary arrays
!
!-----------------------------------------------------------------------
!
  INTEGER, PARAMETER :: mxradvr=10
  INTEGER :: iradvr(mxradvr)

  INTEGER (KIND=selected_int_kind(4)), ALLOCATABLE :: itmp(:,:)

  INTEGER, ALLOCATABLE :: coli(:), colj(:), colk(:,:), numlev(:)
  REAL,    ALLOCATABLE :: radcolref(:,:), radcolvel(:,:)  
  REAL,    ALLOCATABLE :: radcolnyq(:,:), radcoltim(:,:)  

!-----------------------------------------------------------------------
!
!  Statistics variables
!
!-----------------------------------------------------------------------
!
  INTEGER, ALLOCATABLE :: kntref(:)
  INTEGER, ALLOCATABLE :: kntvel(:)
  INTEGER, ALLOCATABLE :: kntzdr(:)
  INTEGER, ALLOCATABLE :: kntkdp(:)
  INTEGER, ALLOCATABLE :: kntrhv(:)
!
  REAL :: refint,velint,nyqint,timint
  PARAMETER (refint=0.,                                                 &
             velint=0.,                                                 &
             nyqint=0.,                                                 &
             timint=0.)
!
!-----------------------------------------------------------------------
!
!  NEXRAD station location variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: nnextab
  INTEGER, PARAMETER :: mxnextab=200
  CHARACTER (LEN=4) :: nextstn(mxnextab)
  REAL :: nextlat(mxnextab)
  REAL :: nextlon(mxnextab)
  REAL :: nextelv(mxnextab)
  CHARACTER (LEN=256) :: nextfile
!
!-----------------------------------------------------------------------
!
!  Background map variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: ovrmap,mapgrid,nmapfile,lmapfile
  INTEGER, PARAMETER  :: mxmapfile= 10
  CHARACTER (LEN=256) :: mapfile(mxmapfile)

  REAL :: latgrid,longrid
  COMMON /mappar2/ latgrid,longrid

!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=3)   :: chplt
  CHARACTER (LEN=18)  :: timplt
  CHARACTER (LEN=256) :: fname
  CHARACTER (LEN=256) :: charg
  INTEGER :: gridtilt_opt
  REAL :: rdummy
  REAL :: x,y,xrd,yrd,gridlat,gridlon,elev
  REAL :: yloc,labsize,mrksize
  REAL :: latnot(2)
  REAL :: ctrx,ctry,swx,swy,dxkm,dykm,hlfdxkm,hlfdykm
  REAL :: x1,x2,y1,y2,pl,pr,pb,pt,px,py,pxc,pyc,xs,ys
  REAL :: latrad,lonrad,chsize,xname,yname
  REAL :: xstn,xstnorg,ystn,ystnorg
  REAL :: xlenchk,xoffset,yoffset
  INTEGER :: i,j,k,kk,klev,ipt,itab,lindex,icol,jcol,munit
  INTEGER :: idummy, iyr, imon, idy, ihr, imin, isec
  INTEGER :: narg,npltlvl
  INTEGER :: dmpfmt,lens,sd_id,iradfmt,ierr,istatus
  LOGICAL :: found,fexist

  REAL, PARAMETER :: misval=-999.

  INTEGER :: typelev, numradcol, nummaxlev
  REAL    :: xmin, xmax, ymin, ymax

  CHARACTER(LEN=256) :: outfilename

  INTEGER :: nextframe = .false.
  INTEGER :: numvar, idp
  INTEGER :: iargc

!-----------------------------------------------------------------------
! 
!  Misc local arrays
! 
!-----------------------------------------------------------------------
!
  INTEGER,ALLOCATABLE :: tilt_time(:)
  REAL,   ALLOCATABLE :: tilt_angle(:)
!
!-----------------------------------------------------------------------
!
!  NAMELIST declarations
!
!-----------------------------------------------------------------------
!
  NAMELIST /file_name/ fname, gridtilt_opt, numvar
  NAMELIST /var_plot/ pltref,pltvel,pltnyq,plttim,pltzdr,pltkdp,pltrhv
  NAMELIST /level_plot/ npltlvl,pltlvl
  NAMELIST /map_plot/ ovrmap,mapgrid,latgrid,longrid,nmapfile,mapfile

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of Executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
!  Initializations
!
!-----------------------------------------------------------------------
!
  pltref=1
  pltvel=0
  pltnyq=0
  plttim=0
  npltlvl=1
  pltlvl=0
  latrad=0.
  lonrad=0.
  gridtilt_opt = 0
  numvar = 2
  pltzdr=0
  pltkdp=0
  pltrhv=0
!
!------------------------------------------------------------------------
!
! First check to see if help command line argument is specified.
!
!------------------------------------------------------------------------
!
  istatus=0
  narg=iargc()
  IF(narg > 0 ) THEN
    CALL getarg(1,charg)
    IF(charg(1:5) == '-help' .OR. charg(1:6) == '--help') THEN
      CALL getarg(0,charg)
      WRITE(6,'(a)') ' Usage: ',TRIM(charg),' [fname]'
      STOP
    END IF
  END IF
!  
!-----------------------------------------------------------------------
!
!  Get user input from input file
!
!-----------------------------------------------------------------------
!
  READ(5,file_name,END=100)
  WRITE(6,'(/a)') 'Namelist block file_name successfully read.'

  READ(5,var_plot,END=100)
  WRITE(6,'(/a)') 'Namelist block var_plot successfully read.'

  READ(5,level_plot,END=100)
  WRITE(6,'(/a)') 'Namelist block lvl_plot successfully read.'

  READ(5,map_plot,ERR=100)
  WRITE(6,'(/a)') 'Namelist block map_plot successfully read.'

  GOTO 101

  100 CONTINUE
  
  WRITE(6,'(/a,a)') 'Error reading NAMELIST file. The program will abort.'
  
  STOP

  101 CONTINUE

  IF(numvar <= 2) THEN
    pltzdr = 0; pltkdp = 0; pltrhv = 0
  ENDIF
!
!-----------------------------------------------------------------------
!
! If there is a variable in the command line, use it to replace
! the filename.
!
!-----------------------------------------------------------------------
!
  IF(narg > 0) THEN
    CALL getarg(1,fname)
    WRITE(6,'(3a)') 'Filename ',TRIM(fname),                           &
                     ' obtained from command line.'
  END IF
!
!-----------------------------------------------------------------------
!
! Read in nexrad site location table
!
!-----------------------------------------------------------------------
!
  nextfile='radarinfo.dat'
  INQUIRE(FILE=TRIM(nextfile),EXIST=fexist)
  IF(.not.fexist) THEN
    nextfile='data/adas/radarinfo.dat'
    INQUIRE(FILE=TRIM(nextfile),EXIST=fexist)
  END IF
  IF(.not.fexist) THEN
    nextfile='../data/adas/radarinfo.dat'
    INQUIRE(FILE=TRIM(nextfile),EXIST=fexist)
  END IF
  IF(fexist) CALL rdnextab(mxnextab,nnextab,nextfile,          &
                           nextstn,nextlat,nextlon,nextelv)

!-----------------------------------------------------------------------
!
! Determine file format
!
!-----------------------------------------------------------------------

  CALL asnctl ('NEWLOCAL', 1, ierr)
  CALL asnfile(fname, '-F f77 -N ieee', ierr)

  lens=LEN(trim(fname))
  IF (fname(lens-3:lens) == 'hdf4') THEN
    dmpfmt=3
  ELSE
    dmpfmt=1
  ENDIF
  write(6,'(a,i4,a,a)')                                                 &
          ' rdradcol: dmpfmt=',dmpfmt,'  fname=',TRIM(fname)

!-----------------------------------------------------------------------
!
! Get dimensions dynamically
!
!-----------------------------------------------------------------------

  IF (dmpfmt == 1) THEN

    OPEN(11,FILE=fname,FORM='unformatted',STATUS='old')
!
!-----------------------------------------------------------------------
!
! Read header data
!
!-----------------------------------------------------------------------
!
    READ(11) radid
    READ(11) ireftim,itime,vcpnum,isource,idummy,                       &
             idummy,idummy,nx,ny,nz

    READ(11) runname
    READ(11) hdmpfmt,strhopt,mapproj,irngmin,irngmax,                   &
             typelev,numradcol,nummaxlev,idummy,idummy
    READ(11) dx,dy,dz,dzmin,ctrlat,                                     &
             ctrlon,trulat1,trulat2,trulon,sclfct,                      &
             latrad,lonrad,elvrad,refelvmin,refelvmax
    READ(11) nradvr,iradvr

    READ(11) xmin, xmax, ymin, ymax

    IF(gridtilt_opt == 1) THEN
      READ(11) idummy,idummy,idummy,idummy,idummy,idummy
      READ(11) rdummy,rdummy,rdummy
      ALLOCATE (tilt_angle(nz),STAT=istatus)
      ALLOCATE (tilt_time(nz),STAT=istatus)
      READ(11) tilt_angle(1:nz)
      READ(11) tilt_time(1:nz)
    ENDIF

  ELSE
    CALL hdfopen(trim(fname), 1, sd_id)

    IF (sd_id < 0) THEN
      WRITE (6,'(3a)') 'RDRADCOL: ERROR opening hdf file:',         &
                trim(fname),' for reading.'
      istatus = 1
      STOP
    END IF

    CALL hdfrdc(sd_id, 4, 'radid', radid, istatus)
    CALL hdfrdi(sd_id, 'ireftim',  ireftim, istatus)
    CALL hdfrdi(sd_id, 'itime',    itime, istatus)
    CALL hdfrdi(sd_id, 'vcpnum',   vcpnum, istatus)
    CALL hdfrdi(sd_id, 'isource',  isource, istatus)

    CALL hdfrdc(sd_id, 40, 'runname', runname, istatus)
    CALL hdfrdi(sd_id, 'iradfmt', iradfmt, istatus)
    CALL hdfrdi(sd_id, 'strhopt', strhopt, istatus)
    CALL hdfrdi(sd_id, 'mapproj', mapproj, istatus)
    CALL hdfrdi(sd_id, 'nx', nx, istatus)
    CALL hdfrdi(sd_id, 'ny', ny, istatus)
    CALL hdfrdi(sd_id, 'nz', nz, istatus)
    CALL hdfrdr(sd_id, 'dx', dx, istatus)
    CALL hdfrdr(sd_id, 'dy', dy, istatus)
    CALL hdfrdr(sd_id, 'dz', dz, istatus)
    CALL hdfrdr(sd_id, 'dzmin', dzmin, istatus)
    CALL hdfrdr(sd_id, 'ctrlat', ctrlat, istatus)
    CALL hdfrdr(sd_id, 'ctrlon', ctrlon, istatus)
    CALL hdfrdr(sd_id, 'trulat1', trulat1, istatus)
    CALL hdfrdr(sd_id, 'trulat2', trulat2, istatus)
    CALL hdfrdr(sd_id, 'trulon', trulon, istatus)
    CALL hdfrdr(sd_id, 'sclfct', sclfct, istatus)
    CALL hdfrdr(sd_id, 'latrad', latrad, istatus)
    CALL hdfrdr(sd_id, 'lonrad', lonrad, istatus)
    CALL hdfrdr(sd_id, 'elvrad', elvrad, istatus)
    CALL hdfrdi(sd_id, 'irngmin', irngmin, istatus)
    CALL hdfrdi(sd_id, 'irngmax', irngmax, istatus)
    CALL hdfrdr(sd_id, 'refelvmin', refelvmin, istatus)
    CALL hdfrdr(sd_id, 'refelvmax', refelvmax, istatus)
    CALL hdfrdi(sd_id, 'nradvr', nradvr, istatus)
    CALL hdfrd1di(sd_id,'iradvr', mxradvr,iradvr,istatus)

    CALL hdfrdi(sd_id, 'typelev', typelev, istatus)
    CALL hdfrdr(sd_id, 'xmin', xmin, istatus)
    CALL hdfrdr(sd_id, 'xmax', xmax, istatus)
    CALL hdfrdr(sd_id, 'ymin', ymin, istatus)
    CALL hdfrdr(sd_id, 'ymax', ymax, istatus)

    CALL hdfrdi(sd_id, 'numradcol', numradcol, istatus)
    CALL hdfrdi(sd_id, 'nummaxelv', nummaxlev, istatus)
  END IF

  IF (nx < 4 .OR. ny < 4 .OR. nz < 4) THEN
    WRITE(6,'(1x,3(a,I4),a)') 'ERROR: dimension sizes are (nx,ny,nz) = ', &
                            nx,',',ny,',',nz,'.'
    WRITE(6,'(8x,2a,/2a,/)') 'The file may be in old radcol format.',   &
       ' You should use pltradcol from ','early ARPS package ',         &
       '(earlier than arps5.2.8) to read the file.'
    STOP
  END IF

!-----------------------------------------------------------------------
!
!  Allocate and Initializations
!
!-----------------------------------------------------------------------
!
  WRITE(6,'(a,3i6)') ' nx,ny,nz= ',nx,ny,nz
  ALLOCATE(gridref(nx,ny,nz),STAT=istatus)
  ALLOCATE(gridvel(nx,ny,nz),STAT=istatus)
  ALLOCATE(gridnyq(nx,ny,nz),STAT=istatus)
  ALLOCATE(gridtim(nx,ny,nz),STAT=istatus)
  gridref=-999999.
  gridvel=-999999.
  gridnyq=-999999.
  gridtim=-999999.
  IF(numvar > 2) THEN
    ALLOCATE(gridzdr(nx,ny,nz),STAT=istatus)
    ALLOCATE(gridkdp(nx,ny,nz),STAT=istatus)
    ALLOCATE(gridrhv(nx,ny,nz),STAT=istatus)
    gridzdr=-999999.
    gridkdp=-999999.
    gridrhv=-999999.
  ENDIF
  
  ALLOCATE(kntref(nz),STAT=istatus)
  ALLOCATE(kntvel(nz),STAT=istatus)
  ALLOCATE(kntzdr(nz),STAT=istatus)
  ALLOCATE(kntkdp(nz),STAT=istatus)
  ALLOCATE(kntrhv(nz),STAT=istatus)
  kntref=0
  kntvel=0
  kntzdr=0
  kntkdp=0
  kntrhv=0
!
!-----------------------------------------------------------------------
!
!  Read radar columns
!
!-----------------------------------------------------------------------
!
  IF (dmpfmt == 1) THEN

    ALLOCATE(readk(nz) ,STAT=istatus)
    ALLOCATE(readhgt(nz),STAT=istatus)
    IF(gridtilt_opt == 1) ALLOCATE(readrng(nz),STAT=istatus)
    ALLOCATE(readref(nz),STAT=istatus)
    ALLOCATE(readvel(nz),STAT=istatus)
    ALLOCATE(readnyq(nz),STAT=istatus)
    ALLOCATE(readtim(nz),STAT=istatus)
    IF(numvar > 2) THEN
      ALLOCATE(readpol(nz),STAT=istatus)
      ALLOCATE(readzdr(nz),STAT=istatus)
      ALLOCATE(readkdp(nz),STAT=istatus)
      ALLOCATE(readrhv(nz),STAT=istatus)
      readpol=0
      readzdr=0
      readkdp=0
      readrhv=0
    ENDIF
    readk=0
    readhgt=0
    IF(gridtilt_opt == 1) readrng=0
    readref=0
    readvel=0
    readnyq=0
    readtim=0

    DO ipt=1,numradcol
      READ(11,END=51) i,j,xrd,yrd,                                        &
                     gridlat,gridlon,elev,klev
      IF(gridtilt_opt == 1) THEN
        READ(11,END=52) (readhgt(kk),kk=1,klev)
        READ(11,END=52) (readrng(kk),kk=1,klev)
        READ(11,END=52) (readvel(kk),kk=1,klev)
        READ(11,END=52) (readref(kk),kk=1,klev)
        IF(numvar > 2) THEN
          READ(11,END=52) (readpol(kk),kk=1,klev)
          READ(11,END=52) (readzdr(kk),kk=1,klev)
          READ(11,END=52) (readpol(kk),kk=1,klev)
          READ(11,END=52) (readkdp(kk),kk=1,klev)
          READ(11,END=52) (readrhv(kk),kk=1,klev)
        ENDIF
      ELSE
        READ(11,END=52) (readk(kk),kk=1,klev)
        READ(11,END=52) (readhgt(kk),kk=1,klev)
        READ(11,END=52) (readref(kk),kk=1,klev)
        READ(11,END=52) (readvel(kk),kk=1,klev)
        READ(11,END=52) (readnyq(kk),kk=1,klev)
        READ(11,END=52) (readtim(kk),kk=1,klev)
      ENDIF
      DO kk=1,klev
        IF(gridtilt_opt == 1) THEN
          k = kk
        ELSE
          k=readk(kk)
          gridnyq(i,j,k)=readnyq(kk)
          gridtim(i,j,k)=readtim(kk)
        ENDIF
        gridref(i,j,k)=readref(kk)
        gridvel(i,j,k)=readvel(kk)
        IF (numvar > 2) THEN
          gridzdr(i,j,k)=readzdr(kk)
          gridkdp(i,j,k)=readkdp(kk)
          gridrhv(i,j,k)=readrhv(kk)
        ENDIF

        IF(gridref(i,j,k) > -200. .AND. gridref(i,j,k) < 200.)  &
           kntref(k)=kntref(k)+1
        IF(gridvel(i,j,k) > -200. .AND. gridvel(i,j,k) < 200.)  &
           kntvel(k)=kntvel(k)+1

        IF (numvar > 2) THEN
          IF(gridzdr(i,j,k) > -200. .AND. gridzdr(i,j,k) < 200.)  &
            kntzdr(k)=kntzdr(k)+1
          IF(gridkdp(i,j,k) > -200. .AND. gridkdp(i,j,k) < 200.)  &
            kntkdp(k)=kntkdp(k)+1
          IF(gridrhv(i,j,k) > -200. .AND. gridrhv(i,j,k) < 200.)  &
            kntrhv(k)=kntrhv(k)+1
        ENDIF

      END DO
    END DO
    51 CONTINUE
    ipt=ipt-1
    WRITE(6,'(a,i6,a)') ' End of file reached after reading',             &
                       ipt,' columns'
    GO TO 55
    52 CONTINUE
    WRITE(6,'(a,i6,a)') ' End of file reached while reading',             &
                       ipt,' column'
    55 CONTINUE
    CLOSE(11)
  ELSE
   !
   !  Allocate hdf temporary arrays
   !
    ALLOCATE (itmp(nummaxlev,numradcol),stat=istatus)
    CALL check_alloc_status(istatus, "pltradcol:itmp")
    
    ALLOCATE(coli(numradcol), STAT = istatus)
    CALL check_alloc_status(istatus, "pltradcol:coli")

    ALLOCATE(colj(numradcol), STAT = istatus)
    CALL check_alloc_status(istatus, "pltradcol:colj")

    ALLOCATE(colk(nummaxlev,numradcol), STAT = istatus)
    CALL check_alloc_status(istatus, "pltradcol:colk")

    ALLOCATE(numlev(numradcol), STAT = istatus)
    CALL check_alloc_status(istatus, "pltradcol:numlev")

    ALLOCATE(radcolref(nummaxlev,numradcol), STAT = istatus)
    CALL check_alloc_status(istatus, "pltradcol:radcolref")

    ALLOCATE(radcolvel(nummaxlev,numradcol), STAT = istatus)
    CALL check_alloc_status(istatus, "pltradcol:radcolvel")

    ALLOCATE(radcolnyq(nummaxlev,numradcol), STAT = istatus)
    CALL check_alloc_status(istatus, "pltradcol:radcolnyq")

    ALLOCATE(radcoltim(nummaxlev,numradcol), STAT = istatus)
    CALL check_alloc_status(istatus, "pltradcol:radcoltim")

    CALL hdfrd1di(sd_id,'numelev',  numradcol,numlev,istatus)
    
    CALL hdfrd1di(sd_id,'radcoli',  numradcol,coli,istatus)
    CALL hdfrd1di(sd_id,'radcolj',  numradcol,colj,istatus)
    CALL hdfrd2di(sd_id,'radcolk',  nummaxlev,numradcol,colk,istatus)

    CALL hdfrd2d (sd_id,'radcolref',nummaxlev,numradcol,radcolref,istatus,itmp)
    CALL hdfrd2d (sd_id,'radcolvel',nummaxlev,numradcol,radcolvel,istatus,itmp)
    CALL hdfrd2d (sd_id,'radcolnyq',nummaxlev,numradcol,radcolnyq,istatus,itmp)
    CALL hdfrd2d (sd_id,'radcoltim',nummaxlev,numradcol,radcoltim,istatus,itmp)

    CALL hdfclose(sd_id,istatus)

    DO ipt=1,numradcol
      i = coli(ipt)
      j = colj(ipt)
      DO kk=1,numlev(ipt)
        k = colk(kk,ipt)

        gridref(i,j,k)=radcolref(kk,ipt)
        gridvel(i,j,k)=radcolvel(kk,ipt)
        gridnyq(i,j,k)=radcolnyq(kk,ipt)
        gridtim(i,j,k)=radcoltim(kk,ipt)

        IF(gridref(i,j,k) > -200. .AND. gridref(i,j,k) < 200.)  &
            kntref(k)=kntref(k)+1
        IF(gridvel(i,j,k) > -200. .AND. gridvel(i,j,k) < 200.)  &
            kntvel(k)=kntvel(k)+1

      END DO
    END DO
  END IF
!
!-----------------------------------------------------------------------
!
!  Write statistics
!
!-----------------------------------------------------------------------
!
  WRITE(6,'(a)') '  k     n ref     n vel'
  DO k = 1, nz
    WRITE(6,'(i3,2i10)') k, kntref(k), kntvel(k)
  END DO
!
!-----------------------------------------------------------------------
!
! Match name to table
!
!-----------------------------------------------------------------------
!
  found=.TRUE.
  IF (latrad == 0. .AND. lonrad == 0.) THEN
    found = .FALSE.
    IF(fexist) THEN
      DO itab=1,nnextab
        IF(nextstn(itab)(1:4) == radid) THEN
          latrad=nextlat(itab)
          lonrad=nextlon(itab)
          found=.TRUE.
          EXIT
        END IF
      END DO
    END IF
    IF(found) THEN
      WRITE(6,'(a,a,a,f10.4,a,f10.4)')                                 &
        'Radar ',radid,' lat:',latrad,' lon:',lonrad
    ELSE
      WRITE(6,'(a,a)') ' Could not find match for ',radid
    END IF
  END IF
!
  CALL abss2ctim(itime, iyr, imon, idy, ihr, imin, isec )
  iyr=MOD(iyr,100)
  WRITE(timplt,815) imon,idy,iyr,ihr,imin
  815 FORMAT(i2.2,'/',i2.2,'/',i2.2,1X,i2.2,':',i2.2,' UTC')
!
!-----------------------------------------------------------------------
!
!  Set up map
!
!-----------------------------------------------------------------------
!
  IF(ctrlon > 180.) ctrlon=ctrlon-360.
  IF(trulon > 180.) trulon=trulon-360.
  latnot(1)=trulat1
  latnot(2)=trulat2
  CALL setmapr(mapproj,sclfct,latnot,trulon)
  CALL lltoxy(1,1,ctrlat,ctrlon,ctrx,ctry)
  swx=ctrx-(0.5*(nx-3)*dx)
  swy=ctry-(0.5*(ny-3)*dy)
  CALL setorig(1,swx,swy)
!
!-----------------------------------------------------------------------
!
!  Initialize ZXPLOT plotting package
!
!-----------------------------------------------------------------------
!
  CALL xdevic_new(1,outfilename,0,0)
  CALL xctrbadv( 1 )  ! Turn on missing value checking for contouring
  CALL xvtrbadv( 1 )  ! Turn on missing value checking for vector plotting
  CALL xsetclrs_new(6,0)
  CALL xaxfmt('(I6)')
  CALL xbadval ( misval ) ! Set the missing value flag to misval
  CALL xcolor(1)
  CALL xintmkr()

  dxkm = dx*0.001
  dykm = dy*0.001
  hlfdxkm = 0.5*dxkm
  hlfdykm = 0.5*dykm
!
! For a grid much larger than the diameter of the radar coverage,
! display only the radar coverage area, centered on the radar.
!
  xorig = 0.0
  yorig = 0.0
  xl = (nx-3)*dxkm
  yl = (ny-3)*dykm
  IF( irngmax == 0) irngmax=230000
  xlenchk=0.004*irngmax
!
! Setup device and define plotting space.
!
  CALL xstpjgrd(mapproj,trulat1,trulat2,trulon,                         &
                ctrlat,ctrlon,xl,yl,xorig,yorig)
!
  CALL xlltoxy(1,1,latrad,lonrad,xstn,ystn)
  xstn=xstn*0.001
  ystn=ystn*0.001

  IF( xl > xlenchk .OR. yl > xlenchk) THEN
    WRITE(6,'(a)')  'Large grid, setting radar-centric domain'
    WRITE(6,'(a,f10.1,a,f10.1)') ' xl: ',xl,'  rngmax: ',(0.001*irngmax)
    xl = 0.002*irngmax
    yl = xl
    xstnorg=xstn
    ystnorg=ystn
    CALL xstpjgrd(mapproj,trulat1,trulat2,trulon,                       &
                  latrad,lonrad,xl,yl,xorig,yorig)
    CALL xlltoxy(1,1,latrad,lonrad,xstn,ystn)
    xstn=xstn*0.001
    ystn=ystn*0.001
    xoffset=xstnorg-xstn
    yoffset=ystnorg-ystn
  ELSE
    WRITE(6,'(a)')  'Small grid, setting grid-centric domain'
    WRITE(6,'(2(a,f12.1))') ' xl:',xl,'  rngmax:',(0.001*irngmax)
    xoffset=0.0
    yoffset=0.0
  END IF
  print *, ' xoffset: ',xoffset
  print *, ' yoffset: ',yoffset
  CALL xmap(0.0,xl, 0.0, yl)
  chsize=yl*0.01
  CALL xmrksz(0.002)
  yloc=1.03*yl 
  labsize=0.03*yl
  mrksize=dxkm*0.00015
!
  IF( xl >= yl) THEN
    CALL xpspac(0.1,0.9, 0.5-yl/xl*0.4,0.5+yl/xl*0.4)
  ELSE
    CALL xpspac(0.5-xl/yl*0.4,0.5+xl/yl*0.4,0.1, 0.9)
  END IF
  xname=xstn-xoffset
  yname=ystn-(yl*0.04)-yoffset
!
!-----------------------------------------------------------------------
!
!  Vertical level loop
!
!-----------------------------------------------------------------------
!
  DO lindex=1,npltlvl
    IF(pltlvl(lindex) < 1 .OR. pltlvl(lindex) > nz) CYCLE
    k=pltlvl(lindex)
    PRINT *, ' plotting k= ',k
    WRITE(lvlid,'(a,i3)') 'k=',k
!
!-----------------------------------------------------------------------
!
!  Plot reflectivity
!
!-----------------------------------------------------------------------
!
    IF( pltref > 0 .AND. kntref(k) > 0 ) THEN
      IF (nextframe) CALL xframe
      CALL xcolor(1)
      CALL xchsiz(labsize)
      CALL xwdwof
      CALL xcharc((xorig+(0.1*xl)),yloc,radid)
      CALL xcharc((xorig+(0.3*xl)),yloc,lvlid)
      CALL xcharc((xorig+(0.6*xl)),yloc,timplt)
      CALL xcharc((xorig+(0.9*xl)),yloc,'Reflectivity')
      CALL xchsiz(chsize)
      IF(pltref == 1) THEN
        CALL xcolor(1)
        DO j=1,ny
          DO i=1,nx
            IF(gridref(i,j,k) > -200. .AND. gridref(i,j,k) < 200.) THEN
              WRITE(chplt,'(i3)') nint(gridref(i,j,k))
              x=(((i-1)*dxkm)-hlfdxkm)-xoffset
              y=(((j-1)*dykm)-hlfdykm)-yoffset
              CALL xcharc(x,y,chplt)
            END IF
          END DO
        END DO
      ELSE IF(pltref == 2) THEN
        DO i=8,23
          CALL xcolor(i)
          WRITE(chplt,'(i3)') ((i-9)*5)
          CALL xcharc((0.03*i*xl),(yl*1.07),chplt)
        END DO
        CALL xcolor(1)
        jcol=1
        DO j=1,ny
          DO i=1,nx
            IF(gridref(i,j,k) > -200. .AND. gridref(i,j,k) < 200.) THEN
              WRITE(chplt,'(i3)') nint(gridref(i,j,k))
              icol=nint(0.2*gridref(i,j,k))+9
              icol=min(max(icol,8),23)
              IF(icol /= jcol) call xcolor(icol) 
              jcol=icol
              x=(((i-1)*dxkm)-hlfdxkm)-xoffset
              y=(((j-1)*dykm)-hlfdykm)-yoffset
!             print *, ' plotting at x,y:',x,y
              CALL xcharc(x,y,chplt)
            END IF
          END DO
        END DO
      ELSE
        CALL xmrksz(mrksize)
        DO i=8,23
          CALL xcolor(i)
          WRITE(chplt,'(i3)') ((i-9)*5)
          CALL xcharc((0.03*i*xl),(yl*1.07),chplt)
        END DO
        CALL xcolor(1)
        jcol=1
        DO j=1,ny
          DO i=1,nx
            IF(gridref(i,j,k) > -200. .AND. gridref(i,j,k) < 200.) THEN
              icol=nint(0.2*gridref(i,j,k))+9
              icol=min(max(icol,8),23)
              IF(icol /= jcol) call xcolor(icol) 
              jcol=icol
              x=(((i-1)*dxkm)-hlfdxkm)-xoffset
              y=(((j-1)*dykm)-hlfdykm)-yoffset
              CALL xmarker(x,y,9)
            END IF
          END DO
        END DO
      END IF
      CALL xcolor(1)
      CALL xmrksz(0.002)
      IF(found) THEN
        CALL xmarker(xstn,ystn,9)
        CALL xcharc(xname,yname,radid)
      END IF
      CALL xaxsca(0.0,xl, 0.001*dx, 0.0, yl, 0.001*dy )
      CALL pltmapx(mxmapfile,nmapfile,mapfile,latgrid,longrid)
      nextframe = .TRUE.
    END IF
!
!-----------------------------------------------------------------------
!
!  Plot velocity
!
!-----------------------------------------------------------------------
!
    IF( pltvel > 0 .AND. kntvel(k) > 0 ) THEN
      IF (nextframe) CALL xframe
      CALL xcolor(1)
      CALL xchsiz(labsize)
      CALL xwdwof
      CALL xcharc((0.1*xl),yloc,radid)
      CALL xcharc((0.3*xl),yloc,lvlid)
      CALL xcharc((0.6*xl),yloc,timplt)
      CALL xcharc((0.9*xl),yloc,'Radial Velocity')
      CALL xchsiz(chsize)
      IF(pltvel == 1) THEN
        CALL xcolor(1)
        DO j=1,ny
          DO i=1,nx
            IF(gridvel(i,j,k) > -200. .AND. gridvel(i,j,k) < 200.) THEN
              WRITE(chplt,'(i3)') nint(gridvel(i,j,k))
              x=(((i-1)*dxkm)-hlfdxkm)-xoffset
              y=(((j-1)*dykm)-hlfdykm)-yoffset
              CALL xcharc(x,y,chplt)
            END IF
          END DO
        END DO
      ELSE IF( pltvel == 2 ) THEN
        DO i=9,23
          CALL xcolor(i)
          WRITE(chplt,'(i3)') ((i-16)*5)
          CALL xcharc((0.03*i*xl),(yl*1.07),chplt)
        END DO
        CALL xcolor(9)
        jcol=9
        DO j=1,ny
          DO i=1,nx
            IF(gridvel(i,j,k) > -200. .AND. gridvel(i,j,k) < 200.) THEN
              WRITE(chplt,'(i3)') nint(gridvel(i,j,k))
              icol=nint(0.2*gridvel(i,j,k))+16
              icol=min(max(icol,9),23)
              IF(icol /= jcol) call xcolor(icol) 
              jcol=icol
              x=(((i-1)*dxkm)-hlfdxkm)-xoffset
              y=(((j-1)*dykm)-hlfdykm)-yoffset
              CALL xcharc(x,y,chplt)
            END IF
          END DO
        END DO
      ELSE
        CALL xmrksz(mrksize)
        DO i=9,23
          CALL xcolor(i)
          WRITE(chplt,'(i3)') ((i-16)*5)
          CALL xcharc((0.03*i*xl),(yl*1.07),chplt)
        END DO
        CALL xcolor(9)
        jcol=9
        DO j=1,ny
          DO i=1,nx
            IF(gridvel(i,j,k) > -200. .AND. gridvel(i,j,k) < 200.) THEN
              icol=nint(0.2*gridvel(i,j,k))+16
              icol=min(max(icol,9),23)
              IF(icol /= jcol) call xcolor(icol) 
              jcol=icol
              x=(((i-1)*dxkm)-hlfdxkm)-xoffset
              y=(((j-1)*dykm)-hlfdykm)-yoffset
              CALL xmarker(xstn,ystn,9)
            END IF
          END DO
        END DO
      END IF
      CALL xcolor(1)
      CALL xmrksz(0.002)
      IF(found) THEN
        CALL xmarker(xstn,ystn,1)
        CALL xcharc(xname,yname,radid)
      END IF
      CALL xaxsca(0.0,xl, 0.001*dx, 0.0, yl, 0.001*dy )
      CALL pltmapx(mxmapfile,nmapfile,mapfile,latgrid,longrid)
      nextframe = .TRUE.
    END IF
!
!-----------------------------------------------------------------------
!
!  Plot Nyquist velocity
!
!-----------------------------------------------------------------------
!
    IF(pltnyq > 0 .AND. kntvel(k) > 0) THEN
      IF (nextframe) CALL xframe
      CALL xcolor(1)
      CALL xchsiz(labsize)
      CALL xwdwof
      CALL xcharc((0.1*xl),yloc,radid)
      CALL xcharc((0.3*xl),yloc,lvlid)
      CALL xcharc((0.6*xl),yloc,timplt)
      CALL xcharc((0.9*xl),yloc,'Nyquist')
      CALL xcolor(15)
      CALL xchsiz(chsize)
      DO j=1,ny
        DO i=1,nx
          IF(gridnyq(i,j,k) > 0. .AND. gridnyq(i,j,k) < 900.) THEN
            WRITE(chplt,'(i3)') nint(10.*gridnyq(i,j,k))
            x=(((i-1)*dxkm)-hlfdxkm)-xoffset
            y=(((j-1)*dykm)-hlfdykm)-yoffset
            CALL xcharc(x,y,chplt)
          END IF
        END DO
      END DO
      CALL xcolor(1)
      CALL xmrksz(0.002)
      IF(found) THEN
        CALL xmarker(xstn,ystn,1)
        CALL xcharc(xname,yname,radid)
      END IF
      CALL xaxsca(0.0,xl, 0.001*dx, 0.0, yl, 0.001*dy )
      CALL pltmapx(mxmapfile,nmapfile,mapfile,latgrid,longrid)
      nextframe = .TRUE.
    END IF
!
!-----------------------------------------------------------------------
!
!  Plot data time
!
!-----------------------------------------------------------------------
!
    IF(plttim > 0 .AND. (kntvel(k) > 0 .OR. kntref(k) > 0)) THEN
      IF (nextframe) CALL xframe
      CALL xcolor(1)
      CALL xchsiz(labsize)
      CALL xwdwof
      CALL xcharc((0.1*xl),yloc,radid)
      CALL xcharc((0.3*xl),yloc,lvlid)
      CALL xcharc((0.6*xl),yloc,timplt)
      CALL xcharc((0.9*xl),yloc,'Time (secs/10)')
      CALL xcolor(15)
      CALL xchsiz(chsize)
      DO j=1,ny
        DO i=1,nx
          IF(gridtim(i,j,k) > -990. .AND. gridtim(i,j,k) < 9999.) THEN
            WRITE(chplt,'(i3)') nint(0.1*gridtim(i,j,k))
            x=(((i-1)*dxkm)-hlfdxkm)-xoffset
            y=(((j-1)*dykm)-hlfdykm)-yoffset
            CALL xcharc(x,y,chplt)
          END IF
        END DO
      END DO
      CALL xcolor(1)
      CALL xmrksz(0.002)
      IF(found) THEN
        CALL xmarker(xstn,ystn,1)
        CALL xcharc(xname,yname,radid)
      END IF
      CALL xaxsca(0.0,xl, 0.001*dx, 0.0, yl, 0.001*dy )
      CALL pltmapx(mxmapfile,nmapfile,mapfile,latgrid,longrid)
      nextframe = .TRUE.
    END IF
!
!-----------------------------------------------------------------------
!
!  Plot differential reflectivity
!
!-----------------------------------------------------------------------
!
    IF( pltzdr > 0 .AND. kntzdr(k) > 0 ) THEN
      IF (nextframe) CALL xframe
      CALL xcolor(1)
      CALL xchsiz(labsize)
      CALL xwdwof
      CALL xcharc((0.1*xl),yloc,radid)
      CALL xcharc((0.3*xl),yloc,lvlid)
      CALL xcharc((0.6*xl),yloc,timplt)
      CALL xcharc((0.9*xl),yloc,'Zdr')
      CALL xchsiz(chsize)
      IF(pltzdr == 1) THEN
        CALL xcolor(1)
        DO j=1,ny
          DO i=1,nx
            IF(gridzdr(i,j,k) > -200. .AND. gridzdr(i,j,k) < 200.) THEN
              WRITE(chplt,'(i3)') nint(gridzdr(i,j,k))
              x=(((i-1)*dxkm)-hlfdxkm)-xoffset
              y=(((j-1)*dykm)-hlfdykm)-yoffset
              CALL xcharc(x,y,chplt)
            END IF
          END DO
        END DO
      ELSE IF(pltzdr == 2) THEN
        DO i=8,23
          CALL xcolor(i)
          WRITE(chplt,'(i3)') ((i-9)*5)
          CALL xcharc((0.03*i*xl),(yl*1.07),chplt)
        END DO
        CALL xcolor(1)
        jcol=1
        DO j=1,ny
          DO i=1,nx
            IF(gridzdr(i,j,k) > -200. .AND. gridzdr(i,j,k) < 200.) THEN
              WRITE(chplt,'(i3)') nint(gridzdr(i,j,k))
              icol=nint(0.2*gridzdr(i,j,k))+9
              icol=min(max(icol,8),23)
              IF(icol /= jcol) call xcolor(icol) 
              jcol=icol
              x=(((i-1)*dxkm)-hlfdxkm)-xoffset
              y=(((j-1)*dykm)-hlfdykm)-yoffset
              CALL xcharc(x,y,chplt)
            END IF
          END DO
        END DO
      ELSE
        CALL xmrksz(mrksize)
        DO i=8,23
          CALL xcolor(i)
          WRITE(chplt,'(i3)') ((i-9)*5)
          CALL xcharc((0.03*i*xl),(yl*1.07),chplt)
        END DO
        CALL xcolor(1)
        jcol=1
        DO j=1,ny
          DO i=1,nx
            IF(gridzdr(i,j,k) > -200. .AND. gridzdr(i,j,k) < 200.) THEN
              icol=nint(0.2*gridzdr(i,j,k))+9
              icol=min(max(icol,8),23)
              IF(icol /= jcol) call xcolor(icol) 
              jcol=icol
              x=(((i-1)*dxkm)-hlfdxkm)-xoffset
              y=(((j-1)*dykm)-hlfdykm)-yoffset
              CALL xmarker(x,y,9)
            END IF
          END DO
        END DO
      END IF
      CALL xcolor(1)
      CALL xmrksz(0.002)
      IF(found) THEN
        CALL xmarker(xstn,ystn,9)
        CALL xcharc(xname,yname,radid)
      END IF
      CALL xaxsca(0.0,xl, 0.001*dx, 0.0, yl, 0.001*dy )
      CALL pltmapx(mxmapfile,nmapfile,mapfile,latgrid,longrid)
      nextframe = .TRUE.
    END IF

  END DO

  CALL xgrend

  WRITE(6,'(/,a,/)') '   ==== PLTRADCOL: Normal Termination ===='

  STOP
END PROGRAM pltradcol

SUBROUTINE pltmapx(maxmapfile,nmapfile,mapfile,latgrid,longrid)
  IMPLICIT NONE
  INTEGER :: maxmapfile
  INTEGER :: nmapfile
  CHARACTER (LEN=256) :: mapfile(maxmapfile)
  REAL :: latgrid
  REAL :: longrid
!
! Misc Local Variables
!
  INTEGER :: i,lmapfile
  LOGICAL :: fexist

  DO i=1,nmapfile

    lmapfile=LEN(mapfile(i))
    CALL strlnth(mapfile(i), lmapfile)

    INQUIRE(FILE=mapfile(i)(1:lmapfile), EXIST = fexist )
    IF( .NOT.fexist) THEN
      WRITE(6,'(a)') 'Map file '//mapfile(i)(1:lmapfile)//              &
          ' not found. Program stopped in PLTMAPX.'
      STOP 001
    END IF

    CALL xdrawmap(10,mapfile(i)(1:lmapfile),latgrid,longrid)

  END DO

  RETURN 
END SUBROUTINE pltmapx

SUBROUTINE rdnextab(mxnextab,nnextab,nextfile,nextstn,nextlat,nextlon,nextelv)
  IMPLICIT NONE
  INTEGER mxnextab
  INTEGER :: nnextab
  CHARACTER (LEN=*) :: nextfile
  CHARACTER (LEN=4) :: nextstn(mxnextab)
  REAL :: nextlat(mxnextab)
  REAL :: nextlon(mxnextab)
  REAL :: nextelv(mxnextab)
!
  INTEGER :: latdeg,latmin,latsec,londeg,lonmin,lonsec,ielev
  INTEGER :: istn,istat
  CHARACTER (LEN=1) :: dummy
!
  OPEN(41,file=TRIM(nextfile),status='old')
  READ(41,'(a)') dummy
  DO istn=1,mxnextab
    READ(41,'(a4,21x,i3,6x,i2,6x,i2,7x,i3,5x,i2,6x,i2,3x,i7)',iostat=istat) &
         nextstn(istn),latdeg,latmin,latsec,londeg,lonmin,lonsec,ielev
    IF(istat /= 0) EXIT
    nextlat(istn)=float(latdeg)+(float(latmin)/60.)+(float(latsec)/3600.)
    nextlon(istn)=-(float(londeg)+(float(lonmin)/60.)+(float(lonsec)/3600.))
    nextelv(istn)=float(ielev)
  END DO
  nnextab=istn-1
  close(41)
END SUBROUTINE rdnextab
