!--------------------------------------------------------------------
!  ppiplt.f90
!--------------------------------------------------------------------
!
!  Purpose: plot radar reflectiviy for each scan
!
!--------------------------------------------------------------------
!
! AUTHOR:
!
! Mingjing Tong (9/21/2005)
!
! MODIFICATIONS:
!
! 12 Jan 2012  Youngsun Jung
! Reformatted the program and updated.
!
!------------------------------------------------------------------------

MODULE VAR_DEF
!--------------------------------------------------------------------
!
!  Declaration of variables
!
!--------------------------------------------------------------------
  IMPLICIT NONE
  SAVE

  INTEGER, PARAMETER :: nprio = 4
  INTEGER :: nx, ny, nz, ntilt, nplt
  INTEGER :: i, j, itilt, l, ii
  INTEGER :: nen3, itime
  INTEGER :: istatus
  REAL, PARAMETER :: mv = -111.0
  REAL, allocatable :: x(:)
  REAL, allocatable :: y(:)
  REAL, allocatable :: xs(:)
  REAL, allocatable :: ys(:)
  REAL, allocatable :: zps(:,:,:)
  REAL, allocatable :: azmsc     (:,:)
  REAL, allocatable :: elvsc     (:,:,:)
  REAL, allocatable :: rngsc     (:,:,:)
  REAL, allocatable :: elvobs(:)
  REAL, allocatable :: hgtobs(:,:,:)    ! high of observation point (Vr)
  REAL, allocatable :: rngoelv(:,:,:)
  REAL, allocatable :: vroelv(:,:,:)    ! Vr obs
  REAL, allocatable :: vrotmp(:,:,:)
  REAL, allocatable :: refotmp(:,:,:)
  REAL, allocatable :: vatmp(:,:,:)
  REAL, allocatable :: refoelv(:,:,:)   ! Global Z observation
  REAL, allocatable :: vphie(:,:,:)     ! Cross-beam wind component
  REAL, allocatable :: var3d(:,:,:)
  REAL :: rdrlat, rdrlon, radarx, radary, rdralt, dazim, rngmin, rngmax
  REAL :: time
  INTEGER :: timeset,iyr,imon,iday,ihour,imin,isec

  INTEGER :: ncyc
  CHARACTER (LEN=10) :: obsnum, timnam
  CHARACTER (LEN=180) :: string, string2
  CHARACTER(LEN=180) :: filename
  CHARACTER (LEN=4) :: rdrnam

  LOGICAL :: fexist
  INTEGER :: lmapfile

  INTEGER :: M, N
  REAL :: dx, dy
  REAL :: xscale, yscale
  REAL :: swx, swy
  REAL :: xl,xr,yb,yt
  REAL :: x1, x2, y1, y2
  REAL :: xbgn, xend, ybgn, yend
  REAL, allocatable :: xc(:,:), yc(:,:)
!
  INTEGER, PARAMETER :: maxnum = 130
  INTEGER :: nfile, nn
  INTEGER :: datanamopt, datafmt
  CHARACTER (LEN=180) :: rdrdtadir, vadtadir, vardtadir
  CHARACTER (LEN=100) :: rdrdtanam(maxnum), vadtanam(maxnum), vardtanam(maxnum)
  CHARACTER (LEN=30) :: radarnam
  CHARACTER (LEN=50) :: vanam, varnam
  INTEGER :: t0, tn, dt

  CHARACTER (LEN=80) :: psfile

  NAMELIST /output/ psfile

  REAL :: dxr, dyr
  INTEGER :: coordopt

  NAMELIST /grid/ dxr, dyr, coordopt

  INTEGER :: ctrlbopt  ! Contour labelling option
  INTEGER :: ctrstyle
  INTEGER :: ctrlbfrq
  INTEGER :: lbmaskopt
  REAL :: ctrlbsiz

  CHARACTER(LEN=256) :: coltabfn_ref, coltabfn_vr
  CHARACTER(LEN=256) :: coltabfn_va, coltabfn_var
  INTEGER :: refpltopt, vrpltopt, vapltopt, varpltopt
  INTEGER :: coltab_ref, ibgncol_ref, iendcol_ref
  INTEGER :: coltab_vr, ibgncol_vr, iendcol_vr
  INTEGER :: coltab_va, ibgncol_va, iendcol_va
  INTEGER :: coltab_var, ibgncol_var, iendcol_var
  INTEGER :: tiltbgn, tiltend
  REAL :: dxscale, dyscale
  REAL :: refminc, refmaxc, refinc
  REAL :: vrminc, vrmaxc, vrinc
  REAL :: vaminc, vamaxc, vainc
  REAL :: varminc, varmaxc, varinc
  INTEGER :: refhlf, vrhlf, vahlf, varhlf
  INTEGER :: refovr, vrovr, vaovr, varovr
  INTEGER :: refzro, vrzro, vazro, varzro
  INTEGER :: refprio, vrprio, vaprio, varprio

  INTEGER, PARAMETER :: mxmapfile= 10
  INTEGER :: mapproj
  REAL :: trulat1
  REAL :: trulat2
  REAL :: trulon
  REAL :: ctrlat
  REAL :: ctrlon
  REAL :: xorig,yorig

  INTEGER :: ovrmap
  INTEGER :: mapgrid
  REAL :: latgrid
  REAL :: longrid
  INTEGER :: nmapfile
  CHARACTER (LEN=180) :: mapfile(mxmapfile)

  INTEGER :: first_fram, ip, k, ovrid
  INTEGER, DIMENSION(nprio) :: iprio, iovr, iplt

END MODULE VAR_DEF

PROGRAM ppiplt
!
!--------------------------------------------------------------------
!
!  Purpose: plot radar reflectiviy for each scan
!
!--------------------------------------------------------------------
!
  USE VAR_DEF
  IMPLICIT NONE

  NAMELIST /input_data/ datanamopt, datafmt, rdrdtadir, vadtadir, vardtadir,    &
                        nfile, rdrdtanam, vadtanam, vardtanam,                  &
                        radarnam, vanam, varnam, t0, tn, dt

  NAMELIST /plot_setup/ ctrlbopt,ctrstyle,ctrlbfrq,ctrlbsiz,lbmaskopt,    &
                        refpltopt,coltab_ref,coltabfn_ref,ibgncol_ref,    &
                        iendcol_ref,refminc,refmaxc,refinc,refovr,refhlf, &
                        refzro,refprio,vrpltopt,coltab_vr,coltabfn_vr,    &
                        ibgncol_vr,iendcol_vr,vrminc,vrmaxc,vrinc,vrovr,  &
                        vrhlf,vrzro,vrprio,vapltopt,coltab_va,coltabfn_va,&
                        ibgncol_va,iendcol_va,vaminc,vamaxc,vainc,vaovr,  &
                        vahlf,vazro,vaprio,varpltopt,coltab_var,coltabfn_var,&
                        ibgncol_var,iendcol_var,varminc,varmaxc,varinc,   &
                        varovr,varhlf,varzro,varprio,tiltbgn,tiltend,x1,  &
                        x2,y1,xbgn,xend,ybgn,yend,dxscale,dyscale

  NAMELIST /projection/ mapproj, trulat1, trulat2, trulon, ctrlat, ctrlon, xorig, yorig

  NAMELIST /map_plot/ ovrmap,mapgrid,latgrid,longrid,nmapfile,mapfile

  INTEGER, PARAMETER :: unum=5          ! unit number for reading in namelist
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  Now we begin to read in the values of parameters:
!
!-----------------------------------------------------------------------
!
  READ(unum,input_data,END=100)
  WRITE(6,'(a)') 'Namelist input_data was sucessfully read in.'

  READ(unum,output,END=100)
  WRITE(6,'(a)') 'Namelist output was sucessfully read in.'

  READ(unum,grid,END=100)
  WRITE(6,'(a)') 'Namelist grid was sucessfully read in.'

  READ (unum,plot_setup,END=100)
  WRITE(6,'(a)')'Namelist block  plot_setup sucessfully read.'

  READ (unum,projection,END=100)
  WRITE(6,'(a)') 'Namelist projection was sucessfully read in.'

  READ(unum,map_plot,END=100)
  WRITE(6,'(a)') 'Namelist map_plot sucessfully read in.'

!-----------------------------------------------------------------------
!
! READ IN DATA
!
!-----------------------------------------------------------------------
  obsnum='xx'
  timnam = 'xxxxxx'
  IF(datanamopt == 1)then
    ncyc=INT((tn-t0)/dt)+1
    write(*,*)'ncyc',ncyc
    nn = ncyc
  ELSE
    nn = nfile
  ENDIF
  DO l=1,nn
    IF(datanamopt == 1)then
      itime=INT(t0+(l-1)*dt)
      time=1.0*itime
      timnam(1:1)=ACHAR( 48+int(itime/100000) )
      nen3=itime-100000*int(itime/100000)
      timnam(2:2)=ACHAR( 48+int(nen3/10000) )
      nen3=nen3-10000*int(nen3/10000)
      timnam(3:3)=ACHAR( 48+int(nen3/1000) )
      nen3=nen3-1000*int(nen3/1000)
      timnam(4:4)=ACHAR( 48+int(nen3/100) )
      nen3=nen3-100*int(nen3/100)
      timnam(5:5)=ACHAR( 48+int(nen3/10) )
      nen3=nen3-10*int(nen3/10)
      timnam(6:6)=ACHAR( 48+nen3 )
    ENDIF

    IF(datafmt == 1)THEN
      IF(refpltopt > 0 .or. vrpltopt > 0)THEN
        IF(datanamopt == 1)then
          filename=trim(rdrdtadir)//'/'//trim(radarnam)//trim(timnam)
          write(*,*)trim(filename)
        ELSE
          filename=trim(rdrdtadir)//'/'//trim(rdrdtanam(l))
          write(*,*)trim(filename)
        ENDIF
        OPEN(UNIT=20,FILE=trim(filename),STATUS='old', &
             FORM='unformatted',IOSTAT= istatus )
        IF( istatus /= 0) GO TO 999
        READ(20,ERR=110,END=120) timeset,iyr,imon,iday,ihour,imin,isec
        READ(20,ERR=110,END=120) ntilt,nx,ny
        READ(20,ERR=110,END=120) rdrnam(1:4)
        READ(20,ERR=110,END=120) rdrlat,rdrlon,radarx,radary,rdralt
        READ(20,ERR=110,END=120) dazim,rngmin,rngmax

        ALLOCATE(hgtobs(nx, ny, ntilt),stat=istatus)
        ALLOCATE(rngoelv(nx, ny, ntilt),stat=istatus)
        ALLOCATE(elvobs(ntilt),stat=istatus)
        ALLOCATE(vroelv(nx, ny, ntilt),stat=istatus)
        ALLOCATE(refoelv(nx, ny, ntilt),stat=istatus)
        vroelv = -999.0
        refoelv = -999.0
        ALLOCATE(vrotmp(nx, ny, ntilt),stat=istatus)
        ALLOCATE(refotmp(nx, ny, ntilt),stat=istatus)

        READ(20,ERR=110,END=120) elvobs
        READ(20,ERR=110,END=120) hgtobs
        READ(20,ERR=110,END=120) rngoelv
        READ(20,ERR=110,END=120) vroelv
        READ(20,ERR=110,END=120) refoelv

        CLOSE(UNIT=20)

      ENDIF
    ELSE
      IF(refpltopt > 0 .or. vrpltopt > 0)THEN
        IF(datanamopt == 1)then
          filename=trim(rdrdtadir)//'/'//trim(radarnam)//trim(timnam)
          write(*,*)trim(filename)
        ELSE
          filename=trim(rdrdtadir)//'/'//trim(rdrdtanam(l))
          write(*,*)trim(filename)
        ENDIF
        OPEN(UNIT=20,FILE=trim(filename),STATUS='old', &
             FORM='unformatted',IOSTAT= istatus )
        IF( istatus /= 0) GO TO 999
        READ(20,ERR=110,END=120)nx,ny,nz,ntilt
        ALLOCATE(elvobs(ntilt),stat=istatus)
        ALLOCATE(xs(nx),stat=istatus)
        ALLOCATE(ys(ny),stat=istatus)
        ALLOCATE(zps(nx,ny,nz),stat=istatus)
        ALLOCATE(azmsc(nx,ny),stat=istatus)
        ALLOCATE(elvsc(nx,ny,nz),stat=istatus)
        ALLOCATE(rngsc(nx,ny,nz),stat=istatus)
        ALLOCATE(hgtobs(nx, ny, ntilt),stat=istatus)
        ALLOCATE(rngoelv(nx, ny, ntilt),stat=istatus)
        ALLOCATE(vroelv(nx, ny, ntilt),stat=istatus)
        ALLOCATE(refoelv(nx, ny, ntilt),stat=istatus)
        READ(20,ERR=110,END=120)elvobs
        READ(20,ERR=110,END=120)xs
        READ(20,ERR=110,END=120)ys
        READ(20,ERR=110,END=120)zps
        READ(20,ERR=110,END=120)azmsc
        READ(20,ERR=110,END=120)elvsc
        READ(20,ERR=110,END=120)rngsc
        READ(20,ERR=110,END=120)hgtobs
        READ(20,ERR=110,END=120)rngoelv
        READ(20,ERR=110,END=120)vroelv
        READ(20,ERR=110,END=120)refoelv
        DEALLOCATE(xs,stat=istatus)
        DEALLOCATE(ys,stat=istatus)
        DEALLOCATE(zps,stat=istatus)
        DEALLOCATE(azmsc,stat=istatus)
        DEALLOCATE(elvsc,stat=istatus)
        DEALLOCATE(rngsc,stat=istatus)
        ALLOCATE(vrotmp(nx, ny, ntilt),stat=istatus)
        ALLOCATE(refotmp(nx, ny, ntilt),stat=istatus)
        CLOSE(UNIT=20)
      ENDIF
    ENDIF

    IF(vapltopt > 0)THEN
      ALLOCATE(vphie(nx, ny, ntilt),stat=istatus)
      ALLOCATE(vatmp(nx, ny, ntilt),stat=istatus)
      IF(datanamopt == 1)THEN
        filename=trim(vadtadir)//'/'//trim(vanam)//trim(timnam)
      ELSE
        filename=trim(vadtadir)//'/'//trim(vadtanam(l))
      write(*,*)trim(filename)
      ENDIF

      OPEN(UNIT=30,FILE=trim(filename),STATUS='old', &
           FORM='unformatted',IOSTAT= istatus )
      IF( istatus /= 0) GO TO 999
      READ(30,ERR=110,END=120) vphie
      CLOSE(UNIT=30)

      DO itilt=1,ntilt
        DO j=1,ny
          DO i=1,nx
            vatmp(i,j,itilt)=vphie(i,j,itilt)
            if(vphie(i,j,itilt) < mv)then
              vphie(i,j,itilt)=-999.0
              vatmp(i,j,itilt)=0.0
            endif
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    IF(varpltopt > 0)THEN
      ALLOCATE(var3d(nx,ny,ntilt),stat=istatus)
      IF(datanamopt == 1)THEN
        filename=trim(vardtadir)//'/'//trim(varnam)//trim(timnam)
      write(*,*)trim(filename)
      ELSE
        filename=trim(vardtadir)//'/'//trim(vardtanam(l))
      write(*,*)trim(filename)
      ENDIF
      OPEN(UNIT=40,FILE=trim(filename),STATUS='old', &
           FORM='unformatted',IOSTAT= istatus )
      IF( istatus /= 0) GO TO 999
      READ(40,ERR=110,END=120)var3d
      CLOSE(UNIT=40)
    ENDIF

    IF(refpltopt > 0 .or. vrpltopt > 0)THEN
      DO itilt=1,ntilt
        DO j=1,ny
          DO i=1,nx
            vrotmp(i,j,itilt)=vroelv(i,j,itilt)
            refotmp(i,j,itilt)=refoelv(i,j,itilt)
            if(vroelv(i,j,itilt) < mv)then
              vroelv(i,j,itilt)=-999.0
              vrotmp(i,j,itilt)=0.0
            endif
            if(refoelv(i,j,itilt) < mv)then
              refoelv(i,j,itilt)=-999.0
              refotmp(i,j,itilt)=0.0
            endif
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    M=nx-3
    N=ny-3

    ALLOCATE(xc(nx,ny),stat=istatus)
    ALLOCATE(yc(nx,ny),stat=istatus)
    ALLOCATE(x(nx),stat=istatus)
    ALLOCATE(y(ny),stat=istatus)
    dx=dxr
    dy=dyr
    if(coordopt == 1)then
      swx=0.0-(nx-3)*dx*0.5  ! -dx/2.
      swy=0.0-(ny-3)*dx*0.5  ! -dx/2.
    else
      swx=0.0
      swy=0.0
    endif
    DO i=1,nx
      x(i)=(i-1)*dx-dx/2.+swx
    ENDDO
    DO j=1,ny
      y(j)=(j-1)*dy-dy/2.+swy
    ENDDO

    IF(xbgn == 0. .and. xend == 0.)THEN
      xl=(x(2)-dx/2.)*0.001
      xr=(x(nx-1)-dx/2.)*0.001
    ELSE
      xl=xbgn
      xr=xend
    ENDIF

    IF(ybgn == 0. .and. xend == 0.)THEN
      yb=(y(2)-dy/2.)*0.001
      yt=(y(ny-1)-dy/2.)*0.001
    ELSE
      yb=ybgn
      yt=yend
    ENDIF

!    write(*,*)'x', x(1), x(2), 'xl', xl, 'xr', xr

    dx=dx*0.001
    dy=dy*0.001
    xscale=(nx-3)*dx
    yscale=(ny-3)*dy

    y2=y1+((x2-x1)*yscale)/xscale

!#
!# Set coordinates at the grid volume center
!#
    xc=0.0
    yc=0.0

    DO j=1,ny
      DO i=1,nx-1
        xc(i,j) = x(i) *0.001
      ENDDO
    ENDDO

    DO j=1,ny-1
      DO i=1,nx
        yc(i,j) = y(j) *0.001
      ENDDO
    ENDDO

    iprio(1)=refprio
    iprio(2)=vrprio
    iprio(3)=vaprio
    iprio(4)=varprio

    iovr(1) = refovr
    iovr(2) = vrovr
    iovr(3) = vaovr
    iovr(4) = varovr

    iplt(1) = refpltopt
    iplt(2) = vrpltopt
    iplt(3) = vapltopt
    iplt(4) = varpltopt

    nplt=0
    DO i=1,nprio
      if(iplt(i)>0)then
        nplt=nplt+1
      endif
    ENDDO

    IF(datanamopt == 2)THEN
      WRITE(string2,'(i4,a,i2,a,i2,i4,a,i2,a,i2,a)') &
            iyr,'/',imon,'/',iday,ihour,':',imin,':',isec
    ELSE
      WRITE(string2,'(''T='',F8.1,A)')time, 's'
    ENDIF
!-----------------------------------------------------------------------
!
! Start to plot
!
!-----------------------------------------------------------------------
    write(*,*)'start to plot'
!#
!# setup device and initialize ZXPLOT
!#
    IF(l == 1) THEN
      CALL xpsfn(trim(psfile),5)
      CALL XDEVIC
      CALL XCFONT(2)
      CALL XPSPAC(x1,x2,y1,y2)
      CALL XMAP(xl,xr,yb,yt)
      CALL XLABMASK(1)

      CALL xlabmask(lbmaskopt)
      IF(ctrlbopt == 0) CALL xcltyp(0)
      IF(ctrlbopt == 1) CALL xclfmt('(f10.1)')
      IF(ctrlbopt == 2) CALL xclfmt('(I10)')

      CALL xcmixl
      IF( ctrstyle == 2) THEN
        CALL xcfull
      ELSE IF( ctrstyle == 3) THEN
        CALL xcdash
      END IF

      CALL xclfrq( ctrlbfrq )

      CALL xctrbadv( 1 )  ! Turn on missing value checking for contouring
!
!  Default value of -999.0 for the flag is used.
!
      CALL xbadval ( -999.0 ) ! Set the missing value flag to -9999.0
    END IF

    DO itilt=tiltbgn,tiltend
!#
!# Plot axes and labels
!#
      IF(l == 1) first_fram = 0
      ovrid = 0
      DO ip = 1, nplt

        IF(refprio == ip)THEN
          IF(refpltopt > 0)THEN
            first_fram = 1
            DO k=1,nprio
              if(k /= 1 .and. iprio(k) == ip+1 .and. iplt(k) > 0 .and. &
                 iovr(k) == 1)then
                first_fram = 0
                GOTO 200
              endif
            ENDDO
        200 CONTINUE
            WRITE(string,'(a,f5.2)') 'WS88d Z at elevation ',elvobs(itilt)
            if(refovr > 0)then
              ovrid=ovrid+1
            else
              ovrid=0
            endif
            CALL ctr2d(refoelv(:,:,itilt),refpltopt,coltab_ref,coltabfn_ref, &
                       ibgncol_ref,iendcol_ref,refminc,refmaxc,refinc,       &
                       refhlf,refzro,refovr,'Z(dBZ)',refotmp(:,:,itilt))

            IF(refovr == 0)THEN
              CALL xcharc(xscale*0.5, yscale+7*dy, trim(string2) )
            ENDIF
            CALL xcharc(xscale*0.5, yscale+4*dy, trim(string) )
          ENDIF ! IF(refpltopt > 0)THEN
        ENDIF ! IF(refprio == ip)THEN

        IF(vrprio == ip)THEN
          IF(vrpltopt > 0)THEN
            first_fram = 1
            DO k=1,nprio
              if(k /= 2 .and. iprio(k) == ip+1 .and. iplt(k) > 0 .and. &
                 iovr(k) == 1)then
                first_fram = 0
                GOTO 300
              endif
            ENDDO
        300 CONTINUE
            if(vrovr > 0)then
              ovrid=ovrid+1
            else
              ovrid=0
            endif
            WRITE(string,'(a,f5.2)') 'WS88d Vr at elevation ',elvobs(itilt)
            CALL ctr2d(vroelv(:,:,itilt),vrpltopt,coltab_vr,coltabfn_vr,        &
                       ibgncol_vr,iendcol_vr,vrminc,vrmaxc,vrinc,               &
                       vrhlf,vrzro,vrovr,'Vr(m/s)',vrotmp(:,:,itilt))

            IF(vrovr == 0)THEN
              CALL xcharc(xscale*0.5, yscale+7*dy, trim(string2) )
            ENDIF
            CALL xcharc(xscale*0.5, yscale+4*dy, trim(string) )

          ENDIF ! IF(vrpltopt > 0)THEN
        ENDIF ! IF(vrprio == ip)THEN

        IF(vaprio == ip)THEN
          IF(vapltopt > 0)THEN
            first_fram = 1
            DO k=1,nprio
              if(k /= 3 .and. iprio(k) == ip+1 .and. iplt(k) > 0 .and. &
                 iovr(k) == 1)then
                first_fram = 0
                GOTO 400
              endif
            ENDDO
        400 CONTINUE
            WRITE(string,'(a,f5.2)') 'WS88d cross-beam component at elevation ',elvobs(itilt)
            if(vaovr > 0)then
              ovrid=ovrid+1
            else
              ovrid=0
            endif
            CALL ctr2d(vphie(:,:,itilt),vapltopt,coltab_va,coltabfn_va,         &
                       ibgncol_va,iendcol_va,vaminc,vamaxc,vainc,               &
                       vahlf,vazro,vaovr,'Va(m/s)',vatmp(:,:,itilt))

            IF(vaovr == 0)THEN
              CALL xcharc(xscale*0.5, yscale+7*dy, trim(string2) )
            ENDIF
            CALL xcharc(xscale*0.5, yscale+4*dy, trim(string) )

          ENDIF
        ENDIF

        IF(varprio == ip)THEN
          IF(varpltopt > 0)THEN
            first_fram = 1
            DO k=1,nprio
              if(k /= 4 .and. iprio(k) == ip+1 .and. iplt(k) > 0 .and. &
                 iovr(k) == 1)then
                first_fram = 0
                GOTO 500
              endif
            ENDDO
        500 CONTINUE
            WRITE(string,'(a,f5.2)') 'cor(Vr, Va) at elevation ',elvobs(itilt)
            if(varovr > 0)then
              ovrid=ovrid+1
            else
              ovrid=0
            endif
            CALL ctr2d(var3d(:,:,itilt),varpltopt,coltab_var,coltabfn_var,      &
                       ibgncol_var,iendcol_var,varminc,varmaxc,varinc,          &
                       varhlf,varzro,varovr,'cor(Vr,Va)',var3d(:,:,itilt))

            IF(varovr == 0)THEN
              CALL xcharc(xscale*0.5, yscale+7*dy, trim(string2) )
            ENDIF
            CALL xcharc(xscale*0.5, yscale+2*dy, trim(string) )

          ENDIF
        ENDIF

        IF(first_fram == 1)THEN

          CALL XSTPJGRD(mapproj,trulat1,trulat2,trulon,ctrlat,ctrlon,       &
                        xscale,yscale,xorig,yorig)
          CALL xwindw(xl,xr,yb,yt)

          IF( ovrmap /= 0 ) THEN
            DO i=1,MIN(mxmapfile,nmapfile)
              lmapfile=len_trim(mapfile(i))
              WRITE(6,'(1x,a,a)') 'Input was ',trim(mapfile(i))

              INQUIRE(FILE=trim(mapfile(i)), EXIST = fexist )
              IF( .NOT.fexist) THEN
                WRITE(6,'(a)') 'Map file '//mapfile(i)(1:lmapfile)//         &
                     ' not found. Corresponding map no plotted.'
              ELSE
                CALL xdrawmap(10,mapfile(i)(1:lmapfile),latgrid,longrid)
              END IF
            END DO
          END IF

          CALL xwdwof
          CALL XFRAME
        ENDIF

      ENDDO ! DO ip = 1, nplt
    ENDDO  ! DO itilt=tiltbgn, tiltend

    IF(ALLOCATED(hgtobs))  DEALLOCATE(hgtobs,stat=istatus)
    IF(ALLOCATED(rngoelv)) DEALLOCATE(rngoelv,stat=istatus)
    IF(ALLOCATED(vroelv))  DEALLOCATE(vroelv,stat=istatus)
    IF(ALLOCATED(refoelv)) DEALLOCATE(refoelv,stat=istatus)
    IF(ALLOCATED(vrotmp))  DEALLOCATE(vrotmp,stat=istatus)
    IF(ALLOCATED(refotmp)) DEALLOCATE(refotmp,stat=istatus)
    IF(ALLOCATED(elvobs))  DEALLOCATE(elvobs,stat=istatus)

    DEALLOCATE(xc,stat=istatus)
    DEALLOCATE(yc,stat=istatus)
    DEALLOCATE(x,stat=istatus)
    DEALLOCATE(y,stat=istatus)

    IF(varpltopt > 0)THEN
      DEALLOCATE(var3d,stat=istatus)
    ENDIF

  ENDDO  ! DO l=1,(nn-n0+1)

  CALL XGREND  ! End graphics
  STOP

  999   CONTINUE
  WRITE(6,'(1x,a,a,a,/1x,i3)')                        &
      'Error occured when opening file ',filename
  STOP

  110   CONTINUE
  WRITE(6,'(/a/)') ' Error reading data in BN2READ'
  STOP

  120   CONTINUE
  WRITE(6,'(/a/)') ' End of file reached in BN2READ'
  STOP

  100  CONTINUE
  WRITE(6,'(/a/)') ' Error reading namelist file'
  STOP

END PROGRAM ppiplt

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CTR2D                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE ctr2d(a,apltopt,coltab_a,coltabfn_a,ibgncol_a,iendcol_a,     &
                 aminc,amaxc,ainc,ahlf,azro,aovr,obsnam,atmp)


!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    Generate contour plots of 2-d field A given its coordinates
!      using ZXPLOT package..
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Mingjing
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  USE VAR_DEF
  IMPLICIT NONE

  REAL :: a(nx,ny),atmp(nx,ny)
  INTEGER :: coltab_a,ibgncol_a,iendcol_a
  CHARACTER (LEN=256) :: coltabfn_a
  CHARACTER (LEN=*) :: obsnam
  INTEGER :: apltopt
  REAL :: aminc,amaxc,ainc
  INTEGER :: ahlf,azro,aovr

  INTEGER :: NCL
  INTEGER :: IW(M,N)
  REAL :: CL(100)
  REAL :: xw(M*N), yw(M*N)
  REAL :: ZMIN, ZMAX, ZINC
  INTEGER :: LEN0
  INTEGER :: MODE
  CHARACTER (LEN=180) :: ch

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF( coltab_a == -1 ) CALL xstctfn(coltabfn_a)
  CALL xsetclrs(coltab_a)
  CALL xcolor(1)
  CALL xcfont(3)
  CALL xchsiz(0.03* yscale)  ! Set character size
  IF(aovr == 0)THEN
    CALL XAXSCA(xl,xr,dx,yb,yt,dy)
  ENDIF

  IF(apltopt == 1 .or. apltopt == 3)THEN
    CL(1)=0.0
    CL(2)=ainc       ! Set tentative contour interval
    CALL xnctrs(2,50)  ! Set lower and upper limits
    MODE=1
    CALL xctrclr(ibgncol_a,iendcol_a)
    CALL xctrlim(aminc, amaxc)
    CALL xwindw(xl,xr,yb,yt)
    CALL XColfil(a(2:nx-2,2:ny-2),xc(2:nx-2,2:ny-2),  &
                 yc(2:nx-2,2:ny-2),IW,xw,yw,M,M,N, CL, NCL,MODE)
    CALL xcpalet(2)            ! Plot vertical   color palette
    CALL xwdwof
    ZMAX=MAXVAL(a(2:nx-2,2:ny-2))
    ZMIN=MINVAL(atmp(2:nx-2,2:ny-2))
    ZINC=CL(2)-CL(1)
    WRITE(ch,'(a,a,a,f6.2,a,f6.2)')obsnam, '      ', 'MIN=', ZMIN, '  MAX=', ZMAX
    LEN0=LEN_TRIM(ch)
    IF(aovr == 1)THEN
      CALL xcharr(xr, yb-(yt-yb)*0.14-yb-(yt-yb)*0.04*ovrid ,ch(1:LEN0))
    ELSE
      CALL xcharr(xr, yb-(yt-yb)*0.14,ch(1:LEN0))
      CALL xcharc(xscale*0.5, yb-(yt-yb)*0.09, 'x (km)')
      CALL XCHORI( 90.0 )
      CALL xcharc(xl-(xr-xl)*0.12, yscale*0.5, 'y (km)')
      CALL XCHORI( 0.0 )
    ENDIF
  ENDIF

  IF(apltopt >= 2)THEN
    CL(1)=0.0
    CL(2)=ainc     ! Set tentative contour interval
    call xnctrs(2,50)  ! Set lower and upper limits
    MODE=1
    CALL xhlfrq(ahlf)
    CALL xczero (azro)
    CALL xctrclr(1,1)
    CALL xctrlim(aminc, amaxc)
    CALL XCLTYP(0)
    CALL XCONTA(a(2:nx-2,2:ny-2),xc(2:nx-2,2:ny-2),  &
                yc(2:nx-2,2:ny-2),IW,M,M,N, CL, NCL,MODE)
    ZMAX=MAXVAL(a(2:nx-2,2:ny-2))
    ZMIN=MINVAL(atmp(2:nx-2,2:ny-2))
    ZINC=CL(2)-CL(1)
    WRITE(ch,'(a,a,a,f6.2,a,f6.2,a,f6.2)')obsnam, '      ',  &
          'MIN=', ZMIN, '  MAX=', ZMAX, '  INC', ZINC
    LEN0=LEN_TRIM(ch)
    IF(aovr == 1)THEN
      CALL xcharr(xr, yb-(yt-yb)*0.14-(yt-yb)*0.04*ovrid ,ch(1:LEN0))
    ELSE
      CALL xcharr(xr, yb-(yt-yb)*0.14,ch(1:LEN0))
      CALL xcharc(xscale*0.5, yb-(yt-yb)*0.09, 'x (km)')
      CALL XCHORI( 90.0 )
      CALL xcharc(xl-(xr-xl)*0.12, yscale*0.5, 'y (km)')
      CALL XCHORI( 0.0 )
    ENDIF
  ENDIF

  RETURN
END SUBROUTINE ctr2d
