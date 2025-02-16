PROGRAM  pltgrid
!#######################################################################
!
!  Program to plot model grids
!
!  Commands:
!
!    zxncarf77 -o pltgrid pltgrid.f
!    pltgrid < pltgrid.input
!
!  or in ARPS root directory, enter
!
!    makearps pltgrid
!    bin/pltgrid  < input/pltgrid.input
!
! MODIFICATION HISTORY:
!
! 11/27/2001 (K. Brewster)
!
! Added plotting of station locations to plot, correcting some
! previous test code.
!
! 03/03/2003 (K. Brewster)
!
! Modified scaling of station markers and name sizes.
!
! 09/04/2007 (K. Brewster)
!
! Added plotting of radars and radar range rings.
!
! 08/06/2008 (Fanyou Kong)
! Add ctrlat1 and ctrlon1 for specifying center locations (along with xy_or_ll)
!
! 03/18/2013 (Guoqing Ge)
! Add plot_polys for plotting polylines (hurricane tracks, tornado paths, etc)
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: nx,ny
  INTEGER :: mapproj
  REAL :: dx,dy,ctrlat,ctrlon,trulat1,trulat2,trulon,sclfct,xl,yl

  INTEGER :: nnwgrd, nnwgrdmax
  PARAMETER (nnwgrdmax = 10) ! maximum number of new grids

  INTEGER :: nx1(nnwgrdmax),ny1(nnwgrdmax)
  INTEGER :: xy_or_ll
  REAL :: xctr1(nnwgrdmax) , yctr1(nnwgrdmax)
  REAL :: ctrlat1(nnwgrdmax) , ctrlon1(nnwgrdmax)
  REAL :: dx1(nnwgrdmax),dy1(nnwgrdmax)

  NAMELIST /grid/ nx,ny,dx,dy,ctrlat,ctrlon

  NAMELIST /projection/ mapproj,trulat1,trulat2,trulon,sclfct

  NAMELIST /newgrid/ nnwgrd,nx1,ny1,xy_or_ll,xctr1,yctr1,ctrlat1,ctrlon1,dx1,dy1

  INTEGER :: ovrmap,mapgrid,mapgridcol,nmapfile
  INTEGER, PARAMETER :: mxmapfile= 10
  CHARACTER (LEN=256) :: mapfile(mxmapfile)
  INTEGER :: mapcol(mxmapfile)
  INTEGER :: mapline_style(mxmapfile)

  REAL :: latgrid,longrid, xorig, yorig
  COMMON /mappar/ ovrmap
  COMMON /mappar1/ nmapfile,mapcol,mapline_style,mapfile
  COMMON /mappar2/ mapgrid,mapgridcol,latgrid,longrid
  NAMELIST /map_plot/ ovrmap,mapgrid,mapgridcol,latgrid,longrid, &
            nmapfile,mapcol,mapline_style,mapfile

  INTEGER :: pltstn
  CHARACTER (LEN=256) :: stnfile
  CHARACTER (LEN=2)   :: statelist(60)
  NAMELIST /station_plot/ pltstn, stnfile, statelist

  CHARACTER (LEN=20) :: stnname
  CHARACTER (LEN=2) :: state
  CHARACTER (LEN=4) :: stnsymbol
  CHARACTER (LEN=80) :: runname
  CHARACTER (LEN=23) :: timestr
  REAL :: lat,long, x1,x2,y1,y2
  INTEGER :: elv,TYPE,length
  REAL :: xstn,ystn

  INTEGER :: mxstalo   ! maximum number of observation stations
  PARAMETER (mxstalo=500)

  INTEGER :: ovrstaopt
  INTEGER :: ovrstam,staset,ovrstan,ovrstav,stacol,markprio,wrtstax
  INTEGER :: nsta_typ,sta_typ(10),sta_marktyp(10),sta_markcol(10)
  REAL :: sta_marksz(10)
  REAL :: wrtstad

  INTEGER :: nsta,nstapro(mxstalo),nstatyp(mxstalo)
  REAL :: latsta(mxstalo), lonsta(mxstalo)
  CHARACTER (LEN=2) :: s_state(mxstalo)
  CHARACTER (LEN=5) :: s_name(mxstalo)
  CHARACTER (LEN=20) :: s_site(mxstalo)
  INTEGER :: s_elev(mxstalo)

  CHARACTER (LEN=256) :: stalofl
  INTEGER :: lstalofl

  NAMELIST /plot_sta/ ovrstaopt,ovrstam,ovrstan,ovrstav,wrtstax,        &
      wrtstad, stacol, markprio, nsta_typ, sta_typ, sta_marktyp,        &
      sta_markcol,sta_marksz,stalofl

!
  INTEGER, PARAMETER :: maxrtab=300
  CHARACTER (LEN=12) :: radtab(maxrtab)
  REAL :: rtablat(maxrtab)
  REAL :: rtablon(maxrtab)

  INTEGER, PARAMETER :: maxrad=100
  REAL :: radlat(maxrad)
  REAL :: radlon(maxrad)

  INTEGER :: pltradopt,nrad
  REAL    :: rngring(maxrad)
  CHARACTER (LEN=12) :: radname(maxrad)

  NAMELIST /plot_rad/ pltradopt,nrad,rngring,radname

  INTEGER, PARAMETER :: max_polys=20
  INTEGER :: nPolys, poly_color(max_polys), iFormat
  INTEGER :: poly_pltline(max_polys), poly_pltmark(max_polys)
  INTEGER :: poly_marktyp(max_polys)
  REAL    :: poly_marksz(max_polys) 
  CHARACTER (LEN=256) :: poly_fname(max_polys)
  NAMELIST /plot_polys/ nPolys, poly_color,poly_fname, poly_pltline,   &
                            poly_pltmark, poly_marksz, poly_marktyp
!
  CHARACTER (LEN=1) :: dummy
  REAL :: swx,swy,ctrx,ctry
  LOGICAL :: fexist
  REAL :: truelat(2)
  REAL :: xl1, yl1, chsize,yoffset
  REAL :: angle,rlat,rlon,xpos,ypos
  REAL :: cornerlat,cornerlon
  REAL :: tlat,tlon
  INTEGER :: ilat,ilatmin,ilatsec,ilon,ilonmin,ilonsec
  INTEGER :: i,k,iang,istat,itype,nradtab,lmapfile,imkrfil
  INTEGER :: iscore,iostatus
  LOGICAL :: found

  CHARACTER(LEN=256) :: outfilename
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  imkrfil=1
!
!-----------------------------------------------------------------------
!
!  Input control parameters map plotting
!
!-----------------------------------------------------------------------
!
  WRITE(6,'(/a/)')                                                      &
      ' Please input control parameters for map plotting.'

  READ(5,grid,ERR=100)

  WRITE(6,'(/5x,a,i8)')'Input nx was',nx
  WRITE(6,'(/5x,a,i8)')'Input ny was',ny
  WRITE(6,'(/5x,a,f12.6,a)')'Input dx was',dx,' meters'
  WRITE(6,'(/5x,a,f12.6,a)')'Input dy was',dy,' meters'


  WRITE(6,'(/5x,a,f12.6,a)')                                            &
       'Input ctrlat was ',ctrlat,' degree North'

  WRITE(6,'(/5x,a,f12.6,a)')                                            &
       'Input ctrlon was ',ctrlon,' degree East'

!-----------------------------------------------------------------------
!
!  Input map projection parameters
!
!-----------------------------------------------------------------------
!
  READ (5,projection,ERR=100)

  WRITE(6,'(/5x,a,i4)') 'Input mapproj was ',mapproj

  WRITE(6,'(/5x,a,f12.6,a)')                                            &
       'Input trulat1 was ',trulat1,' degree North'

  WRITE(6,'(/5x,a,f12.6,a)')                                            &
       'Input trulat2 was ',trulat2,' degree North'

  WRITE(6,'(/5x,a,f12.6)')                                              &
      'The latitude of the center of the model domain was ',ctrlat

  WRITE(6,'(/5x,a,f12.6)')                                              &
      'The longitude of the center of the model domain was ',ctrlon

  WRITE(6,'(/5x,a,f12.6,a)')                                            &
       'Input trulon was ',trulon,' degree East'

  IF ( mapproj == 0 ) THEN
    trulat1 = ctrlat
    trulat2 = ctrlat
    trulon  = ctrlon
  END IF
!
!-----------------------------------------------------------------------
!
!  Set the output grid and the variable control parameters
!
!-----------------------------------------------------------------------
!
  READ (5,newgrid)

  PRINT*,' Input nnwgrd=',nnwgrd

  IF( nnwgrd > nnwgrdmax)  THEN
    PRINT*,'Number of new grids more than maximum allowed.'
    PRINT*,'Increase the size of nnwgrdmax in the program.'
    STOP
  END IF

  PRINT*,' Input nx1=',(nx1(i),i=1,nnwgrd)
  PRINT*,' Input ny1=',(ny1(i),i=1,nnwgrd)
  PRINT*,' Input dx1=',(dx1(i),i=1,nnwgrd)
  PRINT*,' Input dy1=',(dy1(i),i=1,nnwgrd)
  if(xy_or_ll == 1) then
  PRINT*,' Input xctr1=',(xctr1(i),i=1,nnwgrd)
  PRINT*,' Input yctr1=',(yctr1(i),i=1,nnwgrd)
  else
  PRINT*,' Input ctrlat1=',(ctrlat1(i),i=1,nnwgrd)
  PRINT*,' Input ctrlon1=',(ctrlon1(i),i=1,nnwgrd)
  endif

  READ(5,map_plot,ERR=100)

! READ(5,station_plot,ERR=100)

  GO TO 10
  100   WRITE(6,'(a)') 'Error reading NAMELIST file. Program stopped.'
  STOP

  10    CONTINUE

!
!-----------------------------------------------------------------------
!
!  Initialize ZXPLOT plotting package
!
!-----------------------------------------------------------------------
!
  CALL xdevic_new(1,outfilename,0,0)       ! use default
  CALL xafstyl(1)
  CALL xsetclrs_new(1,0)
  CALL xcolor(1)
  CALL xartyp(2)

  CALL xdspac(0.9)

  CALL xaxfmt( '(i10)' )

!-----------------------------------------------------------------------
!
!  Set up map projection.
!
!-----------------------------------------------------------------------
!
  xl = (nx-3)*dx * 0.001
  yl = (ny-3)*dy * 0.001
  xorig = 0.0
  yorig = 0.0

  CALL xstpjgrd(mapproj,trulat1,trulat2,trulon,                         &
                ctrlat,ctrlon,xl,yl,xorig,yorig)

  CALL xxytoll(1,1,xorig*1000.0,yorig*1000.0,cornerlat,cornerlon)

  WRITE(6,'(/1x,a,2f12.6)')                                             &
      'Lat/lon at the SW corner of base grid=',cornerlat,cornerlon

  CALL xxytoll(1,1,(xorig+xl)*1000.0,(yorig+yl)*1000.0,                 &
               cornerlat,cornerlon)
  WRITE(6,'(1x,a,2f12.6/)')                                             &
      'Lat/lon at the NE corner of base grid=',cornerlat,cornerlon

  IF( xl >= yl) THEN
    CALL xpspac(0.1,0.9, 0.5-yl/xl*0.4,0.5+yl/xl*0.4)
  ELSE
    CALL xpspac(0.5-xl/yl*0.4,0.5+xl/yl*0.4,0.1, 0.9)
  END IF
!
  CALL xmap(0.0,xl, 0.0, yl)
  CALL xbox(0.0,xl, 0.0, yl)
!  CALL xaxsca(0.0,xl, dx/1000.0, 0.0, yl, dy/1000.0 )

!  CALL xcolor(9)

  IF( ovrmap /= 0 ) THEN

    DO i=1,MIN(mxmapfile,nmapfile)

      lmapfile=LEN_TRIM(mapfile(i))
      WRITE(6,'(1x,a,a)') 'Input was ',mapfile(i)(1:lmapfile)

      INQUIRE(FILE=mapfile(i)(1:lmapfile), EXIST = fexist )
      IF( .NOT.fexist) THEN
        WRITE(6,'(a)') 'Map file '//mapfile(i)(1:lmapfile)//            &
            ' not found. Corresponding map not plotted.'
      ELSE
        CALL xcolor(mapcol(i))
        IF(mapline_style(i) == 1) THEN
          CALL xthick(1)
          CALL xbrokn(6,3,6,3)
        ELSE IF(mapline_style(i) == 2) THEN
          CALL xthick(1)
        ELSE IF(mapline_style(i) == 3) THEN
          CALL xthick(3)
          CALL xfull
        END IF
        CALL xdrawmap_new(10,mapfile(i)(1:lmapfile),latgrid,longrid, i, mapgridcol)
      END IF

    END DO

  END IF

  CALL xthick(3)
!  CALL xthick(12)
  DO i=1,nnwgrd

    CALL xcolor(18+i)
    xl1 = (nx1(i)-3)*dx1(i) * 0.001
    yl1 = (ny1(i)-3)*dy1(i) * 0.001

  IF (xy_or_ll == 1) THEN
    CALL xxytoll(1,1,xctr1(i),yctr1(i),ctrlat1(i),ctrlon1(i))
  ELSE IF (xy_or_ll == 2) THEN
    CALL xlltoxy(1,1,ctrlat1(i),ctrlon1(i),xctr1(i),yctr1(i))
  ELSE
    WRITE (*,*) "ERROR: xy_or_ll =",xy_or_ll," not a valid value."
    STOP
  END IF
    WRITE(6,'(/1x,a,2f15.4,2f12.6)')                                    &
           'xctr1,yctr1,ctrlat1,ctrlon1=',                              &
            xctr1(i),yctr1(i),ctrlat1(i),ctrlon1(i)

    CALL xxytoll(1,1,xctr1(i)-xl1*500.0,yctr1(i)-yl1*500.0,             &
         cornerlat,cornerlon)
    WRITE(6,'(1x,a,i3,a,2f12.6)')                                       &
           'Lat/lon at the SW corner of new grid no.',i,                &
           '=', cornerlat,cornerlon

    CALL xxytoll(1,1,xctr1(i)+xl1*500.0,yctr1(i)+yl1*500.0,             &
         cornerlat,cornerlon)
    WRITE(6,'(1x,a,i3,a,2f12.6/)')                                      &
          'Lat/lon at the NE corner of new grid no.',i,                 &
          '=', cornerlat,cornerlon

    CALL xbox(xctr1(i)/1000.0-xl1/2, xctr1(i)/1000.0+xl1/2,             &
              yctr1(i)/1000.0-yl1/2, yctr1(i)/1000.0+yl1/2)

  END DO

  ovrstam=0
  ovrstan=0
  ovrstav=0
  wrtstax=0

  print *, ' reading plot_sta'

  READ(5,plot_sta, ERR=72)

  WRITE(6,'(a)')                                                        &
      'Namelist plot_sta was successfully read.'

  lstalofl=LEN_TRIM(stalofl)
  WRITE(6,'(1x,a,a)') 'Station file name was ',stalofl(1:lstalofl)

  72 CONTINUE

  nsta=0

  IF(ovrstaopt /= 0) THEN
    INQUIRE(FILE=stalofl(1:lstalofl), EXIST = fexist )
    IF( .NOT.fexist) THEN
      WRITE(6,'(a)') 'WARNING: Surface station list file '             &
          //stalofl(1:lstalofl)//                                      &
          ' not found. Program continue in PLTGRID.'
    ELSE
      CALL read_station(stalofl,mxstalo,latsta,lonsta,nstatyp,          &
                  nstapro,nsta,s_name,s_state,s_site,s_elev)
      IF(nsta > 0) staset=1
    END IF

    print *, ' Done reading, nsta= ',nsta
    CALL xcolor(1)

    DO i=1,nsta
      found=.false.
      DO itype=1,nsta_typ
        IF(nstatyp(i) == sta_typ(itype)) THEN
          found=.true.
          EXIT
        END IF
      END DO
      IF(.not. found) CYCLE
      print *, ' nstatyp, itype = ',nstatyp(i),itype
      CALL xlltoxy(1,1,latsta(i),lonsta(i), xstn,ystn)
      print *, ' Plotting ',s_name(i),latsta(i),lonsta(i),.001*xstn,.001*ystn
      xstn=xstn*0.001
      ystn=ystn*0.001
      call xcolor(sta_markcol(itype))
      chsize=sta_marksz(itype)*yl*7.
      yoffset=sta_marksz(itype)*yl*10.
      IF( ovrstam > 0) THEN
        print *, 'plotting marker'
        CALL xmrksz(sta_marksz(itype))
        CALL xmarker(xstn,ystn,sta_marktyp(itype))
      END IF
      IF( ovrstan > 0) THEN
        CALL xchsiz(chsize)
        length=LEN_TRIM(s_name(i))
        IF( ovrstam > 0) THEN
          CALL xcharc(xstn, (ystn-yoffset), s_name(i)(1:length))
        ELSE
          CALL xcharc(xstn, ystn, s_name(i)(1:length))
        END IF
      END IF
    END DO
  END IF

  READ(5,plot_rad)
  WRITE(6,'(a)') 'Namelist plot_rad was successfully read.'
  READ(5,plot_polys)
  WRITE(6,'(a)') 'Namelist plot_polys was successfully read.'

  IF(pltradopt > 0 ) THEN

    open(31,file='radarinfo.dat',status='old',iostat=istat)
    IF(istat /= 0) THEN
      WRITE(6,'(a)') ' Did not find radarinfo.dat, trying data/adas'
      open(31,file='./data/adas/radarinfo.dat',status='old',iostat=istat)
    END IF
    IF(istat /= 0) THEN
      WRITE(6,'(a)') ' Did not find radarinfo.dat, trying ../data/adas'
      open(31,file='../data/adas/radarinfo.dat',status='old')
    END IF
!
! Read station table data
!
    read(31,'(a1)') dummy
    DO i=1,maxrtab
      read(31,'(a4,20x,i4,2i8,i10,i7,i8)',iostat=istat) radtab(i), &
           ilat,ilatmin,ilatsec,ilon,ilonmin,ilonsec
      IF(istat /= 0) EXIT
      rtablat(i)=float(ilat)+(float(ilatmin)/60.)+(float(ilatsec)/3600.)
      rtablon(i)=float(ilon)+(float(ilonmin)/60.)+(float(ilonsec)/3600.)
      rtablon(i)=-rtablon(i)
    END DO
    nradtab=i-1
    WRITE(6,'(a,i6,a)') ' Read',nradtab,' NEXRAD radars from table'
!
! Match radar name and table name
!
    DO i=1,nrad
      radlat(i)=0.0
      radlon(i)=0.0
      DO k=1,nradtab
        IF( radname(i) == radtab(k) ) THEN
          WRITE(6,'(a,a)') ' Found ',radname(i)
          radlat(i)=rtablat(k)
          radlon(i)=rtablon(k)
          EXIT
        END IF
      END DO
    END DO
!
! Plot location and max radius for each radar
!
    DO i=1,nrad
      CALL xlltoxy(1,1,radlat(i),radlon(i), xstn,ystn)
      WRITE(*,'(1x,a,i4,a3,a4,4F10.2,a,F10.2)') 'Plotting ',i,' - ',radname(i),  &
             radlat(i),radlon(i),.001*xstn,.001*ystn,', rngring= ',rngring(i)
      xstn=xstn*0.001
      ystn=ystn*0.001
      call xcolor(sta_markcol(1))
      chsize=sta_marksz(1)*yl*6.
      yoffset=sta_marksz(1)*yl*10.
      CALL xmrksz(sta_marksz(1))
      CALL xmarker(xstn,ystn,sta_marktyp(1))
      length=LEN_TRIM(radname(i))
      CALL xchsiz(chsize)
      CALL xcharc(xstn, (ystn-yoffset), radname(i)(1:length))

      CALL xfull
      CALL xthick(1)
      CALL gcircle(radlat(i),radlon(i),0.,rngring(i),rlat,rlon)
      CALL xlltoxy(1,1,rlat,rlon,xpos,ypos)
      CALL xpenup((0.001*xpos),(0.001*ypos))
      DO iang=5,360,5
        angle=float(iang)
        CALL gcircle(radlat(i),radlon(i),angle,rngring(i),rlat,rlon)
        CALL xlltoxy(1,1,rlat,rlon,xpos,ypos)
        CALL xpendn((0.001*xpos),(0.001*ypos))
      END DO
    END DO
  END IF

!!!--------------------------------------------------
!!!  Begin: Plot ploy-lines
!!!  e.g. to plot hurricane tracks or tornado paths
!!!
  DO i = 1, nPolys
    OPEN(99,FILE=trim(poly_fname(i)),status='old',iostat=iostatus)
    IF(iostatus == 0) THEN
      CALL xfull
      CALL xthick(1)
      CALL xcolor(poly_color(i))
      PRINT *, 'Read polyline file'//trim(poly_fname(i))
      READ(99,*) iFormat
      IF (iFormat == 0) THEN 
         PRINT *, 'The file is in format (lat,lon), the corresponding (x,y)(km) is:'
      ELSE
         PRINT *, 'The file is in format (x,y)(km), the corresponding (lat,lon) is:'
      END IF

      READ(99,*,IOSTAT=iostatus, ERR=888) rlat,rlon
      IF (iFormat==0) then 
         CALL xlltoxy(1,1,rlat,rlon,xpos,ypos)
         IF ( poly_pltline(i) > 0 )                           &
                    CALL xpenup((0.001*xpos),(0.001*ypos))
         IF ( poly_pltmark(i) > 0 ) THEN
           CALL xmrksz(poly_marksz(i))
           CALL xmarker((0.001*xpos),(0.001*ypos),poly_marktyp(i))
         END IF
         PRINT *, 0.001*xpos, 0.001*ypos
      ELSE !(rlat,rlon) is now (x,y) (km)
         IF ( poly_pltline(i) > 0 ) CALL xpenup(rlat,rlon)
         IF ( poly_pltmark(i) > 0 ) THEN
           CALL xmrksz(poly_marksz(i))
           CALL xmarker(rlat,rlon,poly_marktyp(i))
         END IF
         CALL xxytoll(1,1,rlat*1000, rlon*1000,xpos,ypos)
         PRINT *, xpos, ypos
      END IF
      DO WHILE ( iostatus == 0)
        READ(99,*,IOSTAT=iostatus, ERR=888) rlat,rlon
        IF (iFormat==0) then 
           CALL xlltoxy(1,1,rlat,rlon,xpos,ypos)
           IF ( poly_pltline(i) > 0 )                           &         
                      CALL xpendn((0.001*xpos),(0.001*ypos))
           IF ( poly_pltmark(i) > 0 ) THEN
             CALL xmrksz(poly_marksz(i))
             CALL xmarker((0.001*xpos),(0.001*ypos),poly_marktyp(i))
           END IF
           PRINT *, 0.001*xpos, 0.001*ypos
        ELSE !(rlat,rlon) is now (x,y) (km)
           IF ( poly_pltline(i) > 0 ) CALL xpendn(rlat,rlon)
           IF ( poly_pltmark(i) > 0 ) THEN
             CALL xmrksz(poly_marksz(i))
             CALL xmarker(rlat,rlon,poly_marktyp(i))
           END IF
           CALL xxytoll(1,1,rlat*1000, rlon*1000,xpos,ypos)
           PRINT *, xpos, ypos
        END IF
      ENDDO
    END IF
    888  CLOSE(99)
  END DO
!!! 
!!! END OF plotting poly-lines
!!!----------------------------------!!!

  OPEN(41,file='tornscore.txt',status='old',iostat=iostatus)
  IF(iostatus == 0) THEN
    print *, ' Opened tornado location file: tornscore.txt'
    read(41,'(a80)') runname
    read(41,'(a23)') timestr
    read(41,'(a1)') dummy
    chsize=sta_marksz(1)*yl*2.5
    CALL xchsiz(chsize)
    DO
      read(41,'(12x,2f12.4,i5)',iostat=iostatus) tlat,tlon,iscore
      IF(iostatus /= 0) EXIT
      print *, ' torn lat,lon,score: ',tlat,tlon,iscore
      CALL xlltoxy(1,1,tlat,tlon,xstn,ystn)
      xstn=xstn*0.001
      ystn=ystn*0.001
      IF( iscore == 1) THEN
        call xcolor(9)
        CALL xcharc(xstn,ystn,'X')
      ELSE IF( iscore == 0) THEN
        call xcolor(20)
        CALL xcharc(xstn,ystn,'0')
      ELSE IF( iscore == -1) THEN
        call xcolor(3)
        CALL xcharc(xstn,ystn,'.')
      END IF
    END DO

    call xcolor(1)
    chsize=sta_marksz(1)*yl*8.
    CALL xchsiz(chsize)
    yoffset=sta_marksz(1)*yl*12.
    CALL xcharc((xl*0.5),(yl-yoffset),timestr)
  END IF

  CALL xgrend

  STOP
END PROGRAM  pltgrid

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE READ_STATION               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE read_station(infile,mxstalo,latsta,lonsta,                   &
           nstatyp,nstapro,nsta,sname,state,sitena,nelev)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  This subroutine will read external staion information.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!    Min Zou (6/1/97)
!
!  Modification history:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=*) :: infile
  INTEGER :: nsta,nstapro(mxstalo),nstatyp(mxstalo)
  REAL    :: latsta(mxstalo), lonsta(mxstalo)
  CHARACTER (LEN=5) :: sname(mxstalo)
  CHARACTER (LEN=2) :: state(mxstalo)
  CHARACTER (LEN=20) :: sitena(mxstalo)
  INTEGER :: nelev(mxstalo)
  CHARACTER (LEN=256) :: line

  OPEN(1,IOSTAT=ios,FILE=infile,STATUS='old',FORM='formatted')
  IF(ios /= 0) THEN     ! error during read
    istatus = -1
    WRITE(6,650) infile
    650       FORMAT(' +++ ERROR opening: ',a70,' +++')
    WRITE(6,651) ios
    651       FORMAT('     IOS code = ',i5)
    RETURN
  END IF

  nsta = 0

! Read only lines that begin with A-Z, a-z, or 0-9 -- treat the rest as
! comments

  DO i=1,mxstalo
    READ(1,'(A)',END=999) line
    IF ( ( line(:1) >= 'A' .AND. line(:1) <= 'Z' ) .OR.                 &
           ( line(:1) >= 'a' .AND. line(:1) <= 'z' ) .OR.               &
           ( line(:1) >= '0' .AND. line(:1) <= '9' ) ) THEN
      nsta = nsta + 1
      READ(line,101,ERR=999)sname(nsta),state(nsta),sitena(nsta),       &
          latsta(nsta),lonsta(nsta),nelev(nsta),nstatyp(nsta),          &
          nstapro(nsta)
    END IF
    101   FORMAT(a5,2X,a2,1X,a20,4X,f8.3,1X,f8.3,1X,i5,2X,i2,i1)
    102   FORMAT(a5,2X,a2,1X,a20,4X,f8.3,1X,f8.3,1X,i5,2X,i2,1X,i1)
  END DO
  999   CONTINUE
  CLOSE(1)

  RETURN
END SUBROUTINE read_station
