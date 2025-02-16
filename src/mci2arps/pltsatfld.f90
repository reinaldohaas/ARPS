PROGRAM pltsatfld
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 PROGRAM PLTSATFLD                    ######
!######                                                      ######
!######                Copyright (c) 1996                    ######
!######    Center for Analysis and Prediction of Storms      ######
!######    University of Oklahoma.  All rights reserved.     ######
!######                                                      ######
!##################################################################
!##################################################################
!
!  PURPOSE:
!
!  Plots data written by remapsat.
!
!  AUTHOR:
!
!  Keith Brewster, CAPS, September, 1997
!
!  LINKING:
!
!  ncargf77 -o pltsatmap pltsatmap.f pltmap.f maproj3d.f timelib3d.f
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nfield
  PARAMETER (nfield=1)

  INTEGER :: nx, ny
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
!
!-----------------------------------------------------------------------
!
!  Read-in data
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=6) :: satname
  REAL :: latsat
  REAL :: lonsat
  INTEGER :: itime
  INTEGER :: isource
  CHARACTER (LEN=6) :: fldname(nfield)
  REAL, ALLOCATABLE :: satfld(:,:,:)
!
!-----------------------------------------------------------------------
!
!  Map plotting variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: maxpts
  PARAMETER (maxpts = 1000)
  REAL :: latmap(maxpts),lonmap(maxpts)
  REAL :: xmap(maxpts),ymap(maxpts)
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
  CHARACTER (LEN=256) :: mapfile
  PARAMETER (mapfile='uscounty.mapdata')
  
  REAL :: x,y
  REAL :: latnot(2)
  REAL :: ctrx,ctry,swx,swy,nex,ney
  INTEGER :: i,j,k,iplt

  INTEGER :: istatus

!-----------------------------------------------------------------------
!
! NAMELIST declaration
!
!-----------------------------------------------------------------------

  NAMELIST /grid_dims/ nx, ny
  NAMELIST /fn/ fname


!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!-----------------------------------------------------------------------
!
!  Get user input
!
!-----------------------------------------------------------------------
!
  READ(5,grid_dims,END=100)
  WRITE(6,'(/a,a)')' Namelist block grid_dims successfully read.'
  
  READ(5,fn,END=100)
  WRITE(6,'(/a,a)')' Namelist block fn successfully read.'

  ALLOCATE( satfld(nx,ny,nfield), STAT=istatus)
  satfld = 0
!
!-----------------------------------------------------------------------
!
!  Read data
!
!-----------------------------------------------------------------------
!
  CALL rdsatfld(nx,ny,nfield,                                           &
                fname,satname,latsat,lonsat,                            &
                itime,isource,fldname,satfld)
!
  CALL opngks
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
  swx=ctrx-((nx-3)/2)*dx
  swy=ctry-((ny-3)/2)*dy
  nex=swx+(nx-3)*dx
  ney=swy+(ny-3)*dy
!
!-----------------------------------------------------------------------
!
!  Field loop
!
!-----------------------------------------------------------------------
!
  DO k=1,nfield
    PRINT *, ' plotting field: ',fldname(k)
!
    IF( fldname(k) == 'cttemp' ) THEN
!
!-----------------------------------------------------------------------
!
!  Plot cloud-top temperature contours
!
!-----------------------------------------------------------------------
!
      CALL set( 0.,1.0,0.0,1.0,                                         &
             -.1,1.1,-.1,1.1,1)
      CALL wtstr(0.05,1.05,satname,16,0,0)
      CALL wtstr(0.6,1.05,timplt,16,0,0)
      CALL wtstr(0.9,1.05,fldname(k),16,0,0)
      CALL set( 0.05,.95,0.05,0.95,                                     &
                swx,nex,swy,ney,1)
      CALL conrec(satfld(1,1,k),                                        &
                nx,nx,ny,0.,0.,20.,-1,-1,-127)
      CALL pltmap(maxpts,mapfile,latmap,lonmap,xmap,ymap)
      CALL frame
!
!-----------------------------------------------------------------------
!
!  Plot cloud-top temperature values
!
!-----------------------------------------------------------------------
!
!      CALL SET( 0.,1.0,0.0,1.0,
!    +           -.1,1.1,-.1,1.1,1)
!      call wtstr(0.05,1.05,satname,16,0,0)
!      call wtstr(0.6,1.05,timplt,16,0,0)
!      call wtstr(0.9,1.05,'fldname(k)',16,0,0)
!      CALL SET( 0.,1.0,0.0,1.0,
!    +           swx,nex,swy,ney,1)
!      DO 150 j=1,ny
!      DO 150 i=1,nx
!        IF(satfld(i,j,k).gt.-200. .and.
!    +         satfld(i,j,k).lt.500.) THEN
!          write(chplt,820) nint(satfld(i,j,k))
      820         FORMAT(i3)
!          x=swx+float(i-1)*dx
!          y=swy+float(j-1)*dy
!          call wtstr(x,y,chplt,8,0,0)
!        END IF
! 150     CONTINUE
!      CALL pltmap(maxpts,mapfile,latmap,lonmap,xmap,ymap)
!      CALL FRAME
    END IF
!
    IF( fldname(k) == 'albedo' ) THEN
!
!-----------------------------------------------------------------------
!
!  Plot albedo contours
!
!-----------------------------------------------------------------------
!
      CALL set( 0.,1.0,0.0,1.0,                                         &
           -.1,1.1,-.1,1.1,1)
      CALL wtstr(0.1,1.05,satname,16,0,0)
      CALL wtstr(0.4,1.05,timplt,16,0,0)
      CALL wtstr(0.8,1.05,fldname(k),16,0,0)
      CALL set( 0.05,0.95,0.05,0.95,                                    &
               swx,nex,swy,ney,1)
      CALL conrec(satfld(1,1,k),                                        &
                  nx,nx,ny,0.0,0.5,0.05,-1,-1,-127)
      CALL pltmap(maxpts,mapfile,latmap,lonmap,xmap,ymap)
      CALL frame
!
!-----------------------------------------------------------------------
!
!  Plot albedo values
!
!-----------------------------------------------------------------------
!
      CALL set( 0.,1.0,0.0,1.0,                                         &
           -.1,1.1,-.1,1.1,1)
      CALL wtstr(0.1,1.05,satname,16,0,0)
      CALL wtstr(0.4,1.05,timplt,16,0,0)
      CALL wtstr(0.8,1.05,fldname(k),16,0,0)
      CALL set( 0.,1.0,0.0,1.0,                                         &
               swx,nex,swy,ney,1)
      DO j=1,ny,5
        DO i=1,nx,5
          IF(satfld(i,j,k) > 0.0 .AND. satfld(i,j,k) < 2.0 ) THEN
            iplt=nint(100.*satfld(i,j,k))
            IF(iplt > 3 .AND. iplt < 10) THEN
              WRITE(chplt,820) iplt
              x=swx+FLOAT(i-1)*dx
              y=swy+FLOAT(j-1)*dy
              CALL wtstr(x,y,chplt,8,0,0)
            END IF
          END IF
        END DO
      END DO
      CALL pltmap(maxpts,mapfile,latmap,lonmap,xmap,ymap)
      CALL frame
    END IF
!
  END DO
!
  CALL clsgks
!
  GOTO 101
  
  100 CONTINUE
  
  WRITE(6,'(/a,a)') 'Error reading NAMELIST file. The program will abort.'
  
  101 CONTINUE
  
  STOP
END PROGRAM pltsatfld

SUBROUTINE satread(nx,ny,varname,DATA)
  IMPLICIT NONE
  INTEGER :: nx,ny
  REAL :: DATA(nx,ny)
  CHARACTER (LEN=6)   :: varname
  CHARACTER (LEN=256) :: fname
  INTEGER :: nxin,nyin
  INTEGER :: i,j

  fname='jian.9804201745.goes08.albedo'
  OPEN(3,FILE=fname,STATUS="old",FORM="unformatted")
  READ (3) varname
  PRINT *, ' varname=',varname
  READ (3) nxin,nyin
  PRINT *, ' nxin,nyin: ',nxin,nyin
  DO j=1,ny
!    print *, 'j=',j
    READ (3) (DATA(i,j),i=1,nx)
  END DO
  CLOSE (3)
  PRINT *, ' data( 1, 1): ',DATA( 1, 1)
  PRINT *, ' data(nx,ny): ',DATA(nx,ny)
  RETURN
END SUBROUTINE satread
