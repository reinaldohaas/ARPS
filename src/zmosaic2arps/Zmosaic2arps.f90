PROGRAM Zmosaic2arps
!
!
!-----------------------------------------------------------------------
!
!  Program for generating a reflectivity file in ARPS grid from
!        NSSL Mosaica reflectivity 
!
!  Author: Ming Hu, CAPS. University of Oklahma.
!  First written: 04/06/2006.
!
!  History:
!
!    5/15/2007 Fanyou Kong
!      Rewrite to remove mosaicfile, and to include 2D field read
!        (e.g., cref, 1h_rad_hsr, etc)
!      Also remove NAMELIST jobname (not used)
!      Add NetCDF file existence check, and uncompress if applicable
!
!    O5/16/2007 Y. Wang
!    Added check for "/" at the end of mosaicPath.
!    Changed wrtvar1 to wrtvar2.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!
!  ARPS grid
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'mp.inc'
  INCLUDE 'netcdf.inc'

!--------------------------------------------------------------
!  ARPS grid 
!--------------------------------------------------------------
  REAL, allocatable :: lath(:,:)    ! Latitude of each terrain point
  REAL, allocatable :: lonh(:,:)    ! Longitude of each terrain point
  REAL, allocatable :: zp(:,:,:)    !  Physical height coordinate defined at
                                    ! w-point of the staggered grid.

  REAL, allocatable :: ref3d(:,:,:) !3D reflectivity in mosaic vertical grid
  REAL, allocatable :: ref_arps_3d(:,:,:)  ! 3D reflectivity in arps grid

  REAL, allocatable :: var2d(:,:)   ! 2D variable
  REAL, allocatable :: compref(:,:) ! composite reflectivity
  REAL, allocatable :: hsr_1hr(:,:) ! 1-h accumulate precipitation

!--------------------------------------------------------------
!  Namelist
!--------------------------------------------------------------
!
! grid 
!
  INTEGER :: nx, ny, nz
!  NAMELIST /jobname/ runname

  NAMELIST /grid/ nx,ny, nz, dx,dy,dz,                              &
            strhopt,dzmin,zrefsfc,dlayer1,dlayer2,strhtune,zflat,   &
            ctrlat,ctrlon
  NAMELIST /projection/ mapproj,trulat1,trulat2,trulon,sclfct
  NAMELIST /terrain/ ternopt,terndta,ternfmt
!
!  For reflectiivty mosaic
!
  CHARACTER*256 mosaicPath
  CHARACTER*15 mosaictime(50)
  CHARACTER*256 mosaicTile(14)
  INTEGER :: tversion,mosaic_opt,ifcompositeRef,numvolume,ntiles
  NAMELIST /mosaicpara/mosaicPath,tversion,mosaic_opt, &
                       ifcompositeRef,numvolume,       &
                       mosaictime,ntiles,mosaicTile

  INTEGER :: dmpfmt
  NAMELIST /output/ dirname, dmpfmt

!--------------------------------------------------------------
!  Grid of NSSL mosaic 
!--------------------------------------------------------------
  INTEGER ::   mscNlon   ! number of longitude of mosaic data
  INTEGER ::   mscNlat   ! number of latitude of mosaic data
  INTEGER ::   mscNlev   ! number of vertical levels of mosaic data
  REAL, allocatable :: msclon(:)        ! longitude of mosaic data
  REAL, allocatable :: msclat(:)        ! latitude of mosaic data
  REAL, allocatable :: msclev(:)        ! level of mosaic data
  REAL, allocatable :: mscValue(:,:)   ! reflectivity

  REAL :: lonMin,latMin,lonMax,latMax,dlon,dlat

!-----------------------------------------------------------------------
!
!  LOCAL VARIABLES:
!
!-----------------------------------------------------------------------
!
  CHARACTER*256 tempfile
  INTEGER :: NCID,iSTATUS

  LOGICAL :: fexist

  INTEGER :: maxlvl

  REAL :: swx,swy,ctrx,ctry    ! Temporary variables.

  INTEGER :: i,j,k,ifl,im,ifn

  REAL :: rlon,rlat
  INTEGER  :: ip,jp
  REAL ::  rip,rjp
  REAL ::  dip,djp
  REAL ::  w1,w2,w3,w4
  REAL ::  ref1,ref2,ref3,ref4
  INTEGER  :: iyr,imon,iday,ihr,imin,isec

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  CALL mpinit_var          ! calling of writtrn needs some variables
!
!-----------------------------------------------------------------------
!
!  Set up the default values for all the variables to be read in
!  using the namelist method. In case the user does not specify a
!  particular value, this value will be used.
!
!-----------------------------------------------------------------------
!
  hdfcompr = 5

  nx = 67
  ny = 67
  nz = 53
!  runname = 'may20'
  dx = 1000.0
  dy = 1000.0
  dz = 400.0

  ctrlat  =  35.0
  ctrlon  = -100.0

  mapproj  = 0
  trulat1  =   30.0
  trulat2  =   60.0
  trulon   = -100.0
  sclfct   =    1.0
!
!-----------------------------------------------------------------------
!
!  Read and Print of Namelist variables and parameters.
!
!-----------------------------------------------------------------------
!
!  READ (5,jobname,ERR=980,END=990)
!  WRITE(6,'(/a,a)') 'The name of this run is: ', runname

  READ (5,mosaicpara,ERR=980,END=990)
  WRITE(6,*) 'NetCDF data type: ', mosaic_opt
  WRITE(6,*) 'Number of times: ', numvolume 
  WRITE(6,*) 'Number of tiles: ', ntiles

  READ(5,grid,ERR=980,END=990)

  READ(5,projection,ERR=980,END=990)

  READ(5,terrain,ERR=980,END=990)

  READ(5,output,ERR=980,END=990)

  lfnkey = LEN_TRIM(mosaicpath)
  IF (mosaicpath(lfnkey:lfnkey) /= '/') mosaicpath = TRIM(mosaicpath) // '/'
!
!  CALL gtlfnkey( runname, lfnkey )
!
!  get the latitude and longitude of ARPS domain
!
  allocate(lath(0:nx,0:ny))
  allocate(lonh(0:nx,0:ny))
  allocate(compref(nx,ny))
  allocate(hsr_1hr(nx,ny))

  CALL getcoordinate(nx,ny,lath,lonh)

  CALL a3dmax0(lath,0,nx,0,nx,0,ny,0,ny,1,1,1,1,latmax,latmin)
  CALL a3dmax0(lonh,0,nx,0,nx,0,ny,0,ny,1,1,1,1,lonmax,lonmin)
!  write(*,*) lonmax,lonmin
!  write(*,*) latmax,latmin
  IF(mosaic_opt == 1) THEN
    allocate(zp(nx,ny,nz))
    allocate(ref_arps_3d(nx,ny,nz))
    CALL getVertcoordinate(nx,ny,nz,zp)
  END IF
!
!  begin to read in mosaic file
!
nv: DO ifl=1, numvolume
!
!  get dimension of mosaic field and allocate arrays
!
    DO im=1, ntiles
      tempfile=trim(mosaicPath)//trim(mosaictime(ifl))//trim(mosaicTile(im))
      write(*,*) 'process file No.',im,'with name:',trim(tempfile)

      INQUIRE(FILE=trim(tempfile), EXIST = fexist )
      IF( .NOT.fexist) THEN
        INQUIRE(FILE=trim(tempfile)//'.gz', EXIST = fexist )
          IF(fexist) THEN
            CALL uncmprs(trim(tempfile)//'.gz')
          ELSE
            PRINT *,trim(tempfile),' or its .gz file does not exist, SKIP'
            CYCLE nv
          END IF
      END IF

       IF( tversion == 14 ) then
         call GET_DIM_ATT_Mosaic(tempfile,mscNlon,mscNlat,mscNlev, &
                   lonMin,latMin,lonMax,latMax,dlon,dlat)
       ELSEIF( tversion == 8 ) then
         call GET_DIM_ATT_Mosaic8(tempfile,mscNlon,mscNlat,mscNlev, &
                   lonMin,latMin,lonMax,latMax,dlon,dlat,mosaic_opt)

       ELSE
         write(*,*) ' unknown tile version !!!'
         stop 123
       ENDIF

       allocate(msclon(mscNlon))
       allocate(msclat(mscNlat))
       allocate(msclev(mscNlev))
       allocate(mscValue(mscNlon,mscNlat))

       if(im == 1) then
         maxlvl=mscNlev
         allocate(ref3d(nx,ny,maxlvl))
         ref3d=-99999.0
       else
         if(maxlvl .ne. mscNlev) then
            write(*,*) ' The vertical dimensions are different in two tiles'
            stop
         endif
       endif

       msclon(1)=lonMin
       DO i=1,mscNlon-1
         msclon(i+1)=msclon(i)+dlon
       ENDDO
       msclat(1)=latMin
       DO i=1,mscNlat-1
         msclat(i+1)=msclat(i)+dlat
       ENDDO
!
!  ingest mosaic file and interpolation
!
       call OPEN_Mosaic(tempfile, NCID)

       if(tversion == 14 ) then
         call Check_DIM_ATT_Mosaic(NCID,mscNlon,mscNlat,mscNlev,  &
               lonMin,latMin,lonMax,latMax,dlon,dlat)
       elseif(tversion == 8 ) then
         call Check_DIM_ATT_Mosaic8(NCID,mscNlon,mscNlat,mscNlev,  &
               lonMin,latMin,lonMax,latMax,dlon,dlat,mosaic_opt)
       endif

       DO k=1, mscNlev

         IF(mosaic_opt == 1) THEN
         call  GET_Mosaic_sngl_Mosaic(NCID,mscNlon,mscNlat,k,mscValue)
         ELSE IF(mosaic_opt == 2) THEN
         call  GET_Mosaic_cref_Mosaic(NCID,mscNlon,mscNlat,mscValue)
         ELSE IF(mosaic_opt == 3) THEN
         call  GET_Mosaic_hsr_Mosaic(NCID,mscNlon,mscNlat,mscValue)
         ELSE
           print *,'mosaic_opt=',mosaic_opt,' is invalid, STOP!!!'
           stop
         END IF
        
!         DO j=0,ny
!         DO i=0,nx
         DO j=1,ny
         DO i=1,nx
            rlat=lath(i,j)
            rlon=lonh(i,j)
          
            if(tversion == 14 ) then
               rip=(rlon-lonMin)/dlon+1
               rjp=(rlat-latMin)/dlat+1
               ip=int(rip)
               jp=int(rjp)
               dip=rip-ip
               djp=rjp-jp
            elseif(tversion == 8 ) then
               rip=(rlon-lonMin)/dlon+1
               rjp=(latMax-rlat)/dlat+1
               ip=int(rip)
               jp=int(rjp)
               dip=rip-ip
               djp=rjp-jp
            else
               write(*,*) ' Unknown Mosaic format !!'
               stop 123
            endif

            if( ip >= 1 .and. ip < mscNlon ) then
            if( jp >= 1 .and. jp < mscNlat ) then
! inside mosaic domain
              w1=(1.0-dip)*(1.0-djp)
              w2=dip*(1.0-djp)
              w3=dip*djp
              w4=(1.0-dip)*djp
              ref1=mscValue(ip,jp)
              ref2=mscValue(ip+1,jp)
              ref3=mscValue(ip+1,jp+1)
              ref4=mscValue(ip,jp+1)
              if(ref1 > -500.0 .and. ref2 > -500.0 .and.  &
               ref3 > -500.0 .and. ref4 > -500.0 ) then
                 ref3d(i,j,k)=(ref1*w1+ref2*w2+ref3*w3+ref4*w4)/10.0
              endif
            endif
            endif
         ENDDO
         ENDDO
       ENDDO ! mscNlev
!
       call CLOSE_Mosaic(NCID)
    
      deallocate(msclon)
      deallocate(msclat)
      deallocate(msclev)
      deallocate(mscValue)

    ENDDO   ! ntiles 

    write(*,*)
    write(*,*) 'successfully process Mosaic data at ', &
               trim(mosaictime(ifl))
!    write(*,*)
!

    IF(mosaic_opt == 1) THEN
      CALL vert_interp_ref(nx,ny,nz,mscNlev,ref3d,ref_arps_3d,zp)
      WRITE(*,*) 'successfully get reflectivity in ARPS grid at ', &
                 trim(mosaictime(ifl))
      !
      !  dump out colume data for data analysis
      !
      READ(mosaictime(ifl),'(i4,i2,i2,1x,i2,i2)') iyr,imon,iday,ihr,imin
      isec = 0
      call wtradcol_Mosaic(nx,ny,nz,                       &
             mosaictime(ifl),iyr,imon,iday,ihr,imin,isec,  &
             lath,lonh,zp,ref_arps_3d)

      ref_arps_3d = -99999.0
      call readadcol_Mosaic(nx,ny,nz,mosaictime(ifl),ref_arps_3d)
      !
      ! get composite reflectivity
      !
      IF (ifcompositeRef==1) THEN
        compref=-99999.0
        DO k=2, nz-1
          DO j=1,ny
            DO i=1,nx
              compref(i,j)=max(compref(i,j),ref_arps_3d(i,j,k))
            ENDDO
          ENDDO
        ENDDO
        !
        !  dump out composite reflectivity
        !
        DO j=1,ny
          DO i=1,nx
            if( compref(i,j) < 0.0 ) compref(i,j)=0
            if( compref(i,j) > 100.0 )  write(*,*) i,j,compref(i,j)
          ENDDO
        ENDDO

        call wrtvar2(nx,ny,1,compref, 'compst','composite reflectivity', &
             'dBZ',0,mosaictime(ifl),dirname,dmpfmt,hdfcompr,0,iSTATUS)

      END IF  ! ifcompositeRef=1

    ELSE IF(mosaic_opt == 2) THEN
      DO j=1,ny
        DO i=1,nx
          compref(i,j) = max(ref3d(i,j,1), 0.0)
        ENDDO
      ENDDO
      call wrtvar2(nx,ny,1,compref, 'compst','composite reflectivity', &
           'dBZ',0,mosaictime(ifl),dirname,dmpfmt,hdfcompr,0,iSTATUS)
    ELSE IF(mosaic_opt == 3) THEN
      DO j=1,ny
        DO i=1,nx
          !hsr_1hr(i,j) = max(ref3d(i,j,1)/25.4, 0.0)  ! convert to inch
          hsr_1hr(i,j) = max(ref3d(i,j,1), 0.0) ! keep mm unit
        ENDDO
      ENDDO
      call wrtvar2(nx,ny,1,hsr_1hr, 'hsr1hr','1-h precipitation', &
           !'in',0,mosaictime(ifl),dirname,dmpfmt,0,0,iSTATUS)
           'mm',0,mosaictime(ifl),dirname,dmpfmt,hdfcompr,0,iSTATUS)
    ELSE
      print *,'mosaic_opt=',mosaic_opt,' is invalid, STOP!!!'
      stop
    END IF

    deallocate(ref3d)

  222 CONTINUE

  ENDDO nv    !numvolume
!  close(ifn)

  DEALLOCATE(lath)
  DEALLOCATE(lonh)
  DEALLOCATE(compref)
  DEALLOCATE(hsr_1hr)
  IF(mosaic_opt == 1) THEN
    DEALLOCATE(zp)
    DEALLOCATE(ref_arps_3d)
  END IF

  CALL exit(0)

  980   CONTINUE
  WRITE(6,'(/1x,a,a)')                                                  &
      'Error occured when reading namelist input file. Program stopped.'

  CALL exit(980)

  990   CONTINUE
  WRITE(6,'(/1x,a,a)')                                                  &
      'End of file reached when reading namelist input file. ',         &
      'Program stopped.'

  CALL exit(990)

  981   CONTINUE
  WRITE(6,'(/1x,a,a)')                                                  &
      'Error occured when reading mosaic file name. Program stopped.'

  CALL exit(981)

  991   CONTINUE
  WRITE(6,'(/1x,a,a)')                                                  &
      'End of file reached when reading mosaic file path. ',            &
      'Program stopped.'
  CALL exit(991)

  STOP
END PROGRAM Zmosaic2arps
