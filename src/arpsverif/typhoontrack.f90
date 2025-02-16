PROGRAM typhoontrack
!#######################################################################
!
! PURPOSE:
!
!   Calculate typhone track based on 2D arbitrary outputs from arpspost.
!   The output file can be plot using GrADS script developed by Yuanbing.
!
!-----------------------------------------------------------------------
!
! AUTHOR:
!
! Lingkun Ran (10/06/2010)
! It was hardcoded for typhoon Morakot.
!
! MODIFICATIONS:
!
! Yunheng Wang (10/07/2010)
! Rewrote for general application and following the ARPS coding conventions.
! Added comments and reorganized for readability.
! Incorporated into the ARPS system formally.
!
!-----------------------------------------------------------------------
!
!
  USE module_arbitrary_vario

  IMPLICIT NONE

  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'

  CHARACTER(LEN=256) :: grdbasfn
  INTEGER :: nprocx_in, nprocy_in
  INTEGER :: hinfmt

  NAMELIST /base_file/ grdbasfn, nprocx_in, nprocy_in, hinfmt

  INTEGER, PARAMETER :: max_path = 6
  INTEGER, PARAMETER :: max_typh = 5       ! maximum number of typhoon
  INTEGER, PARAMETER :: max_time = 200

  INTEGER :: numPath, timeBgn, timeIntvl, timeEnd
  CHARACTER(LEN=256) :: dirname2d
  INTEGER :: finfmt2d

  INTEGER :: pathOpts(max_path)

  NAMELIST /path_option/numPath, timeBgn, timeIntvl, timeEnd, pathOpts, &
                        dirname2d,finfmt2d

  INTEGER :: obsNumb(max_typh)
  REAL    :: obsLats(max_time*max_typh), obsLons(max_time*max_typh),    &
             obsMinp(max_time*max_typh), obsMaxw(max_time*max_typh)

  NAMELIST /obs_option/obsNumb,obsLats, obsLons, obsMinp, obsMaxw

  CHARACTER(LEN=256) :: outdir, outfile, outmslp, outmaxw,arpsdir

  NAMELIST /output/outdir, outfile, outmslp, outmaxw, arpsdir, lvldbg

!-----------------------------------------------------------------------

  CHARACTER(LEN=256) :: filename

  INTEGER :: nx,ny,nz, nzsoil, nstyps, nt

  REAL, ALLOCATABLE :: x(:), y(:), z(:), xs(:),ys(:)
  REAL, ALLOCATABLE :: lat(:,:), lon(:,:)
  REAL, ALLOCATABLE :: hterain(:,:)

  REAL, ALLOCATABLE :: tem1(:,:,:)

  REAL    :: latnot(2)
  REAL    :: ctrx, ctry, swx, swy
  INTEGER :: lengbf

  REAL    :: latmin, latmax, lonmin, lonmax

  INTEGER :: istatus

!-----------------------------------------------------------------------

  INTEGER :: ordMSLP, ordHGT7, ordU10m, ordV10m, ordU850, ordV850,      &
             ordU700, ordV700
  INTEGER :: ordW10m, ordW850, ordW700

  CHARACTER(LEN=6)   :: var2dnames(8)
  CHARACTER (LEN=40) :: var_name,varunits

  INTEGER :: ivar, numvar, numwind, iter

  INTEGER :: timelength

  real ::    dskm             ! the horizontal grid space of arps domain (Unit: km)
  real ::    time

  INTEGER :: i,j,k,it,ii,jj

  INTEGER :: irad, imin, imax, jmin, jmax,ifail,     &
             irad0,imin0,imax0,jmin0,jmax0

  integer :: ix, jy, ix1,jy1,ix2,jy2

  REAL    :: slpmin,slpmin0,slpmax0
  REAL    :: ght700min,ght700min0,ght700max0
  REAL    :: wind10min,wind10min0, wind10max0
  REAL    :: wind850min,wind850min0,wind850max0
  REAL    :: wind700min,wind700min0,wind700max0

  REAL, ALLOCATABLE :: var2d(:,:,:)

  REAL, ALLOCATABLE :: lat0(:,:),lon0(:,:),op(:,:),ow(:,:)

  REAL, ALLOCATABLE :: wind(:,:,:)
                             ! nz = 1, the wind speed at the first level about ground
                             ! nz = 2, the wind speed at 850 hPa level
                             ! nz = 3, the wind speed at 700 hPa level
  REAL, ALLOCATABLE :: swmax(:)
  REAL, ALLOCATABLE :: rlkslp(:)

  REAL, ALLOCATABLE :: tclat(:,:)
                             !tcmax = 1, the typhoon track derived from mean sea level pressure
                             !tcmax = 2, the typhoon track derived from hgt at 700 hPa
                             !tcmax = 3, the typhoon track derived from min wind at 10 m AGL
                             !tcmax = 4, the typhoon track derived from min wind at 850 hPa
                             !tcmax = 5, the typhoon track derived from min wind at 700 hPa
                             !tcmax = 6, the typhoon track derived from the combination of tcmax = 1 and tcmax = 2
  REAL, ALLOCATABLE :: tclon(:,:)


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

 CALL mpinit_var  ! initialize MPI variable for other part of the ARPS subroutines

!-----------------------------------------------------------------------
!
! Reading namelist variables
!
!-----------------------------------------------------------------------

  grdbasfn  = './'
  hinfmt    = 3
  nprocx_in = 1
  nprocy_in = 1

  READ(5,NML=base_file)
  WRITE(6,'(1x,a)') 'Namelist block base_file successfully read.'

  numPath = 6
  timeBgn = 0
  timeIntvl = 3600
  timeEnd = 79200
  pathOpts = (/1,2,3,4,5,6/)

  READ(5,NML=path_option)
  WRITE(6,'(1x,a)') 'Namelist block path_option successfully read.'

  obsNumb = 0
  obsLats = 0.0
  obsLons = 0.0
  obsMinp = 0.0
  obsMaxw = 0.0

  !obsNumb(1) = 17
  !obsLats = (/23.5,23.4,23.4,23.5,23.6,23.9,24.0,24.3,24.6,24.9,25.5,25.3,25.4,25.6,25.7,25.8,25.9/)
  !obsLons = (/123.2,122.6,122.4,122.3,122.1,121.9,121.6,121.5,121.5,121.3,121.0,120.8,120.8,120.8,120.7,120.6,120.5/)
  !obsMinp = (/955.,955.,955.,955.,955.,955.,965.,965.,970.,975.,975.,975.,975.,975.,975.,975.,975./)
  !obsMaxw = (/40.,40.,40.,40.,40.,40.,38.,38.,38.,35.,35.,35.,35.,35.,35.,35.,33./)

  READ(5,NML=obs_option)
  WRITE(6,'(1x,a)') 'Namelist block obs_option successfully read.'

  outdir = './'
  outfile = 'path.dat'
  outmslp = ' '
  outmaxw = ' '
  lvldbg = 0

  READ(5,NML=output)
  WRITE(6,'(1x,a)') 'Namelist block output successfully read.'

!-----------------------------------------------------------------------
!
! Get Observations
!
!-----------------------------------------------------------------------

  ALLOCATE(lat0(MAXVAL(obsNumb),SIZE(obsNumb)), STAT = istatus)
  ALLOCATE(lon0(MAXVAL(obsNumb),SIZE(obsNumb)), STAT = istatus)
  ALLOCATE(op  (MAXVAL(obsNumb),SIZE(obsNumb)), STAT = istatus)
  ALLOCATE(ow  (MAXVAL(obsNumb),SIZE(obsNumb)), STAT = istatus)

  jj = 0
  DO iter = 1, SIZE(obsNumb)
    DO ii = 1, obsNumb(iter)
      jj = jj + 1
      lat0(ii,iter) = obsLats(jj)
      lon0(ii,iter) = obsLons(jj)
      op(ii,iter) = obsMinp(jj)
      ow(ii,iter) = obsMaxw(jj)
    END DO
  END DO

!-----------------------------------------------------------------------
!
! Get arps dimensions
!
!-----------------------------------------------------------------------

  nt = (timeEnd-timeBgn)/timeIntvl + 1

  IF(nprocx_in > 1 .OR. nprocy_in > 1) THEN
    CALL gtsplitfn(grdbasfn,1,1,1,1,1,1,0,0,1,0,filename,istatus)
  ELSE
    WRITE(filename,'(a)') TRIM(grdbasfn)
  END IF
  lengbf = len_trim(filename)

  CALL get_dims_from_data(hinfmt,filename(1:lengbf),       &
          nx,ny,nz,nzsoil,nstyps,istatus)

  IF( istatus /= 0 ) THEN
    PRINT*,'Problem occured when trying to get dimensions from data.'
    PRINT*,'Program stopped.'
    CALL arpsstop('get_dims_from_data error.',1)
  END IF

  IF(nprocx_in > 1 .OR. nprocy_in > 1) THEN
    nx = (nx-3)*nprocx_in + 3
    ny = (ny-3)*nprocy_in + 3
  END IF

  ALLOCATE(x(nx), STAT = istatus)
  ALLOCATE(y(ny), STAT = istatus)
  ALLOCATE(z(nz), STAT = istatus)

  ALLOCATE(hterain(nx,ny), STAT = istatus)
  ALLOCATE(lat(nx,ny), STAT = istatus)
  ALLOCATE(lon(nx,ny), STAT = istatus)

  ALLOCATE(xs(nx), STAT = istatus)
  ALLOCATE(ys(ny), STAT = istatus)

  ALLOCATE(tem1(nx,ny,nz), STAT = istatus)

!-----------------------------------------------------------------------
!
! Get map projection and grid information as well terrain height
!
!-----------------------------------------------------------------------

  CALL get_gridxyzzp(nx,ny,nz,filename,hinfmt,nprocx_in,nprocy_in,      &
                     x,y,z,tem1,istatus)

  hterain(:,:) = tem1(:,:,2)
  dx = x(2)-x(1)
  dy = y(2)-y(1)
!
!-----------------------------------------------------------------------
!
!  Establish coordinate for ARPS scalar fields.
!
!-----------------------------------------------------------------------
!

  DO i=1,nx-1
    xs(i)=0.5*(x(i)+x(i+1))
  END DO
  xs(nx)=2.*xs(nx-1)-xs(nx-2)

  DO j=1,ny-1
    ys(j)=0.5*(y(j)+y(j+1))
  END DO
  ys(ny)=2.*ys(ny-1)-ys(ny-2)

!-----------------------------------------------------------------------
!
!  Set ARPS map projection
!
!-----------------------------------------------------------------------

  latnot(1) = trulat1
  latnot(2) = trulat2
  !print *, mapproj, sclfct,latnot,trulon,ctrlat,ctrlon
  print *,nx,ny,dx,dy
  CALL setmapr(mapproj,sclfct,latnot,trulon)
  CALL lltoxy( 1,1, ctrlat,ctrlon, ctrx, ctry )
  swx = ctrx - (FLOAT((nx-3))/2.) * dx*sclfct
  swy = ctry - (FLOAT((ny-3))/2.) * dy*sclfct
  CALL setorig(1,swx,swy)

  CALL xytoll(nx,ny,xs,ys,lat,lon)

!-----------------------------------------------------------------------
!
! Check variables to be read
!
!-----------------------------------------------------------------------

  ordMSLP = 0
  ordHGT7 = 0
  ordU10m = 0
  ordV10m = 0
  ordU850 = 0
  ordV850 = 0
  ordU700 = 0
  ordV700 = 0

  numvar = 0
  DO ivar = 1, numPath
    SELECT CASE (pathOpts(ivar))
    CASE (1,2,6)
      IF (ordMSLP < 1) THEN
        numvar  = numvar+1
        ordMSLP = numvar
        var2dnames(numvar) = 'mspres'
      END IF
      IF (ordHGT7 < 1) THEN
        numvar  = numvar+1
        ordHGT7 = numvar
        var2dnames(numvar) = 'hgt700'
      END IF
    CASE (3)
      IF (ordU10m < 1) THEN
        numvar  = numvar+1
        ordU10m = numvar
        var2dnames(numvar) = 'u10m__'
      END IF

      IF (ordV10m < 1) THEN
        numvar  = numvar+1
        ordV10m = numvar
        var2dnames(numvar) = 'v10m__'
      END IF

    CASE (4)
      IF (ordU850 < 1) THEN
        numvar  = numvar+1
        ordU850 = numvar
        var2dnames(numvar) = 'u850__'
      END IF

      IF (ordV850 < 1) THEN
        numvar  = numvar+1
        ordV850 = numvar
        var2dnames(numvar) = 'v850__'
      END IF

    CASE (5)
      IF (ordU700 < 1) THEN
        numvar  = numvar+1
        ordU700 = numvar
        var2dnames(numvar) = 'u700__'
      END IF

      IF (ordV700 < 1) THEN
        numvar  = numvar+1
        ordV700 = numvar
        var2dnames(numvar) = 'v700__'
      END IF
    CASE DEFAULT
      WRITE(6,'(a,I2)') 'Wrong option ',pathOpts(ivar)
      CALL arpsstop('Wrong option',1)
    END SELECT
  END DO

!-----------------------------------------------------------------------

  ALLOCATE(var2d(nx,ny,numvar), STAT = istatus)
  ALLOCATE(wind(nx,ny,3),       STAT = istatus)
  ALLOCATE(tclat(max_path,nt),  STAT = istatus)
  ALLOCATE(tclon(max_path,nt),  STAT = istatus)

  ALLOCATE(rlkslp(nt), STAT = istatus)
  ALLOCATE(swmax(nt),  STAT = istatus)

  dskm = MAX(dx,dy) / 1000.0

  lengbf = INDEX(grdbasfn,'.')
  !runname = grdbasfn(1:lengbf-1)

  DO it=1,nt

    timelength = (it-1)*timeIntvl + timeBgn
    time = float(timelength)

    DO ivar=1,numvar

      CALL arbvar_init(0,finfmt2d,iter,1,1,'',                          &
                 dirname2d,runname,var2dnames(ivar),time,lvldbg,istatus)

      CALL arbvar_read(nx,ny,1,var2d(:,:,ivar),var2dnames(ivar),        &
                       var_name,varunits,lvldbg,istatus)

      CALL arbvar_exit(lvldbg,istatus)

    END DO

    ordW10m = 0
    ordW850 = 0
    ordW700 = 0

    numwind = 0
    IF (ordU10m > 0) THEN
      numwind = numwind + 1
      ordW10m = numwind
      DO j=1,ny
        DO i=1,nx
          wind(i,j,ordW10m) = SQRT(var2d(i,j,ordU10m)**2 + var2d(i,j,ordV10m)**2)
        END DO
      END DO
    END IF

    IF (ordU850 > 0) THEN
      numwind = numwind + 1
      ordW850 = numwind
      DO j=1,ny
        DO i=1,nx
          wind(i,j,ordW850) = SQRT(var2d(i,j,ordU850)**2 + var2d(i,j,ordV850)**2)
        END DO
      END DO
    END IF

    IF (ordU700> 0) THEN
      numwind = numwind + 1
      ordW700 = numwind
      DO j=1,ny
        DO i=1,nx
          wind(i,j,ordW700) = SQRT(var2d(i,j,ordU700)**2 + var2d(i,j,ordV700)**2)
        END DO
      END DO
    END IF

    WRITE(6,'(8f10.3)')var2d(nx/2,ny/2,ordMSLP),var2d(nx/2,ny/2,ordHGT7),  &
                       var2d(nx/2,ny/2,ordU700),var2d(nx/2,ny/2,ordV700),  &
                       wind(nx/2,ny/2,1:numwind)

    DO i=1,5
      CALL smooth9p3d(var2d(:,:,ordMSLP),nx,ny)
    END DO

    slpmin     = 1004.0
    ght700min  = 5000.0
    wind10min  = 20.0
    wind850min = 20.0
    wind700min = 20.0

    slpmin0 =   1004.0
    slpmax0 =   900.0

    ght700min0 = 5000.
    ght700max0 = 0.

    wind10min0 = 100.
    wind10max0 = 0.

    wind850min0 = 100.
    wind850max0 = 0.

    wind700min0 = 100.
    wind700max0 = 0.

    irad = max(1, nint (7.5/dskm))          ! define the scan radius (grid number)
    !write(6,'(1x,a,7I4)') '===',numPath,pathOpts
    DO iter = 1, numPath
      SELECT CASE (pathOpts(iter))
      CASE (1)
        CALL find_center(nx,ny,irad,var2d(:,:,ordMSLP),5.,slpmin,       &
                         wind,numwind,slpmin0,slpmax0,                  &
                         ix1,jy1,istatus)

        tclat(1,it)=lat(ix1,jy1)
        tclon(1,it)=lon(ix1,jy1)
      CASE (2)
        CALL find_center(nx,ny,irad,var2d(:,:,ordHGT7),15.,ght700min,   &
                         wind,numwind,ght700min0,ght700max0,            &
                         ix2,jy2,istatus)

        tclat(2,it)=lat(ix2,jy2)
        tclon(2,it)=lon(ix2,jy2)

      CASE (6)   ! it must be done after 1 & 2
        ix = INT( (ix1+ix2)/2 )
        jy = INT( (jy1+jy2)/2 )
        tclat(6,it) = lat(ix,jy )
        tclon(6,it) = lon(ix,jy )

      CASE (3)
        CALL find_center(nx,ny,irad,wind(:,:,ordW10m),4.,wind10min,     &
                         wind,numwind,wind10min0,wind10max0,            &
                         ix,jy,istatus)
        tclat(3,it)=lat(ix,jy)
        tclon(3,it)=lon(ix,jy)

      CASE (4)
        CALL find_center(nx,ny,irad,wind(:,:,ordW850),4.,wind850min,    &
                         wind,numwind,wind850min0,wind850max0,          &
                         ix,jy,istatus)
        tclat(4,it)=lat(ix,jy)
        tclon(4,it)=lon(ix,jy)

      CASE (5)
        CALL find_center(nx,ny,irad,wind(:,:,ordW700),4.,wind700min,    &
                         wind,numwind,wind700min0,wind700max0,          &
                         ix,jy,istatus)
        Tclat(5,it)=lat(ix,jy)
        Tclon(5,it)=lon(ix,jy)
      END SELECT

    END DO

    irad0 = max(1, nint (300./dskm))         ! define the scan radius for calculating maximum wind spped in typhoon eye wall.
    imin0 = max(2,ix1-irad0)
    imax0 = min(nx-1,ix1+irad0)
    jmin0 = max(2,jy1-irad0)
    jmax0 = min(ny-1,jy1+irad0)

    swmax(it)= 0.0

    do ii= imin0,imax0
      do jj= jmin0,jmax0

      if(wind(ii,jj,ordW10m)> swmax(it)) then
        swmax(it)=wind(ii,jj,ordW10m)
      end if

      end do
    end do

    write(*,*) 'the max wind spped =',swmax(it)

    rlkslp(it)= var2d(ix1,jy1,ordMSLP)

  END DO   ! end of it  loop

!-----------------------------------------------------------------------
!
! write down the typhoon track
!
!-----------------------------------------------------------------------

  latmax = MAXVAL(lat)
  latmin = MINVAL(lat)
  lonmax = MAXVAL(lon)
  lonmin = MINVAL(lon)

  OPEN(20,FILE=TRIM(outdir)//TRIM(outfile),STATUS='UNKNOWN')

  WRITE(20,'(2a,4f8.3)') '# ',TRIM(runname),lonmin+2,lonmax-2,latmin+2,latmax-2

  numvar = 0
  IF (obsNumb(1) > 0) THEN
    numvar = numvar + 1
    !#1 Legeng marker(0=typhoon) color style'
    WRITE(20,'(a,I2.2,3(a,I2))') '#',numvar,' Observation ',            &
                         numvar-1,' ',numvar,' ',numvar
    write(20,'(1000(f8.3))')(lon0(it,1),it=1,obsNumb(1))    !observation of typhoon track
    write(20,'(1000(f8.3))')(lat0(it,1),it=1,obsNumb(1))
  END IF

  IF (ordMSLP > 0) THEN
    numvar = numvar + 1
    !#2 Legeng marker(0=typhoon) color style'
    WRITE(20,'(a,I2.2,3(a,I2))') '#',numvar,' MSLP ',                   &
                         numvar-1,' ',numvar,' ',numvar
    write(20,'(1000(f8.3))')(Tclon(1,it),it=1,nt) !typhoon track derived from slp
    write(20,'(1000(f8.3))')(Tclat(1,it),it=1,nt)
  END IF

  IF (ordHGT7 > 0) THEN
    numvar = numvar + 1
    !#3 Legeng marker(0=typhoon) color style'
    WRITE(20,'(a,I2.2,3(a,I2))') '#',numvar,' HGT700 ',                 &
                         numvar-1,' ',numvar,' ',numvar
    write(20,'(1000(f8.3))')(Tclon(2,it),it=1,nt) !typhoon track derived from ght at 700 hPa
    write(20,'(1000(f8.3))')(Tclat(2,it),it=1,nt)
  END IF

  IF (ordMSLP > 0 .AND. ordHGT7 > 0) THEN
    numvar = numvar + 1
    !#4 Legeng marker(0=typhoon) color style'
    WRITE(20,'(a,I2.2,3(a,I2))') '#',numvar,' Combination ',            &
                         numvar-1,' ',numvar,' ',numvar
    write(20,'(1000(f8.3))')(Tclon(6,it),it=1,nt) !typhoon track derived from combination of slp and ght at 700 hPa
    write(20,'(1000(f8.3))')(Tclat(6,it),it=1,nt)
  END IF

  CLOSE(20)

!-----------------------------------------------------------------------
!
! Output Geopotential Height and GrADS control files
!
!-----------------------------------------------------------------------

  CALL output_ght(nx,ny,hterain,outdir,                                 &
                  x,y,latmax,latmin,lonmax,lonmin,istatus)

  DEALLOCATE(x,y,z)
  DEALLOCATE(xs,ys)

  CALL run_plot(runname,arpsdir,outfile,outdir,istatus)

!-----------------------------------------------------------------------
!
! Output minimum sea level pressure if required
!
!-----------------------------------------------------------------------

  IF (LEN_TRIM(outmslp) > 0 .AND. outmslp /= ' ') THEN

    open(21,file=TRIM(outdir)//TRIM(outmslp),form='binary',status='unknown')

    do it =2,nt,3
      write(21)rlkslp(it)
    enddo

    IF (obsNumb(1) > 0) THEN
      do it =1,obsNumb(1)
        write(21) op(it,1)
      enddo
    END IF

    CLOSE(21)
  END IF

!-----------------------------------------------------------------------
!
! Output maximum wind speed if required
!
!-----------------------------------------------------------------------

  IF (LEN_TRIM(outmaxw) > 0 .AND. outmaxw /= ' ') THEN

    open(23,file='swmax.dat', form='binary',status='unknown')

    do it =2,nt,3
      write(23) swmax(it)
    enddo

    IF (obsNumb(1) > 0) THEN
      do it =1,obsNumb(1)
        write(23) ow(it,1)
       enddo
    END IF

    CLOSE(23)
  END IF

END PROGRAM typhoontrack

!#######################################################################

SUBROUTINE smooth9p3d(arr,ix,jx)

  IMPLICIT NONE
  INTEGER :: ix,jx
  REAL    :: arr(ix,jx)

!-----------------------------------------------------------------------

  REAL, ALLOCATABLE :: tem(:,:)

  INTEGER :: i,j, istatus
  REAL    :: wtf,wtfb,wtfc

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ALLOCATE(tem(ix,jx), STAT = istatus)

  wtf = 1.0/16.0

  DO j=1,jx
    DO i=1,ix
      tem(i,j)=arr(i,j)
    END DO
  END DO

  DO  j = 2, jx-2
    DO  i = 2, ix-2
      tem(i,j) = wtf*(arr(i-1,j-1)+2.*arr(i,j-1)+arr(i+1,j-1)           &
                        +2.*arr(i-1,j)+4.*arr(i,j)+2.*arr(i+1,j)        &
                        +arr(i-1,j+1)+2.*arr(i,j+1)+arr(i+1,j+1))
    END DO
  END DO

  DO  j = 2, jx-2
    DO  i = 2, ix-2
      arr(i,j) = tem(i,j)
    END DO
  END DO

  CALL edgfill(arr,1,ix,2,ix-2, 1,jx,1,jx-2, 1,1,1,1)

  DEALLOCATE(tem)

  RETURN
END SUBROUTINE smooth9p3d

!#######################################################################

SUBROUTINE find_center(nx,ny,irad,varin,thredhold,vartype,              &
                       wind, numwind,                                   &
                       inmin,inmax,ix,jy,istatus)
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: nx, ny
  INTEGER, INTENT(IN)  :: irad
  REAL,    INTENT(IN)  :: varin(nx,ny)
  INTEGER, INTENT(IN)  :: numwind
  REAL,    INTENT(IN)  :: wind(nx,ny,numwind)
  REAL,    INTENT(IN)  :: thredhold
  REAL,    INTENT(IN)  :: vartype   ! typical value of this field
  REAL,    INTENT(IN)  :: inmin, inmax
  INTEGER, INTENT(OUT) :: ix, jy
  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
  INTEGER :: i,j, ii, jj
  INTEGER :: imin, imax, jmin, jmax
  REAL    :: sumvar
  INTEGER :: icnt, ivar
  REAL    :: varmin, varmax, typemin
  LOGICAL :: realcenter

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0
  ix = 0
  jy = 0

  typemin = vartype
  DO j= 2,ny-1
    DO i= 2,nx-1

      imin = max(2,i-irad)
      imax = min(nx-1,i+irad)
      jmin = max(2,j-irad)
      jmax = min(ny-1,j+irad)

      varmin = inmin
      varmax = inmax
!write(6,*) i,j,varin(i,j),imin,imax,jmin,jmax
      do ii= imin,imax
        do jj= jmin,jmax

          if (varin(ii,jj) < varmin) then
            varmin = varin(ii,jj)
          endif
          if (varin(ii,jj) > varmax) then
            varmax = varin(ii,jj)
          endif

        enddo
      enddo
!write(6,*) varmin, varmax
      sumvar = 0.0
      icnt   = 0

      do ii= imin,imax
        do jj= jmin,jmax

          if (varmax - varmin <= thredhold ) then
            sumvar = sumvar+varin(ii,jj)
            icnt   = icnt+1
          endif

       enddo
      enddo
!write(6,*) icnt
      IF (icnt > 0) sumvar = sumvar / icnt

      realcenter = .TRUE.
      DO ivar = 1, numwind

        IF (wind(i,j,ivar) > 6.0) THEN
          realcenter = .FALSE.
          !WRITE(6,'(1x,a)') 'WARNING: Wind speed is > 6.0 in the center.'
          istatus = -1
        END IF
      END DO

      IF ( realcenter ) THEN
        IF (icnt > 0 .AND. sumvar < typemin) THEN
          typemin = sumvar
          ix=i
          jy=j
          !write(6,*) '***',ix,jy,wind(i,j,1:3),icnt,sumvar,vartype
        END IF
      END IF
!      write(6,'(a,2I4)') '***',ix,jy
     END DO
  END DO

! write(6,'(a,2I4)') '***',ix,jy
  RETURN
END SUBROUTINE find_center

!#######################################################################

SUBROUTINE output_ght(nx,ny,ht,outdir,                                  &
                      x,y,latmax,latmin,lonmax,lonmin,                  &
                      istatus)
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: nx, ny
  REAL,    INTENT(IN)  :: ht(nx,ny)
  REAL,    INTENT(IN)  :: x(nx), y(ny)
  REAL,    INTENT(IN)  :: latmax, latmin, lonmax, lonmin
  CHARACTER(LEN=*), INTENT(IN) :: outdir

  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'

!-----------------------------------------------------------------------
  CHARACTER (LEN=3) :: monnam(12)            ! Name of months
  DATA monnam/'jan', 'feb', 'mar', 'apr', 'may', 'jun',                 &
              'jul', 'aug', 'sep', 'oct', 'nov', 'dec'/

  INTEGER :: nchout0
  CHARACTER(LEN=256) :: filename
  CHARACTER(LEN=80)  :: chrstr

  REAL :: xbgn, ybgn, xinc, yinc
  REAL :: latinc, loninc
  REAL :: lat11, lon11

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  xbgn = (x(1) + x(2))/2.
  ybgn = (y(1) + y(2))/2.

  xinc = (x(2) - x(1))
  yinc = (y(2) - y(1))

  CALL xytoll(1,1,xbgn,ybgn,lat11,lon11)

  latinc = (latmax-latmin)/(ny-1)
  loninc = (lonmax-lonmin)/(nx-1)

!-----------------------------------------------------------------------

  filename = TRIM(outdir)//TRIM(runname)//'_ght'

  CALL getunit (nchout0)

  OPEN(nchout0,file=TRIM(filename)//'.dat',form='binary',status='unknown')
  WRITE(nchout0) ht
  CLOSE(nchout0)

!-----------------------------------------------------------------------

  WRITE (6,'(1x,a)') 'The GrADS data control file for GHT is '          &
                    //TRIM(filename)//'.ctl'

  OPEN (nchout0, FILE = TRIM(filename)//'.ctl', STATUS = 'unknown')

  WRITE (nchout0,'(2a)') 'DSET    ',TRIM(filename)//'.dat'
  WRITE (nchout0,'(a)')                                                 &
        'TITLE   ARPS model typhoon track output for '//TRIM(runname)
  !WRITE (nchout0,'(a)')                                                 &
  !      'OPTIONS sequential cray_32bit_ieee big_endian'
  WRITE (nchout0,'(a)')                                                 &
        'OPTIONS cray_32bit_ieee big_endian'
  !WRITE (nchout0,'(a,i10)') 'FILEHEADER ',0
  WRITE (nchout0,'(a)')   'UNDEF   -9.e+33'

  IF ( mapproj == 2 ) THEN
    WRITE (nchout0,'(a,1x,i8,1x,i3,a,2f12.6,2i3,3f12.6,2f12.2)')        &
        'PDEF',nx,ny,' LCC',lat11,lon11,1,1,                            &
            trulat1,trulat2,trulon,xinc,yinc
    WRITE (nchout0,'(a,1x,i8,a,f10.4,1x,f10.4)')                        &
        'XDEF',nx,'  LINEAR  ',lonmin,loninc
    WRITE (nchout0,'(a,1x,i8,a,f10.4,1x,f10.4)')                        &
        'YDEF',ny,'  LINEAR  ',latmin,latinc
  ELSE

    WRITE (nchout0,'(a)')                                             &
        '* For i-j-k display, uncomment the following 3 lines.'

    WRITE (nchout0,'(a,1x,i8,a,2i10)') '* XDEF',nx,'  LINEAR  ',1,1

    WRITE (nchout0,'(a,1x,i8,a,2i10)') '* YDEF',ny,'  LINEAR  ',1,1

    WRITE (nchout0,'(a)')                                             &
        '* For x-y-z display, uncomment the following 3 lines.'

    WRITE (nchout0,'(a,1x,i8,a,2f15.4)')                              &
        'XDEF',nx,'  LINEAR  ',xbgn/1000.,xinc/1000.

    WRITE (nchout0,'(a,1x,i8,a,2f15.4)')                              &
        'YDEF',ny,'  LINEAR  ',ybgn/1000.,yinc/1000.

  END IF

  WRITE (nchout0,'(a,2i10)')  'ZDEF        1  LINEAR  ',1,1

  WRITE (chrstr,'(i2.2,a,i2.2,a,i2.2,a3,i4.4)')                       &
          hour,':',minute,'Z',day,monnam(month),year

  WRITE (nchout0,'(3a)')'TDEF        1  LINEAR  ',TRIM(chrstr),' 1mn'
  WRITE (nchout0,'(a)') 'VARS 1'
  WRITE (nchout0,'(a)')                                               &
           'hgt        1 0 Geopotential height above mean sea level (m)'
  WRITE (nchout0,'(a)') 'ENDVARS'

  CLOSE (nchout0)
  CALL retunit(nchout0)

END SUBROUTINE output_ght

!#######################################################################

SUBROUTINE run_plot(runname,arpsdir,pathfile,outdir,istatus)

  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN)  :: runname, arpsdir,pathfile,outdir
  INTEGER,          INTENT(OUT) :: istatus
!-----------------------------------------------------------------------
  CHARACTER(LEN=512) :: chrstr

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  WRITE(chrstr,'(4a)') 'cp ',TRIM(arpsdir),'/scripts/lujing.gs ',TRIM(outdir)
  WRITE(6,'(1x,3a)') 'Running command ''',TRIM(chrstr),''' ...'
  CALL unixcmd(TRIM(chrstr))

  WRITE(chrstr,'(3a)') 'gradsc -lbc "run lujing.gs ',TRIM(pathfile),' 4"'
  WRITE(6,'(1x,3a)') 'Running command ''',TRIM(chrstr),''' ...'
  CALL unixcmd(TRIM(chrstr)//'> /dev/null')

!  CALL unixcmd('rm lujing.gmf > /dev/null')

  RETURN
END SUBROUTINE run_plot
