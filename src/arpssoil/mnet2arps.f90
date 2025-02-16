!
!##################################################################
!##################################################################
!######                                                      ######
!######                PROGRAM ARPSSOIL                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

PROGRAM arpssoil
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Process soil temperature and other soil obs and create
!  soil input file for ARPS.
!
!
!-----------------------------------------------------------------------
!
!  AUTHORS: Min Zou and Keith Brewster
!  May, 1994
!
!  MODIFICATION HISTORY:
!  Oct, 1994 (KB)
!  Changed input times to be more flexible.
!
!  Nov 14, 1994 (KB)
!  Write data in format expected for ARPS subroutine INITSFC
!
!-----------------------------------------------------------------------
!
!  Sizing parameters
!  Note in the following, nvar_max is the max of
!  nvar_ob and nvar_anx, and is used for the size of the
!  Barnes work arrays.
!  ntime is the number of time levels of obs read-in.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: nx, ny, nz               !  dimensions at x, y, z diretions.
  INTEGER :: ntime,nvar_ob,nvar_anx,nvar_max,nvar_out
  PARAMETER (ntime=1,                                                   &
             nvar_ob=3, nvar_anx=3, nvar_max=3, nvar_out=3)
  INTEGER :: mxstn,nsrc_sfc
  PARAMETER (mxstn=500,nsrc_sfc=4)
!
!-----------------------------------------------------------------------
!
!  Barnes analysis parameters
!
!-----------------------------------------------------------------------
!
  INTEGER :: npass
  REAL :: defrange,wlim
  INTEGER :: iqspr
  REAL :: sprdist,derrtol
  PARAMETER (npass = 3,                                                 &
             wlim = 1.e-06,          & ! minimum weight -- establishes max radius
         iqspr= 7,                   & ! Use qobs of pstn to combine x,y,elev
         sprdist = 15.,              & ! km -- closer pairs are made into superob
         derrtol = 20.)                ! km
  INTEGER :: klim(npass)
  DATA klim /1,2,2/
  REAL :: kappa(npass)
  DATA kappa /3.,1.,1./
  INTEGER :: irngsel(nvar_max)
  DATA irngsel /1,1,1/
!
!-----------------------------------------------------------------------
!
!  Barnes scratch space
!
!-----------------------------------------------------------------------
!
  INTEGER :: maxrng
  PARAMETER(maxrng=1)
  INTEGER :: knt(nvar_max)
  REAL :: rngwgt(maxrng)
  REAL :: wgtsum(nvar_max),zsum(nvar_max)
!
  REAL, ALLOCATABLE :: xgrid(:),ygrid(:)
!
!-----------------------------------------------------------------------
!
!  Barnes scale range
!
!-----------------------------------------------------------------------
!
  REAL :: rngobs(mxstn)
  REAL, ALLOCATABLE :: range(:,:)
!
!-----------------------------------------------------------------------
!
!  Primary analysis variables
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: anx(:,:,:)
!
  REAL, ALLOCATABLE :: tsfc(:,:)
  REAL, ALLOCATABLE :: tsoil(:,:)
!
!-----------------------------------------------------------------------
!
!  Primary observation variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: isrc(mxstn)
  REAL :: obs(mxstn,nvar_anx)
  REAL :: odiff(mxstn,nvar_anx)
  REAL :: obanx(mxstn,nvar_anx)

  CHARACTER (LEN=6) :: var_ob(nvar_ob)
  DATA var_ob / 'tsod05','tsod30','solrad'/
  CHARACTER (LEN=6) :: var_anx(nvar_anx)
  DATA var_anx/ 'tsod05','tsod30','solrad'/
!
!-----------------------------------------------------------------------
!
!  Quality Control Parameters and Variables
!
!-----------------------------------------------------------------------
!
  REAL :: rmiss
  PARAMETER (rmiss=-99.)
!  real qclim(nsrc_sfc,nvar_ob)
  INTEGER :: rely(mxstn,nvar_ob,ntime)
!  real qobsrd(mxstn,nvar_ob)
  REAL :: qobs(mxstn,nvar_anx)
!
!-----------------------------------------------------------------------
!
!  Climin is the minimum possible value of each observed variable.
!  Climax is the maximum possible value of each observed variable.
!
!-----------------------------------------------------------------------
!
  REAL :: climin(nvar_ob),climax(nvar_ob)
  DATA climin /  -50.,-50.,0./
  DATA climax /   50., 50.,1500./
!
!-----------------------------------------------------------------------
!
!  Source-dependent parameters
!  Qobs is the expected square error it is used for
!  setting the QC threshold (qclim) and for relative
!  weighting of observations.
!  Zero is not a valid qobs.
!
!-----------------------------------------------------------------------
!
  REAL :: qcmult
  PARAMETER(qcmult=3.0)             ! mult sqrt(qobs) by 3.0 to get qclim
  CHARACTER (LEN=8) :: name_src(nsrc_sfc)
  REAL :: qsrc(nsrc_sfc,nvar_ob)
  DATA name_src /'MESO','MOBLS','MOBLM','SA'/
  DATA qsrc     /   2.,    2.,    2.,    2.,                            &
                    2.,    2.,    2.,    2.,                            &
                  100.,  100.,  100.,  100./
!
!-----------------------------------------------------------------------
!
!  Read_obs parameter
!
!-----------------------------------------------------------------------
!
  INTEGER :: k_sta
  CHARACTER (LEN=4) :: xtnam(mxstn)
  REAL :: xtlat(mxstn),xtlon(mxstn),xtelev(mxstn),xts(mxstn,nvar_anx)
  REAL :: xtx(mxstn),xty(mxstn)

  CHARACTER (LEN=9) :: c9time

  CHARACTER (LEN=15) :: gradsctl
  CHARACTER (LEN=15) :: gradsfile

  CHARACTER (LEN=3) :: monnam(12)
  DATA monnam/'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',                 &
              'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'/
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
!
!-----------------------------------------------------------------------
!
!  Misc. internal variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,initsec,myr,istat,istatus,ierr
  INTEGER :: tsfcout,tsoilout,wsfcout,wdpout,wcanpout,idummy
  REAL :: lllat,lllon,lrlat,lrlon,urlat,urlon,ullat,ullon,rdummy
  REAL :: latnot(2)
  REAL :: x0,y0

  CHARACTER(LEN=256) :: namelist_filename
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  Call initpara to read in ARPS parameters of namelists
!
!-----------------------------------------------------------------------
!
  namelist_filename = ' '
  CALL initpara(nx,ny,nz,nstyps,namelist_filename)
!
  CALL ctim2abss(year,month,day,hour,minute,second,initsec)

!-----------------------------------------------------------------------
!
! Allocate the nx, ny, nz dependent arrays and initialize to zero
!
!-----------------------------------------------------------------------
  ALLOCATE( xgrid(nx), STAT=istatus)
  CALL check_alloc_status(istatus, "arpssoil:xgrid")
  xgrid=0
  ALLOCATE( ygrid(ny), STAT=istatus)
  CALL check_alloc_status(istatus, "arpssoil:ygrid")
  ygrid=0
  ALLOCATE( range(nx,ny), STAT=istatus)
  CALL check_alloc_status(istatus, "arpssoil:range")
  range=0
  ALLOCATE( anx(nx,ny,nvar_anx), STAT=istatus)
  CALL check_alloc_status(istatus, "arpssoil:anx")
  anx=0

  ALLOCATE( tsfc(nx,ny), STAT=istatus)
  CALL check_alloc_status(istatus, "arpssoil:tsfc")
  tsfc=0
  ALLOCATE( tsoil(nx,ny), STAT=istatus)
  CALL check_alloc_status(istatus, "arpssoil:tsoil")
  tsoil=0

  tsfc=anx(:,:,1)
  tsoil=anx(:,:,2)
!
!-----------------------------------------------------------------------
!
!  Some initializations
!
!-----------------------------------------------------------------------

  latnot(1)=trulat1
  latnot(2)=trulat2
  CALL setmapr(mapproj,sclfct,latnot,trulon)
!
  CALL lltoxy(1,1,ctrlat,ctrlon,x0,y0)
  DO i=1,nx
    xgrid(i)=(FLOAT(i)-(FLOAT(nx-3)/2.)-1.5)*dx + x0
  END DO
  DO j=1,ny
    ygrid(j)=(FLOAT(j)-(FLOAT(ny-3)/2.)-1.5)*dy + y0
  END DO
!
  CALL xytoll(1,1, xgrid(1), ygrid(1),lllat,lllon)
  CALL xytoll(1,1,xgrid(nx), ygrid(1),lrlat,lrlon)
  CALL xytoll(1,1,xgrid(nx),ygrid(ny),urlat,urlon)
  CALL xytoll(1,1, xgrid(1),ygrid(ny),ullat,ullon)
  PRINT *, ' lower  left lat,lon = ',lllat,lllon
  PRINT *, ' lower right lat,lon = ',lrlat,lrlon
  PRINT *, ' upper right lat,lon = ',urlat,urlon
  PRINT *, ' upper  left lat,lon = ',ullat,ullon
!
  defrange=6400.e06     ! meters**2
  DO j=1,ny
    DO i=1,nx
      range(i,j)=defrange
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Use initime in the ARPS input file as the time.
!
!-----------------------------------------------------------------------
!
  CALL julday( year, month, day, jday)
  myr=MOD(year,100)
  WRITE(c9time,'(i2.2,i3.3,i2.2,i2.2)') myr,jday,hour,minute
!
!-----------------------------------------------------------------------
!
!  Read in observations. match with station list
!
!-----------------------------------------------------------------------
!
  CALL read_obs(c9time,                                                 &
                xtnam,xtlat,xtlon,xtelev,xts,k_sta)
!
  CALL lltoxy(k_sta,1,xtlat,xtlon,xtx,xty)
!
  DO j= 1,nvar_ob
    DO i= 1,k_sta
      isrc(i)   =1
      obs(i,j)  =xts(i,j)
      qobs(i,j) =qsrc(isrc(i),j)
      rely(i,j,1)   =10
    END DO
  END DO

  CALL barnes_cntrl(nx,ny,mxstn,2,maxrng,                               &
           k_sta,xgrid,ygrid,                                           &
           obs,odiff,obanx,xtx,xty,isrc,qobs,rely,                      &
           var_anx,npass,kappa,irngsel,range,rngobs,                    &
           wlim,klim,knt,rngwgt,wgtsum,zsum,istat,                      &
           anx,istatus)

!
!-----------------------------------------------------------------------
!
!  Calculate again for solar radiation  defrange=1600.
!
!-----------------------------------------------------------------------
!
  defrange=1600.e06    ! meters**2
  DO j=1,ny
    DO i=1,nx
      range(i,j)=defrange
    END DO
  END DO

  CALL barnes_cntrl(nx,ny,mxstn,1,maxrng,                               &
           k_sta,xgrid,ygrid,                                           &
           obs(1,3),odiff(1,3),obanx(1,3),xtx,xty,                      &
           isrc,qobs(1,3),rely(1,3,1),                                  &
           var_anx(3),npass,kappa,irngsel(3),range,rngobs,              &
           wlim,klim,knt(3),rngwgt,wgtsum(3),zsum(3),istat,             &
           anx(1,1,3),istatus)

!
!-----------------------------------------------------------------------
!
!  Make adjustments to output data.
!    1) Convert to Kelvin
!    2) Set wetness values to constants from the input file.
!
!-----------------------------------------------------------------------
!
  DO j=1,ny
    DO i=1,nx
      tsfc(i,j)=tsfc(i,j)+273.15
      tsoil(i,j)=tsoil(i,j)+273.15
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Open a file.
!  Write the var_anx array to the file which was
!  specified as soilinfl in the ARPS input file.
!
!-----------------------------------------------------------------------
!
  idummy=0
  rdummy=0.
!
  tsfcout=1
  tsoilout=1
  wsfcout=0
  wdpout=0
  wcanpout=0
!
  CALL getunit( sfcunit )

  CALL asnctl ('NEWLOCAL', 1, ierr)
  CALL asnfile(soilinfl, '-F f77 -N ieee', ierr)

  OPEN(UNIT=sfcunit,FILE=soilinfl,FORM='unformatted',                   &
         STATUS='unknown',IOSTAT=istat)

  WRITE (sfcunit) nx,ny
  WRITE (sfcunit) tsfcout,tsoilout,wsfcout,wdpout,wcanpout,             &
                  idummy,idummy,idummy,idummy,idummy,                   &
                  idummy,idummy,idummy,idummy,idummy,                   &
                  idummy,idummy,idummy,idummy,idummy

  WRITE (sfcunit) dx,dy,rdummy,rdummy,rdummy,                           &
                  rdummy,rdummy,rdummy,rdummy,rdummy,                   &
                  rdummy,rdummy,rdummy,rdummy,rdummy,                   &
                  rdummy,rdummy,rdummy,rdummy,rdummy

  WRITE (sfcunit) tsfc
  WRITE (sfcunit) tsoil
  CLOSE (sfcunit)
!
!-----------------------------------------------------------------------
!
!  Grads output
!
!-----------------------------------------------------------------------
!
  gradsfile=c9time(1:9)//'.grads'
  gradsctl=c9time(1:9)//'.ctl'
  OPEN(11,FILE=gradsfile,FORM='unformatted')
  DO k = 1,nvar_ob
    WRITE(11) ((anx(i,j,k),i=1,nx),j=1,ny)
  END DO
  CLOSE(11)
  OPEN(11,FILE=gradsctl,STATUS='unknown')
  WRITE(11,'(a,a)') 'TITLE   OLAPS soil temperatures and ',             &
                    'surface solar radiation flux'
  WRITE(11,'(a,a)') 'DSET    ',gradsfile
  WRITE(11,'(a)')   'UNDEF   -999.0'
  WRITE(11,'(a)')   'OPTIONS sequential'
  WRITE(11,'(a,i4,a,2e16.8)')                                           &
                    'XDEF    ',nx,' LINEAR ',-0.5*dx/1000.,dx/1000.
  WRITE(11,'(a,i4,a,2e16.8)')                                           &
                    'YDEF    ',ny,' LINEAR ',-0.5*dy/1000.,dy/1000.
  WRITE(11,'(a)')   'ZDEF    1    LINEAR  0   1'
  WRITE(11,'(a,i2.2,a,i2.2,a,i2.2,a3,i4.4,1x,a)')                       &
      'TDEF    1    LINEAR  ',hour,':',minute,'Z',                      &
      day,monnam(month),year,'1HR'
  WRITE(11,'(a)')   'VARS     3'
  WRITE(11,'(a)')   'ts05    0  99  Soil temperature (K) in 5 cm'
  WRITE(11,'(a)')   'ts30    0  99  Soil temperature (K) in 30 cm'
  WRITE(11,'(a)')   'srad    0  99  Incoming solar radiation'
  WRITE(11,'(a)')   'ENDVARS'
  CLOSE(11)
END PROGRAM arpssoil
