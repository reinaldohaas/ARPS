!
!##################################################################
!##################################################################
!######                                                      ######
!######                  PROGRAM ARPSINTRP_LS                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

PROGRAM arpsintrp_ls
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  This program interpolates gridded data from one ARPS grid to another.
!  It can be used to prepare data for running ARPS in a one-way nested
!  mode. It's exepcted to replace ARPSR2H in this capacity.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!
!  3/27/1997. Written based on ARPSR2H and ARPSCVT.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!  nx, ny, nz: Dimensions of input grid.
!  nx1, ny1, nz1: Dimensions of output grid.
!-----------------------------------------------------------------------
!
  INTEGER :: nx       ! Number of input grid points in the x-direction
  INTEGER :: ny       ! Number of input grid points in the y-direction
  INTEGER :: nz       ! Number of input grid points in the z-direction
  INTEGER :: nzsoil   ! Number of input grid points in the z-direction

  INTEGER :: nx1      ! Number of output grid points in the x-direction
  INTEGER :: ny1      ! Number of output grid points in the y-direction
  INTEGER :: nz1      ! Number of output grid points in the z-direction
  INTEGER :: nzsoil1  ! Number of output grid points in the z-direction

  INTEGER :: nxyz,nxy,nxyz1,nxy1
  INTEGER :: nxyzsoil, nxyzsoil1

!
!-----------------------------------------------------------------------
!
!  Parameters for the utput grid.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nstyps
  PARAMETER (nstyps=1)

  INTEGER :: samgrd   ! Are the output and the input grids the same?
                      ! =1, the grids are the same
                      ! =0, the grids are different

  REAL :: xorig1, yorig1, zorig1
  REAL :: xctr1 , yctr1

  REAL :: dx1,dy1    ! Grid intervals of the refined grid.

  INTEGER :: strhopt1 ! Vertical grid stretching option.
                      ! = 0, no stretching in vertical.
                      ! >= 1, with stretching in vertical.
  REAL :: dz1        ! Average grid spacing in vertical direction in
                     ! transformed computational space (m).
  REAL :: dzmin1     ! Minimun grid spacing in vertical direction in
                     ! physcal space (m).

  REAL :: zrefsfc1   ! The reference height of the surface
                     ! (ground level) (m)

  REAL :: dlayer11   ! The depth of the lower layer with uniform
                     ! (dz=dzmin) vertical spacing (m)

  REAL :: dlayer21   ! The depth of the mid layer with stetched
                     ! vertical spacing (m)

  REAL :: strhtune1  ! Tuning parameter for stretching option 2
                     ! A Value between 0.2 and 5.0 is appropriate.
                     ! A larger value gives a more linear stretching.

  REAL :: zflat1     ! The height at which the grid levels
                     ! becomes flat in the terrain-following
                     ! coordinate transformation (m).
!
!-----------------------------------------------------------------------
!
!  Note:
!
!  Given nx1, ny1 and nz1, the physical doma3din size of the refined
!  grid will be xl1=(nx1-3)*dx1 by yl1=(ny1-3)*dy1 by zh1=(nz1-3)*dz1.
!  Dx1, dy1 and dz1 are the grid intervals of the refined grid.
!
!-----------------------------------------------------------------------
!

  REAL, ALLOCATABLE :: tsoil_tmp(:,:,:,:)  ! 4-D real    array read in
  REAL, ALLOCATABLE :: qsoil_tmp(:,:,:,:)  ! 4-D real    array read in

  REAL, ALLOCATABLE :: a4din(:,:,:,:)  ! 4-D real    array read in
  REAL, ALLOCATABLE :: tsoilin(:,:,:,:)  ! 4-D real    array read in
  REAL, ALLOCATABLE :: qsoilin(:,:,:,:)  ! 4-D real    array read in
  REAL, ALLOCATABLE :: wetcanpin(:,:,:)  ! 3-D real    array read in

  REAL, ALLOCATABLE :: a3dsoilin(:,:,:)  ! 3-D real    array read in
  REAL, ALLOCATABLE :: a3din(:,:,:)  ! 3-D real    array read in

  INTEGER, ALLOCATABLE :: i2din(:,:) ! 2-D integer array read in
  REAL, ALLOCATABLE :: a2din(:,:)    ! 2-D real    array read in
  REAL, ALLOCATABLE :: x     (:)     ! The x-coord. of the physical and
                                     ! computational grid. Defined at u-point.
  REAL, ALLOCATABLE :: y     (:)     ! The y-coord. of the physical and
                                     ! computational grid. Defined at v-point.
  REAL, ALLOCATABLE :: z     (:)     ! The z-coord. of the computational grid.
                                     ! Defined at w-point on the staggered grid.
  REAL, ALLOCATABLE :: zp    (:,:,:) ! The physical height coordinate defined at
                                     ! w-point on the staggered grid.
  REAL, ALLOCATABLE :: zpsoil(:,:,:) ! The physical height coordinate defined at
                                     ! w-point on the staggered grid.
                                     ! w-point on the staggered grid.
  REAL, ALLOCATABLE :: hterain(:,:)  ! The height of terrain.

  REAL, ALLOCATABLE :: tem3d(:,:,:)  ! Work array
!
!-----------------------------------------------------------------------
!
!  Arrays on the new grid:
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: a4dout(:,:,:,:)! 4-D real array interpolated from
                                      ! a4din.
  REAL, ALLOCATABLE :: tsoilout(:,:,:,:)
  REAL, ALLOCATABLE :: qsoilout(:,:,:,:)
  REAL, ALLOCATABLE :: wetcanpout(:,:,:)  ! 3-D real    array read in
  REAL, ALLOCATABLE :: a3dsoilout(:,:,:) ! 3-D real array interpolated
                                         ! from a3dsoilin.
  REAL, ALLOCATABLE :: a3dout(:,:,:) ! 3-D real array interpolated from a3din.

  INTEGER, ALLOCATABLE :: i2dout(:,:)! 2-D integer array derived from i2din
  REAL, ALLOCATABLE :: a2dout(:,:)   ! 2-D real array interpolated from a2din

  REAL, ALLOCATABLE :: x1    (:)     ! New grid x-coord. on the original grid
  REAL, ALLOCATABLE :: x11   (:)     ! New grid x-coord. set for new grid

  REAL, ALLOCATABLE :: y1    (:)     ! New grid y-coord. on the original grid
  REAL, ALLOCATABLE :: y11   (:)     ! New grid y-coord. set for new grid

  REAL, ALLOCATABLE :: z1    (:)     ! The z-coord. of the computational grid.
                                     ! Defined at w-point on the staggered grid.
  REAL, ALLOCATABLE :: zp1   (:,:,:) ! The physical height coordinate defined at
                                     ! w-point on the staggered grid.
  REAL, ALLOCATABLE :: zpsoil1(:,:,:)! The physical height coordinate defined at
                                     ! w-point on the staggered grid.
  REAL, ALLOCATABLE :: j3soil1(:,:,:)! Jacobian
  REAL, ALLOCATABLE :: j3soilinv1(:,:,:)! Jacobian inverse
  REAL, ALLOCATABLE :: hterain1(:,:) ! Terrain height (m)

  REAL, ALLOCATABLE :: tem3d1(:,:,:) ! Work array

  REAL, ALLOCATABLE :: zp1d1 (:)     ! Temporary array
  REAL, ALLOCATABLE :: dzp1d1(:)     ! Temporary array

  REAL :: zflat11,za,zb

!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k
  INTEGER :: is, js
  REAL :: s1,s2,s3,s4
  REAL :: amin, amax
  REAL :: xs1,ys1

  CHARACTER (LEN=256) :: basdmpfn
  INTEGER :: lbasdmpf
  CHARACTER (LEN=256) :: ternfn,sfcoutfl,soiloutfl,temchar
  INTEGER :: lternfn,lfn
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'bndry.inc'
  INCLUDE 'indtflg.inc'

  INTEGER :: hinfmt,nhisfile_max,nhisfile,lengbf,nf,lenfil
  PARAMETER (nhisfile_max=200)
  CHARACTER (LEN=256) :: grdbasfn,hisfile(nhisfile_max)
  INTEGER :: ireturn

  INTEGER :: houtfmt
  CHARACTER (LEN=256) :: filename

  INTEGER :: grdbas

  REAL :: time

  INTEGER :: nfilemax
  PARAMETER (nfilemax=200)

  CHARACTER (LEN=256) :: exbcfn

  CHARACTER (LEN=80) :: timsnd
  CHARACTER (LEN=80) :: new_runname
  INTEGER :: tmstrln
  CHARACTER (LEN=15) :: ctime

  INTEGER :: nfile, length
  REAL :: xeps, yeps
  REAL :: ctrx,ctry,swx,swy,alatpro(2),sclf,dxscl,dyscl
  REAL :: ctrlat1,ctrlon1,latitud1

  CHARACTER (LEN=40) :: fmtver
  PARAMETER (fmtver='004.10 Binary Data')
  CHARACTER (LEN=40) :: fmtverin

  CHARACTER (LEN=10) :: tmunit
  CHARACTER (LEN=12) :: label
!  INTEGER :: nxin,nyin,nzin
  INTEGER :: nxin,nyin,nzin,nzsoilin
  INTEGER :: stgr,oldver
  INTEGER :: inch,nchanl,exbchanl,sfchanl,soilchanl,trnchanl
  INTEGER :: iyr, imon, idy, ihr, imin, isec
  INTEGER :: ierr, itema
  INTEGER :: varout1,mstout1,rainout1,prcout1,iceout1,tkeout1,trbout1
  INTEGER :: sfcout1,landout1

  REAL :: xgrdorg1,ygrdorg1

  INTEGER :: idummy
  REAL :: rdummy

  INTEGER :: istat
  LOGICAL :: fexist, dmpexbc, lsfcdmp, lsoildmp
!
!-----------------------------------------------------------------------
!
!  namelist Declarations:
!
!-----------------------------------------------------------------------
!
  NAMELIST /INPUT/ hinfmt,nhisfile, grdbasfn, hisfile
  NAMELIST /output_dims/ nx1,ny1,nz1

  NAMELIST /jobname/ runname

  NAMELIST /newgrid/ samgrd,strhopt1,xctr1,yctr1,                       &
            dx1,dy1,dz1,dzmin1,                                         &
            zrefsfc1,dlayer11,dlayer21,strhtune1,zflat1

  NAMELIST /output/ dirname,exbcdmp,hdmpfmt,grbpkbit,                   &
            grdout,basout,varout,mstout,rainout,prcout,iceout,          &
            tkeout, trbout,sfcout,landout,                              &
            qcexout,qrexout,qiexout,qsexout,qhexout,                    &
            totout,filcmprs

  REAL :: dzsoil1,zrefsoil1
!  INTEGER :: soilstrhopt
!  REAL :: soildzmin,soildlayer1,soildlayer2,soilstrhtune

  NAMELIST /newgrid_soil/ nzsoil1,dzsoil1,zrefsoil1,soilstrhopt,        &
                          soildzmin,soildlayer1,soildlayer2,            &
                          soilstrhtune

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  mgrid = 1
  nestgrd = 0

  WRITE(6,'(/9(/5x,a)/)')                                               &
     '###############################################################', &
     '###############################################################', &
     '###                                                         ###', &
     '###                Welcome to ARPSINTRP                     ###', &
     '###      This program converts the history dump data        ###', &
     '###      sets generated by ARPS, between various formats.   ###', &
     '###                                                         ###', &
     '###############################################################', &
     '###############################################################'

!
!-----------------------------------------------------------------------
!  Get the names of the input data files.
!-----------------------------------------------------------------------
!
  CALL get_input_file_names(5,hinfmt,grdbasfn,hisfile,nhisfile)

  lengbf = len_trim(grdbasfn)

!  CALL get_dims_from_data(hinfmt,grdbasfn(1:lengbf),                    &
!       nx,ny,nz,nstyps, ireturn)
  CALL get_dims_from_data(hinfmt,grdbasfn(1:lengbf),                    &
       nx,ny,nz,nzsoil,nstyps, ireturn)

  IF( ireturn /= 0 ) THEN
    PRINT*,'Problem occured when trying to get dimensions from data.'
    PRINT*,'Program stopped.'
    STOP
  END IF

!  WRITE(6,'(3(a,i5))') 'nx =',nx,', ny=',ny,', nz=',nz
  WRITE(6,'(4(a,i5))') 'nx =',nx,', ny=',ny,', nz=',nz,', nzsoil=',nzsoil

  IF( nhisfile > nhisfile_max) THEN
    WRITE(6,'(1x,a,i3,/1x,a,/1x,a,/1x,a)')                              &
        'The number of history files to be processed exceeded ',        &
        nhisfile_max,' please reduce the number of files to be processed',&
        'in a single job, or edit the program and re-set parameter ',   &
        'nhisfile_max to a larger value. '
    STOP 9010
  END IF

  length = LEN_trim( grdbasfn )
  WRITE(6,'(1x,a,a)')  ' grdbasfn    =', grdbasfn(1:length)

  DO i=1,nhisfile
    length = LEN_trim( hisfile(i) )
    WRITE(6,'(1x,a,i3,a,a)') ' hisfile(',i,')=',hisfile(i)(1:length)
  END DO

!-----------------------------------------------------------------------
! Read output grid dimensions
!-----------------------------------------------------------------------

  READ(5,output_dims, END=105)
  WRITE(6,'(/a,a/)')                                                    &
      'NAMELIST block output_dims successfully read.'

  READ(5,newgrid_soil, END=105)
  WRITE(6,'(/a,a/)')                                                    &
      'NAMELIST block newgrid_soil successfully read.'

!  WRITE(6,'(3(a,i5))') 'nx1=',nx1,', ny1=',ny1,', nz1=',nz1
  WRITE(6,'(4(a,i5))') 'nx1=',nx1,', ny1=',ny1,', nz1=',nz1, &
                       ', nzsoil1=',nzsoil1

  GO TO 10

105   WRITE(6,'(/a,a/)')                                                &
          'Error reading NAMELIST block intrp_dims. ',                  &
          'Program ARPSINTRP stopped.'
  STOP 1

10 CONTINUE

  READ (5,jobname,ERR=9000)

  WRITE(6,'(/5x,a,a)') 'The name of this run is: ', runname
  new_runname = runname
  CALL gtlfnkey(new_runname, lfnkey)
!
!-----------------------------------------------------------------------
!
!  Set the output grid and the variable control parameters
!
!-----------------------------------------------------------------------
!
  READ (5,newgrid,ERR=9000)

  PRINT*
  PRINT*,' Input parameters for the new refined grid:'
  PRINT*

  PRINT*,' Input samgrd:'
  PRINT*,' Input was ',samgrd

  PRINT*,' Input dx1:'
  PRINT*,' Input was ',dx1

  PRINT*,' Input dy1:'
  PRINT*,' Input was ',dy1

  PRINT*,' Input strhopt1:'
  PRINT*,' Input was ',strhopt1

  PRINT*,' Input dz1:'
  PRINT*,' Input was ',dz1

  PRINT*,' Input dzmin1:'
  PRINT*,' Input was ',dzmin1

  PRINT*,' Input xctr1:'
  PRINT*,' Input was ',xctr1

  PRINT*,' Input yctr1:'
  PRINT*,' Input was ',yctr1
!
!-----------------------------------------------------------------------
!
!  Set the control parameters for output:
!
!-----------------------------------------------------------------------
!
  WRITE(6,'(/a/)')                                                      &
      ' Reading in control parameters for the output data files..'

  READ (5,output,ERR=9000)

  houtfmt = hdmpfmt

  IF( houtfmt /= 1 ) THEN
    WRITE(6,'(/1x,a,a/)') 'Output format is not 1. Reset it to 1.'
    houtfmt = 1
  END IF

  totout = 1
  basout = 0

  ldirnam=LEN_trim(dirname)

  lengbf=LEN_trim(grdbasfn)
  WRITE(6,'(/a,a)')' The grid/base name is ', grdbasfn(1:lengbf)

  INQUIRE(FILE=grdbasfn(1:lengbf), EXIST = fexist )
  IF( fexist ) GO TO 110

  INQUIRE(FILE=grdbasfn(1:lengbf)//'.Z', EXIST = fexist )
  IF( fexist ) THEN
    CALL uncmprs( grdbasfn(1:lengbf)//'.Z' )
    GO TO 110
  END IF

  INQUIRE(FILE=grdbasfn(1:lengbf)//'.gz', EXIST = fexist )
  IF( fexist ) THEN
    CALL uncmprs( grdbasfn(1:lengbf)//'.gz' )
    GO TO 110
  END IF

  WRITE(6,'(/1x,a,/1x,a/)')                                             &
      'File '//grdbasfn(1:lengbf)//                                     &
      ' or its compressed version not found.',                          &
      'Program stopped in ARPSINTRP.'
  STOP 9011

  110   CONTINUE


!  ALLOCATE(a4din(nx,ny,nzsoil,0:nstyps))
!  a4din = 0.0
  ALLOCATE(tsoilin(nx,ny,nzsoil,0:nstyps))
  tsoilin = 0.0
  ALLOCATE(qsoilin(nx,ny,nzsoil,0:nstyps))
  qsoilin = 0.0
  ALLOCATE(wetcanpin(nx,ny,nzsoil))
  wetcanpin = 0.0
  ALLOCATE(a3dsoilin(nx,ny,nzsoil))
  a3dsoilin = 0.0
  ALLOCATE(a3din(nx,ny,nz))
  a3din = 0.0
  ALLOCATE(i2din(nx,ny))
  i2din = 0
  ALLOCATE(a2din(nx,ny))
  a2din =0.0
  ALLOCATE(x     (nx))
  ALLOCATE(y     (ny))
  ALLOCATE(z     (nz))
  x = 0.0
  y = 0.0
  z = 0.0
  ALLOCATE(zp    (nx,ny,nz))
  zp = 0.0
  ALLOCATE(hterain(nx,ny))
  hterain=0.0
  ALLOCATE(tem3d(nx,ny,nz))
  tem3d=0.0
!

  allocate(tsoil_tmp(nx,ny,nzsoil1,0:nstyps))
  tsoil_tmp = 0.0
  allocate(qsoil_tmp(nx,ny,nzsoil1,0:nstyps))
  qsoil_tmp = 0.0

  ALLOCATE(tsoilout(nx1,ny1,nzsoil1,0:nstyps))
  tsoilout = 0.0
  ALLOCATE(qsoilout(nx1,ny1,nzsoil1,0:nstyps))
  qsoilout = 0.0
  ALLOCATE(wetcanpout(nx1,ny1,nzsoil1))
  wetcanpout = 0.0
  ALLOCATE(a3dout(nx1,ny1,nz1))
  a3dout=0.0
  ALLOCATE(i2dout(nx1,ny1))
  i2dout=0.0
  ALLOCATE(a2dout(nx1,ny1))
  a2dout=0.0
  ALLOCATE(x1    (nx1))
  x1=0.0
  ALLOCATE(x11   (nx1))
  x11=0.0
  ALLOCATE(y1    (ny1))
  y1=0.0
  ALLOCATE(y11   (ny1))
  y11=0.0
  ALLOCATE(z1    (nz1))
  z1=0.0
  ALLOCATE(zp1   (nx1,ny1,nz1))
  zp1=0.0
  ALLOCATE(zpsoil1   (nx1,ny1,nzsoil1))
  zpsoil1=0.0
  ALLOCATE(hterain1(nx1,ny1))
  hterain1=0.0
  ALLOCATE(tem3d1(nx1,ny1,nz1))
  tem3d1=0.0
  ALLOCATE(zp1d1 (nz1))
  zp1d1=0.0
  ALLOCATE(dzp1d1(nz1))
  dzp1d1=0.0

!
!-----------------------------------------------------------------------
!
!  Get the IO unit numbers for input grid and base state file
!
!-----------------------------------------------------------------------
!
  CALL getunit( inch )
!
!-----------------------------------------------------------------------
!
!  Cray routines to force binary data file to be in the IEEE format
!
!-----------------------------------------------------------------------
!
  CALL asnctl ('NEWLOCAL', 1, ierr)
  CALL asnfile(grdbasfn(1:lengbf), '-F f77 -N ieee', ierr)
!
!-----------------------------------------------------------------------
!
!  Open the input grdbas file
!
!-----------------------------------------------------------------------
!
  grdbas = 1

  OPEN(UNIT=inch,FILE=grdbasfn(1:lengbf),                               &
           STATUS='old',FORM='unformatted',IOSTAT=istat)
  IF( istat /= 0 ) GO TO 9001

  READ(inch,ERR=9110,END=9120) fmtverin
  READ(inch,ERR=9110,END=9120) runname
  READ(inch,ERR=9110,END=9120) nocmnt

  IF( nocmnt > 0 ) THEN
    DO i=1,nocmnt
      READ(inch,ERR=9110,END=9120) cmnt(i)
    END DO
  END IF

  READ(inch,ERR=9110,END=9120) time,tmunit

  curtim = time
!
!-----------------------------------------------------------------------
!
!  Get dimensions of data in binary file and check against
!  the dimensions defined in dims.inc
!
!-----------------------------------------------------------------------
!
  READ(inch,ERR=9110,END=9120) nxin, nyin, nzin, nzsoilin

  IF( nxin /= nx .OR. nyin /= ny .OR. nzin /= nz .OR. &
      nzsoilin /= nzsoil) THEN
    WRITE(6,'(1x,a)')                                                   &
         ' Dimensions in ARPSINTRP inconsistent with data.'
    WRITE(6,'(1x,a,4I15)') ' Read were: ', nxin, nyin, nzin, nzsoilin
    WRITE(6,'(1x,a)')                                                   &
         ' Program aborted in ARPSINTRP.'
    STOP
  END IF

  WRITE(6,'(1x,a,f8.1,a,f8.3,a/)')                                      &
         'To read grid and base state data at time ', time,             &
         ' secs = ',(time/60.),' mins.'

  READ(inch,ERR=9110,END=9120)                                          &
         grdin,basin,varin,mstin,icein,                                 &
         trbin,idummy,idummy,landin,totin,                              &
         tkein,idummy,idummy,mapproj,month,                             &
         day, year, hour, minute, second

  IF ( grdin /= 1 .OR. basin /= 1 ) THEN
    WRITE(6,'(1x,a,a,a/a)')                                             &
        'File '//grdbasfn(1:lengbf)//' is not .bingrdbas file',         &
        'A .bingrdbas file is required. Program stoped in ARPSINTRP'
    STOP 9012
  END IF

  IF ( varin == 1 ) THEN
    varout1  = varout
  ELSE
    varout1  = 0
  END IF

  IF ( mstin == 1 ) THEN
    mstout1  = mstout
  ELSE
    mstout1 = 0
  END IF

  IF ( icein == 1 ) THEN
    iceout1  = iceout
  ELSE
    iceout1  = 0
  END IF

  IF ( tkein == 1 ) THEN
    tkeout1  = tkeout
  ELSE
    tkeout1  = 0
  END IF

  IF ( trbin == 1 ) THEN
    trbout1  = trbout
  ELSE
    trbout1  = 0
  END IF

  IF ( sfcin == 1 ) THEN
    sfcout1  = sfcout
  ELSE
    sfcout1  = 0
  END IF

  IF ( landin == 1 ) THEN
    landout1 = landout
  ELSE
    landout1 = 0
  END IF

  READ(inch,ERR=9110,END=9120)                                          &
         umove,vmove,xgrdorg,ygrdorg,trulat1,                           &
         trulat2,trulon,sclfct,rdummy,rdummy,                           &
         rdummy,rdummy,rdummy,rdummy,rdummy,                            &
         tstop,thisdmp,latitud,ctrlat,ctrlon

  IF ( totin /= 0 ) THEN
!
!-----------------------------------------------------------------------
!
!  Read in additional parameters for ARPS history dump 4.0 or later
!  version.
!
!-----------------------------------------------------------------------
!
    READ(inch,ERR=9110,END=9120)                                        &
         nstyp, prcout,idummy,idummy,idummy,                            &
         idummy,idummy,idummy,idummy,idummy,                            &
         idummy,idummy,idummy,idummy,idummy,                            &
         idummy,idummy,idummy,idummy,idummy

    IF ( nstyp < 1 ) THEN
      nstyp = 1
    END IF

    IF ( prcin == 1 ) THEN
      prcout1 = prcout
    ELSE
      prcout1 = 0
    END IF

    READ(inch,ERR=9110,END=9120)                                        &
         rdummy,rdummy,rdummy,rdummy,rdummy,                            &
         rdummy,rdummy,rdummy,rdummy,rdummy,                            &
         rdummy,rdummy,rdummy,rdummy,rdummy,                            &
         rdummy,rdummy,rdummy,rdummy,rdummy

  END IF
!
!-----------------------------------------------------------------------
!
!  Read in x,y and z at grid cell centers (scalar points).
!
!----------------------------------------------------------------------
!
  IF( grdin == 1 .OR. grdbas == 1 ) THEN
    READ(inch,ERR=9110,END=9120) label
    READ(inch,ERR=9110,END=9120) x
    WRITE(6,'(1x,a,a12,a)')                                             &
        'Field ',label,' was read into array x.'

    READ(inch,ERR=9110,END=9120) label
    READ(inch,ERR=9110,END=9120) y
    WRITE(6,'(1x,a,a12,a)')                                             &
        'Field ',label,' was read into array y.'

    READ(inch,ERR=9110,END=9120) label
    READ(inch,ERR=9110,END=9120) z
    WRITE(6,'(1x,a,a12,a)')                                             &
        'Field ',label,' was read into array z.'

    READ(inch,ERR=9110,END=9120) label
    READ(inch,ERR=9110,END=9120) zp
    WRITE(6,'(1x,a,a12,a)')                                             &
        'Field ',label,' was read into array zp.'

    READ(inch,ERR=9110,END=9120) label
    READ(inch,ERR=9110,END=9120) zpsoil
    WRITE(6,'(1x,a,a12,a)')                                             &
        'Field ',label,' was read into array zpsoil.'
  END IF  ! grdin
!
!-----------------------------------------------------------------------
!
!  Set up the map projection of the input grid.
!
!-----------------------------------------------------------------------
!
  alatpro(1)=trulat1
  alatpro(2)=trulat2

  dx = x(2)-x(1)
  dy = y(2)-y(1)
  IF( sclfct /= 1.0) THEN
    sclf  = 1.0/sclfct
    dxscl = dx*sclf
    dyscl = dy*sclf
  ELSE
    sclf  = 1.0
    dxscl = dx
    dyscl = dy
  END IF

  PRINT*,mapproj,sclf,alatpro,trulon,ctrlat,ctrlon

  CALL setmapr( mapproj,sclf,alatpro,trulon )
  CALL lltoxy( 1,1, ctrlat,ctrlon, ctrx, ctry )
  swx = ctrx - (FLOAT(nx-3)/2.)*dxscl
  swy = ctry - (FLOAT(ny-3)/2.)*dyscl
  CALL setorig( 1, swx, swy) ! set up the model origin to the coord.

  DO j=1,ny
    DO i=1,nx
      hterain(i,j) = zp(i,j,2)
    END DO
  END DO

  dx = x(2) - x(1)
  dy = y(2) - y(1)
  dz = z(2) - z(1)

  xorig = x(2)
  yorig = y(2)
  zorig = z(2)

  WRITE(6,'(1x,a,3f15.3)') 'xorig, yorig, zorig =',                     &
       xorig, yorig, zorig
!
!-----------------------------------------------------------------------
!
!  If the new grid is same as the original one, no need to
!  interpolate
!
!-----------------------------------------------------------------------
!
  IF ( samgrd == 1 ) THEN
    DO i=1,nx1
      x1(i)=x(i)
      x11(i)=x(i)
    END DO

    DO j=1,ny1
      y1(j)=y(j)
      y11(j)=y(j)
    END DO

    DO k=1,nz1
      z1(k)=z(k)
    END DO

    DO j=1,ny1
      DO i=1,nx1
        hterain1(i,j) = hterain(i,j)
      END DO
    END DO

    DO k=1,nz1
      DO j=1,ny1
        DO i=1,nx1
          zp1(i,j,k) = zp(i,j,k)
        END DO
      END DO
    END DO

    DO k=1,nzsoil1
      DO j=1,ny1
        DO i=1,nx1
          zpsoil1(i,j,k) = zpsoil(i,j,k)
        END DO
      END DO
    END DO

    ctrlat1  = ctrlat
    ctrlon1  = ctrlon
    latitud1 = latitud

    xgrdorg1 = xgrdorg
    ygrdorg1 = ygrdorg

  ELSE

    ebc = 3
    wbc = 3
    sbc = 3
    nbc = 3
!
!-----------------------------------------------------------------------
!
!  If output grid differs from the input grid, perform spatial
!  interpolations.
!
!-----------------------------------------------------------------------
!
    xorig1 = xctr1 - (nx1-3)*dx1*0.5
    yorig1 = yctr1 - (ny1-3)*dy1*0.5
    zorig1 = zorig

    DO i=1,nx1
      x1(i)=xorig1+(i-2)*dx1
    END DO

    DO j=1,ny1
      y1(j)=xorig1+(j-2)*dy1
    END DO

    DO k=1,nz1
      z1(k)=zorig1+(k-2)*dz1
    END DO

    xeps = 0.01*dx
    yeps = 0.01*dy

    IF(x1(1) < x(1)-xeps.OR.x1(nx1) > x(nx)+xeps.OR.                    &
          y1(1) < y(1)-yeps.OR.y1(ny1) > y(ny)+yeps) THEN
      WRITE(6,'(3(/5x,a),/5x,2(a,f12.4),2(a,i6),2(a,f12.4)/5x,a)')      &
          'Sorry, at least part of your new grid is outside the border of', &
          'the original grid, please check input parameters',           &
          'dx1,dy1, nx1, ny1, xctr1 and yctr1. Currently,',             &
          'dx1=',dx1,', dy1=',dy1,', nx1=',nx1,', ny1=',ny1,            &
          ', xctr1=',xctr1,', yctr1=',yctr1,                            &
          'Job stopped in ARPSINTRP.'
      WRITE(6,'(1x,4(a,f12.4)/1x,4(a,f12.4))')                          &
          'x1(1) =',x1(1),' x1(nx1) =',x1(nx1),                         &
          'y1(1) =',y1(1),' y1(ny1) =',y1(ny1),                         &
          'x (1) =',x (1),' x  (nx) =',x  (nx),                         &
          'y (1) =',y (1),' y  (ny) =',y  (ny)
      STOP 9013
    END IF

    CALL xytoll(1,1,xctr1,yctr1,ctrlat1,ctrlon1)
    PRINT*,'ctrlat1=',ctrlat1,', ctrlon1=',ctrlon1

    latitud1 = ctrlat1

    CALL xytoll(1,1,xorig1,yorig1,swx,swy)  ! Find the lat/lon of
                                            ! the new grid origin
    CALL setorig( 2, swx,swy )              ! Set the new origin

    xgrdorg1 = 0.0
    ygrdorg1 = 0.0

    DO i=1,nx1
      x11(i)=x1(i)-xorig1
    END DO

    DO j=1,ny1
      y11(j)=y1(j)-yorig1
    END DO

!
!-----------------------------------------------------------------------
!
!  Intepolate terrain to the new grid
!
!-----------------------------------------------------------------------
!

    DO i=1,nx1-1
      DO j=1,ny1-1

        xs1= (x1(i)+x1(i+1))*0.5
        ys1= (y1(j)+y1(j+1))*0.5

        is = MAX(1, MIN(nx-2, INT((xs1-(x(1)+x(2))*0.5)/dx)+1 ))
        js = MAX(1, MIN(ny-2, INT((ys1-(y(1)+y(2))*0.5)/dy)+1 ))

        s1=ABS((xs1-(x(is  )+x(is+1))*0.5)*(ys1-(y(js  )+y(js+1))*0.5))
        s2=ABS((xs1-(x(is+1)+x(is+2))*0.5)*(ys1-(y(js  )+y(js+1))*0.5))
        s3=ABS((xs1-(x(is+1)+x(is+2))*0.5)*(ys1-(y(js+1)+y(js+2))*0.5))
        s4=ABS((xs1-(x(is  )+x(is+1))*0.5)*(ys1-(y(js+1)+y(js+2))*0.5))

        hterain1(i,j) =                                                 &
            (hterain(is  ,js  )*s3+hterain(is+1,js  )*s4                &
            +hterain(is+1,js+1)*s1+hterain(is  ,js+1)*s2)               &
            /(s1+s2+s3+s4)

      END DO
    END DO

!
!-----------------------------------------------------------------------
!
!  Set up a stretched vertical grid for the output grid.
!
!  For strhopt1=1, function y = a+b*x**3 is used to specify dz as a
!                              function of k.
!  For strhopt1=2, function y = c + a*tanh(b*x) is used to specify dz
!                              as a function of k.
!
!-----------------------------------------------------------------------
!
    IF ( strhopt1 == 0 ) THEN

      DO k=1,nz1
        zp1d1(k) = z1(k)
      END DO

    ELSE IF ( strhopt1 == 1 .OR.strhopt1 == 2 ) THEN

      za = zrefsfc1 + MAX(0.0, MIN(dlayer11, z1(nz1-2)-zrefsfc1 ))
      zb = za       + MAX(0.0, MIN(dlayer21, z1(nz1-1)-za      ))

      IF( dlayer11 >= (nz1-3)*dzmin1 ) THEN
        WRITE(6,'(/1x,a,f13.3,/a,f13.3,a,a)')                           &
            'Can not setup a vertical grid with uniform dz=',dzmin1,    &
            ' over the depth of ',dlayer11,' please specify a smaller ', &
            'value of dlayer11. Program stopped ARPSINTRP.'
        STOP 9014
      END IF

      CALL strhgrd(nz1,strhopt1,zrefsfc1,za,zb,z1(nz1-1),               &
                   dzmin1,strhtune1, zp1d1,dzp1d1)

    ELSE

      WRITE(6,'(1x,a,i3,a/)')                                           &
          'Invalid vertical grid stretching option, strhopt1 was ',strhopt1, &
          '. Program stopped in ARPSINTRP.'
      STOP 9015

    END IF
!
!-----------------------------------------------------------------------
!
!  Physical height of computational grid defined as
!
!  Zp=(z-zrefsfc)*(Zm-hterain)/(Zm-zrefsfc)+hterain for z=<Zm.
!  ZP=z for z>Zm
!
!  where Zm the height at which the grid levels becomes flat.
!  Hm < Zm =< Ztop, hm is the height of mounta3din and Ztop the height
!  of model top.
!
!-----------------------------------------------------------------------
!
    DO k=nz1-1,2,-1
      IF(zp1d1(k) <= zflat1) THEN
        zflat11 = zp1d1(k)
        EXIT
      END IF
    END DO

    300   CONTINUE

    zflat11=MAX(MIN(z1(nz1-1),zflat11),zrefsfc1)

    DO k=2,nz1-1

      IF(zp1d1(k) > zflat11) THEN
        DO j=1,ny1-1
          DO i=1,nx1-1
            zp1(i,j,k)=zp1d1(k)
          END DO
        END DO
      ELSE
        DO j=1,ny1-1
          DO i=1,nx1-1
            zp1(i,j,k)=(zp1d1(k)-zrefsfc1)*(zflat11-hterain1(i,j))      &
                       /(zflat11-zrefsfc1)+hterain1(i,j)
          END DO
        END DO
      END IF

    END DO

    DO j=1,ny1-1
      DO i=1,nx1-1
        zp1(i,j,2)=hterain1(i,j)
        zp1(i,j,1)=2.0*zp1(i,j,2)-zp1(i,j,3)
        zp1(i,j,nz1)=2.0*zp1(i,j,nz1-1)-zp1(i,j,nz1-2)
      END DO
    END DO

!  CALL jacob(nx1,ny1,nz1,x11,y11,z1,zp1,j11,j21,j31)

    CALL inisoilgrd(nx1,ny1,nzsoil1,hterain1,zpsoil1, &
                        j3soil1,j3soilinv1)

  END IF    ! samgrd
!
!-----------------------------------------------------------------------
!
!  Write out terrain data
!
!-----------------------------------------------------------------------
!
  CALL getunit( trnchanl )

  ternfn = new_runname(1:lfnkey)//".trndata"
  lternfn = lfnkey + 8

  IF( dirname /= ' ' ) THEN

    temchar = ternfn
    ternfn = dirname(1:ldirnam)//'/'//temchar
    lternfn  = lternfn + ldirnam + 1

  END IF

  CALL fnversn(ternfn, lternfn )

  PRINT *, 'Write terrain data to ',ternfn(1:lternfn)

  CALL asnctl ('NEWLOCAL', 1, ierr)
  CALL asnfile(ternfn, '-F f77 -N ieee', ierr)

  OPEN(trnchanl,FILE=ternfn,FORM='unformatted',STATUS='new')
  WRITE(trnchanl) nx1,ny1
  WRITE(trnchanl) idummy,mapproj,idummy,idummy,idummy,                  &
                  idummy,idummy, idummy,idummy,idummy,                  &
                  idummy,idummy, idummy,idummy,idummy,                  &
                  idummy,idummy, idummy,idummy,idummy

  WRITE(trnchanl)   dx1,    dy1,ctrlat1,ctrlon1,rdummy,                 &
                 rdummy,trulat1,trulat2, trulon,sclfct,                 &
                 rdummy, rdummy, rdummy, rdummy,rdummy,                 &
                 rdummy, rdummy, rdummy, rdummy,rdummy

  WRITE(trnchanl) hterain1

  CLOSE ( trnchanl )
  CALL retunit ( trnchanl )
!
!-----------------------------------------------------------------------
!
!  Open the surface property data file
!
!-----------------------------------------------------------------------
!
  lsfcdmp = .false.
  IF ( landin == 1 ) THEN        ! take care of soil and veg data
    sfcoutfl = new_runname(1:lfnkey)//".sfcdata"
    lfn = lfnkey + 8

    IF( dirname /= ' ' ) THEN

      temchar = sfcoutfl
      sfcoutfl = dirname(1:ldirnam)//'/'//temchar
      lfn  = lfn + ldirnam + 1

    END IF

    CALL fnversn(sfcoutfl, lfn)

    PRINT *, 'Write surface property data in ',sfcoutfl(1:lfn)

    CALL getunit ( sfchanl )
    CALL asnctl ('NEWLOCAL', 1, ierr)
    CALL asnfile(sfcoutfl, '-F f77 -N ieee', ierr)

    OPEN (UNIT=sfchanl,FILE=sfcoutfl(1:lfn), STATUS = 'new',            &
          FORM='unformatted', ACCESS = 'sequential')

    WRITE(sfchanl) nx1,ny1

    WRITE(sfchanl) mapproj,     1,     1,     1,     1,                 &
                         1,idummy, nstyp,idummy,idummy,                 &
                    idummy,idummy,idummy,idummy,idummy,                 &
                    idummy,idummy,idummy,idummy,idummy

    WRITE(sfchanl)     dx1,   dy1,ctrlat1,ctrlon1,trulat1,              &
                   trulat2,trulon, sclfct, rdummy, rdummy,              &
                    rdummy,rdummy, rdummy, rdummy, rdummy,              &
                    rdummy,rdummy, rdummy, rdummy, rdummy

    lsfcdmp = .true.
  END IF
!
!-----------------------------------------------------------------------
!
!  Open the grid base state data file for the new grid
!
!-----------------------------------------------------------------------
!
  CALL gtbasfn(new_runname(1:lfnkey),dirname,ldirnam,hdmpfmt,           &
               1,0,basdmpfn,lbasdmpf)

  PRINT*
  PRINT*,'Output grid/base state file is ', basdmpfn(1:lbasdmpf)
!
!-----------------------------------------------------------------------
!
!  Get the IO unit numbers for input grid and base state file
!
!-----------------------------------------------------------------------
!
  CALL getunit( nchanl )

  CALL asnctl ('NEWLOCAL', 1, ierr)
  CALL asnfile(basdmpfn(1:lbasdmpf), '-F f77 -N ieee', ierr)
!
!-----------------------------------------------------------------------
!
!  Open the output grdbas file
!
!-----------------------------------------------------------------------
!
  OPEN(UNIT=nchanl,FILE=basdmpfn(1:lbasdmpf),                           &
           STATUS='new',FORM='unformatted',IOSTAT=istat)
  IF( istat /= 0 ) GO TO 9002
!
!-----------------------------------------------------------------------
!
!  Read in other base state arrays and interpolate them into new grid
!
!-----------------------------------------------------------------------
!
!  WRITE (cmnt(nocmnt),'(a,i4,a,i4,a,i4)')                               &
!      ' nx =',nx1,', ny =',ny1,', nz =',nz1
  WRITE (cmnt(nocmnt),'(a,i4,a,i4,a,i4,a,i4)')                           &
      ' nx =',nx1,', ny =',ny1,', nz =',nz1,', nzsoil =',nzsoil1

  WRITE(nchanl) fmtver
  WRITE(nchanl) new_runname
  WRITE(nchanl) nocmnt

  IF( nocmnt > 0 ) THEN
    DO i=1,nocmnt
      WRITE(nchanl) cmnt(i)
      WRITE(6,'(1x,a)') cmnt(i)
    END DO
  END IF

  WRITE(nchanl) time,tmunit
  WRITE(nchanl) nx1,ny1,nz1

  WRITE(nchanl)      1,      1,      0,mstout1,      0,                 &
                     0,      0,      0,landout1,totout,                 &
                idummy, idummy, idummy, mapproj, month,                 &
                   day,   year,   hour, minute, second

  WRITE(nchanl) umove,vmove,xgrdorg1,ygrdorg1,trulat1,                  &
                trulat2,trulon,sclfct,rdummy,rdummy,                    &
                rdummy,rdummy,rdummy,rdummy,rdummy,                     &
                tstop,thisdmp,latitud1,ctrlat1,ctrlon1

  IF ( totout /= 0 ) THEN
    WRITE(nchanl) nstyp,  idummy, idummy, idummy, idummy,               &
                  idummy, idummy, idummy, idummy, idummy,               &
                  idummy, idummy, idummy, idummy, idummy,               &
                  idummy, idummy, idummy, idummy, idummy

    WRITE(nchanl) rdummy, rdummy, rdummy, rdummy, rdummy,               &
                  rdummy, rdummy, rdummy, rdummy, rdummy,               &
                  rdummy, rdummy, rdummy, rdummy, rdummy,               &
                  rdummy, rdummy, rdummy, rdummy, rdummy
  END IF

  WRITE(6,'(/1x,a/)')                                                   &
      'Min. and max. of input data interpolated to the new grid:'

  WRITE(nchanl) 'x coord r1d1'
  WRITE(nchanl) x11
  CALL a3dmax0(x11,1,nx1,1,nx1,1,1,1,1,1,1,1,1, amax,amin)
  WRITE(6,'(/1x,2(a,e13.6))') 'xmin    = ', amin,',  xmax    =',amax

  WRITE(nchanl) 'y coord r1d2'
  WRITE(nchanl) y11
  CALL a3dmax0(y11,1,ny1,1,ny1,1,1,1,1,1,1,1,1, amax,amin)
  WRITE(6,'(1x,2(a,e13.6))') 'ymin    = ', amin,',  ymax    =',amax

  WRITE(nchanl) 'z coord r1d3'
  WRITE(nchanl) z1
  CALL a3dmax0(z1,1,nz1,1,nz1,1,1,1,1,1,1,1,1, amax,amin)
  WRITE(6,'(1x,2(a,e13.6))') 'zmin    = ', amin,',  zmax    =',amax

  WRITE(nchanl) 'zp coor r3d0'
  WRITE(nchanl) zp1
  CALL a3dmax0(zp1,1,nx1,1,nx1-1,1,ny1,1,ny1-1,1,nz1,1,nz1,             &
              amax,amin)
  WRITE(6,'(1x,2(a,e13.6))') 'zpmin   = ', amin,', zpmax    =',amax

  WRITE(nchanl) 'zpsoil coor r3d0'
  WRITE(nchanl) zpsoil1
  CALL a3dmax0(zpsoil1,1,nx1,1,nx1-1,1,ny1,1,ny1-1,1,nzsoil1,1,nzsoil1, &
              amax,amin)
  WRITE(6,'(1x,2(a,e13.6))') 'zpmin   = ', amin,', zpmax    =',amax
!
!-----------------------------------------------------------------------
!
!  Read in base state arrays, interpolate them into new grid, and
!  dump them into the new history file
!
!-----------------------------------------------------------------------
!
  READ(inch,ERR=9110,END=9120) label
  READ(inch,ERR=9110,END=9120) a3din
  WRITE(6,'(1x,a,a12,a)')                                               &
      'Field ',label,' was read into array a3din.'
  stgr = 1              ! U-points
  CALL intrp3d( nx,ny,nz,nx1,ny1,nz1,x,y,z,zp,x1,y1,z1,zp1,             &
                samgrd,stgr,a3din,a3dout,tem3d,tem3d1 )
  WRITE(nchanl) 'ubar    r3d1'
  WRITE(nchanl) a3dout
  CALL a3dmax0(a3dout,1,nx1,1,nx1,1,ny1,1,ny1-1,1,nz1,1,nz1-1,          &
               amax,amin)
  WRITE(6,'(1x,2(a,e13.6))') 'ubarmin = ', amin,',  ubarmax =',amax


  READ(inch,ERR=9110,END=9120) label
  READ(inch,ERR=9110,END=9120) a3din
  WRITE(6,'(1x,a,a12,a)')                                               &
      'Field ',label,' was read into array a3din.'
  stgr = 2              ! V-points
  CALL intrp3d( nx,ny,nz,nx1,ny1,nz1,x,y,z,zp,x1,y1,z1,zp1,             &
                samgrd,stgr,a3din,a3dout,tem3d,tem3d1 )
  WRITE(nchanl) 'vbar    r3d1'
  WRITE(nchanl) a3dout
  CALL a3dmax0(a3dout,1,nx1,1,nx1-1,1,ny1,1,ny1,1,nz1,1,nz1-1,          &
               amax,amin)
  WRITE(6,'(1x,2(a,e13.6))') 'vbarmin = ', amin,',  vbarmax =',amax


  READ(inch,ERR=9110,END=9120) label
  READ(inch,ERR=9110,END=9120) a3din
  WRITE(6,'(1x,a,a12,a)')                                               &
      'Field ',label,' was read into array a3din.'
  stgr = 3              ! W-points
  CALL intrp3d( nx,ny,nz,nx1,ny1,nz1,x,y,z,zp,x1,y1,z1,zp1,             &
                samgrd,stgr,a3din,a3dout,tem3d,tem3d1 )
  WRITE(nchanl) 'wbar    r3d1'
  WRITE(nchanl) a3dout

  READ(inch,ERR=9110,END=9120) label
  READ(inch,ERR=9110,END=9120) a3din
  WRITE(6,'(1x,a,a12,a)')                                               &
      'Field ',label,' was read into array a3din.'
  stgr = 4              ! S-points
  CALL intrp3d( nx,ny,nz,nx1,ny1,nz1,x,y,z,zp,x1,y1,z1,zp1,             &
                samgrd,stgr,a3din,a3dout,tem3d,tem3d1 )
  WRITE(nchanl) 'ptbar   r3d1'
  WRITE(nchanl) a3dout
  CALL a3dmax0(a3dout,1,nx1,1,nx1-1,1,ny1,1,ny1-1,1,nz1,1,nz1-1,        &
              amax,amin)
  WRITE(6,'(1x,2(a,e13.6))') 'ptbarmin= ', amin,',  ptbarmax=',amax

  READ(inch,ERR=9110,END=9120) label
  READ(inch,ERR=9110,END=9120) a3din
  WRITE(6,'(1x,a,a12,a)')                                               &
      'Field ',label,' was read into array a3din.'
  stgr = 4              ! S-points
  CALL intrp3d( nx,ny,nz,nx1,ny1,nz1,x,y,z,zp,x1,y1,z1,zp1,             &
                samgrd,stgr,a3din,a3dout,tem3d,tem3d1 )
  WRITE(nchanl) 'pbar    r3d1'
  WRITE(nchanl) a3dout
  CALL a3dmax0(a3dout,1,nx1,1,nx1-1,1,ny1,1,ny1-1,1,nz1,1,nz1-1,        &
              amax,amin)
  WRITE(6,'(1x,2(a,e13.6))') 'pbarmin = ', amin,',  pbarmax =',amax

  IF ( mstin == 1 ) THEN
    READ(inch,ERR=9110,END=9120) label
    READ(inch,ERR=9110,END=9120) a3din
    WRITE(6,'(1x,a,a12,a)')                                             &
        'Field ',label,' was read into array a3din.'

    IF ( mstout1 == 1 ) THEN
      stgr = 4              ! S-points
      CALL intrp3d( nx,ny,nz,nx1,ny1,nz1,x,y,z,zp,x1,y1,z1,zp1,         &
                    samgrd,stgr,a3din,a3dout,tem3d,tem3d1 )
      WRITE(nchanl) 'qvbar   r3d1'
      WRITE(nchanl) a3dout
      CALL a3dmax0(a3dout,1,nx1,1,nx1-1,1,ny1,1,ny1-1,                  &
                   1,nz1,1,nz1-1,amax,amin)
      WRITE(6,'(1x,2(a,e13.6))')                                        &
          'qvbarmin= ', amin,',  qvbarmax=',amax
    END IF
  END IF

  IF ( landin == 1 ) THEN

    IF (nstyp == 1) THEN
      READ(inch,ERR=9110,END=9120) label
      READ(inch,ERR=9110,END=9120) ((i2din(i,j),i=1,nx),j=1,ny)
      WRITE(6,'(1x,a,a12,a)')                                           &
          'Field ',label,' was read into array i2din'
      CALL dist2d(nx,ny,nx1,ny1,x,y,x1,y1,samgrd,i2din,i2dout)
      WRITE(nchanl) 'soiltyp i2d0'
      WRITE(nchanl) ((i2dout(i,j),i=1,nx1),j=1,ny1)
      IF ( lsfcdmp ) WRITE(sfchanl) ((i2dout(i,j),i=1,nx1),j=1,ny1)

    ELSE

      DO is=1,nstyp
        READ(inch,ERR=9110,END=9120) label
        READ(inch,ERR=9110,END=9120)                                    &
              ((i2din(i,j),i=1,nx),j=1,ny)
        WRITE(6,'(1x,a,a12,a)')                                         &
            'Field ',label,' was read into array i2din'
        CALL dist2d(nx,ny,nx1,ny1,x,y,x1,y1,                            &
                    samgrd,i2din,i2dout)
        WRITE(nchanl) 'soiltyp i2d0'
        WRITE(nchanl) ((i2dout(i,j),i=1,nx1),j=1,ny1)
        IF( lsfcdmp ) WRITE(sfchanl) ((i2dout(i,j),i=1,nx1),j=1,ny1)

        READ(inch,ERR=9110,END=9120) label
        READ(inch,ERR=9110,END=9120)                                    &
              ((a2din(i,j),i=1,nx),j=1,ny)
        WRITE(6,'(1x,a,a12,a,i2)')                                      &
            'Field ',label,' was read into array a2din for is=',is
        CALL intrp2d( nx,ny,nx1,ny1,x,y,x1,y1,                          &
                      samgrd,a2din,a2dout )
        WRITE(nchanl) 'stypfrct r2d0'
        WRITE(nchanl) ((a2dout(i,j),i=1,nx1),j=1,ny1)
        IF( lsfcdmp ) WRITE(sfchanl) ((a2dout(i,j),i=1,nx1),j=1,ny1)
      END DO

    END IF

    READ(inch,ERR=9110,END=9120) label
    READ(inch,ERR=9110,END=9120) i2din
    WRITE(6,'(1x,a,a12,a)')                                             &
        'Field ',label,' was read into array i2din'
    CALL dist2d(nx,ny,nx1,ny1,x,y,x1,y1,samgrd,i2din,i2dout)
    WRITE(nchanl) 'vegtyp  i2d0'
    WRITE(nchanl) i2dout
    IF ( lsfcdmp ) WRITE(sfchanl) i2dout

    READ(inch,ERR=9110,END=9120) label
    READ(inch,ERR=9110,END=9120) a2din
    WRITE(6,'(1x,a,a12,a)')                                             &
        'Field ',label,' was read into array a2din.'
    CALL intrp2d( nx,ny,nx1,ny1,x,y,x1,y1, samgrd,a2din,a2dout )
    WRITE(nchanl) 'lai     r2d0'
    WRITE(nchanl)  a2dout
    IF ( lsfcdmp ) WRITE(sfchanl) a2dout

    READ(inch,ERR=9110,END=9120) label
    READ(inch,ERR=9110,END=9120) a2din
    WRITE(6,'(1x,a,a12,a)')                                             &
        'Field ',label,' was read into array a2din.'
    CALL intrp2d( nx,ny,nx1,ny1,x,y,x1,y1, samgrd,a2din,a2dout )
    WRITE(nchanl) 'roufns  r2d0'
    WRITE(nchanl)  a2dout
    IF ( lsfcdmp ) WRITE(sfchanl) a2dout

    READ(inch,ERR=9110,END=9120) label
    READ(inch,ERR=9110,END=9120) a2din
    WRITE(6,'(1x,a,a12,a)')                                             &
        'Field ',label,' was read into array a2din.'
    CALL intrp2d( nx,ny,nx1,ny1,x,y,x1,y1, samgrd,a2din,a2dout )
    WRITE(nchanl) 'veg     r2d0'
    WRITE(nchanl)  a2dout
    IF ( lsfcdmp ) WRITE(sfchanl) a2dout

  END IF

  CLOSE (inch)
  CLOSE (nchanl)
  IF ( lsfcdmp ) CLOSE (sfchanl)

  CALL retunit( inch )
  CALL retunit( nchanl )
  IF ( lsfcdmp ) CALL retunit (sfchanl)
!
!-----------------------------------------------------------------------
!
!  Loop over time dependent data files
!
!-----------------------------------------------------------------------
!
  grdbas = 0

  DO nfile = 1,nhisfile

    filename = hisfile(nfile)

    lenfil=LEN_trim(filename)
    WRITE(6,'(/a,a,a)')                                                 &
        ' Data set ', filename(1:lenfil) ,' to be processed.'

    INQUIRE(FILE=filename(1:lenfil), EXIST = fexist )
    IF( fexist ) THEN
      GO TO 370
    END IF

    WRITE(6,'(a/a)')                                                    &
        'File '//filename(1:lenfil)//' does not exist.',                &
        'Check if compressed file '//filename(1:lenfil)//'.Z'           &
        //' exists.'

    INQUIRE(FILE=filename(1:lenfil)//'.Z',EXIST = fexist )
    IF( fexist ) THEN
      CALL uncmprs( filename(1:lenfil)//'.Z' )
      GO TO 370
    END IF

    WRITE(6,'(a/a)')                                                    &
        'File '//filename(1:lenfil)//'.Z'//' does not exist.',          &
        'Check if compressed file '//filename(1:lenfil)//'.gz'          &
        //' exists.'

    INQUIRE(FILE=filename(1:lenfil)//'.gz',EXIST = fexist )
    IF( .NOT.fexist ) THEN
      CALL uncmprs( filename(1:lenfil)//'.gz' )
      GO TO 370
    END IF

    WRITE(6,'(a/a)')                                                    &
        'File '//filename(1:lenfil)                                     &
        //' or its compressed version not found.',                      &
        'Program stopped.'
    STOP

    370   CONTINUE         ! also continue to read another time recode
                           ! from GrADS file
!
!-----------------------------------------------------------------------
!
!  Get the IO unit numbers for input files
!
!-----------------------------------------------------------------------
!
    CALL getunit( inch )
!
!-----------------------------------------------------------------------
!
!  Cray routines to force binary data file to be in the IEEE format
!
!-----------------------------------------------------------------------
!
    CALL asnctl ('NEWLOCAL', 1, ierr)
    CALL asnfile(filename(1:lenfil), '-F f77 -N ieee', ierr)
!
!-----------------------------------------------------------------------
!
!  Open the input data file
!
!-----------------------------------------------------------------------
!
    OPEN(UNIT=inch,FILE=filename(1:lenfil),                             &
             STATUS='old',FORM='unformatted',IOSTAT=istat)
    IF( istat /= 0 ) GO TO 9003
!
!-----------------------------------------------------------------------
!
!  Read all input data header
!
!-----------------------------------------------------------------------
!
    READ(inch,ERR=9115,END=9125) fmtverin

    IF ( fmtverin == fmtver ) THEN
      oldver = 0
    ELSE
      oldver = 1
    END IF

    READ(inch,ERR=9115,END=9125) runname
    READ(inch,ERR=9115,END=9125) nocmnt

    IF( nocmnt > 0 ) THEN
      DO i=1,nocmnt
        READ(inch,ERR=9115,END=9125) cmnt(i)
      END DO
    END IF

    READ(inch,ERR=9115,END=9125) time,tmunit
!    READ(inch,ERR=9115,END=9125) nxin, nyin, nzin
    READ(inch,ERR=9115,END=9125) nxin, nyin, nzin, nzsoilin

!    IF( nxin /= nx .OR. nyin /= ny .OR. nzin /= nz ) THEN
!      WRITE(6,'(1x,a)')                                                 &
!           ' Dimensions in ARPSINTRP inconsistent with data.'
!      WRITE(6,'(1x,a,3I15)') ' Read were: ', nxin, nyin, nzin
!      WRITE(6,'(1x,a)')                                                 &
!           ' Program aborted in ARPSINTRP.'
!      STOP
!    END IF
    IF( nxin /= nx .OR. nyin /= ny .OR. nzin /= nz .OR. &
        nzsoilin /= nzsoil ) THEN
      WRITE(6,'(1x,a)')                                                 &
           ' Dimensions in ARPSINTRP inconsistent with data.'
      WRITE(6,'(1x,a,4I15)') ' Read were: ', nxin, nyin, nzin, nzsoilin
      WRITE(6,'(1x,a)')                                                 &
           ' Program aborted in ARPSINTRP.'
      STOP
    END IF

    WRITE(6,'(1x,a,f8.1,a,f8.3,a/)')                                    &
           'To read data at time ', time,                               &
           ' secs = ',(time/60.),' mins.'

    READ(inch,ERR=9115,END=9125)                                        &
           grdin,basin,varin,mstin,icein,                               &
           trbin, sfcin,rainin,landin,totin,                            &
           tkein,idummy,idummy,mapproj, month,                          &
           day, year, hour, minute, second

    IF ( varin == 1 ) THEN
      varout1  = varout
    ELSE
      varout1  = 0
    END IF

    IF ( mstin == 1 ) THEN
      mstout1  = mstout
    ELSE
      mstout1 = 0
    END IF

    IF ( icein == 1 ) THEN
      iceout1  = iceout
    ELSE
      iceout1  = 0
    END IF

    IF ( tkein == 1 ) THEN
      tkeout1  = tkeout
    ELSE
      tkeout1  = 0
    END IF

    IF ( trbin == 1 ) THEN
      trbout1  = trbout
    ELSE
      trbout1  = 0
    END IF

    IF ( rainin == 1 ) THEN
      rainout1 = rainout
    ELSE
      rainout1 = 0
    END IF

    IF ( sfcin == 1 ) THEN
      sfcout1  = sfcout
    ELSE
      sfcout1  = 0
    END IF

    IF ( landin == 1 ) THEN
      landout1 = landout
    ELSE
      landout1 = 0
    END IF

    READ(inch,ERR=9115,END=9125)                                        &
                    umove,vmove,xgrdorg,ygrdorg,trulat1,                &
                    trulat2,trulon,sclfct,rdummy,rdummy,                &
                    rdummy,rdummy,rdummy,rdummy,rdummy,                 &
                    tstop,thisdmp,latitud,ctrlat,ctrlon

    IF ( totin /= 0 ) THEN
!
!-----------------------------------------------------------------------
!
!  Read in additional parameters for ARPS history dump 4.0 or later
!  version.
!
!-----------------------------------------------------------------------
!
      READ(inch,ERR=9115,END=9125)                                      &
           nstyp,  prcin,idummy,idummy,idummy,                          &
           idummy,idummy,idummy,idummy,idummy,                          &
           idummy,idummy,idummy,idummy,idummy,                          &
           idummy,idummy,idummy,idummy,idummy

      IF ( nstyp < 1 ) THEN
        nstyp = 1
      END IF

      IF ( prcin == 1 ) THEN
        prcout1 = prcout
      ELSE
        prcout1 = 0
      END IF

      READ(inch,ERR=9115,END=9125)                                      &
           rdummy,rdummy,rdummy,rdummy,rdummy,                          &
           rdummy,rdummy,rdummy,rdummy,rdummy,                          &
           rdummy,rdummy,rdummy,rdummy,rdummy,                          &
           rdummy,rdummy,rdummy,rdummy,rdummy
    END IF

    curtim = time
!
!-----------------------------------------------------------------------
!
!  Get the history file name for the new grid at time curtim
!
!-----------------------------------------------------------------------
!
    runname = new_runname

    CALL gtdmpfn(runname(1:lfnkey),dirname,                             &
                 ldirnam,curtim,houtfmt,1,0, hdmpfn, ldmpf)

    CALL getunit( nchanl )
    CALL asnctl ('NEWLOCAL', 1, ierr)
    CALL asnfile(hdmpfn(1:ldmpf), '-F f77 -N ieee', ierr)
    OPEN(UNIT=nchanl,FILE=hdmpfn(1:ldmpf),                              &
             STATUS='new',FORM='unformatted',IOSTAT=istat)
    IF( istat /= 0 ) GO TO 9004
!
!-----------------------------------------------------------------------
!
!  Open the surface property data file
!
!-----------------------------------------------------------------------
!
    lsoildmp = .false.
    IF ( sfcin == 1 ) THEN        ! take care of soil variables
      CALL cvttsnd( curtim, timsnd, tmstrln )

      soiloutfl = runname(1:lfnkey)//".soilvar."//timsnd(1:tmstrln)
      lfn = lfnkey + 9 + tmstrln

      IF( dirname /= ' ' ) THEN

        temchar = soiloutfl
        soiloutfl = dirname(1:ldirnam)//'/'//temchar
        lfn  = lfn + ldirnam + 1

      END IF

      CALL fnversn(soiloutfl, lfn)

      PRINT *, 'Write soil initial data to ',soiloutfl(1:lfn)

      CALL getunit ( soilchanl )
      CALL asnctl ('NEWLOCAL', 1, ierr)
      CALL asnfile(soiloutfl, '-F f77 -N ieee', ierr)

      OPEN (UNIT=soilchanl,FILE=soiloutfl(1:lfn), STATUS = 'new',       &
            FORM='unformatted', ACCESS = 'sequential')

      WRITE(soilchanl) nx1,ny1

!      WRITE(soilchanl) mapproj,     1,     1,     1,     1,             &
!                           1,idummy, nstyp,idummy,idummy,               &
!                      idummy,idummy,idummy,idummy,idummy,               &
!                      idummy,idummy,idummy,idummy,idummy
      WRITE(soilchanl) mapproj,     1,     1,                           &
                           1,idummy, nstyp,idummy,idummy,               &
                      idummy,idummy,idummy,idummy,idummy,               &
                      idummy,idummy,idummy,idummy,idummy

      WRITE(soilchanl)     dx1,   dy1,ctrlat1,ctrlon1,trulat1,          &
                     trulat2,trulon, sclfct, rdummy, rdummy,            &
                      rdummy,rdummy, rdummy, rdummy, rdummy,            &
                      rdummy,rdummy, rdummy, rdummy, rdummy

      lsoildmp = .true.
    END IF

    IF ( exbcdmp == 1 ) THEN
!
!-----------------------------------------------------------------------
!
!  Open the EXBC file to dump the data
!
!-----------------------------------------------------------------------
!
      dmpexbc = .true.
      IF ( varin == 0 ) THEN
        WRITE(6,'(1x,a,a)')                                             &
            'The data file does not contain the variables needed ',     &
            'for EXBC data file'
        dmpexbc = .false.
        GO TO 385
      END IF

      IF ( mstout1 == 0 ) THEN
        qcexout = 0
        qrexout = 0
        qiexout = 0
        qsexout = 0
        qhexout = 0
      ELSE IF ( iceout1 == 0 ) THEN
        qiexout = 0
        qsexout = 0
        qhexout = 0
      END IF

      CALL ctim2abss( year,month,day,hour,minute,second, itema )

      itema = itema + INT(curtim)

      CALL abss2ctim( itema,iyr,imon,idy,ihr,imin,isec )

      WRITE (ctime,'(i4.4,2i2.2,a,3i2.2)')                              &
            iyr,imon,idy,'.',ihr,imin,isec

      exbcfn = runname(1:lfnkey)//'.'//ctime
      length = lfnkey + 16
      IF( dirname /= ' ' ) THEN
        temchar = exbcfn
        exbcfn  = dirname(1:ldirnam)//'/'//temchar
        length  = length + ldirnam + 1
      END IF

      CALL fnversn(exbcfn,length)

      WRITE(6,'(1x,a,a)')                                               &
           'The external boundary data file is ',exbcfn(1:length)

      CALL getunit( exbchanl )

      CALL asnctl ('NEWLOCAL', 1, ierr)
      CALL asnfile(exbcfn(1:length), '-F f77 -N ieee', ierr)
      OPEN(UNIT=exbchanl,FILE=exbcfn(1:length),                         &
           STATUS='new',FORM='unformatted',IOSTAT=istat)
      IF( istat /= 0 ) GO TO 9005

      WRITE(exbchanl) nx1,ny1,nz1,dx1,dy1,dz1,ctrlat1,ctrlon1
      WRITE(exbchanl) ctime

      WRITE(exbchanl) varin,varin,varin,varin,varin,                    &
                      mstout1,qcexout,qrexout,qiexout,qsexout,          &
                      qhexout,idummy,idummy,idummy,idummy,              &
                      idummy,idummy,idummy,idummy,idummy,               &
                      idummy,idummy,idummy,idummy,idummy,               &
                      idummy,idummy,idummy,idummy,idummy,               &
                      idummy,idummy,idummy,idummy,idummy,               &
                      idummy,idummy,idummy,idummy,idummy
    END IF

    385   CONTINUE

!
!-----------------------------------------------------------------------
!
!  Write out the header information for new grid
!
!-----------------------------------------------------------------------
!
    WRITE(nchanl) fmtver
    WRITE(nchanl) new_runname
    WRITE(nchanl) nocmnt

    WRITE (cmnt(nocmnt),'(a,i4,a,i4,a,i4)')                             &
        ' nx =',nx1,', ny =',ny1,', nz =',nz1

    IF( nocmnt > 0 ) THEN
      DO i=1,nocmnt
        WRITE(nchanl) cmnt(i)
        WRITE(6,'(1x,a)') cmnt(i)
      END DO
    END IF

    WRITE(nchanl) time,tmunit

    WRITE(nchanl) nx1,ny1,nz1
!
!-----------------------------------------------------------------------
!
!  Write the flags for different data groups.
!
!-----------------------------------------------------------------------
!
    WRITE(nchanl) grdout,      0, varout1, mstout1, iceout1,            &
                  trbout1, sfcout1,rainout1,landout1,totout,            &
                  tkeout1, idummy, idummy, mapproj, month,              &
                     day,   year,   hour, minute, second
    rdummy = 0.0
    WRITE(nchanl)   umove,   vmove,xgrdorg1,ygrdorg1, trulat1,          &
                  trulat2,  trulon,  sclfct,  rdummy,  rdummy,          &
                   rdummy,  rdummy,  rdummy,  rdummy,  rdummy,          &
                    tstop, thisdmp,latitud1, ctrlat1, ctrlon1
!
!-----------------------------------------------------------------------
!
!  If totout=1, write additional parameters to history dump files.
!  This is for ARPS version 4.1.2 or later.
!
!-----------------------------------------------------------------------
!
    IF ( totout == 1 ) THEN
      WRITE(nchanl) nstyp, prcout1, idummy, idummy, idummy,             &
                    idummy, idummy, idummy, idummy, idummy,             &
                    idummy, idummy, idummy, idummy, idummy,             &
                    idummy, idummy, idummy, idummy, idummy

      WRITE(nchanl) rdummy, rdummy, rdummy, rdummy, rdummy,             &
                    rdummy, rdummy, rdummy, rdummy, rdummy,             &
                    rdummy, rdummy, rdummy, rdummy, rdummy,             &
                    rdummy, rdummy, rdummy, rdummy, rdummy
    END IF
!
!-----------------------------------------------------------------------
!
!  Read in x,y and z at grid cell centers (scalar points).
!
!----------------------------------------------------------------------
!
    IF( grdin == 1 .OR. grdbas == 1 ) THEN
      READ(inch,ERR=9115,END=9125) label
      READ(inch,ERR=9115,END=9125) x
      WRITE(6,'(1x,a,a12,a)')                                           &
          'Field ',label,' was read into array x.'

      READ(inch,ERR=9115,END=9125) label
      READ(inch,ERR=9115,END=9125) y
      WRITE(6,'(1x,a,a12,a)')                                           &
          'Field ',label,' was read into array y.'

      READ(inch,ERR=9115,END=9125) label
      READ(inch,ERR=9115,END=9125) z
      WRITE(6,'(1x,a,a12,a)')                                           &
          'Field ',label,' was read into array z.'

      READ(inch,ERR=9115,END=9125) label
      READ(inch,ERR=9115,END=9125) zp
      WRITE(6,'(1x,a,a12,a)')                                           &
          'Field ',label,' was read into array zp.'

      READ(inch,ERR=9115,END=9125) label
      READ(inch,ERR=9115,END=9125) zpsoil
      WRITE(6,'(1x,a,a12,a)')                                           &
          'Field ',label,' was read into array zpsoil.'
    END IF  ! grdin
!
!-----------------------------------------------------------------------
!
!  If grdout=1 or grdbas=1, write out grid variables
!
!-----------------------------------------------------------------------
!
    IF(grdout == 1 .OR. grdbas == 1 ) THEN

      WRITE(nchanl) 'x coord r1d1'
      WRITE(nchanl) x11

      WRITE(nchanl) 'y coord r1d2'
      WRITE(nchanl) y11

      WRITE(nchanl) 'z coord r1d3'
      WRITE(nchanl) z1

      WRITE(nchanl) 'zp coor r3d0'
      WRITE(nchanl) zp1

      WRITE(nchanl) 'zpsoil coor r3d0'
      WRITE(nchanl) zpsoil1

    END IF    ! grdout
!
!-----------------------------------------------------------------------
!
!  Read in base state fields
!
!----------------------------------------------------------------------
!
    IF( basin == 1 .OR. grdbas == 1 ) THEN

      READ(inch,ERR=9115,END=9125) label
      READ(inch,ERR=9115,END=9125) a3din
      WRITE(6,'(1x,a,a12,a)')                                           &
          'Field ',label,' was read into array a3din.'

      READ(inch,ERR=9115,END=9125) label
      READ(inch,ERR=9115,END=9125) a3din
      WRITE(6,'(1x,a,a12,a)')                                           &
          'Field ',label,' was read into array a3din.'

      READ(inch,ERR=9115,END=9125) label
      READ(inch,ERR=9115,END=9125) a3din
      WRITE(6,'(1x,a,a12,a)')                                           &
          'Field ',label,' was read into array a3din.'


      READ(inch,ERR=9115,END=9125) label
      READ(inch,ERR=9115,END=9125) a3din
      WRITE(6,'(1x,a,a12,a)')                                           &
          'Field ',label,' was read into array a3din.'


      READ(inch,ERR=9115,END=9125) label
      READ(inch,ERR=9115,END=9125) a3din
      WRITE(6,'(1x,a,a12,a)')                                           &
          'Field ',label,' was read into array a3din.'

      IF( mstin == 1) THEN
        READ(inch,ERR=9115,END=9125) label
        READ(inch,ERR=9115,END=9125) a3din
        WRITE(6,'(1x,a,a12,a)')                                         &
            'Field ',label,' was read into array a3din.'

      END IF

      IF (landin == 1) THEN

        IF (nstyp == 1) THEN

          READ(inch,ERR=9115,END=9125) label
          READ(inch,ERR=9115,END=9125)                                  &
                  ((i2din(i,j),i=1,nx),j=1,ny)
          WRITE(6,'(1x,a,a12,a,i2)')                                    &
              'Field ',label,' was read into array soiltyp for is=',is

        ELSE

          DO is=1,nstyp
            READ(inch,ERR=9115,END=9125) label
            READ(inch,ERR=9115,END=9125)                                &
                  ((i2din(i,j),i=1,nx),j=1,ny)
            WRITE(6,'(1x,a,a12,a,i2)')                                  &
                'Field ',label,' was read into array i2din for is=',is

            READ(inch,ERR=9115,END=9125) label
            READ(inch,ERR=9115,END=9125)                                &
                  ((a2din(i,j),i=1,nx),j=1,ny)
            WRITE(6,'(1x,a,a12,a,i2)')                                  &
                'Field ',label,' was read into array a2din for is=',is
          END DO

        END IF

        READ(inch,ERR=9115,END=9125) label
        READ(inch,ERR=9115,END=9125) i2din
        WRITE(6,'(1x,a,a12,a,i2)')                                      &
            'Field ',label,' was read into array i2din for is=',is

        READ(inch,ERR=9115,END=9125) label
        READ(inch,ERR=9115,END=9125) a2din
        WRITE(6,'(1x,a,a12,a,i2)')                                      &
            'Field ',label,' was read into array a2din for is=',is

        READ(inch,ERR=9115,END=9125) label
        READ(inch,ERR=9115,END=9125) a2din
        WRITE(6,'(1x,a,a12,a,i2)')                                      &
            'Field ',label,' was read into array a2din for is=',is

        READ(inch,ERR=9115,END=9125) label
        READ(inch,ERR=9115,END=9125) a2din
        WRITE(6,'(1x,a,a12,a,i2)')                                      &
            'Field ',label,' was read into array a2din for is=',is

      END IF

    END IF
!
!-----------------------------------------------------------------------
!
!  Read in arrays for interpolations
!
!----------------------------------------------------------------------
!
    READ(inch,ERR=9115,END=9125) label
    READ(inch,ERR=9115,END=9125) a3din
    WRITE(6,'(1x,a,a12,a)')                                             &
        'Field ',label,' was read into array a3din.'
    stgr = 1              ! U-points
    CALL intrp3d( nx,ny,nz,nx1,ny1,nz1,x,y,z,zp,x1,y1,z1,zp1,           &
                  samgrd,stgr,a3din,a3dout,tem3d,tem3d1 )
    WRITE(nchanl) 'u       r3d1'
    WRITE(nchanl) a3dout
    CALL a3dmax0(a3dout,1,nx1,1,nx1,1,ny1,1,ny1-1,1,nz1,1,nz1-1,        &
                 amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'umin    = ', amin,',  umax    =',amax
    IF ( dmpexbc ) WRITE(exbchanl) a3dout

    READ(inch,ERR=9115,END=9125) label
    READ(inch,ERR=9115,END=9125) a3din
    WRITE(6,'(1x,a,a12,a)')                                             &
        'Field ',label,' was read into array a3din.'
    stgr = 2              ! V-points
    CALL intrp3d( nx,ny,nz,nx1,ny1,nz1,x,y,z,zp,x1,y1,z1,zp1,           &
                  samgrd,stgr,a3din,a3dout,tem3d,tem3d1 )
    WRITE(nchanl) 'v       r3d1'
    WRITE(nchanl) a3dout
    CALL a3dmax0(a3dout,1,nx1,1,nx1-1,1,ny1,1,ny1,1,nz1,1,nz1-1,        &
                 amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'vmin    = ', amin,',  vmax    =',amax
    IF ( dmpexbc ) WRITE(exbchanl) a3dout

    READ(inch,ERR=9115,END=9125) label
    READ(inch,ERR=9115,END=9125) a3din
    WRITE(6,'(1x,a,a12,a)')                                             &
        'Field ',label,' was read into array a3din.'
    stgr = 3              ! W-points
    CALL intrp3d( nx,ny,nz,nx1,ny1,nz1,x,y,z,zp,x1,y1,z1,zp1,           &
                  samgrd,stgr,a3din,a3dout,tem3d,tem3d1 )
    WRITE(nchanl) 'w       r3d1'
    WRITE(nchanl) a3dout
    CALL a3dmax0(a3dout,1,nx1,1,nx1-1,1,ny1,1,ny1,1,nz1,1,nz1-1,        &
                 amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'wmin    = ', amin,',  wmax    =',amax
    IF ( dmpexbc ) WRITE(exbchanl) a3dout

    READ(inch,ERR=9115,END=9125) label
    READ(inch,ERR=9115,END=9125) a3din
    WRITE(6,'(1x,a,a12,a)')                                             &
        'Field ',label,' was read into array a3din.'
    stgr = 4              ! S-points
    CALL intrp3d( nx,ny,nz,nx1,ny1,nz1,x,y,z,zp,x1,y1,z1,zp1,           &
                  samgrd,stgr,a3din,a3dout,tem3d,tem3d1 )
    WRITE(nchanl) 'pt      r3d1'
    WRITE(nchanl) a3dout
    CALL a3dmax0(a3dout,1,nx1,1,nx1-1,1,ny1,1,ny1-1,1,nz1,1,nz1-1,      &
                 amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'ptmin   = ', amin,',  ptmax   =',amax
    IF ( dmpexbc ) WRITE(exbchanl) a3dout

    READ(inch,ERR=9115,END=9125) label
    READ(inch,ERR=9115,END=9125) a3din
    WRITE(6,'(1x,a,a12,a)')                                             &
        'Field ',label,' was read into array a3din.'
    stgr = 4              ! S-points
    CALL intrp3d( nx,ny,nz,nx1,ny1,nz1,x,y,z,zp,x1,y1,z1,zp1,           &
                  samgrd,stgr,a3din,a3dout,tem3d,tem3d1 )
    WRITE(nchanl) 'p       r3d1'
    WRITE(nchanl) a3dout
    CALL a3dmax0(a3dout,1,nx1,1,nx1-1,1,ny1,1,ny1-1,1,nz1,1,nz1-1,      &
                 amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'pmin    = ', amin,',  pmax    =',amax
    IF ( dmpexbc ) WRITE(exbchanl) a3dout

    IF( mstin == 1 ) THEN

      READ(inch,ERR=9115,END=9125) label
      READ(inch,ERR=9115,END=9125) a3din
      WRITE(6,'(1x,a,a12,a)')                                           &
          'Field ',label,' was read into array a3din.'
      IF ( mstout == 1 ) THEN
        stgr = 4              ! S-points
        CALL intrp3d( nx,ny,nz,nx1,ny1,nz1,x,y,z,zp,x1,y1,z1,zp1,       &
                      samgrd,stgr,a3din,a3dout,tem3d,tem3d1 )
        WRITE(nchanl) 'qv      r3d1'
        WRITE(nchanl) a3dout
        CALL a3dmax0(a3dout,1,nx1,1,nx1-1,1,ny1,1,ny1-1,                &
                     1,nz1,1,nz1-1,amax,amin)
        WRITE(6,'(1x,2(a,e13.6))')                                      &
            'qvmin   = ', amin,',  qvmax   =',amax
        IF ( dmpexbc ) WRITE(exbchanl) a3dout
      END IF

      READ(inch,ERR=9115,END=9125) label
      READ(inch,ERR=9115,END=9125) a3din
      WRITE(6,'(1x,a,a12,a)')                                           &
          'Field ',label,' was read into array a3din.'
      IF ( mstout == 1 ) THEN
        stgr = 4              ! S-points
        CALL intrp3d( nx,ny,nz,nx1,ny1,nz1,x,y,z,zp,x1,y1,z1,zp1,       &
                      samgrd,stgr,a3din,a3dout,tem3d,tem3d1 )
        WRITE(nchanl) 'qc      r3d1'
        WRITE(nchanl) a3dout
        CALL a3dmax0(a3dout,1,nx1,1,nx1-1,1,ny1,1,ny1-1,                &
                     1,nz1,1,nz1-1,amax,amin)
        WRITE(6,'(1x,2(a,e13.6))')                                      &
            'qcmin   = ', amin,',  qcmax   =',amax
        IF ( dmpexbc .AND. qcexout == 1 ) WRITE(exbchanl) a3dout
      END IF

      READ(inch,ERR=9115,END=9125) label
      READ(inch,ERR=9115,END=9125) a3din
      WRITE(6,'(1x,a,a12,a)')                                           &
          'Field ',label,' was read into array a3din.'
      IF ( mstout == 1 ) THEN
        stgr = 4              ! S-points
        CALL intrp3d( nx,ny,nz,nx1,ny1,nz1,x,y,z,zp,x1,y1,z1,zp1,       &
                      samgrd,stgr,a3din,a3dout,tem3d,tem3d1 )
        WRITE(nchanl) 'qr      r3d1'
        WRITE(nchanl) a3dout
        CALL a3dmax0(a3dout,1,nx1,1,nx1-1,1,ny1,1,ny1-1,                &
                     1,nz1,1,nz1-1,amax,amin)
        WRITE(6,'(1x,2(a,e13.6))')                                      &
            'qrmin   = ', amin,',  qrmax   =',amax
        IF ( dmpexbc .AND. qrexout == 1 ) WRITE(exbchanl) a3dout
      END IF

      IF( rainin == 1 ) THEN
        READ(inch,ERR=9115,END=9125) label
        READ(inch,ERR=9115,END=9125) a2din
        WRITE(6,'(1x,a,a12,a)')                                         &
            'Field ',label,' was read into array a2din.'
        IF ( mstout == 1 .AND. rainout == 1 ) THEN
          CALL intrp2d( nx,ny,nx1,ny1,x,y,x1,y1,                        &
                        samgrd,a2din,a2dout )
          WRITE(nchanl) 'raing   r3d1'
          WRITE(nchanl) a2dout
          CALL a3dmax0(a2dout,1,nx1,1,nx1-1,1,ny1,1,ny1-1,              &
                       1,1,1,1,amax,amin)
          WRITE(6,'(1x,2(a,e13.6))')                                    &
              'raingmin= ', amin,',  raingmax=',amax
        END IF

        READ(inch,ERR=9115,END=9125) label
        READ(inch,ERR=9115,END=9125) a2din
        WRITE(6,'(1x,a,a12,a)')                                         &
            'Field ',label,' was read into array a2din.'
        IF ( mstout == 1 .AND. rainout == 1 ) THEN
          CALL intrp2d( nx,ny,nx1,ny1,x,y,x1,y1,                        &
                        samgrd,a2din,a2dout )
          WRITE(nchanl) 'rainc   r3d1'
          WRITE(nchanl) a2dout
          CALL a3dmax0(a2dout,1,nx1,1,nx1-1,1,ny1,1,ny1-1,              &
                       1,1,1,1,amax,amin)
          WRITE(6,'(1x,2(a,e13.6))')                                    &
              'raincmin= ', amin,',  raincmax=',amax
        END IF
      END IF

      IF( prcin == 1 ) THEN
        READ(inch,ERR=9115,END=9125) label
        READ(inch,ERR=9115,END=9125) a2din
        WRITE(6,'(1x,a,a12,a)')                                         &
            'Field ',label,' was read into array a2din.'
        IF ( mstout == 1 .AND. prcout == 1 ) THEN
          CALL intrp2d( nx,ny,nx1,ny1,x,y,x1,y1,                        &
                        samgrd,a2din,a2dout )
          WRITE(nchanl) 'prcrt1  r3d1'
          WRITE(nchanl) a2dout
          CALL a3dmax0(a2dout,1,nx1,1,nx1-1,1,ny1,1,ny1-1,              &
                       1,1,1,1,amax,amin)
          WRITE(6,'(1x,2(a,e13.6))')                                    &
              'prcr1min= ', amin,',  prcr1max=',amax
        END IF

        READ(inch,ERR=9115,END=9125) label
        READ(inch,ERR=9115,END=9125) a2din
        WRITE(6,'(1x,a,a12,a)')                                         &
            'Field ',label,' was read into array a2din.'
        IF ( mstout == 1 .AND. prcout == 1 ) THEN
          CALL intrp2d( nx,ny,nx1,ny1,x,y,x1,y1,                        &
                        samgrd,a2din,a2dout )
          WRITE(nchanl) 'prcrt2  r3d1'
          WRITE(nchanl) a2dout
          CALL a3dmax0(a2dout,1,nx1,1,nx1-1,1,ny1,1,ny1-1,              &
                       1,1,1,1,amax,amin)
          WRITE(6,'(1x,2(a,e13.6))')                                    &
              'prcr2min= ', amin,',  prcr1max=',amax
        END IF

        READ(inch,ERR=9115,END=9125) label
        READ(inch,ERR=9115,END=9125) a2din
        WRITE(6,'(1x,a,a12,a)')                                         &
            'Field ',label,' was read into array a2din.'
        IF ( mstout == 1 .AND. prcout == 1 ) THEN
          CALL intrp2d( nx,ny,nx1,ny1,x,y,x1,y1,                        &
                        samgrd,a2din,a2dout )
          WRITE(nchanl) 'prcrt3  r3d1'
          WRITE(nchanl) a2dout
          CALL a3dmax0(a2dout,1,nx1,1,nx1-1,1,ny1,1,ny1-1,              &
                       1,1,1,1,amax,amin)
          WRITE(6,'(1x,2(a,e13.6))')                                    &
              'prcr3min= ', amin,',  prcr1max=',amax
        END IF

        READ(inch,ERR=9115,END=9125) label
        READ(inch,ERR=9115,END=9125) a2din
        WRITE(6,'(1x,a,a12,a)')                                         &
            'Field ',label,' was read into array a2din.'
        IF ( mstout == 1 .AND. prcout == 1 ) THEN
          CALL intrp2d( nx,ny,nx1,ny1,x,y,x1,y1,                        &
                        samgrd,a2din,a2dout )
          WRITE(nchanl) 'prcrt4  r3d1'
          WRITE(nchanl) a2dout
          CALL a3dmax0(a2dout,1,nx1,1,nx1-1,1,ny1,1,ny1-1,              &
                       1,1,1,1,amax,amin)
          WRITE(6,'(1x,2(a,e13.6))')                                    &
              'prcr4min= ', amin,',  prcr1max=',amax
        END IF

      END IF

      IF( icein == 1 ) THEN

        READ(inch,ERR=9115,END=9125) label
        READ(inch,ERR=9115,END=9125) a3din
        WRITE(6,'(1x,a,a12,a)')                                         &
            'Field ',label,' was read into array a3din.'
        IF ( mstout == 1 .AND. iceout == 1 ) THEN
          stgr = 4              ! S-points
          CALL intrp3d( nx,ny,nz,nx1,ny1,nz1,x,y,z,zp,x1,y1,z1,zp1,     &
                        samgrd,stgr,a3din,a3dout,tem3d,tem3d1 )
          WRITE(nchanl) 'qi      r3d1'
          WRITE(nchanl) a3dout
          CALL a3dmax0(a3dout,1,nx1,1,nx1-1,1,ny1,1,ny1-1,              &
                       1,nz1,1,nz1-1,amax,amin)
          WRITE(6,'(1x,2(a,e13.6))')                                    &
              'qimin   = ', amin,',  qimax   =',amax
          IF ( dmpexbc .AND. qiexout == 1 ) WRITE(exbchanl) a3dout
        END IF

        READ(inch,ERR=9115,END=9125) label
        READ(inch,ERR=9115,END=9125) a3din
        WRITE(6,'(1x,a,a12,a)')                                         &
            'Field ',label,' was read into array a3din.'
        IF ( mstout == 1 .AND. iceout == 1 ) THEN
          stgr = 4              ! S-points
          CALL intrp3d( nx,ny,nz,nx1,ny1,nz1,x,y,z,zp,x1,y1,z1,zp1,     &
                        samgrd,stgr,a3din,a3dout,tem3d,tem3d1 )
          WRITE(nchanl) 'qs      r3d1'
          WRITE(nchanl) a3dout
          CALL a3dmax0(a3dout,1,nx1,1,nx1-1,1,ny1,1,ny1-1,              &
                       1,nz1,1,nz1-1,amax,amin)
          WRITE(6,'(1x,2(a,e13.6))')                                    &
              'qsmin   = ', amin,',  qsmax   =',amax
          IF ( dmpexbc .AND. qsexout == 1 ) WRITE(exbchanl) a3dout
        END IF

        READ(inch,ERR=9115,END=9125) label
        READ(inch,ERR=9115,END=9125) a3din
        WRITE(6,'(1x,a,a12,a)')                                         &
            'Field ',label,' was read into array a3din.'
        IF ( mstout == 1 .AND. iceout == 1 ) THEN
          stgr = 4              ! S-points
          CALL intrp3d( nx,ny,nz,nx1,ny1,nz1,x,y,z,zp,x1,y1,z1,zp1,     &
                        samgrd,stgr,a3din,a3dout,tem3d,tem3d1 )
          WRITE(nchanl) 'qh      r3d1'
          WRITE(nchanl) a3dout
          CALL a3dmax0(a3dout,1,nx1,1,nx1-1,1,ny1,1,ny1-1,              &
                       1,nz1,1,nz1-1,amax,amin)
          WRITE(6,'(1x,2(a,e13.6))')                                    &
              'qhmin   = ', amin,',  qhmax   =',amax
          IF ( dmpexbc .AND. qhexout == 1 ) WRITE(exbchanl) a3dout
        END IF

      END IF

    END IF

    IF( tkein == 1 ) THEN

      READ(inch,ERR=9115,END=9125) label
      READ(inch,ERR=9115,END=9125) a3din
      WRITE(6,'(1x,a,a12,a)')                                           &
          'Field ',label,' was read into array a3din.'
      IF ( mstout == 1 .AND. iceout == 1 ) THEN
        stgr = 4              ! S-points
        CALL intrp3d( nx,ny,nz,nx1,ny1,nz1,x,y,z,zp,x1,y1,z1,zp1,       &
                      samgrd,stgr,a3din,a3dout,tem3d,tem3d1 )
        WRITE(nchanl) 'tke     r3d1'
        WRITE(nchanl) a3dout
        CALL a3dmax0(a3dout,1,nx1,1,nx1-1,1,ny1,1,ny1-1,                &
                     1,nz1,1,nz1-1,amax,amin)
        WRITE(6,'(1x,2(a,e13.6))')                                      &
            'tkemin  = ', amin,',  tkemax  =',amax
      END IF

    END IF

    IF( trbin == 1 ) THEN

      READ(inch,ERR=9115,END=9125) label
      READ(inch,ERR=9115,END=9125) a3din
      WRITE(6,'(1x,a,a12,a)')                                           &
          'Field ',label,' was read into array a3din.'
      IF ( mstout == 1 .AND. iceout == 1 ) THEN
        stgr = 4              ! S-points
        CALL intrp3d( nx,ny,nz,nx1,ny1,nz1,x,y,z,zp,x1,y1,z1,zp1,       &
                      samgrd,stgr,a3din,a3dout,tem3d,tem3d1 )
        WRITE(nchanl) 'kmh     r3d1'
        WRITE(nchanl) a3dout
        CALL a3dmax0(a3dout,1,nx1,1,nx1-1,1,ny1,1,ny1-1,                &
                     1,nz1,1,nz1-1,amax,amin)
        WRITE(6,'(1x,2(a,e13.6))')                                      &
            'kmhmin  = ', amin,',  kmhmax  =',amax
      END IF

      IF ( oldver == 0 ) THEN
        READ(inch,ERR=9115,END=9125) label
        READ(inch,ERR=9115,END=9125) a3din
        WRITE(6,'(1x,a,a12,a)')                                         &
            'Field ',label,' was read into array a3din.'
        IF ( mstout == 1 .AND. iceout == 1 ) THEN
          stgr = 4              ! S-points
          CALL intrp3d( nx,ny,nz,nx1,ny1,nz1,x,y,z,zp,x1,y1,z1,zp1,     &
                      samgrd,stgr,a3din,a3dout,tem3d,tem3d1 )
          WRITE(nchanl) 'kmv     r3d1'
          WRITE(nchanl) a3dout
          CALL a3dmax0(a3dout,1,nx1,1,nx1-1,1,ny1,1,ny1-1,              &
                     1,nz1,1,nz1-1,amax,amin)
          WRITE(6,'(1x,2(a,e13.6))')                                    &
              'kmvmin  = ', amin,',  kmvmax  =',amax
        END IF
      END IF

    END IF

    IF( sfcin == 1 ) THEN

      IF (nstyp == 1) THEN

!        READ(inch,ERR=9115,END=9125) label
!        READ(inch,ERR=9115,END=9125) ((a2din(i,j),i=1,nx),j=1,ny)
!        WRITE(6,'(1x,a,a12,a)')                                         &
!            'Field ',label,' was read into array a2din.'
!        CALL intrp2d(nx,ny,nx1,ny1,x,y,x1,y1,samgrd,a2din,a2dout)
!        WRITE(nchanl) 'tsfc    i2d0'
!        WRITE(nchanl) ((a2dout(i,j),i=1,nx1),j=1,ny1)
!        IF ( lsoildmp ) WRITE(soilchanl) a2dout
!
!        READ(inch,ERR=9115,END=9125) label
!        READ(inch,ERR=9115,END=9125) ((a2din(i,j),i=1,nx),j=1,ny)
!        WRITE(6,'(1x,a,a12,a)')                                         &
!            'Field ',label,' was read into array a2din.'
!        CALL intrp2d(nx,ny,nx1,ny1,x,y,x1,y1,samgrd,a2din,a2dout)
!        WRITE(nchanl) 'tsoil   i2d0'
!        WRITE(nchanl) ((a2dout(i,j),i=1,nx1),j=1,ny1)
!       IF ( lsoildmp ) WRITE(soilchanl) a2dout
!
!        READ(inch,ERR=9115,END=9125) label
!        READ(inch,ERR=9115,END=9125) ((a2din(i,j),i=1,nx),j=1,ny)
!        WRITE(6,'(1x,a,a12,a)')                                         &
!            'Field ',label,' was read into array a2din.'
!        CALL intrp2d(nx,ny,nx1,ny1,x,y,x1,y1,samgrd,a2din,a2dout)
!        WRITE(nchanl) 'wetsfc  i2d0'
!        WRITE(nchanl) ((a2dout(i,j),i=1,nx1),j=1,ny1)
!        IF ( lsoildmp ) WRITE(soilchanl) a2dout
!
!        READ(inch,ERR=9115,END=9125) label
!        READ(inch,ERR=9115,END=9125) ((a2din(i,j),i=1,nx),j=1,ny)
!        WRITE(6,'(1x,a,a12,a)')                                         &
!            'Field ',label,' was read into array a2din.'
!        CALL intrp2d(nx,ny,nx1,ny1,x,y,x1,y1,samgrd,a2din,a2dout)
!        WRITE(nchanl) 'wetdp   i2d0'
!        WRITE(nchanl) ((a2dout(i,j),i=1,nx1),j=1,ny1)
!        IF ( lsoildmp ) WRITE(soilchanl) a2dout
!
!        READ(inch,ERR=9115,END=9125) label
!        READ(inch,ERR=9115,END=9125) ((a2din(i,j),i=1,nx),j=1,ny)
!        WRITE(6,'(1x,a,a12,a)')                                         &
!            'Field ',label,' was read into array a2din.'
!        CALL intrp2d(nx,ny,nx1,ny1,x,y,x1,y1,samgrd,a2din,a2dout)
!        WRITE(nchanl) 'wetcanp i2d0'
!        WRITE(nchanl) ((a2dout(i,j),i=1,nx1),j=1,ny1)
!        IF ( lsoildmp ) WRITE(soilchanl) a2dout

        READ(inch,ERR=9115,END=9125) label
        READ(inch,ERR=9115,END=9125)                                    &
           (((tsoilin(i,j,k,0),i=1,nx),j=1,ny),k=1,nzsoil)
        WRITE(nchanl) 'tsoil  rs3d0'
        DO k = 1,nzsoil
          a2din(:,:) = tsoilin(:,:,k,0)
          CALL intrp2d(nx,ny,nx1,ny1,x,y,x1,y1,samgrd,a2din,a2dout)
          WRITE(nchanl) ((a2dout(i,j),i=1,nx1),j=1,ny1)
          IF ( lsoildmp ) WRITE(soilchanl) a2dout
        END DO

        READ(inch,ERR=9115,END=9125) label
        READ(inch,ERR=9115,END=9125)                                    &
           (((qsoilin(i,j,k,0),i=1,nx),j=1,ny),k=1,nzsoil)
        WRITE(nchanl) 'qsoil  rs3d0'
        DO k = 1,nzsoil
          a2din(:,:) = qsoilin(:,:,k,0)
          CALL intrp2d(nx,ny,nx1,ny1,x,y,x1,y1,samgrd,a2din,a2dout)
          WRITE(nchanl) ((a2dout(i,j),i=1,nx1),j=1,ny1)
          IF ( lsoildmp ) WRITE(soilchanl) a2dout
        END DO

        READ(inch,ERR=9115,END=9125) label
        READ(inch,ERR=9115,END=9125) ((a2din(i,j),i=1,nx),j=1,ny)
        WRITE(6,'(1x,a,a12,a)')                                         &
            'Field ',label,' was read into array a2din.'
        CALL intrp2d(nx,ny,nx1,ny1,x,y,x1,y1,samgrd,a2din,a2dout)
        WRITE(nchanl) 'wetcanp i2d0'
        WRITE(nchanl) ((a2dout(i,j),i=1,nx1),j=1,ny1)
        IF ( lsoildmp ) WRITE(soilchanl) a2dout

      ELSE

        DO is=0,nstyp
!          READ(inch,ERR=9115,END=9125) label
!          READ(inch,ERR=9115,END=9125) ((a2din(i,j),i=1,nx),j=1,ny)
!          WRITE(6,'(1x,a,a12,a,i2)')                                    &
!              'Field ',label,' was read into array a2din for is=',is
!          CALL intrp2d(nx,ny,nx1,ny1,x,y,x1,y1, samgrd,a2din,a2dout)
!          WRITE(nchanl) 'tsfc    i2d0'
!          WRITE(nchanl) ((a2dout(i,j),i=1,nx1),j=1,ny1)
!          IF ( lsoildmp .AND. is == 0 ) WRITE(soilchanl) a2dout
!
!          READ(inch,ERR=9115,END=9125) label
!          READ(inch,ERR=9115,END=9125) ((a2din(i,j),i=1,nx),j=1,ny)
!          WRITE(6,'(1x,a,a12,a,i2)')                                    &
!              'Field ',label,' was read into array a2din for is=',is
!          CALL intrp2d(nx,ny,nx1,ny1,x,y,x1,y1, samgrd,a2din,a2dout)
!          WRITE(nchanl) 'tsoil   i2d0'
!          WRITE(nchanl) ((a2dout(i,j),i=1,nx1),j=1,ny1)
!          IF ( lsoildmp .AND. is == 0 ) WRITE(soilchanl) a2dout
!
!          READ(inch,ERR=9115,END=9125) label
!          READ(inch,ERR=9115,END=9125) ((a2din(i,j),i=1,nx),j=1,ny)
!          WRITE(6,'(1x,a,a12,a,i2)')                                    &
!              'Field ',label,' was read into array a2din for is=',is
!          CALL intrp2d(nx,ny,nx1,ny1,x,y,x1,y1, samgrd,a2din,a2dout)
!          WRITE(nchanl) 'wetsfc  i2d0'
!          WRITE(nchanl) ((a2dout(i,j),i=1,nx1),j=1,ny1)
!          IF ( lsoildmp .AND. is == 0 ) WRITE(soilchanl) a2dout
!
!          READ(inch,ERR=9115,END=9125) label
!          READ(inch,ERR=9115,END=9125) ((a2din(i,j),i=1,nx),j=1,ny)
!          WRITE(6,'(1x,a,a12,a,i2)')                                    &
!              'Field ',label,' was read into array a2din for is=',is
!          CALL intrp2d(nx,ny,nx1,ny1,x,y,x1,y1, samgrd,a2din,a2dout)
!          WRITE(nchanl) 'wetdp   i2d0'
!          WRITE(nchanl) ((a2dout(i,j),i=1,nx1),j=1,ny1)
!          IF ( lsoildmp .AND. is == 0 ) WRITE(soilchanl) a2dout
!
!          READ(inch,ERR=9115,END=9125) label
!          READ(inch,ERR=9115,END=9125) ((a2din(i,j),i=1,nx),j=1,ny)
!          WRITE(6,'(1x,a,a12,a,i2)')                                    &
!              'Field ',label,' was read into array a2din for is=',is
!          CALL intrp2d(nx,ny,nx1,ny1,x,y,x1,y1, samgrd,a2din,a2dout)
!          WRITE(nchanl) 'wetcanp i2d0'
!          WRITE(nchanl) ((a2dout(i,j),i=1,nx1),j=1,ny1)
!          IF ( lsoildmp .AND. is == 0 ) WRITE(soilchanl) a2dout

          IF (is <= nstyps) THEN
            READ(inch,ERR=9115,END=9125) label
            READ(inch,ERR=9115,END=9125)                                &
                (((tsoilin(i,j,k,is),i=1,nx),j=1,ny),k=1,nzsoil)
            WRITE(nchanl) 'tsoil  rs3d0'
            DO k = 1,nzsoil
              a2din(:,:) = tsoilin(:,:,k,is)
              CALL intrp2d(nx,ny,nx1,ny1,x,y,x1,y1,samgrd,a2din,a2dout)
              WRITE(nchanl) ((a2dout(i,j),i=1,nx1),j=1,ny1)
              IF ( lsoildmp ) WRITE(soilchanl) a2dout
            END DO


            READ(inch,ERR=9115,END=9125) label
            READ(inch,ERR=9115,END=9125)                                &
                (((qsoilin(i,j,k,is),i=1,nx),j=1,ny),k=1,nzsoil)
            WRITE(nchanl) 'qsoil  rs3d0'
            DO k = 1,nzsoil
              a2din(:,:) = qsoilin(:,:,k,is)
              CALL intrp2d(nx,ny,nx1,ny1,x,y,x1,y1,samgrd,a2din,a2dout)
              WRITE(nchanl) ((a2dout(i,j),i=1,nx1),j=1,ny1)
              IF ( lsoildmp ) WRITE(soilchanl) a2dout
            END DO

            READ(inch,ERR=9115,END=9125) label
            READ(inch,ERR=9115,END=9125) ((a2din(i,j),i=1,nx),j=1,ny)
            WRITE(6,'(1x,a,a12,a,i2)')                                 &
              'Field ',label,' was read into array a2din for is=',is
            CALL intrp2d(nx,ny,nx1,ny1,x,y,x1,y1, samgrd,a2din,a2dout)
            WRITE(nchanl) 'wetcanp i2d0'
            WRITE(nchanl) ((a2dout(i,j),i=1,nx1),j=1,ny1)
            IF ( lsoildmp .AND. is == 0 ) WRITE(soilchanl) a2dout

          ENDIF

        END DO

      END IF
    END IF

    CLOSE ( inch )
    CLOSE ( nchanl )
    IF ( dmpexbc ) CLOSE ( exbchanl )
    IF ( lsoildmp ) CLOSE ( soilchanl )

    CALL retunit ( inch )
    CALL retunit ( nchanl )
    IF ( dmpexbc ) CALL retunit ( exbchanl )
    IF ( lsoildmp ) CALL retunit ( soilchanl )

  END DO

  STOP

  9000  CONTINUE

  WRITE(6,'(1x,a,i2,/1x,a)')                                            &
      'Namelist read error. Job stopped in ARPSINTRP.'
  STOP 9000
!
!-----------------------------------------------------------------------
!
!  Error during read
!
!----------------------------------------------------------------------
!

  9001  CONTINUE

  WRITE(6,'(/a,a/)')                                                    &
      ' Error in openning file ',grdbasfn(1:lengbf)
  STOP 9001

  9002  CONTINUE

  WRITE(6,'(/a,a/)')                                                    &
      ' Error in openning file ',basdmpfn(1:lbasdmpf)
  STOP 9002

  9003  CONTINUE

  WRITE(6,'(/a,a/)')                                                    &
      ' Error in openning file ',filename(1:lenfil)
  STOP 9003

  9004  CONTINUE

  WRITE(6,'(/a,a/)')                                                    &
      ' Error in openning file ',hdmpfn(1:ldmpf)
  STOP 9004

  9005  CONTINUE

  WRITE(6,'(/a,a/)')                                                    &
      ' Error in openning file ',exbcfn(1:length)
  STOP 9005

  9110  CONTINUE

  WRITE(6,'(/a,a/)') ' Error in reading file ',basdmpfn(1:lbasdmpf)

  STOP 9110

  9115  CONTINUE

  WRITE(6,'(/a,a/)') ' Error in reading file ',filename(1:lenfil)

  STOP 9115
!
!-----------------------------------------------------------------------
!
!  End-of-file during read
!
!----------------------------------------------------------------------
!

  9120  CONTINUE

  WRITE(6,'(/a,a/)')                                                    &
      ' Unexpected End of File reached in reading file ',               &
      basdmpfn(1:lbasdmpf)

  STOP 9120

  9125  CONTINUE

  WRITE(6,'(/a,a/)')                                                    &
      ' Unexpected End of File reached in reading file ',               &
      filename(1:lenfil)

  STOP 9125

END PROGRAM arpsintrp_ls
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE INTRP3D                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE intrp3d( nx,ny,nz,nx1,ny1,nz1,x,y,z,zp,x1,y1,z1,zp1,         &
           samgrd,stgr,a3din,a3dout,tem3d,tem3d1 )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Interpolate a 3-d array from an ARPS grid into a new ARPS grid.
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  08/05/1997
!
!  MODIFICATION HISTORY:
!
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Parameters for the utput grid.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz
  INTEGER :: nx1,ny1,nz1

  REAL :: x(nx)
  REAL :: y(ny)
  REAL :: z(nz)
  REAL :: zp(nx,ny,nz)

  REAL :: x1(nx1)
  REAL :: y1(ny1)
  REAL :: z1(nz1)
  REAL :: zp1(nx1,ny1,nz1)

  REAL :: a3din(nx,ny,nz)
  REAL :: a3dout(nx1,ny1,nz1)

  REAL :: tem3d(nx,ny,nz)
  REAL :: tem3d1(nx1,ny1,nz1)

  INTEGER :: samgrd,stgr
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k
  REAL :: s1,s2,s3,s4,sgrdinv
  REAL :: vgridinv,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8
  REAL :: xu1,yu1,zu1,xv1,yv1,zv1,xw1,yw1,zw1,xs1,ys1,zs1
  INTEGER :: iu,ju,ku,iv,jv,kv,iw,jw,kw,is,js,ks
  INTEGER :: kk
  REAL :: zupper,zlower
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'bndry.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF ( samgrd == 1 ) THEN

    DO k=1,nz1
      DO j=1,ny1
        DO i=1,nx1
          a3dout(i,j,k) = a3din(i,j,k)
        END DO
      END DO
    END DO

    RETURN

  END IF

  IF ( stgr == 1 ) THEN   ! for u-points

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=2,nx-1
          tem3d(i,j,k) = 0.25*(zp(i  ,j,k)+zp(i  ,j,k+1)                &
                              +zp(i-1,j,k)+zp(i-1,j,k+1))
        END DO
      END DO
    END DO

    CALL bcsu(nx,ny,nz,1,ny,1,nz,ebc,wbc,tem3d)

    DO k=1,nz1-1
      DO j=1,ny1-1
        DO i=2,nx1-1
          tem3d1(i,j,k) = 0.25*(zp1(i  ,j,k)+zp1(i  ,j,k+1)             &
                               +zp1(i-1,j,k)+zp1(i-1,j,k+1))
        END DO
      END DO
    END DO

    CALL bcsu(nx1,ny1,nz1,1,ny1,1,nz1,ebc,wbc,tem3d1)

    DO i=1,nx1
      DO j=1,ny1-1
        DO k=1,nz1-1
          xu1=x1(i)
          yu1=(y1(j+1)+y1(j))*0.5
          zu1= tem3d1(i,j,k)

          iu = MAX(1, MIN(nx-1, INT((xu1-x(1))/dx)+1 ))
          ju = MAX(1, MIN(ny-2, INT((yu1-(y(1)+y(2))*0.5)/dy)+1 ))

          s1 = ABS((xu1-x(iu  ))*(yu1-(y(ju  )+y(ju+1))*0.5))
          s2 = ABS((xu1-x(iu+1))*(yu1-(y(ju  )+y(ju+1))*0.5))
          s3 = ABS((xu1-x(iu+1))*(yu1-(y(ju+1)+y(ju+2))*0.5))
          s4 = ABS((xu1-x(iu  ))*(yu1-(y(ju+1)+y(ju+2))*0.5))

          sgrdinv = 1.0/(s1+s2+s3+s4)

          ku = 1
          DO kk=nz-2,1,-1
            IF(zu1 >= (tem3d(iu  ,ju  ,kk)*s3+tem3d(iu+1,ju  ,kk)*s4    &
                        +tem3d(iu+1,ju+1,kk)*s1+tem3d(iu  ,ju+1,kk)*s2) &
                        *sgrdinv) THEN
              ku = kk
              EXIT
            END IF
          END DO
          130       CONTINUE

          zlower = (tem3d(iu  ,ju  ,ku)*s3+tem3d(iu+1,ju  ,ku)*s4       &
                   +tem3d(iu+1,ju+1,ku)*s1+tem3d(iu  ,ju+1,ku)*s2)      &
                   *sgrdinv
          zupper = (tem3d(iu  ,ju  ,ku+1)*s3+tem3d(iu+1,ju  ,ku+1)*s4   &
                   +tem3d(iu+1,ju+1,ku+1)*s1+tem3d(iu  ,ju+1,ku+1)*s2)  &
                   *sgrdinv

          tmp1 = ABS(s1*(zu1-zlower))
          tmp2 = ABS(s2*(zu1-zlower))
          tmp3 = ABS(s3*(zu1-zlower))
          tmp4 = ABS(s4*(zu1-zlower))
          tmp5 = ABS(s1*(zu1-zupper))
          tmp6 = ABS(s2*(zu1-zupper))
          tmp7 = ABS(s3*(zu1-zupper))
          tmp8 = ABS(s4*(zu1-zupper))

          vgridinv=1.0/(tmp1+tmp2+tmp3+tmp4+tmp5+tmp6+tmp7+tmp8)

          a3dout(i,j,k) =                                               &
              (a3din(iu  ,ju  ,ku  )*tmp7+a3din(iu+1,ju  ,ku  )*tmp8    &
              +a3din(iu+1,ju+1,ku  )*tmp5+a3din(iu  ,ju+1,ku  )*tmp6    &
              +a3din(iu  ,ju  ,ku+1)*tmp3+a3din(iu+1,ju  ,ku+1)*tmp4    &
              +a3din(iu+1,ju+1,ku+1)*tmp1+a3din(iu  ,ju+1,ku+1)*tmp2)   &
              *vgridinv

        END DO
      END DO
    END DO

  ELSE IF ( stgr == 2 ) THEN

    DO k=1,nz-1
      DO j=2,ny-1
        DO i=1,nx-1
          tem3d(i,j,k) = 0.25*(zp(i,j  ,k)+zp(i,j  ,k+1)                &
                              +zp(i,j-1,k)+zp(i,j-1,k+1))
        END DO
      END DO
    END DO

    CALL bcsu(nx,ny,nz,1,nx,1,nz,nbc,sbc,tem3d)

    DO k=1,nz1-1
      DO j=2,ny1-1
        DO i=1,nx1-1
          tem3d1(i,j,k) = 0.25*(zp1(i,j  ,k)+zp1(i,j  ,k+1)             &
                               +zp1(i,j-1,k)+zp1(i,j-1,k+1))
        END DO
      END DO
    END DO

    CALL bcsv(nx1,ny1,nz1,1,nx1,1,nz1,nbc,sbc,tem3d1)

    DO i=1,nx1-1
      DO j=1,ny1
        DO k=1,nz1-1

          xv1=(x1(i+1)+x1(i))*0.5
          yv1= y1(j)
          zv1= tem3d1(i,j,k)

          iv = MAX(1, MIN(nx-2, INT((xv1-(x(1)+x(2))*0.5)/dx)+1 ))
          jv = MAX(1, MIN(ny-1, INT((yv1- y(1))/dy)+1 ))

          s1 = ABS((xv1-(x(iv  )+x(iv+1))*0.5)*(yv1-y(jv  )))
          s2 = ABS((xv1-(x(iv+1)+x(iv+2))*0.5)*(yv1-y(jv  )))
          s3 = ABS((xv1-(x(iv+1)+x(iv+2))*0.5)*(yv1-y(jv+1)))
          s4 = ABS((xv1-(x(iv  )+x(iv+1))*0.5)*(yv1-y(jv+1)))

          sgrdinv = 1.0/(s1+s2+s3+s4)

          kv = 1
          DO kk=nz-2,1,-1
            IF(zv1 >= (tem3d(iv  ,jv  ,kk)*s3+tem3d(iv+1,jv  ,kk)*s4    &
                        +tem3d(iv+1,jv+1,kk)*s1+tem3d(iv  ,jv+1,kk)*s2) &
                        *sgrdinv) THEN
              kv = kk
              EXIT
            END IF
          END DO
          230       CONTINUE

          zlower = (tem3d(iv  ,jv  ,kv)*s3+tem3d(iv+1,jv  ,kv)*s4       &
                   +tem3d(iv+1,jv+1,kv)*s1+tem3d(iv  ,jv+1,kv)*s2)      &
                   *sgrdinv
          zupper = (tem3d(iv  ,jv  ,kv+1)*s3+tem3d(iv+1,jv  ,kv+1)*s4   &
                   +tem3d(iv+1,jv+1,kv+1)*s1+tem3d(iv  ,jv+1,kv+1)*s2)  &
                   *sgrdinv

          tmp1 = ABS(s1*(zv1-zlower))
          tmp2 = ABS(s2*(zv1-zlower))
          tmp3 = ABS(s3*(zv1-zlower))
          tmp4 = ABS(s4*(zv1-zlower))
          tmp5 = ABS(s1*(zv1-zupper))
          tmp6 = ABS(s2*(zv1-zupper))
          tmp7 = ABS(s3*(zv1-zupper))
          tmp8 = ABS(s4*(zv1-zupper))

          vgridinv=1.0/(tmp1+tmp2+tmp3+tmp4+tmp5+tmp6+tmp7+tmp8)

          a3dout(i,j,k) =                                               &
              (a3din(iv  ,jv  ,kv  )*tmp7+a3din(iv+1,jv  ,kv  )*tmp8    &
              +a3din(iv+1,jv+1,kv  )*tmp5+a3din(iv  ,jv+1,kv  )*tmp6    &
              +a3din(iv  ,jv  ,kv+1)*tmp3+a3din(iv+1,jv  ,kv+1)*tmp4    &
              +a3din(iv+1,jv+1,kv+1)*tmp1+a3din(iv  ,jv+1,kv+1)*tmp2)   &
              *vgridinv

        END DO
      END DO
    END DO

  ELSE IF ( stgr == 3 ) THEN

    DO i=1,nx1-1
      DO j=1,ny1-1
        DO k=1,nz1

          xw1= (x1(i)+x1(i+1))*0.5
          yw1= (y1(j)+y1(j+1))*0.5
          zw1= zp1(i,j,k)

          iw = MAX(1, MIN(nx-2, INT((xw1-(x(1)+x(2))*0.5)/dx)+1 ))
          jw = MAX(1, MIN(ny-2, INT((yw1-(y(1)+y(2))*0.5)/dy)+1 ))

          s1=ABS((xw1-(x(iw  )+x(iw+1))*0.5)                            &
                *(yw1-(y(jw  )+y(jw+1))*0.5))
          s2=ABS((xw1-(x(iw+1)+x(iw+2))*0.5)                            &
                *(yw1-(y(jw  )+y(jw+1))*0.5))
          s3=ABS((xw1-(x(iw+1)+x(iw+2))*0.5)                            &
                *(yw1-(y(jw+1)+y(jw+2))*0.5))
          s4=ABS((xw1-(x(iw  )+x(iw+1))*0.5)                            &
                *(yw1-(y(jw+1)+y(jw+2))*0.5))

          sgrdinv = 1.0/(s1+s2+s3+s4)

          kw = 1
          DO kk=nz-1,1,-1
            IF(zw1 >= (zp(iw  ,jw  ,kk)*s3+zp(iw+1,jw  ,kk)*s4          &
                      +zp(iw+1,jw+1,kk)*s1+zp(iw  ,jw+1,kk)*s2)         &
                      *sgrdinv) THEN
              kw = kk
              EXIT
            END IF
          END DO
          320       CONTINUE

          zlower = (zp(iw  ,jw  ,kw)*s3+zp(iw+1,jw  ,kw)*s4             &
                 +zp(iw+1,jw+1,kw)*s1+zp(iw  ,jw+1,kw)*s2)              &
                 *sgrdinv
          zupper = (zp(iw  ,jw  ,kw+1)*s3+zp(iw+1,jw  ,kw+1)*s4         &
                 +zp(iw+1,jw+1,kw+1)*s1+zp(iw  ,jw+1,kw+1)*s2)          &
                 *sgrdinv

          tmp1 = ABS(s1*(zw1-zlower))
          tmp2 = ABS(s2*(zw1-zlower))
          tmp3 = ABS(s3*(zw1-zlower))
          tmp4 = ABS(s4*(zw1-zlower))
          tmp5 = ABS(s1*(zw1-zupper))
          tmp6 = ABS(s2*(zw1-zupper))
          tmp7 = ABS(s3*(zw1-zupper))
          tmp8 = ABS(s4*(zw1-zupper))

          vgridinv = 1.0/(tmp1+tmp2+tmp3+tmp4+tmp5+tmp6+tmp7+tmp8)

          a3dout(i,j,k) =                                               &
              (a3din(iw  ,jw  ,kw  )*tmp7+a3din(iw+1,jw  ,kw  )*tmp8    &
              +a3din(iw+1,jw+1,kw  )*tmp5+a3din(iw  ,jw+1,kw  )*tmp6    &
              +a3din(iw  ,jw  ,kw+1)*tmp3+a3din(iw+1,jw  ,kw+1)*tmp4    &
              +a3din(iw+1,jw+1,kw+1)*tmp1+a3din(iw  ,jw+1,kw+1)*tmp2)   &
              *vgridinv

        END DO
      END DO
    END DO

  ELSE IF ( stgr == 4 ) THEN

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem3d(i,j,k) = 0.5*(zp(i,j,k)+zp(i,j,k+1))
        END DO
      END DO
    END DO

    DO k=1,nz1-1
      DO j=1,ny1-1
        DO i=1,nx1-1
          tem3d1(i,j,k) = 0.5*(zp1(i,j,k)+zp1(i,j,k+1))
        END DO
      END DO
    END DO

    DO i=1,nx1-1
      DO j=1,ny1-1
        DO k=1,nz1-1

          xs1= (x1(i)+x1(i+1))*0.5
          ys1= (y1(j)+y1(j+1))*0.5
          zs1= tem3d1(i,j,k)

          is = MAX(1, MIN(nx-2, INT((xs1-(x(1)+x(2))*0.5)/dx)+1 ))
          js = MAX(1, MIN(ny-2, INT((ys1-(y(1)+y(2))*0.5)/dy)+1 ))

          s1 =ABS((xs1-(x(is  )+x(is+1))*0.5)                           &
                 *(ys1-(y(js  )+y(js+1))*0.5))
          s2 =ABS((xs1-(x(is+1)+x(is+2))*0.5)                           &
                 *(ys1-(y(js  )+y(js+1))*0.5))
          s3 =ABS((xs1-(x(is+1)+x(is+2))*0.5)                           &
                 *(ys1-(y(js+1)+y(js+2))*0.5))
          s4 =ABS((xs1-(x(is  )+x(is+1))*0.5)                           &
                 *(ys1-(y(js+1)+y(js+2))*0.5))

          sgrdinv = 1.0/(s1+s2+s3+s4)

          ks = 1
          DO kk=nz-2,1,-1
            IF(zs1 >= (tem3d(is  ,js  ,kk)*s3+tem3d(is+1,js  ,kk)*s4    &
                      +tem3d(is+1,js+1,kk)*s1+tem3d(is  ,js+1,kk)*s2)   &
                      *sgrdinv ) THEN
              ks = kk
              EXIT
            END IF
          END DO
          430       CONTINUE

          zlower = (tem3d(is  ,js  ,ks)*s3+tem3d(is+1,js  ,ks)*s4       &
                 +tem3d(is+1,js+1,ks)*s1+tem3d(is  ,js+1,ks)*s2)        &
                 *sgrdinv
          zupper = (tem3d(is  ,js  ,ks+1)*s3+tem3d(is+1,js  ,ks+1)*s4   &
                 +tem3d(is+1,js+1,ks+1)*s1+tem3d(is  ,js+1,ks+1)*s2)    &
                 *sgrdinv

          tmp1 = ABS(s1*(zs1-zlower))
          tmp2 = ABS(s2*(zs1-zlower))
          tmp3 = ABS(s3*(zs1-zlower))
          tmp4 = ABS(s4*(zs1-zlower))
          tmp5 = ABS(s1*(zs1-zupper))
          tmp6 = ABS(s2*(zs1-zupper))
          tmp7 = ABS(s3*(zs1-zupper))
          tmp8 = ABS(s4*(zs1-zupper))

          vgridinv=1.0/(tmp1+tmp2+tmp3+tmp4+tmp5+tmp6+tmp7+tmp8)

          a3dout(i,j,k) =                                               &
              (a3din(is  ,js  ,ks  )*tmp7+a3din(is+1,js  ,ks  )*tmp8    &
              +a3din(is+1,js+1,ks  )*tmp5+a3din(is  ,js+1,ks  )*tmp6    &
              +a3din(is  ,js  ,ks+1)*tmp3+a3din(is+1,js  ,ks+1)*tmp4    &
              +a3din(is+1,js+1,ks+1)*tmp1+a3din(is  ,js+1,ks+1)*tmp2)   &
              *vgridinv

        END DO
      END DO
    END DO

  END IF

  RETURN
END SUBROUTINE intrp3d
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE INTRP2D                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE intrp2d( nx,ny,nx1,ny1,x,y,x1,y1,                            &
           samgrd,a2din,a2dout )
!
!-----------------------------------------------------------------------
!
!  Interpolate a 2-d array from an ARPS grid into a new ARPS grid.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  08/05/1997
!
!  MODIFICATION HISTORY:
!
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Parameters for the utput grid.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny
  INTEGER :: nx1,ny1

  REAL :: x(nx)
  REAL :: y(ny)

  REAL :: x1(nx1)
  REAL :: y1(ny1)

  REAL :: a2din(nx,ny)
  REAL :: a2dout(nx1,ny1)

  INTEGER :: samgrd
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j
  REAL :: s1,s2,s3,s4,sgrdinv
  REAL :: xs1,ys1
  INTEGER :: is,js
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF ( samgrd == 1 ) THEN
    DO j=1,ny1
      DO i=1,nx1
        a2dout(i,j) = a2din(i,j)
      END DO
    END DO

    RETURN

  END IF

  DO i=1,nx1-1
    DO j=1,ny1-1

      xs1= (x1(i)+x1(i+1))*0.5
      ys1= (y1(j)+y1(j+1))*0.5

      is = MAX(1, MIN(nx-2, INT((xs1-(x(1)+x(2))*0.5)/dx)+1 ))
      js = MAX(1, MIN(ny-2, INT((ys1-(y(1)+y(2))*0.5)/dy)+1 ))

      s1=ABS((xs1-(x(is  )+x(is+1))*0.5)*(ys1-(y(js  )+y(js+1))*0.5))
      s2=ABS((xs1-(x(is+1)+x(is+2))*0.5)*(ys1-(y(js  )+y(js+1))*0.5))
      s3=ABS((xs1-(x(is+1)+x(is+2))*0.5)*(ys1-(y(js+1)+y(js+2))*0.5))
      s4=ABS((xs1-(x(is  )+x(is+1))*0.5)*(ys1-(y(js+1)+y(js+2))*0.5))

      sgrdinv = 1.0/(s1+s2+s3+s4)

      a2dout(i,j) =                                                     &
          (a2din(is  ,js  )*s3+a2din(is+1,js  )*s4                      &
          +a2din(is+1,js+1)*s1+a2din(is  ,js+1)*s2)                     &
          *sgrdinv

    END DO
  END DO

  RETURN
END SUBROUTINE intrp2d
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE DIST2D                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE dist2d( nx,ny,nx1,ny1,x,y,x1,y1,                             &
           samgrd,i2din,i2dout )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Re-distribute a 2-d integer array from an ARPS grid into another
!  grid.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  08/05/1997
!
!  MODIFICATION HISTORY:
!
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Parameters for the utput grid.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny
  INTEGER :: nx1,ny1

  REAL :: x(nx)
  REAL :: y(ny)

  REAL :: x1(nx1)
  REAL :: y1(ny1)

  INTEGER :: i2din(nx,ny)
  INTEGER :: i2dout(nx1,ny1)

  INTEGER :: samgrd
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j
  REAL :: xs1,ys1
  INTEGER :: is,js
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF ( samgrd == 1 ) THEN
    DO j=1,ny1
      DO i=1,nx1
        i2dout(i,j) = i2din(i,j)
      END DO
    END DO

    RETURN
  END IF

  DO i=1,nx1-1
    DO j=1,ny1-1

      xs1= (x1(i)+x1(i+1))*0.5
      ys1= (y1(j)+y1(j+1))*0.5

      is = MAX(1, MIN(nx-2, INT((xs1-(x(1)+x(2))*0.5)/dx)+1 ))
      js = MAX(1, MIN(ny-2, INT((ys1-(y(1)+y(2))*0.5)/dy)+1 ))

      i2dout(i,j) = i2din(is,js)
    END DO
  END DO

  RETURN
END SUBROUTINE dist2d
