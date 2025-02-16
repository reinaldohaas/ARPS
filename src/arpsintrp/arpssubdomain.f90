PROGRAM arpssubdomain
!
!##################################################################
!##################################################################
!######                                                      ######
!######                   PROGRAM ARPSSUBDOMAIN              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  (12/1/2002)
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'exbc.inc'
  INCLUDE 'mp.inc'

  REAL, allocatable :: x (:),y (:),z (:)
  REAL, allocatable :: x1(:),y1(:),z1(:)

  REAL, allocatable :: xs (:), ys (:) ! x,y coord for scalar points
  REAL, allocatable :: xs1(:), ys1(:) ! x,y coord for scalar points
  INTEGER, allocatable :: isx(:),jsy(:),iux(:),jvy(:)
  REAL, allocatable:: wgtsx(:),wgtsy(:),wgtux(:),wgtvy(:),temx1yz(:,:)

  REAL, allocatable    :: stypfrct(:,:,:)
  INTEGER, allocatable :: soiltyp(:,:,:)

  REAL, allocatable :: tem1(:,:,:)
  REAL, allocatable :: tem11(:,:,:)

  INTEGER(2), allocatable :: tem2(:,:,:)
  INTEGER(2), allocatable :: tem21(:,:,:)

  INTEGER :: nx,ny,nz,nzsoil   ! Grid dimensions.
  INTEGER :: nstyps            ! Maximum number of soil types.

  INTEGER :: hinfmt,nhisfile_max,nhisfile,lengbf,lenfil
  PARAMETER (nhisfile_max=1000)
  CHARACTER (LEN=256) :: grdbasfn,hisfile(nhisfile_max)

  INTEGER :: nxbgn, nxend, nybgn, nyend, nzbgn, nzend, intrp_opt
  INTEGER :: nx1,ny1,nz1,nzsoil1,nstyp1
  integer :: savespace
  REAL :: xctr1,yctr1,dx1,dy1,xorig1,yorig1

  NAMELIST /output_dims/ intrp_opt,   &
            nxbgn, nxend, nybgn, nyend, nzbgn, nzend, &
            nx1,ny1,xctr1,yctr1,dx1,dy1, savespace,  &
            dz,strhopt,dzmin,zrefsfc,dlayer1,dlayer2,strhtune,zflat

  REAL :: tintv_dmpin, tbgn_dmpin, tend_dmpin
  INTEGER :: nextrafiles,subsample_main,subsample_extra
  CHARACTER (LEN=6) :: extra_field(100)
  CHARACTER (LEN=256) :: dirname_extra_in, dirname_extra_out
  NAMELIST /extra_files/ subsample_main,subsample_extra,tintv_dmpin,tbgn_dmpin,tend_dmpin,extra_field,dirname_extra_in,dirname_extra_out

  CHARACTER (LEN=80) :: runname_input
  NAMELIST /jobname/ runname_input

  INTEGER :: istatus

  NAMELIST /output/ dirname,hdmpfmt,grbpkbit,hdfcompr,          &
            grdout,basout,varout,mstout,rainout,prcout,iceout,          &
            tkeout,trbout,sfcout,landout,totout,filcmprs,readyfl,  &
            exbcdmp,qcexout,qrexout,qiexout,qsexout,qgexout,qhexout,ngbrz,rayklow, &
            terndmp

  INTEGER :: i,j,k,kk,ireturn

  REAL :: time
  LOGICAL :: iexist

  INTEGER :: nfile

  CHARACTER (LEN=30) :: substring

  CHARACTER (LEN=40 ) :: varunits
  CHARACTER (LEN=40 ) :: varname
  CHARACTER (LEN=6 ) :: varid

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  WRITE(6,'(/9(/5x,a)/)')                                               &
     '###############################################################', &
     '###############################################################', &
     '###                                                         ###', &
     '###                Welcome to ARPSSUBDOMAIN                 ###', &
     '###      This program read, interpolate then write out      ###', &
     '###      variables in ARPS history for files one variable   ###', &
     '###      at a time to minimize memory requirement.          ###', &
     '###                                                         ###', &
     '###############################################################', &
     '###############################################################'

! DTD: Assume that this is a single processor job: calls to readvar
! and wrtvar still need the nproc_x and nproc_y variables to be set

  nproc_x = 1
  nproc_y = 1
  max_fopen = 1

!
!-----------------------------------------------------------------------
!  Get the names of the input data files.
!-----------------------------------------------------------------------
!
  hinfmt = 3

  CALL get_input_file_names(5,hinfmt,grdbasfn,hisfile,nhisfile)

  lengbf = len_trim(grdbasfn)
  WRITE(6,'(/a,a)')' The grid/base name is ', trim(grdbasfn)

  lenfil = len_trim(hisfile(1))

  CALL get_dims_from_data(hinfmt,hisfile(1)(1:lenfil),                    &
       nx,ny,nz,nzsoil,nstyps, ireturn)

  Print*,'nx,ny,nz of input data were ', nx,ny,nz

  nzsoil1 = nzsoil
  nstyp1 = nstyp

  ALLOCATE(x (nx ),  STAT=istatus)
  ALLOCATE(y (ny ),  STAT=istatus)
  ALLOCATE(xs(nx ),  STAT=istatus)
  ALLOCATE(ys(ny ),  STAT=istatus)
  ALLOCATE(z (nz ),  STAT=istatus)
  x = 0.0
  y = 0.0
  xs = 0.0
  ys = 0.0
  z = 0.0
!  allocate(zpsoil(nx,ny,nzsoil),stat=istatus)
!  zpsoil = 0.0


  CALL get_gridinfo_from_hdf(grdbasfn(1:lengbf),nx,ny,nz,x,y,z,ireturn)

! Print*,'x,y,z,zp of input data read in.'
! print*,'x(1 )=',x(1)
! print*,'x(nx)=',x(nx)
! print*,'y(1 )=',y(1)
! print*,'y(ny)=',y(ny)

  dx = x(2) - x(1)
  dy = y(2) - y(1)

  IF (nstyps <= 0) nstyps = 1
  nstyp = nstyps

  IF( ireturn /= 0 ) THEN
    PRINT*,'Problem occured when trying to get dimensions from data.'
    PRINT*,'Program stopped.'
    STOP
  END IF

  WRITE(6,'(4(a,i5))') 'nx =',nx,', ny=',ny,', nz=',nz , &
            'nzsoil=',nzsoil

  intrp_opt = 0
  savespace = 1

  READ(5,output_dims, END=100)
  WRITE(6,'(/a,a/)')                                                    &
      'NAMELIST block output_dims successfully read.'

  IF( strhopt == 0.AND.dzmin /= dz ) THEN
    WRITE(6,'(5x,a)')                                                   &
         'For non-stretched case, dzmin was reset to dz.'
    dzmin = dz
  END IF

  if( intrp_opt == 0 ) then
    nx1 = nxend - nxbgn +1
    ny1 = nyend - nybgn +1
  else
    nxbgn = 1
    nxend = nx1
    nybgn = 1
    nyend = ny1
  endif

  nz1 = nzend - nzbgn +1

  WRITE(6,'(4(a,i5))') 'nx1 =',nx1,', ny1=',ny1,', nz1=',nz1

  ALLOCATE(x1 (nx1), STAT=istatus)
  ALLOCATE(y1 (ny1), STAT=istatus)

  ALLOCATE(xs1(nx1), STAT=istatus)
  ALLOCATE(ys1(ny1), STAT=istatus)
  allocate(z1(nz1),stat=istatus)
  x1  = 0.0
  y1  = 0.0
  xs1 = 0.0
  ys1 = 0.0
  z1 = 0.0

!  allocate(zpsoil1(nx1,ny1,nzsoil),stat=istatus)
!  zpsoil1 = 0.0

  allocate(stypfrct(nx,ny,nstyps), stat=istatus)
  allocate(soiltyp (nx,ny,nstyps), stat=istatus)

  READ (5,jobname,END=100)
  WRITE(6,'(/a/)') 'Sucessfully read namelist block JOBNAME.'

  extra_field = ''

  READ (5,extra_files,END=100)
  WRITE(6,'(/a/)') 'Sucessfully read namelist block EXTRA_FILES.'


  WRITE(6,'(/2x,a,a)') 'The name of this run is: ', runname_input
!
!-----------------------------------------------------------------------
!  Set the control parameters for output:
!-----------------------------------------------------------------------
!
  READ (5,output,END=100)

  ldirnam=LEN_trim(dirname)
  print*,'ldirnam, dirname', ldirnam, trim(dirname)

  IF( dirname(1:ldirnam) /= ' ') THEN
    CALL inquiredir(dirname(1:ldirnam),iexist)

    IF( .NOT.iexist ) THEN
      WRITE(6,'(5x,a,2(/5x,a))')                                        &
          'Specified output directory '//dirname(1:ldirnam)//           &
          ' not found.','It will be created by the program.'
      CALL unixcmd( 'mkdir -p '//dirname(1:ldirnam) )
    END IF
  END IF

  WRITE(6,'(5x,a)')                                                   &
      'Output files will be in directory '//dirname(1:ldirnam)//'.'

  WRITE(6,'(/a/)') 'Output control parameters read in are:'
  WRITE(6,output)

  ALLOCATE(tem1(nx,ny,nz),stat=istatus)
  CALL alloc_status_accounting(istatus,nx*ny*nz,'tem1')
  tem1 = 0.0
  ALLOCATE(tem2(nx,ny,nz),stat=istatus)
  CALL alloc_status_accounting(istatus,nx*ny*nz,'tem2')
  tem2 = 0.0

  ALLOCATE(tem11(nx1,ny1,nz1),stat=istatus)
  CALL alloc_status_accounting(istatus,nx1*ny1*nz1,'tem11')
  tem11 = 0.0
  ALLOCATE(tem21(nx1,ny1,nz1),stat=istatus)
  CALL alloc_status_accounting(istatus,nx1*ny1*nz1,'tem21')
  tem21 = 0.0

  if( intrp_opt == 0 ) then
     DO i=1,nx1
       x1(i) = x(nxbgn+i-1)
     END DO
     xctr1 = 0.5*(x(nxbgn)+x(nxend))

     DO j=1,ny1
       y1(j) =y(nybgn+j-1)
     END DO
     yctr1 = 0.5*(y(nybgn)+y(nyend))

   else

     xorig1 = xctr1 - (nx1-3)*dx1*0.5
     yorig1 = yctr1 - (ny1-3)*dy1*0.5

     DO i=1,nx1
       x1(i) =xorig1+(i-2.0)*dx1
       xs1(i)=xorig1+(i-1.5)*dx1
     END DO

     DO j=1,ny1
       y1 (j) =yorig1+(j-2.0)*dy1
       ys1(j) =yorig1+(j-1.5)*dy1
     END DO

     DO i=1,nx-1
       xs(i)=0.5*(x(i)+x(i+1))
     END DO

     DO j=1,ny-1
       ys(j)=0.5*(y(j)+y(j+1))
     END DO

     print*,'xorig1,xs1(1),xs1(nx1-1)=',xorig1,xs1(1),xs1(nx1-1)
     print*,'yorig1,ys1(1),ys1(ny1-1)=',yorig1,ys1(1),ys1(ny1-1)

     allocate(isx(nx1),stat=istatus)
     isx = 0.0
     allocate(jsy(ny1),stat=istatus)
     jsy = 0.0
     allocate(iux(nx1),stat=istatus)
     iux = 0.0
     allocate(jvy(ny1),stat=istatus)
     jvy = 0.0
     allocate(wgtsx(nx1),stat=istatus)
     wgtsx = 0.0
     allocate(wgtsy(ny1),stat=istatus)
     wgtsy = 0.0
     allocate(wgtux(nx1),stat=istatus)
     wgtux = 0.0
     allocate(wgtvy(ny1),stat=istatus)
     wgtvy = 0.0

     allocate(temx1yz(nx1,ny),stat=istatus)
     temx1yz =0.0

     DO i=1,nx1-1
       isx(i) = MAX(1, MIN(nx-2, INT((xs1(i)-xs(1))/dx)+1 ))
       wgtsx(i)= (xs(isx(i)+1)-xs1(i))/(xs(isx(i)+1)-xs(isx(i)))
     END DO

     DO j=1,ny1-1
       jsy(j) = MAX(1, MIN(ny-2, INT((ys1(j)-ys(1))/dy)+1 ))
       wgtsy(j)= (ys(jsy(j)+1)-ys1(j))/(ys(jsy(j)+1)-ys(jsy(j)))
     END DO

     print*,'isx(1),jsy(1)=', isx(1),jsy(1)

     DO i=1,nx1
       iux(i) = MAX(1, MIN(nx-1, INT((x1 (i)-x (1))/dx)+1 ))
       wgtux(i)= (x (iux(i)+1)-x1 (i))/(x (iux(i)+1)-x (iux(i)))
     END DO

     DO j=1,ny1
       jvy(j) = MAX(1, MIN(ny-1, INT((y1 (j)-y (1))/dy)+1 ))
       wgtvy(j)= (y (jvy(j)+1)-y1 (j))/(y (jvy(j)+1)-y (jvy(j)))
     END DO

   end if

   print*,'nx1, ny1 =', nx1, ny1
   print*,'x1(1), x1(nx1), y1(1), y1(ny1) for output grid are :'
   print*,x1(1), x1(nx1), y1(1), y1(ny1)

   DO k=1,nz1
     z1(k)=z(nzbgn+k-1)
   END DO

   ! DTD: added condition for subsampling main history files

   IF(subsample_main == 1) THEN

     DO nfile = 1,nhisfile

      WRITE(6,'(/a,a,a)') ' Data set ', trim(hisfile(nfile)),' to be subsampled.'

      CALL dtareaddump(nx,ny,nz,nzsoil,x,y,z,nstyps,                    &
           trim(grdbasfn),trim(hisfile(nfile)),                         &
           nx1,ny1,nz1,x1,y1,z1,soiltyp,stypfrct,xctr1,yctr1,           &
           nxbgn,nxend,nybgn,nyend,nzbgn,nzend,                         &
           tem1,tem2,tem11,tem21,intrp_opt,savespace,                   &
           isx,jsy,iux,jvy,wgtsx,wgtsy,wgtux,wgtvy,temx1yz)

     END DO

   END IF

   ! DTD:  Added code below to subsample individual files (written by wrtvar subroutine)
   ! No interpolation is performed for now, so only works with the intrpopt == 0
   ! Assumes files were written in HDF for now

   IF(subsample_extra == 1) THEN

     nextrafiles = INT((tend_dmpin-tbgn_dmpin)/tintv_dmpin+1)
     time=tbgn_dmpin

     DO nfile=1,nextrafiles
       IF(nfile == 1) THEN ! Construct the trailer string

         if( nxbgn>999.or.nxend>999.or.nybgn>999.or.nyend>999 &
                               .or.nzbgn>999.or.nzend>999) then
           write(substring,'(6(a1,i4.4))') &
                 '.',nxend,'x',nyend,'x',nzend
         else
           write(substring,'(6(a1,i3.3))') &
                 '.',nxend,'x',nyend,'x',nzend
         endif
       END IF

       ! For each field requested, read in the appropriate file and
       ! subsample the field, then write out to the new file
       tem1 = 0.0
       tem11 = 0.0
       k=1
       DO WHILE(extra_field(k) /= '')
         varid=extra_field(k)
         CALL readvar2(nx,ny,nz,tem1,varid,varname,varunits,time,runname_input,  &
                       dirname_extra_in,3,0,istatus)

         print*,'max val = ',maxval(tem1)
         call copyarray(tem1,nx,ny,nz,tem11, nx1,ny1,nz1, &
              nxbgn,nxend,nybgn,nyend,nzbgn,nzend)
         print*,'max val copied = ',maxval(tem11)
         CALL wrtvar3(nx1,ny1,nz1,tem11,varid,varname,varunits,time,runname_input,  &
                     dirname_extra_out,substring,3,2,0,istatus)
         k=k+1
       END DO

       time=time+tintv_dmpin

     END DO

   END IF

  STOP

  100 WRITE(6,'(a)') 'Error reading NAMELIST file. Program ARPSSUBDOMAIN stopped.'
  STOP

END PROGRAM arpssubdomain
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE DTAREADDUMP                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE dtareaddump(nx,ny,nz,nzsoil,x,y,z,nstyps,    &
           grdbasfn,datafn,   &
           nx1,ny1,nz1,x1,y1,z1,soiltyp,stypfrct,xctr1,yctr1,  &
           nxbgn,nxend,nybgn,nyend,nzbgn,nzend, &
           tem1,tem2,tem11,tem21,intrp_opt,savespace,  &
           isx,jsy,iux,jvy,wgtsx,wgtsy,wgtux,wgtvy,temx1yz)

!
!-----------------------------------------------------------------------
!  Variable Declarations.
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of grid points in the soil
  INTEGER :: nstyps                    ! Number of soil type

  INTEGER :: hinfmt            ! The format of the history data dump
  CHARACTER(LEN=*) :: grdbasfn ! Name of the grid/base state array file
  CHARACTER(LEN=*) :: datafn   ! Name of the other time dependent data file

  REAL :: x(nx),y(ny),z(nz)

  INTEGER :: nx1,ny1,nz1
  REAL :: x1(nx1),y1(ny1),z1(nz1)

  REAL    :: stypfrct(nx,ny,nstyps)
  INTEGER :: soiltyp(nx,ny,nstyps)

  INTEGER :: nxbgn,nxend,nybgn,nyend,nzbgn,nzend

  REAL :: tem1(nx,ny,nz) ! Temporary array
  INTEGER(2) :: tem2(nx,ny,nz) ! Temporary array

  REAL :: tem11(nx1,ny1,nz1) ! Temporary array
  INTEGER(2) :: tem21(nx1,ny1,nz1) ! Temporary array

  integer :: savespace,INTRP_OPT

  integer :: isx(nx1),jsy(ny1),iux(nx1),jvy(ny1)
  real :: wgtsx(nx1),wgtsy(ny1),wgtux(nx1),wgtvy(ny1)
  real :: temx1yz(nx1,ny)

  INTEGER :: grdread,iread
  SAVE grdread

  INTEGER :: istat
  INTEGER :: ireturn           ! Return status indicator
  INTEGER :: grdbas            ! Wether this is a grid/base state
                               ! array dump
  INTEGER :: i,j,k
  LOGICAL :: fexist
  INTEGER :: is

  REAL :: xctr1,yctr1
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'indtflg.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'mp.inc'            ! mpi parameters.
  INCLUDE 'exbc.inc'
  INCLUDE 'phycst.inc'

  DATA grdread /0/
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ireturn = 0
  hinfmt = 3
!
!-----------------------------------------------------------------------
!
!  Open and read grid and base state data file depending on the
!  values of parameters grdin and basin, which are read in from the
!  time dependent data set. If grdin or basin is zero, the grid and
!  base state arrays have to be read in from a separate file.
!
!-----------------------------------------------------------------------
!
  IF( grdread == 0 ) THEN

!   print*,'grdread inside if block=', grdread

    grdbas = 1

    INQUIRE(FILE=grdbasfn, EXIST = fexist )
    IF( fexist ) GO TO 200

    INQUIRE(FILE=grdbasfn//'.Z', EXIST = fexist )
    IF( fexist ) THEN
      CALL uncmprs( grdbasfn//'.Z' )
      GO TO 200
    END IF

    INQUIRE(FILE=grdbasfn//'.gz', EXIST = fexist )
    IF( fexist ) THEN
      CALL uncmprs( grdbasfn//'.gz' )
      GO TO 200
    END IF

    WRITE(6,'(/1x,a,/1x,a/)')                                           &
        'File '//grdbasfn//                             &
        ' or its compressed version not found.',                        &
        'Program stopped in DTAREAD.'
    CALL arpsstop('arpsstop called from dtareaddump during base state read',1)

    200 CONTINUE
!
!-----------------------------------------------------------------------
!  Read grid and base state fields.
!-----------------------------------------------------------------------
!
    CALL hdfreaddump(nx,ny,nz,nzsoil,x,y,z,nstyps,grdbas,trim(grdbasfn), &
         nx1,ny1,nz1,x1,y1,z1,soiltyp,stypfrct,xctr1,yctr1,  &
         nxbgn,nxend,nybgn,nyend,nzbgn,nzend, &
         tem1,tem2,tem11,tem21,intrp_opt,savespace,      &
         isx,jsy,iux,jvy,wgtsx,wgtsy,wgtux,wgtvy,temx1yz,&
         0, 1,1,1,1,1,1, qcexout,qrexout,qiexout,qsexout,qgexout,qhexout,ngbrz,rayklow)

    grdread = 1

  END IF
!
!
!-----------------------------------------------------------------------
!  Read time dependent data fields.
!-----------------------------------------------------------------------
!
!
  grdbas = 0

  INQUIRE(FILE=trim(datafn), EXIST = fexist )
  IF( fexist ) GO TO 100

  INQUIRE(FILE=trim(datafn)//'.Z', EXIST = fexist )
  IF( fexist ) THEN
    CALL uncmprs( trim(datafn)//'.Z' )
    GO TO 100
  END IF

  INQUIRE(FILE=trim(datafn)//'.gz', EXIST = fexist )
  IF( fexist ) THEN
    CALL uncmprs( trim(datafn)//'.gz' )
    GO TO 100
  END IF

  WRITE(6,'(/1x,a,/1x,a/)')                                        &
     'File '//trim(datafn)                               &
     //' or its compressed version not found.',                    &
     'Program stopped in DTAREADDUMP.'
  CALL arpsstop('arpsstop called from dtareaddump during base read-2',1)

  100 CONTINUE

  CALL hdfreaddump(nx,ny,nz,nzsoil,x,y,z,nstyps,grdbas, trim(datafn), &
       nx1,ny1,nz1,x1,y1,z1,soiltyp,stypfrct,xctr1,yctr1,    &
       nxbgn,nxend,nybgn,nyend,nzbgn,nzend, &
       tem1,tem2,tem11,tem21,intrp_opt,savespace,      &
       isx,jsy,iux,jvy,wgtsx,wgtsy,wgtux,wgtvy,temx1yz,&
         exbcdmp, 1,1,1,1,1,1, qcexout,qrexout,qiexout,qsexout,qgexout,qhexout,ngbrz,rayklow)

  RETURN
END SUBROUTINE dtareaddump

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE HDFREAD                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE hdfreaddump(nx,ny,nz,nzsoil,x,y,z,nstyps,grdbas,filename,   &
           nx1,ny1,nz1,x1,y1,z1,soiltyp,stypfrct,xctr1,yctr1,          &
           nxbgn,nxend,nybgn,nyend,nzbgn,nzend, &
           tem1,itmp,tem1o,itmpo,intrp_opt,savespace, &
           isx,jsy,iux,jvy,wgtsx,wgtsy,wgtux,wgtvy,temx1yz,&
           exbcdmpopt, ubcdmp,vbcdmp,wbcdmp,ptbcdmp,prbcdmp,           &
           qvbcdmp,qcbcdmp,qrbcdmp,qibcdmp,qsbcdmp,qgbcdmp,qhbcdmp,ngbrz,rayklowest)


!-----------------------------------------------------------------------
!  PURPOSE:
!  Read in history data in the NCSA HDF4 format.
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  2000/04/15
!
!  MODIFICATION HISTORY:
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nzsoil   Number of grid points in the soil
!
!    grdbas   Data read flag.
!             =1, only grid and base state arrays will be read
!             =0, all arrays will be read based on data
!                 parameter setting.
!    filename  Character variable nhming the input HDF file

!-----------------------------------------------------------------------
!  Variable Declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: nx,ny,nz
  INTEGER :: nzsoil

  INTEGER :: grdbas
  CHARACTER (LEN=*) :: filename
  CHARACTER (LEN=256) :: filename_out
  INTEGER :: lfilename_out

  INTEGER :: itema, lenstr
  INTEGER :: iyr, imon, idy, ihr, imin, isec
  CHARACTER (LEN=256) :: exbcfn
  CHARACTER (LEN=256) :: temchar

  CHARACTER (LEN=15) :: ctime

  INTEGER :: exbcdmpopt

  INTEGER :: ubcdmp,vbcdmp,wbcdmp,ptbcdmp,prbcdmp,                        &
             qvbcdmp,qcbcdmp,qrbcdmp,qibcdmp,qsbcdmp,qgbcdmp,qhbcdmp
  INTEGER :: ngbrz, rayklowest

  INTEGER :: exbccompr

  REAL :: x     (nx)           ! x coord.
  REAL :: y     (ny)           ! y coord.
  REAL :: z     (nz)           ! z coord.
!  REAL :: zpsoil(nx,ny,nzsoil) ! physical x coord. for soil (m)
  REAL, ALLOCATABLE :: zpsoil(:,:,:) ! physical x coord. for soil (m)
  INTEGER :: nstyps

  INTEGER :: nx1,ny1,nz1
  REAL :: x1(nx1),y1(ny1),z1(nz1)
  REAL :: x1_out(nx1),y1_out(ny1)
!  REAL :: zpsoil1(nx1,ny1,nzsoil)
  REAL, ALLOCATABLE :: zpsoil1(:,:,:) ! physical x coord. for soil (m)

  REAL    :: stypfrct(nx,ny,nstyps)
  INTEGER :: soiltyp(nx,ny,nstyps)

  INTEGER :: nxbgn,nxend,nybgn,nyend,nzbgn,nzend

  integer :: intrp_opt

  integer :: isx(nx1),jsy(ny1),iux(nx1),jvy(ny1)
  real :: wgtsx(nx1),wgtsy(ny1),wgtux(nx1),wgtvy(ny1)
  real :: temx1yz(nx1,ny)

  REAL :: tem1(nx,ny,nz) ! Temporary array

  INTEGER (KIND=selected_int_kind(4)) :: itmp(nx,ny,nz) ! Temporary array
  REAL, ALLOCATABLE :: hmax(:), hmin(:) ! Temporary array

  REAL, ALLOCATABLE :: hmaxsoil(:), hminsoil(:) ! Temporary array
  INTEGER (KIND=selected_int_kind(4)), ALLOCATABLE :: itmpsoil(:,:,:)
  INTEGER (KIND=selected_int_kind(4)), ALLOCATABLE :: itmpsoil4(:,:,:,:)

  REAL :: tem1o(nx1,ny1,nz1)              ! Temporary work array
  INTEGER(2) :: itmpo(nx1,ny1,nz1)        ! Temporary array
  REAL, allocatable :: hmaxo(:), hmino(:)        ! Temporary array

  INTEGER(2), allocatable :: itmpsoilo(:,:,:)  ! Temporary array
  INTEGER(2), allocatable :: itmpsoilo4(:,:,:,:)  ! Temporary array
  REAL, allocatable :: hmaxsoilo(:), hminsoilo(:)! Temporary array

  INTEGER :: ireturn
  REAL :: time, dx1, dy1, temvar
  REAL :: ptmin, ptmax

! surface and soil variables
! by mhu
!  INTEGER, allocatable :: soiltyp(:,:,:)
!  REAL, allocatable :: stypfrct(:,:,:)
  INTEGER, allocatable :: vegtyp(:,:)

  INTEGER, allocatable :: soiltyp1(:,:,:)
  REAL, allocatable :: stypfrct1(:,:,:)
  INTEGER, allocatable :: vegtyp1(:,:)

  REAL, allocatable :: rsoil(:,:,:,:)
  REAL, allocatable :: rsoil1(:,:,:,:)
  REAL, allocatable :: wetcsoil(:,:,:)
  REAL, allocatable :: wetcsoil1(:,:,:)

!-----------------------------------------------------------------------
!  Parameters describing routine that wrote the gridded data
!-----------------------------------------------------------------------
!
! 06/28/2002 Zuwen He
!
! fmtver??: to label each data a version.
! intver??: an integer to allow faster comparison than fmtver??,
!           which are strings.
!
! Verion 5.00: significant change in soil variables since version 4.10.
!
!-----------------------------------------------------------------------

  CHARACTER (LEN=40) :: fmtver410,fmtver500,fmtver530
  INTEGER  :: intver,intver410,intver500,intver530

  PARAMETER (fmtver410='004.10 HDF4 Coded Data',intver410=410)
  PARAMETER (fmtver500='005.00 HDF4 Coded Data',intver500=500)
  PARAMETER (fmtver530='005.30 HDF4 Coded Data',intver530=530)

  CHARACTER (LEN=40) :: fmtverin

  CHARACTER (LEN=10) :: tmunit

!-----------------------------------------------------------------------
!  Map
!-----------------------------------------------------------------------

  REAL :: alatpro(2)
  REAL :: ctrlat1,ctrlon1
  REAL :: ctrx, ctry, dxscl, dyscl
  REAL :: xctr1,yctr1
  REAL :: swx,swy

!-----------------------------------------------------------------------
!  Terrain
!-----------------------------------------------------------------------

  REAL, allocatable  :: hterain(:,:)
  CHARACTER(len=256) :: ternfn
  INTEGER :: nunit,lternfn

!-----------------------------------------------------------------------
!  Misc. local variables
!-----------------------------------------------------------------------

  INTEGER :: lchanl
  PARAMETER (lchanl=6)      ! Channel number for formatted printing.

  INTEGER :: i,j,k,is
  INTEGER :: nxin,nyin,nzin,nzsoilin

  INTEGER :: bgrdin,bbasin,bvarin,bicein,btrbin,btkein

  INTEGER :: istat, sd_id,sd_id1, sd_id2
  INTEGER :: nstyp1,nstypin
  CHARACTER (LEN=30) :: substring
  INTEGER :: savespace
  INTEGER :: varflg

  CHARACTER (LEN=80) :: runname_old

  INTEGER :: nq, nqscalarin(nscalar)
  INTEGER :: hdfcomprtmp, exbccomprtmp

!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------

  CHARACTER(LEN=4) :: upcase

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'        ! Grid parameters
  INCLUDE 'indtflg.inc'
  INCLUDE 'alloc.inc'       ! allocation parameters & declarations
  INCLUDE 'phycst.inc'
  INCLUDE 'mp.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF (exbcdmpopt /= 0) THEN
    IF (exbcdmpopt == 3) THEN
      exbccompr = hdfcompr
    ELSE IF (exbcdmpopt == 4) THEN
      exbccompr = 5
    ELSE
      exbccompr = 0
    END IF
  ENDIF

!  if( nxbgn>999.or.nxend>999.or.nybgn>999.or.nyend>999 &
!                            .or.nzbgn>999.or.nzend>999) then
!    write(substring,'(6(a1,i4.4))') &
!    '.',nxbgn,'-',nxend,'x',nybgn,'-',nyend,'x',nzbgn,'-',nzend
!  else
!    write(substring,'(6(a1,i3.3))') &
!    '.',nxbgn,'-',nxend,'x',nybgn,'-',nyend,'x',nzbgn,'-',nzend
!  endif
  if( nxbgn>999.or.nxend>999.or.nybgn>999.or.nyend>999 &
                            .or.nzbgn>999.or.nzend>999) then
    write(substring,'(6(a1,i4.4))') &
    '.',nxend,'x',nyend,'x',nzend
  else
    write(substring,'(6(a1,i3.3))') &
    '.',nxend,'x',nyend,'x',nzend
  endif

  WRITE(*,*) 'HDFREAD: Reading HDF file: ', trim(filename)

  ALLOCATE (hmax(nz),stat=istat)
  IF (istat /= 0) THEN
    WRITE (6,*) "HDFREAD: ERROR allocating hmax, returning"
    RETURN
  END IF
  ALLOCATE (hmin(nz),stat=istat)
  IF (istat /= 0) THEN
    WRITE (6,*) "HDFREAD: ERROR allocating hmin, returning"
    RETURN
  END IF

!-----------------------------------------------------------------------
! Open file for reading
!-----------------------------------------------------------------------

  CALL hdfopen(filename,1,sd_id)
  IF (sd_id < 0) THEN
    WRITE (6,*) "HDFREAD: ERROR opening ",                              &
                 trim(filename)," for reading."
    GO TO 110
  END IF

!-----------------------------------------------------------------------
!  Create file for writing
!-----------------------------------------------------------------------


  lfilename_out = len_trim( filename )
  DO i=lfilename_out,1,-1
    if( filename(i:i) .eq. '/') exit
  enddo
  filename_out = trim(dirname)//filename(i+1:lfilename_out)//trim(substring)
  lfilename_out = len_trim( filename_out )

  CALL fnversn( filename_out, lfilename_out )

  if( hdmpfmt /= 0 ) then

  CALL hdfopen(filename_out(1:lfilename_out),2,sd_id1)
  IF (sd_id1 < 0) THEN
    WRITE (6,*) "HDFDUMP: ERROR creating HDF4 file: ",filename_out(1:lfilename_out)
    WRITE (6,*) "Program stopped."
    STOP
  END IF

  endif

  CALL hdfrdc(sd_id,40,"fmtver",fmtverin,istat)

  WRITE(6,'(/1x,a,a/)')                        &
      'Incoming data format, fmtverin=',fmtverin

  CALL hdfrdc(sd_id,80,"runname",runname_old,istat)
  runname=runname_old
  CALL hdfrdi(sd_id,"nocmnt",nocmnt,istat)
  IF( nocmnt > 0 ) THEN
    CALL hdfrdc(sd_id,80*nocmnt,"cmnt",cmnt,istat)
  END IF

  WRITE(6,'(//''  THE NAME OF THIS RUN IS:  '',A//)') trim(runname)

  WRITE (6,*) "Comments:"
  IF( nocmnt > 0 ) THEN
    DO i=1,nocmnt
      WRITE(6,'(1x,a)') cmnt(i)
    END DO
  END IF


  if( hdmpfmt /= 0 ) then

  CALL hdfwrtc(sd_id1, 40, 'fmtver', fmtverin, istat)

  IF (fmtverin == fmtver410) THEN
    intver=intver410
  ELSE IF (fmtverin == fmtver500) THEN
    intver=intver500
  ELSE IF (fmtverin == fmtver530) THEN
    intver=intver530
  ELSE
    intver=intver500
!    IF (myproc == 0) WRITE(6,'(/1x,a,a,a/)')                        &
!        'Incoming data format, fmtverin=',fmtverin,                 &
!        ', not found. The Job stopped.'
!    CALL arpsstop('arpstop called from HDFREAD. ',1)
  END IF

  IF (myproc == 0) WRITE(6,'(/1x,a,a/)')                        &
      'Incoming data format, fmtverin=',fmtverin

  CALL hdfwrtc(sd_id1, 80, 'runname', runname, istat)
  CALL hdfwrti(sd_id1, 'nocmnt', nocmnt, istat)
  IF( nocmnt > 0 ) THEN
    CALL hdfwrtc(sd_id1, 80*nocmnt, 'cmnt', cmnt, istat)
  END IF

  endif

  WRITE (6,*) " "

  CALL hdfrdc(sd_id,10,"tmunit",tmunit,istat)
  CALL hdfrdr(sd_id,"time",time,istat)


  if( hdmpfmt /= 0 ) then
    WRITE(6,'(1x,a,f13.3,a,a,a/)')                                          &
    'Writing HDF4 data at time=', curtim,' into file ',filename_out(1:lfilename_out)
  else
    WRITE(6,'(1x,a)') 'hdmpfmt =0, no history dump will be created.'
  endif

  if( hdmpfmt /= 0 ) then
  CALL hdfwrtc(sd_id1, 7, 'tmunit', 'seconds', istat)
  CALL hdfwrtr(sd_id1, 'time', time, istat)
  endif

  dx1 = x1(2)-x1(1)
  dy1 = y1(2)-y1(1)

   DO i=1,nx1
        x1_out(i) = x1(i)-x1(2)
   END DO

   DO j=1,ny1
        y1_out(j) = y1(j)-y1(2)
   END DO

!-----------------------------------------------------------------------
!  Get dimensions of data in binary file and check against
!  the dimensions passed to HDFREAD
!-----------------------------------------------------------------------

  CALL hdfrdi(sd_id,"nx",nxin,istat)
  CALL hdfrdi(sd_id,"ny",nyin,istat)
  CALL hdfrdi(sd_id,"nz",nzin,istat)

  IF ( nxin /= nx .OR. nyin /= ny .OR. nzin /= nz ) THEN
    WRITE(6,'(1x,a)') ' Dimensions in HDFREAD inconsistent with data.'
    WRITE(6,'(1x,a,3I15)') ' Read were: ', nxin, nyin, nzin
    WRITE(6,'(1x,a,3I15)') ' Expected:  ', nx, ny, nz
    WRITE(6,'(1x,a)') ' Program aborted in HDFREAD.'
    CALL arpsstop('arpsstop called from HDFREAD due to nxin...',1)
  END IF

!-----------------------------------------------------------------------
!  Read in x,y and z at grid cell centers (scalar points).
!-----------------------------------------------------------------------

  IF( grdin == 1 .OR. grdbas == 1 ) THEN
    CALL hdfrd1d(sd_id,"x",nx,x,istat)
    IF (istat /= 0) GO TO 110
    CALL hdfrd1d(sd_id,"y",ny,y,istat)
    IF (istat /= 0) GO TO 110
    CALL hdfrd1d(sd_id,"z",nz,z,istat)
    IF (istat /= 0) GO TO 110
!    CALL hdfrd3d(sd_id,"zpsoil",nx,ny,nzsoil,zpsoil,istat, &
!                 itmpsoil,hmaxsoil,hminsoil)
!    IF (istat /= 0) GO TO 110
  END IF  ! grdin


  if( hdmpfmt /= 0 ) then
  CALL hdfwrti(sd_id1, 'nx', nx1, istat)
  CALL hdfwrti(sd_id1, 'ny', ny1, istat)
  CALL hdfwrti(sd_id1, 'nz', nz1, istat)
  endif

  CALL hdfrdi(sd_id,"nzsoil",nzsoilin,istat)

  if( hdmpfmt /= 0 ) &
  CALL hdfwrti(sd_id1, 'nzsoil', nzsoilin, istat)

  IF (nzsoilin /= nzsoil) THEN

    WRITE(6,'(1x,a)') &
              ' Dimensions in HDFREADDUMP inconsistent with data.'
    WRITE(6,'(1x,a,I15)') ' Read were: ', nzsoilin
    WRITE(6,'(1x,a,I15)') ' Expected:  ', nzsoil
    WRITE(6,'(1x,a)') ' Program aborted in HDFREADDUMP.'

    CALL arpsstop('arpsstop called from HDFREADDUMP due to nzsoilin...',1)
  END IF

  IF (hdfcompr > 3 .or. exbccompr > 3 ) THEN
    ALLOCATE (hmaxo(nz1),stat=istat)
    IF (istat /= 0) THEN
      WRITE (6,*) "HDFDUMP: ERROR allocating hmaxo, returning"
      RETURN
    END IF
    ALLOCATE (hmino(nz1),stat=istat)
    IF (istat /= 0) THEN
      WRITE (6,*) "HDFDUMP: ERROR allocating hmino, returning"
      RETURN
    END IF

  END IF

!-----------------------------------------------------------------------
!
!  Read in flags for different data groups
!
!-----------------------------------------------------------------------
  IF ( grdbas == 1 ) THEN   ! Read grid and base state arrays

    WRITE(lchanl,'(1x,a,f8.1,a,f8.3,a/)')                               &
         'To read grid and base state data at time ', time,             &
         ' secs = ',(time/60.),' mins.'

    CALL hdfrdi(sd_id,"grdflg",bgrdin,istat)
    CALL hdfrdi(sd_id,"basflg",bbasin,istat)
    CALL hdfrdi(sd_id,"varflg",bvarin,istat)
    CALL hdfrdi(sd_id,"mstflg",mstin,istat)
    CALL hdfrdi(sd_id,"iceflg",bicein,istat)
    CALL hdfrdi(sd_id,"trbflg",btrbin,istat)
    CALL hdfrdi(sd_id,"landflg",landin,istat)
    CALL hdfrdi(sd_id,"totflg",totin,istat)
    CALL hdfrdi(sd_id,"tkeflg",btkein,istat)

    if( hdmpfmt /= 0 ) then

! mhu    landout = 0 ! writing code not implemented in this version

    CALL hdfwrti(sd_id1, 'grdflg', 1, istat)
    CALL hdfwrti(sd_id1, 'basflg', 1, istat)
    CALL hdfwrti(sd_id1, 'varflg', 0, istat)
    CALL hdfwrti(sd_id1, 'mstflg', 1, istat)
    CALL hdfwrti(sd_id1, 'iceflg', 0, istat)
    CALL hdfwrti(sd_id1, 'trbflg', 0, istat)
    CALL hdfwrti(sd_id1, 'sfcflg', 0, istat)
    CALL hdfwrti(sd_id1, 'rainflg', 0, istat)
    CALL hdfwrti(sd_id1, 'landflg', landin*landout, istat)
    CALL hdfwrti(sd_id1, 'totflg', 1, istat)
    CALL hdfwrti(sd_id1, 'tkeflg', 0, istat)
    endif

  ELSE ! Normal data reading

    WRITE(lchanl,'(1x,a,f8.1,a,f8.3,a/)')'To read data for time:',      &
         time,' secs = ',(time/60.),' mins.'

    CALL hdfrdi(sd_id,"grdflg",grdin,istat)
    CALL hdfrdi(sd_id,"basflg",basin,istat)
    CALL hdfrdi(sd_id,"varflg",varin,istat)
    CALL hdfrdi(sd_id,"mstflg",mstin,istat)
    CALL hdfrdi(sd_id,"iceflg",icein,istat)
    CALL hdfrdi(sd_id,"trbflg",trbin,istat)
    CALL hdfrdi(sd_id,"sfcflg",sfcin,istat)
    CALL hdfrdi(sd_id,"rainflg",rainin,istat)
    CALL hdfrdi(sd_id,"landflg",landin,istat)
    CALL hdfrdi(sd_id,"totflg",totin,istat)
    CALL hdfrdi(sd_id,"tkeflg",tkein,istat)

    IF (intver >= intver530) THEN       ! New version

      CALL hdfrdi(sd_id,"nscalar",nscalarin,istat)
      DO nq = 1,nscalar
        CALL hdfrdi(sd_id,'P_'//upcase(qnames(nq)),nqscalarin(nq),istat)
      END DO

    ELSE                        ! maybe old version
      nscalarin = 0
      nqscalarin(:) = 0

      IF (mstin == 1) THEN

        IF (P_QC > 0) THEN
           nqscalarin(P_QC) = 1
           nscalarin = nscalarin + 1
        END IF
        IF (P_QR > 0) THEN
           nqscalarin(P_QR) = 2
           nscalarin = nscalarin + 1
        END IF

        IF (icein == 1) THEN
          IF (P_QI > 0) THEN
             nqscalarin(P_QI) = 3
             nscalarin = nscalarin + 1
          END IF
          IF (P_QS > 0) THEN
             nqscalarin(P_QS) = 4
             nscalarin = nscalarin + 1
          END IF
          IF (P_QH > 0) THEN
             nqscalarin(P_QH) = 5
             nscalarin = nscalarin + 1
          END IF
        END IF

      END IF

    END IF                   ! new version or old version


    if( hdmpfmt /= 0 ) then

!    landout=0   ! no land data for normal dump file

    CALL hdfwrti(sd_id1, 'grdflg', grdin*grdout, istat)
    CALL hdfwrti(sd_id1, 'basflg', basin*basout, istat)

    if( varin == 0 ) then
      varflg = 0
    elseif( varin == 1 .and. varout == 1 ) then
      varflg = 1
    elseif( varin == 1 .and. (varout == 2 .or. varout == 3) ) then
      varflg = varout
    elseif( varin == 2 .and. varout == 2 ) then
      varflg = 2
    elseif( varin == 3 .and. varout == 3 ) then
      varflg = 3
    else
      varflg = 0
    endif
    CALL hdfwrti(sd_id1, 'varflg', varflg , istat)

    CALL hdfwrti(sd_id1, 'mstflg', mstin*mstout, istat)
    CALL hdfwrti(sd_id1, 'iceflg', icein*iceout, istat)
    CALL hdfwrti(sd_id1, 'trbflg', trbin*trbout, istat)
    CALL hdfwrti(sd_id1, 'sfcflg', sfcin*sfcout, istat)
    CALL hdfwrti(sd_id1, 'rainflg',rainin*rainout, istat)

    CALL hdfwrti(sd_id1, 'landflg',landin*landout, istat)
    CALL hdfwrti(sd_id1, 'totflg', 1, istat)
    CALL hdfwrti(sd_id1, 'tkeflg', tkein*tkeout, istat)

    CALL hdfwrti(sd_id1, 'nscalar', nscalar, istat)
    CALL hdfwrti(sd_id1, 'P_QC',    P_QC, istat)
    CALL hdfwrti(sd_id1, 'P_QR',    P_QR, istat)
    CALL hdfwrti(sd_id1, 'P_QI',    P_QI, istat)
    CALL hdfwrti(sd_id1, 'P_QS',    P_QS, istat)
    CALL hdfwrti(sd_id1, 'P_QG',    P_QG, istat)
    CALL hdfwrti(sd_id1, 'P_QH',    P_QH, istat)

    CALL hdfwrti(sd_id1, 'P_NC',    P_NC, istat)
    CALL hdfwrti(sd_id1, 'P_NR',    P_NR, istat)
    CALL hdfwrti(sd_id1, 'P_NI',    P_NI, istat)
    CALL hdfwrti(sd_id1, 'P_NS',    P_NS, istat)
    CALL hdfwrti(sd_id1, 'P_NG',    P_NG, istat)
    CALL hdfwrti(sd_id1, 'P_NH',    P_NH, istat)

    CALL hdfwrti(sd_id1, 'P_ZR',    P_ZR, istat)
    CALL hdfwrti(sd_id1, 'P_ZI',    P_ZI, istat)
    CALL hdfwrti(sd_id1, 'P_ZS',    P_ZS, istat)
    CALL hdfwrti(sd_id1, 'P_ZG',    P_ZG, istat)
    CALL hdfwrti(sd_id1, 'P_ZH',    P_ZH, istat)

    endif

  END IF

  CALL hdfrdi(sd_id,"nstyp",nstyp1,istat)

  IF ( nstyp1 < 1 ) THEN
    nstyp1 = 1
  END IF

  if( hdmpfmt /= 0 ) &
  CALL hdfwrti(sd_id1, 'nstyp', nstyp1, istat)

  IF (nstyp1 > nstyp) THEN
    WRITE (6,*) "HDFREAD: WARNING, nstyp in file (",nstyp1,  &
       ") greater than that specified in input file (",nstyp,  &
       "), using only",nstyp
    nstypin = nstyp
  ELSE
    nstypin = nstyp1
  ENDIF

  CALL hdfrdi(sd_id,"prcflg",prcin,istat)
  CALL hdfrdi(sd_id,"radflg",radin,istat)
  CALL hdfrdi(sd_id,"flxflg",flxin,istat)
  CALL hdfrdi(sd_id,"snowflg",snowin,istat)

  CALL hdfrdi(sd_id,"month",month,istat)
  CALL hdfrdi(sd_id,"day",day,istat)
  CALL hdfrdi(sd_id,"year",year,istat)
  CALL hdfrdi(sd_id,"hour",hour,istat)
  CALL hdfrdi(sd_id,"minute",minute,istat)
  CALL hdfrdi(sd_id,"second",second,istat)

  CALL hdfrdr(sd_id,"umove",umove,istat)
  CALL hdfrdr(sd_id,"vmove",vmove,istat)
  CALL hdfrdr(sd_id,"xgrdorg",xgrdorg,istat)
  CALL hdfrdr(sd_id,"ygrdorg",ygrdorg,istat)

  CALL hdfrdi(sd_id,"mapproj",mapproj,istat)
  CALL hdfrdr(sd_id,"trulat1",trulat1,istat)
  CALL hdfrdr(sd_id,"trulat2",trulat2,istat)
  CALL hdfrdr(sd_id,"trulon",trulon,istat)
  CALL hdfrdr(sd_id,"sclfct",sclfct,istat)
  CALL hdfrdr(sd_id,"tstop",tstop,istat)
  CALL hdfrdr(sd_id,"thisdmp",thisdmp,istat)
  CALL hdfrdr(sd_id,"latitud",latitud,istat)
  CALL hdfrdr(sd_id,"ctrlat",ctrlat,istat)
  CALL hdfrdr(sd_id,"ctrlon",ctrlon,istat)

    alatpro(1)=trulat1
    alatpro(2)=trulat2
    CALL setmapr( mapproj,sclfct,alatpro,trulon )
    CALL lltoxy( 1,1, ctrlat,ctrlon, ctrx, ctry )
    dxscl=x(2)-x(1)
    dyscl=y(2)-y(1)
    swx = ctrx - (FLOAT(nx-3)/2.)*dxscl
    swy = ctry - (FLOAT(ny-3)/2.)*dyscl
    CALL setorig( 1, swx, swy)
    CALL xytoll(1,1,xctr1,yctr1, ctrlat1,ctrlon1)
    write(*,*) 'The new domain has a center point at:'
    write(*,*) 'lat=',ctrlat1,'lon=',ctrlon1

  CALL hdfrdr(sd_id,"ntcloud", ntcloud,istat)
  CALL hdfrdr(sd_id,"n0rain", n0rain,istat)
  CALL hdfrdr(sd_id,"n0snow", n0snow,istat)
  CALL hdfrdr(sd_id,"n0grpl", n0grpl,istat)
  CALL hdfrdr(sd_id,"n0hail", n0hail,istat)
  CALL hdfrdr(sd_id,"rhoice",rhoice,istat)
  CALL hdfrdr(sd_id,"rhosnow",rhosnow,istat)
  CALL hdfrdr(sd_id,"rhogrpl",rhogrpl,istat)
  CALL hdfrdr(sd_id,"rhohail",rhohail,istat)
  CALL hdfrdr(sd_id,"alpharain",alpharain,istat)
  CALL hdfrdr(sd_id,"alphaice",alphaice,istat)
  CALL hdfrdr(sd_id,"alphasnow",alphasnow,istat)
  CALL hdfrdr(sd_id,"alphagrpl",alphagrpl,istat)
  CALL hdfrdr(sd_id,"alphahail",alphahail,istat)

  if( hdmpfmt /= 0 ) then
  CALL hdfwrti(sd_id1, 'prcflg', prcin*prcout, istat)
  CALL hdfwrti(sd_id1, 'radflg', radin*radout, istat)
  CALL hdfwrti(sd_id1, 'flxflg', flxout*flxout, istat)
  CALL hdfwrti(sd_id1, 'snowflg', snowin*snowout, istat)

  CALL hdfwrti(sd_id1, 'day', day, istat)
  CALL hdfwrti(sd_id1, 'year', year, istat)
  CALL hdfwrti(sd_id1, 'month', month, istat)
  CALL hdfwrti(sd_id1, 'hour', hour, istat)
  CALL hdfwrti(sd_id1, 'minute', minute, istat)
  CALL hdfwrti(sd_id1, 'second', second, istat)

  CALL hdfwrtr(sd_id1, 'umove', umove, istat)
  CALL hdfwrtr(sd_id1, 'vmove', vmove, istat)
  CALL hdfwrtr(sd_id1, 'xgrdorg', xgrdorg, istat)
  CALL hdfwrtr(sd_id1, 'ygrdorg', ygrdorg, istat)

  CALL hdfwrti(sd_id1, 'mapproj', mapproj, istat)
  CALL hdfwrtr(sd_id1, 'trulat1', trulat1, istat)
  CALL hdfwrtr(sd_id1, 'trulat2', trulat2, istat)
  CALL hdfwrtr(sd_id1, 'trulon', trulon, istat)
  CALL hdfwrtr(sd_id1, 'sclfct', sclfct, istat)
  CALL hdfwrtr(sd_id1, 'tstop', tstop, istat)
  CALL hdfwrtr(sd_id1, 'thisdmp', thisdmp, istat)
  CALL hdfwrtr(sd_id1, 'latitud', latitud, istat)
  CALL hdfwrtr(sd_id1, 'ctrlat', ctrlat1, istat)
  CALL hdfwrtr(sd_id1, 'ctrlon', ctrlon1, istat)

  CALL hdfwrtr(sd_id1, 'ntcloud',  ntcloud,  istat)
  CALL hdfwrtr(sd_id1, 'n0rain',  n0rain,  istat)
  CALL hdfwrtr(sd_id1, 'n0snow',  n0snow,  istat)
  CALL hdfwrtr(sd_id1, 'n0grpl',  n0grpl,  istat)
  CALL hdfwrtr(sd_id1, 'n0hail',  n0hail,  istat)
  CALL hdfwrtr(sd_id1, 'rhoice', rhoice, istat)
  CALL hdfwrtr(sd_id1, 'rhosnow', rhosnow, istat)
  CALL hdfwrtr(sd_id1, 'rhogrpl', rhogrpl, istat)
  CALL hdfwrtr(sd_id1, 'rhohail', rhohail, istat)
  CALL hdfwrtr(sd_id1, 'alpharain', alpharain, istat)
  CALL hdfwrtr(sd_id1, 'alphaice', alphaice, istat)
  CALL hdfwrtr(sd_id1, 'alphasnow', alphasnow, istat)
  CALL hdfwrtr(sd_id1, 'alphagrpl', alphagrpl, istat)
  CALL hdfwrtr(sd_id1, 'alphahail', alphahail, istat)

  CALL hdfwrtr(sd_id1, 'dx', dx1, istat)
  CALL hdfwrtr(sd_id1, 'dy', dy1, istat)
  endif

!-----------------------------------------------------------------------
!  Write exbc file
!-----------------------------------------------------------------------

  IF (exbcdmpopt /= 0) THEN

    CALL ctim2abss( year,month,day,hour,minute,second, itema )
    itema = itema + INT(time)
    CALL abss2ctim( itema,iyr,imon,idy,ihr,imin,isec )
    WRITE (ctime,'(i4.4,2i2.2,a,3i2.2)')iyr,imon,idy,'.',ihr,imin,isec

    lfnkey = len_trim(runname)
    CALL gtlfnkey( runname, lfnkey )
    exbcfn = runname(1:lfnkey)//'.'//ctime
    lenstr = lfnkey + 16

    IF( dirname /= ' ' ) THEN
      temchar = exbcfn
      exbcfn  = trim(dirname)//'/'//temchar
      lenstr  = len_trim(exbcfn)
    END IF

    !CALL fnversn(exbcfn,lenstr)

    WRITE(6,'(1x,a,a)')                                                 &
         'Dumping to the external boundary format file: ',exbcfn(1:lenstr)

    CALL hdfopen(exbcfn(1:lenstr), 2, sd_id2)
    IF (sd_id2 < 0) THEN
      WRITE (6,*) "HDFREADDUMP: ERROR creating HDF4 file: ",exbcfn(1:lenstr)
      STOP
    END IF

    CALL hdfwrtc(sd_id2, 40, 'fmtver', fmtver500, istat)
    CALL hdfwrtc(sd_id2, 15, 'ctime', ctime, istat)

    CALL hdfwrti(sd_id2, 'nx', nx1, istat)
    CALL hdfwrti(sd_id2, 'ny', ny1, istat)
    CALL hdfwrti(sd_id2, 'nz', nz1, istat)
    CALL hdfwrtr(sd_id2, 'dx', dx1, istat)
    CALL hdfwrtr(sd_id2, 'dy', dy1, istat)
    CALL hdfwrtr(sd_id2, 'dz', dz,  istat)

    CALL hdfwrtr(sd_id2, 'dzmin', dzmin, istat)            ! no known yet
    CALL hdfwrti(sd_id2, 'strhopt', strhopt, istat)
    CALL hdfwrtr(sd_id2, 'zrefsfc', zrefsfc, istat)
    CALL hdfwrtr(sd_id2, 'dlayer1', dlayer1, istat)
    CALL hdfwrtr(sd_id2, 'dlayer2', dlayer2, istat)
    CALL hdfwrtr(sd_id2, 'zflat', zflat, istat)
    CALL hdfwrtr(sd_id2, 'strhtune', strhtune, istat)

    CALL hdfwrti(sd_id2, 'mapproj', mapproj, istat)
    CALL hdfwrtr(sd_id2, 'trulat1', trulat1, istat)
    CALL hdfwrtr(sd_id2, 'trulat2', trulat2, istat)
    CALL hdfwrtr(sd_id2, 'trulon', trulon, istat)
    CALL hdfwrtr(sd_id2, 'sclfct', sclfct, istat)
    CALL hdfwrtr(sd_id2, 'ctrlat', ctrlat1, istat)
    CALL hdfwrtr(sd_id2, 'ctrlon', ctrlon1, istat)

    CALL hdfwrti(sd_id2,"ubcflg",ubcdmp,istat)
    CALL hdfwrti(sd_id2,"vbcflg",vbcdmp,istat)
    CALL hdfwrti(sd_id2,"wbcflg",wbcdmp,istat)
    CALL hdfwrti(sd_id2,"ptbcflg",ptbcdmp,istat)
    CALL hdfwrti(sd_id2,"prbcflg",prbcdmp,istat)
    CALL hdfwrti(sd_id2,"qvbcflg",qvbcdmp,istat)
    CALL hdfwrti(sd_id2,"qcbcflg",qcbcdmp,istat)
    CALL hdfwrti(sd_id2,"qrbcflg",qrbcdmp,istat)
    CALL hdfwrti(sd_id2,"qibcflg",qibcdmp,istat)
    CALL hdfwrti(sd_id2,"qsbcflg",qsbcdmp,istat)
    CALL hdfwrti(sd_id2,"qgbcflg",qgbcdmp,istat)
    CALL hdfwrti(sd_id2,"qhbcflg",qhbcdmp,istat)

    IF (exbcdmpopt == 4) THEN
      CALL hdfwrti(sd_id2, 'clipxy', ngbrz, istat)
      CALL hdfwrti(sd_id2, 'clipz', rayklowest, istat)
    END IF

  ENDIF

!-----------------------------------------------------------------------
!  If grdout=1 or grdbas=1, write out grid variables
!-----------------------------------------------------------------------

  IF((grdout == 1 .OR. grdbas == 1).and.(hdmpfmt /= 0) ) THEN

    CALL hdfwrt1d(x1_out,nx1,sd_id1,'x','x coordinate','m')
    CALL hdfwrt1d(y1_out,ny1,sd_id1,'y','y coordinate','m')
    CALL hdfwrt1d(z1,nz1,sd_id1,'z','z coordinate','m')

    IF( intver >= intver500 .and. (grdin*grdout == 1 .or. grdbas == 1) ) then
      ALLOCATE (zpsoil(nx,ny,nzsoil),stat=istat)
      IF (istat /= 0) THEN
          WRITE (6,*) "HDFDUMP: ERROR allocating zpsoil, returning"
          RETURN
      END IF
      ALLOCATE (zpsoil1(nx1,ny1,nzsoil),stat=istat)
      IF (istat /= 0) THEN
          WRITE (6,*) "HDFDUMP: ERROR allocating zpsoil1, returning"
          RETURN
      END IF
      ALLOCATE (itmpsoil(nx,ny,nzsoil),stat=istat)
      IF (istat /= 0) THEN
        WRITE (6,*) "HDFDUMP: ERROR allocating itmpsoil, returning"
        RETURN
      END IF
      ALLOCATE (hmaxsoil(nzsoil),stat=istat)
      IF (istat /= 0) THEN
        WRITE (6,*) "HDFDUMP: ERROR allocating hmaxsoil, returning"
        RETURN
      END IF
      ALLOCATE (hminsoil(nzsoil),stat=istat)
      IF (istat /= 0) THEN
        WRITE (6,*) "HDFDUMP: ERROR allocating hminisoil, returning"
        RETURN
      END IF
      CALL hdfrd3d(sd_id,"zpsoil",nx,ny,nzsoil,zpsoil,istat, &
                 itmpsoil,hmaxsoil,hminsoil)
      DEALLOCATE(itmpsoil)
      DEALLOCATE(hmaxsoil)
      DEALLOCATE(hminsoil)
      IF (istat /= 0) GO TO 110

      if( intrp_opt == 0 ) then
        call copyarray(zpsoil,nx,ny,nzsoil, zpsoil1, nx1,ny1,nzsoil,    &
            nxbgn,nxend,nybgn,nyend,1,nzsoil)
      else
        call intrpxy3d(zpsoil,nx,1,nx-1,ny,1,ny-1,nzsoil,1,nzsoil,      &
             wgtsx,isx,wgtsy,jsy, intrp_opt,                            &
             zpsoil1,nx1,1,nx1-1,ny1,1,ny1-1, temx1yz)
      endif
      DEALLOCATE(zpsoil)

      ALLOCATE (itmpsoilo(nx1,ny1,nzsoil),stat=istat)
      IF (istat /= 0) THEN
        WRITE (6,*) "HDFDUMP: ERROR allocating itmpsoilo, returning"
        RETURN
      END IF
      ALLOCATE (hmaxsoilo(nzsoil),stat=istat)
      IF (istat /= 0) THEN
        WRITE (6,*) "HDFDUMP: ERROR allocating hmaxsoilo, returning"
        RETURN
      END IF
      ALLOCATE (hminsoilo(nzsoil),stat=istat)
      IF (istat /= 0) THEN
        WRITE (6,*) "HDFDUMP: ERROR allocating hminisoilo, returning"
        RETURN
      END IF
      CALL hdfwrt3d(zpsoil1,nx1,ny1,nzsoil,sd_id1,0,hdfcompr,             &
                'zpsoil','Physical height coordinate (soil)','m',   &
                 itmpsoilo,hmaxsoilo,hminsoilo)
      DEALLOCATE(zpsoil1)
      DEALLOCATE(itmpsoilo)
      DEALLOCATE(hmaxsoilo)
      DEALLOCATE(hminsoilo)
    ENDIF

    if(grdbas == 1 .or. grdin*grdout == 1) then
    CALL hdfrd3d(sd_id,"zp",nx,ny,nz,tem1,istat,itmp,hmax,hmin)
    IF (istat /= 0) GO TO 110

    if( intrp_opt == 0 ) then
      call copyarray(tem1,nx,ny,nz, tem1o, nx1,ny1,nz1, &
           nxbgn,nxend,nybgn,nyend,nzbgn,nzend)
    else
      call intrpxy3d(tem1,nx,1,nx-1,ny,1,ny-1,nz,nzbgn,nzend,           &
                wgtsx,isx,wgtsy,jsy,intrp_opt,                          &
                tem1o,nx1,1,nx1-1,ny1,1,ny1-1, temx1yz)
    endif
    endif

    if(grdbas == 1 .or. grdin*grdout == 1) THEN
    CALL hdfwrt3d(tem1o,nx1,ny1,nz1,sd_id1,0,hdfcompr,                 &
                 'zp','Physical height coordinate','m',              &
                  itmpo,hmaxo,hmino)
    endif

    if(grdbas == 1 .and. terndmp >= 1 ) then
!
!-----------------------------------------------------------------------
!
!  Write out terrain data
!
!-----------------------------------------------------------------------
!
        ALLOCATE (hterain(nx1,ny1),stat=istat)
        IF (istat /= 0) THEN
          WRITE (6,*) "HDFDUMP: ERROR allocating hterain, returning"
          RETURN
        END IF
        hterain(:,:)=tem1o(:,:,2)

        CALL getunit( nunit )

        lfilename_out = len_trim( runname )
        DO i=1,lfilename_out
           if( runname(i:i) .eq.',' .or. runname(i:i) .eq.' ' ) exit
        enddo
        filename_out = trim(dirname)//runname(1:i-1)

        ternfn = trim(filename_out)//".trndata"
        lternfn = len_trim( ternfn )
        CALL fnversn(ternfn, lternfn )

        PRINT *, 'Write terrain data to ',trim(ternfn)

        CALL writtrn(nx1,ny1,ternfn(1:lternfn), dx1,dy1,                &
                   mapproj,trulat1,trulat2,trulon,sclfct,               &
                   ctrlat1,ctrlon1,hterain)

    endif  ! terrain

  ENDIF

!-----------------------------------------------------------------------
!  Read in base state fields
!-----------------------------------------------------------------------

! Print*,'start doing 3d arrays'

  IF( (basout== 1 .OR. grdbas == 1)  .and. hdmpfmt /= 0) THEN

    CALL hdfrd3d(sd_id,"ubar",nx,ny,nz,tem1,istat,itmp,hmax,hmin)
    IF (istat /= 0) GO TO 110

    if( intrp_opt == 0 ) then
    call copyarray(tem1,nx,ny,nz, tem1o, nx1,ny1,nz1, &
         nxbgn,nxend,nybgn,nyend,nzbgn,nzend)
    else
      call intrpxy3d(tem1,nx,1,nx,ny,1,ny-1,nz,nzbgn,nzend,  &
           wgtux,iux,wgtsy,jsy,intrp_opt,tem1o,nx1,1,nx1,ny1,1,ny1-1, temx1yz)
    endif

    if( hdmpfmt /= 0 )  &
    CALL hdfwrt3d(tem1o,nx1,ny1,nz1,sd_id1,1,hdfcompr,            &
                  'ubar','Base state u-velocity','m/s',                 &
                   itmpo,hmaxo,hmino)

    CALL hdfrd3d(sd_id,"vbar",nx,ny,nz,tem1,istat,itmp,hmax,hmin)
    IF (istat /= 0) GO TO 110
    if( intrp_opt == 0 ) then
    call copyarray(tem1,nx,ny,nz, tem1o, nx1,ny1,nz1, &
         nxbgn,nxend,nybgn,nyend,nzbgn,nzend)
    else
      call intrpxy3d(tem1,nx,1,nx-1,ny,1,ny,nz,nzbgn,nzend,  &
           wgtsx,isx,wgtvy,jvy,intrp_opt,tem1o,nx1,1,nx1-1,ny1,1,ny1, temx1yz)
    endif

    if( hdmpfmt /= 0 )  &
    CALL hdfwrt3d(tem1o,nx1,ny1,nz1,sd_id1,1,hdfcompr,            &
                  'vbar','Base state v-velocity','m/s',                 &
                   itmpo,hmaxo,hmino)

    CALL hdfrd3d(sd_id,"wbar",nx,ny,nz,tem1,istat,itmp,hmax,hmin)
    IF (istat /= 0) GO TO 110
    tem1o = 0.0
    if( hdmpfmt /= 0 )  &
    CALL hdfwrt3d(tem1o,nx1,ny1,nz1,sd_id1,1,hdfcompr,            &
                  'wbar','Base state w-velocity','m/s',                 &
                   itmpo,hmaxo,hmino)

    CALL hdfrd3d(sd_id,"ptbar",nx,ny,nz,tem1,istat,itmp,hmax,hmin)
    IF (istat /= 0) GO TO 110
    if( intrp_opt == 0 ) then
      call copyarray(tem1,nx,ny,nz, tem1o, nx1,ny1,nz1, &
           nxbgn,nxend,nybgn,nyend,nzbgn,nzend)
    else
      call intrpxy3d(tem1,nx,1,nx-1,ny,1,ny-1,nz,nzbgn,nzend,  &
           wgtsx,isx,wgtsy,jsy,intrp_opt,tem1o,nx1,1,nx1-1,ny1,1,ny1-1, temx1yz)
    endif

    if( hdmpfmt /= 0 )  &
    CALL hdfwrt3d(tem1o,nx1,ny1,nz1,sd_id1,1,hdfcompr,            &
                  'ptbar','Base state potential temperature','K',       &
                   itmpo,hmaxo,hmino)

    print*,'ptbar(1,1,1)=', tem1o(1,1,1)

    CALL hdfrd3d(sd_id,"pbar",nx,ny,nz,tem1,istat,itmp,hmax,hmin)
    IF (istat /= 0) GO TO 110
    if( intrp_opt == 0 ) then
    call copyarray(tem1,nx,ny,nz, tem1o, nx1,ny1,nz1, &
         nxbgn,nxend,nybgn,nyend,nzbgn,nzend)
    else
      call intrpxy3d(tem1,nx,1,nx-1,ny,1,ny-1,nz,nzbgn,nzend,  &
           wgtsx,isx,wgtsy,jsy,intrp_opt,tem1o,nx1,1,nx1-1,ny1,1,ny1-1, temx1yz)
    endif
    if( hdmpfmt /= 0 )  &
    CALL hdfwrt3d(tem1o,nx1,ny1,nz1,sd_id1,1,hdfcompr,            &
                  'pbar','Base state pressure','Pascal',                &
                   itmpo,hmaxo,hmino)

    IF( mstin == 1 ) THEN
      CALL hdfrd3d(sd_id,"qvbar",nx,ny,nz,tem1,istat,itmp,hmax,hmin)
      IF (istat /= 0) GO TO 110
    ELSE
      tem1 = 0.0
    END IF
!   IF(mstout == 1 ) THEN
! always output qvbar

      if( intrp_opt == 0 ) then
        call copyarray(tem1,nx,ny,nz, tem1o, nx1,ny1,nz1, &
           nxbgn,nxend,nybgn,nyend,nzbgn,nzend)
      else
        call intrpxy3d(tem1,nx,1,nx-1,ny,1,ny-1,nz,nzbgn,nzend,  &
             wgtsx,isx,wgtsy,jsy,intrp_opt,tem1o,nx1,1,nx1-1,ny1,1,ny1-1, temx1yz)
      endif
      if( hdmpfmt /= 0 )  &
      CALL hdfwrt3d(tem1o,nx1,ny1,nz1,sd_id1,1,hdfcompr,            &
            'qvbar','Base state water vapor specific humidity','kg/kg',   &
                     itmpo,hmaxo,hmino)
!   END IF

    IF(landin == 1 .AND. landout == 1  .and. hdmpfmt /= 0 ) THEN

!       allocate(soiltyp(nx,ny,nstyp),stat=istat)
!       IF (istat /= 0) THEN
!          WRITE (6,*) "HDFDUMP: ERROR allocating soiltyp, returning"
!          RETURN
!       END IF
       allocate(soiltyp1(nx1,ny1,nstyp1),stat=istat)
       IF (istat /= 0) THEN
          WRITE (6,*) "HDFDUMP: ERROR allocating soiltyp1, returning"
          RETURN
       END IF
!!       allocate(stypfrct(nx,ny,nstyp),stat=istat)
!       IF (istat /= 0) THEN
!          WRITE (6,*) "HDFDUMP: ERROR allocating stypfrct, returning"
!          RETURN
!       END IF
       allocate(stypfrct1(nx1,ny1,nstyp1),stat=istat)
       IF (istat /= 0) THEN
          WRITE (6,*) "HDFDUMP: ERROR allocating stypfrct1, returning"
          RETURN
       END IF
       allocate(vegtyp(nx,ny),stat=istat)
       IF (istat /= 0) THEN
          WRITE (6,*) "HDFDUMP: ERROR allocating vegtyp, returning"
          RETURN
       END IF
       allocate(vegtyp1(nx1,ny1),stat=istat)
       IF (istat /= 0) THEN
          WRITE (6,*) "HDFDUMP: ERROR allocating vegtyp1, returning"
          RETURN
       END IF
       allocate(itmpsoil(nx,ny,nstyp),stat=istat)
       IF (istat /= 0) THEN
          WRITE (6,*) "HDFDUMP: ERROR allocating itmpsoil, returning"
          RETURN
       END IF
       allocate(hmaxsoil(nstyp),stat=istat)
       IF (istat /= 0) THEN
          WRITE (6,*) "HDFDUMP: ERROR allocating hmaxsoil, returning"
          RETURN
       END IF
       allocate(hminsoil(nstyp),stat=istat)
       IF (istat /= 0) THEN
          WRITE (6,*) "HDFDUMP: ERROR allocating hminsoil, returning"
          RETURN
       END IF

       CALL hdfrd3di(sd_id,"soiltyp",nx,ny,nstyp,soiltyp,istat)
       CALL hdfrd3d(sd_id,"stypfrct",nx,ny,nstyp,stypfrct,istat,     &
                 itmpsoil,hmaxsoil,hminsoil)
       CALL hdfrd2di(sd_id,"vegtyp",nx,ny,vegtyp,istat)

       deallocate(itmpsoil)
       deallocate(hmaxsoil)
       deallocate(hminsoil)

       if( intrp_opt == 0 ) then
         call copyarrayi(soiltyp,nx,ny,nstyp,soiltyp1,nx1,ny1,nstyp1,  &
                        nxbgn,nxend,nybgn,nyend,1,nstyp1)
         call copyarray(stypfrct,nx,ny,nstyp,stypfrct1,nx1,ny1,nstyp1,  &
                        nxbgn,nxend,nybgn,nyend,1,nstyp1)
         call copyarrayi(vegtyp,nx,ny,1,vegtyp1,nx1,ny1,1,  &
                        nxbgn,nxend,nybgn,nyend,1,1)
       else
         CALL intrp_soil_int(nx,ny,nx1,ny1,nstyp,nstyp1,wgtsx,wgtsy,isx,jsy, &
                        soiltyp,stypfrct,vegtyp,                    &
                        soiltyp1,stypfrct1,vegtyp1)
       endif

       allocate(itmpsoilo(nx1,ny1,nstyp1),stat=istat)
       IF (istat /= 0) THEN
          WRITE (6,*) "HDFDUMP: ERROR allocating hminsoilo, returning"
          RETURN
       END IF
       allocate(hmaxsoilo(nstyp1),stat=istat)
       IF (istat /= 0) THEN
          WRITE (6,*) "HDFDUMP: ERROR allocating hminsoilo, returning"
          RETURN
       END IF
       allocate(hminsoilo(nstyp1),stat=istat)
       IF (istat /= 0) THEN
          WRITE (6,*) "HDFDUMP: ERROR allocating hminsoilo, returning"
          RETURN
       END IF

       CALL hdfwrt3di(soiltyp1,nx1,ny1,nstyp1,sd_id1,0,0,             &
                  'soiltyp','Soil type','index')
       CALL hdfwrt3d(stypfrct1,nx1,ny1,nstyp1,sd_id1,0,hdfcompr,      &
                  'stypfrct','Soil type fractional coverage','fraction',  &
                  itmpsoilo,hmaxsoilo,hminsoilo)
       CALL hdfwrt2di(vegtyp1,nx1,ny1,sd_id1,0,0,'vegtyp',             &
                         'Vegetation type','index')

       deallocate(itmpsoilo,stat=istat)
       deallocate(hmaxsoilo,stat=istat)
       deallocate(hminsoilo,stat=istat)
       deallocate(soiltyp1,stat=istat)
       deallocate(stypfrct1,stat=istat)
       deallocate(vegtyp,stat=istat)
       deallocate(vegtyp1,stat=istat)

! lai
       CALL hdfrd2d(sd_id,"lai",nx,ny,tem1,istat,itmp)
       IF (istat /= 0) GO TO 110
       if( intrp_opt == 0 ) then
         call copyarray(tem1,nx,ny,1, tem1o, nx1,ny1,1, &
                     nxbgn,nxend,nybgn,nyend,1,1)
       else
         call intrpxy3d(tem1,nx,1,nx-1,ny,1,ny-1,1,1,1, &
             wgtsx,isx,wgtsy,jsy,intrp_opt,tem1o,nx1,1,nx1-1,ny1,1,ny1-1, temx1yz)
       endif
       if( hdmpfmt /= 0 )  &
          CALL hdfwrt2d(tem1o,nx1,ny1,sd_id1,0,hdfcompr,            &
                  'lai','Leaf Area Index','index',itmpo)
! roufns
       CALL hdfrd2d(sd_id,"roufns",nx,ny,tem1,istat,itmp)
       IF (istat /= 0) GO TO 110
       if( intrp_opt == 0 ) then
         call copyarray(tem1,nx,ny,1, tem1o, nx1,ny1,1, &
                     nxbgn,nxend,nybgn,nyend,1,1)
       else
         call intrpxy3d(tem1,nx,1,nx-1,ny,1,ny-1,1,1,1, &
             wgtsx,isx,wgtsy,jsy,intrp_opt,tem1o,nx1,1,nx1-1,ny1,1,ny1-1, temx1yz)
       endif
       if( hdmpfmt /= 0 )  &
          CALL hdfwrt2d(tem1o,nx1,ny1,sd_id1,0,hdfcompr,            &
                  'roufns','Surface roughness','0-1',itmpo)
! veg
       CALL hdfrd2d(sd_id,"veg",nx,ny,tem1,istat,itmp)
       IF (istat /= 0) GO TO 110
       if( intrp_opt == 0 ) then
         call copyarray(tem1,nx,ny,1, tem1o, nx1,ny1,1, &
                     nxbgn,nxend,nybgn,nyend,1,1)
       else
         call intrpxy3d(tem1,nx,1,nx-1,ny,1,ny-1,1,1,1, &
             wgtsx,isx,wgtsy,jsy,intrp_opt,tem1o,nx1,1,nx1-1,ny1,1,ny1-1, temx1yz)
       endif
       if( hdmpfmt /= 0 )  &
          CALL hdfwrt2d(tem1o,nx1,ny1,sd_id1,0,hdfcompr,            &
                  'veg','Vegetation fraction','fraction',itmpo)

    END IF  ! landin ==1

  END IF  ! grdbas == 1
  IF( grdbas == 1 ) GO TO 930

  IF( varflg == 1 .or. varflg == 2 ) then

!-----------------------------------------------------------------------
!  Read in total values of variables from history dump
!-----------------------------------------------------------------------

    CALL hdfrd3d(sd_id,"u",nx,ny,nz,tem1,istat,itmp,hmax,hmin)
    IF (istat /= 0) GO TO 110
      if( intrp_opt == 0 ) then
      call copyarray(tem1,nx,ny,nz, tem1o, nx1,ny1,nz1, &
           nxbgn,nxend,nybgn,nyend,nzbgn,nzend)
      else
        call intrpxy3d(tem1,nx,1,nx,ny,1,ny-1,nz,nzbgn,nzend,  &
             wgtux,iux,wgtsy,jsy,intrp_opt,tem1o,nx1,1,nx1-1,ny1,1,ny1-1, temx1yz)
      endif
    if( savespace == 1 ) then
      tem1o(:,:,nz1)=0.0
      tem1o(:,:,1)=0.0
      tem1o(:,:,nz1-1)=0.0
    endif
    if( hdmpfmt /= 0 )  &
    CALL hdfwrt3d(tem1o,nx1,ny1,nz1,sd_id1,1,hdfcompr,            &
                  'u','u-velocity','m/s', itmpo,hmaxo,hmino)

    IF( exbcdmpopt /= 0 .and. ubcdmp == 1) THEN
      IF (exbcdmpopt == 4) THEN

        do k=1,min(nz1, rayklowest )
          tem1o((2+ngbrz):(nx1-2-ngbrz),(1+ngbrz):(ny1-1-ngbrz),k:k)=tem1o(2,2,k)
        enddo

      ENDIF
      CALL hdfwrt3d(tem1o,nx1,ny1,nz1,sd_id2,1,exbccompr,                  &
                'u','u-velocity','m/s',itmpo,hmaxo,hmino)
    END IF

    CALL hdfrd3d(sd_id,"v",nx,ny,nz,tem1,istat,itmp,hmax,hmin)
    IF (istat /= 0) GO TO 110
      if( intrp_opt == 0 ) then
      call copyarray(tem1,nx,ny,nz, tem1o, nx1,ny1,nz1, &
           nxbgn,nxend,nybgn,nyend,nzbgn,nzend)
      else
        call intrpxy3d(tem1,nx,1,nx-1,ny,1,ny,nz,nzbgn,nzend,  &
             wgtsx,isx,wgtvy,jvy,intrp_opt,tem1o,nx1,1,nx1-1,ny1,1,ny1-1, temx1yz)
      endif
    if( savespace == 1 ) then
      tem1o(:,:,nz1)=0.0
      tem1o(:,:,1)=0.0
      tem1o(:,:,nz1-1)=0.0
    endif
    if( hdmpfmt /= 0 )  &
    CALL hdfwrt3d(tem1o,nx1,ny1,nz1,sd_id1,1,hdfcompr,            &
                  'v','v-velocity','m/s', itmpo,hmaxo,hmino)

    IF( exbcdmpopt /= 0 .and. vbcdmp == 1) THEN
      IF (exbcdmpopt == 4) THEN
        do k=1,min(nz1,rayklowest)
          tem1o(1+ngbrz:nx1-1-ngbrz,2+ngbrz:ny1-2-ngbrz,k:k)=tem1o(2,2,k)
        enddo
      ENDIF
      CALL hdfwrt3d(tem1o,nx1,ny1,nz1,sd_id2,1,exbccompr,                  &
                'v','v-velocity','m/s',itmpo,hmaxo,hmino)
    END IF

    CALL hdfrd3d(sd_id,"w",nx,ny,nz,tem1,istat,itmp,hmax,hmin)
    IF (istat /= 0) GO TO 110
      if( intrp_opt == 0 ) then
        call copyarray(tem1,nx,ny,nz, tem1o, nx1,ny1,nz1, &
           nxbgn,nxend,nybgn,nyend,nzbgn,nzend)
      else
        call intrpxy3d(tem1,nx,1,nx-1,ny,1,ny-1,nz,nzbgn,nzend,  &
             wgtsx,isx,wgtsy,jsy,intrp_opt,tem1o,nx1,1,nx1-1,ny1,1,ny1-1, temx1yz)
      endif
    if( savespace == 1 ) then
      tem1o(:,:,  1)=0.0
      tem1o(:,:,nz1)=0.0
    endif
    if( hdmpfmt /= 0 )  &
    CALL hdfwrt3d(tem1o,nx1,ny1,nz1,sd_id1,1,hdfcompr,            &
                  'w','w-velocity','m/s', itmpo,hmaxo,hmino)

    IF( exbcdmpopt /= 0 .and. wbcdmp == 1) THEN
      IF (exbcdmpopt == 4) THEN
        do k=1,min(nz1,rayklowest)
          tem1o(1+ngbrz:nx1-1-ngbrz,1+ngbrz:ny1-1-ngbrz,k:k)=tem1o(2,2,k)
        enddo
      ENDIF
      CALL hdfwrt3d(tem1o,nx1,ny1,nz1,sd_id2,1,exbccompr,                  &
                'w','w-velocity','m/s',itmpo,hmaxo,hmino)
    END IF

  ENDIF

  IF( varflg == 1 .or. varflg == 3 ) then

    CALL hdfrd3d(sd_id,"pt",nx,ny,nz,tem1,istat,itmp,hmax,hmin)

      DO k=1,nz-1
        ptmin = tem1(1,1,k)
        ptmax = tem1(1,1,k)
        DO i=1,nx-1
        DO j=1,ny-1
          ptmin=min(ptmin,tem1(i,j,k))
          ptmax=max(ptmax,tem1(i,j,k))
        ENDDO
        ENDDO
        print*,'k, ptmin, ptmax=', k, ptmin, ptmax
      enddo


    IF (istat /= 0) GO TO 110
      if( intrp_opt == 0 ) then
      call copyarray(tem1,nx,ny,nz, tem1o, nx1,ny1,nz1, &
           nxbgn,nxend,nybgn,nyend,nzbgn,nzend)
      else
        call intrpxy3d(tem1,nx,1,nx-1,ny,1,ny-1,nz,nzbgn,nzend,  &
             wgtsx,isx,wgtsy,jsy,intrp_opt,tem1o,nx1,1,nx1-1,ny1,1,ny1-1, temx1yz)
      endif

    if( savespace == 1 ) then
      tem1o(:,:,nz1)=0.0
      tem1o(:,:,1)=0.0
      tem1o(:,:,nz1-1)=0.0
    endif
    if( hdmpfmt /= 0 )  &
    CALL hdfwrt3d(tem1o,nx1,ny1,nz1,sd_id1,1,hdfcompr,            &
                  'pt','Potential temperature','K', itmpo,hmaxo,hmino)

    IF( exbcdmpopt /= 0 .and. ptbcdmp == 1) THEN
      IF (exbcdmpopt == 4) THEN
        do k=1,min(nz1,rayklowest)
          tem1o(1+ngbrz:nx1-1-ngbrz,1+ngbrz:ny1-1-ngbrz,k:k)=tem1o(2,2,k)
        enddo
      ENDIF
      CALL hdfwrt3d(tem1o,nx1,ny1,nz1,sd_id2,1,exbccompr,                  &
                  'pt','Potential temperature','K', itmpo,hmaxo,hmino)
    END IF

    CALL hdfrd3d(sd_id,"p",nx,ny,nz,tem1,istat,itmp,hmax,hmin)

    IF (istat /= 0) GO TO 110
      if( intrp_opt == 0 ) then
        call copyarray(tem1,nx,ny,nz, tem1o, nx1,ny1,nz1, &
             nxbgn,nxend,nybgn,nyend,nzbgn,nzend)
      else
        call intrpxy3d(tem1,nx,1,nx-1,ny,1,ny-1,nz,nzbgn,nzend,  &
             wgtsx,isx,wgtsy,jsy,intrp_opt,tem1o,nx1,1,nx1-1,ny1,1,ny1-1, temx1yz)
      endif
    if( savespace == 1 ) then
      tem1o(:,:,nz1)=0.0
      tem1o(:,:,1)=0.0
      tem1o(:,:,nz1-1)=0.0
    endif
    if( hdmpfmt /= 0 )  &
    CALL hdfwrt3d(tem1o,nx1,ny1,nz1,sd_id1,1,hdfcompr,            &
                  'p','Pressure','Pascal', itmpo,hmaxo,hmino)

    IF( exbcdmpopt /= 0 .and. prbcdmp == 1) THEN
      IF (exbcdmpopt == 4) THEN
        do k=1,min(nz1,rayklowest)
          tem1o(1+ngbrz:nx1-1-ngbrz,1+ngbrz:ny1-1-ngbrz,k:k)=tem1o(2,2,k)
        enddo
      ENDIF
      CALL hdfwrt3d(tem1o,nx1,ny1,nz1,sd_id2,1,exbccompr,                  &
                  'p','Pressure','Pascal', itmpo,hmaxo,hmino)
    END IF
  END IF

!-----------------------------------------------------------------------
!  Read in moisture variables
!-----------------------------------------------------------------------

  IF( mstin == 1 .AND. mstout == 1) THEN

    CALL hdfrd3d(sd_id,"qv",nx,ny,nz,tem1,istat,itmp,hmax,hmin)
    IF (istat /= 0) GO TO 110
      if( intrp_opt == 0 ) then
        call copyarray(tem1,nx,ny,nz, tem1o, nx1,ny1,nz1, &
           nxbgn,nxend,nybgn,nyend,nzbgn,nzend)
      else
        call intrpxy3d(tem1,nx,1,nx-1,ny,1,ny-1,nz,nzbgn,nzend,  &
             wgtsx,isx,wgtsy,jsy,intrp_opt,tem1o,nx1,1,nx1-1,ny1,1,ny1-1, temx1yz)
      endif
    if( savespace == 1 ) then
      tem1o(:,:,nz1)=0.0
      tem1o(:,:,1)=0.0
      tem1o(:,:,nz1-1)=0.0
    endif
    if( hdmpfmt /= 0 )  &
    CALL hdfwrt3d(tem1o,nx1,ny1,nz1,sd_id1,1,hdfcompr,            &
         'qv','Water vapor specific humidity','kg/kg', itmpo,hmaxo,hmino)

    IF( exbcdmpopt /= 0 .and. qvbcdmp == 1) THEN
      IF (exbcdmpopt == 4) THEN
        do k=1,min(nz1,rayklowest)
          tem1o(1+ngbrz:nx1-1-ngbrz,1+ngbrz:ny1-1-ngbrz,k:k)=tem1o(2,2,k)
        enddo
      ENDIF
      CALL hdfwrt3d(tem1o,nx1,ny1,nz1,sd_id2,1,exbccompr,                  &
         'qv','Water vapor specific humidity','kg/kg', itmpo,hmaxo,hmino)
    END IF

    IF (intver < intver530) THEN

      IF (P_QC > 0 )   &
        CALL hdfrd3d(sd_id,"qc",nx,ny,nz,tem1,istat,itmp,hmax,hmin)
      IF (istat /= 0) GO TO 110
      if( intrp_opt == 0 ) then
      call copyarray(tem1,nx,ny,nz, tem1o, nx1,ny1,nz1, &
           nxbgn,nxend,nybgn,nyend,nzbgn,nzend)
        else
          call intrpxy3d(tem1,nx,1,nx-1,ny,1,ny-1,nz,nzbgn,nzend,  &
               wgtsx,isx,wgtsy,jsy,intrp_opt,tem1o,nx1,1,nx1-1,ny1,1,ny1-1, temx1yz)
        endif
      if( savespace == 1 ) then
        tem1o(:,:,nz1)=0.0
        tem1o(:,:,1)=0.0
        tem1o(:,:,nz1-1)=0.0
      endif
      if( hdmpfmt /= 0 )  &
      CALL hdfwrt3d(tem1o,nx1,ny1,nz1,sd_id1,1,hdfcompr,            &
           'qc','Cloud water mixing ratio','kg/kg',itmpo,hmaxo,hmino)

      IF( exbcdmpopt /= 0 .and. qcbcdmp == 1) THEN
        IF (exbcdmpopt == 4) THEN
          do k=1,min(nz1,rayklowest)
            tem1o(1+ngbrz:nx1-1-ngbrz,1+ngbrz:ny1-1-ngbrz,k:k)=tem1o(2,2,k)
          enddo
        ENDIF
        CALL hdfwrt3d(tem1o,nx1,ny1,nz1,sd_id2,1,exbccompr,                  &
            'qc','Cloud water mixing ratio','kg/kg',itmpo,hmaxo,hmino)
      END IF


      IF (P_QR > 0 )   &
        CALL hdfrd3d(sd_id,"qr",nx,ny,nz,tem1,istat,itmp,hmax,hmin)
      IF (istat /= 0) GO TO 110
      if( intrp_opt == 0 ) then
      call copyarray(tem1,nx,ny,nz, tem1o, nx1,ny1,nz1, &
           nxbgn,nxend,nybgn,nyend,nzbgn,nzend)
      else
        call intrpxy3d(tem1,nx,1,nx-1,ny,1,ny-1,nz,nzbgn,nzend,  &
             wgtsx,isx,wgtsy,jsy,intrp_opt,tem1o,nx1,1,nx1-1,ny1,1,ny1-1, temx1yz)
      endif
      if( savespace == 1 ) then
        tem1o(:,:,nz1)=0.0
        tem1o(:,:,1)=0.0
        tem1o(:,:,nz1-1)=0.0
      endif
      if( hdmpfmt /= 0 )  &
      CALL hdfwrt3d(tem1o,nx1,ny1,nz1,sd_id1,1,hdfcompr,            &
           'qr','Rain water mixing ratio','kg/kg',itmpo,hmaxo,hmino)

      IF( exbcdmpopt /= 0 .and. qrbcdmp == 1) THEN
        IF (exbcdmpopt == 4) THEN
          do k=1,min(nz1,rayklowest)
            tem1o(1+ngbrz:nx1-1-ngbrz,1+ngbrz:ny1-1-ngbrz,k:k)=tem1o(2,2,k)
          enddo
        ENDIF
        CALL hdfwrt3d(tem1o,nx1,ny1,nz1,sd_id2,1,exbccompr,                  &
             'qr','Rain water mixing ratio','kg/kg',itmpo,hmaxo,hmino)
      END IF

    ELSE

      DO nq = 1,nscalar
        IF (nqscalarin(nq) > 0 )  THEN
          CALL hdfrd3d(sd_id,qnames(nq),nx,ny,nz,tem1,istat,itmp,hmax,hmin)
          IF (istat /= 0) GO TO 110
          if( intrp_opt == 0 ) then
            call copyarray(tem1,nx,ny,nz, tem1o, nx1,ny1,nz1, &
                    nxbgn,nxend,nybgn,nyend,nzbgn,nzend)
          else
            call intrpxy3d(tem1,nx,1,nx-1,ny,1,ny-1,nz,nzbgn,nzend,  &
                    wgtsx,isx,wgtsy,jsy,intrp_opt,tem1o,nx1,1,nx1-1,ny1,1,ny1-1, temx1yz)
          endif
          if( savespace == 1 ) then
            tem1o(:,:,nz1)=0.0
            tem1o(:,:,1)=0.0
            tem1o(:,:,nz1-1)=0.0
          endif

          ! DTD: Turn off bit-packing for Z-array

          IF(nq >= 13 .and. hdfcompr >= 4) THEN
            hdfcomprtmp = hdfcompr - 3
            exbccomprtmp = exbccompr - 3
          ELSE
            hdfcomprtmp = hdfcompr
            exbccomprtmp = exbccompr
          END IF

          if( hdmpfmt /= 0 )  &
            CALL hdfwrt3d(tem1o,nx1,ny1,nz1,sd_id1,1,hdfcomprtmp,            &
                     TRIM(qnames(nq)),TRIM(qdescp(nq)),'kg/kg',itmpo,hmaxo,hmino)

          IF( exbcdmpopt /= 0 .and. qrbcdmp == 1) THEN
            IF (exbcdmpopt == 4) THEN
              do k=1,min(nz1,rayklowest)
                tem1o(1+ngbrz:nx1-1-ngbrz,1+ngbrz:ny1-1-ngbrz,k:k)=tem1o(2,2,k)
              enddo
            ENDIF
            CALL hdfwrt3d(tem1o,nx1,ny1,nz1,sd_id2,1,exbccomprtmp,                  &
                      TRIM(qnames(nq)),TRIM(qdescp(nq)),'kg/kg',itmpo,hmaxo,hmino)
          END IF

        END IF
      END DO

    END IF

    IF( rainin == 1 .AND. rainout == 1 .and. hdmpfmt /= 0) THEN
      CALL hdfrd2d(sd_id,"raing",nx,ny,tem1,istat,itmp)
      IF (istat /= 0) GO TO 110
      if( intrp_opt == 0 ) then
      call copyarray(tem1,nx,ny,1, tem1o, nx1,ny1,1, &
           nxbgn,nxend,nybgn,nyend,1,1)
      else
        call intrpxy3d(tem1,nx,1,nx-1,ny,1,ny-1,1,1,1, &
             wgtsx,isx,wgtsy,jsy,intrp_opt, tem1o,nx1,1,nx1-1,ny1,1,ny1-1, temx1yz)
      endif
      if( hdmpfmt /= 0 )  &
      CALL hdfwrt2d(tem1o,nx1,ny1,sd_id1,0,hdfcompr,            &
                  'raing','Grid supersaturation rain','mm',itmpo)

      CALL hdfrd2d(sd_id,"rainc",nx,ny,tem1,istat,itmp)
      IF (istat /= 0) GO TO 110
      if( intrp_opt == 0 ) then
      call copyarray(tem1,nx,ny,1, tem1o, nx1,ny1,1, &
           nxbgn,nxend,nybgn,nyend,1,1)
      else
        call intrpxy3d(tem1,nx,1,nx-1,ny,1,ny-1,1,1,1, &
             wgtsx,isx,wgtsy,jsy,intrp_opt, tem1o,nx1,1,nx1-1,ny1,1,ny1-1, temx1yz)
      endif
      if( hdmpfmt /= 0 )  &
      CALL hdfwrt2d(tem1o,nx1,ny1,sd_id1,0,hdfcompr,            &
                  'rainc','Cumulus convective rain','mm',itmpo)
    END IF

    IF( prcin == 1 .AND. prcout ==1 .and. hdmpfmt /= 0) THEN
      CALL hdfrd2d(sd_id,"prcrate1",nx,ny,tem1,istat,itmp)
      IF (istat /= 0) GO TO 110
      if( intrp_opt == 0 ) then
      call copyarray(tem1,nx,ny,1, tem1o, nx1,ny1,1, &
           nxbgn,nxend,nybgn,nyend,1,1)
      else
        call intrpxy3d(tem1,nx,1,nx-1,ny,1,ny-1,1,1,1, &
             wgtsx,isx,wgtsy,jsy,intrp_opt, tem1o,nx1,1,nx1-1,ny1,1,ny1-1, temx1yz)
      endif
      if( hdmpfmt /= 0 )  &
      CALL hdfwrt2d(tem1o,nx1,ny1,sd_id1,0,hdfcompr,            &
          'prcrate1','Total precip. rate','kg/(m**2*s)',itmpo)

      CALL hdfrd2d(sd_id,"prcrate2",nx,ny,tem1,istat,itmp)
      IF (istat /= 0) GO TO 110
      if( intrp_opt == 0 ) then
      call copyarray(tem1,nx,ny,1, tem1o, nx1,ny1,1, &
           nxbgn,nxend,nybgn,nyend,1,1)
      else
        call intrpxy3d(tem1,nx,1,nx-1,ny,1,ny-1,1,1,1, &
             wgtsx,isx,wgtsy,jsy,intrp_opt, tem1o,nx1,1,nx1-1,ny1,1,ny1-1, temx1yz)
      endif
      if( hdmpfmt /= 0 )  &
      CALL hdfwrt2d(tem1o,nx1,ny1,sd_id1,0,hdfcompr,            &
          'prcrate2','Grid scale precip. rate','kg/(m**2*s)',itmpo)

      CALL hdfrd2d(sd_id,"prcrate3",nx,ny,tem1,istat,itmp)
      IF (istat /= 0) GO TO 110
      if( intrp_opt == 0 ) then
      call copyarray(tem1,nx,ny,1, tem1o, nx1,ny1,1, &
           nxbgn,nxend,nybgn,nyend,1,1)
      else
        call intrpxy3d(tem1,nx,1,nx-1,ny,1,ny-1,1,1,1, &
             wgtsx,isx,wgtsy,jsy,intrp_opt, tem1o,nx1,1,nx1-1,ny1,1,ny1-1, temx1yz)
      endif
      if( hdmpfmt /= 0 )  &
      CALL hdfwrt2d(tem1o,nx1,ny1,sd_id1,0,hdfcompr,            &
          'prcrate3','Cumulative precip. rate','kg/(m**2*s)',itmpo)

      CALL hdfrd2d(sd_id,"prcrate4",nx,ny,tem1,istat,itmp)
      IF (istat /= 0) GO TO 110
      if( intrp_opt == 0 ) then
      call copyarray(tem1,nx,ny,1, tem1o, nx1,ny1,1, &
           nxbgn,nxend,nybgn,nyend,1,1)
      else
        call intrpxy3d(tem1,nx,1,nx-1,ny,1,ny-1,1,1,1, &
             wgtsx,isx,wgtsy,jsy,intrp_opt, tem1o,nx1,1,nx1-1,ny1,1,ny1-1, temx1yz)
      endif
      if( hdmpfmt /= 0 )  &
      CALL hdfwrt2d(tem1o,nx1,ny1,sd_id1,0,hdfcompr,            &
          'prcrate4','Microphysics precip. rate','kg/(m**2*s)',itmpo)
    END IF

  END IF

  IF( tkein == 1 .AND. tkeout == 1 .and. hdmpfmt /= 0) THEN

    CALL hdfrd3d(sd_id,"tke",nx,ny,nz,tem1,istat,itmp,hmax,hmin)
    IF (istat /= 0) GO TO 110
    if( intrp_opt == 0 ) then
      call copyarray(tem1,nx,ny,nz, tem1o, nx1,ny1,nz1, &
           nxbgn,nxend,nybgn,nyend,nzbgn,nzend)
      else
        call intrpxy3d(tem1,nx,1,nx-1,ny,1,ny-1,1,1,1, &
             wgtsx,isx,wgtsy,jsy,intrp_opt, tem1o,nx1,1,nx1-1,ny1,1,ny1-1, temx1yz)
      endif
      if( savespace == 1 ) then
        tem1o(:,:,nz1)=0.0
        tem1o(:,:,1)=0.0
        tem1o(:,:,nz1-1)=0.0
      endif
      if( hdmpfmt /= 0 )  &
      CALL hdfwrt3d(tem1o,nx1,ny1,nz1,sd_id1,1,hdfcompr,            &
                  'tke','Turbulent Kinetic Energy','(m/s)**2',          &
                     itmpo,hmaxo,hmino)

  END IF

  IF( trbin == 1 .AND. trbout == 1 .and. hdmpfmt /= 0 ) THEN

    CALL hdfrd3d(sd_id,"kmh",nx,ny,nz,tem1,istat,itmp,hmax,hmin)
    IF (istat /= 0) GO TO 110
    if( intrp_opt == 0 ) then
    call copyarray(tem1,nx,ny,nz, tem1o, nx1,ny1,nz1, &
         nxbgn,nxend,nybgn,nyend,nzbgn,nzend)
      else
        call intrpxy3d(tem1,nx,1,nx-1,ny,1,ny-1,1,1,1, &
             wgtsx,isx,wgtsy,jsy,intrp_opt, tem1o,nx1,1,nx1-1,ny1,1,ny1-1, temx1yz)
      endif
    if( savespace == 1 ) then
      tem1o(:,:,nz1)=0.0
      tem1o(:,:,1)=0.0
      tem1o(:,:,nz1-1)=0.0
    endif
    if( hdmpfmt /= 0 )  &
    CALL hdfwrt3d(tem1o,nx1,ny1,nz1,sd_id1,1,hdfcompr,            &
        'kmh','Hori. turb. mixing coef. for momentum','m**2/s',         &
                   itmpo,hmaxo,hmino)

    CALL hdfrd3d(sd_id,"kmv",nx,ny,nz,tem1,istat,itmp,hmax,hmin)
    IF (istat /= 0) GO TO 110
    if( intrp_opt == 0 ) then
    call copyarray(tem1,nx,ny,nz, tem1o, nx1,ny1,nz1, &
         nxbgn,nxend,nybgn,nyend,nzbgn,nzend)
    else
      call intrpxy3d(tem1,nx,1,nx-1,ny,1,ny-1,1,1,1, &
           wgtsx,isx,wgtsy,jsy,intrp_opt, tem1o,nx1,1,nx1-1,ny1,1,ny1-1, temx1yz)
    endif
    if( savespace == 1 ) then
      tem1o(:,:,nz1)=0.0
      tem1o(:,:,1)=0.0
      tem1o(:,:,nz1-1)=0.0
    endif
    if( hdmpfmt /= 0 )  &
    CALL hdfwrt3d(tem1o,nx1,ny1,nz1,sd_id1,1,hdfcompr,            &
        'kmv','Vert. turb. mixing coef. for momentum','m**2/s',         &
                   itmpo,hmaxo,hmino)
  END IF

!
!-----------------------------------------------------------------------
!
!  If sfcout = 1, write out the surface variables,
!  tsoil, qsoil, and wetcanp.
!  working here
!-----------------------------------------------------------------------

  IF( sfcin == 1 .AND. sfcout == 1 .and. hdmpfmt /= 0 ) THEN

      ALLOCATE (itmpsoil4(nx,ny,nzsoil,nstyp+1),stat=istat)
      IF (istat /= 0) THEN
        WRITE (6,*) "HDFDUMP: ERROR allocating itmpsoil, returning"
        RETURN
      END IF
      ALLOCATE (hmaxsoil(nzsoil),stat=istat)
      IF (istat /= 0) THEN
        WRITE (6,*) "HDFDUMP: ERROR allocating hmaxsoil, returning"
        RETURN
      END IF
      ALLOCATE (hminsoil(nzsoil),stat=istat)
      IF (istat /= 0) THEN
        WRITE (6,*) "HDFDUMP: ERROR allocating hminisoil, returning"
        RETURN
      END IF
      allocate(itmpsoilo4(nx1,ny1,nzsoil,nstyp1+1),stat=istat)
      IF (istat /= 0) THEN
          WRITE (6,*) "HDFDUMP: ERROR allocating hminsoilo, returning"
          RETURN
      END IF
      allocate(hmaxsoilo(nzsoil),stat=istat)
      IF (istat /= 0) THEN
          WRITE (6,*) "HDFDUMP: ERROR allocating hminsoilo, returning"
          RETURN
      END IF
      allocate(hminsoilo(nzsoil),stat=istat)
      IF (istat /= 0) THEN
          WRITE (6,*) "HDFDUMP: ERROR allocating hminsoilo, returning"
          RETURN
      END IF

      ALLOCATE (rsoil(nx,ny,nzsoil,nstyp+1),stat=istat)
      IF (istat /= 0) THEN
        WRITE (6,*) "HDFDUMP: ERROR allocating rsoil, returning"
        RETURN
      END IF
      ALLOCATE (rsoil1(nx1,ny1,nzsoil,nstyp1+1),stat=istat)
      IF (istat /= 0) THEN
        WRITE (6,*) "HDFDUMP: ERROR allocating rsoil1, returning"
        RETURN
      END IF

      CALL hdfrd4d(sd_id,"tsoil",nx,ny,nzsoil,nstyp1+1,rsoil,istat,&
                 itmpsoil4,hmaxsoil,hminsoil)

      if( intrp_opt == 0 ) then
         call copyarray4(rsoil,nx,ny,nstyp,rsoil1,nx1,ny1,nstyp1,nzsoil,  &
                        nxbgn,nxend,nybgn,nyend,1,nstyp1)
      else
        call intrp_soil_real(nx,ny,nx1,ny1,nstyp,nstyp1,nzsoil,wgtsx,wgtsy,isx,jsy,     &
                      rsoil,soiltyp,stypfrct,rsoil1)

      endif

      CALL hdfwrt4d(rsoil1,nx1,ny1,nzsoil,nstyp1+1,sd_id1,0,       &
                  hdfcompr,'tsoil','Soil temperature','K',              &
                  itmpsoilo4,hmaxsoilo,hminsoilo)
!
      CALL hdfrd4d(sd_id,"qsoil",nx,ny,nzsoil,nstyp+1,rsoil,istat,&
                 itmpsoil4,hmaxsoil,hminsoil)

      if( intrp_opt == 0 ) then
         call copyarray4(rsoil,nx,ny,nstyp,rsoil1,nx1,ny1,nstyp1,nzsoil,  &
                        nxbgn,nxend,nybgn,nyend,1,nstyp1)
      else
        call intrp_soil_real(nx,ny,nx1,ny1,nstyp,nstyp1,nzsoil,wgtsx,wgtsy,isx,jsy,     &
                      rsoil,soiltyp,stypfrct,rsoil1)

      endif

      CALL hdfwrt4d(rsoil1,nx1,ny1,nzsoil,nstyp1+1,sd_id1,0,       &
                  hdfcompr,'qsoil','Soil moisture','fraction',   &
                  itmpsoilo4,hmaxsoilo,hminsoilo)

      DEALLOCATE(itmpsoil4)
      DEALLOCATE(hmaxsoil)
      DEALLOCATE(hminsoil)
      DEALLOCATE(itmpsoilo4)
      DEALLOCATE(hmaxsoilo)
      DEALLOCATE(hminsoilo)
      DEALLOCATE(rsoil)
      DEALLOCATE(rsoil1)

      ALLOCATE (itmpsoil(nx,ny,nstyp+1),stat=istat)
      IF (istat /= 0) THEN
        WRITE (6,*) "HDFDUMP: ERROR allocating itmpsoil, returning"
        RETURN
      END IF
      ALLOCATE (hmaxsoil(nstyp+1),stat=istat)
      IF (istat /= 0) THEN
        WRITE (6,*) "HDFDUMP: ERROR allocating hmaxsoil, returning"
        RETURN
      END IF
      ALLOCATE (hminsoil(nstyp+1),stat=istat)
      IF (istat /= 0) THEN
        WRITE (6,*) "HDFDUMP: ERROR allocating hminisoil, returning"
        RETURN
      END IF

      allocate(itmpsoilo(nx1,ny1,nstyp1+1),stat=istat)
      IF (istat /= 0) THEN
          WRITE (6,*) "HDFDUMP: ERROR allocating hminsoilo, returning"
          RETURN
      END IF
      allocate(hmaxsoilo(nstyp1+1),stat=istat)
      IF (istat /= 0) THEN
          WRITE (6,*) "HDFDUMP: ERROR allocating hminsoilo, returning"
          RETURN
      END IF
      allocate(hminsoilo(nstyp1+1),stat=istat)
      IF (istat /= 0) THEN
          WRITE (6,*) "HDFDUMP: ERROR allocating hminsoilo, returning"
          RETURN
      END IF

      ALLOCATE (wetcsoil(nx,ny,nstyp+1),stat=istat)
      IF (istat /= 0) THEN
        WRITE (6,*) "HDFDUMP: ERROR allocating wetcsoil, returning"
        RETURN
      END IF
      ALLOCATE (wetcsoil1(nx1,ny1,nstyp1+1),stat=istat)
      IF (istat /= 0) THEN
        WRITE (6,*) "HDFDUMP: ERROR allocating wetcsoil1, returning"
        RETURN
      END IF

      CALL hdfrd3d(sd_id,"wetcanp",nx,ny,nstyp+1,wetcsoil,istat,        &
                   itmpsoil,hmaxsoil,hminsoil)
      if( intrp_opt == 0 ) then
         call copyarray4(wetcsoil,nx,ny,nstyp,wetcsoil1,nx1,ny1,nstyp1,1,  &
                        nxbgn,nxend,nybgn,nyend,1,nstyp1)
      else
        call intrp_soil_real(nx,ny,nx1,ny1,nstyp,nstyp1,1,wgtsx,wgtsy,isx,jsy,     &
                      wetcsoil,soiltyp,stypfrct,wetcsoil1)
      endif

      CALL hdfwrt3d(wetcsoil1,nx1,ny1,nstyp1+1,sd_id1,0,hdfcompr,        &
                 'wetcanp','Canopy water amount','fraction',            &
                 itmpsoilo,hmaxsoilo,hminsoilo)

      DEALLOCATE(itmpsoil)
      DEALLOCATE(hmaxsoil)
      DEALLOCATE(hminsoil)
      DEALLOCATE(itmpsoilo)
      DEALLOCATE(hmaxsoilo)
      DEALLOCATE(hminsoilo)
      DEALLOCATE(wetcsoil)
      DEALLOCATE(wetcsoil1)

    IF (snowout == 1) THEN

         CALL hdfrd2d(sd_id,"snowdpth",nx,ny,tem1,istat,itmp)
         IF (istat /= 0) GO TO 110
         if( intrp_opt == 0 ) then
           call copyarray(tem1,nx,ny,1, tem1o, nx1,ny1,1, &
                     nxbgn,nxend,nybgn,nyend,1,1)
         else
           call intrpxy3d(tem1,nx,1,nx-1,ny,1,ny-1,1,1,1, &
               wgtsx,isx,wgtsy,jsy,intrp_opt,tem1o,nx1,1,nx1-1,ny1,1,ny1-1, temx1yz)
         endif
         if( hdmpfmt /= 0 )  &
            CALL hdfwrt2d(tem1o,nx1,ny1,sd_id1,0,hdfcompr,            &
                    'snowdpth','Snow depth','m',itmpo)

    END IF

  END IF

!-----------------------------------------------------------------------
!
!  Friendly exit message
!
!-----------------------------------------------------------------------

  930   CONTINUE

  WRITE(6,'(/a,F8.1,a/)')                                               &
  ' Data at time=', time/60,' (min) were successfully read.'

  ireturn = 0

  GO TO 130

!-----------------------------------------------------------------------
!
!  Error during read
!
!-----------------------------------------------------------------------

  110   CONTINUE
  WRITE(6,'(/a/)') ' Error reading data in HDFREAD'
  ireturn=1

  GO TO 130

!-----------------------------------------------------------------------
!
!  End-of-file during read
!
!-----------------------------------------------------------------------

!  120   CONTINUE
!  WRITE(6,'(/a/)') ' End of file reached in HDFREAD'
!  ireturn=2

  130   CONTINUE

  CALL hdfclose(sd_id,istat)

  IF (ireturn == 0) THEN
    IF (istat == 0) THEN
      WRITE(6,'(/a/a)') &
      "HDFREADDUMP: Successfully read ", trim(filename)
    ELSE
      WRITE(6,'(/a,i3,a/,a)') &
      "HDFREADDUMP: ERROR (status=", istat, ") closing ", trim(filename)
    END IF
  END IF

  if( hdmpfmt /= 0 ) then
    CALL hdfclose(sd_id1,istat)
    IF (istat == 0) THEN
      WRITE(6,'(/a/,a,a)') &
      "HDFREADDUMP: Successfully wrote ", filename_out(1:lfilename_out)
    endif
  endif

  IF (exbcdmpopt /= 0) THEN
    CALL hdfclose(sd_id2,istat)
  ENDIF

! DEALLOCATE (itmp,stat=istat)
  DEALLOCATE (hmax,stat=istat)
  DEALLOCATE (hmin,stat=istat)

  IF( hdfcompr > 3 .or. exbccompr > 3) then
! DEALLOCATE (itmpo,stat=istat)
  DEALLOCATE (hmaxo,stat=istat)
  DEALLOCATE (hmino,stat=istat)
  ENDIF

  RETURN
END SUBROUTINE hdfreaddump

!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE GET_GRIDINFO_FROM_HDF         ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE get_gridinfo_from_hdf(filename,nx,ny,nz,x,y,z,ireturn)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Read in grid dimensions from base state/grid history data.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  7/17/2000.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    filename Channel number for binary reading.
!
!  OUTPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: stat, sd_id
  CHARACTER (LEN=*) :: filename

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  REAL :: x(nx),y(ny),z(nz)

  INTEGER (KIND=selected_int_kind(4)), ALLOCATABLE :: itmp(:,:,:) ! Temporary array
  REAL, ALLOCATABLE :: hmax(:), hmin(:) ! Temporary array

  INTEGER :: ireturn           ! Return status indicator

  INTEGER istat

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL hdfopen(filename,1,sd_id)

  IF (sd_id < 0) THEN
    WRITE (6,*) "get_gridinfo_from_hdf: ERROR opening ",                 &
                 trim(filename)," for reading."
    GO TO 110
  ELSE
    WRITE(6,*) 'File ',filename,' openned.'
  END IF

  ALLOCATE (itmp(nx,ny,nz),stat=istat)
  ALLOCATE (hmax(nz),stat=istat)
  ALLOCATE (hmin(nz),stat=istat)

! print*,'sd_id, nx =', sd_id, nx

  CALL hdfrd1d(sd_id,"x",nx,x,istat)

! print*,'istat after reading x =', istat

  IF (istat /= 0) GO TO 110
  CALL hdfrd1d(sd_id,"y",ny,y,istat)
  IF (istat /= 0) GO TO 110
  CALL hdfrd1d(sd_id,"z",nz,z,istat)
  IF (istat /= 0) GO TO 110

  ireturn = 0
  GO TO 130

!-----------------------------------------------------------------------
!
!  Error during read
!
!-----------------------------------------------------------------------

  110   CONTINUE
  WRITE(6,'(/a/)') ' Error reading data in GET_GRIDINFO_FROM_HDF.'
  ireturn=1

130 CONTINUE

!tmp  stat = sfendacc(sd_id)   ! is this necessary?
  CALL hdfclose(sd_id,stat)
  DEALLOCATE (itmp)
  DEALLOCATE (hmax)
  DEALLOCATE (hmin)

  RETURN
END SUBROUTINE get_gridinfo_from_hdf

SUBROUTINE copyarray(arrayin,nx,ny,nz, arrayout, nx1,ny1,nz1, &
           nxbgn,nxend,nybgn,nyend,nzbgn,nzend)

INTEGER :: nx,ny,nz, nx1,ny1,nz1,nxbgn,nxend,nybgn,nyend,nzbgn,nzend
REAL :: arrayin(nx,ny,nz)
REAL :: arrayout(nx1,ny1,nz1)

  DO k=1,nz1
  DO j=1,ny1
  DO i=1,nx1
    arrayout(i,j,k)=arrayin(nxbgn+i-1,nybgn+j-1,nzbgn+k-1)
  ENDDO
  ENDDO
  ENDDO

return

END SUBROUTINE copyarray

SUBROUTINE copyarray4(arrayin,nx,ny,nz, arrayout, nx1,ny1,nz1,nzsoil, &
           nxbgn,nxend,nybgn,nyend,nzbgn,nzend)

INTEGER :: nx,ny,nz, nx1,ny1,nz1,nxbgn,nxend,nybgn,nyend,nzbgn,nzend
INTEGER :: nzsoil
REAL :: arrayin(nx,ny,nzsoil,nz)
REAL :: arrayout(nx1,ny1,nzsoil,nz1)

  DO k=1,nz1
  DO n=1,nzsoil
  DO j=1,ny1
  DO i=1,nx1
    arrayout(i,j,n,k)=arrayin(nxbgn+i-1,nybgn+j-1,n,nzbgn+k-1)
  ENDDO
  ENDDO
  ENDDO
  ENDDO

return

END SUBROUTINE copyarray4


SUBROUTINE copyarrayi(arrayin,nx,ny,nz, arrayout, nx1,ny1,nz1, &
           nxbgn,nxend,nybgn,nyend,nzbgn,nzend)

INTEGER :: nx,ny,nz, nx1,ny1,nz1,nxbgn,nxend,nybgn,nyend,nzbgn,nzend
INTEGER :: arrayin(nx,ny,nz)
INTEGER :: arrayout(nx1,ny1,nz1)

  DO k=1,nz1
  DO j=1,ny1
  DO i=1,nx1
    arrayout(i,j,k)=arrayin(nxbgn+i-1,nybgn+j-1,nzbgn+k-1)
  ENDDO
  ENDDO
  ENDDO

return

END SUBROUTINE copyarrayi
