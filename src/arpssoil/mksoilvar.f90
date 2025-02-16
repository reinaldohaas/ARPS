!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE MKSOILVAR                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE mksoilvar(nx,ny,nz,nzsoil,nstyps,                            &
           mxdays,mxstns,                                               &
           xs,ys,                                                       &
           hterain,zsc,                                                 &
           ptbar,ptprt,                                                 &
           pbar,pprt,                                                   &
           zpsoil,tsoil,qsoil,wetcanp,snowdpth,                         &
           soiltyp,stypfrct,                                            &
           precfile,staid,                                              &
           iptstn,jptstn,iptapi0,jptapi0,                               &
           api0,obsprec,difprec,totpreca,prec,                          &
           itema,initapi,api1,api2,kk2dep,                              &
           totwt,dprec,tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the surface variables for ARPS model.
!  Includes calculation of API from precip data to initialize
!  soil moisture.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: John Mewes (precip and API) and Keith Brewster (temps)
!  March, 1997
!
!  MODIFICATION HISTORY:
!
!  04/11/1997 (Keith Brewster)
!  Added processing of NCEP real-time precip data to API calculation.
!  at the same time delete the table data array veg(14).
!
!  04/15/1997 (Keith Brewster)
!  Added option of reading in pre-calculated soil moisture data,
!  since the real-time precip data are typically only once per day.
!  Minor reorganization of parameters. Bells and whistles added.
!  Renamed to version 2.0
!
!  01/05/1998 (Donghai Wang)
!  Fixed a problem in the case of k1< 0.0.
!
!  12/09/1998 (Donghai Wang)
!  Added the snow cover.
!
!  2000/01/03 (Gene Bassett)
!  Renamed mstinit to soilinit2 and added capability to update soil
!  data over water and snow cover with data in a soil data file.
!
!  05/25/2002  (J. Brotzge)
!  Added additional soil variables for new soil scheme.
!
!  1 June 2002 Eric Kemp
!  Soil variable updates.
!
!  09/2004 (J.Case, ENSCO Inc.)
!  Included provisions for interpolating RUC SSTs to the ARPS
!  grid at water grid points.
!
!-----------------------------------------------------------------------
!
!  Set various parameters used in API calculation
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations....
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz,nstyps
  INTEGER :: nzsoil 
  INTEGER :: mxdays,mxstns
!
  REAL :: xs    (nx)        ! The x-coord. of the physical and
                            ! computational grid. Defined at u-point.
  REAL :: ys    (ny)        ! The y-coord. of the physical and
                            ! computational grid. Defined at v-point.

  REAL :: ptprt  (nx,ny,nz)    ! Perturbation potential temperature (K)
  REAL :: pprt   (nx,ny,nz)    ! Perturbation pressure (Pascal)
                               ! humidity (kg/kg)
  REAL :: ptbar  (nx,ny,nz)    ! Base state potential temperature (K)
  REAL :: pbar   (nx,ny,nz)    ! Base state pressure (Pascal)

  INTEGER :: soiltyp(nx,ny,nstyps)    ! Soil type
  REAL :: stypfrct  (nx,ny,nstyps)

  INTEGER :: soilinit2         ! Moisture initialization option
  INTEGER :: prdata            ! Precip data option
  CHARACTER (LEN=256) :: apifile      ! Name of file containing names of the
                                      ! precip data files.
  CHARACTER (LEN=256) :: apiinit      ! File containing initial API's and
                                      !their corresponding locations
  CHARACTER (LEN=256) :: prcpdir      ! Disk directory containing NCEP precip files
  CHARACTER (LEN=256) :: prcplst      ! File containing list of NCEP precip stns

! J.Case, ENSCO Inc, 9/16/2004 -- Changed to include hour in API dates.
  CHARACTER (LEN=13) :: inidate      ! Day for start of precip files
  CHARACTER (LEN=13) :: fnldate      ! Day for stop of precip files

  REAL :: respapi              ! Desired response of the initial API
!                              analysis to wavelengths equal to
!                              2 x (mean initAPI stn sep).
  INTEGER :: ndays             ! Number of days between initial date
!                              and history file date
  REAL :: k1,k2                ! API depletion minimum coefficients for
!                              the ground surface layer and the root
!                              zone layer
  REAL :: range                ! Range parameter in Barnes' wt. fcn.
!                              denominator
  REAL :: respprec             ! Desired response of the final daily
!                              precip analyses at wavelengths equal
!                              to 2 x (mean stn sep).
  REAL :: gamma(2)             ! Range multipliers for 1st and 2nd
!                              passes
!
!-----------------------------------------------------------------------
!
!  API generation arrays:
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=80) :: precfile(mxstns)
  CHARACTER (LEN=7)  :: staid(mxstns)
  INTEGER :: iptstn(mxstns)
  INTEGER :: jptstn(mxstns)
  INTEGER :: iptapi0(mxstns)
  INTEGER :: jptapi0(mxstns)
  REAL :: api0(mxstns)
  REAL :: obsprec(mxdays,mxstns)
  REAL :: difprec(mxdays,mxstns)
  REAL :: totpreca(mxstns)
  REAL :: prec(nx,ny)
!
!-----------------------------------------------------------------------
!
!  Computed variables
!
!-----------------------------------------------------------------------
!
  REAL :: zpsoil(nx,ny,nzsoil)          ! Soil layer height (m)
  REAL :: tsoil (nx,ny,nzsoil,0:nstyps) ! Soil temperature (K)
  REAL :: qsoil (nx,ny,nzsoil,0:nstyps) ! Soil moisture (m**3/m**3)
  REAL :: wetcanp(nx,ny,0:nstyps)       ! Canopy water amount
  REAL :: snowdpth(nx,ny)               ! Snow depth (m)

  REAL :: hterain(nx,ny)
  REAL :: zsc(nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  Work Arrays
!
!-----------------------------------------------------------------------
!
  INTEGER :: itema(nx,ny)
  REAL :: initapi(nx,ny,nstyps)
  REAL :: api1(nx,ny,nstyps)
  REAL :: api2(nx,ny,nstyps)
  REAL :: kk2dep(nx,ny)
  REAL :: totwt(nx,ny)
  REAL :: dprec(nx,ny)
  REAL :: tem1(nx,ny,nz)

  ! temparary arrays for soilinit2 option 5:
  REAL, ALLOCATABLE :: tsoil2 (:,:,:,:)
  REAL, ALLOCATABLE :: qsoil2 (:,:,:,:)
  REAL, ALLOCATABLE :: wetcanp2(:,:,:)
  INTEGER, ALLOCATABLE :: soiltyp2(:,:,:)
!
!-----------------------------------------------------------------------
!
!  API parameters
!
!-----------------------------------------------------------------------
!
  COMMON /apicom1/ apiinit,apifile,prcpdir,prcplst,inidate,fnldate
  COMMON /apicom2/ soilinit2,prdata
  COMMON /apicom3/ respapi,k1,k2,range,respprec,gamma
!
!-----------------------------------------------------------------------
!
!  Map projection variables
!
!-----------------------------------------------------------------------
!
  REAL :: latnot(2)
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'soilcst.inc'
  INCLUDE 'indtflg.inc'
!
!-----------------------------------------------------------------------
!
!  Miscellaneous local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: sec24hr
  PARAMETER( sec24hr = 86400 )
!
  INTEGER :: i,j,ii,ij,is,is2,istn,jstn,iday,ix,iy,ireturn
  INTEGER :: ipt,jpt,out
  INTEGER :: tstns,totapi0,minday,maxday,length

! J.Case, ENSCO Inc., 9/16/2004 -- Changed variable name for times extracted from fnldate.
  INTEGER :: iniyr,inimo,inidy,inihr,initi,nowti,jlday
  INTEGER :: fnlyr,fnlmo,fnldy,fnlhr,fnlmin,fnlsec
  INTEGER :: itime,iyr,imon,idy,ihr,imin,isec
  INTEGER :: idel,jdel,istr,iend,jstr,jend
  INTEGER :: i1,i2,i3,i4,i5,j1,j2,j3,j4,j5
!
  REAL :: tmp2,tmp3,prs2,prs3,w2,w3,ptk
  REAL :: p0inv,tmpk,pres
  REAL :: xpt,ypt,xctr,yctr,xll,yll
  REAL :: lat,lon,ddx,ddy
  REAL :: wmin,wmax,kk1,kk2,adj1,adj2
  REAL :: apirng,apir2,r2,dist,const,thrdst
  REAL :: sumwt,sumapi,wgt,precip,maxprec
  REAL :: bias,rms,wsfcscl,wdpscl
  REAL :: mindist,totdist,total
  REAL :: d_0,d_1,d_1star,depth1
!
!  INTEGER :: tsfcout,tsoilout,wsfcout,wdpout,wcanpout,snowdout
!  INTEGER :: tsfcin,tsoilin,wsfcin,wdpin,wcanpin,snowdin
!  COMMON /intgcom/tsfcout,tsoilout,wsfcout,wdpout,wcanpout,snowdout

  INTEGER :: zpsoilout,tsoilout,qsoilout,wcanpout,snowdout
  INTEGER :: zpsoilin,tsoilin,qsoilin,wcanpin,snowdin
  COMMON /intgcom/ zpsoilin,tsoilout,qsoilout,wcanpout,snowdout

  REAL :: twater
  COMMON /realcom/twater

  REAL :: amax, amin

  INTEGER :: nbwater,nbw_max
  PARAMETER (nbw_max=128)
  REAL :: tbwater(nbw_max),blat1(nbw_max),blat2(nbw_max)
  REAL :: blon1(nbw_max),blon2(nbw_max)
  COMMON /bwcomi/ nbwater
  COMMON /bwcomr/ tbwater,blat1,blat2,blon1,blon2

  INTEGER :: ib
  INTEGER :: k

!
! J.Manobianco (3/15/2002) Add capability to read/set SST
! J.Manobianco/J.Case (2/6/2003)  Add capability to read RUC SST, soil moisture,
!               and soil temperature on LCC RUCH grid
!
! J.Case
! Other variables
!
  CHARACTER (LEN=10) :: infilename
  INTEGER :: nx_ext,ny_ext,ios
  INTEGER :: iscl(nx,ny),jscl(nx,ny)
  INTEGER :: iproj_ext
  INTEGER :: nisst, njsst, isst, jsst

!  40-km NAM model

  PARAMETER (nx_ext=185,ny_ext=129)

  REAL :: lat_arps(nx,ny),lon_arps(nx,ny)
  REAL :: lat_ruc(nx_ext,ny_ext),lon_ruc(nx_ext,ny_ext)
  REAL :: xs2d(nx,ny),ys2d(nx,ny)
  REAL :: dxfld(nx_ext),dyfld(ny_ext)
  REAL :: rdxfld(nx_ext),rdyfld(ny_ext)
  REAL :: slopey(nx_ext,ny_ext)
  REAL :: alphay(nx_ext,ny_ext)
  REAL :: betay(nx_ext,ny_ext)
  REAL :: blats, blons, res, sstlat, sstlon
  REAL :: sltk_sfc(nx_ext,ny_ext)
  REAL :: sltk_deep(nx_ext,ny_ext)
  REAL :: soim_sfc(nx_ext,ny_ext)
  REAL :: soim_deep(nx_ext,ny_ext)

  REAL :: lat_ext(nx_ext,ny_ext),lon_ext(nx_ext,ny_ext)
  REAL :: x0,y0,x0_ext,y0_ext
  REAL :: x_ext(nx_ext),y_ext(ny_ext)
  REAL :: swlat_ext, swlon_ext
  REAL :: dx_ext, dy_ext

  REAL :: sltk_sfc_interp(nx,ny)
  REAL :: sltk_deep_interp(nx,ny)
  REAL :: soim_sfc_interp(nx,ny)
  REAL :: soim_deep_interp(nx,ny)
  REAL :: goes_sst(nx,ny)
  REAL :: scale_ext, trlon, latnot_ext(2), trlon_ext

  LOGICAL :: sst_flag

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (soilinit2 == 5) THEN
    ALLOCATE(tsoil2  (nx,ny,nzsoil,0:nstyps))
    ALLOCATE(qsoil2  (nx,ny,nzsoil,0:nstyps))
    ALLOCATE(wetcanp2(nx,ny,0:nstyps))
    ALLOCATE(soiltyp2(nx,ny,nstyps))
  END IF

! J.Case, ENSCO Inc., 9/16/2004 -- corrected this bug in assigning API date to soil output.
!  as a temporary fix, read year, month, day, AND hour from the arpssoil.input file.
!
  READ(inidate,'(i4,1x,i2,1x,i2,1x,i2)') iniyr,inimo,inidy,inihr
  CALL ctim2abss(iniyr,inimo,inidy,inihr,0,0, initi)
  READ(fnldate,'(i4,1x,i2,1x,i2,1x,i2)') fnlyr,fnlmo,fnldy,fnlhr
  fnlmin = 0
  fnlsec = 0
  curtim = 0

  CALL ctim2abss(fnlyr, fnlmo, fnldy, fnlhr, fnlmin, fnlsec, nowti)
  nowti=nowti+curtim

  ndays=((nowti-initi)/sec24hr)+1

  WRITE(6,'(a,i6)') ' Number of days of data for API : ',ndays

  IF(ndays > mxdays.AND.(soilinit2 == 2 .OR. soilinit2 == 3)) THEN
    WRITE(6,'(a,i6,/,a)')                                               &
        ' Number of days of data for API exceeds mxdays dimension: ',   &
        mxdays,'   Increase mxdays or change inidate.'
    STOP
  END IF

!
!-----------------------------------------------------------------------
!
!  Find the center lat,lon and the location (in map coord.) of the
!  lower left corner (xll,yll).
!
!-----------------------------------------------------------------------
!
  latnot(1)=trulat1
  latnot(2)=trulat2
  CALL setmapr(mapproj,sclfct,latnot,trulon)
  CALL lltoxy(1,1,ctrlat,ctrlon,xctr,yctr)
  WRITE(6,*) ' dx= ',dx,' dy= ',dy
  xll=xctr-(0.5*(nx-3)*dx)
  yll=yctr-(0.5*(ny-3)*dy)
!
  WRITE(6,'(2x,a,f10.2,a,f10.2//)')                                     &
      ' Grid center:  lat= ',ctrlat,' lon= ',ctrlon

!
  IF( sfcdat == 1 ) THEN

    DO j=1,ny
      DO i=1,nx
        soiltyp(i,j,1) = styp
        stypfrct(i,j,1) = 1.0
      END DO
    END DO

    DO is=2,nstyps
      DO j=1,ny
        DO i=1,nx
          soiltyp (i,j,is) = 0
          stypfrct(i,j,is) = 0.0
        END DO
      END DO
    END DO

  ELSE IF (sfcdat == 2 .OR. (sfcdat == 3.AND.landin /= 1) ) THEN

    DO is=1,nstyps
      DO j=1,ny
        DO i=1,nx
          soiltyp (i,j,is) = 0
          stypfrct(i,j,is) = 0.0
        END DO
      END DO
    END DO

    CALL readsfcdt(nx,ny,nstyps,sfcdtfl,dx,dy,                          &
             mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,       &
                   soiltyp,stypfrct,itema,tem1,tem1,tem1,tem1)
    WRITE(6,*) ' nstyp = ',nstyp

    IF( nstyp == 1 ) THEN

      DO j=1,ny-1
        DO i=1,nx-1
          stypfrct(i,j,1) = 1.0
        END DO
      END DO

    END IF

  ELSE IF (sfcdat == 3 .AND. landin == 1) THEN

    WRITE(6,'(1x,a/)')                                                  &
        'Surface property data in the history file was used.'

  ELSE

    WRITE(6,'(1x,a,i3,a/)')                                             &
        'Invalid surface data input option. sfcdat =',sfcdat,           &
        '. Program stopped in MKSOILVAR.'
    STOP

  END IF

  DO is=0,nstyp

    PRINT*,'In MKSOILVAR for is =', is

    CALL a3dmax0(tsoil(1,1,1,is),                                       &
         1,nx,1,nx-1,1,ny,1,ny-1,1,nzsoil,1,nzsoil, amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'tsoilmin= ', amin,', tsoilmax =',amax

    CALL a3dmax0(qsoil(1,1,1,is),                                       &
         1,nx,1,nx-1,1,ny,1,ny-1,1,nzsoil,1,nzsoil, amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'qsoilmin= ', amin,', qsoilmax =',amax


  END DO


  !
  ! ENSCO inc. J. Case begin
  !
  sst_flag = .true.

  IF(sst_flag) THEN

!
! J.Manobianco (3/15/02) Add setorig call to use xytoll for SST interpolation

    CALL setorig (1,xll,yll)

!   Read in external data

    infilename='sltk.5cm'
    CALL RDEXTSST (nx_ext,ny_ext,infilename,                            &
                   sltk_sfc,blats,blons,res,nisst,njsst,sst_flag)
    print *,'SLTK5 VARS after RDEXT',nisst,njsst,sst_flag
    print *,sltk_sfc(1,1),sltk_sfc(nx_ext,ny_ext),                      &
            sltk_sfc(1,ny_ext),sltk_sfc(nx_ext,1)

    infilename='sltk.40cm'
    CALL RDEXTSST (nx_ext,ny_ext,infilename,                            &
                   sltk_deep,blats,blons,res,nisst,njsst,sst_flag)
    print *,'SLTK40 VARS after RDEXT',nisst,njsst,sst_flag
    print *,sltk_deep(1,1),sltk_deep(nx_ext,ny_ext),                    &
            sltk_deep(1,ny_ext),sltk_deep(nx_ext,1)

!    infilename='soim.1cm'
!    CALL RDEXTSST (nx_ext,ny_ext,infilename,                            &
!                  soim_sfc,blats,blons,res,nisst,njsst,sst_flag)
!    print *,'SOIM1 VARS after RDEXT',nisst,njsst,sst_flag
!    print *,soim_sfc(1,1),soim_sfc(nx_ext,ny_ext),                      &
!           soim_sfc(1,ny_ext),soim_sfc(nx_ext,1)
!
!    infilename='soim.100cm'
!    CALL RDEXTSST (nx_ext,ny_ext,infilename,                            &
!                 soim_deep,blats,blons,res,nisst,njsst,sst_flag)
!    print *,'SOIM100 VARS after RDEXT',nisst,njsst,sst_flag
!    print *,soim_deep(1,1),soim_deep(nx_ext,ny_ext),                    &
!           soim_deep(1,ny_ext),soim_deep(nx_ext,1)

! J.Case, ENSCO Inc. (2/11/2005) -- Reading in GOES SST data from FIT.

     goes_sst=-9999.0
     print *,'Reading GOES SST Product.'
     infilename='sst.data'
     OPEN (unit=45,file=infilename,form='FORMATTED',                     &
           status='OLD',iostat=ios,err=499)
     DO i=1,nx
     DO j=1,ny
       read (45,'(29x,f9.2)') goes_sst(i,j)
       if (goes_sst(i,j) > 0) then
         goes_sst(i,j)=goes_sst(i,j)+273.15
       endif
     ENDDO
     ENDDO
  
     print *,'Normal read of GOES SST data.'
     GOTO 498 

499  print *,'Error reading GOES SST file. IOSTAT = ',ios
498  CONTINUE
                                                                                                                                   
!
! =================================================================
!
!       Compute latitude and longitude of ARPS grid so that x and y
!       locations can be computed relative to the external grid.
!
                                                                                                                                   
    CALL xytoll (nx,ny,xs,ys,lat_arps,lon_arps)
    print *,lat_arps(1,1),lon_arps(1,1),                                 &
           lat_arps(nx-1,ny-1),lon_arps(nx-1,ny-1)
                                                                                                                                   
! J.Case, ENSCO Inc. (9/22/2004)
! Set up external grid information
!
                                                                                                                                   
    iproj_ext = 2        ! LCC
    scale_ext = 1.
    trlon_ext = 265.     ! NAM-40 grid projection
    latnot_ext(1) = 25.  ! NAM-40 grid projection
    latnot_ext(2) = 25.  ! NAM-40 grid projection

!  -- NAM 40 km model
    swlat_ext = 12.190   ! NAM-40 grid projection
    swlon_ext = -133.459 ! NAM-40 grid projection
    dx_ext = 40.63525E+03    ! NAM-40 grid projection
    dy_ext = 40.63525E+03    ! NAM-40 grid projection
                                                                                                                                   
! J.Case, ENSCO Inc. (9/22/2004)
! Set map projection to RUC.
!
                                                                                                                                   
    call setmapr (iproj_ext,scale_ext,latnot_ext,trlon_ext)
    call lltoxy (1,1,swlat_ext,swlon_ext,x0_ext,y0_ext)
                                                                                                                                   
! J.Case, ENSCO Inc.
! x,y locations of NAM:
!
                                                                                                                                   
    x_ext(1)=x0_ext
    y_ext(1)=y0_ext
                                                                                                                                   
    do i=2,nx_ext
      x_ext(i)=x_ext(i-1) + dx_ext
    enddo
    do j=2,ny_ext
      y_ext(j)=y_ext(j-1) + dy_ext
    enddo
                                                                                                                                   
! J.Case, ENSCO Inc.
! lat,lon of external grid:
!
                                                                                                                                   
    call xytoll(nx_ext,ny_ext,x_ext,y_ext,lat_ruc,lon_ruc)
                                                                                                                                   
! J.Case, cont.
! Find x,y locations of external grid.
!
                                                                                                                                   
    call setmapr(iproj_ext,scale_ext,latnot_ext,trlon_ext)
    call setorig(1,x0_ext,y0_ext)
    DO j=1,ny_ext
      call lltoxy(1,1,lat_ruc(1,j),lon_ruc(1,j),                         &
                        x_ext(1),y_ext(j))
    ENDDO
    DO i=1,nx_ext
      call lltoxy(1,1,lat_ruc(i,1),lon_ruc(i,1),                         &
                        x_ext(i),y_ext(1))
    ENDDO
                                                                                                                                   
! J.Case, cont.
!  Find x,y locations of ARPS grid points in terms of external grid.
!
                                                                                                                                   
    CALL lltoxy   (nx,ny,lat_arps,lon_arps,xs2d,ys2d)
    CALL setijloc (nx,ny,nx_ext,ny_ext,xs2d,ys2d,                        &
                     x_ext,y_ext,iscl,jscl)
                                                                                                                                   
!    write (22,*) xs2d
!    write (23,*) ys2d
!    write (24,*) iscl
!    write (25,*) jscl
!    write (26,*) lat_arps
!    write (27,*) lon_arps
                                                                                                                                   
! J.Case, cont.
! Restore map projection to ARPS.
                                                                                                                                   
    CALL setmapr(mapproj,sclfct,latnot,trulon)
    CALL setorig (1,xll,yll)
                                                                                                                                   
! J.Case, cont.
! Interpolate external grid data to ARPS grid points.
                                                                                                                                   
    IF (sst_flag) THEN
      CALL setdxdy(nx_ext,ny_ext,1,nx_ext,1,ny_ext,                      &
                   x_ext,y_ext,dxfld,dyfld,rdxfld,rdyfld)
! J.Case, cont.
! Interpolate shallow soil temperature.
!
      CALL fldint2d(nx,ny,nx_ext,ny_ext,                                 &
                    1,nx,1,ny,1,nx_ext,1,ny_ext,                         &
                    2,xs2d,ys2d,sltk_sfc,                                &
                    x_ext,y_ext,iscl,jscl,                               &
                    dxfld,dyfld,rdxfld,rdyfld,                           &
                    slopey,alphay,betay,sltk_sfc_interp)
! J.Case, cont.
! Interpolate deep soil temperature.
!
      CALL fldint2d(nx,ny,nx_ext,ny_ext,                                 &
                    1,nx,1,ny,1,nx_ext,1,ny_ext,                         &
                    2,xs2d,ys2d,sltk_deep,                               &
                    x_ext,y_ext,iscl,jscl,                               &
                    dxfld,dyfld,rdxfld,rdyfld,                           &
                    slopey,alphay,betay,sltk_deep_interp)
! J.Case, cont.
! Interpolate 1-cm soil moisture.
!
!      CALL fldint2d(nx,ny,nx_ext,ny_ext,                                 &
!                    1,nx,1,ny,1,nx_ext,1,ny_ext,                         &
!                    2,xs2d,ys2d,soim_sfc,                                &
!                    x_ext,y_ext,iscl,jscl,                               &
!                    dxfld,dyfld,rdxfld,rdyfld,                           &
!                    slopey,alphay,betay,soim_sfc_interp)
! J.Case, cont.
! Interpolate 100-cm soil moisture.
!
!      CALL fldint2d(nx,ny,nx_ext,ny_ext,                                 &
!                    1,nx,1,ny,1,nx_ext,1,ny_ext,                         &
!                    2,xs2d,ys2d,soim_deep,                               &
!                    x_ext,y_ext,iscl,jscl,                               &
!                    dxfld,dyfld,rdxfld,rdyfld,                           &
!                    slopey,alphay,betay,soim_deep_interp)                  
        ENDIF   ! sst_flag
                                                                                                                                   
!    write (20,*) sltk_sfc
!    write (21,*) sltk_sfc_interp
                                                                                                                                   
! J.Case, cont.
!  Restore origin back to North Pole (0,0).
!
                                                                                                                                   
    CALL setorig (1,0.,0.)
                                                                                                                               
  ENDIF  ! sst_flag .eq. .true.
                                                                                                                                   
! =================================================================
!

  IF(soilinit2 /= 5) THEN
    p0inv=1./p0
    DO j=1,ny-1
      DO i=1,nx-1
        tmp3 = ptbar(i,j,3) + ptprt(i,j,3)
        prs3 = ALOG(pbar (i,j,3) + pprt (i,j,3))
        tmp2 = ptbar(i,j,2) + ptprt(i,j,2)
        prs2 = ALOG(pbar (i,j,2) + pprt (i,j,2))

        w2=1.+(zsc(i,j,2)-hterain(i,j))/(zsc(i,j,3)-zsc(i,j,2))
        w3=1.-w2

        ptk=w2*tmp2 + w3*tmp3
        pres=EXP((w2*prs2 + w3*prs3))

        tmpk = ptk * (p0inv*pres)**rddcp

        IF (nbwater > 0) CALL xytoll(1,1,xll+dx*(i+0.5),yll+dy*(j+0.5),lat,lon)
!
        DO is=1,nstyp
          IF(soiltyp(i,j,is) /= 0) THEN
            IF(  soiltyp(i,j,is) == 13 ) THEN
              DO k=1,nzsoil 
! J.Case, ENSCO (9/22/2004)
                IF (sst_flag) THEN                
! J.Case, ENSCO (2/11/2005) -- FIT SST product
                  IF (goes_sst(i,j) > 0.0) THEN
                    print *,'Using GOES SST at grid point I/J = ',i,j
                    print *,'Original = ',tsoil(i,j,k,is)
                    print *,'GOES SST = ',goes_sst(i,j)
                    print *,'External = ',sltk_deep_interp(i,j)
                    tsoil (i,j,k,is) = goes_sst(i,j)
                  ELSE
                    print *,'Using external model SST at grid point I/J = ',i,j
                    print *,'Original = ',tsoil(i,j,k,is)
                    print *,'External = ',sltk_deep_interp(i,j)
                    tsoil (i,j,k,is) = sltk_deep_interp(i,j)
                  END IF
                ELSE
                  tsoil (i,j,k,is) = twater + 273.15
                END IF
              END DO 
              DO ib = 1,nbwater
                IF ((lat >= blat1(ib)) .AND. (lat <= blat2(ib))         &
                      .AND. (lon >= blon1(ib)) .AND. (lon <= blon2(ib))) THEN
                  DO k=1,nzsoil 
                    tsoil  (i,j,k,is) = tbwater(ib) + 273.15
                  END DO 
                END IF
              END DO
            ELSE
!TODO:  Replace tsrat and t2rat with array of values for each soil level.
              tsoil (i,j,1,is) = tmpk + ttprt
              DO k=2,nzsoil 
                tsoil  (i,j,k,is) = tmpk + tbprt
              END DO 
            END IF

          END IF
        END DO
      END DO
    END DO

  END IF

  IF ( snowin /= 1 ) THEN
    DO j=1,ny
      DO i=1,nx
        snowdpth(i,j)=snowdpth0
      END DO
    END DO
  END IF

!
!-----------------------------------------------------------------------
!
!  Modify or initialize soil properties using soilinit2 option
!
!-----------------------------------------------------------------------
!
  IF( soilinit2 == 1 ) THEN
!
    DO is=1,nstyp
      DO j=1,ny-1
        DO i=1,nx-1
!
!TODO:  Replace wgrat and w2rat with array of values for each soil level.
          IF( soiltyp(i,j,is) /= 0 ) THEN
            qsoil (i,j,1,is) = wwlt(soiltyp(i,j,is)) + wgrat             &
                 *( wsat(soiltyp(i,j,is)) - wwlt(soiltyp(i,j,is)) )
            DO k=1,nzsoil 
              qsoil (i,j,k,is) = wwlt(soiltyp(i,j,is)) + w2rat           &
                 *( wsat(soiltyp(i,j,is)) - wwlt(soiltyp(i,j,is)) )
            END DO 
            wetcanp(i,j,is) = wetcanp0
          END IF

        END DO
      END DO
    END DO

!
  ELSE IF( soilinit2 == 2 .OR. soilinit2 == 3 ) THEN
!
!-----------------------------------------------------------------------
!
!  Read the initial API data and find its location in the grid
!
!-----------------------------------------------------------------------
!
    totapi0=0
    OPEN(UNIT=3,ERR=205,FILE=apiinit,STATUS='unknown')
!
    out=0
    WRITE(6,*)
    DO istn=1,mxstns
      jstn=istn-out
      READ(3,*,END=201) lat,lon,api0(jstn)
      CALL lltoxy(1,1,lat,lon,xpt,ypt)
      xpt=xpt-xll
      ypt=ypt-yll
      CALL findlc(nx,ny,xs,ys,xpt,ypt,ipt,jpt,ireturn)
      IF( ireturn == 0 ) THEN
        iptapi0(jstn)=ipt
        jptapi0(jstn)=jpt
        WRITE(6,*) 'API obs:',jstn,' found near ',ipt,jpt,              &
            ' with initial API of ',api0(jstn)
        WRITE(6,*) ' --> main soil type:',soiltyp(ipt,jpt,1),           &
                       ' wilting:',wwlt(soiltyp(ipt,jpt,1)),            &
                       ' saturation:',wsat(soiltyp(ipt,jpt,1))
        WRITE(6,*)
      ELSE
        WRITE(6,*)'Initial API station #',jstn,                         &
                 ' is outside the grid..'
        WRITE(6,*)
        out=out+1
      END IF
    END DO
    201   CONTINUE
    totapi0=istn-1-out
    CLOSE(3)
    205   CONTINUE
!
!-----------------------------------------------------------------------
!
!  Read the observed precip data and find its location in the grid
!
!-----------------------------------------------------------------------
!
    IF( prdata == 1 ) THEN

      out=0
      WRITE(6,*)
      OPEN(UNIT=1,ERR=301,FILE=apifile,STATUS='old')
      DO istn=1,mxstns
        jstn=istn-out
        READ(1,1010,END=301) precfile(jstn)
        1010       FORMAT(a80)
        length=LEN(precfile(jstn))
        CALL strlnth(precfile(jstn),length)
        WRITE(6,*) 'File:',precfile(jstn)(1:length)
        OPEN(UNIT=2,FILE=precfile(jstn),STATUS='old')
        READ(2,*) lat,lon
        CALL lltoxy(1,1,lat,lon,xpt,ypt)
        xpt=xpt-xll
        ypt=ypt-yll
        CALL findlc(nx,ny,xs,ys,xpt,ypt,ipt,jpt,ireturn)
        IF( ireturn == 0 ) THEN
!
          iptstn(jstn)=ipt
          jptstn(jstn)=jpt
          totpreca(jstn)=0.0
!
          WRITE(6,*)
          WRITE(6,*) 'Precip stn # ',jstn,' was found near i,j:',       &
              ipt,jpt
          IF( jpt >= nint(ny*0.5) ) THEN
            WRITE(6,*) 'North of Center'
          ELSE
            WRITE(6,*) 'South of Center'
          END IF
          IF( ipt >= nint(nx*0.5)) THEN
            WRITE(6,*) 'East of Center'
          ELSE
            WRITE(6,*) 'West of Center'
          END IF
!
          WRITE(6,*)
          WRITE(6,*)
          DO iday=1,ndays
            READ(2,*) obsprec(iday,jstn) ! in inches...
            obsprec(iday,jstn)=obsprec(iday,jstn)*2.54
!
!  Assume (arbitrarily) that the ground can only absorb 4.0 cm of
!  rainfall per day, with the rest being runoff.
!
            obsprec(iday,jstn)=AMIN1(obsprec(iday,jstn),4.0 )
          END DO
          WRITE(6,*)
        ELSE
          PRINT *,                                                      &
              'Precip station #',jstn,' is outside of the grid..'
          PRINT *
          out=out+1
        END IF
        CLOSE(2)
      END DO
      301     CONTINUE
      CLOSE(1)
      CLOSE(2)
      tstns=istn-1-out

    ELSE IF (prdata == 2) THEN

! J.Case, ENSCO Inc. (9/22/2004)

      CALL rdnceppr(nx,ny,mxstns,mxdays,xs,ys,xll,yll,                  &
                      fnlyr,fnlmo,fnldy,ndays,                          &
                      prcpdir,prcplst,                                  &
                      staid,iptstn,jptstn,obsprec,totpreca,             &
                      tstns)

    ELSE

      WRITE(6,'(i6,a)') prdata,' is not a valid prdata option!'
      STOP


    END IF
!
!---------------------------------------------------------------------
!
!  Analyze the initial root layer API to the grid.  Determine the
!  range parameter using respapi.
!
!---------------------------------------------------------------------
!
    i1=1
    i2=nint(nx*.25)
    i3=nint(nx*.50)
    i4=nint(nx*.75)
    i5=nx-1
!
    j1=1
    j2=nint(ny*.25)
    j3=nint(ny*.50)
    j4=nint(ny*.75)
    j5=ny-1
!
    WRITE(6,*)
    WRITE(6,*) 'Analyzing initial API to grid...'
    WRITE(6,*)
!
    IF( totapi0 > 1 ) THEN
!
      totdist=SQRT((nx*dx*0.001)*(ny*dy*0.001))*                        &
          (1.+SQRT(FLOAT(totapi0)))/FLOAT(totapi0-1) ! using separation that
!        would result if obs were randomly scattered
      apirng=                                                           &
          SQRT(-((2.*totdist/3.14159)**2)*LOG(respapi))
      apir2=apirng*apirng
      const=1.0E-06/apir2
!
      WRITE(6,*)
      WRITE(6,*) 'Using range:',apirng,' km, for initializing API'
      WRITE(6,*)

      DO iy=1,ny-1
        DO ix=1,nx-1
          sumwt=0.
          sumapi=0.
          DO istn=1,totapi0
            ddx=xs(iptapi0(istn))-xs(ix)
            ddy=ys(jptapi0(istn))-ys(iy)
            dist=(ddx*ddx+ddy*ddy)
            wgt=EXP(-dist*const)
            sumwt=sumwt+wgt
            sumapi=sumapi+wgt*api0(istn)
          END DO
          IF( sumwt > 0. ) sumapi=sumapi/sumwt

          DO is=1,nstyp
            IF( soiltyp(ix,iy,is) /= 0 ) initapi(ix,iy,is)=             &
                AMIN1(sumapi,(wsat(soiltyp(ix,iy,is))*d2*100.))
          END DO

        END DO
      END DO

    ELSE ! if no initial values available, set to 0.5*wsat

      DO is=1,nstyp
        DO iy=1,ny-1
          DO ix=1,nx-1
            IF( soiltyp(ix,iy,is) /= 0)                                 &
                initapi(ix,iy,is)=0.5*wsat(soiltyp(ix,iy,is))*d2*100.0
          END DO
        END DO
      END DO

    END IF
!
    WRITE(6,*)
    WRITE(6,*) 'Initial API analysis:'
    WRITE(6,*) '---------------------'
    WRITE(6,*)
!
    WRITE(6,1040) j5,initapi(i1,j5,1),initapi(i2,j5,1),                 &
                     initapi(i3,j5,1),initapi(i4,j5,1),                 &
                     initapi(i5,j5,1)

    WRITE(6,1040) j4,initapi(i1,j4,1),initapi(i2,j4,1),                 &
                     initapi(i3,j4,1),initapi(i4,j4,1),                 &
                     initapi(i5,j4,1)

    WRITE(6,1040) j3,initapi(i1,j3,1),initapi(i2,j3,1),                 &
                     initapi(i3,j3,1),initapi(i4,j3,1),                 &
                     initapi(i5,j3,1)

    WRITE(6,1040) j2,initapi(i1,j2,1),initapi(i2,j2,1),                 &
                     initapi(i3,j2,1),initapi(i4,j2,1),                 &
                     initapi(i5,j2,1)

    WRITE(6,1040) j1,initapi(i1,j1,1),initapi(i2,j1,1),                 &
                     initapi(i3,j1,1),initapi(i4,j1,1),                 &
                     initapi(i5,j1,1)

    WRITE(6,*)
    WRITE(6,1035) i1,i2,i3,i4,i5
    WRITE(6,*)
!
    minday=15
    depth1=0.01
!
    DO iy=1,ny-1
      DO ix=1,nx-1
        kk2dep(ix,iy)=1.0
      END DO
    END DO

    DO is=1,nstyp
      DO iy=1,ny-1
        DO ix=1,nx-1
          api1(ix,iy,is)=initapi(ix,iy,is)*depth1/d2
          api2(ix,iy,is)=initapi(ix,iy,is)
        END DO
      END DO
    END DO
!
    maxprec=0.0
    DO iday=1,ndays
      precip=0.0
      DO istn=1,tstns
        IF( obsprec(iday,istn) >= 0. ) precip=precip+obsprec(iday,istn)
      END DO
      precip=precip/tstns
      IF( precip > maxprec ) THEN
        maxprec=precip
        maxday=iday
      END IF
    END DO
!
!----------------------------------------------------------------------
!
!  If the range is not specified in the input file (i.e. = 0.0)
!  determine the own range based on the input desired response.
!
!----------------------------------------------------------------------
!
    totdist=0.0
    DO ii=1,tstns
      mindist=1.0E30
      DO ij=1,tstns
        IF( ii /= ij ) THEN
          ddx=(xs(iptstn(ii))-xs(iptstn(ij)))
          ddy=(ys(jptstn(ii))-ys(jptstn(ij)))
          mindist=AMIN1(mindist,(ddx*ddx+ddy*ddy))
        END IF
      END DO
      totdist=totdist+0.001*SQRT(mindist)
    END DO
!
    WRITE(6,*) 'Average station separation (km):',totdist/tstns
    WRITE(6,*)
!
    IF( range == 0. ) THEN
      DO ii=1,1000000
        d_0=.000001*ii
        d_1=d_0**(gamma(2)**1.0)
        d_1star=d_0+(1-d_0)*d_1
        IF( d_1star >= respprec ) EXIT
      END DO
      396     CONTINUE
!
      range=SQRT(-((2*(totdist/(tstns*1.0))/3.14159)**2)*LOG(d_0))
      WRITE(6,*)
      WRITE(6,*) 'Range not specified, using range:',range
      WRITE(6,1026) ' -- 1st pass response at 2 x (stn sep):',d_0
      WRITE(6,1026) ' -- 2nd pass response at 2 x (stn sep):',          &
           d_1star
      1026     FORMAT(a40,f10.8)
      WRITE(6,*)
    END IF
!
    r2=range*range
!
!----------------------------------------------------------------------
!
!  Begin main day loop
!  Objectively analyze the precipitation data on 'ndays' different
!  days using a 2-pass analysis.
!
!----------------------------------------------------------------------
!
    DO iday=1,ndays
!
      DO iy=1,ny-1
        DO ix=1,nx-1
          prec(ix,iy)=0.
          totwt(ix,iy)=0.
        END DO
      END DO
      itime=initi+(iday-1)*sec24hr
      CALL abss2ctim( itime,iyr,imon,idy,ihr,imin,isec )
!
!    1st pass....
!
      const=1.0E-06/(r2*gamma(1))
      thrdst=1.0E06*(6.*totdist/tstns)*(6.*totdist/tstns)
      idel=INT(SQRT(thrdst)/dx)+1
      jdel=INT(SQRT(thrdst)/dy)+1
!
      DO istn=1,tstns
        IF( obsprec(iday,istn) >= 0. ) THEN
          istr =MAX(1,(iptstn(istn)-idel))
          iend =MIN((nx-1),(iptstn(istn)+idel))
          jstr =MAX(1,(jptstn(istn)-jdel))
          jend =MIN((ny-1),(jptstn(istn)+jdel))
          DO iy=jstr,jend
            DO ix=istr,iend
              ddx=(xs(iptstn(istn))-xs(ix))
              ddy=(ys(jptstn(istn))-ys(iy))
              dist=(ddx*ddx+ddy*ddy)
              IF( dist < thrdst ) THEN
                wgt=EXP(-dist*const)
                prec(ix,iy)=prec(ix,iy)+wgt*obsprec(iday,istn)
                totwt(ix,iy)=totwt(ix,iy)+wgt
              END IF
            END DO
          END DO
        END IF
      END DO

      DO iy=1,ny-1
        DO ix=1,nx-1
          IF( totwt(ix,iy) > 0. ) prec(ix,iy)=prec(ix,iy)/totwt(ix,iy)
        END DO
      END DO
!
      IF( iday == maxday ) THEN
!
        WRITE(6,'(a,i4,a,i2.2,a,i2.2)')                                 &
                'Example precip analysis on ',                          &
                iyr,'-',imon,'-',idy
        WRITE(6,'(a)') '-------------------------------------'
        WRITE(6,*)
!
        1035       FORMAT(7X,i3,5X,i3,5X,i3,5X,i3,5X,i3)
        1040       FORMAT(i3,3X,f5.2,3X,f5.2,3X,f5.2,3X,f5.2,3X,f5.2)
        WRITE(6,1040) j5,prec(i1,j5),prec(i2,j5),                       &
                         prec(i3,j5),prec(i4,j5),                       &
                         prec(i5,j5)

        WRITE(6,1040) j4,prec(i1,j4),prec(i2,j4),                       &
                         prec(i3,j4),prec(i4,j4),                       &
                         prec(i5,j4)

        WRITE(6,1040) j3,prec(i1,j3),prec(i2,j3),                       &
                         prec(i3,j3),prec(i4,j3),                       &
                         prec(i5,j3)

        WRITE(6,1040) j2,prec(i1,j2),prec(i2,j2),                       &
                         prec(i3,j2),prec(i4,j2),                       &
                         prec(i5,j2)

        WRITE(6,1040) j1,prec(i1,j1),prec(i2,j1),                       &
                         prec(i3,j1),prec(i4,j1),                       &
                         prec(i5,j1)
        WRITE(6,*)
        WRITE(6,1035) i1,i2,i3,i4,i5
        WRITE(6,*)
!
        rms=0.0
        bias=0.0
        total=0.0
        DO istn=1,tstns
          IF( obsprec(iday,istn) >= 0. ) THEN
            1045           FORMAT(a4,1X,f6.3,1X,a6,1X,f6.3)
            WRITE(6,1045) 'obs:',obsprec(iday,istn),' anal:',           &
                    prec(iptstn(istn),jptstn(istn))
            rms=rms+(obsprec(iday,istn)-                                &
                    prec(iptstn(istn),jptstn(istn)))**2
            bias=bias+obsprec(iday,istn)-                               &
                prec(iptstn(istn),jptstn(istn))
            total=total+1.
          END IF
        END DO
        IF( total > 0. ) THEN
          rms=SQRT(rms/total)
          bias=-bias/total
          WRITE(6,*)
          WRITE(6,*) 'RMS Error:',rms,' Bias:',bias
          WRITE(6,*)
        END IF

      END IF
!
!  2nd pass....
!
      DO istn=1,tstns
        IF( obsprec(iday,istn) >= 0. ) THEN
          difprec(iday,istn)=obsprec(iday,istn)-                        &
                      prec(iptstn(istn),jptstn(istn))
        ELSE
          difprec(iday,istn)=-900.
        END IF
      END DO

      DO iy=1,ny-1
        DO ix=1,nx-1
          dprec(ix,iy)=0.
          totwt(ix,iy)=0.
        END DO
      END DO
!
      const=1.0E-06/(r2*gamma(2))
      thrdst=1.0E06*(2.5*totdist/tstns)*(2.5*totdist/tstns)
      idel=INT(SQRT(thrdst)/dx)+1
      jdel=INT(SQRT(thrdst)/dy)+1
      DO istn=1,tstns
        IF( difprec(iday,istn) > -100. ) THEN
          istr =MAX(1,(iptstn(istn)-idel))
          iend =MIN((nx-1),(iptstn(istn)+idel))
          jstr =MAX(1,(jptstn(istn)-jdel))
          jend =MIN((ny-1),(jptstn(istn)+jdel))
          DO iy=jstr,jend
            DO ix=istr,iend
              ddx=xs(iptstn(istn))-xs(ix)
              ddy=ys(jptstn(istn))-ys(iy)
              dist=(ddx*ddx+ddy*ddy)
              IF( dist < thrdst ) THEN
                wgt=EXP(-dist*const)
                dprec(ix,iy)=dprec(ix,iy)+wgt*difprec(iday,istn)
                totwt(ix,iy)=totwt(ix,iy)+wgt
              END IF
            END DO
          END DO
        END IF
      END DO
!
!
!  Sum difference and impose limits on daily precip values...
!
      DO iy=1,ny-1
        DO ix=1,nx-1
          IF( totwt(ix,iy) > 0. )                                       &
              prec(ix,iy)=prec(ix,iy)+(dprec(ix,iy)/totwt(ix,iy))
          prec(ix,iy)=AMAX1(prec(ix,iy),0.0)
          prec(ix,iy)=AMIN1(prec(ix,iy),4.0)
        END DO
      END DO
!
      IF( iday == maxday) THEN
!
        itime=initi+(iday-1)*sec24hr
        CALL abss2ctim( itime,iyr,imon,idy,ihr,imin,isec )
        WRITE(6,'(a,i4,a,i2.2,a,i2.2)')                                 &
                'Example precip analysis on ',                          &
                iyr,'-',imon,'-',idy
        WRITE(6,*) '-------------------------------------'
        WRITE(6,*)
!
        WRITE(6,1040) j5,prec(i1,j5),prec(i2,j5),                       &
                         prec(i3,j5),prec(i4,j5),                       &
                         prec(i5,j5)

        WRITE(6,1040) j4,prec(i1,j4),prec(i2,j4),                       &
                         prec(i3,j4),prec(i4,j4),                       &
                         prec(i5,j4)

        WRITE(6,1040) j3,prec(i1,j3),prec(i2,j3),                       &
                         prec(i3,j3),prec(i4,j3),                       &
                         prec(i5,j3)

        WRITE(6,1040) j2,prec(i1,j2),prec(i2,j2),                       &
                         prec(i3,j2),prec(i4,j2),                       &
                         prec(i5,j2)

        WRITE(6,1040) j1,prec(i1,j1),prec(i2,j1),                       &
                         prec(i3,j1),prec(i4,j1),                       &
                         prec(i5,j1)
        WRITE(6,*)
        WRITE(6,1035) i1,i2,i3,i4,i5
        WRITE(6,*)
        rms=0.0
        bias=0.0
        total=0.0
        DO istn=1,tstns
          IF(obsprec(iday,istn) >= 0.) THEN
            WRITE(6,1045) 'obs:',obsprec(iday,istn),' anal:',           &
                prec(iptstn(istn),jptstn(istn))
            rms=rms+(obsprec(iday,istn)-                                &
                prec(iptstn(istn),jptstn(istn)))**2
            bias=bias+obsprec(iday,istn)-                               &
                prec(iptstn(istn),jptstn(istn))
            total=total+1.
          END IF
        END DO
        IF( total > 0. ) THEN
          rms=SQRT(rms/total)
          bias=-bias/total
          WRITE(6,*)
          WRITE(6,*) 'RMS Error:',rms,' Bias:',bias
          WRITE(6,*)
        END IF
      END IF
!
!---------------------------------------------------------------------
!
!  Calculate the API - assume that API will be relatively
!  insensitive to the initial API after about 90 days.
!  Recall d2=100.0cm and depth1 (m) is set below.  Note that
!  depth1 is not necessarilly equal to d1 in soilcst.inc, which
!  is only a normalizing depth.  In Noilhan-Planton (1989) they
!  suggest that depth1 is only a few millimeters.
!
!  k1 and k2 are minimum depletion coefficients for the topsoil
!  and root layer.  kk1 and kk2 are the time-dependent depletion
!  coefficients where kk? is at a maximum on minday (for minimum
!  moisture depletion - set to January 15, or Julian day 15) and
!  kk? is at a minimum 6 months later (July 15)
!
!  The time-dependent depletion coefficients can THEN be modified
!  up or down according to the apparent volumetric soil water
!  content (i.e. the API) relative to wilting and saturation using
!  the adj1 and adj2 terms, the values of which can be derived
!  empirically using ARM data.  With the limited data available
!  at this time, the adj1 and adj2 terms were simply set to 1.0
!  due to inability to show value in adding this term.  However,
!  intuitively there should be some relationship.
!
!  The API depletion coefficient multiplies the difference
!  between the calculated API and the wilting point (or the
!  initial API, whichever is lower) instead of the multiplying
!  the prior API as in the common formulation.  The new method
!  was derived empirically and appeared to give better correlation
!  with observed data as well as a lower dependence on the initial
!  API using the depletion coefficient formula below.
!
!---------------------------------------------------------------------
!
!
      adj1=1.0
      adj2=1.0
      CALL julday( iyr, imon, idy, jlday )

      kk1=1.0-adj1*                                                     &
             (1.0-k1)*SIN((3.14159*(jlday-minday))                      &
             /365.0)
      kk2=1.0-adj2*                                                     &
             (1.0-k2)*SIN((3.14159*(jlday-minday))                      &
             /365.0)

      DO iy=1,ny-1
        DO ix=1,nx-1
          kk2dep(ix,iy)=kk2dep(ix,iy)*kk2
        END DO
      END DO

      DO is=1,nstyp
        DO iy=1,ny-1
          DO ix=1,nx-1
!
            IF(soiltyp(ix,iy,is) > 0) THEN

              IF( k1 > 0. ) THEN
                api1(ix,iy,is)=prec(ix,iy)+kk1*api1(ix,iy,is)
              ELSE
                api1(ix,iy,is)=0.0
              END IF
              api1(ix,iy,is)=                                           &
                  AMIN1(api1(ix,iy,is),                                 &
                      (wsat(soiltyp(ix,iy,is))*depth1*100.))
!
              api2(ix,iy,is)=prec(ix,iy)+                               &
                  kk2*(api2(ix,iy,is)-AMIN1(initapi(ix,iy,is),          &
                    wwlt(soiltyp(ix,iy,is))*(d2*100.)))+                &
                  AMIN1(initapi(ix,iy,is),                              &
                   wwlt(soiltyp(ix,iy,is))*(d2*100.))
              api2(ix,iy,is)=                                           &
                  AMIN1(api2(ix,iy,is),(wsat(soiltyp(ix,iy,is))*d2*100.))
!
!  Carry out 1/2 of the daily depletion on the final day (the
!  day of the run) without adding in a precip term.
!
              IF( iday == ndays) THEN
                IF( k1 > 0. )                                           &
                    api1(ix,iy,is)=(api1(ix,iy,is)+api1(ix,iy,is)*kk1)/2.0
                api2(ix,iy,is)=(api2(ix,iy,is)+api2(ix,iy,is)*kk2)/2.0
              END IF

            END IF
!
          END DO
        END DO
      END DO

      DO istn=1,tstns
!      write(6,1020) 'Stn #:',istn,'Day:',iyr,'-',imon,'-',idy,
!    :       'Obs. precip (cm):',obsprec(iday,istn),
!    :       'Anal. prec (cm):',prec(iptstn(istn),jptstn(istn))
        totpreca(istn)=totpreca(istn)+                                  &
                    prec(iptstn(istn),jptstn(istn))
      END DO
      1020   FORMAT(a6,i4,3X,a5,i4,a,i2.2,a,i2.2,3X,a17,f5.2,3X,a16,f5.2)

    END DO
!
    WRITE(6,*)
    WRITE(6,*) 'Top layer API analysis:'
    WRITE(6,*) '-----------------------'
    WRITE(6,*)
!
    WRITE(6,1040) j5,api1(i1,j5,1),api1(i2,j5,1),                       &
                     api1(i3,j5,1),api1(i4,j5,1),                       &
                     api1(i5,j5,1)

    WRITE(6,1040) j4,api1(i1,j4,1),api1(i2,j4,1),                       &
                     api1(i3,j4,1),api1(i4,j4,1),                       &
                     api1(i5,j4,1)

    WRITE(6,1040) j3,api1(i1,j3,1),api1(i2,j3,1),                       &
                     api1(i3,j3,1),api1(i4,j3,1),                       &
                     api1(i5,j3,1)

    WRITE(6,1040) j2,api1(i1,j2,1),api1(i2,j2,1),                       &
                     api1(i3,j2,1),api1(i4,j2,1),                       &
                     api1(i5,j2,1)

    WRITE(6,1040) j1,api1(i1,j1,1),api1(i2,j1,1),                       &
                     api1(i3,j1,1),api1(i4,j1,1),                       &
                     api1(i5,j1,1)

    WRITE(6,*)
    WRITE(6,1035) i1,i2,i3,i4,i5
    WRITE(6,*)
!
    WRITE(6,*)
    WRITE(6,*) 'Root zone API analysis:'
    WRITE(6,*) '-----------------------'
    WRITE(6,*)
!
    WRITE(6,1040) j5,api2(i1,j5,1),api2(i2,j5,1),                       &
                     api2(i3,j5,1),api2(i4,j5,1),                       &
                     api2(i5,j5,1)

    WRITE(6,1040) j4,api2(i1,j4,1),api2(i2,j4,1),                       &
                     api2(i3,j4,1),api2(i4,j4,1),                       &
                     api2(i5,j4,1)

    WRITE(6,1040) j3,api2(i1,j3,1),api2(i2,j3,1),                       &
                     api2(i3,j3,1),api2(i4,j3,1),                       &
                     api2(i5,j3,1)

    WRITE(6,1040) j2,api2(i1,j2,1),api2(i2,j2,1),                       &
                     api2(i3,j2,1),api2(i4,j2,1),                       &
                     api2(i5,j2,1)

    WRITE(6,1040) j1,api2(i1,j1,1),api2(i2,j1,1),                       &
                     api2(i3,j1,1),api2(i4,j1,1),                       &
                     api2(i5,j1,1)

    WRITE(6,*)
    WRITE(6,1035) i1,i2,i3,i4,i5
    WRITE(6,*)
!
    WRITE(6,*)
    WRITE(6,*) 'Dependence of root zone API on initial API: '
    WRITE(6,*) '(most useful when adj2 .ne. 1.0 -- see code)'
    WRITE(6,*) '--------------------------------------------'
    WRITE(6,*)
!
    WRITE(6,1040) j5,kk2dep(i1,j5),kk2dep(i2,j5),                       &
                     kk2dep(i3,j5),kk2dep(i4,j5),                       &
                     kk2dep(i5,j5)

    WRITE(6,1040) j4,kk2dep(i1,j4),kk2dep(i2,j4),                       &
                     kk2dep(i3,j4),kk2dep(i4,j4),                       &
                     kk2dep(i5,j4)

    WRITE(6,1040) j3,kk2dep(i1,j3),kk2dep(i2,j3),                       &
                     kk2dep(i3,j3),kk2dep(i4,j3),                       &
                     kk2dep(i5,j3)

    WRITE(6,1040) j2,kk2dep(i1,j2),kk2dep(i2,j2),                       &
                     kk2dep(i3,j2),kk2dep(i4,j2),                       &
                     kk2dep(i5,j2)

    WRITE(6,1040) j1,kk2dep(i1,j1),kk2dep(i2,j1),                       &
                     kk2dep(i3,j1),kk2dep(i4,j1),                       &
                     kk2dep(i5,j1)

    WRITE(6,*)
    WRITE(6,1035) i1,i2,i3,i4,i5
    WRITE(6,*)
!
!--------------------------------------------------------------------
!
!  Convert the API's to wetsfc and wetdp.  Limit both 'wetsfc' and
!  'wetdp' to ~'wsat' as a maximum.  Limit 'wetsfc' to 'wgeq' as
!  a minimum and limit 'wetdp' to ~'wwlt' as a minimum.  If k1=0.0
!  set 'wetsfc' equal to 'wgeq' and if k1<0.0 set 'wetsfc' ~ equal
!  to 0.8*'wgeq'.  See input file for more details.
!
!  Also, convert API from scalar points back to x,y grid.
!
!--------------------------------------------------------------------
!
    IF(soilinit2 == 3) THEN
      wsfcscl=wgrat
      wdpscl=w2rat
    ELSE
      wsfcscl=1.0
      wdpscl=1.0
    END IF
!
    DO is=1,nstyp
      DO iy=1,ny-1
        DO ix=1,nx-1
          IF( soiltyp(ix,iy,is) /= 0 ) THEN

            qsoil(ix,iy,1,is) =wsfcscl*api1(ix,iy,is)/(depth1*100.0)
            DO k = 2,nzsoil
              qsoil(ix,iy,k,is) =wdpscl*api2(ix,iy,is)/(d2*100.0)
            END DO

            wetcanp(ix,iy,is)=wetcanp0
!
            IF( soiltyp(ix,iy,is) /= 13 ) THEN
              wmin=0.95*wwlt(soiltyp(ix,iy,is)) ! arbitrary min,max values
              wmax=0.95*wsat(soiltyp(ix,iy,is)) ! for qsoil(nzsoil)...
            ELSE  ! if soiltyp is water
              wmin=wwlt(soiltyp(ix,iy,is))
              wmax=wsat(soiltyp(ix,iy,is))
            END IF
            DO k = 2,nzsoil
              IF( qsoil(ix,iy,k,is) < wmin ) THEN
                qsoil(ix,iy,k,is)=wmin ! keep qsoil(k) close to wwlt
              ELSE IF( qsoil(ix,iy,k,is) > wmax ) THEN
                qsoil(ix,iy,k,is)=wmax ! keep qsoil(k) slightly less than saturated
              END IF
            END DO
!            wmin=qsoil(ix,iy,nzsoil,is)/wsat(soiltyp(ix,iy,is))
            wmin=qsoil(ix,iy,2,is)/wsat(soiltyp(ix,iy,is))
            wmin=wmin - ( awgeq(soiltyp(ix,iy,is)) *                  &
              wmin**pwgeq(soiltyp(ix,iy,is)) )  *                     &
              ( 1.0 - wmin**(8.0*pwgeq(soiltyp(ix,iy,is))) )
            wmin=wmin*wsat(soiltyp(ix,iy,is))
            IF( qsoil(ix,iy,1,is) > wmax ) THEN
              qsoil(ix,iy,1,is)=wsat(soiltyp(ix,iy,is))
            ELSE IF( qsoil(ix,iy,1,is) < wmin ) THEN
              qsoil(ix,iy,1,is)=wmin
            END IF
!
            IF( k1 == 0. ) THEN
              qsoil(ix,iy,1,is)=wmin
            ELSE IF(  k1 < 0. .AND. soiltyp(ix,iy,is) /= 13) THEN
              qsoil(ix,iy,1,is)=wmin-AMIN1(0.2*wmin,0.05)
            END IF

          END IF

        END DO
      END DO
    END DO

    DO j=1,ny-1
      DO i=1,nx-1
        DO k=1,nzsoil 
          tsoil (i,j,k,0) =  0.0
          qsoil (i,j,k,0) =  0.0
        END DO 
        wetcanp(i,j,0) =  0.0
      END DO
    END DO

    DO is = 1,nstyp
      DO j=1,ny-1
        DO i=1,nx-1
          DO k=1,nzsoil 
            tsoil (i,j,k,0) = tsoil  (i,j,k,0)                           &
                         + tsoil (i,j,k,is) * stypfrct(i,j,is)
            qsoil (i,j,k,0) = qsoil (i,j,k,0)                            &
                         + qsoil (i,j,k,is) * stypfrct(i,j,is)

          END DO
 
          wetcanp(i,j,0) = wetcanp(i,j,0)                               &
                         + wetcanp(i,j,is) * stypfrct(i,j,is)
        END DO
      END DO
    END DO

    DO is=0,nstyp
      CALL edgfill(tsoil (1,1,1,is),1,nx,1,nx-1,1,ny,1,ny-1,1,nzsoil, &
               1,nzsoil)
      CALL edgfill(qsoil (1,1,1,is),1,nx,1,nx-1,1,ny,1,ny-1,nzsoil,1, &
               1,nzsoil)
      CALL edgfill(wetcanp(1,1,is),1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
    END DO
!
!-------------------------------------------------------------------
!
!    Print out the calculated values of wetsfc, wetdp, and wetcanp.
!
!-------------------------------------------------------------------
!
    i5=nx
    j5=ny
!
    WRITE(6,*)
    WRITE(6,*) 'wetsfc analysis:'
    WRITE(6,*) '-----------------------'
    WRITE(6,*)
!
    WRITE(6,1040) j5,qsoil(i1,j5,1,0),qsoil(i2,j5,1,0),                 &
                 qsoil(i3,j5,1,0),qsoil(i4,j5,1,0),                     &
                 qsoil(i5,j5,1,0)

    WRITE(6,1040) j4,qsoil(i1,j4,1,0),qsoil(i2,j4,1,0),                 &
                 qsoil(i3,j4,1,0),qsoil(i4,j4,1,0),                     &
                 qsoil(i5,j4,1,0)

    WRITE(6,1040) j3,qsoil(i1,j3,1,0),qsoil(i2,j3,1,0),                 &
                 qsoil(i3,j3,1,0),qsoil(i4,j3,1,0),                     &
                 qsoil(i5,j3,1,0)

    WRITE(6,1040) j2,qsoil(i1,j2,1,0),qsoil(i2,j2,1,0),                 &
                 qsoil(i3,j2,1,0),qsoil(i4,j2,1,0),                     &
                 qsoil(i5,j2,1,0)

    WRITE(6,1040) j1,qsoil(i1,j1,1,0),qsoil(i2,j1,1,0),                 &
                 qsoil(i3,j1,1,0),qsoil(i4,j1,1,0),                     &
                 qsoil(i5,j1,1,0)

    WRITE(6,*)
    WRITE(6,1035) i1,i2,i3,i4,i5
    WRITE(6,*)
!
    WRITE(6,*)
    WRITE(6,*) 'wetdp analysis:'
    WRITE(6,*) '-----------------------'
    WRITE(6,*)
!
    WRITE(6,1040) j5,qsoil(i1,j5,nzsoil,0),qsoil(i2,j5,nzsoil,0),       &
                 qsoil(i3,j5,nzsoil,0),qsoil(i4,j5,nzsoil,0),           &
                 qsoil(i5,j5,nzsoil,0)

    WRITE(6,1040) j4,qsoil(i1,j4,nzsoil,0),qsoil(i2,j4,nzsoil,0),       &
                 qsoil(i3,j4,nzsoil,0),qsoil(i4,j4,nzsoil,0),           &
                 qsoil(i5,j4,nzsoil,0)

    WRITE(6,1040) j3,qsoil(i1,j3,nzsoil,0),qsoil(i2,j3,nzsoil,0),       &
                 qsoil(i3,j3,nzsoil,0),qsoil(i4,j3,nzsoil,0),           &
                 qsoil(i5,j3,nzsoil,0)

    WRITE(6,1040) j2,qsoil(i1,j2,nzsoil,0),qsoil(i2,j2,nzsoil,0),       &
                 qsoil(i3,j2,nzsoil,0),qsoil(i4,j2,nzsoil,0),           &
                 qsoil(i5,j2,nzsoil,0)

    WRITE(6,1040) j1,qsoil(i1,j1,nzsoil,0),qsoil(i2,j1,nzsoil,0),       &
                 qsoil(i3,j1,nzsoil,0),qsoil(i4,j1,nzsoil,0),           &
                 qsoil(i5,j1,nzsoil,0)

    WRITE(6,*)
    WRITE(6,1035) i1,i2,i3,i4,i5
    WRITE(6,*)
!
!  DO 720 istn=1,tstns
!    totpreco=0.0
!    DO 710 iday=1,ndays
!      IF(obsprec(iday,istn).gt.0.)
!    :    totpreco=totpreco+obsprec(iday,istn)
! 710   CONTINUE
!    write(6,*)
!    write(6,1021) 'Total recorded precip:',totpreco
!    write(6,1021) 'Total analyzed precip:',totpreca(istn)
!1021   format(a22,f6.2)
!    write(6,1022) 'API:',api2(iptstn(istn),jptstn(istn),1)
!1022   format(a4,18x,f6.2)
!    write(6,1023) 'Qsoilsfc:',qsoil(iptstn(istn),jptstn(istn),1,0)
!1023   format(a9,15x,f6.4)
!    write(6,1024) 'Qsoildp:',qsoil(iptstn(istn),jptstn(istn),nzsoil,0)
!1024   format(a8,16x,f6.4)
!    write(6,*)
!    write(6,*)
! 720 CONTINUE
!
    IF( totapi0 > 0 ) THEN
      WRITE(6,*)
      WRITE(6,*) 'Root zone volumetric water contents'                  &
          //' at initial API stations:'
      WRITE(6,*)
      DO i=1,totapi0
        WRITE(6,1042) 'Stn #:',i,' has wetdp of:',                      &
            qsoil(iptapi0(i),jptapi0(i),nzsoil,0)
        1042     FORMAT(a6,1X,i2,a15,1X,f5.3)
      END DO
      WRITE(6,*)
    END IF
!
    WRITE(6,*)
    WRITE(6,*) 'wetcanp analysis:'
    WRITE(6,*) '-----------------------'
    WRITE(6,*)
!
    WRITE(6,1040) j5,wetcanp(i1,j5,0),wetcanp(i2,j5,0),                 &
                 wetcanp(i3,j5,0),wetcanp(i4,j5,0),                     &
                 wetcanp(i5,j5,0)

    WRITE(6,1040) j4,wetcanp(i1,j4,0),wetcanp(i2,j4,0),                 &
                 wetcanp(i3,j4,0),wetcanp(i4,j4,0),                     &
                 wetcanp(i5,j4,0)

    WRITE(6,1040) j3,wetcanp(i1,j3,0),wetcanp(i2,j3,0),                 &
                 wetcanp(i3,j3,0),wetcanp(i4,j3,0),                     &
                 wetcanp(i5,j3,0)

    WRITE(6,1040) j2,wetcanp(i1,j2,0),wetcanp(i2,j2,0),                 &
                 wetcanp(i3,j2,0),wetcanp(i4,j2,0),                     &
                 wetcanp(i5,j2,0)

    WRITE(6,1040) j1,wetcanp(i1,j1,0),wetcanp(i2,j1,0),                 &
                 wetcanp(i3,j1,0),wetcanp(i4,j1,0),                     &
                 wetcanp(i5,j1,0)

    WRITE(6,*)
    WRITE(6,1035) i1,i2,i3,i4,i5
    WRITE(6,*)
!
    1041 FORMAT(i3,3X,f5.1,3X,f5.1,3X,f5.1,3X,f5.1,3X,f5.1)
!
    WRITE(6,*)
    WRITE(6,*) 'tsoil analysis:'
    WRITE(6,*) '-----------------------'
    WRITE(6,*)
!
    DO k=1,nzsoil

    WRITE(6,1041) j5,tsoil(i1,j5,k,0),tsoil(i2,j5,k,0),                   &
                 tsoil(i3,j5,k,0),tsoil(i4,j5,k,0),                       &
                 tsoil(i5,j5,k,0)

    WRITE(6,1041) j4,tsoil(i1,j4,k,0),tsoil(i2,j4,k,0),                   &
                 tsoil(i3,j4,k,0),tsoil(i4,j4,k,0),                       &
                 tsoil(i5,j4,k,0)

    WRITE(6,1041) j3,tsoil(i1,j3,k,0),tsoil(i2,j3,k,0),                   &
                 tsoil(i3,j3,k,0),tsoil(i4,j3,k,0),                       &
                 tsoil(i5,j3,k,0)

    WRITE(6,1041) j2,tsoil(i1,j2,k,0),tsoil(i2,j2,k,0),                   &
                 tsoil(i3,j2,k,0),tsoil(i4,j2,k,0),                       &
                 tsoil(i5,j2,k,0)

    WRITE(6,1041) j1,tsoil(i1,j1,k,0),tsoil(i2,j1,k,0),                   &
                 tsoil(i3,j1,k,0),tsoil(i4,j1,k,0),                       &
                 tsoil(i5,j1,k,0)
    END DO 

    WRITE(6,*)
    WRITE(6,1035) i1,i2,i3,i4,i5
    WRITE(6,*)
!
!--------------------------------------------------------------------
!
!    Read in soil moisture data from a file.
!
!--------------------------------------------------------------------
!
  ELSE IF( soilinit2 == 4 ) THEN

    WRITE(6,'(a,a)') 'Reading soil moisture from ',soilinfl
    CALL readsoil(nx,ny,nzsoil,nstyps,soilinfl,dx,dy,zpsoil,            &
             mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,       &
             zpsoilin,tsoilin,qsoilin,wcanpin,snowdin,          &
                  tem1,qsoil,wetcanp,snowdpth,soiltyp)

    DO is = 1,nstyp
      DO j=1,ny-1
        DO i=1,nx-1
          wetcanp(i,j,is) = wetcanp(i,j,0)

          DO k=1,nzsoil

          qsoil (i,j,k,is) = qsoil (i,j,k,0)
          tsoil (i,j,k,0) = tsoil  (i,j,k,0)                          &
                       + tsoil  (i,j,k,is) * stypfrct(i,j,is)
          END DO 

        END DO
      END DO
    END DO

!
!--------------------------------------------------------------------
!
!    Merge in soil data from soilinfl into soil data provided by history
!    dump where ever there is water, also use snowdpth in soilinfl
!    if present.
!
!--------------------------------------------------------------------
!
  ELSE IF( soilinit2 == 5 ) THEN

    WRITE(6,'(a,a)') 'Merging water soil data from ',soilinfl

    ! Note that nstyp=1 is assumed (or at least only the composite
    ! type is used) for data in soilinfl
    soiltyp2 = soiltyp
    CALL readsoil(nx,ny,nstyp,soilinfl,dx,dy,                          &
             mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,       &
                  zpsoilin,tsoilin,qsoilin,wcanpin,snowdin,          &
                  tsoil2,qsoil2,wetcanp2,snowdpth,soiltyp2)

!wdt Copyright (c) 2001 Weather Decision Technologies, Inc.
!   DO is=1,nstyp
!    DO k=1,nzsoil 
!     DO j=1,ny-1
!       DO i=1,nx-1
!         IF (soiltyp(i,j,is) == 13) THEN
!           IF (tsoilin /= 0)  tsoil(i,j,k,is)  = tsoil2(i,j,k,0)
!           IF (qsoilin /= 0)  qsoil(i,j,k,is)  = qsoil2(i,j,k,0)
!           IF (wcanpin /= 0)  wetcanp(i,j,is) = wetcanp2(i,j,0)
!         END IF
!       END DO
!      END DO 
!     END DO
!   END DO

    DO is=1,nstyp
      DO j=1,ny-1
        DO i=1,nx-1
          IF (soiltyp(i,j,is) == 13) THEN
            DO is2=1,nstyp
              IF (soiltyp2(i,j,is2) == 13) THEN
                DO k=1,nzsoil
                  IF (tsoilin /= 0) tsoil(i,j,k,is) = tsoil2(i,j,k,is2)
                  IF (qsoilin /= 0) qsoil(i,j,k,is) = qsoil2(i,j,k,is2)
                END DO 
                IF (wcanpin /= 0)  wetcanp(i,j,is) = wetcanp2(i,j,is2)
                EXIT
              END IF
            END DO
          END IF
        END DO
      END DO
    END DO

        ! is=0 filled in when program returns to main program

!  ELSE IF( soilinit2 == 6 ) THEN
!
!    DO is=1,nstyps
!    DO j=1,ny-1
!    DO i=1,nx-1
!
!      IF( soiltyp(i,j,is) /= 0 ) THEN
!        IF ( soiltyp(i,j,is).eq.13 ) THEN
!          wetsfc (i,j,is) = 1.0
!          wetdp  (i,j,is) = 1.0
!          wetcanp(i,j,is) = 0.0
!        ELSE
!         wetsfc (i,j,is) = soim1_interp(i,j) + 0.05  ! RUC offset
!!          wetsfc (i,j,is) = soim1_interp(i,j)
!          wmin=0.95*(wwlt(soiltyp(i,j,is)))
!          wmax=0.95*(wsat(soiltyp(i,j,is)))
!                                                                                                                                   
!          IF (wetsfc(i,j,is).gt.wmax) THEN
!      print *,'Resetting wetsfc to 0.95*wsat at I/J/IS= ',i,j,is
!      print *,'WETSFC/0.95*WSAT = ',wetsfc(i,j,is),wmax
!            wetsfc(i,j,is)= wmax
!          ENDIF
!          IF (wetsfc(i,j,is).lt.(wmin-amin1(0.2*wmin,0.05))) THEN
!      print *,'Resetting wetsfc to wmin at I/J/IS= ',i,j,is
!      print *,'WETSFC/WMIN = ',wetsfc(i,j,is),(wmin-amin1(0.2*wmin,0.05))
!            wetsfc(i,j,is)= wmin-amin1(0.2*wmin,0.05)
!          ENDIF
!                                                                                                                                   
!          wetdp  (i,j,is) = soim100_interp(i,j) + 0.05  ! RUC offset
!!          wetdp  (i,j,is) = soim100_interp(i,j)
!                                                                                                                                   
!          IF (wetdp(i,j,is).gt.wmax) THEN
!      print *,'Resetting wetdp to 0.95*wsat at I/J/IS= ',i,j,is
!      print *,'WETDP/0.95*WSAT = ',wetdp(i,j,is),wmax
!            wetdp(i,j,is)= wmax
!          ENDIF
!          IF (wetdp(i,j,is).lt.wmin) THEN
!      print *,'Resetting wetdp to wmin at I/J/IS= ',i,j,is
!      print *,'WETDP/WMIN = ',wetdp(i,j,is),wmin
!            wetdp(i,j,is)= wmin
!          ENDIF
!          wetcanp(i,j,is) = wetcanp0
!        ENDIF
!      ENDIF
!
!    ENDDO
!    ENDDO
!    ENDDO
  
  ELSE

    WRITE(6,'(i6,a)') soilinit2,' is not a valid soilinit2 option!'
    STOP

  END IF

  !wdt update
  IF (soilinit2 == 5) THEN
    DEALLOCATE(tsoil2)
    DEALLOCATE(qsoil2)
    DEALLOCATE(wetcanp2)
    DEALLOCATE(soiltyp2)
  END IF
!
!-------------------------------------------------------------------
!
!    End of subroutine
!
!-------------------------------------------------------------------
!
  PRINT *,'Finishes mksoilvar OK!!'
!
  RETURN
!
END SUBROUTINE mksoilvar
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE FINDLC                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE findlc(nx,ny,xs,ys,xpt,ypt,ipt,jpt,ireturn)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Searches in x and y to find i,j location of xpt, ypt.
!
!  X and Y do not have to be on a regular grid, however it is
!  assumed that x and y are monotonically increasing as i and j
!  indices increase.
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  April 1992.
!
!  MODIFICATION HISTORY:
!
!  February, 1993 (K. Brewster)
!  Additional documentation for ARPS 3.1 release
!
!  October, 1994 (K. Brewster)
!  Changed to reference scalar points.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    xs       x coordinate of scalar points in physical/comp. space (m)
!    ys       y coordinate of scalar points in physical/comp. space (m)
!
!    xpt      location to find in x coordinate (m)
!    ypt      location to find in y coordinate (m)
!
!  OUTPUT:
!
!    ipt      i index to the west of desired location
!    jpt      j index to the south of desired location
!
!-----------------------------------------------------------------------
!
!  Arguments
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny          ! Dimensions of ARPS grids
  REAL :: xs(nx)            ! x coordinate of scalar grid points in
                            ! physical/comp. space (m)
  REAL :: ys(ny)            ! y coordinate of grid points in
                            ! physical/comp. space (m)

  REAL :: xpt               ! location to find in x coordinate
  REAL :: ypt               ! location to find in y coordinate
  INTEGER :: ipt            ! i index to the west of desired
                            ! location
  INTEGER :: jpt            ! j index to the south of desired
                            ! location
  INTEGER :: ireturn
!
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ireturn=0
!
  DO i=2,nx
    IF(xpt < xs(i)) EXIT
  END DO
  101 CONTINUE
  ipt=i-1
  IF( xpt > xs(nx-1) .OR. xpt < xs(1) ) ireturn=-1
  DO j=2,ny
    IF( ypt < ys(j)) EXIT
  END DO
  201 CONTINUE
  jpt=j-1
  IF( ypt > ys(ny-1) .OR. ypt < ys(1) ) ireturn=-2

  RETURN
END SUBROUTINE findlc
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE RDNCEPPR                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE rdnceppr(nx,ny,mxstns,mxdays,xs,ys,xll,yll,                  &
           iendyr,iendmo,ienddy,ndays,                                  &
           prcpdir,prcplst,                                             &
           staid,iptstn,jptstn,obsprec,totpreca,                        &
           tstns)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Program to use history data dump and surface property data to
!  create the soil variables
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!    04/11/1995
!
!  MODIFICATION HISTORY:
!
!  05/07/1998 (Yuhe Liu)
!  Hub-CAPS's modification to the precipitation format change in
!  extracting the information from each line.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Arguments
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,mxstns,mxdays
!
  REAL :: xs(nx)
  REAL :: ys(ny)
  REAL :: xll
  REAL :: yll
  INTEGER :: iendyr,iendmo,ienddy,ndays
  CHARACTER (LEN=*) :: prcpdir
  CHARACTER (LEN=*) :: prcplst
  CHARACTER (LEN=7) :: staid(mxstns)
  INTEGER :: iptstn(mxstns)
  INTEGER :: jptstn(mxstns)
  REAL :: obsprec(mxdays,mxstns)
  REAL :: totpreca(mxstns)
  INTEGER :: tstns
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: mxstafile
  PARAMETER( mxstafile = 99000)
  INTEGER :: sec24hr
  PARAMETER( sec24hr = 86400 )

  REAL :: in2cm,maxprcm
! J.Case, ENSCO Inc. (10/26/2004)
!  PARAMETER (in2cm=2.54,maxprcm=4.00)
  PARAMETER (in2cm=2.54,maxprcm=8.00)
!
  CHARACTER (LEN=7) :: rdstaid
  CHARACTER (LEN=7) :: rdst(3)
  CHARACTER (LEN=7) :: dummy
  CHARACTER (LEN=256) :: rfcname
  CHARACTER (LEN=80)  :: line
  INTEGER :: iunit,istn,jstn,kstn,iday,out
  INTEGER :: ilat(3),ilon(3)
  INTEGER :: ipt,jpt,ios,ireturn,lendir
  INTEGER :: istart,itime,iyr,imon,idy,ihr,imin,isec
  INTEGER :: nread,nmatch
  REAL :: lat,lon,xpt,ypt
  REAL :: rdlat,rdlon,rdobpr
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-------------------------------------------------------------------
!
!  Initialize precipitation observations to missing
!
!-------------------------------------------------------------------
!
  DO istn=1,mxstns
    staid(istn)='     '
    DO iday=1,mxdays
      obsprec(iday,istn)=-9.0
    END DO
  END DO
!
!-------------------------------------------------------------------
!
!  Get list of precip stations
!  and keep those that are within the domain.
!
!-------------------------------------------------------------------
!
  iunit=38
  OPEN(iunit,ERR=910,FILE=prcplst,IOSTAT=ios,                           &
             STATUS='old',FORM='formatted')
! J.Case, ENSCO Inc. (10/26/2004)
!  READ(iunit,'(a5)') dummy
  out=0
  nread=0

! J.Case -- commented out read statement -- rewrote.
! J.Case, ENSCO Inc. (10/26/2004)

!  DO istn=1,mxstns,3
!    READ(iunit,'(a7,i5,i6,8x,a7,i5,i6,8x,a7,i5,i6)',                    &
!              ERR=101,END=101)                                          &
!              rdst(1),ilat(1),ilon(1),                                  &
!              rdst(2),ilat(2),ilon(2),                                  &
!              rdst(3),ilat(3),ilon(3)
!    nread=nread+3
!    DO kstn=1,3
!      jstn=istn-out+(kstn-1)
!      lat = 0.01*FLOAT(ilat(kstn))
!      lon =-0.01*FLOAT(ilon(kstn))
!

  DO istn=1,mxstns
    READ(iunit,'(1x,f5.2,3x,f5.2,8x,a7)',ERR=101,END=101)               &
              lat,lon,rdst(1)                                           
    nread=nread+1
    DO kstn=1,1
      jstn=istn-out
      lon = -lon
      CALL lltoxy(1,1,lat,lon,xpt,ypt)
      xpt=xpt-xll
      ypt=ypt-yll
      CALL findlc(nx,ny,xs,ys,xpt,ypt,ipt,jpt,ireturn)
      IF( ireturn == 0 ) THEN

        staid(jstn)=rdst(kstn)
        iptstn(jstn)=ipt
        jptstn(jstn)=jpt
        totpreca(jstn)=0.0
        WRITE(6,*) ' Precip stn ',rdst(kstn),' was found near i,j:',    &
                   ipt,jpt

      ELSE

        write(6,*) ' Precip stn ',rdst(kstn),' is outside of the grid.'
        out=out+1

      END IF
    END DO
  END DO

  101 CONTINUE
  tstns=istn-out-1
  CLOSE(iunit)
  WRITE(6,'(a,i6,a)') ' Read ',nread,' precip station locations.'
  WRITE(6,'(i12,a)') tstns,' inside grid.'
!
!-------------------------------------------------------------------
!
!  Begin main day loop, reading precip file for each day.
!
!-------------------------------------------------------------------
!
  lendir = LEN( prcpdir )
  CALL strlnth( prcpdir, lendir )

  CALL ctim2abss(iendyr,iendmo,ienddy,12,0,0, istart)
  istart=istart-((ndays-1)*sec24hr)
!
  DO iday=1,ndays
    itime = istart + ((iday-1)*sec24hr)
    CALL abss2ctim( itime, iyr, imon, idy, ihr, imin, isec )
    WRITE(rfcname,'(a,a,i4.4,i2.2,i2.2)')                               &
        prcpdir(1:lendir),'/usa-dlyprcp-',iyr,imon,idy
    WRITE(6,'(/a,a)') ' Opening: ',rfcname
! J.Case, ENSCO Inc. (10/26/2004) -- Fixed bug that exited out of loop entirely if 
! one day of data is missing.
    OPEN(iunit,ERR=351,FILE=rfcname,                                    &
         STATUS='old',FORM='formatted')
    nread=0
    nmatch=0
    READ(iunit,'(a5)',ERR=351,END=351) dummy
    DO kstn=1,mxstafile
      READ(iunit,'(a80)',ERR=350,END=351) line
      READ(line,*,ERR=350,END=350) rdlat,rdlon,rdobpr
      IF (line(23:23) == ' ') THEN
        rdstaid = line(24:31)
      ELSE
        rdstaid = line(23:30)
      END IF

      nread=nread+1
      DO istn=1,tstns
        IF( rdstaid == staid(istn) ) THEN
          nmatch=nmatch+1
          obsprec(iday,istn)=in2cm*rdobpr
          obsprec(iday,istn)=AMIN1(maxprcm,obsprec(iday,istn))
          EXIT
        END IF
      END DO
      301     CONTINUE
    END DO
    350 CONTINUE
    351 CONTINUE
    CLOSE(iunit)
    WRITE(6,'(a,i6,a)') ' Read ',nread,' precip obs.'
    WRITE(6,'(i12,a)') nmatch,' obs matched inside grid.'
  END DO
!
!-------------------------------------------------------------------
!
!    Normal end
!
!-------------------------------------------------------------------
!
  RETURN
!
!-------------------------------------------------------------------
!
!    Error destinations
!
!-------------------------------------------------------------------
!
910 CONTINUE
  WRITE(6,'(a,a,/a,i6)')                                                &
  'Error opening list of precip stations,',prcplst,' iostat = ',ios
  STOP

END SUBROUTINE rdnceppr

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
SUBROUTINE rdextsst(nx,ny,infilename,                                   &
                     grid,blats,blons,res,nisst,njsst,sst_flag)
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  INTEGER :: nx,ny
  REAL :: grid(nx,ny), blats, blons, res

  INTEGER :: nisst, njsst
  LOGICAL :: sst_flag

  CHARACTER (LEN=80) dummy
  CHARACTER (LEN=10) infilename

  blats = 0.
  blons = -180.
  res   = 1.

  OPEN (unit=50,file=infilename,form='FORMATTED',                       &
        status='OLD',iostat=ios,err=100)

  100 sst_flag = .true.
!  print *,ios,' ',infilename,nx,ny
  print *,'INFILE = ',infilename
  IF (ios .ne. 0) THEN
    sst_flag = .false.
    RETURN
  END IF
!
  read (50,'(//////,71x,i4,1x,i4)') nisst,njsst
  print *,'Total grid dimensions = ',nisst,njsst
  read (50,'(10x,2I5,11x,2I5)') ix1, ix2, iy1, iy2
  print *,'Subset domain = ',ix1,ix2,iy1,iy2
  read (50,'(///)')
  read (50,110) (kk,k = ix1, ix2 )
  110 format (9x,8(3X,I4,2X),/(8x,8(3X,I4,2X)))
!
!     Loop through rows from the top down.
!
  DO j = iy2, iy1, -1
    read (50,120) jj,(grid(k,j),k=ix1,ix2)
  END DO
  120 format (4x,I3,1X,8F9.4,/(8X,8F9.4))

  RETURN
END SUBROUTINE rdextsst

