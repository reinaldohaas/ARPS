PROGRAM arpsprt
!
!##################################################################
!##################################################################
!######                                                      ######
!######                   PROGRAM ARPSPRT                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Program to examine ARPS grids that have been saved as
!  history files.  This prints out the variables to a file or
!  screen.
!
!  It shares with the model the include files  'globcst.inc'
!  for storage parameters.
!
!  It reads in a history file produced by ARPS 3.0 in a user
!  specified format.
!
!  Parameters grdin,basin,mstin,icein,trbin are read in from the
!  data file itself, therefore are determined internally.
!  Arrays that are not read in retain their initial zero values.
!  These parameters are passed among subroutines through
!  a common block defined in 'indtflg.inc'.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster  OU School of Meteorology. April 1992
!
!  MODIFICATION HISTORY:
!   14 May 1992  (KB) changed from arps2.5 to arps3.0
!   03 Aug 1992  (KB) updated to account for changes in arps3.0
!
!   8/27/1992 (M. Xue)
!   To call dtaread to read new data format.
!
!   8/30/1992 (K. Brewster)
!   Moved label 101 so that grid file name is not reentered for
!   second and subsequent data reads.
!
!   9/1/94 (Y. Lu)
!   Cleaned up documentation.
!
!   10/11/1994 (K. Brewster)
!   Further update for dtaread.
!
!-----------------------------------------------------------------------
!
!  DATA ARRAYS READ IN:
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (km)
!    zp       z coordinate of grid points in computational space (m)
!    zpsoil   z coordinate of grid points of the soil model
!
!    uprt     x component of perturbation velocity (m/s)
!    vprt     y component of perturbation velocity (m/s)
!    wprt     Vertical component of perturbation velocity in Cartesian
!             coordinates (m/s).
!
!    ptprt    Perturbation potential temperature (K)
!    pprt     Perturbation pressure (Pascal)
!
!    qvprt    Perturbation water vapor mixing ratio (kg/kg)
!    qscalar  Hydrometeor scalars
!
!    ubar     Base state x velocity component (m/s)
!    vbar     Base state y velocity component (m/s)
!    wbar     Base state z velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    rhobar   Base state air density (kg/m**3)
!    qvbar    Base state water vapor mixing ratio (kg/kg)
!
!  CALCULATED DATA ARRAYS:
!
!    u        x component of velocity (m/s)
!    v        y component of velocity (m/s)
!    w        z component of velocity (m/s)
!    pt       Potential temperature (K)
!    qv       Water vapor mixing ratio (kg/kg)
!
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  INTEGER :: nzsoil

!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'indtflg.inc'
!
!-----------------------------------------------------------------------
!
!  Arrays to be read in:
!
!-----------------------------------------------------------------------
!

  REAL, ALLOCATABLE :: x     (:)         ! The x-coord. of the physical and
                                         ! computational grid. Defined at u-point.
  REAL, ALLOCATABLE :: y     (:)         ! The y-coord. of the physical and
                                         ! computational grid. Defined at v-point.
  REAL, ALLOCATABLE :: z     (:)         ! The z-coord. of the computational grid.
                                         ! Defined at w-point on the staggered grid.
  REAL, ALLOCATABLE :: zp    (:,:,:)     ! The physical height coordinate defined at
                                         ! w-point of the staggered grid.
  REAL, ALLOCATABLE :: zpsoil  (:,:,:)   ! The physical height coordinate
                                         ! of the soil levels
  REAL, ALLOCATABLE :: u     (:,:,:)     ! Total u-velocity (m/s)
  REAL, ALLOCATABLE :: v     (:,:,:)     ! Total v-velocity (m/s)
  REAL, ALLOCATABLE :: w     (:,:,:)     ! Total w-velocity (m/s)
  REAL, ALLOCATABLE :: pt    (:,:,:)     ! Total potential temperature (K)

  REAL, ALLOCATABLE :: qv    (:,:,:)     ! Water vapor specific humidity (kg/kg)
  REAL, ALLOCATABLE :: qscalar (:,:,:,:)
  REAL, ALLOCATABLE :: tke   (:,:,:)     ! Turbulent Kinetic Energy ((m/s)**2)
  REAL, ALLOCATABLE :: kmh   (:,:,:)     ! Horizontal turb. mixing coef. for
                                         ! momentum. ( m**2/s )
  REAL, ALLOCATABLE :: kmv   (:,:,:)     ! Vertical turb. mixing coef. for
                                         ! momentum. ( m**2/s )
  REAL, ALLOCATABLE :: ubar  (:,:,:)     ! Base state u-velocity (m/s)
  REAL, ALLOCATABLE :: vbar  (:,:,:)     ! Base state v-velocity (m/s)
  REAL, ALLOCATABLE :: wbar  (:,:,:)     ! Base state w-velocity (m/s)
  REAL, ALLOCATABLE :: ptbar (:,:,:)     ! Base state potential temperature (K)
  REAL, ALLOCATABLE :: rhobar(:,:,:)     ! Base state air density (kg/m**3)
  REAL, ALLOCATABLE :: pbar  (:,:,:)     ! Base state pressure (Pascal)
  REAL, ALLOCATABLE :: qvbar (:,:,:)     ! Base state water vapor specific humidity
                                         ! (kg/kg)
  INTEGER :: nstyps
!  PARAMETER (nstyps=4)

  INTEGER, ALLOCATABLE :: soiltyp (:,:,:)! Soil type
  REAL, ALLOCATABLE :: stypfrct(:,:,:)   ! Fraction of soil types
  INTEGER , ALLOCATABLE:: vegtyp  (:,:)  ! Vegetation type
  REAL, ALLOCATABLE :: lai     (:,:)     ! Leaf Area Index
  REAL, ALLOCATABLE :: roufns  (:,:)     ! Surface roughness
  REAL, ALLOCATABLE :: veg     (:,:)     ! Vegetation fraction

  REAL, ALLOCATABLE :: tsoil  (:,:,:,:)  ! soil temperature (K)
  REAL, ALLOCATABLE :: qsoil  (:,:,:,:)  ! soil moisture
  REAL, ALLOCATABLE :: wetcanp(:,:,:)    ! Canopy water amount
  REAL , ALLOCATABLE:: snowdpth(:,:)     ! Snow depth (m)

  REAL, ALLOCATABLE :: raing  (:,:)      ! Grid supersaturation rain
  REAL, ALLOCATABLE :: rainc  (:,:)      ! Cumulus convective rain
  REAL, ALLOCATABLE :: prcrate(:,:,:)    ! precipitation rate (kg/(m**2*s))
                                         ! prcrate(1,1,1) = total precip. rate
                                         ! prcrate(1,1,2) = grid scale precip. rate
                                         ! prcrate(1,1,3) = cumulus precip. rate
                                         ! prcrate(1,1,4) = microphysics precip. rate

  REAL, ALLOCATABLE :: radfrc(:,:,:)     ! Radiation forcing (K/s)
  REAL, ALLOCATABLE :: radsw (:,:)       ! Solar radiation reaching the surface
  REAL, ALLOCATABLE :: rnflx (:,:)       ! Net radiation flux absorbed by surface
  REAL, ALLOCATABLE :: radswnet (:,:) ! Net solar radiation, SWin - SWout
  REAL, ALLOCATABLE :: radlwin  (:,:) ! Incoming longwave radiation


  REAL, ALLOCATABLE :: usflx (:,:)       ! Surface flux of u-momentum (kg/(m*s**2))
  REAL, ALLOCATABLE :: vsflx (:,:)       ! Surface flux of v-momentum (kg/(m*s**2))
  REAL, ALLOCATABLE :: ptsflx(:,:)       ! Surface heat flux (K*kg/(m*s**2))
  REAL, ALLOCATABLE :: qvsflx(:,:)       ! Surface moisture flux (kg/(m**2*s))

  REAL, ALLOCATABLE :: uprt  (:,:,:)     ! Perturbation u-velocity (m/s)
  REAL, ALLOCATABLE :: vprt  (:,:,:)     ! Perturbation v-velocity (m/s)
  REAL, ALLOCATABLE :: wprt  (:,:,:)     ! Perturbation w-velocity (m/s)
  REAL, ALLOCATABLE :: ptprt (:,:,:)     ! Perturbation potential temperature (K)
  REAL, ALLOCATABLE :: pprt  (:,:,:)     ! Perturbation pressure (Pascal)
  REAL, ALLOCATABLE :: qvprt (:,:,:)     ! Perturbation water vapor specific
                                         ! humidity (kg/kg)

!
!-----------------------------------------------------------------------
!
!  Other data variables
!
!-----------------------------------------------------------------------
!
  REAL :: time
!
!-----------------------------------------------------------------------
!
!  Temporary working arrays:
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: tem1(:,:,:)
  REAL, ALLOCATABLE :: tem2(:,:,:)
  REAL, ALLOCATABLE :: tem3(:,:,:)
!
!-----------------------------------------------------------------------
!
!  User request stuff
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=80) :: inline(5)
  INTEGER :: maxpar
  PARAMETER (maxpar=10)
  CHARACTER (LEN=2) :: params(maxpar)
  INTEGER :: islice(5),jslice(5),kslice(5)
  INTEGER :: nparams,ip
  INTEGER :: ireply(5)
  INTEGER :: nreply
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,ireturn,mode,nm
!
  CHARACTER (LEN=256) :: filename
  CHARACTER (LEN=256) :: grdbasfn
  CHARACTER (LEN=25)  :: title
  INTEGER :: lenfil,lengbf
!
  INTEGER :: iplot(3,5)
  INTEGER :: ib(3),ie(3),jb(3),je(3),kb(3),ke(3)
  INTEGER :: nchin

  INTEGER :: istatus

!-----------------------------------------------------------------------
!
!  NAMELIST definition
!
!-----------------------------------------------------------------------

!  NAMELIST /grid_dims/ nx, ny, nz, nzsoil
  NAMELIST /input_fn/ hdmpfmt, grdbasfn, filename
  NAMELIST /field_number/ nreply
  NAMELIST /plot_options/ ireply, inline, iplot, kslice, jslice, islice
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  WRITE(6,'(/11(/5x,a)/)')                                              &
     '###############################################################', &
     '###############################################################', &
     '###                                                         ###', &
     '###                  Welcome to ARPSPRT                     ###', &
     '###      This program reads in the history dump data        ###', &
     '###      sets generated by ARPS, and prints the array       ###', &
     '###      contents as 2-D arrays tables at user specified    ###', &
     '###      slices.                                            ###', &
     '###                                                         ###', &
     '###############################################################', &
     '###############################################################'

!
!-----------------------------------------------------------------------
!
!  Get the name of the input data set.
!
!-----------------------------------------------------------------------
!
  READ(5,input_fn,END=100)
  WRITE(6,'(/a,a)')' Namelist block input_fn successfully read.'

  lengbf=LEN_trim(grdbasfn)
  lenfil=LEN_trim(filename)
  WRITE(6,'(/a,a)')' The data set name is ', filename(1:lenfil)

  CALL get_dims_from_data(hdmpfmt,filename(1:lenfil),nx,ny,nz,nzsoil,   &
                          nstyps,ireturn)

  nstyp = nstyps ! Copy to global variable

!-----------------------------------------------------------------------
!
!  Allocate variables and initialize to zero
!
!-----------------------------------------------------------------------

  ALLOCATE(x(nx),STAT=istatus)
  x=0
  ALLOCATE(y(ny),STAT=istatus)
  y=0
  ALLOCATE(z(nz),STAT=istatus)
  z=0
  ALLOCATE(zp(nx,ny,nz),STAT=istatus)
  zp=0
  ALLOCATE(zpsoil(nx,ny,nzsoil),STAT=istatus)
  zpsoil=0
  ALLOCATE(u(nx,ny,nz),STAT=istatus)
  u=0
  ALLOCATE(v(nx,ny,nz),STAT=istatus)
  v=0
  ALLOCATE(w(nx,ny,nz),STAT=istatus)
  w=0
  ALLOCATE(pt(nx,ny,nz),STAT=istatus)
  pt=0
  ALLOCATE(qv(nx,ny,nz),STAT=istatus)
  qv=0
  ALLOCATE(qscalar(nx,ny,nz,nscalar),STAT=istatus)
  qscalar = 0.0
  ALLOCATE(tke(nx,ny,nz),STAT=istatus)
  tke=0
  ALLOCATE(kmh(nx,ny,nz),STAT=istatus)
  kmh=0
  ALLOCATE(kmv(nx,ny,nz),STAT=istatus)
  kmv=0
  ALLOCATE(ubar(nx,ny,nz),STAT=istatus)
  ubar=0
  ALLOCATE(vbar(nx,ny,nz),STAT=istatus)
  vbar=0
  ALLOCATE(wbar(nx,ny,nz),STAT=istatus)
  wbar=0
  ALLOCATE(ptbar(nx,ny,nz),STAT=istatus)
  ptbar=0
  ALLOCATE(rhobar(nx,ny,nz),STAT=istatus)
  rhobar=0
  ALLOCATE(pbar(nx,ny,nz),STAT=istatus)
  pbar=0
  ALLOCATE(qvbar(nx,ny,nz),STAT=istatus)
  qvbar=0
  ALLOCATE(soiltyp(nx,ny,nstyps),STAT=istatus)
  soiltyp=0
  ALLOCATE(stypfrct(nx,ny,nstyps),STAT=istatus)
  stypfrct=0
  ALLOCATE(vegtyp(nx,ny),STAT=istatus)
  vegtyp=0
  ALLOCATE(lai(nx,ny),STAT=istatus)
  lai=0
  ALLOCATE(roufns(nx,ny),STAT=istatus)
  roufns=0
  ALLOCATE(veg(nx,ny),STAT=istatus)
  veg=0
  ALLOCATE(tsoil(nx,ny,nzsoil,0:nstyps),STAT=istatus)
  tsoil=0
  ALLOCATE(qsoil(nx,ny,nzsoil,0:nstyps),STAT=istatus)
  qsoil=0
  ALLOCATE(wetcanp(nx,ny,0:nstyps),STAT=istatus)
  wetcanp=0
  ALLOCATE(snowdpth(nx,ny),STAT=istatus)
  snowdpth=0
  ALLOCATE(raing(nx,ny),STAT=istatus)
  raing=0
  ALLOCATE(rainc(nx,ny),STAT=istatus)
  rainc=0
  ALLOCATE(prcrate(nx,ny,4),STAT=istatus)
  prcrate=0
  ALLOCATE(radfrc(nx,ny,nz),STAT=istatus)
  radfrc=0
  ALLOCATE(radsw(nx,ny),STAT=istatus)
  radsw=0
  ALLOCATE(rnflx(nx,ny),STAT=istatus)
  rnflx=0
  ALLOCATE(radswnet(nx,ny),STAT=istatus)
  radswnet=0
  ALLOCATE(radlwin(nx,ny),STAT=istatus)
  radlwin=0
  ALLOCATE(usflx(nx,ny),STAT=istatus)
  usflx=0
  ALLOCATE(vsflx(nx,ny),STAT=istatus)
  vsflx=0
  ALLOCATE(ptsflx(nx,ny),STAT=istatus)
  ptsflx=0
  ALLOCATE(qvsflx(nx,ny),STAT=istatus)
  qvsflx=0
  ALLOCATE(uprt(nx,ny,nz),STAT=istatus)
  uprt=0
  ALLOCATE(vprt(nx,ny,nz),STAT=istatus)
  vprt=0
  ALLOCATE(wprt(nx,ny,nz),STAT=istatus)
  wprt=0
  ALLOCATE(ptprt(nx,ny,nz),STAT=istatus)
  ptprt=0
  ALLOCATE(pprt(nx,ny,nz),STAT=istatus)
  pprt=0
  ALLOCATE(qvprt(nx,ny,nz),STAT=istatus)
  qvprt=0
  ALLOCATE(tem1(nx,ny,nz),STAT=istatus)
  tem1=0
  ALLOCATE(tem2(nx,ny,nz),STAT=istatus)
  tem2=0
  ALLOCATE(tem3(nx,ny,nz),STAT=istatus)
  tem3=0
!
!-----------------------------------------------------------------------
!
!  A few initializations
!
!-----------------------------------------------------------------------

  DO mode=1,3
    ib(mode)=2
    ie(mode)=nx-2
    jb(mode)=2
    je(mode)=ny-2
    kb(mode)=2
    ke(mode)=nz-2
  END DO
!
!-----------------------------------------------------------------------
!
!  Read all input data arrays
!
!-----------------------------------------------------------------------
!
  CALL dtaread(nx,ny,nz,nzsoil,nstyps,                                  &
               hdmpfmt,nchin,grdbasfn(1:lengbf),lengbf,                 &
               filename(1:lenfil),lenfil,time,                          &
               x,y,z,zp,zpsoil,uprt,vprt,wprt,ptprt,pprt,               &
               qvprt, qscalar, tke,kmh,kmv,                             &
               ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,            &
               soiltyp,stypfrct,vegtyp,lai,roufns,veg,                  &
               tsoil,qsoil,wetcanp,snowdpth,                            &
               raing,rainc,prcrate,                                     &
               radfrc,radsw,rnflx,radswnet,radlwin,                     &
               usflx,vsflx,ptsflx,qvsflx,                               &
               ireturn, tem1,tem2,tem3)

  curtim = time
!
!-----------------------------------------------------------------------
!
!  ireturn = 0 for a successful read
!
!-----------------------------------------------------------------------
!
  IF( ireturn == 0 ) THEN   ! successful read

!-----------------------------------------------------------------------
!
! Determine how many variables to process
!
!-----------------------------------------------------------------------

  READ(5,field_number,END=100)
  WRITE(6,'(/a,a)') 'Namelist block field_number successfully read.'

!
!-----------------------------------------------------------------------
!
!  Determine variables to plot
!
!-----------------------------------------------------------------------
!

  READ(5,plot_options,END=100)
  WRITE(6,'(/a,a)') 'Namelist block plot_options successfully read.'


  DO nm=1,nreply
    IF( ireply(nm) /= 1) THEN

      CALL strcap(inline(nm),inline(nm),80)
      WRITE(6,'(a,/a80)') '  inline: ',inline(nm)
      CALL parsln(inline(nm),params,80,maxpar,nparams)

      WRITE(6,'(a,i4,a,20a)')                                           &
          ' read ',nparams,' params: ', (params(ip),ip=1,nparams)
      WRITE(6,'(a)') ' '

      IF( kslice(nm) <= 0) kslice(nm) = (nz-2)/2+1
      IF( jslice(nm) <= 0) jslice(nm) = (ny-2)/2+1
      IF( islice(nm) <= 0) islice(nm) = (nx-2)/2+1
    ELSE
      nparams=5
      params(1)='uu'
      params(2)='vv'
      params(3)='ww'
      params(4)='TH'
      params(5)='PR'
!
!-----------------------------------------------------------------------
!
!  Set the printing parameter for each mode
!
!-----------------------------------------------------------------------
!
      iplot(1,nm)=1       ! control for plotting x-y slice
      iplot(2,nm)=1       ! control for plotting x-z slice
      iplot(3,nm)=0       ! control for plotting y-z slice

      kslice(nm) = 3
      jslice(nm) = (ny-2)/2+1
      islice(nm) = 8
    END IF
!
!-----------------------------------------------------------------------
!
!  Transfer the islice info into the ib,ie vectors
!  Likewise for j and k
!
!-----------------------------------------------------------------------
!
    kb(1)=kslice(nm)
    ke(1)=kslice(nm)

    jb(2)=jslice(nm)
    je(2)=jslice(nm)

    ib(3)=islice(nm)
    ie(3)=islice(nm)
!
!-----------------------------------------------------------------------
!
!  Calculate total fields from that for base state and perturbations
!
!-----------------------------------------------------------------------
!
    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          u(i,j,k)=uprt(i,j,k)+ubar(i,j,k)
          v(i,j,k)=vprt(i,j,k)+vbar(i,j,k)
          w(i,j,k)=wprt(i,j,k)+wbar(i,j,k)
          pt(i,j,k)=ptprt(i,j,k)+ptbar(i,j,k)
          qv(i,j,k)=qvprt(i,j,k)+qvbar(i,j,k)
        END DO
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Loop through three plotting modes
!  producing x-y, x-z and y-z slices of 2-d plotting each time.
!  If iplot(mode)=0, plotting is skipped.
!
!-----------------------------------------------------------------------
!
    DO mode = 1, 3
      WRITE(6,'(a,i3)') 'mode =',mode
      IF( iplot(mode,nm) /= 0 ) THEN
        DO ip=1,nparams
          WRITE(6,'(a,i4,a,a2)')                                        &
              '  ip= ',ip,'  param= ',params(ip)
          IF (params(ip) == 'uu') THEN
            WRITE(title,'(a6,1x,a)') runname(1:6),'u'
            CALL wrigar(u,                                              &
                        1,nx,1,ny,1,nz,ib(mode),ie(mode),               &
                        jb(mode),je(mode),kb(mode),ke(mode),            &
                        title,0.,mode)

          ELSE IF (params(ip) == 'vv') THEN
            WRITE(title,'(a6,1x,a)') runname(1:6),'v'
            CALL wrigar(v,                                              &
                        1,nx,1,ny,1,nz,ib(mode),ie(mode),               &
                        jb(mode),je(mode),kb(mode),ke(mode),            &
                        title,0.,mode)

          ELSE IF (params(ip) == 'ww') THEN
            WRITE(title,'(a6,1x,a)') runname(1:6),'w'
            CALL wrigar(w,                                              &
                        1,nx,1,ny,1,nz,ib(mode),ie(mode),               &
                        jb(mode),je(mode),kb(mode),ke(mode),            &
                        title,0.,mode)

          ELSE IF (params(ip) == 'pt') THEN
            WRITE(title,'(a6,1x,a)') runname(1:6),'pt'
            CALL wrigar(pt,                                             &
                        1,nx,1,ny,1,nz,ib(mode),ie(mode),               &
                        jb(mode),je(mode),kb(mode),ke(mode),            &
                        title,0.,mode)

          ELSE IF (params(ip) == 'pp') THEN
            WRITE(title,'(a6,1x,a)') runname(1:6),'pprt '
            CALL wrigar(pprt ,                                          &
                        1,nx,1,ny,1,nz,ib(mode),ie(mode),               &
                        jb(mode),je(mode),kb(mode),ke(mode),            &
                        title,0.,mode)

          ELSE IF (params(ip) == 'qv') THEN
            WRITE(title,'(a6,1x,a)') runname(1:6),'qv'
            CALL wrigar(qv,                                             &
                        1,nx,1,ny,1,nz,ib(mode),ie(mode),               &
                        jb(mode),je(mode),kb(mode),ke(mode),            &
                        title,0.,mode)

          ELSE IF (params(ip) == 'qc' .AND. P_QC > 0) THEN
            WRITE(title,'(a6,1x,a)') runname(1:6),'qc'
            CALL wrigar(qscalar(:,:,:,P_QC),                            &
                        1,nx,1,ny,1,nz,ib(mode),ie(mode),               &
                        jb(mode),je(mode),kb(mode),ke(mode),            &
                        title,0.,mode)

          ELSE IF (params(ip) == 'qr' .AND. P_QR > 0) THEN
            WRITE(title,'(a6,1x,a)') runname(1:6),'qr'
            CALL wrigar(qscalar(:,:,:,P_QR),                            &
                        1,nx,1,ny,1,nz,ib(mode),ie(mode),               &
                        jb(mode),je(mode),kb(mode),ke(mode),            &
                        title,0.,mode)

          ELSE IF (params(ip) == 'qi' .AND. P_QI > 0) THEN
            WRITE(title,'(a6,1x,a)') runname(1:6),'qi'
            CALL wrigar(qscalar(:,:,:,P_QI),                            &
                        1,nx,1,ny,1,nz,ib(mode),ie(mode),               &
                        jb(mode),je(mode),kb(mode),ke(mode),            &
                        title,0.,mode)

          ELSE IF (params(ip) == 'qs'  .AND. P_QS > 0) THEN
            WRITE(title,'(a6,1x,a)') runname(1:6),'qs'
            CALL wrigar(qscalar(:,:,:,P_QS),                            &
                        1,nx,1,ny,1,nz,ib(mode),ie(mode),               &
                        jb(mode),je(mode),kb(mode),ke(mode),            &
                        title,0.,mode)

          ELSE IF (params(ip) == 'qh' .AND. P_QH > 0) THEN
            WRITE(title,'(a6,1x,a)') runname(1:6),'qh'
            CALL wrigar(qscalar(:,:,:,P_QH),                            &
                        1,nx,1,ny,1,nz,ib(mode),ie(mode),               &
                        jb(mode),je(mode),kb(mode),ke(mode),            &
                        title,0.,mode)
!
!-----------------------------------------------------------------------
!
!  Perturbation quantities
!
!-----------------------------------------------------------------------
!

          ELSE IF (params(ip) == 'up') THEN
            WRITE(title,'(a6,1x,a)') runname(1:6),'uprt '
            CALL wrigar(uprt ,                                          &
                        1,nx,1,ny,1,nz,ib(mode),ie(mode),               &
                        jb(mode),je(mode),kb(mode),ke(mode),            &
                        title,0.,mode)

          ELSE IF (params(ip) == 'vp') THEN
            WRITE(title,'(a6,1x,a)') runname(1:6),'vprt '
            CALL wrigar(vprt ,                                          &
                        1,nx,1,ny,1,nz,ib(mode),ie(mode),               &
                        jb(mode),je(mode),kb(mode),ke(mode),            &
                        title,0.,mode)

          ELSE IF (params(ip) == 'wp') THEN
            WRITE(title,'(a6,1x,a)') runname(1:6),'wprt '
            CALL wrigar(wprt ,                                          &
                        1,nx,1,ny,1,nz,ib(mode),ie(mode),               &
                        jb(mode),je(mode),kb(mode),ke(mode),            &
                        title,0.,mode)

          ELSE IF (params(ip) == 'tp') THEN
            WRITE(title,'(a6,1x,a)') runname(1:6),'ptprt'
            CALL wrigar(ptprt,                                          &
                        1,nx,1,ny,1,nz,ib(mode),ie(mode),               &
                        jb(mode),je(mode),kb(mode),ke(mode),            &
                        title,0.,mode)

          ELSE IF (params(ip) == 'qp') THEN
            WRITE(title,'(a6,1x,a)') runname(1:6),'qvprt '
            CALL wrigar(qvprt ,                                         &
                        1,nx,1,ny,1,nz,ib(mode),ie(mode),               &
                        jb(mode),je(mode),kb(mode),ke(mode),            &
                        title,0.,mode)

          ELSE IF (params(ip) == 'ub') THEN
            WRITE(title,'(a6,1x,a)') runname(1:6),'u base'
            CALL wrigar(ubar,                                           &
                        1,nx,1,ny,1,nz,ib(mode),ie(mode),               &
                        jb(mode),je(mode),kb(mode),ke(mode),            &
                        title,0.,mode)

          ELSE IF (params(ip) == 'vb') THEN
            WRITE(title,'(a6,1x,a)') runname(1:6),'v base'
            CALL wrigar(vbar,                                           &
                        1,nx,1,ny,1,nz,ib(mode),ie(mode),               &
                        jb(mode),je(mode),kb(mode),ke(mode),            &
                        title,0.,mode)

          ELSE IF (params(ip) == 'wb') THEN
            WRITE(title,'(a6,1x,a)') runname(1:6),'w base'
            CALL wrigar(wbar,                                           &
                        1,nx,1,ny,1,nz,ib(mode),ie(mode),               &
                        jb(mode),je(mode),kb(mode),ke(mode),            &
                        title,0.,mode)

          ELSE IF (params(ip) == 'tb') THEN
            WRITE(title,'(a6,1x,a)') runname(1:6),'thbase'
            CALL wrigar(ptbar,                                          &
                        1,nx,1,ny,1,nz,ib(mode),ie(mode),               &
                        jb(mode),je(mode),kb(mode),ke(mode),            &
                        title,0.,mode)

          ELSE IF (params(ip) == 'pb') THEN
            WRITE(title,'(a6,1x,a)') runname(1:6),'pbar'
            CALL wrigar(pbar,                                           &
                        1,nx,1,ny,1,nz,ib(mode),ie(mode),               &
                        jb(mode),je(mode),kb(mode),ke(mode),            &
                        title,0.,mode)

          ELSE IF (params(ip) == 'qb') THEN
            WRITE(title,'(a6,1x,a)') runname(1:6),'qvbase'
            CALL wrigar(qvbar,                                          &
                        1,nx,1,ny,1,nz,ib(mode),ie(mode),               &
                        jb(mode),je(mode),kb(mode),ke(mode),            &
                        title,0.,mode)
          ELSE IF (params(ip) == 'st') THEN
            WRITE(title,'(a6,1x,a)') runname(1:6),'tsoil_sfc  t'
            CALL wrigar( tsoil(1,1,1,0),                                          &
                 1,nx,1,ny,1,1,ib(mode),ie(mode),                       &
                 jb(mode),je(mode),1,1,                                 &
                 title,0.,mode)

          ELSE IF (params(ip) == 'LT') THEN
            WRITE(title,'(a6,1x,a)') runname(1:6),'tsoil_soil t'
            CALL wrigar( tsoil(1,1,2,0),                                         &
                 1,nx,1,ny,1,1,ib(mode),ie(mode),                       &
                 jb(mode),je(mode),1,1,                                 &
                 title,0.,mode)

          ELSE IF (params(ip) == 'sm') THEN
            WRITE(title,'(a6,1x,a)') runname(1:6),'qsoil_wetsfc'
            CALL wrigar( qsoil(1,1,1,0),                                        &
                 1,nx,1,ny,1,1,ib(mode),ie(mode),                       &
                 jb(mode),je(mode),1,1,                                 &
                 title,0.,mode)

          ELSE IF (params(ip) == 'lm') THEN
            WRITE(title,'(a6,1x,a)') runname(1:6),'qsoil_wetdp '
            CALL wrigar( qsoil(1,1,2,0),                                &
                 1,nx,1,ny,1,1,ib(mode),ie(mode),                       &
                 jb(mode),je(mode),1,1,                                 &
                 title,0.,mode)

          ELSE IF (params(ip) == 'cm') THEN
            WRITE(title,'(a6,1x,a)') runname(1:6),'wetcan'
            CALL wrigar( wetcanp,                                       &
                 1,nx,1,ny,1,1,ib(mode),ie(mode),                       &
                 jb(mode),je(mode),1,1,                                 &
                 title,0.,mode)

          END IF
        END DO
      END IF                                  ! iplot(mode).NE.0
    END DO

  END DO  !  nm

  END IF                                       ! successful read

  GOTO 101

  100 CONTINUE

  WRITE(6,'(/a,a)') 'Error reading NAMELIST file. The program will abort.'

  101 CONTINUE

  STOP
END PROGRAM arpsprt
