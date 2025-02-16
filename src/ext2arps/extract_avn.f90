
!##################################################################
!##################################################################
!######                                                      ######
!######                   PROGRAM EXTRACT_AVN                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

PROGRAM extract_avn
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!  Ming Xue (August 2000)
!
!  MODIFICATION HISTORY:
!
!  Nov. 5, 2000. (Ming Xue)
!  Added variables at sigma = 0.995
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Input and output grid dimensions
!
!-----------------------------------------------------------------------

  INTEGER :: nx_ext, ny_ext, nz_ext ! dimensions of external data grid
  INCLUDE 'globcst.inc'
!
!-----------------------------------------------------------------------
!
!  NAMELIST variables for input data
!
!-----------------------------------------------------------------------
!

  CHARACTER (LEN=256) :: dir_extd       ! directory of external data

  CHARACTER (LEN=80) :: extdname       ! Prefix string of external file name.
                                       ! Used ONLY for ARPS data input
                                       ! The file name should look like:
                                       !      extdname.time_string

  INTEGER :: nextdfil            ! number of external data files to process
  CHARACTER (LEN=29) :: extdtime(50)   ! external data times
                                       ! format mm-dd-yyyy:hh:mm:ss/HHH:MM:SS
                                       !        _initial time_____ /forecast

  REAL :: latbgn,latend,lonbgn,lonend, del_lat, del_lon
  INTEGER :: ibgn,iend,jbgn,jend

  CHARACTER (LEN=80)  :: jobname
  CHARACTER (LEN=256) :: outdir
!
  NAMELIST /extdfile/ jobname,dir_extd,extdname,                        &
                      nextdfil,extdtime,                                &
                      latbgn,latend,lonbgn,lonend,outdir,filcmprs

  REAL, allocatable :: lat_ext(:,:)    ! external data latidude
  REAL, allocatable :: lon_ext(:,:)    ! external data longitude

  REAL, allocatable :: p_ext(:,:,:)    ! Pressure (Pascals)
  REAL, allocatable :: hgt_ext(:,:,:)  ! Height (m)
  REAL, allocatable :: t_ext(:,:,:)    ! Temperature (K)
  REAL, allocatable :: u_ext(:,:,:)    ! Eastward wind component
  REAL, allocatable :: v_ext(:,:,:)    ! Northward wind component
  REAL, allocatable :: qv_ext(:,:,:)   ! Specific humidity (kg/kg)
!  REAL, allocatable :: qc_ext(:,:,:)   ! Cloud H2O mixing ratio (kg/kg)
!  REAL, allocatable :: qr_ext(:,:,:)   ! Rain  H2O mixing ratio (kg/kg)
!  REAL, allocatable :: qi_ext(:,:,:)   ! Ice   H2O mixing ratio (kg/kg)
!  REAL, allocatable :: qs_ext(:,:,:)   ! Snow  H2O mixing ratio (kg/kg)
!  REAL, allocatable :: qh_ext(:,:,:)   ! Hail  H2O mixing ratio (kg/kg)

  REAL, allocatable :: qscalar_ext(:,:,:,:)

  REAL, allocatable :: tsfc_ext   (:,:)   ! Temperature at surface (K)
  REAL, allocatable :: tsoil_ext  (:,:)   ! Deep soil temperature (K)
  REAL, allocatable :: wetsfc_ext (:,:)   ! Surface soil moisture
  REAL, allocatable :: wetdp_ext  (:,:)   ! Deep soil moisture
  REAL, allocatable :: wetcanp_ext(:,:)   ! Canopy water amount
  REAL, allocatable :: snowdpth_ext(:,:)  ! Snow depth (m)

  REAL, allocatable :: trn_ext    (:,:)   ! External terrain (m)
  REAL, allocatable :: psfc_ext   (:,:)   ! Surface pressure (Pa)

  REAL, allocatable :: ugrd_ext   (:,:)   ! u at sigma=0.995 (m/s)
  REAL, allocatable :: vgrd_ext   (:,:)   ! v at sigma=0.995 (m/s)
  REAL, allocatable :: tgrd_ext   (:,:)   ! T at sigma=0.995 (K)
  REAL, allocatable :: rhgrd_ext  (:,:)   ! relative humidity at sigma=0.995 (K)
  REAL, allocatable :: ptgrd_ext  (:,:)   ! PT at sigma=0.995 (K)
  REAL, allocatable :: pmsl_ext   (:,:)   ! MSL pressure (Pa)
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=19) :: extdinit
  CHARACTER (LEN=9) :: extdfcst
  CHARACTER (LEN=9) :: julfname

  INTEGER :: i,istatus
  INTEGER :: iyr,imo,iday,ihr,imin,isec,jldy
  INTEGER :: ifhr,ifmin,ifsec,mfhr
  INTEGER :: myr,iabssec,jabssec
  INTEGER :: ifile
  INTEGER :: iextmn,iextmx,jextmn,jextmx
  INTEGER :: idiag,jdiag

  REAL :: latnot(2)
  REAL :: amin,amax
  REAL :: qvmin,qvmax,qvval
  REAL :: csconst,pconst
  REAL :: deltaz,tv_ext
  REAL :: pres,temp,qvsat,rh,tvbar,qvprt,qtot
  REAL :: xdiag,ydiag,dd,dmin,latd,lond
  REAL :: ppasc,pmb,tc,tdc,theta,smix,e,bige,alge,dir,spd

  CHARACTER (LEN=80) :: timsnd
  INTEGER :: lfn, tmstrln
  INTEGER :: tsfcout,tsoilout,wetsout,wetdout,wetcout,snowdout
  INTEGER :: isnow,jsnow,ii,jj
  REAL :: xumin,xumax,yvmin,yvmax

  INTEGER :: strlen, ireturn
  CHARACTER (LEN=3) :: FMT
  CHARACTER (LEN=100) :: tmp_ch
  INTEGER :: iproj_ext, nunit, idummy
  REAL :: scale_ext,trlon_ext,latnot_ext(2),x0_ext,y0_ext,rdummy
  INTEGER :: iuout, ivout, ipout, ihout,itout,                          &
          iqvout, itsfcout,itsoilout,iwsfcout,iwdpout,                  &
          iwcnpout,isnowout,itrnout,ipsfcout,                           &
          iugrdout,ivgrdout,itgrdout,iptgrdout,irhgrdout,ipmslout

  CHARACTER (LEN=13)  :: outtime
  CHARACTER (LEN=256) :: outfile
  INTEGER :: len1,outflen,lenrun,ierr
  LOGICAL :: iexist
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!
!-----------------------------------------------------------------------
!
!  Read in additional namelists for external file specifications.
!
!-----------------------------------------------------------------------
!
  dir_extd = './'
  extdname = 'may20'
  nextdfil = 1
  extdtime(1) = '1977-05-20.21:00:00+000:00:00'
  filcmprs = 1

  READ (5,extdfile)

  WRITE (6, '(3x,a,a,a)')'dir_extd   = ''', trim(dir_extd), ''','
  WRITE (6, '(3x,a,a,a)')'extdname   = ''', trim(extdname), ''','
  WRITE (6, '(3x,a,i4,a)')'nextdfil  = ', nextdfil, ','

  DO i=1,nextdfil
    WRITE (6, '(3x,a,i2.2,a,a,a)')                                      &
        'extdtime(',i,') = ''', trim(extdtime(i)), ''','
  END DO

  2  CONTINUE

  nx_ext = 360
  ny_ext = 181
  nz_ext = 26

  PRINT*,'nx_ext, ny_ext, nz_ext = ', nx_ext, ny_ext, nz_ext

  IF( lonbgn <= 0.0 ) THEN
    lonbgn = lonbgn + 360.0
  END IF

  IF( lonend <= 0.0) THEN
    lonend = lonend + 360.0
  END IF

  IF( lonbgn > lonend ) lonend = lonend + 360.0

  latbgn = latbgn + 90.0
  latend = latend + 90.0

  ibgn = nint(lonbgn/1.0)+1
  iend = nint(lonend/1.0)+1
  jbgn = nint(latbgn/1.0)+1
  jend = nint(latend/1.0)+1

  PRINT *,'ibgn= ',ibgn,' iend=', iend
  PRINT *,'jbgn= ',jbgn,' jend=', jend

  allocate(lat_ext(ibgn:iend,jbgn:jend),stat=istatus)
  allocate(lon_ext(ibgn:iend,jbgn:jend),stat=istatus)
!
  allocate(p_ext(ibgn:iend,jbgn:jend,nz_ext),stat=istatus)
  allocate(hgt_ext(ibgn:iend,jbgn:jend,nz_ext),stat=istatus)
  allocate(t_ext(ibgn:iend,jbgn:jend,nz_ext),stat=istatus)
  allocate(u_ext(ibgn:iend,jbgn:jend,nz_ext),stat=istatus)
  allocate(v_ext(ibgn:iend,jbgn:jend,nz_ext),stat=istatus)
  allocate(qv_ext(ibgn:iend,jbgn:jend,nz_ext),stat=istatus)
!  allocate(qc_ext(ibgn:iend,jbgn:jend,nz_ext),stat=istatus)
!  allocate(qr_ext(ibgn:iend,jbgn:jend,nz_ext),stat=istatus)
!  allocate(qi_ext(ibgn:iend,jbgn:jend,nz_ext),stat=istatus)
!  allocate(qs_ext(ibgn:iend,jbgn:jend,nz_ext),stat=istatus)
!  allocate(qh_ext(ibgn:iend,jbgn:jend,nz_ext),stat=istatus)

!  allocate(qscalar_ext(ibgn:iend,jbgn:jend,nz_ext,nscalar),stat=istatus)

  allocate(tsfc_ext   (ibgn:iend,jbgn:jend),stat=istatus)
  allocate(tsoil_ext  (ibgn:iend,jbgn:jend),stat=istatus)
  allocate(wetsfc_ext (ibgn:iend,jbgn:jend),stat=istatus)
  allocate(wetdp_ext  (ibgn:iend,jbgn:jend),stat=istatus)
  allocate(wetcanp_ext(ibgn:iend,jbgn:jend),stat=istatus)
  allocate(snowdpth_ext(ibgn:iend,jbgn:jend),stat=istatus)

  allocate(trn_ext    (ibgn:iend,jbgn:jend),stat=istatus)
  allocate(psfc_ext   (ibgn:iend,jbgn:jend),stat=istatus)

  allocate(ugrd_ext   (ibgn:iend,jbgn:jend),stat=istatus)
  allocate(vgrd_ext   (ibgn:iend,jbgn:jend),stat=istatus)
  allocate(tgrd_ext   (ibgn:iend,jbgn:jend),stat=istatus)
  allocate(rhgrd_ext  (ibgn:iend,jbgn:jend),stat=istatus)
  allocate(ptgrd_ext  (ibgn:iend,jbgn:jend),stat=istatus)
  allocate(pmsl_ext   (ibgn:iend,jbgn:jend),stat=istatus)

  lat_ext=-999.0
  lon_ext=-999.0
  p_ext  =-999.0
  hgt_ext=-999.0
  t_ext  =-999.0
  u_ext  =-999.0
  v_ext  =-999.0
  qv_ext =-999.0
!  qc_ext =-999.0
!  qr_ext =-999.0
!  qi_ext =-999.0
!  qs_ext =-999.0
!  qh_ext =-999.0
!  qscalar_ext=-999.0
  tsfc_ext    =-999.0
  tsoil_ext   =-999.0
  wetsfc_ext  =-999.0
  wetdp_ext   =-999.0
  wetcanp_ext =-999.0
  snowdpth_ext=-999.0
  trn_ext     =-999.0
  psfc_ext    =-999.0

  ugrd_ext    =-999.0
  vgrd_ext    =-999.0
  tgrd_ext    =-999.0
  rhgrd_ext   =-999.0
  ptgrd_ext   =-999.0
  pmsl_ext    =-999.0
!
!-----------------------------------------------------------------------
!
!  Loop through the data times provided via NAMELIST.
!
!-----------------------------------------------------------------------
!
  DO ifile=1,nextdfil
!
!-----------------------------------------------------------------------
!
!  Time conversions.
!  Formats:  extdtime='1994-05-06.18:00:00+000:00:00'
!            julfname='941261800'
!
!-----------------------------------------------------------------------
!
    READ(extdtime(ifile),'(a19,1x,a9)') extdinit,extdfcst
    IF(extdfcst == '         ') extdfcst='000:00:00'
    READ(extdinit,                                                      &
        '(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)',ERR=920,END=920)           &
        iyr,imo,iday,ihr,imin,isec
    CALL julday(iyr,imo,iday,jldy)
    myr=MOD(iyr,100)
    ifhr=0
    ifmin=0
    ifsec=0
    READ(extdfcst,                                                      &
        '(i3,1x,i2,1x,i2)',ERR=4,END=4) ifhr,ifmin,ifsec
    4    CONTINUE
    mfhr=MOD(ifhr,24)
    jldy = jldy + ifhr/24
    WRITE(julfname,                                                     &
         '(i2.2,i3.3,i2.2,i2.2)') myr,jldy,ihr,mfhr
    CALL ctim2abss(iyr,imo,iday,ihr,imin,isec,iabssec)
    jabssec=(ifhr*3600) + (ifmin*60) + ifsec + iabssec

!    write(6,'(a,a9,a,/19x,a,a19,a/a,a/a,i16,a,/,i26,a)')
!    :     ' Calling rdextfil, looking for ',
!    :       extdfcst,' hour forecast '
!
!-----------------------------------------------------------------------
!
!  Get NCEP AVN GRIB #3
!
!-----------------------------------------------------------------------
!

    CALL get_avn_grb(nx_ext,ny_ext,nz_ext,ibgn,iend,jbgn,jend,          &
         dir_extd,extdname,extdinit,extdfcst,julfname,                  &
         iproj_ext,scale_ext,trlon_ext,latnot_ext,x0_ext,y0_ext,        &
         lat_ext,lon_ext,p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,        &
         tsfc_ext,tsoil_ext,wetsfc_ext,wetdp_ext,wetcanp_ext,           &
         snowdpth_ext,trn_ext,psfc_ext,                                 &
         ugrd_ext,vgrd_ext,tgrd_ext,rhgrd_ext,ptgrd_ext,pmsl_ext,       &
         istatus)

    IF(istatus /= 1) GO TO 999

    PRINT*,' '
    CALL a3dmax0(lat_ext,ibgn,iend,ibgn,iend,jbgn,jend,jbgn,jend,       &
               1,1,1,1,amax,amin)
    WRITE(6,'(1x,2(a,e13.6))')                                          &
          'lat_ext_min= ', amin,', lat_ext_max=',amax
    CALL a3dmax0(lon_ext,ibgn,iend,ibgn,iend,jbgn,jend,jbgn,jend,       &
               1,1,1,1,amax,amin)
    WRITE(6,'(1x,2(a,e13.6))')                                          &
          'lon_ext_min= ', amin,', lon_ext_max=',amax
!
    CALL a3dmax0(p_ext  ,ibgn,iend,ibgn,iend,jbgn,jend,jbgn,jend,       &
               1,nz_ext,1,nz_ext,amax,amin)
    WRITE(6,'(1x,2(a,e13.6))')                                          &
          'p_ext_min  = ', amin,', p_ext_max  =',amax
    CALL a3dmax0(hgt_ext,ibgn,iend,ibgn,iend,jbgn,jend,jbgn,jend,       &
               1,nz_ext,1,nz_ext,amax,amin)
    WRITE(6,'(1x,2(a,e13.6))')                                          &
          'hgt_ext_min= ', amin,', hgt_ext_max=',amax
    CALL a3dmax0(t_ext  ,ibgn,iend,ibgn,iend,jbgn,jend,jbgn,jend,       &
               1,nz_ext,1,nz_ext,amax,amin)
    WRITE(6,'(1x,2(a,e13.6))')                                          &
          't_ext_min  = ', amin,', t_ext_max  =',amax
    CALL a3dmax0(u_ext  ,ibgn,iend,ibgn,iend,jbgn,jend,jbgn,jend,       &
               1,nz_ext,1,nz_ext,amax,amin)
    WRITE(6,'(1x,2(a,e13.6))')                                          &
          'u_ext_min  = ', amin,', u_ext_max  =',amax
    CALL a3dmax0(v_ext  ,ibgn,iend,ibgn,iend,jbgn,jend,jbgn,jend,       &
               1,nz_ext,1,nz_ext,amax,amin)
    WRITE(6,'(1x,2(a,e13.6))')                                          &
          'v_ext_min  = ', amin,', v_ext_max  =',amax
    CALL a3dmax0(qv_ext ,ibgn,iend,ibgn,iend,jbgn,jend,jbgn,jend,       &
               1,nz_ext,1,nz_ext,amax,amin)
    WRITE(6,'(1x,2(a,e13.6))')                                          &
          'qv_ext_min = ', amin,', qv_ext_max =',amax

!    DO nq=1,nscalar
!      CALL a3dmax0(qscalar_ext(1,1,1,nq) ,ibgn,iend,ibgn,iend,jbgn,jend,jbgn,jend,       &
!                 1,nz_ext,1,nz_ext,amax,amin)
!      WRITE(6,'(1x,a,e13.6,2a,e13.6)')                                          &
!            TRIM(qnames(nq))//'_ext_min = ', amin,', ',TRIM(qnames(nq))//'_ext_max =',amax
!
!    END DO

    CALL a3dmax0(tsfc_ext,ibgn,iend,ibgn,iend,jbgn,jend,jbgn,jend,      &
               1,1,1,1,amax,amin)
    WRITE(6,'(1x,2(a,e13.6))')                                          &
          'tsfc_ext_min   = ', amin,', tsfc_ext_max   =',amax
    CALL a3dmax0(tsoil_ext,ibgn,iend,ibgn,iend,jbgn,jend,jbgn,jend,     &
               1,1,1,1,amax,amin)
    WRITE(6,'(1x,2(a,e13.6))')                                          &
          'tsoil_ext_min  = ', amin,', tsoil_ext_max  =',amax
    CALL a3dmax0(wetsfc_ext,ibgn,iend,ibgn,iend,jbgn,jend,jbgn,jend,    &
               1,1,1,1,amax,amin)
    WRITE(6,'(1x,2(a,e13.6))')                                          &
          'wetsfc_ext_min = ', amin,', wetsfc_ext_max =',amax
    CALL a3dmax0(wetdp_ext,ibgn,iend,ibgn,iend,jbgn,jend,jbgn,jend,     &
               1,1,1,1,amax,amin)
    WRITE(6,'(1x,2(a,e13.6))')                                          &
          'wetdp_ext_min  = ', amin,', wetdp_ext_max  =',amax
    CALL a3dmax0(wetcanp_ext,                                           &
                 ibgn,iend,ibgn,iend,jbgn,jend,jbgn,jend,               &
                 1,1,1,1,amax,amin)
    WRITE(6,'(1x,2(a,e13.6))')                                          &
          'wetcanp_ext_min= ', amin,', wetcanp_ext_max=',amax
    CALL a3dmax0(snowdpth_ext,                                          &
               ibgn,iend,ibgn,iend,jbgn,jend,jbgn,jend,                 &
               1,1,1,1,amax,amin)
    WRITE(6,'(1x,2(a,e13.6))')                                          &
          'snowd_ext_min  = ', amin,', snow_ext_max   =',amax
    CALL a3dmax0(psfc_ext,                                              &
               ibgn,iend,ibgn,iend,jbgn,jend,jbgn,jend,                 &
               1,1,1,1,amax,amin)
    WRITE(6,'(1x,2(a,e13.6))')                                          &
          'psfc_ext_min   = ', amin,', psfc_ext_max   =',amax

    CALL a3dmax0(trn_ext,ibgn,iend,ibgn,iend,jbgn,jend,jbgn,jend,       &
               1,1,1,1,amax,amin)
    WRITE(6,'(1x,2(a,e13.6))')                                          &
          'trn_ext_min    = ', amin,', trn_ext_max    =',amax

    CALL a3dmax0(ugrd_ext,ibgn,iend,ibgn,iend,jbgn,jend,jbgn,jend,      &
               1,1,1,1,amax,amin)
    WRITE(6,'(1x,2(a,e13.6))')                                          &
          'ugrd_ext_min   = ', amin,', ugrd_ext_max   =',amax
    CALL a3dmax0(vgrd_ext,ibgn,iend,ibgn,iend,jbgn,jend,jbgn,jend,      &
               1,1,1,1,amax,amin)
    WRITE(6,'(1x,2(a,e13.6))')                                          &
          'vgrd_ext_min   = ', amin,', vgrd_ext_max   =',amax
    CALL a3dmax0(tgrd_ext,ibgn,iend,ibgn,iend,jbgn,jend,jbgn,jend,      &
               1,1,1,1,amax,amin)
    WRITE(6,'(1x,2(a,e13.6))')                                          &
          'tgrd_ext_min   = ', amin,', tgrd_ext_max   =',amax
    CALL a3dmax0(rhgrd_ext,ibgn,iend,ibgn,iend,jbgn,jend,jbgn,jend,     &
               1,1,1,1,amax,amin)
    WRITE(6,'(1x,2(a,e13.6))')                                          &
          'rhgrd_ext_min  = ', amin,', rhgrd_ext_max  =',amax
    CALL a3dmax0(ptgrd_ext,ibgn,iend,ibgn,iend,jbgn,jend,jbgn,jend,     &
               1,1,1,1,amax,amin)
    WRITE(6,'(1x,2(a,e13.6))')                                          &
          'ptgrd_ext_min  = ', amin,', ptgrd_ext_max  =',amax
    CALL a3dmax0( pmsl_ext,ibgn,iend,ibgn,iend,jbgn,jend,jbgn,jend,     &
               1,1,1,1,amax,amin)
    WRITE(6,'(1x,2(a,e13.6))')                                          &
          ' pmsl_ext_min  = ', amin,',  pmsl_ext_max  =',amax

!
!-----------------------------------------------------------------------
!
!  Write out the data patch
!
!-----------------------------------------------------------------------
!
    READ (extdinit,'(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)')                &
         iyr,imo,iday,ihr,imin,isec

    myr=MOD(iyr,100)
    ifhr=0
    ifmin=0
    ifsec=0

    READ(extdfcst,'(i3)',ERR=5,END=5) ifhr
    5     CONTINUE

    WRITE (outtime,'(i4.4,i2.2,i2.2,i2.2,a1,i2.2)')                     &
           iyr,imo,iday,ihr,'f',ifhr

    len1=len_trim(outdir)
    outflen=len1

    IF( outflen == 0 .OR. outdir(1:outflen) == ' ' ) THEN
      outdir = '.'
      outflen=1
    END IF

    outflen = len_trim( outdir )
    IF( outdir(outflen:outflen) /= '/') THEN
      outflen=outflen+1
      outdir(outflen:outflen)='/'
    END IF

!  print*,'outdir =', trim(outdir)

    CALL INQUIREDIR(trim(outdir),iexist)

    IF( .NOT.iexist ) THEN

      WRITE(6,'(5x,a,2(/5x,a))')                                        &
          'Specified output directory '//trim(outdir)//' not found.',   &
          'It was created by the program.'
      CALL unixcmd( 'mkdir -p '//trim(outdir) )

    END IF

    lenrun = len_trim( jobname )

    outfile = outdir(1:outflen)//jobname(1:lenrun)                      &
                                  //'.'//outtime(1:13)
    outflen = outflen + lenrun + 14
    CALL fnversn(outfile, outflen)

    CALL getunit( nunit)
    CALL asnctl ('NEWLOCAL', 1, ierr)
    CALL asnfile(outfile(1:outflen), '-F f77 -N ieee', ierr)

    OPEN(UNIT=nunit,FILE=outfile(1:outflen),                            &
         STATUS='unknown',FORM='unformatted',IOSTAT=istatus)

    IF( istatus /= 0 ) THEN
      WRITE(6,'(1x,a,a,/1x,i3,a)')                                      &
          'Error occured when opening file ',outfile(1:outflen),        &
          'using FORTRAN unit ',nunit,' Program stopped.'
      STOP
    END IF

    PRINT*,'To creat file ',trim(outfile)

    del_lat = 1.0
    del_lon = 1.0

    WRITE(nunit) iend-ibgn+1, jend-jbgn+1,nz_ext
    WRITE(nunit) lonbgn,lonend,latbgn,latend
    WRITE(nunit) del_lon, del_lat

    iuout    = 1
    ivout    = 1
    ipout    = 1
    ihout    = 1
    itout    = 1
    iqvout   = 1
    itsfcout = 1
    itsoilout= 1
    iwsfcout = 1
    iwdpout  = 1
    iwcnpout = 0
    isnowout = 0
    itrnout  = 1
    ipsfcout = 1

    iugrdout = 1
    ivgrdout = 1
    itgrdout = 1
    iptgrdout= 1
    ipmslout = 1
    irhgrdout= 1

    idummy = 0.0
    WRITE(nunit) iuout, ivout, ipout, ihout,itout,                      &
          iqvout, itsfcout,itsoilout,iwsfcout,iwdpout,                  &
          iwcnpout,isnowout,itrnout,ipsfcout,iugrdout,                  &
          ivgrdout,itgrdout,iptgrdout,irhgrdout,ipmslout,               &
          iproj_ext,idummy,idummy, idummy, idummy,                      &
          idummy, idummy,  idummy, idummy, idummy

    rdummy = 0.0
    WRITE(nunit) scale_ext,trlon_ext,latnot_ext(1),latnot_ext(2),       &
          x0_ext, y0_ext, rdummy, rdummy, rdummy, rdummy,               &
          rdummy, rdummy, rdummy, rdummy, rdummy,                       &
          rdummy, rdummy, rdummy, rdummy, rdummy,                       &
          rdummy, rdummy, rdummy, rdummy, rdummy,                       &
          rdummy, rdummy, rdummy, rdummy, rdummy

    IF( iuout==1 ) WRITE(nunit) 'u(m/s)....'
    IF( iuout==1 ) WRITE(nunit) u_ext
    IF( ivout==1 ) WRITE(nunit) 'v(m/s)....'
    IF( ivout==1 ) WRITE(nunit) v_ext
    IF( ipout==1 ) WRITE(nunit) 'p(P)......'
    IF( ipout==1 ) WRITE(nunit) p_ext
    IF( ihout==1 ) WRITE(nunit) 'hgt(m)....'
    IF( ihout==1 ) WRITE(nunit) hgt_ext
    IF( itout==1 ) WRITE(nunit) 't(K)......'
    IF( itout==1 ) WRITE(nunit) t_ext
    IF( iqvout==1) WRITE(nunit) 'qv(g/g)...'
    IF( iqvout==1) WRITE(nunit) qv_ext

    IF( itsfcout==1 ) WRITE(nunit) 'tsfc(K)...'
    IF( itsfcout==1 ) WRITE(nunit) tsfc_ext
    IF( itsoilout==1) WRITE(nunit) 'tsoil(K)..'
    IF( itsoilout==1) WRITE(nunit) tsoil_ext
    IF( iwsfcout==1 ) WRITE(nunit) 'wetsfc....'
    IF( iwsfcout==1 ) WRITE(nunit) wetsfc_ext
    IF( iwdpout==1  ) WRITE(nunit) 'wetdp.....'
    IF( iwdpout==1  ) WRITE(nunit) wetdp_ext

    IF( iwcnpout==1 ) WRITE(nunit) 'wetcanp...'
    IF( iwcnpout==1 ) WRITE(nunit) wetcanp_ext
    IF( isnowout==1 ) WRITE(nunit) 'snowdp(m).'
    IF( isnowout==1 ) WRITE(nunit) snowdpth_ext
    IF( itrnout ==1 ) WRITE(nunit) 'trn(m)....'
    IF( itrnout ==1 ) WRITE(nunit) trn_ext
    IF( ipsfcout==1 ) WRITE(nunit) 'psfc(Pa)..'
    IF( ipsfcout==1 ) WRITE(nunit) psfc_ext
    IF( iugrdout==1 ) WRITE(nunit) 'ugrd(m/s).'
    IF( iugrdout==1 ) WRITE(nunit)  ugrd_ext
    IF( ivgrdout==1 ) WRITE(nunit) 'vgrd(m/s).'
    IF( ivgrdout==1 ) WRITE(nunit)  vgrd_ext
    IF( itgrdout==1 ) WRITE(nunit) 'tgrd(K)...'
    IF( itgrdout==1 ) WRITE(nunit)  tgrd_ext
    IF( irhgrdout==1) WRITE(nunit) 'rhgrd.....'
    IF( irhgrdout==1) WRITE(nunit)  rhgrd_ext
    IF( iptgrdout==1) WRITE(nunit) 'ptgrd(K)..'
    IF( iptgrdout==1) WRITE(nunit)  ptgrd_ext
    IF( ipmslout==1 ) WRITE(nunit) 'pmsl(Pa)..'
    IF( ipmslout==1 ) WRITE(nunit)  pmsl_ext

    CLOSE (UNIT=nunit)
    CALL retunit( nunit)
    IF( filcmprs .eq. 1 ) CALL cmprs( outfile(1:outflen) )

    PRINT*,'Finished writing file ',trim(outfile)
    PRINT*,' '

  END DO
!
!-----------------------------------------------------------------------
!
!  Friendly exit message.
!
!-----------------------------------------------------------------------
!
  WRITE(6,'(a/a,i4,a)')                                                 &
          ' Normal succesful completion of EXT2ARPS.',                  &
          ' Processed',nextdfil,' file(s)'
  STOP
!
!-----------------------------------------------------------------------
!
!  Problem doing time conversions.
!
!-----------------------------------------------------------------------
!
  920 CONTINUE
  WRITE(6,'(a,/,a,i4,a,i4/,a,a19)')                                     &
          ' Aborting, error in time format for external file',          &
          ' File number:',ifile,' of',nextdfil,                         &
          ' External file time provided:',extdtime(ifile)
  STOP
!
!-----------------------------------------------------------------------
!
!  Error status returned from rdextfil
!
!-----------------------------------------------------------------------
!
  999 CONTINUE
  WRITE(6,'(a,i6)')                                                     &
          ' Aborting, error reading external file. istatus=',           &
            istatus
  STOP
END PROGRAM extract_avn
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GET_AVN_GRB                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE get_avn_grb(nx_ext,ny_ext,nz_ext,ibgn,iend,jbgn,jend,        &
           dir_extd,extdname,extdinit,extdfcst,julfname,                &
           iproj_ext,scale_ext,trlon_ext,latnot_ext,x0_ext,y0_ext,      &
           lat_ext,lon_ext,p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,      &
           tsfc_ext,tsoil_ext,wetsfc_ext,wetdp_ext,wetcanp_ext,         &
           snowdpth_ext,trn_ext,psfc_ext,                               &
           ugrd_ext,vgrd_ext,tgrd_ext,rhgrd_ext,ptgrd_ext,pmsl_ext,     &
           istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Reads and pass out a section of NCEP AVN GRIB
!  (Grid #3, 1x1 degree) data file.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: M. Xue
!  07/25/2000 Based on GETNMCAVN3.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    dir_extd      Directory name for external file
!    extdname      Prefix string of external file name
!    extdinit      Initialized time in mm-dd-yyyy:hh:mm:ss format
!    extdfcst      Forecast hour in HHH:MM:SS format
!    julfname      File name in yyjjjhhmm format
!
!  OUTPUT:
!
!    iproj_ext     Map projection number of external data
!    scale_ext     Scale factor of external data
!    trlon_ext     True longitude of external data (degrees E)
!    latnot_ext(2) True latitude(s) of external data (degrees N)
!    x0_ext        x coordinate of origin of external data
!    y0_ext        y coordinate of origin of external data
!    lat_ext       latitude of external data points (degrees N)
!    lon_ext       longitude of external data points (degrees E)
!    p_ext         pressure (Pascal)
!    hgt_ext       height (m)
!    t_ext         temperature (K)
!    qv_ext        specific humidity (kg/kg)
!    u_ext         u wind component (m/s)
!    v_ext         v wind component (m/s)
!    qc_ext        Cloud water mixing ratio (kg/kg)
!    qr_ext        Rain water mixing ratio (kg/kg)
!    qi_ext        Ice mixing ratio (kg/kg)
!    qs_ext        Snow mixing ratio (kg/kg)
!    qh_ext        Hail mixing ratio (kg/kg)
!
!    tsfc_ext      Surface temperature
!    tsoil_ext     Soil temperature
!    wetsfc_ext    Top layer soil moisture
!    wetdp_ext     Deep soil moisture
!    wetcanp_ext   Water content on canopy
!
!    trn_ext       External terrain (m)
!    psfc_ext      Surface pressure (Pa)
!    ugrd_ext	   u at sigma=0.995 (m/s)
!    vgrd_ext	   v at sigma=0.995 (m/s)
!    tgrd_ext	   T at sigma=0.995 (K)
!    rhgrd_ext     relative humidity at sigma=0.995 (K)
!    ptgrd_ext     PT at sigma=0.995 (K)
!    pmsl_ext      MSL pressure (Pa)
!
!    istatus       status indicator
!
!  WORK ARRAYS:
!
!    var_grb3d     Arrays to store the GRIB 3-D variables:
!                  var_grb3d(nxgrb,nygrb,nzgrb,1,1) - Temperature (K)
!                  var_grb3d(nxgrb,nygrb,nzgrb,2,1) - Specific humidity
!                                                     (kg/kg)
!                  var_grb3d(nxgrb,nygrb,nzgrb,3,1) - u wind (m/s)
!                  var_grb3d(nxgrb,nygrb,nzgrb,4,1) - v wind (m/s)
!                  var_grb3d(nxgrb,nygrb,nzgrb,5,1) - Geopotential
!                                                     height (gpm)
!                  var_grb3d(nxgrb,nygrb,nzgrb,6,1) - Pressure vertical
!                                                     velocity (Pa/s)
!                                                     (if applied)
!                  var_grb3d(nxgrb,nygrb,nzgrb,1,2) - soil temp. (K)
!                  var_grb3d(nxgrb,nygrb,nzgrb,2,2) - vol. soil moist.
!                                                     (m**3/m**3)
!
!    var_grb2d     Arrays to store the GRIB 2-D variables:
!                  var_grb2d(nxgrb,nygrb,1,1) - Surface pressure (Pa)
!                  var_grb2d(nxgrb,nygrb,2,1) - Geopotential height (gpm)
!                  var_grb2d(nxgrb,nygrb,3,1) - Surface temperature (K)
!                  var_grb2d(nxgrb,nygrb,4,1) - Plant canopy surface water (kg/m**2)
!                  var_grb2d(nxgrb,nygrb,6,1) - Snow depth
!
!                  var_grb2d(nxgrb,nygrb,1,2) - temperature at sigma = 0.995 (K)
!                  var_grb2d(nxgrb,nygrb,2,2) - u, east-west velocity at sigma = 0.995 (m/s)
!                  var_grb2d(nxgrb,nygrb,3,2) - v, north-south velocity at sigma = 0.995 (m/s)
!                  var_grb2d(nxgrb,nygrb,4,2) - relative humidity at sigma = 0.995
!                  var_grb2d(nxgrb,nygrb,5,2) - potential temperature at sigma = 0.995 (K)
!                  var_grb2d(nxgrb,nygrb,6,2) - undefined, unused
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INCLUDE 'gribcst.inc'

  CHARACTER (LEN=*) :: dir_extd
  CHARACTER (LEN=*) :: extdname

  CHARACTER (LEN=19) :: extdinit
  CHARACTER (LEN=9) :: extdfcst
  CHARACTER (LEN=9) :: julfname
!
!-----------------------------------------------------------------------
!
!  External grid variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iproj_ext
  REAL :: scale_ext,trlon_ext
  REAL :: latnot_ext(2)
  REAL :: x0_ext,y0_ext
  REAL :: dx_ext,dy_ext
!
!-----------------------------------------------------------------------
!
!  Output external variable arrays
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx_ext,ny_ext,nz_ext
  INTEGER :: ibgn,iend,jbgn,jend

  REAL :: lat_ext(ibgn:iend,jbgn:jend)
  REAL :: lon_ext(ibgn:iend,jbgn:jend)
  REAL :: p_ext  (ibgn:iend,jbgn:jend,nz_ext)   ! Pressure (Pascals)
  REAL :: hgt_ext(ibgn:iend,jbgn:jend,nz_ext)   ! Height (m)
  REAL :: t_ext  (ibgn:iend,jbgn:jend,nz_ext)   ! Temperature (K)
  REAL :: qv_ext (ibgn:iend,jbgn:jend,nz_ext)   ! Specific humidity (kg/kg)
  REAL :: u_ext  (ibgn:iend,jbgn:jend,nz_ext)   ! Eastward wind component
  REAL :: v_ext  (ibgn:iend,jbgn:jend,nz_ext)   ! Northward wind component
!  REAL :: qc_ext (ibgn:iend,jbgn:jend,nz_ext)   ! Cloud H2O mixing ratio (kg/kg)
!  REAL :: qr_ext (ibgn:iend,jbgn:jend,nz_ext)   ! Rain  H2O mixing ratio (kg/kg)
!  REAL :: qi_ext (ibgn:iend,jbgn:jend,nz_ext)   ! Ice   mixing ratio (kg/kg)
!  REAL :: qs_ext (ibgn:iend,jbgn:jend,nz_ext)   ! Snow  mixing ratio (kg/kg)
!  REAL :: qh_ext (ibgn:iend,jbgn:jend,nz_ext)   ! Hail  mixing ratio (kg/kg)

  REAL :: tsfc_ext   (ibgn:iend,jbgn:jend)      ! Temperature at surface (K)
  REAL :: tsoil_ext  (ibgn:iend,jbgn:jend)      ! Deep soil temperature (K)
  REAL :: wetsfc_ext (ibgn:iend,jbgn:jend)      ! Surface soil moisture
  REAL :: wetdp_ext  (ibgn:iend,jbgn:jend)      ! Deep soil moisture
  REAL :: wetcanp_ext(ibgn:iend,jbgn:jend)      ! Canopy water amount
  REAL :: snowdpth_ext(ibgn:iend,jbgn:jend)     ! Snow depth (m)

  REAL :: trn_ext    (ibgn:iend,jbgn:jend)      ! External terrain (m)
  REAL :: psfc_ext   (ibgn:iend,jbgn:jend)      ! Surface pressure (Pa)

  REAL :: ugrd_ext   (ibgn:iend,jbgn:jend)      ! u at sigma=0.995 (m/s)
  REAL :: vgrd_ext   (ibgn:iend,jbgn:jend)      ! v at sigma=0.995 (m/s)
  REAL :: tgrd_ext   (ibgn:iend,jbgn:jend)      ! T at sigma=0.995 (K)
  REAL :: rhgrd_ext  (ibgn:iend,jbgn:jend)      ! relative humidity at sigma=0.995 (K)
  REAL :: ptgrd_ext  (ibgn:iend,jbgn:jend)      ! PT at sigma=0.995 (K)
  REAL :: pmsl_ext   (ibgn:iend,jbgn:jend)      ! MSL pressure (Pa)
!
!-----------------------------------------------------------------------
!
!  Other  external variable arrays
!
!-----------------------------------------------------------------------
!
!  real x_ext(nx_ext)
!  real y_ext(ny_ext)

  INTEGER :: istatus
!
!-----------------------------------------------------------------------
!
!  Work arrays for storing grib data
!
!-----------------------------------------------------------------------
!
  REAL, allocatable :: var_grb2d(:,:,:,:)   ! GRIB variables
  REAL, allocatable :: var_grb3d(:,:,:,:,:) ! GRIB 3-D variables
  INTEGER, allocatable :: var_lev3d(:,:,:)  ! Levels (hybrid) for
                                            ! each 3-D variable
  REAL, allocatable :: rcdata(:)            ! temporary data array
!
!-----------------------------------------------------------------------
!
!  Original grid variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iproj
  REAL :: scale,trlon,x0,y0
  REAL :: latnot(2)
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: gribfile
  CHARACTER (LEN=13)  :: gribtime
  INTEGER :: i,j,k,kk
  INTEGER :: iyr,imo,iday,myr, jldy
  INTEGER :: ihr,imin,isec
  INTEGER :: ifhr,ifmin,ifsec
  INTEGER :: grbflen, len1, lenrun

  INTEGER :: m,n,nz1,max_nr2d,max_nr3d,min_nr3d,nz2

  REAL :: govrd

  INTEGER :: chklev, lvscan

  INTEGER :: iret             ! Return flag
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
!
!-----------------------------------------------------------------------
!
!  GRIB grid information
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=42) :: gridesc ! Grid description

  INTEGER :: iproj_grb    ! Map projection indicator
  INTEGER :: gthin        ! Indicator of whether the grid is "thinned"

  INTEGER :: ni_grb       ! Number of points along x-axis
  INTEGER :: nj_grb       ! Number of points along y-axis
  INTEGER :: np_grb       ! Total number of horizontal grid points

  INTEGER :: nk_grb       ! Number of vertical parameters
  REAL :: zk_grb(nz_ext)  ! Vertical coordinate parameters

  INTEGER :: npeq         ! Number of lat circles from pole to equator
  INTEGER :: nit(nz_ext)  ! Number of x-points for thinned grid

  REAL :: pi_grb          ! x-coordinate of pole point
  REAL :: pj_grb          ! y-coordinate of pole point
  INTEGER :: ipole        ! Projection center flag

  REAL :: di_grb          ! x-direction increment or grid length
  REAL :: dj_grb          ! y-direction increment or grid length

  REAL :: latsw           ! Latitude  of South West corner point
  REAL :: lonsw           ! Longitude of South West corner point
  REAL :: latne           ! Latitude  of North East corner point
  REAL :: lonne           ! Longitude of North East corner point

  REAL :: lattru1         ! Latitude (1st) at which projection is true
  REAL :: lattru2         ! Latitude (2nd) at which projection is true
  REAL :: lontrue         ! Longitude      at which projection is true

  REAL :: latrot          ! Latitude  of southern pole of rotation
  REAL :: lonrot          ! Longitude of southern pole of rotation
  REAL :: angrot          ! Angle of rotation

  REAL :: latstr          ! Latitude  of the pole of stretching
  REAL :: lonstr          ! Longitude of the pole of stretching
  REAL :: facstr          ! Stretching factor

  INTEGER :: scanmode     ! Scanning indicator
  INTEGER :: iscan        ! x-direction   scanning indicator
  INTEGER :: jscan        ! y-direction   scanning indicator
  INTEGER :: kscan        ! FORTRAN index scanning indicator

  INTEGER :: ires         ! Resolution direction increments indicator
  INTEGER :: iearth       ! Earth shape indicator: spherical or oblate?
  INTEGER :: icomp        ! (u,v) components decomposition indicator

  INTEGER :: jpenta       ! J-Pentagonal resolution parameter
  INTEGER :: kpenta       ! K-Pentagonal resolution parameter
  INTEGER :: mpenta       ! M-Pentagonal resolution parameter
  INTEGER :: ispect       ! Spectral representation type
  INTEGER :: icoeff       ! Spectral coefficient storage mode

  REAL :: xp_grb          ! X coordinate of sub-satellite point
  REAL :: yp_grb          ! Y coordinate of sub-satellite point
  REAL :: xo_grb          ! X coordinate of image sector origin
  REAL :: yo_grb          ! Y coordinate of image sector origin
  REAL :: zo_grb          ! Camera altitude from center of Earth

  INTEGER :: isrc
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  allocate(var_grb2d(nx_ext,ny_ext,n2dvs,n2dlvt))
  allocate(var_grb3d(nx_ext,ny_ext,nz_ext,n3dvs,n3dlvt))
  allocate(rcdata(nx_ext*ny_ext))
  allocate(var_lev3d(nz_ext,n3dvs,n3dlvt))

  IF(extdfcst == '         ') extdfcst='000:00:00'

  READ (extdinit,'(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)')                  &
       iyr,imo,iday,ihr,imin,isec

  CALL julday(iyr,imo,iday,jldy)

  myr=MOD(iyr,100)
  ifhr=0
  ifmin=0
  ifsec=0

  READ(extdfcst,'(i3)',ERR=4,END=4) ifhr

  4     CONTINUE

  WRITE (gribtime,'(i4.4,i2.2,i2.2,i2.2,a1,i2.2)')                      &
         iyr,imo,iday,ihr,'f',ifhr

  len1=LEN(dir_extd)
  grbflen=len1

  CALL strlnth( dir_extd, grbflen )

  IF( grbflen == 0 .OR. dir_extd(1:grbflen) == ' ' ) THEN
    dir_extd = '.'
    grbflen=1
  END IF

  IF( dir_extd(grbflen:grbflen) /= '/'.AND.grbflen < len1 ) THEN
    grbflen=grbflen+1
    dir_extd(grbflen:grbflen)='/'
  END IF

  lenrun = LEN( extdname )
  CALL strlnth( extdname, lenrun )

  gribfile = dir_extd(1:grbflen)//extdname(1:lenrun)                    &
                                //'.'//gribtime(1:13)
  grbflen = grbflen + lenrun + 14
!
!-----------------------------------------------------------------------
!
!  RDNMCGRB reads NMC GRIB data
!
!-----------------------------------------------------------------------
!
  gridtyp  = avn3grid
  mproj_grb = avn3proj

  n2dvars  = avn3nvs2d
  n2dlvtps = avn3nlvt2d

  DO k=1,n2dlvtps
    DO n=1,n2dvars
      var_id2d(n,k) = avn3var_id2d(n,k)
    END DO
    levtyp2d(k) = avn3levs2d(k)
  END DO

  n3dvars  = avn3nvs3d
  n3dlvtps = avn3nlvt3d

  DO m=1,n3dlvtps
    DO n=1,n3dvars
      var_id3d(n,m) = avn3var_id3d(n,m)
    END DO
    levtyp3d(m) = avn3levs3d(m)
  END DO

  CALL rdnmcgrb(nx_ext,ny_ext,nz_ext,gribfile,grbflen, gribtime,        &
                gridesc, iproj_grb, gthin,                              &
                ni_grb,nj_grb,np_grb, nk_grb,zk_grb, npeq,nit,          &
                pi_grb,pj_grb,ipole, di_grb,dj_grb,                     &
                latsw,lonsw, latne,lonne,                               &
                latrot,lonrot,angrot,                                   &
                latstr,lonstr,facstr,                                   &
                lattru1,lattru2,lontrue,                                &
                scanmode, iscan,jscan,kscan,                            &
                ires,iearth,icomp,                                      &
                jpenta,kpenta,mpenta,ispect,icoeff,                     &
                xp_grb,yp_grb, xo_grb,yo_grb,zo_grb,                    &
                rcdata,var_grb2d,var_grb3d,var_lev3d,lvldbg,iret)

  IF (iret /= 0)  THEN
    ISTATUS = -888
    GOTO 999
  END IF

  max_nr2d = 0
  DO n=1,n2dvars
    DO m=1,n2dlvtps
      max_nr2d = MAX( max_nr2d, var_nr2d(n,m) )
    END DO
  END DO

  max_nr3d = 0
  min_nr3d = nz_ext
  DO n=1,n3dvars
    max_nr3d = MAX( max_nr3d, var_nr3d(n,1) )
    min_nr3d = MIN( min_nr3d, var_nr3d(n,1) )
  END DO

  IF ( max_nr3d == 0 ) THEN
    WRITE (6,'(a)')                                                     &
        'No 3-D variable was found in the GRIB file',                   &
        'Program stopped in GET_AVN_GRB.'
!   STOP
    ISTATUS = -888
    GOTO 999
  END IF

  IF ( max_nr2d == 0 ) THEN
    WRITE (6,'(a)')                                                     &
        'No 2-D variables was found in the GRIB file'
  END IF

!  write (6,'(/a7,2x,6(i7))')
!    :    'Lev\\VID',(var_id3d(n,1),n=1,n3dvars)
!  DO 60 k=1,max_nr3d
!    var_lev3d(k,5,1) = var_lev3d(k,1,1)
!    var_lev3d(k,6,1) = var_lev3d(k,1,1)
!    write (6,'(/i5,4x,6(i7))')
!    :    k,(var_lev3d(k,n,1),n=1,n3dvars)
  60    CONTINUE

  DO k=1,min_nr3d
    DO n=2,n3dvars
      IF ( var_lev3d(k,1,1) /= var_lev3d(k,n,1) ) THEN
        WRITE (6,'(a)')                                                 &
            'Variables were not at the same level.',                    &
            'Program stopped in GET_AVN_GRB.'
!       STOP
        ISTATUS = -888
        GOTO 999
      END IF
    END DO
  END DO

  IF ( iproj_grb == 5 .AND. ipole == 0 ) THEN      ! Center in NP
    iproj_ext = 1
  ELSE IF ( iproj_grb == 5 .AND. ipole == 1 ) THEN  ! Center in SP
    iproj_ext = -1
  ELSE IF ( iproj_grb == 3 .AND. ipole == 0 ) THEN  ! Center in NP
    iproj_ext = 2
  ELSE IF ( iproj_grb == 3 .AND. ipole == 1 ) THEN  ! Center in SP
    iproj_ext = -2
  ELSE IF ( iproj_grb == 1 ) THEN
    iproj_ext = 3
  ELSE IF ( iproj_grb == 0 ) THEN
    iproj_ext = 4
  ELSE
    WRITE (6,'(a)')                                                     &
        'Unknown map projection. Set to non-projection.'
    iproj_ext = 0
  END IF

  scale_ext = 1.0
  latnot_ext(1) = lattru1
  latnot_ext(2) = lattru2
  trlon_ext = lontrue

  dx_ext = di_grb
  dy_ext = dj_grb

  DO i=ibgn,iend
    DO j=jbgn,jend
      lon_ext(i,j)= lonsw + (i-1) * dx_ext
      lat_ext(i,j)= latsw + (j-1) * dy_ext
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Retrieve 2-D variables
!
!-----------------------------------------------------------------------
!
  DO j=jbgn,jend
    DO i=ibgn,iend

      isrc = MOD(i,360)
      IF( isrc == 0) isrc = 360

      IF ( var_nr2d(1,1) == 0 ) THEN
        psfc_ext   (i,j) = -999.0
      ELSE
        psfc_ext   (i,j) = var_grb2d(isrc,j,1,1)  ! Pa
      END IF

      IF ( var_nr2d(2,1) == 0 ) THEN
        trn_ext    (i,j) = -999.0
      ELSE
        trn_ext    (i,j) = var_grb2d(isrc,j,2,1) ! gpm (same as geometric meter?)
      END IF

      IF ( var_nr2d(1,2) == 0 ) THEN
        tgrd_ext    (i,j) = -999.0
      ELSE
        tgrd_ext    (i,j) = var_grb2d(isrc,j,1,2)
      END IF
      IF ( var_nr2d(2,2) == 0 ) THEN
        ugrd_ext    (i,j) = -999.0
      ELSE
        ugrd_ext    (i,j) = var_grb2d(isrc,j,2,2)
      END IF
      IF ( var_nr2d(3,2) == 0 ) THEN
        vgrd_ext    (i,j) = -999.0
      ELSE
        vgrd_ext    (i,j) = var_grb2d(isrc,j,3,2)
      END IF
      IF ( var_nr2d(4,2) == 0 ) THEN
        rhgrd_ext    (i,j) = -999.0
      ELSE
        rhgrd_ext    (i,j) = var_grb2d(isrc,j,4,2)*0.01 ! 0-1.0
      END IF
      IF ( var_nr2d(5,2) == 0 ) THEN
        ptgrd_ext    (i,j) = -999.0
      ELSE
        ptgrd_ext    (i,j) = var_grb2d(isrc,j,5,2)
      END IF
      IF ( var_nr2d(1,3) == 0 ) THEN
        pmsl_ext     (i,j) = -999.0
      ELSE
        pmsl_ext     (i,j) = var_grb2d(isrc,j,1,3)
      END IF

      IF ( var_nr3d(1,2) == 0 ) THEN
        tsfc_ext  (i,j) = -999.0
        tsoil_ext (i,j) = -999.0
        wetsfc_ext(i,j) = -999.0
        wetdp_ext (i,j) = -999.0
      ELSE
        tsfc_ext  (i,j) = var_grb2d(i,j,3,1)       ! sfc temp.

        IF ( nint(var_grb2d(i,j,5,1)) == 1 ) THEN  ! soil temp over land
          tsoil_ext (i,j) = var_grb3d(isrc,j,1,1,2)

          IF ( tsoil_ext (i,j) <= 200. ) THEN
            tsoil_ext (i,j) = tsfc_ext(i,j)
          END IF

          wetsfc_ext(i,j) = var_grb3d(isrc,j,2,2,2)
          wetdp_ext(i,j)  = var_grb3d(isrc,j,1,2,2)
        ELSE                                       ! sfc temp over sea

          tsoil_ext (i,j) = tsfc_ext(i,j)

          wetsfc_ext(i,j) = 1.0
          wetdp_ext(i,j)  = 1.0
        END IF

      END IF

      IF ( var_nr2d(4,1) == 0 ) THEN
        wetcanp_ext(i,j) = -999.0
      ELSE
        wetcanp_ext(i,j) = var_grb2d(isrc,j,4,1)*1.e-3     ! in meter
      END IF

      IF ( var_nr2d(6,1) == 0 ) THEN
        snowdpth_ext(i,j) = -999.
      ELSE

!      Convert water equiv. of accum. snow depth (kg/m**2) to meters
!      (where 1 meter liquid water is set equivqlent to 10 meters snow).
!          0.01 = 10. (m snow/m liquid) / (1000 kg/m**3)

        snowdpth_ext(i,j) = 0.01 * var_grb2d(isrc,j,6,1)  ! in meters

      END IF

    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Retrieve 3-D variables
!
!-----------------------------------------------------------------------
!
  nz1 = MIN(var_nr3d(1,1),nz_ext)

  IF ( var_lev3d(1,1,1) > var_lev3d(nz1,1,1) ) THEN  ! 1st level at sfc
    chklev = 1
    lvscan = 0
  ELSE
    chklev = -1
    lvscan = nz1+1
  END IF

  DO k=1,nz1
    kk = chklev * k + lvscan
    DO j=jbgn,jend
      DO i=ibgn,iend
        isrc = MOD(i,360)
        IF( isrc == 0) isrc = 360
        p_ext  (i,j,kk) = 100.0 * FLOAT(var_lev3d(k,1,1)) ! Pressure
        hgt_ext(i,j,kk) = var_grb3d(isrc,j,k,1,1)
        u_ext  (i,j,kk) = var_grb3d(isrc,j,k,2,1)    ! u wind (m/s)
        v_ext  (i,j,kk) = var_grb3d(isrc,j,k,3,1)    ! v wind (m/s)

        t_ext  (i,j,kk) = var_grb3d(isrc,j,k,4,1)    ! Temperature (K)

!        qc_ext (i,j,kk) = -999.
!        qr_ext (i,j,kk) = -999.
!        qi_ext (i,j,kk) = -999.
!        qs_ext (i,j,kk) = -999.
!        qh_ext (i,j,kk) = -999.
      END DO
    END DO
  END DO

  CALL getqvs(iend-ibgn+1,jend-jbgn+1,nz1,                              &
              1,iend-ibgn+1,1,jend-jbgn+1,1,nz1,                        &
              p_ext, t_ext, qv_ext )

  nz2 = MIN( nz1, min_nr3d )

  DO k=1,nz2
    kk = chklev * k + lvscan
    DO j=jbgn,jend
      DO i=ibgn,iend
        isrc = MOD(i,360)
        IF( isrc == 0) isrc = 360
        qv_ext(i,j,kk) = 0.01*var_grb3d(isrc,j,k,5,1)*qv_ext(i,j,kk)
      END DO
    END DO
  END DO

  IF ( nz2 < nz1 ) THEN
    DO k=nz2+1,nz1
      kk = chklev * k + lvscan
      DO j=jbgn,jend
        DO i=ibgn,iend
          qv_ext(i,j,kk) = 0.0
          var_lev3d(k,5,1) = var_lev3d(k,1,1)
        END DO
      END DO
    END DO
  END IF

  istatus = 1

  999 CONTINUE
  deallocate(var_grb2d,var_grb3d,rcdata,var_lev3d)

  RETURN
END SUBROUTINE get_avn_grb
