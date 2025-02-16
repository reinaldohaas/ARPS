!
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GETGEMRUC                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE getgemruc(nx_ext, ny_ext, nz_ext,dir_extd,extdname,          &
           extdinit,extdfcst,julfname,i4time,                           &
           iproj_ext,scale_ext,                                         &
           trlon_ext,latnot_ext,x0_ext,y0_ext,                          &
           lat_ext,lon_ext,                                             &
           p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,                      &
           istatus, tem1_ext)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Reads external file for processing by ext2arps, a program
!  that converts external files to ARPS variables and format.
!
!  This version of rdextfil reads RUC (MAPS) data in GEMPAK format
!
!  The script Gemenviron must be run by the process running
!  this program.  It defines GEMPAK symbolic links.
!
!  Because of the GEMPAK parameter include file, implicit none
!  is not used in this program.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  April, 1995
!
!  MODIFICATIONS:
!  Changed external pressure array to be 3-d.to be compatible
!  with new ext2arps.    9 August 1995  Keith Brewster
!
!  Added code to allow for read creating data forecast hours
!  other than those stored in the file by linear time interpolation,
!  when necessary.  Added tem1_ext to the argument list.
!  20 March 1996  Keith Brewster
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    extdinit      Initialized time in mm-dd-yyyy:hh:mm:ss format
!    extdfcst      Forecast hour in HHH:MM:SS format
!    julfname      File name in yyjjjhhmm format
!    i4time        Absolute time in seconds (from 1 Jan 1970) of desired data
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
!    qc_ext        Cloud water mixing ratio (kg/kg)
!    qr_ext        Rain  water mixing ratio (kg/kg)
!    qi_ext        Ice         mixing ratio (kg/kg)
!    qs_ext        Snow        mixing ratio (kg/kg)
!    qh_ext        Hail        mixing ratio (kg/kg)
!    u_ext         u wind component (m/s)
!    v_ext         v wind component (m/s)
!    istatus       status indicator
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INCLUDE 'GEMINC:GEMPRM.PRM'
  INTEGER :: nx_ext, ny_ext, nz_ext
!
  CHARACTER (LEN=* ) :: dir_extd
  CHARACTER (LEN=* ) :: extdname

  CHARACTER (LEN=* ) :: extdinit
  CHARACTER (LEN=* ) :: extdfcst
  CHARACTER (LEN=* ) :: julfname
  INTEGER :: i4time
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
  REAL :: x_ext(nx_ext),y_ext(ny_ext)
!
!------------------------------------------------------------------------
!
!  MAPS/RUC Grid definition parameters
!
!------------------------------------------------------------------------
!
  INTEGER :: iproj_ruc
  REAL :: latnot1_ruc,latnot2_ruc,trlon_ruc,scale_ruc
  PARAMETER (iproj_ruc  =1,  & ! Polar Stereographic
         latnot1_ruc=40.,                                               &
             latnot2_ruc=0.,                                            &
             trlon_ruc  =-105.,                                         &
             scale_ruc  =1.0)
  REAL :: dx_ruc,dy_ruc,swlat_ruc,swlon_ruc
  PARAMETER (dx_ruc = 60000., dy_ruc=60000.,                            &
             swlat_ruc = 22.8373,                                       &
             swlon_ruc = -120.491)
  INTEGER :: fhrinc
  PARAMETER (fhrinc=3)

  INTEGER :: nz_ext_mx
  PARAMETER(nz_ext_mx=41)
!
!------------------------------------------------------------------------
!
!  MAPS/RUC variables
!
!------------------------------------------------------------------------
!
  INTEGER :: ipr_ext(nz_ext_mx)
  DATA ipr_ext / 1100,1075,1050,1025,1000,975,950,925,                  &
                  900, 875, 850, 825, 800,775,750,725,                  &
                  700, 675, 650, 625, 600,575,550,525,                  &
                  500, 475, 450, 425, 400,375,350,325,                  &
                  300, 275, 250, 225, 200,175,150,125,                  &
                  100/
!
!-----------------------------------------------------------------------
!
!  Output external variable arrays
!
!-----------------------------------------------------------------------
!
  REAL :: lat_ext(nx_ext,ny_ext)
  REAL :: lon_ext(nx_ext,ny_ext)
  REAL :: p_ext(nx_ext,ny_ext,nz_ext)
  REAL :: hgt_ext(nx_ext,ny_ext,nz_ext)
  REAL :: t_ext(nx_ext,ny_ext,nz_ext)
  REAL :: qv_ext(nx_ext,ny_ext,nz_ext)
  REAL :: u_ext(nx_ext,ny_ext,nz_ext)
  REAL :: v_ext(nx_ext,ny_ext,nz_ext)
!  REAL :: qc_ext(nx_ext,ny_ext,nz_ext)    ! Cloud H2O mixing ratio (kg/kg)
!  REAL :: qr_ext(nx_ext,ny_ext,nz_ext)    ! Rain  H2O mixing ratio (kg/kg)
!  REAL :: qi_ext(nx_ext,ny_ext,nz_ext)    ! Ice   H2O mixing ratio (kg/kg)
!  REAL :: qs_ext(nx_ext,ny_ext,nz_ext)    ! Snow  H2O mixing ratio (kg/kg)
!  REAL :: qh_ext(nx_ext,ny_ext,nz_ext)    ! Hail  H2O mixing ratio (kg/kg)
  REAL :: tem1_ext(nx_ext,ny_ext,nz_ext)

  INTEGER :: istatus
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
!------------------------------------------------------------------------
!
!  GEMPAK variables
!
!------------------------------------------------------------------------
!
  INTEGER :: inocord,iprcord,ithcord,ihtcord
  PARAMETER (inocord =0, iprcord =1,                                    &
             ithcord =2, ihtcord =3)
!
  REAL :: grid(llnanl)
  REAL :: rnav(llnnav)
  INTEGER :: nxgem,nygem,iret
  INTEGER :: level(2)
  INTEGER :: ighdr(llgdhd)
  CHARACTER (LEN=20)  :: lastim,timep(2),timef(2)
  CHARACTER (LEN=72)  :: gdcur
  CHARACTER (LEN=256) :: gdfile
!
  CHARACTER (LEN=12) :: rucvar
!
!------------------------------------------------------------------------
!
!  Physcial parameters
!
!------------------------------------------------------------------------
!
  REAL :: rd,g
  PARAMETER (rd=287.053, g=9.81)
!
  REAL :: gamma,rddg,xconst
  PARAMETER ( gamma = 0.0065,    & ! degrees/m  lapse rate
          rddg  = (rd/g),                                               &
              xconst = (rd*gamma/g) )
!
!------------------------------------------------------------------------
!
!  Misc internal variables
!
!------------------------------------------------------------------------
!
  CHARACTER (LEN=8) :: gmpktm
  INTEGER :: i,j,k,klev,kstart,itime
  INTEGER :: fhr,fhpast,fhfutr
  INTEGER :: iflno,numgrd,maxgrd,len1,len2
  REAL :: pratio,dln,tbar,const,qvsat,wgtp,wgtf
  INTEGER :: iyr,imo,iday,ihr,imin

  LOGICAL :: tmint,wrtflg,newfil
  INTEGER :: myr
!
  LOGICAL :: init_called

  REAL, ALLOCATABLE :: utmp(:,:), vtmp(:,:)
!
!-----------------------------------------------------------------------
!
!  Function f_qvsat and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
  REAL :: f_qvsatl

!fpp$ expand (f_qvsatl)
!dir$ inline always f_qvsatl

  SAVE init_called
  DATA init_called /.false./
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ALLOCATE(utmp(nx_ext,ny_ext))
  ALLOCATE(vtmp(nx_ext,ny_ext))

  iproj_ext=iproj_ruc
  scale_ext=scale_ruc       ! report lengths in m
  trlon_ext=trlon_ruc      ! orientation of external data grids
  latnot_ext(1)=latnot1_ruc
  latnot_ext(2)=latnot2_ruc
!
  DO k=1,nz_ext
    DO j=1,ny_ext
      DO i=1,nx_ext
        p_ext(i,j,k)=100.*FLOAT(ipr_ext(k))
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Get the lat,lon of the MAPS/RUC grid points
!
!-----------------------------------------------------------------------
!
  CALL getmapr(iproj,scale,latnot,trlon,x0,y0)
  CALL setmapr(iproj_ext,scale_ext,latnot_ext,trlon_ext)
  CALL lltoxy(1,1,swlat_ruc,swlon_ruc,x0_ext,y0_ext)
!
  DO i=1,nx_ext
    x_ext(i)=x0_ext+(i-1)*dx_ruc
  END DO
  DO j=1,ny_ext
    y_ext(j)=y0_ext+(j-1)*dy_ruc
  END DO

  CALL xytoll(nx_ext,ny_ext,x_ext,y_ext,lat_ext,lon_ext)
!
  PRINT *, ' maps point 34,17: ',lat_ext(34,17),lon_ext(34,17)
  PRINT *, ' maps point nx,ny: ',lat_ext(nx_ext,ny_ext),                &
                                 lon_ext(nx_ext,ny_ext)
!
!-----------------------------------------------------------------------
!
!  Initialize GEMPAK sans TAE
!
!-----------------------------------------------------------------------
!
  IF( .NOT.init_called) THEN
    CALL in_bdta(iret)
    init_called=.true.
  END IF
!
!-----------------------------------------------------------------------
!
!  Build RUC file name
!
!-----------------------------------------------------------------------
!
  READ(extdinit,'(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)')                   &
                  iyr,imo,iday,ihr,imin
  myr=MOD(iyr,100)
  len1=LEN(dir_extd)
  len2=len1
  CALL strlnth( dir_extd, len2 )
  IF( dir_extd(len2:len2) /= '/' .AND.  len2 < len1) THEN
    len2=len2+1
    dir_extd(len2:len2)='/'
  END IF

  len1 = LEN( extdname )
  CALL strlnth( extdname, len1 )

  WRITE(gdfile,'(a,a,a,i4.4,i2.2,i2.2,i2.2)')                           &
                 dir_extd(1:len2),extdname(1:len1),'.',                 &
                 iyr,imo,iday,ihr
!
!-----------------------------------------------------------------------
!
!  Open the grid file
!
!-----------------------------------------------------------------------
!
  PRINT *, ' Opening gdfile= ',TRIM(gdfile)
  CALL gr_open  ( gdfile, wrtflg , gdcur, iflno, lastim,                &
                  grid, rnav, numgrd, maxgrd, newfil, iret)
  IF  ( iret == 0 )  THEN
    PRINT *, ' gdcur = ',gdcur
    PRINT *, ' iflno = ',iflno
    PRINT *, ' lastim = ',lastim
    PRINT *, ' numgrd = ',numgrd
    PRINT *, ' maxgrd = ',maxgrd
    PRINT *, ' newfil = ',newfil
    PRINT *, ' iret = ',iret
  ELSE
    WRITE(6,'(a,a,/,a,i4)') '  Error opening file ',gdfile,             &
                 '  GEMPAK GR_OPEN return status: ',iret
    istatus=iret
    RETURN
  END IF
!
!-----------------------------------------------------------------------
!
!  Build the  GEMPAK grid time string
!  It has format yymodd/hhmnFHHH
!  yy: year      mo: month  dd: GMT day
!  hh: GMT hour  mn: minute
!  F: separation charcter
!  HHH: forecast hour (000 = analysis)
!  example  time(1)='950126/1200F000'
!
!-----------------------------------------------------------------------
!
  timep(1)='                    '
  timep(2)='                    '
  timef(1)='                    '
  timef(2)='                    '
  READ(extdfcst,'(i3)') fhr
  IF(MOD(fhr,fhrinc) /= 0) THEN
    tmint=.true.
    fhpast=(fhr/fhrinc)*fhrinc
    fhfutr=fhpast+fhrinc
    wgtp=FLOAT(fhfutr-fhr)/FLOAT(fhrinc)
    wgtf=1.-wgtp
  ELSE
    tmint=.false.
    fhpast=fhr
    fhfutr=fhr
    wgtp=1.0
    wgtf=0.0
  END IF

  WRITE(timep(1),'(i2.2,i2.2,i2.2,a1,i2.2,i2.2,a1,i3.3)')               &
        myr,imo,iday,'/',ihr,imin,'F',fhpast
  WRITE(6,'(a,a)') '  Past   GEMPAK time string ',timep(1)
  WRITE(timef(1),'(i2.2,i2.2,i2.2,a1,i2.2,i2.2,a1,i3.3)')               &
        myr,imo,iday,'/',ihr,imin,'F',fhfutr
  WRITE(6,'(a,a)') '  Future GEMPAK time string ',timef(1)
!
!-----------------------------------------------------------------------
!
!  Data in the RUC files are only available from 100-1000 mb.
!  Thus this process does not start at pr_ext vertical level one,
!  which is at 1100 mb.  Find 1000 mb in the ipr_ext vector.
!
!-----------------------------------------------------------------------
!
  DO klev=1,nz_ext-1
    IF(ipr_ext(klev) <= 1000) EXIT
  END DO
  kstart=klev
  PRINT *, ' kstart = ',kstart
!
!-----------------------------------------------------------------------
!
!  Go through each of the RUC variables collecting those
!  interpolated to pressure surfaces at NMC from the
!  original model coordinates (hybrid sigma-isentropic)
!
!-----------------------------------------------------------------------
!
  DO klev=kstart,nz_ext
    PRINT *, ' Reading level ',ipr_ext(klev)
    level(1)=ipr_ext(klev)
    level(2)=-1
!
!-----------------------------------------------------------------------
!
!  Heights
!
!-----------------------------------------------------------------------
!
    rucvar='HGHT'
    CALL gd_rdat  ( iflno, timep,                                       &
                    level, iprcord, rucvar,                             &
                    hgt_ext(1,1,klev), nxgem, nygem, ighdr, iret )
    IF  ( iret == 0 )  THEN
      PRINT *, ' tpast hght(40,40) = ',hgt_ext(40,40,klev)
      IF( nxgem /= nx_ext .OR. nygem /= ny_ext ) THEN
        WRITE(6,'(a,/,a,2i12,/,a,2i12)')                                &
            '  Error in grid dimensions.',                              &
            '  GEMPAK returned nx,ny:',nxgem,nygem,                     &
            '  Expected        nx,ny:',nx_ext,ny_ext
        istatus=-31
        RETURN
      END IF
      IF(tmint) THEN
        CALL gd_rdat  ( iflno, timef,                                   &
                        level, iprcord, rucvar,                         &
                        tem1_ext, nxgem, nygem, ighdr, iret )
        IF  ( iret == 0 )  THEN
          PRINT *, ' tfutr hght(40,40) = ',tem1_ext(40,40,1)
          DO j=1,ny_ext
            DO i=1,nx_ext
              hgt_ext(i,j,klev)=wgtp*hgt_ext(i,j,klev)+                 &
                                wgtf*tem1_ext(i,j,1)
            END DO
          END DO
        END IF
      END IF
    ELSE
      WRITE(6,'(a,a,/,a,i4)') '  Error reading file ',gdfile,           &
              '  GEMPAK GD_RDAT return status: ',iret
      istatus=iret
      RETURN
    END IF
!
!-----------------------------------------------------------------------
!
!  Temperatures
!
!-----------------------------------------------------------------------
!
    rucvar='TMPK'
    CALL gd_rdat  ( iflno, timep,                                       &
                    level, iprcord, rucvar,                             &
                    t_ext(1,1,klev), nxgem, nygem, ighdr, iret )
    IF  ( iret == 0 )  THEN
      PRINT *, ' tpast t(40,40) = ',t_ext(40,40,klev)
      IF( nxgem /= nx_ext .OR. nygem /= ny_ext ) THEN
        WRITE(6,'(a,/,a,2i12,/,a,2i12)')                                &
            '  Error in grid dimensions.',                              &
            '  GEMPAK returned nx,ny:',nxgem,nygem,                     &
            '  Expected        nx,ny:',nx_ext,ny_ext
        istatus=-31
        RETURN
      END IF
      IF(tmint) THEN
        CALL gd_rdat  ( iflno, timef,                                   &
                        level, iprcord, rucvar,                         &
                        tem1_ext, nxgem, nygem, ighdr, iret )
        IF  ( iret == 0 )  THEN
          PRINT *, ' tfutr t(40,40) = ',tem1_ext(40,40,1)
          DO j=1,ny_ext
            DO i=1,nx_ext
              t_ext(i,j,klev)=wgtp*t_ext(i,j,klev)+                     &
                              wgtf*tem1_ext(i,j,1)
            END DO
          END DO
        END IF
      END IF
    ELSE
      WRITE(6,'(a,a,/,a,i4)') '  Error reading file ',gdfile,           &
              '  GEMPAK GD_RDAT return status: ',iret
      istatus=iret
      RETURN
    END IF
!
!-----------------------------------------------------------------------
!
!  Relative humidity
!
!-----------------------------------------------------------------------
!
    rucvar='RELH'
    CALL gd_rdat  ( iflno, timep,                                       &
                    level, iprcord, rucvar,                             &
                    qv_ext(1,1,klev), nxgem, nygem, ighdr, iret )
    IF  ( iret == 0 )  THEN
      PRINT *, ' tpast rh(40,40) = ',qv_ext(40,40,klev)
      IF( nxgem /= nx_ext .OR. nygem /= ny_ext ) THEN
        WRITE(6,'(a,/,a,2i12,/,a,2i12)')                                &
            '  Error in grid dimensions.',                              &
            '  GEMPAK returned nx,ny:',nxgem,nygem,                     &
            '  Expected        nx,ny:',nx_ext,ny_ext
        istatus=-31
        RETURN
      END IF
      IF(tmint) THEN
        CALL gd_rdat  ( iflno, timef,                                   &
                        level, iprcord, rucvar,                         &
                        tem1_ext, nxgem, nygem, ighdr, iret )
        IF  ( iret == 0 )  THEN
          PRINT *, ' tfutr rh(40,40) = ',tem1_ext(40,40,1)
          DO j=1,ny_ext
            DO i=1,nx_ext
              qv_ext(i,j,klev)=wgtp*qv_ext(i,j,klev)+                   &
                               wgtf*tem1_ext(i,j,1)
            END DO
          END DO
        END IF
      END IF
    ELSE
      WRITE(6,'(a,a,/,a,i4)') '  Error reading file ',gdfile,           &
              '  GEMPAK GD_RDAT return status: ',iret
      istatus=iret
      RETURN
    END IF
!
!-----------------------------------------------------------------------
!
!  U - Grid relative velocities
!
!-----------------------------------------------------------------------
!
    rucvar='UREL'
    CALL gd_rdat  ( iflno, timep,                                       &
                    level, iprcord, rucvar,                             &
                    u_ext(1,1,klev), nxgem, nygem, ighdr, iret )
    IF  ( iret == 0 )  THEN
      PRINT *, ' tpast u(40,40) = ',u_ext(40,40,klev)
      IF( nxgem /= nx_ext .OR. nygem /= ny_ext ) THEN
        WRITE(6,'(a,/,a,2i12,/,a,2i12)')                                &
            '  Error in grid dimensions.',                              &
            '  GEMPAK returned nx,ny:',nxgem,nygem,                     &
            '  Expected        nx,ny:',nx_ext,ny_ext
        istatus=-31
        RETURN
      END IF
      IF(tmint) THEN
        CALL gd_rdat  ( iflno, timef,                                   &
                        level, iprcord, rucvar,                         &
                        tem1_ext, nxgem, nygem, ighdr, iret )
        IF  ( iret == 0 )  THEN
          PRINT *, ' tfutr u(40,40) = ',tem1_ext(40,40,1)
          DO j=1,ny_ext
            DO i=1,nx_ext
              u_ext(i,j,klev)=wgtp*u_ext(i,j,klev)+                     &
                              wgtf*tem1_ext(i,j,1)
            END DO
          END DO
        END IF
      END IF
    ELSE
      WRITE(6,'(a,a,/,a,i4)') '  Error reading file ',gdfile,           &
              '  GEMPAK GD_RDAT return status: ',iret
      istatus=iret
      RETURN
    END IF
!
!-----------------------------------------------------------------------
!
!  V - Grid relative velocities
!
!-----------------------------------------------------------------------
!
    rucvar='VREL'
    CALL gd_rdat  ( iflno, timep,                                       &
                    level, iprcord, rucvar,                             &
                    v_ext(1,1,klev), nxgem, nygem, ighdr, iret )
    IF  ( iret == 0 )  THEN
      PRINT *, ' tpast v(40,40) = ',v_ext(40,40,klev)
      IF( nxgem /= nx_ext .OR. nygem /= ny_ext ) THEN
        WRITE(6,'(a,/,a,2i12,/,a,2i12)')                                &
            '  Error in grid dimensions.',                              &
            '  GEMPAK returned nx,ny:',nxgem,nygem,                     &
            '  Expected        nx,ny:',nx_ext,ny_ext
        istatus=-31
        RETURN
      END IF
      IF(tmint) THEN
        CALL gd_rdat  ( iflno, timef,                                   &
                        level, iprcord, rucvar,                         &
                        tem1_ext, nxgem, nygem, ighdr, iret )
        IF  ( iret == 0 )  THEN
          PRINT *, ' tfutr v(40,40) = ',tem1_ext(40,40,1)
          DO j=1,ny_ext
            DO i=1,nx_ext
              v_ext(i,j,klev)=wgtp*v_ext(i,j,klev)+                     &
                              wgtf*tem1_ext(i,j,1)
            END DO
          END DO
        END IF
      END IF
    ELSE
      WRITE(6,'(a,a,/,a,i4)') '  Error reading file ',gdfile,           &
              '  GEMPAK GD_RDAT return status: ',iret
      istatus=iret
      RETURN
    END IF
  END DO
!
!-----------------------------------------------------------------------
!
!  Extrapolate data to fill in data below 1000 mb
!  To begin, set all to be equal to values at lowest
!  available MAPS/RUC level.
!
!-----------------------------------------------------------------------
!
  DO klev=1,kstart-1
    DO j=1,ny_ext
      DO i=1,nx_ext
        hgt_ext(i,j,klev)=hgt_ext(i,j,kstart)
        t_ext(i,j,klev)  =t_ext(i,j,kstart)
        qv_ext(i,j,klev) =qv_ext(i,j,kstart)
        u_ext(i,j,klev)  =u_ext(i,j,kstart)
        v_ext(i,j,klev)  =v_ext(i,j,kstart)
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Extrapolate temperature using standard atmospheric
!  lapse rate.
!
!-----------------------------------------------------------------------
!
  DO klev=1,kstart-1
    DO j=1,ny_ext
      DO i=1,nx_ext
        t_ext(i,j,klev)=t_ext(i,j,kstart)*                              &
               ( p_ext(i,j,klev) /                                      &
                 p_ext(i,j,kstart) )**xconst
      END DO
    END DO
    PRINT *, ' pr,tks,text = ',ipr_ext(klev),                           &
               t_ext(60,36,kstart),t_ext(60,36,klev)
  END DO
!
!-----------------------------------------------------------------------
!
!  Change RUC Relative Humidity to Mixing Ratio for ARPS.
!  RH is 0-100% so it is multipled by 0.01
!
!-----------------------------------------------------------------------
!
  DO klev=1,nz_ext
    DO j=1,ny_ext
      DO i=1,nx_ext
        qvsat = f_qvsatl( p_ext(i,j,klev), t_ext(i,j,klev) )
        qv_ext(i,j,klev)=0.01*qv_ext(i,j,klev)*qvsat
      END DO
    END DO
  END DO
  PRINT *, ipr_ext(2),' mixing ratio ',qv_ext(60,36,2)
  PRINT *, ipr_ext(4),' mixing ratio ',qv_ext(60,36,4)
  PRINT *, ipr_ext(6),' mixing ratio ',qv_ext(60,36,6)
  PRINT *, ipr_ext(8),' mixing ratio ',qv_ext(60,36,8)
!
!-----------------------------------------------------------------------
!
!  Set height field
!  by integrating T down from kstart level
!
!-----------------------------------------------------------------------
!
  DO klev=kstart-1,1,-1
    DO j=1,ny_ext
      DO i=1,nx_ext
        dln=ALOG(p_ext(i,j,klev)/p_ext(i,j,klev+1))
        const=dln*rddg
        tbar=0.5*(t_ext(i,j,klev)+t_ext(i,j,klev+1))
        hgt_ext(i,j,klev)=hgt_ext(i,j,klev+1)-const*tbar
      END DO
    END DO
    PRINT *, ' pr,hgt(ks),hgt = ',ipr_ext(klev),                        &
               hgt_ext(60,36,kstart),hgt_ext(60,36,klev)
  END DO
!
!-----------------------------------------------------------------------
!
!  Rotate winds to be relative to true north.
!  The RUC data are sent as grid-relative.
!
!-----------------------------------------------------------------------
!
  DO klev=1,nz_ext
!2001-05-16 GMB: Having umap & uear (or vmap & vear) point to
!the same array causes numerical errors when optimizing.
    CALL uvmptoe(nx_ext,ny_ext,u_ext(1,1,klev),v_ext(1,1,klev),         &
                 lon_ext,utmp,vtmp)
    u_ext(:,:,klev) = utmp(:,:)
    v_ext(:,:,klev) = vtmp(:,:)
  END DO
!
!-----------------------------------------------------------------------
!
!  Fill qc,qr,qi,qs,qh arrays with missing value.
!
!-----------------------------------------------------------------------
!
!  DO k=1,nz_ext
!    DO j=1,ny_ext
!      DO i=1,nx_ext
!        qc_ext(i,j,k)=-999.
!        qr_ext(i,j,k)=-999.
!        qi_ext(i,j,k)=-999.
!        qs_ext(i,j,k)=-999.
!        qh_ext(i,j,k)=-999.
!      END DO
!    END DO
!  END DO
!
!-----------------------------------------------------------------------
!
!  Reset map projection to previous values
!
!-----------------------------------------------------------------------
!
  CALL setmapr(iproj,scale,latnot,trlon)
  CALL setorig(1,x0,y0)
!
!-----------------------------------------------------------------------
!
!  Set good status
!
!-----------------------------------------------------------------------
!
  istatus=1

  DEALLOCATE(utmp)
  DEALLOCATE(vtmp)

  RETURN
END SUBROUTINE getgemruc
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GETGEMRUC2                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE getgemruc2(nx_ext, ny_ext, nz_ext,dir_extd,extdname,         &
           extdinit,extdfcst,julfname,i4time,                           &
           iproj_ext,scale_ext,                                         &
           trlon_ext,latnot_ext,x0_ext,y0_ext,                          &
           lat_ext,lon_ext,                                             &
           p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,                      &
           istatus, tem1_ext)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Reads external file for processing by ext2arps, a program
!  that converts external files to ARPS variables and format.
!
!  This version of rdextfil reads RUC (MAPS) data in GEMPAK format
!
!  The script Gemenviron must be run by the process running
!  this program.  It defines GEMPAK symbolic links.
!
!  Because of the GEMPAK parameter include file, implicit none
!  is not used in this program.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  April, 1995
!
!  MODIFICATIONS:
!  Changed external pressure array to be 3-d.to be compatible
!  with new ext2arps.    9 August 1995  Keith Brewster
!
!  Added code to allow for read creating data forecast hours
!  other than those stored in the file by linear time interpolation,
!  when necessary.  Added tem1_ext to the argument list.
!  20 March 1996  Keith Brewster
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    extdinit      Initialized time in mm-dd-yyyy:hh:mm:ss format
!    extdfcst      Forecast hour in HHH:MM:SS format
!    julfname      File name in yyjjjhhmm format
!    i4time        Absolute time in seconds (from 1 Jan 1970) of desired data
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
!    qc_ext        Cloud water mixing ratio (kg/kg)
!    qr_ext        Rain  water mixing ratio (kg/kg)
!    qi_ext        Ice         mixing ratio (kg/kg)
!    qs_ext        Snow        mixing ratio (kg/kg)
!    qh_ext        Hail        mixing ratio (kg/kg)
!    u_ext         u wind component (m/s)
!    v_ext         v wind component (m/s)
!    istatus       status indicator
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INCLUDE 'GEMINC:GEMPRM.PRM'
  INTEGER :: nx_ext, ny_ext, nz_ext
!
  CHARACTER (LEN=* ) :: dir_extd
  CHARACTER (LEN=* ) :: extdname

  CHARACTER (LEN=* ) :: extdinit
  CHARACTER (LEN=* ) :: extdfcst
  CHARACTER (LEN=* ) :: julfname
  INTEGER :: i4time
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
  REAL :: x_ext(nx_ext),y_ext(ny_ext)
!
!------------------------------------------------------------------------
!
!  MAPS/RUC Grid definition parameters
!
!------------------------------------------------------------------------
!
  INTEGER :: iproj_ruc
  REAL :: latnot1_ruc,latnot2_ruc,trlon_ruc,scale_ruc
  PARAMETER (iproj_ruc  =2,  & ! Lambert Conformal
         latnot1_ruc=25.0,                                              &
             latnot2_ruc=25.0,                                          &
             trlon_ruc  =-95.0,                                         &
             scale_ruc  =1.0)
  REAL :: dx_ruc,dy_ruc,swlat_ruc,swlon_ruc
  PARAMETER (dx_ruc = 40.63525E3, dy_ruc=dx_ruc,                        &
             swlat_ruc = 16.2810,                                       &
             swlon_ruc = -126.1378)
  INTEGER :: fhrinc
  PARAMETER (fhrinc=3)      ! 1998/02/03 file contains hours 0-3,6,9,12.
!
!------------------------------------------------------------------------
!
!  MAPS/RUC variables
!
!------------------------------------------------------------------------
!
  INTEGER :: nz_ext_mx
  PARAMETER(nz_ext_mx=41)
  INTEGER :: ipr_ext(nz_ext_mx)
  DATA ipr_ext / 1100,1075,1050,1025,1000,975,950,925,                  &
                  900, 875, 850, 825, 800,775,750,725,                  &
                  700, 675, 650, 625, 600,575,550,525,                  &
                  500, 475, 450, 425, 400,375,350,325,                  &
                  300, 275, 250, 225, 200,175,150,125,                  &
                  100/
!
!-----------------------------------------------------------------------
!
!  Output external variable arrays
!
!-----------------------------------------------------------------------
!
  REAL :: lat_ext(nx_ext,ny_ext)
  REAL :: lon_ext(nx_ext,ny_ext)
  REAL :: p_ext(nx_ext,ny_ext,nz_ext)
  REAL :: hgt_ext(nx_ext,ny_ext,nz_ext)
  REAL :: t_ext(nx_ext,ny_ext,nz_ext)
  REAL :: qv_ext(nx_ext,ny_ext,nz_ext)
  REAL :: u_ext(nx_ext,ny_ext,nz_ext)
  REAL :: v_ext(nx_ext,ny_ext,nz_ext)
!  REAL :: qc_ext(nx_ext,ny_ext,nz_ext)    ! Cloud H2O mixing ratio (kg/kg)
!  REAL :: qr_ext(nx_ext,ny_ext,nz_ext)    ! Rain  H2O mixing ratio (kg/kg)
!  REAL :: qi_ext(nx_ext,ny_ext,nz_ext)    ! Ice   H2O mixing ratio (kg/kg)
!  REAL :: qs_ext(nx_ext,ny_ext,nz_ext)    ! Snow  H2O mixing ratio (kg/kg)
!  REAL :: qh_ext(nx_ext,ny_ext,nz_ext)    ! Hail  H2O mixing ratio (kg/kg)
  REAL :: tem1_ext(nx_ext,ny_ext,nz_ext)

  INTEGER :: istatus
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
!------------------------------------------------------------------------
!
!  GEMPAK variables
!
!------------------------------------------------------------------------
!
  INTEGER :: inocord,iprcord,ithcord,ihtcord
  PARAMETER (inocord =0, iprcord =1,                                    &
             ithcord =2, ihtcord =3)
!
  REAL :: grid(llnanl)
  REAL :: rnav(llnnav)
  INTEGER :: nxgem,nygem,iret
  INTEGER :: level(2)
  INTEGER :: ighdr(llgdhd)
  CHARACTER (LEN=20)  :: lastim,timep(2),timef(2)
  CHARACTER (LEN=72)  :: gdcur
  CHARACTER (LEN=256) :: gdfile
!
  CHARACTER (LEN=12) :: rucvar
!
!------------------------------------------------------------------------
!
!  Physcial parameters
!
!------------------------------------------------------------------------
!
  REAL :: rd,g
  PARAMETER (rd=287.053, g=9.81)
!
  REAL :: gamma,rddg,xconst
  PARAMETER ( gamma = 0.0065,    & ! degrees/m  lapse rate
          rddg  = (rd/g),                                               &
              xconst = (rd*gamma/g) )
!
!------------------------------------------------------------------------
!
!  Misc internal variables
!
!------------------------------------------------------------------------
!
  CHARACTER (LEN=8) :: gmpktm
  INTEGER :: i,j,k,klev,kstart,itime
  INTEGER :: fhr,fhpast,fhfutr
  INTEGER :: iflno,numgrd,maxgrd,len1,len2
  REAL :: pratio,dln,tbar,const,qvsat,wgtp,wgtf
  INTEGER :: iyr,imo,iday,ihr,imin
  LOGICAL :: tmint,wrtflg,newfil
  INTEGER :: myr

  REAL, ALLOCATABLE :: utmp(:,:), vtmp(:,:)
!
!-----------------------------------------------------------------------
!
!  Function f_qvsatl and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
  REAL :: f_qvsatl

!fpp$ expand (f_qvsatl)
!dir$ inline always f_qvsatl

  LOGICAL :: init_called
  SAVE init_called
  DATA init_called /.false./
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ALLOCATE(utmp(nx_ext,ny_ext))
  ALLOCATE(vtmp(nx_ext,ny_ext))

  iproj_ext=iproj_ruc
  scale_ext=scale_ruc       ! report lengths in m
  trlon_ext=trlon_ruc      ! orientation of external data grids
  latnot_ext(1)=latnot1_ruc
  latnot_ext(2)=latnot2_ruc
!
  DO k=1,nz_ext
    DO j=1,ny_ext
      DO i=1,nx_ext
        p_ext(i,j,k)=100.*FLOAT(ipr_ext(k))
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Get the lat,lon of the MAPS/RUC grid points
!
!-----------------------------------------------------------------------
!
  CALL getmapr(iproj,scale,latnot,trlon,x0,y0)
  CALL setmapr(iproj_ext,scale_ext,latnot_ext,trlon_ext)
  CALL lltoxy(1,1,swlat_ruc,swlon_ruc,x0_ext,y0_ext)
!
  DO i=1,nx_ext
    x_ext(i)=x0_ext+(i-1)*dx_ruc
  END DO
  DO j=1,ny_ext
    y_ext(j)=y0_ext+(j-1)*dy_ruc
  END DO

  CALL xytoll(nx_ext,ny_ext,x_ext,y_ext,lat_ext,lon_ext)
!
  PRINT *, ' maps point 34,17: ',lat_ext(34,17),lon_ext(34,17)
  PRINT *, ' maps point nx,ny: ',lat_ext(nx_ext,ny_ext),                &
                                 lon_ext(nx_ext,ny_ext)
!
!-----------------------------------------------------------------------
!
!  Initialize GEMPAK sans TAE
!
!-----------------------------------------------------------------------
!
  IF( .NOT.init_called) THEN
    CALL in_bdta(iret)
    init_called=.true.
  END IF
!
!-----------------------------------------------------------------------
!
!  Build RUC file name
!
!-----------------------------------------------------------------------
!
  READ(extdinit,'(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)')                   &
                  iyr,imo,iday,ihr,imin
  myr=MOD(iyr,100)
  len1=LEN(dir_extd)
  len2=len1
  CALL strlnth( dir_extd, len2 )
  IF( dir_extd(len2:len2) /= '/' .AND.  len2 < len1) THEN
    len2=len2+1
    dir_extd(len2:len2)='/'
  END IF

  len1 = LEN( extdname )
  CALL strlnth( extdname, len1 )

  WRITE(gdfile,'(a,a,a,i4.4,i2.2,i2.2,i2.2)')                           &
                 dir_extd(1:len2),extdname(1:len1),'.',                 &
                 iyr,imo,iday,ihr
!
!-----------------------------------------------------------------------
!
!  Open the grid file
!
!-----------------------------------------------------------------------
!
  PRINT *, ' Opening gdfile= ',gdfile
  CALL gr_open  ( gdfile, wrtflg , gdcur, iflno, lastim,                &
                  grid, rnav, numgrd, maxgrd, newfil, iret)
  IF  ( iret == 0 )  THEN
    PRINT *, ' gdcur = ',gdcur
    PRINT *, ' iflno = ',iflno
    PRINT *, ' lastim = ',lastim
    PRINT *, ' numgrd = ',numgrd
    PRINT *, ' maxgrd = ',maxgrd
    PRINT *, ' newfil = ',newfil
    PRINT *, ' iret = ',iret
  ELSE
    WRITE(6,'(a,a,/,a,i4)') '  Error opening file ',gdfile,             &
                 '  GEMPAK GR_OPEN return status: ',iret
    istatus=iret
    RETURN
  END IF
!
!-----------------------------------------------------------------------
!
!  Build the  GEMPAK grid time string
!  It has format yymodd/hhmnFHHH
!  yy: year      mo: month  dd: GMT day
!  hh: GMT hour  mn: minute
!  F: separation charcter
!  HHH: forecast hour (000 = analysis)
!  example  time(1)='950126/1200F000'
!
!-----------------------------------------------------------------------
!
  timep(1)='                    '
  timep(2)='                    '
  timef(1)='                    '
  timef(2)='                    '
  READ(extdfcst,'(i3)') fhr
  IF(MOD(fhr,fhrinc) /= 0) THEN
    tmint=.true.
    fhpast=(fhr/fhrinc)*fhrinc
    fhfutr=fhpast+fhrinc
    wgtp=FLOAT(fhfutr-fhr)/FLOAT(fhrinc)
    wgtf=1.-wgtp
  ELSE
    tmint=.false.
    fhpast=fhr
    fhfutr=fhr
    wgtp=1.0
    wgtf=0.0
  END IF

  WRITE(timep(1),'(i2.2,i2.2,i2.2,a1,i2.2,i2.2,a1,i3.3)')               &
        myr,imo,iday,'/',ihr,imin,'F',fhpast
  WRITE(6,'(a,a)') '  Past   GEMPAK time string ',timep(1)
  WRITE(timef(1),'(i2.2,i2.2,i2.2,a1,i2.2,i2.2,a1,i3.3)')               &
        myr,imo,iday,'/',ihr,imin,'F',fhfutr
  WRITE(6,'(a,a)') '  Future GEMPAK time string ',timef(1)
!
!-----------------------------------------------------------------------
!
!  Data in the RUC files are only available from 100-1000 mb.
!  Thus this process does not start at pr_ext vertical level one,
!  which is at 1100 mb.  Find 1000 mb in the ipr_ext vector.
!
!-----------------------------------------------------------------------
!
  DO klev=1,nz_ext-1
    IF(ipr_ext(klev) <= 1000) EXIT
  END DO
  kstart=klev
  PRINT *, ' kstart = ',kstart
!
!-----------------------------------------------------------------------
!
!  Go through each of the RUC variables collecting those
!  interpolated to pressure surfaces at NMC from the
!  original model coordinates (hybrid sigma-isentropic)
!
!-----------------------------------------------------------------------
!
  DO klev=kstart,nz_ext
    PRINT *, ' Reading level ',ipr_ext(klev)
    level(1)=ipr_ext(klev)
    level(2)=-1
!
!-----------------------------------------------------------------------
!
!  Heights
!
!-----------------------------------------------------------------------
!
    rucvar='HGHT'
    CALL gd_rdat  ( iflno, timep,                                       &
                    level, iprcord, rucvar,                             &
                    hgt_ext(1,1,klev), nxgem, nygem, ighdr, iret )
    IF  ( iret == 0 )  THEN
      PRINT *, ' tpast hght(40,40) = ',hgt_ext(40,40,klev)
      IF( nxgem /= nx_ext .OR. nygem /= ny_ext ) THEN
        WRITE(6,'(a,/,a,2i12,/,a,2i12)')                                &
            '  Error in grid dimensions.',                              &
            '  GEMPAK returned nx,ny:',nxgem,nygem,                     &
            '  Expected        nx,ny:',nx_ext,ny_ext
        istatus=-31
        RETURN
      END IF
      IF(tmint) THEN
        CALL gd_rdat  ( iflno, timef,                                   &
                        level, iprcord, rucvar,                         &
                        tem1_ext, nxgem, nygem, ighdr, iret )
        IF  ( iret == 0 )  THEN
          PRINT *, ' tfutr hght(40,40) = ',tem1_ext(40,40,1)
          DO j=1,ny_ext
            DO i=1,nx_ext
              hgt_ext(i,j,klev)=wgtp*hgt_ext(i,j,klev)+                 &
                                wgtf*tem1_ext(i,j,1)
            END DO
          END DO
        END IF
      END IF
    ELSE
      WRITE(6,'(a,a,/,a,i4)') '  Error reading file ',gdfile,           &
              '  GEMPAK GD_RDAT return status: ',iret
      istatus=iret
      RETURN
    END IF
!
!-----------------------------------------------------------------------
!
!  Temperatures
!
!-----------------------------------------------------------------------
!
    rucvar='TMPK'
    CALL gd_rdat  ( iflno, timep,                                       &
                    level, iprcord, rucvar,                             &
                    t_ext(1,1,klev), nxgem, nygem, ighdr, iret )
    IF  ( iret == 0 )  THEN
      PRINT *, ' tpast t(40,40) = ',t_ext(40,40,klev)
      IF( nxgem /= nx_ext .OR. nygem /= ny_ext ) THEN
        WRITE(6,'(a,/,a,2i12,/,a,2i12)')                                &
            '  Error in grid dimensions.',                              &
            '  GEMPAK returned nx,ny:',nxgem,nygem,                     &
            '  Expected        nx,ny:',nx_ext,ny_ext
        istatus=-31
        RETURN
      END IF
      IF(tmint) THEN
        CALL gd_rdat  ( iflno, timef,                                   &
                        level, iprcord, rucvar,                         &
                        tem1_ext, nxgem, nygem, ighdr, iret )
        IF  ( iret == 0 )  THEN
          PRINT *, ' tfutr t(40,40) = ',tem1_ext(40,40,1)
          DO j=1,ny_ext
            DO i=1,nx_ext
              t_ext(i,j,klev)=wgtp*t_ext(i,j,klev)+                     &
                              wgtf*tem1_ext(i,j,1)
            END DO
          END DO
        END IF
      END IF
    ELSE
      WRITE(6,'(a,a,/,a,i4)') '  Error reading file ',gdfile,           &
              '  GEMPAK GD_RDAT return status: ',iret
      istatus=iret
      RETURN
    END IF
!
!-----------------------------------------------------------------------
!
!  Relative humidity
!
!-----------------------------------------------------------------------
!
    rucvar='RELH'
    CALL gd_rdat  ( iflno, timep,                                       &
                    level, iprcord, rucvar,                             &
                    qv_ext(1,1,klev), nxgem, nygem, ighdr, iret )
    IF  ( iret == 0 )  THEN
      PRINT *, ' tpast rh(40,40) = ',qv_ext(40,40,klev)
      IF( nxgem /= nx_ext .OR. nygem /= ny_ext ) THEN
        WRITE(6,'(a,/,a,2i12,/,a,2i12)')                                &
            '  Error in grid dimensions.',                              &
            '  GEMPAK returned nx,ny:',nxgem,nygem,                     &
            '  Expected        nx,ny:',nx_ext,ny_ext
        istatus=-31
        RETURN
      END IF
      IF(tmint) THEN
        CALL gd_rdat  ( iflno, timef,                                   &
                        level, iprcord, rucvar,                         &
                        tem1_ext, nxgem, nygem, ighdr, iret )
        IF  ( iret == 0 )  THEN
          PRINT *, ' tfutr rh(40,40) = ',tem1_ext(40,40,1)
          DO j=1,ny_ext
            DO i=1,nx_ext
              qv_ext(i,j,klev)=wgtp*qv_ext(i,j,klev)+                   &
                               wgtf*tem1_ext(i,j,1)
            END DO
          END DO
        END IF
      END IF
    ELSE
      WRITE(6,'(a,a,/,a,i4)') '  Error reading file ',gdfile,           &
              '  GEMPAK GD_RDAT return status: ',iret
      istatus=iret
      RETURN
    END IF
!
!-----------------------------------------------------------------------
!
!  U - Grid relative velocities
!
!-----------------------------------------------------------------------
!
    rucvar='UREL'
    CALL gd_rdat  ( iflno, timep,                                       &
                    level, iprcord, rucvar,                             &
                    u_ext(1,1,klev), nxgem, nygem, ighdr, iret )
    IF  ( iret == 0 )  THEN
      PRINT *, ' tpast u(40,40) = ',u_ext(40,40,klev)
      IF( nxgem /= nx_ext .OR. nygem /= ny_ext ) THEN
        WRITE(6,'(a,/,a,2i12,/,a,2i12)')                                &
            '  Error in grid dimensions.',                              &
            '  GEMPAK returned nx,ny:',nxgem,nygem,                     &
            '  Expected        nx,ny:',nx_ext,ny_ext
        istatus=-31
        RETURN
      END IF
      IF(tmint) THEN
        CALL gd_rdat  ( iflno, timef,                                   &
                        level, iprcord, rucvar,                         &
                        tem1_ext, nxgem, nygem, ighdr, iret )
        IF  ( iret == 0 )  THEN
          PRINT *, ' tfutr u(40,40) = ',tem1_ext(40,40,1)
          DO j=1,ny_ext
            DO i=1,nx_ext
              u_ext(i,j,klev)=wgtp*u_ext(i,j,klev)+                     &
                              wgtf*tem1_ext(i,j,1)
            END DO
          END DO
        END IF
      END IF
    ELSE
      WRITE(6,'(a,a,/,a,i4)') '  Error reading file ',gdfile,           &
              '  GEMPAK GD_RDAT return status: ',iret
      istatus=iret
      RETURN
    END IF
!
!-----------------------------------------------------------------------
!
!  V - Grid relative velocities
!
!-----------------------------------------------------------------------
!
    rucvar='VREL'
    CALL gd_rdat  ( iflno, timep,                                       &
                    level, iprcord, rucvar,                             &
                    v_ext(1,1,klev), nxgem, nygem, ighdr, iret )
    IF  ( iret == 0 )  THEN
      PRINT *, ' tpast v(40,40) = ',v_ext(40,40,klev)
      IF( nxgem /= nx_ext .OR. nygem /= ny_ext ) THEN
        WRITE(6,'(a,/,a,2i12,/,a,2i12)')                                &
            '  Error in grid dimensions.',                              &
            '  GEMPAK returned nx,ny:',nxgem,nygem,                     &
            '  Expected        nx,ny:',nx_ext,ny_ext
        istatus=-31
        RETURN
      END IF
      IF(tmint) THEN
        CALL gd_rdat  ( iflno, timef,                                   &
                        level, iprcord, rucvar,                         &
                        tem1_ext, nxgem, nygem, ighdr, iret )
        IF  ( iret == 0 )  THEN
          PRINT *, ' tfutr v(40,40) = ',tem1_ext(40,40,1)
          DO j=1,ny_ext
            DO i=1,nx_ext
              v_ext(i,j,klev)=wgtp*v_ext(i,j,klev)+                     &
                              wgtf*tem1_ext(i,j,1)
            END DO
          END DO
        END IF
      END IF
    ELSE
      WRITE(6,'(a,a,/,a,i4)') '  Error reading file ',gdfile,           &
              '  GEMPAK GD_RDAT return status: ',iret
      istatus=iret
      RETURN
    END IF
  END DO
!
!-----------------------------------------------------------------------
!
!  Extrapolate data to fill in data below 1000 mb
!  To begin, set all to be equal to values at lowest
!  available MAPS/RUC level.
!
!-----------------------------------------------------------------------
!
  DO klev=1,kstart-1
    DO j=1,ny_ext
      DO i=1,nx_ext
        hgt_ext(i,j,klev)=hgt_ext(i,j,kstart)
        t_ext(i,j,klev)  =t_ext(i,j,kstart)
        qv_ext(i,j,klev) =qv_ext(i,j,kstart)
        u_ext(i,j,klev)  =u_ext(i,j,kstart)
        v_ext(i,j,klev)  =v_ext(i,j,kstart)
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Extrapolate temperature using standard atmospheric
!  lapse rate.
!
!-----------------------------------------------------------------------
!
  DO klev=1,kstart-1
    DO j=1,ny_ext
      DO i=1,nx_ext
        t_ext(i,j,klev)=t_ext(i,j,kstart)*                              &
               ( p_ext(i,j,klev) /                                      &
                 p_ext(i,j,kstart) )**xconst
      END DO
    END DO
    PRINT *, ' pr,tks,text = ',ipr_ext(klev),                           &
               t_ext(60,36,kstart),t_ext(60,36,klev)
  END DO
!
!-----------------------------------------------------------------------
!
!  Change RUC Relative Humidity to Mixing Ratio for ARPS
!  RH is 0-100% so it is multipled by 0.01
!
!-----------------------------------------------------------------------
!
  DO klev=1,nz_ext
    DO j=1,ny_ext
      DO i=1,nx_ext
        qvsat = f_qvsatl( p_ext(i,j,klev), t_ext(i,j,klev) )
        qv_ext(i,j,klev)=0.01*qv_ext(i,j,klev)*qvsat
      END DO
    END DO
  END DO
  PRINT *, ipr_ext(2),' mixing ratio ',qv_ext(60,36,2)
  PRINT *, ipr_ext(4),' mixing ratio ',qv_ext(60,36,4)
  PRINT *, ipr_ext(6),' mixing ratio ',qv_ext(60,36,6)
  PRINT *, ipr_ext(8),' mixing ratio ',qv_ext(60,36,8)
!
!-----------------------------------------------------------------------
!
!  Set height field
!  by integrating T down from kstart level
!
!-----------------------------------------------------------------------
!
  DO klev=kstart-1,1,-1
    DO j=1,ny_ext
      DO i=1,nx_ext
        dln=ALOG(p_ext(i,j,klev)/p_ext(i,j,klev+1))
        const=dln*rddg
        tbar=0.5*(t_ext(i,j,klev)+t_ext(i,j,klev+1))
        hgt_ext(i,j,klev)=hgt_ext(i,j,klev+1)-const*tbar
      END DO
    END DO
    PRINT *, ' pr,hgt(ks),hgt = ',ipr_ext(klev),                        &
               hgt_ext(60,36,kstart),hgt_ext(60,36,klev)
  END DO
!
!-----------------------------------------------------------------------
!
!  Rotate winds to be relative to true north.
!  The RUC data are sent as grid-relative.
!
!-----------------------------------------------------------------------
!
  DO klev=1,nz_ext
!2001-05-16 GMB: Having umap & uear (or vmap & vear) point to
!the same array causes numerical errors when optimizing.
    CALL uvmptoe(nx_ext,ny_ext,u_ext(1,1,klev),v_ext(1,1,klev),         &
                 lon_ext,utmp,vtmp)
    u_ext(:,:,klev) = utmp(:,:)
    v_ext(:,:,klev) = vtmp(:,:)
  END DO
!
!-----------------------------------------------------------------------
!
!  Fill qc,qr,qi,qs,qh arrays with missing value.
!
!-----------------------------------------------------------------------
!
!  DO k=1,nz_ext
!    DO j=1,ny_ext
!      DO i=1,nx_ext
!        qc_ext(i,j,k)=-999.
!        qr_ext(i,j,k)=-999.
!        qi_ext(i,j,k)=-999.
!        qs_ext(i,j,k)=-999.
!        qh_ext(i,j,k)=-999.
!      END DO
!    END DO
!  END DO
!
!-----------------------------------------------------------------------
!
!  Reset map projection to previous values
!
!-----------------------------------------------------------------------
!
  CALL setmapr(iproj,scale,latnot,trlon)
  CALL setorig(1,x0,y0)
!
!-----------------------------------------------------------------------
!
!  Set good status
!
!-----------------------------------------------------------------------
!
  istatus=1

  DEALLOCATE(utmp)
  DEALLOCATE(vtmp)

  RETURN
END SUBROUTINE getgemruc2
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GETGEMETA                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE getgemeta(nx_ext,ny_ext,nz_ext,dir_extd,extdname,            &
           extdinit,extdfcst,julfname,i4time,                           &
           iproj_ext,scale_ext,                                         &
           trlon_ext,latnot_ext,x0_ext,y0_ext,                          &
           lat_ext,lon_ext,                                             &
           p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,                      &
           istatus, tem1_ext)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Reads external file for processing by ext2arps, a program
!  that converts external files to ARPS variables and format.
!
!  This version of rdextfil reads NCEP Eta data in GEMPAK format
!
!  The script Gemenviron must be run by the process running
!  this program.  It defines GEMPAK symbolic links.
!
!  Because of the GEMPAK parameter include file, implicit none
!  is not used in this program.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  Based on GETGEMRUC
!  June, 1996
!
!  MODIFICATIONS:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    extdinit      Initialized time in mm-dd-yyyy:hh:mm:ss format
!    extdfcst      Forecast hour in HHH:MM:SS format
!    julfname      File name in yyjjjhhmm format
!    i4time        Absolute time in seconds (from 1 Jan 1970) of desired data
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
!    qc_ext        Cloud water mixing ratio (kg/kg)
!    qr_ext        Rain  water mixing ratio (kg/kg)
!    qi_ext        Ice         mixing ratio (kg/kg)
!    qs_ext        Snow        mixing ratio (kg/kg)
!    qh_ext        Hail        mixing ratio (kg/kg)
!    u_ext         u wind component (m/s)
!    v_ext         v wind component (m/s)
!    istatus       status indicator
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INCLUDE 'GEMINC:GEMPRM.PRM'
  INTEGER :: nx_ext, ny_ext, nz_ext
!
  CHARACTER (LEN=* ) :: dir_extd
  CHARACTER (LEN=* ) :: extdname

  CHARACTER (LEN=* ) :: extdinit
  CHARACTER (LEN=* ) :: extdfcst
  CHARACTER (LEN=* ) :: julfname
  INTEGER :: i4time
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
  REAL :: x_ext(nx_ext),y_ext(ny_ext)
!
!------------------------------------------------------------------------
!
!  Eta Grid definition parameters
!
!------------------------------------------------------------------------
!
  INTEGER :: iproj_eta
  REAL :: latnot1_eta,latnot2_eta,trlon_eta,scale_eta
  PARAMETER (iproj_eta  =1,  & ! Polar Stereographic
         latnot1_eta=40.,                                               &
             latnot2_eta=0.,                                            &
             trlon_eta  =-105.,                                         &
             scale_eta  =1.0)
  REAL :: dx_eta,dy_eta,swlat_eta,swlon_eta
  PARAMETER (dx_eta = 80000., dy_eta=80000.,                            &
             swlat_eta = 17.53,                                         &
             swlon_eta = -129.30)
  INTEGER :: fhrinc
  PARAMETER (fhrinc=6)
!
!------------------------------------------------------------------------
!
!  Eta variables
!
!------------------------------------------------------------------------
!
  INTEGER :: nz_ext_mx
  PARAMETER(nz_ext_mx=41)
  INTEGER :: ipr_ext(nz_ext_mx)
  DATA ipr_ext / 1000,975,950,                                          &
                  900,875,850,825,800,775,750,725,                      &
                  700,675,650,625,600,575,550,525,                      &
                  500,475,450,425,400,375,350,325,                      &
                  300,275,250,225,200,175,150,125,                      &
                  100, 75, 50,  0,  0,  0/
!
!-----------------------------------------------------------------------
!
!  Output external variable arrays
!
!-----------------------------------------------------------------------
!
  REAL :: lat_ext(nx_ext,ny_ext)
  REAL :: lon_ext(nx_ext,ny_ext)
  REAL :: p_ext(nx_ext,ny_ext,nz_ext)
  REAL :: hgt_ext(nx_ext,ny_ext,nz_ext)
  REAL :: t_ext(nx_ext,ny_ext,nz_ext)
  REAL :: qv_ext(nx_ext,ny_ext,nz_ext)
  REAL :: u_ext(nx_ext,ny_ext,nz_ext)
  REAL :: v_ext(nx_ext,ny_ext,nz_ext)
!  REAL :: qc_ext(nx_ext,ny_ext,nz_ext)    ! Cloud H2O mixing ratio (kg/kg)
!  REAL :: qr_ext(nx_ext,ny_ext,nz_ext)    ! Rain  H2O mixing ratio (kg/kg)
!  REAL :: qi_ext(nx_ext,ny_ext,nz_ext)    ! Ice   H2O mixing ratio (kg/kg)
!  REAL :: qs_ext(nx_ext,ny_ext,nz_ext)    ! Snow  H2O mixing ratio (kg/kg)
!  REAL :: qh_ext(nx_ext,ny_ext,nz_ext)    ! Hail  H2O mixing ratio (kg/kg)
  REAL :: tem1_ext(nx_ext,ny_ext,nz_ext)

  INTEGER :: istatus
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
!------------------------------------------------------------------------
!
!  GEMPAK variables
!
!------------------------------------------------------------------------
!
  INTEGER :: inocord,iprcord,ithcord,ihtcord
  PARAMETER (inocord =0, iprcord =1,                                    &
             ithcord =2, ihtcord =3)
!
  REAL :: grid(llnanl)
  REAL :: rnav(llnnav)
  INTEGER :: nxgem,nygem,iret
  INTEGER :: level(2)
  INTEGER :: ighdr(llgdhd)
  CHARACTER (LEN=20)  :: lastim,timep(2),timef(2)
  CHARACTER (LEN=72)  :: gdcur
  CHARACTER (LEN=256) :: gdfile
!
  CHARACTER (LEN=12) :: etavar
!
!------------------------------------------------------------------------
!
!  Physcial parameters
!
!------------------------------------------------------------------------
!
  REAL :: rd,g
  PARAMETER (rd=287.053, g=9.81)
!
  REAL :: gamma,rddg,xconst
  PARAMETER ( gamma = 0.0065,    & ! degrees/m  lapse rate
          rddg  = (rd/g),                                               &
              xconst = (rd*gamma/g) )
!
!------------------------------------------------------------------------
!
!  Misc internal variables
!
!------------------------------------------------------------------------
!
  CHARACTER (LEN=8) :: gmpktm
  INTEGER :: i,j,k,klev,kstart,itime
  INTEGER :: fhr,fhpast,fhfutr
  INTEGER :: iflno,numgrd,maxgrd,len1,len2
  REAL :: pratio,dln,tbar,const,qvsat,wgtp,wgtf
  INTEGER :: iyr,imo,iday,ihr,imin

  LOGICAL :: tmint,wrtflg,newfil
  LOGICAL :: init_called

  REAL, ALLOCATABLE :: utmp(:,:), vtmp(:,:)
!
!-----------------------------------------------------------------------
!
!  Function f_qvsat and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
  REAL :: f_qvsat

!fpp$ expand (f_qvsat)
!dir$ inline always f_qvsat

  SAVE init_called
  DATA init_called /.false./
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ALLOCATE(utmp(nx_ext,ny_ext))
  ALLOCATE(vtmp(nx_ext,ny_ext))

  iproj_ext=iproj_eta
  scale_ext=scale_eta       ! report lengths in m
  trlon_ext=trlon_eta      ! orientation of external data grids
  latnot_ext(1)=latnot1_eta
  latnot_ext(2)=latnot2_eta
!
  DO k=1,nz_ext
    DO j=1,ny_ext
      DO i=1,nx_ext
        p_ext(i,j,k)=100.*FLOAT(ipr_ext(k))
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Get the lat,lon of the Eta grid points
!
!-----------------------------------------------------------------------
!
  CALL getmapr(iproj,scale,latnot,trlon,x0,y0)
  CALL setmapr(iproj_ext,scale_ext,latnot_ext,trlon_ext)
  CALL lltoxy(1,1,swlat_eta,swlon_eta,x0_ext,y0_ext)
!
  DO i=1,nx_ext
    x_ext(i)=x0_ext+(i-1)*dx_eta
  END DO
  DO j=1,ny_ext
    y_ext(j)=y0_ext+(j-1)*dy_eta
  END DO

  CALL xytoll(nx_ext,ny_ext,x_ext,y_ext,lat_ext,lon_ext)
!
  PRINT *, ' eta point 34,17: ',lat_ext(34,17),lon_ext(34,17)
  PRINT *, ' eta point nx,ny: ',lat_ext(nx_ext,ny_ext),                 &
                                lon_ext(nx_ext,ny_ext)
!
!-----------------------------------------------------------------------
!
!  Initialize GEMPAK sans TAE
!
!-----------------------------------------------------------------------
!
  IF( .NOT.init_called) THEN
    CALL in_bdta(iret)
    init_called=.true.
  END IF
!
!-----------------------------------------------------------------------
!
!  Build Eta file name
!
!-----------------------------------------------------------------------
!
  READ(extdinit,'(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)')                   &
                  iyr,imo,iday,ihr,imin
  myr=MOD(iyr,100)
  len1=LEN(dir_extd)
  len2=len1
  CALL strlnth( dir_extd, len2 )
  IF( dir_extd(len2:len2) /= '/' .AND.  len2 < len1) THEN
    len2=len2+1
    dir_extd(len2:len2)='/'
  END IF

  len1 = LEN( extdname )
  CALL strlnth( extdname, len1 )

  WRITE(gdfile,'(a,a,a,i4.4,i2.2,i2.2,i2.2)')                           &
                 dir_extd(1:len2),extdname(1:len1),'.',                 &
                 iyr,imo,iday,ihr
!
!-----------------------------------------------------------------------
!
!  Open the grid file
!
!-----------------------------------------------------------------------
!
  PRINT *, ' Opening gdfile= ',gdfile
  CALL gr_open  ( gdfile, wrtflg , gdcur, iflno, lastim,                &
                  grid, rnav, numgrd, maxgrd, newfil, iret)
  IF  ( iret == 0 )  THEN
    PRINT *, ' gdcur = ',gdcur
    PRINT *, ' iflno = ',iflno
    PRINT *, ' lastim = ',lastim
    PRINT *, ' numgrd = ',numgrd
    PRINT *, ' maxgrd = ',maxgrd
    PRINT *, ' newfil = ',newfil
    PRINT *, ' iret = ',iret
  ELSE
    WRITE(6,'(a,a,/,a,i4)') '  Error opening file ',gdfile,             &
                 '  GEMPAK GR_OPEN return status: ',iret
    istatus=iret
    RETURN
  END IF
!
!-----------------------------------------------------------------------
!
!  Build the  GEMPAK grid time string
!  It has format yymodd/hhmnFHHH
!  yy: year      mo: month  dd: GMT day
!  hh: GMT hour  mn: minute
!  F: separation charcter
!  HHH: forecast hour (000 = analysis)
!  example  time(1)='950126/1200F000'
!
!-----------------------------------------------------------------------
!
  timep(1)='                    '
  timep(2)='                    '
  timef(1)='                    '
  timef(2)='                    '
  READ(extdfcst,'(i3)') fhr
  IF(MOD(fhr,fhrinc) /= 0) THEN
    tmint=.true.
    fhpast=(fhr/fhrinc)*fhrinc
    fhfutr=fhpast+fhrinc
    wgtp=FLOAT(fhfutr-fhr)/FLOAT(fhrinc)
    wgtf=1.-wgtp
  ELSE
    tmint=.false.
    fhpast=fhr
    fhfutr=fhr
    wgtp=1.0
    wgtf=0.0
  END IF

  WRITE(timep(1),'(i2.2,i2.2,i2.2,a1,i2.2,i2.2,a1,i3.3)')               &
        myr,imo,iday,'/',ihr,imin,'F',fhpast
  WRITE(6,'(a,a)') '  Past   GEMPAK time string ',timep(1)
  WRITE(timef(1),'(i2.2,i2.2,i2.2,a1,i2.2,i2.2,a1,i3.3)')               &
        myr,imo,iday,'/',ihr,imin,'F',fhfutr
  WRITE(6,'(a,a)') '  Future GEMPAK time string ',timef(1)
!
!-----------------------------------------------------------------------
!
!  Data in the RUC files are only available from 100-1000 mb.
!  Thus this process does not start at pr_ext vertical level one,
!  which is at 1100 mb.  Find 1000 mb in the ipr_ext vector.
!
!-----------------------------------------------------------------------
!
  DO klev=1,nz_ext-1
    IF(ipr_ext(klev) <= 1000) EXIT
  END DO
  kstart=klev
  PRINT *, ' kstart = ',kstart
!
!-----------------------------------------------------------------------
!
!  Go through each of the RUC variables collecting those
!  interpolated to pressure surfaces at NMC from the
!  original model coordinates (hybrid sigma-isentropic)
!
!-----------------------------------------------------------------------
!
  DO klev=kstart,nz_ext
    PRINT *, ' Reading level ',ipr_ext(klev)
    level(1)=ipr_ext(klev)
    level(2)=-1
!
!-----------------------------------------------------------------------
!
!  Heights
!
!-----------------------------------------------------------------------
!
    etavar='HGHT'
    CALL gd_rdat  ( iflno, timep,                                       &
                    level, iprcord, etavar,                             &
                    hgt_ext(1,1,klev), nxgem, nygem, ighdr, iret )
    IF  ( iret == 0 )  THEN
      PRINT *, ' tpast hght(40,40) = ',hgt_ext(40,40,klev)
      IF( nxgem /= nx_ext .OR. nygem /= ny_ext ) THEN
        WRITE(6,'(a,/,a,2i12,/,a,2i12)')                                &
            '  Error in grid dimensions.',                              &
            '  GEMPAK returned nx,ny:',nxgem,nygem,                     &
            '  Expected        nx,ny:',nx_ext,ny_ext
        istatus=-31
        RETURN
      END IF
      IF(tmint) THEN
        CALL gd_rdat  ( iflno, timef,                                   &
                        level, iprcord, etavar,                         &
                        tem1_ext, nxgem, nygem, ighdr, iret )
        IF  ( iret == 0 )  THEN
          PRINT *, ' tfutr hght(40,40) = ',tem1_ext(40,40,1)
          DO j=1,ny_ext
            DO i=1,nx_ext
              hgt_ext(i,j,klev)=wgtp*hgt_ext(i,j,klev)+                 &
                                wgtf*tem1_ext(i,j,1)
            END DO
          END DO
        END IF
      END IF
    ELSE
      WRITE(6,'(a,a,/,a,i4)') '  Error reading file ',gdfile,           &
              '  GEMPAK GD_RDAT return status: ',iret
      istatus=iret
      RETURN
    END IF
!
!-----------------------------------------------------------------------
!
!  Temperatures
!
!-----------------------------------------------------------------------
!
    etavar='TMPK'
    CALL gd_rdat  ( iflno, timep,                                       &
                    level, iprcord, etavar,                             &
                    t_ext(1,1,klev), nxgem, nygem, ighdr, iret )
    IF  ( iret == 0 )  THEN
      PRINT *, ' tpast t(40,40) = ',t_ext(40,40,klev)
      IF( nxgem /= nx_ext .OR. nygem /= ny_ext ) THEN
        WRITE(6,'(a,/,a,2i12,/,a,2i12)')                                &
            '  Error in grid dimensions.',                              &
            '  GEMPAK returned nx,ny:',nxgem,nygem,                     &
            '  Expected        nx,ny:',nx_ext,ny_ext
        istatus=-31
        RETURN
      END IF
      IF(tmint) THEN
        CALL gd_rdat  ( iflno, timef,                                   &
                        level, iprcord, etavar,                         &
                        tem1_ext, nxgem, nygem, ighdr, iret )
        IF  ( iret == 0 )  THEN
          PRINT *, ' tfutr t(40,40) = ',tem1_ext(40,40,1)
          DO j=1,ny_ext
            DO i=1,nx_ext
              t_ext(i,j,klev)=wgtp*t_ext(i,j,klev)+                     &
                              wgtf*tem1_ext(i,j,1)
            END DO
          END DO
        END IF
      END IF
    ELSE
      WRITE(6,'(a,a,/,a,i4)') '  Error reading file ',gdfile,           &
              '  GEMPAK GD_RDAT return status: ',iret
      istatus=iret
      RETURN
    END IF
!
!-----------------------------------------------------------------------
!
!  Relative humidity
!
!-----------------------------------------------------------------------
!
    etavar='SPFH'
    CALL gd_rdat  ( iflno, timep,                                       &
                    level, iprcord, etavar,                             &
                    qv_ext(1,1,klev), nxgem, nygem, ighdr, iret )
    IF  ( iret == 0 )  THEN
      PRINT *, ' tpast qv(40,40) = ',qv_ext(40,40,klev)
      IF( nxgem /= nx_ext .OR. nygem /= ny_ext ) THEN
        WRITE(6,'(a,/,a,2i12,/,a,2i12)')                                &
            '  Error in grid dimensions.',                              &
            '  GEMPAK returned nx,ny:',nxgem,nygem,                     &
            '  Expected        nx,ny:',nx_ext,ny_ext
        istatus=-31
        RETURN
      END IF
      IF(tmint) THEN
        CALL gd_rdat  ( iflno, timef,                                   &
                        level, iprcord, etavar,                         &
                        tem1_ext, nxgem, nygem, ighdr, iret )
        IF  ( iret == 0 )  THEN
          PRINT *, ' tfutr qv(40,40) = ',tem1_ext(40,40,1)
          DO j=1,ny_ext
            DO i=1,nx_ext
              qv_ext(i,j,klev)=wgtp*qv_ext(i,j,klev)+                   &
                               wgtf*tem1_ext(i,j,1)
            END DO
          END DO
        END IF
      END IF
    ELSE
      WRITE(6,'(a,a,/,a,i4)') '  Error reading file ',gdfile,           &
              '  GEMPAK GD_RDAT return status: ',iret
      istatus=iret
      RETURN
    END IF
!
!-----------------------------------------------------------------------
!
!  U - Grid relative velocities
!
!-----------------------------------------------------------------------
!
    etavar='UREL'
    CALL gd_rdat  ( iflno, timep,                                       &
                    level, iprcord, etavar,                             &
                    u_ext(1,1,klev), nxgem, nygem, ighdr, iret )
    IF  ( iret == 0 )  THEN
      PRINT *, ' tpast u(40,40) = ',u_ext(40,40,klev)
      IF( nxgem /= nx_ext .OR. nygem /= ny_ext ) THEN
        WRITE(6,'(a,/,a,2i12,/,a,2i12)')                                &
            '  Error in grid dimensions.',                              &
            '  GEMPAK returned nx,ny:',nxgem,nygem,                     &
            '  Expected        nx,ny:',nx_ext,ny_ext
        istatus=-31
        RETURN
      END IF
      IF(tmint) THEN
        CALL gd_rdat  ( iflno, timef,                                   &
                        level, iprcord, etavar,                         &
                        tem1_ext, nxgem, nygem, ighdr, iret )
        IF  ( iret == 0 )  THEN
          PRINT *, ' tfutr u(40,40) = ',tem1_ext(40,40,1)
          DO j=1,ny_ext
            DO i=1,nx_ext
              u_ext(i,j,klev)=wgtp*u_ext(i,j,klev)+                     &
                              wgtf*tem1_ext(i,j,1)
            END DO
          END DO
        END IF
      END IF
    ELSE
      WRITE(6,'(a,a,/,a,i4)') '  Error reading file ',gdfile,           &
              '  GEMPAK GD_RDAT return status: ',iret
      istatus=iret
      RETURN
    END IF
!
!-----------------------------------------------------------------------
!
!  V - Grid relative velocities
!
!-----------------------------------------------------------------------
!
    etavar='VREL'
    CALL gd_rdat  ( iflno, timep,                                       &
                    level, iprcord, etavar,                             &
                    v_ext(1,1,klev), nxgem, nygem, ighdr, iret )
    IF  ( iret == 0 )  THEN
      PRINT *, ' tpast v(40,40) = ',v_ext(40,40,klev)
      IF( nxgem /= nx_ext .OR. nygem /= ny_ext ) THEN
        WRITE(6,'(a,/,a,2i12,/,a,2i12)')                                &
            '  Error in grid dimensions.',                              &
            '  GEMPAK returned nx,ny:',nxgem,nygem,                     &
            '  Expected        nx,ny:',nx_ext,ny_ext
        istatus=-31
        RETURN
      END IF
      IF(tmint) THEN
        CALL gd_rdat  ( iflno, timef,                                   &
                        level, iprcord, etavar,                         &
                        tem1_ext, nxgem, nygem, ighdr, iret )
        IF  ( iret == 0 )  THEN
          PRINT *, ' tfutr v(40,40) = ',tem1_ext(40,40,1)
          DO j=1,ny_ext
            DO i=1,nx_ext
              v_ext(i,j,klev)=wgtp*v_ext(i,j,klev)+                     &
                              wgtf*tem1_ext(i,j,1)
            END DO
          END DO
        END IF
      END IF
    ELSE
      WRITE(6,'(a,a,/,a,i4)') '  Error reading file ',gdfile,           &
              '  GEMPAK GD_RDAT return status: ',iret
      istatus=iret
      RETURN
    END IF
  END DO
!
!-----------------------------------------------------------------------
!
!  Extrapolate data to fill in data below 1000 mb
!  To begin, set all to be equal to values at lowest
!  available MAPS/RUC level.
!
!-----------------------------------------------------------------------
!
  DO klev=1,kstart-1
    DO j=1,ny_ext
      DO i=1,nx_ext
        hgt_ext(i,j,klev)=hgt_ext(i,j,kstart)
        t_ext(i,j,klev)  =t_ext(i,j,kstart)
        qv_ext(i,j,klev) =qv_ext(i,j,kstart)
        u_ext(i,j,klev)  =u_ext(i,j,kstart)
        v_ext(i,j,klev)  =v_ext(i,j,kstart)
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Extrapolate temperature using standard atmospheric
!  lapse rate.
!
!-----------------------------------------------------------------------
!
  DO klev=1,kstart-1
    DO j=1,ny_ext
      DO i=1,nx_ext
        t_ext(i,j,klev)=t_ext(i,j,kstart)*                              &
               ( p_ext(i,j,klev) /                                      &
                 p_ext(i,j,kstart) )**xconst
      END DO
    END DO
    PRINT *, ' pr,tks,text = ',ipr_ext(klev),                           &
               t_ext(60,36,kstart),t_ext(60,36,klev)
  END DO
!
!-----------------------------------------------------------------------
!
!  Change RUC Relative Humidity to Mixing Ratio for ARPS
!  RH is 0-100% so it is multipled by 0.01
!
!-----------------------------------------------------------------------
!
!    DO 425 klev=1,nz_ext
!    DO 425 j=1,ny_ext
!    DO 425 i=1,nx_ext
!      qvsat = f_qvsat( p_ext(i,j,klev), t_ext(i,j,klev) )
!      qv_ext(i,j,klev)=0.01*qv_ext(i,j,klev)*qvsat
! 425   CONTINUE
!    print *, ipr_ext(2),' mixing ratio ',qv_ext(60,36,2)
!    print *, ipr_ext(4),' mixing ratio ',qv_ext(60,36,4)
!    print *, ipr_ext(6),' mixing ratio ',qv_ext(60,36,6)
!    print *, ipr_ext(8),' mixing ratio ',qv_ext(60,36,8)
!
!-----------------------------------------------------------------------
!
!  Set height field
!  by integrating T down from kstart level
!
!-----------------------------------------------------------------------
!
  DO klev=kstart-1,1,-1
    DO j=1,ny_ext
      DO i=1,nx_ext
        dln=ALOG(p_ext(i,j,klev)/p_ext(i,j,klev+1))
        const=dln*rddg
        tbar=0.5*(t_ext(i,j,klev)+t_ext(i,j,klev+1))
        hgt_ext(i,j,klev)=hgt_ext(i,j,klev+1)-const*tbar
      END DO
    END DO
    PRINT *, ' pr,hgt(ks),hgt = ',ipr_ext(klev),                        &
               hgt_ext(60,36,kstart),hgt_ext(60,36,klev)
  END DO
!
!-----------------------------------------------------------------------
!
!  Rotate winds to be relative to true north.
!  The RUC data are sent as grid-relative.
!
!-----------------------------------------------------------------------
!
  DO klev=1,nz_ext
!2001-05-16 GMB: Having umap & uear (or vmap & vear) point to
!the same array causes numerical errors when optimizing.
    CALL uvmptoe(nx_ext,ny_ext,u_ext(1,1,klev),v_ext(1,1,klev),         &
                 lon_ext,utmp,vtmp)
    u_ext(:,:,klev) = utmp(:,:)
    v_ext(:,:,klev) = vtmp(:,:)
  END DO
!
!-----------------------------------------------------------------------
!
!  Fill qc,qr,qi,qs,qh arrays with missing value.
!
!-----------------------------------------------------------------------
!
!  DO k=1,nz_ext
!    DO j=1,ny_ext
!      DO i=1,nx_ext
!        qc_ext(i,j,k)=-999.
!        qr_ext(i,j,k)=-999.
!        qi_ext(i,j,k)=-999.
!        qs_ext(i,j,k)=-999.
!        qh_ext(i,j,k)=-999.
!      END DO
!    END DO
!  END DO
!
!-----------------------------------------------------------------------
!
!  Reset map projection to previous values
!
!-----------------------------------------------------------------------
!
  CALL setmapr(iproj,scale,latnot,trlon)
  CALL setorig(1,x0,y0)
!
!-----------------------------------------------------------------------
!
!  Set good status
!
!-----------------------------------------------------------------------
!
  istatus=1

  DEALLOCATE(utmp)
  DEALLOCATE(vtmp)

  RETURN
END SUBROUTINE getgemeta
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GETGEMETA2                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE getgemeta2(nx_ext,ny_ext,nz_ext,dir_extd,extdname,           &
           extdinit,extdfcst,julfname,i4time,                           &
           iproj_ext,scale_ext,                                         &
           trlon_ext,latnot_ext,x0_ext,y0_ext,                          &
           lat_ext,lon_ext,                                             &
           p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,                      &
           istatus, tem1_ext)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Reads external file for processing by ext2arps, a program
!  that converts external files to ARPS variables and format.
!
!  This version of rdextfil reads NCEP Eta data in GEMPAK format
!
!  The script Gemenviron must be run by the process running
!  this program.  It defines GEMPAK symbolic links.
!
!  Because of the GEMPAK parameter include file, implicit none
!  is not used in this program.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  Based on GETGEMRUC
!  June, 1996
!
!  MODIFICATIONS:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    extdinit      Initialized time in mm-dd-yyyy:hh:mm:ss format
!    extdfcst      Forecast hour in HHH:MM:SS format
!    julfname      File name in yyjjjhhmm format
!    i4time        Absolute time in seconds (from 1 Jan 1970) of desired data
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
!    qc_ext        Cloud water mixing ratio (kg/kg)
!    qr_ext        Rain  water mixing ratio (kg/kg)
!    qi_ext        Ice         mixing ratio (kg/kg)
!    qs_ext        Snow        mixing ratio (kg/kg)
!    qh_ext        Hail        mixing ratio (kg/kg)
!    u_ext         u wind component (m/s)
!    v_ext         v wind component (m/s)
!    istatus       status indicator
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INCLUDE 'GEMINC:GEMPRM.PRM'
  INTEGER :: nx_ext, ny_ext, nz_ext
!
  CHARACTER (LEN=* ) :: dir_extd
  CHARACTER (LEN=* ) :: extdname

  CHARACTER (LEN=* ) :: extdinit
  CHARACTER (LEN=* ) :: extdfcst
  CHARACTER (LEN=* ) :: julfname
  INTEGER :: i4time
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
  REAL :: x_ext(nx_ext),y_ext(ny_ext)
!
!------------------------------------------------------------------------
!
!  Eta Grid definition parameters
!
!------------------------------------------------------------------------
!
  INTEGER :: iproj_eta
  REAL :: latnot1_eta,latnot2_eta,trlon_eta,scale_eta
  PARAMETER (iproj_eta  =1,  & ! Polar Stereographic
         latnot1_eta=40.15,                                             &
             latnot2_eta=0.,                                            &
             trlon_eta  =-105.,                                         &
             scale_eta  =1.0)
  REAL :: dx_eta,dy_eta,swlat_eta,swlon_eta
  PARAMETER (dx_eta = 80000., dy_eta=80000.,                            &
             swlat_eta = -0.27,                                         &
             swlon_eta = -139.48)
  INTEGER :: fhrinc
  PARAMETER (fhrinc=6)
!
!------------------------------------------------------------------------
!
!  Eta variables
!
!------------------------------------------------------------------------
!
  INTEGER :: nz_ext_mx
  PARAMETER(nz_ext_mx=41)

  INTEGER :: ipr_ext(nz_ext_mx)
  DATA ipr_ext / 1000, 975, 950, 925,                                   &
                  900, 875, 850, 825, 800,775,750,725,                  &
                  700, 675, 650, 625, 600,575,550,525,                  &
                  500, 475, 450, 425, 400,375,350,325,                  &
                  300, 275, 250, 225, 200,175,150,125,                  &
                  100,  75,  50,   0,   0/
!
!-----------------------------------------------------------------------
!
!  Output external variable arrays
!
!-----------------------------------------------------------------------
!
  REAL :: lat_ext(nx_ext,ny_ext)
  REAL :: lon_ext(nx_ext,ny_ext)
  REAL :: p_ext(nx_ext,ny_ext,nz_ext)
  REAL :: hgt_ext(nx_ext,ny_ext,nz_ext)
  REAL :: t_ext(nx_ext,ny_ext,nz_ext)
  REAL :: qv_ext(nx_ext,ny_ext,nz_ext)
  REAL :: u_ext(nx_ext,ny_ext,nz_ext)
  REAL :: v_ext(nx_ext,ny_ext,nz_ext)
!  REAL :: qc_ext(nx_ext,ny_ext,nz_ext)    ! Cloud H2O mixing ratio (kg/kg)
!  REAL :: qr_ext(nx_ext,ny_ext,nz_ext)    ! Rain  H2O mixing ratio (kg/kg)
!  REAL :: qi_ext(nx_ext,ny_ext,nz_ext)    ! Ice   H2O mixing ratio (kg/kg)
!  REAL :: qs_ext(nx_ext,ny_ext,nz_ext)    ! Snow  H2O mixing ratio (kg/kg)
!  REAL :: qh_ext(nx_ext,ny_ext,nz_ext)    ! Hail  H2O mixing ratio (kg/kg)
  REAL :: tem1_ext(nx_ext,ny_ext,nz_ext)

  INTEGER :: istatus
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
!------------------------------------------------------------------------
!
!  GEMPAK variables
!
!------------------------------------------------------------------------
!
  INTEGER :: inocord,iprcord,ithcord,ihtcord
  PARAMETER (inocord =0, iprcord =1,                                    &
             ithcord =2, ihtcord =3)
!
  REAL :: grid(llnanl)
  REAL :: rnav(llnnav)
  INTEGER :: nxgem,nygem,iret
  INTEGER :: level(2)
  INTEGER :: ighdr(llgdhd)
  CHARACTER (LEN=20) :: lastim,timep(2),timef(2)
  CHARACTER (LEN=72) :: gdcur
  CHARACTER (LEN=256) :: gdfile
!
  CHARACTER (LEN=12) :: etavar
!
!------------------------------------------------------------------------
!
!  Physcial parameters
!
!------------------------------------------------------------------------
!
  REAL :: rd,g
  PARAMETER (rd=287.053, g=9.81)
!
  REAL :: gamma,rddg,xconst
  PARAMETER ( gamma = 0.0065,    & ! degrees/m  lapse rate
          rddg  = (rd/g),                                               &
              xconst = (rd*gamma/g) )
!
!------------------------------------------------------------------------
!
!  Misc internal variables
!
!------------------------------------------------------------------------
!
  CHARACTER (LEN=8) :: gmpktm
  INTEGER :: i,j,k,klev,kstart,itime
  INTEGER :: fhr,fhpast,fhfutr
  INTEGER :: iflno,numgrd,maxgrd,len1,len2
  REAL :: pratio,dln,tbar,const,qvsat,wgtp,wgtf
  INTEGER :: iyr,imo,iday,ihr,imin

  LOGICAL :: tmint,wrtflg,newfil
  LOGICAL :: init_called

  REAL, ALLOCATABLE :: utmp(:,:), vtmp(:,:)
!
!-----------------------------------------------------------------------
!
!  Function f_qvsat and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
  REAL :: f_qvsat

!fpp$ expand (f_qvsat)
!dir$ inline always f_qvsat

  SAVE init_called
  DATA init_called /.false./
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ALLOCATE(utmp(nx_ext,ny_ext))
  ALLOCATE(vtmp(nx_ext,ny_ext))

  iproj_ext=iproj_eta
  scale_ext=scale_eta       ! report lengths in m
  trlon_ext=trlon_eta      ! orientation of external data grids
  latnot_ext(1)=latnot1_eta
  latnot_ext(2)=latnot2_eta
!
  DO k=1,nz_ext
    DO j=1,ny_ext
      DO i=1,nx_ext
        p_ext(i,j,k)=100.*FLOAT(ipr_ext(k))
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Get the lat,lon of the Eta grid points
!
!-----------------------------------------------------------------------
!
  CALL getmapr(iproj,scale,latnot,trlon,x0,y0)
  CALL setmapr(iproj_ext,scale_ext,latnot_ext,trlon_ext)
  CALL lltoxy(1,1,swlat_eta,swlon_eta,x0_ext,y0_ext)
!
  DO i=1,nx_ext
    x_ext(i)=x0_ext+(i-1)*dx_eta
  END DO
  DO j=1,ny_ext
    y_ext(j)=y0_ext+(j-1)*dy_eta
  END DO

  CALL xytoll(nx_ext,ny_ext,x_ext,y_ext,lat_ext,lon_ext)
!
  PRINT *, ' eta point 34,17: ',lat_ext(34,17),lon_ext(34,17)
  PRINT *, ' eta point nx,ny: ',lat_ext(nx_ext,ny_ext),                 &
                                lon_ext(nx_ext,ny_ext)
  PRINT *, ' above should be: ', 32.75,  -14.60
!
!-----------------------------------------------------------------------
!
!  Initialize GEMPAK sans TAE
!
!-----------------------------------------------------------------------
!
  IF( .NOT.init_called) THEN
    CALL in_bdta(iret)
    init_called=.true.
  END IF
!
!-----------------------------------------------------------------------
!
!  Build Eta file name
!
!-----------------------------------------------------------------------
!
  READ(extdinit,'(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)')                   &
                  iyr,imo,iday,ihr,imin
  myr=MOD(iyr,100)
  len1=LEN(dir_extd)
  len2=len1
  CALL strlnth( dir_extd, len2 )
  IF( dir_extd(len2:len2) /= '/' .AND.  len2 < len1) THEN
    len2=len2+1
    dir_extd(len2:len2)='/'
  END IF

  len1 = LEN( extdname )
  CALL strlnth( extdname, len1 )

  WRITE(gdfile,'(a,a,a,i4.4,i2.2,i2.2,i2.2)')                           &
                 dir_extd(1:len2),extdname(1:len1),'.',                 &
                 iyr,imo,iday,ihr
!
!-----------------------------------------------------------------------
!
!  Open the grid file
!
!-----------------------------------------------------------------------
!
  PRINT *, ' Opening gdfile= ',gdfile
  CALL gr_open  ( gdfile, wrtflg , gdcur, iflno, lastim,                &
                  grid, rnav, numgrd, maxgrd, newfil, iret)
  IF  ( iret == 0 )  THEN
    PRINT *, ' gdcur = ',gdcur
    PRINT *, ' iflno = ',iflno
    PRINT *, ' lastim = ',lastim
    PRINT *, ' numgrd = ',numgrd
    PRINT *, ' maxgrd = ',maxgrd
    PRINT *, ' newfil = ',newfil
    PRINT *, ' iret = ',iret
  ELSE
    WRITE(6,'(a,a,/,a,i4)') '  Error opening file ',gdfile,             &
                 '  GEMPAK GR_OPEN return status: ',iret
    istatus=iret
    RETURN
  END IF
!
!-----------------------------------------------------------------------
!
!  Build the  GEMPAK grid time string
!  It has format yymodd/hhmnFHHH
!  yy: year      mo: month  dd: GMT day
!  hh: GMT hour  mn: minute
!  F: separation charcter
!  HHH: forecast hour (000 = analysis)
!  example  time(1)='950126/1200F000'
!
!-----------------------------------------------------------------------
!
  timep(1)='                    '
  timep(2)='                    '
  timef(1)='                    '
  timef(2)='                    '
  READ(extdfcst,'(i3)') fhr
  IF(MOD(fhr,fhrinc) /= 0) THEN
    tmint=.true.
    fhpast=(fhr/fhrinc)*fhrinc
    fhfutr=fhpast+fhrinc
    wgtp=FLOAT(fhfutr-fhr)/FLOAT(fhrinc)
    wgtf=1.-wgtp
  ELSE
    tmint=.false.
    fhpast=fhr
    fhfutr=fhr
    wgtp=1.0
    wgtf=0.0
  END IF

  WRITE(timep(1),'(i2.2,i2.2,i2.2,a1,i2.2,i2.2,a1,i3.3)')               &
        myr,imo,iday,'/',ihr,imin,'F',fhpast
  WRITE(6,'(a,a)') '  Past   GEMPAK time string ',timep(1)
  WRITE(timef(1),'(i2.2,i2.2,i2.2,a1,i2.2,i2.2,a1,i3.3)')               &
        myr,imo,iday,'/',ihr,imin,'F',fhfutr
  WRITE(6,'(a,a)') '  Future GEMPAK time string ',timef(1)

!
!-----------------------------------------------------------------------
!
!  Data in the RUC files are only available from 100-1000 mb.
!  Thus this process does not start at pr_ext vertical level one,
!  which is at 1100 mb.  Find 1000 mb in the ipr_ext vector.
!
!-----------------------------------------------------------------------
!
  DO klev=1,nz_ext-1
    IF(ipr_ext(klev) <= 1000) EXIT
  END DO
  kstart=klev
  PRINT *, ' kstart = ',kstart
!
!-----------------------------------------------------------------------
!
!  Go through each of the RUC variables collecting those
!  interpolated to pressure surfaces at NMC from the
!  original model coordinates (hybrid sigma-isentropic)
!
!-----------------------------------------------------------------------
!
  DO klev=kstart,nz_ext
    PRINT *, ' Reading level ',ipr_ext(klev)
    level(1)=ipr_ext(klev)
    level(2)=-1
!
!-----------------------------------------------------------------------
!
!  Heights
!
!-----------------------------------------------------------------------
!
    etavar='HGHT'
    CALL gd_rdat  ( iflno, timep,                                       &
                    level, iprcord, etavar,                             &
                    hgt_ext(1,1,klev), nxgem, nygem, ighdr, iret )
    IF  ( iret == 0 )  THEN
      PRINT *, ' tpast hght(40,40) = ',hgt_ext(40,40,klev)
      IF( nxgem /= nx_ext .OR. nygem /= ny_ext ) THEN
        WRITE(6,'(a,/,a,2i12,/,a,2i12)')                                &
            '  Error in grid dimensions.',                              &
            '  GEMPAK returned nx,ny:',nxgem,nygem,                     &
            '  Expected        nx,ny:',nx_ext,ny_ext
        istatus=-31
        RETURN
      END IF
      IF(tmint) THEN
        CALL gd_rdat  ( iflno, timef,                                   &
                        level, iprcord, etavar,                         &
                        tem1_ext, nxgem, nygem, ighdr, iret )
        IF  ( iret == 0 )  THEN
          PRINT *, ' tfutr hght(40,40) = ',tem1_ext(40,40,1)
          DO j=1,ny_ext
            DO i=1,nx_ext
              hgt_ext(i,j,klev)=wgtp*hgt_ext(i,j,klev)+                 &
                                wgtf*tem1_ext(i,j,1)
            END DO
          END DO
        END IF
      END IF
    ELSE
      WRITE(6,'(a,a,/,a,i4)') '  Error reading file ',gdfile,           &
              '  GEMPAK GD_RDAT return status: ',iret
      istatus=iret
      RETURN
    END IF
!
!-----------------------------------------------------------------------
!
!  Temperatures
!
!-----------------------------------------------------------------------
!
    etavar='TMPK'
    CALL gd_rdat  ( iflno, timep,                                       &
                    level, iprcord, etavar,                             &
                    t_ext(1,1,klev), nxgem, nygem, ighdr, iret )
    IF  ( iret == 0 )  THEN
      PRINT *, ' tpast t(40,40) = ',t_ext(40,40,klev)
      IF( nxgem /= nx_ext .OR. nygem /= ny_ext ) THEN
        WRITE(6,'(a,/,a,2i12,/,a,2i12)')                                &
            '  Error in grid dimensions.',                              &
            '  GEMPAK returned nx,ny:',nxgem,nygem,                     &
            '  Expected        nx,ny:',nx_ext,ny_ext
        istatus=-31
        RETURN
      END IF
      IF(tmint) THEN
        CALL gd_rdat  ( iflno, timef,                                   &
                        level, iprcord, etavar,                         &
                        tem1_ext, nxgem, nygem, ighdr, iret )
        IF  ( iret == 0 )  THEN
          PRINT *, ' tfutr t(40,40) = ',tem1_ext(40,40,1)
          DO j=1,ny_ext
            DO i=1,nx_ext
              t_ext(i,j,klev)=wgtp*t_ext(i,j,klev)+                     &
                              wgtf*tem1_ext(i,j,1)
            END DO
          END DO
        END IF
      END IF
    ELSE
      WRITE(6,'(a,a,/,a,i4)') '  Error reading file ',gdfile,           &
              '  GEMPAK GD_RDAT return status: ',iret
      istatus=iret
      RETURN
    END IF
!
!-----------------------------------------------------------------------
!
!  Relative humidity
!
!-----------------------------------------------------------------------
!
    etavar='SPFH'
    CALL gd_rdat  ( iflno, timep,                                       &
                    level, iprcord, etavar,                             &
                    qv_ext(1,1,klev), nxgem, nygem, ighdr, iret )
    IF  ( iret == 0 )  THEN
      PRINT *, ' tpast rh(40,40) = ',qv_ext(40,40,klev)
      IF( nxgem /= nx_ext .OR. nygem /= ny_ext ) THEN
        WRITE(6,'(a,/,a,2i12,/,a,2i12)')                                &
            '  Error in grid dimensions.',                              &
            '  GEMPAK returned nx,ny:',nxgem,nygem,                     &
            '  Expected        nx,ny:',nx_ext,ny_ext
        istatus=-31
        RETURN
      END IF
      IF(tmint) THEN
        CALL gd_rdat  ( iflno, timef,                                   &
                        level, iprcord, etavar,                         &
                        tem1_ext, nxgem, nygem, ighdr, iret )
        IF  ( iret == 0 )  THEN
          PRINT *, ' tfutr rh(40,40) = ',tem1_ext(40,40,1)
          DO j=1,ny_ext
            DO i=1,nx_ext
              qv_ext(i,j,klev)=wgtp*qv_ext(i,j,klev)+                   &
                               wgtf*tem1_ext(i,j,1)
            END DO
          END DO
        END IF
      END IF
    ELSE
      WRITE(6,'(a,a,/,a,i4)') '  Error reading file ',gdfile,           &
              '  GEMPAK GD_RDAT return status: ',iret
      istatus=iret
      RETURN
    END IF
!
!-----------------------------------------------------------------------
!
!  U - Grid relative velocities
!
!-----------------------------------------------------------------------
!
    etavar='UREL'
    CALL gd_rdat  ( iflno, timep,                                       &
                    level, iprcord, etavar,                             &
                    u_ext(1,1,klev), nxgem, nygem, ighdr, iret )
    IF  ( iret == 0 )  THEN
      PRINT *, ' tpast u(40,40) = ',u_ext(40,40,klev)
      IF( nxgem /= nx_ext .OR. nygem /= ny_ext ) THEN
        WRITE(6,'(a,/,a,2i12,/,a,2i12)')                                &
            '  Error in grid dimensions.',                              &
            '  GEMPAK returned nx,ny:',nxgem,nygem,                     &
            '  Expected        nx,ny:',nx_ext,ny_ext
        istatus=-31
        RETURN
      END IF
      IF(tmint) THEN
        CALL gd_rdat  ( iflno, timef,                                   &
                        level, iprcord, etavar,                         &
                        tem1_ext, nxgem, nygem, ighdr, iret )
        IF  ( iret == 0 )  THEN
          PRINT *, ' tfutr u(40,40) = ',tem1_ext(40,40,1)
          DO j=1,ny_ext
            DO i=1,nx_ext
              u_ext(i,j,klev)=wgtp*u_ext(i,j,klev)+                     &
                              wgtf*tem1_ext(i,j,1)
            END DO
          END DO
        END IF
      END IF
    ELSE
      WRITE(6,'(a,a,/,a,i4)') '  Error reading file ',gdfile,           &
              '  GEMPAK GD_RDAT return status: ',iret
      istatus=iret
      RETURN
    END IF
!
!-----------------------------------------------------------------------
!
!  V - Grid relative velocities
!
!-----------------------------------------------------------------------
!
    etavar='VREL'
    CALL gd_rdat  ( iflno, timep,                                       &
                    level, iprcord, etavar,                             &
                    v_ext(1,1,klev), nxgem, nygem, ighdr, iret )
    IF  ( iret == 0 )  THEN
      PRINT *, ' tpast v(40,40) = ',v_ext(40,40,klev)
      IF( nxgem /= nx_ext .OR. nygem /= ny_ext ) THEN
        WRITE(6,'(a,/,a,2i12,/,a,2i12)')                                &
            '  Error in grid dimensions.',                              &
            '  GEMPAK returned nx,ny:',nxgem,nygem,                     &
            '  Expected        nx,ny:',nx_ext,ny_ext
        istatus=-31
        RETURN
      END IF
      IF(tmint) THEN
        CALL gd_rdat  ( iflno, timef,                                   &
                        level, iprcord, etavar,                         &
                        tem1_ext, nxgem, nygem, ighdr, iret )
        IF  ( iret == 0 )  THEN
          PRINT *, ' tfutr v(40,40) = ',tem1_ext(40,40,1)
          DO j=1,ny_ext
            DO i=1,nx_ext
              v_ext(i,j,klev)=wgtp*v_ext(i,j,klev)+                     &
                              wgtf*tem1_ext(i,j,1)
            END DO
          END DO
        END IF
      END IF
    ELSE
      WRITE(6,'(a,a,/,a,i4)') '  Error reading file ',gdfile,           &
              '  GEMPAK GD_RDAT return status: ',iret
      istatus=iret
      RETURN
    END IF
  END DO
!
!-----------------------------------------------------------------------
!
!  Extrapolate data to fill in data below 1000 mb
!  To begin, set all to be equal to values at lowest
!  available MAPS/RUC level.
!
!-----------------------------------------------------------------------
!
  DO klev=1,kstart-1
    DO j=1,ny_ext
      DO i=1,nx_ext
        hgt_ext(i,j,klev)=hgt_ext(i,j,kstart)
        t_ext(i,j,klev)  =t_ext(i,j,kstart)
        qv_ext(i,j,klev) =qv_ext(i,j,kstart)
        u_ext(i,j,klev)  =u_ext(i,j,kstart)
        v_ext(i,j,klev)  =v_ext(i,j,kstart)
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Extrapolate temperature using standard atmospheric
!  lapse rate.
!
!-----------------------------------------------------------------------
!
  DO klev=1,kstart-1
    DO j=1,ny_ext
      DO i=1,nx_ext
        t_ext(i,j,klev)=t_ext(i,j,kstart)*                              &
               ( p_ext(i,j,klev) /                                      &
                 p_ext(i,j,kstart) )**xconst
      END DO
    END DO
    PRINT *, ' pr,tks,text = ',ipr_ext(klev),                           &
               t_ext(60,36,kstart),t_ext(60,36,klev)
  END DO
!
!-----------------------------------------------------------------------
!
!  Change RUC Relative Humidity to Mixing Ratio for ARPS
!  RH is 0-100% so it is multipled by 0.01
!
!-----------------------------------------------------------------------
!
!     DO 425 klev=1,nz_ext
!     DO 425 j=1,ny_ext
!     DO 425 i=1,nx_ext
!       qvsat = f_qvsat( p_ext(i,j,klev), t_ext(i,j,klev) )
!       qv_ext(i,j,klev)=0.01*qv_ext(i,j,klev)*qvsat
!  425   CONTINUE
  PRINT *, ipr_ext(2),' mixing ratio ',qv_ext(60,36,2)
  PRINT *, ipr_ext(4),' mixing ratio ',qv_ext(60,36,4)
  PRINT *, ipr_ext(6),' mixing ratio ',qv_ext(60,36,6)
  PRINT *, ipr_ext(8),' mixing ratio ',qv_ext(60,36,8)
!
!-----------------------------------------------------------------------
!
!  Set height field
!  by integrating T down from kstart level
!
!-----------------------------------------------------------------------
!
  DO klev=kstart-1,1,-1
    DO j=1,ny_ext
      DO i=1,nx_ext
        dln=ALOG(p_ext(i,j,klev)/p_ext(i,j,klev+1))
        const=dln*rddg
        tbar=0.5*(t_ext(i,j,klev)+t_ext(i,j,klev+1))
        hgt_ext(i,j,klev)=hgt_ext(i,j,klev+1)-const*tbar
      END DO
    END DO
    PRINT *, ' pr,hgt(ks),hgt = ',ipr_ext(klev),                        &
               hgt_ext(60,36,kstart),hgt_ext(60,36,klev)
  END DO
!
!-----------------------------------------------------------------------
!
!  Rotate winds to be relative to true north.
!  The RUC data are sent as grid-relative.
!
!-----------------------------------------------------------------------
!
  DO klev=1,nz_ext
!2001-05-16 GMB: Having umap & uear (or vmap & vear) point to
!the same array causes numerical errors when optimizing.
    CALL uvmptoe(nx_ext,ny_ext,u_ext(1,1,klev),v_ext(1,1,klev),         &
                 lon_ext,utmp,vtmp)
    u_ext(:,:,klev) = utmp(:,:)
    v_ext(:,:,klev) = vtmp(:,:)
  END DO
!
!-----------------------------------------------------------------------
!
!  Fill qc,qr,qi,qs,qh arrays with missing value.
!
!-----------------------------------------------------------------------
!
!  DO k=1,nz_ext
!    DO j=1,ny_ext
!      DO i=1,nx_ext
!        qc_ext(i,j,k)=-999.
!        qr_ext(i,j,k)=-999.
!        qi_ext(i,j,k)=-999.
!        qs_ext(i,j,k)=-999.
!        qh_ext(i,j,k)=-999.
!      END DO
!    END DO
!  END DO
!
!-----------------------------------------------------------------------
!
!  Reset map projection to previous values
!
!-----------------------------------------------------------------------
!
  CALL setmapr(iproj,scale,latnot,trlon)
  CALL setorig(1,x0,y0)
!
!-----------------------------------------------------------------------
!
!  Set good status
!
!-----------------------------------------------------------------------
!
  istatus=1

  DEALLOCATE(utmp)
  DEALLOCATE(vtmp)

  RETURN
END SUBROUTINE getgemeta2
