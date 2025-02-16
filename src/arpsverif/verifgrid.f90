!########################################################################
!########################################################################
!######                                                            ######
!######                   PROGRAM VERIFGRID                        ######
!######                                                            ######
!######                      Developed by                          ######
!######         Center for Analysis and Prediction of Storms       ######
!######                   University of Oklahoma                   ######
!######                                                            ######
!########################################################################
!########################################################################

PROGRAM verifgrid

!----------------------------------------------------------------------- 
!
! PURPOSE:
!
! Reads in two verification files, performs horizontal grid 
! interpolation if desired by the user, and outputs bias and RMSE
! scores.
!
! AUTHOR:  Eric Kemp, November 1999
!
! MODIFICATION HISTORY
! Eric Kemp, 31 March 2000
! Added lat/lon coordinates of corners for NCL plotting.  Also replaced
! some modules with FORTRAN 77 include files, modified character
! variable declarations, and corrected bug in wind differences 
! calculations.
!
! Eric Kemp, 7 July 2000
! Added allocation statement for arrays iloc and jloc when grids match.
!
!-----------------------------------------------------------------------
!
! Use modules
!
!-----------------------------------------------------------------------

! USE extdims2
  USE verif
  
!-----------------------------------------------------------------------
!
! Variable declarations
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INCLUDE 'globcst.inc'

  INTEGER :: fdate(6),adate(6)
  INTEGER,ALLOCATABLE :: ftimesec(:),atimesec(:),timesec(:)

  INTEGER :: fmapproj,amapproj
  REAL :: fscale,ascale
  REAL :: ftrulat(2),ftrulon,atrulat(2),atrulon
  REAL :: fdx,fdy,adx,ady
  REAL :: fscswlat,fscswlon,ascswlat,ascswlon
  REAL :: fmapswlat,fmapswlon,amapswlat,amapswlon
  REAL :: fmapnwlat,fmapnwlon,amapnwlat,amapnwlon
  REAL :: fmapselat,fmapselon,amapselat,amapselon
  REAL :: fmapnelat,fmapnelon,amapnelat,amapnelon

  REAL :: corner_lat(2,2), corner_lon(2,2)
  REAL :: imapswlat,imapswlon,imapnwlat,imapnwlon, &
          imapselat,imapselon,imapnelat,imapnelon

  CHARACTER*4 :: fmodel,fgrid,amodel,agrid, fid_p, aid_p, fid_sfc, aid_sfc
  CHARACTER*16 :: fname_p, fname_sfc, aname_p, aname_sfc
  CHARACTER*3 :: funit_p, funit_sfc, aunit_p, aunit_sfc

  LOGICAL :: comcoord

  REAL,ALLOCATABLE :: fx_ext(:),ax_ext(:),fxs(:),axs(:)
  REAL,ALLOCATABLE :: fy_ext(:),ay_ext(:),fys(:),ays(:)
  REAL :: fx0,ax0,fy0,ay0
  REAL,ALLOCATABLE :: flat_ext(:,:),alat_ext(:,:)
  REAL,ALLOCATABLE :: flon_ext(:,:),alon_ext(:,:)

  REAL,ALLOCATABLE :: dxfld(:),dyfld(:),rdxfld(:),rdyfld(:)

  REAL,ALLOCATABLE,DIMENSION(:,:) :: slopey, alphay, betay

  REAL,ALLOCATABLE :: ftem(:,:),item(:,:),atem(:,:),dtem(:,:)
  REAL :: avgitem,avgatem
  INTEGER :: itemcount,atemcount

  INTEGER,ALLOCATABLE :: iloc(:,:),jloc(:,:)

  REAL,ALLOCATABLE :: fx2d(:,:),fy2d(:,:)
  REAL,ALLOCATABLE :: ax2d(:,:),ay2d(:,:)

  REAL,ALLOCATABLE, DIMENSION(:,:,:,:) :: &
      fvar_sfc, ivar_sfc, avar_sfc, diff_sfc
  REAL,ALLOCATABLE, DIMENSION(:,:,:,:,:) :: &
      fvar_p, ivar_p, avar_p, diff_p

  REAL,ALLOCATABLE :: fvar2d(:,:)

  REAL :: vtrulat(2),vscale
  REAL,ALLOCATABLE :: vxs(:),vys(:)
  REAL,ALLOCATABLE :: vx(:),vy(:)
  REAL,ALLOCATABLE :: vx2d(:,:),vy2d(:,:)

!-----------------------------------------------------------------------
!
! Functions
!
!-----------------------------------------------------------------------

  REAL,EXTERNAL :: pntint2d2
!fpp$ expand (pntint2d2)
!dir$ inline always pntint2d2

!-----------------------------------------------------------------------
!
! Miscellaneous variables
!
!-----------------------------------------------------------------------

  REAL,PARAMETER :: missing = -9999.
  INTEGER,PARAMETER :: imissing = -9999

  INTEGER :: i,j,k,l,ii,jj,kk,ll,m,n
  INTEGER :: fn_datasets,an_datasets
  INTEGER :: maxntime

  INTEGER :: imin,imax,jmin,jmax
  INTEGER :: aimin,aimax,ajmin,ajmax

  REAL :: maxy2d,miny2d,maxx2d,minx2d
  REAL :: maxay2d,minay2d,maxax2d,minax2d

  REAL :: vxctr,vyctr

  REAL :: vxmin,vxmax,vymin,vymax

  INTEGER,ALLOCATABLE :: counter_p(:,:,:), counter_sfc(:,:)

  REAL,ALLOCATABLE, DIMENSION(:,:,:) :: bias_p, rms_p
  REAL,ALLOCATABLE, DIMENSION(:,:) :: bias_sfc, rms_sfc

  REAL :: vctrx,vctry,vswx,vswy

  INTEGER :: status

  INTEGER :: fistart,fiend,fjstart,fjend
  INTEGER :: fimid,fjmid
  REAL :: fxmid,fymid

  REAL :: axpnt,aypnt
  
  INTEGER :: mastercounter
  CHARACTER*72 :: line72= &
    '------------------------------------------------------------------------'

  INTEGER :: fflag_p,fflag_sfc
  INTEGER :: aflag_p,aflag_sfc
  INTEGER :: vflag_p,vflag_sfc

!----------------------------------------------------------------------- 
!
! Namelists
!
!-----------------------------------------------------------------------

  INTEGER :: fnx,fny,fntime, fextdopt
  CHARACTER*256 :: ffilename
  NAMELIST /forecast/ fnx,fny,fntime,fextdopt,ffilename 

  INTEGER :: anx,any,antime, aextdopt
  CHARACTER*256 :: afilename
  NAMELIST /analysis/ anx,any,antime,aextdopt,afilename

  INTEGER :: korder
  NAMELIST /interp/ korder

  INTEGER :: vmapproj, vnx,vny, vbuffer, vimin,vimax,vjmin,vjmax
  REAL :: vdx,vdy,vctrlat,vctrlon, vtrulat1,vtrulat2,vtrulon,vsclfct
  NAMELIST /verification/ vdx,vdy,vctrlat,vctrlon,vmapproj,vtrulat1,vtrulat2,  &
                   vtrulon,vsclfct,vnx,vny,vbuffer,vimin,vimax,vjmin,   &
                   vjmax

  CHARACTER*256 :: statsfile,hfilename,hdiffile
  INTEGER :: if_diff=0
  NAMELIST /stats/ statsfile,hfilename,hdiffile, if_diff

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  fmodel=''
  fgrid=''
  amodel=''
  agrid=''

  WRITE(6,'(//5x,a)')                                                   &
 &'###################################################################'
  WRITE(6,'(5x,a,/5x,a)')                                               &
 &'#                                                                 #',&
 &'# Welcome to VERIFGRID, a program that reads in verification      #'
  WRITE(6,'(5x,a)')                                                     &
 &'# HDF files from CVT2VERIF and calculates statistics in text      #'
  WRITE(6,'(5x,a)')                                                     &
 &'# and HDF format.                                                 #',&
 &'#                                                                 #',&
 &'###################################################################'
  WRITE(6,*)

!-----------------------------------------------------------------------
!
! Read in namelists
!
!-----------------------------------------------------------------------

  PRINT *, 'Reading NAMELIST forecast'
  READ(5,forecast)
  PRINT *, 'Reading NAMELIST analysis'
  READ(5,analysis)
  PRINT *, 'Reading NAMELIST interp'
  READ(5,interp)
  PRINT *, 'Reading NAMELIST verification'
  READ(5,verification)
  PRINT *, 'Reading NAMELIST stats'
  READ(5,stats)

!-----------------------------------------------------------------------
!
! Read forecast data.
!
!-----------------------------------------------------------------------

  ALLOCATE (ftimesec(fntime), &
            fvar_p(fnx,fny,nlevel,fntime,nvar_sfc), &
            fvar_sfc(fnx,fny,fntime,nvar_sfc), &
            STAT=status)

  IF (status /= 0) CALL alloc_fail (status, 'f1')

  CALL rdverif ( fnx,fny,nlevel,fntime,nvar_p,nvar_sfc,		&
                 ffilename,fmodel,fgrid,fdate,ftimesec,pressure,	&
                 fmapproj, fscale, ftrulat(1),ftrulat(2),ftrulon, &
                 fdx,fdy, fscswlat,fscswlon, &
                 fmapswlat,fmapswlon,fmapnwlat,fmapnwlon, &  
                 fmapselat,fmapselon,fmapnelat,fmapnelon, &  
                 fflag_p,fvar_p, fid_p, fname_p, funit_p, &
		 fflag_sfc,fvar_sfc, fid_sfc, fname_sfc, funit_sfc )

  CALL setmapr(fmapproj,fscale,ftrulat,ftrulon)
  CALL lltoxy(1,1,fscswlat,fscswlon,fx0,fy0)
  CALL setorig(1,fx0,fy0)
 
  ALLOCATE (fx_ext(fnx), fxs(fnx), fy_ext(fny), fys(fny), &
            flat_ext(fnx,fny), flon_ext(fnx,fny), STAT=status)

  IF (status /= 0) CALL alloc_fail (status, 'f2')

  DO i = 1,fnx
    fx_ext(i) = (i-1)*fdx
  END DO

  DO j = 1,fny
    fy_ext(j) = (j-1)*fdy
  END DO

  fxs = fx_ext
  fys = fy_ext

  CALL xytoll(fnx,fny,fxs,fys,flat_ext,flon_ext)

!-----------------------------------------------------------------------
!
! Read analysis data
!
!-----------------------------------------------------------------------

  ALLOCATE (atimesec(antime),&
            avar_p(anx,any,nlevel,antime,nvar_sfc), &
            avar_sfc(anx,any,antime,nvar_sfc), &
            STAT=status)

  IF (status /= 0) CALL alloc_fail (status, 'a1')
  avar_sfc(:,:,:,:) = missing  

  CALL rdverif ( anx,any,nlevel,antime,nvar_p,nvar_sfc,		&
                 afilename,amodel,agrid,adate,atimesec,pressure,	&
                 amapproj, ascale, atrulat(1),atrulat(2),atrulon, &
                 adx,ady,ascswlat,ascswlon, & 
                 amapswlat,amapswlon,amapnwlat,amapnwlon, &  
                 amapselat,amapselon,amapnelat,amapnelon, &  
                 aflag_p,avar_p, aid_p, aname_p, aunit_p, &
		 aflag_sfc,avar_sfc, aid_sfc, aname_sfc, aunit_sfc )

  CALL setmapr(amapproj,ascale,atrulat,atrulon)
  CALL lltoxy(1,1,ascswlat,ascswlon,ax0,ay0)
  CALL setorig(1,ax0,ay0)
 
  ALLOCATE (ax_ext(anx), axs(anx), ay_ext(any), ays(any), &
            alat_ext(anx,any), alon_ext(anx,any), STAT=status)
  IF (status /= 0) CALL alloc_fail (status, 'a2')

  DO i = 1,anx
    ax_ext(i) = (i-1)*adx
  END DO

  DO j = 1,any
    ay_ext(j) = (j-1)*ady
  END DO

  axs = ax_ext
  ays = ay_ext

  CALL xytoll(anx,any,axs,ays,alat_ext,alon_ext)

!-----------------------------------------------------------------------
!
! Determine if start dates for forecast and analysis datasets match.
!
!-----------------------------------------------------------------------

  DO i=1,6
    IF (fdate(i) /= adate(i)) THEN
      PRINT *, 'FATAL: Forecast and analysis start dates do not match'
      PRINT *, 'fdate= ', fdate
      PRINT *, 'adate= ', adate
    STOP
    END IF
  END DO

!-----------------------------------------------------------------------
!
! Determine if grids match.
!
!-----------------------------------------------------------------------

  PRINT *, 'Checking to see if grids match'
  ALLOCATE (fx2d(fnx,fny), fy2d(fnx,fny), &
            ax2d(anx,any), ay2d(anx,any), STAT=status)
  IF (status /= 0) CALL alloc_fail (status, 'fa')

  IF (fmapproj == amapproj .AND. fscale == ascale .AND.                 &
      ftrulat(1) == atrulat(1) .AND. ftrulat(2) == atrulat(2) .AND.     &
      ftrulon == atrulon .AND. fdx == adx .AND. fdy == ady .AND.        &
      fnx == anx .AND. fny == any .AND. fscswlat == ascswlat .AND.      &
      fscswlon == ascswlon ) THEN
    comcoord = .TRUE.
    WRITE(6,*)' Grids share a common coordinate system'
    
    DO j = 1,any
      ax2d(:,j) = axs(:)
    END DO
    DO i = 1,anx
      ay2d(i,:) = ays(:)
    END DO

  ELSE
    comcoord = .FALSE.
    WRITE(6,*)'Grid coordinate systems differ, converting lat/lon'

    CALL setmapr(fmapproj,fscale,ftrulat,ftrulon)
    CALL setorig(1,fx0,fy0)
    CALL lltoxy(anx,any,alat_ext,alon_ext,ax2d,ay2d)
    CALL lltoxy(fnx,fny,flat_ext,flon_ext,fx2d,fy2d)
  ENDIF

!-----------------------------------------------------------------------
!
! If grids are different, horizontally interpolate the forecast
! to the analysis.
!
!-----------------------------------------------------------------------

  ALLOCATE (ivar_p(anx,any,nlevel,fntime,nvar_p), &
            ivar_sfc(anx,any,fntime,nvar_sfc), STAT=status)
  IF (status /= 0) CALL alloc_fail (status, 'ivar')
  ivar_p(:,:,:,:,:) = missing
  ivar_sfc(:,:,:,:) = missing

  CALL FLUSH (6)

  IF (.NOT.comcoord) THEN

    ALLOCATE (dxfld(fnx), dyfld(fny), rdxfld(fnx), rdyfld(fny), &
              slopey(fnx,fny), alphay(fnx,fny), betay(fnx,fny), &
              ftem(fnx,fny), item(anx,any), atem(anx,any), dtem(anx,any), &
              iloc(anx,any), jloc(anx,any), fvar2d(fnx,fny), STAT=status)
    IF (status /= 0) CALL alloc_fail (status, 'x')

    item(:,:) = missing
    atem(:,:) = missing
    dtem(:,:) = missing 

    iloc(:,:) = imissing
    jloc(:,:) = imissing
 
    PRINT *, 'Calling setdxdy'
    CALL setdxdy(fnx,fny,1,fnx-1,1,fny-1,fxs,fxs,dxfld,dyfld,            &
                 rdxfld,rdyfld)
    PRINT *, 'Calling setijloc2'
    CALL setijloc2(fnx,fny,anx,any,iloc,jloc,fx2d,fy2d,ax2d,ay2d)

    PRINT *, 'Interpolating sfc vars'
    DO l = 1,nvar_sfc
      DO k = 1,fntime
        CALL setdrvy2(fnx,fny,1,2,fnx-1,2,fny-1,1,1,dyfld,               &
                      rdyfld,fvar_sfc(1,1,k,l),slopey,alphay,betay)
        DO n = 2,any-1
          DO m = 2,anx-1
            IF (iloc(m,n) /= imissing .AND. jloc(m,n) /= imissing) THEN
              i = iloc(m,n)
              j = jloc(m,n)
              IF (fvar_sfc(i,j,k,l) /= missing .AND.                   &
                  fvar_sfc(i-1,j,k,l) /= missing .AND.                 &
                  fvar_sfc(i+1,j,k,l) /= missing .AND.                 &
                  fvar_sfc(i,j-1,k,l) /= missing .AND.                 &
                  fvar_sfc(i,j+1,k,l) /= missing) THEN

                ivar_sfc(m,n,k,l) = pntint2d2(fnx,fny,3,fnx-3,3,fny-2, &
                                           korder,fxs,fys,ax2d(m,n), &
                                           ay2d(m,n),i,j,            &
                                           fvar_sfc(1,1,k,l),          &
                                           dxfld,dyfld,rdxfld,rdyfld,&
					   slopey,alphay,betay)
              ENDIF
            ENDIF
          END DO
        END DO
      END DO
    END DO

    PRINT *, 'Interpolating pressure vars'
    DO l = 1,nvar_p
      DO k = 1,fntime
	DO kk = 1,nlevel
	  CALL setdrvy2(fnx,fny,1,2,fnx-1,2,fny-1,1,1,dyfld,               &
			rdyfld,fvar_p(1,1,kk,k,l),slopey,alphay,betay)
	  DO n = 2,any-1
	    DO m = 2,anx-1
	      IF (iloc(m,n) /= imissing .AND. jloc(m,n) /= imissing) THEN
		i = iloc(m,n)
		j = jloc(m,n)
		IF (fvar_p(i,j,kk,k,l) /= missing .AND.                   &
		    fvar_p(i-1,j,kk,k,l) /= missing .AND.                 &
		    fvar_p(i+1,j,kk,k,l) /= missing .AND.                 &
		    fvar_p(i,j-1,kk,k,l) /= missing .AND.                 &
		    fvar_p(i,j+1,kk,k,l) /= missing) THEN

		  ivar_p(m,n,kk,k,l) = pntint2d2(fnx,fny,3,fnx-3,3,fny-3, &
					     korder,fxs,fys,ax2d(m,n), &
					     ay2d(m,n),i,j,            &
					     fvar_p(1,1,kk,k,l),        &
					     dxfld,dyfld,rdxfld,rdyfld,&
					     slopey,alphay,betay)
		END IF
	      END IF
	    END DO	!m
	  END DO	!n
        END DO	!kk
      END DO	!k
    END DO	!l

  ELSE

    ALLOCATE (iloc(anx,any), jloc(anx,any), STAT=status) ! EMK 7/7/00
    ivar_p = fvar_p
    ivar_sfc = fvar_sfc

  ENDIF

!-----------------------------------------------------------------------
!
! Set limits of verification domain.
!
!-----------------------------------------------------------------------

  IF (vbuffer >= 0) THEN
    vimin = 1+vbuffer
    vjmin = 1+vbuffer
    vimax = vnx-vbuffer
    vjmax = vny-vbuffer
  ENDIF

  IF (vimin < 1) THEN    
    vimin = 1
    WRITE(6,*)'ERROR:  vimin < 1.  Resetting vimin to ',vimin
  ENDIF
  IF (vjmin < 1) THEN    
    vjmin = 1
    WRITE(6,*)'ERROR:  vjmin < 1.  Resetting vjmin to ',vjmin
  ENDIF
  IF (vimax > vnx) THEN    
    vimax = vnx
    WRITE(6,*)'ERROR:  vimax > vnx.  Resetting vimax to ',vimax
  ENDIF
  IF (vjmax > vny) THEN    
    vjmax = vny
    WRITE(6,*)'ERROR:  vjmax > vny.  Resetting vjmax to ',vjmax
  ENDIF

  vtrulat(1) = vtrulat1
  vtrulat(2) = vtrulat2
  vscale = vsclfct

  ALLOCATE (vxs(vnx),vys(vny), vx(vnx),vy(vny), vx2d(vnx,vny),vy2d(vnx,vny), &
            STAT=status)
  IF (status /= 0) CALL alloc_fail (status, 'z')

  iloc(:,:) = imissing
  jloc(:,:) = imissing

  CALL setmapr(vmapproj,vscale,vtrulat,vtrulon)
  CALL lltoxy(1,1,vctrlat,vctrlon,vctrx,vctry)
  vswx = vctrx - (FLOAT(vnx-3)/2.) * (vdx*vsclfct)
  vswy = vctry - (FLOAT(vny-3)/2.) * (vdy*vsclfct)
  CALL setorig (1,vswx,vswy)

  DO i = 1,vnx
    vx(i) = (i-2)*vdx
  END DO
  DO j = 1,vny
    vy(j) = (j-2)*vdy
  END DO

  DO i = 1,vnx-1
    vxs(i) = 0.5*(vx(i)+vx(i+1))
  END DO
  vxs(vnx) = (2.*vx(vnx))-vx(vnx-1)

  DO j = 1,vny-1
    vys(j) = 0.5*(vy(j)+vy(j+1))
  END DO
  vys(vny) = (2.*vy(vny))-vy(vny-1)
 
  DO j = 1,vny
    vx2d(:,j) = vxs(:)
  END DO
  DO i = 1,vnx
    vy2d(i,:) = vys(:)
  END DO


  CALL lltoxy(anx,any,alat_ext,alon_ext,ax2d,ay2d)
  CALL setijloc2(vnx,vny,anx,any,iloc,jloc,vx2d,vy2d,ax2d,ay2d)

  CALL xytoll ( 1, 1, vx2d(vimin,vjmin), vy2d(vimin,vjmin), &
                corner_lat(1,1), corner_lon(1,1) )
  CALL xytoll ( 1, 1, vx2d(vimin,vjmax), vy2d(vimin,vjmax), &
                corner_lat(1,2), corner_lon(1,2) )
  CALL xytoll ( 1, 1, vx2d(vimax,vjmin), vy2d(vimax,vjmin), &
                corner_lat(2,1), corner_lon(2,1) )
  CALL xytoll ( 1, 1, vx2d(vimax,vjmax), vy2d(vimax,vjmax), &
                corner_lat(2,2), corner_lon(2,2) )

!-----------------------------------------------------------------------
!
! Calculate differences for each variable
!
!-----------------------------------------------------------------------

  maxntime = MIN(fntime,antime)
  WRITE(6,*)'maxntime = ',maxntime

  ALLOCATE (timesec(maxntime), &
      diff_p(anx,any,nlevel,maxntime,nvar_sfc_stats), &
      diff_sfc(anx,any,maxntime,nvar_sfc_stats), &
      STAT=status)
  IF (status /= 0) CALL alloc_fail (status, 'd')

  diff_sfc(:,:,:,:) = missing
  diff_p(:,:,:,:,:) = missing

  m = 0

  DO k = 1,fntime
    DO kk = 1,antime
      IF (ftimesec(k) == atimesec(kk)) THEN
        m = m + 1
        timesec(m) = ftimesec(k)
      END IF
    END DO
  END DO
 
  maxntime = m
  PRINT *, 'Times in fcst file  (', m, '):', ftimesec(1:fntime)
  PRINT *, 'Times in anl  file  (', m, '):', atimesec(1:antime)
  PRINT *, 'Times in both files (', m, '):', timesec(1:m)

  PRINT *, 'Computing differences of pressure vars'
  DO l = 1,nvar_p
    DO k = 1, fntime
      DO kk = 1,antime
	DO m = 1,maxntime
	  IF (timesec(m) == ftimesec(k) .AND.                   &
	      timesec(m) == atimesec(kk)) THEN
	    mastercounter = 0
	    DO ll = 1,nlevel
	      DO j = 1,any-1
		DO i = 1,anx-1
		  IF (iloc(i,j) /= imissing .AND. iloc(i,j) >= vimin     &
		      .AND. iloc(i,j) <= vimax .AND.                  &
		      jloc(i,j) /= imissing .AND. jloc(i,j) >= vjmin     &
		      .AND. jloc(i,j) <= vjmax) THEN
		    mastercounter = mastercounter + 1

		    IF (ivar_p(i,j,ll,k,l) /= missing .AND.              &
			avar_p(i,j,ll,kk,l) /= missing) THEN
		      diff_p(i,j,ll,m,l) = ivar_p(i,j,ll,k,l) - &
		                              avar_p(i,j,ll,kk,l)
		    ENDIF

		  ENDIF
		END DO	!i
	      END DO	!j
	    END DO	!ll
	  ENDIF
	END DO	!m
      END DO	!kk
    END DO	!k
  END DO	!l

  PRINT *, 'Computing differences of sfc vars'
  DO l = 1,nvar_sfc
    DO k = 1, fntime
      DO kk = 1,antime
	DO m = 1,maxntime
	  IF (timesec(m) == ftimesec(k) .AND.                   &
	      timesec(m) == atimesec(kk)) THEN
	    mastercounter = 0

	    DO j = 2,any-1
	      DO i = 2,anx-1
		IF (iloc(i,j) /= imissing .AND. iloc(i,j) >= vimin     &
		    .AND. iloc(i,j) <= vimax .AND.                  &
		    jloc(i,j) /= imissing .AND. jloc(i,j) >= vjmin     &
		    .AND. jloc(i,j) <= vjmax) THEN
		  mastercounter = mastercounter + 1

		  IF (ivar_sfc(i,j,k,l) /= missing .AND.              &
		      avar_sfc(i,j,kk,l) /= missing) THEN
	    diff_sfc(i,j,m,l) = ivar_sfc(i,j,k,l) - avar_sfc(i,j,kk,l)
		  ENDIF

		ENDIF
	      END DO	!i
	    END DO	!j

	  ENDIF
	END DO	!m
      END DO	!kk
    END DO	!k
  END DO	!l

!-----------------------------------------------------------------------
!
! Calculate differences for total winds
!
!-----------------------------------------------------------------------

  PRINT *, 'Computing differences of total winds at sfc'
  n = 0		! index of sfc
  DO k = 1,fntime
    DO kk = 1,antime
      DO m = 1,maxntime
	IF (timesec(m) == ftimesec(k) .AND.  &
	    timesec(m) == atimesec(kk)) THEN
	  DO j = 1,any-1  
	    DO i = 1,anx-1
	      IF (iloc(i,j) /= imissing .AND. iloc(i,j) >= vimin .AND. & 
		  iloc(i,j) <= vimax .AND. jloc(i,j) /= imissing .AND. &
		  jloc(i,j) >= vjmin .AND. jloc(i,j) <= vjmax) THEN
		IF (avar_sfc(i,j,kk,id_u) /= missing .AND.           &
		    avar_sfc(i,j,kk,id_v) /= missing .AND.            &
		    ivar_sfc(i,j,k,id_u) /= missing .AND.          &
		    ivar_sfc(i,j,k,id_u) /= missing) THEN
		  diff_sfc(i,j,m,id_wspd) =  &
	   SQRT( ivar_sfc(i,j,k,id_u)**2  + ivar_sfc(i,j,k,id_v)**2 ) - &
	   SQRT( avar_sfc(i,j,kk,id_u)**2 + avar_sfc(i,j,kk,id_v)**2 )
		END IF                       
	      END IF
	    END DO
	  END DO
	END IF
      END DO
    END DO
  END DO

  PRINT *, 'Computing differences of total winds at pressure levels'
  DO k = 1,fntime
    DO kk = 1,antime
      DO m = 1,maxntime
	IF (timesec(m) == ftimesec(k) .AND.               &
	    timesec(m) == atimesec(kk)) THEN
	  DO n=1,nlevel
	    DO j = 1,any-1  
	      DO i = 1,anx-1
		IF (iloc(i,j) /= imissing .AND. iloc(i,j) >= vimin .AND. & 
		    iloc(i,j) <= vimax .AND. jloc(i,j) /= imissing .AND. &
		    jloc(i,j) >= vjmin .AND. jloc(i,j) <= vjmax) THEN
		  IF (avar_p(i,j,n,kk,id_u) /= missing .AND.           &
		      avar_p(i,j,n,kk,id_v) /= missing .AND.            &
		      ivar_p(i,j,n,k,id_u) /= missing .AND.          &
		      ivar_p(i,j,n,k,id_u) /= missing) THEN
		    diff_p(i,j,n,m,id_wspd) =  &
	     SQRT( ivar_p(i,j,n,k,id_u)**2  + ivar_p(i,j,n,k,id_v)**2 ) - &
	     SQRT( avar_p(i,j,n,kk,id_u)**2 + avar_p(i,j,n,kk,id_v)**2 )
		  END IF                       
		END IF
	      END DO
	    END DO
	  END DO
	END IF
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
! Calculate statistics
!
!-----------------------------------------------------------------------

  ALLOCATE ( bias_p(nlevel,maxntime,nvar_p_stats),		&
             rms_p(nlevel,maxntime,nvar_p_stats),		&
             bias_sfc(maxntime,nvar_sfc_stats),			&
	     rms_sfc(maxntime,nvar_sfc_stats),			&
             counter_sfc(maxntime,nvar_sfc_stats), 		&
             counter_p(nlevel,maxntime,nvar_p_stats),		&
             STAT=status)
  IF (status /= 0) CALL alloc_fail (status, 'stats')

  bias_p(:,:,:)   = 0.0
  rms_p(:,:,:)    = 0.0
  bias_sfc(:,:)   = 0.0
  rms_sfc(:,:)    = 0.0
  counter_sfc(:,:)  = 0
  counter_p(:,:,:)= 0
  
  ! 5D vars

  PRINT *, 'Computing statistics on pressure levels'
  DO l = 1,nvar_p
    DO k = 1,maxntime
      DO m = 1,nlevel
	DO j = 1,any
	  DO i = 1,anx
	    
	    IF (diff_p(i,j,m,k,l) /= missing) THEN
	      bias_p(m,k,l) = bias_p(m,k,l) + diff_p(i,j,m,k,l)
	      rms_p(m,k,l) = rms_p(m,k,l) + diff_p(i,j,m,k,l)**2
	      counter_p(m,k,l) = counter_p(m,k,l) + 1
	    ENDIF

	  END DO
	END DO

	IF (counter_p(m,k,l) > 0) THEN
	  bias_p(m,k,l) = bias_p(m,k,l)/counter_p(m,k,l)
	  rms_p(m,k,l) = SQRT( rms_p(m,k,l)/counter_p(m,k,l) )
	ELSE
          bias_p(m,k,l)   = missing
          rms_p(m,k,l)    = missing
	ENDIF

      END DO	!m
    END DO	!k
  END DO	!l

  ! 4D vars

  PRINT *, 'Computing statistics at sfc'
  DO l = 1,nvar_sfc
    DO k = 1,maxntime
      DO j = 1,any
        DO i = 1,anx
          
          IF (diff_sfc(i,j,k,l) /= missing) THEN
            bias_sfc(k,l) = bias_sfc(k,l) + diff_sfc(i,j,k,l)
            rms_sfc(k,l) = rms_sfc(k,l) + diff_sfc(i,j,k,l)**2
            counter_sfc(k,l) = counter_sfc(k,l) + 1
          ENDIF

	  END DO
	END DO

      IF (counter_sfc(k,l) > 0) THEN
        bias_sfc(k,l) = bias_sfc(k,l)/counter_sfc(k,l)
        rms_sfc(k,l) = SQRT( rms_sfc(k,l)/counter_sfc(k,l) )
      ELSE
	bias_sfc(k,l)   = missing
	rms_sfc(k,l)    = missing
      ENDIF

    END DO	!k
  END DO	!l

  PRINT *, 'Computing wind statistics at p'
  DO k = 1,maxntime
    DO m = 1,nlevel
      DO j = 1,any
	DO i = 1,anx
	  IF (diff_p(i,j,m,k,id_wspd) /= missing) THEN
	    bias_p(m,k,id_wspd) = bias_p(m,k,id_wspd) + diff_p(i,j,m,k,id_wspd)
	    rms_p(m,k,id_wspd) = rms_p(m,k,id_wspd) + diff_p(i,j,m,k,id_wspd)**2
	    rms_p(m,k,id_wvec) = rms_p(m,k,id_wvec) + &
			 diff_p(i,j,m,k,id_u)**2 + diff_p(i,j,m,k,id_v)**2
	    counter_p(m,k,id_wspd) = counter_p(m,k,id_wspd) + 1
	  ENDIF
	END DO
      END DO

      counter_p(m,k,id_wvec) = counter_p(m,k,id_wspd)

      IF (counter_p(m,k,id_wspd) > 0) THEN
	bias_p(m,k,id_wspd) = bias_p(m,k,id_wspd)/counter_p(m,k,id_wspd)
	rms_p(m,k,id_wspd) = SQRT( rms_p(m,k,id_wspd)/counter_p(m,k,id_wspd) )
	rms_p(m,k,id_wvec) = SQRT( rms_p(m,k,id_wvec)/counter_p(m,k,id_wvec) )
      ELSE
	bias_p(m,k,id_wspd) = missing
	rms_p(m,k,id_wspd)  = missing
	rms_p(m,k,id_wvec)  = missing
      ENDIF
      
    END DO
  END DO

  PRINT *, 'Computing wind statistics at sfc'
  DO k = 1,maxntime
    DO j = 1,any
      DO i = 1,anx
	IF (diff_sfc(i,j,k,id_wspd) /= missing) THEN
	  bias_sfc(k,id_wspd) = bias_sfc(k,id_wspd) + diff_sfc(i,j,k,id_wspd)
	  rms_sfc(k,id_wspd) = rms_sfc(k,id_wspd) + diff_sfc(i,j,k,id_wspd)**2
	  rms_sfc(k,id_wvec) = rms_sfc(k,id_wvec) + &
		       diff_sfc(i,j,k,id_u)**2 + diff_sfc(i,j,k,id_v)**2
	  counter_sfc(k,id_wspd) = counter_sfc(k,id_wspd) + 1
	ENDIF
      END DO
    END DO

    counter_sfc(k,id_wvec) = counter_sfc(k,id_wspd)

    IF (counter_sfc(k,id_wspd) > 0) THEN
      bias_sfc(k,id_wspd) = bias_sfc(k,id_wspd)/counter_sfc(k,id_wspd)
      rms_sfc(k,id_wspd) = SQRT( rms_sfc(k,id_wspd)/counter_sfc(k,id_wspd) )
      rms_sfc(k,id_wvec) = SQRT( rms_sfc(k,id_wvec)/counter_sfc(k,id_wvec) )
    ELSE
      bias_sfc(k,id_wspd) = missing
      rms_sfc(k,id_wspd)  = missing
      rms_sfc(k,id_wvec)  = missing
    ENDIF
    
  END DO

!-----------------------------------------------------------------------
!
! Output statistics
!
!-----------------------------------------------------------------------

  PRINT *, 'Writing text output to ', TRIM(statsfile)
  OPEN(UNIT=13,FILE=statsfile,STATUS="REPLACE")

  WRITE(13,600) 'Forecast initialized ',fdate(1),'-',fdate(2),         &
             '-',fdate(3),'.',fdate(4),':',fdate(5),':',fdate(6)
  WRITE(13,610) 'Forecast model/grid:', TRIM(fmodel), TRIM(fgrid)
  WRITE(13,625) NINT((fnx+1)*fdx/1000.0), NINT((fny+1)*fdx/1000.0), &
                fdx/1000.0, fnx,fny
  WRITE(13,*)
  WRITE(13,610) 'Analysis model/grid:', TRIM(amodel), TRIM(agrid)
  WRITE(13,625) NINT((anx+1)*adx/1000.0), NINT((any+1)*adx/1000.0), &
                adx/1000.0, anx,any
  WRITE(13,*)
  WRITE(13,*) 'Verification Region:'
  WRITE(13,625) NINT((vnx+1)*vdx/1000.0), NINT((vny+1)*vdx/1000.0), &
                vdx/1000.0, vnx,vny
  WRITE(13,650) 'NW/NE', &
      corner_lat(1,2),corner_lon(1,2),corner_lat(2,2),corner_lon(2,2)
  WRITE(13,650) 'SW/SE', &
      corner_lat(1,1),corner_lon(1,1),corner_lat(2,1),corner_lon(2,1)
  WRITE(13,*) 'Buffer: vimin/max, vjmin/max: ', vimin, vimax, vjmin, vjmax
  WRITE(13,*) 'Number of grid points in verif region: ', mastercounter
  WRITE(13,*)

  DO l = 1,nvar_sfc_stats
    WRITE(13,'(A)') TRIM(line72)
    WRITE(13,810) TRIM(varname_sfc(l)), TRIM(varunit_sfc(l))
    WRITE(13,*) 'Time(hr)  Grid Points      Bias  RMS Error'
    DO k = 1,maxntime
      WRITE(13,700) timesec(k)/3600, counter_sfc(k,l), &
          bias_sfc(k,l), rms_sfc(k,l)
    END DO
  END DO

  DO l = 1,nvar_p_stats
  DO m = 1,nlevel
    WRITE(13,'(A)') TRIM(line72)
    WRITE(13,810) TRIM(varname_p(l)), TRIM(varunit_p(l))
    WRITE(13,*) 'Level = ', pressure(m) / 100.0, ' mb'
    WRITE(13,*) 'Time(hr)  Grid Points      Bias  RMS Error'
    DO k = 1,maxntime
      WRITE(13,700) timesec(k)/3600, counter_p(m,k,l), &
          bias_p(m,k,l), rms_p(m,k,l)
    END DO
  END DO
  END DO

  PRINT *, 'Writing statistics to HDF file ', TRIM(hfilename)

  CALL wrt_verif_stats_diff (hfilename, if_diff, hdiffile, fdate, &
                     amodel,agrid,adx,anx,any, &
                     fmodel,fgrid,fdx, &
                     vdx, vdy, vnx, vny, vmapproj, &
		     vtrulat1, vtrulat2, vtrulon, vsclfct, vctrlat, vctrlon, &
		     corner_lat, corner_lon, &
		     nlevel, maxntime, nvar_p_stats, nvar_sfc_stats, &
		     timesec, pressure, missing, mastercounter, &
		     counter_p, bias_p, rms_p, diff_p, &
                     counter_sfc, bias_sfc, rms_sfc, diff_sfc, &
	             varid_p, varname_p, varunit_p, &
	             varid_sfc, varname_sfc, varunit_sfc)

  WRITE(6,*)'Program verifgrid successfully completed.'

  600 FORMAT (/1x,a21,i4.4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2)
  610 FORMAT (1X,A,1X,A,' / ',A)
  625 FORMAT (' Domain:',I5, ' x', I5, ' km, dx=', F5.1, &
              ' km (',I4,' x',I4,')' )
  650 FORMAT (1X,A,' corner lat/lon:', 2F10.3, 2X, 2F10.3)
  700 FORMAT (1x,i8,2x,i11,2x,f8.2,3x,f8.2)
  800 FORMAT (1x,i8,2x,i11,4x,f8.2,4x,f8.2,10x,f8.2)
  810 FORMAT ( 1X,A," (",A,")" )

END PROGRAM verifgrid
