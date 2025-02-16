Program plt1dtrajc
!
! Program to plotting 1D (time v.s. plotted quantity)
! curves of quantities calculated along trajectories by arpscalctrajc
! program.
!
! To compile and run the program:
!
! zxncarf90/zxncarpgf90/zxncarifort -o plt1dtrajc plt1dtrajc.f90
! plt1dtrajc < plt1dtrajc.input
! The input file can be found in src/input directory in ARPS root directory.
!

  Implicit NONE
  REAL :: dx, dy, dz
  REAL :: xlow, xhigh, ylow, yhigh, zlow, zhigh
  REAL :: xtra
  REAL, allocatable :: xtrajc(:,:),ytrajc(:,:),ztrajc(:,:),ttrajc(:)
  REAL, allocatable :: utrajc(:,:),vtrajc(:,:),wtrajc(:,:)
  REAL, allocatable :: vortx_trajc(:,:),vorty_trajc(:,:),vortz_trajc(:,:), &
                       vorts_trajc(:,:),vortc_trajc(:,:)
  REAL, allocatable :: buoy_trajc(:,:),buoyq_trajc(:,:),frcs_trajc(:,:),pprt_trajc(:,:)
  REAL, allocatable :: upgrad_trajc(:,:),vpgrad_trajc(:,:),wpgrad_trajc(:,:)
  REAL, allocatable :: vorts_gen_trajc(:,:),vortc_gen_trajc(:,:)
  REAL, allocatable :: vortz_stch_trajc(:,:),vortz_tilt_trajc(:,:),vhtrajc(:,:)
  REAL, allocatable :: tem_trajc(:)


  INTEGER :: ntrajcs, npoints, npoints_indata
  INTEGER :: i,j,k, istat, nunit, nstart, nend, ntrajcs_in
  CHARACTER(LEN=256) :: trajc_fn, string
  CHARACTER(LEN=80)  :: runname
  CHARACTER(LEN=256) :: outfilename
  CHARACTER(LEN=512) :: string1
  COMMON /runname_common/runname
  integer trajc_color

  REAL :: tstart, tzero, tend, tstart_plt, tend_plt, tref
  NAMELIST /input/ trajc_fn, tstart_plt, tend_plt, tref
  REAL :: missing_value, eps, xlimit, ylimit
  INTEGER :: first_point, i1, i2, istatus, i_in, npoints_in, lfnkey, zxout_unit
  REAL :: x1,x2,y1,y2

  REAL :: umin,umax,vmin,vmax,wmin,wmax,vsmin,vsmax,vcmin,vcmax,  &
          vxmin,vxmax,vymin,vymax,vzmin,vzmax
  REAL :: xmin,xmax,ysmin,ysmax,zxmin,zxmax,pmin,pmax,fmin,fmax, &
          upgmin,upgmax,vpgmin,vpgmax,wpgmin,wpgmax,vsgmin,vsgmax,vcgmin,vcgmax
  REAL :: trajc_no_bgn,trajc_no_end, xo, yo

  missing_value = -99999.0
  eps = 0.001

  tstart_plt = 0.0
  tend_plt = 0.0
  trajc_fn = 'default_file_name'
  READ(5,input)

  nunit = 15

  OPEN(UNIT=nunit,FILE=trim(trajc_fn),STATUS='old',   &
        FORM='formatted',IOSTAT= istat )

  READ(nunit,'(a)') runname
  READ(nunit,'(6e17.6)') xlow, xhigh, ylow, yhigh, zlow, zhigh
  write(6,'(6e17.6)') xlow, xhigh, ylow, yhigh, zlow, zhigh


  READ(nunit,'(3e17.6)') dx, dy, dz
  write(6,'(3e17.6)') dx, dy, dz

  READ(nunit,'(3e17.6)') tstart, tzero, tend
  READ(nunit,'(i10)') npoints
  READ(nunit,'(i10)') ntrajcs
  READ(nunit,'(a)') string1

  print*,'npoints, ntrajcs=', npoints, ntrajcs

  allocate(xtrajc(npoints,ntrajcs))
  allocate(ytrajc(npoints,ntrajcs))
  allocate(ztrajc(npoints,ntrajcs))
  allocate(ttrajc(npoints))

  allocate(vortx_trajc(npoints,ntrajcs),stat=istatus)
  allocate(vorty_trajc(npoints,ntrajcs),stat=istatus)
  allocate(vortz_trajc(npoints,ntrajcs),stat=istatus)
  allocate(vorts_trajc(npoints,ntrajcs),stat=istatus)
  allocate(vortc_trajc(npoints,ntrajcs),stat=istatus)
  allocate(buoy_trajc(npoints,ntrajcs),stat=istatus)
  allocate(buoyq_trajc(npoints,ntrajcs),stat=istatus)
  allocate(vorts_gen_trajc(npoints,ntrajcs),stat=istatus)
  allocate(vortc_gen_trajc(npoints,ntrajcs),stat=istatus)
  allocate(frcs_trajc(npoints,ntrajcs),stat=istatus)
  allocate(pprt_trajc(npoints,ntrajcs),stat=istatus)
  allocate(upgrad_trajc(npoints,ntrajcs),stat=istatus)
  allocate(vpgrad_trajc(npoints,ntrajcs),stat=istatus)
  allocate(wpgrad_trajc(npoints,ntrajcs),stat=istatus)
  allocate(utrajc(npoints,ntrajcs),stat=istatus)
  allocate(vtrajc(npoints,ntrajcs),stat=istatus)
  allocate(wtrajc(npoints,ntrajcs),stat=istatus)
  allocate(vhtrajc(npoints,ntrajcs),stat=istatus)
  allocate(vortz_tilt_trajc(npoints,ntrajcs),stat=istatus)
  allocate(vortz_stch_trajc(npoints,ntrajcs),stat=istatus)
  allocate(tem_trajc(npoints),stat=istatus)

  print*,'npoints = ', npoints

  DO j = 1, npoints
    read(nunit,*,err=1002,end=1005) ttrajc(j)
    print*,'ttrajc(j), j =', ttrajc(j), j
    DO i = 1,ntrajcs

      read(nunit,'(f8.1,6f10.2,22e14.6)',end=1001,err=1002) ttrajc(j),  &
          xtrajc(j,i),ytrajc(j,i),ztrajc(j,i), &
          utrajc(j,i),vtrajc(j,i),wtrajc(j,i), pprt_trajc(j,i), &
          vortx_trajc(j,i),vorty_trajc(j,i),vortz_trajc(j,i),   &
          vorts_trajc(j,i), vortc_trajc(j,i), &
          buoy_trajc(j,i),buoyq_trajc(j,i),frcs_trajc(j,i),     &
          upgrad_trajc(j,i),vpgrad_trajc(j,i),wpgrad_trajc(j,i),  &
          vorts_gen_trajc(j,i),vortc_gen_trajc(j,i),    &
          vortz_tilt_trajc(j,i),vortz_stch_trajc(j,i)

      write(6,'(i5,f8.1,6f10.2,22e14.6)') i, ttrajc(j)

!     write(6,'(i5,f8.1,6f10.2,22e14.6)') i, ttrajc(j),  &
!         xtrajc(j,i),ytrajc(j,i),ztrajc(j,i), &
!         utrajc(j,i),vtrajc(j,i),wtrajc(j,i), pprt_trajc(j,i), &
!         vortx_trajc(j,i),vorty_trajc(j,i),vortz_trajc(j,i),   &
!         vorts_trajc(j,i), vortc_trajc(j,i), &
!         buoy_trajc(j,i),buoyq_trajc(j,i),frcs_trajc(j,i),     &
!         upgrad_trajc(j,i),vpgrad_trajc(j,i),wpgrad_trajc(j,i),  &
!         vorts_gen_trajc(j,i),vortc_gen_trajc(j,i),    &
!         vortz_tilt_trajc(j,i),vortz_stch_trajc(j,i)

      IF( xtrajc(j,i)==missing_value .or. ytrajc(j,i)==missing_value .or. &
          ztrajc(j,i)==missing_value) then
          utrajc(j,i)=missing_value
          vtrajc(j,i)=missing_value
          wtrajc(j,i)=missing_value
          pprt_trajc (j,i)=missing_value
          vortx_trajc(j,i)=missing_value
          vorty_trajc(j,i)=missing_value
          vortz_trajc(j,i)=missing_value
          vorts_trajc(j,i)=missing_value
          vortc_trajc(j,i)=missing_value
          buoy_trajc (j,i)=missing_value
          buoyq_trajc(j,i)=missing_value
          frcs_trajc (j,i)=missing_value
          upgrad_trajc(j,i)=missing_value
          vpgrad_trajc(j,i)=missing_value
          wpgrad_trajc(j,i)=missing_value
          vorts_gen_trajc(j,i)=missing_value
          vortc_gen_trajc(j,i)=missing_value
          vortz_tilt_trajc(j,i)=missing_value
          vortz_stch_trajc(j,i)=missing_value
      ENDIF
    ENDDO
  ENDDO


    goto 1003
1001 continue
    print*,'Error encountered when reading trajectory file.'
    print*,'Number of points in trajectory reset to ', j-1
    npoints_indata = j-1
    goto 1003

1004 continue
    print*,'Error encountered when reading trajectory file.'
    goto 1003

1002 continue
    print*,'EOF encountered when reading trajectory file.'
    print*,'Number of points in trajectory reset to ', j-1
    npoints_indata = j-1
    goto 1003

1005 continue
    print*,'EOF encountered when reading trajectory file.'
    npoints_indata = j-1

1003 continue

  CLOSE(UNIT=nunit)

  if(.false.) then

  OPEN(UNIT=nunit,FILE=trim(trajc_fn)//'_reordered',STATUS='unknown',   &
        FORM='formatted',IOSTAT= istat )

  write(nunit,'(a)') trim(string1)
  write(nunit,'(2i10)') ntrajcs, npoints_indata
  DO i = 1,ntrajcs
    write(nunit,'(i10)') i
    DO j = 1, npoints_indata
      write(nunit,'(f8.1,6f10.2,22e14.6)') ttrajc(j),  &
          xtrajc(j,i),ytrajc(j,i),ztrajc(j,i), &
          utrajc(j,i),vtrajc(j,i),wtrajc(j,i), pprt_trajc(j,i), &
          vortx_trajc(j,i),vorty_trajc(j,i),vortz_trajc(j,i),   &
          vorts_trajc(j,i), vortc_trajc(j,i), &
          buoy_trajc(j,i),buoyq_trajc(j,i),frcs_trajc(j,i),     &
          upgrad_trajc(j,i),vpgrad_trajc(j,i),wpgrad_trajc(j,i),  &
          vorts_gen_trajc(j,i),vortc_gen_trajc(j,i),    &
          vortz_tilt_trajc(j,i),vortz_stch_trajc(j,i)
    ENDDO
  ENDDO

  CLOSE(UNIT=nunit)

  endif


  IF( tstart_plt /= tend_plt ) then
    nstart = npoints_indata
    DO j=1,npoints_indata
      IF( ttrajc(j) >= tstart_plt ) then
        nstart = j
        exit
      ENDIF
    ENDDO

    nend   = 1
    DO j=npoints_indata,1,-1
      IF( ttrajc(j) <= tend_plt ) then
        nend = j
        exit
      ENDIF
    ENDDO
  ELSE
    nstart = 1
    nend   = npoints_indata
  ENDIF

  DO j=1,ntrajcs
    do i=1,npoints_indata
      if( utrajc(i,j) == missing_value .or. vtrajc(i,j) == missing_value ) then
        vhtrajc(i,j)= missing_value
      else
        vhtrajc(i,j)=sqrt( utrajc(i,j)**2 + vtrajc(i,j)**2 )
      endif
    enddo
  ENDDO

  DO i=1,2
  DO j=1,ntrajcs
    CALL smooth3pmv(utrajc(1,j),npoints,nstart,min(nend,npoints_indata-1), missing_value, tem_trajc )
    CALL smooth3pmv(vtrajc(1,j),npoints,nstart,min(nend,npoints_indata-1), missing_value, tem_trajc )
    CALL smooth3pmv(wtrajc(1,j),npoints,nstart,min(nend,npoints_indata-1), missing_value, tem_trajc )
    CALL smooth3pmv(vhtrajc(1,j),npoints,nstart,min(nend,npoints_indata-1), missing_value, tem_trajc )
    CALL smooth3pmv(pprt_trajc(1,j),npoints,nstart,min(nend,npoints_indata-1), missing_value, tem_trajc )
    CALL smooth3pmv(vortx_trajc(1,j),npoints,nstart,min(nend,npoints_indata-1), missing_value, tem_trajc )
    CALL smooth3pmv(vorty_trajc(1,j),npoints,nstart,min(nend,npoints_indata-1), missing_value, tem_trajc )
    CALL smooth3pmv(vortz_trajc(1,j),npoints,nstart,min(nend,npoints_indata-1), missing_value, tem_trajc )
    CALL smooth3pmv(vorts_trajc(1,j),npoints,nstart,min(nend,npoints_indata-1), missing_value, tem_trajc )
    CALL smooth3pmv(vortc_trajc(1,j),npoints,nstart,min(nend,npoints_indata-1), missing_value, tem_trajc )
    CALL smooth3pmv( buoy_trajc(1,j),npoints,nstart,min(nend,npoints_indata-1), missing_value, tem_trajc )
    CALL smooth3pmv(buoyq_trajc(1,j),npoints,nstart,min(nend,npoints_indata-1), missing_value, tem_trajc )
    CALL smooth3pmv( frcs_trajc(1,j),npoints,nstart,min(nend,npoints_indata-1), missing_value, tem_trajc )
    CALL smooth3pmv(upgrad_trajc(1,j),npoints,nstart,min(nend,npoints_indata-1), missing_value, tem_trajc )
    CALL smooth3pmv(vpgrad_trajc(1,j),npoints,nstart,min(nend,npoints_indata-1), missing_value, tem_trajc )
    CALL smooth3pmv(wpgrad_trajc(1,j),npoints,nstart,min(nend,npoints_indata-1), missing_value, tem_trajc )
    CALL smooth3pmv(vortz_tilt_trajc(1,j),npoints,nstart,min(nend,npoints_indata-1), missing_value, tem_trajc )
    CALL smooth3pmv(vortz_stch_trajc(1,j),npoints,nstart,min(nend,npoints_indata-1), missing_value, tem_trajc )
    CALL smooth3pmv(vorts_gen_trajc(1,j),npoints,nstart,min(nend,npoints_indata-1), missing_value, tem_trajc )
    CALL smooth3pmv(vortc_gen_trajc(1,j),npoints,nstart,min(nend,npoints_indata-1), missing_value, tem_trajc )
  ENDDO
  ENDDO

  print*,'tstart_plt, tend_plt, nstart, nend =', tstart_plt, tend_plt, nstart, nend

! tref = 13250.0

! tref = tref/60.0
! tstart_plt =  tstart_plt/60.0
! tend_plt   = tend_plt /60.0
! ttrajc     = ttrajc/60.0

  CALL gtlfnkey ( runname, lfnkey )

  zxout_unit=2
  CALL xpsfn(runname(1:lfnkey)//'.ps', zxout_unit)

  call xdevic(1,outfilename,0,0)
  CALL xdspac(0.95)
  call xsetclrs(3)
  CALL xspace( 1, 6 , 0.0 , xlimit,ylimit)
  CALL xpmagn(0.05, 0.02)

  xlow = 0.0
  xhigh = 48.0
  ylow = 0.0
  yhigh = 48.0
  zlow  = 0.0
  zhigh = 14.0001

!-----------------------------------------------------------------------
! Loop over trajectories
!-----------------------------------------------------------------------

  DO j=1,ntrajcs

  CALL xdspac(0.95)

  CALL xspace( 1, 4 , 0.0 , xlimit,ylimit)
  CALL xpmagn(0.05, 0.03)

!-----------------------------------------------------------------------
! Plot x,y
!-----------------------------------------------------------------------

  call xcolor(1)
  call xnwpic
  call xmap(tstart_plt, tend_plt, xlow, xhigh)
  call xqmap(x1,x2,y1,y2)

  call xaxnsz( (y2-y1)*0.05 )
  xo = 12600.0
  yo = 0.0
  call xaxfmt( '(I5)')
  call XAXTIK(1,1)
  CALL XAXANT(-1,-1)
  call XAXISX1(XO,y1, 60.0,600.0)
  call XAXISY1(x1,YO, 1.0,4.0)
  CALL xbox(x1,x2,y1,y2)

  call xpenup(x1,0.0)
  call xpendn(x2,0.0)
  call xpenup(tref,y1)
  call xpendn(tref,y2)
  call xpenup(tref,y1)
  call xpendn(tref,y2)

  call xchsiz( (y2-y1)*0.06 )
  call xchori(0.0)
  call xcharc( (x1+x2)*0.5, y1-0.15*(y2-y1), 'Time (min)')
  call xchori(90.0)
  call xcharc( x1-(x2-x1)*0.05, (y1+y2)*0.5, 'X and Y (m)')
  call xchori(0.0)
  call xwindw(x1,x2,y1,y2)

  print*,'tstart_plt, tend_plt, zlow, zhigh=', tstart_plt, tend_plt, zlow, zhigh

  call xcolor( mod( j+2, 253 ) )
  call xcolor( 3 )
  first_point = 1
  do i=nstart,nend
    if( xtrajc(i,j) == missing_value ) cycle
    if( first_point == 1) then
      call xpenup( ttrajc(i),xtrajc(i,j)*0.001 )
      first_point = 0
    else
      call xpendn( ttrajc(i),xtrajc(i,j)*0.001 )
    endif
  enddo
  call xchsiz( (y2-y1)*0.04 )
  call xpenup( x1+(x2-x1)*0.10, y2-(y2-y1)*0.1)
  call xpendn( x1+(x2-x1)*0.03, y2-(y2-y1)*0.1)
  call xcharl( x1+(x2-x1)*0.12, y2-(y2-y1)*0.1, 'X')

  call xcolor( 6 )
  first_point = 1
  do i=nstart,nend
    if( ytrajc(i,j) == missing_value ) cycle
    if( first_point == 1) then
      call xpenup( ttrajc(i),ytrajc(i,j)*0.001 )
      first_point = 0
    else
      call xpendn( ttrajc(i),ytrajc(i,j)*0.001 )
    endif
  enddo
  call xchsiz( (y2-y1)*0.04 )
  call xpenup( x1+(x2-x1)*0.10, y2-(y2-y1)*0.15)
  call xpendn( x1+(x2-x1)*0.03, y2-(y2-y1)*0.15)
  call xcharl( x1+(x2-x1)*0.12, y2-(y2-y1)*0.15, 'Y')
  call xwdwof

!-----------------------------------------------------------------------
! Plot Z
!-----------------------------------------------------------------------

  call xcolor(1)
  call xmap(tstart_plt, tend_plt, zlow, zhigh)
  call xqmap(x1,x2,y1,y2)
  call xaxnsz( (y2-y1)*0.05 )
  xo = 12600.0
  yo = 0.0
  call xaxfmt( '(I5)')

  call XAXTIK(-1,-1)
  CALL XAXANT(0,1)

  call XAXISX1(XO,y2, 60.0,600.0)
  call XAXISY1(x2,YO, 0.2,2.0)

  call XAXTIK(1,1)
  CALL XAXANT(-1,-1)

  CALL xbox(x1,x2,y1,y2)
  call xchori(90.0)
  call xcharc( x2+(x2-x1)*0.05, (y1+y2)*0.5, 'Z (m)')
  call xchori( 0.0)

  call xwindw(x1,x2,y1,y2)
  call xcolor( 12 )
  first_point = 1
  do i=nstart,nend
    if( ztrajc(i,j) == missing_value ) cycle
    if( first_point == 1) then
      call xpenup( ttrajc(i),ztrajc(i,j)*0.001 )
      first_point = 0
    else
      call xpendn( ttrajc(i),ztrajc(i,j)*0.001 )
    endif
  enddo
  call xchsiz( (y2-y1)*0.04 )
  call xpenup( x1+(x2-x1)*0.10, y2-(y2-y1)*0.20)
  call xpendn( x1+(x2-x1)*0.03, y2-(y2-y1)*0.20)
  call xcharl( x1+(x2-x1)*0.12, y2-(y2-y1)*0.20, 'Z')
  call xwdwof

!-----------------------------------------------------------------------
! Plot W and Vh
!-----------------------------------------------------------------------

  call xnwpic
  call xcolor(1)
  wmin = -20.0
  wmax = 80.0
  call xmap(tstart_plt, tend_plt, wmin, wmax)
  call xqmap(x1,x2,y1,y2)
  call xchsiz( (y2-y1)*0.06 )
  call xchori(0.0)
  call xcharc( (x1+x2)*0.5, y1-0.15*(y2-y1), 'Time (min)')
  call xchori(90.0)
  call xcharc( x1-(x2-x1)*0.05, (y1+y2)*0.5, 'W and Vh(m/s)')
  call xchori(0.0)

  call xaxnsz( (y2-y1)*0.05 )
  xo = 12600.0
  yo = 0.0
  call xaxfmt( '(I5)')

  call XAXTIK(1,1)
  CALL XAXANT(-1,-1)

  call XAXISX1(XO,y1, 60.0,600.0)
  call XAXISY1(x1,YO, 0.05*(y2-y1),0.10*(y2-y1))

  CALL xbox(x1,x2,y1,y2)

  call xpenup(x1,0.0)
  call xpendn(x2,0.0)
  call xpenup(tref,y1)
  call xpendn(tref,y2)
  call xwindw(x1,x2,y1,y2)

  call xcolor( mod( j+2, 253 ) )
  call xcolor( 3 )
  first_point = 1
  do i=nstart,nend
    if( wtrajc(i,j) == missing_value ) cycle
    if( first_point == 1) then
      call xpenup( ttrajc(i),wtrajc(i,j) )
      first_point = 0
    else
      call xpendn( ttrajc(i),wtrajc(i,j) )
    endif
  enddo
  call xchsiz( (y2-y1)*0.04 )
  call xpenup( x1+(x2-x1)*0.10, y2-(y2-y1)*0.10)
  call xpendn( x1+(x2-x1)*0.03, y2-(y2-y1)*0.10)
  call xcharl( x1+(x2-x1)*0.12, y2-(y2-y1)*0.10, 'W')

  call xwdwof

  call xcolor(1)
! vmin = 0.0
! vmax = 80.0
! call xmap(tstart_plt, tend_plt, vmin, vmax)
! call xqmap(x1,x2,y1,y2)
! call xchsiz( (y2-y1)*0.06 )
! call xchori(90.0)
! call xcharc( x2+(x2-x1)*0.05, (y1+y2)*0.5, 'Vh (m/s)')
! call xchori(0.0)

  call xaxnsz( (y2-y1)*0.05 )

  call XAXTIK(-1,-1)
  CALL XAXANT( 0,0)

  call xaxfmt( '(I5)')
  call XAXISX1(XO,y2, 60.0,600.0)
  call XAXISY1(x2,YO, 0.05*(y2-y1),0.10*(y2-y1))

  call XAXTIK(1,1)
  CALL XAXANT(-1,-1)

  call xpenup(x1,0.0)
  call xpendn(x2,0.0)
  call xpenup(tref,y1)
  call xpendn(tref,y2)

  call xwindw(x1,x2,y1,y2)
  call xcolor( mod( j+2, 253 ) )
  call xcolor( 6 )
  first_point = 1
  do i=nstart,nend
    if( vhtrajc(i,j) == missing_value ) cycle
    if( first_point == 1) then
      call xpenup( ttrajc(i),vhtrajc(i,j))
      first_point = 0
    else
      call xpendn( ttrajc(i),vhtrajc(i,j))
    endif
  enddo
  call xchsiz( (y2-y1)*0.04 )
  call xpenup( x1+(x2-x1)*0.10, y2-(y2-y1)*0.15)
  call xpendn( x1+(x2-x1)*0.03, y2-(y2-y1)*0.15)
  call xcharl( x1+(x2-x1)*0.12, y2-(y2-y1)*0.15, 'Vh')
  call xwdwof

!-----------------------------------------------------------------------
! Plot vorts, vortc, vorth
!-----------------------------------------------------------------------

  call xnwpic
  call xcolor(1)
  vsmin = -0.5
  vsmax =  1.0
  call xmap(tstart_plt, tend_plt, vsmin, vsmax)
  call xqmap(x1,x2,y1,y2)
  call xchsiz( (y2-y1)*0.06 )
  call xchori(0.0)
  call xcharc( (x1+x2)*0.5, y1-0.15*(y2-y1), 'Time (min)')
  call xchori(90.0)
  call xcharc( x1-(x2-x1)*0.05, (y1+y2)*0.5, 'Vorts, Vortc and Vorth (1/s)')
  call xchori(0.0)

  call xaxnsz( (y2-y1)*0.05 )
  xo = 12600.0
  yo = 0.0

  call XAXTIK(1,1)
  CALL XAXANT(-1,-1)

  call xaxfmt( '(I5)')
  call XAXISX1(XO,y1, 60.0,600.0)

  call xaxfmt( '(f10.2)')
  call XAXISY1(x1,YO, 0.05*(y2-y1), 0.1*(y2-y1))

  call XAXTIK(-1,-1)
  CALL XAXANT(0,0)

  call XAXISX1(XO,y2, 60.0,600.0)
  call XAXISY1(x2,YO, 0.05*(y2-y1), 0.1*(y2-y1))

  call XAXTIK(1,1)
  CALL XAXANT(-1,-1)

  CALL xbox(x1,x2,y1,y2)

! call xaxsca(tstart_plt, tend_plt, 60.0, vsmin, vsmax, 0.05*(vsmax-vsmin) )
  call xpenup(x1,0.0)
  call xpendn(x2,0.0)
  call xpenup(tref,y1)
  call xpendn(tref,y2)
  call xwindw(x1,x2,y1,y2)

!-----------------------------------------------------------------------
! Plot streamwise Vort.
!-----------------------------------------------------------------------

  call xcolor( 3 )
  first_point = 1
  do i=nstart,nend
    if( vorts_trajc(i,j) == missing_value ) cycle
    if( first_point == 1) then
      call xpenup( ttrajc(i),vorts_trajc(i,j) )
      first_point = 0
    else
      call xpendn( ttrajc(i),vorts_trajc(i,j) )
    endif
  enddo
  call xchsiz( (y2-y1)*0.04 )
  call xpenup( x1+(x2-x1)*0.10, y2-(y2-y1)*0.10)
  call xpendn( x1+(x2-x1)*0.03, y2-(y2-y1)*0.10)
  call xcharl( x1+(x2-x1)*0.12, y2-(y2-y1)*0.10, 'Streamwise Vort.')

!-----------------------------------------------------------------------
! Plot cross-stream Vort.
!-----------------------------------------------------------------------

  call xcolor( 6 )
  first_point = 1
  do i=nstart,nend
    if( vortc_trajc(i,j) == missing_value ) cycle
    if( first_point == 1) then
      call xpenup( ttrajc(i),vortc_trajc(i,j) )
      first_point = 0
    else
      call xpendn( ttrajc(i),vortc_trajc(i,j) )
    endif
  enddo
  call xchsiz( (y2-y1)*0.04 )
  call xpenup( x1+(x2-x1)*0.10, y2-(y2-y1)*0.15)
  call xpendn( x1+(x2-x1)*0.03, y2-(y2-y1)*0.15)
  call xcharl( x1+(x2-x1)*0.12, y2-(y2-y1)*0.15, 'Cross-stream Vort.')

!-----------------------------------------------------------------------
! Plot Horizontal Vort.
!-----------------------------------------------------------------------

  call xcolor( 12 )
  first_point = 1
  do i=nstart,nend
    if( vortc_trajc(i,j) == missing_value .or. &
        vorts_trajc(i,j) == missing_value ) cycle
    if( first_point == 1) then
      call xpenup( ttrajc(i),sqrt(vorts_trajc(i,j)**2+vortc_trajc(i,j)**2) )
      first_point = 0
    else
      call xpendn( ttrajc(i),sqrt(vorts_trajc(i,j)**2+vortc_trajc(i,j)**2) )
    endif
  enddo
  call xchsiz( (y2-y1)*0.04 )
  call xpenup( x1+(x2-x1)*0.10, y2-(y2-y1)*0.20)
  call xpendn( x1+(x2-x1)*0.03, y2-(y2-y1)*0.20)
  call xcharl( x1+(x2-x1)*0.12, y2-(y2-y1)*0.20, 'Horizontal Vort.')

  call xwdwof

!-----------------------------------------------------------------------
! Plot Vert. and Total . Vort.
!-----------------------------------------------------------------------

  call xnwpic
  call xcolor(1)
  vzmin = -0.5
  vzmax =  2.0
  call xmap(tstart_plt, tend_plt, vzmin, vzmax)
  call xqmap(x1,x2,y1,y2)
  call xchsiz( (y2-y1)*0.06 )
  call xchori(0.0)
  call xcharc( (x1+x2)*0.5, y1-0.15*(y2-y1), 'Time (min)')
  call xchori(90.0)
  call xcharc( x1-(x2-x1)*0.05, (y1+y2)*0.5, 'Vertical & Total Vort (1/s)')
  call xchori(0.0)
  call xaxfmt( '(f10.2)')
  call xaxnsz( (y2-y1)*0.05 )

  xo = 12600.0
  yo = 0.0

  call XAXTIK(1,1)
  CALL XAXANT(-1,-1)

  call xaxfmt( '(I5)')
  call XAXISX1(XO,y1, 60.0,600.0)

  call xaxfmt( '(f10.2)')
  call XAXISY1(x1,YO, 0.05*(y2-y1), 0.1*(y2-y1))

  call XAXTIK(-1,-1)
  CALL XAXANT(0,0)

  call XAXISX1(XO,y2, 60.0,600.0)
  call XAXISY1(x2,YO, 0.05*(y2-y1), 0.1*(y2-y1))

  call XAXTIK(1,1)
  CALL XAXANT(-1,-1)

  CALL xbox(x1,x2,y1,y2)

! call xaxsca(tstart_plt, tend_plt, 60.0, vzmin, vzmax, 0.05*(vzmax-vzmin) )
  call xpenup(x1,0.0)
  call xpendn(x2,0.0)
  call xpenup(tref,y1)
  call xpendn(tref,y2)

!-----------------------------------------------------------------------
! Plot Vert. Vort.
!-----------------------------------------------------------------------

  call xwindw(tstart_plt, tend_plt, vzmin, vzmax )
  call xcolor( mod( j+2, 253 ) )
  call xcolor( 3 )
  first_point = 1
  do i=nstart,nend
    if( vortz_trajc(i,j) == missing_value ) cycle
    if( first_point == 1) then
      call xpenup( ttrajc(i),vortz_trajc(i,j) )
      first_point = 0
    else
      call xpendn( ttrajc(i),vortz_trajc(i,j) )
    endif
  enddo
  call xchsiz( (y2-y1)*0.04 )
  call xpenup( x1+(x2-x1)*0.10, y2-(y2-y1)*0.10)
  call xpendn( x1+(x2-x1)*0.03, y2-(y2-y1)*0.10)
  call xcharl( x1+(x2-x1)*0.12, y2-(y2-y1)*0.10, 'Vertical Vort.')

!-----------------------------------------------------------------------
! Plot Total Vort.
!-----------------------------------------------------------------------

  call xcolor( 6 )
  first_point = 1
  do i=nstart,nend
    if( vortx_trajc(i,j) == missing_value .or. &
        vorty_trajc(i,j) == missing_value .or. &
        vortz_trajc(i,j) == missing_value ) cycle
    if( first_point == 1) then
      call xpenup( ttrajc(i), sqrt(vortx_trajc(i,j)**2+vorty_trajc(i,j)**2+vortz_trajc(i,j)**2) )
      first_point = 0
    else
      call xpendn( ttrajc(i), sqrt(vortx_trajc(i,j)**2+vorty_trajc(i,j)**2+vortz_trajc(i,j)**2) )
    endif
  enddo
  call xchsiz( (y2-y1)*0.04 )
  call xpenup( x1+(x2-x1)*0.10, y2-(y2-y1)*0.15)
  call xpendn( x1+(x2-x1)*0.03, y2-(y2-y1)*0.15)
  call xcharl( x1+(x2-x1)*0.12, y2-(y2-y1)*0.15, 'Total Vort.')
  call xwdwof

!-----------------------------------------------------------------------
! End of time series plotting
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Start plotting 3D trajectories
!-----------------------------------------------------------------------

  call xcolor(11)

  dx = 25.0
  dy = 25.0
  dz = 20.0
  trajc_no_bgn = j
  trajc_no_end = j

  trajc_color = 11

  CALL ZX3DSTRM(dx,dy,dz,xtrajc,ytrajc,ztrajc,ttrajc, npoints, ntrajcs,  &
           xlow*1000, xhigh*1000, ylow*1000, yhigh*1000, zlow*1000, zhigh*1000, nstart, nend , &
           trajc_no_bgn, trajc_no_end, trajc_color )

  trajc_color = 15
  CALL ZX3DSTRM(dx,dy,dz,xtrajc,ytrajc,ztrajc,ttrajc, npoints, ntrajcs,  &
           xlow*1000, xhigh*1000, ylow*1000, yhigh*1000, zlow*1000, zhigh*1000, 1 , nstart, &
           trajc_no_bgn, trajc_no_end, trajc_color )

  trajc_color =  9
  CALL ZX3DSTRM(dx,dy,dz,xtrajc,ytrajc,ztrajc,ttrajc, npoints, ntrajcs,  &
           xlow*1000, xhigh*1000, ylow*1000, yhigh*1000, zlow*1000, zhigh*1000, nend , npoints_indata, &
           trajc_no_bgn, trajc_no_end, trajc_color )

  call xnwfrm

  enddo  ! loop over trajectories J

  call xgrend

  stop
end


SUBROUTINE gtlfnkey( runname, lfnkey )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Find out the number of characters to be used to construct file
!  names.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  03/15/93
!
!  MODIFICATION HISTORY:
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!  12/12/1996 (Yuhe Liu)
!  Removed the restrict of 6 characters to runname.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!

  IMPLICIT NONE

  CHARACTER (LEN=* ) :: runname ! Input
  INTEGER :: lfnkey         ! Output

  INTEGER :: lenstr, firstb, firstc
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  lenstr = LEN( runname )

  firstb = INDEX( runname, ' ')
  firstc = INDEX( runname, ',')

  IF( firstb == 0) firstb = lenstr+1
  IF( firstc == 0) firstc = lenstr+1

  lfnkey = MAX(1, MIN( lenstr, firstb-1, firstc-1 ) )

  RETURN
END SUBROUTINE gtlfnkey

! call xdevic
! CALL xdspac(0.9)
! call xsetclrs(1)

! call xframe

! call xpspac(0.05, 0.95, 0.05, 0.95)

! call xmap(tstart_plt/60.0, tend_plt/60.0, zlow, zhigh)
! call xcolor(1)
! call xcharc( 0.5*(tstart_plt+tend_plt)/60.0, zlow, trim(trajc_fn) )

! write(string,'(6(a,f5.1))') 'x1=',xlow*0.001,', x2=',xhigh*0.001, &
!      ', y1=',ylow*0.001,' y2=',yhigh*0.001,', z1=',zlow*0.001,', z2=',zhigh*0.001

! call xcharc( 0.5*(tstart_plt+tend_plt)/60.0, zlow+0.05*(zhigh-zlow), trim(string) )
! call xcharc( 0.5*(tstart_plt+tend_plt)/60.0, zlow+0.00*(zhigh-zlow), trim(trajc_fn) )

! CALL ZX3DSTRM(dx,dy,dz,xtrajc,ytrajc,ztrajc,ttrajc, npoints, ntrajcs,  &
!          xlow, xhigh, ylow, yhigh, zlow, zhigh, nstart, nend )

! call xgrend


SUBROUTINE ZX3DSTRM(dx,dy,dz,xga,yga,zga,tga, npoints, ntrajcs,  &
           xlow, xhigh, ylow, yhigh, zlow, zhigh, nstart, nend,&
           trajc_no_bgn, trajc_no_end, trajc_color )

  implicit none

  REAL :: dx, dy, dz

  INTEGER :: npoints,ntrajcs
  real xga(npoints,ntrajcs)
  real zga(npoints,ntrajcs)
  real yga(npoints,ntrajcs)
  real tga(npoints)

  REAL :: xlow, xhigh, ylow, yhigh, zlow, zhigh, theta
  REAL :: x1,x2,y1,y2
  REAL :: trajc_no_bgn, trajc_no_end

  character (len=10) :: view

  real ueye,veye,weye,uori,vori,wori,uthi,vthi,wthi

  real normalize, missing_value, eps, zga_normalized1,zga_normalized2

  integer i,j,k,l,n, nstart, nend, trajc_color
  CHARACTER(LEN=80) :: runname
  COMMON /runname_common/runname

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

   print *,'entering zx3dplot plotting routine...'

  missing_value = -99999.0
  eps = 0.001


!-----------------------------------draw lines------------------

  normalize = 1.0*(xhigh-xlow) / (zhigh-zlow)

  ueye = xlow - 2.0*(xhigh-xlow)
  veye = ylow - 2.0*(yhigh-ylow)
  weye = (0.5 *((zlow+zhigh)/2.0))*normalize

  uori = (xlow+xhigh)/2.0
  vori = (ylow+yhigh)/2.0
  wori = zlow*normalize

  uthi = uori
  vthi = vori
  wthi = 2.0*wori

  veye = (ylow+yhigh)/2.0
  weye = ((zlow+zhigh)/4.0 + 2.0*(zhigh-zlow))*normalize

  uori = (xlow+xhigh)/2.0
  vori = (ylow+yhigh)/2.0
  wori = (zlow+zhigh)/4.0*normalize


! do k = 0, 8 ! elevation angle
  do k = 7, 7 ! elevation angle
! do k = 0, 0 ! elevation angle 0 = top
! do k = 8, 8 ! elevation angle

  view = 'south'

  if ( view == 'west' ) then

  theta = 2.0*atan(1.0)*k/8.0
  veye = (ylow+yhigh)/2.0
  ueye = (xlow+xhigh)/2.0 + sin(-theta ) * 4.0*(zhigh-zlow)*normalize
  weye = ((zlow+zhigh)/4.0 + cos(-theta )* 4.0*(zhigh-zlow))*normalize
  uthi = xhigh
  vthi = vori
  wthi = (zlow+zhigh)/2.0*normalize

  elseif( view == 'south' ) then  ! from south

  theta = 2.0*atan(1.0)*k/8.0
  ueye = (xlow+xhigh)/2.0
  veye = (ylow+yhigh)/2.0 + sin(-theta ) * 4.0*(zhigh-zlow)*normalize
  weye = ((zlow+zhigh)/4.0 + cos(-theta )* 4.0*(zhigh-zlow))*normalize
  uthi = uori
  vthi = yhigh
  wthi = (zlow+zhigh)/2.0*normalize

  elseif( view == 'north' )then

  theta = 2.0*atan(1.0)*k/8.0
  ueye = (xlow+xhigh)/2.0
  veye = (ylow+yhigh)/2.0 - sin(-theta ) * 4.0*(zhigh-zlow)*normalize
  weye = ((zlow+zhigh)/4.0 + cos(-theta )* 4.0*(zhigh-zlow))*normalize
  uthi = uori
  vthi = -yhigh
  wthi = (zlow+zhigh)/2.0*normalize

  elseif( view == 'east' )then

  theta = 2.0*atan(1.0)*k/8.0
  veye = (ylow+yhigh)/2.0
  ueye = (xlow+xhigh)/2.0 - sin(-theta ) * 4.0*(zhigh-zlow)*normalize
  weye = ((zlow+zhigh)/4.0 + cos(-theta )* 4.0*(zhigh-zlow))*normalize
  uthi = -xhigh
  vthi = vori
  wthi = (zlow+zhigh)/2.0*normalize

  endif

  L = 20
! do L = 0, 160 ! athimuthal angle
  theta = 8.0*atan(1.0)*L/160.0 - 2.0*atan(1.0)

! do L = 0, 80 ! athimuthal angle
! theta = 8.0*atan(1.0)*L/80.0 - 2.0*atan(1.0)

! do L = 0, 0 ! athimuthal angle
! theta = 8.0*atan(1.0)*L/160.0 - 2.0*atan(1.0)

! ueye = (xlow+xhigh)/2.0 + cos(theta ) * 4.0*(zhigh-zlow)*normalize
! veye = (ylow+yhigh)/2.0 + sin(theta ) * 4.0*(zhigh-zlow)*normalize
! weye = ((zlow+zhigh)/2.0)*normalize
! uthi = uori
! vthi = vori
! wthi = 3*(zlow+zhigh)/4.0*normalize

! print *,' calling tdinit'
  CALL TDINIT(ueye,veye,weye,uori,vori,wori,uthi,vthi,wthi,0)

  call xcolor(1)
! print *,' calling tdgrds start'
! CALL TDGRDS(xlow,ylow,zlow*normalize,xhigh,yhigh,  &
!             zhigh*normalize,dx*10,dy*10, dz*normalize,2,1)

  call tdline( xlow, ylow , zhigh*normalize, xhigh, ylow , zhigh*normalize)
  call tdline( xlow, yhigh, zhigh*normalize, xhigh, yhigh, zhigh*normalize)
  call tdline( xlow, ylow , zhigh*normalize, xlow , yhigh, zhigh*normalize)
  call tdline( xhigh,ylow , zhigh*normalize, xhigh, yhigh, zhigh*normalize)
  call tdline( xlow, ylow , zlow *normalize, xhigh, ylow , zlow *normalize)
  call tdline( xlow, yhigh, zlow *normalize, xhigh, yhigh, zlow *normalize)
  call tdline( xlow, ylow , zlow *normalize, xlow , yhigh, zlow *normalize)
  call tdline( xhigh,ylow , zlow *normalize, xhigh, yhigh, zlow *normalize)
  call tdline( xlow, ylow, zlow*normalize, xlow , ylow, zhigh*normalize)
  call tdline( xlow, yhigh,zlow*normalize, xlow , yhigh,zhigh*normalize)
  call tdline( xhigh,ylow, zlow*normalize, xhigh, ylow, zhigh*normalize)
  call tdline( xhigh,yhigh,zlow*normalize, xhigh, yhigh,zhigh*normalize)

  CALL TDLINE( (xlow+xhigh)*0.5, yhigh,zlow*normalize, &
               (xlow+xhigh)*0.5, yhigh+0.05*(yhigh-ylow),zlow*normalize)

! DO j=nstart+1,nend,5    ! have trajectories grow with frames
  DO j=nend,nend        ! plot all points in one frame

! DO i=1,ntrajcs,2
! DO i=1,73,2
! DO i=4*73+1,5*73,2
! DO i=73,ntrajcs,73
! print *, 'starting to make plot'
! DO i=1,ntrajcs

  DO i=trajc_no_bgn,trajc_no_end

!   call xcolor( mod( i+2, 253 ) )
    call xcolor( trajc_color )

    print*,'plotting trajectory No. ', j, 'from points ', nstart,' to ',j

    DO n=nstart, j-1

!   DO n=361, j-1

      if( .not. (                                  &
          abs(xga(n  ,i)-missing_value).lt.eps.or. &
          abs(yga(n  ,i)-missing_value).lt.eps.or. &
          abs(zga(n  ,i)-missing_value).lt.eps.or. &
          abs(xga(n+1,i)-missing_value).lt.eps.or. &
          abs(yga(n+1,i)-missing_value).lt.eps.or. &
          abs(zga(n+1,i)-missing_value).lt.eps ) ) then

        zga_normalized1 = zga(n  ,i)*normalize
        zga_normalized2 = zga(n+1,i)*normalize

        CALL TDLINE(xga(n  ,i),yga(n  ,i),zga_normalized1, &
                    xga(n+1,i),yga(n+1,i),zga_normalized2)
      endif
    ENDDO
! ENDDO	 ! i=1,ntrajcs

! print *,' calling TDGRDS end'

! print*,'xlow,ylow,zlow,xhigh,yhigh,zhigh*normalize=',xlow,ylow,zlow,xhigh,yhigh,zhigh*normalize
! print*,'dx,dy,dz=', dx,dy,dz

  ENDDO	 ! i=1,ntrajcs

  CALL GSPLCI(1)
! CALL TDGRDS(xlow,ylow,zlow*normalize,xhigh,yhigh,  &
!              zhigh*normalize,dx*10,dy*10,dz*normalize,2,0)


! CALL TDPLCH( 0.5*(xlow+xhigh), ylow-0.05*(yhigh-ylow), trim(runname),0.01*(xhigh-xlow),0.0, 0)
! CALL TDPLCH( 0.5*(xlow+xhigh), 0.5*(ylow+yhigh), 'Caption Here',(xhigh-xlow),0.0, 0)

! call xqmap(x1, x2, y1, y2)
! print*,'x1, x2, y1, y2=', x1, x2, y1, y2
! call xchsiz( 0.04*(y2-y1) )
! call xcharc( 0.5*(x1+x2), y1+0.05*(y2-y1), trim(runname) )

! CALL xFRAME

  ENDDO  ! j - number of points

  enddo ! k - viewing angle

  RETURN
  END

SUBROUTINE ZX3D_cir(dx,dy,dz,xga,yga,zga,tga, npoints, ntrajcs,  &
           xlow, xhigh, ylow, yhigh, zlow, zhigh)

  implicit none

  REAL :: dx, dy, dz

  INTEGER :: npoints,ntrajcs
  real xga(npoints,ntrajcs)
  real zga(npoints,ntrajcs)
  real yga(npoints,ntrajcs)
  real tga(npoints)

  REAL :: xlow, xhigh, ylow, yhigh, zlow, zhigh, theta

  real ueye,veye,weye,uori,vori,wori,uthi,vthi,wthi

  real normalize, missing_value, eps, zga_normalized1,zga_normalized2

  integer i,j,k, n, color_counter, i1
  CHARACTER(LEN=80) :: runname
  COMMON /runname_common/runname

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

   print *,'entering zx3dplot plotting routine...'

  missing_value = -99999.0
  eps = 0.001


!-----------------------------------draw lines------------------

  normalize = 2.0 * (xhigh-xlow) / (zhigh-zlow)

  ueye = xlow - 2.0*(xhigh-xlow)
  veye = ylow - 2.0*(yhigh-ylow)
  weye = (0.5 *((zlow+zhigh)/2.0))*normalize

  uori = (xlow+xhigh)/2.0
  vori = (ylow+yhigh)/2.0
  wori = zlow*normalize

  uthi = uori
  vthi = vori
  wthi = 2.0*wori

  veye = (ylow+yhigh)/2.0
  weye = ((zlow+zhigh)/4.0 + 2.0*(zhigh-zlow))*normalize

  uori = (xlow+xhigh)/2.0
  vori = (ylow+yhigh)/2.0
  wori = (zlow+zhigh)/4.0*normalize


  do k = 0, 20 ! viewing angle   ! elevation angle
  theta = 2.0*atan(1.0)*k/20.0
  ueye = (xlow+xhigh)/2.0
  veye = (ylow+yhigh)/2.0 + sin(-theta ) * 2.0*(zhigh-zlow)*normalize
  weye = ((zlow+zhigh)/4.0 + cos(-theta )* 2.0*(zhigh-zlow))*normalize
  uthi = uori
  vthi = yhigh
  wthi = (zlow+zhigh)/2.0*normalize

! do k = 0, 20 ! athimuthal angle
! theta = 8.0*atan(1.0)*k/20.0
! ueye = (xlow+xhigh)/2.0 + cos(theta ) * 2.0*(zhigh-zlow)*normalize
! veye = (ylow+yhigh)/2.0 + sin(theta ) * 2.0*(zhigh-zlow)*normalize
! weye = ((zlow+zhigh)/4.0)*normalize
! uthi = uori
! vthi = vori
! wthi = (zlow+zhigh)/2.0*normalize

  CALL TDINIT(ueye,veye,weye,uori,vori,wori,uthi,vthi,wthi,0)

  CALL TDGRDS(xlow,ylow,zlow*normalize,xhigh,yhigh,  &
              zhigh*normalize,dx*10,dy*10, dz*normalize,2,1)

  CALL TDLINE( (xlow+xhigh)*0.5, yhigh,zlow*normalize, &
               (xlow+xhigh)*0.5, yhigh+0.01*(yhigh-ylow),zlow*normalize)

  color_counter = 2

! DO n=1,npoints, 50

! DO n=361-100,361+100,5
! DO n=361-100,361,1
  DO n=361,361+10,1

  color_counter = color_counter + 3

  call xcolor( mod( color_counter, 253 ) )

  DO i=1,ntrajcs

      i1 = i + 1
      IF( i == ntrajcs ) i1 = 1

      if( .not. (                                  &
          abs(xga(n  ,i)-missing_value).lt.eps.or. &
          abs(yga(n  ,i)-missing_value).lt.eps.or. &
          abs(zga(n  ,i)-missing_value).lt.eps.or. &
          abs(xga(n,i1)-missing_value).lt.eps.or. &
          abs(yga(n,i1)-missing_value).lt.eps.or. &
          abs(zga(n,i1)-missing_value).lt.eps ) ) then

        zga_normalized1 = zga(n  ,i)*normalize
        zga_normalized2 = zga(n,i1)*normalize

        CALL TDLINE(xga(n  ,i),yga(n  ,i),zga_normalized1, &
                    xga(n,i1),yga(n,i1),zga_normalized2)
      endif

  ENDDO	 ! i=1,ntrajcs-1

  ENDDO	 ! n -  points

  CALL GSPLCI(1)
  CALL TDGRDS(xlow,ylow,zlow*normalize,xhigh,yhigh,  &
               zhigh*normalize,dx*10,dy*10,dz*normalize,2,0)

  CALL FRAME

  enddo ! k - viewing angle

  RETURN
  END
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SMOOTH3PMV                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE smooth3pmv( arr,nx,ibgn,iend,missing_value, tem1 )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Smooth a 1-D array by the filter of { 1 2 1 }
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!
!  5/26/2004.
!
!  Modification History
!  Based on smooth9pmv.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  nx       Number of grid points in the x-direction
!  ibgn     First index in x-direction in the soomthing region.
!  iend     Last  index in x-direction in the soomthing region.
!
!  arr    1-D array
!  missing_value
!
!  OUTPUT:
!
!  arr    1-D array
!
!  TEMPORARY:
!
!  tem1     Temporary 1-D array
!
!-----------------------------------------------------------------------
!  Variable Declarations.
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: nx         ! Number of grid points in the x-direction
  INTEGER :: ibgn
  INTEGER :: iend
!
  REAL :: arr (nx)   ! 1-D array
  REAL :: tem1(nx)   ! Temporary array
  REAL :: missing_value, mv
!
!-----------------------------------------------------------------------
!  Misc. local variables:
!-----------------------------------------------------------------------
!
  INTEGER :: i,ip,im
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  Beginning of executable code...
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  mv = missing_value

  DO i=1,nx
    IF( ABS(arr(i)-mv) <= 1.0e-10 ) arr(i)=mv
  END DO

  DO i = ibgn,iend
    ip=MIN(nx,i+1)
    im=MAX( 1,i-1)
    IF(.not. (arr(im)==mv.OR.arr(i)==mv.OR.arr(ip)==mv)) then
      tem1(i) =0.25*(arr(im)+2.*arr(i)+arr(ip))
    END IF
  END DO

  DO i = ibgn,iend
    arr(i) = tem1(i)
  END DO

  RETURN
END SUBROUTINE smooth3pmv

