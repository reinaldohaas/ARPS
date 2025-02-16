PROGRAM ARPSTRAJC
!
!##################################################################
!##################################################################
!######                                                      ######
!######                   PROGRAM ARPSTRAJC                  ######
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
!  (4/08/2004)
!
!  MODIFICATION HISTORY:
!
!  9/28/2005
!  Modified to handle terrain-following grid. 3D zp is defined.
!  Removed used dignositc codes and arrays
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  REAL, ALLOCATABLE :: x(:),y(:),z(:)
  REAL, ALLOCATABLE :: xs (:), ys (:) ! x,y coord for scalar points

  REAL, ALLOCATABLE :: u(:,:,:)
  REAL, ALLOCATABLE :: v(:,:,:)
  REAL, ALLOCATABLE :: w(:,:,:)

  REAL, ALLOCATABLE :: zp(:,:,:)
  REAL, ALLOCATABLE :: tem1(:,:,:)
  REAL, ALLOCATABLE :: tem2(:,:,:)
  REAL, ALLOCATABLE :: tem3(:,:,:)
  REAL, ALLOCATABLE :: tem4(:,:,:)
  REAL, ALLOCATABLE :: tem5(:,:,:)
  REAL, ALLOCATABLE :: tem6(:,:,:)

! REAL, allocatable :: ptbar(:),pbar(:),rhobar(:),qvbar(:)

  INTEGER (KIND=selected_int_kind(4)), ALLOCATABLE :: itmp(:,:,:) ! Temporary array

  INTEGER :: nx,ny,nz,nzsoil   ! Grid dimensions.
  INTEGER :: nstyps            ! Maximum number of soil types.

  INTEGER :: hinfmt,nhisfile_max,nhisfile,lengbf
  PARAMETER (nhisfile_max=2000)
  CHARACTER (LEN=256) :: grdbasfn,hisfile(nhisfile_max)

  INTEGER :: hdmpinopt, nf
  CHARACTER (LEN=256) :: hdmpftrailer
  CHARACTER (LEN=256) :: hdmpfheader
  REAL :: tintv_dmpin, tbgn_dmpin, tend_dmpin

  NAMELIST /history_data/ hinfmt, nhisfile, grdbasfn, hisfile,          &
           hdmpinopt,hdmpfheader,hdmpftrailer,                          &
           tintv_dmpin, tbgn_dmpin, tend_dmpin

  INTEGER, PARAMETER :: nmax_trajcs=1000, nmax_times=20
  REAL :: xtrajc0(nmax_trajcs,nmax_times),                              &
          ytrajc0(nmax_trajcs,nmax_times),                              &
          ztrajc0(nmax_trajcs,nmax_times),ttrajc0
  INTEGER :: ntrajc0,ntrajc0_sub,ntrajcs, npoints
!
! initrajc = 1 !  specified locations
!          = 2 !  circle of radius radius and athimuthal angle increments of theta_inc
!          = 3 !  square of radius radius and athimuthal angle increments of theta_inc
!          = 4 !  read in from trajectory data file
!
  CHARACTER(LEN=256) :: trajc_fn_in(nmax_times)

  REAL :: radius,xwidth,ywidth,xctr,yctr,zctr(100),theta_inc
  INTEGER :: initrajc,nzctr, ntimes,subinterval
  CHARACTER (LEN=256) :: trajcfn_header
  REAL :: reftime(nmax_times)

  NAMELIST /trajectories/ trajcfn_header,ntimes,subinterval, reftime, initrajc, ntrajcs, &
  xtrajc0,ytrajc0,ztrajc0,ttrajc0,ntrajc0, &
  radius,xwidth,ywidth,xctr,yctr,zctr,nzctr,theta_inc,dirname,trajc_fn_in

  REAL, allocatable :: xtrajc(:,:,:),ytrajc(:,:,:),ztrajc(:,:,:),ttrajc(:)
  REAL, allocatable :: utrajc(:,:,:),vtrajc(:,:,:),wtrajc(:,:,:)

! REAL, allocatable :: xweight(:,:),yweight(:,:),zweight(:,:)
! INTEGER, allocatable :: itrajc(:,:),jtrajc(:,:),ktrajc(:,:)

  INTEGER :: istatus, cur_level
  REAL :: tinc

  REAL :: min_pte
  INTEGER :: min_pte_i, min_pte_j
  INTEGER :: ibgn,iend,jbgn,jend

  REAL :: f_tdew
!
!-----------------------------------------------------------------------
!  Include files:
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'grid.inc'

  INTEGER :: i,j,k,ns,ireturn

  REAL :: time, time1, time2, amin, amax, degree, degree2radian, theta
  LOGICAL :: iexist


  REAL :: u1_changed, utem,vtem,wtem, factor
  INTEGER :: nfile,nsub

  REAL :: ttrajc0_in
  INTEGER :: ntrajcs_in, ntrajcs_each_j

  CHARACTER (LEN=80) :: timsnd1,timsnd2,timsnd0,timsnd3
  INTEGER :: tmstrln1,tmstrln2,tmstrln0, nunit(nmax_times),ltrajc_fn,istat,tmstrln3
  INTEGER :: kl, notinteg, nunit1(nmax_times)

  CHARACTER (LEN=256) :: trajc_fn

  CHARACTER (LEN=40) :: varunits, var_name

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
     '###                Welcome to ARPSTRAJC                     ###', &
     '###                                                         ###', &
     '###############################################################', &
     '###############################################################'

  CALL mpinit_var

!
!-----------------------------------------------------------------------
!  Get the names of the input data files.
!-----------------------------------------------------------------------
!

  READ(5,history_data,ERR=100)

  WRITE(6,'(a,i3)') 'Input hdmpinopt=', hdmpinopt
  WRITE(6,'(a,i3)') 'Input hinfmt   =', hinfmt
  WRITE(6,'(a,i3)') 'Input nhisfile =', nhisfile

  IF( hinfmt /= 3 ) then
    Print*,'Only HDF format (3) is supported.'
    STOP
  ENDIF

  WRITE(6,'(a)')'Namelist history_data was successfully read.'
  IF( hdmpinopt == 1 ) THEN

    CALL gthinfns(hdmpfheader,hdmpftrailer,hinfmt,                      &
                  tintv_dmpin, tbgn_dmpin, tend_dmpin,                  &
                  grdbasfn,hisfile,nhisfile)

  END IF

  IF( nhisfile > nhisfile_max ) THEN

    WRITE(6,'(/a,a,I5,/a/)')                                            &
        'The number of history files to be processed',                  &
        ' exceeded the maximum number ',nhisfile_max,                   &
        'Please increase the size of array hisfile.'

  END IF

  lengbf =len_trim( grdbasfn )
  WRITE(6,'(/a,a)')'The grid base-state file name is ',                 &
                    grdbasfn(1:lengbf)

  DO nfile=1,nhisfile
    WRITE(6,'(a,i3,a,a)')                                               &
        'History file No.',nfile,' is ',trim(hisfile(nfile))
  END DO


  initrajc = 1

  reftime = 0.0
  trajcfn_header = ' '
  nzctr = 1
  ntimes = 1

  ttrajc0=0.0
  ntrajc0=1

  xtrajc0(:,:)=0.0
  ytrajc0(:,:)=0.0
  ztrajc0(:,:)=0.0

  READ(5,trajectories, END=100)
  WRITE(6,'(/a,a/)') 'NAMELIST block trajectories successfully read.'

  trajc_fn_in(:)='default_trajecotry_filename'

!  write(6,*) 'Print out namelist block trajectories parameters read in.'
!  write(6,trajectories)

  IF( ntimes > nmax_times ) then
    print*,'ntimes in input file exceeded max allowed value of ', nmax_times
    STOP
  ENDIF

  IF( initrajc == 1 .and. ntrajcs > nmax_trajcs ) then
    Print*,'ntrajcs larger than maximum allowed. Increase nmax_trajcs in program to ',ntrajcs
    STOP
  endif

  degree2radian = 4.0*atan(1.0)/180.0

  IF( initrajc == 2 ) then  ! circle

    ntrajcs_each_j = nint( 360.0/theta_inc ) + 1
    ntrajcs = nzctr * ntrajcs_each_j

    IF( ntrajcs > nmax_trajcs ) then
      Print*,'ntrajcs larger than maximum allowed. Increase nmax_trajcs in program to ',ntrajcs
      STOP
    endif

    print*,'ntrajcs,nzctr,initrajc=',ntrajcs,nzctr,initrajc

    do k=1,ntimes
    do j=1,nzctr

      xtrajc0((j-1)*(ntrajcs_each_j)+1,k)= xctr
      ytrajc0((j-1)*(ntrajcs_each_j)+1,k)= yctr
      ztrajc0((j-1)*(ntrajcs_each_j)+1,k)= zctr(j)

      do i=1,ntrajcs
        theta = (i-1)*theta_inc * degree2radian
        xtrajc0((j-1)*(ntrajcs_each_j)+(i+1),k)= xctr+radius*sin( theta )
        ytrajc0((j-1)*(ntrajcs_each_j)+(i+1),k)= yctr+radius*cos( theta )
        ztrajc0((j-1)*(ntrajcs_each_j)+(i+1),k)= zctr(j)
      enddo
    enddo
    enddo

  endif

  IF( initrajc == 3 ) then  ! box
!    ntrajcs = nint( 360.0/theta_inc )

    ! DTD: changed to a box that has trajectories at even increments within it:
    ! In this case theta_inc is interpreted as the distance in x and y between
    ! trajectories inside the box.

    ntrajcs = ((2.*radius)/theta_inc+1.)**2.

    IF( ntrajcs > nmax_trajcs ) then
      Print*,'ntrajcs larger than maximum allowed. Increase nmax_trajcs in program to ',ntrajcs
      STOP
    endif

    do k=1,ntimes
    do i=1,ntrajcs

    xtrajc0(i,k) = xctr-radius+MOD(i-1,INT(2.*radius/theta_inc+1.))*theta_inc
    ytrajc0(i,k) = yctr-radius+INT((i-1)/(2.*radius/theta_inc+1.))*theta_inc

!      degree = (i-1)*theta_inc
!      theta = degree * degree2radian

!      if( degree <= 45.0 ) then
!        xtrajc0(i,k)= xctr+radius
!        ytrajc0(i,k)= yctr+radius*tan( theta )
!      elseif( degree <= 135.0 ) then
!        xtrajc0(i,k)= xctr+radius*cos( theta )
!        ytrajc0(i,k)= yctr+radius
!      elseif( degree <= 225.0 ) then
!        xtrajc0(i,k)= xctr-radius
!        ytrajc0(i,k)= yctr+radius*tan( theta )
!      elseif( degree <= 315.0 ) then
!        xtrajc0(i,k)= xctr+radius*cos( theta )
!        ytrajc0(i,k)= yctr-radius
!      else
!        xtrajc0(i,k)= xctr+radius
!        ytrajc0(i,k)= yctr+radius*tan( theta )
!      endif

      ztrajc0(i,k)=zctr(1)
    enddo
    enddo
  endif

! print*,'I am here. initrajc, ntrajcs=',  initrajc, ntrajcs

  IF( initrajc == 4 ) then  ! read from file

  DO k=1,ntimes

  CALL getunit(nunit(k))
  OPEN(UNIT=nunit(k),FILE=trim(trajc_fn_in(k)),STATUS='old',FORM='formatted')

  READ(nunit(k),*)
  READ(nunit(k),*)
  READ(nunit(k),*)
  READ(nunit(k),*)
  READ(nunit(k),*)
  READ(nunit(k),*)

  READ(nunit(k),'(4e17.6)') ttrajc0_in
  READ(nunit(k),'(i10)') ntrajcs_in
  print*,'ntrajcs_in=', ntrajcs_in

  IF( ntrajcs_in /= ntrajcs ) then
    print*,'ntrajcs read in .ne. ntrajcs in program. Value in data used.'
    ntrajcs = ntrajcs_in
  ENDIF
  READ(nunit(k),'(6e17.6)') ((xtrajc0(i,k),ytrajc0(i,k),ztrajc0(i,k)),i=1,ntrajcs)

  CLOSE(nunit(k))
  call retunit(nunit(k))

  ENDDO

  ENDIF

! DO k=1,ntimes
! print*,'I am here now. ntrajcs=',ntrajcs,' No. of times=', k
! write(6,'(6e17.6)') ((xtrajc0(i,k),ytrajc0(i,k),ztrajc0(i,k)),i=1,ntrajcs)
! ENDDO

  IF( hdmpinopt == 1 ) THEN

    CALL gthinfns(hdmpfheader,hdmpftrailer,hinfmt,                      &
                  tintv_dmpin, tbgn_dmpin, tend_dmpin,                  &
                  grdbasfn,hisfile,nhisfile)

  END IF

  IF( nhisfile > nhisfile_max ) THEN

    WRITE(6,'(/a,a,I5,/a/)')                                            &
        'The number of history files to be processed',                  &
        ' exceeded the maximum number ',nhisfile_max,                   &
        'Please increase the size of array hisfile.'

  END IF

  lengbf =len_trim( grdbasfn )
  WRITE(6,'(/a,a)')'The grid base-state file name is ',                 &
                    grdbasfn(1:lengbf)

  DO nf=1,nhisfile
    WRITE(6,'(a,i3,a,a)')                                               &
        'History file No.',nf,' is ',trim(hisfile(nf))
  END DO


  npoints = nhisfile*subinterval-(subinterval-1)

  lengbf = len_trim(grdbasfn)
  WRITE(6,'(/a,a)')' The grid/base name is ', trim(grdbasfn)

  CALL get_dims_from_data(hinfmt,grdbasfn(1:lengbf),                    &
       nx,ny,nz,nzsoil,nstyps, ireturn)

  Print*,'nx,ny,nz of input data were ', nx,ny,nz

  allocate(x(nx ),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:x")
  x = 0.0

  allocate(y(ny ),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:y")
  y = 0.0

  allocate(z(nz ),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:z")
  z = 0.0

  allocate(xs(nx),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:xs")
  xs = 0.0
  allocate(ys(ny),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:ys")
  ys = 0.0

  allocate(u(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:u")
  u=0.0
  allocate(v(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:v")
  v=0.0
  allocate(w(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:w")
  w=0.0

  allocate(zp(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:zp")
  zp=0.0

  allocate(itmp(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:itmp")
  itmp=0

  allocate(tem1(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:tem1")

  allocate(tem2(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:tem2")

  allocate(tem3(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:tem3")

  allocate(tem4(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:tem4")

  allocate(tem5(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:tem5")

  allocate(tem6(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:tem6")

  if( hdmpinopt == 1 ) then ! determine ntrajc0  from ttrajc0
!   ttrajc0 is known
    ntrajc0 = nint( (ttrajc0 - tbgn_dmpin)/tintv_dmpin ) + 1
  else
!   ntrajc0 is known from namelist input
!   Get ttrajc0  from NO. ntrajc0 history data
  endif

  CALL get_gridinfo_from_hdf(grdbasfn(1:lengbf),nx,ny,nz,x,y,z,ireturn)

  Print*,'x,y,z of input data read in.'
  print*,'x(1 )=',x(1)
  print*,'x(nx)=',x(nx)
  print*,'y(1 )=',y(1)
  print*,'y(ny)=',y(ny)

  dx = x(2) - x(1)
  dy = y(2) - y(1)

  IF(initrajc == 5) THEN  ! DTD: Trajectories initialized near surface in low-level cold pool
                          ! by choosing points at regular intervals below a certain pte threshold

    ! Note: this option assumes there is only one trajectory initial time for now

    k = 0
    nfile = ntrajc0
    CALL hdfreadvar(nx,ny,nz,trim(grdbasfn),0,"zp", zp, itmp )

! DTD: commented out below for now and read in file that contains pteprt
! explicitly

    time=reftime(1)
    CALL readvar2(nx,ny,nz,tem5,"pteprt",var_name,    &
                      varunits, time, trajcfn_header, dirname,            &
                      3, 0, istatus)

!    CALL hdfreadvar(nx,ny,nz,trim(grdbasfn),time,'pbar', tem1, itmp )
!    CALL hdfreadvar(nx,ny,nz,trim(grdbasfn),time,'ptbar', tem2, itmp )
!    CALL hdfreadvar(nx,ny,nz,trim(grdbasfn),time,'qvbar', tem3, itmp )

!    CALL pt2pte(nx,ny,nz,1,nx-1,1,ny-1,1,nz-1,tem1,tem2,tem3,tem4)

!    CALL hdfreadvar(nx,ny,nz,trim(hisfile(nfile)),time,'p', tem1, itmp )
!    CALL hdfreadvar(nx,ny,nz,trim(hisfile(nfile)),time,'pt', tem2, itmp )
!    CALL hdfreadvar(nx,ny,nz,trim(hisfile(nfile)),time,'qv', tem3, itmp )

!    CALL pt2pte(nx,ny,nz,1,nx-1,1,ny-1,1,nz-1,tem1,tem2,tem3,tem5)

    ibgn=2
    iend=nx-2
    jbgn=2
    jend=ny-2

    DO i=1,nx-1
      xs(i)=0.5*(x(i)+x(i+1))
      IF((xctr-xwidth/2.) >= xs(i)) THEN
        ibgn = i
      END IF
      IF((xctr+xwidth/2.) >= xs(i)) THEN
        iend = i
      END IF
    END DO

    DO j=1,ny-1
      ys(j)=0.5*(y(j)+y(j+1))
      IF((yctr-ywidth/2.) >= ys(j)) THEN
        jbgn = j
      END IF
      IF((yctr+ywidth/2.) >= ys(j)) THEN
        jend = j
      END IF
    END DO

    print*,'ibgn,iend',ibgn,iend
    print*,'jbgn,jend',jbgn,jend

    DO i=ibgn,iend
      DO j=jbgn,jend
!    DO i=2,nx-2,1
!      DO j=2,ny-2,1
!        tem5(i,j,2) = tem5(i,j,2)-tem4(i,j,2)
        !print*,tem1(i,j,2)
        !WRITE(0,*)'thetae_prt at i,j,2',i,j,tem5(i,j,2)
        IF(tem5(i,j,2) < -4.0) THEN
          k = k + 1
          xtrajc0(k,1) = 0.5*(x(i)+x(i+1))
          ytrajc0(k,1) = 0.5*(y(j)+y(j+1))
          ztrajc0(k,1) = (zp(i,j,2)+zp(i,j,3))*0.5-zp(i,j,2)
        END IF
      END DO
    END DO

    ntrajcs = k

  END IF

  IF(initrajc == 6) THEN  ! DTD: Start a trajectory from the minimum surface theta_e
                          ! Also from minimum surface T and minimum surface Td
                          ! In the future I will identify areas surrounding these points
                          ! and start other trajectories which can then be averaged.

    ! Note: this option assumes there is only one trajectory initial time for now

    nfile = ntrajc0
    ntrajcs = 9

    ! Minimum theta_e

    CALL hdfreadvar(nx,ny,nz,trim(grdbasfn),0,"zp", zp, itmp )

    CALL hdfreadvar(nx,ny,nz,trim(hisfile(nfile)),time,'p', tem1, itmp )
    CALL hdfreadvar(nx,ny,nz,trim(hisfile(nfile)),time,'pt', tem2, itmp )
    CALL hdfreadvar(nx,ny,nz,trim(hisfile(nfile)),time,'qv', tem3, itmp )

    CALL pt2pte(nx,ny,nz,1,nx-1,1,ny-1,1,nz-1,tem1,tem2,tem3,tem5)

    min_pte = tem5(2,2,2)

    ibgn=2
    iend=nx-2
    jbgn=2
    jend=ny-2

    DO i=1,nx-1
      xs(i)=0.5*(x(i)+x(i+1))
      IF((xctr-xwidth/2.) >= xs(i)) THEN
        ibgn = i
      END IF
      IF((xctr+xwidth/2.) >= xs(i)) THEN
        iend = i
      END IF
    END DO

    DO j=1,ny-1
      ys(j)=0.5*(y(j)+y(j+1))
      IF((yctr-ywidth/2.) >= ys(j)) THEN
        jbgn = j
      END IF
      IF((yctr+ywidth/2.) >= ys(j)) THEN
        jend = j
      END IF
    END DO

    print*,'ibgn,iend',ibgn,iend
    print*,'jbgn,jend',jbgn,jend

    ! Find the grid point with minimum surface theta_e
    ! within the constraints of the given bounding box
    DO i=ibgn,iend
      DO j=jbgn,jend
        IF(tem5(i,j,2) < min_pte) THEN
          min_pte = tem5(i,j,2)
          min_pte_i = i
          min_pte_j = j
        END IF
      END DO
    END DO

    xtrajc0(1,1) = 0.5*(x(min_pte_i)+x(min_pte_i+1))
    ytrajc0(1,1) = 0.5*(y(min_pte_j)+y(min_pte_j+1))
    DO i=1,9
      ztrajc0(i,1) = (zp(min_pte_i,min_pte_j,2)+zp(min_pte_i,min_pte_j,3))*0.5-zp(min_pte_i,min_pte_j,2)    ! initial trajectory at height of first scalar point above ground
    END DO
    ! Also start trajectories at 8 surrounding scalar points

!    xtrajc0(2,1) = 0.5*(x(min_pte_i-1)+x(min_pte_i))
!    ytrajc0(2,1) = 0.5*(y(min_pte_j)+y(min_pte_j+1))

!    xtrajc0(3,1) = 0.5*(x(min_pte_i-1)+x(min_pte_i))
!    ytrajc0(3,1) = 0.5*(y(min_pte_j-1)+y(min_pte_j))

!    xtrajc0(4,1) = 0.5*(x(min_pte_i)+x(min_pte_i+1))
!    ytrajc0(4,1) = 0.5*(y(min_pte_j-1)+y(min_pte_j))

!    xtrajc0(5,1) = 0.5*(x(min_pte_i+1)+x(min_pte_i+2))
!    ytrajc0(5,1) = 0.5*(y(min_pte_j-1)+y(min_pte_j))

!    xtrajc0(6,1) = 0.5*(x(min_pte_i+1)+x(min_pte_i+2))
!    ytrajc0(6,1) = 0.5*(y(min_pte_j)+y(min_pte_j+1))

!    xtrajc0(7,1) = 0.5*(x(min_pte_i+1)+x(min_pte_i+2))
!    ytrajc0(7,1) = 0.5*(y(min_pte_j+1)+y(min_pte_j+2))

!    xtrajc0(8,1) = 0.5*(x(min_pte_i)+x(min_pte_i+1))
!    ytrajc0(8,1) = 0.5*(y(min_pte_j+1)+y(min_pte_j+2))

!    xtrajc0(9,1) = 0.5*(x(min_pte_i-1)+x(min_pte_i))
!    ytrajc0(9,1) = 0.5*(y(min_pte_j+1)+y(min_pte_j+2))

    ! Start trajectories at 8 points in a box surrounding central point

    xtrajc0(2,1) = xtrajc0(1,1)-0.50*dx
    ytrajc0(2,1) = ytrajc0(1,1)

    xtrajc0(3,1) = xtrajc0(1,1)-0.50*dx
    ytrajc0(3,1) = ytrajc0(1,1)-0.50*dy

    xtrajc0(4,1) = xtrajc0(1,1)
    ytrajc0(4,1) = ytrajc0(1,1)-0.50*dy

    xtrajc0(5,1) = xtrajc0(1,1)+0.50*dx
    ytrajc0(5,1) = ytrajc0(1,1)-0.50*dy

    xtrajc0(6,1) = xtrajc0(1,1)+0.50*dx
    ytrajc0(6,1) = ytrajc0(1,1)

    xtrajc0(7,1) = xtrajc0(1,1)+0.50*dx
    ytrajc0(7,1) = ytrajc0(1,1)+0.50*dy

    xtrajc0(8,1) = xtrajc0(1,1)
    ytrajc0(8,1) = ytrajc0(1,1)+0.50*dy

    xtrajc0(9,1) = xtrajc0(1,1)-0.50*dx
    ytrajc0(9,1) = ytrajc0(1,1)+0.50*dy

    ! Minimum T

!    DO k=1,nz-1
!      DO j=1,ny-1
!        DO i=1,nx-1
!          tem5(i,j,k) = tem2(i,j,k)*(tem1(i,j,k)/p0)**rddcp   ! temperature is in tem5
!        END DO
!      END DO
!    END DO

!    !Re-use min_pte for temperature here
!    min_pte = tem5(2,2,2)

    ! Find the grid point with minimum surface temperature
    ! within the constraints of the given bounding box
!    DO i=MAX(2,INT((xctr-radius)/dx)),MIN(nx-2,INT((xctr+radius)/dx))
!      DO j=MAX(2,INT((yctr-radius)/dy)),MIN(nx-2,INT((yctr+radius)/dy))
!        IF(tem5(i,j,2) < min_pte) THEN
!          min_pte = tem5(i,j,2)
!          min_pte_i = i
!          min_pte_j = j
!        END IF
!      END DO
!    END DO

!    xtrajc0(2,1) = x(min_pte_i)
!    ytrajc0(2,1) = y(min_pte_j)
!    ztrajc0(2,1) = (zp(i,j,2)+zp(i,j,3))*0.5    ! initial trajectory at height of first scalar point above ground

    ! Minimum Td

!    DO k=1,nz-1
!      DO j=1,ny-1
!        DO i=1,nx-1
!          tem5(i,j,k) = f_tdew(tem1(i,j,k),tem5(i,j,k),tem3(i,j,k))  ! dewpoint is in tem5
!        END DO
!      END DO
!    END DO

    !Re-use min_pte for dewpoint here
!    min_pte = tem5(2,2,2)

    ! Find the grid point with minimum surface dewpoint
    ! within the constraints of the given bounding box
!    DO i=MAX(2,INT((xctr-radius)/dx)),MIN(nx-2,INT((xctr+radius)/dx))
!      DO j=MAX(2,INT((yctr-radius)/dy)),MIN(nx-2,INT((yctr+radius)/dy))
!        IF(tem5(i,j,2) < min_pte) THEN
!          min_pte = tem5(i,j,2)
!          min_pte_i = i
!          min_pte_j = j
!        END IF
!      END DO
!    END DO

!    xtrajc0(3,1) = x(min_pte_i)
!    ytrajc0(3,1) = y(min_pte_j)
!    ztrajc0(3,1) = (zp(i,j,2)+zp(i,j,3))*0.5    ! initial trajectory at height of first scalar point above ground

  END IF

  IF (initrajc == 7) THEN ! Start trajectories in and near a point at the minimum w at a certain
                          ! k-level within a given bounding box

    ! Note: this option assumes there is only one trajectory initial time for now

    nfile = ntrajc0
    ntrajcs = 9

    ! Minimum w

    CALL hdfreadvar(nx,ny,nz,trim(grdbasfn),0,"zp", zp, itmp )

    CALL hdfreadvar(nx,ny,nz,trim(hisfile(nfile)),time,'w', tem1, itmp )

    ibgn=2
    iend=nx-2
    jbgn=2
    jend=ny-2

    DO i=1,nx-1
      xs(i)=0.5*(x(i)+x(i+1))
      IF((xctr-xwidth/2.) >= xs(i)) THEN
        ibgn = i
      END IF
      IF((xctr+xwidth/2.) >= xs(i)) THEN
        iend = i
      END IF
    END DO

    DO j=1,ny-1
      ys(j)=0.5*(y(j)+y(j+1))
      IF((yctr-ywidth/2.) >= ys(j)) THEN
        jbgn = j
      END IF
      IF((yctr+ywidth/2.) >= ys(j)) THEN
        jend = j
      END IF
    END DO

    print*,'ibgn,iend',ibgn,iend
    print*,'jbgn,jend',jbgn,jend

    ! "Re-use" min_pte here for min_w
    min_pte = tem1(2,2,6)

    ! Find the grid point with minimum w
    ! within the constraints of the given bounding box
    DO i=ibgn,iend
      DO j=jbgn,jend
        IF(tem1(i,j,6) < min_pte) THEN
          min_pte = tem1(i,j,6)
          min_pte_i = i
          min_pte_j = j
        END IF
      END DO
    END DO

    xtrajc0(1,1) = 0.5*(x(min_pte_i)+x(min_pte_i+1))
    ytrajc0(1,1) = 0.5*(y(min_pte_j)+y(min_pte_j+1))
    DO i=1,9
      ztrajc0(i,1) = zp(min_pte_i,min_pte_j,6)
    END DO
    ! Also start trajectories at 8 surrounding scalar points

!    xtrajc0(2,1) = 0.5*(x(min_pte_i-1)+x(min_pte_i))
!    ytrajc0(2,1) = 0.5*(y(min_pte_j)+y(min_pte_j+1))

!    xtrajc0(3,1) = 0.5*(x(min_pte_i-1)+x(min_pte_i))
!    ytrajc0(3,1) = 0.5*(y(min_pte_j-1)+y(min_pte_j))

!    xtrajc0(4,1) = 0.5*(x(min_pte_i)+x(min_pte_i+1))
!    ytrajc0(4,1) = 0.5*(y(min_pte_j-1)+y(min_pte_j))

!    xtrajc0(5,1) = 0.5*(x(min_pte_i+1)+x(min_pte_i+2))
!    ytrajc0(5,1) = 0.5*(y(min_pte_j-1)+y(min_pte_j))

!    xtrajc0(6,1) = 0.5*(x(min_pte_i+1)+x(min_pte_i+2))
!    ytrajc0(6,1) = 0.5*(y(min_pte_j)+y(min_pte_j+1))

!    xtrajc0(7,1) = 0.5*(x(min_pte_i+1)+x(min_pte_i+2))
!    ytrajc0(7,1) = 0.5*(y(min_pte_j+1)+y(min_pte_j+2))

!    xtrajc0(8,1) = 0.5*(x(min_pte_i)+x(min_pte_i+1))
!    ytrajc0(8,1) = 0.5*(y(min_pte_j+1)+y(min_pte_j+2))

!    xtrajc0(9,1) = 0.5*(x(min_pte_i-1)+x(min_pte_i))
!    ytrajc0(9,1) = 0.5*(y(min_pte_j+1)+y(min_pte_j+2))

    ! Start trajectories at 8 points in a box surrounding central point

    xtrajc0(2,1) = xtrajc0(1,1)-0.50*dx
    ytrajc0(2,1) = ytrajc0(1,1)

    xtrajc0(3,1) = xtrajc0(1,1)-0.50*dx
    ytrajc0(3,1) = ytrajc0(1,1)-0.50*dy

    xtrajc0(4,1) = xtrajc0(1,1)
    ytrajc0(4,1) = ytrajc0(1,1)-0.50*dy

    xtrajc0(5,1) = xtrajc0(1,1)+0.50*dx
    ytrajc0(5,1) = ytrajc0(1,1)-0.50*dy

    xtrajc0(6,1) = xtrajc0(1,1)+0.50*dx
    ytrajc0(6,1) = ytrajc0(1,1)

    xtrajc0(7,1) = xtrajc0(1,1)+0.50*dx
    ytrajc0(7,1) = ytrajc0(1,1)+0.50*dy

    xtrajc0(8,1) = xtrajc0(1,1)
    ytrajc0(8,1) = ytrajc0(1,1)+0.50*dy

    xtrajc0(9,1) = xtrajc0(1,1)-0.50*dx
    ytrajc0(9,1) = ytrajc0(1,1)+0.50*dy
  END IF

  DO k=1,ntimes
  print*,'I am here now. ntrajcs=',ntrajcs,' No. of times=', k
  write(6,'(6e17.6)') ((xtrajc0(i,k),ytrajc0(i,k),ztrajc0(i,k)),i=1,ntrajcs)
  ENDDO

  print*,'npoints,ntrajcs=', npoints,ntrajcs

  allocate(xtrajc(npoints,ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:xtrajc")

  allocate(ytrajc(npoints,ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:ytrajc")

  allocate(ztrajc(npoints,ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:ztrajc")

  allocate(utrajc(npoints,ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:utrajc")

  allocate(vtrajc(npoints,ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:vtrajc")

  allocate(wtrajc(npoints,ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:wtrajc")

! allocate(xweight(npoints,ntrajcs),stat=istatus)
! CALL check_alloc_status(istatus, "arpstrajc:xweight")

! allocate(yweight(npoints,ntrajcs),stat=istatus)
! CALL check_alloc_status(istatus, "arpstrajc:yweight")

! allocate(zweight(npoints,ntrajcs),stat=istatus)
! CALL check_alloc_status(istatus, "arpstrajc:zweight")

! allocate(itrajc(npoints,ntrajcs),stat=istatus)
! CALL check_alloc_status(istatus, "arpstrajc:itrajc")

! allocate(jtrajc(npoints,ntrajcs),stat=istatus)
! CALL check_alloc_status(istatus, "arpstrajc:jtrajc")

! allocate(ktrajc(npoints,ntrajcs),stat=istatus)
! CALL check_alloc_status(istatus, "arpstrajc:ktrajc")


  allocate(ttrajc(npoints),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:ttrajc")

  ntrajc0_sub = ntrajc0*subinterval-(subinterval-1)

  DO k=1,ntimes
  DO j=1,ntrajcs
    xtrajc(ntrajc0_sub,j,k)=xtrajc0(j,k)
    ytrajc(ntrajc0_sub,j,k)=ytrajc0(j,k)
    ztrajc(ntrajc0_sub,j,k)=ztrajc0(j,k)
  ENDDO
  ENDDO

  IF( ireturn /= 0 ) THEN
    PRINT*,'Problem occured when trying to get dimensions from data.'
    PRINT*,'Program stopped.'
    STOP
  END IF

  WRITE(6,'(4(a,i5))') 'nx =',nx,', ny=',ny,', nz=',nz ,'nzsoil=',nzsoil

!-----------------------------------------------------------------------
!  Find start and end times for trajectories from data
!  for hdmpinopt =2.
!-----------------------------------------------------------------------

   IF( hdmpinopt == 1 ) then

     ttrajc(1) = tbgn_dmpin
     ttrajc(npoints) = tend_dmpin

   ELSE

     nfile = nhisfile
     CALL dtareadwind(nx,ny,nz,x,y,z,zp, trim(grdbasfn),trim(hisfile(nfile)),   &
          time,u,v,w, itmp)
     ttrajc(npoints) = time

     nfile = 1
     print*,'calling dtareadwind to read file ', trim(hisfile(nfile))
     CALL dtareadwind(nx,ny,nz,x,y,z,zp, trim(grdbasfn),trim(hisfile(nfile)),   &
          time,u,v,w, itmp)
     ttrajc(1) = time

   ENDIF

!-----------------------------------------------------------------------
!  Read history data at the starting (reference) time of trajectories.
!-----------------------------------------------------------------------

   nfile = ntrajc0
   nsub = nfile*subinterval-(subinterval-1)
   WRITE(6,'(/a,a,a)') ' Data set ', trim(hisfile(nfile)),' to be read.'
   CALL dtareadwind(nx,ny,nz,x,y,z,zp, trim(grdbasfn),trim(hisfile(nfile)),   &
        time,u,v,w, itmp)

   ttrajc0 = time
   time1 = time
   ttrajc(nsub)=time

   print*,'ntrajc0, ttrajc0=', ntrajc0, ttrajc0

   print*,'At time=', ttrajc0, 'xtrajc(1),ytrajc(1),ztrajc(1) =',  &
           xtrajc(nfile,1,1),ytrajc(nfile,1,1),ztrajc(nfile,1,1)

   print*,'Data at Time =', time, ' read.'
!write(0,*) 'before a3dmax0u'
   CALL a3dmax0(u,1,nx,1,nx,1,ny,1,ny-1,1,nz,1,nz-1, amax,amin)
   WRITE(6,'(1x,2(a,e13.6))') 'umin    = ', amin,',  umax    =',amax
!write(0,*) 'before a3dmax0v'
   CALL a3dmax0(v,1,nx,1,nx-1,1,ny,1,ny,1,nz,1,nz-1, amax,amin)
   WRITE(6,'(1x,2(a,e13.6))') 'vmin    = ', amin,',  vmax    =',amax
!write(0,*) 'before a3dmax0w'
   CALL a3dmax0(w,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz, amax,amin)
   WRITE(6,'(1x,2(a,e13.6))') 'wmin    = ', amin,',  wmax    =',amax


   print*,'Data at Time =', time, ' read.'

!write(6,*) 'here'
   DO i=1,nx-1
     xs(i) = 0.5*(x(i)+x(i+1))
   ENDDO

   DO j=1,ny-1
     ys(j) = 0.5*(y(j)+y(j+1))
   ENDDO

!write(6,*) 'ntimes = ',ntimes
!-----------------------------------------------------------------------
!  Open trajectory data file and write header information
!-----------------------------------------------------------------------

   IF( trajcfn_header == ' ' ) trajcfn_header = runname

   DO k=1,ntimes
   CALL cvttsnd( ttrajc0        , timsnd0, tmstrln0 )
   CALL cvttsnd( ttrajc(1)      , timsnd1, tmstrln1 )
   CALL cvttsnd( ttrajc(npoints), timsnd2, tmstrln2 )
   CALL cvttsnd( reftime(k),         timsnd3, tmstrln3 )

   trajc_fn = trim(trajcfn_header)//'.trajc_'//trim(timsnd1)//'-'//trim(timsnd2)//'_'//trim(timsnd3)
   CALL getunit( nunit(k))

   ltrajc_fn = len_trim(trajc_fn)
   CALL fnversn( trajc_fn, ltrajc_fn )
   OPEN(UNIT=nunit(k),FILE=trim(trajc_fn),STATUS='unknown',   &
        FORM='formatted',IOSTAT= istat )

   CALL getunit( nunit1(k))
   OPEN(UNIT=nunit1(k),FILE=trim(trajc_fn)//'.backward',STATUS='unknown',   &
        FORM='formatted',IOSTAT= istat )

   WRITE(nunit(k),'(a)') trim(trajcfn_header)
   WRITE(nunit(k),'(6e17.6)') x(2),x(nx-1),y(2),y(ny-1),z(2),z(nz-1)
   WRITE(nunit(k),'(3e17.6)') dx, dy, z(2)-z(1)

   WRITE(nunit(k),'(3e17.6)') ttrajc(1),ttrajc0,ttrajc(npoints)
   WRITE(nunit(k),'(i10)') npoints
   WRITE(nunit(k),'(i10)') ntrajcs

   WRITE(nunit1(k),'(a)') trim(trajcfn_header)
   WRITE(nunit1(k),'(6e17.6)') x(2),x(nx-1),y(2),y(ny-1),z(2),z(nz-1)
   WRITE(nunit1(k),'(3e17.6)') dx, dy, z(2)-z(1)

   WRITE(nunit1(k),'(3e17.6)') ttrajc(1),ttrajc0,ttrajc(npoints)
   WRITE(nunit1(k),'(i10)') npoints
   WRITE(nunit1(k),'(i10)') ntrajcs

   notinteg = 1
   tinc = 1.0 ! should not be used because notinteg = 1
   cur_level = ntrajc0_sub
   CALL calc_trajc(nx,ny,nz,xs,ys,zp,u,v,w,tinc,cur_level, &
           xtrajc(1,1,k),ytrajc(1,1,k),ztrajc(1,1,k),  &
           utrajc(1,1,k),vtrajc(1,1,k),wtrajc(1,1,k),      &
!          itrajc,jtrajc,ktrajc,xweight,yweight,zweight,   &
           npoints,ntrajcs, notinteg)

   j=ntrajc0_sub
   write(nunit1(k),'(4e17.6)') ttrajc(j)
   write(nunit1(k),'(i10)') ntrajcs
   write(nunit1(k),'(6e17.6)') ((xtrajc(j,i,k),ytrajc(j,i,k),ztrajc(j,i,k)),i=1,ntrajcs)

   ENDDO

!-----------------------------------------------------------------------
!  Calculate backward trajectories.
!-----------------------------------------------------------------------

   time1 = ttrajc(ntrajc0_sub)

   DO nfile = ntrajc0-1,1,-1

     WRITE(6,'(/a,a,a)') ' Data set ', trim(hisfile(nfile)),' to be read.'
     CALL dtareadwind(nx,ny,nz,x,y,z,zp, trim(grdbasfn),trim(hisfile(nfile)),   &
                      time,u,v,w, itmp)
     time2 = time

     print*,'Data at Time =', time, ' read.'

     CALL a3dmax0(u,1,nx,1,nx,1,ny,1,ny-1,1,nz,1,nz-1, amax,amin)
     WRITE(6,'(1x,2(a,e13.6))') 'umin    = ', amin,',  umax    =',amax
     CALL a3dmax0(v,1,nx,1,nx-1,1,ny,1,ny,1,nz,1,nz-1, amax,amin)
     WRITE(6,'(1x,2(a,e13.6))') 'vmin    = ', amin,',  vmax    =',amax
     CALL a3dmax0(w,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz, amax,amin)
     WRITE(6,'(1x,2(a,e13.6))') 'wmin    = ', amin,',  wmax    =',amax

     ! DTD: get winds at the current time as well, to allow
     ! for temporal interpolation. Store in tem1, tem2, tem3

     WRITE(6,'(/a,a,a)') ' Data set ', trim(hisfile(nfile+1)),' to be read.'
     CALL dtareadwind(nx,ny,nz,x,y,z,zp, trim(grdbasfn),trim(hisfile(nfile+1)),   &
          time,tem1,tem2,tem3, itmp)

     print*,'Data at Time =', time, ' read.'

     CALL a3dmax0(tem1,1,nx,1,nx,1,ny,1,ny-1,1,nz,1,nz-1, amax,amin)
     WRITE(6,'(1x,2(a,e13.6))') 'umin    = ', amin,',  umax    =',amax
     CALL a3dmax0(tem2,1,nx,1,nx-1,1,ny,1,ny,1,nz,1,nz-1, amax,amin)
     WRITE(6,'(1x,2(a,e13.6))') 'vmin    = ', amin,',  vmax    =',amax
     CALL a3dmax0(tem3,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz, amax,amin)
     WRITE(6,'(1x,2(a,e13.6))') 'wmin    = ', amin,',  wmax    =',amax

!write(0,*) 'here'

!-----------------------------------------------------------------------
!    Calculate trajectories
!-----------------------------------------------------------------------
!write(0,*) 'calling calc_trajc ',time2, time1
     tinc = (time2-time1)/subinterval
!     cur_level = nfile+1
     notinteg = 0

     DO ns = 1,subinterval
       cur_level = (nfile+1)*subinterval-(subinterval-1)-(ns-1)
       nsub = cur_level-1
       ttrajc(nsub)=time1+ns*tinc
       write(0,*) 'time1,time2,nsub,ttrajc(nsub)',time1,time2,nsub,ttrajc(nsub)

       !DTD: temporally interpolate winds to current step
       ! wind components in tem4,tem5,tem6
       DO k=1,nz
        DO j=1,ny
          DO i=1,nx
            tem4(i,j,k) = tem1(i,j,k)+(ttrajc(nsub)-time1)*(u(i,j,k)-tem1(i,j,k))/(time2-time1)
            tem5(i,j,k) = tem2(i,j,k)+(ttrajc(nsub)-time1)*(v(i,j,k)-tem2(i,j,k))/(time2-time1)
            tem6(i,j,k) = tem3(i,j,k)+(ttrajc(nsub)-time1)*(w(i,j,k)-tem3(i,j,k))/(time2-time1)
          END DO
        END DO
       END DO

       DO k=1,ntimes
!  write(0,*) 'calling calc_trajc'
       CALL calc_trajc(nx,ny,nz,xs,ys,zp,tem4,tem5,tem6,tinc,cur_level, &
            xtrajc(1,1,k),ytrajc(1,1,k),ztrajc(1,1,k),  &
            utrajc(1,1,k),vtrajc(1,1,k),wtrajc(1,1,k),      &
!           itrajc,jtrajc,ktrajc,xweight,yweight,zweight,   &
            npoints,ntrajcs, notinteg)

!       print*,'At time=', ttrajc(nsub), 'xtrajc(1),ytrajc(1) =', xtrajc(nsub,1,1),ytrajc(nsub,1,1)

!       time1 = time2

       j = nsub
       write(nunit1(k),'(4e17.6)') ttrajc(j)
       write(nunit1(k),'(i10)') ntrajcs
       write(nunit1(k),'(6e17.6)') ((xtrajc(j,i,k),ytrajc(j,i,k),ztrajc(j,i,k)),i=1,ntrajcs)
       ENDDO

     ENDDO ! subinterval

     time1 = time2

   ENDDO

   DO k=1,ntimes

!-----------------------------------------------------------------------
!  Write out the backward trajectory points into the combined trajectory file
!-----------------------------------------------------------------------

     do j=1,ntrajc0_sub
       write(nunit(k),'(4e17.6)') ttrajc(j)
       write(nunit(k),'(i10)') ntrajcs
       write(nunit(k),'(6e17.6)') ((xtrajc(j,i,k),ytrajc(j,i,k),ztrajc(j,i,k)),i=1,ntrajcs)
     enddo

     CLOSE(nunit1(k),status='delete') ! No need for this file anymore
     CALL retunit(nunit1(k))

   ENDDO

!-----------------------------------------------------------------------
!  Calculate forward trajectories
!  First re-read data at starting (reference) time if necessary.
!-----------------------------------------------------------------------

   time1 = ttrajc(ntrajc0_sub)

   DO nfile = ntrajc0+1,nhisfile,+1

     WRITE(6,'(/a,a,a)') ' Data set ', trim(hisfile(nfile)),' to be read.'
     CALL dtareadwind(nx,ny,nz,x,y,z,zp, trim(grdbasfn),trim(hisfile(nfile)),   &
          time,u,v,w, itmp)
     time2 = time
!     ttrajc(nfile)=time
     print*,'Data at Time =', time, ' read.'

     CALL a3dmax0(u,1,nx,1,nx,1,ny,1,ny-1,1,nz,1,nz-1, amax,amin)
     WRITE(6,'(1x,2(a,e13.6))') 'umin    = ', amin,',  umax    =',amax
     CALL a3dmax0(v,1,nx,1,nx-1,1,ny,1,ny,1,nz,1,nz-1, amax,amin)
     WRITE(6,'(1x,2(a,e13.6))') 'vmin    = ', amin,',  vmax    =',amax
     CALL a3dmax0(w,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz, amax,amin)
     WRITE(6,'(1x,2(a,e13.6))') 'wmin    = ', amin,',  wmax    =',amax

     ! DTD: get winds at the current time as well, to allow
     ! for temporal interpolation. Store in tem1, tem2, tem3

     WRITE(6,'(/a,a,a)') ' Data set ', trim(hisfile(nfile-1)),' to be read.'
     CALL dtareadwind(nx,ny,nz,x,y,z,zp, trim(grdbasfn),trim(hisfile(nfile-1)),   &
          time,tem1,tem2,tem3, itmp)

     print*,'Data at Time =', time, ' read.'

     CALL a3dmax0(tem1,1,nx,1,nx,1,ny,1,ny-1,1,nz,1,nz-1, amax,amin)
     WRITE(6,'(1x,2(a,e13.6))') 'umin    = ', amin,',  umax    =',amax
     CALL a3dmax0(tem2,1,nx,1,nx-1,1,ny,1,ny,1,nz,1,nz-1, amax,amin)
     WRITE(6,'(1x,2(a,e13.6))') 'vmin    = ', amin,',  vmax    =',amax
     CALL a3dmax0(tem3,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz, amax,amin)
     WRITE(6,'(1x,2(a,e13.6))') 'wmin    = ', amin,',  wmax    =',amax


!write(0,*) 'here'
!-----------------------------------------------------------------------
!    Calculate trajectories
!-----------------------------------------------------------------------
!write(0,*) 'calling calc_trajc ',time2, time1
     tinc = (time2-time1)/subinterval
!     cur_level = nfile-1
     notinteg = 0

!    print*,'cur_level =', cur_level
!    print*,'xtrajc(cur_level,:) =', xtrajc(cur_level,1), xtrajc(cur_level,2)

     DO ns=1,subinterval

       cur_level = (nfile-1)*subinterval-(subinterval-1)+(ns-1)
       nsub = cur_level+1
       ttrajc(nsub)=time1+ns*tinc

       !DTD: temporally interpolate winds to current step
       ! wind components in tem4,tem5,tem6
       DO k=1,nz
        DO j=1,ny
          DO i=1,nx
            tem4(i,j,k) = tem1(i,j,k)+(ttrajc(nsub)-time1)*(u(i,j,k)-tem1(i,j,k))/(time2-time1)
            tem5(i,j,k) = tem2(i,j,k)+(ttrajc(nsub)-time1)*(v(i,j,k)-tem2(i,j,k))/(time2-time1)
            tem6(i,j,k) = tem3(i,j,k)+(ttrajc(nsub)-time1)*(w(i,j,k)-tem3(i,j,k))/(time2-time1)
          END DO
        END DO
       END DO

       DO k=1,ntimes

!  write(0,*) 'calling calc_trajc'
       CALL calc_trajc(nx,ny,nz,xs,ys,zp,tem4,tem5,tem6,tinc,cur_level, &
            xtrajc(1,1,k),ytrajc(1,1,k),ztrajc(1,1,k),  &
            utrajc(1,1,k),vtrajc(1,1,k),wtrajc(1,1,k),      &
  !           itrajc,jtrajc,ktrajc,xweight,yweight,zweight,   &
            npoints,ntrajcs, notinteg)

  !     time1 = time2

!       print*,'At time=', time2, 'xtrajc(1),ytrajc(1),ztrajc(1) =',  &
!       xtrajc(nfile,1,k),ytrajc(nfile,1,k),ztrajc(nfile,1,k)

!-----------------------------------------------------------------------
!  Write out the forward trajectory points
!-----------------------------------------------------------------------
       j = nsub
       write(nunit(k),'(4e17.6)') ttrajc(j)
       write(nunit(k),'(i10)') ntrajcs
       write(nunit(k),'(6e17.6)') ((xtrajc(j,i,k),ytrajc(j,i,k),ztrajc(j,i,k)),i=1,ntrajcs)

       ENDDO

    ENDDO ! subinterval

    time1 = time2

   ENDDO


!-----------------------------------------------------------------------
!  Close the trajectory file
!-----------------------------------------------------------------------

   DO k=1,ntimes
     CLOSE(UNIT=nunit(k))
     CALL retunit( nunit(k))
   ENDDO

  STOP

100 WRITE(6,'(a)') 'Error reading NAMELIST file. Program ARPSTRAJC stopped.'
  STOP

END PROGRAM ARPSTRAJC
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE DTAREADWIND                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE dtareadwind(nx,ny,nz,x,y,z,zp, grdbasfn,datafn, time,u,v,w, itmp)

!
!-----------------------------------------------------------------------
!  Variable Declarations.
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  INTEGER :: hinfmt            ! The format of the history data dump
  CHARACTER(LEN=*) :: grdbasfn ! Name of the grid/base state array file
  CHARACTER(LEN=*) :: datafn   ! Name of the other time dependent data file

  REAL :: x(nx),y(ny),z(nz)
  REAL :: time

  REAL :: u(nx,ny,nz)
  REAL :: v(nx,ny,nz)
  REAL :: w(nx,ny,nz)
  REAL :: zp(nx,ny,nz)

  INTEGER (KIND=selected_int_kind(4)) :: itmp(nx,ny,nz) ! Temporary array

  INTEGER :: grdread,iread
  SAVE grdread

  INTEGER :: istat
  INTEGER :: ireturn           ! Return status indicator
  INTEGER :: grdbas            ! Wether this is a grid/base state
                               ! array dump
  INTEGER :: i,j,k
  LOGICAL :: fexist
  INTEGER :: is
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

  Print*,'In dtareadwind, grdbasfn =', trim(grdbasfn)
  Print*,'In dtareadwind, datafn =', trim(datafn)

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
        'Program stopped in DTAREADWIND location 1.'
    CALL arpsstop('arpsstop called from dtareadwind during base state read',1)

    200 CONTINUE
!
!-----------------------------------------------------------------------
!  Read grid and base state fields.
!-----------------------------------------------------------------------
!
    CALL hdfreadwind(nx,ny,nz,x,y,z,zp,grdbas,trim(grdbasfn), &
                     time,u,v,w, itmp)

    grdread = 1

  END IF
!
!-----------------------------------------------------------------------
!  Read data fields.
!-----------------------------------------------------------------------
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
     'Program stopped in DTAREADWIND location 2.'
  CALL arpsstop('arpsstop called from dtareadwind during base read-2',1)

  100 CONTINUE

  CALL hdfreadwind(nx,ny,nz,x,y,z,zp,grdbas,trim(datafn), &
                   time,u,v,w, itmp)

  RETURN
END SUBROUTINE dtareadwind

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

SUBROUTINE hdfreadwind(nx,ny,nz,x,y,z,zp,grdbas,filename,   &
           time,u,v,w, itmp)

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

  INTEGER :: grdbas
  CHARACTER (LEN=*) :: filename

  REAL :: x     (nx)           ! x coord.
  REAL :: y     (ny)           ! y coord.
  REAL :: z     (nz)           ! z coord.

  REAL :: time
  REAL :: u(nx,ny,nz)
  REAL :: v(nx,ny,nz)
  REAL :: w(nx,ny,nz)
  REAL :: zp(nx,ny,nz)

  INTEGER (KIND=selected_int_kind(4)) :: itmp(nx,ny,nz) ! Temporary array
  REAL :: hmax(nz), hmin(nz) ! Temporary array

  INTEGER :: ireturn

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

  CHARACTER (LEN=40) :: fmtver410,fmtver500
  INTEGER  :: intver,intver410,intver500

  PARAMETER (fmtver410='004.10 HDF4 Coded Data',intver410=410)
  PARAMETER (fmtver500='005.00 HDF4 Coded Data',intver500=500)

  CHARACTER (LEN=40) :: fmtverin
  CHARACTER (LEN=10) :: tmunit

!-----------------------------------------------------------------------
!  Misc. local variables
!-----------------------------------------------------------------------

  INTEGER :: lchanl
  PARAMETER (lchanl=6)      ! Channel number for formatted printing.

  INTEGER :: i,j,k,is
  INTEGER :: nxin,nyin,nzin

  INTEGER :: bgrdin,bbasin,bvarin,bicein,btrbin,btkein

  INTEGER :: istat, sd_id
  INTEGER :: varflg, istatus

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'        ! Grid parameters
  INCLUDE 'indtflg.inc'


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  Beginning of executable code...
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


  WRITE(*,*) 'HDFREAD: Reading HDF file: ', trim(filename)

!-----------------------------------------------------------------------
! Open file for reading
!-----------------------------------------------------------------------

  CALL hdfopen(filename,1,sd_id)
  IF (sd_id < 0) THEN
    WRITE (6,*) "HDFREAD: ERROR opening ", trim(filename)," for reading."
    GO TO 110
  END IF

  fmtverin = fmtver500

  WRITE(6,'(/1x,a,a/)') 'Incoming data format, fmtverin=',fmtverin

  CALL hdfrdc(sd_id,80,"runname",runname,istat)
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

  WRITE (6,*) " "

  CALL hdfrdc(sd_id,10,"tmunit",tmunit,istat)
  CALL hdfrdr(sd_id,"time",time,istat)

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
  END IF  ! grdin

!-----------------------------------------------------------------------
!  Read in flags for different data groups
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
    CALL hdfrdi(sd_id,"totflg",totin,istat)
    CALL hdfrdi(sd_id,"tkeflg",tkein,istat)

    print*,'done reading parameters'

  END IF

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

  IF( grdin == 1 .OR. grdbas == 1 ) THEN
    CALL hdfrd3d(sd_id,"zp",nx,ny,nz,zp,istat,itmp,hmax,hmin)
    IF (istat /= 0) GO TO 110
  ENDIF

!-----------------------------------------------------------------------
!  Read in base state fields
!-----------------------------------------------------------------------

! Print*,'start doing 3d arrays'

  IF( (basin == 1 .OR. grdbas == 1)  .and. hdmpfmt /= 0) THEN

!   CALL hdfrd3d(sd_id,"ubar",nx,ny,nz,ubar,istat,itmp,hmax,hmin)
!   IF (istat /= 0) GO TO 110

!   CALL hdfrd3d(sd_id,"vbar",nx,ny,nz,vbar,istat,itmp,hmax,hmin)
!   IF (istat /= 0) GO TO 110

!   CALL hdfrd3d(sd_id,"ptbar",nx,ny,nz,ptbar,istat,itmp,hmax,hmin)
!   IF (istat /= 0) GO TO 110

!   CALL hdfrd3d(sd_id,"pbar",nx,ny,nz,pbar,istat,itmp,hmax,hmin)
!   IF (istat /= 0) GO TO 110

!   IF( mstin == 1 ) THEN
!     CALL hdfrd3d(sd_id,"qvbar",nx,ny,nz,qvbar,istat,itmp,hmax,hmin)
!     IF (istat /= 0) GO TO 110
!   ELSE
!     qvbar= 0.0
!   END IF

  ENDIF

  landout = 0 ! skipping landout for know

  IF( grdbas == 1 ) GO TO 930

  IF( varin == 1 ) then

!-----------------------------------------------------------------------
!  Read in total values of variables from history dump
!-----------------------------------------------------------------------

    CALL hdfrd3d(sd_id,"u",nx,ny,nz,u,istat,itmp,hmax,hmin)
    IF (istat /= 0) GO TO 110

    CALL hdfrd3d(sd_id,"v",nx,ny,nz,v,istat,itmp,hmax,hmin)
    IF (istat /= 0) GO TO 110

    CALL hdfrd3d(sd_id,"w",nx,ny,nz,w,istat,itmp,hmax,hmin)
    IF (istat /= 0) GO TO 110

!   CALL hdfrd3d(sd_id,"pt",nx,ny,nz,pt,istat,itmp,hmax,hmin)
!   IF (istat /= 0) GO TO 110

!   CALL hdfrd3d(sd_id,"p",nx,ny,nz,p,istat,itmp,hmax,hmin)
!   IF (istat /= 0) GO TO 110

!   IF( mstin == 1 ) THEN

!     CALL hdfrd3d(sd_id,"qv",nx,ny,nz,qv,istat,itmp,hmax,hmin)
!     IF (istat /= 0) GO TO 110
!
!     CALL hdfrd3d(sd_id,"qc",nx,ny,nz,qc,istat,itmp,hmax,hmin)
!     IF (istat /= 0) GO TO 110
!
!     CALL hdfrd3d(sd_id,"qr",nx,ny,nz,qr,istat,itmp,hmax,hmin)
!     IF (istat /= 0) GO TO 110


!     IF( icein == 1 ) THEN
!       CALL hdfrd3d(sd_id,"qi",nx,ny,nz,qi,istat,itmp,hmax,hmin)
!       IF (istat /= 0) GO TO 110
!
!       CALL hdfrd3d(sd_id,"qs",nx,ny,nz,qs,istat,itmp,hmax,hmin)
!       IF (istat /= 0) GO TO 110
!
!       CALL hdfrd3d(sd_id,"qh",nx,ny,nz,qh,istat,itmp,hmax,hmin)
!       IF (istat /= 0) GO TO 110

!     END IF

!   END IF

  END IF

! IF( tkein == 1 ) THEN

!   CALL hdfrd3d(sd_id,"tke",nx,ny,nz,tke,istat,itmp,hmax,hmin)
!   IF (istat /= 0) GO TO 110

! END IF

! IF( trbin == 1 ) THEN

!   CALL hdfrd3d(sd_id,"kmh",nx,ny,nz,kmh,istat,itmp,hmax,hmin)
!   IF (istat /= 0) GO TO 110

!   CALL hdfrd3d(sd_id,"kmv",nx,ny,nz,,istat,itmp,hmax,hmin)
!   IF (istat /= 0) GO TO 110

! END IF

!
! Surface arrays skipped
!

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
  STOP
  ireturn=1

  130   CONTINUE

  CALL hdfclose(sd_id,istat)

  IF (ireturn == 0) THEN
    IF (istat == 0) THEN
      WRITE(6,'(/a/a)') &
      "HDFREADWIND: Successfully read ", trim(filename)
    ELSE
      WRITE(6,'(/a,i3,a/,a)') &
      "HDFREADWIND: ERROR (status=", istat, ") closing ", trim(filename)
      STOP
    END IF
  END IF

  RETURN
END SUBROUTINE hdfreadwind

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
  print*,'Program stopped in GET_GRIDINFO_FROM_HDF.'
  STOP
  ireturn=1

130 CONTINUE

!tmp  stat = sfendacc(sd_id)   ! is this necessary?
  CALL hdfclose(sd_id,stat)

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

SUBROUTINE intrpx3d(ain,nx,is,ie, ny,js,je, nz,ks,ke, wgtx,ix,          &
           aout,nx1,is1,ie1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Perform interpolation in the first dimension
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  4/1/1999.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: nx,is,ie, ny,js,je, nz,ks,ke
  INTEGER :: nx1,is1,ie1
  REAL :: ain (nx ,ny,nz)
  REAL :: aout(nx1,ny,nz)
  REAL :: wgtx(nx1)
  INTEGER :: ix(nx1)
  INTEGER :: i,j,k

    DO k=ks ,ke
      DO j=js ,je
        DO i=is1,ie1
          aout(i,j,k)=      wgtx(i) *ain(ix(i)  ,j,k)                 &
                      +(1.0-wgtx(i))*ain(ix(i)+1,j,k)
        END DO
      END DO
    END DO

  RETURN
END SUBROUTINE intrpx3d

SUBROUTINE intrpy3d(ain,nx,is,ie, ny,js,je, nz,ks,ke, wgty,jy,          &
           aout,ny1,js1,je1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Perform interpolation in the second dimension
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  4/1/1999.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,is,ie, ny,js,je, nz,ks,ke
  INTEGER :: ny1,js1,je1
  REAL :: ain (nx,ny ,nz)
  REAL :: aout(nx,ny1,nz)
  REAL :: wgty(ny1)
  INTEGER :: jy(ny1)
  INTEGER :: i,j,k

    DO k=ks ,ke
      DO j=js1,je1
        DO i=is ,ie
          aout(i,j,k)=      wgty(j) *ain(i,jy(j)  ,k)                 &
                      +(1.0-wgty(j))*ain(i,jy(j)+1,k)
        END DO
      END DO
    END DO

  RETURN
END SUBROUTINE intrpy3d

SUBROUTINE intrpxy3d(ain,nx,is,ie, ny,js,je, nz,ks,ke,                  &
           wgtx,ix,wgty,jy,                                             &
           aout,nx1,is1,ie1, ny1,js1,je1,                               &
           temx1yz)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Perform interpolation in the first and second dimensions
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  4/1/1999.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,is,ie, ny,js,je, nz,ks,ke
  INTEGER :: nx1,is1,ie1, ny1,js1,je1
  REAL :: ain (nx ,ny ,nz)
  REAL :: aout(nx1,ny1,ks:ke)
  REAL :: wgtx(nx1),wgty(ny1)
  INTEGER :: ix(nx1),jy(ny1)
  INTEGER :: k

  REAL :: temx1yz(nx1,ny)

  DO k=ks,ke

  CALL intrpx3d(ain(1,1,k),nx,is,ie, ny,js,je, 1,1,1,                &
                wgtx,ix,temx1yz,nx1,is1,ie1)

  CALL intrpy3d(temx1yz,nx1,is1,ie1, ny,js,je, 1,1,1,                &
                wgty,jy,aout(1,1,k),ny1,js1,je1)

  ENDDO

  CALL edgfill(aout,1,nx1,is1,ie1,1,ny1,js1,je1,ks,ke,ks,ke)

  RETURN
END SUBROUTINE intrpxy3d

SUBROUTINE calc_trajc(nx,ny,nz,xs,ys,zp,u,v,w,tinc,clevel, &
           xtrajc,ytrajc,ztrajc,utrajc,vtrajc,wtrajc,      &
!          itrajc,jtrajc,ktrajc,xweight,yweight,zweight,   &
           ntrajc_points,ntrajcs, notinteg)
!
!-----------------------------------------------------------------------
!  PURPOSE:
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  (4/08/2004)
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!  Variable Declarations:
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz
  REAL :: xs(nx),ys(ny),zp(nx,ny,nz)

  REAL :: u(nx,ny,nz)
  REAL :: v(nx,ny,nz)
  REAL :: w(nx,ny,nz)

  INTEGER :: ntrajc_points,ntrajcs

  REAL :: xtrajc(ntrajc_points,ntrajcs),ytrajc(ntrajc_points,ntrajcs)
  REAL :: ztrajc(ntrajc_points,ntrajcs)

  REAL :: utrajc(ntrajc_points,ntrajcs),vtrajc(ntrajc_points,ntrajcs),wtrajc(ntrajc_points,ntrajcs)
! REAL :: xweight(ntrajc_points,ntrajcs),yweight(ntrajc_points,ntrajcs),zweight(ntrajc_points,ntrajcs)
! INTEGER :: itrajc(ntrajc_points,ntrajcs),jtrajc(ntrajc_points,ntrajcs),ktrajc(ntrajc_points,ntrajcs)

  REAL :: tinc
  INTEGER :: clevel, flevel,get_ktrajc
  INTEGER :: niteration,niteration_max,ntrajc,level_inc

  REAL :: dx, dy, utrajc1,vtrajc1,wtrajc1
  REAL :: uy1z1,uy2z1,uy1z2,uy2z2,vy1z1,vy2z1,vy1z2,vy2z2,wy1z1,wy2z1,wy1z2,wy2z2
  REAL :: uz1,uz2,vz1,vz2,wz1,wz2
  REAL :: missing_value

  INTEGER :: notinteg   ! for the first time step where xtrajc,ytrajc,ztrajc
                        ! are known and no time integration of trajectory is performed.
  INTEGER :: niterat,itrajc1,jtrajc1,ktrajc1

  REAL :: xweight1,yweight1,zweight1

  REAL :: zs1d(nz),zy1(nz),zy2(nz) ! automatic work array

  INTEGER :: i,j,k, vlevel
!
!-----------------------------------------------------------------------
!  Include files:
!-----------------------------------------------------------------------
!

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  write(0,*) 'entering calc_trajc'
  dx = xs(2)-xs(1)
  dy = ys(2)-ys(1)

  IF( tinc == 0 ) then
    Print*,'tinc cannot be equal 0, job stopped in calc_trajc.'
    STOP
  ENDIF

! print*,'mark 1'

  IF( tinc > 0 ) level_inc = 1
  IF( tinc < 0 ) level_inc = -1

  flevel = clevel + level_inc
  if( notinteg == 1 ) flevel = clevel

  niteration_max = 7
  IF( notinteg == 1 ) THEN
    niterat = 1
  ELSE
    niterat = niteration_max
  ENDIF

! print*,'ntrajc_points,ntrajcs=', ntrajc_points,ntrajcs
! print*,'mark 2, flevel =', flevel, clevel, level_inc

  xtrajc(flevel,:) = xtrajc(clevel,:)
  ytrajc(flevel,:) = ytrajc(clevel,:)
  ztrajc(flevel,:) = ztrajc(clevel,:)

! print*,'xtrajc(clevel,:), xtrajc(flevel,:)=', xtrajc(clevel,:), xtrajc(flevel,:)
! print*,'ytrajc(clevel,:), ytrajc(flevel,:)=', ytrajc(clevel,:), ytrajc(flevel,:)
! print*,'ztrajc(clevel,:), ztrajc(flevel,:)=', ztrajc(clevel,:), ztrajc(flevel,:)

! print*,' I am here'

  missing_value = -99999.0

!write(0,*) 'inside calc_trajc ',ntrajcs
  DO ntrajc = 1,ntrajcs

!write(0,*) 'ntrajc = ',ntrajc
    IF( xtrajc(clevel,ntrajc) == missing_value .or.      &
        ytrajc(clevel,ntrajc) == missing_value .or.      &
        ztrajc(clevel,ntrajc) == missing_value ) cycle

    DO niteration = 1, niterat

!write(0,*) 'niteration = ',niteration, niterat
      IF( niteration == 1 .and. notinteg /= 1 ) THEN
        utrajc1 = utrajc(clevel,ntrajc)
        vtrajc1 = vtrajc(clevel,ntrajc)
        wtrajc1 = wtrajc(clevel,ntrajc)
      ELSE

      itrajc1 = MAX(1, MIN(nx-2, INT((xtrajc(flevel,ntrajc)-xs(1))/dx)+1 ))
      jtrajc1 = MAX(1, MIN(ny-2, INT((ytrajc(flevel,ntrajc)-ys(1))/dy)+1 ))

!write(0,*) 'xtrajc = ',flevel,xtrajc(flevel,ntrajc)

      xweight1 = (xtrajc(flevel,ntrajc)-xs(itrajc1))/dx
      yweight1 = (ytrajc(flevel,ntrajc)-ys(jtrajc1))/dy

      ! DTD test initial trajectory height

      IF(notinteg == 1) THEN
        print*,'Before: itrajc1,jtrajc1,ztrajc,zp(k=2)',itrajc1,jtrajc1,ztrajc(clevel,ntrajc),zp(itrajc1,jtrajc1,2)
        ! adjust initial trajectory height to be consistent with zp array, which is in ASL (user input should
        ! be interpreted as being AGL)
        ztrajc(clevel,ntrajc)=ztrajc(clevel,ntrajc)+zp(itrajc1,jtrajc1,2)
        print*,'After: itrajc1,jtrajc1,ztrajc,zp(k=1)',itrajc1,jtrajc1,ztrajc(clevel,ntrajc),zp(itrajc1,jtrajc1,2)

      END IF

      DO k=1,nz-1

        zy1(k) = (1.0-xweight1)*(zp(itrajc1  ,jtrajc1  ,k)+zp(itrajc1  ,jtrajc1  ,k+1))*0.5 &
                     +xweight1 *(zp(itrajc1+1,jtrajc1  ,k)+zp(itrajc1+1,jtrajc1  ,k+1))*0.5
        zy2(k) = (1.0-xweight1)*(zp(itrajc1  ,jtrajc1+1,k)+zp(itrajc1  ,jtrajc1+1,k+1))*0.5 &
                     +xweight1 *(zp(itrajc1+1,jtrajc1+1,k)+zp(itrajc1+1,jtrajc1+1,k+1))*0.5

        zs1d(k) = ( 1.0-yweight1 )*zy1(k) + yweight1*zy2(k)

      END DO

!write(0,*) ' here 2'
      ktrajc1 = get_ktrajc(ztrajc(flevel,ntrajc),zs1d,nz-1)

!     print*,'niteration =', niteration
!     print*,'xtrajc, ytrajc, ztrajc=', xtrajc(flevel,ntrajc),ytrajc(flevel,ntrajc),ztrajc(flevel,ntrajc)
!     print*,'itrajc,jtrajc  ,ktrajc =',itrajc1,jtrajc1  ,ktrajc1

      uy1z1 = (1.0-xweight1)*(u(itrajc1  ,jtrajc1  ,ktrajc1  )+u(itrajc1+1,jtrajc1  ,ktrajc1  ))*0.5 &
                   +xweight1*(u(itrajc1+1,jtrajc1  ,ktrajc1  )+u(itrajc1+2,jtrajc1  ,ktrajc1  ))*0.5
      uy2z1 = (1.0-xweight1)*(u(itrajc1  ,jtrajc1+1,ktrajc1  )+u(itrajc1+1,jtrajc1+1,ktrajc1  ))*0.5 &
                  +xweight1 *(u(itrajc1+1,jtrajc1+1,ktrajc1  )+u(itrajc1+2,jtrajc1+1,ktrajc1  ))*0.5
      uy1z2 = (1.0-xweight1)*(u(itrajc1  ,jtrajc1  ,ktrajc1+1)+u(itrajc1+1,jtrajc1  ,ktrajc1+1))*0.5 &
                  +xweight1* (u(itrajc1+1,jtrajc1  ,ktrajc1+1)+u(itrajc1+2,jtrajc1  ,ktrajc1+1))*0.5
      uy2z2 = (1.0-xweight1)*(u(itrajc1  ,jtrajc1+1,ktrajc1+1)+u(itrajc1+1,jtrajc1+1,ktrajc1+1))*0.5 &
                  +xweight1* (u(itrajc1+1,jtrajc1+1,ktrajc1+1)+u(itrajc1+2,jtrajc1+1,ktrajc1+1))*0.5

      vy1z1 = (1.0-xweight1)*(v(itrajc1  ,jtrajc1  ,ktrajc1  )+v(itrajc1  ,jtrajc1+1,ktrajc1  ))*0.5 &
                  +xweight1* (v(itrajc1+1,jtrajc1  ,ktrajc1  )+v(itrajc1+1,jtrajc1+1,ktrajc1  ))*0.5
      vy2z1 = (1.0-xweight1)*(v(itrajc1  ,jtrajc1+1,ktrajc1  )+v(itrajc1  ,jtrajc1+2,ktrajc1  ))*0.5 &
                  +xweight1* (v(itrajc1+1,jtrajc1+1,ktrajc1  )+v(itrajc1+1,jtrajc1+2,ktrajc1  ))*0.5
      vy1z2 = (1.0-xweight1)*(v(itrajc1  ,jtrajc1  ,ktrajc1+1)+v(itrajc1  ,jtrajc1+1,ktrajc1+1))*0.5 &
                  +xweight1* (v(itrajc1+1,jtrajc1  ,ktrajc1+1)+v(itrajc1+1,jtrajc1+1,ktrajc1+1))*0.5
      vy2z2 = (1.0-xweight1)*(v(itrajc1  ,jtrajc1+1,ktrajc1+1)+v(itrajc1  ,jtrajc1+2,ktrajc1+1))*0.5 &
                  +xweight1* (v(itrajc1+1,jtrajc1+1,ktrajc1+1)+v(itrajc1+1,jtrajc1+2,ktrajc1+1))*0.5

      wy1z1 = (1.0-xweight1)*(w(itrajc1  ,jtrajc1  ,ktrajc1  )+w(itrajc1  ,jtrajc1  ,ktrajc1+1))*0.5 &
                  +xweight1* (w(itrajc1+1,jtrajc1  ,ktrajc1  )+w(itrajc1+1,jtrajc1  ,ktrajc1+1))*0.5
      wy2z1 = (1.0-xweight1)*(w(itrajc1  ,jtrajc1+1,ktrajc1  )+w(itrajc1  ,jtrajc1+1,ktrajc1+1))*0.5 &
                  +xweight1* (w(itrajc1+1,jtrajc1+1,ktrajc1  )+w(itrajc1+1,jtrajc1+1,ktrajc1+1))*0.5
      wy1z2 = (1.0-xweight1)*(w(itrajc1  ,jtrajc1  ,ktrajc1+1)+w(itrajc1  ,jtrajc1  ,ktrajc1+2))*0.5 &
                  +xweight1* (w(itrajc1+1,jtrajc1  ,ktrajc1+1)+w(itrajc1+1,jtrajc1  ,ktrajc1+2))*0.5
      wy2z2 = (1.0-xweight1)*(w(itrajc1  ,jtrajc1+1,ktrajc1+1)+w(itrajc1  ,jtrajc1+1,ktrajc1+2))*0.5 &
                  +xweight1* (w(itrajc1+1,jtrajc1+1,ktrajc1+1)+w(itrajc1+1,jtrajc1+1,ktrajc1+2))*0.5

!write(0,*) 'wylz2,wy2z2',wy1z2,wy2z2, xweight1
!write(0,*) 'here 3 ',itrajc1, jtrajc1, ktrajc1, w(itrajc1  ,jtrajc1  ,ktrajc1+1)
      uz1 = ( 1.0-yweight1 )*uy1z1 + yweight1*uy2z1
      uz2 = ( 1.0-yweight1 )*uy1z2 + yweight1*uy2z2
      vz1 = ( 1.0-yweight1 )*vy1z1 + yweight1*vy2z1
      vz2 = ( 1.0-yweight1 )*vy1z2 + yweight1*vy2z2
      wz1 = ( 1.0-yweight1 )*wy1z1 + yweight1*wy2z1
      wz2 = ( 1.0-yweight1 )*wy1z2 + yweight1*wy2z2

!write(0,*) 'here 4 ',ktrajc1,zs1d(ktrajc1+1), zs1d(ktrajc1)
      zweight1 = (ztrajc(flevel,ntrajc)-zs1d(ktrajc1))/(zs1d(ktrajc1+1)-zs1d(ktrajc1))
!write(0,*) 'here 4.1'
      utrajc1 = (1.0-zweight1)*uz1+zweight1*uz2
!write(0,*) 'here 4.2'
      vtrajc1 = (1.0-zweight1)*vz1+zweight1*vz2
!write(0,*) 'here 4.3', wz1, wz2
      wtrajc1 = (1.0-zweight1)*wz1+zweight1*wz2
!write(0,*) 'here 5'

      ENDIF

!     print*,'xweight,yweight,zweight=', xweight1,yweight1,zweight1
!     print*,'tinc =', tinc


!write(0,*) 'notinteg = ',notinteg
      IF( notinteg /= 1 ) THEN

        IF( xtrajc(clevel,ntrajc) .ne. missing_value ) then
          xtrajc(flevel,ntrajc) = xtrajc(clevel,ntrajc) + tinc * 0.5*(utrajc(clevel,ntrajc)+utrajc1)
        ELSE
          xtrajc(flevel,ntrajc) = missing_value
        ENDIF

        IF( ytrajc(clevel,ntrajc) .ne. missing_value ) then
          ytrajc(flevel,ntrajc) = ytrajc(clevel,ntrajc) + tinc * 0.5*(vtrajc(clevel,ntrajc)+vtrajc1)
        ELSE
          ytrajc(flevel,ntrajc) = missing_value
        ENDIF

        ztrajc(flevel,ntrajc) = ztrajc(clevel,ntrajc) + tinc * 0.5*(wtrajc(clevel,ntrajc)+wtrajc1)
!       print*,'xtrajc, ytrajc, ztrajc=', xtrajc(flevel,ntrajc),ytrajc(flevel,ntrajc),ztrajc(flevel,ntrajc)
      ENDIF

!     print*,'clevel,flevel,niteration,ntrajc=',clevel,flevel,niteration,ntrajc
!     print*,'utrajc0,vtrajc0,wtrajc0=', utrajc(clevel,ntrajc),vtrajc(clevel,ntrajc),wtrajc(clevel,ntrajc)
!     print*,'utrajc1,vtrajc1,wtrajc1=', utrajc1,vtrajc1,wtrajc1
!     print*,'xtrajc, ytrajc, ztrajc=', xtrajc(flevel,ntrajc),ytrajc(flevel,ntrajc),ztrajc(flevel,ntrajc)

    ENDDO  ! niteration

!write(0,*) 'here 1'
    utrajc(flevel,ntrajc) = utrajc1
    vtrajc(flevel,ntrajc) = vtrajc1
    wtrajc(flevel,ntrajc) = wtrajc1

!   xweight(flevel,ntrajc) = xweight1  ! save the calculated weights in array for reuse
!   yweight(flevel,ntrajc) = yweight1  ! save the calculated weights in array for reuse
!   zweight(flevel,ntrajc) = zweight1  ! save the calculated weights in array for reuse

!   itrajc(flevel,ntrajc) = itrajc1  ! save the calculated weights in array for reuse
!   jtrajc(flevel,ntrajc) = jtrajc1  ! save the calculated weights in array for reuse
!   ktrajc(flevel,ntrajc) = ktrajc1  ! save the calculated weights in array for reuse

    IF( notinteg /= 1 ) THEN
      ztrajc(flevel,ntrajc) = max(ztrajc(flevel,ntrajc), (zs1d(1)+zs1d(2))*0.5+0.1*(zs1d(2)-zs1d(1)))
      ztrajc(flevel,ntrajc) = min(ztrajc(flevel,ntrajc), (zs1d(nz-1)+zs1d(nz-2))*0.5 )

      if( xtrajc(flevel,ntrajc) .lt. xs(1) .or. xtrajc(flevel,ntrajc).gt. xs(nx-1) ) then
        xtrajc(flevel,ntrajc) = missing_value
      endif

      if( ytrajc(flevel,ntrajc) .lt. ys(1) .or. ytrajc(flevel,ntrajc).gt. ys(ny-1) ) then
        ytrajc(flevel,ntrajc) = missing_value
      endif
    ENDIF

  ENDDO

  RETURN
END SUBROUTINE calc_trajc

INTEGER FUNCTION get_ktrajc(z, zs1d, nz)

  IMPLICIT NONE
  INTEGER :: nz, k
  REAL :: z, zs1d(nz)

  IF( z < zs1d(1) ) then
    get_ktrajc = 1
    RETURN
  ENDIF

  IF( z >= zs1d(nz) ) then
    get_ktrajc = nz-1
    RETURN
  ENDIF

  DO k=1,nz-1

    IF( z >= zs1d(k) .and. z < zs1d(k+1) ) then
      get_ktrajc = k
      EXIT
    endif

  ENDDO

  RETURN

END FUNCTION get_ktrajc


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

SUBROUTINE hdfreadvar(nx,ny,nz,filename, time, varname, var, itmp )

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

  CHARACTER (LEN=*) :: varname
  CHARACTER (LEN=*) :: filename

  REAL :: time
  REAL :: var(nx,ny,nz)

  INTEGER (KIND=selected_int_kind(4)) :: itmp(nx,ny,nz) ! Temporary array
  REAL :: hmax(nz), hmin(nz) ! Temporary array

  INTEGER :: ireturn

  INTEGER :: i,j,k,is
  INTEGER :: nxin,nyin,nzin

  INTEGER :: istat, sd_id
  INTEGER :: varflg, istatus

  REAL :: timein

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'        ! Grid parameters
  INCLUDE 'indtflg.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  Beginning of executable code...
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  WRITE(*,*) 'HDFREAD: Reading HDF file: ', trim(filename)

!-----------------------------------------------------------------------
! Open file for reading
!-----------------------------------------------------------------------

  CALL hdfopen(filename,1,sd_id)
  IF (sd_id < 0) THEN
    WRITE (6,*) "HDFREAD: ERROR opening ", trim(filename)," for reading."
    GO TO 110
  END IF

  CALL hdfrdr(sd_id,"time",timein,istat)

  IF( timein /= time ) then
    print*,'Warning: time in data does not match time passed into READHDFVAR.'
    Print*,'time in program =', time, ', time in data =', timein
    print*,'time in program reset to ', timein
    time = timein
  ENDIF

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
!  Read in total values of variables from history dump
!-----------------------------------------------------------------------

    CALL hdfrd3d(sd_id,varname,nx,ny,nz,var,istat,itmp,hmax,hmin)
    IF (istat /= 0) GO TO 110

!-----------------------------------------------------------------------
!  Friendly exit message
!-----------------------------------------------------------------------

  930   CONTINUE

  WRITE(6,'(/a,a,a,F8.1,a/)')                                               &
  ' Variable ', varname,' at time=', time/60,' (min) were successfully read.'

  ireturn = 0

  GO TO 130

!-----------------------------------------------------------------------
!  Error during read
!-----------------------------------------------------------------------

  110   CONTINUE
  WRITE(6,'(/a/)') ' Error reading data in HDFREADVAR'
!  STOP
  ireturn=1
  RETURN

  130   CONTINUE

  CALL hdfclose(sd_id,istat)

  IF (ireturn == 0) THEN
    IF (istat == 0) THEN
      WRITE(6,'(/a/a)') &
      "HDFREADVAR: Successfully read ", trim(filename)
    ELSE
      WRITE(6,'(/a,i3,a/,a)') &
      "HDFREADVAR: ERROR (status=", istat, ") closing ", trim(filename)
    END IF
  END IF

  RETURN
END SUBROUTINE hdfreadvar
