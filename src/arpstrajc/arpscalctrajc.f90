PROGRAM arpscalctrajc
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 PROGRAM ARPSCALCTRAJC                ######
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
!  09/05/2012 (Y. Wang)
!  Merged Brett J. Roberts', Youngsun Jung's and Dan Dawson's works
!  together.
!
!  Cleared code and document. Reorganized the code structure.
!
!  Added run-time options to replace the hard-coded options.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  USE arps_precision
  IMPLICIT NONE

  INTEGER, PARAMETER :: nmax_times=20
  INTEGER, PARAMETER :: nscalar_max = 20

  INTEGER :: nx,ny,nz,nzsoil   ! Grid dimensions.
  INTEGER :: nstyps            ! Maximum number of soil types.

  INTEGER :: hinfmt
  CHARACTER (LEN=256) :: grdbasfn,hisfile

  REAL(P),  ALLOCATABLE :: x (:),y (:),z (:)
  REAL(P),  ALLOCATABLE :: xs (:), ys (:), zs1d(:) ! x,y coord for scalar points

  REAL(P),  ALLOCATABLE :: zp(:,:,:), ptbar(:,:,:), pbar(:,:,:),qvbar(:,:,:),rhobar(:,:,:)

  REAL(P),  ALLOCATABLE :: xtrajc(:,:,:),ytrajc(:,:,:),ztrajc(:,:,:),ttrajc(:)

  REAL(P),  ALLOCATABLE :: xtrajc1(:,:),ytrajc1(:,:),ztrajc1(:,:)

  REAL(P),  ALLOCATABLE :: utrajc(:,:),vtrajc(:,:),wtrajc(:,:)
  REAL(P),  ALLOCATABLE :: pt_trajc(:,:),p_trajc(:,:),qv_trajc(:,:),q_trajc(:,:),ptprt_trajc(:,:)
  REAL(P),  ALLOCATABLE :: qvcalc_trajc(:,:),qvcalcmp_trajc(:,:),qvcalcmix_trajc(:,:),qvmix_trajc(:,:)
  !REAL(P),  ALLOCATABLE :: xweight(:,:),yweight(:,:),zweight(:,:)
  !INTEGER,  ALLOCATABLE :: itrajc(:,:),jtrajc(:,:),ktrajc(:,:)

  REAL(P),  ALLOCATABLE :: vortx_trajc(:,:),vorty_trajc(:,:),vortz_trajc(:,:), &
                           vorts_trajc(:,:),vortc_trajc(:,:)
  REAL(P),  ALLOCATABLE :: vortx_stch_trajc(:,:),vorty_stch_trajc(:,:),vortz_stch_trajc(:,:)
  REAL(P),  ALLOCATABLE :: vortx_tilt_trajc(:,:),vorty_tilt_trajc(:,:),vortz_tilt_trajc(:,:)
  REAL(P),  ALLOCATABLE :: vortx_gen_trajc(:,:),vorty_gen_trajc(:,:)
  REAL(P),  ALLOCATABLE :: vortx_mix_trajc(:,:),vorty_mix_trajc(:,:),vortz_mix_trajc(:,:)
  REAL(P),  ALLOCATABLE :: vorts_stch_trajc(:,:),vortc_stch_trajc(:,:)
  REAL(P),  ALLOCATABLE :: vorts_tilt_trajc(:,:),vortc_tilt_trajc(:,:)
  REAL(P),  ALLOCATABLE :: vorts_gen_trajc(:,:),vortc_gen_trajc(:,:)
  REAL(P),  ALLOCATABLE :: vorts_mix_trajc(:,:),vortc_mix_trajc(:,:)
  REAL(P),  ALLOCATABLE :: vorts_exhg_trajc(:,:),vortc_exhg_trajc(:,:)

  REAL(P),  ALLOCATABLE :: psi_trajc(:,:),vh_trajc(:,:)
  REAL(P),  ALLOCATABLE :: buoy_trajc(:,:),buoyq_trajc(:,:)
  REAL(P),  ALLOCATABLE :: frcs_trajc(:,:),pprt_trajc(:,:)
  REAL(P),  ALLOCATABLE :: upgrad_trajc(:,:),vpgrad_trajc(:,:),wpgrad_trajc(:,:)

  REAL(P),  ALLOCATABLE :: mpcool_trajc(:,:,:),mpheat_trajc(:,:,:), mptrate_trajc(:,:)
  REAL(P),  ALLOCATABLE :: mpcoolproctot_trajc(:,:,:),mpheatproctot_trajc(:,:,:)
  REAL(P),  ALLOCATABLE :: mpcoolrate_trajc(:,:,:),mpheatrate_trajc(:,:,:)
  REAL(P),  ALLOCATABLE :: mpcooltot_trajc(:,:,:),mpheattot_trajc(:,:,:)
  REAL(P),  ALLOCATABLE :: temp_trajc(:,:),ptcalc_trajc(:,:),ptcalcmix_trajc(:,:),ptmix_trajc(:,:)
  REAL(P),  ALLOCATABLE :: ptcalcmp_trajc(:,:),ppi_trajc(:,:)
  REAL(P),  ALLOCATABLE :: qscalar(:,:,:,:), qscalar_trajc(:,:,:)
  REAL(P),  ALLOCATABLE :: tke_trajc(:,:), kmh_trajc(:,:), kmv_trajc(:,:), &
                           pte_trajc(:,:), ptw_trajc(:,:), mse_trajc(:,:)
  REAL(P),  ALLOCATABLE :: ptecalc_trajc(:,:), ptecalcmp_trajc(:,:),ptecalcmix_trajc(:,:)
  REAL(P),  ALLOCATABLE :: Dm_trajc(:,:,:),rho_trajc(:,:)
  REAL(DP), ALLOCATABLE :: N0_trajc(:,:,:),alpha_trajc(:,:,:),N0eff_trajc(:,:,:)

  REAL(P),  ALLOCATABLE :: xtrajcavg(:,:),ytrajcavg(:,:),ztrajcavg(:,:)
  REAL(P),  ALLOCATABLE :: utrajcavg(:,:),vtrajcavg(:,:),wtrajcavg(:,:)
  REAL(P),  ALLOCATABLE :: pt_trajcavg(:,:),ptcalc_trajcavg(:,:)
  REAL(P),  ALLOCATABLE :: ptcalcmp_trajcavg(:,:),ptcalcmix_trajcavg(:,:)
  REAL(P),  ALLOCATABLE :: qv_trajcavg(:,:),qvcalc_trajcavg(:,:)
  REAL(P),  ALLOCATABLE :: qvcalcmp_trajcavg(:,:),qvcalcmix_trajcavg(:,:)
  REAL(P),  ALLOCATABLE :: pte_trajcavg(:,:),ptecalc_trajcavg(:,:)
  REAL(P),  ALLOCATABLE :: ptecalcmp_trajcavg(:,:),ptecalcmix_trajcavg(:,:)
  REAL(P),  ALLOCATABLE :: ntrajcsavg(:,:)
  LOGICAL,  ALLOCATABLE :: trajc_missing(:,:,:)
  REAL(P),  ALLOCATABLE :: qscalar_trajcavg(:,:,:)

!-----------------------------------------------------------------------
!
! Working arrays
!
!-----------------------------------------------------------------------

  REAL, ALLOCATABLE :: tem1 (:,:,:), tem2(:,:,:), tem3(:,:,:), tem4(:,:,:)
  REAL, ALLOCATABLE :: tem5 (:,:,:), tem6(:,:,:), tem7(:,:,:)

  INTEGER (KIND=INT16), ALLOCATABLE :: itmp(:,:,:) ! Temporary array

  REAL, ALLOCATABLE :: tem1_trajc(:,:), tem2_trajc(:,:), tem3_trajc(:,:)  ! Temporary arrays along trajectories

  INTEGER :: istatus, cur_level

  REAL :: cx(6)     ! (PI/6)*rho_qx

  REAL    :: xlow, xhigh, ylow, yhigh, zlow, zhigh, pttem
  INTEGER :: nstart, nend, npnt, ntrj, ninc
  REAL    :: tzero, tend
!
!-----------------------------------------------------------------------
!  Include files:
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'phycst.inc'

!-----------------------------------------------------------------------
!
! NAMELIST Definitions
!
!-----------------------------------------------------------------------

  CHARACTER (LEN=256) :: hdmpftrailer
  CHARACTER (LEN=256) :: hdmpfheader

  REAL :: tstart_calc, tend_calc, tinc_data, tinc_calc
  REAL :: reftime(nmax_times)
  CHARACTER(LEN=256) :: trajc_fn_in(nmax_times), trajc_fn_out(10)

  INTEGER :: ntimes

  NAMELIST /input/ hdmpfheader,hdmpftrailer,mphyopt,dirname,            &
                   tstart_calc,tend_calc,tinc_data,tinc_calc,           &
                   reftime,trajc_fn_in,ntimes

  INTEGER :: do_diagnostics, have_mixing_terms
  INTEGER :: avg_trajc_flag,therm_intg_flag, DSD_flag
  INTEGER :: mprate_flag, misc_fields_flag

  NAMELIST /options/ do_diagnostics, have_mixing_terms,                 &
                     avg_trajc_flag,therm_intg_flag, DSD_flag,          &
                     mprate_flag, misc_fields_flag

  REAL :: Ntcfix,N0rfix,N0sfix,N0gfix,N0hfix
  REAL :: alpharfix,alphaifix,alphasfix,alphagfix,alphahfix
  REAL :: rhor,rhoi,rhos,rhog,rhoh
  NAMELIST /microph_param/ Ntcfix,N0rfix,N0sfix,N0gfix,N0hfix,          &
                     alpharfix,alphaifix,alphasfix,alphagfix,alphahfix, &
                     rhor,rhoi,rhos,rhog,rhoh


!-----------------------------------------------------------------------
!
! Misc. local working variables
!
!-----------------------------------------------------------------------

  INTEGER :: i,j,k,ii,kk,n,ireturn
  INTEGER :: nq
  INTEGER :: dirstrpos

  CHARACTER(LEN=128) :: fmtstr

  REAL    :: time, time1, time2, amin, amax, degree, degree2radian, theta
  LOGICAL :: iexist

  REAL    :: u1_changed, utem,vtem,wtem, factor
  INTEGER, ALLOCATABLE :: npnt_trajc_start(:,:)
  INTEGER, ALLOCATABLE :: npnt_trajc_start_flag(:,:)

  INTEGER :: ntrajcs, npoints, npoints_in

  REAL    :: ttrajc0_in
  INTEGER :: ntrajcs_in, ntrajcs_each_j

  CHARACTER (LEN=80) :: timsnd
  INTEGER :: tmstrln, nunit(nmax_times,10),ltrajc_fn(10),istat,tmstrln3
  INTEGER :: kl, counter

  REAL    :: missing_value, pi
  REAL    :: temp1, temp2, temp3
  REAL    :: del_x, del_y, psi_prev, psi_next, del_psi

  INTEGER :: tmstrln0,tmstrln1,tmstrln2, idot
  CHARACTER (LEN=6) :: timsnd0,timsnd1,timsnd2

  CHARACTER (LEN=80) :: varunits(17)
  CHARACTER (LEN=80) :: varname(17)
  CHARACTER (LEN=6)  :: varid
  CHARACTER (LEN=6)  :: varid_mpheat(17)
  CHARACTER (LEN=6)  :: varid_mptrate
  CHARACTER (LEN=6)  :: varid_mpcool(10)
  CHARACTER (LEN=6)  :: varid_qscalar(17)

  CHARACTER (LEN=256) :: mixdirname,mixrunname

!-----------------------------------------------------------------------
!
! Working constants and functions
!
!-----------------------------------------------------------------------

  REAL, PARAMETER :: LV = 2.501e6     ! Latent heat of vaporization/condensation (J/kg)
  REAL, PARAMETER :: LF = 3.34e5       ! Latent heat of melting/freezing (J/kg)
  REAL, PARAMETER :: LS = LV+LF       ! Latent heat of sublimation/deposition (J/kg)
!  REAL, PARAMETER :: CP = 1005.46     ! Specific heat at constant pressure for dry air (J/kg/K)
!  REAL, PARAMETER :: RD = 287.05       ! Gas constant for dry air (J/kg/K)
!  REAL, PARAMETER :: p0 = 1.0e5       ! Reference surface pressure (Pa)
  REAL, PARAMETER :: p0inv = 1.0/p0
!  REAL, PARAMETER :: rddcp = RD/CP

  REAL :: f_pt2pte

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  CALL mpinit_var

  missing_value = -99999.0

!-----------------------------------------------------------------------
!
! Read run-time parameters
!
!-----------------------------------------------------------------------

  hinfmt = 3
  reftime = 0.0

  !
  ! input
  !
  READ(5,input,ERR=100)
  WRITE(6,'(a)')'Namelist input was successfully read.'
  WRITE(6,input)

  !
  ! options
  !
  do_diagnostics = 10
  have_mixing_terms = 0

  ! Set flags to turn on/off various calculations
  avg_trajc_flag =     0      ! Perform averaging on multiple trajectories
  therm_intg_flag =    0      ! Integrate thermodynamic energy equation
                              ! along trajectories (requires certain special
                              ! arrays containing mixing and microphysical
                              ! heating/cooling terms)
  mprate_flag =        0      ! Calculate various microphysics rates along trajectories
                              ! Note, this is automatically turned on if therm_intg_flag = 1
  DSD_flag =           0      ! Calculate various DSD parameters along trajectories

  misc_fields_flag =   0      ! Calculate various miscellaneous fields along trajectories

  READ(5,options,ERR=100)
  WRITE(6,'(a)')'Namelist options was successfully read.'
  WRITE(6,options)

  IF(therm_intg_flag == 1) mprate_flag = 1

  !
  ! microph_param
  !
  READ(5,microph_param,ERR=100)
  WRITE(6,'(a)')'Namelist microph_param was successfully read.'
  WRITE(6,microph_param)

!-----------------------------------------------------------------------
!
! Read trajectory files to initialize the trajectories
!
!-----------------------------------------------------------------------

  DO k=1,ntimes

    CALL getunit (nunit(k,1))

    OPEN(UNIT=nunit(k,1),FILE=trim(trajc_fn_in(k)),STATUS='old',   &
          FORM='formatted',IOSTAT= istat )

    READ(nunit(k,1),'(a)') runname
    READ(nunit(k,1),'(6e17.6)') xlow, xhigh, ylow, yhigh, zlow, zhigh

    WRITE(6,'(6e17.6)') xlow, xhigh, ylow, yhigh, zlow, zhigh
    READ(nunit(k,1),'(3e17.6)') dx, dy, dz
    WRITE(6,'(3e17.6)') dx, dy, dz

    READ(nunit(k,1),'(3e17.6)') tstart, tzero, tend
    READ(nunit(k,1),'(i10)') npoints
    READ(nunit(k,1),'(i10)') ntrajcs

  END DO

  ALLOCATE(xtrajc(npoints,ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:xtrajc")

  ALLOCATE(ytrajc(npoints,ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:ytrajc")

  ALLOCATE(ztrajc(npoints,ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:ztrajc")

  ALLOCATE(ttrajc(npoints),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:ttrajc")

  ALLOCATE(xtrajc1(ntrajcs,ntimes),stat=istatus)
  ALLOCATE(ytrajc1(ntrajcs,ntimes),stat=istatus)
  ALLOCATE(ztrajc1(ntrajcs,ntimes),stat=istatus)

  IF(avg_trajc_flag == 1) THEN
    ALLOCATE(xtrajcavg(npoints,ntimes),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:xtrajcavg")
    xtrajcavg=0.0
    ALLOCATE(ytrajcavg(npoints,ntimes),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:ytrajcavg")
    ytrajcavg=0.0
    ALLOCATE(ztrajcavg(npoints,ntimes),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:ztrajcavg")
    ztrajcavg=0.0
    ALLOCATE(utrajcavg(npoints,ntimes),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:utrajcavg")
    utrajcavg=0.0
    ALLOCATE(vtrajcavg(npoints,ntimes),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:vtrajcavg")
    vtrajcavg=0.0
    ALLOCATE(wtrajcavg(npoints,ntimes),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:wtrajcavg")
    wtrajcavg=0.0
    ALLOCATE(pt_trajcavg(npoints,ntimes),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:pt_trajcavg")
    pt_trajcavg=0.0
    ALLOCATE(ptcalc_trajcavg(npoints,ntimes),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:ptcalc_trajcavg")
    ptcalc_trajcavg=0.0
    ALLOCATE(ptcalcmp_trajcavg(npoints,ntimes),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:ptcalcmp_trajcavg")
    ptcalcmp_trajcavg=0.0
    ALLOCATE(ptcalcmix_trajcavg(npoints,ntimes),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:ptcalcmix_trajcavg")
    ptcalcmix_trajcavg=0.0
    ALLOCATE(qv_trajcavg(npoints,ntimes),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:qv_trajcavg")
    qv_trajcavg=0.0
    ALLOCATE(qvcalc_trajcavg(npoints,ntimes),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:qvcalc_trajcavg")
    qvcalc_trajcavg=0.0
    ALLOCATE(qvcalcmp_trajcavg(npoints,ntimes),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:qvcalcmp_trajcavg")
    qvcalcmp_trajcavg=0.0
    ALLOCATE(qvcalcmix_trajcavg(npoints,ntimes),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:qvcalcmix_trajcavg")
    qvcalcmix_trajcavg=0.0
    ALLOCATE(pte_trajcavg(npoints,ntimes),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:pte_trajcavg")
    pte_trajcavg=0.0
    ALLOCATE(ptecalc_trajcavg(npoints,ntimes),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:ptecalc_trajcavg")
    ptecalc_trajcavg=0.0
    ALLOCATE(ptecalcmp_trajcavg(npoints,ntimes),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:ptecalcmp_trajcavg")
    ptecalcmp_trajcavg=0.0
    ALLOCATE(ptecalcmix_trajcavg(npoints,ntimes),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:ptecalcmix_trajcavg")
    ptecalcmix_trajcavg=0.0
    ALLOCATE(ntrajcsavg(npoints,ntimes),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:ntrajcsavg")
    ntrajcsavg=0.0
    ALLOCATE(qscalar_trajcavg(npoints,ntimes,nscalar_max),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:qscalar_trajcavg")
    qscalar_trajcavg = 0.0
  END IF

  ALLOCATE(trajc_missing(npoints,ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:trajc_missing")
  trajc_missing = .false.
  ALLOCATE(npnt_trajc_start(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:npnt_trajc_start")
  npnt_trajc_start = 0
  ALLOCATE(npnt_trajc_start_flag(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:npnt_trajc_start_flag")
  npnt_trajc_start_flag = 0

  DO k=1,ntimes
  DO j=1,npoints
    READ(nunit(k,1),'(4e17.6)',err=115,end=115) ttrajc(j)
    READ(nunit(k,1),'(i10)',err=115,end=115) ntrajcs_in
    IF( ntrajcs_in /= ntrajcs ) then
      print*,'ntrajcs read in .ne. ntrajcs in program.'
      print*,'Job stopped'
      STOP
    ENDIF
    READ(nunit(k,1),'(6e17.6)',err=115,end=115) ((xtrajc(j,i,k),ytrajc(j,i,k),ztrajc(j,i,k)),i=1,ntrajcs)
  ENDDO

  CLOSE(UNIT=nunit(k,1))
  CALL retunit(nunit(k,1))
  ENDDO

  npoints_in = npoints

  DO j=1,npoints
    DO k=1,ntimes
      DO i=1,ntrajcs
        IF(xtrajc(j,i,k) /= missing_value .and. ytrajc(j,i,k) /= missing_value .and. ztrajc(j,i,k) /= missing_value) THEN
          IF(avg_trajc_flag == 1) THEN
            ntrajcsavg(j,k)=ntrajcsavg(j,k)+1
          END IF
          trajc_missing(j,i,k) = .false.
        ELSE
          trajc_missing(j,i,k) = .true.
        END IF
        IF(trajc_missing(j,i,k) == .false. .and. npnt_trajc_start_flag(i,k) == 0) THEN
          npnt_trajc_start(i,k) = j
          npnt_trajc_start_flag(i,k) = 1
        END IF
      END DO
    END DO
  END DO
  ! Calculate an "average" trajectory by averaging the x, y, and z coordinates at
  ! each valid trajectory point.
  IF(avg_trajc_flag == 1) THEN
    DO j=1,npoints
      DO k=1,ntimes
        DO i=1,ntrajcs
          IF( .not. trajc_missing(j,i,k)) THEN
            xtrajcavg(j,k)=xtrajcavg(j,k)+xtrajc(j,i,k)/ntrajcsavg(j,k)
            ytrajcavg(j,k)=ytrajcavg(j,k)+ytrajc(j,i,k)/ntrajcsavg(j,k)
            ztrajcavg(j,k)=ztrajcavg(j,k)+ztrajc(j,i,k)/ntrajcsavg(j,k)
          END IF
        END DO
      END DO
    END DO
  END IF

  GOTO 125
  115 continue

  npoints_in = max(1,j-1)

  125 continue

  print*,'done reading trajectory data'

! IF( tstart_calc /= tend_calc ) then
    nstart = npoints_in
    DO j=1,npoints_in
      IF( ttrajc(j) >= tstart_calc ) then
        nstart = j
        exit
      END IF
    END DO

    nend   = 1
    DO j=npoints_in,1,-1
      IF( ttrajc(j) <= tend_calc ) then
        nend = j
        exit
      END IF
    END DO
! ELSE
!   nstart = 1
!   nend   = npoints
! END IF


  PRINT *,'ttrajc(nstart),ttrajc(nend), tinc_calc=',                    &
          ttrajc(nstart),ttrajc(nend), tinc_calc

  ninc = max(1, nint( tinc_calc/(ttrajc(nstart+1)-ttrajc(nstart)) ))

  PRINT *,'tstart_calc, tend_calc, nstart, nend =',                     &
          tstart_calc, tend_calc, nstart, nend

  grdbasfn = trim(hdmpfheader)//'.hdfgrdbas'//trim(hdmpftrailer)

  CALL get_dims_from_data(hinfmt,trim(grdbasfn),                        &
                          nx,ny,nz,nzsoil,nstyps, ireturn)

  PRINT *,'nx,ny,nz of input data were ', nx,ny,nz

  ALLOCATE(x(nx ),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:x")
  x = 0.0

  ALLOCATE(y(ny ),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:y")
  y = 0.0

  ALLOCATE(z(nz ),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:z")
  z = 0.0

  ALLOCATE(xs(nx),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:xs")
  xs = 0.0

  ALLOCATE(ys(ny),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:ys")
  ys = 0.0

  ALLOCATE(zs1d(nz),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:zs")
  zs1d = 0.0

  ALLOCATE(zp(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:zp")
  zp=0.0

  ALLOCATE(ptbar(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:ptbar")
  ptbar=0.0

  ALLOCATE(pbar(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:pbar")
  pbar=0.0

  ALLOCATE(qvbar(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:qvbar")
  qvbar=0.0

  ALLOCATE(rhobar(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:rhobar")
  rhobar=0.0

  ALLOCATE(tem1 (nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:tem1 ")
  tem1 =0.0

  ALLOCATE(tem2 (nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:tem2 ")
  tem2 =0.0

  ALLOCATE(tem3 (nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:tem3 ")
  tem3 =0.0

  ALLOCATE(tem4 (nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:tem4 ")
  tem4 =0.0

  ALLOCATE(itmp(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:itmp")
  itmp=0

  allocate(tem5 (nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:tem5 ")
  tem5 =0.0

  allocate(tem6 (nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:tem6 ")
  tem6 =0.0

  allocate(tem7 (nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:tem7 ")
  tem7 =0.0

  allocate(psi_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:psi_trajc")
  psi_trajc = missing_value

  allocate(vh_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:vh_trajc")
  vh_trajc = missing_value

  ALLOCATE(vortx_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:vortx_trajc")
  vortx_trajc = missing_value

  ALLOCATE(vorty_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:vorty_trajc")
  vorty_trajc = missing_value

  ALLOCATE(vortz_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:vortz_trajc")
  vortz_trajc = missing_value

  ALLOCATE(vorts_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:vorts_trajc")
  vorts_trajc= missing_value

  ALLOCATE(vortc_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:vortc_trajc")
  vortc_trajc= missing_value

  allocate(vortx_stch_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:vortx_stch_trajc")
  vortx_stch_trajc = missing_value

  allocate(vorty_stch_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:vorty_stch_trajc")
  vorty_stch_trajc = missing_value

  allocate(vortz_stch_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:vortz_stch_trajc")
  vortz_stch_trajc= missing_value

  allocate(vorts_stch_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:vorts_stch_trajc")
  vorts_stch_trajc= missing_value

  allocate(vortc_stch_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:vortc_stch_trajc")
  vortc_stch_trajc= missing_value

  allocate(vortx_tilt_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:vortx_tilt_trajc")
  vortx_tilt_trajc = missing_value

  allocate(vorty_tilt_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:vorty_tilt_trajc")
  vorty_tilt_trajc = missing_value

  allocate(vortz_tilt_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:vortz_tilt_trajc")
  vortz_tilt_trajc= missing_value

  allocate(vorts_tilt_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:vorts_tilt_trajc")
  vorts_tilt_trajc= missing_value

  allocate(vortc_tilt_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:vortc_tilt_trajc")
  vortc_tilt_trajc= missing_value

  allocate(vortx_gen_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:vortx_gen_trajc")
  vortx_gen_trajc = missing_value

  ALLOCATE(vorts_gen_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:vorts_gen_trajc")
  vorts_gen_trajc= missing_value

  ALLOCATE(vortc_gen_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:vortc_gen_trajc")
  vortc_gen_trajc= missing_value

  allocate(vortx_mix_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:vortx_mix_trajc")
  vortx_mix_trajc = missing_value

  allocate(vorty_mix_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:vorty_mix_trajc")
  vorty_mix_trajc = missing_value

  allocate(vortz_mix_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:vortz_mix_trajc")
  vortz_mix_trajc = missing_value

  allocate(vorts_mix_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:vorts_mix_trajc")
  vorts_mix_trajc = missing_value

  allocate(vortc_mix_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:vortc_mix_trajc")
  vortc_mix_trajc = missing_value

  allocate(vorts_exhg_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:vorts_exhg_trajc")
  vorts_exhg_trajc = missing_value

  allocate(vortc_exhg_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:vortc_exhg_trajc")
  vortc_exhg_trajc = missing_value

  ALLOCATE(buoy_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:buoy_trajc")
  buoy_trajc = missing_value

  ALLOCATE(buoyq_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:buoya_trajc")
  buoyq_trajc= missing_value

  ALLOCATE(frcs_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:frcs_trajc")
  frcs_trajc = missing_value

  ALLOCATE(pprt_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:pprt_trajc")
  pprt_trajc = missing_value

  ALLOCATE(ptprt_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:ptprt_trajc")
  ptprt_trajc = missing_value

  ALLOCATE(pt_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:pt_trajc")
  pt_trajc = missing_value

  ALLOCATE(qv_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:qv_trajc")
  qv_trajc = missing_value

  ALLOCATE(p_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:p_trajc")
  p_trajc = missing_value

  ALLOCATE(q_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:q_trajc")
  q_trajc = missing_value

  ALLOCATE(upgrad_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:upgrad_trajc")
  upgrad_trajc = missing_value

  ALLOCATE(vpgrad_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:vpgrad_trajc")
  vpgrad_trajc = missing_value

  ALLOCATE(wpgrad_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:wpgrad_trajc")
  wpgrad_trajc = missing_value

  ALLOCATE(utrajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:utrajc")
  utrajc = missing_value

  ALLOCATE(vtrajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:vtrajc")
  vtrajc = missing_value

  ALLOCATE(wtrajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:wtrajc")
  wtrajc = missing_value

  ALLOCATE(temp_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:temp_trajc")
  temp_trajc = 0.0

  ALLOCATE(ppi_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:ppi_trajc")
  ppi_trajc = missing_value


  ALLOCATE(qscalar_trajc(ntrajcs,ntimes,nscalar_max),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:qscalar_trajc")
  qscalar_trajc = 0.0

  IF ( misc_fields_flag > 0 ) THEN

    ALLOCATE(tke_trajc(ntrajcs,ntimes),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:tke_trajc")
    tke_trajc = 0.0

    ALLOCATE(kmh_trajc(ntrajcs,ntimes),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:kmh_trajc")
    kmh_trajc = 0.0

    ALLOCATE(kmv_trajc(ntrajcs,ntimes),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:kmv_trajc")
    kmv_trajc = 0.0

    ALLOCATE(ptw_trajc(ntrajcs,ntimes),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:ptw_trajc")
    ptw_trajc = 0.0

    ALLOCATE(mse_trajc(ntrajcs,ntimes),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:mse_trajc")
    mse_trajc = 0.0

  END IF

  ALLOCATE(pte_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:pte_trajc")
  pte_trajc = 0.0

  ALLOCATE(qscalar(nx,ny,nz,nscalar_max),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:qscalar")
  qscalar = 0.0

  ALLOCATE(rho_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:rho_trajc")
  rho_trajc = 0.0

  ALLOCATE(tem1_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:tem1_trajc")
  tem1_trajc = 0.0

  ALLOCATE(tem2_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:tem2_trajc")
  tem2_trajc = 0.0

  ALLOCATE(tem3_trajc(ntrajcs,ntimes),stat=istatus)
  CALL check_alloc_status(istatus, "arpstrajc:tem3_trajc")
  tem3_trajc = 0.0

  !ALLOCATE(xweight(ntrajcs,ntimes),stat=istatus)
  !CALL check_alloc_status(istatus, "arpstrajc:xweight")
  !
  !ALLOCATE(yweight(ntrajcs,ntimes),stat=istatus)
  !CALL check_alloc_status(istatus, "arpstrajc:yweight")
  !
  !ALLOCATE(zweight(ntrajcs,ntimes),stat=istatus)
  !CALL check_alloc_status(istatus, "arpstrajc:zweight")
  !
  !ALLOCATE(itrajc(ntrajcs,ntimes),stat=istatus)
  !CALL check_alloc_status(istatus, "arpstrajc:itrajc")
  !
  !ALLOCATE(jtrajc(ntrajcs,ntimes),stat=istatus)
  !CALL check_alloc_status(istatus, "arpstrajc:jtrajc")
  !
  !ALLOCATE(ktrajc(ntrajcs,ntimes),stat=istatus)
  !CALL check_alloc_status(istatus, "arpstrajc:ktrajc")

  IF(mprate_flag == 1) THEN

    ALLOCATE(mpcool_trajc(ntrajcs,ntimes,10),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:mpcool_trajc")
    mpcool_trajc = 0.0

    ALLOCATE(mpheat_trajc(ntrajcs,ntimes,17),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:mpheat_trajc")
    mpheat_trajc = 0.0

    ALLOCATE(mptrate_trajc(ntrajcs,ntimes),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:mptrate_trajc")
    mptrate_trajc = 0.0

    ALLOCATE(mpcoolproctot_trajc(ntrajcs,ntimes,10),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:mpcoolproctot_trajc")
    mpcoolproctot_trajc = 0.0

    ALLOCATE(mpheatproctot_trajc(ntrajcs,ntimes,17),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:mpheatproctot_trajc")
    mpheatproctot_trajc = 0.0

    ALLOCATE(mpcoolrate_trajc(ntrajcs,ntimes,10),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:mpcoolrate_trajc")
    mpcoolrate_trajc = 0.0

    ALLOCATE(mpheatrate_trajc(ntrajcs,ntimes,17),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:mpheatrate_trajc")
    mpheatrate_trajc = 0.0

    ALLOCATE(mpcooltot_trajc(ntrajcs,ntimes,10),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:mpcooltot_trajc")
    mpcooltot_trajc = 0.0

    ALLOCATE(mpheattot_trajc(ntrajcs,ntimes,17),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:mpheattot_trajc")
    mpheattot_trajc = 0.0

  END IF

  IF(therm_intg_flag == 1) THEN

    ALLOCATE(ptcalc_trajc(ntrajcs,ntimes),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:ptcalc_trajc")
    ptcalc_trajc = 0.0

    ALLOCATE(ptmix_trajc(ntrajcs,ntimes),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:ptmix_trajc")
    ptmix_trajc = 0.0

    ALLOCATE(ptcalcmp_trajc(ntrajcs,ntimes),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:ptcalcmp_trajc")
    ptcalcmp_trajc = 0.0

    ALLOCATE(ptcalcmix_trajc(ntrajcs,ntimes),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:ptcalcmix_trajc")
    ptcalcmix_trajc = 0.0

    ALLOCATE(qvcalc_trajc(ntrajcs,ntimes),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:qvcalc_trajc")
    qvcalc_trajc = 0.0

    ALLOCATE(qvmix_trajc(ntrajcs,ntimes),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:qvmix_trajc")
    qvmix_trajc = 0.0

    ALLOCATE(qvcalcmp_trajc(ntrajcs,ntimes),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:qvcalcmp_trajc")
    qvcalcmp_trajc = 0.0

    ALLOCATE(qvcalcmix_trajc(ntrajcs,ntimes),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:qvcalcmix_trajc")
    qvcalcmix_trajc = 0.0

    ALLOCATE(ptecalc_trajc(ntrajcs,ntimes),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:ptecalc_trajc")
    ptecalc_trajc = 0.0

    ALLOCATE(ptecalcmp_trajc(ntrajcs,ntimes),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:ptecalcmp_trajc")
    ptecalcmp_trajc = 0.0

    ALLOCATE(ptecalcmix_trajc(ntrajcs,ntimes),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:ptecalcmix_trajc")
    ptecalcmix_trajc = 0.0


  END IF

  IF(DSD_flag == 1) THEN

    ALLOCATE(N0_trajc(ntrajcs,ntimes,6),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:N0_trajc")
    N0_trajc = 0.d0

    ALLOCATE(alpha_trajc(ntrajcs,ntimes,6),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:alpha_trajc")
    alpha_trajc = 0.d0

    ALLOCATE(Dm_trajc(ntrajcs,ntimes,6),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:Dm_trajc")
    Dm_trajc = 0.0

    ALLOCATE(N0eff_trajc(ntrajcs,ntimes,6),stat=istatus)
    CALL check_alloc_status(istatus, "arpstrajc:N0eff_trajc")
    N0eff_trajc = 0.d0

  END IF

  nscalar  = 0
  nscalarq = 0
  IF (mphyopt == 1) THEN
    varid_qscalar(1) = 'qc'
    varid_qscalar(2) = 'qr'
    nscalar = 2
    nscalarq = 2

  ELSE IF(mphyopt >= 2 .and. mphyopt <= 7) THEN
    varid_qscalar(1) = 'qc'
    varid_qscalar(2) = 'qr'
    varid_qscalar(3) = 'qi'
    varid_qscalar(4) = 'qs'
    IF(mphyopt == 2) THEN
      varid_qscalar(5) = 'qh'
    ELSE
      varid_qscalar(5) = 'qg'
    END IF
    IF(DSD_flag == 1) THEN
      nscalar = 10
      nscalarq = 5
      varid_qscalar(6) = 'nc'
      varid_qscalar(7) = 'nr'
      varid_qscalar(8) = 'ni'
      varid_qscalar(9) = 'ns'
      IF(mphyopt == 2) THEN
        varid_qscalar(10) = 'nh'
      ELSE
        varid_qscalar(10) = 'ng'
      END IF
    ELSE
      nscalar = 5
      nscalarq = 5
    END IF
  ELSE IF(mphyopt == 8) THEN
    varid_qscalar(1) = 'qc'
    varid_qscalar(2) = 'qr'
    varid_qscalar(3) = 'qi'
    varid_qscalar(4) = 'qs'
    varid_qscalar(5) = 'qg'
    varid_qscalar(6) = 'qh'
    IF(DSD_flag == 1) THEN
      nscalar = 12
      nscalarq = 6
      varid_qscalar(7) = 'nc'
      varid_qscalar(8) = 'nr'
      varid_qscalar(9) = 'ni'
      varid_qscalar(10) = 'ns'
      varid_qscalar(11) = 'ng'
      varid_qscalar(12) = 'nh'
    ELSE
      nscalar = 6
      nscalarq = 6
    END IF
  ELSE IF(mphyopt == 9 .or. mphyopt == 10) THEN
    varid_qscalar(1) = 'qc'
    varid_qscalar(2) = 'qr'
    varid_qscalar(3) = 'qi'
    varid_qscalar(4) = 'qs'
    varid_qscalar(5) = 'qg'
    varid_qscalar(6) = 'qh'
    varid_qscalar(7) = 'nc'
    varid_qscalar(8) = 'nr'
    varid_qscalar(9) = 'ni'
    varid_qscalar(10) = 'ns'
    varid_qscalar(11) = 'ng'
    varid_qscalar(12) = 'nh'
    nscalar = 12
    nscalarq = 6
  ELSE IF(mphyopt == 11) THEN
    varid_qscalar(1) = 'qc'
    varid_qscalar(2) = 'qr'
    varid_qscalar(3) = 'qi'
    varid_qscalar(4) = 'qs'
    varid_qscalar(5) = 'qg'
    varid_qscalar(6) = 'qh'
    varid_qscalar(7) = 'nc'
    varid_qscalar(8) = 'nr'
    varid_qscalar(9) = 'ni'
    varid_qscalar(10) = 'ns'
    varid_qscalar(11) = 'ng'
    varid_qscalar(12) = 'nh'
    varid_qscalar(13) = 'zr'
    varid_qscalar(14) = 'zi'
    varid_qscalar(15) = 'zs'
    varid_qscalar(16) = 'zg'
    varid_qscalar(17) = 'zh'
    nscalar = 17
    nscalarq = 6
  ELSE
    WRITE(*,'(1x,I0,a)') 'ERROR: unsupported microphisics schems, mphyopt = ',mphyopt,'.'
    CALL arpsstop('ERROR: wrong value of mphyopt.',1)
  END IF

  ! Calculate cx (constant in Nt,N0,etc. calculations) for each species

  pi = 4.0 * atan(1.0)

  cx(1) = (pi/6.)*rhor
  cx(2) = (pi/6.)*rhor
  IF(mphyopt >= 8 .and. mphyopt <= 11) THEN
    cx(3) = 440.0      ! Constant value used in MY scheme for ice (bullet rosettes)
  ELSE
    cx(3) = (pi/6.)*rhoi
  END IF
  cx(4) = (pi/6.)*rhos
  IF(mphyopt > 2) THEN
    cx(5) = (pi/6.)*rhog
  ELSE
    cx(5) = (pi/6.)*rhoh
  END IF
  IF(mphyopt >=8) THEN
    cx(6) = (pi/6.)*rhoh
  END IF

!-----------------------------------------------------------------------
! Read x,y,z
!-----------------------------------------------------------------------

  CALL hdfreadxyz(nx,ny,nz,x,y,z, trim(grdbasfn), time)

  Print*,'x,y,z of input data read in.'
  print*,'x(1 )=',x(1)
  print*,'x(nx)=',x(nx)
  print*,'y(1 )=',y(1)
  print*,'y(ny)=',y(ny)

  dx = x(2) - x(1)
  dy = y(2) - y(1)

  IF( ireturn /= 0 ) THEN
    PRINT*,'Problem occured when trying to get dimensions from data.'
    PRINT*,'Program stopped.'
    STOP
  END IF

  WRITE(6,'(4(a,i5))') 'nx =',nx,', ny=',ny,', nz=',nz ,'nzsoil=',nzsoil

!-----------------------------------------------------------------------
!  First re-read data at starting (reference) time if necessary.
!-----------------------------------------------------------------------

  CALL hdfreadvar(nx,ny,nz,trim(grdbasfn),time,'zp', zp, itmp )

  CALL hdfreadvar(nx,ny,nz,trim(grdbasfn),time,'ptbar', ptbar, itmp )

  CALL hdfreadvar(nx,ny,nz,trim(grdbasfn),time,'pbar', pbar, itmp )

  CALL hdfreadvar(nx,ny,nz,trim(grdbasfn),time,'qvbar', qvbar, itmp )

  DO i=1,nx-1
    xs(i) = 0.5*(x(i)+x(i+1))
  END DO
  DO j=1,ny-1
    ys(j) = 0.5*(y(j)+y(j+1))
   END DO

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        rhobar(i,j,k) = pbar(i,j,k)/(rd*ptbar(i,j,k)*(pbar(i,j,k)/p0)**rddcp)
      END DO
    END DO
  END DO

  PRINT *,'nstart, nend, ninc =', nstart, nend, ninc

!-----------------------------------------------------------------------
!  Write out data along trajectories
!-----------------------------------------------------------------------

  DO k=1,ntimes

    idot = index( trim(trajc_fn_in(k)) , '.')

    CALL cvttsnd( ttrajc(nstart)  ,timsnd0, tmstrln0 )
    CALL cvttsnd( ttrajc(nend)    ,timsnd1, tmstrln1 )
    CALL cvttsnd( reftime(k),      timsnd2, tmstrln2 )

    fmtstr = trajc_fn_in(k)(1:idot-1)//'.trajc_'//                      &
        & trim(timsnd0)//'-'//trim(timsnd1)//'_'//trim(timsnd2)

    trajc_fn_out(1) = TRIM(fmtstr)//'.data'
    trajc_fn_out(2) = TRIM(fmtstr)//'.mprate1'
    trajc_fn_out(3) = TRIM(fmtstr)//'.mprate2'
    trajc_fn_out(4) = TRIM(fmtstr)//'.mprate3'
    trajc_fn_out(5) = TRIM(fmtstr)//'.mprate4'
    trajc_fn_out(6) = TRIM(fmtstr)//'.DSD'
    trajc_fn_out(7) = TRIM(fmtstr)//'.thermintg'
    trajc_fn_out(8) = TRIM(fmtstr)//'.trajcavg'

    DO i=1,8

      IF((i == 1) .or. (i >=2 .and. i <= 5 .and. mprate_flag == 1) .or.   &
         (i == 6 .and. DSD_flag == 1) .or. (i == 7 .and.                  &
          therm_intg_flag == 1) .or. (i == 8 .and. avg_trajc_flag == 1)) THEN

        CALL getunit(nunit(k,i))

        ltrajc_fn(i) = len_trim(trajc_fn_out(i))
        CALL fnversn( trajc_fn_out(i), ltrajc_fn(i) )

        OPEN(UNIT=nunit(k,i),FILE=trim(trajc_fn_out(i)),STATUS='unknown',   &
             FORM='formatted',IOSTAT= istat )

        counter = 0
        DO j = nstart, nend, ninc
          counter = counter + 1
        ENDDO

        WRITE(nunit(k,i),'(a)') trim(runname)
        WRITE(nunit(k,i),'(6e17.6)') xs(1),xs(nx-1),ys(1),ys(ny-1),zp(1,1,1),zp(1,1,nz-1)
        WRITE(nunit(k,i),'(3e17.6)') dx, dy, z(2)-z(1)

        WRITE(nunit(k,i),'(3e17.6)') ttrajc(nstart),ttrajc(nstart),ttrajc(nend)
        WRITE(nunit(k,i),'(i10)') counter
        WRITE(nunit(k,i),'(i10)') ntrajcs
      END IF

      IF(i == 1) WRITE(nunit(k,1),'(10a)')                              &
    't,      x,        y,        z,           ',                        &
    'u,        v,          w,        pt,            ',                  &
    'p,        qv,         ptprt,    pprt,          ',                  &
    'vortx,  vorty,  vortz,  vorts,  vortc,  ',                         &
    'buoy, buoyq,  frcs, upgrad,  vpgrad,  wpgrad, ',                   &
    'vortx_stch, vorty_stch, vorts_stch, vortc_stch, vortz_stch, ',     &
    'vortx_tilt, vorty_tilt, vorts_tilt, vortc_tilt, vortz_tilt, ',     &
    'vortx_mix,  vorty_mix,      vorts_mix,      vortc_mix,      vortz_mix,  ', &
    'vortx_gen,  vorty_gen,      vorts_gen,      vortc_gen,      ',     &
    'vorts_exhg, vortc_exhg'

    END DO

    IF(mprate_flag == 1) THEN

      WRITE(nunit(k,2),'(a,21(a,","),3a)')   &
      't,x,y,z,u,v,w,',                                                 &
      (varid_qscalar(i),i=1,nscalar),                                   &
      'evapqc,evapqr,sublqi,sublqs,sublqg',                             &
      'sublqh,meltqi,meltqs,meltqg,meltqh,condqc,condqr',               &
      'nuclqi,depoqi,depoqs,depoqg,depoqh,frzqci,colqci',               &
      'colqcs,colqcg,colqch,frzqrh,colqri,colqrs,colqrg,colqrh'

      WRITE(nunit(k,3),'(6a)')   &
      't,x,y,z,u,v,w,evapqctot,evapqrtot,sublqitot,sublqstot,',         &
      'sublqgtot,sublqhtot,meltqitot,meltqstot,meltqgtot,',             &
      'meltqhtot,condqctot,condqrtot,nuclqitot,depoqitot,',             &
      'depoqstot,depoqgtot,depoqhtot,frzqcitot,colqcitot,',             &
      'colqcstot,colqcgtot,colqchtot,frzqrhtot,colqritot,',             &
      'colqrstot,colqrgtot,colqrhtot'

      WRITE(nunit(k,4),'(6a)')   &
      't,x,y,z,u,v,w,evapqctrate,evapqrtrate,sublqitrate,sublqstrate,', &
      'sublqgtrate,sublqhtrate,meltqitrate,meltqstrate,meltqgtrate,',   &
      'meltqhtrate,condqctrate,condqrtrate,nuclqitrate,depoqitrate,',   &
      'depoqstrate,depoqgtrate,depoqhtrate,frzqcitrate,colqcitrate,',   &
      'colqcstrate,colqcgtrate,colqchtrate,frzqrhtrate,colqritrate,',   &
      'colqrstrate,colqrgtrate,colqrhtrate'

      WRITE(nunit(k,5),'(6a)')   &
      't,x,y,z,u,v,w,evapqcctot,evapqrctot,sublqictot,sublqsctot,',     &
      'sublqgctot,sublqhctot,meltqictot,meltqsctot,meltqgctot,',        &
      'meltqhctot,condqchtot,condqrhtot,nuclqihtot,depoqihtot,',        &
      'depoqshtot,depoqghtot,depoqhhtot,frzqcihtot,colqcihtot,',        &
      'colqcshtot,colqcghtot,colqchhtot,frzqrhhtot,colqrihtot,',        &
      'colqrshtot,colqrghtot,colqrhhtot'
    END IF

    IF(DSD_flag == 1) THEN
      WRITE(nunit(k,6),'(a,20(a,","),a,a)') 't,x,y,z,u,v,w,',           &
      (varid_qscalar(i),i=1,nscalar),                                   &
      'N0c,N0r,N0i,N0s,N0g,N0h,N0effc,N0effr',                          &
      'N0effi,N0effs,N0effg,N0effh,alphac,alphar,alphai,alphas',        &
      'alphag,alphah,Dmc,Dmr,Dmi,Dms,Dmg,Dmh'

    END IF

    IF(therm_intg_flag == 1) THEN

      WRITE(nunit(k,7),'(4a)')                                          &
      't,x,y,z,u,v,w,pt,ptcalcmix,ptcalcmp,ptcalc,qv,qvcalcmix,',       &
      'qvcalcmp,qvcalc,pte,ptecalcmix,ptecalcmp,ptecalc'

    END IF

    IF(avg_trajc_flag == 1) THEN
      WRITE(nunit(k,8),'(a,17(a,","),a,a)') 't,x,y,z,u,v,w,',           &
      (varid_qscalar(i),i=1,nscalar),                                   &
      'pt,ptcalcmix,ptcalcmp,ptcalc,qv,qvcalcmix,',                     &
      'qvcalcmp,qvcalc,pte,ptecalcmix,ptecalcmp,ptecalc'
    END IF

   END DO

!-----------------------------------------------------------------------
!  Start calculations along trajectories
!-----------------------------------------------------------------------

  DO npnt = nstart, nend, ninc

    DO k=1,ntimes
      DO j=1,ntrajcs
        xtrajc1(j,k)=xtrajc(npnt,j,k)
        ytrajc1(j,k)=ytrajc(npnt,j,k)
        ztrajc1(j,k)=ztrajc(npnt,j,k)
      END DO
    END DO

    ! DTD: only read history dump at history dump time interval

    IF(MOD((ttrajc(npnt)-tstart_calc),tinc_data) == 0) THEN

      time=ttrajc(npnt)

      ! print*,'Inside npnt loop, npnt = ', npnt

      ! May have ninc and hdmpfheader, hdmpftrailer depends on time

      CALL cvttsnd( ttrajc(npnt), timsnd, tmstrln )
      hisfile = trim(hdmpfheader)//'.hdf'//timsnd(1:tmstrln)//trim(hdmpftrailer)

      WRITE(6,'(/a,a,a)') ' Data set ', trim(hisfile),' to be read.'

    END IF

    CALL hdfreadvar(nx,ny,nz,trim(hisfile),time,'u', tem1, itmp )
    CALL hdfreadvar(nx,ny,nz,trim(hisfile),time,'v', tem2, itmp )
    CALL hdfreadvar(nx,ny,nz,trim(hisfile),time,'w', tem3, itmp )

    CALL a3dmax0(tem1,1,nx,1,nx,1,ny,1,ny-1,1,nz,1,nz-1, amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'umin    = ', amin,',  umax    =',amax
    CALL a3dmax0(tem2,1,nx,1,nx-1,1,ny,1,ny,1,nz,1,nz-1, amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'vmin    = ', amin,',  vmax    =',amax
    CALL a3dmax0(tem3,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz, amax,amin)
    WRITE(6,'(1x,2(a,e13.6))') 'wmin    = ', amin,',  wmax    =',amax

!-----------------------------------------------------------------------
!  Do diagnostics
!-----------------------------------------------------------------------

    IF( do_diagnostics == 1 .or. do_diagnostics == 10 ) THEN

      CALL cal_xvort_flat(tem2,tem3,y,zp,nx,ny,nz, tem4)
      CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem4,                        &
           xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, vortx_trajc )

      CALL cal_yvort_flat(tem1,tem3,x,zp,nx,ny,nz, tem5)
      CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem5,                        &
           xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, vorty_trajc )

      ! vertical vorticity
      CALL cal_zvort_flat(tem1,tem2,x,y ,nx,ny,nz, tem6)
      CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem6,                        &
           xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, vortz_trajc )

!-----------------------------------------------------------------------
! Stretching of vorticity
!-----------------------------------------------------------------------

      ! Stretching of vortx
      ! vortx_stch = vortx * du/dx
      !DO k=2,nz-1
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=2,nx-2
            tem7(i,j,k) = tem4(i,j,k) * (tem1(i+1,j,k)-tem1(i,j,k))/(dx)
          END DO
        END DO
      END DO
      !CALL edgfill(tem7,1,nx,2,nx-2, 1,ny,1,ny-1, 1,nz,2,nz-1)
      CALL edgfill(tem7,1,nx,2,nx-2, 1,ny,1,ny-1, 1,nz,1,nz-1)
      CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem7,                          &
           xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, vortx_stch_trajc )

      ! Stretching of vorty
      ! vorty_stch = vorty * dv/dy
      !DO k=2,nz-1
      DO k=1,nz-1
        DO j=2,ny-2
          DO i=1,nx-1
            tem7(i,j,k) = tem5(i,j,k) * (tem2(i,j+1,k)-tem2(i,j,k))/(dy)
          END DO
        END DO
      END DO
      !CALL edgfill(tem7,1,nx,1,nx-1, 1,ny,2,ny-2, 1,nz,2,nz-1)
      CALL edgfill(tem7,1,nx,1,nx-1, 1,ny,2,ny-2, 1,nz,1,nz-1)
      CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem7,                          &
           xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, vorty_stch_trajc )

      ! Stretching of vortz
      ! vortz_stch = vortz * dw/dz
      DO k=2,nz-2
        DO j=1,ny-1
          DO i=1,nx-1
            tem7(i,j,k) = tem6(i,j,k) * (tem3(i,j,k+1)-tem3(i,j,k))/(zp(i,j,k+1)-zp(i,j,k))
          END DO
        END DO
      END DO
      CALL edgfill(tem7,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,2,nz-2)
      CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem7,                          &
           xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, vortz_stch_trajc )

!-----------------------------------------------------------------------
! Tilting of vorticity
!-----------------------------------------------------------------------

      ! Tilting of vorticity into x-dir
      ! vortx_tilt = (vorty * du/dy) + (vortz * du/dz)
      DO k=2,nz-2
      !DO k=1,nz-2
        DO j=2,ny-2
          DO i=1,nx-1
            IF( k <= 2 ) THEN
              kl = k
              factor = 1.0
            ELSE
              kl = k-1
              factor = 0.5
            END IF

           !IF(k == 1) THEN
           !  kl = 2
           !ELSE
           !  kl = k
           !END IF

            tem7(i,j,k) = tem5(i,j,k) * 0.5*(tem1(i,j+1,k)+tem1(i+1,j+1,k)-tem1(i,j-1,k)-tem1(i+1,j-1,k))/(2*dy)  &
                        + tem6(i,j,k) * factor*0.5*(tem1(i,j,k+1)+tem1(i+1,j,k+1)-tem1(i,j,kl)-tem1(i+1,j,kl))/((zp(i,j,k+1)-zp(i,j,k)))
          END DO
        END DO
      END DO
      CALL edgfill(tem7,1,nx,1,nx-1, 1,ny,2,ny-2, 1,nz,2,nz-2)
      !CALL edgfill(tem7,1,nx,1,nx-1, 1,ny,2,ny-2, 1,nz,1,nz-2)
      CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem7,                          &
           xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, vortx_tilt_trajc )

      ! Tilting of vorticity into y-dir
      ! vorty_tilt = (vortx * dv/dx) + (vortz * dv/dz)
      DO k=2,nz-2
      !DO k=1,nz-2
        DO j=1,ny-1
          DO i=2,nx-2
            IF( k <= 2 ) THEN
              kl = k
              factor = 1.0
            ELSE
              kl = k-1
              factor = 0.5
            END IF

           !IF(k == 1) THEN
           !   kl = 2
           !ELSE
           !   kl = k
           !END IF

            tem7(i,j,k) = tem4(i,j,k) * 0.5*(tem2(i+1,j,k)+tem2(i+1,j+1,k)-tem2(i-1,j,k)-tem2(i-1,j+1,k))/(2*dx)  &
                         + tem6(i,j,k) * factor*0.5*(tem2(i,j,k+1)+tem2(i,j+1,k+1)-tem2(i,j,kl)-tem2(i,j+1,kl))/((zp(i,j,k+1)-zp(i,j,k)))
          END DO
        END DO
      END DO
      !CALL edgfill(tem7,1,nx,2,nx-2, 1,ny,1,ny-1, 1,nz,2,nz-2)
      CALL edgfill(tem7,1,nx,2,nx-2, 1,ny,1,ny-1, 1,nz,1,nz-2)
      CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem7,                        &
           xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, vorty_tilt_trajc )

      ! Tilting of vorticity into z-dir
      ! vortz_tilt = (vortx * dw/dx) + (vorty * dw/dy)
      DO k=1,nz-1
        DO j=2,ny-2
          DO i=2,nx-2
            tem7(i,j,k) = tem4(i,j,k) * 0.5*(tem3(i+1,j,k)+tem3(i+1,j,k+1)-tem3(i-1,j,k)-tem3(i-1,j,k+1))/(2*dx)  &
                         + tem5(i,j,k) * 0.5*(tem3(i,j+1,k)+tem3(i,j+1,k+1)-tem3(i,j-1,k)-tem3(i,j-1,k+1))/(2*dy)
          END DO
        END DO
      END DO
      CALL edgfill(tem7,1,nx,2,nx-2, 1,ny,2,ny-2, 1,nz,1,nz-1)
      CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem7,                        &
           xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, vortz_tilt_trajc )

!-----------------------------------------------------------------------
!
! Wind trajectories and averages
!
!-----------------------------------------------------------------------

      CALL avgx(tem1,0, nx,ny,nz,1,nx-1,1,ny-1,1,nz-1, tem4 )
      CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem4,                        &
           xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, utrajc )
      tem1 = tem4

      CALL avgy(tem2,0, nx,ny,nz,1,nx-1,1,ny-1,1,nz-1, tem4 )
      CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem4,                        &
           xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, vtrajc )
      tem2 = tem4

      CALL avgz(tem3,0, nx,ny,nz,1,nx-1,1,ny-1,1,nz-1, tem4 )
      CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem4,                        &
           xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, wtrajc )
      tem3 = tem4

      IF(avg_trajc_flag == 1) THEN
        DO k=1,ntimes
          DO i=1,ntrajcs
            IF( .not. trajc_missing(npnt,i,k)) THEN
              utrajcavg(npnt,k)=utrajcavg(npnt,k)+utrajc(i,k)/ntrajcsavg(npnt,k)
              vtrajcavg(npnt,k)=vtrajcavg(npnt,k)+vtrajc(i,k)/ntrajcsavg(npnt,k)
              wtrajcavg(npnt,k)=wtrajcavg(npnt,k)+wtrajc(i,k)/ntrajcsavg(npnt,k)
            END IF
          END DO
        END DO
      END IF

!-----------------------------------------------------------------------
!
! Angle/direction of horizontal wind: psi = atan(v/u)
! tem5 = psi
!
!-----------------------------------------------------------------------

      PI = 4.0 * atan(1.0)

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            IF (tem1(i,j,k).eq.0.0 .and. tem2(i,j,k).gt.0.0) THEN
              tem5 (i,j,k) = (pi/2.0)
            ELSE IF (tem1(i,j,k).eq.0.0 .and. tem2(i,j,k).lt.0.0) THEN
              tem5 (i,j,k) = 3.0*pi/2.0
            ELSE IF (tem1(i,j,k).gt.0.0 .and. tem2(i,j,k).eq.0.0) THEN
              tem5 (i,j,k) = 0.0
            ELSE IF (tem1(i,j,k).lt.0.0 .and. tem2(i,j,k).eq.0.0) THEN
              tem5 (i,j,k) = pi
            ELSE IF (tem1(i,j,k).gt.0.0 .and. tem2(i,j,k).gt.0.0) THEN
              tem5 (i,j,k) = atan(tem2(i,j,k)/tem1(i,j,k))
            ELSE IF (tem1(i,j,k).gt.0.0 .and. tem2(i,j,k).lt.0.0) THEN
              tem5 (i,j,k) = atan(tem2(i,j,k)/tem1(i,j,k)) + (2.0*pi)
            ELSE IF (tem1(i,j,k).lt.0.0 .and. tem2(i,j,k).lt.0.0) THEN
              tem5 (i,j,k) = atan(tem2(i,j,k)/tem1(i,j,k)) + (pi)
            ELSE IF (tem1(i,j,k).lt.0.0 .and. tem2(i,j,k).gt.0.0) THEN
              tem5 (i,j,k) = atan(tem2(i,j,k)/tem1(i,j,k)) + (pi)
            ELSE
              tem5 (i,j,k) = 0.0
            ENDIF
          END DO
        END DO
      END DO
      !
      ! horizontal wind speed
      !
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            tem6(i,j,k) = sqrt( (tem1(i,j,k))**2 + (tem2(i,j,k))**2 )
          END DO
        END DO
      END DO

      CALL intrp_trajc(nx,ny,nz, xs,ys,zp, tem5,                       &
           xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, psi_trajc )
      CALL intrp_trajc(nx,ny,nz, xs,ys,zp, tem6,                       &
           xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, vh_trajc )

!-----------------------------------------------------------------------
! Horizontal streamwise and crosswise vorticity
!-----------------------------------------------------------------------

      DO k=1,ntimes
        DO ntrj=1,ntrajcs
          utem = utrajc(ntrj,k)
          vtem = vtrajc(ntrj,k)
          wtem = wtrajc(ntrj,k)

          ! Streamwise vorticity
          ! vorts = (vortx*u + vorty*v) / V_h
          vorts_trajc(ntrj,k)=( vortx_trajc(ntrj,k)*utem                &
                               +vorty_trajc(ntrj,k)*vtem ) /            &
                               (1.0e-10 + sqrt(utem**2+vtem**2))

          ! Crosswise vorticity
          ! vortc = (-vortx*v + vorty*u) / V_h
          vortc_trajc(ntrj,k)=(-vortx_trajc(ntrj,k)*vtem                &
                               +vorty_trajc(ntrj,k)*utem) /             &
                               (1.0e-10 + sqrt(utem**2+vtem**2))

          ! Stretching of streamwise vorticity
          ! vorts_stch = (vortx_stch*u + vorty_stch*v) / V_h
          vorts_stch_trajc(ntrj,k)=( vortx_stch_trajc(ntrj,k)*utem     &
                               +vorty_stch_trajc(ntrj,k)*vtem)         &
              /(1.0e-10 + sqrt(utem*utem + vtem*vtem))


          ! Stretching of crosswise vorticity
          ! vortc_stch = (-vortx_stch*v + vorty_stch*u) / V_h
          vortc_stch_trajc(ntrj,k)=(-vortx_stch_trajc(ntrj,k)*vtem     &
                               +vorty_stch_trajc(ntrj,k)*utem)         &
              /(1.0e-10 + sqrt(utem*utem + vtem*vtem))

          ! Tilting of vorticity into streamwise direction
          ! vorts_tilt = (vortx_tilt*u + vorty_tilt*v) / V_h
          vorts_tilt_trajc(ntrj,k)=( vortx_tilt_trajc(ntrj,k)*utem     &
                               +vorty_tilt_trajc(ntrj,k)*vtem)         &
              /(1.0e-10 + sqrt(utem*utem + vtem*vtem))

          ! Tilting of vorticity into crosswise direction
          ! vortc_tilt = (-vortx_tilt*v + vorty_tilt*u) / V_h
          vortc_tilt_trajc(ntrj,k)=(-vortx_tilt_trajc(ntrj,k)*vtem     &
                                +vorty_tilt_trajc(ntrj,k)*utem)        &
               /(1.0e-10 + sqrt(utem*utem + vtem*vtem))

        END DO
      END DO

    ENDIF

    IF( do_diagnostics == 2 .or. do_diagnostics == 10 ) then

!-----------------------------------------------------------------------
! Buoyancy term
!-----------------------------------------------------------------------

       CALL hdfreadvar(nx,ny,nz,trim(hisfile),time,'qc', tem1, itmp )
       CALL hdfreadvar(nx,ny,nz,trim(hisfile),time,'qr', tem2, itmp )

       CALL edgfill(tem1,1,nx,2,nx-2, 1,ny,2,ny-2, 1,nz,2,nz-2)
       CALL edgfill(tem2,1,nx,2,nx-2, 1,ny,2,ny-2, 1,nz,2,nz-2)
       tem4 = tem1 + tem2

       !DTD: attempt to read in other ice variables and add to the loading
       tem1 = 0.0
       tem2 = 0.0
       tem3 = 0.0
       CALL hdfreadvar(nx,ny,nz,trim(hisfile),time,'qi', tem1, itmp )
       CALL hdfreadvar(nx,ny,nz,trim(hisfile),time,'qs', tem2, itmp )
       IF(mphyopt >= 8 .and. mphyopt <= 11) THEN
         CALL hdfreadvar(nx,ny,nz,trim(hisfile),time,'qg', tem3, itmp )
       ELSE
         tem3 = 0.0
       END IF

       tem4 = tem4 + tem1 + tem2 + tem3

       tem1 = 0.0
       CALL hdfreadvar(nx,ny,nz,trim(hisfile),time,'qh', tem1, itmp )
       tem4 = tem4 + tem1

       DO k=1,nz-1
         DO j=1,ny-1
           DO i=1,nx-1
             tem4(i,j,k) = -g*tem4(i,j,k)/(1+qvbar(i,j,k))
           END DO
         END DO
       END DO

       CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem4,                        &
            xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, buoyq_trajc )

       CALL hdfreadvar(nx,ny,nz,trim(hisfile),time,'pt', tem1, itmp )
       CALL hdfreadvar(nx,ny,nz,trim(hisfile),time,'p',  tem2, itmp )
       CALL hdfreadvar(nx,ny,nz,trim(hisfile),time,'qv', tem3, itmp )

       CALL edgfill(tem1,1,nx,2,nx-2, 1,ny,2,ny-2, 1,nz,2,nz-2)
       CALL edgfill(tem2,1,nx,2,nx-2, 1,ny,2,ny-2, 1,nz,2,nz-2)
       CALL edgfill(tem3,1,nx,2,nx-2, 1,ny,2,ny-2, 1,nz,2,nz-2)

       CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem1,                        &
            xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, pt_trajc )
       CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem2,                        &
            xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, p_trajc )
       CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem3,                        &
            xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, qv_trajc )

       DO k=1,nz
         DO j=1,ny
           DO i=1,nx
             tem1(i,j,k) = tem1(i,j,k)-ptbar(i,j,k)
             tem2(i,j,k) = tem2(i,j,k)- pbar(i,j,k)
             tem3(i,j,k) = tem3(i,j,k)-qvbar(i,j,k)
           ENDDO
         ENDDO
       ENDDO

       CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem1,                        &
            xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, ptprt_trajc )
       CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem2,                        &
            xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, pprt_trajc )
!
!-----------------------------------------------------------------------
!  wbuoy = g ( ptprt/ptbar-pprt/(cpdcv*pbar)+        ! Dan Dawson's fix
!          qvprt/(rddrv+qvbar)-(qvprt+qc+qr+qs+qi+qh)/(1+qvbar)
!          -(ptprt*ptprt)/(ptbar*ptbar)                        !2nd-order
!          +0.5*(ptprt*pprt)/(cpdcv*ptbar*pbar))               !2nd-order
!-----------------------------------------------------------------------

       DO k=1,nz-1
         DO j=1,ny-1
           DO i=1,nx-1
             pttem = tem1(i,j,k)/ptbar(i,j,k)
             tem1(i,j,k) = g*( pttem - tem2(i,j,k)/((cpdcv)*pbar(i,j,k))&
                -pttem*pttem+0.5*pttem*tem2(i,j,k)/(cpdcv*pbar(i,j,k))  &
                + tem3(i,j,k)/(rddrv+qvbar(i,j,k)))                     &
              + tem4(i,j,k) - g*tem3(i,j,k)/(1.+qvbar(i,j,k))   ! Dan Dawson's fix
           END DO
         END DO
       END DO

       CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem1,                        &
            xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, buoy_trajc )

!-----------------------------------------------------------------------
! Baroclinic vorticity generation
!-----------------------------------------------------------------------

      ! tem2 = theta,  tem3 = p
      CALL hdfreadvar(nx,ny,nz,trim(hisfile),time,'pt', tem2, itmp )
      CALL hdfreadvar(nx,ny,nz,trim(hisfile),time,'p',  tem3, itmp )

      CALL edgfill(tem2,1,nx,2,nx-2, 1,ny,2,ny-2, 1,nz,2,nz-2)
      CALL edgfill(tem3,1,nx,2,nx-2, 1,ny,2,ny-2, 1,nz,2,nz-2)

      ! tem5 = qv
      CALL hdfreadvar(nx,ny,nz,trim(hisfile),time,'qv', tem5, itmp )

      CALL edgfill(tem5,1,nx,2,nx-2, 1,ny,2,ny-2, 1,nz,2,nz-2)

      tem6(:,:,:) = 0.0  ! the total hydrometeor/condensate mixing ratio
      DO nq = 1, nscalarq

        CALL hdfreadvar(nx,ny,nz,trim(hisfile),time,varid_qscalar(nq), tem7, itmp )
        CALL edgfill(tem7,1,nx,2,nx-2, 1,ny,2,ny-2, 1,nz,2,nz-2)

        tem6(:,:,:) = tem6(:,:,:) + tem7(:,:,:)
      END DO

      ! tem2 = theta_rho (density potential temperature)
      ! Incorporates effects of virtual temp. correction (qv) and condensate loading (qc, qr)
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            ! theta_rho = theta * (1 + 0.61*qv - qht) , where qht is the total hydrometeor/condensate mixing ratio
            tem2(i,j,k) = tem2(i,j,k) * (1 + 0.61*tem5(i,j,k) - tem6(i,j,k))
          END DO
        END DO
      END DO

      ! tem4 = rho
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            ! rho = p / (R_d * theta_rho * (p/p0)^(R_d/c_p))
            tem4(i,j,k) = tem3(i,j,k) / (rd * tem2(i,j,k) * (tem3(i,j,k)/p0)**rddcp)
          END DO
        END DO
      END DO

      ! Baroclinic generation of vortx
      ! vortx_gen = (1 / rho^2) * (-(d(rho)/dz * dp/dy) + (d(rho)/dy * dp/dz))
      DO k=2,nz-2
        DO j=2,ny-2
          DO i=2,nx-2
            tem5(i,j,k) = (1 / (tem4(i,j,k)**2)) * (                   &
               -(tem4(i,j,k+1)-tem4(i,j,k-1))/(zp(i,j,k+1)-zp(i,j,k-1))&
                * (tem3(i,j+1,k)-tem3(i,j-1,k))/(2*dy)                 &
               +(tem3(i,j,k+1)-tem3(i,j,k-1))/(zp(i,j,k+1)-zp(i,j,k-1))&
                * (tem4(i,j+1,k)-tem4(i,j-1,k))/(2*dy)                 &
               )
            !OLD BUOYANCY FORMULATION:
            !tem5(i,j,k) = (tem1(i,j+1,k) - tem1(i,j-1,k)) / (2*dy)
          END DO
        END DO
      END DO
      CALL edgfill(tem5,1,nx,2,nx-2, 1,ny,2,ny-2, 1,nz,2,nz-2)
      CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem5,                        &
           xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, vortx_gen_trajc )

      ! Baroclinic generation of vorty
      ! vorty_gen = (1 / rho^2) * ((d(rho)/dz * dp/dx) - (d(rho)/dx * dp/dz))
      DO k=2,nz-2
        DO j=2,ny-2
          DO i=2,nx-2
             tem5(i,j,k) = (1 / (tem4(i,j,k)*tem4(i,j,k))) * (         &
                (tem4(i,j,k+1)-tem4(i,j,k-1))/(zp(i,j,k+1)-zp(i,j,k-1))&
                * (tem3(i+1,j,k)-tem3(i-1,j,k))/(2*dx)                 &
               -(tem3(i,j,k+1)-tem3(i,j,k-1))/(zp(i,j,k+1)-zp(i,j,k-1))&
                * (tem4(i+1,j,k)-tem4(i-1,j,k))/(2*dx)                 &
               )
            !OLD BUOYANCY FORMULATION:
            !tem5(i,j,k) = -(tem1(i+1,j,k) - tem1(i-1,j,k)) / (2*dx)
          END DO
        END DO
      END DO
      CALL edgfill(tem5,1,nx,2,nx-2, 1,ny,2,ny-2, 1,nz,2,nz-2)
      CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem5,                        &
           xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, vorty_gen_trajc )

      DO k=1,ntimes
        DO ntrj=1,ntrajcs
          utem = utrajc(ntrj,k)
          vtem = vtrajc(ntrj,k)

          ! Baroclinic generation of streamwise vorticity
          ! vorts_gen = (vortx_gen*u + vorty_gen*v) / V_h
          vorts_gen_trajc(ntrj,k)=( vortx_gen_trajc(ntrj,k)*utem       &
                                   +vorty_gen_trajc(ntrj,k)*vtem)      &
            /(1.0e-10 + sqrt(utem*utem + vtem*vtem))

          ! Baroclinic generation of crosswise vorticity
          ! vortc_gen = (-vortx_gen*v + vorty_gen*u) / V_h
          vortc_gen_trajc(ntrj,k)=( -vortx_gen_trajc(ntrj,k)*vtem      &
                                    +vorty_gen_trajc(ntrj,k)*utem)     &
            /(1.0e-10 + sqrt(utem*utem + vtem*vtem))
        END DO
      END DO

!-----------------------------------------------------------------------
! Mixing term (of vorticity equation)
!-----------------------------------------------------------------------

       IF (have_mixing_terms > 0) THEN

         dirstrpos  = SCAN(hdmpfheader,'/',BACK=.TRUE.)
         mixdirname = hdmpfheader(:dirstrpos-1)
         mixrunname = hdmpfheader(dirstrpos+1:)

         varunits(1)= '(m/s**2)'
         varunits(2)= '(m/s**2)'
         varunits(3)= '(m/s**2)'
         varname(1) = 'Total mixing of u'
         varname(2) = 'Total mixing of v'
         varname(3) = 'Total mixing of w'

         ! Set tem1 = umix, tem2 = vmix, tem3 = wmix
         varid = 'mix_u'
         CALL readvar2(nx,ny,nz,tem1,varid,varname(1),varunits(1),time, &
                       mixrunname,mixdirname,3,0,istatus)
         varid = 'mix_v'
         CALL readvar2(nx,ny,nz,tem2,varid,varname(2),varunits(2),time, &
                       mixrunname,mixdirname,3,0,istatus)
         varid = 'mix_w'
         CALL readvar2(nx,ny,nz,tem3,varid,varname(3),varunits(3),time, &
                       mixrunname,mixdirname,3,0,istatus)

         ! tem4 = vortx_mix
         CALL cal_xvortmix_flat(tem2,tem3,y,zp,nx,ny,nz, tem4)
         CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem4,                      &
              xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, vortx_mix_trajc )

         ! tem5 = vorty_mix
         CALL cal_yvortmix_flat(tem1,tem3,x,zp,nx,ny,nz, tem5)
         CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem5,                      &
              xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, vorty_mix_trajc )

         ! tem6 = vortz_mix
         CALL cal_zvortmix_flat(tem1,tem2,x,y ,nx,ny,nz, tem6)
         CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem6,                      &
              xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, vortz_mix_trajc )

         DO k=1,ntimes
           DO ntrj=1,ntrajcs
             utem = utrajc(ntrj,k)
             vtem = vtrajc(ntrj,k)

             ! Streamwise mixing generation
             ! vorts_mix = (vortx_mix*u + vorty_mix*v) / V_h
             vorts_mix_trajc(ntrj,k)=( vortx_mix_trajc(ntrj,k)*utem     &
                                      +vorty_mix_trajc(ntrj,k)*vtem)    &
               /(1.0e-10 + sqrt(utem*utem + vtem*vtem))

             ! Crosswise mixing generation
             ! vortc_mix = (-vortx_mix*v + vorty_mix*u) / V_h
             vortc_mix_trajc(ntrj,k)=( -vortx_mix_trajc(ntrj,k)*vtem    &
                                       +vorty_mix_trajc(ntrj,k)*utem)   &
               /(1.0e-10 + sqrt(utem*utem + vtem*vtem))
           END DO
         END DO

       END IF

!-----------------------------------------------------------------------
! Exchange terms between streamwise and crosswise vorticity
!
! NOTE: These terms require finite differencing in time, unlike the
!       other terms. Therefore, they depend upon the trajectory position
!       at the previous and/or next time step.
!-----------------------------------------------------------------------

       DO k=1,ntrajcs
         DO j=1,ntimes

           ! Correction applied if this is the first or last timestep
           ! in the arpstrajc output file. In these cases, we'll simply
           ! assume dpsi/dt is the same value as at the next/previous
           ! timestep, since we don't have the information to do a proper
           ! finite-difference in time.
           IF (npnt == 1) THEN
              i = npnt + 1
           ELSE IF (npnt == (tend-tstart+1)) THEN
              i = npnt - 1
           ELSE
              i = npnt
           END IF

           ! Calculate average psi between previous and current timestep
           !  i.e., centered at t = n - (1/2)
           del_x = xtrajc(i,k,j) - xtrajc(i-1,k,j)
           del_y = ytrajc(i,k,j) - ytrajc(i-1,k,j)

           IF (del_x.eq.0.0 .and. del_y.gt.0.0) THEN
              psi_prev = (pi/2.0)
           ELSE IF (del_x.eq.0.0 .and. del_y.lt.0.0) THEN
              psi_prev = 3.0*pi/2.0
           ELSE IF (del_x.gt.0.0 .and. del_y.eq.0.0) THEN
              psi_prev = 0.0
           ELSE IF (del_x.lt.0.0 .and. del_y.eq.0.0) THEN
              psi_prev = pi
           ELSE IF (del_x.gt.0.0 .and. del_y.gt.0.0) THEN
              psi_prev = atan(del_y/del_x)
           ELSE IF (del_x.gt.0.0 .and. del_y.lt.0.0) THEN
              psi_prev = atan(del_y/del_x) + (2.0*pi)
           ELSE IF (del_x.lt.0.0 .and. del_y.lt.0.0) THEN
              psi_prev = atan(del_y/del_x) + (pi)
           ELSE IF (del_x.lt.0.0 .and. del_y.gt.0.0) THEN
              psi_prev = atan(del_y/del_x) + (pi)
           ELSE
              psi_prev = 0.0
           END IF

           ! Calculate average psi between current and next timestep
           !  i.e., centered at t = n + (1/2)
           del_x = xtrajc(i+1,k,j) - xtrajc(i,k,j)
           del_y = ytrajc(i+1,k,j) - ytrajc(i,k,j)

           IF (del_x.eq.0.0 .and. del_y.gt.0.0) THEN
              psi_next = (pi/2.0)
           ELSE IF (del_x.eq.0.0 .and. del_y.lt.0.0) THEN
              psi_next = 3.0*pi/2.0
           ELSE IF (del_x.gt.0.0 .and. del_y.eq.0.0) THEN
              psi_next = 0.0
           ELSE IF (del_x.lt.0.0 .and. del_y.eq.0.0) THEN
              psi_next = pi
           ELSE IF (del_x.gt.0.0 .and. del_y.gt.0.0) THEN
              psi_next = atan(del_y/del_x)
           ELSE IF (del_x.gt.0.0 .and. del_y.lt.0.0) THEN
              psi_next = atan(del_y/del_x) + (2.0*pi)
           ELSE IF (del_x.lt.0.0 .and. del_y.lt.0.0) THEN
              psi_next = atan(del_y/del_x) + (pi)
           ELSE IF (del_x.lt.0.0 .and. del_y.gt.0.0) THEN
              psi_next = atan(del_y/del_x) + (pi)
           ELSE
              psi_next = 0.0
           END IF

           ! Temporal change in psi at current timestep (centered difference)
           del_psi = psi_next - psi_prev

           ! If abs(del_psi) > pi, correct the value by taking periodicity into consideration
           IF (del_psi.gt.pi) THEN
              del_psi = (2 * pi) - del_psi
           ELSE IF (del_psi.lt.(-pi)) THEN
              del_psi = del_psi + (2 * pi)
           END IF

           ! vorts_exhg = vortc * d(psi)/dt
           vorts_exhg_trajc(k,j) = vortc_trajc(k,j) * del_psi / tinc_calc
           ! vortc_exhg = -vorts * d(psi)/dt
           vortc_exhg_trajc(k,j) = -(vorts_trajc(k,j) * del_psi) / tinc_calc
         END DO
       END DO

    END IF  ! do_diagnostics == 2

    IF( do_diagnostics == 3 .or. do_diagnostics == 10 ) then

!-----------------------------------------------------------------------
! Calculate pressure gradient force components
!-----------------------------------------------------------------------

       CALL hdfreadvar(nx,ny,nz,trim(hisfile),time,'p',  tem2, itmp )

       CALL edgfill(tem2,1,nx,2,nx-2, 1,ny,2,ny-2, 1,nz,2,nz-2)
       tem2(:,:,:) = tem2(:,:,:) - pbar(:,:,:)

       DO k=2,nz-2

         IF( k == 2 ) THEN
           kl = 2
           factor = 1.0
         ELSE
           kl = k-1
           factor = 2.0
         ENDIF

         DO j=1,ny-1
           DO i=1,nx-1
             !rhobar(i,j,k) = pbar(i,j,k)/(rd*ptbar(i,j,k)*(pbar(i,j,k)/p0)**rddcp)
             tem1(i,j,k)= -(tem2(i,j,k+1)-tem2(i,j,kl))/                &
                          (factor*rhobar(i,j,k)*(zp(i,j,k+1)-zp(i,j,k)))
           END DO
         END DO

!-----------------------------------------------------------------------
! Note that we are not including the correction term due to non-flat
! coordinate surfaces at this time.
!-----------------------------------------------------------------------

         DO j=1,ny-1
           DO i=2,nx-2
             tem3(i,j,k)=-(tem2(i+1,j,k)-tem2(i-1,j,k))/(2*dx*rhobar(i,j,k))
           END DO
         END DO

         DO j=2,ny-2
           DO i=1,nx-1
             tem4(i,j,k)=-(tem2(i,j+1,k)-tem2(i,j,k))/(2*dy*rhobar(i,j,k))
           END DO
         END DO
       END DO

       CALL edgfill(tem1,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,2,nz-2)
       CALL edgfill(tem3,1,nx,2,nx-2, 1,ny,1,ny-1, 1,nz,2,nz-2)
       CALL edgfill(tem4,1,nx,1,nx-1, 1,ny,2,ny-2, 1,nz,2,nz-2)

       CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem1,                        &
            xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, wpgrad_trajc )
       CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem3,                        &
            xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, upgrad_trajc )
       CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem4,                        &
            xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, vpgrad_trajc )

       DO k=1,ntimes
         DO ntrj=1,ntrajcs
           utem = utrajc(ntrj,k)
           vtem = vtrajc(ntrj,k)
           wtem = wtrajc(ntrj,k)
           frcs_trajc(ntrj,k)=(upgrad_trajc(ntrj,k)*utem                &
                              +vpgrad_trajc(ntrj,k)*vtem+               &
               (wpgrad_trajc(ntrj,k)+buoy_trajc(ntrj,k))*wtem)          &
               /(1.0e-10+sqrt(utem**2+vtem**2+wtem**2))
         ENDDO
       ENDDO

    END IF

!-----------------------------------------------------------------------
!
! Misc fields (tke, kmh, kmv etc.)
!
!-----------------------------------------------------------------------

     IF(misc_fields_flag == 1) THEN
       !DTD: turbulent kinetic energy

       CALL hdfreadvar(nx,ny,nz,trim(hisfile),time,'tke', tem1, itmp )
       CALL edgfill(tem1,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
       CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem1,                   &
            xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, tke_trajc )

       !DTD: horizontal mixing coefficient (kmh)

       CALL hdfreadvar(nx,ny,nz,trim(hisfile),time,'kmh', tem1, itmp )
       CALL edgfill(tem1,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
       CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem1,                   &
            xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, kmh_trajc )

       !DTD: vertical mixing coefficient (kmv)

       CALL hdfreadvar(nx,ny,nz,trim(hisfile),time,'kmv', tem1, itmp )
       CALL edgfill(tem1,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
       CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem1,                   &
            xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, kmv_trajc )

       !DTD: equivalent potential temperature

       CALL hdfreadvar(nx,ny,nz,trim(hisfile),time,'p', tem1, itmp )
       CALL hdfreadvar(nx,ny,nz,trim(hisfile),time,'pt', tem2, itmp )
       CALL hdfreadvar(nx,ny,nz,trim(hisfile),time,'qv', tem3, itmp )

       CALL pt2pte(nx,ny,nz,1,nx-1,1,ny-1,1,nz-1,tem1,tem2,tem3,tem4)

       CALL edgfill(tem4,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
       CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem4,                   &
            xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, pte_trajc )

       IF(avg_trajc_flag == 1) THEN
         DO k=1,ntimes
           DO i=1,ntrajcs
             IF( .not. trajc_missing(npnt,i,k)) THEN
               pte_trajcavg(npnt,k)=pte_trajcavg(npnt,k)+pte_trajc(i,k)/ntrajcsavg(npnt,k)
             END IF
           END DO
         END DO
       END IF

       !DTD: wet-bulb potential temperature

       CALL edgfill(tem3,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
       CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem3,                   &
            xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, tem1_trajc )

       DO k=1,ntimes
         DO i=1,ntrajcs
           ptw_trajc(i,k) = pte_trajc(i,k) - (LV/CP)*tem1_trajc(i,k)
         END DO
       END DO

       !DTD: Moist static energy

       DO k=1,nz-1
         DO j=1,ny-1
           DO i=1,nx-1
             ! Temperature (in tem2)
             tem2(i,j,k) = tem2(i,j,k)*(tem1(i,j,k)/p0)**rddcp

             ! Moist static energy (in tem1)
             tem1(i,j,k) = cp*tem2(i,j,k) + 0.5*g*(zp(i,j,k)+zp(i,j,k+1)) + LV*tem3(i,j,k)

           END DO
         END DO
       END DO

       CALL edgfill(tem1,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
       CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem1,                   &
            xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, mse_trajc )

     END IF

     IF(DSD_flag == 1 .or. mprate_flag == 1) THEN
       !DTD: Calculate various microphysical source/sink terms along trajectories

       qscalar = 0.0

       !mphyopt = 0
       !print*,'nx,ny,nz',nx,ny,nz

       ! Read in needed variables to calculate air density (rho) which is used in several of the
       ! following calculations.

       CALL hdfreadvar(nx,ny,nz,trim(hisfile),time,'p', tem1, itmp )
       CALL hdfreadvar(nx,ny,nz,trim(hisfile),time,'pt', tem2, itmp )

       CALL edgfill(tem1,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
       CALL edgfill(tem2,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)

       DO k=2,nz-2
         DO j=2,ny-2
           DO i=2,nx-2
             tem3(i,j,k) = tem1(i,j,k)/(rd*tem2(i,j,k)*(tem1(i,j,k)/p0)**rddcp)  ! rho is in tem3
           END DO
         END DO
       END DO

       CALL edgfill(tem3,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
       ! Interpolate rho to trajectory point
       CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem3,                   &
          xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, rho_trajc)

       ! Calculate temperature along trajectory

       DO k=2,nz-2
         DO j=2,ny-2
           DO i=2,nx-2
             tem3(i,j,k) = tem2(i,j,k)*(tem1(i,j,k)/p0)**rddcp   ! temperature is in tem3
           END DO
         END DO
       END DO

       CALL edgfill(tem3,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
       ! Interpolate temperature to trajectory point
       CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem3,                        &
          xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, temp_trajc)

       ! Interpolate microphysics scalars to trajectory
       IF(mphyopt < 9) THEN
         DO nq=1,nscalarq
           CALL hdfreadvar(nx,ny,nz,trim(hisfile),time,varid_qscalar(nq), tem1, itmp )
           CALL edgfill(tem1,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
           CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem1,                    &
               xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, qscalar_trajc(:,:,nq) )
           IF(avg_trajc_flag == 1) THEN
             DO k=1,ntimes
               DO i=1,ntrajcs
                 IF( .not. trajc_missing(npnt,i,k)) THEN
                   qscalar_trajcavg(npnt,k,nq)=qscalar_trajcavg(npnt,k,nq) &
                               +qscalar_trajc(i,k,nq)/ntrajcsavg(npnt,k)
                 END IF
               END DO
             END DO
           END IF
         END DO
       ELSE
         DO nq=1,nscalar
           CALL hdfreadvar(nx,ny,nz,trim(hisfile),time,varid_qscalar(nq), tem1, itmp )
           CALL edgfill(tem1,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
           CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem1,                    &
               xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, qscalar_trajc(:,:,nq) )
           IF(avg_trajc_flag == 1) THEN
             DO k=1,ntimes
               DO i=1,ntrajcs
                 IF( .not. trajc_missing(npnt,i,k)) THEN
                   qscalar_trajcavg(npnt,k,nq)=qscalar_trajcavg(npnt,k,nq) &
                               +qscalar_trajc(i,k,nq)/ntrajcsavg(npnt,k)
                 END IF
               END DO
             END DO
           END IF
         END DO
       END IF

     END IF

     IF(DSD_flag == 1) THEN

       IF(mphyopt == 2) THEN
         DO nq=1,nscalarq
           ! Calculate N0 along trajectories (in this case constant everywhere)

           IF(nq == 1) THEN
             N0_trajc(:,:,nq) = 0.d0   ! Not applicable to monodisperse cloud
           ELSE IF(nq == 2) THEN
             N0_trajc(:,:,nq) = dble(N0rfix)
           ELSE IF(nq == 3) THEN
             N0_trajc(:,:,nq) = 0.d0   ! Not applicable to monodisperse ice
           ELSE IF(nq == 4) THEN
             N0_trajc(:,:,nq) = dble(N0sfix)
           ELSE IF(nq == 5) THEN
             N0_trajc(:,:,nq) = dble(N0hfix)
           END IF

           ! Calculate alpha along trajectories (in this case zero everywhere)

           alpha_trajc(:,:,nq) = 0.d0

           ! Calculate Nt along trajectory from values of q and N0

           IF(nq == 1) THEN
             qscalar_trajc(:,:,nq+5) = Ntcfix  ! Fixed number concentration for cloud droplets
           ELSE IF(nq == 3) THEN
             qscalar_trajc(:,:,nq+5) = 0.0  ! Fletcher equation for ice number concentration: TODO: need to determine this for LIN scheme
           ELSE
             CALL cal_Nt(ntrajcs,ntimes,rho_trajc,qscalar_trajc(:,:,nq),&
                         N0_trajc(:,:,nq),cx(nq),alpha_trajc(:,:,nq),   &
                         qscalar_trajc(:,:,nq+6))
           END IF

           ! Calculate mean-mass diameter along trajectories

           CALL cal_Dm(ntrajcs,ntimes,rho_trajc,qscalar_trajc(:,:,nq),  &
                       qscalar_trajc(:,:,nq+6),cx,Dm_trajc(:,:,nq))
         END DO
       ELSE IF(mphyopt >= 5 .and. mphyopt <= 7) THEN
         print*,'WSM6 not supported yet for DSD calculations'
       ELSE IF(mphyopt == 8) THEN
         DO nq=1,nscalarq
           ! Calculate alpha along trajectories

           IF(nq == 1) THEN
             alpha_trajc(:,:,nq) = 1.d0  ! Fixed alpha = 1 for cloud droplets
           ELSE IF(nq == 2) THEN
             alpha_trajc(:,:,nq) = dble(alpharfix)
           ELSE IF(nq == 3) THEN
             alpha_trajc(:,:,nq) = dble(alphaifix)
           ELSE IF(nq == 4) THEN
             alpha_trajc(:,:,nq) = dble(alphasfix)
           ELSE IF(nq == 5) THEN
             alpha_trajc(:,:,nq) = dble(alphagfix)
           ELSE IF(nq == 6) THEN
             alpha_trajc(:,:,nq) = dble(alphahfix)
           END IF

           ! Calculate N0 along trajectories

           IF(nq == 1) THEN
             ! Nt constant for cloud
             qscalar_trajc(:,:,nq+6) = 1.0e8
             !CALL cal_N0(ntrajcs,ntimes,rho_trajc,qscalar_trajc(:,:,nq),qscalar_trajc(:,:,nq+6),cx(nq),alpha_trajc(:,:,nq),N0_trajc(:,:,nq),N0eff_trajc(:,:,nq))
             ! Set N0 for cloud to zero for now for testing
             N0_trajc(:,:,nq) = 0.d0
           ELSE IF(nq == 2) THEN
             N0_trajc(:,:,nq) = dble(N0rfix)
           ELSE IF(nq == 3) THEN
             N0_trajc(:,:,nq) = 0.d0 ! Not applicable (Nt calculated from Cooper eqn)
           ELSE IF(nq == 4) THEN
             N0_trajc(:,:,nq) = dble(N0sfix)
           ELSE IF(nq == 5) THEN
             N0_trajc(:,:,nq) = dble(N0gfix)
           ELSE IF(nq == 6) THEN
             N0_trajc(:,:,nq) = dble(N0hfix)
           END IF

           ! Calculate Nt along trajectories

           IF(nq == 3) THEN
             DO k=1,ntimes
               DO i=1,ntrajcs
                 qscalar_trajc(i,k,nq+6) = 5.*exp(0.304*(273.15-max(233.,temp_trajc(i,k))))  ! Cooper eqn for Nt for ice
               END DO
             END DO
           ELSE IF (nq /= 1) THEN
             CALL cal_Nt(ntrajcs,ntimes,rho_trajc,qscalar_trajc(:,:,nq),&
                         N0_trajc(:,:,nq),cx(nq),alpha_trajc(:,:,nq),   &
                         qscalar_trajc(:,:,nq+6))
           END IF

           ! Calculate mean-mass diameter along trajectories

           CALL cal_Dm(ntrajcs,ntimes,rho_trajc,qscalar_trajc(:,:,nq),  &
                       qscalar_trajc(:,:,nq+6),cx(nq),Dm_trajc(:,:,nq))
         END DO

       ELSE IF(mphyopt == 9 .or. mphyopt == 10) THEN

         DO nq=1,nscalarq
           ! Calculate Dm along trajectories

           CALL cal_Dm(ntrajcs,ntimes,rho_trajc,qscalar_trajc(:,:,nq),  &
                       qscalar_trajc(:,:,nq+6),cx(nq),Dm_trajc(:,:,nq))

           ! Calculate alpha along trajectories

           IF(mphyopt == 9) THEN  ! Fixed alpha
             IF(nq == 1) THEN
               alpha_trajc(:,:,nq) = 1.d0  ! Fixed alpha = 1 for cloud droplets
             ELSE IF(nq == 2) THEN
               alpha_trajc(:,:,nq) = dble(alpharfix)
             ELSE IF(nq == 3) THEN
               alpha_trajc(:,:,nq) = dble(alphaifix)
             ELSE IF(nq == 4) THEN
               alpha_trajc(:,:,nq) = dble(alphasfix)
             ELSE IF(nq == 5) THEN
               alpha_trajc(:,:,nq) = dble(alphagfix)
             ELSE IF(nq == 6) THEN
               alpha_trajc(:,:,nq) = dble(alphahfix)
             END IF
           ELSE IF(mphyopt == 10) THEN ! Diagnostic alpha
             IF(nq == 1) THEN
               alpha_trajc(:,:,nq) = 1.d0  ! Fixed alpha = 1 for cloud droplets
             ELSE
               CALL diag_alpha(ntrajcs,ntimes,rho_trajc,varid_qscalar(nq), &
                               Dm_trajc(:,:,nq),alpha_trajc(:,:,nq))
             END IF
           END IF

           ! Calculate N0 (and "effective" N0) along trajectories

           CALL cal_N0(ntrajcs,ntimes,rho_trajc,qscalar_trajc(:,:,nq),    &
                       qscalar_trajc(:,:,nq+6),cx(nq),alpha_trajc(:,:,nq),&
                       N0_trajc(:,:,nq),N0eff_trajc(:,:,nq))

         END DO
       ELSE IF(mphyopt == 11) THEN
         DO nq=1,nscalarq
           ! Calculate Dm along trajectories
           CALL cal_Dm(ntrajcs,ntimes,rho_trajc,qscalar_trajc(:,:,nq),  &
                       qscalar_trajc(:,:,nq+6),cx(nq),Dm_trajc(:,:,nq))

           ! Calculate alpha along trajectories
           IF(nq >= 2) THEN
             CALL solve_alpha(ntrajcs,ntimes,rho_trajc,cx(nq),               &
                              qscalar_trajc(:,:,nq),qscalar_trajc(:,:,nq+6), &
                              qscalar_trajc(:,:,nq+11),alpha_trajc(:,:,nq))
           ELSE
             alpha_trajc(:,:,nq) = 1.d0  ! Fixed alpha for cloud droplets
           END IF

           ! Calculate N0 (and "effective" N0) along trajectories

           CALL cal_N0(ntrajcs,ntimes,rho_trajc,qscalar_trajc(:,:,nq),     &
                       qscalar_trajc(:,:,nq+6),cx(nq),alpha_trajc(:,:,nq), &
                       N0_trajc(:,:,nq),N0eff_trajc(:,:,nq))

         END DO
       ELSE
         WRITE(6,'(a)') 'This microphysics option is not supported yet.'
       END IF
     END IF ! DSD_flag == 1

     IF(mprate_flag == 1 .and. npnt /= nend) THEN
       varid_mpcool(1) = 'evapqc'
       varid_mpcool(2) = 'evapqr'
       varid_mpcool(3) = 'sublqi'
       varid_mpcool(4) = 'sublqs'
       varid_mpcool(5) = 'sublqg'
       varid_mpcool(6) = 'sublqh'
       varid_mpcool(7) = 'meltqi'
       varid_mpcool(8) = 'meltqs'
       varid_mpcool(9) = 'meltqg'
       varid_mpcool(10) = 'meltqh'

       DO nq=1,10

         CALL readvar3(nx,ny,nz, tem1, varid_mpcool(nq), varname(nq),   &
                 varunits(nq), time, runname, dirname, hdmpftrailer,3, 0, istatus)

         CALL edgfill(tem1,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,2,nz-1)
         CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem1,                      &
            xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, mpcool_trajc(:,:,nq) )

         ! Calculate cumulative quantities (total evap/subl/melt and cooling from
         ! these processes, per unit parcel mass)

         DO k=1,ntimes
           DO i=1,ntrajcs

             mpcoolproctot_trajc(i,k,nq) = mpcoolproctot_trajc(i,k,nq)+ &
                                           mpcool_trajc(i,k,nq)*tinc_calc

             IF(nq >= 1 .and. nq <= 2) THEN
               ! Instantaneous rate and Cumulative total evaporational cooling (J/kg/s,J/kg)
               mpcoolrate_trajc(i,k,nq) = LV*mpcool_trajc(i,k,nq)
               mpcooltot_trajc(i,k,nq) = mpcooltot_trajc(i,k,nq) +      &
                                      mpcoolrate_trajc(i,k,nq)*tinc_calc
             ELSE IF(nq >= 3 .and. nq <= 6) THEN
               ! Instantaneous rate and Cumulative total sublimational cooling (J/kg/s,J/kg)
               mpcoolrate_trajc(i,k,nq) = LS*mpcool_trajc(i,k,nq)
               mpcooltot_trajc(i,k,nq) = mpcooltot_trajc(i,k,nq) +      &
                                      mpcoolrate_trajc(i,k,nq)*tinc_calc
             ELSE
               !Instantaneous rate and Cumulative total melting cooling (J/kg/s,J/kg)
               mpcoolrate_trajc(i,k,nq) = LF*mpcool_trajc(i,k,nq)
               mpcooltot_trajc(i,k,nq) = mpcooltot_trajc(i,k,nq) +      &
                                      mpcoolrate_trajc(i,k,nq)*tinc_calc
             END IF

           END DO
         END DO

       END DO

       varid_mpheat(1)='condqc'
       varid_mpheat(2)='condqr'
       varid_mpheat(3)='nuclqi'
       varid_mpheat(4)='depoqi'
       varid_mpheat(5)='depoqs'
       varid_mpheat(6)='depoqg'
       varid_mpheat(7)='depoqh'
       varid_mpheat(8)='frzqci'
       varid_mpheat(9)='colqci'
       varid_mpheat(10)='colqcs'
       varid_mpheat(11)='colqcg'
       varid_mpheat(12)='colqch'
       varid_mpheat(13)='frzqrh'
       varid_mpheat(14)='colqri'
       varid_mpheat(15)='colqrs'
       varid_mpheat(16)='colqrg'
       varid_mpheat(17)='colqrh'

       DO nq=1,17

         CALL readvar3(nx,ny,nz, tem1, varid_mpheat(nq), varname(nq),   &
                 varunits(nq), time, runname, dirname, hdmpftrailer,3, 0, istatus)

         CALL edgfill(tem1,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,2,nz-1)
         CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem1,                      &
            xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, mpheat_trajc(:,:,nq) )

         ! Calculate cumulative quantities (total frz/depo/cond and heating from
         ! these processes, per unit parcel mass)

         DO k=1,ntimes
           DO i=1,ntrajcs

             mpheatproctot_trajc(i,k,nq) = mpheatproctot_trajc(i,k,nq)+ &
                                           mpheat_trajc(i,k,nq)*tinc_calc

             ! Instantaneous rate and Cumulative total condensational heating (J/kg/s,J/kg)
             IF(nq >= 1 .and. nq <= 2) THEN
               mpheatrate_trajc(i,k,nq) = LV*mpheat_trajc(i,k,nq)
               mpheattot_trajc(i,k,nq) = mpheattot_trajc(i,k,nq) +      &
                                      mpheatrate_trajc(i,k,nq)*tinc_calc
             ELSE IF(nq >= 3 .and. nq <= 7) THEN ! Instantaneous rate and Cumulative total depositional heating (J/kg)
               mpheatrate_trajc(i,k,nq) = LS*mpheat_trajc(i,k,nq)
               mpheattot_trajc(i,k,nq) = mpheattot_trajc(i,k,nq) +      &
                                      mpheatrate_trajc(i,k,nq)*tinc_calc
             ELSE ! Instantaneous rate and Cumulative total freezing heating (J/kg)
               mpheatrate_trajc(i,k,nq) = LF*mpheat_trajc(i,k,nq)
               mpheattot_trajc(i,k,nq) = mpheattot_trajc(i,k,nq) +      &
                                      mpheatrate_trajc(i,k,nq)*tinc_calc
             END IF

           END DO
         END DO


       END DO

       varid_mptrate = 'mptrat'

       CALL readvar3(nx,ny,nz, tem1, varid_mptrate, varname(1),         &
                 varunits(1), time, runname, dirname, hdmpftrailer, 3, 0, istatus)

       CALL edgfill(tem1,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,2,nz-1)
       CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem1,                        &
            xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, mptrate_trajc )


       ! Perform integration of thermodynamic energy (potential temperature)
       ! Equation if wanted.
       IF(therm_intg_flag == 1) THEN

         CALL hdfreadvar(nx,ny,nz,trim(hisfile),time,'pt', tem1, itmp )
         CALL hdfreadvar(nx,ny,nz,trim(hisfile),time,'p', tem2, itmp )

         DO k=1,nz-1
           DO j=1,ny-1
             DO i=1,nx-1
                 ! Inverse exner function in tem3
                 tem3(i,j,k) = (p0 / tem2(i,j,k)) ** rddcp
             END DO
           END DO
         END DO

         !WRITE(0,*) 'p at i=3,j=3,k=1 and 2',tem2(3,3,1),tem2(3,3,2)

         CALL edgfill(tem3,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
         CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem3,                      &
              xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, ppi_trajc )


         ! Save initial potential temperature and mixing ratio of each trajectory

         DO k=1,ntimes
           DO i=1,ntrajcs
             IF(npnt == npnt_trajc_start(i,k)) THEN
               ! Save initial potential temperature of each trajectory
               !WRITE(0,*) 'k,i,npnt,pt_trajc,ppi_trajc',k,i,npnt,pt_trajc(i,k),ppi_trajc(i,k)
               ptcalc_trajc(i,k) = pt_trajc(i,k)
               ptcalcmix_trajc(i,k) = pt_trajc(i,k)
               ptcalcmp_trajc(i,k) = pt_trajc(i,k)
             END IF
           END DO
         END DO

         DO k=1,ntimes
           DO i=1,ntrajcs
             IF(npnt == npnt_trajc_start(i,k)) THEN
               ! Save initial qv of each trajectory
               qvcalc_trajc(i,k) = qv_trajc(i,k)
               qvcalcmix_trajc(i,k) = qv_trajc(i,k)
               qvcalcmp_trajc(i,k) = qv_trajc(i,k)
             END IF
           END DO
         END DO

         !DTD: integrated equivalent potential temperature

         DO k=1,ntimes
           DO i=1,ntrajcs
              IF( .not. trajc_missing(npnt,i,k)) THEN
                temp1 = p_trajc(i,k)
                temp2 = ptcalc_trajc(i,k)
                temp3 = qvcalc_trajc(i,k)
                !WRITE(0,*) 'k,i,p_trajc,ptcalc_trajc,qvcalc_trajc',k,i,p_trajc(i,k),ptcalc_trajc(i,k),qvcalc_trajc(i,k)
                ptecalc_trajc(i,k) = f_pt2pte(temp1,temp2,temp3)
                temp1 = p_trajc(i,k)
                temp2 = ptcalcmix_trajc(i,k)
                temp3 = qvcalcmix_trajc(i,k)
                ptecalcmix_trajc(i,k) = f_pt2pte(temp1,temp2,temp3)
                temp1 = p_trajc(i,k)
                temp2 = ptcalcmp_trajc(i,k)
                temp3 = qvcalcmp_trajc(i,k)
                ptecalcmp_trajc(i,k) = f_pt2pte(temp1,temp2,temp3)
              END IF

    !         ptecalc_trajc(i,k) = f_pt2pte(p_trajc(i,k),ptcalc_trajc(i,k),qvcalc_trajc(i,k))
    !         ptecalcmix_trajc(i,k) = f_pt2pte(p_trajc(i,k),ptcalcmix_trajc(i,k),qvcalcmix_trajc(i,k))
    !         CALL pt2pte(1,1,1,1,1,1,1,1,1,p_trajc(i,k),ptcalc_trajc(i,k),     &
    !                     qvcalc_trajc(i,k),ptecalc_trajc(i,k))
    !         CALL pt2pte(1,1,1,1,1,1,1,1,1,p_trajc(i,k),ptcalcmix_trajc(i,k),     &
    !                     qvcalcmix_trajc(i,k),ptecalcmix_trajc(i,k))
           END DO
         END DO

         j=npnt

         IF(avg_trajc_flag == 1) THEN
           DO k=1,ntimes
             DO i=1,ntrajcs
               IF( .not. trajc_missing(npnt,i,k)) THEN
                 pt_trajcavg(npnt,k)=pt_trajcavg(npnt,k)+pt_trajc(i,k)/ntrajcsavg(npnt,k)
                 ptcalc_trajcavg(npnt,k)=ptcalc_trajcavg(npnt,k)+ptcalc_trajc(i,k)/ntrajcsavg(npnt,k)
                 ptcalcmp_trajcavg(npnt,k)=ptcalcmp_trajcavg(npnt,k)+ptcalcmp_trajc(i,k)/ntrajcsavg(npnt,k)
                 ptcalcmix_trajcavg(npnt,k)=ptcalcmix_trajcavg(npnt,k)+ptcalcmix_trajc(i,k)/ntrajcsavg(npnt,k)
                 qv_trajcavg(npnt,k)=qv_trajcavg(npnt,k)+qv_trajc(i,k)/ntrajcsavg(npnt,k)
                 qvcalc_trajcavg(npnt,k)=qvcalc_trajcavg(npnt,k)+qvcalc_trajc(i,k)/ntrajcsavg(npnt,k)
                 qvcalcmp_trajcavg(npnt,k)=qvcalcmp_trajcavg(npnt,k)+qvcalcmp_trajc(i,k)/ntrajcsavg(npnt,k)
                 qvcalcmix_trajcavg(npnt,k)=qvcalcmix_trajcavg(npnt,k)+qvcalcmix_trajc(i,k)/ntrajcsavg(npnt,k)
                 ptecalc_trajcavg(npnt,k)=ptecalc_trajcavg(npnt,k)+ptecalc_trajc(i,k)/ntrajcsavg(npnt,k)
                 ptecalcmp_trajcavg(npnt,k)=ptecalcmp_trajcavg(npnt,k)+ptecalcmp_trajc(i,k)/ntrajcsavg(npnt,k)
                 ptecalcmix_trajcavg(npnt,k)=ptecalcmix_trajcavg(npnt,k)+ptecalcmix_trajc(i,k)/ntrajcsavg(npnt,k)
               END IF
             END DO
           END DO

           DO k=1,ntimes
            !WRITE(nunit_thermintg(k),'(f8.1)') ttrajc(j)
             WRITE(nunit(k,8),'(f8.1,",",6(f10.2,","),29(e14.6,","))') ttrajc(j),  &
                   xtrajcavg(j,k),ytrajcavg(j,k),ztrajcavg(j,k), &
                   utrajcavg(j,k),vtrajcavg(j,k),wtrajcavg(j,k),    &
                   (qscalar_trajcavg(j,k,nq),nq=1,nscalar),   &
                   pt_trajcavg(j,k),ptcalcmix_trajcavg(j,k),ptcalcmp_trajcavg(j,k),ptcalc_trajcavg(j,k), &
                   qv_trajcavg(j,k),qvcalcmix_trajcavg(j,k),qvcalcmp_trajcavg(j,k),qvcalc_trajcavg(j,k), &
                   pte_trajcavg(j,k),ptecalcmix_trajcavg(j,k),ptecalcmp_trajcavg(j,k), &
                   ptecalc_trajcavg(j,k)
           ENDDO

         END IF

         ! Dump out the values to file here, since the integration is forward in time,
         ! so the new calculated values will be valid for the *following* time step.

         DO k=1,ntimes
           !WRITE(nunit_thermintg(k),'(f8.1)') ttrajc(j)
           DO i = 1,ntrajcs
             WRITE(nunit(k,7),'(f8.1,",",6(f10.2,","),12(e14.6,","))') ttrajc(j),  &
                   xtrajc(j,i,k),ytrajc(j,i,k),ztrajc(j,i,k), &
                   utrajc(i,k),vtrajc(i,k),wtrajc(i,k),    &
                   pt_trajc(i,k),ptcalcmix_trajc(i,k),ptcalcmp_trajc(i,k),ptcalc_trajc(i,k), &
                   qv_trajc(i,k),qvcalcmix_trajc(i,k),qvcalcmp_trajc(i,k),qvcalc_trajc(i,k), &
                   pte_trajc(i,k),ptecalcmix_trajc(i,k),ptecalcmp_trajc(i,k), &
                   ptecalc_trajc(i,k)
           ENDDO
         ENDDO

         ! Integrate theta conservation equation along trajectories

         ! Calculate the turbulent mixing term for potential temperature along the trajectories
         ! and add to the calculated ptcalc_trajc above

         varid = 'mix_pt'

         CALL readvar3(nx,ny,nz, tem1, varid, varname(1),               &
                       varunits(1),time,runname,dirname, hdmpftrailer,  &
                       3, 0, istatus)

         ! Calculate the coordinate transformation jacobian tem2
         DO k=1,nz-1
          DO j=1,ny-1
            DO i=1,nx-1
              tem2(i,j,k) = (zp(i,j,k+1)-zp(i,j,k))/(z(k+1)-z(k))
            END DO
          END DO
         END DO

         DO k=1,nz-1
           DO j=1,ny-1
             tem2(nx,j,k) = tem2(nx-1,j,k)
           END DO
         END DO

         DO k=1,nz-1
           DO i=1,nx
             tem2(i,ny,k) = tem2(i,ny-1,k)
           END DO
         END DO

         DO j=1,ny
           DO i=1,nx
             tem2(i,j,nz) = tem2(i,j,nz-1)
           END DO
         END DO

          DO k=1,nz-1
            DO j=1,ny-1
              DO i=1,nx-1
                tem1(i,j,k) = tem1(i,j,k)/(rhobar(i,j,k)*tem2(i,j,k))
              END DO
            END DO
          END DO

!         CALL edgfill(tem1,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
         CALL edgfill(tem1,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,2,nz-1)
         CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem1,                      &
              xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, ptmix_trajc )

         DO k=1,ntimes
           DO i=1,ntrajcs
             ptcalc_trajc(i,k) = ptcalc_trajc(i,k) + ptmix_trajc(i,k)*tinc_calc
             ! also keep track of theta change due just to mixing
             ptcalcmix_trajc(i,k) = ptcalcmix_trajc(i,k) + ptmix_trajc(i,k)*tinc_calc
           END DO
         END DO

         ! Calculate new potential temperature based on heating/cooling rates above

         DO k=1,ntimes
           DO i=1,ntrajcs
             DO nq=1,10
               ptcalc_trajc(i,k) = ptcalc_trajc(i,k) - ppi_trajc(i,k)*mpcoolrate_trajc(i,k,nq)/CP*tinc_calc
               ptcalcmp_trajc(i,k) = ptcalcmp_trajc(i,k) - ppi_trajc(i,k)*mpcoolrate_trajc(i,k,nq)/CP*tinc_calc
             END DO
             DO nq=1,17
               ptcalc_trajc(i,k) = ptcalc_trajc(i,k) + ppi_trajc(i,k)*mpheatrate_trajc(i,k,nq)/CP*tinc_calc
               ptcalcmp_trajc(i,k) = ptcalcmp_trajc(i,k) + ppi_trajc(i,k)*mpheatrate_trajc(i,k,nq)/CP*tinc_calc
             END DO
           END DO
         END DO

         ! Integrate qv conservation equation along trajectories

         ! Calculate the turbulent mixing term for qv along the trajectories
         ! and add to the calculated qvcalc_trajc above

         varid = 'mix_qv'

         CALL readvar3(nx,ny,nz, tem1, varid, varname(1),               &
                       varunits(1),time,runname, dirname, hdmpftrailer, &
                       3, 0, istatus)

          DO k=1,nz-1
            DO j=1,ny-1
              DO i=1,nx-1
                tem1(i,j,k) = tem1(i,j,k)/(rhobar(i,j,k)*tem2(i,j,k))
              END DO
            END DO
          END DO

!         IF(npnt == nstart) THEN
!         DO j=1,ny
!           DO i=1,nx
!             WRITE(0,*) 'mix_qv at k = 1: i,j,mix_qv',i,j,tem1(i,j,1)
!             WRITE(0,*) 'mix_qv at k = 2: i,j,mix_qv',i,j,tem1(i,j,2)
!           END DO
!         END DO
!         END IF

         ! Ensure bottom boundary is handled
!         DO i=1,nx
!           DO j=1,ny
!             ! Zero-gradient for mixing term appropriate?
!             tem1(i,j,1)=tem(i,j,2)
!           END DO
!         END DO

!         CALL edgfill(tem1,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
         CALL edgfill(tem1,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,2,nz-1)

         CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem1,                      &
                          xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes,     &
                          qvmix_trajc )

         DO k=1,ntimes
           DO i=1,ntrajcs
             qvcalc_trajc(i,k) = qvcalc_trajc(i,k) + qvmix_trajc(i,k)*tinc_calc
             ! also keep track of changes to qv due just to mixing
             qvcalcmix_trajc(i,k) = qvcalcmix_trajc(i,k) + qvmix_trajc(i,k)*tinc_calc
           END DO
         END DO

         ! Add contribution from microphysics

         DO k=1,ntimes
           DO i=1,ntrajcs
             DO nq=1,6
               qvcalc_trajc(i,k) = qvcalc_trajc(i,k) + mpcool_trajc(i,k,nq)*tinc_calc
               qvcalcmp_trajc(i,k) = qvcalcmp_trajc(i,k) + mpcool_trajc(i,k,nq)*tinc_calc
             END DO
             DO nq=1,7
               qvcalc_trajc(i,k) = qvcalc_trajc(i,k) - mpheat_trajc(i,k,nq)*tinc_calc
               qvcalcmp_trajc(i,k) = qvcalcmp_trajc(i,k) - mpheat_trajc(i,k,nq)*tinc_calc
             END DO
           END DO
         END DO

         ! Add contribution from "clipping" (reuse qmix_trajc)

!         varid = 'clipqv'

!         CALL readvar2(nx,ny,nz, tem1, varid, varname(1),     &
!                   varunits(1), time, runname, dirname, 3, 0, istatus)

!         CALL edgfill(tem1,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,2,nz-1)

!         CALL intrp_trajc(nx,ny,nz,xs,ys,zp, tem1,                   &
!              xtrajc1,ytrajc1,ztrajc1, ntrajcs, ntimes, qvmix_trajc )

!         DO k=1,ntimes
!           DO i=1,ntrajcs
!             qvcalc_trajc(i,k) = qvcalc_trajc(i,k) + qvmix_trajc(i,k)*tinc_calc
!             qvcalcmp_trajc(i,k) = qvcalcmp_trajc(i,k) + qvmix_trajc(i,k)*tinc_calc
!           END DO
!         END DO

       END IF ! therm_intg_flag == 1

     END IF ! mprate_flag == 1

!-----------------------------------------------------------------------
!  Write out data along trajectories
!-----------------------------------------------------------------------

    j = npnt
    DO k=1,ntimes
      !WRITE(nunit(k),'(f8.1)') ttrajc(j)
      DO i = 1,ntrajcs
        WRITE(nunit(k,1),'(f8.1,a,6(f10.2,","),37(e14.6,","))')        &
              ttrajc(j),',',                                           &
              xtrajc(j,i,k),ytrajc(j,i,k),ztrajc(j,i,k),               &
              utrajc(i,k),vtrajc(i,k),wtrajc(i,k),                     &
              pt_trajc(i,k),p_trajc(i,k),qv_trajc(i,k),                &
              ptprt_trajc(i,k),pprt_trajc(i,k),                        &
              vortx_trajc(i,k),vorty_trajc(i,k),vortz_trajc(i,k),      &
              vorts_trajc(i,k), vortc_trajc(i,k),                      &
              buoy_trajc(i,k),buoyq_trajc(i,k),frcs_trajc(i,k),        &
              upgrad_trajc(i,k),vpgrad_trajc(i,k),wpgrad_trajc(i,k),   &
              vortx_stch_trajc(i,k), vorty_stch_trajc(i,k),            &
              vorts_stch_trajc(i,k), vortc_stch_trajc(i,k),            &
              vortz_stch_trajc(i,k),                                   &
              vortx_tilt_trajc(i,k), vorty_tilt_trajc(i,k),            &
              vorts_tilt_trajc(i,k), vortc_tilt_trajc(i,k),            &
              vortz_tilt_trajc(i,k),                                   &
              vortx_mix_trajc(i,k), vorty_mix_trajc(i,k),              &
              vorts_mix_trajc(i,k), vortc_mix_trajc(i,k),              &
              vortz_mix_trajc(i,k),                                    &
              vortx_gen_trajc(i,k), vorty_gen_trajc(i,k),              &
              vorts_gen_trajc(i,k), vortc_gen_trajc(i,k),              &
              vorts_exhg_trajc(i,k), vortc_exhg_trajc(i,k)
      END DO
    END DO

    IF(mprate_flag == 1) THEN
      DO k=1,ntimes
        !WRITE(nunit(k,2),'(f8.1)') ttrajc(j)
        DO i = 1,ntrajcs
          WRITE(nunit(k,2),'(f8.1,a,6(f10.2,","),44(e14.6,","))')      &
                ttrajc(j),',',                                         &
                xtrajc(j,i,k),ytrajc(j,i,k),ztrajc(j,i,k),             &
                utrajc(i,k),vtrajc(i,k),wtrajc(i,k),                   &
                (qscalar_trajc(i,k,nq),nq=1,nscalar),                  &
                (mpcool_trajc(i,k,nq),nq=1,10),                        &
                (mpheat_trajc(i,k,nq),nq=1,17)
        ENDDO
      ENDDO

      DO k=1,ntimes
        !WRITE(nunit(k,3),'(f8.1)') ttrajc(j)
        DO i = 1,ntrajcs
          WRITE(nunit(k,3),'(f8.1,a,6(f10.2,","),27(e14.6,","))')      &
                ttrajc(j),',',                                         &
                xtrajc(j,i,k),ytrajc(j,i,k),ztrajc(j,i,k),             &
                utrajc(i,k),vtrajc(i,k),wtrajc(i,k),                   &
                (mpcoolproctot_trajc(i,k,nq),nq=1,10),                 &
                (mpheatproctot_trajc(i,k,nq),nq=1,17)
        ENDDO
      ENDDO

      DO k=1,ntimes
        !WRITE(nunit_mp3(k),'(f8.1)') ttrajc(j)
        DO i = 1,ntrajcs
          WRITE(nunit(k,4),'(f8.1,a,6(f10.2,","),27(e14.6,","))')      &
                ttrajc(j),',',                                         &
                xtrajc(j,i,k),ytrajc(j,i,k),ztrajc(j,i,k),             &
                utrajc(i,k),vtrajc(i,k),wtrajc(i,k),                   &
                (mpcoolrate_trajc(i,k,nq),nq=1,10),                    &
                (mpheatrate_trajc(i,k,nq),nq=1,17)
        ENDDO
      ENDDO

      DO k=1,ntimes
        !WRITE(nunit_mp4(k),'(f8.1)') ttrajc(j)
        DO i = 1,ntrajcs
          WRITE(nunit(k,5),'(f8.1,a,6(f10.2,","),27(e14.6,","))')      &
                ttrajc(j),',',                                         &
                xtrajc(j,i,k),ytrajc(j,i,k),ztrajc(j,i,k),             &
                utrajc(i,k),vtrajc(i,k),wtrajc(i,k),                   &
                (mpcooltot_trajc(i,k,nq),nq=1,10),                     &
                (mpheattot_trajc(i,k,nq),nq=1,17)
        ENDDO
      ENDDO
    END IF

    IF(DSD_flag == 1) THEN
      DO k=1,ntimes
        !WRITE(nunit_mp5(k),'(f8.1)') ttrajc(j)
        DO i = 1,ntrajcs
          WRITE(nunit(k,6),'(f8.1,a,6(f10.2,","),41(e14.6,","))')      &
                ttrajc(j),',',                                         &
                xtrajc(j,i,k),ytrajc(j,i,k),ztrajc(j,i,k),             &
                utrajc(i,k),vtrajc(i,k),wtrajc(i,k),                   &
                (qscalar_trajc(i,k,nq),nq=1,nscalar),                  &
                (N0_trajc(i,k,nq),nq=1,nscalarq),                      &
                (N0eff_trajc(i,k,nq),nq=1,nscalarq),                   &
                (alpha_trajc(i,k,nq),nq=1,nscalarq),                   &
                (Dm_trajc(i,k,nq),nq=1,nscalarq)
        ENDDO
      ENDDO
    END IF

  END DO

!-----------------------------------------------------------------------
!  Close the trajectory file
!-----------------------------------------------------------------------

  DO i=1,8
    DO k=1,ntimes
      CALL flush(nunit(k,i))
      CLOSE(UNIT=nunit(k,i))
      CALL retunit( nunit(k,i))
    ENDDO
  END DO
  PRINT *, 'done writing trajectory data'

  STOP

  100 WRITE(6,'(a)') 'Error reading NAMELIST file. Program ARPSTRAJC stopped.'
  STOP

END PROGRAM arpscalctrajc

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
    WRITE (6,*) "get_gridinfo_from_hdf: ERROR opening ",                &
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
  END DO
  END DO
  END DO

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
          aout(i,j,k)=      wgtx(i) *ain(ix(i)  ,j,k)                   &
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
          aout(i,j,k)=      wgty(j) *ain(i,jy(j)  ,k)                   &
                      +(1.0-wgty(j))*ain(i,jy(j)+1,k)
        END DO
      END DO
    END DO

  RETURN
END SUBROUTINE intrpy3d

SUBROUTINE intrpxy3d(ain,nx,is,ie, ny,js,je, nz,ks,ke,                  &
                     wgtx,ix,wgty,jy,                                   &
                     aout,nx1,is1,ie1, ny1,js1,je1,                     &
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

  CALL intrpx3d(ain(1,1,k),nx,is,ie, ny,js,je, 1,1,1,                   &
                wgtx,ix,temx1yz,nx1,is1,ie1)

  CALL intrpy3d(temx1yz,nx1,is1,ie1, ny,js,je, 1,1,1,                   &
                wgty,jy,aout(1,1,k),ny1,js1,je1)

  END DO

  CALL edgfill(aout,1,nx1,is1,ie1,1,ny1,js1,je1,ks,ke,ks,ke)

  RETURN
END SUBROUTINE intrpxy3d

INTEGER FUNCTION get_ktrajc(z, zs, nz)

  IMPLICIT NONE
  INTEGER :: nz, k
  REAL :: z, zs(nz)

  IF( z < zs(1) ) then
    get_ktrajc = 1
    RETURN
  END IF

  IF( z >= zs(nz) ) then
    get_ktrajc = nz-1
    RETURN
  END IF

  DO k=1,nz-1

    IF( z >= zs(k) .and. z < zs(k+1) ) then
      get_ktrajc = k
      EXIT
    END IF

  END DO

  RETURN

END FUNCTION get_ktrajc

FUNCTION gamma(xx)

!  Modified from "Numerical Recipes"

  IMPLICIT NONE

! PASSING PARAMETERS:
  DOUBLE PRECISION, INTENT(IN) :: xx

! LOCAL PARAMETERS:
  DOUBLE PRECISION  :: gamma
  INTEGER  :: j
  DOUBLE PRECISION  :: ser,stp,tmp,x,y,cof(6)


  SAVE cof,stp
  DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,                &
       24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,   &
       -.5395239384953d-5,2.5066282746310005d0/
  x=xx
  y=x
  tmp=x+5.5d0
  tmp=(x+0.5d0)*log(tmp)-tmp
  ser=1.000000000190015d0
! do j=1,6   !original
  do j=1,4
!!do j=1,3   !gives result to within ~ 3 %
     y=y+1.d0
     ser=ser+cof(j)/y
  enddo
  gamma=tmp+log(stp*ser/x)
  gamma= exp(gamma)

END FUNCTION gamma


SUBROUTINE cal_Nt(ntrajcs,ntimes,rho_trajc,q_trajc,N0_trajc,cx,alpha_trajc,Nt_trajc)

!
!-----------------------------------------------------------------------
!  PURPOSE:  Calculates number concentration on a trajectory
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Dawson
!  (11/24/2007)
!
!  MODIFICATION HISTORY:
!  (04/21/2008) Dan Dawson
!  Modified to be consistent with the main dsdlib package (changed some
!  variables to double precision)
!-----------------------------------------------------------------------
!  Variable Declarations:
!-----------------------------------------------------------------------
!

  INTEGER :: ntrajcs, ntimes
  REAL :: rho_trajc(ntrajcs,ntimes),q_trajc(ntrajcs,ntimes)
  REAL*8 :: N0_trajc(ntrajcs,ntimes),alpha_trajc(ntrajcs,ntimes)
  REAL :: cx
  REAL :: Nt_trajc(ntrajcs,ntimes)
  REAL*8 :: gamma1,gamma4

  REAL*8:: gamma

  INTEGER i,k

  DO k=1,ntimes
    DO i=1,ntrajcs
       gamma1 = gamma(1.d0+alpha_trajc(i,k))
       gamma4 = gamma(4.d0+alpha_trajc(i,k))

       !print*,'gamma1,gamma4,cx',gamma1,gamma4,cx

       Nt_trajc(i,k) = sngl((N0_trajc(i,k)*gamma1)**(3.d0/(4.d0+alpha_trajc(i,k)))*   &
                         ((gamma1/gamma4)*dble(rho_trajc(i,k))* &
                         dble(q_trajc(i,k))/dble(cx))**((1.d0+alpha_trajc(i,k))/(4.d0+alpha_trajc(i,k))))
    END DO
  END DO

END SUBROUTINE cal_Nt

SUBROUTINE cal_N0(ntrajcs,ntimes,rho_trajc,q_trajc,Nt_trajc,cx,alpha_trajc,N0_trajc,N0eff_trajc)

!
!-----------------------------------------------------------------------
!  PURPOSE:  Calculates intercept parameter and "effective" intercept parameter
!            on a trajectory
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Dawson
!  (11/24/2007)
!
!  MODIFICATION HISTORY:
!  (04/21/2008) Dan Dawson
!  Modified to be consistent with the main dsdlib package (changed some
!  variables to double precision)
!-----------------------------------------------------------------------
!  Variable Declarations:
!-----------------------------------------------------------------------
!

  INTEGER :: ntrajcs, ntimes
  REAL :: rho_trajc(ntrajcs,ntimes),q_trajc(ntrajcs,ntimes)
  REAL*8 :: N0_trajc(ntrajcs,ntimes),alpha_trajc(ntrajcs,ntimes),N0eff_trajc(ntrajcs,ntimes)
  REAL :: cx
  REAL :: Nt_trajc(ntrajcs,ntimes)
  REAL*8 :: gamma1, gamma4

  REAL*8:: gamma

  REAL*8:: lamda

  INTEGER i,k

  DO k=1,ntimes
    DO i=1,ntrajcs

      gamma1 = gamma(1.d0+alpha_trajc(i,k))
      gamma4 = gamma(4.d0+alpha_trajc(i,k))

      IF(rho_trajc(i,k) > 0.0 .and. q_trajc(i,k) > 0.0) THEN
        lamda = ((gamma4/gamma1)*dble(cx)*dble(Nt_trajc(i,k))/(dble(rho_trajc(i,k))*  &
                dble(q_trajc(i,k))))**(1.d0/3.d0)
      ELSE
        lamda = 0.d0
      END IF

      N0_trajc(i,k) = dble(Nt_trajc(i,k))*lamda**(0.5d0*(1.d0+alpha_trajc(i,k)))*     &
                      (1.d0/gamma1)*lamda**(0.5d0*(1.d0+alpha_trajc(i,k)))

      IF(lambda /= 0.d0) THEN
        N0eff_trajc(i,k) = N0_trajc(i,k)*(((4.d0+alpha_trajc(i,k))/lamda)**alpha_trajc(i,k))*        &
                           gamma4*(128.d0/3.d0)/((4.d0+alpha_trajc(i,k))**(4.d0+alpha_trajc(i,k)))
      ELSE
        N0eff_trajc(i,k) = 0.d0
      END IF

    END DO
  END DO

END SUBROUTINE cal_N0

SUBROUTINE cal_Dm(ntrajcs,ntimes,rho_trajc,q_trajc,Nt_trajc,cx,Dm_trajc)
!
!-----------------------------------------------------------------------
!  PURPOSE:  Calculates mean-mass diameter on a trajectory
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Dawson
!  (11/24/2007)
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!  Variable Declarations:
!-----------------------------------------------------------------------
!

  INTEGER :: ntrajcs, ntimes
  REAL :: rho_trajc(ntrajcs,ntimes),q_trajc(ntrajcs,ntimes)
  REAL :: cx
  REAL :: Nt_trajc(ntrajcs,ntimes),Dm_trajc(ntrajcs,ntimes)

  INTEGER i,k

  DO k=1,ntimes
    DO i=1,ntrajcs
      IF(Nt_trajc(i,k) /= 0.0) THEN
        Dm_trajc(i,k) = (rho_trajc(i,k)*q_trajc(i,k)/(cx*Nt_trajc(i,k)))**(1./3.)
      ELSE
        Dm_trajc(i,k) = 0.0
      END IF
    END DO
  END DO

END SUBROUTINE cal_Dm

SUBROUTINE diag_alpha(ntrajcs,ntimes,rho_trajc,varid_qscalar,Dm_trajc,alpha_trajc)
!
!-----------------------------------------------------------------------
!  PURPOSE:  Calculates shape parameter alpha on a trajectory
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Dawson
!  (11/24/2007)
!
!  MODIFICATION HISTORY:
!  (04/21/2008) Dan Dawson
!  Modified to be consistent with the main dsdlib package (changed some
!  variables to double precision)
!-----------------------------------------------------------------------
!  Variable Declarations:
!-----------------------------------------------------------------------
!

  INTEGER :: ntrajcs, ntimes
  REAL :: rho_trajc(ntrajcs,ntimes),Dm_trajc(ntrajcs,ntimes)
  REAL*8 :: alpha_trajc(ntrajcs,ntimes)
  REAL*8 :: ddalpha
  REAL*8 :: ddm

  REAL*8 :: diagAlpha

  CHARACTER (LEN=6) :: varid_qscalar

  INTEGER i,k,nq

  IF(varid_qscalar == 'qr') THEN
    nq = 1
  ELSE IF(varid_qscalar == 'qi') THEN
    nq = 2
  ELSE IF(varid_qscalar == 'qs') THEN
    nq = 3
  ELSE IF(varid_qscalar == 'qg') THEN
    nq = 4
  ELSE IF(varid_qscalar == 'qh') THEN
    nq = 5
  END IF

  DO k=1,ntimes
    DO i=1,ntrajcs
      ddm = dble(Dm_trajc(i,k))
      alpha_trajc(i,k) = diagAlpha(ddm,nq)
    END DO
  END DO

END SUBROUTINE diag_alpha

FUNCTION diagAlpha(Dm,x)

  IMPLICIT NONE

  integer :: x
  real*8  :: diagAlpha,Dm
  real, dimension(5) :: c1,c2,c3,c4
  real, parameter    :: pi = 3.14159265
  real*8, parameter  :: alphaMAX= 80.d0
  data c1 /19.0, 12.0, 4.5, 5.5, 3.7/
  data c2 / 0.6,  0.7, 0.5, 0.7, 0.3/
  data c3 / 1.8,  1.7, 5.0, 4.5, 9.0/
  data c4 /17.0, 11.0, 5.5, 8.5, 6.5/
  diagAlpha= c1(x)*tanh(c2(x)*(1.e3*Dm-c3(x)))+c4(x)
  if (x==5.and.Dm>0.008) diagAlpha= 1.e3*Dm-2.6
  diagAlpha= min(diagAlpha, alphaMAX)

END function diagAlpha

SUBROUTINE solve_alpha(ntrajcs,ntimes,rho_trajc,cx,q_trajc,Nt_trajc,Z_trajc,alpha_trajc)
!
!-----------------------------------------------------------------------
!  PURPOSE:  Calculates shape parameter alpha on a trajectory
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Dawson
!  (11/24/2007)
!
!  MODIFICATION HISTORY:
!  (04/21/2008) Dan Dawson
!  Modified to be consistent with the main dsdlib package (changed some
!  variables to double precision)
!-----------------------------------------------------------------------
!  Variable Declarations:
!-----------------------------------------------------------------------
!

  INTEGER :: ntrajcs, ntimes
  REAL :: rho_trajc(ntrajcs,ntimes),q_trajc(ntrajcs,ntimes),Nt_trajc(ntrajcs,ntimes),Z_trajc(ntrajcs,ntimes)
  REAL*8 :: alpha_trajc(ntrajcs,ntimes)

  REAL*8 :: solveAlpha
  REAL*8 :: dsA

  REAL :: cx

  INTEGER i,k

  DO k=1,ntimes
    DO i=1,ntrajcs
      IF(q_trajc(i,k) > 0.0 .and. Nt_trajc(i,k) > 0.0 .and. Z_trajc(i,k) > 0.0) THEN
        alpha_trajc(i,k) = solveAlpha(q_trajc(i,k),Nt_trajc(i,k),Z_trajc(i,k),cx,rho_trajc(i,k))
      ELSE
        alpha_trajc(i,k) = 0.d0
      END IF
    END DO
  END DO


END SUBROUTINE solve_alpha

FUNCTION solveAlpha(Q,N,Z,Cx,rho)

 IMPLICIT NONE

! PASSING PARAMETERS:
  real, INTENT(IN) :: Q, N, Z, Cx, rho

! LOCAL PARAMETERS:
  real*8 :: solveAlpha
  real   :: a,g,a1,g1,g2,tmp1
  integer :: i
  real, parameter :: alphaMax= 40.
  real, parameter :: epsQ    = 1.e-14
  real, parameter :: epsN    = 1.e-3
  real, parameter :: epsZ    = 1.e-32

!  Q         mass mixing ratio
!  N         total concentration
!  Z         reflectivity
!  Cx        (pi/6)*RHOx
!  rho       air density
!  a         alpha (returned as solveAlpha)
!  g         function g(a)= [(6+a)(5+a)(4+a)]/[(3+a)(2+a)(1+a)],
!              where g = (Cx/(rho*Q))**2.*(Z*N)


  if (Q==0. .or. N==0. .or. Z==0. .or. Cx==0. .or. rho==0.) then
  ! For testing/debugging only; this module should never be called
  ! if the above condition is true.
    print*,'*** STOPPED in MODULE ### solveAlpha *** '
    print*,'*** : ',Q,N,Z,Cx*1.9099,rho
    stop
  endif

  IF (Q>epsQ .and. N>epsN .and. Z>epsZ ) THEN

     tmp1= Cx/(rho*Q)
     g   = tmp1*Z*tmp1*N    ! g = (Z*N)*[Cx / (rho*Q)]^2

 !Note: The above order avoids OVERFLOW, since tmp1*tmp1 is very large

!----------------------------------------------------------!
! !Solve alpha numerically: (brute-force; for testing only)
!      a= 0.
!      g2= 999.
!      do i=0,4000
!         a1= i*0.01
!         g1= (6.+a1)*(5.+a1)*(4.+a1)/((3.+a1)*(2.+a1)*(1.+a1))
!         if(abs(g-g1)<abs(g-g2)) then
!            a = a1
!            g2= g1
!         endif
!      enddo
!----------------------------------------------------------!

!Piecewise-polynomial approximation of g(a) to solve for a:
     if (g>=20.) then
       a= 0.
     else
       g2= g*g
       if (g<20.  .and.g>=13.31) a= 3.3638e-3*g2 - 1.7152e-1*g + 2.0857e+0
       if (g<13.31.and.g>=7.123) a= 1.5900e-2*g2 - 4.8202e-1*g + 4.0108e+0
       if (g<7.123.and.g>=4.200) a= 1.0730e-1*g2 - 1.7481e+0*g + 8.4246e+0
       if (g<4.200.and.g>=2.946) a= 5.9070e-1*g2 - 5.7918e+0*g + 1.6919e+1
       if (g<2.946.and.g>=1.793) a= 4.3966e+0*g2 - 2.6659e+1*g + 4.5477e+1
       if (g<1.793.and.g>=1.405) a= 4.7552e+1*g2 - 1.7958e+2*g + 1.8126e+2
       if (g<1.405.and.g>=1.230) a= 3.0889e+2*g2 - 9.0854e+2*g + 6.8995e+2
       if (g<1.230) a= alphaMax
     endif

     solveAlpha= max(0.,min(a,alphaMax))

  ELSE

     solveAlpha= 0.

  ENDIF

END FUNCTION solveAlpha

SUBROUTINE intrp_trajc(nx,ny,nz,xs,ys,zp, u,                &
           xtrajc,ytrajc,ztrajc, ntrajcs,ntimes, utrajc )

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
  REAL :: u(nx,ny,nz) ! input - field to be interpolated to trajectory points

  INTEGER :: ntrajcs, ntimes

  REAL :: xtrajc(ntrajcs,ntimes),ytrajc(ntrajcs,ntimes),ztrajc(ntrajcs,ntimes)
  REAL :: utrajc(ntrajcs,ntimes)  ! output - u interpolated to trajectory points

  INTEGER :: get_ktrajc

  REAL :: dx, dy, utrajc1,vtrajc1,wtrajc1
  REAL :: uy1z1,uy2z1,uy1z2,uy2z2,vy1z1,vy2z1,vy1z2,vy2z2,wy1z1,wy2z1,wy1z2,wy2z2
  REAL :: uz1,uz2,vz1,vz2,wz1,wz2
  REAL :: missing_value

  INTEGER :: ntrajc,itrajc1,jtrajc1,ktrajc1

  REAL :: xweight1,yweight1,zweight1
!

  REAL :: zs1d(nz),zy1(nz),zy2(nz) ! automatic work array
  INTEGER :: i,k,l

  dx = xs(2)-xs(1)
  dy = ys(2)-ys(1)

  missing_value = -99999.0

  DO k=1,ntimes
  DO i= 1,ntrajcs

    IF( xtrajc(i,k) == missing_value .or.      &
        ytrajc(i,k) == missing_value .or.      &
        ztrajc(i,k) == missing_value ) then

      utrajc(i,k) = missing_value

    ELSE

    itrajc1 = MAX(1, MIN(nx-2, INT((xtrajc(i,k)-xs(1))/dx)+1 ))
    jtrajc1 = MAX(1, MIN(ny-2, INT((ytrajc(i,k)-ys(1))/dy)+1 ))
!    ktrajc1 = get_ktrajc(ztrajc(i,k),zs,nz-1)

    xweight1 = (xtrajc(i,k)-xs(itrajc1))/dx
    yweight1 = (ytrajc(i,k)-ys(jtrajc1))/dy

    DO l=1,nz-1

      zy1(l) = (1.0-xweight1)*(zp(itrajc1  ,jtrajc1  ,l)+zp(itrajc1  ,jtrajc1  ,l+1))*0.5 &
                   +xweight1 *(zp(itrajc1+1,jtrajc1  ,l)+zp(itrajc1+1,jtrajc1  ,l+1))*0.5
      zy2(l) = (1.0-xweight1)*(zp(itrajc1  ,jtrajc1+1,l)+zp(itrajc1  ,jtrajc1+1,l+1))*0.5 &
                   +xweight1 *(zp(itrajc1+1,jtrajc1+1,l)+zp(itrajc1+1,jtrajc1+1,l+1))*0.5

      zs1d(l) = ( 1.0-yweight1 )*zy1(l) + yweight1*zy2(l)

    END DO

    ktrajc1 = get_ktrajc(ztrajc(i,k),zs1d,nz-1)

    uy1z1 = (1.0-xweight1)*u(itrajc1,jtrajc1  ,ktrajc1  )+xweight1*u(itrajc1+1,jtrajc1  ,ktrajc1  )
    uy2z1 = (1.0-xweight1)*u(itrajc1,jtrajc1+1,ktrajc1  )+xweight1*u(itrajc1+1,jtrajc1+1,ktrajc1  )
    uy1z2 = (1.0-xweight1)*u(itrajc1,jtrajc1  ,ktrajc1+1)+xweight1*u(itrajc1+1,jtrajc1  ,ktrajc1+1)
    uy2z2 = (1.0-xweight1)*u(itrajc1,jtrajc1+1,ktrajc1+1)+xweight1*u(itrajc1+1,jtrajc1+1,ktrajc1+1)

    uz1 = ( 1.0-yweight1 )*uy1z1 + yweight1*uy2z1
    uz2 = ( 1.0-yweight1 )*uy1z2 + yweight1*uy2z2

    zweight1 = (ztrajc(i,k)-zs1d(ktrajc1))/(zs1d(ktrajc1+1)-zs1d(ktrajc1))
    utrajc(i,k) = (1.0-zweight1)*uz1+zweight1*uz2

    END IF
  END DO
  END DO

  RETURN
END SUBROUTINE intrp_trajc

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
!
!  Beginning of executable code...
!
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
  END IF

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

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE HDFREADXYZ                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE hdfreadxyz(nx,ny,nz,x,y,z, filename, time)

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
!-----------------------------------------------------------------------
!  Variable Declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: nx,ny,nz

  CHARACTER (LEN=*) :: filename

  REAL :: x(nx)           ! x coord.
  REAL :: y(ny)           ! y coord.
  REAL :: z(nz)           ! z coord.

  REAL :: time
  INTEGER :: ireturn

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

  INTEGER :: istat, sd_id
  INTEGER :: varflg, istatus

!-----------------------------------------------------------------------
!  Include files:
!-----------------------------------------------------------------------

! INCLUDE 'globcst.inc'
! INCLUDE 'grid.inc'        ! Grid parameters
! INCLUDE 'indtflg.inc'

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

  CALL hdfrd1d(sd_id,"x",nx,x,istat)
  IF (istat /= 0) GO TO 110
  CALL hdfrd1d(sd_id,"y",ny,y,istat)
  IF (istat /= 0) GO TO 110
  CALL hdfrd1d(sd_id,"z",nz,z,istat)
  IF (istat /= 0) GO TO 110

!-----------------------------------------------------------------------
!  Friendly exit message
!-----------------------------------------------------------------------

  930   CONTINUE

  WRITE(6,'(/a,F8.1,a/)')                                               &
  ' Data at time=', time/60,' (min) were successfully read.'

  ireturn = 0

  GO TO 130

!-----------------------------------------------------------------------
!  Error during read
!-----------------------------------------------------------------------

  110   CONTINUE
  WRITE(6,'(/a/)') ' Error reading data in HDFREAD'
  ireturn=1

  130   CONTINUE

  CALL hdfclose(sd_id,istat)

  IF (ireturn == 0) THEN
    IF (istat == 0) THEN
      WRITE(6,'(/a/a)') &
      "HDFREADXYZ: Successfully read ", trim(filename)
    ELSE
      WRITE(6,'(/a,i3,a/,a)') &
      "HDFREADXYZ: ERROR (status=", istat, ") closing ", trim(filename)
    END IF
  END IF

  RETURN
END SUBROUTINE hdfreadxyz

SUBROUTINE cal_xvort(v,w,y,zp,nx,ny,nz, xvort)

!
! Rewritten (DTD 10-10-05) to include terrain-following case
!

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: xvort(nx,ny,nz)
  REAL :: v(nx,ny,nz), w(nx,ny,nz)
  REAL :: y(ny),zp(nx,ny,nz)
  REAL :: zweightl,zweightr,zs(nx,ny,nz)
  REAL :: ws(nx,ny,nz),wr,wl
  REAL :: dwdy, dvdz


  INTEGER :: i,j,k,kl,k1,k2
  REAL :: factor
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO j=2,ny-2
    DO i=1,nx-1
      DO k=2,nz-2
        ! Determine value of z and w at center scalar point

        zs(i,j,k)=0.5*(zp(i,j,k)+zp(i,j,k+1))
        ws(i,j,k)=0.5*(w(i,j,k)+w(i,j,k+1))
      END DO
    END DO
  END DO

  DO j=2,ny-2
    DO i=1,nx-1
      DO k=2,nz-2
        k1=k
        k2=k

        ! Determine k levels for vertical interpolation on either side of center point
        IF(zs(i,j,k) < zp(i,j-1,k)) THEN
          DO WHILE (zs(i,j,k) < zp(i,j-1,k1))
            k1=k1-1
          END DO
        ELSE IF(zs(i,j,k) > zp(i,j-1,k+1)) THEN
          DO WHILE (zs(i,j,k) > zp(i,j-1,k1+1))
            k1=k1+1
          END DO
        END IF

        IF(zs(i,j,k) < zp(i,j+1,k)) THEN
          DO WHILE (zs(i,j,k) < zp(i,j+1,k2))
            k2=k2-1
          END DO
        ELSE IF(zs(i,j,k) > zp(i,j+1,k+1)) THEN
          DO WHILE (zs(i,j,k) > zp(i,j+1,k2+1))
            k2=k2+1
          END DO
        END IF

        ! Precalculate left-side and right-side interpolated w values

        IF(k1 >= 2 .and. k1 <= nz-2) THEN
          zweightl=(zs(i,j,k)-zp(i,j-1,k1))/(zp(i,j-1,k1+1)-zp(i,j-1,k1))
          wl=(1.0-zweightl)*zp(i,j-1,k1)+zweightl*zp(i,j-1,k1+1)
        END IF

        IF(k2 >= 2 .and. k2 <= nz-2) THEN
          zweightr=(zs(i,j,k)-zp(i,j+1,k2))/(zp(i,j+1,k2+1)-zp(i,j+1,k2))
          wr=(1.0-zweightr)*zp(i,j+1,k2)+zweightr*zp(i,j+1,k2+1)
        END IF


        if( k == 2 ) then
          kl = 2
          factor = 2.0
        else
          kl = k-1
          factor = 1.0
        END IF


        ! Make sure the k levels are not below ground or above the model top
        IF((k1 < 2 .and. k2 < 2) .or. (k1 > nz-2 .and. k2 > nz-2)) THEN
          ! No calculation can be performed!  Set vorticity to zero at this point.
          xvort(i,j,k) = 0.0
        ELSE IF(k1 < 2) THEN       ! Do a right-sided difference
          dwdy=(wr-ws(i,j,k))/(y(j+1)-y(j))
        ELSE IF(k2 < 2) THEN       ! Do a left-sided difference
          dwdy=(ws(i,j,k)-wl)/(y(j)-y(j-1))
        ELSE                  ! Do a normal centered difference
          dwdy=(wr-wl)/(y(j+1)-y(j-1))
        END IF

        if( k == 2 ) then
          kl = 2
          factor = 2.0
        else
          kl = k-1
          factor = 1.0
        END IF

        ! Calculate dvdz (there will be some error here due to assumption of uniform grid in the vertical)

        dvdz=0.25*factor*(v(i,j,k+1)+v(i,j+1,k+1)-v(i,j,kl )-v(i,j+1,kl ))/(zs(i,j,k+1)-zs(i,j,k))

        ! Calculate xvort

        xvort(i,j,k)=dwdy-dvdz
      END DO
    END DO
  END DO

  DO  k=2,nz-2
    DO i=1,nx-1
      xvort(i,   1,k)=xvort(i,   2,k)
      xvort(i,ny-1,k)=xvort(i,ny-2,k)
    END DO
  END DO

  DO j=1,ny-1
    DO i=1,nx-1
      xvort(i,j,   1)=xvort(i,j,   2)
      xvort(i,j,nz-1)=xvort(i,j,nz-2)
    END DO
  END DO

  RETURN
END SUBROUTINE cal_xvort

SUBROUTINE cal_yvort(u,w,x,zp,nx,ny,nz, yvort)

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: yvort(nx,ny,nz)
  REAL :: zp(nx,ny,nz)
  REAL :: u(nx,ny,nz), w(nx,ny,nz)
  REAL :: x(ny),zs(nx,ny,nz)
  REAL :: wl,wr,ws(nx,ny,nz)
  REAL :: dwdx,dudz,zweightl,zweightr
  INTEGER :: k1,k2

  INTEGER :: i,j,k, kl
  REAL :: factor
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!


  DO j=1,ny-1
    DO i=2,nx-2
      DO k=2,nz-2
        ! Determine value of z and w at center scalar point

        zs(i,j,k)=0.5*(zp(i,j,k)+zp(i,j,k+1))
        ws(i,j,k)=0.5*(w(i,j,k)+w(i,j,k+1))
      END DO
    END DO
  END DO


  DO j=1,ny-1
    DO i=2,nx-2
      DO k=2,nz-2
        k1=k
        k2=k

        ! Determine k levels for vertical interpolation on either side of center point
        IF(zs(i,j,k) < zp(i-1,j,k)) THEN
          DO WHILE (zs(i,j,k) < zp(i-1,j,k1))
            k1=k1-1
          END DO
        ELSE IF(zs(i,j,k) > zp(i-1,j,k+1)) THEN
          DO WHILE (zs(i,j,k) > zp(i-1,j,k1+1))
            k1=k1+1
          END DO
        END IF

        IF(zs(i,j,k) < zp(i+1,j,k)) THEN
          DO WHILE (zs(i,j,k) < zp(i+1,j,k2))
            k2=k2-1
          END DO
        ELSE IF(zs(i,j,k) > zp(i+1,j,k+1)) THEN
          DO WHILE (zs(i,j,k) > zp(i+1,j,k2+1))
            k2=k2+1
          END DO
        END IF

        ! Precalculate left-side and right-side interpolated w values

        IF(k1 >= 2 .and. k1 <= nz-2) THEN
          zweightl=(zs(i,j,k)-zp(i-1,j,k1))/(zp(i-1,j,k1+1)-zp(i-1,j,k1))
          wl=(1.0-zweightl)*zp(i-1,j,k1)+zweightl*zp(i-1,j,k1+1)
        END IF

        IF(k2 >= 2 .and. k2 <= nz-2) THEN
          zweightr=(zs(i,j,k)-zp(i+1,j,k2))/(zp(i+1,j,k2+1)-zp(i+1,j,k2))
          wr=(1.0-zweightr)*zp(i+1,j,k2)+zweightr*zp(i+1,j,k2+1)
        END IF

        if( k == 2 ) then
          kl = 2
          factor = 2.0
        else
          kl = k-1
          factor = 1.0
        END IF

        ! Make sure the k levels are not below ground or above the model top
        IF((k1 < 2 .and. k2 < 2) .or. (k1 > nz-2 .and. k2 > nz-2)) THEN
          ! No calculation can be performed!  Set vorticity to zero at this point.
          yvort(i,j,k) = 0.0
        ELSE IF(k1 < 2) THEN       ! Do a right-sided difference
          dwdx=(wr-ws(i,j,k))/(x(i+1)-x(i))
        ELSE IF(k2 < 2) THEN       ! Do a left-sided difference
          dwdx=(ws(i,j,k)-wl)/(x(i)-x(i-1))
        ELSE                  ! Do a normal centered difference
          dwdx=(wr-wl)/(x(i+1)-x(i-1))
        END IF

        if( k == 2 ) then
          kl = 2
          factor = 2.0
        else
          kl = k-1
          factor = 1.0
        END IF
        ! Calculate dudz

        dudz=0.25*factor*(u(i,j,k+1)+u(i+1,j,k+1)-u(i,j,kl)-u(i+1,j,kl))/(zs(i,j,k+1)-zs(i,j,k))

        ! Calculate xvort

        yvort(i,j,k)=dudz-dwdx
      END DO
    END DO
  END DO

  DO  k=2,nz-2
    DO j=1,ny-1
      yvort(   1,j,k)=yvort(   2,j,k)
      yvort(nx-1,j,k)=yvort(nx-2,j,k)
    END DO
  END DO

  DO j=1,ny-1
    DO i=1,nx-1
      yvort(i,j,   1)=yvort(i,j,   2)
      yvort(i,j,nz-1)=yvort(i,j,nz-2)
    END DO
  END DO


  RETURN
END SUBROUTINE cal_yvort
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINES cal_zvort                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE cal_zvort(u,v,x,y,nx,ny,nz,zvort,zp)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate Vort*10^5 (1/s) value.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  3/2/98
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: zvort(nx,ny,nz),zp(nx,ny,nz)
  REAL :: u(nx,ny,nz), v(nx,ny,nz)
  REAL :: x(nx), y(ny)
  INTEGER :: kw, ke, ks, kn
  REAL :: zs(nx,ny,nz)
  REAL :: us(nx,ny,nz),vs(nx,ny,nz)
  REAL :: viw,vie,uin,uis,zweightw,zweighte,zweightn,zweights
  REAL :: dudy,dvdx

  INTEGER :: i,j,k
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  ! Precalculate z, u, and v, at the scalar points

  DO j=2,ny-2
    DO i=2,nx-2
      DO k=2,nz-2
        zs(i,j,k)=0.5*(zp(i,j,k)+zp(i,j,k+1))
        us(i,j,k)=0.5*(u(i+1,j,k)+u(i,j,k))
        vs(i,j,k)=0.5*(v(i,j+1,k)+v(i,j,k))
      END DO
    END DO
  END DO

  DO j=2,ny-2
    DO i=2,nx-2
      DO k=2,nz-2
        kw=k
        ke=k
        ks=k
        kn=k

        ! Determine k levels for vertical interpolation of u and v on either side of center point
        IF(zs(i,j,k) < zp(i-1,j,k)) THEN
          DO WHILE (zs(i,j,k) < zp(i-1,j,kw))
            kw=kw-1
          END DO
        ELSE IF(zs(i,j,k) > zp(i-1,j,k+1)) THEN
          DO WHILE (zs(i,j,k) > zp(i-1,j,kw+1))
            kw=kw+1
          END DO
        END IF

        IF(zs(i,j,k) < zp(i+1,j,k)) THEN
          DO WHILE (zs(i,j,k) < zp(i+1,j,ke))
            ke=ke-1
          END DO
        ELSE IF(zs(i,j,k) > zp(i+1,j,k+1)) THEN
          DO WHILE (zs(i,j,k) > zp(i+1,j,ke+1))
            ke=ke+1
          END DO
        END IF

        IF(zs(i,j,k) < zp(i,j-1,k)) THEN
          DO WHILE (zs(i,j,k) < zp(i,j-1,ks))
            ks=ks-1
          END DO
        ELSE IF(zs(i,j,k) > zp(i,j-1,k+1)) THEN
          DO WHILE (zs(i,j,k) > zp(i,j-1,ks+1))
            ks=ks+1
          END DO
        END IF

        IF(zs(i,j,k) < zp(i,j+1,k)) THEN
          DO WHILE (zs(i,j,k) < zp(i,j+1,kn))
            kn=kn-1
          END DO
        ELSE IF(zs(i,j,k) > zp(i,j+1,k+1)) THEN
          DO WHILE (zs(i,j,k) > zp(i,j+1,kn+1))
            kn=kn+1
          END DO
        END IF

        ! Precalculate vertically-interpolated u and v values

        IF(kw >= 2 .and. kw <= nz-2) THEN
          zweightw=(zs(i,j,k)-zp(i-1,j,kw))/(zp(i-1,j,kw+1)-zp(i-1,j,kw))
          viw=(1.0-zweightw)*vs(i-1,j,kw)+zweightw*vs(i-1,j,kw+1)
        END IF

        IF(ke >= 2 .and. ke <= nz-2) THEN
          zweighte=(zs(i,j,k)-zp(i+1,j,ke))/(zp(i+1,j,ke+1)-zp(i+1,j,ke))
          vie=(1.0-zweighte)*vs(i+1,j,ke)+zweighte*vs(i,j,ke+1)
        END IF

        IF(ks >= 2 .and. ks <= nz-2) THEN
          zweights=(zs(i,j,k)-zp(i,j-1,ks))/(zp(i,j-1,ks+1)-zp(i,j-1,ks))
          uis=(1.0-zweights)*us(i,j-1,ks)+zweights*us(i,j-1,ks+1)
        END IF

        IF(kn >= 2 .and. kn <= nz-2) THEN
          zweightn=(zs(i,j,k)-zp(i,j+1,kn))/(zp(i,j+1,kn+1)-zp(i,j+1,kn))
          uin=(1.0-zweightn)*us(i,j+1,kn)+zweightn*us(i,j+1,kn+1)
        END IF

        ! Make sure the k levels are not below ground or above the model top and calculate dudy and dvdx
        IF((ks < 2 .and. kn < 2) .or. (ks > nz-2 .and. kn > nz-2)) THEN
          dudy = 0.0
        ELSE IF(ks < 2) THEN       ! Do a right-sided difference
          dudy=(uin-us(i,j,k))/(y(j+1)-y(j))
        ELSE IF(kn < 2) THEN       ! Do a left-sided difference
          dudy=(us(i,j,k)-uis)/(y(j)-y(j-1))
        ELSE                  ! Do a normal centered difference
          dudy=(uin-uis)/(y(j+1)-y(j-1))
        END IF

        IF((kw < 2 .and. ke < 2) .or. (kw > nz-2 .and. ke > nz-2)) THEN
          dvdx = 0.0
        ELSE IF(kw < 2) THEN       ! Do a right-sided difference
          dvdx=(vie-vs(i,j,k))/(x(i+1)-x(i))
        ELSE IF(ke < 2) THEN       ! Do a left-sided difference
          dvdx=(vs(i,j,k)-viw)/(x(i)-x(i-1))
        ELSE                  ! Do a normal centered difference
          dvdx=(vie-viw)/(x(i+1)-x(i-1))
        END IF

        ! Finally, calculate the vorticity
        zvort(i,j,k)=dvdx-dudy

     END DO
   END DO
  END DO

  DO j=2,ny-2
    DO i=2,nx-2
      zvort(i,j,   1)=zvort(i,j,   2)
      zvort(i,j,nz-1)=zvort(i,j,nz-2)
    END DO
  END DO

  DO k=1,nz-1
    DO j=2,ny-2
      zvort(   1,j,k)=zvort(   2,j,k)
      zvort(nx-1,j,k)=zvort(nx-2,j,k)
    END DO
  END DO

  DO  k=1,nz-1
    DO i=1,nx-1
      zvort(i,   1,k)=zvort(i,   2,k)
      zvort(i,ny-1,k)=zvort(i,ny-2,k)
    END DO
  END DO

  RETURN
END SUBROUTINE cal_zvort
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINES cal_xvort_flat              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
! Only valid for Cartesian grid, flat terrain case
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate Vort*10^5 (1/s) value.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  3/2/98
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
SUBROUTINE cal_xvort_flat(v,w,y,zp,nx,ny,nz, xvort)

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: xvort(nx,ny,nz)
  REAL :: v(nx,ny,nz), w(nx,ny,nz)
  REAL :: y(ny),zp(nx,ny,nz)

  INTEGER :: i,j,k,kl
  REAL :: factor
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  !DO k=2,nz-2
  DO k=1,nz-2

    !IF( k == 2 ) THEN
    IF( k <= 2 ) THEN
      kl = k
      factor = 2.0
    ELSE
      kl = k-1
      factor = 1.0
    END IF

    DO  j=2,ny-2
      DO  i=1,nx-1
        xvort(i,j,k)=0.25*(                                             &
         (w(i,j+1,k)+w(i,j+1,k+1)-w(i,j-1,k)-w(i,j-1,k+1))/(y(j+1)-y(j))&
 -factor*(v(i,j,k+1)+v(i,j+1,k+1)-v(i,j,kl )-v(i,j+1,kl ))/(zp(i,j,k+1)-zp(i,j,k)) )
      END DO
    END DO
  END DO

  DO  k=2,nz-2
    DO i=1,nx-1
      xvort(i,   1,k)=xvort(i,   2,k)
      xvort(i,ny-1,k)=xvort(i,ny-2,k)
    END DO
  END DO


  DO j=1,ny-1
    DO i=1,nx-1
      !xvort(i,j,   1)=xvort(i,j,   2)
      xvort(i,j,nz-1)=xvort(i,j,nz-2)
    END DO
  END DO

  RETURN
END SUBROUTINE cal_xvort_flat

SUBROUTINE cal_yvort_flat(u,w,x,zp,nx,ny,nz, yvort)

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: yvort(nx,ny,nz)
  REAL :: u(nx,ny,nz), w(nx,ny,nz)
  REAL :: x(ny),zp(nx,ny,nz)

  INTEGER :: i,j,k, kl
  REAL :: factor
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  !DO k=2,nz-2
  DO k=1,nz-2
    !IF( k == 2 ) THEN
    IF( k <= 2 ) THEN
      kl = k
      factor = 2.0
    ELSE
      kl = k-1
      factor = 1.0
    END IF
    DO  j=1,ny-1
      DO  i=2,nx-2
        yvort(i,j,k)=0.25*(                                                &
        -(w(i+1,j,k)+w(i+1,j,k+1)-w(i-1,j,k)-w(i-1,j,k+1))/(x(i+1)-x(i))   &
 +factor*(u(i,j,k+1)+u(i+1,j,k+1)-u(i,j,kl )-u(i+1,j,kl ))   &
        /(zp(i,j,k+1)-zp(i,j,k)) )
      END DO
    END DO
  END DO

  DO  k=2,nz-2
    DO j=1,ny-1
      yvort(   1,j,k)=yvort(   2,j,k)
      yvort(nx-1,j,k)=yvort(nx-2,j,k)
    END DO
  END DO

  DO j=1,ny-1
    DO i=1,nx-1
      !yvort(i,j,   1)=yvort(i,j,   2)
      yvort(i,j,nz-1)=yvort(i,j,nz-2)
    END DO
  END DO

  RETURN
END SUBROUTINE cal_yvort_flat

SUBROUTINE cal_zvort_flat(u,v,x,y,nx,ny,nz,zvort)

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: zvort(nx,ny,nz)
  REAL :: u(nx,ny,nz), v(nx,ny,nz)
  REAL :: x(nx), y(ny)

  INTEGER :: i,j,k
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  DO k=1,nz-2
    DO  j=2,ny-2
      DO  i=2,nx-2
        zvort(i,j,k)= (                                            &
               (v(i+1,j,k)-v(i-1,j,k)+v(i+1,j+1,k)-v(i-1,j+1,k))/       &
               (4*(x(i+1)-x(i)))-                                       &
               (u(i,j+1,k)-u(i,j-1,k)+u(i+1,j+1,k)-u(i+1,j-1,k))/       &
               (4*(y(j+1)-y(j))) )
      END DO
    END DO
  END DO

  DO j=2,ny-2
    DO i=2,nx-2
      !zvort(i,j,   1)=zvort(i,j,   2)
      zvort(i,j,nz-1)=zvort(i,j,nz-2)
    END DO
  END DO

  DO k=1,nz-1
    DO j=2,ny-2
      zvort(   1,j,k)=zvort(   2,j,k)
      zvort(nx-1,j,k)=zvort(nx-2,j,k)
    END DO
  END DO

  DO k=1,nz-1
    DO i=1,nx-1
      zvort(i,   1,k)=zvort(i,   2,k)
      zvort(i,ny-1,k)=zvort(i,ny-2,k)
    END DO
  END DO

  RETURN
END SUBROUTINE cal_zvort_flat

!
!##################################################################
!##################################################################
!######                                                      ######
!######           SUBROUTINES cal_xvortmix_flat              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
! NEW SUBROUTINES FOR MIXING TERMS
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate mixing vorticity x-component value.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  3/2/98
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
! Only valid for Cartesian grid, flat terrain case
!
SUBROUTINE cal_xvortmix_flat(v,w,y,zp,nx,ny,nz, xvort)

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: xvort(nx,ny,nz)
  REAL :: v(nx,ny,nz), w(nx,ny,nz)
  REAL :: y(ny),zp(nx,ny,nz)

  INTEGER :: i,j,k,kl
  REAL :: factor
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  DO k=2,nz-2

    IF( k == 2 ) then
      kl = 2
      factor = 2.0
    ELSE
      kl = k-1
      factor = 1.0
    END IF

    DO  j=2,ny-2
      DO  i=1,nx-1
        xvort(i,j,k)=0.25*(                                             &
         (w(i,j+1,k)+w(i,j+1,k+1)-w(i,j-1,k)-w(i,j-1,k+1))/(y(j+1)-y(j))&
 -factor*(v(i,j,k+1)+v(i,j+1,k+1)-v(i,j,kl )-v(i,j+1,kl ))/(zp(i,j,k+1)-zp(i,j,k)) )
      END DO
    END DO
  END DO

  DO  k=2,nz-2
    DO i=1,nx-1
      xvort(i,   1,k)=xvort(i,   2,k)
      xvort(i,ny-1,k)=xvort(i,ny-2,k)
    END DO
  END DO

  DO j=1,ny-1
    DO i=1,nx-1
      xvort(i,j,   1)=xvort(i,j,   2)
      xvort(i,j,nz-1)=xvort(i,j,nz-2)
    END DO
  END DO

  RETURN
END SUBROUTINE cal_xvortmix_flat

SUBROUTINE cal_yvortmix_flat(u,w,x,zp,nx,ny,nz, yvort)

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: yvort(nx,ny,nz)
  REAL :: u(nx,ny,nz), w(nx,ny,nz)
  REAL :: x(ny),zp(nx,ny,nz)

  INTEGER :: i,j,k, kl
  REAL :: factor
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  DO k=2,nz-2
    IF( k == 2 ) then
      kl = 2
      factor = 2.0
    ELSE
      kl = k-1
      factor = 1.0
    END IF
    DO  j=1,ny-1
      DO  i=2,nx-2
        yvort(i,j,k)=0.25*(                                             &
        -(w(i+1,j,k)+w(i+1,j,k+1)-w(i-1,j,k)-w(i-1,j,k+1))/(x(i+1)-x(i))&
 +factor*(u(i,j,k+1)+u(i+1,j,k+1)-u(i,j,kl )-u(i+1,j,kl ))   &
        /(zp(i,j,k+1)-zp(i,j,k)) )
      END DO
    END DO
  END DO

  DO  k=2,nz-2
    DO j=1,ny-1
      yvort(   1,j,k)=yvort(   2,j,k)
      yvort(nx-1,j,k)=yvort(nx-2,j,k)
    END DO
  END DO

  DO j=1,ny-1
    DO i=1,nx-1
      yvort(i,j,   1)=yvort(i,j,   2)
      yvort(i,j,nz-1)=yvort(i,j,nz-2)
    END DO
  END DO

  RETURN
END SUBROUTINE cal_yvortmix_flat

SUBROUTINE cal_zvortmix_flat(u,v,x,y,nx,ny,nz,zvort)

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: zvort(nx,ny,nz)
  REAL :: u(nx,ny,nz), v(nx,ny,nz)
  REAL :: x(nx), y(ny)

  INTEGER :: i,j,k
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  DO k=2,nz-2
    DO  j=2,ny-2
      DO  i=2,nx-2
        zvort(i,j,k)= (                                            &
               (v(i+1,j,k)-v(i-1,j,k)+v(i+1,j+1,k)-v(i-1,j+1,k))/       &
               (4*(x(i+1)-x(i)))-                                       &
               (u(i,j+1,k)-u(i,j-1,k)+u(i+1,j+1,k)-u(i+1,j-1,k))/       &
               (4*(y(j+1)-y(j))) )
      END DO
    END DO
  END DO

  DO j=2,ny-2
    DO i=2,nx-2
      zvort(i,j,   1)=zvort(i,j,   2)
      zvort(i,j,nz-1)=zvort(i,j,nz-2)
    END DO
  END DO

  DO k=1,nz-1
    DO j=2,ny-2
      zvort(   1,j,k)=zvort(   2,j,k)
      zvort(nx-1,j,k)=zvort(nx-2,j,k)
    END DO
  END DO

  DO  k=1,nz-1
    DO i=1,nx-1
      zvort(i,   1,k)=zvort(i,   2,k)
      zvort(i,ny-1,k)=zvort(i,ny-2,k)
    END DO
  END DO

  RETURN
END SUBROUTINE cal_zvortmix_flat
