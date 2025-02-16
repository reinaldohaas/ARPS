!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE BARNES_R5                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE barnes_r5 (nx,ny,nz,tor,t,wt_p,cf_modelfg,                   &
           l_stn_clouds,default_clear_cover,cld_snd,wt_snd,             &
           spcng,i_snd,j_snd,n_cld_snd,                                 &
           istatus)
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Horizontally interpolate the SAO cloud soundings to ADAS grid.
!  The interpolation scheme is Barnes-type, the weight given to a
!  data point is:
!
!  wt(m) = [(100.0*bias)/(1.533*dgrd(m)**2 + 1.0)]**(5.0/2.0)
!  wt(m) ~ 34367.237 /(dgrd(m)**5)
!
!  where
!
!    m  is the index of a data point
!    bias=1.0        ! it can be a function of data location
!    dgrd(m) is the distance from the mth data point to
!          the grid point, in terms of the number of grid spacings.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:  (Jian Zhang)
!  07/95   Based on the LAPS cloud analysis code
!
!  MODIFICATION HISTORY:
!  02/06/96  J. Zhang
!            Modified for ADAS. Added documents.
!  05/10/96  J. Zhang
!            Added the interpolation for the station grid column
!            (skipped during the first step).
!  03/14/97  J. Zhang
!            Clean up the code and implemented for the offcial
!            arps4.2.4 version
!  07/30/97  J. Zhang
!            Added minimum weight threshold for Barnes interpolation.
!            This can avoid the interpolations at the grid points
!            where there is no station nearby.
!  08/06/97  J. Zhang
!            Change adascld24.inc to adascld25.inc.
!            Set weight_modelfg=1.0 for the cases when there is
!            no SAO cloud soundings.
!  09/09/97  J. Zhang
!            Using the different influence cutoff radius for
!            cloud soundingd with different cloud coverage
!  09/10/97  J. Zhang
!            Awoke the procedure which uses background RH field to
!            derive cloud cover field when there is no SAO cloud
!            sounding.
!            Change adascld25.inc to adascld26.inc.
!  11/18/97  J. Zhang
!            Moved "max_obs" to adascld26.inc.
!  11/19/97  J. Zhang
!            Added test for the case when ncnt > max_obs.
!  11/05/07  K.W. Thomas
!            Set "iskip=1" (no skipping) as the default, as it really
!            really doesn't take more than a few seconds.  Allow users
!            to undo this if they have a really old and slow system.
!            Don't enter the bilinear interpolation code when there
!            isn't interpolation to be done (iskip=1).
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!
!  OUTPUT:
!
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INCLUDE 'adas.inc'
  INCLUDE 'bndry.inc'
  INCLUDE 'mp.inc'
!
!   INPUT from adas.inc:
!   real r_missing      ! missing value flag
!   real range_cld      ! cutoff radius for interpolation
!   integer max_cld_snd ! max. possible # of cloudy soundings
!   integer max_obs     ! max. possible # of cldy data points in 3D
                        ! domain (should be about max_cld_snd*nz)
!   integer i_perimeter ! an extension of the ADAS domain for
                        ! searching SAOs.
!
! INPUT:
!
  INTEGER :: nx,ny,nz           ! grid size
  INTEGER :: n_cld_snd          ! actual # of cld snd

  REAL :: tor(nx,ny,nz)         ! 1st guess of gridded cloud cover
  REAL :: cf_modelfg(nx,ny,nz)  ! 1st guess for cloud cover field
  REAL :: wt_p(nx,ny,nz)        ! weights for first guess field

  LOGICAL :: l_stn_clouds       ! using SAO stns'cld sndings?
  REAL :: default_clear_cover   ! defaule cld value for clear sky

  REAL :: cld_snd(max_cld_snd,nz) ! obs. cloud cover sounding
  REAL :: wt_snd(max_cld_snd,nz)  ! weights for obs. cld sounding

  INTEGER :: i_snd(max_cld_snd)   ! i-loc of cloud sounding stn
  INTEGER :: j_snd(max_cld_snd)   ! j-loc of cloud sounding stn
!
! OUTPUT:
!
  REAL :: t(nx,ny,nz)       ! the gridded analysis of cld cvr
  REAL :: tem1(nx,ny,nz)    ! temporary array
  INTEGER :: istatus
!
! LOCAL:
!
  INTEGER :: nx_lut,ny_lut      ! dims for look up tbl of intp wgt.
  INTEGER :: ix_low,ix_high,iy_low,iy_high
  INTEGER :: ix_low_lg,ix_high_lg,iy_low_lg,iy_high_lg

  INTEGER :: n_fnorm
!
  INTEGER :: iskip        ! # of grid pts skipped when performing
                         ! Barnes intp. (for time-saving)
  INTEGER :: lowi_lut(nx) ! index for skipped grid points.
  INTEGER :: lowj_lut(ny) ! index for skipped grid points.
!
  REAL, allocatable :: fnorm(:)    ! normalized weights
!
  INTEGER :: iob(max_obs)   ! i-loc of each cld data pt. in ADAS grid
  INTEGER :: job(max_obs)   ! j-loc of each cld data pt. in ADAS grid
  INTEGER :: kob(max_obs)   ! k-lvl of each cld data pt. in ADAS grid
  INTEGER :: nob(max_obs)   ! sequential index of each cldy data pt.

  REAL, allocatable :: iiilut(:,:)
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: nlast (nz)        ! the index for the cloud sounding data
                             ! point at each level.
  LOGICAL :: l_analyze (nz)
!
  REAL :: weight_modelfg       ! the weight given to the 1st guess
                             ! cloud cover field
  REAL :: rden_ratio,radm2,rfac
  PARAMETER ( rfac=10.0, radm2=1.533/rfac**2)
!   parameter ( radm2=1.533)
!
  REAL :: exp_dist_wt          ! power of distance in weight
  PARAMETER( exp_dist_wt = 5.0)

  REAL :: spcng                ! grid spacing in m

  REAL :: rr,rr_max
  INTEGER :: iii,iii_max

  INTEGER :: iiizero,ncnt
  REAL :: bias_iii             ! a bias may be a func of locatns

  LOGICAL :: cldprt
  INTEGER :: i,j,k,n,nobs,nstart,nstop,ii,jj,nn,n1,jmid,ibeg
  REAL :: weight,sum,sumwt

  INTEGER :: il,ih,jl,jh,km1
  REAL :: z1,z2,z3,z4,fraci,fracj
  INTEGER :: astat
  INTEGER :: nxlg, nylg
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  nx_lut = nx + i_perimeter - 1
  ny_lut = ny + i_perimeter - 1
  ix_low = 1  - i_perimeter
  ix_high= nx + i_perimeter
  iy_low = 1  - i_perimeter
  iy_high= ny + i_perimeter

  IF (mp_opt > 0 ) THEN

     nxlg = (nx-3) * (nproc_x) + 3
     nylg = (ny-3) * (nproc_y) + 3
     nx_lut = nxlg + i_perimeter - 1
     ny_lut = nylg + i_perimeter - 1
!
!-----------------------------------------------------------------------
!
!  The MPI case is a bit trickier.  All processors need to have the same
!  list of stations that they do computations before, and match the
!  non-MPI case.  The catch is that "ix/iy_low/high" above is on the
!  small grid, we must compute "ix/iy_low/high_lg" on the large grid,
!  but relative to the small grid values of the current processor.
!
!-----------------------------------------------------------------------
!
      ix_low_lg = ix_low - (loc_x-1) * (nx-3)
      ix_high_lg = ix_high + (nproc_x-loc_x) * (nx - 3)
      iy_low_lg = iy_low - (loc_y-1) * (ny-3)
      iy_high_lg = iy_high + (nproc_y-loc_y) * (ny - 3)

  ENDIF

  n_fnorm = 1.6 * ((nx_lut*nx_lut)+(ny_lut*ny_lut))

  allocate (fnorm(n_fnorm),stat=astat)
  IF (astat /= 0) THEN
    WRITE (6,'(a)') "BARNES_R5: ERROR allocating fnorm, exiting"
    STOP
  END IF
  allocate (iiilut(-nx_lut:nx_lut,-ny_lut:ny_lut),stat=astat)
  IF (astat /= 0) THEN
    WRITE (6,'(a)') "BARNES_R5: ERROR allocating fnorm, exiting"
    STOP
  END IF

  cldprt=(clddiag == 1)
  IF (myproc /= 0) cldprt = .FALSE.
  jmid=ny/2
  ibeg=MAX(1,nx-10)
  IF (myproc == 0) WRITE(6,'(1x,a)') ' Barnes_r5 called'

  istatus = 0
  IF(cldprt) THEN
    WRITE(6,'(a,i8)') ' # of cloud soundings = ',n_cld_snd
    WRITE(6,'(/a)') ' ==== BARNES cloud cover first guess===='
    WRITE(6,'(/1X,3X,11(1X,i4,1X))') (i,i=nx-10,nx)
    WRITE(6,'(1X,i3,11F6.1)')                                           &
         (k,(cf_modelfg(i,jmid,k),i=ibeg,nx),k=nz,1,-1)
  END IF
!
!-----------------------------------------------------------------------
!
!  Obtain and/or iterate for value of iskip
!
!-----------------------------------------------------------------------
!
  iskip = nint(20000. / spcng)           ! # of grid pts skipped
                                         ! during calcul'n
  iskip = MAX(iskip,1)

!
!-----------------------------------------------------------------------
!
!  The purpose of "iskip" is to make the computations a bit faster.  This
!  may at one time have led to faster runs.  Limit tested shows that
!  processing all the data (iskip=1) adds only a small amount of time.
!  In addition, on complex problems, MPI runs are better.
!
!  In the interests of producing the same answers for MPI and non-MPI
!  runs, "iskip" is now hardwired to be 1.  On slower system, commenting
!  out the next line will make the run a bit faster, so it can be
!  uncommented out if necessary.
!
!-----------------------------------------------------------------------
!

  iskip = 1

  100   rden_ratio = ((FLOAT(nx)-1.)/FLOAT(iskip))

  IF( ABS(rden_ratio - nint(rden_ratio)) > 0.001) THEN
    IF (myproc == 0) WRITE(6,*)' Bad value of iskip'
    IF(iskip > 1)THEN
      iskip = iskip - 1
      GO TO 100
    ELSE IF( iskip == 1) THEN
      IF (myproc == 0) WRITE(6,*)' Code error - stop'
      CALL arpsstop("problem with 'iskip'",1)
    END IF
  ELSE
    IF (myproc == 0) WRITE(6,*)' Good value of iskip = ',iskip
  END IF
!
!-----------------------------------------------------------------------
!
!  Create an array (lookup table) of weights as function of the
!  distance iii (in the unit of grid spacings)
!
!-----------------------------------------------------------------------
!
  IF( l_stn_clouds .AND. myproc == 0 )                                 &
    PRINT*, ' # of cloud soundings = ',n_cld_snd

  iiizero = 0

  DO iii = 1,n_fnorm
    bias_iii = 1E0    ! It can be func of location iii later

    fnorm(iii) = (100.*bias_iii/FLOAT(iii)) ** (exp_dist_wt/2.0)
    IF( fnorm(iii) == 0. .AND. iiizero == 0) THEN
      iiizero = iii
      IF (myproc == 0)                                                 &
        WRITE(6,*)' WARNING: fnorm array reached zero, iii=',iiizero
    END IF
  END DO  ! III

  IF (myproc == 0) THEN
    WRITE(6,*)' Min possible weight = ',fnorm(n_fnorm)
    WRITE(6,*)' Max possible weight = ',fnorm(1)
  END IF

  IF (n_cld_snd <= 0) THEN
    weight_modelfg = 1.0
  ELSE
    weight_modelfg = 0.0
  END IF
!
!-----------------------------------------------------------------------
!
!  Count the number of cloudy data points on each cloud height
!  level, and mark each cloudy data point with iob- and job-index.
!
!-----------------------------------------------------------------------
!
  ncnt=0        ! # of cloudy data points on each cloud grid level

  DO k = 1, nz

    IF (.NOT. l_stn_clouds) THEN   ! using gridded non-model 1st
                                   ! guess field
      nlast(k) = ncnt
      DO j=1,ny
        DO i=1,nx
          IF( tor(i,j,k) == r_missing) GO TO 223
          ncnt=ncnt+1
          iob(ncnt)=i
          job(ncnt)=j
          kob(ncnt)=k
          nlast(k) = ncnt
          223         CONTINUE
        END DO ! i
      END DO ! j

    ELSE ! l_stn_clouds = .TRUE.,  use cloud soundings
      nlast(k) = ncnt

      IF (n_cld_snd >= 1) THEN
        DO n = 1,n_cld_snd
          IF(cld_snd(n,k) < default_clear_cover) GO TO 233
!
!-----------------------------------------------------------------------
!
!  Test if out of the bounds of the ADAS_plus domain
!
!-----------------------------------------------------------------------
!
          IF (mp_opt > 0) THEN
            IF(i_snd(n) < ix_low_lg .OR. i_snd(n) > ix_high_lg) GO TO 233
            IF(j_snd(n) < iy_low_lg .OR. j_snd(n) > iy_high_lg) GO TO 233
          ELSE
            IF(i_snd(n) < ix_low .OR. i_snd(n) > ix_high) GO TO 233
            IF(j_snd(n) < iy_low .OR. j_snd(n) > iy_high) GO TO 233
          END IF

          ncnt=ncnt+1
          IF (ncnt > max_obs) THEN
            IF (myproc == 0) THEN
              PRINT*,'# of cldy grid pts exceeds the max. # allowed'
              PRINT*,'# of cldy grid pts:',ncnt,' # allowed:',max_obs
              PRINT*,'Plz increase MAX_OBS in adas.inc.'
              PRINT*,' Aborting...'
            END IF
            CALL arpsstop("MAX_OBS too small in 'adas.inc'",1)
          END IF
          iob(ncnt)=i_snd(n)
          job(ncnt)=j_snd(n)
          kob(ncnt)=k
          nob(ncnt)=n
          nlast(k) = ncnt
          233           CONTINUE
        END DO ! n
      END IF

    END IF ! l_stn_clouds

    IF( k > 1) THEN
      km1 = k - 1
      l_analyze(k) = .false.

      IF(.NOT. l_stn_clouds) THEN
        DO j=1,ny
          DO i=1,nx
            IF( tor(i,j,k) /= tor(i,j,km1)) THEN
              l_analyze(k) = .true.
              GO TO 250
            END IF
          END DO ! i
        END DO ! j

      ELSE
        DO n=1,n_cld_snd
          IF( cld_snd(n,k) /= cld_snd(n,km1)) THEN
            l_analyze(k) = .true.
            goto 250
          END IF
        END DO ! n

      END IF ! l_stn_clouds

      250       CONTINUE

    ELSE ! k .eq. 1
      l_analyze(k) = .true.
    END IF

  END DO ! K

  IF(ncnt == 0) THEN
    IF (myproc == 0) THEN
      PRINT*, 'No SAO cloud sounding data for Barnes'
      PRINT*, 'Using the background cloud cover field'
    END IF

    DO k=2,nz-2
      DO j=1,ny
        DO i=1,nx
          t(i,j,k) = cf_modelfg(i,j,k)
        END DO
      END DO
    END DO
    istatus = 1
    RETURN
  ELSE
    IF (myproc == 0) WRITE(6,*)' Ncnt = ',ncnt
  END IF
!
!-----------------------------------------------------------------------
!
!  Create a lookup table for fnorm(iii): the weights as a function
!  of i- and j- distances.
!
!-----------------------------------------------------------------------
!
  rr_max = (nx_lut-1)**2 + (ny_lut-1)**2
  iii_max = radm2*100.*rr_max + 1.
  IF (myproc == 0)                                                      &
    PRINT*,' radm2*100,iii_max,n_fnorm: ',radm2*100.,iii_max,n_fnorm

  IF( iii_max > n_fnorm) THEN
    IF (myproc == 0)                                                    &
      PRINT*,' iii_max is too large, increase n_fnorm',iii_max,n_fnorm
    istatus = 0
    RETURN
  END IF
!
  DO i = -nx_lut,nx_lut
    DO j = -ny_lut,ny_lut
      rr=i*i+j*j
      iii=radm2*100.*rr+1.
      IF( iii > n_fnorm) iii=n_fnorm
      iiilut(i,j) = fnorm(iii)            ! lookup table for weights
                                          ! fnorm as a function of
                                          ! distance iii.
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Loop through each cloud height level.
!
!-----------------------------------------------------------------------
!
  DO k = 1, nz

    IF(k == 1)THEN
      nstart = 1
    ELSE
      nstart = nlast(k-1) + 1
    END IF

    nstop = nlast(k)
    nobs = nstop-nstart+1       ! # of obs. in one level


    IF (nobs <= 0) THEN  ! No Obs; Set level to model 1st guess
!
      IF(cldprt) WRITE(6,'(a,4i5,a)') ' lvl,nstart,nstop,nobs=',        &
           k,nstart,nstop,nobs,' No Obs; Set level to model fg'
      DO j=1,ny
        DO i=1,nx
          t(i,j,k) = cf_modelfg(i,j,k)
        END DO ! i
      END DO ! j

    ELSE IF( weight_modelfg > 1.0E-15 .OR. l_analyze(k) ) THEN

      IF(cldprt) WRITE(6,'(a,4i5)') ' lvl,nstart,nstop,nobs=',          &
        k,nstart,nstop,nobs
!
!-----------------------------------------------------------------------
!
!    Analyze every iskip-th  grid point
!
!-----------------------------------------------------------------------
!
      DO j=1,ny,iskip
        DO i=1,nx,iskip
          sum=0.
          sumwt=0.
          IF(.NOT. l_stn_clouds)THEN
            DO n=nstart,nstop
              ii=iob(n)
              jj=job(n)
              iii = (i-ii)*(i-ii) + (j-jj)*(j-jj)
              rr = FLOAT (iii)
              rr = spcng*SQRT (rr)
              IF (rr <= range_cld) THEN
!
!-----------------------------------------------------------------------
!
!c              gridded obs are being weighted
!
!-----------------------------------------------------------------------
!
                weight = iiilut(i-ii,j-jj) * wt_p(ii,jj,k)
                sum=weight*tor(ii,jj,k)+sum
                sumwt=sumwt+weight
              END IF
            END DO  ! n
          ELSE
            DO n=nstart,nstop
              ii=iob(n)
              jj=job(n)
              nn=nob(n)
              iii = (i-ii)*(i-ii) + (j-jj)*(j-jj)
              rr = FLOAT (iii)
              rr = spcng*SQRT (rr)
              IF (rr <= range_cld) THEN
!
!-----------------------------------------------------------------------
!
!c              Station obs are being weighted
!
!-----------------------------------------------------------------------
!
                weight = iiilut(i-ii,j-jj) * wt_snd(nn,k)
                sum=weight*cld_snd(nn,k)+sum
                sumwt=sumwt+weight
              END IF
            END DO ! n
          END IF

          sum   = sum   + weight_modelfg * cf_modelfg(i,j,k)
          sumwt = sumwt + weight_modelfg

          IF (ABS(sumwt) <= 1.0E-15) THEN
!9/3/97        t(i,j,k) = r_missing
!              istatus = 0
!              print*,'(00): ',i,j,k,' tot sum:',sum,'  tot sumwt:'
!     :                ,sumwt,'  t:',t(i,j,k)
            t(i,j,k) = 0.01
          ELSE
            t(i,j,k)=sum/sumwt
          END IF

        END DO ! i
      END DO ! j
!
!-----------------------------------------------------------------------
!
!      Bilinearly interpolate to fill in rest of domain
!
!-----------------------------------------------------------------------
!
      IF (iskip > 1) THEN
        DO i = 1,nx
          lowi_lut(i) = (i-1)/iskip*iskip + 1
          IF( i == nx) lowi_lut(i) = lowi_lut(i) - iskip
        END DO ! i

        DO j = 1,ny
          lowj_lut(j) = (j-1)/iskip*iskip + 1
          IF( j == ny) lowj_lut(j) = lowj_lut(j) - iskip
        END DO ! i

        DO j=1,ny
          jl = lowj_lut(j)
          jh = jl + iskip
          IF (jh > ny) jh = ny
          fracj = FLOAT(j-jl)/FLOAT(iskip)

          DO i=1,nx
            il = lowi_lut(i)
            ih = il + iskip
            IF (ih > nx) ih = nx
            fraci = FLOAT(i-il)/FLOAT(iskip)

            z1=t(il,jl,k)
            z2=t(ih,jl,k)
            z3=t(ih,jh,k)
            z4=t(il,jh,k)

            t(i,j,k) =  z1+(z2-z1)*fraci+(z4-z1)*fracj                    &
                      - (z2+z4-z3-z1)*fraci*fracj
          END DO ! i
        END DO ! j
      END IF ! iskip > 1

      IF (l_stn_clouds .AND. n_cld_snd >= 1) THEN
!
!-----------------------------------------------------------------------
!
!  Reinterpolate at the station grid points which were skipped
!  previously.
!
!-----------------------------------------------------------------------
!
        DO n1=1,n_cld_snd
          i = i_snd(n1)
          j = j_snd(n1)
          IF(i < 1.OR.i > nx.OR.j < 1.OR.j > ny) GO TO 252
          IF( (i-1) /= (i-1)/iskip*iskip .OR. (j-1) /= (j-1)/iskip*iskip) THEN
            sum=0.
            sumwt=0.
            DO n=nstart,nstop
              ii=iob(n)
              jj=job(n)
              nn=nob(n)
              iii = (i-ii)*(i-ii) + (j-jj)*(j-jj)
              rr = FLOAT (iii)
              rr = spcng*SQRT (rr)
              IF (rr <= range_cld) THEN
!
!-----------------------------------------------------------------------
!
!c                Station obs are being weighted
!
!-----------------------------------------------------------------------
!
                weight = iiilut(i-ii,j-jj) * wt_snd(nn,k)
                sum=weight*cld_snd(nn,k)+sum
                sumwt=sumwt+weight
              END IF
            END DO ! n

            sum = sum + weight_modelfg * cf_modelfg(i,j,k)
            sumwt = sumwt + weight_modelfg

            IF (ABS(sumwt) <= 1.0E-15) THEN
!                  t(i,j,k) = r_missing
!                  istatus = 0
!                  print*,'(00*): ',i,j,k,' tot sum:',sum,
!     :                 '  tot sumwt:',sumwt,'  t:',t(i,j,k)
              t(i,j,k) = 0.01
            ELSE
              t(i,j,k)=sum/sumwt
            END IF

          END IF  ! Stn (i,j) was skipped during the 1st-pass intp.
          252           CONTINUE

        END DO ! n1
      END IF  ! l_stn_clouds = .true. and n_cld_snd >0

    ELSE       ! no 1st guess and obs. is identical to lower lvl
!
!-----------------------------------------------------------------------
!
!c        Obs are identical; Copy analysis from 1 level below
!
!-----------------------------------------------------------------------
!
      IF(cldprt) WRITE(6,'(a,4i5,a)') ' lvl,nstart,nstop,nobs=',        &
        k,nstart,nstop,nobs,' Identical Obs; Copy from 1 lvl down'

      km1 = k - 1

      DO j=1,ny
        DO i=1,nx
          t(i,j,k) = t(i,j,km1)
        END DO ! i
      END DO ! j

    END IF

  END DO ! K

  IF (mp_opt > 0) THEN
    CALL mpsendrecv2dew(t,nx,ny,nz,ebc,wbc,1,tem1)
    CALL mpsendrecv2dns(t,nx,ny,nz,nbc,sbc,2,tem1)
  END IF

  istatus = 1
!
  DEALLOCATE (fnorm,stat=astat)
  DEALLOCATE (iiilut,stat=astat)
  RETURN
END SUBROUTINE barnes_r5
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE HGT_TO_ZK                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE hgt_to_zk (z,zk,nk,zs_1d,istatus)
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  To find the k-position of a given height in the model grid
!  If the height is above the model highest level, the largest
!  level index (nk) will be returned. If the height is below
!  the lowest model level, the lowestlevel index (1) will be
!  returned.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (Jian Zhang)
!  05/01/96
!
!  MODIFICATION HISTORY:
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
!  INPUT:
  REAL :: z         ! the height (m) to be coverted to k-index
  INTEGER :: nk     ! number of vertical levels in model grid
  REAL :: zs_1d(nk) ! the physical height of each model level
!
!  OUTPUT:
  INTEGER :: istatus
  REAL :: zk        ! the k-position of the ht. in the model grid
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: k,nkm1
  REAL :: height_above,thickness,frac
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  zk = FLOAT(nk-1)
  height_above = zs_1d(nk-1)

  IF( z > height_above) THEN
    nkm1=nk-1
!        write(6,*)
!        write(6,601) z,zk,nkm1,zs_1d(nkm1)
!601     format(' Warning: above model domain in HGT_TO_ZK'/
!     :        1x,' z=',f9.0,' zk=',f6.2,' nk=',i2,'  zs(nk)=',f9.0)
!        write(6,*) 'Top level index will be assigned to ZK.'
    istatus=0
!
    GO TO 999
  END IF

  DO k = nk-2,1,-1

    IF( height_above >= z .AND. zs_1d(k) <= z) THEN
      thickness = height_above - zs_1d(k)
      frac = (z - zs_1d(k))/thickness
      zk = FLOAT(k) + frac
      istatus=1
      GO TO 999
    END IF

    height_above = zs_1d(k)
  END DO ! K

  zk = 1.0

!      write(6,*)
!      write(6,602) z,zk,zs_1d(1)
!602   format(' Warning: below model level_1 in HGT_TO_ZK'/
!     :      1x,' z=',f9.0,' zk=',f6.2,' k=1  zs(1)=',f9.0)
!      write(6,*) 'Bottom level index will be assigned to ZK.'
  istatus=0

  999   CONTINUE
  RETURN
END SUBROUTINE hgt_to_zk


SUBROUTINE wmix (p,t,dd,w,nl)
!
!            07-jun-84       source code from wang
!                            (not on u. of w. source tape)
!
!            dd      - dewpoint depression
  DIMENSION p(*), t(*), dd(*), w(*)
  DO i = 1, nl
    td = t(i) - dd(i)
    IF (td > 253.0) GO TO 110
    es = vpice(td)
    GO TO 120
    110     es = satvap(td)
    120     w(i) = 622.0 * es / p(i)
  END DO
  RETURN
END SUBROUTINE wmix
!
!##################################################################
!##################################################################
!######                                                      ######
!######              FUNCTION WSAT                           ######
!######                                                      ######
!##################################################################
!##################################################################
!

  FUNCTION wsat(p, t)
!
!         1 4-may-84 added from csu vas code
!                            (not on u. of w. source tape)
!
! jz 4/29/98 to avoid warning message when compile on origin
  DIMENSION p(1),t(1),dd(1),w(1)
  dd(1)  = 0.0
! jz 4/29/98 end
  CALL wmix(p,t,dd,w,1)
  wsat=w(1)
!      wsat=w
!
!         1 4-may-84 return added
!
  RETURN
  END FUNCTION wsat
!
!##################################################################
!##################################################################
!######                                                      ######
!######              FUNCTION VPICE                          ######
!######                                                      ######
!##################################################################
!##################################################################
!

  FUNCTION vpice(temp)
!
!            07-jun-84       source code from wong
!                            (not on u. of w. source tape)
  REAL(KIND=8) :: c(6)
  REAL(KIND=8) :: t
!
  DATA c/ 0.7859063157D+00,  0.3579242320D-01,                          &
       -0.1292820828D-03,  0.5937519208D-06,                            &
        0.4482949133D-09,  0.2176664827D-10/
!
  t = temp - 273.16
  vplog = c(1) + t * (c(2) + t * (c(3) + t * (c(4) + t * (c(5) +        &
               t *  c(6)))))
  vpice = 10.0 ** vplog
  RETURN
  END FUNCTION vpice
!
!##################################################################
!##################################################################
!######                                                      ######
!######              FUNCTION SATVAP                         ######
!######                                                      ######
!##################################################################
!##################################################################
!

  FUNCTION satvap(temp)
!
!            07-jun-84       source code from wang
!                            (not on u. of w. source tape)
!
  REAL(KIND=8) :: c(10)
  REAL(KIND=8) :: s, t
!
  DATA c/ 0.9999968760D-00, -0.9082695004D-02,                          &
        0.7873616869D-04, -0.6111795727D-06,                            &
        0.4388418740D-08, -0.2988388486D-10,                            &
        0.2187442495D-12, -0.1789232111D-14,                            &
        0.1111201803D-16, -0.3099457145D-19/
  t = temp - 273.16
  s = 6.1078000000D+00 / (c(1) + t * (c(2) + t * (c(3) +                &
                               t * (c(4) + t * (c(5) +                  &
                               t * (c(6) + t * (c(7) +                  &
                               t * (c(8) + t * (c(9) +                  &
                               t * c(10)))))))))) ** 8
  satvap = s
  RETURN
  END FUNCTION satvap

SUBROUTINE precw(p,w,u,np)
!
!            02-may-84       sequence numbers removed
!
  DIMENSION p(*),w(*),u(*)
  DATA f/1961.33/
  w1=w(1)
  p1=p(1)
  s=0.
  u(1)=s
!
!            07-jun-84       correction of tape read error
!  do 100 i=2,n
  DO i = 2, np
    w2=w(i)
    p2=p(i)
    dp=ABS(p2-p1)
    s=s+(w1+w2)*dp/f
    u(i)=s
    w1=w2
    p1=p2
  END DO
  RETURN
END SUBROUTINE precw
!
!##################################################################
!##################################################################
!######                                                      ######
!######              FUNCTION VSKINT                         ######
!######                                                      ######
!##################################################################
!##################################################################
!

  FUNCTION vskint(tbb,jcof,juse,jact)
!
! * obtain sfc-skin-temp from vas brightness temperatures
!     jcof=0 for empirical coefficients
!     jcof=1 for climatological(theoretical) coeffs.
!
!     juse=0 to use 'day'(2-chan) or 'night'(3-chan) eqn.
!        according to solar zenith angle
! *******  currently forced out; will use day only if chan 12 missing.
!     juse=1 to use 'day' eqn.
!     juse=2 to use 'night' eqn.
!
!     jact=0 if nothing done
!     jact=1 if 'day' eqn. actually used
!     jact=2 if 'night' eqn. actually used
  COMMON/nav/vlat,vlon,vzen,szen,il,ie,iras,ipic,itime,jtime,jday
  COMMON/tsurfc/tscc(4,2)
  COMMON/tsurfe/tsce(4,2)
  DIMENSION tsc(4,2),tx(3),ic(3),tbb(*)
  DATA ic/7,8,12/,nx/3/,vmisg/999999./
  DO i=1,nx
    j=ic(i)
    tx(i)=tbb(j)
  END DO
  jact=0
  ts=tx(2)
  IF(tx(1) == vmisg) GO TO 180
  IF(tx(2) == vmisg) GO TO 180
  IF(jcof /= 0) GO TO 110
  DO j=1,2
    DO i=1,4
      tsc(i,j)=tsce(i,j)
    END DO
  END DO
  GO TO 130
  110 DO j=1,2
    DO i=1,4
      tsc(i,j)=tscc(i,j)
    END DO
  END DO
  130 j=juse
  IF(j /= 0) GO TO 150
  j=2                       !changed 1 --> 2
!  if(szen.ge.90.) j=2       !deleted line
  150 IF(j == 1) GO TO 160
  IF(tx(3) == vmisg) j=1
  IF(tx(3) < tx(2)-4.) j=1 !changed tx(2) --> tx(2)-4.
  160 jact=j
  ts=tsc(4,j)
  DO i=1,nx
    ts=ts+tx(i)*tsc(i,j)
  END DO
  IF(ts == 0.) ts=vmisg
  180 vskint=ts
  RETURN
  END FUNCTION vskint
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE PCP_MXR                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE pcp_mxr (nx,ny,nz,t_3d,p_3d,zs,hterain,ref_3d,               &
                    cldpcp_type_3d,                                     &
                    qr_3d,qs_3d,qh_3d,                                  &
                    cloudopt,refthr1,refthr2,hgtrefthr,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Perform 3D precipitation mixing ratio (in g/kg) analysis using
!  radar reflectivity data. For rain water, using Kessler (1969)
!  formula:
!            qr(g/kg) = a*(rho*arg)**b                  (1)
!
!  Here arg = Z (mm**6/m**3), and dBZ = 10log10 (arg).
!  Coeffcients a=17300.0, and b=7/4.
!  rho represents the air density.
!
!  For snow and hail, using Rogers and Yau (1989) formula:
!
!            qs(g/kg) = c*(rho*arg)**d                  (2)
!
!  where, c=38000.0,  d=2.2
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (Jian Zhang)
!  06/13/96
!
!  MODIFICATION HISTORY:
!  07/30/97 (J. Zhang)
!           Added precipitation type in the argument list so that
!           mixing ratios of different precip. types can be computed.
!  09/04/97 (J. Zhang)
!            Changed the radar echo thresholds for inserting precip.
!            from radar reflectivities.
!  03/17/2008  K. Brewster
!            Added code for cloudopt=4, cycling, that uses the
!            background ratios of precipitation types to distribute
!            precip from observed reflectivity.
!  04/23/2008  K. Brewster
!            Updated reflectivity thresholding to be consistent
!            with cloud insertion, using refthr1 and refthr2
!
!-----------------------------------------------------------------------
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
  INTEGER, INTENT(IN) :: nx,ny,nz      ! Model grid size
!
  REAL, INTENT(IN) :: t_3d(nx,ny,nz)   ! Temperature (deg. Kelvin)
  REAL, INTENT(IN) :: p_3d(nx,ny,nz)   ! Pressure (Pascal)
  REAL, INTENT(IN) :: zs(nx,ny,nz)     ! Height of scalar points (m)
  REAL, INTENT(IN) :: hterain(nx,ny)   ! Height of scalar points
  REAL, INTENT(IN) :: ref_3d(nx,ny,nz) ! radar reflectivity (dBZ)
  INTEGER, INTENT(IN) :: cldpcp_type_3d(nx,ny,nz)  ! cloud/precip type field

  REAL, INTENT(INOUT) :: qr_3d(nx,ny,nz)  ! rain mixing ratio in (g/kg)
  REAL, INTENT(INOUT) :: qs_3d(nx,ny,nz)  ! snow/sleet/frz-rain mixing ratio
                                          ! in (g/kg)
  REAL, INTENT(INOUT) :: qh_3d(nx,ny,nz)  ! hail mixing ratio in (g/kg)

  INTEGER, INTENT(IN) :: cloudopt
  REAL, INTENT(IN) :: refthr1
  REAL, INTENT(IN) :: refthr2
  REAL, INTENT(IN) :: hgtrefthr
  INTEGER, INTENT(OUT) :: istatus
!
! LOCAL:
! Coefficiants for Z - qr relationship

  REAL, PARAMETER :: a=17300.0
  REAL, PARAMETER :: b=(7.0/4.0)
  REAL, PARAMETER :: c=38000.0
  REAL, PARAMETER :: d=2.2
  REAL, PARAMETER :: br=(1.0/b)
  REAL, PARAMETER :: dr=(1.0/d)
!
  REAL, PARAMETER :: rair = 287.04    ! Gas constant (J/deg/kg)
  REAL, PARAMETER :: thresh_miss = -90.
  INTEGER :: pcptype
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k, iarg
  REAL :: rhobar
  REAL :: zrbkg,zsbkg,zfbkg,zbkg
  REAL :: qw_tot,qf_tot,ztot,zqr,zqf,qsratio
  REAL :: qrold,qsold,qhold,refz,refobs,refbkg,reftest
  REAL :: ref_thresh
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
  istatus=0
!
!-----------------------------------------------------------------------
!
!  Compute the precip mixing ratio in g/kg from radar reflectivity
!  factor following Kessler (1969) or Rogers and Yau (1989).
!
!-----------------------------------------------------------------------
!
  IF(cloudopt < 3) THEN
    DO k = 1,nz-1
      DO j = 1,ny-1
        DO i = 1,nx-1
          IF ((zs(i,j,k)-hterain(i,j)) <= hgtrefthr) THEN
            ref_thresh = refthr1
          ELSE
            ref_thresh = refthr2
          END IF
          IF ( ref_3d(i,j,k) >= ref_thresh ) THEN    ! valid radar refl.
            rhobar = p_3d(i,j,k)/rair/t_3d(i,j,k)
            ztot = 10.0**(0.1*ref_3d(i,j,k))
            iarg = cldpcp_type_3d(i,j,k)
            pcptype = iarg/16              ! precip. type

            IF (pcptype == 1.OR.pcptype == 3) THEN   ! rain or Z R
              qr_3d(i,j,k) = (ztot/a)**br/rhobar
            ELSE IF (pcptype == 2) THEN                   ! snow
              qs_3d(i,j,k) = (ztot/c)**dr/rhobar
            ELSE IF (pcptype == 4.OR.pcptype == 5) THEN   ! hail or sleet
              qh_3d(i,j,k) = (ztot/c)**dr/rhobar
            ELSE                                          ! unknown
              IF(t_3d(i,j,k) > 273.15) THEN
                qr_3d(i,j,k) = (ztot/a)**br/rhobar
                qs_3d(i,j,k) = 0.0
                qh_3d(i,j,k) = 0.0
              ELSE IF( ref_3d(i,j,k) > 45.) THEN
                qr_3d(i,j,k) = 0.0
                qs_3d(i,j,k) = 0.0
                qh_3d(i,j,k) = (ztot/c)**dr/rhobar
              ELSE
                qr_3d(i,j,k) = 0.0
                qs_3d(i,j,k) = (ztot/c)**dr/rhobar
                qh_3d(i,j,k) = 0.0
              END IF
            END IF
          ELSE IF( ref_3d(i,j,k) > thresh_miss ) THEN
            qr_3d(i,j,k) = 0.
            qs_3d(i,j,k) = 0.
            qh_3d(i,j,k) = 0.
          END IF
        END DO ! k
      END DO ! i
    END DO ! j
  ELSE IF (cloudopt == 3) THEN ! warm rain only version
    DO k = 1,nz-1
      DO j = 1,ny-1
        DO i = 1,nx-1
          IF ((zs(i,j,k)-hterain(i,j)) <= hgtrefthr) THEN
            ref_thresh = refthr1
          ELSE
            ref_thresh = refthr2
          END IF
          IF ( ref_3d(i,j,k) >= ref_thresh ) THEN    ! valid radar refl.
            rhobar = p_3d(i,j,k)/rair/t_3d(i,j,k)
            ztot = 10.0**(0.1*ref_3d(i,j,k))
            qr_3d(i,j,k) = (ztot/a)**br/rhobar
            qs_3d(i,j,k) = 0.
            qh_3d(i,j,k) = 0.
          ELSE IF( ref_3d(i,j,k) > thresh_miss ) THEN
            qr_3d(i,j,k) = 0.
            qs_3d(i,j,k) = 0.
            qh_3d(i,j,k) = 0.
          END IF
        END DO ! k
      END DO ! i
    END DO ! j
  ELSE
    DO k = 1,nz-1
      DO j = 1,ny-1
        DO i = 1,nx-1
          IF ( ref_3d(i,j,k) > thresh_miss ) THEN    ! valid radar refl.
            qw_tot = qr_3d(i,j,k)+qs_3d(i,j,k)+qh_3d(i,j,k)
            qf_tot = qs_3d(i,j,k)+qh_3d(i,j,k)
            rhobar = p_3d(i,j,k)/(rair*t_3d(i,j,k))
            IF ((zs(i,j,k)-hterain(i,j)) <= hgtrefthr) THEN
              ref_thresh = refthr1
            ELSE
              ref_thresh = refthr2
            END IF
            IF ( ref_3d(i,j,k) >= ref_thresh ) THEN    ! valid radar refl.
              refobs=ref_3d(i,j,k)
            ELSE IF (ref_3d(i,j,k) > 0. .AND. qw_tot > 0.) THEN
              refz=0.
              IF(qr_3d(i,j,k) > 0.) THEN
                refz=a*((rhobar*qr_3d(i,j,k))**b)
              END IF
              IF(qf_tot > 0.) THEN
                refz=refz+c*((rhobar*qf_tot)**d)
              END IF
              refbkg=10.0*LOG10(max(refz,1.0))
              refobs=ref_3d(i,j,k)
            ELSE
              refobs=0.
            END IF
            IF(refobs > 0.) THEN
              ztot = 10.0**(0.1*refobs)
              iarg = cldpcp_type_3d(i,j,k)
              pcptype = iarg/16              ! precip. type
              qrold=qr_3d(i,j,k)
              qsold=qs_3d(i,j,k)
              qhold=qh_3d(i,j,k)
              IF(qw_tot == 0.) THEN
                IF (pcptype == 1.OR.pcptype == 3) THEN   ! rain or Z R
                  qr_3d(i,j,k) = ((ztot/a)**br)/rhobar
                  qs_3d(i,j,k) = 0.0
                  qh_3d(i,j,k) = 0.0
                ELSE IF (pcptype == 2) THEN                   ! snow
                  qr_3d(i,j,k) = 0.0
                  qs_3d(i,j,k) = ((ztot/c)**dr)/rhobar
                  qh_3d(i,j,k) = 0.0
                ELSE IF (pcptype == 4.OR.pcptype == 5) THEN   ! hail or sleet
                  qr_3d(i,j,k) = 0.0
                  qs_3d(i,j,k) = 0.0
                  qh_3d(i,j,k) = ((ztot/c)**dr)/rhobar
                ELSE                                          ! unknowA
                  IF(t_3d(i,j,k) > 273.15) THEN
                    qr_3d(i,j,k) = ((ztot/a)**br)/rhobar
                    qs_3d(i,j,k) = 0.0
                    qh_3d(i,j,k) = 0.0
                  ELSE IF( ref_3d(i,j,k) > 45.) THEN
                    qr_3d(i,j,k) = 0.0
                    qs_3d(i,j,k) = 0.0
                    qh_3d(i,j,k) = ((ztot/c)**dr)/rhobar
                  ELSE
                    qr_3d(i,j,k) = 0.0
                    qs_3d(i,j,k) = ((ztot/c)**dr)/rhobar
                    qh_3d(i,j,k) = 0.0
                  END IF
                END IF
              ELSE  ! use ratio of background precip types
                IF(qf_tot > 0.0 ) THEN  ! frozen exists
                  IF(qr_3d(i,j,k) > 0.0) THEN  ! mixed
                    zrbkg=a*((rhobar*qr_3d(i,j,k))**b)
                    zfbkg=c*((rhobar*qf_tot)**d)
                    zbkg=zrbkg+zfbkg
                    zqr=ztot*(zrbkg/zbkg)
                    qr_3d(i,j,k) = ((zqr/a)**br)/rhobar
                    zqf=max(0.,(ztot-zqr))
                    IF(qh_3d(i,j,k) == 0.) THEN  ! frozen only snow
                      qs_3d(i,j,k) = ((zqf/c)**dr)/rhobar
                    ELSE IF (qs_3d(i,j,k) == 0.) THEN ! frozen only hail
                      qh_3d(i,j,k) = ((zqf/c)**dr)/rhobar
                    ELSE ! mixed frozen
                      qsratio = qs_3d(i,j,k)/qf_tot
                      qf_tot = ((zqf/c)**dr)/rhobar
                      qs_3d(i,j,k) = qf_tot*qsratio
                      qh_3d(i,j,k) = qf_tot*(1.0-qsratio)
                    END IF
                  ELSE ! frozen only
                    qr_3d(i,j,k) = 0.
                    IF(qh_3d(i,j,k) == 0.) THEN  ! frozen only snow
                      qs_3d(i,j,k) = ((ztot/c)**dr)/rhobar
                    ELSE IF (qs_3d(i,j,k) == 0.) THEN ! frozen only hail
                      qh_3d(i,j,k) = ((ztot/c)**dr)/rhobar
                    ELSE ! mixed frozen
                      qsratio = qs_3d(i,j,k)/qf_tot
                      qf_tot = ((ztot/c)**dr)/rhobar
                      qs_3d(i,j,k) = qf_tot*qsratio
                      qh_3d(i,j,k) = qf_tot*(1.0-qsratio)
                    END IF
                  END IF
                ELSE ! liquid only
                  qr_3d(i,j,k) = ((ztot/a)**br)/rhobar
                  qs_3d(i,j,k) = 0.
                  qh_3d(i,j,k) = 0.
                END IF
              END IF
            ELSE ! refobs = 0.
              qr_3d(i,j,k)=0.
              qs_3d(i,j,k)=0.
              qh_3d(i,j,k)=0.
            END IF
!           IF(j == 264) THEN
!             refz=0.
!             qf_tot=qs_3d(i,j,k)+qh_3d(i,j,k)
!             IF(qr_3d(i,j,k) > 0.) THEN
!               refz=a*((rhobar*qr_3d(i,j,k))**b)
!             END IF
!             IF(qf_tot > 0.) THEN
!               refz=refz+c*((rhobar*qf_tot)**d)
!             END IF
!             reftest=10.0*LOG10(max(refz,1.0))
!             WRITE(47,'(4i4,7f8.3)') &
!               i,j,k,pcptype,qw_tot,qrold,qr_3d(i,j,k), &
!                   qsold,qs_3d(i,j,k),qhold,qh_3d(i,j,k)
!             WRITE(47,'(16x,2f8.3,f8.4)') ref_3d(i,j,k),reftest, &
!                  (ref_3d(i,j,k)-reftest)
!           END IF
          ELSE IF( ref_3d(i,j,k) > thresh_miss ) THEN
            qr_3d(i,j,k) = 0.
            qs_3d(i,j,k) = 0.
            qh_3d(i,j,k) = 0.
          END IF
        END DO ! k
      END DO ! i
    END DO ! j

  END IF
!
!-----------------------------------------------------------------------
!
  istatus = 1
!
  RETURN
END SUBROUTINE pcp_mxr
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE PCP_MXR_FERRIER             ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE pcp_mxr_ferrier (nx,ny,nz,t_3d,p_3d,zs,hterain,ref_3d,      &
                            cldpcp_type_3d,                            &
                            qr_3d,qs_3d,qh_3d,                         &
                            cloudopt,refthr1,refthr2,hgtrefthr,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Perform 3D precipitation mixing ratio (in g/kg) analysis using
!  radar reflectivity data. For rain water, using Ferrier et al (1995)
!  formulation:
!
!
!     For rain water:
!
!          18
!        10   * 720                              1.75
!  Zer = --------------------------- * (rho * qr)
!          1.75      0.75       1.75
!        pi     * N0r     * rhor
!
!
!     For dry snow (t <= 0 C):
!
!
!          18           2                     0.25
!        10  * 720 * |K|                * rhos
!                       ice                                    1.75
!  Zes = ----------------------------------------- * (rho * qs)     t <= 0 C
!          1.75         2          0.75       2
!        pi        * |K|      * N0s     * rhoi
!                     water
!
!
!     For wet snow (t >= 0 C):
!
!
!              18
!            10   * 720                                 1.75
!  Zes =     ----------------------------   * (rho * qs)            t >  0 C
!              1.75      0.75       1.75
!            pi     * N0s     * rhos
!
!
!     For hail water:
!
!
!          /   18                       \ 0.95
!         /  10   * 720                  \              1.6625
!  Zeh = |   ---------------------------- | * (rho * qh)
!         \    1.75      0.75       1.75 /
!          \ pi     * N0h     * rhoh    /
!
!  Here Zx (mm**6/m**3, x=r,s,h), and dBZ = 10log10 (Zx).
!  rho represents the air density, rhor,rhos,rhoh are the density of
!  rain, snow and hail respectively. Other variables are all constants
!  for this scheme, see below.
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (Donghai Wang and Eric Kemp)
!  07/20/2000
!
!  MODIFICATION HISTORY:
!
!  11/09/2000 Keith Brewster
!  Moved some parameters with real-valued exponentiation to be
!  computed at runtime due to compiler complaint.
!
!  04/07/2003 Keith Brewster
!  Restructured code to make more tractable.and consistent with
!  the reflec_ferrier subroutine.
!
!  04/23/2008  K. Brewster
!  Updated reflectivity thresholding to be consistent
!  with cloud insertion, using refthr1 and refthr2
!
!  04/05/2010  K. Brewster
!  Updated code to use n0rain,n0snow,n0hail,rhosnow and rhohail
!  from the input file or background data instead of constant PARAMETER values.
!  Defaults are preserved here as PARAMETERs.
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
  INTEGER, INTENT(IN) :: nx,ny,nz              ! Model grid size
!
  REAL, INTENT(IN) :: t_3d(nx,ny,nz)    ! Temperature (deg. Kelvin)
  REAL, INTENT(IN) :: p_3d(nx,ny,nz)    ! Pressure (Pascal)
  REAL, INTENT(IN) :: zs(nx,ny,nz)      ! Height of scalars (m)
  REAL, INTENT(IN) :: hterain(nx,ny)    ! Terrain height (m)
  REAL, INTENT(IN) :: ref_3d(nx,ny,nz)  ! radar reflectivity (dBZ)
  INTEGER, INTENT(IN) :: cldpcp_type_3d(nx,ny,nz) ! cloud/precip type field

  REAL, INTENT(INOUT) :: qr_3d(nx,ny,nz)  ! rain mixing ratio in (g/kg)
  REAL, INTENT(INOUT) :: qs_3d(nx,ny,nz)  ! snow/sleet/frz-rain mixing ratio
                                          ! in (g/kg)
  REAL, INTENT(INOUT) :: qh_3d(nx,ny,nz)  ! hail mixing ratio in (g/kg)

  INTEGER, INTENT(IN) :: cloudopt
  REAL, INTENT(IN) :: refthr1
  REAL, INTENT(IN) :: refthr2
  REAL, INTENT(IN) :: hgtrefthr
  INTEGER, INTENT(OUT) :: istatus

  REAL,PARAMETER :: ki2 = 0.176 ! Dielectric factor for ice if other
                                !   than melted drop diameters are used.
  REAL,PARAMETER :: kw2=0.93 ! Dielectric factor for water.

  REAL,PARAMETER :: m3todBZ=1.0E+18 ! Conversion factor from m**3 to
                                    !   mm**6 m**-3.
  REAL,PARAMETER :: Zefact=720.0 ! Multiplier for Ze components.
  REAL,PARAMETER :: lg10div=0.10 ! Log10 multiplier (1/10)

  REAL,PARAMETER :: pi=3.1415926 ! Pi.
  REAL,PARAMETER :: N0r_def=8.0E+06 ! Intercept parameter in 1/(m^4) for rain.
  REAL,PARAMETER :: N0s_def=3.0E+06 ! Intercept parameter in 1/(m^4) for snow.
  REAL,PARAMETER :: N0h_def=4.0E+04 ! Intercept parameter in 1/(m^4) for hail.

  REAL,PARAMETER :: N0xpowf=3.0/7.0 ! Power to which N0r,N0s & N0h are
                                    !   raised.
  REAL,PARAMETER :: K2powf=4.0/7.0  ! Power to which K-squared
                                    !  of ice, water are raised
  REAL,PARAMETER :: zkpowf=4.0/7.0  ! Power to which Zk is raised
  REAL,PARAMETER :: zepowf=4.0/7.0  ! Power to which Ze is raised
  REAL,PARAMETER :: zehpowf=(4.0/7.0)*1.0526  ! Power to which Zeh is raised

  REAL,PARAMETER :: rhoi=917.  ! Density of ice (kg m**-3)
  REAL,PARAMETER :: rhor=1000. ! Density of rain (kg m**-3)
  REAL,PARAMETER :: rhos_def=100.  ! Density of snow (kg m**-3)
  REAL,PARAMETER :: rhoh_def=913.  ! Density of hail (kg m**-3)

  REAL,PARAMETER :: rhoipowf=8.0/7.0  ! Power to which rhoi is raised.
  REAL,PARAMETER :: rhospowf=1.0/7.0  ! Power to which rhos is raised.
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  REAL, PARAMETER :: thresh_miss = -90.
  INTEGER :: i,j,k, iarg
  INTEGER :: pcptype
  REAL :: N0r,N0s,N0h,rhos,rhoh
  REAL :: zkconst,zerf,zesnegf,zesposf,zehf,rfract
  REAL :: ze,zer,zeh,zes,rhoinv,tc
  REAL :: ref_thresh
!
!-----------------------------------------------------------------------
!
! Include Files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
  INCLUDE 'mp.inc'
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
!  Intiailize constant factors in the Ze terms for rain, snow and hail,
!  respectively, in Ferrier.
!
!  These are the inverse of those presented in the reflec_ferrier function.
!
!-----------------------------------------------------------------------
!
  istatus=0

  N0r=N0r_def
  N0s=N0s_def
  N0h=N0h_def
  rhos=rhos_def
  rhoh=rhoh_def
  IF (n0rain > 0) N0r=n0rain
  IF (n0snow > 0) N0s=n0snow
  IF (n0hail > 0) N0h=n0hail
  IF (rhosnow > 0) rhos=rhosnow
  IF (rhohail > 0) rhoh=rhohail

  zkconst = (Zefact*m3todBZ) ** zkpowf

  zerf=1000.*(pi * (N0r**N0xpowf) * rhor )/zkconst

  zesnegf=1000.*(pi*(kw2**k2powf)*(N0s**N0xpowf)*(rhoi**rhoipowf)) /    &
          ( zkconst * (ki2**k2powf) * (rhos**rhospowf) )

  zesposf=1000.*( pi * (N0s**N0xpowf) * rhos) / zkconst

  zehf=1000.*( pi * (N0h**N0xpowf) * rhoh) / zkconst

!-----------------------------------------------------------------------
!
!  Compute the precip mixing ratio in g/kg from radar reflectivity
!  factor following Ferrier et al (1995).
!
!-----------------------------------------------------------------------
!
  IF(cloudopt < 3) THEN
  DO k = 2,nz-1
    DO j = 1,ny-1
      DO i = 1,nx-1
        IF ((zs(i,j,k)-hterain(i,j)) <= hgtrefthr) THEN
          ref_thresh = refthr1
        ELSE
          ref_thresh = refthr2
        END IF
        IF ( ref_3d(i,j,k) >= ref_thresh ) THEN    ! valid radar refl.
          rhoinv = (rd*t_3d(i,j,k))/p_3d(i,j,k)
          ze = 10.0**(0.1*ref_3d(i,j,k))
          iarg = cldpcp_type_3d(i,j,k)
          pcptype = iarg/16              ! precip. type
          tc = t_3d(i,j,k) - 273.15

          IF (pcptype == 1) THEN   ! rain
            qr_3d(i,j,k) = rhoinv * zerf * (ze**zepowf)
          ELSE IF (pcptype == 2) THEN                   ! snow
            IF (tc <= 0.0) THEN     !dry snow
              qs_3d(i,j,k) = rhoinv * zesnegf * (ze**zepowf)
            ELSE IF (tc < 5.0) THEN     !wet snow
              rfract=0.20*tc
              zer=rfract*ze
              zes=(1.-rfract)*ze
              qs_3d(i,j,k) = rhoinv * zesposf * (zes**zepowf)
              qr_3d(i,j,k) = rhoinv * zerf * (zer**zepowf)
            ELSE
              qr_3d(i,j,k) = rhoinv * zerf * (ze**zepowf)
            END IF
          ELSE IF (pcptype == 3) THEN   ! ZR
            qr_3d(i,j,k) = rhoinv * zerf * (ze**zepowf)
          ELSE IF (pcptype == 4) THEN   ! sleet
            IF (tc <= 0.0) THEN     ! hail category
              qh_3d(i,j,k) = rhoinv * zehf * (ze**zehpowf)
            ELSE IF( tc < 10. ) THEN
              rfract=0.10*tc
              zer=rfract*ze
              zeh=(1.-rfract)*ze
              qr_3d(i,j,k) = rhoinv * zerf * (zer**zepowf)
              qh_3d(i,j,k) = rhoinv * zehf * (zeh**zehpowf)
            ELSE
              qr_3d(i,j,k) = rhoinv * zerf * (ze**zepowf)
            END IF
          ELSE IF (pcptype == 5) THEN   ! hail
            qh_3d(i,j,k) = rhoinv * zehf * (ze**zehpowf)
          ELSE                                          ! unknown
            IF (tc <= 0.0) THEN     !dry snow
              qs_3d(i,j,k) = rhoinv * zesnegf * (ze**zepowf)
            ELSE IF ( tc < 5.0 ) THEN     !wet snow
              rfract=0.20*tc
              zer=rfract*ze
              zes=(1.-rfract)*ze
              qs_3d(i,j,k) = rhoinv * zesposf * (zes**zepowf)
              qr_3d(i,j,k) = rhoinv * zerf * (zer**zepowf)
            ELSE ! rain
              qr_3d(i,j,k) = rhoinv * zerf * (ze**zepowf)
            END IF
          END IF
        ELSE IF( ref_3d(i,j,k) > thresh_miss) THEN
          qr_3d(i,j,k) = 0.
          qs_3d(i,j,k) = 0.
          qh_3d(i,j,k) = 0.
        END IF
      END DO ! k
    END DO ! i
  END DO ! j
  ELSE  ! warm rain option, no ice, use rain equation
    DO k = 2,nz-1
      DO j = 1,ny-1
        DO i = 1,nx-1
          IF ((zs(i,j,k)-hterain(i,j)) <= hgtrefthr) THEN
            ref_thresh = refthr1
          ELSE
            ref_thresh = refthr2
          END IF
          IF ( ref_3d(i,j,k) >= ref_thresh ) THEN    ! valid radar refl.
            rhoinv = (rd*t_3d(i,j,k))/p_3d(i,j,k)
            ze = 10.0**(0.1*ref_3d(i,j,k))
            qr_3d(i,j,k) = rhoinv * zerf * (ze**zepowf)
            qs_3d(i,j,k) = 0.
            qh_3d(i,j,k) = 0.
          ELSE IF( ref_3d(i,j,k) > thresh_miss) THEN
            qr_3d(i,j,k) = 0.
            qs_3d(i,j,k) = 0.
            qh_3d(i,j,k) = 0.
          END IF
        END DO ! i
      END DO ! j
    END DO ! k
  END IF
  IF (myproc == 0) PRINT*,'Finish Ferrier ...'
!
!-----------------------------------------------------------------------
!
  istatus = 1
!
  RETURN
END SUBROUTINE pcp_mxr_ferrier
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE REFMOSAIC                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
SUBROUTINE refmosaic(nradfil,nx,ny,nz,mx_rad,                           &
           xs,ys,zs,radfname,lvldbg,ref_mos_3d,rhinf,rvinf,             &
           refl,tem1,tem2,istatus)
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  This routine patches together reflectivity fields observed from
!  multiple radars to create one 3D radar reflectivity field.
!  It also fills the data gaps between the radar beams and elevations).
!  These gaps are possible because that the real-time radar data
!  processing only put data to the grid point nearest to the radar gate.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (Jian Zhang)
!  06/04/96.
!
!  MODIFICATION HISTORY:
!  11/01/01 (Keith Brewster)
!  Created refmosaic from rad_patch_fill saving lots of memory when
!  there are more than a couple of radars.
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
  INTEGER :: nradfil              ! # of radars currently available
  INTEGER :: nx,ny,nz             ! model grid size
  INTEGER :: mx_rad
  REAL    :: xs(nx)               ! x-coor at scalar points
  REAL    :: ys(ny)               ! y-coor at scalar points
  REAL    :: zs(nx,ny,nz)         ! z-coor at scalar points
  CHARACTER(LEN=256) radfname(mx_rad)
  INTEGER, INTENT(IN) :: lvldbg
  REAL    :: rhinf
  REAL    :: rvinf
!
!  OUTPUT:
!
  REAL    :: ref_mos_3d (nx,ny,nz) ! combined radar refl. field
  INTEGER :: istatus
!
!  Temporary working array
!
  REAL :: refl(nx,ny,nz)
  REAL :: tem1(nx,ny,nz)
  REAL :: tem2(nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,irad
  INTEGER :: im1,ip1,jm1,jp1,km1,kp1
  INTEGER :: isrcrad,istat_rad
  REAL    :: latrad,lonrad,elvrad
  REAL    :: x_m1,x_p1,y_m1,y_p1,hgt_km1,hgt_kp1,frac
  REAL    :: refl_ip1,refl_im1,refl_jm1,refl_jp1,refl_km1,refl_kp1
  LOGICAL :: got_im1,got_ip1,got_jm1,got_jp1,got_km1,got_kp1
  INTEGER :: strhopt_rad,mapproj_rad
  INTEGER :: nxlg,nylg
  REAL    :: dx_rad,dy_rad,dz_rad,dzmin_rad
  REAL    :: ctrlat_rad,ctrlon_rad,tlat1_rad,tlat2_rad,tlon_rad
  REAL    :: sclfct_rad
  CHARACTER (LEN=5) stnrad
  CHARACTER (LEN=128) warn_string
  INTEGER :: istat

  REAL, ALLOCATABLE :: temrad(:,:,:)
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
   INCLUDE 'grid.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  istatus=0

  warn_string = 'WARNING: Radar data are not used ' //                  &
                   'due to grid or map projection inconsistencies'

  IF (myproc == 0) WRITE(6,'(a)')                                       &
    ' Mosaiking reflectivity from multiple radar observations'

  nxlg = (nx-3) * (nproc_x) + 3
  nylg = (ny-3) * (nproc_y) + 3

!  ALLOCATE(temrad(nxlg,nylg,nz),stat=istat)
!  CALL check_alloc_status(istatus,"refmosaic:temrad")
!
!-----------------------------------------------------------------------
!
!  Combining the multiple radar obs. to produce one set of
!  3D reflectivity field.
!
!-----------------------------------------------------------------------
!
  DO irad = 1,nradfil
!    CALL fill1rad_old(nx,ny,nz,nxlg,nylg,radfname(irad),lvldbg,xs,ys,zs,    &
!                  strhopt_rad,mapproj_rad,dx_rad,dy_rad,dz_rad,         &
!                  dzmin_rad,ctrlat_rad,ctrlon_rad,                      &
!                  tlat1_rad,tlat2_rad,tlon_rad,sclfct_rad,              &
!                  isrcrad,stnrad,latrad,lonrad,elvrad,                  &
!                  tem1,refl,tem2,tem2,temrad,                           &
!                  istat_rad)

    CALL fill1rad(nx,ny,nz,radfname(irad),lvldbg,xs,ys,zs,              &
                  strhopt_rad,mapproj_rad,dx_rad,dy_rad,dz_rad,         &
                  dzmin_rad,ctrlat_rad,ctrlon_rad,                      &
                  tlat1_rad,tlat2_rad,tlon_rad,sclfct_rad,              &
                  isrcrad,stnrad,latrad,lonrad,elvrad,                  &
                  tem1,refl,tem2,tem2,                                  &
                  istat_rad)
!
!-----------------------------------------------------------------------
!
!  For successful read of radar data.
!
!  Check IF the radar data are remapped on the same grid as
!  the one we are using for ADAS.
!
!-----------------------------------------------------------------------
!
    IF (istat_rad == 0 ) THEN
      IF (myproc == 0)                                                  &
        WRITE(6,'(a,a)') ' Successfully read radar datafile ',radfname(irad)
      IF (ABS(strhopt_rad-strhopt) > 1.0E-4 .OR.                        &
              mapproj_rad /= mapproj) THEN
        IF (myproc == 0) THEN
          WRITE(6,*) warn_string
          WRITE(6,*)                                                    &
          'RADAR: strhopt=',strhopt_rad, '  mapproj=',mapproj_rad
          WRITE(6,*) 'ADAS: strhopt=',strhopt, '  mapproj=',mapproj
        ENDIF
        CYCLE
      END IF
!
      IF (ABS(dx_rad-dx) > 1.0E-4 .OR. ABS(dy_rad-dy) > 1.0E-4 .OR.     &
          ABS(dz_rad-dz) > 1.0E-4 .OR.                                  &
          ABS(dzmin_rad-dzmin) > 1.0E-4) THEN
        IF (myproc == 0) THEN
          WRITE(6,*) warn_string
          WRITE(6,*) 'RADAR: dx=',dx_rad, ' dy=',dy_rad,' dz=',dz_rad
          WRITE(6,*) 'ADAS: dx=',dx, ' dy=',dy,' dz=',dz
        ENDIF
        CYCLE
      END IF
!
      IF (mapproj /= 0) THEN
        IF (ABS(ctrlat_rad-ctrlat) > 1.0E-4 .OR.                        &
              ABS(ctrlon_rad-ctrlon) > 1.0E-4) THEN
          IF (myproc == 0) THEN
            WRITE(6,*) warn_string
            WRITE(6,*) 'RADAR: ctrlat=',ctrlat_rad,' ctrlon=',ctrlon_rad
            WRITE(6,*) 'ADAS: ctrlat=',ctrlat, ' ctrlon=',ctrlon
          ENDIF
          CYCLE
        END IF
!
        IF(ABS(sclfct_rad-sclfct) > 1.0E-4) THEN
          IF (myproc == 0) THEN
            WRITE(6,*) warn_string
            WRITE(6,*) 'RADAR: sclfct=',sclfct_rad
            WRITE(6,*) 'ADAS: sclfct=',sclfct
          ENDIF
          CYCLE
        END IF
      END IF
!
      IF (mapproj == 1 .OR. mapproj == 3) THEN
        IF (ABS(tlat1_rad-trulat1) > 1.0E-4 .OR.                        &
              ABS(tlon_rad-trulon) > 1.0E-4) THEN
          IF (myproc == 0) THEN
            WRITE(6,*) warn_string
            WRITE(6,*) 'RADAR: trulat1=',tlat1_rad, ' trulon=',tlon_rad
            WRITE(6,*) 'ADAS: trulat1=',trulat1, ' trulon=',trulon
          ENDIF
          CYCLE
        END IF
      ELSE IF (mapproj == 2) THEN
        IF (ABS(tlat1_rad-trulat1) > 1.0E-4 .OR.                        &
            ABS(tlat2_rad-trulat2) > 1.0E-4 &
            .OR. ABS(tlon_rad-trulon) > 1.0E-4) THEN
          IF (myproc == 0) THEN
            WRITE(6,*) warn_string
            WRITE(6,*) 'RADAR: trulat1=',tlat1_rad,                     &
                   ' trulat2=',tlat2_rad,' trulon=',tlon_rad
            WRITE(6,*) 'ADAS: trulat1=',trulat1,                        &
                   ' trulat2=',trulat2,' trulon=',trulon
          ENDIF
          CYCLE
        END IF
      END IF
!
!-----------------------------------------------------------------------
!
!     For now, mosaic by taking the maximum value
!
!-----------------------------------------------------------------------
!
      IF (myproc == 0)                                                  &
        WRITE(6,'(a,a,a)') ' Bringing radar ',stnrad,' into mosaic '

! OpenMP changed loop order to j,k,i:
!$OMP PARALLEL DO PRIVATE (j,k,i,irad)
      DO j=1,ny
        DO k=1,nz-1
          DO i=1,nx
            ref_mos_3d(i,j,k) = MAX(ref_mos_3d(i,j,k),refl(i,j,k))
          END DO
        END DO
      END DO
    END IF ! read ok
  END DO  !irad

!
!-----------------------------------------------------------------------
!
!  Vertical interpolation to fill the gaps between different
!  elevations.
!
!-----------------------------------------------------------------------
!
! WRITE(6,'(a)') ' Vertically interpolating radar refl.'
!! OpenMP changed loop order to j,i,k:
!!$OMP PARALLEL DO PRIVATE(j,i,k,kp1,got_kp1,refl_kp1,hgt_kp1,km1,got_km1,refl_km1,hgt_km1,frac)
! DO j=1,ny
!   DO i=1,nx
!     DO k=2,nz-2
!       tem1(i,j,k) = ref_mos_3d(i,j,k)
!       IF(ref_mos_3d(i,j,k) < -10.) THEN
!         got_kp1=.FALSE.
!         DO kp1=k+1,nz-1
!           IF (zs(i,j,kp1) > (zs(i,j,k)+rvinf)) EXIT
!           IF (ref_mos_3d(i,j,kp1) >= -10.) THEN
!             refl_kp1 = ref_mos_3d(i,j,kp1)
!             hgt_kp1 = zs(i,j,kp1)
!             got_kp1=.TRUE.
!             EXIT
!           END IF
!         END DO

!         got_km1=.FALSE.
!         DO km1=k-1,2,-1
!           IF(zs(i,j,km1) < (zs(i,j,k)-rvinf)) EXIT
!           IF(ref_mos_3d(i,j,km1) >= -10) THEN
!             refl_km1 = ref_mos_3d(i,j,km1)
!             hgt_km1 = zs(i,j,km1)
!             got_km1=.TRUE.
!             EXIT
!           END IF
!         END DO

!         IF(got_km1 .AND. got_kp1 .AND.                             &
!            (hgt_kp1-hgt_km1) < rvinf ) THEN
!           frac = (zs(i,j,k)-hgt_km1)/(hgt_kp1-hgt_km1)
!           tem1(i,j,k) = refl_kp1*frac + refl_km1*(1.-frac)
!         END IF
!       END IF
!     END DO
!   END DO
! END DO

!
!-----------------------------------------------------------------------
!
!  Horizontal interpolation to further fill gaps between
!  the beams.
!
!-----------------------------------------------------------------------
!
! WRITE(6,'(a)') ' Horizontally interpolating radar refl.'

! DO k=1,nz-1
!   DO i=2,nx-2
!     DO j=2,ny-2
!       ref_mos_3d(i,j,k)=tem1(i,j,k)
!       IF(tem1(i,j,k) < -10.) THEN
!         got_ip1=.FALSE.
!         DO ip1=i+1,nx-1
!           IF (xs(ip1) > (xs(i)+rhinf)) EXIT
!           IF(tem1(ip1,j,k) >= -10.) THEN
!             refl_ip1=tem1(ip1,j,k)
!             x_p1=xs(ip1)
!             got_ip1=.TRUE.
!             EXIT
!           END IF
!         END DO

!         got_im1=.FALSE.
!         DO im1=i-1,1,-1
!           IF (xs(im1) < (xs(i)-rhinf)) EXIT
!           IF(tem1(im1,j,k) >= -10.) THEN
!             refl_im1=tem1(im1,j,k)
!             x_m1=xs(im1)
!             got_im1=.TRUE.
!             EXIT
!           END IF
!         END DO

!         got_jp1=.FALSE.
!         DO jp1=j+1,ny-1
!           IF (ys(jp1) > (ys(j)+rhinf)) EXIT
!           IF(tem1(i,jp1,k) >= -10.) THEN
!             refl_jp1=tem1(i,jp1,k)
!             y_p1=ys(jp1)
!             got_jp1=.TRUE.
!             EXIT
!           END IF
!         END DO

!         got_jm1=.FALSE.
!         DO jm1=j-1,1,-1
!           IF (ys(jm1) < (ys(j)-rhinf)) EXIT
!           IF(tem1(i,jm1,k) >= -10.) THEN
!             refl_jm1=tem1(i,jm1,k)
!             y_m1=ys(jm1)
!             got_jm1=.TRUE.
!             EXIT
!           END IF
!         END DO

!         IF(got_im1 .AND. got_ip1 .AND. (x_p1-x_m1) < rhinf .AND.  &
!            got_jm1 .AND. got_jp1 .AND. (y_p1-y_m1) < rhinf ) THEN
!           frac = (xs(i)-x_m1)/(x_p1-x_m1)
!           ref_mos_3d(i,j,k) = refl_ip1*frac + refl_im1*(1-frac)
!           frac = (ys(j)-y_m1)/(y_p1-y_m1)
!           ref_mos_3d(i,j,k) = (refl_jp1*frac + refl_jm1*(1-frac)    &
!                                    + ref_mos_3d(i,j,k))*0.5
!         ELSE IF (got_im1 .AND. got_ip1 .AND.                      &
!                  (x_p1-x_m1) < rhinf) THEN
!           frac = (xs(i)-x_m1)/(x_p1-x_m1)
!           ref_mos_3d(i,j,k) = refl_ip1*frac + refl_im1*(1-frac)
!         ELSE IF (got_jm1 .AND. got_jp1 .AND.                      &
!                  (y_p1-y_m1) < rhinf) THEN
!           frac = (ys(j)-y_m1)/(y_p1-y_m1)
!           ref_mos_3d(i,j,k) = refl_jp1*frac + refl_jm1*(1-frac)
!         END IF
!       END IF
!     END DO
!   END DO
! END DO

!
!-----------------------------------------------------------------------
!
! Done.
!
!-----------------------------------------------------------------------
!
  istatus=1
!  DEALLOCATE(temrad)

  RETURN
END SUBROUTINE refmosaic
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE READRAD                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE readrad(nx,ny,nz,isrcrad,stnrad                              &
           ,latrad,lonrad,elvrad                                        &
           ,gridvel,gridref,gridnyq,gridtim                             &
           ,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Reads radar data remapped on the ARPS grid.
!  This routine requires the remapping to occur on the same grid
!  as the analysis.
!
!  THIS CODE IS NOT CALLED ANYWHERE!!!
!
!-----------------------------------------------------------------------
!
!  AUTHOR:  Jian Zhang
!  05/1996  Read the remapped radar data which was written by the
!           corresponding output routine "wrtrad" in remaplib.f.
!
!  MODIFICATION HISTORY:
!  03/19/97  J. Zhang
!            Added a line of error message when there is trouble
!            reading a radar file.
!  04/03/97  J. Zhang
!            Added the option of reading the data file created
!            from "WTRADCOL".  Added output for the remapping
!            parameters in the radar file (e.g., strhopt,mapproj,
!            dx,dy,dz,dzmin,ctrlat,ctrlon,tlat1,tlat2,tlon,scale)
!  04/07/97  J. Zhang
!            Added  the QC for the case when i,j,k outside the model
!            domain
!  04/09/97  J. Zhang
!            Added the Initializations for gridref, girdvel...
!  04/11/97  J. Zhang
!            Include dims.inc for nx,ny,nz
!  04/14/97  J. Zhang
!            Added message output for the case when actual # of
!            radar files exceeds the maximum allowed number in the
!            ADAS include file.  When that happens, the program will
!            stop.
!  03/31/98  J. Zhang
!            Deleted the option for reading the radar data file
!            created from "WRTRAD".
!
!-----------------------------------------------------------------------
!
!  INCLUDE:  (from dims.inc and adas.inc)
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    mx_rad     maximum number of radars
!
!    nradfil    number of radar files
!    radfname   file name for radar datasets
!
!  OUTPUT:
!
!    isrcrad  index of radar source
!    stnrad   radar site name    character*4
!    latrad   latitude of radar  (degrees N)
!    lonrad   longitude of radar (degrees E)
!    elvrad   elevation of feed horn of radar (m MSL)
!
!    gridvel  radial velocity on ARPS grid
!    gridref  reflectivity on ARPS grid
!    gridnyq  nyquist velocity on ARPS grid
!    gridtim  observation time at ARPS grid
!
!    istatus  status indicator
!
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
  INCLUDE 'adas.inc'     ! ADAS parameters
!
!-----------------------------------------------------------------------
!
!  INPUT:
  INTEGER :: nx,ny,nz    ! the ARPS grid size
!
!  LOCAL:
!
  REAL :: readk(nz)
  REAL :: readhgt(nz)
  REAL :: readref(nz)
  REAL :: readvel(nz)
  REAL :: readnyq(nz)
  REAL :: readtim(nz)
!
  INTEGER :: kntref(nz)
  INTEGER :: kntvel(nz)
  INTEGER :: iradvr
  INTEGER :: nradvr
!
!  OUTPUT:
  INTEGER :: istatus
!
!  OUTPUT:  ARPS radar arrays
  REAL :: gridvel(nx,ny,nz,mx_rad)
  REAL :: gridref(nx,ny,nz,mx_rad)
  REAL :: gridnyq(nx,ny,nz,mx_rad)
  REAL :: gridtim(nx,ny,nz,mx_rad)
!
!  OUTPUT:  Radar site variables
  INTEGER :: isrcrad(0:mx_rad)
  CHARACTER (LEN=5) :: stnrad(mx_rad)
  REAL :: latrad(mx_rad)
  REAL :: lonrad(mx_rad)
  REAL :: elvrad(mx_rad)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=4)   :: stn
  CHARACTER (LEN=80)  :: runname
  CHARACTER (LEN=256) :: fname
  INTEGER :: ireftim,itime,vcpnum,idummy
  INTEGER :: iradfmt,strhopt,mapprin
  INTEGER :: nchanl,ierr
  INTEGER :: iyr, imon, idy, ihr, imin, isec
  INTEGER :: i,j,k,krad,kk,ipt,klev

  REAL :: dxin,dyin,dzin,dzminin,ctrlatin
  REAL :: ctrlonin,tlat1in,tlat2in,tlonin,scalin,rdummy
  REAL :: xrd,yrd,gridlat,gridlon,elev
!
!-----------------------------------------------------------------------
!
!  Common block that stores remapping parameters for the radar
!  data file.
!
!-----------------------------------------------------------------------
!
  COMMON/remapfactrs_rad/strhopt,mapprin
  COMMON/remapfactrs_rad2/dxin,dyin,dzin,dzminin,                       &
           ctrlatin,ctrlonin,tlat1in,tlat2in,tlonin,scalin
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
  istatus=0
  IF (nradfil > mx_rad) THEN
    WRITE(6,'(a,i3,a,i3/a)')                                            &
        ' ERROR: nradfil ',nradfil,' exceeds mx_rad dimension',         &
        mx_rad,' please increase MX_RAD in the .inc file'
    PRINT*,' ABORTING from READRAD......'
    STOP
  END IF
  IF(nradfil < 1) THEN
    WRITE(6,*) 'No radar data available. Returning from READRAD...'
    RETURN
  END IF
!
!-----------------------------------------------------------------------
!
!  Initializations
!
!-----------------------------------------------------------------------
!
! OpenMP changed loop order to j,krad,k,i:
!$OMP PARALLEL DO PRIVATE(j,krad,k,i)
  DO j=1,ny
    DO krad = 1, mx_rad
      DO k=1,nz
        DO i=1,nx
          gridref(i,j,k,krad)=-9999.
          gridvel(i,j,k,krad)=-9999.
          gridnyq(i,j,k,krad)=-9999.
          gridtim(i,j,k,krad)=-9999.
        END DO
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Loop through all radars
!
!-----------------------------------------------------------------------
!
  DO krad = 1, nradfil

    fname=radfname(krad)
    CALL asnctl ('NEWLOCAL', 1, ierr)
    CALL asnfile(fname, '-F f77 -N ieee', ierr)

    CALL getunit( nchanl )
    OPEN(UNIT=nchanl,FILE=trim(fname),ERR=399,                          &
         FORM='unformatted',STATUS='old')
!
!-----------------------------------------------------------------------
!
!  Read radar description variables
!
!-----------------------------------------------------------------------
!
    istatus=1
    isrcrad(krad)=1

    READ(nchanl) stn
    stnrad(krad)=stn
    READ(nchanl) ireftim,itime,vcpnum,idummy,idummy,                    &
               idummy,idummy,idummy,idummy,idummy

    CALL abss2ctim(itime, iyr, imon, idy, ihr, imin, isec )
    iyr=MOD(iyr,100)
    WRITE(6,'(/a,i2.2,a,i2.2,a,i2.2,1X,i2.2,a,i2.2,a)')                 &
       'Reading remapped raw radar data for: ',                         &
       imon,'/',idy,'/',iyr,ihr,':',imin,' UTC'

    READ(nchanl) runname
    READ(nchanl) iradfmt,strhopt,mapprin,idummy,idummy,                 &
               idummy,idummy,idummy,idummy,idummy

    READ(nchanl) dxin,dyin,dzin,dzminin,ctrlatin,                       &
           ctrlonin,tlat1in,tlat2in,tlonin,scalin,                      &
           latrad(krad),lonrad(krad),elvrad(krad),                      &
           rdummy,rdummy

    DO k=1,nz
      kntref(k) = 0
      kntvel(k) = 0
    END DO

    READ(nchanl) nradvr,iradvr

    DO ipt=1,(nx*ny)

      READ(nchanl,END=51) i,j,xrd,yrd,                                  &
                     gridlat,gridlon,elev,klev
      READ(nchanl,END=52) (readk(kk),kk=1,klev)
      READ(nchanl,END=52) (readhgt(kk),kk=1,klev)
      READ(nchanl,END=52) (readref(kk),kk=1,klev)
      READ(nchanl,END=52) (readvel(kk),kk=1,klev)
      READ(nchanl,END=52) (readnyq(kk),kk=1,klev)
      READ(nchanl,END=52) (readtim(kk),kk=1,klev)

      IF(i <= nx.AND.i >= 1 .AND. j <= ny.AND.j >= 1) THEN
        DO kk=1,klev
          k=nint(readk(kk))
          IF(k <= nz.AND.k >= 1) THEN
            gridref(i,j,k,krad)=readref(kk)
            gridvel(i,j,k,krad)=readvel(kk)
            gridnyq(i,j,k,krad)=readnyq(kk)
            gridtim(i,j,k,krad)=readtim(kk)
            IF (gridref(i,j,k,krad) > -200. .AND. gridref(i,j,k,krad) < 200.) &
                kntref(k)=kntref(k)+1
            IF (gridvel(i,j,k,krad) > -200. .AND. gridvel(i,j,k,krad) < 200.) &
                kntvel(k)=kntvel(k)+1
          END IF  ! 1 < k < nz
        END DO  ! kk = 1, klev
      END IF  ! 1 < i < nx  & 1 < j < ny

    END DO  ! ipt = 1, nx*ny

    51      CONTINUE
    ipt=ipt-1
    WRITE(6,'(a,i6,a)') ' End of file reached after reading',           &
                       ipt,' columns'
    GO TO 55
    52      CONTINUE
    WRITE(6,'(a,i6,a)') ' End of file reached while reading',           &
                       ipt,' column'
    55      CONTINUE
!
!-----------------------------------------------------------------------
!
!  Write statistics
!
!-----------------------------------------------------------------------
!
    WRITE(6,'(a)') '  k   n ref    n vel'
    DO k=1,nz
      WRITE(6,'(i3,2i10)') k,kntref(k),kntvel(k)
    END DO

    CLOSE(nchanl)
    CALL retunit( nchanl )

    GO TO 400

    399     CONTINUE
    PRINT*,'Error reading the radar file:',fname

    400     CONTINUE
  END DO  ! KRAD = 1, nradfil

  RETURN
END SUBROUTINE readrad
!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE READ1RAD                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE read1rad(nx,ny,nz,nxlg,nylg,radfname,lvldbg,                 &
                  strhopt_rad,mapproj_rad,dx_rad,dy_rad,dz_rad,         &
                  dzmin_rad,ctrlat_rad,ctrlon_rad,                      &
                  tlat1_rad,tlat2_rad,tlon_rad,sclfct_rad,              &
                  isrcrad,stnrad,latrad,lonrad,elvrad,                  &
                  gridvel,gridref,gridnyq,gridtim,tem1,                 &
                  istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Reads radar data remapped on the ARPS grid.
!  This routine requires the remapping to occur on the same grid
!  as the analysis.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:  Jian Zhang
!  05/1996  Read the remapped radar data which was written by the
!           corresponding output routine "wrtrad" in remaplib.f.
!
!  MODIFICATION HISTORY:
!  03/19/97  J. Zhang
!            Added a line of error message when there is trouble
!            reading a radar file.
!  04/03/97  J. Zhang
!            Added the option of reading the data file created
!            from "WTRADCOL".  Added output for the remapping
!            parameters in the radar file (e.g., strhopt,mapproj,
!            dx,dy,dz,dzmin,ctrlat,ctrlon,tlat1,tlat2,tlon,scale)
!  04/07/97  J. Zhang
!            Added  the QC for the case when i,j,k outside the model
!            domain
!  04/09/97  J. Zhang
!            Added the Initializations for gridref, girdvel...
!  04/11/97  J. Zhang
!            Include dims.inc for nx,ny,nz
!  04/14/97  J. Zhang
!            Added message output for the case when actual # of
!            radar files exceeds the maximum allowed number in the
!            ADAS include file.  When that happens, the program will
!            stop.
!  03/31/98  J. Zhang
!            Deleted the option for reading the radar data file
!            created from "WRTRAD".
!  11/01/01  K. Brewster
!            Modified readrad to read1rad to read only one radar file
!            at a time to save memory.
!  02/02/04  K. Thomas
!            Change so a missing HDF4 isn't a fatal error.
!  02/18/04  K. Thomas
!            Allocated I/O unit not returned when binary file doesn't exit.
!            This update fixes the problem.
!  03/05/07  K. W. Thomas
!            MPI.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nxlg     Number of grid points in the x-direction (large grid/MPI)
!    nylg     Number of grid points in the y-direction (large grid/MPI)
!    radfname   file name for radar datasets
!
!  OUTPUT:
!
!    isrcrad  index of radar source
!    stnrad   radar site name    character*4
!    latrad   latitude of radar  (degrees N)
!    lonrad   longitude of radar (degrees E)
!    elvrad   elevation of feed horn of radar (m MSL)
!
!    gridvel  radial velocity on ARPS grid
!    gridref  reflectivity on ARPS grid
!    gridnyq  nyquist velocity on ARPS grid
!    gridtim  observation time at ARPS grid
!
!    istatus  status indicator
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INCLUDE 'mp.inc'
!
  INTEGER :: nx,ny,nz    ! the ARPS grid size
  INTEGER :: nxlg,nylg
  CHARACTER(LEN=*)    :: radfname
  INTEGER, INTENT(IN) :: lvldbg   ! Added by Y. Wang
  INTEGER :: strhopt_rad,mapproj_rad
  REAL    :: dx_rad,dy_rad,dz_rad,dzmin_rad
  REAL    :: ctrlat_rad,ctrlon_rad,tlat1_rad,tlat2_rad,tlon_rad
  REAL    :: sclfct_rad
!
!  OUTPUT:  Radar site variables
!
  INTEGER :: isrcrad
  CHARACTER (LEN=5) :: stnrad
  REAL :: latrad
  REAL :: lonrad
  REAL :: elvrad
!
!  OUTPUT:  ARPS radar arrays
!
  REAL :: gridvel(nx,ny,nz)
  REAL :: gridref(nx,ny,nz)
  REAL :: gridnyq(nx,ny,nz)
  REAL :: gridtim(nx,ny,nz)
  REAL :: tem1(nxlg,nylg,nz)
  INTEGER :: istatus
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  REAL :: readk(nz)
  REAL :: readhgt(nz)
  REAL :: readref(nz)
  REAL :: readvel(nz)
  REAL :: readnyq(nz)
  REAL :: readtim(nz)
!
  INTEGER, ALLOCATABLE :: kntref(:)
  INTEGER, ALLOCATABLE :: kntvel(:)
!
  INTEGER :: mxradvr,nradvr
  PARAMETER(mxradvr=10)
  INTEGER :: iradvr(mxradvr)

  CHARACTER (LEN=4) :: stn
  CHARACTER (LEN=80) :: runname
  INTEGER :: ireftim,itime,vcpnum,idummy
  INTEGER :: iradfmt,strhopt,mapprin
  INTEGER :: nchanl,ierr
  INTEGER :: iyr, imon, idy, ihr, imin, isec
  INTEGER :: i,j,k,krad,kk,ipt,klev
  INTEGER :: istat

  REAL :: xrd,yrd,gridlat,gridlon,elev,rdummy

  LOGICAL :: verbose = .FALSE.
  INTEGER :: ips, ipe, jps, jpe

  REAL,  ALLOCATABLE :: tem2(:,:,:), tem3(:,:,:), tem4(:,:,:)
!
!-----------------------------------------------------------------------
!
! Temporary  hdf arrays
!
!-----------------------------------------------------------------------
!
  INTEGER (KIND=selected_int_kind(4)), ALLOCATABLE :: itmp(:,:,:)
  REAL,                                ALLOCATABLE :: hmax(:), hmin(:)
  INTEGER :: lens,dmpfmt,sd_id,isource
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  istatus=0

  IF (lvldbg > 90) THEN
    verbose = .TRUE.

    IF (mp_opt > 0) THEN
      ips = 2           ! Patch start index
      ipe = nx-2        ! Patch end index
      jps = 2
      jpe = ny-2
      IF (loc_x == 1)       ips = 1    ! To cover the same domain as non-MPI runs
      IF (loc_x == nproc_x) ipe = nx
      IF (loc_y == 1)       jps = 1
      IF (loc_y == nproc_y) jpe = ny
    ELSE
      ips = 1             ! Domain start index
      ipe = nx            ! Domain end index
      jps = 1
      jpe = ny
    END IF

    ALLOCATE(kntref(nz))
    ALLOCATE(kntvel(nz))

    DO k=1,nz
      kntref(k)=0
      kntvel(k)=0
    END DO
  END IF
!
!-----------------------------------------------------------------------
!
!  Initializations
!
!-----------------------------------------------------------------------
!
!
! OpenMP changed loop order to j,k,i:
!$OMP PARALLEL DO PRIVATE(j,k,i)
  DO j=1,ny
    DO k=1,nz
      DO i=1,nx
        gridref(i,j,k)=-9999.
        gridvel(i,j,k)=-9999.
        gridnyq(i,j,k)=-9999.
        gridtim(i,j,k)=-9999.
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
! Open file
!
!-----------------------------------------------------------------------
!
  CALL asnctl ('NEWLOCAL', 1, ierr)
  CALL asnfile(radfname, '-F f77 -N ieee', ierr)

  lens=LEN(trim(radfname))
  IF(radfname(lens-3:lens)=='hdf4')THEN
    dmpfmt=2
  ELSE
    dmpfmt=1
  ENDIF
  IF (myproc == 0) WRITE(6,'(1x,3a,I2)') 'radfname = ',TRIM(radfname),', dmpfmt = ',dmpfmt

  IF(dmpfmt==1)THEN
    isrcrad=1

    ALLOCATE (tem2(nxlg,nylg,nz),stat=istatus)
    CALL check_alloc_status(istatus,"read1rad:tem2")
    ALLOCATE (tem3(nxlg,nylg,nz),stat=istatus)
    CALL check_alloc_status(istatus,"read1rad:tem3")
    ALLOCATE (tem4(nxlg,nylg,nz),stat=istatus)
    CALL check_alloc_status(istatus,"read1rad:tem4")

    IF (myproc == 0) THEN
      CALL getunit( nchanl )
      OPEN(UNIT=nchanl,FILE=trim(radfname),ERR=399,                     &
           FORM='unformatted',STATUS='old')
!
!-----------------------------------------------------------------------
!
!  Read radar description variables
!
!-----------------------------------------------------------------------
!
      READ(nchanl) stn
      READ(nchanl) ireftim,itime,vcpnum,idummy,idummy,                  &
                   idummy,idummy,idummy,idummy,idummy

      CALL abss2ctim(itime, iyr, imon, idy, ihr, imin, isec )
      iyr=MOD(iyr,100)
      WRITE(6,'(/a,i2.2,a,i2.2,a,i2.2,1X,i2.2,a,i2.2,a)')               &
              'Reading remapped raw radar data for: ',                  &
              imon,'/',idy,'/',iyr,ihr,':',imin,' UTC'

      READ(nchanl) runname
      READ(nchanl) iradfmt,strhopt_rad,mapproj_rad,idummy,idummy,       &
                   idummy,idummy,idummy,idummy,idummy

      READ(nchanl) dx_rad,dy_rad,dz_rad,dzmin_rad,ctrlat_rad,           &
                   ctrlon_rad,tlat1_rad,tlat2_rad,tlon_rad,sclfct_rad,  &
                   latrad,lonrad,elvrad,                                &
                   rdummy,rdummy

      READ(nchanl) nradvr,iradvr

      DO ipt=1,(nxlg*nylg)

        READ(nchanl,END=51) i,j,xrd,yrd,                                &
                            gridlat,gridlon,elev,klev
        READ(nchanl,END=52) (readk(kk),kk=1,klev)
        READ(nchanl,END=52) (readhgt(kk),kk=1,klev)
        READ(nchanl,END=52) (readref(kk),kk=1,klev)
        READ(nchanl,END=52) (readvel(kk),kk=1,klev)
        READ(nchanl,END=52) (readnyq(kk),kk=1,klev)
        READ(nchanl,END=52) (readtim(kk),kk=1,klev)

        IF(i <= nxlg.AND.i >= 1 .AND. j <= nylg.AND.j >= 1) THEN
          DO kk=1,klev
            k=nint(readk(kk))
            IF(k <= nz.AND.k >= 1) THEN
              tem1(i,j,k)=readref(kk)
              tem2(i,j,k)=readvel(kk)
              tem3(i,j,k)=readnyq(kk)
              tem4(i,j,k)=readtim(kk)
            END IF  ! 1 < k < nz
          END DO  ! kk = 1, klev
        END IF  ! 1 < i < nx  & 1 < j < ny

      END DO  ! ipt = 1, nx*ny

      51  CONTINUE
      ipt=ipt-1
      WRITE(6,'(a,i6,a)') ' End of file reached after reading',         &
                          ipt,' columns'
      GO TO 55
      52  CONTINUE
      WRITE(6,'(a,i6,a)') ' End of file reached while reading',         &
                          ipt,' column'
      55  CONTINUE
!
!-----------------------------------------------------------------------
!
!  Write statistics
!
!-----------------------------------------------------------------------
!
      PRINT*,'finishing reading binary radar file'
      CLOSE(nchanl)
      CALL retunit( nchanl )

      GO TO 400

      399     CONTINUE
      PRINT*,'Error reading the radar file:',radfname
      CALL retunit( nchanl )
      istatus = -1

      400     CONTINUE
    END IF

    CALL mpupdatec(istatus,1)
    IF (istatus /= 0) RETURN

    CALL mpupdatec(stn,4)
    stnrad=stn
    CALL mpupdatei(strhopt_rad,1)
    CALL mpupdatei(mapproj_rad,1)
    CALL mpupdater(dx_rad,1)
    CALL mpupdater(dy_rad,1)
    CALL mpupdater(dz_rad,1)
    CALL mpupdater(dzmin_rad,1)
    CALL mpupdater(ctrlat_rad,1)
    CALL mpupdater(ctrlon_rad,1)
    CALL mpupdater(tlat1_rad,1)
    CALL mpupdater(tlat2_rad,1)
    CALL mpupdater(tlon_rad,1)
    CALL mpupdater(sclfct_rad,1)
    CALL mpupdater(latrad,1)
    CALL mpupdater(lonrad,1)
    CALL mpupdater(elvrad,1)

    CALL mpisplit3d(tem1,nx,ny,nz,gridref)
    CALL mpisplit3d(tem2,nx,ny,nz,gridvel)
    CALL mpisplit3d(tem3,nx,ny,nz,gridnyq)
    CALL mpisplit3d(tem4,nx,ny,nz,gridtim)

    DEALLOCATE(tem2,tem3,tem4)

  ELSE !HDF4 file

    IF (myproc == 0) THEN
      ALLOCATE (itmp(nxlg,nylg,nz),stat=istatus)
      CALL check_alloc_status(istatus,"read1rad:itmp")
      ALLOCATE (hmax(nz),stat=istatus)
      CALL check_alloc_status(istatus,"read1rad:hmax")
      ALLOCATE (hmin(nz),stat=istatus)
      CALL check_alloc_status(istatus,"read1rad:hmin")

      CALL hdfopen(trim(radfname), 1, sd_id)
      IF (sd_id < 0) THEN
        WRITE (6,*) "READ1RAD: ERROR opening ",trim(radfname)," for reading."
        istatus = -1
      ELSE
        istatus = 0
      END IF
    END IF
    CALL mpupdatei(istatus,1)

    IF (istatus == 0) THEN   ! Open file sucessfully

      IF (myproc == 0) THEN
        CALL hdfrdc(sd_id, 4, 'radid', stn, istatus)
        stnrad=stn
        CALL hdfrdi(sd_id, 'ireftim', ireftim, istatus)
        CALL hdfrdi(sd_id, 'itime', itime, istatus)
        CALL hdfrdi(sd_id, 'vcpnum', vcpnum, istatus)
        CALL hdfrdi(sd_id, 'isource', isource, istatus)
        CALL hdfrdc(sd_id, 40, 'runname', runname, istatus)
        CALL hdfrdi(sd_id, 'iradfmt', iradfmt, istatus)
        CALL hdfrdi(sd_id, 'strhopt', strhopt_rad, istatus)
        CALL hdfrdi(sd_id, 'mapproj', mapproj_rad, istatus)
        CALL hdfrdr(sd_id, 'dx', dx_rad, istatus)
        CALL hdfrdr(sd_id, 'dy', dy_rad, istatus)
        CALL hdfrdr(sd_id, 'dz', dz_rad, istatus)
        CALL hdfrdr(sd_id, 'dzmin', dzmin_rad, istatus)
        CALL hdfrdr(sd_id, 'ctrlat', ctrlat_rad, istatus)
        CALL hdfrdr(sd_id, 'ctrlon', ctrlon_rad, istatus)
        CALL hdfrdr(sd_id, 'trulat1', tlat1_rad, istatus)
        CALL hdfrdr(sd_id, 'trulat2', tlat2_rad, istatus)
        CALL hdfrdr(sd_id, 'trulon', tlon_rad, istatus)
        CALL hdfrdr(sd_id, 'sclfct', sclfct_rad, istatus)
        CALL hdfrdr(sd_id, 'latrad', latrad, istatus)
        CALL hdfrdr(sd_id, 'lonrad', lonrad, istatus)
        CALL hdfrdr(sd_id, 'elvrad', elvrad, istatus)
        CALL hdfrdi(sd_id, 'nradvr', nradvr, istatus)
        CALL hdfrd1di(sd_id,'iradvr', mxradvr,iradvr,istatus)

        IF (verbose) PRINT *, ' Got nradvr,iradvr: ',nradvr,iradvr
        CALL abss2ctim(itime, iyr, imon, idy, ihr, imin, isec )
        iyr=MOD(iyr,100)
        WRITE(6,'(/1x,a,i2.2,a,i2.2,a,i2.2,1X,i2.2,a,i2.2,a)')          &
          'Reading remapped raw radar data for: ',                      &
          imon,'/',idy,'/',iyr,ihr,':',imin,' UTC'
      ENDIF
      CALL mpupdatei(istatus,1)

      IF (istatus == 0 ) THEN  ! read header sucessfully
        CALL mpupdatei(isrcrad,1)
        CALL mpupdatec(stnrad,5)
        CALL mpupdatei(strhopt_rad,1)
        CALL mpupdatei(mapproj_rad,1)
        CALL mpupdater(dx_rad,1)
        CALL mpupdater(dy_rad,1)
        CALL mpupdater(dz_rad,1)
        CALL mpupdater(dzmin_rad,1)
        CALL mpupdater(ctrlat_rad,1)
        CALL mpupdater(ctrlon_rad,1)
        CALL mpupdater(tlat1_rad,1)
        CALL mpupdater(tlat2_rad,1)
        CALL mpupdater(tlon_rad,1)
        CALL mpupdater(sclfct_rad,1)
        CALL mpupdater(latrad,1)
        CALL mpupdater(lonrad,1)
        CALL mpupdater(elvrad,1)

        IF (myproc == 0) THEN
          CALL hdfrd3d(sd_id,"gridref",nxlg,nylg,nz,tem1,istatus,itmp,hmax,hmin)
        END IF
        CALL mpupdatei(istatus,1)
        IF (istatus /= 0) GO TO 115
        CALL mpisplit3d(tem1,nx,ny,nz,gridref)

        IF (myproc == 0) THEN
          CALL hdfrd3d(sd_id,"gridvel",nxlg,nylg,nz,tem1,istatus,itmp,hmax,hmin)
        END IF
        CALL mpupdatei(istatus,1)
        IF (istatus /= 0) GO TO 115
        CALL mpisplit3d(tem1,nx,ny,nz,gridvel)

        IF (myproc == 0) THEN
          CALL hdfrd3d(sd_id,"gridnyq",nxlg,nylg,nz,tem1,istatus,itmp,hmax,hmin)
        END IF
        CALL mpupdatei(istatus,1)
        IF (istatus /= 0) GO TO 115
        CALL mpisplit3d(tem1,nx,ny,nz,gridnyq)

        IF (myproc == 0) THEN
          CALL hdfrd3d(sd_id,"gridtim",nxlg,nylg,nz,tem1,istatus,itmp,hmax,hmin)
          CALL hdfclose(sd_id,istatus)
        END IF
        CALL mpupdatei(istatus,1)
        IF (istatus /= 0) GO TO 115
        CALL mpisplit3d(tem1,nx,ny,nz,gridtim)

        istatus=0
      END IF
    END IF

    IF (myproc == 0) DEALLOCATE(hmin,hmax,itmp)

  END IF

  IF (verbose .AND. istatus == 0 ) THEN

    DO j=jps,jpe
      DO k=1,nz
        DO i=ips,ipe
          IF (gridref(i,j,k) > -200. .AND. gridref(i,j,k) < 200.)     &
            kntref(k)=kntref(k)+1
          IF (gridvel(i,j,k) > -200. .AND. gridvel(i,j,k) < 200.)     &
            kntvel(k)=kntvel(k)+1
        END DO
      END DO
    END DO

    DO k=1,nz
      CALL mptotali(kntref(k))
      CALL mptotali(kntvel(k))
    END DO

    IF (myproc == 0) THEN
      WRITE(6,'(a)') '  k       n ref      n vel'
      DO k = 1,nz
        WRITE(6,'(i3,2i10)') k,kntref(k),kntvel(k)
      END DO
    END IF

    DEALLOCATE(kntref)
    DEALLOCATE(kntvel)

  END IF

  RETURN

  !
  !  Destination for hdf read error
  !

  115   CONTINUE
  WRITE(6,'(/a/)') ' Error reading data in HDFREAD'

  istatus=-11

  IF (myproc == 0) DEALLOCATE(hmin,hmax, itmp)

  DEALLOCATE(kntref)
  DEALLOCATE(kntvel)

  RETURN

END SUBROUTINE read1rad
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GET_VIS                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE process_vis(solar_alt,nx,ny,                                 &
           cloud_frac_vis_a,albedo,r_missing,lvldbg,istatus)
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read satellite visible imagery data from LAPS ".lvd" file.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (Jian Zhang)
!   04/1996   Based on the LAPS cloud analysis code of 1995.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INCLUDE 'mp.inc'
!  INPUT:
  INTEGER, INTENT(IN) :: nx,ny
  REAL,    INTENT(IN) :: solar_alt(nx,ny)
  REAL,    INTENT(IN) :: r_missing
  INTEGER, INTENT(IN) :: lvldbg
!
!  OUTPUT:
!
  REAL :: albedo(nx,ny)
  REAL :: cloud_frac_vis_a(nx,ny)
!
!  LOCAL
!
  INTEGER :: ihist_alb(-10:20)
  INTEGER :: ihist_frac_sat(-10:20)
  INTEGER :: ih_cf_sat(-10:20)
!
!-----------------------------------------------------------------------
!
!  Misc. variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,iscr_alb
  INTEGER :: ibeg,iend,jbeg,jend
  INTEGER :: n_missing_albedo
  REAL :: albedo_to_cf,cloud_frac_vis
  REAL :: frac,term1,term2
  INTEGER :: istatus,iscr_frac_sat
  LOGICAL :: verbose = .FALSE.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (lvldbg > 90) verbose = .TRUE.
!
!-----------------------------------------------------------------------
!
!  Initialize all histogram arrays.
!
!-----------------------------------------------------------------------
!
  DO i=-10,20
    ihist_alb(i) = 0
    ihist_frac_sat(i) = 0
    ih_cf_sat(i) = 0
  END DO
!
!-----------------------------------------------------------------------
!
!  Read the interpolated satellite albedo data.
!  These data were created by laps2arps.f (based on ter2arps)
!  by interpolating OLAPS (vortex95) grid satellite data onto
!  the ARPS grid.
!
!-----------------------------------------------------------------------
!
  ibeg=MAX(1,nx/2-5)
  iend=MIN(nx,nx/2+5)
  jbeg=MAX(1,ny/2-5)
  jend=MIN(ny,ny/2+5)
!
  IF (verbose) THEN
    WRITE(6,'(a)') ' Satellite albedo value examples in the mid-domain:'

    WRITE(6,'(/4x,11(2x,i2,2x))') (i,i=ibeg,iend)
    WRITE(6,'(1x,i3,11f6.2)')                                           &
       (j,(albedo(i,j),i=ibeg,iend),j=jend,jbeg,-1)
  END IF
!
  n_missing_albedo = 0
!
!-----------------------------------------------------------------------
!
!  Horizontal array loop
!
!-----------------------------------------------------------------------
!
  DO j = 1,ny
    DO i = 1,nx
!
!-----------------------------------------------------------------------
!
!  We now only use the VIS data if the solar alt exceeds 15 deg
!
!-----------------------------------------------------------------------
!
      IF (solar_alt(i,j) > 15.0 .AND. albedo(i,j) > 0.) THEN
!
!-----------------------------------------------------------------------
!
!  Translate the albedo into cloud fraction
!  Store histogram information for satellite data
!
!-----------------------------------------------------------------------
!
        iscr_alb  = nint(albedo(i,j)*10.)
        iscr_alb  = MIN(MAX(iscr_alb,-10),20)
        ihist_alb(iscr_alb) = ihist_alb(iscr_alb) + 1

        cloud_frac_vis = albedo_to_cf(albedo(i,j))
!
!-----------------------------------------------------------------------
!
!  Fudge the frac at low solar elevation angles
!  Note the ramp extrapolates down to 9 deg to account for slight
!  errors in determining the solar elevation
!
!-----------------------------------------------------------------------
!
        IF(solar_alt(i,j) < 20. .AND. solar_alt(i,j) >= 9.) THEN
          frac = (20. - solar_alt(i,j)) / 10.
          term1 = .13 * frac
          term2 = 1. + term1
          cloud_frac_vis = (cloud_frac_vis + term1) * term2
        END IF

        iscr_frac_sat = nint(cloud_frac_vis*10.)
        iscr_frac_sat = MIN(MAX(iscr_frac_sat,-10),20)
        ihist_frac_sat(iscr_frac_sat) =                                 &
                   ihist_frac_sat(iscr_frac_sat) + 1
!
!-----------------------------------------------------------------------
!
!  Make sure satellite cloud fraction is between 0 and 1
!
!-----------------------------------------------------------------------
!
        cloud_frac_vis=max(min(cloud_frac_vis,1.0),0.0)
        cloud_frac_vis_a(i,j) = cloud_frac_vis

        iscr_frac_sat = nint(cloud_frac_vis*10.)
        iscr_frac_sat = MIN(MAX(iscr_frac_sat,-10),20)
        ih_cf_sat(iscr_frac_sat) =                                      &
                   ih_cf_sat(iscr_frac_sat) + 1

      ELSE   ! albedo(i,j) = r_missing
        n_missing_albedo =  n_missing_albedo + 1
        cloud_frac_vis_a(i,j) = r_missing
      END IF  ! albedo(i,j).ne.r_missing

    END DO ! i
  END DO ! j

  IF(n_missing_albedo == nx*ny) THEN ! Return with status = 0
    IF (myproc == 0) WRITE(6,'(a)')' All albedos were missing, return from process_vis'
    istatus = 0
    RETURN
  END IF

  IF (verbose) THEN
    CALL mptotali(n_missing_albedo)
    IF (myproc == 0) THEN

      WRITE(6,'(/a,i8/)')' N_MISSING_ALBEDO =',n_missing_albedo

      WRITE(6,'(a)')'       HISTOGRAMS'
      WRITE(6,'(a)')'  I   Alb  CFsat  CFsato'
      DO i = -5,15
        WRITE(6,'(i4,i5,i7,i8)')                                        &
          i,ihist_alb(i),ihist_frac_sat(i),ih_cf_sat(i)
      END DO ! i

      WRITE(6,'(a)') ' ==== process_vis: Albedo derived cld cvr'
      WRITE(6,'(/4x,11(2x,i2,2x))') (i,i=ibeg,iend)
      WRITE(6,'(1X,i3,11f6.2)')                                         &
           (j,(cloud_frac_vis_a(i,j),i=ibeg,iend),j=jend,jbeg,-1)
    END IF
  END IF
  istatus = 1

  RETURN
END SUBROUTINE process_vis

!
!##################################################################
!##################################################################
!######                                                      ######
!######                FUNCTION ALBEDO_TO_CF                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

FUNCTION albedo_to_cf( albedo )
!
!-----------------------------------------------------------------------
!
!  This version assumes 0.15 land albedo and 0.80 total cloud albedo.
!  Albedo fields are now calibrated before this call using the new
!  visible calibration files.
!
!  cloud_frac_vis = (albedo - 0.15) / (0.80 - 0.15)
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  REAL, INTENT(IN) :: albedo
  REAL :: albedo_to_cf
!
! Local Variables
!
  REAL, PARAMETER :: land_albedo=0.15
  REAL, PARAMETER :: max_albedo=0.80
  REAL, PARAMETER :: alb_denom=(1./(max_albedo-land_albedo))
  REAL :: cf_vis

  cf_vis = (albedo-land_albedo)*alb_denom

  albedo_to_cf = MIN(MAX(cf_vis,0.),1.0)

  RETURN
END FUNCTION albedo_to_cf
!
!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE SPREAD2                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE spread2 (cld_snd,wt_snd,i_snd,j_snd,n_cld_snd                &
           ,max_cld_snd,nz,iadas,jadas,k,cover,wt)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  This subroutine inserts the cloud sounding into the analysis arrays
!  at one point location.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!  07/95
!
!  MODIFICATION HISTORY:
!  03/20/96. (Jian Zhang)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: nz           ! number of cloud analysis levels
  INTEGER :: max_cld_snd  ! max. possible # of cloud soundings
                        ! in the domain

  INTEGER :: n_cld_snd    ! sequential index of cloud sounding sta.

  INTEGER :: iadas        ! i-loc of cld sounding stn in ADAS grid
  INTEGER :: jadas        ! j-loc of cld sounding stn in ADAS grid
  INTEGER :: k            ! the level where obs. cloud locates

  INTEGER :: i_snd (max_cld_snd)    ! array for i-location of
                                  ! cloud sounding stations
  INTEGER :: j_snd (max_cld_snd)    ! array for j-location of
                                  ! cloud sounding stations

  REAL :: cld_snd(max_cld_snd,nz)   ! cloud cover sounding array
  REAL :: wt_snd(max_cld_snd,nz)    ! weight array for cloud sounding

  REAL :: cover                     ! obs. cloud coverage
  REAL :: wt                        ! weight for the cloud cover obs.

  cld_snd(n_cld_snd,k) = cover
  wt_snd(n_cld_snd,k) = wt
  i_snd(n_cld_snd) = iadas
  j_snd(n_cld_snd) = jadas

  RETURN
END SUBROUTINE spread2
!
!##################################################################
!##################################################################
!######                                                      ######
!######              FUNCTION T_GROUND_K                     ######
!######                                                      ######
!##################################################################
!##################################################################
!

  FUNCTION t_ground_k (t_sfc_k,solar_alt,solar_ha                       &
                ,solar_dec,rlat,cvr_snow,r_missing_data,i,j,nx,ny)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate the ground skin temperature as a function of
!  the surface air temperature and solar positions.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:  (Jian Zhang)
!          Based on the LAPS cloud analysis code of 09/95
!
!  MODIFICATION HISTORY:
!  04/29/96  J. Zhang
!            Added ARPS format documents.
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!  Variable Declarition
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INCLUDE 'mp.inc'
!
!  INPUT:
  INTEGER :: i,j,nx,ny
  REAL :: t_sfc_k            ! surface air temperature
  REAL :: solar_alt          ! solar altitude angle
  REAL :: solar_ha           ! solar hour angle
  REAL :: solar_dec          ! solar declination angle
  REAL :: rlat               ! latitude of the grid point
  REAL :: cvr_snow           ! ground snow cover at the grid point
  REAL :: r_missing_data     ! flag for data hole
!
!  OUTPUT:
  REAL :: t_ground_k         ! ground skin temperature
!
!  LOCAL:
  INTEGER :: imid,jmid
  REAL :: high_alt,low_alt,corr_low,corr_high,corr_ramp
  REAL :: t_ground_fullcorr,t_snow_k,corr_negonly
  REAL :: solar_transit_alt,corr
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  imid=(nx/2)+1
  jmid=(ny/2)+1
  solar_transit_alt = 90. - ABS(rlat - solar_dec)

  IF(solar_ha < 0.) THEN ! morning
    high_alt = solar_transit_alt * .66
    low_alt  = solar_transit_alt * .36
!       low_alt  = high_alt - 11.
  ELSE                     ! afternoon
    high_alt = solar_transit_alt * .79
    low_alt  = solar_transit_alt * .49
!       low_alt  = high_alt - 11.
  END IF

  corr_high = 4.
  corr_low = -4.
  corr_ramp = (corr_high - corr_low) / (high_alt - low_alt)
!
!-----------------------------------------------------------------------
!
!  Warmer at day, colder at night
!
!-----------------------------------------------------------------------
!
  IF(solar_alt < low_alt)THEN
    corr = corr_low

  ELSE IF(solar_alt > high_alt)THEN
    corr = corr_high

  ELSE ! ramp the function
    corr = corr_low + (solar_alt - low_alt) * corr_ramp

  END IF

  IF(i == imid .AND. j == jmid .AND. myproc == 0) THEN
    WRITE(6,'(1x,a,2i4,a)') 'Sample solar position parms for point (',  &
          i,j,')'
    WRITE(6,'(1x,a,a)') '  hangle  altitude transalt  ratio  ',         &
                   ' highalt   lowalt  correctn'
    WRITE(6,'(1X,7F9.2)') solar_ha,solar_alt,solar_transit_alt,         &
          solar_alt/solar_transit_alt,high_alt,low_alt,corr
  END IF

  corr_negonly = MIN(corr,0.0)            ! Only colder at night

  IF(cvr_snow /= r_missing_data) THEN
    t_ground_fullcorr = t_sfc_k + corr  ! Warmer at day, colder at night
    t_snow_k = t_sfc_k + corr_negonly   ! Only colder at night
    t_snow_k = MIN(t_snow_k,273.15)     ! Snow can't be above 0C
    t_ground_k = t_ground_fullcorr * (1.-cvr_snow)                      &
                          + t_snow_k * cvr_snow
  ELSE
    t_ground_k = t_sfc_k + corr_negonly ! Only colder at night
  END IF

  RETURN
  END FUNCTION t_ground_k
!
!##################################################################
!##################################################################
!######                                                      ######
!######                   SUBROUTINE FILL1RAD                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE fill1rad_old(nx,ny,nz,nxlg,nylg,radfname,lvldbg,xs,ys,zs,        &
                  strhopt_rad,mapproj_rad,dx_rad,dy_rad,dz_rad,         &
                  dzmin_rad,ctrlat_rad,ctrlon_rad,                      &
                  tlat1_rad,tlat2_rad,tlon_rad,sclfct_rad,              &
                  isrcrad,stnrad,latrad,lonrad,elvrad,                  &
                  gridvel,gridref,gridnyq,gridtim,tem1,                 &
                  istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Reads radar data remapped on the ARPS grid.
!  This routine requires the remapping to occur on the same grid
!  as the analysis.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:  Jian Zhang
!  05/1996  Read the remapped radar data which was written by the
!           corresponding output routine "wrtrad" in remaplib.f.
!
!  MODIFICATION HISTORY:
!  03/19/97  J. Zhang
!            Added a line of error message when there is trouble
!            reading a radar file.
!  04/03/97  J. Zhang
!            Added the option of reading the data file created
!            from "WTRADCOL".  Added output for the remapping
!            parameters in the radar file (e.g., strhopt,mapproj,
!            dx,dy,dz,dzmin,ctrlat,ctrlon,tlat1,tlat2,tlon,scale)
!  04/07/97  J. Zhang
!            Added  the QC for the case when i,j,k outside the model
!            domain
!  04/09/97  J. Zhang
!            Added the Initializations for gridref, girdvel...
!  04/11/97  J. Zhang
!            Include dims.inc for nx,ny,nz
!  04/14/97  J. Zhang
!            Added message output for the case when actual # of
!            radar files exceeds the maximum allowed number in the
!            ADAS include file.  When that happens, the program will
!            stop.
!  03/31/98  J. Zhang
!            Deleted the option for reading the radar data file
!            created from "WRTRAD".
!  11/01/01  K. Brewster
!            Modified readrad to read1rad to read only one radar file
!            at a time to save memory.
!  04/23/07  K. Brewster
!            Added MPI handling code patterned after read1rad.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nxlg     Number of grid points in the x-direction (large grid/MPI)
!    nylg     Number of grid points in the y-direction (large grid/MPI)
!    radfname   file name for radar datasets
!    lvldbg   Debug level
!
!  OUTPUT:
!
!    isrcrad  index of radar source
!    stnrad   radar site name    character*4
!    latrad   latitude of radar  (degrees N)
!    lonrad   longitude of radar (degrees E)
!    elvrad   elevation of feed horn of radar (m MSL)
!
!    gridvel  radial velocity on ARPS grid
!    gridref  reflectivity on ARPS grid
!    gridnyq  nyquist velocity on ARPS grid
!    gridtim  observation time at ARPS grid
!    temrad   temporary radar array
!
!    istatus  status indicator
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INCLUDE 'mp.inc'
!
  INTEGER :: nx,ny,nz    ! the ARPS local grid size
  INTEGER :: nxlg,nylg   ! total domain size
  REAL :: xs(nx)
  REAL :: ys(ny)
  REAL :: zs(nx,ny,nz)
  CHARACTER (LEN=132) :: radfname
  INTEGER :: lvldbg
  INTEGER :: strhopt_rad,mapproj_rad
  REAL    :: dx_rad,dy_rad,dz_rad,dzmin_rad
  REAL    :: ctrlat_rad,ctrlon_rad,tlat1_rad,tlat2_rad,tlon_rad
  REAL    :: sclfct_rad
!
!  OUTPUT:  Radar site variables
!
  INTEGER :: isrcrad
  CHARACTER (LEN=5) :: stnrad
  REAL :: latrad
  REAL :: lonrad
  REAL :: elvrad
!
!  OUTPUT:  ARPS radar arrays
!
  REAL :: gridvel(nx,ny,nz)
  REAL :: gridref(nx,ny,nz)
  REAL :: gridnyq(nx,ny,nz)
  REAL :: gridtim(nx,ny,nz)
  REAL :: tem1(nxlg,nylg,nz)
  INTEGER :: istatus
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  REAL :: readk(nz)
  REAL :: readhgt(nz)
  REAL :: readref(nz)
  REAL :: readvel(nz)
  REAL :: readnyq(nz)
  REAL :: readtim(nz)
!
  INTEGER, ALLOCATABLE :: kntref(:)
  INTEGER, ALLOCATABLE :: kntvel(:)
!
  INTEGER :: mxradvr,nradvr
  PARAMETER(mxradvr=10)
  INTEGER :: iradvr(mxradvr)

  CHARACTER (LEN=4) :: stn
  CHARACTER (LEN=80) :: runname
  INTEGER :: ireftim,itime,vcpnum,idummy
  INTEGER :: iradfmt,strhoptin,mapprin
  INTEGER :: nchanl,ierr,iostatus
  INTEGER :: iyr, imon, idy, ihr, imin, isec
  INTEGER :: i,j,k,irad,jrad,krad,kk,ipt,klev

  INTEGER :: irngmin,irngmax
  INTEGER :: ilen,jlen,istart,iend,jstart,jend
  REAL :: height,rmin,rmax,rmax2,xrad,yrad,dist,dist2,elev,range
  REAL :: refelvmin,refelvmax,refrngmin,refrngmax

  REAL :: xrd,yrd,gridlat,gridlon,rdummy

  LOGICAL :: verbose = .FALSE.
  INTEGER :: ips, ipe, jps, jpe
!
  REAL, PARAMETER :: defelvmin = 0.5
  REAL, PARAMETER :: defelvmax = 19.5
  REAL, PARAMETER :: defrmin = 10.E03
  REAL, PARAMETER :: defrmax = 230.E03
!
!-----------------------------------------------------------------------
!
! Temporary arrays
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: tem2(:,:,:)
  REAL, ALLOCATABLE :: tem3(:,:,:)
  REAL, ALLOCATABLE :: tem4(:,:,:)

!-----------------------------------------------------------------------
!
! hdf arrays
!
!-----------------------------------------------------------------------
!
  INTEGER (KIND=selected_int_kind(4)), ALLOCATABLE :: itmp(:,:,:)
  REAL, ALLOCATABLE :: hmax(:), hmin(:) ! Temporary array
  INTEGER :: lens,dmpfmt,sd_id,isource
!
!-----------------------------------------------------------------------
!
! Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'grid.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  istatus=0
  refelvmin=-999.
  refelvmax=-999.
  refrngmin=-999.
  refrngmax=-999.
!
!-----------------------------------------------------------------------
!
!  Initializations
!
!-----------------------------------------------------------------------
!
  ALLOCATE(kntref(nz))
  ALLOCATE(kntvel(nz))

  IF (lvldbg > 90) THEN
    verbose = .TRUE.

    IF (mp_opt > 0) THEN
      ips = 2           ! Patch start index
      ipe = nx-2        ! Patch end index
      jps = 2
      jpe = ny-2
      IF (loc_x == 1)       ips = 1    ! To cover the same domain as non-MPI runs
      IF (loc_x == nproc_x) ipe = nx
      IF (loc_y == 1)       jps = 1
      IF (loc_y == nproc_y) jpe = ny
    ELSE
      ips = 1             ! Domain start index
      ipe = nx            ! Domain end index
      jps = 1
      jpe = ny
    END IF

  END IF

!
! OpenMP changed loop order to j,k,i:
!$OMP PARALLEL DO PRIVATE(j,k,i)
  DO j=1,ny
    DO k=1,nz
      DO i=1,nx
        gridref(i,j,k)=-9999.
        gridvel(i,j,k)=-9999.
        gridnyq(i,j,k)=-9999.
        gridtim(i,j,k)=-9999.
      END DO
    END DO
  END DO
!
  DO k=1,nz
    kntref(k)=0
    kntvel(k)=0
  END DO
!
!-----------------------------------------------------------------------
!
! Open file
!
!-----------------------------------------------------------------------
!
  CALL asnctl ('NEWLOCAL', 1, ierr)
  CALL asnfile(radfname, '-F f77 -N ieee', ierr)

  lens=LEN(trim(radfname))
  IF(radfname(lens-3:lens)=='hdf4')THEN
    dmpfmt=2
  ELSE
    dmpfmt=1
  ENDIF
  IF(myproc == 0) print*,'radfname=',TRIM(radfname),'dmpfmt=',dmpfmt

  IF( dmpfmt == 1 ) THEN
    ALLOCATE (tem2(nxlg,nylg,nz),stat=istatus)
    CALL check_alloc_status(istatus,"read1rad:tem2")
    ALLOCATE (tem3(nxlg,nylg,nz),stat=istatus)
    CALL check_alloc_status(istatus,"read1rad:tem3")
    ALLOCATE (tem4(nxlg,nylg,nz),stat=istatus)
    CALL check_alloc_status(istatus,"read1rad:tem4")

    IF(myproc == 0) THEN

      CALL getunit( nchanl )
      OPEN(UNIT=nchanl,FILE=trim(radfname),iostat=iostatus,            &
         FORM='unformatted',STATUS='old')
      IF(iostatus /= 0) THEN
        WRITE(6,'(a,a)') 'Error opening the radar file:',              &
          TRIM(radfname)
        istatus=-1
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!
!  Read radar description variables
!
!-----------------------------------------------------------------------
!
!      istatus=1
!
      READ(nchanl) stn
      stnrad=stn
      READ(nchanl) ireftim,itime,vcpnum,isrcrad,idummy,                &
                   idummy,idummy,idummy,idummy,idummy
!
      CALL abss2ctim(itime, iyr, imon, idy, ihr, imin, isec )
      iyr=MOD(iyr,100)
      WRITE(6,'(/a,i2.2,a,i2.2,a,i2.2,1X,i2.2,a,i2.2,a)')              &
         'Reading remapped raw radar data for: ',                      &
         imon,'/',idy,'/',iyr,ihr,':',imin,' UTC'
!
      READ(nchanl) runname
      READ(nchanl) iradfmt,strhopt_rad,mapproj_rad,irngmin,irngmax,    &
                   idummy,idummy,idummy,idummy,idummy

      READ(nchanl) dx_rad,dy_rad,dz_rad,dzmin_rad,ctrlat_rad,          &
           ctrlon_rad,tlat1_rad,tlat2_rad,tlon_rad,sclfct_rad,         &
           latrad,lonrad,elvrad,refelvmin,refelvmax
!
      READ(nchanl) nradvr,iradvr
!
! For NEXRAD radar source determine radar coverage area
! Use default values if these were not written by 88d2arps (old version)
! Fill radar coverage area with ref=0 as a default "no-echo" region.
!
      IF( isrcrad == 1 .OR. isrcrad == 2) THEN
        IF( refelvmin < 0.) refelvmin = defelvmin
        IF( refelvmax < refelvmin) refelvmax = defelvmax
        rmin=float(irngmin)
        rmax=float(irngmax)
        IF( rmin < 0. ) rmin = defrmin
        IF( rmax <= rmin ) rmax = defrmax
        WRITE(6,'(a,f9.1,a,f9.1,a)')                                   &
        ' Radar range limits: ',(0.001*rmin),' to ',(0.001*rmax),' km'
        WRITE(6,'(a,f9.2,a,f9.2,a)')                                   &
        ' Radar eleva limits: ',refelvmin,' to ',refelvmax,' degrees'
        rmax2=rmax*rmax
!
        CALL lltoxy(1,1,latrad,lonrad,xrad,yrad)
        irad=1+INT((xrad-xs(1))/dx)
        jrad=1+INT((yrad-ys(1))/dy)
        ilen=1+INT(rmax/dx)
        jlen=1+INT(rmax/dy)
        istart=max(1,(irad-ilen))
        iend=min(nx,(irad+ilen))
        jstart=max(1,(jrad-jlen))
        jend=min(ny,(jrad+jlen))
        DO j=jstart,jend
          DO i=istart,iend
            dist2=(xs(i)-xrad)*(xs(i)-xrad)+(ys(j)-yrad)*(ys(j)-yrad)
            IF(dist2 <= rmax2) THEN
              dist=sqrt(dist2)
              DO k=2,nz-1
                height=zs(i,j,k)-elvrad
                CALL beamelv(height,dist,elev,range)
                IF(elev > refelvmin .AND. elev < refelvmax .AND.         &
                 range > rmin .AND. range < rmax) gridref(i,j,k)=0.0
              END DO
            END IF
          END DO
        END DO
      END IF
!
      DO ipt=1,(nxlg*nylg)

        READ(nchanl,iostat=iostatus) i,j,xrd,yrd,                          &
                   gridlat,gridlon,elev,klev
        IF(iostatus /= 0) EXIT
        READ(nchanl,iostat=iostatus) (readk(kk),kk=1,klev)
        IF(iostatus /= 0) EXIT
        READ(nchanl,iostat=iostatus) (readhgt(kk),kk=1,klev)
        IF(iostatus /= 0) EXIT
        READ(nchanl,iostat=iostatus) (readref(kk),kk=1,klev)
        IF(iostatus /= 0) EXIT
        READ(nchanl,iostat=iostatus) (readvel(kk),kk=1,klev)
        IF(iostatus /= 0) EXIT
        READ(nchanl,iostat=iostatus) (readnyq(kk),kk=1,klev)
        IF(iostatus /= 0) EXIT
        READ(nchanl,iostat=iostatus) (readtim(kk),kk=1,klev)
        IF(iostatus /= 0) EXIT

        IF(i <= nxlg .AND. i >= 1 .AND. j <= nylg .AND. j >= 1) THEN
          DO kk=1,klev
            k=nint(readk(kk))
            IF(k <= nz.AND.k >= 1) THEN
              gridref(i,j,k)=max(readref(kk),gridref(i,j,k))
              gridvel(i,j,k)=readvel(kk)
              gridnyq(i,j,k)=readnyq(kk)
              gridtim(i,j,k)=readtim(kk)
              IF (gridref(i,j,k) > -200. .AND. gridref(i,j,k) < 200.)    &
                kntref(k)=kntref(k)+1
              IF (gridvel(i,j,k) > -200. .AND. gridvel(i,j,k) < 200.)    &
                kntvel(k)=kntvel(k)+1
            END IF  ! 1 < k < nz
          END DO  ! kk = 1, klev
        END IF  ! 1 < i < nxlg  & 1 < j < nylg

      END DO  ! ipt = 1, nx*ny

      ipt=ipt-1
      WRITE(6,'(a,i6,a)') ' End of file reached after reading',          &
                       ipt,' columns'
!
!-----------------------------------------------------------------------
!
!  Write statistics
!
!-----------------------------------------------------------------------
!
      WRITE(6,'(a)') ' Finished reading binary radar file'
      CLOSE(nchanl)
      CALL retunit( nchanl )

    END IF

    CALL mpupdatei(istatus,1)
    IF (istatus /= 0) RETURN

    CALL mpupdatec(stn,4)
    stnrad=stn
    CALL mpupdatei(strhopt_rad,1)
    CALL mpupdatei(mapproj_rad,1)
    CALL mpupdater(dx_rad,1)
    CALL mpupdater(dy_rad,1)
    CALL mpupdater(dz_rad,1)
    CALL mpupdater(dzmin_rad,1)
    CALL mpupdater(ctrlat_rad,1)
    CALL mpupdater(ctrlon_rad,1)
    CALL mpupdater(tlat1_rad,1)
    CALL mpupdater(tlat2_rad,1)
    CALL mpupdater(tlon_rad,1)
    CALL mpupdater(sclfct_rad,1)
    CALL mpupdater(latrad,1)
    CALL mpupdater(lonrad,1)
    CALL mpupdater(elvrad,1)

    CALL mpisplit3d(tem1,nx,ny,nz,gridref)
    CALL mpisplit3d(tem2,nx,ny,nz,gridvel)
    CALL mpisplit3d(tem3,nx,ny,nz,gridnyq)
    CALL mpisplit3d(tem4,nx,ny,nz,gridtim)

    DEALLOCATE(tem2,tem3,tem4)

  ELSE !HDF4 file

    IF (myproc == 0) THEN
      ALLOCATE (itmp(nxlg,nylg,nz),stat=istatus)
      CALL check_alloc_status(istatus,"read1rad:itmp")
      ALLOCATE (hmax(nz),stat=istatus)
      CALL check_alloc_status(istatus,"read1rad:hmax")
      ALLOCATE (hmin(nz),stat=istatus)
      CALL check_alloc_status(istatus,"read1rad:hmin")

      CALL hdfopen(trim(radfname), 1, sd_id)
      IF (sd_id < 0) THEN
        WRITE (6,*) "FILL1RAD: ERROR opening ",                        &
                trim(radfname)," for reading."
        istatus = -1
      ELSE
        istatus = 0
      END IF
    END IF
    CALL mpupdatei(istatus,1)

    IF (istatus == 0) THEN   ! Open file sucessfully

      IF (myproc == 0) THEN

        CALL hdfrdc(sd_id, 4, 'radid', stn, istatus)
        stnrad=stn
        CALL hdfrdi(sd_id, 'ireftim', ireftim, istatus)
        CALL hdfrdi(sd_id, 'itime', itime, istatus)
        CALL hdfrdi(sd_id, 'vcpnum', vcpnum, istatus)
        CALL hdfrdi(sd_id, 'isource', isource, istatus)
        CALL hdfrdc(sd_id, 40, 'runname', runname, istatus)
        CALL hdfrdi(sd_id, 'iradfmt', iradfmt, istatus)
        CALL hdfrdi(sd_id, 'strhopt', strhopt_rad, istatus)
        CALL hdfrdi(sd_id, 'mapproj', mapproj_rad, istatus)
        CALL hdfrdr(sd_id, 'dx', dx_rad, istatus)
        CALL hdfrdr(sd_id, 'dy', dy_rad, istatus)
        CALL hdfrdr(sd_id, 'dz', dz_rad, istatus)
        CALL hdfrdr(sd_id, 'dzmin', dzmin_rad, istatus)
        CALL hdfrdr(sd_id, 'ctrlat', ctrlat_rad, istatus)
        CALL hdfrdr(sd_id, 'ctrlon', ctrlon_rad, istatus)
        CALL hdfrdr(sd_id, 'trulat1', tlat1_rad, istatus)
        CALL hdfrdr(sd_id, 'trulat2', tlat2_rad, istatus)
        CALL hdfrdr(sd_id, 'trulon', tlon_rad, istatus)
        CALL hdfrdr(sd_id, 'sclfct', sclfct_rad, istatus)
        CALL hdfrdr(sd_id, 'latrad', latrad, istatus)
        CALL hdfrdr(sd_id, 'lonrad', lonrad, istatus)
        CALL hdfrdr(sd_id, 'elvrad', elvrad, istatus)
        CALL hdfrdi(sd_id, 'irngmin', irngmin, istatus)
        CALL hdfrdi(sd_id, 'irngmax', irngmax, istatus)
        CALL hdfrdr(sd_id, 'refelvmin', refelvmin, istatus)
        CALL hdfrdr(sd_id, 'refelvmax', refelvmax, istatus)
        CALL hdfrdi(sd_id, 'nradvr', nradvr, istatus)
        CALL hdfrd1di(sd_id,'iradvr', mxradvr,iradvr,istatus)
        IF (verbose) PRINT *, ' Got nradvr,iradvr: ',nradvr,iradvr
        CALL abss2ctim(itime, iyr, imon, idy, ihr, imin, isec )
        iyr=MOD(iyr,100)
        WRITE(6,'(/a,i2.2,a,i2.2,a,i2.2,1X,i2.2,a,i2.2,a)')            &
           'Reading remapped raw radar data for: ',                    &
           imon,'/',idy,'/',iyr,ihr,':',imin,' UTC'
      END IF
      CALL mpupdatei(istatus,1)

      IF (istatus == 0 ) THEN  ! read header sucessfully
        CALL mpupdatei(isrcrad,1)
        CALL mpupdatec(stnrad,5)
        CALL mpupdatei(strhopt_rad,1)
        CALL mpupdatei(mapproj_rad,1)
        CALL mpupdater(dx_rad,1)
        CALL mpupdater(dy_rad,1)
        CALL mpupdater(dz_rad,1)
        CALL mpupdater(dzmin_rad,1)
        CALL mpupdater(ctrlat_rad,1)
        CALL mpupdater(ctrlon_rad,1)
        CALL mpupdater(tlat1_rad,1)
        CALL mpupdater(tlat2_rad,1)
        CALL mpupdater(tlon_rad,1)
        CALL mpupdater(sclfct_rad,1)
        CALL mpupdater(latrad,1)
        CALL mpupdater(lonrad,1)
        CALL mpupdater(elvrad,1)
        CALL mpupdatei(irngmin,1)
        CALL mpupdatei(irngmax,1)
        CALL mpupdater(refelvmin,1)
        CALL mpupdater(refelvmax,1)

        IF (myproc == 0) THEN
          CALL hdfrd3d(sd_id,"gridref",nxlg,nylg,nz,tem1,istatus,itmp,hmax,hmin)
        END IF
        CALL mpupdatei(istatus,1)
        IF (istatus /= 0) GO TO 115
        CALL mpisplit3d(tem1,nx,ny,nz,gridref)

        IF (myproc == 0) THEN
          CALL hdfrd3d(sd_id,"gridvel",nxlg,nylg,nz,tem1,istatus,itmp,hmax,hmin)
        END IF
        CALL mpupdatei(istatus,1)
        IF (istatus /= 0) GO TO 115
        CALL mpisplit3d(tem1,nx,ny,nz,gridvel)

        IF (myproc == 0) THEN
          CALL hdfrd3d(sd_id,"gridnyq",nxlg,nylg,nz,tem1,istatus,itmp,hmax,hmin)
        END IF
        CALL mpupdatei(istatus,1)
        IF (istatus /= 0) GO TO 115
        CALL mpisplit3d(tem1,nx,ny,nz,gridnyq)

        IF (myproc == 0) THEN
          CALL hdfrd3d(sd_id,"gridtim",nxlg,nylg,nz,tem1,istatus,itmp,hmax,hmin)
          CALL hdfclose(sd_id,istatus)
        END IF
        CALL mpupdatei(istatus,1)
        IF (istatus /= 0) GO TO 115
        CALL mpisplit3d(tem1,nx,ny,nz,gridtim)

        istatus=0
      END IF
    END IF

    IF (myproc == 0) DEALLOCATE(hmin,hmax,itmp)

  END IF

  IF (verbose .AND. istatus == 0 ) THEN

    DO j=jps,jpe
      DO k=1,nz
        DO i=ips,ipe
          IF (gridref(i,j,k) > -200. .AND. gridref(i,j,k) < 200.)     &
            kntref(k)=kntref(k)+1
          IF (gridvel(i,j,k) > -200. .AND. gridvel(i,j,k) < 200.)     &
            kntvel(k)=kntvel(k)+1
        END DO
      END DO
    END DO

    DO k=1,nz
      CALL mptotali(kntref(k))
      CALL mptotali(kntvel(k))
    END DO

    IF (myproc == 0) THEN
      WRITE(6,'(a)') '  k       n ref      n vel'
      DO k = 1,nz
        WRITE(6,'(i3,2i10)') k,kntref(k),kntvel(k)
      END DO
    END IF

  END IF

  RETURN
!
!  Destination for hdf read error
!

  115   CONTINUE
  WRITE(6,'(/a/)') ' Error reading data in HDFREAD'

  istatus=-11

  DEALLOCATE(kntref)
  DEALLOCATE(kntvel)

  RETURN

END SUBROUTINE fill1rad_old
!
!##################################################################
!##################################################################
!######                                                      ######
!######                   SUBROUTINE FILL1RAD                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE fill1rad(nx,ny,nz,radfname,lvldbg,xs,ys,zs,                  &
                  strhopt_rad,mapproj_rad,dx_rad,dy_rad,dz_rad,         &
                  dzmin_rad,ctrlat_rad,ctrlon_rad,                      &
                  tlat1_rad,tlat2_rad,tlon_rad,sclfct_rad,              &
                  isrcrad,stnrad,latrad,lonrad,elvrad,                  &
                  gridvel,gridref,gridnyq,gridtim,                      &
                  istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Reads radar data remapped on the ARPS grid.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:  Yunheng Wang 07/02/2007
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    radfname   file name for radar datasets
!    lvldbg   Debug level
!
!  OUTPUT:
!
!    isrcrad  index of radar source
!    stnrad   radar site name    character*4
!    latrad   latitude of radar  (degrees N)
!    lonrad   longitude of radar (degrees E)
!    elvrad   elevation of feed horn of radar (m MSL)
!
!    gridvel  radial velocity on ARPS grid
!    gridref  reflectivity on ARPS grid
!    gridnyq  nyquist velocity on ARPS grid
!    gridtim  observation time at ARPS grid
!    temrad   temporary radar array
!
!    istatus  status indicator
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INCLUDE 'mp.inc'
!
  INTEGER :: nx,ny,nz    ! the ARPS local grid size
  REAL :: xs(nx)
  REAL :: ys(ny)
  REAL :: zs(nx,ny,nz)
  CHARACTER (LEN=132) :: radfname
  INTEGER :: lvldbg
  INTEGER :: strhopt_rad,mapproj_rad
  REAL    :: dx_rad,dy_rad,dz_rad,dzmin_rad
  REAL    :: ctrlat_rad,ctrlon_rad,tlat1_rad,tlat2_rad,tlon_rad
  REAL    :: sclfct_rad
!
!  OUTPUT:  Radar site variables
!
  INTEGER :: isrcrad
  CHARACTER (LEN=5) :: stnrad
  REAL :: latrad
  REAL :: lonrad
  REAL :: elvrad
!
!  OUTPUT:  ARPS radar arrays
!
  REAL :: gridvel(nx,ny,nz)
  REAL :: gridref(nx,ny,nz)
  REAL :: gridnyq(nx,ny,nz)
  REAL :: gridtim(nx,ny,nz)
  INTEGER :: istatus
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  REAL :: readk(nz)
  REAL :: readhgt(nz)
  REAL :: readref(nz)
  REAL :: readvel(nz)
  REAL :: readnyq(nz)
  REAL :: readtim(nz)
!
  INTEGER, ALLOCATABLE :: kntref(:)
  INTEGER, ALLOCATABLE :: kntvel(:)
!
  INTEGER :: mxradvr,nradvr
  PARAMETER(mxradvr=10)
  INTEGER :: iradvr(mxradvr)

  CHARACTER (LEN=4) :: stn
  CHARACTER (LEN=80) :: runname
  INTEGER :: ireftim,itime,vcpnum,idummy
  INTEGER :: iradfmt,strhoptin,mapprin
  INTEGER :: nchanl,ierr,iostatus
  INTEGER :: iyr, imon, idy, ihr, imin, isec
  INTEGER :: i,j,k,irad,jrad,kradf,kk,ipt,klev

  INTEGER :: irngmin,irngmax
  INTEGER :: ilen,jlen,istart,iend,jstart,jend
  REAL :: height,rmin,rmax,rmax2,xrad,yrad,dist,dist2,elev,range
  REAL :: refelvmin,refelvmax,refrngmin,refrngmax

  REAL :: xrd,yrd,gridlat,gridlon,rdummy

  LOGICAL :: verbose = .FALSE.
  INTEGER :: ips, ipe, jps, jpe
!
  REAL, PARAMETER :: defelvmin = 0.5
  REAL, PARAMETER :: defelvmax = 19.5
  REAL, PARAMETER :: defrmin = 10.E03
  REAL, PARAMETER :: defrmax = 230.E03
!
!-----------------------------------------------------------------------
!
! Temporary arrays
!
!-----------------------------------------------------------------------
!
  INTEGER :: typelev
  REAL    :: xmin, xmax, ymin, ymax
  INTEGER :: numradcol, nummaxlev

  INTEGER, ALLOCATABLE :: coli(:), colj(:), colk(:,:)
  INTEGER, ALLOCATABLE :: numlev(:)
  REAL,    ALLOCATABLE :: collat(:), collon(:)
  REAL,    ALLOCATABLE :: radcolhgt(:,:), radcolref(:,:), radcolvel(:,:), &
                          radcolnyq(:,:), radcoltim(:,:)
  INTEGER :: kcol
  INTEGER :: iproc, jproc  ! not related to iproc/jproc used elsewhere

!-----------------------------------------------------------------------
!
! hdf arrays
!
!-----------------------------------------------------------------------
!
  INTEGER (KIND=selected_int_kind(4)), ALLOCATABLE :: itmp(:,:)
  INTEGER :: lens,dmpfmt,sd_id,isource
!
!-----------------------------------------------------------------------
!
! Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'grid.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  istatus=0
  refelvmin=-999.
  refelvmax=-999.
  refrngmin=-999.
  refrngmax=-999.
  numradcol=0
  nummaxlev=0
!
!-----------------------------------------------------------------------
!
!  Initializations
!
!-----------------------------------------------------------------------
!
  IF (lvldbg > 90) THEN
    verbose = .TRUE.

    IF (mp_opt > 0) THEN
      ips = 2           ! Patch start index
      ipe = nx-2        ! Patch end index
      jps = 2
      jpe = ny-2
      IF (loc_x == 1)       ips = 1    ! To cover the same domain as non-MPI runs
      IF (loc_x == nproc_x) ipe = nx
      IF (loc_y == 1)       jps = 1
      IF (loc_y == nproc_y) jpe = ny
    ELSE
      ips = 1             ! Domain start index
      ipe = nx            ! Domain end index
      jps = 1
      jpe = ny
    END IF

    ALLOCATE(kntref(nz))
    ALLOCATE(kntvel(nz))

    DO k=1,nz
      kntref(k)=0
      kntvel(k)=0
    END DO

  END IF

!
! OpenMP changed loop order to j,k,i:
!$OMP PARALLEL DO PRIVATE(j,k,i)
  DO j=1,ny
    DO k=1,nz
      DO i=1,nx
        gridref(i,j,k)=-999.
        gridvel(i,j,k)=-999.
        gridnyq(i,j,k)=-999.
        gridtim(i,j,k)=-999.
      END DO
    END DO
  END DO
!
!
!-----------------------------------------------------------------------
!
! Open file
!
!-----------------------------------------------------------------------
!
  CALL asnctl ('NEWLOCAL', 1, ierr)
  CALL asnfile(radfname, '-F f77 -N ieee', ierr)

  lens=LEN(trim(radfname))
  IF(radfname(lens-3:lens)=='hdf4')THEN
    dmpfmt=3
  ELSE
    dmpfmt=1
  ENDIF
  IF (myproc == 0) WRITE(6,'(1x,a,I5,3a,i2)') 'Processor: ',myproc,     &
                   ' reading ',TRIM(radfname),', dmpfmt = ',dmpfmt

  IF( dmpfmt == 1 ) THEN

    CALL getunit( nchanl )
    OPEN(UNIT=nchanl,FILE=trim(radfname),iostat=iostatus,               &
         FORM='unformatted',STATUS='old')
    IF(iostatus /= 0) THEN
      WRITE(6,'(a,a)') 'Error opening the radar file:',TRIM(radfname)
      istatus=-1
      GO TO 999
    END IF
!
!-----------------------------------------------------------------------
!
!  Read radar description variables
!
!-----------------------------------------------------------------------
!
    istatus=1

    READ(nchanl) stn
    stnrad=stn
    READ(nchanl) ireftim,itime,vcpnum,isrcrad,idummy,                   &
                 idummy,idummy,idummy,idummy,idummy

    READ(nchanl) runname
    READ(nchanl) iradfmt,strhopt_rad,mapproj_rad,irngmin,irngmax,       &
!                 idummy,idummy,idummy,idummy,idummy
                 typelev,numradcol,nummaxlev,idummy,idummy

    READ(nchanl) dx_rad,dy_rad,dz_rad,dzmin_rad,ctrlat_rad,             &
         ctrlon_rad,tlat1_rad,tlat2_rad,tlon_rad,sclfct_rad,            &
         latrad,lonrad,elvrad,refelvmin,refelvmax

    READ(nchanl) nradvr,iradvr

    READ(nchanl) xmin, xmax, ymin, ymax
  ELSE
!   IF (.NOT. ALLOCATED(itmp)) THEN
!     ALLOCATE (itmp(nx,ny,nz),stat=istatus)
!     CALL check_alloc_status(istatus, "rdradcol:itmp")
!   END IF

    CALL hdfopen(trim(radfname), 1, sd_id)
    IF (sd_id < 0) THEN
      IF (myproc == 0) WRITE (6,*) "FILL1RAD: ERROR opening ",          &
              trim(radfname)," for reading."
      istatus = -1
      GO TO 999
    ELSE
      istatus = 0
    END IF

    CALL hdfrdc(sd_id, 4, 'radid', stn, istatus)

    !
    !  Check for corrupt file.
    !
    IF (istatus /= 0 ) THEN
      WRITE(6,'(1x,a,1x,a)') 'FILL1RAD:  ERROR corrupt file: ',    &
                trim(radfname),' for reading.'
      CALL hdfclose(sd_id,istatus)
      istatus = -1
      GOTO 999
    END IF

    stnrad=stn
    CALL hdfrdi(sd_id, 'ireftim', ireftim, istatus)
    CALL hdfrdi(sd_id, 'itime', itime, istatus)
    CALL hdfrdi(sd_id, 'vcpnum', vcpnum, istatus)
    CALL hdfrdi(sd_id, 'isource', isrcrad, istatus)
    CALL hdfrdc(sd_id, 40, 'runname', runname, istatus)
    CALL hdfrdi(sd_id, 'iradfmt', iradfmt, istatus)
    CALL hdfrdi(sd_id, 'strhopt', strhopt_rad, istatus)
    CALL hdfrdi(sd_id, 'mapproj', mapproj_rad, istatus)
    CALL hdfrdr(sd_id, 'dx', dx_rad, istatus)
    CALL hdfrdr(sd_id, 'dy', dy_rad, istatus)
    CALL hdfrdr(sd_id, 'dz', dz_rad, istatus)
    CALL hdfrdr(sd_id, 'dzmin', dzmin_rad, istatus)
    CALL hdfrdr(sd_id, 'ctrlat', ctrlat_rad, istatus)
    CALL hdfrdr(sd_id, 'ctrlon', ctrlon_rad, istatus)
    CALL hdfrdr(sd_id, 'trulat1', tlat1_rad, istatus)
    CALL hdfrdr(sd_id, 'trulat2', tlat2_rad, istatus)
    CALL hdfrdr(sd_id, 'trulon', tlon_rad, istatus)
    CALL hdfrdr(sd_id, 'sclfct', sclfct_rad, istatus)
    CALL hdfrdr(sd_id, 'latrad', latrad, istatus)
    CALL hdfrdr(sd_id, 'lonrad', lonrad, istatus)
    CALL hdfrdr(sd_id, 'elvrad', elvrad, istatus)
    CALL hdfrdi(sd_id, 'irngmin', irngmin, istatus)
    CALL hdfrdi(sd_id, 'irngmax', irngmax, istatus)
    CALL hdfrdr(sd_id, 'refelvmin', refelvmin, istatus)
    CALL hdfrdr(sd_id, 'refelvmax', refelvmax, istatus)
    CALL hdfrdi(sd_id, 'nradvr', nradvr, istatus)
    CALL hdfrd1di(sd_id,'iradvr', mxradvr,iradvr,istatus)
    IF (verbose) PRINT *, ' Got nradvr,iradvr: ',nradvr,iradvr

    CALL hdfrdi(sd_id, 'typelev', typelev, istatus)

    CALL hdfrdi(sd_id, 'numradcol', numradcol, istatus)
    CALL hdfrdi(sd_id, 'nummaxelv', nummaxlev, istatus)

    CALL hdfrdr(sd_id, 'xmin',   xmin,    istatus)
    CALL hdfrdr(sd_id, 'xmax',   xmax,    istatus)
    CALL hdfrdr(sd_id, 'ymin',   xmin,    istatus)
    CALL hdfrdr(sd_id, 'ymax',   ymax,    istatus)
  END IF

  CALL abss2ctim(itime, iyr, imon, idy, ihr, imin, isec )
  iyr=MOD(iyr,100)
  IF (myproc == 0) WRITE(6,'(/a,i2.2,a,i2.2,a,i2.2,1X,i2.2,a,i2.2,a)') &
     'Reading remapped raw radar data for: ',                      &
     imon,'/',idy,'/',iyr,ihr,':',imin,' UTC'

  IF (typelev /= 1)  THEN
    WRITE(6,'(1x,a,/,1x,2a,/)')                                         &
      'Sub fill1rad was rewritten to support new radar data format. ',  &
      'In order to read data in old format, please change all calls of',&
      ' fill1rad into fill1rad_old in the source code.'
!   CALL arpsstop('Maybe in old data format?',1)
    istatus = -1
    GO TO 888
  END IF

!
!-----------------------------------------------------------------------
!
!  Make sure we don't have a corrupt file.
!
!-----------------------------------------------------------------------
!

  IF (numradcol < 1 .OR. nummaxlev < 1) THEN
    WRITE(6,'(//,a)') 'WARNING: Data corruption.  File not processed.'
    istatus = -1
    GO TO 888
  END IF

  ALLOCATE(coli(numradcol), STAT = istatus)
  CALL check_alloc_status(istatus, "fill1rad:coli")

  ALLOCATE(colj(numradcol), STAT = istatus)
  CALL check_alloc_status(istatus, "fill1rad:colj")

  ALLOCATE(colk(nummaxlev,numradcol), STAT = istatus)
  CALL check_alloc_status(istatus, "fill1rad:colk")

  ALLOCATE(numlev(numradcol), STAT = istatus)
  CALL check_alloc_status(istatus, "fill1rad:numlev")

  ALLOCATE(collat(numradcol), STAT = istatus)
  CALL check_alloc_status(istatus, "fill1rad:collat")

  ALLOCATE(collon(numradcol), STAT = istatus)
  CALL check_alloc_status(istatus, "fill1rad:collon")

  ALLOCATE(radcolhgt(nummaxlev,numradcol), STAT = istatus)
  CALL check_alloc_status(istatus, "fill1rad:radcolhgt")

  ALLOCATE(radcolref(nummaxlev,numradcol), STAT = istatus)
  CALL check_alloc_status(istatus, "fill1rad:radcolref")

  ALLOCATE(radcolvel(nummaxlev,numradcol), STAT = istatus)
  CALL check_alloc_status(istatus, "fill1rad:radcolvel")

  ALLOCATE(radcolnyq(nummaxlev,numradcol), STAT = istatus)
  CALL check_alloc_status(istatus, "fill1rad:radcolnyq")

  ALLOCATE(radcoltim(nummaxlev,numradcol), STAT = istatus)
  CALL check_alloc_status(istatus, "fill1rad:radcoltim")

  radcolhgt(:,:) = -999.
  radcolref(:,:) = -999.
  radcolvel(:,:) = -999.
  radcolnyq(:,:) = -999.
  radcoltim(:,:) = -999.

  IF (dmpfmt == 1) THEN ! READ binary files

    kcol  = 0
    DO kcol=1,numradcol
      READ(nchanl,END=201) i,j,xrd,yrd,                           &
                           latrad,lonrad,elev,klev

      coli(kcol) = i
      colj(kcol) = j

      collat(kcol)  = latrad
      collon(kcol)  = lonrad
      numlev(kcol)  = klev

      READ(nchanl,END=202)      (colk(kk,kcol),kk=1,klev)
      READ(nchanl,END=202) (radcolhgt(kk,kcol),kk=1,klev)
      READ(nchanl,END=202) (radcolref(kk,kcol),kk=1,klev)
      READ(nchanl,END=202) (radcolvel(kk,kcol),kk=1,klev)
      READ(nchanl,END=202) (radcolnyq(kk,kcol),kk=1,klev)
      READ(nchanl,END=202) (radcoltim(kk,kcol),kk=1,klev)
    END DO

    201   CONTINUE
    WRITE(6,'(a,i6,a)') ' End of file reached after reading',     &
                        kcol,' columns'
    GO TO 205

    202   CONTINUE
    WRITE(6,'(a,i6,a)') ' End of file reached while reading',     &
                        kcol,' column'

    205   CONTINUE

!   CLOSE(nchanl)
!   CALL retunit( nchanl )

    nummaxlev = MAX(nummaxlev,klev)

  ELSE                  ! READ hdf file

    CALL hdfrd1d (sd_id,'radcollat',numradcol,collat,istatus)
    CALL hdfrd1d (sd_id,'radcollon',numradcol,collon,istatus)
    CALL hdfrd1di(sd_id,'numelev',  numradcol,numlev,istatus)

    CALL hdfrd1di(sd_id,'radcoli',  numradcol,coli,istatus)
    CALL hdfrd1di(sd_id,'radcolj',  numradcol,colj,istatus)
    CALL hdfrd2di(sd_id,'radcolk',  nummaxlev,numradcol,colk,istatus)

    IF (.NOT. ALLOCATED(itmp)) THEN
      ALLOCATE (itmp(nummaxlev,numradcol),stat=istatus)
      CALL check_alloc_status(istatus, "rdradcol:itmp")
    END IF

    CALL hdfrd2d (sd_id,'radcolhgt',nummaxlev,numradcol,radcolhgt,istatus,itmp)
    CALL hdfrd2d (sd_id,'radcolref',nummaxlev,numradcol,radcolref,istatus,itmp)
    CALL hdfrd2d (sd_id,'radcolvel',nummaxlev,numradcol,radcolvel,istatus,itmp)
    CALL hdfrd2d (sd_id,'radcolnyq',nummaxlev,numradcol,radcolnyq,istatus,itmp)
    CALL hdfrd2d (sd_id,'radcoltim',nummaxlev,numradcol,radcoltim,istatus,itmp)

!   CALL hdfclose(sd_id,istatus)
  END IF

  DO kcol = 1, numradcol
     iproc = (coli(kcol)-2)/(nx-3) + 1
     jproc = (colj(kcol)-2)/(ny-3) + 1

     IF (loc_x == iproc .AND. loc_y == jproc) THEN
       i =  MOD((coli(kcol)-2),(nx-3)) + 2
       j =  MOD((colj(kcol)-2),(ny-3)) + 2
       klev = numlev(kcol)

       DO kk=1,klev
         k=colk(kk,kcol)
         IF(k <= nz.AND.k >= 1) THEN
           gridref(i,j,k)=max(radcolref(kk,kcol),gridref(i,j,k))
           gridvel(i,j,k)=radcolvel(kk,kcol)
           gridnyq(i,j,k)=radcolnyq(kk,kcol)
           gridtim(i,j,k)=radcoltim(kk,kcol)
         END IF  ! 1 < k < nz
       END DO  ! kk = 1, klev
    END IF  ! 1 < i < nxlg  & 1 < j < nylg
  END DO

  IF (verbose) THEN

    DO j=jps,jpe
      DO k=1,nz
        DO i=ips,ipe
          IF (gridref(i,j,k) > -200. .AND. gridref(i,j,k) < 200.)     &
            kntref(k)=kntref(k)+1
          IF (gridvel(i,j,k) > -200. .AND. gridvel(i,j,k) < 200.)     &
            kntvel(k)=kntvel(k)+1
        END DO
      END DO
    END DO

    DO k=1,nz
      CALL mptotali(kntref(k))
      CALL mptotali(kntvel(k))
    END DO

    IF (myproc == 0) THEN
      WRITE(6,'(a)') '  k       n ref      n vel'
      DO k = 1,nz
        WRITE(6,'(i3,2i10)') k,kntref(k),kntvel(k)
      END DO
    END IF

  END IF

!
! Label 999 is for opens that failed.
! Label 888 is for all successes and those failures after open works.
!

  888 CONTINUE

  IF ( dmpfmt == 1 ) THEN
    CLOSE(nchanl)
    CALL retunit( nchanl )
  ELSE
    CALL hdfclose(sd_id,istatus)
  END IF

  999 CONTINUE

  IF (verbose) THEN
    DEALLOCATE(kntref)
    DEALLOCATE(kntvel)
  END IF

  IF (ALLOCATED(itmp)) DEALLOCATE(itmp)

  IF (ALLOCATED(coli)) THEN
    DEALLOCATE(coli, colj, colk, numlev)
    DEALLOCATE(collat, collon)
    DEALLOCATE(radcolhgt, radcolref, radcolvel, radcolnyq, radcoltim)
  END IF

  RETURN
END SUBROUTINE fill1rad
