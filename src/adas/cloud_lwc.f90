!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE CLOUD_LWC                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE cloud_lwc (nx,ny,nz,time,dirname,runname                     &
           ,clouds_3d                                                   &
           ,temp_3d,rh_3d,p_3d,zs_3d                                    &
           ,istat_radar,radar_3d                                        &
           ,iflag_slwc,slwc_3d,cice_3d,ctmp_3d                          &
           ,l_flag_incld_w,w_3d                                         &
           ,l_mask_pcptype,l_flag_cldtyp,cldpcp_type_3d                 &
           ,l_flag_icing,icing_index_3d, tem1,                          &
           istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  This routine calculates SLWC, Cloud Type, and Icing Index.
!  This routine also does the in-cloud w-field. These have been
!  combined to make more efficient use of looping through large
!  arrays. The input flags can be adjusted to allow only some of
!  these to be returned.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:  Jian Zhang
!  02/96    Based on LAPS code by Steve Albers,    04/94
!
!  MODIFICATION HISTORY:
!
!  05/10/96  J. Zhang
!            Modified for ADAS format. Added full documentation.
!  08/21/96  J. Zhang
!            Added subroutine cloud_type_qc, a quality check
!            scheme to assure the consistency of cumulonimbus clouds.
!  04/11/97  J. Zhang
!            Included the dims.inc file for nx,ny,nz.
!  08/05/97  J. Zhang
!            Included the adascld25.inc file for L_Cb_qc,
!            r_missing, wmhr_cu, wmhr_sc, and wc_st.
!  09/11/97  J. Zhang
!            Change adascld25.inc to adascld26.inc.
!  11/02/97  J. Zhang
!            Change the input argument "slwc_3d(i,j,k)" (unit: g/kg)
!            in the call to the icing index algorithm subroutine
!            to "lwc_g_m3" (unit: g/m3).
!  11/18/97  J. Zhang
!            Added flag 'clddiag' in adascld26.inc for whether to
!            output the special cloud fields.
!  07/10/01  K. Brewster
!            Cloud temperature added as an output array for use in
!            optional latent heating adjustment to temperature.
!  07/10/01  K. Brewster
!            Removed radar-based depletion of ice.
!  07/10/01  K. Brewster
!            Added 1-d cloud temp array to get_sfm_1d call.
!            Cleaned up calculation of entrainment dilution to make
!            code more readable.
!  04/17/2007 Y. Wang
!            Changed the logic for diagnostic outputs based on the MPI
!            changes to ensure the same outputs between MPI and non-MPI
!            runs, and the order of i,j loops for efficiency.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INCLUDE 'adas.inc'
  INCLUDE 'mp.inc'
  INCLUDE 'bndry.inc'
!
!-----------------------------------------------------------------------
!
!  LWC Flags:
!  iflag_slwc =  1 : Adiabatic LWC
!  iflag_slwc =  2 : Adjusted  LWC
!  iflag_slwc =  3 : Adjusted  SLWC
!  iflag_slwc = 13 : New Smith-Feddes LWC
!
!  INCLUDE:  (from adas.inc)
!  clddiag           ! flag for whether to output the special cloud
!
!  INPUT:
  INTEGER :: nx,ny,nz             ! model grid size
  REAL :: time               ! analysis time in seconds from the model
                           ! initial time
  CHARACTER (LEN=*) :: dirname
  CHARACTER (LEN=*) :: runname

!  INCLUDE:  (from adas.inc)
!  real r_missing               ! indicator for missing data
!  real wmhr_Cu      ! coef_wmax for Cu-clouds (m/s/m)
!  real wmhr_Sc      ! coef_wmax for Sc-clouds (m/s/m)
!  real wc_St
!  real thresh_cvr  ! cld fractn threshld for "cloudy"

  REAL :: clouds_3d(nx,ny,nz)     ! 3D cloud fraction analysis
  REAL :: temp_3d(nx,ny,nz)       ! 3D temp.(K) field on model grid
  REAL :: p_3d(nx,ny,nz)          ! 3D pres.(pa) field on model grid
  REAL :: rh_3d(nx,ny,nz)         ! 3D relative humidity (%) on model grid
  REAL :: zs_3d(nx,ny,nz)         ! physical hgts (m) at model scalar g.p.
  REAL :: tem1(nx,ny,nz)

  INTEGER :: istat_radar          ! Valid 3D radar data? (1=yes)
  REAL :: radar_3d(nx,ny,nz)      ! 3D reflectivity interp. on model grid

  INTEGER :: iflag_slwc           ! indicator for CLWC schemes

  LOGICAL :: l_flag_incld_w       ! Analize w-field in clouds?
  LOGICAL :: l_flag_cldtyp        ! Analize cloud/precip types?
  LOGICAL :: l_flag_icing         ! Analyze icing severity?
!
!  OUTPUT:
  INTEGER :: istatus
  REAL :: slwc_3d(nx,ny,nz)       ! liquid water content field (g/kg)
  REAL :: cice_3d(nx,ny,nz)       ! ice water content (g/kg)
  REAL :: ctmp_3d(nx,ny,nz)       ! cloud temperature (K)
  REAL :: w_3d(nx,ny,nz)          ! vertical velocity (m/s) in clouds
  REAL :: w_tem(nx,ny,nz)          ! vertical velocity (m/s) in clouds
  INTEGER :: cldpcp_type_3d(nx,ny,nz) ! cloud/previp type field
  INTEGER :: icing_index_3d(nx,ny,nz) ! icing severity index
  LOGICAL :: l_mask_pcptype(nx,ny)    ! Analyze 2D precip type using
                                    ! simulated radar data ??
!  LOCAL:

  INTEGER :: ibase_array(nx,ny)   ! lowest cloud base hgt index
  INTEGER :: itop_array(nx,ny)    ! highest cloud top hgt index
!
!  LOCAL:
  REAL :: t_1d(nz)
  REAL :: zs_1d(nz)
  REAL :: zp_1d(nz)
  REAL :: p_mb_1d(nz)
  REAL :: p_pa_1d(nz)
  REAL :: slwc_1d(nz)
  REAL :: w_1d(nz)
  REAL :: cice_1d(nz)
  REAL :: ctmp_1d(nz)
  REAL :: dte_dz_1d(nz)
  INTEGER :: cldtyp_1d(nz)
  INTEGER :: istat_1
  INTEGER :: is_dup
  LOGICAL :: l_cb_qc   ! Check the spatial continuity of the Cb cloud?
  PARAMETER (l_cb_qc = .true.)
!
!  Temporary array:
!  REAL :: tem1(nz)
!
!  Constant:
  REAL :: rair
  PARAMETER (rair=287.04)             ! J/deg/kg
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,k1,kb,kt,mode,kr, k_highest,k_base,k_top
  INTEGER :: imid,jmid,ibeg,iend,kbeg,kend
  INTEGER :: n_cld_columns,n_slwc_call
  REAL    :: cld_base_m,cld_top_m,cld_base_qc_m,cld_top_qc_m
  REAL    :: xindex,ramp,rho,lwc_g_m3
  LOGICAL :: l_cloud, l_prt
  INTEGER :: itype,i_precip_type,i_cldtyp, iarg,istatwrt
  CHARACTER (LEN=2) :: c2_type

  CHARACTER (LEN=6)  :: varid
  CHARACTER (LEN=20) :: varname
  CHARACTER (LEN=20) :: varunits

  INTEGER :: idiag, jdiag, idiagproc, jdiagproc
  INTEGER :: ips, ipe, jps, jpe
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (myproc == 0) THEN
    WRITE(6,*)' Start LWC subroutine.'
    PRINT*,'wmhr_Cu=',wmhr_cu,' wmhr_Sc=',wmhr_sc,' wc_St=', wc_st
    PRINT*,'thresh_cvr=',thresh_cvr
  END IF
  istatus = 0
!
!-----------------------------------------------------------------------
!
!  Print out information about the data passed in.
!
!-----------------------------------------------------------------------
!
  imid=nx/2
  ibeg=MAX(1,(nx/2-5))
  iend=MIN(nx,(nx/2+5))
  jmid=ny/2
  kbeg=MAX(1,(nz/2-5))
  kend=MIN(nz,(nz/2+5))
  IF(clddiag == 1 .AND. myproc == 0) THEN
    WRITE(6,'(a)')'grid point heights (m) at j=jmid,i=imid-5,+5'
    DO k=kbeg,kend
      WRITE(6,622) k,(zs_3d(i,jmid,k),i=ibeg,iend)
      622       FORMAT (1X,'k=',i2,11F7.0)
    END DO
  END IF
!
!-----------------------------------------------------------------------
!
!  Check consistency of input flags
!
!-----------------------------------------------------------------------
!
  IF(l_flag_icing) THEN
    iflag_slwc = 13
  END IF

  IF(iflag_slwc >= 1) THEN
    l_flag_cldtyp = .true.
  END IF
!
!-----------------------------------------------------------------------
!
!  Generate Lowest Base and Highest Top Arrays
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) WRITE(6,*)' Generating Lowest Base and Highest Top Arrays'
  DO j = 1,ny
    DO i = 1,nx
      ibase_array(i,j) = nz
      itop_array(i,j) = 0
    END DO
  END DO

  DO j = 1,ny-1
    DO i = 1,nx-1
      DO k = nz-1,1,-1
        IF(clouds_3d(i,j,k) >= thresh_cvr) THEN
          ibase_array(i,j) = k
        END IF
      END DO
    END DO
  END DO

  DO j = 1,ny-1
    DO i = 1,nx-1
      DO k = 1,nz-1
        IF(clouds_3d(i,j,k) >= thresh_cvr) THEN
          itop_array(i,j) = k
        END IF
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Initialize liquid water content array
!
!-----------------------------------------------------------------------
!
  IF(iflag_slwc /= 0) THEN
    IF (myproc == 0) WRITE(6,'(/,2x,a)') 'Initializing SLWC array'
    DO k = 1,nz
      DO j = 1,ny
        DO i = 1,nx
          slwc_3d(i,j,k) = 0.0
          cice_3d(i,j,k) = 0.0
          ctmp_3d(i,j,k) = temp_3d(i,j,k)
        END DO
      END DO
    END DO
  END IF
!
!-----------------------------------------------------------------------
!
!  Initialize cloud vertical velicity array
!
!-----------------------------------------------------------------------
!
  IF(l_flag_incld_w)THEN
    IF (myproc == 0) WRITE(6,*)' Initializing w array'

    DO k = 1,nz
      DO j = 1,ny
        DO i = 1,nx
          w_3d(i,j,k) = r_missing
          w_tem(i,j,k) = 0.0
        END DO
      END DO
    END DO
  END IF
!
!-----------------------------------------------------------------------
!
!  Initialize cloud type array
!
!-----------------------------------------------------------------------
!
  IF(l_flag_cldtyp) THEN
    IF (myproc == 0) WRITE(6,*)' Initializing Cloud Type Array'
    DO k = 1,nz
      DO j = 1,ny
        DO i = 1,nx
          cldpcp_type_3d(i,j,k) = 0
        END DO
      END DO
    END DO
  END IF
!
  n_cld_columns = 0
  n_slwc_call = 0
!
!-----------------------------------------------------------------------
!
!  Find Cloud Layers and Computing Output Field(s)
!  The procedure works column by column.
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0)                                                     &
    WRITE(6,*)' Finding Cloud Layers and Computing Output Field(s)'

!
! Don't let duplicate (overlap) points get into the computations.
!
  IF (mp_opt > 0) THEN
    ips = 2             ! Patch index
    ipe = nx-2
    jps = 2
    jpe = ny-2
    IF (loc_x == 1)       ips = 1         ! to cover the same domain as non-MPI runs below
    IF (loc_x == nproc_x) ipe = nx-1
    IF (loc_y == 1)       jps = 1
    IF (loc_y == nproc_y) jpe = ny-1
  ELSE
    ips = 1             ! Memory index for the whole domain
    ipe = nx-1
    jps = 1
    jpe = ny-1
  END IF

  idiag = 14            ! diagnostic output point
  jdiag = 14

  idiagproc = (idiag-2) / (nx-3) + 1   ! We want to output 1 diagnostic output only
  jdiagproc = (jdiag-2) / (ny-3) + 1   ! and it should be the same point as non-MPI runs

  IF (loc_x == idiagproc) THEN
    idiag = MOD((idiag-2),(nx-3)) + 2   ! Local index for global middle point
  ELSE
    jdiag = -999
  END IF

  IF (loc_y == jdiagproc) THEN
    jdiag = MOD((jdiag-2),(ny-3)) + 2   ! Local index for global middle point
  ELSE
    jdiag = -999
  END IF

!  DO i = 1,nx-1
!    is_dup = 0
!    IF ((mp_opt > 0) .AND. (i >= nx-2) .AND. (loc_x .NE. nproc_x)) is_dup=1
!    DO j = 1,ny-1
!      IF ((mp_opt > 0) .AND. (j >= ny-2) .AND. (loc_y .NE. nproc_y)) is_dup=1

  DO j = jps, jpe
    DO i = ips, ipe

      l_prt=.false.
!      IF( i == idiag .and. j == jdiag) l_prt=.true.

      l_cloud = .false.
!
!-----------------------------------------------------------------------
!
!  Generate vertical sounding at this grid point and
!  call SLWC routine
!
!-----------------------------------------------------------------------
!
      IF( ibase_array(i,j) < nz) THEN ! At least one lyr exists

        l_cloud = .true.
        n_cld_columns = n_cld_columns + 1

        DO k = 1,nz-1                      ! Initialize
          t_1d(k) = temp_3d(i,j,k)
          zs_1d(k) = zs_3d(i,j,k)
          p_pa_1d(k) = p_3d(i,j,k)
          p_mb_1d(k) = 0.01*p_3d(i,j,k)
          cldtyp_1d(k) = 0
        END DO
        IF(l_prt) THEN
          WRITE(6,*)
          PRINT*,i,j
          WRITE(6,'(1X,a,33F7.1)') 'T:',(t_1d(k),k=2,34)
          WRITE(6,'(1X,a,33F7.1)') 'P:',(p_mb_1d(k),k=2,34)
          WRITE(6,'(1X,a,33F7.0)') 'z:',(zs_1d(k),k=2,34)
          WRITE(6,'(1X,a,33F6.1)') 'rh:',(rh_3d(i,j,k),k=2,34)
          WRITE(6,'(1X,a,33F6.2)') 'cldcvr:',(clouds_3d(i,j,k),k=2,34)
          WRITE(6,'(1X,a,33F6.2)') 'radar:',(radar_3d(i,j,k),k=2,34)
        END IF
!
!-----------------------------------------------------------------------
!
!  Get Base and Top
!
!-----------------------------------------------------------------------
!
        k = MAX(ibase_array(i,j)-1, 1)
        k_highest = MIN(itop_array(i,j)+1,(nz-1))

        DO WHILE (k <= k_highest)

          IF(clouds_3d(i,j,k+1) >= thresh_cvr                           &
                .AND. clouds_3d(i,j,k) < thresh_cvr                     &
                .OR. k == 1 .AND. clouds_3d(i,j,k) >= thresh_cvr) THEN
            cld_base_m = 0.5 * (zs_1d(k) + zs_1d(k+1))
            k_base = k + 1

            k = k + 1
            DO WHILE (k <= nz-1)
              IF(clouds_3d(i,j,k) > thresh_cvr                          &
                    .AND. clouds_3d(i,j,k+1) <= thresh_cvr) THEN

                cld_top_m = 0.5 * (zs_1d(k) + zs_1d(k+1))
                k_top = k
!
!-----------------------------------------------------------------------
!
!  We have now defined a cloud base and top
!
!-----------------------------------------------------------------------
!
                kb = k_base
                kt = k_top
!
!-----------------------------------------------------------------------
!
!  Make sure cloud base and top stay in the model domain
!
!-----------------------------------------------------------------------
!
                kb = MIN(kb,nz-1)
                kt = MIN(kt,nz-1)
!
                IF(l_prt)                             &
                  PRINT*,i,j,' Cld bs/top index',kb,kt

                IF(l_flag_cldtyp) THEN ! Get 1d Stab.& Cld type
                  CALL get_stability (nz,t_1d,zs_1d,p_mb_1d             &
                                           ,kb,kt,dte_dz_1d,tem1)

                  IF(l_prt) THEN
                    WRITE(6,*)
                    PRINT*,i,j
                    WRITE(6,641) (1000.0*dte_dz_1d(kr),kr=2,34)
                    641           FORMAT(1X,'dTdz(km):',11F6.2)
                  END IF

                  DO k1 = kb,kt
                    CALL get_cloudtype(t_1d(k1),dte_dz_1d(k1)           &
                             ,cld_base_m,cld_top_m,itype,c2_type)
!
                    IF(radar_3d(i,j,k1) > 45.) THEN
                      itype = 10 ! CB
                    END IF

                    IF(l_prt) THEN
                      WRITE(6,*)
                      WRITE(6,642) i,j,k1,t_1d(k1),dte_dz_1d(k1)        &
                                  ,0.001*cld_base_m,0.001*cld_top_m     &
                                  ,radar_3d(i,j,k1),itype,c2_type
                      642 FORMAT(1X,3I3,' t=',f6.1,' dtz=',f6.2,' base/top=' &
                                   ,2F7.3,' refl=',e7.2,' typ=',i2,a3)
                    END IF

                    cldtyp_1d(k1) = itype
                    cldpcp_type_3d(i,j,k1) = itype
                  END DO  !k1
                END IF  ! l_flag_cldtyp.eq.true
!
!-----------------------------------------------------------------------
!
!  We have now defined a cloud base and top
!
!-----------------------------------------------------------------------
!
                IF(iflag_slwc /= 0) THEN

                  n_slwc_call = n_slwc_call + 1

                  IF(iflag_slwc < 10) THEN ! simple adiabatc scheme
                    CALL get_slwc1d (nz,cld_base_m,cld_top_m,kb,kt      &
                        ,zs_1d,t_1d,p_pa_1d,iflag_slwc,slwc_1d)

                  ELSE ! iflag_slwc > 10, new Smith-Feddes scheme
                    mode = 1
                    DO k1 = 1,nz ! Initialize
                      slwc_1d(k1) = 0.0
                      cice_1d(k1) = 0.0
                      ctmp_1d(k1) = temp_3d(i,j,k1)
                    END DO
!
!-----------------------------------------------------------------------
!
!  QC the data going into SMF
!
!-----------------------------------------------------------------------
!
                    IF(cld_top_m > zs_1d(nz-1) - 110.) THEN
                      cld_top_qc_m = zs_1d(nz-1) - 110.
                      cld_base_qc_m = MIN(cld_base_m,cld_top_qc_m - 110.)
                      WRITE(6,'(1x,a,2I4,3I6)')                         &
                            'QC data for SMF, ht/bs/tp',i,j,            &
                            nint(zs_1d(nz-1)),nint(cld_base_qc_m),      &
                            nint(cld_top_qc_m)

                    ELSE ! normal case
                      cld_top_qc_m = cld_top_m
                      cld_base_qc_m = cld_base_m
                    END IF
!
                    IF(l_prt) THEN
                      WRITE(6,*) '== cloud_lwc: Calling GET_SFM_1D =='
                      WRITE(6,643) i,j,0.001*cld_top_qc_m,              &
                                   0.001*cld_base_qc_m
                      643 FORMAT(1X,2I3,' ctop/base_qc=',2F7.3)
                    END IF

                    CALL get_sfm_1d(nz,cld_base_qc_m,cld_top_qc_m       &
                                     ,zs_1d,p_mb_1d,t_1d                &
                                     ,slwc_1d,cice_1d,ctmp_1d,l_prt)
!
                    IF(l_prt) THEN
                      WRITE(6,*) ' After get_sfm_1d: ',i,j
                      WRITE(6,644) (slwc_1d(kr),kr=1,35)
                      644 FORMAT(1X,' lwc=',7E10.3)
                      WRITE(6,645) (cice_1d(kr),kr=1,35)
                      645 FORMAT(1X,' ice=',7E10.3)
                    END IF

                  END IF ! iflag_slwc < 10
                END IF ! iflag_slwc .ne. 0

                DO k1 = kb,kt ! Loop through the cloud layer
!
!-----------------------------------------------------------------------
!
                  IF(iflag_slwc /= 0) THEN
!                   IF(slwc_1d(k1) > 0.) THEN
!                     IF(istat_radar == 1 .AND. temp_3d(i,j,k1) < 273.15) THEN
!
!-----------------------------------------------------------------------
!
!  Apply Radar Depletion when temperature is below 0 degree (C)
!
!-----------------------------------------------------------------------
!
!                       IF(cldtyp_1d(k1) /= 10) THEN ! Not CB
!                         IF(radar_3d(i,j,k1) <= 8.) THEN ! No Depletion
!                           CONTINUE
!                         ELSE IF(radar_3d(i,j,k1) > 11.) THEN
!                                                  ! Total Depltion
!                           cice_1d(k1) = cice_1d(k1) + slwc_1d(k1)
!                           slwc_1d(k1) = 0.0
!
!                           IF(l_prt) THEN
!                             WRITE(6,*) ' **** Radar Depeletion--total (>11dBZ)'
!                             WRITE(6,646) i,j,k1,radar_3d(i,j,k1),cice_1d(k1) &
!                                         ,slwc_1d(k1)
!                           END IF
!                         ELSE                        ! Ramped Depletion
!                           ramp = 1.0 - ((radar_3d(i,j,k1) - 8.) / 3.)
!                           cice_1d(k1) = cice_1d(k1) + slwc_1d(k1)     &
!                                                            * (1.-ramp)
!                           slwc_1d(k1) = slwc_1d(k1) * ramp
!
!                            IF(l_prt) THEN
!                              WRITE(6,*) ' **** Radar Depeletion--partial ( 9-11dBZ)'
!                              WRITE(6,646) i,j,k1,radar_3d(i,j,k1),cice_1d(k1) &
!                                          ,slwc_1d(k1)
!                            END IF
!
!                          END IF  ! Dbz Sufficient for Depletion
!                       END IF  ! Not a CB
!                     END IF ! istat_radar=1 & temp < 273.15
!                   END IF ! slwc_1d(k1) .gt. 0.
!
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Pull lwc/ice from 1d array to 3d array
!
!-----------------------------------------------------------------------
!
                    IF(slwc_1d(k1) > 0.) slwc_3d(i,j,k1)=slwc_1d(k1)
                    IF(cice_1d(k1) > 0.) cice_3d(i,j,k1)=cice_1d(k1)
                    ctmp_3d(i,j,k1)=ctmp_1d(k1)
!
                    IF(l_prt) THEN
                      WRITE(6,*) ' cloud_lwc: final lwc and ice content:'
                      WRITE(6,646) i,j,k1,radar_3d(i,j,k1),cice_3d(i,j,k1) &
                                  ,slwc_3d(i,j,k1)
                      646 FORMAT(1x,3i3,' refl=',e7.2,' ice=',e10.3,' lwc=',e10.3)
                    END IF
                  END IF ! iflag_slwc .ne. 0

                END DO ! k1

                goto 1000

              END IF ! Found Cloud Top, cvr(k) > 0.1 & cvr(k+1) < 0.1

              k = k + 1

            END DO ! k,  while (k .le. nz-1)

          END IF ! Found Cld_Bs, cld(k+1)>.1&cld(k)<.1 OR k=1&cld(k)>.1

          1000 k = k + 1

        END DO ! k,  while (k .le. k_highest)

      END IF ! ibase_array.NE.nz, (At least one layer exists)

    END DO ! j
  END DO ! i
!
!-----------------------------------------------------------------------
!
!  Quality check cloud type field
!
!-----------------------------------------------------------------------
!
  IF (mp_opt > 0) THEN  ! cldpcp_type_3d is used inside cloud_type_qc
                        ! Actually it is an INTEGER. We assume both
                        ! INTEGER and REAL are the same size here.
!    CALL acct_interrupt(mp_acct)
    CALL mpsendrecv2dew(cldpcp_type_3d,nx,ny,nz,ebc,wbc,0,tem1)
    CALL mpsendrecv2dns(cldpcp_type_3d,nx,ny,nz,nbc,sbc,0,tem1)
  END IF

  IF(l_cb_qc) THEN
    CALL cloud_type_qc (nx,ny,nz,zs_3d,clouds_3d,cldpcp_type_3d)
  END IF
!
!-----------------------------------------------------------------------
!
!  Calculate in-cloud w-field.
!
!-----------------------------------------------------------------------
!
  DO i = 1,nx-1
    DO j = 1,ny-1
      l_cloud = .false.
!
      IF( ibase_array(i,j) /= nz) THEN ! At least one lyr exists
        l_cloud = .true.

        DO k = 1,nz-1                      ! Initialize
          t_1d(k) = temp_3d(i,j,k)
          zs_1d(k) = zs_3d(i,j,k)
          p_pa_1d(k) = p_3d(i,j,k)
          p_mb_1d(k) = 0.01*p_3d(i,j,k)
          cldtyp_1d(k) = cldpcp_type_3d(i,j,k)
        END DO
      END IF ! ibase_array.NE.nz, (At least one layer exists)

      IF(l_flag_incld_w) THEN
        IF(l_cloud) THEN
          IF(l_prt) THEN
            WRITE(6,*) ' cloud_lwc: Calling incld_w',i,j
            WRITE(6,650) (cldtyp_1d(kr),kr=2,34)
            650         FORMAT('CTP',11I7)
            WRITE(6,652) (zs_1d(kr),kr=2,34)
            652         FORMAT('zs:',11F7.0)
          END IF

          DO k = 2,nz-1
            zp_1d(k) = 0.5*(zs_1d(k-1)+zs_1d(k))
          END DO
          zp_1d(1) = zs_1d(1) - (zp_1d(2)-zs_1d(1))

          CALL cloud_w(nz,cldtyp_1d,zp_1d,w_1d                          &
                       ,r_missing,wmhr_cu,wmhr_sc,wc_st)

          IF(l_prt) THEN
            WRITE(6,653) (w_1d(kr),kr=2,34)
            653 FORMAT('w:',11F7.2)
            WRITE(6,654) (p_mb_1d(kr),kr=2,34)
            654 FORMAT('p:',11F7.1)
          END IF

          DO k = 1,nz-1 ! Transfer the 1D w into the output w_3d array
            IF(w_1d(k) /= r_missing) THEN
              w_3d(i,j,k) = w_1d(k)
              w_tem(i,j,k) = w_1d(k)
            END IF
          END DO

        END IF ! l_cloud
      END IF ! l_flag_incld_w

    END DO ! j
  END DO ! i

  IF(clddiag == 1) THEN
    IF (cld_files == 1) THEN
      varid='cloudw'
      varname='Cloud Vert Veloc'
      varunits='m/s'
      CALL wrtvar1(nx,ny,nz,w_tem,varid,varname,varunits,               &
                   time,runname,dirname,istatwrt)
    END IF
  END IF
!
!-----------------------------------------------------------------------
!
  CALL mptotali(n_cld_columns)
  CALL mptotali(n_slwc_call)
  IF (myproc == 0) WRITE(6,'(/1X,2(a,i8)/)')    &
    ' n_cld_columns = ', n_cld_columns,', n_slwc_call = ',n_slwc_call
!
!-----------------------------------------------------------------------
!
!  Simulate Radar Data
!
!-----------------------------------------------------------------------
!
  IF(.true.) THEN
    DO i = 1,nx
      DO j = 1,ny
        l_mask_pcptype(i,j) = .false.
      END DO ! j
    END DO ! i
  END IF
!
!-----------------------------------------------------------------------
!
!  Get 3D precip type
!
!-----------------------------------------------------------------------
!
  IF(istat_radar == 1) THEN ! Get 3D precip type

    IF (myproc == 0) WRITE(6,*) ' Analyzing cloud/precipitation type.'

    CALL pcp_type_3d(nx,ny,nz,temp_3d,rh_3d,p_3d                        &
                    ,radar_3d,l_mask_pcptype                            &
                    ,cldpcp_type_3d,istat_1)
    IF(istat_1 /= 1) THEN
      WRITE(6,*) 'Bad status returned by PCP_TYPE_3D, aborting'
      RETURN
    END IF

  END IF ! Valid Radar Data

  IF(clddiag == 1) THEN
    IF (cld_files == 1) THEN
      varid='cptype'
      varname='Precipitation Type'
      varunits='Category Number'
      CALL wrtvar1(nx,ny,nz,cldpcp_type_3d,varid,varname,varunits,       &
                  time,runname,dirname,istatwrt)
    END IF
  END IF
!
!-----------------------------------------------------------------------
!
!  Computing Icing Severity Index field
!
!-----------------------------------------------------------------------
!
  IF(l_flag_icing) THEN

    IF (myproc == 0) WRITE(6,*)' Computing Icing Severity Index field'

    DO i = 1,nx-1
      DO j = 1,ny-1
        DO k = 1,nz-1
          iarg = cldpcp_type_3d(i,j,k)
          i_precip_type = iarg/16
          i_cldtyp  = iarg - i_precip_type*16

          IF(iarg > 0) THEN ! clouds or precip are present
            rho = p_3d(i,j,k)/(rair*temp_3d(i,j,k))
            lwc_g_m3 = slwc_3d(i,j,k)*rho      ! g/m3 for icing
            CALL isi3(lwc_g_m3,temp_3d(i,j,k)-273.15                    &
                                ,i_cldtyp,i_precip_type,xindex)
            iarg = xindex
            icing_index_3d(i,j,k) = iarg

          ELSE
            icing_index_3d(i,j,k) = 0

          END IF

        END DO ! k
      END DO ! j
    END DO ! i

    IF(clddiag == 1) THEN
      IF (cld_files == 1) THEN
        varid='icingd'
        varname='Icing Index'
        varunits='Index Number'
        CALL wrtvar1(nx,ny,nz, icing_index_3d, varid,varname,varunits,  &
                    time,runname,dirname,istatwrt)
      END IF
    END IF
!
  END IF ! l_flag_icing =.true.
!
!-----------------------------------------------------------------------
!
  istatus = 1
  RETURN
END SUBROUTINE cloud_lwc
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE GET_STABILITY               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE get_stability (nz,t_1d,zs_1d,p_mb_1d,kbtm,ktop               &
           ,dte_dz_1d,thetae_1d)
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  This routine returns stability at a given level given
!  1D temperature and pressure array inputs
!
!-----------------------------------------------------------------------
!
!  AUTHOR:  Jian Zhang
!  05/96    Based on LAPS cloud analysis code of 07/95
!
!  MODIFICATION HISTORY:
!
!  05/11/96  (J. Zhang)
!            Modified for ADAS format. Added full documentation.
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
!  INPUT:
  INTEGER :: nz         ! number of vertical model levels
  REAL :: t_1d(nz)      ! temperature (degree Kelvin) profile
  REAL :: zs_1d(nz)     ! heights (m MSL) of each level
  REAL :: p_mb_1d(nz)   ! pressure (mb) at each level
  INTEGER :: kbtm,ktop  ! indices of the bottom and top cloud layer
!
!  OUTPUT:
  REAL :: dte_dz_1d(nz) ! stability array
!
!  LOCAL:
  REAL :: thetae_1d(nz) ! (equivalent) potential temperature.
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: k,km1,kp1,klow,khigh
  REAL :: os_fast
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
!  Calculate Stability
!
!-----------------------------------------------------------------------
!
  klow  = MAX(kbtm-1,1)
  khigh = MIN(ktop+1,nz-1)

  DO k = klow,khigh
    thetae_1d(k)  = os_fast(t_1d(k), p_mb_1d(k))
  END DO ! k

  dte_dz_1d=0.

  DO k = kbtm,ktop
    km1  = MAX(k-1,1)
    kp1  = MIN(k+1,nz-1)

    IF( (zs_1d(kp1) - zs_1d(km1)) <= 0.) THEN
      print *, ' Error in get_stability '
      print *, ' k,kp1,km1 = ',k,kp1,km1
      print *, ' zs_1d(kp1),zs_1d(km1)= ',zs_1d(kp1),zs_1d(km1),        &
                 (zs_1d(kp1) - zs_1d(km1))
      STOP
    ELSE
      dte_dz_1d(k) = (thetae_1d(kp1) - thetae_1d(km1))                  &
                           / (zs_1d(kp1) - zs_1d(km1))
    END IF
  END DO ! k

  RETURN
END SUBROUTINE get_stability


!
!##################################################################
!##################################################################
!######                                                      ######
!######            FUNCTION OS_FAST                          ######
!######                                                      ######
!##################################################################
!##################################################################
!

  FUNCTION os_fast(tk,p)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  THIS FUNCTION RETURNS THE EQUIVALENT POTENTIAL TEMPERATURE OS
!  (K) FOR A PARCEL OF AIR SATURATED AT TEMPERATURE T (K)
!  AND PRESSURE P (MILLIBARS).
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (BAKER,SCHLATTER)
!  05/17/1982
!
!
!  MODIFICATION HISTORY:
!  05/11/96 (Jian Zhang)
!  Modified for ADAS grid. Add document stuff.
!
!-----------------------------------------------------------------------
!
!  Variables declaration
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  INPUT:
  REAL :: tk     ! temperature in kelvin
  REAL :: p      ! pressure in mb
!
!  OUTPUT:
  REAL :: os_fast  ! equivalent potential temperature
!
!  LOCAL:
  REAL :: b            ! empirical const. approx.= latent heat of
                     ! vaporiz'n for water devided by the specific
                     ! heat at const. pressure for dry air.
  DATA b/2.6518986/

  REAL :: tc,x,w
  REAL :: eslo
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  tc = tk - 273.15
!
!-----------------------------------------------------------------------
!
!  From W routine
!
!-----------------------------------------------------------------------
!
  x= eslo(tc)
  w= 622.*x/(p-x)

  os_fast= tk*((1000./p)**.286)*(EXP(b*w/tk))

  RETURN
  END FUNCTION os_fast



!
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE GET_CLOUDTYPE               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE get_cloudtype(temp_k,dte_dz,cbase_m,ctop_m                   &
           ,itype,c2_type)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  This routine returns cloud type at a given point given
!  temperature and stability inputs
!
!-----------------------------------------------------------------------
!
!  AUTHOR:  Jian Zhang
!  05/96    Based on the LAPS cloud analysis code of 05/1995
!
!  MODIFICATION HISTORY:
!
!  05/11/96  (J. Zhang)
!            Modified for ADAS format. Added full documentation.
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
!  INPUT:
  REAL :: temp_k       ! temperature
  REAL :: dte_dz       ! stability factor
  REAL :: cbase_m      ! height at cloud base level
  REAL :: ctop_m       ! height at cloud top level
!
!  OUTPUT:
  INTEGER :: itype     ! cloud type index
  CHARACTER (LEN=2) :: c2_type
!
!  LOCAL:
  CHARACTER (LEN=2) :: c2_cldtyps(10)

  DATA c2_cldtyps /'St','Sc','Cu','Ns','Ac'                             &
                  ,'As','Cs','Ci','Cc','Cb'/
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  REAL :: depth_m,temp_c
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  temp_c = temp_k - 273.15
  depth_m = ctop_m - cbase_m
!
!-----------------------------------------------------------------------
!
!  Go from Stability to Cloud Type
!
!-----------------------------------------------------------------------
!
  IF ( temp_c >= -10.) THEN
    IF (dte_dz >= +.001) THEN
      itype = 1      ! St
    ELSE IF (dte_dz < +.001 .AND. dte_dz >= -.001)  THEN
      itype = 2      ! Sc
    ELSE IF (dte_dz < -.001 .AND. dte_dz >= -.005)  THEN
      itype = 3      ! Cu
    ELSE ! dte_dz .lt. -.005
      IF(depth_m > 5000) THEN
        itype = 10   ! Cb
      ELSE  ! depth < 5km
        itype = 3    ! Cu
      END IF
    END IF

  ELSE IF (temp_c < -10. .AND. temp_c >= -20.) THEN

    IF (dte_dz < 0.) THEN
      IF(depth_m > 5000) THEN
        itype = 10   ! Cb
      ELSE
        itype = 5    ! Ac
      END IF
    ELSE
      itype = 6      ! As
    END IF

  ELSE               ! temp_c.lt.-20.

    IF (dte_dz >= +.0005) THEN
      itype = 7      ! Cs
    ELSE IF (dte_dz < +.0005 .AND. dte_dz >= -.0005) THEN
      itype = 8      ! Ci
    ELSE             ! dte_dz .lt. -.0005
      itype = 9      ! Cc
    END IF

    IF(depth_m > 5000 .AND. dte_dz < -.0000) THEN
      itype = 10     ! Cb
    END IF

  END IF

  c2_type = c2_cldtyps(itype)

  RETURN
END SUBROUTINE get_cloudtype






!
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE GET_SLWC1D                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE get_slwc1d (nk,cbase_m,ctop_m,kbase,ktop                     &
           ,zs_1d,t_1d,p_pa_1d,iflag_slwc,slwc_1d)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  This routine calls a subroutine "lwc_rep" which calculates
!  the adiabatic liquid water content.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:  Jian Zhang
!  05/96    Based on the LAPS cloud analysis code of 07/1995
!
!  MODIFICATION HISTORY:
!
!  05/13/96  (Jian Zhang)
!            Modified for ADAS format. Added full documentation.
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
!  INPUT:
  INTEGER :: iflag_slwc     ! indicator for LWC scheme option
  INTEGER :: nk             ! number of model vertical levels
  REAL :: t_1d(nk)          ! temperature (k) in one model column
  REAL :: zs_1d(nk)         ! heights (m) at grd pts in one model column
  REAL :: p_pa_1d(nk)       ! pressure (pa) in one model column
  REAL :: cbase_m,ctop_m    ! heights (m) of cloud base and top levels
  INTEGER :: kbase,ktop        ! vertical index of cloud base and top levels
!
!  OUTPUT:
  REAL :: slwc_1d(nk)       ! estimated adiabatic liquid water
!
!  LOCAL:
  INTEGER :: i_status1,i_status2 ! flag for subroutine calling
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: k
  REAL :: p_low,p_high,cbase_pa,cbase_k,ctop_k,frac_k                   &
       ,grid_top_pa,grid_top_k
  REAL :: fraction,thickness,dlog_space
  REAL :: adiabatic_lwc,adjusted_lwc,adjusted_slwc
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
!  Initialize
!
!-----------------------------------------------------------------------
!
  DO k = 1,nk
    slwc_1d(k) = 0.0
  END DO

  IF(ctop_m > cbase_m) THEN
!
!-----------------------------------------------------------------------
!
!  Determine Lowest and Highest Grid Points within the cloud
!
!-----------------------------------------------------------------------
!
    IF(ktop >= kbase .AND. kbase >= 2) THEN
!
!-----------------------------------------------------------------------
!
!  Get cloud base pressure and temperature
!
!-----------------------------------------------------------------------
!
      cbase_pa = -999.         ! Default value is off the grid
      DO k = 1,nk-2
        IF(zs_1d(k+1) > cbase_m .AND. zs_1d(k) <= cbase_m) THEN
          thickness = zs_1d(k+1) - zs_1d(k)
          fraction = (cbase_m - zs_1d(k))/thickness
          p_low = p_pa_1d(k)
          p_high = p_pa_1d(k+1)
          dlog_space = LOG(p_high/p_low)
          cbase_pa = p_low * EXP(dlog_space*fraction)
        END IF
      END DO ! k

      frac_k=(cbase_m-zs_1d(kbase-1))/(zs_1d(kbase)-zs_1d(kbase-1))
      IF(frac_k /= fraction)                                            &
          PRINT*,' **GET_SLWC1D**  frac=',fraction,' frac_k=',frac_k

      cbase_k = t_1d(kbase-1)*(1.0-frac_k) + t_1d(kbase)*frac_k
!
!-----------------------------------------------------------------------
!
!  Get cloud top temperature
!
!-----------------------------------------------------------------------
!
      frac_k = (ctop_m-zs_1d(ktop-1)) / (zs_1d(ktop)-zs_1d(ktop-1))
      ctop_k = t_1d(ktop-1)*(1.0 - frac_k) + t_1d(ktop) * frac_k
!
!-----------------------------------------------------------------------
!
!  Calculate SLWC at each vertical grid point. For each level
!  we use an assumed cloud extending from the actual cloud base
!  to the height of the grid point in question.
!
!-----------------------------------------------------------------------
!
      DO k=kbase,ktop
        grid_top_pa = p_pa_1d(k)
        grid_top_k = t_1d(k)

        CALL slwc_revb(cbase_pa,cbase_k                                 &
                  ,grid_top_pa,grid_top_k,ctop_k                        &
                  ,adiabatic_lwc,adjusted_lwc,adjusted_slwc             &
                  ,i_status1,i_status2)
!
        IF(i_status2 == 1) THEN
          IF(iflag_slwc == 1) THEN
            slwc_1d(k) = adiabatic_lwc
          ELSE IF(iflag_slwc == 2) THEN
            slwc_1d(k) = adjusted_lwc
          ELSE IF(iflag_slwc == 3) THEN
            slwc_1d(k) = adjusted_slwc
          END IF
        ELSE
          WRITE(6,*)' Error Detected in SLWC'
        END IF
      END DO ! k
    END IF ! ktop > kbase & kbase > 2,  thick enough cloud exists
  END IF ! ctop_m > cbase_m,  cloud exists

  RETURN
END SUBROUTINE get_slwc1d






SUBROUTINE slwc_revb(cb_pa,cb_k,gt_pa,gt_k,ct_k,                        &
           adiabatic_lwc,adjusted_lwc,adjusted_slwc,                    &
           i_status1,i_status2)
!
!.......................HISTORY.............................
!
!     WRITTEN: CA. 1982 BY W. A. COOPER IN HP FORTRAN 4
!
!....... CALCULATES TEMPERATURE T AND LIQUID WATER CONTENT FROM
!..      CLOUD BASE PRESSURE P0 AND TEMPERATURE T0, FOR ADIABATIC
!..      ASCENT TO THE PRESSURE P.
!..     ->  INPUT:  CLOUD BASE PRESSURE P0 AND TEMPERATURE T0
!..                 PRESSURE AT OBSERVATION LEVEL P
!..     ->  OUTPUT: "ADIABATIC" TEMPERATURE T AND LIQUID WATER CONTENT
!
!     MODIFIED: November 1989 by Paul Lawson for LAPS/WISP.  Routine
!               now calculates adiabatic liquid water content
!               (ADIABATIC_LWC) using cloud base pressure and grid-top
!               temperature and pressure.  Also calculated are ADJUSTED_LWC,
!               which adjusts ADIABATIC_LWC using an empirical cloud
!               water depletion algorithm, and ADJUSTED_SLWC, which is
!               ADIABATIC_LWC in regions where T < 0 C adjusted
!               using an empirical algorithm by Marcia Politovich.
!
!               Subroutine is now hardwired for stratiform cloud only.
!               Can be modified to include Cu with input from LAPS main.
!
!               revb: ca 12/89 Calculate adiabatic lwc by going from cloud
!                     base to LAPS grid level instead to cloud top, thus
!                     helping to better calculate in layer clouds.
!                     Add TG (grid temperature) to calcualtion.
!
!               revc: 2/27/90 Correct error in code.  Zero-out slwc when grid
!                     temperature (GT) > 0.
!
!               J.Z.: 4/7/97 Correct error in code
!                     Grid temperature should be TG, not GT.
!
!
!     OUTPUTS:  ADIABATIC_LWC
!               ADJUSTED_LWC
!               ADJUSTED_SLWC
!               I_STATUS1 - 1 when -20 < cld_top_temp < 0 for Stratus
!                           0 Otherwise
!               I_STATUS2 - 1 when valid input data provided from main
!
  DATA eps/0.622/,cpd/1.0042E3/,cw/4.218E3/,rd/287.05/,alhv/2.501E6/
  INTEGER :: cty
  INTEGER :: i_status1, i_status2
  i_status1=1
  i_status2=1
!   2 Print *,'ENTER: P-BASE(mb), T-BASE(C), P-TOP, T-TOP, CLD TYPE'
!  READ(5,*) P0, T0, P, CTT, CTY
!  If(CTY.ne.0.and.CTY.ne.1) Go to 2
!
!  Hardwire cloud type (CTY) for stratus for now
!
  cty=0
!
!.....Convert Pa to mb and Kelvin to Celcius
!
  p0 = cb_pa/100.
  p  = gt_pa/100.
  t0 = cb_k - 273.15
  tg = gt_k - 273.15
  ctt= ct_k - 273.15
!  Print *, 'CTT in Sub = ', CTT
!
!  Check for valid input data...
!
  IF(p0 > 1013..OR.p0 < 50.) THEN
    i_status2=0
    RETURN
  ELSE
  END IF
!
!
  IF(t0 > 50..OR.t0 < -70.) THEN
    i_status2=0
    RETURN
  ELSE
  END IF
!
!
  IF(p > 1013..OR.p < 50.) THEN
    i_status2=0
    RETURN
  ELSE
  END IF
!
!  Set I_STATUS1 = F if 0 < cld top < -20 C (for stratus).
!
  IF(tg >= 0..OR.ctt < -20.) i_status1=0
!
  tk=t0+273.15
  e=vapor(t0)
  r=eps*e/(p0-e)
  cpt=cpd+r*cw
  thetaq=tk*(1000./(p0-e))**(rd/cpt)*EXP(alhv*r/(cpt*tk))
! 1ST APPROX
  t1=tk
  e=vapor(t1-273.15)
  rv=eps*e/(p-e)
  t1=thetaq/((1000./(p-e))**(rd/cpt)*EXP(alhv*rv/(cpt*t1)))
! SUCCESSIVE APPROXIMATIONS
  DO i=1,10
    e=vapor(t1-273.15)
    rv=eps*e/(p-e)
    t1=(thetaq/((1000./(p-e))**(rd/cpt)*EXP(alhv*rv/(cpt*t1)))          &
        +t1)/2.
    t=t1-273.15
!  Print *, P0,T0,P,T,E,RV,THETAQ
  END DO
! GET LWC
  e=vapor(t)
  rv=eps*e/(p-e)
  tw=r-rv
  adiabatic_lwc=tw*p*28.9644/(8.314E7*t1)*1.e9
  IF(adiabatic_lwc < 0.) adiabatic_lwc=0.
!  Print *, 'Adiabtic LWC = ', ADIABATIC_LWC
  IF(tg >= 0.) THEN
!
    adjusted_slwc=0.                                          ! Added 2/27/90
!

    IF(cty == 0.) THEN
      IF(ctt < -20.) THEN
        adjusted_lwc=0.
      ELSE IF(ctt < -15..AND.ctt >= -20.) THEN
        adjusted_lwc=adiabatic_lwc/8.
      ELSE IF(ctt < -10..AND.ctt >= -15.) THEN
        adjusted_lwc=adiabatic_lwc/4.
      ELSE
        adjusted_lwc=adiabatic_lwc/2.
      END IF
    ELSE
      IF(ctt < -25.) THEN
        adjusted_lwc=0.
      ELSE IF(ctt < -15..AND.ctt >= -25.) THEN
        adjusted_lwc=adiabatic_lwc/8.
      ELSE IF(ctt < -10..AND.ctt >= -15.) THEN
        adjusted_lwc=adiabatic_lwc/4.
      ELSE
        adjusted_lwc=adiabatic_lwc/2.
      END IF
    END IF
  ELSE
    IF(cty == 0.) THEN
      IF(ctt < -20.) THEN
        adjusted_lwc=0.
        adjusted_slwc=0.
      ELSE IF(ctt < -15..AND.ctt >= -20.) THEN
        adjusted_lwc=adiabatic_lwc/8.
        adjusted_slwc=adiabatic_lwc/8.
      ELSE IF(ctt < -10..AND.ctt >= -15.) THEN
        adjusted_lwc=adiabatic_lwc/4.
        adjusted_slwc=adiabatic_lwc/4.
      ELSE
        adjusted_lwc=adiabatic_lwc/2.
        adjusted_slwc=adiabatic_lwc/2.
      END IF
    ELSE
      IF(ctt < -25.) THEN
        adjusted_lwc=0.
        adjusted_slwc=0.
      ELSE IF(ctt < -15..AND.ctt >= -25.) THEN
        adjusted_lwc=adiabatic_lwc/8.
        adjusted_slwc=adiabatic_lwc/8.
      ELSE IF(ctt < -10..AND.ctt >= -15.) THEN
        adjusted_lwc=adiabatic_lwc/4.
        adjusted_slwc=adiabatic_lwc/4.
      ELSE
        adjusted_lwc=adiabatic_lwc/2.
        adjusted_slwc=adiabatic_lwc/2.
      END IF
    END IF
  END IF
!  Print *,'Adjusted LWC = ', ADJUSTED_LWC
!  Print *,'Adjusted SLWC = ', ADJUSTED_SLWC
END SUBROUTINE slwc_revb





!  FUNCTION TO CALCULATE VAPOR PRESSURE:
!

  FUNCTION vapor(tfp)
! INPUT IS IN DEGREES C.  IF GT 0, ASSUMED TO BE DEW POINT.  IF
! LESS THAN 0, ASSUMED TO BE FROST POINT.
! ROUTINE CODES GOFF-GRATCH FORMULA
  tvap=273.16+tfp
  IF(tfp > 0.) GO TO 1
! THIS IS ICE SATURATION VAPOR PRESSURE
  IF(tvap <= 0) tvap=1E-20
  e=-9.09718*(273.16/tvap-1.)-3.56654*ALOG10(273.16/tvap)               &
      +0.876793*(1.-tvap/273.16)
  vapor=6.1071*10.**e
  RETURN
  1    CONTINUE
! THIS IS WATER SATURATION VAPOR PRESSURE
  IF(tvap <= 0) tvap=1E-20
  e=-7.90298*(373.16/tvap-1.)+5.02808*ALOG10(373.16/tvap)               &
      -1.3816E-7*(10.**(11.344*(1.-tvap/373.16))-1.)                    &
      +8.1328E-3*(10.**(3.49149*(1-373.16/tvap))-1)
  vapor=1013.246*10.**e
  RETURN
  END FUNCTION vapor
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE GET_SFM_1D                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE get_sfm_1d (nz,zcb,zctop,zs_1d,p_mb_1d,t_1d,ql,qi,cldt,      &
                       l_prt)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!c-----------------------------------------------------------------
!c
!c    This is the streamlined version of the Smith-Feddes
!c    and Temperature Adjusted LWC calculation methodologies
!c    produced at Purdue University under sponsorship
!c    by the FAA Technical Center.
!c
!c    Currently, this subroutine will only use the Smith-
!c    Feddes and will only do so as if there are solely
!c    stratiform clouds present, however, it is very easy
!c    to switch so that only the Temperature Adjusted
!c    method is used.
!c
!c    Dilution by glaciation is also included, it is a
!c    linear function of in cloud temperature going from
!c    all liquid water at -10 C to all ice at -30 C
!c    as such the amount of ice is also calculated
!
!-----------------------------------------------------------------------
!
!  AUTHOR:  Jian Zhang
!  05/96    Based on the LAPS cloud analysis code of 07/1995
!
!  MODIFICATION HISTORY:
!
!  05/16/96 (Jian Zhang)
!           Modified for ADAS format. Added full documentation.
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
!
!-----------------------------------------------------------------------
!
!  INPUT:
  INTEGER :: nz             ! number of model vertical levels
  REAL :: zs_1d(nz)         ! physical height (m) at each scalar level
  REAL :: p_mb_1d(nz)       ! pressure (mb) at each level
  REAL :: t_1d(nz)          ! temperature (K) at each level

  REAL :: zcb               ! cloud base height (m)
  REAL :: zctop             ! cloud top height (m)
!
!  OUTPUT:
  REAL :: ql(nz)            ! liquid water content (g/kg)
  REAL :: qi(nz)            ! ice water content (g/kg)
  REAL :: cldt(nz)
!
!  LOCAL:
  REAL :: calw(300)
  REAL :: cali(300)
  REAL :: catk(300)
  REAL :: entr(300)
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  REAL :: dz,rv,rair,grav,cp,rlvo,rlso,dlvdt,eso
  REAL :: c,a1,b1,c1,a2,b2,c2
  REAL :: delz,delt,cldbtm,cldbp,cldtpt,tbar
  REAL :: arg,fraclw,tlwc
  REAL :: temp,press,zbase,alw,zht,ht,y
  REAL :: rl,es,qvs1,p,des,dtz,es2,qvs2
  INTEGER :: i,j,k,nlevel,nlm1,ip,kctop,kctop1,kcb,kcb1
  REAL :: dtdz,dttdz,zcloud,entc,tmpk
  LOGICAL :: l_prt
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
!  Initialize 1d liquid water and ice arrays (for 100m layers)
!
!-----------------------------------------------------------------------
!
  DO i=1,200
    calw(i)=0.0
    cali(i)=0.0
  END DO
!      if(i_prt.le.20) then
!        i_prt=i_prt+1
!        l_prt=.true.
!      else
!        l_prt=.false.
!      endif
!
!-----------------------------------------------------------------------
!
!  Preset some constants and coefficients.
!
!-----------------------------------------------------------------------
!
  dz=100.0                ! m
  rv=461.5                ! J/deg/kg
  rair=287.04             ! J/deg/kg
  grav=9.81               ! m/s2
  cp=1004.                ! J/deg/kg
  rlvo=2.5003E+6          ! J/kg
  rlso=2.8339E+6          ! J/kg
  dlvdt=-2.3693E+3        ! J/kg/K
  eso=610.78              ! pa
  c=0.01
  a1=8.4897
  b1=-13.2191
  c1=4.7295
  a2=10.357
  b2=-28.2416
  c2=8.8846
!
!-----------------------------------------------------------------------
!
!  Calculate indices of cloud top and base
!
!-----------------------------------------------------------------------
!
  DO k=1,nz-1
    IF(zs_1d(k) < zcb .AND. zs_1d(k+1) > zcb) THEN
      kcb=k
      kcb1=kcb+1
    END IF
    IF(zs_1d(k) < zctop .AND. zs_1d(k+1) > zctop) THEN
      kctop=k
      kctop1=kctop+1
    END IF
  END DO
  IF(l_prt .AND. myproc == 0) THEN
    WRITE(6,*) ' get_sfm_1d: input at cloud base:'
    WRITE(6,600) zcb,kcb,zs_1d(kcb),t_1d(kcb)                           &
                 ,kcb1,zs_1d(kcb1),t_1d(kcb1)
    600     FORMAT(1X,' base=',f8.0,' kcb=',i2,' ht=',f8.0,' T=',f6.1,  &
                     ' kcb1=',i2,' ht=',f8.0,' T=',f6.1)
    WRITE(6,*) ' get_sfm_1d: input at cloud top:'
    WRITE(6,601) zctop,kctop,zs_1d(kctop),t_1d(kctop)                   &
                 ,kctop1,zs_1d(kctop1),t_1d(kctop1)
    601     FORMAT(1X,' top=',f8.0,' kctop=',i2,' ht=',f8.0,' T=',f6.1, &
                     ' kctop1=',i2,' ht=',f8.0,' T=',f6.1)
  END IF
!
!-----------------------------------------------------------------------
!
!  Obtain cloud base and top conditions
!
!-----------------------------------------------------------------------
!
  delz   = zs_1d(kcb+1)-zs_1d(kcb)
  delt   = t_1d(kcb+1)-t_1d(kcb)
  cldbtm = delt*(zcb-zs_1d(kcb))/delz+t_1d(kcb)
  tbar   = (cldbtm+t_1d(kcb))/2.
  arg    = -grav*(zcb-zs_1d(kcb))/rair/tbar
  cldbp  = p_mb_1d(kcb)*EXP(arg)
  delz   = zs_1d(kctop+1)-zs_1d(kctop)
  delt   = t_1d(kctop+1)-t_1d(kctop)
  cldtpt = delt*(zctop-zs_1d(kctop))/delz+t_1d(kctop)
!
!-----------------------------------------------------------------------
!
!  Calculate cloud lwc profile for cloud base/top pair
!
!-----------------------------------------------------------------------
!
  temp   = cldbtm
  press  = cldbp*100.0
  zbase  = zcb
  nlevel = ((zctop-zcb)/100.0)+1
  IF(nlevel <= 0) nlevel=1
  alw    = 0.0
  calw(1)= 0.0
  cali(1)= 0.0
  catk(1)= temp
  entr(1)= 1.0
  nlm1   = nlevel-1
  IF(nlm1 < 1) nlm1=1
! Prevent array out of bounds
  IF(nlm1 >= 200) nlm1=199
  zht    = zbase

  DO j=1,nlm1
    rl   = rlvo+(273.15-temp)*dlvdt
    arg  = rl*(temp-273.15)/273.15/temp/rv
    es   = eso*EXP(arg)
    qvs1 = 0.622*es/(press-es)
!        rho1 = press/(rair*temp)
    arg  = -grav*dz/rair/temp
    p    = press*EXP(arg)

    IF(l_prt .AND. myproc == 0) THEN
      WRITE(6,605) j,zht,temp,press,1000.0*qvs1,es,rl
      605       FORMAT(1X,i2,' ht=',f8.0,' T=',f6.1,' P=',f9.1,' qvs=', &
                       f7.3,' es=',f6.1,' Lv=',e8.3)
    END IF
!
!-----------------------------------------------------------------------
!
!  Calculate saturated adiabatic lapse rate
!
!-----------------------------------------------------------------------
!
    des   = es*rl/temp/temp/rv
    dtz   = -grav*((1.0+0.621*es*rl/(press*rair*temp))/                 &
                 (cp+0.621*rl*des/press))
    zht   = zht+dz
    press = p
    temp  = temp+dtz*dz
    rl    = rlvo+(273.15-temp)*dlvdt
    arg   = rl*(temp-273.15)/273.15/temp/rv
    es2   = eso*EXP(arg)
    qvs2  = 0.622*es2/(press-es2)

    IF(l_prt .AND. myproc == 0) THEN
      WRITE(6,605) j+1,zht,temp,press,1000.0*qvs2,es2,rl
!605       format(1x,i2,' ht=',f8.0,' T=',f6.1,' P=',f9.1,' qvs=',
!  :           f7.3,' es=',f6.1,' Lv=',e8.3)
    END IF

!        rho2  = press/(rair*temp)
!        alw   = alw+(qvs1-qvs2)*(rho1+rho2)/2.   ! kg/m3

    alw   = alw+(qvs1-qvs2)                   ! kg/kg
    calw(j+1) = alw

    IF (l_prt .AND. myproc == 0) THEN
      WRITE(6,9015) j,1000.0*calw(j+1),zht
      9015      FORMAT(1X,'j=',i3,'  adiab.lwc =',f7.3,'  alt =',f8.0)
    END IF
!
!-----------------------------------------------------------------------
!
!  Reduction of lwc by entrainment
!
!-----------------------------------------------------------------------
!
    ht = (zht-zbase)*.001
!
!c   ------------------------------------------------------------------
!c
!c                          skatskii's curve(convective)
!c
!c   ------------------------------------------------------------------
!c      if(ht.lt.0.3) then
!c        y    = -1.667*(ht-0.6)
!c      elseif(ht.lt.1.0) then
!c        arg1 = b1*b1-4.0*a1*(c1-ht)
!c        y    = (-b1-sqrt(arg1))/(2.0*a1)
!c      elseif(ht.lt.2.9) then
!c        arg2 = b2*b2-4.0*a2*(c2-ht)
!c        y    = (-b2-sqrt(arg2))/(2.0*a2)
!c      else
!c        y    = 0.26
!c      endif
!c
!c   ------------------------------------------------------------------
!c
!c                         warner's curve(stratiform)
!c
!c   ------------------------------------------------------------------
    IF(ht < 0.032) THEN
      y = -11.0*ht+1.0           ! y(ht=0.032) = 0.648
    ELSE IF(ht <= 0.177) THEN
      y = -1.4*ht+0.6915         ! y(ht=0.177) = 0.4437
    ELSE IF(ht <= 0.726) THEN
      y = -0.356*ht+0.505        ! y(ht=0.726) = 0.2445
    ELSE IF(ht <= 1.5) THEN
      y = -0.0608*ht+0.2912      ! y(ht=1.5) = 0.2
    ELSE
      y = 0.20
    END IF
!
!-----------------------------------------------------------------------
!
!  Calculate reduced lwc by entrainment and dilution
!
!  Note at -5 C and warmer, all liquid.   ! changed from -10 KB
!       at -25 C and colder, all ice      ! changed from -30 KB
!       Linear ramp between.
!
!-----------------------------------------------------------------------
!
    IF(temp < 268.15) THEN
      IF(temp > 248.15) THEN
        fraclw=0.05*(temp-248.15)
      ELSE
        fraclw=0.0
      END IF
    ELSE
      fraclw=1.0
    END IF

    tlwc=1000.*y*calw(j+1)                ! g/kg
    calw(j+1)=tlwc*fraclw
    cali(j+1)=tlwc*(1.-fraclw)
    catk(j+1)=temp
    entr(j+1)=y

    IF(l_prt .AND. myproc == 0) THEN
      WRITE(6,*) ' Get_sfm_1d: entrainment dilution'
      WRITE(6,608) j+1,ht,1000.0*tlwc,calw(j+1),cali(j+1)
      608       FORMAT(1X,i2,' ht=',f8.3,' alw=',f8.3,' lwc=',f7.3,     &
                       ' ice=',f7.3)
    END IF

  END DO
!
!-----------------------------------------------------------------------
!
!  Alternative calculation procedure using the observed or
!  inferred in cloud temperature profile
!
!-----------------------------------------------------------------------
!
  IF(.true.) GO TO 455        ! forced goto, diseffect the following
  nlevel = (zctop-zcb)/100.0
  temp   = cldbtm
  press  = cldbp*100.0
  IF(nlevel <= 0) nlevel=0
  alw     = 0.0
  calw(1) = 0.0
  nlm1    = nlevel-1
  IF(nlm1 < 1) nlm1=1
  dtdz = (cldtpt-cldbtm)/(zctop-zcb)
  zht  = zbase
  DO j=1,nlm1
    rl   = rlvo+(temp-273.15)*dlvdt
    arg  = rl*(273.15-temp)/273.15/temp/rv
    es   = eso*EXP(arg)
    qvs1 = 0.622*es/(press-es)
!        rho1 = press/(rair*temp)
    arg  = -grav*dz/rair/temp
    p    = press*EXP(arg)
    des  = es*rl/temp/temp/rv
    dtz  = -grav*((1.0+0.621*es*rl/(press*rair*temp))/                  &
                  (cp+0.621*rl*des/press))
    IF(dtdz < dtz) THEN
      dttdz = dtz-(dtdz-dtz)
    ELSE
      dttdz = dtdz
    END IF
    zht   = zht+dz
    press = p
    temp  = temp+dttdz*dz
    rl    = rlvo+(273.15-temp)*dlvdt
    arg   = rl*(temp-273.15)/273.15/temp/rv
    es2   = eso*EXP(arg)
    qvs2  = 0.622*es2/(press-es2)
!        rho2  = press/(rair*temp)
!        alw   = alw+(qvs1-qvs2)*(rho1+rho2)/2.   ! kg/m3

    alw   = alw+(qvs1-qvs2)                   ! kg/kg
    IF(alw < 0.0) alw=0.0
!
!-----------------------------------------------------------------------
!
!  Application of a simple linear glaciation
!c   ---------------------------------------------------------------
!c
!c         all liquid T > -15 C
!c         partially liquid -15 C > T > -25 C
!c         all ice    T < -25 C
!c
!-----------------------------------------------------------------------
!
    IF(cldtpt < 258.15) THEN
      IF(cldtpt > 248.15) THEN
        fraclw = 0.1*(cldtpt-248.15)
      ELSE
        fraclw = 0.0
      END IF
    ELSE
      fraclw = 1.0
    END IF
    calw(j+1) = alw*fraclw*1000.0
    WRITE(6,9015) j,calw(j+1),zht
    !9015   format(1x,'j=',i3,'  adiab.lwc =',f9.5,'  alt =',f8.0)
  END DO
  455   CONTINUE
!c
!-----------------------------------------------------------------------
!
!  Obtain profile of LWCs at the given grid point
!
!-----------------------------------------------------------------------
!
  DO ip=2,nz-1
    IF(zs_1d(ip) <= zcb .OR. zs_1d(ip) > zctop) THEN
      ql(ip)=0.0
      qi(ip)=0.0
      cldt(ip)=t_1d(ip)
    ELSE
      DO j=2,nlevel
        zcloud = zcb+(j-1)*dz
        IF(zcloud >= zs_1d(ip)) THEN
          ql(ip) = (zs_1d(ip)-zcloud+100.)*(calw(j)-calw(j-1))*0.01     &
                       +calw(j-1)
          qi(ip) = (zs_1d(ip)-zcloud+100.)*(cali(j)-cali(j-1))*0.01     &
                       +cali(j-1)
          tmpk = (zs_1d(ip)-zcloud+100.)*(catk(j)-catk(j-1))*0.01     &
                       +catk(j-1)
          entc = (zs_1d(ip)-zcloud+100.)*(entr(j)-entr(j-1))*0.01     &
                       +entr(j-1)
          cldt(ip) = (1.-entc)*t_1d(ip) + entc*tmpk

          IF(l_prt .AND. myproc == 0) THEN
            WRITE(6,*) ' Get_sfm_1d: assigning ql(ip),qi(ip)'
            WRITE(6,609) ip,zs_1d(ip),zcb,zctop
            609             FORMAT(1X,' ip=',i2,' ht=',f8.0,' zcb='     &
                                    ,f8.0,' zctop=',f8.0)
            WRITE(6,610) j,zcloud,j-1,zcloud-dz
            610             FORMAT(1X,' j=',i2,' z=',f8.0,' j-1=',i2,' z=',f8.0)
            WRITE(6,611) calw(j),calw(j-1),ql(ip)                       &
                       , cali(j),cali(j-1),qi(ip)
            611             FORMAT(1X,' lwc_j=',f7.3,' lwc_1=',f7.3,' ql=',f7.3, &
                                       ' ice_j=',f7.3,' ice_1=',f7.3,' qi=',f7.3)
          END IF
          EXIT
        END IF
      END DO
!      475      CONTINUE
    END IF
  END DO
!
!-----------------------------------------------------------------------
!
!  Write out file of lwc comparisons
!
!-----------------------------------------------------------------------
!
!  9001 FORMAT(i7)
!  9002 FORMAT(1X,3I2,1X,14F8.2,i2,i3)
!  9004 FORMAT(1X,2E15.8,'ihr, imin, isec=',3(i8,1X))
!  9005 FORMAT(2X,'Predicted LWC',8X,'Observed LWC',/)
!  9014 FORMAT(1X,8E15.8)
  RETURN
END SUBROUTINE get_sfm_1d
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE CLOUD_TYPE_QC               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE cloud_type_qc (nx,ny,nz,zs_3d,cloud_cvr_3d,                  &
                          cloud_type_3d)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  This routine checks and modifies the analyzed cloud type field
!  to assure the spatial consistency of cumulonimbus clouds.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (Jian Zhang)
!  08/96
!
!  MODIFICATION HISTORY:
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
  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!
!  INPUT:
  INTEGER :: nx,ny,nz
  REAL :: zs_3d(nx,ny,nz)        ! physical heights

!  INPUT/OUTPUT:
  INTEGER :: cloud_type_3d(nx,ny,nz)     ! cloud/precip type field
  REAL    :: cloud_cvr_3d(nx,ny,nz)      ! cloud cover analysis
!
!  LOCAL:
!
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,k1,k_cb_top,k_cb_base
  REAL    :: depth_cb,depth_cu,depth_sc
  REAL    :: cldcvr_cb,cldcvr_below,cldcvr_above
  LOGICAL :: cld_below_cb,cld_above_cb,l_prt
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (myproc == 0)                                                     &
    WRITE(6,*) 'Quality check consistency of comulunimbus cloud.'
!
!-----------------------------------------------------------------------
!
!  Quality check the horizontal consistency of Cb cloud
!
!-----------------------------------------------------------------------
!
  l_prt=.false.
  DO k = 1,nz-1
    DO j = 1,ny-1
      DO i = 2,nx-2
!          l_prt = j.eq.29.and.i.lt.15.and.k.lt.25
        IF(l_prt .AND. myproc == 0) THEN
          PRINT*,i,j,k,'iii00ctyp:',cloud_type_3d(i,j,k),               &
                 ' cvr:',cloud_cvr_3d(i,j,k),                           &
                 'ctym1:',cloud_type_3d(i-1,j,k),                       &
                 'ctyp1:',cloud_type_3d(i+1,j,k)
        END IF
        IF(cloud_type_3d(i,j,k) == 10 ) CYCLE
        IF(cloud_type_3d(i-1,j,k) == 10 .AND.cloud_type_3d(i+1,j,k) == 10) THEN
          cloud_type_3d(i,j,k) = 10
          IF(cloud_cvr_3d(i,j,k) < 0.1)                                 &
              cloud_cvr_3d(i,j,k) = 0.5*(cloud_cvr_3d(i-1,j,k)          &
                                       + cloud_cvr_3d(i+1,j,k))
        END IF
        IF(l_prt .AND. myproc == 0) THEN
          PRINT*,i,j,k,'iii11ctyp:',cloud_type_3d(i,j,k),               &
                 ' cvr:',cloud_cvr_3d(i,j,k),                           &
                 'ctym1:',cloud_type_3d(i-1,j,k),                       &
                 'ctyp1:',cloud_type_3d(i+1,j,k)
        END IF
      END DO
    END DO  !j

    l_prt=.false.
    DO j = 2,ny-2
      DO i = 1,nx-1
!          l_prt = j.eq.29.and.i.lt.15.and.k.lt.25
        IF(l_prt .AND. myproc == 0) THEN
          PRINT*,i,j,k,'jjj00ctyp:',cloud_type_3d(i,j,k),               &
                 ' cvr:',cloud_cvr_3d(i,j,k),                           &
                 'ctym1:',cloud_type_3d(i,j-1,k),                       &
                 'ctyp1:',cloud_type_3d(i,j+1,k)
        END IF
        IF(cloud_type_3d(i,j,k) == 10 ) CYCLE
        IF(cloud_type_3d(i,j-1,k) == 10 .AND.cloud_type_3d(i,j+1,k) == 10) THEN
          cloud_type_3d(i,j,k) = 10
          IF(cloud_cvr_3d(i,j,k) < 0.1)                                 &
              cloud_cvr_3d(i,j,k) = 0.5*(cloud_cvr_3d(i,j-1,k)          &
                                       + cloud_cvr_3d(i,j+1,k))
        END IF
        IF(l_prt .AND. myproc == 0) THEN
          PRINT*,i,j,k,'jjj11ctyp:',cloud_type_3d(i,j,k),               &
                 ' cvr:',cloud_cvr_3d(i,j,k),                           &
                 'ctym1:',cloud_type_3d(i,j-1,k),                       &
                 'ctyp1:',cloud_type_3d(i,j+1,k)
        END IF
      END DO
    END DO  !i
  END DO  !k
!
!-----------------------------------------------------------------------
!
!  Quality check the vertical consistency of Cb cloud
!
!-----------------------------------------------------------------------
!
  DO j = 1,ny-1
    DO i = 1,nx-1
      DO k = 3, nz-2
!          l_prt = j.eq.29.and.i.lt.15.and.k.lt.25
        cld_above_cb =.false.
        cld_below_cb =.false.
        cldcvr_above = 0.01
        cldcvr_below = 0.01
        IF(l_prt .AND. myproc == 0) THEN
          PRINT*,i,j,k,'kkk00ctyp:',cloud_type_3d(i,j,k),               &
                 'ctyp2:',cloud_type_3d(i,j,k+2),                       &
                 'ctyp1:',cloud_type_3d(i,j,k+1),                       &
                 ' ctyp:',cloud_type_3d(i,j,k),                         &
                 'ctym1:',cloud_type_3d(i,j,k-1),                       &
                 'ctym2:',cloud_type_3d(i,j,k-2)
          PRINT*,i,j,k,'kkk00cvr:',cloud_cvr_3d(i,j,k),                 &
                 'cvrp2:',cloud_cvr_3d(i,j,k+2),                        &
                 'cvrp1:',cloud_cvr_3d(i,j,k+1),                        &
                 ' ccvr:',cloud_cvr_3d(i,j,k),                          &
                 'cvrm1:',cloud_cvr_3d(i,j,k-1),                        &
                 'cvrm2:',cloud_cvr_3d(i,j,k-2)
        END IF
        IF(cloud_type_3d(i,j,k) == 3.OR.cloud_type_3d(i,j,k) == 10) CYCLE
        IF(cloud_type_3d(i,j,k+1) == 10) THEN
          cld_above_cb =.true.
          cldcvr_above = cloud_cvr_3d(i,j,k+1)
          GO TO 14
        ELSE IF(cloud_type_3d(i,j,k+2) == 10) THEN
          cld_above_cb =.true.
          cldcvr_above = cloud_cvr_3d(i,j,k+2)
        END IF

        14        CONTINUE
        IF(cloud_type_3d(i,j,k-1) == 10 .OR.cloud_type_3d(i,j,k-1) == 3) THEN
          cld_below_cb =.true.
          cldcvr_below = cloud_cvr_3d(i,j,k-1)
          GO TO 15
        ELSE IF (cloud_type_3d(i,j,k-2) == 10                           &
                  .OR. cloud_type_3d(i,j,k-2) == 3) THEN
          cld_below_cb =.true.
          cldcvr_below = cloud_cvr_3d(i,j,k-2)
        END IF

        15        CONTINUE
        IF(l_prt .AND. myproc == 0) THEN
          PRINT*,i,j,k,'kkk111:',cld_below_cb,cld_above_cb,             &
                 'cvr_below:',cldcvr_below,                             &
                 'cvr_above:',cldcvr_above
        END IF
        IF(cld_below_cb .AND. cld_above_cb) THEN
          cloud_type_3d(i,j,k) = 10
          IF(cloud_cvr_3d(i,j,k) < 0.1) THEN
            cloud_cvr_3d(i,j,k) = 0.5*(cldcvr_above+cldcvr_below)
          END IF
        END IF
        IF(l_prt .AND. myproc == 0) THEN
          PRINT*,i,j,k,'kkk22ctyp:',cloud_type_3d(i,j,k),               &
                 'ctyp2:',cloud_type_3d(i,j,k+2),                       &
                 'ctyp1:',cloud_type_3d(i,j,k+1),                       &
                 ' ctyp:',cloud_type_3d(i,j,k),                         &
                 'ctym1:',cloud_type_3d(i,j,k-1),                       &
                 'ctym2:',cloud_type_3d(i,j,k-2),                       &
                 'cvr:',cloud_cvr_3d(i,j,k)
          PRINT*,i,j,k,'kkk22cvr:',cloud_cvr_3d(i,j,k),                 &
                 'cvrp2:',cloud_cvr_3d(i,j,k+2),                        &
                 'cvrp1:',cloud_cvr_3d(i,j,k+1),                        &
                 ' ccvr:',cloud_cvr_3d(i,j,k),                          &
                 'cvrm1:',cloud_cvr_3d(i,j,k-1),                        &
                 'cvrm2:',cloud_cvr_3d(i,j,k-2)
        END IF
      END DO
    END DO  !j
  END DO  !i
!
!-----------------------------------------------------------------------
!
!  Final quality check the vertical consistency of Cb cloud
!
!-----------------------------------------------------------------------
!
  DO j = 1,ny-1
    DO i = 1,nx-1
!
      depth_cb = 0.0
      depth_cu = 0.0
      depth_sc = 0.0

      DO k = 2, nz-2
!          l_prt = j.eq.29.and.i.lt.15.and.k.lt.25
        IF(l_prt .AND. myproc == 0) THEN
          PRINT*,i,j,k,'fff00ctyp:',cloud_type_3d(i,j,k),               &
                 'ccvr:',cloud_cvr_3d(i,j,k)
        END IF
        IF(cloud_type_3d(i,j,k) == 10)                                  &
            depth_cb = depth_cb + 0.5*(zs_3d(i,j,k+1)-zs_3d(i,j,k-1))
        IF(cloud_type_3d(i,j,k) == 3)                                   &
            depth_cu = depth_cu + 0.5*(zs_3d(i,j,k+1)-zs_3d(i,j,k-1))
        IF(cloud_type_3d(i,j,k) == 2)                                   &
            depth_sc = depth_sc + 0.5*(zs_3d(i,j,k+1)-zs_3d(i,j,k-1))
      END DO  !k
!          l_prt = j.eq.29.and.i.lt.15
      IF(l_prt .AND. myproc == 0) THEN
        PRINT*,i,j,'fff11depth_cb:',depth_cb,                           &
               ' depth_Cu:',depth_cu,' depth_Sc:',depth_sc
      END IF

      IF(depth_cb >= 2500.0 .OR.depth_cb >= 1500.0.AND.depth_cu >= 5000.0 &
            .OR.depth_cu >= 7000.0                                      &
            .OR.depth_cb > 200.0.AND.(depth_cu+depth_sc) >= 9000.) THEN

        DO k = 2, nz-2
          IF(cloud_type_3d(i,j,k) /= 0) THEN
            k_cb_base = k
            GO TO 17
          END IF
        END DO
        17        CONTINUE
        DO k = nz-1,2,-1
          IF(cloud_type_3d(i,j,k) /= 0) THEN
            k_cb_top = k
            GO TO 18
          END IF
        END DO
        18        CONTINUE
!          l_prt = j.eq.29.and.i.lt.15
        IF(l_prt .AND. myproc == 0) THEN
          PRINT*,i,j,'fff22k_Cb_base:',k_cb_base,                       &
                 ' k_Cb_top:',k_cb_top
        END IF

        cldcvr_cb = cloud_cvr_3d(i,j,k_cb_base)
        DO k1 = k_cb_base, k_cb_top
          IF(l_prt .AND. myproc == 0) THEN
            PRINT*,i,j,k1,'fff33cldcvr_Cb_below:',cldcvr_cb,            &
                   ' cvr:',cloud_cvr_3d(i,j,k1),                        &
                   ' ctyp:',cloud_type_3d(i,j,k1)
          END IF
          cloud_type_3d(i,j,k1) = 10
          IF(cloud_cvr_3d(i,j,k1) < 0.1) THEN
            cloud_cvr_3d(i,j,k1) = cldcvr_cb
          ELSE
            cldcvr_cb = cloud_cvr_3d(i,j,k1)
          END IF
          IF(l_prt .AND. myproc == 0) THEN
            PRINT*,i,j,k1,'fff44cldcvr_Cb_below:',cldcvr_cb,            &
                   ' cvr:',cloud_cvr_3d(i,j,k1),                        &
                  ' ctyp:',cloud_type_3d(i,j,k1)
          END IF
        END DO
      END IF
    END DO  !i
  END DO  !j

  RETURN
END SUBROUTINE cloud_type_qc
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE CLOUD_W                       ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!##################################################################
!##################################################################
!

SUBROUTINE cloud_w (nz,cloud_type,zs_1d,w                               &
           ,r_missing, wmhr_cu,wmhr_sc,wc_st)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Calculate in-cloud w-field as a function of cloud thickness
!  and cloud type using prespecified profile shape.
!  (e.g., parabolic for cumulus and constant for stratus clouds).
!
!-----------------------------------------------------------------------
!
!  AUTHOR:  Jian Zhang
!  05/96    Based on the LAPS cloud analysis code of 07/1995
!
!  MODIFICATION HISTORY:
!  05/16/96 (Jian Zhang)
!           Modified for ADAS format. Added full documentation.
!  05/01/97 (Jian Zhang)
!           Change the way that the values of "wmhr"s
!           are defined. Previously their values are hard coded
!           in this subroutine. Now they are defined in
!           "adascld.inc" file.
!  08/05/97 (Jian Zhang)
!           Further changes to the way that the values of "wmhr"s
!           are defined. Now they are input arguments.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

!  Input:
  REAL :: r_missing   ! indicator for missing data
  REAL :: wmhr_cu      ! coef_wmax for Cu-clouds (m/s/m)
  REAL :: wmhr_sc      ! coef_wmax for Sc-clouds (m/s/m)
  REAL :: wc_st
!
!  INPUT:
  INTEGER :: nz             ! number of model vertical levels
  INTEGER :: cloud_type(nz) ! clour type index, see following table
!
!-----------------------------------------------------------------------
!
!   Cloud Type:  '  ','St','Sc','Cu','Ns','Ac','As','Cs','Ci','Cc','Cb'
!   Integer Value: 0    1    2    3    4    5    6    7    8    9   10
!
!-----------------------------------------------------------------------
  REAL :: zs_1d(nz)        ! phusical heights at each level
!
!  OUTPUT
  REAL :: w(nz)             ! estimated in-cloud w values at each level
!
!  LOCAL:
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  REAL :: ratio, wcld, parabolic_w_profile

  INTEGER :: k, k1, kbase, ktop
  REAL :: zbase, ztop
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  Zero out return vector
!
!-----------------------------------------------------------------------
!
  DO k = 1, nz
    w(k) = 0.
  END DO
!
!-----------------------------------------------------------------------
!
!  Put in the wcld's for cumuliform clouds (Cu or Cb) first.
!
!-----------------------------------------------------------------------
!
  ratio = wmhr_cu
  DO k = 1, nz-1
    IF (cloud_type(k) == 3  .OR.  cloud_type(k) == 10) THEN
      kbase = k
      GO TO 10
    END IF
  END DO
  GO TO 100

  10    DO k = kbase, nz-1
    IF (cloud_type(k) == 3  .OR.  cloud_type(k) == 10) THEN
      ktop = k
    ELSE
      GO TO 20
    END IF
  END DO

  20    k1 = k          ! save our place in the column
  zbase = zs_1d(kbase)
  ztop  = zs_1d(ktop)
  DO k = 1, nz-1
    wcld = parabolic_w_profile (zbase, ztop, ratio, zs_1d(k))
    IF (wcld > 0.) THEN
      w(k) = wcld
    ELSE
      w(k) = 0.
    END IF
  END DO

  k1 = k1 + 1
  IF (k1 >= nz-1) GO TO 100
!
!-----------------------------------------------------------------------
!
!  Try for another level of Cu.
!
!-----------------------------------------------------------------------
!
  DO k = k1, nz-1
    IF (cloud_type(k) == 3  .OR.  cloud_type(k) == 10) THEN
      kbase = k
      GO TO 10
    END IF
  END DO
!
!-----------------------------------------------------------------------
!
!  Now do the stratocumulus or similar clouds (Sc, Ac, Cc, Ns).
!
!-----------------------------------------------------------------------
!
  100   ratio = wmhr_sc
  DO k = 1, nz-1
    IF (cloud_type(k) == 2 .OR. cloud_type(k) == 4                      &
          .OR. cloud_type(k) == 5 .OR. cloud_type(k) == 9) THEN
      kbase = k
      GO TO 110
    END IF
  END DO
  GO TO 200

  110   DO k = kbase, nz-1
    IF (cloud_type(k) == 2 .OR. cloud_type(k) == 4                      &
          .OR. cloud_type(k) == 5 .OR. cloud_type(k) == 9) THEN
      ktop = k
    ELSE
      GO TO 120
    END IF
  END DO

  120   k1 = k          ! save our place in the column
  zbase = zs_1d(kbase)
  ztop  = zs_1d(ktop)
  DO k = 1, nz-1
    wcld = parabolic_w_profile (zbase, ztop, ratio, zs_1d(k))
    IF (wcld > w(k)) w(k) = wcld
  END DO
  k1 = k1 + 1
  IF (k1 >= nz) GO TO 200       ! try for stratiform clouds
!
!-----------------------------------------------------------------------
!
!  Try for another level of Sc.
!
!-----------------------------------------------------------------------
!
  DO k = k1, nz-1
    IF (cloud_type(k) == 2 .OR. cloud_type(k) == 4                      &
          .OR. cloud_type(k) == 5 .OR. cloud_type(k) == 9) THEN
      kbase = k
      GO TO 110
    END IF
  END DO

!
!-----------------------------------------------------------------------
!
!  Make sure there is non-zero wcld wherever there are clouds of
!  any kind.  Also, return missing-data value for any place where
!  wcld value is not diagnosed.
!
!-----------------------------------------------------------------------
!
  200   DO k = 1, nz-1
    IF(cloud_type(k) /= 0 .AND. w(k) < wc_st) w(k) = wc_st
    IF (w(k) == 0.) w(k) = r_missing
  END DO

  RETURN
END SUBROUTINE cloud_w






!
!##################################################################
!##################################################################
!######                                                      ######
!######         FUNCTION PARABOLIC_W_PROFILE                 ######
!######                                                      ######
!##################################################################
!##################################################################
!

  FUNCTION parabolic_w_profile (zbase, ztop, ratio, z)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Define a parabolic-shaped w-profile.
!  The vertical velocity is zero at cloud top, peaks one third of
!  the way up from the base, and extends below the base by one third
!  of the cloud depth.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!  07/95
!
!
!  MODIFICATION HISTORY:
!  05/16/96 (Jian Zhang)
!  Added ARPS format document
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
!  INPUT:
  REAL :: zbase ! height at cloud base
  REAL :: ztop  ! height at cloud top
  REAL :: ratio ! specifies the max. w in two cld typs as func of z
  REAL :: z     ! height at the given grid point
!
!  OUTPUT
  REAL :: parabolic_w_profile ! w value at the given height
!
!  LOCAL:
  REAL :: depth, vvmax, vvspan, halfspan, height_vvmax, x
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  depth = ztop - zbase
  IF (depth <= 0.) THEN
    parabolic_w_profile = 0.
    RETURN
  END IF

  vvmax = ratio * depth
  vvspan = depth * 4. / 3.
  halfspan = vvspan / 2.
  height_vvmax = ztop - halfspan
  x = -vvmax/(halfspan*halfspan)

  parabolic_w_profile = x * (z-height_vvmax)*(z-height_vvmax) + vvmax

  RETURN
  END FUNCTION parabolic_w_profile




!
!
!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE PCP_TYPE_3D                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE pcp_type_3d (nx,ny,nz,temp_3d,rh_3d,p_pa_3d,                 &
                        radar_3d,l_mask,cldpcp_type_3d,istatus)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  This routine returns 3D cloud and precipitation type field.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:  Jian Zhang
!  05/1996  Based on the LAPS cloud analysis code developed by
!           Steve Albers.
!
!  This program modifies the most significant 4 bits of the integer
!  array by inserting multiples of 16.
!
!  MODIFICATION HISTORY:
!
!  05/16/96 (J. Zhang)
!           Modified for ADAS format. Added full documentation.
!  01/20/98 (J. Zhang)
!           Fixed a bug that no precip. type was assigned for a
!           grid point at the top of the radar echo with Tw
!           falling in the range of 0 to 1.3 degree C.
!  01/21/98 (J. Zhang)
!           Fixed a bug that does the freezing/refreezing test
!           on ice precipitates.
!  02/17/98 (J. Zhang)
!           Change the hail diagnose procedure.
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
!  INPUT:
  INTEGER :: nx,ny,nz                  ! Model grid size
  REAL :: temp_3d(nx,ny,nz)            ! temperature (K)
  REAL :: rh_3d(nx,ny,nz)              ! relative humudity
  REAL :: p_pa_3d(nx,ny,nz)            ! pressure (Pascal)
  REAL :: radar_3d(nx,ny,nz)           ! radar refl. (dBZ)
!
!  OUTPUT:
  INTEGER :: istatus
  INTEGER :: cldpcp_type_3d(nx,ny,nz)! cld/precip type
  INTEGER :: itype                   ! cld/precip type index
  LOGICAL :: l_mask(nx,ny)             ! "Potential" Precip Type
!
!  LOCAL functions:
  REAL :: dwpt                         ! for dew point calcl'n
  REAL :: wb_melting_thres             ! define melting temp. thresh.
  REAL :: tw                           ! for wet-bulb temp calcl'n
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,k_upper
  REAL    :: t_c,td_c,t_wb_c,temp_lower_c,temp_upper_c,tbar_c,          &
             p_mb,thickns,frac_below_zero
  INTEGER :: iprecip_type,iprecip_type_last,iflag_melt,                 &
             iflag_refreez
  REAL    :: zero_c,rlayer_refreez_max,rlayer_refreez
  INTEGER :: n_zr,n_sl,n_last
  REAL    :: tmelt_c

  LOGICAL :: verbose
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

  verbose = .FALSE.
!
!-----------------------------------------------------------------------
!
!  Stuff precip type into cloud type array
!  0 - No Precip
!  1 - Rain
!  2 - Snow
!  3 - Freezing Rain
!  4 - Sleet
!  5 - Hail
!
!-----------------------------------------------------------------------
!
  zero_c = 273.15
  rlayer_refreez_max = 0.0

  n_zr = 0
  n_sl = 0
  n_last = 0

  DO j = 1, ny-1
    DO i = 1, nx-1

      iflag_melt = 0
      iflag_refreez = 0
      rlayer_refreez = 0.0

      iprecip_type_last = 0

      DO k = nz-1,1,-1

        IF(radar_3d(i,j,k) >= 0. .OR. l_mask(i,j)) THEN
!
!-----------------------------------------------------------------------
!
!  Set refreezing flag
!
!-----------------------------------------------------------------------
!
          t_c  = temp_3d(i,j,k) - zero_c
          td_c = dwpt(t_c,rh_3d(i,j,k))
          p_mb = 0.01*p_pa_3d(i,j,k)

          tmelt_c = wb_melting_thres(t_c,radar_3d(i,j,k))
          t_wb_c = tw(t_c,td_c,p_mb)

          IF(t_wb_c < 0.) THEN
            IF(iflag_melt == 1) THEN
!
!-----------------------------------------------------------------------
!
!  Integrate below freezing temperature times column thickness
!   - ONLY for portion of layer below freezing
!
!-----------------------------------------------------------------------
!
              temp_lower_c = t_wb_c
              k_upper = MIN(k+1,nz-1)
!
!-----------------------------------------------------------------------
!
!  For simplicity and efficiency, the assumption is here made that
!  the wet bulb depression is constant throughout the level.
!
!-----------------------------------------------------------------------
!
              temp_upper_c = t_wb_c + ( temp_3d(i,j,k_upper)            &
                                          - temp_3d(i,j,k))
              IF(temp_upper_c <= 0.) THEN
                frac_below_zero = 1.0
                tbar_c = 0.5 * (temp_lower_c + temp_upper_c)

              ELSE ! Layer straddles the freezing level
                frac_below_zero = temp_lower_c                          &
                                        / (temp_lower_c - temp_upper_c)
                tbar_c = 0.5 * temp_lower_c

              END IF

              thickns = p_pa_3d(i,j,k_upper) - p_pa_3d(i,j,k)
              rlayer_refreez = rlayer_refreez                           &
                   + ABS(tbar_c * thickns * frac_below_zero)

              IF(rlayer_refreez >= 25000.) THEN
                iflag_refreez = 1
              END IF

              rlayer_refreez_max =                                      &
                                MAX(rlayer_refreez_max,rlayer_refreez)

            END IF ! iflag_melt = 1

          ELSE ! Temp > 0C
            iflag_refreez = 0
            rlayer_refreez = 0.0

          END IF ! T < 0.0c, Temp is below freezing
!
!-----------------------------------------------------------------------
!
!  Set melting flag
!
!-----------------------------------------------------------------------
!
          IF(t_wb_c >= tmelt_c) THEN
            iflag_melt = 1
          END IF

          IF(t_wb_c >= tmelt_c) THEN  ! Melted to Rain
            iprecip_type = 1

          ELSE ! Check if below zero_c (Refrozen Precip or Snow)
            IF(t_wb_c < 0.0) THEN
              IF(iflag_melt == 1) THEN
                IF(iprecip_type_last == 1 .OR.iprecip_type_last == 3) THEN
                                   ! test if rain or zr freeze
                  IF(iflag_refreez == 0) THEN ! Freezing Rain
                    n_zr = n_zr + 1
                    IF(n_zr < 30) THEN
                      WRITE(6,5)i,j,k,t_wb_c,temp_3d(i,j,k)             &
                          ,rh_3d(i,j,k)
                      5 FORMAT('zr',3I3,2F8.2,f8.1)
                    END IF
                    iprecip_type = 3

                  ELSE  ! (iflag_refreez = 1)  ! Sleet
                    n_sl = n_sl + 1
                    iprecip_type = 4
                  END IF ! iflag_refreez .eq. 0
                ELSE
                  iprecip_type = iprecip_type_last  ! Unchanged
                  n_last = n_last + 1
                  IF(n_last < 5 .AND. verbose) THEN
                    WRITE(6,*)'Unchanged Precip',i,j,k,t_wb_c
                  END IF
                END IF      ! liquid precip. at upper level?

              ELSE    ! iflag_melt =0        ! Snow
                iprecip_type = 2

              END IF   ! iflag_melt = 1?
            ELSE ! t_wb_c >= 0c, and t_wb_c < tmelt_c

              IF (iprecip_type_last == 0) THEN        !   1/20/98
                iprecip_type = 1    ! rain:at echo top and 0<Tw<1.3C
                iflag_melt = 1
              ELSE
                iprecip_type = iprecip_type_last
                n_last = n_last + 1
                IF(n_last < 5 .AND. verbose) THEN
                  WRITE(6,*)'Unchanged Precip',i,j,k,t_wb_c
                END IF
              END IF

            END IF ! t_wb_c < 0c
          END IF  ! t_wb_c >= tmelt_c

        ELSE ! radar_3d < 0dBZ;  No Radar Echo
          iprecip_type = 0
          iflag_melt = 0
          iflag_refreez = 0
          rlayer_refreez = 0.0

        END IF ! radar_3d(i,j,k).ge.0. .or. l_mask(i,j);  Radar Echo?
!
!-----------------------------------------------------------------------
!
!  Insert most sig 4 bits into array
!
!-----------------------------------------------------------------------
!
        itype = cldpcp_type_3d(i,j,k)
        itype = itype - (itype/16)*16     ! Initialize the 4 bits
        itype = itype + iprecip_type * 16 ! Add in the new value
        cldpcp_type_3d(i,j,k) = itype

        iprecip_type_last = iprecip_type

      END DO ! k
    END DO ! j
  END DO ! i

  DO j = 1, ny-1
    DO i = 1, nx-1
      DO k = 1,nz-1
        IF(radar_3d(i,j,k) >= 50.) THEN
          iprecip_type = 5
          itype = cldpcp_type_3d(i,j,k)
          itype = itype - (itype/16)*16     ! Initialize the 4 bits
          itype = itype + iprecip_type * 16 ! Add in the new value
          cldpcp_type_3d(i,j,k) = itype
        END IF
      END DO ! k
    END DO ! j
  END DO ! i

  IF (myproc == 0) THEN
    WRITE(6,*)' rlayer_refreez_max = ',rlayer_refreez_max
    WRITE(6,*)' n_frz_rain/n_sleet = ',n_zr,n_sl
  END IF
  istatus=1

  RETURN
END SUBROUTINE pcp_type_3d
!
!##################################################################
!##################################################################
!######                                                      ######
!######            FUNCTION WB_MELTING_THRES                 ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

  FUNCTION wb_melting_thres(t_c,dbz)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  This function calculates the wet-bulb threshold for melting snow
!  into rain as a function of dbz and t_c.
!
!  Currently it's simply set to a constant. Later it may be defined
!  as function of temperature and radar echo intensity.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!  07/95
!
!  MODIFICATION HISTORY:
!  05/17/96 (Jian Zhang)
!
!-----------------------------------------------------------------------

  wb_melting_thres = 1.3  ! Units are C

  RETURN
  END FUNCTION wb_melting_thres



!
!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE ISI3                         ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE isi3 (slw,temp,ictype,iptype,xindex)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  This subroutine is to calculate icing severity index
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (M. POLITOVICH, RAP, NCAR)
!  30 Dec 1991
!
!  MODIFICATION HISTORY:
!
!  05/16/96  (Jian Zhang)
!            Modified for ADAS format. Added full documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT:  slw     liquid water content in g/m3
!          temp    temperature in deg C
!          ictype  cloud type
!          iptype  precip type
!  OUTPUT: index   severity index - real number from 0-6
!
!  NOTE: there are really 3 categories, with
!         continuous or intermittent designation
!        1-3 = cats 1, 2 and 3 continuous
!        4-6 = cats 1, 2 and 3 intermittent
!        continuous:   cloud types st,as,cs
!        intermittent: cloud types cu,cs,cb,ac,sc
!
!  precip type: if zr, category 3 is designated
!
!  scat and tcat are arrays of thresholds for slw and temp
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  PARAMETER (islw=3, itemp=3, ityp=2)
  DIMENSION scat(islw),tcat(itemp)       ! lw and temp categories
  DATA scat /0.01,0.1,0.5/
  DATA tcat /-10.,-5.,0./
!
  DIMENSION iaray(islw+1, itemp+1, ityp) ! 3D severity index matrix
!
  DATA iaray/  0, 1, 2, 3, 0, 2, 2, 3,                                  &
               0, 2, 3, 3, 0, 0, 0, 0,                                  &
               0, 4, 5, 6, 0, 5, 5, 6,                                  &
               0, 5, 6, 6, 0, 0, 0, 0/
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  Sort input values
!  **NOTE that - test is for "ge" rather than "gt"
!
!-----------------------------------------------------------------------
!
  i = 1
  DO ii = 1, islw
    IF (slw > scat(ii)) i = ii + 1
  END DO

  j = 1
  DO ii = 1, itemp
    IF (temp > tcat(ii)) j = ii+1
  END DO
!
!-----------------------------------------------------------------------
!
!  Cloud type section (continuous/intermittent)
!   Types: 0/none 1/St 2/Sc 3/Cu 4/Ns 5/Ac
!          6/As   7/Cs 8/Ci 9/Cc 10/Cb
!   Continuous, k=1
!   Intermittent k=2
!  ** NOTE that - if there is no cloud type (=0) but there
!   IS slw calculated, the slw overrides**
!
!-----------------------------------------------------------------------
!
  k = 1
  IF (ictype == 2) k = 2
  IF (ictype == 3) k = 2
  IF (ictype == 5) k = 2
  IF (ictype >= 9) k = 2
!
!-----------------------------------------------------------------------
!
!  Precip type section (freezing rain/drizzle=cat 3)
!   Types: 0/none 1/rn 2/sn 3/zr 4/sl 5/ha 6-10/no assignment
!  ** NOTE that - freezing rain sets slw to category 4,
!   and, sets temperature to category 3 (even if it is
!   somehow diagnosed to "above freezing"
!   Also - if no cloud type given (eg, below cloud), zr
!    is considered "continuous"***
!
!-----------------------------------------------------------------------
!
  IF (iptype == 3) j = 3
  IF (iptype == 3) i = 4
!
!     write (6,5995) slw,temp,type
!  15995  FORMAT (' l,t,ty:',2F7.2, a2)
!     write (6,5996) i,j,k
!  15996  FORMAT (' i,j,k:',3I3)
!
!-----------------------------------------------------------------------
!
!  assign severity index from sorted values
!
!-----------------------------------------------------------------------
!
  INDEX = iaray(i,j,k)
  xindex = INDEX
!
  RETURN
END SUBROUTINE isi3
