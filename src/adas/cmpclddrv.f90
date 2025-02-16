SUBROUTINE cmpclddrv(nx,ny,nz,i4timanx,                                 &
           xs,ys,zs,j3,hterain,latgr,longr,                             &
           nobsng,indexsng,stnsng,isrcsng,csrcsng,xsng,ysng,            &
           nxlg,nylg,xslg,yslg,                                         &
           timesng,latsng,lonsng,hgtsng,                                &
           kloud,store_amt,store_hgt,                                   &
           stnrad,isrcrad,latrad,lonrad,elvrad,ixrad,jyrad,             &
           pprt,ptprt,qv,qscalar,w,                                     &
           pbar,ptbar,qvbar,rhostr,                                     &
           istat_radar,ref_mos_3d,w_cld,cldpcp_type_3d,                 &
           icing_index_3d,l_mask_pcptype,                               &
           p_3d,t_3d,rh_3d,qvprt_old,qw_old,                            &
           qr_cld,qs_cld,qh_cld,rho,tem1,tem2)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Driver for complex cloud analysis code.
!  This subroutine serves as an interface between the ADAS
!  code and the cloud code.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!  Keith Brewster, CAPS
!  Based on code originally installed directly into adas.f
!  by Jian Zhang, CAPS
!  February, 1998
!
!  MODIFICATION HISTORY:
!  05/05/1998 Jian Zhang
!  Abandoned the cloud grid, using the ARPS grid instead.
!
!  08/30/2000 D.H. Wang
!  Added option for Ferrier formulation of cloud water.
!
!  11/02/2000 K. Brewster
!  Modifications to allow tiny, insignificant differences in
!  mapping parameters between radar and ADAS.  Minor clean-up.
!
!  02/05/2001 Richard Carpenter, WDT
!  Precipitate determined from radar enhances, rather than replaces, the
!  background field.
!
!-----------------------------------------------------------------------
!  July, 2001 (K. Brewster)
!  Changed wrtvar calls to wrtvar1 to match new routine for writing
!  3-D arrays to be read by arpsplt as an arbitary array.
!
!  Some clean-up of comments and lines that had been commented-out.
!
!  Added new latent heating options, cldptopt = 4 and 5.
!  cldptop = 4 replaces temperature with cloud parcel temperature
!  (with entrainment) everywhere there is cloud.
!  cldptopt =5 replaces adusts temperature to cloud parcel temperature
!  beginning where w=-0.2 m/s ramping to full application of cloud parcel
!  temperature where w>=0.
!
!  11/01/2001 (K. Brewster)
!  Restructured the building of the reflectivity mosaic and moved
!  it into new subroutine, refmosaic, replacing rad_patch_fill.
!
!  02/07/2006 (K. W. Thomas)
!  MPI update.
!
!  04/27/2010 (K. Brewster)
!  Restructured handling of ref_mos_3d to make it an input to this
!  routine, avoid the necessity of re-reading the radar files.
!
!-----------------------------------------------------------------------
!
  USE my3mom_fncs_mod, ONLY: gammaDP

  IMPLICIT NONE
  INCLUDE 'adas.inc'
  INCLUDE 'bndry.inc'
  INCLUDE 'mp.inc'
  INCLUDE 'globcst.inc'

  INTEGER :: i4timanx        ! analysis time in seconds from 00:00 UTC
                             ! Jan. 1, 1960
!
!-----------------------------------------------------------------------
!
!  Grid variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz
!
  REAL :: xs(nx)
  REAL :: ys(ny)
  REAL :: zs(nx,ny,nz)
  REAL :: j3(nx,ny,nz)
  REAL :: hterain(nx,ny)
  REAL :: latgr(nx,ny)
  REAL :: longr(nx,ny)

  REAL :: pprt(nx,ny,nz)
  REAL :: ptprt(nx,ny,nz)
  REAL :: w(nx,ny,nz)
  REAL :: qv(nx,ny,nz)
  REAL :: qscalar(nx,ny,nz,nscalar)
  REAL :: pbar(nx,ny,nz)
  REAL :: ptbar(nx,ny,nz)
  REAL :: qvbar(nx,ny,nz)
  REAL :: rhostr(nx,ny,nz)

  INTEGER :: nxlg, nylg
  REAL ::  xslg(nxlg)
  REAL ::  yslg(nylg)

!
!-----------------------------------------------------------------------
!
!  Single-level data
!
!-----------------------------------------------------------------------
!
  INTEGER :: nobsng
  INTEGER :: indexsng(nobsng)
  REAL :: latsng(mx_sng),lonsng(mx_sng)
  REAL :: hgtsng(mx_sng)
  REAL :: xsng(mx_sng),ysng(mx_sng)
  INTEGER :: timesng(mx_sng)
  CHARACTER (LEN=5) :: stnsng(mx_sng)
  INTEGER :: isrcsng(mx_sng)
  CHARACTER (LEN=8) :: csrcsng(mx_sng)
  CHARACTER (LEN=4) :: store_amt(mx_sng,5)
  INTEGER :: kloud(mx_sng)
  REAL :: store_hgt(mx_sng,5)
!
!-----------------------------------------------------------------------
!
!  Radar site variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=5) :: stnrad(mx_rad)
  INTEGER :: isrcrad(0:mx_rad)
  REAL :: latrad(mx_rad),lonrad(mx_rad)
  REAL :: elvrad(mx_rad)
  INTEGER :: ixrad(mx_rad),jyrad(mx_rad)
!
!-----------------------------------------------------------------------
!
!  3D radar array
!
!-----------------------------------------------------------------------
!
  INTEGER :: istat_radar
  REAL :: ref_mos_3d(nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  Output variables
!
!-----------------------------------------------------------------------
!
  REAL :: clouds_3d (nx,ny,nz)  ! final 3D fractnl cloud cover analysis.

  REAL :: cloud_ceiling (nx,ny) ! cloud ceiling heights( m AGL)
  REAL :: vil (nx,ny) ! vertically intregrated cloud liquid/ice water

  REAL :: slwc_3d (nx,ny,nz)    ! cloud liquid water content (kg/kg)
  REAL :: cice_3d (nx,ny,nz)    ! cloud ice content (kg/kg)
  REAL :: ctmp_3d (nx,ny,nz)    ! cloud temperature (K)
  REAL :: w_cld (nx,ny,nz)      ! vertical velocity (m/s) in clouds

  INTEGER :: cldpcp_type_3d(nx,ny,nz) ! cloud/precip type field
  INTEGER :: icing_index_3d(nx,ny,nz) ! icing severity index
  LOGICAL :: l_mask_pcptype(nx,ny)    ! analyze precip type using simulated
                                      ! radar data?
!
!-----------------------------------------------------------------------
!
!  Temporary arrays
!
!-----------------------------------------------------------------------
!
  REAL :: p_3d     (nx,ny,nz)
  REAL :: t_3d     (nx,ny,nz)
  REAL :: rh_3d    (nx,ny,nz)
  REAL :: qvprt_old(nx,ny,nz)
  REAL :: qw_old   (nx,ny,nz)
  REAL :: qr_cld   (nx,ny,nz)
  REAL :: qs_cld   (nx,ny,nz)
  REAL :: qh_cld   (nx,ny,nz)
  REAL :: rho      (nx,ny,nz)
  REAL :: tem1     (nx,ny,nz)
  REAL :: tem2     (nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  REAL :: f_qvsat
  REAL, PARAMETER :: rhsubsat=0.950
  REAL, PARAMETER :: smfct1=0.5

  REAL, PARAMETER :: epsdif=0.0001

!
  CHARACTER (LEN=6) :: varid
  CHARACTER (LEN=20) :: varname
  CHARACTER (LEN=20) :: varunits
  INTEGER :: i,j,k,nq,irad
  INTEGER :: istatus_cld,istatus_lwc,istatus_pcp,istatwrt
  REAL :: xrad,yrad,qvsat,radri,radrj,p0inv,ppi,qpcp,rh2qv
  REAL :: tgrid,dqv_prt,dqw,qw_new,arg,rh,rhobar,wratio,ptdiff,ptcld
!
  INTEGER :: strhopt_rad                  ! streching option
  INTEGER :: mapproj_rad                  ! map projection indicator
  REAL :: dx_rad,dy_rad,dz_rad,dzmin_rad  ! grid spcngs
  REAL :: ctrlat_rad,ctrlon_rad           ! central lat and lon
  REAL :: tlat1_rad,tlat2_rad,tlon_rad    ! true lat and lon
  REAL :: scale_rad                       ! map scale factor
  REAL :: max_pt_adj,extrnz,extrnz1,extr1,relh,frac
  REAL :: min_pt_adj
  REAL :: cldpcpfrc,onemcpf,cldtot,pcptot,cldlim
  CHARACTER (LEN=72) :: warn_string

! constant in equation for radar reflectivity
  REAL :: Gr, Gi, Gs, Gg, Gh

  REAL, PARAMETER :: pi = 3.14159265


!
!-----------------------------------------------------------------------
!
!  Include file
!
!-----------------------------------------------------------------------
!
!  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'grid.inc'        ! Grid parameters
!
!-----------------------------------------------------------------------
!
!  Specify the scheme and products wanted.
!
!-----------------------------------------------------------------------
!
  INTEGER :: iflag_slwc
  LOGICAL :: l_flag_incld_w ,l_flag_cldtyp,l_flag_icing
!
!-----------------------------------------------------------------------
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of the cloud analysis procedure.
!  Some initializations.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  iflag_slwc = 13 ! New Smith-Fedds model for Liquid Water Content
  l_flag_incld_w = .true.  ! Analyze in-cloud w-field
  l_flag_cldtyp = .true.   ! Diagnose cloud type field
  l_flag_icing = .true.  ! Analyzing icing severity
!
!-----------------------------------------------------------------------
!
!  Compute state variables on grid including k=1, nz-1, nz, via
!  linear extrapolation, except RH is taken to be constant.
!  Save the old buoyancy field.
!
!-----------------------------------------------------------------------
!
! OpenMP changed loop order to j,k,i:
!$OMP PARALLEL DO PRIVATE(i,j,k,qvsat)
  DO j=1,ny-1
    DO k=2,nz-1
      DO i=1,nx-1
        p_3d(i,j,k) = pprt(i,j,k)+pbar(i,j,k)
        t_3d(i,j,k) = (ptprt(i,j,k)+ptbar(i,j,k))*                      &
                      ((pprt(i,j,k)+pbar(i,j,k))/p0)**rddcp
        qvsat=f_qvsat(p_3d(i,j,k),t_3d(i,j,k))
        rh_3d(i,j,k)=100.*MIN(1.,MAX(0.,(qv(i,j,k)/qvsat)))
        qvprt_old(i,j,k) = qv(i,j,k) - qvbar(i,j,k)
        qw_old(i,j,k) = 0.0
        DO nq=1,nscalarq
          qw_old(i,j,k) = qw_old(i,j,k) + qscalar(i,j,k,nq)
        END DO
      END DO
    END DO
  END DO

! OpenMP:
!$OMP PARALLEL DO PRIVATE(j,i,extrnz1,extrnz,extr1)
  DO j = 1, ny-1
    DO i = 1, nx-1
      extrnz1 = (zs(i,j,nz-1)-zs(i,j,nz-3))                             &
             / (zs(i,j,nz-2)-zs(i,j,nz-3))
      t_3d(i,j,nz-1) = t_3d(i,j,nz-2)                                   &
                       + (t_3d(i,j,nz-2)-t_3d(i,j,nz-3))*extrnz1
      p_3d(i,j,nz-1) = p_3d(i,j,nz-2)                                   &
                       + (p_3d(i,j,nz-2)-p_3d(i,j,nz-3))*extrnz1
      rh_3d(i,j,nz-1) = rh_3d(i,j,nz-2)

      extrnz = (zs(i,j,nz)-zs(i,j,nz-3))                                &
             / (zs(i,j,nz-2)-zs(i,j,nz-3))
      t_3d(i,j,nz) = t_3d(i,j,nz-3)                                     &
                       + (t_3d(i,j,nz-2)-t_3d(i,j,nz-3))*extrnz
      p_3d(i,j,nz) = p_3d(i,j,nz-3)                                     &
                       + (p_3d(i,j,nz-2)-p_3d(i,j,nz-3))*extrnz
      rh_3d(i,j,nz) = rh_3d(i,j,nz-2)

      extr1 = (zs(i,j,1)-zs(i,j,3))                                     &
             / (zs(i,j,2)-zs(i,j,3))
      t_3d(i,j,1) = t_3d(i,j,3)                                         &
                       + (t_3d(i,j,2)-t_3d(i,j,3))*extr1
      p_3d(i,j,1) = p_3d(i,j,3)                                         &
                       + (p_3d(i,j,2)-p_3d(i,j,3))*extr1
      rh_3d(i,j,1) = rh_3d(i,j,2)
    END DO
  END DO

  CALL edgfill( t_3d ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz)
  CALL edgfill( p_3d ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz)
  CALL edgfill(rh_3d ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz)
  CALL edgfill(   zs ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz)
!
!-----------------------------------------------------------------------
!
!  Initialize the 3D radar reflectivity array.
!  Reorganized loops to eliminate IF's inside loop.
!
!-----------------------------------------------------------------------
!
   qr_cld = 0.0
   qs_cld = 0.0
   qh_cld = 0.0

  IF( P_QR > 0 ) THEN
! OpenMP changed loop order to j,k,i:
!$OMP PARALLEL DO PRIVATE(j,i,k)
    DO j=1,ny
      DO k=1,nz
        DO i=1,nx
          qr_cld(i,j,k)=1000.0*qscalar(i,j,k,P_QR)
        END DO
      END DO
    END DO
  END IF
!
  IF (P_QS > 0) THEN
!$OMP PARALLEL DO PRIVATE(j,i,k)
    DO j=1,ny
      DO k=1,nz
        DO i=1,nx
          qs_cld(i,j,k)=1000.0*qscalar(i,j,k,P_QS)
        END DO
      END DO
    END DO
  END IF
!
  IF (P_QG > 0) THEN
!$OMP PARALLEL DO PRIVATE(j,i,k)
    DO j=1,ny
      DO k=1,nz
        DO i=1,nx
          qh_cld(i,j,k)=1000.0*qscalar(i,j,k,P_QG)
        END DO
      END DO
    END DO
  END IF
!
  IF (P_QH > 0) THEN 
!$OMP PARALLEL DO PRIVATE(j,i,k)
    DO j=1,ny
      DO k=1,nz
        DO i=1,nx
          qh_cld(i,j,k)=1000.0*qscalar(i,j,k,P_QH)
        END DO
      END DO
    END DO
  END IF
!
!-----------------------------------------------------------------------
!
!  Analyze 3D fractional cloud cover field by using
!  SAO, radar reflectivity, IR and VIS satellite data.
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) WRITE(6,*) ' Calling Cloud Cover Analysis: cloud_cv'
!
  CALL cloud_cv (nx,ny,nz,i4timanx,curtim,dirname,runname(1:lfnkey),    &
            hdmpfmt,hdfcompr,lvldbg,                                    &
            xs,ys,zs,dx,dy,hterain,latgr,longr,                         &
            nxlg,nylg,xslg,yslg,                                        &
            p_3d,t_3d,rh_3d,                                            &
            nobsng,indexsng,stnsng,isrcsng,csrcsng,xsng,ysng,           &
            timesng,latsng,lonsng,hgtsng,                               &
            kloud,store_amt,store_hgt,                                  &
            ref_mos_3d,istat_radar,                                     &
            clouds_3d,cloud_ceiling,tem1,tem2,                          &
            istatus_cld)

  IF(istatus_cld /= 1) THEN
    WRITE(6,*) 'Bad status returned from cloud_cv, Aborting...'
    GO TO 999
  END IF

  IF (mp_opt > 0) THEN    ! Since clouds_3d will be used later
    CALL acct_interrupt(mp_acct)
    CALL mpsendrecv2dew(clouds_3d,nx,ny,nz,ebc,wbc,0,tem1)
    CALL mpsendrecv2dns(clouds_3d,nx,ny,nz,nbc,sbc,0,tem1)
  END IF

!
!-----------------------------------------------------------------------
!
!  Analyze 3D cloud liquid/ice water field.
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) WRITE(6,*) ' Calculating cloud liquid/ice water.'
!
!-----------------------------------------------------------------------
!
!  Analyzing cloud liquid water field
!
!-----------------------------------------------------------------------
!
  CALL cloud_lwc (nx,ny,nz,curtim,dirname,runname(1:lfnkey),            &
                  clouds_3d,                                            &
                  t_3d,rh_3d,p_3d,zs,                                   &
                  istat_radar,ref_mos_3d,                               &
                  iflag_slwc,slwc_3d,cice_3d,ctmp_3d,                   &
                  l_flag_incld_w,w_cld,                                 &
                  l_mask_pcptype,l_flag_cldtyp,cldpcp_type_3d,          &
                  l_flag_icing,icing_index_3d,tem1,                     &
                  istatus_lwc)

  IF(istatus_lwc /= 1) THEN
    IF (myproc == 0) WRITE(6,*)' Bad status returned from cloud_lwc'
    GO TO 999
  END IF

  IF (mp_opt > 0) THEN    ! Note the computation in subroutine cloud_type_qc
    CALL acct_interrupt(mp_acct)
    CALL mpsendrecv2dew(clouds_3d,nx,ny,nz,ebc,wbc,0,tem1)
    CALL mpsendrecv2dns(clouds_3d,nx,ny,nz,nbc,sbc,0,tem1)
  END IF

! OK TO HERE without radars.
!
!-----------------------------------------------------------------------
!
!  Calculate 3D Precipitation mixing ratio in g/kg.
!  Note that qr_cld, qs_cld, and qh_cld are diagnosed
!  qr, qs and qh in g/kg, respectively.
!
!-----------------------------------------------------------------------
!
  IF(istat_radar == 1) THEN
    IF (cldqropt == 1) THEN
      !
      ! Kessler's scheme
      !
      IF (myproc == 0) THEN
        WRITE(6,'(a)') ' Computing Precip mixing ratio.'
        WRITE(6,'(a)') ' Using Kessler radar reflectivity equations...'
      END IF
      CALL pcp_mxr (nx,ny,nz,t_3d,p_3d,zs,hterain,ref_mos_3d,           &
                    cldpcp_type_3d,                                     &
                    qr_cld,qs_cld,qh_cld,                               &
                    cloudopt,refthr1,refthr2,hgtrefthr,istatus_pcp)

    ELSE IF (cldqropt == 2) THEN
      !
      ! Ferrier's scheme
      !
      IF (myproc == 0) THEN
        WRITE(6,'(a)') ' Computing Precip mixing ratio.'
        WRITE(6,'(a)') ' Using Ferrier radar reflectivity equations...'
      END IF
      CALL pcp_mxr_ferrier (nx,ny,nz,t_3d,p_3d,zs,hterain,ref_mos_3d,   &
                   cldpcp_type_3d,                                      &
                   qr_cld,qs_cld,qh_cld,                                &
                   cloudopt,refthr1,refthr2,hgtrefthr,istatus_pcp)
    END IF   !cldqropt=1 or 2
  END IF
!
!-----------------------------------------------------------------------
!
!  Calculate the vertically integrated condensates
!  (in unit of kg/kg)
!
!-----------------------------------------------------------------------
!
! OpenMP:
!$OMP PARALLEL DO PRIVATE(j,i,k,arg,rhobar)
  DO j=1,ny
    DO i=1,nx
      vil(i,j) = 0.0
      DO k=2,nz-1
        arg = slwc_3d(i,j,k)+cice_3d(i,j,k)                             &
               +qr_cld(i,j,k)+qs_cld(i,j,k)+qh_cld(i,j,k)   ! g/kg
        rhobar = rhostr(i,j,k)/j3(i,j,k)
        rho(i,j,k) = (pbar(i,j,k)+pprt(i,j,k))/(rd*t_3d(i,j,k))
        arg = 0.001* 0.5*(zs(i,j,k+1)-zs(i,j,k-1))*arg /rhobar  ! kg/m**2
        vil(i,j) = vil(i,j) +arg
      END DO
    END DO
  END DO
  IF(cld_files == 1) THEN
    varid='colvil'
    varname='Cloud VIL'
    varunits='kg/m**2'
    CALL wrtvar1(nx,ny,1,vil,varid,varname,varunits,                    &
                curtim,runname(1:lfnkey),dirname,istatwrt)
  END IF
  IF (myproc == 0) THEN
    WRITE(6,'(a)')' Cloud options: '
    WRITE(6,'(a,i3)') ' cldqvopt=',cldqvopt
    WRITE(6,'(a,i3)') ' cldqcopt=',cldqcopt
    WRITE(6,'(a,i3)') ' cldqropt=',cldqropt
    WRITE(6,'(a,i3)') ' cldptopt=',cldptopt
    WRITE(6,'(a,i3)') ' cldwopt =',cldwopt
  END IF
!
!-----------------------------------------------------------------------
!
!  Enhance the cloud liquid/ice water mixing ratio fields (convert
!   g/kg to kg/kg) and THEN DO smoothing.
!
!-----------------------------------------------------------------------
!
  IF (cldqcopt == 1) THEN
    IF (myproc == 0) WRITE(6,'(a)')' Enhancing qc and qi-fields'
    DO k=2,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
!         IF (qc(i,j,k) < 0.001*slwc_3d(i,j,k)) qc(i,j,k) = 0.001*slwc_3d(i,j,k)
!         IF (qi(i,j,k) < 0.001*cice_3d(i,j,k)) qi(i,j,k) = 0.001*cice_3d(i,j,k)
!          qc(i,j,k) = max( 0.0, 0.001*slwc_3d(i,j,k))
!          qi(i,j,k) = max( 0.0, 0.001*cice_3d(i,j,k))
           IF (P_QC > 0) qscalar(i,j,k,P_QC) = max( 0.0, 0.001*slwc_3d(i,j,k))
           IF (P_QI > 0) qscalar(i,j,k,P_QI) = max( 0.0, 0.001*cice_3d(i,j,k))
! Limit total cloud water plus ice to local qvsat.
          qvsat=f_qvsat( p_3d(i,j,k), t_3d(i,j,k))
          qvsat = qvslimit_2_qc *qvsat
!          arg = qc(i,j,k)+qi(i,j,k)
          arg = 0.0
          IF (P_QC > 0) arg = arg + qscalar(i,j,k,P_QC)
          IF (P_QI > 0) arg = arg + qscalar(i,j,k,P_QI)
          IF (arg > 1.0E-10 .AND. arg > qvsat) THEN
!            relh = qc(i,j,k)/arg
!            qc(i,j,k) = relh*qvsat
!            qi(i,j,k) = (1.0-relh)*qvsat
           relh = qscalar(i,j,k,P_QC)/arg
           IF (P_QC > 0) qscalar(i,j,k,P_QC) = relh*qvsat
           IF (P_QI > 0) qscalar(i,j,k,P_QI) = (1.0-relh)*qvsat
          END IF
          IF(mphyopt >= 9) THEN
            qscalar(i,j,k,P_NC) = ntcloud
            qscalar(i,j,k,P_NI) = 5.*exp(0.304*(273.15-max(233.,t_3d(i,j,k))))
          END IF
          IF(mphyopt == 11) THEN
            Gi = ((6.+alphaice)*(5+alphaice)*(4+alphaice))/((3.+alphaice)*(2+alphaice)*(1+alphaice))
            IF(qscalar(i,j,k,P_NI) /= 0.0) THEN
              qscalar(i,j,k,P_ZI) = Gi/((440.0)**2.)*((rho(i,j,k)*qscalar(i,j,k,P_QI))**2.)/qscalar(i,j,k,P_NI)
            ELSE
              qscalar(i,j,k,P_ZI) = 0.0
            END IF
          END IF
        END DO
      END DO
    END DO

    IF (mp_opt > 0) THEN  ! Before smoothing, the arrays must be MPI valid
      CALL acct_interrupt(mp_acct)
      IF (P_QC > 0) THEN
        CALL mpsendrecv2dew(qscalar(1,1,1,P_QC),nx,ny,nz,ebc,wbc,0,tem1)
        CALL mpsendrecv2dns(qscalar(1,1,1,P_QC),nx,ny,nz,nbc,sbc,0,tem1)
        IF(mphyopt >= 9) THEN
          CALL mpsendrecv2dew(qscalar(1,1,1,P_NC),nx,ny,nz,ebc,wbc,0,tem1)
          CALL mpsendrecv2dns(qscalar(1,1,1,P_NC),nx,ny,nz,nbc,sbc,0,tem1)
        END IF
      END IF
      IF (P_QI > 0) THEN
        CALL mpsendrecv2dew(qscalar(1,1,1,P_QI),nx,ny,nz,ebc,wbc,0,tem1)
        CALL mpsendrecv2dns(qscalar(1,1,1,P_QI),nx,ny,nz,nbc,sbc,0,tem1)
        IF(mphyopt >= 9) THEN
          CALL mpsendrecv2dew(qscalar(1,1,1,P_NI),nx,ny,nz,ebc,wbc,0,tem1)
          CALL mpsendrecv2dns(qscalar(1,1,1,P_NI),nx,ny,nz,nbc,sbc,0,tem1)
        END IF
        IF(mphyopt == 11) THEN
          CALL mpsendrecv2dew(qscalar(1,1,1,P_ZI),nx,ny,nz,ebc,wbc,0,tem1)
          CALL mpsendrecv2dns(qscalar(1,1,1,P_ZI),nx,ny,nz,nbc,sbc,0,tem1)
        END IF
      END IF
    END IF

    IF (smth_opt == 1) THEN
      DO k=2,nz-1
        IF (P_QC > 0) CALL smooth9p(qscalar(1,1,k,P_QC),nx,ny,1,nx-1,1,ny-1,0,tem1(1,1,1))
        IF (P_QI > 0) CALL smooth9p(qscalar(1,1,k,P_QI),nx,ny,1,nx-1,1,ny-1,0,tem1(1,1,1))
        IF(mphyopt >= 9) THEN
          CALL smooth9p(qscalar(1,1,k,P_NC),nx,ny,1,nx-1,1,ny-1,0,tem1(1,1,1))
          CALL smooth9p(qscalar(1,1,k,P_NI),nx,ny,1,nx-1,1,ny-1,0,tem1(1,1,1))
        END IF
        IF(mphyopt == 11) THEN
          CALL smooth9p(qscalar(1,1,k,P_ZI),nx,ny,1,nx-1,1,ny-1,0,tem1(1,1,1))
        END IF
      END DO
    ELSE IF(smth_opt == 2) THEN
      IF (P_QC > 0) CALL smooth3d(nx,ny,nz, 1,nx-1,1,ny-1,1,nz-1,0,smfct1,zs,         &
                    qscalar(1,1,1,P_QC),tem1,qscalar(1,1,1,P_QC))
      IF (P_QI > 0) CALL smooth3d(nx,ny,nz, 1,nx-1,1,ny-1,1,nz-1,0,smfct1,zs,         &
                    qscalar(1,1,1,P_QI),tem1,qscalar(1,1,1,P_QI))
      IF(mphyopt >= 9) THEN
        CALL smooth3d(nx,ny,nz, 1,nx-1,1,ny-1,1,nz-1,0,smfct1,zs,         &
                      qscalar(1,1,1,P_NC),tem1,qscalar(1,1,1,P_NC))
        CALL smooth3d(nx,ny,nz, 1,nx-1,1,ny-1,1,nz-1,0,smfct1,zs,         &
                      qscalar(1,1,1,P_NI),tem1,qscalar(1,1,1,P_NI))
      END IF
      IF(mphyopt == 11) THEN
        CALL smooth3d(nx,ny,nz, 1,nx-1,1,ny-1,1,nz-1,0,smfct1,zs,         &
                      qscalar(1,1,1,P_ZI),tem1,qscalar(1,1,1,P_ZI))
      END IF
    END IF
  END IF  ! cldqcopt.eq.1?
!
!-----------------------------------------------------------------------
!
!  Enhance and THEN smooth the precipitate mixing ratio fields.
!
!-----------------------------------------------------------------------
!
  IF (cldqropt == 1 .or. cldqropt==2) THEN
    IF (myproc == 0) WRITE(6,'(a)')' Enhancing qr, qs, and qh-fields'
    DO k=2,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          qr_cld(i,j,k) = MIN(0.001*qr_cld(i,j,k),qrlimit)
          qs_cld(i,j,k) = MIN(0.001*qs_cld(i,j,k),qrlimit)
          qh_cld(i,j,k) = MIN(0.001*qh_cld(i,j,k),qrlimit)
        END DO
      END DO
    END DO

    IF (mp_opt > 0) THEN  ! Before smoothing, the arrays must be MPI valid
      CALL acct_interrupt(mp_acct)
      CALL mpsendrecv2dew(qr_cld,nx,ny,nz,ebc,wbc,0,tem1)
      CALL mpsendrecv2dns(qr_cld,nx,ny,nz,nbc,sbc,0,tem1)
      CALL mpsendrecv2dew(qs_cld,nx,ny,nz,ebc,wbc,0,tem1)
      CALL mpsendrecv2dns(qs_cld,nx,ny,nz,nbc,sbc,0,tem1)
      CALL mpsendrecv2dew(qh_cld,nx,ny,nz,ebc,wbc,0,tem1)
      CALL mpsendrecv2dns(qh_cld,nx,ny,nz,nbc,sbc,0,tem1)
    END IF

    IF (smth_opt == 1) THEN
      DO k=2,nz-1
        CALL smooth9p(qr_cld(1,1,k),nx,ny,1,nx-1,1,ny-1,0,tem1(1,1,1))
        CALL smooth9p(qs_cld(1,1,k),nx,ny,1,nx-1,1,ny-1,0,tem1(1,1,1))
        CALL smooth9p(qh_cld(1,1,k),nx,ny,1,nx-1,1,ny-1,0,tem1(1,1,1))
      END DO
    ELSE IF(smth_opt == 2) THEN
      CALL smooth3d(nx,ny,nz, 1,nx-1,1,ny-1,2,nz-1,0,smfct1,zs,         &
                    qr_cld(1,1,1),tem1,qr_cld(1,1,1))
      CALL smooth3d(nx,ny,nz, 1,nx-1,1,ny-1,2,nz-1,0,smfct1,zs,         &
                    qs_cld(1,1,1),tem1,qs_cld(1,1,1))
      CALL smooth3d(nx,ny,nz, 1,nx-1,1,ny-1,2,nz-1,0,smfct1,zs,         &
                    qh_cld(1,1,1),tem1,qh_cld(1,1,1))
    END IF
    frac = 1.0-frac_qr_2_qc
    DO k=2,nz-1
      DO j=1,ny-1
        DO i=1,nx-1

          IF (P_QC > 0) qscalar(i,j,k,P_QC) = qscalar(i,j,k,P_QC) + frac_qr_2_qc*qr_cld(i,j,k)
          IF (P_QI > 0) qscalar(i,j,k,P_QI) = qscalar(i,j,k,P_QI) +  &
                                frac_qr_2_qc * (qs_cld(i,j,k) + qh_cld(i,j,k))
          IF (P_QR > 0) qscalar(i,j,k,P_QR) = MAX( 0.0, frac*qr_cld(i,j,k) )
          IF (P_QS > 0) qscalar(i,j,k,P_QS) = MAX( 0.0, frac*qs_cld(i,j,k) )
          IF (P_QH > 0) qscalar(i,j,k,P_QH) = MAX( 0.0, frac*qh_cld(i,j,k) )
          IF (P_QG > 0) qscalar(i,j,k,P_QG) = MAX( 0.0, frac*qh_cld(i,j,k) )

          IF(mphyopt >= 9) THEN
            qscalar(i,j,k,P_NC) = ntcloud
            qscalar(i,j,k,P_NI) = 5.*exp(0.304*(273.15-max(233.,t_3d(i,j,k))))
            qscalar(i,j,k,P_NR) = sngl(gammaDP(1.d0+dble(alpharain)))   &
                                 *(n0rain**(3./(4.+alpharain)))         &
                       *(rho(i,j,k)*qscalar(i,j,k,P_QR)/((pi/6.)*1000.* &
                       sngl(gammaDP(4.d0+dble(alpharain)))))**((1.+alpharain)/(4.+alpharain))
            qscalar(i,j,k,P_NS) = sngl(gammaDP(1.d0+dble(alphasnow)))   &
                                 *(n0snow**(3./(4.+alphasnow)))*        &
                       (rho(i,j,k)*qscalar(i,j,k,P_QS)/((pi/6.)*rhosnow* &
                       sngl(gammaDP(4.d0+dble(alphasnow)))))**((1.+alphasnow)/(4.+alphasnow))
            qscalar(i,j,k,P_NH) = sngl(gammaDP(1.d0+dble(alphahail)))   &
                                 *(n0hail**(3./(4.+alphahail)))*        &
                       (rho(i,j,k)*qscalar(i,j,k,P_QH)/((pi/6.)*rhohail* &
                       sngl(gammaDP(4.d0+dble(alphahail)))))**((1.+alphahail)/(4.+alphahail))
            IF(mp_opt > 0) THEN  ! DTD: is the following needed?
              CALL acct_interrupt(mp_acct)
              CALL mpsendrecv2dew(qscalar(1,1,1,P_NR),nx,ny,nz,ebc,wbc,0,tem1)
              CALL mpsendrecv2dns(qscalar(1,1,1,P_NR),nx,ny,nz,nbc,sbc,0,tem1)
              CALL mpsendrecv2dew(qscalar(1,1,1,P_NS),nx,ny,nz,ebc,wbc,0,tem1)
              CALL mpsendrecv2dns(qscalar(1,1,1,P_NS),nx,ny,nz,nbc,sbc,0,tem1)
              CALL mpsendrecv2dew(qscalar(1,1,1,P_NH),nx,ny,nz,ebc,wbc,0,tem1)
              CALL mpsendrecv2dns(qscalar(1,1,1,P_NH),nx,ny,nz,nbc,sbc,0,tem1)
            END IF
          END IF
          IF(mphyopt == 11) THEN
            Gr = ((6.+alpharain)*(5.+alpharain)*(4.+alpharain))/((3.+alpharain)*(2.+alpharain)*(1.+alpharain))
            Gi = ((6.+alphaice)*(5.+alphaice)*(4.+alphaice))/((3.+alphaice)*(2.+alphaice)*(1.+alphaice))
            Gs = ((6.+alphasnow)*(5.+alphasnow)*(4.+alphasnow))/((3.+alphasnow)*(2.+alphasnow)*(1.+alphasnow))
            Gh = ((6.+alphahail)*(5.+alphahail)*(4.+alphahail))/((3.+alphahail)*(2.+alphahail)*(1.+alphahail))
            IF(qscalar(i,j,k,P_NI) /= 0.0) THEN
              qscalar(i,j,k,P_ZI) = Gi/((440.0)**2.)*((rho(i,j,k)*qscalar(i,j,k,P_QI))**2.)/qscalar(i,j,k,P_NI)
            ELSE
              qscalar(i,j,k,P_ZI) = 0.0
            END IF
            IF(qscalar(i,j,k,P_NR) /= 0.0) THEN
              qscalar(i,j,k,P_ZR) = Gr/(((pi/6.)*1000.)**2.)*((rho(i,j,k)*qscalar(i,j,k,P_QR))**2.)/qscalar(i,j,k,P_NR)
            ELSE
              qscalar(i,j,k,P_ZR) = 0.0
            END IF
            IF(qscalar(i,j,k,P_NS) /= 0.0) THEN
              qscalar(i,j,k,P_ZS) = Gs/(((pi/6.)*rhosnow)**2.)*((rho(i,j,k)*qscalar(i,j,k,P_QS))**2.)/qscalar(i,j,k,P_NS)
            ELSE
              qscalar(i,j,k,P_ZS) = 0.0
            END IF
            IF(qscalar(i,j,k,P_NH) /= 0.0) THEN
               qscalar(i,j,k,P_ZH) = Gh/(((pi/6.)*rhohail)**2.)*((rho(i,j,k)*qscalar(i,j,k,P_QH))**2.)/qscalar(i,j,k,P_NH)
            ELSE
               qscalar(i,j,k,P_ZH) = 0.0
            END IF
            IF(mp_opt > 0) THEN ! DTD: is the following needed?
              CALL acct_interrupt(mp_acct)
              CALL mpsendrecv2dew(qscalar(1,1,1,P_ZR),nx,ny,nz,ebc,wbc,0,tem1)
              CALL mpsendrecv2dns(qscalar(1,1,1,P_ZR),nx,ny,nz,nbc,sbc,0,tem1)
              CALL mpsendrecv2dew(qscalar(1,1,1,P_ZS),nx,ny,nz,ebc,wbc,0,tem1)
              CALL mpsendrecv2dns(qscalar(1,1,1,P_ZS),nx,ny,nz,nbc,sbc,0,tem1)
              CALL mpsendrecv2dew(qscalar(1,1,1,P_ZH),nx,ny,nz,ebc,wbc,0,tem1)
              CALL mpsendrecv2dns(qscalar(1,1,1,P_ZH),nx,ny,nz,nbc,sbc,0,tem1)
            END IF
          END IF
        END DO
      END DO
    END DO
  END IF   ! cldqropt.eq.1?

!-----------------------------------------------------------------------
!
!  Adjust the cloud water fields for possible double counting from
!  the combination of the Smith-Feddes cloud scheme and the
!  reflectivity-to-precipitation algorithm.
!
!  For now a simple approach is taken.  If clouds and precipation
!  are present, limit cloud water to 5% of the precipitation value.
!  Then adjust precipitation fields to remove water that already
!  is represented by cloud fields.
!
!  Note possible conflict with frac_qr_2_qc above, so for best
!  results set frac_qr_2_qc to ZERO.
!
!-----------------------------------------------------------------------

  cldpcpfrc=0.05
  DO k=2,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        cldtot = 0.0
        pcptot = 0.0
        IF (P_QC > 0) cldtot = cldtot + qscalar(i,j,k,P_QC)
        IF (P_QI > 0) cldtot = cldtot + qscalar(i,j,k,P_QI)
        IF (P_QR > 0) pcptot = pcptot + qscalar(i,j,k,P_QR)
        IF (P_QS > 0) pcptot = pcptot + qscalar(i,j,k,P_QS)
        IF (P_QG > 0) pcptot = pcptot + qscalar(i,j,k,P_QG)
        IF (P_QH > 0) pcptot = pcptot + qscalar(i,j,k,P_QH)
        IF(cldtot > 0. .AND. pcptot > 0.) THEN
          IF( pcptot < cldtot ) THEN
            cldlim=cldtot-pcptot
          ELSE
            cldlim=cldpcpfrc*pcptot
          END IF
          IF(cldtot > cldlim) THEN
            IF (P_QC > 0) qscalar(i,j,k,P_QC)=qscalar(i,j,k,P_QC)*cldlim/cldtot
            IF (P_QI > 0) qscalar(i,j,k,P_QI)=qscalar(i,j,k,P_QI)*cldlim/cldtot
            IF(mphyopt >= 9) THEN
              qscalar(i,j,k,P_NC) = ntcloud
              qscalar(i,j,k,P_NI) = 5.*exp(0.304*(273.15-max(233.,t_3d(i,j,k))))
            END IF
            IF(mphyopt == 11) THEN
              Gi = ((6.+alphaice)*(5+alphaice)*(4+alphaice))/((3.+alphaice)*(2+alphaice)*(1+alphaice))
              IF(qscalar(i,j,k,P_NI) /= 0.0) THEN
                qscalar(i,j,k,P_ZI) = Gi/((440.0)**2.)*((rho(i,j,k)*qscalar(i,j,k,P_QI))**2.)/qscalar(i,j,k,P_NI)
              ELSE
                qscalar(i,j,k,P_ZI) = 0.0
              END IF
            END IF
          END IF
        END IF
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
!  Adjust the purterbation potential temperature field to account
!  for the latent heating release.
!
!  Corrected by D.H. Wang for ppi term in temperature adjustment
!  Adjustment is to temperature then converted to potential
!  temperature.
!
!-----------------------------------------------------------------------
!
  IF (cldptopt == 3) THEN
    IF (myproc == 0) THEN
      WRITE(6,'(a)')' Adjusting ptprt to account for latent heating.'
      WRITE(6,'(a,f10.4,a,f10.4)')                                      &
         ' frac of qc:',frac_qc_2_lh,' adj_lim:',max_lh_2_pt
    END IF
    p0inv=1./p0
    max_pt_adj = 0.0
    DO k=2,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          ppi = ((pbar(i,j,k)+pprt(i,j,k))*p0inv) ** rddcp
!          arg = lathv*frac_qc_2_lh*(qc(i,j,k)+qi(i,j,k))/(cp*ppi)
          arg = 0.0
          IF (P_QC > 0) arg = arg + qscalar(i,j,k,P_QC)
          IF (P_QI > 0) arg = arg + qscalar(i,j,k,P_QI)
          arg = lathv*frac_qc_2_lh*arg/(cp*ppi)
          max_pt_adj = MAX(max_pt_adj,arg)
          ptprt(i,j,k) = ptprt(i,j,k) + MIN(arg,max_lh_2_pt)
          max_pt_adj = MAX(max_pt_adj,arg)
        END DO
      END DO
    END DO
!
!   In this call calls below, "min_pt_adj" is just a dummy argument, as the
!   subroutine called returns a max and a min.
!
    min_pt_adj = 0
    CALL mpmax0(max_pt_adj,min_pt_adj)
    IF (myproc == 0) PRINT*,'max_adj=',max_pt_adj
  ELSE IF (cldptopt == 4) THEN
    IF (myproc == 0) THEN
      WRITE(6,'(a)')' Adjusting ptprt to account for latent heating in w.'
      PRINT*,'frac of qc:',frac_qc_2_lh,' adj_lim:',max_lh_2_pt
    END IF
    max_pt_adj = 0.0
    DO k=2,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          !
          ! Ningzhu's change on 3/23/2011
          !
          !IF(w(i,j,k) > 0. ) THEN
          !  wratio=1.0
          !  ptcld=ctmp_3d(i,j,k)*(p0/p_3d(i,j,k))**rddcp
          !  ptdiff=ptcld-(ptbar(i,j,k)+ptprt(i,j,k))
          !  IF(ptdiff > 0.) THEN
          !    arg = frac_qc_2_lh*wratio*ptdiff
          !    ptprt(i,j,k) = ptprt(i,j,k) + MIN(arg,max_lh_2_pt)
          !    max_pt_adj = MAX(max_pt_adj,arg)
          !  END IF
          !END IF
          IF(w(i,j,k) > 0. ) THEN
            ppi = ((pbar(i,j,k)+pprt(i,j,k))*p0inv) ** rddcp
            !arg = lathv*frac_qc_2_lh*(qc(i,j,k)+qi(i,j,k))/(cp*ppi)
            arg = 0.0
            IF (P_QC > 0) arg = arg + qscalar(i,j,k,P_QC)
            IF (P_QI > 0) arg = arg + qscalar(i,j,k,P_QI)
            arg = lathv*frac_qc_2_lh*arg/(cp*ppi)
            ptprt(i,j,k) = ptprt(i,j,k) + MIN(arg,max_lh_2_pt)
            max_pt_adj = MAX(max_pt_adj,arg)
          END IF
        END DO
      END DO
    END DO
    min_pt_adj = 0
    CALL mpmax0(max_pt_adj,min_pt_adj)
    IF (myproc == 0) PRINT*,'max_adj=',max_pt_adj
  ELSE IF (cldptopt == 5) THEN
    IF (myproc == 0) THEN
      WRITE(6,'(a)')' Adjusting ptprt to moist-adiab cloud temp for w>-0.2'
      PRINT*,'frac of qc:',frac_qc_2_lh,' adj_lim:',max_lh_2_pt
    END IF
    max_pt_adj = 0.0
    DO k=2,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          wratio=min(max(0.,(5.0*(w(i,j,k)+0.2))),1.0)
          ptcld=ctmp_3d(i,j,k)*(p0/p_3d(i,j,k))**rddcp
          ptdiff=ptcld-(ptbar(i,j,k)+ptprt(i,j,k))
          IF(ptdiff > 0.) THEN
            arg = frac_qc_2_lh*wratio*ptdiff
            ptprt(i,j,k) = ptprt(i,j,k) + MIN(arg,max_lh_2_pt)
            max_pt_adj = MAX(max_pt_adj,arg)
          END IF
        END DO
      END DO
    END DO
    min_pt_adj = 0
    CALL mpmax0(max_pt_adj,min_pt_adj)
    IF (myproc == 0) PRINT*,'max_adj=',max_pt_adj

  ELSE IF (cldptopt == 6) THEN
    IF (myproc == 0) THEN
      WRITE(6,'(a)')' Adjusting ptprt to moist-adiab cloud temp for w>0.0'
      PRINT*,'frac of qc:',frac_qc_2_lh,' adj_lim:',max_lh_2_pt
    END IF
    max_pt_adj = 0.0
    DO k=2,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          IF(w(i,j,k) > 0.) THEN
            ptcld=ctmp_3d(i,j,k)*(p0/p_3d(i,j,k))**rddcp
            ptdiff=ptcld-(ptbar(i,j,k)+ptprt(i,j,k))
            IF(ptdiff > 0.) THEN
              arg = frac_qc_2_lh*ptdiff
              ptprt(i,j,k) = ptprt(i,j,k) + MIN(arg,max_lh_2_pt)
              max_pt_adj = MAX(max_pt_adj,arg)
            END IF
          END IF
        END DO
      END DO
    END DO
    min_pt_adj = 0
    CALL mpmax0(max_pt_adj,min_pt_adj)
    IF (myproc == 0) PRINT*,'max_adj=',max_pt_adj

  END IF   ! cldptopt=3?
!
!-----------------------------------------------------------------------
!
!  Enhance rh*-field in the cloudy area with the cloud cover
!  greater than thresh_cvr and then smooth the enhanced field.
!
!-----------------------------------------------------------------------
!
  IF (cldqvopt == 1) THEN
    IF (myproc == 0) WRITE(6,'(a,f5.2,a,f5.2,a,f5.1,a,f5.1,a)')         &
            ' Enhancing RH-field for cldcvr between'                    &
             ,cvr2rh_thr1,' and ', cvr2rh_thr2                          &
             ,'with linear ramp from ',rh_thr1*100.0                    &
             ,'% to ',rh_thr2*100.0,'%'

    DO k=2,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          IF(clouds_3d(i,j,k) > cvr2rh_thr1) THEN
            rh = (clouds_3d(i,j,k)-cvr2rh_thr1)*(rh_thr2-rh_thr1)       &
                    /(cvr2rh_thr2-cvr2rh_thr1) + rh_thr1
            rh = MIN(rh,rh_thr2)
            rh_3d(i,j,k) = MAX((100.*rh),rh_3d(i,j,k))
          END IF
        END DO
      END DO
    END DO

    IF (mp_opt > 0) THEN  ! Before smoothing, the arrays must be MPI valid
      CALL acct_interrupt(mp_acct)
      CALL mpsendrecv2dew(rh_3d,nx,ny,nz,ebc,wbc,0,tem1)
      CALL mpsendrecv2dns(rh_3d,nx,ny,nz,nbc,sbc,0,tem1)
    END IF

    IF (smth_opt == 1) THEN
      DO k=2,nz-1
        CALL smooth9p(rh_3d(1,1,k),nx,ny,1,nx-1,1,ny-1,0,tem1(1,1,1))
      END DO
    ELSE IF(smth_opt == 2) THEN
      CALL smooth3d(nx,ny,nz, 1,nx-1,1,ny-1,2,nz-1,0,smfct1,zs,         &
                    rh_3d,tem1,rh_3d)

    END IF
!
!-----------------------------------------------------------------------
!
!  Convert the rh* field back into qv-field
!  Where precip was removed, set qv according to rhsubsat
!
!-----------------------------------------------------------------------
!
    print *, ' Here QV adjustment loop 1'
    DO k=2,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          qvsat=f_qvsat( p_3d(i,j,k), t_3d(i,j,k))
          rh2qv=0.01*rh_3d(i,j,k)*qvsat
          qpcp=qr_cld(i,j,k)+qs_cld(i,j,k)+qh_cld(i,j,k)
          IF ( qpcp == 0.0 ) THEN  ! no precip here
            qw_new = 0.0
            DO nq=1,nscalarq
              qw_new = qw_new + qscalar(i,j,k,nq)
            END DO
            dqw = qw_new - qw_old(i,j,k)
            IF( dqw < 0.0 ) THEN   ! total water was reduced
              qv(i,j,k)=MIN(qv(i,j,k),rh2qv)
              qv(i,j,k)=MIN(qv(i,j,k),(rhsubsat*qvsat))
            ELSE
              qv(i,j,k)=MAX(qv(i,j,k),rh2qv)
            END IF
          ELSE
            qv(i,j,k)=MAX(qv(i,j,k),rh2qv)
          END IF
        END DO
      END DO
    END DO

  END IF          ! cldqvopt=1
!
!-----------------------------------------------------------------------
!
!  Adjust the purterbation potential temperature field to preserve
!  the previous buoyancy field.
!
!-----------------------------------------------------------------------
!
  IF (cldptopt > 0.AND.cldptopt < 3) THEN
    IF (myproc == 0) THEN
      WRITE(6,'(a)')' Adjusting ptprt-field to preserve buoyancy'
      PRINT*,'frac of qw:',frac_qw_2_pt
    END IF
    DO k=2,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          qw_new = 0.0
          DO nq=1,nscalarq
            qw_new = qw_new + qscalar(i,j,k,nq)
          END DO
          dqw = qw_new - qw_old(i,j,k)
          IF (cldptopt == 1) THEN
            arg =  dqw/(1.0+qvbar(i,j,k))
          ELSE
            dqv_prt = qv(i,j,k) - qvbar(i,j,k) - qvprt_old(i,j,k)
            arg = (dqv_prt + dqw)/(1.0+qvbar(i,j,k))                    &
                    - dqv_prt/(0.622+qvbar(i,j,k))
          END IF
          ptprt(i,j,k) = ptprt(i,j,k) + ptbar(i,j,k)*arg*frac_qw_2_pt
        END DO
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Re-adjust the qv field to account for the temp. adjustment
!  Where the water fields were reduced, limit rh to rhsubsat
!
!-----------------------------------------------------------------------
!
    IF (cldqvopt == 1) THEN

      print *, ' Here QV adjustment loop 2'

      DO k=2,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            tgrid =(ptprt(i,j,k)+ptbar(i,j,k))*((p_3d(i,j,k)/p0)**rddcp)
            qvsat=f_qvsat( p_3d(i,j,k), tgrid )
            rh2qv=0.01*rh_3d(i,j,k)*qvsat
            qpcp=qr_cld(i,j,k)+qs_cld(i,j,k)+qh_cld(i,j,k)
            IF ( qpcp == 0.0 ) THEN  ! no precip here
              qw_new = 0.0
              DO nq=1,nscalarq
                qw_new = qw_new + qscalar(i,j,k,nq)
              END DO
              dqw = qw_new - qw_old(i,j,k)
              IF( dqw < 0.0 ) THEN   ! total water was reduced
                qv(i,j,k)=MIN(qv(i,j,k),rh2qv)
                qv(i,j,k)=MIN(qv(i,j,k),(rhsubsat*qvsat))
              ELSE
                qv(i,j,k)=MAX(qv(i,j,k),rh2qv)
              END IF
            ELSE
              qv(i,j,k)=MAX(qv(i,j,k),rh2qv)
            END IF
          END DO
        END DO
      END DO

    END IF          ! cldqvopt=1
  END IF     ! cldptopt.eq.1?
!
!
!-----------------------------------------------------------------------
!
!  Enhance and THEN smooth w-field using analyzed in-cloud w-field.
!
!-----------------------------------------------------------------------
!
  IF (cldwopt == 1) THEN
    IF (myproc == 0) WRITE(6,'(a,f8.4)')  &
      ' Enhancing w-field for areas w/ cldcvr >', thresh_cvr

    DO k=2,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          IF(w_cld(i,j,k) > w(i,j,k)) THEN
            w(i,j,k) = w_cld(i,j,k)
          END IF
        END DO
      END DO
    END DO

    IF (mp_opt > 0) THEN  ! Before smoothing, the arrays must be MPI valid
      CALL acct_interrupt(mp_acct)
      CALL mpsendrecv2dew(w,nx,ny,nz,ebc,wbc,3,tem1)
      CALL mpsendrecv2dns(w,nx,ny,nz,nbc,sbc,3,tem1)
    END IF

    IF (smth_opt == 1) THEN
      DO k=2,nz-1
        CALL smooth9p(w(1,1,k),nx,ny,1,nx-1,1,ny-1,0,tem1(1,1,1))
      END DO
    ELSE IF(smth_opt == 2) THEN
      CALL smooth3d(nx,ny,nz, 1,nx-1,1,ny-1,2,nz-1,3,smfct1,zs,         &
                    w(1,1,1),tem1,w(1,1,1))
    END IF

  END IF   ! cldwopt = 1
!
!-----------------------------------------------------------------------
!
  999   CONTINUE

  RETURN
END SUBROUTINE cmpclddrv
!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE INICLDGRD                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE inicldgrd (nx,ny,nz,zs,                                      &
           bgqcopt,default_clear_cover,                                 &
           topo,t_3d,p_3d,rh_3d,cf_modelfg,t_sfc_k,psfc_pa,             &
           cldcv,wtcldcv,z_ref_lcl,z_lcl,rh_modelfg,                    &
           r_missing,tem1,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Create background temperature and cloud fractional
!  cover fields on the cloud analysis grid (which has the same
!  horizontal grid as the ADAS, but the vertical levels are
!  different from ADAS).
!
!  This should probably be free of very small scale horizontal
!  structures for best results when combining with the satellite
!  data (cloud analysis uses surface observations, radar echo,
!  and satellite imagery data).
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (Jian Zhang)
!  03/1996
!
!  MODIFICATION HISTORY
!
!  03/14/97  J. Zhang
!            Cleaning up the code and implemented for the official
!            arps4.2.4 version
!  09/10/97  J. Zhang
!            Using a quadratic relationship for deriving cloud
!            fractional cover field from gridded relative humidity
!            analysis. Added calculation for lifting condensation
!            levels.
!  09/15/97  J. Zhang
!            Fixed a bug when calling function rh_to_cldcv
!  05/06/98  J. Zhang
!            Abandoned the cloud grid, using the ARPS grid instead.
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  nx, ny, nz              ! ADAS grid size.
!
!c Analysis variables
!
!  zs (nx,ny,nz)           ! The physical height coordinate
                           ! defined at w-point of staggered grid.
!  t_3d (nx,ny,nz)         ! Temperature field
!  p_3d (nx,ny,nz)         ! Pressure field
!  rh_3d(nx,ny,nz)         ! Relative humidity field
!
!  OUTPUT:
!
!   REAL*4 cf_modelfg(nx,ny,nz)  ! model first guess for cld cv.
!   REAL*4 rh_modelfg(nx,ny,nz)  ! model first guess for RH
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INCLUDE 'mp.inc'
  INCLUDE 'bndry.inc'
!
!-----------------------------------------------------------------------
!
!  Variables declaration
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
  INTEGER :: nx,ny,nz              ! ADAS grid size
!
  INTEGER :: bgqcopt
  REAL :: z_ref_lcl                ! ref. level for computing LCL
  REAL :: default_clear_cover      ! default value for clear sky
  REAL :: r_missing                ! bad or missing data flag
!
  REAL :: zs (nx,ny,nz)            ! Hgts of each ADAS gird pt.
  REAL :: topo (nx,ny)             ! Terrain height
  REAL :: t_3d (nx,ny,nz)          ! Temperature field (K)
  REAL :: p_3d (nx,ny,nz)          ! Pressure field (pa)
  REAL :: rh_3d (nx,ny,nz)         ! relative humidity (0-100) field
!
!  OUTPUT:
!
  REAL :: cf_modelfg(nx,ny,nz)     ! Output, 1st guess of cloud
                                 ! cover on cloud height grid
  REAL :: rh_modelfg(nx,ny,nz)     ! Output, 1st guess of relative
                                 ! humidity on cloud height grid
  REAL :: t_sfc_k(nx,ny)           ! Air temp. at sfc
  REAL :: psfc_pa(nx,ny)           ! pressure (Pascal) at sfc
  REAL :: z_lcl(nx,ny)             ! lifting condensatn lvl (MSL)
!
  REAL :: cldcv (nx,ny,nz)         ! 3D gridded frac cld cv analysis.
  REAL :: wtcldcv (nx,ny,nz)       ! wgt assigned to cld cvr analysis
!
  INTEGER :: istatus               ! flag indicate process status
!
!   WORK variable:
!
  REAL :: tem1 (nx,ny,nz)          ! MPI temporary variable
!
!  FUNCTIONS:
  REAL :: rh_to_cldcv,dwpt
!
!  CONSTANTS:
  REAL :: gamma_d   ! dry adiabatic lapse rate (K/m)
  PARAMETER (gamma_d = 9.8/1004.0)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  REAL :: arg,frac_z
  REAL :: t_ref_k,t_ref_c,rh_ref,td_ref_c,z_ref
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (myproc == 0) WRITE(6,*) 'Initializing the cloud analysis grid'
  istatus = 0
!
!-----------------------------------------------------------------------
!
!  Initialize surface temperature fields with missing flag
!  Initialize background cloud cover and the weight arrays.
!  Initialize model first guess fields with default values
!
!-----------------------------------------------------------------------
!
  DO j=1,ny
    DO i=1,nx
      t_sfc_k(i,j) = r_missing
    END DO
  END DO

  DO k=1,nz
    DO j=1,ny
      DO i=1,nx
        cldcv(i,j,k) = r_missing
        wtcldcv(i,j,k) = r_missing
        cf_modelfg(i,j,k) = default_clear_cover
        rh_modelfg(i,j,k) = 0.0
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Calculate the surface temperature using standard atmosphere
!  profile.
!
!-----------------------------------------------------------------------
!
  DO j = 1, ny-1
    DO i = 1, nx-1
      psfc_pa(i,j) = 2.*p_3d(i,j,2)-p_3d(i,j,3)
      t_sfc_k(i,j) = 2.*t_3d(i,j,2)-t_3d(i,j,3)
    END DO
  END DO
  DO i = 1, nx-1
    psfc_pa(i,ny) = psfc_pa(i,ny-1)
    t_sfc_k(i,ny) = t_sfc_k(i,ny-1)
  END DO
  DO j = 1, ny
    psfc_pa(nx,j) = psfc_pa(nx-1,j)
    t_sfc_k(nx,j) = t_sfc_k(nx-1,j)
  END DO
!
!-----------------------------------------------------------------------
!
!  Find the lifting condensation level
!
!-----------------------------------------------------------------------
!
  DO j = 1,ny-1
    DO i = 1,nx-1
      z_ref = z_ref_lcl + topo(i,j)
      z_ref = min(max(z_ref,(zs(i,j,2)+1.0)),(zs(i,j,nz-1)-1.0))
!
!-----------------------------------------------------------------------
!
!  Find the model temperature and dewpoint at the reference level.
!
!-----------------------------------------------------------------------
!
      DO k = 3,nz-1
        IF (z_ref >= zs(i,j,k-1)) THEN
          frac_z = (z_ref-zs(i,j,k-1))/(zs(i,j,k)-zs(i,j,k-1))
          t_ref_k = t_3d(i,j,k-1)                                       &
                     + frac_z*(t_3d(i,j,k)-t_3d(i,j,k-1))
          t_ref_c = t_ref_k - 273.15
!
          rh_ref = rh_3d(i,j,k-1)                                       &
                 + frac_z*(rh_3d(i,j,k)-rh_3d(i,j,k-1))
          td_ref_c = dwpt(t_ref_c,rh_ref)
        END IF
      END DO  ! k = 2,nz-1
!
      z_lcl(i,j) = z_ref + (t_ref_c - td_ref_c)/gamma_d
      z_lcl(i,j) = min(zs(i,j,nz-1),max(z_lcl(i,j),zs(i,j,2)))

    END DO  ! I
  END DO  ! J
  DO j = 1, ny-1
    z_lcl(nx,j) = z_lcl(nx-1,j)
  END DO
  DO i = 1, nx
    z_lcl(i,ny) = z_lcl(i,ny-1)
  END DO
!
!-----------------------------------------------------------------------
!
!  Remap model temperature and rh fields to cloud height grid
!  and convert model rh to cloud cover
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) WRITE(6,*)
!
  DO k = 2,nz-1
    DO j = 1,ny-1
      DO i = 1,nx-1
        rh_modelfg(i,j,k) = 0.01*rh_3d(i,j,k)
        IF (zs(i,j,k) >= z_lcl(i,j) .AND. bgqcopt > 0) THEN
          arg = zs(i,j,k) - topo(i,j)
          cf_modelfg(i,j,k) = rh_to_cldcv(rh_modelfg(i,j,k),arg)
        END IF
      END DO  ! I
    END DO  ! J
  END DO  ! K

  DO i = 1, nx-1
    DO j = 1, ny-1
      rh_modelfg(i,j,1) = rh_modelfg(i,j,2)
      rh_modelfg(i,j,nz) = rh_modelfg(i,j,nz-1)
    END DO
  END DO

  DO k = 1,nz
    DO j = 1,ny-1
      rh_modelfg(nx,j,k) = rh_modelfg(nx-1,j,k)
      cf_modelfg(nx,j,k) = cf_modelfg(nx-1,j,k)
    END DO  ! J
    DO i = 1,nx
      rh_modelfg(i,ny,k) = rh_modelfg(i,ny-1,k)
      cf_modelfg(i,ny,k) = cf_modelfg(i,ny-1,k)
    END DO  ! I

  END DO  ! K

  IF (mp_opt > 0) THEN
    CALL mpsendrecv2dew(cf_modelfg,nx,ny,nz,ebc,wbc,1,tem1)
    CALL mpsendrecv2dns(cf_modelfg,nx,ny,nz,nbc,sbc,2,tem1)
  END IF

  istatus=1

  RETURN
END SUBROUTINE inicldgrd

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE READRAD_JZ                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE readrad_jz(nx,ny,nz,isrcrad,stnrad                           &
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
!  08/06/97  J. Zhang
!            Change adascld24.inc to adascld25.inc.
!  09/11/97  J. Zhang
!            Change adascld25.inc to adascld26.inc.
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
!    fradname   file name for radar datasets
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
!    tem1     Temporary work array.
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
!-----------------------------------------------------------------------

  INCLUDE 'adas.inc'     ! ADAS parameters
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!-----------------------------------------------------------------------

  INTEGER :: nx,ny,nz    ! the ARPS grid size
!
!  INCLUDE:   (from adas.inc)
!  integer mx_rad
!  INPUT:   (from namelist in the .input file)
!  integer nradfil
!  character*132 fradname(mx_rad)
!
!  LOCAL:
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
  INTEGER :: iopt_wrtrad
  PARAMETER (iopt_wrtrad=2)
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
!  Temporary working array
!
!-----------------------------------------------------------------------
!
  REAL :: tem1(nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=4)   :: stn
  CHARACTER (LEN=6)   :: runname
  CHARACTER (LEN=256) :: fname
  INTEGER :: ireftim,itime,vcpnum,idummy
  INTEGER :: hdmpfmt,strhopt,mapprin
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
  PRINT*,'Reading radar data from',nradfil,' sites.'
  IF (nradfil > mx_rad) THEN
    WRITE(6,'(a,i3,a,i3/a)')                                            &
        ' ERROR: nradfil ',nradfil,' exceeds mx_rad dimension',         &
        mx_rad,' please increase MX_RAD in the .inc file'
    PRINT*,' ABORTING from READRAD......'
    CALL arpsstop("MX_RAD array too small",1)
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
  DO krad = 1, mx_rad
    DO k=1,nz
      DO j=1,ny
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
!
    READ(nchanl) stn
    stnrad(krad)=stn
    PRINT*,'Reading ',stnrad(krad),' radar data from ',                 &
              radfname(krad)
    READ(nchanl) ireftim,itime,vcpnum,idummy,idummy,                    &
               idummy,idummy,idummy,idummy,idummy
!
    CALL abss2ctim(itime, iyr, imon, idy, ihr, imin, isec )
    iyr=MOD(iyr,100)
    WRITE(6,815) imon,idy,iyr,ihr,imin
    815   FORMAT(i2.2,'/',i2.2,'/',i2.2,1X,i2.2,':',i2.2,' UTC')
!
    READ(nchanl) runname
    READ(nchanl) hdmpfmt,strhopt,mapprin,idummy,idummy,                 &
               idummy,idummy,idummy,idummy,idummy

    READ(nchanl) dxin,dyin,dzin,dzminin,ctrlatin,                       &
             ctrlonin,tlat1in,tlat2in,tlonin,scalin,                    &
             latrad(krad),lonrad(krad),elvrad(krad),                    &
             rdummy,rdummy
!
    IF (iopt_wrtrad == 2) THEN
!
!-----------------------------------------------------------------------
!
!  Read the data file created from subroutine "WRTRAD"
!
!-----------------------------------------------------------------------
!
      READ(nchanl) tem1    ! Reflectivity
      DO i=1,nx
        DO j=1,ny
          DO k=1,nz
            gridref(i,j,k,krad) = tem1(i,j,k)
          END DO
        END DO
      END DO
!
      READ(nchanl) tem1    ! Radial Velocity
      DO i=1,nx
        DO j=1,ny
          DO k=1,nz
            gridvel(i,j,k,krad) = tem1(i,j,k)
          END DO
        END DO
      END DO

      READ(nchanl) tem1    ! Nyquist Velocity
      DO i=1,nx
        DO j=1,ny
          DO k=1,nz
            gridnyq(i,j,k,krad) = tem1(i,j,k)
          END DO
        END DO
      END DO

      READ(nchanl) tem1    ! Time (scnds) from the reference time
      DO i=1,nx
        DO j=1,ny
          DO k=1,nz
            gridtim(i,j,k,krad) = tem1(i,j,k)
          END DO
        END DO
      END DO
!
    ELSE IF (iopt_wrtrad == 1) THEN
!
!-----------------------------------------------------------------------
!
!  Read the data file created from subroutine "WTRADCOL"
!
!-----------------------------------------------------------------------
!
      DO k=1,nz
        kntref(k) = 0
        kntvel(k) = 0
      END DO

      READ(nchanl) nradvr,iradvr
!
      DO ipt=1,(nx*ny)

        READ(nchanl,END=51) i,j,xrd,yrd,                                &
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
              IF (gridref(i,j,k,krad) > -200.                           &
                  .AND. gridref(i,j,k,krad) < 200.)                     &
                  kntref(k)=kntref(k)+1
              IF (gridvel(i,j,k,krad) > -200.                           &
                  .AND. gridvel(i,j,k,krad) < 200.)                     &
                  kntvel(k)=kntvel(k)+1
            END IF  ! 1 < k < nz
          END DO  ! kk = 1, klev
        END IF  ! 1 < i < nx  & 1 < j < ny

      END DO  ! ipt = 1, nx*ny

      51        CONTINUE
      ipt=ipt-1
      WRITE(6,'(a,i6,a)') ' End of file reached after reading',         &
                         ipt,' columns'
      GO TO 55
      52        CONTINUE
      WRITE(6,'(a,i6,a)') ' End of file reached while reading',         &
                         ipt,' column'
      55        CONTINUE
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
!
      CLOSE(nchanl)
      CALL retunit( nchanl )

    END IF  ! iopt_wrtrad
    GO TO 400

    399     CONTINUE
    PRINT*,'Error reading the radar file:',fname

    400     CONTINUE
  END DO  ! KRAD = 1, nradfil

  RETURN
END SUBROUTINE readrad_jz
