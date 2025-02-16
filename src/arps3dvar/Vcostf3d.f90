!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE  COSTF                    ######
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
!  Define the total costfunction for analysis
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!
!  Jidong Gao, CAPS, July, 2000
!
!-----------------------------------------------------------------------

SUBROUTINE costf(ount,numctr,ctrv, cfun_single,                         &
           gdu_err,gdv_err,gdp_err,gdt_err,gdq_err,gdw_err,             &
           u_ctr,v_ctr,p_ctr,t_ctr,q_ctr,w_ctr, psi, phi,               &
           gdscal, nx,ny,nz,                                            &
           nvar,nvarradin,nvarrad,nzua,nzrdr,nzret,                     &
           mapfct,j1,j2,j3,aj3x,aj3y,aj3z,j3inv,rhostr,                 &
           rhostru, rhostrv, rhostrw, div3,                             &
           mxsng,mxua,mxrad,mxcolrad,mxret,mxcolret,                    &
           nsrcsng,nsrcua,nsrcrad,nsrcret,ncat,                         &
           mxpass,npass,iwstat,xs,ys,zs,x,y,z,zp,hterain,               &
           icatg,xcor,nam_var,                                          &
           ownsng,usesng,xsng,ysng,hgtsng,thesng,                       &
           obsng,odifsng,qobsng2,qualsng,isrcsng,icatsng,nobsng,        &
           ownua,useua,xua,yua,hgtua,theua,                             &
           obsua,odifua,qobsua2,qualua,isrcua,nlevsua,nobsua,           &
           ownrad,userad,elvrad,xradc,yradc,                            &
           distrad,uazmrad,vazmrad,hgtradc,theradc,dsdr,dhdr,           &
           obsrad,odifrad,qobsrad2,qualrad,                             &
           irad,isrcrad,nlevrad,ncolrad,                                &
           xretc,yretc,hgtretc,theretc,                                 &
           obsret,odifret,qobsret,qualret,                              &
           iret,isrcret,nlevret,ncolret,                                &
           srcsng,srcua,srcrad,srcret,                                  &
           ianxtyp,iusesng,iuseua,iuserad,iuseret,                      &
           xyrange,kpvrsq,wlim,zrange,zwlim,                            &
           thrng,rngsqi,knt,wgtsum,zsum,                                &
           corsng,corua,corrad,corret,                                  &
           xsng_p,ysng_p,ihgtsng,xua_p,yua_p,ihgtua,                    &
           xradc_p,yradc_p,ihgtradc,zsng_1,zsng_2,                      &
           zua_1,zua_2,zradc_1,zradc_2,                                 &
           oanxsng,oanxua,oanxrad,oanxret,                              &
           sngsw, uasw, radsw, retsw,                                   &
           ipass_filt,hradius,radius_z,                                &
           div_opt,cntl_var, smth_opt,                                  &
           wgt_div_h,wgt_div_v,wgt_smth,                                &
           thermo_opt, wgt_thermo, sinlat,ffu,ffv,                      &
           anx,tem1,tem2,tem3,tem4,smcu,smcv,smcw,                      &
           u,v,w,wcont,istatus)

!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  USE arps_precision

  IMPLICIT NONE
!
!  INCLUDE 'varpara.inc'
!
  INCLUDE 'phycst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'bndry.inc'     ! to include ebc, nbc, sbc, wbc etc.
  INCLUDE 'globcst.inc'   ! to include mp_acct etc.
  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  Input Sizing Arguments
!
!-----------------------------------------------------------------------
!
  INTEGER, INTENT(IN) :: ount    ! output unit
  INTEGER, INTENT(IN) :: nx,ny,nz,sngsw, uasw, radsw, retsw
  INTEGER, INTENT(IN) :: nvar,nvarradin,nvarrad
  INTEGER, INTENT(IN) :: nzua,nzrdr,nzret
  INTEGER, INTENT(IN) :: mxsng,mxua,mxrad,mxcolrad,mxret,mxcolret
  INTEGER, INTENT(IN) :: nsrcsng,nsrcua,nsrcrad,nsrcret,ncat
  INTEGER, INTENT(IN) :: mxpass,npass
  INTEGER, INTENT(IN) :: ipass_filt,div_opt,cntl_var
  INTEGER, INTENT(IN) :: smth_opt, thermo_opt
  REAL,    INTENT(IN) :: radius_z(nx,ny,nz)
  REAL    :: hradius
  REAL    :: wgt_div_h, wgt_div_v, wgt_smth, wgt_thermo

  LOGICAL, INTENT(IN) :: ownsng(mxsng), ownua(mxua), ownrad(mxcolrad)
  LOGICAL, INTENT(IN) :: usesng(mxsng), useua(mxua), userad(mxcolrad)
!
!-----------------------------------------------------------------------
!
!  input grid arguments
!
!-----------------------------------------------------------------------
!
!
  REAL(P) :: x     (nx)        ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL(P) :: y     (ny)        ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL(P) :: z     (nz)        ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL(P) :: zp    (nx,ny,nz)  ! The physical height coordinate defined at
                               ! w-point of the staggered grid.
  REAL(P) :: hterain(nx,ny)    ! The height of the terrain.

  REAL(P) :: mapfct(nx,ny,8)   ! Map factors at scalar, u and v points

  REAL(P) :: j1    (nx,ny,nz)  ! Coordinate transformation Jacobian
                               ! defined as - d( zp )/d( x ).
  REAL(P) :: j2    (nx,ny,nz)  ! Coordinate transformation Jacobian
                               ! defined as - d( zp )/d( y ).
  REAL(P) :: j3    (nx,ny,nz)  ! Coordinate transformation Jacobian
                               ! defined as d( zp )/d( z ).
  REAL(P) :: aj3x  (nx,ny,nz)  ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE X-DIR.
  REAL(P) :: aj3y  (nx,ny,nz)  ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Y-DIR.
  REAL(P) :: aj3z  (nx,ny,nz)  ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Z-DIR.
  REAL(P) :: j3inv (nx,ny,nz)     ! Inverse of j3
  REAL(P) :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3

  REAL(P) :: rhostru(nx,ny,nz)    ! Averaged rhostr at u points (kg/m**3).
  REAL(P) :: rhostrv(nx,ny,nz)    ! Averaged rhostr at v points (kg/m**3).
  REAL(P) :: rhostrw(nx,ny,nz)    ! Averaged rhostr at w points (kg/m**3).
!
  REAL(P) :: xs(nx)
  REAL(P) :: ys(ny)
  REAL(P) :: zs(nx,ny,nz)
  INTEGER :: icatg(nx,ny)
  REAL(P) :: xcor(ncat,ncat)
!
!-----------------------------------------------------------------------
!
!  Input Observation Arguments
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=6) :: nam_var(nvar)
  REAL(P) :: xsng(mxsng)
  REAL(P) :: ysng(mxsng)
  REAL(P) :: hgtsng(mxsng)
  REAL(P) :: thesng(mxsng)
  REAL(P) :: obsng(nvar,mxsng)
  REAL(P) :: odifsng(nvar,mxsng)
  REAL(P) :: qobsng2(nvar,mxsng)
  INTEGER :: qualsng(nvar,mxsng)
  INTEGER :: isrcsng(mxsng)
  INTEGER :: icatsng(mxsng)
  INTEGER :: nobsng

!  REAL(P) :: xsng_p(mxsng),ysng_p(mxsng)
  REAL(DP):: xsng_p(mxsng),ysng_p(mxsng)
  REAL(P) :: zsng_1(mxsng),zsng_2(mxsng)
  INTEGER :: ihgtsng(mxsng)
!
  REAL(P) :: xua(mxua)
  REAL(P) :: yua(mxua)
  REAL(P) :: hgtua(nzua,mxua)
  REAL(P) :: theua(nzua,mxua)
  REAL(P) :: obsua(nvar,nzua,mxua)
  REAL(P) :: odifua(nvar,nzua,mxua)
  REAL(P) :: qobsua2(nvar,nzua,mxua)
  INTEGER :: qualua(nvar,nzua,mxua)
  INTEGER :: nlevsua(mxua)
  INTEGER :: isrcua(mxua)
  INTEGER :: nobsua
!
!  REAL(P) :: xua_p(mxua),yua_p(mxua)
  double precision :: xua_p(mxua),yua_p(mxua)
  REAL(P) :: zua_1(nzua,mxua),zua_2(nzua,mxua)
  INTEGER :: ihgtua(nzua,mxua)
!
  REAL(P) :: elvrad(mxrad)
  REAL(P) :: xradc(mxcolrad)
  REAL(P) :: yradc(mxcolrad)
  REAL(P) :: distrad(mxcolrad)
  REAL(P) :: uazmrad(mxcolrad)
  REAL(P) :: vazmrad(mxcolrad)
  REAL(P) :: hgtradc(nzrdr,mxcolrad)
  REAL(P) :: theradc(nzrdr,mxcolrad)
  REAL(P) :: dsdr(nzrdr,mxcolrad)
  REAL(P) :: dhdr(nzrdr,mxcolrad)
  REAL(P) :: obsrad(nvarradin,nzrdr,mxcolrad)
  REAL(P) :: odifrad(nvarrad,nzrdr,mxcolrad)
  REAL(P) :: qobsrad2(nvarrad,nzrdr,mxcolrad)
  INTEGER :: qualrad(nvarrad,nzrdr,mxcolrad)
  INTEGER :: nlevrad(mxcolrad)
  INTEGER :: irad(mxcolrad)
  INTEGER :: isrcrad(0:mxrad)
  INTEGER :: ncolrad
!
!  REAL(P) :: xradc_p(mxcolrad),yradc_p(mxcolrad)
  double precision :: xradc_p(mxcolrad),yradc_p(mxcolrad)
  REAL(P) :: zradc_1(nzrdr,mxcolrad),zradc_2(nzrdr,mxcolrad)
  INTEGER :: ihgtradc(nzrdr,mxcolrad)
!
  REAL(P) :: xretc(mxcolret)
  REAL(P) :: yretc(mxcolret)
  REAL(P) :: hgtretc(nzret,mxcolret)
  REAL(P) :: theretc(nzret,mxcolret)
  REAL(P) :: obsret(nvar,nzret,mxcolret)
  REAL(P) :: odifret(nvar,nzret,mxcolret)
  REAL(P) :: qobsret(nvar,nzret,mxcolret)
  INTEGER :: qualret(nvar,nzret,mxcolret)
  INTEGER :: nlevret(mxcolret)
  INTEGER :: iret(mxcolret)
  INTEGER :: isrcret(0:mxret)
  INTEGER :: ncolret
!
!-----------------------------------------------------------------------
!
!  Input Analysis Control Variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=8) :: srcsng(nsrcsng)
  CHARACTER (LEN=8) :: srcua (nsrcua )
  CHARACTER (LEN=8) :: srcrad(nsrcrad)
  CHARACTER (LEN=8) :: srcret(nsrcret)

  INTEGER :: ianxtyp(mxpass)
  INTEGER :: iusesng(0:nsrcsng)
  INTEGER :: iuseua(0:nsrcua)
  INTEGER :: iuserad(0:nsrcrad)
  INTEGER :: iuseret(0:nsrcret)

  REAL(P) :: xyrange(mxpass)
  REAL(P) :: kpvrsq(nvar)
  REAL(P) :: wlim
  REAL(P) :: zrange(mxpass)
  REAL(P) :: zwlim
  REAL(P) :: thrng(mxpass)
  INTEGER :: iwstat
!
!-----------------------------------------------------------------------
!
!  Scratch Space
!
!-----------------------------------------------------------------------
!
  REAL(P) :: rngsqi(nvar)
  INTEGER :: knt(nvar,nz)
  REAL(P) :: wgtsum(nvar,nz)
  REAL(P) :: zsum(nvar,nz)
!
!-----------------------------------------------------------------------
!
!  Output Variables at Observation Locations
!
!-----------------------------------------------------------------------
!
  REAL(P) :: corsng(mxsng,nvar)
  REAL(P) :: corua(nzua,mxua,nvar)
  REAL(P) :: corrad(nzrdr,mxcolrad,nvarrad)
  REAL(P) :: corret(nzret,mxcolret,nvar)

  REAL(P) :: oanxsng(nvar,mxsng)
  REAL(P) :: oanxua(nvar,nzua,mxua)
  REAL(P) :: oanxrad(nvarrad,nzrdr,mxcolrad)
  REAL(P) :: oanxret(nvar,nzret,mxcolret)
!
!-----------------------------------------------------------------------
!
!  Output Grid
!
!-----------------------------------------------------------------------
!
  REAL(P) :: anx(nx,ny,nz,nvar)
!
!-----------------------------------------------------------------------
!
!  Work arrays
!
!-----------------------------------------------------------------------
!
  REAL(P) :: tem1(nx,ny,nz)
  REAL(P) :: tem2(nx,ny,nz)
  REAL(P) :: tem3(nx,ny,nz)
  REAL(P) :: tem4(nx,ny,nz)
  REAL(P) :: smcu(nx,ny,nz)
  REAL(P) :: smcv(nx,ny,nz)
  REAL(P) :: smcw(nx,ny,nz)

  REAL(P) :: u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz), wcont(nx,ny,nz)
  ! working arrays here
!  INTEGER :: ibegin(ny)
!  INTEGER :: iend(ny)
!
!-----------------------------------------------------------------------
!
!  Return status
!
!-----------------------------------------------------------------------
!
  INTEGER :: istatus
!
!-----------------------------------------------------------------------
!
!  Misc.local variables
!
!-----------------------------------------------------------------------
!
  REAL(P) :: ftabinv,setexp
  INTEGER :: i,j,k,isrc
  REAL(P) :: rpass,zrngsq,thrngsq

  INTEGER :: ibgn, iend, jbgn, jend
!
!   output:
!   ------
!    cfun:  value of cost function
!
!-----------------------------------------------------------------------
!
! argument
! --------
  INTEGER :: numctr,iflag,num
  REAL(P) :: ctrv(numctr)
  REAL(P) :: gdu_err(nx,ny,nz)
  REAL(P) :: gdv_err(nx,ny,nz)
  REAL(P) :: gdp_err(nx,ny,nz)
  REAL(P) :: gdt_err(nx,ny,nz)
  REAL(P) :: gdq_err(nx,ny,nz)
  REAL(P) :: gdw_err(nx,ny,nz)
  REAL(P) ::  gdscal(nx,ny,nz)


  REAL(P) ::   u_ctr(nx,ny,nz)
  REAL(P) ::   v_ctr(nx,ny,nz)
  REAL(P) ::   p_ctr(nx,ny,nz)
  REAL(P) ::   t_ctr(nx,ny,nz)
  REAL(P) ::   q_ctr(nx,ny,nz)
  REAL(P) ::   w_ctr(nx,ny,nz)
  REAL(P) ::     psi(nx,ny,nz)
  REAL(P) ::     phi(nx,ny,nz)
  REAL(P) ::     ffu(nx,ny,nz)
  REAL(P) ::     ffv(nx,ny,nz)
  REAL(P) ::  sinlat(nx,ny,nz)

!additional memory in model space
!
!  REAL(P), DIMENSION (:,:,:), allocatable :: u,v,w,wcont
!  REAL(P) :: zuu,zvv,zpp,ztt,zqq,zww
  DOUBLE PRECISION :: zuu,zvv,zpp,ztt,zqq,zww
!
!
!  REAL(P) :: cfun,f_b,f_bu,f_bv,f_bp,f_bt,f_bq,f_bw,ugrid,vgrid,vr
!  REAL(P) :: wgrid,sum1
  REAL(P) :: cfun,ugrid,vgrid,vr
  DOUBLE PRECISION :: f_b,f_bu,f_bv,f_bp,f_bt,f_bq,f_bw
  REAL(P) :: wgrid
!  REAL(P) :: f_osng,f_ousng,f_ovsng,f_opsng,f_otsng,f_oqsng,f_owsng
!  REAL(P) :: f_oua,f_ouua,f_ovua,f_opua,f_otua,f_oqua,f_owua
!  REAL(P) :: f_orad,f_ourad,f_ovrad,f_oprad,f_otrad,f_oqrad,f_owrad
  DOUBLE PRECISION :: f_osng,f_ousng,f_ovsng,f_opsng,f_otsng,f_oqsng,f_owsng
  DOUBLE PRECISION :: f_oua,f_ouua,f_ovua,f_opua,f_otua,f_oqua,f_owua
  DOUBLE PRECISION :: f_orad,f_ourad,f_ovrad,f_oprad,f_otrad,f_oqrad,f_owrad
  DOUBLE PRECISION :: f_smth_u, f_smth_v, f_smth_w
  DOUBLE PRECISION :: sum1

  DOUBLE PRECISION :: f_div, f_thermo
  REAL(P) :: div3(nx,ny,nz)
  REAL(P) :: maxdiv3, idiv,jdiv,kdiv
  REAL(P) :: tdiv1, tdiv2,tdiv3,tdiv4

! array save coast function of each variable, the nvar+1 is divergence and radar
  DOUBLE PRECISION :: cfun_single(nvar+1), cfun_total

!  REAL(P), ALLOCATABLE :: outlg(:,:,:)
!  INTEGER :: nxlg, nylg

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  nxlg = (nx-3)*nproc_x+3
!  nylg = (ny-3)*nproc_y+3

  ibgn = 1
  iend = nx-1
  jbgn = 1
  jend = ny-1

  IF (loc_x > 1)       ibgn = 2
  IF (loc_x < nproc_x) iend = nx-2

  IF (loc_y > 1)       jbgn = 2
  IF (loc_y < nproc_y) jend = ny-2

!
!
!  allocate control variable arrays
!  -------------------------------------------------------
!
  !ALLOCATE (     u(nx,ny,nz), STAT = istatus )
  !CALL check_alloc_status(istatus, "gradt:u")
  u(:,:,:) = 0.0
  !
  !ALLOCATE (     v(nx,ny,nz), STAT = istatus )
  !CALL check_alloc_status(istatus, "gradt:v")
  v(:,:,:) = 0.0
  !
  !ALLOCATE (     w(nx,ny,nz), STAT = istatus )
  !CALL check_alloc_status(istatus, "gradt:w")
  w(:,:,:) = 0.0
  !
  !ALLOCATE ( wcont(nx,ny,nz), STAT = istatus )
  !CALL check_alloc_status(istatus, "gradt:wcont")
  wcont(:,:,:) = 0.0

  DO i = 1, nvar+1
    cfun_single(i) = 0.0
  END DO

  !
  ! sngsw,uasw,radsw should not be changed from now on
  !
!  IF( sngsw > 0 ) THEN
!   DO j = 1, mxsng
!     DO i = 1, nvar
!       oanxsng(i,j) = 0.0
!     END DO
!   END DO
!  END IF
!
!
!  IF( uasw > 0 ) THEN
!    DO k = 1, mxua
!      DO j = 1, nzua
!        DO i = 1, nvar
!          oanxua(i,j,k) = 0.0
!        END DO
!      END DO
!    END DO
!  END IF

  IF( radsw > 0 ) THEN
    DO k = 1, mxcolrad
      DO j = 1, nzrdr
        DO i = 1, nvarrad
          oanxrad(i,j,k) = 0.0
        END DO
      END DO
    END DO
  END IF

  DO k = 1, nz
    DO j = 1, ny
      DO i = 1, nx
        u_ctr(i,j,k) = 0.
        v_ctr(i,j,k) = 0.
          psi(i,j,k) = 0.
          phi(i,j,k) = 0.
        p_ctr(i,j,k) = 0.
        t_ctr(i,j,k) = 0.
        q_ctr(i,j,k) = 0.
        w_ctr(i,j,k) = 0.
      END DO
    END DO
  END DO
!
! From ctrv to U, V, W, p_ctr, t_ctr, q_ctr etc.
!------------ctrv-----control variables--------one dimensional------------------
!

  CALL adtrans(numctr,nx,ny,nz,psi,phi,p_ctr,t_ctr,q_ctr,w_ctr,ctrv,tem4)
!
!-----------from here until div_opt all operations are on the scalar points----
!
  f_b = 0.

  f_bu = 0.
  f_bv = 0.
  f_bp = 0.
  f_bt = 0.
  f_bq = 0.
  f_bw = 0.

  DO k = 1,nz-1
    DO j = jbgn,jend             ! local sum
      DO i = ibgn,iend
        zuu =   psi(i,j,k) *   psi(i,j,k)
        zvv =   phi(i,j,k) *   phi(i,j,k)
        zpp = p_ctr(i,j,k) * p_ctr(i,j,k)
        ztt = t_ctr(i,j,k) * t_ctr(i,j,k)
        zqq = q_ctr(i,j,k) * q_ctr(i,j,k)
        zww = w_ctr(i,j,k) * w_ctr(i,j,k)

        f_bu = f_bu + zuu
        f_bv = f_bv + zvv
        f_bp = f_bp + zpp
        f_bt = f_bt + ztt
        f_bq = f_bq + zqq
        f_bw = f_bw + zww
      END DO
    END DO
  END DO

  f_b  = f_bu+f_bv+f_bp+f_bt+f_bq+f_bw   ! local sum.
  IF( 1 == 2 ) THEN
    f_b=0.0
    cfun_single(1)=0.0
    cfun_single(2)=0.0
    cfun_single(3)=0.0
    cfun_single(4)=0.0
    cfun_single(5)=0.0
    cfun_single(6)=0.0
  ELSE
    cfun_single(1)=f_bu
    cfun_single(2)=f_bv
    cfun_single(3)=f_bp
    cfun_single(4)=f_bt
    cfun_single(5)=f_bq
    cfun_single(6)=f_bw
  ENDIF
!
!  print*,'f_b===',f_b
!
! --------------------------------------------------------
  tem1=0.0;tem2=0.0;tem3=0.0;tem4=0.0      !working array
  CALL ctr_to_vbl(ipass_filt,hradius,radius_z,nx,ny,nz,gdu_err,        &
                  gdscal,  psi, tem1,tem2,tem3,tem4)

  tem1=0.0;tem2=0.0;tem3=0.0;tem4=0.0        
  CALL ctr_to_vbl(ipass_filt,hradius,radius_z,nx,ny,nz,gdv_err,        &
                  gdscal,  phi, tem1,tem2,tem3,tem4)

  tem1=0.0;tem2=0.0;tem3=0.0;tem4=0.0    
  CALL ctr_to_vbl(ipass_filt,hradius,radius_z,nx,ny,nz,gdp_err,        &
                  gdscal,p_ctr, tem1,tem2,tem3,tem4)

  tem1=0.0;tem2=0.0;tem3=0.0;tem4=0.0
  CALL ctr_to_vbl(ipass_filt,hradius,radius_z,nx,ny,nz,gdt_err,        &
                  gdscal,t_ctr, tem1,tem2,tem3,tem4)

  tem1=0.0;tem2=0.0;tem3=0.0;tem4=0.0
  CALL ctr_to_vbl(ipass_filt,hradius,radius_z,nx,ny,nz,gdq_err,        &
                  gdscal,q_ctr, tem1,tem2,tem3,tem4)

  tem1=0.0;tem2=0.0;tem3=0.0;tem4=0.0
  CALL ctr_to_vbl(ipass_filt,hradius,radius_z,nx,ny,nz,gdw_err,        &
                  gdscal,w_ctr, tem1,tem2,tem3,tem4)
  tem1=0.0;tem2=0.0;tem3=0.0;tem4=0.0
!
!
!  Option (cntl_var ==0),  u, v as control variables
!
!  u_ctr, v_ctr are scalars
!
  IF(cntl_var == 0) THEN

    DO k = 1, nz
      DO i = 1, nx
        DO j = 1, ny
          u_ctr(i,j,k) = psi(i,j,k)
          v_ctr(i,j,k) = phi(i,j,k)
        END DO
      END DO
    END DO

  ELSE

    DO k = 1, nz
      DO i = 2, nx-1
        DO j = 2, ny-1
          u_ctr(i,j,k) = ( psi(i-1,j+1,k)+psi(i,j+1,k)                  &
                          -psi(i-1,j-1,k)-psi(i,j-1,k) )/dy/4.          &
                        +( phi(i,  j,  k)-phi(i-1,j,k) )/dx
          u_ctr(i,j,k) = u_ctr(i,j,k)*mapfct(i,j,2)

          v_ctr(i,j,k) = ( psi(i+1,j-1,k)+psi(i+1,j,k)                  &
                          -psi(i-1,j-1,k)-psi(i-1,j,k) )/dx/4.          &
                        +( phi(i,  j,  k)-phi(i,j-1,k) )/dy
          v_ctr(i,j,k) = v_ctr(i,j,k)*mapfct(i,j,3)
        END DO
      END DO

      DO j=2,ny-1
        u_ctr( 1,j,k)=u_ctr( 2,j,k)+u_ctr( 2,j,k)-u_ctr( 3,j,k)
        v_ctr( 1,j,k)=v_ctr( 2,j,k)+v_ctr( 2,j,k)-v_ctr( 3,j,k)
        u_ctr(nx,j,k)=u_ctr(nx-1,j,k)+u_ctr(nx-1,j,k)-u_ctr(nx-2,j,k)
        v_ctr(nx,j,k)=v_ctr(nx-1,j,k)+v_ctr(nx-1,j,k)-v_ctr(nx-2,j,k)
      END DO

      DO i=2,nx-1
        u_ctr(i, 1,k)=u_ctr(i, 2,k)+u_ctr(i, 2,k)-u_ctr(i, 3,k)
        v_ctr(i, 1,k)=v_ctr(i, 2,k)+v_ctr(i, 2,k)-v_ctr(i, 3,k)
        u_ctr(i,ny,k)=u_ctr(i,ny-1,k)+u_ctr(i,ny-1,k)-u_ctr(i,ny-2,k)
        v_ctr(i,ny,k)=v_ctr(i,ny-1,k)+v_ctr(i,ny-1,k)-v_ctr(i,ny-2,k)
      END DO

      u_ctr(1,1 ,k)=0.5*( u_ctr(2,1,k)+u_ctr(1,2,k) )
      v_ctr(1,1 ,k)=0.5*( v_ctr(2,1,k)+v_ctr(1,2,k) )
      u_ctr(1,ny,k)=0.5*( u_ctr(1,ny-1,k)+u_ctr(2,ny,k) )
      v_ctr(1,ny,k)=0.5*( v_ctr(1,ny-1,k)+v_ctr(2,ny,k) )

      u_ctr(nx,1,k)=0.5*( u_ctr(nx-1,1,k)+u_ctr(nx,2,k) )
      v_ctr(nx,1,k)=0.5*( v_ctr(nx-1,1,k)+v_ctr(nx,2,k) )
      u_ctr(nx,ny,k)=0.5*( u_ctr(nx,ny-1,k)+u_ctr(nx-1,ny,k) )
      v_ctr(nx,ny,k)=0.5*( v_ctr(nx,ny-1,k)+v_ctr(nx-1,ny,k) )

    END DO

    IF (mp_opt > 0) THEN
      CALL acct_interrupt(mp_acct)
      CALL mpsendrecv2dew(u_ctr, nx, ny, nz, ebc, wbc, 0, tem4)
      CALL mpsendrecv2dns(u_ctr, nx, ny, nz, nbc, sbc, 0, tem4)

      CALL mpsendrecv2dew(v_ctr, nx, ny, nz, ebc, wbc, 0, tem4)
      CALL mpsendrecv2dns(v_ctr, nx, ny, nz, nbc, sbc, 0, tem4)
      CALL acct_stop_inter
    END IF

  END IF
!
!  loading single level data
! --------------------------------------------------------
!
  f_osng  = 0.
  f_ousng = 0.
  f_ovsng = 0.
  f_opsng = 0.
  f_otsng = 0.
  f_oqsng = 0.
  f_owsng = 0.

!  IF (myproc == 0) WRITE(ount,'(1x,a,i2)') '==== sngsw ===== ',sngsw

  IF(sngsw > 0) THEN

    CALL linearint_3d(nx,ny,nz,u_ctr(1,1,1),xs(1),ys(1),zs(1,1,1),      &
           usesng,1, mxsng, icatsng, nobsng,                            &
           1,xsng_p,ysng_p,zsng_1,zsng_2,hgtsng,ihgtsng,corsng(1,1) )

    CALL linearint_3d(nx,ny,nz,v_ctr(1,1,1),xs(1), ys(1),zs(1,1,1),     &
           usesng,1, mxsng, icatsng, nobsng,                            &
           2,xsng_p,ysng_p,zsng_1,zsng_2,hgtsng,ihgtsng,corsng(1,2) )

    CALL linearint_3d(nx,ny,nz,p_ctr(1,1,1),xs(1),ys(1),zs(1,1,1),      &
           usesng,1, mxsng, icatsng, nobsng,                            &
           3,xsng_p,ysng_p,zsng_1,zsng_2,hgtsng,ihgtsng,corsng(1,3) )

    CALL linearint_3d(nx,ny,nz,t_ctr(1,1,1),xs(1),ys(1),zs(1,1,1),      &
           usesng,1, mxsng, icatsng, nobsng,                            &
           4,xsng_p,ysng_p,zsng_1,zsng_2,hgtsng,ihgtsng,corsng(1,4) )

    CALL linearint_3d(nx,ny,nz,q_ctr(1,1,1),xs(1),ys(1),zs(1,1,1),      &
           usesng,1, mxsng, icatsng, nobsng,                            &
           5,xsng_p,ysng_p,zsng_1,zsng_2,hgtsng,ihgtsng,corsng(1,5) )

!     call linearint_3d(nx,ny,nz,w_ctr(1,1,1),xs(1),ys(1),zs(1,1,1),
!        1, mxsng, icatsng, nobsng,                                   &
!        6,xsng_p,ysng_p,zsng_1,zsng_2,hgtsng,ihgtsng,corsng(1,6) )
!
    num=0
    DO i = 1, nobsng
      IF ( usesng(i) ) THEN
        iflag = 1
        IF(qualsng(1,i) <= 0 .or. ihgtsng(i)<0 ) iflag = 0
        zuu  = corsng(i,1) - odifsng(1,i)
        IF(iflag /= 0) THEN
          corsng(i,1) = iflag*zuu / ( qobsng2(1,i) )
        ELSE
          corsng(i,1) = 0.0
        END IF

        zuu     = zuu * corsng(i,1)
        IF (ownsng(i)) f_ousng = f_ousng +  zuu
!
!
!    if(iflag.eq.1) num=num+1
!    print*,' num==',num,nobsng
!    print*,'qobsng2==', qobsng2(1,i),odifsng(1,i),iflag
!
!
        iflag = 1
        IF(qualsng(2,i) <= 0 .or. ihgtsng(i)<0 ) iflag = 0
        zvv = corsng(i,2) - odifsng(2,i)
        IF(iflag /= 0) THEN
          corsng(i,2) = iflag*zvv / ( qobsng2(2,i) )
        ELSE
          corsng(i,2) = 0.0
        END IF
        zvv     = zvv * corsng(i,2)
        IF (ownsng(i)) f_ovsng = f_ovsng +  zvv
!
!    if(iflag.eq.1) num=num+1
!    print*,' num==',num,nobsng
!    print*,'corsng==', corsng(i,2),odifsng(2,i)
!
!
        iflag = 1
        IF(qualsng(3,i) <= 0 .or. ihgtsng(i)<0) iflag = 0
        zpp        = corsng(i,3) - odifsng(3,i)
        IF(iflag /= 0) THEN
          corsng(i,3) = iflag*zpp / ( qobsng2(3,i) )
        ELSE
          corsng(i,3) = 0.0
        END IF
        zpp        = zpp * corsng(i,3)
        IF (ownsng(i)) f_opsng    = f_opsng +  zpp
!
!    if(iflag.eq.1) num=num+1
!    print*,' num==',num,nobsng
!    print*,'corsng==', corsng(i,3),odifsng(3,i)
!    print*,'qobsng2==', qobsng2(3,i),odifsng(3,i),iflag
!    print*,'f_opsng====',qobsng2(1,i),qobsng2(2,i),
!    :                       qobsng2(3,i),
!    :                       qobsng2(4,i),qobsng2(5,i)
!
        iflag = 1
        IF(qualsng(4,i) <= 0 .or. ihgtsng(i)<0) iflag = 0
        ztt        = corsng(i,4) - odifsng(4,i)
        IF(iflag /= 0) THEN
          corsng(i,4) = iflag*ztt / ( qobsng2(4,i) )
        ELSE
          corsng(i,4) = 0.0
        END IF
        ztt        = ztt * corsng(i,4)
        IF (ownsng(i)) f_otsng    = f_otsng +  ztt
!
!    if(iflag.eq.1) num=num+1
!    print*,' num==',num,nobsng
!    print*,'corsng==', corsng(i,4),odifsng(4,i)
!
        iflag = 1
        IF(qualsng(5,i) <= 0 .or. ihgtsng(i)<0) iflag = 0
        zqq        = corsng(i,5) - odifsng(5,i)
        IF(iflag /= 0) THEN
          corsng(i,5) = iflag*zqq / ( qobsng2(5,i) )
        ELSE
          corsng(i,5) = 0.0
        END IF
        zqq        = zqq * corsng(i,5)
        IF (ownsng(i)) f_oqsng    = f_oqsng +  zqq
!
!
!!   if(iflag.eq.1) num=num+1
!!   print*,'corsng==', corsng(i,5),odifsng(5,i)
!!   print*,' qobsng2 odifsng=',qobsng2(5,i),odifsng(5,i)
!    :            ,qualsng(5,i),iflag
!
!
!
!    if(1.eq.0) THEN
!
!       iflag = 1
!    if(qualsng(6,i).le.0 .and. ihgtsng<0) iflag = 0
!    zww        = corsng(i,6) - odifsng(6,i)
!    if(iflag.eq.1) num=num+1
!    print*,' num==',num,nobsng
!    if(iflag.ne.0) then
!      corsng(i,6) = iflag*zww / ( qobsng2(6,i) )
!    else
!      corsng(i,6) = 0.0
!    end if
!    zww        = zww * corsng(i,6)
!    f_owsng    = f_owsng +  zww
!    print*,'obsng==', corsng(i,6),odifsng(6,i)
!      end if
!
!  observation cost function for single level data
!  --------------------------------------------------------------
       if (ownsng(i)) num = num + 1

      END IF  ! usesng(i)

    END DO

    f_osng=f_ousng+f_ovsng+f_opsng+f_otsng+f_oqsng+f_owsng
    cfun_single(1)=cfun_single(1)+f_ousng
    cfun_single(2)=cfun_single(2)+f_ovsng
    cfun_single(3)=cfun_single(3)+f_opsng
    cfun_single(4)=cfun_single(4)+f_otsng
    cfun_single(5)=cfun_single(5)+f_oqsng
    cfun_single(6)=cfun_single(6)+f_owsng

!    call mpsumdp(f_ousng,1)
!    call mpsumdp(f_ovsng,1)
!    call mpsumdp(f_opsng,1)
!    call mpsumdp(f_otsng,1)
!    call mpsumdp(f_oqsng,1)
!    call mpsumdp(f_owsng,1)
!    call mpsumdp( f_osng,1)
!
!    IF (myproc == 0) THEN
!      write(ount,*) '======================='
!      write(ount,*) '  npass=',npass,' f_ousng=', f_ousng
!      write(ount,*) '  npass=',npass,' f_ovsng=', f_ovsng
!      write(ount,*) '  npass=',npass,' f_opsng=', f_opsng
!      write(ount,*) '  npass=',npass,' f_otsng=', f_otsng
!      write(ount,*) '  npass=',npass,' f_oqsng=', f_oqsng
!      write(ount,*) '  npass=',npass,' f_owsng=', f_owsng
!
!      write(ount,*) '  npass=',npass,' f_osng =', f_osng
!      write(ount,*) '======================='
!      call flush(ount)
!    END IF

!    WRITE(ount,*) myproc,':  npass=',npass,' f_ousng=', f_ousng
!    WRITE(ount,'(5x,a,I2.2,a,I6,a)') 'Domain ',myproc,' processed ',num,' SNG data.'
!    call mpbarrier

  END IF
!
!  loading upper level data
! --------------------------------------------------------
!
  num = 0
  f_oua  = 0.
  f_ouua = 0.
  f_ovua = 0.
  f_opua = 0.
  f_otua = 0.
  f_oqua = 0.
  f_owua = 0.

!  IF (myproc == 0) WRITE(ount,'(1x,a,i2)') '==== uasw ===== ',uasw

  IF(uasw > 0) THEN

    CALL linearint_3d(nx,ny,nz,u_ctr(1,1,1),xs(1),ys(1),zs(1,1,1),      &
           useua, nzua, mxua, nlevsua, nobsua,                          &
           1,xua_p,yua_p,zua_1,zua_2,hgtua,ihgtua,corua(1,1,1) )

    CALL linearint_3d(nx,ny,nz,v_ctr(1,1,1),xs(1),ys(1),zs(1,1,1),      &
           useua, nzua, mxua, nlevsua, nobsua,                          &
           2,xua_p,yua_p,zua_1,zua_2,hgtua,ihgtua,corua(1,1,2) )

    CALL linearint_3d(nx,ny,nz,p_ctr(1,1,1),xs(1),ys(1),zs(1,1,1),      &
           useua, nzua, mxua, nlevsua, nobsua,                          &
           3,xua_p,yua_p,zua_1,zua_2,hgtua,ihgtua,corua(1,1,3) )

    CALL linearint_3d(nx,ny,nz,t_ctr(1,1,1),xs(1),ys(1),zs(1,1,1),      &
           useua, nzua, mxua, nlevsua, nobsua,                          &
           4,xua_p,yua_p,zua_1,zua_2,hgtua,ihgtua,corua(1,1,4) )

    CALL linearint_3d(nx,ny,nz,q_ctr(1,1,1),xs(1),ys(1),zs(1,1,1),      &
           useua, nzua, mxua, nlevsua, nobsua,                          &
           5,xua_p,yua_p,zua_1,zua_2,hgtua,ihgtua,corua(1,1,5) )
!
!   CALL linearint_3d(nx,ny,nz,w_ctr(1,1,1),xs(1),ys(1),zs(1,1,1),      &
!          nzua, mxua, nlevsua, nobsua,                                 &
!          6,xua_p,yua_p,zua_1,zua_2,hgtua,ihgtua,corua(1,1,6) )
!
!
    DO i = 1, nobsua
      IF ( useua(i) ) THEN
        DO j = 1, nlevsua(i)

          iflag = 1
          IF(qualua(1,j,i) <= 0 .OR. ihgtua(j,i)<0) iflag = 0
          zuu          = corua(j,i,1) - odifua(1,j,i)
          IF(iflag /= 0) THEN
            corua(j,i,1) = iflag*zuu /( qobsua2(1,j,i) )
          ELSE
            corua(j,i,1) = 0.0
          END IF
          zuu          = zuu * corua(j,i,1)
          IF (ownua(i)) f_ouua = f_ouua +  zuu

          iflag = 1
          IF(qualua(2,j,i) <= 0 .or. ihgtua(j,i)<0) iflag = 0
          zvv          = corua(j,i,2) - odifua(2,j,i)
          IF(iflag /= 0) THEN
            corua(j,i,2) = iflag*zvv/( qobsua2(2,j,i) )
          ELSE
            corua(j,i,2) = 0.0
          END IF
          zvv          = zvv * corua(j,i,2)
          IF (ownua(i)) f_ovua = f_ovua +  zvv

          iflag = 1
          IF(qualua(3,j,i) <= 0 .or. ihgtua(j,i)<0) iflag = 0
          zpp          = corua(j,i,3) - odifua(3,j,i)
          IF(iflag /= 0) THEN
            corua(j,i,3) = iflag*zpp /( qobsua2(3,j,i) )
          ELSE
            corua(j,i,3) = 0.0
          END IF
          zpp          = zpp * corua(j,i,3)
          IF (ownua(i)) f_opua = f_opua +  zpp

          iflag = 1
          IF(qualua(4,j,i) <= 0 .or. ihgtua(j,i)<0) iflag = 0
          ztt          = corua(j,i,4) - odifua(4,j,i)
          IF(iflag /= 0) THEN
            corua(j,i,4) = iflag*ztt /( qobsua2(4,j,i) )
          ELSE
            corua(j,i,4) = 0.0
          END IF
          ztt          = ztt * corua(j,i,4)
          IF (ownua(i)) f_otua = f_otua +  ztt


          iflag = 1
          IF(qualua(5,j,i) <= 0 .or. ihgtua(j,i)<0) iflag = 0
          zqq          = corua(j,i,5) - odifua(5,j,i)
          IF(iflag /= 0) THEN
            corua(j,i,5) = iflag*zqq / ( qobsua2(5,j,i) )
          ELSE
            corua(j,i,5) = 0.0
          END IF
          zqq          = zqq * corua(j,i,5)
          IF (ownua(i)) f_oqua  = f_oqua +  zqq

!
!    print*,' qobsua2 odifua=',qobsua2(5,k,i),odifua(5,k,i)
!   :            ,qualua(5,k,i),iflag
!
!         iflag = 1
!      if(qualua(6,k,i).le.0 .or. ihgtua(j,i)<0) iflag = 0
!      zww          = corua(j,i,6) - odifua(6,j,i)
!      if(iflag.ne.0) then
!        corua(j,i,6) = iflag*zww / ( qobsua2(6,j,i) )
!      else
!        corua(j,i,6) = 0.0
!      end if
!      zww          = zww * corua(j,i,6)
!      f_owua       = f_owua +  zww
!

        END DO
        IF (ownua(i)) num = num + 1
      END IF
    END DO

    f_oua = f_ouua+f_ovua+f_opua+f_otua+f_oqua+f_owua
    cfun_single(1)=cfun_single(1)+f_ouua
    cfun_single(2)=cfun_single(2)+f_ovua
    cfun_single(3)=cfun_single(3)+f_opua
    cfun_single(4)=cfun_single(4)+f_otua
    cfun_single(5)=cfun_single(5)+f_oqua
    cfun_single(6)=cfun_single(6)+f_owua

!    IF (myproc == 0) THEN
!    WRITE(ount,*) myproc, 'f_oua = ',f_oua
!    print*,'f_ouua = ',f_ouua
!    print*,'f_ovua = ',f_ovua
!    print*,'f_opua = ',f_opua
!    print*,'f_otua = ',f_otua
!    print*,'f_oqua = ',f_oqua
!    print*,'f_owua = ',f_owua
!      call flush(ount)
!    END IF
!    WRITE(ount,'(5x,a,I2.2,a,I6,a)') 'Domain ',myproc,' processed ',num,' UA data.'


  END IF
!
!  loading radar data
! --------------------------------------------------------
!
!
  num = 0
  f_orad  = 0.
  f_ovrad = 0.
  f_oqrad = 0.
  f_div   = 0.

!  IF (myproc == 0) WRITE(ount,'(1x,a,i2)') '==== radsw ===== ',radsw

  IF(radsw > 0) THEN

    CALL linearint_3d(nx,ny,nz,u_ctr(1,1,1),xs(1),ys(1),zs(1,1,1),      &
           userad, nzrdr, mxcolrad, nlevrad, ncolrad,                   &
           1,xradc_p,yradc_p,zradc_1,zradc_2,hgtradc,                   &
                                            ihgtradc,corrad(1,1,1) )

    CALL linearint_3d(nx,ny,nz,v_ctr(1,1,1),xs(1),ys(1),zs(1,1,1),      &
           userad, nzrdr, mxcolrad, nlevrad, ncolrad,                   &
           2,xradc_p,yradc_p,zradc_1,zradc_2,hgtradc,                   &
                                            ihgtradc,corrad(1,1,2) )

    CALL linearint_3d(nx,ny,nz,w_ctr(1,1,1),xs(1),ys(1),zs(1,1,1),      &
           userad, nzrdr, mxcolrad, nlevrad, ncolrad,                   &
           6,xradc_p,yradc_p,zradc_1,zradc_2,hgtradc,                   &
                                            ihgtradc,corrad(1,1,3) )

    DO i = 1, ncolrad
      IF(iuserad(isrcrad(irad(i))) > 0 .AND. userad(i) ) THEN

        DO k=1,nlevrad(i)

          IF(qualrad(2,k,i) > 0 .AND.ihgtradc(k,i)>=0 ) THEN
!           abs(odifrad(2,k,i))< 25.0 ) THEN
!           abs(odifrad(2,k,i))< 20.0 .AND. abs(obsrad(2,k,i))> 0.05 ) THEN

            vr = ( uazmrad(i)*corrad(k,i,1)                             &
                +vazmrad(i)*corrad(k,i,2) ) * dsdr(k,i)                 &
                           +corrad(k,i,3)  * dhdr(k,i)

            zvv        = vr - odifrad(2,k,i)
!
!         corrad(k,i,2)=zvv /( qobsrad2(2,k,i) )
!         zvv        = zvv * corrad(k,i,2)
!
            oanxrad(2,k,i)=zvv /( qobsrad2(2,k,i) )
            zvv        = zvv * oanxrad(2,k,i)

            IF (ownrad(i)) f_ovrad = f_ovrad +  zvv

          END IF  ! (qualrad(2,k,i) > 0 )
!!
!!      IF(1 == 0) THEN
!!
!!         CALL linearint_3d(nx,ny,nz,q_ctr(1,1,1),xs(1),ys(1),zs(1,1,1), &
!!              5,xradc(i),yradc(i),hgtradc(k,i),oanxrad(1,k,i),iflag )
!!         IF(qualrad(3,k,i) <= 0) iflag = 0
!!_gao   zqq        = oanxrad(1,k,i) - odifrad(1,k,i)
!!         zqq        = oanxrad(1,k,i) - odifrad(1,k,i)
!!         IF(iflag /= 0) THEN
!!_gao     oanxrad(1,k,i)=iflag*zqq /(qobsrad2(1,k,i)*qobsrad2(1,k,i))
!!           oanxrad(1,k,i)=iflag*zqq /( qobsrad2(1,k,i) )
!!         ELSE
!!           oanxrad(1,k,i)=0.0
!!         END IF
!!         zqq        = zqq * oanxrad(1,k,i)
!!         f_oqrad    = f_oqrad +  zqq
!!       END IF
!!
        END DO  ! nlevrad(i)

        IF (ownrad(i)) num = num +1
      END IF
    END DO

!    WRITE(ount,'(5x,a,I2.2,a,I6,a)') 'Domain ',myproc,' processed ',num,' RAD data.'
!
!-----------------------------------------------------------------------
!
!  Compute the momentum divergence term, defined as
!
!  div = d(u*rhostr)/dx + d(v*rhostr)/dy + d(wcont*rhostr)/dz.
!
!  NOTE: u, v & w should be on vector points when using in divergence constraint
!        and all other arps model constraints in the future.
!
!-----------------------------------------------------------------------
!
    IF( div_opt== 1 ) THEN

      ! First, convert u_ctr, v_ctr, w_ctr from scalar grid to vector grid
      ! of u, v, w

      DO k= 1,nz-1
        DO j= 1,ny-1
          DO i= 2,nx-1
            u(i,j,k) = 0.5*(anx(i,j,k,1)+u_ctr(i,j,k)+anx(i-1,j,k,1)+u_ctr(i-1,j,k))
          END DO
        END DO
      END DO

      DO k= 1,nz-1
        DO j= 1,ny-1
          u(1,j,k)  = u(2,j,k)
          u(nx,j,k) = u(nx-1,j,k)
        END DO
      END DO

      DO k= 1,nz-1
        DO j= 2,ny-1
          DO i= 1,nx-1
            v(i,j,k) = 0.5*(anx(i,j,k,2)+v_ctr(i,j,k)+anx(i,j-1,k,2)+v_ctr(i,j-1,k))
          END DO
        END DO
      END DO

      DO k= 1,nz-1
        DO i= 1,nx-1
          v(i,1,k)  = v(i,2,k)
          v(i,ny,k) = v(i,ny-1,k)
        END DO
      END DO

      DO k=2,nz-1
        DO j= 1,ny-1
          DO i= 1,nx-1
            w(i,j,k) = 0.5*(anx(i,j,k,6)+w_ctr(i,j,k) + anx(i,j,k-1,6)+w_ctr(i,j,k-1))
          END DO
        END DO
      END DO

      DO j= 1,ny-1
        DO i= 1,nx-1
          w(i,j,1)  = w(i,j,2)
          w(i,j,nz) = w(i,j,nz-1)
        END DO
      END DO

      IF (mp_opt > 0) THEN
        !CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(u, nx, ny, nz, ebc, wbc, 1, tem4)
        CALL mpsendrecv2dns(u, nx, ny, nz, nbc, sbc, 1, tem4)

        CALL mpsendrecv2dew(v, nx, ny, nz, ebc, wbc, 2, tem4)
        CALL mpsendrecv2dns(v, nx, ny, nz, nbc, sbc, 2, tem4)

        !CALL mpsendrecv2dew(w, nx, ny, nz, ebc, wbc, 3, tem4)
        !CALL mpsendrecv2dns(w, nx, ny, nz, nbc, sbc, 3, tem4)
        !CALL acct_stop_inter
      END IF

      CALL wcontra(nx,ny,nz,u,v,w,mapfct,j1,j2,j3,aj3z,                 &
                   rhostr,rhostru,rhostrv,rhostrw,wcont,tem1,tem2)
      !
      ! tem1 is on U grid, tem2 on V grid, tem3 on W grid, and div3 is on scalar grid
      !
      IF(wgt_div_h > 0.0 ) THEN
        DO k=1,nz-1
          DO j=1,ny-1
            DO i=1,nx
              tem1(i,j,k)=u(i,j,k)*rhostru(i,j,k)*mapfct(i,j,5)
            END DO
          END DO
        END DO

        DO k=1,nz-1
          DO j=1,ny
            DO i=1,nx-1
              tem2(i,j,k)=v(i,j,k)*rhostrv(i,j,k)*mapfct(i,j,6)
            END DO
          END DO
        END DO
      END IF

      IF(wgt_div_v > 0.0 ) THEN
        DO k=1,nz
          DO j=1,ny-1
            DO i=1,nx-1
             tem3(i,j,k)=wcont(i,j,k)*rhostrw(i,j,k)
            END DO
          END DO
        END DO
      END IF

      IF( (wgt_div_h > 0.0) .AND. (wgt_div_v > 0.0)) THEN
        DO k=1,nz-1
          DO j=jbgn,jend
            DO i=ibgn,iend
              div3(i,j,k) = 1./wgt_div_h * j3inv(i,j,k)                 &
                         * ( mapfct(i,j,7)                              &
                         * ((tem1(i+1,j,k)-tem1(i,j,k))*dxinv           &
                          +(tem2(i,j+1,k)-tem2(i,j,k))*dyinv))          &
                         +1./wgt_div_v * j3inv(i,j,k)                   &
                         * ((tem3(i,j,k+1)-tem3(i,j,k))*dzinv )
              f_div = f_div + div3(i,j,k) * div3(i,j,k)
            END DO
          END DO
        END DO
      ELSE IF( (wgt_div_h > 0.0) .AND. (wgt_div_v < 0.0)) THEN
        DO k=1,nz-1
          DO j=jbgn,jend
            DO i=ibgn,iend
              div3(i,j,k) = 1./wgt_div_h * j3inv(i,j,k) * ( mapfct(i,j,7) &
                            * ((tem1(i+1,j,k)-tem1(i,j,k))*dxinv          &
                              +(tem2(i,j+1,k)-tem2(i,j,k))*dyinv))
              f_div = f_div + div3(i,j,k) * div3(i,j,k)
            END DO
          END DO
        END DO
      END IF

      IF (mp_opt > 0) THEN
        !CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(div3, nx, ny, nz, ebc, wbc, 0, tem4)
        CALL mpsendrecv2dns(div3, nx, ny, nz, nbc, sbc, 0, tem4)
        !CALL acct_stop_inter
      END IF

    END IF   ! end of the divergence constraint

  END IF  ! end of  the radar data loop

  f_orad = f_ovrad+f_oqrad+f_div
  cfun_single(nvar+1)=f_orad
!  print*,'f_ovrad f_oqrad f_div==',f_ovrad,f_oqrad,f_div,wgt_div_h,wgt_div_v
!
!
!
!-----------------------------------------------------------------------
!
!  Compute the thermowind constraint term, defined as
!
!  first: du/dz +(g/f/anx(5))*d(anx(5)/dy -u/anx(5)*d(anx(5))/dz
!  first: dv/dz -(g/f/anx(5))*d(anx(5)/dx -v/anx(5)*d(anx(5))/dz
!
!
!-----------------------------------------------------------------------
!
!
  f_thermo = 0.0

  IF ( thermo_opt==1 ) then
    !
    ! u, v, w used as temporary arrays on scalar grid, where w is used for temperature here
    !
    DO k= 1,nz-1
      DO j= 1,ny-1
        DO i= 1,nx-1
          u(i,j,k) =  anx(i,j,k,1)+ u_ctr(i,j,k)
          v(i,j,k) =  anx(i,j,k,2)+ v_ctr(i,j,k)
          w(i,j,k) = (anx(i,j,k,4)+ t_ctr(i,j,k) )                     &
                   * ((anx(i,j,k,3)+p_ctr(i,j,k) )/p0 )**rddcp
        END DO
      END DO
    END DO

    ! Added by WYH for application in MPI mode
    DO j = 1,ny-1          ! Set top boundary (W points) since (k+1) is used below
      DO i = 1,nx-1        ! Assume zero normal gradient condition.
        u(i,j,nz) = u(i,j,nz-1)
        v(i,j,nz) = v(i,j,nz-1)
        w(i,j,nz) = w(i,j,nz-1)
      END DO
    END DO

    DO k = 1,nz            ! Set east boundary (U points) since (i+1) is used below
      DO j = 1,ny-1        ! Assume zero normal gradient condition.
        w(nx,j,k) = w(nx-2,j,k)
      END DO
    END DO

    DO k = 1,nz            ! Set north boundary (V points) since (j+1) is used below
      DO i = 1,nx          ! Assume zero normal gradient condition.
        w(i,ny,k) = w(i,ny-2,k)
      END DO
    END DO
    ! End of addition from WYH

    !IF (mp_opt > 0) THEN
    !  !CALL acct_interrupt(mp_acct)
    !  CALL mpsendrecv2dew(u, nx, ny, nz, ebc, wbc, 0, tem4)
    !  CALL mpsendrecv2dns(u, nx, ny, nz, nbc, sbc, 0, tem4)
    !
    !  CALL mpsendrecv2dew(v, nx, ny, nz, ebc, wbc, 0, tem4)
    !  CALL mpsendrecv2dns(v, nx, ny, nz, nbc, sbc, 0, tem4)
    !
    !  CALL mpsendrecv2dew(w, nx, ny, nz, ebc, wbc, 0, tem4)
    !  CALL mpsendrecv2dns(w, nx, ny, nz, nbc, sbc, 0, tem4)
    !  !CALL acct_stop_inter
    !END IF

    tem1(:,:,:) = 0.0
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          ! Note that w is temperature here
          ffu(i,j,k) = (u(i,j,k+1)-u(i,j,k))*dzinv*j3inv(i,j,k)              &
                      +g/sinlat(i,j,k)/w(i,j,k)*(w(i,j+1,k)-w(i,j,k))*dyinv  &
                      -u(i,j,k)/w(i,j,k)*(w(i,j,k+1)-w(i,j,k))*dzinv*j3inv(i,j,k)

          ffv(i,j,k) = (v(i,j,k+1)-v(i,j,k))*dzinv*j3inv(i,j,k)              &
                      -g/sinlat(i,j,k)/w(i,j,k)*(w(i+1,j,k)-w(i,j,k))*dxinv  &
                      -v(i,j,k)/w(i,j,k)*(w(i,j,k+1)-w(i,j,k))*dzinv*j3inv(i,j,k)
        END DO
      END DO
    END DO

    IF (mp_opt > 0) THEN
      !CALL acct_interrupt(mp_acct)
      !CALL mpsendrecv2dew(ffu, nx, ny, nz, ebc, wbc, 0, tem4)
      CALL mpsendrecv2dns(ffu, nx, ny, nz, nbc, sbc, 0, tem4)
      CALL mpsendrecv2dew(ffv, nx, ny, nz, ebc, wbc, 0, tem4)
      !CALL mpsendrecv2dns(ffv, nx, ny, nz, nbc, sbc, 0, tem4)
      !CALL acct_stop_inter
    END IF

    DO k = 1,nz-1
      DO j = jbgn,jend
        DO i = ibgn,iend
          f_thermo = f_thermo                                           &
                    +wgt_thermo*(ffu(i,j,k)*ffu(i,j,k)+ffv(i,j,k)*ffv(i,j,k))
        END DO
      END DO
    END DO
    cfun_single(nvar+1)=cfun_single(nvar+1)+f_thermo
  END IF

  f_smth_u = 0.0
  f_smth_v = 0.0
  f_smth_w = 0.0
  IF(smth_opt == 1) THEN

    ! Added by WYH to define and set boundary for u,v & w, otherwise,
    ! the definitions of u, v & w are ambiguous depending on values of thermo_opt, radsw, div_opt etc.
    DO k= 1,nz-1
      DO j= 1,ny-1
        DO i= 1,nx-1
          u(i,j,k) = anx(i,j,k,1)+ u_ctr(i,j,k)
          v(i,j,k) = anx(i,j,k,2)+ v_ctr(i,j,k)
          w(i,j,k) = anx(i,j,k,6)+ w_ctr(i,j,k)
        END DO
      END DO
    END DO

    DO j = 1,ny-1          ! Set top boundary (W points) since (k+1) is used below
      DO i = 1,nx-1        ! Assume zero normal gradient condition.
        u(i,j,nz) = u(i,j,nz-1)
        v(i,j,nz) = v(i,j,nz-1)
        w(i,j,nz) = w(i,j,nz-1)
      END DO
    END DO

    DO k = 1,nz            ! Set east boundary (U points) since (i+1) is used below
      DO j = 1,ny-1        ! Assume zero normal gradient condition.
        u(nx,j,k) = u(nx-2,j,k)
        v(nx,j,k) = v(nx-2,j,k)
        w(nx,j,k) = w(nx-2,j,k)
      END DO
    END DO

    DO k = 1,nz            ! Set north boundary (V points) since (j+1) is used below
      DO i = 1,nx          ! Assume zero normal gradient condition.
        u(i,ny,k) = u(i,ny-2,k)
        v(i,ny,k) = v(i,ny-2,k)
        w(i,ny,k) = w(i,ny-2,k)
      END DO
    END DO
    ! End of Adding from WYH

    DO k=2,nz-1
      DO j=2,ny-1
        DO i=2,nx-1
          smcu(i,j,k)= (u(i+1,j,k)+u(i-1,j,k)+u(i,j+1,k)+u(i,j-1,k)     &
                       +u(i,j,k+1)+u(i,j,k-1)-6.*u(i,j,k))
          smcv(i,j,k)= (v(i+1,j,k)+v(i-1,j,k)+v(i,j+1,k)+v(i,j-1,k)     &
                       +v(i,j,k+1)+v(i,j,k-1)-6.*v(i,j,k))
          smcw(i,j,k)= (w(i+1,j,k)+w(i-1,j,k)+w(i,j+1,k)+w(i,j-1,k)     &
                       +w(i,j,k+1)+w(i,j,k-1)-6.*w(i,j,k))
        END DO
      END DO
    END DO

    DO j=2,ny-1
      DO i=2,nx-1
        smcu(i,j,1)= (u(i+1,j,1)+u(i-1,j,1)+u(i,j+1,1)+u(i,j-1,1)       &
                     -4.*u(i,j,1))
        smcv(i,j,1)= (v(i+1,j,1)+v(i-1,j,1)+v(i,j+1,1)+v(i,j-1,1)       &
                     -4.*v(i,j,1))
        smcw(i,j,1)= (w(i,j,1)+w(i-1,j,1)+w(i,j+1,1)+w(i,j-1,1)         &
                     -4.*w(i,j,1))
      END DO
    END DO

    IF (mp_opt > 0) THEN
      !CALL acct_interrupt(mp_acct)
      CALL mpsendrecv2dew(smcu, nx, ny, nz, ebc, wbc, 0, tem4)
      CALL mpsendrecv2dns(smcu, nx, ny, nz, nbc, sbc, 0, tem4)

      CALL mpsendrecv2dew(smcv, nx, ny, nz, ebc, wbc, 0, tem4)
      CALL mpsendrecv2dns(smcv, nx, ny, nz, nbc, sbc, 0, tem4)

      CALL mpsendrecv2dew(smcw, nx, ny, nz, ebc, wbc, 0, tem4)
      CALL mpsendrecv2dns(smcw, nx, ny, nz, nbc, sbc, 0, tem4)
      !CALL acct_stop_inter
    END IF

    DO k=1,nz-1
      DO j=2,jend
        DO i=2,iend
          f_smth_u= f_smth_u + wgt_smth*smcu(i,j,k)**2
          f_smth_v= f_smth_v + wgt_smth*smcv(i,j,k)**2
          f_smth_w= f_smth_w + wgt_smth*smcw(i,j,k)**2
        END DO
      END DO
    END DO

    cfun_single(1) = cfun_single(1) + f_smth_u
    cfun_single(2) = cfun_single(2) + f_smth_v
    cfun_single(6) = cfun_single(6) + f_smth_w

  END IF

  DO i=1,nvar+1
    cfun_single(i)=cfun_single(i)/2.0
  END DO

!  call mpsumdp(cfun_single,nvar+1)
!  call mpsumdp(f_ovrad,1)
!  call mpsumdp(f_oqrad,1)
!  call mpsumdp(f_div,1)
!  call mpsumdp(f_orad,1)
!  IF (myproc == 0) THEN
!    DO i=1,nvar+1
!      WRITE(6,*) 'cfun_single(',i,') = ',cfun_single(i)
!    END DO
!    write(6,*) 'f_orad  = ',f_orad
!    write(6,*) 'f_ovrad = ',f_ovrad
!    write(6,*) 'f_oqrad = ',f_oqrad
!    write(6,*) 'f_div   = ',f_div
!    call flush(ount)
!  END IF

!  cfun = (f_b+f_osng+f_oua+f_orad)/2.0
!
!  cfun_total=0
!  DO i=1,nvar+1
!    cfun_total=cfun_total+cfun_single(i)
!  ENDDO
!  write(*,*) 'check cost function==',cfun_total, cfun
!
!
!  PRINT*,' cfun b o===',f_b,f_osng,f_oua,f_orad,f_div
!  PRINT*,' cfun b o===',f_b,2*cfun-f_b, f_div
!
!  print *, ' cfun===',cfun, ' fb=',f_b,' f_osng=',f_osng     &
!          ,' f_oua==',f_oua,' f_orad=',f_orad
!
!  if (icall .eq. 0)  then
!      write (911, '(/,a,a,a,a)')
!    $ '  call    u+v+t+q+p             u                v',
!    $                 '                t                q',
!    $                                 '                 p'
!  endif
!
!  write (911, '(2x,i3,2x,1pe15.8,2x,a,2x,5(1pe15.8,2x))')
!    $       icall, f_b/2., 'f_b', f_bu/2., f_bv/2., f_bt/2.,
!    $                                f_bq/2., f_bp/2.
!
!  write (911, '(2x,i3,2x,1pe15.8,2x,a,2x,5(1pe15.8,2x))')
!    $       icall, f_o/2., 'f_o', f_ou/2., f_ov/2., f_ot/2.,
!    $                                f_oq/2., f_op/2.
!
!  write (911, '(2x,i3,2x,1pe15.8,2x,a,/))')
!    $       icall, cfun,'f_cost'
!
!
!  DEALLOCATE( u, v, w, wcont )

  RETURN
END SUBROUTINE costf
