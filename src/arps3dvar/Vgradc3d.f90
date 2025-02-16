!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE GRADT                     ######
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
!  Calculating the gradient of costfunction.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!
!  Jidong Gao, CAPS, July, 2000
!
!-----------------------------------------------------------------------
!
SUBROUTINE gradt(numctr,ctrv,grad,                                      &
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
           xsng,ysng,hgtsng,thesng,                                     &
           obsng,odifsng,qobsng,qualsng,isrcsng,icatsng,nobsng,         &
           xua,yua,hgtua,theua,                                         &
           obsua,odifua,qobsua,qualua,isrcua,nlevsua,nobsua,            &
           elvrad,xradc,yradc,                                          &
           distrad,uazmrad,vazmrad,hgtradc,theradc,dsdr,dhdr,           &
           obsrad,odifrad,qobsrad,qualrad,                              &
           irad,isrcrad,nlevrad,ncolrad,                                &
           xretc,yretc,hgtretc,theretc,                                 &
           obsret,odifret,qobsret,qualret,                              &
           usesng,useua, userad,                                        &
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
           sngsw,uasw,radsw,retsw,                                      &
           ipass_filt,hradius,radius_z,                                &
           div_opt,cntl_var,smth_opt,                                   &
           wgt_div_h,wgt_div_v, wgt_smth,                               &
           thermo_opt,wgt_thermo, sinlat,ffu,ffv,                       &
           anx,tem1,tem2,tem3,tem4,smcu,smcv,smcw,                      &
           u_grd,v_grd,p_grd,t_grd,q_grd,w_grd,istatus )
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
!  INCLUDE 'varpara.inc'
!
   INCLUDE 'phycst.inc'
   INCLUDE 'grid.inc'
   INCLUDE 'mp.inc'
   INCLUDE 'bndry.inc'

!-----------------------------------------------------------------------
!
!  Input Sizing Arguments
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz,sngsw,uasw,radsw,retsw
  INTEGER :: nvar,nvarradin,nvarrad
  INTEGER :: nzua,nzrdr,nzret
  INTEGER :: mxsng,mxua,mxrad,mxcolrad,mxret,mxcolret
  INTEGER :: nsrcsng,nsrcua,nsrcrad,nsrcret,ncat
  INTEGER :: mxpass,npass
  INTEGER :: ipass_filt,div_opt,cntl_var,smth_opt
  INTEGER :: thermo_opt
  REAL    :: hradius
  REAL    :: wgt_div_h, wgt_div_v, wgt_smth, wgt_thermo
  REAL    :: radius_z(nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  input grid arguments
!
!-----------------------------------------------------------------------
!
  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.
  REAL :: hterain(nx,ny)       ! The height of the terrain.

  REAL :: mapfct(nx,ny,8)      ! Map factors at scalar, u and v points

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! defined as - d( zp )/d( x ).
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! defined as - d( zp )/d( y ).
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian
                               ! defined as d( zp )/d( z ).
  REAL :: aj3x  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE X-DIR.
  REAL :: aj3y  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Y-DIR.
  REAL :: aj3z  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Z-DIR.
  REAL :: j3inv (nx,ny,nz)     ! Inverse of j3
  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3

  REAL :: rhostru(nx,ny,nz)    ! Averaged rhostr at u points (kg/m**3).
  REAL :: rhostrv(nx,ny,nz)    ! Averaged rhostr at v points (kg/m**3).
  REAL :: rhostrw(nx,ny,nz)    ! Averaged rhostr at w points (kg/m**3).

  REAL    :: xs(nx)
  REAL    :: ys(ny)
  REAL    :: zs(nx,ny,nz)
  INTEGER :: icatg(nx,ny)
  REAL    :: xcor(ncat,ncat)
!
!-----------------------------------------------------------------------
!
!  Input Observation Arguments
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=6) :: nam_var(nvar)
  REAL :: xsng(mxsng)
  REAL :: ysng(mxsng)
  REAL :: hgtsng(mxsng)
  REAL :: thesng(mxsng)
  REAL :: obsng(nvar,mxsng)
  REAL :: odifsng(nvar,mxsng)
  REAL :: qobsng(nvar,mxsng)
  INTEGER :: qualsng(nvar,mxsng)
  INTEGER :: isrcsng(mxsng)
  INTEGER :: icatsng(mxsng)
  INTEGER :: nobsng

!  REAL :: xsng_p(mxsng),ysng_p(mxsng)
  DOUBLE PRECISION :: xsng_p(mxsng),ysng_p(mxsng)
  REAL :: zsng_1(mxsng),zsng_2(mxsng)
  INTEGER :: ihgtsng(mxsng)

  REAL :: xua(mxua)
  REAL :: yua(mxua)
  REAL :: hgtua(nzua,mxua)
  REAL :: theua(nzua,mxua)
  REAL :: obsua(nvar,nzua,mxua)
  REAL :: odifua(nvar,nzua,mxua)
  REAL :: qobsua(nvar,nzua,mxua)
  INTEGER :: qualua(nvar,nzua,mxua)
  INTEGER :: nlevsua(mxua)
  INTEGER :: isrcua(mxua)
  INTEGER :: nobsua

!  REAL :: xua_p(mxua),yua_p(mxua)
  double precision :: xua_p(mxua),yua_p(mxua)
  REAL :: zua_1(nzua,mxua),zua_2(nzua,mxua)
  INTEGER :: ihgtua(nzua,mxua)

  REAL :: elvrad(mxrad)
  REAL :: xradc(mxcolrad)
  REAL :: yradc(mxcolrad)
  REAL :: distrad(mxcolrad)
  REAL :: uazmrad(mxcolrad)
  REAL :: vazmrad(mxcolrad)
  REAL :: hgtradc(nzrdr,mxcolrad)
  REAL :: theradc(nzrdr,mxcolrad)
  REAL :: dsdr(nzrdr,mxcolrad)
  REAL :: dhdr(nzrdr,mxcolrad)
  REAL :: obsrad(nvarradin,nzrdr,mxcolrad)
  REAL :: odifrad(nvarrad,nzrdr,mxcolrad)
  REAL :: qobsrad(nvarrad,nzrdr,mxcolrad)
  INTEGER :: qualrad(nvarrad,nzrdr,mxcolrad)
  INTEGER :: nlevrad(mxcolrad)
  INTEGER :: irad(mxcolrad)
  INTEGER :: isrcrad(0:mxrad)
  INTEGER :: ncolrad

!  REAL :: xradc_p(mxcolrad),yradc_p(mxcolrad)
  double precision :: xradc_p(mxcolrad),yradc_p(mxcolrad)
  REAL :: zradc_1(nzrdr,mxcolrad),zradc_2(nzrdr,mxcolrad)
  INTEGER :: ihgtradc(nzrdr,mxcolrad)

  LOGICAL, INTENT(IN) :: usesng(mxsng), useua(mxua), userad(mxcolrad)

  REAL :: xretc(mxcolret)
  REAL :: yretc(mxcolret)
  REAL :: hgtretc(nzret,mxcolret)
  REAL :: theretc(nzret,mxcolret)
  REAL :: obsret(nvar,nzret,mxcolret)
  REAL :: odifret(nvar,nzret,mxcolret)
  REAL :: qobsret(nvar,nzret,mxcolret)
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

  REAL :: xyrange(mxpass)
  REAL :: kpvrsq(nvar)
  REAL :: wlim
  REAL :: zrange(mxpass)
  REAL :: zwlim
  REAL :: thrng(mxpass)
  INTEGER :: iwstat
!
!-----------------------------------------------------------------------
!
!  Scratch Space
!
!-----------------------------------------------------------------------
!
  REAL    :: rngsqi(nvar)
  INTEGER :: knt(nvar,nz)
  REAL    :: wgtsum(nvar,nz)
  REAL    :: zsum(nvar,nz)
!
!-----------------------------------------------------------------------
!
!  Output Variables at Observation Locations
!
!-----------------------------------------------------------------------
!
  REAL :: corsng(mxsng,nvar)
  REAL :: corua(nzua,mxua,nvar)
  REAL :: corrad(nzrdr,mxcolrad,nvarrad)
  REAL :: corret(nzret,mxcolret,nvar)

  REAL :: oanxsng(nvar,mxsng)
  REAL :: oanxua(nvar,nzua,mxua)
  REAL :: oanxrad(nvarrad,nzrdr,mxcolrad)
  REAL :: oanxret(nvar,nzret,mxcolret)
!
!-----------------------------------------------------------------------
!
!  Output Grid
!
!-----------------------------------------------------------------------
!
  REAL :: anx(nx,ny,nz,nvar)
!
!-----------------------------------------------------------------------
!
!  Work arrays
!
!-----------------------------------------------------------------------
!
  REAL ::   tem1(nx,ny,nz)
  REAL ::   tem2(nx,ny,nz)
  REAL ::   tem3(nx,ny,nz)
  REAL ::   tem4(nx,ny,nz)
  REAL ::   smcu(nx,ny,nz)
  REAL ::   smcv(nx,ny,nz)
  REAL ::   smcw(nx,ny,nz)
  REAL ::    ffu(nx,ny,nz)
  REAL ::    ffv(nx,ny,nz)
  REAL :: sinlat(nx,ny,nz)
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
  REAL    :: ftabinv,setexp
  INTEGER :: i,j,k,isrc
  REAL    :: rpass,zrngsq,thrngsq
!
!----------------------------------------------------------------------
!
!    input:
!    -----
!      numctr: dimension of control variables
!      ctrv:   current perturbation (1d arrays)
!
!    output:
!    ------
!     gradient of costfunction
!
!-----------------------------------------------------------------------
!
!
  INTEGER :: numctr,iflag,num
  REAL :: ctrv (numctr), grad (numctr)
  REAL :: gdu_err(nx,ny,nz)
  REAL :: gdv_err(nx,ny,nz)
  REAL :: gdp_err(nx,ny,nz)
  REAL :: gdt_err(nx,ny,nz)
  REAL :: gdq_err(nx,ny,nz)
  REAL :: gdw_err(nx,ny,nz)
  REAL ::  gdscal(nx,ny,nz)
  REAL :: ugrid,vgrid,wgrid,vr,zv

  REAL ::   u_ctr(nx,ny,nz)
  REAL ::   v_ctr(nx,ny,nz)
  REAL ::   p_ctr(nx,ny,nz)
  REAL ::   t_ctr(nx,ny,nz)
  REAL ::   q_ctr(nx,ny,nz)
  REAL ::   w_ctr(nx,ny,nz)
  REAL ::     psi(nx,ny,nz)
  REAL ::     phi(nx,ny,nz)

  REAL, INTENT(OUT) :: u_grd(nx,ny,nz), v_grd(nx,ny,nz), w_grd(nx,ny,nz)
  REAL, INTENT(OUT) :: p_grd(nx,ny,nz), t_grd(nx,ny,nz), q_grd(nx,ny,nz)

  REAL, DIMENSION (:,:,:), ALLOCATABLE :: u,v,w,wcont
  REAL, DIMENSION (:,:,:), ALLOCATABLE :: ut,vt,wt

!  REAL,DIMENSION (:,:,:), ALLOCATABLE :: u_grd,v_grd,p_grd,t_grd,q_grd,w_grd

  REAL    :: sum1,sum2
  INTEGER :: isum

  REAL    :: f_div
  REAL    :: div3(nx,ny,nz)

  INTEGER :: ivar
  REAL    :: temp1,temp2

!  REAL, ALLOCATABLE :: outlg(:,:,:)
!  INTEGER :: nxlg, nylg
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  nxlg = (nx-3)*nproc_x+3
!  nylg = (ny-3)*nproc_y+3

  ALLOCATE ( u(nx,ny,nz), STAT = istatus )
  CALL check_alloc_status(istatus, "gradt:u")
  u = 0.0
  ALLOCATE ( v(nx,ny,nz), STAT = istatus )
  CALL check_alloc_status(istatus, "gradt:v")
  v = 0.0
  ALLOCATE ( w(nx,ny,nz), STAT = istatus )
  CALL check_alloc_status(istatus, "gradt:w")
  w = 0.0

  ALLOCATE ( wcont(nx,ny,nz), STAT = istatus )
  CALL check_alloc_status(istatus, "gradt:wcont")
  wcont = 0.0

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
!----same as cost function--define u,v,w at vector points, all other-----
!    variables are defined in scalar points
!
!-----temp arrays u_ctr, v_ctr, w_ctr are defined in vector points, these
!-----arrays are only used for div_opt== 1
!
  CALL adtrans(numctr,nx,ny,nz,u,v,p_ctr,t_ctr,q_ctr,w,ctrv,tem4)

!
! All the following arrays are on scalar grid
!
!  allocate (u_grd(nx,ny,nz), STAT = istatus)
!  CALL check_alloc_status(istatus, "gradt:u_grd")
!  allocate (v_grd(nx,ny,nz), STAT = istatus)
!  CALL check_alloc_status(istatus, "gradt:v_grd")
!  allocate (p_grd(nx,ny,nz), STAT = istatus)
!  CALL check_alloc_status(istatus, "gradt:p_grd")
!  allocate (t_grd(nx,ny,nz), STAT = istatus)
!  CALL check_alloc_status(istatus, "gradt:t_grd")
!  allocate (q_grd(nx,ny,nz), STAT = istatus)
!  CALL check_alloc_status(istatus, "gradt:q_grd")
!  allocate (w_grd(nx,ny,nz), STAT = istatus)
!  CALL check_alloc_status(istatus, "gradt:w_grd")

  DO k = 1, nz
    DO j = 1, ny
      DO i = 1, nx
        u_grd(i,j,k) = 0.
        v_grd(i,j,k) = 0.
        p_grd(i,j,k) = 0.
        t_grd(i,j,k) = 0.
        q_grd(i,j,k) = 0.
        w_grd(i,j,k) = 0.
      END DO
    END DO
  END DO
!
!  observation for single level data
!--------------------------------------------------------------
!

  IF(sngsw > 0) THEN

    num = 0

    CALL adlinearint_3d(nx,ny,nz,u_grd(1,1,1),xs(1),ys(1),zs(1,1,1),    &
           usesng, 1, mxsng, icatsng, nobsng, 1,                        &
           xsng_p,ysng_p,zsng_1,zsng_2,hgtsng,ihgtsng,corsng(1,1),tem4 )

    CALL adlinearint_3d(nx,ny,nz,v_grd(1,1,1),xs(1),ys(1),zs(1,1,1),    &
           usesng, 1, mxsng, icatsng, nobsng, 2,                        &
           xsng_p,ysng_p,zsng_1,zsng_2,hgtsng,ihgtsng,corsng(1,2),tem4 )

    CALL adlinearint_3d(nx,ny,nz,p_grd(1,1,1),xs(1),ys(1),zs(1,1,1),    &
           usesng, 1, mxsng, icatsng, nobsng, 3,                        &
           xsng_p,ysng_p,zsng_1,zsng_2,hgtsng,ihgtsng,corsng(1,3),tem4 )

    CALL adlinearint_3d(nx,ny,nz,t_grd(1,1,1),xs(1),ys(1),zs(1,1,1),    &
           usesng, 1, mxsng, icatsng, nobsng, 4,                        &
           xsng_p,ysng_p,zsng_1,zsng_2,hgtsng,ihgtsng,corsng(1,4),tem4 )

    CALL adlinearint_3d(nx,ny,nz,q_grd(1,1,1),xs(1),ys(1),zs(1,1,1),    &
           usesng, 1, mxsng, icatsng, nobsng, 5,                        &
           xsng_p,ysng_p,zsng_1,zsng_2,hgtsng,ihgtsng,corsng(1,5),tem4 )

    iflag = 0
!
!  call alinearint_3d(nx,ny,nz,w_grd(1,1,1),xs(1),ys(1),zs(1,1,1),
!    :                6,xsng(i),ysng(i),hgtsng(i),obsng(6,i),iflag )
!
!
      IF(iflag == 1) num = num+1
!
!    print*, 'num of valid data in grad cal===',num
!
  END IF
!
!  observation cost function for upper level data
!  --------------------------------------------------------------
!
  IF(uasw > 0) THEN

    CALL adlinearint_3d(nx,ny,nz,u_grd(1,1,1),xs(1),ys(1),zs(1,1,1),    &
             useua,nzua, mxua, nlevsua, nobsua, 1,                      &
             xua_p,yua_p,zua_1,zua_2,hgtua,ihgtua,corua(1,1,1), tem4 )

    CALL adlinearint_3d(nx,ny,nz,v_grd(1,1,1),xs(1),ys(1),zs(1,1,1),    &
             useua,nzua, mxua, nlevsua, nobsua, 2,                      &
             xua_p,yua_p,zua_1,zua_2,hgtua,ihgtua,corua(1,1,2), tem4 )

    CALL adlinearint_3d(nx,ny,nz,p_grd(1,1,1),xs(1),ys(1),zs(1,1,1),    &
             useua,nzua, mxua, nlevsua, nobsua, 3,                      &
             xua_p,yua_p,zua_1,zua_2,hgtua,ihgtua,corua(1,1,3), tem4 )

    CALL adlinearint_3d(nx,ny,nz,t_grd(1,1,1),xs(1),ys(1),zs(1,1,1),    &
             useua,nzua, mxua, nlevsua, nobsua, 4,                      &
             xua_p,yua_p,zua_1,zua_2,hgtua,ihgtua,corua(1,1,4), tem4 )

    CALL adlinearint_3d(nx,ny,nz,q_grd(1,1,1),xs(1),ys(1),zs(1,1,1),    &
             useua,nzua, mxua, nlevsua, nobsua, 5,                      &
             xua_p,yua_p,zua_1,zua_2,hgtua,ihgtua,corua(1,1,5), tem4 )
!
!       CALL alinearint_3d(nx,ny,nz,w_grd(1,1,1),xs(1),ys(1),zs(1,1,1), &
!          nzua, mxua, nlevsua, nobsua,                                 &
!          6,xua_p,yua_p,zua_1,zua_2,hgtua,ihgtua,corua(1,1,6) )
!
  END IF
!
!  loading radar data
! --------------------------------------------------------
!
  sum1 = 0.0
  sum2 = 0.0

  IF(radsw > 0) THEN

    corrad=0.0
    sum1= 0.
    sum2= 0.
    isum = 0
    DO i = 1, ncolrad
!   if(iuserad(isrcrad(irad(i))).GT.0) THEN
!
      IF (userad(i)) THEN
        DO k=1,nlevrad(i)

          IF(qualrad(2,k,i) > 0 .AND. ihgtradc(k,i)>=0 ) THEN
! .AND.
!           abs(odifrad(2,k,i))< 25.0 ) THEN
!           abs(odifrad(2,k,i))< 20.0 .AND. abs(obsrad(2,k,i))> 0.05 ) THEN


!         corrad(k,i,1) = uazmrad(i)*corrad(k,i,2)*dsdr(k,i)
!         corrad(k,i,3) =            corrad(k,i,2)*dhdr(k,i)
!         corrad(k,i,2) = vazmrad(i)*corrad(k,i,2)*dsdr(k,i)
!
!           sum2 = sum2 + oanxrad(2,k,i)*oanxrad(2,k,i)
            corrad(k,i,1) = uazmrad(i)*oanxrad(2,k,i)*dsdr(k,i)
            corrad(k,i,3) =            oanxrad(2,k,i)*dhdr(k,i)
            corrad(k,i,2) = vazmrad(i)*oanxrad(2,k,i)*dsdr(k,i)
            isum = isum +1

            sum1=sum1+corrad(k,i,1)*corrad(k,i,1)
            sum1=sum1+corrad(k,i,2)*corrad(k,i,2)
            sum1=sum1+corrad(k,i,3)*corrad(k,i,3)
!         ugrid = uazmrad(i)*oanxrad(2,k,i)*dsdr(k,i)
!         wgrid =            oanxrad(2,k,i)*dhdr(k,i)
!         vgrid = vazmrad(i)*oanxrad(2,k,i)*dsdr(k,i)
!         sum1=sum1+ugrid*ugrid
!         sum1=sum1+wgrid*wgrid
!         sum1=sum1+vgrid*vgrid
!         sum1=sum1+dsdr(k,i)*dsdr(k,i)+dhdr(k,i)*dhdr(k,i)
!         sum2=sum2+uazmrad(i)*uazmrad(i)+vazmrad(i)*vazmrad(i)

          END IF

        END DO

      END IF   ! userad(i)

    END DO

    CALL adlinearint_3d(nx,ny,nz,u_grd(1,1,1),xs(1),ys(1),zs(1,1,1),    &
           userad,nzrdr, mxcolrad, nlevrad, ncolrad,1,xradc_p,yradc_p,  &
           zradc_1,zradc_2,hgtradc,ihgtradc,corrad(1,1,1), tem4 )

    CALL adlinearint_3d(nx,ny,nz,v_grd(1,1,1),xs(1),ys(1),zs(1,1,1),    &
           userad,nzrdr, mxcolrad, nlevrad, ncolrad, 2,xradc_p,yradc_p, &
           zradc_1,zradc_2,hgtradc,ihgtradc,corrad(1,1,2), tem4 )

    CALL adlinearint_3d(nx,ny,nz,w_grd(1,1,1),xs(1),ys(1),zs(1,1,1),    &
           userad,nzrdr, mxcolrad, nlevrad, ncolrad, 6,xradc_p,yradc_p, &
           zradc_1,zradc_2,hgtradc,ihgtradc,corrad(1,1,3), tem4 )

    IF( div_opt== 1 ) THEN

      ! tem1 is on U grid, tem2 on V grid, tem3 on W grid
      DO k=1,nz
        DO j=1,ny
          DO i=1,nx
            tem1(i,j,k)  = 0.0
            tem2(i,j,k)  = 0.0
            tem3(i,j,k)  = 0.0
          END DO
        END DO
      END DO

      IF( (wgt_div_h > 0.0) .and. (wgt_div_v > 0.0)) THEN
        DO k=1,nz-1
          DO j=1,ny-1
            DO i=1,nx-1
              tem1(i+1,j,k)= tem1(i+1,j,k) + 1./wgt_div_h*j3inv(i,j,k)*    &
                                                mapfct(i,j,7)*div3(i,j,k)*dxinv
              tem1(i  ,j,k)= tem1(i  ,j,k) - 1./wgt_div_h*j3inv(i,j,k)*    &
                                                mapfct(i,j,7)*div3(i,j,k)*dxinv

              tem2(i,j+1,k)= tem2(i,j+1,k) + 1./wgt_div_h*j3inv(i,j,k)*    &
                                                mapfct(i,j,7)*div3(i,j,k)*dyinv
              tem2(i  ,j,k)= tem2(i  ,j,k) - 1./wgt_div_h*j3inv(i,j,k)*    &
                                                mapfct(i,j,7)*div3(i,j,k)*dyinv

              tem3(i,j,k+1)= tem3(i,j,k+1) + 1./wgt_div_v*j3inv(i,j,k)*    &
                                                div3(i,j,k)*dzinv
              tem3(i  ,j,k)= tem3(i  ,j,k) - 1./wgt_div_v*j3inv(i,j,k)*    &
                                                div3(i,j,k)*dzinv
            END DO
          END DO
        END DO
      ELSE IF( (wgt_div_h > 0.0) .AND. (wgt_div_v < 0.0)) THEN
        DO k=1,nz-1
          DO j=1,ny-1
            DO i=1,nx-1
              tem1(i+1,j,k)= tem1(i+1,j,k) + 1./wgt_div_h*j3inv(i,j,k)*   &
                                                mapfct(i,j,7)*div3(i,j,k)*dxinv
              tem1(i  ,j,k)= tem1(i  ,j,k) - 1./wgt_div_h*j3inv(i,j,k)*   &
                                                mapfct(i,j,7)*div3(i,j,k)*dxinv

              tem2(i,j+1,k)= tem2(i,j+1,k) + 1./wgt_div_h*j3inv(i,j,k)*   &
                                                mapfct(i,j,7)*div3(i,j,k)*dyinv
              tem2(i  ,j,k)= tem2(i  ,j,k) - 1./wgt_div_h*j3inv(i,j,k)*   &
                                                mapfct(i,j,7)*div3(i,j,k)*dyinv
            END DO
          END DO
        END DO
      END IF

      IF (mp_opt > 0) THEN
        !CALL acct_interrupt(mp_acct)
        CALL mpsendrecv2dew(tem1, nx, ny, nz, ebc, wbc, 1, tem4)
        CALL mpsendrecv2dns(tem1, nx, ny, nz, nbc, sbc, 1, tem4)

        CALL mpsendrecv2dew(tem2, nx, ny, nz, ebc, wbc, 2, tem4)
        CALL mpsendrecv2dns(tem2, nx, ny, nz, nbc, sbc, 2, tem4)
        !CALL acct_stop_inter
      END IF

      !
      ! u_ctr, v_ctr, wcont are vectors here
      !
      IF(wgt_div_v > 0.0 ) THEN
        DO k=1,nz
          DO j=1,ny-1
            DO i=1,nx-1
              wcont(i,j,k)=tem3(i,j,k)*rhostrw(i,j,k)
              tem3(i,j,k)=0.0
            END DO
          END DO
        END DO
      ENDIF

      IF(wgt_div_h > 0.0 ) THEN
        DO k=1,nz-1
          DO j=1,ny
            DO i=1,nx-1
              v_ctr(i,j,k)=tem2(i,j,k)*rhostrv(i,j,k)*mapfct(i,j,6)
              tem2(i,j,k)=0.0
            END DO
          END DO
        END DO

        DO k=1,nz-1
          DO j=1,ny-1
            DO i=1,nx
              u_ctr(i,j,k)=tem1(i,j,k)*rhostru(i,j,k)*mapfct(i,j,5)
              tem1(i,j,k) = 0.0
            END DO
          END DO
        END DO
      ENDIF

      ! No. 2 message passing inside which is commented out.
      CALL adwcontra(nx,ny,nz,u_ctr,v_ctr,w_ctr,mapfct,j1,j2,j3,aj3z,   &
                     rhostr,rhostru,rhostrv,rhostrw,wcont,tem1,tem2,tem4)

      !zhaokun add the gradient of mass contiunity into the total grd. Note the
      !u_ctr,v_ctr_w_ctr is at vector point, and should be converting to mass points
      DO k = 1, nz-1
        DO j = 1, ny-1
          DO i = 1, nx-1
            u_grd(i,j,k) = u_grd(i,j,k) + 0.5*(u_ctr(i,j,k)+u_ctr(i+1,j,k))
            v_grd(i,j,k) = v_grd(i,j,k) + 0.5*(v_ctr(i,j,k)+v_ctr(i,j+1,k))
            w_grd(i,j,k) = w_grd(i,j,k) + 0.5*(w_ctr(i,j,k)+w_ctr(i,j,k+1))
            u_ctr(i,j,k) = 0.0
            v_ctr(i,j,k) = 0.0
            w_ctr(i,j,k) = 0.0
          END DO
        END DO
      END DO

    END IF  !end of div_opt

  END IF

  IF( thermo_opt == 1 ) THEN

    ALLOCATE ( ut(nx,ny,nz), STAT = istatus )
    CALL check_alloc_status( istatus, "gradt:ut" )
    ut = 0.0
    ALLOCATE ( vt(nx,ny,nz), STAT = istatus )
    CALL check_alloc_status( istatus, "gradt:vt" )
    vt = 0.0
    ALLOCATE ( wt(nx,ny,nz), STAT = istatus )
    CALL check_alloc_status( istatus, "gradt:wt" )
    wt = 0.0

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          ut(i,j,k) = anx(i,j,k,1)+u_ctr(i,j,k)
          vt(i,j,k) = anx(i,j,k,2)+v_ctr(i,j,k)
          ! w used for temperature here temporarily
          wt(i,j,k) = (  anx(i,j,k,4)+t_ctr(i,j,k) )                    &
                     *(( anx(i,j,k,3)+p_ctr(i,j,k) )/p0 )**rddcp
          wcont(i,j,k) = 0.0
          ffu(i,j,k)=wgt_thermo*ffu(i,j,k)
          ffv(i,j,k)=wgt_thermo*ffv(i,j,k)
        END DO
      END DO
    END DO

    tem1(:,:,:) = 0.0
    tem2(:,:,:) = 0.0
    tem3(:,:,:) = 0.0

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem1(i,j,k+1)=tem1(i,j,k+1) +ffu(i,j,k)*dzinv*j3inv(i,j,k)
          tem1(i,j,k  )=tem1(i,j,k  ) -ffu(i,j,k)*dzinv*j3inv(i,j,k)

          tem3(i,j,k  )=tem3(i,j,k  ) -g/sinlat(i,j,k)*ffu(i,j,k)        &
                            /wt(i,j,k)**2*(wt(i,j+1,k)-wt(i,j,k))*dyinv
          tem3(i,j+1,k)=tem3(i,j+1,k) +g/sinlat(i,j,k)/wt(i,j,k)         &
                                      *ffu(i,j,k)*dyinv
          tem3(i,j  ,k)=tem3(i,j  ,k) -g/sinlat(i,j,k)/wt(i,j,k)         &
                                      *ffu(i,j,k)*dyinv

          tem1(i,j,k  )=tem1(i,j,k  ) -ffu(i,j,k)/wt(i,j,k)              &
                            *(wt(i,j,k+1)-wt(i,j,k))*dzinv*j3inv(i,j,k)
          tem3(i,j,k  )=tem3(i,j,k  ) +ut(i,j,k)*ffu(i,j,k)/wt(i,j,k)**2 &
                            *(wt(i,j,k+1)-wt(i,j,k))*dzinv*j3inv(i,j,k)
          tem3(i,j,k+1)=tem3(i,j,k+1) -ut(i,j,k)/wt(i,j,k)*ffu(i,j,k)    &
                            *dzinv*j3inv(i,j,k)
          tem3(i,j,k  )=tem3(i,j,k  ) +ut(i,j,k)/wt(i,j,k)*ffu(i,j,k)    &
                            *dzinv*j3inv(i,j,k)
          ffu(i,j,k) = 0.0
!
          tem2(i,j,k+1)=tem2(i,j,k+1) +ffv(i,j,k)*dzinv*j3inv(i,j,k)
          tem2(i,j,k  )=tem2(i,j,k  ) -ffv(i,j,k)*dzinv*j3inv(i,j,k)
!
          tem3(i,j,k  )=tem3(i,j,k  ) +g/sinlat(i,j,k)*ffv(i,j,k)        &
                            /wt(i,j,k)**2*(wt(i+1,j,k)-wt(i,j,k))*dxinv
          tem3(i+1,j,k)=tem3(i+1,j,k) -g/sinlat(i,j,k)/wt(i,j,k)         &
                                      *ffv(i,j,k)*dxinv
          tem3(i,j  ,k)=tem3(i,j  ,k) +g/sinlat(i,j,k)/w(i,j,k)          &
                                      *ffv(i,j,k)*dxinv
!
          tem2(i,j,k  )=tem2(i,j,k  ) -ffv(i,j,k)/wt(i,j,k)              &
                            *(wt(i,j,k+1)-wt(i,j,k))*dzinv*j3inv(i,j,k)
          tem3(i,j,k  )=tem3(i,j,k  ) +vt(i,j,k)*ffv(i,j,k)/wt(i,j,k)**2 &
                            *(wt(i,j,k+1)-wt(i,j,k))*dzinv*j3inv(i,j,k)
          tem3(i,j,k+1)=tem3(i,j,k+1) -vt(i,j,k)/wt(i,j,k)*ffv(i,j,k)    &
                            *dzinv*j3inv(i,j,k)
          tem3(i,j,k  )=tem3(i,j,k  ) +vt(i,j,k)/wt(i,j,k)*ffv(i,j,k)    &
                            *dzinv*j3inv(i,j,k)
          ffv(i,j,k) = 0.0
        END DO
      END DO
    END DO

    IF (mp_opt > 0) THEN
      !CALL acct_interrupt(mp_acct)
      CALL mpsendrecv2dew(tem1, nx, ny, nz, ebc, wbc, 0, tem4)
      CALL mpsendrecv2dns(tem1, nx, ny, nz, nbc, sbc, 0, tem4)

      CALL mpsendrecv2dew(tem2, nx, ny, nz, ebc, wbc, 0, tem4)
      CALL mpsendrecv2dns(tem2, nx, ny, nz, nbc, sbc, 0, tem4)

      CALL mpsendrecv2dew(tem3, nx, ny, nz, ebc, wbc, 0, tem4)
      CALL mpsendrecv2dns(tem3, nx, ny, nz, nbc, sbc, 0, tem4)
      !CALL acct_stop_inter
    END IF
!
    DO k= 1,nz-1
      DO j= 1,ny-1
        DO i= 1,nx-1
          u_grd(i,j,k) = u_grd(i,j,k) + tem1(i,j,k)
          v_grd(i,j,k) = v_grd(i,j,k) + tem2(i,j,k)
          t_grd(i,j,k) = t_grd(i,j,k) + tem3(i,j,k)                      &
                  *(( anx(i,j,k,3)+p_ctr(i,j,k) )/p0 )**rddcp
          p_grd(i,j,k) = p_grd(i,j,k) + tem3(i,j,k)                      &
                  *rddcp*( anx(i,j,k,4)+t_ctr(i,j,k) )/p0                &
                  *(( anx(i,j,k,3)+p_ctr(i,j,k) )/p0 )**(rddcp-1)
          !tem1(i,j,k) = 0.0
          !tem2(i,j,k) = 0.0
          !tem3(i,j,k) = 0.0
        END DO
      END DO
    END DO

    DEALLOCATE( ut, vt, wt )
  END IF

  IF( smth_opt== 1 ) THEN

    DO  k=1,nz
      DO  j=1,ny
        DO  i=1,nx
          wcont(i,j,k)=0.0   !!!! used as temporary array !!!!
        END DO
      END DO
    END DO

    DO i=nx-1,1,-1
      DO j=ny-1,1,-1
        DO k=nz-1,1,-1
          wcont(i,j,k)=smcu(i,j,k)*wgt_smth
        END DO
      END DO
    END DO

    DO i=nx-1,2,-1
      DO j=ny-1,2,-1
        u_grd(i,j,1)   =u_grd(i,j,1)-wcont(i,j,1)*4.0
        u_grd(i,j+1,1) =u_grd(i,j+1,1) +wcont(i,j,1)
        u_grd(i,j-1,1) =u_grd(i,j-1,1) +wcont(i,j,1)
        u_grd(i+1,j,1) =u_grd(i+1,j,1) +wcont(i,j,1)
        u_grd(i-1,j,1) =u_grd(i-1,j,1) +wcont(i,j,1)
      END DO
    END DO
!
    DO i=nx-1,2,-1
      DO j=ny-1,2,-1
        DO k=nz-1,2,-1
          u_grd(i,j,k)=u_grd(i,j,k) -wcont(i,j,k)*6.0
          u_grd(i,j,k+1)=u_grd(i,j,k+1) +wcont(i,j,k)
          u_grd(i,j,k-1)=u_grd(i,j,k-1) +wcont(i,j,k)
          u_grd(i,j+1,k)=u_grd(i,j+1,k) +wcont(i,j,k)
          u_grd(i,j-1,k)=u_grd(i,j-1,k) +wcont(i,j,k)
          u_grd(i+1,j,k)=u_grd(i+1,j,k) +wcont(i,j,k)
          u_grd(i-1,j,k)=u_grd(i-1,j,k) +wcont(i,j,k)
        END DO
      END DO
    END DO

    IF (loc_x > 1) THEN  ! count column 1's affect on column 2
      DO j=ny-1,2,-1
        DO k=nz-1,1,-1
          u_grd(2,j,k)=u_grd(2,j,k) +wcont(1,j,k)
        END DO
      END DO
    END IF

    IF (loc_y > 1) THEN  ! Count row 1's affect on column 2
      DO i=nx-1,2,-1
        DO k=nz-1,1,-1
          u_grd(i,2,k)=u_grd(i,2,k) +wcont(i,1,k)
        END DO
      END DO
    END IF

    DO i=nx-1,2,-1
      DO j=ny-1,2,-1
        DO k=nz-1,1,-1
          wcont(i,j,k)=smcv(i,j,k)*wgt_smth
        END DO
      END DO
    END DO

    DO i=nx-1,2,-1
      DO j=ny-1,2,-1
        v_grd(i,j,1)=v_grd(i,j,1)-wcont(i,j,1)*4.0
        v_grd(i,j+1,1)=v_grd(i,j+1,1)+wcont(i,j,1)
        v_grd(i,j-1,1)=v_grd(i,j-1,1)+wcont(i,j,1)
        v_grd(i+1,j,1)=v_grd(i+1,j,1)+wcont(i,j,1)
        v_grd(i-1,j,1)=v_grd(i-1,j,1)+wcont(i,j,1)
      END DO
    END DO
!
    DO i=nx-1,2,-1
      DO j=ny-1,2,-1
        DO k=nz-1,2,-1
          v_grd(i,j,k)=v_grd(i,j,k)-wcont(i,j,k)*6.0
          v_grd(i,j,k+1)=v_grd(i,j,k+1)+wcont(i,j,k)
          v_grd(i,j,k-1)=v_grd(i,j,k-1)+wcont(i,j,k)
          v_grd(i,j+1,k)=v_grd(i,j+1,k)+wcont(i,j,k)
          v_grd(i,j-1,k)=v_grd(i,j-1,k)+wcont(i,j,k)
          v_grd(i+1,j,k)=v_grd(i+1,j,k)+wcont(i,j,k)
          v_grd(i-1,j,k)=v_grd(i-1,j,k)+wcont(i,j,k)
        END DO
      END DO
    END DO

    IF (loc_x > 1) THEN
      DO j=ny-1,2,-1
        DO k=nz-1,1,-1
          v_grd(2,j,k)=v_grd(2,j,k) +wcont(1,j,k)
        END DO
      END DO
    END IF

    IF (loc_y > 1) THEN
      DO i=nx-1,2,-1
        DO k=nz-1,1,-1
          v_grd(i,2,k)=v_grd(i,2,k) +wcont(i,1,k)
        END DO
      END DO
    END IF

    DO i=nx-1,2,-1
      DO j=ny-1,2,-1
        DO k=nz-1,1,-1
          wcont(i,j,k)=smcw(i,j,k)*wgt_smth
        END DO
      END DO
    END DO

    DO i=nx-1,2,-1
      DO j=ny-1,2,-1
        w_grd(i,j,1)=w_grd(i,j,1)-wcont(i,j,1)*4.0
        w_grd(i,j+1,1)=w_grd(i,j+1,1)+wcont(i,j,1)
        w_grd(i,j-1,1)=w_grd(i,j-1,1)+wcont(i,j,1)
        w_grd(i+1,j,1)=w_grd(i+1,j,1)+wcont(i,j,1)
        w_grd(i-1,j,1)=w_grd(i-1,j,1)+wcont(i,j,1)
      END DO
    END DO

    DO i=nx-1,2,-1
      DO j=ny-1,2,-1
        DO k=nz-1,2,-1
          w_grd(i,j,k)=w_grd(i,j,k)-wcont(i,j,k)*6.0
          w_grd(i,j,k+1)=w_grd(i,j,k+1)+wcont(i,j,k)
          w_grd(i,j,k-1)=w_grd(i,j,k-1)+wcont(i,j,k)
          w_grd(i,j+1,k)=w_grd(i,j+1,k)+wcont(i,j,k)
          w_grd(i,j-1,k)=w_grd(i,j-1,k)+wcont(i,j,k)
          w_grd(i+1,j,k)=w_grd(i+1,j,k)+wcont(i,j,k)
          w_grd(i-1,j,k)=w_grd(i-1,j,k)+wcont(i,j,k)
        END DO
      END DO
    END DO

    IF (loc_x > 1) THEN  ! count column 1's affect on column 2
      DO j=ny-1,2,-1
        DO k=nz-1,1,-1
          w_grd(2,j,k)=w_grd(2,j,k) +wcont(1,j,k)
        END DO
      END DO
    END IF

    IF (loc_y > 1) THEN  ! Count row 1's affect on column 2
      DO i=nx-1,2,-1
        DO k=nz-1,1,-1
          w_grd(i,2,k)=w_grd(i,2,k) +wcont(i,1,k)
        END DO
      END DO
    END IF

    IF (mp_opt > 0) THEN
      !CALL acct_interrupt(mp_acct)
      CALL mpsendrecv2dew(u_grd, nx, ny, nz, ebc, wbc, 0, tem4)
      CALL mpsendrecv2dns(u_grd, nx, ny, nz, nbc, sbc, 0, tem4)

      CALL mpsendrecv2dew(v_grd, nx, ny, nz, ebc, wbc, 0, tem4)
      CALL mpsendrecv2dns(v_grd, nx, ny, nz, nbc, sbc, 0, tem4)

      CALL mpsendrecv2dew(w_grd, nx, ny, nz, ebc, wbc, 0, tem4)
      CALL mpsendrecv2dns(w_grd, nx, ny, nz, nbc, sbc, 0, tem4)
      !CALL acct_stop_inter
    END IF
  END IF
!
!
!--IF cntl_var == 0 u & v are control varibles---------------------
!
!
  IF( cntl_var == 0 ) THEN

    DO k = 1, nz-1
      DO j = 1, ny-1
        DO i = 1, nx-1
          psi(i,j,k)= u_grd(i,j,k)
          phi(i,j,k)= v_grd(i,j,k)
        END DO
      END DO
    END DO

  ELSE    ! Should be checked for MPI

    DO k = 1, nz
      IF (loc_x == nproc_x .AND. loc_y == nproc_y) THEN
        v_grd(nx-1,ny,k)=v_grd(nx-1,ny,k)+0.5*v_grd(nx,ny,k)
        v_grd(nx,ny-1,k)=v_grd(nx,ny-1,k)+0.5*v_grd(nx,ny,k)
        v_grd(nx,ny,k)=0.0
        u_grd(nx-1,ny,k)=u_grd(nx-1,ny,k)+0.5*u_grd(nx,ny,k)
        u_grd(nx,ny-1,k)=u_grd(nx,ny-1,k)+0.5*u_grd(nx,ny,k)
        u_grd(nx,ny,k)=0.0
      END IF

      IF (loc_x == nproc_x .AND. loc_y == 1) THEN
        v_grd(nx,2,k)=v_grd(nx,2,k)+0.5*v_grd(nx,1,k)
        v_grd(nx-1,1,k)=v_grd(nx-1,1,k)+0.5*v_grd(nx,1,k)
        v_grd(nx,1,k)=0.0
        u_grd(nx,2,k)=u_grd(nx,2,k)+0.5*u_grd(nx,1,k)
        u_grd(nx-1,1,k)=u_grd(nx-1,1,k)+0.5*u_grd(nx,1,k)
        u_grd(nx,1,k)=0.0
      END IF

      IF (loc_x == 1 .AND. loc_y == nproc_y) THEN
        v_grd(2,ny,k)=v_grd(2,ny,k)+0.5*v_grd(1,ny,k)
        v_grd(1,ny-1,k)=v_grd(1,ny-1,k)+0.5*v_grd(1,ny,k)
        v_grd(1,ny,k)=0.0
        u_grd(2,ny,k)=u_grd(2,ny,k)+0.5*u_grd(1,ny,k)
        u_grd(1,ny-1,k)=u_grd(1,ny-1,k)+0.5*u_grd(1,ny,k)
        u_grd(1,ny,k)=0.0
      END IF

      IF (loc_x == 1 .AND. loc_y == 1) THEN
        v_grd(1,2,k)=v_grd(1,2,k)+0.5*v_grd(1,1,k)
        v_grd(2,1,k)=v_grd(2,1,k)+0.5*v_grd(1,1,k)
        v_grd(1,1,k)=0.0
        u_grd(1,2,k)=u_grd(1,2,k)+0.5*u_grd(1,1,k)
        u_grd(2,1,k)=u_grd(2,1,k)+0.5*u_grd(1,1,k)
        u_grd(1,1,k)=0.0
      END IF

      IF (loc_y == nproc_y) THEN
        DO i=2,nx-1
          v_grd(i,ny-2,k)=v_grd(i,ny-2,k)-v_grd(i,ny,k)
          v_grd(i,ny-1,k)=v_grd(i,ny-1,k)+v_grd(i,ny,k)
          v_grd(i,ny-1,k)=v_grd(i,ny-1,k)+v_grd(i,ny,k)
          v_grd(i,ny,k)=0.0
          u_grd(i,ny-2,k)=u_grd(i,ny-2,k)-u_grd(i,ny,k)
          u_grd(i,ny-1,k)=u_grd(i,ny-1,k)+u_grd(i,ny,k)
          u_grd(i,ny-1,k)=u_grd(i,ny-1,k)+u_grd(i,ny,k)
          u_grd(i,ny,k)=0.0
        END DO
      END IF

      IF (loc_y == 1) THEN
        DO i=2,nx-1
          v_grd(i,   3,k)=v_grd(i,   3,k)-v_grd(i, 1,k)
          v_grd(i,   2,k)=v_grd(i,   2,k)+v_grd(i, 1,k)
          v_grd(i,   2,k)=v_grd(i,   2,k)+v_grd(i, 1,k)
          v_grd(i, 1,k)=0.0
          u_grd(i,   3,k)=u_grd(i,   3,k)-u_grd(i, 1,k)
          u_grd(i,   2,k)=u_grd(i,   2,k)+u_grd(i, 1,k)
          u_grd(i,   2,k)=u_grd(i,   2,k)+u_grd(i, 1,k)
          u_grd(i, 1,k)=0.0
        END DO
      END IF

      IF (loc_x == nproc_x) THEN
        DO j=2,ny-1
          v_grd(nx-2,j,k)=v_grd(nx-2,j,k)-v_grd(nx,j,k)
          v_grd(nx-1,j,k)=v_grd(nx-1,j,k)+v_grd(nx,j,k)
          v_grd(nx-1,j,k)=v_grd(nx-1,j,k)+v_grd(nx,j,k)
          v_grd(nx,j,k)=0.0
          u_grd(nx-2,j,k)=u_grd(nx-2,j,k)-u_grd(nx,j,k)
          u_grd(nx-1,j,k)=u_grd(nx-1,j,k)+u_grd(nx,j,k)
          u_grd(nx-1,j,k)=u_grd(nx-1,j,k)+u_grd(nx,j,k)
          u_grd(nx,j,k)=0.0
        END DO
      END IF

      IF (loc_x == 1) THEN
        DO j=2,ny-1
          v_grd(3,j,k)=v_grd(   3,j,k)-v_grd( 1,j,k)
          v_grd(2,j,k)=v_grd(   2,j,k)+v_grd( 1,j,k)
          v_grd(2,j,k)=v_grd(   2,j,k)+v_grd( 1,j,k)
          v_grd(1,j,k)=0.0
          u_grd(3,j,k)=u_grd(   3,j,k)-u_grd( 1,j,k)
          u_grd(2,j,k)=u_grd(   2,j,k)+u_grd( 1,j,k)
          u_grd(2,j,k)=u_grd(   2,j,k)+u_grd( 1,j,k)
          u_grd(1,j,k)=0.0
        END DO
      END IF

      DO j = 2, ny-1
        DO i = 2, nx-1

          u_grd(i,j,k)   = u_grd(i,j,k)*mapfct(i,j,2)
          psi(i-1,j+1,k) = psi(i-1,j+1,k) + u_grd(i,j,k)/dy/4.
          psi(i,  j+1,k) = psi(i,  j+1,k) + u_grd(i,j,k)/dy/4.
          psi(i-1,j-1,k) = psi(i-1,j-1,k) - u_grd(i,j,k)/dy/4.
          psi(i,  j-1,k) = psi(i,  j-1,k) - u_grd(i,j,k)/dy/4.

          phi(i,  j,k)  = phi(i,  j,k) + u_grd(i,j,k)/dx
          phi(i-1,j,k)  = phi(i-1,j,k) - u_grd(i,j,k)/dx
          u_grd(i,j,k)  = 0.0   ! unnecessary - WYH

          v_grd(i,j,k)   = v_grd(i,j,k)*mapfct(i,j,3)
          psi(i+1,j-1,k) = psi(i+1,j-1,k) + v_grd(i,j,k)/dx/4.
          psi(i+1,  j,k) = psi(i+1,  j,k) + v_grd(i,j,k)/dx/4.
          psi(i-1,j-1,k) = psi(i-1,j-1,k) - v_grd(i,j,k)/dx/4.
          psi(i-1,j,  k) = psi(i-1,  j,k) - v_grd(i,j,k)/dx/4.

          phi(i,  j,k)  = phi(i,  j,k) + v_grd(i,j,k)/dy
          phi(i,j-1,k)  = phi(i,j-1,k) - v_grd(i,j,k)/dy
          v_grd(i,j,k)  = 0.0   ! unnecessary - WYH

        END DO
      END DO

      IF (loc_x > 1) THEN ! count 1st column's effect on 2nd column
        DO j = 2, ny-1
          v_grd(1,j,k) = v_grd(1,j,k)*mapfct(1,j,3)
          psi(2,j-1,k) = psi(2,j-1,k) + v_grd(1,j,k)/dx/4.
          psi(2,  j,k) = psi(2,  j,k) + v_grd(1,j,k)/dx/4.
          v_grd(1,j,k) = 0.0
        END DO
      END IF

      IF (loc_y > 1) THEN ! count 1st row's effect on 2nd row
        DO i = 2, nx-1
          u_grd(i,1,k) = u_grd(i,1,k)*mapfct(i,1,2)
          psi(i-1,2,k) = psi(i-1,2,k) + u_grd(i,1,k)/dy/4.
          psi(i,  2,k) = psi(i,  2,k) + u_grd(i,1,k)/dy/4.
          u_grd(i,1,k) = 0.0
        END DO
      END IF
    END DO

    !
    ! Msg for u_grd, v_grd, psi, phi
    !
    IF (mp_opt > 0) THEN   ! Maybe not exactly the same, should be checked later
!      CALL acct_interrupt(mp_acct)
      CALL mpsendrecv2dew(psi, nx, ny, nz, ebc, wbc, 0, tem4)
      CALL mpsendrecv2dns(psi, nx, ny, nz, nbc, sbc, 0, tem4)

      CALL mpsendrecv2dew(phi, nx, ny, nz, ebc, wbc, 0, tem4)
      CALL mpsendrecv2dns(phi, nx, ny, nz, nbc, sbc, 0, tem4)
!      CALL acct_stop_inter
    END IF

  END IF

  tem1=0.0;tem2=0.0;tem3=0.0;tem4=0.0
  CALL  vbl_to_ctr(ipass_filt,hradius,radius_z,nx,ny,nz,gdu_err,       &
                   gdscal,  psi, tem1,tem2,tem3,tem4)

  tem1=0.0;tem2=0.0;tem3=0.0;tem4=0.0
  CALL  vbl_to_ctr(ipass_filt,hradius,radius_z,nx,ny,nz,gdv_err,       &
                   gdscal,  phi, tem1,tem2,tem3,tem4)

  tem1=0.0;tem2=0.0;tem3=0.0;tem4=0.0
  CALL  vbl_to_ctr(ipass_filt,hradius,radius_z,nx,ny,nz,gdp_err,       &
                   gdscal,p_grd, tem1,tem2,tem3,tem4)

  tem1=0.0;tem2=0.0;tem3=0.0;tem4=0.0
  CALL  vbl_to_ctr(ipass_filt,hradius,radius_z,nx,ny,nz,gdt_err,       &
                   gdscal,t_grd, tem1,tem2,tem3,tem4)

  tem1=0.0;tem2=0.0;tem3=0.0;tem4=0.0
  CALL  vbl_to_ctr(ipass_filt,hradius,radius_z,nx,ny,nz,gdq_err,       &
                   gdscal,q_grd, tem1,tem2,tem3,tem4)

  tem1=0.0;tem2=0.0;tem3=0.0;tem4=0.0
  CALL  vbl_to_ctr(ipass_filt,hradius,radius_z,nx,ny,nz,gdw_err,       &
                   gdscal,w_grd, tem1,tem2,tem3,tem4)

!
!------------return gradients from scalars to vector points--------------
!
  !w_ctr = 0.0  ! w_ctr will be actually w_grd below
  !DO k=1, nz-1
  !  DO j=1, ny-1
  !    DO i=1, nx-1
  !      w_ctr(i,j,k)   = w_ctr(i,j,k)  + 0.5*w_grd(i,j,k)
  !      w_ctr(i,j,k+1) = w_ctr(i,j,k+1)+ 0.5*w_grd(i,j,k)
  !    END DO
  !  END DO
  !END DO

  DO k = 1, nz-1
    DO j = 1, ny-1
      DO i = 1, nx-1
        u(i,j,k) = psi(i,j,k) + u(i,j,k)
      END DO
    END DO
  END DO


  DO k = 1, nz-1
    DO j = 1, ny-1
      DO i = 1, nx-1
        v(i,j,k) = phi(i,j,k) + v(i,j,k)
      END DO
    END DO
  END DO

  DO k = 1, nz
    DO j = 1, ny-1
      DO i = 1, nx-1
        !w(i,j,k) = w_ctr(i,j,k) + w(i,j,k)
        w(i,j,k) = w_grd(i,j,k) + w(i,j,k)
      END DO
    END DO
  END DO

  DO k = 1, nz-1
    DO j = 1, ny-1
      DO i = 1, nx-1
        p_grd(i,j,k) = p_grd(i,j,k) + p_ctr(i,j,k)
        t_grd(i,j,k) = t_grd(i,j,k) + t_ctr(i,j,k)
        q_grd(i,j,k) = q_grd(i,j,k) + q_ctr(i,j,k)
      END DO
    END DO
  END DO

  CALL trans(numctr,nx,ny,nz,u,v,p_grd,t_grd,q_grd,w,grad)

  DEALLOCATE( u, v, w, wcont )
!  DEALLOCATE( u_grd, v_grd, p_grd, t_grd, q_grd, w_grd )

  RETURN
END SUBROUTINE gradt

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ADWCONTRA                  ######
!######                                                      ######
!######                Copyright (c) 1994                    ######
!######  Center for the Analysis and Prediction of Storms    ######
!######  University of Oklahoma.  All rights reserved.       ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE ADWCONTRA(nx,ny,nz,u,v,w,mapfct,j1,j2,j3,aj3z,               &
                     rhostr,rhostru,rhostrv,rhostrw,wcont,ustr,vstr,tem1)
!
!-----------------------------------------------------------------------
!
! PURPOSE:
!
!   Perform the adjoint operations on WCONTRA. WCONTRA
!   calculates wcont, the contravariant vertical velocity (m/s)
!
!-----------------------------------------------------------------------
!
! AUTHOR:   Jidong Gao
! 5/20/96 converted to ARPS4.0.22
!
!-----------------------------------------------------------------------
!
! INPUT:
!
!       nx       Number of grid points in the x-direction (east/west)
!       ny       Number of grid points in the y-direction (north/south)
!       nz       Number of grid points in the vertical
!
!       u        x component of velocity at all time levels (m/s)
!       v        y component of velocity at all time levels (m/s)
!       w        Vertical component of Cartesian velocity
!                at all time levels (m/s)
!
!    mapfct   Map factors at scalar, u and v points
!
!       j1       Coordinate transform Jacobian -d(zp)/dx
!       j2       Coordinate transform Jacobian -d(zp)/dy
!       j3       Coordinate transform Jacobian  d(zp)/dz
!    aj3z     Avgz of the coordinate transformation Jacobian  d(zp)/dz
!
!       rhostr   j3 times base state density rhobar(kg/m**3).
!       rhostru  Average rhostr at u points (kg/m**3).
!       rhostrv  Average rhostr at v points (kg/m**3).
!       rhostrw  Average rhostr at w points (kg/m**3).
!
! OUTPUT:
!
!       wcont    Vertical component of contravariant velocity in
!                computational coordinates (m/s)
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
! Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER, INTENT(IN) :: nx,ny,nz          ! The number of grid points in 3
                                           ! directions

  REAL, INTENT(INOUT) :: u     (nx,ny,nz)  ! Total u-velocity (m/s)
  REAL, INTENT(INOUT) :: v     (nx,ny,nz)  ! Total v-velocity (m/s)
  REAL, INTENT(INOUT) :: w     (nx,ny,nz)  ! Total w-velocity (m/s)

  REAL,    INTENT(IN) :: mapfct(nx,ny,8)   ! Map factors at scalar, u and v points

  REAL,    INTENT(IN) :: j1    (nx,ny,nz)  ! Coordinate transform Jacobian
                                           ! defined as - d( zp )/d( x ).
  REAL,    INTENT(IN) :: j2    (nx,ny,nz)  ! Coordinate transform Jacobian
                                           ! defined as - d( zp )/d( y ).
  REAL,    INTENT(IN) :: j3    (nx,ny,nz)  ! Coordinate transform Jacobian
                                           ! defined as d( zp )/d( z ).
  REAL,    INTENT(IN) :: aj3z  (nx,ny,nz)  ! Coordinate transformation Jacobian defined
                                           ! as d( zp )/d( z ) AVERAGED IN THE Z-DIR.

  REAL,    INTENT(IN) :: rhostr(nx,ny,nz)  ! j3 times base state density rhobar
                                           ! (kg/m**3).
  REAL,    INTENT(IN) :: rhostru(nx,ny,nz) ! Average rhostr at u points (kg/m**3).
  REAL,    INTENT(IN) :: rhostrv(nx,ny,nz) ! Average rhostr at v points (kg/m**3).
  REAL,    INTENT(IN) :: rhostrw(nx,ny,nz) ! Average rhostr at w points (kg/m**3).

  REAL, INTENT(OUT)   :: wcont (nx,ny,nz)  ! Vertical velocity in computational
                                           ! coordinates (m/s)

  REAL, INTENT(INOUT) :: ustr  (nx,ny,nz)  ! temporary work array
  REAL, INTENT(INOUT) :: vstr  (nx,ny,nz)  ! temporary work array

  REAL, INTENT(INOUT) :: tem1  (nx,ny,nz)  ! temporary work array
!
!-----------------------------------------------------------------------
!
! Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'bndry.inc'
  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
!  INTEGER :: onvf
!  REAL    :: tem2  (nx,ny,nz)     ! temporary work array
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!     Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  CALL advbcwcont(nx,ny,nz,wcont)

  IF( crdtrns == 0 ) THEN  ! No coord. transformation case.

    DO k= 2,nz-1
      DO j= 1,ny-1
        DO i= 1,nx-1
          w(i,j,k)=w(i,j,k) + wcont(i,j,k)
          wcont(i,j,k)=0.0
        END DO
      END DO
    END DO

  ELSEIF( ternopt == 0) THEN

    DO k= 2,nz-1
      DO j= 1,ny-1
        DO i= 1,nx-1
          w(i,j,k)=w(i,j,k) + wcont(i,j,k)/aj3z(i,j,k)
          wcont(i,j,k)=0.0
        END DO
      END DO
    END DO

  ELSE

    ! ustr is on U grid, vstr on V grid below
    DO k= 2,nz-1
      DO j= 1,ny-1
        DO i= 1,nx-1
!         wcont(i,j,k)= (                                               &
!             ((ustr(i  ,j,k)+ustr(i  ,j,k-1))*j1(i  ,j,k)              &
!             +(ustr(i+1,j,k)+ustr(i+1,j,k-1))*j1(i+1,j,k)              &
!             +(vstr(i  ,j,k)+vstr(i  ,j,k-1))*j2(i  ,j,k)              &
!             +(vstr(i,j+1,k)+vstr(i,j+1,k-1))*j2(i,j+1,k))             &
!             * mapfct(i,j,8)                                           &
!             / rhostrw(i,j,k) + w(i,j,k)                               &
!             ) /aj3z(i,j,k)
!
         ustr(i,j,k  )=ustr(i,j,k  ) + wcont(i,j,k)*j1(i,j,k)           &
            *mapfct(i,j,8)/rhostrw(i,j,k)/aj3z(i,j,k)
         ustr(i,j,k-1)=ustr(i,j,k-1) + wcont(i,j,k)*j1(i,j,k)           &
            *mapfct(i,j,8)/rhostrw(i,j,k)/aj3z(i,j,k)
         ustr(i+1,j,k  )=ustr(i+1,j,k  ) + wcont(i,j,k)*j1(i+1,j,k)     &
            *mapfct(i,j,8)/rhostrw(i,j,k)/aj3z(i,j,k)
         ustr(i+1,j,k-1)=ustr(i+1,j,k-1) + wcont(i,j,k)*j1(i+1,j,k)     &
            *mapfct(i,j,8)/rhostrw(i,j,k)/aj3z(i,j,k)
         vstr(i  ,j,k  )=vstr(i  ,j,k  ) + wcont(i,j,k)*j2(i,j  ,k)     &
            *mapfct(i,j,8)/rhostrw(i,j,k)/aj3z(i,j,k)
         vstr(i  ,j,k-1)=vstr(i  ,j,k-1) + wcont(i,j,k)*j2(i,j  ,k)     &
            *mapfct(i,j,8)/rhostrw(i,j,k)/aj3z(i,j,k)
         vstr(i,j+1,k  )=vstr(i,j+1,k  ) + wcont(i,j,k)*j2(i,j+1,k)     &
            *mapfct(i,j,8)/rhostrw(i,j,k)/aj3z(i,j,k)
         vstr(i,j+1,k-1)=vstr(i,j+1,k-1) + wcont(i,j,k)*j2(i,j+1,k)     &
            *mapfct(i,j,8)/rhostrw(i,j,k)/aj3z(i,j,k)
         w(i,j,k)=w(i,j,k) + wcont(i,j,k)/aj3z(i,j,k)
         wcont(i,j,k) = 0.0
        END DO
      END DO
    END DO

! No. 2
!
!    IF (mp_opt > 0) THEN
!      !CALL acct_interrupt(mp_acct)
!      CALL mpsendrecv2dew(ustr, nx, ny, nz, ebc, wbc, 1, tem1)
!      CALL mpsendrecv2dns(ustr, nx, ny, nz, nbc, sbc, 1, tem1)
!
!      CALL mpsendrecv2dew(vstr, nx, ny, nz, ebc, wbc, 2, tem1)
!      CALL mpsendrecv2dns(vstr, nx, ny, nz, nbc, sbc, 2, tem1)
!      !CALL acct_stop_inter
!    END IF

    DO k= 1,nz-1
      DO j= 1,ny
        DO i= 1,nx-1
!         vstr(i,j,k)=v(i,j,k)*rhostrv(i,j,k)

          v(i,j,k)=v(i,j,k)+vstr(i,j,k)*rhostrv(i,j,k)
          vstr(i,j,k) = 0.0
        END DO
      END DO
    END DO

    DO k= 1,nz-1
      DO j= 1,ny-1
        DO i= 1,nx
!         ustr(i,j,k)=u(i,j,k)*rhostru(i,j,k)

          u(i,j,k)=u(i,j,k)+ustr(i,j,k)*rhostru(i,j,k)
          ustr(i,j,k)=0.0
        END DO
      END DO
    END DO

  END IF

  RETURN
END SUBROUTINE ADWCONTRA
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE ADVBCWCONT               ######
!######                                                      ######
!######                Copyright (c) 1992-1994               ######
!######  Center for the Analysis and Prediction of Storms    ######
!######  University of Oklahoma.  All rights reserved.       ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE ADVBCWCONT(nx,ny,nz,wcont)
!
!-----------------------------------------------------------------------
!
! PURPOSE:
!
!     Perform the adjoint of
!     Set the top and bottom boundary conditions for wcont.
!
!-----------------------------------------------------------------------
!
! AUTHOR:
!     5/22/96 Jidong Gao
!
!     MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
! INPUT :
!
!       nx       Number of grid points in the x-direction (east/west)
!       ny       Number of grid points in the y-direction (north/south)
!       nz       Number of grid points in the vertical
!
!       wcont    Contravariant vertical velocity (m/s)
!
! OUTPUT:
!
!       wcont    Top and bottom values of wcont.
!
!-----------------------------------------------------------------------
!
! Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE                 ! Force explicit declarations

  INTEGER, INTENT(IN) ::  nx,ny,nz  ! Number of grid points in x, y and z
                                    ! directions

  REAL, INTENT(INOUT) :: wcont (nx,ny,nz)   ! Contravariant vertical velocity (m/s)
!
!-----------------------------------------------------------------------
!
! Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k
!
!-----------------------------------------------------------------------
!
!     Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'bndry.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
!-----------------------------------------------------------------------
!
! Set the top boundary condition
!
!-----------------------------------------------------------------------
!
  IF(tbc == 1) THEN             ! Rigid lid boundary condition

    DO j=1,ny-1
      DO i=1,nx-1
!         wcont(i,j,nz)=-wcont(i,j,nz-2)
!         wcont(i,j,nz-1)=0.0

          wcont(i,j,nz-2)=wcont(i,j,nz-2)-wcont(i,j,nz)
          wcont(i,j,nz)=0.0
          wcont(i,j,nz-1)=0.0
      END DO
    END DO

  ELSE IF(tbc == 2) THEN         ! Periodic boundary condition.

    DO j=1,ny-1
      DO i=1,nx-1
!         wcont(i,j,nz)=wcont(i,j,3)

          wcont(i,j,3)=wcont(i,j,3)+wcont(i,j,nz)
          wcont(i,j,nz)=0.0
      END DO
    END DO

  ELSE IF(tbc == 3.OR.tbc == 4) THEN ! Zero normal gradient condition.

    DO j=1,ny-1
      DO i=1,nx-1
!         wcont(i,j,nz)=wcont(i,j,nz-1)

          wcont(i,j,nz-1)=wcont(i,j,nz-1)+wcont(i,j,nz)
          wcont(i,j,nz) = 0.0
      END DO
    END DO

  ELSE

    WRITE(6,900) 'VBCWCONT', tbc
    CALL arpsstop ("",1)
!    stop

  END IF
!
!-----------------------------------------------------------------------
!
! Set the bottom boundary condition
!
!-----------------------------------------------------------------------
!
  IF(bbc == 1) THEN             ! Non-penetrative ground condition

    DO j=1,ny-1
      DO i=1,nx-1
!         wcont(i,j,1)=-wcont(i,j,3)
!         wcont(i,j,2)=0.0

          wcont(i,j,3)=wcont(i,j,3)-wcont(i,j,1)
          wcont(i,j,1)=0.0
          wcont(i,j,2)=0.0
      END DO
    END DO

  ELSE IF(bbc == 2) THEN         ! Periodic boundary condition.

    DO j=1,ny-1
      DO i=1,nx-1
!         wcont(i,j,1)=wcont(i,j,nz-2)

          wcont(i,j,nz-2)=wcont(i,j,nz-2)+wcont(i,j,1)
          wcont(i,j,1) = 0.0
      END DO
    END DO

  ELSE IF(bbc == 3) THEN         ! Zero normal gradient condition.

    DO j=1,ny-1
      DO i=1,nx-1
!         wcont(i,j,1)=wcont(i,j,2)

          wcont(i,j,2)=wcont(i,j,2)+wcont(i,j,1)
          wcont(i,j,1)= 0.0
      END DO
    END DO

  ELSE

    write(6,900) 'VBCWCONT', bbc
    CALL arpsstop('Wrong bbc in VBCWCONT.',1)
!    STOP
  END IF

  900   format(1x,'Invalid boundary condition option found in ',a,      &
              /1x,'The option was ',i3,' Job stopped.')

  RETURN
END SUBROUTINE ADVBCWCONT

