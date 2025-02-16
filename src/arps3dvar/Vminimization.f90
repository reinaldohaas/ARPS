!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE MINIMIZATION               ######
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
!  Main driver for gradient check and minimization.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!
!  Jidong Gao, CAPS, July, 2000
!
!  MODIFICATIONS:
!
!  Yunheng Wang (11/01/2007)
!  Added two more temporary arrays for parallelizing recursive filter.
!
!-----------------------------------------------------------------------
!
SUBROUTINE minimization(nx,ny,nz,                                       &
           nvar,nvarradin,nvarrad,nzua,nzrdr,nzret,                     &
           mapfct,j1,j2,j3,aj3x,aj3y,aj3z,j3inv,rhostr,                 &
           mxsng,mxua,mxrad,mxcolrad,mxret,mxcolret,                    &
           nsrcsng,nsrcua,nsrcrad,nsrcret,ncat,                         &
           mxpass,ipass,iwstat,                                         &
           xs,ys,zs,x,y,z,zp,hterain,                                   &
           icatg,xcor,nam_var,                                          &
           indexsng,usesng,xsng,ysng,hgtsng,thesng,                     &
           odifsng,qobsng,qualsng,isrcsng,icatsng,nobsng,               &
           indexua,useua,xua,yua,hgtua,theua,                           &
           odifua,qobsua,qualua,isrcua,nlevsua,nobsua,                  &
           indexrad,userad,elvrad,xradc,yradc,                          &
           distrad,uazmrad,vazmrad,hgtradc,theradc,dsdr,dhdr,           &
           obsrad,odifrad,qobsrad,qualrad,                              &
           irad,isrcrad,nlevrad,ncolrad,                                &
           xretc,yretc,hgtretc,theretc,                                 &
           odifret,qobsret,qualret,                                     &
           iret,isrcret,nlevret,ncolret,                                &
           srcsng,srcua,srcrad,srcret,                                  &
           ianxtyp,iusesng,iuseua,iuserad,iuseret,                      &
           xyrange,kpvrsq,wlim,zrange,zwlim,                            &
           thrng,rngsqi,knt,wgtsum,zsum,                                &
           corsng,corua,corrad,corret,                                  &
           oanxsng,oanxua,oanxrad,oanxret,                              &
           nz_tab,hqback,qback,                                         &
           nsngfil,nuafil,nradfil,raduvobs,nretfil,                     &
           anx, smcu,smcv,smcw, sinlat, ffu,ffv, tem2d,                 &
           tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,tem9,tem10,          &
           radius_z,istatus)
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx         Number of grid points in the x-direction (east/west)
!    ny         Number of grid points in the y-direction (north/south)
!    nz         Number of grid points in the vertical
!    nvar       Number of analysis variables
!    nvarradin  Number of variables in the obsrad array
!    nvarrad    Number of variables in the odifrad array
!    nzua       Maximumnumber of levels in UA observations
!    nzrdr      Maximum number of levels in a radar column
!    nzret      Maximum number of levels in a retrieval column
!    mxsng      Maximum number of single level data
!    mxua       Maximum number of upper air observations
!    mxrad      Maximum number of radars
!    mxcolrad   Maximum number of radar columns
!    mxcolret   Maximum number of retrieval columns
!    mxpass     Maximum number of iterations
!    npass      Number of iterations
!    iwstat     Status indicator for writing statistics
!
!    xs         x location of scalar pts (m)
!    ys         y location of scalar pts (m)
!    zs         z location of scalar pts (m)
!
!    xsng       x location of single-level data
!    ysng       y location of single-level data
!    hgtsng     elevation of single-level data
!    thesng     theta (potential temperature) of single-level data
!
!    obsng      single-level observations
!    odifsng    difference between single-level obs and analysis
!    qobsng     normalized observation error
!    qualsng    single-level data quality indicator
!    isrcsng    index of single-level data source
!    nobsng     number of single-level observations
!
!    xua        x location of upper air data
!    yua        y location of upper air data
!    hgtua      elevation of upper air data
!    theua      theta (potential temperature) of upper air data
!
!    obsua      upper air observations
!    odifua     difference between upper air obs and analysis
!    qobsua     normalized observation error
!    qualua     upper air data quality indicator
!    isrcua     index of upper air data source
!    nlevsua    number of levels of data for each upper air location
!    nobsua     number of upper air observations
!
!    anx        Analyzed variables (or first guess)
!    qback      Background (first guess) error
!
!    nradfil    number of radar files
!    fradname   file name for radar dataset
!    srcrad     name of radar sources
!
!    latrad   latitude of radar  (degrees N)
!    lonrad   longitude of radar (degrees E)
!    elvrad   elevation of feed horn of radar (m MSL)
!    xradc    x location of radar column
!    yradc    y location of radar column
!    irad     radar number
!    nlevrad  number of levels of radar data in each column
!    distrad  distance of radar column from source radar
!    uazmrad  u-component of radar beam for each column
!    vazmrad  v-component of radar beam for each column
!    hgtradc  height (m MSL) of radar observations
!    obsrad   radar observations
!    oanxrad  analysis (first guess) value at radar data location
!    odifrad  difference between radar observation and analysis
!    qobsrad  normalized observation error
!    qualrad  radar data quality indicator
!    ncolrad  number of radar columns read-in (local number for MPI)
!    istatus  status indicator
!
!    latret   latitude of retrieval radar  (degrees N)
!    lonret   longitude of retrieval radar (degrees E)
!    xretc    x location of retrieval column
!    yretc    y location of retrieval column
!    iret     retrieval number
!    nlevret  number of levels of retrieval data in each column
!    hgtretc  height (m MSL) of retrieval observations
!    obsret   retrieval observations
!    odifret  difference between retrieval observation and analysis
!    qobsret  normalized observation error
!    qualret  retrieval data quality indicator
!    ncolret  number of retr columns read-in
!    istatus  status indicator
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  USE arps_precision
  USE Vmodule_3dvar

  IMPLICIT NONE

  INCLUDE 'phycst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'varpara.inc'

  !INCLUDE '3dvarcputime.inc'

  INCLUDE 'bndry.inc'   ! to include web, ebc, nbc, sbc etc.
  INCLUDE 'globcst.inc' ! to include mp_acct
  INCLUDE 'mp.inc'      ! to include mp_opt, myproc

  REAL :: cputimearray(10)
  COMMON /cputime3dvar/ cputimearray

!-----------------------------------------------------------------------
!
!  Input Sizing Arguments
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz
  INTEGER :: nvar,nvarradin,nvarrad
  INTEGER :: nzua,nzrdr,nzret
  INTEGER :: mxsng,mxua,mxrad,mxcolrad,mxret,mxcolret
  INTEGER :: nsrcsng,nsrcua,nsrcrad,nsrcret,ncat
  INTEGER :: mxpass,ipass
  INTEGER :: nsngfil,nuafil,nradfil,nretfil
  INTEGER :: raduvobs    ! extra flag for radar data inherit from adas

  INTEGER, INTENT(IN) :: indexsng(mxsng), indexua(mxua), indexrad(mxcolrad)
  LOGICAL, INTENT(IN) :: usesng(mxsng), useua(mxua), userad(mxcolrad)
!
!-----------------------------------------------------------------------
!
!  Input Grid Arguments
!
!-----------------------------------------------------------------------
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
  REAL(P) :: j3inv (nx,ny,nz)  ! Inverse of j3
  REAL(P) :: rhostr(nx,ny,nz)  ! Base state density rhobar times j3
!
  REAL(P) :: xs(nx)
  REAL(P) :: ys(ny)
  REAL(P) :: zs(nx,ny,nz)
!
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
  DOUBLE PRECISION,    DIMENSION(:), ALLOCATABLE :: xsng_p, ysng_p
  REAL(P),             DIMENSION(:), ALLOCATABLE:: zsng_1, zsng_2
  INTEGER,             DIMENSION(:), ALLOCATABLE:: ihgtsng
  REAL(P) :: hgtsng(mxsng)
  REAL(P) :: thesng(mxsng)
  REAL(P) :: obsng(nvar,mxsng)
  REAL(P) :: odifsng(nvar,mxsng)
  REAL(P) :: qobsng(nvar,mxsng)
  INTEGER :: qualsng(nvar,mxsng)
  INTEGER :: isrcsng(mxsng)
  INTEGER :: icatsng(mxsng)
  INTEGER :: nobsng
!
  REAL(P) :: xua(mxua)
  REAL(P) :: yua(mxua)
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: xua_p, yua_p
  REAL(P),          DIMENSION(:,:), ALLOCATABLE:: zua_1, zua_2
  INTEGER,          DIMENSION(:,:), ALLOCATABLE:: ihgtua
  REAL(P) :: hgtua(nzua,mxua)
  REAL(P) :: theua(nzua,mxua)
  REAL(P) :: obsua(nvar,nzua,mxua)
  REAL(P) :: odifua(nvar,nzua,mxua)
  REAL(P) :: qobsua(nvar,nzua,mxua)
  INTEGER :: qualua(nvar,nzua,mxua)
  INTEGER :: nlevsua(mxua)
  INTEGER :: isrcua(mxua)
  INTEGER :: nobsua
!
  REAL(P) :: elvrad(mxrad)
  REAL(P) :: xradc(mxcolrad)
  REAL(P) :: yradc(mxcolrad)
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE:: xradc_p, yradc_p
  REAL(P),          DIMENSION(:,:), ALLOCATABLE:: zradc_1, zradc_2
  INTEGER,          DIMENSION(:,:), ALLOCATABLE:: ihgtradc
  REAL(P) :: distrad(mxcolrad)
  REAL(P) :: uazmrad(mxcolrad)
  REAL(P) :: vazmrad(mxcolrad)
  REAL(P) :: hgtradc(nzrdr,mxcolrad)
  REAL(P) :: theradc(nzrdr,mxcolrad)
  REAL(P) :: dsdr(nzrdr,mxcolrad)
  REAL(P) :: dhdr(nzrdr,mxcolrad)
  REAL(P) :: obsrad(nvarradin,nzrdr,mxcolrad)
  REAL(P) :: odifrad(nvarrad,nzrdr,mxcolrad)
  REAL(P) :: qobsrad(nvarrad,nzrdr,mxcolrad)
  INTEGER :: qualrad(nvarrad,nzrdr,mxcolrad)
  INTEGER :: nlevrad(mxcolrad)
  INTEGER :: irad(mxcolrad)
  INTEGER :: isrcrad(0:mxrad)
  INTEGER :: ncolrad
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
  INTEGER :: iusesng(0:nsrcsng,mxpass)
  INTEGER :: iuseua(0:nsrcua,mxpass)
  INTEGER :: iuserad(0:nsrcrad,mxpass)
  INTEGER :: iuseret(0:nsrcret,mxpass)

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
!
  REAL(P) :: oanxsng(nvar,mxsng)
  REAL(P) :: oanxua(nvar,nzua,mxua)
  REAL(P) :: oanxrad(nvarrad,nzrdr,mxcolrad)
  REAL(P) :: oanxret(nvar,nzret,mxcolret)
!
  INTEGER :: nz_tab
  REAL(P) :: hqback(nz_tab)
  REAL(P) :: qback(nvar,nz_tab)
!
  REAL(P),DIMENSION(nz_tab) :: qback_wk
  REAL(P),DIMENSION(nz_tab) :: hqback_wk
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
!  variable declarations for varirational analysis: by gao
!
!-----------------------------------------------------------------------
!
  REAL(P) :: eps,gtol,ftol
  REAL(P) :: cfun
  DOUBLE PRECISION :: cfun_single(nvar+1),cfun_total
  INTEGER :: mp,lp,iprint(2),iflag,point
  LOGICAL :: diagco
  COMMON /va15dd/mp,lp, gtol

  INTEGER :: icall
  REAL(P) :: rchek,gxnn
  REAL(P) :: cfun1

  REAL(P), ALLOCATABLE :: fa(:), fr(:)
!
!-----------------------------------------------------------------------
!
!  end of variable declarations for variational anaylsis:
!
!-----------------------------------------------------------------------
!
  REAL(P), DIMENSION (nx,ny,nz) :: smcu,smcv,smcw
  REAL(P), DIMENSION (nx,ny,nz) :: sinlat
  REAL(P), DIMENSION (nx,ny,nz) :: ffu,ffv
  REAL(P), DIMENSION (nx,ny)    :: tem2d

!
!-----------------------------------------------------------------------
!
!  Work arrays
!
!-----------------------------------------------------------------------
!
  REAL(P) :: tem1(nx,ny,nz), tem2(nx,ny,nz), tem3(nx,ny,nz), tem4(nx,ny,nz)
  REAL(P) :: tem5(nx,ny,nz), tem6(nx,ny,nz), tem7(nx,ny,nz), tem8(nx,ny,nz)
  REAL(P) :: tem9(nx,ny,nz), tem10(nx,ny,nz)
  REAL(P) :: radius_z(nx,ny,nz), rL,rU,rrr

  LOGICAL :: ownsng(mxsng), ownua(mxua)          ! with global index
  LOGICAL :: ownrad(mxcolrad)                    ! with local index
!
!-----------------------------------------------------------------------
!
!  Return status
!
!-----------------------------------------------------------------------
!
  INTEGER :: istatus,ierr
!
!-----------------------------------------------------------------------
!
!  Misc.local variables
!
!-----------------------------------------------------------------------
!
  REAL(P) :: ftabinv,setexp
  INTEGER :: i,j,k,isrc,sngsw,uasw,radsw,retsw, nL,nU,nV
  INTEGER :: num, numsingle
  REAL(P) :: rpass,zrngsq,thrngsq
  REAL(P) :: sum1,sum2,sum3,sum4,sum5,sum6
  REAL(P) :: f_cputime
  REAL(P) :: cpuin,cpuout
  INTEGER :: ount

  !REAL(P),DIMENSION (:,:,:), ALLOCATABLE :: smcu,smcv,smcw
  !REAL(P),DIMENSION (:,:,:), ALLOCATABLE :: sinlat
  !REAL(P),DIMENSION (:,:,:), ALLOCATABLE :: ffu,ffv
  !REAL(P),DIMENSION (:,:),   ALLOCATABLE :: tem2d

  REAL(P) :: omega2, sinclat
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ount = 6
!
!-----------------------------------------------------------------------
!
!  allocate control variable arrays
!
!-----------------------------------------------------------------------
!
  !ALLOCATE ( smcu(nx,ny,nz), STAT = istatus )
  !CALL check_alloc_status(istatus, "minimization:smcu")
  smcu(:,:,:) = 0.0
  !
  !ALLOCATE ( smcv(nx,ny,nz), STAT = istatus )
  !CALL check_alloc_status(istatus, "minimization:smcv")
  smcv(:,:,:) = 0.0
  !
  !ALLOCATE ( smcw(nx,ny,nz), STAT = istatus )
  !CALL check_alloc_status(istatus, "minimization:smcw")
  smcw(:,:,:) = 0.0
  !
  !IF( ALLOCATED(sinlat) ) DEALLOCATE ( sinlat )
  !ALLOCATE ( sinlat (nx,ny,nz) )
  !CALL check_alloc_status(istatus, "minimization:sinlat")
  !
  !IF( ALLOCATED(ffu) ) DEALLOCATE ( ffu )
  !ALLOCATE ( ffu(nx,ny,nz), STAT = istatus )
  !CALL check_alloc_status(istatus, "minimization:ffu")
  ffu(:,:,:) = 0.0
  !
  !IF( ALLOCATED(ffv) ) DEALLOCATE ( ffv )
  !ALLOCATE ( ffv (nx,ny,nz), STAT = istatus )
  !CALL check_alloc_status(istatus, "minimization:ffv")
  ffv(:,:,:) = 0.0
  !
  !IF( ALLOCATED(tem2d) ) DEALLOCATE ( tem2d )
  !ALLOCATE ( tem2d(nx,ny), STAT = istatus )
  !CALL check_alloc_status(istatus, "minimization:tem2d")


!
!------------------------calculate coriolis costant-----------------------
!
   CALL gtsinlat(nx,ny,x,y, tem2d, tem1,tem2, tem3)
!
! Added the effects of spatial gradient of map factor on the coriolis force
!
  IF( coriopt == 3 .OR. coriopt == 4) THEN
    DO j = 1,ny
      DO i = 1,nx
        tem1(i,j,1) = tem2d(i,j)/SQRT(1-tem2d(i,j)**2)  ! tan(lat)
      END DO
    END DO

    DO k = 1,nz
      DO j = 1,ny-1
        DO i = 1,nx-1
          !
          ! fm = U*My - V*Mx + U*TAN(lat)/eradius
          !
          tem2(i,j,k) = 0.5*( (anx(i,j,k,1)+anx(i+1,j,k,1))*                &
                              ((mapfct(i,j+1,3)-mapfct(i,j,3))/dy)          &
                             -(anx(i,j,k,2)+anx(i,j+1,k,2))*                &
                              ((mapfct(i+1,j,2)-mapfct(i,j,2))/dx)          &
                             +((anx(i,j,k,1)+anx(i+1,j,k,1))*tem1(i,j,1))/  &
                              eradius)
        END DO
      END DO
    END DO

  END IF
!
  omega2 = 2.0* omega
  IF( coriopt == 1 ) THEN

    sinclat = SIN( ATAN(1.0)/45.0 * ctrlat )
    DO k = 1,nz
      DO j=1,ny
        DO i=1,nx
          sinlat(i,j,k) = omega2* sinclat
        END DO
      END DO
    END DO
  ELSE IF( coriopt == 2 ) THEN

    sinclat = SIN( ATAN(1.0)/45.0 * ctrlat )
    DO k = 1,nz
      DO j=1,ny
        DO i=1,nx
          sinlat(i,j,k) = omega2* sinclat
        END DO
      END DO
    END DO

  ELSE IF( coriopt == 3 ) THEN

    DO k = 1,nz
      DO j=1,ny-1
        DO i=1,nx-1
          sinlat(i,j,k) = omega2*tem2d(i,j) +tem2(i,j,k)
        END DO
      END DO
    END DO

  ELSE IF( coriopt == 4 ) THEN

    DO k = 1,nz
      DO j=1,ny-1
        DO i=1,nx-1
          sinlat(i,j,k) = omega2*tem2d(i,j) +tem2(i,j,k)
        END DO
      END DO
    END DO

  END IF
!
!-------------------------------------------------------
!
  CALL rhouvw(nx,ny,nz,rhostr,rhostru,rhostrv,rhostrw)

!
!    PRINT*,'qbackin=',(qback(1,i),i=1,nz_tab)
!
   ! Background error for U
   num = 0
   DO isrc = 1, nz_tab
     IF( qback(1,isrc) > -900.0) THEN
       num = num+1
       qback_wk (num) =  qback(1,isrc)
       hqback_wk(num) =  hqback(isrc)
     END IF
   END DO

   DO k = 1,nz
     DO j = 1,ny
       DO i = 1,nx
         CALL linear(num,hqback_wk,qback_wk,zs(i,j,k),gdu_err(i,j,k) )
       END DO
     END DO
   END DO

   ! Background error for V
   num = 0
   DO isrc = 1, nz_tab
     IF( qback(2,isrc) > -900.0) THEN
       num = num+1
       qback_wk(num)  = qback(2,isrc)
       hqback_wk(num) = hqback(isrc)
     end if
   END DO

   DO k = 1,nz
     DO j = 1,ny
       DO i = 1,nx
         CALL linear(num,hqback_wk,qback_wk,zs(i,j,k),gdv_err(i,j,k) )
       END DO
     END DO
   END DO

   ! Background error for P
   num = 0
   DO isrc = 1, nz_tab
      IF( qback(3,isrc) > -900.0) THEN
        num = num+1
        qback_wk(num)  = qback(3,isrc)
        hqback_wk(num) = hqback(isrc)
      END IF
   END DO

   DO k = 1,nz
     DO j = 1,ny
       DO i = 1,nx
         CALL linear(num,hqback_wk,qback_wk,zs(i,j,k),gdp_err(i,j,k) )
       END DO
     END DO
   END DO

   ! Background error for T
   num = 0
   DO isrc = 1, nz_tab
     IF( qback(4,isrc) > -900.0) THEN
       num = num+1
       qback_wk(num) =  qback(4,isrc)
       hqback_wk(num) = hqback(isrc)
     END IF
   END DO

   DO k = 1,nz
     DO j = 1,ny
       DO i = 1,nx
         CALL linear(num,hqback_wk,qback_wk,zs(i,j,k),gdt_err(i,j,k) )
       END DO
     END DO
   END DO

   ! Background error for QV
   num = 0
   DO isrc = 1, nz_tab
     IF( qback(5,isrc) > -900.0) THEN
       num = num+1
       qback_wk(num) =  qback(5,isrc)
       hqback_wk(num) = hqback(isrc)
     END IF
   END DO

   DO k = 1,nz
     DO j = 1,ny
       DO i = 1,nx
         CALL linear(num,hqback_wk,qback_wk,zs(i,j,k),gdq_err(i,j,k) )
       END DO
     END DO
   END DO

   ! Background error for W
!
!  num = 0
!  DO isrc = 1, nz_tab
!     if( qback(1,isrc) > -900.0) then
!     num = num+1
!     qback_wk(num) =  qback(1,isrc)
!    hqback_wk(num) = hqback(isrc)
!     end if
!  END DO
!
!  DO k = 1,nz
!   DO j = 1,ny
!    DO i = 1,nx
!      call linear(num,hqback_wk,qback_wk,zs(i,j,k),gdw_err(i,j,k))
!    END DO
!   END DO
!  END DO
!
   DO k = 1,nz
     DO j = 1,ny
       DO i = 1,nx
         gdw_err(i,j,k) = gdu_err(i,j,k)
       END DO
     END DO
   END DO

!
!----------------if ( cntl_var==0) u & v as control variable-------------
!
  IF(cntl_var == 0) THEN

!    IF( .FALSE. ) THEN   ! Fixed error matrices
!      gdu_err(:,:,:) = 1.0
!      gdv_err(:,:,:) = 1.0
!      gdp_err(:,:,:) = 1.0
!      gdt_err(:,:,:) = 1.0
!      gdq_err(:,:,:) = 1.0
!      gdw_err(:,:,:) = 1.0
!    ELSE
      DO k = 1,nz
        DO j = 1,ny
          DO i = 1,nx
            gdu_err(i,j,k) = SQRT( gdu_err(i,j,k) )
            gdv_err(i,j,k) = SQRT( gdv_err(i,j,k) )
            gdp_err(i,j,k) = SQRT( gdp_err(i,j,k) )
            gdt_err(i,j,k) = SQRT( gdt_err(i,j,k) )
            gdq_err(i,j,k) = SQRT( gdq_err(i,j,k) )
            gdw_err(i,j,k) = SQRT( gdw_err(i,j,k) )
          END DO
        END DO
      END DO
!    ENDIF

  ELSE   ! Stream function & potential as control variables

    DO k = 1,nz
      DO j = 1,ny
        DO i = 1,nx
        gdu_err(i,j,k) = SQRT( gdu_err(i,j,k) )*(dx+dy)/2.
        gdv_err(i,j,k) = SQRT( gdv_err(i,j,k) )*(dx+dy)/2.
        gdp_err(i,j,k) = SQRT( gdp_err(i,j,k) )
        gdt_err(i,j,k) = SQRT( gdt_err(i,j,k) )
        gdq_err(i,j,k) = SQRT( gdq_err(i,j,k) )
        gdw_err(i,j,k) = SQRT( gdw_err(i,j,k) )
        END DO
      END DO
    END DO

  END IF

!-----------------------------------------------------------------------
!
! Allocate arrays and set flags
!
!-----------------------------------------------------------------------

  ownsng(:)  = .FALSE.
  ownua (:)  = .FALSE.
  ownrad(:)  = .FALSE.

  DO i=1,numctr
    ctrv(i)=0.0
    xgus(i)=0.0
    grad(i)=0.0
  END DO

  sngsw = 0
  uasw  = 0
  radsw = 0
  retsw = 0
  IF( nsngfil /= 0 ) sngsw = 1
  IF( nuafil  /= 0 ) uasw  = 1
  IF( nradfil /= 0 .AND. raduvobs /= 0) radsw = 1
  IF( nretfil /= 0 ) retsw = 1

!
!  caculating some parameters for interpolation ! added by Jidong Gao
!

  !
  ! SNG
  !
  DO i = 1, nobsng
    icatsng(i) = 1
    IF (indexsng(i) == myproc) ownsng(i) = .TRUE.
  END DO

  IF( ALLOCATED(xsng_p) ) DEALLOCATE (xsng_p)
  ALLOCATE ( xsng_p(mxsng), STAT=istatus )
  CALL check_alloc_status(istatus, "minimization:xsng_p")

  IF( ALLOCATED(ysng_p) ) DEALLOCATE (ysng_p)
  ALLOCATE ( ysng_p(mxsng), STAT=istatus )
  CALL check_alloc_status(istatus, "minimization:ysng_p")

  IF( ALLOCATED(zsng_1) ) DEALLOCATE (zsng_1)
  ALLOCATE ( zsng_1(mxsng), STAT=istatus )
  CALL check_alloc_status(istatus, "minimization:zsng_l")

  IF( ALLOCATED(zsng_2) ) DEALLOCATE (zsng_2)
  ALLOCATE ( zsng_2(mxsng), STAT=istatus )
  CALL check_alloc_status(istatus, "minimization:zsng_2")

  IF( ALLOCATED(ihgtsng) ) DEALLOCATE (ihgtsng)
  ALLOCATE ( ihgtsng(mxsng), STAT=istatus )
  CALL check_alloc_status(istatus, "minimization:ihgtsng")

  CALL map_to_mod2(nx,ny,mxsng,nobsng,usesng,xs,ys,xsng,ysng,xsng_p,ysng_p)

  CALL map_to_modz(1, mxsng, icatsng, nobsng,nx,ny,nz,zs,               &
                   xsng_p,ysng_p,hgtsng,ihgtsng,zsng_1,zsng_2)

  DO i = 1, nobsua
    IF (indexua(i) == myproc) ownua(i) = .TRUE.
  END DO

  IF( allocated(xua_p) ) deallocate (xua_p)
  ALLOCATE ( xua_p(mxua), STAT=istatus )
  CALL check_alloc_status(istatus, "minimization:xua_p")

  IF( allocated(yua_p) ) deallocate (yua_p)
  ALLOCATE ( yua_p(mxua), STAT=istatus )
  CALL check_alloc_status(istatus, "minimization:yua_pua_p")

  IF( allocated(zua_1) ) deallocate (zua_1)
  ALLOCATE ( zua_1(nzua,mxua), STAT=istatus )
  CALL check_alloc_status(istatus, "minimization:zua_1")

  IF( allocated(zua_2) ) deallocate (zua_2)
  ALLOCATE ( zua_2(nzua,mxua), STAT=istatus )
  CALL check_alloc_status(istatus, "minimization:zua_2")

  IF( allocated( ihgtua) ) deallocate (ihgtua)
  ALLOCATE ( ihgtua(nzua,mxua), STAT=istatus )
  CALL check_alloc_status(istatus, "minimization:ihgtua")

  CALL map_to_mod2(nx,ny,mxua,nobsua,useua,xs,ys,xua,yua,xua_p,yua_p)
  CALL map_to_modz(nzua, mxua, nlevsua, nobsua, nx,ny,nz,zs,            &
                   xua_p,yua_p, hgtua, ihgtua,zua_1,zua_2)

  DO i = 1, ncolrad   ! local number
    IF (indexrad(i) == myproc) ownrad(i)  = .TRUE.
  END DO

  IF( allocated(xradc_p) ) deallocate (xradc_p)
  ALLOCATE ( xradc_p(mxcolrad), STAT=istatus )
  CALL check_alloc_status(istatus, "minimization:xradc_p")

  IF( allocated(yradc_p) ) deallocate (yradc_p)
  ALLOCATE ( yradc_p(mxcolrad), STAT=istatus )
  CALL check_alloc_status(istatus, "minimization:yradc_p")

  IF( allocated(zradc_1) ) deallocate (zradc_1)
  ALLOCATE ( zradc_1(nzrdr,mxcolrad), STAT=istatus )
  CALL check_alloc_status(istatus, "minimization:zradc_1")

  IF( allocated(zradc_2) ) deallocate (zradc_2)
  ALLOCATE ( zradc_2(nzrdr,mxcolrad), STAT=istatus )
  CALL check_alloc_status(istatus, "minimization:zrad_2")

  IF( allocated(ihgtradc) ) deallocate (ihgtradc)
  ALLOCATE ( ihgtradc(nzrdr,mxcolrad), STAT=istatus )
  CALL check_alloc_status(istatus, "minimization:ihgtradc")

  CALL map_to_mod2(nx,ny,mxcolrad,ncolrad,userad,xs,ys,xradc,yradc,     &
                                                       xradc_p,yradc_p)
  CALL map_to_modz(nzrdr, mxcolrad, nlevrad, ncolrad,nx,ny,nz,zs,       &
                   xradc_p,yradc_p,hgtradc,ihgtradc,zradc_1,zradc_2)

  !BEGIN: to specifiy vertical influence radius
  IF(vradius_opt ==1) THEN  ! vradius is in grid points
    DO j =1, ny-1
     DO i =1, nx-1
      DO k =1, nz-1
        radius_z(i,j,k)=vradius(ipass)
      END DO
     END DO
    END DO
  ELSE IF(vradius_opt ==2) THEN  ! vradius in in unit of km
    DO j =1, ny-1
     DO i =1, nx-1
      rrr=vradius(ipass)*1000.0  !user-specified vertical radius in meters
      DO k =2, nz-1
        nL=0; nU=0; rL =0.0; rU=0.0
        ! count upward
        DO nV =k, nz-1
          rU = zp(i,j,nV) - zp(i,j,k)
          IF ( rU >= rrr ) THEN
            nU = nV - k
            exit ! found nU, exit the do loop
          END IF
        END DO
        IF (nU == 0 .and. nV >= nz-1 ) nU = nz -1 - k
        ! count downward
        DO nV = k, 2, -1
          rL = zp(i,j,k) - zp(i,j,nV) 
          IF ( rL >= rrr ) THEN
            nL = k - nV
            exit ! found nL, exit the do loop
          END IF
        END DO
        IF (nL == 0 .and. nV <= 2 ) nL = k - 2 

        radius_z(i,j,k) = ( nL * max(rrr, rL) + nU * max(rrr, rU) ) / ( rU + rL )
      END DO
      radius_z(i,j,1) = radius_z(i,j,2)
      radius_z(i,j,nz) = radius_z(i,j,nz-1)
     END DO
    END DO
  END IF ! IF (vradius_opt == 1 or 2)
  IF (myproc == 0 ) THEN
    WRITE(*,*) 'The vertical influence radius is specified in grid points as:'
    WRITE(*,*) '---- The followings are k, zp, dzp, vradius at i=3,j=3 ----'
    i=3; j =3  
    DO k =1, nz-1
      WRITE(*,*) k, zp(i,j,k), zp(i,j,k+1)-zp(i,j,k), radius_z(i,j,k)
    END DO
  END IF
  !END  : to specifiy vertical influence radius

  IF(chk_opt == 1) THEN
    ipass = 1

    CALL mpmaxi(sngsw)    ! All processors must have the same flags
    CALL mpmaxi(uasw)
    CALL mpmaxi(radsw)

    tem7=0.0;tem8=0.0;tem9=0.0;tem4=0.0;tem5=0.0
    CALL scale_factor(nx,ny,nz,gdscal,ipass_filt(1),                    &
            hradius(1),radius_z(1,1,1),tem7,tem8,tem9,tem4,tem5)

    rchek=1.0E-1
    ALLOCATE (fa(20))
    ALLOCATE (fr(20))
!
!-----------------------------------------------------------------------
!
!      OBTAIN THE COST FUNCTION AND THE GRADIENT.
!
!-----------------------------------------------------------------------
!
    CALL costf(ount,numctr,ctrv,cfun_single,                            &
          gdu_err,gdv_err,gdp_err,gdt_err,gdq_err,gdw_err,              &
          u_ctr,v_ctr,p_ctr,t_ctr,q_ctr,w_ctr, psi, phi,                &
          gdscal, nx,ny,nz,                                             &
          nvar,nvarradin,nvarrad,nzua,nzrdr,nzret,                      &
          mapfct,j1,j2,j3,aj3x,aj3y,aj3z,j3inv,rhostr,                  &
          rhostru, rhostrv, rhostrw, div3,                              &
          mxsng,mxua,mxrad,mxcolrad,mxret,mxcolret,                     &
          nsrcsng,nsrcua,nsrcrad,nsrcret,ncat,                          &
          mxpass,ipass,iwstat,xs,ys,zs,x,y,z,zp,hterain,                &
          icatg,xcor,nam_var,                                           &
          ownsng,usesng,xsng,ysng,hgtsng,thesng,                        &
          obsng,odifsng,qobsng,qualsng,isrcsng,icatsng,nobsng,          &
          ownua,useua,xua,yua,hgtua,theua,                              &
          obsua,odifua,qobsua,qualua,isrcua,nlevsua,nobsua,             &
          ownrad,userad,elvrad,xradc,yradc,                             &
          distrad,uazmrad,vazmrad,hgtradc,theradc,dsdr,dhdr,            &
          obsrad,odifrad,qobsrad,qualrad,                               &
          irad,isrcrad,nlevrad,ncolrad,                                 &
          xretc,yretc,hgtretc,theretc,                                  &
          obsret,odifret,qobsret,qualret,                               &
          iret,isrcret,nlevret,ncolret,                                 &
          srcsng,srcua,srcrad,srcret,                                   &
          ianxtyp,iusesng(0,ipass),iuseua(0,ipass),                     &
          iuserad(0,ipass),iuseret(0,ipass),                            &
          xyrange,kpvrsq,wlim,zrange,zwlim,                             &
          thrng,rngsqi,knt,wgtsum,zsum,                                 &
          corsng,corua,corrad,corret,                                   &
          xsng_p,ysng_p,ihgtsng,xua_p,yua_p,ihgtua,                     &
          xradc_p,yradc_p,ihgtradc,zsng_1,zsng_2,                       &
          zua_1,zua_2,zradc_1,zradc_2,                                  &
          oanxsng,oanxua,oanxrad,oanxret,                               &
          sngsw,uasw,radsw,retsw,                                       &
          ipass_filt(ipass),hradius(ipass),radius_z(1,1,1),             &
          div_opt,cntl_var,smth_flag,                                   &
          wgt_div_h(ipass),wgt_div_v(ipass),wgt_smth(ipass),            &
          thermo_opt,wgt_thermo(ipass),sinlat,ffu,ffv,                  &
          anx,tem1,tem2,tem3,tem4,smcu,smcv,smcw,                       &
          tem5,tem6,tem7,tem8,istatus)

    cfun_total=0.0
    DO i=1,nvar+1
      cfun_total=cfun_total+cfun_single(i)
    ENDDO

    CALL mpsumdp(cfun_total,1)

    cfun = cfun_total
    IF (myproc == 0) PRINT *, 'cfun  ==== ',cfun

    CALL gradt(numctr,ctrv,grad,                                        &
          gdu_err,gdv_err,gdp_err,gdt_err,gdq_err,gdw_err,              &
          u_ctr,v_ctr,p_ctr,t_ctr,q_ctr,w_ctr, psi, phi,                &
          gdscal, nx,ny,nz,                                             &
          nvar,nvarradin,nvarrad,nzua,nzrdr,nzret,                      &
          mapfct,j1,j2,j3,aj3x,aj3y,aj3z,j3inv,rhostr,                  &
          rhostru, rhostrv, rhostrw, div3,                              &
          mxsng,mxua,mxrad,mxcolrad,mxret,mxcolret,                     &
          nsrcsng,nsrcua,nsrcrad,nsrcret,ncat,                          &
          mxpass,ipass,iwstat,xs,ys,zs,x,y,z,zp,hterain,                &
          icatg,xcor,nam_var,xsng,ysng,hgtsng,thesng,                   &
          obsng,odifsng,qobsng,qualsng,isrcsng,icatsng,nobsng,          &
          xua,yua,hgtua,theua,                                          &
          obsua,odifua,qobsua,qualua,isrcua,nlevsua,nobsua,             &
          elvrad,xradc,yradc,                                           &
          distrad,uazmrad,vazmrad,hgtradc,theradc,dsdr,dhdr,            &
          obsrad,odifrad,qobsrad,qualrad,                               &
          irad,isrcrad,nlevrad,ncolrad,                                 &
          xretc,yretc,hgtretc,theretc,                                  &
          obsret,odifret,qobsret,qualret,                               &
          usesng,useua,userad,                                          &
          iret,isrcret,nlevret,ncolret,                                 &
          srcsng,srcua,srcrad,srcret,                                   &
          ianxtyp,iusesng(0,ipass),iuseua(0,ipass),                     &
          iuserad(0,ipass),iuseret(0,ipass),                            &
          xyrange,kpvrsq,wlim,zrange,zwlim,                             &
          thrng,rngsqi,knt,wgtsum,zsum,                                 &
          corsng,corua,corrad,corret,                                   &
          xsng_p,ysng_p,ihgtsng,xua_p,yua_p,ihgtua,                     &
          xradc_p,yradc_p,ihgtradc,zsng_1,zsng_2,                       &
          zua_1,zua_2,zradc_1,zradc_2,                                  &
          oanxsng,oanxua,oanxrad,oanxret,                               &
          sngsw,uasw,radsw,retsw,                                       &
          ipass_filt(ipass),hradius(ipass),radius_z(1,1,1),             &
          div_opt,cntl_var,smth_flag,                                   &
          wgt_div_h(ipass),wgt_div_v(ipass),wgt_smth(ipass),            &
          thermo_opt,wgt_thermo(ipass),sinlat,ffu,ffv,                  &
          anx,tem1,tem2,tem3,tem4,smcu,smcv,smcw,                       &
          tem5,tem6,tem7,tem8,tem9,tem10,istatus)

    gxnn=0.
    DO i=1,numctr
      gxnn=gxnn+grad(i)*grad(i)
    END DO
    CALL mptotal(gxnn)

    IF (myproc == 0) PRINT *,'gxnn  ==== ',gxnn

    DO j=1,20
      rchek=rchek*1.0E-1
      DO i=1,numctr
        xgus(i)=ctrv(i)+rchek*grad(i)
      END DO
!
!-----------------------------------------------------------------------
!
!      OBTAIN THE COST FUNCTION AND THE GRADIENT.
!
!-----------------------------------------------------------------------
!
      CALL costf(ount,numctr,xgus,cfun_single,                          &
            gdu_err,gdv_err,gdp_err,gdt_err,gdq_err,gdw_err,            &
            u_ctr,v_ctr,p_ctr,t_ctr,q_ctr,w_ctr, psi, phi,              &
            gdscal, nx,ny,nz,                                           &
            nvar,nvarradin,nvarrad,nzua,nzrdr,nzret,                    &
            mapfct,j1,j2,j3,aj3x,aj3y,aj3z,j3inv,rhostr,                &
            rhostru, rhostrv, rhostrw, div3,                            &
            mxsng,mxua,mxrad,mxcolrad,mxret,mxcolret,                   &
            nsrcsng,nsrcua,nsrcrad,nsrcret,ncat,                        &
            mxpass,ipass,iwstat,xs,ys,zs,x,y,z,zp,hterain,              &
            icatg,xcor,nam_var,                                         &
            ownsng,usesng, xsng,ysng,hgtsng,thesng,                     &
            obsng,odifsng,qobsng,qualsng,isrcsng,icatsng,nobsng,        &
            ownua,useua,xua,yua,hgtua,theua,                            &
            obsua,odifua,qobsua,qualua,isrcua,nlevsua,nobsua,           &
            ownrad,userad,elvrad,xradc,yradc,                           &
            distrad,uazmrad,vazmrad,hgtradc,theradc,dsdr,dhdr,          &
            obsrad,odifrad,qobsrad,qualrad,                             &
            irad,isrcrad,nlevrad,ncolrad,                               &
            xretc,yretc,hgtretc,theretc,                                &
            obsret,odifret,qobsret,qualret,                             &
            iret,isrcret,nlevret,ncolret,                               &
            srcsng,srcua,srcrad,srcret,                                 &
            ianxtyp,iusesng(0,ipass),iuseua(0,ipass),                   &
            iuserad(0,ipass),iuseret(0,ipass),                          &
            xyrange,kpvrsq,wlim,zrange,zwlim,                           &
            thrng,rngsqi,knt,wgtsum,zsum,                               &
            corsng,corua,corrad,corret,                                 &
            xsng_p,ysng_p,ihgtsng,xua_p,yua_p,ihgtua,                   &
            xradc_p,yradc_p,ihgtradc,zsng_1,zsng_2,                     &
            zua_1,zua_2,zradc_1,zradc_2,                                &
            oanxsng,oanxua,oanxrad,oanxret,                             &
            sngsw,uasw,radsw,retsw,                                     &
            ipass_filt(ipass),hradius(ipass),radius_z(1,1,1),           &
            div_opt,cntl_var,smth_flag,                                 &
            wgt_div_h(ipass),wgt_div_v(ipass),wgt_smth(ipass),          &
            thermo_opt,wgt_thermo(ipass),sinlat,ffu,ffv,                &
            anx,tem1,tem2,tem3,tem4,smcu,smcv,smcw,                     &
            tem5,tem6,tem7,tem8,istatus)

      cfun_total=0
      DO i=1,nvar+1
        cfun_total=cfun_total+cfun_single(i)
      END DO

      CALL mpsumdp(cfun_total,1)

      cfun1=cfun_total

      IF (myproc == 0) PRINT *, 'cfun1 ==== ',cfun1,cfun1-cfun

      fa(j)=rchek
      fr(j)=(cfun1-cfun)/(gxnn*rchek)
    END DO

    IF (myproc == 0) THEN
      DO j=1,16
        WRITE(*,*) '  Rchek=',fa(j),'  fr==',fr(j)
      END DO
    END IF

    DEALLOCATE (fa, fr)
    CALL arpsstop(' ',1)

  END IF

  IF(assim_opt == 1) THEN
!
!-----------------------------------------------------------------------
!
!  Data assimilation begins here.
!
!-----------------------------------------------------------------------
!
    DO i=1,numctr
      ctrv(i)=0.0
    END DO
!
!-----------------------------------------------------------------------
!
!  Loop through ipass analysis iterations
!
!-----------------------------------------------------------------------
!
!    DO ipass=1,npass
      icall  = 0
      eps    = 1.0E-12
      ftol   = 1.0E-4
      !iprint = 0
      iprint(:) = -1
      iflag  = 0
      diagco = .false.
      grad   = 0.0

      tem7=0.0;tem8=0.0;tem9=0.0;tem4=0.0;tem5=0.0
      CALL scale_factor(nx,ny,nz,gdscal,ipass_filt(ipass),              &
             hradius(ipass),radius_z(1,1,1),tem7,tem8,tem9,tem4,tem5)
!
!-----------------------------------------------------------------------
!
!  Set single-level usage switch based in iusesng
!
!-----------------------------------------------------------------------
!
      IF (myproc == 0) WRITE(6,'(/a,i4)') ' Source usage switches for pass: ',ipass
      IF (myproc == 0) WRITE(6,'(/a,i4)') '   Single level data'
      sngsw=0

      IF( nsngfil /= 0 ) THEN
        DO isrc=1,nsrcsng
          IF(iusesng(isrc,ipass) > 0) THEN
            sngsw=1
            IF (myproc == 0) WRITE(6,'(a,a)') '        Using          ',srcsng(isrc)
          END IF
        END DO
      END IF
      IF(sngsw == 0 .AND. myproc == 0) WRITE(6,'(a)') '         none'
!
!-----------------------------------------------------------------------
!
!  Set multiple-level usage switch based in iuseua
!
!-----------------------------------------------------------------------
!
      IF (myproc == 0) WRITE(6,'(/a,i4)') '   Multiple level data'
      uasw=0
      IF( nuafil /= 0 ) THEN
        DO isrc=1,nsrcua
          IF(iuseua(isrc,ipass) > 0) THEN
            uasw=1
            IF (myproc == 0) WRITE(6,'(a,a)') '        Using          ',srcua(isrc)
          END IF
        END DO
      END IF
      IF(uasw == 0 .AND. myproc == 0) WRITE(6,'(a)') '        none'
!
!-----------------------------------------------------------------------
!
!  Set radar-data usage switch based in iuserad
!
!-----------------------------------------------------------------------
!
      IF (myproc == 0) WRITE(6,'(/a,i4)') '   Raw radar data'
      radsw=0
      IF( nradfil /= 0 .AND. raduvobs /= 0 ) THEN
        DO isrc=1,nsrcrad
          IF(iuserad(isrc,ipass) > 0) THEN
            radsw=1
            IF (myproc == 0) WRITE(6,'(a,a)') '        Using          ',srcrad(isrc)
          END IF
        END DO
      END IF
      IF(radsw == 0 .AND. myproc == 0) WRITE(6,'(a)') '         none'
!
!-----------------------------------------------------------------------
!
!  Set ret usage switch based in iuseret
!
!-----------------------------------------------------------------------
!
      IF (myproc == 0) WRITE(6,'(/3x,a,i4)') 'Retrieval radar data'
      retsw=0
      IF( nretfil /= 0 ) THEN
        DO isrc=1,nsrcret
          IF(iuseret(isrc,ipass) > 0) THEN
            retsw=1
            IF (myproc == 0) WRITE(6,'(7x,a,9x,a)') ' Using ',srcret(isrc)
          END IF
        END DO
      END IF
      IF(retsw == 0 .AND. myproc == 0) WRITE(6,'(9x,a)') 'none'

      CALL mpmaxi(sngsw)    ! All processors must have the same flags
      CALL mpmaxi(uasw)
      CALL mpmaxi(radsw)
!
!=======================================================================
!
      20    CONTINUE        ! iteration begin below

!-----------------------------------------------------------------------
!
!      OBTAIN THE COST FUNCTION AND THE GRADIENT.
!
!-----------------------------------------------------------------------
!
      IF (myproc == 0 )  WRITE(ount,'(1x,a,I3,a,/,1x,3(a,I2))')         &
               '=========== iteration - ',icall, ' ==========',         &
               'sngsw = ',sngsw,', uasw = ',uasw,', radsw = ',radsw

      !if (icall > 2) call arpsstop(" ",1)
      cpuin = f_cputime()

      CALL costf(ount,numctr,ctrv,cfun_single,                          &
            gdu_err,gdv_err,gdp_err,gdt_err,gdq_err,gdw_err,            &
            u_ctr,v_ctr,p_ctr,t_ctr,q_ctr,w_ctr, psi, phi,              &
            gdscal, nx,ny,nz,                                           &
            nvar,nvarradin,nvarrad,nzua,nzrdr,nzret,                    &
            mapfct,j1,j2,j3,aj3x,aj3y,aj3z,j3inv,rhostr,                &
            rhostru, rhostrv, rhostrw, div3,                            &
            mxsng,mxua,mxrad,mxcolrad,mxret,mxcolret,                   &
            nsrcsng,nsrcua,nsrcrad,nsrcret,ncat,                        &
            mxpass,ipass,iwstat,xs,ys,zs,x,y,z,zp,hterain,              &
            icatg,xcor,nam_var,                                         &
            ownsng,usesng,xsng,ysng,hgtsng,thesng,                      &
            obsng,odifsng,qobsng,qualsng,isrcsng,icatsng,nobsng,        &
            ownua,useua,xua,yua,hgtua,theua,                            &
            obsua,odifua,qobsua,qualua,isrcua,nlevsua,nobsua,           &
            ownrad,userad,elvrad,xradc,yradc,                           &
            distrad,uazmrad,vazmrad,hgtradc,theradc,dsdr,dhdr,          &
            obsrad,odifrad,qobsrad,qualrad,                             &
            irad,isrcrad,nlevrad,ncolrad,                               &
            xretc,yretc,hgtretc,theretc,                                &
            obsret,odifret,qobsret,qualret,                             &
            iret,isrcret,nlevret,ncolret,                               &
            srcsng,srcua,srcrad,srcret,                                 &
            ianxtyp,iusesng(0,ipass),iuseua(0,ipass),                   &
            iuserad(0,ipass),iuseret(0,ipass),                          &
            xyrange,kpvrsq,wlim,zrange,zwlim,                           &
            thrng,rngsqi,knt,wgtsum,zsum,                               &
            corsng,corua,corrad,corret,                                 &
            xsng_p,ysng_p,ihgtsng,xua_p,yua_p,ihgtua,                   &
            xradc_p,yradc_p,ihgtradc,zsng_1,zsng_2,                     &
            zua_1,zua_2,zradc_1,zradc_2,                                &
            oanxsng,oanxua,oanxrad,oanxret,                             &
            sngsw,uasw,radsw,retsw,                                     &
            ipass_filt(ipass),hradius(ipass),radius_z(1,1,1),          &
            div_opt,cntl_var,smth_flag,                                 &
            wgt_div_h(ipass),wgt_div_v(ipass),wgt_smth(ipass),          &
            thermo_opt,wgt_thermo(ipass),sinlat,ffu,ffv,                &
            anx,tem1,tem2,tem3,tem4,smcu,smcv,smcw,                     &
            tem5,tem6,tem7,tem8,istatus)

      cpuout=f_cputime()
      cputimearray(1)=cputimearray(1)+(cpuout-cpuin)

      cfun_total=0
      DO i=1,nvar+1
        cfun_total=cfun_total+cfun_single(i)
      END DO

!-----------------------------------------------------------------------
!
!  MPI messages
!
!-----------------------------------------------------------------------

      CALL mpsumdp(cfun_total,1)

      IF (myproc == 0) WRITE(ount,'(/,1x,a,F20.2,/)') 'After COSTF, cfun_total = ',cfun_total

      cfun=cfun_total

!-----------------------------------------------------------------------
!
!  Calculate the gradient of the cost function
!
!-----------------------------------------------------------------------
      cpuin=f_cputime()

      CALL gradt(numctr,ctrv,grad,                                      &
            gdu_err,gdv_err,gdp_err,gdt_err,gdq_err,gdw_err,            &
            u_ctr,v_ctr,p_ctr,t_ctr,q_ctr,w_ctr, psi, phi,              &
            gdscal, nx,ny,nz,                                           &
            nvar,nvarradin,nvarrad,nzua,nzrdr,nzret,                    &
            mapfct,j1,j2,j3,aj3x,aj3y,aj3z,j3inv,rhostr,                &
            rhostru, rhostrv, rhostrw, div3,                            &
            mxsng,mxua,mxrad,mxcolrad,mxret,mxcolret,                   &
            nsrcsng,nsrcua,nsrcrad,nsrcret,ncat,                        &
            mxpass,ipass,iwstat,xs,ys,zs,x,y,z,zp,hterain,              &
            icatg,xcor,nam_var,                                         &
            xsng,ysng,hgtsng,thesng,                                    &
            obsng,odifsng,qobsng,qualsng,isrcsng,icatsng,nobsng,        &
            xua,yua,hgtua,theua,                                        &
            obsua,odifua,qobsua,qualua,isrcua,nlevsua,nobsua,           &
            elvrad,xradc,yradc,                                         &
            distrad,uazmrad,vazmrad,hgtradc,theradc,dsdr,dhdr,          &
            obsrad,odifrad,qobsrad,qualrad,                             &
            irad,isrcrad,nlevrad,ncolrad,                               &
            xretc,yretc,hgtretc,theretc,                                &
            obsret,odifret,qobsret,qualret,                             &
            usesng,useua,userad,                                        &
            iret,isrcret,nlevret,ncolret,                               &
            srcsng,srcua,srcrad,srcret,                                 &
            ianxtyp,iusesng(0,ipass),iuseua(0,ipass),                   &
            iuserad(0,ipass),iuseret(0,ipass),                          &
            xyrange,kpvrsq,wlim,zrange,zwlim,                           &
            thrng,rngsqi,knt,wgtsum,zsum,                               &
            corsng,corua,corrad,corret,                                 &
            xsng_p,ysng_p,ihgtsng,xua_p,yua_p,ihgtua,                   &
            xradc_p,yradc_p,ihgtradc,zsng_1,zsng_2,                     &
            zua_1,zua_2,zradc_1,zradc_2,                                &
            oanxsng,oanxua,oanxrad,oanxret,                             &
            sngsw,uasw,radsw,retsw,                                     &
            ipass_filt(ipass),hradius(ipass),radius_z(1,1,1),          &
            div_opt,cntl_var,smth_flag,                                 &
            wgt_div_h(ipass),wgt_div_v(ipass),wgt_smth(ipass),          &
            thermo_opt,wgt_thermo(ipass),sinlat,ffu,ffv,                &
            anx,tem1,tem2,tem3,tem4,smcu,smcv,smcw,                     &
            tem5,tem6,tem7,tem8,tem9,tem10,istatus)

      cpuout=f_cputime()
      cputimearray(2)=cputimearray(2)+(cpuout-cpuin)
!
!-----------------------------------------------------------------------
!
! Call the LBFGS minimization algorithm.
!
!-----------------------------------------------------------------------
!
      cpuin=f_cputime()

      CALL va15ad(ount,numctr,mgra,ctrv,cfun,grad,diagco,diag,iprint,   &
                  eps,swork,ywork,point,work,iflag,ftol)

      cpuout=f_cputime()
      cputimearray(3)=cputimearray(3)+(cpuout-cpuin)

      IF(iflag <= 0) GO TO 50
      icall=icall + 1
      IF(icall > maxin(ipass) ) GO TO 50
      GO TO 20
      50  CONTINUE

!    END DO  ! end of do loop over ipass

    IF (myproc == 0) PRINT *, '-----------------after minimization!---------------------'
!
!-----------------------------------------------------------------------
!
! allocate the arrays for storing the optimal perturbation
!
!
    DO k = 1,nz
      DO j = 1,ny
        DO i = 1,nx
          u_ctr(i,j,k) = 0.0
          v_ctr(i,j,k) = 0.0
            psi(i,j,k) = 0.0
            phi(i,j,k) = 0.0
          p_ctr(i,j,k) = 0.0
          t_ctr(i,j,k) = 0.0
          q_ctr(i,j,k) = 0.0
          w_ctr(i,j,k) = 0.0
           psi(i,j,k)  = 0.0
           phi(i,j,k)  = 0.0
        END DO
      END DO
    END DO

    CALL adtrans(numctr,nx,ny,nz, psi,phi, p_ctr,t_ctr,q_ctr,w_ctr,ctrv,tem4)
!
! ----------------------------------------------------------------------
!
!   U_CTR (psi) & V_CTR (phi) are already on scalar grid
!
!-----------------------------------------------------------------------
!

    tem7=0.0; tem8=0.0; tem9=0.0; tem4=0.0
    CALL ctr_to_vbl(ipass_filt(ipass),hradius(ipass),radius_z(1,1,1),        &
                    nx,ny,nz,gdu_err, gdscal,  psi, tem7,tem8,tem9,tem4)

    tem7=0.0; tem8=0.0; tem9=0.0; tem4=0.0
    CALL ctr_to_vbl(ipass_filt(ipass),hradius(ipass),radius_z(1,1,1),        &
                    nx,ny,nz,gdv_err, gdscal,  phi, tem7,tem8,tem9,tem4)

    tem7=0.0; tem8=0.0; tem9=0.0; tem4=0.0
    CALL ctr_to_vbl(ipass_filt(ipass),hradius(ipass),radius_z(1,1,1),        &
                    nx,ny,nz,gdp_err, gdscal,p_ctr, tem7,tem8,tem9,tem4)

    tem7=0.0; tem8=0.0; tem9=0.0; tem4=0.0
    CALL ctr_to_vbl(ipass_filt(ipass),hradius(ipass),radius_z(1,1,1),        &
                    nx,ny,nz,gdt_err, gdscal,t_ctr, tem7,tem8,tem9,tem4)

    tem7=0.0; tem8=0.0; tem9=0.0; tem4=0.0
    CALL ctr_to_vbl(ipass_filt(ipass),hradius(ipass),radius_z(1,1,1),        &
                    nx,ny,nz,gdq_err, gdscal,q_ctr, tem7,tem8,tem9,tem4)

    tem7=0.0; tem8=0.0; tem9=0.0; tem4=0.0
    CALL ctr_to_vbl(ipass_filt(ipass),hradius(ipass),radius_z(1,1,1),        &
                    nx,ny,nz,gdw_err, gdscal,w_ctr,  tem7,tem8,tem9,tem4)

!
!-----------------------------------------------------------------------
!
! Option cntl_var = 0  U,V as control variables
! Option cntl_var = 1  psi, phi as control variable
!
!-----------------------------------------------------------------------
!
    IF(cntl_var  == 0 ) THEN

      DO k = 1, nz-1
        DO j = 1, ny-1
          DO i = 1, nx-1
            u_ctr(i,j,k) =  psi(i,j,k)
            v_ctr(i,j,k) =  phi(i,j,k)
          END DO
        END DO
      END DO

    ELSE

      DO k = 1, nz-1
        DO j = 2, ny-1
          DO i = 2, nx-1
            u_ctr(i,j,k) = ( psi(i-1,j+1,k)+psi(i,j+1,k)                &
                            -psi(i-1,j-1,k)-psi(i,j-1,k) )/dy/4.        &
                         + ( phi(i,  j,  k)-phi(i-1,j,k) )/dx
            u_ctr(i,j,k) = u_ctr(i,j,k)*mapfct(i,j,2)

            v_ctr(i,j,k) = ( psi(i+1,j-1,k)+psi(i+1,j,k)                &
                            -psi(i-1,j-1,k)-psi(i-1,j,k) )/dx/4.        &
                         + ( phi(i,  j,  k)-phi(i,j-1,k) )/dy
            v_ctr(i,j,k) = v_ctr(i,j,k)*mapfct(i,j,3)
          END DO
        END DO

        DO j=2,ny-1
          u_ctr( 1,j,k)=u_ctr( 2,j,k)+u_ctr( 2,j,k)-u_ctr( 3,j,k)
          v_ctr( 1,j,k)=v_ctr( 2,j,k)+v_ctr( 2,j,k)-v_ctr( 3,j,k)
          u_ctr(nx,j,k)=u_ctr(nx-1,j,k)+u_ctr(nx-1,j,k)-u_ctr(nx-2,j,k) ! Not necessary since analysis is
          v_ctr(nx,j,k)=v_ctr(nx-1,j,k)+v_ctr(nx-1,j,k)-v_ctr(nx-2,j,k) ! on SCALAR grid
        END DO

        DO i=2,nx-1
          u_ctr(i, 1,k)=u_ctr(i, 2,k)+u_ctr(i, 2,k)-u_ctr(i, 3,k)
          v_ctr(i, 1,k)=v_ctr(i, 2,k)+v_ctr(i, 2,k)-v_ctr(i, 3,k)
          u_ctr(i,ny,k)=u_ctr(i,ny-1,k)+u_ctr(i,ny-1,k)-u_ctr(i,ny-2,k) ! Not necessary since analysis is
          v_ctr(i,ny,k)=v_ctr(i,ny-1,k)+v_ctr(i,ny-1,k)-v_ctr(i,ny-2,k) ! on SCALAR grid
        END DO

        u_ctr(1,1 ,k)  = 0.5*( u_ctr(2,1,k)+u_ctr(1,2,k) )
        v_ctr(1,1 ,k)  = 0.5*( v_ctr(2,1,k)+v_ctr(1,2,k) )
        u_ctr(1,ny,k)  = 0.5*( u_ctr(1,ny-1,k)+u_ctr(2,ny,k) )      ! Not necessary since analysis is
        v_ctr(1,ny,k)  = 0.5*( v_ctr(1,ny-1,k)+v_ctr(2,ny,k) )      ! on SCALAR grid

        u_ctr(nx,1,k)  = 0.5*( u_ctr(nx-1,1,k)+u_ctr(nx,2,k) )      ! Not necessary since analysis is
        v_ctr(nx,1,k)  = 0.5*( v_ctr(nx-1,1,k)+v_ctr(nx,2,k) )      ! on SCALAR grid
        u_ctr(nx,ny,k) = 0.5*( u_ctr(nx,ny-1,k)+u_ctr(nx-1,ny,k) )
        v_ctr(nx,ny,k) = 0.5*( v_ctr(nx,ny-1,k)+v_ctr(nx-1,ny,k) )

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

!  call smooth ( u_ctr(1,1,1), NX, NY, NZ, NUM_SMTH, 1)
!  call smooth ( v_ctr(1,1,1), NX, NY, NZ, NUM_SMTH, 1)
!  call smooth ( p_ctr(1,1,1), NX, NY, NZ, NUM_SMTH, 0)
!  call smooth ( t_ctr(1,1,1), NX, NY, NZ, NUM_SMTH, 0)
!  call smooth ( q_ctr(1,1,1), NX, NY, NZ, NUM_SMTH, 0)
!  call smooth ( w_ctr(1,1,1), NX, NY, NZ, NUM_SMTH, 0)

!
!-----------------------------------------------------------------------
!
!  add the optimal perturbation to background
!
!-----------------------------------------------------------------------
!
!
    DO k = 1, nz
      DO j = 1, ny
        DO i = 1, nx

          anx(i,j,k,1) = anx(i,j,k,1) + u_ctr(i,j,k)
          anx(i,j,k,2) = anx(i,j,k,2) + v_ctr(i,j,k)
          anx(i,j,k,3) = anx(i,j,k,3) + p_ctr(i,j,k)
          anx(i,j,k,4) = anx(i,j,k,4) + t_ctr(i,j,k)
          anx(i,j,k,5) = anx(i,j,k,5) + q_ctr(i,j,k)
          anx(i,j,k,6) = anx(i,j,k,6) + w_ctr(i,j,k)

          IF ( anx(i,j,k,5) <= 0.) anx(i,j,k,5) = 0.0

        END DO
      END DO
    END DO

  END IF  ! assim_opt == 1

!----------------------------------------------------------------------
!
! Clear up before Return
!
!----------------------------------------------------------------------

  istatus=0

  !DEALLOCATE( smcu, smcv, smcw )
  !DEALLOCATE( sinlat, ffu, ffv, tem2d )

  DEALLOCATE( xsng_p,  ysng_p,  zsng_1,  zsng_2,  ihgtsng )
  DEALLOCATE( xua_p,   yua_p,   zua_1,   zua_2,   ihgtua )
  DEALLOCATE( xradc_p, yradc_p, zradc_1, zradc_2, ihgtradc )

  RETURN
END SUBROUTINE minimization

!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE LINEAR                    ######
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
!
!  Linear interplation in one direction.
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!
!  Jidong Gao, CAPS, July, 2000
!
!-----------------------------------------------------------------------
!
SUBROUTINE linear(n,xx,yy,u,f)

  USE arps_precision

  INTEGER :: n
  REAL(P) :: xx(n),yy(n),u,f

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF( u <= xx(1)) THEN

    a1=( u-xx(2))/(xx(1)-xx(2))
    a2=( u-xx(1))/(xx(2)-xx(1))
    f=a1*yy(1)+a2*yy(2)

  ELSE

    DO i=2,n
      IF( u <= xx(i) )  EXIT
    END DO
    i = min(i,n)

    a1=( u-xx(i)  )/(xx(i-1)-xx(i))
    a2=( u-xx(i-1))/(xx(i)-xx(i-1))
    f=a1*yy(i-1)+a2*yy(i)

  END IF

  RETURN
END SUBROUTINE linear
