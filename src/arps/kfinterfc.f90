!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE KFINTERFC                  ######
!######                                                      ######
!######                     Developed by                     ######
!######          Coastal Meteorology Research Program        ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE kfinterfc(nx,ny,nz,u,v,w,pprt,ptprt,qv,pbar,ptbar,zp,        &
           ptcumsrc,qcumsrc,rainc,prcrate,nca,raincv,                   &
           ptop,psfc,psb,                                               &
           sigma,hfsig,ub,vb,w0,tb,qvb,tem1,tem2)

!
!-----------------------------------------------------------------------
!   PURPOSE:
!
!    To calculate the necessary inputs for the Kain-Fritsch
!    cumulus parameterization scheme. Note that the inputs
!    to the k-f scheme are arranged in such a way that the
!    k index is counted from top down, i.e., k=1 is top and
!    k=kx is the bottom.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Zonghui Huo
!  08/01/97
!
!  MODIFICATION HISTORY:
!
!  1/6/98 (Z. Huo and J. Kain)
!  Fixed a bug in the calculation of hfsig. Optimize the code the
!  call to kfpara.f.
!
!  2/3/98 Zonghui Huo
!  A bug fix in loop 563. tem2(i,j,k) = 0.5*(tem1(i,j,k)-tem1(i,j,k-1))
!  should have been       tem2(i,j,k) = 0.5*(tem1(i,j,k)+tem1(i,j,k-1))
!
!  2/9/1998 (M. Xue)
!  Made a few additional optimization.
!
!  4/15/1998 (Donghai Wang)
!  Added the source terms for qc,qr,qi and qs equations.
!  Added the running average vertical velocity instead of the 2-level
!  time average.
!  Passed two options mphyopt and kffbfct from arps.input to K-F scheme.
!
!  4/18/2002 (Zuwen He)
!  Passed sub-saturation option kfsubsattrig from arps.input to K-F scheme.
!
!----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    u        x component of velocity at a given time level (m/s)
!    v        y component of velocity at a given time level (m/s)
!    w        time-average vertical velocity (m/s)
!    pprt     Perturbation pressure (Pascal)
!    ptprt    Perturbation potential temperature (K)
!    qv       Water vapor specific humidity (kg/kg)
!    pbar     Base state pressure (Pascal)
!    ptbar    Base state potential temperature (K)
!    zp       The physical height coordinate defined at w-point
!
!  OUTPUT:
!
!    ptcumsrc Potential temperature source term.
!    qcumsrc Water specific humidity source term
!    rainc    Accumulated precipitation by cumulus convection (mm)
!    prcrate  precipitatioon rate (mm/s)
!    raincv   K-F convective rainfall (cm)
!    nca      K-F counter for CAPE release
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  INTEGER :: nx,ny,nz

  REAL :: u      (nx,ny,nz) ! Total u-velocity (m/s)
  REAL :: v      (nx,ny,nz) ! Total v-velocity (m/s)
  REAL :: w      (nx,ny,nz) ! Total w-velocity (m/s)

  REAL :: pprt   (nx,ny,nz) ! Perturbation pressure (Pascal)
  REAL :: ptprt  (nx,ny,nz) ! Perturbation potential temperature (K)
  REAL :: qv     (nx,ny,nz) ! Water vapor specific humidity (kg/kg)
  REAL :: pbar   (nx,ny,nz) ! Base state pressure (Pascal)
  REAL :: ptbar  (nx,ny,nz) ! Base state potential temperature (K)

  REAL :: zp     (nx,ny,nz) ! The physical height coordinate defined
                            ! at w-point
!
! Output array
!

  REAL :: ptcumsrc(nx,ny,nz) ! Source term in pt equation.
  REAL :: qcumsrc(nx,ny,nz,5)  ! Source term in water equations due
                               ! to cumulus parameterization:
                               ! qcumsrc(1,1,1,1) for qv equation
                               ! qcumsrc(1,1,1,2) for qc equation
                               ! qcumsrc(1,1,1,3) for qr equation
                               ! qcumsrc(1,1,1,4) for qi equation
                               ! qcumsrc(1,1,1,5) for qs equation
  REAL :: rainc(nx,ny)       ! Accumulated precipitation by convec(mm)
  REAL :: prcrate(nx,ny)     ! precipitation rate (kg/(m**2*s))
  REAL :: raincv  (nx,ny)    ! K-F output (cm)
  INTEGER :: nca     (nx,ny)    ! counter for CAPE release

!
! Local array
!

  REAL :: ptop (nx,ny)    ! Pressure at the model top (cb)
  REAL :: psfc (nx,ny)    ! Surface pressure (cb=1000Pa)
  REAL :: psb  (nx,ny)    ! Psfc - Ptop, used by mm5
  REAL :: sigma(nx,ny,nz) ! sigma levels where w is defined
  REAL :: hfsig(nx,ny,nz) ! half sigma levels where all variables defined
  REAL :: ub   (nx,ny,nz) ! u at the mass points (times pstr)
  REAL :: vb   (nx,ny,nz) ! v at the mass points
  REAL :: w0   (nx,ny,nz) ! w at the mass points
  REAL :: tb   (nx,ny,nz) ! Temperature (K) at the mass points
  REAL :: qvb  (nx,ny,nz) ! Water vapor specific humidity (kg/kg) at
                          ! mass points of the physical domain
!
! More local array related to the optimization to call the kfpatra.f
!
  INTEGER :: ncuyes,k300,kmax,kmin,nzmax,nxmax
  REAL :: rovg,psfck,p300,ttmp,qtmp,es,qs,semmx,semmin
  REAL :: cp,xlv,g,aliq,bliq,cliq,dliq
  PARAMETER(nzmax=100,nxmax=500)
  PARAMETER(cp=1005.7,xlv=2.5E6,g=9.8,rovg=29.29)
  PARAMETER(aliq=613.3,bliq=17.502,cliq=4780.8,dliq=32.19)
  INTEGER :: icuyes(nxmax)
  INTEGER :: lsb   (nxmax)
  REAL :: sem   (nzmax)
  REAL :: sems  (nzmax)
  REAL :: prs   (nzmax)
  REAL :: tv0   (nzmax)
  REAL :: z0    (nzmax)

  REAL :: tem1   (nx,ny,nz) ! temporary
  REAL :: tem2   (nx,ny,nz) ! temporary

!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,ix,jx,kx,nk
  REAL :: delz,tv, dtbig10
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Set the loop index for the physical space mass point since
!  the K-F scheme is only performed in physical space on mass
!  points
!
!-----------------------------------------------------------------------
!

  ix=nx-1
  jx=ny-1
  kx=nz-3
!
!-----------------------------------------------------------------------
!
! Calculate the pressure at model top (ptop) and model bottom (psfc)
! The temprature (K) is writen in tv.
!
!-----------------------------------------------------------------------
!

! Surface pressure

  DO j=1,ny-1
    DO i=1,nx-1
      delz = (zp(i,j,3)-zp(i,j,2))/2.0
      tv = (ptbar(i,j,2)+ptprt(i,j,2))                                  &
          * (1.0E-5*(pbar(i,j,2)+pprt(i,j,2)))**0.286
      psfc(i,j) = (1.0E-3)*(pbar(i,j,2)+pprt(i,j,2))                    &
                  *EXP(9.8*delz/(287.0*tv))

    END DO
  END DO

! Model top pressure

  DO j=1,ny-1
    DO i=1,nx-1
      delz = (zp(i,j,nz-1)-zp(i,j,nz-2))/2.0
      tv = (ptbar(i,j,nz-2)+ptprt(i,j,nz-2))                            &
          * (1.0E-5*(pbar(i,j,nz-2)+pprt(i,j,nz-2)))**0.286
      ptop(i,j) = (1.0E-3)*(pbar(i,j,nz-2)+pprt(i,j,nz-2))              &
                  *EXP(-9.8*delz/(287.0*tv))

    END DO
  END DO

  DO j=1,jx
    DO i=1,ix
      psb(i,j) = psfc(i,j) -  ptop(i,j)
    END DO
  END DO

!
! Calculate the sigma levels corresponding to the Gal-Chen levels
! and the half sigma levels where all other variables are defined
!

  DO k=1,kx
    DO j=1,jx
      DO i=1,ix
        tem1(i,j,k)=((pbar(i,j,k+1)+pprt(i,j,k+1))*(1.0E-3)             &
                  -ptop(i,j))/(psfc(i,j)-ptop(i,j))
      END DO
    END DO
  END DO

  DO j=1,jx
    DO i=1,ix
      tem2(i,j,1) = 1.0
      tem2(i,j,kx+1) = 0.0
    END DO
  END DO

  DO k=2,kx
    DO j=1,jx
      DO i=1,ix
        tem2(i,j,k) = 0.5*(tem1(i,j,k)+tem1(i,j,k-1))
      END DO
    END DO
  END DO

  DO k=1,kx
    nk=kx-k+1
    DO j=1,jx
      DO i=1,ix
        hfsig(i,j,nk) = tem1(i,j,k)
      END DO
    END DO
  END DO

  DO k=1,kx+1
    nk=kx-k+2
    DO j=1,jx
      DO i=1,ix
        sigma(i,j,nk) = tem2(i,j,k)
      END DO
    END DO
  END DO

!
! Get u and v winds on mass points (Arakawa c grid in ARPS)
!

  DO k=1,kx
    nk=kx-k+1
    DO j=1,jx
      DO i=1,ix
        ub(i,j,nk) =  (u(i,j,k+1)+u(i+1,j,k+1))*0.5
        vb(i,j,nk) =  (v(i,j,k+1)+v(i,j+1,k+1))*0.5
      END DO
    END DO
  END DO

!
! Calculate the w on the mass points (vertical half levels)
!

  DO k=1,kx
    nk=kx-k+1
    DO j=1,jx
      DO i=1,ix
        w0(i,j,nk) =  (w(i,j,k+1)+w(i,j,k+2))*0.5
      END DO
    END DO
  END DO

!
! Calculate tb (in K) and put qvb on mass levels, Rd/Cp=2858566
!

  DO k=1,kx
    nk=kx-k+1
    DO j=1,jx
      DO i=1,ix
        tb(i,j,nk) = (ptbar(i,j,k+1)+ptprt(i,j,k+1))*                   &
                     ((pbar(i,j,k+1) +pprt(i,j,k+1))*                   &
                      1.0E-5)**0.2858566
      END DO
    END DO
  END DO

  DO k=1,kx
    nk=kx-k+1
    DO j=1,jx
      DO i=1,ix
        qvb(i,j,nk) =  qv(i,j,k+1)
      END DO
    END DO
  END DO


!-----------------------------------------------------------------------
!
  DO k=nz-3,1,-1
    nk=nz-k
    DO j=1,ny-1
      DO i=1,nx-1
        tem1(i,j,k) = qcumsrc(i,j,nk-1,2)
        tem2(i,j,k) = qcumsrc(i,j,nk-1,3)
      END DO
    END DO
  END DO

  DO k=nz-3,1,-1
    DO j=1,ny-1
      DO i=1,nx-1
        qcumsrc(i,j,k,2) = tem1(i,j,k)
        qcumsrc(i,j,k,3) = tem2(i,j,k)
      END DO
    END DO
  END DO

  DO k=nz-3,1,-1
    nk=nz-k
    DO j=1,ny-1
      DO i=1,nx-1
        tem1(i,j,k) = qcumsrc(i,j,nk-1,4)
        tem2(i,j,k) = qcumsrc(i,j,nk-1,5)
      END DO
    END DO
  END DO

  DO k=nz-3,1,-1
    DO j=1,ny-1
      DO i=1,nx-1
        qcumsrc(i,j,k,4) = tem1(i,j,k)
        qcumsrc(i,j,k,5) = tem2(i,j,k)
      END DO
    END DO
  END DO


!
!-----------------------------------------------------------------------
!
! Call the kain-Fritsch scheme. It is checked by east-west i slice.
!
! Initialize temperary arrays, which are used to hold dt/dt and dq/dt.
!
!-----------------------------------------------------------------------
!
  DO k=1,nz
    DO j=1,ny
      DO i=1,nx
        tem1(i,j,k) = 0.0
        tem2(i,j,k) = 0.0
      END DO
    END DO
  END DO

!
!-----------------------------------------------------------------------
!...make a quick for check for these things that can eliminate
!...grid points from the possibility of convective initiation:
!
!     1.) convection already active at this point
!     2.) point close to grid-domain boundary
!     3.) downward motion at all levels in lowest 300 mb
!     4.) no CAPE
!
!-----------------------------------------------------------------------
!
 !!!! Peter, the variable lists in the multitasking directives below
 !!!! (in SGI compiler format) should be looked at carefully since they
 !!!! were used with earlier versions of the code.  But, the idea that
 !!!! they put into effect, which is one of sending the KF call for individual
 !!!! j-slices to different processors, should be fairly easy to activate.
!c$doacross
!c$& share(ix,jx,kx,w0,nca,psfc,a,psb,ptop,
!c$&       tb,qvb,rovg),
!c$& local(i,j,k,nk,ncuyes,icuyes,psfck,p300,k300,prs,ttmp,qtmp,
!c$&       tv0,z0,es,qs,sem,sems,semmx,semmin,kmax,kmin)
!
!
!-----------------------------------------------------------------------
!
  DO j=1,jx
    ncuyes=0
    DO i = 1,ix
      icuyes(i)=0
      lsb(i) = 0
    END DO
    DO i = 1,ix
!
!...if parameterization already active at this point, go on to next point...
!
      IF (nca(i,j) > 0) CYCLE
      psfck = psfc(i,j)*1.e3
!
!...find highest model level k-value in lowest 300 mb...
!
      p300 = psfck-3.e4
      k300 = 1
      DO k = 2,kx
        nk = kx-k+1
        prs(k)=1.e3*(hfsig(i,j,nk)*psfc(i,j)+ptop(i,j))
        IF(prs(k) < p300) EXIT
        k300 = k
      END DO
!      15     CONTINUE
!
!...vertical velocity must be upward at some level in the
!...lowest 300 mb...
!
      DO k = 1,k300
        IF(w0(i,j,k) > 0.) GO TO 25
      END DO
      CYCLE
      25     CONTINUE
!
!...calculate moist static energy and saturation moist static energy
!...at each vertical level...
!
      DO k = kx,1,-1
        ttmp = tb(i,j,k)
        qtmp = qvb(i,j,k)
        prs(k)=1.e3*(hfsig(i,j,k)*psfc(i,j)+ptop(i,j))
        tv0(k) = ttmp*(1.+.608*qtmp)
        IF(k == kx) THEN
          z0(k) = 0.
        ELSE
          z0(k) = z0(k+1)-rovg*0.5*(tv0(k)+tv0(k+1))*                   &
                ALOG(prs(k)/prs(k+1))
        END IF
        es=aliq*EXP((bliq*ttmp-cliq)/(ttmp-dliq))
        qs=0.622*es/(prs(k)-es)
        sems(k)=cp*ttmp+xlv*qs+g*z0(k)
        sem(k)=cp*ttmp+xlv*qtmp+g*z0(k)
      END DO
!
!...determine level of maximum MSE in lowest 300 mb...
!
      semmx = 0.
      DO nk =1,k300
        k = kx-nk+1
        semmx=AMAX1(semmx,sem(k))
        IF(semmx == sem(k)) kmax=k
      END DO
!
!...determine level of minimum saturated MSE above kmax...
!
      semmin = 1.e6
      DO k = kmax,1,-1
        semmin=AMIN1(semmin,sems(k))
        IF(semmin == sems(k)) kmin=k
      END DO
      IF(semmx > semmin)THEN
        ncuyes=ncuyes+1
        icuyes(ncuyes)=i
        lsb(i) = kx-kmax+1
      END IF
    END DO
    IF(ncuyes > 0)THEN

      CALL kfpara(nx,ny,nz,j,dtbig,dx,mphyopt,kffbfct,                  &
                  kfsubsattrig,                                         & 
                  ptop,psb,sigma,hfsig,ub,vb,w0,tb,qvb,                 &
                  tem1,tem2,qcumsrc(1,1,1,2),qcumsrc(1,1,1,3),          &
                  qcumsrc(1,1,1,4),qcumsrc(1,1,1,5),                    &
                  nca,raincv,ncuyes,icuyes,lsb)
    END IF
  END DO

!
!-----------------------------------------------------------------------
!
! Put the K-F results into the forcing terms (note the index change,
! and the lateral boundary values, the k-f temperature rate should be
! converted back to potential temperature rate)) Rd/Cp=2858566
!
!-----------------------------------------------------------------------
!
  DO k=2,nz-2
    nk=nz-k
    DO j=1,ny-1
      DO i=1,nx-1
        ptcumsrc(i,j,k) = ptcumsrc(i,j,k)                               &
                        + tem1(i,j,nk-1)*(1.0E5/(pbar(i,j,k)            &
                        + pprt(i,j,k)))**0.2858566
        qcumsrc(i,j,k,1) = qcumsrc(i,j,k,1) + tem2(i,j,nk-1)
      END DO
    END DO
  END DO

  DO k=2,nz-2
    nk=nz-k
    DO j=1,ny-1
      DO i=1,nx-1
        tem1(i,j,k) = qcumsrc(i,j,nk-1,2)
        tem2(i,j,k) = qcumsrc(i,j,nk-1,3)
      END DO
    END DO
  END DO

  DO k=2,nz-2
    DO j=1,ny-1
      DO i=1,nx-1
        qcumsrc(i,j,k,2) = tem1(i,j,k)
        qcumsrc(i,j,k,3) = tem2(i,j,k)
      END DO
    END DO
  END DO

  DO k=2,nz-2
    nk=nz-k
    DO j=1,ny-1
      DO i=1,nx-1
        tem1(i,j,k) = qcumsrc(i,j,nk-1,4)
        tem2(i,j,k) = qcumsrc(i,j,nk-1,5)
      END DO
    END DO
  END DO

  DO k=2,nz-2
    DO j=1,ny-1
      DO i=1,nx-1
        qcumsrc(i,j,k,4) = tem1(i,j,k)
        qcumsrc(i,j,k,5) = tem2(i,j,k)
      END DO
    END DO
  END DO

  dtbig10 = 10.0/dtbig

  DO j=1,ny-1
    DO i=1,nx-1
      prcrate(i,j) = raincv(i,j)*dtbig10
    END DO
  END DO
  RETURN
END SUBROUTINE kfinterfc
