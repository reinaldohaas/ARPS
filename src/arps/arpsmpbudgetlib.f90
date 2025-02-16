
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
  DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,               &
       24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,  &
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

SUBROUTINE cal_Nt(nx,ny,nz,rhoa,q,N0,cx,alpha,Ntx)

!   
!-----------------------------------------------------------------------
!  PURPOSE:  Calculates number concentration at scalar points
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Dawson 
!  (02/06/2008) 
!   
!  MODIFICATION HISTORY:
!  
!  03/31/08 - converted intermediate calculations to double precision
!             as well as a few of the input arguments.
! 
!-----------------------------------------------------------------------
!  Variable Declarations:
!-----------------------------------------------------------------------
!     

  INTEGER :: nx,ny,nz
  REAL :: rhoa(nx,ny,nz),q(nx,ny,nz)
  REAL*8 :: N0(nx,ny,nz),alpha(nx,ny,nz)
  REAL :: cx
  REAL :: Ntx(nx,ny,nz)
  REAL*8 :: gamma1,gamma4

  REAL*8 :: gamma

  INTEGER i,j,k

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        gamma1 = gamma(1.d0+alpha(i,j,k))
        gamma4 = gamma(4.d0+alpha(i,j,k))

       !print*,'gamma1,gamma4,cx',gamma1,gamma4,cx

       Ntx(i,j,k) = sngl((N0(i,j,k)*gamma1)**(3.d0/(4.d0+alpha(i,j,k)))*   &
                         ((gamma1/gamma4)*dble(rhoa(i,j,k))* &
                         dble(q(i,j,k))/dble(cx))**((1.d0+alpha(i,j,k))/(4.d0+alpha(i,j,k))))
      END DO
    END DO
  END DO

END SUBROUTINE cal_Nt

SUBROUTINE cal_lamda(nx,ny,nz,rhoa,q,Ntx,cx,alpha,lamda)

!
!-----------------------------------------------------------------------
!  PURPOSE:  Calculates slope parameter lamda
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Dawson
!  (02/06/2008)
!
!  MODIFICATION HISTORY:
!  (03/31/2008)
!  Converted intermediate calculations and arrays alpha and lamda to
!  double precision. 
!-----------------------------------------------------------------------
!  Variable Declarations:
!-----------------------------------------------------------------------
!

  INTEGER :: nx,ny,nz
  REAL :: rhoa(nx,ny,nz),q(nx,ny,nz)
  REAL*8 :: lamda(nx,ny,nz),alpha(nx,ny,nz)
  REAL :: cx
  REAL :: Ntx(nx,ny,nz)
  REAL*8 :: gamma1, gamma4

  REAL*8 :: gamma


  INTEGER i,j,k

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1

        gamma1 = gamma(1.d0+alpha(i,j,k))
        gamma4 = gamma(4.d0+alpha(i,j,k))

        IF(rhoa(i,j,k) > 0.0 .and. q(i,j,k) > 0.0) THEN
          lamda(i,j,k) = ((gamma4/gamma1)*dble(cx)*dble(Ntx(i,j,k))/(dble(rhoa(i,j,k))*  &
                dble(q(i,j,k))))**(1.d0/3.d0)
        ELSE
          lamda(i,j,k) = 0.d0
        END IF

      END DO
    END DO
  END DO

END SUBROUTINE cal_lamda


SUBROUTINE cal_N0(nx,ny,nz,rhoa,q,Ntx,cx,alpha,N0,N0_eff)

!   
!-----------------------------------------------------------------------
!  PURPOSE:  Calculates intercept parameter and "effective" intercept parameter
!            
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Dawson 
!  (02/06/2008) 
!   
!  MODIFICATION HISTORY:
!   
!  (03/26/2008)
!  Recast N0 as a double precision variable, and used double precision for
!  all intermediate calculations.  The calling subroutine should
!  also define it as double precision.  For situations with large alpha,
!  N0 can become very large, and loss of precision can result.
!  Also tweaked the calculation of N0 a bit to avoid overflow, in keeping
!  With Jason Milbrandt's calculation of N0 just before evaporation in 
!  the multi-moment code.
!
!-----------------------------------------------------------------------
!  Variable Declarations:
!-----------------------------------------------------------------------
!     

!  USE MM_FUNCTIONS

  INTEGER :: nx,ny,nz
  REAL :: rhoa(nx,ny,nz),q(nx,ny,nz),Ntx(nx,ny,nz)
  REAL*8 :: N0(nx,ny,nz),alpha(nx,ny,nz),N0_eff(nx,ny,nz)
  REAL :: cx
  REAL*8 :: gamma1, gamma4

  REAL*8 :: gamma

  DOUBLE PRECISION :: lamda

  INTEGER i,j,k

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1

        gamma1 = gamma(1.d0+alpha(i,j,k))
        gamma4 = gamma(4.d0+alpha(i,j,k))

        IF(rhoa(i,j,k) > 0.0 .and. q(i,j,k) > 0.0) THEN
          lamda = ((gamma4/gamma1)*dble(cx)*dble(Ntx(i,j,k))/(dble(rhoa(i,j,k))*  &
              dble(q(i,j,k))))**(1.d0/3.d0)
        ELSE
          lamda = 0.d0
        END IF

        N0(i,j,k) = dble(Ntx(i,j,k))*lamda**(0.5d0*(1.d0+alpha(i,j,k)))*    &
                    (1.d0/gamma1)*lamda**(0.5d0*(1.d0+alpha(i,j,k)))

        IF(lamda /= 0.d0) THEN
          N0_eff(i,j,k) = N0(i,j,k)*(((4.d0+alpha(i,j,k))/lamda)**    &
                           alpha(i,j,k))*gamma4*(128.d0/3.d0)/         &
                           ((4.d0+alpha(i,j,k))**(4.d0+alpha(i,j,k)))
        ELSE
          N0_eff(i,j,k) = 0.d0
        END IF
      END DO
    END DO
  END DO

END SUBROUTINE cal_N0

SUBROUTINE cal_Dm(nx,ny,nz,rhoa,q,Ntx,cx,Dm)
!   
!-----------------------------------------------------------------------
!  PURPOSE:  Calculates mean-mass diameter
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Dawson 
!  (02/06/2008) 
!   
!  MODIFICATION HISTORY:
!   
!-----------------------------------------------------------------------
!  Variable Declarations:
!-----------------------------------------------------------------------
!     

  INTEGER :: nx,ny,nz
  REAL :: rhoa(nx,ny,nz),q(nx,ny,nz)
  REAL :: cx
  REAL :: Ntx(nx,ny,nz),Dm(nx,ny,nz)

  INTEGER i,j,k

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        IF(Ntx(i,j,k) /= 0.0) THEN
          Dm(i,j,k) = (rhoa(i,j,k)*q(i,j,k)/(cx*Ntx(i,j,k)))**(1./3.)
        ELSE
          Dm(i,j,k) = 0.0
        END IF
      END DO
    END DO
  END DO

END SUBROUTINE cal_Dm

SUBROUTINE diag_alpha(nx,ny,nz,rhoa,varid_qscalar,Dm,alpha)
!   
!-----------------------------------------------------------------------
!  PURPOSE:  Calculates shape parameter alpha
!-----------------------------------------------------------------------
!     
!  AUTHOR: Dan Dawson   
!  (02/06/2008)         
!   
!  MODIFICATION HISTORY:
!  (03/31/2008)
!  Changed alpha array to double precision.
!-----------------------------------------------------------------------
!  Variable Declarations:
!-----------------------------------------------------------------------
!     

  INTEGER :: nx,ny,nz
  REAL :: rhoa(nx,ny,nz),Dm(nx,ny,nz)
  REAL*8 :: alpha(nx,ny,nz)
  REAL*8 :: ddm

  REAL*8 :: diagAlpha

  CHARACTER (LEN=2) :: varid_qscalar

  INTEGER i,j,k,nq

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

  print*,'nq',nq

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        ddm = dble(Dm(i,j,k))
        alpha(i,j,k) = diagAlpha(ddm,nq)
      END DO
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

SUBROUTINE solve_alpha(nx,ny,nz,rhoa,cx,q,Ntx,Z,alpha)
!   
!-----------------------------------------------------------------------
!  PURPOSE:  Calculates shape parameter alpha
!-----------------------------------------------------------------------
!     
!  AUTHOR: Dan Dawson   
!  (02/06/2008)         
!   
!  MODIFICATION HISTORY:
!  (03/31/2008)
!  Changed alpha array to double precision
!-----------------------------------------------------------------------
!  Variable Declarations:
!-----------------------------------------------------------------------
!     

  INTEGER :: nx,ny,nz
  REAL :: rhoa(nx,ny,nz),q(nx,ny,nz),Ntx(nx,ny,nz),Z(nx,ny,nz)
  REAL*8 :: alpha(nx,ny,nz)

  REAL*8 :: solveAlpha
  REAL*8 :: dsA

  REAL :: cx

  INTEGER i,j,k

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        IF(q(i,j,k) > 0.0 .and. Ntx(i,j,k) > 0.0 .and. Z(i,j,k) > 0.0) THEN

          alpha(i,j,k) = solveAlpha(q(i,j,k),Ntx(i,j,k),Z(i,j,k),cx,rhoa(i,j,k))

        ELSE
          alpha(i,j,k) = 0.d0
        END IF
      END DO
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
