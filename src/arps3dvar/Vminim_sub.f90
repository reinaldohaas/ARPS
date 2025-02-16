!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE VA15AD                    ######
!######                                                      ######
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
!  minimization code.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!
!  unknown, copied from Florida state University. 
!
!-----------------------------------------------------------------------
!
!
!  ----------------------------------------------------------------------
!
!  CALL VA15AD(Numctr,MGRA,ctrv,CFUN,grad,DIAGCO,DIAG,IPRINT,
!    :                  EPS,SWORK,YWORK,POINT,WORK,IFLAG,FTOL)
!
!

SUBROUTINE va15ad(ount,n,m,x,f,g,diagco,diag,iprint,eps,s,y,            &
                  point,w,iflag,ftol)
  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: ount
  INTEGER, INTENT(IN)    :: n, m
  REAL,    INTENT(IN)    :: x(n),g(n)
  REAL,    INTENT(INOUT) :: diag(n), s(m*n), y(m*n), w(n+2*m)
  REAL,    INTENT(IN)    :: f, eps, ftol
  LOGICAL, INTENT(IN)    :: diagco
  INTEGER, INTENT(IN)    :: iprint(2)
  INTEGER, INTENT(INOUT) :: iflag
  INTEGER, INTENT(OUT)   :: point

  REAL    :: xtol,stpmin,stpmax,stp,ys,sq,yr,beta,xnorm,gnorm,yy,stp1
  INTEGER :: bound,iter,nfun,nfev,cp

  LOGICAL :: finish

  INTEGER :: mp, lp
  REAL    :: gtol
  COMMON /va15dd/mp,lp, gtol

  REAL    :: one, zero
  SAVE
  DATA one,zero/1.0E+0,0.0E+0/

  INTEGER :: i, j, k, ii
  INTEGER :: maxfev, info, npt
!
! Functions
!
  REAL :: ddot

!  INCLUDE 'mp.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!
!  ------------------------------------------------------------
!  INITIALIZE
!  ------------------------------------------------------------
!
  IF(iflag == 0) GO TO 1
  GO TO (72,10) iflag
  1  iter= 0
  IF(n <= 0.OR.m <= 0) GO TO 96
  IF(gtol <= 1.d-04) THEN
    IF(lp > 0) WRITE(lp,145)
    gtol=1.d-02
  END IF
  nfun= 1
  point= 0
  finish= .false.
  IF(diagco) THEN
    DO i=1,n
      IF (diag(i) <= zero) GO TO 95
    END DO
  ELSE
    DO i=1,n
      diag(i)= 1.0D0
    END DO
  END IF
  DO i=1,n
    s(i)= -g(i)*diag(i)
  END DO
  gnorm= SQRT(ddot(n,g,g))
  IF( gnorm > 0.0 ) THEN
    stp1= one/gnorm
  ELSE
    iflag=-13 
    RETURN
  ENDIF
!
!  PARAMETERS FOR LINE SEARCH ROUTINE
!  ----------------------------------
  xtol   = 1.0D-17
  stpmin = 1.0D-20
  stpmax = 1.0D+20
  maxfev = 20

!  WRITE(ount,*)'iprint=',iprint,iter,nfun,n,m,f,stp,finish
  CALL va15bd(iprint,iter,nfun,n,m,x,f,g,stp,finish)
!
!    ------------------------------------------------------------
!  MAIN ITERATION LOOP
!    --------------------------------------------------------
!
  8    iter= iter+1
  info=0
  bound=iter-1
  IF (iter >  m) bound=m
  IF (iter == 1) GO TO 65
!
!  ------------------------------------------------------------
!  COMPUTE -HG AND THE DIAGONAL SCALING MATRIX IN DIAG
!  ------------------------------------------------------------
!
  IF(.NOT.diagco) THEN
    DO i=1,n
      diag(i)= ys/yy
    END DO
  ELSE
    iflag=2
    RETURN
  END IF
  10  CONTINUE
  DO i=1,n
    IF (diag(i) <= zero) GO TO 95
  END DO

  cp= point
  IF (point == 0) cp=m
  w(n+cp)= one/ys
  DO i=1,n
    w(i)= -g(i)
  END DO
  cp= point
  DO ii= 1,bound
    cp=cp-1
    IF (cp == -1)cp=m-1
    sq= ddot(n,s(cp*n+1),w)
    w(n+m+cp+1)= w(n+cp+1)*sq
    DO k=1,n
      w(k)= w(k)-w(n+m+cp+1)*y(cp*n+k)
    END DO
  END DO
!
  DO i=1,n
    w(i)=diag(i)*w(i)
  END DO
  DO ii=1,bound
    yr= ddot(n,y(cp*n+1),w)
    beta= w(n+cp+1)*yr
    DO k=1,n
      w(k)= w(k)+s(cp*n+k)*(w(n+m+cp+1)-beta)
    END DO
    cp=cp+1
    IF (cp == m)cp=0
  END DO
!
!  ------------------------------------------------------------
!  STORE THE NEW DIRECTION IN S
!  ------------------------------------------------------------
!
  DO j=1,n
    s(point*n+j)= w(j)
  END DO
!
!  ------------------------------------------------------------
!  OBTAIN THE MINIMIZER OF THE FUNCTION ALONG THE
!  DIRECTION S BY USING THE LINE SEARCH ROUTINE OF VD05AD
!  ------------------------------------------------------------
  65  nfev=0
  stp=one
  IF (iter == 1) stp=stp1
  DO i=1,n
    w(i)=g(i)
  END DO
  72  CONTINUE
!
!  write(ount,'(4x,I2,a,I2,F12.5)') myproc,' before vd05ad',info,stp
  CALL vd05ad(n,x,f,g,s(point*n+1),stp,ftol,gtol,                       &
              xtol,stpmin,stpmax,maxfev,info,nfev,diag)
!  write(ount,'(4x,I2,a,I2,F12.5)') myproc,' after vd05ad',info,stp
  call mpbarrier
!
  IF (info == -1) THEN
    iflag=1
    RETURN
  END IF
  IF (info /= 1) GO TO 90
  nfun= nfun + nfev
!
!  ------------------------------------------------------------
!  COMPUTE THE NEW S AND Y
!  ------------------------------------------------------------
!
  npt=point*n
  DO i=1,n
    s(npt+i)= stp*s(npt+i)
    y(npt+i)= g(i)-w(i)
  END DO
  ys= ddot(n,y(npt+1),s(npt+1))
  yy= ddot(n,y(npt+1),y(npt+1))
  point=point+1
  IF (point == m)point=0
!
!  ------------------------------------------------------------
!  CONVERGENCE CHECK
!  ------------------------------------------------------------
!
  gnorm= ddot(n,g,g)
  gnorm=SQRT(gnorm)
  xnorm= ddot(n,x,x)
  xnorm=SQRT(xnorm)
! xnorm= MAX1(1.0,xnorm)
  xnorm= MAX(1.0,xnorm)

  IF (gnorm/xnorm <= eps) finish=.true.

  CALL va15bd(iprint,iter,nfun,n,m,x,f,g,stp,finish)

  IF (finish) THEN
    iflag=0
    RETURN
  END IF
  GO TO 8
!
!  ------------------------------------------------------------
!  END OF MAIN ITERATION LOOP. ERROR EXITS.
!  ------------------------------------------------------------
!
  90  IF(lp <= 0) RETURN
  IF (info == 0) THEN
    iflag= -1
    WRITE(lp,100)iflag
  ELSE IF (info == 2) THEN
    iflag= -2
    WRITE(lp,105)iflag
  ELSE IF (info == 3) THEN
    iflag= -3
    WRITE(lp,110)iflag
  ELSE IF (info == 4) THEN
    iflag= -4
    WRITE(lp,115)iflag
  ELSE IF (info == 5) THEN
    iflag= -5
    WRITE(lp,120)iflag
  ELSE IF (info == 6) THEN
    iflag= -6
    WRITE(lp,125)iflag
  END IF
  RETURN
!
  95  iflag= -7
  IF(lp > 0) WRITE(lp,135)iflag,i
  RETURN
  96  iflag= -8
  IF(lp > 0) WRITE(lp,140)iflag
!
!  ------------------------------------------------------------
!  FORMATS
!  ------------------------------------------------------------
!
  100  FORMAT(/' IFLAG= ',i2,/' IMPROPER INPUT PARAMETERS DURING',      &
              ' THE LINE SEARCH.')
  105  FORMAT(/' IFLAG= ',i2,/' RELATIVE WIDTH OF THE INTERVAL OF',     &
              ' UNCERTAINTY IN THE LINE SEARCH'                         &
              /'IS OF THE ORDER OF machine roundoff.')
  110  FORMAT(/' IFLAG= ',i2,/' NUMBER OF CALLS TO FUNCTION IN THE',    &
              ' LINE SEARCH HAS REACHED 20.')
  115  FORMAT(/' IFLAG= ',i2,/' THE STEP IN THE LINE SEARCH IS',        &
              ' TOO SMALL.')
  120  FORMAT(/' IFLAG= ',i2,/' THE STEP IN THE LINE SEARCH IS',        &
              ' TOO LARGE.')
  125  FORMAT(/' IFLAG= ',i2,/' ROUNDING ERRORS PREVENT FURTHER',       &
              ' PROGRESS IN THE LINE SEARCH.')
  135  FORMAT(/' IFLAG= ',i2,/' THE',i5,'-TH DIAGONAL ELEMENT OF THE',  &
              ' INVERSE HESSIAN APPROXIMATION IS NOT POSITIVE')
  140  FORMAT(/' IFLAG= ',i2,/' IMPROPER INPUT PARAMETERS (N OR M',     &
              ' ARE NOT POSITIVE)')
  145  FORMAT(/'  GTOL IS LESS THAN OR EQUAL TO 1.D-04',                &
              / 'IT HAS BEEN RESET TO 1.D-02')
  RETURN
END SUBROUTINE va15ad

!
!#######################################################################
!
SUBROUTINE va15bd(iprint,iter,nfun,n,m,x,f,g,stp,finish)

  IMPLICIT NONE
!
!  ---------------------------------------------------------------------
!  THIS ROUTINE PRINTS MONITORING INFORMATION. THE FREQUENCY AND AMOUNT
!  OF OUTPUT ARE SPECIFIED AS FOLLOWS:
!
!  IPRINT(1) < 0 : NO OUTPUT IS GENERATED
!  IPRINT(1) = 0 : OUTPUT ONLY AT FIRST AND LAST ITERATION
!  IPRINT(1) 0 : OUTPUT EVERY IPRINT(1) ITERATION
!  IPRINT(2) = 0 : ITERATION COUNT, FUNCTION VALUE, NORM OF THE GRADIENT
!                  ,NUMBER OF FUNCTION CALLS AND STEP LENGTH
!  IPRINT(2) = 1 : + VECTOR OF VARIABLES AND GRADIENT VECTOR AT THE
!                    INITIAL POINT
!  IPRINT(2) = 2 : + VECTOR OF VARIABLES
!  IPRINT(2) = 3 : + GRADIENT VECTOR
!  ---------------------------------------------------------------------
!
  INTEGER, INTENT(IN) :: iprint(2),iter,nfun, n, m
  REAL,    INTENT(IN) :: f,stp
  REAL,    INTENT(IN) :: x(n),g(n)
  LOGICAL, INTENT(OUT):: finish

  REAL    :: gnorm,factor,gtol
  INTEGER :: prob,mp, lp
  COMMON /set/ factor,prob
  COMMON /va15dd/mp,lp, gtol

  INTEGER :: i

  REAL    :: ddot
!
  IF (iprint(1) < 0) RETURN
  gnorm= ddot(n,g,g)
  gnorm= SQRT(gnorm)
  IF (iter == 0)THEN
    WRITE(mp,10)
    WRITE(mp,20) prob,n,m
    WRITE(mp,30) f,gnorm
    IF (iprint(2) >= 1)THEN
      WRITE(mp,40)
      WRITE(mp,50) (x(i),i=1,n)
      WRITE(mp,60)
      WRITE(mp,50) (g(i),i=1,n)
    END IF
    WRITE(mp,10)
    WRITE(mp,70)
  ELSE
    IF ((iprint(1) == 0).AND.(iter /= 1.AND..NOT.finish)) RETURN
    IF (iprint(1) /= 0) THEN
      IF(MOD(iter-1,iprint(1)) == 0 .OR. finish) THEN
        WRITE(mp,80)iter,nfun,f,gnorm,stp
      ELSE
        RETURN
      END IF
    ELSE
      WRITE(mp,80)iter,nfun,f,gnorm,stp
    END IF
    IF (iprint(2) == 2.OR.iprint(2) == 3)THEN
      IF (finish)THEN
        WRITE(mp,90)
      ELSE
        WRITE(mp,40)
      END IF
      WRITE(mp,50)(x(i),i=1,n)
      IF (iprint(2) == 3)THEN
        WRITE(mp,60)
        WRITE(mp,50)(g(i),i=1,n)
      END IF
    END IF
    IF (finish) WRITE(mp,100)
  END IF
!
  10   FORMAT('*************************************************')
  20   FORMAT(' PROB=',i3,'   N=',i9,'   NUMBER OF CORRECTIONS=',i2)
  30   FORMAT(' F= ',1PD10.3,'   GNORM= ',1PD10.3)
  40   FORMAT(' VECTOR X= ')
  50   FORMAT(6(2X,1PD10.3))
  60   FORMAT(' GRADIENT VECTOR G= ')
  70   FORMAT(/'   I   NFN',4X,'FUNC',8X,'GNORM',7X,'STEPLENGTH'/)
  80   FORMAT(2(i4,1X),3X,3(1PD10.3,2X))
  90   FORMAT(' FINAL POINT X= ')
  100  FORMAT(/' THE MINIMIZATION TERMINATED WITHOUT DETECTING ERRORS.', &
              /' IFLAG = 0')
!
  RETURN
END SUBROUTINE va15bd
!
!   ----------------------------------------------------------
!   DATA BLOCK
!   ----------------------------------------------------------
!
  BLOCK DATA va15cd
  COMMON /va15dd/mp,lp, gtol
  INTEGER :: lp
  REAL :: gtol
  DATA mp,lp,gtol/6,6,9.0E-01/
  END

!
!#######################################################################
!
SUBROUTINE vd05ad(n,x,f,g,s,stp,ftol,gtol,xtol,                         &
                  stpmin,stpmax,maxfev,info,nfev,wa)

  IMPLICIT NONE

!  **********
!
!  SUBROUTINE VD05AD
!
!  THE PURPOSE OF VD05AD IS TO FIND A STEP WHICH SATISFIES
!  A SUFFICIENT DECREASE CONDITION AND A CURVATURE CONDITION.
!  THE USER MUST PROVIDE A SUBROUTINE WHICH CALCULATES THE
!  FUNCTION AND THE GRADIENT.
!
!  AT EACH STAGE THE SUBROUTINE UPDATES AN INTERVAL OF
!  UNCERTAINTY WITH ENDPOINTS STX AND STY. THE INTERVAL OF
!  UNCERTAINTY IS INITIALLY CHOSEN SO THAT IT CONTAINS A
!  MINIMIZER OF THE MODIFIED FUNCTION
!
!       F(X+STP*S) - F(X) - FTOL*STP*(GRADF(X)'S).
!
!  IF A STEP IS OBTAINED FOR WHICH THE MODIFIED FUNCTION
!  HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE DERIVATIVE,
!  THEN THE INTERVAL OF UNCERTAINTY IS CHOSEN SO THAT IT
!  CONTAINS A MINIMIZER OF F(X+STP*S).
!
!  THE ALGORITHM IS DESIGNED TO FIND A STEP WHICH SATISFIES
!  THE SUFFICIENT DECREASE CONDITION
!
!        F(X+STP*S) .LE. F(X) + FTOL*STP*(GRADF(X)'S),
!
!  AND THE CURVATURE CONDITION
!
!        ABS(GRADF(X+STP*S)'S)) .LE. GTOL*ABS(GRADF(X)'S).
!
!  IF FTOL IS LESS THAN GTOL AND IF, FOR EXAMPLE, THE FUNCTION
!  IS BOUNDED BELOW, THEN THERE IS ALWAYS A STEP WHICH SATISFIES
!  BOTH CONDITIONS. IF NO STEP CAN BE FOUND WHICH SATISFIES BOTH
!  CONDITIONS, THEN THE ALGORITHM USUALLY STOPS WHEN ROUNDING
!  ERRORS PREVENT FURTHER PROGRESS. IN THIS CASE STP ONLY
!  SATISFIES THE SUFFICIENT DECREASE CONDITION.
!
!  THE SUBROUTINE STATEMENT IS
!
!     SUBROUTINE VD05AD(N,X,F,G,S,STP,FTOL,GTOL,XTOL,
!                       STPMIN,STPMAX,MAXFEV,INFO,NFEV,WA)
!  WHERE
!
!    N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!      OF VARIABLES.
!
!    X IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE
!      BASE POINT FOR THE LINE SEARCH. ON OUTPUT IT CONTAINS
!      X + STP*S.
!
!    F IS A VARIABLE. ON INPUT IT MUST CONTAIN THE VALUE OF F
!      AT X. ON OUTPUT IT CONTAINS THE VALUE OF F AT X + STP*S.
!
!    G IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE
!      GRADIENT OF F AT X. ON OUTPUT IT CONTAINS THE GRADIENT
!      OF F AT X + STP*S.
!
!    S IS AN INPUT ARRAY OF LENGTH N WHICH SPECIFIES THE
!      SEARCH DIRECTION.
!
!    STP IS A NONNEGATIVE VARIABLE. ON INPUT STP CONTAINS AN
!      INITIAL ESTIMATE OF A SATISFACTORY STEP. ON OUTPUT
!      STP CONTAINS THE FINAL ESTIMATE.
!
!    FTOL AND GTOL ARE NONNEGATIVE INPUT VARIABLES. TERMINATION
!      OCCURS WHEN THE SUFFICIENT DECREASE CONDITION AND THE
!      DIRECTIONAL DERIVATIVE CONDITION ARE SATISFIED.
!
!    XTOL IS A NONNEGATIVE INPUT VARIABLE. TERMINATION OCCURS
!      WHEN THE RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY
!      IS AT MOST XTOL.
!
!    STPMIN AND STPMAX ARE NONNEGATIVE INPUT VARIABLES WHICH
!      SPECIFY LOWER AND UPPER BOUNDS FOR THE STEP.
!
!    MAXFEV IS A POSITIVE INTEGER INPUT VARIABLE. TERMINATION
!      OCCURS WHEN THE NUMBER OF CALLS TO FCN IS AT LEAST
!      MAXFEV BY THE END OF AN ITERATION.
!
!    INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS:
!
!      INFO = 0  IMPROPER INPUT PARAMETERS.
!
!      INFO =-1  A RETURN IS MADE TO COMPUTE THE FUNCTION AND GRADIENT.
!
!      INFO = 1  THE SUFFICIENT DECREASE CONDITION AND THE
!                DIRECTIONAL DERIVATIVE CONDITION HOLD.
!
!      INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY
!                IS AT MOST XTOL.
!
!      INFO = 3  NUMBER OF CALLS TO FCN HAS REACHED MAXFEV.
!
!      INFO = 4  THE STEP IS AT THE LOWER BOUND STPMIN.
!
!      INFO = 5  THE STEP IS AT THE UPPER BOUND STPMAX.
!
!      INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS.
!                THERE MAY NOT BE A STEP WHICH SATISFIES THE
!                SUFFICIENT DECREASE AND CURVATURE CONDITIONS.
!                TOLERANCES MAY BE TOO SMALL.
!
!    NFEV IS AN INTEGER OUTPUT VARIABLE SET TO THE NUMBER OF
!      CALLS TO FCN.
!
!    WA IS A WORK ARRAY OF LENGTH N.
!
!  SUBPROGRAMS CALLED
!
!    HARWELL-SUPPLIED...VD05BD
!
!    FORTRAN-SUPPLIED...ABS,MAX,MIN
!
!  ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
!  JORGE J. MORE', DAVID J. THUENTE
!
!  **********

  INTEGER, INTENT(IN)    :: n,maxfev
  REAL,    INTENT(IN)    :: ftol,gtol,xtol,stpmin,stpmax
  REAL,    INTENT(IN)    :: s(n)

  REAL,    INTENT(INOUT) :: x(n),g(n),wa(n)
  INTEGER, INTENT(INOUT) :: info
  REAL,    INTENT(INOUT) :: f, stp

  INTEGER, INTENT(OUT)   :: nfev

  SAVE

  INTEGER :: infoc,j
  LOGICAL :: brackt,stage1
  REAL    :: dg,dgm,dginit,dgtest,dgx,dgxm,dgy,dgym,                    &
             finit,ftest1,fm,fx,fxm,fy,fym,p5,p66,stx,sty,              &
             stmin,stmax,width,width1,xtrapf,zero
  DATA p5,p66,xtrapf,zero /0.5E0,0.66E0,4.0E0,0.0E0/

!
! external function
!
  REAL :: ddot

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF(info == -1) GO TO 45
  infoc = 1
!
!  CHECK THE INPUT PARAMETERS FOR ERRORS.
!
  IF (n <= 0 .OR. stp <= zero .OR. ftol < zero .OR.                     &
      gtol < zero .OR. xtol < zero .OR. stpmin < zero                   &
      .OR. stpmax < stpmin .OR. maxfev <= 0) RETURN
!
!  COMPUTE THE INITIAL GRADIENT IN THE SEARCH DIRECTION
!  AND CHECK THAT S IS A DESCENT DIRECTION.
!
!  dginit = zero
!  DO j = 1, n
!    dginit = dginit + g(j)*s(j)
!  END DO
  dginit = ddot(n,g,s)
  IF (dginit >= zero) RETURN
!
!  INITIALIZE LOCAL VARIABLES.
!
  brackt = .false.
  stage1 = .true.
  nfev = 0
  finit = f
  dgtest = ftol*dginit
  width = stpmax - stpmin
  width1 = width/p5
  DO j = 1, n
    wa(j) = x(j)
  END DO
!
!  THE VARIABLES STX, FX, DGX CONTAIN THE VALUES OF THE STEP,
!  FUNCTION, AND DIRECTIONAL DERIVATIVE AT THE BEST STEP.
!  THE VARIABLES STY, FY, DGY CONTAIN THE VALUE OF THE STEP,
!  FUNCTION, AND DERIVATIVE AT THE OTHER ENDPOINT OF
!  THE INTERVAL OF UNCERTAINTY.
!  THE VARIABLES STP, F, DG CONTAIN THE VALUES OF THE STEP,
!  FUNCTION, AND DERIVATIVE AT THE CURRENT STEP.
!
  stx = zero
  fx  = finit
  dgx = dginit
  sty = zero
  fy  = finit
  dgy = dginit
!
!  START OF ITERATION.
!
  30 CONTINUE
!
!     SET THE MINIMUM AND MAXIMUM STEPS TO CORRESPOND
!     TO THE PRESENT INTERVAL OF UNCERTAINTY.
!
  IF (brackt) THEN
    stmin = MIN(stx,sty)
    stmax = MAX(stx,sty)
  ELSE
    stmin = stx
    stmax = stp + xtrapf*(stp - stx)
  END IF
!
!     FORCE THE STEP TO BE WITHIN THE BOUNDS STPMAX AND STPMIN.
!
  stp = MAX(stp,stpmin)
  stp = MIN(stp,stpmax)
!
!  WRITE(0,*)'stp =================',stp
!
!     IF AN UNUSUAL TERMINATION IS TO OCCUR THEN LET
!     STP BE THE LOWEST POINT OBTAINED SO FAR.
!
  IF ((brackt .AND. (stp <= stmin .OR. stp >= stmax))                   &
      .OR. nfev >= maxfev-1 .OR. infoc == 0                             &
      .OR. (brackt .AND. stmax-stmin <= xtol*stmax)) stp = stx
!
!     EVALUATE THE FUNCTION AND GRADIENT AT STP
!     AND COMPUTE THE DIRECTIONAL DERIVATIVE.
!
  DO j = 1, n
    x(j) = wa(j) + stp*s(j)
  END DO
  info=-1
  RETURN
!
  45    info=0
  nfev = nfev + 1
  dg = ddot(n,g,s)
  
  ftest1 = finit + stp*dgtest
!
!     TEST FOR CONVERGENCE.
!
  IF ((brackt .AND. (stp <= stmin .OR. stp >= stmax)) .OR. infoc == 0) info = 6
  IF (stp == stpmax .AND. f <= ftest1 .AND. dg <= dgtest) info = 5
  IF (stp == stpmin .AND. (f > ftest1 .OR. dg >= dgtest)) info = 4
  IF (nfev >= maxfev) info = 3
  IF (brackt .AND. stmax-stmin <= xtol*stmax) info = 2
  IF (f <= ftest1 .AND. ABS(dg) <= gtol*(-dginit)) info = 1
!
!     CHECK FOR TERMINATION.
!
  IF (info /= 0) RETURN
!
!     IN THE FIRST STAGE WE SEEK A STEP FOR WHICH THE MODIFIED
!     FUNCTION HAS A NONPOSITIVE VALUE AND NONNEGATIVE DERIVATIVE.
!
  IF (stage1 .AND. f <= ftest1 .AND.                                    &
      dg >= MIN(ftol,gtol)*dginit) stage1 = .false.
!
!     A MODIFIED FUNCTION IS USED TO PREDICT THE STEP ONLY IF
!     WE HAVE NOT OBTAINED A STEP FOR WHICH THE MODIFIED
!     FUNCTION HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE
!     DERIVATIVE, AND IF A LOWER FUNCTION VALUE HAS BEEN
!     OBTAINED BUT THE DECREASE IS NOT SUFFICIENT.
!
  IF (stage1 .AND. f <= fx .AND. f > ftest1) THEN
!
!        DEFINE THE MODIFIED FUNCTION AND DERIVATIVE VALUES.
!
    fm = f - stp*dgtest
    fxm = fx - stx*dgtest
    fym = fy - sty*dgtest
    dgm = dg - dgtest
    dgxm = dgx - dgtest
    dgym = dgy - dgtest
!
!        CALL CSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
!        AND TO COMPUTE THE NEW STEP.
!
    CALL vd05bd(stx,fxm,dgxm,sty,fym,dgym,stp,fm,dgm,                   &
               brackt,stmin,stmax,infoc)
!
!        RESET THE FUNCTION AND GRADIENT VALUES FOR F.
!
    fx = fxm + stx*dgtest
    fy = fym + sty*dgtest
    dgx = dgxm + dgtest
    dgy = dgym + dgtest
  ELSE
!
!        CALL VD05BD TO UPDATE THE INTERVAL OF UNCERTAINTY
!        AND TO COMPUTE THE NEW STEP.
!
    CALL vd05bd(stx,fx,dgx,sty,fy,dgy,stp,f,dg,                         &
                brackt,stmin,stmax,infoc)
  END IF
!
!     FORCE A SUFFICIENT DECREASE IN THE SIZE OF THE
!     INTERVAL OF UNCERTAINTY.
!
  IF (brackt) THEN
    IF (ABS(sty-stx) >= p66*width1) stp = stx + p5*(sty - stx)
    width1 = width
    width = ABS(sty-stx)
  END IF
!
!     END OF ITERATION.
!
  GO TO 30
END SUBROUTINE vd05ad
!
!#######################################################################
!
SUBROUTINE vd05bd(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,                 &
                  stpmin,stpmax,info)
  IMPLICIT NONE
!  **********
!
!  SUBROUTINE VD05BD
!
!  THE PURPOSE OF VD05BD IS TO COMPUTE A SAFEGUARDED STEP FOR
!  A LINESEARCH AND TO UPDATE AN INTERVAL OF UNCERTAINTY FOR
!  A MINIMIZER OF THE FUNCTION.
!
!  THE PARAMETER STX CONTAINS THE STEP WITH THE LEAST FUNCTION
!  VALUE. THE PARAMETER STP CONTAINS THE CURRENT STEP. IT IS
!  ASSUMED THAT THE DERIVATIVE AT STX IS NEGATIVE IN THE
!  DIRECTION OF THE STEP. IF BRACKT IS SET TRUE THEN A
!  MINIMIZER HAS BEEN BRACKETED IN AN INTERVAL OF UNCERTAINTY
!  WITH ENDPOINTS STX AND STY.
!
!  THE SUBROUTINE STATEMENT IS
!
!    SUBROUTINE VD05BD(STX,FX,DX,STY,FY,DY,STP,FP,DP,BRACKT,
!                     STPMIN,STPMAX,INFO)
!
!  WHERE
!
!    STX, FX, AND DX ARE VARIABLES WHICH SPECIFY THE STEP,
!      THE FUNCTION, AND THE DERIVATIVE AT THE BEST STEP OBTAINED
!      SO FAR. THE DERIVATIVE MUST BE NEGATIVE IN THE DIRECTION
!      OF THE STEP, THAT IS, DX AND STP-STX MUST HAVE OPPOSITE
!      SIGNS. ON OUTPUT THESE PARAMETERS ARE UPDATED APPROPRIATELY.
!
!    STY, FY, AND DY ARE VARIABLES WHICH SPECIFY THE STEP,
!      THE FUNCTION, AND THE DERIVATIVE AT THE OTHER ENDPOINT OF
!      THE INTERVAL OF UNCERTAINTY. ON OUTPUT THESE PARAMETERS ARE
!      UPDATED APPROPRIATELY.
!
!    STP, FP, AND DP ARE VARIABLES WHICH SPECIFY THE STEP,
!      THE FUNCTION, AND THE DERIVATIVE AT THE CURRENT STEP.
!      IF BRACKT IS SET TRUE THEN ON INPUT STP MUST BE
!      BETWEEN STX AND STY. ON OUTPUT STP IS SET TO THE NEW STEP.
!
!    BRACKT IS A LOGICAL VARIABLE WHICH SPECIFIES IF A MINIMIZER
!      HAS BEEN BRACKETED. IF THE MINIMIZER HAS NOT BEEN BRACKETED
!      THEN ON INPUT BRACKT MUST BE SET FALSE. IF THE MINIMIZER
!      IS BRACKETED THEN ON OUTPUT BRACKT IS SET TRUE.
!
!    STPMIN AND STPMAX ARE INPUT VARIABLES WHICH SPECIFY LOWER
!      AND UPPER BOUNDS FOR THE STEP.
!
!    INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS:
!      IF INFO = 1,2,3,4,5, THEN THE STEP HAS BEEN COMPUTED
!      ACCORDING TO ONE OF THE FIVE CASES BELOW. OTHERWISE
!      INFO = 0, AND THIS INDICATES IMPROPER INPUT PARAMETERS.
!
!  SUBPROGRAMS CALLED
!
!    FORTRAN-SUPPLIED ... ABS,MAX,MIN,SQRT
!
!  ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
!  JORGE J. MORE', DAVID J. THUENTE
!
!  **********
  REAL,    INTENT(INOUT) :: stx,fx,dx,sty,fy,dy,stp,fp,dp
  REAL,    INTENT(IN)    :: stpmin,stpmax
  LOGICAL, INTENT(INOUT) :: brackt
  INTEGER, INTENT(OUT)   :: info

  LOGICAL :: bound
  REAL    :: gamma,p,q,r,s,sgnd,stpc,stpf,stpq,theta

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  info = 0
!
!  CHECK THE INPUT PARAMETERS FOR ERRORS.
!
  IF ((brackt .AND. (stp <= MIN(stx,sty) .OR. stp >= MAX(stx,sty))) .OR. &
      dx*(stp-stx) >= 0.0 .OR. stpmax < stpmin) RETURN
!
!  DETERMINE IF THE DERIVATIVES HAVE OPPOSITE SIGN.
!
  sgnd = dp*(dx/ABS(dx))
!
!  FIRST CASE. A HIGHER FUNCTION VALUE.
!  THE MINIMUM IS BRACKETED. IF THE CUBIC STEP IS CLOSER
!  TO STX THAN THE QUADRATIC STEP, THE CUBIC STEP IS TAKEN,
!  ELSE THE AVERAGE OF THE CUBIC AND QUADRATIC STEPS IS TAKEN.
!
  IF (fp > fx) THEN
    info = 1
    bound = .true.
    theta = 3*(fx - fp)/(stp - stx) + dx + dp
    s = MAX(ABS(theta),ABS(dx),ABS(dp))
    gamma = s*SQRT((theta/s)**2 - (dx/s)*(dp/s))
    IF (stp < stx) gamma = -gamma
    p = (gamma - dx) + theta
    q = ((gamma - dx) + gamma) + dp
    r = p/q
    stpc = stx + r*(stp - stx)
    stpq = stx + ((dx/((fx-fp)/(stp-stx)+dx))/2)*(stp - stx)
    IF (ABS(stpc-stx) < ABS(stpq-stx)) THEN
      stpf = stpc
    ELSE
      stpf = stpc + (stpq - stpc)/2
    END IF
    brackt = .true.
!
!  SECOND CASE. A LOWER FUNCTION VALUE AND DERIVATIVES OF
!  OPPOSITE SIGN. THE MINIMUM IS BRACKETED. IF THE CUBIC
!  STEP IS CLOSER TO STX THAN THE QUADRATIC (SECANT) STEP,
!  THE CUBIC STEP IS TAKEN, ELSE THE QUADRATIC STEP IS TAKEN.
!
  ELSE IF (sgnd < 0.0) THEN
    info = 2
    bound = .false.
    theta = 3*(fx - fp)/(stp - stx) + dx + dp
    s = MAX(ABS(theta),ABS(dx),ABS(dp))
    gamma = s*SQRT((theta/s)**2 - (dx/s)*(dp/s))
    IF (stp > stx) gamma = -gamma
    p = (gamma - dp) + theta
    q = ((gamma - dp) + gamma) + dx
    r = p/q
    stpc = stp + r*(stx - stp)
    stpq = stp + (dp/(dp-dx))*(stx - stp)
    IF (ABS(stpc-stp) > ABS(stpq-stp)) THEN
      stpf = stpc
    ELSE
      stpf = stpq
    END IF
    brackt = .true.
!
!  THIRD CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
!  SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DECREASES.
!  THE CUBIC STEP IS ONLY USED IF THE CUBIC TENDS TO INFINITY
!  IN THE DIRECTION OF THE STEP OR IF THE MINIMUM OF THE CUBIC
!  IS BEYOND STP. OTHERWISE THE CUBIC STEP IS DEFINED TO BE
!  EITHER STPMIN OR STPMAX. THE QUADRATIC (SECANT) STEP IS ALSO
!  COMPUTED AND IF THE MINIMUM IS BRACKETED THEN THE THE STEP
!  CLOSEST TO STX IS TAKEN, ELSE THE STEP FARTHEST AWAY IS TAKEN.
!
  ELSE IF (ABS(dp) < ABS(dx)) THEN
    info = 3
    bound = .true.
    theta = 3*(fx - fp)/(stp - stx) + dx + dp
    s = MAX(ABS(theta),ABS(dx),ABS(dp))
!
!     THE CASE GAMMA = 0 ONLY ARISES IF THE CUBIC DOES NOT TEND
!     TO INFINITY IN THE DIRECTION OF THE STEP.
!
    gamma = s*SQRT(MAX(0.0,(theta/s)**2 - (dx/s)*(dp/s)))
    IF (stp > stx) gamma = -gamma
    p = (gamma - dp) + theta
    q = (gamma + (dx - dp)) + gamma
    r = p/q
    IF (r < 0.0 .AND. gamma /= 0.0) THEN
      stpc = stp + r*(stx - stp)
    ELSE IF (stp > stx) THEN
      stpc = stpmax
    ELSE
      stpc = stpmin
    END IF
    stpq = stp + (dp/(dp-dx))*(stx - stp)
    IF (brackt) THEN
      IF (ABS(stp-stpc) < ABS(stp-stpq)) THEN
        stpf = stpc
      ELSE
        stpf = stpq
      END IF
    ELSE
      IF (ABS(stp-stpc) > ABS(stp-stpq)) THEN
        stpf = stpc
      ELSE
        stpf = stpq
      END IF
    END IF
!
!  FOURTH CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
!  SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DOES
!  NOT DECREASE. IF THE MINIMUM IS NOT BRACKETED, THE STEP
!  IS EITHER STPMIN OR STPMAX, ELSE THE CUBIC STEP IS TAKEN.
!
  ELSE
    info = 4
    bound = .false.
    IF (brackt) THEN
      theta = 3*(fp - fy)/(sty - stp) + dy + dp
      s = MAX(ABS(theta),ABS(dy),ABS(dp))
      gamma = s*SQRT((theta/s)**2 - (dy/s)*(dp/s))
      IF (stp > sty) gamma = -gamma
      p = (gamma - dp) + theta
      q = ((gamma - dp) + gamma) + dy
      r = p/q
      stpc = stp + r*(sty - stp)
      stpf = stpc
    ELSE IF (stp > stx) THEN
      stpf = stpmax
    ELSE
      stpf = stpmin
    END IF
  END IF
!
!  UPDATE THE INTERVAL OF UNCERTAINTY. THIS UPDATE DOES NOT
!  DEPEND ON THE NEW STEP OR THE CASE ANALYSIS ABOVE.
!
  IF (fp > fx) THEN
    sty = stp
    fy = fp
    dy = dp
  ELSE
    IF (sgnd < 0.0) THEN
      sty = stx
      fy = fx
      dy = dx
    END IF
    stx = stp
    fx = fp
    dx = dp
  END IF
!
!  COMPUTE THE NEW STEP AND SAFEGUARD IT.
!
  stpf = MIN(stpmax,stpf)
  stpf = MAX(stpmin,stpf)
  stp = stpf
  IF (brackt .AND. bound) THEN
    IF (sty > stx) THEN
      stp = MIN(stx+0.66*(sty-stx),stp)
    ELSE
      stp = MAX(stx+0.66*(sty-stx),stp)
    END IF
  END IF

  RETURN
END SUBROUTINE vd05bd

!
!#######################################################################
!
REAL FUNCTION ddot(n,x,y)

  IMPLICIT NONE
!
!   -------------------------------------------------------
!   THIS FUNCTION COMPUTES THE INNER PRODUCT OF TWO VECTORS
!   -------------------------------------------------------
!
  INTEGER, INTENT(IN) :: n
  REAL,    INTENT(IN) :: x(n),y(n)

!  REAL    :: prod
  DOUBLE PRECISION :: prod
  INTEGER          :: m,i,mp1

  ddot = 0.0
  IF (n < 0) RETURN

  prod = 0.0
  m = MOD(n,5)
  IF ( m /= 0 ) THEN
    DO i=1,m
      prod= prod + x(i)*y(i)
    END DO
  END IF
  IF ( n >= 5 ) THEN
    mp1 = m + 1
    DO i = mp1,n,5
      prod = prod + x(i    )*y(i    ) + x(i + 1)*y(i + 1) +             &
                    x(i + 2)*y(i + 2) + x(i + 3)*y(i + 3) +             &
                    x(i + 4)*y(i + 4)
    END DO
  END IF

  CALL mpsumdp(prod,1)

  ddot = prod

! CALL mptotal(ddot)

  RETURN
END FUNCTION ddot
