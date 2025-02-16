  FUNCTION pimach(dum)
!***BEGIN PROLOGUE  PIMACH
!
!  This subprogram supplies the value of the constant PI correct to
!  machine precision where
!
!  PI=3.1415926535897932384626433832795028841971693993751058209749446
!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  PIMACH
!
!***FIRST EXECUTABLE STATEMENT  PIMACH
  pimach = 3.14159265358979
  RETURN
  END FUNCTION pimach


SUBROUTINE vcost(m,n,x,xt,mdimx,wsave)

!***BEGIN PROLOGUE  VCOST
!***DATE WRITTEN   860701   (YYMMDD)
!***REVISION DATE  900509   (YYMMDD)
!***CATEGORY NO.  J1A3
!***KEYWORDS  FAST FOURIER TRANSFORM, COSINE TRANSFORM, MULTIPLE
!          SEQUENCES
!***AUTHOR  BOISVERT, R. F. (NIST)
!***PURPOSE  Cosine transform of one or more real, even sequences.
!***DESCRIPTION
!
!  Subroutine VCOST computes the discrete Fourier cosine transform
!  of M even sequences X(J,I), J=1,...,M.  The transform is defined
!  below at output parameter X.
!
!  The array WSAVE which is used by subroutine VCOST must be
!  initialized by calling subroutine VCOSTI(N,WSAVE).
!
!  Input Parameters
!
!  M       the number of sequences to be transformed.
!
!  N       the length of the sequence to be transformed.  N must be
!       greater than 1.  The method is most efficient when N-1 is
!       is a product of small primes.
!
!  X       an array of size at least X(MDIMX,N) which contains the
!       the sequences to be transformed.  The sequences are stored
!       in the ROWS of X.  Thus, the Jth sequence is stored in
!       X(J,I), I=1,..,N.
!
!  XT      a work array of size at least XT(MDIMX,N-1).
!
!  MDIMX   the first dimension of the array X exactly as it appears in
!       the calling program.
!
!  WSAVE   a work array which must be dimensioned at least 3*N+15
!       in the program that calls VCOST.  The WSAVE array must be
!       initialized by calling subroutine VCOSTI(N,WSAVE), and a
!       different WSAVE array must be used for each different
!       value of N.  This initialization does not have to be
!       repeated so long as N remains unchanged.  Thus subsequent
!       transforms can be obtained faster than the first.
!
!  Output Parameters
!
!  X       For I=1,...,N and J=1,...,M
!
!          X(J,I) = ( X(J,1)+(-1)**(I-1)*X(J,N)
!
!            + the sum from K=2 to K=N-1
!
!              2*X(J,K)*COS((K-1)*(I-1)*PI/(N-1)) )/SQRT(2*(N-1))
!
!  WSAVE   contains initialization calculations which must not be
!       destroyed between calls of VCOST.
!
!  -----------------------------------------------------------------
!
!  NOTE  -  A call of VCOST followed immediately by another call
!        of VCOST will return the original sequences X.  Thus,
!        VCOST is the correctly normalized inverse of itself.
!
!  -----------------------------------------------------------------
!
!  VCOST is a straightforward extension of the subprogram COST to
!  handle M simultaneous sequences.  The scaling of the sequences
!  computed by VCOST is different than that of COST.  COST was
!  originally developed by P. N. Swarztrauber of NCAR.
!
!***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
!            Computations, (G. Rodrigue, ed.), Academic Press, 1982,
!            pp. 51-83.
!***ROUTINES CALLED  VRFFTF
!***END PROLOGUE  VCOST
  DIMENSION       x(mdimx,*), xt(mdimx,*), wsave(*)
!***FIRST EXECUTABLE STATEMENT  VCOST
  IF (m <= 0)  GO TO 900
  IF (n <= 1)  GO TO 900
  IF (n > 3)  GO TO 400
  IF (n == 3)  GO TO 300
!
!  CASE  N = 2
!
  scale = SQRT(0.50E0)
  DO j=1,m
    x1h = scale*(x(j,1)+x(j,2))
    x(j,2) = scale*(x(j,1)-x(j,2))
    x(j,1) = x1h
  END DO
  GO TO 900
!
!  CASE  N = 3
!
  300 CONTINUE
  scale = 0.50E0
  DO j=1,m
    x1p3 = x(j,1)+x(j,3)
    tx2 = x(j,2)+x(j,2)
    x(j,2) = scale*(x(j,1)-x(j,3))
    x(j,1) = scale*(x1p3+tx2)
    x(j,3) = scale*(x1p3-tx2)
  END DO
  GO TO 900
!
!  CASE  N .GT. 3
!
!  ... PREPROCESSING
!
  400 CONTINUE
  nm1 = n-1
  np1 = n+1
  ns2 = n/2
  DO j=1,m
    xt(j,1) = x(j,1)-x(j,n)
    x(j,1) = x(j,1)+x(j,n)
  END DO
  DO k=2,ns2
    kc = np1-k
    DO j=1,m
      t1 = x(j,k)+x(j,kc)
      t2 = x(j,k)-x(j,kc)
      xt(j,1) = xt(j,1)+wsave(kc)*t2
      t2 = wsave(k)*t2
      x(j,k) = t1-t2
      x(j,kc) = t1+t2
    END DO
  END DO
  modn = MOD(n,2)
  IF (modn /= 0) THEN
    DO j=1,m
      x(j,ns2+1) = x(j,ns2+1)+x(j,ns2+1)
    END DO
  END IF
  DO j=1,m
    x(j,n) = xt(j,1)
  END DO
!
!  ... REAL PERIODIC TRANSFORM
!
  CALL vrfftf (m,nm1,x,xt,mdimx,wsave(np1))
!
!  ... POSTPROCESSING
!
  factor = 1.0/SQRT(REAL(nm1))
  DO j=1,m
    xt(j,1) = x(j,2)
    x(j,2) = factor*x(j,n)
  END DO
  DO i=4,n,2
    DO j=1,m
      xi = x(j,i)
      x(j,i) = x(j,i-2)-x(j,i-1)
      x(j,i-1) = xt(j,1)
      xt(j,1) = xi
    END DO
  END DO
  IF (modn /= 0) THEN
    DO j=1,m
      x(j,n) = xt(j,1)
    END DO
  END IF
!
!  ... NORMALIZATION
!
  scale = SQRT(0.5)
  DO i=1,n
    DO j=1,m
      x(j,i) = scale*x(j,i)
    END DO
  END DO
!
!  EXIT
!
  900 CONTINUE
  RETURN
END SUBROUTINE vcost

SUBROUTINE vcosti(n,wsave)
!***BEGIN PROLOGUE  VCOSTI
!***DATE WRITTEN   860701   (YYMMDD)
!***REVISION DATE  900509   (YYMMDD)
!***CATEGORY NO.  J1A3
!***KEYWORDS  FAST FOURIER TRANSFORM, COSINE TRANSFORM, MULTIPLE
!          SEQUENCES
!***AUTHOR  BOISVERT, R. F. (NIST)
!***PURPOSE  Initialize for VCOST.
!***DESCRIPTION
!
!  Subroutine VCOSTI initializes the array WSAVE which is used in
!  subroutine VCOST.  The prime factorization of N together with
!  a tabulation of the trigonometric functions are computed and
!  stored in WSAVE.
!
!  Input Parameter
!
!  N       the length of the sequence to be transformed.  The method
!       is most efficient when N-1 is a product of small primes.
!
!  Output Parameter
!
!  WSAVE   a work array which must be dimensioned at least 3*N+15.
!       Different WSAVE arrays are required for different values
!       of N.  The contents of WSAVE must not be changed between
!       calls of VCOST.
!
!  -----------------------------------------------------------------
!
!  VCOSTI is a straightforward extension of the subprogram COSTI to
!  handle M simultaneous sequences.  COSTI was originally developed
!  by P. N. Swarztrauber of NCAR.
!
!***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
!            Computations, (G. Rodrigue, ed.), Academic Press, 1982,
!            pp. 51-83.
!***ROUTINES CALLED  VRFFTI
!***END PROLOGUE  VCOSTI
  DIMENSION       wsave(*)
!***FIRST EXECUTABLE STATEMENT  VCOSTI
  pi = pimach(1.0)
  IF (n <= 3) RETURN
  nm1 = n-1
  np1 = n+1
  ns2 = n/2
  dt = pi/REAL(nm1)
  fk = 0.
  DO k=2,ns2
    fk = fk+1.
    wsave(k) = 2.*SIN(fk*dt)
  END DO
  fk = 0.
  DO k=2,ns2
    kc = np1-k
    fk = fk+1.
    wsave(kc) = 2.*COS(fk*dt)
  END DO
  CALL vrffti (nm1,wsave(n+1))
  RETURN
END SUBROUTINE vcosti

SUBROUTINE vradb2 (mp,ido,l1,cc,ch,mdimc,wa1)
!
!  VRFFTPK, VERSION 1, AUGUST 1985
!
  DIMENSION  cc(mdimc,ido,2,l1)    ,ch(mdimc,ido,l1,2),                 &
                  wa1(ido)
  DO k=1,l1
    DO m=1,mp
      ch(m,1,k,1) = cc(m,1,1,k)+cc(m,ido,2,k)
      ch(m,1,k,2) = cc(m,1,1,k)-cc(m,ido,2,k)
    END DO
  END DO
  IF (ido-2 < 0) THEN
    GO TO   107
  ELSE IF (ido-2 == 0) THEN
    GO TO   105
  END IF
  idp2 = ido+2
  DO k=1,l1
    DO i=3,ido,2
      ic = idp2-i
      DO m=1,mp
        ch(m,i-1,k,1) = cc(m,i-1,1,k)+cc(m,ic-1,2,k)
        ch(m,i,k,1) = cc(m,i,1,k)-cc(m,ic,2,k)
        ch(m,i-1,k,2) = wa1(i-2)*(cc(m,i-1,1,k)-cc(m,ic-1,2,k))         &
            -wa1(i-1)*(cc(m,i,1,k)+cc(m,ic,2,k))
        ch(m,i,k,2) = wa1(i-2)*(cc(m,i,1,k)+cc(m,ic,2,k))+wa1(i-1)      &
            *(cc(m,i-1,1,k)-cc(m,ic-1,2,k))
      END DO
    END DO
  END DO
  IF (MOD(ido,2) == 1) RETURN
  105 DO k=1,l1
    DO m=1,mp
      ch(m,ido,k,1) = cc(m,ido,1,k)+cc(m,ido,1,k)
      ch(m,ido,k,2) = -(cc(m,1,2,k)+cc(m,1,2,k))
    END DO
  END DO
  107 RETURN
END SUBROUTINE vradb2

SUBROUTINE vradb3 (mp,ido,l1,cc,ch,mdimc,wa1,wa2)
!
!  VRFFTPK, VERSION 1, AUGUST 1985
!
  DIMENSION  cc(mdimc,ido,3,l1)    ,ch(mdimc,ido,l1,3),                 &
                  wa1(ido)   ,wa2(ido)
  arg=2.*pimach(1.0)/3.
  taur=COS(arg)
  taui=SIN(arg)
  DO k=1,l1
    DO m=1,mp
      ch(m,1,k,1) = cc(m,1,1,k)+2.*cc(m,ido,2,k)
      ch(m,1,k,2) = cc(m,1,1,k)+(2.*taur)*cc(m,ido,2,k)                 &
          -(2.*taui)*cc(m,1,3,k)
      ch(m,1,k,3) = cc(m,1,1,k)+(2.*taur)*cc(m,ido,2,k)                 &
          +2.*taui*cc(m,1,3,k)
    END DO
  END DO
  IF (ido == 1) RETURN
  idp2 = ido+2
  DO k=1,l1
    DO i=3,ido,2
      ic = idp2-i
      DO m=1,mp
        ch(m,i-1,k,1) = cc(m,i-1,1,k)+(cc(m,i-1,3,k)+cc(m,ic-1,2,k))
        ch(m,i,k,1) = cc(m,i,1,k)+(cc(m,i,3,k)-cc(m,ic,2,k))
        ch(m,i-1,k,2) = wa1(i-2)*                                       &
            ((cc(m,i-1,1,k)+taur*(cc(m,i-1,3,k)+cc(m,ic-1,2,k)))-       &
            (taui*(cc(m,i,3,k)+cc(m,ic,2,k))))                          &
                     -wa1(i-1)*                                         &
            ((cc(m,i,1,k)+taur*(cc(m,i,3,k)-cc(m,ic,2,k)))+             &
            (taui*(cc(m,i-1,3,k)-cc(m,ic-1,2,k))))
        ch(m,i,k,2) = wa1(i-2)*                                         &
            ((cc(m,i,1,k)+taur*(cc(m,i,3,k)-cc(m,ic,2,k)))+             &
            (taui*(cc(m,i-1,3,k)-cc(m,ic-1,2,k))))                      &
                    +wa1(i-1)*                                          &
            ((cc(m,i-1,1,k)+taur*(cc(m,i-1,3,k)+cc(m,ic-1,2,k)))-       &
            (taui*(cc(m,i,3,k)+cc(m,ic,2,k))))
        ch(m,i-1,k,3) = wa2(i-2)*                                       &
            ((cc(m,i-1,1,k)+taur*(cc(m,i-1,3,k)+cc(m,ic-1,2,k)))+       &
            (taui*(cc(m,i,3,k)+cc(m,ic,2,k))))                          &
                      -wa2(i-1)*                                        &
            ((cc(m,i,1,k)+taur*(cc(m,i,3,k)-cc(m,ic,2,k)))-             &
            (taui*(cc(m,i-1,3,k)-cc(m,ic-1,2,k))))
        ch(m,i,k,3) = wa2(i-2)*                                         &
            ((cc(m,i,1,k)+taur*(cc(m,i,3,k)-cc(m,ic,2,k)))-             &
            (taui*(cc(m,i-1,3,k)-cc(m,ic-1,2,k))))                      &
                   +wa2(i-1)*                                           &
            ((cc(m,i-1,1,k)+taur*(cc(m,i-1,3,k)+cc(m,ic-1,2,k)))+       &
            (taui*(cc(m,i,3,k)+cc(m,ic,2,k))))
      END DO
    END DO
  END DO
  RETURN
END SUBROUTINE vradb3

SUBROUTINE vradb4 (mp,ido,l1,cc,ch,mdimc,wa1,wa2,wa3)
!
!  VRFFTPK, VERSION 1, AUGUST 1985
!
  DIMENSION  cc(mdimc,ido,4,l1)  ,ch(mdimc,ido,l1,4)    ,               &
                  wa1(ido)  ,wa2(ido)  ,wa3(ido)
  sqrt2=SQRT(2.)
  DO k=1,l1
    DO m=1,mp
      ch(m,1,k,3) = (cc(m,1,1,k)+cc(m,ido,4,k))                         &
          -(cc(m,ido,2,k)+cc(m,ido,2,k))
      ch(m,1,k,1) = (cc(m,1,1,k)+cc(m,ido,4,k))                         &
          +(cc(m,ido,2,k)+cc(m,ido,2,k))
      ch(m,1,k,4) = (cc(m,1,1,k)-cc(m,ido,4,k))                         &
          +(cc(m,1,3,k)+cc(m,1,3,k))
      ch(m,1,k,2) = (cc(m,1,1,k)-cc(m,ido,4,k))                         &
          -(cc(m,1,3,k)+cc(m,1,3,k))
    END DO
  END DO
  IF (ido-2 < 0) THEN
    GO TO   107
  ELSE IF (ido-2 == 0) THEN
    GO TO   105
  END IF
  idp2 = ido+2
  DO k=1,l1
    DO i=3,ido,2
      ic = idp2-i
      DO m=1,mp
        ch(m,i-1,k,1) = (cc(m,i-1,1,k)+cc(m,ic-1,4,k))                  &
            +(cc(m,i-1,3,k)+cc(m,ic-1,2,k))
        ch(m,i,k,1) = (cc(m,i,1,k)-cc(m,ic,4,k))                        &
            +(cc(m,i,3,k)-cc(m,ic,2,k))
        ch(m,i-1,k,2)=wa1(i-2)*((cc(m,i-1,1,k)-cc(m,ic-1,4,k))          &
            -(cc(m,i,3,k)+cc(m,ic,2,k)))-wa1(i-1)                       &
            *((cc(m,i,1,k)+cc(m,ic,4,k))+(cc(m,i-1,3,k)-cc(m,ic-1,2,k)))
        ch(m,i,k,2)=wa1(i-2)*((cc(m,i,1,k)+cc(m,ic,4,k))                &
            +(cc(m,i-1,3,k)-cc(m,ic-1,2,k)))+wa1(i-1)                   &
            *((cc(m,i-1,1,k)-cc(m,ic-1,4,k))-(cc(m,i,3,k)+cc(m,ic,2,k)))
        ch(m,i-1,k,3)=wa2(i-2)*((cc(m,i-1,1,k)+cc(m,ic-1,4,k))          &
            -(cc(m,i-1,3,k)+cc(m,ic-1,2,k)))-wa2(i-1)                   &
            *((cc(m,i,1,k)-cc(m,ic,4,k))-(cc(m,i,3,k)-cc(m,ic,2,k)))
        ch(m,i,k,3)=wa2(i-2)*((cc(m,i,1,k)-cc(m,ic,4,k))                &
            -(cc(m,i,3,k)-cc(m,ic,2,k)))+wa2(i-1)                       &
            *((cc(m,i-1,1,k)+cc(m,ic-1,4,k))-(cc(m,i-1,3,k)             &
            +cc(m,ic-1,2,k)))
        ch(m,i-1,k,4)=wa3(i-2)*((cc(m,i-1,1,k)-cc(m,ic-1,4,k))          &
            +(cc(m,i,3,k)+cc(m,ic,2,k)))-wa3(i-1)                       &
            *((cc(m,i,1,k)+cc(m,ic,4,k))-(cc(m,i-1,3,k)-cc(m,ic-1,2,k)))
        ch(m,i,k,4)=wa3(i-2)*((cc(m,i,1,k)+cc(m,ic,4,k))                &
            -(cc(m,i-1,3,k)-cc(m,ic-1,2,k)))+wa3(i-1)                   &
            *((cc(m,i-1,1,k)-cc(m,ic-1,4,k))+(cc(m,i,3,k)+cc(m,ic,2,k)))
      END DO
    END DO
  END DO
  IF (MOD(ido,2) == 1) RETURN
  105 CONTINUE
  DO k=1,l1
    DO m=1,mp
      ch(m,ido,k,1) = (cc(m,ido,1,k)+cc(m,ido,3,k))                     &
          +(cc(m,ido,1,k)+cc(m,ido,3,k))
      ch(m,ido,k,2) = sqrt2*((cc(m,ido,1,k)-cc(m,ido,3,k))              &
          -(cc(m,1,2,k)+cc(m,1,4,k)))
      ch(m,ido,k,3) = (cc(m,1,4,k)-cc(m,1,2,k))                         &
          +(cc(m,1,4,k)-cc(m,1,2,k))
      ch(m,ido,k,4) = -sqrt2*((cc(m,ido,1,k)-cc(m,ido,3,k))             &
          +(cc(m,1,2,k)+cc(m,1,4,k)))
    END DO
  END DO
  107 RETURN
END SUBROUTINE vradb4

SUBROUTINE vradb5 (mp,ido,l1,cc,ch,mdimc,wa1,wa2,wa3,wa4)
!
!  VRFFTPK, VERSION 1, AUGUST 1985
!
  DIMENSION  cc(mdimc,ido,5,l1)    ,ch(mdimc,ido,l1,5),                 &
               wa1(ido)     ,wa2(ido)     ,wa3(ido)     ,wa4(ido)
  arg=2.*pimach(1.0)/5.
  tr11=COS(arg)
  ti11=SIN(arg)
  tr12=COS(2.*arg)
  ti12=SIN(2.*arg)
  DO k=1,l1
    DO m=1,mp
      ch(m,1,k,1) = cc(m,1,1,k)+2.*cc(m,ido,2,k)+2.*cc(m,ido,4,k)
      ch(m,1,k,2) = (cc(m,1,1,k)+tr11*2.*cc(m,ido,2,k)                  &
          +tr12*2.*cc(m,ido,4,k))-(ti11*2.*cc(m,1,3,k)                  &
          +ti12*2.*cc(m,1,5,k))
      ch(m,1,k,3) = (cc(m,1,1,k)+tr12*2.*cc(m,ido,2,k)                  &
          +tr11*2.*cc(m,ido,4,k))-(ti12*2.*cc(m,1,3,k)                  &
          -ti11*2.*cc(m,1,5,k))
      ch(m,1,k,4) = (cc(m,1,1,k)+tr12*2.*cc(m,ido,2,k)                  &
          +tr11*2.*cc(m,ido,4,k))+(ti12*2.*cc(m,1,3,k)                  &
          -ti11*2.*cc(m,1,5,k))
      ch(m,1,k,5) = (cc(m,1,1,k)+tr11*2.*cc(m,ido,2,k)                  &
          +tr12*2.*cc(m,ido,4,k))+(ti11*2.*cc(m,1,3,k)                  &
          +ti12*2.*cc(m,1,5,k))
    END DO
  END DO
  IF (ido == 1) RETURN
  idp2 = ido+2
  DO k=1,l1
    DO i=3,ido,2
      ic = idp2-i
      DO m=1,mp
        ch(m,i-1,k,1) = cc(m,i-1,1,k)+(cc(m,i-1,3,k)+cc(m,ic-1,2,k))    &
            +(cc(m,i-1,5,k)+cc(m,ic-1,4,k))
        ch(m,i,k,1) = cc(m,i,1,k)+(cc(m,i,3,k)-cc(m,ic,2,k))            &
            +(cc(m,i,5,k)-cc(m,ic,4,k))
        ch(m,i-1,k,2) = wa1(i-2)*((cc(m,i-1,1,k)+tr11*                  &
            (cc(m,i-1,3,k)+cc(m,ic-1,2,k))+tr12                         &
            *(cc(m,i-1,5,k)+cc(m,ic-1,4,k)))-(ti11*(cc(m,i,3,k)         &
            +cc(m,ic,2,k))+ti12*(cc(m,i,5,k)+cc(m,ic,4,k))))            &
            -wa1(i-1)*((cc(m,i,1,k)+tr11*(cc(m,i,3,k)-cc(m,ic,2,k))     &
            +tr12*(cc(m,i,5,k)-cc(m,ic,4,k)))+(ti11*(cc(m,i-1,3,k)      &
            -cc(m,ic-1,2,k))+ti12*(cc(m,i-1,5,k)-cc(m,ic-1,4,k))))
        ch(m,i,k,2) = wa1(i-2)*((cc(m,i,1,k)+tr11*(cc(m,i,3,k)          &
            -cc(m,ic,2,k))+tr12*(cc(m,i,5,k)-cc(m,ic,4,k)))             &
            +(ti11*(cc(m,i-1,3,k)-cc(m,ic-1,2,k))+ti12                  &
            *(cc(m,i-1,5,k)-cc(m,ic-1,4,k))))+wa1(i-1)                  &
            *((cc(m,i-1,1,k)+tr11*(cc(m,i-1,3,k)                        &
            +cc(m,ic-1,2,k))+tr12*(cc(m,i-1,5,k)+cc(m,ic-1,4,k)))       &
            -(ti11*(cc(m,i,3,k)+cc(m,ic,2,k))+ti12                      &
            *(cc(m,i,5,k)+cc(m,ic,4,k))))
        ch(m,i-1,k,3) = wa2(i-2)                                        &
            *((cc(m,i-1,1,k)+tr12*(cc(m,i-1,3,k)+cc(m,ic-1,2,k))        &
            +tr11*(cc(m,i-1,5,k)+cc(m,ic-1,4,k)))-(ti12*(cc(m,i,3,k)    &
            +cc(m,ic,2,k))-ti11*(cc(m,i,5,k)+cc(m,ic,4,k))))            &
            -wa2(i-1)                                                   &
            *((cc(m,i,1,k)+tr12*(cc(m,i,3,k)-                           &
            cc(m,ic,2,k))+tr11*(cc(m,i,5,k)-cc(m,ic,4,k)))              &
            +(ti12*(cc(m,i-1,3,k)-cc(m,ic-1,2,k))-ti11                  &
            *(cc(m,i-1,5,k)-cc(m,ic-1,4,k))))
        ch(m,i,k,3) = wa2(i-2)                                          &
            *((cc(m,i,1,k)+tr12*(cc(m,i,3,k)-                           &
            cc(m,ic,2,k))+tr11*(cc(m,i,5,k)-cc(m,ic,4,k)))              &
            +(ti12*(cc(m,i-1,3,k)-cc(m,ic-1,2,k))-ti11                  &
            *(cc(m,i-1,5,k)-cc(m,ic-1,4,k))))                           &
            +wa2(i-1)                                                   &
            *((cc(m,i-1,1,k)+tr12*(cc(m,i-1,3,k)+cc(m,ic-1,2,k))        &
            +tr11*(cc(m,i-1,5,k)+cc(m,ic-1,4,k)))-(ti12*(cc(m,i,3,k)    &
            +cc(m,ic,2,k))-ti11*(cc(m,i,5,k)+cc(m,ic,4,k))))
        ch(m,i-1,k,4) = wa3(i-2)                                        &
            *((cc(m,i-1,1,k)+tr12*(cc(m,i-1,3,k)+cc(m,ic-1,2,k))        &
            +tr11*(cc(m,i-1,5,k)+cc(m,ic-1,4,k)))+(ti12*(cc(m,i,3,k)    &
            +cc(m,ic,2,k))-ti11*(cc(m,i,5,k)+cc(m,ic,4,k))))            &
            -wa3(i-1)                                                   &
            *((cc(m,i,1,k)+tr12*(cc(m,i,3,k)-                           &
            cc(m,ic,2,k))+tr11*(cc(m,i,5,k)-cc(m,ic,4,k)))              &
            -(ti12*(cc(m,i-1,3,k)-cc(m,ic-1,2,k))-ti11                  &
            *(cc(m,i-1,5,k)-cc(m,ic-1,4,k))))
        ch(m,i,k,4) = wa3(i-2)                                          &
            *((cc(m,i,1,k)+tr12*(cc(m,i,3,k)-                           &
            cc(m,ic,2,k))+tr11*(cc(m,i,5,k)-cc(m,ic,4,k)))              &
            -(ti12*(cc(m,i-1,3,k)-cc(m,ic-1,2,k))-ti11                  &
            *(cc(m,i-1,5,k)-cc(m,ic-1,4,k))))                           &
            +wa3(i-1)                                                   &
            *((cc(m,i-1,1,k)+tr12*(cc(m,i-1,3,k)+cc(m,ic-1,2,k))        &
            +tr11*(cc(m,i-1,5,k)+cc(m,ic-1,4,k)))+(ti12*(cc(m,i,3,k)    &
            +cc(m,ic,2,k))-ti11*(cc(m,i,5,k)+cc(m,ic,4,k))))
        ch(m,i-1,k,5) = wa4(i-2)                                        &
            *((cc(m,i-1,1,k)+tr11*(cc(m,i-1,3,k)+cc(m,ic-1,2,k))        &
            +tr12*(cc(m,i-1,5,k)+cc(m,ic-1,4,k)))+(ti11*(cc(m,i,3,k)    &
            +cc(m,ic,2,k))+ti12*(cc(m,i,5,k)+cc(m,ic,4,k))))            &
            -wa4(i-1)                                                   &
            *((cc(m,i,1,k)+tr11*(cc(m,i,3,k)-cc(m,ic,2,k))              &
            +tr12*(cc(m,i,5,k)-cc(m,ic,4,k)))-(ti11*(cc(m,i-1,3,k)      &
            -cc(m,ic-1,2,k))+ti12*(cc(m,i-1,5,k)-cc(m,ic-1,4,k))))
        ch(m,i,k,5) = wa4(i-2)                                          &
            *((cc(m,i,1,k)+tr11*(cc(m,i,3,k)-cc(m,ic,2,k))              &
            +tr12*(cc(m,i,5,k)-cc(m,ic,4,k)))-(ti11*(cc(m,i-1,3,k)      &
            -cc(m,ic-1,2,k))+ti12*(cc(m,i-1,5,k)-cc(m,ic-1,4,k))))      &
            +wa4(i-1)                                                   &
            *((cc(m,i-1,1,k)+tr11*(cc(m,i-1,3,k)+cc(m,ic-1,2,k))        &
            +tr12*(cc(m,i-1,5,k)+cc(m,ic-1,4,k)))+(ti11*(cc(m,i,3,k)    &
            +cc(m,ic,2,k))+ti12*(cc(m,i,5,k)+cc(m,ic,4,k))))
      END DO
    END DO
  END DO
  RETURN
END SUBROUTINE vradb5

SUBROUTINE vradbg (mp,ido,ip,l1,idl1,cc,c1,c2,ch,ch2,                   &
!
!  VRFFTPK, VERSION 1, AUGUST 1985
! &
         mdimc,wa)
  DIMENSION    ch(mdimc,ido,l1,ip)    ,cc(mdimc,ido,ip,l1) ,            &
             c1(mdimc,ido,l1,ip)     ,c2(mdimc,idl1,ip),                &
                  ch2(mdimc,idl1,ip)       ,wa(ido)
  tpi=2.*pimach(1.0)
  arg = tpi/FLOAT(ip)
  dcp = COS(arg)
  dsp = SIN(arg)
  idp2 = ido+2
  nbd = (ido-1)/2
  ipp2 = ip+2
  ipph = (ip+1)/2
  IF (ido < l1) GO TO 103
  DO k=1,l1
    DO i=1,ido
      DO m=1,mp
        ch(m,i,k,1) = cc(m,i,1,k)
      END DO
    END DO
  END DO
  GO TO 106
  103 DO i=1,ido
    DO k=1,l1
      DO m=1,mp
        ch(m,i,k,1) = cc(m,i,1,k)
      END DO
    END DO
  END DO
  106 DO j=2,ipph
    jc = ipp2-j
    j2 = j+j
    DO k=1,l1
      DO m=1,mp
        ch(m,1,k,j) = cc(m,ido,j2-2,k)+cc(m,ido,j2-2,k)
        ch(m,1,k,jc) = cc(m,1,j2-1,k)+cc(m,1,j2-1,k)
      END DO
    END DO
  END DO
  IF (ido == 1) GO TO 116
  IF (nbd < l1) GO TO 112
  DO j=2,ipph
    jc = ipp2-j
    DO k=1,l1
      DO i=3,ido,2
        ic = idp2-i
        DO m=1,mp
          ch(m,i-1,k,j) = cc(m,i-1,2*j-1,k)+cc(m,ic-1,2*j-2,k)
          ch(m,i-1,k,jc) = cc(m,i-1,2*j-1,k)-cc(m,ic-1,2*j-2,k)
          ch(m,i,k,j) = cc(m,i,2*j-1,k)-cc(m,ic,2*j-2,k)
          ch(m,i,k,jc) = cc(m,i,2*j-1,k)+cc(m,ic,2*j-2,k)
        END DO
      END DO
    END DO
  END DO
  GO TO 116
  112 DO j=2,ipph
    jc = ipp2-j
    DO i=3,ido,2
      ic = idp2-i
      DO k=1,l1
        DO m=1,mp
          ch(m,i-1,k,j) = cc(m,i-1,2*j-1,k)+cc(m,ic-1,2*j-2,k)
          ch(m,i-1,k,jc) = cc(m,i-1,2*j-1,k)-cc(m,ic-1,2*j-2,k)
          ch(m,i,k,j) = cc(m,i,2*j-1,k)-cc(m,ic,2*j-2,k)
          ch(m,i,k,jc) = cc(m,i,2*j-1,k)+cc(m,ic,2*j-2,k)
        END DO
      END DO
    END DO
  END DO
  116 ar1 = 1.
  ai1 = 0.
  DO l=2,ipph
    lc = ipp2-l
    ar1h = dcp*ar1-dsp*ai1
    ai1 = dcp*ai1+dsp*ar1
    ar1 = ar1h
    DO ik=1,idl1
      DO m=1,mp
        c2(m,ik,l) = ch2(m,ik,1)+ar1*ch2(m,ik,2)
        c2(m,ik,lc) = ai1*ch2(m,ik,ip)
      END DO
    END DO
    dc2 = ar1
    ds2 = ai1
    ar2 = ar1
    ai2 = ai1
    DO j=3,ipph
      jc = ipp2-j
      ar2h = dc2*ar2-ds2*ai2
      ai2 = dc2*ai2+ds2*ar2
      ar2 = ar2h
      DO ik=1,idl1
        DO m=1,mp
          c2(m,ik,l) = c2(m,ik,l)+ar2*ch2(m,ik,j)
          c2(m,ik,lc) = c2(m,ik,lc)+ai2*ch2(m,ik,jc)
        END DO
      END DO
    END DO
  END DO
  DO j=2,ipph
    DO ik=1,idl1
      DO m=1,mp
        ch2(m,ik,1) = ch2(m,ik,1)+ch2(m,ik,j)
      END DO
    END DO
  END DO
  DO j=2,ipph
    jc = ipp2-j
    DO k=1,l1
      DO m=1,mp
        ch(m,1,k,j) = c1(m,1,k,j)-c1(m,1,k,jc)
        ch(m,1,k,jc) = c1(m,1,k,j)+c1(m,1,k,jc)
      END DO
    END DO
  END DO
  IF (ido == 1) GO TO 132
  IF (nbd < l1) GO TO 128
  DO j=2,ipph
    jc = ipp2-j
    DO k=1,l1
      DO i=3,ido,2
        DO m=1,mp
          ch(m,i-1,k,j) = c1(m,i-1,k,j)-c1(m,i,k,jc)
          ch(m,i-1,k,jc) = c1(m,i-1,k,j)+c1(m,i,k,jc)
          ch(m,i,k,j) = c1(m,i,k,j)+c1(m,i-1,k,jc)
          ch(m,i,k,jc) = c1(m,i,k,j)-c1(m,i-1,k,jc)
        END DO
      END DO
    END DO
  END DO
  GO TO 132
  128 DO j=2,ipph
    jc = ipp2-j
    DO i=3,ido,2
      DO k=1,l1
        DO m=1,mp
          ch(m,i-1,k,j) = c1(m,i-1,k,j)-c1(m,i,k,jc)
          ch(m,i-1,k,jc) = c1(m,i-1,k,j)+c1(m,i,k,jc)
          ch(m,i,k,j) = c1(m,i,k,j)+c1(m,i-1,k,jc)
          ch(m,i,k,jc) = c1(m,i,k,j)-c1(m,i-1,k,jc)
        END DO
      END DO
    END DO
  END DO
  132 CONTINUE
  IF (ido == 1) RETURN
  DO ik=1,idl1
    DO m=1,mp
      c2(m,ik,1) = ch2(m,ik,1)
    END DO
  END DO
  DO j=2,ip
    DO k=1,l1
      DO m=1,mp
        c1(m,1,k,j) = ch(m,1,k,j)
      END DO
    END DO
  END DO
  IF (nbd > l1) GO TO 139
  is = -ido
  DO j=2,ip
    is = is+ido
    idij = is
    DO i=3,ido,2
      idij = idij+2
      DO k=1,l1
        DO m=1,mp
          c1(m,i-1,k,j) = wa(idij-1)*ch(m,i-1,k,j)-wa(idij)*            &
              ch(m,i,k,j)
          c1(m,i,k,j) = wa(idij-1)*ch(m,i,k,j)+wa(idij)*                &
              ch(m,i-1,k,j)
        END DO
      END DO
    END DO
  END DO
  GO TO 143
  139 is = -ido
  DO j=2,ip
    is = is+ido
    DO k=1,l1
      idij = is
      DO i=3,ido,2
        idij = idij+2
        DO m=1,mp
          c1(m,i-1,k,j) = wa(idij-1)*ch(m,i-1,k,j)-wa(idij)*            &
              ch(m,i,k,j)
          c1(m,i,k,j) = wa(idij-1)*ch(m,i,k,j)+wa(idij)*                &
              ch(m,i-1,k,j)
        END DO
      END DO
    END DO
  END DO
  143 RETURN
END SUBROUTINE vradbg

SUBROUTINE vradf2 (mp,ido,l1,cc,ch,mdimc,wa1)
!
!  VRFFTPK, VERSION 1, AUGUST 1985
!
  DIMENSION   ch(mdimc,ido,2,l1)  ,cc(mdimc,ido,l1,2)     ,             &
                  wa1(ido)
  DO k=1,l1
    DO m=1,mp
      ch(m,1,1,k) = cc(m,1,k,1)+cc(m,1,k,2)
      ch(m,ido,2,k) = cc(m,1,k,1)-cc(m,1,k,2)
    END DO
  END DO
  IF (ido-2 < 0) THEN
    GO TO   107
  ELSE IF (ido-2 == 0) THEN
    GO TO   105
  END IF
  idp2 = ido+2
  DO k=1,l1
    DO i=3,ido,2
      ic = idp2-i
      DO m=1,mp
        ch(m,i,1,k) = cc(m,i,k,1)+(wa1(i-2)*cc(m,i,k,2)-                &
            wa1(i-1)*cc(m,i-1,k,2))
        ch(m,ic,2,k) = (wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*                  &
            cc(m,i-1,k,2))-cc(m,i,k,1)
        ch(m,i-1,1,k) = cc(m,i-1,k,1)+(wa1(i-2)*cc(m,i-1,k,2)+          &
            wa1(i-1)*cc(m,i,k,2))
        ch(m,ic-1,2,k) = cc(m,i-1,k,1)-(wa1(i-2)*cc(m,i-1,k,2)+         &
            wa1(i-1)*cc(m,i,k,2))
      END DO
    END DO
  END DO
  IF (MOD(ido,2) == 1) RETURN
  105 DO k=1,l1
    DO m=1,mp
      ch(m,1,2,k) = -cc(m,ido,k,2)
      ch(m,ido,1,k) = cc(m,ido,k,1)
    END DO
  END DO
  107 RETURN
END SUBROUTINE vradf2

SUBROUTINE vradf3 (mp,ido,l1,cc,ch,mdimc,wa1,wa2)
!
!  VRFFTPK, VERSION 1, AUGUST 1985
!
  DIMENSION   ch(mdimc,ido,3,l1)  ,cc(mdimc,ido,l1,3)     ,             &
                  wa1(ido)     ,wa2(ido)
  arg=2.*pimach(1.0)/3.
  taur=COS(arg)
  taui=SIN(arg)
  DO k=1,l1
    DO m=1,mp
      ch(m,1,1,k) = cc(m,1,k,1)+(cc(m,1,k,2)+cc(m,1,k,3))
      ch(m,1,3,k) = taui*(cc(m,1,k,3)-cc(m,1,k,2))
      ch(m,ido,2,k) = cc(m,1,k,1)+taur*                                 &
          (cc(m,1,k,2)+cc(m,1,k,3))
    END DO
  END DO
  IF (ido == 1) RETURN
  idp2 = ido+2
  DO k=1,l1
    DO i=3,ido,2
      ic = idp2-i
      DO m=1,mp
        ch(m,i-1,1,k) = cc(m,i-1,k,1)+((wa1(i-2)*cc(m,i-1,k,2)+         &
            wa1(i-1)*cc(m,i,k,2))+(wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*     &
            cc(m,i,k,3)))
        ch(m,i,1,k) = cc(m,i,k,1)+((wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*      &
            cc(m,i-1,k,2))+(wa2(i-2)*cc(m,i,k,3)-wa2(i-1)*              &
            cc(m,i-1,k,3)))
        ch(m,i-1,3,k) = (cc(m,i-1,k,1)+taur*((wa1(i-2)*                 &
            cc(m,i-1,k,2)+wa1(i-1)*cc(m,i,k,2))+(wa2(i-2)*              &
            cc(m,i-1,k,3)+wa2(i-1)*cc(m,i,k,3))))+(taui*((wa1(i-2)*     &
            cc(m,i,k,2)-wa1(i-1)*cc(m,i-1,k,2))-(wa2(i-2)*              &
            cc(m,i,k,3)-wa2(i-1)*cc(m,i-1,k,3))))
        ch(m,ic-1,2,k) = (cc(m,i-1,k,1)+taur*((wa1(i-2)*                &
            cc(m,i-1,k,2)+wa1(i-1)*cc(m,i,k,2))+(wa2(i-2)*              &
            cc(m,i-1,k,3)+wa2(i-1)*cc(m,i,k,3))))-(taui*((wa1(i-2)*     &
            cc(m,i,k,2)-wa1(i-1)*cc(m,i-1,k,2))-(wa2(i-2)*              &
            cc(m,i,k,3)-wa2(i-1)*cc(m,i-1,k,3))))
        ch(m,i,3,k) = (cc(m,i,k,1)+taur*((wa1(i-2)*cc(m,i,k,2)-         &
            wa1(i-1)*cc(m,i-1,k,2))+(wa2(i-2)*cc(m,i,k,3)-wa2(i-1)*     &
            cc(m,i-1,k,3))))+(taui*((wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*   &
            cc(m,i,k,3))-(wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*              &
            cc(m,i,k,2))))
        ch(m,ic,2,k) = (taui*((wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*         &
            cc(m,i,k,3))-(wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*              &
            cc(m,i,k,2))))-(cc(m,i,k,1)+taur*((wa1(i-2)*cc(m,i,k,2)-    &
            wa1(i-1)*cc(m,i-1,k,2))+(wa2(i-2)*cc(m,i,k,3)-wa2(i-1)*     &
            cc(m,i-1,k,3))))
      END DO
    END DO
  END DO
  RETURN
END SUBROUTINE vradf3

SUBROUTINE vradf4 (mp,ido,l1,cc,ch,mdimc,wa1,wa2,wa3)
!
!  VRFFTPK, VERSION 1, AUGUST 1985
!
  DIMENSION    cc(mdimc,ido,l1,4)   ,ch(mdimc,ido,4,l1)     ,           &
                  wa1(ido)     ,wa2(ido)     ,wa3(ido)
  hsqt2=SQRT(2.)/2.
  DO k=1,l1
    DO m=1,mp
      ch(m,1,1,k) = (cc(m,1,k,2)+cc(m,1,k,4))                           &
          +(cc(m,1,k,1)+cc(m,1,k,3))
      ch(m,ido,4,k) = (cc(m,1,k,1)+cc(m,1,k,3))                         &
          -(cc(m,1,k,2)+cc(m,1,k,4))
      ch(m,ido,2,k) = cc(m,1,k,1)-cc(m,1,k,3)
      ch(m,1,3,k) = cc(m,1,k,4)-cc(m,1,k,2)
    END DO
  END DO
  IF (ido-2 < 0) THEN
    GO TO   107
  ELSE IF (ido-2 == 0) THEN
    GO TO   105
  END IF
  idp2 = ido+2
  DO k=1,l1
    DO i=3,ido,2
      ic = idp2-i
      DO m=1,mp
        ch(m,i-1,1,k) = ((wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*              &
            cc(m,i,k,2))+(wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)*              &
            cc(m,i,k,4)))+(cc(m,i-1,k,1)+(wa2(i-2)*cc(m,i-1,k,3)+       &
            wa2(i-1)*cc(m,i,k,3)))
        ch(m,ic-1,4,k) = (cc(m,i-1,k,1)+(wa2(i-2)*cc(m,i-1,k,3)+        &
            wa2(i-1)*cc(m,i,k,3)))-((wa1(i-2)*cc(m,i-1,k,2)+            &
            wa1(i-1)*cc(m,i,k,2))+(wa3(i-2)*cc(m,i-1,k,4)+              &
            wa3(i-1)*cc(m,i,k,4)))
        ch(m,i,1,k) = ((wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*                  &
            cc(m,i-1,k,2))+(wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*              &
            cc(m,i-1,k,4)))+(cc(m,i,k,1)+(wa2(i-2)*cc(m,i,k,3)-         &
            wa2(i-1)*cc(m,i-1,k,3)))
        ch(m,ic,4,k) = ((wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*                 &
            cc(m,i-1,k,2))+(wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*              &
            cc(m,i-1,k,4)))-(cc(m,i,k,1)+(wa2(i-2)*cc(m,i,k,3)-         &
            wa2(i-1)*cc(m,i-1,k,3)))
        ch(m,i-1,3,k) = ((wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*                &
            cc(m,i-1,k,2))-(wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*              &
            cc(m,i-1,k,4)))+(cc(m,i-1,k,1)-(wa2(i-2)*cc(m,i-1,k,3)+     &
            wa2(i-1)*cc(m,i,k,3)))
        ch(m,ic-1,2,k) = (cc(m,i-1,k,1)-(wa2(i-2)*cc(m,i-1,k,3)+        &
            wa2(i-1)*cc(m,i,k,3)))-((wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*     &
            cc(m,i-1,k,2))-(wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*              &
            cc(m,i-1,k,4)))
        ch(m,i,3,k) = ((wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)*                &
            cc(m,i,k,4))-(wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*              &
            cc(m,i,k,2)))+(cc(m,i,k,1)-(wa2(i-2)*cc(m,i,k,3)-           &
            wa2(i-1)*cc(m,i-1,k,3)))
        ch(m,ic,2,k) = ((wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)*               &
            cc(m,i,k,4))-(wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*              &
            cc(m,i,k,2)))-(cc(m,i,k,1)-(wa2(i-2)*cc(m,i,k,3)-wa2(i-1)*  &
            cc(m,i-1,k,3)))
      END DO
    END DO
  END DO
  IF (MOD(ido,2) == 1) RETURN
  105 CONTINUE
  DO k=1,l1
    DO m=1,mp
      ch(m,ido,1,k) = (hsqt2*(cc(m,ido,k,2)-cc(m,ido,k,4)))+            &
          cc(m,ido,k,1)
      ch(m,ido,3,k) = cc(m,ido,k,1)-(hsqt2*(cc(m,ido,k,2)-              &
          cc(m,ido,k,4)))
      ch(m,1,2,k) = (-hsqt2*(cc(m,ido,k,2)+cc(m,ido,k,4)))-             &
          cc(m,ido,k,3)
      ch(m,1,4,k) = (-hsqt2*(cc(m,ido,k,2)+cc(m,ido,k,4)))+             &
          cc(m,ido,k,3)
    END DO
  END DO
  107 RETURN
END SUBROUTINE vradf4

SUBROUTINE vradf5 (mp,ido,l1,cc,ch,mdimc,wa1,wa2,wa3,wa4)
!
!  VRFFTPK, VERSION 1, AUGUST 1985
!
  DIMENSION  cc(mdimc,ido,l1,5)    ,ch(mdimc,ido,5,l1)     ,            &
             wa1(ido)     ,wa2(ido)     ,wa3(ido)     ,wa4(ido)
  arg=2.*pimach(1.0)/5.
  tr11=COS(arg)
  ti11=SIN(arg)
  tr12=COS(2.*arg)
  ti12=SIN(2.*arg)
  DO k=1,l1
    DO m=1,mp
      ch(m,1,1,k) = cc(m,1,k,1)+(cc(m,1,k,5)+cc(m,1,k,2))+              &
          (cc(m,1,k,4)+cc(m,1,k,3))
      ch(m,ido,2,k) = cc(m,1,k,1)+tr11*(cc(m,1,k,5)+cc(m,1,k,2))+       &
          tr12*(cc(m,1,k,4)+cc(m,1,k,3))
      ch(m,1,3,k) = ti11*(cc(m,1,k,5)-cc(m,1,k,2))+ti12*                &
          (cc(m,1,k,4)-cc(m,1,k,3))
      ch(m,ido,4,k) = cc(m,1,k,1)+tr12*(cc(m,1,k,5)+cc(m,1,k,2))+       &
          tr11*(cc(m,1,k,4)+cc(m,1,k,3))
      ch(m,1,5,k) = ti12*(cc(m,1,k,5)-cc(m,1,k,2))-ti11*                &
          (cc(m,1,k,4)-cc(m,1,k,3))
    END DO
  END DO
  IF (ido == 1) RETURN
  idp2 = ido+2
  DO k=1,l1
    DO i=3,ido,2
      ic = idp2-i
      DO m=1,mp
        ch(m,i-1,1,k) = cc(m,i-1,k,1)+((wa1(i-2)*cc(m,i-1,k,2)+         &
            wa1(i-1)*cc(m,i,k,2))+(wa4(i-2)*cc(m,i-1,k,5)+wa4(i-1)*     &
            cc(m,i,k,5)))+((wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*            &
            cc(m,i,k,3))+(wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)*cc(m,i,k,4)))
        ch(m,i,1,k) = cc(m,i,k,1)+((wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*      &
            cc(m,i-1,k,2))+(wa4(i-2)*cc(m,i,k,5)-wa4(i-1)*              &
            cc(m,i-1,k,5)))+((wa2(i-2)*cc(m,i,k,3)-wa2(i-1)*            &
            cc(m,i-1,k,3))+(wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*              &
            cc(m,i-1,k,4)))
        ch(m,i-1,3,k) = cc(m,i-1,k,1)+tr11*                             &
            ( wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*cc(m,i,k,2)               &
            +wa4(i-2)*cc(m,i-1,k,5)+wa4(i-1)*cc(m,i,k,5))+tr12*         &
            ( wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*cc(m,i,k,3)               &
            +wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)*cc(m,i,k,4))+ti11*         &
            ( wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*cc(m,i-1,k,2)               &
            -(wa4(i-2)*cc(m,i,k,5)-wa4(i-1)*cc(m,i-1,k,5)))+ti12*       &
            ( wa2(i-2)*cc(m,i,k,3)-wa2(i-1)*cc(m,i-1,k,3)               &
            -(wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*cc(m,i-1,k,4)))
        ch(m,ic-1,2,k) = cc(m,i-1,k,1)+tr11*                            &
            ( wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*cc(m,i,k,2)               &
            +wa4(i-2)*cc(m,i-1,k,5)+wa4(i-1)*cc(m,i,k,5))+tr12*         &
            ( wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*cc(m,i,k,3)               &
            +wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)*cc(m,i,k,4))-(ti11*        &
            ( wa1(i-2)*cc(m,i,k,2)-wa1(i-1)*cc(m,i-1,k,2)               &
            -(wa4(i-2)*cc(m,i,k,5)-wa4(i-1)*cc(m,i-1,k,5)))+ti12*       &
            ( wa2(i-2)*cc(m,i,k,3)-wa2(i-1)*cc(m,i-1,k,3)               &
            -(wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*cc(m,i-1,k,4))))
        ch(m,i,3,k) = (cc(m,i,k,1)+tr11*((wa1(i-2)*cc(m,i,k,2)-         &
            wa1(i-1)*cc(m,i-1,k,2))+(wa4(i-2)*cc(m,i,k,5)-wa4(i-1)*     &
            cc(m,i-1,k,5)))+tr12*((wa2(i-2)*cc(m,i,k,3)-wa2(i-1)*       &
            cc(m,i-1,k,3))+(wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*              &
            cc(m,i-1,k,4))))+(ti11*((wa4(i-2)*cc(m,i-1,k,5)+            &
            wa4(i-1)*cc(m,i,k,5))-(wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*     &
            cc(m,i,k,2)))+ti12*((wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)*       &
            cc(m,i,k,4))-(wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*              &
            cc(m,i,k,3))))
        ch(m,ic,2,k) = (ti11*((wa4(i-2)*cc(m,i-1,k,5)+wa4(i-1)*         &
            cc(m,i,k,5))-(wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*              &
            cc(m,i,k,2)))+ti12*((wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)*       &
            cc(m,i,k,4))-(wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*              &
            cc(m,i,k,3))))-(cc(m,i,k,1)+tr11*((wa1(i-2)*cc(m,i,k,2)-    &
            wa1(i-1)*cc(m,i-1,k,2))+(wa4(i-2)*cc(m,i,k,5)-wa4(i-1)*     &
            cc(m,i-1,k,5)))+tr12*((wa2(i-2)*cc(m,i,k,3)-wa2(i-1)*       &
            cc(m,i-1,k,3))+(wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*              &
            cc(m,i-1,k,4))))
        ch(m,i-1,5,k) = (cc(m,i-1,k,1)+tr12*((wa1(i-2)*                 &
            cc(m,i-1,k,2)+wa1(i-1)*cc(m,i,k,2))+(wa4(i-2)*              &
            cc(m,i-1,k,5)+wa4(i-1)*cc(m,i,k,5)))+tr11*((wa2(i-2)*       &
            cc(m,i-1,k,3)+wa2(i-1)*cc(m,i,k,3))+(wa3(i-2)*              &
            cc(m,i-1,k,4)+wa3(i-1)*cc(m,i,k,4))))+(ti12*((wa1(i-2)*     &
            cc(m,i,k,2)-wa1(i-1)*cc(m,i-1,k,2))-(wa4(i-2)*cc(m,i,k,5)-  &
            wa4(i-1)*cc(m,i-1,k,5)))-ti11*((wa2(i-2)*cc(m,i,k,3)-       &
            wa2(i-1)*cc(m,i-1,k,3))-(wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*     &
            cc(m,i-1,k,4))))
        ch(m,ic-1,4,k) = (cc(m,i-1,k,1)+tr12*((wa1(i-2)*                &
            cc(m,i-1,k,2)+wa1(i-1)*cc(m,i,k,2))+(wa4(i-2)*              &
            cc(m,i-1,k,5)+wa4(i-1)*cc(m,i,k,5)))+tr11*((wa2(i-2)*       &
            cc(m,i-1,k,3)+wa2(i-1)*cc(m,i,k,3))+(wa3(i-2)*              &
            cc(m,i-1,k,4)+wa3(i-1)*cc(m,i,k,4))))-(ti12*((wa1(i-2)*     &
            cc(m,i,k,2)-wa1(i-1)*cc(m,i-1,k,2))-(wa4(i-2)*cc(m,i,k,5)-  &
            wa4(i-1)*cc(m,i-1,k,5)))-ti11*((wa2(i-2)*cc(m,i,k,3)-       &
            wa2(i-1)*cc(m,i-1,k,3))-(wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*     &
            cc(m,i-1,k,4))))
        ch(m,i,5,k) = (cc(m,i,k,1)+tr12*((wa1(i-2)*cc(m,i,k,2)-         &
            wa1(i-1)*cc(m,i-1,k,2))+(wa4(i-2)*cc(m,i,k,5)-wa4(i-1)*     &
            cc(m,i-1,k,5)))+tr11*((wa2(i-2)*cc(m,i,k,3)-wa2(i-1)*       &
            cc(m,i-1,k,3))+(wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*              &
            cc(m,i-1,k,4))))+(ti12*((wa4(i-2)*cc(m,i-1,k,5)+            &
            wa4(i-1)*cc(m,i,k,5))-(wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*     &
            cc(m,i,k,2)))-ti11*((wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)*       &
            cc(m,i,k,4))-(wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*              &
            cc(m,i,k,3))))
        ch(m,ic,4,k) = (ti12*((wa4(i-2)*cc(m,i-1,k,5)+wa4(i-1)*         &
            cc(m,i,k,5))-(wa1(i-2)*cc(m,i-1,k,2)+wa1(i-1)*              &
            cc(m,i,k,2)))-ti11*((wa3(i-2)*cc(m,i-1,k,4)+wa3(i-1)*       &
            cc(m,i,k,4))-(wa2(i-2)*cc(m,i-1,k,3)+wa2(i-1)*              &
            cc(m,i,k,3))))-(cc(m,i,k,1)+tr12*((wa1(i-2)*cc(m,i,k,2)-    &
            wa1(i-1)*cc(m,i-1,k,2))+(wa4(i-2)*cc(m,i,k,5)-wa4(i-1)*     &
            cc(m,i-1,k,5)))+tr11*((wa2(i-2)*cc(m,i,k,3)-wa2(i-1)*       &
            cc(m,i-1,k,3))+(wa3(i-2)*cc(m,i,k,4)-wa3(i-1)*              &
            cc(m,i-1,k,4))))
      END DO
    END DO
  END DO
  RETURN
END SUBROUTINE vradf5

SUBROUTINE vradfg (mp,ido,ip,l1,idl1,cc,c1,c2,ch,ch2,mdimc,wa)
!
!  VRFFTPK, VERSION 1, AUGUST 1985
!
  DIMENSION     ch(mdimc,ido,l1,ip)   ,cc(mdimc,ido,ip,l1)  ,           &
              c1(mdimc,ido,l1,ip)    ,c2(mdimc,idl1,ip),                &
                  ch2(mdimc,idl1,ip)           ,wa(ido)
  tpi=2.*pimach(1.0)
  arg = tpi/FLOAT(ip)
  dcp = COS(arg)
  dsp = SIN(arg)
  ipph = (ip+1)/2
  ipp2 = ip+2
  idp2 = ido+2
  nbd = (ido-1)/2
  IF (ido == 1) GO TO 119
  DO ik=1,idl1
    DO m=1,mp
      ch2(m,ik,1) = c2(m,ik,1)
    END DO
  END DO
  DO j=2,ip
    DO k=1,l1
      DO m=1,mp
        ch(m,1,k,j) = c1(m,1,k,j)
      END DO
    END DO
  END DO
  IF (nbd > l1) GO TO 107
  is = -ido
  DO j=2,ip
    is = is+ido
    idij = is
    DO i=3,ido,2
      idij = idij+2
      DO k=1,l1
        DO m=1,mp
          ch(m,i-1,k,j) = wa(idij-1)*c1(m,i-1,k,j)+wa(idij)             &
              *c1(m,i,k,j)
          ch(m,i,k,j) = wa(idij-1)*c1(m,i,k,j)-wa(idij)                 &
              *c1(m,i-1,k,j)
        END DO
      END DO
    END DO
  END DO
  GO TO 111
  107 is = -ido
  DO j=2,ip
    is = is+ido
    DO k=1,l1
      idij = is
      DO i=3,ido,2
        idij = idij+2
        DO m=1,mp
          ch(m,i-1,k,j) = wa(idij-1)*c1(m,i-1,k,j)+wa(idij)             &
              *c1(m,i,k,j)
          ch(m,i,k,j) = wa(idij-1)*c1(m,i,k,j)-wa(idij)                 &
              *c1(m,i-1,k,j)
        END DO
      END DO
    END DO
  END DO
  111 IF (nbd < l1) GO TO 115
  DO j=2,ipph
    jc = ipp2-j
    DO k=1,l1
      DO i=3,ido,2
        DO m=1,mp
          c1(m,i-1,k,j) = ch(m,i-1,k,j)+ch(m,i-1,k,jc)
          c1(m,i-1,k,jc) = ch(m,i,k,j)-ch(m,i,k,jc)
          c1(m,i,k,j) = ch(m,i,k,j)+ch(m,i,k,jc)
          c1(m,i,k,jc) = ch(m,i-1,k,jc)-ch(m,i-1,k,j)
        END DO
      END DO
    END DO
  END DO
  GO TO 121
  115 DO j=2,ipph
    jc = ipp2-j
    DO i=3,ido,2
      DO k=1,l1
        DO m=1,mp
          c1(m,i-1,k,j) = ch(m,i-1,k,j)+ch(m,i-1,k,jc)
          c1(m,i-1,k,jc) = ch(m,i,k,j)-ch(m,i,k,jc)
          c1(m,i,k,j) = ch(m,i,k,j)+ch(m,i,k,jc)
          c1(m,i,k,jc) = ch(m,i-1,k,jc)-ch(m,i-1,k,j)
        END DO
      END DO
    END DO
  END DO
  GO TO 121
  119 DO ik=1,idl1
    DO m=1,mp
      c2(m,ik,1) = ch2(m,ik,1)
    END DO
  END DO
  121 DO j=2,ipph
    jc = ipp2-j
    DO k=1,l1
      DO m=1,mp
        c1(m,1,k,j) = ch(m,1,k,j)+ch(m,1,k,jc)
        c1(m,1,k,jc) = ch(m,1,k,jc)-ch(m,1,k,j)
      END DO
    END DO
  END DO
!
  ar1 = 1.
  ai1 = 0.
  DO l=2,ipph
    lc = ipp2-l
    ar1h = dcp*ar1-dsp*ai1
    ai1 = dcp*ai1+dsp*ar1
    ar1 = ar1h
    DO ik=1,idl1
      DO m=1,mp
        ch2(m,ik,l) = c2(m,ik,1)+ar1*c2(m,ik,2)
        ch2(m,ik,lc) = ai1*c2(m,ik,ip)
      END DO
    END DO
    dc2 = ar1
    ds2 = ai1
    ar2 = ar1
    ai2 = ai1
    DO j=3,ipph
      jc = ipp2-j
      ar2h = dc2*ar2-ds2*ai2
      ai2 = dc2*ai2+ds2*ar2
      ar2 = ar2h
      DO ik=1,idl1
        DO m=1,mp
          ch2(m,ik,l) = ch2(m,ik,l)+ar2*c2(m,ik,j)
          ch2(m,ik,lc) = ch2(m,ik,lc)+ai2*c2(m,ik,jc)
        END DO
      END DO
    END DO
  END DO
  DO j=2,ipph
    DO ik=1,idl1
      DO m=1,mp
        ch2(m,ik,1) = ch2(m,ik,1)+c2(m,ik,j)
      END DO
    END DO
  END DO
!
  IF (ido < l1) GO TO 132
  DO k=1,l1
    DO i=1,ido
      DO m=1,mp
        cc(m,i,1,k) = ch(m,i,k,1)
      END DO
    END DO
  END DO
  GO TO 135
  132 DO i=1,ido
    DO k=1,l1
      DO m=1,mp
        cc(m,i,1,k) = ch(m,i,k,1)
      END DO
    END DO
  END DO
  135 DO j=2,ipph
    jc = ipp2-j
    j2 = j+j
    DO k=1,l1
      DO m=1,mp
        cc(m,ido,j2-2,k) = ch(m,1,k,j)
        cc(m,1,j2-1,k) = ch(m,1,k,jc)
      END DO
    END DO
  END DO
  IF (ido == 1) RETURN
  IF (nbd < l1) GO TO 141
  DO j=2,ipph
    jc = ipp2-j
    j2 = j+j
    DO k=1,l1
      DO i=3,ido,2
        ic = idp2-i
        DO m=1,mp
          cc(m,i-1,j2-1,k) = ch(m,i-1,k,j)+ch(m,i-1,k,jc)
          cc(m,ic-1,j2-2,k) = ch(m,i-1,k,j)-ch(m,i-1,k,jc)
          cc(m,i,j2-1,k) = ch(m,i,k,j)+ch(m,i,k,jc)
          cc(m,ic,j2-2,k) = ch(m,i,k,jc)-ch(m,i,k,j)
        END DO
      END DO
    END DO
  END DO
  RETURN
  141 DO j=2,ipph
    jc = ipp2-j
    j2 = j+j
    DO i=3,ido,2
      ic = idp2-i
      DO k=1,l1
        DO m=1,mp
          cc(m,i-1,j2-1,k) = ch(m,i-1,k,j)+ch(m,i-1,k,jc)
          cc(m,ic-1,j2-2,k) = ch(m,i-1,k,j)-ch(m,i-1,k,jc)
          cc(m,i,j2-1,k) = ch(m,i,k,j)+ch(m,i,k,jc)
          cc(m,ic,j2-2,k) = ch(m,i,k,jc)-ch(m,i,k,j)
        END DO
      END DO
    END DO
  END DO
  RETURN
END SUBROUTINE vradfg

SUBROUTINE vrfftb(m,n,r,rt,mdimr,wsave)
!***BEGIN PROLOGUE  VRFFTB
!***DATE WRITTEN   850801   (YYMMDD)
!***REVISION DATE  900509   (YYMMDD)
!***CATEGORY NO.  J1A1
!***KEYWORDS  FAST FOURIER TRANSFORM, REAL PERIODIC TRANSFORM,
!          FOURIER SYNTHESIS, BACKWARD TRANSFORM, MULTIPLE SEQUENCES
!***AUTHOR  SWEET, R.A. (NIST) AND LINDGREN, L.L. (NIST)
!***PURPOSE  Backward real periodic transform, M sequences.
!***DESCRIPTION
!
!  Subroutine VRFFTB computes the synthesis (backward transform) of a
!  number of real periodic sequences from their Fourier coefficients.
!  Specifically, for each set of independent Fourier coefficients
!  F(K), the corresponding real periodic sequence is computed.
!
!  The array WSAVE which is used by subroutine VRFFTB must be
!  initialized by calling subroutine VRFFTI(N,WSAVE).
!
!
!  Input Parameters
!
!  M       the number of sets of coefficients.
!
!  N       the length of the sequences of coefficients to be
!       transformed.  The method is most efficient when N is a
!       product of small primes, however n may be any positive
!       integer.
!
!  R       areal two-dimensional array of size MDIMX x N containing the
!       coefficients to be transformed.  Each set of coefficients
!       F(K), K\0,1,..,N-1, is stored as a ROW of R.  Specifically,
!       the I-th set of independent Fourier coefficients is stored
!
!             R(I,1) = REAL( F(I,0) ),
!
!             R(I,2*K) = REAL( F(I,K) )
!
!             R(I,2*K+1) = IMAG( F(I,K) )
!
!                for K = 1, 2, . . . , M-1,
!
!             and, when N is even,
!
!             R(I,N) = REAL( F(I,N/2) ).
!
!  RT      a real two-dimensional work array of size MDIMX x N.
!
!  MDIMR   the row (or first) dimension of the arrays R and RT exactly
!       as they appear in the calling program.  This parameter is
!       used to specify the variable dimension of these arrays.
!
!  WSAVE   a real one-dimensional work array which must be dimensioned
!       at least N+15.  The WSAVE array must be initialized by
!       calling subroutine VRFFTI.  A different WSAVE array must be
!       used for each different value of N.  This initialization does
!       not have to be repeated so long as N remains unchanged.  The
!       same WSAVE array may be used by VRFFTB and VRFFTB.
!
!  Output Parameters
!
!  R       contains M real periodic sequences corresponding to the given
!       coefficients.  Specifically, the I-th row of R contains the
!       real periodic sequence corresponding to the I-th set of
!       independent Fourier coefficients F(I,K) stored as
!
!            R(I,J) = X(I,J-1) ,   J = 1, 2, . . . , N, where
!
!            X(I,J) = SQRT(1/N)* F(I,0) + (-1)**J*F(I,N/2)
!                     + 2*SUM(K=1,M)[ REAL(F(I,2K))*COS(2K*J*PI/N)
!                     - IMAG(F(I,2K+1))*SIN(2K*J*PI/N) ]  ,
!
!              when N is even, and
!
!            X(I,J) = SQRT(1/N)* F(I,0) +
!                     2*SUM(K=1,M)[ REAL(F(I,2K))*COS(2K*J*PI/N)
!                     - IMAG(F(I,2K+1))*SIN(2K*J*PI/N) ]  ,
!
!              when N is odd.
!
!  WSAVE   contains results which must not be destroyed between calls
!       to VRFFTF or VRFFTB.
!
!  -----------------------------------------------------------------
!
!  NOTE  -  A call of VRFFTF followed immediately by a call of
!        of VRFFTB will return the original sequences R.  Thus,
!        VRFFTB is the correctly normalized inverse of VRFFTF.
!
!  -----------------------------------------------------------------
!
!  VRFFTB is a straightforward extension of the subprogram RFFTB to
!  handle M simultaneous sequences.  RFFTB was originally developed
!  by P. N. Swarztrauber of NCAR.
!
!
!           * * * * * * * * * * * * * * * * * * * * *
!           *                                       *
!           *         PROGRAM SPECIFICATIONS        *
!           *                                       *
!           * * * * * * * * * * * * * * * * * * * * *
!
!
!  DIMENSION OF    R(MDIMR,N), RT(MDIMR,N), WSAVE(N+15)
!  ARGUMENTS
!
!  LATEST          AUGUST 1, 1985
!  REVISION
!
!  SUBPROGRAMS     VRFFTI, VRFTI1, VRFFTF, VRFTF1, VRADF2, VRADF3,
!  REQUIRED        VRADF4, VRADF5, VRADFG, VRFFTB, VRFTB1, VRADB2,
!                  VRADB3, VRADB4, VRADB5, VRADBG, PIMACH
!
!  SPECIAL         NONE
!  CONDITIONS
!
!  COMMON          NONE
!  BLOCKS
!
!  I/O             NONE
!
!  PRECISION       SINGLE
!
!  SPECIALIST      ROLAND SWEET
!
!  LANGUAGE        FORTRAN
!
!  HISTORY         WRITTEN BY LINDA LINDGREN AND ROLAND SWEET AT THE
!                  NATIONAL BUREAU OF STANDARDS (BOULDER).
!
!  ALGORITHM       A REAL VARIANT OF THE STOCKHAM AUTOSORT VERSION
!                  OF THE COOLEY-TUKEY FAST FOURIER TRANSFORM.
!
!  PORTABILITY     AMERICAN NATIONAL STANDARDS INSTITUTE FORTRAN 77.
!                  THE ONLY MACHINE DEPENDENT CONSTANT IS LOCATED IN
!                  THE FUNCTION PIMACH.
!
!  REQUIRED        COS,SIN
!  RESIDENT
!  ROUTINES
!
!
!***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
!            Computations, (G. Rodrigue, ed.), Academic Press, 1982,
!            pp. 51-83.
!***ROUTINES CALLED  VRFTB1
!***END PROLOGUE  VRFFTB
!
!  VRFFTPK, VERSION 1, AUGUST 1985
!
  DIMENSION     r(mdimr,n),rt(mdimr,n),wsave(n+15)
  IF (n == 1) RETURN
  CALL vrftb1 (m,n,r,rt,mdimr,wsave(1),wsave(n+1))
  RETURN
END SUBROUTINE vrfftb

SUBROUTINE vrfftf (m,n,r,rt,mdimr,wsave)
!***BEGIN PROLOGUE  VRFFTF
!***DATE WRITTEN   850801   (YYMMDD)
!***REVISION DATE  900509   (YYMMDD)
!***CATEGORY NO.  J1A1
!***KEYWORDS  FAST FOURIER TRANSFORM, REAL PERIODIC TRANSFORM,
!          FOURIER ANALYSIS, FORWARD TRANSFORM, MULTIPLE SEQUENCES
!***AUTHOR  SWEET, R.A. (NIST) AND LINDGREN, L.L. (NIST)
!***PURPOSE  Forward real periodic transform, M sequences.
!***DESCRIPTION
!
!  Subroutine VRFFTF computes the Fourier coefficients (forward
!  transform) of a number of real periodic sequences.  Specifically,
!  for each sequence the subroutine claculates the independent
!  Fourier coefficients described below at output parameter R.
!
!  The array WSAVE which is used by subroutine VRFFTF must be
!  initialized by calling subroutine VRFFTI(N,WSAVE).
!
!
!  Input Parameters
!
!  M       the number of sequences to be transformed.
!
!  N       the length of the sequences to be transformed.  The method
!       is most efficient when N is a product of small primes,
!       however n may be any positive integer.
!
!  R       areal two-dimensional array of size MDIMX x N containing the
!       the sequences to be transformed.  The sequences are stored
!       in the ROWS of R.  Thus, the I-th sequence to be transformed,
!       X(I,J), J=0,1,...,N-1, is stored as
!
!            R(I,J) = X(I,J-1) , J=1, 2, . . . , N.
!
!  RT      a real two-dimensional work array of size MDIMX x N.
!
!  MDIMR   the row (or first) dimension of the arrays R and RT exactly
!       as they appear in the calling program.  This parameter is
!       used to specify the variable dimension of these arrays.
!
!  WSAVE   a real one-dimensional work array which must be dimensioned
!       at least N+15.  The WSAVE array must be initialized by
!       calling subroutine VRFFTI.  A different WSAVE array must be
!       used for each different value of N.  This initialization does
!       not have to be repeated so long as N remains unchanged.  The
!       same WSAVE array may be used by VRFFTF and VRFFTB.
!
!  Output Parameters
!
!  R       contains the Fourier coefficients F(K) for each of the M
!       input sequences.  Specifically, row I of R, R(I,J),
!       J=1,2,..,N, contains the independent Fourier coefficients
!       F(I,K), for the I-th input sequence stored as
!
!          R(I,1) = REAL( F(I,0) ),
!                 = SQRT(1/N)*SUM(J=0,N-1)[ X(I,J) ],
!
!          R(I,2*K) = REAL( F(I,K) )
!                   = SQRT(1/N)*SUM(J=0,N-1)[X(I,J)*COS(2J*K*PI/N)]
!
!          R(I,2*K+1) = IMAG( F(I,K) )
!                     =-SQRT(1/N)*SUM(J=0,N-1)[X(I,J)*SIN(2J*K*PI/N)]
!
!                for K = 1, 2, . . . , M-1,
!
!           and, when N is even,
!
!           R(I,N) = REAL( F(I,N/2) ).
!                  = SQRT(1/N)*SUM(J=0,N-1)[ (-1)**J*X(I,J) ].
!
!  WSAVE   contains results which must not be destroyed between calls
!       to VRFFTF or VRFFTB.
!
!  -----------------------------------------------------------------
!
!  NOTE  -  A call of VRFFTF followed immediately by a call of
!        of VRFFTB will return the original sequences R.  Thus,
!        VRFFTB is the correctly normalized inverse of VRFFTF.
!
!  -----------------------------------------------------------------
!
!  VRFFTF is a straightforward extension of the subprogram RFFTF to
!  handle M simultaneous sequences.  RFFTF was originally developed
!  by P. N. Swarztrauber of NCAR.
!
!
!           * * * * * * * * * * * * * * * * * * * * *
!           *                                       *
!           *         PROGRAM SPECIFICATIONS        *
!           *                                       *
!           * * * * * * * * * * * * * * * * * * * * *
!
!
!  DIMENSION OF    R(MDIMR,N), RT(MDIMR,N), WSAVE(N+15)
!  ARGUMENTS
!
!  LATEST          AUGUST 1, 1985
!  REVISION
!
!  SUBPROGRAMS     VRFFTI, VRFTI1, VRFFTF, VRFTF1, VRADF2, VRADF3,
!  REQUIRED        VRADF4, VRADF5, VRADFG, VRFFTB, VRFTB1, VRADB2,
!                  VRADB3, VRADB4, VRADB5, VRADBG, PIMACH
!
!  SPECIAL         NONE
!  CONDITIONS
!
!  COMMON          NONE
!  BLOCKS
!
!  I/O             NONE
!
!  PRECISION       SINGLE
!
!  SPECIALIST      ROLAND SWEET
!
!  LANGUAGE        FORTRAN
!
!  HISTORY         WRITTEN BY LINDA LINDGREN AND ROLAND SWEET AT THE
!                  NATIONAL BUREAU OF STANDARDS (BOULDER).
!
!  ALGORITHM       A REAL VARIANT OF THE STOCKHAM AUTOSORT VERSION
!                  OF THE COOLEY-TUKEY FAST FOURIER TRANSFORM.
!
!  PORTABILITY     AMERICAN NATIONAL STANDARDS INSTITUTE FORTRAN 77.
!                  THE ONLY MACHINE DEPENDENT CONSTANT IS LOCATED IN
!                  THE FUNCTION PIMACH.
!
!  REQUIRED        COS,SIN
!  RESIDENT
!  ROUTINES
!
!
!***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
!            Computations, (G. Rodrigue, ed.), Academic Press, 1982,
!            pp. 51-83.
!***ROUTINES CALLED  VRFTF1
!***END PROLOGUE  VRFFTF
!
!  VRFFTPK, VERSION 1, AUGUST 1985
!
  DIMENSION       r(mdimr,n)  ,rt(mdimr,n)    ,wsave(n+15)
!***FIRST EXECUTABLE STATEMENT  VRFFTF
  IF (n == 1) RETURN
  CALL vrftf1 (m,n,r,rt,mdimr,wsave(1),wsave(n+1))
  RETURN
END SUBROUTINE vrfftf

SUBROUTINE vrffti (n,wsave)
!***BEGIN PROLOGUE  VRFFTI
!***DATE WRITTEN   860701   (YYMMDD)
!***REVISION DATE  900509   (YYMMDD)
!***CATEGORY NO.  J1A1
!***KEYWORDS  FAST FOURIER TRANSFORM, REAL PERIODIC TRANSFORM,
!          MULTIPLE SEQUENCES
!***AUTHOR  SWEET, R.A. (NIST) AND LINDGREN, L.L. (NIST)
!***PURPOSE  Initialization for VRFFTF and VRFFTB.
!***DESCRIPTION
!
!  Subroutine VRFFTI initializes the array WSAVE which is used in
!  both VRFFTF and VRFFTB.  The prime factorization of N together with
!  a tabulation of certain trigonometric functions are computed and
!  stored in the array WSAVE.
!
!  Input Parameter
!
!  N       the length of the sequence to be transformed.  There is no
!       restriction on N.
!
!  Output Parameter
!
!  WSAVE   a work array which must be dimensioned at least N+15.
!       The same work array can be used for both VRFFTF and VRFFTB
!       as long as N remains unchanged.  Different WSAVE arrays
!       are required for different values of N.  The contents of
!       WSAVE must not be changed between calls of VRFFTF or VRFFTB.
!
!
!           * * * * * * * * * * * * * * * * * * * * *
!           *                                       *
!           *         PROGRAM SPECIFICATIONS        *
!           *                                       *
!           * * * * * * * * * * * * * * * * * * * * *
!
!
!  DIMENSION OF    R(MDIMR,N), RT(MDIMR,N), WSAVE(N+15)
!  ARGUMENTS
!
!  LATEST          AUGUST 1, 1985
!  REVISION
!
!  SUBPROGRAMS     VRFFTI, VRFTI1, VRFFTF, VRFTF1, VRADF2, VRADF3,
!  REQUIRED        VRADF4, VRADF5, VRADFG, VRFFTB, VRFTB1, VRADB2,
!                  VRADB3, VRADB4, VRADB5, VRADBG, PIMACH
!
!  SPECIAL         NONE
!  CONDITIONS
!
!  COMMON          NONE
!  BLOCKS
!
!  I/O             NONE
!
!  PRECISION       SINGLE
!
!  SPECIALIST      ROLAND SWEET
!
!  LANGUAGE        FORTRAN
!
!  HISTORY         WRITTEN BY LINDA LINDGREN AND ROLAND SWEET AT THE
!                  NATIONAL BUREAU OF STANDARDS (BOULDER).
!
!  ALGORITHM       A REAL VARIANT OF THE STOCKHAM AUTOSORT VERSION
!                  OF THE COOLEY-TUKEY FAST FOURIER TRANSFORM.
!
!  PORTABILITY     AMERICAN NATIONAL STANDARDS INSTITUTE FORTRAN 77.
!                  THE ONLY MACHINE DEPENDENT CONSTANT IS LOCATED IN
!                  THE FUNCTION PIMACH.
!
!  REQUIRED        COS,SIN
!  RESIDENT
!  ROUTINES
!
!
!***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
!            Computations, (G. Rodrigue, ed.), Academic Press, 1982,
!            pp. 51-83.
!***ROUTINES CALLED  VRFTI1
!***END PROLOGUE  VRFFTI
!
!  VRFFTPK, VERSION 1, AUGUST 1985
!
  DIMENSION       wsave(n+15)
!***FIRST EXECUTABLE STATEMENT  VRFFTI
  IF (n == 1) RETURN
  CALL vrfti1 (n,wsave(1),wsave(n+1))
  RETURN
END SUBROUTINE vrffti

SUBROUTINE vrftb1 (m,n,c,ch,mdimc,wa,fac)
!
!  VRFFTPK, VERSION 1, AUGUST 1985
!
  DIMENSION       ch(mdimc,n), c(mdimc,n), wa(n) ,fac(15)
  nf = fac(2)
  na = 0
  l1 = 1
  iw = 1
  DO k1=1,nf
    ip = fac(k1+2)
    l2 = ip*l1
    ido = n/l2
    idl1 = ido*l1
    IF (ip /= 4) GO TO 103
    ix2 = iw+ido
    ix3 = ix2+ido
    IF (na /= 0) GO TO 101
    CALL vradb4 (m,ido,l1,c,ch,mdimc,wa(iw),wa(ix2),wa(ix3))
    GO TO 102
    101    CALL vradb4 (m,ido,l1,ch,c,mdimc,wa(iw),wa(ix2),wa(ix3))
    102    na = 1-na
    GO TO 115
    103    IF (ip /= 2) GO TO 106
    IF (na /= 0) GO TO 104
    CALL vradb2 (m,ido,l1,c,ch,mdimc,wa(iw))
    GO TO 105
    104    CALL vradb2 (m,ido,l1,ch,c,mdimc,wa(iw))
    105    na = 1-na
    GO TO 115
    106    IF (ip /= 3) GO TO 109
    ix2 = iw+ido
    IF (na /= 0) GO TO 107
    CALL vradb3 (m,ido,l1,c,ch,mdimc,wa(iw),wa(ix2))
    GO TO 108
    107    CALL vradb3 (m,ido,l1,ch,c,mdimc,wa(iw),wa(ix2))
    108    na = 1-na
    GO TO 115
    109    IF (ip /= 5) GO TO 112
    ix2 = iw+ido
    ix3 = ix2+ido
    ix4 = ix3+ido
    IF (na /= 0) GO TO 110
    CALL vradb5 (m,ido,l1,c,ch,mdimc,wa(iw),wa(ix2),wa(ix3),wa(ix4))
    GO TO 111
    110 CALL vradb5 (m,ido,l1,ch,c,mdimc,wa(iw),wa(ix2),wa(ix3),wa(ix4))
    111    na = 1-na
    GO TO 115
    112    IF (na /= 0) GO TO 113
    CALL vradbg (m,ido,ip,l1,idl1,c,c,c,ch,ch,mdimc,wa(iw))
    GO TO 114
    113    CALL vradbg (m,ido,ip,l1,idl1,ch,ch,ch,c,c,mdimc,wa(iw))
    114    IF (ido == 1) na = 1-na
    115    l1 = l2
    iw = iw+(ip-1)*ido
  END DO
  scale=SQRT(1./n)
  IF (na == 0) GO TO 118
  DO j=1,n
    DO i=1,m
      c(i,j) = scale*ch(i,j)
    END DO
  END DO
  RETURN
  118 DO j=1,n
    DO i=1,m
      c(i,j)=scale*c(i,j)
    END DO
  END DO
  RETURN
END SUBROUTINE vrftb1

SUBROUTINE vrftf1 (m,n,c,ch,mdimc,wa,fac)
!
!  VRFFTPK, VERSION 1, AUGUST 1985
!
  DIMENSION       ch(mdimc,n) ,c(mdimc,n)  ,wa(n)   ,fac(15)
  nf = fac(2)
  na = 1
  l2 = n
  iw = n
  DO k1=1,nf
    kh = nf-k1
    ip = fac(kh+3)
    l1 = l2/ip
    ido = n/l2
    idl1 = ido*l1
    iw = iw-(ip-1)*ido
    na = 1-na
    IF (ip /= 4) GO TO 102
    ix2 = iw+ido
    ix3 = ix2+ido
    IF (na /= 0) GO TO 101
    CALL vradf4 (m,ido,l1,c,ch,mdimc,wa(iw),wa(ix2),wa(ix3))
    GO TO 110
    101    CALL vradf4 (m,ido,l1,ch,c,mdimc,wa(iw),wa(ix2),wa(ix3))
    GO TO 110
    102    IF (ip /= 2) GO TO 104
    IF (na /= 0) GO TO 103
    CALL vradf2 (m,ido,l1,c,ch,mdimc,wa(iw))
    GO TO 110
    103    CALL vradf2 (m,ido,l1,ch,c,mdimc,wa(iw))
    GO TO 110
    104    IF (ip /= 3) GO TO 106
    ix2 = iw+ido
    IF (na /= 0) GO TO 105
    CALL vradf3 (m,ido,l1,c,ch,mdimc,wa(iw),wa(ix2))
    GO TO 110
    105    CALL vradf3 (m,ido,l1,ch,c,mdimc,wa(iw),wa(ix2))
    GO TO 110
    106    IF (ip /= 5) GO TO 108
    ix2 = iw+ido
    ix3 = ix2+ido
    ix4 = ix3+ido
    IF (na /= 0) GO TO 107
    CALL vradf5(m,ido,l1,c,ch,mdimc,wa(iw),wa(ix2),wa(ix3),wa(ix4))
    GO TO 110
    107 CALL vradf5 (m,ido,l1,ch,c,mdimc,wa(iw),wa(ix2),wa(ix3),wa(ix4))
    GO TO 110
    108    IF (ido == 1) na = 1-na
    IF (na /= 0) GO TO 109
    CALL vradfg (m,ido,ip,l1,idl1,c,c,c,ch,ch,mdimc,wa(iw))
    na = 1
    GO TO 110
    109    CALL vradfg (m,ido,ip,l1,idl1,ch,ch,ch,c,c,mdimc,wa(iw))
    na = 0
    110    l2 = l1
  END DO
  scale=SQRT(1./n)
  IF (na == 1) GO TO 113
  DO j=1,n
    DO i=1,m
      c(i,j) = scale*ch(i,j)
    END DO
  END DO
  RETURN
  113 DO j=1,n
    DO i=1,m
      c(i,j)=scale*c(i,j)
    END DO
  END DO
  RETURN
END SUBROUTINE vrftf1

SUBROUTINE vrfti1 (n,wa,fac)
!
!  VRFFTPK, VERSION 1, AUGUST 1985
!
  DIMENSION       wa(n)      ,fac(15)    ,ntryh(4)
  DATA ntryh(1),ntryh(2),ntryh(3),ntryh(4)/4,2,3,5/
  nl = n
  nf = 0
  j = 0
  101 j = j+1
  IF (j-4 > 0) THEN
    GO TO   103
  END IF
  ntry = ntryh(j)
  GO TO 104
  103 ntry = ntry+2
  104 nq = nl/ntry
  nr = nl-ntry*nq
  IF (nr == 0) THEN
    GO TO   105
  ELSE
    GO TO   101
  END IF
  105 nf = nf+1
  fac(nf+2) = ntry
  nl = nq
  IF (ntry /= 2) GO TO 107
  IF (nf == 1) GO TO 107
  DO i=2,nf
    ib = nf-i+2
    fac(ib+2) = fac(ib+1)
  END DO
  fac(3) = 2
  107 IF (nl /= 1) GO TO 104
  fac(1) = n
  fac(2) = nf
  tpi = 2.*pimach(1.0)
  argh = tpi/FLOAT(n)
  is = 0
  nfm1 = nf-1
  l1 = 1
  IF (nfm1 == 0) RETURN
  DO k1=1,nfm1
    ip = fac(k1+2)
    ld = 0
    l2 = l1*ip
    ido = n/l2
    ipm = ip-1
    DO j=1,ipm
      ld = ld+l1
      i = is
      argld = FLOAT(ld)*argh
      fi = 0.
      DO ii=3,ido,2
        i = i+2
        fi = fi+1.
        arg = fi*argld
        wa(i-1) = COS(arg)
        wa(i) = SIN(arg)
      END DO
      is = is+ido
    END DO
    l1 = l2
  END DO
  RETURN
END SUBROUTINE vrfti1
