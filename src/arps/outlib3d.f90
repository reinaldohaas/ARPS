!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE A3DMAX0                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE a3dmax0(a,m1,m2,i1,i2,n1,n2,j1,j2,l1,l2,k1,k2,amax,amin)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Find the maximum and minimum of a 3-D array, a, in a specified
!  subdomain.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/01/91.
!
!  MODIFICATION HISTORY:
!
!  5/06/92 (M. Xue)
!  Added full documentation.
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    a        A 3-D array whose max. and min. values are sought here.
!
!    m1,m2    i-index of array, a.
!    i1,i2    i-index defining a subdomain of array, a, where
!             the max. and min. of the array is sought
!
!    n1,n2    j-index of array, a.
!    j1,j2    j-index defining a subdomain of array, a, where
!             the max. and min. of the array is sought
!
!    l1,l2    k-index of array, a.
!    k1,k2    k-index defining a subdomain of array, a, where
!             the max. and min. of the array is sought
!
!
!  OUTPUT:
!
!    amax     The maximum value of an array, a, in a specified
!             subdomain.
!    amin     The minimum value of an array, a, in a specified
!             subdomain.
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: m1,n1,l1,m2,n2,l2
  REAL :: a(m1:m2,n1:n2,l1:l2)

  REAL :: amin, amax
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i1,i2,j1,j2,k1,k2,i,j,k
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  amax=a(i1,j1,k1)

  DO k=k1,k2
    DO j=j1,j2
      DO i=i1,i2
        amax = MAX(amax, a(i,j,k))
      END DO
    END DO
  END DO

  amin=a(i1,j1,k1)

  DO k=k1,k2
    DO j=j1,j2
      DO i=i1,i2
        amin = MIN(amin, a(i,j,k))
      END DO
    END DO
  END DO

  IF (mp_opt > 0) CALL mpmax0(amax,amin)

  RETURN
END SUBROUTINE a3dmax0

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE A3DMAX0LCL                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE a3dmax0lcl(a,m1,m2,i1,i2,n1,n2,j1,j2,l1,l2,k1,k2,amax,amin)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Find the maximum and minimum of a 3-D array, a, in a specified
!  subdomain but only on the local processor when running message
!  passing version.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/01/91.
!
!  MODIFICATION HISTORY:
!
!  5/06/92 (M. Xue)
!  Added full documentation.
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    a        A 3-D array whose max. and min. values are sought here.
!
!    m1,m2    i-index of array, a.
!    i1,i2    i-index defining a subdomain of array, a, where
!             the max. and min. of the array is sought
!
!    n1,n2    j-index of array, a.
!    j1,j2    j-index defining a subdomain of array, a, where
!             the max. and min. of the array is sought
!
!    l1,l2    k-index of array, a.
!    k1,k2    k-index defining a subdomain of array, a, where
!             the max. and min. of the array is sought
!
!
!  OUTPUT:
!
!    amax     The maximum value of an array, a, in a specified
!             subdomain.
!    amin     The minimum value of an array, a, in a specified
!             subdomain.
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: m1,n1,l1,m2,n2,l2
  REAL :: a(m1:m2,n1:n2,l1:l2)

  REAL :: amin, amax
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i1,i2,j1,j2,k1,k2,i,j,k
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  amax=a(i1,j1,k1)

  DO k=k1,k2
    DO j=j1,j2
      DO i=i1,i2
        amax = MAX(amax, a(i,j,k))
      END DO
    END DO
  END DO

  amin=a(i1,j1,k1)

  DO k=k1,k2
    DO j=j1,j2
      DO i=i1,i2
        amin = MIN(amin, a(i,j,k))
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE a3dmax0lcl
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE A3DMAX                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE a3dmax(a,m1,m2,i1,i2,n1,n2,j1,j2,l1,l2,k1,k2,                &
           amax,amin,imax,jmax,kmax, imin,jmin,kmin)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Find the maximum, minimum and the index locations of an array, a,
!  in a specified subdomain.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/01/91.
!
!  MODIFICATION HISTORY:
!
!  5/06/92 (M. Xue)
!  Added full documentation.
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    a        A 3-D array whose max. and min. values are sought here.
!
!    m1,m2    i-index of array, a.
!    i1,i2    i-index defining a subdomain of array, a, where
!             the max. and min. of the array is sought
!
!    n1,n2    j-index of array, a.
!    j1,j2    j-index defining a subdomain of array, a, where
!             the max. and min. of the array is sought
!
!    l1,l2    k-index of array, a.
!    k1,k2    k-index defining a subdomain of array, a, where
!             the max. and min. of the array is sought
!
!
!  OUTPUT:
!
!    amax     The maximum value of an array, a, in a specified
!             subdomain.
!    amin     The minimum value of an array, a, in a specified
!             subdomain.

!    imax,jmax,kmax  The index location of the maximum value (i,j,k)
!    imin,jmin,kmin  The index location of the minimum value (i,j,k)
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: m1,n1,l1,m2,n2,l2
  REAL :: a(m1:m2,n1:n2,l1:l2)

  REAL :: amin, amax
  INTEGER :: imax,jmax,kmax, imin,jmin,kmin
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i1,i2,j1,j2,k1,k2,i,j,k
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  imax = i1
  jmax = j1
  kmax = k1

  amax=a(i1,j1,k1)

  DO k=k1,k2
    DO j=j1,j2
      DO i=i1,i2
        IF( a(i,j,k) > amax ) THEN
          amax = a(i,j,k)
          imax = i
          jmax = j
          kmax = k
        END IF
      END DO
    END DO
  END DO

  imin = i1
  jmin = j1
  kmin = k1

  amin=a(i1,j1,k1)

  DO k=k1,k2
    DO j=j1,j2
      DO i=i1,i2
        IF( a(i,j,k) < amin ) THEN
          amin = a(i,j,k)
          imin = i
          jmin = j
          kmin = k
        END IF
      END DO
    END DO
  END DO

  IF (mp_opt > 0) CALL mpmax(amax,amin,m2,n2,l2,imax,jmax,kmax,imin,jmin,kmin)

  RETURN
END SUBROUTINE a3dmax
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRIGAR                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE wrigar(a,nx0,nx1,ny0,ny1,nz0,nz1,i0,i1,j0,j1,k0,k1,          &
           title,tkoff,mode)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    Orchestrate the formatted printing of select 2-d slices from 3-D
!    arrays. A subdomain (i0:i1,j0:j1,k0:k1) is selected for printing.
!    A maximum width of 32 numbers can be printed in a horizontal
!    direction. A portion of the array table may be truncated. In the
!    case when the  arrays are large, the data are printed for every
!    other data point.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  12/10/89.
!
!  MODIFICATION HISTORY:
!
!  5/20/92 (M. Xue)
!  Added full documentation.
!
!  9/20/93 (A. Sathye)
!  Fixed problem with late nx,ny.
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    a        A 3-D array whose content will be printed.
!
!    nx0,nx1  i-index of array a.
!    ny0,ny1  j-index of array a.
!    nz0,nz1  k-index of array a.
!
!    i0,i1    i-index defining a subdomain of array, a, where
!             the max. and min. of the array is sought
!    j0,j1    j-index defining a subdomain of array, a, where
!             the max. and min. of the array is sought
!    k0,k1    k-index defining a subdomain of array, a, where
!             the max. and min. of the array is sought
!
!    title    The name of the field (character string).
!    tkoff    A value to be subtracted from array, a (for printing).
!
!    Mode     A print control option for the selection of the slice
!             orientation.
!             = 0, all slices will be printed.
!             = 1, x-y slices will be printed.
!             = 2, x-z slices will be printed.
!             = 3, y-z slices will be printed.
!
!  OUTPUT:
!
!    None.
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx0,nx1,ny0,ny1,nz0,nz1
  INTEGER :: i0,i1,j0,j1,k0,k1

  REAL :: a(nx0:nx1,ny0:ny1,nz0:nz1)
  CHARACTER (LEN=*) :: title
  REAL :: tkoff
  INTEGER :: mode
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,modd,nap,nexp
  REAL :: armin,armax,scale,factor

  INTEGER :: nplan, mstep
  PARAMETER(nplan=1000,mstep=10)

  INTEGER :: iap(nplan)

  INTEGER :: is,js,ks,ndigit,nunit,leng,nxmi,nymi,nzmi

  PARAMETER(ndigit=2 ,nunit=6,leng=132)
  PARAMETER(nxmi=33,nymi=5,nzmi=20)

  REAL :: eps
  PARAMETER(eps=1.0E-30)

  INTEGER :: ixp(9),iyp(9),izp(9)
  SAVE ixp,iyp,izp

  DATA ixp /nxmi,-1,-1,-1,-1,-1,-1,-1,-1/
  DATA iyp /-1,-1,-1,-1,-1,-1,-1,-1,-1/
  DATA izp /-1,-1,-1,-1,-1,-1,-1,-1,-1/
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  DO i = 1, nplan

    iap(i) = i

  END DO
!
!-----------------------------------------------------------------------
!
!  Calculate the max and min for scaling:
!
!-----------------------------------------------------------------------
!
!
  armin=a(i0,j0,k0)
  armax=armin

  DO k=k0,k1
    DO j=j0,j1
      DO i=i0,i1
        armin=MIN(a(i,j,k),armin)
        armax=MAX(a(i,j,k),armax)
      END DO
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Determine the scaling factor:
!
!-----------------------------------------------------------------------
!
!
  WRITE(nunit,1001) title,armin,armax
  IF(armax == armin .OR. mode == 4) RETURN

  armax=armax-tkoff
  armin=armin-tkoff

  scale=1.0E-30+MAX(ABS(armax),ABS(armin))
  nexp=INT(ALOG10(scale))-ndigit

  IF(nint(scale*10.**(-nexp)) == 10**(-nexp+1)) nexp=nexp+1

  IF(ALOG10(scale) <= 0.)  nexp=nexp-1

  factor=10.**(-nexp)

  IF(k0 == k1) THEN

    modd=1

  ELSE IF(j0 == j1) THEN

    modd=2

  ELSE IF(i0 == i1) THEN

    modd=3

  ELSE

    modd=mode

  END IF

  is = MAX(1, nint( (i1-i0+1)/32.0 ))
  js = MAX(1, nint( (j1-j0+1)/32.0 ))
  ks = 1

  IF(modd == 0) THEN

    CALL outarr(a,nx0,nx1,i0,i1,is,ixp,ny0,ny1,j0,j1,js,iyp,            &
                nz0,nz1,k0,k1,ks,izp,nplan,leng,title,nunit,            &
                scale,nexp,factor,ndigit,tkoff,modd)

  ELSE

    IF(modd == 1) THEN

      nap=k1-k0+1

    ELSE IF(modd == 2) THEN

      nap=j1-j0+1

    ELSE

      nap=i1-i0+1

    END IF

    CALL outarr(a,nx0,nx1,i0,i1,is,iap,ny0,ny1,j0,j1,js,                &
                iap,                                                    &
                nz0,nz1,k0,k1,ks,iap,nap,leng,title,nunit,              &
                scale,nexp,factor,ndigit,tkoff,modd)
  END IF

  1001  FORMAT(/t2,a8,t30,'Min=',e14.7,' Max=',e14.7)

  RETURN
END SUBROUTINE wrigar

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE OUTARR                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE outarr (a,l0,l1,i0,i1,is,ixp,m0,m1,j0,j1,js,iyp,             &
           n0,n1,k0,k1,ks,isp,nplan,leng,title,nchan,scale,nexp,        &
           factor,ndigit,tkoff,mode)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Print a formatted table of a 2-d array. This subroutine is called
!  by WRIGRA.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  12/10/89.
!
!  MODIFICATION HISTORY:
!
!  5/06/92 (M. Xue)
!  Added full documentation.
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: l0,l1,m0,m1,n0,n1
  REAL :: a(l0:l1,m0:m1,n0:n1)
  REAL :: tkoff, factor, scale
  INTEGER :: mode, ndigit

  CHARACTER (LEN=*) :: title
  CHARACTER (LEN=9) :: cmod(3)
  CHARACTER (LEN=2) :: numero(39)
  CHARACTER (LEN=20) :: forma1
  CHARACTER (LEN=10) :: forma2
  CHARACTER (LEN=13) :: forma3

  INTEGER :: nplan
  INTEGER :: ixp(nplan),iyp(nplan),isp(nplan)

  INTEGER :: leng,lengl,i0,i1,j0,j1,k0,k1,is,js,ks
  INTEGER :: nchan,nexp,lfield,line,ndec,ndig
  INTEGER :: i11,i12,li,imax,kp,j11,j12,mi,jmax,jp,ip
  INTEGER :: i,j,k

  SAVE cmod,forma1,numero,forma2,forma3
  DATA cmod/'plane is=','plane iy=','plane ix='/
  DATA numero/'01','02','03','04','05','06','07','08','09'              &
      ,'10','11','12','13','14','15','16','17','18','19'                &
      ,'20','21','22','23','24','25','26','27','28','29'                &
      ,'30','31','32','33','34','35','36','37','38','39'/
  DATA forma1/'(1x,i2,''*'',T6,31I04)'/
  DATA forma2/'(t6,31i04)'/
  DATA forma3/'(i4,04e30.23)'/

!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  lengl=leng-6
  lfield=ndigit+2
  line=lengl/lfield
  forma1(18:19)=numero(ndigit+2)
  forma1(15:16)=numero(line)
  forma2(8:9)=numero(ndigit+2)
  forma2(5:6)=numero(line)
  ndec=MAX(1,ndigit+2-7)
  ndig=ndec+7
  forma3(8:9)=numero(ndig)
  forma3(11:12)=numero(ndec)
  forma3(5:6)=numero(line)

  IF(mode == 0 .OR. mode == 1) THEN

    i11=i0+MIN(i1-i0,(line-1)*is)
    imax=MIN(lengl,(lfield*((i1-i0)/is+1)))
    li=i0+(line-1)*is
    i12=MIN(i1,li)

    DO kp=1,nplan

!      IF(isp(kp).ge.n0) THEN
      k=k0-1+isp(kp)
      WRITE(nchan,10) title,tkoff,cmod(1),k,nexp
      WRITE(nchan,forma2) (i,i=i0,i11,is)
      WRITE(nchan,4) ('*',i=2,imax)

      DO j=j1,j0,-js

        IF(ndigit < 15) THEN
          WRITE(nchan,forma1) j,( nint((a(i,j,k)-tkoff)*factor)         &
                ,i=i0,i12,is)
        ELSE
          WRITE(nchan,forma3) j,(a(i,j,k),i=i0,i12,is)
        END IF

      END DO

      WRITE(nchan,4) ('*',i=2,imax)
!      ENDIF

    END DO

  END IF

  IF(mode == 0 .OR. mode == 2) THEN

    i11=i0+MIN(i1-i0,(line-1)*is)
    imax=MIN(lengl,(lfield*((i1-i0)/is+1)))
    li=i0+(line-1)*is
    i12=MIN(i1,li)

    DO jp=1,nplan

!      IF(iyp(jp).ge.j0) THEN
      j=j0-1+iyp(jp)
      WRITE(nchan,10) title,tkoff,cmod(2),j,nexp
      WRITE(nchan,forma2) (i,i=i0,i11,is)
      WRITE(nchan,4) ('*',i=2,imax)

      DO k=k1,k0,-ks

        IF(ndigit < 15) THEN
          WRITE(nchan,forma1) k,(nint((a(i,j,k)-tkoff)*factor)          &
                ,i=i0,i12,is)
        ELSE
          WRITE(nchan,forma3) k,(a(i,j,k),i=i0,i12,is)
        END IF

      END DO
      WRITE(nchan,4) ('*',i=2,imax)

!      ENDIF

    END DO

  END IF

  IF(mode == 0 .OR. mode == 3) THEN

    j11=j0+MIN(j1-j0,(line-1)*js)
    jmax=MIN(lengl,(lfield*((j1-j0)/js+1)))
    mi=j0+(line-1)*js
    j12=MIN(j1,mi)

    DO ip=1,nplan

!      IF(ixp(ip).ge.i0) THEN
      i=i0-1+ixp(ip)
      WRITE(nchan,10) title,tkoff,cmod(3),i,nexp
      WRITE(nchan,forma2) (j,j=j11,j0,-js)
      WRITE(nchan,4) ('*',j=2,jmax)

      DO k=k1,k0,-ks

        IF(ndigit < 15) THEN
          WRITE(nchan,forma1) k,(nint((a(i,j,k)-tkoff)*factor)          &
                ,j=j12,j0,-js)
        ELSE
          WRITE(nchan,forma3) k,(a(i,j,k),j=j12,j0,-js)
        END IF

      END DO

      WRITE(nchan,4) ('*',j=2,jmax)
!      ENDIF

    END DO

  END IF

  RETURN

!  2     FORMAT(t30,'Units of 10**',i3)
!  3     FORMAT(t6,31I4)
  4     FORMAT(t7,125A1)
!  5     FORMAT(1X,i2,'*',t6,31I4)
  10    FORMAT(//t7,a,'(add:',e10.3,')',3X,a,i4,5X,'Units of 10**'      &
           ,i3,/)

END SUBROUTINE outarr
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GTBASFN                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE gtbasfn(fnkey,dirname,ldirnam,hdmpfmt,mgrid,nestgrd,         &
           basdmpfn,lbasdmpf )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Return a unique name for the grid and base state array dump file.
!  The naming convention of the history data dump is:
!
!      fnkey.bingrdbas.aa unformatted binary data set.
!      fnkey.ascgrdbas.aa formatted ASCII data set.
!      fnkey.hdfgrdbas.aa HDF data set
!      fnkey.pakgrdbas.aa packed binary data set
!
!  where fnkey is a string for name construction and aa is a two
!  digit number appended to the data set name in case a data file
!  called fnkey.hdfgrdbas.(aa-1) already exists on the disk.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  11/20/91
!
!  MODIFICATION HISTORY:
!
!  5/05/92 (M. Xue)
!  Added full documentation.
!
!  7/30/92 (M. Xue)
!  Added option for packed data set.
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    fnkey    A character string used to name the files
!    dirname  A file directory name to be attached to the file name.
!    ldirnam  Length of the directory string name
!    hdmpfmt  Parameter specifying the format of output data set.
!    mgrid    The grid number
!    nestgrd  Flag for nested grid run.
!
!-----------------------------------------------------------------------
!
!  OUTPUT:
!
!    basdmpfn The name of the data dump file.
!    lbasdmpf The length of character string basdmpfn.
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  CHARACTER (LEN=* ) :: fnkey   ! A character string used to name the
                                ! files
  CHARACTER (LEN=* ) :: dirname ! A string giving the file directory
                                ! name
  INTEGER :: ldirnam         ! The length of the directory name
  INTEGER :: hdmpfmt         ! Parameter specifying the format of
                             ! output data set.
  INTEGER :: mgrid           ! The grid number
  INTEGER :: nestgrd         ! Flag for nested grid run.
  CHARACTER (LEN=* ) :: basdmpfn ! Name of grid & base state data file.
  INTEGER :: lbasdmpf        ! Length of character string basdmpfn.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: lfnkey
  CHARACTER (LEN=256) :: temchar
  INTEGER :: istatus
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  lfnkey = LEN( fnkey)

  IF( hdmpfmt == 1.OR.hdmpfmt == 0) THEN       ! Unformatted binary dump
    basdmpfn = fnkey(1:lfnkey)//'.bingrdbas'
  ELSE IF( hdmpfmt == 2 ) THEN   ! Formatted ASCII dump
    basdmpfn = fnkey(1:lfnkey)//'.ascgrdbas'
  ELSE IF( hdmpfmt == 3 ) THEN   ! HDF data dump
    basdmpfn = fnkey(1:lfnkey)//'.hdfgrdbas'
  ELSE IF( hdmpfmt == 4 ) THEN   ! Packed binary dump
    basdmpfn = fnkey(1:lfnkey)//'.pakgrdbas'
  ELSE IF( hdmpfmt == 5 ) THEN   ! Savi3D data dump

!-----------------------------------------------------------------------
!    For Savi3D data dump, the grid and base state information is
!    always written together with the other fields.
!-----------------------------------------------------------------------

  ELSE IF( hdmpfmt == 6 ) THEN   ! Binary with skipping
    basdmpfn = fnkey(1:lfnkey)//'.bn2grdbas'
  ELSE IF( hdmpfmt == 7 ) THEN   ! NetCDF format
    basdmpfn = fnkey(1:lfnkey)//'.netgrdbas'
  ELSE IF( hdmpfmt == 8 ) THEN   ! Packed NetCDF format
    basdmpfn = fnkey(1:lfnkey)//'.ncgrdbas'
  ELSE IF( hdmpfmt == 9 ) THEN   ! GrADS data dump
!-----------------------------------------------------------------------
!    For GrADS data dump, the grid and base state information is
!    always written together with the other fields.
!-----------------------------------------------------------------------
  ELSE IF( hdmpfmt == 10 ) THEN  ! GRIB format
    basdmpfn = fnkey(1:lfnkey)//'.grbgrdbas'
  END IF

  lbasdmpf = 10 + LEN(fnkey)

  IF(nestgrd == 1) THEN

    WRITE(basdmpfn((lbasdmpf+1):(lbasdmpf+4)),'(a,i2.2)')'.g',mgrid
    lbasdmpf = lbasdmpf + 4

  END IF

  IF (mp_opt > 0 .AND. joindmp(FINDX_H) == 0) THEN
    temchar = basdmpfn
    CALL gtsplitfn(temchar,1,1,loc_x,loc_y,1,1,0,0,0,2,basdmpfn,istatus)
    lbasdmpf  = LEN_TRIM(basdmpfn)
  END IF

  IF( dirname /= ' ' ) THEN

    temchar = basdmpfn
    basdmpfn = dirname(1:ldirnam)//'/'//temchar
    lbasdmpf  = lbasdmpf + ldirnam + 1

  END IF

  CALL fnversn(basdmpfn, lbasdmpf)

  RETURN
END SUBROUTINE gtbasfn

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GTDMPFN                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE gtdmpfn(fnkey,dirname,ldirnam,curtim,hdmpfmt,                &
                   mgrid,nestgrd, hdmpfn, ldmpf )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Return a unique name for the history data dump at time 'curtim'.
!  The naming convention of the history data dump is:
!
!      fnkey.hdfnnnnnn.aa for hdf data set
!      fnkey.binnnnnnn.aa for unformatted binary data set.
!
!  where fnkey is a string for name construction, nnnnnn indicates
!  the time of the data set in hour/minute/second format, and aa is
!  a two digit number appended to the data set name in case a data
!  file called fnkey.hdfnnnnnn.(aa-1) already exists on the disk.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  11/20/91
!
!  MODIFICATION HISTORY:
!
!  5/05/92 (M. Xue)
!  Added full documentation.
!
!  7/30/92 (M. Xue)
!  Added option for packed data set.
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    fnkey    A character string for the name of this run.
!    dirname  A file directory name to be attached to the file name.
!    ldirnam  Length of the directory string name
!    curtim   The model time in seconds.
!    hdmpfmt  Parameter specifying the format of output data set.
!    mgrid    The grid number
!    nestgrd  Flag for nested grid run.
!
!-----------------------------------------------------------------------
!
!  OUTPUT:
!
!    hdmpfn   The name of the history data dump file.
!    ldmpf    The length of character string hdmpfn.
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  CHARACTER (LEN=* ) :: fnkey   ! A character string used to name the
                                ! files
  CHARACTER (LEN=* ) :: dirname ! A string giving the file directory
                                ! name.
  INTEGER            :: ldirnam ! The length of the directory name.
  REAL               :: curtim  ! Current model time.
  INTEGER            :: hdmpfmt ! Parameter specifying the format of
                                ! output data set.
  INTEGER            :: mgrid    ! The grid number.
  INTEGER            :: nestgrd  ! Flag for nested grid run.
  CHARACTER (LEN=* ) :: hdmpfn   ! Name of the history data dump file.
  INTEGER            :: ldmpf    ! Length of character string hdmpfn.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=80 ) :: timsnd
  CHARACTER (LEN=256) :: temchar
  INTEGER :: tmstrln
  INTEGER :: istatus
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  CALL cvttsnd( curtim, timsnd, tmstrln )

  IF( hdmpfmt == 1.OR.hdmpfmt == 0) THEN       ! Unformatted binary dump
    WRITE(hdmpfn,'(a,a)') fnkey, '.bin'//timsnd(1:tmstrln)
  ELSE IF( hdmpfmt == 2 ) THEN   ! Formatted ASCII dump
    WRITE(hdmpfn,'(a,a)') fnkey, '.asc'//timsnd(1:tmstrln)
  ELSE IF( hdmpfmt == 3 ) THEN   ! HDF data dump
    WRITE(hdmpfn,'(a,a)') fnkey, '.hdf'//timsnd(1:tmstrln)
  ELSE IF( hdmpfmt == 4 ) THEN   ! Packed binary dump
    WRITE(hdmpfn,'(a,a)') fnkey, '.pak'//timsnd(1:tmstrln)
  ELSE IF( hdmpfmt == 5 ) THEN   ! Data dump for Savi3D
    WRITE(hdmpfn,'(a,a)') fnkey, '.svi'
  ELSE IF( hdmpfmt == 6 ) THEN   ! Binary with skipping
    WRITE(hdmpfn,'(a,a)') fnkey, '.bn2'//timsnd(1:tmstrln)
  ELSE IF( hdmpfmt == 7 ) THEN   ! NetCDF format
    WRITE(hdmpfn,'(a,a)') fnkey, '.net'//timsnd(1:tmstrln)
  ELSE IF( hdmpfmt == 8 ) THEN   ! NetCDF one file
    WRITE(hdmpfn,'(a,a)') fnkey, '.nc'
  ELSE IF( hdmpfmt == 9 ) THEN   ! Data dump for GrADS
    WRITE(hdmpfn,'(a,a)') fnkey, '.gad'
  ELSE IF( hdmpfmt == 10 ) THEN   ! Data dump for GrADS
    WRITE(hdmpfn,'(a,a)') fnkey, '.grb'//timsnd(1:tmstrln)
  ELSE IF( hdmpfmt == 11 ) THEN   ! Data dump for Vis5D
    WRITE(hdmpfn,'(a,a)') fnkey, '.v5d'//timsnd(1:tmstrln)
  END IF

  IF(  hdmpfmt == 5 .OR. hdmpfmt == 9 ) THEN
    ldmpf = LEN(fnkey) + 4
  ELSE IF( hdmpfmt == 8 ) THEN
    ldmpf = LEN(fnkey) + 3
  ELSE
    ldmpf = LEN(fnkey) + 4 + tmstrln
  END IF

  IF(nestgrd == 1) THEN

    WRITE(hdmpfn((ldmpf+1):(ldmpf+4)), '(a,i2.2)') '.g',mgrid
    ldmpf = ldmpf + 4

  END IF

  IF (mp_opt > 0 .AND. joindmp(FINDX_H) == 0) THEN
    temchar = hdmpfn
    CALL gtsplitfn(temchar,1,1,loc_x,loc_y,1,1,0,0,0,2,hdmpfn,istatus)
  END IF

  IF( dirname /= ' ' ) THEN
    temchar = hdmpfn
    hdmpfn  = dirname(1:ldirnam)//temchar
  END IF

  ldmpf = LEN_TRIM(hdmpfn)

  CALL fnversn(hdmpfn, ldmpf)

  RETURN
END SUBROUTINE gtdmpfn
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GTRSTFN                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE gtrstfn(fnkey,dirname,ldirnam,curtim,                        &
           mgrid,nestgrd, rstoutf, lrstof )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Return a unique name for the restart data dump at time 'curtim'.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  03/12/93
!
!  MODIFICATION HISTORY:
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    fnkey    A character string used to name the files
!    dirname  A file directory name to be attached to the file name.
!    ldirnam  Length of the directory string name
!    curtim   The model time in seconds.
!    mgrid    The grid number
!    nestgrd  Flag for nested grid run.
!
!-----------------------------------------------------------------------
!
!  OUTPUT:
!
!    rstoutf  The name of the restart dump file.
!    lrstof   The length of the file name.
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  CHARACTER (LEN=* ) :: fnkey   ! A character string used to name the
                                ! files
  CHARACTER (LEN=* ) :: dirname ! A string giving the file directory
                                ! name.
  INTEGER            :: ldirnam ! The length of the directory name.
  REAL               :: curtim  ! Current model time.
  INTEGER            :: mgrid   ! The grid number.
  INTEGER            :: nestgrd ! Flag for nested grid run.
  CHARACTER (LEN=* ) :: rstoutf ! Name of the history data dump file.
  INTEGER            :: lrstof  ! Length of character string hdmpfn.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=80 ) :: timsnd
  CHARACTER (LEN=256) :: temchar
  INTEGER :: tmstrln
  INTEGER :: istatus
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'mp.inc'            ! Message passing parameters.
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
!  Construct the file name:
!
!  The format of the restart file names is
!
!    fnkey.rstnnnnnn
!
!  where fnkey are characters used to name the files.
!  nnnnnn is a 6 digit integer number indicating the time of
!  the data in seconds.
!
!-----------------------------------------------------------------------
!

  CALL cvttsnd( curtim, timsnd, tmstrln )

  WRITE(rstoutf,'(a,a)') fnkey, '.rst'//timsnd(1:tmstrln)
  lrstof = LEN(fnkey) + 4 + tmstrln

  IF( nestgrd == 1 ) THEN

!-----------------------------------------------------------------------
!
!  Attach the grid number to the file name
!
!-----------------------------------------------------------------------

    WRITE(rstoutf((lrstof+1):(lrstof+4)), '(a,i2.2)') '.g',mgrid
    lrstof = lrstof + 4

  END IF

  IF( dirname /= ' ' ) THEN

    temchar = rstoutf(1:lrstof)
    rstoutf = dirname(1:ldirnam)//'/'//temchar
    lrstof  = lrstof + ldirnam + 1

  END IF

  IF (mp_opt > 0 .AND. joindmp(FINDX_R) <= 0) THEN
    temchar = rstoutf
    CALL gtsplitfn(temchar,1,1,loc_x,loc_y,1,1,0,0,0,2,rstoutf,istatus)
    lrstof  = LEN_TRIM(rstoutf)
  ELSE
!
!-----------------------------------------------------------------------
!
!  Append a version number to the file name if the named file already
!  exists.
!
!-----------------------------------------------------------------------
!
    CALL fnversn(rstoutf, lrstof)
  END IF

  RETURN
END SUBROUTINE gtrstfn

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GTLOGFN                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE gtlogfn( fnkey, mgrid, nestgrd, logfn, llogfn )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Return a unique name for the input log file. The name convention
!  of the log file is:
!
!      fnkey.log.aa
!
!  where fnkey is a character string for naming the files and aa is
!  a two digit number appended to the data set name in case the file
!  fnkey.log.(aa-1) already exists.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  3/17/1991.
!
!  MODIFICATION HISTORY:
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    fnkey    A character string used to name the files
!    mgrid    The grid number
!    nestgrd  Flag for nested grid run.
!
!  OUTPUT:
!
!    logfn    Log file filename.
!    llogfn   The length of character string logfn.
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  CHARACTER (LEN=*  ) :: fnkey ! A character string used to name the
                               ! files
  INTEGER :: mgrid             ! The grid number
  INTEGER :: nestgrd           ! Flag for nested grid run.
  CHARACTER (LEN=256) :: logfn  ! Log file filename.
  INTEGER :: llogfn             ! The length of character string logfn.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  WRITE(logfn,'(a,a)') fnkey, '.log'
  llogfn = 4 + LEN(fnkey)

  IF( nestgrd == 1 ) THEN
    WRITE(logfn((llogfn+1):(llogfn+4)), '(a,i2.2)') '.g',mgrid
    llogfn = llogfn + 4
  END IF

  CALL fnversn(logfn, llogfn )

  RETURN
END SUBROUTINE gtlogfn

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE FNVERSN                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE fnversn( filename, fnlen )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Append the file version number to a file name if the named
!  file already exists.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  3/17/1991.
!
!  MODIFICATION HISTORY:
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    filename name of the file on input
!    fnlen    file name length on input
!
!  OUTPUT:
!
!    filename name of the file on output
!    fnlen    file name length on output
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  CHARACTER (LEN=* ) :: filename ! name of the file on input/output
  INTEGER :: fnlen               ! file name length on input/output
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256)  :: temchar
  LOGICAL :: iexist, iexist1, iexist2
  INTEGER :: nnn, fnlen_old, fnlen_tem
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  temchar = filename

  fnlen_old = fnlen
  fnlen_tem = fnlen

  nnn = 0

  200 CONTINUE

  INQUIRE(FILE=temchar(1:fnlen_tem),        EXIST=iexist)
  INQUIRE(FILE=temchar(1:fnlen_tem)//'.Z' , EXIST=iexist1)
  INQUIRE(FILE=temchar(1:fnlen_tem)//'.gz', EXIST=iexist2)

  IF( iexist .OR. iexist1 .OR. iexist2 ) THEN

    nnn = nnn+1

    IF( nnn > 99) THEN
      WRITE(6,'(/1x,a,/1x,a/)')                                         &
          'An alternative name could not be found for ',                &
          temchar(1:fnlen_old),' Job stopped in FNVERSN.'
      CALL arpsstop('File name not found.',1)
    END IF

    WRITE(temchar((fnlen_old+1):(fnlen_old+3)),'(a,i2.2)')'.',nnn
    fnlen_tem = fnlen_old + 3
    GO TO 200

  END IF

  fnlen = fnlen_tem
  filename(1:fnlen) = temchar(1:fnlen_tem)

  RETURN
END SUBROUTINE fnversn
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE STRLNTH                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE strlnth( string, length )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Return the length of the non-blank part of a character string.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  11/20/91
!
!  MODIFICATION HISTORY:
!
!  5/05/92 (M. Xue)
!  Added full documentation.
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    string   A character string
!    length   The declared length of the character string 'string'.
!
!  OUTPUT:
!
!    length   The length of the non-blank part of the string.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN) :: string ! A character string for the name of
                                         ! this run.
  INTEGER,       INTENT(INOUT) :: length ! The length of the non-blank part
                                         ! of a string.

  INTEGER :: i
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  IF (length == 0) length = LEN(string)

  DO i = length,1,-1

    IF(string(i:i) /= ' ') EXIT

  END DO

  length = MAX(1,i)

  RETURN
END SUBROUTINE strlnth


!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE STRMIN                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE strmin( string, length )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Minimize a string length by removing consecutive blank spaces.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  1/15/93
!
!  MODIFICATION HISTORY:
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    string   A character string
!    length   The declared length of the character string 'string'.
!
!  OUTPUT:
!
!    length   The length of string with consecutive blank spaces
!             removed.
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  CHARACTER (LEN=*     ) :: string ! A character string for the name of
                                   ! this run.
  INTEGER :: length                ! The length of the non-blank part
                                   ! of a string.

  CHARACTER (LEN=1  ) ::  str_1
  CHARACTER (LEN=256) :: str

  INTEGER :: i,len_old
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF( length > 256) THEN
    PRINT*,'Work string defined in STRMIN was too small.'
    PRINT*,'The output from this subroutine may not be correct.'
    length=256
  END IF

  len_old = length
  length = 1

  str = string
  DO i = 2,len_old

    str_1 = str(i-1:i-1)
    IF(.NOT.(str(i:i) == ' '.AND.                                       &
          (str_1 == ' '.OR.str_1 == '('.OR.str_1 == '='))) THEN
      length=length+1
      string(length:length)=str(i:i)
    END IF

  END DO

!  200   CONTINUE

  RETURN
END SUBROUTINE strmin

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE COMLNTH                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE comlnth( string, length )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Return the length of the non-blank part of a character string.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  11/20/91
!
!  MODIFICATION HISTORY:
!
!  5/05/92 (M. Xue)
!  Added full documentation.
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    string   A character string
!    length   The declared length of the character string 'string'.
!
!  OUTPUT:
!
!    length   The length of the non-blank part of the string.
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  CHARACTER (LEN=*     ) :: string ! A character string for the name of
                                   ! this run.
  INTEGER :: length            ! The length of the non-blank part
                               ! of a string.
  INTEGER :: i
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  DO i = 1, 80

    IF(string(i:i) == ' ') THEN
      IF (string(i+1:i+1) == ' ') EXIT
    END IF

  END DO

  length = MAX(1,i-1)

  RETURN
END SUBROUTINE comlnth
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CVTTIM                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE cvttim(timsnd, timhms)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Convert time, given in seconds, into a 6 character string
!  containing time in the hour/minute/second format.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  11/02/92
!
!  MODIFICATION HISTORY:
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    timsnd   Time in seconds
!
!  OUTPUT:
!
!    timhms    string contain time in hour:minute:second format
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!

  IMPLICIT NONE

  REAL :: timsnd               ! Time in seconds
  CHARACTER (LEN=6      ) :: timhms  ! string contain time in
                                     ! hour:minute:second format
  INTEGER :: h,m,s
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  h = INT(timsnd/3600.0 )
  m = INT((timsnd-h*3600.0)/60.0)
  s = nint(timsnd-h*3600.0-m*60.0)

  IF( s == 60) THEN
    m = m+1
    s = 0
  END IF

  WRITE(timhms,'(3i2.2)') h,m,s

  RETURN
END SUBROUTINE cvttim
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CVTTSND                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE cvttsnd(time, timsnd, tmstrln)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Convert time given in second into a character string
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  07/17/2000
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    time     Time in seconds
!
!  OUTPUT:
!
!    timsnd   Time string in seconds
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!

  IMPLICIT NONE
  REAL,              INTENT(IN)  :: time       ! Time in seconds
  CHARACTER (LEN=*), INTENT(OUT) :: timsnd     ! Time string in seconds
  INTEGER,           INTENT(OUT) :: tmstrln    ! Length of time string

  INTEGER :: itime
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  itime = nint(time)

  IF ( itime < 999999 ) THEN
    tmstrln=6
    IF( LEN(timsnd) < tmstrln) GO TO 10
    WRITE(timsnd,'(i6.6)') itime
  ELSE IF ( itime < 9999999 ) THEN
    tmstrln=7
    IF( LEN(timsnd) < tmstrln) GO TO 10
    WRITE(timsnd,'(i7.7)') itime
  ELSE IF ( itime < 99999999 ) THEN
    tmstrln=8
    IF( LEN(timsnd) < tmstrln) GO TO 10
    WRITE(timsnd,'(i8.8)') itime
  ELSE IF ( itime < 999999999 ) THEN
    tmstrln=9
    IF( LEN(timsnd) < tmstrln) GO TO 10
    WRITE(timsnd,'(i9.9)') itime
  ELSE
    WRITE (6,'(a/a,i4/a,e16.8,a)')                                      &
        'WARNING: The time is too large to fit in 9 characters',        &
        '         time    = ',time, ' seconds.'
    CALL arpsstop('arpsstop called from CVTTSND time too large ',1)
  END IF

  RETURN

  10    CONTINUE

  WRITE(6,'(a/a/a)')                                                    &
      'String timsnd passed into CVTTSND not long enough.',             &
      'Need to be at least ',tmstrln,' charactere long.',               &
      'Job stopped in CVTTSND.'
  CALL arpsstop('arpsstop called from CVTTSND string insufficient       &
       & length',1)

  RETURN
END SUBROUTINE cvttsnd
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GTLFNKEY                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE gtlfnkey( runname, lfnkey )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Find out the number of characters to be used to construct file
!  names.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  03/15/93
!
!  MODIFICATION HISTORY:
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!  12/12/1996 (Yuhe Liu)
!  Removed the restrict of 6 characters to runname.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!

  IMPLICIT NONE

  CHARACTER (LEN=* ) :: runname ! Input
  INTEGER :: lfnkey         ! Output

  INTEGER :: lenstr, firstb, firstc
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  lenstr = LEN( runname )

  firstb = INDEX( runname, ' ')
  firstc = INDEX( runname, ',')

  IF( firstb == 0) firstb = lenstr+1
  IF( firstc == 0) firstc = lenstr+1

  lfnkey = MAX(1, MIN( lenstr, firstb-1, firstc-1 ) )

  RETURN
END SUBROUTINE gtlfnkey
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GETUNIT                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE getunit( nunit )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Find a free FORTRAN I/O unit from a list and return that unit.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  4/21/93
!
!  MODIFICATION HISTORY:
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  OUTPUT:
!
!  nunit   A free fortran I/O unit number.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nunit

  LOGICAL :: used

  INTEGER :: list(30), nfree
  SAVE list, nfree
  DATA list /10,11,12,13,14,15,16,17,18,19,                             &
             20,21,22,23,24,25,26,27,28,29,                             &
             30,31,32,33,34,35,36,37,38,39/
  DATA nfree /30/

  INCLUDE 'mp.inc'          ! Message passing parameters.

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  5 CONTINUE
  IF( nfree < 1) THEN
    IF (myproc == 0) WRITE(6,'(1x,a,a)')  &
        'No more unit number is available from the list. ',             &
        'Job stopped in GETUNIT.'
    CALL arpsstop('arpsstop called from GETINIT out of file numbers',1)
  END IF

  nunit = list(nfree)

  nfree = nfree-1

  INQUIRE( UNIT=nunit, OPENED=used)
  IF( used ) GO TO 5

  IF (myproc == 0) WRITE(6,'(1x,a,i3,a)') 'Fortran I/O unit ',nunit,  &
       ' picked from the free list.'

  RETURN
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                   ENTRY RETUNIT                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
  ENTRY retunit( nunit )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Return a freed unit to the list.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  4/21/93
!
!  MODIFICATION HISTORY:
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  nunit   A freed fortran I/O unit number to be returned.
!
!-----------------------------------------------------------------------
!


  INQUIRE( UNIT=nunit, OPENED=used)
  IF( used ) RETURN

  nfree = nfree + 1
  IF( nfree <= 30 ) list( nfree ) = nunit
  nfree = MIN( nfree, 30)
  IF (myproc == 0 ) WRITE(6,'(1x,a,i3,a)') 'Fortran I/O unit ',nunit,  &
       ' returned to the free list.'

  RETURN
END SUBROUTINE getunit
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SECTHRZ                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE secthrz(nx,ny,nz,s,z,ss1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Interpolate 3-D data to a given horizontal level.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue & Hao Jin
!  12/18/92.
!
!  MODIFICATION HISTORY:
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    s        3-dimensional array of data to contour
!    s1       2-dimensional array of data to contour
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in physical space (m)
!
!    z01      value of x for first i grid point to plot
!
!  OUTPUT:
!    ss1      interpolated 3-D data to a given horizontal level
!
!-----------------------------------------------------------------------
!
!  Parameters of output
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: s(nx,ny,nz)          ! 3-dimensional array of data to contour
  REAL :: z(nx,ny,nz)          ! z coordinate of grid points
                               ! in physical space (m)
  REAL :: ss1(nx,ny)           ! interpolated 3-D data to a
                               ! given horizontal level
!
!-----------------------------------------------------------------------
!
!  Common blocks for plotting control parameters
!
!-----------------------------------------------------------------------
!
  REAL :: z01                      ! the given height of the slice
  COMMON /sliceh/z01

!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
!
!-----------------------------------------------------------------------
!
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
!  Find index for interpolation
!
!-----------------------------------------------------------------------
!
  DO i=1,nx-1
    DO j=1,ny-1
      IF(z01 <= z(i,j,1)) GO TO 11
      IF(z01 >= z(i,j,nz-1)) GO TO 12
      DO k=2,nz-2
        IF(z01 >= z(i,j,k).AND.z01 < z(i,j,k+1)) GO TO 15
      END DO

      11    k=1
      GO TO 15
      12    k=nz-1
      GO TO 15

      15    ss1(i,j)=s(i,j,k)+(s(i,j,k+1)-s(i,j,k))*                    &
                     (z01-z(i,j,k))/(z(i,j,k+1)-z(i,j,k))

!-----------------------------------------------------------------------
!
!  If the data point is below the ground level, set the
!  data value to the missing value.
!
!-----------------------------------------------------------------------

      IF( z01 < z(i,j,2) ) ss1(i,j) = -9999.0

    END DO
  END DO

  RETURN
END SUBROUTINE secthrz
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SECTVRT                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE sectvrt(nx,ny,nz,s,x,y,z,dx,dy,ss1,zs1,n,xp,yp)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Interpolate a 3-D data to 2-d vectical plane.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue & Hao Jin
!  12/18/92.
!
!  MODIFICATION HISTORY:
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!  01/07/2003 (Yunheng Wang)
!  Added message passing code. After the call, all processors will
!  get the same output (ss1,zs1).
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    s        3-dimensional array of data to contour
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    dx, dy   Grid spacing in x and y directions
!
!    n
!    xp
!    yp
!
!  OUTPUT:
!    ss1      interpolated a 3-D data to 2-d vectical plane
!    zs1      interpolated a 3-D data to 2-d vectical plane
!
!
!-----------------------------------------------------------------------
!
!  Parameters of output
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER, INTENT(IN) :: n

  REAL, TARGET, INTENT(IN) :: s(nx,ny,nz)  ! 3-dimensional array of data to contour
  REAL, TARGET, INTENT(IN) :: x(nx,ny,nz)  ! x coordinate of grid points
                                           ! in physical/comp. space (m)
  REAL, TARGET, INTENT(IN) :: y(nx,ny,nz)  ! y coordinate of grid points
                                           ! in physical/comp. space (m)
  REAL, TARGET, INTENT(IN) :: z(nx,ny,nz)  ! z coordinate of grid points
                                           ! in computational space (m)
  REAL, INTENT(IN) :: dx, dy

  REAL, INTENT(OUT) :: ss1(n,nz)
  REAL, INTENT(OUT) :: zs1(n,nz)
  REAL, INTENT(IN)  :: xp(n)
  REAL, INTENT(IN)  :: yp(n)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,k
  INTEGER :: is,js
  REAL :: s1,s2,s3,s4,sgrid,xs1,ys1

  INTEGER :: nxlg, nylg
  REAL, POINTER :: xmp(:,:,:), ymp(:,:,:), zmp(:,:,:), smp(:,:,:)
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'mp.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  nxlg = (nx-3)*nproc_x + 3
  nylg = (ny-3)*nproc_y + 3

  IF(mp_opt > 0) THEN
    ALLOCATE(xmp(nxlg,nylg,nz))
    ALLOCATE(ymp(nxlg,nylg,nz))
    ALLOCATE(zmp(nxlg,nylg,nz))
    ALLOCATE(smp(nxlg,nylg,nz))

    CALL mpimerge3d(x,nx,ny,nz,xmp)
    CALL mpimerge3d(y,nx,ny,nz,ymp)
    CALL mpimerge3d(z,nx,ny,nz,zmp)
    CALL mpimerge3d(s,nx,ny,nz,smp)
  ELSE
    xmp => x
    ymp => y
    zmp => z
    smp => s
  END IF

!
!-----------------------------------------------------------------------
!
!  Intepolate the date to the given point.
!
!-----------------------------------------------------------------------
!
  IF(myproc == 0) THEN

    DO k=1,nz-1
      DO i=1,n
        xs1=xp(i)
        ys1=yp(i)

        is = MAX(1, MIN(nxlg-2, INT( (xs1-xmp(1,2,2))/dx + 1 ) ))
        js = MAX(1, MIN(nylg-2, INT( (ys1-ymp(2,1,2))/dy + 1 ) ))

        s1 = (xs1-xmp(is  ,js  ,k))*(ys1-ymp(is  ,js  ,k))
        s2 =-(xs1-xmp(is+1,js  ,k))*(ys1-ymp(is+1,js  ,k))
        s3 = (xs1-xmp(is+1,js+1,k))*(ys1-ymp(is+1,js+1,k))
        s4 =-(xs1-xmp(is  ,js+1,k))*(ys1-ymp(is  ,js+1,k))
        sgrid = (xmp(is+1,js,k)-xmp(is,js,k))*(ymp(is,js+1,k)-ymp(is,js,k))

        zs1(i,k) =(zmp(is  ,js  ,k)*s3+zmp(is+1,js  ,k)*s4              &
                  +zmp(is+1,js+1,k)*s1+zmp(is  ,js+1,k)*s2)/sgrid
        ss1(i,k) =(smp(is  ,js  ,k)*s3+smp(is+1,js  ,k)*s4              &
                  +smp(is+1,js+1,k)*s1+smp(is  ,js+1,k)*s2)/sgrid

      END DO
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Make sure all PEs have the same data and deallocate working arrays
!
!-----------------------------------------------------------------------

  IF (mp_opt > 0) THEN
    CALL mpupdater(zs1,n*nz)
    CALL mpupdater(ss1,n*nz)

    DEALLOCATE(xmp,ymp,zmp,smp)
  END IF

  RETURN
END SUBROUTINE sectvrt
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE REFLEC                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE reflec(nx,ny,nz, rhobar, qr, qs, qh, reflc )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Compute the radar reflectivity factor following Kessler (1969).
!  Here, arg=Z (mm**6/m**3), and dBz = 10log10 (arg).
!
!-----------------------------------------------------------------------
!
!  AUTHOR: K. Droegemeier and M.Xue
!  4/19/93
!
!  MODIFICATION HISTORY:
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!  12/6/95 (J. Zong and M. Xue)
!  Added qs and qh to the argument list of this subroutine to
!  facilitate inclusion of the contributions of qs and qh to reflec-
!  tivity. A relation between radar reflectivity factor and snow
!  content is adopted from Rogers and Yau (1989) and extended to
!  represent the effects of snow and graupel/hail on the
!  reflectivity. globcst.inc is included to pass the value of ice.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    rhobar   Base state density (kg/m**3)
!    qr       Rainwater mixing ratio (kg/kg)
!    qs       Snow mixing ratio (kg/kg)
!    qh       Hail mixing ratio (kg/kg)
!
!  OUTPUT:
!
!    reflc    Radar reflectivity factor.
!
!-----------------------------------------------------------------------
!

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
!  Include files:
!
!-----------------------------------------------------------------------
!
!
  INCLUDE 'globcst.inc'
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz

  REAL :: rhobar(nx,ny,nz)     ! Base state air density (kg/m**3)
  REAL :: qr    (nx,ny,nz)     ! Rain water mixing ratio (kg/kg)
  REAL :: qs    (nx,ny,nz)     ! Snow mixing ratio (kg/kg)
  REAL :: qh    (nx,ny,nz)     ! Hail mixing ratio (kg/kg)

  REAL :: reflc (nx,ny,nz)     ! Radar reflectivity (dBZ)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  REAL    :: arg
  REAL    :: coef, svnfrth
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
!  Compute the radar reflectivity factor following Kessler (1969).
!  Here, arg=Z (mm**6/m**3), and dBz = 10log10 (arg).
!
!-----------------------------------------------------------------------
!
!  IF ( mphyopt >= 5 .AND. mphyopt <= 7 ) THEN
!    svnfrth = 1.24
!    coef    = 20740.0
!  ELSE
    svnfrth = 7./4.
    coef    = 17300.0
!  END IF

  DO k = 1,nz-1
    DO j = 1,ny-1
      DO i = 1,nx-1

        arg = coef*( rhobar(i,j,k)*1000.0 * MAX(0.0,qr(i,j,k)) )**svnfrth

        IF (ice == 1) THEN

          arg = arg + 38000.0*( rhobar(i,j,k)*1000.0                    &
                              * MAX(0.0,qs(i,j,k)+qh(i,j,k)) )**2.2

        END IF

        reflc(i,j,k) = 10.0*ALOG10( MAX(arg,1.0) )

      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE reflec

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE REFLEC_WR                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE reflec_wr(nx,ny,nz, rhobar, qr, reflc )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Compute the radar reflectivity factor following Kessler (1969).
!  Here, arg=Z (mm**6/m**3), and dBz = 10log10 (arg).
!
!  This version is identical to subroutine REFLEC except for the removal
!  of ice arrays and calculations
!
!-----------------------------------------------------------------------
!
!  AUTHOR: K. Droegemeier and M.Xue
!  4/19/93
!
!  MODIFICATION HISTORY:
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!  12/6/95 (J. Zong and M. Xue)
!  Added qs and qh to the argument list of this subroutine to
!  facilitate inclusion of the contributions of qs and qh to reflec-
!  tivity. A relation between radar reflectivity factor and snow
!  content is adopted from Rogers and Yau (1989) and extended to
!  represent the effects of snow and graupel/hail on the
!  reflectivity. globcst.inc is included to pass the value of ice.
!
!  8/29/07 (D. Dawson)
!  Split off from subroutine REFLEC to only include rainwater fields
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    rhobar   Base state density (kg/m**3)
!    qr       Rainwater mixing ratio (kg/kg)
!
!  OUTPUT:
!
!    reflc    Radar reflectivity factor.
!
!-----------------------------------------------------------------------
!

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
!  Include files:
!
!-----------------------------------------------------------------------
!
!
  INCLUDE 'globcst.inc'
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz

  REAL :: rhobar(nx,ny,nz)     ! Base state air density (kg/m**3)
  REAL :: qr    (nx,ny,nz)     ! Rain water mixing ratio (kg/kg)

  REAL :: reflc (nx,ny,nz)     ! Radar reflectivity (dBZ)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  REAL    :: arg
  REAL    :: coef, svnfrth
  REAL    :: n0r
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
!  Compute the radar reflectivity factor following Kessler (1969).
!  Here, arg=Z (mm**6/m**3), and dBz = 10log10 (arg).
!
!-----------------------------------------------------------------------
!
!  IF ( mphyopt >= 5 .AND. mphyopt <= 7 ) THEN
!    svnfrth = 1.24
!    coef    = 20740.0
!  ELSE
    svnfrth = 7./4.
    coef    = 17300.0
!  END IF

  DO k = 1,nz-1
    DO j = 1,ny-1
      DO i = 1,nx-1

        ! DTD: updated to allow calculation for diagnostic N0 version of Zhang

        !IF(mphyopt == 12 .or. mphyopt == 13) THEN
        !  n0r = 7.84e6*(rhobar(i,j,k)*1000.0 * MAX(0.0,qr(i,j,k)))**0.681
        !  arg = 3.07e9*MAX(0.0,(n0r**(-3./4.)))*( rhobar(i,j,k)*1000.0 * MAX(0.0,qr(i,j,k)) )**svnfrth
        !ELSE

          arg = coef*( rhobar(i,j,k)*1000.0 * MAX(0.0,qr(i,j,k)) )**svnfrth

        !END IF

        reflc(i,j,k) = 10.0*ALOG10( MAX(arg,1.0) )

      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE reflec_wr

!
!########################################################################
!########################################################################
!#########                                                      #########
!#########               SUBROUTINE REFLEC_FERRIER              #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE reflec_ferrier(nx,ny,nz, rho, qscalar, t, rff)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This subroutine estimates logarithmic radar reflectivity using
! equations customized for use with the ARPS Lin-Tao microphysics
! package.
!
! The equations (after some algebra to optimize code) are as follows:
!
! dBZ = 10*LOG10(Ze)
!
! Ze = Zer + Zes + Zeh (contributions from rain, snow, and hail).
!
!         k * 720                                 1.75
!   Zer = --------------------------- * (rho * qr)
!           1.75      0.75       1.75
!         pi     * N0r     * rhor
!
!   Zes = Zesnegf for "dry" snow (T < 0 C), or Zes = Zesposf for "wet"
!   snow (T > 0 C)
!
!                          2                     0.25
!             k * 720 * |K|                * rhos
!                          ice                                    1.75
!   Zesnegf = --------------------------------------- * (rho * qs)
!               1.75       2          0.75         2
!             pi     *  |K|      * N0s     * rhoice
!                          water
!
!             k * 720                                 1.75
!   Zesposf = --------------------------- * (rho * qs)
!               1.75      0.75       1.75
!             pi     * N0s     * rhos
!
!          /  k * 720                     \ 0.95             1.6625
!   Zeh = |   ---------------------------  |     * (rho * qh)
!          \    1.75      0.75       1.75 /
!           \ pi     * N0h     * rhoh    /
!
!-----------------------------------------------------------------------
!
! REFERENCES:
!
! Jahn, D., D. Weber, E. Kemp, and H. Neeman, 2000:  Evidence of
!   convective-induced turbulence outside the immediate storm region:
!   Part III.  CAPS report submitted to AlliedSignal/Honeywell, 37pp.
!
! Ferrier, B. S., W.-K. Tao, and J. Simpson, 1995:  A double-moment
!   multiple-phase four-class bulk ice scheme.  Part II:  Simulations
!   of convective storms in different large-scale environments and
!   comparisons with other bulk parameterizations.  J. Atmos. Sci.,
!   45, 3846-3879.
!
! McCumber, M., W.-K. Tao, J. Simpson, R. Penc, and S.-T. Soong, 1991:
!   Comparison of ice-phase microphysical parameterization schemes using
!   numerical simulations of tropical convection.  J. Appl. Meteor.,
!   30, 985-1004.
!
! Smith, P. L., 1984:  Equivalent radar reflectivity factors for snow
!   and ice particle.  J. Climate Appl. Meteor., 23, 1258-1260.
!
! Smith, P. L., Jr., C. G. Myers, and H. D. Orville, 1975:  Radar
!   reflectivity factor calculations in numerical cloud models using
!   bulk parameterization of precipitation.  J. Appl. Meteor., 14,
!   1156-1165.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Henry Neeman, Spring 2000.
!
! MODIFICATION HISTORY:
!
! Eric Kemp, 8 October 2001
! Reformatted code for ARPS Fortran 90 standard.
!
! Ming Xue, 16 Oct. 2001
!   Changed ni,nj,nk to nx,ny,nz. Change loop bounds. Changed the order
!   of argument list. Removed IF( ice == 0 ) check.
!
! 04/05/03 (Keith Brewster)
!   Clean-up of some unused parameters for clarity.
!
! 09/28/06 (Dan Dawson)
!   Allowed intercept parameters for rain/snow/hail to vary
!   based on values in history dump for the Lin scheme.
!   Added switches to distinguish between Lin scheme and
!   the original, diagnostic N0, and SCG versions of the WSM6 scheme.
!   The WSM6 scheme uses qg instead of qh, and the intercept parameters
!   and densities are different from the Lin scheme.  Also, the snow
!   intercept parameter is diagnosed in the WSM6 scheme.
!
!-----------------------------------------------------------------------
!
! Force explicit declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Include files.
!-----------------------------------------------------------------------

  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
! Declare arguments.
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nx,ny,nz ! Dimensions of grid

  REAL, INTENT(IN) :: rho(nx,ny,nz) ! Air density (kg m**-3)
!  REAL, INTENT(IN) :: qr(nx,ny,nz) ! Rain mixing ratio (kg kg**-1)
!  REAL, INTENT(IN) :: qs(nx,ny,nz) ! Snow mixing ratio (kg kg**-1)
!  REAL, INTENT(IN) :: qg(nx,ny,nz) ! Graupel mixing ratio (kg kg**-1)
!  REAL, INTENT(IN) :: qh(nx,ny,nz) ! Hail mixing ratio (kg kg**-1)
  REAL, INTENT(IN) :: qscalar(nx,ny,nz,nscalar) ! Microphysics scalar array
  REAL, INTENT(IN) :: t(nx,ny,nz) ! Temperature (K)

  REAL, INTENT(OUT) :: rff(nx,ny,nz) ! Reflectivity (dBZ)

!-----------------------------------------------------------------------
! Declare local parameters.
!-----------------------------------------------------------------------

  REAL,PARAMETER :: ki2 = 0.176 ! Dielectric factor for ice if other
                                !   than melted drop diameters are used.
  REAL,PARAMETER :: kw2=0.93 ! Dielectric factor for water.

  REAL,PARAMETER :: degKtoC=273.15 ! Conversion factor from degrees K to
                                   !   degrees C

  REAL,PARAMETER :: m3todBZ=1.0E+18 ! Conversion factor from m**3 to
                                    !   mm**6 m**-3.

  REAL,PARAMETER :: pi=3.1415926 ! Pi.

  REAL,PARAMETER :: pipowf=7.0/4.0 ! Power to which pi is raised.

!  REAL,PARAMETER :: N0r=8.0E+06 ! Intercept parameter in 1/(m^4) for rain.
!  REAL,PARAMETER :: N0s=3.0E+06 ! Intercept parameter in 1/(m^4) for snow.
!  REAL,PARAMETER :: N0h=4.0E+04 ! Intercept parameter in 1/(m^4) for hail.
  REAL :: N0r, N0s, N0h

  REAL,PARAMETER :: N0xpowf=3.0/4.0 ! Power to which N0r,N0s & N0h are
                                    !   raised.

  REAL,PARAMETER :: approxpow=0.95 ! Approximation power for hail
                                   !   integral.

  REAL,PARAMETER :: rqrpowf=7.0/4.0 ! Power to which product rho * qr
                                    !   is raised.
  REAL,PARAMETER :: rqsnpowf=7.0/4.0 ! Power to which product rho * qs
                                     !   is raised (dry snow).
  REAL,PARAMETER :: rqsppowf=7.0/4.0 ! Power to which product rho * qs
                                     !   is raised (wet snow).

  REAL,PARAMETER :: rqhpowf=(7.0/4.0)*approxpow ! Power to which product
                                                !   rho * qh is raised.

  REAL,PARAMETER :: rhoi=917.  ! Density of ice (kg m**-3)
  REAL,PARAMETER :: rhor=1000. ! Density of rain (kg m**-3)
!  REAL,PARAMETER :: rhos=100.  ! Density of snow (kg m**-3)
!  REAL,PARAMETER :: rhoh=913.  ! Density of hail (kg m**-3)
  REAL           :: rhos
  REAL           :: rhoh

  REAL,PARAMETER :: rhoipowf=2.0     ! Power to which rhoi is raised.
  REAL,PARAMETER :: rhospowf=1.0/4.0 ! Power to which rhos is raised.
  REAL,PARAMETER :: rhoxpowf=7.0/4.0 ! Power to which rhoh is raised.

  REAL,PARAMETER :: Zefact=720.0 ! Multiplier for Ze components.

  REAL,PARAMETER :: lg10mul=10.0 ! Log10 multiplier

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------

  REAL    :: rcomp,scomp,hcomp,sumcomp
  INTEGER :: i,j,k
  REAL    :: Zerf,Zesnegf,Zesposf,Zehf

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  rhos = rhosnow
  rhoh = rhohail

  IF (rhos <= 0.0 .OR. rhoh <= 0.0) THEN
    IF (myproc == 0) WRITE(6,'(/3(1x,a,/))')                            &
          'WARNING: rhosnow and rhohail are not initialized.',          &
          '         You may be reading earlier history files.',         &
          '         rhos and rhoh were reset as rhos = 100, rhoh = 913.'

    rhos = 100.
    rhoh = 913.

!    WRITE(6,'(/a,2(a,F7.2),a,/)')    &
!                  'Either rhosnow or rhohail is not initialized, ',     &
!                  'rhosnow = ',rhosnow,' rhohail = ',rhohail,           &
!                  'Program aborting ......'
!    CALL arpsstop('Wrong parameters, rhosnow or rhohail.',1)
  END IF

!-----------------------------------------------------------------------
! Determine which microphysics scheme is being used and set parameters
! accordingly.
!-----------------------------------------------------------------------

  IF (n0rain <= 0.0) THEN
    IF (myproc == 0) WRITE(6,'(/3(1x,a,/))')                            &
          'WARNING: n0rain is not initialized.',                        &
          '         You may be reading earlier history files.',         &
          '         n0rain was reset to the default value of 8.0e6'

    n0rain = 8.0e6
  END IF

  IF (n0snow <= 0.0) THEN
    IF (myproc == 0) WRITE(6,'(/3(1x,a,/))')                            &
          'WARNING: n0snow is not initialized.',                        &
          '         You may be reading earlier history files.',         &
          '         n0snow was reset to the default value of 3.0e6'

    n0snow = 3.0e6
  END IF

  IF (n0hail <= 0.0) THEN
    IF (myproc == 0) WRITE(6,'(/3(1x,a,/))')                            &
          'WARNING: n0hail is not initialized.',                        &
          '         You may be reading earlier history files.',         &
          '         n0hail was reset to the default value of 4.0e4'

    n0hail = 4.0e4
  END IF

  IF (rhosnow <= 0.0) THEN
    IF (myproc == 0) WRITE(6,'(/3(1x,a,/))')                            &
          'WARNING: rhosnow is not initialized.',                       &
          '         You may be reading earlier history files.',         &
          '         rhosnow was reset to the default value of 100.'

    rhosnow = 100.
  END IF
  IF (rhogrpl <= 0.0) THEN
    IF (myproc == 0) WRITE(6,'(/3(1x,a,/))')                            &
          'WARNING: rhogrpl is not initialized.',                       &
          '         You may be reading earlier history files.',         &
          '         rhogrpl was reset to the default value of 400.'

    rhogrpl = 100.
  END IF
  IF (rhohail <= 0.0) THEN
    IF (myproc == 0) WRITE(6,'(/3(1x,a,/))')                            &
          'WARNING: rhohail is not initialized.',                       &
          '         You may be reading earlier history files.',         &
          '         rhohail was reset to the default value of 913.'

    rhohail = 913.
  END IF

  IF(mphyopt == 2) THEN    ! Lin scheme
    N0r = n0rain
    N0s = n0snow
    N0h = n0hail
  ELSE IF(mphyopt == 5) THEN    ! original WSM6 scheme
    N0r = 8.0e6
    N0h = 4.0e6
    rhoh = 500.
    rhos = 100.
  ELSE IF(mphyopt == 7) THEN    ! diagnostic N0r WSM6 scheme
    N0h = 4.0e6
    rhoh = 500.
    rhos = 100.
  ELSE                     ! Schultz, Straka schemes : assumed same as default Lin for now
    N0r = 8.0e6
    N0s = 3.0e6
    N0h = 4.0e4
  END IF

!-----------------------------------------------------------------------
! First gather all the constants together.  (These are treated as
! variables because Fortran 90 does not allow real exponents when
! calculating parameters).
!-----------------------------------------------------------------------

  Zerf = (m3todBZ * Zefact) /  &
                  ((pi ** pipowf) * (N0r ** N0xpowf) *  &
                   (rhor ** rhoxpowf))

  IF(mphyopt == 2 .or. mphyopt == 3 .or. mphyopt == 4) THEN
  ! These schemes have constant N0 for snow
    Zesnegf = ((m3todBZ * Zefact   * Ki2 * (rhos ** rhospowf)) /  &
                   ((pi ** pipowf) * Kw2 * (N0s ** N0xpowf) *  &
                    (rhoi ** rhoipowf)))
    Zesposf = ((m3todBZ * Zefact) /  &
                   ((pi ** pipowf) * (N0s ** N0xpowf) *  &
                    (rhos ** rhoxpowf)))
  END IF

  Zehf = (((m3todBZ * Zefact) /  &
                    ((pi ** pipowf) * (N0h ** N0xpowf) *  &
                     (rhoh ** rhoxpowf))) ** approxpow)

!-----------------------------------------------------------------------
! Now loop through the scalar grid points.
!-----------------------------------------------------------------------

!  DO k = 2,nz-1      ! original
  DO k = 1,nz-1       ! Eric 8/8/03
    DO j = 1,ny-1
      DO i = 1,nx-1

!-----------------------------------------------------------------------
! Check for bad air density value.
!-----------------------------------------------------------------------

        IF (rho(i,j,k) <= 0.0) THEN
          rff(i,j,k) = 0.0
        ELSE

!-----------------------------------------------------------------------
! Calculate reflectivity contribution from rain.
!-----------------------------------------------------------------------

          IF (P_QR <= 0) THEN
            rcomp = 0.0
          ELSE IF(qscalar(i,j,k,P_QR) <= 0.0) THEN
            rcomp = 0.0
          ELSE IF(P_QR > 0) THEN
            IF(mphyopt >= 2 .and. mphyopt <= 5) THEN
              rcomp = Zerf*((qscalar(i,j,k,P_QR)*rho(i,j,k)) ** rqrpowf)
            ELSE IF(mphyopt == 6) THEN   ! WSM6 SCG for rain
              rcomp = -0.2645*(LOG10(qscalar(i,j,k,P_QR)*rho(i,j,k)*1000.) ** 4)     &
                      -1.43*(LOG10(qscalar(i,j,k,P_QR)*rho(i,j,k)*1000.) ** 3)       &
                      -0.8737*(LOG10(qscalar(i,j,k,P_QR)*rho(i,j,k)*1000.) ** 2)     &
                      +16.413*LOG10(qscalar(i,j,k,P_QR)*rho(i,j,k)*1000.)+43.1076
            ELSE IF(mphyopt == 7) THEN   ! WSM6 diagnostic N0 for rain
              rcomp = 2.074e4*((qscalar(i,j,k,P_QR)*rho(i,j,k)*1000.) ** 1.24)
            END IF
          END IF

!-----------------------------------------------------------------------
! Calculate reflectivity contribution from snow (dry or wet).
!-----------------------------------------------------------------------

          IF (P_QS <= 0) THEN
            scomp = 0.0
          ELSE IF(qscalar(i,j,k,P_QS) <= 0.0) THEN
            scomp = 0.0
          ELSE IF (P_QS > 0) THEN
            IF (mphyopt >= 5 .and. mphyopt <= 7) THEN   ! WSM6 schemes
              ! Calculate temperature dependent N0s based on Hong et al. (2004,MWR)
              ! eqn. 6

              N0s = min(2e6*exp(0.12*(273.15-t(i,j,k))),1e11)

              ! Now calculate constant in snow reflectivity equation

              Zesnegf = ((m3todBZ * Zefact   * Ki2 * (rhos ** rhospowf)) /  &
                   ((pi ** pipowf) * Kw2 * (N0s ** N0xpowf) *  &
                    (rhoi ** rhoipowf)))
              Zesposf = ((m3todBZ * Zefact) /  &
                   ((pi ** pipowf) * (N0s ** N0xpowf) *  &
                    (rhos ** rhoxpowf)))

            END IF

            IF (t(i,j,k) <= degKtoC) THEN
              scomp = Zesnegf*((qscalar(i,j,k,P_QS)*rho(i,j,k)) ** rqsnpowf)
            ELSE
              scomp = Zesposf*((qscalar(i,j,k,P_QS)*rho(i,j,k)) ** rqsppowf)
            END IF
          END IF

!-----------------------------------------------------------------------
! Calculate reflectivity contribution from hail.
!-----------------------------------------------------------------------

          IF (mphyopt >= 5 .and. mphyopt <= 7) THEN ! WSM6 scheme uses graupel
            ! First check to make sure the graupel array exists.  If reading
            ! a history file from an earlier version, the graupel array for the
            ! WSM6 scheme is sometimes found in the hail array instead
            IF(P_QH > 0) THEN ! Assume graupel array is in hail array
              IF (qscalar(i,j,k,P_QH) <= 0.0) THEN
                hcomp = 0.0
              ELSE
                hcomp = Zehf*((qscalar(i,j,k,P_QH)*rho(i,j,k)) ** rqhpowf)
              END IF
            ELSE IF (P_QG > 0) THEN ! Assume graupel array is defined (i.e. P_QG > 0)
              IF (qscalar(i,j,k,P_QG) <= 0.0) THEN
                hcomp = 0.0
              ELSE
                hcomp = Zehf*((qscalar(i,j,k,P_QG)*rho(i,j,k)) ** rqhpowf)
              END IF
            END IF
          ELSE
            IF (P_QH < 0) THEN
              hcomp = 0.0
            ELSE IF (qscalar(i,j,k,P_QH) <= 0.0) THEN
              hcomp = 0.0
            ELSE
              hcomp = Zehf*((qscalar(i,j,k,P_QH)*rho(i,j,k)) ** rqhpowf)
            END IF
          END IF

!-----------------------------------------------------------------------
! Now add the contributions and convert to logarithmic reflectivity
! factor dBZ.
!-----------------------------------------------------------------------

          IF (mphyopt /= 6) THEN
            sumcomp = rcomp + scomp + hcomp
            rff(i,j,k) = lg10mul * LOG10(MAX(sumcomp,1.0))
          ELSE
            sumcomp = scomp + hcomp
            rff(i,j,k) = lg10mul * LOG10(MAX(sumcomp,1.0))+rcomp
          END IF

        END IF !  IF (rho(i,j,k) <= 0.0) ... ELSE ...

      END DO ! DO i
    END DO ! DO j
  END DO ! DO k

  RETURN
END SUBROUTINE reflec_ferrier

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE REFLEC_MM                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE reflec_MM(nx,ny,nz, rho, qscalar, t, rff)
!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This subroutine estimates logarithmic radar reflectivity using
! equations based on the formulation described in Milbrandt and Yau (2005 Part I)
! It is designed for use with any of the 1, 2, or 3-moment versions of
! of the MY multi-moment microphysics scheme within ARPS.
!-----------------------------------------------------------------------
!
! AUTHOR:  Dan Dawson, Fall 2006.
!
! MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
! Force explicit declarations.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Include files.
!-----------------------------------------------------------------------

  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'mp.inc'
!-----------------------------------------------------------------------
! Declare arguments.
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nx,ny,nz ! Dimensions of grid

  REAL, INTENT(IN) :: rho(nx,ny,nz) ! Air density (kg m**-3)
  REAL, INTENT(IN) :: qscalar(nx,ny,nz,nscalar) ! Scalar microphysics array
  REAL, INTENT(IN) :: t(nx,ny,nz) ! Temperature (K)

  REAL, INTENT(OUT) :: rff(nx,ny,nz) ! Reflectivity (dBZ)

!-----------------------------------------------------------------------
! Declare local parameters.
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: nqscalar=6       ! Number of microphysics species (6 for MY)
  REAL,    PARAMETER :: pi = 3.14159265
  REAL,    PARAMETER :: piovr6 = pi/6.0
  REAL,    PARAMETER :: deratio = 0.224     ! Ratio of dielectric constants for ice and water
  REAL,    PARAMETER :: cr = (piovr6)*1000.  ! constant in mass power law relation for water particle
  REAL,    PARAMETER :: onethird = 1./3.

!  REAL, PARAMETER :: rhor = 1000        ! Density of liquid water
!  REAL, PARAMETER :: rhoi = 500         ! Density of ice particles
!  REAL, PARAMETER :: rhos = 100         ! Density of snow
!  REAL, PARAMETER :: rhog = 400         ! Density of graupel
!  REAL, PARAMETER :: rhoh = 900         ! Density of hail

  REAL, PARAMETER :: rhor = 1000        ! Density of liquid water
  REAL :: rhoi         ! Density of ice particles
  REAL :: rhos         ! Density of snow
  REAL :: rhog         ! Density of graupel
  REAL :: rhoh         ! Density of hail

!  REAL, PARAMETER :: N0rfix = 1.0e6
!  REAL, PARAMETER :: N0sfix = 1.0e7
!  REAL, PARAMETER :: N0gfix = 4.0e5
!  REAL, PARAMETER :: N0hfix = 1.0e5

  REAL :: N0rfix, N0sfix, N0gfix, N0hfix

  ! Constants in diagnostic alpha relations

  REAL, PARAMETER :: c1r = 19.0, c2r = 0.6, c3r = 1.8, c4r = 17.0
  REAL, PARAMETER :: c1i = 12.0, c2i = 0.7, c3i = 1.7, c4i = 11.0
  REAL, PARAMETER :: c1s = 4.5,  c2s = 0.5, c3s = 5.0, c4s = 5.5
  REAL, PARAMETER :: c1g = 5.5,  c2g = 0.7, c3g = 4.5, c4g = 8.5
  REAL, PARAMETER :: c1h = 3.7,  c2h = 0.3, c3h = 9.0, c4h = 6.5
  REAL, PARAMETER :: c5h = 1.0,  c6h = 6.5


  REAL :: Ntr
  REAL :: Nti
  REAL :: Nts
  REAL :: Ntg
  REAL :: Nth
  REAL :: Zr                 ! rain reflectivity
  REAL :: Zi                 ! ice reflectivity
  REAL :: Zs                 ! snow reflectivity
  REAL :: Zsd                ! dry snow reflectivity
  REAL :: Zsw                ! wet snow reflectivity
  REAL :: Zg                 ! graupel reflectivity
  REAL :: Zgd                ! dry graupel reflectivity
  REAL :: Zgw                ! wet graupel reflectivity
  REAL :: Zh                 ! hail reflectivity
  REAL :: Zhd                ! dry hail reflectivity
  REAL :: Zhw                ! wet hail reflectivity
  REAL :: Zt                 ! total reflectivity
  REAL :: alpha              ! shape parameter in microphysical gamma distribution
  REAL :: dmx                ! mean-mass diameter of distribution
  REAL :: Gx                 ! parameter in radar equation

  REAL, PARAMETER :: epsQ = 1.0e-14
  REAL, PARAMETER :: epsN = 1.0e-3
  REAL, PARAMETER :: epsZ = 1.0e-32

  ! Note: for now, these must be set to the same as the values in the original
  ! run of the 2-moment version, or the reflectivity values will be incorrect.
  ! EDIT 08/04/08: values now passed in through common block in phycst.inc

!  REAL, PARAMETER :: alpharfix = 0.0
!  REAL, PARAMETER :: alphaifix = 0.0
!  REAL, PARAMETER :: alphasfix = 0.0
!  REAL, PARAMETER :: alphagfix = 0.0
!  REAL, PARAMETER :: alphahfix = 0.0

  REAL :: alpharfix, alphaifix, alphasfix, alphagfix, alphahfix

  INTEGER :: i,j,k,nq

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  print*,'mphyopt = ',mphyopt

  IF (ntcloud <= 0.0) THEN
    IF (myproc == 0) WRITE(6,'(/3(1x,a,/))')                            &
          'WARNING: ntcloud is not initialized.',          &
          '         You may be reading earlier history files.',         &
          '         ntcloud was reset to the default value of 1.0e8'

    ntcloud = 1.0e8
  END IF
  IF (n0rain <= 0.0) THEN
    IF (myproc == 0) WRITE(6,'(/3(1x,a,/))')                            &
          'WARNING: n0rain is not initialized.',          &
          '         You may be reading earlier history files.',         &
          '         n0rain was reset to the default value of 8.0e6'

    n0rain = 8.0e6
  END IF
  IF (n0snow <= 0.0) THEN
    IF (myproc == 0) WRITE(6,'(/3(1x,a,/))')                            &
          'WARNING: n0snow is not initialized.',          &
          '         You may be reading earlier history files.',         &
          '         n0snow was reset to the default value of 3.0e6'

    n0snow = 3.0e6
  END IF
  IF (n0grpl <= 0.0) THEN
    IF (myproc == 0) WRITE(6,'(/3(1x,a,/))')                            &
          'WARNING: n0grpl is not initialized.',          &
          '         You may be reading earlier history files.',         &
          '         n0grpl was reset to the default value of 4.0e5'

    n0grpl = 4.0e5
  END IF
  IF (n0hail <= 0.0) THEN
    IF (myproc == 0) WRITE(6,'(/3(1x,a,/))')                            &
          'WARNING: n0hail is not initialized.',          &
          '         You may be reading earlier history files.',         &
          '         n0hail was reset to the default value of 4.0e4'

    n0hail = 4.0e4
  END IF
  IF (rhoice <= 0.0) THEN
    IF (myproc == 0) WRITE(6,'(/3(1x,a,/))')                            &
          'WARNING: rhoice is not initialized.',          &
          '         You may be reading earlier history files.',         &
          '         rhoice was reset to the default value of 500.'

    rhoice = 500.
  END IF
  IF (rhosnow <= 0.0) THEN
    IF (myproc == 0) WRITE(6,'(/3(1x,a,/))')                            &
          'WARNING: rhosnow is not initialized.',          &
          '         You may be reading earlier history files.',         &
          '         rhosnow was reset to the default value of 100.'

    rhosnow = 100.
  END IF
  IF (rhogrpl <= 0.0) THEN
    IF (myproc == 0) WRITE(6,'(/3(1x,a,/))')                            &
          'WARNING: rhogrpl is not initialized.',          &
          '         You may be reading earlier history files.',         &
          '         rhogrpl was reset to the default value of 400.'

    rhogrpl = 100.
  END IF
  IF (rhohail <= 0.0) THEN
    IF (myproc == 0) WRITE(6,'(/3(1x,a,/))')                            &
          'WARNING: rhohail is not initialized.',          &
          '         You may be reading earlier history files.',         &
          '         rhohail was reset to the default value of 913.'

    rhohail = 913.
  END IF

  rhoi = rhoice
  rhos = rhosnow
  rhog = rhogrpl
  rhoh = rhohail
  N0rfix = n0rain
  N0sfix = n0snow
  N0gfix = n0grpl
  N0hfix = n0hail
  alpharfix = alpharain
  alphaifix = alphaice
  alphasfix = alphasnow
  alphagfix = alphagrpl
  alphahfix = alphahail

!  print*,'rhoi,rhos,rhog,rhoh',rhoi,rhos,rhog,rhoh
!  print*,'N0rfix,N0sfix,N0gfix,N0hfix',N0rfix,N0sfix,N0gfix,N0hfix
!  print*,'alpharfix,alphaifix,alphagfix,alphahfix',alpharfix,alphaifix,alphasfix,alphagfix,alphahfix

  Zr = 0.0
  Zi = 0.0
  Zs = 0.0
  Zg = 0.0
  Zh = 0.0

  rff = 0.0

  DO k = 1,nz-1
    DO j = 1,ny-1
      DO i = 1,nx-1

      Zr = 0.0
      Zi = 0.0
      Zs = 0.0
      Zg = 0.0
      Zh = 0.0

        DO nq = 1,nqscalar
          IF(mphyopt == 11) THEN   ! 3-moment case: reflectivity predicted
            IF(qscalar(i,j,k,nq+2*nqscalar-1) >= epsZ) THEN
              IF(nq == P_QR) THEN
                Zr = (((pi/6.)*rhor/cr)**2.)*qscalar(i,j,k,nq+2*nqscalar-1)
                !print*,'rain block,qr,nr,Zr_raw,Zr',qscalar(i,j,k,nq),qscalar(i,j,k,nq+nqscalar),qscalar(i,j,k,nq+2*nqscalar-1),Zr
              ELSE IF(nq == P_QI) THEN
                Zi = deratio*((440./cr)**2.)*qscalar(i,j,k,nq+2*nqscalar-1)
              ELSE IF(nq == P_QS) THEN
                !IF(t(i,j,k) <= 273.15) THEN
                Zs = deratio*(((pi/6.)*rhos/cr)**2.)*qscalar(i,j,k,nq+2*nqscalar-1)
                !ELSE
                !  Zs = (((pi/6.)*rhos/cr)**2.)*qscalar(i,j,k,nq+2*nqscalar-1)
                !END IF
              ELSE IF(nq == P_QG) THEN
                Zgd = deratio*(((pi/6.)*rhog/cr)**2.)*qscalar(i,j,k,nq+2*nqscalar-1)
                !Zgw = (((pi/6.)*rhog/cr)**2.)*qscalar(i,j,k,nq+2*nqscalar-1)
                !IF(t(i,j,k) <= 275.65 .and. t(i,j,k) >= 270.65) THEN
                !  Zg = Zgw*(t(i,j,k)-270.65)/5.0 + Zgd*(275.65-t(i,j,k))/5.0
                !ELSE IF (t(i,j,k) < 270.65) THEN
                Zg = Zgd
                !ELSE
                !  Zg = Zgw
                !END IF
              ELSE IF(nq == P_QH) THEN
                Zhd = deratio*(((pi/6.)*rhoh/cr)**2.)*qscalar(i,j,k,nq+2*nqscalar-1)
                !Zhw = (((pi/6.)*rhoh/cr)**2.)*qscalar(i,j,k,nq+2*nqscalar-1)
                !IF(t(i,j,k) <= 275.65 .and. t(i,j,k) >= 270.65) THEN
                !  Zh = Zhw*(t(i,j,k)-270.65)/5.0 + Zhd*(275.65-t(i,j,k))/5.0
                !ELSE IF (t(i,j,k) < 270.65) THEN
                 Zh = Zhd
                !ELSE
                !  Zh = Zhw
                !END IF
              END IF
            END IF
          ELSE IF(mphyopt == 10 .OR. mphyopt == 12) THEN     ! 2-moment with diagnostic alpha
            IF(qscalar(i,j,k,nq) >= epsQ .and. qscalar(i,j,k,nq+nqscalar) >= epsN) THEN
              IF(nq == P_QR) THEN
                dmx = (rho(i,j,k)*qscalar(i,j,k,nq)                     &
                       /((pi/6.)*rhor*qscalar(i,j,k,nq+nqscalar)))**(1./3.)
                alpha = c1r*TANH(c2r*(dmx-c3r))+c4r
                Gx = (6+alpha)*(5+alpha)*(4+alpha)/((3+alpha)*(2+alpha)*(1+alpha))
                Zr = ((1./cr)**2.)*Gx                                   &
                   * ((rho(i,j,k)*qscalar(i,j,k,nq))**2.)/qscalar(i,j,k,nq+nqscalar)
              ELSE IF(nq == P_QI) THEN
                dmx = (rho(i,j,k)*qscalar(i,j,k,nq)/((pi/6.)*rhoi*qscalar(i,j,k,nq+nqscalar)))**(1./3.)
                alpha = c1i*TANH(c2i*(dmx-c3i))+c4i
                Gx = (6+alpha)*(5+alpha)*(4+alpha)/((3+alpha)*(2+alpha)*(1+alpha))
                Zi = deratio*((1/cr)**2.)*Gx                            &
                   * ((rho(i,j,k)*qscalar(i,j,k,nq))**2.)/qscalar(i,j,k,nq+nqscalar)
              ELSE IF(nq == P_QS) THEN
                dmx = (rho(i,j,k)*qscalar(i,j,k,nq)/((pi/6.)*rhos*qscalar(i,j,k,nq+nqscalar)))**(1./3.)
                alpha = c1s*TANH(c2s*(dmx-c3s))+c4s
                Gx = (6+alpha)*(5+alpha)*(4+alpha)/((3+alpha)*(2+alpha)*(1+alpha))
                Zs = deratio*((1/cr)**2.)*Gx                            &
                   * ((rho(i,j,k)*qscalar(i,j,k,nq))**2.)/qscalar(i,j,k,nq+nqscalar)
              ELSE IF(nq == P_QG) THEN
                dmx = (rho(i,j,k)*qscalar(i,j,k,nq)/((pi/6.)*rhog*qscalar(i,j,k,nq+nqscalar)))**(1./3.)
                alpha = c1g*TANH(c2g*(dmx-c3g))+c4g
                Gx = (6+alpha)*(5+alpha)*(4+alpha)/((3+alpha)*(2+alpha)*(1+alpha))
                Zg = deratio*((1/cr)**2.)*Gx                            &
                   * ((rho(i,j,k)*qscalar(i,j,k,nq))**2.)/qscalar(i,j,k,nq+nqscalar)
              ELSE IF(nq == P_QH) THEN
                dmx = (rho(i,j,k)*qscalar(i,j,k,nq)/((pi/6.)*rhoh*qscalar(i,j,k,nq+nqscalar)))**(1./3.)
                IF(dmx < 8e-3) THEN
                  alpha = c1h*TANH(c2h*(dmx-c3h))+c4h
                ELSE
                  alpha = c5h*dmx-c6h
                END IF
                Gx = (6+alpha)*(5+alpha)*(4+alpha)/((3+alpha)*(2+alpha)*(1+alpha))
                Zh = deratio*((1/cr)**2.)*Gx                            &
                   * ((rho(i,j,k)*qscalar(i,j,k,nq))**2.)/qscalar(i,j,k,nq+nqscalar)
              END IF
            END IF
          ELSE IF(mphyopt == 9) THEN     ! 2-moment scheme with fixed alpha
            IF(qscalar(i,j,k,nq) >= epsQ .and. qscalar(i,j,k,nq+nqscalar) >= epsN) THEN
              IF(nq == P_QR) THEN
                alpha = alpharfix
                Gx = (6+alpha)*(5+alpha)*(4+alpha)/((3+alpha)*(2+alpha)*(1+alpha))
                Zr = ((1./cr)**2.)*Gx                                   &
                   * ((rho(i,j,k)*qscalar(i,j,k,nq))**2.)/qscalar(i,j,k,nq+nqscalar)
                !print*,'rain block,i,j,k,alpha,Gx,rho,qr,nr',i,j,k,alpha,Gx,rho(i,j,k),qscalar(i,j,k,nq),qscalar(i,j,k,nq+nqscalar),10.0*LOG10(1.0e18*Zr)
              ELSE IF(nq == P_QI) THEN
                alpha = alphaifix
                Gx = (6+alpha)*(5+alpha)*(4+alpha)/((3+alpha)*(2+alpha)*(1+alpha))
                Zi = deratio*((1./cr)**2.)*Gx                           &
                   * ((rho(i,j,k)*qscalar(i,j,k,nq))**2.)/qscalar(i,j,k,nq+nqscalar)
                !  print*,'ice block'
              ELSE IF(nq == P_QS) THEN
                alpha = alphasfix
                Gx = (6+alpha)*(5+alpha)*(4+alpha)/((3+alpha)*(2+alpha)*(1+alpha))
                Zs = deratio*((1./cr)**2.)*Gx                           &
                   * ((rho(i,j,k)*qscalar(i,j,k,nq))**2.)/qscalar(i,j,k,nq+nqscalar)
                !  print*,'snow block'
              ELSE IF(nq == P_QG) THEN
                alpha = alphagfix
                Gx = (6+alpha)*(5+alpha)*(4+alpha)/((3+alpha)*(2+alpha)*(1+alpha))
                Zg = deratio*((1./cr)**2.)*Gx                           &
                   * ((rho(i,j,k)*qscalar(i,j,k,nq))**2.)/qscalar(i,j,k,nq+nqscalar)
                !  print*,'graupel block'
              ELSE IF(nq == P_QH) THEN
                alpha = alphahfix
                Gx = (6+alpha)*(5+alpha)*(4+alpha)/((3+alpha)*(2+alpha)*(1+alpha))
                Zh = deratio*((1./cr)**2.)*Gx                           &
                   * ((rho(i,j,k)*qscalar(i,j,k,nq))**2.)/qscalar(i,j,k,nq+nqscalar)
                !print*,'hail block'
              END IF
            END IF
          ELSE IF(mphyopt == 8) THEN     ! 1-moment scheme with fixed N0x (and alpha)
            IF(qscalar(i,j,k,nq) >= epsQ) THEN
              IF(nq == P_QR) THEN
                Ntr = N0rfix**(3./4.)*(rho(i,j,k)*qscalar(i,j,k,nq)/(pi*rhor))**(1./4.)
                alpha = alpharfix
                Gx = (6+alpha)*(5+alpha)*(4+alpha)/((3+alpha)*(2+alpha)*(1+alpha))
                IF(Ntr >= epsN) THEN
                  Zr = ((1/cr)**2.)*Gx*((rho(i,j,k)*qscalar(i,j,k,nq))**2.)/Ntr
                END IF
              ELSE IF(nq == P_QI) THEN
                Nti = 5.*exp(0.304*(273.15-max(233.,t(i,j,k))))
                alpha = alphaifix
                Gx = (6+alpha)*(5+alpha)*(4+alpha)/((3+alpha)*(2+alpha)*(1+alpha))
                IF(Nti >= epsN) THEN
                  Zi = deratio*((1/cr)**2.)*Gx*((rho(i,j,k)*qscalar(i,j,k,nq))**2.)/Nti
                END IF
              ELSE IF(nq == P_QS) THEN
                Nts = N0sfix**(3./4.)*(rho(i,j,k)*qscalar(i,j,k,nq)/(pi*rhos))**(1./4.)
                alpha = alphasfix
                Gx = (6+alpha)*(5+alpha)*(4+alpha)/((3+alpha)*(2+alpha)*(1+alpha))
                IF(Nts >= epsN) THEN
                  Zs = deratio*((1/cr)**2.)*Gx*((rho(i,j,k)*qscalar(i,j,k,nq))**2.)/Nts
                END IF
              ELSE IF(nq == P_QG) THEN
                Ntg = N0gfix**(3./4.)*(rho(i,j,k)*qscalar(i,j,k,nq)/(pi*rhog))**(1./4.)
                alpha = alphagfix
                Gx = (6+alpha)*(5+alpha)*(4+alpha)/((3+alpha)*(2+alpha)*(1+alpha))
                IF(Ntg >= epsN) THEN
                  Zg = deratio*((1/cr)**2.)*Gx*((rho(i,j,k)*qscalar(i,j,k,nq))**2.)/Ntg
                END IF
              ELSE IF(nq == P_QH) THEN
                Nth = N0hfix**(3./4.)*(rho(i,j,k)*qscalar(i,j,k,nq)/(pi*rhoh))**(1./4.)
                alpha = alphahfix
                Gx = (6+alpha)*(5+alpha)*(4+alpha)/((3+alpha)*(2+alpha)*(1+alpha))
                IF(Nth >= epsN) THEN
                  Zh = deratio*((1/cr)**2.)*Gx*((rho(i,j,k)*qscalar(i,j,k,nq))**2.)/Nth
                END IF
              END IF
            END IF
          END IF
        END DO  ! nq loop

        ! Sum up contributions from all hydrometeor species and convert
        ! to logarithmic reflectivity

        Zt = 1.0e18*(Zr+Zi+Zs+Zg+Zh)   ! Now in units of mm^6*m^-3

        !IF(Zt >= epsZ) THEN
        !  print*,'Zt = ',Zt
        !END IF
        rff(i,j,k) = 10.0*LOG10(MAX(Zt,1.0))   ! Now in units of dBZ


        !IF(rff(i,j,k) >= 0.0) THEN
        !  print*,rff(i,j,k)
        !END IF

      END DO  ! k loop
    END DO  ! j loop
  END DO  ! i loop

  RETURN
END SUBROUTINE reflec_MM
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SETCORNERLL                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE setcornerll( nx,ny, x, y )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the latitudes and longitudes set all corner points.
!
!  Before calling this subroutine, the map projection should have
!  been set up.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!          09/30/1997
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points for the model
!             grid in the east-west direction.
!    ny       Number of grid points for the model
!             grid in the north-south direction.
!
!    x        Analysis grid points in the e-w direction
!             (in grid units)
!    y        Analysis grid points in the n-s direction
!             (in grid units)
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx                ! Number of model grid points
                               ! in the east-west direction.
  INTEGER :: ny                ! Number of model grid points
                               ! in the north-south direction
  REAL :: x(nx)                ! 2-D model grid points east-west
                               ! direction (model grid units)
  REAL :: y(ny)                ! 2-D model grid points north-south
                               ! direction (model grid units)
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  REAL :: tema1, tema2, temb1, temb2
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  tema1 = 0.5*(x(1)+x(2))
  tema2 = 0.5*(x(nx-1)+x(nx)) + dx  ! The nx'th scalar point

  temb1 = 0.5*(y(1)+y(2))
  temb2 = 0.5*(y(ny-1)+y(ny)) + dy  ! The ny'th scalar point

  CALL xytoll(1,1, tema1,temb1, swlats,swlons) ! for scalar grid
  CALL xytoll(1,1, tema2,temb2, nelats,nelons) ! for scalar grid
  CALL xytoll(1,1, x(1), temb1, swlatu,swlonu) ! for u-wind grid
  CALL xytoll(1,1, x(nx),temb2, nelatu,nelonu) ! for u-wind grid
  CALL xytoll(1,1, tema1,y(1),  swlatv,swlonv) ! for v-wind grid
  CALL xytoll(1,1, tema2,y(ny), nelatv,nelonv) ! for v-wind grid

  RETURN
END SUBROUTINE setcornerll

SUBROUTINE wrtvar(nx,ny,nz, var, varnam, time, runname,dirname)
!
!-----------------------------------------------------------------------
!
!  Purpose:
!  To write an array 'var' into a file in binary format.
!
!  Author: Ming Xue
!
!  Input:
!
!  nx, ny, nz   Dimensions of input array 'var'.
!  var          Input array to be written out.
!  varnam       String of length 6 to designate the input array
!  time         The time of the input data array in seconds.
!
!  Output:
!
!  A disk file containing array 'var'.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz
  REAL :: var(nx,ny,nz)
  CHARACTER (LEN=6) :: varnam
  CHARACTER (LEN=*) :: runname
  CHARACTER (LEN=*) :: dirname

  CHARACTER (LEN=80 ) :: timsnd
  CHARACTER (LEN=256) :: vfnam
  CHARACTER (LEN=256) :: temchar
  INTEGER :: lvfnam, ierr, istat, nunit
  INTEGER :: tmstrln
  INTEGER :: lrunnam,ldirnam
  REAL :: time

  INTEGER :: i,j,k
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  CALL cvttsnd( time, timsnd, tmstrln )
  CALL gtlfnkey( runname, lrunnam )

  vfnam = runname(1:lrunnam)//'.'//varnam//timsnd(1:tmstrln)

  lvfnam = len_trim(vfnam)
  ldirnam = LEN_trim(dirname)

  IF( ldirnam == 0 ) THEN
    dirname = '.'
    ldirnam = 1
  END IF

  PRINT*,'ldirnam,dirname=',ldirnam,dirname

  IF( dirname /= ' ' ) THEN

    temchar = vfnam
    vfnam = dirname(1:ldirnam)//'/'//temchar
    lvfnam = lvfnam + ldirnam + 1

  END IF

  PRINT*,'lvfnam,vfnam=',lvfnam,vfnam
  PRINT*,'lvfnam,vfnam=',lvfnam,vfnam(1:lvfnam)

  CALL getunit( nunit)

  CALL asnctl ('NEWLOCAL', 1, ierr)
  CALL asnfile(vfnam, '-F f77 -N ieee', ierr)

  OPEN(UNIT=nunit,FILE=trim(vfnam(1:lvfnam)),STATUS='unknown',          &
         FORM='unformatted',IOSTAT= istat )

  WRITE(nunit) varnam
  WRITE(nunit) nx,ny, nz
  DO k=1,nz
    WRITE(nunit) ((var(i,j,k),i=1,nx),j=1,ny)
  END DO

  CLOSE(UNIT=nunit)
  CALL retunit(nunit)

  RETURN
END SUBROUTINE wrtvar



SUBROUTINE readvar(nx,ny,nz, varnam, time, var, runname)

!-----------------------------------------------------------------------
!
!  Purpose:
!  To read in array 'var' from a file.
!
!  Author: Ming Xue
!
!  Modifications:
!  2/8/1998 (M.Xue)
!
!  Added check on the arrays bounds read in from data.
!
!  Input:
!
!  nx, ny, nz   Dimensions of array 'var' to be read in.
!  varnam       String of length 6 to designate the array
!  time         The time of the array in seconds.
!
!  Output:
!
!  var          Array readin from a disk file.
!
!-----------------------------------------------------------------------


  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz
  CHARACTER (LEN=6) :: varnam_in
  CHARACTER (LEN=*) :: varnam
  CHARACTER (LEN=*) :: runname
  REAL :: time

  REAL :: var(nx,ny,nz)

  CHARACTER (LEN=80 ) :: timsnd
  CHARACTER (LEN=256) :: vfnam
  INTEGER :: ierr, istat, nunit, lvar
  INTEGER :: tmstrln
  INTEGER :: lrunnam

  INTEGER :: nx_in,ny_in,nz_in

  INTEGER :: i,j,k
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  CALL cvttsnd( time, timsnd, tmstrln )
  CALL gtlfnkey( runname, lrunnam )

  lvar = LEN(varnam)
  vfnam = runname(1:lrunnam)//'.'//varnam(1:lvar)//timsnd(1:tmstrln)

  CALL getunit( nunit)

  CALL asnctl ('NEWLOCAL', 1, ierr)
  CALL asnfile(vfnam, '-F f77 -N ieee', ierr)

  WRITE(6,'(''READING ARRAY '',A,'' FROM FILE '',A)')                   &
        varnam(1:lvar),vfnam

  OPEN(UNIT=nunit,FILE=trim(vfnam),STATUS='old',                        &
         FORM='unformatted',IOSTAT= istat )

  READ(nunit) varnam_in
  READ(nunit) nx_in,ny_in,nz_in

  IF(nx_in /= nx .OR. ny_in /= ny .OR. nz_in /= nz) THEN
    WRITE(6,'(a,/a,a,/a,3I5,/a,3I5)')                                   &
        'Warning in subroutine READVAR: Dimensions of data file ',      &
        vfnam,' do not agree with the expected dimensions.',            &
        'nx, ny and nz in the data are ',nx_in,ny_in,nz_in,             &
        'nx, ny and nz expected    are ',nx,ny,nz
  END IF

  DO k=1,nz
    READ (nunit) ((var(i,j,k),i=1,nx),j=1,ny)
  END DO

  CLOSE(UNIT=nunit)
  CALL retunit(nunit)

  WRITE(6,'(''ARRAY '',A,'' READ FROM FILE '',A)')                      &
        varnam_in,vfnam

  RETURN
END SUBROUTINE readvar
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRTVAR1                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!


SUBROUTINE wrtvar1(nx,ny,nz,var, varid,varname,varunits,                &
           time,runname,dirname, istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!    To write an array 'var' into a file in binary format.
!
!  AUTHOR: Ming Xue
!
!  MODIFICATION HISTORY:
!
!  1998/03/17 (R. Carpenter)
!  Introduced I/O status variable and array name and units.
!
!  INPUT:
!
!  nx, ny, nz   Dimensions of input array 'var'.
!  var          Array to be written out.
!  varid        String of length 6 (padded with _ as necessary) to
!               designate the input array. (E.g., 'w_____')
!  varname      String describing the field (e.g., 'Vertical Velocity')
!  varunits     String describing the units (e.g., 'm/s')
!  time         The model time in seconds.
!  runname      Run name
!  dirname      Output directory (use '.' for current directory)
!
!  OUTPUT:
!
!  istatus       Exit status (0=okay, 1=warning, 2=error)
!
!  I/O:
!
!  An unformatted binary file named  dirname/runname.varid{time}  is created:
!             nx, ny, nz
!             var
!             varname
!             varunits
!
!-----------------------------------------------------------------------
!

  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz, istatus
  REAL :: var(nx,ny,nz)
  CHARACTER (LEN=* ) :: varunits
  CHARACTER (LEN=* ) :: varname
  CHARACTER (LEN=40) :: varunits1
  CHARACTER (LEN=40) :: varname1
  CHARACTER (LEN=6 ) :: varid

  CHARACTER (LEN=256) :: tmpdirnam
  CHARACTER (LEN=256) :: vfnam
  INTEGER :: ierr, IOSTAT, nunit
  REAL :: time
  CHARACTER (LEN=*  ) :: dirname
  CHARACTER (LEN=*  ) :: runname
  CHARACTER (LEN=256) :: timestring
  INTEGER :: lfnkey
  INTEGER :: ltimestring, lvfnam
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IOSTAT = 0
  iSTATUS = 0

  CALL gtlfnkey( runname, lfnkey )

  tmpdirnam = dirname
  IF ( LEN(tmpdirnam) == 0 .OR. tmpdirnam == ' ' ) THEN
    tmpdirnam = '.'
  END IF

  WRITE (vfnam,'(4A,A6)') trim(tmpdirnam),'/',runname(1:lfnkey),'.',varid

  CALL cvttsnd( time, timestring, ltimestring )
  lvfnam = len_trim(vfnam)

  vfnam(lvfnam+1:lvfnam+ltimestring) = timestring(1:ltimestring)

  PRINT *, 'WRTVAR1: Writing array to ', trim(vfnam)

  CALL getunit( nunit)

  CALL asnctl ('NEWLOCAL', 1, ierr)
  CALL asnfile(vfnam, '-F f77 -N ieee', ierr)

  OPEN(UNIT=nunit,FILE=vfnam,STATUS='unknown',FORM='unformatted',       &
      ERR=9000, IOSTAT=IOSTAT)

  WRITE(nunit, ERR=9000, IOSTAT=IOSTAT) nx,ny,nz
  WRITE(nunit, ERR=9000, IOSTAT=IOSTAT) var

  varname1 = varname
  varunits1= varunits
  WRITE(nunit, ERR=9000, IOSTAT=IOSTAT) varname1
  WRITE(nunit, ERR=9000, IOSTAT=IOSTAT) varunits1

  CLOSE(UNIT=nunit)
  CALL retunit(nunit)

  RETURN        ! Normal return

! I/O error handling
! Note that IOSTAT < 0 should not occur in this subroutine.

  9000 CONTINUE        ! I/O errors

  CLOSE(UNIT=nunit)
  CALL retunit(nunit)

  IF (IOSTAT < 0) THEN
    iSTATUS = 2
    PRINT *,                                                            &
        'WRTVAR1: I/O ERRORS OCCURRED ',                                &
        '(possible end of record or file): ',                           &
        IOSTAT, iSTATUS
  ELSE IF (IOSTAT > 0) THEN
    iSTATUS = 2
    PRINT *, 'WRTVAR1: UNRECOVERABLE I/O ERRORS OCCURRED: ',            &
        IOSTAT, iSTATUS
  END IF

  RETURN
END SUBROUTINE wrtvar1
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE READVAR1                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE readvar1(nx,ny,nz,var, varid,varname,varunits,               &
           time,runname,dirname, iSTATUS)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!  To read in array 'var' from a file.
!
!  AUTHOR: Ming Xue
!
!  MODIFICATION HISTORY:
!
!  2/8/1998 (M.Xue)
!  Added check on the arrays bounds read in from data.
!
!  1998/03/17 (R. Carpenter)
!  Introduced I/O istatus variable and array name and units.
!
!  INPUT:
!
!  nx, ny, nz   Dimensions of input array 'var'.
!  varid        String of length 6 (padded with _ as necessary) to
!               designate the input array. (E.g., 'w_____')
!  time         The model time in seconds.
!  runname      Run name
!  dirname      Input directory (use '.' for current directory)
!
!  OUTPUT:
!
!  var          Array to be written out.
!  varname      String describing the field (e.g., 'Vertical Velocity')
!  varunits     String describing the units (e.g., 'm/s')
!  istatus       Exit status (0=okay, 1=warning, 2=error)
!
!  I/O:
!
!  An unformatted binary file named  dirname/runname.varid{time}  is read:
!             nx, ny, nz
!             var
!             varname
!             varunits
!
!-----------------------------------------------------------------------


  IMPLICIT NONE

  INTEGER :: nx,ny,nz, iSTATUS
  CHARACTER (LEN=6) :: varid
  CHARACTER (LEN=*) :: varname, varunits
  CHARACTER (LEN=40) :: varname1, varunits1
  CHARACTER (LEN=*) :: runname, dirname
  REAL :: time
  REAL :: var(nx,ny,nz)

  CHARACTER (LEN=256) :: tmpdirnam
  CHARACTER (LEN=256) :: vfnam
  CHARACTER (LEN=256 ) :: timestring
  INTEGER :: ierr, IOSTAT, nunit

  INTEGER :: nx_in,ny_in,nz_in, lfnkey

  INTEGER :: lvfnam, ltimestring
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  iSTATUS = 0

  CALL gtlfnkey( runname, lfnkey )

  tmpdirnam = dirname
  IF ( LEN(tmpdirnam) == 0 .OR. tmpdirnam == ' ' ) THEN
    tmpdirnam = '.'
  END IF

  WRITE (vfnam,'(4A,A6)') trim(tmpdirnam),'/',runname(1:lfnkey),'.',varid

  CALL cvttsnd( time, timestring, ltimestring )
  lvfnam = len_trim(vfnam)
  vfnam(lvfnam+1:lvfnam+ltimestring) = timestring(1:ltimestring)

  PRINT *, 'READVAR1: Reading array from ',trim(vfnam)

  CALL getunit( nunit)

  CALL asnctl ('NEWLOCAL', 1, ierr)
  CALL asnfile(vfnam, '-F f77 -N ieee', ierr)

  OPEN(UNIT=nunit,FILE=trim(vfnam),STATUS='old', FORM='unformatted',    &
      ERR=9000, IOSTAT=IOSTAT)

  READ(nunit, ERR=9000, END=9000, IOSTAT=IOSTAT) nx_in,ny_in,nz_in

  IF(nx_in /= nx .OR. ny_in /= ny .OR. nz_in /= nz) THEN
    WRITE(6,'(a,/a,a,/a,3I5,/a,3I5)')                                   &
        'Warning in subroutine READVAR1: Dimensions of data file ',     &
        vfnam,' do not agree with the expected dimensions.',            &
        'nx, ny and nz in the data are ',nx_in,ny_in,nz_in,             &
        'nx, ny and nz expected    are ',nx,ny,nz
  END IF

  READ(nunit, ERR=9000, END=9000, IOSTAT=IOSTAT) var
  READ(nunit, ERR=9000,END=9000,IOSTAT=IOSTAT) varname1
  READ(nunit, ERR=9000,END=9000,IOSTAT=IOSTAT) varunits1

  varname = varname1
  varunits = varunits1

  CLOSE(UNIT=nunit)
  CALL retunit(nunit)

  RETURN        ! Normal return

! I/O error handling
! Note that warning (istatus=1) (e.g., end of record or file) is implemented as
! error (status=2) because of ambiguities in IOSTAT.

  9000 CONTINUE        ! I/O errors

  CLOSE(UNIT=nunit)
  CALL retunit(nunit)

  IF (IOSTAT < 0) THEN
    iSTATUS = 2
    PRINT *, 'READVAR1: I/O ERRORS OCCURRED ',                          &
        '(possible end of record or file): ',IOSTAT, iSTATUS
  ELSE IF (IOSTAT > 0) THEN
    iSTATUS = 2
    PRINT *, 'READVAR1: UNRECOVERABLE I/O ERRORS OCCURRED: ',           &
        IOSTAT,iSTATUS
  END IF

  RETURN
END SUBROUTINE readvar1
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRTVAR2                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE wrtvar2(nx,ny,nz,var,varid,varname,varunits,time,runname,    &
                   dirname,foutfmt,hdfcompr,mpiflag,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!    To write an array 'var' into a file. The file can be in three format
!    based on the argument, foutfmt:
!       foutfmt = 1,      Binary file
!       foutfmt = 3,      HDF4 file,
!         hdfcompr = 0-7, Specify HDF4 compression option, see arps.input
!       foutfmt = 7,      NetCDF file
!       others,           Error message
!
!  NOTE:
!    This subroutine includes "mp.inc", So all variables in mp.inc must be
!    set before calling this subroutine. It also needs to link with the
!    following objects:
!
!           alloclib.o
!           mpisubs.o          or       nompsubs.o
!           hdfio3d.o                   nohdfio3d.o
!           netio3d.o                   nonetio3d.o
!
!  AUTHOR: Yunheng Wang (2005/06/14)
!  Based on wrtvar1. The interface is the same as wrtvar1 except that
!  two more arguments, foutfmt, hdfcompr were added.
!
!  MODIFICATION HISTORY:
!
!  INPUT:
!
!    nx, ny, nz   Dimensions of input array 'var'.
!    var          Array to be written out.
!    varid        String of length 6 (padded with _ as necessary) to
!                 designate the input array. (E.g., 'w_____')
!    varname      String describing the field (e.g., 'Vertical Velocity')
!    varunits     String describing the units (e.g., 'm/s')
!    time         The model time in seconds.
!    runname      Run name
!    dirname      Output directory (use '.' for current directory)
!    foutfmt      File output format, 1,3, or 7
!    hdfcompr     HDF file compression option, only used when foutfmt=3
!    mpiflag      Flag to indicate whether to merge the array or not
!                 = 0, Do not merge
!                 = 1, merge if the program is in MPI mode.
!
!  OUTPUT:
!
!    istatus       Exit status (0=okay, 1=unknown file format, 2=error)
!
!  I/O:
!
!  An file named  dirname/runname.xxxvarid{time}  is created:
!             nx, ny, nz
!             var
!             varname
!             varunits
!
!  where xxx is
!          bin        for foutfmt = 1
!          hdf        for foutfmt = 3
!          net        for foutfmt = 7
!
!-----------------------------------------------------------------------
!

  IMPLICIT NONE

  INTEGER,          INTENT(IN)  :: nx,ny,nz
  REAL,             INTENT(IN)  :: var(nx,ny,nz)
  CHARACTER(LEN=6), INTENT(IN)  :: varid
  CHARACTER(LEN=*), INTENT(IN)  :: varunits
  CHARACTER(LEN=*), INTENT(IN)  :: varname
  CHARACTER(LEN=*), INTENT(IN)  :: dirname
  CHARACTER(LEN=*), INTENT(IN)  :: runname
  INTEGER,          INTENT(IN)  :: foutfmt
  INTEGER,          INTENT(IN)  :: hdfcompr
  INTEGER,          INTENT(IN)  :: mpiflag
  REAL,             INTENT(IN)  :: time
  INTEGER,          INTENT(OUT) :: istatus

  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  CHARACTER(LEN=40)  :: varunits1
  CHARACTER(LEN=40)  :: varname1
  CHARACTER(LEN=256) :: tmpstr
  CHARACTER(LEN=256) :: vfnam
  CHARACTER(LEN=256) :: timestring
  INTEGER            :: ierr, nunit
  INTEGER            :: lfnkey
  INTEGER            :: ltimestring, lvfnam

  CHARACTER(LEN=3), PARAMETER :: fmtstr(7) =                            &
                          (/'bin','xxx','hdf','xxx','xxx','xxx','net'/)

  INTEGER               :: nxlg, nylg
  REAL,    ALLOCATABLE  :: varout(:,:,:)

  INTEGER(KIND=selected_int_kind(4)), ALLOCATABLE  :: itmp(:,:,:)
  REAL,                               ALLOCATABLE  :: hmax(:), hmin(:)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  istatus = 0
  vfnam = ' '

  CALL gtlfnkey( runname, lfnkey )

  tmpstr = dirname
  IF ( LEN_TRIM(tmpstr) == 0 .OR. tmpstr == ' ' ) THEN
    tmpstr = '.'
  END IF

  WRITE (vfnam,'(4A,A3,A6)') trim(tmpstr),'/',runname(1:lfnkey),'.', &
                             fmtstr(foutfmt),varid

  CALL cvttsnd( time, timestring, ltimestring )
  lvfnam = len_trim(vfnam)

  vfnam(lvfnam+1:lvfnam+ltimestring) = timestring(1:ltimestring)
  lvfnam = lvfnam + ltimestring

!
! Merge variables if needed
!

  IF (mp_opt > 0 .AND. mpiflag > 0) THEN
    nxlg = (nx-3)*nproc_x + 3
    nylg = (ny-3)*nproc_y + 3
  ELSE
   IF (mp_opt > 0) THEN
     tmpstr = vfnam
     CALL gtsplitfn(tmpstr,1,1,loc_x,loc_y,1,1,0,0,0,2,vfnam,istatus)
     lvfnam = LEN_TRIM(vfnam)
   END IF
   nxlg = nx
   nylg = ny
  END IF

  CALL fnversn(vfnam,lvfnam)

  ALLOCATE(varout(nxlg,nylg,nz), STAT= istatus)

  IF (mp_opt > 0 .AND. mpiflag > 0) THEN
    CALL mpimerge3d(var,nx,ny,nz,varout)
  ELSE
    varout(:,:,:) = var(:,:,:)
  END IF

  !
  ! Write variable in file
  !
  IF (mpiflag == 0 .OR. myproc == 0) PRINT *, 'WRTVAR2: Writing array to ', trim(vfnam)

  varname1 = varname
  varunits1= varunits

  IF (myproc == 0 .OR. mpiflag == 0) THEN
    IF (foutfmt == 1) THEN
      CALL getunit( nunit)

      CALL asnctl ('NEWLOCAL', 1, ierr)
      CALL asnfile(vfnam, '-F f77 -N ieee', ierr)

      OPEN(UNIT=nunit,FILE=vfnam,STATUS='unknown',FORM='unformatted',       &
                      ERR=9000, IOSTAT=istatus)

      WRITE(nunit, ERR=9000, IOSTAT=istatus) nxlg,nylg,nz
      WRITE(nunit, ERR=9000, IOSTAT=istatus) varout

      WRITE(nunit, ERR=9000, IOSTAT=istatus) varname1
      WRITE(nunit, ERR=9000, IOSTAT=istatus) varunits1

      CLOSE(UNIT=nunit)
      CALL retunit(nunit)

      GO TO 8000

      ! I/O error handling
      ! Note that IOSTAT < 0 should not occur in this subroutine.

      9000 CONTINUE        ! I/O errors

      CLOSE(UNIT=nunit)
      CALL retunit(nunit)

      IF (istatus < 0) THEN
        PRINT *, 'WRTVAR1: I/O ERRORS OCCURRED ',                         &
                 '(possible end of record or file): ', istatus
      ELSE IF (istatus > 0) THEN
        PRINT *, 'WRTVAR1: UNRECOVERABLE I/O ERRORS OCCURRED: ',istatus
      END IF
      istatus = 2

      8000 CONTINUE

    ELSE IF (foutfmt == 3) THEN

      IF (hdfcompr > 3) THEN
        ALLOCATE (itmp(nxlg,nylg,nz),stat=istatus)
        CALL check_alloc_status(istatus, 'wrtvar2:itmp')
        ALLOCATE (hmax(nz),stat=istatus)
        CALL check_alloc_status(istatus, 'wrtvar2:hmax')
        ALLOCATE (hmin(nz),stat=istatus)
        CALL check_alloc_status(istatus, 'wrtvar2:hmin')
      END IF

      CALL hdfopen(vfnam,2,nunit)
      CALL hdfwrti(nunit, 'nx', nxlg, istatus)
      CALL hdfwrti(nunit, 'ny', nylg, istatus)
      CALL hdfwrti(nunit, 'nz', nz,   istatus)
      CALL hdfwrt3d(varout,nxlg,nylg,nz,nunit,0,hdfcompr,               &
                    varid,varname1,varunits1,itmp,hmax,hmin)
      CALL hdfclose(nunit,istatus)

      IF (hdfcompr > 3) THEN
        DEALLOCATE(itmp)
        DEALLOCATE(hmax,hmin)
      END IF

    ELSE IF (foutfmt == 7) THEN

      CALL netopen(vfnam,'C',nunit)
      CALL net_define_onevar(nunit,nxlg,nylg,nz,varid,varname1,varunits1,istatus)
      CALL netwrt3d(nunit,0,0,varid,varout,nxlg,nylg,nz)
      CALL netclose(nunit)

    ELSE
      WRITE(6,'(1x,a,I2,a,/,a)') 'Unsupported file format - ',foutfmt,    &
               ' inside wrtvar2.','Please check the subroutine arguments.'
      istatus = 1
    END IF
  END IF

  CALL mpupdatei(istatus, 1)

  DEALLOCATE(varout)

  RETURN        ! Normal return
END SUBROUTINE wrtvar2
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRTCOL2GRD                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
 SUBROUTINE wrtcol2grd(nx,ny,nz,nxny,ncolp,                             &
                       colvar,icolp,jcolp,                              &
                       varid,varname,varunits,varmiss,time,runname,     &
                       dirname,foutfmt,hdfcompr,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!    To write a radar column data array 'var' into a file as a 3d grid.
!    The file can be written in one three formats
!     based on the argument, foutfmt:
!       foutfmt = 1,      Binary file
!       foutfmt = 3,      HDF4 file,
!         hdfcompr = 0-7, Specify HDF4 compression option, see arps.input
!       foutfmt = 7,      NetCDF file
!       others,           Error message
!
!  NOTE:
!    This subroutine includes "mp.inc", So all variables in mp.inc must be
!    set before calling this subroutine. It also needs to link with the
!    following objects:
!
!           alloclib.o
!           mpisubs.o          or       nompsubs.o
!           hdfio3d.o                   nohdfio3d.o
!           netio3d.o                   nonetio3d.o
!
!  AUTHOR: Keith Brewster (03/09/2010)
!  Based on wrtvar2.
!
!  MODIFICATION HISTORY:
!
!  INPUT:
!
!    nx, ny, nz   Dimensions of the ARPS 3d grid
!    ncolp        Number of columns (for this processor) of input array 'var'.
!    var          Array to be written out.
!    varid        String of length 6 (padded with _ as necessary) to
!                 designate the input array. (E.g., 'w_____')
!    varname      String describing the field (e.g., 'Vertical Velocity')
!    varunits     String describing the units (e.g., 'm/s')
!    varmiss      Missing value
!    time         The model time in seconds.
!    runname      Run name
!    dirname      Output directory (use '.' for current directory)
!    foutfmt      File output format, 1,3, or 7
!    hdfcompr     HDF file compression option, only used when foutfmt=3
!
!  OUTPUT:
!
!    istatus       Exit status (0=okay, 1=unknown file format, 2=error)
!
!  I/O:
!
!  An file named  dirname/runname.xxxvarid{time}  is created:
!             nx, ny, nz
!             var
!             varname
!             varunits
!
!  where xxx is
!          bin        for foutfmt = 1
!          hdf        for foutfmt = 3
!          net        for foutfmt = 7
!
!-----------------------------------------------------------------------
!

  IMPLICIT NONE

  INTEGER,          INTENT(IN)  :: nx,ny,nz
  INTEGER,          INTENT(IN)  :: nxny,ncolp
  REAL,             INTENT(IN)  :: colvar(nz,nxny)
  INTEGER,          INTENT(IN)  :: icolp(nxny)
  INTEGER,          INTENT(IN)  :: jcolp(nxny)
  CHARACTER(LEN=6), INTENT(IN)  :: varid
  CHARACTER(LEN=*), INTENT(IN)  :: varunits
  REAL,             INTENT(IN)  :: varmiss
  CHARACTER(LEN=*), INTENT(IN)  :: varname
  CHARACTER(LEN=*), INTENT(IN)  :: dirname
  CHARACTER(LEN=*), INTENT(IN)  :: runname
  INTEGER,          INTENT(IN)  :: foutfmt
  INTEGER,          INTENT(IN)  :: hdfcompr
  REAL,             INTENT(IN)  :: time
  INTEGER,          INTENT(OUT) :: istatus

  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  CHARACTER(LEN=40)  :: varunits1
  CHARACTER(LEN=40)  :: varname1
  CHARACTER(LEN=256) :: tmpstr
  CHARACTER(LEN=256) :: vfnam
  CHARACTER(LEN=256) :: timestring
  INTEGER, PARAMETER :: itagbas = 4100
  INTEGER            :: ierr, nunit
  INTEGER            :: lfnkey
  INTEGER            :: ltimestring, lvfnam
  INTEGER            :: idat,iloc,jloc,k,jproc,itag
  INTEGER            :: ncprc,nsize,kntcol

  CHARACTER(LEN=3), PARAMETER :: fmtstr(7) =                            &
                          (/'bin','xxx','hdf','xxx','xxx','xxx','net'/)

  INTEGER               :: nxlg, nylg

  REAL,    ALLOCATABLE  :: varout(:,:,:)
  REAL,    ALLOCATABLE  :: colvar_tmp(:,:)
  INTEGER, ALLOCATABLE  :: icolp_tmp(:)
  INTEGER, ALLOCATABLE  :: jcolp_tmp(:)

  INTEGER(KIND=selected_int_kind(4)), ALLOCATABLE  :: itmp(:,:,:)
  REAL,                               ALLOCATABLE  :: hmax(:), hmin(:)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  istatus = 0
  vfnam = ' '

  CALL gtlfnkey( runname, lfnkey )

  tmpstr = dirname
  IF ( LEN_TRIM(tmpstr) == 0 .OR. tmpstr == ' ' ) THEN
    tmpstr = '.'
  END IF

  print *, ' wrtcol2grd: varid: ',varid

  WRITE (vfnam,'(4A,A3,A6)') trim(tmpstr),'/',runname(1:lfnkey),'.', &
                             fmtstr(foutfmt),varid

  CALL cvttsnd( time, timestring, ltimestring )
  lvfnam = len_trim(vfnam)

  vfnam(lvfnam+1:lvfnam+ltimestring) = timestring(1:ltimestring)
  lvfnam = lvfnam + ltimestring

!
! Merge variables if needed
!

  IF (nprocs > 1) THEN
    nxlg = (nx-3)*nproc_x + 3
    nylg = (ny-3)*nproc_y + 3
  ELSE
   nxlg = nx
   nylg = ny
  END IF

  CALL fnversn(vfnam,lvfnam)

  ALLOCATE(varout(nxlg,nylg,nz), STAT= istatus)
  CALL check_alloc_status(istatus,'wrtcol2grd:varout')
  varout=varmiss

  IF (nprocs > 1 ) THEN
    IF( myproc > 0 ) THEN ! non-zero processors send data

      itag=itagbas+1
      CALL mpsendi(ncolp,1,0,itag,ierr)

      itag=itagbas+2
      CALL mpsendi(icolp,ncolp,0,itag,ierr)
      print *, ' myproc:',myproc,' send icolp_tmp err:',ierr

      itag=itagbas+3
      CALL mpsendi(jcolp,ncolp,0,itag,ierr)
      print *, ' myproc:',myproc,' send jcolp_tmp err:',ierr
      print *, ' myproc:',myproc,' sent icolp(1),jcolp(1):',icolp(1),jcolp(1)
      print *, ' myproc:',myproc,' sent icolp(n),jcolp(n):',icolp(ncolp),jcolp(ncolp)

      nsize=nz*ncolp
      itag=itagbas+4
      CALL mpsendr(colvar,nsize,0,itag,ierr)

    ELSE   ! processor zero

      ALLOCATE(icolp_tmp(nxny),stat=istatus)
      CALL check_alloc_status(istatus,'wrtcol2grd:icolp_tmp')
      ALLOCATE(jcolp_tmp(nxny),stat=istatus)
      CALL check_alloc_status(istatus,'wrtcol2grd:jcolp_tmp')
      ALLOCATE(colvar_tmp(nz,nxny),stat=istatus)
      CALL check_alloc_status(istatus,'wrtcol2grd:colvar_tmp')


      DO idat=1,ncolp
        iloc=icolp(idat)
        jloc=jcolp(idat)
        DO k=1,nz
          varout(iloc,jloc,k)=colvar(k,idat)
        END DO
      END DO
      kntcol=ncolp

      print *, ' Processor zero xferred ',kntcol,' columns.'
!
!   Loop through all other processors
!
      DO jproc=1,(nprocs-1)

        icolp_tmp=-199
        jcolp_tmp=-199
        colvar_tmp=varmiss

        WRITE(6,'(a,i6)')  &
          ' wrtcol2grd:Gathering column data from:',jproc
        itag=itagbas+1
        CALL mprecvi(ncprc,1,jproc,itag,ierr)
        WRITE(6,'(a,i12,a,i6)') ' Recvd ncolp:',ncprc,' from:',jproc
        IF( ncprc > 0 ) THEN

          itag=itagbas+2
          CALL mprecvi(icolp_tmp,ncprc,jproc,itag,ierr)
          print *, 'from jproc:',jproc,' recv icolp_tmp err:',ierr

          itag=itagbas+3
          CALL mprecvi(jcolp_tmp,ncprc,jproc,itag,ierr)
          print *, 'from jproc:',jproc,' recv icolp_tmp err:',ierr
      print *, ' jproc:',jproc,' recv icolp(1),jcolp(1):',icolp_tmp(1),jcolp_tmp(1)
      print *, ' jproc:',jproc,' recv icolp(n),jcolp(n):',icolp_tmp(ncprc),jcolp_tmp(ncprc)

          nsize=nz*ncprc
          itag=itagbas+4
          CALL mprecvr(colvar_tmp,nsize,jproc,itag,ierr)

          DO idat=1,ncprc
            iloc=icolp_tmp(idat)
            jloc=jcolp_tmp(idat)
            DO k=1,nz
              varout(iloc,jloc,k)=colvar_tmp(k,idat)
            END DO
          END DO
          kntcol=kntcol+ncprc

        END IF  ! non-zero ncprc
      END DO  ! other processors loop

      WRITE(6,'(a,i12,a)') ' Gathered',kntcol,' non-missing columns'

      DEALLOCATE(colvar_tmp,icolp_tmp,jcolp_tmp)

    END IF

  ELSE

    DO idat=1,ncolp
      iloc=icolp(idat)
      jloc=jcolp(idat)
      DO k=1,nz
        varout(iloc,jloc,k)=colvar(k,idat)
      END DO
    END DO

  END IF

  !
  ! Write variable in file
  !
  IF (myproc == 0) PRINT *, 'WRTVAR2COL: Writing array to ', trim(vfnam)

  varname1 = varname
  varunits1= varunits

  IF (myproc == 0) THEN
    IF (foutfmt == 1) THEN
      CALL getunit( nunit)

      CALL asnctl ('NEWLOCAL', 1, ierr)
      CALL asnfile(vfnam, '-F f77 -N ieee', ierr)

      OPEN(UNIT=nunit,FILE=vfnam,STATUS='unknown',FORM='unformatted',       &
                      ERR=9000, IOSTAT=istatus)

      WRITE(nunit, ERR=9000, IOSTAT=istatus) nxlg,nylg,nz
      WRITE(nunit, ERR=9000, IOSTAT=istatus) varout

      WRITE(nunit, ERR=9000, IOSTAT=istatus) varname1
      WRITE(nunit, ERR=9000, IOSTAT=istatus) varunits1

      CLOSE(UNIT=nunit)
      CALL retunit(nunit)

      GO TO 8000

      ! I/O error handling
      ! Note that IOSTAT < 0 should not occur in this subroutine.

      9000 CONTINUE        ! I/O errors

      CLOSE(UNIT=nunit)
      CALL retunit(nunit)

      IF (istatus < 0) THEN
        PRINT *, 'WRTVAR1: I/O ERRORS OCCURRED ',                         &
                 '(possible end of record or file): ', istatus
      ELSE IF (istatus > 0) THEN
        PRINT *, 'WRTVAR1: UNRECOVERABLE I/O ERRORS OCCURRED: ',istatus
      END IF
      istatus = 2

      8000 CONTINUE

    ELSE IF (foutfmt == 3) THEN

      IF (hdfcompr > 3) THEN
        ALLOCATE (itmp(nxlg,nylg,nz),stat=istatus)
        CALL check_alloc_status(istatus, 'wrtvar2:itmp')
        ALLOCATE (hmax(nz),stat=istatus)
        CALL check_alloc_status(istatus, 'wrtvar2:hmax')
        ALLOCATE (hmin(nz),stat=istatus)
        CALL check_alloc_status(istatus, 'wrtvar2:hmin')
      END IF

      CALL hdfopen(vfnam,2,nunit)
      CALL hdfwrti(nunit, 'nx', nxlg, istatus)
      CALL hdfwrti(nunit, 'ny', nylg, istatus)
      CALL hdfwrti(nunit, 'nz', nz,   istatus)
      CALL hdfwrt3d(varout,nxlg,nylg,nz,nunit,0,hdfcompr,               &
                    varid,varname1,varunits1,itmp,hmax,hmin)
      CALL hdfclose(nunit,istatus)

      IF (hdfcompr > 3) THEN
        DEALLOCATE(itmp)
        DEALLOCATE(hmax,hmin)
      END IF

    ELSE IF (foutfmt == 7) THEN

      CALL netopen(vfnam,'C',nunit)
      CALL net_define_onevar(nunit,nxlg,nylg,nz,varid,varname1,varunits1,istatus)
      CALL netwrt3d(nunit,0,0,varid,varout,nxlg,nylg,nz)
      CALL netclose(nunit)

    ELSE
      WRITE(6,'(1x,a,I2,a,/,a)') 'Unsupported file format - ',foutfmt,    &
               ' inside wrtvar2.','Please check the subroutine arguments.'
      istatus = 1
    END IF
  END IF

  CALL mpupdatei(istatus, 1)

  DEALLOCATE(varout)

  RETURN        ! Normal return
END SUBROUTINE wrtcol2grd
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE READVAR2                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE readvar2(nx,ny,nz,var, varid,varname,varunits,               &
                    time,runname,dirname, finfmt, mpiflag, istatus)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!    To read in array 'var' from a file.
!
!  NOTE:
!    1. Two new parameters were added, finfmt and mpiflag.
!      finfmt = 1,  input file is in binary format
!      finfmt = 3,  HDF 4 file
!      finfmt = 7,  NetCDF file
!
!      mpiflag = 1, Read in the array and them split it if
!                   It is in MPI mode.
!      mpiflag = 0, Do not split.
!
!    2. Need to link with the following objects:
!
!        alloclib.o
!        mpisubs.o          or        nompsubs.o
!        netio3d.o                    nonetio3d.o
!        hdfio3d.o                    nohdfio3d.o
!
!    3. Includes 'mp.inc', should set at least the following variables
!       before call this subroutine, even for no MPI mode
!
!         mp_opt, myproc,  nproc_x, nproc_y, loc_x, loc_y
!
!  AUTHOR:
!    Yunheng Wang (2005/06/14)
!    Based on readvar1.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!  INPUT:
!
!    nx, ny, nz   Dimensions of input array 'var'.
!    varid        String of length 6 (padded with _ as necessary) to
!                 designate the input array. (E.g., 'w_____')
!    time         The model time in seconds.
!    runname      Run name
!    dirname      Input directory (use '.' for current directory)
!    finfmt       Input file format
!    mpiflag      Split or not? (=0, no split, = 1 split)
!
!  OUTPUT:
!
!    var          Array to be written out.
!    varname      String describing the field (e.g., 'Vertical Velocity')
!    varunits     String describing the units (e.g., 'm/s')
!    istatus       Exit status (0=okay, 1=warning, 2=error)
!
!  I/O:
!
!  An file named  dirname/runname.xxxvarid{time}  is read:
!             nx, ny, nz
!             var
!             varname
!             varunits
!
!  where 'xxx' is bin, hdf, or net
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER,          INTENT(IN)  :: nx,ny,nz
  CHARACTER(LEN=6), INTENT(IN)  :: varid
  CHARACTER(LEN=*), INTENT(IN)  :: runname, dirname
  REAL,             INTENT(IN)  :: time
  INTEGER,          INTENT(IN)  :: finfmt, mpiflag
  CHARACTER(LEN=*), INTENT(OUT) :: varname, varunits
  REAL,             INTENT(OUT) :: var(nx,ny,nz)

  INCLUDE  'mp.inc'

  CHARACTER(LEN=40)  :: varname1, varunits1
  CHARACTER(LEN=256) :: tmpdirnam
  CHARACTER(LEN=256) :: vfnam
  CHARACTER(LEN=256) :: timestring
  INTEGER            :: ierr, istatus, nunit

  INTEGER            :: nx_in,ny_in,nz_in, lfnkey
  INTEGER            :: lvfnam, ltimestring

  CHARACTER(LEN=3), PARAMETER :: fmtstr(7) =                            &
                          (/'bin','xxx','hdf','xxx','xxx','xxx','net'/)

  INTEGER                     :: nxlg, nylg
  REAL, ALLOCATABLE           :: varin(:,:,:)

  INTEGER(KIND=selected_int_kind(4)), ALLOCATABLE  :: itmp(:,:,:)
  REAL,                               ALLOCATABLE  :: hmax(:), hmin(:)

  INTEGER            :: i, j, k
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  istatus = 0

  CALL gtlfnkey( runname, lfnkey )

  tmpdirnam = dirname
  IF ( LEN(tmpdirnam) == 0 .OR. tmpdirnam == ' ' ) tmpdirnam = '.'

  WRITE (vfnam,'(4A,A3,A6)') trim(tmpdirnam),'/',runname(1:lfnkey),'.', &
                             fmtstr(finfmt),varid

  CALL cvttsnd( time, timestring, ltimestring )
  lvfnam = len_trim(vfnam)
  vfnam(lvfnam+1:lvfnam+ltimestring) = timestring(1:ltimestring)
  IF (mp_opt > 0 .AND. mpiflag == 0) THEN   ! MPI mode, all processes read.
     tmpdirnam = vfnam
     CALL gtsplitfn(tmpdirnam,1,1,loc_x,loc_y,1,1,0,0,1,2,vfnam,istatus)
     lvfnam = LEN_TRIM(vfnam)
     nxlg = nx
     nylg = ny
  ELSE                    ! no MPI mode or MPI mode only root process read.
     nxlg = (nx-3)*nproc_x + 3
     nylg = (ny-3)*nproc_y + 3
  END IF

  ALLOCATE(varin(nxlg,nylg,nz), STAT=istatus)
  CALL check_alloc_status(istatus, 'readvar2:varin')

  varunits  = ' '
  varunits1 = ' '
  varname1  = ' '
  IF (myproc == 0 .OR. mpiflag == 0) THEN  ! no MPI mode, myproc == 0
                                           ! MPI mode,  mpiflag = 1 --- Only root process reads
                                           !            mpiflag = 0 --- all processes read
    WRITE(6,'(/,1x,2a)') 'READVAR2: Reading array from ',trim(vfnam)

    IF (finfmt == 1) THEN
      CALL getunit( nunit)

      CALL asnctl ('NEWLOCAL', 1, ierr)
      CALL asnfile(vfnam, '-F f77 -N ieee', ierr)

      OPEN(UNIT=nunit,FILE=trim(vfnam),STATUS='old', FORM='unformatted',&
                      ERR=9000, IOSTAT=istatus)

      READ(nunit, ERR=9000, END=9000, IOSTAT=istatus) nx_in,ny_in,nz_in

      IF(nx_in /= nxlg .OR. ny_in /= nylg .OR. nz_in /= nz) THEN
        WRITE(6,'(a,/a,a,/a,3I5,/a,3I5)')                               &
            'Warning in subroutine READVAR2: Dimensions of data file ', &
            vfnam,' do not agree with the expected dimensions.',        &
            'nx, ny and nz in the data are ',nx_in,ny_in,nz_in,         &
            'nx, ny and nz expected    are ',nxlg,nylg,nz
        istatus = 2
        GO TO 9000
      END IF

      READ(nunit, ERR=9000,END=9000,IOSTAT=istatus) varin
      READ(nunit, ERR=9000,END=9000,IOSTAT=istatus) varname1
      READ(nunit, ERR=9000,END=9000,IOSTAT=istatus) varunits1

      CLOSE(UNIT=nunit)
      CALL retunit(nunit)

      GO TO 8000

      ! I/O error handling
      ! Note that warning (istatus=1) (e.g., end of record or file) is implemented as
      ! error (status=2) because of ambiguities in IOSTAT.

      9000 CONTINUE        ! I/O errors

      CLOSE(UNIT=nunit)
      CALL retunit(nunit)

      IF (istatus < 0) THEN
        PRINT *, 'READVAR2: I/O ERRORS OCCURRED ',                          &
                 '(possible end of record or file): ', istatus
      ELSE IF (istatus > 0) THEN
        PRINT *, 'READVAR2: UNRECOVERABLE I/O ERRORS OCCURRED: ', istatus
      END IF

      8000 CONTINUE

    ELSE IF (finfmt == 3) THEN

      ALLOCATE (itmp(nxlg,nylg,nz),stat=istatus)
      CALL check_alloc_status(istatus, "READVAR2:itmp")
      ALLOCATE (hmax(nz),stat=istatus)
      CALL check_alloc_status(istatus, "READVAR2:hmax")
      ALLOCATE (hmin(nz),stat=istatus)
      CALL check_alloc_status(istatus, "READVAR2:hmin")

      CALL hdfopen(vfnam,1,nunit)
      CALL hdfrd3d(nunit,varid,nxlg,nylg,nz,varin,istatus,itmp,hmax,hmin)
      CALL get_var_attr_from_hdf(nunit,varid,'comment',varname1, 40,istatus)
      CALL get_var_attr_from_hdf(nunit,varid,'units',  varunits1,40,istatus)
      CALL hdfclose(nunit,istatus)
      !print *, varname1, varunits1, nxlg, nylg,nz

      DEALLOCATE(itmp)
      DEALLOCATE(hmax,hmin)

    ELSE IF (finfmt == 7) THEN

      CALL netopen(vfnam,'R',nunit)
      CALL net_get_onevar(nunit,nx_in,ny_in,nz_in,varid,varname1,varunits1,istatus)
      CALL netread3d(nunit,0,0,varid,nx_in,ny_in,nz_in,varin)
      CALL netclose(nunit)

    ELSE
      WRITE(6,'(a,I2)') 'Unknown input file format - ',finfmt
      istatus = 1
    END IF

  END IF
  CALL mpupdatei(istatus,1)

  IF (istatus == 0) THEN                     ! successful read
    IF (mp_opt > 0 .AND. mpiflag > 0) THEN   ! split data
      CALL mpisplit3d(varin,nx,ny,nz,var)
    ELSE             ! no MPI mode, or all processes read its own data
      DO k = 1,nz
        DO j = 1,ny
          DO i = 1,nx
            var(i,j,k) = varin(i,j,k)
          END DO
        END DO
      END DO
    END IF
    varname = varname1
    varunits = varunits1
  END IF

  DEALLOCATE(varin)

  RETURN
END SUBROUTINE readvar2

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRTVAR3                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE wrtvar3(nx,ny,nz,var,varid,varname,varunits,time,runname,    &
                   dirname,postfix,foutfmt,hdfcompr,mpiflag,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!    To write an array 'var' into a file. The file can be in three format
!    based on the argument, foutfmt:
!       foutfmt = 1,      Binary file
!       foutfmt = 3,      HDF4 file,
!         hdfcompr = 0-7, Specify HDF4 compression option, see arps.input
!       foutfmt = 7,      NetCDF file
!       others,           Error message
!
!  NOTE:
!    This subroutine includes "mp.inc", So all variables in mp.inc must be
!    set before calling this subroutine. It also needs to link with the
!    following objects:
!
!           alloclib.o
!           mpisubs.o          or       nompsubs.o
!           hdfio3d.o                   nohdfio3d.o
!           netio3d.o                   nonetio3d.o
!
!  AUTHOR: Yunheng Wang (2005/06/14)
!  Based on wrtvar1. The interface is the same as wrtvar1 except that
!  two more arguments, foutfmt, hdfcompr were added.
!
!  MODIFICATION HISTORY:
!
!  INPUT:
!
!    nx, ny, nz   Dimensions of input array 'var'.
!    var          Array to be written out.
!    varid        String of length 6 (padded with _ as necessary) to
!                 designate the input array. (E.g., 'w_____')
!    varname      String describing the field (e.g., 'Vertical Velocity')
!    varunits     String describing the units (e.g., 'm/s')
!    time         The model time in seconds.
!    runname      Run name
!    dirname      Output directory (use '.' for current directory)
!    foutfmt      File output format, 1,3, or 7
!    hdfcompr     HDF file compression option, only used when foutfmt=3
!    mpiflag      Flag to indicate whether to merge the array or not
!                 = 0, Do not merge
!                 = 1, merge if the program is in MPI mode.
!
!  OUTPUT:
!
!    istatus       Exit status (0=okay, 1=unknown file format, 2=error)
!
!  I/O:
!
!  An file named  dirname/runname.xxxvarid{time}  is created:
!             nx, ny, nz
!             var
!             varname
!             varunits
!
!  where xxx is
!          bin        for foutfmt = 1
!          hdf        for foutfmt = 3
!          net        for foutfmt = 7
!
!-----------------------------------------------------------------------
!

  IMPLICIT NONE

  INTEGER,          INTENT(IN)  :: nx,ny,nz
  REAL,             INTENT(IN)  :: var(nx,ny,nz)
  CHARACTER(LEN=6), INTENT(IN)  :: varid
  CHARACTER(LEN=*), INTENT(IN)  :: postfix
  CHARACTER(LEN=*), INTENT(IN)  :: varunits
  CHARACTER(LEN=*), INTENT(IN)  :: varname
  CHARACTER(LEN=*), INTENT(IN)  :: dirname
  CHARACTER(LEN=*), INTENT(IN)  :: runname
  INTEGER,          INTENT(IN)  :: foutfmt
  INTEGER,          INTENT(IN)  :: hdfcompr
  INTEGER,          INTENT(IN)  :: mpiflag
  REAL,             INTENT(IN)  :: time
  INTEGER,          INTENT(OUT) :: istatus

  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  CHARACTER(LEN=40)  :: varunits1
  CHARACTER(LEN=40)  :: varname1
  CHARACTER(LEN=256) :: tmpdirnam
  CHARACTER(LEN=256) :: vfnam
  CHARACTER(LEN=256) :: timestring
  INTEGER            :: ierr, nunit
  INTEGER            :: lfnkey
  INTEGER            :: ltimestring, lvfnam

  CHARACTER(LEN=3), PARAMETER :: fmtstr(7) =                            &
                          (/'bin','xxx','hdf','xxx','xxx','xxx','net'/)

  INTEGER               :: nxlg, nylg
  REAL,    ALLOCATABLE  :: varout(:,:,:)

  INTEGER(KIND=selected_int_kind(4)), ALLOCATABLE  :: itmp(:,:,:)
  REAL,                               ALLOCATABLE  :: hmax(:), hmin(:)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  istatus = 0

  CALL gtlfnkey( runname, lfnkey )

  tmpdirnam = dirname
  IF ( LEN_TRIM(tmpdirnam) == 0 .OR. tmpdirnam == ' ' ) THEN
    tmpdirnam = '.'
  END IF

  WRITE (vfnam,'(4A,A3,A6)') trim(tmpdirnam),'/',runname(1:lfnkey),'.', &
                             fmtstr(foutfmt),varid

  CALL cvttsnd( time, timestring, ltimestring )
  lvfnam = len_trim(vfnam)

  vfnam(lvfnam+1:lvfnam+ltimestring) = timestring(1:ltimestring)
  vfnam=trim(vfnam)//trim(postfix)

!
! Merge variables if needed
!

  IF (mp_opt > 0 .AND. mpiflag > 0) THEN
    nxlg = (nx-3)*nproc_x + 3
    nylg = (ny-3)*nproc_y + 3
  ELSE
   IF (mp_opt > 0) THEN
     tmpdirnam = vfnam
     CALL gtsplitfn(tmpdirnam,1,1,loc_x,loc_y,1,1,0,0,0,2,vfnam,istatus)
     lvfnam = LEN_TRIM(vfnam)
   END IF
   nxlg = nx
   nylg = ny
  END IF

  ALLOCATE(varout(nxlg,nylg,nz), STAT= istatus)

  IF (mp_opt > 0 .AND. mpiflag > 0) THEN
    CALL mpimerge3d(var,nx,ny,nz,varout)
  ELSE
    varout(:,:,:) = var(:,:,:)
  END IF

!
! Write variable in file
!
  PRINT *, 'WRTVAR3: Writing array to ', trim(vfnam)

  varname1 = varname
  varunits1= varunits

  IF (myproc == 0 .OR. mpiflag == 0) THEN
    IF (foutfmt == 1) THEN
      CALL getunit( nunit)

      CALL asnctl ('NEWLOCAL', 1, ierr)
      CALL asnfile(vfnam, '-F f77 -N ieee', ierr)

      OPEN(UNIT=nunit,FILE=vfnam,STATUS='unknown',FORM='unformatted',       &
                      ERR=9000, IOSTAT=istatus)

      WRITE(nunit, ERR=9000, IOSTAT=istatus) nxlg,nylg,nz
      WRITE(nunit, ERR=9000, IOSTAT=istatus) varout

      WRITE(nunit, ERR=9000, IOSTAT=istatus) varname1
      WRITE(nunit, ERR=9000, IOSTAT=istatus) varunits1

      CLOSE(UNIT=nunit)
      CALL retunit(nunit)

      GO TO 8000

      ! I/O error handling
      ! Note that IOSTAT < 0 should not occur in this subroutine.

      9000 CONTINUE        ! I/O errors

      CLOSE(UNIT=nunit)
      CALL retunit(nunit)

      IF (istatus < 0) THEN
        PRINT *, 'WRTVAR1: I/O ERRORS OCCURRED ',                         &
                 '(possible end of record or file): ', istatus
      ELSE IF (istatus > 0) THEN
        PRINT *, 'WRTVAR1: UNRECOVERABLE I/O ERRORS OCCURRED: ',istatus
      END IF
      istatus = 2

      8000 CONTINUE

    ELSE IF (foutfmt == 3) THEN

      IF (hdfcompr > 3) THEN
        ALLOCATE (itmp(nxlg,nylg,nz),stat=istatus)
        CALL check_alloc_status(istatus, 'wrtvar3:itmp')
        ALLOCATE (hmax(nz),stat=istatus)
        CALL check_alloc_status(istatus, 'wrtvar3:hmax')
        ALLOCATE (hmin(nz),stat=istatus)
        CALL check_alloc_status(istatus, 'wrtvar3:hmin')
      END IF

      CALL hdfopen(vfnam,2,nunit)
      CALL hdfwrti(nunit, 'nx', nxlg, istatus)
      CALL hdfwrti(nunit, 'ny', nylg, istatus)
      CALL hdfwrti(nunit, 'nz', nz,   istatus)
      CALL hdfwrt3d(varout,nxlg,nylg,nz,nunit,0,hdfcompr,               &
                    varid,varname1,varunits1,itmp,hmax,hmin)
      CALL hdfclose(nunit,istatus)

      IF (hdfcompr > 3) THEN
        DEALLOCATE(itmp)
        DEALLOCATE(hmax,hmin)
      END IF

    ELSE IF (foutfmt == 7) THEN

      CALL netopen(vfnam,'C',nunit)
      CALL net_define_onevar(nunit,nxlg,nylg,nz,varid,varname1,varunits1,istatus)
      CALL netwrt3d(nunit,0,0,varid,varout,nxlg,nylg,nz)
      CALL netclose(nunit)

    ELSE
      WRITE(6,'(1x,a,I2,a,/,a)') 'Unsupported file format - ',foutfmt,    &
               ' inside wrtvar3.','Please check the subroutine arguments.'
      istatus = 1
    END IF
  END IF

  CALL mpupdatei(istatus, 1)

  DEALLOCATE(varout)

  RETURN        ! Normal return
END SUBROUTINE wrtvar3
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE READVAR3                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE readvar3(nx,ny,nz,var,varid,varname,varunits,time,runname,  &
                    dirname,postfix,finfmt,mpiflag,istatus)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!    To read in array 'var' from a file.
!
!  NOTE:
!    1. Two new parameters were added, finfmt and mpiflag.
!      finfmt = 1,  input file is in binary format
!      finfmt = 3,  HDF 4 file
!      finfmt = 7,  NetCDF file
!
!      mpiflag = 1, Read in the array and them split it if
!                   It is in MPI mode.
!      mpiflag = 0, Do not split.
!
!    2. Need to link with the following objects:
!
!        alloclib.o
!        mpisubs.o          or        nompsubs.o
!        netio3d.o                    nonetio3d.o
!        hdfio3d.o                    nohdfio3d.o
!
!    3. Includes 'mp.inc', should set at least the following variables
!       before call this subroutine, even for no MPI mode
!
!         mp_opt, myproc,  nproc_x, nproc_y, loc_x, loc_y
!
!  AUTHOR:
!    Yunheng Wang (2005/06/14)
!    Based on readvar1.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!  INPUT:
!
!    nx, ny, nz   Dimensions of input array 'var'.
!    varid        String of length 6 (padded with _ as necessary) to
!                 designate the input array. (E.g., 'w_____')
!    time         The model time in seconds.
!    runname      Run name
!    dirname      Input directory (use '.' for current directory)
!    finfmt       Input file format
!    mpiflag      Split or not? (=0, no split, = 1 split)
!
!  OUTPUT:
!
!    var          Array to be written out.
!    varname      String describing the field (e.g., 'Vertical Velocity')
!    varunits     String describing the units (e.g., 'm/s')
!    istatus       Exit status (0=okay, 1=warning, 2=error)
!
!  I/O:
!
!  An file named  dirname/runname.xxxvarid{time}  is read:
!             nx, ny, nz
!             var
!             varname
!             varunits
!
!  where 'xxx' is bin, hdf, or net
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER,          INTENT(IN)  :: nx,ny,nz
  CHARACTER(LEN=6), INTENT(IN)  :: varid
  CHARACTER(LEN=*), INTENT(IN)  :: postfix
  CHARACTER(LEN=*), INTENT(IN)  :: runname, dirname
  REAL,             INTENT(IN)  :: time
  INTEGER,          INTENT(IN)  :: finfmt, mpiflag
  CHARACTER(LEN=*), INTENT(OUT) :: varname, varunits
  REAL,             INTENT(OUT) :: var(nx,ny,nz)

  INCLUDE  'mp.inc'

  CHARACTER(LEN=40)  :: varname1, varunits1
  CHARACTER(LEN=256) :: tmpdirnam
  CHARACTER(LEN=256) :: vfnam
  CHARACTER(LEN=256) :: timestring
  INTEGER            :: ierr, istatus, nunit

  INTEGER            :: nx_in,ny_in,nz_in, lfnkey
  INTEGER            :: lvfnam, ltimestring

  CHARACTER(LEN=3), PARAMETER :: fmtstr(7) =                            &
                          (/'bin','xxx','hdf','xxx','xxx','xxx','net'/)

  INTEGER                     :: nxlg, nylg
  REAL, ALLOCATABLE           :: varin(:,:,:)

  INTEGER(KIND=selected_int_kind(4)), ALLOCATABLE  :: itmp(:,:,:)
  REAL,                               ALLOCATABLE  :: hmax(:), hmin(:)

  INTEGER            :: i, j, k
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  istatus = 0

  CALL gtlfnkey( runname, lfnkey )

  tmpdirnam = dirname
  IF ( LEN(tmpdirnam) == 0 .OR. tmpdirnam == ' ' ) tmpdirnam = '.'

  WRITE (vfnam,'(4A,A3,A6)') trim(tmpdirnam),'/',runname(1:lfnkey),'.', &
                             fmtstr(finfmt),varid

  CALL cvttsnd( time, timestring, ltimestring )
  vfnam=trim(vfnam)//trim(postfix)
  lvfnam = len_trim(vfnam)
  vfnam(lvfnam+1:lvfnam+ltimestring) = timestring(1:ltimestring)

  IF (mp_opt >0 .AND. mpiflag == 0) THEN   ! MPI mode, all processes read.
     tmpdirnam = vfnam
     CALL gtsplitfn(tmpdirnam,1,1,loc_x,loc_y,1,1,0,0,1,2,vfnam,istatus)
     lvfnam = LEN_TRIM(vfnam)
     nxlg = nx
     nylg = ny
  ELSE                                     ! no MPI mode or MPI mode only root process read
     nxlg = (nx-3)*nproc_x + 3
     nylg = (ny-3)*nproc_y + 3
  END IF

  ALLOCATE(varin(nxlg,nylg,nz), STAT=istatus)
  CALL check_alloc_status(istatus, 'readvar3:varin')

  PRINT *, 'READVAR3: Reading array from ',trim(vfnam)

  IF (myproc == 0 .OR. mpiflag == 0) THEN  ! no MPI mode, myproc == 0
                                           ! MPI mode,  mpiflag = 1 --- Only root process reads
                                           !            mpiflag = 0 --- all processes read
    IF (finfmt == 1) THEN
      CALL getunit( nunit)

      CALL asnctl ('NEWLOCAL', 1, ierr)
      CALL asnfile(vfnam, '-F f77 -N ieee', ierr)

      OPEN(UNIT=nunit,FILE=trim(vfnam),STATUS='old', FORM='unformatted',&
                      ERR=9000, IOSTAT=istatus)

      READ(nunit, ERR=9000, END=9000, IOSTAT=istatus) nx_in,ny_in,nz_in

      IF(nx_in /= nxlg .OR. ny_in /= nylg .OR. nz_in /= nz) THEN
        WRITE(6,'(a,/a,a,/a,3I5,/a,3I5)')                                   &
            'Warning in subroutine READVAR3: Dimensions of data file ',     &
            vfnam,' do not agree with the expected dimensions.',            &
            'nx, ny and nz in the data are ',nx_in,ny_in,nz_in,             &
            'nx, ny and nz expected    are ',nxlg,nylg,nz
        istatus = 2
        GO TO 9000
      END IF

      READ(nunit, ERR=9000,END=9000,IOSTAT=istatus) varin
      READ(nunit, ERR=9000,END=9000,IOSTAT=istatus) varname1
      READ(nunit, ERR=9000,END=9000,IOSTAT=istatus) varunits1

      CLOSE(UNIT=nunit)
      CALL retunit(nunit)

      GO TO 8000

      ! I/O error handling
      ! Note that warning (istatus=1) (e.g., end of record or file) is implemented as
      ! error (status=2) because of ambiguities in IOSTAT.

      9000 CONTINUE        ! I/O errors

      CLOSE(UNIT=nunit)
      CALL retunit(nunit)

      IF (istatus < 0) THEN
        PRINT *, 'READVAR3: I/O ERRORS OCCURRED ',                          &
                 '(possible end of record or file): ', istatus
      ELSE IF (istatus > 0) THEN
        PRINT *, 'READVAR3: UNRECOVERABLE I/O ERRORS OCCURRED: ', istatus
      END IF

      8000 CONTINUE

    ELSE IF (finfmt == 3) THEN

      ALLOCATE (itmp(nxlg,nylg,nz),stat=istatus)
      CALL check_alloc_status(istatus, "READVAR3:itmp")
      ALLOCATE (hmax(nz),stat=istatus)
      CALL check_alloc_status(istatus, "READVAR3:hmax")
      ALLOCATE (hmin(nz),stat=istatus)
      CALL check_alloc_status(istatus, "READVAR3:hmin")

      CALL hdfopen(vfnam,1,nunit)
      CALL hdfrd3d(nunit,varid,nxlg,nylg,nz,varin,istatus,itmp,hmax,hmin)
      CALL get_var_attr_from_hdf(nunit,varid,'comment',varname1, 40,istatus)
      CALL get_var_attr_from_hdf(nunit,varid,'units',  varunits1,40,istatus)
      CALL hdfclose(nunit,istatus)
      print *, varname1, varunits1, nxlg, nylg,nz

      DEALLOCATE(itmp)
      DEALLOCATE(hmax,hmin)

    ELSE IF (finfmt == 7) THEN

      CALL netopen(vfnam,'R',nunit)
      CALL net_get_onevar(nunit,nx_in,ny_in,nz_in,varid,varname1,varunits1,istatus)
      CALL netread3d(nunit,0,0,varid,nxlg,nylg,nz,varin)
      CALL netclose(nunit)

    ELSE
      WRITE(6,'(a,I2)') 'Unknown input file format - ',finfmt
      istatus = 1
    END IF

  END IF

  IF (istatus == 0) THEN                     ! successful read
    IF (mp_opt > 0 .AND. mpiflag > 0) THEN   ! split data
      CALL mpisplit3d(varin,nx,ny,nz,var)
    ELSE                                     ! no MPI mode, or all processes read its own data
      DO k = 1,nz
        DO j = 1,ny
          DO i = 1,nx
            var(i,j,k) = varin(i,j,k)
          END DO
        END DO
      END DO
    END IF
    varname = varname1
    varunits = varunits1
  END IF

  DEALLOCATE(varin)

  RETURN
END SUBROUTINE readvar3
!
!##################################################################
!##################################################################
!######                                                      ######
!######          SUBROUTINE GET_INPUT_FILE_NAMES             ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE get_input_file_names(unum,hinfmt,grdbasfn,hisfile,nhisfile)
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!    7/17/2000.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    hinfmt
!
!  OUTPUT:
!
!    grdbasfn,hisfile,nhisfile
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER, INTENT(IN) :: unum
  INTEGER :: hinfmt ! input data format
  INTEGER :: hdmpinopt

  CHARACTER (LEN=256) :: hdmpftrailer
  CHARACTER (LEN=256) :: hdmpfheader
  REAL :: tintv_dmpin, tbgn_dmpin, tend_dmpin

  INTEGER, PARAMETER  :: nhisfile_max = 200
  INTEGER             :: nhisfile
  CHARACTER (LEN=256) :: grdbasfn
  CHARACTER (LEN=256) :: hisfile(nhisfile_max)

  NAMELIST /history_data/ hinfmt, nhisfile, grdbasfn, hisfile,          &
           hdmpinopt,hdmpfheader,hdmpftrailer,                          &
           tintv_dmpin, tbgn_dmpin, tend_dmpin

  INTEGER :: lengbf,nf,length
!
!-----------------------------------------------------------------------
!
!  Get the names of the input data files.
!
!-----------------------------------------------------------------------
!
  READ(unum,history_data,ERR=100)

  WRITE(6,'(a)')'Namelist history_data was successfully read.'
  WRITE(6,'(2x,a,i3)') 'Input hdmpinopt=', hdmpinopt
  WRITE(6,'(2x,a,i3)') 'Input hinfmt   =', hinfmt
  WRITE(6,'(2x,a,i3)') 'Input nhisfile =', nhisfile

  IF( hdmpinopt == 1 ) THEN

    CALL gthinfns(hdmpfheader,hdmpftrailer,hinfmt,                      &
                  tintv_dmpin, tbgn_dmpin, tend_dmpin,                  &
                  grdbasfn,hisfile,nhisfile)

  END IF

  IF( nhisfile > nhisfile_max ) THEN

    WRITE(6,'(/a,a,I5,/a/)')                                            &
        'The number of history files to be processed',                  &
        ' exceeded the maximum number ',nhisfile_max,                   &
        'Please increase the size of array hisfile.'

  END IF

  lengbf =len_trim( grdbasfn )
  WRITE(6,'(/2x,a,a)')'The grid base-state file name is ',              &
                    grdbasfn(1:lengbf)

  DO nf=1,nhisfile
    WRITE(6,'(2x,a,i5,a,a)')                                            &
       'History file      No. ',nf,'   is ',trim(hisfile(nf))
  END DO
!
  GO TO 10

  100   WRITE(6,'(a,a)') 'Error reading NAMELIST file. ',               &
                       'Job stopped in GET_INPUT_FILE_NAMES.'
  CALL arpsstop('arpsstop called from GET_INPUT_FILE_NAMES error'//     &
                'reading namelist',1)

  10    CONTINUE

  RETURN
END SUBROUTINE get_input_file_names

!
!##################################################################
!##################################################################
!######                                                      ######
!######          SUBROUTINE GET_OUTPUT_DIRNAME               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE get_output_dirname(mpigroup,indir,time,mmbid,dirname,istatus)

!#######################################################################
!
! PURPOSE: Find the proper output directory name at the specified
!          forecast time
!
!#######################################################################

  IMPLICIT NONE

  INTEGER,            INTENT(IN)  :: mpigroup   ! 0: only 1 process calls
                                                ! 1: all process members call
  CHARACTER(LEN=256), INTENT(IN)  :: indir      ! got from namelist
  REAL,               INTENT(IN)  :: time       ! forecast time
  INTEGER,            INTENT(IN)  :: mmbid      ! ensemble id
  CHARACTER(LEN=256), INTENT(OUT) :: dirname    ! specific directory name at this time

  INTEGER, INTENT(OUT) :: istatus

  INCLUDE 'mp.inc'
!-----------------------------------------------------------------------

  INTEGER :: indx1, indx2, indx3
  INTEGER :: digit
  CHARACTER(LEN=1)   :: letter,trail
  CHARACTER(LEN=20)  :: fmtstr
  CHARACTER(LEN=256) :: tmpstr
  LOGICAL :: iexist

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  trail = ' '
  indx2 = LEN_TRIM(indir)
  IF (indir(indx2:indx2) /= '/') trail = '/'

  WRITE(dirname,'(2a)') TRIM(indir),trail
  indx2 = LEN_TRIM(dirname)

  tmpstr = dirname
  indx1 = INDEX(tmpstr,'%')   ! find the first '%'

  IF (mpigroup == 0) THEN    ! Require the parent directoy only

    IF (indx1 > 0) THEN
      indx2 = INDEX(indir(1:indx1),'/',.TRUE.)
      IF (indx2 <= 0) THEN
        WRITE(dirname,'(a)') './'
      ELSE
        WRITE(dirname,'(a)') indir(1:indx2)
      END IF
    END IF

  ELSE        ! Subsitute all patterns

    DO WHILE (indx1 > 0)

      digit = -1
      READ(tmpstr(indx1+1:indx1+2),'(I1,a1)') digit,letter
      IF (digit > 9 .OR. digit < 1) THEN
        WRITE(*,'(1x,a,I0,2a1,3a)') 'ERROR: Unknown format (%',digit,',', &
        letter,') in the string ',TRIM(indir),' for dirname.'
        CALL arpsstop('ERROR: Wrong namelist paramter,',1)
      END IF

      IF (letter == 'T') THEN
        IF (time < 0.0 .AND. time > -1.0E-10) THEN  ! grid and base file
          WRITE(dirname,'(3a)') tmpstr(1:indx1-1),'grdbas',tmpstr(indx1+3:indx2)
          digit = 6
        ELSE
          WRITE(fmtstr,'(a,I1,a,I1,a)') '(a,I',digit,'.',digit,',a)'
          WRITE(dirname,fmt=fmtstr) tmpstr(1:indx1-1),NINT(time),tmpstr(indx1+3:indx2)
        END IF
      ELSE IF (letter == 'N') THEN
        IF (mmbid < 0) THEN         ! background mean
          WRITE(dirname,'(3a)') tmpstr(1:indx1-1),'fmean',tmpstr(indx1+3:indx2)
          digit = 5
        ELSE IF (mmbid == 0)  THEN  ! analysis mean
          WRITE(dirname,'(3a)') tmpstr(1:indx1-1),'amean',tmpstr(indx1+3:indx2)
          digit = 5
        ELSE
          WRITE(fmtstr,'(a,I1,a,I1,a)') '(a,I',digit,'.',digit,',a)'
          WRITE(dirname,fmt=fmtstr) tmpstr(1:indx1-1),mmbid,tmpstr(indx1+3:indx2)
        END IF
      ELSE
        WRITE(*,'(1x,a,a1,3a)') 'ERROR: Unsupported template "',          &
        letter,'" in the string ',TRIM(indir),' for dirname.'
        CALL arpsstop('ERROR: Wrong namelist paramter.',2)
      END IF

      indx3 = INDEX(tmpstr(indx1+3:indx2),'%')   ! find next '%'
      IF (indx3 <= 0) EXIT

      indx1 = indx1-1 + digit + indx3
      tmpstr = dirname
      indx2 = LEN_TRIM(tmpstr)

    END DO

  END IF
  !PRINT *, 'dirname = ', TRIM(dirname)
!
!-----------------------------------------------------------------------
!
!  Check if the specified output directory exists, if not, create it.
!
!-----------------------------------------------------------------------

  IF( mpigroup == 0 .OR. (mpigroup == 1 .AND. myproc == 0) ) THEN
    CALL inquiredir(TRIM(dirname),iexist)

    IF( .NOT. iexist ) THEN
      !CALL unixcmd( 'mkdir -p '//dirname(1:ldirnam) )
      CALL makedir(dirname,istatus)
      WRITE(6,'(5x,a,2(/5x,a))') 'Required output directory <'          &
        //TRIM(dirname)//'> not found.', 'It was created by the program.'
    END IF

    WRITE(6,'(5x,a,F9.2,3a)')                                           &
      'Output files for forecast time = ',time,' will be in directory <',TRIM(dirname),'>.'

  END IF

  IF (mpigroup > 0)  CALL mpbarrier  ! wait until the directory is created.

  RETURN
END SUBROUTINE get_output_dirname
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE GTHINFNS                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE gthinfns(hdmpfheader,hdmpftrailer,hinfmt,                    &
                    tintv_dmpin, tbgn_dmpin, tend_dmpin,                &
                    grdbasfn,hisfile,nhisfile)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Retrun a list of history file names given the start, end times
!  and time interval, as well as the name header and trailing string.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!    4/7/2000
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    hdmpfheader,hdmpftrailer,hinfmt
!    tintv_dmpin,tbgn_dmpin,tend_dmpin
!
!  OUTPUT:
!
!    grdbasfn,hisfile,nhisfile
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: hinfmt

  CHARACTER (LEN=256) :: hdmpfheader,hdmpftrailer
  REAL :: tintv_dmpin, tbgn_dmpin, tend_dmpin

  INTEGER, PARAMETER  :: nhisfile_max = 200
  INTEGER             :: nhisfile
  CHARACTER (LEN=256) :: grdbasfn
  CHARACTER (LEN=256) :: hisfile(nhisfile_max)

  REAL :: time
  CHARACTER (LEN=80 ) :: timsnd
  INTEGER :: tmstrln
  INTEGER :: lheader,ltrailer
  INTEGER :: length,i
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  lheader = len_trim(hdmpfheader)

  IF( hinfmt == 1.OR.hinfmt == 0) THEN
    grdbasfn = hdmpfheader(1:lheader)//'.bingrdbas'
  ELSE IF( hinfmt == 2 ) THEN   ! Formatted ASCII dump
    grdbasfn = hdmpfheader(1:lheader)//'.ascgrdbas'
  ELSE IF( hinfmt == 3 ) THEN   ! HDF data dump
    grdbasfn = hdmpfheader(1:lheader)//'.hdfgrdbas'
  ELSE IF( hinfmt == 4 ) THEN   ! Packed binary dump
    grdbasfn = hdmpfheader(1:lheader)//'.pakgrdbas'
  ELSE IF( hinfmt == 5 ) THEN   ! Savi3D data dump

!-----------------------------------------------------------------------
!  For Savi3D data dump, the grid and base state information is
!  always written together with the other fields.
!-----------------------------------------------------------------------

  ELSE IF( hinfmt == 6 ) THEN   ! Binary with skipping
    grdbasfn = hdmpfheader(1:lheader)//'.bn2grdbas'
  ELSE IF( hinfmt == 7 ) THEN   ! NetCDF format
    grdbasfn = hdmpfheader(1:lheader)//'.netgrdbas'
  ELSE IF( hinfmt == 8 ) THEN   ! Packed NetCDF format
    grdbasfn = hdmpfheader(1:lheader)//'.nc'
  ELSE IF( hinfmt == 9 ) THEN   ! GrADS data dump

!-----------------------------------------------------------------------
!  For GrADS data dump, the grid and base state information is
!  always written together with the other fields.
!-----------------------------------------------------------------------

  ELSE IF( hinfmt == 10 ) THEN  ! GRIB format
    grdbasfn = hdmpfheader(1:lheader)//'.grbgrdbas'
  END IF

  ltrailer= len_trim(hdmpftrailer)

  IF( ltrailer > 0 ) THEN
    length = LEN( grdbasfn)
    CALL strlnth( grdbasfn, length )
    grdbasfn(length+1:length+ltrailer)=hdmpftrailer(1:ltrailer)
  END IF

  time = tbgn_dmpin
  nhisfile = 0
  DO i=1,100000000

    IF( time > tend_dmpin + 0.01*tintv_dmpin ) GO TO 105
    nhisfile = nhisfile + 1

    IF( nhisfile > nhisfile_max) THEN
      WRITE(6,'(1x,a,i3,/1x,a,/1x,a,/1x,a)')                            &
          'The number of history files to be processed exceeded ',      &
          nhisfile_max,                                                 &
          ' please reduce the number of files to be processed',         &
          'in a single job, or edit the program and re-set parameter',  &
          ' nhisfile_max to a larger value. '
      CALL arpsstop('arpsstop called from GTHINFNS number of files'//   &
          'exceeded',1)
    END IF

    CALL cvttsnd( time , timsnd, tmstrln )

    hisfile(i)=' '

    IF( hinfmt == 1.OR.hinfmt == 0) THEN       ! Unformatted binary dump
      hisfile(i)(1:lheader+4+tmstrln)                                   &
          =hdmpfheader(1:lheader)//'.bin'//timsnd(1:tmstrln)
    ELSE IF( hinfmt == 2 ) THEN   ! Formatted ASCII dump
      hisfile(i)(1:lheader+4+tmstrln)                                   &
          =hdmpfheader(1:lheader)//'.asc'//timsnd(1:tmstrln)
    ELSE IF( hinfmt == 3 ) THEN   ! HDF data dump
      hisfile(i)(1:lheader+4+tmstrln)                                   &
          =hdmpfheader(1:lheader)//'.hdf'//timsnd(1:tmstrln)
    ELSE IF( hinfmt == 4 ) THEN   ! Packed binary dump
      hisfile(i)(1:lheader+4+tmstrln)                                   &
          =hdmpfheader(1:lheader)//'.pak'//timsnd(1:tmstrln)
    ELSE IF( hinfmt == 5 ) THEN   ! Data dump for Savi3D
      hisfile(i)(1:lheader+4+tmstrln)                                   &
          =hdmpfheader(1:lheader)//'.svi'
    ELSE IF( hinfmt == 6 ) THEN   ! Binary with skipping
      hisfile(i)(1:lheader+4+tmstrln)                                   &
          =hdmpfheader(1:lheader)//'.bn2'//timsnd(1:tmstrln)
    ELSE IF( hinfmt == 7 ) THEN   ! NetCDF format
      hisfile(i)(1:lheader+4+tmstrln)                                   &
          =hdmpfheader(1:lheader)//'.net'//timsnd(1:tmstrln)
    ELSE IF( hinfmt == 8 ) THEN   ! Packed NetCDF format
      hisfile(i)(1:lheader+4+tmstrln)                                   &
          =hdmpfheader(1:lheader)//'.nc'
    ELSE IF( hinfmt == 9 ) THEN   ! Data dump for GrADS
      hisfile(i)(1:lheader+4+tmstrln)                                   &
          =hdmpfheader(1:lheader)//'.gad'
    ELSE IF( hinfmt == 10 ) THEN   ! Data dump for GrADS
      hisfile(i)(1:lheader+4+tmstrln)                                   &
          =hdmpfheader(1:lheader)//'.grb'//timsnd(1:tmstrln)
    ELSE IF( hinfmt == 11 ) THEN   ! Data dump for Vis5D
      hisfile(i)(1:lheader+4+tmstrln)                                   &
          =hdmpfheader(1:lheader)//'.v5d'//timsnd(1:tmstrln)
    END IF

    IF( ltrailer > 0 ) THEN
      length = len_trim( hisfile(i) )
      hisfile(i)(length+1:length+ltrailer)                              &
                =hdmpftrailer(1:ltrailer)
    END IF

    time = time + tintv_dmpin

  END DO

  105   CONTINUE

  RETURN
END SUBROUTINE gthinfns

SUBROUTINE gtsplitfn(fbase,nxpatch,nypatch,iorig,jorig,iloc,jloc,       &
                     ioffset,joffset,rdwrtflag,dbglvl,filename,istatus)

!#######################################################################
! Get split file name
!#######################################################################

  IMPLICIT NONE

  CHARACTER(LEN=*),   INTENT(IN)  :: fbase     ! Filename base
                                               ! i.e. the normal ARPS file name no in split format
  INTEGER,            INTENT(IN)  :: nxpatch, nypatch
                                   ! number of patches each process will handle
  INTEGER,            INTENT(IN)  :: iorig, jorig
                                   ! Processor rank
  INTEGER,            INTENT(IN)  :: iloc, jloc
                                   ! local rank within this process, must within (1, nxpatch) and (1, nypatch)
  INTEGER,            INTENT(IN)  :: ioffset, joffset
                                   ! Offset value to be added. Both are 0 if all patches are handled
  INTEGER,            INTENT(IN)  :: rdwrtflag
                                   ! Read write flag,
                                   ! = 1, To read, file must already exist.
                                   ! = 0, To write
  INTEGER,            INTENT(IN)  :: dbglvl

  CHARACTER(LEN=256), INTENT(OUT) :: filename
  INTEGER,            INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  INTEGER :: inum, jnum
  CHARACTER(LEN=256) :: filename_saved
  LOGICAL :: fexist
  INTEGER :: flen

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  inum = (iorig-1)*nxpatch+iloc+ioffset
  jnum = (jorig-1)*nypatch+jloc+joffset

  filename = ' '
  WRITE(filename,'(2a,2I3.3)') TRIM(fbase),'_',inum,jnum

  INQUIRE( FILE = TRIM(filename), EXIST = fexist )

  IF (rdwrtflag > 0 ) THEN      ! File must already exist

    IF (.NOT. fexist ) THEN

      IF (dbglvl > 100) WRITE(6,'(1x,3a,/,10x,a)')                      &
       'WARNING: file "',TRIM(filename),'" does not exist.',            &
       'Trying old file name convention with 4 characters for patch number ...'

      filename_saved = filename

      filename = ' '
      WRITE(filename,'(2a,2I2.2)') TRIM(fbase),'_',inum,jnum

      INQUIRE( FILE = TRIM(filename), EXIST = fexist )

      IF (.NOT. fexist) THEN
        istatus = -1
        IF (dbglvl >=0 ) THEN  ! Skip error message when dbglvl < 0.
        WRITE(6,'(1x,3a,/,8x,3a,/,8x,a,I2)')                            &
         'ERROR: Both file  "',TRIM(filename_saved),'" and file with old name', &
         'convention "',TRIM(filename),'" do not exist.', &
         'Subroutine gtsplitfn returning with istatus = ',istatus
        END IF
      END IF

    END IF

  ELSE          ! If already exists, add version number

    flen = LEN_TRIM(filename)
    CALL fnversn(filename, flen)

  END IF

  RETURN
END SUBROUTINE gtsplitfn

!#######################################################################

SUBROUTINE dbgwrt(unt,var,nx,ny,nz,ibgn,iend,jbgn,jend,kbgn,kend,istatus)

!-----------------------------------------------------------------------
!
! Purpose:
!   Write ASCII data to an opened file with unit "unt" of the value of a
!   3D array. Note that the index written will be global index instead of
!   local index passed in.
!
!----------------------------------------------------------------------
! AUTHOR: Yunheng Wang (2009)
!
!----------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: unt
  INTEGER, INTENT(IN)  :: nx,ny,nz
  REAL,    INTENT(IN)  :: var(nx,ny,nz)
  INTEGER, INTENT(IN)  :: ibgn,iend,jbgn,jend,kbgn,kend

  INTEGER, INTENT(OUT) :: istatus

  INCLUDE 'mp.inc'

  INTEGER :: i,j,k
  INTEGER :: ilg, jlg

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  DO k = kbgn,kend
    DO j = jbgn,jend
      DO i = ibgn, iend
        ilg = (loc_x-1)*(nx-3)+i
        jlg = (loc_y-1)*(ny-3)+j
        IF ( MOD((ilg-2),5) == 0 ) WRITE(unt,FMT='(/,3I4,a)',ADVANCE='NO') ilg,jlg,k,':'
        WRITE(unt,FMT='(F12.5)',ADVANCE='NO') var(i,j,k)
      END DO
      WRITE(unt,'(/,a)') '------------------------------------'
    END DO
    WRITE(unt,'(/,a)') '===================================='
  END DO

  istatus=0

  RETURN
END SUBROUTINE dbgwrt

!#######################################################################

SUBROUTINE min_set_scalars(anscalar,anscalarq,aqnames,aqdescp,           &
                       aP_QC,aP_QR,aP_QI,aP_QS,aP_QG,aP_QH,             &
                       aP_NC,aP_NR,aP_NI,aP_NS,aP_NG,aP_NH,             &
                             aP_ZR,aP_ZI,aP_ZS,aP_ZG,aP_ZH,             &
                       bnscalar,bnscalarq,bqnames,bqdescp,              &
                       bP_QC,bP_QR,bP_QI,bP_QS,bP_QG,bP_QH,             &
                       bP_NC,bP_NR,bP_NI,bP_NS,bP_NG,bP_NH,             &
                             bP_ZR,bP_ZI,bP_ZS,bP_ZG,bP_ZH,             &
                       cnscalar,cnscalarq,cqnames,cqdescp,              &
                       cP_QC,cP_QR,cP_QI,cP_QS,cP_QG,cP_QH,             &
                       cP_NC,cP_NR,cP_NI,cP_NS,cP_NG,cP_NH,             &
                             cP_ZR,cP_ZI,cP_ZS,cP_ZG,cP_ZH,             &
                       aQindex,bQindex, istatus )

!-----------------------------------------------------------------------
!
! Purpose:
!   Give two sets (prefixed with "a" and "b" respectively) of scalar indices,
!   find the minimum common set (intersection) and assign it to set "c".
!   aQindex, bQindex are also set for convention.
!
!----------------------------------------------------------------------
! AUTHOR: Yunheng Wang (08/12/2009)
!
!----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: anscalar, anscalarq
  INTEGER, INTENT(IN) :: aP_QC, aP_QR, aP_QI, aP_QS, aP_QH, aP_QG
  INTEGER, INTENT(IN) :: aP_NC, aP_NR, aP_NI, aP_NS, aP_NH, aP_NG
  INTEGER, INTENT(IN) ::        aP_ZR, aP_ZI, aP_ZS, aP_ZH, aP_ZG

  CHARACTER(LEN=40), INTENT(IN) :: aqnames(20)
  CHARACTER(LEN=40), INTENT(IN) :: aqdescp(20)

  INTEGER, INTENT(IN) :: bnscalar, bnscalarq
  INTEGER, INTENT(IN) :: bP_QC, bP_QR, bP_QI, bP_QS, bP_QH, bP_QG
  INTEGER, INTENT(IN) :: bP_NC, bP_NR, bP_NI, bP_NS, bP_NH, bP_NG
  INTEGER, INTENT(IN) ::        bP_ZR, bP_ZI, bP_ZS, bP_ZH, bP_ZG

  CHARACTER(LEN=40), INTENT(IN) :: bqnames(20)
  CHARACTER(LEN=40), INTENT(IN) :: bqdescp(20)

  INTEGER, INTENT(OUT) :: cnscalar, cnscalarq
  INTEGER, INTENT(OUT) :: cP_QC, cP_QR, cP_QI, cP_QS, cP_QH, cP_QG
  INTEGER, INTENT(OUT) :: cP_NC, cP_NR, cP_NI, cP_NS, cP_NH, cP_NG
  INTEGER, INTENT(OUT) ::        cP_ZR, cP_ZI, cP_ZS, cP_ZH, cP_ZG

  CHARACTER(LEN=40), INTENT(OUT) :: cqnames(20)
  CHARACTER(LEN=40), INTENT(OUT) :: cqdescp(20)

  INTEGER, INTENT(OUT) :: aQindex(anscalar), bQindex(anscalar)

  INTEGER, INTENT(OUT) :: istatus

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  cnscalar  = 0
  cnscalarq = 0
  cP_QC = -1; cP_QR = -1; cP_QI = -1; cP_QS = -1; cP_QH = -1; cP_QG = -1
  cP_NC = -1; cP_NR = -1; cP_NI = -1; cP_NS = -1; cP_NG = -1; cP_NH = -1
              cP_ZR = -1; cP_ZI = -1; cP_ZS = -1; cP_ZG = -1; cP_ZH = -1;
  cqnames(:)= ' '; cqdescp(:)= ' '

  IF (aP_QC > 0 .AND. bP_QC > 0) THEN
    cnscalar = cnscalar + 1
    cP_QC    = cnscalar
    cqnames(cP_QC) = aqnames(aP_QC)
    cqdescp(cP_QC) = aqdescp(aP_QC)
    aQindex(cP_QC) = aP_QC
    bQindex(cP_QC) = bP_QC
  END IF

  IF (aP_QR > 0 .AND. bP_QR > 0) THEN
    cnscalar = cnscalar + 1
    cP_QR    = cnscalar
    cqnames(cP_QR) = aqnames(aP_QR)
    cqdescp(cP_QR) = aqdescp(aP_QR)
    aQindex(cP_QR) = aP_QR
    bQindex(cP_QR) = bP_QR
  END IF

  IF (aP_QI > 0 .AND. bP_QI > 0) THEN
    cnscalar = cnscalar + 1
    cP_QI    = cnscalar
    cqnames(cP_QI) = aqnames(aP_QI)
    cqdescp(cP_QI) = aqdescp(aP_QI)
    aQindex(cP_QI) = aP_QI
    bQindex(cP_QI) = bP_QI
  END IF

  IF (aP_QS > 0 .AND. bP_QS > 0) THEN
    cnscalar = cnscalar + 1
    cP_QS    = cnscalar
    cqnames(cP_QS) = aqnames(aP_QS)
    cqdescp(cP_QS) = aqdescp(aP_QS)
    aQindex(cP_QS) = aP_QS
    bQindex(cP_QS) = bP_QS
  END IF

  IF (aP_QG > 0 .AND. bP_QG > 0) THEN
    cnscalar = cnscalar + 1
    cP_QG    = cnscalar
    cqnames(cP_QG) = aqnames(aP_QG)
    cqdescp(cP_QG) = aqdescp(aP_QG)
    aQindex(cP_QG) = aP_QG
    bQindex(cP_QG) = bP_QG
  END IF

  IF (aP_QH > 0 .AND. bP_QH > 0) THEN
    cnscalar = cnscalar + 1
    cP_QH    = cnscalar
    cqnames(cP_QH) = aqnames(aP_QH)
    cqdescp(cP_QH) = aqdescp(aP_QH)
    aQindex(cP_QH) = aP_QH
    bQindex(cP_QH) = bP_QH
  END IF

  cnscalarq = cnscalar

  IF (aP_NC > 0 .AND. bP_NC > 0) THEN
    cnscalar = cnscalar + 1
    cP_NC    = cnscalar
    cqnames(cP_NC) = aqnames(aP_NC)
    cqdescp(cP_NC) = aqdescp(aP_NC)
    aQindex(cP_NC) = aP_NC
    bQindex(cP_NC) = bP_NC
  END IF

  IF (aP_NR > 0 .AND. bP_NR > 0) THEN
    cnscalar = cnscalar + 1
    cP_NR    = cnscalar
    cqnames(cP_NR) = aqnames(aP_NR)
    cqdescp(cP_NR) = aqdescp(aP_NR)
    aQindex(cP_NR) = aP_NR
    bQindex(cP_NR) = bP_NR
  END IF

  IF (aP_NI > 0 .AND. bP_NI > 0) THEN
    cnscalar = cnscalar + 1
    cP_NI    = cnscalar
    cqnames(cP_NI) = aqnames(aP_NI)
    cqdescp(cP_NI) = aqdescp(aP_NI)
    aQindex(cP_NI) = aP_NI
    bQindex(cP_NI) = bP_NI
  END IF

  IF (aP_NS > 0 .AND. bP_NS > 0) THEN
    cnscalar = cnscalar + 1
    cP_NS    = cnscalar
    cqnames(cP_NS) = aqnames(aP_NS)
    cqdescp(cP_NS) = aqdescp(aP_NS)
    aQindex(cP_NS) = aP_NS
    bQindex(cP_NS) = bP_NS
  END IF

  IF (aP_NG > 0 .AND. bP_NG > 0) THEN
    cnscalar = cnscalar + 1
    cP_NG    = cnscalar
    cqnames(cP_NG) = aqnames(aP_NG)
    cqdescp(cP_NG) = aqdescp(aP_NG)
    aQindex(cP_NG) = aP_NG
    bQindex(cP_NG) = bP_NG
  END IF

  IF (aP_NH > 0 .AND. bP_NH > 0) THEN
    cnscalar = cnscalar + 1
    cP_NH    = cnscalar
    cqnames(cP_NH) = aqnames(aP_NH)
    cqdescp(cP_NH) = aqdescp(aP_NH)
    aQindex(cP_NH) = aP_NH
    bQindex(cP_NH) = bP_NH
  END IF

  IF (aP_ZR > 0 .AND. bP_ZR > 0) THEN
    cnscalar = cnscalar + 1
    cP_ZR    = cnscalar
    cqnames(cP_ZR) = aqnames(aP_ZR)
    cqdescp(cP_ZR) = aqdescp(aP_ZR)
    aQindex(cP_ZR) = aP_ZR
    bQindex(cP_ZR) = bP_ZR
  END IF

  IF (aP_ZI > 0 .AND. bP_ZI > 0) THEN
    cnscalar = cnscalar + 1
    cP_ZI    = cnscalar
    cqnames(cP_ZI) = aqnames(aP_ZI)
    cqdescp(cP_ZI) = aqdescp(aP_ZI)
    aQindex(cP_ZI) = aP_ZI
    bQindex(cP_ZI) = bP_ZI
  END IF

  IF (aP_ZS > 0 .AND. bP_ZS > 0) THEN
    cnscalar = cnscalar + 1
    cP_ZS    = cnscalar
    cqnames(cP_ZS) = aqnames(aP_ZS)
    cqdescp(cP_ZS) = aqdescp(aP_ZS)
    aQindex(cP_ZS) = aP_ZS
    bQindex(cP_ZS) = bP_ZS
  END IF

  IF (aP_ZG > 0 .AND. bP_ZG > 0) THEN
    cnscalar = cnscalar + 1
    cP_ZG    = cnscalar
    cqnames(cP_ZG) = aqnames(aP_ZG)
    cqdescp(cP_ZG) = aqdescp(aP_ZG)
    aQindex(cP_ZG) = aP_ZG
    bQindex(cP_ZG) = bP_ZG
  END IF

  IF (aP_ZH > 0 .AND. bP_ZH > 0) THEN
    cnscalar = cnscalar + 1
    cP_ZH    = cnscalar
    cqnames(cP_ZH) = aqnames(aP_ZH)
    cqdescp(cP_ZH) = aqdescp(aP_ZH)
    aQindex(cP_ZH) = aP_ZH
    bQindex(cP_ZH) = bP_ZH
  END IF

  RETURN
END SUBROUTINE min_set_scalars
