!##################################################################
!######                                                      ######
!######                 PLOT_COLTAB                          ######
!######                                                      ######
!######                Copyright (c) 1996                    ######
!######    Center for Analysis and Prediction of Storms      ######
!######    University of Oklahoma.  All rights reserved.     ######
!######                                                      ######
!##################################################################
!##################################################################
!##################################################################
!

PROGRAM plot_coltab
  IMPLICIT NONE
  INTEGER :: coltab
  CHARACTER (LEN=132) :: ch
  CHARACTER (LEN=256) :: coltabfn
  INTEGER :: kcout, LEN
  PARAMETER(kcout=25)

  INTEGER :: i, num,  lind(kcout)
  REAL :: cline(kcout)
  REAL :: lblmag
  COMMON /labmag/ lblmag
  INTEGER :: icolor,lbcolor                ! required color
  COMMON /recolor/icolor,lbcolor
  REAL :: xl, xr, yb,yt

  CHARACTER(LEN=256) :: outfilename

  lblmag=1.0
  lbcolor=1

  CALL xdevic_new(1,outfilename,0,0)

  xl=0.10
  xr=0.90
  yb=0.01
  yt=0.99
  CALL xdspac(1.0)
  CALL xwindw(xl, xr, yb, yt)

  PRINT*,'enter color table (-1 or 1 or 2 or 4 or 5)'
  PRINT*,'       ps:  1 to 5 is same as arpsplt color map '
  PRINT*,'            -1 user self defined color map'
  READ(5,*) coltab
  IF(coltab == -1) THEN
    PRINT*,'Please enter color file name (quote the name) '
    READ(5,*) coltabfn
    CALL xstctfn(coltabfn)
    LEN = LEN_TRIM(coltabfn)
    WRITE(ch,'(''COLOR MAP FILE IS: '',A)') coltabfn(1:LEN)
    PRINT*,'The color map file is:',coltabfn(1:LEN)
  ELSE
    WRITE(ch,'(''COLOR MAP '',I2)') coltab
  END IF

  CALL setcolors ( coltab )
  CALL color(lbcolor)
  CALL xafstyl(1)
  CALL xartyp(2)
  CALL xhlfrq(20)

  CALL xchsiz(0.035*(yt-yb)* lblmag )
  IF(LEN > 50) CALL xchsiz(0.029*(yt-yb)* lblmag )
  LEN = 132
  CALL strmin(ch, LEN)
  CALL xcharc(0.5, 0.935, ch(1:LEN))

  CALL xchsiz(0.025*(yt-yb)* lblmag )
  yb=0.85
  yt=0.90
  DO num=0,9
    DO i=1,kcout
      lind(i)= num*kcout+i-1
      cline(i) = FLOAT(num*kcout+i-1)
    END DO
    CALL label(lind,cline,0.1,0.9,yb, yt ,kcout)
    yb=yb-0.08
    yt=yt-0.08
  END DO

  DO i = 1,6
    lind(i)= 249+i
    cline(i) = FLOAT (249+i)
  END DO
  CALL label(lind,cline,0.1,0.9,yb, yt ,6)
  CALL xframe
  CALL xgrend
  STOP
END PROGRAM plot_coltab

SUBROUTINE label(lind,cline,xl,xr,yb,yt, kcolor)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    Generate color label plots of 2-d field A given its
!      coordinate using ZXPLOT and ncar package..
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!    15/08/92
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!
!    lind     color index
!    kcout    number of color index
!    xl       Left bound of the physical domain
!    xr       Right bound of the physical domain
!    yb      label box bottom bound of the physical domain.
!    yt      label box top bound of the physical domain.
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: k,l
  INTEGER :: kcout
  PARAMETER (kcout=25)
  REAL :: xl,xr,yb,yt   ! the physical domain space left, right, bottom, top
  INTEGER :: lind(kcout)
  CHARACTER (LEN=20) :: llbs
  CHARACTER (LEN=20) :: ctmp
  REAL :: xra(5),yra(5)     ! array for single color box
  REAL :: dtx
  REAL :: x,y             ! color values position
  REAL :: xs,ys
  REAL :: cline(kcout)
  REAL :: lblmag
  COMMON /labmag/ lblmag
  INTEGER :: icolor,lbcolor                ! required color
  COMMON /recolor/icolor,lbcolor
  INTEGER :: kcolor
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  xs= xr-xl
  ys= yt-yb
  CALL xcfont(3)
  CALL xchsiz(0.25*ys* lblmag )
  dtx=xs/REAL(kcout)
  y=yb-0.37*ys
  DO k=1,kcolor
    xra(1)=xl+dtx*REAL(k-1)
    xra(2)=xl+dtx*REAL(k)
    xra(3)=xra(2)
    xra(4)=xra(1)
    xra(5)=xra(1)
    yra(1)=yt
    yra(2)=yra(1)
    yra(3)=yb
    yra(4)=yra(3)
    yra(5)=yra(1)
    CALL color(lind(k))
    CALL zfilln(xra,yra,5)
    IF (kcout < 20 .OR. (kcout >= 20 .AND. MOD(k-1,2) == 0) ) THEN
      CALL xrch(cline(k),ctmp,l)
      llbs=ctmp(1:l-2)
      CALL color(lbcolor)
      x = xra(1)
      CALL xcharl(x,y,llbs(1:l-2))
    END IF
  END DO
  CALL color(lbcolor)
  CALL xwdwof

  RETURN
END SUBROUTINE label

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

  CHARACTER (LEN=*) :: string ! A character string for the name of
                              ! this run.
  INTEGER :: length            ! The length of the non-blank part
                               ! of a string.

  CHARACTER (LEN=1) ::  str_1
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
    IF( .NOT. (str(i:i) == ' ' .AND.                                    &
          (str_1 == ' ' .OR. str_1 == '(' .OR. str_1 == '='))) THEN
      length=length+1
      string(length:length)=str(i:i)
    END IF

  END DO

  RETURN
END SUBROUTINE strmin

