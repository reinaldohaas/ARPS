!#######################################################################
!
!     This file contains subroutine that interface ZXPLOT main
!     library package with the NCAR graphics low-level premitives.
!
!     Written by Ming Xue, CAPS, University of Oklahoma,
!     100 E. Boyd, Norman, OK 73019. mxue@uoknor.edu.
!
!     Copyright to Ming Xue. All rights reserved. 1986-1996.
!
!#######################################################################
!
      SUBROUTINE xdevic
c-----------------------------------------------------------------------
c To be compatible with old programs. Newly created code should call
c xdevic_new directly.
c-----------------------------------------------------------------------
      CHARACTER(LEN=256) filename
      INTEGER iwtype
      iwtype = 1
      CALL xdevic_new(iwtype,filename,0,13)
!          filename is un-initialized here so that it will behavior the same
!          as before (file name is "gmeta" or that set by environment variable NCARG_GKS_OUTPUT).
      RETURN
      END SUBROUTINE xdevic

      SUBROUTINE XDEVIC_NEW(IWTYPE,filename,len_filename,iounitin)
      INTEGER IWTYPE
      CHARACTER(LEN=256) filename  ! output graphical file name (base), the extension
                                   ! will be determined by IWTYPE
      INTEGER len_filename         ! valid character length in filename
                                   ! = 0 to use the default "gmeta" for backward compatibility.
      INTEGER iounitin             ! A file unit for output, used only by PS interface.
                                   ! It is ignored here.
      SAVE NZXCAL, MBORDR
      DATA NZXCAL /0/ ,MBORDR/0/
      COMMON /XPSIZE/ PSIZE
      COMMON /XPSD01/ XSIDE, YSIDE
      integer ncunique
      integer icoltable
      common /xcltbl/ ncunique,icoltable

      integer iasf(13)
      data iasf /13*1/

      PARAMETER (IERRF=6, LUNIT=2,  IWKID=1)

      CHARACTER(LEN=6) extstr
!
!#######################################################################
!
!     Set up divice by calling NCAR GKS package.
!
!#######################################################################
!
!     CALL opngks

      IF (IWTYPE < 1) IWTYPE = 1

      CALL GOPKS (IERRF, ISZDM)

      IF (len_filename > 0) THEN      ! else, use default file name 'gmeta' or set
                                      ! by environment variable NCARG_GKS_OUTPUT.
        SELECT CASE (IWTYPE)
        CASE (1)
          extstr = '.gmeta'
        CASE (11, 12)
          extstr = '.pdf  '
        CASE (20,23,26,29)
          extstr = '.ps   '
        CASE (21,24,27,30)
          extstr = '.eps  '
        CASE (22,25,28,31)
          extstr = '.epsi '
        CASE DEFAULT
          WRITE(extstr,'(a1,I5.5)') '.',IWTYPE
        END SELECT

        WRITE(6,'(/,1x,3a,/)') '==== Creating Graphical file: ',
     +          filename(1:len_filename)//TRIM(extstr),' ===='

        CALL NGSETC('ME',filename(1:len_filename)//TRIM(extstr))
!c     CALL NGSETI('Workstation',IWKID)
!c     CALL NGSETI('Full background',1)

      END IF

      CALL GOPWK (IWKID, LUNIT, IWTYPE)
      CALL GACWK (IWKID)

c
c#######################################################################
c
c     Turn off the clipping indicator (GKS routine)
c
c#######################################################################
c
      call gsclip(0)
c
c#######################################################################
c
c     Set all aspect source flags to "individual" (GKS routine)
c
c#######################################################################
c
      call gsasf (iasf)
c
c#######################################################################
C
C     Define 16 different color indices, for indices 0 through 15.
c     The color corresponding to index 0 is black and the color
c     corresponding to index 1 is white.
c
c#######################################################################
C
      CALL XICHAR
c
c#######################################################################
c
C     XDSPAC must prceeds all other calls to ZXPLOT routines.
c
c#######################################################################
c
      CALL XDSPAC( 1.0 )
      CALL XMINIT
      NZXCAL=1

      CALL XSETCLRS(icoltable,1)
      CALL xcolor(1)
c
      RETURN

      ENTRY XFRAME
c
c#######################################################################
c
C     Advance plotting onto next frame
c
c#######################################################################
c
        CALL XQMAP (  X1,X2,Y1,Y2 )
        CALL XLPNUP(X1,Y1)
        CALL FRAME
      RETURN

      ENTRY XGREND
c
c#######################################################################
c
C     Terminate graphic plotting
c
c#######################################################################
c
      IF( NZXCAL.EQ.0)  THEN
        PRINT*,' XGREND called before device is set up.'
        RETURN
      ENDIF
        CALL XQMAP (  X1,X2,Y1,Y2 )
        CALL XLPNUP(X1,Y1)
        call plotif(0.0, 0.0, 2)

        CALL FRAME
        CALL GDAWK (IWKID)
        CALL GCLWK (IWKID)
        CALL GCLKS
c        CALL clsgks
      RETURN

      ENTRY XFBORD(MBORD)
      MBORDR=MBORD
      RETURN

      END SUBROUTINE XDEVIC_NEW


      SUBROUTINE XDSPAC(PSIZE)
!
!#######################################################################
!
!     Define a normalized device (ND) space (0.0,XSIDE,0.0,YSIDE) on
!     device provided. The size of this space is 'PSIZE' times of the
!     max device  space.
!     YSIDE should be 1.0, XSIDE can be bigger or smaller that 1.0
!     depending on device. XSIDE, YSIDE should be passed to ZXPLOT
!     through common block PSIDES.
!
!     This routine sets up the ND space based on NCAR graphics.
!
!#######################################################################
!
!
      COMMON /XPSIZE/ PSIZE1
      COMMON /XPSD01/ XSIDE, YSIDE
      save xp1,xp2,yp1,yp2,x1,x2,y1,y2

      PSIZE1=PSIZE
C
      XSIDE=1.0
      YSIDE=1.0
      IF( PSIZE.LE.0.0) THEN
        PRINT*,'Zero or negative values for PSIZE not permitted!'
        PRINT*,'It should be in range of 0.1 to 1.,please reset value'
        STOP
      ENDIF

      EDGEX=0.5*XSIDE*(1.0/PSIZE-1)
      EDGEY=EDGEX
      X1=-EDGEX
      X2= EDGEX+XSIDE
      Y1=-EDGEY
      Y2= EDGEY+YSIDE

      xp1 = 0.0
      xp2 = 1.0
      yp1 = 0.0
      yp2 = 1.0

      CALL Set(xp1,xp2,yp1,yp2, X1,X2,Y1,Y2, 1)

      return
c
c#######################################################################
c
      entry xqset(xp1a,xp2a,yp1a,yp2a,x1a,x2a,y1a,y2a)

      xp1a = xp1
      xp2a = xp2
      yp1a = yp1
      yp2a = yp2

      x1a = x1
      x2a = x2
      y1a = y1
      y2a = y2

      RETURN
c
      entry xzx2ncar( xpos, ypos )
c
c#######################################################################
c
c     Find the corrdinate in the ncargraphic device space
c     given its corrdinate in the zxplot non-dimensional
c     device space.
c
c#######################################################################
c
      xpos = 0.0 + (xpos-x1)/(x2-x1)*1.0
      ypos = 0.0 + (ypos-y1)/(y2-y1)*1.0

      return
      END SUBROUTINE XDSPAC

      SUBROUTINE PPENUP(X1,Y1)
c
c#######################################################################
c
C     Connect pen up routines
c
c#######################################################################
c
        X=X1
        Y=Y1
        CALL Frstpt(X,Y)
      RETURN
      END SUBROUTINE PPENUP

      SUBROUTINE PPENDN(X1,Y1)
c
c#######################################################################
c
C     Connect pen down routines
c
c#######################################################################
c
        X=X1
        Y=Y1
        CALL Vector(X,Y)
      RETURN
      END SUBROUTINE PPENDN

      SUBROUTINE ZFILLN(X,Y,NP)
      integer np
      real x(np),y(np)
      CALL XFILAREA(x,y,np)
      RETURN
      END SUBROUTINE ZFILLN

      SUBROUTINE XFILAREA(X,Y,NP)
c
c#######################################################################
c
C     To fill a polygon defined by (x(i),y(i),i=1,np) with predefined
c     color (set by CALL COLOR).
c
c     Modified by Ming Xue to use GKS area fill routine GFA instead of
c     Ncar Graphics soft fill routine SFWRLD. (1/17/96).
c
c#######################################################################
c
      REAL X(*),Y(*)
      parameter (npmax=100000)
      REAL xra(npmax),yra(npmax)

      common /xwndw1/ xw1,xw2,yw1,yw2, iwndon
      real x1,x2,y1,y2
c
c#######################################################################
c
c     GSFAIS(INTS) is called first to set the area fill style, where
c     INTS=0 for hollow fill,
c         =1 for solid fill,
c         =2 for pattern fill and
C         =3 for hatch fill.
c
!c     Then GSAFSI(istyle) can be called to set the hatch or pattern style.
c     for INIT=2 and 3.
c     In NCAR graphics, hatch options are 1 for horizontal, 2 for vertical,
c     3 for right slanting, 4 for left slanting line, and 5 to horizontal and
c     vertical lines, 6 for right and left slanting lines.
c
c#######################################################################
c
      if(iwndon.eq.1)then
        nn=0
        do i=1,np-1
          x1=x(i)
          x2=x(i+1)
          y1=y(i)
          y2=y(i+1)
          call xtstwd(x1,y1,x2,y2,idispl)
          if(idispl.ne.0) then
            nn=nn+1
            IF( nn.gt.npmax) GOTO 999
            xra(nn)=x1
            yra(nn)=y1
            nn=nn+1
            IF( nn.gt.npmax) GOTO 999
            xra(nn)=x2
            yra(nn)=y2
          end if
        end do
      else
        IF( np.gt.npmax) then
          write(6,'(1x,a,/1x,a,i6)')
     :   'Work array xra and yra defined in XFILAREA not large enough.',
     :   'Only ',npmax,' number of pointed were used.'
        endif
        do i=1,min(npmax,np)
          xra(i) = x(i)
          yra(i) = y(i)
        enddo
        nn = min(npmax,np)
      end if

      if(nn .ge. 3) then
        DO I=1,nn
          CALL XTRANS(xra(I),yra(I))
        END DO
        CALL GFA(nn, xra, yra)
      endif

      RETURN
999   write(6,'(1x,a,/1x,a,i6)')
     :'Work array xra and yra defined in XFILAREA were too small',
     :'Please increase the array size to more than ',npmax

      RETURN
      END SUBROUTINE XFILAREA

      SUBROUTINE XICHAR
c
c#######################################################################
c
!C     To build an equivalence between READING amdahl CHAR function and
C     LONDON cray  CHAR function.
c
C     CHAR( I ) (london) =  CHAR( ICRAM(I) ) (reading)
C     ICHAR( C ) (reading)= ICRAM( ICHAR( C ) (london) )
c
c#######################################################################
c
      COMMON /XCHR30/ NCRAM
      INTEGER ICRAM(127) ,NCRAM(256)
      DATA ICRAM /
C 1-30
     :   32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32,
     :   32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32,
C 31-60
     :   32, 32, 90,127,123, 91,108, 80,125, 77, 93, 92, 78,107, 96,
     :   75, 97,240,241,242,243,244,245,246,247,248,249,122, 94, 76,
C 61-90
     :  126,110,111, 32,193,194,195,196,197,198,199,200,201,209,210,
     :  211,212,213,214,215,216,217,226,227,228,229,230,231,232,233,
C 91-120
     :  173,224,189,113,109,121,129,130,131,132,133,134,135,136,137,
     :  145,146,147,148,149,150,151,152,153,162,163,164,165,166,167,
C 121-127
     :  168,169,192, 32,208, 95, 32  /
      DO I=1,127
        NCRAM(I)=ICRAM(I)
      END DO
      DO  I=128, 256
        NCRAM(I)= 32
      END DO

      RETURN
      END SUBROUTINE XICHAR

      subroutine psftnm(a)
      return
      end subroutine psftnm

      subroutine psgray(a)
      return
      end subroutine psgray

      SUBROUTINE XDFCLRS
c
c#######################################################################
C
C     Define a set of RGB color triples for colors 1 through 15.
c
c#######################################################################
C
      DIMENSION RGBV(3,0:15)
c
c#######################################################################
C
C     Define the RGB color triples needed below.
c
c#######################################################################
C
      DATA RGBV /
     +               0.0 , 0.0  , 0.0  ,
     +              1.00 , 1.00 , 1.00 ,
     +              0.70 , 0.70 , 0.70 ,
     +              0.75 , 0.50 , 1.00 ,
     +              0.50 , 0.00 , 1.00 ,
     +              0.00 , 0.00 , 1.00 ,
     +              0.00 , 0.50 , 1.00 ,
     +              0.00 , 1.00 , 1.00 ,
     +              0.00 , 1.00 , 0.60 ,
     +              0.00 , 1.00 , 0.00 ,
     +              0.70 , 1.00 , 0.00 ,
     +              1.00 , 1.00 , 0.00 ,
     +              1.00 , 0.75 , 0.00 ,
     +              1.00 , 0.38 , 0.38 ,
     +              1.00 , 0.00 , 0.38 ,
     +              1.00 , 0.00 , 0.00 /
c
c#######################################################################
C
C     Define 16 different color indices, for indices 0 through 15.  The
C     color corresponding to index 0 is black and the color corresponding
C     to index 1 is white.
c
c#######################################################################
C
      DO I=0,15
        CALL GSCR (1,I,RGBV(1,I),RGBV(2,I),RGBV(3,I))
      END DO

      RETURN
      END SUBROUTINE XDFCLRS

      subroutine xafstyl(nstyle)
c
c#######################################################################
c
c     Set area fill style
c
c     nstyle = 0 hollow fill
c            = 1 solid fill
c            = 2 pattern fill
c            = 3 hatch fill
c
c#######################################################################
c
      nt = nstyle
      if( nstyle.lt.0 .or. nstyle.gt.3 ) nt = 0
c
      call gsfais( nt )

      return
      end subroutine xafstyl

      subroutine xafpatn(npat)
c
c#######################################################################
c
c     Set the hatch pattern when the fill style is hatch fill
c
c     npat = 1 horizontal lines
c          = 2 vertical lines
c          = 3 lines of positive slope
c          = 4 lines of negative slope
c          = 5 horizontal and vertical lines
c          = 6 lines of postive and negative slope
c
c#######################################################################
c
      np = npat
      if( npat.lt.1 .or. npat.gt.6 ) np = 1

      call gsfasi(np)

      return
      end subroutine xafpatn

      subroutine xlncinx(ind)
c
c#######################################################################
c
c     Set the line color index defined in xdfclrs
c
c#######################################################################
c
      indx = ind
      if( ind .lt.0 .or. ind.gt.15 ) indx = 0

      call gsplci(indx)
      return
      end subroutine xlncinx

      subroutine xafcinx(ind)
c
c#######################################################################
c
c     Set the area fill color index that is defined in xdfclrs
c
c#######################################################################
c
      indx = ind
      if( ind .lt.0 .or. ind.gt.15 ) indx = 0

      call gsfaci(indx)
      return
      end subroutine xafcinx

C#BEGINCOMMENTBLOCK
!
!     Ncargraphic strmline routine with modefications so that
!     originally hard wired parameters can be altered through
!     common blocks.
!
!
!
!      SUBROUTINE STRMLN (U,V,WORK,IMAX,IPTSX,JPTSY,NSET,IER)
!C
!C +-----------------------------------------------------------------+
!C |                                                                 |
!C |                Copyright (C) 1989 by UCAR                       |
!C |        University Corporation for Atmospheric Research          |
!C |                    All Rights Reserved                          |
!C |                                                                 |
!C |                 NCARGRAPHICS  Version 3.00                      |
!C |                                                                 |
!C +-----------------------------------------------------------------+
!C
!C SUBROUTINE STRMLN (U,V,WORK,IMAX,IPTSX,JPTSY,NSET,IER)
!C
!C DIMENSION OF           U(IMAX,JPTSY) , V(IMAX,JPTSY) ,
!C ARGUMENTS              WORK(2*IMAX*JPTSY)
!C
!C PURPOSE                STRMLN draws a streamline representation of
!C                        the flow field. The representation is
!C                        independent of the flow speed.
!C
!C USAGE                  If the following assumptions are met, use
!C
!C                            CALL EZSTRM  (U,V,WORK,IMAX,JMAX)
!C
!C                          Assumptions:
!C                            --The whole array is to be processed.
!C                            --The arrays are dimensioned
!C                              U(IMAX,JMAX) , V(IMAX,JMAX) and
!C                              WORK(2*IMAX*JMAX).
!C                            --Window and viewport are to be chosen
!C                              by STRMLN.
!C                            --PERIM is to be called.
!C
!C                        If these assumptions are not met, use
!C
!C                            CALL STRMLN (U,V,WORK,IMAX,IPTSX,JPTSY,
!C                                         NSET,IER)
!C
!C                        The user must call FRAME in the calling
!C                        routine.
!C
!C                        The user may change various internal
!C                        parameters via common blocks. See below.
!C
!C ARGUMENTS
!C
!C ON INPUT               U, V
!C                          Two dimensional arrays containing the
!C                          velocity fields to be plotted.
!C
!C                          Note:  If the U AND V components
!C                          are, for example, defined in Cartesian
!C                          coordinates and the user wishes to plot them
!C                          on a different projection (i.e., stereo-
!C                          graphic), then the appropriate
!C                          transformation must be made to the U and V
!C                          components via the functions FU and FV
!C                          (located in DRWSTR).
!C
!C                        WORK
!C                          User provided work array.  The dimension
!C                          of this array must be .GE. 2*IMAX*JPTSY.
!C
!C                          Caution:  This routine does not check the
!C                          size of the work array.
!C
!C                        IMAX
!C                          The first dimension of U and V in the
!C                          calling program. (X-direction)
!C
!C                        IPTSX
!C                          The number of points to be plotted in the
!C                          first subscript direction.  (X-direction)
!C
!C                        JPTSY
!C                          The number of points to be plotted in the
!C                          second subscript direction. (Y-direction)
!C
!C                        NSET
!C                          Flag to control scaling
!C                          > 0  STRMLN assumes that the window
!C                               and viewport have been set by the
!C                               user in such a way as to properly
!C                               scale the plotting instructions
!C                               generated by STRMLN. PERIM is not
!C                               called.
!C                          = 0  STRMLN will establish the window and
!C                               viewport to properly scale the
!C                               plotting instructions to the standard
!C                               configuration. PERIM is called to draw
!C                               the border.
!C                          < 0  STRMLN establishes the window
!C                               and viewport so as to place the
!C                               streamlines within the limits
!C                               of the user's window.  PERIM is
!C                               not called.
!C
!C ON OUTPUT              Only the IER argument may be changed. All
!C                        other arguments are unchanged.
!C
!C
!C                        IER
!C                          =  0 when no errors are detected
!C                          = -1 when the routine is called with ICYC
!C                               .NE. 0  and the data are not cyclic
!C                               (ICYC is an internal parameter
!C                               described below); in this case the
!C                               routine will draw the
!C                               streamlines with the non-cyclic
!C                               interpolation formulas.
!C
!C ENTRY POINTS           STRMLN, DRWSTR, EZSTRM, GNEWPT, CHKCYC
!C
!C COMMON BLOCKS          STR01, STR02, STR03, STR04
!C
!C REQUIRED LIBRARY       GRIDAL, GBYTES, and the SPPS
!C ROUTINES
!C
!C REQUIRED GKS LEVEL     0A
!C
!C I/O                    None
!C
!C PRECISION              Single
!C
!C LANGUAGE               FORTRAN 77
!C
!C HISTORY                Written and standardized in November 1973.
!C
!C                        Converted to FORTRAN 77 and GKS in June, 1984.
!C
!C
!C PORTABILITY            FORTRAN 77
!C
!C ALGORITHM              Wind components are normalized to the value
!C                        of DISPL. The least significant two
!C                        bits of the work array are
!C                        utilized as flags for each grid box. Flag 1
!C                        indicates whether any streamline has
!C                        previously passed through this box.  Flag 2
!C                        indicates whether a directional arrow has
!C                        already appeared in a box. Judicious use
!C                        of these flags prevents overcrowding of
!C                        streamlines and directional arrows.
!C                        Experience indicates that a final pleasing
!C                        picture is produced when streamlines are
!C                        initiated in the center of a grid box. The
!C                        streamlines are drawn in one direction then
!C                        in the opposite direction.
!C
!C REFERENCE              The techniques utilized here are described
!C                        in an article by Thomas Whittaker (U. of
!C                        Wisconsin) which appeared in the notes and
!C                        correspondence section of Monthly Weather
!C                        Review, June 1977.
!C
!C TIMING                 Highly variable
!C                          It depends on the complexity of the
!C                          flow field and the parameters:  DISPL,
!C                          DISPC , CSTOP , INITA , INITB , ITERC ,
!C                          and IGFLG. (See below for a discussion
!C                          of these parameters.) If all values
!C                          are default, then a simple linear
!C                          flow field for a 40 x 40 grid will
!C                          take about 0.4 seconds on the CRAY1-A;
!C                          a fairly complex flow field will take about
!C                          1.5 seconds on the CRAY1-A.
!C
!C
!C INTERNAL PARAMETERS
!C
!C                        NAME     DEFAULT            FUNCTION
!C                        ----     -------            --------
!C
!C                        EXT       0.25      Lengths of the sides of the
!C                                            plot are proportional to
!C                                            IPTSX and JPTSY except in
!C                                            the case when MIN(IPTSX,JPT)
!C                                            / MAX(IPTSX,JPTSY) .LT. EXT;
!C                                            in that case a square
!C                                            graph is plotted.
!C
!C                        SIDE      0.90      Length of longer edge of
!C                                            plot. (See also EXT.)
!C
!C                        XLT       0.05      Left hand edge of the plot.
!C                                            (0.0 = left edge of frame)
!C                                            (1.0 = right edge of frame)
!C
!C                        YBT       0.05      Bottom edge of the plot.
!C                                            (0.0 = bottom ; 1.0 = top)
!C
!C                                            (YBT+SIDE and XLT+SIDE must
!C                                            be .LE. 1. )
!C
!C                        INITA     2         Used to precondition grid
!C                                            boxes to be eligible to
!C                                            start a streamline.
!C                                            For example, a value of 4
!C                                            means that every fourth
!C                                            grid box is eligible ; a
!C                                            value of 2 means that every
!C                                            other grid box is eligible.
!C                                            (see INITB)
!C
!C                        INITB     2         Used to precondition grid
!C                                            boxes to be eligible for
!C                                            direction arrows.
!C                                            If the user changes the
!C                                            default values of INITA
!C                                            and/or INITB, it should
!C                                            be done such that
!C                                            MOD(INITA,INITB) = 0 .
!C                                            For a dense grid try
!C                                            INITA=4 and INITB=2 to
!C                                            reduce the CPU time.
!C
!C                        AROWL     0.33      Length of direction arrow.
!C                                            For example, 0.33 means
!C                                            each directional arrow will
!C                                            take up a third of a grid
!C                                            box.
!C
!C                        ITERP     35        Every 'ITERP' iterations
!C                                            the streamline progress
!C                                            is checked.
!C
!C                        ITERC     -99       The default value of this
!C                                            parameter is such that
!C                                            it has no effect on the
!C                                            code. When set to some
!C                                            positive value, the program
!C                                            will check for streamline
!C                                            crossover every 'ITERC'
!C                                            iterations. (The routine
!C                                            currently does this every
!C                                            time it enters a new grid
!C                                            box.)
!C                                            Caution:  When this
!C                                            parameter is activated,
!C                                            CPU time will increase.
!C
!C                        IGFLG     0         A value of zero means that
!C                                            the sixteen point Bessel
!C                                            Interpolation Formula will
!C                                            be utilized where possible;
!C                                            when near the grid edges,
!C                                            quadratic and bi-linear
!C                                            interpolation  will be
!C                                            used. This mixing of
!C                                            interpolation schemes can
!C                                            sometimes cause slight
!C                                            raggedness near the edges
!C                                            of the plot.  If IGFLG.NE.0,
!C                                            then only the bilinear
!C                                            interpolation formula
!C                                            is used; this will generally
!C                                            result in slightly faster
!C                                            plot times but a less
!C                                            pleasing plot.
!C
!C                        IMSG      0         If zero, then no missing
!C                                            U and V components are
!C                                            present.
!C                                            If .NE. 0, STRMLN will
!C                                            utilize the
!C                                            bi-linear interpolation
!C                                            scheme and terminate if
!C                                            any data points are missing.
!C
!C                        UVMSG     1.E+36    Value assigned to a missing
!C                                            point.
!C
!C                        ICYC      0         Zero means the data are
!C                                            non-cyclic in the X
!C                                            direction.
!C                                            If .NE 0, the
!C                                            cyclic interpolation
!C                                            formulas will be used.
!C                                            (Note:  Even if the data
!C                                            are cyclic in X, leaving
!C                                            ICYC = 0 will do no harm.)
!C
!C                        DISPL     0.33      The wind speed is
!C                                            normalized to this value.
!C                                            (See the discussion below.)
!C
!C                        DISPC     0.67      The critical displacement.
!C                                            If after 'ITERP' iterations
!C                                            the streamline has not
!C                                            moved this distance, the
!C                                            streamline will be
!C                                            terminated.
!C
!C                        CSTOP     0.50      This parameter controls
!C                                            the spacing between
!C                                            streamlines.  The checking
!C                                            is done when a new grid
!C                                            box is entered.
!C
!C DISCUSSION OF          Assume a value of 0.33 for DISPL.  This
!C DISPL,DISPC            means that it will take three steps to move
!C AND CSTOP              across one grid box if the flow was all in the
!C                        X direction. If the flow is zonal, then a
!C                        larger value of DISPL is in order.
!C                        If the flow is highly turbulent, then
!C                        a smaller value is in order.  The smaller
!C                        DISPL, the more the CPU time.  A value
!C                        of 2 to 4 times DISPL is a reasonable value
!C                        for DISPC.  DISPC should always be greater
!C                        than DISPL. A value of 0.33 for CSTOP would
!C                        mean that a maximum of three stream-
!C                        lines will be drawn per grid box. This max
!C                        will normally only occur in areas of singular
!C                        points.
!C
!C                                            ***************************
!C                                            Any or all of the above
!C                                            parameters may be changed
!C                                            by utilizing common blocks
!C                                            STR02 and/or STR03
!C                                            ***************************
!C
!C                        UXSML               A number which is small
!C                                            compared to the average
!C                                            normalized u component.
!C                                            Set automatically.
!C
!C                        NCHK      750       This parameter is located
!C                                            in DRWSTR. It specifies the
!C                                            length of the circular
!C                                            lists  used for checking
!C                                            for STRMLN crossovers.
!C                                            For most plots this number
!C                                            may be reduced to 500
!C                                            or less and the plots will
!C                                            not be altered.
!C
!C                        ISKIP               Number of bits to be
!C                                            skipped to get to the
!C                                            least two significant bits
!C                                            in a floating point number.
!C                                            The default value is set to
!C                                            I1MACH(5) - 2 . This value
!C                                            may have to be changed
!C                                            depending on the target
!C                                            computer; see subroutine
!C                                            DRWSTR.
!C
!C
!C
!      DIMENSION       U(IMAX,JPTSY)         ,V(IMAX,JPTSY)           ,
!     1                WORK(IMAX*JPTSY+1)
!      DIMENSION       WNDW(4)               ,VWPRT(4)
!C
!      COMMON /STR01/  IS         ,IEND      ,JS        ,JEND
!     1             ,  IEND1      ,JEND1     ,I         ,J
!     2             ,  X          ,Y         ,DELX      ,DELY
!     3             ,  ICYC1      ,IMSG1     ,IGFL1
!      COMMON /STR02/  EXT , SIDE , XLT , YBT
!      COMMON /STR03/  INITA , INITB , AROWL , ITERP , ITERC , IGFLG
!     1             ,  IMSG , UVMSG , ICYC , DISPL , DISPC , CSTOP
!C
!      SAVE
!C
!C THE FOLLOWING CALL IS FOR MONITORING LIBRARY USE AT NCAR
!C
!      CALL Q8QST4 ( 'GRAPHX', 'STRMLN', 'STRMLN', 'VERSION 01')
!C
!      IER = 0
!C
!C LOAD THE COMMUNICATION COMMON BLOCK WITH PARAMETERS
!C
!      IS = 1
!      IEND = IPTSX
!      JS = 1
!      JEND = JPTSY
!      IEND1 = IEND-1
!      JEND1 = JEND-1
!      IEND2 = IEND-2
!      JEND2 = JEND-2
!      XNX = FLOAT(IEND-IS+1)
!      XNY = FLOAT(JEND-JS+1)
!      ICYC1 = ICYC
!      IGFL1 = IGFLG
!      IMSG1 = 0
!C
!C IF ICYC .NE. 0 THEN CHECK TO MAKE SURE THE CYCLIC CONDITION EXISTS.
!C
!      IF (ICYC1.NE.0) CALL CHKCYC  (U,V,IMAX,JPTSY,IER)
!C
!C SAVE ORIGINAL  NORMALIZATION TRANSFORMATION NUMBER
!C
!      CALL GQCNTN ( IERR,NTORIG )
!C
!C SET UP SCALING
!C
!      IF (NSET) 10 , 20 , 60
!   10 CALL GETUSV ( 'LS' , ITYPE )
!      CALL GQNT ( NTORIG,IERR,WNDW,VWPRT )
!      CALL GETUSV('LS',IOLLS)
!      X1 = VWPRT(1)
!      X2 = VWPRT(2)
!      Y1 = VWPRT(3)
!      Y2 = VWPRT(4)
!      X3 = IS
!      X4 = IEND
!      Y3 = JS
!      Y4 = JEND
!      GO TO  55
!C
!   20 ITYPE = 1
!      X1 = XLT
!      X2 = (XLT+SIDE)
!      Y1 = YBT
!      Y2 = (YBT+SIDE)
!      X3 = IS
!      X4 = IEND
!      Y3 = JS
!      Y4 = JEND
!      IF (AMIN1(XNX,XNY)/AMAX1(XNX,XNY).LT.EXT) GO TO  50
!      IF (XNX-XNY)  30, 50, 40
!   30 X2 = (SIDE*(XNX/XNY) + XLT)
!      GO TO  50
!   40 Y2 = (SIDE*(XNY/XNX) + YBT)
!   50 CONTINUE
!C
!C CENTER THE PLOT
!C
!      DX = 0.25*( 1. - (X2-X1) )
!      DY = 0.25*( 1. - (Y2-Y1) )
!      X1 = (XLT+DX)
!      X2 = (X2+DX )
!      Y1 = (YBT+DY)
!      Y2 = (Y2+DY )
!C
!   55 CONTINUE
!C
!C SAVE NORMALIZATION TRANSFORMATION 1
!C
!      CALL GQNT ( 1,IERR,WNDW,VWPRT )
!C
!C DEFINE AND SELECT NORMALIZATION TRANS, SET LOG SCALING
!C
!      CALL SET(X1,X2,Y1,Y2,X3,X4,Y3,Y4,ITYPE)
!C
!      IF (NSET.EQ.0) CALL PERIM (1,0,1,0)
!C
!   60 CONTINUE
!C
!C DRAW THE STREAMLINES
!C .   BREAK THE WORK ARRAY INTO TWO PARTS.  SEE DRWSTR FOR FURTHER
!C .   COMMENTS ON THIS.
!C
!      CALL DRWSTR (U,V,WORK(1),WORK(IMAX*JPTSY+1),IMAX,JPTSY)
!C
!C RESET NORMALIATION TRANSFORMATION 1 TO ORIGINAL VALUES
!C
!      IF (NSET .LE. 0) THEN
!        CALL SET(VWPRT(1),VWPRT(2),VWPRT(3),VWPRT(4),
!     -           WNDW(1),WNDW(2),WNDW(3),WNDW(4),IOLLS)
!      ENDIF
!      CALL GSELNT (NTORIG)
!C
!      RETURN
!      END
!
!      BLOCK DATA
!c
!c#######################################################################
!c
!c     Initialize parameters for STRMLN.
!c
!c#######################################################################
!c
!      COMMON /STR01/  IS         ,IEND      ,JS        ,JEND
!     1             ,  IEND1      ,JEND1     ,I         ,J
!     2             ,  X          ,Y         ,DELX      ,DELY
!     3             ,  ICYC1      ,IMSG1     ,IGFL1
!      COMMON /STR02/  EXT , SIDE , XLT , YBT
!      COMMON /STR03/  INITA , INITB , AROWL , ITERP , ITERC , IGFLG
!     1             ,  IMSG , UVMSG , ICYC , DISPL , DISPC , CSTOP
!C
!      DATA EXT  / 0.25/,
!     :     SIDE / 0.90/,
!     :     XLT  / 0.05/,
!     :     YBT  / 0.05/,
!     :     INITA/ 6/,
!     :     INITB/ 6/,
!     :     AROWL/ 0.33/,
!     :     ITERP/ 35/,
!     :     ITERC/ -99/,
!     :     IGFLG/ 0/,
!     :     ICYC / 0/,
!     :     IMSG / 0/,
!     :     UVMSG/ 1.E+36/,
!     :     DISPL/ 0.33/,
!     :     DISPC/ 0.67/,
!     :     CSTOP/ 0.50/
!
!      END
C#ENDCOMMENTBLOCK

      SUBROUTINE strmln_new (u,v,work,imax,iptsx,jptsy,nset,
     :                       magmax,magmin,ier)
!#######################################################################
!
! PURPOSE: Plot stream line using NCARG STINIT/STREAM.
!          It is a replacement of the obsolete STRMLN entry point.
!
!#######################################################################
!
! AUTHOR: Yunheng Wang (11/05/2009)
!
!----------------------------------------------------------------------
!
        IMPLICIT NONE
        INTEGER, INTENT(IN)  :: imax, iptsx, jptsy
        INTEGER, INTENT(IN)  :: nset   ! whether to color the stream lines
                                       ! with magitude and the number of color levels
        REAL,    INTENT(IN)  :: u(imax,jptsy), v(imax,jptsy)
        REAL,    INTENT(OUT) :: work(3*iptsx*jptsy)
        REAL,    INTENT(OUT) :: magmax, magmin
        INTEGER, INTENT(OUT) :: ier
!----------------------------------------------------------------------

        INTEGER :: iwrk, nlv
        INTEGER :: i, icolor

!
!     Parameters for color pallete plotting.
!
        INTEGER, PARAMETER :: nctrlvls_max=1000    ! Max. number of contour values
        REAL    :: ctrlvls(nctrlvls_max)     ! contour values dividing the filled areas
        INTEGER :: clrindx(nctrlvls_max)     ! plot color index bar color index
        INTEGER :: nctrlvls                  ! Number of contour levels

        COMMON /xcflvls/nctrlvls,ctrlvls,clrindx

        REAL    :: tvl

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        CALL STSETI('SET -- Set Call Flag', 0)

        CALL STSETR('USV -- U Array Special Value', -9999.0)
        CALL STSETR('VSV -- V Array Special Value', -9999.0)
        CALL STSETI('SVF -- Special Value Flag', 1)


        IF (nset > 1) THEN
          nlv = MIN(nset,255)

          CALL STSETI('CTV -- Color Threshold Value Control', 1)
          CALL STSETI('NLV -- Number of Color Levels', nlv)

          CALL XQCOLOR(icolor)
          icolor = icolor - 1

          DO i = 1,nlv
            CALL STSETI('PAI -- Parameter Array Index',i)
            CALL STSETI('CLR -- GKS Color Index',i+icolor)
          END DO

        ELSE
          CALL STSETI('CTV -- Color Threshold Value Control', 0)
        END IF

        iwrk = 3*iptsx*jptsy
        CALL STINIT(u,imax,v,imax,work,imax,iptsx,jptsy,work,iwrk)

        CALL STREAM(u,v,work,ier,ier,work)

!----------------------------------------------------------------------
!
!  Save values for plotting color palette
!
!----------------------------------------------------------------------

        IF (nset > 1) THEN

          nctrlvls = nlv + 1

          DO i = 1,nlv
            CALL STSETI('PAI -- Parameter Array Index',i)
            CALL STGETR('TVL -- GKS Color Index',tvl)

            clrindx(min(i,  nctrlvls_max-1)) = icolor+i
            ctrlvls(min(i+1,nctrlvls_max-1)) = tvl
          END DO
          ctrlvls(1) = 2*ctrlvls(2) - ctrlvls(3)

          magmax = ctrlvls(nlv+1)
          magmin = ctrlvls(1)

        END IF

        RETURN
      END SUBROUTINE strmln_new

      SUBROUTINE XWRTCTBL(rgb_table,nc_max,exchange)

c#######################################################################
c
c     Activate pre-defined color table in rgb_table
c
c#######################################################################

      integer exchange
      integer nc_max,j,jj
      real rgb_table(3,nc_max)
      integer ncoltable

      CALL XQCLRTBL( ncoltable )

      if( ncoltable.lt.0 ) goto 100

      if( ncoltable.eq.0 ) then

        rgb_table(1,1) = 0.
        rgb_table(2,1) = 0.
        rgb_table(3,1) = 0.
        DO j=2,nc_max
          rgb_table(1,j) = 1.
          rgb_table(2,j) = 1.
          rgb_table(3,j) = 1.
        END DO

      else if( ncoltable.ne.4. ) then

        rgb_table(1,1) = 0.
        rgb_table(2,1) = 0.
        rgb_table(3,1) = 0.
        rgb_table(1,2) = 1.
        rgb_table(2,2) = 1.
        rgb_table(3,2) = 1.

      end if

      IF (exchange == 1) THEN
        rgb_table(1,1) = 1.
        rgb_table(2,1) = 1.
        rgb_table(3,1) = 1.
        rgb_table(1,2) = 0.
        rgb_table(2,2) = 0.
        rgb_table(3,2) = 0.
      END IF

100   continue

      DO j=1,nc_max
        jj = j -1
        CALL gscr(1,jj,rgb_table(1,j),rgb_table(2,j),rgb_table(3,j))
      END DO

      CALL gsplci(1)
      CALL gspmci(1)
      CALL gstxci(1)

      CALL xafstyl(1)  ! Set default area fill pattern to solid

      RETURN
      END SUBROUTINE XWRTCTBL

      SUBROUTINE COLOR (icolor)
      integer icolor
        call xcolor(icolor)
      RETURN
      END SUBROUTINE COLOR

      SUBROUTINE XCOLOR(icolor)
c
c#######################################################################
c
c     PURPOSE:
c
c     Set the color index for line, text, marker and area fill.
c
c#######################################################################
c
c     INPUT:
c
c     icolor  color index corresponding to the color table set by
c             setcolors.
c
c#######################################################################
c
c     Variable Declarations.
c
c#######################################################################
c
      implicit none
c
      integer icolor,ii
      integer kolor,lcolor
      save kolor
      data kolor /1/
c
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
C     Beginning of executable code...
C
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
c
c   To set the polyline color index, use

      kolor = icolor

      CALL SFLUSH            ! flush line segments before changing color
      ii = Min (icolor,299)
      ii = Max (ii,0)
      CALL GSPLCI (ii)

c   To set the polymarker color index, use
      CALL GSPMCI (ii)

c   To set the text color index, use
      CALL GSTXCI (ii)

c   To set the fill area color index, use
      CALL GSFACI (ii)

      RETURN

      ENTRY XQCOLOR(lcolor)
      lcolor = kolor
      RETURN
      END SUBROUTINE XCOLOR

      SUBROUTINE xpsfn(fn)
      character fn*(*)
      return
      END SUBROUTINE xpsfn

      SUBROUTINE xpaprlnth(yside)
      real yside
      return
      END SUBROUTINE xpaprlnth

      SUBROUTINE xpscmnt(ch)
      character*(*) ch
      return
      END SUBROUTINE xpscmnt

