!c
!c ZXPLOT plotting package developed by Ming Xue with
!c contributions from Zuojun Zhang.
!c The author reserves all rights to the package.
!c
!c Change history:
!c
!c Version 2.2: Window clipping and multilevel rotated rectangular
!c masking added. May, 1991.
!c
!c Version 2.2a:
!c All printing output are routed to fortran unit NCH. 11/30/91.
!c Subroutine xstchn(nch) added.
!c
!c Version 2.2b.
!c subroutine xqspac added. (4/3/1992)
!c
!c Version 2.2c  (10/25/1993)
!c Added bad value skipping capability in the contour routines.
!c Zhang Zuojun and Ming Xue.
!c
!c This file contains subroutines common to all versions (PS and NCAR graphics
!c version)
!c
!c Correction made in xaxinc (4/20/94)
!c
!c 8/31/1995
!c Wind vector unit is now subject to the limit of vmax also.
!c
!c 1/24/96.
!c Fixed a problem in XVECTU with the first guess of umax,umin,vmax,vmin
!c when the first value is missing.
!c
!c 2/3/96 (M.Xue)
!c fixed a problem in xnwpic when the rotation angle is 90.0 for xspace.
!c
!c 10/13/1998 (M. Xue)
!c Added cross-hatching routines (xhatcha etc.) at arbitary angles.
!c
!c Added routine XCTRHL to label H and L contour centers. It is
!c called within XCONTA when XHLLABL is called with argument 1 or 2.
!c H and L labeling is off by default.
!c
!c Reorganized the subroutine into three files. zxplot3.f contains
!c common subroutines for all versions (e.g. PS and NCARgraphics versions),
!c zxpslib3.f, zxnglib3.f and zxgenlib3.f contains version dependent
!c subroutines, and xncar3.f and xpost3.f contain package dependent
!c drivers.
!c
!c ZXPLOT is upgrade to Version 3.0
!c
!c 11/19/1998 (Ming Xue)
!c Added subroutine XSETCLRS and other color related routines.
!c
!c Streamline routines for color filled maps (XCOLFIL),
!c variable colored contours (mode=4 for XCONTA).
!c Included map projection setup routines in zxplot3.f.
!c Also include map plotting XDRAWMAP.
!c
!c Added color map number five.
!c Added wind barb plotting routines, XBARB and XBARBS.
!c
!c
!c To do list:
!c Use buffer to store line segments. Remove redundant movetos.
!c Streamline color palette plotting routine.
!c Wind barbs scaling still an issue.
!c Tune fond sizes for PS
!c
!c Work on ZXPLOT Version 3.0 User's Guide.
!c Document new routines.
!c
!c Known problems: Missing value handling is not quite right for
!c xhatcha.
!c
!c Known problem: xcontc sometimes run into infinite loop
!c and overwrite the arrays for storing the line segments
!c
!c Added subroutine xcontc1, a simpler but less efficient version
!c of contour color mapping routine. Entry xcontcopt added to
!c set the option for using xcontc or xcontc1. One can also
!c set the option via common block xcontc_opt.
!c xcontc1 does not share the problem with xcontc.
!c
!c Added an option for XCOLFIL to plot pixdel type color fill
!c of a 2-D field.  Call xcontcopt(3) before call xcolfil to
!c activate this option. XPIXELFIL is called by xcolfil in this
!c case instead of xcontc or xcontc1.
!c
!c Added subroutine XCTR_THICK_THIN_RATIO to change to default
!c line thinkness ratio between highlighted thick contours to
!c thin contours. The default is 2.
!c
!c Fix dead-loop problem with contouring routines (XCONTR and XCONTJ) (YHW)
!c
!c Keith Brewster, CAPS, April, 2009
!c Modified to add marker types 5-10 for filled symbols
!c
      SUBROUTINE XMINIT
!C To initialize ZXPLOT package (Called by XDEVIC when setting up device)
      CHARACTER CLABEL*20
      COMMON /XPSD01/ XSIDE, YSIDE
      COMMON /XPHY01/ PL,PR,PB,PT,XRANGE,YRANGE
      COMMON /XPHO03/ DXPO,DYPO
      COMMON /XMAP04/ X1,X2,Y1,Y2,XSCALE,YSCALE
      COMMON /XMAO05/ DXMOP,DYMOP
      COMMON /XFTR06/ XFACTR,YFACTR
      COMMON /XPRF07/ XPREF,YPREF
      COMMON /XMRF08/ XMREF,YMREF,XMPREF,YMPREF
      COMMON /XAGS09/ DRANG,CRANG,XANGLE,XSYMAN,SRANG,KSR,XA,YA
      COMMON /XSCS10/ SINDRA,COSDRA,SINMRA,COSMRA,SINSRA,COSSRA
     :               ,SINXA,COSXA,SINYA,COSYA ,CHSIN,CHCOS
      COMMON /XPEN11/ XPEN,YPEN,FLEN,BLEN,NPD,XMPEN,YMPEN
      COMMON /XCHA20/ HCTR,SCTR,CRATIO, KFONT,NUNDLN
      COMMON /XASC12/ IASCII(300)
      COMMON /XCHP21/ XCHPEN, YCHPEN ,XCHMO,YCHMO,XCHPO,YCHPO
      COMMON /XLPN13/ HF1,HB1,HF2,HB2,LFULL,lfull0,LTHICK,DTHICK
      COMMON /XLAB14/ DLABEL,WLABEL,HLABEL,SIZLB,KLBTYP,ICLI,ICLON
      COMMON /XLAB15/ CLABEL
      COMMON /XLAB16/ LCLAB
      COMMON /XLBA33/ LABROT
      COMMON /XCRF17/ CLREF,LCPTN,LABTYP,ICLF,LHILIT,IHLF,KCT0
      COMMON /XAXS18/ KANX,KANY, KTKX,KTKY
      COMMON /XCLM19/ NMIN, NMAX
      COMMON /XCDV23/ NSUBDV
      COMMON /XCMD24/ MTD
      COMMON /XCIR25/ XCIR(9) ,YCIR(9) , RPOINT
      COMMON /XPRJ26/ KPROJC
      COMMON /XCHR31/ CHDATA
      COMMON /XCHR32/ ICDATA
      COMMON /XFMT33/ LBFMT, AXFMT
      COMMON /XFMT34/ LLBFMT,LAXFMT
      COMMON /XHCH35/ DH
      COMMON /XART36/ KARTYP,KVMODE,VSC
      COMMON /XAXM35/ NTMAG, ANMAG, ANSIZ
      COMMON /XHLL36/ hllabel
      COMMON /ZCHOLE/ NHOLE,SPECIA,nvtrbadv
      common /xwndw1/ xw1,xw2,yw1,yw2, iwndon
      common /xmask1/ xm1(99),xm2(99),ym1(99),ym2(99),rmangl(99),
     :                cosmsa(99),sinmsa(99),lvlmsk
      common /xcwndw/ icwndw, xcpen, ycpen
      common /xoutch/ nch
      integer ncunique
      integer icoltable
      common /xcltbl/ ncunique,icoltable
      integer iclrbgn,iclrend  ! Beginning and ending colors of contours
      common /xctrclor/iclrbgn,iclrend
      real ctrmin, ctrmax
      common /xctrmx/ ctrmin, ctrmax

      integer nctrlvls_max
      parameter(nctrlvls_max=1000)    ! Max. number of contour values
      real ctrlvls(nctrlvls_max)     ! contour values dividing the filled areas
      integer clrindx(nctrlvls_max)  ! plot color index bar color index
      integer nctrlvls      ! Number of contour levels
      common /xcflvls/nctrlvls,ctrlvls,clrindx
      common /xfctr1/ fctr

      CHARACTER CHDATA(127)*300,LBFMT*50,AXFMT*10
      INTEGER ICDATA (0:150, 32:127)
      INTEGER NASCII(300)

      character cpalnfmt*15
      common /xcplnfmt/ cpalnfmt

      integer labmask
      common /labmask1/ labmask

      common /xlimzf/ limzf, zfmax, zfmin
      integer icontcopt
      common /xcontc_opt/ icontcopt

      integer ictr_thick_thin_ratio
      common /ctr_thick_thin_ratio/ ictr_thick_thin_ratio

      integer icplswitch
      common /xcplswitch/ icplswitch

C Note XSIDE and YSIDE are the length of the ND-space x and y sides.
C They should be defined outside this package passing through common
C block PSIDES. The prefered values are XSIDE=1.5, YSIDE=1.0
      PL=0.0
      PR=XSIDE
      PB=0.0
      PT=YSIDE
      XRANGE=PR-PL
      YRANGE=PT-PB
      fctr = sqrt( abs(xrange * yrange) )
      X1=0.0
      X2=1.0
      Y1=0.0
      Y2=1.0
      XSCALE=X2-X1
      YSCALE=Y2-Y1
      XFACTR=XRANGE/XSCALE
      YFACTR=YRANGE/YSCALE
      XMREF=X1
      YMREF=Y1
      XMPREF=PL
      YMPREF=PB
      XPREF=0.5*(PL+PR)
      XPREF=0.5*(PT+PB)
      DXPO=0.0
      DYPO=0.0
      DXMOP=0.0
      DYMOP=0.0
      XMPEN=XMREF
      YMPEN=YMREF
      HCTR=0.02
      SCTR=0.02/YFACTR
      CRATIO=0.75
      KFONT=2
      NUNDLN=0
      DRANG=0.0
      CRANG=0.0
      SRANG=90.0
      KSR=0
      XA=0.0
      YA=0.0
      SINDRA=0.0
      COSDRA=1.0
      SINMRA=0.0
      COSMRA=1.0
      SINSRA=1.0
      COSSRA=0.0
      SINXA=0.0
      COSXA=1.0
      SINYA=0.0
      COSYA=1.0

      XANGLE=CRANG+DRANG+(90-SRANG)*KSR+XA
      XSYMAN=0.0
      CHSIN=0.0
      CHCOS=1.0

      HHH  =0.001*YSIDE
      HF1=HHH*10
      HB1=HHH*5
      HF2=HHH*10
      HB2=HHH*5
      LFULL =1
      LFULL0=1
      LTHICK=1
      DTHICK=0.0007*YSIDE

      ICLI=3
      CLABEL=' '
      LCLAB=1
      DLABEL=XRANGE/ICLI
      HLABEL=0.015
      SIZLB =0.015
      LABROT=1
      KLBTYP=-1
      WLABEL=HLABEL*LCLAB*0.77
      ICLON=0

      CLREF=0.0
      LCPTN=0
      LABTYP=2
      ICLF=2
      LHILIT=1
      IHLF=4
      KCT0=1

      KANX=-1
      KTKX=1
      KANY=-1
      KTKY=1
      AXFMT= '*'
      LBFMT= '*'
      LLBFMT=1
      LAXFMT=1
      DH=0.015
      KARTYP=2
      KVMODE=1
      VSC=1.0
      NTMAG=0
      hllabel=0

      NMIN=8
      NMAX=20

      NSUBDV=4
      MTD=0
      iwndon=0
      icwndw=0
      lvlmsk=0

      nch = 6

      NHOLE=0
      nvtrbadv=0
      SPECIA=-9999.

      ctrmin = 0.0
      ctrmax = 0.0

      icoltable = 1
      iclrbgn = 1
      iclrend = 1

      nctrlvls=0

      cpalnfmt = '*'

      labmask = 0

      limzf=0
      zfmin=-999.0
      zfmax= 999.0

      icontcopt=1

      ictr_thick_thin_ratio = 2

      Print*,'icplswitch set to 1 in XMINIT'
      icplswitch = 1

      PI=4*ATAN(1.)
      DO I=1,9
        XCIR(I)=COS((I-1)*0.25*PI)
        YCIR(I)=SIN((I-1)*0.25*PI)
      END DO
      RPOINT=0.0010*YSIDE
      KPROJC=0
      DO I=1,300
        IASCII(I)=NASCII(I)
      END DO

      CALL XINTMKR

      GOTO ( 501, 502, 503, 504 ) KFONT
 501  CALL XCSETB(CHDATA)
      GOTO 505
 502  CALL XCSETC(CHDATA)
      GOTO 505
 503  CALL XCSETA(CHDATA)
      GOTO 505
 504  CALL XCSETD(CHDATA)
 505  CONTINUE
      DO I=32,127
        ICDATA(0,I)=0
      END DO
      RETURN

      DATA NASCII /
*  1-16
     +000,000,000,000,000,000,000,000,000,000,000,000,000,000,000,000,
*  17-32
     +000,000,000,000,000,000,000,000,000,000,000,000,000,000,000,000,
*  33-48
     +000,000,000,000,000,000,000,000,000,000,000,000,000,000,000,000,
*  49-64
     +000,000,000,000,000,000,000,000,000,000,000,000,000,000,000,32 ,
*  65-80
     +000,000,000,000,000,000,000,000,000,000,46 ,60 ,40 ,43 ,124,38 ,
*  81-96
     +000,000,000,000,000,000,000,000,000,33 ,36 ,42 ,41 ,59 ,126,45 ,
*  97-112
     +47 ,000,000,000,000,000,000,000,000,000,44 ,37 ,95 ,62 ,63 ,000,
*  113-128
     +94 ,000,000,000,000,000,000,000,96 ,58 ,35 ,64 ,000,61 ,34 ,000,
*  129-144
     +97 ,98 ,99 ,100,101,102,103,104,105,000,000,000,000,000,000,000,
*  145-160
     +106,107,108,109,110,111,112,113,114,000,000,000,000,000,000,000,
*  161-176
     +000,115,116,117,118,119,120,121,122,000,000,000,91 ,000,000,000,
*  177-192
     +000,000,000,000,000,000,000,000,000,000,000,000,93 ,000,000,123,
*  193-208
     +65 ,66 ,67 ,68 ,69 ,70 ,71 ,72 ,73 ,000,000,000,000,000,000,125,
*  209-224
     +74 ,75 ,76 ,77 ,78 ,79 ,80 ,81 ,82 ,000,000,000,000,000,000,92 ,
*  225-240
     +000,83 ,84 ,85 ,86 ,87 ,88 ,89 ,90 ,000,000,000,000,000,000,48 ,
*  241-256
     +49 ,50 ,51 ,52 ,53 ,54 ,55 ,56 ,57 ,000,000,000,000,000,000,000,
*  257-272
     +000,000,000,000,000,000,000,000,000,000,000,000,000,000,000,000,
*  273-288
     +000,000,000,000,000,000,000,000,000,000,000,000,000,000,000,000,
*  289-300
     +000,000,000,000,000,000,000,000,000,000,000,000             /
      END SUBROUTINE XMINIT
c
      SUBROUTINE XSETCLRS(col_tab)
c-----------------------------------------------------------------------
c To be back-compatible with old programs. Newly created code should call
c xsetclrs_new directly.
c-----------------------------------------------------------------------
      INTEGER col_tab
      CALL xsetclrs_new(col_tab,0)
      RETURN
      END SUBROUTINE XSETCLRS

      SUBROUTINE XSETCLRS_new(col_tab,exchange)
c
c#######################################################################
c
c     PURPOSE:
c
c     Setup the color tables for ZXPLOT.
c
c#######################################################################
c
c     AUTHOR: Min Zou
c     8/28/94
c
c     1/17/96 (Ming Xue).
c     Added call to xafsty to set the default style of area fill
c     that uses GFA.
c
c     1/20/96 (Min Zou)
c     Added grayscale color table (col_tab=4, and user-specfied
c     color table (read in from a file) options.
c
c     4/15/96 (Zuojun Zhang)
c     Added multi-color scales (col_tab=5) and eliminate loading
c     unreferenced color table definition in the previous version.
c
c     11/16/1998 (Ming Xue)
c     Unified the NCAR graphics and PS versions.
c
c#######################################################################
c
c
c     INPUT:
c
c     col_tab  = -1, color table read in from file coltabfn, which
c                    can be set using entry XSTCTFN.
c              = 0,  Black and white plot. All lines are black.
c              = 1,  Predefined color table No. 1.
c              = 2,  Predefined color table No. 2,
c                    which is No.1 in reversed order.
c              = 3,  Predefined color table No. 3.
c              = 4,  Gray scale color table.
c              = 5,  predefined 200 multi-spectum color table (zj).
c              = 6,  Color table containing reflectivity palet etc.
c
c#######################################################################
c
c     Variable Declarations.
c
c#######################################################################
c
      IMPLICIT NONE
c
      integer exchange
      integer i,nc_max,col_table_no,j,jj
      integer col_tab,sum_colors,lenfil
      parameter (nc_max=256,col_table_no=20)
      real rgbv(3,25),rgbf(3,18),rgbg(3,18),rgb_zj(3,200),rgb6(3,139)
      real rgb_table(3,nc_max)
      logical fexist
      data rgb_table /768*1./  ! initial set rgb

      character*256 coltabfn
      character col_tab_fn*(*)

      save coltabfn
      data coltabfn /'zx_color.tbl'/

      integer ncunique , ncoltable
      integer icoltable
      common /xcltbl/ ncunique,icoltable

c     Index: for color table 1
c 0=bLack 1=white 2=yellow 3=dark 4=royal 5=light 6=sky 7=turquoise
c 8=aqua 9=olive 10=yellow-green 11=light 12=bright 13=kelly 14=green
c 15=yellow 16=maize 17=orange 18=red-orange 19=bright 20=red 21=dark
c 22=brown c 23=violet 24= mauve
C
      data ((rgbv(i,j),i=1,3),j=1,25)/
     :         1.000, 1.000, 1.000,  0.000, 0.000, 0.000,
     :         1.000, 0.804, 0.000,  0.000, 0.000, 0.702,
     :         0.000, 0.353, 1.000,  0.000, 0.553, 1.000,
     :         0.000, 0.753, 1.000,  0.000, 1.000, 1.000,
     :         0.631, 1.000, 1.000,  0.073, 0.612, 0.015,
     :         0.087, 0.737, 0.018,  0.095, 0.799, 0.019,
     :         0.102, 0.862, 0.021,  0.110, 0.924, 0.022,
     :         0.119, 1.000, 0.024,  1.000, 1.000, 0.000,
     :         1.000, 0.804, 0.000,  1.000, 0.604, 0.000,
     :         1.000, 0.400, 0.000,  1.000, 0.000, 0.000,
     :         0.804, 0.000, 0.000,  0.604, 0.000, 0.000,
     :         0.400, 0.000, 0.000,  0.400, 0.000, 0.400,
     :         0.559, 0.085, 0.433 /

C Index for color table 3:
c 0=black 1=white 2=yellow 3=turquoise 4=carolina 5=blue 6=bright
c 7=green 8=dark 9=yellow 10=dark 11=orange 12=bright 13=red 14=dark
c 15=magenta 16=purple 17= white

      data ((rgbf(i,j),i=1,3),j=1,18)/
     :          1.000, 1.000, 1.000,  0.000, 0.000, 0.000,
     :          1.000, 0.804, 0.000,  0.000, 0.925, 0.925,
     :          0.004, 0.627, 0.961,  0.000, 0.000, 0.965,
     :          0.000, 1.000, 0.000,  0.000, 0.784, 0.000,
     :          0.000, 0.565, 0.000,  1.000, 1.000, 0.000,
     :          0.906, 0.753, 0.000,  1.000, 0.565, 0.000,
     :          1.000, 0.000, 0.000,  0.839, 0.000, 0.000,
     :          0.753, 0.000, 0.000,  1.000, 0.000, 1.000,
     :          0.600, 0.333, 0.788,  1.000, 1.000, 1.000 /
c
C     Index for color table 5:
      data ((rgb_zj(i,j),i=1,3),j=1,38)/
     :          1.000, 1.000, 1.000, 0.000, 0.000, 0.000,
     :          1.000, 0.000, 0.000, 0.000, 1.000, 0.000,
     :          0.000, 0.000, 1.000, 1.000, 1.000, 0.000,
     :          0.000, 1.000, 1.000, 1.000, 0.000, 1.000,
     :          0.631, 1.000, 1.000, 0.000, 0.000, 0.702,
     :          0.300, 0.000, 0.000, 0.478, 0.000, 0.000,
     :          0.656, 0.000, 0.000, 0.834, 0.000, 0.000,
     :          1.000, 0.006, 0.006, 0.992, 0.198, 0.000,
     :          1.000, 0.356, 0.012, 1.000, 0.517, 0.029,
     :          1.000, 0.688, 0.036, 1.000, 0.883, 0.019,
     :          0.902, 1.000, 0.178, 0.818, 1.000, 0.440,
     :          0.813, 1.000, 0.623, 0.848, 1.000, 0.766,
     :          0.930, 1.000, 0.883, 0.883, 1.000, 0.930,
     :          0.766, 1.000, 0.848, 0.623, 1.000, 0.813,
     :          0.440, 1.000, 0.818, 0.178, 1.000, 0.902,
     :          0.019, 0.883, 1.000, 0.036, 0.688, 1.000,
     :          0.029, 0.517, 1.000, 0.012, 0.356, 1.000,
     :          0.000, 0.198, 0.992, 0.006, 0.006, 1.000,
     :          0.000, 0.000, 0.834, 0.000, 0.000, 0.656/
      data ((rgb_zj(i,j),i=1,3),j=39,76)/
     :          0.000, 0.000, 0.478, 0.000, 0.000, 0.300,
     :          0.300, 0.000, 0.000, 0.567, 0.000, 0.000,
     :          0.834, 0.000, 0.000, 1.000, 0.051, 0.051,
     :          1.000, 0.184, 0.184, 1.000, 0.358, 0.277,
     :          1.000, 0.529, 0.373, 1.000, 0.698, 0.471,
     :          1.000, 0.859, 0.577, 1.000, 1.000, 0.702,
     :          0.926, 1.000, 0.777, 0.795, 1.000, 0.641,
     :          0.644, 1.000, 0.525, 0.483, 1.000, 0.419,
     :          0.317, 1.000, 0.317, 0.184, 1.000, 0.184,
     :          0.051, 1.000, 0.051, 0.000, 0.834, 0.000,
     :          0.000, 0.567, 0.000, 0.000, 0.300, 0.000,
     :          0.000, 0.000, 0.300, 0.000, 0.000, 0.567,
     :          0.000, 0.000, 0.834, 0.050, 0.050, 1.000,
     :          0.184, 0.184, 1.000, 0.353, 0.282, 1.000,
     :          0.520, 0.382, 1.000, 0.680, 0.489, 1.000,
     :          0.831, 0.605, 1.000, 0.958, 0.745, 1.000,
     :          1.000, 0.745, 0.958, 1.000, 0.605, 0.831,
     :          1.000, 0.489, 0.680, 1.000, 0.382, 0.520,
     :          1.000, 0.282, 0.353, 1.000, 0.184, 0.184/
      data ((rgb_zj(i,j),i=1,3),j=77,114)/
     :          1.000, 0.051, 0.051, 0.834, 0.000, 0.000,
     :          0.567, 0.000, 0.000, 0.300, 0.000, 0.000,
     :          0.000, 0.300, 0.000, 0.000, 0.567, 0.000,
     :          0.000, 0.834, 0.000, 0.051, 1.000, 0.051,
     :          0.184, 1.000, 0.184, 0.277, 1.000, 0.358,
     :          0.373, 1.000, 0.529, 0.471, 1.000, 0.698,
     :          0.577, 1.000, 0.859, 0.702, 1.000, 1.000,
     :          0.777, 0.926, 1.000, 0.641, 0.795, 1.000,
     :          0.525, 0.644, 1.000, 0.419, 0.483, 1.000,
     :          0.317, 0.317, 1.000, 0.184, 0.184, 1.000,
     :          0.050, 0.050, 1.000, 0.000, 0.000, 0.834,
     :          0.000, 0.000, 0.567, 0.000, 0.000, 0.300,
     :          0.901, 0.901, 0.901, 0.812, 0.812, 0.812,
     :          0.723, 0.723, 0.723, 0.634, 0.634, 0.634,
     :          0.545, 0.545, 0.545, 0.456, 0.456, 0.456,
     :          0.367, 0.367, 0.367, 0.278, 0.278, 0.278,
     :          0.189, 0.189, 0.189, 0.100, 0.100, 0.100,
     :          0.200, 0.000, 1.000, 0.349, 0.000, 1.000,
     :          0.503, 0.000, 1.000, 0.676, 0.000, 1.000/
      data ((rgb_zj(i,j),i=1,3),j=115,152)/
     :          0.880, 0.000, 1.000, 1.000, 0.100, 0.880,
     :          1.000, 0.200, 0.676, 1.000, 0.300, 0.503,
     :          1.000, 0.200, 0.349, 1.000, 0.100, 0.200,
     :          1.000, 0.200, 0.000, 1.000, 0.349, 0.000,
     :          1.000, 0.503, 0.000, 1.000, 0.676, 0.000,
     :          1.000, 0.880, 0.000, 0.880, 1.000, 0.000,
     :          0.676, 1.000, 0.000, 0.503, 1.000, 0.000,
     :          0.349, 1.000, 0.000, 0.200, 1.000, 0.000,
     :          0.100, 1.000, 0.200, 0.200, 1.000, 0.349,
     :          0.300, 1.000, 0.503, 0.200, 1.000, 0.676,
     :          0.100, 1.000, 0.880, 0.000, 0.880, 1.000,
     :          0.000, 0.676, 1.000, 0.000, 0.503, 1.000,
     :          0.000, 0.349, 1.000, 0.000, 0.200, 1.000,
     :          0.300, 0.000, 0.000, 0.567, 0.000, 0.000,
     :          0.834, 0.000, 0.000, 1.000, 0.051, 0.051,
     :          1.000, 0.184, 0.184, 1.000, 0.317, 0.317,
     :          1.000, 0.451, 0.451, 1.000, 0.585, 0.585,
     :          1.000, 0.718, 0.718, 1.000, 0.851, 0.851,
     :          0.702, 1.000, 1.000, 0.436, 1.000, 1.000/
      data ((rgb_zj(i,j),i=1,3),j=153,190)/
     :          0.169, 1.000, 1.000, 0.000, 0.951, 0.951,
     :          0.000, 0.817, 0.817, 0.000, 0.684, 0.684,
     :          0.000, 0.551, 0.551, 0.000, 0.417, 0.417,
     :          0.000, 0.284, 0.284, 0.000, 0.150, 0.150,
     :          0.000, 0.300, 0.000, 0.000, 0.567, 0.000,
     :          0.000, 0.834, 0.000, 0.051, 1.000, 0.051,
     :          0.184, 1.000, 0.184, 0.317, 1.000, 0.317,
     :          0.451, 1.000, 0.451, 0.585, 1.000, 0.585,
     :          0.718, 1.000, 0.718, 0.851, 1.000, 0.851,
     :          0.851, 0.851, 1.000, 0.718, 0.718, 1.000,
     :          0.585, 0.585, 1.000, 0.451, 0.451, 1.000,
     :          0.317, 0.317, 1.000, 0.184, 0.184, 1.000,
     :          0.050, 0.050, 1.000, 0.000, 0.000, 0.834,
     :          0.000, 0.000, 0.567, 0.000, 0.000, 0.300,
     :          0.150, 0.150, 0.000, 0.284, 0.284, 0.000,
     :          0.417, 0.417, 0.000, 0.551, 0.551, 0.000,
     :          0.684, 0.684, 0.000, 0.817, 0.817, 0.000,
     :          0.951, 0.951, 0.000, 1.000, 1.000, 0.169,
     :          1.000, 1.000, 0.436, 1.000, 1.000, 0.702/
      data ((rgb_zj(i,j),i=1,3),j=191,200)/
     :          1.000, 0.702, 1.000, 1.000, 0.436, 1.000,
     :          1.000, 0.169, 1.000, 0.951, 0.000, 0.951,
     :          0.817, 0.000, 0.817, 0.684, 0.000, 0.684,
     :          0.551, 0.000, 0.551, 0.417, 0.000, 0.417,
     :          0.284, 0.000, 0.284, 0.150, 0.000, 0.150/

C     Index for color table 6:
      data ((rgb6(i,j),i=1,3),j=1,57)/
     : 1.000,1.000,1.000,0.000,0.000,0.000,0.200,0.200,0.200,
     : 0.350,0.350,0.350,0.500,0.500,0.500,0.650,0.650,0.650,
     : 0.800,0.800,0.800,1.000,1.000,0.000,0.000,0.553,1.000,
     : 0.000,0.753,1.000,0.000,1.000,1.000,0.631,1.000,1.000,
     : 0.100,0.800,0.100,0.200,0.900,0.200,0.500,1.000,0.200,
     : 0.500,0.900,0.500,1.000,1.000,0.000,1.000,0.804,0.000,
     : 1.000,0.604,0.000,1.000,0.400,0.000,1.000,0.000,0.000,
     : 0.900,0.300,0.400,0.800,0.000,0.900,0.900,0.000,1.000,
     : 0.800,0.700,0.800,0.700,0.600,0.700,0.600,0.500,0.600,
     : 0.500,0.400,0.500,0.000,0.925,0.925,0.004,0.627,0.961,
     : 0.000,0.000,0.965,0.000,1.000,0.000,0.000,0.784,0.000,
     : 0.000,0.565,0.000,1.000,1.000,0.000,0.906,0.753,0.000,
     : 1.000,0.565,0.000,1.000,0.000,0.000,0.839,0.000,0.000,
     : 0.753,0.000,0.000,1.000,0.000,1.000,0.600,0.333,0.788,
     : 0.900,0.900,0.900,1.000,0.000,1.000,0.749,0.000,1.000,
     : 0.498,0.000,1.000,0.000,0.000,1.000,0.000,0.349,1.000,
     : 0.000,0.549,1.000,0.000,0.749,1.000,0.000,1.000,1.000,
     : 0.000,0.902,0.800,0.000,0.800,0.498,0.000,0.702,0.000,
     : 0.498,0.800,0.000,0.800,0.902,0.000,1.000,1.000,0.000/
      data ((rgb6(i,j),i=1,3),j=58,114)/
     : 1.000,0.800,0.000,1.000,0.600,0.000,1.000,0.400,0.000,
     : 1.000,0.000,0.000,0.800,0.000,0.000,0.600,0.000,0.000,
     : 0.400,0.000,0.000,0.400,0.000,0.400,0.600,0.000,0.600,
     : 0.800,0.000,0.800,1.000,0.000,1.000,0.749,0.000,1.000,
     : 0.498,0.000,1.000,0.050,0.050,0.050,0.100,0.100,0.100,
     : 0.900,0.900,0.900,0.950,0.950,0.950,1.000,0.000,0.800,
     : 1.000,0.000,1.000,0.800,0.000,1.000,0.600,0.000,1.000,
     : 0.400,0.000,1.000,0.000,0.000,1.000,0.000,0.400,1.000,
     : 0.000,0.600,1.000,0.000,0.800,1.000,0.000,1.000,1.000,
     : 0.000,1.000,0.800,0.000,1.000,0.600,0.000,1.000,0.000,
     : 0.600,1.000,0.000,0.800,1.000,0.000,1.000,1.000,0.000,
     : 1.000,0.800,0.000,1.000,0.600,0.000,1.000,0.400,0.000,
     : 1.000,0.000,0.000,1.000,0.000,0.400,1.000,0.000,0.600,
     : 1.000,0.000,0.800,1.000,0.000,1.000,0.800,0.000,1.000,
     : 0.600,0.000,1.000,0.400,0.000,1.000,0.000,0.000,1.000,
     : 0.000,0.400,1.000,0.000,0.600,1.000,0.000,0.800,1.000,
     : 0.800,0.600,0.400,1.000,0.800,0.600,1.000,1.000,0.600,
     : 0.600,1.000,0.400,0.400,1.000,0.600,0.000,1.000,0.000,
     : 0.200,0.800,0.600,0.200,0.600,0.400,0.100,0.400,0.200/
      data ((rgb6(i,j),i=1,3),j=115,139)/
     : 1.000,0.804,0.000,1.000,0.604,0.000,1.000,0.400,0.000,
     : 0.000,1.000,1.000,0.000,0.500,1.000,0.000,0.000,1.000,
     : 0.400,0.400,1.000,0.600,0.600,1.000,0.800,0.800,1.000,
     : 1.000,1.000,1.000,1.000,0.800,0.800,1.000,0.600,0.600,
     : 1.000,0.400,0.400,1.000,0.000,0.000,1.000,0.500,0.000,
     : 1.000,1.000,0.000,1.000,0.878,0.584,0.898,0.780,0.427,
     : 0.761,0.643,0.271,0.624,0.545,0.075,0.486,0.369,0.000,
     : 1.000,0.804,0.000,1.000,0.604,0.000,1.000,0.400,0.000,
     : 1.000,0.000,0.000/
C
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
C     Beginning of executable code...
C
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
      icoltable = col_tab

      IF (col_tab .eq.0 ) THEN   ! all lines are black

        DO i=1,3
          rgb_table(i,1) = 1.
        END DO
        DO j=2,nc_max
          DO i=1,3
            rgb_table(i,j) = 0.
          END DO
        END DO
        ncunique = nc_max

      ELSE IF (col_tab .eq. 1 ) THEN

        ncunique = 25
        DO j=1,nc_max
          jj=j
          IF (j.ge.4)  jj = mod(j-4,ncunique-3)+4
          DO i=1,3
            rgb_table(i,j) = rgbv(i,jj)
          END DO
        END DO

      ELSE IF (col_tab .eq.2 ) THEN

        ncunique = 25
        DO j=1,nc_max
          jj=j
          IF (j.ge.4) jj = mod(24-mod(j-1,22),22)+4
          DO i=1,3
            rgb_table(i,j) = rgbv(i,jj)
          END DO
        END DO

      ELSE IF( col_tab.eq.3) THEN

        ncunique = 18
        DO j=1,nc_max
          jj=j
          IF (j.ge.4) jj = mod(j-4,ncunique-3)+4
          DO i=1,3
            rgb_table(i,j) = rgbf(i,jj)
          END DO
        END DO

      ELSE IF( col_tab.eq.4) THEN    ! gray shade

        DO j=1,3
          DO i=1,3
            rgbg(i,j)=rgbv(i,j)
          END DO
        END DO

c use logarithmic

        rgbg(1,4)=0.995
        rgbg(2,4)=rgbg(1,4)
        rgbg(3,4)=rgbg(1,4)
        DO j=5,17
          rgbg(1,j) = log10((15.- real(j) + 4.)/15.*10.)
          rgbg(2,j)=rgbg(1,j)
          rgbg(3,j)=rgbg(1,j)
        END DO
        rgbg(1,18)=0.0
        rgbg(2,18)=rgbg(1,18)
        rgbg(3,18)=rgbg(1,18)
c
        ncunique = 18
        DO j=1,nc_max
          jj=j
          IF (j.ge.4) jj = mod(j-4,ncunique-3)+4
          DO i=1,3
            rgb_table(i,j) = rgbg(i,jj)
          END DO
        END DO


      ELSE IF( col_tab.eq.5 ) THEN ! define 200-element color table

        DO j=1,nc_max
          jj=j
          IF (j.ge.4)  jj = mod(j-4,197)+4
          DO i=1,3
            rgb_table(i,j)=rgb_zj(i,jj)
          END DO
        END DO
        ncunique = 200

      ELSE IF( col_tab.eq.6 ) THEN ! define 139-element color table

        ncunique = 139
        DO j=1,nc_max
          jj=j
          IF (j.ge.4)  jj = mod(j-4,ncunique-3)+4
          DO i=1,3
            rgb_table(i,j)=rgb6(i,jj)
          END DO
        END DO

      ELSE IF( col_tab.eq.7) THEN    ! gray shade

        DO j=1,3
          DO i=1,3
            rgbg(i,j)=rgbv(i,j)
          END DO
        END DO

c use logarithmic

        rgbg(1,18)=0.995
        rgbg(2,18)=rgbg(1,18)
        rgbg(3,18)=rgbg(1,18)
        DO jj=5,17
          j = 17 - jj + 5
          rgbg(1,j) = log10((15.- real(jj) + 4.)/15.*10.)
          rgbg(2,j)=rgbg(1,j)
          rgbg(3,j)=rgbg(1,j)
        END DO
        rgbg(1,4)=0.0
        rgbg(2,4)=rgbg(1,4)
        rgbg(3,4)=rgbg(1,4)
c
        ncunique = 18
        DO j=1,nc_max
          jj=j
          IF (j.ge.4) jj = mod(j-4,ncunique-3)+4
          DO i=1,3
            rgb_table(i,j) = rgbg(i,jj)
          END DO
        END DO

      ELSE IF (col_tab.eq.-1) THEN  ! user specifies own color table

        lenfil = max(1, index(coltabfn, ' ')-1 )

        inquire(file=coltabfn(1:lenfil),exist=fexist)
        IF(.not.fexist) THEN
          write(6,'(1x,a,a,a/1x,a/1x,a,a)')
     :    'Color table file ',coltabfn(1:lenfil),' does not exist.',
     :    'Please respecify the file name (this file is required when',
     :    'color table option -1 is chosen. ',
     :    'The default is zx_color.tbl).'
           STOP 101
        ENDIF

        open(1,file=coltabfn(1:lenfil),form='formatted',status='old')
        sum_colors=0
        DO j=1,nc_max
          read(1,*,end=100) (rgb_table(i,j),i=1,3)
          sum_colors=sum_colors+1
        END DO
100     CONTINUE
        IF(sum_colors.lt.nc_max ) THEN
          DO j=sum_colors+1,nc_max
            jj=mod(j, sum_colors)
            IF(jj.eq.0) jj = sum_colors
            DO i=1,3
              rgb_table(i,j) = rgb_table(i,jj)
            END DO
          END DO
        ENDIF
        ncunique = sum_colors

        CLOSE(1)

      END IF
c
c#######################################################################
c
c     Setup the color index (color table)
c
c#######################################################################
c
      CALL XWRTCTBL(rgb_table,nc_max,exchange)

      RETURN

      ENTRY XSTCTFN(col_tab_fn)
c
c#######################################################################
c
c
c     PURPOSE:
c
c     To be called before SETCOLORS to reset the default filename
c     (zx_color.tbl) for the case of user-specified color table.
c
c#######################################################################
c

      coltabfn = col_tab_fn

      RETURN

      ENTRY XQCLRTBL( ncoltable )

c#######################################################################
c
c     Return current color table number
c
c#######################################################################

      ncoltable = icoltable

      RETURN
      END SUBROUTINE XSETCLRS_new

      SUBROUTINE SETCOLORS(col_tab)
        integer col_tab
        call xsetclrs(col_tab,0)
      RETURN
      END SUBROUTINE SETCOLORS

      SUBROUTINE XPSPAC( PL0,PR0,PB0,PT0)
C Define the individual picture plotting area in ND-space
C Its arguments should be in the range of (0.0,1.5,0.0,1.0)
C By default this area covers the whole ND-space .
C All transformations in vector space are reset as default (cancelled)
C when XPSPAC is called.
      COMMON /XPHY01/PL,PR,PB,PT,XRANGE,YRANGE
      COMMON /XMAP04/ X1,X2,Y1,Y2,XSCALE,YSCALE
      COMMON /XFTR06/ XFACTR,YFACTR
      COMMON /XPHO03/ DXPO,DYPO
      COMMON /XMAO05/ DXMOP,DYMOP
      COMMON /XPRF07/XPREF,YPREF

      XRANGE0=PR-PL
      YRANGE0=PT-PB

      PL=PL0
      PR=PR0
      PB=PB0
      PT=PT0
      XRANGE=PR-PL
      YRANGE=PT-PB
      XFACTR=XRANGE/XSCALE
      YFACTR=YRANGE/YSCALE
      CALL XPSCOF
      CALL XMREFP(X1,Y1)
      CALL XUNMLC
      CALL XMROFF
      CALL XSROFF
      CALL XOBOFF

      if( abs(xrange-xrange0).gt.0.001.or.
     :    abs(yrange-yrange0).gt.0.001) then
        call xqthik(lthick)
        call xthick(lthick)
        call xbrokn0
      endif

      RETURN
      END SUBROUTINE XPSPAC

      subroutine xstchn(nch0)
      common /xoutch/ nch
      nch = nch0
      return
      end subroutine xstchn

      SUBROUTINE XQPSPC(PX1,PX2,PY1, PY2)
C Return the current picture space parameters defined by XPSPAC
C subject to no picture scaling.
      COMMON /XPHY01/PL,PR,PB,PT,XRANGE,YRANGE
      PX1=PL
      PX2=PR
      PY1=PB
      PY2=PT
      RETURN

      ENTRY XQRANG( RANGEX,RANGEY)
C Return the actual length of picture sides measured in ND-space
C subject to picture scaling. ( X and Y denote direction of axes
C before coordinate rotation but subject to overall picture ratation.)
      RANGEX=XRANGE
      RANGEY=YRANGE
      RETURN
      END SUBROUTINE XQPSPC

      SUBROUTINE XSPACE(NUMPH,NUMPV,ROTANG,XLIM, YLIM)
!C
!C  SUBROUTINE TO SET A GRAPHIC SPACE CONTAINING  NUMPH*NUMPV
!C  PICTURES IN ONE FRAME OF FILM BY MOVING AND ROTATING COORDINATES.
!C  INPUT : NUMPH,NUMPV- NUMBER OF PICTURES IN HORIZONATL AND VERTICAL
!C                       IN EACH FRAME
!C          ROTANG- THE ANGLE ATHAT EACH PICTURE IS ROTATED THROUGH
!C  OUTPUT: XLIMIT,YLIMIT--
!C          DEFINE THE MAXIMUM PLOTTING AREA FOR EACH PICTURE
!C          (-XLIMIT/2,XLIMIT/2,-YLIMIT/2,YLIMIT/2),
!C  ENTRIES: XNWPIC, XNWFRM, XPMAGN
!C
!C Option to switch off annotation for certain sub-pictures are included
!C This is controled by Entry XFAUTO for automatic frame setting.
      SAVE NCALLS, NOPIC , KFAUTO
      SAVE NCOUNT,NUMPX,NUMPY,NUMPIC,XLIMIT,YLIMIT,XMAGIN,YMAGIN,PANGLE
      DATA NCOUNT,NUMPX,NUMPY,NUMPIC,XLIMIT,YLIMIT,XMAGIN,YMAGIN,PANGLE
     ;     /  0,    1,    1,    1,     1.5,   1.0,  0.0,   0.0 ,   0.0/
      DATA KFAUTO /0/

      integer icplswitch
      common /xcplswitch/ icplswitch

      COMMON /XPSD01/ XSIDE, YSIDE
      COMMON /XAXS18/ KANX,KANY, KTKX,KTKY
      common /xoutch/ nch
      DATA NCALLS /0/


      NUMPX=NUMPH
      NUMPY=NUMPV
      NUMPIC=NUMPX*NUMPY
      NCOUNT=0
      XLIMIT=XSIDE/NUMPX
      YLIMIT=YSIDE/NUMPY
      PANGLE=ROTANG
      IF(PANGLE.EQ.90.0) THEN
      XLIM=YLIMIT
      YLIM=XLIMIT
      ELSE
      XLIM=XLIMIT
      YLIM=YLIMIT
      ENDIF
      RETURN

      ENTRY XNWPIC
C Used in relating to XSPACE to define the picture plotting space
C for next picture.
      NCALLS=NCALLS+1
      NCOUNT=NCOUNT+1
      IORIGN= MOD( NCOUNT, NUMPIC)
      IF((IORIGN.EQ. 1.OR.NUMPIC.EQ.1).AND.NCALLS.GT.1) CALL XFRAME
      NOPIC=IORIGN
      IF(IORIGN.EQ.0) NOPIC=NUMPIC
      WRITE(NCH,*) 'Picture No. ', NOPIC,' in the frame.'

      IF(PANGLE.EQ.90.0) THEN
        XOR=XLIMIT*(INT((NOPIC-1)/NUMPY)+0.5)
        YOR=YLIMIT*(MOD(NOPIC-1, NUMPY)+0.5)
      ELSE
        XOR=XLIMIT*(MOD(NOPIC-1,NUMPX)+0.5)
        YOR=YLIMIT*(INT(NUMPY-(NOPIC-1)/NUMPX)-0.5)
      ENDIF
      XRANGE=XLIMIT-2*XMAGIN
      YRANGE=YLIMIT-2*YMAGIN
      XC=XRANGE/2
      YC=YRANGE/2
C
      IF(PANGLE.EQ.90.0) THEN
        CALL XPSPAC( XOR-YC,XOR+YC,YOR-XC,YOR+XC)
      ELSE
        CALL XPSPAC( XOR-XC,XOR+XC,YOR-YC,YOR+YC)
      ENDIF
C
      IF( PANGLE.NE.0.0)  THEN
        CALL XDREFP( XOR,YOR)
        CALL XDRANG( PANGLE)
      ENDIF
      PPANG=PANGLE

      IF( KFAUTO.EQ.0) RETURN
       CALL XAXANT(-1,-1)
       IF(PANGLE.NE.90.0) THEN
         NX=NUMPX
         NY=NUMPY
       ELSE
         NX=NUMPY
         NY=NUMPX
       ENDIF
       NSEQH=MOD(NOPIC,NX)
       IF( NSEQH.EQ.0) NSEQH=NX
       IF( NSEQH.NE.1) KANY=0
       IF( NOPIC.LE. (NY-1)*NX ) KANX=0

       icplswitch = 1
       if( MOD(NOPIC,NX).ne.0) icplswitch = 0

       print*,'NOPIC,NX, MOD(NOPIC,NX),icplswitch=',
     : NOPIC,NX, MOD(NOPIC,NX),icplswitch

      RETURN

      ENTRY XPMAGN( XM,YM)
C Used in XSPACE to set the margins of graghic in the picture space
C provided. If not called ,default values of zero are provided.
      XMAGIN=XM
      YMAGIN=YM
      RETURN

      ENTRY XNWFRM
C Used in XSPACE to terminate the current picture frame and move on
C to the next page
      IFRAME=1
      NCOUNT=0
      RETURN

      ENTRY XQNPIC(NPIC)
      NPIC=NOPIC
      RETURN

      ENTRY XFAUTO(KFAU)
C* ADDED IN ZXPLOTI *
      KFAUTO=KFAU
      RETURN

      ENTRY XQSPAC(NPICH, NPICV, RANGLE, XLIM, YLIM)
c
      NPICH = NUMPX
      NPICV = NUMPY
      RANGLE = PANGLE
      XLIM = XLIMIT
      YLIM = YLIMIT
      RETURN

      END SUBROUTINE XSPACE

      SUBROUTINE XDRSET(X1,Y1)
C Perform rotation around device reference point.
      COMMON /XPHY01/PL,PR,PB,PT,XRANGE,YRANGE
      COMMON /XPHO03/ DXPO,DYPO
      COMMON /XPRF07/ XPREF,YPREF
      COMMON /XAGS09/ DRANG,CRANG,XANGLE,XSYMAN,SRANG,KSR,XA,YA
      COMMON /XSCS10/ SINDRA,COSDRA,SINMRA,COSMRA,SINSRA,COSSRA
     :              ,SINXA,COSXA,SINYA,COSYA ,CHSIN,CHCOS

      IF( DRANG.EQ.0.0) RETURN
      X2=X1-XPREF
      Y2=Y1-YPREF
      X1=X2*COSDRA-Y2*SINDRA+XPREF
      Y1=X2*SINDRA+Y2*COSDRA+YPREF
      RETURN

      ENTRY XDREFP(XP,YP)
C Define the device reference point in ND-space for overall picture
C rotation.
      XPREF=XP
      YPREF=YP
      RETURN

      ENTRY XQDREF(XP1,YP1)
      XP1=XPREF
      YP1=YPREF
      RETURN

      ENTRY XDRANG(ANG)
C Set the angle the overall picture is rotated through around the
C device reference point. (Defined by XDREFP )
      DRANG=ANG
      XANGLE=CRANG+DRANG+(90-SRANG)*KSR+XA
      RADANG= ATAN(1.)/45.0*ANG
      SINDRA= SIN( RADANG)
      COSDRA= COS( RADANG)
      RETURN

      ENTRY XDROFF
C Switch off rotation around the device reference point.
      DRANG=0.0
      XANGLE=CRANG+DRANG+(90-SRANG)*KSR+XA
      SINDRA=0.0
      COSDRA=1.0
      RETURN

      ENTRY XQDRAG(ANG1)
      ANG1=DRANG
      RETURN

      ENTRY XDLOCA(XPLOC,YPLOC)
C Define the position in ND-space to which the device reference point
C is moved. ( Picture translation in ND-space )
      DXPO=XPLOC-XPREF
      DYPO=YPLOC-YPREF
      RETURN

      ENTRY XUNDLC
C Cancel picture translation in ND-space.
      DXPO=0.0
      DYPO=0.0
      RETURN

      ENTRY XQDLOC(XPLOC1,YPLOC1)
      XPLOC1=XPREF+DXPO
      YPLOC1=YPREF+DYPO
      RETURN
      END

      SUBROUTINE XMAP(XL,XR,YB,YT)
C Map the picture space. (Define maths coordinates on the picture space)
C Transformations in vector space are reset as default when remaped.
      COMMON /XPHY01/PL,PR,PB,PT,XRANGE,YRANGE
      COMMON /XMAP04/ X1,X2,Y1,Y2,XSCALE,YSCALE
      COMMON /XFTR06/ XFACTR,YFACTR

      XSCALE=XR-XL
      YSCALE=YT-YB
      X1=XL
      X2=XR
      Y1=YB
      Y2=YT
      XFACTR=XRANGE/XSCALE
      YFACTR=YRANGE/YSCALE
      CALL XMREFP(X1,Y1)
      CALL XPSCOF
      CALL XUNMLC
      CALL XMROFF
      CALL XSROFF
      CALL XOBOFF
      RETURN

      ENTRY XQMAP (XL0,XR0,YB0,YT0)
C Return the range of the current mapping space
      XL0=X1
      XR0=X2
      YB0=Y1
      YT0=Y2
      RETURN
      END

      SUBROUTINE XMREFP(XREF,YREF)
C Define the picture reference point in mapped vector space as the
C center of picture scaling, deformation, rotation, and translation.
C All these transformations are cancelled when either XPSPAC or XMAP
C is called  (But the transformations defined in ND-space remain in
C effect).
      COMMON /XPHY01/PL,PR,PB,PT,XRANGE,YRANGE
      COMMON /XMAP04/ X1,X2,Y1,Y2,XSCALE,YSCALE
      COMMON /XMAO05/ DXMOP,DYMOP
      COMMON /XFTR06/ XFACTR,YFACTR
      COMMON /XMRF08/ XMREF,YMREF,XMPREF,YMPREF
      XMREF=XREF
      YMREF=YREF
      XMPREF=PL+(XMREF-X1)*(PR-PL)/(X2-X1)
      YMPREF=PB+(YMREF-Y1)*(PT-PB)/(Y2-Y1)
      RETURN

      ENTRY XQMREF(XREF1,YREF1)
      XREF1=XMREF
      YREF1=YMREF
      RETURN

      ENTRY XPSCAL(SCALEX,SCALEY)
C Define scaling factors
      XRANGE=(PR-PL)*SCALEX
      YRANGE=(PT-PB)*SCALEY
      XFACTR=XRANGE/XSCALE
      YFACTR=YRANGE/YSCALE
      RETURN

      ENTRY XPSCOF
C Switch off scaling
      XRANGE= PR-PL
      YRANGE= PT-PB
      XFACTR=XRANGE/XSCALE
      YFACTR=YRANGE/YSCALE
      RETURN

      ENTRY XQPSCL(SCALX1,SCALY1)
      SCALX1=XRANGE/(PR-PL)
      SCALY1=YRANGE/(PT-PB)
      RETURN

      ENTRY XMLOCA(XLOC,YLOC)
C Translate the picture reference point to the location (XLOC,YLOC)
C defined in mapped vetor space.
      DXMOP=(XLOC-XMREF) *XRANGE/XSCALE
      DYMOP=(YLOC-YMREF) *YRANGE/YSCALE
      RETURN

      ENTRY XUNMLC
C Switch off the translation in mapped vector space
      DXMOP=0.0
      DYMOP=0.0
      RETURN

      ENTRY XQMLOC(XLOC1,YLOC1)
      XLOC1=XMREF+DXMOP*XSCALE/XRANGE
      YLOC1=YMREF+DYMOP*YSCALE/YRANGE
      RETURN
      END

      SUBROUTINE XMRSET(X2,Y2)
C Perform rotation around picture reference point (XMREF,YMREF)
      COMMON /XMRF08/ XMREF,YMREF,XMPREF,YMPREF
      COMMON /XAGS09/ DRANG,CRANG,XANGLE,XSYMAN,SRANG,KSR,XA,YA
      COMMON /XSCS10/ SINDRA,COSDRA,SINMRA,COSMRA,SINSRA,COSSRA
     :              ,SINXA,COSXA,SINYA,COSYA ,CHSIN,CHCOS
      IF( CRANG.EQ.0.0) RETURN
      X1=X2-XMPREF
      Y1=Y2-YMPREF
      X2=X1*COSMRA-Y1*SINMRA +XMPREF
      Y2=X1*SINMRA+Y1*COSMRA +YMPREF
      RETURN

      ENTRY XMRANG(ANG)
C Set coordinate rotation angle (It supercedes the previous value.)
      CRANG=ANG
      RADANG= ATAN(1.)/45.0*ANG
      SINMRA= SIN( RADANG)
      COSMRA= COS( RADANG)
      XANGLE=CRANG+DRANG+(90-SRANG)*KSR+XA
      RETURN

      ENTRY XMROFF
C Turn off coordinate rotation
      CRANG=0.0
      XANGLE=CRANG+DRANG+(90-SRANG)*KSR+XA
      SINMRA=0.0
      COSMRA=1.0
      RETURN

      ENTRY XQMRAG(ANG1)
      ANG1=CRANG
      RETURN

      ENTRY XSRSET(X2,Y2)
C Perform shearing of the picture in x or y direction
      IF( SRANG.EQ.90.0) RETURN
        IF( KSR.EQ.0) THEN
          Y1=Y2-YMPREF
          X2=X2+Y1*COSSRA
          Y2=YMPREF+Y1*SINSRA
        ELSE
          X1=X2-XMPREF
          Y2=Y2+X1*COSSRA
          X2=XMPREF+X1*SINSRA
        ENDIF
      RETURN

      ENTRY XSHEAR( XYANGL, KSHEAR)
C Shear the picture in x or y direction.
C XYANGL  The angle between x and y-axis after shearing (default 90.0)
C KSHEAR  =0 when X-axis to be is fixed, =1 when Y-axis is to be fixed
      SRANG=XYANGL
      KSR=KSHEAR
      RAD=ATAN(1.0)/45.0*SRANG
      SINSRA=SIN(RAD)
      COSSRA=COS(RAD)
      XANGLE=CRANG+DRANG+(90-SRANG)*KSR+XA
      RETURN

      ENTRY XSROFF
C Switch off picture shearing
      SRANG=90.0
      KSR=0
      SINSRA=1.0
      COSSRA=0.0
      XANGLE=CRANG+DRANG+(90-SRANG)*KSR+XA
      RETURN
      END

      SUBROUTINE XOBSET(X2,Y2)
C Perform non-orthogonal rotation of coordinate system. (Deformation)
      COMMON /XMRF08/ XMREF,YMREF,XMPREF,YMPREF
      COMMON /XAGS09/ DRANG,CRANG,XANGLE,XSYMAN,SRANG,KSR,XA,YA
      COMMON /XSCS10/ SINDRA,COSDRA,SINMRA,COSMRA,SINSRA,COSSRA
     :              ,SINXA,COSXA,SINYA,COSYA ,CHSIN,CHCOS
      IF( XA.EQ.0.0.AND.YA.EQ.0.0) RETURN
      X1=X2-XMPREF
      Y1=Y2-YMPREF
      X2=X1*COSXA-Y1*SINYA+XMPREF
      Y2=X1*SINXA+Y1*COSYA+YMPREF
      RETURN

      ENTRY XOBANG( XANG, YANG)
C Define angles for non-orthogonal coordiante rotation.
C XANG -- The angle x-axis is rotated through relative to old x-axis.
C YANG -- The angle y-axis is rotated through relative to old y-axis.
      XA=XANG
      YA=YANG
      DR =ATAN(1.0)/45.0
      SINXA=SIN(XA*DR)
      COSXA=COS(XA*DR)
      SINYA=SIN(YA*DR)
      COSYA=COS(YA*DR)
      XANGLE=CRANG+DRANG+(90-SRANG)*KSR +XA
      RETURN

      ENTRY XOBOFF
C Switch off non-orthogonal rotation
      XA=0.0
      YA=0.0
      SINXA=0.0
      COSXA=1.0
      SINYA=0.0
      COSYA=1.0
      XANGLE=CRANG+DRANG+(90-SRANG)*KSR+XA
      RETURN

      ENTRY XQOBAG( XANG1, YANG1)
      XANG1=XA
      YANG1=YA
      RETURN
      END

      FUNCTION XLTRNX(X)
C Perform linear transformation in x-direction
      COMMON /XFTR06/ XFACTR,YFACTR
      COMMON /XMRF08/ XMREF,YMREF,XMPREF,YMPREF
      XLTRNX=(X-XMREF)*XFACTR+XMPREF
      RETURN
      END

      FUNCTION XLTRNY(Y)
C Perform linear transformation in y-direction
      COMMON /XFTR06/ XFACTR,YFACTR
      COMMON /XMRF08/ XMREF,YMREF,XMPREF,YMPREF
      XLTRNY=(Y-YMREF)*YFACTR+YMPREF
      RETURN
      END

      SUBROUTINE XLINVT(X,Y)
C Perform inverse linear transformation (from ND-space to mathe space).
      COMMON /XFTR06/ XFACTR,YFACTR
      COMMON /XMRF08/ XMREF,YMREF,XMPREF,YMPREF
      X=XMREF+(X-XMPREF)/XFACTR
      Y=YMREF+(Y-YMPREF)/YFACTR
      RETURN
      END

      SUBROUTINE XTRANS(X,Y)
C    TRANSFORM POINT (X,Y)  FROM MATHEMATICAL SPACE BACK
C    TO ABSOLUTE PICTURE PLOTTING SPACE
      COMMON /XPHO03/ DXPO,DYPO
      COMMON /XMAO05/ DXMOP,DYMOP
      COMMON /XPRJ26/ KPROJC
      EXTERNAL XLTRNX,XLTRNY
      X1=X
      Y1=Y
      IF(KPROJC.NE.0)  CALL XPROJC(X1,Y1)
      X1=XLTRNX(X1)
      Y1=XLTRNY(Y1)
C     CALL XSRSET(X1,Y1)
      CALL XOBSET(X1,Y1)
      CALL XMRSET(X1,Y1)
      X1=X1+DXMOP
      Y1=Y1+DYMOP
      CALL XDRSET(X1,Y1)
      X=X1+DXPO
      Y=Y1+DYPO
      RETURN
      END

      SUBROUTINE XPROJC(X,Y)
C A dummy routine which can used to define a projection (transformation)
C by user.  XPRJON should called when projection is to be switched on.
      RETURN
      END

      SUBROUTINE XPRJON
      COMMON /XPRJ26/ KPROJC
C Switch on user defined projection through XPROJC.
      KPROJC=1
      RETURN

      ENTRY      XPRJOF
C Switch off user defined projection through XPROJC.
      KPROJC=0
      RETURN
      END

      FUNCTION XPNTSD(X1,Y1,X2,Y2)
C Measure the distance in ND-space between two points
C  (X1,Y1) and (X2,Y2) defined in maths space
      PX1=X1
      PY1=Y1
      CALL XTRANS(PX1,PY1)
      PX2=X2
      PY2=Y2
      CALL XTRANS(PX2,PY2)

      XPNTSD=SQRT((PX2-PX1)*(PX2-PX1)+(PY2-PY1)*(PY2-PY1))
      RETURN
      END

      subroutine xmask(x1,x2,y1,y2)
      common /xmask1/ xm1(99),xm2(99),ym1(99),ym2(99),rmangl(99),
     :                cosmsa(99),sinmsa(99),lvlmsk
      common /xoutch/ nch
      if(lvlmsk.ge.98)then
        write(nch,'(a)')
     :  'Warning: level of masking exceeded 99, it was set to 99.'
      endif
      lvlmsk=min(99,lvlmsk+1)
      xm1(lvlmsk)=x1
      xm2(lvlmsk)=x2
      ym1(lvlmsk)=y1
      ym2(lvlmsk)=y2
      rmangl(lvlmsk)=0.0
      cosmsa(lvlmsk)=1.0
      sinmsa(lvlmsk)=0.0
      return

      entry xqlmsk(level)
      level=lvlmsk
      return
      end

      subroutine xrmask(x0,y0,xl,yl,rangle)
      common /xmask1/ xm1(99),xm2(99),ym1(99),ym2(99),rmangl(99),
     :                cosmsa(99),sinmsa(99),lvlmsk
      common /xoutch/ nch
      if(lvlmsk.ge.98)then
        write(nch,'(a)')
     :  'Warning: level of masking exceeded 99, it was set to 99.'
      endif
      lvlmsk=min(99,lvlmsk+1)
      pi=4.0*atan(1.0)
      rmangl(lvlmsk)=rangle
      cosmsa(lvlmsk)=cos(pi*rangle/180.0)
      sinmsa(lvlmsk)=sin(pi*rangle/180.0)
      xm1(lvlmsk)= x0*cosmsa(lvlmsk)+y0*sinmsa(lvlmsk)
      xm2(lvlmsk)= xm1(lvlmsk)+xl
      ym1(lvlmsk)=-x0*sinmsa(lvlmsk)+y0*cosmsa(lvlmsk)
      ym2(lvlmsk)= ym1(lvlmsk)+yl
      return
      end

      subroutine xmsprj(x,y,level)
      common /xmask1/ xm1(99),xm2(99),ym1(99),ym2(99),rmangl(99),
     :                cosmsa(99),sinmsa(99),lvlmsk
      lvl=level
      x0=x
      y0=y
      x= x0*cosmsa(lvl)+y0*sinmsa(lvl)
      y=-x0*sinmsa(lvl)+y0*cosmsa(lvl)
      return
      end

      subroutine xmsrpr(x,y,level)
      common /xmask1/ xm1(99),xm2(99),ym1(99),ym2(99),rmangl(99),
     :                cosmsa(99),sinmsa(99),lvlmsk
      lvl=level
      x0=x
      y0=y
      x= x0*cosmsa(lvl)-y0*sinmsa(lvl)
      y= x0*sinmsa(lvl)+y0*cosmsa(lvl)
      return
      end

      subroutine xunmsk( level )
      common /xmask1/ xm1(99),xm2(99),ym1(99),ym2(99),rmangl(99),
     :                cosmsa(99),sinmsa(99),lvlmsk
      lvlmsk=max(0,level-1)
      return
      end

      subroutine xrsmsk( level )
      common /xmask1/ xm1(99),xm2(99),ym1(99),ym2(99),rmangl(99),
     :                cosmsa(99),sinmsa(99),lvlmsk
      lvlmsk=level
      return
      end

      subroutine xtsmsk(x1,y1,x2,y2, lnsegs)
      common /xmask1/ xm1(99),xm2(99),ym1(99),ym2(99),rmangl(99),
     :                cosmsa(99),sinmsa(99),lvlmsk
      common /xoutch/ nch
      real x1(199),y1(199),x2(199),y2(199), idispl(199)
      dimension ic1(4),ic2(4)
      logical xinbdr

      if(lvlmsk.eq.0) then
        lnsegs=1
        return
      endif
      do 100 lv=1,lvlmsk

      lines=lnsegs
      do 200 ln=1,lines
      idispl(ln)=0
      xln1 = x1(ln)
      yln1 = y1(ln)
      xln2 = x2(ln)
      yln2 = y2(ln)
      if(rmangl(lv).ne.0.0) then
        call xmsprj(x1(ln),y1(ln),lv)
        call xmsprj(x2(ln),y2(ln),lv)
      endif
      call xmscod(x1(ln),y1(ln),lv,ic1)
      call xmscod(x2(ln),y2(ln),lv,ic2)

      isum1=ic1(1)+ic1(2)+ic1(3)+ic1(4)
      isum2=ic2(1)+ic2(2)+ic2(3)+ic2(4)

      if(isum1+isum2.eq.0) then  ! both ends are inside
        idispl(ln)=0
        x1(ln) = xln1
        y1(ln) = yln1
        x2(ln) = xln2
        y2(ln) = yln2
        goto 200
      endif

      do 20 i=1,4
      if(ic1(i)+ic2(i).eq.2) then   ! the line is obviously out side
        idispl(ln)=1
        x1(ln) = xln1
        y1(ln) = yln1
        x2(ln) = xln2
        y2(ln) = yln2
        goto 200
      endif
20    continue

      if(isum1.eq.0.or.isum2.eq.0) then ! one end is inside

      isw=0
      if(isum1.eq.0)then
        ic01=ic1(1)
        ic02=ic1(2)
        ic03=ic1(3)
        ic04=ic1(4)
        do 30 i=1,4
30      ic1(i)=ic2(i)
        ic2(1)=ic01
        ic2(2)=ic02
        ic2(3)=ic03
        ic2(4)=ic04
        x0=x1(ln)
        y0=y1(ln)
        x1(ln)=x2(ln)
        y1(ln)=y2(ln)
        x2(ln)=x0
        y2(ln)=y0
        isw=1
      endif

      knt = 0
      if(ic1(1).eq.1)then
       y0=y1(ln)+(xm1(lv)-x1(ln))*(y2(ln)-y1(ln))
     :                               /(x2(ln)-x1(ln))
       x0=xm1(lv)
       goto 160
      elseif(ic1(2).eq.1)then
       y0=y1(ln)+(xm2(lv)-x1(ln))*(y2(ln)-y1(ln))
     :                               /(x2(ln)-x1(ln))
       x0=xm2(lv)
       goto 160
      endif
150   if(ic1(3).eq.1)then
       x0=x1(ln)+(ym1(lv)-y1(ln))*(x2(ln)-x1(ln))
     :                               /(y2(ln)-y1(ln))
       y0=ym1(lv)
      elseif(ic1(4).eq.1)then
       x0=x1(ln)+(ym2(lv)-y1(ln))*(x2(ln)-x1(ln))
     :                               /(y2(ln)-y1(ln))
       y0=ym2(lv)
      endif
160   continue
      if(.not.xinbdr(x0,y0,xm1(lv),xm2(lv),ym1(lv),ym2(lv)))then
      knt=knt+1
      if(knt.gt.10)then
        WRITE(NCH,*)'Dead loop encountered in XTSMSK, job stopped.'
        stop 991
      endif
      goto150
      endif
      if(rmangl(lv).ne.0.0) call xmsrpr(x0,y0,lv)
      x2(ln)=x0
      y2(ln)=y0
      if(isw.eq.1)then
        x1(ln) = xln2
        y1(ln) = yln2
      else
        x1(ln) = xln1
        y1(ln) = yln1
      endif
      if(isw.eq.1)then
       x0=x1(ln)
       y0=y1(ln)
       x1(ln)=x2(ln)
       y1(ln)=y2(ln)
       x2(ln)=x0
       y2(ln)=y0
       isum2=0
      endif
      idispl(ln)=1

      else                           ! both ends are outside

      xa=x1(ln)
      ya=y1(ln)
      kount=0
      if(ic1(1).eq.1)then
       yb=y1(ln)+(xm1(lv)-x1(ln))*(y2(ln)-y1(ln))
     :                               /(x2(ln)-x1(ln))
       xb=xm1(lv)
       goto 250
      elseif(ic1(2).eq.1)then
       yb=y1(ln)+(xm2(lv)-x1(ln))*(y2(ln)-y1(ln))
     :                               /(x2(ln)-x1(ln))
       xb=xm2(lv)
       goto 250
      endif
260   if(ic1(3).eq.1)then
       xb=x1(ln)+(ym1(lv)-y1(ln))*(x2(ln)-x1(ln))
     :                               /(y2(ln)-y1(ln))
       yb=ym1(lv)
      elseif(ic1(4).eq.1)then
       xb=x1(ln)+(ym2(lv)-y1(ln))*(x2(ln)-x1(ln))
     :                               /(y2(ln)-y1(ln))
       yb=ym2(lv)
      endif
250   continue
      kount = kount+1
      if(kount.gt.10)then
        WRITE(NCH,*)'Dead loop encountered in XTSMSK, job stopped.'
        stop 992
      endif

      if(.not.xinbdr(xb,yb,xm1(lv),xm2(lv),ym1(lv),ym2(lv))
     :   .and.kount.eq.1)   goto260
      if(.not.xinbdr(xb,yb,xm1(lv),xm2(lv),ym1(lv),ym2(lv)))then
        idispl(ln)=1
        x1(ln) = xln1
        y1(ln) = yln1
        x2(ln) = xln2
        y2(ln) = yln2
        goto 200
      else
        lnsegs=lnsegs+1
        x1(lnsegs)=xln1
        y1(lnsegs)=yln1
        xb0=xb
        yb0=yb
        if(rmangl(lv).ne.0.0) call xmsrpr(xb0,yb0,lv)
        x2(lnsegs)=xb0
        y2(lnsegs)=yb0
        idispl(lnsegs)=1
        x1(ln)=xb
        y1(ln)=yb
      endif

      kount=0
      if(ic2(1).eq.1)then
        y0=y1(ln)+(xm1(lv)-x1(ln))*(y2(ln)-y1(ln))
     :                               /(x2(ln)-x1(ln))
        x0=xm1(lv)
        goto 360
      elseif(ic2(2).eq.1)then
        y0=y1(ln)+(xm2(lv)-x1(ln))*(y2(ln)-y1(ln))
     :                               /(x2(ln)-x1(ln))
        x0=xm2(lv)
        goto 360
      endif
350   continue
      if(ic2(3).eq.1)then
        x0=x1(ln)+(ym1(lv)-y1(ln))*(x2(ln)-x1(ln))
     :                               /(y2(ln)-y1(ln))
        y0=ym1(lv)
      elseif(ic2(4).eq.1)then
        x0=x1(ln)+(ym2(lv)-y1(ln))*(x2(ln)-x1(ln))
     :                               /(y2(ln)-y1(ln))
        y0=ym2(lv)
      endif
360   continuE
      kount=kount+1
      if(kount.gt.10)then
        WRITE(NCH,*)'Dead loop encountered in XTSMSK, job stopped.'
        stop 993
      endif
      if(.not.xinbdr(x0,y0,xm1(lv),xm2(lv),ym1(lv),ym2(lv))
     : .and. kount.eq.1)goto350
      if(rmangl(lv).ne.0.0) call xmsrpr(x0,y0,lv)
      x1(ln)=x0
      y1(ln)=y0
      if(rmangl(lv).ne.0.0) call xmsrpr(x2(ln),y2(ln),lv)
      idispl(ln)=1

      endif
200   continue
      lin = 0
      do 400 i=1,lnsegs
        if(idispl(i).eq.1)then
        lin=lin+1
        x1(lin)=x1(i)
        x2(lin)=x2(i)
        y1(lin)=y1(i)
        y2(lin)=y2(i)
        endif
400   continue
      lnsegs=lin
100   continue
      return
      end

      subroutine xmscod(x,y,level,ic)
      common /xmask1/ xm1(99),xm2(99),ym1(99),ym2(99),rmangl(99),
     :                cosmsa(99),sinmsa(99),lvlmsk
      integer ic(4)
      do 10 i=1,4
10      ic(i)=0
      if(x.le.xm1(level)) ic(1)=1
      if(x.ge.xm2(level)) ic(2)=1
      if(y.le.ym1(level)) ic(3)=1
      if(y.ge.ym2(level)) ic(4)=1
      return
      end

      logical function xinbdr(x,y,xw1,xw2,yw1,yw2)
c check if the point is within the border
      xinbdr=.false.
      if(x.ge.xw1.and.x.le.xw2.and.y.ge.yw1.and.y.le.yw2)xinbdr=.true.
      return
      end

      subroutine xtstwd(x1,y1,x2,y2,idispl)
      common /xwndw1/ xw1,xw2,yw1,yw2, iwndon
      dimension ic1(4),ic2(4)
      common /xoutch/ nch
      knt = 0
5     knt = knt+1
      call xencod(x1,y1,ic1)
      call xencod(x2,y2,ic2)

      isum1=ic1(1)+ic1(2)+ic1(3)+ic1(4)
      isum2=ic2(1)+ic2(2)+ic2(3)+ic2(4)
      idispl=1
      if(isum1+isum2.eq.0) goto 999

      idispl=0
      do 20 i=1,4
20    if(ic1(i)+ic2(i).eq.2) goto 999
c
c make sure (x1,y1) is outside the window
      isw=0
      if(isum1.eq.0)then
        ic01=ic1(1)
        ic02=ic1(2)
        ic03=ic1(3)
        ic04=ic1(4)
        do 30 i=1,4
30      ic1(i)=ic2(i)
        ic2(1)=ic01
        ic2(2)=ic02
        ic2(3)=ic03
        ic2(4)=ic04
        x0=x1
        y0=y1
        x1=x2
        y1=y2
        x2=x0
        y2=y0
        isw=1
      endif

      if(ic1(1).eq.1)then
       y1=y1+(xw1-x1)*(y2-y1)/(x2-x1)
       x1=xw1
      elseif(ic1(2).eq.1)then
       y1=y1+(xw2-x1)*(y2-y1)/(x2-x1)
       x1=xw2
      elseif(ic1(3).eq.1)then
       x1=x1+(yw1-y1)*(x2-x1)/(y2-y1)
       y1=yw1
      elseif(ic1(4).eq.1)then
       x1=x1+(yw2-y1)*(x2-x1)/(y2-y1)
       y1=yw2
      endif
      if(isw.eq.1)then
       x0=x1
       y0=y1
       x1=x2
       y1=y2
       x2=x0
       y2=y0
      endif
      idispl=1
      if(knt.gt.10)then
        WRITE(NCH,*)'Dead loop encountered in XTSTWD, job stopped.'
        stop 991
      endif
      goto 5
999   return
      end

      subroutine xencod(x,y,ic)
      common /xwndw1/ xw1,xw2,yw1,yw2, iwndon
      integer ic(4)
      do 10 i=1,4
10      ic(i)=0
      if(x.lt.xw1) ic(1)=1
      if(x.gt.xw2) ic(2)=1
      if(y.lt.yw1) ic(3)=1
      if(y.gt.yw2) ic(4)=1
      return
      end

      subroutine xwindw(x1,x2,y1,y2)
      common /xwndw1/ xw1,xw2,yw1,yw2, iwndon
      common /xcwndw/ icwndw, xcpen, ycpen
      xw1=x1
      xw2=x2
      yw1=y1
      yw2=y2
      iwndon=1
      icwndw=1
      return
      end

      subroutine xqwdwon(kwndon)
      common /xwndw1/ xw1,xw2,yw1,yw2, iwndon
      kwndon = iwndon
      return
      end

      subroutine xwdwof
      common /xwndw1/ xw1,xw2,yw1,yw2, iwndon
      common /xcwndw/ icwndw, xcpen, ycpen
      iwndon=0
      icwndw=0
      return
      end

      subroutine xqwndw(x1,x2,y1,y2)
      common /xwndw1/ xw1,xw2,yw1,yw2, iwndon
      common /xcwndw/ icwndw, xcpen, ycpen
      x1=xw1
      x2=xw2
      y1=yw1
      y2=yw2
      return
      end

      subroutine xtstchwrt(x,y,wrtch)
      common /xwndw1/ xw1,xw2,yw1,yw2, iwndon
      common /xmask1/ xm1(99),xm2(99),ym1(99),ym2(99),rmangl(99),
     :                cosmsa(99),sinmsa(99),lvlmsk
      integer wrtchw,wrtchm,wrtch

      wrtchw = 0
      if(iwndon.eq.0)then
        wrtchw=1
      elseif(x.ge.xw1.and.x.le.xw2.and.y.ge.yw1.and.y.le.yw2)then
        wrtchw=1
      endif

c     print*,'in xtstchwrt,  xw1,xw2,yw1,yw2,wrtchw=',
c    :            xw1,xw2,yw1,yw2,wrtchw

      wrtchm = 1
      if(lvlmsk.eq.0) then
        wrtchm = 1
      else
        do 100 lv=1,lvlmsk
          x1=x
          y1=y
          call xmsprj(x1,y1,lv)
          if(x1.ge.xm1(lv).and.x1.le.xm2(lv).and.
     :       y1.ge.ym1(lv).and.y1.le.ym2(lv)) then
            wrtchm = 0
          endif
100     continue
      endif

      wrtch=0
      if(wrtchw.eq.1.and.wrtchm.eq.1) wrtch=1

      return
      end

      SUBROUTINE xpenup(x,y)
C position pen at point (x,y) defined in maths space
      common /xpen11/ xpen,ypen,flen,blen,npd,xmpen,ympen
      xpp= x
      ypp= y
      xmpen= x
      ympen= y
      call xtrans(xpp,ypp)
      call xtpnup(xpp,ypp)
      flen=0.0
      blen=0.0
      npd=0
      return
      end

      subroutine xpendn(x,y)
C Join  point (x,y) defined in maths space
      common /xpen11/ xpen,ypen,flen,blen,npd,xmpen,ympen
      common /xlpn13/ hf1,hb1,hf2,hb2,lfull,lfull0,lthick, dthick
      common /xwndw1/ xw1,xw2,yw1,yw2, iwndon
      common /xmask1/ xm1(99),xm2(99),ym1(99),ym2(99),rmangl(99),
     :                cosmsa(99),sinmsa(99),lvlmsk
      real xa(199),ya(199),xb(199),yb(199)

      save hf,hb
      xp1=xpen
      yp1=ypen
      xp2 = x
      yp2 = y
      xmpen0=xmpen
      ympen0=ympen
      xmpen = x
      ympen = y

      x1=xmpen0
      y1=ympen0
      x2=xmpen
      y2=ympen
      if(iwndon.eq.1)then
      call xtstwd(x1,y1,x2,y2,idispl)
      if(idispl.eq.0)then
        CALL XTRANS(X2 ,Y2 )
        CALL XTPNup(X2 ,Y2 )
        return
      endif
      endif

      lnsegs=1
      xa(1)=x1
      ya(1)=y1
      xb(1)=x2
      yb(1)=y2
      if(lvlmsk.gt.0)then
      call xtsmsk(xa,ya,xb,yb,lnsegs)

      if(lnsegs.eq.0)then
        CALL XTRANS(X2 ,Y2 )
        CALL XTPNup(X2 ,Y2 )
        return
      endif
      endif

      do 100 lin=1,lnsegs
      x1=xa(lin)
      y1=ya(lin)
      xp2=xb(lin)
      yp2=yb(lin)
      x2=xp2
      y2=yp2

      if(x1.ne.xmpen0.or.y1.ne.ympen0.or.lin.ne.1)then
        xp1=x1
        yp1=y1
        CALL XTRANS(XP1 ,YP1 )
        CALL XTPNup(XP1 ,YP1 )
      endif

      CALL XTRANS(XP2 ,YP2 )
      IF(LFULL0.EQ.1) THEN
        CALL XTPNDN(XP2 ,YP2 )
        goto 15
      endif

      ZL=SQRT((XP2-XP1)*(XP2-XP1)+(YP2-YP1)*(YP2-YP1))
      IF( ZL.LT.1.0E-20  ) GO TO 16
      XR=(XP2-XP1)/ZL
      YR=(YP2-YP1)/ZL
      IF(MOD(NPD,2).EQ.0)THEN
        HF=HF1
        HB=HB1
      ELSE
        HF=HF2
        HB=HB2
      ENDIF
      IF(BLEN.NE.0.0 ) GOTO 28
 20   IF(ZL-(HF-FLEN)) 22,21,21
 21   XP1=XP1+(HF-FLEN)*XR
      YP1=YP1+(HF-FLEN)*YR
      IF( HF.LT.1.0E-10) THEN
        CALL XPPONT(XP1,YP1)
      ELSE
        CALL XTPNDN (XP1,YP1)
      ENDIF
      ZL=ZL-(HF-FLEN)
      FLEN=0.0
 28   IF(ZL-(HB-BLEN)) 26,25,25
 25   XP1=XP1+(HB-BLEN)*XR
      YP1=YP1+(HB-BLEN)*YR
      CALL XTPNUP(XP1,YP1)
      ZL=ZL-(HB-BLEN)
      BLEN=0.0
      NPD=NPD+1
      IF(MOD(NPD,2).EQ.0)THEN
        HF=HF1
        HB=HB1
      ELSE
        HF=HF2
        HB=HB2
      ENDIF
      GO TO 20
 22   FLEN=FLEN+ZL
      BLEN=0.0
      IF( HF.GE.1.0E-10)CALL XTPNDN(XP2,YP2)
      GO TO 15
 26   BLEN=BLEN+ZL
      FLEN=0.0
 16   CALL XTPNUP(XP2,YP2)
 15   CONTINUE
100   continue

      if(iwndon.eq.0.and.lvlmsk.le.0) return
      if(x2.ne.xmpen.or.y2.ne.ympen)then
        XP2 = Xmpen
        YP2 = Ympen
        CALL XTRANS(Xp2 ,Yp2 )
        CALL XTPNup(Xp2 ,Yp2 )
      endif
      RETURN


      entry xqmpen( xmp, ymp )
      xmp=xmpen
      ymp=ympen
      return

      entry xqppen( xpp, ypp )
      xpp=xpen
      ypp=ypen
      return
      end

      SUBROUTINE XLPNUP(XP,YP)
C Used in contouring routine in the place of XPENUP to incorperate
C contour labeling.
      COMMON /XFTR06/ XFACTR,YFACTR
      COMMON /XMAP04/ XL,XR,YB,YT,XSCALE,YSCALE
      COMMON /XLAB14/ DLABEL,WLABEL,HLABEL,SIZLB,KLBTYP,ICLI,ICLON
      CHARACTER CLABEL*20
      COMMON /CLABEL/XP1,YP1,XPP1,YPP1,DL,WL
      COMMON /XLAB15/CLABEL
      COMMON /XLAB16/ LCLAB
      COMMON /XLBA33/ LABROT
      common /xwndw1/ xw1,xw2,yw1,yw2, iwndon
      REAL XBUF(0:31), YBUF(0:31)
      SAVE NCALLS ,XLB1, YLB1 ,XP2,YP2 ,NCALDN,XBUF,YBUF,NBUF
      DATA NCALLS ,NCALDN /0, 0/, NBUF/0/

      REAL XLBP1, YLBP1
      SAVE XLBP1, YLBP1
      DATA XLBP1, YLBP1 /0, 0/

      integer labmask
      common /labmask1/ labmask

      real xe(4),ye(4)

      IF(NBUF.NE.0.AND. NCALDN.NE.0 ) THEN
        CALL XPENUP(XBUF(0), YBUF(0))
        DO 200 NN=0,NBUF
  200   CALL XPENDN(XBUF(NN),YBUF(NN))
        NBUF=0
        WL=0.0
      ENDIF

      NCALLS=NCALLS+1
      XP1=XP
      YP1=YP
      CALL XPENUP(XP1,YP1)
      IF( ICLI .EQ.0.OR.ICLON.EQ.0) RETURN
      XPP1=XLTRNX(XP)
      YPP1=XLTRNY(YP)
      DL=DLABEL*ABS(SIN(137.0*NCALLS))*0.5
      WL=0.0
      NBUF=0
      IF( KLBTYP ) 1,2,3
 1      SIZLB=0.02*(YT-YB)
 3      HLABEL=SIZLB*YFACTR
        WLABEL= HLABEL*LCLAB*0.77
 2    CONTINUE
      RETURN

      ENTRY XLPNDN (XP,YP)
C Used in contouring routine in the place of XPENDN to incorperate
C contour labeling.
      XP2=XP
      YP2=YP
      IF( ICLI .EQ.0.OR.ICLON.EQ.0) GOTO 9101
      XPP2=XLTRNX(XP)
      YPP2=XLTRNY(YP)
      D21=SQRT((XPP2-XPP1)*(XPP2-XPP1)+(YPP2-YPP1)*(YPP2-YPP1))
      IF(D21.LT.1.0E-10)GOTO 9101
      DL=DL+D21
      IF(DL-DLABEL) 9101,9102,9102
 9102 IF(WL.EQ.0.0) THEN
        DRATIO=(DLABEL-DL+D21)/D21
        XLB1=XP1+DRATIO*(XP2-XP1)
        YLB1=YP1+DRATIO*(YP2-YP1)
        XLBP1=XLTRNX(XLB1)
        YLBP1=XLTRNY(YLB1)
        CALL XPENDN( XLB1,YLB1 )
        NCALDN=1
        XP1=XLB1
        YP1=YLB1
        XPP1=XLBP1
        YPP1=YLBP1
        XBUF(0)=XLB1
        YBUF(0)=YLB1
        NBUF=0
      ENDIF
      D21=SQRT((XPP2-XPP1)*(XPP2-XPP1)+(YPP2-YPP1)*(YPP2-YPP1))
      WL=SQRT( (XPP2-XLBP1)*(XPP2-XLBP1)+(YPP2-YLBP1)*(YPP2-YLBP1))
      IF(WL    -WLABEL) 9111,9112,9112
 9112 DRATIO=(WLABEL-WL+D21) /D21
      XLB2=XP1+DRATIO*(XP2-XP1)
      YLB2=YP1+DRATIO*(YP2-YP1)
      XLBP2=XLTRNX(XLB2)
      YLBP2=XLTRNY(YLB2)
      IF( LABROT.EQ.0) THEN
        ANG=0.0
        DPX=0.0
        DPY=HLABEL*0.4
        XLB =(XLBP1+ XLBP2)*0.5
        YLB =(YLBP1+ YLBP2)*0.5-DPY
        GOTO 130
      ENDIF
      IF(ABS( XLBP2-XLBP1).LE. 1.0E-10) THEN
        ANG=90.0
      ELSE
        ANG=ATAN( (YLBP2-YLBP1)/(XLBP2-XLBP1) )*180/3.1415926535
      ENDIF
        DX=XLBP2-XLBP1
        DY=YLBP2-YLBP1
        IF( ABS(DX).LT.1.0E-6) THEN
          DPX=HLABEL*0.4
          DPY=0.0
        ELSE
          AK=DY/DX
          A=SQRT(1.0+AK*AK)
          DPX=HLABEL*AK/A  *0.4
          DPY=HLABEL/A     *0.4
        ENDIF
        XLB =(XLBP1+ XLBP2)*0.5  +ABS(DPX) *SIGN(1.0, ANG)
        YLB =(YLBP1+ YLBP2)*0.5  -ABS(DPY)
 130    CONTINUE

      CALL XLINVT( XLB, YLB )

      IF(iwndon.eq.0 .or. ((xlb-xw1)*(xlb-xw2).le.0.0 .and.
     :                     (ylb-yw1)*(ylb-yw2).le.0.0) ) THEN

C Draw boxes  around labels . Usful when wish to blank the labeled area.

        IF( labmask.ne.0) then

        call xqcolor(kcolor)
        call xcolor(0)

        DBX=DPX*1.5
        DBY=DPY*1.5
        XE(1)=XLBP1+DBX
        YE(1)=YLBP1-DBY
        XE(2)=XLBP1-DBX
        YE(2)=YLBP1+DBY
        XE(3)=XLBP2-DBX
        YE(3)=YLBP2+DBY
        XE(4)=XLBP2+DBX
        YE(4)=YLBP2-DBY
        CALL XLINVT( XE(1),YE(1))
        CALL XLINVT( XE(2),YE(2))
        CALL XLINVT( XE(3),YE(3))
        CALL XLINVT( XE(4),YE(4))
        CALL XFILAREA(xe,ye,4)

c       CALL XPENUP( XE1,YE1)
c       CALL XPENDN( XE2,YE2)
c       CALL XPENDN( XE3,YE3)
c       CALL XPENDN( XE4,YE4)
c       CALL XPENDN( XE1,YE1)
C
        call xcolor(kcolor)

        endif
c
C write the label
        CALL XQCHOR( ANGSYM )
        CALL XCHORI( ANG)
        CALL XQCHMG( CMAG )
        CALL XCHMAG( HLABEL)
        CALL XCHARC( XLB,YLB, CLABEL(1:LCLAB) )
        CALL XCHORI( ANGSYM )
        CALL XCHMAG( CMAG )
      ENDIF

      DL=MIN( WL    -WLABEL , DLABEL)
      WL=0.0
      NBUF=0
      CALL XPENUP( XLB2, YLB2)
      CALL XPENDN( XP2,YP2 )
      GOTO 9150
 9111 WL=WL
C     CALL XPENUP(XP2,YP2)
      NBUF=NBUF+1
      XBUF(NBUF)=XP2
      YBUF(NBUF)=YP2
      IF(NBUF.GE.30) THEN
         DO 210 NN=0,NBUF,2
         XBUF(NN/2)=XBUF(NN)
  210    YBUF(NN/2)=YBUF(NN)
         NBUF=NN/2
      ENDIF
      GOTO 9150
 9101 CALL XPENDN(XP2,YP2)
 9150 XP1=XP2
      YP1=YP2
      XPP1=XPP2
      YPP1=YPP2
      CONTINUE
      RETURN

      ENTRY XLABMASK(lbmsk)
      labmask = lbmsk

      RETURN
      END

      SUBROUTINE XLBINT( NCLI)
      CHARACTER CLABEL*20 ,CLABL*(*), LABEL*20, LBFORM*(*)
      COMMON /XLAB14/ DLABEL,WLABEL,HLABEL,SIZLB,KLBTYP,ICLI,ICLON
      COMMON /XLAB15/CLABEL
      COMMON /XLAB16/ LCLAB
      COMMON /XLBA33/ LABROT
        CALL XQRANG( XRG,YRG)
        DLABEL=XRG/NCLI
        ICLI=NCLI
      RETURN

      ENTRY XLBSIZ( H1 )
        KLBTYP=1
        SIZLB=abs(H1)
      RETURN

      ENTRY XLBMAG( H)
        KLBTYP=0
        HLABEL=abs(H)
        WLABEL= HLABEL*LCLAB*0.77
      RETURN

      ENTRY XLBROT(KROT)
        LABROT=KROT
      RETURN

      ENTRY XLBON
        ICLON=1
      RETURN

      ENTRY XLBOFF
        ICLON=0
      RETURN

      ENTRY XQLBON(KLBON)
        KLBON=ICLON
      RETURN

      ENTRY XLABEL( CLABL )
        CLABEL=CLABL
        LCLAB=   LEN ( CLABL )
        WLABEL=HLABEL*LCLAB*0.77
      RETURN

      ENTRY XQLABL( LABEL , LCH)
        LABEL= CLABEL
        LCH=LCLAB
      RETURN

      ENTRY XLBFM ( LBFORM )
      RETURN
      END

      SUBROUTINE XINUMB(X,Y,I, FORM )
      CHARACTER CH*132 , FORM*(*)
      IF( FORM.EQ. '*' ) THEN
        CALL XICH(I,CH,LCH)
      ELSE
        WRITE( CH, FORM )  I
        LCH= ICLENG ( CH )
      ENDIF
      CALL XCHARL(X,Y,CH(1:LCH) )
      END
      SUBROUTINE XRNUMB(X,Y,R, FORM )
      CHARACTER CH*132 , FORM*(*)
      IF( FORM.EQ. '*' ) THEN
        CALL XRCH(R,CH,LCH)
      ELSE
        WRITE( CH, FORM )  R
        LCH= ICLENG ( CH )
      ENDIF
      CALL XCHARL(X,Y,CH(1:LCH) )
      END
      SUBROUTINE XICH( I,CH, LCH)
      CHARACTER CH*20
      WRITE(CH,'( I20  )')  I
      LCH=20
      CALL XCHLJ( CH, LCH )
      END
      SUBROUTINE XRCH( R,CH,LCH)
C Return real number R as a character string in automatically set format
      CHARACTER CH*20
      ABSR=ABS(R)
        IF(ABSR.GE.1.0E5.OR.(ABSR.GT.0.0.AND.ABSR.LT.1.0E-2))THEN
          WRITE(CH,'(1P,E20.2)') R
        ELSEIF( ABSR.LT.0.1.AND. ABSR.NE.0.0) THEN
          WRITE(CH,'(F20.2)') R
        ELSE
          WRITE(CH,'(F20.1)') R
        ENDIF
      LCH=20
      CALL XCHLJ( CH, LCH)
      END
      SUBROUTINE XCHLJ( CH,LCH)
C Left justify a character string.
      CHARACTER CH*(*) , CH1*20
      K=1
      LCH=LEN( CH )
      DO 1 L=1,LCH
        IF( CH(L:L).NE.' ') THEN
          K=L
          GOTO 2
        ENDIF
 1    CONTINUE
 2    CH1=CH
      CH=' '
      CH(1:LCH-K+1)=CH1(K:LCH)
      LCH=LCH-K+1
      RETURN
      END


      SUBROUTINE XLETER(XO,YO,STRING, IPOS )
      COMMON /XCHP21/ XCHPEN, YCHPEN ,XCHMO,YCHMO,XCHPO,YCHPO
      COMMON /XPHY01/ PL,PR,PB,PT,XRANGE,YRANGE
      COMMON /XMAP04/ XL,XR,YB,YT,XSCALE,YSCALE
      COMMON /XFTR06/ XFACTR,YFACTR
      COMMON /XCHA20/ HCTR,SCTR,CRATIO, KFONT,NUNDLN
      COMMON /XASC12/ IASCII(300)
      COMMON /XCHR30/ ICRAM(256)
      COMMON /XCHR31/ CHDATA
      COMMON /XCHR32/ ICDATA
      INTEGER ICDATA (0:150, 32:127)
      CHARACTER     CTEMP*1, CH*5,CHDATA(127)*300
      CHARACTER*(*) STRING
      LOGICAL MODE
      common /xoutch/ nch

      CH=CHDATA(2)
      READ(CH  ,103)IXCHR,IYCHR
      CXY=0.75
      SY = HCTR /( YFACTR *IYCHR)
      SX = HCTR /( XFACTR *IYCHR) *CRATIO/CXY
      IF( KFONT.EQ.2) THEN
      FACTOR=6.0/4.2
      SX=SX*FACTOR
      SY=SY*FACTOR
      ENDIF

      XCHMO=XO
      YCHMO=YO
      XCHPO=XLTRNX(XO)
      YCHPO=XLTRNY(YO)
      N =   LEN (STRING)
      IF( IPOS.LT.0) GOTO 600
      ITX=0
      DO 8 ICHR=1,N
      CTEMP = STRING (ICHR:ICHR)
      I = ICHAR(CTEMP)
c     I = ICRAM(I)
c     I = IASCII(I)
      IF( I.EQ.0) I=32
      IF( ICDATA(0,I).NE.KFONT) THEN
        CALL XCHDEC(ICDATA,CHDATA,I)
        ICDATA(0,I)=KFONT
      ENDIF
      NCD=ICDATA(1,I)
      IX= ICDATA(NCD-1,I)
      IF( IX.GE.50)   IX=IX-50
      ITX=ITX+IX
  8   CONTINUE
      XWIDTH= ITX* SX
 600  IF( IPOS ) 601,602,603
 602     XSPOS=XO-0.5*XWIDTH
      GOTO 300
 603     XSPOS=XO- XWIDTH
      GOTO 300
 601     XSPOS=XO
 300  CONTINUE
      YSPOS = YO
      XSP = XSPOS
      YSP = YSPOS
      XTPOS = XSPOS
      YTPOS = YSPOS
      CALL XCPNUP(XSP, YSP)
      DO 1 ICHR=1,N
         CTEMP = STRING (ICHR:ICHR)
      I = ICHAR(CTEMP)
c     I = ICRAM(I)
c     I = IASCII(I)
      IF (I .EQ. 0) THEN
        I=32
        IF( CTEMP.NE.' ')
     :  WRITE(NCH,*)' Can not draw character ',CTEMP,' it was replaced'
     : ,' by a blank by ZXPLOT.'
      ENDIF

      IF( ICDATA(0,I).NE.KFONT) THEN
        CALL XCHDEC(ICDATA,CHDATA,I)
        ICDATA(0,I)=KFONT
      ENDIF
      NCD=ICDATA(1,I)

      DO 3 ICD=2,NCD,2
         IX= ICDATA(ICD,I)
         JY= ICDATA(ICD+1,I)
           MODE=.TRUE.
         IF( IX.GE.50) THEN
           MODE=.FALSE.
           IX=IX-50
         ENDIF
         XTPOS = XSP + FLOAT(IX)*SX
         YTPOS = YSP + FLOAT(JY)*SY
C        IF (XTPOS.GT.XR.OR.XTPOS.LT.XL) WRITE(NCH,*)'Out of bound in x-dir.'
C        IF (YTPOS.GT.YT.OR.YTPOS.LT.YB) WRITE(NCH,*)'Out of bound in y-dir.'
         IF (MODE)  THEN
            CALL XCPNDN(XTPOS, YTPOS )
         ELSE
            CALL XCPNUP(XTPOS, YTPOS )
         ENDIF
  3   CONTINUE
         XSP=XTPOS
         YSP=YTPOS
  1   CONTINUE

      DO 2 N=1,NUNDLN
         XST = XSPOS - 5.0*SX
         YST = YSPOS - 15.0*SY
         XFI = XTPOS + 5.0*SX
         YFI = YTPOS - 15.0*SY
         CALL XCPNUP (XST,YST)
         CALL XCPNDN (XFI,YFI)
    2 CONTINUE
C     XCHPEN=XTPOS
C     YCHPEN=YTPOS
      RETURN
  103 FORMAT(I2,1X,I2)
      END

      FUNCTION XCHLEN(STRING)
      COMMON /XPHY01/ PL,PR,PB,PT,XRANGE,YRANGE
      COMMON /XMAP04/ XL,XR,YB,YT,XSCALE,YSCALE
      COMMON /XFTR06/ XFACTR,YFACTR
      COMMON /XCHA20/ HCTR,SCTR,CRATIO, KFONT,NUNDLN
      COMMON /XASC12/ IASCII(300)
      COMMON /XCHR30/ ICRAM(256)
      COMMON /XCHR31/ CHDATA
      COMMON /XCHR32/ ICDATA
      INTEGER ICDATA (0:150, 32:127)
      CHARACTER     CTEMP*1, CH*5,CHDATA(127)*300
      CHARACTER*(*) STRING

      CH=CHDATA(2)
      READ(CH  ,103)IXCHR,IYCHR
      CXY=0.75
      CALL XQPSCL( XSC, YSC )
      SY = 1.0/ ( IYCHR* YFACTR)*YSC
      SX = 1.0/ ( IYCHR* XFACTR)*XSC      *CRATIO/CXY
      IF( KFONT.EQ.2) THEN
      FACTOR=6.0/4.2
      SX=SX*FACTOR
      SY=SY*FACTOR
      ENDIF
      N =   LEN (STRING)
      ITX=0
      DO 8 ICHR=1,N
      CTEMP = STRING (ICHR:ICHR)
      I = ICHAR(CTEMP)
c     I = ICRAM(I)
c     I = IASCII(I)
      IF( I.EQ.0) I=32
      IF( ICDATA(0,I).NE.KFONT) THEN
        CALL XCHDEC(ICDATA,CHDATA,I)
        ICDATA(0,I)=KFONT
      ENDIF
      NCD=ICDATA(1,I)
      IX= ICDATA(NCD-1,I)
      IF( IX.GE.50)   IX=IX-50
      ITX=ITX+IX
  8   CONTINUE
      XWIDTH= ITX* SX
      XP1= 0.0
      YP1= 0.0
      XP2= XWIDTH
      YP2= 0.0
      CALL XCTRAN(XP1,YP1)
      CALL XCTRAN(XP2,YP2)
      XCHLEN=SQRT( (XP2-XP1)*(XP2-XP1)+ (YP2-YP1)*(YP2-YP1))
      RETURN
  103 FORMAT(I2,1X,I2)
      END

      FUNCTION ICLENG( CH )
      CHARACTER*(*) CH
      ICLENG=0
      IC=LEN( CH )
      DO 5 L=1,IC
 5    IF( CH(L:L).NE.' ') ICLENG=L
      RETURN
      END

      SUBROUTINE XCHOBL( CTROBL )
      COMMON /XCHA20/  HCTR,SCTR,CRATIO, KFONT,NUNDLN
      CRATIO= CTROBL
      RETURN

      ENTRY XQCHOB( COBL )
      COBL=CRATIO
      RETURN
      END

      SUBROUTINE XCHLIN( N )
      COMMON /XCHA20/ HCTR,SCTR,CRATIO, KFONT,NUNDLN
      NUNDLN=N
      RETURN

      ENTRY XQCHLN ( NN )
      NN=NUNDLN
      RETURN
      END

      SUBROUTINE XSTRLNTH( string, length )
c     Return the length of the non-blank part of a character string.
c     INPUT:
c       string   A character string
c       length   The declared length of the character string 'string'.
c     OUTPUT:
c       length   The length of the non-blank part of the string.
c
      implicit none
      character string*(*)
      integer length
      integer i

      DO 100 i = length,1,-1
        IF(string(i:i) .ne. ' '.and.string(i:i).ne.' ') GOTO 200
100   continue
200   CONTINUE
      length = max(1,i)
      RETURN
      END

      SUBROUTINE XCHRST(X2,Y2)
C Perform rotation around picture reference point (XMREF,YMREF)
      COMMON /XCHP21/ XCHPEN, YCHPEN ,XCHMO,YCHMO,XCHPO,YCHPO
      COMMON /XAGS09/ DRANG,CRANG,XANGLE,XSYMAN,SRANG,KSR,XA,YA
      COMMON /XSCS10/ SINDRA,COSDRA,SINMRA,COSMRA,SINSRA,COSSRA
     :              ,SINXA,COSXA,SINYA,COSYA ,CHSIN,CHCOS
      IF( XSYMAN.EQ.0.0) RETURN
      X1=X2-XCHPO
      Y1=Y2-YCHPO
      X2=X1*CHCOS-Y1*CHSIN +XCHPO
      Y2=X1*CHSIN+Y1*CHCOS +YCHPO
      RETURN

      ENTRY XCHORI(CHANG)
      XSYMAN=CHANG
      IF( CHANG.EQ.0) GOTO 3
      RADANG= ATAN(1.)/45.0*XSYMAN
      CHSIN= SIN( RADANG)
      CHCOS= COS( RADANG)
      XANGLE=CRANG+DRANG+(90-SRANG)*KSR+XA
      RETURN

 3    CHSIN= 0.0
      CHCOS= 1.0
      RETURN

      ENTRY XQCHOR(SYMANG)
      SYMANG=XSYMAN
      RETURN
      END

      SUBROUTINE XCPNUP(X,Y)
      COMMON /XCHP21/ XCHPEN, YCHPEN ,XCHMO,YCHMO,XCHPO,YCHPO
      COMMON /XPEN11/ XPEN,YPEN,FLEN,BLEN,NPD,XMPEN,YMPEN
      common /xcwndw/ icwndw, xcpen, ycpen
      common /xmask1/ xm1(99),xm2(99),ym1(99),ym2(99),rmangl(99),
     :                cosmsa(99),sinmsa(99),lvlmsk
      XCHPEN=X
      YCHPEN=Y
      X1=X
      Y1=Y
      CALL XCTRAN( X1,Y1)
      xpen=x1
      ypen=y1
      CALL PPENUP( X1,Y1)
      if(icwndw.eq.1.or.lvlmsk.ge.1)then
      xcpen=xchpen
      ycpen=ychpen
      Xcpen=XLTRNX(Xcpen)
      Ycpen=XLTRNY(Ycpen)
      CALL XCHRST(Xcpen,ycpen)
      call xlinvt(xcpen,ycpen)
      XmPEN=xcpen
      YmPEN=ycpen
      endif
      RETURN
      END

      subroutine xcpndn(x,y)
      common /xchp21/ xchpen, ychpen ,xchmo,ychmo,xchpo,ychpo
      COMMON /XPEN11/ XPEN,YPEN,FLEN,BLEN,NPD,XMPEN,YMPEN
      common /xcwndw/ icwndw, xcpen, ycpen
      common /xmask1/ xm1(99),xm2(99),ym1(99),ym2(99),rmangl(99),
     :                cosmsa(99),sinmsa(99),lvlmsk
      real xa(199),ya(199),xb(199),yb(199)
      xchpen=x
      ychpen=y
      if(icwndw.eq.0.and.lvlmsk.eq.0)then
        x2=x
        y2=y
        call xctran(x2,y2)
        xpen=x2
        ypen=y2
        call ppendn(x2,y2)
        xcpen=x
        ycpen=y
        Xcpen=XLTRNX(Xcpen)
        Ycpen=XLTRNY(Ycpen)
        CALL XCHRST(Xcpen,ycpen)
        call xlinvt(xcpen,ycpen)
        XmPEN=xcpen
        YmPEN=ycpen
        goto 999
      endif

      xcpen0=xcpen
      ycpen0=ycpen
      x1=xcpen0
      y1=ycpen0
      X2=X
      Y2=Y
      X2=XLTRNX(X2)
      Y2=XLTRNY(Y2)
      CALL XCHRST(X2,Y2)
      call xlinvt(x2,y2)
      xcpen=x2
      ycpen=y2
      xmpen=xcpen
      ympen=ycpen

      if(icwndw.eq.1)then
      call xtstwd(x1,y1,x2,y2,idispl)
      x2a=x2
      y2a=y2
      if(idispl.ne.1)then
        call xtrans(x2,y2)
        xpen=x2
        ypen=y2
        call ppenup(x2,y2)
        goto 999
      endif
      endif

      lnsegs=1
      xa(1)=x1
      ya(1)=y1
      xb(1)=x2
      yb(1)=y2
      if(lvlmsk.gt.0)then
      call xtsmsk(xa,ya,xb,yb,lnsegs)

      if(lnsegs.eq.0)then
        call xtrans(x2 ,y2 )
        xpen=x2
        ypen=y2
        call ppenup(x2 ,y2 )
        goto 999
      endif
      endif

      do 100 lin=1,lnsegs
      x1=xa(lin)
      y1=ya(lin)
      x2=xb(lin)
      y2=yb(lin)
      x2a=x2
      y2a=y2
      if(x1.ne.xcpen0.or.y1.ne.ycpen0.or.lin.ne.1)then
        call xtrans(x1 ,y1 )
        xpen=x1
        ypen=y1
        call ppenup(x1 ,y1 )
      endif
      call xtrans(x2,y2)
      xpen=x2
      ypen=y2
      call ppendn(x2,y2)
100   continue

      if(x2a.ne.xcpen.or.y2a.ne.ycpen)then
        x2=xcpen
        y2=ycpen
        call xtrans(x2,y2)
        xpen=x2
        ypen=y2
        call ppenup(x2,y2)
      endif
999   continue
      RETURN

      ENTRY XQCPEN( XCHP, YCHP )
      XCHP=XCHPEN
      YCHP=YCHPEN
      RETURN
      END

      SUBROUTINE XCTRAN(X,Y)
      COMMON /XPHO03/ DXPO,DYPO
      COMMON /XMAO05/ DXMOP,DYMOP
      COMMON /XPRJ26/ KPROJC
      EXTERNAL XLTRNX,XLTRNY
      X1=X
      Y1=Y
      IF(KPROJC.NE.0)  CALL XPROJC(X1,Y1)
      X1=XLTRNX(X1)
      Y1=XLTRNY(Y1)
      CALL XCHRST(X1,Y1)
C     CALL XSRSET(X1,Y1)
      CALL XOBSET(X1,Y1)
      CALL XMRSET(X1,Y1)
      X1=X1+DXMOP
      Y1=Y1+DYMOP
      CALL XDRSET(X1,Y1)
      X=X1+DXPO
      Y=Y1+DYPO
      RETURN
      END

C UTILITY ROUTINES:
      SUBROUTINE XPOINT(X,Y)
C Plot a point at position (X,Y) of mathematical space with predefined
C size. ( The size can be defined by XRPONT).
      COMMON /XCIR25/ XCIR(9) ,YCIR(9) , RPOINT
      X1=X
      Y1=Y
      CALL XTRANS(X1,Y1)
      CALL PPENUP(X1+RPOINT*XCIR(1), Y1+RPOINT*YCIR(1) )
      DO 6 I=1,9
 6    CALL PPENDN(X1+RPOINT*XCIR(I), Y1+RPOINT*YCIR(I) )
      RETURN

      ENTRY XPPONT(XP,YP)
      CALL PPENUP(XP+RPOINT*XCIR(1), YP+RPOINT*YCIR(1) )
      DO 5 I=3,9,2
 5    CALL PPENDN(XP+RPOINT*XCIR(I), YP+RPOINT*YCIR(I) )
      RETURN

      ENTRY XPNTSZ(R)
C Define the size of points to be plotted by XPOINT by their radius
C in ND-space.  By default R=0.0005
      RPOINT=R
      RETURN
      END

      SUBROUTINE XINTMKR
c
c     This subroutine will define (initialize) some additional marker
c     shapes for zxplot Library.
c
c     Keith Brewster, CAPS, April, 2009
c     Modified to add marker types 5-10 for filled symbols
c
      parameter (mxmrkty=10, mxmrkp=10)
      COMMON /XPSD01/ XSIDE, YSIDE
      COMMON /XMRK25/ IMKRFIL(mxmrkty),
     :       XMRK(mxmrkp,mxmrkty),YMRK(mxmrkp,mxmrkty),
     :       MDX(mxmrkty),RMARKER

      PI = 4.0*ATAN(1.0)
      DO 72 J=1,mxmrkty
        DO 71 I=1,mxmrkp
          XMRK(I,J)=0.0
          YMRK(I,J)=0.0
71      CONTINUE
        MDX(J)=0
        IMKRFIL(J)=0
72    CONTINUE

      DO 80 I=1,9
        XMRK(I,1)=COS((I-1)*0.25*PI)
        YMRK(I,1)=SIN((I-1)*0.25*PI)
80    CONTINUE
      mdx(1)=9

      XMRK(1,2) = 1.0
      YMRK(1,2) = 1.0
      XMRK(2,2) = 0.
      YMRK(2,2) = -1.0
      XMRK(3,2) = -1.0
      YMRK(3,2) = 1.0
      XMRK(4,2) = 1.0
      YMRK(4,2) = 1.0
      MDX(2)=4
      XMRK(1,3) = 0.0
      YMRK(1,3) = 2.0
      XMRK(2,3) = 1.0
      YMRK(2,3) = 0.0
      XMRK(3,3) = -1.0
      YMRK(3,3) = 0.0
      XMRK(4,3) = 0.0
      YMRK(4,3) = 2.0
      MDX(3)=4
      XMRK(1,4) = 1.0
      YMRK(1,4) = 1.0
      XMRK(2,4) = 1.0
      YMRK(2,4) = -1.0
      XMRK(3,4) = -1.0
      YMRK(3,4) = -1.0
      XMRK(4,4) = -1.0
      YMRK(4,4) = 1.0
      XMRK(5,4) = 1.0
      YMRK(5,4) = 1.0
      MDX(4)=5
      XMRK(1,5) = 1.0
      YMRK(1,5) = 0.0
      XMRK(2,5) = 0.0
      YMRK(2,5) = -1.0
      XMRK(3,5) = -1.0
      YMRK(3,5) = 0.0
      XMRK(4,5) = 0.0
      YMRK(4,5) = 1.0
      XMRK(5,5) = 1.0
      YMRK(5,5) = 0.0
      MDX(5)=5

      DO 92 J=6,10
        DO 91 I=1,mxmrkp
          XMRK(I,J)=XMRK(I,J-5)
          YMRK(I,J)=YMRK(I,J-5)
91      CONTINUE
        MDX(J)=MDX(J-5)
        IMKRFIL(J)=1
92    CONTINUE

      RMARKER=0.0010*YSIDE

      RETURN
      END

      SUBROUTINE XMARKER(X,Y,ITY)
c
c     This subroutine will draw markers
c
c     x,y  marker''s coordination
c     ity  Marker number
c
c     Keith Brewster, CAPS, April, 2009
c     Added processing for filled symbols through IMKRFIL
c
      parameter (mxmrkty=10, mxmrkp=10)
      COMMON /XMRK25/ IMKRFIL(mxmrkty),
     :       XMRK(mxmrkp,mxmrkty),YMRK(mxmrkp,mxmrkty),
     :       MDX(mxmrkty),RMARKER
      common /xwndw1/ xw1,xw2,yw1,yw2, iwndon
      real Xw(mxmrkp), Yw(mxmrkp)

      IF(iwndon.ne.0) THEN
        IF((x-xw1)*(x-xw2).gt.0.0 .or.((y-yw1)*(y-yw2).gt.0.0))
     :  return
      ENDIF

      X1=X
      Y1=Y
      CALL XTRANS(X1,Y1)

      IF( IMKRFIL(ITY).eq.0 ) THEN
        CALL PPENUP(X1+RMARKER*XMRK(1,ITY),Y1+RMARKER*YMRK(1,ITY) )
        DO 7 I=1,MDX(ITY)
 7      CALL PPENDN(X1+RMARKER*XMRK(I,ITY), Y1+RMARKER*YMRK(I,ITY))
      ELSE
        CALL xafstyl(IMKRFIL(ITY))
        DO 9 I=1,MDX(ITY)
          Xw(I) = X1+RMARKER*XMRK(I,ITY)
          Yw(I) = Y1+RMARKER*YMRK(I,ITY)
          CALL XLINVT (Xw(I),Yw(I))
 9      CONTINUE
        CALL XFILAREA(Xw,Yw,MDX(ITY))
      ENDIF

      RETURN

      ENTRY XMRKSZ(R)
        RMARKER=R
      RETURN

      END

      SUBROUTINE XBOX(X1,X2,Y1,Y2)
      CALL XPENUP( X1,Y1)
      CALL XPENDN ( X2,Y1)
      CALL XPENDN ( X2,Y2)
      CALL XPENDN ( X1,Y2)
      CALL XPENDN ( X1,Y1)
      RETURN
      END

      SUBROUTINE XBORDR
C DRAW A BORDER AROUND MAPPED AERA
      COMMON /XMAP04/ X1,X2,Y1,Y2,XSCALE,YSCALE
      CALL XBOX(X1,X2,Y1,Y2)
      RETURN
      END

      SUBROUTINE XAXES(XO,XSTEP,YO,YSTEP)
C Draw X and Y axis through (XO,YO) with tick interval XSTEP and YSTEP
C If XSTEP or YSTEP=0.0,the intervals are set automatically
      CALL XAXISX1(XO,YO,XSTEP, 0.0)
      CALL XAXISY1(XO,YO,YSTEP, 0.0)
      RETURN
      END

      SUBROUTINE XAXISX(XO,YO,XSTEP)
      CALL XAXISX1(XO,YO,XSTEP, 0.0)
      RETURN
      END

      SUBROUTINE XAXISY(XO,YO,YSTEP)
      CALL XAXISY1(XO,YO,YSTEP, 0.0)
      RETURN
      END

      SUBROUTINE XAXISX1(XO,YO,XSTEP_in,XMJSTEP)
C  To draw X-AXIS through (XO,YO) with tick interval of XSTEP.
C  If XSTEP=0.0,the interval is set automatically
c 2/17/1999 (M.Xue) Wrote this XAXISX and added YMJSTEP.

      PARAMETER( JUMP =2 )
      COMMON /XMAP04/ XL,XR,YB,YT,XSCALE,YSCALE
      COMMON /XPHY01/ PL,PR,PB,PT,XRANGE,YRANGE
      COMMON /XFTR06/ XFACTR,YFACTR
      COMMON /XFMT33/ LBFMT, AXFMT
      COMMON /XFMT34/ LLBFMT,LAXFMT
      COMMON /XAXS18/ KANX,KANY, KTKX,KTKY
      CHARACTER LBFMT*50, AXFMT*10
      CHARACTER CH*20
      EXTERNAL XAXINC
      COMMON /XAXM35/ NTMAG, ANMAG, ANSIZ
      real xmjstep
      integer ijump, passed, ifold

      IF( KTKX.EQ.0 .and.kanx.eq.0) RETURN
      IF( xr.eq.xl ) return

      UNITH=SQRT( ABS(XRANGE*YRANGE))*0.01
      UH= MIN( ABS(XRANGE), ABS(YRANGE) )*0.03
      IF( NTMAG.EQ.0) ANMAG=UH
      IF( NTMAG.EQ.2) ANMAG=ANSIZ*YFACTR
      HX= ANMAG/XFACTR
      HY= ANMAG/YFACTR
      CALL XQCHMG( HOLD )
      CALL XCHMAG( ANMAG )

      CALL XPENUP( XL, YO)
      CALL XPENDN( XR, YO)

      xstep = xstep_in

      IF( XSTEP.EQ.0.0) THEN
        XSTEPJ=XAXINC(XR-XL)
        xstep = xstepj
      ELSEIF(xmjstep.eq.0.0) THEN
        XSTEPJ=XSTEP
 5      IF( abs(NINT((XR-XL)/XSTEPJ)).GT.6  ) THEN
          XSTEPJ=XSTEPJ*2
          GOTO 5
        ENDIF
      ELSE
        ijump=nint( xmjstep/xstep)
        XSTEPJ=XSTEP*ijump
      ENDIF

      IF( KTKX.EQ.0) GOTO 110
      HTICK=UNITH /YFACTR*KTKX
      AXL=XO+ NINT((XL-XO)/XSTEP)*XSTEP

      AXL=XO+ NINT((XL-XO)/XSTEPJ)*XSTEPJ - XSTEPJ

      tem=XSTEPJ/XSTEP

      IFOLD=NINT(XSTEPJ/XSTEP)
      IFOLD=NINT( tem )
      ifold = max(1, ifold)

      epsx = (xr-xl)*0.001

      passed =0
      DO 150 i=0,5000

      X=AXL+I*XSTEP

      IF ((X-(XL-epsx))*(X-(XR+epsx)).le.0.0) then
        passed=1

        IF( KTKX.NE.0) THEN
          XHT=HTICK
          IF( MOD(I, IFOLD).EQ.0) XHT=HTICK+HTICK
          CALL XPENUP(X,YO)
          CALL XPENDN(X,YO+XHT)
        ENDIF

        IF( KANX.NE.0 .and. MOD(I, IFOLD).EQ.0) THEN

          IF( AXFMT(1:LAXFMT).EQ.'*') THEN
            CALL XRCH(X, CH, LCH)
          ELSE
            DO 503 KCH=1,LAXFMT
              IF((AXFMT(KCH:KCH).EQ.'I')
     :             .or.(AXFMT(KCH:KCH).EQ.'i'))THEN
              WRITE(CH,AXFMT(1:LAXFMT)) NINT(X)
              GOTO 504
            ENDIF
 503        CONTINUE
              WRITE(CH,AXFMT(1:LAXFMT)) X
 504        LCH=ICLENG( CH )
            CALL XCHLJ( CH(1:LCH), LCH)
          ENDIF

          YSHIFT=HY*1.5*KANX
          CALL XCHARC(X, YO+YSHIFT, CH(1:LCH) )
        ENDIF

      ELSE
        IF( passed.eq.1) GOTO 110
      ENDIF

150   CONTINUE

110   CONTINUE

      CALL XCHMAG( HOLD )
      RETURN
      END

      SUBROUTINE XAXISY1(XO,YO,YSTEP_in,YMJSTEP)
C  To draw Y-AXIS through (XO,YO) with tick interval of YSTEP.
C  If YSTEP=0.0,the interval is set automatically
c 2/17/1999 (M.Xue) Wrote this XAXISY and added XMJSTEP.
      implicit none
      real xo,yo,ystep_in,ymjstep
      real XL,XR,YB,YT,XSCALE,YSCALE
      COMMON /XMAP04/ XL,XR,YB,YT,XSCALE,YSCALE
      real PL,PR,PB,PT,XRANGE,YRANGE
      COMMON /XPHY01/ PL,PR,PB,PT,XRANGE,YRANGE
      real XFACTR,YFACTR
      COMMON /XFTR06/ XFACTR,YFACTR
      CHARACTER LBFMT*50, AXFMT*10
      COMMON /XFMT33/ LBFMT, AXFMT
      integer LLBFMT,LAXFMT
      COMMON /XFMT34/ LLBFMT,LAXFMT
      integer KANX,KANY, KTKX,KTKY
      COMMON /XAXS18/ KANX,KANY, KTKX,KTKY
      CHARACTER CH*20
      real xaxinc, ystep
      EXTERNAL XAXINC
      real anmag,ansiz
      integer ntmag,jfold,lch,kch,icleng
      COMMON /XAXM35/ NTMAG, ANMAG, ANSIZ
      real y,eps,unith,uh,hx,hy,hold,ystepj,AYB,yht,htick
      integer j,jjump, passed

      IF( KTKY.EQ.0 .and.kanY.eq.0) RETURN
      IF( yt.eq.yb ) return

      UNITH=SQRT( ABS(XRANGE*YRANGE))*0.01
      UH= MIN( ABS(XRANGE), ABS(YRANGE) )*0.03
      IF( NTMAG.EQ.0) ANMAG=UH
      IF( NTMAG.EQ.2) ANMAG=ANSIZ*YFACTR
      HX= ANMAG/XFACTR
      HY= ANMAG/YFACTR
      CALL XQCHMG( HOLD )
      CALL XCHMAG( ANMAG )

      CALL XPENUP( XO, YB)
      CALL XPENDN( XO, YT)

      ystep = ystep_in

      IF( YSTEP.EQ.0.0) THEN
        YSTEPJ=XAXINC(YT-YB)
        ystep = ystepj
      ELSEIF(ymjstep.eq.0.0) THEN
        YSTEPJ=YSTEP
 5      IF( abs( NINT((YT-YB)/YSTEPJ)).GT.6  ) THEN
          YSTEPJ=YSTEPJ*2
          GOTO 5
        ENDIF
      ELSE
        jjump=nint( ymjstep/YSTEP)
        YSTEPJ=YSTEP*jjump
      ENDIF

      IF( KTKY.EQ.0) GOTO 110
      HTICK=UNITH/XFACTR*KTKY
      AYB=YO+ NINT((YB-YO)/YSTEP)*YSTEP

      AYB=YO+ NINT((YB-YO)/YSTEPJ)*YSTEPJ - YSTEPJ

      JFOLD=NINT(YSTEPJ/YSTEP)
      eps = (YT-YB)*0.001

      passed =0
      DO 150 j=0,5000

      Y=AYB+j*YSTEP
      IF ((Y-(YB-eps))*(Y-(YT+eps)).le.0.0) then
        passed = 1

        IF( KTKX.NE.0) THEN
          YHT=HTICK
          IF( MOD(j, JFOLD).EQ.0) YHT=HTICK+HTICK
          CALL XPENUP(XO,Y)
          CALL XPENDN(XO+YHT,Y)
        ENDIF

        IF( KANY.NE.0 .and. MOD(J, JFOLD).EQ.0) THEN

          IF( AXFMT(1:LAXFMT).EQ.'*') THEN
            CALL XRCH(Y, CH, LCH)
          ELSE
            DO 503 KCH=1,LAXFMT
              IF((AXFMT(KCH:KCH).EQ.'I')
     :             .or.(AXFMT(KCH:KCH).EQ.'i'))THEN
              WRITE(CH,AXFMT(1:LAXFMT)) NINT(Y)
              GOTO 504
            ENDIF
 503        CONTINUE
              WRITE(CH,AXFMT(1:LAXFMT)) Y
 504        LCH=ICLENG( CH )
            CALL XCHLJ( CH(1:LCH), LCH)
          ENDIF

          IF(KANY.EQ.-1)
     :      CALL XCHARR(XO-HX*0.7, Y-0.4*HY, CH(1:LCH) )
          IF(KANY.EQ. 1)
     :      CALL XCHARL(XO+HX*0.7, Y-0.4*HY, CH(1:LCH) )

        ENDIF

      ELSE
        IF(passed.eq.1) GOTO 110
      ENDIF
150   CONTINUE

110   CONTINUE

      CALL XCHMAG( HOLD )
      RETURN
      END
      SUBROUTINE XAXANT(KANTX,KANTY)
      COMMON /XAXS18/ KANX,KANY, KTKX,KTKY
      COMMON /XFMT33/ LBFMT, AXFMT
      COMMON /XFMT34/ LLBFMT,LAXFMT
      CHARACTER LBFMT*50, AXFMT*10,AXFM*(*)
      COMMON /XAXM35/ NTMAG, ANMAG, ANSIZ
C KANTX KANTY--  Axis annotation parameters.
C KANTX=1  annotation lacated above x-axis
C KANTX=-1 annotation lacated below x-axis
C KANTX=0  annotation on x-axis is suppressed
C KANTY=1  annotation lacated to the right of y-axis
C KANTY=-1 annotation lacated to the left  of y-axis
C KANTY=0  annotation on y-axis is suppressed
C Default: KANTX=-1, KANTY=-1
      KANX=KANTX
      KANY=KANTY
      RETURN

      ENTRY XAXTIK(KTIKX,KTIKY)
C KTIKX KTIKY--  Axis ticking parameters.
C KTIKX=1  ticking lacated above x-axis
C KTIKX=-1 ticking lacated below x-axis
C KTIKX=0  ticking on x-axis is suppressed
C KTIKY=1  ticking lacated to the right of y-axis
C KTIKY=-1 ticking lacated to the left  of y-axis
C KTIKY=0  ticking on y-axis is suppressed
C Default: KTIKX= 1, KTIKY= 1
      KTKX=KTIKX
      KTKY=KTIKY
      RETURN

      ENTRY XAXDEF
C* ZXPLOTI *
C To restore the default values of parameters for axis annotation
C and ticking.
      KTKX=1
      KTKY=1
      KANX=-1
      KANY=-1
      RETURN

      ENTRY XAXFMT( AXFM )
C* MODIFIED IN ZXPLOTI, INTEGER FORMAT ALLOWED. *
      LAXFMT=LEN(AXFM)
      AXFMT=AXFM
      RETURN
      END

      SUBROUTINE XXAXIS(XCOOR,XVALUE,N,YO)
C  To draw an X-AXIS through (XCOOR(1),YO) and tickmark the axis
C  at x=xcoor(i) with value xvalue(i) for i=1,n.
      real xcoor(n), xvalue(n)
      COMMON /XMAP04/ XL,XR,YB,YT,XSCALE,YSCALE
      COMMON /XPHY01/ PL,PR,PB,PT,XRANGE,YRANGE
      COMMON /XFTR06/ XFACTR,YFACTR
      COMMON /XFMT33/ LBFMT, AXFMT
      COMMON /XFMT34/ LLBFMT,LAXFMT
      COMMON /XAXS18/ KANX,KANY, KTKX,KTKY
      CHARACTER LBFMT*50, AXFMT*10
      CHARACTER CH*20
      COMMON /XAXM35/ NTMAG, ANMAG, ANSIZ
      UNITH=SQRT( ABS(XRANGE*YRANGE))*0.01
      UH= MIN( ABS(XRANGE), ABS(YRANGE) )*0.03
      IF( NTMAG.EQ.0) ANMAG=UH
      IF( NTMAG.EQ.2) ANMAG=ANSIZ*YFACTR
      HX= ANMAG/XFACTR
      HY= ANMAG/YFACTR
      CALL XQCHMG( HOLD )
      CALL XCHMAG( ANMAG )
      CALL XPENUP( Xcoor(1), YO)
      CALL XPENDN( Xcoor(n), YO)
      YSHIFT=HY*1.0*KANX
      HTICK=UNITH*1.5/YFACTR*KTKX
      DO 150 I=1,n
         X=xcoor(i)
         IF( AXFMT(1:LAXFMT).EQ.'*') THEN
           CALL XRCH(Xvalue(i), CH, LCH)
         ELSE
           DO 503 KCH=1,LAXFMT
             IF((AXFMT(KCH:KCH).EQ.'I')
     :            .or.(AXFMT(KCH:KCH).EQ.'i'))THEN
             WRITE(CH,AXFMT(1:LAXFMT)) NINT(Xvalue(i))
             GOTO 504
             ENDIF
 503       CONTINUE
             WRITE(CH,AXFMT(1:LAXFMT)) Xvalue(i)
 504       LCH=ICLENG( CH )
           CALL XCHLJ( CH(1:LCH), LCH)
         ENDIF

         IF(KANX.EQ. 1) CALL XCHARC(X, YO+0.5*HY, CH(1:LCH))
         IF(KANX.EQ.-1) CALL XCHARC(X, YO-1.0*HY, CH(1:LCH))

         IF(KTKX.NE.0)THEN
           CALL XPENUP(X,YO)
           CALL XPENDN(X,YO+HTICK)
         ENDIF
 150  CONTINUE
      CALL XCHMAG( HOLD )
      RETURN
      END

      SUBROUTINE XYAXIS(XO,YCOOR,YVALUE,N)
C  To draw Y-axis through (XO,YCOOR(1)) and tickmark at y=ycoord(j)
C  with value yvalue(j) for j=1,n.
      REAL YCOOR(N),YVALUE(N)
      COMMON /XMAP04/ XL,XR,YB,YT,XSCALE,YSCALE
      COMMON /XPHY01/PL,PR,PB,PT,XRANGE,YRANGE
      COMMON /XFTR06/ XFACTR,YFACTR
      COMMON /XAXS18/ KANX,KANY, KTKX,KTKY
      COMMON /XFMT33/ LBFMT, AXFMT
      COMMON /XFMT34/ LLBFMT,LAXFMT
      CHARACTER LBFMT*50, AXFMT*10
      CHARACTER CH*20
      COMMON /XAXM35/ NTMAG, ANMAG, ANSIZ
      UNITH=SQRT( ABS(XRANGE*YRANGE))*0.01
      UH= MIN( ABS(XRANGE), ABS(YRANGE) )*0.03
      IF( NTMAG.EQ.0) ANMAG=UH
      IF( NTMAG.EQ.2) ANMAG=ANSIZ*YFACTR
      HX= ANMAG/XFACTR
      HY= ANMAG/YFACTR
      CALL XQCHMG( HOLD )
      CALL XCHMAG( ANMAG )

      CALL XPENUP( XO, Ycoor(1))
      CALL XPENDN( XO, Ycoor(n))
      HTICK=UNITH*1.5/XFACTR*KTKY
      DO 250 J=1,n
         Y=ycoor(j)
         IF( AXFMT(1:LAXFMT).EQ.'*') THEN
           CALL XRCH(Yvalue(j), CH, LCH)
         ELSE
           DO 503 KCH=1,LAXFMT
           IF((AXFMT(KCH:KCH).EQ.'I')
     :       .or.(axfmt(kch:kch).eq.'i')) THEN
             WRITE(CH,AXFMT(1:LAXFMT)) NINT(Yvalue(j))
             GOTO 504
           ENDIF
 503       CONTINUE
             WRITE(CH,AXFMT(1:LAXFMT)) Yvalue(j)
 504       LCH=ICLENG( CH )
           CALL XCHLJ(CH(1:LCH), LCH)
         ENDIF
         IF( KANY ) 301,300,302
 301     CALL XCHARR(XO-HX*0.3,Y-0.2*HY, CH(1:LCH) )
         GOTO 300
 302     CALL XCHARL(XO+HX*0.3,Y-0.2*HY, CH(1:LCH) )
 300     IF(KTKY.NE.0)THEN
           CALL XPENUP(XO,Y)
           CALL XPENDN (XO+HTICK,Y)
         ENDIF
 250  CONTINUE
      CALL XCHMAG( HOLD )
      RETURN
      END

      REAL FUNCTION XAXINC(X)
c
C TO SET ANNOTATION INCREMENT (ANNOTATIONS >=4 AND =<16 FOR FOLD=1.0)
c
c     Corrected version. 4/20/1994 Ming Xue.
c
      integer D
      real xlog

      IF(x.eq.0.0) THEN
        xaxinc = 1.0
        return
      ENDIF

      xlog = log10(x)

      IPOWER=INT(xlog )
      if( xlog.lt.0.0 ) ipower = ipower-1

      D= INT(X/10.0**IPOWER)

      FOLD=1.0
      IF(D.GE.1.AND.D.LT.3) THEN
        XAXINC=2.0*10.0**(IPOWER-1)*FOLD
      ELSEIF(D.GE.3.AND.D.LT.7) THEN
        XAXINC=5.0*10.0**(IPOWER-1)*FOLD
      ELSEIF(D.GE.7.AND.D.LT.10) THEN
        XAXINC=1.0*10.0** IPOWER*FOLD
      ELSEIF( d.eq.0) then
        XAXINC=1.0*10.0** IPOWER*FOLD
      ENDIF

      IF(XAXINC .EQ.0.0) XAXINC=X*0.1

      RETURN
      END

      SUBROUTINE XAXSCA(XL,XR,XSTEP, YB,YT,YSTEP)
C To draw ticks on border defined by (xl,xr,yb,yt) and annotate the
C ticks.  Modifications: Options for annotation and ticking included
C just as those in axis plotting routines.

      CALL XAXSCA1(XL,XR,XSTEP,0.0,YB,YT,YSTEP,0.0)

      RETURN
      END

      SUBROUTINE XAXSCA1(XL,XR,XSTEP,XMJSTEP,YB,YT,YSTEP,YMJSTEP)
c
c add new variable XJUMP, YJUMP by Min , others almost like XAXSCA
C To draw ticks on border defined by (xl,xr,yb,yt) and annotate the
C ticks.  Modifications: Options for annotation and ticking included
C just as those in axis plotting routines.
c
c Changed made by Ming Xue, 2/16/1998
c
c MX: 2/5/1999.
c XMJSTEP and YMJSTEP are now used to define major tick mark steps.
c 2/17/1999 (M.Xue) Rewrote this subroutine based on origin XAXSCA.
c
      PARAMETER( JUMP=2 )
      COMMON /XPHY01/PL,PR,PB,PT,XRANGE,YRANGE
      COMMON /XFTR06/ XFACTR,YFACTR
      COMMON /XFMT33/ LBFMT, AXFMT
      COMMON /XFMT34/ LLBFMT,LAXFMT
      COMMON /XAXS18/ KANX,KANY, KTKX,KTKY
      CHARACTER LBFMT*50, AXFMT*10
      CHARACTER CH*20
      integer xjump, yjump, passed
      real xmjstep, ymjstep
      COMMON /XAXM35/ NTMAG, ANMAG, ANSIZ
      SAVE KOR, X0, Y0
      DATA KOR /0/
      real eps

      IF( XSTEP.EQ.0.0.OR. YSTEP.EQ.0.0) RETURN

      UNITH=SQRT(ABS( XRANGE*YRANGE))*0.01
      UH= MIN( ABS(XRANGE), ABS(YRANGE) )*0.03
      IF( NTMAG.EQ.0) ANMAG=UH
      IF( NTMAG.EQ.2) ANMAG=ANSIZ*YFACTR
      HX= ANMAG/XFACTR
      HY= ANMAG/YFACTR
      CALL XQCHMG( HOLD )
      CALL XCHMAG( ANMAG )

      IF(xmjstep.eq.0.0) THEN
        XSTEPJ=XSTEP*JUMP
 5      IF( NINT((XR-XL)/XSTEPJ).GT.6  ) THEN
          XSTEPJ=XSTEPJ*2
          GOTO 5
        ENDIF
      ELSE
        xjump=nint( xmjstep/xstep)
        XSTEPJ=XSTEP*xjump
      ENDIF

      IF(ymjstep.eq.0.0) THEN
        YSTEPJ=YSTEP*JUMP
 6      IF( NINT((YT-YB)/YSTEPJ).GT.6  ) THEN
          YSTEPJ=YSTEPJ*2
          GOTO 6
        ENDIF
      ELSE
        yjump=nint( ymjstep/ystep)
        YSTEPJ=YSTEP*yjump
      ENDIF

      IF( KOR.EQ.1) THEN
        XO=X0
        YO=Y0
      ELSE
        XO=XL
        YO=YB
      ENDIF

      CALL XBOX(XL,XR,YB,YT)

      IF( KANX.EQ.0 .AND. KTKX.EQ.0 ) GOTO 160

      HTICK=UNITH/YFACTR*KTKX
      AXL=XO+NINT((XL-XO)/XSTEP)*XSTEP
      AXL=XO+NINT((XL-XO)/XSTEPJ)*XSTEPJ - XSTEPJ
      IFOLD=NINT(XSTEPJ/XSTEP)

      eps = 0.001*(xr-xl)

      passed =0
      DO 100 i=0,5000

      X=AXL+I*XSTEP

      IF ((x-(xl-eps))*(x-(xr+eps)).le.0.0) then
        passed =1

        IF( KTKX.NE.0) THEN
          XHT = HTICK
          IF( MOD(I, IFOLD).EQ.0) XHT=HTICK+HTICK
          CALL XPENUP(X,YB)
          CALL XPENDN(X,YB+XHT)
          CALL XPENUP(X,YT)
          CALL XPENDN(X,YT-XHT)
        endif

        IF( KANX.NE.0 .and.  MOD(I,IFOLD).EQ.0) THEN
          IF( AXFMT(1:LAXFMT).EQ.'*') THEN
            CALL XRCH(X, CH, LCH)
          ELSE
            DO 501 KCH=1,LAXFMT
              IF((AXFMT(KCH:KCH).EQ.'I')
     :          .or.(axfmt(kch:kch).eq.'i')) THEN
                WRITE(CH,AXFMT(1:LAXFMT)) NINT(X)
                GOTO 502
              ENDIF
 501        CONTINUE
            WRITE(CH,AXFMT(1:LAXFMT)) X
 502        LCH=ICLENG( CH )
            CALL XCHLJ(CH(1:LCH), LCH)
          ENDIF
          IF(KANX.EQ. 1)CALL XCHARC(X,YT+1.5*HY,CH(1:LCH))
          IF(KANX.EQ.-1)CALL XCHARC(X,YB-1.5*HY,CH(1:LCH))
        ENDIF

      ELSE
        IF(passed.eq.1) GOTO 160
      ENDIF
100   CONTINUE

160   CONTINUE

      IF( KANY.EQ.0 .AND. KTKY.EQ.0 ) GOTO 260

      HTICK=UNITH/XFACTR*KTKY
      AYB=YO+NINT((YB-YO)/YSTEP)*YSTEP
      AYB=YO+NINT((YB-YO)/YSTEPJ)*YSTEPJ - YSTEPJ

      JFOLD = NINT(YSTEPJ/YSTEP)

      eps = 0.001*(YT-YB)

      passed=0
      DO 200 j=0,5000

      Y=AYB+J*YSTEP

      IF ((y-(yb-eps))*(y-(yt+eps)).le.0.0)then
        passed =1

        IF(KTKY.NE.0)THEN
          YHT = HTICK
          IF( MOD(J, JFOLD).EQ.0) YHT=HTICK+HTICK
          CALL XPENUP(XL,Y)
          CALL XPENDN(XL+YHT,Y)
          CALL XPENUP(XR,Y)
          CALL XPENDN(XR-YHT,Y)
        endif

        IF( KANY.NE.0 .and. MOD(J, JFOLD).EQ.0) THEN
          IF( AXFMT(1:LAXFMT).EQ.'*') THEN
            CALL XRCH(Y, CH, LCH)
          ELSE
            DO 503 KCH=1,LAXFMT
            IF((AXFMT(KCH:KCH).EQ.'I')
     :        .or.(axfmt(kch:kch).eq.'i')) THEN
              WRITE(CH,AXFMT(1:LAXFMT)) NINT(Y)
              GOTO 504
            ENDIF
 503        CONTINUE
            WRITE(CH,AXFMT(1:LAXFMT)) Y
 504        LCH=ICLENG( CH )
            CALL XCHLJ(CH(1:LCH), LCH)
          ENDIF
          IF(KANY.EQ.-1)
     :      CALL XCHARR(XL-HX*0.7, Y-0.4*HY, CH(1:LCH) )
          IF(KANY.EQ. 1)
     :      CALL XCHARL(XR+HX*0.7, Y-0.4*HY, CH(1:LCH) )
        ENDIF

      ELSE
        IF( passed.eq.1) GOTO 260
      ENDIF

200   CONTINUE

260   CONTINUE

      CALL XCHMAG( HOLD )

      RETURN

      ENTRY XAXSOR(X1, Y1)
      KOR=1
      X0=X1
      Y0=Y1
      RETURN

      END

      SUBROUTINE XAXNMG(A)
      COMMON /XAXM35/ NTMAG, ANMAG, ANSIZ
      ANMAG = abs(A)
      NTMAG=1
      RETURN

      ENTRY XAXNSZ(B)
      NTMAG=2
      ANSIZ=abs(B)
      END

      SUBROUTINE XCLEVL(Z,MD, M,N,ZZMAX,ZZMIN,ZZINC,CL,NCNT)
!   TO DETERMINE CONTOUR INCRMENT AND CONTOUR VALUES FOR Z(M,N)
!      REAL Z(MD,1 ),CL(*)                        ! original
      REAL Z(MD,N ),CL(*)
      COMMON /XCLM19/ NMIN, NMAX
      COMMON /XCRF17/CLREF,LCPTN,LABTYP,ICLF,LHILIT,IHLF,KCT0
      COMMON /ZCHOLE/ NHOLE,SPECIA,nvtrbadv
      common /xoutch/ nch
      integer mxset
      NCMIN=NMIN
      NCMAX=NMAX
      ZINC=ZZINC

      mxset = 0
      DO 20 J=1,N
      DO 20 I=1,M

        IF(NHOLE.EQ.1.AND.abs(Z(I,J)-SPECIA).lt.1.0e-6)GOTO 20

        IF( mxset.eq.0) THEN
          ZMAX1= Z(I,J)
          ZMIN1= Z(I,J)
          mxset = 1
        ELSE
          ZMAX1= MAX (ZMAX1,Z(I,J))
          ZMIN1= MIN (ZMIN1,Z(I,J))
        ENDIF
20    CONTINUE

      DIFF=ZMAX1-ZMIN1
      IF(DIFF.le.ABS( ZINC)*1.0E-6) THEN
        WRITE(NCH,'(1x,a,/1x,a)')
     :  'Bad first guess of contour increment or field is constant',
     :  'number of contours is one.'
        NCNT=1
        CL(1)= ZMIN1
        ZZMIN= ZMIN1
        ZZMAX= ZMAX1
        ZZINC= 0.0
        RETURN
      ENDIF

 4    KCOUNT=0
 1    CONTINUE
      EPS=0.001*ZINC
      KCOUNT=KCOUNT+1
      IF( KCOUNT.GT.20) GOTO 998
      KZINC=(ZMIN1-CLREF)/ZINC
      ZMIN=KZINC*ZINC+CLREF
      KZINC=(ZMAX1-CLREF)/ZINC
      ZMAX=KZINC*ZINC+CLREF
      IF(ZMIN1-CLREF.GT.0.0) ZMIN=ZMIN+ZINC
      IF(ZMAX1-CLREF.LT.0.0) ZMAX=ZMAX-ZINC
C
      CLV=ZMIN-ZINC
      NCNT=0
 6    CLV=CLV+ZINC
      IF(CLV-ZMAX-EPS.gt.0.0) GOTO 8
      NCNT=NCNT+1

      IF(NCNT.GT.NCMAX) THEN
        ZINC=ZINC*2
        WRITE(nch,1000) NCMAX, ZINC
 1000   FORMAT(' Number of contours > ',I3,' ,Zinc is doubled. Zinc='
     :   ,E10.3)
        GO TO 1
      ENDIF

      IF( ABS( CLV-CLREF ).LT.EPS ) CLV=CLREF
      CL(NCNT)=CLV
      GOTO 6

 8    CONTINUE

      IF( NCNT.LT.NCMIN) THEN
        ZINC=ZINC/2
        WRITE(nch,2000) NCMIN,ZINC
 2000   FORMAT(' Number of contours < ',I3,' ,Zinc is halved. Zinc='
     :  ,E10.3)
        GO TO 1
      ENDIF

      WRITE(nch,'('' * Number of contours= '',I5,''  MIN='',E12.4,
     : '' MAX='', E12.4,'' INC='',E12.5 )')
     ;    NCNT,ZMIN1,ZMAX1,ZINC
      ZZMAX=ZMAX
      ZZMIN=ZMIN
      ZZINC=ZINC

      RETURN

 998  WRITE(NCH,*)' Contour levels can not be selected by XCNTLV.'
      WRITE(NCH,*)
     :' Plz alter input contour interval or limits of contour number'
      RETURN

      ENTRY XCTREF( CREF)
C Set reference contour level. Default is 0.0 .
        CLREF=CREF
      RETURN

      ENTRY XNCTRS(NMIN1, NMAX1)
C Set upper and lower limit of the number of contours
      NMAX=NMAX1
      NMIN=NMIN1
      RETURN
      END

      SUBROUTINE XCTRHL(Z,X,Y,MD,M,N)
c
c This routine put H,L labels at the maximum and minium
c centers of a contour field.
c Written Oct 13, 1998 by Ming Xue
c
      implicit none

      integer hllabel,hllabel0
      integer llbfmt,laxfmt
      integer NHOLE,nvtrbadv
      real SPECIA
      COMMON /XFMT33/ LBFMT, AXFMT
      COMMON /XFMT34/ LLBFMT,LAXFMT
      COMMON /XHLL36/ hllabel
      COMMON /ZCHOLE/ NHOLE,SPECIA,nvtrbadv
      CHARACTER LBFMT*50, AXFMT*10
      CHARACTER CH*20
      integer lch,i,j,kch,icleng
      real hch,zmin,zmax
c
c Input through argument list
c
      integer md,m,n
      REAL X(MD,*),Y(MD,*), Z(MD,*)

      if(hllabel.eq.0) return

      do j=1,n
        do i=1,m
          IF(.not.(NHOLE.EQ.1.AND.abs(Z(I,J)-SPECIA).lt.1.0e-6))then
          if(hllabel.eq.1) then
            zmax=max(z(max(1,i-1),j),z(i,min(n,j+1)),
     :               z(min(m,i+1),j),z(i,max(1,j-1)))
            zmin=min(z(max(1,i-1),j),z(i,min(n,j+1)),
     :               z(min(m,i+1),j),z(i,max(1,j-1)))
          else
            zmax=max(z(max(1,i-1),j),z(i,min(n,j+1)),
     :               z(max(1,i-2),j),z(i,min(n,j+2)),
     :               z(min(m,i+1),j),z(i,max(1,j-1)),
     :               z(min(m,i+2),j),z(i,max(1,j-2)),
     :      z(max(1,i-1),max(1,j-1)),z(max(1,i-1),min(n,j+1)),
     :      z(min(m,i+1),min(n,j+1)),z(min(m,i+1),max(1,j-1)))

            zmin=min(z(max(1,i-1),j),z(i,min(n,j+1)),
     :               z(max(1,i-2),j),z(i,min(n,j+2)),
     :               z(min(m,i+1),j),z(i,max(1,j-1)),
     :               z(min(m,i+2),j),z(i,max(1,j-2)),
     :      z(max(1,i-1),max(1,j-1)),z(max(1,i-1),min(n,j+1)),
     :      z(min(m,i+1),min(n,j+1)),z(min(m,i+1),max(1,j-1)))
          endif

          if(z(i,j).gt.zmax) then
            call xqchsz(hch)
            call xchsiz(2*hch)
            call xcharc(x(i,j),y(i,j),'H')
            call xchsiz(hch)

            IF( LBFMT(1:LLBFMT).eq.'*') THEN
              CALL XRCH( z(i,j), CH, LCH)
            ELSE
              DO 507 KCH=1,LLBFMT
                IF((LBFMT(KCH:KCH).EQ.'I')
     :              .or.(lbfmt(kch:kch).eq.'i')) THEN
                  WRITE(CH,LBFMT(1:LLBFMT)) NINT(z(i,j))
                  GOTO 508
                ENDIF
 507          CONTINUE

              WRITE( CH, LBFMT(1:LLBFMT) ) z(i,j)
 508          LCH=ICLENG( CH )
              CALL XCHLJ( CH(1:LCH), LCH)
            ENDIF

            call xqchsz(hch)
            CALL Xcharc(x(i,j),y(i,j)-hch,ch(1:lch))
          endif

          if(z(i,j).lt.zmin) then
            call xqchsz(hch)
            call xchsiz(2*hch)
            call xcharc(x(i,j),y(i,j),'L')
            call xchsiz(hch)

            IF( LBFMT(1:LLBFMT).eq.'*') THEN
              CALL XRCH( z(i,j), CH, LCH)
            ELSE
              DO 503 KCH=1,LLBFMT
                IF((LBFMT(KCH:KCH).EQ.'I')
     :              .or.(lbfmt(kch:kch).eq.'i')) THEN
                  WRITE(CH,LBFMT(1:LLBFMT)) NINT(z(i,j))
                  GOTO 504
                ENDIF
 503          CONTINUE

              WRITE( CH, LBFMT(1:LLBFMT) ) z(i,j)
 504          LCH=ICLENG( CH )
              CALL XCHLJ( CH(1:LCH), LCH)
            ENDIF

            call xqchsz(hch)
            CALL Xcharc(x(i,j),y(i,j)-hch,ch(1:lch))
          endif

        endif
        enddo
      enddo

      RETURN

      ENTRY XHLLABL(hllabel0)

      hllabel=hllabel0

      RETURN
      END

      SUBROUTINE XCOLFIL(a,x,y,iwrk,xw,yw,md,m,n, cl0,ncl, mode)
c
c#######################################################################
c
c     Generate color filled contour plots of 2-d field A given its
c     coordinates x and y.
c
c#######################################################################
c
c     INPUT:
c
c       a        2-dimensional slice of data to contour
c       x        x coordinate of grid points in plot space (over on page)
c       y        y coordinate of grid points in plot space (up on page)
c       md       first dimension of a
c       iwrk,xw,yw Work arrays
c       m        number of points in the first dimension of a to be plotted
c       n        second dimension of a
c       cl0      contour levels
c       ncl      Number of contour levels
c       mode     =1,2,3,4. As in XCONTA.
c
c#######################################################################
c
      implicit none
      integer md,m,n
      real a(md,n)
      real x(md,n)
      real y(md,n)
      integer iwrk(*)
      real xw(*),yw(*)  ! dimension for color routine zcontc at least 8*m
      real cl0(*), cl(0:500)

      integer ncl, mode
c
      real zinc             ! contour interval

      real zmax, zmin       ! The real max and min for the field
      real ctrmin, ctrmax

      integer iclrbgn, iclrend
      integer nmin, nmax
      COMMON /XCLM19/ NMIN, NMAX
      integer NHOLE,nvtrbadv
      real SPECIA
      COMMON /ZCHOLE/ NHOLE,SPECIA,nvtrbadv
      integer nch
      common /xoutch/ nch

      integer LCPTN,LABTYP,ICLF,LHILIT,IHLF,KCT0
      real CLREF
      COMMON /XCRF17/ CLREF,LCPTN,LABTYP,ICLF,LHILIT,IHLF,KCT0
c
c     Parameters for color pallete plotting.
c
      integer nctrlvls_max
      parameter(nctrlvls_max=1000)    ! Max. number of contour values
      real ctrlvls(nctrlvls_max)     ! contour values dividing the filled areas
      integer clrindx(nctrlvls_max)  ! plot color index bar color index
      integer nctrlvls      ! Number of contour levels

      common /xcflvls/nctrlvls,ctrlvls,clrindx

      integer icontcopt
      common /xcontc_opt/ icontcopt
c
c     Local variables
c
      integer icol, kolor
      integer mxset,nn
      integer i,j,nclmin,nclmax,kcl, ncl0,ncl1, kcl0
      integer ctrmin_set,ctrmax_set
      real eps,clv,clv1,clv2
      real tem,zmin1,zmax1


c     print*,'inside xcolfil '

      call xqcolor(kolor)

      IF( MODE.LT.1.OR.MODE.GT.4) THEN
        WRITE(NCH,*)
     : ' Input MODE for XCOLFIL not between 1 and 4, job stoped.'
        STOP 999
      ENDIF

      mxset = 0
      DO 2 J=1,N
      DO 2 I=1,M

        IF(NHOLE.EQ.1.AND.abs(a(I,J)-SPECIA).lt.1.0e-6)GOTO 2

        IF( mxset.eq.0) THEN
          ZMAX= a(I,J)
          ZMIN= a(I,J)
          mxset = 1
        ELSE
          ZMAX= MAX (ZMAX,a(I,J))
          ZMIN= MIN (ZMIN,a(I,J))
        ENDIF
 2    CONTINUE
      IF( mxset.eq.0) RETURN

      IF( mode.eq.1) then
        ZINC=CL0(2)-CL0(1)
        CALL XCLEVL(a,MD, M,N,ZMAX1,ZMIN1,ZINC,CL0,NCL)
        IF( NCL.EQ.1 ) RETURN
      ELSEIF( mode.eq.2) then
        NCLMAX=NMAX
        NCLMIN=NMIN
        ZINC=CL0(2)-CL0(1)
        CALL XNCTRS( 0, 500 )
        CALL XCLEVL(a,MD, M,N,ZMAX1,ZMIN1,ZINC,CL0,NCL)
        CALL XNCTRS( NCLMIN, NCLMAX )
        IF( NCL.EQ.1 ) RETURN
      ELSEIF( mode.eq.3) then
        IF( ncl.le.0) return
        ZINC=CL0(2)-CL0(1)
        EPS=0.001*ZINC
        CLV=CL0(1)-ZINC
        kcl = 0
50      continue
        CLV=CLV+ZINC
        IF(CLV-ZMAX.gt.0.0) goto 150
        IF( ABS( CLV-CLREF ).LT. EPS ) CLV=CLREF
        kcl = kcl + 1
        CL0(KCL)=CLV
        GOTO 50
150     CONTINUE
c       ncl = kcl
      ELSEIF( mode.eq.4) then
        ZINC=0.0  ! Undetermined - unequal intervals
      ENDIF

      CALL xqctrlim(ctrmin, ctrmax)

      ctrmax_set = 1
      ctrmin_set = 1
      IF(ctrmax.eq.-9999.0 ) ctrmax_set = 0
      IF(ctrmin.eq.-9999.0 ) ctrmin_set = 0
      IF( ctrmax.eq.0.0.and.ctrmin.eq.0.0 ) THEN
        ctrmin_set = 0
        ctrmax_set = 0
      ENDIF

c     print*,'inside xcolfil 1'
c     print*,'mode=',mode

      IF( MODE.ne.4) THEN

        IF(ctrmin_set.eq.1.and.ctrmax_set.eq.1) then
          ncl0=1
          ncl1=1
          cl(1)=ctrmin
45        clv=cl(ncl1)+zinc
          if(clv.gt.ctrmax+1.0e-5*zinc) goto 450
          ncl1=ncl1+1
          cl(ncl1)=clv
          goto 45
450       continue
        ELSEIF(ctrmin_set.eq.1.and.ctrmax_set.eq.0) then

          ncl0=1
          ncl1=1
          if( ctrmin.gt.zmax) then
            ncl1=1
            goto 550
          endif

          cl(1)=ctrmin
65        clv=cl(ncl1)+zinc
          if(clv.gt.zmax) goto 250
          ncl1=ncl1+1
          cl(ncl1)=clv
          goto 65
250       continue
          if(cl(ncl1).lt.zmax-1.0e-5*zinc)then
            ncl1 =  ncl1+1
            cl(ncl1) = zmax
          endif
        ELSEIF(ctrmin_set.eq.0.and.ctrmax_set.eq.1) then
          ncl0=1
          ncl1=1
          if( ctrmax.lt.zmin) then
            ncl1=1
            goto 550
          endif
          nn = int((ctrmax-zmin)/zinc)
          cl(1)=ctrmax-nn*zinc
75        clv=cl(ncl1)+zinc
          if(clv.gt.ctrmax+1.0e-5*zinc) goto 350
          ncl1=ncl1+1
          cl(ncl1)=clv
          goto 75
350       continue
          if(cl(1).gt.zmin+1.0e-5*zinc)then
            ncl0 = 0
            cl(ncl0) = zmin
          endif
        ELSE

          ncl0 = 1
          ncl1 = ncl
          do i = 1,ncl
            cl(i)=cl0(i)
          enddo
          if(cl(1).gt.zmin)then
            ncl0 = 0
            cl(ncl0) = zmin
          endif
          if(cl(ncl1).lt.zmax)then
            ncl1 =  ncl1+1
            cl(ncl1) = zmax
          endif

        ENDIF

      ELSE ! mode =4

        ncl0 = 1
        ncl1 = ncl
        do i = 1,ncl
          cl(i)=cl0(i)
        enddo

      ENDIF

550   continue

      CALL XQCTRCLR(iclrbgn, iclrend)

      kcl0 = 0
      DO 100 KCL=ncl0, ncl1-1

        CLV1=CL(KCL)
        CLV2=CL(kcl+1)

        if(clv2.lt.clv1) then
          tem = clv2
          clv2 = clv1
          clv1 = tem
        endif

c     print*,'ncl0,ncl1,kcl,clv1,clv2=',ncl0,ncl1,kcl,clv1,clv2

        kcl0 = kcl0+1

        IF(iclrbgn.eq.iclrend) THEN
          icol = iclrbgn
        ELSEIF(iclrbgn.le.iclrend) THEN
          icol= iclrbgn + mod(kcl0-1, iclrend-iclrbgn+1)
        ELSE IF(iclrbgn.gt.iclrend) THEN
          icol= iclrbgn - mod(kcl0-1, iclrbgn-iclrend+1)
        END IF

        call xcolor(icol)
        IF( clv1.lt.zmax+1.0e-10*(clv2-clv1).and.
     :      clv2.gt.zmin-1.0e-10*(clv2-clv1)) then

          if( icontcopt.eq.1) then
            CALL XCONTC(a,x,y,iwrk,xw,yw,md,m,n,clv1,clv2)
          else if( icontcopt.eq.2) then
            CALL XCONTC1(a,x,y,md,m,n,clv1,clv2)
          else if( icontcopt.eq.3) then
            CALL XPIXELFIL(a,x,y,md,m,n,clv1,clv2)
          else
            Print*,'Wrong option for color fill.'
          endif

        endif
c
c       Save values for plotting color palette
c
        clrindx(min(kcl0,nctrlvls_max-1))=icol
        ctrlvls(min(kcl0,nctrlvls_max-1))=CL(KCL)
        ctrlvls(min(kcl0+1,nctrlvls_max))=CL(kcl+1)

100   CONTINUE
      nctrlvls = kcl0+1

      call xcolor(kolor)

      RETURN
      END

      subroutine xpixelfil(z,x,y,md,m,n,c1,c2)
!
! This routine fills pixels with values between c1 and c2
! with one predefined color.
!
! To do: add missing value skipping capability
!
      dimension z(md,*),x(md,*),y(md,*)
      real xcell(4), ycell(4)

      DO j=1,n-1
      DO i=1,m-1

        xcell(1) = x(i,j)
        xcell(2) = x(i+1,j)
        xcell(3) = x(i+1,j+1)
        xcell(4) = x(i,j+1)

        ycell(1) = y(i,j)
        ycell(2) = y(i+1,j)
        ycell(3) = y(i+1,j+1)
        ycell(4) = y(i,j+1)

        zmean = 0.25*(z(i,j)+z(i+1,j)+z(i,j+1)+z(i+1,j+1))

        IF( zmean <= c2 .and. zmean >= c1 ) then
          CALL XFILAREA(Xcell,Ycell,4)
        ENDIF

      ENDDO
      ENDDO

      RETURN

      END


      SUBROUTINE XCTRLIM(ctrmin1, ctrmax1)
c-----------------------------------------------------------------------
c     Set lower and upper limits (the range) of the values beyond which no
c     contour is plotted. Used by XCONTA and XCOLFIL.
c     IF set to -9999.0, then the min or max in the field is used.
c     e.g., CALL XCTRLIM(0.0, -9999.0) will plot all positive contours.
c-----------------------------------------------------------------------
      real ctrmin1, ctrmax1,ctrmin2, ctrmax2
      real ctrmin, ctrmax
      common /xctrmx/ ctrmin, ctrmax
      ctrmin = ctrmin1
      ctrmax = ctrmax1
      RETURN

      ENTRY XQCTRLIM(ctrmin2, ctrmax2)
      ctrmin2 = ctrmin
      ctrmax2 = ctrmax
      RETURN
      END

      SUBROUTINE XCONTA(Z,X,Y,IWRK,MD, M,N, CL0, NCL, MODE )
      IMPLICIT NONE
      INTEGER md,lcptn,labtyp,iclf,lhilit,ihlf,kct0,nmin,nmax
      INTEGER llbfmt,laxfmt,nhole,nvtrbadv,nch
      REAL X(MD,*),Y(MD,*),Z(MD,*), CL0(*), cl(0:900)
      INTEGER IWRK(*)
      COMMON /XCRF17/ CLREF,LCPTN,LABTYP,ICLF,LHILIT,IHLF,KCT0
      COMMON /XCLM19/ NMIN, NMAX
      COMMON /XFMT33/ LBFMT, AXFMT
      COMMON /XFMT34/ LLBFMT,LAXFMT
      COMMON /ZCHOLE/ NHOLE,SPECIA,nvtrbadv
      CHARACTER LBFMT*50, AXFMT*10, LBFM*(*)
      common /xoutch/ nch
      CHARACTER CH*20
      real ctrmin, ctrmax
      common /xctrmx/ ctrmin, ctrmax
      real rmin, rmax

      integer ictr_thick_thin_ratio
      common /ctr_thick_thin_ratio/ ictr_thick_thin_ratio

      integer icol
      integer mxset
      real specia,clref

      integer m,n,ncl,mode,nn,if0,ib0,id1,id0,klbon,kfull,kf
      integer kb1,kf2,kb2,kthick,idref,nocl,nctr_thick_thin_ratio
      integer kf1,lch,kch,icleng,iclrbgn,iclrend,ico,kzero
      real zinc,zmax,zmin,xrag,yrag,xu

      integer i,j,nclmin,nclmax,kcl,ncl0,ncl1,kcl0
      integer ctrmin_set,ctrmax_set
      real eps,clv,clv1,clv2
      real zmin1,zmax1

      IF( MODE.LT.1.OR.MODE.GT.4) THEN
        WRITE(NCH,*)
     : ' Input MODE for XCONTA not between 1 and 3, job stoped.'
        STOP 999
      ENDIF

      mxset = 0
      DO 2 J=1,N
      DO 2 I=1,M

        IF(NHOLE.EQ.1.AND.abs(Z(I,J)-SPECIA).lt.1.0e-6)GOTO 2

        IF( mxset.eq.0) THEN
          ZMAX1= Z(I,J)
          ZMIN1= Z(I,J)
          mxset = 1
        ELSE
          ZMAX1= MAX (ZMAX1,Z(I,J))
          ZMIN1= MIN (ZMIN1,Z(I,J))
        ENDIF
 2    CONTINUE
      IF( mxset.eq.0) RETURN ! All values missing

      IF( mode.eq.1) then
        ZINC=CL0(2)-CL0(1)
        CALL XCLEVL(Z,MD, M,N,ZMAX,ZMIN,ZINC,CL0,NCL)
        IF( NCL.EQ.1 ) RETURN
      ELSEIF( mode.eq.2) then
        NCLMAX=NMAX
        NCLMIN=NMIN
        ZINC=CL0(2)-CL0(1)
        CALL XNCTRS( 0, 500 )
        CALL XCLEVL(Z,MD, M,N,ZMAX,ZMIN,ZINC,CL0,NCL)
        CALL XNCTRS( NCLMIN, NCLMAX )
        IF( NCL.EQ.1 ) RETURN
      ELSEIF( mode.eq.3) then
!       IF( ncl.le.0) return
!       ZINC=CL(2)-CL(1)

        IF( ncl.le.0) return
        ZINC=CL0(2)-CL0(1)
        EPS=0.001*ZINC
        CLV=CL0(1)-ZINC
        kcl = 0
50      continue
        CLV=CLV+ZINC
        IF(CLV-ZMAX1.gt.0.0) goto 150
        IF( ABS( CLV-CLREF ).LT. EPS ) CLV=CLREF
        kcl = kcl + 1
        CL0(KCL)=CLV
        GOTO 50
150     CONTINUE

      ELSEIF( mode.eq.4) then
        ZINC=1.0  ! Undetermined - unequal intervals
      ENDIF

      rmax=zmax1
      rmin=zmin1
      IF( ctrmax.eq.0.0.and.ctrmin.eq.0.0 ) THEN
        rmax=zmax1
        rmin=zmin1
      ELSE
        rmax=ctrmax
        rmin=ctrmin
        IF(ctrmax.eq.-9999.0 ) rmax=zmax1
        IF(ctrmin.eq.-9999.0 ) rmin=zmin1
      END IF

      CALL xqctrlim(ctrmin, ctrmax)

      ctrmax_set = 1
      ctrmin_set = 1
      IF(ctrmax == -9999.0 ) ctrmax_set = 0
      IF(ctrmin == -9999.0 ) ctrmin_set = 0
      IF( ctrmax == 0.0 .AND. ctrmin == 0.0 ) THEN
        ctrmin_set = 0
        ctrmax_set = 0
      ENDIF

c     print*,'inside xcolfil 1'
c     print*,'mode=',mode

      IF( MODE.ne.4) THEN

        IF(ctrmin_set == 1 .AND. ctrmax_set == 1) then
          ncl0=1
          ncl1=1
          cl(1)=ctrmin
45        clv=cl(ncl1)+zinc
          IF ( clv > ctrmax+1.0e-5*zinc .OR. ncl1 >= 900) GOTO 450
          ncl1=ncl1+1
          cl(ncl1)=clv
!          print *, ncl1, clv, zinc, ctrmax
          GOTO 45
450       CONTINUE
        ELSE IF(ctrmin_set.eq.1.and.ctrmax_set.eq.0) THEN

          ncl0=1
          ncl1=1
          IF( ctrmin.gt.zmax1) THEN
            ncl1=1
            GOTO 550
          ENDIF

          cl(1)=ctrmin
65        clv=cl(ncl1)+zinc
          if(clv.gt.zmax1) goto 250
          ncl1=ncl1+1
          cl(ncl1)=clv
          goto 65
250       continue
          if(cl(ncl1).lt.zmax1-1.0e-5*zinc)then
            ncl1 =  ncl1+1
            cl(ncl1) = zmax1
          endif
        ELSEIF(ctrmin_set.eq.0.and.ctrmax_set.eq.1) then
          ncl0=1
          ncl1=1
          if( ctrmax.lt.zmin1) then
            ncl1=1
            goto 550
          endif
          nn = int((ctrmax-zmin1)/zinc)
          cl(1)=ctrmax-nn*zinc
75        clv=cl(ncl1)+zinc
          if(clv.gt.ctrmax+1.0e-5*zinc) goto 350
          ncl1=ncl1+1
          cl(ncl1)=clv
          goto 75
350       continue
          if(cl(1).gt.zmin1+1.0e-5*zinc)then
            ncl0 = 0
            cl(ncl0) = zmin1
          endif
        ELSE

          ncl0 = 1
          ncl1 = ncl
          do i = 1,ncl
            cl(i)=cl0(i)
          enddo
          if(cl(1).gt.zmin1)then
            ncl0 = 0
            cl(ncl0) = zmin1
          endif
          if(cl(ncl1).lt.zmax1)then
            ncl1 =  ncl1+1
            cl(ncl1) = zmax1
          endif

        ENDIF

      ELSE ! mode =4

        ncl0 = 1
        ncl1 = ncl
        do i = 1,ncl
          cl(i)=cl0(i)
        enddo

      ENDIF

550   continue


      CALL XQRANG( XRAG, YRAG )
      XU= SQRT( ABS(XRAG*YRAG) )
      IF0= 12*XU
      IB0=  7*XU
      ID1=  3*XU
      ID0=  3*XU
      CALL XQLBON( KLBON)
      CALL XQFULL( KFULL )
      IF( KFULL.EQ.0 ) CALL XQBRKN( KF1,KB1,KF2,KB2)
      CALL XQTHIK( KTHICK )

      EPS=0.001*ZINC
      CLV=CL(1)-ZINC

!     if( mode.eq.4) then
!       print*,'mode=4', mode
!       print*,'ncl0,ncl1,cl=',ncl0,ncl1,(cl(i),i=ncl0,ncl1)
!       print*,'IHLF=',IHLF
!       print*,'iclf=',iclf
!       iclf = 1
!     endif

      kcl0 = 0
      DO 100 KCL=ncl0,ncl1

        IF( MODE.EQ.3) THEN
          CLV=CLV+ZINC
          IF(CLV-rmax.gt.0.0) goto 100
          IF( ABS( CLV-CLREF ).LT. EPS ) CLV=CLREF
          CL(KCL )=CLV
        ENDIF

        CLV=CL(KCL)

!       print*,'clv, cl(kcl), kcl=', clv, cl(kcl), kcl

        IF(mode.ne.4)then
          IF(clv.lt.rmin .or. clv.gt.rmax) GOTO 100
        ENDIF

        kcl0 = kcl0+1
c
c     Set labeling option for each contour
c
        IDREF=NINT((CL(KCL)-CLREF)/ZINC)
        IF( mode.eq.4) IDREF = KCL-1

        IF((MOD(IDREF,ICLF)==0 .OR. (ncl1==ncl0)) .AND. LABTYP/=0) THEN
          IF( LABTYP.LT.0) THEN
            CALL XLBON
          ELSEIF( LABTYP.eq.1) THEN
            IF( NCL1.EQ.ncl0) THEN
              NOCL=ncl0
            ELSE
              NOCL=IDREF
            ENDIF
            CALL XICH( NOCL  , CH, LCH)
            CALL XLABEL( CH(1:LCH) )
            CALL XLBON
          ELSEIF( LABTYP.eq.2) THEN
            IF( LBFMT(1:LLBFMT).NE.'*') THEN
              DO 503 KCH=1,LLBFMT
                IF((LBFMT(KCH:KCH).EQ.'I')
     :            .or.(lbfmt(kch:kch).eq.'i')) THEN
                  WRITE(CH,LBFMT(1:LLBFMT)) NINT(CL(KCL))
                  GOTO 504
                ENDIF
 503          CONTINUE
              WRITE( CH, LBFMT(1:LLBFMT) ) CL(KCL)
 504          LCH=ICLENG( CH )
              CALL XCHLJ( CH(1:LCH), LCH)
            ELSE
              CALL XRCH( CL(KCL)  , CH, LCH)
            ENDIF
            CALL XLABEL( CH(1:LCH) )
            CALL XLBON
          ENDIF
        ELSE
          CALL XLBOFF
        ENDIF
c
c     Set highlighting option for each contour
c
        IF((MOD(IDREF,IHLF).EQ.0.OR.NCL1.EQ.ncl0).AND.LHILIT.NE.0) THEN
          CALL XTHICK(ictr_thick_thin_ratio)
        ELSE
          CALL XTHICK(1)
        ENDIF

        IF( LCPTN.eq.0 ) THEN
          IF( CL(KCL).lt.0.0 ) THEN
            CALL XBROKN(IF0,IB0,IF0,IB0)
          ELSEIF( abs(CL(KCL)).lt.eps ) THEN
            IF(KCT0.EQ.1) CALL XBROKN(ID1, ID0,ID1,ID0 )
            IF(KCT0.EQ.2) CALL XBROKN(ID1,ID0,IF0,ID0)
            IF(KCT0.EQ.3) THEN
              CALL XFULL
              CALL XTHICK(ictr_thick_thin_ratio)
            ENDIF
          ELSEIF( CL(KCL).gt.0.0 ) THEN
            CALL XFULL
          ENDIF
        ELSEIF( LCPTN.eq.1 ) THEN
          CALL XFULL
        ELSEIF( LCPTN.eq.2 ) THEN
          CALL XBROKN(IF0,IB0,IF0,IB0)
        ELSEIF( LCPTN.eq.4 ) THEN
          CALL XBROKN(ID1, ID0,ID1,ID0 )
        ENDIF

        IF( KCT0.EQ.0.AND.ABS(CLV).LT.1.0e-10*ZINC) GOTO 100

        CALL XQCTRCLR(iclrbgn, iclrend)
        IF(iclrbgn.eq.iclrend) THEN
          icol = iclrbgn
        ELSEIF(iclrbgn.le.iclrend) THEN
          icol= iclrbgn + mod(KCL0-1, iclrend-iclrbgn+1)
        ELSE IF(iclrbgn.gt.iclrend) THEN
          icol= iclrbgn - mod(KCL0-1, iclrbgn-iclrend+1)
        END IF

        call xcolor(icol)

        IF(clv.lt.zmin1.or. clv.gt.zmax1) GOTO 100

        IF( NHOLE.EQ.1 ) THEN
          CALL XCONTJ(Z,X,Y,IWRK,MD,M,N,CLV,SPECIA)
        ELSE
          CALL XCONTR(Z,X,Y,IWRK,MD, M,N,CLV )
        ENDIF

100   CONTINUE

      call XCTRHL(Z,X,Y,MD,M,N)

      IF( KFULL.EQ.1) CALL XFULL
      IF( KFULL.NE.1) CALL XBROKN( KF1,KB1,KF2,KB2)
      IF( KLBON.EQ.1) CALL XLBON
      IF( KLBON.NE.1) CALL XLBOFF
      CALL XTHICK( KTHICK )
      RETURN

      ENTRY XCLFMT( LBFM )
        LLBFMT=LEN(LBFM)
        LBFMT=LBFM
      RETURN

      ENTRY XQCZRO(KZERO)
        KZERO=KCT0
      RETURN

      ENTRY XCTR_THICK_THIN_RATIO(nctr_thick_thin_ratio)
        ictr_thick_thin_ratio = nctr_thick_thin_ratio
      RETURN

      END

      SUBROUTINE XCTRBADV(MHOLE)
        COMMON /ZCHOLE/ NHOLE,SPECIA,nvtrbadv
        NHOLE=MHOLE
      RETURN

      ENTRY XBADVAL(SPECM)
        SPECIA=SPECM
      RETURN

      ENTRY XVTRBADV(MHOLE)
        nvtrbadv = MHOLE
      RETURN

      END

      SUBROUTINE XCTRCLR(klrbgn, klrend)
        implicit none
        integer klrbgn,klrend    ! Beginning and ending colors of contours
        integer iclrbgn,iclrend  ! Beginning and ending colors of contours
        common /xctrclor/iclrbgn,iclrend
        iclrbgn = klrbgn
        iclrend = klrend
      RETURN
      END

      SUBROUTINE XQCTRCLR(klrbgn, klrend)
        implicit none
        integer klrbgn,klrend    ! Beginning and ending colors of contours
        integer iclrbgn,iclrend  ! Beginning and ending colors of contours
        common /xctrclor/iclrbgn,iclrend
        klrbgn = iclrbgn
        klrend = iclrend
      RETURN
      END

      SUBROUTINE ZCONTA(Z,ZG,IWRK,MD,M ,N ,CL,NCL, MODE)
      DIMENSION ZG(MD ,*),Z(MD ,*),IWRK(M ,*),CL(*)
      COMMON /XCRF17/ CLREF,LCPTN,LABTYP,ICLF,LHILIT,IHLF,KCT0
      COMMON /XCLM19/ NMIN, NMAX
      COMMON /XFMT33/ LBFMT, AXFMT
      COMMON /XFMT34/ LLBFMT,LAXFMT
      COMMON /ZCHOLE/ NHOLE,SPECIA,nvtrbadv
      CHARACTER LBFMT*50, AXFMT*10
      CHARACTER CH*20
      COMPLEX   ZG
      IF( MODE.LT.1.OR.MODE.GT.3) THEN
        PRINT*,' Input MODE for XCONTB not between 1 and 3, job stoped.'
        STOP 999
      ENDIF
      GOTO ( 50,51,52 ) MODE
 50   ZINC=CL(2)-CL(1)
      CALL XCLEVL(Z,MD, M,N,ZMAX,ZMIN,ZINC,CL,NCL)
      IF( NCL.EQ.1 ) RETURN
      GOTO 55
 51   NCLMAX=NMAX
      NCLMIN=NMIN
      ZINC=CL(2)-CL(1)
      CALL XNCTRS( 0, 500 )
      CALL XCLEVL(Z,MD, M,N,ZMAX,ZMIN,ZINC,CL,NCL)
      IF( NCL.EQ.1 ) RETURN
      CALL XNCTRS( NCLMIN, NCLMAX )
      GOTO 55
 52   IF( NCL-1 ) 101, 102, 103
 101  RETURN
 102  ZINC=1.0
      GOTO 104
 103  ZINC=CL(2)-CL(1)
 104  CONTINUE

 55   CONTINUE

      mxset = 0
      DO 2 J=1,N
        DO 2 I=1,M

        IF(NHOLE.EQ.1.AND.abs(Z(I,J)-SPECIA).lt.1.0e-6)GOTO 2

        IF( mxset.eq.0) THEN
          ZMAX1= Z(I,J)
          ZMIN1= Z(I,J)
          mxset = 1
        ELSE
          ZMAX1= MAX (ZMAX1,Z(I,J))
          ZMIN1= MIN (ZMIN1,Z(I,J))
        ENDIF
 2    CONTINUE

      CALL XQRANG( XRAG, YRAG )
      XU= MIN( XRAG, YRAG)
      IF0= 20*XU
      IB0=  7*XU
      ID1=  7*XU
      ID0= 15*XU
      CALL XQLBON( KLBON)
      CALL XQFULL( KFULL )
      IF( KFULL.EQ.0 ) CALL XQBRKN( KF1,KB1,KF2,KB2)
      CALL XQTHIK( KTHICK )

      EPS=0.001*ZINC
      CLV=CL(1)-ZINC
      DO 10 KCL=1,NCL
      IF( MODE.EQ.3) THEN
        CLV=CLV+ZINC
        IF(CLV-ZMAX1   ) 4,4,10
 4      IF( ABS( CLV-CLREF ).LT. EPS        ) CLV=CLREF
        CL(KCL )=CLV
      ENDIF
        IDREF=NINT((CL(KCL)-CLREF)/ZINC)
        IF((MOD(IDREF,ICLF).EQ.0.OR.NCL.EQ.1).AND. LABTYP.NE.0) THEN
          IF( LABTYP.LT.0) GOTO 46
          GOTO  (41,42) LABTYP
 41       IF( NCL.EQ. 1) THEN
            NOCL=1
          ELSE
            NOCL=IDREF
          ENDIF
          CALL XICH( NOCL,CH,LCH)
          GOTO 43
 42       CONTINUE
C         IF( FLOAT( INT( ZINC ) ).EQ. ZINC .AND. FLOAT( INT(CLREF))
C    :        .EQ. CLREF) THEN
C           CALL XICH( INT(CL(KCL)), CH, LCH )
          IF( LBFMT(1:LLBFMT).NE.'*') THEN
            WRITE( CH, LBFMT(1:LLBFMT) ) CL(KCL)
            LCH=ICLENG( CH )
          ELSE
            CALL XRCH( CL(KCL),CH,LCH)
          ENDIF
 43       CONTINUE
          CALL XLABEL( CH(1:LCH) )
 46       CALL XLBON
        ELSE
          CALL XLBOFF
        ENDIF

        IF((MOD(IDREF,IHLF).EQ.0.OR.NCL.EQ.1).AND. LHILIT.NE.0 )THEN
          CALL XTHICK(2)
        ELSE
          CALL XTHICK(1)
        ENDIF

        GOTO ( 30,31,32,33 ) LCPTN+1
 30     IF( ABS( CL(KCL) ).LT. EPS ) GOTO 22
          IF( CL(KCL)) 21,22,23
 21         CALL XBROKN(IF0,IB0,IF0,IB0)
            GOTO 24
 22         CALL XBROKN(ID1,ID0,ID1,ID0 )
            GOTO 24
 23         CALL XFULL
 24       CONTINUE
          GOTO 35
 31         CALL XFULL
          GOTO 35
 32         CALL XBROKN(IF0,IB0,IF0,IB0)
          GOTO 35
 33         CALL XBROKN(ID1,ID0,ID1,ID0 )
 35       CONTINUE
        CLV=CL(KCL)
        IF( KCT0.EQ.0.AND.ABS(CLV).LT.0.1*ZINC) GOTO 10
        IF( NHOLE.EQ.1 ) THEN
          CALL ZCONTJ(Z,ZG,IWRK,MD,M,N,CLV,SPECIA)
        ELSE
            CALL ZCONTR(Z,ZG ,IWRK,MD, M,N,CLV )
        ENDIF
 10   CONTINUE
        IF( KFULL.EQ.1) CALL XFULL
        IF( KFULL.NE.1) CALL XBROKN( KF1,KB1,KF2,KB2)
        IF( KLBON.EQ.1) CALL XLBON
        IF( KLBON.NE.1) CALL XLBOFF
        CALL XTHICK( KTHICK )
      RETURN
      END

      SUBROUTINE XCMIXL
      COMMON/XCRF17/CLREF,LCPTN,LABTYP,ICLF,LHILIT,IHLF,KCT0
C Contour plotting pattern is set so that lines are dash,dotted,solid
C for negative ,zero, positve values respectively. This is default.
        LCPTN=0
      RETURN

      ENTRY XCFULL
        LCPTN=1
      RETURN

      ENTRY XCDASH
C Set contour plotting pattern as dash   lines.
        LCPTN=2
      RETURN

      ENTRY XCDOT
        LCPTN=3
      RETURN

      ENTRY XCLTYP( LTYPE)
C Define type of labels on contours.
C LTYPE-- parameter controling contour labeling.
C LTYPE <0,  label is specified by user through XLABEL,
C       =0,  no labeling is done.
C       =1,  label the contour number,  number=0 for zero contour.
C       =2,  label the contour values.
C By default LTYPE=2.
C Note setting LTYPE=0 is the only way to suppress labels outside
C routine XCONTA as XLBON and XLBOFF are called inside XCONTA.
      LABTYP=LTYPE
      RETURN

      ENTRY XCLFRQ( NCLF)
C Set contour labeling frequency so that every NCLFth contour relative
C to reference contour is labeled. Default NCLF=2.
      ICLF=NCLF
      RETURN

      ENTRY XHILIT( KHILIT )
      LHILIT=KHILIT
      RETURN

      ENTRY XHLFRQ( NHLF )
      IHLF=NHLF
      RETURN

      ENTRY XCZERO( KCZERO )
C Option of zero contour plotting.
C KCZERO=0, zero line is suppressed, by default KCZERO=1.
      KCT0=KCZERO
      RETURN
      END

      SUBROUTINE XCONTR(ZG,X,Y,IWRK,MD,MG,JG,CV)
      DIMENSION ZG(MD ,*),X(MD ,*),Y(MD,*),IWRK(MG ,*)
C* The final edition of the contouring package  2nd ed
C*  Zhang Zuojun, Jan. 1988

      DOUBLE PRECISION CVn, normscl
      DOUBLE PRECISION H5n
      INTEGER          normexp

C     IFUN1(K)=K+MG*((MGP-K)/MGP-K/MGP)
      D(P1,P2,B1,B2   )=B1+(CV-P1)*(B2-B1)/(P2-P1)

C     Normalize CV and H5
      IF (CV == 0) THEN
        normexp = 0
      ELSE
        normexp = ANINT(LOG10(ABS(CV)))
      END IF
      normscl = 10**(-1.*normexp)
      cvn = dble(cv)*normscl
c     write(0,*) cv,normexp, normscl, cvn

      MGP=MG+1
      JGP=JG+1
      DO 4 J=1,JG
      DO 4 I=1,MG
    4 IWRK(I,J)=0
      DO 1 JJ=1,2*(MG+JG-2)
      IF(JJ.LT.MG) THEN
          I4=JJ
          J4=1
          ISW=1
      ELSEIF(JJ.LT.MG+JG-1) THEN
          I4=MG
          J4=JJ-MG+1
          ISW=4
      ELSEIF(JJ.LT.MG+MG+JG-2) THEN
          I4=MG+MG+JG-JJ-1
          J4=JG
          ISW=3
      ELSEIF(JJ.LT.MG+MG+JG+JG-3) THEN
          I4=1
          J4=MG+MG+JG+JG-2-JJ
          ISW=2
      ENDIF
      INI=MOD(ISW  ,2)*(1-2*(MOD(ISW,4)/2))
      INJ=MOD(ISW+1,2)*(1-2*(MOD(ISW,4)/2))
      I1=I4+INI
      J1=J4+INJ
      IF(I1.EQ.0.OR.I1.EQ.MGP.OR.J1.EQ.0.OR.J1.EQ.JGP)GOTO 1
      H1=ZG(I1,J1)
      H4=ZG(I4,J4)
      IF(H1.GE.CV.OR.H4.LT.CV ) GOTO 1
      X1= X(I1,J1)
      X4= X(I4,J4)
      Y1= Y(I1,J1)
      Y4= Y(I4,J4)
      XA=D(H4,H1,X4,X1)
      YA=D(H4,H1,Y4,Y1)
      CALL XCURUP( XA, YA )
      I2=I1+MOD(ISW-1,2)*(1-2*((ISW-1)/2))
      J2=J1+MOD(ISW  ,2)*(1-2*((ISW-1)/2))
      I3=I4+MOD(ISW-1,2)*(1-2*((ISW-1)/2))
      J3=J4+MOD(ISW  ,2)*(1-2*((ISW-1)/2))
  201 H1=ZG(I1,J1)
      H2=ZG(I2,J2)
      H3=ZG(I3,J3)
      H4=ZG(I4,J4)
      H5n=0.25*(dble(H1)+dble(H2)+dble(H3)+dble(H4))
      H5 = H5n
      H5n= normscl*H5n
      X1= X(I1,J1)
      X2= X(I2,J2)
      X3= X(I3,J3)
      X4= X(I4,J4)
      Y1= Y(I1,J1)
      Y2= Y(I2,J2)
      Y3= Y(I3,J3)
      Y4= Y(I4,J4)
      IF(H2-CV) 52,53,53
   52 IF(H3-CV) 63,62,62
   53 IF(H3-CV) 54,61,61
c  54 IF(H5-CV) 63,61,61
   54 if ( (H5n-CVn) > -1.0E-5) then
        go to 61
      else
        go to 63
      end if
   61 ISA=1
      XB=D(H1,H2,X1,X2)
      YB=D(H1,H2,Y1,Y2)
      I4=I2
      J4=J2
      GOTO 60
   62 ISA=2
      XB=D(H2,H3,X2,X3)
      YB=D(H2,H3,Y2,Y3)
      I1=I2
      J1=J2
      I4=I3
      J4=J3
      GOTO 60
   63 ISA=3
      XB=D(H3,H4,X3,X4)
      YB=D(H3,H4,Y3,Y4)
      I1=I3
      J1=J3
   60 ISW=MOD(ISW-ISA+5,4)+1
      I2=I1+MOD(ISW-1,2)*(1-2*((ISW-1)/2))
      J2=J1+MOD(ISW  ,2)*(1-2*((ISW-1)/2))
      I3=I4+MOD(ISW-1,2)*(1-2*((ISW-1)/2))
      J3=J4+MOD(ISW  ,2)*(1-2*((ISW-1)/2))
      IF( I2.EQ.0.OR.I3.EQ.0.OR.I2.EQ.MGP.OR.I3.EQ.MGP  .OR.
     :    J2.EQ.0.OR.J3.EQ.0.OR.J2.EQ.JGP.OR.J3.EQ.JGP) THEN
          CALL XCURDN( XB, YB,0 , 1 )
      ELSE
          IF(XB.NE.XA. OR.YB.NE.YA) CALL XCURDN( XB , YB, 0 ,0)
          XA=XB
          YA=YB
          IWRK(I1,J1)=1
          IWRK(I4,J4)=1
          GOTO 201
      ENDIF
   1  CONTINUE
      DO 2 J=2,JG-1
      DO 2 I=1,MG-1
      ISW=1
      I10=I+1
      J10=J
      I40=I
      J40=J
      IF(IWRK(I10,J10).EQ.1.AND.IWRK(I40,J40).EQ.1) GOTO 2
      H1=ZG(I10,J10)
      H4=ZG(I40,J40)
      IF(H1.GE.CV.OR.H4.LT.CV ) GOTO 2
      I1=I10
      J1=J10
      I4=I40
      J4=J40
      X1= X(I1,J1)
      X4= X(I4,J4)
      Y1= Y(I1,J1)
      Y4= Y(I4,J4)
      XA=D(H4,H1,X4,X1)
      YA=D(H4,H1,Y4,Y1)
      CALL XCURUP( XA, YA )
  101 I2=I1+MOD(ISW-1,2)*(1-2*((ISW-1)/2))
      J2=J1+MOD(ISW  ,2)*(1-2*((ISW-1)/2))
      I3=I4+MOD(ISW-1,2)*(1-2*((ISW-1)/2))
      J3=J4+MOD(ISW  ,2)*(1-2*((ISW-1)/2))
      H1=ZG(I1,J1)
      H2=ZG(I2,J2)
      H3=ZG(I3,J3)
      H4=ZG(I4,J4)
      H5n=0.25*(dble(H1)+dble(H2)+dble(H3)+dble(H4))
      H5 = H5n
      H5n= normscl*H5n
      X1= X(I1,J1)
      X2= X(I2,J2)
      X3= X(I3,J3)
      X4= X(I4,J4)
      Y1= Y(I1,J1)
      Y2= Y(I2,J2)
      Y3= Y(I3,J3)
      Y4= Y(I4,J4)
      IF(H2-CV) 12,13,13
   12 IF(H3-CV) 23,22,22
   13 IF(H3-CV) 14,21,21
c  14 IF(H5-CV) 23,21,21
   14 if ( (H5n-CVn) > -1.0E-5) then
        go to 21
      else
        go to 23
      end if
   21 ISA=1
      XB=D(H1,H2,X1,X2)
      YB=D(H1,H2,Y1,Y2)
      I4=I2
      J4=J2
      GOTO 30
   22 ISA=2
      XB=D(H2,H3,X2,X3)
      YB=D(H2,H3,Y2,Y3)
      I1=I2
      J1=J2
      I4=I3
      J4=J3
      GOTO 30
   23 ISA=3
      XB=D(H3,H4,X3,X4)
      YB=D(H3,H4,Y3,Y4)
      I1=I3
      J1=J3
   30 IF( I1.EQ.I10.AND.J1.EQ.J10.AND.I4.EQ.I40.AND.J4.EQ.J40) THEN
          CALL XCURDN(XB,YB,1,1)
      ELSE
          IF(XB.NE.XA. OR.YB.NE.YA) CALL XCURDN( XB , YB, 1 ,0)
          XA=XB
          YA=YB
          IWRK(I1,J1)=1
          IWRK(I4,J4)=1
          ISW=MOD(ISW-ISA+5,4)+1
          GOTO 101
      ENDIF
    2 CONTINUE
      CALL XLPNUP( X(1,1) ,Y(1,1) )
      RETURN
      END

      SUBROUTINE XCONTJ(ZG,X,Y,IWRK,MD,MG,JG,CV,SPEC)
      DIMENSION ZG(MD ,*),X(MD ,*),Y(MD,*),IWRK(MG ,*)
C*   New update for contouring allowing special value holes (SPEC)
C*   The second edition of the contour tracing
C*   Zhang Zuojun, Jan. 1988
C*   New update including contouring on triagle grids
C*   When MODE=0 contouring perform on retangular grids (default)
C*   When MODE=1 contouring perform on triangular grids .
C     IFUN1(K)=K+MG*((MGP-K)/MGP-K/MGP)
c
c    Converted by Ming Xue, Oct. 1993 to use real arrays for
c    grid coordinates.
c
      DOUBLE PRECISION CVn, normscl
      DOUBLE PRECISION H5n
      INTEGER          normexp

      D(P1,P2,B1,B2   )=B1+(CV-P1)*(B2-B1)/(P2-P1)

C     Normalize CV and H5
      IF (CV == 0) THEN
        normexp = 0
      ELSE
        normexp = ANINT(LOG10(ABS(CV)))
      END IF
      normscl = 10**(-1.*normexp)
      cvn = dble(cv)*normscl
c     write(0,*) cv,normexp, normscl, cvn

      CALL ZQCONM(MODE)
      DUM=SPEC
      MGP=MG+1
      JGP=JG+1
      DO 4 J=1,JG
      DO 4 I=1,MG
    4 IWRK(I,J)=0
      DO 1 JJ=1,2*(MG+JG-2)
      IF(JJ.LT.MG) THEN
          I4=JJ
          J4=1
          ISW=1
      ELSEIF(JJ.LT.MG+JG-1) THEN
          I4=MG
          J4=JJ-MG+1
          ISW=4
      ELSEIF(JJ.LT.MG+MG+JG-2) THEN
          I4=MG+MG+JG-JJ-1
          J4=JG
          ISW=3
      ELSEIF(JJ.LT.MG+MG+JG+JG-3) THEN
          I4=1
          J4=MG+MG+JG+JG-2-JJ
          ISW=2
      ENDIF
      INI=MOD(ISW  ,2)*(1-2*(MOD(ISW,4)/2))
      INJ=MOD(ISW+1,2)*(1-2*(MOD(ISW,4)/2))
      I1=I4+INI
      J1=J4+INJ
      IF(I1.EQ.0.OR.I1.EQ.MGP.OR.J1.EQ.0.OR.J1.EQ.JGP)GOTO 1
      H1=ZG(I1,J1)
      H4=ZG(I4,J4)
      IF(H1.GE.CV.OR.H4.LT.CV ) GOTO 1
      X1= X(I1,J1)
      X4= X(I4,J4)
      Y1= Y(I1,J1)
      Y4= Y(I4,J4)
      XA=D(H4,H1,X4,X1)
      YA=D(H4,H1,Y4,Y1)
      CALL XCURUP( XA, YA )

      I2=I1+MOD(ISW-1,2)*(1-2*((ISW-1)/2))
      J2=J1+MOD(ISW  ,2)*(1-2*((ISW-1)/2))
      I3=I4+MOD(ISW-1,2)*(1-2*((ISW-1)/2))
      J3=J4+MOD(ISW  ,2)*(1-2*((ISW-1)/2))
  201 H1=ZG(I1,J1)
      H2=ZG(I2,J2)
      H3=ZG(I3,J3)
      H4=ZG(I4,J4)
      H5n=0.25*(dble(H1)+dble(H2)+dble(H3)+dble(H4))
      H5 = H5n
      H5n= normscl*H5n
      IF(H1.EQ.DUM.OR.H2.EQ.DUM.OR.H3.EQ.DUM.OR.H4.EQ.DUM)THEN
        ISWTCH=0
      ELSE
        ISWTCH=1
      ENDIF
      X1= X(I1,J1)
      X2= X(I2,J2)
      X3= X(I3,J3)
      X4= X(I4,J4)
      Y1= Y(I1,J1)
      Y2= Y(I2,J2)
      Y3= Y(I3,J3)
      Y4= Y(I4,J4)
      IF(MODE.EQ.1) THEN
        X5=0.25*(X1+X2+X3+X4)
        Y5=0.25*(Y1+Y2+Y3+Y4)
      ENDIF
      IF(H2-CV) 52,53,53
   52 IF(H3-CV) 63,62,62
   53 IF(H3-CV) 54,61,61
c* 54 IF(H5-CV) 63,61,61
   54 if ( (H5n-CVn) > -1.0E-5) then
        go to 61
      else
        go to 63
      end if
   61 ISA=1
      IF(MODE.EQ.1.AND.ISWTCH.EQ.1) THEN
          IF(H5.LT.CV.AND.H3.GE.CV) THEN
              XC=D(H4,H5,X4,X5)
              YC=D(H4,H5,Y4,Y5)
              IF(XC.NE.XA.or.YC.NE.YA)
     :           CALL XCURDN(XC,YC,0,0)
              XA=XC
              YA=YC
              XC=D(H3,H5,X3,X5)
              YC=D(H3,H5,Y3,Y5)
              IF(XC.NE.XA.or.YC.NE.YA)
     :           CALL XCURDN(XC,YC,0,0)
              XA=XC
              YA=YC
              XC=D(H2,H5,X2,X5)
              YC=D(H2,H5,Y2,Y5)
              IF(XC.NE.XA.or.YC.NE.YA)
     :           CALL XCURDN(XC,YC,0,0)
              XA=XC
              YA=YC
          ELSE
              XC=D(H1,H5,X1,X5)
              YC=D(H1,H5,Y1,Y5)
              IF(XC.NE.XA.or.YC.NE.YA) CALL XCURDN(XC,YC,0,0)
              XA=XC
              YA=YC
          ENDIF
      ENDIF
      XB=D(H1,H2,X1,X2)
      YB=D(H1,H2,Y1,Y2)
      I4=I2
      J4=J2
      GOTO 60
   62 ISA=2
      IF(MODE.EQ.1.AND.ISWTCH.EQ.1) THEN
          IF(H5.LT.CV) THEN
              XC=D(H4,H5,X4,X5)
              YC=D(H4,H5,Y4,Y5)
              IF(XC.NE.XA.or.YC.NE.YA) CALL XCURDN(XC,YC,0,0)
              XA=XC
              YA=YC
              XC=D(H3,H5,X3,X5)
              YC=D(H3,H5,Y3,Y5)
              IF(XC.NE.XA.or.YC.NE.YA) CALL XCURDN(XC,YC,0,0)
              XA=XC
              YA=YC
          ELSE
              XC=D(H1,H5,X1,X5)
              YC=D(H1,H5,Y1,Y5)
              IF(XC.NE.XA.or.YC.NE.YA) CALL XCURDN(XC,YC,0,0)
              XA=XC
              YA=YC
              XC=D(H2,H5,X2,X5)
              YC=D(H2,H5,Y2,Y5)
              IF(XC.NE.XA.or.YC.NE.YA) CALL XCURDN(XC,YC,0,0)
              XA=XC
              YA=YC
          ENDIF
      ENDIF
      XB=D(H2,H3,X2,X3)
      YB=D(H2,H3,Y2,Y3)
      I1=I2
      J1=J2
      I4=I3
      J4=J3
      GOTO 60
   63 ISA=3
      IF(MODE.EQ.1.AND.ISWTCH.EQ.1) THEN
          IF(H5.GE.CV.AND.H2.LT.CV) THEN
              XC=D(H1,H5,X1,X5)
              YC=D(H1,H5,Y1,Y5)
              IF(XC.NE.XA.or.YC.NE.YA) CALL XCURDN(XC,YC,0,0)
              XA=XC
              YA=YC
              XC=D(H2,H5,X2,X5)
              YC=D(H2,H5,Y2,Y5)
              IF(XC.NE.XA.or.YC.NE.YA) CALL XCURDN(XC,YC,0,0)
              XA=XC
              YA=YC
              XC=D(H3,H5,X3,X5)
              YC=D(H3,H5,Y3,Y5)
              IF(XC.NE.XA.or.YC.NE.YA) CALL XCURDN(XC,YC,0,0)
              XA=XC
              YA=YC
          ELSE
              XC=D(H4,H5,X4,X5)
              YC=D(H4,H5,Y4,Y5)
              IF(XC.NE.XA.or.YC.NE.YA) CALL XCURDN(XC,YC,0,0)
              XA=XC
              YA=YC
          ENDIF
      ENDIF
      XB=D(H3,H4,X3,X4)
      YB=D(H3,H4,Y3,Y4)
      I1=I3
      J1=J3
   60 ISW=MOD(ISW-ISA+5,4)+1
      I2=I1+MOD(ISW-1,2)*(1-2*((ISW-1)/2))
      J2=J1+MOD(ISW  ,2)*(1-2*((ISW-1)/2))
      I3=I4+MOD(ISW-1,2)*(1-2*((ISW-1)/2))
      J3=J4+MOD(ISW  ,2)*(1-2*((ISW-1)/2))
      IF( I2.EQ.0.OR.I3.EQ.0.OR.I2.EQ.MGP.OR.I3.EQ.MGP  .OR.
     :    J2.EQ.0.OR.J3.EQ.0.OR.J2.EQ.JGP.OR.J3.EQ.JGP )THEN
          IF(ISWTCH.EQ.1) THEN
            CALL XCURDN(XB,YB,0,1)
          ELSE
            CALL XCURUP(XB,YB)
          ENDIF
      ELSE
        IF(ISWTCH.EQ.1.AND.(XB.NE.XA.or.YB.NE.YA))THEN
          IF(ZG(I2,J2).EQ.DUM.OR.ZG(I3,J3).EQ.DUM)THEN
              CALL XCURDN(XB,YB,0,0)
          ELSE
              CALL XCURDN(XB,YB,0,1)
          ENDIF
        ELSE
            CALL XCURUP(XB,YB)
        ENDIF
        XA=XB
        YA=YB
        IWRK(I1,J1)=1
        IWRK(I4,J4)=1
        GOTO 201
      ENDIF
   1  CONTINUE
      DO 2 J=2,JG-1
      DO 2 I=1,MG-1
      ISW=1
      I10=I+1
      J10=J
      I40=I
      J40=J
      IF(IWRK(I10,J10).EQ.1.AND.IWRK(I40,J40).EQ.1) GOTO 2
      H1=ZG(I10,J10)
      H4=ZG(I40,J40)
      IF(H1.GE.CV.OR.H4.LT.CV ) GOTO 2
      I1=I10
      J1=J10
      I4=I40
      J4=J40
      X1= X(I1,J1)
      Y1= Y(I1,J1)
      X4= X(I4,J4)
      Y4= Y(I4,J4)
      XA=D(H4,H1,X4,X1)
      YA=D(H4,H1,Y4,Y1)
      CALL XCURUP(XA,YA)
      I2=I1+MOD(ISW-1,2)*(1-2*((ISW-1)/2))
      J2=J1+MOD(ISW  ,2)*(1-2*((ISW-1)/2))
      I3=I4+MOD(ISW-1,2)*(1-2*((ISW-1)/2))
      J3=J4+MOD(ISW  ,2)*(1-2*((ISW-1)/2))
  101 H1=ZG(I1,J1)
      H2=ZG(I2,J2)
      H3=ZG(I3,J3)
      H4=ZG(I4,J4)
      H5n=0.25*(dble(H1)+dble(H2)+dble(H3)+dble(H4))
      H5=H5n
      H5n = normscl*H5n
      IF(H1.EQ.DUM.OR.H2.EQ.DUM.OR.H3.EQ.DUM.OR.H4.EQ.DUM)THEN
        ISWTCH=0
      ELSE
        ISWTCH=1
      ENDIF
      X1= X(I1,J1)
      X2= X(I2,J2)
      X3= X(I3,J3)
      X4= X(I4,J4)
      Y1= Y(I1,J1)
      Y2= Y(I2,J2)
      Y3= Y(I3,J3)
      Y4= Y(I4,J4)

      IF(MODE.EQ.1) THEN
        X5=0.25*(X1+X2+X3+X4)
        Y5=0.25*(Y1+Y2+Y3+Y4)
      ENDIF

      IF(H2-CV) 12,13,13
   12 IF(H3-CV) 23,22,22
   13 IF(H3-CV) 14,21,21
c  14 IF(H5-CV) 23,21,21
   14 IF ( (H5n-CVn) > -1.0E-5) THEN
        GO TO 21
      ELSE
        GO TO 23
      END IF
   21 ISA=1
      IF(MODE.EQ.1.AND.ISWTCH.EQ.1) THEN
          IF(H5.LT.CV.AND.H3.GE.CV) THEN
              XC=D(H4,H5,X4,X5)
              YC=D(H4,H5,Y4,Y5)
              IF(XC.NE.XA.or.YC.NE.YA) CALL XCURDN(XC,YC,1,0)
              XA=XC
              YA=YC
              XC=D(H3,H5,X3,X5)
              YC=D(H3,H5,Y3,Y5)
              IF(XC.NE.XA.or.YC.NE.YA) CALL XCURDN(XC,YC,1,0)
              XA=XC
              YA=YC
              XC=D(H2,H5,X2,X5)
              YC=D(H2,H5,Y2,Y5)
              IF(XC.NE.XA.or.YC.NE.YA) CALL XCURDN(XC,YC,1,0)
              XA=XC
              YA=YC
          ELSE
              XC=D(H1,H5,X1,X5)
              YC=D(H1,H5,Y1,Y5)
              IF(XC.NE.XA.or.YC.NE.YA) CALL XCURDN(XC,YC,1,0)
              XA=XC
              YA=YC
          ENDIF
      ENDIF
      XB=D(H1,H2,X1,X2)
      YB=D(H1,H2,Y1,Y2)
      I4=I2
      J4=J2
      GOTO 30
   22 ISA=2
      IF(MODE.EQ.1.AND.ISWTCH.EQ.1) THEN
          IF(H5.LT.CV) THEN
              XC=D(H4,H5,X4,X5)
              YC=D(H4,H5,Y4,Y5)
              IF(XC.NE.XA.or.YC.NE.YA) CALL XCURDN(XC,YC,1,0)
              XA=XC
              YA=YC
              XC=D(H3,H5,X3,X5)
              YC=D(H3,H5,Y3,Y5)
              IF(XC.NE.XA.or.YC.NE.YA) CALL XCURDN(XC,YC,1,0)
              XA=XC
              YA=YC
          ELSE
              XC=D(H1,H5,X1,X5)
              YC=D(H1,H5,Y1,Y5)
              IF(XC.NE.XA.or.YC.NE.YA) CALL XCURDN(XC,YC,1,0)
              XA=XC
              YA=YC
              XC=D(H2,H5,X2,X5)
              YC=D(H2,H5,Y2,Y5)
              IF(XC.NE.XA.or.YC.NE.YA) CALL XCURDN(XC,YC,1,0)
              XA=XC
              YA=YC
          ENDIF
      ENDIF
      XB=D(H2,H3,X2,X3)
      YB=D(H2,H3,Y2,Y3)
      I1=I2
      J1=J2
      I4=I3
      J4=J3
      GOTO 30
   23 ISA=3
      IF(MODE.EQ.1.AND.ISWTCH.EQ.1) THEN
          IF(H5.GE.CV.AND.H2.LT.CV) THEN
              XC=D(H1,H5,X1,X5)
              YC=D(H1,H5,Y1,Y5)
              IF(XC.NE.XA.or.YC.NE.YA) CALL XCURDN(XC,YC,1,0)
              XA=XC
              YA=YC
              XC=D(H2,H5,X2,X5)
              YC=D(H2,H5,Y2,Y5)
              IF(XC.NE.XA.or.YC.NE.YA) CALL XCURDN(XC,YC,1,0)
              XA=XC
              YA=YC
              XC=D(H3,H5,X3,X5)
              YC=D(H3,H5,Y3,Y5)
              IF(XC.NE.XA.or.YC.NE.YA) CALL XCURDN(XC,YC,1,0)
              XA=XC
              YA=YC
          ELSE
              XC=D(H4,H5,X4,X5)
              YC=D(H4,H5,Y4,Y5)
              IF(XC.NE.XA.or.YC.NE.YA) CALL XCURDN(XC,YC,1,0)
              XA=XC
              YA=YC
          ENDIF
      ENDIF
      XB=D(H3,H4,X3,X4)
      YB=D(H3,H4,Y3,Y4)
      I1=I3
      J1=J3
   30 IF (I1.EQ.I10.AND.J1.EQ.J10.AND.I4.EQ.I40.AND.J4.EQ.J40) THEN
        IF(ISWTCH.EQ.1)THEN
              CALL XCURDN(XB,YB,0,1)
        ELSE
            CALL XCURUP(XB,YB)
        ENDIF
      ELSE
        IWRK(I1,J1)=1
        IWRK(I4,J4)=1
        ISW=MOD(ISW-ISA+5,4)+1
        I2=I1+MOD(ISW-1,2)*(1-2*((ISW-1)/2))
        J2=J1+MOD(ISW  ,2)*(1-2*((ISW-1)/2))
        I3=I4+MOD(ISW-1,2)*(1-2*((ISW-1)/2))
        J3=J4+MOD(ISW  ,2)*(1-2*((ISW-1)/2))

        IF (I2.EQ.0.OR.I3.EQ.0.OR.I2.EQ.MGP.OR.I3.EQ.MGP  .OR.
     :     J2.EQ.0.OR.J3.EQ.0.OR.J2.EQ.JGP.OR.J3.EQ.JGP ) THEN
          IF(ISWTCH.EQ.1)THEN
              CALL XCURDN(XB,YB,0,1)
          ELSE
            CALL XCURUP(XB,YB)
          ENDIF
        ELSE
          IF(ZG(I2,J2).EQ.DUM.OR.ZG(I3,J3).EQ.DUM)THEN
            IF(ISWTCH.EQ.1) THEN
              CALL XCURDN(XB,YB,0,1)
            ELSE
              CALL XCURUP(XB,YB)
            ENDIF
          ELSE
            IF(ISWTCH.EQ.1) THEN
                 IF(XB.NE.XA.or.YB.NE.YA) CALL XCURDN(XB,YB,0,1)
            ELSE
              CALL XCURUP(XB,YB)
            ENDIF
          ENDIF
          XA=XB
          YA=YB
          GOTO 101
        END IF
      ENDIF
    2 CONTINUE
      CALL XLPNUP( X(1,1), Y(1,1) )
      RETURN
      END

      SUBROUTINE ZCONTR(ZG,Z,IWRK,MD,MG,JG,CV)
      DIMENSION ZG(MD ,*),Z(MD ,*),IWRK(MG ,*)
C*   The second edition of the contour tracing
C*   Zhang Zuojun, Jan. 1988
C*   New update including contouring on triagle grids
C*   When MODE=0 contouring perform on retangular grids (default)
C*   When MODE=1 contouring perform on triangular grids .
      COMPLEX   Z,B1,B2,ZA,ZB,ZC,Z1,Z2,Z3,Z4,Z5,D
C     IFUN1(K)=K+MG*((MGP-K)/MGP-K/MGP)
      D(P1,P2,B1,B2   )=B1+(CV-P1)*(B2-B1)/(P2-P1)
      CALL ZQCONM(MODE)
      MGP=MG+1
      JGP=JG+1
      DO 4 J=1,JG
      DO 4 I=1,MG
    4 IWRK(I,J)=0
      DO 1 JJ=1,2*(MG+JG-2)
      IF(JJ.LT.MG) THEN
          I4=JJ
          J4=1
          ISW=1
      ELSEIF(JJ.LT.MG+JG-1) THEN
          I4=MG
          J4=JJ-MG+1
          ISW=4
      ELSEIF(JJ.LT.MG+MG+JG-2) THEN
          I4=MG+MG+JG-JJ-1
          J4=JG
          ISW=3
      ELSEIF(JJ.LT.MG+MG+JG+JG-3) THEN
          I4=1
          J4=MG+MG+JG+JG-2-JJ
          ISW=2
      ENDIF
      INI=MOD(ISW  ,2)*(1-2*(MOD(ISW,4)/2))
      INJ=MOD(ISW+1,2)*(1-2*(MOD(ISW,4)/2))
      I1=I4+INI
      J1=J4+INJ
      IF(I1.EQ.0.OR.I1.EQ.MGP.OR.J1.EQ.0.OR.J1.EQ.JGP)GOTO 1
      H1=ZG(I1,J1)
      H4=ZG(I4,J4)
      IF(H1.GE.CV.OR.H4.LT.CV ) GOTO 1
      Z1= Z(I1,J1)
      Z4= Z(I4,J4)
      ZA=D(H4,H1,Z4,Z1)
      CALL XCURUP(REAL(ZA),AIMAG(ZA))
      I2=I1+MOD(ISW-1,2)*(1-2*((ISW-1)/2))
      J2=J1+MOD(ISW  ,2)*(1-2*((ISW-1)/2))
      I3=I4+MOD(ISW-1,2)*(1-2*((ISW-1)/2))
      J3=J4+MOD(ISW  ,2)*(1-2*((ISW-1)/2))
  201 H1=ZG(I1,J1)
      H2=ZG(I2,J2)
      H3=ZG(I3,J3)
      H4=ZG(I4,J4)
      H5=0.25*(H1+H2+H3+H4)
      Z1= Z(I1,J1)
      Z2= Z(I2,J2)
      Z3= Z(I3,J3)
      Z4= Z(I4,J4)
      IF(MODE.EQ.1) Z5=0.25*(Z1+Z2+Z3+Z4)
      IF(H2-CV) 52,53,53
   52 IF(H3-CV) 63,62,62
   53 IF(H3-CV) 54,61,61
   54 IF(H5-CV) 63,61,61
   61 ISA=1
      IF(MODE.EQ.1) THEN
          IF(H5.LT.CV.AND.H3.GE.CV) THEN
              ZC=D(H4,H5,Z4,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),0,0)
              ZA=ZC
              ZC=D(H3,H5,Z3,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),0,0)
              ZA=ZC
              ZC=D(H2,H5,Z2,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),0,0)
              ZA=ZC
          ELSE
              ZC=D(H1,H5,Z1,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),0,0)
              ZA=ZC
          ENDIF
      ENDIF
      ZB=D(H1,H2,Z1,Z2)
      I4=I2
      J4=J2
      GOTO 60
   62 ISA=2
      IF(MODE.EQ.1) THEN
          IF(H5.LT.CV) THEN
              ZC=D(H4,H5,Z4,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),0,0)
              ZA=ZC
              ZC=D(H3,H5,Z3,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),0,0)
              ZA=ZC
          ELSE
              ZC=D(H1,H5,Z1,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),0,0)
              ZA=ZC
              ZC=D(H2,H5,Z2,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),0,0)
              ZA=ZC
          ENDIF
      ENDIF
      ZB=D(H2,H3,Z2,Z3)
      I1=I2
      J1=J2
      I4=I3
      J4=J3
      GOTO 60
   63 ISA=3
      IF(MODE.EQ.1) THEN
          IF(H5.GE.CV.AND.H2.LT.CV) THEN
              ZC=D(H1,H5,Z1,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),0,0)
              ZA=ZC
              ZC=D(H2,H5,Z2,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),0,0)
              ZA=ZC
              ZC=D(H3,H5,Z3,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),0,0)
              ZA=ZC
          ELSE
              ZC=D(H4,H5,Z4,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),0,0)
              ZA=ZC
          ENDIF
      ENDIF
      ZB=D(H3,H4,Z3,Z4)
      I1=I3
      J1=J3
   60 ISW=MOD(ISW-ISA+5,4)+1
      I2=I1+MOD(ISW-1,2)*(1-2*((ISW-1)/2))
      J2=J1+MOD(ISW  ,2)*(1-2*((ISW-1)/2))
      I3=I4+MOD(ISW-1,2)*(1-2*((ISW-1)/2))
      J3=J4+MOD(ISW  ,2)*(1-2*((ISW-1)/2))
      IF( I2.EQ.0.OR.I3.EQ.0.OR.I2.EQ.MGP.OR.I3.EQ.MGP  .OR.
     :    J2.EQ.0.OR.J3.EQ.0.OR.J2.EQ.JGP.OR.J3.EQ.JGP) THEN
          CALL XCURDN(REAL(ZB),AIMAG(ZB),0,1)
      ELSE
          IF(ZB.NE.ZA) CALL XCURDN(REAL(ZB),AIMAG(ZB),0,0)
          ZA=ZB
          IWRK(I1,J1)=1
          IWRK(I4,J4)=1
          GOTO 201
      ENDIF
   1  CONTINUE
      DO 2 J=2,JG-1
      DO 2 I=1,MG-1
      ISW=1
      I10=I+1
      J10=J
      I40=I
      J40=J
      IF(IWRK(I10,J10).EQ.1.AND.IWRK(I40,J40).EQ.1) GOTO 2
      H1=ZG(I10,J10)
      H4=ZG(I40,J40)
      IF(H1.GE.CV.OR.H4.LT.CV ) GOTO 2
      I1=I10
      J1=J10
      I4=I40
      J4=J40
      Z1= Z(I1,J1)
      Z4= Z(I4,J4)
      ZA=D(H4,H1,Z4,Z1)
      CALL XCURUP(REAL(ZA),AIMAG(ZA))
  101 I2=I1+MOD(ISW-1,2)*(1-2*((ISW-1)/2))
      J2=J1+MOD(ISW  ,2)*(1-2*((ISW-1)/2))
      I3=I4+MOD(ISW-1,2)*(1-2*((ISW-1)/2))
      J3=J4+MOD(ISW  ,2)*(1-2*((ISW-1)/2))
      H1=ZG(I1,J1)
      H2=ZG(I2,J2)
      H3=ZG(I3,J3)
      H4=ZG(I4,J4)
      H5=0.25*(H1+H2+H3+H4)
      Z1= Z(I1,J1)
      Z2= Z(I2,J2)
      Z3= Z(I3,J3)
      Z4= Z(I4,J4)
      IF(MODE.EQ.1) Z5=0.25*(Z1+Z2+Z3+Z4)
      IF(H2-CV) 12,13,13
   12 IF(H3-CV) 23,22,22
   13 IF(H3-CV) 14,21,21
   14 IF(H5-CV) 23,21,21
   21 ISA=1
      IF(MODE.EQ.1) THEN
          IF(H5.LT.CV.AND.H3.GE.CV) THEN
              ZC=D(H4,H5,Z4,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),1,0)
              ZA=ZC
              ZC=D(H3,H5,Z3,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),1,0)
              ZA=ZC
              ZC=D(H2,H5,Z2,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),1,0)
              ZA=ZC
          ELSE
              ZC=D(H1,H5,Z1,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),1,0)
              ZA=ZC
          ENDIF
      ENDIF
      ZB=D(H1,H2,Z1,Z2)
      I4=I2
      J4=J2
      GOTO 30
   22 ISA=2
      IF(MODE.EQ.1) THEN
          IF(H5.LT.CV) THEN
              ZC=D(H4,H5,Z4,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),1,0)
              ZA=ZC
              ZC=D(H3,H5,Z3,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),1,0)
              ZA=ZC
          ELSE
              ZC=D(H1,H5,Z1,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),1,0)
              ZA=ZC
              ZC=D(H2,H5,Z2,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),1,0)
              ZA=ZC
          ENDIF
      ENDIF
      ZB=D(H2,H3,Z2,Z3)
      I1=I2
      J1=J2
      I4=I3
      J4=J3
      GOTO 30
   23 ISA=3
      IF(MODE.EQ.1) THEN
          IF(H5.GE.CV.AND.H2.LT.CV) THEN
              ZC=D(H1,H5,Z1,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),1,0)
              ZA=ZC
              ZC=D(H2,H5,Z2,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),1,0)
              ZA=ZC
              ZC=D(H3,H5,Z3,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),1,0)
              ZA=ZC
          ELSE
              ZC=D(H4,H5,Z4,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),1,0)
              ZA=ZC
          ENDIF
      ENDIF
      ZB=D(H3,H4,Z3,Z4)
      I1=I3
      J1=J3
   30 IF( I1.EQ.I10.AND.J1.EQ.J10.AND.I4.EQ.I40.AND.J4.EQ.J40) THEN
          CALL XCURDN(REAL(ZB),AIMAG(ZB),1,1)
      ELSE
          IF(ZB.NE.ZA) CALL XCURDN(REAL(ZB),AIMAG(ZB),1,0)
          ZA=ZB
          IWRK(I1,J1)=1
          IWRK(I4,J4)=1
          ISW=MOD(ISW-ISA+5,4)+1
          GOTO 101
      ENDIF
    2 CONTINUE
      CALL XLPNUP( REAL(Z(1,1)), AIMAG(Z(1,1)) )
      RETURN
      END

      SUBROUTINE ZCONTJ(ZG,Z,IWRK,MD,MG,JG,CV,SPEC)
      DIMENSION ZG(MD ,*),Z(MD ,*),IWRK(MG ,*)
C*   New update for contouring allowing special value holes (SPEC)
C*   The second edition of the contour tracing
C*   Zhang Zuojun, Jan. 1988
C*   New update including contouring on triagle grids
C*   When MODE=0 contouring perform on retangular grids (default)
C*   When MODE=1 contouring perform on triangular grids .
      COMPLEX   Z,B1,B2,ZA,ZB,ZC,Z1,Z2,Z3,Z4,Z5,D
C     IFUN1(K)=K+MG*((MGP-K)/MGP-K/MGP)
      D(P1,P2,B1,B2   )=B1+(CV-P1)*(B2-B1)/(P2-P1)
      CALL ZQCONM(MODE)
      DUM=SPEC
      MGP=MG+1
      JGP=JG+1
      DO 4 J=1,JG
      DO 4 I=1,MG
    4 IWRK(I,J)=0
      DO 1 JJ=1,2*(MG+JG-2)
      IF(JJ.LT.MG) THEN
          I4=JJ
          J4=1
          ISW=1
      ELSEIF(JJ.LT.MG+JG-1) THEN
          I4=MG
          J4=JJ-MG+1
          ISW=4
      ELSEIF(JJ.LT.MG+MG+JG-2) THEN
          I4=MG+MG+JG-JJ-1
          J4=JG
          ISW=3
      ELSEIF(JJ.LT.MG+MG+JG+JG-3) THEN
          I4=1
          J4=MG+MG+JG+JG-2-JJ
          ISW=2
      ENDIF
      INI=MOD(ISW  ,2)*(1-2*(MOD(ISW,4)/2))
      INJ=MOD(ISW+1,2)*(1-2*(MOD(ISW,4)/2))
      I1=I4+INI
      J1=J4+INJ
      IF(I1.EQ.0.OR.I1.EQ.MGP.OR.J1.EQ.0.OR.J1.EQ.JGP)GOTO 1
      H1=ZG(I1,J1)
      H4=ZG(I4,J4)
      IF(H1.GE.CV.OR.H4.LT.CV ) GOTO 1
      Z1= Z(I1,J1)
      Z4= Z(I4,J4)
      ZA=D(H4,H1,Z4,Z1)
      CALL XCURUP(REAL(ZA),AIMAG(ZA))
      I2=I1+MOD(ISW-1,2)*(1-2*((ISW-1)/2))
      J2=J1+MOD(ISW  ,2)*(1-2*((ISW-1)/2))
      I3=I4+MOD(ISW-1,2)*(1-2*((ISW-1)/2))
      J3=J4+MOD(ISW  ,2)*(1-2*((ISW-1)/2))
  201 H1=ZG(I1,J1)
      H2=ZG(I2,J2)
      H3=ZG(I3,J3)
      H4=ZG(I4,J4)
      H5=0.25*(H1+H2+H3+H4)
      IF(H1.EQ.DUM.OR.H2.EQ.DUM.OR.H3.EQ.DUM.OR.H4.EQ.DUM)THEN
        ISWTCH=0
      ELSE
        ISWTCH=1
      ENDIF
      Z1= Z(I1,J1)
      Z2= Z(I2,J2)
      Z3= Z(I3,J3)
      Z4= Z(I4,J4)
      IF(MODE.EQ.1) Z5=0.25*(Z1+Z2+Z3+Z4)
      IF(H2-CV) 52,53,53
   52 IF(H3-CV) 63,62,62
   53 IF(H3-CV) 54,61,61
   54 IF(H5-CV) 63,61,61
   61 ISA=1
      IF(MODE.EQ.1.AND.ISWTCH.EQ.1) THEN
          IF(H5.LT.CV.AND.H3.GE.CV) THEN
              ZC=D(H4,H5,Z4,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),0,0)
              ZA=ZC
              ZC=D(H3,H5,Z3,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),0,0)
              ZA=ZC
              ZC=D(H2,H5,Z2,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),0,0)
              ZA=ZC
          ELSE
              ZC=D(H1,H5,Z1,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),0,0)
              ZA=ZC
          ENDIF
      ENDIF
      ZB=D(H1,H2,Z1,Z2)
      I4=I2
      J4=J2
      GOTO 60
   62 ISA=2
      IF(MODE.EQ.1.AND.ISWTCH.EQ.1) THEN
          IF(H5.LT.CV) THEN
              ZC=D(H4,H5,Z4,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),0,0)
              ZA=ZC
              ZC=D(H3,H5,Z3,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),0,0)
              ZA=ZC
          ELSE
              ZC=D(H1,H5,Z1,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),0,0)
              ZA=ZC
              ZC=D(H2,H5,Z2,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),0,0)
              ZA=ZC
          ENDIF
      ENDIF
      ZB=D(H2,H3,Z2,Z3)
      I1=I2
      J1=J2
      I4=I3
      J4=J3
      GOTO 60
   63 ISA=3
      IF(MODE.EQ.1.AND.ISWTCH.EQ.1) THEN
          IF(H5.GE.CV.AND.H2.LT.CV) THEN
              ZC=D(H1,H5,Z1,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),0,0)
              ZA=ZC
              ZC=D(H2,H5,Z2,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),0,0)
              ZA=ZC
              ZC=D(H3,H5,Z3,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),0,0)
              ZA=ZC
          ELSE
              ZC=D(H4,H5,Z4,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),0,0)
              ZA=ZC
          ENDIF
      ENDIF
      ZB=D(H3,H4,Z3,Z4)
      I1=I3
      J1=J3
   60 ISW=MOD(ISW-ISA+5,4)+1
      I2=I1+MOD(ISW-1,2)*(1-2*((ISW-1)/2))
      J2=J1+MOD(ISW  ,2)*(1-2*((ISW-1)/2))
      I3=I4+MOD(ISW-1,2)*(1-2*((ISW-1)/2))
      J3=J4+MOD(ISW  ,2)*(1-2*((ISW-1)/2))
      IF( I2.EQ.0.OR.I3.EQ.0.OR.I2.EQ.MGP.OR.I3.EQ.MGP  .OR.
     :    J2.EQ.0.OR.J3.EQ.0.OR.J2.EQ.JGP.OR.J3.EQ.JGP )THEN
        IF(ISWTCH.EQ.1) THEN
              CALL XCURDN(REAL(ZB),AIMAG(ZB),0,1)
        ELSE
            CALL XCURUP(REAL(ZB),AIMAG(ZB))
        ENDIF
      ELSE
        IF(ISWTCH.EQ.1.AND.ZB.NE.ZA)THEN
            IF(ZG(I2,J2).EQ.DUM.OR.ZG(I3,J3).EQ.DUM)THEN
                  CALL XCURDN(REAL(ZB),AIMAG(ZB),0,0)
            ELSE
                  CALL XCURDN(REAL(ZB),AIMAG(ZB),0,1)
            ENDIF
        ELSE
            CALL XCURUP(REAL(ZB),AIMAG(ZB))
        ENDIF
          ZA=ZB
          IWRK(I1,J1)=1
          IWRK(I4,J4)=1
          GOTO 201
      ENDIF
   1  CONTINUE
      DO 2 J=2,JG-1
      DO 2 I=1,MG-1
      ISW=1
      I10=I+1
      J10=J
      I40=I
      J40=J
      IF(IWRK(I10,J10).EQ.1.AND.IWRK(I40,J40).EQ.1) GOTO 2
      H1=ZG(I10,J10)
      H4=ZG(I40,J40)
      IF(H1.GE.CV.OR.H4.LT.CV ) GOTO 2
      I1=I10
      J1=J10
      I4=I40
      J4=J40
      Z1= Z(I1,J1)
      Z4= Z(I4,J4)
      ZA=D(H4,H1,Z4,Z1)
      CALL XCURUP(REAL(ZA),AIMAG(ZA))
      I2=I1+MOD(ISW-1,2)*(1-2*((ISW-1)/2))
      J2=J1+MOD(ISW  ,2)*(1-2*((ISW-1)/2))
      I3=I4+MOD(ISW-1,2)*(1-2*((ISW-1)/2))
      J3=J4+MOD(ISW  ,2)*(1-2*((ISW-1)/2))
  101 H1=ZG(I1,J1)
      H2=ZG(I2,J2)
      H3=ZG(I3,J3)
      H4=ZG(I4,J4)
      H5=0.25*(H1+H2+H3+H4)
      IF(H1.EQ.DUM.OR.H2.EQ.DUM.OR.H3.EQ.DUM.OR.H4.EQ.DUM)THEN
        ISWTCH=0
      ELSE
        ISWTCH=1
      ENDIF
      Z1= Z(I1,J1)
      Z2= Z(I2,J2)
      Z3= Z(I3,J3)
      Z4= Z(I4,J4)
      IF(MODE.EQ.1) Z5=0.25*(Z1+Z2+Z3+Z4)
      IF(H2-CV) 12,13,13
   12 IF(H3-CV) 23,22,22
   13 IF(H3-CV) 14,21,21
   14 IF(H5-CV) 23,21,21
   21 ISA=1
      IF(MODE.EQ.1.AND.ISWTCH.EQ.1) THEN
          IF(H5.LT.CV.AND.H3.GE.CV) THEN
              ZC=D(H4,H5,Z4,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),1,0)
              ZA=ZC
              ZC=D(H3,H5,Z3,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),1,0)
              ZA=ZC
              ZC=D(H2,H5,Z2,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),1,0)
              ZA=ZC
          ELSE
              ZC=D(H1,H5,Z1,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),1,0)
              ZA=ZC
          ENDIF
      ENDIF
      ZB=D(H1,H2,Z1,Z2)
      I4=I2
      J4=J2
      GOTO 30
   22 ISA=2
      IF(MODE.EQ.1.AND.ISWTCH.EQ.1) THEN
          IF(H5.LT.CV) THEN
              ZC=D(H4,H5,Z4,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),1,0)
              ZA=ZC
              ZC=D(H3,H5,Z3,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),1,0)
              ZA=ZC
          ELSE
              ZC=D(H1,H5,Z1,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),1,0)
              ZA=ZC
              ZC=D(H2,H5,Z2,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),1,0)
              ZA=ZC
          ENDIF
      ENDIF
      ZB=D(H2,H3,Z2,Z3)
      I1=I2
      J1=J2
      I4=I3
      J4=J3
      GOTO 30
   23 ISA=3
      IF(MODE.EQ.1.AND.ISWTCH.EQ.1) THEN
          IF(H5.GE.CV.AND.H2.LT.CV) THEN
              ZC=D(H1,H5,Z1,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),1,0)
              ZA=ZC
              ZC=D(H2,H5,Z2,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),1,0)
              ZA=ZC
              ZC=D(H3,H5,Z3,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),1,0)
              ZA=ZC
          ELSE
              ZC=D(H4,H5,Z4,Z5)
              IF(ZC.NE.ZA) CALL XCURDN(REAL(ZC),AIMAG(ZC),1,0)
              ZA=ZC
          ENDIF
      ENDIF
      ZB=D(H3,H4,Z3,Z4)
      I1=I3
      J1=J3
   30 IF( I1.EQ.I10.AND.J1.EQ.J10.AND.I4.EQ.I40.AND.J4.EQ.J40) THEN
        IF(ISWTCH.EQ.1)THEN
              CALL XCURDN(REAL(ZB),AIMAG(ZB),0,1)
        ELSE
            CALL XCURUP(REAL(ZB),AIMAG(ZB))
        ENDIF
      ELSE
          IWRK(I1,J1)=1
          IWRK(I4,J4)=1
          ISW=MOD(ISW-ISA+5,4)+1
          I2=I1+MOD(ISW-1,2)*(1-2*((ISW-1)/2))
          J2=J1+MOD(ISW  ,2)*(1-2*((ISW-1)/2))
          I3=I4+MOD(ISW-1,2)*(1-2*((ISW-1)/2))
          J3=J4+MOD(ISW  ,2)*(1-2*((ISW-1)/2))
        IF(ZG(I2,J2).EQ.DUM.OR.ZG(I3,J3).EQ.DUM)THEN
           IF(ISWTCH.EQ.1) THEN
             CALL XCURDN(REAL(ZB),AIMAG(ZB),0,1)
           ELSE
             CALL XCURUP(REAL(ZB),AIMAG(ZB))
           ENDIF
        ELSE
           IF(ISWTCH.EQ.1) THEN
                 IF(ZB.NE.ZA) CALL XCURDN(REAL(ZB),AIMAG(ZB),0,1)
           ELSE
             CALL XCURUP(REAL(ZB),AIMAG(ZB))
           ENDIF
        ENDIF
          ZA=ZB
          GOTO 101
      ENDIF
    2 CONTINUE
      CALL XLPNUP( REAL(Z(1,1)), AIMAG(Z(1,1)) )
      RETURN
      END

      SUBROUTINE ZCONTM(MODES)
      SAVE MODE
      MODE=MODES
      RETURN

C***************
      ENTRY      ZQCONM(MODET)
      MODET=MODE
      RETURN

      DATA MODE/0/
      END
      SUBROUTINE ZRCNTR(ZG,Z,MD,MG,NG,CV)
C*     The final edition of the contouring routine
C*                                 12.6.1987
      DIMENSION ZG(MD,*),Z(MD,*)
      COMPLEX   Z,B1,B2,ZA,ZB,Z1,Z2,Z3,Z4,D
                D(P1,P2,B1,B2   )=B1+(CV-P1)*(B2-B1)/(P2-P1)
      JG=ABS(NG)
      IDUB=0
      DO 1 J=1,MG-1
      DO 2 I=1,JG-1
      H1=ZG(J  ,I  )
      H2=ZG(J+1,I  )
      Z1= Z(J  ,I  )
      Z2= Z(J+1,I  )
      H3=ZG(J+1,I+1)
      H4=ZG(  J,I+1)
      Z3= Z(J+1,I+1)
      Z4= Z(  J,I+1)
      IF(H1-CV)11,20,20
   11 IF(H2-CV)12,14,14
   12 IF(H3-CV)13,15,15
   13 IF(H4-CV) 2,30,30
   14 IF(H3-CV)16,17,17
   15 IF(H4-CV)31,32,32
   16 IF(H4-CV)33,34,34
   17 IF(H4-CV)35,36,36
   20 IF(H2-CV)21,23,23
   21 IF(H3-CV)22,24,24
   22 IF(H4-CV)36,35,35
   23 IF(H3-CV)25,26,26
   24 IF(H4-CV)37,33,33
   25 IF(H4-CV)32,31,31
   26 IF(H4-CV)30, 2, 2
   30 ZA=D(H1,H4,Z1,Z4)
      ZB=D(H4,H3,Z4,Z3)
      GOTO 40
   31 ZA=D(H2,H3,Z2,Z3)
      ZB=D(H4,H3,Z4,Z3)
      GOTO 40
   32 ZA=D(H1,H4,Z1,Z4)
      ZB=D(H2,H3,Z2,Z3)
      GOTO 40
   33 ZA=D(H1,H2,Z1,Z2)
      ZB=D(H2,H3,Z2,Z3)
      IDUB=0
      GOTO 40
   34 IDUB=1
      H5=0.25*(H1+H2+H3+H4)
      IF(H5.GT.CV) IDUB=-1
      GOTO (30,31) 2-(1+IDUB)/2
   35 ZA=D(H1,H2,Z1,Z2)
      ZB=D(H4,H3,Z4,Z3)
      GOTO 40
   36 ZA=D(H1,H2,Z1,Z2)
      ZB=D(H1,H4,Z1,Z4)
      IDUB=0
      GOTO 40
   37 IDUB=-1
      H5=0.25*(H1+H2+H3+H4)
      IF(H5.GT.CV) IDUB=1
      GOTO (30,31) 2-(1+IDUB)/2
      GOTO 31
   40 CONTINUE
      IF(NG.GT.0) THEN
         CALL XPENUP(REAL(ZA),AIMAG(ZA))
         CALL XPENDN(REAL(ZB),AIMAG(ZB))
      ELSE
         ZB=ZA+0.7*(ZB-ZA)
         CALL XPENUP(REAL(ZA),AIMAG(ZA))
         CALL XPENDN(REAL(ZB),AIMAG(ZB))
      ENDIF
      IF(IDUB)36,2,33
    2 CONTINUE
    1 CONTINUE
      RETURN
      END

      SUBROUTINE ZRCNTA(ZG,Z,MD,MG,JG,CVL,NC)
      DIMENSION ZG(MD,*),Z(MD,*),CVL(*)
      COMPLEX   Z
      DO 100 K=1,NC
      CV=CVL(K)
  100 CALL ZRCNTR(ZG,Z,MD,MG,JG,CV)
      RETURN
      END

      SUBROUTINE ZRCNTB(ZG,Z,MD,MG,NG,CVL,NC)
C Contouring on triangular grid
      DIMENSION ZG(MD,*),Z(MD,*),CVL(*)
      COMPLEX   Z,B1,B2,ZA,ZB,Z1,Z2,Z3,D
      D(P1,P2,B1,B2   )=B1+(CV-P1)*(B2-B1)/(P2-P1)
      JG=ABS(NG)
      DO  50 I=1,MG-1
      DO  50 J=1,JG-1
      H3=0.25*(ZG(I,J)+ZG(I,J+1)+ZG(I+1,J+1)+ZG(I+1,J))
      Z3=0.25*( Z(I,J)+ Z(I,J+1)+ Z(I+1,J+1)+ Z(I+1,J))
      DO  50 M=1,4
      H1=ZG(I+MOD(M  ,4)/2,J+MOD(M-1,4)/2)
      Z1=Z (I+MOD(M  ,4)/2,J+MOD(M-1,4)/2)
      H2=ZG(I+MOD(M+1,4)/2,J+MOD(M  ,4)/2)
      Z2=Z (I+MOD(M+1,4)/2,J+MOD(M  ,4)/2)
      DO  50 K=1,NC
      CV=CVL(K)
      IF(H1-CV) 1, 2, 2
    1 IF(H2-CV) 3, 4, 4
    3 IF(H3-CV)50,30,30
    4 IF(H3-CV)20,10,10
    2 IF(H2-CV) 5, 6, 6
    5 IF(H3-CV)10,20,20
    6 IF(H3-CV)30,50,50
   10 ZA=D(H3,H1,Z3,Z1)
      ZB=D(H1,H2,Z1,Z2)
      GOTO 40
   20 ZA=D(H1,H2,Z1,Z2)
      ZB=D(H2,H3,Z2,Z3)
      GOTO 40
   30 ZA=D(H2,H3,Z2,Z3)
      ZB=D(H3,H1,Z3,Z1)
   40 IF(ZA.EQ.ZB) GOTO 50
      IF(NG.GT.0) THEN
         CALL XPENUP(REAL(ZA),AIMAG(ZA))
         CALL XPENDN(REAL(ZB),AIMAG(ZB))
      ELSE
         ZB=ZA+0.7*(ZB-ZA)
         CALL XPENUP(REAL(ZA),AIMAG(ZA))
         CALL XPENDN(REAL(ZB),AIMAG(ZB))
      ENDIF
   50 CONTINUE
      RETURN
      END

      SUBROUTINE XHATCH(Z,X,Y,MD,M,N,CL1,CL2, MODE)
      REAL X(MD,*),Y(MD,*),Z(MD,*)
      IF( MODE.EQ.0) THEN
        CALL XHATCX(Z,X,Y,MD,M,N,CL1,CL2)
        CALL XHATCY(Z,X,Y,MD,M,N,CL1,CL2)
        RETURN
      ENDIF
      IF( MODE.EQ.1)   CALL XHATCX(Z,X,Y,MD,M,N,CL1,CL2)
      IF( MODE.EQ.-1)  CALL XHATCY(Z,X,Y,MD,M,N,CL1,CL2)
      RETURN
      END

      SUBROUTINE XHATCX(Z,X,Y,MD,M,N,CL1,CL2)
      REAL X(MD,*),Y(MD,*), Z(MD,*) ,XP(10),YP(10)
      COMMON /XFTR06/ XFACTR,YFACTR
      COMMON /XHCH35/ DH
      D(P1,P2,B1,B2)=B1+(CV-P1)*(B2-B1)/(P2-P1)
      IF( CL1.EQ.CL2)RETURN
C
      MM=M-1
      NM=N-1
      XSIGN=SIGN( 1.0, X(M,1)-X(1,1) )
      DXP = DH*XSIGN
      DX=DXP/XFACTR
      XS= X(1,1)
      IS=1
  1   CONTINUE
      XS=XS+DX
      IF((XS-X(M,1))*XSIGN. GT. 0.0) RETURN
  2   IF((XS-X(IS+1,1))*XSIGN.GT.0.0) THEN
        IS=IS+1
        IF(IS.GT.M-1)RETURN
        GOTO 2
      ENDIF
      HS=Z(IS,1)+(Z(IS+1,1)-Z(IS,1))/(X(IS+1,1)-X(IS,1))*(XS-X(IS,1))
      IF( HS.GE.CL1.AND.HS.LE.CL2) THEN
         YPP= Y(IS,1)+(Y(IS+1,1)-Y(IS,1))
     :   /(X(IS+1,1)-X(IS,1))*(XS-X(IS,1))
         CALL XPENUP(XS,YPP)
         MODEP =1
      ELSE
         MODEP =0
      ENDIF
      IP=IS
      DO 540 J=1,NM
      X1=X(IP,J)
      Y1=Y(IP,J)
      H1=Z(IP,J)
      X4=X(IP+1,J)
      Y4=Y(IP+1,J)
      H4=Z(IP+1,J)
      X2=X(IP,J+1)
      Y2=Y(IP,J+1)
      H2=Z(IP,J+1)
      X3=X(IP+1,J+1)
      Y3=Y(IP+1,J+1)
      H3=Z(IP+1,J+1)
      CV=CL1
      NP=0
      DO 300 IK=1,2
      IDUB=0
      IF(H1-CV)11,20,20
 11   IF(H2-CV)12,14,14
 12   IF(H3-CV)13,15,15
 13   IF(H4-CV)250,250,30
 14   IF(H3-CV)16,17,17
 15   IF(H4-CV)31,32,32
 16   IF(H4-CV)33,39,39
 17   IF(H4-CV)35,36,36
 20   IF(H2-CV)21,23,23
 21   IF(H3-CV)22,24,24
 22   IF(H4-CV)36,35,35
 23   IF(H3-CV)25,26,26
 24   IF(H4-CV)38,33,33
 25   IF(H4-CV)32,31,31
 26   IF(H4-CV)30,250,250
 38   IF((H1+H2+H3+H4)*0.25-CV) 34,34,37
 39   IF((H1+H2+H3+H4)*0.25-CV) 37,37,34
  30  XA=D(H1,H4,X1,X4)
      YA=D(H1,H4,Y1,Y4)
      XB=D(H4,H3,X4,X3)
      YB=D(H4,H3,Y4,Y3)
      GOTO 40
  31  XA=D(H2,H3,X2,X3)
      YA=D(H2,H3,Y2,Y3)
      XB=D(H4,H3,X4,X3)
      YB=D(H4,H3,Y4,Y3)
      GOTO 40
  32  XA=D(H1,H4,X1,X4)
      YA=D(H1,H4,Y1,Y4)
      XB=D(H2,H3,X2,X3)
      YB=D(H2,H3,Y2,Y3)
      GOTO 40
  33  XA=D(H1,H2,X1,X2)
      YA=D(H1,H2,Y1,Y2)
      XB=D(H2,H3,X2,X3)
      YB=D(H2,H3,Y2,Y3)
      IDUB=0
      GOTO 40
  34  IDUB=1
      GOTO 31
  35  XA=D(H1,H2,X1,X2)
      YA=D(H1,H2,Y1,Y2)
      XB=D(H4,H3,X4,X3)
      YB=D(H4,H3,Y4,Y3)
      GOTO 40
  36  XA=D(H1,H2,X1,X2)
      YA=D(H1,H2,Y1,Y2)
      XB=D(H1,H4,X1,X4)
      YB=D(H1,H4,Y1,Y4)
      IDUB=0
      GOTO 40
  37  IDUB=-1
      GOTO 30
  40  CONTINUE
      IF(XA.EQ.XB) GOTO 245
      XC=MIN(XA,XB)
      XD=MAX(XA,XB)
      IF( XS.GT.XC.AND.XS.LE.XD) THEN
        NP=NP+1
        XP(NP)=XS
        YP(NP)=YA+(YB-YA)/(XB-XA)*(XS-XA)
      ENDIF
 245  IF(IDUB)33,250,36
 250  CV=CL2
 300  CONTINUE

      IF( NP.GT.2) CALL XHAT01(YP,XP,NP)
      DO 350 NPL=1,NP
        IF( MODEP.EQ.0) THEN
          CALL XPENUP(XP(NPL),YP(NPL))
          MODEP=1
        ELSE
          CALL XPENDN(XP(NPL),YP(NPL))
          MODEP=0
        ENDIF
 350  CONTINUE
 540  CONTINUE
      IF( MODEP.EQ.1) THEN
        YPP=Y(IS,N)+(Y(IS+1,N)-Y(IS,N))/(X(IS+1,N)-X(IS,N))*(XS-X(IS,N))
        CALL XPENDN(XS,YPP)
        MODP=0
      ENDIF
      IF( IS.LE.M-1)GOTO 1
      RETURN
      END

      SUBROUTINE XHATCY(Z,X,Y,MD,M,N,CL1,CL2)
      REAL X(MD,*),Y(MD,*), Z(MD,*) ,XP(10),YP(10)
      COMMON /XFTR06/ XFACTR,YFACTR
      COMMON /XHCH35/ DH
      D(P1,P2,B1,B2)=B1+(CV-P1)*(B2-B1)/(P2-P1)
      IF( CL1.EQ.CL2)RETURN
C
      MM=M-1
      NM=N-1
      YSIGN=SIGN( 1.0, Y(1,N)-Y(1,1))
      DYP = DH*YSIGN
      DY=DYP/YFACTR
      YS= Y(1,1)
      JS=1
  1   CONTINUE
      YS=YS+DY
      IF((YS-Y(1,N))*YSIGN. GT. 0.0) RETURN
  2   IF((YS-Y(1,JS+1  ))*YSIGN.GT.0.0) THEN
        JS=JS+1
        IF( JS.GT. N-1      ) RETURN
        GOTO 2
      ENDIF
      HS=Z(1,JS)+(Z(1,JS+1)-Z(1,JS))/(Y(1,JS+1)-Y(1,JS))*(YS-Y(1,JS))
      IF( HS.GE.CL1.AND.HS.LE.CL2) THEN
      XPP=X(1,JS)+(X(1,JS+1)-X(1,JS))/(Y(1,JS+1)-Y(1,JS))*(YS-Y(1,JS))
         CALL XPENUP(XPP,YS)
         MODEP =1
      ELSE
         MODEP =0
      ENDIF
      DO 540 I=1,M-1
      JP=JS
      X1=X(I,JP)
      Y1=Y(I,JP)
      H1=Z(I,JP)
      X4=X(I+1,JP)
      Y4=Y(I+1,JP)
      H4=Z(I+1,JP)
      X2=X(I,JP+1)
      Y2=Y(I,JP+1)
      H2=Z(I,JP+1)
      X3=X(I+1,JP+1)
      Y3=Y(I+1,JP+1)
      H3=Z(I+1,JP+1)
      CV=CL1
      NP=0
      DO 300 IK=1,2
      IDUB=0
 10   IF(H1-CV)11,20,20
 11   IF(H2-CV)12,14,14
 12   IF(H3-CV)13,15,15
 13   IF(H4-CV)250,250,30
 14   IF(H3-CV)16,17,17
 15   IF(H4-CV)31,32,32
 16   IF(H4-CV)33,39,39
 17   IF(H4-CV)35,36,36
 20   IF(H2-CV)21,23,23
 21   IF(H3-CV)22,24,24
 22   IF(H4-CV)36,35,35
 23   IF(H3-CV)25,26,26
 24   IF(H4-CV)38,33,33
 25   IF(H4-CV)32,31,31
 26   IF(H4-CV)30,250,250
 38   IF((H1+H2+H3+H4)*0.25-CV) 34,34,37
 39   IF((H1+H2+H3+H4)*0.25-CV) 37,37,34
  30  XA=D(H1,H4,X1,X4)
      YA=D(H1,H4,Y1,Y4)
      XB=D(H4,H3,X4,X3)
      YB=D(H4,H3,Y4,Y3)
      GOTO 40
  31  XA=D(H2,H3,X2,X3)
      YA=D(H2,H3,Y2,Y3)
      XB=D(H4,H3,X4,X3)
      YB=D(H4,H3,Y4,Y3)
      GOTO 40
  32  XA=D(H1,H4,X1,X4)
      YA=D(H1,H4,Y1,Y4)
      XB=D(H2,H3,X2,X3)
      YB=D(H2,H3,Y2,Y3)
      GOTO 40
  33  XA=D(H1,H2,X1,X2)
      YA=D(H1,H2,Y1,Y2)
      XB=D(H2,H3,X2,X3)
      YB=D(H2,H3,Y2,Y3)
      IDUB=0
      GOTO 40
  34  IDUB=1
      GOTO 31
  35  XA=D(H1,H2,X1,X2)
      YA=D(H1,H2,Y1,Y2)
      XB=D(H4,H3,X4,X3)
      YB=D(H4,H3,Y4,Y3)
      GOTO 40
  36  XA=D(H1,H2,X1,X2)
      YA=D(H1,H2,Y1,Y2)
      XB=D(H1,H4,X1,X4)
      YB=D(H1,H4,Y1,Y4)
      IDUB=0
      GOTO 40
  37  IDUB=-1
      GOTO 30
  40  CONTINUE
 50   IF(YA.EQ.YB) GOTO 245
      YC=MIN(YA,YB)
      YD=MAX(YA,YB)
      IF( YS.GT.YC.AND.YS.LE.YD) THEN
        NP=NP+1
        YP(NP)=YS
        XP(NP)=XA+(XB-XA)/(YB-YA)*(YS-YA)
      ENDIF
 245  IF(IDUB)33,250,36
 250  CV=CL2
 300  CONTINUE

      IF( NP.GT.2) CALL XHAT01(XP,YP,NP)
      DO 350 NPL=1,NP
        IF( MODEP.EQ.0) THEN
          CALL XPENUP(XP(NPL),YP(NPL))
          MODEP=1
        ELSE
          CALL XPENDN(XP(NPL),YP(NPL))
          MODEP=0
        ENDIF
 350  CONTINUE
 540  CONTINUE
      IF( MODEP.EQ.1) THEN
      XPP=X(M,JS)+(X(M,1+JS)-X(M,JS))/(Y(M,1+JS)-Y(M,JS))*(YS-Y(M,JS))
        CALL XPENDN(XPP,YS)
        MODEP=0
      ENDIF
      IF( JS.LE.(N-1)) GOTO 1
      RETURN
      END

      SUBROUTINE XDHTCH( DD )
      COMMON /XHCH35/ DH
      DH=DD
      RETURN
      END

      SUBROUTINE XHAT01(XP,YP,N)
C arrange data sequence in xp and yp in ascending order of XP element
      REAL XP(N), YP(N), XPT(20),YPT(20)
      INTEGER IN(20)


      IF( n.gt. 20 ) then
        print*,'work arrray defined in XHAT01 not big enough.'
        print*,'plotting job stopped'
        stop
      endif

      XPT(1)=XP(1)
      YPT(1)=YP(1)
      DO 4 I=1,N
 4      IN(I)=0
      K=1
      KK=1
 3    DO 1 I=1,N
        IF( XP(I).LT.XPT(K).AND.IN(I).EQ.0) THEN
          KK=I
          XPT(K)=XP(I)
          YPT(K)=YP(I)
        ENDIF
  1   CONTINUE
      IN(KK)=1
      K=K+1
      IF(K.GT.N) GOTO 6
      DO 5 I=1,N
        IF( IN(I).EQ.0) THEN
          YPT(K)=YP(I)
          XPT(K)=XP(I)
          KK=I
          GOTO 3
        ENDIF
 5    CONTINUE
      GOTO 3
 6    DO 2 I=1,N
        XP(I)=XPT(I)
 2      YP(I)=YPT(I)
      END

      SUBROUTINE XHATCHA(Z,X,Y,xwk,ywk,MD,M,N,CL1a,CL2a,
     :           hatch_angle)
c
c This routine does hatching of arbitary orientation
c between two contour values
c Written Oct 13, 1998 by Ming Xue
c
      implicit none
c
c Input through argument list
c
      integer md,m,n
      REAL X(MD,*),Y(MD,*), Z(MD,*)
      real cl1a,cl2a,cl1,cl2

      REAL Xwk(MD,*),Ywk(MD,*)
c
c Input through common blocks
c
      real XFACTR,YFACTR,DH
      COMMON /XFTR06/ XFACTR,YFACTR
      COMMON /XHCH35/ DH

      real hatch_angle,sinhagl,coshagl
c     common /xhatch_angle/ hatch_angle,sinhagl,coshagl
      real sinhagla,coshagla

c
c Miselaneous local variables
c
      integer i,j,np,idub,npl,ik,i1,j1,i2,j2
      integer ledgefound,uedgefound
      real XP(20),YP(20) ! Work arrays
      real d,p1,p2,b1,b2,cv
      real dxp,dx,hs,xs,xc,xd,tem
      real x1,y1,h1,x2,y2,h2,x3,y3,h3,x4,y4,h4,xa,xb,ya,yb
      real xwkmin,xwkmax,ywkmin,ywkmax,xmin,xmax
      real ys1,ys2
      integer iedge1,jedge1,iedge2,jedge2
      integer iedge3,jedge3,iedge4,jedge4
      real fxtrns,fytrns,fxorig,fyorig
      real xtrns,ytrns,xorig,yorig
      integer i3,j3,i4,j4
      real ys1_new,ys2_new,yxratio
c
c Inline functions
c
      D(P1,P2,B1,B2)=B1+(CV-P1)*(B2-B1)/(P2-P1)

      fxtrns(xorig,yorig)= xorig*coshagl+yorig*sinhagl
      fytrns(xorig,yorig)=-xorig*sinhagl+yorig*coshagl

      fxorig(xtrns,ytrns)= xtrns*coshagl-ytrns*sinhagl
      fyorig(xtrns,ytrns)= xtrns*sinhagl+ytrns*coshagl
c
c Start of executable statements
c
      IF(CL1a.EQ.CL2a) RETURN ! Then there is nothing to do.
      cl1=min(cl1a,cl2a)
      cl2=max(cl1a,cl2a)

      tem = (hatch_angle-90.0)*atan(1.0)/45.0
      sinhagla= sin( tem )
      coshagla = cos( tem )

      yxratio = yfactr/xfactr
      tem=1.0/sqrt(coshagla**2+(yxratio*sinhagla)**2)
      sinhagl = yxratio*sinhagla*tem
      coshagl = coshagla*tem
c
c Transform into a rotated coordinate. Which angle between
c the new y axis and the old x axis is hatch_angle
C
      DO i=1,m
        DO j=1,n
          xwk(i,j)=fxtrns(x(i,j),y(i,j))
          ywk(i,j)=fytrns(x(i,j),y(i,j))
        ENDDO
      ENDDO

      xwkmin = xwk(1,1)
      xwkmax = xwk(1,1)
      ywkmin = ywk(1,1)
      ywkmax = ywk(1,1)
      DO i=1,m
        DO j=1,n
          xwkmin=min(xwk(i,j),xwkmin)
          xwkmax=max(xwk(i,j),xwkmax)
          ywkmin=min(ywk(i,j),ywkmin)
          ywkmax=max(ywk(i,j),ywkmax)
        ENDDO
      ENDDO

      DXP = DH
      DX=DXP/XFACTR
      XS= Xwkmin+dx

100   CONTINUE ! Come back for another hatching line

c     call xpenup(fxorig(xs,ywkmin),fyorig(xs,ywkmin))
c     call xpendn(fxorig(xs,ywkmax),fyorig(xs,ywkmax))

C Scan boxes on the edges, and found the edge of the box with the
c smallest intercepting y with the hatch line.
c
      ledgefound = 0
      ys1=ywkmax
      uedgefound = 0
      ys2=ywkmin

      DO i=1,m-1
      DO j=1,n-1
        IF(i.eq.1.or.i.eq.m-1.or.j.eq.1.or.j.eq.n-1)then

        xmin=min(xwk(i,j),xwk(i+1,j),xwk(i,j+1),xwk(i+1,j+1))
        xmax=max(xwk(i,j),xwk(i+1,j),xwk(i,j+1),xwk(i+1,j+1))
        if(xs.ge.xmin.and.xs.lt.xmax) then

          call getlnsgmnt(xwk,ywk,md,i,j,xs,ys1_new,i1,j1,i2,j2,
     :                    ys2_new,i3,j3,i4,j4)

          if(ledgefound.eq.0.or.
     :      (ledgefound.ne.0.and.ys1_new.lt.ys1))then
            iedge1=i1
            jedge1=j1
            iedge2=i2
            jedge2=j2
            ys1=ys1_new
            ledgefound=1
          endif

          if(uedgefound.eq.0.or.
     :      (uedgefound.ne.0.and.ys2_new.gt.ys2))then
            iedge3=i3
            jedge3=j3
            iedge4=i4
            jedge4=j4
            ys2=ys2_new
            uedgefound=1
          endif

        endif

        endif
      ENDDO
      ENDDO

c starting or ending box not found. No very likely, just in case
      IF(ledgefound.eq.0.or.uedgefound.eq.0) RETURN

      HS=Z(iedge1,jedge1)+(Z(iedge2,jedge2)-Z(iedge1,jedge1))/
     :   (Xwk(iedge2,jedge2)-Xwk(iedge1,jedge1))
     :   *(XS-Xwk(iedge1,jedge1))

      NP=0

      IF( HS.GE.CL1.AND.HS.LE.CL2) THEN
        np=np+1
        xp(np)=XS
        yp(np)=ys1
      ENDIF

      HS=Z(iedge3,jedge3)+(Z(iedge4,jedge4)-Z(iedge3,jedge3))/
     :   (Xwk(iedge4,jedge4)-Xwk(iedge3,jedge3))
     :   *(XS-Xwk(iedge3,jedge3))

      IF( HS.GE.CL1.AND.HS.LE.CL2) THEN
        np=np+1
        xp(np)=XS
        yp(np)=ys2
      ENDIF

      DO 540 I=1,m-1
      DO 540 J=1,N-1

        X1=Xwk(I,J)    ! low-left
        Y1=Ywk(I,J)
        H1=Z(I,J)

        X4=Xwk(I+1,J)  ! low-right
        Y4=Ywk(I+1,J)
        H4=Z(I+1,J)

        X2=Xwk(I,J+1)  ! upper-left
        Y2=Ywk(I,J+1)
        H2=Z(I,J+1)

        X3=Xwk(I+1,J+1) ! upper-right
        Y3=Ywk(I+1,J+1)
        H3=Z(I+1,J+1)

        xmin=min(x1,x2,x3,x4)
        xmax=max(x1,x2,x3,x4)
        if(.not.(xs.ge.xmin.and.xs.lt.xmax)) goto 540

        CV=CL1

        DO 300 IK=1,2  ! Test CL1 and CL2, hence 2 here

          IDUB=0

 10       IF(H1-CV)11,20,20
 11       IF(H2-CV)12,14,14
 12       IF(H3-CV)13,15,15
 13       IF(H4-CV)250,250,30
 14       IF(H3-CV)16,17,17
 15       IF(H4-CV)31,32,32
 16       IF(H4-CV)33,39,39
 17       IF(H4-CV)35,36,36
 20       IF(H2-CV)21,23,23
 21       IF(H3-CV)22,24,24
 22       IF(H4-CV)36,35,35
 23       IF(H3-CV)25,26,26
 24       IF(H4-CV)38,33,33
 25       IF(H4-CV)32,31,31
 26       IF(H4-CV)30,250,250
 38       IF((H1+H2+H3+H4)*0.25-CV) 34,34,37
 39       IF((H1+H2+H3+H4)*0.25-CV) 37,37,34
  30      XA=D(H1,H4,X1,X4)
          YA=D(H1,H4,Y1,Y4)
          XB=D(H4,H3,X4,X3)
          YB=D(H4,H3,Y4,Y3)
          GOTO 40
  31      XA=D(H2,H3,X2,X3)
          YA=D(H2,H3,Y2,Y3)
          XB=D(H4,H3,X4,X3)
          YB=D(H4,H3,Y4,Y3)
          GOTO 40
  32      XA=D(H1,H4,X1,X4)
          YA=D(H1,H4,Y1,Y4)
          XB=D(H2,H3,X2,X3)
          YB=D(H2,H3,Y2,Y3)
          GOTO 40
  33      XA=D(H1,H2,X1,X2)
          YA=D(H1,H2,Y1,Y2)
          XB=D(H2,H3,X2,X3)
          YB=D(H2,H3,Y2,Y3)
          IDUB=0
          GOTO 40
  34      IDUB=1
          GOTO 31
  35      XA=D(H1,H2,X1,X2)
          YA=D(H1,H2,Y1,Y2)
          XB=D(H4,H3,X4,X3)
          YB=D(H4,H3,Y4,Y3)
          GOTO 40
  36      XA=D(H1,H2,X1,X2)
          YA=D(H1,H2,Y1,Y2)
          XB=D(H1,H4,X1,X4)
          YB=D(H1,H4,Y1,Y4)
          IDUB=0
          GOTO 40
  37      IDUB=-1
          GOTO 30
  40      CONTINUE
 50       IF(XA.EQ.XB) GOTO 245
          XC=MIN(XA,XB)
          XD=MAX(XA,XB)
          IF( XS.GT.XC.AND.XS.LE.XD) THEN
            NP=NP+1
            XP(NP)=XS
            YP(NP)=YA+(YB-YA)/(XB-XA)*(XS-XA)
          ENDIF
 245      IF(IDUB)33,250,36
 250      CV=CL2
 300    CONTINUE

 540  CONTINUE

      IF( NP.GE.2) then
        CALL Xsortxyp(xp,yp,NP)
        do npl=1,np,2
          CALL XPENUP(fxorig(XP(npl),YP(npl))
     :               ,fyorig(XP(npl),YP(npl)))
          CALL XPENDN(fxorig(XP(NPL+1),YP(NPL+1))
     :               ,fyorig(XP(NPL+1),YP(NPL+1)))
        enddo
      endif

      XS=XS+DX
      IF(XS.lt.Xwkmax) GOTO 100

      RETURN
      END

      SUBROUTINE getlnsgmnt(xwk,ywk,md,i,j,xs,
     :           ys1,ie1,je1,ie2,je2,ys2,ie3,je3,ie4,je4)
c
c This version does hatching of arbitary orientation
c between two contour values
c
      implicit none
c
c Input through argument list
c
      integer md
      REAL Xwk(MD,*),Ywk(MD,*)
      integer i,j
      integer iedge1(4),iedge2(4),jedge1(4),jedge2(4)
      real xs,ys1,ys2
      integer ie1,je1,ie2,je2,ie3,je3,ie4,je4
      integer nch
      common /xoutch/ nch
c
c Miselaneous local variables
c
      integer i1,j1,i2,j2,icount
      real tem,ys(4)

      icount = 0

      i1=i
      j1=j
      i2=i
      j2=j+1
      IF( (xwk(i1,j1)-xs)*(xwk(i2,j2)-xs).le.0.0 .and.
     :      (xwk(i1,j1).ne.xwk(i2,j2)) ) then
        tem=1.0/(xwk(i2,j2)-xwk(i1,j1))*(xs-xwk(i1,j1))
        icount = icount+1
        ys(icount)=ywk(i1,j1)+(ywk(i2,j2)-ywk(i1,j1))*tem
        iedge1(icount)=i1
        jedge1(icount)=j1
        iedge2(icount)=i2
        jedge2(icount)=j2
      endif

      i1=i
      j1=j+1
      i2=i+1
      j2=j+1
      IF( (xwk(i1,j1)-xs)*(xwk(i2,j2)-xs).le.0.0 .and.
     :      (xwk(i1,j1).ne.xwk(i2,j2))  ) then
        tem=1.0/(xwk(i2,j2)-xwk(i1,j1))*(xs-xwk(i1,j1))
        icount = icount+1
        ys(icount)=ywk(i1,j1)+(ywk(i2,j2)-ywk(i1,j1))*tem
        iedge1(icount)=i1
        jedge1(icount)=j1
        iedge2(icount)=i2
        jedge2(icount)=j2
      endif

      i1=i+1
      j1=j+1
      i2=i+1
      j2=j
      IF( (xwk(i1,j1)-xs)*(xwk(i2,j2)-xs).le.0.0 .and.
     :      (xwk(i1,j1).ne.xwk(i2,j2))  ) then
        tem=1.0/(xwk(i2,j2)-xwk(i1,j1))*(xs-xwk(i1,j1))
        icount = icount+1
        ys(icount)=ywk(i1,j1)+(ywk(i2,j2)-ywk(i1,j1))*tem
        iedge1(icount)=i1
        jedge1(icount)=j1
        iedge2(icount)=i2
        jedge2(icount)=j2
      endif

      i1=i
      j1=j
      i2=i+1
      j2=j
      IF( (xwk(i1,j1)-xs)*(xwk(i2,j2)-xs).le.0.0 .and.
     :      (xwk(i1,j1).ne.xwk(i2,j2))  ) then
        tem=1.0/(xwk(i2,j2)-xwk(i1,j1))*(xs-xwk(i1,j1))
        icount = icount+1
        ys(icount)=ywk(i1,j1)+(ywk(i2,j2)-ywk(i1,j1))*tem
        iedge1(icount)=i1
        jedge1(icount)=j1
        iedge2(icount)=i2
        jedge2(icount)=j2
      endif

      IF( icount.eq.0 .or. icount.gt.2 ) then
        write(nch,'(3(1x,a/))')
     :  'No or more than two intercepting side(s) found.',
     :  'Something is wrong.',
     :  'Program stopped in subroutine gtlnsgmnt.'
        Stop
      ENDIF

      if( ys(2).ge.ys(1)) then
        ie1=iedge1(1)
        je1=jedge1(1)
        ie2=iedge2(1)
        je2=jedge2(1)
        ie3=iedge1(2)
        je3=jedge1(2)
        ie4=iedge2(2)
        je4=jedge2(2)
        ys1=ys(1)
        ys2=ys(2)
      else
        ie1=iedge1(2)
        je1=jedge1(2)
        ie2=iedge2(2)
        je2=jedge2(2)
        ie3=iedge1(1)
        je3=jedge1(1)
        ie4=iedge2(1)
        je4=jedge2(1)
        ys1=ys(2)
        ys2=ys(1)
      endif

      RETURN
      END

      SUBROUTINE xsortxyp(xp,yp,n)
c sort xp and yp in ascending order in yp
      real xp(n),yp(n)

      do j=2,n
        a=yp(j)
        b=xp(j)
        do i=j-1,1,-1
          if(yp(i).le.a) goto 10
          yp(i+1)=yp(i)
          xp(i+1)=xp(i)
        enddo
        i=0
10      continue
        yp(i+1)=a
        xp(i+1)=b
      enddo

      return
      end

      SUBROUTINE XVECTU(U,V,MD,M,ISTEP,N,JSTEP,XLENG,UUNIT)
C
C Assess ranges if U,V values, and set length XLENG at which the unit
C vector UUNIT is plotted in x-direction.
C The length of vector in y-direction is scaled according to mapping.
C (Notice the non-isotropicity.)
C XLENG was set as XSCALE/(M-1)*ISTEP where XSCALE is the horizontal
C scale of mapped area.
C UUNIT was set that the longest vector falls between 0.75*XLENG
C and 1.5*XLENG in length.
c
c Fixed a problem with the first guess of umax,umin,vmax,vmin
c when the first value is missing.
c Jan. 24, 1995.
c
      REAL U(MD,N),V(MD,N)
      COMMON /XART36/ KARTYP,KVMODE,VSC
      COMMON /ZCHOLE/ NHOLE,SPECIA,nvtrbadv
      character*(*) sumax,sumin,svmax,svmin
      CHARACTER CH1*80
      SAVE UMAX, UMIN, VMAX, VMIN
      common /xoutch/ nch
      integer mxset

      mxset = 0

      DO 5 J=1,N,JSTEP
      DO 5 I=1,M,ISTEP
        IF(nvtrbadv.eq.1.and.(U(I,J).eq.SPECIA.or.V(I,J).eq.SPECIA))
     :  goto 5
        IF(mxset.eq.0) THEN
          UMAX=U(i,j)
          VMAX=V(i,j)
          UMIN=UMAX
          VMIN=VMAX
          mxset = 1
        ELSE
          UMAX=MAX(UMAX,U(I,J))
          UMIN=MIN(UMIN,U(I,J))
          VMAX=MAX(VMAX,V(I,J))
          VMIN=MIN(VMIN,V(I,J))
        ENDIF
5     CONTINUE

      WRITE(nch,'(4(a,F7.2),a)') ' Umin=',UMIN,' Umax=',UMAX,
     : ' Vmin=',VMIN,' Vmax=',VMAX,', -- ZXPLOT.'
      IF (UMAX == UMIN .AND. VMAX == VMIN ) GO TO 500

      UNIT=UUNIT
      IF( KVMODE.EQ.2) GOTO 105

 25   IF( MAX( ABS(UMAX), ABS(UMIN),abs(vmax),abs(vmin))
     :    .LT. UNIT*0.75 ) THEN
      UNIT=UNIT/2
      WRITE(nch,100) UNIT
 100  FORMAT(' Max vector < 0.75*UNIT , UNIT is halved. UNIT='
     :   ,3F9.4)
      GO TO 25
      ENDIF

 30   IF( MAX( ABS(UMAX),ABS(UMIN),abs(vmax),abs(vmin))
     :    .GT.UNIT*1.5  ) THEN
      UNIT=UNIT*2
      WRITE(nch,200) UNIT
 200  FORMAT(' Max vector > 1.5 *UNIT ,UNIT is doubled. UNIT='
     :   ,3F9.4)
      GO TO 30
      ENDIF
      CONTINUE
 105  UUNIT=UNIT

 500  CALL XQMAP(XL,XR,YB,YT)
      XLENG=(XR-XL)/(M-1)* ISTEP*VSC

      RETURN

      ENTRY XVSCAL( VSC0)
      VSC=VSC0
      RETURN

      ENTRY XVMODE(KVM)
      KVMODE=KVM
      RETURN

      ENTRY XVLMT(  CSIZE)
C Call XVLIMT with UMAX,UMIN,VMAX,VMIN saved in XVECTU.
C CSIZE set the size of characters. Default is about 0.012.
C Default value is assumed if CSIZE is 0.0.
      CALL XVLIMT(UMAX,UMIN,VMAX,VMIN, CSIZE)
      RETURN

      ENTRY XVLIMIT(x,y,sUMAX,sUMIN,sVMAX,sVMIN)

      WRITE(CH1,'(a,a,G9.3E2,3(a,a,a,G9.3E2))')
     :      sumax,'=',UMAX,',',sumin,'=',UMIN,',',
     :      svmax,'=',VMAX,',',svmin,'=',VMIN
      lch = 80
      CALL xstrmin (CH1,LCH)
      CALL XCHARC(x,y,ch1(1:lch))
      RETURN
      END

      SUBROUTINE XVECTU_MP( xl, xr, yb, yt,Umax, Umin, Vmax, Vmin,
     :                      m, istep, XLENG, UUNIT)
C
C Based on input ranges if U,V values, and set length XLENG at which the unit
C vector UUNIT is plotted in x-direction.
C The length of vector in y-direction is scaled according to mapping.
C (Notice the non-isotropicity.)
C XLENG was set as XSCALE/(M-1)*ISTEP where XSCALE is the horizontal
C scale of mapped area.
C UUNIT was set that the longest vector falls between 0.75*XLENG
C and 1.5*XLENG in length.
c
c It is based on XVECTU, but written specificially for MPI mode when
c global UMAX, UMIN, VMAX, VMIN are known beforehand.
c     -- By Y. Wang 10/30/2009
c
      IMPLICIT NONE
      REAL :: xl, xr, yb, yt
      REAL :: Umax, Umin, Vmax, Vmin
      INTEGER :: m, istep
      REAL :: XLENG, UUNIT

      INTEGER :: kartyp, kvmode
      REAL    :: vsc
      COMMON /XART36/ KARTYP,KVMODE,VSC

      INTEGER :: nhole, nvtrbadv
      REAL    :: specia
      COMMON /ZCHOLE/ NHOLE,SPECIA,nvtrbadv

      INTEGER :: nch
      common /xoutch/ nch

      REAL :: tunit

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!      WRITE(nch,'(4(a,F7.2),a)') ' Umin=',UMIN,' Umax=',UMAX,
!     : ' Vmin=',VMIN,' Vmax=',VMAX,', -- ZXPLOT.'

      IF (UMAX == UMIN .AND. VMAX == VMIN ) GO TO 500

      tUNIT=UUNIT
      IF( KVMODE.EQ.2) GOTO 105

 25   IF( MAX( ABS(UMAX), ABS(UMIN),abs(vmax),abs(vmin))
     :    .LT. tUNIT*0.75 ) THEN

        tUNIT=tUNIT/2
        WRITE(nch,100) tUNIT
 100    FORMAT(' Max vector < 0.75*UNIT , UNIT is halved. UNIT='
     :          ,3F9.4)
        GO TO 25
      ENDIF

 30   IF( MAX( ABS(UMAX),ABS(UMIN),abs(vmax),abs(vmin))
     :    .GT. tUNIT*1.5  ) THEN
        tUNIT=tUNIT*2
        WRITE(nch,200) tUNIT
 200    FORMAT(' Max vector > 1.5 *UNIT ,UNIT is doubled. UNIT='
     :          ,3F9.4)
        GO TO 30
      ENDIF
      CONTINUE
 105  UUNIT=tUNIT

 500  CALL XQMAP(XL,XR,YB,YT)
      XLENG=(XR-XL)/(M-1)* ISTEP*VSC

      RETURN
      END SUBROUTINE XVECTU_MP


      SUBROUTINE XVECTR(U,V,X,Y,MD,M,ISTEP,N,JSTEP,XLENG,UUNIT)
C Plot vector feilds ( U,V)   . Unit X-component UUNIT is plotted with
C length XLENG in mapped maths coordinate. Length of vector in
C Y-direction is scaled according to the mapping.
      REAL U(MD,N),V(MD,N),X(MD,N),Y(MD,N)
      COMMON /ZCHOLE/ NHOLE,SPECIA,nvtrbadv
      common /xwndw1/ xw1,xw2,yw1,yw2, iwndon

      DO 15 J=1,N,JSTEP
      DO 15 I=1,M,ISTEP
        IF(iwndon.eq.0 .or. ((x(i,j)-xw1)*(x(i,j)-xw2).le.0.0 .and.
     :                       (y(i,j)-yw1)*(y(i,j)-yw2).le.0.0) ) THEN
        IF(nvtrbadv.eq.0) THEN
          CALL XARROW( U(I,J),V(I,J),X(I,J),Y(I,J),XLENG,UUNIT)
        ELSE
          IF(u(i,j).ne.SPECIA.and.v(i,j).ne.SPECIA )
     :    CALL XARROW( U(I,J),V(I,J),X(I,J),Y(I,J),XLENG,UUNIT)
        ENDIF
        ENDIF
 15   CONTINUE
      RETURN
      END

      SUBROUTINE XBARBS(U,V,X,Y,MD,M,ISTEP,N,JSTEP,wunits,xleng,barbopt)

C     Plot wind barbs for vector feilds ( U,V).

      REAL U(MD,N),V(MD,N),X(MD,N),Y(MD,N)
      COMMON /ZCHOLE/ NHOLE,SPECIA,nvtrbadv
      common /xwndw1/ xw1,xw2,yw1,yw2, iwndon
      integer wunits,barbopt


      DO 15 J=1,N,JSTEP
      DO 15 I=1,M,ISTEP
        IF(iwndon.eq.0 .or. ((x(i,j)-xw1)*(x(i,j)-xw2).le.0.0 .and.
     :                       (y(i,j)-yw1)*(y(i,j)-yw2).le.0.0) ) THEN
        IF(nvtrbadv.eq.0) THEN
          CALL XBARB(U(I,J),V(I,J),X(I,J),Y(I,J),wunits,XLENG,
     :               barbopt)
        ELSE
          IF(u(i,j).ne.SPECIA.and.v(i,j).ne.SPECIA )
     :    CALL XBARB(U(I,J),V(I,J),X(I,J),Y(I,J),wunits,XLENG,
     :               barbopt)
        ENDIF
        ENDIF
 15   CONTINUE
      RETURN
      END

      SUBROUTINE XBARB(U,V,X0,Y0,wunits,XLENG,barbopt)

c     Plot a wind barb at (x0, y0) for wind vector (u,v).
C     Its length is specified in terms of the length in
c     the x-coordinate direction.

      IMPLICIT NONE

      real u,v,x0,y0
      integer wunits ! Wind units. =1: m/s, =2: knots or miles/per hour
      real xleng
      integer barbopt ! Option for plotting the direction of wind barb
                      ! =1, wind barb direction conforms to the streamlines if plotted,
                      !     i.e., it depends on the grid aspect ration.
                      ! =2, wind bard direction represents the absolute direction,
                      !     regardless the grid aspect ratio.

      real pi, angle1, sina1, cosa1
c     PARAMETER(PI=3.14159,ANGLE1=(120./180+1)*PI)
      PARAMETER(PI=3.14159,ANGLE1=(105./180+1)*PI)
      real speed, sinta,costa,arrow
      integer nhalf, nfull, nfifty,i
      real dx,dy,px0,py0,px1,py1,dpx,dpy,px1a,py1a
      real XPNTSD,XPLENG,DPXY,dx1,dy1

      real xfactr,yfactr
      COMMON /XFTR06/ XFACTR,YFACTR

c     sina1=-0.866027355 ! sin( angle1 )
c     cosa1= 0.499996603 ! cos( ANGLE1 )

c     sina1=sin( angle1 )
c     cosa1=cos( angle1 )

      sina1= -0.965925826
      cosa1=  0.258819045

      speed = sqrt( u*u + v*v)
      IF( speed .lt. 1.0e-10 ) RETURN

      If (wunits.EQ.1) Then
        nhalf   = Nint (speed/2.5)        !for m/s
      ElseIf (wunits.EQ.2) Then
        nhalf   = Nint (speed/5.0)        !for mph/knots
      End If
      nfifty = nhalf / 10
      nhalf  = nhalf - nfifty * 10
      nfull  = nhalf / 2
      nhalf  = nhalf - nfull * 2

      XPLENG =XPNTSD( 0.0,0.0,XLENG,0.0)

      PX0=X0
      PY0=Y0
      CALL XTRANS(PX0,PY0)

      IF( barbopt.eq.1) then
        DX=u/speed*xleng
        DY=v/speed*xleng
        PX1=X0-DX
        PY1=Y0-DY
        CALL XTRANS(PX1,PY1)
        DPX=PX1-PX0
        DPY=PY1-PY0
      else
        DPX=-u*XPLENG/speed
        DPY=-v*XPLENG/speed
        px1=px0+dpx
        py1=py0+dpy
      endif

      DPXY=SQRT( DPX*DPX+DPY*DPY)
      IF(DPXY.le.1.0E-30) RETURN

      SINTA=DPY/DPXY
      COSTA=DPX/DPXY

      ARROW=0.50* XPLENG
      DX1=ARROW*(COSTA*COSA1-SINTA*SINA1)
      DY1=ARROW*(SINTA*COSA1+COSTA*SINA1)
      CALL XTPNUP(PX0,PY0)
      CALL XTPNDN(PX1,PY1)

      px1a = px1
      py1a = py1

      DO i=1,nfifty ! Plot flags
        CALL XTPNUP(PX1a    , PY1a    )
        CALL XTPNDN(PX1a+DX1, PY1a+DY1)
        px1a = px1a+0.30*(px0-px1)
        py1a = py1a+0.30*(py0-py1)
        CALL XTPNDN(PX1a,PY1a)
      ENDDO

      DO i=1,nfull ! Plot full-length barbs
        CALL XTPNUP(PX1a    , PY1a    )
        CALL XTPNDN(PX1a+DX1, PY1a+DY1)
        px1a = px1a+0.15*(px0-px1)
        py1a = py1a+0.15*(py0-py1)
      ENDDO

      if( nhalf.ge.1 .and. (nfifty.eq.0.and.nfull.eq.0)) then
        px1a = px1a+0.20*(px0-px1)
        py1a = py1a+0.20*(py0-py1)
      endif

      DO i=1,nhalf ! Plot half-length barbs
        CALL XTPNUP(PX1a    , PY1a    )
        CALL XTPNDN(PX1a+DX1*0.5, PY1a+DY1*0.5)
        px1a = px1a+0.15*(px0-px1)
        py1a = py1a+0.15*(py0-py1)
      ENDDO

      RETURN
      END SUBROUTINE xbarb

      SUBROUTINE XCONTS(Z,X,Y,MD,ND, ZINC)
      REAL Z(MD,ND),X(MD,ND),Y(MD,ND),CL(150)
      INTEGER IWRK(10000)
      CL(1)=0.0
      CL(2)=ZINC
      MODE=1
      M =MD
      N =ND
      IST=1
      JST=1
      CALL XCONTA(Z(IST,JST),X(IST,JST),Y(IST,JST)
     :           ,IWRK,MD,M,N,CL,NCL,MODE)

      IF( NCL.EQ.0) THEN
           ZMAX= CL(1)
           ZINC1= ZMAX
      ELSE
           ZMAX=CL(NCL)
           ZINC1=CL(2)-CL(1)
       ENDIF
      CALL XCLIMT(ZMAX, CL(1),ZINC1, 0.0)

      RETURN
      END

      SUBROUTINE XVECTS(U,V,X,Y,MD,ISTEP,ND,JSTEP,XLENG1,UUNIT1)
      REAL U(MD,ND),V(MD,ND),X(MD,ND),Y(MD,ND)

      M =MD
      N =ND
      IST=1
      JST=1
      XLENG=XLENG1
      UUNIT=UUNIT1
      CALL XVECTU(U(IST,JST),V(IST,JST),MD,M,ISTEP,N,JSTEP,XLENG,UUNIT)
      CALL XVECTR(U(IST,JST),V(IST,JST),X(IST,JST),Y(IST,JST),
     :            MD,M,ISTEP,N,JSTEP,XLENG,UUNIT)
      CALL XQMAP(XL,XR,YB,YT)
      X0=XL+(XR-XL)*0.75
      Y0=YT+(YT-YB)*0.03
      KEY=0
      AM=1.0
      IF( (M-1)/ISTEP.GT.30) AM=2.0
      CALL XVECTK(X0,Y0,XLENG*AM,UUNIT*AM, KEY)
      CALL XVLMT(0.0)
      RETURN
      END

      SUBROUTINE XCAPTN(TITLES,NUM,CH, LC )
C PLOT CAPTIONS ALONG THE BORDER
      common /xoutch/ nch
      CHARACTER TITLES(NUM)*50, CH*100
      CALL XQMAP(XL,XR,YB,YT)
      CALL XQRANG( XRG, YRG )
      CALL XQCHOR( ANG0 )
      CALL XQCHMG( SIZ0 )
      SIZ=0.05*MIN( XRG, YRG)
      HX=SIZ*(XR-XL)/XRG
      HY=SIZ*(YT-YB)/YRG
      CALL XCHMAG( SIZ )
      IF( NUM.GE.1) THEN
        CALL XQOBAG( XANG, YANG )
        CALL XCHORI( 90.0 +YANG-XANG )
        CALL XCHARC( XL- 2.5*HX,  0.5*(YT+YB), TITLES(1)(1:10) )
        CALL XCHORI(0.0)
        CALL XCHARC( 0.5*(XR+XL),  YB-2*HY , TITLES(1)(11:20))
      ENDIF
      IF( NUM.GE.2) THEN
        CALL XCHARL(XL,YB-4*HY,TITLES(2)(1:20)//'  '//ch(1:lc))
      ENDIF
        YCP=YB-4.5*HY
      DO 10 K=3,NUM
        YCP=YCP-HY
        CALL XCHARL(XL,YCP,TITLES(K))
 10   CONTINUE
      CALL XCHORI( ANG0 )
      CALL XCHMAG( SIZ0 )
      WRITE(NCH,*) TITLES(2)(1:20) ,' is to be plotted..'
      RETURN
      END

      SUBROUTINE XCLIMT(FMAX,FMIN,FINC ,CTRSIZ)
      CHARACTER CH*150
      CALL XQMAP(XL,XR,YB,YT)
      CALL XQRANG( XRG, YRG )
      CALL XQCHMG( SIZ0 )
      IF(ctrsiz.eq.0) THEN
        SIZ=0.03*MIN( XRG, YRG)
        CALL XCHMAG(SIZ)
      ELSE
        CALL XCHSIZ(CTRSIZ)
      ENDIF
      WRITE(CH,'(''Min='',G9.3E2,'' Max='',G9.3E2,
     :    '' Contour interval='',G9.3E2)')FMIN,FMAX,FINC
      lch = 54
      CALL xstrmin ( CH, LCH)
      CALL XCHARC( 0.5*(XL+XR), YT+0.02*(YT-YB),CH(1:lch) )
      CALL XCHMAG( SIZ0)
      RETURN
      END

      SUBROUTINE XVLIMT(UMAX,UMIN,VMAX,VMIN , ctrsiz)
      CHARACTER CH1*80, ch2*80
      CALL XQMAP(XL,XR,YB,YT)
      CALL XQCHMG( SIZ0 )
      CALL XQRANG( XRG, YRG )
      CALL XQCHMG( SIZ0 )
      IF( ctrsiz .eq. 0.0 ) THEN
        SIZ=0.03*MIN( XRG, YRG)
        CALL XCHMAG(SIZ)
      ELSE
        CALL XCHSIZ(CTRSIZ)
      ENDIF
      WRITE(CH1,'('' Umax='',G9.3E2,'' Umin='',G9.3E2)')UMAX,UMIN
      WRITE(CH2,'('' Wmax='',G9.3E2,'' Wmin='',G9.3E2)')VMAX,VMIN
      lch = 30
      CALL xstrmin ( CH1, LCH)
      CALL XCHARL(XL,YT+0.02*(YT-YB),CH1(1:lch))
      lch = 30
      CALL xstrmin ( CH2, LCH)
      CALL XCHARL(XL,YT+0.06*(YT-YB),CH2(1:lch) )
      CALL XCHMAG( SIZ0)
      RETURN
      END

      SUBROUTINE XVECTK(x0,y0, xleng, uunit, key)
c
c#######################################################################
c     PURPOSE:
c     Plot unit vectors  starting at (X0,Y0)
C     KEY=-1, 0, 1, 2, for none,in both X and Y-direction,X only, Y only
c     INPUT: x0 y0 xleng uunit key c
c#######################################################################
c
      implicit none
c
      real x0,y0
      real xleng,uunit
      integer key
c
c#######################################################################
c
c     Misc. local Variables
c
c#######################################################################
c
      real pi,angle1,angle2
      parameter(pi=3.14159,angle1=(10./180+1)*pi,angle2=(-10./180+1)*pi)
      real sina1,cosa1,sina2,cosa2
      parameter(sina1=-.17365,cosa1=-.98481,sina2= -sina1,cosa2=cosa1)
      real xrg,xl,xr,xscale
      real yrg,yb,yt,yscale,yf
      real dx,dx1,dx2
      real dy,dy1,dy2
      real px0,px1,px2
      real py0,py1,py2
      real dph,dpv
      real costa,sinta
      real vunit
      real arrow
      real asym
      real xang
      real yang
      integer lch
      character ch*20

      INTEGER kartyp, KVMODE
      REAL vsc

      REAL :: vfactor
      REAL :: vlenfactor        ! vertical vector factor
      DATA vlenfactor/1.0/

      COMMON /VFACTOR/ vlenfactor    ! vertical vector factor

      COMMON /XART36/ KARTYP,KVMODE,VSC
c
      CALL xqrang(xrg,yrg)
      CALL xqmap(xl,xr,yb,yt)
      xscale=xr-xl
      yscale=yt-yb
      yf=0.4+( min(xrg,yrg)-0.4)*0.5

      vunit=uunit
      dx=xleng

      IF( KARTYP==3 .OR. KARTYP==4 ) then
        dy=xleng* (yt-yb)/(xr-xl) * xrg/yrg
      ELSE
        dy=xleng*vlenfactor
      ENDIF

      px0=x0
      py0=y0
      px1=x0+dx
      py1=y0
      px2=x0
      py2=y0+dy
      CALL xtrans(px0,py0)
      CALL xtrans(px1,py1)
      CALL xtrans(px2,py2)
      dph=sqrt( (px1-px0)**2+(py1-py0)**2  )
      dpv=sqrt( (px2-px0)**2+(py2-py0)**2  )
 5    IF( dpv.gt.(1.5*dph) ) THEN
        dpv=dpv*0.5
        dy=dy*0.5
        vunit=vunit*0.5
        GO TO 5
      ENDIF
 6    IF( dpv.lt.(0.75*dph) ) THEN
        dpv=dpv*2
        dy=dy*2
        vunit=vunit*2
        GO TO 6
      ENDIF
      px2=x0
      py2=y0+dy
      CALL xtrans(px2,py2)

      IF( key.eq.0.or.key.eq.1) THEN
        arrow=0.30*dph
        costa=(px1-px0)/dph
        sinta=(py1-py0)/dph
        dx1=arrow*(costa*cosa1-sinta*sina1)
        dy1=arrow*(sinta*cosa1+costa*sina1)
        dx2=arrow*(costa*cosa2-sinta*sina2)
        dy2=arrow*(sinta*cosa2+costa*sina2)
        CALL xtpnup(px0,py0)
        CALL xtpndn(px1    , py1)
        CALL xtpndn(px1+dx1, py1+dy1)
        CALL xtpnup(px1    , py1    )
        CALL xtpndn(px1+dx2, py1+dy2)
        write(ch,'(f6.1)') uunit
        lch=6
        CALL xchlj (ch,lch)
        call XSTRLNTH(ch,lch)
        CALL xcharl(x0+dx   +0.01*xscale,y0,ch(1:lch) )
      ENDIF
c
      IF( key.eq.0.or.key.eq.2 ) THEN
        arrow=0.30*dph
        costa=(px2-px0)/dpv
        sinta=(py2-py0)/dpv
        dx1=arrow*(costa*cosa1-sinta*sina1)
        dy1=arrow*(sinta*cosa1+costa*sina1)
        dx2=arrow*(costa*cosa2-sinta*sina2)
        dy2=arrow*(sinta*cosa2+costa*sina2)
        CALL xtpnup(px0,py0)
        CALL xtpndn(px2    , py2)
        CALL xtpndn(px2+dx1, py2+dy1)
        CALL xtpnup(px2    , py2    )
        CALL xtpndn(px2+dx2, py2+dy2)
        CALL xqobag( xang, yang )
        CALL xqchor( asym )
        CALL xchori(90.0+ yang- xang)
        write(ch,'(f6.1)') vunit
        lch=6
        CALL xchlj (ch,lch)
        call XSTRLNTH(ch,lch)
        CALL xcharl(x0-.01*xscale ,y0 ,ch(1:lch) )
        CALL xchori( asym )
      ENDIF
      RETURN

      ENTRY set_vertical_factor(vfactor)
      IF (vfactor <= 0.0 .OR. vfactor >= 1.0) THEN
        vlenfactor = 1.0
      ELSE
        vlenfactor = vfactor
      END IF
      RETURN

      END

      SUBROUTINE XGRAPH( X,Y,N )
C Plot single valued curve y=y(x), where for each x,there is a unique y.
C METHOD-- parameter controling curve plotting pattern
C          =0 ,points are joined up by straight lines
C          =1 or 2, two points are joined up using seamed quadratics.
C when =1,slope at data points are calculated using monotonic method
C when =2, using Bessel methed.
      REAL X(*), Y(*)
      IF( N.LE.1) RETURN
      IF( N.EQ.2) THEN
        CALL XLPNUP( X(1),Y(1))
        CALL XLPNDN( X(2),Y(2))
        RETURN
      ENDIF
      CALL XQCVMD( MTD )
      IF( MTD.EQ.0) THEN
        CALL XLPNUP(X(1),Y(1))
      ELSE
        CALL XGRPUP(X(1),Y(1))
      ENDIF
      NEND=0
      DO 10 I=2,N
        IF( MTD.EQ.0) THEN
          CALL XLPNDN(X(I),Y(I))
        ELSE
          IF( I.EQ.N) NEND=1
          CALL XGRPDN(X(I),Y(I),NEND)
        ENDIF
 10   CONTINUE
      CALL XLPNUP( X(N),Y(N) )
      RETURN
      END

      SUBROUTINE XGRPUP(X,Y)
C Position pen at the starting point of a sigle valued curve y=y(x).
      COMMON /SVCURV/ X1,Y1,X2,Y2,X3,Y3,NPTS, SL1,SL2
      COMMON /XCDV23/ NSUBDV
      COMMON /XCMD24/ MTD
      common /xoutch/ nch
      X2=X
      Y2=Y
      NPTS=1
      CALL XLPNUP(X2,Y2)
      RETURN

      ENTRY XGRPDN(X,Y,NEND)
C Join (x,y(x)) with a smooth seamed quadratic curve, where y=y(x) is a
C single valued function.
      KEND=NEND
      IF( MTD.EQ.0.OR .NSUBDV.LE.1) THEN
        CALL XLPNDN(X,Y)
        RETURN
      ENDIF
      NEND1=NEND
      IF( NPTS.LT.3)GOTO 8
      IF( X.EQ.X2)THEN
        IF(Y.NE.Y2)GOTO 999
        IF( NEND1.EQ.1) GOTO 22
        GOTO 5
      ENDIF
      NPTS=NPTS+1
      X1=X2
      Y1=Y2
      X2=X3
      Y2=Y3
      X3=X
      Y3=Y
      GOTO 13
 8    IF(NPTS.EQ.2) GOTO 10
      IF( X.EQ.X2) THEN
        IF(Y.NE.Y2)GOTO 999
        GOTO 5
      ENDIF
         NPTS=NPTS+1
         X3=X
         Y3=Y
         GOTO 5
 10   IF( X.EQ.X3) THEN
        IF(Y.NE.Y3)GOTO 999
        GOTO 5
      ENDIF
        NPTS=NPTS+1
        X4=X
        Y4=Y
      CALL XQUADR(X2,Y2,X3,Y3,X4,Y4,A,B,C)
      X1=X2-(X4-X3)
      Y1=C+(B+A*X1)*X1
 13   IF( (X3-X2)*(X2-X1).LT.0.0) GOTO 997
      KOUNT=1
      IF( NPTS.EQ.3) KOUNT=2
 21   DO 20 J=1,KOUNT
        D21=(Y2-Y1)/(X2-X1)
        D32=(Y3-Y2)/(X3-X2)
        IF(MTD.EQ.2) THEN
          DD31=(D32-D21)/(X3-X1)
          SL2=( D21+D32-DD31*(X1-2*X2+X3))*0.5
          GOTO 16
        ELSEIF( MTD.EQ.1) THEN
            IF( D21*D32.GT. 0.0) THEN
              AA=(1.0+(X3-X2)/(X3-X1))/3.0
              SLINV=AA/D21+(1-AA)/D32
              SL2=1.0/SLINV
            ELSE
              SL2=0.0
            ENDIF
            GOTO 16
        ENDIF
 16   IF( KOUNT.EQ.2.AND.J.EQ.1) THEN
        SL1=SL2
        X1=X2
        Y1=Y2
        X2=X3
        Y2=Y3
        X3=X4
        Y3=Y4
      ENDIF
 20   CONTINUE
      CALL XSEAMQ( X1,Y1,SL1,X2,Y2,SL2,A1,A2,B,C)
      XSUB=X1
      XINC=(X2-X1)/NSUBDV
      XC=(X2+X1)*0.5
      KEND=0
      DO 15 ISUB=2,NSUBDV
        XSUB=XSUB+XINC
        A12=A1
        IF( ISUB.GT.NSUBDV/2 ) A12=A2
        YSUB=C +(B +A12*(XSUB-XC))*(XSUB-XC)
 15   CALL XLPNDN( XSUB,YSUB)
      IF( NEND.EQ.1.AND.NEND1.EQ.0) KEND=1
      CALL XLPNDN( X2,Y2)
        SL1=SL2
 22   IF(NEND1.EQ.0) THEN
        RETURN
      ELSE
        CALL XQUADR(X1,Y1,X2,Y2,X3,Y3,A,B,C)
        X4=X3-(X1-X2)
        Y4=C+(B+A*X4)*X4
        X1=X2
        Y1=Y2
        X2=X3
        Y2=Y3
        X3=X4
        Y3=Y4
        NEND1=0
        KOUNT=1
        GOTO 21
      ENDIF
 5    IF( NEND1.NE.1) RETURN
        IF( NPTS.EQ.2) THEN
          KEND=1
          CALL XLPNDN( X2,Y2)
        ELSE
          RETURN
        ENDIF
      RETURN

 999  WRITE(NCH,*)
     :  ' Input data are controdicting! Curve plotting aborted.'
      RETURN

 997  WRITE(NCH,*)' Input data not in correct order! Plotting aborted.'
      RETURN
      END

      SUBROUTINE XCURVE( X,Y,N, KLOSE)
C Plot multiple-valued curve X(t), Y(t).
C METHOD-- parameter controling curve plotting pattern.
C          =0 ,points are joined up by straight lines
C  =1 or 2, two points are joined up using parametric seamed quadratics.
C when =1,slope at data points are calculated using monotonic method
C when =2, using Bessel methed.
      REAL X(*), Y(*)
      IF( N.LE.1) RETURN
      IF( N.EQ.2) THEN
        CALL XLPNUP( X(1),Y(1))
        CALL XLPNDN( X(2),Y(2))
        RETURN
      ENDIF
      CALL XQCVMD(MTD)
      IF( MTD.EQ.0) THEN
        CALL XLPNUP(X(1),Y(1))
      ELSE
        CALL XCURUP(X(1),Y(1))
      ENDIF
      NEND=0
      DO 10 I=2,N
        IF( MTD.EQ.0) THEN
          CALL XLPNDN(X(I),Y(I))
        ELSE
          IF( I.EQ.N) NEND=1
          CALL XCURDN(X(I),Y(I), KLOSE, NEND)
        ENDIF
 10   CONTINUE
      IF( MTD.EQ.0.AND.KLOSE.EQ.1.AND.(X(1).NE.X(N).OR.Y(1).NE.Y(N)))
     :  CALL XLPNDN(X(1), Y(1))
      CALL XLPNUP( X(1),Y(1) )
      RETURN
      END

      SUBROUTINE XCURUP(X,Y)
      COMMON /MVCURV/ X1,Y1,X2,Y2,X3,Y3,T2,NPTS ,SL1X,SL1Y,SL2X,SL2Y
      SAVE X01,Y01,X02,Y02,X03,Y03,IEND
      COMMON /XCDV23/ NSUBDV
      COMMON /XCMD24/ MTD
      X2=X
      Y2=Y
      NPTS=1
      IEND=0
      CALL XLPNUP(X2,Y2)
      RETURN

      ENTRY XCURDN(X,Y,KCLOSE,NEND)
      IF( MTD.EQ.0.OR .NSUBDV.LE.1) THEN
        CALL XLPNDN(X,Y )
        RETURN
      ENDIF
      NEND1=NEND
      KLOSE=KCLOSE
      IF( NPTS.LT.3)GOTO 8
      IF( X.EQ.X2.AND.Y.EQ.Y2) THEN
        IF( NEND1.EQ.1) GOTO 22
        GOTO 5
      ENDIF
      NPTS=NPTS+1
      T2=T2+1.
      X1=X2
      Y1=Y2
      X2=X3
      Y2=Y3
      X3=X
      Y3=Y
      GOTO 13
 8    IF(NPTS.EQ.2) GOTO 10
      IF( X.EQ.X2.AND.Y.EQ.Y2) GOTO 5
         NPTS=NPTS+1
         X3=X
         Y3=Y
         IF( KLOSE.EQ.1) CALL XLPNUP(X3,Y3)
         GOTO 5
 10   IF( X.EQ.X3.AND.Y.EQ.Y3) GOTO 5
        NPTS=NPTS+1
        X4=X
        Y4=Y
      T2=2.0
      IF( KLOSE.EQ.1) THEN
        X1=X2
        Y1=Y2
        X2=X3
        Y2=Y3
        X3=X4
        Y3=Y4
        X01=X1
        X02=X2
        X03=X3
        Y01=Y1
        Y02=Y2
        Y03=Y3
        KOUNT=1
        GOTO 21
      ENDIF
      T1=0.0
      CALL XQUADR(T2-1,X2,T2,X3,T2+1,X4,A,B,C)
      X1=C
      CALL XQUADR(T2-1,Y2,T2,Y3,T2+1,Y4,A,B,C)
      Y1=C
 13   KOUNT=1
      IF( NPTS.EQ.3.AND.KLOSE.NE.1)  KOUNT=2
 21   DO 20 J=1,KOUNT
        IF(MTD.EQ.2) THEN
          SL2X=(X3-X1)*0.5
          SL2Y=(Y3-Y1)*0.5
          GOTO 16
        ELSEIF( MTD.EQ.1) THEN
          D21X= X2-X1
          D32X= X3-X2
          D21Y= Y2-Y1
          D32Y= Y3-Y2
          IF( D21X*D32X.GT. 0.0) THEN
            SLINV=0.5/D21X+0.5/D32X
            IF(SLINV.EQ.0.0) THEN
              SL2X=0.0
            ELSE
              SL2X=1.0/SLINV
            ENDIF
          ELSE
            SL2X=0.0
          ENDIF
          IF( D21Y*D32Y.GT. 0.0) THEN
            SLINV=0.5/D21Y+0.5/D32Y
            IF(SLINV.EQ.0.0) THEN
              SL2Y=0.0
            ELSE
              SL2Y=1.0/SLINV
            ENDIF
          ELSE
            SL2Y=0.0
          ENDIF
          GOTO 16
        ENDIF
 16   CONTINUE
      IF(KLOSE.EQ.1.AND.NPTS.LE.3) THEN
        SL1X=SL2X
        SL1Y=SL2Y
        RETURN
      ENDIF
      IF( KOUNT.EQ.2.AND.J.EQ.1) THEN
        SL1X=SL2X
        SL1Y=SL2Y
        X1=X2
        Y1=Y2
        X2=X3
        Y2=Y3
        X3=X4
        Y3=Y4
        T2=2.0
      ENDIF
 20   CONTINUE

      CALL XSEAMQ( T2-1.0,X1,SL1X,T2,X2,SL2X,A1X,A2X,BX,CX)
      T=T2-1
      TC=T2-0.5
      XX=CX+(BX+A1X*(T-TC))*(T-TC)
      CALL XSEAMQ( T2-1.0,Y1,SL1Y,T2,Y2,SL2Y,A1Y,A2Y,BY,CY)
      TINC=1.0/NSUBDV
      KEND=0
      DO 15 ISUB=2,NSUBDV
        T=(ISUB-1)*TINC-0.5
        A12X=A1X
        A12Y=A1Y
        IF( ISUB.GT.NSUBDV/2 ) THEN
          A12Y=A2Y
          A12X=A2X
        ENDIF
        XSUB=CX+(BX+A12X*T)*T
        YSUB=CY+(BY+A12Y*T)*T
 15   CALL XLPNDN( XSUB,YSUB)
      IF(KLOSE.EQ.1) THEN
         IF( IEND.EQ.3) KEND=1
      ELSE
         IF( NEND .EQ.1.AND.NEND1.EQ.0  ) KEND=1
      ENDIF
      CALL XLPNDN( X2,Y2)
        SL1X=SL2X
        SL1Y=SL2Y

 22   IF(NEND1.EQ.0) THEN
        RETURN
      ELSEIF(KLOSE.EQ.1) THEN
        IF(IEND.GE.3) RETURN
        X1=X2
        X2=X3
        Y1=Y2
        Y2=Y3
 27     IEND=IEND+1
        IF( IEND.EQ.1) THEN
          IF(X3.EQ.X01.AND.Y3.EQ.Y01) GOTO 27
          X3=X01
          Y3=Y01
        ELSEIF( IEND.EQ.2) THEN
          X3=X02
          Y3=Y02
        ELSEIF( IEND.EQ.3) THEN
          X3=X03
          Y3=Y03
        ENDIF
        NPTS=NPTS+1
        T2=T2+1
        KOUNT=1
          GOTO 21
        ELSE
        T4=NPTS+1.0
        CALL XQUADR(T2-1,X1,T2,X2,T2+1,X3,A,B,C)
        X4=C+(B+A*T4)*T4
        CALL XQUADR(T2-1,Y1,T2,Y2,T2+1,Y3,A,B,C)
        Y4=C+(B+A*T4)*T4
        X1=X2
        Y1=Y2
        X2=X3
        Y2=Y3
        X3=X4
        Y3=Y4
        T2=T2+1
        NEND1=0
        KOUNT=1
        GOTO 21
      ENDIF
 5    IF( NEND1.EQ.1.AND. NPTS.EQ.2) THEN
          KEND=1
          CALL XLPNDN( X3,Y3)
      ENDIF
      RETURN
      END

      SUBROUTINE XCURDV(NDIV )
      COMMON /XCDV23/ NSUBDV
      NSUBDV= NDIV
      RETURN
      END

      SUBROUTINE XCVMTD( METHOD )
C Set the method for curve plotting when using XGRPUP,XGRPDN,XGRAPH,
C XCURUP,XCURUP,XCURDN. By default METHOD=0
      COMMON /XCMD24/ MTD
      MTD= METHOD
      RETURN
      ENTRY XQCVMD( MTD1 )
      MTD1=MTD
      RETURN
      END

      SUBROUTINE XQUADR(X1,Y1,X2,Y2,X3,Y3,A,B,C)
C Return the coefficients A,B,C of a quadratic polynomial fitting
C points (X1,Y1),(X2,Y2),(X3,Y3)
      D21=(Y2-Y1)/(X2-X1)
      D32=(Y3-Y2)/(X3-X2)
      A=(D32-D21)/(X3-X1)
      B=(D21*(X3+X2)-D32*(X1+X2))/(X3-X1)
      C=(X3*Y1-X1*Y3+X1*X3*(D32-D21))/(X3-X1)
      RETURN
      END

      SUBROUTINE XSEAMQ(X1,Y1,SL1,X2,Y2,SL2,A1,A2,B,C)
C Fit points (X1,Y1) and(X2,Y2) with slopes SL1,SL2 at the corresponding
C seamed quadratics so that y=c+(b+a*(x-xc))*(x-xc) where xc=(x1+x2)/2
C and a=a1 for (x2-x1)*(x-xc)<0.0 ,a=a2 for (x2-x1)*(x-xc)>=0.0
      H=(X2-X1)*0.5
      HH=H*H
      TL=Y1+0.5*H*SL1
      TR=Y2-0.5*H*SL2
      C =0.5*(TL+TR)
      B =(TR-TL)/H
      A1=(Y1-1.5*TL+0.5*TR)/HH
      A2=(Y2-1.5*TR+0.5*TL)/HH
      RETURN
      END

C* COLOUR FILLING ROUTINES
      SUBROUTINE ZCONTB(ZG,Z,MD,MG,JG,C1,C2)
C Colour filling of triangular blocks
      DIMENSION ZG(MD,*),Z(MD,*),X(5),Y(5)
      COMPLEX   Z,B1,B2,ZA,ZB,ZC,ZD,Z1,Z2,Z3,Z4,Z5,D
                D(P1,P2,B1,B2   )=B1+(CV-P1)*(B2-B1)/(P2-P1)
      CV1=MIN(C1,C2)
      CV2=MAX(C1,C2)
      DO 101 I=1,MG-1
      DO 102 J=1,JG-1
      H5=0.25*(ZG(I,J)+ZG(I,J+1)+ZG(I+1,J+1)+ZG(I+1,J))
      Z5=0.25*( Z(I,J)+ Z(I,J+1)+ Z(I+1,J+1)+ Z(I+1,J))
      DO 103 M=1,4
      H1=ZG(I+MOD(M  ,4)/2,J+MOD(M-1,4)/2)
      Z1=Z (I+MOD(M  ,4)/2,J+MOD(M-1,4)/2)
      H2=ZG(I+MOD(M+1,4)/2,J+MOD(M  ,4)/2)
      Z2=Z (I+MOD(M+1,4)/2,J+MOD(M  ,4)/2)
      H3=H5
      Z3=Z5
      IF(H2.LT.H1) THEN
          H4=H2
          H2=H1
          H1=H4
          Z4=Z2
          Z2=Z1
          Z1=Z4
      ENDIF
      IF(H3.LT.H1) THEN
          H4=H3
          H3=H1
          H1=H4
          Z4=Z3
          Z3=Z1
          Z1=Z4
      ENDIF
      IF(H3.LT.H2) THEN
          H4=H3
          H3=H2
          H2=H4
          Z4=Z3
          Z3=Z2
          Z2=Z4
      ENDIF
      IWS=0
      DO 104 K=1,2
      CV=CV1*(2-K)+CV2*(K-1)
      IF(H1-CV) 1, 3, 3
    1 IF(H2-CV) 2,10,10
    2 IF(H3-CV)30,20,20
    3 IWA=0
      GOTO 104
   10 IWA=1
      ZA=D(H3,H1,Z3,Z1)
      ZB=D(H1,H2,Z1,Z2)
      IF(K.EQ.1) GOTO 60
      IF(IWS.EQ.0) THEN
          X(1)= REAL(Z1)
          X(2)= REAL(ZA)
          X(3)= REAL(ZB)
          Y(1)=AIMAG(Z1)
          Y(2)=AIMAG(ZA)
          Y(3)=AIMAG(ZB)
          NP=3
      ELSEIF(IWS.EQ.1) THEN
          X(1)= REAL(ZA)
          X(2)= REAL(ZB)
          X(3)= REAL(ZD)
          X(4)= REAL(ZC)
          Y(1)=AIMAG(ZA)
          Y(2)=AIMAG(ZB)
          Y(3)=AIMAG(ZD)
          Y(4)=AIMAG(ZC)
          NP=4
      ENDIF
      GOTO 60
   20 IWA=2
      ZA=D(H2,H3,Z2,Z3)
      ZB=D(H3,H1,Z3,Z1)
      IF(K.EQ.1) GOTO 60
      IF(IWS.EQ.0) THEN
          X(1)= REAL(ZA)
          X(2)= REAL(ZB)
          X(3)= REAL(Z1)
          X(4)= REAL(Z2)
          Y(1)=AIMAG(ZA)
          Y(2)=AIMAG(ZB)
          Y(3)=AIMAG(Z1)
          Y(4)=AIMAG(Z2)
          NP=4
      ELSEIF(IWS.EQ.1) THEN
          X(1)= REAL(ZA)
          X(2)= REAL(ZB)
          X(3)= REAL(ZC)
          X(4)= REAL(ZD)
          X(5)= REAL(Z2)
          Y(1)=AIMAG(ZA)
          Y(2)=AIMAG(ZB)
          Y(3)=AIMAG(ZC)
          Y(4)=AIMAG(ZD)
          Y(5)=AIMAG(Z2)
          NP=5
      ELSEIF(IWS.EQ.2) THEN
          X(1)= REAL(ZA)
          X(2)= REAL(ZB)
          X(3)= REAL(ZD)
          X(4)= REAL(ZC)
          Y(1)=AIMAG(ZA)
          Y(2)=AIMAG(ZB)
          Y(3)=AIMAG(ZD)
          Y(4)=AIMAG(ZC)
          NP=4
      ENDIF
      GOTO 60
   30 IWA=3
      IF(K.EQ.1) GOTO 103
      IF(IWS.EQ.0) THEN
          X(1)= REAL(Z1)
          X(2)= REAL(Z2)
          X(3)= REAL(Z3)
          Y(1)=AIMAG(Z1)
          Y(2)=AIMAG(Z2)
          Y(3)=AIMAG(Z3)
          NP=3
      ELSEIF(IWS.EQ.1) THEN
          X(1)= REAL(ZC)
          X(2)= REAL(ZD)
          X(3)= REAL(Z2)
          X(4)= REAL(Z3)
          Y(1)=AIMAG(ZC)
          Y(2)=AIMAG(ZD)
          Y(3)=AIMAG(Z2)
          Y(4)=AIMAG(Z3)
          NP=4
      ELSEIF(IWS.EQ.2) THEN
          X(1)= REAL(ZC)
          X(2)= REAL(ZD)
          X(3)= REAL(Z3)
          Y(1)=AIMAG(ZC)
          Y(2)=AIMAG(ZD)
          Y(3)=AIMAG(Z3)
          NP=3
      ENDIF
   60 IF(K.NE.1) CALL XFILAREA(X,Y,NP)
      IF(IWA.EQ.3) GOTO 103
      ZC=ZA
      ZD=ZB
      IWS=IWA
  104 CONTINUE
  103 CONTINUE
  102 CONTINUE
  101 CONTINUE
      RETURN
      END

      SUBROUTINE ZCONTC(ZG,Z,IWRK,X,Y,MD,MG,JG,C1,C2)
C   filling colour between two contour values
      DIMENSION ZG(MD,*),Z(MD,*),IWRK(MG,*),X(*),Y(*)
      COMPLEX   Z,B1,B2,ZA,ZB,Z1,Z2,Z3,Z4,D
                D(P1,P2,B1,B2   )=B1+(CV-P1)*(B2-B1)/(P2-P1)
      MGP=MG+1
      CV1=MIN(C1,C2)
      CV2=MAX(C1,C2)
      DO 3 J=1,JG-1
      HMN=MIN(ZG(1,J),ZG(1,J+1))
      HMX=MAX(ZG(1,J),ZG(1,J+1))
      DO 50 I=2,MG
      HMN=MIN(HMN,ZG(I,J),ZG(I,J+1))
   50 HMX=MAX(HMX,ZG(I,J),ZG(I,J+1))
      IF(HMN.GE.CV1.AND.HMX.LE.CV2) THEN
          NP=1
          ZA=Z(1,J)
          X(NP)= REAL(ZA)
          Y(NP)=AIMAG(ZA)
          DO 51 I=2,MG
          IF(Z(I,J).NE.ZA) THEN
              NP=NP+1
              ZA=Z(I,J)
              X(NP)= REAL(ZA)
              Y(NP)=AIMAG(ZA)
          ENDIF
   51     CONTINUE
          DO 52 I=1,MG
          IF(Z(MG-I+1,J+1).NE.ZA) THEN
              NP=NP+1
              ZA=Z(MG-I+1,J+1)
              X(NP)= REAL(ZA)
              Y(NP)=AIMAG(ZA)
          ENDIF
   52     CONTINUE
          CALL XFILAREA(X,Y,NP)
      ELSEIF(HMN.GT.CV2.OR.HMX.LT.CV1) THEN
          GOTO 3
      ENDIF
      DO 4 JJ=1,4
      DO 4 I=1,MG
    4 IWRK(I,JJ)=0
      DO 2 I=1,MG-1
      I1=I
      J1=J
      I2=I+1
      J2=J
      I3=I+1
      J3=J+1
      I4=I
      J4=J+1
      ISIG= 1
   29 CV=0.5*(1+ISIG)*CV1+0.5*(ISIG-1)*CV2
      H1=ZG(I1,J1)*ISIG
      H2=ZG(I2,J2)*ISIG
      H3=ZG(I3,J3)*ISIG
      H4=ZG(I4,J4)*ISIG
      IF(H1-CV)31,36,36
   31 IF(H2-CV)32,34,34
   32 IF(H3-CV)33,35,35
   33 IF(H4-CV)46,42,42
   34 IF(H3-CV)44,35,35
   35 IF(H4-CV)43,42,42
   36 IF(H2-CV)41,37,37
   37 IF(H3-CV)44,38,38
   38 IF(H4-CV)43,46,46
   41 ISW=1
      I10=I2
      J10=J2
      I40=I1
      J40=J1
      GOTO 45
   42 ISW=2
      I10=I1
      J10=J1
      I40=I4
      J40=J4
      GOTO 45
   43 ISW=3
      I10=I4
      J10=J4
      I40=I3
      J40=J3
      GOTO 45
   44 ISW=4
      I10=I3
      J10=J3
      I40=I2
      J40=J2
      GOTO 45
   46 IF(ISIG.EQ.1) THEN
         IF(H1.LT.CV ) GOTO 2
         ISIG=-1
         GOTO 29
      ENDIF
      GOTO 2
   45 J0=J+ISIG-2
      IF(IWRK(I10,J10-J0).EQ.1.AND.IWRK(I40,J40-J0).EQ.1) GOTO 2
      I1=I10
      J1=J10
      I4=I40
      J4=J40
      CV=0.5*(1+ISIG)*CV1+0.5*(ISIG-1)*CV2
      H1=ZG(I1,J1)*ISIG
      H4=ZG(I4,J4)*ISIG
      Z1= Z(I1,J1)
      Z4= Z(I4,J4)
      ZA=D(H4,H1,Z4,Z1)
      NP=1
      X(NP)=REAL(ZA)
      Y(NP)=AIMAG(ZA)
  101 I2=I1+MOD(ISW-1,2)*(1-2*((ISW-1)/2))
      J2=J1+MOD(ISW  ,2)*(1-2*((ISW-1)/2))
      I3=I4+MOD(ISW-1,2)*(1-2*((ISW-1)/2))
      J3=J4+MOD(ISW  ,2)*(1-2*((ISW-1)/2))
      IF(I2.EQ.0.OR.I3.EQ.0.OR.I2.EQ.MGP.OR.I3.EQ.MGP)GOTO 103
      IF(J2.EQ.J-1.OR.J3.EQ.J-1.OR.J2.EQ.J+2.OR.J3.EQ.J+2)GOTO 103
      GOTO 104
  103 ISW=MOD(ISW+1,4)+1
      KORNER=0
  112 INI=MOD(ISW  ,2)*(1-2*(MOD(ISW,4)/2))
      INJ=MOD(ISW+1,2)*(1-2*(MOD(ISW,4)/2))
      CVC=0.5*(ISIG+1)*CV2+0.5*(ISIG-1)*CV1
      H4=ZG(I4,J4)*ISIG
      IF(KORNER.EQ.0.AND.H4.GT.CVC) THEN
          ISIG=-ISIG
          CV=0.5*(1+ISIG)*CV1+0.5*(ISIG-1)*CV2
          I4=I1
          J4=J1
          I1=I4+INI
          J1=J4+INJ
          H1=ZG(I1,J1)*ISIG
          H4=ZG(I4,J4)*ISIG
          Z1= Z(I1,J1)
          Z4= Z(I4,J4)
          ZA=D(H4,H1,Z4,Z1)
          IWRK(I1,J1-J-ISIG+2)=1
          IWRK(I4,J4-J-ISIG+2)=1
          NP=NP+1
          X(NP)=REAL(ZA)
          Y(NP)=AIMAG(ZA)
          IF (I1.NE.I10.OR.J1.NE.J10.OR.I4.NE.I40.OR.J4.NE.J40) GOTO 101
          GOTO 100
      ELSE
          I1=I4
          J1=J4
          NP=NP+1
          X(NP)=REAL(Z(I1,J1))
          Y(NP)=AIMAG(Z(I1,J1))
      ENDIF
  111 I4=I1
      J4=J1
      I1=I4+INI
      J1=J4+INJ
      IF(I1.EQ.0.OR.I1.EQ.MGP.OR.J1.EQ.J-1.OR.J1.EQ.J+2)GOTO 113
      GOTO 114
  113 ISW=MOD(ISW+2,4)+1
      KORNER=1
      GOTO 112
  114 H1=ZG(I1,J1)*ISIG
      IF(H1.GT.CV2*(1+ISIG)*0.5+CV1*(ISIG-1)*0.5) THEN
          ISIG=-ISIG
          H1=-H1
          CV=0.5*(1+ISIG)*CV1+0.5*(ISIG-1)*CV2
      ELSEIF(H1-CV.GE.0.) THEN
          NP=NP+1
          X(NP)=REAL(Z(I1,J1))
          Y(NP)=AIMAG(Z(I1,J1))
          GOTO 111
      ENDIF
      H4=ZG(I4,J4)*ISIG
      Z1= Z(I1,J1)
      Z4= Z(I4,J4)
      ZA=D(H4,H1,Z4,Z1)
      IWRK(I1,J1-J-ISIG+2)=1
      IWRK(I4,J4-J-ISIG+2)=1
      NP=NP+1
      X(NP)=REAL(ZA)
      Y(NP)=AIMAG(ZA)
      IF (I1.NE.I10.OR.J1.NE.J10.OR.I4.NE.I40.OR.J4.NE.J40) GOTO 101
      GOTO 100
  104 H1=ZG(I1,J1)*ISIG
      H2=ZG(I2,J2)*ISIG
      H3=ZG(I3,J3)*ISIG
      H4=ZG(I4,J4)*ISIG
      H5=0.25*(H1+H2+H3+H4)
      Z1= Z(I1,J1)
      Z2= Z(I2,J2)
      Z3= Z(I3,J3)
      Z4= Z(I4,J4)
      IF(H1-CV) 11, 2,15
   11 IF(H2-CV) 12,13,13
   12 IF(H3-CV) 23,22,22
   13 IF(H3-CV) 14,21,21
   14 IF(H5-CV) 23,21,21
   15 IF(H2-CV) 16,16,18
   16 IF(H3-CV) 21,21,17
   17 IF(H5-CV) 21,21,23
   18 IF(H3-CV) 22,22,23
   21 ISA=1
      ZB=D(H1,H2,Z1,Z2)
      I4=I2
      J4=J2
      GOTO 30
   22 ISA=2
      ZB=D(H2,H3,Z2,Z3)
      I1=I2
      J1=J2
      I4=I3
      J4=J3
      GOTO 30
   23 ISA=3
      ZB=D(H3,H4,Z3,Z4)
      I1=I3
      J1=J3
   30 IF(ZB.NE.ZA) THEN
          NP=NP+1
          X(NP)=REAL(ZB)
          Y(NP)=AIMAG(ZB)
      ENDIF
      IWRK(I1,J1-J-ISIG+2)=1
      IWRK(I4,J4-J-ISIG+2)=1
      ZA=ZB
      ISW=MOD(ISW-ISA+5,4)+1
      IF (I1.NE.I10.OR.J1.NE.J10.OR.I4.NE.I40.OR.J4.NE.J40) GOTO 101
  100 CALL XFILAREA(X,Y,NP-1)
    2 CONTINUE
    3 CONTINUE
      RETURN
      END

      subroutine xcfill(z,x,y,iw,xw,yw,md,m,n,cl,ncl,mode)
C
c Produce a colour or gray filling map by calling Xcontc.
C Z : array defining the data surface.
C X,Y: 2-d arrays defining the coordintes of grid points.
C IW: integeter working space of size M*N.
C xw,yw: real working space of size 8*M.
C Md, M, N: dimension of arrays
C CL, NCl: input 1-d array containing the value of contours between
C          shaded areas.
C Mode: =0, no shading.
C if mode>0, colour varies from white to black for lower to higher
C            contour values.
C if mode<0, the order of gray filling is reveresed, i.e. black for
C            minimum, white for maximum.
C ZFmin, ZFmax: Lower and upper limits of the contour values between
C               which the shading is done.
C if zfmin=-999.0, zfmin=min(z), if zfmax=999.0, zfmax=max(z).
C
      dimension z(md,*),x(md,*),y(md,*),iw(md,*),xw(*),yw(*),cl(1)

c     data limzf,zfmin,zfmax /0,-999.0, 999.0/
c     save limzf, zfmax, zfmin

      common /xlimzf/ limzf, zfmax, zfmin

      integer icontcopt
      common /xcontc_opt/ icontcopt

      if(mode.eq.0) return
      zmax1=z(1,1)
      zmin1=zmax1
      do 1 j=1,n
      do 1 i=1,m
        zmax1=max(z(i,j),zmax1)
        zmin1=min(z(i,j),zmin1)
1     continue

      ncl1=ncl
      if(limzf.ne.0) then
        if(zfmin.ne.-999.0)then
        if(zfmin.ge.cl(ncl1))then
          ncl1=1
          return
        elseif(zfmin.lt.cl(1)) then
          do 14 k=ncl1,1,-1
14        cl(k+1)=cl(k)
          cl(1)=zmin1
          ncl1=ncl1+1
          goto 12
        endif
        do 10 k=1,ncl1-1
          if(cl(k).le.zfmin.and.cl(k+1).gt.zfmin) then
            cl(1)=zfmin
            ncl1=ncl1-k+1
            do 11 kk=2,ncl1
11            cl(kk)=cl(k+kk-1)
            goto 12
          endif
10      continue
        endif
12      if(zfmax.ne.999.0)then
        if(zfmax.le.cl(1))then
          ncl1=1
          return
        elseif(zfmax.gt.cl(ncl1)) then
          ncl1=ncl1+1
          cl(ncl1)=zmax1
          goto 22
        endif
        do 20 k=1,ncl1-1
          if(cl(k).lt.zfmax.and.cl(k+1).ge.zfmax) then
            cl(k+1)=zfmax
            ncl1=k+1
            goto 22
          endif
20      continue
        endif
22      continue
      endif
      if((limzf.eq.0.or.zfmin.eq.-999.0).and.(zmin1.lt.cl(1)))then
        do 32 k=ncl1,1,-1
32        cl(k+1)=cl(k)
        cl(1)=zmin1
        ncl1=ncl1+1
      endif
      if((limzf.eq.0.or.zfmax.eq.999.0).and.(zmax1.gt.cl(ncl1)))then
        ncl1=ncl1+1
        cl(ncl1)=zmax1
      endif

      call xqthik(ithick)
      call xthick(0)

      do 50 k=1,ncl1-1


        cl1=cl(k)
        cl2=cl(k+1)
C
C set gray degree:
C
        if(mode.gt.0) gray=(ncl1-k-0.5)/(ncl1-1)
        if(mode.lt.0) gray=(k-0.5)/(ncl1-1)
        call PSgray(gray)

        if( icontcopt.eq.1) then
          call xcontc(z,x,y,iw,xw,yw,md,m,n,cl1,cl2)
        else
          CALL XCONTC1(z,x,y,md,m,n,cl1,cl2)
        endif

 50   continue
      call xthick(ithick)
      call PSgray(0.0)
      return
      end

      subroutine xcflim(zfmi_1, zfma_1)
      common /xlimzf/ limzf, zfmax, zfmin
      limzf=1
      zfmax=zfma_1
      zfmin=zfmi_1
      return
      end

      subroutine xcontc(zg,zr,zi,iwrk,x,y,md,mg,jg,c1,c2)
c
C This is the 'real' version of ZCONTC.F
c Fill in colour between two contour values. developed by zunjun zhang
c at Reading university, england. jan, 1988.
c Modified by Shian-Jiann Lin at University of Oklahoma, Mar. 3,1990.
c external routine called: XFILAREA(x,y,np)
c Polygon is defined by (x(i),y(i),i=1,np)
c
      dimension zg(md,*),Zr(md,*),Zi(md,*)
      integer iwrk(mg,*)
      real x(*),y(*)  ! at least 8*m
c
c Single precision may cause the infinite loops noted above.  It fixes one
c case on LEMIEUX.
c
C Change to double is one way to fix the problem. If we find it does not
C work later, it can be changed to normalied fix just as xcontj does.
C                                            -- Commented by WYH
      double precision h1,h2,h3,h4,h5,cv,cvc
      double precision p1,p2

      dr(p1,p2,b1r,b2r)=b1r+(cv-p1)*(b2r-b1r)/(p2-p1)
      di(p1,p2,b1i,b2i)=b1i+(cv-p1)*(b2i-b1i)/(p2-p1)
      mgp=mg+1
      if(c1.gt.c2) then
      cv1=c2
      cv2=c1
      else
      cv1=c1
      cv2=c2
      endif
      do 3 j=1,jg-1
      hmx=zg(1,j)
      hmn=zg(1,j+1)
      if(hmn.gt.hmx) then
        hmx=hmn
        hmn=zg(1,j)
      endif
C
        do 50 i=2,mg
        a=zg(i,j)
        b=zg(i,j+1)
        if(a.gt.b) then
          if(a.gt.hmx) hmx=a
          if(b.lt.hmn) hmn=b
        else
          if(b.gt.hmx) hmx=b
          if(a.lt.hmn) hmn=a
        endif
50    continue

      if(hmn.ge.cv1.and.hmx.le.cv2) then
        np=1
      x(np)=Zr(1,j)
      y(np)=Zi(1,j)
        do 51 i=2,mg
          if(Zr(i,j).ne.x(np).or.Zi(i,j).ne.y(np)) then
            np=np+1
c           if(np.gt.md*100) then
c             print*,'np =', np, ' exceeding ', md*100
c             stop
c           endif
            x(np)= Zr(i,j)
            y(np)= Zi(i,j)
          endif
   51   continue
        do 52 i=1,mg
          img=mg-i+1
          if(Zr(img,j+1).ne.x(np).or.Zi(img,j+1).ne.y(np)) then
            np=np+1
c           if(np.gt.md*100) then
c             print*,'np =', np, ' exceeding ', md*100
c             stop
c           endif
            x(np) = Zr(img,j+1)
            y(np) = Zi(img,j+1)
          endif
   52   continue
c       print*,'calling xfilarea 1'
        call XFILAREA(x,y,np)
c       print*,'done calling xfilarea'
      elseif(hmn.gt.cv2.or.hmx.lt.cv1) then
        goto 3
      endif

      do 4 jj=1,4
      do 4 i=1,mg
    4 iwrk(i,jj)=0
      do 2 i=1,mg-1
      i1=i
      j1=j
      i2=i+1
      j2=j
      i3=i+1
      j3=j+1
      i4=i
      j4=j+1
      isig= 1
   29 cv=0.5*(1+isig)*cv1+0.5*(isig-1)*cv2
      h1=zg(i1,j1)*isig
      h2=zg(i2,j2)*isig
      h3=zg(i3,j3)*isig
      h4=zg(i4,j4)*isig
      if(h1-cv)31,36,36
   31 if(h2-cv)32,34,34
   32 if(h3-cv)33,35,35
   33 if(h4-cv)46,42,42
   34 if(h3-cv)44,35,35
   35 if(h4-cv)43,42,42
   36 if(h2-cv)41,37,37
   37 if(h3-cv)44,38,38
   38 if(h4-cv)43,46,46
   41 isw=1
      i10=i2
      j10=j2
      i40=i1
      j40=j1
      goto 45
   42 isw=2
      i10=i1
      j10=j1
      i40=i4
      j40=j4
      goto 45
   43 isw=3
      i10=i4
      j10=j4
      i40=i3
      j40=j3
      goto 45
   44 isw=4
      i10=i3
      j10=j3
      i40=i2
      j40=j2
      goto 45
   46 if(isig.eq.1) then
        if(h1.lt.cv ) goto 2
        isig=-1
        goto 29
      endif
      goto 2
   45 j0=j+isig-2
      if(iwrk(i10,j10-j0).eq.1.and.iwrk(i40,j40-j0).eq.1) goto 2
      i1=i10
      j1=j10
      i4=i40
      j4=j40
      cv=0.5*(1+isig)*cv1+0.5*(isig-1)*cv2
      h1=zg(i1,j1)*isig
      h4=zg(i4,j4)*isig

      np=1
      x(np)=dr(h4,h1,Zr(i4,j4),Zr(i1,j1))
      y(np)=di(h4,h1,Zi(i4,j4),Zi(i1,j1))
C
  101 i2=i1+mod(isw-1,2)*(1-2*((isw-1)/2))
      j2=j1+mod(isw  ,2)*(1-2*((isw-1)/2))
      i3=i4+mod(isw-1,2)*(1-2*((isw-1)/2))
      j3=j4+mod(isw  ,2)*(1-2*((isw-1)/2))
      if(i2.eq.0.or.i3.eq.0.or.i2.eq.mgp.or.i3.eq.mgp)goto 103
      if(j2.eq.j-1.or.j3.eq.j-1.or.j2.eq.j+2.or.j3.eq.j+2)goto 103
      goto 104
  103 isw=mod(isw+1,4)+1
      korner=0
  112 ini=mod(isw  ,2)*(1-2*(mod(isw,4)/2))
      inj=mod(isw+1,2)*(1-2*(mod(isw,4)/2))
      cvc=0.5*(isig+1)*cv2+0.5*(isig-1)*cv1
      h4=zg(i4,j4)*isig
      if(korner.eq.0.and.h4.gt.cvc) then
        isig=-isig
        cv=0.5*(1+isig)*cv1+0.5*(isig-1)*cv2
        i4=i1
        j4=j1
        i1=i4+ini
        j1=j4+inj
        h1=zg(i1,j1)*isig
        h4=zg(i4,j4)*isig
        iwrk(i1,j1-j-isig+2)=1
        iwrk(i4,j4-j-isig+2)=1
        np=np+1
c       if(np.gt.md*100) then
c         print*,'np =', np, ' exceeding ', md*100
c         stop
c       endif
      x(np)=dr(h4,h1,Zr(i4,j4),Zr(i1,j1))
      y(np)=di(h4,h1,Zi(i4,j4),Zi(i1,j1))
        if (i1.ne.i10.or.j1.ne.j10.or.i4.ne.i40.or.j4.ne.j40) goto 101
        goto 100
      else
        i1=i4
        j1=j4
        np=np+1
c       if(np.gt.md*100) then
c         print*,'np =', np, ' exceeding ', md*100
c         stop
c       endif
      x(np)=Zr(i1,j1)
      y(np)=Zi(i1,j1)
      endif
  111 i4=i1
      j4=j1
      i1=i4+ini
      j1=j4+inj
      if(i1.eq.0.or.i1.eq.mgp.or.j1.eq.j-1.or.j1.eq.j+2)goto 113
      goto 114
  113 isw=mod(isw+2,4)+1
      korner=1
      goto 112
  114 h1=zg(i1,j1)*isig
      if(h1.gt.cv2*(1+isig)*0.5+cv1*(isig-1)*0.5) then
        isig=-isig
        h1=-h1
        cv=0.5*(1+isig)*cv1+0.5*(isig-1)*cv2
      elseif(h1-cv.ge.0.) then
        np=np+1
c       if(np.gt.md*100) then
c         print*,'np =', np, ' exceeding ', md*100
c         stop
c       endif
      x(np)=Zr(i1,j1)
      y(np)=Zi(i1,j1)
        goto 111
      endif
      h4=zg(i4,j4)*isig
      iwrk(i1,j1-j-isig+2)=1
      iwrk(i4,j4-j-isig+2)=1
      np=np+1
c     if(np.gt.md*100) then
c       print*,'np =', np, ' exceeding ', md*100
c       stop
c     endif
      x(np)=dr(h4,h1,Zr(i4,j4),Zr(i1,j1))
      y(np)=di(h4,h1,Zi(i4,j4),Zi(i1,j1))
      if (i1.ne.i10.or.j1.ne.j10.or.i4.ne.i40.or.j4.ne.j40) goto 101
      goto 100
  104 h1=zg(i1,j1)*isig
      h2=zg(i2,j2)*isig
      h3=zg(i3,j3)*isig
      h4=zg(i4,j4)*isig
      h5=0.25*(h1+h2+h3+h4)
      if(h1-cv) 11, 2,15
   11 if(h2-cv) 12,13,13
   12 if(h3-cv) 23,22,22
   13 if(h3-cv) 14,21,21
   14 if(h5-cv) 23,21,21
   15 if(h2-cv) 16,16,18
   16 if(h3-cv) 21,21,17
   17 if(h5-cv) 21,21,23
   18 if(h3-cv) 22,22,23
   21 isa=1
      zbr=dr(h1,h2,Zr(i1,j1),Zr(i2,j2))
      zbi=di(h1,h2,Zi(i1,j1),Zi(i2,j2))
c
      i4=i2
      j4=j2
      goto 30
   22 isa=2
      zbr=dr(h2,h3,Zr(i2,j2),Zr(i3,j3))
      zbi=di(h2,h3,Zi(i2,j2),Zi(i3,j3))
c
      i1=i2
      j1=j2
      i4=i3
      j4=j3
      goto 30
   23 isa=3
      zbr=dr(h3,h4,Zr(i3,j3),Zr(i4,j4))
      zbi=di(h3,h4,Zi(i3,j3),Zi(i4,j4))
C
      i1=i3
      j1=j3
30    if(zbr.ne.x(np).or.zbi.ne.y(np)) then
        np=np+1
c       if(np.gt.md*100) then
c         print*,'np =', np, ' exceeding ', md*100
c         stop
c       endif
        x(np)=zbr
        y(np)=zbi
      endif
      iwrk(i1,j1-j-isig+2)=1
      iwrk(i4,j4-j-isig+2)=1
      isw=mod(isw-isa+5,4)+1
      if (i1.ne.i10.or.j1.ne.j10.or.i4.ne.i40.or.j4.ne.j40) goto 101
  100 continue
c     print*,'calling xfilarea 2'
      call XFILAREA(x,y,np-1)
c     print*,'done calling xfilarea'
    2 continue
    3 continue
      return
      end

C* SURFACE VIEWING ROUTINES
      SUBROUTINE ZSFPLT(SURFAS,MD,M,N,WORK)
C     isometric surface viewing
      REAL SURFAS(MD,*),WORK(2,*)
      SAVE IQ,ANGISM,SCALE,XRANGE,XX0,YY0,MODE
      DATA IQ/1/,ANGISM/1.043862/,SCALE/1./,XX0,YY0,XRANGE/2*0.,2./
     :    ,MODE/0/
C        MD: the first dimension of array to be viewed
C         M: array dimension of x direction
C         N: array dimension of y direction
C      WORK: working space of total dimension at least 2*max(M,N)
C
      IQUAD=IQ
      EPS=1.E-4
      M1=(MOD(IQUAD,4)/2)*(M-1)+1
      N1=((IQUAD-1)/2)*(N-1)   +1
      DX=XRANGE/(M+N)
      DY=DX/TAN(ANGISM)
      IF( MOD(IQUAD,2).EQ.0) DX=-ABS(DX)
      HS=   SCALE
C     plot parallel along x direction
C     first line
C     loop 101,102 draw line in x-direction
      IF(MODE.EQ.2) GOTO 301
      DO 101 I=1,M
      X=(I-1)*DX+XX0
      Y=(I-1)*DY+YY0
      II=ABS(I-M1)+1
      JJ=ABS(1-N1)+1
      Y=Y+SURFAS(II,JJ)*HS
      WORK(1,I)=Y
      WORK(2,I)=Y
      IF(I.EQ.1) THEN
        CALL XPENUP(X,Y)
      ELSE
        CALL XPENDN(X,Y)
      ENDIF
  101 CONTINUE
      DO 102 J=2,N
C     loop 103 : upper surface
      DO 103 I=M,2,-1
        X= (I-1)*DX-(J-1)*DX +XX0
        Y0= (I-1)*DY+(J-1)*DY+YY0
        II=ABS(I-M1)+1
        JJ=ABS(J-N1)+1
        Y=Y0+SURFAS(II,JJ  )*HS
        PMAX1=WORK(1, I)
        PMAX2=WORK(1, I-1)
        IF(Y.GE.WORK(1,I-1)) THEN
          IPLOT=1
          WORK(1,I)=Y
          IF(I.EQ.M) THEN
            CALL ZSF002(X,Y,IPLOT)
          ELSE
            II1=ABS(I  -M1)+1
            JJ1=ABS(J-1-N1)+1
            Y2=Y0+SURFAS(II1,JJ1)*HS-DY
            IF(ABS(PMAX1-Y2).LT.EPS) Y1=Y2
            CALL ZSF001(X,Y,IPLOT,PMAX2,PMAX1,Y1,DX)
          ENDIF
        ELSE
          IPLOT=0
          IF(I.EQ.M) THEN
            CALL ZSF002(X,Y,0)
          ELSE
            CALL ZSF001(X,Y,IPLOT,PMAX2,PMAX1,Y1,DX)
          ENDIF
          WORK(1,I)=WORK(1,I-1)
        ENDIF
        Y1=Y
  103 CONTINUE
      X=-(J-1)*DX+XX0
      Y=+(J-1)*DY+YY0
      II=ABS(1-M1)+1
      JJ=ABS(J-N1)+1
      Y=Y+SURFAS(II,JJ)*HS
      WORK(1,1)=Y
      CALL ZSF003(X,Y)

C     loop 104 : Lower surface
      DO 104 I=M,2,-1
        X= (I-1)*DX-(J-1)*DX+XX0
        Y0=(I-1)*DY+(J-1)*DY+YY0
        II=ABS(I-M1)+1
        JJ=ABS(J-N1)+1
        Y=Y0+SURFAS(II,JJ)*HS
        PMIN1=WORK(2, I)
        PMIN2=WORK(2, I-1)
        IF(Y.LE.WORK(2,I-1)) THEN
          IPLOT=1
          WORK(2,I)=Y
          IF(I.EQ.M) THEN
            CALL ZSF002(X,Y,IPLOT)
          ELSE
            II1=ABS(I  -M1)+1
            JJ1=ABS(J-1-N1)+1
            Y2=Y0+SURFAS(II1,JJ1)*HS-DY
            IF(ABS(PMIN1-Y2).LT.EPS) Y1=Y2
            CALL ZSF001(X,Y,IPLOT,PMIN2,PMIN1,Y1,DX)
          ENDIF
        ELSE
          IPLOT=0
          IF(I.EQ.M) THEN
            CALL ZSF002(X,Y,IPLOT)
          ELSE
            CALL ZSF001(X,Y,IPLOT,PMIN2,PMIN1,Y1,DX)
          ENDIF
          WORK(2,I)=WORK(2,I-1)
        ENDIF
        Y1=Y
  104 CONTINUE
      X=-(J-1)*DX+XX0
      Y=+(J-1)*DY+YY0
      II=ABS(1-M1)+1
      JJ=ABS(J-N1)+1
      Y=Y+SURFAS(II,JJ)*HS
      WORK(2,1)=Y
      CALL ZSF003(X,Y)
  102 CONTINUE

C     loop 201,202 draw line in y-direction
  301 IF(MODE.EQ.1) RETURN
      DO 201 J=1,N
      X=-(J-1)*DX+XX0
      Y= (J-1)*DY+YY0
      II=ABS(1-M1)+1
      JJ=ABS(J-N1)+1
      Y=Y+SURFAS(II,JJ)*HS
      WORK(1,J)=Y
      WORK(2,J)=Y
      IF(J.EQ.1) THEN
        CALL XPENUP(X,Y)
      ELSE
        CALL XPENDN(X,Y)
      ENDIF
  201 CONTINUE

      DO 202 I=2,M
C     loop 203 : upper surface
      DO 203 J=N,2,-1
        X= (I-1)*DX-(J-1)*DX     +XX0
        Y0= (I-1)*DY+(J-1)*DY    +YY0
        II=ABS(I-M1)+1
        JJ=ABS(J-N1)+1
        Y=Y0+SURFAS(II,JJ  )*HS
        PMAX1=WORK(1, J)
        PMAX2=WORK(1, J-1)
        IF(Y.GE.WORK(1,J-1)) THEN
          IPLOT=1
          WORK(1,J)=Y
          IF(J.EQ.N) THEN
            CALL ZSF002(X,Y,IPLOT)
          ELSE
            II1=ABS(I-1-M1)+1
            JJ1=ABS(J  -N1)+1
            Y2=Y0+SURFAS(II1,JJ1)*HS-DY
            IF(ABS(PMAX1-Y2).LT.EPS) Y1=Y2
            CALL ZSF001(X,Y,IPLOT,PMAX2,PMAX1,Y1,-DX)
          ENDIF
        ELSE
          IPLOT=0
          IF(J.EQ.N) THEN
            CALL ZSF002(X,Y,0)
          ELSE
            CALL ZSF001(X,Y,IPLOT,PMAX2,PMAX1,Y1,-DX)
          ENDIF
          WORK(1,J)=WORK(1,J-1)
        ENDIF
        Y1=Y
  203 CONTINUE
      X=(I-1)*DX+XX0
      Y=(I-1)*DY+YY0
      II=ABS(I-M1)+1
      JJ=ABS(1-N1)+1
      Y=Y+SURFAS(II,JJ  )*HS
      WORK(1,1)=Y
      CALL ZSF003(X,Y)
C     loop 204 : lower surface
      DO 204 J=N,2,-1
        X= (I-1)*DX-(J-1)*DX       +XX0
        Y0= (I-1)*DY+(J-1)*DY      +YY0
        II=ABS(I-M1)+1
        JJ=ABS(J-N1)+1
        Y=Y0+SURFAS(II,JJ  )*HS
        PMIN1=WORK(2, J)
        PMIN2=WORK(2, J-1)
        IF(Y.LE.WORK(2,J-1)) THEN
          IPLOT=1
          WORK(2,J)=Y
          IF(J.EQ.N) THEN
            CALL ZSF002(X,Y,IPLOT)
          ELSE
            II1=ABS(I-1-M1)+1
            JJ1=ABS(J  -N1)+1
            Y2=Y0+SURFAS(II1,JJ1)*HS-DY
            IF(ABS(PMIN1-Y2).LT.EPS) Y1=Y2
            CALL ZSF001(X,Y,IPLOT,PMIN2,PMIN1,Y1,-DX)
          ENDIF
        ELSE
          IPLOT=0
          IF(J.EQ.N) THEN
            CALL ZSF002(X,Y,IPLOT)
          ELSE
            CALL ZSF001(X,Y,IPLOT,PMIN2,PMIN1,Y1,-DX)
          ENDIF
          WORK(2,J)=WORK(2,J-1)
        ENDIF
        Y1=Y
  204 CONTINUE
      X=(I-1)*DX+XX0
      Y=(I-1)*DY+YY0
      II=ABS(I-M1)+1
      JJ=ABS(1-N1)+1
      Y=Y+SURFAS(II,JJ  )*HS
      WORK(2,1)=Y
      CALL ZSF003(X,Y)
  202 CONTINUE
      RETURN

      ENTRY ZSFSTL(MODES)
C  MODE=0 draw surface lines along both axes directions (default).
C  MODE=1 draw surface lines along the x-direction.
C  MODE=2 draw surface lines along the y-direction.
      MODE=MODES
      RETURN

      ENTRY ZSFVEW(IQS)
C  Define the corner through which the surface is viewed
C     IQ=1,2,3 & 4
      IQ=IQS
      RETURN

      ENTRY ZSFSCL(SCALES)
C  Set the scaler which scales the surplot data
C  Decreasing SCALE results in decreasing in size of the plot
      SCALE=SCALES
      RETURN

      ENTRY ZSFANG(ANGS)
C  Set the isometric angle of the surface viewing
      ANGISM=4.*ATAN(1.)*ANGS/180.
      RETURN

      ENTRY ZSFLOC(X0S,Y0S,RANGEX)
C     left to right range in plotting space
C     (X0S,Y0S) : the position of reference point
C     i.e. the closest grid point to the viewer
      XRANGE=RANGEX
      XX0=X0S
      YY0=Y0S
      RETURN
      END

      SUBROUTINE ZSF001(X,Y,IPLOT,PY,PYP,YP,DX)
      SAVE ILAST
      IF(IPLOT.EQ.ILAST) THEN
        IF(IPLOT.GT.0) THEN
          CALL XPENDN(X,Y)
        ENDIF
      ELSE
        XX=(PY-Y)/(YP-PYP + PY-Y)
        YY=Y+(YP-Y)*XX
        XX=X+XX*DX
        IF(IPLOT.EQ.0) THEN
          CALL XPENDN(XX,YY)
          IPLOT=0
        ELSE
          CALL XPENUP(XX,YY)
          CALL XPENDN(X, Y )
          IPLOT=1
        ENDIF
      ENDIF
      ILAST=IPLOT
      RETURN

      ENTRY ZSF002(X,Y,IPLOT)
C      this entry is used to initialize the state of plotting line
          CALL XPENUP(X,Y)
      ILAST=IPLOT
      RETURN

      ENTRY ZSF003(X,Y)
      IF(ILAST.GT.0) THEN
        CALL XPENDN(X,Y)
      ELSE
        CALL XPENUP(X,Y)
      ENDIF
      RETURN
      END

      SUBROUTINE XCHDEC(ICDATA,CHDATA,I)
C To decode character set data for charactere No. I.
      CHARACTER CHDATA(127)*300,XCH*2,YCH*2,STR2*4
      INTEGER ICDATA (0:150, 32:127)
      ILEN = ICLENG(CHDATA(I))
      I0=0
      J0=0
      ICDATA(1,I)=ILEN/2+1
      ICD=1
      DO 5 INUM=1,ILEN,4
         STR2 = CHDATA(I) (INUM:INUM+3)
         READ(STR2,100)XCH,YCH
         JX = XDECOD(XCH)
         JY = XDECOD(YCH)
         I00=0
         IF (JX .GT. 127) THEN
            JX = JX - 128
            I00=1
         ENDIF
         I0=I0+JX-64
         J0=J0+JY-64
         IF( I00.EQ.1) THEN
           ICDATA(ICD+1,I)=I0
         ELSE
           ICDATA(ICD+1,I)=I0+50
         ENDIF
         ICDATA(ICD+2,I)=J0
         ICD=ICD+2
 5    CONTINUE
  100 FORMAT(2A,2A)
      RETURN
      END

      FUNCTION XDECOD(CH)
      CHARACTER   CH*2 ,F*1
      COMMON /XCHR30/ ICRAM(256)
  100 FORMAT(1A)
  101 FORMAT(1X,1A)
      READ(CH,100) F
      J = ICHAR(F)
      J=  ICRAM( J )
      IF (J .GT. 200) THEN
          J=J-240
      ELSE
          J=J-183
      END IF
      READ(CH,101) F
      K = ICHAR(F)
      K=  ICRAM( K )
      IF (K .GT. 200) THEN
         K = K - 240
      ELSE
         K=K-183
      ENDIF
      XDECOD = 16*J+K
      END
C
      SUBROUTINE XCSETA(C)
      CHARACTER*300 C(127)

           C(1) =' ROMAN CHARACTER SET.'
           C(2) ='20 30 Size of characters in x and y direction.'
           C(32)='5540'
           C(33)='4155BF3EC134C14CBF42403EC03A4035BF3FC13FC141BF414A3E'
           C(34)='4258BF3FBF41C141C13FC03EBF3EBF3F4845BF3FBF41C141C13FC0
     +3EBF3EBF3F492D'
           C(35)='4855B9244D5CB9243A51CE40313ACE40453C'
           C(36)  ='4559C023445DC0234456BF3FC13FC141C041BE42BD41BC40BD3F
     +BE3EC03EC13EC13FC23FC63EC23FC23E3249C23EC23FC63EC23FC13FC13EC03CBE
     +3EBD3FBC40BD41BE42C041C141C13FBF3F533D'
           C(37)  ='D2553340C23EC03EBF3EBE3FBE40BE42C042C142C241C240C23F
     +C33FC340C341C2413C32BE3FBF3EC03EC23EC240C241C142C042BE42BE404839'
           C(38)  ='524DBF3FC13FC141C041BF41BF40BF3FBF3EBE3BBE3DBE3EBE3F
     +BD40BD41BF42C043C142C644C242C142C042BF42BE41BE3FBF3EC03EC13DC23DC5
     +39C23EC33FC140C141C041313EBE41BF42C043C142C2424046C13EC835C23EC23F
     +4640'
           C(39)  ='4258BF3FBF41C141C13FC03EBF3EBF3F492D'
           C(40)  ='4759BE3EBE3DBE3CBF3BC03CC13BC23CC23DC23E3E5EBE3CBF3D
     +BF3BC03CC13BC13DC23C4945'
           C(41)  ='4059C23EC23DC23CC13BC03CBF3BBE3CBE3DBE3E425EC23CC13D
     +C13BC03CBF3BBF3DBE3C4C45'
           C(42)  ='4555C0343B49CA3A4046B63A5034'
           C(43)  ='4952C02E3749D2404837'
           C(44)  ='4241BF3FBF41C141C13FC03EBF3EBF3F4944'
           C(45)  ='4049D2404837'
           C(46)  ='4140BF41C141C13F493F'
           C(47)  ='5259AE205647'
           C(48)  ='4655BD3FBE3DBF3BC03DC13BC23DC33FC240C341C243C145C043
     +BF45BE43BD41BE40BE3FBF3FBF3EBF3BC03DC13BC13EC13FC23F4240C241C141C1
     +42C145C043BF45BF42BF414A2C'
           C(49)  ='4051C241C343C02B3F54C02C3C40C9404640'
           C(50)  ='4151C13FBF3FBF41C041C142C141C341C440C33FC13FC13EC03E
     +BF3EBD3EBB3EBE3FBE3EBF3DC03D4955C23FC13FC13EC03EBF3EBD3EBC3E3B39C1
     +41C240C53EC340C241C1413540C53DC440C141C142463D'
           C(51)  ='4151C13FBF3FBF41C041C142C141C341C440C33FC13EC03DBF3E
     +BD3FBD404349C23FC13EC03DBF3EBE3FC23FC23EC13EC03DBF3EBF3FBD3FBC40BD
     +41BF41BF42C041C141C13FBF3F4B46C13DC03DBF3EBF3F483F'
           C(52)  ='4A53C02D4140C055B531D040453A'
           C(53)  ='4255BE36C242C341C340C33FC23EC13DC03EBF3DBE3EBD3FBD40
     +BD41BF41BF42C041C141C13FBF3F474AC23FC23EC13DC03EBF3DBE3EBE3F3A55CA
     +40363FC540C541482B'
           C(54)  ='4C52BF3FC13FC141C041BF42BE41BD40BD3FBE3EBF3EBF3CC03A
     +C13DC23EC33FC240C341C242C143C041BF43BE42BD41BF40BD3FBE3EBF3D464EBE
     +3FBE3EBF3EBF3CC03AC13DC23EC23F4240C241C242C143C041BF43BE424A34'
           C(55  )='4055C03A4042C142C242C240C53DC240C141C142333EC241C240
     +C53E4443C03DBF3DBC3BBF3EBF3DC03B464FBB3BBF3EBF3DC03B4E40'
           C(56  )='4555BD3FBF3EC03DC13EC33FC440C341C142C043BF42BD41BC40
     +BE3FBF3EC03DC13EC23F4440C241C142C043BF42BE413C37BD3FBF3FBF3EC03CC1
     +3EC13FC33FC440C341C141C142C044BF42BF41BD413C40BE3FBF3FBF3EC03CC13E
     +C13FC23F4440C241C141C142C044BF42BF414935'
           C(57  )='4D4EBF3DBE3EBD3FBF40BD41BE42BF43C041C143C242C341C240
     +C33FC23EC13DC03ABF3CBF3EBE3EBD3FBD40BE41BF42C041C141C13FBF3F4445BE
     +41BE42BF43C041C143C242C2414240C23FC23EC13DC03ABF3CBF3EBE3E4B3F'
           C(58  )='414EBF3FC13FC141BF414034BF3FC13FC141BF414E3E'
           C(59  )='414EBF3FC13FC141BF414032BF41C141C13FC03EBF3EBF3F4E44
     +'
           C(60  )='5052B037D0374440'
           C(61  )='404CD2402E3AD240483A'
           C(62  )='4052D037B0375440'
           C(63  )='4151C13FBF3FBF41C041C142C141C241C340C33FC13FC13EC03E
     +BF3EBF3FBC3EC03D414EC23FC13FC13EC03EBF3EBE3E3E37BF3FC13FC141BF414D
     +3E'
           C(64  )='4F4DBF42BE41BD40BE3FBF3FBF3DC03DC13EC23FC340C241C142
     +3B48BE3EBF3DC03DC13EC13F474BBF38C03EC23FC240C242C143C042BF43BF42BE
     +42BE41BD41BD40BD3FBE3FBE3EBF3EBF3DC03DC13DC13EC23EC23FC33FC340C341
     +C241C1413E4DBF38C03EC13F483B'
           C(65  )='4955B92B4755C72B3952C62E3546C940333AC6404640C6404340
     +'
           C(66  )='4340C0554140C02B3C55CC40C33FC13FC13EC03EBF3EBF3FBD3F
     +404AC23FC13FC13EC03EBF3EBF3FBE3F3840C840C33FC13FC13EC03DBF3EBF3FBD
     +3FB4404C4BC23FC13FC13EC03DBF3EBF3FBE3F4A40'
           C(67  )='4E52C13DC046BF3DBE42BD41BE40BD3FBE3EBF3EBF3DC03BC13D
     +C13EC23EC33FC240C341C242C1423850BE3FBE3EBF3EBF3DC03BC13DC13EC23EC2
     +3F4E40'
           C(68  )='4340C0554140C02B3C55CA40C33FC23EC13EC13DC03BBF3DBF3E
     +BE3EBD3FB6404A55C23FC23EC13EC13DC03BBF3DBF3EBE3EBE3F4C40'
           C(69  )='4340C0554140C02B464FC038364ED040C03ABF463536C6403635
     +D040C046BF3A4640'
           C(70  )='4340C0554140C02B464FC038364ED040C03ABF463536C6403635
     +C7404D40'
           C(71  )='4E52C13DC046BF3DBE42BD41BE40BD3FBE3EBF3EBF3DC03BC13D
     +C13EC23EC33FC240C341C2423952BE3FBE3EBF3EBF3DC03BC13DC13EC23EC23F47
     +48C0384148C0383C48C7404538'
           C(72  )='4340C0554140C02B4C40C0554140C02B2F55C7404640C7403036
     +CC403035C7404640C7404440'
           C(73  )='4340C0554140C02B3C55C740392BC7404440'
           C(74  )='4855C02FBF3DBE3FBE40BE41BF42C042C141C13FBF3F4651C02F
     +BF3DBF3F3F55C740442B'
           C(75  )='4340C0554140C02B4D55B3334544C834374CC8343055C7404640
     +C6402D2BC7404640C6404340'
           C(76  )='4340C0554140C02B3C55C740392BCF40C046BF3A4440'
           C(77  )='4340C0554140C62E3952C72BC755C02B4140C0552E40C4404D40
     +C4402B2BC6404840C7404440'
           C(78  )='4340C0554140CC2D3451CC2DC0553040C4404940C6402D2BC640
     +5140'
           C(79  )='4755BD3FBE3EBF3EBF3CC03DC13CC13EC23EC33FC240C341C242
     +C142C144C043BF44BF42BE42BD41BE40BE3FBE3EBF3EBF3CC03DC13CC13EC23EC2
     +3F4240C241C242C142C144C043BF44BF42BE42BE414D2B'
           C(80  )='4340C0554140C02B3C55CC40C33FC13FC13EC03DBF3EBF3FBD3F
     +B840484BC23FC13FC13EC03DBF3EBF3FBE3F3436C7404F40'
           C(81  )='4755BD3FBE3EBF3EBF3CC03DC13CC13EC23EC33FC240C341C242
     +C142C144C043BF44BF42BE42BD41BE40BE3FBE3EBF3EBF3CC03DC13CC13EC23EC2
     +3F4240C241C242C142C144C043BF44BF42BE42BE413B2DC041C142C241C140C23F
     +C13EC139C13FC240C142C0413B45C13CC13EC13FC140C1414643'
           C(82  )='4340C0554140C02B3C55CC40C33FC13FC13EC03EBF3EBF3FBD3F
     +B840484AC23FC13FC13EC03EBF3EBF3FBE3F3435C740424BC23FC13FC339C13FC1
     +40C1413948C13EC239C13FC240C142C041443D'
           C(83  )='4D52C143C03ABF43BE42BD41BD40BD3FBE3EC03EC13EC13FC23F
     +C63EC23FC23E3249C23EC23FC63EC23FC13FC13EC03CBE3EBD3FBD40BD41BE42BF
     +43C03AC143533D'
           C(84  )='4740C0554140C02B3955BF3AC046CF40C03ABF46362BC7404840
     +'
           C(85  )='4355C031C13DC23EC33FC240C341C242C143C04F3340C031C13D
     +C23EC23F3755C7404740C640442B'
           C(86  )='4255C72B3A55C62E4752B92B3755C6404640C640422B'
           C(87  )='4355C42B3D55C3304450BC2B4455C42B3D55C3304450BC2B3155
     +C7404940C640422B'
           C(88  )='4255CD2B3455CD2B4055B22B3E55C6404640C6402E2BC6404640
     +C6404240'
           C(89  )='4255C735C0363A55C735C0364755B935364BC6404740C640332B
     +C7404840'
           C(90  )='CD554140B32B4055BF3AC046CE40322BCE40C046BF3A4740'
           C(91  )='4059C0204140C0603F40C7403920C7404847'
           C(92  )='4059D2204447'
           C(93  )='4659C0204140C0603940C7403920C7404947'
           C(94  )='4052C747C739462E'
           C(95  )='403ED4404442'
           C(96  )='4058C13FC141BF41BF3FC03EC13EC13F472D'
           C(97  )='424CC03FBF40C041C141C241C440C23FC13FC13EC039C13EC13F
     +3D4CC037C13EC23FC1403C4ABF3FBA3FBD3FBF3EC03EC13EC33FC340C241C24239
     +45BE3FBF3EC03EC13EC23F5040'
           C(98  )='4340C0554140C02B404BC242C241C240C33FC23EC13DC03EBF3D
     +BE3EBD3FBE40BE41BE42464BC23FC23EC13DC03EBF3DBE3EBE3F3655C440512B'
           C(99  )='4C4BBF3FC13FC141C041BE42BE41BD40BD3FBE3EBF3DC03EC13D
     +C23EC33FC240C341C242394BBE3FBE3EBF3DC03EC13DC23EC23F4D40'
           C(100 )='4C55C02B4155C02B3F4BBE42BE41BE40BD3FBE3EBF3DC03EC13D
     +C23EC33FC240C241C2423A4BBE3FBE3EBF3DC03EC13DC23EC23F4355C4403F2BC4
     +404540'
           C(101 )='4148CC40C042BF42BF41BE41BD40BD3FBE3EBF3DC03EC13DC23E
     +C33FC240C341C2423F45C043BF423B41BE3FBE3EBF3DC03EC13DC23EC23F4D40'
           C(102 )='4854BF3FC13FC141C041BF41BE40BE3FBF3EC02E4355BF3FBF3E
     +C02E3C4EC8403832C7404640'
           C(103 )='464EBE3FBF3FBF3EC03EC13EC13FC23FC240C241C141C142C042
     +BF42BF41BE41BE403E3FBF3EC03CC13E4640C142C044BF42413FC141C241C03FBE
     +403739BF3FBF3EC03FC13EC33FC540C33FC13F3345C13FC33FC540C33FC13EC03F
     +BF3EBD3FBA40BD41BF42C041C142C3414F40'
           C(104 )='4340C0554140C02B404BC242C341C240C33FC13EC0353C4EC23F
     +C13EC0353255C4403C2BC7404440C7404440'
           C(105 )='4355BF3FC13FC141BF414039C032414EC0323C4EC4403C32C740
     +4440'
           C(106 )='4555BF3FC13FC141BF414139C02EBF3EBE3FBE40BF41C041C141
     +C13FBF3F4454C02EBF3EBF3F3F55C4404532'
           C(107 )='4340C0554140C02B4A4EB6364544C6383948C6383255C4404739
     +C6402F32C7404440C6404440'
           C(108 )='4340C0554140C02B3C55C4403C2BC7404440'
           C(109 )='434EC032414EC032404BC242C341C240C33FC13EC0353C4EC23F
     +C13EC035414BC242C341C240C33FC13EC0353C4EC23FC13EC035274EC4403C32C7
     +404440C7404440C7404440'
           C(110 )='434EC032414EC032404BC242C341C240C33FC13EC0353C4EC23F
     +C13EC035324EC4403C32C7404440C7404440'
           C(111 )='464EBD3FBE3EBF3DC03EC13DC23EC33FC240C341C242C143C042
     +BF43BE42BD41BE40BE3FBE3EBF3DC03EC13DC23EC23F4240C241C242C143C042BF
     +43BE42BE414B32'
           C(112 )='4339C0554140C02B4052C242C241C240C33FC23EC13DC03EBF3D
     +BE3EBD3FBE40BE41BE42464BC23FC23EC13DC03EBF3DBE3EBE3F364EC4403C2BC7
     +404E47'
           C(113 )='4C39C0554140C02B3F52BE42BE41BE40BD3FBE3EBF3DC03EC13D
     +C23EC33FC240C241C2423A4BBE3FBE3EBF3DC03EC13DC23EC23F4339C7404447'
           C(114 )='434EC032414EC0324048C143C242C241C340C13FC03FBF3FBF41
     +C1413441C4403C32C7404A40'
           C(115 )='4A4CC142C03CBF42BF41BE41BC40BE3FBF3FC03EC13FC23FC53E
     +C23FC13F3547C13FC23FC53EC23FC13FC03DBF3FBE3FBC40BE41BF41BF42C03CC1
     +42503E'
           C(116 )='4355C02FC13DC23FC240C241C1423952C02FC13DC13F3A4EC840
     +4732'
           C(117 )='434EC035C13EC33FC240C341C242364BC035C13EC23F474EC032
     +414EC032314EC4404740C4403F32C4404440'
           C(118 )='424EC6323B4EC534464CBA32384EC6404440C6404232'
           C(119 )='434EC4323D4EC335444BBC32444EC4323D4EC335444BBC32314E
     +C7404940C6404232'
           C(120 )='424ECB32364ECB32404EB4323E4EC6404440C6403032C6404440
     +C6404440'
           C(121 )='424EC6323B4EC534464CBA32BE3CBE3EBE3FBF40BF41C141C13F
     +3E54C6404440C6404332'
           C(122 )='4B4EB5324C4EB532404EBF3CC044CC403432CC40C044BF3C4740
     +'
           C(123 )='4559BE3FBF3FBF3EC03EC13EC13FC13EC03EBE3E414EBF3EC03E
     +C13EC13FC13EC03EBF3EBC3EC43EC13EC03EBF3EBF3FBF3EC03EC13E3F4EC23EC0
     +3EBF3EBF3FBF3EC03EC13EC13FC23F4747'
           C(124 )='4059C0204747'
           C(125 )='4059C23FC13FC13EC03EBF3EBF3FBF3EC03EC23E3F4EC13EC03E
     +BF3EBF3FBF3EC03EC13EC43EBC3EBF3EC03EC13EC13FC13EC03EBF3E414EBE3EC0
     +3EC13EC13FC13EC03EBF3EBF3FBE3F4B47'
           C(126 )='4046C042C143C241C240C23FC43DC23FC240C241C1422E3EC142
     +C241C240C23FC43DC23FC240C241C143C0424834'
           C(127 )='5540'
      RETURN
      END

      SUBROUTINE XCSETB(C)
      CHARACTER*300 C(127)

           C(1) =' Small and Simple.'
           C(2) =' 6  8  '
           C(32 )='4640'
           C(33 )='C0414041C0444240403A'
           C(34 )='4044C1424140BF3E403C4340'
           C(35 )='4042C4404042BC404142C03A4240C046433A'
           C(36 )='4041C340C141BF41BE40BF41C141C3403E41C03A4440'
           C(37 )='4046C03FC140C041BF404440BC3A4440BF40C041C140C03F4240'
     +
           C(38 )='4442BE3EBF40BF41C041C242C041BF41BF3FC03FC43C4240'
           C(39 )='4044C142403A4240'
           C(40 )='40464240BE3EC03EC23E4240'
           C(41 )='4046C23EC03EBE3E4440'
           C(42 )='4241C044423EBC404442BC3C4044C43C423F'
           C(43 )='4241C044423EBC40463D'
           C(44 )='4041C03FBF3F41414240'
           C(45 )='4043C4404240403D'
           C(46 )='C041423F'
           C(47 )='C4464240403A'
           C(48 )='4140BF41C044C141C140C13FC03CBF3FBF404440'
           C(49 )='4045C141C03A3F40C2404240'
           C(50 )='4045C141C240C13FC03FBF3FBE40BF3FC03EC4404240'
           C(51 )='4045C141C240C13FC03FBF3FBF404140C13FC03FBF3FBE40BF414
     +03F4640'
           C(52 )='4442BC40C344C03A4340'
           C(53 )='4041C13FC240C141C042BF41BD40C042C4404240403A'
           C(54 )='4043C340C13FC03FBF3FBE40BF41C043C242C1404340403A'
           C(55 )='4046C440BD3A4540'
           C(56 )='4140BF41C041C141C240C141C041BF41BE40BF3FC03FC13F4240C
     +13FC03FBF3FBE404540'
           C(57 )='4140C140C242C043BF41BE40BF3FC03FC13FC3404240403D'
           C(58 )='4044C03F403FC03F423F'
           C(59 )='4044C03F403FC03EBF3F41414240'
           C(60 )='40464240BD3DC33D4240'
           C(61 )='4044C440403EBC404640403E'
           C(62 )='4046C33DBD3D4540'
           C(63 )='4045C141C140C13FC03FBF3FC03F403FC03F4340'
           C(64 )='4343BF3FBF40C041C141C140C03EC141C042BF41BE40BF3FC03CC
     +13FC3404240'
           C(65 )='C042C244C23CC03E3C42C440423E'
           C(66 )='C340C141C041BF41BE404240C141C041BF41BD404140C03A4540'
     +
           C(67 )='44404041BF3FBE40BF41C044C141C240C13F423E403D'
           C(68 )='C340C141C044BF41BD404140C03A4540'
           C(69 )='C046C4403C40403DC2403E40403DC4404240'
           C(70 )='C046C4403C40403DC240403D4440'
           C(71 )='4343C140C03DBD40BF41C044C141C3404240403A'
           C(72 )='C046403DC4404043C03A4240'
           C(73 )='4046C2403F40C03A3F40C2404240'
           C(74 )='4041C13FC240C141C0454240403A'
           C(75 )='C0464440BD3DBF404140C33D4240'
           C(76 )='4046C03AC4404240'
           C(77 )='C046C23CC244C03A4240'
           C(78 )='C046C43AC046403A4240'
           C(79 )='C046C440C03ABC404640'
           C(80 )='C046C340C13FC03FBF3FBD40463D'
           C(81 )='4242C13FBF3FBF40BF41C044C141C240C13FC03DBF3FC13F4240'
     +
           C(82 )='C046C340C13FC03FBF3FBD404140C33D4240'
           C(83 )='4041C13FC240C141BC44C141C240C13F4240403B'
           C(84 )='4046C4403E40C03A4440'
           C(85 )='4046C03BC13FC240C141C0454240403A'
           C(86 )='4046C33AC3464240403A'
           C(87 )='4046C23AC143C13DC2464240403A'
           C(88 )='C4463C40C43A4240'
           C(89 )='4046C23DC03D4043C2434240403A'
           C(90 )='4046C440BC3AC4404240'
           C(91 )='C046C240403ABE404440'
           C(92 )='4046C43A4240'
           C(93 )='4046C240C03ABE404440'
           C(94 )='4044C242C23E423C'
           C(95 )='403FC4404241'
           C(96 )='4046C13E423C'
           C(97 )='4240BF40BF41C042C141C140C13FC03EBF3F4141C13F4240'
           C(98 )='C046403CC242C140C13FC03EBF3FBF40BE42403E4640'
           C(99 )='4444BD40BF3FC03EC13FC3404240'
           C(100)='4442BE3EBF40BF41C042C141C140C23E4044C03A4240'
           C(101)='4042C340C141BF41BE40BF3FC03EC13FC2404340'
           C(102)='4043C3404142BF41BF40BF3FC03B4540'
           C(103)='403FC13FC240C141C044BF41BE40BF3FC03EC13FC3404240'
           C(104)='C046403CC242C140C13FC03D4240'
           C(105)='C0434041C0414240403B'
           C(106)='403FC13FC140C141C0444041C0414240403B'
           C(107)='C046403CC240C2423E3EC23E4240'
           C(108)='4046C03BC13F4240'
           C(109)='C044403FC141C13FC03F4041C141C13FC03D4240'
           C(110)='C044403FC141C140C13FC03D4240'
           C(111)='4340BE40BF41C042C141C240C13FC03EBF3F4340'
           C(112)='403EC046403FC141C240C13FC03EBF3FBD404640'
           C(113)='443EC046403FBF41BE40BF3FC03EC13FC3404240'
           C(114)='C044403EC242C140C13F4240403D'
           C(115)='C340C141BF41BE40BF41C141C340423C'
           C(116)='4044C4403E42C03BC13FC141423F'
           C(117)='4044C03DC13FC140C2424042C03C4240'
           C(118)='4044C23CC244423C'
           C(119)='4044C13CC144C13CC144423C'
           C(120)='C4443C40C43C4240'
           C(121)='4044C23C4244BD3ABF4040424640'
           C(122)='4044C440BC3CC4404240'
           C(123)='40464240BF3FC03FBF3FC13FC03FC13F4240'
           C(124)='C046403A4240'
           C(125)='C141C041C141BF41C041BF41403A4440'
           C(126)='4041C141C240C141423E403F'
           C(127)='4640'

      RETURN
      END

      SUBROUTINE XCSETC(C)
      CHARACTER*300 C(127)

           C(1)  ='  A Simple Large Character Set. '
           C(2)  ='30 40 '
           C(32 )='5340'
           C(33 )='4155BF3EC134C14CBF42403EC03A4035BF3FC13FC141BF414C3E'
     +
           C(34 )='4258BF3FBF41C141C13FC03EBF3EBF3F4845BF3FBF41C141C13FC
     +03EBF3EBF3F492D'
           C(35 )='4855B9244D5CB9243A51CE40313ACE40453C'
           C(36 )='4559C023445DC0234556BE42BD41BC40BD3FBE3EC03EC13EC13FC
     +23FC63EC23FC13FC13EC03DBE3EBD3FBC40BD41BE42543D'
           C(37 )='D2553340C23EC03EBF3EBE3FBE40BE42C042C142C241C240C23FC
     +33FC340C341C2413C32BE3FBF3EC03EC23EC240C241C142C042BE42BE404839'
           C(38 )='524DBF3FC13FC141C041BF41BF40BF3FBF3EBE3BBE3DBE3EBE3FB
     +D40BD41BF42C043C142C644C242C142C042BF42BE41BE3FBF3EC03EC13DC23DC53
     +9C23EC33FC140C141C041313EBE41BF42C043C142C2424046C13EC835C23EC23F4
     +640'
           C(39 )='4258BF3FBF41C141C13FC03EBF3EBF3F492D'
           C(40 )='4759BE3EBE3DBE3CBF3BC03CC13BC23CC23DC23E4747'
           C(41 )='4059C23EC23DC23CC13BC03CBF3BBE3CBE3DBE3E4E47'
           C(42 )='454FC0343B49CA3A4046B63A503A'
           C(43 )='4952C02E3749D2404837'
           C(44 )='4241BF3FBF41C141C13FC03EBF3EBF3F4A44'
           C(45 )='4049D2404837'
           C(46 )='4142BF3FC13FC141483F'
           C(47 )='5259AE205647'
           C(48 )='4655BD3FBE3DBF3BC03DC13BC23DC33FC240C341C243C145C043B
     +F45BE43BD41BE404E2B'
           C(49 )='4051C241C343C02B4B40'
           C(50 )='4150C041C142C141C241C440C23FC13FC13EC03EBF3EBE3DB636C
     +E404640'
           C(51 )='4255CB40BA38C340C23FC13FC13DC03EBF3DBE3EBD3FBD40BD41B
     +F41BF42543C'
           C(52 )='4A55B632CF403B4EC02B4A40'
           C(53 )='4C55B640BF37C141C341C340C33FC23EC13DC03EBF3DBE3EBD3FB
     +D40BD41BF41BF42543C'
           C(54 )='4C52BF42BD41BE40BD3FBE3DBF3BC03BC13CC23EC33FC140C341C
     +242C143C041BF43BE42BD41BF40BD3FBE3EBF3D5439'
           C(55 )='4055CE40B62B5040'
           C(56 )='4555BD3FBF3EC03EC13EC23FC43FC33FC23EC13EC03DBF3EBF3FB
     +D3FBC40BD41BF41BF42C043C142C242C341C441C241C142C042BF42BD41BC404F2
     +B'
           C(57 )='4D4EBF3DBE3EBD3FBF40BD41BE42BF43C041C143C242C341C140C
     +33FC23EC13CC03BBF3BBE3DBD3FBE40BD41BF42533D'
           C(58 )='414EBF3FC13FC141BF414034BF3FC13FC141BF414E3E'
           C(59 )='414EBF3FC13FC141BF414032BF41C141C13FC03EBF3EBF3F4E44'
     +
           C(60 )='5052B037D0374440'
           C(61 )='404CD2402E3AD240483A'
           C(62 )='4052D037B0375440'
           C(63 )='4151C13FBF3FBF41C041C142C141C241C340C33FC13EC03EBF3EB
     +F3FBC3EC03D414E423F413F413E403E3F3E3E3E3E37BF3FC13FC141BF414D3E'
           C(64 )='4F4DBF42BE41BD40BE3FBF3FBF3DC03DC13EC23FC340C241C1423
     +B48BE3EBF3DC03DC13EC13F474BBF38C03EC23FC240C242C143C042BF43BF42BE4
     +2BE41BD41BD40BD3FBE3FBE3EBF3EBF3DC03DC13DC13EC23EC23FC33FC340C341C
     +241C1413E4DBF38C03EC13F483B'
           C(65 )='C855C82B3347CA404939'
           C(66 )='C055C940C33FC13FC13EC03EBF3EBF3FBD3F3740C940C33FC13FC
     +13EC03DBF3EBF3FBD3FB7405540'
           C(67 )='4F50BF42BE42BE41BC40BE3FBE3EBF3EBF3DC03BC13DC13EC23EC
     +23FC440C241C242C142463B'
           C(68 )='C055C740C33FC23EC13EC13DC03BBF3DBF3EBE3EBD3FB9405540'
     +
           C(69 )='C055CD403336C8403835CD404640'
           C(70 )='C055CD403336C8404A35'
           C(71 )='4F50BF42BE42BE41BC40BE3FBE3EBF3EBF3DC03BC13DC13EC23EC
     +23FC440C241C242C142C043BB404B38'
           C(72 )='C0554E40C02B324BCE404835'
           C(73 )='C055482B'
           C(74 )='4A55C030BF3DBF3FBE3FBE40BE41BF41BF43C0425039'
           C(75 )='C0554E40B2324545C9344740'
           C(76 )='4055C02BCC404540'
           C(77 )='C055C82BC855C02B4840'
           C(78 )='C055CE2BC055482B'
           C(79 )='4655BE3FBE3EBF3EBF3DC03BC13DC13EC23EC23FC440C241C242C
     +142C143C045BF43BF42BE42BE41BC40502B'
           C(80 )='C055C940C33FC13FC13EC03DBF3EBF3FBD3FB7405536'
           C(81 )='4655BE3FBE3EBF3EBF3DC03BC13DC13EC23EC23FC440C241C242C
     +142C143C045BF43BF42BE42BE41BC40432FC63A4742'
           C(82 )='C055C940C33FC13FC13EC03EBF3EBF3FBD3FB7404740C7354740'
     +
           C(83 )='4E52BE42BD41BC40BD3FBE3EC03EC13EC13FC23FC63EC23FC13FC
     +13EC03DBE3EBD3FBC40BD41BE42543D'
           C(84 )='4755C02B3955CE40442B'
           C(85 )='4055C031C13DC23EC33FC240C341C242C143C04F482B'
           C(86 )='4055C82BC855442B'
           C(87 )='4055C52BC555C52BC555442B'
           C(88 )='CE553240CE2B4640'
           C(89 )='4055C836C84A3836C0354C40'
           C(90 )='4055CE40B22BCE404640'
           C(91 )='4759B940C020C7404847'
           C(92 )='4059D2204447'
           C(93 )='4059C740C020B9404C47'
           C(94 )='4052C747C739462E'
           C(95 )='403ED4404442'
           C(96 )='4058C13FC141BF41BF3FC03EC13EC13F472D'
           C(97 )='4C4EC032404BBE42BE41BD40BE3FBE3EBF3DC03EC13DC23EC23FC
     +340C241C242473D'
           C(98 )='C0554036C242C241C340C23FC23EC13DC03EBF3DBE3EBE3FBD40B
     +E41BE42533D'
           C(99 )='4C4BBE42BE41BD40BE3FBE3EBF3DC03EC13DC23EC23FC340C241C
     +242463D'
           C(100)='4C55C02B404BBE42BE41BD40BE3FBE3EBF3DC03EC13DC23EC23FC
     +340C241C242473D'
           C(101)='4048CC40C042BF42BF41BE41BD40BE3FBE3EBF3DC03EC13DC23EC
     +23FC340C241C242463D'
           C(102)='4855BE40BE3FBF3DC02F3D4EC7404532'
           C(103)='4C4EC030BF3DBF3FBE3FBD40BE414951BE42BE41BD40BE3FBE3EB
     +F3DC03EC13DC23EC23FC340C241C242473D'
           C(104)='C0554035C343C241C340C23FC13DC0364840'
           C(105)='4054C13FC141BF41BF3F413AC0324740'
           C(106)='4454C13FC141BF41BF3F413AC02FBF3DBE3FBE404A47'
           C(107)='C0554A39B6364444C7384640'
           C(108)='C055482B'
           C(109)='C04E403CC343C241C340C23FC13DC036404AC343C241C340C23FC
     +13DC0364840'
           C(110)='C04E403CC343C241C340C23FC13DC0364840'
           C(111)='454EBE3FBE3EBF3DC03EC13DC23EC23FC340C241C242C143C042B
     +F43BE42BE41BD404E32'
           C(112)='404EC02B4052C242C241C340C23FC23EC13DC03EBF3DBE3EBE3FB
     +D40BE41BE42533D'
           C(113)='4C4EC02B4052BE42BE41BD40BE3FBE3EBF3DC03EC13DC23EC23FC
     +340C241C242473D'
           C(114)='C04E403AC143C242C241C3404532'
           C(115)='4B4BBF42BD41BD40BD3FBF3EC13EC23FC53FC23FC13EC03FBF3EB
     +D3FBD40BD41BF42513D'
           C(116)='4355C02FC13DC23FC240384EC7404732'
           C(117)='404EC036C13DC23FC340C241C343404AC0324840'
           C(118)='404EC632C64E4432'
           C(119)='404EC432C44EC432C44E4632'
           C(120)='CB4E4032B54E5132'
           C(121)='414EC632464EBA32BE3CBE3EBE3FBF405047'
           C(122)='404ECB40B532CB404640'
           C(123)='4559BE3FBF3FBF3EC03EC13EC13FC13EC03EBE3EBE3FC23FC23EC
     +03EBF3EBF3FBF3EC03EC13EC13FC23F4747'
           C(124)='4059C0204747'
           C(125)='4059C23FC13FC13EC03EBF3EBF3FBF3EC03EC23EC23FBE3FBE3EC
     +03EC13EC13FC13EC03EBF3EBF3FBE3F4B47'
           C(126)='4046C042C143C241C240C23FC43DC23FC240C241C1422E3EC142C
     +241C240C23FC43DC23FC240C241C143C0424834'
           C(127)='5340'

      RETURN
      END

      SUBROUTINE XCSETD(C)
      CHARACTER*300 C(127)

           C(1)  =' Italic Character set.'
           C(2)  ='20 30 '
           C(32 )='5540'
           C(33 )='4655BF3FBE34434CBD34434DC13FBC343E3ABF3FC13FC141BF414
     +E3E'
           C(34 )='4857BD40C242C13EBF3EBE3E4B44BD40C242C13EBF3EBE3E4A2D'
     +
           C(35 )='4855B9244D5CB9243A51CE40313ACE40453C'
           C(36 )='4A59B8234D5DB8234955BF3FC13FC141C041BF42BF41BD41BC40B
     +D3FBE3EC03EC13EC13FC73CC23E3549C23EC73CC13FC13EC03DBF3EBF3FBD3FBC4
     +0BD41BF41BF42C041553B'
           C(37 )='D2553340C23EC03EBF3EBE3FBE40BE42C042C142C241C240C23FC
     +33FC340C341C2413C32BE3FBF3EC03EC23EC240C241C142C042BE42BE404839'
           C(38 )='524DBF3FC13FC141C041BF41BF40BF3FBF3EBE3BBE3DBE3EBE3FB
     +D40BD41BF42C043C142C644C242C142C042BF42BE41BE3FBF3EC03EC13DC23DC53
     +9C23EC33FC140C141C041313EBE41BF42C043C142C2424046C13EC835C23EC23F4
     +640'
           C(39 )='4857BD40C242C13EBF3EBE3E4A2D'
           C(40 )='4C59BC3DBD3DBE3DBE3CBF3BC03CC13BC13DC13E455DBD3CBE3CB
     +F3DBF3BC03B4E3F'
           C(41 )='4959C13EC13DC13BC03CBF3BBE3CBE3DBD3DBC3D4960C13DC13BC
     +03BBF3BBF3D473C'
           C(42 )='4655C0343B49CA3A3640CA46462E'
           C(43 )='4D52BB2E3A49D2404637'
           C(44 )='4340BD40C242C13EBF3EBE3E4C44'
           C(45 )='4049D2404837'
           C(46 )='4140BF41C141C13F493F'
           C(47 )='5A59A6205847'
           C(48 )='4955BD3FBE3EBE3DBF3DBF3CC03DC13DC13FC23FC240C341C242C
     +243C143C144C043BF43BF41BE41BE40BE3FBE3EBE3DBF3DBF3CC03DC13DC23E424
     +0C241C242C243C143C144C043BF43482D'
           C(49 )='4140C8553E3CBB2F4140C655BD3DBD3EBE3F4743BC3E4A30'
           C(50 )='454F4142C13FBF3FBF41C041C142C141C341C340C33FC13EC03EB
     +F3EBE3EBD3EBC3EBD3EBE3EBE3C4D55C23FC13EC03EBF3EBE3EBA3C3A3AC141C24
     +0C53EC340C241C142353FC53DC340C241C143463C'
           C(51 )='4551C13FBF3FBF41C041C142C141C341C340C33FC13EC03EBF3EB
     +D3EBD3F434AC23FC13EC03EBF3EBE3E3B3FC240C33FC13FC13EC03DBF3EBF3FBD3
     +FBC40BD41BF41BF42C041C141C13FBF3F4847C23FC13FC13EC03DBF3EBF3FBE3F4
     +C40'
           C(52 )='4E54BA2C4755BA2B4655B131CC40493A'
           C(53 )='5046374FBB36454ACA40363FC540C5413136C141C341C340C33FC
     +13FC13EC03DBF3DBE3EBD3FBD40BD41BF41BF42C041C141C13FBF3F4849C23FC13
     +FC13EC03DBF3DBE3EBE3F4D40'
           C(54 )='4E52BF3FC13FC141C041BF42BE41BD40BD3FBE3EBE3DBF3DBF3CC
     +03CC13EC13FC23FC340C341C242C142C043BF42BF41BE41BD40BE3FBE3EBF3E484
     +EBE3FBE3EBE3DBF3DBF3CC03BC13E453FC241C242C142C0444937'
           C(55 )='4355BE3A4F46BF3DBE3DBB3ABE3DBF3EBF3C494FBA3ABE3DBF3EB
     +F3C3F52C343C240C53D3741C241C240C53EC240C141452D'
           C(56)='51553840BD3FBF3FBF3EC03DC13EC23FC340C441C141C142C043BF
     +42BD41BD40BE3FBF3FBF3EC03DC13EC13F4340C341C141C142C043BF42BE413B36
     +BC3FBE3EBF3EC03DC13EC33FC440C441C141C142C043BF42BF41BE413D40BD3FBE
     +3EBF3EC03DC13EC23F4440C341C141C142C0444838'
           C(57 )='4B4A4344BF3EBE3EBE3FBD40BE41BF41BF42C043C142C242C341C
     +340C23FC13FC13EC03CBF3CBF3DBE3DBE3EBD3FBD40BE41BF42C041C141C13FBF3
     +F4347BF42C044C142C242C241453FC13EC03BBF3CBF3DBE3DBE3E4D3F'
           C(58 )='444EBF3FC13FC141BF413D34BF3FC13FC141BF414D3F'
           C(59 )='454EBF3FC13FC141BF413D32BF41C141C13FC03FBF3EBE3E4D44'
     +
           C(60 )='5352AD37D0374640'
           C(61 )='434CD2402B3AD240483A'
           C(62 )='4252CE37AE375540'
           C(63 )='4151C13FBF3FBF41C041C142C141C341C440C33FC13EC03EBF3EB
     +F3FBA3EBE3FC03EC13FC240434EC23FC13EC03EBF3EBF3FBE3F3A36BF3FC13FC14
     +1BF41553E'
           C(64 )='4F4DBF42BE41BD40BE3FBF3FBF3DC03DC13EC23FC340C241C1423
     +B48BE3EBF3DC03DC13EC13F474BBF38C03EC23FC240C242C143C042BF43BF42BE4
     +2BE41BD41BD40BD3FBE3FBE3EBF3EBF3DC03DC13DC13EC23EC23FC33FC340C341C
     +241C1413E4DBF38C03EC13F483B'
           C(65 )='4F55B32B4D55C12B3E53C12D3746C940313AC6404640C6404540'
     +
           C(66 )='4955BA2B4755BA2B4255CB40C33FC13EC03EBF3DBF3FBD3F414AC
     +23FC13EC03EBF3DBF3FBE3F3740C940C23FC13EC03EBF3DBE3EBC3FB440504BC13
     +FC13EC03EBF3DBE3EBD3F4C40'
           C(67 )='4F53C140C142BF3AC042BF42BF41BE41BD40BD3FBE3EBE3DBF3DB
     +F3CC03DC13DC13FC33FC340C241C242C1423C50BE3FBE3EBE3DBF3DBF3CC03DC13
     +DC13FC23F4E40'
           C(68 )='4955BA2B4755BA2B4255C940C33FC13FC13DC03CBF3CBE3CBE3EB
     +E3FBC3FB7404F55C23FC13FC13DC03CBF3CBE3CBE3EBE3FBD3F4E40'
           C(69 )='4955BA2B4755BA2B4A4FBE383A4ECF40BF3AC0463336C6403335C
     +F40C245BD3B4940'
           C(70 )='4955BA2B4755BA2B4A4FBE383A4ECF40BF3AC0463336C6403335C
     +7405040'
           C(71 )='4F53C140C142BF3AC042BF42BF41BE41BD40BD3FBE3EBE3DBF3DB
     +F3CC03DC13DC13FC33FC240C341C242C2443B4EBE3FBE3EBE3DBF3DBF3CC03DC13
     +DC13FC23F4240C241C242C2443D40C7404439'
           C(72 )='4955BA2B4755BA2B5255BA2B4755BA2B3555C7404640C7402D36C
     +C402D35C7404640C7404640'
           C(73 )='4955BA2B4755BA2B4255C740332BC7404840'
           C(74 )='4E55BB2FBF3EBF3FBE3FBE40BE41BF42C042C141C13FBF3F4C51B
     +B2FBF3EBE3E4555C740412B'
           C(75 )='4955BA2B4755BA2B5355AF334744C4343B4CC4343655C7404640C
     +640272BC7404640C6404440'
           C(76 )='4955BA2B4755BA2B4255C740332BCF40C246BD3A4640'
           C(77 )='4955BA2B4655C12B4055C12D4C53B32B4D55BA2B4755BA2B3455C
     +4404D40C440252BC6404840C7404640'
           C(78 )='4955BA2B4655C72E394FC72E4655BA2B3655C3404A40C640272BC
     +6405340'
           C(79 )='4A55BD3FBE3EBE3DBF3DBF3CC03DC13DC13FC23FC340C341C242C
     +243C143C144C043BF43BF41BE41BD40BE3FBE3EBE3DBF3DBF3CC03DC13DC23E434
     +0C241C242C243C143C144C043BF43BE42482B'
           C(80 )='4955BA2B4755BA2B4255CC40C33FC13EC03EBF3DBE3EBC3FB8404
     +B4BC23FC13EC03EBF3DBE3EBD3F3136C7405040'
          C(81)='4A55BD3FBE3EBE3DBF3DBF3CC03DC13DC13FC23FC340C341C242C24
     +3C143C144C043BF43BF41BE41BD40BE3FBE3EBE3DBF3DBF3CC03DC13DC23E4340C
     +241C242C243C143C144C043BF43BE42362DC041C142C241C140C23FC13EC039C13
     +FC240C142C0413C45C13AC13FC140C1414743'
           C(82 )='4955BA2B4755BA2B4255CB40C33FC13EC03EBF3DBF3FBD3FB7404
     +A4AC23FC13EC03EBF3DBF3FBE3F3C40C23FC13FC138C13FC240C142C0413B46C23
     +9C13FC140C1412C3EC7405140'
           C(83 )='5153C140C142BF3AC042BF42BF41BD41BC40BD3FBE3EC03EC13EC
     +13FC73CC23E3549C23EC73CC13FC13EC03DBF3EBF3FBD3FBC40BD41BF41BF42C04
     +2BF3AC142C140553E'
           C(84 )='4B55BA2B4755BA2B3F55BD3AC246CF40BF3AC046302BC7404C40'
     +
           C(85 )='4555BD35BF3CC03DC13EC33FC440C341C242C143C44F3340BD35B
     +F3CC03DC13EC23F3D55C7404740C640442B'
           C(86 )='4555C12B4055C12D4C53B32B3D55C6404640C640432B'
           C(87 )='4555BE2B4355BE2D4953B62B4A55BE2B4355BE2D4953B62B3755C
     +7404940C640422B'
           C(88 )='4855C72B3A55C72B4655AC2B4455C6404640C640282BC6404640C
     +6404740'
           C(89 )='4555C436BD354055C436BD354D55B636394AC6404740C6402D2BC
     +7404C40'
           C(90 )='5355AD2B5455AD2B4655BD3AC246CE402C2BCE40C246BD3A4940'
     +
           C(91 )='4859B8204960B8204760C7403120C7404C47'
           C(92 )='4959CC204447'
           C(93 )='4E59B8204960B8204160C7403120C7404C47'
           C(94 )='4552C747C739452E'
           C(95 )='403ED4404442'
           C(96 )='4058C13FC141BF41BF3FC03EC13EC13F472D'
           C(97 )='4D4EBE39BF3CC03EC13FC340C242C1423D4ABE39BF3CC03EC13F3
     +F47C043BF43BE41BE40BD3FBE3DBF3DC03DC13EC13FC23FC240C241C243C1433B4
     +7BE3FBE3DBF3DC03CC13E533F'
           C(98 )='4455BC33C03DC13DC13F4354BC33C143C242C241C240C23FC13FC
     +13EC03DBF3DBE3DBD3FBE40BE41BF43C0444945C13EC03CBF3DBE3DBE3F3B55C44
     +04E2B'
           C(99 )='4B4BC03FC140C041BF42BE41BD40BD3FBE3DBF3DC03DC13EC13FC
     +23FC240C341C2433B4ABE3FBE3DBF3DC03CC13E4F3F'
           C(100)='4F55BC32BF3CC03EC13FC340C242C1423F51BC32BF3CC03EC13F3
     +F47C043BF43BE41BE40BD3FBE3DBF3DC03DC13EC13FC23FC240C241C243C1433B4
     +7BE3FBE3DBF3DC03CC13E4A54C440452B'
           C(101)='4145C441C341C342C142BF42BE41BD40BD3FBE3DBF3DC03DC13EC
     +13FC23FC240C341C2423B4BBE3FBE3DBF3DC03CC13E4F3F'
           C(102)='5054BF3FC13FC141C041BF41BE40BE3FBF3FBF3EBF3DBD32BF3CB
     +F3E4A5BBE3EBF3EBF3CBE37BF3CBF3DBF3EBF3FBE3FBE40BF41C041C141C13FBF3
     +F4554CA404032'
           C(103)='504EBC32BF3DBE3DBD3FBD40BE41BF41C041C141C13FBF3F4E53B
     +C32BF3DBE3DBE3F474EC043BF43BE41BE40BD3FBE3DBF3DC03DC13EC13FC23FC24
     +0C241C243C1433B47BE3FBE3DBF3DC03CC13E503F'
           C(104)='4655BA2B4755BA2B4247C244C242C241C240C23FC13FC03EBE3AC
     +03DC13F3E4EC23EC03EBE3AC03DC13FC340C242C1423151C4404E2B'
           C(105)='4955BF3FC13FC141BF413835C142C242C340C13FC03DBE3AC03DC
     +13F3F4EC13FC03DBE3AC03DC13FC340C242C142443C'
           C(106)='4B55BF3FC13FC141BF413835C142C242C340C13FC03DBD36BF3DB
     +F3EBF3FBE3FBE40BF41C041C141C13FBF3F4854C13FC03DBD36BF3DBF3EBE3E4B4
     +7'
           C(107)='4755BA2B4755BA2B4D4DBF3FC13FC141C041BF41BF40BE3FBC3CB
     +E3FBE404240C23FC23AC13F3B48C13FC23AC13FC240C241C2433451C4404D2B'
           C(108)='4555BC32BF3CC03EC13FC340C242C1423F51BC32BF3CC03EC13F4
     +055C440462B'
           C(109)='404AC142C242C340C13FC03EBF3CBE39414EC13FC03EBF3CBE394
     +347C244C242C241C240C23FC13FC03EBD36404EC23EC03EBD364347C244C242C24
     +1C240C23FC13FC03EBE3AC03DC13F3E4EC23EC03EBE3AC03DC13FC340C242C1424
     +33C'
           C(110)='404AC142C242C340C13FC03EBF3CBE39414EC13FC03EBF3CBE394
     +347C244C242C241C240C23FC13FC03EBE3AC03DC13F3E4EC23EC03EBE3AC03DC13
     +FC340C242C142433C'
           C(111)='464EBD3FBE3DBF3DC03DC13EC13FC23FC240C341C243C143C043B
     +F42BF41BE41BE40BE3FBE3DBF3DC03CC13E443FC241C243C143C044BF424733'
           C(112)='424AC142C242C340C13FC03EBF3CBC324355C13FC03EBF3CBC324
     +54EC143C243C241C240C23FC13FC13EC03DBF3DBE3DBD3FBE40BE41BF43C043494
     +6C13EC03CBF3DBE3DBE3F3339C7404F47'
           C(113)='4D4EBA2B4755BA2B434EC043BF43BE41BE40BD3FBE3DBF3DC03DC
     +13EC13FC23FC240C241C243C1433B47BE3FBE3DBF3DC03CC13E4238C7404847'
           C(114)='404AC142C242C340C13FC03EBF3CBE39414EC13FC03EBF3CBE394
     +347C244C242C241C240C13FC03FBF3FBF41C1414333'
           C(115)='4C4CC03FC140C041BF41BD41BD40BD3FBF3FC03EC13FC73CC13F3
     +747C13FC73CC13FC03DBF3FBD3FBD40BD41BF41C041C140C03F503E'
           C(116)='4655BC32BF3CC03EC13FC340C242C1423F51BC32BF3CC03EC13F3
     +D4EC9404532'
           C(117)='404AC142C242C340C13FC03DBE3AC03EC23E3E4EC13FC03DBE3AC
     +03EC13FC23FC240C241C242C2444247BE39BF3CC03EC13FC340C242C1423D4ABE3
     +9BF3CC03EC13F4740'
           C(118)='404AC142C242C340C13FC03DBE3AC03EC23E3E4EC13FC03DBE3AC
     +03EC13FC23FC140C341C242C243C144C044BF40C13E4434'
           C(119)='404AC142C242C340C13FC03DBE3AC03EC23E3E4EC13FC03DBE3AC
     +03EC13FC23FC240C241C242C1424249BE37C03DC13FC23FC240C241C242C142C14
     +4C045BF40C13E3842BE37C03DC23E4C40'
           C(120)='414AC243C241C340C13EC03D3E45C13EC03DBF3CBF3EBE3EBE3FB
     +F40BF41C041C141C13FBF3F4644C03DC13EC340C241C2434049BF3FC13FC141C04
     +1BF41BF40BE3FBE3EBF3EBF3CC03DC13E4B40'
           C(121)='404AC142C242C340C13FC03DBE3AC03EC23E3E4EC13FC03DBE3AC
     +03EC13FC23FC240C241C242C2444347BC32BF3DBE3DBD3FBD40BE41BF41C041C14
     +1C13FBF3F4E53BC32BF3DBE3DBE3F4D47'
           C(122)='4E4EBF3EBE3EB83ABE3EBF3E414AC142C242C340C43E3740C241C
     +340C43FC2403436C240C43FC340C2413740C43EC340C242C142473C'
           C(123)='4C59BD3FBF3FBF3EC03DC13DC03FBF3EBE3E444EBF3DC03EC13DC
     +03EBF3EBF3FBB3EC43EC13FC03EBF3DBD3EBF3FBF3EC03E434EC13FC13EC03FBF3
     +EBE3FBF3FBE3DC03FC13EC23F40474940'
           C(124)='4859B8204C47'
           C(125)='4659C23FC13FC03EBF3DBD3DBF3EC03FC23E414EC13EBF3DBD3DB
     +F3EC03EC13FC43EBB3EBE3DC03EC13DC03D444CBD3EBF3EC03FC13EC03DBF3EBE3
     +E4B47'
           C(126)='4046C042C143C241C240C23FC43DC23FC240C241C1422E3EC142C
     +241C240C23FC43DC23FC240C241C143C0424834'
           C(127)='5540'

      RETURN
      END

      SUBROUTINE XCPALET(mode)
c
c#######################################################################
c
c     PURPOSE:
c
c     Generate color label plots of 2-d field A given its
c     coordinate using ZXPLOT and ncar package..
c
c#######################################################################
c
c     AUTHOR: Min Zou
c     15/08/92
C
c#######################################################################
c
c     INPUT:
c
c     ctrlvls(nctrlvls)     Contour values dividing the filled areas
c     clrindx(nctrlvls-1)   Plot color index bar color index
c     nctrlvls              Number of contour levels
c
c     mode  Option for positioning the color palette
c           = 1, color bar is located below the plotting space
c           = 2, color bar is located to the right of plotting space
c
c#######################################################################
c
c     Variable Declarations.
c
c#######################################################################
c
      implicit none
c
      integer mode

      integer nctrlvls_max
      parameter(nctrlvls_max=1000)    ! Max. number of contour values
      real ctrlvls(nctrlvls_max)     ! contour values dividing the filled areas
      integer clrindx(nctrlvls_max)  ! plot color index bar color index
      integer nctrlvls, nctrlvls_lim ! Number of contour levels
      common /xcflvls/nctrlvls,ctrlvls,clrindx

      character cpalnfmt*15, xtem*15
      common /xcplnfmt/ cpalnfmt

      integer icplswitch
      common /xcplswitch/ icplswitch

      real xl,xr,yb,yt
      character*20 ch
      integer lch
      real xra(5),yra(5)    ! array for single color box

      real dtx,dty,x,y,xs,ys
      real byt,byb,bxl,bxr

      integer k,kwndon,iskip
      real xwd1,xwd2,ywd1,ywd2,hch
c
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
C     Beginning of executable code...
C
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
c     print*,'cpalnfmt=',cpalnfmt

c     print*,'mode, icplswitch=', mode, icplswitch

      IF( mode .eq. 2 .and. icplswitch .eq. 0 ) return

      IF(nctrlvls .gt.nctrlvls_max) THEN
        write(6,'(a,/a,i5,a)')
     :   'The number of contours exceeded maximum allowed.',
     :   'Only ',nctrlvls_max,' contours will be plotted'
      ENDIF
      nctrlvls_lim =min(nctrlvls,nctrlvls_max)

      call xqmap(xl,xr,yb,yt)

      byt=yb-0.13*(yt-yb)
      byb=byt-0.06*(yt-yb)
c
c Find out old window setting
c
      call xqwdwon(kwndon)
      if( kwndon.eq.1) then
        call xqwndw(xwd1,xwd2,ywd1,ywd2)
        call xwdwof
      endif

      IF(nctrlvls_lim.lt.15)then
        iskip=1
      ELSEIF(nctrlvls_lim.lt.30) then
        iskip=2
      else
        iskip=3
      endif

      call xqchsz(hch)

      if(mode.eq.1) then   ! Place color bar below plotting window.

        byt= yb-0.07*(yt-yb)
        byb=byt-0.03*(yt-yb)

        xs= xr-xl
        ys= byt-byb
        dtx=xs/(nctrlvls_lim-1)

        DO k=1,nctrlvls_lim-1

          xra(1)=xl+dtx*(k-1)
          xra(2)=xl+dtx*k
          xra(3)=xra(2)
          xra(4)=xra(1)
          yra(1)=byt
          yra(2)=byt
          yra(3)=byb
          yra(4)=byb
          CALL XCOLOR(clrindx(k))
          CALL XFILAREA(xra,yra,4)

          CALL xcolor(1)
          call xbox(xra(1),xra(2),yra(3),yra(1))

          IF(mod(k-1,iskip).eq.0) then
            if( cpalnfmt(1:1).eq.'*') then
              CALL XRCH_new(ctrlvls(k),ch,lch)
            else
              if( index(cpalnfmt,'I').eq.0 .and.
     :            index(cpalnfmt,'i').eq.0 ) then
                write(ch,cpalnfmt) ctrlvls(k)
              else
                write(ch,cpalnfmt) nint(ctrlvls(k))
              endif
              lch = 15
              call xstrlnth(ch,lch)
            endif
            CALL xcharc(xra(1),byb-1.3*hch,ch(1:lch))
          END IF

        END DO

        CALL xcolor(1)

        CALL XRCH_new(ctrlvls(nctrlvls_lim),ch,lch)
        CALL xcharc(xr,byb-1.3*hch,ch(1:lch))

      else if( mode == 2 ) then ! Place color bar to the right of plotting window.

        bxr = xr+0.07*(xr-xl)
        bxl = xr+0.03*(xr-xl)

        xs = bxr-bxl
        ys = 0.94*(yt-yb)
        dty=ys/(nctrlvls_lim-1)
        x=bxr+0.20*xs

        DO k=1,nctrlvls_lim-1

          yra(1)=0.030*(yt-yb)+yb+dty*(k-1)
          yra(2)=0.030*(yt-yb)+yb+dty*k
          yra(3)=yra(2)
          yra(4)=yra(1)
          xra(1)=bxl
          xra(2)=bxl
          xra(3)=bxr
          xra(4)=bxr

          CALL XCOLOR(clrindx(k))
          CALL XFILAREA(xra,yra,4)
          CALL xcolor(1)
          call xbox(xra(1),xra(3),yra(1),yra(2))

          IF(mod(k-1,iskip).eq.0) then
            if( cpalnfmt(1:1).eq.'*') then
              CALL XRCH_new(ctrlvls(k),ch,lch)
            else

              if( index(cpalnfmt,'I').eq.0 .and.
     :            index(cpalnfmt,'i').eq.0 ) then
                write(ch,cpalnfmt) ctrlvls(k)
              else
                write(ch,cpalnfmt) nint(ctrlvls(k))
              endif
              lch = 15
              call xstrlnth(ch,lch)
            endif
            CALL xcharl(x,yra(1)-0.3*hch,ch(1:lch))
          END IF
        END DO

        CALL xcolor(1)
        y=0.025*(yt-yb)+yb+ys
        CALL XRCH_new(ctrlvls(nctrlvls_lim),ch,lch)
        CALL xcharl(x,y-0.3*hch,ch(1:lch))

      end if
c
c Restore old windin setting
c
      if( kwndon.eq.1) call xwindw(xwd1,xwd2,ywd1,ywd2)

      RETURN

      ENTRY XCPALNFMT( xtem )

      cpalnfmt = xtem

      RETURN
      END SUBROUTINE XCPALET

      SUBROUTINE XRCH_new( R,CH,LCH)

C Return real number R as a character string in automatically set format
      REAL R
      CHARACTER CH*20, STR*20

      CALL XGETFMT(R,STR)
      IF(ABS(R).LT. 1.E-20) THEN
        WRITE(CH,'(F3.1)') R
      ELSE IF(ABS(R).GE. 1.E-2 .AND. ABS(R).LT. 1.0) THEN
        WRITE(CH,'(F4.2)') R
      ELSE
        WRITE(CH,STR) R
      ENDIF

      LCH=20
      CALL xstrlnth( CH, LCH)
      CALL xstrmin ( CH, LCH)
      RETURN
      END SUBROUTINE XRCH_new

      SUBROUTINE xgetfmt(R,CH)
      INTEGER NPOZ
      CHARACTER CH*20,FORM,NDROB
      WRITE(CH,10)R
 10   FORMAT(G11.4)
      DO I=20,1,-1
        IF(CH(I:I).EQ.'0'.OR.CH(I:I).EQ.' ') THEN
           CH(I:I)=' '
        ELSE
          GOTO 1
        END IF
      END DO
 1    CONTINUE
      NPOZ=0
      NDOT=0
      NMANT=0
      NDROB=' '
      FORM='F'
      DO I = 1,20
        IF(CH(I:I).NE.' ' ) NPOZ=NPOZ+1
        IF(CH(I:I).EQ.'E') FORM='E'
        IF(NDROB.EQ.'.'.AND.CH(I:I).NE.' ') NDOT=NDOT+1
        IF(CH(I:I).EQ.'.') NDROB='.'
        IF(FORM.NE.'E') NMANT=NPOZ
      END DO
      NPOZ=NPOZ
      IF(FORM.EQ.'F') THEN
        IF(NDOT.NE.0) THEN
          write(CH,20) '(',FORM,NPOZ,'.',NDOT,')'
        ELSE
          write(CH,20) '(',FORM,NPOZ,'.',NDOT,')'
        END IF
      elseif(FORM.EQ.'E') then
        CH = '(1PE20.2)'
      ELSE
        write(CH,20) '(',FORM,NPOZ,'.',NMANT,')'
      END IF
 20   FORMAT(A1,A1,I1,A1,I1,A1)
      RETURN
      END SUBROUTINE xgetfmt

      SUBROUTINE XSTRMIN( string, length )
c
c#######################################################################
c
c     PURPOSE:
c
c     Minimize a string length by removing consecutive blank spaces.
c
c#######################################################################
c
c     AUTHOR: Ming Xue
c     1/15/93
c
c#######################################################################
c
c     INPUT:
c       string   A character string
c       length   The declared length of the character string 'string'.
c     OUTPUT:
c       length   The length of string with consecutive blank spaces
c                removed.
c
c#######################################################################

      implicit none
      character string*(*)
      integer length
      character str*256, str_1
      integer i,len_old
c
      IF( length.gt.256) THEN
        print*,'Work string defined in XSTRMIN was too small.'
        print*,'The output from this subroutine may not be correct.'
        length=256
      ENDIF

      len_old = length
      length = 1

      str = string
      DO 100 i = 2,len_old

        str_1 = str(i-1:i-1)
        IF(.not.(str(i:i).eq.' '.and.
     :    (str_1.eq.' '.or.str_1.eq.'('.or.str_1.eq.'='))) THEN
          length=length+1
          string(length:length)=str(i:i)
        ENDIF

100   CONTINUE

      DO 200 i = 1,length
        if( string(i:i).ne.' ') goto 300
200   CONTINUE
300   CONTINUE
      IF( i.ne.1) then
        str=string
        string(1:length-i+1)=str(i:length)
        length = length-i+1
      ENDIF

      RETURN
      END

      SUBROUTINE XSTPJGRD(mapproj,trulat1,trulat2,trulon,
     :           ctrlat,ctrlon,xl,yl,xorig,yorig)
c
c     Set up map projection grid
c
      implicit none
      integer mapproj
      real trulat1,trulat2,trulon
      real ctrlat,ctrlon,xl,yl,xorig,yorig
      real swx,swy,ctrx,ctry

      call XSTMPRJ(mapproj,trulat1,trulat2,trulon)

      CALL xlltoxy( 1,1, ctrlat,ctrlon, ctrx, ctry )

      swx = ctrx - (xl*0.5+xorig)*1000.0
      swy = ctry - (yl*0.5+yorig)*1000.0

      CALL xsetorig( 1, swx, swy)

      RETURN
      END
c
c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######            ARPS Map Projection Subsystem.            ######
c     ######                                                      ######
c     ######                     Developed by                     ######
c     ######     Center for Analysis and Prediction of Storms     ######
c     ######                University of Oklahoma                ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################
c
c      General Information
c
c      This set of subroutines allows for transformation between
c      lat-lon coordinates and any one of three map projections: Polar
c      Stereographic, Lambert Conformal or Mercator.
c
c      In order for the transformation subroutines to work, the
c      map projection must first be set up by calling setmapr.  The
c      user may wish to call setorig immediately after setmapr to
c      established an origin (given a lat-long or x-y in the default
c      system) other than the default origin (e.g., the north pole).
c
c      All lat-lons are in degrees (positive north, negative south,
c      positive east and negative west).  Note carefully the dimensions
c      of x,y -- it differs among the subroutines to conform to ARPS usage.
c      x,y coordinates are meters on earth but may be changed using the scale
c      parameter in setmapr to change to km (scale=0.001) or to a different
c      sphere (e.g., scale=mars_radius/earth_radius).
c
c
c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######                SUBROUTINE XSTMPRJ                    ######
c     ######                                                      ######
c     ######                     Developed by                     ######
c     ######     Center for Analysis and Prediction of Storms     ######
c     ######                University of Oklahoma                ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################
c
      SUBROUTINE XSTMPRJ(iproj,trulat1,trulat2,trulon)
c
c#######################################################################
c
c     PURPOSE:
c
c     Set constants for map projections, which are stored in
c     the common block named /xprojcst/.
c
c
c#######################################################################
c
c     AUTHOR: Keith Brewster
c     11/13/93.
c
c     MODIFICATION HISTORY:
c     03/30/1995 (K. Brewster)
c     Corrected error in Lambert Conformal scaling and added code to
c     allow Lambert Tangent projection (lat1=lat2 in Lambert Conformal).
c     Resulted in redefinition of projc1 for option 2.
c
c#######################################################################
c
c     INPUT:
c
c       iproj        Map projection number
c                    1=North Polar Stereographic   (-1 South Pole)
c                    2=Northern Lambert Conformal  (-2 Southern)
c                    3=Mercator
c                    4=Lat,Lon
c
c       scale        Map scale factor,  at latitude=latnot
c                    Distance on map = (Distance on earth) * scale
c                    For ARPS model runs, generally this is 1.0
c                    For ARPS plotting this will depend on window
c                    size and the area to be plotted.
c
c       latnot(2)    Real "True" latitude(s) of map projection
c                    (degrees, positive north)
c                    Except for iproj=1, only latnot(1) is used
c
c       orient       Longitude line that runs vertically on the map.
c                    (degrees, negative west, positive east)
c
c#######################################################################
c
c     Variable Declarations.
c
c#######################################################################
c
      implicit none
      integer iproj
      real trulat1, trulat2, trulon
      real scale                       ! map scale factor
      real latnot(2)                   ! true latitude (degrees N)
      real orient                      ! orientation longitude (degrees E)

      real d2rad,eradius
      parameter (d2rad=3.141592654/180.,
     :           eradius = 6371000. )  ! mean earth radius in m

      integer jproj,jpole
      real trulat(2),rota,scmap,xorig,yorig,
     :     projc1,projc2,projc3,projc4,projc5
      common /xprojcst/ jproj,jpole,trulat,rota,scmap,xorig,yorig,
     :                 projc1,projc2,projc3,projc4,projc5
c
c#######################################################################
c
c     Misc. local variables:
c
c#######################################################################
c
      real denom1,denom2,denom3

c#######################################################################
c
c
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
C     Beginning of executable code...
C
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
      scale = 1.0
      latnot(1) = trulat1
      latnot(2) = trulat2
      orient = trulon

      xorig=0.
      yorig=0.
      jproj=iabs(iproj)
      jpole=isign(1,iproj)
c
c#######################################################################
c
c     No map projection
c
c#######################################################################
c
      IF ( jproj.eq.0 ) THEN
c       write(6,'(a)')
c    :  '  No map projection will be used.'
c
c#######################################################################
c
c     Polar Stereographic projection
c     For this projection:
c         projc1 is the scaled earth's radius, scale times eradius
c         projc2 is the numerator of emfact, the map image scale factor.
c         projc3 is projc2 times the scaled earth's radius.
c
c#######################################################################
c
      ELSEIF( jproj.eq.1 ) THEN
        trulat(1)=latnot(1)
        rota=orient
        scmap=scale
        projc1=scale*eradius
        projc2=(1. + sin(d2rad*jpole*trulat(1)) )
        projc3=projc1*projc2
        IF(jpole.gt.0) THEN
c         write(6,'(a/,a)')
c    :    '  Map projection set to Polar Stereographic',
c    :    '  X origin, Y origin set to 0.,0. at the North Pole.'
        ELSE
c         write(6,'(a/,a)')
c    :    '  Map projection set to Polar Stereographic',
c    :    '  X origin, Y origin set to 0.,0. at the South Pole.'
        END IF
c
c#######################################################################
c
c     Lambert Conformal Conic Projection.
c     For this projection:
c         projc1 is the scaled earth's radius, scale times eradius/n
c         projc2 is cos of trulat(1)
c         projc3 is tan (45. - trulat/2) a const for local map scale
c         projc4 is the cone constant, n
c
c#######################################################################
c
      ELSE IF( jproj.eq.2 ) THEN
        trulat(1)=latnot(1)
        trulat(2)=latnot(2)
        rota=orient
        scmap=scale
        projc2=cos(d2rad*trulat(1))
        projc3=tan(d2rad*(45.-0.5*jpole*trulat(1)))
        denom1=cos(d2rad*trulat(2))
        denom2=tan(d2rad*(45.-0.5*jpole*trulat(2)))
        IF(denom2.ne.0.) THEN
          denom3=alog( projc3/denom2 )
        ELSE
          denom3=0.
        END IF
        IF(denom1.ne.0. and. denom3.ne.0.) THEN
          projc4=alog( projc2/denom1 ) / denom3
c         print *, '  The cone constant is : ',projc4
          IF( projc4.lt.0.) THEN
            write(6,'(a/,a,f9.2,a,f9.2,/a)')
     :    '  Warning in SETMAPR for Lambert Projection',
     :    '  For the true latitudes provided, ',
     :       trulat(1),' and ',trulat(2),
     :    '  projection must be from opposite pole...changing pole.'
            jpole=-jpole
            projc3=tan(d2rad*(45.-0.5*jpole*trulat(1)) )
            denom2=tan(d2rad*(45.-0.5*jpole*trulat(2)))
            IF(denom2.ne.0.) THEN
              denom3=alog( projc3/denom2 )
            ELSE
              denom3=0.
            END IF
            IF(denom1.ne.0. and. denom3.ne.0.) THEN
              projc4=alog( projc2/denom1 ) / denom3
c             print *, '  The revised cone constant is : ',projc4
            ELSE
              write(6,'(a/,a,f9.2,a,f9.2)')
     :      '  Error (1) in SETMAPR for Lambert Projection',
     :      '  Illegal combination of trulats one: ',
     :         trulat(1),' and two: ',trulat(2)
              STOP
            END IF
          END IF
          projc1=scale*eradius/projc4
        ELSE IF(denom3.eq.0. .and. denom2.ne.0.) THEN   ! tangent
          write(6,'(a/,a,f9.2,a,f9.2)')
     :    '  Using Tangent Lambert Projection',
     :    '  Based on input combination of trulats one: ',
     :       trulat(1),' and two: ',trulat(2)
          projc4=sin(d2rad*jpole*trulat(1))
c         print *, '  The cone constant is : ',projc4
          IF( projc4.lt.0.) THEN
            write(6,'(a/,a,f9.2,a,f9.2,/a)')
     :    '  Warning in SETMAPR for Lambert Projection',
     :    '  For the true latitudes provided, ',
     :       trulat(1),' and ',trulat(2),
     :    '  projection must be from opposite pole...changing pole.'
            jpole=-jpole
            projc4=sin(d2rad*jpole*trulat(1))
          END IF
          projc1=scale*eradius/projc4
        ELSE
          write(6,'(a/,a,f9.2,a,f9.2)')
     :    '  Error (1) in SETMAPR for Lambert Projection',
     :    '  Illegal combination of trulats one: ',
     :       trulat(1),' and two: ',trulat(2)
          STOP
        END IF

        IF(jpole.gt.0) THEN
c         write(6,'(a/,a)')
c    :    '  Map projection set to Lambert Conformal',
c    :    '  X origin, Y origin set to 0.,0. at the North Pole.'
        ELSE
c         write(6,'(a/,a)')
c    :    '  Map projection set to Lambert Conformal',
c    :    '  X origin, Y origin set to 0.,0. at the South Pole.'
        END IF
c
c#######################################################################
c
c     Mercator Projection.
c     For this projection:
c         projc1 is the scaled earth's radius, scale times eradius
c         projc2 is cos of trulat(1)
c         projc3 is projc1 times projc2
c
c#######################################################################
c
      ELSE IF( jproj.eq.3 ) THEN
        trulat(1)=latnot(1)
        rota=orient
        scmap=scale
        projc1=scale*eradius
        projc2=cos(d2rad*trulat(1))
        projc3=projc1*projc2
        IF(projc2.le.0.) THEN
          write(6,'(a/,a,f9.2,a,f9.2)')
     :    '  Error (1) in SETMAPR for Mercator Projection',
     :    '  Illegal true latitude provided: ',trulat(1)
          STOP
        END IF
        write(6,'(a/,a,f6.1/,a)')
     :    '  Map projection set to Mercator',
     :    '  X origin, Y origin set to 0.,0. at the equator,',rota,
     :    '  Y positive toward the North Pole.'
c
c#######################################################################
c
c     Lat, Lon Projection.
c     For this projection:
c         projc1 is the scaled earth's radius, scale times eradius
c         projc2 is cos of trulat(1)
c         projc3 is projc1 times projc2 times 180/pi
c
c#######################################################################
c
      ELSE IF( jproj.eq.4 ) THEN
        trulat(1)=latnot(1)
        rota=orient
        scmap=scale
        projc1=scale*eradius
        projc2=cos(d2rad*trulat(1))
        IF(projc2.le.0.) THEN
          write(6,'(a/,a,f9.2,a,f9.2)')
     :    '  Error (1) in SETMAPR for Lat,Lon Projection',
     :    '  Illegal true latitude provided: ',trulat(1)
          STOP
        END IF
        projc3=projc1*projc2/d2rad
        write(6,'(a/,a,/a)')
     :    '  Map projection set to Lat, Lon',
     :    '  X origin, Y origin set to 0.,0. at the equator, 0. long',
     :    '  Y positive toward the North Pole.'
      ELSE
        write(6,'(i4,a)') iproj,' projection is not supported'
        STOP
      END IF

      RETURN
      END
c
c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######                SUBROUTINE GETMAPR                    ######
c     ######                                                      ######
c     ######                     Developed by                     ######
c     ######     Center for Analysis and Prediction of Storms     ######
c     ######                University of Oklahoma                ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################
c
      SUBROUTINE XGETMAPR(iproj,scale,latnot,orient,x0,y0)
c
c#######################################################################
c
c     PURPOSE:
c
c     Get the constants for the current map projection, which are stored
c     in the common block named /xprojcst/.
c
c
c#######################################################################
c
c     AUTHOR: Keith Brewster
c     9/17/94.
c
c     MODIFICATION HISTORY:
c     1/17/96  Corrected retrieval of iproj to assign sign from jpole.
c
c#######################################################################
c
c     OUTPUT:
c
c       iproj        Map projection number
c                    1=North Polar Stereographic   (-1 South Pole)
c                    2=Northern Lambert Conformal  (-2 Southern)
c                    3=Mercator
c                    4=Lat,Lon
c
c       scale        Map scale factor,  at latitude=latnot
c                    Distance on map = (Distance on earth) * scale
c                    For ARPS model runs, generally this is 1.0
c                    For ARPS plotting this will depend on window
c                    size and the area to be plotted.
c
c       latnot(2)    Real "True" latitude(s) of map projection
c                    (degrees, positive north)
c                    Except for iproj=2, only latnot(1) is used
c
c       orient       Longitude line that runs vertically on the map.
c                    (degrees, negative west, positive east)
c
c       x0           x coordinate of origin
c       y0           y coordinate of origin
c
c#######################################################################
c
c     Variable Declarations.
c
c#######################################################################
c
      implicit none
      integer iproj       ! map projection number
      real scale          ! map scale factor
      real latnot(2)      ! true latitude (degrees N)
      real orient         ! orientation longitude (degrees E)
      real x0             ! x coordinate of origin
      real y0             ! y coordinate of origin

      integer jproj,jpole
      real trulat(2),rota,scmap,xorig,yorig,
     :     projc1,projc2,projc3,projc4,projc5
      common /xprojcst/ jproj,jpole,trulat,rota,scmap,xorig,yorig,
     :                 projc1,projc2,projc3,projc4,projc5

c#######################################################################
c
c
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
C     Beginning of executable code...
C
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
      iproj=jproj*jpole
      scale=scmap
      latnot(1)=trulat(1)
      latnot(2)=trulat(2)
      orient=rota
      x0=xorig
      y0=yorig
      RETURN
      END
c
c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######                SUBROUTINE XSETORIG                   ######
c     ######                                                      ######
c     ######                     Developed by                     ######
c     ######     Center for Analysis and Prediction of Storms     ######
c     ######                University of Oklahoma                ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################
c
      SUBROUTINE XSETORIG(iopt,x0,y0)
c
c#######################################################################
c
c     PURPOSE:
c
c     Set the origin for the map projection.
c     This is call after subroutine mapproj if the origin
c     must be moved from the original position, which is the
c     pole for the polar stereographic projection and the
c     Lambert conformal, and the equator for Mercator.
c
c#######################################################################
c
c     AUTHOR: Keith Brewster
c     11/20/93.
c
c     MODIFICATION HISTORY:
c
c#######################################################################
c
c     INPUT:
c
c       iopt        origin setting option
c                   1: origin given in corrdinate x,y
c                   2: origin given in lat,lon on earth
c
c       x0          first coordinate of origin
c       y0          second coordinate of origin
c
c
c#######################################################################
c
c     Variable Declarations.
c
c#######################################################################
c
      implicit none
      integer iopt       ! origin setting option
      real x0            ! first coordinate of origin
      real y0            ! second coordinate of origin

      integer jproj,jpole
      real trulat(2),rota,scmap,xorig,yorig,
     :     projc1,projc2,projc3,projc4,projc5
      common /xprojcst/ jproj,jpole,trulat,rota,scmap,xorig,yorig,
     :                 projc1,projc2,projc3,projc4,projc5
c
c#######################################################################
c
c     Misc. local variables:
c
c#######################################################################
c
      real xnew,ynew,rlat,rlon

c#######################################################################
c
c
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
C     Beginning of executable code...
C
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
c
c#######################################################################
c
c     iopt=1 origin is given in x,y in absolute coordinates.
c
c#######################################################################
c
      IF( iopt.eq.1 ) THEN
        xorig=x0
        yorig=y0
        CALL xxytoll(1,1,0.,0.,rlat,rlon)

c       write(6,'(/a,f18.2,f18.2,/a,f16.2,f16.2/)')
c    : '  Coordinate origin set to absolute x,y =',xorig,yorig,
c    : '    Latitude, longitude= ',rlat,rlon
c
c#######################################################################
c
c     iopt=2 origin is given in lat,lon on earth
c
c#######################################################################
c
c
      ELSE IF( iopt.eq.2 ) THEN
        xorig=0.
        yorig=0.
        CALL xlltoxy(1,1,x0,y0,xnew,ynew)
        xorig=xnew
        yorig=ynew
c       write(6,'(/a,f16.2,f16.2,/a,f16.2,f16.2/)')
c    : '  Coordinate origin set to absolute x,y =',xorig,yorig,
c    : '    Latitude, longitude= ',x0,y0

      ELSE
        CALL xxytoll(1,1,0.,0.,rlat,rlon)
c       write(6,'(/a,i4,a,/a,f16.2,f16.2,/a,f16.2,f16.2)')
c    : ' Setorig option ',iopt,' not supported.',
c    : '    Coordinate origin unchanged at x,y =',xorig,yorig,
c    : '    Latitude, longitude= ',rlat,rlon
      END IF
      RETURN
      END
c
c
c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######                SUBROUTINE XXYTOLL                    ######
c     ######                                                      ######
c     ######                     Developed by                     ######
c     ######     Center for Analysis and Prediction of Storms     ######
c     ######                University of Oklahoma                ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################
c
c
      SUBROUTINE XXYTOLL(idim,jdim,x,y,rlat,rlon)
c
c#######################################################################
c
c     PURPOSE:
c
c     Determine latitude and longitude given X,Y coordinates on
c     map projection.  SETMAPR must be called before this routine
c     to set-up the map projection constants.
c
c#######################################################################
c
c     AUTHOR: Keith Brewster
c     11/13/93.
c
c     MODIFICATION HISTORY:
c     01/17/96  Bug in southern hemisphere for Polar Stereo and
c               Mercator projections fixed.
c
c#######################################################################
c
c     INPUT:
c
c       idim     Number of points in x direction.
c       jdim     Number of points in y direction.
c
c       rlat     Array of latitude.
c                (degrees, negative south, positive north)
c
c       rlon     Array of longitude.
c                (degrees, negative west, positive east)
c
c     OUTPUT:
c
c       x        Vector of x in map coordinates
c       y        Vector of y in map coordinates
c                Units are meters unless the scale parameter is
c                not equal to 1.0
c
c#######################################################################
c
c     Variable Declarations.
c
c#######################################################################
c
      implicit none
      integer idim,jdim
      real x(idim),y(jdim),rlat(idim,jdim),rlon(idim,jdim)

      real r2deg,eradius
      parameter (r2deg=180./3.141592654,
     :           eradius = 6371000. )  ! mean earth radius in m
c
      integer jproj,jpole
      real trulat(2),rota,scmap,xorig,yorig,
     :     projc1,projc2,projc3,projc4,projc5
      common /xprojcst/ jproj,jpole,trulat,rota,scmap,xorig,yorig,
     :                 projc1,projc2,projc3,projc4,projc5
c
c#######################################################################
c
c     Misc. local variables:
c
c#######################################################################
c
      integer i,j
      real xabs,yabs,yjp
      real radius,ratio,dlon

c#######################################################################
c
c
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
C     Beginning of executable code...
C
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
c#######################################################################
c
c     No map projection
c
c#######################################################################
c
      IF ( jproj.eq.0 ) THEN
        ratio=r2deg/eradius
        DO 10 j = 1, jdim
        DO 10 i = 1, idim
          rlat(i,j) = ratio*(y(j)+yorig)
          rlon(i,j) = ratio*(x(i)+xorig)
10      CONTINUE
c
c#######################################################################
c
c     Polar Stereographic projection
c     For this projection:
c         projc1 is the scaled earth's radius, scale times eradius
c         projc2 is the numerator of emfact, the map image scale factor.
c         projc3 is projc2 times the scaled earth's radius.
c
c#######################################################################
c
      ELSEIF( jproj.eq.1 ) THEN
        DO 100 j=1,jdim
        DO 100 i=1,idim
          yabs=y(j)+yorig
          xabs=x(i)+xorig
          radius=sqrt( xabs*xabs + yabs*yabs )/projc3
          rlat(i,j) = jpole*(90. - 2.*r2deg*atan(radius))
          rlat(i,j)=amin1(rlat(i,j), 90.)
          rlat(i,j)=amax1(rlat(i,j),-90.)

          IF((jpole*yabs).gt.0.) THEN
            dlon=180. + r2deg*atan(-xabs/yabs)
          ELSE IF((jpole*yabs).lt.0.) THEN
            dlon=r2deg*atan(-xabs/yabs)
          ELSE IF (xabs.gt.0.) THEN     ! y=0.
            dlon=90.
          ELSE
            dlon=-90.
          END IF
          rlon(i,j)= rota + jpole*dlon
          IF(rlon(i,j).gt. 180) rlon(i,j)=rlon(i,j)-360.
          IF(rlon(i,j).lt.-180) rlon(i,j)=rlon(i,j)+360.
          rlon(i,j)=amin1(rlon(i,j), 180.)
          rlon(i,j)=amax1(rlon(i,j),-180.)
c
 100    CONTINUE
c
c#######################################################################
c
c     Lambert Conformal Conic Projection.
c     For this projection:
c         projc1 is the scaled earth's radius, scale times eradius/n
c         projc2 is cos of trulat(1)
c         projc3 is tan (45. - trulat/2) a const for local map scale
c         projc4 is the cone constant, n
c
c#######################################################################
c
      ELSE IF ( jproj.eq.2 ) THEN
        DO 200 j=1,jdim
        DO 200 i=1,idim
          yabs=y(j)+yorig
          xabs=x(i)+xorig
          radius=sqrt( xabs*xabs+ yabs*yabs )
          ratio=projc3*((radius/(projc1*projc2))**(1./projc4))
          rlat(i,j)=jpole*(90. -2.*r2deg*(atan(ratio)))
          rlat(i,j)=amin1(rlat(i,j), 90.)
          rlat(i,j)=amax1(rlat(i,j),-90.)

          yjp=jpole*yabs
          IF(yjp.gt.0.) THEN
            dlon=180. + r2deg*atan(-xabs/yabs)/projc4
          ELSE IF(yjp.lt.0.) THEN
            dlon=r2deg*atan(-xabs/yabs)/projc4
          ELSE IF (xabs.gt.0.) THEN     ! y=0.
            dlon=90./projc4
          ELSE
            dlon=-90./projc4
          END IF
          rlon(i,j)= rota + jpole*dlon
          IF(rlon(i,j).gt. 180) rlon(i,j)=rlon(i,j)-360.
          IF(rlon(i,j).lt.-180) rlon(i,j)=rlon(i,j)+360.
          rlon(i,j)=amin1(rlon(i,j), 180.)
          rlon(i,j)=amax1(rlon(i,j),-180.)

 200    CONTINUE
c
c#######################################################################
c
c     Mercator Projection.
c     For this projection:
c         projc1 is the scaled earth's radius, scale times eradius
c         projc2 is cos of trulat(1)
c         projc3 is projc1 times projc2
c
c#######################################################################
c
      ELSE IF( jproj.eq.3 ) THEN
        DO 300 j=1,jdim
        DO 300 i=1,idim
          yabs=y(j)+yorig
          xabs=x(i)+xorig
          rlat(i,j)=(90. - 2.*r2deg*atan(exp(-yabs/projc3)))
          rlat(i,j)=amin1(rlat(i,j), 90.)
          rlat(i,j)=amax1(rlat(i,j),-90.)
          dlon=r2deg*(xabs/projc3)
          rlon(i,j)=rota + dlon
          IF(rlon(i,j).gt. 180) rlon(i,j)=rlon(i,j)-360.
          IF(rlon(i,j).lt.-180) rlon(i,j)=rlon(i,j)+360.
 300    CONTINUE
c
c#######################################################################
c
c     Lat, Lon Projection.
c     For this projection:
c         projc1 is the scaled earth's radius, scale times eradius
c         projc2 is cos of trulat(1)
c         projc3 is projc1 times projc2 times 180/pi
c
c#######################################################################
c
      ELSE IF( jproj.eq.4 ) THEN
        DO 400 j=1,jdim
        DO 400 i=1,idim
          rlon(i,j)=x(j)-xorig
          rlat(i,j)=y(j)-yorig
  400   CONTINUE
      ELSE
        write(6,'(i4,a)') jproj,' projection is not supported'
        STOP
      END IF

      RETURN
      END
c
c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######                SUBROUTINE XLLTOXY                    ######
c     ######                                                      ######
c     ######                     Developed by                     ######
c     ######     Center for Analysis and Prediction of Storms     ######
c     ######                University of Oklahoma                ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################
c
      SUBROUTINE XLLTOXY(idim,jdim,rlat,rlon,xloc,yloc)
c
c#######################################################################
c
c     PURPOSE:
c
c     Determine x, y coordinates on map projection from the given latitude
c     and longitude. SETMAPR must be called before this routine to set-up
c     the map projection constants.
c
c#######################################################################
c
c     AUTHOR: Keith Brewster
c     11/11/93.
c
c     MODIFICATION HISTORY:
c
c#######################################################################
c
c     INPUT:
c
c       idim     Array dimension in x direction
c       jdim     Array dimension in y direction
c
c       rlat     Real vector of latitude.
c                (degrees, negative south, positive north)
c
c       rlon     Real vector of longitude.
c                (degrees, negative west, positive east)
c
c     OUTPUT:
c
c       xloc     Real vector of x in map coordinates
c       yloc     Real vector of y in map coordinates
c
c#######################################################################
c
c     Variable Declarations.
c
c#######################################################################
c
      implicit none

      integer idim,jdim
      real rlat(idim,jdim),rlon(idim,jdim)
      real xloc(idim,jdim),yloc(idim,jdim)

      real d2rad,eradius
      parameter (d2rad=3.141592654/180.,
     :           eradius = 6371000. )  ! mean earth radius in m

      integer jproj,jpole
      real trulat(2),rota,scmap,xorig,yorig,
     :     projc1,projc2,projc3,projc4,projc5
      common /xprojcst/ jproj,jpole,trulat,rota,scmap,xorig,yorig,
     :                 projc1,projc2,projc3,projc4,projc5
c
c#######################################################################
c
c     Misc. local variables:
c
c#######################################################################
c
      integer i,j
      real radius,denom,dlon,ratio
      real tem

c#######################################################################
c
c
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
C     Beginning of executable code...
C
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
c#######################################################################
c
c     No map projection
c
c#######################################################################
c
      IF( jproj.eq.0 ) THEN
        ratio=d2rad*eradius
        DO 10 j = 1, jdim
        DO 10 i = 1, idim
          tem = rlon(i,j)
          if( tem.lt.-180.0) tem = 360.0+tem
          if( tem.gt. 180.0) tem = tem-360.0
          xloc(i,j) = ratio*tem       - xorig
          yloc(i,j) = ratio*rlat(i,j) - yorig
10      CONTINUE
c
c#######################################################################
c
c     Polar Stereographic projection
c     For this projection:
c         projc1 is the scaled earth's radius, scale times eradius
c         projc2 is the numerator of emfact, the map image scale factor.
c         projc3 is projc2 times the scaled earth's radius.
c
c#######################################################################
c
      ELSE IF( jproj.eq.1 ) THEN
        DO 100 j=1,jdim
        DO 100 i=1,idim
          denom=(1. + sin(d2rad*jpole*rlat(i,j)))
          IF(denom.eq.0.) denom=1.0E-10
          radius=jpole*projc3*cos(d2rad*rlat(i,j))/denom
          tem = rlon(i,j)-rota
          if( tem.lt.-180.0) tem = 360.0+tem
          if( tem.gt. 180.0) tem = tem-360.0
          dlon=jpole*d2rad*tem
          xloc(i,j)= radius*sin(dlon) - xorig
          yloc(i,j)=-radius*cos(dlon) - yorig
 100    CONTINUE
c
c#######################################################################
c
c     Lambert Conformal Conic Projection.
c     For this projection:
c         projc1 is the scaled earth's radius, scale times eradius/n
c         projc2 is cos of trulat(1)
c         projc3 is tan (45. - trulat/2) a const for local map scale
c         projc4 is the cone constant, n
c
c#######################################################################
c
      ELSE IF( jproj.eq.2 ) THEN
        DO 200 j=1,jdim
        DO 200 i=1,idim
          radius=projc1*projc2
     :          *(tan(d2rad*(45.-0.5*jpole*rlat(i,j)))/projc3)**projc4
c         dlon=projc4*d2rad*(rlon(i,j)-rota)
cmx
          tem = rlon(i,j)-rota
          if( tem.lt.-180.0) tem = 360.0+tem
          if( tem.gt. 180.0) tem = tem-360.0
          dlon=projc4*d2rad*tem

          xloc(i,j)=       radius*sin(dlon) - xorig
          yloc(i,j)=-jpole*radius*cos(dlon) - yorig
 200    CONTINUE
c
c#######################################################################
c
c     Mercator Projection.
c     For this projection:
c         projc1 is the scaled earth's radius, scale times eradius
c         projc2 is cos of trulat(1)
c         projc3 is projc1 times projc2
c
c#######################################################################
c
      ELSE IF(jproj.eq.3) THEN
        DO 300 j=1,jdim
        DO 300 i=1,idim
          dlon=rlon(i,j)-rota
          IF(dlon.lt.-180.) dlon=dlon+360.
          IF(dlon.gt. 180.) dlon=dlon-360.
          xloc(i,j)=projc3*d2rad*dlon - xorig
          denom=tan(d2rad*(45. - 0.5*rlat(i,j)))
          IF( denom.le.0. ) denom=1.0E-10
          yloc(i,j)=-projc3*alog(denom) - yorig
 300    CONTINUE
c
c#######################################################################
c
c     Lat, Lon Projection.
c     For this projection:
c         projc1 is the scaled earth's radius, scale times eradius
c         projc2 is cos of trulat(1)
c         projc3 is projc1 times projc2 times 180/pi
c
c#######################################################################
c
      ELSE IF(jproj.eq.4) THEN
        DO 400 j=1,jdim
        DO 400 i=1,idim
          tem = rlon(i,j)
          if( tem.lt.-180.0) tem = 360.0+tem
          if( tem.gt. 180.0) tem = tem-360.0
          xloc(i,j)=tem      -xorig
          yloc(i,j)=rlat(i,j)-yorig
  400   CONTINUE
      ELSE
        write(6,'(i4,a)') jproj,' projection is not supported'
        STOP
      END IF
      RETURN
      END
c
c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######                SUBROUTINE XLATTOMF                   ######
c     ######                                                      ######
c     ######                     Developed by                     ######
c     ######     Center for Analysis and Prediction of Storms     ######
c     ######                University of Oklahoma                ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################
c
      SUBROUTINE XLATTOMF(idim,jdim,rlat,emfact)
c
c#######################################################################
c
c     PURPOSE:
c
c     Determine the map scale factor, emfact, at a given latitude.
c
c#######################################################################
c
c     AUTHOR: Keith Brewster
c     11/11/93.
c
c     MODIFICATION HISTORY:
c
c#######################################################################
c
c     INPUT:
c
c       idim        Array dimension in x direction
c       jdim        Array dimension in y direction
c
c       rlat        Real vector of latitudes.
c                   (degrees, negative south, positive north)
c
c     OUTPUT:
c
c       emfact      Vector of map scale factors corresponding to the
c                   input latitudes (map scale includes the projection
c                   image scale times the overall scale of the map).
c
c#######################################################################
c
c     Variable Declarations.
c
c#######################################################################
c
      implicit none

      integer idim,jdim         ! dimensions of arrays
      real rlat(idim,jdim)      ! latitude (degrees)
      real emfact(idim,jdim)    ! local map scale factor

      real d2rad
      parameter (d2rad=3.141592654/180.)

      integer jproj,jpole
      real trulat(2),rota,scmap,xorig,yorig,
     :     projc1,projc2,projc3,projc4,projc5
      common /xprojcst/ jproj,jpole,trulat,rota,scmap,xorig,yorig,
     :                 projc1,projc2,projc3,projc4,projc5
c
c#######################################################################
c
c     Misc. local variables:
c
c#######################################################################
c
      integer i,j
      real denom

c#######################################################################
c
c
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
C     Beginning of executable code...
C
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
c#######################################################################
c
c     No map projection
c
c#######################################################################
c
      IF( jproj.eq.0 ) THEN
        DO 10 j=1,jdim
        DO 10 i=1,idim
          emfact(i,j)=1.0
 10     CONTINUE
c
c#######################################################################
c
c     Polar Stereographic projection
c     For this projection:
c         projc1 is the scaled earth's radius, scale times eradius
c         projc2 is the numerator of emfact, the map image scale factor.
c         projc3 is projc2 times the scaled earth's radius.
c
c#######################################################################
c
      ELSE IF( jproj.eq.1 ) THEN
        DO 100 j=1,jdim
        DO 100 i=1,idim
          denom=(1. + sin(d2rad*jpole*rlat(i,j)))
          IF(denom.eq.0.) denom=1.0E-10
          emfact(i,j)=scmap*projc2/denom
 100    CONTINUE
c
c#######################################################################
c
c     Lambert Conformal Conic Projection.
c     For this projection:
c         projc1 is the scaled earth's radius, scale times eradius/n
c         projc2 is cos of trulat(1)
c         projc3 is tan (45. - trulat/2) a const for local map scale
c         projc4 is the cone constant, n
c
c#######################################################################
c
      ELSE IF( jproj.eq.2 ) THEN
        DO 200 j=1,jdim
        DO 200 i=1,idim
          denom=cos( d2rad*rlat(i,j) )
          IF(denom.lt.1.0E-06) THEN
            emfact(i,j)=1.0e+10
          ELSE
            emfact(i,j)=scmap*(projc2/denom)
     :               *(tan(d2rad*(45.-0.5*jpole*rlat(i,j)))
     :               /projc3)**projc4
          END IF
          emfact(i,j)=amax1(emfact(i,j),1.0e-10)
          emfact(i,j)=amin1(emfact(i,j),1.0e+10)
 200    CONTINUE
c
c#######################################################################
c
c     Mercator Projection.
c     For this projection:
c         projc1 is the scaled earth's radius, scale times eradius
c         projc2 is cos of trulat(1)
c
c#######################################################################
c
      ELSE IF(jproj.eq.3) THEN
        DO 300 j=1,jdim
        DO 300 i=1,idim
          denom=cos( d2rad*rlat(i,j) )
          IF(denom.eq.0.) denom=1.0E-10
          emfact(i,j)=projc2/denom
 300    CONTINUE
c
c#######################################################################
c
c     Lat, Lon Projection.
c     For this projection:
c         projc1 is the scaled earth's radius, scale times eradius
c         projc2 is cos of trulat(1)
c         projc3 is projc1 times projc2 times 180/pi
c
c#######################################################################
c
      ELSE IF(jproj.eq.4) THEN
        DO 400 j=1,jdim
        DO 400 i=1,idim
          denom=cos( d2rad*rlat(i,j) )
          IF(denom.eq.0.) denom=1.0E-10
          emfact(i,j)=projc3/denom
  400   CONTINUE
      ELSE
        write(6,'(i4,a)') jproj,' projection is not supported'
        STOP
      END IF
      RETURN
      END
c
c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######                SUBROUTINE XXYTOMF                    ######
c     ######                                                      ######
c     ######                     Developed by                     ######
c     ######     Center for Analysis and Prediction of Storms     ######
c     ######                University of Oklahoma                ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################
c
      SUBROUTINE XXYTOMF(idim,jdim,x,y,emfact)
c
c#######################################################################
c
c     PURPOSE:
c
c     Determine the map scale factor, emfact, given x,y in the projected
c     space.
c
c#######################################################################
c
c     AUTHOR: Keith Brewster
c     11/11/93.
c
c     MODIFICATION HISTORY:
c
c#######################################################################
c
c     INPUT:
c
c       idim     Array dimension in x direction.
c       jdim     Array dimension in y direction.
c
c       x        x coordinate values (meters if scmap=1.0)
c       y        y coordinate values (meters if scmap=1.0)
c
c     OUTPUT:
c
c       emfact    Vector of map scale factors corresponding to the
c                input x,y's.
c
c#######################################################################
c
c     Variable Declarations.
c
c#######################################################################
c
      implicit none

      integer idim,jdim       ! array dimensions
      real x(idim)            ! x map coordinate
      real y(jdim)            ! y map coordinate
      real emfact(idim,jdim)  ! local map scale factor

      real d2rad,r2deg
      parameter (d2rad=3.141592654/180.,
     :           r2deg=180./3.141592654)

      integer jproj,jpole
      real trulat(2),rota,scmap,xorig,yorig,
     :     projc1,projc2,projc3,projc4,projc5
      common /xprojcst/ jproj,jpole,trulat,rota,scmap,xorig,yorig,
     :                 projc1,projc2,projc3,projc4,projc5
c
c#######################################################################
c
c     Misc. local variables:
c
c#######################################################################
c
      integer i,j
      real xabs,yabs,rlat,ratio,radius,denom
c
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
C     Beginning of executable code...
C
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
c#######################################################################
c
c     No map projection
c
c#######################################################################
      IF( jproj.eq.0 ) THEN
        DO 10 j=1,jdim
        DO 10 i=1,idim
          emfact(i,j)=1.0
  10    CONTINUE
c
c#######################################################################
c
c     Polar Stereographic projection
c     For this projection:
c         projc1 is the scaled earth's radius, scale times eradius
c         projc2 is the numerator of emfact, the map image scale factor.
c         projc3 is projc2 times the scaled earth's radius.
c
c#######################################################################
c
      ELSE IF( jproj.eq.1 ) THEN
        DO 100 j=1,jdim
        DO 100 i=1,idim
          xabs=x(i)+xorig
          yabs=y(j)+yorig
          radius=sqrt( xabs*xabs + yabs*yabs )/projc3
          rlat = 90. - 2.*r2deg*atan(radius)
          rlat=amin1(rlat, 90.)
          rlat=amax1(rlat,-90.)
          denom=(1. + sin(d2rad*rlat))
          IF(denom.eq.0.) denom=1.0E-10
          emfact(i,j)=scmap*projc2/denom
 100    CONTINUE
c
c#######################################################################
c
c     Lambert Conformal Conic Projection.
c     For this projection:
c         projc1 is the scaled earth's radius, scale times eradius/n
c         projc2 is cos of trulat(1)
c         projc3 is tan (45. - trulat/2) a const for local map scale
c         projc4 is the cone constant, n
c
c#######################################################################
c
      ELSE IF( jproj.eq.2 ) THEN
        DO 200 j=1,jdim
        DO 200 i=1,idim
          xabs=x(i)+xorig
          yabs=y(j)+yorig
          radius=sqrt( xabs*xabs+ yabs*yabs )
          ratio=projc3*((radius/(projc1*projc2))**(1./projc4))
          rlat=90. -2.*r2deg*(atan(ratio))
          rlat=amin1(rlat, 90.)
          rlat=amax1(rlat,-90.)
          denom=cos( d2rad*rlat )
          IF(denom.eq.0.) denom=1.0E-10
          emfact(i,j)=scmap*(projc2/denom)
     :               *(tan(d2rad*(45.-0.5*rlat))/projc3)**projc4
 200    CONTINUE
c
c#######################################################################
c
c     Mercator Projection.
c     For this projection:
c         projc1 is the scaled earth's radius, scale times eradius
c         projc2 is cos of trulat(1)
c         projc3 is projc1 times projc2
c
c#######################################################################
c
      ELSE IF(jproj.eq.3) THEN
        DO 300 j=1,jdim
          yabs=y(j)+yorig
          rlat=90. - 2.*r2deg*atan(exp(-yabs/projc3))
          rlat=amin1(rlat, 90.)
          rlat=amax1(rlat,-90.)
          denom=cos( d2rad*rlat )
          IF(denom.eq.0.) denom=1.0E-10
          DO 300 i=1,idim
            emfact(i,j)=projc2/denom
 300    CONTINUE
c
c#######################################################################
c
c     Lat, Lon Projection.
c     For this projection:
c         projc1 is the scaled earth's radius, scale times eradius
c         projc2 is cos of trulat(1)
c         projc3 is projc1 times projc2 times 180/pi
c
c#######################################################################
c
      ELSE IF(jproj.eq.4) THEN
        DO 400 j=1,jdim
          yabs=y(j)+yorig
          denom=cos( d2rad*yabs )
          IF(denom.eq.0.) denom=1.0E-10
          DO 400 i=1,idim
            emfact(i,j)=projc3/denom
 400    CONTINUE
      ELSE
        write(6,'(i4,a)') jproj,' projection is not supported'
        STOP
      END IF
      RETURN
      END
c
c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######                SUBROUTINE XDDROTUV                   ######
c     ######                                                      ######
c     ######                     Developed by                     ######
c     ######     Center for Analysis and Prediction of Storms     ######
c     ######                University of Oklahoma                ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################
c
      SUBROUTINE XDDROTUV(nsta,stalon,dd,ff,ddrot,umap,vmap)
c
c#######################################################################
c
c     PURPOSE:
c
c     Rotate wind from earth direction to map orientation.
c
c#######################################################################
c
c     AUTHOR: Keith Brewster
c     11/20/93.
c
c     MODIFICATION HISTORY:
c     03/30/95  (K. Brewster)
c     Removed the map scale factor from the conversion of winds
c     from u,v on the earth to projection u,v.  Affected argument
c     list of ddrotuv.
c
c#######################################################################
c
c     INPUT:
c
c       nsta      array dimension
c
c       stalon    longitude (degrees E)
c
c       dd        wind direction (degrees from north)
c       ff        wind speed
c
c     OUTPUT:
c
c       ddrot     wind direction rotated to map orientation
c
c       umap      u wind component on map (same units as ff)
c       vmap      v wind component on map (same units as ff)
c
c#######################################################################
c
c     Variable Declarations.
c
c#######################################################################
c
      implicit none

      integer nsta               ! array dimension
      real stalon(nsta)          ! longitude (degrees E)
      real dd(nsta)              ! wind direction
      real ff(nsta)              ! speed
      real ddrot(nsta)           ! wind direction rotated to map orientation
      real umap(nsta)            ! u wind component on map
      real vmap(nsta)            ! v wind component on map

      real d2rad,r2deg
      parameter (d2rad=3.141592654/180.,
     :           r2deg=180./3.141592654)

      integer jproj,jpole
      real trulat(2),rota,scmap,xorig,yorig,
     :     projc1,projc2,projc3,projc4,projc5
      common /xprojcst/ jproj,jpole,trulat,rota,scmap,xorig,yorig,
     :                 projc1,projc2,projc3,projc4,projc5
c
c#######################################################################
c
c     Misc. local variables:
c
c#######################################################################
c
      integer i
      real arg

c#######################################################################
c
c
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
C     Beginning of executable code...
C
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
c#######################################################################
c
c     No map projection.
c     Just do conversion from ddff to u,v.
c
c#######################################################################
c
      IF( jproj.eq.0 ) THEN
        DO 50 i=1,nsta
          ddrot(i)=dd(i)
          arg = (ddrot(i) * d2rad)
          umap(i) = -ff(i) * sin(arg)
          vmap(i) = -ff(i) * cos(arg)
  50    CONTINUE
c
c#######################################################################
c
c     Polar Stereographic projection
c     For this projection:
c         projc1 is the scaled earth's radius, scale times eradius
c         projc2 is the numerator of emfact, the map image scale factor.
c         projc3 is projc2 times the scaled earth's radius.
c
c#######################################################################
c
      ELSE IF( jproj.eq.1 ) THEN
        DO 100 i=1,nsta
          ddrot(i)=dd(i) + rota - stalon(i)
          arg = (ddrot(i) * d2rad)
          umap(i) = -ff(i) * sin(arg)
          vmap(i) = -ff(i) * cos(arg)
 100    CONTINUE
c
c#######################################################################
c
c     Lambert Conformal Conic Projection.
c     For this projection:
c         projc1 is the scaled earth's radius, scale times eradius/n
c         projc2 is cos of trulat(1)
c         projc3 is tan (45. - trulat/2) a const for local map scale
c         projc4 is the cone constant, n
c
c#######################################################################
c
      ELSE IF( jproj.eq.2 ) THEN
        DO 200 i=1,nsta
          ddrot(i)=dd(i) + projc4*(rota - stalon(i))
          arg = (ddrot(i) * d2rad)
          umap(i) = -ff(i) * sin(arg)
          vmap(i) = -ff(i) * cos(arg)
 200    CONTINUE
c
c#######################################################################
c
c     Mercator Projection.
c     For this projection:
c         projc1 is the scaled earth's radius, scale times eradius
c         projc2 is cos of trulat(1)
c         projc3 is projc1 times projc2
c
c#######################################################################
c
      ELSE IF(jproj.eq.3) THEN
        DO 300 i=1,nsta
          ddrot(i)=dd(i)
          arg = (ddrot(i) * d2rad)
          umap(i) = -ff(i) * sin(arg)
          vmap(i) = -ff(i) * cos(arg)
 300    CONTINUE
c
c#######################################################################
c
c     Lat, Lon Projection.
c     For this projection:
c         projc1 is the scaled earth's radius, scale times eradius
c         projc2 is cos of trulat(1)
c         projc3 is projc1 times projc2 times 180/pi
c
c#######################################################################
c
      ELSE IF(jproj.eq.4) THEN
        DO 400 i=1,nsta
          ddrot(i)=dd(i)
          arg = (ddrot(i) * d2rad)
          umap(i) = -ff(i) * sin(arg)
          vmap(i) = -ff(i) * cos(arg)
  400   CONTINUE
      ELSE
        write(6,'(i4,a)') jproj,' projection is not supported'
        STOP
      END IF
      RETURN
      END
c
c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######                SUBROUTINE XUVROTDD                   ######
c     ######                                                      ######
c     ######                     Developed by                     ######
c     ######     Center for Analysis and Prediction of Storms     ######
c     ######                University of Oklahoma                ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################
c
      SUBROUTINE XUVROTDD(idim,jdim,elon,umap,vmap,dd,ff)
c
c#######################################################################
c
c     PURPOSE:
c     Convert winds u, v in map coordinates to wind direction and speed
c     in earth coordinates.
c
c#######################################################################
c
c     AUTHOR: Keith Brewster
c     11/20/93.
c
c     MODIFICATION HISTORY:
c     03/30/95  (K. Brewster)
c     Removed the map scale factor from the conversion of winds
c     from u,v on the earth to projection u,v.  Affected argument
c     list of uvrotdd.
c
c#######################################################################
c
c     INPUT:
c       idim       Array dimension in the x direction
c       jdim       Array dimension in the y direction
c
c       elon       Earth longitude (degrees E)
c
c       umap       u wind component on map
c       vmap       v wind component on map
c
c     OUTPUT:
c       dd         wind direction on earth
c       ff         wind speed on earth
c
c#######################################################################
c
c     Variable Declarations.
c
c#######################################################################
c
      implicit none

      integer idim,jdim       ! array dimensions
      real elon(idim,jdim)    ! longitude (degrees E)
      real umap(idim,jdim)    ! u wind component on map
      real vmap(idim,jdim)    ! v wind component on map

      real dd(idim,jdim)      ! direction
      real ff(idim,jdim)      ! wind speed

      real d2rad,r2deg
      parameter (d2rad=3.141592654/180.,
     :           r2deg=180./3.141592654)

      integer jproj,jpole
      real trulat(2),rota,scmap,xorig,yorig,
     :     projc1,projc2,projc3,projc4,projc5
      common /xprojcst/ jproj,jpole,trulat,rota,scmap,xorig,yorig,
     :                 projc1,projc2,projc3,projc4,projc5
c
c#######################################################################
c
c     Misc. local variables:
c
c#######################################################################
c
      integer i,j
      real dlon

c#######################################################################
c
c
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
C     Beginning of executable code...
C
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
c#######################################################################
c
c     No map projection
c
c#######################################################################
c
      IF( jproj.eq.0 ) THEN
        DO 50 j=1,jdim
        DO 50 i=1,idim
          ff(i,j) = sqrt(umap(i,j)*umap(i,j) + vmap(i,j)*vmap(i,j))

          IF(vmap(i,j).gt.0.) THEN
            dlon=r2deg*atan(umap(i,j)/vmap(i,j))
          ELSE IF(vmap(i,j).lt.0.) THEN
            dlon=180. + r2deg*atan(umap(i,j)/vmap(i,j))
          ELSE IF(umap(i,j).ge.0.) THEN
            dlon=90.
          ELSE
            dlon=-90.
          END IF

          dd(i,j)= dlon + 180.
          dd(i,j)= dd(i,j)-360.*(nint(dd(i,j))/360)
  50    CONTINUE
c
c#######################################################################
c
c     Polar Stereographic projection
c     For this projection:
c         projc1 is the scaled earth's radius, scale times eradius
c         projc2 is the numerator of emfact, the map image scale factor.
c         projc3 is projc2 times the scaled earth's radius.
c
c#######################################################################
c
      ELSE IF( jproj.eq.1 ) THEN
        DO 100 j=1,jdim
        DO 100 i=1,idim

          ff(i,j) = sqrt(umap(i,j)*umap(i,j) + vmap(i,j)*vmap(i,j))

          IF(vmap(i,j).gt.0.) THEN
            dlon=r2deg*atan(umap(i,j)/vmap(i,j))
          ELSE IF(vmap(i,j).lt.0.) THEN
            dlon=180. + r2deg*atan(umap(i,j)/vmap(i,j))
          ELSE IF(umap(i,j).ge.0.) THEN
            dlon=90.
          ELSE
            dlon=-90.
          END IF

          dd(i,j)= dlon + 180. + elon(i,j) - rota
          dd(i,j)= dd(i,j)-360.*(nint(dd(i,j))/360)
 100    CONTINUE
c
c#######################################################################
c
c     Lambert Conformal Conic Projection.
c     For this projection:
c         projc1 is the scaled earth's radius, scale times eradius/n
c         projc2 is cos of trulat(1)
c         projc3 is tan (45. - trulat/2) a const for local map scale
c         projc4 is the cone constant, n
c
c#######################################################################
c
      ELSE IF( jproj.eq.2 ) THEN
        DO 200 j=1,jdim
        DO 200 i=1,idim
          ff(i,j) = sqrt(umap(i,j)*umap(i,j) + vmap(i,j)*vmap(i,j))

          IF(vmap(i,j).gt.0.) THEN
            dlon=r2deg*atan(umap(i,j)/vmap(i,j))
          ELSE IF(vmap(i,j).lt.0.) THEN
            dlon=180. + r2deg*atan(umap(i,j)/vmap(i,j))
          ELSE IF(umap(i,j).ge.0.) THEN
            dlon=90.
          ELSE
            dlon=-90.
          END IF

          dd(i,j)= dlon + 180. + projc4*(elon(i,j) - rota)
          dd(i,j)= dd(i,j)-360.*(nint(dd(i,j))/360)
 200    CONTINUE
c
c#######################################################################
c
c     Mercator Projection.
c     For this projection:
c         projc1 is the scaled earth's radius, scale times eradius
c         projc2 is cos of trulat(1)
c         projc3 is projc1 times projc2
c
c#######################################################################
c
      ELSE IF(jproj.eq.3) THEN
        DO 300 j=1,jdim
        DO 300 i=1,idim
          ff(i,j) = sqrt(umap(i,j)*umap(i,j) + vmap(i,j)*vmap(i,j))

          IF(vmap(i,j).gt.0.) THEN
            dlon=r2deg*atan(umap(i,j)/vmap(i,j))
          ELSE IF(vmap(i,j).lt.0.) THEN
            dlon=180. + r2deg*atan(umap(i,j)/vmap(i,j))
          ELSE IF(umap(i,j).ge.0.) THEN
            dlon=90.
          ELSE
            dlon=-90.
          END IF

          dd(i,j)= dlon + 180.
          dd(i,j)= dd(i,j)-360.*(nint(dd(i,j))/360)
 300    CONTINUE
c
c#######################################################################
c
c     Lat, Lon Projection.
c     For this projection:
c         projc1 is the scaled earth's radius, scale times eradius
c         projc2 is cos of trulat(1)
c         projc3 is projc1 times projc2 times 180/pi
c
c#######################################################################
c
      ELSE IF(jproj.eq.4) THEN
        DO 400 j=1,jdim
        DO 400 i=1,idim
          ff(i,j) = sqrt(umap(i,j)*umap(i,j) + vmap(i,j)*vmap(i,j))

          IF(vmap(i,j).gt.0.) THEN
            dlon=r2deg*atan(umap(i,j)/vmap(i,j))
          ELSE IF(vmap(i,j).lt.0.) THEN
            dlon=180. + r2deg*atan(umap(i,j)/vmap(i,j))
          ELSE IF(umap(i,j).ge.0.) THEN
            dlon=90.
          ELSE
            dlon=-90.
          END IF

          dd(i,j)= dlon + 180.
          dd(i,j)= dd(i,j)-360.*(nint(dd(i,j))/360)
  400   CONTINUE
      ELSE
        write(6,'(i4,a)') jproj,' projection is not supported'
        STOP
      END IF
      RETURN
      END
c
c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######                SUBROUTINE XUVETOMP                   ######
c     ######                                                      ######
c     ######                     Developed by                     ######
c     ######     Center for Analysis and Prediction of Storms     ######
c     ######                University of Oklahoma                ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################
c
      SUBROUTINE XUVETOMP(idim,jdim,uear,vear,lon,umap,vmap)
c
c#######################################################################
c
c     PURPOSE:
c
c     Transform u, v wind from earth coordinates to map coordinates.
c
c#######################################################################
c
c     AUTHOR: Keith Brewster
c     04/30/94.
c
c     MODIFICATION HISTORY:
c     03/30/95  (K. Brewster)
c     Removed the map scale factor from the conversion of winds
c     from u,v on the earth to projection u,v.  Affected argument
c     list of uvetomp.
c     04/30/96  (KB)
c     Streamlined the computation for iproj=1 and iproj=2.
c     12/11/96  (KB)
c     Corrected a bug in the computation for iproj=1 and iproj=2.
c
c#######################################################################
c
c     INPUT:
c
c       idim       Array dimension in the x direction
c       jdim       Array dimension in the y direction
c
c       uear       u (eastward) wind component on earth
c       vear       v (northwrd) wind component on earth
c
c       lon        earth longitude
c
c     OUTPUT:
c
c       umap       u wind component on map
c       vmap       v wind component on map
c
c#######################################################################
c
c     Variable Declarations.
c
c#######################################################################
c
      implicit none

      integer idim,jdim        ! array dimensions
      real uear(idim,jdim)     ! u (eastward) wind component on earth
      real vear(idim,jdim)     ! v (northward) wind component on earth
      real lon(idim,jdim)      ! longitude (degrees east)

      real umap(idim,jdim)     ! u wind component on map
      real vmap(idim,jdim)     ! v wind component on map

      real d2rad,r2deg
      parameter (d2rad=3.141592654/180.,
     :           r2deg=180./3.141592654)

      integer jproj,jpole
      real trulat(2),rota,scmap,xorig,yorig,
     :     projc1,projc2,projc3,projc4,projc5
      common /xprojcst/ jproj,jpole,trulat,rota,scmap,xorig,yorig,
     :                 projc1,projc2,projc3,projc4,projc5
c
c#######################################################################
c
c     Misc. local variables:
c
c#######################################################################
c
      integer i,j
      real dlon,arg,dxdlon,dydlon,utmp,vtmp
c
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
C     Beginning of executable code...
C
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
c#######################################################################
c
c     No map projection
c
c#######################################################################
c
      IF( jproj.eq.0 ) THEN
        DO 50 j=1,jdim
        DO 50 i=1,idim
          umap(i,j) = uear(i,j)
          vmap(i,j) = vear(i,j)
  50    CONTINUE
c
c#######################################################################
c
c     Polar Stereographic projection
c     For this projection:
c         projc1 is the scaled earth's radius, scale times eradius
c         projc2 is the numerator of emfact, the map image scale factor.
c         projc3 is projc2 times the scaled earth's radius.
c
c#######################################################################
c
      ELSE IF( jproj.eq.1 ) THEN
        DO 100 j=1,jdim
        DO 100 i=1,idim
          dlon=(lon(i,j)-rota)
          arg=d2rad*dlon
          dxdlon=cos(arg)
          dydlon=sin(arg)
          utmp=uear(i,j)
          vtmp=vear(i,j)
          umap(i,j)=utmp*dxdlon - vtmp*dydlon
          vmap(i,j)=vtmp*dxdlon + utmp*dydlon
 100    CONTINUE
c
c#######################################################################
c
c     Lambert Conformal Conic Projection.
c     For this projection:
c         projc1 is the scaled earth's radius, scale times eradius/n
c         projc2 is cos of trulat(1)
c         projc3 is tan (45. - trulat/2) a const for local map scale
c         projc4 is the cone constant, n
c
c#######################################################################
c
      ELSE IF( jproj.eq.2 ) THEN
        DO 200 j=1,jdim
        DO 200 i=1,idim
          dlon=(lon(i,j)-rota)
          arg=d2rad*projc4*(dlon - 360.*nint(dlon/360.))
          dxdlon=cos(arg)
          dydlon=sin(arg)
          utmp=uear(i,j)
          vtmp=vear(i,j)
          umap(i,j)=utmp*dxdlon - vtmp*dydlon
          vmap(i,j)=vtmp*dxdlon + utmp*dydlon
 200    CONTINUE
c
c#######################################################################
c
c     Mercator Projection.
c     For this projection:
c         projc1 is the scaled earth's radius, scale times eradius
c         projc2 is cos of trulat(1)
c         projc3 is projc1 times projc2
c
c#######################################################################
c
      ELSE IF(jproj.eq.3) THEN
        DO 300 j=1,jdim
        DO 300 i=1,idim
          umap(i,j) = uear(i,j)
          vmap(i,j) = vear(i,j)
 300    CONTINUE
c
c#######################################################################
c
c     Lat, Lon Projection.
c     For this projection:
c         projc1 is the scaled earth's radius, scale times eradius
c         projc2 is cos of trulat(1)
c         projc3 is projc1 times projc2 times 180/pi
c
c#######################################################################
c
      ELSE IF(jproj.eq.4) THEN
        DO 400 j=1,jdim
        DO 400 i=1,idim
          umap(i,j) = uear(i,j)
          vmap(i,j) = vear(i,j)
 400    CONTINUE
      ELSE
        write(6,'(i4,a)') jproj,' projection is not supported'
        STOP
      END IF
      RETURN
      END
c
c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######                SUBROUTINE XUVMPTOE                   ######
c     ######                                                      ######
c     ######                     Developed by                     ######
c     ######     Center for Analysis and Prediction of Storms     ######
c     ######                University of Oklahoma                ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################
c
      SUBROUTINE XUVMPTOE(idim,jdim,umap,vmap,lon,uear,vear)
c
c#######################################################################
c
c     PURPOSE:
c
c     Transform u, v wind from map coordinates to earth coordinates.
c
c#######################################################################
c
c     AUTHOR: Keith Brewster
c     04/30/94.
c
c     MODIFICATION HISTORY:
c     03/30/95  (K. Brewster)
c     Removed the map scale factor from the conversion of winds
c     from u,v on the map to earth u,v.  Affected argument
c     list of uvmptoe.
c     04/30/96  (KB)
c     Streamlined the computation for iproj=1 and iproj=2.
c     12/11/96  (KB)
c     Corrected a bug in the computation for iproj=1 and iproj=2.
c
c#######################################################################
c
c     INPUT:
c
c       idim           Array dimension in x direction
c       jdim           Array dimension in y direction
c
c       umap           u wind component on map
c       vmap           v wind component on map
c
c       lon            Longitude (degrees E)
c
c     OUTPUT:
c
c       uear           u (eastward) wind component on earth
c       vear           v (northward) wind component on earth
c
c#######################################################################
c
c     Variable Declarations.
c
c#######################################################################
c
      implicit none

      integer idim,jdim       ! array dimensions
      real lon(idim,jdim)     ! longitude (degrees E)
      real umap(idim,jdim)    ! u wind component on map
      real vmap(idim,jdim)    ! v wind component on map

      real uear(idim,jdim)    ! u (eastward) wind component on earth
      real vear(idim,jdim)    ! v (northward) wind component on earth

      real d2rad,r2deg
      parameter (d2rad=3.141592654/180.,
     :           r2deg=180./3.141592654)

      integer jproj,jpole
      real trulat(2),rota,scmap,xorig,yorig,
     :     projc1,projc2,projc3,projc4,projc5
      common /xprojcst/ jproj,jpole,trulat,rota,scmap,xorig,yorig,
     :                 projc1,projc2,projc3,projc4,projc5
c
c#######################################################################
c
c     Misc. local variables:
c
c#######################################################################
c
      integer i,j
      real dlon,arg,utmp,vtmp,dxdlon,dydlon
c
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
C     Beginning of executable code...
C
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c
c#######################################################################
c
c     No map projection
c
c#######################################################################
c
      IF( jproj.eq.0 ) THEN
        DO 50 j=1,jdim
        DO 50 i=1,idim
          uear(i,j) = umap(i,j)
          vear(i,j) = vmap(i,j)
  50    CONTINUE
c
c#######################################################################
c
c     Polar Stereographic projection
c     For this projection:
c         projc1 is the scaled earth's radius, scale times eradius
c         projc2 is the numerator of emfact, the map image scale factor.
c         projc3 is projc2 times the scaled earth's radius.
c
c#######################################################################
c
      ELSE IF( jproj.eq.1 ) THEN
        DO 100 j=1,jdim
        DO 100 i=1,idim
          dlon=(lon(i,j)-rota)
          arg=d2rad*dlon
          dxdlon=cos(arg)
          dydlon=sin(arg)
          utmp=umap(i,j)
          vtmp=vmap(i,j)
          uear(i,j)=utmp*dxdlon + vtmp*dydlon
          vear(i,j)=vtmp*dxdlon - utmp*dydlon
 100    CONTINUE
c
c#######################################################################
c
c     Lambert Conformal Conic Projection.
c     For this projection:
c         projc1 is the scaled earth's radius, scale times eradius/n
c         projc2 is cos of trulat(1)
c         projc3 is tan (45. - trulat/2) a const for local map scale
c         projc4 is the cone constant, n
c
c#######################################################################
c
      ELSE IF( jproj.eq.2 ) THEN
        DO 200 j=1,jdim
        DO 200 i=1,idim
          dlon=(lon(i,j)-rota)
          arg=d2rad*projc4*(dlon - 360.*nint(dlon/360.))
          dxdlon=cos(arg)
          dydlon=sin(arg)
          utmp=umap(i,j)
          vtmp=vmap(i,j)
          uear(i,j)=utmp*dxdlon + vtmp*dydlon
          vear(i,j)=vtmp*dxdlon - utmp*dydlon
 200    CONTINUE
c
c#######################################################################
c
c     Mercator Projection.
c     For this projection:
c         projc1 is the scaled earth's radius, scale times eradius
c         projc2 is cos of trulat(1)
c         projc3 is projc1 times projc2
c
c#######################################################################
c
      ELSE IF(jproj.eq.3) THEN
        DO 300 j=1,jdim
        DO 300 i=1,idim
          uear(i,j) = umap(i,j)
          vear(i,j) = vmap(i,j)
 300    CONTINUE
c
c#######################################################################
c
c     Lat, Lon Projection.
c     For this projection:
c         projc1 is the scaled earth's radius, scale times eradius
c         projc2 is cos of trulat(1)
c         projc3 is projc1 times projc2 times 180/pi
c
c#######################################################################
c
      ELSE IF(jproj.eq.4) THEN
        DO 400 j=1,jdim
        DO 400 i=1,idim
          uear(i,j) = umap(i,j)
          vear(i,j) = vmap(i,j)
 400    CONTINUE
      ELSE
        write(6,'(i4,a)') jproj,' projection is not supported'
        STOP
      END IF
      RETURN
      END

      SUBROUTINE XDRAWMAP_old(nunit, mapfile, latgrid, longrid)
c
c-----------------------------------------------------------------------
c     This subroutine will draw a map within a rectagular box in
c     a map projection space. The map projection and plotting space
c     should have been properly set before calling this subroutine.
c-----------------------------------------------------------------------
c
c     nunit     the channel of the mapfile data
c     mapfile   character of map file name
c     latgrid,longrid (degree): the intervals between lat and lon grid lines.
c             < 0.0, no grid lines in the given direction,
c             = 0.0, internally determined,
c             = any real number, typically from 1.0 to 10.0 degrees.
c-----------------------------------------------------------------------
c     Author: Ming Xue
c-----------------------------------------------------------------------
      implicit none
      integer nunit
      character mapfile*(*)
      real latgrid,longrid

      integer nmax
      parameter (nmax = 100)
      real xloc(nmax),yloc(nmax),lat(nmax),long(nmax)
      real lonmin,lonmax,latmin,latmax
      real x1,x2,y1,y2,xw1,xw2,yw1,yw2
      integer iseg,ilast,ndata,i,j
      real pi2deg,x,y
      integer jmax,imax,iwndwon
      integer lsample
      parameter ( lsample = 21)
      real xlat(lsample),ylon(lsample)
     :     ,rlat(lsample,lsample),rlon(lsample,lsample)
      real latgrid1,longrid1

      real trulon
      common /pass/trulon

      CALL xqmap (x1,x2,y1,y2)
      CALL xqwdwon( iwndwon )
      CALL xqwndw(xw1,xw2,yw1,yw2)

      CALL xwindw(x1,x2,y1,y2)

      OPEN(nunit,file=mapfile,form='formatted',status='old')

      ILAST=-999
      ndata = 0
      pi2deg = 180.0/3.1415926535

  200 CONTINUE

      READ (nunit,*,END=900) ISEG,X,Y

      IF( iseg.ne.ilast ) GOTO 300

      ndata = ndata+1
      long(ndata) = x*pi2deg
      lat (ndata) = y*pi2deg

      IF(ndata.eq.nmax ) GOTO 300

      GOTO 200

300   CONTINUE

      IF( ndata.gt.1) THEN
        CALL xlltoxy(ndata,1,lat,long,xloc,yloc)
        CALL xpenup(xloc(1)*0.001 ,yloc(1)*0.001 )
        DO 350 i=2,ndata
          CALL xpendn(xloc(i)*0.001 ,yloc(i)*0.001 )
350     CONTINUE
      ENDIF

      ndata = 1
      long(ndata) = x*pi2deg
      lat (ndata) = y*pi2deg
      ilast = iseg

      IF( iseg.eq.-10) GOTO 900
      GOTO 200

  900 CONTINUE

      REWIND nunit

c To draw grid lines

      do i=1,lsample
      do j=1,lsample
        xlat(i)=(x1+(i-1)*(x2-x1)/(lsample-1))*1000.0
        ylon(j)=(y1+(j-1)*(y2-y1)/(lsample-1))*1000.0
      enddo
      enddo

      CALL xxytoll(lsample,lsample,xlat,ylon,rlat,rlon)

      lonmin = rlon(1,1)
      lonmax = rlon(1,1)
      latmin = rlat(1,1)
      latmax = rlat(1,1)
      do i=1,lsample
      do j=1,lsample
        lonmin=min(lonmin, rlon(i,j))
        lonmax=max(lonmax, rlon(i,j))
        latmin=min(latmin, rlat(i,j))
        latmax=max(latmax, rlat(i,j))
      enddo
      enddo

c     print*,'lonmin, latmin, lonmax, latmax=',
c    : lonmin, latmin, lonmax, latmax

      latgrid1 = latgrid
      longrid1 = longrid

      IF( latgrid.eq.0.0 ) latgrid1=5.0
      IF( longrid.eq.0.0 ) longrid1=5.0

      if(lonmin.lt.0.0) lonmin = lonmin-longrid1
      if(latmin.lt.0.0) latmin = latmin-latgrid1
      if(lonmax.gt.0.0) lonmax = lonmax+longrid1
      if(latmax.gt.0.0) latmax = latmax+latgrid1

      lonmin=int(lonmin/longrid1)*longrid1
      lonmax=int(lonmax/longrid1)*longrid1
      latmin=int(latmin/latgrid1)*latgrid1
      latmax=int(latmax/latgrid1)*latgrid1

c     print*,'lonmin, latmin, lonmax, latmax=',
c    : lonmin, latmin, lonmax, latmax

      CALL xbrokn(6,3,6,3)

      IF( latgrid .lt. 0.0) GOTO 650
      jmax=nint((latmax-latmin)/latgrid1)
      DO 600 j=1,jmax+1
        DO 610 i=1,100
          lat (i) = latmin + (j-1)*latgrid1
          long(i) = lonmin + (i-1)/99.0*(lonmax-lonmin)
610     CONTINUE
        CALL xlltoxy( 100,1, lat, long, xloc, yloc)

        CALL xpenup(xloc(1)*0.001 ,yloc(1)*0.001 )
        DO 620 i=2,100
          xloc(i) = xloc(i)*0.001
          yloc(i) = yloc(i)*0.001
          CALL xpendn(xloc(i),yloc(i))
620     CONTINUE
600   CONTINUE
650   CONTINUE

      IF( longrid .lt. 0.0) GOTO 750
      imax=nint((lonmax-lonmin)/longrid1)
      DO 700 i=1,imax+1
        DO 710 j=1,11
          lat (j) = latmin + (j-1)*(latmax-latmin)/10.0
          long(j) = lonmin + (i-1)*longrid1
710     CONTINUE
        CALL xlltoxy( 11,1, lat, long, xloc, yloc)
        CALL xpenup(xloc(1)*0.001 ,yloc(1)*0.001 )
        DO 720 j=2,11
          xloc(j) = xloc(j)*0.001
          yloc(j) = yloc(j)*0.001
          CALL xpendn(xloc(j),yloc(j))
720     CONTINUE
700    CONTINUE
750    CONTINUE

      CALL XFULL

      IF( iwndwon.eq.1) then
        CALL xwindw(xw1,xw2,yw1,yw2)
      else
        CALL xwdwof
      endif

      RETURN
      END

      SUBROUTINE xdrawmap (nunit,mapfile,latgrid,longrid)
        IMPLICIT NONE
        INTEGER,   INTENT(IN) :: nunit
        CHARACTER, INTENT(IN) :: mapfile*(*)
        REAL,      INTENT(IN) :: latgrid,longrid

        !-----------------------------------------------------------

        INTEGER :: gridcol

        !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        CALL xqcolor(gridcol)
        CALL xdrawmap_new(nunit, mapfile, latgrid, longrid, 1, gridcol )

        RETURN
      END SUBROUTINE xdrawmap

      SUBROUTINE xdrawmap_new(nunit, mapfile, latgrid, longrid,         &
     &                        n, gridcol)
!
!-----------------------------------------------------------------------
!     This subroutine will draw a map within a rectagular box in
!     a map projection space. The map projection and plotting space
!     should have been properly set before calling this subroutine.
!-----------------------------------------------------------------------
!
!     nunit     the channel of the mapfile data
!     mapfile   character of map file name
!     latgrid,longrid (degree): the intervals between lat and lon grid lines.
!             < -90.0, no grid lines in the given direction,
!             = 0.0, internally determined,
!             = any real number, typically from 1.0 to 10.0 degrees.
!             < 0.0, grid lines plus axis annotations
!-----------------------------------------------------------------------
!     Author: Ming Xue
!     1/18/199 (M. Xue)
!     Changed to used new format of map data obtained from NCAR
!
!     08/17/2010 (Y. Wang)
!     Added parameter gridcol for controlling map grid colors and
!     Added general option for latitude/longitude axies labels
!     based on former work by Ningzhu Du.
!
!     11/14/2011 (Y. Wang)
!     Added a new parameter n to denote the order in the map drawing.
!     This is to avoid drawing the lat/lon axies labels several times.
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER,   INTENT(IN) :: nunit
      CHARACTER, INTENT(IN) :: mapfile*(*)
      REAL,      INTENT(IN) :: latgrid,longrid
      INTEGER,   INTENT(IN) :: n
      INTEGER,   INTENT(IN) :: gridcol

      integer nmax
      parameter (nmax = 100)
      real xloc(nmax),yloc(nmax),lat(nmax),long(nmax)
      real lonmin,lonmax,latmin,latmax
      real x1,x2,y1,y2,xw1,xw2,yw1,yw2
      integer ndata,i,j
      integer jmax,imax,iwndwon
      integer lsample
      parameter ( lsample = 21)
      real xlat(lsample),ylon(lsample)
     :     ,rlat(lsample,lsample),rlon(lsample,lsample)
      real latgrid1,longrid1
      integer NPTS,IGID
      real XLATMX,XLATMN,XLONMX,XLONMN

      INTEGER   :: label_on
      character :: ch*6, affix*1

      real trulon
      common /pass/trulon

      CALL xqmap (x1,x2,y1,y2)
      CALL xqwdwon( iwndwon )
      CALL xqwndw(xw1,xw2,yw1,yw2)

      CALL xwindw(x1,x2,y1,y2)

      do i=1,lsample
      do j=1,lsample
        xlat(i)=(x1+(i-1)*(x2-x1)/(lsample-1))*1000.0
        ylon(j)=(y1+(j-1)*(y2-y1)/(lsample-1))*1000.0
      enddo
      enddo

      CALL xxytoll(lsample,lsample,xlat,ylon,rlat,rlon)

      lonmin = rlon(1,1)
      lonmax = rlon(1,1)
      latmin = rlat(1,1)
      latmax = rlat(1,1)
      do i=1,lsample
      do j=1,lsample
        lonmin=min(lonmin, rlon(i,j))
        lonmax=max(lonmax, rlon(i,j))
        latmin=min(latmin, rlat(i,j))
        latmax=max(latmax, rlat(i,j))
      enddo
      enddo

! Give some wigger room so that a line whose both ends are outside
! the plotting box but part of the line is inside may get plotted.
      lonmin = lonmin - 0.25  ! subtract 0.25 degree - increase it if some lines are missing
      lonmax = lonmax + 0.25
      latmin = latmin - 0.25
      latmax = latmax + 0.25

      OPEN(nunit,file=mapfile,form='formatted',status='old')
      read (nunit,'(a)')  ! Skip header line 1
      read (nunit,'(a)')  ! Skip header line 2

  200 CONTINUE

   10 READ (nunit,1001,END=900)
     :      NPTS,IGID,XLATMX,XLATMN,XLONMX,XLONMN
c
c igid=1:  CONTINENTAL OUTLINES
C     =2:  US STATE BOUNDARIES (HIGHER RESOLUTION THAN 1)
C     =3:  INTERNATIONAL POLITICAL BOUNDARIES

 1001 FORMAT(2I8,4F8.3)

      IF( npts .lt.2 ) GOTO 200

      ndata = (npts+1)/2
      READ (nunit,1002,END=900)(lat(I),long(I),I=1,ndata)
 1002 FORMAT(10F8.3)

      if( (xlatmn-latmin)*(xlatmn-latmax).lt.0.0 .or.
     :    (xlatmx-latmin)*(xlatmx-latmax).lt.0.0 .or.
     :    (xlonmn-lonmin)*(xlonmn-lonmax).lt.0.0 .or.
     :    (xlonmx-lonmin)*(xlonmx-lonmax).lt.0.0 ) then

      CALL xlltoxy(ndata,1,lat,long,xloc,yloc)
      CALL xpenup(xloc(1)*0.001 ,yloc(1)*0.001 )
      DO 350 i=2,ndata
        CALL xpendn(xloc(i)*0.001 ,yloc(i)*0.001 )
350   CONTINUE

      endif

      GOTO 200

  900 CONTINUE

      CLOSE (nunit)

! To draw grid lines

      label_on = 0

      latgrid1 = latgrid
      longrid1 = longrid
      IF (latgrid1 < 0.0) THEN
        latgrid1 = ABS(latgrid1)
        IF (n == 1) label_on = label_on + 1
      END IF
      IF (longrid1 < 0.0) THEN
        longrid1 = ABS(longrid1)
        IF (n == 1) label_on = label_on + 1
      END IF

      IF( latgrid .eq. 0.0 ) latgrid1=5.0
      IF( longrid .eq. 0.0 ) longrid1=5.0

      if(lonmin.lt.0.0) lonmin = lonmin-longrid1
      if(latmin.lt.0.0) latmin = latmin-latgrid1
      if(lonmax.gt.0.0) lonmax = lonmax+longrid1
      if(latmax.gt.0.0) latmax = latmax+latgrid1

      lonmin=int(lonmin/longrid1)*longrid1
      lonmax=int(lonmax/longrid1)*longrid1
      latmin=int(latmin/latgrid1)*latgrid1
      latmax=int(latmax/latgrid1)*latgrid1

      CALL xbrokn(6,3,6,3)
      CALL xthick(1)
      CALL xchmag(0.015)

      IF( latgrid < -900.0) GOTO 650
      call xcolor(gridcol)
      jmax=nint((latmax-latmin)/latgrid1)
      DO 600 j=1,jmax+1
        DO 610 i=1,100
          lat (i) = latmin + (j-1)*latgrid1
          long(i) = lonmin + (i-1)/99.0*(lonmax-lonmin)
610     CONTINUE
        CALL xlltoxy( 100,1, lat, long, xloc, yloc)

        IF (label_on > 0) THEN
          CALL xwdwof
          IF (lat(1)>=0.0 .AND. lat(1)<=90.0)THEN
            WRITE(affix,'(a)') 'N'
          ELSE
            WRITE(affix,'(a)') 'S'
          END IF
          DO i = 1, 100
            IF ( xloc(i) >= 0 ) then
              !print*,"lat(1),x1",lat(1),x1
              WRITE(ch,'(I2,a)') INT(lat(1)), affix
              IF ( yloc(i)*0.001 >= y1 .AND. yloc(i)*0.001 <= y2 )then
                CALL XCHARR(x1-0.008*(x2-x1),yloc(i)*0.001,CH(1:3))  ! Fectch only the first 3 characters
              END IF
              EXIT
            END IF
          END DO
          CALL xqwdwon( iwndwon )
          CALL xqwndw(xw1,xw2,yw1,yw2)
          CALL xwindw(x1,x2,y1,y2)
        END IF

        CALL xpenup(xloc(1)*0.001 ,yloc(1)*0.001 )
        DO 620 i=2,100
          xloc(i) = xloc(i)*0.001
          yloc(i) = yloc(i)*0.001
          CALL xpendn(xloc(i),yloc(i))
620     CONTINUE
600   CONTINUE
650   CONTINUE

      IF( longrid < -900.0 ) GOTO 750
      CALL xcolor(gridcol)
      imax=nint((lonmax-lonmin)/longrid1)
      DO i=1,imax+1
        DO j=1,11
          lat (j) = latmin + (j-1)*(latmax-latmin)/10.0
          long(j) = lonmin + (i-1)*longrid1
        END DO
        CALL xlltoxy( 11,1, lat, long, xloc, yloc)

        IF (label_on > 0) THEN
          CALL xwdwof
          IF( long(1)>=0.0 .AND. long(1)<=180.0)THEN
            WRITE(affix,'(a)') 'E'
          ELSE
            WRITE(affix,'(a)') 'W'
          END IF

          IF ( xloc(1)*0.001 >= x1 .and. xloc(1)*0.001 <= x2 ) THEN
            WRITE(CH,'(I3,a)') INT(abs(long(1))),affix
            CALL XCHARC(xloc(1)*0.001,Y1-0.04*(Y2-Y1),CH(1:4))
          ENDIF
          CALL xqwdwon( iwndwon )
          CALL xqwndw(xw1,xw2,yw1,yw2)
          CALL xwindw(x1,x2,y1,y2)
        END IF

        CALL xpenup(xloc(1)*0.001 ,yloc(1)*0.001 )
        DO j=2,11
          xloc(j) = xloc(j)*0.001
          yloc(j) = yloc(j)*0.001
          CALL xpendn(xloc(j),yloc(j))
        END DO
      END DO
 750  CONTINUE

      CALL XFULL

!c     IF( iwndwon.eq.1) then
!c       CALL xwindw(xw1,xw2,yw1,yw2)
!c     else
!c       CALL xwdwof
!c     endif

      RETURN
      END

      subroutine xintsy(a)
      return
      end

      subroutine xcontc1(z,x,y,md,m,n,c1,c2)
!
! This routine has the same functionality of xcontx, but
! use a simpler but less efficient algorithm
!
! To do: add missing value skipping capability
!
      dimension z(md,*),x(md,*),y(md,*)
      real xcell(4), ycell(4),zcell(4)

      DO j=1,n-1
      DO i=1,m-1

        xcell(1) = x(i,j)
        xcell(2) = x(i+1,j)
        xcell(3) = x(i+1,j+1)
        xcell(4) = x(i,j+1)
        ycell(1) = y(i,j)
        ycell(2) = y(i+1,j)
        ycell(3) = y(i+1,j+1)
        ycell(4) = y(i,j+1)
        zcell(1) = z(i,j)
        zcell(2) = z(i+1,j)
        zcell(3) = z(i+1,j+1)
        zcell(4) = z(i,j+1)

        call fillcell(xcell, ycell, zcell, c1, c2)

      ENDDO
      ENDDO

      RETURN

      END

      subroutine fillcell(xc,yc,zc, cl1, cl2)
!
! Fill areas between cl1 and cl2 with a specified color
! within a cell defined by four points
!
      implicit none
      real xc(4), yc(4), zc(4), cl1, cl2

      REAL D,p1,p2,b1,b2,cv
      REAL zmin,zmax,x1,x2,y1,y2,z1,z2,cl,z12min,z12max
      REAL xp(20),yp(20)
      INTEGER np,i,i1,i2,no_cl1_found,no_cl2_found
      INTEGER cl1_already_found,cl2_already_found,first_cl,current_cl

      D(P1,P2,B1,B2,CV)=B1+(CV-P1)*(B2-B1)/(P2-P1)

      zmin = min(zc(1),zc(2),zc(3),zc(4))
      zmax = max(zc(1),zc(2),zc(3),zc(4))
      IF( cl1 .gt. zmax ) return
      IF( cl2 .lt. zmin ) return

      IF( zmax.le.cl2.and.zmin.ge.cl1) then
        call xfilarea(xc,yc,4)
        return
      endif

      np = 0
      no_cl1_found = 0
      no_cl2_found = 0
      first_cl = 0
      current_cl = 0

      DO i=1,4

        i1=i
        i2=i+1
        if(i2.gt.4)i2=i2-4

        x1=xc(i1)
        x2=xc(i2)
        y1=yc(i1)
        y2=yc(i2)
        z1=zc(i1)
        z2=zc(i2)

c     write(6,'(a,2i3,6f7.3)')
c    :  'i1,i2,x1,y1,x2,y2,z1,z2=',i1,i2,x1,y1, x2,y2,z1,z2

        z12min = min(z1,z2)
        z12max = max(z1,z2)

!     IF( z12min.gt.cl2) CYCLE
!     IF( z12max.lt.cl1) CYCLE

        IF( z12min.gt.cl2) goto 121
        IF( z12max.lt.cl1) goto 121

        IF( i.eq.1.and.(z1.le.cl2.and.z1.ge.cl1)) then
          np = np + 1
          xp(np)=x1
          yp(np)=y1
          current_cl = 3 ! corner point
c         print*,'np,xp,yp=',np,xp(np),yp(np),first_cl
        endif

        IF( z12max.le.cl2.and.z12min.ge.cl1) then
c         np = np + 1
c         xp(np)=x1
c         yp(np)=y1
          np = np + 1
          xp(np)=x2
          yp(np)=y2
          current_cl = 3 ! corner point
          !CYCLE
          goto 121
        ENDIF

        IF( z2.gt.z1) THEN

          cl = cl1 ! look for cl1 first

          cl2_already_found = 0
!c         IF( (cl-z1)*(cl-z2).lt.0 ) then
          IF( cl.ge.z1.and.cl.le.z2 ) then
            np = np + 1
            xp(np)=d(z1,z2,x1,x2,cl)
            yp(np)=d(z1,z2,y1,y2,cl)
            no_cl1_found = no_cl1_found + 1
            if( first_cl.eq.0) first_cl = 1
            current_cl = 1

            cl = cl2 ! now look for cl2

!c           IF( (cl-z1)*(cl-z2).lt.0 ) then
            IF( cl.ge.z1.and.cl.le.z2 ) then
              np = np + 1
              xp(np)=d(z1,z2,x1,x2,cl)
              yp(np)=d(z1,z2,y1,y2,cl)
              no_cl2_found = no_cl2_found + 1
              cl2_already_found = 1
              if( first_cl.eq.0) first_cl = 2
              current_cl = 2
            else
              np = np + 1
              xp(np)=x2
              yp(np)=y2
              current_cl = 3 ! corner point
            endif
          endif

          cl = cl2 ! now look for cl2
!c         IF( cl2_already_found.eq.0.and.(cl-z1)*(cl-z2).lt.0) then
          IF(cl2_already_found.eq.0.and.(cl.ge.z1.and.cl.le.z2))then
            np = np + 1
            xp(np)=d(z1,z2,x1,x2,cl)
            yp(np)=d(z1,z2,y1,y2,cl)
            no_cl2_found = no_cl2_found + 1
            if( first_cl.eq.0) first_cl = 2
            current_cl = 2
          endif

        else  ! z2.le.z1 case

          cl = cl2 ! look for cl2 first

          cl1_already_found = 0
!c        IF( (cl-z1)*(cl-z2).lt.0 ) then
          IF( cl.ge.z2.and.cl.le.z1 ) then
            np = np + 1
            xp(np)=d(z1,z2,x1,x2,cl)
            yp(np)=d(z1,z2,y1,y2,cl)
            no_cl2_found = no_cl2_found + 1
            if( first_cl.eq.0) first_cl = 2
            current_cl = 2

            cl = cl1 ! now look for cl1
!c          IF( (cl-z1)*(cl-z2).lt.0 ) then
            IF( cl.ge.z2.and.cl.le.z1 ) then
              np = np + 1
              xp(np)=d(z1,z2,x1,x2,cl)
              yp(np)=d(z1,z2,y1,y2,cl)
              cl1_already_found = 1
              no_cl1_found = no_cl1_found + 1
              if( first_cl.eq.0) first_cl = 1
              current_cl = 1
            else
              np = np + 1
              xp(np)=x2
              yp(np)=y2
              current_cl = 3 ! corner point
            endif
          endif

          cl = cl1 ! now look for cl1
!c        IF(cl1_already_found.eq.0.and.(cl-z1)*(cl-z2).lt.0) then
          IF(cl1_already_found.eq.0.and.(cl.ge.z2.and.cl.le.z1))then
            np = np + 1
            xp(np)=d(z1,z2,x1,x2,cl)
            yp(np)=d(z1,z2,y1,y2,cl)
            no_cl1_found = no_cl1_found + 1
            if( first_cl.eq.0) first_cl = 1
            current_cl = 1
          endif

        endif


        IF( no_cl2_found .eq.2 .or. no_cl1_found.eq.2 ) then
          IF( current_cl .eq. first_cl ) then
            call xfilarea(xp,yp,np)
            no_cl1_found = 0
            no_cl2_found = 0
            first_cl=0
            current_cl=0
            np = 0
          ENDIF
        ENDIF

121     CONTINUE
      END DO ! loop over four sides

111   continue
      IF(np.ne.0) then
        call xfilarea(xp,yp,np)
      endif

      return
      end subroutine fillcell

      subroutine xcontcopt(kcontcopt)
      integer icontcopt
      common /xcontc_opt/ icontcopt

      icontcopt = kcontcopt

      return
      end subroutine xcontcopt
