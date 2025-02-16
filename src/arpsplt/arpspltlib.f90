!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE INITPLTPARA                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE initpltpara(nx,ny,nz,nzsoil,nstyps,outfilename,iwtype,istatus)

!
!-----------------------------------------------------------------------
!
!  This is the subroutine to initilize ARPSPLT parameter from the
!  namelist input file arpsplt.input
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang, CAPS/OU.
!    12/17/2002.
!
!  MODIFICATION HISTORY:
!  05/08/2012 (Y. Wang)
!  Added capability to read command line for namelist file.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz          ! Grid dimensions.
  INTEGER :: nzsoil            ! levels of soil model
  INTEGER :: nstyps            ! Maximum number of soil types.

  INTEGER, INTENT(OUT) :: istatus

  INTEGER, PARAMETER :: nhisfile_max=200
  INTEGER, PARAMETER :: max_dim = 200
  INTEGER, PARAMETER :: fzone = 3
  INTEGER, PARAMETER :: nscalarmax=30
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'indtflg.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'arpsplt.inc'
  INCLUDE 'alloc.inc'

  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
! Variables for mpi jobs
!
!-----------------------------------------------------------------------
  INTEGER :: nprocx_in, nprocy_in   ! number of processors in input data
  INTEGER :: ncompressx, ncompressy ! compression in x and y direction:
                                    ! ncompressx=nprocx_in/nproc_x
                                    ! ncompressy=nprocy_in/nproc_y
  INTEGER :: nproc_node
  INTEGER :: nprocx_lw, nprocy_lw

  NAMELIST /message_passing/ nproc_x, nproc_y, max_fopen, nproc_node,   &
                             nprocx_in, nprocy_in,nprocx_lw,nprocy_lw
  COMMON /init1_mpi/ ncompressx, ncompressy, nproc_node, nprocx_lw, nprocy_lw

!-----------------------------------------------------------------------
!
! Variables in NAMELISTs
!
!-----------------------------------------------------------------------

  INTEGER :: hinfmt, nhisfile
  CHARACTER (LEN=256) :: grdbasfn
  CHARACTER (LEN=256) :: hisfile(nhisfile_max)
                              ! base name, DONOT contain processor info.
  CHARACTER (LEN=256) :: hdmpftrailer
  COMMON /init2_hisf/ hinfmt,nhisfile, grdbasfn, hisfile,hdmpftrailer

  INTEGER :: layout,nxpic,nypic,inwfrm
  REAL    :: paprlnth
  NAMELIST /page_setup/ layout, nxpic, nypic, inwfrm,paprlnth
  COMMON /init3_page/ layout, inwfrm, paprlnth

  INTEGER :: iorig
  REAL    :: xbgn,xend,ybgn,yend,zbgn,zend,zsoilbgn,zsoilend
  REAL    :: yxstrch        ! Stretching factor for x-y plots
  REAL    :: zxstrch        ! Stretching factor for x-z plots
  REAL    :: zystrch        ! Stretching factor for y-z plots
  REAL    :: zhstrch        ! Stretching factor for arbitrary vertical slices
  REAL    :: zsoilxstrch    ! Stretching factor for x-z plots for the soil model
  REAL    :: zsoilystrch    ! Stretching factor for y-z plots for the soil model
  REAL    :: winsiz         ! A global factor for window size
  REAL    :: margnx, margny ! margin
  INTEGER :: pcolbar        ! position of color bar
  INTEGER :: iskip, jskip
  NAMELIST /plotting_setup/ iorig, xorig, yorig,                        &
      xbgn, xend, ybgn, yend, zbgn, zend, zsoilbgn, zsoilend,           &
      yxstrch, zxstrch, zystrch, zhstrch, zsoilxstrch, zsoilystrch,     &
      winsiz, margnx, margny,pcolbar,iskip,jskip
  COMMON /init4_plotset/ iorig,zbgn,zend,zsoilbgn,zsoilend,             &
                       yxstrch,zxstrch,zystrch,                         &
                       zhstrch,zsoilxstrch,zsoilystrch,                 &
                       margnx,margny
  COMMON /pltwdw/ xbgn,xend,ybgn,yend,iskip,jskip
  COMMON /windows/ winsiz

  INTEGER :: col_table
  CHARACTER (LEN=256) :: color_map
  NAMELIST /col_table_cntl/ col_table,color_map
  COMMON /init5_coltab/ color_map
  COMMON /coltable/col_table,pcolbar

  INTEGER :: lnmag,fontopt,lbaxis,axlbfmt
  INTEGER :: haxisu, vaxisu
  INTEGER :: tickopt
  INTEGER :: presaxis_no
  INTEGER :: ctrlbopt, ctrstyle, ctrlbfrq
  INTEGER :: lbmaskopt
  REAL    :: lblmag   ! A global magnification factor for labels.
  REAL    :: ctrlbsiz, axlbsiz
  REAL    :: hmintick,vmajtick,vmintick,hmajtick
  REAL    :: pres_val(20), pres_z(20)
  NAMELIST /style_tuning/ lblmag,lnmag, fontopt,                        &
      lbaxis,axlbfmt,axlbsiz, haxisu, vaxisu,                           &
      tickopt,hmintick,vmajtick,vmintick,hmajtick,                      &
      presaxis_no,pres_val,                                             &
      ctrlbopt,ctrstyle,ctrlbfrq,ctrlbsiz,lbmaskopt
  COMMON /init6_style/ lnmag, ctrlbopt, ctrstyle
  COMMON /labmag/ lblmag, ctrlbsiz, axlbsiz
  COMMON /var_par/ fontopt,haxisu, vaxisu,lbaxis,tickopt,               &
          hmintick,vmajtick,vmintick,hmajtick,axlbfmt
  COMMON /pressbar_par/presaxis_no,pres_val,pres_z
  COMMON /clb_frq/ ctrlbfrq

  INTEGER :: smooth
  NAMELIST /smooth_cntl/ smooth
  COMMON /smoothopt/smooth

  INTEGER :: ntitle,titcol, wpltime
  REAL :: titsiz
  CHARACTER (LEN=256) :: title(3), footer_l, footer_c, footer_r
  NAMELIST /title_setup/ntitle,titcol,titsiz,title
  NAMELIST /footer_setup/wpltime,footer_l,footer_c, footer_r
  COMMON /titpar1/title, footer_l, footer_c, footer_r
  COMMON /titpar2/ntitle,titcol,wpltime, nxpic, nypic
  COMMON /titpar3/titsiz

  INTEGER :: ovrmap,mapgrid,mapgridcol,nmapfile,mapcol(maxmap),         &
             mapline_style(maxmap)
  REAL :: latgrid,longrid
  CHARACTER (LEN=256) :: mapfile(maxmap)
  NAMELIST /map_plot/ ovrmap,mapgrid,latgrid,longrid,mapgridcol,        &
      nmapfile,mapfile,mapcol,mapline_style
  COMMON /mappar / ovrmap
  COMMON /mappar1/ nmapfile,mapcol,mapline_style,mapfile
  COMMON /mappar2/ mapgrid,mapgridcol, latgrid,longrid

  INTEGER :: missfill_opt,missval_colind    ! miss value color index
  NAMELIST /multi_setup/ missfill_opt,missval_colind
  COMMON /multi_value/ missfill_opt, missval_colind

  INTEGER :: nslice_xy, slice_xy(max_dim)
  INTEGER :: nslice_xz, slice_xz(max_dim)
  INTEGER :: nslice_yz, slice_yz(max_dim)
  INTEGER :: nslice_h, nslice_p, nslice_pt, nslice_v
  REAL :: slice_h(max_dim), slice_p(max_dim),  slice_pt(max_dim)
  REAL :: xpnt1(max_dim),ypnt1(max_dim),xpnt2(max_dim),ypnt2(max_dim)

  INTEGER :: nslice_xy_soil, slice_xy_soil(max_dim)
  INTEGER :: nslice_xz_soil, slice_xz_soil(max_dim)
  INTEGER :: nslice_yz_soil, slice_yz_soil(max_dim)

  NAMELIST /xy_slice_cntl/ nslice_xy, slice_xy
  NAMELIST /xz_slice_cntl/ nslice_xz, slice_xz
  NAMELIST /yz_slice_cntl/ nslice_yz, slice_yz
  NAMELIST /h_slice_cntl/  nslice_h, slice_h
  NAMELIST /v_slice_cntl/  nslice_v, xpnt1,ypnt1,xpnt2,ypnt2
  NAMELIST /p_slice_cntl/  nslice_p, slice_p
  NAMELIST /pt_slice_cntl/ nslice_pt, slice_pt
  NAMELIST /xy_soil_slice_cntl/ nslice_xy_soil, slice_xy_soil
  NAMELIST /xz_soil_slice_cntl/ nslice_xz_soil, slice_xz_soil
  NAMELIST /yz_soil_slice_cntl/ nslice_yz_soil, slice_yz_soil
  COMMON /init7_slice/ nslice_xy, nslice_xz, nslice_yz, nslice_h,        &
                       nslice_v,  nslice_p,  nslice_pt,                  &
                       nslice_xy_soil, nslice_xz_soil, nslice_yz_soil,   &
                       slice_xy, slice_xz, slice_yz, slice_h,            &
                       slice_p,  slice_pt, xpnt1, ypnt1, xpnt2, ypnt2,   &
                       slice_xy_soil, slice_xz_soil, slice_yz_soil,      &
                       imove

  INTEGER :: imove
  NAMELIST /domain_move/ imove, umove, vmove

!-----------------------------------------------------------------------
!
!  *inc        -- Contour intervals
!  *minc,*maxc -- Limited variable minimum and maximum for color
!                 contour shade
!  *ovr        -- Overlay control parameters
!  *hlf        -- highlighting frequency for contour parameters
!  *zro        -- define the attributes of zero contour to be plotted
!                 parameters
!  *sty        -- Define the option for contour line stypes.
!
!-----------------------------------------------------------------------

  INTEGER :: hplot, msfplt, thkplt, tplot, uplot, vplot, vhplot, vsplot,&
             wplot, ptplot, pplot,  ipvplt
  REAL :: hinc, msfinc, thkinc, tinc, uinc, vinc, vhinc, vsinc,         &
          winc, ptinc, pinc, ipvinc
  REAL :: hminc, msfminc, thkminc, tminc, uminc, vminc, vhminc, vsminc, &
          wminc, ptminc,  pminc, ipvminc
  REAL :: hmaxc, msfmaxc, thkmaxc, tmaxc, umaxc, vmaxc, vhmaxc, vsmaxc, &
          wmaxc, ptmaxc,  pmaxc, ipvmaxc
  INTEGER :: hovr, msfovr, thkovr, tovr, uovr, vovr, vhovr, vsovr,      &
             wovr , ptovr , povr, ipvovr
  INTEGER :: hcol1, msfcol1, thkcol1, tcol1, ucol1, vcol1, vhcol1, vscol1,      &
             wcol1 , ptcol1 , pcol1, ipvcol1
  INTEGER :: hcol2, msfcol2, thkcol2, tcol2, ucol2, vcol2, vhcol2, vscol2,      &
             wcol2 , ptcol2 , pcol2, ipvcol2
  INTEGER :: hprio, msfprio, thkprio, tprio, uprio, vprio, vhprio, vsprio,      &
             wprio , ptprio , pprio, ipvprio
  INTEGER :: hhlf, msfhlf, thkhlf, thlf, uhlf, vhlf, vhhlf, vshlf,      &
             whlf , pthlf , phlf, ipvhlf
  INTEGER :: hzro, msfzro, thkzro, tzro, uzro, vzro, vhzro, vszro,      &
             wzro , ptzro , pzro, ipvzro
  INTEGER :: hsty, msfsty, thksty, tsty, usty, vsty, vhsty, vssty,      &
             wsty , ptsty , psty, ipvsty
  CHARACTER (LEN=1) :: tunits  ! units for temperature F or C
  INTEGER :: vhunits

  NAMELIST /sclrplt_cntl1/                                              &
      hplot,  hinc,   hminc,   hmaxc,  hovr,   hcol1,hcol2, hprio,      &
            hhlf, hzro, hsty,                                           &
      msfplt,msfinc,msfminc, msfmaxc,msfovr,msfcol1,msfcol2,msfprio,    &
            msfhlf, msfzro,  msfsty,                                    &
      thkplt,thkinc,thkminc, thkmaxc,thkovr,thkcol1,thkcol2,thkprio,    &
            thkhlf, thkzro,  thksty,                                    &
      tplot,  tinc,   tminc,   tmaxc,  tovr,   tcol1,tcol2, tprio,      &
            tunits, thlf, tzro, tsty,                                   &
      uplot,  uinc,   uminc,   umaxc,  uovr,   ucol1,ucol2, uprio,      &
            uhlf, uzro, usty,                                           &
      vplot,  vinc,   vminc,   vmaxc,  vovr,   vcol1,vcol2, vprio,      &
            vhlf, vzro, vsty,                                           &
      vhplot, vhinc,  vhminc,  vhmaxc, vhovr,  vhcol1,vhcol2,vhprio,    &
            vhunits, vhhlf, vhzro, vhsty,                               &
      vsplot, vsinc,  vsminc,  vsmaxc, vsovr,  vscol1,vscol2,vsprio,    &
            vshlf, vszro, vssty,                                        &
      wplot,  winc,   wminc,   wmaxc,  wovr,   wcol1,wcol2, wprio,      &
            whlf, wzro, wsty,                                           &
      ptplot, ptinc,  ptminc,  ptmaxc, ptovr,  ptcol1,ptcol2,ptprio,    &
            pthlf, ptzro, ptsty,                                        &
      pplot , pinc,   pminc,   pmaxc,  povr,   pcol1,pcol2, pprio,      &
            phlf, pzro, psty,                                           &
      ipvplt,ipvinc,ipvminc, ipvmaxc,ipvovr,ipvcol1,ipvcol2,ipvprio,    &
            ipvhlf, ipvzro,  ipvsty

  COMMON /init8_cntl1/                 &
      hplot,  hinc,   hminc,   hmaxc,  hovr,   hcol1,hcol2, hprio,      &
            hhlf, hzro, hsty,                                           &
      msfplt,msfinc,msfminc, msfmaxc,msfovr,msfcol1,msfcol2,msfprio,    &
            msfhlf, msfzro,  msfsty,                                    &
      thkplt,thkinc,thkminc, thkmaxc,thkovr,thkcol1,thkcol2,thkprio,    &
            thkhlf, thkzro,  thksty,                                    &
      tplot,  tinc,   tminc,   tmaxc,  tovr,   tcol1,tcol2, tprio,      &
            thlf, tzro, tsty,                                   &
      uplot,  uinc,   uminc,   umaxc,  uovr,   ucol1,ucol2, uprio,      &
            uhlf, uzro, usty,                                           &
      vplot,  vinc,   vminc,   vmaxc,  vovr,   vcol1,vcol2, vprio,      &
            vhlf, vzro, vsty,                                           &
      vhplot, vhinc,  vhminc,  vhmaxc, vhovr,  vhcol1,vhcol2,vhprio,    &
            vhunits, vhhlf, vhzro, vhsty,                               &
      vsplot, vsinc,  vsminc,  vsmaxc, vsovr,  vscol1,vscol2,vsprio,    &
            vshlf, vszro, vssty,                                        &
      wplot,  winc,   wminc,   wmaxc,  wovr,   wcol1,wcol2, wprio,      &
            whlf, wzro, wsty,                                           &
      ptplot, ptinc,  ptminc,  ptmaxc, ptovr,  ptcol1,ptcol2,ptprio,    &
            pthlf, ptzro, ptsty,                                        &
      pplot , pinc,   pminc,   pmaxc,  povr,   pcol1,pcol2, pprio,      &
            phlf, pzro, psty,                                           &
      ipvplt,ipvinc,ipvminc, ipvmaxc,ipvovr,ipvcol1,ipvcol2,ipvprio,    &
            ipvhlf, ipvzro,  ipvsty

  INTEGER :: qvplot,qscalarplot(nscalarmax),qwplot,qtplot
  REAL :: qvinc,qscalarinc(nscalarmax),qwinc,qtinc
  REAL :: qvminc,qscalarminc(nscalarmax),qwminc,qtminc
  REAL :: qvmaxc,qscalarmaxc(nscalarmax),qwmaxc,qtmaxc
  INTEGER :: qvovr,qscalarovr(nscalarmax),qwovr,qtovr
  INTEGER :: qvcol1,qscalarcol1(nscalarmax),qwcol1,qtcol1
  INTEGER :: qvcol2,qscalarcol2(nscalarmax),qwcol2,qtcol2
  INTEGER :: qvprio,qscalarprio(nscalarmax),qwprio,qtprio
  INTEGER :: qvhlf,qscalarhlf(nscalarmax),qwhlf,qthlf
  INTEGER :: qvzro,qscalarzro(nscalarmax),qwzro,qtzro
  INTEGER :: qvsty,qscalarsty(nscalarmax),qwsty,qtsty

  NAMELIST /sclrplt_cntl2/                                              &
      qvplot, qvinc,  qvminc,  qvmaxc, qvovr,  qvcol1,qvcol2,qvprio,    &
            qvhlf, qvzro, qvsty,                                        &
      qscalarplot,qscalarinc,qscalarminc,qscalarmaxc,qscalarovr,        &
            qscalarcol1,qscalarcol2,qscalarprio,qscalarhlf,             &
            qscalarzro,qscalarsty,                                      &
      qwplot, qwinc,  qwminc,  qwmaxc, qwovr,  qwcol1,qwcol2,qwprio,    &
            qwhlf, qwzro, qwsty,                                        &
      qtplot, qtinc,  qtminc,  qtmaxc, qtovr,  qtcol1,qtcol2,qtprio,    &
            qthlf, qtzro, qtsty
  COMMON /init9_cntl2/            &
      qvplot, qvinc,  qvminc,  qvmaxc, qvovr,  qvcol1,qvcol2,qvprio,    &
            qvhlf, qvzro, qvsty,                                        &
      qscalarplot,qscalarinc,qscalarminc,qscalarmaxc,qscalarovr,        &
            qscalarcol1,qscalarcol2,qscalarprio,qscalarhlf,             &
            qscalarzro,qscalarsty,                                      &
      qwplot, qwinc,  qwminc,  qwmaxc, qwovr,  qwcol1,qwcol2,qwprio,    &
            qwhlf, qwzro, qwsty,                                        &
      qtplot, qtinc,  qtminc,  qtmaxc, qtovr,  qtcol1,qtcol2,qtprio,    &
            qthlf, qtzro, qtsty

  INTEGER :: kmhplt,kmvplt,tkeplt,rhplot,tdplot,rfplot,rfcplt,pteplt,   &
             zdrplt,kdpplt,zdpplt,rhvplt
  INTEGER :: rfopt, dualpol
  CHARACTER (LEN=256) :: rsadir
  REAL :: wavelen
  REAL :: kmhinc,kmvinc,tkeinc,rhinc,tdinc,rfinc,rfcinc,pteinc,         &
          zdrinc,kdpinc,zdpinc,rhvinc
  REAL :: kmhminc,kmvminc,tkeminc,rhminc,tdminc,rfminc,rfcminc,pteminc,&
          zdrminc,kdpminc,zdpminc,rhvminc
  REAL :: kmhmaxc,kmvmaxc,tkemaxc,rhmaxc,tdmaxc,rfmaxc,rfcmaxc,ptemaxc, &
          zdrmaxc,kdpmaxc,zdpmaxc,rhvmaxc
  INTEGER :: kmhovr,kmvovr,tkeovr,rhovr,tdovr,rfovr,rfcovr,pteovr,      &
             zdrovr,kdpovr,zdpovr,rhvovr
  INTEGER :: kmhcol1,kmvcol1,tkecol1,rhcol1,tdcol1,rfcol1,rfccol1,ptecol1,&
             zdrcol1,kdpcol1,zdpcol1,rhvcol1
  INTEGER :: kmhcol2,kmvcol2,tkecol2,rhcol2,tdcol2,rfcol2,rfccol2,ptecol2,&
             zdrcol2,kdpcol2,zdpcol2,rhvcol2
  INTEGER :: kmhprio,kmvprio,tkeprio,rhprio,tdprio,rfprio,rfcprio,pteprio,&
             zdrprio,kdpprio,zdpprio,rhvprio
  INTEGER :: kmhhlf,kmvhlf,tkehlf,rhhlf,tdhlf,rfhlf,rfchlf,ptehlf,      &
             zdrhlf,kdphlf,zdphlf,rhvhlf
  INTEGER :: kmhzro,kmvzro,tkezro,rhzro,tdzro,rfzro,rfczro,ptezro,      &
             zdrzro,kdpzro,zdpzro,rhvzro
  INTEGER :: kmhsty,kmvsty,tkesty,rhsty,tdsty,rfsty,rfcsty,ptesty,      &
             zdrsty,kdpsty,zdpsty,rhvsty
  CHARACTER (LEN=1) :: tdunits ! units for dew-point temp F or C

  NAMELIST /sclrplt_cntl3/                                              &
      kmhplt, kmhinc, kmhminc,kmhmaxc,kmhovr,kmhcol1,kmhcol2,kmhprio,   &
            kmhhlf, kmhzro, kmhsty,                                     &
      kmvplt, kmvinc, kmvminc,kmvmaxc,kmvovr,kmvcol1,kmvcol2,kmvprio,   &
            kmvhlf, kmvzro, kmvsty,                                     &
      tkeplt, tkeinc, tkeminc, tkemaxc,tkeovr,tkecol1,tkecol2,tkeprio,  &
            tkehlf, tkezro, tkesty,                                     &
      rhplot, rhinc,  rhminc,  rhmaxc,  rhovr,   rhcol1,rhcol2,rhprio,  &
            rhhlf, rhzro,   rhsty,                                      &
      tdplot, tdinc,  tdminc,  tdmaxc,  tdovr,   tdcol1,tdcol2,tdprio,  &
            tdunits, tdhlf, tdzro, tdsty,                               &
      dualpol,rfopt,rsadir,wavelen,graupel_ON,hail_ON,                  &
      rfplot, rfinc,  rfminc,  rfmaxc,  rfovr,   rfcol1,rfcol2,rfprio,  &
            rfhlf, rfzro, rfsty,                                        &
      rfcplt, rfcinc, rfcminc, rfcmaxc,rfcovr,rfccol1,rfccol2,rfcprio,  &
            rfchlf, rfczro, rfcsty,                                     &
      pteplt, pteinc, pteminc, ptemaxc,pteovr,ptecol1,ptecol2,pteprio,  &
            ptehlf, ptezro, ptesty,                                     &
      zdrplt, zdrinc, zdrminc, zdrmaxc,zdrovr,zdrcol1,zdrcol2,zdrprio,  &
            zdrhlf, zdrzro, zdrsty,                                     &
      kdpplt, kdpinc, kdpminc, kdpmaxc,kdpovr,kdpcol1,kdpcol2,kdpprio,  &
            kdphlf, kdpzro, kdpsty,                                     &
      zdpplt, zdpinc, zdpminc, zdpmaxc,zdpovr,zdpcol1,zdpcol2,zdpprio,  &
            zdphlf, zdpzro, zdpsty,                                     &
      rhvplt, rhvinc, rhvminc, rhvmaxc,rhvovr,rhvcol1,rhvcol2,rhvprio,  &
            rhvhlf, rhvzro, rhvsty

  COMMON /init10_cntl3/                              &
      kmhplt, kmhinc, kmhminc,kmhmaxc,kmhovr,kmhcol1,kmhcol2,kmhprio,   &
            kmhhlf, kmhzro, kmhsty,                                     &
      kmvplt, kmvinc, kmvminc,kmvmaxc,kmvovr,kmvcol1,kmvcol2,kmvprio,   &
            kmvhlf, kmvzro, kmvsty,                                     &
      tkeplt, tkeinc, tkeminc, tkemaxc,tkeovr,tkecol1,tkecol2,tkeprio,  &
            tkehlf, tkezro, tkesty,                                     &
      rhplot, rhinc,  rhminc,  rhmaxc,  rhovr,   rhcol1,rhcol2,rhprio,  &
            rhhlf, rhzro,   rhsty,                                      &
      tdplot, tdinc,  tdminc,  tdmaxc,  tdovr,   tdcol1,tdcol2,tdprio,  &
            tdhlf, tdzro, tdsty,                                        &
      dualpol, rfopt, rsadir, wavelen,                                  &
      rfplot, rfinc,  rfminc,  rfmaxc,  rfovr,   rfcol1,rfcol2,rfprio,  &
            rfhlf, rfzro, rfsty,                                        &
      rfcplt, rfcinc, rfcminc, rfcmaxc,rfcovr,rfccol1,rfccol2,rfcprio,  &
            rfchlf, rfczro, rfcsty,                                     &
      pteplt, pteinc, pteminc, ptemaxc,pteovr,ptecol1,ptecol2,pteprio,  &
            ptehlf, ptezro, ptesty,                                     &
      zdrplt, zdrinc, zdrminc, zdrmaxc,zdrovr,zdrcol1,zdrcol2,zdrprio,  &
            zdrhlf, zdrzro, zdrsty,                                     &
      kdpplt, kdpinc, kdpminc, kdpmaxc,kdpovr,kdpcol1,kdpcol2,kdpprio,  &
            kdphlf, kdpzro, kdpsty,                                     &
      zdpplt, zdpinc, zdpminc, zdpmaxc,zdpovr,zdpcol1,zdpcol2,zdpprio,  &
            zdphlf, zdpzro, zdpsty,                                     &
      rhvplt, rhvinc, rhvminc, rhvmaxc,rhvovr,rhvcol1,rhvcol2,rhvprio,  &
            rhvhlf, rhvzro, rhvsty

  COMMON /init810_char_units/ tunits, tdunits

  INTEGER :: upplot,vpplot,wpplot,ptpplt,ppplot,qvpplt,                 &
             vorpplt, divpplt, divqplt
  REAL :: upinc,vpinc,wpinc,ptpinc,ppinc,qvpinc,                        &
             vorpinc, divpinc, divqinc
  REAL :: upminc,vpminc,wpminc,ptpminc,ppminc,qvpminc,                  &
             vorpminc, divpminc, divqminc
  REAL :: upmaxc,vpmaxc,wpmaxc,ptpmaxc,ppmaxc,qvpmaxc,                  &
             vorpmaxc, divpmaxc, divqmaxc
  INTEGER :: upovr,vpovr,wpovr,ptpovr,ppovr,qvpovr,                     &
             vorpovr, divpovr, divqovr
  INTEGER :: upcol1,vpcol1,wpcol1,ptpcol1,ppcol1,qvpcol1,               &
             vorpcol1, divpcol1, divqcol1
  INTEGER :: upcol2,vpcol2,wpcol2,ptpcol2,ppcol2,qvpcol2,               &
             vorpcol2, divpcol2, divqcol2
  INTEGER :: upprio,vpprio,wpprio,ptpprio,ppprio,qvpprio,               &
             vorpprio, divpprio, divqprio
  INTEGER :: uphlf,vphlf,wphlf,ptphlf,pphlf,qvphlf,                     &
             vorphlf, divphlf, divqhlf
  INTEGER :: upzro,vpzro,wpzro,ptpzro,ppzro,qvpzro,                     &
             vorpzro, divpzro, divqzro
  INTEGER :: upsty,vpsty,wpsty,ptpsty,ppsty,qvpsty,                     &
             vorpsty, divpsty, divqsty

  NAMELIST /sclrplt_cntl_prt1/                                          &
      upplot, upinc,  upminc,   upmaxc,   upovr,upcol1,upcol2,upprio,   &
            uphlf, upzro, upsty,                                        &
      vpplot, vpinc,  vpminc,   vpmaxc,   vpovr,vpcol1,vpcol2,vpprio,   &
            vphlf, vpzro, vpsty,                                        &
      wpplot, wpinc,  wpminc,   wpmaxc,   wpovr,wpcol1,wpcol2,wpprio,   &
            wphlf, wpzro, wpsty,                                        &
      ptpplt, ptpinc, ptpminc,ptpmaxc,ptpovr,ptpcol1,ptpcol2,ptpprio,   &
            ptphlf, ptpzro, ptpsty,                                     &
      ppplot, ppinc,  ppminc, ppmaxc,  ppovr,   ppcol1,ppcol2,ppprio,   &
            pphlf, ppzro, ppsty,                                        &
      qvpplt, qvpinc, qvpminc,qvpmaxc,qvpovr,qvpcol1,qvpcol2,qvpprio,   &
            qvphlf, qvpzro, qvpsty,                                     &
      vorpplt,vorpinc,vorpminc, vorpmaxc, vorpovr, vorpcol1,vorpcol2,   &
            vorphlf,  vorpprio, vorpzro, vorpsty,                       &
      divpplt,divpinc,divpminc, divpmaxc, divpovr, divpcol1,divpcol2,   &
            divphlf,  divpprio, divpzro, divpsty,                       &
      divqplt,divqinc,divqminc, divqmaxc, divqovr, divqcol1,divqcol2,   &
            divqhlf,divqprio, divqzro,divqsty
  COMMON /init11_cntl_prt1/               &
      upplot, upinc,  upminc,   upmaxc,   upovr,upcol1,upcol2,upprio,   &
            uphlf, upzro, upsty,                                        &
      vpplot, vpinc,  vpminc,   vpmaxc,   vpovr,vpcol1,vpcol2,vpprio,   &
            vphlf, vpzro, vpsty,                                        &
      wpplot, wpinc,  wpminc,   wpmaxc,   wpovr,wpcol1,wpcol2,wpprio,   &
            wphlf, wpzro, wpsty,                                        &
      ptpplt, ptpinc, ptpminc,ptpmaxc,ptpovr,ptpcol1,ptpcol2,ptpprio,   &
            ptphlf, ptpzro, ptpsty,                                     &
      ppplot, ppinc,  ppminc, ppmaxc,  ppovr,   ppcol1,ppcol2,ppprio,   &
            pphlf, ppzro, ppsty,                                        &
      qvpplt, qvpinc, qvpminc,qvpmaxc,qvpovr,qvpcol1,qvpcol2,qvpprio,   &
            qvphlf, qvpzro, qvpsty,                                     &
      vorpplt,vorpinc,vorpminc, vorpmaxc, vorpovr, vorpcol1,vorpcol2,   &
            vorphlf,  vorpprio, vorpzro, vorpsty,                       &
      divpplt,divpinc,divpminc, divpmaxc, divpovr, divpcol1,divpcol2,   &
            divphlf,  divpprio, divpzro, divpsty,                       &
      divqplt,divqinc,divqminc, divqmaxc, divqovr, divqcol1,divqcol2,   &
            divqhlf,divqprio, divqzro,divqsty

  INTEGER :: gricplt, avorplt, rhiplot
  REAl :: gricinc, avorinc, rhiinc
  REAl :: gricminc, avorminc, rhiminc
  REAl :: gricmaxc, avormaxc, rhimaxc
  INTEGER :: gricovr, avorovr, rhiovr
  INTEGER :: griccol1, avorcol1, rhicol1
  INTEGER :: griccol2, avorcol2, rhicol2
  INTEGER :: gricprio, avorprio, rhiprio
  INTEGER :: grichlf, avorhlf, rhihlf
  INTEGER :: griczro, avorzro, rhizro
  INTEGER :: gricsty, avorsty, rhisty
  NAMELIST /sclrplt_cntl_prt2/                                          &
      gricplt,gricinc,gricminc, gricmaxc, gricovr, griccol1,griccol2,   &
            grichlf,gricprio, griczro, gricsty,                         &
      avorplt,avorinc,avorminc, avormaxc, avorovr, avorcol1,avorcol2,   &
            avorhlf,avorprio, avorzro, avorsty,                         &
      rhiplot, rhiinc,  rhiminc,  rhimaxc,  rhiovr,   rhicol1,rhicol2,  &
            rhiprio,  rhihlf, rhizro , rhisty
  COMMON /init12_cntl_prt2/          &
      gricplt,gricinc,gricminc, gricmaxc, gricovr, griccol1,griccol2,   &
            grichlf,gricprio, griczro, gricsty,                         &
      avorplt,avorinc,avorminc, avormaxc, avorovr, avorcol1,avorcol2,   &
            avorhlf,avorprio, avorzro, avorsty,                         &
      rhiplot, rhiinc,  rhiminc,  rhimaxc,  rhiovr,   rhicol1,rhicol2,  &
            rhiprio,  rhihlf, rhizro , rhisty

  INTEGER :: istride,jstride,kstride
  INTEGER :: vtrplt, vtpplt, xuvplt, strmplt, vagplt
  REAL :: vtrunit, vtpunit, xuvunit, strmunit, vagunit
  INTEGER :: vtrovr, vtpovr, xuvovr, strmovr, vagovr
  INTEGER :: vtrcol1, vtpcol1, xuvcol1, strmcol1, vagcol1
  INTEGER :: vtrcol2, vtpcol2, xuvcol2, strmcol2, vagcol2
  INTEGER :: vtrprio, vtpprio, xuvprio, strmprio, vagprio
  INTEGER :: vtrunits, vtpunits, xuvunits, strmunits, vagunits
  INTEGER :: vtrtype, vtptype, xuvtype, strmtype, vagtype

  NAMELIST /vctrplt_cntl/istride,jstride,kstride,                       &
      vtrplt, vtrunit,vtrovr,vtrcol1,vtrcol2,vtrprio,vtrunits,vtrtype,  &
      vtpplt, vtpunit, vtpovr,vtpcol1,vtpcol2,vtpprio,vtpunits,vtptype,  &
      xuvplt, xuvunit,xuvovr,xuvcol1,xuvcol2,xuvprio,xuvunits,xuvtype,  &
      strmplt,strmunit,strmovr,strmcol1,strmcol2,strmprio,strmunits,    &
      strmtype,                                                         &
      vagplt, vagunit,vagovr,vagcol1,vagcol2,vagprio,vagunits,vagtype

  COMMON /init13_cntl_vctr/ istride,jstride,kstride,                    &
      vtrplt, vtrunit,vtrovr,vtrcol1,vtrcol2,vtrprio,vtrunits,vtrtype,  &
      vtpplt, vtpunit, vtpovr,vtpcol1,vtpcol2,vtpprio,vtpunits,vtptype,  &
      xuvplt, xuvunit,xuvovr,xuvcol1,xuvcol2,xuvprio,xuvunits,xuvtype,  &
      strmplt,strmunit,strmovr,strmcol1,strmcol2,strmprio,strmunits,    &
      strmtype,                                                         &
      vagplt, vagunit,vagovr,vagcol1,vagcol2,vagprio,vagunits,vagtype

  INTEGER :: vtrstrm, vtrstmovr, vtrstmcol1, vtrstmcol2, vtrstmprio
  INTEGER :: vtpstrm, vtpstmovr, vtpstmcol1, vtpstmcol2, vtpstmprio
  NAMELIST /strmplt_cntl/                                               &
      vtrstrm, vtrstmovr, vtrstmcol1, vtrstmcol2, vtrstmprio,           &
      vtpstrm, vtpstmovr, vtpstmcol1, vtpstmcol2, vtpstmprio
  COMMON /init14_cntl_strm/                                             &
      vtrstrm, vtrstmovr, vtrstmcol1, vtrstmcol2, vtrstmprio,           &
      vtpstrm, vtpstmovr, vtpstmcol1, vtpstmcol2, vtpstmprio


  INTEGER :: trnplt,wetcanplt,raincplt,raingplt,raintplt
  REAL :: trninc,wcpinc,raincinc,rainginc,raintinc
  REAL :: trnminc,wcpminc,raincminc,raingminc,raintminc
  REAL :: trnmaxc,wcpmaxc,raincmaxc,raingmaxc,raintmaxc
  INTEGER :: trnovr,wcpovr,racovr,ragovr,ratovr
  INTEGER :: trncol1,wcpcol1,raccol1,ragcol1,ratcol1
  INTEGER :: trncol2,wcpcol2,raccol2,ragcol2,ratcol2
  INTEGER :: trnprio,wcpprio,racprio,ragprio,ratprio
  INTEGER :: trnhlf,wcphlf,rachlf,raghlf,rathlf
  INTEGER :: trnzro,wcpzro,raczro,ragzro,ratzro
  INTEGER :: trnsty,wcpsty,racsty,ragsty,ratsty
  INTEGER :: racunit, ragunit, ratunit
  INTEGER :: rainicplt,rainigplt,rainitplt
  REAL    :: rainicinc,rainiginc,rainitinc
  REAL    :: rainicminc, rainigminc, rainitminc
  REAL    :: rainicmaxc, rainigmaxc, rainitmaxc
  INTEGER :: raicovr,raigovr,raitovr
  INTEGER :: raiccol1,raigcol1,raitcol1
  INTEGER :: raiccol2,raigcol2,raitcol2
  INTEGER :: raichlf,raighlf,raithlf
  INTEGER :: raicprio,raigprio,raitprio
  INTEGER :: raiczro,raigzro,raitzro
  INTEGER :: raicsty,raigsty,raitsty
  INTEGER :: raicunit,raigunit,raitunit

  NAMELIST /sfc_plot1/                                                  &
      trnplt,trninc,trnminc, trnmaxc,trnovr,trncol1,trncol2,trnprio,    &
            trnhlf, trnzro, trnsty,                                     &
      wetcanplt,wcpinc,wcpminc,wcpmaxc,wcpovr,wcpcol1,wcpcol2,wcpprio,  &
            wcphlf, wcpzro, wcpsty,                                     &
      raincplt,raincinc,raincminc,raincmaxc,racovr,raccol1,raccol2,     &
            rachlf, racprio, raczro, racsty, racunit,                   &
      raingplt,rainginc,raingminc,raingmaxc,ragovr,ragcol1,ragcol2,     &
            raghlf,  ragprio, ragzro, ragsty, ragunit,                  &
      raintplt,raintinc,raintminc,raintmaxc,ratovr,ratcol1,ratcol2,     &
            rathlf,  ratprio, ratzro, ratsty, ratunit,                  &
      rainicplt,rainicinc,rainicminc,rainicmaxc,raicovr,raiccol1,       &
            raiccol2,raichlf,raicprio,raiczro,raicsty,raicunit,         &
      rainigplt,rainiginc,rainigminc,rainigmaxc,raigovr,raigcol1,       &
            raigcol2,raighlf,raigprio,raigzro,raigsty,raigunit,         &
      rainitplt,rainitinc,rainitminc,rainitmaxc,raitovr,raitcol1,       &
            raitcol2,raithlf,raitprio,raitzro,raitsty,raitunit

  COMMON /init15_sfc/               &
      trnovr,trncol1,trncol2,trnprio,trnhlf, trnzro, trnsty,            &
      wetcanplt,wcpinc,wcpminc,wcpmaxc,wcpovr,wcpcol1,wcpcol2,wcpprio,  &
            wcphlf, wcpzro, wcpsty,                                     &
      raincplt,raincinc,raincminc,raincmaxc,racovr,raccol1,raccol2,     &
            rachlf, racprio, raczro, racsty, racunit,                   &
      raingplt,rainginc,raingminc,raingmaxc,ragovr,ragcol1,ragcol2,     &
            raghlf,  ragprio, ragzro, ragsty, ragunit,                  &
      raintplt,raintinc,raintminc,raintmaxc,ratovr,ratcol1,ratcol2,     &
            rathlf,  ratprio, ratzro, ratsty, ratunit,                  &
      rainicplt,rainicinc,rainicminc,rainicmaxc,raicovr,raiccol1,       &
            raiccol2,raichlf,raicprio,raiczro,raicsty,raicunit,         &
      rainigplt,rainiginc,rainigminc,rainigmaxc,raigovr,raigcol1,       &
            raigcol2,raighlf,raigprio,raigzro,raigsty,raigunit,         &
      rainitplt,rainitinc,rainitminc,rainitmaxc,raitovr,raitcol1,       &
            raitcol2,raithlf,raitprio,raitzro,raitsty,raitunit


  INTEGER :: tsoilplt, qsoilplt
  REAL :: tsoilinc, qsoilinc
  REAL :: tsoilminc, qsoilminc
  REAL :: tsoilmaxc, qsoilmaxc
  INTEGER :: tsoilovr, qsoilovr
  INTEGER :: tsoilcol1, qsoilcol1
  INTEGER :: tsoilcol2, qsoilcol2
  INTEGER :: tsoilhlf, qsoilhlf
  INTEGER :: tsoilprio, qsoilprio
  INTEGER :: tsoilzro, qsoilzro

  NAMELIST /soil_plot/                                                  &
      tsoilplt,tsoilinc,tsoilminc,tsoilmaxc,tsoilovr,                   &
            tsoilcol1,tsoilcol2,tsoilhlf,tsoilprio,tsoilzro,            &
      qsoilplt,qsoilinc,qsoilminc,qsoilmaxc,qsoilovr,                   &
            qsoilcol1,qsoilcol2,qsoilhlf,qsoilprio,qsoilzro
  COMMON /init19_soil/               &
      tsoilplt,tsoilinc,tsoilminc,tsoilmaxc,tsoilovr,                   &
            tsoilcol1,tsoilcol2,tsoilhlf,tsoilprio,tsoilzro,            &
      qsoilplt,qsoilinc,qsoilminc,qsoilmaxc,qsoilovr,                   &
            qsoilcol1,qsoilcol2,qsoilhlf,qsoilprio,qsoilzro

  INTEGER :: pslplt,capeplt,cinplt,thetplt,heliplt,uhplt,               &
             brnplt,brnuplt,srlfplt,srmfplt
  REAL :: uhmnhgt,uhmxhgt
  REAL :: pslinc,capeinc,cininc,thetinc,heliinc,uhinc,                  &
             brninc,brnuinc,srlfinc,srmfinc
  REAL :: pslminc,capeminc,cinminc,thetminc,heliminc,uhminc,            &
             brnminc,bruminc,srlminc,srmminc
  REAL :: pslmaxc,capemaxc,cinmaxc,thetmaxc,helimaxc,uhmaxc,            &
             brnmaxc,brumaxc,srlmaxc,srmmaxc
  INTEGER :: pslovr,capovr,cinovr,theovr,helovr,uhovr,                  &
             brnovr,brnuovr,srlfovr,srmfovr
  INTEGER :: pslcol1,capcol1,cincol1,thecol1,helcol1,uhcol1,            &
             brncol1,brnucol1,srlfcol1,srmfcol1
  INTEGER :: pslcol2,capcol2,cincol2,thecol2,helcol2,uhcol2,            &
             brncol2,brnucol2,srlfcol2,srmfcol2
  INTEGER :: pslprio,capprio,cinprio,theprio,helprio,uhprio,            &
             brnprio,bruprio,srlprio,srmprio
  INTEGER :: pslhlf,caphlf,cinhlf,thehlf,helhlf,uhhlf,                  &
             brnhlf,brnuhlf,srlfhlf,srmfhlf
  INTEGER :: pslzro,capzro,cinzro,thezro,helzro,uhzro,                  &
             brnzro,brnuzro,srlfzro,srmfzro
  INTEGER :: pslsty,capsty,cinsty,thesty,helsty,uhsty,                  &
             brnsty,brnusty,srlfsty,srmfsty

  NAMELIST /sfc_plot2/                                                  &
      pslplt,pslinc, pslminc, pslmaxc,pslovr,pslcol1,pslcol2,pslprio,   &
            pslhlf, pslzro, pslsty,                                     &
      capeplt,capeinc,capeminc,capemaxc,capovr,capcol1,capcol2,capprio, &
            caphlf, capzro, capsty,                                     &
      cinplt, cininc, cinminc, cinmaxc, cinovr,cincol1,cincol2,cinprio, &
            cinhlf, cinzro, cinsty,                                     &
      thetplt,thetinc,thetminc,thetmaxc,theovr,thecol1,thecol2,theprio, &
            thehlf, thezro, thesty,                                     &
      heliplt,heliinc,heliminc,helimaxc,helovr,helcol1,helcol2,helprio, &
            helhlf, helzro, helsty,                                     &
      uhplt,uhinc,uhminc,uhmaxc,uhovr,uhcol1,uhcol2,uhprio,             &
            uhhlf, uhzro, uhsty, uhmnhgt, uhmxhgt,                      &
      brnplt, brninc, brnminc, brnmaxc, brnovr,brncol1,brncol2,brnprio, &
            brnhlf, brnzro, brnsty,                                     &
      brnuplt, brnuinc, bruminc, brumaxc, brnuovr, brnucol1,brnucol2,   &
            brnuhlf,  brnuzro, brnusty, bruprio,                        &
      srlfplt, srlfinc, srlminc, srlmaxc, srlfovr, srlfcol1,srlfcol2,   &
            srlfhlf,  srlfzro, srlfsty, srlprio,                        &
      srmfplt, srmfinc, srmminc, srmmaxc, srmfovr, srmfcol1,srmfcol2,   &
            srmfhlf, srmfzro, srmfsty, srmprio
  COMMON /init16_sfc/               &
      pslplt,pslinc, pslminc, pslmaxc,pslovr,pslcol1,pslcol2,pslprio,   &
            pslhlf, pslzro, pslsty,                                     &
      capeplt,capeinc,capeminc,capemaxc,capovr,capcol1,capcol2,capprio, &
            caphlf, capzro, capsty,                                     &
      cinplt, cininc, cinminc, cinmaxc, cinovr,cincol1,cincol2,cinprio, &
            cinhlf, cinzro, cinsty,                                     &
      thetplt,thetinc,thetminc,thetmaxc,theovr,thecol1,thecol2,theprio, &
            thehlf, thezro, thesty,                                     &
      heliplt,heliinc,heliminc,helimaxc,helovr,helcol1,helcol2,helprio, &
            helhlf, helzro, helsty,                                     &
      uhplt,uhinc,uhminc,uhmaxc,uhovr,uhcol1,uhcol2,uhprio,             &
            uhhlf, uhzro, uhsty, uhmnhgt, uhmxhgt,                      &
      brnplt, brninc, brnminc, brnmaxc, brnovr,brncol1,brncol2,brnprio, &
            brnhlf, brnzro, brnsty,                                     &
      brnuplt, brnuinc, bruminc, brumaxc, brnuovr, brnucol1,brnucol2,   &
            brnuhlf,  brnuzro, brnusty, bruprio,                        &
      srlfplt, srlfinc, srlminc, srlmaxc, srlfovr, srlfcol1,srlfcol2,   &
            srlfhlf,  srlfzro, srlfsty, srlprio,                        &
      srmfplt, srmfinc, srmminc, srmmaxc, srmfovr, srmfcol1,srmfcol2,   &
            srmfhlf, srmfzro, srmfsty, srmprio

  INTEGER :: liplt,capsplt,blcoplt,viqcplt,viqiplt,viqrplt,viqsplt,     &
             viqhplt, vilplt
  REAL :: liinc,capsinc,blcoinc,viqcinc,viqiinc,viqrinc,viqsinc,        &
             viqhinc, vilinc
  REAL :: liminc,capsminc,blcominc,viqcminc,viqiminc,viqrminc,viqsminc, &
             viqhminc, vilminc
  REAL :: limaxc,capsmaxc,blcomaxc,viqcmaxc,viqimaxc,viqrmaxc,viqsmaxc, &
             viqhmaxc, vilmaxc
  INTEGER :: liovr,capsovr,blcoovr,viqcovr,viqiovr,viqrovr,viqsovr,     &
             viqhovr, vilovr
  INTEGER :: licol1,capscol1,blcocol1,viqccol1,viqicol1,viqrcol1,viqscol1, &
             viqhcol1, vilcol1
  INTEGER :: licol2,capscol2,blcocol2,viqccol2,viqicol2,viqrcol2,viqscol2, &
             viqhcol2, vilcol2
  INTEGER :: liprio,capsprio,blcoprio,viqcprio,viqiprio,viqrprio,viqsprio, &
             viqhprio, vilprio
  INTEGER :: lihlf,capshlf,blcohlf,viqchlf,viqihlf,viqrhlf,viqshlf,     &
             viqhhlf, vilhlf
  INTEGER :: lizro,capszro,blcozro,viqczro,viqizro,viqrzro,viqszro,     &
             viqhzro, vilzro
  INTEGER :: listy,capssty,blcosty,viqcsty,viqisty,viqrsty,viqssty,     &
             viqhsty, vilsty

  NAMELIST /sfc_plot3/                                                  &
      liplt, liinc, liminc, limaxc, liovr, licol1,licol2,liprio,        &
            lihlf, lizro, listy,                                        &
      capsplt, capsinc, capsminc, capsmaxc, capsovr, capscol1,capscol2, &
            capshlf, capszro, capssty, capsprio,                        &
      blcoplt, blcoinc, blcominc, blcomaxc, blcoovr, blcocol1,blcocol2, &
            blcohlf, blcozro, blcosty, blcoprio,                        &
      viqcplt, viqcinc, viqcminc, viqcmaxc, viqcovr, viqccol1,viqccol2, &
            viqchlf, viqczro, viqcsty, viqcprio,                        &
      viqiplt, viqiinc, viqiminc, viqimaxc, viqiovr, viqicol1,viqicol2, &
            viqihlf, viqizro, viqisty, viqiprio,                        &
      viqrplt, viqrinc, viqrminc, viqrmaxc, viqrovr, viqrcol1,viqrcol2, &
            viqrhlf, viqrzro, viqrsty, viqrprio,                        &
      viqsplt, viqsinc, viqsminc, viqsmaxc, viqsovr, viqscol1,viqscol2, &
            viqshlf, viqszro, viqssty,viqsprio,                         &
      viqhplt, viqhinc, viqhminc, viqhmaxc, viqhovr, viqhcol1,viqhcol2, &
            viqhhlf, viqhzro, viqhsty,viqhprio,                         &
      vilplt, vilinc, vilminc, vilmaxc, vilovr, vilcol1,vilcol2,        &
            vilhlf, vilzro, vilsty, vilprio
  COMMON /init17_sfc/               &
      liplt, liinc, liminc, limaxc, liovr, licol1,licol2,liprio,        &
            lihlf, lizro, listy,                                        &
      capsplt, capsinc, capsminc, capsmaxc, capsovr, capscol1,capscol2, &
            capshlf, capszro, capssty, capsprio,                        &
      blcoplt, blcoinc, blcominc, blcomaxc, blcoovr, blcocol1,blcocol2, &
            blcohlf, blcozro, blcosty, blcoprio,                        &
      viqcplt, viqcinc, viqcminc, viqcmaxc, viqcovr, viqccol1,viqccol2, &
            viqchlf, viqczro, viqcsty, viqcprio,                        &
      viqiplt, viqiinc, viqiminc, viqimaxc, viqiovr, viqicol1,viqicol2, &
            viqihlf, viqizro, viqisty, viqiprio,                        &
      viqrplt, viqrinc, viqrminc, viqrmaxc, viqrovr, viqrcol1,viqrcol2, &
            viqrhlf, viqrzro, viqrsty, viqrprio,                        &
      viqsplt, viqsinc, viqsminc, viqsmaxc, viqsovr, viqscol1,viqscol2, &
            viqshlf, viqszro, viqssty,viqsprio,                         &
      viqhplt, viqhinc, viqhminc, viqhmaxc, viqhovr, viqhcol1,viqhcol2, &
            viqhhlf, viqhzro, viqhsty,viqhprio,                         &
      vilplt, vilinc, vilminc, vilmaxc, vilovr, vilcol1,vilcol2,        &
            vilhlf, vilzro, vilsty, vilprio

  INTEGER :: viiplt,vicplt,ctcplt,vitplt,pwplt,tprplt,gprplt,cprplt
  REAL :: viiinc,vicinc,ctcinc,vitinc,pwinc,tprinc,gprinc,cprinc
  REAL :: viiminc,vicminc,ctcminc,vitminc,pwminc,tprminc,gprminc,cprminc
  REAL :: viimaxc,vicmaxc,ctcmaxc,vitmaxc,pwmaxc,tprmaxc,gprmaxc,cprmaxc
  INTEGER :: viiovr,vicovr,ctcovr,vitovr,pwovr,tprovr,gprovr,cprovr
  INTEGER :: viicol1,viccol1,ctccol1,vitcol1,pwcol1,tprcol1,gprcol1,cprcol1
  INTEGER :: viicol2,viccol2,ctccol2,vitcol2,pwcol2,tprcol2,gprcol2,cprcol2
  INTEGER :: viihlf,vichlf,ctchlf,vithlf,pwhlf,tprhlf,gprhlf,cprhlf
  INTEGER :: viizro,viczro,ctczro,vitzro,pwzro,tprzro,gprzro,cprzro
  INTEGER :: viisty,vicsty,ctcsty,vitsty,pwsty,tprsty,gprsty,cprsty
  INTEGER :: viiprio,vicprio,ctcprio,vitprio,pwprio,tprprio,gprprio,cprprio
  INTEGER :: tprunits, gprunits, cprunits

  NAMELIST /sfc_plot4/                                                  &
      viiplt, viiinc, viiminc, viimaxc, viiovr, viicol1,viicol2,        &
            viihlf, viizro, viisty,  viiprio,                           &
      vicplt, vicinc, vicminc, vicmaxc, vicovr, viccol1,viccol2,        &
            vichlf, viczro, vicsty, vicprio,                            &
      ctcplt, ctcinc, ctcminc, ctcmaxc, ctcovr, ctccol1,ctccol2,        &
            ctchlf, ctczro, ctcsty, ctcprio,                            &
      vitplt, vitinc, vitminc, vitmaxc, vitovr, vitcol1,vitcol2,        &
            vithlf, vitzro, vitsty, vitprio,                            &
      pwplt, pwinc, pwminc, pwmaxc, pwovr, pwcol1,pwcol2,               &
            pwhlf, pwzro, pwsty, pwprio,                                &
      tprplt, tprinc, tprminc, tprmaxc, tprovr, tprcol1,tprcol2,        &
            tprhlf, tprzro, tprsty, tprprio, tprunits,                  &
      gprplt, gprinc, gprminc, gprmaxc, gprovr, gprcol1,gprcol2,        &
            gprhlf, gprzro, gprsty, gprprio, gprunits,                  &
      cprplt, cprinc, cprminc, cprmaxc, cprovr, cprcol1,cprcol2,        &
            cprhlf, cprzro, cprsty, cprprio, cprunits
  COMMON /init18_sfc/               &
      viiplt, viiinc, viiminc, viimaxc, viiovr, viicol1,viicol2,        &
            viihlf, viizro, viisty,  viiprio,                           &
      vicplt, vicinc, vicminc, vicmaxc, vicovr, viccol1,viccol2,        &
            vichlf, viczro, vicsty, vicprio,                            &
      ctcplt, ctcinc, ctcminc, ctcmaxc, ctcovr, ctccol1,ctccol2,        &
            ctchlf, ctczro, ctcsty, ctcprio,                            &
      vitplt, vitinc, vitminc, vitmaxc, vitovr, vitcol1,vitcol2,        &
            vithlf, vitzro, vitsty, vitprio,                            &
      pwplt, pwinc, pwminc, pwmaxc, pwovr, pwcol1,pwcol2,               &
            pwhlf, pwzro, pwsty, pwprio,                                &
      tprplt, tprinc, tprminc, tprmaxc, tprovr, tprcol1,tprcol2,        &
            tprhlf, tprzro, tprsty, tprprio, tprunits,                  &
      gprplt, gprinc, gprminc, gprmaxc, gprovr, gprcol1,gprcol2,        &
            gprhlf, gprzro, gprsty, gprprio, gprunits,                  &
      cprplt, cprinc, cprminc, cprmaxc, cprovr, cprcol1,cprcol2,        &
            cprhlf, cprzro, cprsty, cprprio, cprunits

  INTEGER :: soiltpplt,vegtpplt,laiplt,rouplt,vegplt,snowdplt
  REAL :: soiltpinc,vegtpinc,laiinc,rouinc,veginc,snowdinc
  REAL :: soiltpminc,vegtpminc,laiminc,rouminc,vegminc,snowdminc
  REAL :: soiltpmaxc,vegtpmaxc,laimaxc,roumaxc,vegmaxc,snowdmaxc
  INTEGER :: styovr,vtyovr,laiovr,rouovr,vegovr,snowdovr
  INTEGER :: stycol1,vtycol1,laicol1,roucol1,vegcol1,snowdcol1
  INTEGER :: stycol2,vtycol2,laicol2,roucol2,vegcol2,snowdcol2
  INTEGER :: styprio,vtyprio,laiprio,rouprio,vegprio,snowdprio
  INTEGER :: styhlf,vtyhlf,laihlf,rouhlf,veghlf,snowdhlf
  INTEGER :: styzro,vtyzro,laizro,rouzro,vegzro,snowdzro
  INTEGER :: stysty,vtysty,laisty,rousty,vegsty,snowdsty
  INTEGER :: soiltpn       ! number of soil type 1 to 4

  NAMELIST /sfc_cha_plot/                                               &
      soiltpplt,soiltpinc,soiltpminc,soiltpmaxc,styovr,stycol1,stycol2, &
            styhlf, styzro, stysty,styprio,soiltpn,                     &
      vegtpplt,vegtpinc,vegtpminc,vegtpmaxc,vtyovr,vtycol1,vtycol2,     &
            vtyhlf, vtyzro, vtysty,vtyprio,                             &
      laiplt,laiinc,laiminc,laimaxc,laiovr,laicol1,laicol2,laiprio,     &
            laihlf, laizro, laisty,                                     &
      rouplt,rouinc,rouminc,roumaxc,rouovr,roucol1,roucol2,rouprio,     &
            rouhlf, rouzro, rousty,                                     &
      vegplt,veginc,vegminc,vegmaxc,vegovr,vegcol1,vegcol2,vegprio,     &
            veghlf, vegzro, vegsty,                                     &
      snowdplt,snowdinc,snowdminc,snowdmaxc,snowdovr,snowdcol1,         &
            snowdcol2, snowdprio,snowdhlf, snowdzro, snowdsty
  COMMON /init20_sfccha/               &
      soiltpplt,soiltpinc,soiltpminc,soiltpmaxc,styovr,stycol1,stycol2, &
            styhlf, styzro, stysty,styprio,soiltpn,                     &
      vegtpplt,vegtpinc,vegtpminc,vegtpmaxc,vtyovr,vtycol1,vtycol2,     &
            vtyhlf, vtyzro, vtysty,vtyprio,                             &
      laiplt,laiinc,laiminc,laimaxc,laiovr,laicol1,laicol2,laiprio,     &
            laihlf, laizro, laisty,                                     &
      rouplt,rouinc,rouminc,roumaxc,rouovr,roucol1,roucol2,rouprio,     &
            rouhlf, rouzro, rousty,                                     &
      vegplt,veginc,vegminc,vegmaxc,vegovr,vegcol1,vegcol2,vegprio,     &
            veghlf, vegzro, vegsty,                                     &
      snowdplt,snowdinc,snowdminc,snowdmaxc,snowdovr,snowdcol1,         &
            snowdcol2, snowdprio,snowdhlf, snowdzro, snowdsty

  INTEGER :: setcontopt ,setcontnum
  CHARACTER (LEN=12) :: setcontvar(maxuneva)
  REAL :: setconts(maxunevm,maxuneva)
  NAMELIST /setcont_cntl/setcontopt,setcontnum,setcontvar,setconts
  COMMON /setcont_var/setcontvar
  COMMON /setcon_par/setcontopt,setcontnum,setconts

  INTEGER :: arbvaropt   ! plot arbitrary variable
  INTEGER :: finfmt3d(maxarbvar), finfmt2d(maxarbvar)
  CHARACTER (LEN=256) :: dirname3d(maxarbvar),dirname2d(maxarbvar)
  CHARACTER (LEN=256) :: filename3d(maxarbvar),filename2d(maxarbvar)
  CHARACTER (LEN=6)   :: var3d(maxarbvar),var2d(maxarbvar)
  INTEGER :: var3dnum, var3dplot(maxarbvar)
  REAL    :: var3dinc(maxarbvar), var3dminc(maxarbvar),                 &
             var3dmaxc(maxarbvar)
  INTEGER :: var3dovr(maxarbvar),var3dcol1(maxarbvar),                  &
          var3dcol2(maxarbvar),var3dprio(maxarbvar),                    &
          var3dhlf(maxarbvar),var3dzro(maxarbvar),                      &
          var3dsty(maxarbvar)
  INTEGER :: var2dnum, var2dplot(maxarbvar)
  REAL    :: var2dinc(maxarbvar), var2dminc(maxarbvar),                 &
             var2dmaxc(maxarbvar)
  INTEGER :: var2dovr(maxarbvar),var2dcol1(maxarbvar),                  &
          var2dcol2(maxarbvar), var2dprio(maxarbvar),                   &
          var2dhlf(maxarbvar),var2dzro(maxarbvar),                      &
          var2dsty(maxarbvar)

  INTEGER :: vtr2dnum
  CHARACTER (LEN=256) :: diruv2d(maxarbvar),filenameu2d(maxarbvar),filenamev2d(maxarbvar)
  INTEGER :: finfmtuv2d(maxarbvar)
  CHARACTER (LEN=6)  :: vtru2d(maxarbvar),vtrv2d(maxarbvar)
  INTEGER :: iastride(maxarbvar),jastride(maxarbvar)
  INTEGER :: vtraplt(maxarbvar), magaplt(maxarbvar)
  REAL :: vtraunit(maxarbvar)
  REAL :: magainc(maxarbvar),magaminc(maxarbvar),magamaxc(maxarbvar)
  INTEGER :: vtraovr(maxarbvar), magaovr(maxarbvar)
  INTEGER :: magahlf(maxarbvar), magazro(maxarbvar)
  INTEGER :: vtracol1(maxarbvar), magacol1(maxarbvar)
  INTEGER :: vtracol2(maxarbvar), magacol2(maxarbvar)
  INTEGER :: vtraprio(maxarbvar), magaprio(maxarbvar)
  INTEGER :: vtraunits(maxarbvar), magaunits(maxarbvar)
  INTEGER :: vtratype(maxarbvar), magasty(maxarbvar)

  NAMELIST /arbvar_cntl/arbvaropt,                                      &
      var3dnum,dirname3d,finfmt3d,filename3d,                           &
      var3d,var3dplot, var3dinc, var3dminc,var3dmaxc,                   &
      var3dovr, var3dhlf, var3dzro,var3dsty,var3dcol1, var3dcol2,       &
      var3dprio, var2dnum,dirname2d,finfmt2d, filename2d,               &
      var2d,var2dplot, var2dinc, var2dminc,var2dmaxc,                   &
      var2dovr, var2dhlf, var2dzro, var2dsty, var2dcol1, var2dcol2,     &
      var2dprio,                                                        &
      vtr2dnum,diruv2d,vtru2d,vtrv2d,finfmtuv2d,filenameu2d,filenamev2d,&
      iastride,jastride,                                                &
      vtraplt,vtraunit,vtraovr,vtracol1,vtracol2,vtraprio,              &
      vtraunits,vtratype,                                               &
      magaplt,magainc,magaminc,magamaxc,magaovr,magahlf,magazro,        &
      magasty,magacol1,magacol2,magaunits

  COMMON /init21_cntl_arbvar/arbvaropt,                                 &
      var3dnum,dirname3d,finfmt3d, filename3d,                          &
      var3d,var3dplot, var3dinc, var3dminc,var3dmaxc,                   &
      var3dovr, var3dhlf, var3dzro,var3dsty,var3dcol1, var3dcol2,       &
      var3dprio, var2dnum,dirname2d,finfmt2d, filename2d,               &
      var2d,var2dplot, var2dinc, var2dminc,var2dmaxc,                   &
      var2dovr, var2dhlf, var2dzro, var2dsty, var2dcol1, var2dcol2,     &
      var2dprio

  COMMON /init26_cntl_vctra/ vtr2dnum,diruv2d,vtru2d,vtrv2d,            &
      finfmtuv2d,filenameu2d,filenamev2d,iastride,jastride,             &
      vtraplt, vtraunit,vtraovr,vtracol1,vtracol2,                      &
      vtraprio,vtraunits,vtratype,                                      &
      magaplt,magainc,magaminc,magamaxc,magaovr,magahlf,magazro,        &
      magasty,magacol1,magacol2,magaunits

  INTEGER :: number_of_boxes, boxcol
  REAL :: bctrx(10),bctry(10),blengx(10),blengy(10)
  REAL :: bx1(10), bx2(10),by1(10),by2(10)
  NAMELIST /plot_boxes/ number_of_boxes,boxcol,                         &
      bctrx,bctry,blengx,blengy
  COMMON /boxesopt/number_of_boxes,boxcol,bx1,bx2,by1,by2

  INTEGER :: number_of_polys, polycol
  REAL :: vertx(max_verts,max_polys), verty(max_verts,max_polys)
  NAMELIST /plot_polylines/ number_of_polys,polycol,vertx,verty
  COMMON /polysopt/number_of_polys,polycol,vertx,verty

  INCLUDE 'arpstrajc.inc'
  NAMELIST /plot_trajectories/ trajc_plt_opt,trajc_plt_bgn_time, &
           trajc_plt_end_time,traj_col,trajc_fn_in,ntimes, &
           trajc_lbl_opt,trajc_lbl_frq,trajc_lbl_siz,trajc_mkr_typ, &
           trajc_mkr_frq,trajc_mkr_siz,trajc_lbl_fmt, &
           ntrajc_start, ntrajc_end, ntrajc_stride

  INTEGER :: ovrlaymulopt, ovrmul_num
  CHARACTER (LEN=12) :: ovrname, ovrmulname(50)
  NAMELIST /ovrlay_mul/ovrlaymulopt,ovrname,ovrmul_num,ovrmulname
  COMMON /init22_ovrlay/ovrlaymulopt,ovrname,ovrmul_num,ovrmulname

  INTEGER :: ovrtrn
  NAMELIST /ovr_terrain/ ovrtrn
  REAL :: ztmin,ztmax
  COMMON /trnpar/trnplt,ovrtrn,trninc,trnminc,trnmaxc,                  &
         ztmin,ztmax

  INTEGER :: w3dplt, q3dplt
  REAL :: wisosf,qisosf
  NAMELIST /wirfrm_plot/ w3dplt, wisosf, q3dplt, qisosf
  COMMON /init23_wirfrm/ w3dplt, wisosf, q3dplt, qisosf

  INTEGER :: ovrobs,nsfcobfl,obsset,obscol,obs_marktyp
  CHARACTER (LEN=256) :: sfcobfl(mxsfcobfl)
  REAL :: obs_marksz,obs_valsz
  NAMELIST /plot_obs/ ovrobs,nsfcobfl,sfcobfl,obscol,obs_marktyp,      &
                      obs_marksz,obs_valsz
  COMMON /init24_obs/ nsfcobfl,sfcobfl
  COMMON /obspar/ ovrobs,obsset,obscol,obs_marktyp,                    &
                  obs_marksz,obs_valsz

  INTEGER :: ovrstaopt
  INTEGER :: ovrstam,staset,ovrstan,ovrstav,wrtstax,stacol,markprio
  INTEGER :: nsta_typ,sta_typ(30),sta_marktyp(30),sta_markcol(30)
  REAL :: sta_marksz(30), wrtstad
  CHARACTER (LEN=256) :: stalofl
  NAMELIST /plot_sta/ ovrstaopt,ovrstam,ovrstan,ovrstav,wrtstax,        &
      wrtstad, stacol, markprio, nsta_typ, sta_typ, sta_marktyp,        &
      sta_markcol,sta_marksz,stalofl
  COMMON /sta_par/ ovrstaopt,ovrstam,staset,ovrstan,ovrstav,stacol,     &
         markprio,nsta_typ,sta_typ,sta_marktyp,sta_markcol,             &
         sta_marksz,stalofl,wrtstax,wrtstad

!---------------------------------------------------------------------
!
!  *min  -- Profile plot lower bound
!  *max  -- Profile plot upper bound
!
!----------------------------------------------------------------------

  INTEGER :: profopt,nprof,npicprof
  REAL    :: xprof(max_dim), yprof(max_dim)
  INTEGER :: uprof, vprof, wprof, ptprof, pprof, qvprof,                &
             qcprof,qrprof,qiprof,qsprof,qhprof,rhprof,                 &
             kmhprof,kmvprof,tkeprof,rfprof,pteprf,                     &
             upprof,vpprof,wpprof,ptpprf,ppprof,qvpprf,                 &
             vorpprf,divpprf, tsoilprof,qsoilprof
  REAL :: uprmin, vprmin, wprmin, ptprmin, pprmin,qvprmin,              &
          qcpmin,qrpmin,qipmin,qspmin,qhpmin,rhpmin,                    &
          kmhpmin,kmvpmin,tkepmin,rfpmin,ptepmin,                       &
          uppmin,vppmin,wppmin,ptppmin,pppmin,qvppmin,                  &
          vorppmin,divppmin,tsoilprofmin,qsoilprofmin
  REAL :: uprmax, vprmax, wprmax, ptprmax, pprmax,qvprmax,              &
          qcpmax,qrpmax,qipmax,qspmax,qhpmax,rhpmax,                    &
          kmhpmax,kmvpmax,tkepmax,rfpmax,ptepmax,                       &
          uppmax,vppmax,wppmax,ptppmax,pppmax,qvppmax,                  &
          vorppmax,divppmax,tsoilprofmax,qsoilprofmax

  REAL :: zprofbgn, zprofend, zsoilprofbgn, zsoilprofend
  INTEGER :: nxprpic, nyprpic

  NAMELIST /profile_cntl/ profopt, nprof, xprof, yprof,                 &
      npicprof, uprof, uprmin, uprmax, vprof, vprmin, vprmax,           &
      wprof,wprmin,wprmax,  ptprof,ptprmin,ptprmax,                     &
      pprof,pprmin,pprmax,  qvprof,qvprmin,qvprmax,                     &
      qcprof,qcpmin,qcpmax, qrprof,qrpmin,qrpmax,                       &
      qiprof,qipmin,qipmax, qsprof,qspmin,qspmax,                       &
      qhprof,qhpmin,qhpmax, rhprof,rhpmin,rhpmax,                       &
      kmhprof,kmhpmin,kmhpmax, kmvprof,kmvpmin,kmvpmax,                 &
      tkeprof,tkepmin,tkepmax,                                          &
      rfprof,rfpmin,rfpmax, pteprf,ptepmin,ptepmax,                     &
      upprof,uppmin,uppmax, vpprof,vppmin,vppmax,                       &
      wpprof,wppmin,wppmax, ptpprf,ptppmin,ptppmax,                     &
      ppprof,pppmin,pppmax, qvpprf,qvppmin,qvppmax,                     &
      vorpprf, vorppmin, vorppmax, divpprf, divppmin, divppmax,         &
      zprofbgn,zprofend,                                                &
      tsoilprof,tsoilprofmin,tsoilprofmax,                              &
      qsoilprof,qsoilprofmin,qsoilprofmax,                              &
      zsoilprofbgn,zsoilprofend,                                        &
      nxprpic, nyprpic
  COMMON /init25_prof/ profopt, nprof, xprof, yprof,                    &
      npicprof, uprof, uprmin, uprmax, vprof, vprmin, vprmax,           &
      wprof,wprmin,wprmax,  ptprof,ptprmin,ptprmax,                     &
      pprof,pprmin,pprmax,  qvprof,qvprmin,qvprmax,                     &
      qcprof,qcpmin,qcpmax, qrprof,qrpmin,qrpmax,                       &
      qiprof,qipmin,qipmax, qsprof,qspmin,qspmax,                       &
      qhprof,qhpmin,qhpmax, rhprof,rhpmin,rhpmax,                       &
      kmhprof,kmhpmin,kmhpmax, kmvprof,kmvpmin,kmvpmax,                 &
      tkeprof,tkepmin,tkepmax,                                          &
      rfprof,rfpmin,rfpmax, pteprf,ptepmin,ptepmax,                     &
      upprof,uppmin,uppmax, vpprof,vppmin,vppmax,                       &
      wpprof,wppmin,wppmax, ptpprf,ptppmin,ptppmax,                     &
      ppprof,pppmin,pppmax, qvpprf,qvppmin,qvppmax,                     &
      vorpprf, vorppmin, vorppmax, divpprf, divppmin, divppmax,         &
      zprofbgn,zprofend,                                                &
      tsoilprof,tsoilprofmin,tsoilprofmax,                              &
      qsoilprof,qsoilprofmin,qsoilprofmax,                              &
      zsoilprofbgn,zsoilprofend,                                        &
      nxprpic, nyprpic

  CHARACTER(LEN=256), INTENT(OUT) :: outfilename
  INTEGER,            INTENT(OUT) :: iwtype

  NAMELIST /output/ dirname,outfilename,iwtype,lvldbg

!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
  INTEGER :: lmapfile, lenstr
  INTEGER :: ireturn
  INTEGER :: lengbf,nf,lenfil
  INTEGER :: indxslic
  LOGICAL :: fexist
  INTEGER :: i,j,k,nq,ifile
  INTEGER :: nxlg, nylg

  CHARACTER(LEN=256) :: strtmp

  INTEGER :: unum         ! unit number for reading in namelist

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  istatus = 0

  CALL mpinit_proc(0)

  IF(myproc == 0) THEN
    WRITE(6,'(/ 16(/5x,a)//)')                                          &
     '###############################################################', &
     '###############################################################', &
     '####                                                       ####', &
     '####                 Welcome to ARPSPLT                    ####', &
     '####                                                       ####', &
     '####       A graphic analysis program for model ARPS 5.3   ####', &
     '####                                                       ####', &
     '####            The graphic plotting is based              ####', &
     '####              on graphic package ZXPLOT                ####', &
     '####               by Ming Xue CAPS/SOM/OU                 ####', &
     '####            (http://www.caps.ou.edu/ZXPLOT)            ####', &
     '####                                                       ####', &
     '###############################################################', &
     '###############################################################'


    unum = COMMAND_ARGUMENT_COUNT()
    IF (unum > 0) THEN
      CALL GET_COMMAND_ARGUMENT(1, strtmp, lenstr, istatus )
      IF ( strtmp(1:1) == ' ' .OR. istatus /= 0 ) THEN  ! Use standard input to be backward-compatible
        unum = 5
      ELSE
        INQUIRE(FILE=TRIM(strtmp),EXIST=fexist)
        IF (.NOT. fexist) THEN
          WRITE(6,'(1x,3a)') 'WARNING: namelist file - ',               &
                TRIM(strtmp),' does not exist. Falling back to standard input.'
          unum = 5
        END IF
      END IF
    ELSE
      unum = 5
    END IF

    IF (unum /= 5) THEN
      CALL getunit( unum )
      OPEN(unum,FILE=TRIM(strtmp),STATUS='OLD',FORM='FORMATTED')
      WRITE(*,'(1x,3a,/,1x,a,/)') 'Reading ARPSPLT namelist from file - ', &
              TRIM(strtmp),' ... ','========================================'
    ELSE
      WRITE(*,'(2(1x,a,/))') 'Waiting namelist from standard input ... ', &
                             '========================================'
    END IF

  END IF
!-----------------------------------------------------------------------
!
! First, initialize MPI jobs
!
!-----------------------------------------------------------------------

  nprocx_lw = 1
  nprocy_lw = 1
  nproc_node = 1

  IF(myproc == 0) THEN
    READ(unum,message_passing,ERR=100)
    WRITE(6,'(a)')'Namelist message_passing was successfully read.'
  END IF
  CALL mpupdatei(nproc_x,1)
  CALL mpupdatei(nproc_y,1)
  CALL mpupdatei(max_fopen,1)
  CALL mpupdatei(nproc_node,1)
!  CALL mpupdatei(readsplit,1)

  CALL mpupdatei(nprocx_in,1)
  CALL mpupdatei(nprocy_in,1)

  CALL mpinit_var

  IF (mp_opt == 0) THEN     ! no-mpi specific
    nproc_node = 1
    ncompressx = nprocx_in
    ncompressy = nprocy_in
  ELSE                      ! mpi specific
    readstride = max_fopen
    dumpstride = nprocs

    IF (nprocx_in == 1 .AND. nprocy_in == 1) THEN
      readsplit(:) = 1
      readstride = nprocs
      nproc_node = 1

      ncompressx = 1
      ncompressy = 1
    ELSE
      readsplit(:) = 0

      ncompressx = nprocx_in/nproc_x
      ncompressy = nprocy_in/nproc_y

      IF ( (MOD(nprocx_in,nproc_x) /= 0 .OR. MOD(nprocy_in,nproc_y) /= 0)   &
           .OR. (ncompressx < 1 .OR. ncompressy < 1)  ) THEN
        IF (myproc == 0) WRITE(6,'(3x,a/,2(3x,2(a,I2)/))')                  &
          'nprocx_in (nprocy_in) must be multiples of nproc_x(nproc_y)',    &
          'nprocx_in = ',nprocx_in, ',nprocy_in = ',nprocy_in,              &
          '  nproc_x = ',nproc_x,  ',  nproc_y = ', nproc_y
        CALL mpexit(1);
        STOP
      END IF

    END IF

    IF (nproc_node <= 1) THEN      ! ignore nproc_node
      nproc_node = 1
    ELSE                           ! ignore max_fopen
      readstride = nprocs
    END IF
  END IF

!
!-----------------------------------------------------------------------
!
!  Get the names of the input data files.
!
!-----------------------------------------------------------------------
!
  IF(myproc == 0) THEN
    CALL get_input_file_names(unum,hinfmt,grdbasfn,hisfile,nhisfile)
    lengbf = len_trim(grdbasfn)

    IF(mp_opt > 0 .AND. readsplit(FINDX_H) <= 0) THEN
      CALL gtsplitfn(hisfile(1),1,1,loc_x,loc_y,1,1,                    &
                     nprocx_lw-1,nprocy_lw-1,1,lvldbg,strtmp,ireturn)
      lenfil = LEN_TRIM(strtmp)
    ELSE IF (mp_opt == 0 .AND. (nprocx_in > 1 .OR. nprocy_in > 1) ) THEN
      CALL gtsplitfn(hisfile(1),1,1,loc_x,loc_y,1,1,                    &
                     nprocx_lw-1,nprocy_lw-1,1,lvldbg,strtmp,ireturn)
      lenfil = LEN_TRIM(strtmp)
    ELSE
      strtmp = hisfile(1)
      lenfil = LEN_TRIM(hisfile(1))
    END IF

    CALL get_dims_from_data(hinfmt,strtmp(1:lenfil),                    &
                            nx,ny,nz,nzsoil,nstyps, ireturn)

    IF( ireturn /= 0 ) THEN
      PRINT*,'Problem occured when trying to get dimensions from data.'
      PRINT*,'Program stopped.'
      !STOP
    END IF

    IF (mp_opt > 0) THEN
      IF( readsplit(FINDX_H) > 0 ) THEN
        IF( MOD(nx-fzone,nproc_x) /= 0 .OR. MOD(ny-fzone,nproc_y) /= 0) THEN
          WRITE(6,'(a/,a/,4(a,i5))')      &
             'The specification of nproc_x or nproc_y is not matched with nx or ny.',&
             'nx-3 and ny-3 must be multiples of nproc_x and nproc_y respectively.', &
             'nx = ', nx, ' ny = ', ny, ' nproc_x = ',nproc_x, ' nproc_y = ',nproc_y
          nx = 0
          ny = 0
        ELSE
          nx = (nx-fzone)/nproc_x + fzone
          ny = (ny-fzone)/nproc_y + fzone
        END IF
      ELSE
        nx = (nx-fzone)*ncompressx + fzone
        ny = (ny-fzone)*ncompressy + fzone
      END IF
    ELSE IF (nprocx_in > 1 .OR. nprocy_in > 1) THEN
      nx = (nx-fzone)*ncompressx + fzone
      ny = (ny-fzone)*ncompressy + fzone
    END IF

  END IF ! myproc == 0
  CALL mpupdatei(ireturn,1)
  IF( ireturn /= 0 ) THEN
    !PRINT*,'Problem occured when trying to get dimensions from data.'
    !PRINT*,'Program stopped.'
    CALL arpsstop('Error when trying to get dimensions from data',1)
  END IF

  CALL mpupdatei(hinfmt,1)
  CALL mpupdatec(grdbasfn,256)
  CALL mpupdatei(nhisfile,1)
  CALL mpupdatec(hisfile,256*nhisfile_max)

  CALL mpupdatei(nx,1)
  CALL mpupdatei(ny,1)
  CALL mpupdatei(nz,1)
  CALL mpupdatei(nzsoil,1)
  CALL mpupdatei(nstyps,1)
  CALL mpupdatei(nscalar,1)
  CALL mpupdatei(nscalarq,1)
  CALL mpupdatei(P_QC,1)
  CALL mpupdatei(P_QR,1)
  CALL mpupdatei(P_QI,1)
  CALL mpupdatei(P_QS,1)
  CALL mpupdatei(P_QH,1)
  CALL mpupdatei(P_QG,1)
  CALL mpupdatei(P_NC,1)
  CALL mpupdatei(P_NR,1)
  CALL mpupdatei(P_NI,1)
  CALL mpupdatei(P_NS,1)
  CALL mpupdatei(P_NH,1)
  CALL mpupdatei(P_NG,1)
  CALL mpupdatei(P_ZR,1)
  CALL mpupdatei(P_ZI,1)
  CALL mpupdatei(P_ZS,1)
  CALL mpupdatei(P_ZH,1)
  CALL mpupdatei(P_ZG,1)
  CALL mpupdatec(qnames,nscalar*40)
  CALL mpupdatec(qdescp,nscalar*40)

  IF( nx <= 0 .OR. ny <= 0 ) CALL mpexit(1);

  nstyp = nstyps            ! Copy to global variable

  IF(myproc == 0) THEN
    WRITE(6,'(2x,5(a,i5))') 'nx = ',nx,', ny = ',ny,', nz = ',nz,       &
                            ', nzsoil = ',nzsoil,', nstyps = ', nstyps
  END IF

  nxlg = (nx-3)*nproc_x + 3
  nylg = (ny-3)*nproc_y + 3

!-----------------------------------------------------------------------
! Set certain defaul options / values
!-----------------------------------------------------------------------

  msfplt = 0
  ipvplt = 0
  vagplt = 0
  thkplt = 0

  paprlnth = 1.5  ! default value
  lnmag = 1       ! default value

  lblmag = 1.0
  winsiz = 1.0
  margnx = 0.1
  margny = 0.1
  pcolbar = 1
  axlbfmt = -1
  axlbsiz = 0.025
  tickopt=0

  ctrlbopt  = 1
  ctrstyle  = 1
  ctrlbfrq  = 2
  ctrlbsiz  = 0.02
  lbmaskopt = 0

  istride = 0
  jstride = 0
  kstride = 0

  iskip = 0
  jskip = 0

  uhmnhgt = 2000.
  uhmxhgt = 5000.

!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!
!  Read in plotting control parameters
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!
!
!-----------------------------------------------------------------------
!
!  Page set-up parameters
!
!-----------------------------------------------------------------------
!
  IF(myproc == 0) THEN
    READ(unum,page_setup,ERR=100)
    WRITE(6,'(a)')'Namelist page_setup was successfully read.'
  END IF
  CALL mpupdatei(layout,1)
  CALL mpupdatei(nxpic,1)
  CALL mpupdatei(nypic,1)
  CALL mpupdatei(inwfrm,1)
  CALL mpupdater(paprlnth,1)

  IF(myproc == 0) THEN
    READ(unum,plotting_setup,ERR=100)
    WRITE(6,'(a)')'Namelist plotting_setup was successfully read.'
  END IF
  CALL mpupdatei(iorig,1)
  CALL mpupdater(xorig,1)
  CALL mpupdater(yorig,1)
  CALL mpupdater(xbgn,1)
  CALL mpupdater(xend,1)
  CALL mpupdater(ybgn,1)
  CALL mpupdater(yend,1)
  CALL mpupdater(zbgn,1)
  CALL mpupdater(zend,1)
  CALL mpupdater(zsoilbgn,1)
  CALL mpupdater(zsoilend,1)
  CALL mpupdater(yxstrch,1)
  CALL mpupdater(zxstrch,1)
  CALL mpupdater(zystrch,1)
  CALL mpupdater(zhstrch,1)
  CALL mpupdater(zsoilxstrch,1)
  CALL mpupdater(zsoilystrch,1)
  CALL mpupdater(winsiz,1)
  CALL mpupdater(margnx,1)
  CALL mpupdater(margny,1)
  CALL mpupdater(pcolbar,1)
  CALL mpupdatei(iskip,1)
  CALL mpupdatei(jskip,1)

  IF(myproc == 0) THEN
    READ(unum,col_table_cntl, ERR=100)
    WRITE(6,'(a)')'namelist col_table_cntl was successfully read.'
    WRITE(6,'(2x,a,i3)') 'Color table is : ',col_table
  END IF
  CALL mpupdatei(col_table,1)
  CALL mpupdatec(color_map,256)

  IF(myproc == 0) THEN
    READ(unum,style_tuning,ERR=100)
    WRITE(6,'(a)')'Namelist style_tuning was successfully read.'
  END IF
  CALL mpupdater(lblmag,1)
  CALL mpupdatei(lnmag,1)
  CALL mpupdatei(fontopt,1)
  CALL mpupdatei(lbaxis,1)
  CALL mpupdatei(axlbfmt,1)
  CALL mpupdater(axlbsiz,1)
  CALL mpupdatei(haxisu,1)
  CALL mpupdatei(vaxisu,1)
  CALL mpupdatei(tickopt,1)
  CALL mpupdater(hmintick,1)
  CALL mpupdater(vmajtick,1)
  CALL mpupdater(vmintick,1)
  CALL mpupdater(hmajtick,1)
  CALL mpupdatei(presaxis_no,1)
  CALL mpupdater(pres_val,20)
  CALL mpupdatei(ctrlbopt,1)
  CALL mpupdatei(ctrstyle,1)
  CALL mpupdatei(ctrlbfrq,1)
  CALL mpupdater(ctrlbsiz,1)
  CALL mpupdatei(lbmaskopt,1)

  IF(myproc == 0) THEN
    READ(unum,smooth_cntl, ERR=100)
    WRITE(6,'(a)')'Namelist smooth_cntl was successfully read.'
    WRITE(6,'(2x,a,i3)') 'Smoothing option is : ',smooth
  END IF
  CALL mpupdatei(smooth,1)


  IF(myproc == 0) THEN
    READ(unum,title_setup, ERR=100)
    WRITE(6,'(a)')'Namelist title_setup was successfully read.'
  END IF
  CALL mpupdatei(ntitle,1)
  CALL mpupdatei(titcol,1)
  CALL mpupdater(titsiz,1)
  CALL mpupdatec(title,3*256)

  IF(myproc == 0) THEN
    READ(unum,footer_setup, ERR=100)
    WRITE(6,'(a)')'Namelist footer_setup was successfully read.'
  END IF
  CALL mpupdatei(wpltime,1)
  CALL mpupdatec(footer_l,256)
  CALL mpupdatec(footer_c,256)
  CALL mpupdatec(footer_r,256)
!
!-----------------------------------------------------------------------
!
!  Input control parameters map plotting
!
!-----------------------------------------------------------------------
!
  mapgrid = 0 ! no longer used.

  IF(myproc == 0) THEN
    READ(unum,map_plot,ERR=100)
    WRITE(6,'(a)')'Namelist map_plot was successfully read.'

    IF(nmapfile > maxmap)                                                 &
        WRITE(6,'(a)')'Warning: the maximum map files should be ',maxmap

    DO i=1,nmapfile
      lmapfile = LEN_TRIM(mapfile(i))
      WRITE(6,'(2x,a,a)') 'Input was ',mapfile(i)(1:lmapfile)

      IF(ovrmap == 1) THEN
        INQUIRE(FILE=mapfile(i)(1:lmapfile), EXIST = fexist )
        IF( .NOT.fexist) THEN
          WRITE(6,'(a)') 'Warning: Map file '//mapfile(i)(1:lmapfile)     &
              //' not found. Program will be continue'
        END IF
      END IF
    END DO
  END IF
  CALL mpupdatei(ovrmap,1)
  CALL mpupdatei(mapgrid,1)
  CALL mpupdater(latgrid,1)
  CALL mpupdater(longrid,1)
  CALL mpupdatei(mapgridcol,1)
  CALL mpupdatei(nmapfile,1)
  CALL mpupdatec(mapfile,256*maxmap)
  CALL mpupdatei(mapcol,maxmap)
  CALL mpupdatei(mapline_style,maxmap)

  IF(myproc == 0) THEN
    READ(unum,multi_setup,ERR=100)
    WRITE(6,'(a)')'Namelist multi_setup was successfully read.'
  END IF
  CALL mpupdatei(missfill_opt,1)
  CALL mpupdatei(missval_colind,1)

!
!-----------------------------------------------------------------------
!
!  Input control parameters plotting type
!
!-----------------------------------------------------------------------
!
  IF(myproc == 0) THEN
    READ(unum,xy_slice_cntl,ERR=100)
    WRITE(6,'(a)')'Namelist xy_slice_cntl was successfully read.'
  END IF
  CALL mpupdatei(nslice_xy,1)
  CALL mpupdatei(slice_xy,max_dim)

  IF(myproc == 0) THEN
    READ(unum,xz_slice_cntl,ERR=100)
    WRITE(6,'(a)')'Namelist xz_slice_cntl was successfully read.'
  END IF
  CALL mpupdatei(nslice_xz,1)
  CALL mpupdatei(slice_xz,max_dim)

  IF(myproc == 0) THEN
    READ(unum,yz_slice_cntl,ERR=100)
    WRITE(6,'(a)')'Namelist yz_slice_cntl was successfully read.'
  END IF
  CALL mpupdatei(nslice_yz,1)
  CALL mpupdatei(slice_yz,max_dim)

  IF(myproc == 0) THEN
    READ(unum,h_slice_cntl,ERR=100)
    WRITE(6,'(a)')'Namelist h_slice_cntl was successfully read.'

    IF(abs(nslice_h) > max_dim ) THEN
      WRITE(6,'(a,i3,a)') 'Please give value smaller than ',          &
          nz,'. Program stopped !'
      CALL mpexit(0)
    END IF
  END IF
  CALL mpupdatei(nslice_h,1)
  CALL mpupdater(slice_h,max_dim)

  IF(myproc == 0) THEN
    READ(unum,xy_soil_slice_cntl,ERR=100)
    WRITE(6,'(a)')'Namelist xy_soil_slice_cntl was successfully read.'
  END IF
  CALL mpupdatei(nslice_xy_soil,1)
  CALL mpupdatei(slice_xy_soil,max_dim)

  IF(myproc == 0) THEN
    READ(unum,xz_soil_slice_cntl,ERR=100)
    WRITE(6,'(a)')'Namelist xz_soil_slice_cntl was successfully read.'
  END IF
  CALL mpupdatei(nslice_xz_soil,1)
  CALL mpupdatei(slice_xz_soil,max_dim)

  IF(myproc == 0) THEN
    READ(unum,yz_soil_slice_cntl,ERR=100)
    WRITE(6,'(a)')'Namelist yz_soil_slice_cntl was successfully read.'
  END IF
  CALL mpupdatei(nslice_yz_soil,1)
  CALL mpupdatei(slice_yz_soil,max_dim)

  IF(myproc == 0) THEN
    READ(unum,v_slice_cntl,ERR=100)
    WRITE(6,'(a)')'Namelist v_slice_cntl was successfully read.'

    IF(nslice_v > nxlg+nylg) THEN
      WRITE(6,'(1x,a,i3,a)') ' Please give a value smaller than ',        &
          nx+ny,'. Program stopped !'
      CALL mpexit(0)
    END IF
  END IF
  CALL mpupdatei(nslice_v,1)
  CALL mpupdater(xpnt1,max_dim)
  CALL mpupdater(ypnt1,max_dim)
  CALL mpupdater(xpnt2,max_dim)
  CALL mpupdater(ypnt2,max_dim)

  IF(myproc == 0) THEN
    READ(unum,p_slice_cntl,ERR=100)
    WRITE(6,'(a)')'Namelist p_slice_cntl was successfully read.'

    IF (nslice_p > 1) THEN
      WRITE(6,'(2x,a)',ADVANCE='NO') 'Pressure slices to be plotted are at p = '
      DO indxslic = 1,nslice_p
        WRITE(6,'(f10.3,a)',ADVANCE='NO') slice_p(indxslic),','
      END DO
      WRITE(*,*)

      IF(nslice_p > nz) THEN
        WRITE(6,'(1x,a,i3,a)') 'Please give value smaller than ',       &
            nz,'. Program stopped !'
        CALL mpexit(0)
      END IF
    END IF

  END IF
  CALL mpupdatei(nslice_p,1)
  CALL mpupdater(slice_p,max_dim)

  IF(myproc == 0) THEN
    READ(unum,pt_slice_cntl,ERR=100)
    WRITE(6,'(a)')'Namelist pt_slice_cntl was successfully read.'

    IF (nslice_pt > 0) THEN
      IF(nslice_pt > max_dim) THEN
        WRITE(6,'(1x,2(a,i3))')                                         &
            'Warning: Maximum number of PT slices allowed is max_dim= ',&
            max_dim,'nslice_pt reset t= ',max_dim
        nslice_pt = max_dim
      END IF

      WRITE(6,'(2x,a)')'Isentropic slices to be plotted are at theta=:'
      DO indxslic = 1,nslice_pt
        WRITE(6,'(1x,f10.3)') slice_pt(indxslic)
      END DO
    END IF
  END IF
  CALL mpupdatei(nslice_pt,1)
  CALL mpupdater(slice_pt,max_dim)

  IF(myproc == 0) THEN
    READ(unum,domain_move,ERR=100)
    WRITE(6,'(a)')'Namelist domain_move was successfully read.'
  END IF
  CALL mpupdatei(imove,1)
  CALL mpupdater(umove,1)
  CALL mpupdater(vmove,1)

  IF(myproc == 0) THEN
    READ(unum,sclrplt_cntl1, ERR=100)
    WRITE(6,'(a)')'Namelist sclrplt_cntl1 was successfully read.'
  END IF

  CALL mpupdatei(hplot,1)
  CALL mpupdater(hinc,1)
  CALL mpupdater(hminc,1)
  CALL mpupdater(hmaxc,1)
  CALL mpupdatei(hovr,1)
  CALL mpupdatei(hcol1,1)
  CALL mpupdatei(hcol2,1)
  CALL mpupdatei(hprio,1)
  CALL mpupdatei(hhlf,1)
  CALL mpupdatei(hzro,1)
  CALL mpupdatei(hsty,1)

  CALL mpupdatei(msfplt,1)
  CALL mpupdater(msfinc,1)
  CALL mpupdater(msfminc,1)
  CALL mpupdater(msfmaxc,1)
  CALL mpupdatei(msfovr,1)
  CALL mpupdatei(msfcol1,1)
  CALL mpupdatei(msfcol2,1)
  CALL mpupdatei(msfprio,1)
  CALL mpupdatei(msfhlf,1)
  CALL mpupdatei(msfzro,1)
  CALL mpupdatei(msfsty,1)

  CALL mpupdatei(thkplt,1)
  CALL mpupdater(thkinc,1)
  CALL mpupdater(thkminc,1)
  CALL mpupdater(thkmaxc,1)
  CALL mpupdatei(thkovr,1)
  CALL mpupdatei(thkcol1,1)
  CALL mpupdatei(thkcol2,1)
  CALL mpupdatei(thkprio,1)
  CALL mpupdatei(thkhlf,1)
  CALL mpupdatei(thkzro,1)
  CALL mpupdatei(thksty,1)

  CALL mpupdatei(tplot,1)
  CALL mpupdater(tinc,1)
  CALL mpupdater(tminc,1)
  CALL mpupdater(tmaxc,1)
  CALL mpupdatei(tovr,1)
  CALL mpupdatei(tcol1,1)
  CALL mpupdatei(tcol2,1)
  CALL mpupdatei(tprio,1)
  CALL mpupdatec(tunits,1)
  CALL mpupdatei(thlf,1)
  CALL mpupdatei(tzro,1)
  CALL mpupdatei(tsty,1)

  CALL mpupdatei(uplot,1)
  CALL mpupdater(uinc,1)
  CALL mpupdater(uminc,1)
  CALL mpupdater(umaxc,1)
  CALL mpupdatei(uovr,1)
  CALL mpupdatei(ucol1,1)
  CALL mpupdatei(ucol2,1)
  CALL mpupdatei(uprio,1)
  CALL mpupdatei(uhlf,1)
  CALL mpupdatei(uzro,1)
  CALL mpupdatei(usty,1)

  CALL mpupdatei(vplot,1)
  CALL mpupdater(vinc,1)
  CALL mpupdater(vminc,1)
  CALL mpupdater(vmaxc,1)
  CALL mpupdatei(vovr,1)
  CALL mpupdatei(vcol1,1)
  CALL mpupdatei(vcol2,1)
  CALL mpupdatei(vprio,1)
  CALL mpupdatei(vhlf,1)
  CALL mpupdatei(vzro,1)
  CALL mpupdatei(vsty,1)

  CALL mpupdatei(vhplot,1)
  CALL mpupdater(vhinc,1)
  CALL mpupdater(vhminc,1)
  CALL mpupdater(vhmaxc,1)
  CALL mpupdatei(vhovr,1)
  CALL mpupdatei(vhcol1,1)
  CALL mpupdatei(vhcol2,1)
  CALL mpupdatei(vhprio,1)
  CALL mpupdatei(vhunits,1)
  CALL mpupdatei(vhhlf,1)
  CALL mpupdatei(vhzro,1)
  CALL mpupdatei(vhsty,1)

  CALL mpupdatei(vsplot,1)
  CALL mpupdater(vsinc,1)
  CALL mpupdater(vsminc,1)
  CALL mpupdater(vsmaxc,1)
  CALL mpupdatei(vsovr,1)
  CALL mpupdatei(vscol1,1)
  CALL mpupdatei(vscol2,1)
  CALL mpupdatei(vsprio,1)
  CALL mpupdatei(vshlf,1)
  CALL mpupdatei(vszro,1)
  CALL mpupdatei(vssty,1)

  CALL mpupdatei(wplot,1)
  CALL mpupdater(winc,1)
  CALL mpupdater(wminc,1)
  CALL mpupdater(wmaxc,1)
  CALL mpupdatei(wovr,1)
  CALL mpupdatei(wcol1,1)
  CALL mpupdatei(wcol2,1)
  CALL mpupdatei(wprio,1)
  CALL mpupdatei(whlf,1)
  CALL mpupdatei(wzro,1)
  CALL mpupdatei(wsty,1)

  CALL mpupdatei(ptplot,1)
  CALL mpupdater(ptinc,1)
  CALL mpupdater(ptminc,1)
  CALL mpupdater(ptmaxc,1)
  CALL mpupdatei(ptovr,1)
  CALL mpupdatei(ptcol1,1)
  CALL mpupdatei(ptcol2,1)
  CALL mpupdatei(ptprio,1)
  CALL mpupdatei(pthlf,1)
  CALL mpupdatei(ptzro,1)
  CALL mpupdatei(ptsty,1)

  CALL mpupdatei(pplot,1)
  CALL mpupdater(pinc,1)
  CALL mpupdater(pminc,1)
  CALL mpupdater(pmaxc,1)
  CALL mpupdatei(povr,1)
  CALL mpupdatei(pcol1,1)
  CALL mpupdatei(pcol2,1)
  CALL mpupdatei(pprio,1)
  CALL mpupdatei(phlf,1)
  CALL mpupdatei(pzro,1)
  CALL mpupdatei(psty,1)

  CALL mpupdatei(ipvplt,1)
  CALL mpupdater(ipvinc,1)
  CALL mpupdater(ipvminc,1)
  CALL mpupdater(ipvmaxc,1)
  CALL mpupdatei(ipvovr,1)
  CALL mpupdatei(ipvcol1,1)
  CALL mpupdatei(ipvcol2,1)
  CALL mpupdatei(ipvprio,1)
  CALL mpupdatei(ipvhlf,1)
  CALL mpupdatei(ipvzro,1)
  CALL mpupdatei(ipvsty,1)

  ! Make sure all plotting control arrays for qscalar are initialized
  ! in the case that less than nscalar are explicitly listed in the
  ! namelist file

  DO nq=1,nscalarmax
    qscalarplot(nq)=0
    qscalarinc(nq)=0
    qscalarminc(nq)=0
    qscalarmaxc(nq)=0
    qscalarovr(nq)=0
    qscalarcol1(nq)=0
    qscalarcol2(nq)=0
    qscalarprio(nq)=0
    qscalarhlf(nq)=0
    qscalarzro(nq)=0
    qscalarsty(nq)=0

  END DO

  IF(myproc == 0) THEN
    READ(unum,sclrplt_cntl2, ERR=100)
    WRITE(6,'(a)')'Namelist sclrplt_cntl2 was successfully read.'
  END IF
  CALL mpupdatei(qvplot,1)
  CALL mpupdater(qvinc,1)
  CALL mpupdater(qvminc,1)
  CALL mpupdater(qvmaxc,1)
  CALL mpupdatei(qvovr,1)
  CALL mpupdatei(qvcol1,1)
  CALL mpupdatei(qvcol2,1)
  CALL mpupdatei(qvprio,1)

  CALL mpupdatei(qvhlf,1)
  CALL mpupdatei(qvzro,1)
  CALL mpupdatei(qvsty,1)

    CALL mpupdatei(qscalarplot,nscalarmax)
    CALL mpupdater(qscalarinc,nscalarmax)
    CALL mpupdater(qscalarminc,nscalarmax)
    CALL mpupdater(qscalarmaxc,nscalarmax)
    CALL mpupdatei(qscalarovr,nscalarmax)
    CALL mpupdatei(qscalarcol1,nscalarmax)
    CALL mpupdatei(qscalarcol2,nscalarmax)
    CALL mpupdatei(qscalarprio,nscalarmax)
    CALL mpupdatei(qscalarhlf,nscalarmax)
    CALL mpupdatei(qscalarzro,nscalarmax)
    CALL mpupdatei(qscalarsty,nscalarmax)

  CALL mpupdatei(qwplot,1)
  CALL mpupdater(qwinc,1)
  CALL mpupdater(qwminc,1)
  CALL mpupdater(qwmaxc,1)
  CALL mpupdatei(qwovr,1)
  CALL mpupdatei(qwcol1,1)
  CALL mpupdatei(qwcol2,1)
  CALL mpupdatei(qwprio,1)

  CALL mpupdatei(qwhlf,1)
  CALL mpupdatei(qwzro,1)
  CALL mpupdatei(qwsty,1)

  CALL mpupdatei(qtplot,1)
  CALL mpupdater(qtinc,1)
  CALL mpupdater(qtminc,1)
  CALL mpupdater(qtmaxc,1)
  CALL mpupdatei(qtovr,1)
  CALL mpupdatei(qtcol1,1)
  CALL mpupdatei(qtcol2,1)
  CALL mpupdatei(qtprio,1)

  CALL mpupdatei(qthlf,1)
  CALL mpupdatei(qtzro,1)
  CALL mpupdatei(qtsty,1)

  IF(myproc == 0) THEN
    READ(unum,sclrplt_cntl3, ERR=100)
    WRITE(6,'(a)')'Namelist sclrplt_cntl3 was successfully read.'
  END IF
  CALL mpupdatei(kmhplt,1)
  CALL mpupdater(kmhinc,1)
  CALL mpupdater(kmhminc,1)
  CALL mpupdater(kmhmaxc,1)
  CALL mpupdatei(kmhovr,1)
  CALL mpupdatei(kmhcol1,1)
  CALL mpupdatei(kmhcol2,1)
  CALL mpupdatei(kmhprio,1)
  CALL mpupdatei(kmhhlf,1)
  CALL mpupdatei(kmhzro,1)
  CALL mpupdatei(kmhsty,1)

  CALL mpupdatei(kmvplt,1)
  CALL mpupdater(kmvinc,1)
  CALL mpupdater(kmvminc,1)
  CALL mpupdater(kmvmaxc,1)
  CALL mpupdatei(kmvovr,1)
  CALL mpupdatei(kmvcol1,1)
  CALL mpupdatei(kmvcol2,1)
  CALL mpupdatei(kmvprio,1)
  CALL mpupdatei(kmvhlf,1)
  CALL mpupdatei(kmvzro,1)
  CALL mpupdatei(kmvsty,1)

  CALL mpupdatei(tkeplt,1)
  CALL mpupdater(tkeinc,1)
  CALL mpupdater(tkeminc,1)
  CALL mpupdater(tkemaxc,1)
  CALL mpupdatei(tkeovr,1)
  CALL mpupdatei(tkecol1,1)
  CALL mpupdatei(tkecol2,1)

  CALL mpupdatei(tkeprio,1)
  CALL mpupdatei(tkehlf,1)
  CALL mpupdatei(tkezro,1)
  CALL mpupdatei(tkesty,1)

  CALL mpupdatei(rhplot,1)
  CALL mpupdater(rhinc,1)
  CALL mpupdater(rhminc,1)
  CALL mpupdater(rhmaxc,1)
  CALL mpupdatei(rhovr,1)
  CALL mpupdatei(rhcol1,1)
  CALL mpupdatei(rhcol2,1)
  CALL mpupdatei(rhprio,1)
  CALL mpupdatei(rhhlf,1)
  CALL mpupdatei(rhzro,1)
  CALL mpupdatei(rhsty,1)

  CALL mpupdatei(tdplot,1)
  CALL mpupdater(tdinc,1)
  CALL mpupdater(tdminc,1)
  CALL mpupdater(tdmaxc,1)
  CALL mpupdatei(tdovr,1)
  CALL mpupdatei(tdcol1,1)
  CALL mpupdatei(tdcol2,1)
  CALL mpupdatei(tdprio,1)
  CALL mpupdatec(tdunits,1)
  CALL mpupdatei(tdhlf,1)
  CALL mpupdatei(tdzro,1)
  CALL mpupdatei(tdsty,1)

  CALL mpupdatei(rfopt,1)
  CALL mpupdatei(dualpol,1)
  CALL mpupdatei(rsadir,256)
  CALL mpupdatei(wavelen,1)
  CALL mpupdatei(graupel_ON,1)
  CALL mpupdatei(hail_ON,1)

!-----------------------------------------------------------------------
!  Read scattering amplitude calculated employing T-matrix application.
!-----------------------------------------------------------------------
  if(dualpol == 2) then
!    CALL read_table (rsadir)
  endif

  CALL mpupdatei(rfplot,1)
  CALL mpupdater(rfinc,1)
  CALL mpupdater(rfminc,1)
  CALL mpupdater(rfmaxc,1)
  CALL mpupdatei(rfovr,1)
  CALL mpupdatei(rfcol1,1)
  CALL mpupdatei(rfcol2,1)
  CALL mpupdatei(rfprio,1)
  CALL mpupdatei(rfhlf,1)
  CALL mpupdatei(rfzro,1)
  CALL mpupdatei(rfsty,1)

  CALL mpupdatei(rfcplt,1)
  CALL mpupdater(rfcinc,1)
  CALL mpupdater(rfcminc,1)
  CALL mpupdater(rfcmaxc,1)
  CALL mpupdatei(rfcovr,1)
  CALL mpupdatei(rfccol1,1)
  CALL mpupdatei(rfccol2,1)
  CALL mpupdatei(rfcprio,1)
  CALL mpupdatei(rfchlf,1)
  CALL mpupdatei(rfczro,1)
  CALL mpupdatei(rfcsty,1)

  CALL mpupdatei(pteplt,1)
  CALL mpupdater(pteinc,1)
  CALL mpupdater(pteminc,1)
  CALL mpupdater(ptemaxc,1)
  CALL mpupdatei(pteovr,1)
  CALL mpupdatei(ptecol1,1)
  CALL mpupdatei(ptecol2,1)
  CALL mpupdatei(pteprio,1)
  CALL mpupdatei(ptehlf,1)
  CALL mpupdatei(ptezro,1)
  CALL mpupdatei(ptesty,1)

  CALL mpupdatei(zdrplt,1)
  CALL mpupdater(zdrinc,1)
  CALL mpupdater(zdrminc,1)
  CALL mpupdater(zdrmaxc,1)
  CALL mpupdatei(zdrovr,1)
  CALL mpupdatei(zdrcol1,1)
  CALL mpupdatei(zdrcol2,1)
  CALL mpupdatei(zdrprio,1)
  CALL mpupdatei(zdrhlf,1)
  CALL mpupdatei(zdrzro,1)
  CALL mpupdatei(zdrsty,1)

  CALL mpupdatei(kdpplt,1)
  CALL mpupdater(kdpinc,1)
  CALL mpupdater(kdpminc,1)
  CALL mpupdater(kdpmaxc,1)
  CALL mpupdatei(kdpovr,1)
  CALL mpupdatei(kdpcol1,1)
  CALL mpupdatei(kdpcol2,1)
  CALL mpupdatei(kdpprio,1)
  CALL mpupdatei(kdphlf,1)
  CALL mpupdatei(kdpzro,1)
  CALL mpupdatei(kdpsty,1)

  CALL mpupdatei(zdpplt,1)
  CALL mpupdater(zdpinc,1)
  CALL mpupdater(zdpminc,1)
  CALL mpupdater(zdpmaxc,1)
  CALL mpupdatei(zdpovr,1)
  CALL mpupdatei(zdpcol1,1)
  CALL mpupdatei(zdpcol2,1)
  CALL mpupdatei(zdpprio,1)
  CALL mpupdatei(zdphlf,1)
  CALL mpupdatei(zdpzro,1)
  CALL mpupdatei(zdpsty,1)

  CALL mpupdatei(rhvplt,1)
  CALL mpupdater(rhvinc,1)
  CALL mpupdater(rhvminc,1)
  CALL mpupdater(rhvmaxc,1)
  CALL mpupdatei(rhvovr,1)
  CALL mpupdatei(rhvcol1,1)
  CALL mpupdatei(rhvcol2,1)
  CALL mpupdatei(rhvprio,1)
  CALL mpupdatei(rhvhlf,1)
  CALL mpupdatei(rhvzro,1)
  CALL mpupdatei(rhvsty,1)

  IF(myproc == 0) THEN
    READ(unum,sclrplt_cntl_prt1, ERR=100)
    WRITE(6,'(a)')'Namelist sclrplt_cntl_prt1 was successfully read.'
  END IF
  CALL mpupdatei(upplot,1)
  CALL mpupdater(upinc,1)
  CALL mpupdater(upminc,1)
  CALL mpupdater(upmaxc,1)
  CALL mpupdatei(upovr,1)
  CALL mpupdatei(upcol1,1)
  CALL mpupdatei(upcol2,1)
  CALL mpupdatei(upprio,1)
  CALL mpupdatei(uphlf,1)
  CALL mpupdatei(upzro,1)
  CALL mpupdatei(upsty,1)

  CALL mpupdatei(vpplot,1)
  CALL mpupdater(vpinc,1)
  CALL mpupdater(vpminc,1)
  CALL mpupdater(vpmaxc,1)
  CALL mpupdatei(vpovr,1)
  CALL mpupdatei(vpcol1,1)
  CALL mpupdatei(vpcol2,1)
  CALL mpupdatei(vpprio,1)
  CALL mpupdatei(vphlf,1)
  CALL mpupdatei(vpzro,1)
  CALL mpupdatei(vpsty,1)

  CALL mpupdatei(wpplot,1)
  CALL mpupdater(wpinc,1)
  CALL mpupdater(wpminc,1)
  CALL mpupdater(wpmaxc,1)
  CALL mpupdatei(wpovr,1)
  CALL mpupdatei(wpcol1,1)
  CALL mpupdatei(wpcol2,1)
  CALL mpupdatei(wpprio,1)
  CALL mpupdatei(wphlf,1)
  CALL mpupdatei(wpzro,1)
  CALL mpupdatei(wpsty,1)

  CALL mpupdatei(ptpplt,1)
  CALL mpupdater(ptpinc,1)
  CALL mpupdater(ptpminc,1)
  CALL mpupdater(ptpmaxc,1)
  CALL mpupdatei(ptpovr,1)
  CALL mpupdatei(ptpcol1,1)
  CALL mpupdatei(ptpcol2,1)
  CALL mpupdatei(ptpprio,1)
  CALL mpupdatei(ptphlf,1)
  CALL mpupdatei(ptpzro,1)
  CALL mpupdatei(ptpsty,1)

  CALL mpupdatei(ppplot,1)
  CALL mpupdater(ppinc,1)
  CALL mpupdater(ppminc,1)
  CALL mpupdater(ppmaxc,1)
  CALL mpupdatei(ppovr,1)
  CALL mpupdatei(ppcol1,1)
  CALL mpupdatei(ppcol2,1)
  CALL mpupdatei(ppprio,1)
  CALL mpupdatei(pphlf,1)
  CALL mpupdatei(ppzro,1)
  CALL mpupdatei(ppsty,1)

  CALL mpupdatei(qvpplt,1)
  CALL mpupdater(qvpinc,1)
  CALL mpupdater(qvpminc,1)
  CALL mpupdater(qvpmaxc,1)
  CALL mpupdatei(qvpovr,1)
  CALL mpupdatei(qvpcol1,1)
  CALL mpupdatei(qvpcol2,1)
  CALL mpupdatei(qvpprio,1)
  CALL mpupdatei(qvphlf,1)
  CALL mpupdatei(qvpzro,1)
  CALL mpupdatei(qvpsty,1)

  CALL mpupdatei(vorpplt,1)
  CALL mpupdater(vorpinc,1)
  CALL mpupdater(vorpminc,1)
  CALL mpupdater(vorpmaxc,1)
  CALL mpupdatei(vorpovr,1)
  CALL mpupdatei(vorpcol1,1)
  CALL mpupdatei(vorpcol2,1)
  CALL mpupdatei(vorphlf,1)
  CALL mpupdatei(vorpprio,1)
  CALL mpupdatei(vorpzro,1)
  CALL mpupdatei(vorpsty,1)

  CALL mpupdatei(divpplt,1)
  CALL mpupdater(divpinc,1)
  CALL mpupdater(divpminc,1)
  CALL mpupdater(divpmaxc,1)
  CALL mpupdatei(divpovr,1)
  CALL mpupdatei(divpcol1,1)
  CALL mpupdatei(divpcol2,1)
  CALL mpupdatei(divphlf,1)
  CALL mpupdatei(divpprio,1)
  CALL mpupdatei(divpzro,1)
  CALL mpupdatei(divpsty,1)

  CALL mpupdatei(divqplt,1)
  CALL mpupdater(divqinc,1)
  CALL mpupdater(divqminc,1)
  CALL mpupdater(divqmaxc,1)
  CALL mpupdatei(divqovr,1)
  CALL mpupdatei(divqcol1,1)
  CALL mpupdatei(divqcol2,1)
  CALL mpupdatei(divqhlf,1)
  CALL mpupdatei(divqprio,1)
  CALL mpupdatei(divqzro,1)
  CALL mpupdatei(divqsty,1)

  IF(myproc == 0) THEN
    READ(unum,sclrplt_cntl_prt2, ERR=100)
    WRITE(6,'(a)')'Namelist sclrplt_cntl_prt2 was successfully read.'
  END IF
  CALL mpupdatei(gricplt,1)
  CALL mpupdater(gricinc,1)
  CALL mpupdater(gricminc,1)
  CALL mpupdater(gricmaxc,1)
  CALL mpupdatei(gricovr,1)
  CALL mpupdatei(griccol1,1)
  CALL mpupdatei(griccol2,1)
  CALL mpupdatei(grichlf,1)
  CALL mpupdatei(gricprio,1)
  CALL mpupdatei(griczro,1)
  CALL mpupdatei(gricsty,1)

  CALL mpupdatei(avorplt,1)
  CALL mpupdater(avorinc,1)
  CALL mpupdater(avorminc,1)
  CALL mpupdater(avormaxc,1)
  CALL mpupdatei(avorovr,1)
  CALL mpupdatei(avorcol1,1)
  CALL mpupdatei(avorcol2,1)
  CALL mpupdatei(avorhlf,1)
  CALL mpupdatei(avorprio,1)
  CALL mpupdatei(avorzro,1)
  CALL mpupdatei(avorsty,1)

  CALL mpupdatei(rhiplot,1)
  CALL mpupdater(rhiinc,1)
  CALL mpupdater(rhiminc,1)
  CALL mpupdater(rhimaxc,1)
  CALL mpupdatei(rhiovr,1)
  CALL mpupdatei(rhicol1,1)
  CALL mpupdatei(rhicol2,1)
  CALL mpupdatei(rhiprio,1)
  CALL mpupdatei(rhihlf,1)
  CALL mpupdatei(rhizro,1)
  CALL mpupdatei(rhisty,1)

  IF(myproc == 0) THEN
    READ(unum,vctrplt_cntl, ERR=100)
    WRITE(6,'(a)')'Namelist vctrplt_cntl was successfully read.'
  END IF
  CALL mpupdatei(istride,1)
  CALL mpupdatei(jstride,1)
  CALL mpupdatei(kstride,1)

  CALL mpupdatei(vtrplt,1)
  CALL mpupdater(vtrunit,1)
  CALL mpupdatei(vtrovr,1)
  CALL mpupdatei(vtrcol1,1)
  CALL mpupdatei(vtrcol2,1)
  CALL mpupdatei(vtrprio,1)
  CALL mpupdatei(vtrunits,1)
  CALL mpupdatei(vtrtype,1)

  CALL mpupdatei(vtpplt,1)
  CALL mpupdater(vtpunit,1)
  CALL mpupdatei(vtpovr,1)
  CALL mpupdatei(vtpcol1,1)
  CALL mpupdatei(vtpcol2,1)
  CALL mpupdatei(vtpprio,1)
  CALL mpupdatei(vtpunits,1)
  CALL mpupdatei(vtptype,1)

  CALL mpupdatei(xuvplt,1)
  CALL mpupdater(xuvunit,1)
  CALL mpupdatei(xuvovr,1)
  CALL mpupdatei(xuvcol1,1)
  CALL mpupdatei(xuvcol2,1)
  CALL mpupdatei(xuvprio,1)
  CALL mpupdatei(xuvunits,1)
  CALL mpupdatei(xuvtype,1)

  CALL mpupdatei(strmplt,1)
  CALL mpupdater(strmunit,1)
  CALL mpupdatei(strmovr,1)
  CALL mpupdatei(strmcol1,1)
  CALL mpupdatei(strmcol2,1)
  CALL mpupdatei(strmprio,1)
  CALL mpupdatei(strmunits,1)
  CALL mpupdatei(strmtype,1)

  CALL mpupdatei(vagplt,1)
  CALL mpupdater(vagunit,1)
  CALL mpupdatei(vagovr,1)
  CALL mpupdatei(vagcol1,1)
  CALL mpupdatei(vagcol2,1)
  CALL mpupdatei(vagprio,1)
  CALL mpupdatei(vagunits,1)
  CALL mpupdatei(vagtype,1)

  IF(myproc == 0) THEN
    READ(unum,strmplt_cntl,ERR=100)
    WRITE(6,'(a)')'Namelist strmplt_cntl was successfully read.'
!    IF ( nprocs > 1 .AND. vtrstrm /= 0 ) THEN
!        WRITE(6,'(a)')                                                         &
!        '*** Streamlines won''t work in an MPI job.  Turning off.  ***'
!        vtrstrm = 0
!    END IF
  END IF
  CALL mpupdatei(vtrstrm,1)
  CALL mpupdatei(vtrstmovr,1)
  CALL mpupdatei(vtrstmcol1,1)
  CALL mpupdatei(vtrstmcol2,1)
  CALL mpupdatei(vtrstmprio,1)

  CALL mpupdatei(vtpstrm,1)
  CALL mpupdatei(vtpstmovr,1)
  CALL mpupdatei(vtpstmcol1,1)
  CALL mpupdatei(vtpstmcol2,1)
  CALL mpupdatei(vtpstmprio,1)

!
!-----------------------------------------------------------------------
!
!  Input control parameters for 2-d surface feild plotting
!
!-----------------------------------------------------------------------
!
  IF(myproc == 0) THEN
    READ(unum,sfc_plot1,ERR=100)
    WRITE(6,'(a)')'Namelist sfc_plot1 was successfully read.'
  END IF
  CALL mpupdatei(trnplt,1)
  CALL mpupdater(trninc,1)
  CALL mpupdater(trnminc,1)
  CALL mpupdater(trnmaxc,1)
  CALL mpupdatei(trnovr,1)
  CALL mpupdatei(trncol1,1)
  CALL mpupdatei(trncol2,1)
  CALL mpupdatei(trnprio,1)
  CALL mpupdatei(trnhlf,1)
  CALL mpupdatei(trnzro,1)
  CALL mpupdatei(trnsty,1)

  CALL mpupdatei(wetcanplt,1)
  CALL mpupdater(wcpinc,1)
  CALL mpupdater(wcpminc,1)
  CALL mpupdater(wcpmaxc,1)
  CALL mpupdatei(wcpovr,1)
  CALL mpupdatei(wcpcol1,1)
  CALL mpupdatei(wcpcol2,1)
  CALL mpupdatei(wcpprio,1)
  CALL mpupdatei(wcphlf,1)
  CALL mpupdatei(wcpzro,1)
  CALL mpupdatei(wcpsty,1)

  CALL mpupdatei(raincplt,1)
  CALL mpupdater(raincinc,1)
  CALL mpupdater(raincminc,1)
  CALL mpupdater(raincmaxc,1)
  CALL mpupdatei(racovr,1)
  CALL mpupdatei(raccol1,1)
  CALL mpupdatei(raccol2,1)
  CALL mpupdatei(rachlf,1)
  CALL mpupdatei(racprio,1)
  CALL mpupdatei(raczro,1)
  CALL mpupdatei(racsty,1)
  CALL mpupdatei(racunit,1)

  CALL mpupdatei(raingplt,1)
  CALL mpupdater(rainginc,1)
  CALL mpupdater(raingminc,1)
  CALL mpupdater(raingmaxc,1)
  CALL mpupdatei(ragovr,1)
  CALL mpupdatei(ragcol1,1)
  CALL mpupdatei(ragcol2,1)
  CALL mpupdatei(raghlf,1)
  CALL mpupdatei(ragprio,1)
  CALL mpupdatei(ragzro,1)
  CALL mpupdatei(ragsty,1)
  CALL mpupdatei(ragunit,1)

  CALL mpupdatei(raintplt,1)
  CALL mpupdater(raintinc,1)
  CALL mpupdater(raintminc,1)
  CALL mpupdater(raintmaxc,1)
  CALL mpupdatei(ratovr,1)
  CALL mpupdatei(ratcol1,1)
  CALL mpupdatei(ratcol2,1)
  CALL mpupdatei(rathlf,1)
  CALL mpupdatei(ratprio,1)
  CALL mpupdatei(ratzro,1)
  CALL mpupdatei(ratsty,1)
  CALL mpupdatei(ratunit,1)

  CALL mpupdatei(rainicplt,1)
  CALL mpupdater(rainicinc,1)
  CALL mpupdater(rainicminc,1)
  CALL mpupdater(rainicmaxc,1)
  CALL mpupdatei(raicovr,1)
  CALL mpupdatei(raiccol1,1)
  CALL mpupdatei(raiccol2,1)
  CALL mpupdatei(raichlf,1)
  CALL mpupdatei(raicprio,1)
  CALL mpupdatei(raiczro,1)
  CALL mpupdatei(raicsty,1)
  CALL mpupdatei(raicunit,1)

  CALL mpupdatei(rainigplt,1)
  CALL mpupdater(rainiginc,1)
  CALL mpupdater(rainigminc,1)
  CALL mpupdater(rainigmaxc,1)
  CALL mpupdatei(raigovr,1)
  CALL mpupdatei(raigcol1,1)
  CALL mpupdatei(raigcol2,1)
  CALL mpupdatei(raighlf,1)
  CALL mpupdatei(raigprio,1)
  CALL mpupdatei(raigzro,1)
  CALL mpupdatei(raigsty,1)
  CALL mpupdatei(raigunit,1)

  CALL mpupdatei(rainitplt,1)
  CALL mpupdater(rainitinc,1)
  CALL mpupdater(rainitminc,1)
  CALL mpupdater(rainitmaxc,1)
  CALL mpupdatei(raitovr,1)
  CALL mpupdatei(raitcol1,1)
  CALL mpupdatei(raitcol2,1)
  CALL mpupdatei(raithlf,1)
  CALL mpupdatei(raitprio,1)
  CALL mpupdatei(raitzro,1)
  CALL mpupdatei(raitsty,1)
  CALL mpupdatei(raitunit,1)

  IF(myproc == 0) THEN
    READ(unum,soil_plot,ERR=100)
    WRITE(6,'(a)')'Namelist soil_plot was successfully read.'
  END IF

  CALL mpupdatei(tsoilplt,1)
  CALL mpupdater(tsoilinc,1)
  CALL mpupdater(tsoilminc,1)
  CALL mpupdater(tsoilmaxc,1)
  CALL mpupdatei(tsoilovr,1)

  CALL mpupdatei(tsoilcol1,1)
  CALL mpupdatei(tsoilcol2,1)
  CALL mpupdatei(tsoilhlf,1)
  CALL mpupdatei(tsoilprio,1)
  CALL mpupdatei(tsoilzro,1)

  CALL mpupdatei(qsoilplt,1)
  CALL mpupdater(qsoilinc,1)
  CALL mpupdater(qsoilminc,1)
  CALL mpupdater(qsoilmaxc,1)
  CALL mpupdatei(qsoilovr,1)

  CALL mpupdatei(qsoilcol1,1)
  CALL mpupdatei(qsoilcol2,1)
  CALL mpupdatei(qsoilhlf,1)
  CALL mpupdatei(qsoilprio,1)
  CALL mpupdatei(qsoilzro,1)

  IF(myproc == 0) THEN
    READ(unum,sfc_plot2,ERR=100)
    WRITE(6,'(a)')'Namelist sfc_plot2 was successfully read.'
  END IF
  CALL mpupdatei(pslplt,1)
  CALL mpupdater(pslinc,1)
  CALL mpupdater(pslminc,1)
  CALL mpupdater(pslmaxc,1)
  CALL mpupdatei(pslovr,1)
  CALL mpupdatei(pslcol1,1)
  CALL mpupdatei(pslcol2,1)
  CALL mpupdatei(pslprio,1)
  CALL mpupdatei(pslhlf,1)
  CALL mpupdatei(pslzro,1)
  CALL mpupdatei(pslsty,1)

  CALL mpupdatei(capeplt,1)
  CALL mpupdater(capeinc,1)
  CALL mpupdater(capeminc,1)
  CALL mpupdater(capemaxc,1)
  CALL mpupdatei(capovr,1)
  CALL mpupdatei(capcol1,1)
  CALL mpupdatei(capcol2,1)
  CALL mpupdatei(capprio,1)
  CALL mpupdatei(caphlf,1)
  CALL mpupdatei(capzro,1)
  CALL mpupdatei(capsty,1)

  CALL mpupdatei(cinplt,1)
  CALL mpupdater(cininc,1)
  CALL mpupdater(cinminc,1)
  CALL mpupdater(cinmaxc,1)
  CALL mpupdatei(cinovr,1)
  CALL mpupdatei(cincol1,1)
  CALL mpupdatei(cincol2,1)
  CALL mpupdatei(cinprio,1)
  CALL mpupdatei(cinhlf,1)
  CALL mpupdatei(cinzro,1)
  CALL mpupdatei(cinsty,1)

  CALL mpupdatei(thetplt,1)
  CALL mpupdater(thetinc,1)
  CALL mpupdater(thetminc,1)
  CALL mpupdater(thetmaxc,1)
  CALL mpupdatei(theovr,1)
  CALL mpupdatei(thecol1,1)
  CALL mpupdatei(thecol2,1)
  CALL mpupdatei(theprio,1)
  CALL mpupdatei(thehlf,1)
  CALL mpupdatei(thezro,1)
  CALL mpupdatei(thesty,1)

  CALL mpupdatei(heliplt,1)
  CALL mpupdater(heliinc,1)
  CALL mpupdater(heliminc,1)
  CALL mpupdater(helimaxc,1)
  CALL mpupdatei(helovr,1)
  CALL mpupdatei(helcol1,1)
  CALL mpupdatei(helcol2,1)
  CALL mpupdatei(helprio,1)
  CALL mpupdatei(helhlf,1)
  CALL mpupdatei(helzro,1)
  CALL mpupdatei(helsty,1)

  CALL mpupdatei(uhplt,1)
  CALL mpupdater(uhinc,1)
  CALL mpupdater(uhminc,1)
  CALL mpupdater(uhmaxc,1)
  CALL mpupdatei(uhovr,1)
  CALL mpupdatei(uhcol1,1)
  CALL mpupdatei(uhcol2,1)
  CALL mpupdatei(uhprio,1)
  CALL mpupdatei(uhhlf,1)
  CALL mpupdatei(uhzro,1)
  CALL mpupdatei(uhsty,1)
  CALL mpupdater(uhmnhgt,1)
  CALL mpupdater(uhmxhgt,1)

  CALL mpupdatei(brnplt,1)
  CALL mpupdater(brninc,1)
  CALL mpupdater(brnminc,1)
  CALL mpupdater(brnmaxc,1)
  CALL mpupdatei(brnovr,1)
  CALL mpupdatei(brncol1,1)
  CALL mpupdatei(brncol2,1)
  CALL mpupdatei(brnprio,1)
  CALL mpupdatei(brnhlf,1)
  CALL mpupdatei(brnzro,1)
  CALL mpupdatei(brnsty,1)

  CALL mpupdatei(brnuplt,1)
  CALL mpupdater(brnuinc,1)
  CALL mpupdater(bruminc,1)
  CALL mpupdater(brumaxc,1)
  CALL mpupdatei(brnuovr,1)
  CALL mpupdatei(brnucol1,1)
  CALL mpupdatei(brnucol2,1)
  CALL mpupdatei(brnuhlf,1)
  CALL mpupdatei(brnuzro,1)
  CALL mpupdatei(brnusty,1)
  CALL mpupdatei(bruprio,1)

  CALL mpupdatei(srlfplt,1)
  CALL mpupdater(srlfinc,1)
  CALL mpupdater(srlminc,1)
  CALL mpupdater(srlmaxc,1)
  CALL mpupdatei(srlfovr,1)
  CALL mpupdatei(srlfcol1,1)
  CALL mpupdatei(srlfcol2,1)
  CALL mpupdatei(srlfhlf,1)
  CALL mpupdatei(srlfzro,1)
  CALL mpupdatei(srlfsty,1)
  CALL mpupdatei(srlprio,1)

  CALL mpupdatei(srmfplt,1)
  CALL mpupdater(srmfinc,1)
  CALL mpupdater(srmminc,1)
  CALL mpupdater(srmmaxc,1)
  CALL mpupdatei(srmfovr,1)
  CALL mpupdatei(srmfcol1,1)
  CALL mpupdatei(srmfcol2,1)
  CALL mpupdatei(srmfhlf,1)
  CALL mpupdatei(srmfzro,1)
  CALL mpupdatei(srmfsty,1)
  CALL mpupdatei(srmprio,1)

  IF(myproc == 0) THEN
    READ(unum,sfc_plot3,ERR=100)
    WRITE(6,'(a)')'Namelist sfc_plot3 was successfully read.'
  END IF

  CALL mpupdatei(liplt,1)
  CALL mpupdater(liinc,1)
  CALL mpupdater(liminc,1)
  CALL mpupdater(limaxc,1)
  CALL mpupdatei(liovr,1)
  CALL mpupdatei(licol1,1)
  CALL mpupdatei(licol2,1)
  CALL mpupdatei(liprio,1)
  CALL mpupdatei(lihlf,1)
  CALL mpupdatei(lizro,1)
  CALL mpupdatei(listy,1)

  CALL mpupdatei(capsplt,1)
  CALL mpupdater(capsinc,1)
  CALL mpupdater(capsminc,1)
  CALL mpupdater(capsmaxc,1)
  CALL mpupdatei(capsovr,1)
  CALL mpupdatei(capscol1,1)
  CALL mpupdatei(capscol2,1)
  CALL mpupdatei(capshlf,1)
  CALL mpupdatei(capszro,1)
  CALL mpupdatei(capssty,1)
  CALL mpupdatei(capsprio,1)

  CALL mpupdatei(blcoplt,1)
  CALL mpupdater(blcoinc,1)
  CALL mpupdater(blcominc,1)
  CALL mpupdater(blcomaxc,1)
  CALL mpupdatei(blcoovr,1)
  CALL mpupdatei(blcocol1,1)
  CALL mpupdatei(blcocol2,1)
  CALL mpupdatei(blcohlf,1)
  CALL mpupdatei(blcozro,1)
  CALL mpupdatei(blcosty,1)
  CALL mpupdatei(blcoprio,1)

  CALL mpupdatei(viqcplt,1)
  CALL mpupdater(viqcinc,1)
  CALL mpupdater(viqcminc,1)
  CALL mpupdater(viqcmaxc,1)
  CALL mpupdatei(viqcovr,1)
  CALL mpupdatei(viqccol1,1)
  CALL mpupdatei(viqccol2,1)
  CALL mpupdatei(viqchlf,1)
  CALL mpupdatei(viqczro,1)
  CALL mpupdatei(viqcsty,1)
  CALL mpupdatei(viqcprio,1)

  CALL mpupdatei(viqiplt,1)
  CALL mpupdater(viqiinc,1)
  CALL mpupdater(viqiminc,1)
  CALL mpupdater(viqimaxc,1)
  CALL mpupdatei(viqiovr,1)
  CALL mpupdatei(viqicol1,1)
  CALL mpupdatei(viqicol2,1)
  CALL mpupdatei(viqihlf,1)
  CALL mpupdatei(viqizro,1)
  CALL mpupdatei(viqisty,1)
  CALL mpupdatei(viqiprio,1)

  CALL mpupdatei(viqrplt,1)
  CALL mpupdater(viqrinc,1)
  CALL mpupdater(viqrminc,1)
  CALL mpupdater(viqrmaxc,1)
  CALL mpupdatei(viqrovr,1)
  CALL mpupdatei(viqrcol1,1)
  CALL mpupdatei(viqrcol2,1)
  CALL mpupdatei(viqrhlf,1)
  CALL mpupdatei(viqrzro,1)
  CALL mpupdatei(viqrsty,1)
  CALL mpupdatei(viqrprio,1)

  CALL mpupdatei(viqsplt,1)
  CALL mpupdater(viqsinc,1)
  CALL mpupdater(viqsminc,1)
  CALL mpupdater(viqsmaxc,1)
  CALL mpupdatei(viqsovr,1)
  CALL mpupdatei(viqscol1,1)
  CALL mpupdatei(viqscol2,1)
  CALL mpupdatei(viqshlf,1)
  CALL mpupdatei(viqszro,1)
  CALL mpupdatei(viqssty,1)
  CALL mpupdatei(viqsprio,1)

  CALL mpupdatei(viqhplt,1)
  CALL mpupdater(viqhinc,1)
  CALL mpupdater(viqhminc,1)
  CALL mpupdater(viqhmaxc,1)
  CALL mpupdatei(viqhovr,1)
  CALL mpupdatei(viqhcol1,1)
  CALL mpupdatei(viqhcol2,1)
  CALL mpupdatei(viqhhlf,1)
  CALL mpupdatei(viqhzro,1)
  CALL mpupdatei(viqhsty,1)
  CALL mpupdatei(viqhprio,1)

  CALL mpupdatei(vilplt,1)
  CALL mpupdater(vilinc,1)
  CALL mpupdater(vilminc,1)
  CALL mpupdater(vilmaxc,1)
  CALL mpupdatei(vilovr,1)
  CALL mpupdatei(vilcol1,1)
  CALL mpupdatei(vilcol2,1)
  CALL mpupdatei(vilhlf,1)
  CALL mpupdatei(vilzro,1)
  CALL mpupdatei(vilsty,1)
  CALL mpupdatei(vilprio,1)

  IF(myproc == 0) THEN
    READ(unum,sfc_plot4,ERR=100)
    WRITE(6,'(a)')'Namelist sfc_plot4 was successfully read.'
  END IF
  CALL mpupdatei(viiplt,1)
  CALL mpupdater(viiinc,1)
  CALL mpupdater(viiminc,1)
  CALL mpupdater(viimaxc,1)
  CALL mpupdatei(viiovr,1)
  CALL mpupdatei(viicol1,1)
  CALL mpupdatei(viicol2,1)
  CALL mpupdatei(viihlf,1)
  CALL mpupdatei(viizro,1)
  CALL mpupdatei(viisty,1)
  CALL mpupdatei(viiprio,1)

  CALL mpupdatei(vicplt,1)
  CALL mpupdater(vicinc,1)
  CALL mpupdater(vicminc,1)
  CALL mpupdater(vicmaxc,1)
  CALL mpupdatei(vicovr,1)
  CALL mpupdatei(viccol1,1)
  CALL mpupdatei(viccol2,1)
  CALL mpupdatei(vichlf,1)
  CALL mpupdatei(viczro,1)
  CALL mpupdatei(vicsty,1)
  CALL mpupdatei(vicprio,1)

  CALL mpupdatei(ctcplt,1)
  CALL mpupdater(ctcinc,1)
  CALL mpupdater(ctcminc,1)
  CALL mpupdater(ctcmaxc,1)
  CALL mpupdatei(ctcovr,1)
  CALL mpupdatei(ctccol1,1)
  CALL mpupdatei(ctccol2,1)
  CALL mpupdatei(ctchlf,1)
  CALL mpupdatei(ctczro,1)
  CALL mpupdatei(ctcsty,1)
  CALL mpupdatei(ctcprio,1)

  CALL mpupdatei(vitplt,1)
  CALL mpupdater(vitinc,1)
  CALL mpupdater(vitminc,1)
  CALL mpupdater(vitmaxc,1)
  CALL mpupdatei(vitovr,1)
  CALL mpupdatei(vitcol1,1)
  CALL mpupdatei(vitcol2,1)
  CALL mpupdatei(vithlf,1)
  CALL mpupdatei(vitzro,1)
  CALL mpupdatei(vitsty,1)
  CALL mpupdatei(vitprio,1)

  CALL mpupdatei(pwplt,1)
  CALL mpupdater(pwinc,1)
  CALL mpupdater(pwminc,1)
  CALL mpupdater(pwmaxc,1)
  CALL mpupdatei(pwovr,1)
  CALL mpupdatei(pwcol1,1)
  CALL mpupdatei(pwcol2,1)
  CALL mpupdatei(pwhlf,1)
  CALL mpupdatei(pwzro,1)
  CALL mpupdatei(pwsty,1)
  CALL mpupdatei(pwprio,1)

  CALL mpupdatei(tprplt,1)
  CALL mpupdater(tprinc,1)
  CALL mpupdater(tprminc,1)
  CALL mpupdater(tprmaxc,1)
  CALL mpupdatei(tprovr,1)
  CALL mpupdatei(tprcol1,1)
  CALL mpupdatei(tprcol2,1)
  CALL mpupdatei(tprhlf,1)
  CALL mpupdatei(tprzro,1)
  CALL mpupdatei(tprsty,1)
  CALL mpupdatei(tprprio,1)
  CALL mpupdatei(tprunits,1)

  CALL mpupdatei(gprplt,1)
  CALL mpupdater(gprinc,1)
  CALL mpupdater(gprminc,1)
  CALL mpupdater(gprmaxc,1)
  CALL mpupdatei(gprovr,1)
  CALL mpupdatei(gprcol1,1)
  CALL mpupdatei(gprcol2,1)
  CALL mpupdatei(gprhlf,1)
  CALL mpupdatei(gprzro,1)
  CALL mpupdatei(gprsty,1)
  CALL mpupdatei(gprprio,1)
  CALL mpupdatei(gprunits,1)

  CALL mpupdatei(cprplt,1)
  CALL mpupdater(cprinc,1)
  CALL mpupdater(cprminc,1)
  CALL mpupdater(cprmaxc,1)
  CALL mpupdatei(cprovr,1)
  CALL mpupdatei(cprcol1,1)
  CALL mpupdatei(cprcol2,1)
  CALL mpupdatei(cprhlf,1)
  CALL mpupdatei(cprzro,1)
  CALL mpupdatei(cprsty,1)
  CALL mpupdatei(cprprio,1)
  CALL mpupdatei(cprunits,1)

!-----------------------------------------------------------------------
!
!  Input control parameters for 2-d surface characteristics plotting
!
!-----------------------------------------------------------------------
!
  IF(myproc == 0) THEN
    READ(unum,sfc_cha_plot,ERR=100)
    WRITE(6,'(a)')                                                      &
      'Namelist sfc_cha_plot was successfully read.'
  END IF
  CALL mpupdatei(soiltpplt,1)
  CALL mpupdater(soiltpinc,1)
  CALL mpupdater(soiltpminc,1)
  CALL mpupdater(soiltpmaxc,1)
  CALL mpupdatei(styovr,1)
  CALL mpupdatei(stycol1,1)
  CALL mpupdatei(stycol2,1)
  CALL mpupdatei(styhlf,1)
  CALL mpupdatei(styzro,1)
  CALL mpupdatei(stysty,1)
  CALL mpupdatei(styprio,1)
  CALL mpupdatei(soiltpn,1)

  CALL mpupdatei(vegtpplt,1)
  CALL mpupdater(vegtpinc,1)
  CALL mpupdater(vegtpminc,1)
  CALL mpupdater(vegtpmaxc,1)
  CALL mpupdatei(vtyovr,1)
  CALL mpupdatei(vtycol1,1)
  CALL mpupdatei(vtycol2,1)
  CALL mpupdatei(vtyhlf,1)
  CALL mpupdatei(vtyzro,1)
  CALL mpupdatei(vtysty,1)
  CALL mpupdatei(vtyprio,1)

  CALL mpupdatei(laiplt,1)
  CALL mpupdater(laiinc,1)
  CALL mpupdater(laiminc,1)
  CALL mpupdater(laimaxc,1)
  CALL mpupdatei(laiovr,1)
  CALL mpupdatei(laicol1,1)
  CALL mpupdatei(laicol2,1)
  CALL mpupdatei(laiprio,1)
  CALL mpupdatei(laihlf,1)
  CALL mpupdatei(laizro,1)
  CALL mpupdatei(laisty,1)

  CALL mpupdatei(rouplt,1)
  CALL mpupdater(rouinc,1)
  CALL mpupdater(rouminc,1)
  CALL mpupdater(roumaxc,1)
  CALL mpupdatei(rouovr,1)
  CALL mpupdatei(roucol1,1)
  CALL mpupdatei(roucol2,1)
  CALL mpupdatei(rouprio,1)
  CALL mpupdatei(rouhlf,1)
  CALL mpupdatei(rouzro,1)
  CALL mpupdatei(rousty,1)

  CALL mpupdatei(vegplt,1)
  CALL mpupdater(veginc,1)
  CALL mpupdater(vegminc,1)
  CALL mpupdater(vegmaxc,1)
  CALL mpupdatei(vegovr,1)
  CALL mpupdatei(vegcol1,1)
  CALL mpupdatei(vegcol2,1)
  CALL mpupdatei(vegprio,1)
  CALL mpupdatei(veghlf,1)
  CALL mpupdatei(vegzro,1)
  CALL mpupdatei(vegsty,1)

  CALL mpupdatei(snowdplt,1)
  CALL mpupdater(snowdinc,1)
  CALL mpupdater(snowdminc,1)
  CALL mpupdater(snowdmaxc,1)
  CALL mpupdatei(snowdovr,1)
  CALL mpupdatei(snowdcol1,1)
  CALL mpupdatei(snowdcol2,1)
  CALL mpupdatei(snowdprio,1)
  CALL mpupdatei(snowdhlf,1)
  CALL mpupdatei(snowdzro,1)
  CALL mpupdatei(snowdsty,1)

!-----------------------------------------------------------------------
!
!    Input control parameter for uneven contour interval
!
!-----------------------------------------------------------------------

  DO i=1,maxuneva
    setcontvar(i)(1:12) = '            '
    DO k=1,maxunevm
      setconts(k,i) = -9999.
    END DO
  END DO

  IF(myproc == 0) THEN
    READ(unum,setcont_cntl,ERR=100)
    WRITE(6,'(a)')'Namelist setcont_cntl was successfully read.'
  END IF
  CALL mpupdatei(setcontopt,1)
  CALL mpupdatei(setcontnum,1)
  CALL mpupdatec(setcontvar,12*maxuneva)
  CALL mpupdater(setconts,maxunevm*maxuneva)

  filename3d = ' '
  filename2d = ' '
  filenameu2d = ' '
  filenamev2d = ' '
  IF(myproc == 0) THEN
    READ(unum,arbvar_cntl,ERR=100)
    WRITE(6,'(a)')'Namelist arbvar_cntl was successfully read.'
  END IF
  CALL mpupdatei(arbvaropt,1)

  CALL mpupdatei(var3dnum,1)
  CALL mpupdatec(dirname3d,256*maxarbvar)
  CALL mpupdatei(finfmt3d,maxarbvar)
  CALL mpupdatec(var3d,6*maxarbvar)
  CALL mpupdatei(var3dplot,maxarbvar)
  CALL mpupdater(var3dinc,maxarbvar)
  CALL mpupdater(var3dminc,maxarbvar)
  CALL mpupdater(var3dmaxc,maxarbvar)
  CALL mpupdatei(var3dovr,maxarbvar)
  CALL mpupdatei(var3dhlf,maxarbvar)
  CALL mpupdatei(var3dzro,maxarbvar)
  CALL mpupdatei(var3dsty,maxarbvar)
  CALL mpupdatei(var3dcol1,maxarbvar)
  CALL mpupdatei(var3dcol2,maxarbvar)
  CALL mpupdatei(var3dprio,maxarbvar)

  CALL mpupdatei(var2dnum,1)
  CALL mpupdatec(dirname2d,256*maxarbvar)
  CALL mpupdatei(finfmt2d,maxarbvar)
  CALL mpupdatec(var2d,6*maxarbvar)
  CALL mpupdatei(var2dplot,maxarbvar)
  CALL mpupdater(var2dinc,maxarbvar)
  CALL mpupdater(var2dminc,maxarbvar)
  CALL mpupdater(var2dmaxc,maxarbvar)
  CALL mpupdatei(var2dovr,maxarbvar)
  CALL mpupdatei(var2dhlf,maxarbvar)
  CALL mpupdatei(var2dzro,maxarbvar)
  CALL mpupdatei(var2dsty,maxarbvar)
  CALL mpupdatei(var2dcol1,maxarbvar)
  CALL mpupdatei(var2dcol2,maxarbvar)
  CALL mpupdatei(var2dprio,maxarbvar)

  CALL mpupdatei(vtr2dnum,1)
  CALL mpupdatec(diruv2d,256*maxarbvar)
  CALL mpupdatec(vtru2d,6*maxarbvar)
  CALL mpupdatec(vtrv2d,6*maxarbvar)
  CALL mpupdatei(finfmtuv2d,maxarbvar)
  CALL mpupdatei(iastride,maxarbvar)
  CALL mpupdatei(jastride,maxarbvar)
  CALL mpupdatei(vtraplt,maxarbvar)
  CALL mpupdater(vtraunit,maxarbvar)
  CALL mpupdatei(vtraovr,maxarbvar)
  CALL mpupdatei(vtracol1,maxarbvar)
  CALL mpupdatei(vtracol2,maxarbvar)
  CALL mpupdatei(vtraprio,maxarbvar)
  CALL mpupdatei(vtraunits,maxarbvar)
  CALL mpupdatei(vtratype,maxarbvar)
  CALL mpupdatei(magaplt,maxarbvar)
  CALL mpupdater(magainc,maxarbvar)
  CALL mpupdater(magaminc,maxarbvar)
  CALL mpupdater(magamaxc,maxarbvar)
  CALL mpupdatei(magaovr,maxarbvar)
  CALL mpupdatei(magahlf,maxarbvar)
  CALL mpupdatei(magazro,maxarbvar)
  CALL mpupdatei(magasty,maxarbvar)
  CALL mpupdatei(magacol1,maxarbvar)
  CALL mpupdatei(magacol2,maxarbvar)
  CALL mpupdatei(magaunits,maxarbvar)

  CALL mpupdatec(filename2d,maxarbvar*256)
  CALL mpupdatec(filename3d,maxarbvar*256)
  CALL mpupdatec(filenameu2d,maxarbvar*256)
  CALL mpupdatec(filenamev2d,maxarbvar*256)

!
!-----------------------------------------------------------------------
!
!  Input control parameters plotting boxes
!
!-----------------------------------------------------------------------
!
  IF(myproc == 0) THEN
    READ(unum,plot_boxes,ERR=100)
    WRITE(6,'(a)')'Namelist plot_box was successfully read.'
  END IF
  CALL mpupdatei(number_of_boxes,1)
  CALL mpupdatei(boxcol,1)
  CALL mpupdater(bctrx,10)
  CALL mpupdater(bctry,10)
  CALL mpupdater(blengx,10)
  CALL mpupdater(blengy,10)

  IF(number_of_boxes /= 0) THEN
    DO k=1,number_of_boxes
      WRITE(6,'(1x,a,i3,a,2f10.5)')                                     &
          'Center of box No.',k,' is at ',bctrx(k),bctry(k)
      WRITE(6,'(1x,a,i3,a,2f10.5)')                                     &
          'The size of box No.',k,' is ',blengx(k),blengy(k)
    END DO

    DO k=1,number_of_boxes
      bx1(k)=bctrx(k) - blengx(k)*0.5
      bx2(k)=bctrx(k) + blengx(k)*0.5
      by1(k)=bctry(k) - blengy(k)*0.5
      by2(k)=bctry(k) + blengy(k)*0.5
    END DO

  END IF
!
!
!-----------------------------------------------------------------------
!
!  Input control parameters plotting polylines
!
!-----------------------------------------------------------------------
!
  DO j=1,max_polys
    DO i=1,max_verts
      vertx(i,j) = -9999.
      verty(i,j) = -9999.
    END DO
  END DO

  IF(myproc == 0) THEN
    READ(unum,plot_polylines,ERR=100)
    WRITE(6,'(a)')'Namelist plot_polylines was successfully read.'

    IF(number_of_polys /= 0) THEN
      DO k=1,number_of_polys
        WRITE(6,'(1x,a,i2)')'The number of polyline is : ',k
        DO j = 1, max_verts
          IF(vertx(j,k) /= -9999. .AND. verty(j,k) /= -9999.)             &
            WRITE(6,'(1x,a,2f10.5)')                                    &
            'The position of vertices are: ',vertx(j,k),verty(j,k)
        END DO
      END DO
    END IF
  END IF  !myproc == 0
  CALL mpupdatei(number_of_polys,1)
  CALL mpupdatei(polycol,1)
  CALL mpupdater(vertx,max_verts*max_polys)
  CALL mpupdater(verty,max_verts*max_polys)

!
!-----------------------------------------------------------------------
!
!  Input control parameters for plotting trajectories
!  -added by Dan Dawson 12/03/04
!
!-----------------------------------------------------------------------
!

  IF(myproc == 0) THEN

    trajc_plt_opt = 0
    ntimes = 1
    trajc_fn_in = 'may20.trajc_001800-009000_006600'
    trajc_plt_bgn_time = 0.0
    trajc_plt_end_time = 0.0

    ntrajc_start = 1
    ntrajc_end = -1
    ntrajc_stride = 1

    traj_col = 1

    trajc_lbl_opt = 3
    trajc_lbl_frq = 10
    trajc_lbl_siz = 0.015
    trajc_lbl_fmt = 0

    trajc_mkr_typ = 6
    trajc_mkr_frq = 5
    trajc_mkr_siz = 0.004

    READ(unum,plot_trajectories,ERR=100)
    WRITE(6,'(a)')'Namelist plot_trajectories was successfully read.'

    IF (trajc_plt_opt > 0) THEN
      DO k = 1,ntimes
        WRITE(6,'(a,a)')'The trajectory file name is : ',trajc_fn_in(k)
      END DO
    END IF

  END IF  !myproc == 0
  CALL mpupdatei(trajc_plt_opt,1)
  CALL mpupdater(trajc_plt_bgn_time, 1)
  CALL mpupdater(trajc_plt_end_time, 1)
  CALL mpupdatei(traj_col          , 1)
  CALL mpupdatec(trajc_fn_in       , 256*nmax_times)
  CALL mpupdatei(ntimes            , 1)
  CALL mpupdatei(trajc_lbl_opt     , 1)
  CALL mpupdatei(trajc_lbl_frq     , 1)
  CALL mpupdater(trajc_lbl_siz     , 1)
  CALL mpupdatei(trajc_mkr_typ     , 1)
  CALL mpupdatei(trajc_mkr_frq     , 1)
  CALL mpupdater(trajc_mkr_siz     , 1)
  CALL mpupdatei(trajc_lbl_fmt     , 1)
  CALL mpupdatei(ntrajc_start      , 1)
  CALL mpupdatei(ntrajc_end        , 1)
  CALL mpupdatei(ntrajc_stride     , 1)
!
!-----------------------------------------------------------------------
!
!  Input control parameters for overlay one filed to many fields
!
!-----------------------------------------------------------------------
!
  IF(myproc == 0) THEN
    READ(unum,ovrlay_mul,ERR=100)
    WRITE(6,'(a)')'Namelist ovrlay_mul was successfully read.'
  END IF
  CALL mpupdatei(ovrlaymulopt,1)
  CALL mpupdatec(ovrname,12)
  CALL mpupdatei(ovrmul_num,1)
  CALL mpupdatec(ovrmulname,50*12)

  IF(ovrlaymulopt == 0) THEN
    ovrname(1:12)='            '
    ovrmul_num = 0
    DO i = 1,50
      ovrmulname(i)(1:12) ='            '
    END DO
  END IF
!
!-----------------------------------------------------------------------
!
!  Input control parameters for terrain overlay
!
!-----------------------------------------------------------------------
!
  IF(myproc == 0) THEN
    READ(unum,ovr_terrain,ERR=100)
    WRITE(6,'(a)')'Namelist ovr_terrain was successfully read.'
  END IF
  CALL mpupdatei(ovrtrn,1)
!
!-----------------------------------------------------------------------
!
!  Input control parameters for 3-D wireframe plotting
!
!-----------------------------------------------------------------------
!
  IF(myproc == 0) THEN
    READ(unum,wirfrm_plot,ERR=100)
    WRITE(6,'(a)')'Namelist wirfrm_plot was successfully read.'
  END IF
  CALL mpupdatei(w3dplt,1)
  CALL mpupdater(wisosf,1)
  CALL mpupdatei(q3dplt,1)
  CALL mpupdater(qisosf,1)
!
!-----------------------------------------------------------------------
!
!  Parameters for overlaying observations
!
!-----------------------------------------------------------------------
!
  ovrobs=0
  obscol=1
  obs_marktyp=4
  obs_marksz=0.002
  obs_valsz=0.01
  nsfcobfl=0
  sfcobfl(:)=' '

  IF(myproc == 0) THEN
    READ(unum,plot_obs,ERR=71)
    WRITE(6,'(a)')'Namelist plot_obs was successfully read.'

    nsfcobfl=min(mxsfcobfl,nsfcobfl)
    IF (ovrobs > 0 .AND. nsfcobfl > 0) THEN
      DO ifile=1,nsfcobfl
        WRITE(6,'(2x,a,i3,2a)') 'Surface obs file(',ifile, &
                '): ',TRIM(sfcobfl(ifile))
      END DO
    END IF
  END IF
  CALL mpupdatei(ovrobs,1)
  CALL mpupdatei(nsfcobfl,1)
  CALL mpupdatec(sfcobfl,mxsfcobfl*256)
  CALL mpupdatei(obscol,1)
  CALL mpupdatei(obs_marktyp,1)
  CALL mpupdater(obs_marksz,1)
  CALL mpupdater(obs_valsz,1)

  71 CONTINUE

!-----------------------------------------------------------------------
!
!  Parameters for overlaying airport location
!
!-----------------------------------------------------------------------
!
  ovrstam=0
  ovrstan=0
  ovrstav=0
  wrtstax=0

  IF(myproc == 0) THEN
    READ(unum,plot_sta, ERR=72)
    WRITE(6,'(a)') 'Namelist plot_sta was successfully read.'

    IF (ovrstaopt > 0) THEN
      WRITE(6,'(2x,a,a)') 'Station file name: ',TRIM(stalofl)
    END IF
  END IF
  CALL mpupdatei(ovrstaopt,1)
  CALL mpupdatei(ovrstam,1)
  CALL mpupdatei(ovrstan,1)
  CALL mpupdatei(ovrstav,1)
  CALL mpupdatei(wrtstax,1)
  CALL mpupdater(wrtstad,1)
  CALL mpupdatei(stacol,1)
  CALL mpupdatei(markprio,1)
  CALL mpupdatei(nsta_typ,1)
  CALL mpupdatei(sta_typ,30)
  CALL mpupdatei(sta_marktyp,30)
  CALL mpupdatei(sta_markcol,30)
  CALL mpupdater(sta_marksz,30)
  CALL mpupdatec(stalofl,256)

  72 CONTINUE

!-----------------------------------------------------------------------
!
!  Input control parameter for profile
!
!-----------------------------------------------------------------------

  xprof = 0.0
  yprof = 0.0

  IF(myproc == 0) THEN
    READ(unum,profile_cntl,ERR=100)
    WRITE(6,'(a)')'Namelist profile_cntl was successfully read.'

    IF (nprof > max_dim) THEN
      WRITE (6,'(1x,a,i4)') 'Too many profiles. Limited to ',max_dim
      nprof = max_dim
    END IF
  END IF
  CALL mpupdatei(profopt,1)
  CALL mpupdatei(nprof,1)
  CALL mpupdater(xprof,max_dim)
  CALL mpupdater(yprof,max_dim)
  CALL mpupdatei(npicprof,1)
  CALL mpupdatei(uprof,1)
  CALL mpupdater(uprmin,1)
  CALL mpupdater(uprmax,1)
  CALL mpupdatei(vprof,1)
  CALL mpupdater(vprmin,1)
  CALL mpupdater(vprmax,1)
  CALL mpupdatei(wprof,1)
  CALL mpupdater(wprmin,1)
  CALL mpupdater(wprmax,1)
  CALL mpupdatei(ptprof,1)
  CALL mpupdater(ptprmin,1)
  CALL mpupdater(ptprmax,1)
  CALL mpupdatei(pprof,1)
  CALL mpupdater(pprmin,1)
  CALL mpupdater(pprmax,1)
  CALL mpupdatei(qvprof,1)
  CALL mpupdater(qvprmin,1)
  CALL mpupdater(qvprmax,1)
  CALL mpupdatei(qcprof,1)
  CALL mpupdater(qcpmin,1)
  CALL mpupdater(qcpmax,1)
  CALL mpupdatei(qrprof,1)
  CALL mpupdater(qrpmin,1)
  CALL mpupdater(qrpmax,1)
  CALL mpupdatei(qiprof,1)
  CALL mpupdater(qipmin,1)
  CALL mpupdater(qipmax,1)
  CALL mpupdatei(qsprof,1)
  CALL mpupdater(qspmin,1)
  CALL mpupdater(qspmax,1)
  CALL mpupdatei(qhprof,1)
  CALL mpupdater(qhpmin,1)
  CALL mpupdater(qhpmax,1)
  CALL mpupdatei(rhprof,1)
  CALL mpupdater(rhpmin,1)
  CALL mpupdater(rhpmax,1)
  CALL mpupdatei(kmhprof,1)
  CALL mpupdater(kmhpmin,1)
  CALL mpupdater(kmhpmax,1)
  CALL mpupdatei(kmvprof,1)
  CALL mpupdater(kmvpmin,1)
  CALL mpupdater(kmvpmax,1)
  CALL mpupdatei(tkeprof,1)
  CALL mpupdater(tkepmin,1)
  CALL mpupdater(tkepmax,1)
  CALL mpupdatei(rfprof,1)
  CALL mpupdater(rfpmin,1)
  CALL mpupdater(rfpmax,1)
  CALL mpupdatei(pteprf,1)
  CALL mpupdater(ptepmin,1)
  CALL mpupdater(ptepmax,1)
  CALL mpupdatei(upprof,1)
  CALL mpupdater(uppmin,1)
  CALL mpupdater(uppmax,1)
  CALL mpupdatei(vpprof,1)
  CALL mpupdater(vppmin,1)
  CALL mpupdater(vppmax,1)
  CALL mpupdatei(wpprof,1)
  CALL mpupdater(wppmin,1)
  CALL mpupdater(wppmax,1)
  CALL mpupdatei(ptpprf,1)
  CALL mpupdater(ptppmin,1)
  CALL mpupdater(ptppmax,1)
  CALL mpupdatei(ppprof,1)
  CALL mpupdater(pppmin,1)
  CALL mpupdater(pppmax,1)
  CALL mpupdatei(qvpprf,1)
  CALL mpupdater(qvppmin,1)
  CALL mpupdater(qvppmax,1)
  CALL mpupdatei(vorpprf,1)
  CALL mpupdater(vorppmin,1)
  CALL mpupdater(vorppmax,1)
  CALL mpupdatei(divpprf,1)
  CALL mpupdater(divppmin,1)
  CALL mpupdater(divppmax,1)

  CALL mpupdater(zprofbgn,1)
  CALL mpupdater(zprofend,1)

  CALL mpupdatei(tsoilprof,1)
  CALL mpupdater(tsoilprofmin,1)
  CALL mpupdater(tsoilprofmax,1)
  CALL mpupdatei(qsoilprof,1)
  CALL mpupdater(qsoilprofmin,1)
  CALL mpupdater(qsoilprofmax,1)
  CALL mpupdater(zsoilprofbgn,1)
  CALL mpupdater(zsoilprofend,1)

  CALL mpupdatei(nxprpic,1)
  CALL mpupdatei(nyprpic,1)

  dirname = './'
  outfilename = ' '
  lvldbg      = 0
  IF(myproc == 0) THEN
    READ(unum,output,ERR=100)
    WRITE(6,'(a)')'Namelist output was successfully read.'

    lenstr = LEN_TRIM(dirname)
    IF(lenstr > 0) THEN
      IF(dirname(lenstr:lenstr) /= '/') THEN
        dirname(lenstr+1:lenstr+1) = '/'
        lenstr = lenstr + 1
      END IF
    ELSE
      dirname = './'
    END IF

  END IF
  CALL mpupdatec(dirname,256)
  CALL mpupdatec(outfilename,256)
  CALL mpupdatei(iwtype,1)
  CALL mpupdatei(lvldbg,1)

  IF (myproc == 0) WRITE(6,'(/,a,/)') 'NAMELIST file was read and parsed succefully.'

  IF (unum /= 5 .AND. myproc == 0) THEN
    CLOSE( unum )
    CALL retunit( unum )
  END IF

  GO TO 10

  100 WRITE(6,'(a)')'Error reading NAMELIST file. Job stopped in ARPSPLT.'
  CALL mpexit(1)
  STOP

  10   CONTINUE

  RETURN

END SUBROUTINE initpltpara
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CTR3D                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE ctr3d(b,x,y,z, x1,x2,dx,y1,y2,dy,z1,z2,dz,                   &
           nx,ibgn,iend, ny,jbgn,jend, nz,kbgn,kend,                    &
           label,time,slicopt, kslice, jslice, islice,                  &
           n,xp,yp,axy2d,av2d,zp, runname, factor,tem1,tem2,tem3,       &
           tem4,bb,tem5,hterain,pltopt)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    Set-up 2-d slices of a 3-d data array to contour with
!      subroutine ctr2d.
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!
!  MODIFICATION HISTORY:
!    6/08/92  Added full documentation (K. Brewster)
!
!  12/25/1992 M. Xue and H. Jin
!    Added capability to plot arbitary cross sections.
!
!  8/28/1994 M. Zou
!    Added color shader to contour plot,add full documentation
!
!  3/25/96 (K. Brewster)
!    Added variables isize,jsize,ksize
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    b        3-dimension array of variable
!    x        x-coord of scalar point (km)
!    y        y-coord of scalar point (km)
!    z        z-coord of scalar point in computation space (km)
!    label    character string describing the contents of a plot
!    time     model runing time step
!    slicopt  slice orientation indicator
!             = 1, x-y slice of at k=kslice is plotted.
!             = 2, x-z slice of at j=jslice is plotted.
!             = 3, y-z slice of at i=islice is plotted.
!             = 4, horizontal slice at z index islice is plotted.
!             = 5, xy-z cross section of wind islice is plotted.
!             = 6, data field on constant p-level is plotted.
!             = 0, all of the three slices above are plotted.
!    axy2d    2d x-y array
!    av2d     2D array for the vertical slice
!    xp       x-coordinate of grid points on arbitary vertical
!               cross-section
!    yp       y-coordinate of grid points on arbitary vertical
!               cross-section
!    zp       z-coordinate of grid points on arbitary vertical
!               cross-section
!    runname  character string describing the model run
!    factor   scaling factor
!    hterain  2-D terrain data for contour
!    trnplt   flag to plot terrain (0/1)
!  WORK ARRAY:
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
!    tem4     Temporary work array.
!    tem5     Temporary work array.
!
!  (These arrays are defined and used locally (i.e. inside this
!   subroutine), they may also be passed into routines called by
!   this one. Exiting the call to this subroutine, these temporary
!   work arrays may be used for other purposes therefore their
!   contents overwritten. Please examine the usage of work arrays
!   before you alter the code.)
!
!   pp01      The pressure (mb) value at the specific p-level
!   ercpl     reciprocal of exponent
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz
  INTEGER :: n

  REAL :: b(nx,ny,nz)
  REAL :: x(nx,ny,nz)
  REAL :: y(nx,ny,nz)
  REAL :: z(nx,ny,nz)

  REAL :: axy2d(nx,ny)
  REAL :: av2d(n,nz),zp(n,nz)
  REAL :: xp(n),yp(n)

  REAL :: x1,x2,dx,y1,y2,dy,z1,z2,dz
  INTEGER :: ibgn,iend,jbgn,jend,kbgn,kend,length

  CHARACTER (LEN=9) :: timhms
  CHARACTER (LEN=*) :: label

  REAL :: time

  INTEGER :: slicopt
  INTEGER :: kslice,jslice,islice

  CHARACTER (LEN=*) :: runname

  REAL :: factor

  INTEGER :: trnplt            ! plot terrain option (0/1/2/3)
  INTEGER :: pltopt            ! plot variable option (0/1/2/3)
  REAL :: hterain(nx,ny)       ! The height of the terrain.

  REAL :: tem1(*)
  REAL :: tem2(*)
  REAL :: tem3(*)
  REAL :: tem4(*)
  REAL :: bb(nx,ny,nz)
  REAL :: tem5(*)          ! size must >= 6*nx*ny
!
!-----------------------------------------------------------------------
!
!  Some constants
!
!-----------------------------------------------------------------------
!
  REAL            :: pp01
  REAL, PARAMETER :: ercpl = 0.3678794              ! exp(-1.0)
!
!-----------------------------------------------------------------------
!
!  Common blocks for plotting control parameters
!
!-----------------------------------------------------------------------
!
  REAL :: x01,y01                  ! the first  point of interpolation
  REAL :: x02,y02                  ! the second point of interpolation
  REAL :: zlevel                   ! the given height of the slice
  REAL :: sinaf,cosaf,dist,sqrtdxy
  COMMON /slicev/x01,y01,x02,y02,sinaf,cosaf,dist,sqrtdxy
  COMMON /sliceh/zlevel

  INTEGER :: ovrtrn               ! overlay terrain option (0/1)
  REAL :: trninc,trnmin, trnmax   ! terrain interval minimum, maximum
  REAL :: ztmin,ztmax
  COMMON /trnpar/ trnplt,ovrtrn,trninc,trnmin, trnmax,ztmin,ztmax

  INTEGER :: smooth
  COMMON /smoothopt/smooth

  INTEGER :: xfont   ! the font of character
  INTEGER :: haxisu, vaxisu
  INTEGER :: lbaxis
  INTEGER :: tickopt
  INTEGER :: axlbfmt
  REAL :: hmintick,vmajtick,vmintick,hmajtick
  COMMON /var_par/ xfont,haxisu,vaxisu,lbaxis,tickopt,hmintick,         &
          vmajtick, vmintick,hmajtick,axlbfmt
  CHARACTER (LEN=4) :: stem2
  CHARACTER (LEN=1) :: stem1
  REAL :: x_tmp
  COMMON /tmphc2/ x_tmp

  REAL :: tmpx, tmpy
  CHARACTER (LEN=20) :: distc
  REAL :: x101, y101, x102,y102
  COMMON /slicev1/x101, y101, x102,y102
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,ij,ik,jk,isize,jsize,ksize, llabel
  CHARACTER (LEN=120) :: label_copy
  CHARACTER (LEN=120) :: title

  INTEGER :: wrtflag
  CHARACTER (LEN=90) :: levlab
  CHARACTER (LEN=50) :: timelab
  CHARACTER (LEN=25) :: timestring
  COMMON /timelev/wrtflag, timelab, levlab, timestring

  INTEGER :: xpbgn,xpend,ypbgn,ypend
  COMMON /processors/ xpbgn, xpend, ypbgn, ypend

  INTEGER :: agl
  COMMON /agl_or_asl/ agl

  INTEGER :: idsize, jdsize, mnsize
  INTEGER :: tinds, tind1,tind2,tind3,tind4,tind5,tind6
                ! temporary arrays index, assume size of tem5 > 6*nx*ny
!-----------------------------------------------------------------------
!
!  Common blocks for including iskip, jskip
!
!-----------------------------------------------------------------------
!
  REAL    :: xw1,xw2,yw1,yw2
  INTEGER :: iskip, jskip
  COMMON /pltwdw/ xw1,xw2,yw1,yw2, iskip, jskip

  REAL, ALLOCATABLE :: var_tr(:,:)
  REAL, ALLOCATABLE :: xptr(:,:)
  REAL, ALLOCATABLE :: zptr(:,:)
  REAL, ALLOCATABLE :: temptr1(:,:),temptr2(:,:),temptr3(:,:)
  REAL, ALLOCATABLE :: temptr4(:,:),temptr5(:,:),temptr6(:,:)
  REAL, ALLOCATABLE :: xptr1d(:),yptr1d(:)
  REAL :: xptr1,xptr2,dxptr,xtr,ytr
  INTEGER :: itrajc_start, itrajc_end, iptr
  INTEGER :: ireturn,j_bgn,j_end
  INTEGER :: j_bgn_last,j_end_last,istatus,kz

  REAL :: yxratio
  COMMON /yratio/ yxratio       ! the scaling factor the y/x ratio.

!----------------------------------------------------------------------
!
! Include files
!
!---------------------------------------------------------------------
  INCLUDE 'mp.inc'
  INCLUDE 'arpstrajc.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  isize = (iend-ibgn)+1
  jsize = (jend-jbgn)+1
  ksize = (kend-kbgn)+1

  idsize = isize            ! global maximum isize
  jdsize = jsize
  CALL mpmaxi(idsize)
  CALL mpmaxi(jdsize)

  mnsize = idsize*jdsize
  mnsize = MAX(mnsize,idsize*ksize,jdsize*ksize)

  tind1 = 1             ! reuse a 3d temporary array 'tem5' as several 2D
  tind2 = tind1+mnsize  ! arrays inside ctr2d
  tind3 = tind2+mnsize
  tind4 = tind3+mnsize
  tind5 = tind4+mnsize
  tind6 = tind5+mnsize

!  tinds = SIZE(tem5)
!  IF (tinds < 6*mnsize) THEN
!    WRITE(6,'(3a)') 'ERROR: temporary array tem5 is too small ',        &
!                    'inside ctr3d while plotting ', label
!    CALL arpsstop('Temporary array too small inside ctr3d.',1)
!  END IF

  label_copy = label
  llabel = 120
  CALL xstrlnth(label_copy, llabel)

  IF(myproc == 0) CALL xpscmnt('Start plotting '//label_copy(1:llabel))
!
!-----------------------------------------------------------------------
!
!  Set up terrain, if needed.
!
!-----------------------------------------------------------------------
!
  IF(trnplt == 1 .OR.trnplt == 2 .OR. ovrtrn == 1)  THEN
    DO j=jbgn,jend
      DO i=ibgn,iend
        ij = i-ibgn+1 + (j-jbgn)*isize
        tem4(ij)=hterain(i,j)
      END DO
    END DO
  END IF

  CALL get_forecast_hms( time, timhms )

  WRITE(timelab,'(a,F8.1,a)') 'T=',time,' s ('//TRIM(timhms)//')'
  CALL get_time_string ( time, timestring,'Z  ',0 )

  IF ( slicopt == 2  .OR. slicopt == 3  .OR. slicopt == 5 .OR.          &
       slicopt == 10 .OR. slicopt == 11 .OR. slicopt == 12 ) THEN
    CALL cal_dist(haxisu,dx,dy,x01,y01,x02,y02,                         &
                  slicopt,tmpx,tmpy,distc)
  END IF

  IF(slicopt == 1 .OR. slicopt == 0 ) THEN

    k = kslice
    DO j=jbgn,jend
      DO i=ibgn,iend
        ij = i-ibgn+1 + (j-jbgn)*isize
        tem1(ij) = -9999.0
        IF(b(i,j,k) /= -9999.0) tem1(ij)=b(i,j,k)*factor
        tem2(ij)=x(i,j,k)
        tem3(ij)=y(i,j,k)
      END DO
    END DO


    IF (k /= 2) THEN
      WRITE(levlab,'(''GRID LEVEL='',I3)')k
      WRITE(title,'(a)') label
    ELSE
      WRITE(levlab,'(''First level above ground (surface)'')')
      WRITE(title,'(a)') label
    END IF

    length = LEN_TRIM(title)
    CALL strmin ( title, length)

    DO i=1,smooth
      CALL smooth9pmv(tem1,isize,jsize,1,isize,1,jsize,tem5)
    END DO

    CALL ctr2d(tem1,tem2,tem3, x1,x2,dx, y1,y2,dy,iskip+1,jskip+1,      &
               isize,jsize,title(1:length),runname,                     &
               tem4,slicopt,pltopt,mnsize,                              &
               tem5(tind1),tem5(tind2),tem5(tind3),                     &
               tem5(tind4),tem5(tind5),tem5(tind6))

!
!-----------------------------------------------------------------------
!
!  slicopt=2   Plot x-z cross-section
!
!-----------------------------------------------------------------------
!
  ELSE IF (slicopt == 2 .OR. slicopt == 0 ) THEN

    x_tmp = y(1,jslice,1)

    j = jslice
    DO k=kbgn,kend
      DO i=ibgn,iend
        ik = i-ibgn+1 + (k-kbgn)*isize
        tem1(ik) = -9999.0
        IF(b(i,j,k) /= -9999.0) tem1(ik)=b(i,j,k)*factor
        tem2(ik)=x(i,j,k)
        tem3(ik)=z(i,j,k)
      END DO
    END DO

    j = j + (ny-3)*(ypbgn-1)
    dist = (j-1.5)*tmpy

    length= LEN(distc)
    CALL strmin ( distc, length)
    WRITE(levlab,'(''X-Z Plane at Y='',F8.1,A)')dist,distc(1:length)

    WRITE(title,'( a )') label

    length = LEN_TRIM(title)
    CALL strmin ( title, length)

    DO i=1,smooth
      CALL smooth9pmv(tem1,isize,ksize,1,isize,1,ksize,tem5)
    END DO

    CALL ctr2d(tem1,tem2,tem3, x1,x2,dx, z1,z2,dz,iskip+1,1,            &
               isize,ksize,title(1:length),runname,                     &
               tem4,slicopt,pltopt,mnsize,                              &
               tem5(tind1),tem5(tind2),tem5(tind3),                     &
               tem5(tind4),tem5(tind5),tem5(tind6))

!
!-----------------------------------------------------------------------
!
!  slicopt=3   Plot y-z cross-section
!
!-----------------------------------------------------------------------
!
  ELSE IF ( slicopt == 3 .OR. slicopt == 0) THEN

    x_tmp = x(islice,1,1)

    i = islice
    DO k=kbgn,kend
      DO j=jbgn,jend
        jk = j-jbgn+1 + (k-kbgn)*jsize
        tem1(jk) = -9999.0
        IF(b(i,j,k) /= -9999.0) tem1(jk)=b(i,j,k)*factor
        tem2(jk)=y(i,j,k)
        tem3(jk)=z(i,j,k)
      END DO
    END DO

    i = i + (nx-3)*(xpbgn-1)
    dist = (i-1.5)*tmpx
    length= LEN_TRIM(distc)
    CALL strmin ( distc, length)
    WRITE(levlab,'(''Y-Z Plane at X='',F8.1,A)')dist,distc(1:length)

    WRITE(title,'( a )' ) label

    length = LEN_TRIM(title)
    CALL strmin ( title, length)

    DO i=1,smooth
      CALL smooth9pmv(tem1,jsize,ksize,1,jsize,1,ksize,tem5)
    END DO

    CALL ctr2d(tem1,tem2,tem3, y1,y2,dy, z1,z2,dz,jskip+1,1,            &
               jsize,ksize,title(1:length),runname,                     &
               tem4,slicopt,pltopt,mnsize,                              &
               tem5(tind1),tem5(tind2),tem5(tind3),                     &
               tem5(tind4),tem5(tind5),tem5(tind6))
!
!-----------------------------------------------------------------------
!
!  slicopt=4   Plot horizontal slice at given height
!  slicopt=6   Plot constant pressure slice at given pressure(mb)
!  slicopt=7   Plot isentropic surfaces
!
!-----------------------------------------------------------------------
!
  ELSE IF( slicopt == 4 .OR. slicopt == 6 .OR. slicopt == 7 ) THEN

    DO k=kbgn,kend
      DO j=jbgn,jend
        DO i=ibgn,iend
          bb(i,j,k) = -9999.0
          IF(b(i,j,k) /= -9999.0) bb(i,j,k)= b(i,j,k)*factor
        END DO
      END DO
    END DO

    IF( agl == 1 ) THEN
!     print*,'AGL: calling hintrp2'
      CALL hintrp2(nx,ny,nz,kbgn,kend,bb,z,zlevel,hterain,axy2d)
    ELSE
      CALL hintrp1(nx,ny,nz,kbgn,kend,bb,z,zlevel,axy2d)
    END IF

    DO j=jbgn,jend
      DO i=ibgn,iend
        ij = i-ibgn+1 + (j-jbgn)*isize
        tem1(ij)=axy2d(i,j)
        tem2(ij)=x(i,j,2)
        tem3(ij)=y(i,j,2)
      END DO
    END DO

    IF(slicopt == 4) THEN
      IF ( agl == 1 ) THEN
        IF( zlevel > 1.0) THEN
          WRITE(levlab,'(''Z='',F7.3,'' km AGL'')') zlevel
        ELSE
          WRITE(levlab,'(''Z='',F7.0,'' m AGL'')') (1000.*zlevel)
        END IF
      ELSE
        IF( zlevel > 1.0) THEN
          WRITE(levlab,'(''Z='',F7.3,'' km MSL'')') zlevel
        ELSE
          WRITE(levlab,'(''Z='',F7.0,'' m MSL'')') (1000.*zlevel)
        END IF
      END IF
    ELSE IF(slicopt == 6) THEN
      pp01 = 0.01*ercpl**zlevel
      WRITE(levlab,'(''P='',F7.2,A)') pp01, ' MB'
    ELSE
      WRITE(levlab,'(''Theta='',F5.1,A)') zlevel, ' (K)'
    END IF

    WRITE(title,'(a)') label

    length = LEN_TRIM(title)
    CALL strmin ( title, length)

    DO i=1,smooth
      CALL smooth9pmv(tem1,isize,jsize,1,isize,1,jsize,tem5)
    END DO

    CALL ctr2d(tem1,tem2,tem3, x1,x2,dx, y1,y2,dy,iskip+1,jskip+1,      &
               isize,jsize,title(1:length),runname,                     &
               tem4,slicopt,pltopt,mnsize,                              &
               tem5(tind1),tem5(tind2),tem5(tind3),                     &
               tem5(tind4),tem5(tind5),tem5(tind6))

!-----------------------------------------------------------------------
!  slicopt= 12   Plot vectical slice following trajectories
!-----------------------------------------------------------------------

  ELSE IF( slicopt == 12 .AND. trajc_plt_opt == 2 ) THEN

    k = 1            ! supporting only ntimes=1 for slice plotting
    i = itrajc_index ! i=itrajc_start, itrajc_end, ntrajc_stride

    j_bgn = npoints_bgn(i,k)
    j_end = npoints_end(i,k)

    ALLOCATE(xptr(j_bgn:j_end,kbgn:kend),STAT=istatus)
    ALLOCATE(xptr1d(j_bgn:j_end),STAT=istatus)
    ALLOCATE(yptr1d(j_bgn:j_end),STAT=istatus)

    ALLOCATE(temptr1(j_bgn:j_end,kbgn:kend),STAT=istatus)
    ALLOCATE(temptr2(j_bgn:j_end,kbgn:kend),STAT=istatus)
    ALLOCATE(temptr3(j_bgn:j_end,kbgn:kend),STAT=istatus)
    ALLOCATE(temptr4(j_bgn:j_end,kbgn:kend),STAT=istatus)
    ALLOCATE(temptr5(j_bgn:j_end,kbgn:kend),STAT=istatus)
    ALLOCATE(temptr6(j_bgn:j_end,kbgn:kend),STAT=istatus)

    ALLOCATE(zptr(j_bgn:j_end,nz),STAT=istatus)
    ALLOCATE(var_tr(j_bgn:j_end,1:nz),STAT=istatus)

    xptr1d(j_bgn:j_end)=xtrajc(j_bgn:j_end,i,k)*0.001
    yptr1d(j_bgn:j_end)=ytrajc(j_bgn:j_end,i,k)*0.001

    PRINT*,'plotting vertical slice along trajectory No. ', i

    CALL sectvrt(nx,ny,nz,b,x,y,z,dx,dy,var_tr,zptr,j_end-j_bgn+1,xptr1d,yptr1d)

    DO kz=kbgn,kend
      xptr(j_bgn,kz)=0.0
      DO j=j_bgn+1,j_end
        xptr(j,kz) = xptr(j-1,kz)+sqrt((xptr1d(j)-xptr1d(j-1))**2+(yptr1d(j)-yptr1d(j-1))**2)
      END DO
    END DO

    xptr1 = xptr(j_bgn,kbgn)
    xptr2 = xptr(j_end,kbgn)
    dxptr=(xptr2-xptr1)/(j_end-j_bgn)

    length=LEN_TRIM(distc)
    CALL strmin ( distc, length)

    IF(axlbfmt == -1 .OR. axlbfmt == 1 ) THEN
      WRITE(levlab,                                                     &
      '(''Vert. plane along traj. '',4(A,F8.1),A,A)')                   &
      '(',xptr1d(j_bgn),',',yptr1d(j_bgn),') to (',xptr1d(j_end)        &
      ,',',yptr1d(j_end),') ', distc(1:length)
      WRITE(title,'(a)') label
    ELSE IF(axlbfmt == 0) THEN
      WRITE(levlab,                                                     &
      '(''Vert. plane along traj. '',4(A,I5),A,A)')                     &
      '(',nint(xptr1d(j_bgn)),',',nint(yptr1d(j_bgn)),') through (',    &
      nint(xptr1d(j_end)),',',nint(yptr1d(j_end)),') ', distc(1:length)
      WRITE(title,'(a)') label
    ELSE
      WRITE(title,'(''V-W '',A)') label
      WRITE(levlab,                                                     &
      '(''Vert. plane along traj. '',4(A,F8.2),A,A)')                   &
      '(',xptr1d(j_bgn),',',yptr1d(j_bgn),') through (',xptr1d(j_end)   &
      ,',',yptr1d(j_end),') ', distc(1:length)
    END IF

    length = LEN_TRIM(title)
    CALL strmin ( title, length)

    DO iptr=1,smooth
      CALL smooth9pmv(var_tr(j_bgn:j_end,kbgn:kend),j_end-j_bgn+1,      &
           kend-kbgn+1,1,j_end-j_bgn+1,1,kend-kbgn+1, temptr1)
    END DO

    yxratio = (xptr2-xptr1)/(z2-z1)

    CALL ctr2d(var_tr(j_bgn:j_end,kbgn:kend),xptr(j_bgn:j_end,kbgn:kend), &
         zptr(j_bgn:j_end,kbgn:kend), xptr1,xptr2,dxptr, z1,z2,dz,1,1,    &
         j_end-j_bgn+1, kend-kbgn+1, title(1:length),runname,             &
         tem4,slicopt,pltopt,(j_end-j_bgn+1)*(kend-kbgn+1),               &
         temptr1,temptr2,temptr3,temptr4,temptr5,temptr6)

    DEALLOCATE(xptr)
    DEALLOCATE(xptr1d)
    DEALLOCATE(yptr1d)

    DEALLOCATE(temptr1)
    DEALLOCATE(temptr2)
    DEALLOCATE(temptr3)
    DEALLOCATE(temptr4)
    DEALLOCATE(temptr5)
    DEALLOCATE(temptr6)

    DEALLOCATE(zptr)
    DEALLOCATE(var_tr)

!-----------------------------------------------------------------------
!  slicopt= 5   Plot vectical slice through two given points
!-----------------------------------------------------------------------

  ELSE IF( slicopt == 5 ) THEN

    CALL sectvrt(nx,ny,nz,b,x,y,z,dx,dy,av2d,zp,n,xp,yp)

    DO k=kbgn,kend
      DO i=ibgn,iend
        ik = i-ibgn+1 + (k-kbgn)*isize
        tem1(ik) = -9999.0
        IF(av2d(i,k) /= -9999.0) tem1(ik)=av2d(i,k)*factor
        tem2(ik)=x1+(i-ibgn)* sqrtdxy
        tem3(ik)=zp(i,k)
      END DO
    END DO

    IF(axlbfmt == -1 .OR. axlbfmt == 1 ) THEN
      length=LEN_TRIM(distc)
      CALL strmin ( distc, length)
      WRITE(levlab,                                                     &
          '(''Vertical Plane From '',4(A,F8.1),A,A)')                   &
          '(',x101,',',y101,') through (',x102,',',y102,') ',           &
          distc(1:length)
      WRITE(title,'(a)') label
    ELSE IF(axlbfmt == 0) THEN
      length= LEN_TRIM(distc)
      CALL strmin ( distc, length)
      WRITE(levlab,                                                     &
          '(''Vertical Plane From '',4(A,I5),A,A)')                     &
          '(',NINT(x101),',',NINT(y101),') through (',NINT(x102),','    &
          ,NINT(y102),') ', distc(1:length)
      WRITE(title,'(a)') label
    ELSE
      length=LEN_TRIM(distc)
      CALL strmin ( distc, length)
!     WRITE(stem1,'(i1)')axlbfmt
!     WRITE(stem2,'(a3,a1)')'f8.',stem1

      WRITE(title,'(''V-W '',A)') label
      WRITE(levlab,                                                     &
          '(''Vertical Plane From '',4(A,F8.2),A,A)')                   &
          '(',x101,',',y101,') through (',x102,',',y102,') ',           &
          distc(1:length)
    END IF

    length = LEN_TRIM(title)
    CALL strmin ( title, length)

    DO i=1,smooth
      CALL smooth9pmv(tem1,isize,ksize,1,isize,1,ksize,tem5)
    END DO

    CALL ctr2d(tem1,tem2,tem3, x1,x2,sqrtdxy, z1,z2,dz,1,1,             &
               isize,ksize,title(1:length),runname,                     &
               tem4,slicopt,pltopt,mnsize,                              &
               tem5(tind1),tem5(tind2),tem5(tind3),                     &
               tem5(tind4),tem5(tind5),tem5(tind6))
!
!-----------------------------------------------------------------------
!
!  slicopt=9   Plot x-y cross-section of the soil model
!
!  06/03/2002 Zuwen He
!
!  slicopt (9) is the same as slicopt (1), except that
!  the labels.
!
!-----------------------------------------------------------------------
!
  ELSE IF(slicopt == 9) THEN

    k = kslice
    DO j=jbgn,jend
      DO i=ibgn,iend
        ij = i-ibgn+1 + (j-jbgn)*isize
        tem1(ij) = -9999.0
        IF(b(i,j,k) /= -9999.0) tem1(ij)=b(i,j,k)*factor
        tem2(ij)=x(i,j,k)
        tem3(ij)=y(i,j,k)
      END DO
    END DO

    WRITE(levlab,'(''GRID LEVEL (SOIL) ='',I3)')k
    WRITE(title,'(a)') label

    length = LEN_TRIM(title)
    CALL strmin ( title, length)

    DO i=1,smooth
      CALL smooth9pmv(tem1,isize,jsize,1,isize,1,jsize,tem5)
    END DO

    CALL ctr2d(tem1,tem2,tem3, x1,x2,dx, y1,y2,dy,iskip+1,jskip+1,      &
               isize,jsize,title(1:length),runname,                     &
               tem4,slicopt,pltopt,mnsize,                              &
               tem5(tind1),tem5(tind2),tem5(tind3),                     &
               tem5(tind4),tem5(tind5),tem5(tind6))
!
!-----------------------------------------------------------------------
!
! Zuwen He, 06/06/2002
!
!  slicopt=10  Plot x-z cross-section of the soil model.
!
!-----------------------------------------------------------------------
!
  ELSE IF (slicopt == 10) THEN

    x_tmp = y(1,jslice,1)

    j = jslice
    DO k=kbgn,kend
      DO i=ibgn,iend
        ik = i-ibgn+1 + (k-kbgn)*isize
        tem1(ik) = -9999.0
        IF(b(i,j,k) /= -9999.0) tem1(ik)=b(i,j,k)*factor
        tem2(ik)=x(i,j,k)
        tem3(ik)=z(i,j,k)
      END DO
    END DO

    j = j + (ny-3)*(ypbgn-1)
    dist = (j-1.5)*tmpy
    length=LEN_TRIM(distc)
    CALL strmin ( distc, length)
    WRITE(levlab,'(''X-Z Plane (Soil) at Y='',F8.1,A)')dist,distc(1:length)

    WRITE(title,'( a )') label

    length = LEN_TRIM(title)
    CALL strmin ( title, length)

    DO i=1,smooth
      CALL smooth9pmv(tem1,isize,ksize,1,isize,1,ksize,tem5)
    END DO

    CALL ctr2d(tem1,tem2,tem3, x1,x2,dx, z1,z2,dz,iskip+1,1,            &
               isize,ksize,title(1:length),runname,                     &
               tem4,slicopt,pltopt,mnsize,                              &
               tem5(tind1),tem5(tind2),tem5(tind3),                     &
               tem5(tind4),tem5(tind5),tem5(tind6))
!
!-----------------------------------------------------------------------
!
!  slicopt=11   Plot y-z cross-section of the soil model
!
!-----------------------------------------------------------------------
!
  ELSE IF ( slicopt == 11) THEN

    x_tmp = x(islice,1,1)

    i = islice
    DO k=kbgn,kend
      DO j=jbgn,jend
        jk = j-jbgn+1 + (k-kbgn)*jsize
        tem1(jk) = -9999.0
        IF(b(i,j,k) /= -9999.0) tem1(jk)=b(i,j,k)*factor
        tem2(jk)=y(i,j,k)
        tem3(jk)=z(i,j,k)
      END DO
    END DO

    i = i + (nx-3)*(xpbgn-1)
    dist = (i-1.5)*tmpx
    length=LEN_TRIM(distc)
    CALL strmin ( distc, length)
    WRITE(levlab,'(''Y-Z Plane (Soil) at X='',F8.1,A)')dist,distc(1:length)

    write (*,*) "levlab", levlab

    WRITE(title,'( a )' ) label

    length = LEN_TRIM(title)
    CALL strmin ( title, length)

    DO i=1,smooth
      CALL smooth9pmv(tem1,jsize,ksize,1,jsize,1,ksize,tem5)
    END DO

    CALL ctr2d(tem1,tem2,tem3, y1,y2,dy, z1,z2,dz,jskip+1,1,            &
               jsize,ksize,title(1:length),runname,                     &
               tem4,slicopt,pltopt,mnsize,                              &
               tem5(tind1),tem5(tind2),tem5(tind3),                     &
               tem5(tind4),tem5(tind5),tem5(tind6))

  END IF

  IF(myproc == 0) CALL xpscmnt('End plotting '//label_copy(1:llabel))

  RETURN
END SUBROUTINE ctr3d
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CTR2D                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE ctr2d(a,x,y,xl,xr,dx,yb,yt,dy,istep,jstep,m,n,title,runname, &
                 hterain,slicopt,pltopt,mnsize,                         &
                 plota,plotx,ploty,iwrk,xwk,ywk)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    Generate contour plots of 2-d field A given its coordinates
!      using ZXPLOT package..
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!
!  MODIFICATION HISTORY:
!
!  6/08/92 (K. Brewster)
!  Added full documentation.
!
!  8/28/94 (M. Zou)
!  Added color routing , overlay terrain.
!
!  1/24/96 (J. Zong and M. Xue)
!  Fixed a problem related to finding the minimum and maximum of the
!  2D array, a, when there exist missing data. Initial min. and max.
!  should be set to values other than the missing value, -9999.0.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    a        2-dimensional slice of data to contour
!
!    x        x coordinate of grid points in plot space (over on page)
!    y        y coordinate of grid points in plot space (up on page)
!
!
!    xl       Left bound of the physical domain
!    xr       Right bound of the physical domain
!    dx    Spacing between x-axis tick marks
!    yb       Bottom bound of the physical domain.
!    yt       Top bound of the physical domain.
!    dy    Spacing between y-axis tick marks
!
!    m        first dimension of a
!    n        second dimension of a
!
!    title    character string describing the contents of a
!    runname  character string describing the model run
!
!    hterain  2-D terrain data to contour
!    slicopt  slice orientation indicator
!             = 1, x-y slice of at k=kslice is plotted.
!             = 2, x-z slice of at j=jslice is plotted.
!             = 3, y-z slice of at i=islice is plotted.
!             = 4, horizontal slice at z index islice is plotted.
!             = 5, xy-z cross section of wind islice is plotted.
!             = 6, data field on constant p-level is plotted.
!             = 0, all of the three slices above are plotted.
!    plot      variable plot option (0/1/2/3)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INCLUDE 'arpsplt.inc'

  INTEGER, INTENT(IN) :: m,n

  REAL,    INTENT(IN) :: a(m,n)
  REAL,    INTENT(IN) :: x(m,n)
  REAL,    INTENT(IN) :: y(m,n)
  REAL,    INTENT(IN) :: xl,xr,dx,yb,yt,dy
  INTEGER, INTENT(IN) :: istep, jstep
  REAL,    INTENT(IN) :: hterain(m,n)

  CHARACTER(LEN=*), INTENT(IN) :: runname
  CHARACTER(LEN=*), INTENT(IN) :: title

  INTEGER, INTENT(IN) :: pltopt     ! variavle plot option (0/1/2/3)
  INTEGER, INTENT(IN) :: slicopt

  INTEGER, INTENT(IN)    :: mnsize  ! maximum m*n among all processors
  REAL,    INTENT(INOUT) :: plota(mnsize), plotx(mnsize), ploty(mnsize)
  INTEGER, INTENT(INOUT) :: iwrk(mnsize)
  REAL,    INTENT(INOUT) :: xwk(mnsize),   ywk(mnsize)
!
!-----------------------------------------------------------------------
!
!  Plotting control common blocks
!
!-----------------------------------------------------------------------
!
  INTEGER :: layover
  COMMON /laypar/ layover

  INTEGER :: ovrobs,obsset,obscol,obs_marktyp
  REAL :: obs_marksz,obs_valsz
  COMMON /obspar/ ovrobs,obsset,obscol,obs_marktyp,                    &
                  obs_marksz,obs_valsz

  INTEGER :: ovrstaopt
  INTEGER :: ovrstam,staset,ovrstan,ovrstav,stacol,markprio,wrtstax
  INTEGER :: nsta_typ,sta_typ(30),sta_marktyp(30),sta_markcol(30)
  REAL    :: sta_marksz(30),wrtstad
  CHARACTER (LEN=256) :: stalofl
  COMMON /sta_par/ ovrstaopt,ovrstam,staset,ovrstan,ovrstav,stacol,     &
         markprio,nsta_typ,sta_typ,sta_marktyp,sta_markcol,             &
         sta_marksz,stalofl,wrtstax,wrtstad

  REAL :: ctinc,ctmin,ctmax,vtunt  ! contour interval and vector unit
  COMMON /incunt/ ctinc,ctmin,ctmax,vtunt

  INTEGER :: icolor,icolor1,lbcolor,trcolor                ! required color
  COMMON /recolor/icolor,icolor1,lbcolor,trcolor

  INTEGER :: flag
  INTEGER :: xfont   ! the font of character
  INTEGER :: haxisu, vaxisu
  INTEGER :: lbaxis, tickopt, axlbfmt
  REAL    :: hmintick,vmajtick,vmintick,hmajtick
  COMMON /var_par/ xfont,haxisu,vaxisu,lbaxis,tickopt,hmintick,         &
          vmajtick, vmintick,hmajtick,axlbfmt

  REAL :: yxratio
  COMMON /yratio/ yxratio       ! the scaling factor the y/x ratio.

  INTEGER :: ntitle,titcol, nxpic, nypic, wpltime
  REAL    :: titsiz
  CHARACTER (LEN=256) :: ptitle(3), footer_l, footer_c, footer_r

  COMMON /titpar1/ptitle, footer_l, footer_c, footer_r
  COMMON /titpar2/ntitle,titcol,wpltime, nxpic, nypic
  COMMON /titpar3/titsiz
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=300) :: ch1, ch

  INTEGER :: istatus

  INTEGER :: i,j
  REAL :: cl(900)       ! contour levels
  REAL :: pl,pr,pb,pt   ! plot space left, right, bottom, top coordinate
  REAL :: px,py         ! plot space left-right length and up-down height
  REAL :: pxc,pyc       ! plot space left-right center and
                        !            up-down    center
  REAL :: xs,ys         ! real space left-right length and up-down height
  REAL :: zinc          ! contour interval
  REAL :: zmin,zmax     ! max and min of data array
  INTEGER :: ncl,mode1

  REAL :: zlevel
  COMMON/sliceh/zlevel

  INTEGER :: timeovr
  COMMON /timover/ timeovr

  REAL :: lblmag, ctrlbsiz, axlbsiz
  COMMON /labmag/ lblmag, ctrlbsiz, axlbsiz
  REAL :: xfinc

  INTEGER :: col_table,pcolbar
  COMMON /coltable/col_table,pcolbar

  INTEGER :: LEN0,len1

  CHARACTER (LEN=12) :: varname
  COMMON /varplt1/ varname

  CHARACTER (LEN=150) :: f_ch

  INTEGER :: setcontopt, setcontnum
  REAL :: setconts(maxunevm,maxuneva)
  COMMON /setcon_par/setcontopt,setcontnum,setconts
  INTEGER :: ncont
  REAL :: tcont(maxunevm)

  INTEGER :: wrtflag
  CHARACTER (LEN=25) :: timestring
  CHARACTER (LEN=90) :: levlab
  CHARACTER (LEN=50) :: timelab
  COMMON /timelev/wrtflag,timelab, levlab, timestring

  CHARACTER (LEN=80) :: prestr
  INTEGER :: preflag
  COMMON /preinfo/ prestr,preflag

  REAL :: x101, y101, x102,y102
  COMMON /slicev1/x101, y101, x102,y102

  REAL :: xttmp,yltmp,yttmp     !! local temporary variable
  REAL :: f_cputime,cpu1, cpu2
  DOUBLE PRECISION :: f_walltime,second1,second2
  REAL :: hatch_angle

  INTEGER :: missval_colind, missfill_opt    ! miss value color index
  COMMON /multi_value/ missfill_opt,missval_colind
  INTEGER :: missfill
  DATA missfill/0/

  INTEGER :: xnwpic_called
  COMMON /callnwpic/xnwpic_called

  INTEGER :: iclfrq
  INTEGER :: ctrlbfrq
  COMMON /clb_frq/ ctrlbfrq

!----------------------------------------------------------------------
!
! Message passing version
!
!---------------------------------------------------------------------

  INTEGER :: xpbgn,xpend,ypbgn,ypend
  COMMON /processors/ xpbgn, xpend, ypbgn, ypend

  INCLUDE 'mp.inc'

  INTEGER :: ii,jj
  INTEGER :: mm,nn  ! temporay varaible only useful for processor 0
  INTEGER :: mnn
  INTEGER            :: ierr, itags, itagr
  INTEGER, PARAMETER :: destination = 0
  INTEGER            :: source
!  CHARACTER(LEN=4) :: sourcechar

  REAL    :: clsaved(900)
  INTEGER :: nclsaved, nminctr, nmaxctr

  REAL, PARAMETER :: eps = 1.0E-6
  LOGICAL :: adjust_low, adjust_high
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ch  = ' '   ! initialize those string is necessary
  ch1 = ' '

!-----------------------------------------------------------------------
!
!  Check for adequate room in work array
!
!-----------------------------------------------------------------------
!
  second1= f_walltime()
  cpu1 = f_cputime()

  IF(myproc == 0) THEN

    WRITE(6,'(/,1x,a,a)') 'Plotting ',title

    IF( layover == 0  .OR. xnwpic_called == 0) THEN
      CALL xnwpic
      xnwpic_called = 1
      timeovr = 0             ! set overlayer terrain again
      wrtflag = 0             !
      preflag = 0
      prestr = levlab
      len1=LEN_TRIM(prestr)
      CALL strmin(prestr,len1)
    ELSE
      timeovr=1
      wrtflag = wrtflag + 1
    END IF
!
!-----------------------------------------------------------------------
!
!  Get plotting space variables
!
!-----------------------------------------------------------------------
!

    CALL xqpspc( pl, pr, pb, pt)
    px = pr - pl
    py = pt - pb
    pxc = (pr+pl)/2
    pyc = (pb+pt)/2

    xs = xr-xl
    ys = yt-yb
!
!-----------------------------------------------------------------------
!
!  Let the longest lenth determine size scaling of plot
!
!-----------------------------------------------------------------------
!
    IF( py/px >= (ys*yxratio)/xs ) THEN
      py = (ys*yxratio)/xs*px
      CALL xpspac(pl, pr, pyc-py/2, pyc+py/2 )
    ELSE
      px = xs/(ys*yxratio)*py
      CALL xpspac(pxc-px/2, pxc+px/2, pb, pt)
    END IF
!
!-----------------------------------------------------------------------
!
!  Set the real distance to plot distance scaling
!
!-----------------------------------------------------------------------
!
    CALL xmap( xl, xr, yb, yt )

    CALL xlbsiz( ctrlbsiz*(yt-yb)*lblmag )

  END IF  ! myproc == 0
!
!-----------------------------------------------------------------------
!
!  Find max and min of data array
!
!-----------------------------------------------------------------------
!
  missfill = 0
  zmax  = -999999999999.0
  zmin  =  999999999999.0
  DO j=1,n
    DO i=1,m
      IF(ABS(a(i,j)-(-9999.0)) < 1.0E-6) THEN
        missfill = 1
        CYCLE
      END IF
      zmax= MAX (zmax,a(i,j))
      zmin= MIN (zmin,a(i,j))
    END DO
  END DO

  IF (loc_x < xpbgn .OR. loc_x > xpend .OR. loc_y < ypbgn .OR. loc_y > ypend) THEN
    zmax = -9999.9
    zmin =  9999.9
  END IF
  CALL mpmax0 (zmax, zmin)      !? only inside xpbgn-xpend, ypbgn-ypend
  CALL mpmax0i(i,missfill)      ! Ensure missfill = 1 only when it
                                ! is 1 in all processors

  IF (missfill == 1 .AND. missfill_opt == 1)         &
      CALL fillmissval ( m,n,xl, xr, yb,yt )
!
!-----------------------------------------------------------------------
!
!  Find proper contour interval and then contour field
!  using ZXPLOT routine xconta
!
!-----------------------------------------------------------------------

  ncont = 0
  ncl   = 1
  cl(:) = 0.0

  IF( zmax-zmin > 1.0E-20 ) THEN
!
!-----------------------------------------------------------------------
!
!    Check to see if user defined contour levels is available for the
!    current variable.
!
!-----------------------------------------------------------------------

    IF(myproc == 0) CALL xcolor(lbcolor)

    CALL get_contour (ncont, tcont)
    iclfrq = ctrlbfrq
    IF(setcontopt > 0 .AND. ncont > 0) THEN
      ch1(1:11)=' contours: '
      len0=11
      DO i =1,ncont
        CALL xrch1(tcont(i),f_ch,len1)
        WRITE(ch1(len0+1:),'(2a)') f_ch(1:len1), ', '
        len0=len0+len1+2
      END DO
      DO i=1,ncont
        cl(i)=tcont(i)
      END DO
      ncl = ncont
      mode1 = 4
      iclfrq = 1
      GO TO 150
    END IF

    IF( ctinc == 0.0) THEN
      cl(2)=cl(1)+ xfinc(zmax-zmin)/2
      IF(cl(2)-cl(1) == 0.0) cl(2)=cl(1)+1.0
      nminctr = 8
      nmaxctr = 20
      !CALL xnctrs( 8,20)
      mode1=1
    ELSE IF ( ctinc == -9999.) THEN
      CALL set_interval(zmin,zmax,ctmin,ctmax,cl)
      ctinc = cl(2)-cl(1)
      zinc = ctinc
      nminctr = 8
      nmaxctr = 20
      !CALL xnctrs( 8,20)
      mode1=1
    ELSE
      cl(2)=cl(1)+ctinc
      ! Not one, as one ends up giving a division by zero later.
      nminctr = 2
      nmaxctr = 900
      !CALL xnctrs(1,900)
      mode1=1
    END IF

    ! new subroutine call for MPI mode. NOTE that mode1 reset to 4
    IF( mp_opt > 0 .AND. mode1 == 1) THEN
    	adjust_low  = .FALSE.
    	adjust_high = .FALSE.
      IF(ABS(ctmax-ctmin) < eps) THEN
        ctmax = zmax
        ctmin = zmin
        adjust_low  = .TRUE.
        adjust_high = .TRUE.
      ELSE IF (ctmax < -9990) THEN
        ctmax = zmax
        adjust_high = .TRUE.
      ELSE IF (ctmin < -9990) THEN
        ctmin = zmin
        adjust_low = .TRUE.
      END IF

      IF (myproc == 0) THEN
        IF (.NOT. adjust_low) THEN
          cl(2) = cl(2)-cl(1)+ctmin
          cl(1) = ctmin
        END IF
        CALL setcontr(ctmin,ctmax,nminctr,nmaxctr,cl,ncl)
        IF (cl(1) > zmin .AND. adjust_low) THEN
          DO i = ncl,1,-1
            cl(i+1) = cl(i)
          END DO
          cl(1) = zmin
          ncl   = ncl+1
        END IF
        IF (cl(ncl) < zmax .AND. adjust_high) THEN
          ncl = ncl+1
          cl(ncl) = zmax
        END IF
      END IF
!      CALL mpupdatei(ncl,1)
!      CALL mpupdater(cl, 900)
      mode1 = 4
    END IF

    150     CONTINUE

    nclsaved       = ncl
    clsaved(1:ncl) = cl(1:ncl)
    CALL xnctrs(nminctr,nmaxctr)

    zinc = cl(2)-cl(1)

!-----------------------------------------------------------------------
!
!  Plot contour or color filled contour fields
!
!-----------------------------------------------------------------------
    IF(myproc == 0) THEN

      CALL xwindw(xl, xr, yb, yt)

      CALL xctrlim(ctmin,ctmax)
      CALL xclfrq(iclfrq)

    END IF  ! myproc == 0

    ii = 0
    DO j = 1,n,jstep    ! to pack data to plota, to be passed to processor 0
      DO i = 1,m,istep
        ii = ii+1
        plota(ii) = a(i,j)
        plotx(ii) = x(i,j)
        ploty(ii) = y(i,j)
      END DO
      IF ( MOD(m-1,istep) /= 0) THEN   ! Plot the last index MUST
        ii = ii+1
        plota(ii) = a(m,j)
        plotx(ii) = x(m,j)
        ploty(ii) = y(m,j)
      END IF
    END DO

    IF ( MOD(n-1,jstep) /= 0) THEN     ! Must plot the last index
      j = n
      DO i = 1,m,istep
        ii = ii+1
        plota(ii) = a(i,j)
        plotx(ii) = x(i,j)
        ploty(ii) = y(i,j)
      END DO
      IF ( MOD(m-1,istep) /= 0) THEN
        ii = ii+1
        plota(ii) = a(m,j)
        plotx(ii) = x(m,j)
        ploty(ii) = y(m,j)
      END IF
    END IF

    mm = (m-1)/istep+1      ! mm*nn is the valid data to be ploted, this assignmen
    nn = (n-1)/jstep+1      ! is for processor 0 only. All other processors will
    mnn = ii             ! will pass their dimensions to processor 0 later
    IF (MOD(m-1,istep) /= 0) mm = mm + 1
    IF (MOD(n-1,jstep) /= 0) nn = nn + 1

    IF (mnn /= mm*nn) THEN
      WRITE(6,'(1x,a,4(I7,a))')                                              &
      'ERROR: Wrong size for arrays is detected in ctr2d, mnsize = ',mnn, &
      ', but it should be mm*nn = ',mm,'*',nn,' = ',mm*nn,'.'
      CALL arpsstop('Wrong size detected in ctr2d.',1)
    END IF

    DO jj = ypbgn,ypend
      DO ii = xpbgn, xpend

        source = (ii+(jj-1)*nproc_x-1)
        IF (source == 0) GOTO 600

        CALL inctag
        IF (myproc == source ) THEN
          itags = gentag + 4
          CALL mpsendi(mm,1,destination,itags,ierr)
          itags = gentag + 5
          CALL mpsendi(nn,1,destination,itags,ierr)

          itags = gentag
          CALL mpsendr(plota,mnn,destination,itags,ierr)
          itags = gentag + 1
          CALL mpsendr(plotx,mnn,destination,itags,ierr)
          itags = gentag + 2
          CALL mpsendr(ploty,mnn,destination,itags,ierr)
        END IF

        IF (myproc == 0) THEN
          itagr = gentag + 4
          CALL mprecvi(mm,1,source,itagr,ierr)
          itagr = gentag + 5
          CALL mprecvi(nn,1,source,itagr,ierr)

          itagr = gentag
          CALL mprecvr(plota,mm*nn,source,itagr,ierr)
          itagr = gentag + 1
          CALL mprecvr(plotx,mm*nn,source,itagr,ierr)
          itagr = gentag+2
          CALL mprecvr(ploty,mm*nn,source,itagr,ierr)

          ncl = nclsaved
          cl(1:nclsaved) = clsaved(1:nclsaved)
        END IF

        600 CONTINUE

        !WRITE(sourcechar,'(I04)') source

        IF (myproc == 0) THEN

          !CALL xpscmnt('Begin plotting processor :: '//sourcechar)

          IF(pltopt == 1) THEN
            CALL xctrclr(icolor, icolor)
            CALL xconta(plota,plotx,ploty,iwrk,mm,mm,nn,cl,ncl,mode1)
          ELSE IF( pltopt == 2) THEN
            CALL xctrclr(icolor, icolor1)
            CALL xcolfil(plota,plotx,ploty,iwrk,xwk,ywk,mm,mm,nn,cl,ncl,mode1)
            CALL xchmag(0.025*sqrt(px*py))
            CALL xcpalet(pcolbar)
          ELSE IF(pltopt == 4) THEN
            CALL xctrclr(icolor, icolor1)
            CALL xconta(plota,plotx,ploty,iwrk,mm,mm,nn,cl,ncl,mode1)
          ELSE IF(pltopt == 5) THEN
            CALL xctrclr(icolor, icolor1)
            CALL xcolfil(plota,plotx,ploty,iwrk,xwk,ywk,mm,mm,nn,cl,ncl,mode1)
            CALL xchmag(0.025*sqrt(px*py))
            CALL xcpalet(pcolbar)
            CALL xctrclr(lbcolor, lbcolor)
            CALL xconta(plota,plotx,ploty,iwrk,mm,mm,nn,cl,ncl,mode1)
          ELSE IF(pltopt == 6) THEN
            CALL xctrclr(icolor, icolor)
            CALL xdhtch(0.003)
            CALL xctrclr(icolor, icolor)
            ncl = 2
            mode1 = 4
            cl(1) = ctmin
            cl(2) = ctmax
            CALL xclfrq(1)
            CALL xhilit(0)
            CALL xconta(plota,plotx,ploty,iwrk,mm,mm,nn,cl,ncl,mode1)
            CALL xhilit(1)

            hatch_angle = 45.0
            CALL xdhtch(0.004)
            CALL xhatcha(plota,plotx,ploty,xwk,ywk,mm,mm,nn,ctmin,1.0E10,hatch_angle)

            CALL xdhtch(0.002)
            CALL xhatcha(plota,plotx,ploty,xwk,ywk,mm,mm,nn,ctmax,1.0E10,hatch_angle)

          END IF

          CALL xclfrq(2)

          !CALL xpscmnt('End plotting processor ::'//sourcechar)

        END IF  ! myproc == 0

        CALL mpbarrier         ! sync the processors

      END DO
    END DO

  ELSE
    cl(2)=1.0
    ncl=2
  END IF  ! zmax-zmin > 1.0E-20

  IF(ctinc == 0.0) THEN
    zinc = cl(2) - cl(1)
  ELSE
    zinc = ctinc
  END IF
!
!-----------------------------------------------------------------------
!
!  Plot map, boxes and polygons.
!
!-----------------------------------------------------------------------
!
  IF(myproc == 0) CALL pltextra(slicopt, pltopt)
!
!-----------------------------------------------------------------------
!
!  Plot terrain etc.
!
!-----------------------------------------------------------------------
!
    ii = 0
    DO j = 1,n,jstep             ! again useful for processor 0 only
      DO i = 1,m,istep
        ii = ii+1
        plota(ii) = hterain(i,j)
        plotx(ii) = x(i,j)
        ploty(ii) = y(i,j)
      END DO
      IF ( MOD(m-1,istep) /= 0) THEN
        ii = ii+1
        plota(ii) = a(m,j)
        plotx(ii) = x(m,j)
        ploty(ii) = y(m,j)
      END IF
    END DO

    IF ( MOD(n-1,jstep) /= 0) THEN
      j = n
      DO i = 1,m,istep
        ii = ii+1
        plota(ii) = a(i,j)
        plotx(ii) = x(i,j)
        ploty(ii) = y(i,j)
      END DO
      IF ( MOD(m-1,istep) /= 0) THEN
        ii = ii+1
        plota(ii) = a(m,j)
        plotx(ii) = x(m,j)
        ploty(ii) = y(m,j)
      END IF
    END IF

    mm = (m-1)/istep + 1
    nn = (n-1)/jstep + 1
    mnn = ii
    IF (MOD(m-1,istep) /= 0) mm = mm + 1
    IF (MOD(n-1,jstep) /= 0) nn = nn + 1

    DO jj = ypbgn,ypend
      DO ii = xpbgn, xpend

        source = (ii+(jj-1)*nproc_x-1)
        IF (source == 0) GOTO 602

        CALL inctag
        IF (myproc == source ) THEN
          itags = gentag + 4
          CALL mpsendi(mm,1,destination,itags,ierr)
          itags = gentag + 5
          CALL mpsendi(nn,1,destination,itags,ierr)


          itags = gentag
          CALL mpsendr(plota,mnn,destination,itags,ierr)
          itags = gentag + 1
          CALL mpsendr(plotx,mnn,destination,itags,ierr)
          itags = gentag + 2
          CALL mpsendr(ploty,mnn,destination,itags,ierr)
        END IF

        IF (myproc == 0) THEN
          itagr = gentag + 4
          CALL mprecvi(mm,1,source,itagr,ierr)
          itagr = gentag + 5
          CALL mprecvi(nn,1,source,itagr,ierr)

          itagr = gentag
          CALL mprecvr(plota,mm*nn,source,itagr,ierr)
          itagr = gentag + 1
          CALL mprecvr(plotx,mm*nn,source,itagr,ierr)
          itagr = gentag+2
          CALL mprecvr(ploty,mm*nn,source,itagr,ierr)
        END IF

        602 CONTINUE

        IF (myproc == 0) THEN

!-----------------------------------------------------------------------
!
!  Terrain outline in vertical slices.
!
!-----------------------------------------------------------------------
          IF(slicopt == 2 .OR. slicopt == 3 .OR. slicopt == 5 .OR.  &
             slicopt == 10 .OR. slicopt == 11) THEN

            CALL xcolor(trcolor)
            CALL xthick(3)
            CALL xpenup( plotx(1), ploty(1)-0.5*(ploty(mm+1)-ploty(1)) )
            DO i=2,mm
              CALL xpendn( plotx(i), ploty(i)-0.5*(ploty(i+mm)-ploty(i)) )
            END DO
            CALL xthick(1)

          END IF   ! slicopt
!
!-----------------------------------------------------------------------
!
!  Overlay terrain contour if required in x-y level
!  or Plot terrain outline in slice zlevel
!
!-----------------------------------------------------------------------
!
           IF ( timeovr == 0 ) CALL plttrn(plota,plotx,ploty,mm,nn,     &
                                           slicopt,iwrk,xwk,ywk)

        END IF  ! myproc == 0

        CALL mpbarrier         ! sync the processors

      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Plot station labels
!
!-----------------------------------------------------------------------
!
  IF (ovrstaopt == 1 .AND. (wrtstax == 1 .OR. staset == 1) ) THEN
    ii = 0
    DO j = 1,n,jstep
      DO i = 1,m,istep
        ii = ii+1
        plota(ii) = a(i,j)
        plotx(ii) = x(i,j)
        ploty(ii) = y(i,j)
      END DO
      IF ( MOD(m-1,istep) /= 0) THEN
        ii = ii+1
        plota(ii) = a(m,j)
        plotx(ii) = x(m,j)
        ploty(ii) = y(m,j)
      END IF
    END DO

    IF ( MOD(n-1,jstep) /= 0) THEN
      j = n
      DO i = 1,m,istep
        ii = ii+1
        plota(ii) = a(i,j)
        plotx(ii) = x(i,j)
        ploty(ii) = y(i,j)
      END DO
      IF ( MOD(m-1,istep) /= 0) THEN
        ii = ii+1
        plota(ii) = a(m,j)
        plotx(ii) = x(m,j)
        ploty(ii) = y(m,j)
      END IF
    END IF

    mm = (m-1)/istep + 1
    nn = (n-1)/jstep + 1
    mnn = ii
    IF (MOD(m-1,istep) /= 0) mm = mm + 1
    IF (MOD(n-1,jstep) /= 0) nn = nn + 1


    DO jj = ypbgn,ypend
      DO ii = xpbgn, xpend

        source = (ii+(jj-1)*nproc_x-1)
        IF (source == 0) GOTO 603

        CALL inctag
        IF (myproc == source ) THEN
          itags = gentag + 4
          CALL mpsendi(mm,1,destination,itags,ierr)
          itags = gentag + 5
          CALL mpsendi(nn,1,destination,itags,ierr)


          itags = gentag
          CALL mpsendr(plota,mnn,destination,itags,ierr)
          itags = gentag + 1
          CALL mpsendr(plotx,mnn,destination,itags,ierr)
          itags = gentag + 2
          CALL mpsendr(ploty,mnn,destination,itags,ierr)
        END IF

        IF (myproc == 0) THEN
          itagr = gentag + 4
          CALL mprecvi(mm,1,source,itagr,ierr)
          itagr = gentag + 5
          CALL mprecvi(nn,1,source,itagr,ierr)

          itagr = gentag
          CALL mprecvr(plota,mm*nn,source,itagr,ierr)
          itagr = gentag + 1
          CALL mprecvr(plotx,mm*nn,source,itagr,ierr)
          itagr = gentag+2
          CALL mprecvr(ploty,mm*nn,source,itagr,ierr)
        END IF

        603 CONTINUE

        IF (myproc == 0) THEN

          IF( wrtstax == 1 .AND. (timeovr == 0 .OR.                     &
                                  (timeovr== 1 .AND. pltopt == 2)) .AND.&
              (slicopt == 2  .OR. slicopt == 3 .OR. slicopt == 5 .OR.   &
               slicopt == 10 .OR. slicopt == 11) ) THEN
             CALL xchmag(0.025*px * lblmag)
             flag=1
             CALL pltsta(plota,plota,plotx,ploty,mm,nn,flag,slicopt)
           END IF

           IF( staset == 1 .AND. (ovrstam == 1 .OR. ovrstan == 1        &
                                  .OR. ovrstav == 1)   .AND.            &
               (slicopt == 1 .OR. slicopt == 4 .OR. slicopt == 6 .OR.   &
                slicopt == 7 .OR. slicopt == 8 .OR. slicopt == 9) .AND. &
               (timeovr == 0 .OR. (timeovr == 1.AND.pltopt == 2) )) THEN
             CALL xchmag(0.025*px * lblmag)
             flag=0
             CALL pltsta(plota,plota,plotx,ploty,mm,nn,flag,slicopt)
           END IF

        END IF  ! myproc == 0

        CALL mpbarrier         ! sync the processors

      END DO
    END DO

  END IF  ! ovrstaopt == 1

  IF (myproc == 0)  CALL xwdwof
!
!-----------------------------------------------------------------------
!
!  Plot observations
!
!-----------------------------------------------------------------------
!
  IF(obsset == 1 .AND. ovrobs == 1 .AND.  &
     (slicopt == 1 .OR. slicopt == 4 .OR. slicopt == 6 .OR. &
      slicopt == 7 .OR. slicopt == 8 .OR. slicopt ==  9) ) THEN

    IF (myproc == 0) THEN
      CALL xchmag(0.025*px * lblmag)
      CALL pltobs(1)
    END IF

    obsset=0

  END IF
!
!-----------------------------------------------------------------------
!
!  Plot axes with tick marks
!
!-----------------------------------------------------------------------
!
  IF(myproc == 0) THEN

    CALL pltaxes(slicopt,dx,dy)

    IF(ntitle>0 .AND. nxpic==1 .AND. nypic==1 .AND. timeovr == 0 ) THEN
      CALL xcolor(titcol)
      CALL xchmag(0.025*px * titsiz)
      DO i=1,ntitle
        LEN0=256
        CALL strlnth(ptitle(i),LEN0)
        CALL xchori(0.)
        CALL xcharc( xl+xs/2,yt+(0.25-(i-1)*0.06)*ys*px/py,ptitle(i)(1:LEN0))
      END DO
      CALL xcolor(lbcolor)
    END IF

    CALL xchmag( 0.025*px * lblmag )

    ! plot time and level label
    IF ( layover < 1 ) THEN
      IF(levlab /= ' ') THEN
        len1=LEN_TRIM(levlab)
        CALL strmin(levlab,len1)
        CALL xcharc((xl+xr)*0.5,yt+0.06*ys*px/py, levlab(1:len1))
        preflag = 1
      END IF
      len1=LEN_TRIM(timelab)
      CALL strmin(timelab,len1)
      CALL xcharc((xl+xr)*0.5,yt+0.015*ys*px/py,                            &
                        timestring(1:25)//' '//timelab(1:len1))
    END IF

    IF(preflag == 0 .AND. levlab /= ' ') THEN
      len1=LEN_TRIM(levlab)
      CALL strmin(levlab,len1)
      CALL xcharc((xl+xr)*0.5,yt+0.06*ys*px/py, levlab(1:len1))
      preflag = 1
    END IF

    LEN0 = LEN_TRIM(title)
    CALL strmin(title, LEN0)

    IF( title(LEN0:LEN0) == ')' ) LEN0 = max(1,LEN0-1)

    IF(pltopt == 2) THEN
      WRITE(f_ch, '(2a)') title(1:LEN0), ', Shaded)'
    ELSE IF( pltopt == 5 ) THEN
      WRITE(f_ch, '(2a)') title(1:LEN0), ', Shaded/Contour)'
    ELSE
      WRITE(f_ch, '(2a)' )title(1:LEN0), ', contour)'
    END IF

! if first levlab is not equal second levlab then attatch levlab on f_ch
    LEN0=LEN_TRIM(f_ch)
    CALL strmin(f_ch, LEN0)
    len1=LEN_TRIM(levlab)
    CALL strmin(levlab,len1)

    IF(pltopt == 1) CALL xcolor(icolor)
    ! plot variable name
    IF (preflag == 1 .AND. prestr /= levlab .AND. prestr /= ' '           &
          .AND.layover /= 0 .AND. levlab /= ' ') THEN
      CALL xchmag( 0.018*px * lblmag )
    ELSE
      CALL xchmag( 0.028*px * lblmag )
    END IF
    IF(prestr(1:1) == ' ' .AND. layover /= 0 ) prestr=levlab
                                              ! save for next time use

    IF(lbaxis == 1 ) THEN
      IF( wrtstax == 0) THEN
        yltmp = 0.08
      ELSE
        yltmp = 0.14
      END IF
    ELSE
      yltmp = 0.12
    END IF

    LEN0=LEN_TRIM(f_ch)
    CALL strmin(f_ch,LEN0)

    CALL xchmag( 0.025*px * lblmag )
    CALL xcolor(lbcolor)

    xttmp = xl-0.05*(xr-xl)
    yttmp = yb-(yltmp+wrtflag*0.030)*ys*px/py
    CALL xcharl(xttmp, yttmp, f_ch(1:LEN0))

    IF ( pltopt == 1 .OR. pltopt == 3 .OR. pltopt == 4 .OR. pltopt == 5 ) THEN
      IF ( ABS(zmin-zmax) <= 1.e-15 .OR. ncont > 0)  THEN
        WRITE(ch,'(2(a,G9.3E2))' ) 'Min=',zmin,' Max=',zmax
      ELSE
        WRITE(ch,'(3(a,G10.4E2))') 'Min=',zmin,' Max=',zmax,' inc=',zinc
      END IF
    ELSE IF( pltopt == 2 ) THEN
        WRITE(ch,'(2(a,G9.3E2))' ) 'Min=', zmin,' Max=',zmax
    END IF

    xttmp = xr+0.05*(xr-xl)
    yttmp = yb-(yltmp+wrtflag*0.030)*ys*px/py
    LEN0=LEN_TRIM(ch)
    CALL strmin(ch,LEN0)
    CALL xcharr(xttmp, yttmp, ch(1:LEN0))
    IF (ncont > 1 .AND. (pltopt == 1 .OR. pltopt == 4) ) THEN
      wrtflag = wrtflag+1
      xttmp = xr+0.05*(xr-xl)
      yttmp = yb-(yltmp+wrtflag*0.030)*ys*px/py
      len1=LEN_TRIM(ch1)
      CALL strmin(ch1,len1)
      CALL xcharr(xttmp, yttmp, ch1(1:len1))
    END IF

!-----------------------------------------------------------------------
!
!  Plot additional text below the figure
!
!-----------------------------------------------------------------------

    CALL label2d(runname)

    CALL xpspac(pl, pr, pb, pt)  ! set frame back

  END IF  ! myproc == 0

  cpu2 = f_cputime()
  second2 = f_walltime()

!  write(6,*) '!!!!  total cpu time for one CTR2D  :',                   &
!             cpu2-cpu1,' PLOT:',varname

  RETURN
END SUBROUTINE ctr2d
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CTRINC                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE ctrinc( ctinc0, ctmin0, ctmax0 )

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    Set the contour interval for field to plotted by CTR2D.
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!
!  MODIFICATION HISTORY:
!    6/08/92  Added full documentation (K. Brewster)
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    ctinc0    Contour interval
!              If CTINC0 = 0.0, the interval is internally determined.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  REAL :: ctinc0,ctmin0,ctmax0
!
!-----------------------------------------------------------------------
!
!  Plotting control common blocks
!
!-----------------------------------------------------------------------
!
  REAL :: ctinc,ctmin,ctmax,vtunt   ! contour interval and vector unit
  COMMON /incunt/ ctinc,ctmin,ctmax,vtunt
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  ctinc = ctinc0
  ctmin = ctmin0
  ctmax = ctmax0

  RETURN
END SUBROUTINE ctrinc
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE LABEL2D                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE label2d(runname)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Plot certain text labels for VTR2D and CTR2d.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!
!  MODIFICATION HISTORY:
!  Taked from former CTR2D.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    xl       Left bound of the physical domain
!    xr       Right bound of the physical domain
!    yb       Bottom bound of the physical domain.
!    yt       Top bound of the physical domain.
!
!    runname  character string describing the model run
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  CHARACTER (LEN=*) :: runname
  INTEGER :: layover
  COMMON /laypar/ layover

  INTEGER :: icolor,icolor1,lbcolor,trcolor                ! required color
  COMMON /recolor/icolor,icolor1,lbcolor,trcolor
!
  INTEGER :: xfont   ! the font of character
  INTEGER :: haxisu, vaxisu
  INTEGER :: lbaxis
  INTEGER :: tickopt
  INTEGER :: axlbfmt
  REAL :: hmintick,vmajtick,vmintick,hmajtick
  COMMON /var_par/ xfont,haxisu,vaxisu,lbaxis,tickopt,hmintick,         &
          vmajtick, vmintick,hmajtick,axlbfmt
!
  INTEGER :: ntitle,titcol, nxpic, nypic, wpltime
  REAL :: titsiz
  CHARACTER (LEN=256) :: ptitle(3), footer_l, footer_c, footer_r

  COMMON /titpar1/ptitle, footer_l, footer_c, footer_r
  COMMON /titpar2/ntitle,titcol,wpltime, nxpic, nypic
  COMMON /titpar3/titsiz

  REAL :: xl,xr,yb,yt
  REAL :: pl, pr, pb, pt
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: nopic
  REAL :: xlimit, ylimit, rotang
  INTEGER :: nhpic, nvpic,ifont

  INTEGER :: ovrtrn ,trnplt       ! overlay terrain option (0/1)
  REAL :: trninc,trnmin, trnmax   ! terrain interval minimum, maximum
  REAL :: ztmin,ztmax
  COMMON /trnpar/ trnplt,ovrtrn,trninc,trnmin, trnmax,ztmin,ztmax
  INTEGER :: timeovr
  COMMON /timover/ timeovr

  REAL :: lblmag, ctrlbsiz, axlbsiz
  COMMON /labmag/ lblmag, ctrlbsiz, axlbsiz
!
  INTEGER :: col_table,pcolbar
  COMMON /coltable/col_table,pcolbar
!
  CHARACTER (LEN=24)  :: tzstring
  CHARACTER (LEN=24)  :: tz
  CHARACTER (LEN=256) :: datetimestr

  INTEGER :: lnblnk, len1, len2, len3
  CHARACTER (LEN=256) :: string_l, string_c, string_r

  CHARACTER (LEN=8) :: tzone
  CHARACTER (LEN=10) :: cur_time
  CHARACTER (LEN=8) :: cur_date
  INTEGER :: t_values(8)

  REAL :: ytmp, hch
  REAL :: px, py
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  CALL xqmap(xl,xr,yb,yt)

  CALL xqpspc( pl, pr, pb, pt)
  px = pr - pl
  py = pt - pb

  CALL xcolor(lbcolor)

  CALL xqnpic(nopic)
  CALL xqspac(nhpic, nvpic, rotang, xlimit, ylimit)

  CALL xqpspc( pl, pr, pb, pt)

  CALL xchmag( 0.021*(pr-pl) * lblmag )

  IF(timeovr == 0) THEN
    IF(nopic == nhpic*(nvpic-1)+1 ) THEN

      IF ( wpltime == 1) THEN

        CALL date_and_time(cur_date,cur_time,tzone,t_values)

        IF(t_values(4) == 0) THEN
          tzstring = ' UTC'
        ELSE
          tzstring = ' Local Time'
        END IF

        WRITE (datetimestr,999) 'Plotted ',                             &
            t_values(1),t_values(2),t_values(3),                        &
            t_values(5),t_values(6),tzstring
        999        FORMAT (a, i4.4,'/',i2.2,'/',i2.2,' ',i2.2,':',i2.2,a)
      END IF

      IF ( footer_l == ' ') THEN
        string_l = 'ARPSPLT/ZXPLOT '
      ELSE
        string_l = footer_l
      END IF

      IF( footer_c == '  ') THEN
        string_c = runname
      ELSE
        string_c = footer_c
      END IF

      IF(wpltime == 1 ) THEN
        string_r = datetimestr(:lnblnk(datetimestr))
      ELSE
        string_r = footer_r
      END IF

      CALL xqcfnt(ifont)
      CALL xcfont(xfont)

      ytmp = 0.29

      CALL xqchsz(hch)

      IF ( layover < 1) THEN

        len1=LEN_TRIM(string_l)
        CALL strmin(string_l, len1)
        len2=LEN_TRIM(string_c)
        CALL strmin(string_c, len2)
        len3=LEN_TRIM(string_r)
        CALL strmin(string_r, len3)

        CALL xcharc(xl+0.5*(xr-xl),                                     &
             yb-(ytmp+layover*0.03)*(yt-yb)*px/py,                      &
             string_l(1:len1)//'  '//string_c(1:len2)//'  '//           &
             string_r(1:len3))

      END IF
      CALL xcfont(ifont)
    END IF
    timeovr=1
  END IF

  RETURN
END SUBROUTINE label2d
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE VTR3D                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE vtr3d(u,v,w, x,y,z, xw,xe,dx, ys,yn,dy, zb,zt,dz,            &
           nx,ibgn,iend,ist, ny,jbgn,jend,jst, nz,kbgn,kend,kst,        &
           kslice, jslice, islice, label,time, runname, factor,         &
           slicopt,n,xp,yp,zp,u1,v1,u2,v2,w2,                           &
           tem1,tem2,tem3,tem4,                                         &
           tem5,tem6,hterain)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    Plot vector fields in 2-d slices
!
!  AUTHOR: Ming Xue
!
!  MODIFICATION HISTORY:
!    6/08/92  Added full documentation (K. Brewster)
!
!  3/25/96 (K. Brewster)
!    Added variables isize,jsize,ksize
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    u        3-dimensional array of u wind components (m/s)
!    v        3-dimensional array of v wind components (m/s)
!    w        3-dimensional array of w wind components (m/s)
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in physical space (m)
!
!    xw       value of x for first i grid point to plot
!    xe       value of x for last i grid point to plot
!    ys       value of y for first j grid point to plot
!    yn       value of y for last j grid point to plot
!    zb       value of z for first k grid point to plot
!    zt       value of z for last k grid point to plot
!
!    nx       first dimension of b
!    ibgn     index of first i grid point to plot
!    iend     index of last  i grid point to plot
!
!    ny       second dimension of b
!    jbgn     index of first j grid point to plot
!    jend     index of last  j grid point to plot
!
!    nz       third dimension of b
!    kbgn     index of first k grid point to plot
!    kend     index of last  k grid point to plot
!
!    ist      step size in x direction
!    jst      step size in y direction
!    kst      step size in z direction
!
!    time     time of data in seconds
!
!    kslice   k index of plane for slicopt=1 x-y slice
!    jslice   j index of plane for slicopt=2 x-z slice
!    islice   i index of plane for slicopt=1 y-z slice
!
!    runname  character string decribing run
!
!    factor   scaling factor for winds
!             V*factor wind vectors are plotted
!
!    slicopt  slice orientation indicator
!             = 1, x-y slice of at k=kslice is plotted.
!             = 2, x-z slice of at j=jslice is plotted.
!             = 3, y-z slice of at i=islice is plotted.
!             = 4, horizontal slice at z index islice is plotted.
!             = 5, xy-z cross section of wind islice is plotted.
!             = 6, data field on constant p-level is plotted.
!             = 0, all of the three slices above are plotted.
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
!    tem4     Temporary work array.
!    tem5     Temporary work array.
!
!  (These arrays are defined and used locally (i.e. inside this
!   subroutine), they may also be passed into routines called by
!   this one. Exiting the call to this subroutine, these temporary
!   work arrays may be used for other purposes therefore their
!   contents overwritten. Please examine the usage of work arrays
!   before you alter the code.)
!
!   pp01      The pressure (mb) value at the specific p-level
!   ercpl     reciprocal of exponent
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz
  INTEGER :: n

  REAL :: u(nx,ny,nz)
  REAL :: v(nx,ny,nz)
  REAL :: w(nx,ny,nz)

  REAL :: x(nx,ny,nz)
  REAL :: y(nx,ny,nz)
  REAL :: z(nx,ny,nz)

  REAL :: u1(nx,ny),v1(nx,ny)
  REAL :: u2(n,nz),v2(n,nz),w2(n,nz),zp(n,nz)
  REAL :: xp(n),yp(n)

  REAL :: hterain(nx,ny)           ! The height of the terrain.

  INTEGER :: kslice,jslice,islice
  CHARACTER (LEN=*) :: runname
  CHARACTER (LEN=*) :: label

  REAL :: xw,xe,dx,ys,yn,dy,zb,zt,dz
  INTEGER :: ibgn,iend,ist, jbgn,jend,jst, kbgn,kend,kst

  REAL :: time,factor
  INTEGER :: slicopt

  INTEGER :: iunits, itype
  COMMON /windvtr/iunits, itype

  CHARACTER (LEN=12) :: varname
  COMMON /varplt1/ varname

  REAL :: xw1,xe1,ys1,yn1
  COMMON /xuvpar/xw1,xe1,ys1,yn1

!
!-----------------------------------------------------------------------
!
!  Some constants
!
!-----------------------------------------------------------------------
!
  REAL :: pp01
  REAL, PARAMETER :: ercpl=0.3678794              ! exp(-1.0)
!
!-----------------------------------------------------------------------
!
!  Work arrays: tem1,tem2,tem3,tem4,tem5 of size at least
!          max( nx*ny, nx*nz, ny*nz).
!
!-----------------------------------------------------------------------
!
  REAL :: tem1(*),tem2(*),tem3(*),tem4(*),tem5(*)
  REAL :: tem6(*)
!
!-----------------------------------------------------------------------
!
!  Common blocks for plotting control parameters
!
!-----------------------------------------------------------------------
!
  REAL :: x01,y01                  ! the first  point of interpolation
  REAL :: x02,y02                  ! the second point of interpolation
  REAL :: zlevel                   ! the given height of the slice
  REAL :: sinaf,cosaf,dist,sqrtdxy
  COMMON /slicev/x01,y01,x02,y02,sinaf,cosaf,dist,sqrtdxy
  COMMON /sliceh/zlevel

  INTEGER :: ovrobs,obsset,obscol,obs_marktyp
  REAL :: obs_marksz,obs_valsz
  COMMON /obspar/ ovrobs,obsset,obscol,obs_marktyp,                    &
                  obs_marksz,obs_valsz

  INTEGER :: icolor,icolor1,lbcolor,trcolor        ! required color
  COMMON /recolor/icolor,icolor1,lbcolor,trcolor

  INTEGER :: trnplt                ! flag to plot terain (1 or 0)
  INTEGER :: ovrtrn         ! overlay terrain option (0/1)
  REAL :: trninc,trnmin, trnmax    ! terrain interval minimum, maximum
  REAL :: ztmin,ztmax
  COMMON /trnpar/ trnplt,ovrtrn,trninc,trnmin, trnmax,ztmin,ztmax

!
!-----------------------------------------------------------------------
!
!  Misc. local Variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,ij,ik,jk,istep,jstep,length,isize,jsize,ksize
  REAL :: uunit
  CHARACTER (LEN=9)   :: timhms
  CHARACTER (LEN=120) :: title

  INTEGER :: xfont   ! the font of character
  INTEGER :: haxisu, vaxisu
  INTEGER :: lbaxis
  INTEGER :: tickopt
  INTEGER :: axlbfmt
  REAL :: hmintick,vmajtick,vmintick,hmajtick
  COMMON /var_par/ xfont,haxisu,vaxisu,lbaxis,tickopt,hmintick,         &
          vmajtick, vmintick,hmajtick,axlbfmt
  CHARACTER (LEN=6) :: stem2
  CHARACTER (LEN=1) :: stem1

  INTEGER :: smooth
  COMMON /smoothopt/smooth

  INTEGER :: id

  REAL :: x_tmp
  COMMON /tmphc2/ x_tmp

  INTEGER :: wrtflag
  CHARACTER (LEN=90) :: levlab
  CHARACTER (LEN=50) :: timelab
  CHARACTER (LEN=25) :: timestring
  COMMON /timelev/wrtflag, timelab, levlab, timestring

  REAL :: tmpx, tmpy
  CHARACTER (LEN=20) :: distc
  REAL :: x101, y101, x102,y102
  COMMON /slicev1/x101, y101, x102,y102
  INTEGER :: llabel
  CHARACTER (LEN=120) :: label_copy

  INTEGER :: xpbgn,xpend,ypbgn,ypend
  COMMON /processors/ xpbgn, xpend, ypbgn, ypend

  INTEGER :: agl
  COMMON /agl_or_asl/ agl

  INTEGER :: idsize, jdsize, mnsize
  INTEGER :: tinds, tind1,tind2,tind3,tind4,tind5,tind6,tind7,tind8

  REAL, ALLOCATABLE :: u2_tr(:,:)
  REAL, ALLOCATABLE :: w2_tr(:,:)
  REAL, ALLOCATABLE :: xptr(:,:)
  REAL, ALLOCATABLE :: zptr(:,:)
  REAL, ALLOCATABLE :: temptr1(:,:),temptr2(:,:),temptr3(:,:)
  REAL, ALLOCATABLE :: temptr4(:,:),temptr5(:,:),temptr6(:,:)
  REAL, ALLOCATABLE :: temptr7(:,:),temptr8(:,:)
  REAL, ALLOCATABLE :: xptr1d(:),yptr1d(:)
  REAL :: xptr1,xptr2,dxptr,xtr,ytr
  INTEGER :: itrajc_start, itrajc_end, iptr
  INTEGER :: ireturn,j_bgn,j_end
  INTEGER :: j_bgn_last,j_end_last,istatus,kz
  REAL :: sinaf_ptr,cosaf_ptr, delta_x,delta_y

!----------------------------------------------------------------------
!
! Include files
!
!---------------------------------------------------------------------

  INCLUDE 'mp.inc'
  INCLUDE 'arpstrajc.inc'

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  isize=(iend-ibgn)+1
  jsize=(jend-jbgn)+1
  ksize=(kend-kbgn)+1

  idsize = isize            ! global maximum isize
  jdsize = jsize
  CALL mpmaxi(idsize)
  CALL mpmaxi(jdsize)

  mnsize = idsize*jdsize
  mnsize = MAX(mnsize,idsize*ksize,jdsize*ksize)

  tind1 = 1             ! reuse a 3d temporary array 'tem5' as several 2D
  tind2 = tind1+mnsize  ! arrays inside vtr2d
  tind3 = tind2+mnsize
  tind4 = tind3+mnsize
  tind5 = tind4+mnsize
  tind6 = tind5+mnsize
  tind7 = tind6+mnsize
  tind8 = tind7+mnsize

!  tinds = SIZE(tem6)
!  IF (tinds < 5*mnsize) THEN
!    WRITE(6,'(3a)') 'ERROR: temporary array tem6 is too small ',        &
!                    'inside vtr3d while plotting ',label
!    CALL arpsstop('Temporary array too small inside vtr3d.',1)
!  END IF

  label_copy = label
  llabel = 120
  CALL xstrlnth(label_copy, llabel)
  IF(myproc == 0)CALL xpscmnt('Start plotting '//label_copy(1:llabel))
!
!-----------------------------------------------------------------------
!
!  slicopt=1   Plot u-v field
!
!-----------------------------------------------------------------------
!
  CALL get_forecast_hms( time, timhms )

  WRITE(timelab,'(a,F8.1,a)') 'T=',time,' s ('//TRIM(timhms)//')'

  CALL get_time_string ( time, timestring,'Z  ',0 )

  IF ( slicopt == 2 .OR. slicopt == 3  .OR. slicopt == 5) THEN
    CALL cal_dist(haxisu,dx,dy,x01,y01,x02,y02,slicopt,                 &
                  tmpx,tmpy,distc)
  END IF

!
!-----------------------------------------------------------------------
!
!  Set up terrain, if needed.
!
!-----------------------------------------------------------------------
!
  IF(trnplt == 1 .OR.trnplt == 2 .OR. ovrtrn == 1)  THEN
    DO j=jbgn,jend
      DO i=ibgn,iend
        ij = i-ibgn+1 + (j-jbgn)*isize
        tem5(ij)=hterain(i,j)
      END DO
    END DO
  END IF

  IF( slicopt == 1 .OR. slicopt == 0 ) THEN

    k = kslice
    DO j=jbgn,jend
      DO i=ibgn,iend
        ij = i-ibgn+1 + (j-jbgn)*isize
        tem1(ij) = -9999.0
        tem2(ij) = -9999.0
        IF(u(i,j,k) /= -9999.0) tem1(ij)=u(i,j,k)*factor
        IF(v(i,j,k) /= -9999.0) tem2(ij)=v(i,j,k)*factor
        tem3(ij)=x(i,j,k)
        tem4(ij)=y(i,j,k)
      END DO
    END DO

    IF (k /= 2) THEN
      WRITE(levlab,'(''GRID LEVEL='',I3)')k
      WRITE(title,'(''U-V '',A)')label
    ELSE
      WRITE(levlab,'(''First level above ground (surface)'')')
      WRITE(title,'(''U-V '',A)') label
    END IF

    length = 120
    CALL strlnth( title, length )
    CALL strmin ( title, length)

    uunit = 10.0
    CALL xvmode(1)
    istep = ist
    jstep = jst

    DO i=1,smooth
      CALL smooth9pmv(tem1,isize,jsize,1,isize,1,jsize,tem6)
      CALL smooth9pmv(tem2,isize,jsize,1,isize,1,jsize,tem6)
    END DO

    CALL vtr2d(tem1,tem2,tem3,tem4, uunit, xw,xe,dx,ys,yn,dy,           &
               isize,istep,jsize,jstep,title(1:length),runname, 1,      &
               tem5,slicopt,mnsize,tem6(tind1),tem6(tind2),tem6(tind3), &
               tem6(tind4),tem6(tind5),tem6(tind6),tem6(tind7),         &
               tem6(tind8))
!
!-----------------------------------------------------------------------
!
!  slicopt=2   Plot u-w field
!
!-----------------------------------------------------------------------
!
  ELSE IF( slicopt == 2 .OR. slicopt == 0 ) THEN

    x_tmp = y(1,jslice,1)

    j = jslice

    j = j + (ypbgn-1)*(ny-3)
    dist = (j-1.5)*tmpy
    length=LEN_TRIM(distc)
    CALL strmin ( distc, length)
    WRITE(levlab,'(''X-Z Plane at Y='',F8.1,A)')dist,distc(1:length)

    IF(varname(1:6) == 'xuvplt') THEN
      xw1=xw
      xe1=xe
      ys1=ys
      yn1=yn
      id=4
      DO k=kbgn,kend
        DO i=ibgn,iend
          ik = i-ibgn+1 + (k-kbgn)*isize
          tem1(ik) = -9999.0
          tem2(ik) = -9999.0
          IF(u(i,jslice,k) /= -9999.0) tem1(ik)=u(i,jslice,k)*factor
          IF(v(i,jslice,k) /= -9999.0) tem2(ik)=v(i,jslice,k)*factor
          tem3(ik)=x(i,jslice,k)
          tem4(ik)=z(i,jslice,k)
        END DO
      END DO
      WRITE(title,'(''U-V '',A)')label
    ELSE
      !CALL set_vertical_factor ( (zt-zb)/(xe-xw) )
      id=2
      DO k=kbgn,kend
        DO i=ibgn,iend
          ik = i-ibgn+1 + (k-kbgn)*isize
          tem1(ik) = -9999.0
          tem2(ik) = -9999.0
          IF(u(i,jslice,k) /= -9999.0) tem1(ik)=u(i,jslice,k)*factor
          IF(w(i,jslice,k) /= -9999.0) tem2(ik)=w(i,jslice,k)*factor
          tem3(ik)=x(i,jslice,k)
          tem4(ik)=z(i,jslice,k)
        END DO
      END DO
      WRITE(title,'(''U-W '',A)')label
    END IF

    length = 120
    CALL strlnth( title, length )
    CALL strmin ( title, length)

    DO i=1,smooth
      CALL smooth9pmv(tem1,isize,ksize,1,isize,1,ksize,tem6)
      CALL smooth9pmv(tem2,isize,ksize,1,isize,1,ksize,tem6)
    END DO

    uunit = 10.0
    CALL xvmode(1)
    istep = ist
    jstep = kst
    CALL vtr2d(tem1,tem2,tem3,tem4,uunit, xw,xe,dx,zb,zt,dz,            &
               isize,istep,ksize,jstep,title(1:length),runname, id,     &
               tem5,slicopt,mnsize,tem6(tind1),tem6(tind2),tem6(tind3), &
               tem6(tind4),tem6(tind5),tem6(tind6),tem6(tind7),         &
               tem6(tind8))

!    IF(varname(1:6) /= 'xuvplt') CALL set_vertical_factor( 1.0 )
!
!-----------------------------------------------------------------------
!
!  slicopt=3   Plot v-w field
!
!-----------------------------------------------------------------------
!
  ELSE IF( slicopt == 3 .OR. slicopt == 0 ) THEN

!    x_tmp = y(1,jslice,1)
    x_tmp = x(islice,1,1)

    i = islice

    i = i+ (xpbgn-1)*(nx-3)
    dist = (i-1.5)*tmpx
    length=LEN_TRIM(distc)
    CALL strmin ( distc, length)
    WRITE(levlab,'(''Y-Z Plane at X='',F8.1,A)')dist,distc(1:length)

    IF(varname(1:6) == 'xuvplt') THEN
      xw1=xw
      xe1=xe
      ys1=ys
      yn1=yn
      id=4
      DO k=kbgn,kend
        DO j=jbgn,jend
          jk = j-jbgn+1 + (k-kbgn)*jsize
          tem1(jk) = -9999.0
          tem2(jk) = -9999.0
          IF(u(islice,j,k) /= -9999.0) tem1(jk)=u(islice,j,k)*factor
          IF(v(islice,j,k) /= -9999.0) tem2(jk)=v(islice,j,k)*factor
          tem3(jk)=y(islice,j,k)
          tem4(jk)=z(islice,j,k)
        END DO
      END DO
      WRITE(title,'(''U-V '',A)')label
    ELSE
!      CALL set_vertical_factor ( (zt-zb)/(yn-ys) )
      id=3
      DO k=kbgn,kend
        DO j=jbgn,jend
          jk = j-jbgn+1 + (k-kbgn)*jsize
          tem1(jk) = -9999.0
          tem2(jk) = -9999.0
          IF(v(islice,j,k) /= -9999.0) tem1(jk)=v(islice,j,k)*factor
          IF(w(islice,j,k) /= -9999.0) tem2(jk)=w(islice,j,k)*factor
          tem3(jk)=y(islice,j,k)
          tem4(jk)=z(islice,j,k)
        END DO
      END DO
      WRITE(title,'(''V-W '',A)')label
    END IF

    length = 120
    CALL strlnth( title, length )
    CALL strmin ( title, length)

    DO i=1,smooth
      CALL smooth9pmv(tem1,jsize,ksize,1,jsize,1,ksize,tem6)
      CALL smooth9pmv(tem2,jsize,ksize,1,jsize,1,ksize,tem6)
    END DO

    uunit = 10.0
    CALL xvmode(1)
    istep = jst
    jstep = kst
    CALL vtr2d(tem1,tem2,tem3,tem4,uunit, ys,yn,dy,zb,zt,dz,            &
               jsize,istep,ksize,jstep,title(1:length),runname, id,     &
               tem5,slicopt,mnsize,tem6(tind1),tem6(tind2),tem6(tind3), &
               tem6(tind4),tem6(tind5),tem6(tind6),tem6(tind7),         &
               tem6(tind8))

!    IF(varname(1:6) /= 'xuvplt') CALL set_vertical_factor( 1.0 )
!
!-----------------------------------------------------------------------
!
!  slicopt=4   Plot u-v field on constant z levels
!  slicopt=6   Plot u-v field on constant pressure levels
!  slicopt=7   Plot u-v field on constant PT levels
!
!-----------------------------------------------------------------------
!
  ELSE IF( slicopt == 4 .OR. slicopt == 6 .OR. slicopt == 7 ) THEN

!    CALL hintrp(nx,ny,nz,u,z,zlevel,u1)
!    CALL hintrp(nx,ny,nz,v,z,zlevel,v1)
    IF( agl == 1 ) THEN
      print*,'AGL: calling hintrp2'
      CALL hintrp2(nx,ny,nz,kbgn,kend,u,z,zlevel,hterain,u1)
      CALL hintrp2(nx,ny,nz,kbgn,kend,v,z,zlevel,hterain,v1)
    ELSE
      CALL hintrp1(nx,ny,nz,kbgn,kend,u,z,zlevel,u1)
      CALL hintrp1(nx,ny,nz,kbgn,kend,v,z,zlevel,v1)
    END IF

    DO j=jbgn,jend
      DO i=ibgn,iend
        ij = i-ibgn+1 + (j-jbgn)*isize
        tem1(ij) = -9999.0
        tem2(ij) = -9999.0
        IF(u1(i,j) /= -9999.0) tem1(ij)=u1(i,j)*factor
        IF(v1(i,j) /= -9999.0) tem2(ij)=v1(i,j)*factor
        tem3(ij)=x(i,j,2)
        tem4(ij)=y(i,j,2)
      END DO
    END DO

    IF( slicopt == 4) THEN
      IF ( agl == 1 ) THEN
        IF( zlevel > 1.0) THEN
          WRITE(levlab,'(''Z='',F7.3,'' km AGL'')') zlevel
        ELSE
          WRITE(levlab,'(''Z='',F7.0,'' m AGL'')') (1000.*zlevel)
        END IF
      ELSE
        IF( zlevel > 1.0) THEN
          WRITE(levlab,'(''Z='',F7.3,'' km MSL'')') zlevel
        ELSE
          WRITE(levlab,'(''Z='',F7.0,'' m MSL'')') (1000.*zlevel)
        END IF
      END IF
    ELSE IF( slicopt == 6) THEN
      pp01 = 0.01*ercpl**zlevel
      WRITE(levlab,'(''P='',F7.2,A)') pp01, ' MB'
    ELSE
      WRITE(levlab,'(''THETA='',F5.1,A)') zlevel, ' (K)'
    END IF

    WRITE(title,'(''U-V '',A)') label

    length = 120
    CALL strlnth( title, length )
    CALL strmin ( title, length)

    uunit = 10.0
    CALL xvmode(1)
    istep = ist
    jstep = jst

    DO i=1,smooth
      CALL smooth9pmv(tem1,isize,jsize,1,isize,1,jsize,tem6)
      CALL smooth9pmv(tem2,isize,jsize,1,isize,1,jsize,tem6)
    END DO

    CALL vtr2d(tem1,tem2,tem3,tem4, uunit, xw,xe,dx,ys,yn,dy,           &
               isize,istep,jsize,jstep,title(1:length),runname, 1,      &
               tem5,slicopt,mnsize,tem6(tind1),tem6(tind2),tem6(tind3), &
               tem6(tind4),tem6(tind5),tem6(tind6),tem6(tind7),         &
               tem6(tind8))
!
!-----------------------------------------------------------------------
!  slicopt=12 Plot u-v field in vertical cross section along a trajectory
!-----------------------------------------------------------------------
!
  ELSE IF( slicopt == 12 .and. trajc_plt_opt == 2 ) THEN

    k = 1  ! supporting only ntimes=1 for slice plotting
    i = itrajc_index ! i=itrajc_start, itrajc_end, ntrajc_stride

    j_bgn = npoints_bgn(i,k)
    j_end = npoints_end(i,k)

!   print*,'trajectory No. ',i,', j_bgn,j_end=',j_bgn,j_end

    ALLOCATE(xptr(j_bgn:j_end,kbgn:kend),STAT=istatus)
    ALLOCATE(xptr1d(j_bgn:j_end),STAT=istatus)
    ALLOCATE(yptr1d(j_bgn:j_end),STAT=istatus)

    ALLOCATE(temptr1(j_bgn:j_end,kbgn:kend),STAT=istatus)
    ALLOCATE(temptr2(j_bgn:j_end,kbgn:kend),STAT=istatus)
    ALLOCATE(temptr3(j_bgn:j_end,kbgn:kend),STAT=istatus)
    ALLOCATE(temptr4(j_bgn:j_end,kbgn:kend),STAT=istatus)
    ALLOCATE(temptr5(j_bgn:j_end,kbgn:kend),STAT=istatus)
    ALLOCATE(temptr6(j_bgn:j_end,kbgn:kend),STAT=istatus)
    ALLOCATE(temptr7(j_bgn:j_end,kbgn:kend),STAT=istatus)
    ALLOCATE(temptr8(j_bgn:j_end,kbgn:kend),STAT=istatus)

    ALLOCATE(zptr(j_bgn:j_end,nz),STAT=istatus)
    ALLOCATE(u2_tr(j_bgn:j_end,1:nz),STAT=istatus)
    ALLOCATE(w2_tr(j_bgn:j_end,1:nz),STAT=istatus)

    xptr1d(j_bgn:j_end)=xtrajc(j_bgn:j_end,i,k)*0.001
    yptr1d(j_bgn:j_end)=ytrajc(j_bgn:j_end,i,k)*0.001

    DO kz=kbgn,kend
      xptr(j_bgn,kz)=0.0
      DO j=j_bgn+1,j_end
        xptr(j,kz) = xptr(j-1,kz)+sqrt((xptr1d(j)-xptr1d(j-1))**2+ &
                                       (yptr1d(j)-yptr1d(j-1))**2)
      END DO
    END DO

    xptr1 = xptr(j_bgn,kbgn)
    xptr2 = xptr(j_end,kbgn)
    dxptr=(xptr2-xptr1)/(j_end-j_bgn)

    CALL sectvrt(nx,ny,nz,u,x,y,z,dx,dy,u2_tr,zptr,j_end-j_bgn+1,xptr1d,yptr1d)
    CALL sectvrt(nx,ny,nz,v,x,y,z,dx,dy,w2_tr,zptr,j_end-j_bgn+1,xptr1d,yptr1d)

! Project (u,v) to the local tangent of the trajectory plane.

    id=2
    DO kz=kbgn,kend
      DO j=j_bgn,j_end

        IF(u2_tr(j,kz) /= -9999.0 .AND. w2_tr(j,kz) /= -9999.0) THEN

          delta_x = xptr1d(min(j+1,j_end))-xptr1d(max(j_bgn,j-1))
          delta_y = yptr1d(min(j+1,j_end))-yptr1d(max(j_bgn,j-1))
          dist = sqrt(delta_x**2+delta_y**2)

          sinaf_ptr=delta_y/(dist+0.000001)
          cosaf_ptr=delta_x/(dist+0.000001)

          u2_tr(j,kz)=(u2_tr(j,kz)*cosaf_ptr+w2_tr(j,kz)*sinaf_ptr)*factor
        ELSE
          u2_tr(j,kz)=-9999.0
        ENDIF

      END DO
    END DO

    CALL sectvrt(nx,ny,nz,w,x,y,z,dx,dy,w2_tr,zptr,j_end-j_bgn+1,xptr1d,yptr1d)

    DO kz=kbgn,kend
      DO j=j_bgn,j_end
        IF(w2_tr(j,kz) /= -9999.0) w2_tr(j,kz)=w2_tr(j,kz)*factor
      END DO
    END DO

    length=LEN_TRIM(distc)
    CALL strmin ( distc, length)

    IF(axlbfmt == -1 .OR. axlbfmt == 1 ) THEN
      length=LEN_TRIM(distc)
      CALL strmin(distc,length)
      IF(varname(1:6) == 'xuvplt') THEN
        length=LEN_TRIM(distc)
        CALL strmin(distc,length)
        WRITE(title,'(''U-V '',A)') label
        WRITE(levlab,'(''XY-Z PLOT FROM '',4(A,F5.1),A,A)')             &
            '(',x101,',',y101,') through (',x102,',',y102,') ',         &
            distc(1:length)
      ELSE
        WRITE(title,'(''UV-W '',A)') label
        WRITE(levlab,                                                   &
            '(''Vertical Plane From '',4(A,F8.1),A,A)')                 &
            '(',x101,',',y101,') through (',x102,',',y102,') ',         &
            distc(1:length)
      END IF

    ELSE IF(axlbfmt == 0) THEN
      length=LEN_TRIM(distc)
      CALL strmin ( distc, length)
      IF(varname(1:6) == 'xuvplt') THEN
        WRITE(title,'(''U-V '',A)') label
        WRITE(levlab,'(''XY-Z PLOT FROM '',4(A,I5),A,A)')               &
            '(',NINT(x101),',',NINT(y101),') through (',                &
            NINT(x102),',',NINT(y102),') ',distc(1:length)
      ELSE
        WRITE(title,'(''UV-W '',A)') label
        WRITE(levlab,                                                   &
            '(''Vertical Plane From '',4(A,I5),A,A)')                   &
            '(',NINT(x101),',',NINT(y101),') through (',                &
            NINT(x102),',',NINT(y102),') ',distc(1:length)
      END IF
    ELSE
      length=LEN_TRIM(distc)
      CALL strmin ( distc, length)
!     WRITE(stem1,'(i1)')axlbfmt
!     WRITE(stem2,'(a3,a1)')'f8.',stem1

      IF(varname(1:6) == 'xuvplt') THEN
        WRITE(title,'(''U-V '',A)') label
        WRITE(levlab,'(''XY-Z PLOT FROM '',4(A,f8.2),A,A)')            &
            '(',x101,',',y101,') through (',x102,',',y102,') ',         &
            distc(1:length)
      ELSE
        WRITE(title,'(''UV-W '',A)') label
        WRITE(levlab,                                                   &
            '(''Vertical Plane From '',4(A,f8.2),A,A)')                &
            '(',x101,',',y101,') through (',x102,',',y102,') ',         &
            distc(1:length)
      END IF
    END IF

    length = 120
    CALL strlnth( title, length )
    CALL strmin ( title, length)

    DO kz=1,smooth
      CALL smooth9pmv(u2_tr(j_bgn:j_end,kbgn:kend),j_end-j_bgn+1, &
           kend-kbgn+1,1,j_end-j_bgn+1,1,kend-kbgn+1, temptr1)
      CALL smooth9pmv(w2_tr(j_bgn:j_end,kbgn:kend),j_end-j_bgn+1, &
             kend-kbgn+1,1,j_end-j_bgn+1,1,kend-kbgn+1, temptr1)
    END DO

    uunit = 10.0
    CALL xvmode(1)
    istep = ist
    jstep = kst
    CALL vtr2d(u2_tr(j_bgn:j_end,kbgn:kend),w2_tr(j_bgn:j_end,kbgn:kend), &
               xptr(j_bgn:j_end,kbgn:kend), zptr(j_bgn:j_end,kbgn:kend), &
               uunit, xptr1,xptr2,dxptr,zb,zt,dz,j_end-j_bgn+1, &
               istep,kend-kbgn+1,jstep,title(1:length),runname, id,  &
               tem5,slicopt,(j_end-j_bgn+1)*(kend-kbgn+1), &
               temptr1,temptr2,temptr3,temptr4,temptr5,temptr6,temptr7,temptr8)

    DEALLOCATE(xptr)
    DEALLOCATE(xptr1d)
    DEALLOCATE(yptr1d)

    DEALLOCATE(temptr1)
    DEALLOCATE(temptr2)
    DEALLOCATE(temptr3)
    DEALLOCATE(temptr4)
    DEALLOCATE(temptr5)
    DEALLOCATE(temptr6)
    DEALLOCATE(temptr7)
    DEALLOCATE(temptr8)

    DEALLOCATE(zptr)
    DEALLOCATE(u2_tr)
    DEALLOCATE(w2_tr)

!-----------------------------------------------------------------------
!  slicopt=5 Plot u-v field in vertical cross section through 2 points
!-----------------------------------------------------------------------

  ELSE IF( slicopt == 5 ) THEN

    CALL sectvrt(nx,ny,nz,u,x,y,z,dx,dy,u2,zp,n,xp,yp)
    CALL sectvrt(nx,ny,nz,v,x,y,z,dx,dy,v2,zp,n,xp,yp)
    CALL sectvrt(nx,ny,nz,w,x,y,z,dx,dy,w2,zp,n,xp,yp)

    IF(varname(1:6) == 'xuvplt') THEN
      xw1=xw
      xe1=xe
      ys1=ys
      yn1=yn
      id=4
      DO k=kbgn,kend
        DO i=ibgn,iend
          ik = i-ibgn+1 + (k-kbgn)*isize
          tem1(ik) = -9999.0
          tem2(ik) = -9999.0
          IF(u2(i,k) /= -9999.0) tem1(ik)= u2(i,k)*factor
          IF(v2(i,k) /= -9999.0) tem2(ik)= v2(i,k)*factor
          tem3(ik)=xw+(i-ibgn)* sqrtdxy
          tem4(ik)=zp(i,k)
        END DO
      END DO
    ELSE
      id=2
      DO k=kbgn,kend
        DO i=ibgn,iend
          ik = i-ibgn+1 + (k-kbgn)*isize
          tem1(ik) = -9999.0
          tem2(ik) = -9999.0
          IF(u2(i,k) /= -9999.0 .AND. v2(i,k) /= -9999.0)               &
               tem1(ik)=(u2(i,k)*cosaf+v2(i,k)*sinaf)*factor
          IF(w2(i,k) /= -9999.0) tem2(ik)=w2(i,k)*factor
          tem3(ik)=xw+(i-ibgn)* sqrtdxy
          tem4(ik)=zp(i,k)
        END DO
      END DO
    END IF

    IF(axlbfmt == -1 .OR. axlbfmt == 1 ) THEN
      length=LEN_TRIM(distc)
      CALL strmin(distc,length)
      IF(varname(1:6) == 'xuvplt') THEN
        length=LEN_TRIM(distc)
        CALL strmin(distc,length)
        WRITE(title,'(''U-V '',A)') label
        WRITE(levlab,'(''XY-Z PLOT FROM '',4(A,F5.1),A,A)')             &
            '(',x101,',',y101,') through (',x102,',',y102,') ',         &
            distc(1:length)
      ELSE
        WRITE(title,'(''UV-W '',A)') label
        WRITE(levlab,                                                   &
            '(''Vertical Plane From '',4(A,F8.1),A,A)')                 &
            '(',x101,',',y101,') through (',x102,',',y102,') ',         &
            distc(1:length)
      END IF
    ELSE IF(axlbfmt == 0) THEN
      length=LEN_TRIM(distc)
      CALL strmin ( distc, length)
      IF(varname(1:6) == 'xuvplt') THEN
        WRITE(title,'(''U-V '',A)') label
        WRITE(levlab,'(''XY-Z PLOT FROM '',4(A,I5),A,A)')               &
            '(',NINT(x101),',',NINT(y101),') through (',                &
            NINT(x102),',',NINT(y102),') ',distc(1:length)
      ELSE
        WRITE(title,'(''UV-W '',A)') label
        WRITE(levlab,                                                   &
            '(''Vertical Plane From '',4(A,I5),A,A)')                   &
            '(',NINT(x101),',',NINT(y101),') through (',                &
            NINT(x102),',',NINT(y102),') ',distc(1:length)
      END IF
    ELSE
      length=LEN_TRIM(distc)
      CALL strmin ( distc, length)
!     WRITE(stem1,'(i1)')axlbfmt
!     WRITE(stem2,'(a3,a1)')'f8.',stem1

      IF(varname(1:6) == 'xuvplt') THEN
        WRITE(title,'(''U-V '',A)') label
        WRITE(levlab,'(''XY-Z PLOT FROM '',4(A,f8.2),A,A)')            &
            '(',x101,',',y101,') through (',x102,',',y102,') ',         &
            distc(1:length)
      ELSE
        WRITE(title,'(''UV-W '',A)') label
        WRITE(levlab,                                                   &
            '(''Vertical Plane From '',4(A,f8.2),A,A)')                &
            '(',x101,',',y101,') through (',x102,',',y102,') ',         &
            distc(1:length)
      END IF
    END IF

    length = 120
    CALL strlnth( title, length )
    CALL strmin ( title, length)

    DO i=1,smooth
      CALL smooth9pmv(tem1,isize,ksize,1,isize,1,ksize,tem6)
      CALL smooth9pmv(tem2,isize,ksize,1,isize,1,ksize,tem6)
    END DO

    uunit = 10.0
    CALL xvmode(1)
    istep = ist
    jstep = kst
    CALL vtr2d(tem1,tem2,tem3,tem4,uunit, xw,xe,sqrtdxy,zb,zt,dz,       &
               isize,istep,ksize,jstep,title(1:length),runname, id,     &
               tem5,slicopt,mnsize,tem6(tind1),tem6(tind2),tem6(tind3), &
               tem6(tind4),tem6(tind5),tem6(tind6),tem6(tind7),         &
               tem6(tind8))

  END IF

  IF(myproc == 0) CALL xpscmnt('End plotting '//label_copy(1:llabel))

  RETURN
END SUBROUTINE vtr3d
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE VTR2D                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE vtr2d(u,v,x,y,uunit1, xl,xr,dx,yb,yt,dy,                     &
                 m,istep,n,jstep,char1,char2, vpltmod,                  &
                 hterain,slicopt,mnsize,                                &
                 plotu,plotv,plota,plotx,ploty,iwrk,xwk,ywk)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    Plot 2-d wind (u1,u2) vector field defined on grid points (x,y)
!    using ZXPLOT package..
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!
!  MODIFICATION HISTORY:
!
!  1/24/96 (J. Zong and M. Xue)
!  Fixed a problem related to finding the minima and maxima of u & v
!  when there exist missing data. The initial min. and max. should be
!  set to values other than the missing value, -9999.0.
!
!-----------------------------------------------------------------------
!
!
!  INPUT:
!    u        m by n 2-dimensional array of u (left-to-right)
!               wind components (m/s)
!    v        m by n 2-dimensional array of v (down-to-up)
!               wind components (m/s)
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!
!    uunit1
!
!    xl,xr    The left and right bound of the physical domain.
!    dx       Spacing between the x-axis tick marks
!    yb,yt    Bottom and top bound of the physical domain.
!    dy       Spacing between the y-axis tick marks
!
!    m        First dimension of vector component array
!    istep    Step increment for plotting in x direction
!
!    n        Second dimension of vector component array
!    jstep    Step increment for plotting in y direction
!
!    char1    First character string to plot (title)
!    char2    Second character string to plot (runname)
!
!    vpltmod  vpltmod = 1 for u-v vector   (u=u, v=v in model space)
!             vpltmod = 2 for u-w vector   (u=u, v=w in model space)
!             vpltmod = 3 for v-w vector   (u=v, v=w in model space)
!    hterain  the height of terrain
!    slicopt  slice orientation indicator
!             = 1, x-y slice of at k=kslice is plotted.
!             = 2, x-z slice of at j=jslice is plotted.
!             = 3, y-z slice of at i=islice is plotted.
!             = 4, horizontal slice at z index islice is plotted.
!             = 5, xy-z cross section of wind islice is plotted.
!             = 6, data field on constant p-level is plotted.
!             = 0, all of the three slices above are plotted.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INCLUDE 'arpsplt.inc'

  INTEGER, INTENT(IN) :: m,n

  REAL,    INTENT(IN) :: u(m,n)
  REAL,    INTENT(IN) :: v(m,n)
  REAL,    INTENT(IN) :: x(m,n)
  REAL,    INTENT(IN) :: y(m,n)

  REAL,    INTENT(IN) :: uunit1
  REAL,    INTENT(IN) :: xl,xr,dx,yb,yt,dy
  INTEGER, INTENT(IN) :: istep,jstep

  CHARACTER(LEN=*), INTENT(IN)    :: char2
  CHARACTER(LEN=*), INTENT(INOUT) :: char1

  INTEGER, INTENT(IN) :: slicopt,vpltmod

  REAL,    INTENT(IN) :: hterain(m,n)             ! The height of the terrain.

  INTEGER, INTENT(IN)    :: mnsize
  REAL,    INTENT(INOUT) :: plotu(mnsize)
  REAL,    INTENT(INOUT) :: plotv(mnsize)
  REAL,    INTENT(INOUT) :: plota(mnsize)
  REAL,    INTENT(INOUT) :: plotx(mnsize)
  REAL,    INTENT(INOUT) :: ploty(mnsize)
  INTEGER, INTENT(INOUT) :: iwrk(mnsize)
  REAL,    INTENT(INOUT) :: xwk(mnsize), ywk(mnsize)

!-----------------------------------------------------------------------
!
!  Plotting control common blocks
!
!-----------------------------------------------------------------------
!
  INTEGER :: layover
  REAL    :: ctinc,ctmin,ctmax,vtunt    !contour interval and vector unit
  REAL    :: xleng,vunit
  REAL    :: yxratio                  !the scaling factor the y/x ratio.
  INTEGER :: iunits, itype

  COMMON /laypar/  layover
  COMMON /incunt/  ctinc,ctmin,ctmax,vtunt
  COMMON /vecscl/  xleng,vunit
  COMMON /yratio/  yxratio
  COMMON /windvtr/ iunits, itype

  INTEGER :: ovrstaopt
  INTEGER :: ovrstam,staset,ovrstan,ovrstav,stacol,markprio,wrtstax
  INTEGER :: nsta_typ,sta_typ(30),sta_marktyp(30),sta_markcol(30)
  REAL    :: sta_marksz(30),wrtstad
  CHARACTER (LEN=256) :: stalofl
  COMMON /sta_par/ ovrstaopt,ovrstam,staset,ovrstan,ovrstav,stacol,     &
         markprio,nsta_typ,sta_typ,sta_marktyp,                         &
         sta_markcol,sta_marksz,stalofl,wrtstax,wrtstad

  INTEGER :: icolor,icolor1,lbcolor,trcolor                ! required color
  COMMON /recolor/icolor,icolor1,lbcolor,trcolor

  INTEGER :: ovrobs,obsset,obscol,obs_marktyp
  REAL :: obs_marksz,obs_valsz
  COMMON /obspar/ ovrobs,obsset,obscol,obs_marktyp,                    &
                  obs_marksz,obs_valsz

  REAL :: lblmag, ctrlbsiz, axlbsiz
  COMMON /labmag/ lblmag, ctrlbsiz, axlbsiz

  INTEGER :: flag, haxisu, vaxisu, lbaxis, tickopt, axlbfmt
  INTEGER :: xfont   ! the font of character
  REAL    :: hmintick,vmajtick,vmintick,hmajtick
  COMMON /var_par/ xfont,haxisu,vaxisu,lbaxis,tickopt,hmintick,         &
         vmajtick, vmintick,hmajtick,axlbfmt

  REAL :: ubarb(200,200), vbarb(200,200)
  COMMON /windtmp/ubarb, vbarb

  REAL :: zlevel
  COMMON /sliceh/ zlevel

  INTEGER :: timeovr
  COMMON /timover/ timeovr

  INTEGER :: ntitle,titcol, nxpic, nypic, wpltime
  REAL    :: titsiz
  CHARACTER (LEN=256) :: ptitle(3), footer_l, footer_c, footer_r
  COMMON /titpar1/ptitle, footer_l, footer_c, footer_r
  COMMON /titpar2/ntitle,titcol,wpltime, nxpic, nypic
  COMMON /titpar3/titsiz

  INTEGER :: col_table,pcolbar
  COMMON /coltable/ col_table,pcolbar

  CHARACTER (LEN=12) :: varname
  COMMON /varplt1/ varname

  REAL :: xw1,xe1,ys1,yn1
  COMMON /xuvpar/ xw1,xe1,ys1,yn1

  INTEGER :: wrtflag
  CHARACTER (LEN=90) :: levlab
  CHARACTER (LEN=50) :: timelab
  CHARACTER (LEN=25) :: timestring
  COMMON /timelev/wrtflag, timelab, levlab, timestring

  CHARACTER (LEN=80) :: prestr
  INTEGER            :: preflag
  COMMON /preinfo/ prestr,preflag

  REAL :: x101, y101, x102,y102
  COMMON /slicev1/x101, y101, x102,y102

  INTEGER :: xnwpic_called
  COMMON /callnwpic/xnwpic_called

  INTEGER :: xpbgn,xpend,ypbgn,ypend        ! for MPI jobs
  COMMON /processors/ xpbgn,xpend,ypbgn,ypend
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,key
  REAL :: pl,pr,pb,pt   ! plot space left, right, bottom, top coordinate
  REAL :: px,py         ! plot space left-right length and up-down height
  REAL :: xs,ys         ! real space left-right length and up-down height
  REAL :: pxc,pyc       ! plot space left-right center and
                        ! up-down    center
  REAL :: x0,y0
  REAL :: umax,umin     ! max and min of u component
  REAL :: vmax,vmin     ! max and min of v component
  REAL :: uunit, uunit0
  REAL :: am


  INTEGER :: len0, len1
  REAL    :: xleng0,istand
  INTEGER :: iunits0

  CHARACTER (LEN=15)  :: ichar2
  CHARACTER (LEN=150) :: f_char1
  CHARACTER (LEN=150) :: ch

  REAL :: ytmp   !!local temporary variable

  REAL :: f_cputime,cpu1,cpu2
  DOUBLE PRECISION :: f_walltime,second1,second2

  REAL, ALLOCATABLE :: var_tr(:,:)
  REAL, ALLOCATABLE :: xptr(:,:)
  REAL, ALLOCATABLE :: zptr(:,:)
  REAL, ALLOCATABLE :: temptr1(:,:),temptr2(:,:),temptr3(:,:)
  REAL, ALLOCATABLE :: temptr4(:,:),temptr5(:,:),temptr6(:,:)
  REAL, ALLOCATABLE :: xptr1d(:),yptr1d(:)
  REAL :: xptr1,xptr2,dxptr,xtr,ytr
  INTEGER :: itrajc_start, itrajc_end, iptr
  INTEGER :: j_bgn,j_end
  INTEGER :: istatus,kz


  INCLUDE 'mp.inc'
  INCLUDE 'arpstrajc.inc'

  INTEGER :: ii,jj,mm,nn
  INTEGER :: ierr, itags, itagr

  INTEGER, PARAMETER :: destination = 0
  INTEGER            :: source
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  second1= f_walltime()
  cpu1 = f_cputime()

  IF(myproc == 0) THEN
    WRITE(6,'(/1x,a,a)') 'Plotting ',char1

    IF( layover == 0 .OR. xnwpic_called == 0) THEN
      CALL xnwpic
      xnwpic_called=1
      timeovr=0
      wrtflag = 0
      preflag = 0
      prestr = levlab
      len1=LEN_TRIM(prestr)
      CALL strmin(prestr,len1)
    ELSE
      timeovr=1
      wrtflag = wrtflag + 1
    END IF
!
!-----------------------------------------------------------------------
!
!  Get plotting space variables
!
!-----------------------------------------------------------------------
!
    CALL xqpspc( pl, pr, pb, pt)
    px = pr - pl
    py = pt - pb

    xs = xr-xl
    ys = yt-yb

    pxc = (pr+pl)/2
    pyc = (pb+pt)/2
!
!-----------------------------------------------------------------------
!
!  Let the longest lenth determine size scaling of plot
!
!-----------------------------------------------------------------------
!
    IF( py/px >= ys*yxratio/xs ) THEN
      py = ys*yxratio/xs*px
      CALL xpspac(pl, pr, pyc-py/2, pyc+py/2 )
    ELSE
      px = xs/(ys*yxratio)*py
      CALL xpspac(pxc-px/2, pxc+px/2, pb, pt)
    END IF
!
!-----------------------------------------------------------------------
!
!  Set the real distance to plot distance scaling
!
!-----------------------------------------------------------------------
!
    CALL xmap( xl, xr, yb,yt)
!
!-----------------------------------------------------------------------
!
!  Plot maps, boxes, and polygons
!
!-----------------------------------------------------------------------
!
    CALL xcolor(lbcolor)

    CALL pltextra(slicopt, 1 )

  END IF  ! myproc == 0
!
!-----------------------------------------------------------------------
!
!  Find max and min of data array
!
!-----------------------------------------------------------------------
!
  DO j=1,n
    DO i=1,m
      IF(u(i,j) == -9999.0 .OR. v(i,j) == -9999.0) CYCLE
      umin = u(i,j)
      vmin = v(i,j)
      GO TO 110
    END DO
  END DO
  110   CONTINUE

  umax=umin
  vmax=vmin

  DO j=1,n
    DO i=1,m
      IF(u(i,j) > umax .AND. u(i,j) /= -9999.0) umax=u(i,j)
      IF(u(i,j) < umin .AND. u(i,j) /= -9999.0) umin=u(i,j)
      IF(v(i,j) > vmax .AND. v(i,j) /= -9999.0) vmax=v(i,j)
      IF(v(i,j) < vmin .AND. v(i,j) /= -9999.0) vmin=v(i,j)
    END DO
  END DO

  IF (loc_x < xpbgn .OR. loc_x > xpend .OR. loc_y < ypbgn .OR. loc_y > ypend) THEN
    umax = -9999.9
    vmax = -9999.9
    umin =  9999.9
    vmin =  9999.9
  END IF
  CALL mpmax0(umax, umin) !? only inside xpbgn-xpend, ypbgn-ypend
  CALL mpmax0(vmax, vmin) !? only inside xpbgn-xpend, ypbgn-ypend
!
!-----------------------------------------------------------------------
!
!  Fill various labels
!
!-----------------------------------------------------------------------
!
  IF(myproc == 0) THEN

    CALL xchmag( 0.025*px  * lblmag )
    CALL xcolor(lbcolor)

    IF ( layover < 1 ) THEN

      len1=LEN_TRIM(timelab)
      CALL strmin(timelab,len1)
      CALL xcharc((xl+xr)*0.5,yt+0.015*ys*px/py,                            &
                        timestring(1:25)//' '//timelab(1:len1))

      IF(levlab /= ' ') THEN
        len1=LEN_TRIM(levlab)
        CALL strmin(levlab,len1)
        CALL xcharc(xl+xs*0.5,yt+0.06*ys*px/py, levlab(1:len1))
        preflag = 1
      END IF
!     len1=80
!     CALL strmin(levlab,len1)
!     CALL xcharc(xl+xs*0.5,yt+0.06*ys*px/py, levlab(1:len1))
    END IF
    IF(preflag == 0 .AND. levlab /= ' ') THEN
      len1=LEN_TRIM(levlab)
      CALL strmin(levlab,len1)
      CALL xcharc(xl+xs*0.5,yt+0.06*ys*px/py, levlab(1:len1))
      preflag = 1
    END IF

    IF( vpltmod == 1 .OR. vpltmod == 4 ) THEN
      WRITE(ch,'(4(a,F7.2))')                                           &
          'Umin=',umin,' Umax=',umax,' Vmin=',vmin,' Vmax=',vmax
    ELSE IF( vpltmod == 2 ) THEN
      WRITE(ch,'(4(a,F7.2))')                                           &
          'Umin=',umin,' Umax=',umax,' Wmin=',vmin,' Wmax=',vmax
    ELSE
      WRITE(ch,'(4(a,F7.2))')                                           &
          'Vmin=',umin,' Vmax=',umax,' Wmin=',vmin,' Wmax=',vmax
    END IF

    LEN0= LEN_TRIM(char1)
    CALL strmin(char1,LEN0)
    IF( char1(LEN0:LEN0) == ')' ) char1(LEN0:LEN0)=','
    IF(itype == 1) THEN
      WRITE(f_char1, '(a,'' Vector)'')') char1(1:LEN0)
    ELSE IF(itype == 2) THEN
      WRITE(f_char1, '(a, '' Barb)'')')char1(1:LEN0)
    END IF

! if first levlab is not equal second levlab then attatch levlab on f_ch
!
    LEN0=LEN_TRIM(f_char1)
    CALL strmin(f_char1,LEN0)
    len1=LEN_TRIM(levlab)
    CALL strmin(levlab,len1)
!  IF (preflag.eq.1 .and. prestr(1:len1).ne.levlab(1:len1)
!    :   .and. prestr(1:1).ne.' '
!    :   .and.layover.ne.0 .and. levlab(1:1).ne.' ') THEN
!    write(f_char1,'(a,a)') f_char1(1:len0),levlab(1:len1)
!  ENDIF

    WRITE(6,'(1x,a51)') ch(1:51)
    CALL xcolor(icolor)

    IF(lbaxis == 1) THEN
      IF(wrtstax == 0) THEN
        ytmp = 0.08
      ELSE
        ytmp =0.14
      END IF
    ELSE
      ytmp = 0.12
    END IF
    LEN0=LEN_TRIM(f_char1)
    CALL strmin(f_char1,LEN0)

    CALL xchmag(0.025*px * lblmag )
    CALL xcharl(xl-0.05*(xr-xl), yb-(yt-yb)*px/py*(ytmp+wrtflag*0.030),         &
             f_char1(1:LEN0))

    len1=LEN_TRIM(ch)
    CALL strmin(ch,len1)
    CALL xcharr(xr+0.05*(xr-xl), yb-(yt-yb)*px/py*(ytmp+wrtflag*0.030),         &
             ch(1:len1))
!
!-----------------------------------------------------------------------
!
!  Set vector unit and plot vectors.
!
!-----------------------------------------------------------------------
!
  ! Set parameter for barb

    xleng0 = (pr-pl)/(m-1) * istep * 0.65
    IF(iunits == 1 .AND. itype == 2) THEN
      iunits0=1
      istand = 5.
      WRITE(ichar2,'(a15)')'5 m/s'
    ELSE IF(iunits == 2 .AND. itype == 2) THEN
      iunits0=2
      istand = 10.
      WRITE(ichar2,'(a15)')'10 knots'
    ELSE IF (iunits == 3 .AND. itype == 2) THEN
      iunits0=2
      istand = 10.
      WRITE(ichar2,'(a15)')'10 MPH'
    END IF

    IF(layover >= 1) CALL xcolor(icolor)

    CALL xcolor(icolor)
    CALL xwindw(xl, xr, yb, yt)

  END IF  ! myproc == 0

  uunit=uunit1
  IF( vtunt /= 0.0 ) THEN
    uunit=vtunt
    CALL xvmode(2)
  END IF
  CALL xmap(xl,xr, yb,yt)

  IF (mp_opt > 0) THEN
    IF (myproc == 0) THEN
      CALL xvectu_mp(xl,xr,yb,yt,umax,umin,vmax,vmin,m,istep,xleng,uunit)
!      xleng = (xleng-2.0)/(xpend-xpbgn+1)  ! m is smaller in mpi mode than serial mode
      xleng = (xleng)/(xpend-xpbgn+1)  ! m is smaller in mpi mode than serial mode
    END IF
  ELSE
    CALL xvectu(u,v,m,m,istep,n,jstep,xleng,uunit)
  END IF

  DO j = 1, n
    DO i = 1, m
      ii = i+(j-1)*m
      plotu(ii) = u(i,j)
      plotv(ii) = v(i,j)
      plotx(ii) = x(i,j)
      ploty(ii) = y(i,j)
    END DO
  END DO
  mm = m
  nn = n

  DO jj = ypbgn,ypend
    DO ii = xpbgn, xpend

      source = (ii+(jj-1)*nproc_x-1)
      IF (source == 0) GOTO 600

      CALL inctag
      IF (myproc == source ) THEN

        itags = gentag + 4
        CALL mpsendi(m,1,destination,itags,ierr)
        itags = gentag + 5
        CALL mpsendi(n,1,destination,itags,ierr)

        itags = gentag
        CALL mpsendr(u,m*n,destination,itags,ierr)
        itags = gentag+3
        CALL mpsendr(v,m*n,destination,itags,ierr)
        itags = gentag + 1
        CALL mpsendr(x,m*n,destination,itags,ierr)
        itags = gentag + 2
        CALL mpsendr(y,m*n,destination,itags,ierr)
      END IF


      IF (myproc == 0) THEN

        plotu = 0.0
        plotv = 0.0
        plotx = 0.0
        ploty = 0.0
        mm = 0
        nn = 0

        itagr = gentag + 4
        CALL mprecvi(mm,1,source,itagr,ierr)
        itagr = gentag + 5
        CALL mprecvi(nn,1,source,itagr,ierr)

        itagr = gentag
        CALL mprecvr(plotu,mm*nn,source,itagr,ierr)
        itagr = gentag + 3
        CALL mprecvr(plotv,mm*nn,source,itagr,ierr)
        itagr = gentag + 1
        CALL mprecvr(plotx,mm*nn,source,itagr,ierr)
        itagr = gentag + 2
        CALL mprecvr(ploty,mm*nn,source,itagr,ierr)

      END IF

      600 CONTINUE

      IF (myproc == 0) THEN
        IF(itype == 1) THEN
          CALL xvectr(plotu,plotv,plotx,ploty,mm,mm,istep,nn,jstep,xleng,uunit)
        ELSE IF(itype == 2) THEN
          CALL xbarbs(plotu,plotv,plotx,ploty,mm,mm,istep,nn,jstep,iunits0,xleng*0.65,2)
        END IF

        CALL xwdwof

      END IF  ! myproc == 0

      CALL mpbarrier         ! sync the processors

    END DO
  END DO

  IF (myproc == 0) THEN
!
!-----------------------------------------------------------------------
!
!  Plot axes with tick marks
!
!-----------------------------------------------------------------------
!
    CALL pltaxes(slicopt,dx,dy)

    vunit=uunit
    x0=xl-(xr-xl)*0.05
!   y0=yb+(yt-yb)*0.07
    y0=yt+0.020*ys*px/py

    key=0
    am=0.5
    IF( ((m-1)/istep) > 30 ) am=1.0
    IF(itype == 1) THEN
      IF(varname(1:6) == 'xuvplt') CALL xmap(xl, xr,yb,yt)
      CALL xvectk(x0,y0,xleng*am,uunit*am, key)

      CALL xmap(xl, xr, yb, yt)
    END IF

  END IF ! myproc == 0

!
!-----------------------------------------------------------------------
!
! Plot terrain etc.
!
!-----------------------------------------------------------------------
!

  DO j = 1,n
    DO i = 1,m
        ii = i+ (j-1)*m
        plota(ii) = hterain(i,j)
        plotu(ii) = u(i,j)
        plotv(ii) = v(i,j)
        plotx(ii) = x(i,j)
        ploty(ii) = y(i,j)
    END DO
  END DO
  mm = m
  nn = n

  DO jj = ypbgn,ypend
    DO ii = xpbgn, xpend

      source = (ii+(jj-1)*nproc_x-1)
      IF (source == 0) GOTO 602

      CALL inctag
      IF (myproc == source ) THEN
        itags = gentag + 4
        CALL mpsendi(m,1,destination,itags,ierr)
        itags = gentag + 5
        CALL mpsendi(n,1,destination,itags,ierr)


        itags = gentag
        CALL mpsendr(hterain,m*n,destination,itags,ierr)
        itags = gentag + 3
        CALL mpsendr(u,m*n,destination,itags,ierr)
        itags = gentag + 4
        CALL mpsendr(v,m*n,destination,itags,ierr)
        itags = gentag + 1
        CALL mpsendr(x,m*n,destination,itags,ierr)
        itags = gentag + 2
        CALL mpsendr(y,m*n,destination,itags,ierr)
      END IF

      plota = 0.0
      plotu = 0.0
      plotv = 0.0
      plotx = 0.0
      ploty = 0.0
      mm = 0
      nn = 0

      IF (myproc == 0) THEN
        itagr = gentag + 4
        CALL mprecvi(mm,1,source,itagr,ierr)
        itagr = gentag + 5
        CALL mprecvi(nn,1,source,itagr,ierr)

        itagr = gentag
        CALL mprecvr(plota,mm*nn,source,itagr,ierr)
        itagr = gentag + 3
        CALL mprecvr(plotu,mm*nn,source,itagr,ierr)
        itagr = gentag + 4
        CALL mprecvr(plotv,mm*nn,source,itagr,ierr)
        itagr = gentag + 1
        CALL mprecvr(plotx,mm*nn,source,itagr,ierr)
        itagr = gentag+2
        CALL mprecvr(ploty,mm*nn,source,itagr,ierr)
      END IF

      602 CONTINUE

      IF (myproc == 0) THEN

!-----------------------------------------------------------------------
!
!  Plot terrain profile in vertical slices
!
!-----------------------------------------------------------------------

        IF(slicopt == 2 .OR. slicopt == 3 .OR.slicopt == 5) THEN
          CALL xcolor(trcolor)
          CALL xthick(2)
          CALL xpenup( plotx(1), ploty(1)-0.5*(ploty(1+mm)-ploty(1)) )
          DO i=2,mm
            CALL xpendn(plotx(i), ploty(i)-0.5*(ploty(i+mm)-ploty(i)) )
          END DO
          CALL xthick(1)
        END IF
!
!-----------------------------------------------------------------------
!
!  Overlay terrain contour if required in x-y level
!      or Plot terrain outline in this slice zlevel .
!
!-----------------------------------------------------------------------
!
        IF(timeovr == 0) CALL plttrn(plota,plotx,ploty,mm,nn,slicopt,   &
                                     iwrk,xwk,ywk)

        CALL xcolor(lbcolor)

        CALL xwindw(xl, xr, yb, yt)

!
!-----------------------------------------------------------------------
!
!  Plot station labels
!
!-----------------------------------------------------------------------
!
        CALL xcolor(lbcolor)
        IF(ovrstaopt == 1 .AND. staset == 1 .AND.                       &
           (ovrstam == 1 .OR. ovrstan == 1 .OR. ovrstav == 1) .AND.     &
           (slicopt == 1 .OR. slicopt == 4 .OR. slicopt == 6 .OR.       &
            slicopt == 7 .OR. slicopt == 8 ) .AND. timeovr == 0  ) THEN

          CALL xchmag(0.025*px * lblmag)
          CALL pltsta(plotu,plotv,plotx,ploty,mm,nn,0,slicopt)
          !staset=0
        END IF
        IF (ovrstaopt == 1 .AND. wrtstax == 1 .AND. timeovr == 0 .AND.  &
            (slicopt == 2 .OR. slicopt == 3 .OR. slicopt == 5) ) THEN
          CALL xchmag(0.025*px * lblmag)
          flag=1
          CALL pltsta(plotu,plotv,plotx,ploty,mm,nn,flag,slicopt)
        END IF

        CALL xwdwof

      END IF  ! myproc == 0

      CALL mpbarrier         ! sync the processors

    END DO
  END DO

  IF(myproc == 0) THEN
!-----------------------------------------------------------------------
!
!  Plot observations
!
!-----------------------------------------------------------------------
!
    IF(ovrobs == 1 .AND. obsset == 1 .AND.                              &
       (slicopt == 1 .OR. slicopt == 4 .OR.                             &
        slicopt == 6 .OR. slicopt == 7 .OR. slicopt == 8)) THEN
      CALL pltobs(3)
      obsset=0
    END IF

!-----------------------------------------------------------------------
!
!  Plot additional text below the figure
!
!-----------------------------------------------------------------------

    CALL label2d(char2)

  END IF  ! myproc == 0

  cpu2 = f_cputime()
  second2 = f_walltime()

!  write(6,*) '!!!!  total cpu time for one VTR2D  :',                  &
!             cpu2-cpu1,' PLOT:',varname
  RETURN
END SUBROUTINE vtr2d
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE VTRUNT                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE vtrunt( vtunt0 )

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    Set the wind vector unit for wind field to be plotted by VTR2D.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!
!  MODIFICATION HISTORY:
!    6/08/92  Added full documentation (K. Brewster)
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    vtunt0    Unit vector
!              If VTUNT0 = 0.0, the unit is internally determined.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  REAL :: vtunt0
!
!-----------------------------------------------------------------------
!
!  Plotting control common blocks
!
!-----------------------------------------------------------------------
!
  REAL :: ctinc,ctmin,ctmax,vtunt   ! contour interval and vector unit
  COMMON /incunt/ ctinc,ctmin,ctmax,vtunt
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  vtunt = vtunt0

  RETURN
END SUBROUTINE vtrunt
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE STRM3D                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE strm3d(u,v,w, x,y,z, xw,xe,dx, ys,yn,dy, zb,zt,dz,           &
           nx,ibgn,iend,ist, ny,jbgn,jend,jst, nz,kbgn,kend,kst,        &
           kslice, jslice, islice, time, runname,factor,slicopt,        &
           n,xp,yp,zp,u1,v1,u2,v2,w2,                                   &
           tem1,tem2,tem3,tem4,tem5,                                    &
           tem6,hterain)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    Plot a streamline field in 2-d slices
!
!  AUTHOR: Ming Xue
!    1/16/1992
!
!  MODIFICATION HISTORY:
!
!  3/25/96 (K. Brewster)
!    Added variables isize,jsize,ksize
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    u        3-dimensional array of u wind components (m/s)
!    v        3-dimensional array of v wind components (m/s)
!    w        3-dimensional array of w wind components (m/s)
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in physcal space (m)
!
!    xw       value of x for first i grid point to plot
!    xe       value of x for last i grid point to plot
!    ys       value of y for first j grid point to plot
!    yn       value of y for last j grid point to plot
!    zb       value of z for first k grid point to plot
!    zt       value of z for last k grid point to plot
!
!    nx       first dimension of b
!    ibgn     index of first i grid point to plot
!    iend     index of last  i grid point to plot
!
!    ny       second dimension of b
!    jbgn     index of first j grid point to plot
!    jend     index of last  j grid point to plot
!
!    nz       third dimension of b
!    kbgn     index of first k grid point to plot
!    kend     index of last  k grid point to plot
!
!    ist      step size in x direction
!    jst      step size in y direction
!    kst      step size in z direction
!
!    time     time of data in seconds
!
!    kslice   k index of plane for slicopt=1 x-y slice
!    jslice   j index of plane for slicopt=2 x-z slice
!    islice   i index of plane for slicopt=1 y-z slice
!
!    runname  character string decribing run
!
!    factor   scaling factor for winds
!             V*factor wind vectors are plotted
!
!    slicopt  slice orientation indicator
!             = 1, x-y slice of at k=kslice is plotted.
!             = 2, x-z slice of at j=jslice is plotted.
!             = 3, y-z slice of at i=islice is plotted.
!             = 4, horizontal slice at z index islice is plotted.
!             = 5, xy-z cross section of wind islice is plotted.
!             = 6, data field on constant p-level is plotted.
!             = 0, all of the three slices above are plotted.
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
!    tem4     Temporary work array.
!    tem5     Temporary work array.
!    tem6     Temporary work array.
!
!  (These arrays are defined and used locally (i.e. inside this
!   subroutine), they may also be passed into routines called by
!   this one. Exiting the call to this subroutine, these temporary
!   work arrays may be used for other purposes therefore their
!   contents overwritten. Please examine the usage of work arrays
!   before you alter the code.)
!
!   pp01      The pressure (mb) value at the specific p-level
!   ercpl     reciprocal of exponent
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: nx,ny,nz
  INTEGER :: n
!
  REAL :: u(nx,ny,nz)
  REAL :: v(nx,ny,nz)
  REAL :: w(nx,ny,nz)

  REAL :: x(nx,ny,nz)
  REAL :: y(nx,ny,nz)
  REAL :: z(nx,ny,nz)
!
  REAL :: u1(nx,ny),v1(nx,ny)
  REAL :: u2(n,nz),v2(n,nz),w2(n,nz),zp(n,nz)
  REAL :: xp(n),yp(n)

  INTEGER :: kslice,jslice,islice
  CHARACTER (LEN=*) :: runname

  REAL :: xw,xe,dx,ys,yn,dy,zb,zt,dz
  INTEGER :: ibgn,iend,ist, jbgn,jend,jst, kbgn,kend,kst

  REAL :: time,factor
  INTEGER :: slicopt

  REAL :: x_tmp
  COMMON /tmphc2/ x_tmp

  INTEGER :: isource

!
!-----------------------------------------------------------------------
!
!  Some constants
!
!-----------------------------------------------------------------------
!
  REAL :: pp01, ercpl
  PARAMETER (ercpl=0.3678794)              ! exp(-1.0)
!
!-----------------------------------------------------------------------
!
!  Work arrays: tem1,tem2,tem3 of size at least
!               max( nx*ny, nx*nz, ny*nz).
!
!-----------------------------------------------------------------------
!
  REAL :: tem1(*),tem2(*),tem3(*),tem4(*),tem5(*),tem6(*)

  INTEGER, PARAMETER :: nzmax = 300
  REAL :: fdata(nzmax),zdata(nzmax),fprof(nzmax),zprof(nzmax)
!
!-----------------------------------------------------------------------
!
!  Common blocks for plotting control parameters
!
!-----------------------------------------------------------------------
!
  REAL :: x01,y01                  ! the first  point of interpolation
  REAL :: x02,y02                  ! the second point of interpolation
  REAL :: zlevel                   ! the given height of the slice
  REAL :: sinaf,cosaf,dist,sqrtdxy
  COMMON /slicev/x01,y01,x02,y02,sinaf,cosaf,dist,sqrtdxy
  COMMON /sliceh/zlevel
!
!-----------------------------------------------------------------------
!
!  Misc. local Variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,ij,ik,jk,length,isize,jsize,ksize
  CHARACTER (LEN=9)   :: timhms
  CHARACTER (LEN=120) :: title
!
  INTEGER :: icolor,icolor1,lbcolor,trcolor                ! required color
  COMMON /recolor/icolor,icolor1,lbcolor,trcolor
!
  INTEGER :: trnplt                ! flag to plot terain (1 or 0)
  REAL :: hterain(nx,ny)           ! The height of the terrain.

  INTEGER :: ovrtrn                ! overlay terrain option (0/1)
  REAL :: trninc,trnmin, trnmax    ! terrain interval minimum, maximum
  REAL :: ztmin,ztmax
  COMMON /trnpar/ trnplt,ovrtrn,trninc,trnmin, trnmax,ztmin,ztmax
!
  INTEGER :: xfont   ! the font of character
  INTEGER :: haxisu, vaxisu
  INTEGER :: lbaxis
  INTEGER :: tickopt
  INTEGER :: axlbfmt
  REAL :: hmintick,vmajtick,vmintick,hmajtick
  COMMON /var_par/ xfont,haxisu,vaxisu,lbaxis,tickopt,hmintick,         &
          vmajtick, vmintick,hmajtick,axlbfmt
  CHARACTER (LEN=4) :: stem2
  CHARACTER (LEN=1) :: stem1

  INTEGER :: smooth
  COMMON /smoothopt/smooth

  INTEGER :: wrtflag
  CHARACTER (LEN=90) :: levlab
  CHARACTER (LEN=50) :: timelab
  CHARACTER (LEN=25) :: timestring
  COMMON /timelev/wrtflag, timelab, levlab, timestring

  REAL :: tmpx, tmpy
  CHARACTER (LEN=20) :: distc
  REAL :: x101, y101, x102,y102
  COMMON /slicev1/x101, y101, x102,y102

  INTEGER :: xpbgn,xpend,ypbgn,ypend
  COMMON /processors/ xpbgn, xpend, ypbgn, ypend

!----------------------------------------------------------------------
!
! Include files
!
!---------------------------------------------------------------------

  INCLUDE 'mp.inc'

  REAL, ALLOCATABLE :: hterain_lg(:,:)
  REAL, ALLOCATABLE :: u_lg(:,:), v_lg(:,:), w_lg(:,:)
  REAL, ALLOCATABLE :: x_lg(:,:), y_lg(:,:), z_lg(:,:)

  INTEGER :: nxlg, nylg
  INTEGER :: istatus

  REAL,    ALLOCATABLE :: temht(:,:)
  REAL,    ALLOCATABLE :: temu(:,:), temv(:,:), temw(:,:)
  REAL,    ALLOCATABLE :: temx(:,:), temy(:,:), temz(:,:)
  INTEGER, ALLOCATABLE :: temwrk1(:,:)
  REAL,    ALLOCATABLE :: temwrk2(:,:,:), temwrk3(:,:)

  INTEGER :: iwrk
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  IF (mp_opt >0) THEN
    !WRITE(6,'(2a/,a/)') 'Sorry, since subroutine strmln is an ',        &
    !         'internal procedure of NCARG package. It is not MPI''ed.', &
    !         'No streamline field is plotted.'
    !RETURN
    nxlg = (nx-3)*nproc_x + 3
    nylg = (ny-3)*nproc_y + 3

!  END IF

  isize=(iend-ibgn)/ist+1
  jsize=(jend-jbgn)/jst+1
  ksize=(kend-kbgn)/kst+1

  iwrk = 2                  ! Does not color magnitude
  IF (icolor1-icolor > 0) iwrk = 3

!
!-----------------------------------------------------------------------
!
!  setup  time label
!
!-----------------------------------------------------------------------
!
  CALL get_forecast_hms( time, timhms )
  WRITE(timelab,'(a,F8.1,a)') 'T=',time,' s ('//TRIM(timhms)//')'

  CALL get_time_string ( time, timestring,'Z  ',0 )
!
!-----------------------------------------------------------------------
!
!  Set up terrain, if needed.
!
!-----------------------------------------------------------------------
!
  IF(trnplt == 1 .OR.trnplt == 2 .OR. ovrtrn == 1)  THEN

    ALLOCATE(hterain_lg(nxlg,nylg),   STAT = istatus)
    ALLOCATE(temht     (isize,jsize), STAT = istatus)

    CALL mpimerge2d(hterain,nx,ny,hterain_lg)

    jk = 0
    DO j=jbgn,jend,jst
      jk = jk+1
      ik = 0
      DO i=ibgn,iend,ist
        ik = ik+1
        temht(ik,jk)=hterain_lg(i,j)
      END DO
    END DO
  END IF

  IF ( slicopt == 2 .OR. slicopt == 3  .OR. slicopt == 5) THEN
    CALL cal_dist(haxisu,dx,dy,x01,y01,x02,y02,slicopt,                 &
                  tmpx,tmpy,distc)
  END IF

!
!-----------------------------------------------------------------------
!
!  slicopt=1   Plot u-v field
!
!-----------------------------------------------------------------------
!
  IF( slicopt == 1 .OR. slicopt == 0 ) THEN

    ALLOCATE(u_lg (nxlg, nylg),  STAT = istatus)
    ALLOCATE(v_lg (nxlg, nylg),  STAT = istatus)
    ALLOCATE(x_lg (nxlg, nylg),  STAT = istatus)
    ALLOCATE(y_lg (nxlg, nylg),  STAT = istatus)

    ALLOCATE(temu (isize,jsize), STAT = istatus)
    ALLOCATE(temv (isize,jsize), STAT = istatus)
    ALLOCATE(temx (isize,jsize), STAT = istatus)
    ALLOCATE(temy (isize,jsize), STAT = istatus)

    ALLOCATE(temwrk1 (isize,jsize),      STAT = istatus)
    ALLOCATE(temwrk2 (isize,jsize,iwrk), STAT = istatus)
    ALLOCATE(temwrk3 (isize,jsize),      STAT = istatus)

    k = kslice
    CALL mpimerge2d(u(:,:,k),nx,ny,u_lg)
    CALL mpimerge2d(v(:,:,k),nx,ny,v_lg)
    CALL mpimerge2d(x(:,:,k),nx,ny,x_lg)
    CALL mpimerge2d(y(:,:,k),nx,ny,y_lg)

    IF (myproc == 0) THEN
      jk = 0
      DO j=jbgn,jend,jst
        jk = jk + 1
        ik = 0
        DO i=ibgn,iend,ist
          ik = ik + 1
          temu(ik,jk) = -9999.0
          temv(ik,jk) = -9999.0
          IF(u_lg(i,j) /= -9999.0) temu(ik,jk)=u_lg(i,j)*factor
          IF(v_lg(i,j) /= -9999.0) temv(ik,jk)=v_lg(i,j)*factor
          temx(ik,jk)=x_lg(i,j)
          temy(ik,jk)=y_lg(i,j)
        END DO
      END DO

      IF (k /= 2) THEN
        WRITE(title,'(''U-V Streamline'')')
        WRITE(levlab,'(''X-Y cross-section through k='',I3)')k
      ELSE
        WRITE(title,'(''U-V Streamline'')')
        WRITE(levlab,'(''X-Y cross-section through k=2 (surface)'')')
      END IF

      length = 120
      CALL strlnth( title, length )
      CALL strmin ( title, length )

      DO i=1,smooth
        CALL smooth9pmv(temu,isize,jsize,1,isize,1,jsize,temwrk2)
        CALL smooth9pmv(temv,isize,jsize,1,isize,1,jsize,temwrk2)
      END DO

      CALL strm2d(temu,temv, xw,xe,ys,yn, dx, dy, isize,jsize,            &
           title(1:length),runname, temx,temy,                            &
           temht,slicopt,temwrk1,temwrk2,temwrk3)
    END IF

!
!-----------------------------------------------------------------------
!
!  slicopt=2   Plot u-w streamline
!
!-----------------------------------------------------------------------
!
  ELSE IF( slicopt == 2 .OR. slicopt == 0 ) THEN

    ALLOCATE(u_lg(nxlg,nz),       STAT = istatus)
    ALLOCATE(w_lg(nxlg,nz),       STAT = istatus)
    ALLOCATE(x_lg(nxlg,nz),       STAT = istatus)
    ALLOCATE(z_lg(nxlg,nz),       STAT = istatus)

    ALLOCATE(temu (isize,ksize),    STAT = istatus)
    ALLOCATE(temw (isize,ksize),    STAT = istatus)
    ALLOCATE(temx (isize,ksize),    STAT = istatus)
    ALLOCATE(temz (isize,ksize),    STAT = istatus)

    ALLOCATE(temwrk1 (isize,ksize),      STAT = istatus)
    ALLOCATE(temwrk2 (isize,ksize,iwrk), STAT = istatus)
    ALLOCATE(temwrk3 (isize,ksize),      STAT = istatus)

    j = jslice   ! local index

    DO k=1,nz
      DO i=1,nx
        ik = i + (k-1)*nx
        tem1(ik) = u(i,j,k)
        tem2(ik) = w(i,j,k)
        tem3(ik) = x(i,j,k)
        tem4(ik) = z(i,j,k)
      END DO
    END DO

    CALL mpimerge2dx(tem1,nx,nz,ypbgn,u_lg,istatus)
    CALL mpimerge2dx(tem2,nx,nz,ypbgn,w_lg,istatus)
    CALL mpimerge2dx(tem3,nx,nz,ypbgn,x_lg,istatus)
    CALL mpimerge2dx(tem4,nx,nz,ypbgn,z_lg,istatus)

    x_tmp = y(1,jslice,1)
    isource = (ypbgn-1)*nproc_x + xpbgn-1
    IF (isource /= 0) THEN
      IF ( myproc == isource ) CALL mpsendr(x_tmp,1,      0,123,istatus)
      IF ( myproc == 0)        CALL mprecvr(x_tmp,1,isource,123,istatus)
    END IF

    IF (myproc == 0) THEN

      jk = 0
      DO k=kbgn,kend,kst
        jk = jk + 1
        ik = 0
        DO i=ibgn,iend,ist
          ik = ik + 1
          temu(ik,jk) = -9999.0
          temw(ik,jk) = -9999.0
          IF(u_lg(i,k) /= -9999.0) temu(ik,jk)=u_lg(i,k)*factor
          IF(w_lg(i,k) /= -9999.0) temw(ik,jk)=w_lg(i,k)*factor
          temx(ik,jk)=x_lg(i,k)
          temz(ik,jk)=z_lg(i,k)
        END DO
      END DO

      IF( nzmax < ksize) THEN
        WRITE(6,'(1x,a)')                                                 &
            'nzmax given in STRM3D too small. Job stopped.'
        istatus = -1
        GO TO 999
      END IF

      DO k=1,ksize
        zprof(k)= zb+(k-1)*(zt-zb)/(kend-kbgn)
      END DO

      CALL unigrid(isize,ksize,temu,temz,fdata,zdata,fprof,zprof)
      CALL unigrid(isize,ksize,temw,temz,fdata,zdata,fprof,zprof)

      WRITE(title,'(''U-W Streamline '')')
      j = j + (ypbgn-1)*(ny-3)    ! Changed back to global index
      dist = (j-1)*tmpy
      length=LEN_TRIM(distc)
      CALL strmin ( distc, length)
      WRITE(levlab,'(''X-Z cross-section through j='',I3,                 &
                &  '' (y = '',f8.1,a,a)')j,dist,distc(1:length),')'

      length = 120
      CALL strlnth( title, length )
      CALL strmin ( title, length)

      DO i=1,smooth
        CALL smooth9pmv(temu,isize,ksize,1,isize,1,ksize,temwrk2)
        CALL smooth9pmv(temw,isize,ksize,1,isize,1,ksize,temwrk2)
      END DO

      CALL strm2d(temu,temw, xw,xe,zb,zt, dx, dz, isize,ksize,            &
           title(1:length),runname, temx ,temz,                           &
           temht,slicopt,temwrk1,temwrk2,temwrk3)

      999 CONTINUE
    END IF
    CALL mpupdatei(istatus,1)
    IF (istatus < 0) CALL arpsstop('ERROR: inside strm3d.',1)

!
!-----------------------------------------------------------------------
!
!  slicopt=3   Plot v-w field
!
!-----------------------------------------------------------------------
!
  ELSE IF( slicopt == 3 .OR. slicopt == 0 ) THEN

    ALLOCATE(v_lg(nylg,nz),       STAT = istatus)
    ALLOCATE(w_lg(nylg,nz),       STAT = istatus)
    ALLOCATE(y_lg(nylg,nz),       STAT = istatus)
    ALLOCATE(z_lg(nylg,nz),       STAT = istatus)

    ALLOCATE(temv (jsize,ksize),    STAT = istatus)
    ALLOCATE(temw (jsize,ksize),    STAT = istatus)
    ALLOCATE(temy (jsize,ksize),    STAT = istatus)
    ALLOCATE(temz (jsize,ksize),    STAT = istatus)

    ALLOCATE(temwrk1 (jsize,ksize),      STAT = istatus)
    ALLOCATE(temwrk2 (jsize,ksize,iwrk), STAT = istatus)
    ALLOCATE(temwrk3 (jsize,ksize),      STAT = istatus)

    i = islice            ! local index

    DO k=1,nz
      DO j=1,ny
        jk = j + (k-1)*ny
        tem1(jk) = v(i,j,k)
        tem2(jk) = w(i,j,k)
        tem3(jk) = y(i,j,k)
        tem4(jk) = z(i,j,k)
      END DO
    END DO

    CALL mpimerge2dy(tem1,ny,nz,xpbgn,v_lg,istatus)
    CALL mpimerge2dy(tem2,ny,nz,xpbgn,w_lg,istatus)
    CALL mpimerge2dy(tem3,ny,nz,xpbgn,y_lg,istatus)
    CALL mpimerge2dy(tem4,ny,nz,xpbgn,z_lg,istatus)

    x_tmp = x(islice,1,1)
    isource = (ypbgn-1)*nproc_x + xpbgn-1
    IF (isource /= 0) THEN
      IF ( myproc == isource ) CALL mpsendr(x_tmp,1,      0,123,istatus)
      IF ( myproc == 0)        CALL mprecvr(x_tmp,1,isource,123,istatus)
    END IF

    IF (myproc == 0) THEN

      jk = 0
      DO k=kbgn,kend,kst
        jk = jk + 1
        ik = 0
        DO j=jbgn,jend,jst
          ik = ik + 1
          temv(ik,jk) = -9999.0
          temw(ik,jk) = -9999.0
          IF(v_lg(j,k) /= -9999.0) temv(ik,jk)=v_lg(j,k)*factor
          IF(w_lg(j,k) /= -9999.0) temw(ik,jk)=w_lg(j,k)*factor
          temy(ik,jk)=y_lg(j,k)
          temz(ik,jk)=z_lg(j,k)
        END DO
      END DO

      IF( nzmax < ksize) THEN
        WRITE(6,'(1x,a)')                                                 &
            'nzmax given in STRM3D too small. Job stopped.'
        istatus = -1
        GO TO 998
      END IF

      DO k=1,ksize
        zprof(k)= zb+(k-1)*(zt-zb)/(kend-kbgn)
      END DO

      CALL unigrid(jsize,ksize,temv,temz,fdata,zdata,fprof,zprof)
      CALL unigrid(jsize,ksize,temw,temz,fdata,zdata,fprof,zprof)

      i = i + (xpbgn-1)*(nx-3)   ! Changed back to global index
      dist = (i-1)*tmpx
      length=LEN_TRIM(distc)
      CALL strmin ( distc, length)
      WRITE(levlab,'(''Y-Z cross-section through i='',I3,                 &
      &   '' ( x='',f8.1,a,a)')i,dist, distc(1:length),')'
      WRITE(title,'(''V-W Streamline'')')

      length = 120
      CALL strlnth( title, length )
      CALL strmin ( title, length )

      DO i=1,smooth
        CALL smooth9pmv(temv,jsize,ksize,1,jsize,1,ksize,temwrk2)
        CALL smooth9pmv(temw,jsize,ksize,1,jsize,1,ksize,temwrk2)
      END DO

      CALL strm2d(temv,temw, ys,yn,zb,zt, dy, dz,jsize,ksize,             &
           title(1:length),runname, temy ,temz,                           &
           temht,slicopt,temwrk1,temwrk2,temwrk3)

      998 CONTINUE
    END IF
    CALL mpupdatei(istatus,1)
    IF (istatus < 0) CALL arpsstop('ERROR: inside strm3d.',1)

!
!-----------------------------------------------------------------------
!
!  slicopt=4   Plot u-v streamlines on constant z levels
!  slicopt=6   Plot u-v streamlines on constant pressure levels
!  slicopt=7   Plot u-v streamlines on constant PT levels
!
!-----------------------------------------------------------------------
!
  ELSE IF( slicopt == 4 .OR. slicopt == 6 .OR. slicopt == 7 ) THEN

!    CALL hintrp(nx,ny,nz,u,z,zlevel,u1)
!    CALL hintrp(nx,ny,nz,v,z,zlevel,v1)

    CALL hintrp1(nx,ny,nz,kbgn,kend,u,z,zlevel,tem1)
    CALL hintrp1(nx,ny,nz,kbgn,kend,v,z,zlevel,tem2)

    ALLOCATE(u_lg(nxlg,nylg),       STAT = istatus)
    ALLOCATE(v_lg(nxlg,nylg),       STAT = istatus)
    ALLOCATE(x_lg(nxlg,nylg),       STAT = istatus)
    ALLOCATE(y_lg(nxlg,nylg),       STAT = istatus)

    ALLOCATE(temu (isize,jsize),    STAT = istatus)
    ALLOCATE(temv (isize,jsize),    STAT = istatus)
    ALLOCATE(temx (isize,jsize),    STAT = istatus)
    ALLOCATE(temy (isize,jsize),    STAT = istatus)

    ALLOCATE(temwrk1 (isize,jsize),      STAT = istatus)
    ALLOCATE(temwrk2 (isize,jsize,iwrk), STAT = istatus)
    ALLOCATE(temwrk3 (isize,jsize),      STAT = istatus)

    CALL mpimerge2d(tem1,nx,ny,u_lg)
    CALL mpimerge2d(tem2,nx,ny,v_lg)
    CALL mpimerge2d(x(:,:,2),nx,ny,x_lg)
    CALL mpimerge2d(y(:,:,2),nx,ny,y_lg)

    IF (myproc == 0) THEN

      jk = 0
      DO j=jbgn,jend,jst
        jk = jk + 1
        ik = 0
        DO i=ibgn,iend,ist
          ik = ik + 1
          temu(ik,jk) = -9999.0
          temv(ik,jk) = -9999.0
          IF(u_lg(i,j) /= -9999.0) temu(ik,jk)=u_lg(i,j)*factor
          IF(v_lg(i,j) /= -9999.0) temv(ik,jk)=v_lg(i,j)*factor
          temx(ik,jk)=x_lg(i,j)
          temy(ik,jk)=y_lg(i,j)
        END DO
      END DO

      IF( slicopt == 4) THEN
        WRITE(levlab,'(''Z='',F7.3,A,'' MSL'')')zlevel,' km'
      ELSE IF( slicopt == 6) THEN
        pp01=0.01*ercpl**zlevel
        WRITE(levlab,'(''P='',F7.2,a)') pp01, ' MB'
      ELSE
        WRITE(levlab,'(''Theta='',F5.1,a)') zlevel, ' (K)'
      END IF

      WRITE(title,'(''U-V Streamline'')')

      length = 120
      CALL strlnth( title, length )
      CALL strmin ( title, length)

      DO i=1,smooth
        CALL smooth9pmv(temu,isize,jsize,1,isize,1,jsize,temwrk2)
        CALL smooth9pmv(temv,isize,jsize,1,isize,1,jsize,temwrk2)
      END DO

      CALL strm2d(temu,temv, xw,xe,ys,yn, dx, dy, isize,jsize,            &
           title(1:length),runname, temx ,temy,                           &
           temht,slicopt,temwrk1,temwrk2,temwrk3)
    END IF
!
!-----------------------------------------------------------------------
!
!  slicopt=5   Plot V-w field
!
!-----------------------------------------------------------------------
!
  ELSE IF( slicopt == 5 ) THEN

    CALL sectvrt(nx,ny,nz,u,x,y,z,dx,dy,u2,zp,n,xp,yp)
    CALL sectvrt(nx,ny,nz,v,x,y,z,dx,dy,v2,zp,n,xp,yp)
    CALL sectvrt(nx,ny,nz,w,x,y,z,dx,dy,w2,zp,n,xp,yp)

    ALLOCATE(temu (n,ksize),    STAT = istatus)
    ALLOCATE(temw (n,ksize),    STAT = istatus)
    ALLOCATE(temx (n,ksize),    STAT = istatus)
    ALLOCATE(temz (n,ksize),    STAT = istatus)

    ALLOCATE(temwrk1 (n,ksize),      STAT = istatus)
    ALLOCATE(temwrk2 (n,ksize,iwrk), STAT = istatus)
    ALLOCATE(temwrk3 (n,ksize),      STAT = istatus)

    IF (myproc == 0) THEN

      jk = 0
      DO k=kbgn,kend,kst
        jk = jk + 1
        ik = 0
        DO i=1,n
          ik = ik + 1
          temu(ik,jk) = -9999.0
          temw(ik,jk) = -9999.0
          IF(u2(i,k) /= -9999.0 .AND. v2(i,k) /= -9999.0)               &
            temu(ik,jk) = (u2(i,k)*cosaf+v2(i,k)*sinaf)*factor
          IF(w2(i,k) /= -9999.0)                                        &
            temw(ik,jk) = w2(i,k)*factor
          temx(ik,jk) = xw+(i-ibgn)*sqrtdxy
          temz(ik,jk) = zp(i,k)
        END DO
      END DO

      IF( nzmax < ksize) THEN
        WRITE(6,'(1x,a)') 'nzmax given in STRM3D too small. Job stopped.'
        istatus = -1
        GO TO 997
      END IF

      DO k=1,ksize
        zprof(k)= zb+(k-1)*(zt-zb)/(kend-kbgn)
      END DO

      CALL unigrid(n,ksize,temu,temz,fdata,zdata,fprof,zprof)
      CALL unigrid(n,ksize,temw,temz,fdata,zdata,fprof,zprof)

      IF(axlbfmt == -1 .OR. axlbfmt == 1 ) THEN
        length=LEN_TRIM(distc)
        CALL strmin(distc,length)
        WRITE(title,'(''V-W streamline'')')
        WRITE(levlab,                                                   &
            '(''Vert. cross-section through '',4(A,F8.1),A,A)')         &
            '(',x101,',',y101,') (',x102,',',y102,')',distc(1:length)
      ELSE IF(axlbfmt == 0 ) THEN
        length=LEN_TRIM(distc)
        CALL strmin(distc,length)
        WRITE(title,'(''V-W streamline'')')
        WRITE(levlab,                                                   &
            '(''Vert. cross-section through '',4(A,I5),A,A)')           &
            '(',NINT(x101),',',NINT(y101),') (',NINT(x102),',',         &
                NINT(y102),')',distc(1:length)
      ELSE
!       WRITE(stem1,'(i1)')axlbfmt
!       WRITE(stem2,'(a3,a1)')'f8.',stem1
        WRITE(title,'(''V-W streamline'')')
        WRITE(levlab,                                                   &
            '(''Vert. cross-section through '',4(A,f8.2),A,A)')         &
            '(',x101,',',y101,') (',x102,',',y102,')',distc(1:length)

      END IF

      length = 120
      CALL strlnth( title, length )
      CALL strmin ( title, length)

      DO i=1,smooth
        CALL smooth9pmv(temu,n,ksize,1,n,1,ksize,temwrk2)
        CALL smooth9pmv(temw,n,ksize,1,n,1,ksize,temwrk2)
      END DO

      CALL strm2d(temu,temw, xw,xe,zb,zt, sqrtdxy, dz, n,ksize,         &
           title(1:length),runname, temx,temz,                          &
           temht,slicopt,temwrk1,temwrk2,temwrk3)

      997 CONTINUE
    END IF
    CALL mpupdatei(istatus,1)
    IF (istatus < 0) CALL arpsstop('ERROR: inside strm3d.',1)

  END IF

!-----------------------------------------------------------------------
!
! Before returning
!
!-----------------------------------------------------------------------

  IF (ALLOCATED(hterain_lg)) DEALLOCATE(hterain_lg)
  IF (ALLOCATED(temht))      DEALLOCATE(temht)

  IF (ALLOCATED(u_lg))     DEALLOCATE(u_lg)
  IF (ALLOCATED(v_lg))     DEALLOCATE(v_lg)
  IF (ALLOCATED(w_lg))     DEALLOCATE(w_lg)
  IF (ALLOCATED(x_lg))     DEALLOCATE(x_lg)
  IF (ALLOCATED(y_lg))     DEALLOCATE(y_lg)
  IF (ALLOCATED(z_lg))     DEALLOCATE(z_lg)

  IF (ALLOCATED(temu))     DEALLOCATE(temu)
  IF (ALLOCATED(temv))     DEALLOCATE(temv)
  IF (ALLOCATED(temw))     DEALLOCATE(temw)
  IF (ALLOCATED(temx))     DEALLOCATE(temx)
  IF (ALLOCATED(temy))     DEALLOCATE(temy, STAT = istatus)
  IF (ALLOCATED(temz))     DEALLOCATE(temz)

  IF (ALLOCATED(temwrk1))  DEALLOCATE(temwrk1)
  IF (ALLOCATED(temwrk2))  DEALLOCATE(temwrk2)
  IF (ALLOCATED(temwrk3))  DEALLOCATE(temwrk3)

  RETURN
END SUBROUTINE strm3d
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE STRM2D                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE strm2d(u,v,xl,xr,yb,yt,dx,dy,m,n,char1,char2, x,y,           &
                  hterain,slicopt,iwrk,xwk,ywk)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    Plot streamlines of a 2-d wind (u1,u2) field using ncargraphic
!    subroutine strmln
!
!  INPUT:
!    u        m by n 2-dimensional array of u (left-to-right)
!               wind components (m/s)
!    v        m by n 2-dimensional array of v (down-to-up)
!               wind components (m/s)
!
!    xl,xr    The left and right bound of the physical domain.
!    yb,yt    Bottom and top bound of the physical domain.
!    dx,dy    Grid interval in x and y direction (km)
!
!    m        First dimension of vector component array
!    n        Second dimension of vector component array
!
!    char1    First character string to plot (title)
!    char2    Second character string to plot (runname)
!
!    x        x coordinate of grid points in plot space (over on page)
!    y        y coordinate of grid points in plot space (up on page)
!
!    hterain  the height of terrain
!    slicopt  slice orientation indicator
!       slicopt = 1, x-y slice of u,v at z index kslice is plotted.
!       slicopt = 2, x-z slice of u,w at y index jslice is plotted.
!       slicopt = 3, y-z slice of v,w at x index islice is plotted.
!       slicopt = 4, x-y slice of u,v at z index islice is plotted.
!       slicopt = 5, xy-z cross section of wind islice is plotted.
!       slicopt = 6, data field on constant p-level is plotted.
!       slicopt = 0, all of the three slices above are plotted.
!
!  WORK ARRAY
!    iwrk      A work array of size at least m*n*2
!    xwk       A work array of size at least m*n*2
!    ywk
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INCLUDE 'arpsplt.inc'

  INTEGER :: m,n
  INTEGER :: i

  REAL :: u(m,n)
  REAL :: v(m,n)

  REAL :: x(m,n)
  REAL :: y(m,n)
  REAL :: xl,xr,yb,yt,dx,dy

  INTEGER, INTENT(INOUT) :: iwrk(m,n)
  REAL,    INTENT(INOUT) :: xwk(m,n,2), ywk(m,n)

  CHARACTER (LEN=*) :: char2
  CHARACTER (LEN=*) :: char1
  INTEGER :: ierror
!
!-----------------------------------------------------------------------
!
!  Plotting control common blocks
!
!-----------------------------------------------------------------------
!
  INTEGER :: layover
  COMMON /laypar/ layover

  REAL :: ctinc,ctmin,ctmax,vtunt  !contour interval and vector unit
  COMMON /incunt/ ctinc,ctmin,ctmax,vtunt

  INTEGER :: icolor,icolor1,lbcolor,trcolor                ! required color
  COMMON /recolor/icolor,icolor1,lbcolor,trcolor

  INTEGER :: flag
  INTEGER :: xfont   ! the font of character
  INTEGER :: haxisu, vaxisu
  INTEGER :: lbaxis
  INTEGER :: tickopt
  INTEGER :: axlbfmt
  REAL :: hmintick,vmajtick,vmintick,hmajtick
  COMMON /var_par/ xfont,haxisu,vaxisu,lbaxis,tickopt,hmintick,         &
         vmajtick, vmintick,hmajtick,axlbfmt

  INTEGER :: ovrstaopt
  INTEGER :: ovrstam,staset,ovrstan,ovrstav,stacol,markprio,wrtstax
  INTEGER :: nsta_typ,sta_typ(30),sta_marktyp(30),sta_markcol(30)
  REAL :: sta_marksz(30),wrtstad
  CHARACTER (LEN=256) :: stalofl
  COMMON /sta_par/ ovrstaopt,ovrstam,staset,ovrstan,ovrstav,stacol,     &
         markprio,nsta_typ,sta_typ,sta_marktyp,                         &
         sta_markcol,sta_marksz,stalofl,wrtstax,wrtstad

  REAL :: yxratio                  !the scaling factor the y/x ratio.
  COMMON /yratio/ yxratio

  INTEGER :: col_table,pcolbar
  COMMON /coltable/col_table,pcolbar
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: nopic,nhpic,nvpic,ifont
  REAL :: pl,pr,pb,pt   ! plot space left, right, bottom, top coordinate
  REAL :: px,py         ! plot space left-right length and up-down height
  REAL :: xs,ys         ! real space left-right length and up-down height
  REAL :: pxc,pyc       ! plot space left-right center and
                        ! up-down    center
  REAL :: xlimit,ylimit
  REAL :: rotang
  REAL :: xp1,xp2,yp1,yp2
  REAL :: xd1,xd2,yd1,yd2,xpos1,xpos2,ypos1,ypos2

  REAL :: zlevel
  COMMON/sliceh/zlevel

  REAL :: hterain(m,n)       ! The height of the terrain.

  INTEGER :: slicopt

  INTEGER :: timeovr
  COMMON /timover/ timeovr

  REAL :: lblmag, ctrlbsiz, axlbsiz
  COMMON /labmag/ lblmag, ctrlbsiz, axlbsiz

  INTEGER :: len1

  INTEGER :: wrtflag
  CHARACTER (LEN=90) :: levlab
  CHARACTER (LEN=50) :: timelab
  CHARACTER (LEN=25) :: timestring
  COMMON /timelev/wrtflag, timelab, levlab, timestring

  REAL :: x101, y101, x102,y102
  COMMON /slicev1/x101, y101, x102,y102

  INTEGER :: xnwpic_called
  COMMON /callnwpic/xnwpic_called

  REAL    :: yltmp   ! local temporary variable
  INTEGER :: nset

  CHARACTER(LEN=80) :: f_char1
  INTEGER           :: LEN0

  REAL :: xttmp, yttmp
  REAL :: magmax, magmin
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  WRITE(6,'(/1x,a,a)') 'Plotting ',char1

  IF( layover == 0 .OR. xnwpic_called == 0) THEN
    CALL xnwpic
    xnwpic_called=1
    timeovr=0
    wrtflag = 0
  ELSE
    timeovr=1
    wrtflag = wrtflag + 1
  END IF
!
!-----------------------------------------------------------------------
!
!  Get plotting space variables
!
!-----------------------------------------------------------------------
!
  CALL xqpspc( pl, pr, pb, pt)
  px = pr - pl
  py = pt - pb

  xs = xr-xl
  ys = yt-yb

  pxc = (pr+pl)/2
  pyc = (pb+pt)/2
!
!-----------------------------------------------------------------------
!
!  Let the longest lenth determine size scaling of plot
!
!-----------------------------------------------------------------------
!
  IF( py/px >= ys*yxratio/xs ) THEN
    py = ys*yxratio/xs*px
    CALL xpspac(pl, pr, pyc-py/2, pyc+py/2 )
  ELSE
    px = xs/(ys*yxratio)*py
    CALL xpspac(pxc-px/2, pxc+px/2, pb, pt)
  END IF
!
!-----------------------------------------------------------------------
!
!  Set the real distance to plot distance scaling
!
!-----------------------------------------------------------------------
!
  CALL xmap( xl, xr, yb,yt)
!
!-----------------------------------------------------------------------
!
!  Plot map, boxes and polygons.
!
!-----------------------------------------------------------------------
!
  CALL xcolor(lbcolor)

  CALL pltextra(slicopt, 1 )

  xpos1 = xl
  xpos2 = xr
  ypos1 = yb
  ypos2 = yt

  CALL xtrans(xpos1,ypos1)
  CALL xtrans(xpos2,ypos2)
  CALL xzx2ncar(xpos1,ypos1)
  CALL xzx2ncar(xpos2,ypos2)
!
  IF(slicopt == 2 .OR. slicopt == 3 .OR.slicopt == 5) THEN
    CALL xcolor(trcolor)
    CALL xthick(3)
    CALL xpenup( x(1,1), y(1,1)-0.5*(y(1,2)-y(1,1)) )
    DO i=2,m
      CALL xpendn( x(i,1), y(i,1)-0.5*(y(i,2)-y(i,1)) )
    END DO
    CALL xthick(1)
  END IF
!
!-----------------------------------------------------------------------
!
!  Overlay terrain contour if required in x-y level
!  or Plot terrain outline in this slice zlevel .
!
!-----------------------------------------------------------------------
!
  IF( timeovr == 0) CALL plttrn(hterain,x,y,m,n,slicopt,iwrk,xwk,ywk)

  CALL xcolor(lbcolor)
!
!-----------------------------------------------------------------------
!
!  Plot station labels
!
!-----------------------------------------------------------------------
!
  IF(ovrstaopt == 1 .AND. staset == 1 .AND.                             &
        (ovrstam == 1.OR.ovrstan == 1.OR.ovrstav == 1).AND.             &
        (slicopt == 1.OR.slicopt == 4.OR.slicopt == 6                   &
        .OR.slicopt == 7.OR.slicopt == 8)                               &
        .AND.timeovr == 0 ) THEN
    CALL xchmag(0.025*px * lblmag)
    CALL pltsta(u,v,x,y,m,n,0,slicopt)
!    staset=0
  END IF
!
!-----------------------------------------------------------------------
!
!  Plot observations
!
!-----------------------------------------------------------------------
!
  IF( ovrstaopt == 1 .AND. wrtstax == 1 .AND. timeovr == 0              &
        .AND.(slicopt == 2.OR.slicopt == 3.OR. slicopt == 5) ) THEN
    CALL xchmag(0.025*px * lblmag)
    flag=1
    CALL pltsta(u,v,x,y,m,n,flag,slicopt)
  END IF
!
!-----------------------------------------------------------------------
!
!  Plot streamlines
!
!-----------------------------------------------------------------------
!
  CALL xcolor(lbcolor)
!
  CALL xqset(xp1,xp2,yp1,yp2, xd1,xd2,yd1,yd2)
  CALL set(xpos1,xpos2,ypos1,ypos2, 1.0, FLOAT(m), 1.0, FLOAT(n),1)

  CALL xcolor(icolor)

  nset = (icolor1-icolor)+1
  CALL strmln_new(u,v, xwk, m, m, n, nset, magmax, magmin, ierror)
  CALL set(xp1,xp2,yp1,yp2, xd1,xd2,yd1,yd2, 1)
!
!-----------------------------------------------------------------------
!
!  Plot axes with tick marks
!
!-----------------------------------------------------------------------
!
  CALL xcolor(lbcolor)

  CALL pltaxes(slicopt,dx,dy)

!-----------------------------------------------------------------------
!
!  Plot labels
!
!-----------------------------------------------------------------------

  CALL xcolor(lbcolor)

  CALL xqnpic(nopic)
  CALL xqspac(nhpic, nvpic, rotang, xlimit, ylimit)

! write time and level
  CALL xchmag( 0.025*px  * lblmag )
  IF ( layover < 1 ) THEN

    len1=LEN_TRIM(timelab)
    CALL strmin(timelab,len1)
    CALL xcharc((xl+xr)*0.5,yt+0.015*ys*px/py,                            &
                      timestring(1:25)//' '//timelab(1:len1))

    len1=LEN_TRIM(levlab)
    CALL strmin(levlab,len1)
    CALL xcharc(xl+xs*0.5,yt+0.06*ys*px/py, levlab(1:len1))
  END IF

! write variable label
  IF (nset <= 1) CALL xcolor(icolor)
  IF(lbaxis == 1) THEN
    IF(wrtstax == 0) THEN
      yltmp = 0.08
    ELSE
      yltmp = 0.14
    END IF
  ELSE
    yltmp =0.12
  END IF

  CALL xchmag( 0.025*px  * lblmag )
  CALL xcharl(xl-0.05*(xr-xl), yb-(yt-yb)*px/py*(yltmp+layover*0.030),char1)
!  CALL xcharc(0.5*(xr+xl), yb-(yt-yb)*px/py*(yltmp+layover*0.030),char1)

  CALL xcolor(lbcolor)
  IF (nset > 1) THEN
    LEN0 = LEN_TRIM(char1)
    CALL strmin(char1, LEN0)

    IF( char1(LEN0:LEN0) == ')' ) LEN0 = max(1,LEN0-1)
    WRITE(f_char1, '(a,a)') char1(1:LEN0),' (m/s, magnitude coloring)'

    LEN0=LEN_TRIM(f_char1)
    CALL strmin(f_char1,LEN0)

    CALL xchmag( 0.025*px * lblmag )
    CALL xcolor(lbcolor)

    xttmp = xl-0.05*(xr-xl)
    yttmp = yb-(yltmp+wrtflag*0.030)*ys*px/py
    CALL xcharl(xttmp, yttmp, f_char1(1:LEN0))

    f_char1 = ' '
    WRITE(f_char1,'(a,G9.3E2,a,G9.3E2)') 'Min=',magmin,' Max=',magmax

    xttmp = xr+0.05*(xr-xl)
    yttmp = yb-(yltmp+wrtflag*0.030)*ys*px/py
    LEN0=LEN_TRIM(f_char1)
    CALL strmin(f_char1,LEN0)
    CALL xcharr(xttmp, yttmp, f_char1(1:LEN0))
  END IF

  IF (timeovr == 0) THEN
    IF(nopic == nhpic*(nvpic-1)+1 ) THEN
      !ytmp =0.25
      !IF(layover < 1) CALL xcharl(xl,yb-(ytmp+layover*0.03)*(yt-yb)*px/py, char2 )

      !CALL xqcfnt(ifont)
      !CALL xcfont(xfont)
      !ytmp = 0.20
!     !IF(layover < 1) CALL xcharl(xl,yb-(0.20+layover*0.03)*(yt-yb)*px/py,'CAPS/ARPS ' )
      !
      !CALL xcfont(ifont)

      IF (layover < 1) CALL label2d(char2)
    END IF
    timeovr=1
  END IF

!----------------------------------------------------------------------
!
!  plot colorbar for pregressive coloring
!
!----------------------------------------------------------------------

   IF (nset > 1)  CALL xcpalet(pcolbar)

  RETURN
END SUBROUTINE strm2d
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CTRSFC                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE ctrsfc(a,x,y,x1,x2,dx,y1,y2,dy,                              &
           nx,ibgn,iend,ny,jbgn,jend,                                   &
           label,time,runname,factor,tem1,tem2,tem3,                    &
           tem4,tem5,hterain,slicopt,pltopt)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    To plot a contour map for a 2-d surface array.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!    4/20/1994
!
!  MODIFICATION HISTORY:
!
!    9/27/95 (Yuhe Liu)
!    Fixed a bug in call of smth. Added the temporary array tem5 to
!    the argument list.
!
!    3/25/96 (Keith Brewster)
!    Added variables isize,jsize and replaced smth with smooth9p
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    a        2-d surface array.
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!
!    x1       value of x for first i grid point to plot
!    x2       value of x for last i grid point to plot
!    dx
!    y1       value of y for first j grid point to plot
!    y2       value of y for last j grid point to plot
!    dy
!
!    nx       first dimension of a
!    ibgn     index of first i grid point to plot
!    iend     index of last  i grid point to plot
!
!    ny       second dimension of a
!    jbgn     index of first j grid point to plot
!    jend     index of last  j grid point to plot
!
!    label    character string describing the contents of a
!
!    time     time of data in seconds
!
!    runname  character string decribing run
!
!    factor   scaling factor for data
!             contours are labelled a*factor
!    slicopt  slice orientation indicator
!       slicopt = 1, x-y slice of u,v at z index kslice is plotted.
!       slicopt = 2, x-z slice of u,w at y index jslice is plotted.
!       slicopt = 3, y-z slice of v,w at x index islice is plotted.
!       slicopt = 4, x-y slice of u,v at z index islice is plotted.
!       slicopt = 5, xy-z cross section of wind islice is plotted.
!       slicopt = 6, data field on constant p-level is plotted.
!       slicopt = 0, all of the three slices above are plotted.
!    plot     variable plot option (0/1/2/3)
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
!    tem4     Temporary work array.
!    tem5     Temporary work array.
!
!
!    hterain  The height of the terrain.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny

  REAL :: a(nx,ny)
  REAL :: x(nx,ny)
  REAL :: y(nx,ny)

  REAL :: x1,x2,dx,y1,y2,dy
  INTEGER :: ibgn,iend,jbgn,jend,length

  CHARACTER (LEN=9) :: timhms
  CHARACTER (LEN=*) :: label
  CHARACTER (LEN=*) :: runname

  REAL :: time
  REAL :: factor

  REAL :: tem1(*)
  REAL :: tem2(*)
  REAL :: tem3(*)
  REAL :: tem4(*)
  REAL :: tem5(*)

  REAL :: hterain(nx,ny)

  INTEGER :: slicopt
  INTEGER :: pltopt       ! variable plot option (0/1/2/3)

  INTEGER :: ovrtrn,trnplt               ! overlay terrain option (0/1)
  REAL :: trninc,trnmin, trnmax   ! terrain interval minimum, maximum
  REAL :: ztmin,ztmax
  COMMON /trnpar/ trnplt,ovrtrn,trninc,trnmin, trnmax,ztmin,ztmax

  INTEGER :: smooth
  COMMON /smoothopt/smooth

  INTEGER :: wrtflag
  CHARACTER (LEN=120) :: label_copy
  CHARACTER (LEN=90) :: levlab
  CHARACTER (LEN=50) :: timelab
  CHARACTER (LEN=25) :: timestring
  COMMON /timelev/wrtflag, timelab, levlab, timestring

  INTEGER :: xpbgn,xpend,ypbgn,ypend
  COMMON /processors/ xpbgn, xpend, ypbgn, ypend

!-----------------------------------------------------------------------
!
!  Common blocks for including iskip, jskip
!
!-----------------------------------------------------------------------
!
  REAL    :: xw1,xw2,yw1,yw2
  INTEGER :: iskip, jskip
  COMMON /pltwdw/ xw1,xw2,yw1,yw2, iskip, jskip

!----------------------------------------------------------------------
!
! Include files
!
!---------------------------------------------------------------------

  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,ij,isize,jsize,llabel
  CHARACTER (LEN=120) :: title

  INTEGER :: idsize, jdsize, mnsize
  INTEGER :: tinds, tind1,tind2,tind3,tind4,tind5,tind6
                ! temporary arrays index, assume size of tem5 > 6*nx*ny
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (slicopt == 2 .OR. slicopt == 3 .OR. slicopt == 5 .OR. slicopt >9) RETURN

  isize=(iend-ibgn)+1
  jsize=(jend-jbgn)+1

  idsize = isize            ! global maximum isize
  jdsize = jsize
  CALL mpmaxi(idsize)
  CALL mpmaxi(jdsize)

  mnsize = idsize*jdsize

  tind1 = 1             ! reuse a 3d temporary array 'tem5' as several 2D
  tind2 = tind1+mnsize  ! arrays inside ctr2d
  tind3 = tind2+mnsize
  tind4 = tind3+mnsize
  tind5 = tind4+mnsize
  tind6 = tind5+mnsize

!  tinds = SIZE(tem5)
!  IF (tinds < 6*mnsize) THEN
!    WRITE(6,*) 'ERROR: temporary array tem5 is too small.'
!    WRITE(6,*) '       Inside ctrsfc: isize = ',isize,' jsize = ',jsize, &
!               ' size(tem5) = ',tinds
!    CALL arpsstop('Temporary array too small inside ctrsfc.',1)
!  END IF

  label_copy = label
  llabel = 120
  CALL xstrlnth(label_copy, llabel)
  IF(myproc == 0)CALL xpscmnt('Start plotting '//label_copy(1:llabel))
!
!-----------------------------------------------------------------------
!
!  Set up terrain, if needed.
!
!-----------------------------------------------------------------------
!
  IF(ovrtrn == 1)  THEN
    DO j=jbgn,jend
      DO i=ibgn,iend
        ij = i-ibgn+1 + (j-jbgn)*isize
        tem4(ij)=hterain(i,j)
      END DO
    END DO
  END IF

  CALL get_forecast_hms( time, timhms )

  WRITE(timelab,'(a,F8.1,a)') 'T=',time,' s ('//TRIM(timhms)//')'
  CALL get_time_string ( time, timestring,'Z  ',0 )

  DO j=jbgn,jend
    DO i=ibgn,iend
      ij = i-ibgn+1 + (j-jbgn)*isize
      tem1(ij) = -9999.0
      IF(a(i,j) /= -9999.0) tem1(ij)=a(i,j)*factor
      tem2(ij)=x(i,j)
      tem3(ij)=y(i,j)
    END DO
  END DO

  levlab=' '
  WRITE(title,'(a)') label

  length = 120
  CALL strlnth( title, length)
  CALL strmin ( title, length)

  DO i=1,smooth
    CALL smooth9pmv(tem1,isize,jsize,1,isize,1,jsize,tem5)
  END DO

  CALL ctr2d(tem1,tem2,tem3, x1,x2,dx, y1,y2,dy,iskip+1,jskip+1,        &
             isize,jsize,title(1:length),runname,                       &
             tem4,slicopt,pltopt,mnsize,                                &
             tem5(tind1),tem5(tind2),tem5(tind3),                       &
             tem5(tind4),tem5(tind5),tem5(tind6))

  IF(myproc == 0) CALL xpscmnt('End plotting '//label_copy(1:llabel))

  RETURN
END SUBROUTINE ctrsfc
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE OVERLAY                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE overlay (layovr)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    Set the layover counter parameter in the laypar common block
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!
!  MODIFICATION HISTORY:
!    6/08/92  Added full documentation (K. Brewster)
!
!    8/08/93 (MX)
!    Automatically set the overlay parameter when input is not zero.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    layovr   The 'overlay' parameter.
!             If layover .ne. 0, the following 2-d contour plot will be
!             superimposed on the previous plot.
!             layover =1, 2, ... indicating this is the
!                layover'th (1st or 2nd ...) plot to be overlayed
!                on the previous one.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: layovr
!
!-----------------------------------------------------------------------
!
!  Plotting control common blocks
!
!-----------------------------------------------------------------------
!
  INTEGER :: layover, first_frame
  COMMON /laypar/ layover
  COMMON /frstfrm/ first_frame
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF( first_frame == 1 .OR. layovr == 0 ) THEN
    layover = 0
  ELSE
    layover = layover + 1
  END IF

  first_frame = 0

  RETURN
END SUBROUTINE overlay

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE STYXRT                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE styxrt( yxrt )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    Set the scaling factor of the y/x ratio of the plot.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!
!  MODIFICATION HISTORY:
!    6/08/92  Added full documentation (K. Brewster)
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    yxrt      Ratio of height to length of plot space
!              Default is set in the main program to 1.0
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  REAL :: yxrt
!
!-----------------------------------------------------------------------
!
!  Plotting control common blocks
!
!-----------------------------------------------------------------------
!
  REAL :: yxratio
  COMMON /yratio/ yxratio       ! the scaling factor the y/x ratio.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  yxratio = yxrt

  RETURN
END SUBROUTINE styxrt
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 FUNCTION XFINC                       ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

  REAL FUNCTION xfinc(x)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    Automatically divide domain (0,x) to a number of subdomain
!    with interval xfinc which is >=4 and =<16 for fold=1.0
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!    sometime
!
!  MODIFICATIONS:
!    6/09/92  Added full documentation (K. Brewster)
!
!-----------------------------------------------------------------------
!
!  INPUT:
!     x       not sure
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  REAL :: x
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: ipower
  REAL :: d,fold
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ipower= INT( ALOG10(x) )
  d= INT(x/(10.0**ipower))
  fold=1.0
  xfinc=0.1*x
  IF( d >= 0.0 .AND. d < 3.0 ) THEN
    xfinc=2.0*10.0**(ipower-1)
  ELSE IF( d >= 3.0 .AND. d < 7.0 ) THEN
    xfinc=5.0*10.0**(ipower-1)*fold
  ELSE IF( d >= 7.0 .AND. d < 10. ) THEN
    xfinc=1.0*10.0** ipower*fold
  END IF
  IF(xfinc == 0.0) xfinc=x*0.1
  RETURN
  END FUNCTION xfinc
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CLIPWD                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE clipwd(x1,y1,x2,y2,idispl)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Return the portion of a line that is within a given window
!  (xw1,xw2,yw1,yw2)
!
!  If the given line is completely outside the window,
!  idispl=0, otherwise, idispl=1.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  3/6/93
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    x1       value of x for first i grid point to plot
!    x2       value of x for last i grid point to plot
!    y1       value of y for first j grid point to plot
!    y2       value of y for last j grid point to plot
!
!    idispl   line orientation indicator
!       idispl = 0, the given line is completely outside the window
!       idispl = 1, the given line is partly inside the window
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
  REAL :: x1,x2,y1,y2
  INTEGER :: idispl
!
!-----------------------------------------------------------------------
!
!  Common blocks for plotting control parameters
!
!-----------------------------------------------------------------------
!
  REAL    :: xw1,xw2,yw1,yw2
  INTEGER :: iskip, jskip
  COMMON /pltwdw/ xw1,xw2,yw1,yw2, iskip, jskip
  INTEGER :: ic1(4),ic2(4)
!
!-----------------------------------------------------------------------
!
!  Misc. local Variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,knt,isw
  REAL :: x0,y0
  REAL :: isum1,isum2,ic01,ic02,ic03,ic04
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  knt = 0
  5     knt = knt+1
  CALL encodwd(x1,y1,ic1)
  CALL encodwd(x2,y2,ic2)

  isum1=ic1(1)+ic1(2)+ic1(3)+ic1(4)
  isum2=ic2(1)+ic2(2)+ic2(3)+ic2(4)
  idispl=1
  IF(isum1+isum2 == 0) GO TO 999

  idispl=0
  DO i=1,4
    IF(ic1(i)+ic2(i) == 2) GO TO 999
  END DO
!
!-----------------------------------------------------------------------
!
!  Make sure (x1,y1) is outside the window
!
!-----------------------------------------------------------------------
!
  isw=0
  IF(isum1 == 0) THEN
    ic01=ic1(1)
    ic02=ic1(2)
    ic03=ic1(3)
    ic04=ic1(4)
    DO i=1,4
      ic1(i)=ic2(i)
    END DO
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
  END IF

  IF(ic1(1) == 1) THEN
    y1=y1+(xw1-x1)*(y2-y1)/(x2-x1)
    x1=xw1
  ELSE IF(ic1(2) == 1) THEN
    y1=y1+(xw2-x1)*(y2-y1)/(x2-x1)
    x1=xw2
  ELSE IF(ic1(3) == 1) THEN
    x1=x1+(yw1-y1)*(x2-x1)/(y2-y1)
    y1=yw1
  ELSE IF(ic1(4) == 1) THEN
    x1=x1+(yw2-y1)*(x2-x1)/(y2-y1)
    y1=yw2
  END IF

  IF(isw == 1) THEN
    x0=x1
    y0=y1
    x1=x2
    y1=y2
    x2=x0
    y2=y0
  END IF

  idispl=1

  IF(knt > 10) THEN
    WRITE(6,*)'Dead loop encountered in CLIPWD, job stopped.'
    STOP 991
  END IF

  GO TO 5

  999   RETURN
END SUBROUTINE clipwd
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ENCODWD                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE encodwd(x,y,ic)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Encode a line section for window clipping purpose.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  3/6/93
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    x       value of x
!    y       value of y
!    ic
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
  REAL :: x,y
  INTEGER :: ic(4)
!
!-----------------------------------------------------------------------
!
!  Common blocks for plotting control parameters
!
!-----------------------------------------------------------------------
!
  REAL :: xw1,xw2,yw1,yw2
  INTEGER :: iskip, jskip
  COMMON /pltwdw/ xw1,xw2,yw1,yw2, iskip, jskip
!
!-----------------------------------------------------------------------
!
!  Misc. local Variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO i=1,4
    ic(i)=0
  END DO

  IF(x < xw1) ic(1)=1
  IF(x > xw2) ic(2)=1
  IF(y < yw1) ic(3)=1
  IF(y > yw2) ic(4)=1

  RETURN
END SUBROUTINE encodwd
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CTRCOL                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE ctrcol (icol,icol0)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    Set the color  for field to plotted by CTR2D.
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    icol    begin color index
!    icol0   end color index
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: icol,icol0
!
!-----------------------------------------------------------------------
!
!  Plotting control common blocks
!
!-----------------------------------------------------------------------
!
  INTEGER :: icolor,icolor1,lbcolor,trcolor
  COMMON /recolor/icolor,icolor1,lbcolor,trcolor

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  icolor  = icol
  icolor1 = icol0

  RETURN
END SUBROUTINE ctrcol


!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CTRVTR                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE ctrvtr (units0,type0)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    Set the units and type for plot wind
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    units0   the units of wind
!    type0    the type of wind
!    wcolor0  the color index
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: units0,type0
!
!-----------------------------------------------------------------------
!
!  Plotting control common blocks
!
!-----------------------------------------------------------------------
!
  INTEGER :: iunits, itype
  COMMON /windvtr/iunits, itype
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  iunits = units0
  itype = type0

  RETURN
END SUBROUTINE ctrvtr


!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE  VARPLT                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE varplt( var )

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    Set the variable plot name for xconta.
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!
!  MODIFICATION HISTORY:
!    3/28/96
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    var   variable plot name
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  CHARACTER (LEN=*)  :: var
  CHARACTER (LEN=12) :: varname

  COMMON /varplt1/ varname

  varname=var

  RETURN
END SUBROUTINE varplt

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE VTRSFC                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE vtrsfc(u,v, x,y, xw,xe,dx, ys,yn,dy,                         &
           nx,ibgn,iend,ist, ny,jbgn,jend,jst,                          &
           label,time, runname, factor, slicopt,                        &
           tem1,tem2,tem3,tem4,                                         &
           tem5,tem6,hterain)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    Plot vector fields in 2-d array
!
!  AUTHOR: Min Zou
!  4/28/97
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    u        2-dimensional array of u wind components (m/s)
!    v        2-dimensional array of v wind components (m/s)
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in physical space (m)
!
!    xw       value of x for first i grid point to plot
!    xe       value of x for last i grid point to plot
!    ys       value of y for first j grid point to plot
!    yn       value of y for last j grid point to plot
!
!    nx       first dimension of b
!    ibgn     index of first i grid point to plot
!    iend     index of last  i grid point to plot
!
!    ny       second dimension of b
!    jbgn     index of first j grid point to plot
!    jend     index of last  j grid point to plot
!
!
!    time     time of data in seconds
!
!    runname  character string decribing run
!
!    factor   scaling factor for winds
!             V*factor wind vectors are plotted
!    slicopt  slice orientation indicator
!       slicopt = 1, x-y slice of u,v at z index kslice is plotted.
!       slicopt = 2, x-z slice of u,w at y index jslice is plotted.
!       slicopt = 3, y-z slice of v,w at x index islice is plotted.
!       slicopt = 4, x-y slice of u,v at z index islice is plotted.
!       slicopt = 5, xy-z cross section of wind islice is plotted.
!       slicopt = 6, data field on constant p-level  is plotted.
!       slicopt = 0, all of the three slices above are plotted.
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
!    tem4     Temporary work array.
!    tem5     Temporary work array.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny

  REAL :: u(nx,ny)
  REAL :: v(nx,ny)

  REAL :: x(nx,ny)
  REAL :: y(nx,ny)

  CHARACTER (LEN=*) :: runname
  CHARACTER (LEN=*) :: label

  REAL :: xw,xe,dx,ys,yn,dy
  INTEGER :: ibgn,iend,ist, jbgn,jend,jst

  REAL :: time,factor
  INTEGER :: slicopt

  INTEGER :: iunits, itype
  COMMON /windvtr/iunits, itype

  CHARACTER (LEN=12) :: varname
  COMMON /varplt1/ varname

  REAL :: xw1,xe1,ys1,yn1
  COMMON /xuvpar/xw1,xe1,ys1,yn1
!
!-----------------------------------------------------------------------
!
!  Work arrays: tem1,tem2,tem3,tem4,tem5 of size at least
!          max( nx*ny, nx*nz, ny*nz).
!
!-----------------------------------------------------------------------
!
  REAL :: tem1(*),tem2(*),tem3(*),tem4(*),tem5(*)
  REAL :: tem6(*)
!
!-----------------------------------------------------------------------
!
!  Common blocks for plotting control parameters
!
!-----------------------------------------------------------------------
!
  REAL :: x01,y01                  ! the first  point of interpolation
  REAL :: x02,y02                  ! the second point of interpolation
  REAL :: zlevel                   ! the given height of the slice
  REAL :: sinaf,cosaf,dist,sqrtdxy
  COMMON /slicev/x01,y01,x02,y02,sinaf,cosaf,dist,sqrtdxy
  COMMON /sliceh/zlevel

  INTEGER :: ovrobs,obsset,obscol,obs_marktyp
  REAL :: obs_marksz,obs_valsz
  COMMON /obspar/ ovrobs,obsset,obscol,obs_marktyp,                    &
                  obs_marksz,obs_valsz
!
!-----------------------------------------------------------------------
!
!  Misc. local Variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,ij,istep,jstep,length,isize,jsize
  REAL :: uunit
  CHARACTER (LEN=9)   :: timhms
  CHARACTER (LEN=120) :: title

  INTEGER :: icolor,icolor1,lbcolor,trcolor        ! required color
  COMMON /recolor/icolor,icolor1,lbcolor,trcolor

  INTEGER :: trnplt                ! flag to plot terain (1 or 0)
  REAL :: hterain(nx,ny)           ! The height of the terrain.

  INTEGER :: ovrtrn         ! overlay terrain option (0/1)

  REAL :: trninc,trnmin, trnmax    ! terrain interval minimum, maximum
  REAL :: ztmin,ztmax
  COMMON /trnpar/ trnplt,ovrtrn,trninc,trnmin, trnmax,ztmin,ztmax

  INTEGER :: smooth
  COMMON /smoothopt/smooth

  INTEGER :: wrtflag, llabel
  CHARACTER (LEN=90) :: levlab
  CHARACTER (LEN=50) :: timelab
  CHARACTER (LEN=25) :: timestring
  COMMON /timelev/wrtflag,timelab, levlab, timestring
  CHARACTER (LEN=120) :: label_copy

  INTEGER :: xpbgn,xpend,ypbgn,ypend
  COMMON /processors/ xpbgn, xpend, ypbgn, ypend

  INTEGER :: idsize, jdsize, mnsize
  INTEGER :: tinds, tind1,tind2,tind3,tind4,tind5,tind6,tind7,tind8

!----------------------------------------------------------------------
!
! Include files
!
!---------------------------------------------------------------------
  INCLUDE 'mp.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (slicopt == 2 .OR. slicopt == 3 .OR. slicopt == 5 .OR. slicopt >9) RETURN

  isize=(iend-ibgn)+1
  jsize=(jend-jbgn)+1

  idsize = isize            ! global maximum isize
  jdsize = jsize
  CALL mpmaxi(idsize)
  CALL mpmaxi(jdsize)

  mnsize = idsize*jdsize

  tind1 = 1             ! reuse a 3d temporary array 'tem6' as several 2D
  tind2 = tind1+mnsize  ! arrays inside vtr2d
  tind3 = tind2+mnsize
  tind4 = tind3+mnsize
  tind5 = tind4+mnsize
  tind6 = tind5+mnsize
  tind7 = tind6+mnsize
  tind8 = tind7+mnsize

!  tinds = SIZE(tem6)
!  IF (tinds < 5*mnsize) THEN
!    WRITE(6,*) 'ERROR: temporary array tem6 is too small.'
!    CALL arpsstop('Temporary array too small inside vtrsfc.',1)
!  END IF

!
!-----------------------------------------------------------------------
!
!  Set up terrain, if needed.
!
!-----------------------------------------------------------------------
!
  label_copy = label
  llabel = 120
  CALL xstrlnth(label_copy, llabel)

  IF(myproc ==0) CALL xpscmnt('Start plotting '//label_copy(1:llabel))

  IF(trnplt == 1 .OR.trnplt == 2 .OR. ovrtrn == 1)  THEN
    DO j=jbgn,jend
      DO i=ibgn,iend
        ij = i-ibgn+1 + (j-jbgn)*isize
        tem5(ij)=hterain(i,j)
      END DO
    END DO
  END IF

  CALL get_forecast_hms( time, timhms )
  WRITE(timelab,'(a,F8.1,a)') 'T=',time,' s ('//TRIM(timhms)//')'

  CALL get_time_string ( time, timestring,'Z  ',0 )

!   length=50
!   CALL strmin(timelab,length)
!   write(timelab,'(a,'' '',a)') timestring(1:21), timelab(1:length)
!   print*,'in vtrsfc', timelab


  DO j=jbgn,jend
    DO i=ibgn,iend
      ij = i-ibgn+1 + (j-jbgn)*isize
      tem1(ij) = -9999.0
      tem2(ij) = -9999.0
      IF(u(i,j) /= -9999.0) tem1(ij)=u(i,j)*factor
      IF(v(i,j) /= -9999.0) tem2(ij)=v(i,j)*factor
      tem3(ij)=x(i,j)
      tem4(ij)=y(i,j)
    END DO
  END DO

  levlab = 'First level above ground (surface)'
  WRITE(title,'(2A)') 'U-V ',label

!  length = 120
!  CALL strlnth( title, length )
!  CALL strmin ( title, length)
  length = LEN_TRIM(title)

  uunit = 10.0
  CALL xvmode(1)
  istep = ist
  jstep = jst

  DO i=1,smooth
    CALL smooth9pmv(tem1,isize,jsize,1,isize,1,jsize,tem6)
    CALL smooth9pmv(tem2,isize,jsize,1,isize,1,jsize,tem6)
  END DO

  CALL vtr2d(tem1,tem2,tem3,tem4, uunit, xw,xe,dx,ys,yn,dy,             &
             isize,istep,jsize,jstep,title(1:length),runname, 1,        &
             tem5,slicopt,mnsize,tem6(tind1),tem6(tind2),tem6(tind3),   &
             tem6(tind4),tem6(tind5),tem6(tind6),tem6(tind7),tem6(tind8))

  IF(myproc==0) CALL xpscmnt('End plotting '//label_copy(1:llabel))

  RETURN
END SUBROUTINE vtrsfc
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SET_INTERVAL               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE set_interval(zmin1,zmax1,ctmin,ctmax,cl)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    Limited contour interval when uinc = -9999.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    zmin1    the minimum value of 2-D array
!    zmax1    The maximum value of 2-D array
!    ctmin    the input minimum value
!    ctmax    the input maximum value
!    cl       the intervals
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  REAL :: cl(*)
  REAL :: ctmin, ctmax
  REAL :: zmin1, zmax1

  REAL :: zmin, zmax
  REAL :: zinc

  INTEGER :: nmin, nmax
  COMMON /xclm19/ nmin, nmax

  REAL    :: clref, lcptn
  INTEGER :: labtyp, iclf, lhilit, ihlf, kct0
  COMMON /xcrf17/clref,lcptn,labtyp,iclf,lhilit,ihlf,kct0

  INTEGER :: nhole, nvtrbadv
  REAL    :: specia
  COMMON /zchole/ nhole,specia,nvtrbadv

  INTEGER :: nch
  COMMON /xoutch/ nch

  REAL    :: diff, kzinc, clv, eps
  INTEGER :: ncmin, ncmax, ncnt, kcount

  REAL :: xfinc
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF(ctmin == 0.0 .AND. ctmax == 0.0) THEN
    cl(2)=cl(1)+ xfinc(zmax1-zmin1)/2
  ELSE
    cl(2)=cl(1)+ xfinc(ctmax-ctmin)/2
  END IF
  IF(cl(2)-cl(1) == 0.0) cl(2)=cl(1)+1.0
  zinc = cl(2)-cl(1)

  ncmin=nmin
  ncmax=nmax
  diff=ctmax-ctmin
  IF( diff - ABS(zinc)*1.0E-6  > 0.0) THEN
    GOTO 4
  END IF
  WRITE(nch,'(a,a)')                                               &
    ' Bad first guess of contour increment or field is constant',  &
    ', number of contours is one.'
  ncnt=1
  cl(1)= ctmin
  RETURN

  4 kcount=0
  1 CONTINUE
  eps=0.1*zinc
  kcount=kcount+1
  IF( kcount > 20) GO TO 998
  kzinc=(ctmin-clref)/zinc
  zmin=kzinc*zinc+clref
  kzinc=(ctmax-clref)/zinc
  zmax=kzinc*zinc+clref
  IF(ctmin-clref > 0.0) zmin=zmin+zinc
  IF(ctmax-clref < 0.0) zmax=zmax-zinc
!
  clv=zmin-zinc
  ncnt=0
  6    clv=clv+zinc
  IF(clv-zmax-eps > 0.0) GO TO     8
  ncnt=ncnt+1
  IF(ncnt > ncmax) THEN
    zinc=zinc*2
    WRITE(nch,1000) ncmax, zinc
    1000 FORMAT(' Number of contours > ',i3,' ,Zinc is doubled. Zinc='  &
            ,e10.3)
    GO TO 1
  END IF
  IF( ABS( clv-clref ) < eps ) clv=clref
  cl(ncnt)=clv
  GO TO 6
  8    CONTINUE

  IF( ncnt < ncmin) THEN
    zinc=zinc/2
    WRITE(nch,2000) ncmin,zinc
    2000 FORMAT(' Number of contours < ',i3,' ,Zinc is halved. Zinc='   &
           ,e10.3)
    GO TO 1
  END IF
  WRITE(nch,'('' * NUMBER OF CONTOURS= '',I5,''  Min='',E12.4,          &
  &   '' Max='', e12.4,'' Inc='',e12.5 )')                              &
      ncnt,ctmin,ctmax,zinc

  IF( zmin1 >= ctmin .AND. zmax1 <= ctmax) THEN
    zinc = cl(2) - cl(1)
    WRITE(nch,'(''SET MINIMUM CONTOUR INTERVAL IS'',E12.4,              &
    &   '' ctmin='',e12.4,'' ctmax='',e12.4 )')zinc,ctmin,ctmax
    CALL xctref(zinc)
    CALL xnctrs( 1,900)
  ELSE
    WRITE(nch,'(''NO NEED SET MINIMUM CONTOUR INTERVAL'')' )
    WRITE(nch,'(''CNTOUR INTERVAL IS SET AUTOMATICALLY'')' )
    cl(2)=cl(1)+ xfinc(zmax1-zmin1)/2
    IF(cl(2)-cl(1) == 0.0) cl(2)=cl(1)+1.0
  END IF
  RETURN
  998  WRITE(nch,*)' Contour levels can not be selected by XCNTLV.'
  WRITE(nch,*)                                                          &
      ' Plz alter input contour interval or limits of contour number'

END SUBROUTINE set_interval

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE XRCH1                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE xrch1( r,ch,lch)

! Return real number R as a character string in automatically set format
  REAL :: r
  CHARACTER (LEN=20) :: str
  CHARACTER (LEN=*) :: ch

  CALL get_format(r,str)
  IF(ABS(r-0.0) < 1.e-20) THEN
    WRITE(ch,'(F3.1)') r
  ELSE
    WRITE(ch,str) r
  END IF
  lch=20
  CALL strlnth( ch, lch)
  CALL strmin ( ch, lch)
  RETURN
END SUBROUTINE xrch1

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE  GET_FORMAT                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE get_format(r,ch)
  INTEGER :: npoz
  CHARACTER (LEN=1) :: FORM,ndrob
  CHARACTER (LEN=20) :: ch
  WRITE(ch,10)r
  10   FORMAT(g11.4)
  DO i=20,1,-1
    IF(ch(i:i) == '0'.OR.ch(i:i) == ' ') THEN
      ch(i:i)=' '
    ELSE
      GO TO 1
    END IF
  END DO
  1    CONTINUE
  npoz=0
  ndot=0
  nmant=0
  ndrob=' '
  FORM='F'
  DO i = 1,20
    IF(ch(i:i) /= ' ' ) npoz=npoz+1
    IF(ch(i:i) == 'E') FORM='E'
    IF(ndrob == '.'.AND.ch(i:i) /= ' ') ndot=ndot+1
    IF(ch(i:i) == '.') ndrob='.'
    IF(FORM /= 'E') nmant=npoz
  END DO
  npoz=npoz
  IF(FORM == 'F') THEN
    IF(ndot /= 0) THEN
      WRITE(ch,20) '(',FORM,npoz,'.',ndot,')'
    ELSE
      WRITE(ch,20) '(',FORM,npoz,'.',ndot,')'
    END IF
  ELSE IF(FORM == 'E') THEN
    ch = '(1PE20.2)'
  ELSE
    WRITE(ch,20) '(',FORM,npoz,'.',nmant,')'
  END IF
  20   FORMAT(a1,a1,i1,a1,i1,a1)
  RETURN
END SUBROUTINE get_format

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE DRAWMAP                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE drawmap(nunit)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  This subroutine will plot the map
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!
!  MODIFICATION HISTORY:
!    6/2/97 Min Zou
!    Read multiple mapfile only once. Using differnt line style to
!    plot mapdata.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nunit     the channel of the mapfile data
!    mapfile   character of map file name
!
!-----------------------------------------------------------------------
!
!  Variable Declarations
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INCLUDE 'arpsplt.inc'

  INTEGER :: nunit,i
  INTEGER :: lmapfile

  CHARACTER (LEN=256) :: mapfile(maxmap)
  INTEGER :: mapgrid,mapgridcol, kolor
  REAL :: latgrid,longrid
  INTEGER :: nmapfile,mapcol(maxmap),mapline_style(maxmap)
  COMMON /mappar1/nmapfile,mapcol,mapline_style,mapfile
  COMMON /mappar2/mapgrid,mapgridcol,latgrid,longrid

  REAL :: x1,x2,y1,y2

  CALL xpscmnt('Start of map plotting ')

  CALL xqmap (x1,x2,y1,y2)
  CALL xwindw(x1,x2,y1,y2)
  CALL xqcolor(kolor)

  DO i=1,nmapfile
    CALL xcolor(mapcol(i))
    IF(mapline_style(i) == 1) THEN
      CALL xthick(1)
      CALL xbrokn(6,3,6,3)
    ELSE IF(mapline_style(i) == 2) THEN
      CALL xthick(1)
    ELSE IF(mapline_style(i) == 3) THEN
      CALL xthick(3)
      CALL xfull
    END IF

    lmapfile=256
    CALL xstrlnth(mapfile(i), lmapfile)

    CALL xdrawmap_new(nunit,mapfile(i)(1:lmapfile),latgrid,longrid,     &
                      i,mapgridcol)
  END DO

  CALL xcolor(kolor)
  CALL xfull
  CALL xwdwof

  CALL xpscmnt('End of map plotting ')

  RETURN
END SUBROUTINE drawmap
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE PLTOBS                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE pltobs(obopt)
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Plots observations on an arpsplt contour map.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!  obopt    Plotting option
!           1  Plot data in obs1 as characters
!           2  Plot data in obs1 and obs2 as characters
!           3  Plot wind arrows with obs1 as u and obs2 as v.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!

  IMPLICIT NONE

  INCLUDE 'arpsplt.inc'
!
!  Arguments
!
  INTEGER :: obopt
!
  INTEGER :: nobs
  REAL :: latob(mxsfcob)
  REAL :: lonob(mxsfcob)
  REAL :: obs1(mxsfcob)
  REAL :: obs2(mxsfcob)
!
  COMMON /sfc_obs1/ nobs
  COMMON /sfc_obs2/ latob,lonob,obs1,obs2
!
!  Plotting parameters
!
  CHARACTER (LEN=1) :: cross
  PARAMETER(cross='+')
!
  INTEGER :: ovrobs,obsset,obscol,obs_marktyp
  REAL :: obs_marksz,obs_valsz
  COMMON /obspar/ ovrobs,obsset,obscol,obs_marktyp,                    &
                  obs_marksz,obs_valsz
!
!  Plotting common blocks
!
  INTEGER :: icolor,icolor1,lbcolor,trcolor       ! required color
  COMMON /recolor/icolor,icolor1,lbcolor,trcolor

  REAL :: ctinc,ctmin,ctmax,vtunt  ! contour interval and vector unit
  COMMON /incunt/ ctinc,ctmin,ctmax,vtunt

  REAL :: xleng,vunit
  COMMON /vecscl/ xleng,vunit

  INTEGER :: iunits, itype
  COMMON /windvtr/iunits, itype
!
!  Misc local variables
!
  INTEGER :: iob
  REAL :: orgmag,yoff,yoff2
  REAL :: x1,x2,y1,y2
  REAL :: xob,yob
  CHARACTER (LEN=4) :: chplot
  INTEGER :: imkrfil,intobs
!
!  Set-up plotting space and zxplot variables
!
  CALL xcolor(lbcolor)
  CALL xqmap (x1,x2,y1,y2)
  CALL xwindw(x1,x2,y1,y2)
  CALL xchori(0.0)
  CALL xqchsz(orgmag)

  yoff=0.5*orgmag
  yoff2=2.*yoff
  CALL xchmag(obs_valsz)

  IF(obopt == 1) THEN
    CALL xcolor(obscol)
    CALL xmrksz(obs_marksz)
    DO iob=1,nobs
      IF(obs1(iob) > -98. .AND. obs1(iob) < 500.) THEN
        CALL xlltoxy(1,1,latob(iob),lonob(iob),xob,yob)
        intobs=max(min(nint(obs1(iob)),999),-999)
        WRITE(chplot,'(i4)') intobs
        CALL xcharc((0.001*xob),(0.001*yob+yoff),chplot)
        CALL xmarker((0.001*xob),(0.001*yob),obs_marktyp)
      END IF
    END DO
  ELSE IF(obopt == 2) THEN
    CALL xcolor(obscol)
    CALL xmrksz(obs_marksz)
    DO iob=1,nobs
      CALL xlltoxy(1,1,latob(iob),lonob(iob),xob,yob)
      IF(obs1(iob) > -98. .AND. obs1(iob) < 500.) THEN
        intobs=max(min(nint(obs1(iob)),999),-999)
        WRITE(chplot,'(i4)') intobs
        CALL xcharc((0.001*xob),(0.001*yob+yoff),chplot)
        CALL xmarker((0.001*xob),(0.001*yob),MOD(obs_marktyp,5))
      END IF
      IF(obs2(iob) > -98. .AND. obs2(iob) < 500.) THEN
        intobs=max(min(nint(obs1(iob)),999),-999)
        WRITE(chplot,'(i4)') intobs
        CALL xcharc((0.001*xob),(0.001*yob+yoff),chplot)
        CALL xmarker((0.001*xob),(0.001*yob),obs_marktyp)
      END IF
    END DO
  ELSE IF(obopt == 3) THEN
    CALL xcolor(obscol)
    CALL xmrksz(obs_marksz)
    DO iob=1,nobs
      IF(obs1(iob) > -98. .AND. obs1(iob) < 500. .AND.                  &
         obs2(iob) > -98. .AND. obs2(iob) < 500.) THEN
        CALL xlltoxy(1,1,latob(iob),lonob(iob),xob,yob)
        xob=0.001*xob
        yob=0.001*yob
        IF( xob > x1 .AND. xob < x2 .AND. yob > y1 .AND. yob < y2 ) THEN
          CALL xarrow(obs1(iob),obs2(iob),xob,yob,xleng,vunit)
          !CALL XBARB(U,V,X0,Y0,wunits,XLENG,barbopt)
          CALL xmarker(xob,yob,obs_marktyp)
        END IF
      END IF
    END DO
  END IF

  CALL xcolor(lbcolor)
  CALL xchsiz(orgmag)
  CALL xfull
  CALL xwdwof

  RETURN
END SUBROUTINE pltobs
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE PLTSTA                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE pltsta(a,b,x,y,m,n,flag,slicopt)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  This subroutine will plot some station information.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!    Min Zou (6/1/97)
!
!  Modification history:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    a, b    2-dimension array of Variable
!    m, n    Array dimensions
!    x, y    x-coord and y-coord of the staions
!    flag    a flag for different plot
!    slicopt  slice orientation indicator
!
!-----------------------------------------------------------------------
!
!  Variable Declarations
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INCLUDE 'arpsplt.inc'

  INTEGER :: m,n
  REAL :: a(m,n)
  REAL :: b(m,n)
  REAL :: x(m,n)
  REAL :: y(m,n)
  INTEGER :: nsta,nstapro(mxstalo),nstatyp(mxstalo)
  REAL :: latsta(mxstalo), lonsta(mxstalo)
  CHARACTER (LEN=5) :: s_name(mxstalo)
  INTEGER :: ovrstaopt
  INTEGER :: ovrstam,staset,ovrstan,ovrstav,stacol,markprio,wrtstax
  INTEGER :: nsta_typ,sta_typ(30),sta_marktyp(30), sta_markcol(30)
  REAL :: sta_marksz(30)
  REAL :: wrtstad
  CHARACTER (LEN=256) :: stalofl
  COMMON /sta_par/ ovrstaopt,ovrstam,staset,ovrstan,ovrstav,stacol,     &
         markprio, nsta_typ,sta_typ,sta_marktyp,                        &
         sta_markcol,sta_marksz,stalofl,wrtstax,wrtstad
  COMMON /sta_loc/latsta,lonsta,nstatyp,nstapro,nsta
  COMMON /sta_loc1/s_name
  REAL :: xob(mxstalo), yob(mxstalo),aob(mxstalo),bob(mxstalo)
  COMMON /xob_yob/xob, yob
  INTEGER :: LEN,i,j

  REAL :: x01,x02,y01,y02
  REAL :: sinaf,cosaf,dist,sqrtdxy
  COMMON /slicev/x01,y01,x02,y02,sinaf,cosaf,dist,sqrtdxy

  INTEGER :: icolor,icolor1,lbcolor,trcolor       ! required color
  COMMON /recolor/icolor,icolor1,lbcolor,trcolor

  INTEGER :: layover
  COMMON /laypar/ layover

  CHARACTER (LEN=12) :: varname
  COMMON /varplt1/ varname

  REAL :: xori1,xori2,yori1,yori2,zori1,zori2
  COMMON /tmphc1/xori1,xori2,yori1,yori2,zori1,zori2

  REAL :: xleng,vunit
  COMMON /vecscl/ xleng,vunit

  INTEGER :: iunits, itype
  COMMON /windvtr/iunits, itype

  REAL :: x_tmp
  COMMON /tmphc2/ x_tmp

  REAL :: x1,x2,y1,y2
  REAL :: orgmag,yoff,yoff2,xoff, xoff2
  CHARACTER (LEN=30) :: ctmp

  INTEGER :: flag,slicopt,fg
  REAL :: xdist, ydist,xd0,yd0,xa,xb
  SAVE fg

  REAL :: xleng0, spd, dir, istand
  INTEGER :: iunits0, imkrfil
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! calculate xob and yob on the xy plane
  CALL xwindw(xori1,xori2,yori1,yori2)
  IF(fg == 0) THEN
    CALL xlltoxy(nsta,1,latsta,lonsta,xob,yob)
    DO i= 1,nsta
      xob(i) = xob(i)*0.001
      yob(i) = yob(i)*0.001
    END DO
    fg=1
  END IF

  CALL xqmap (x1,x2,y1,y2)
  CALL xwindw(x1,x2,y1,y2)
  CALL xchori(0.0)
  CALL xqchsz(orgmag)

  yoff=0.5*orgmag
  yoff2=3.*yoff
  xoff = 0.001*(x2-x1)
  xoff2 = 4.*xoff
  CALL xchsiz(0.8*orgmag)

  IF(ovrstav == 1 ) THEN  ! interpolation
    CALL intepo (nsta,xob,yob,aob,m,n,x,y,a)
    IF(varname(1:6) == 'vtrplt' .OR. varname(1:6) == 'vtpplt')          &
        CALL intepo (nsta,xob,yob,bob,m,n,x,y,b)
  END IF

  IF (flag == 1) THEN
    CALL xcolor(stacol)
    IF(slicopt == 5) THEN
      xa = (y02-y01)/(x02-x01)
      xb = y01 - xa*x01
    END IF
    DO i = 1,nsta
      IF(nstapro(i) <= markprio) THEN
        IF(wrtstax == 1 )THEN
          CALL xwindw(x1,x2-0.005*(x2-x1), y1-0.3*(y2-y1),y1)
          LEN=5
          CALL strlnth(s_name(i),LEN)
          CALL xchori(90.)
          IF(slicopt == 2 .OR. slicopt == 10) THEN
            CALL xwindw(xori1,xori2-0.005*(xori2-xori1),                &
                        y1-0.3*(y2-y1),y1)
            IF ( (xob(i) <= xori2.AND.xob(i) >= xori1) .AND.            &
                 (yob(i) <= yori2.AND.yob(i) >= yori1) ) THEN
              IF( ABS(yob(i)-x_tmp) <= wrtstad ) THEN
!               CALL XCHARR((xob(i)),y1-1.75*yoff2,
                CALL xcharr((xob(i)),y1-1.50*yoff2,                     &
                             s_name(i)(1:LEN))
              END IF
            END IF
          ELSE IF(slicopt == 3 .OR. slicopt == 11) THEN
            CALL xwindw(yori1,yori2-0.005*(yori2-yori1),                &
                    y1-0.3*(y2-y1),y1)
            IF ( (xob(i) <= xori2.AND.xob(i) >= xori1)                  &
                  .AND. (yob(i) <= yori2.AND.yob(i) >= yori1) ) THEN
              IF( ABS(xob(i)-x_tmp) <= wrtstad ) THEN
                CALL xcharr((yob(i)),y1-1.75*yoff2,                     &
                    s_name(i)(1:LEN))
              END IF
            END IF
          ELSE IF( slicopt == 5) THEN
            CALL xwindw(x1,x2-0.005*(x2-x1),                            &
                    y1-0.3*(y2-y1),y1)
            IF ( (xob(i) <= xori2.AND.xob(i) >= xori1)                  &
                  .AND. (yob(i) <= yori2.AND.yob(i) >= yori1) ) THEN
              xd0 = 1./(xa*xa+1.0)*((yob(i)-xb)*xa+xob(i))
              yd0 = xa*xd0+xb
              xdist = SQRT((x01-xd0)*(x01-xd0) + (y01-yd0)*(y01-yd0))
              ydist = SQRT((xob(i)-xd0)*(xob(i)-xd0)+                   &
                         (yob(i)-yd0)*(yob(i)-yd0))
              xdist = xdist+x1
              IF(ABS(ydist) <= wrtstad ) THEN
                CALL xcharr(xdist,y1-1.75*yoff2,                        &
                    s_name(i)(1:LEN))
              END IF
            END IF
          END IF
          CALL xchori(0.)
        END IF
      END IF
    END DO

  ELSE IF(flag == 0) THEN
    CALL xwindw(x1,x2,y1,y2)
    DO i = 1,nsta
      IF( (xob(i) >= x1.AND.xob(i) <= x2) .AND.                         &
            (yob(i) >= y1.AND.yob(i) <= y2) ) THEN
        IF(nstapro(i) <= markprio) THEN
          IF(ovrstan == 1) THEN
            LEN=5
            CALL strlnth(s_name(i),LEN)
            CALL xcharc((xob(i)),(yob(i)-yoff2),                        &
                  s_name(i)(1:LEN))
          END IF
          IF(ovrstam == 1) THEN
            DO j=1,nsta_typ
              CALL xmrksz(sta_marksz(j))
              CALL xcolor(sta_markcol(j))
              IF(nstatyp(i) == sta_typ(j)) THEN
                CALL xmarker((xob(i)),(yob(i)),                         &
                     sta_marktyp(j))
                IF(ovrstan == 1) THEN
                  LEN=5
                  CALL strlnth(s_name(i),LEN)
                  CALL xcharc((xob(i)),(yob(i)-yoff2),                  &
                      s_name(i)(1:LEN))
                END IF
              END IF
            END DO
          END IF
          IF(ovrstav == 1) THEN
            CALL xcolor(stacol)
            IF(varname(1:6) == 'vtrplt' .OR. varname(1:6) == 'vtpplt') THEN
!          IF(i.eq.1) THEN
              xleng0=xleng*0.0004
              IF(iunits == 2 ) THEN
                iunits0=2
                istand = 10.
                WRITE(ctmp,'(a30)')'10 knots'
              ELSE IF (iunits == 3) THEN
                iunits0=2
                istand = 10.
                WRITE(ctmp,'(a30)')'10 MPH'
              ELSE IF(iunits == 1) THEN
                iunits0=1
                istand = 5.
                WRITE(ctmp,'(a30)')'5 m/s'
              END IF
!          ENDIF
              IF(aob(i) /= -9999. .AND. bob(i) /= -9999.) THEN
                spd = SQRT(aob(i)*aob(i)+bob(i)*bob(i))
                dir = ATAN2(-1.*aob(i),-1.*bob(i))*180./3.1415926
                IF(dir <= 0.) dir = 360.+dir
!            CALL barb((xob(i)),
!    :                   (yob(i)),dir,spd,iunits0-1, xleng0)
                CALL xbarb(aob(i),bob(i),xob(i),yob(i),                 &
                     iunits0,xleng*0.65,2)

              END IF
            ELSE
              IF(aob(i) /= -9999.) THEN
                CALL xrch(aob(i),ctmp,LEN)
                IF(layover == 0) CALL xcharr((xob(i)-xoff2),            &
                    (yob(i)+yoff),ctmp(1:LEN))
                IF(layover == 1) CALL xcharl((xob(i)+xoff2),            &
                    (yob(i)+yoff) ,ctmp(1:LEN))
                IF(layover == 2) CALL xcharc((xob(i)+xoff2),            &
                    (yob(i)-yoff),ctmp(1:LEN))
                IF(layover == 3) CALL xcharl((xob(i))+xoff2,            &
                    (yob(i)-yoff) ,ctmp(1:LEN))
              END IF
            END IF
          END IF

        END IF
      END IF
    END DO
  END IF

  CALL xcolor(lbcolor)
  CALL xchsiz(orgmag)
  CALL xfull
  CALL xwdwof

  RETURN
END SUBROUTINE pltsta
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE RUNLAB                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE runlab(runname)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    Plot a run label at the lower left cornor of the picture frame.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    runname      character string of run label
!
!-----------------------------------------------------------------------
!
!  Variables Declarations
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  CHARACTER (LEN=*) :: runname
  REAL :: xl, xr, yb, yt, rotang, xlimit, ylimit
  INTEGER :: nopic, nxpic, nypic

  REAL :: ytmp, hch
  REAL :: pl, pr, pb, pt, px, py
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  CALL xqmap(xl,xr,yb,yt)

  CALL xqpspc( pl, pr, pb, pt)
  px = pr - pl
  py = pt - pb

  CALL xqnpic(nopic)
  CALL xqspac(nxpic, nypic, rotang, xlimit, ylimit)

  IF( rotang == 0.0 ) THEN

    IF(nopic == nxpic*nypic -(nxpic-1)) THEN
      CALL xcharl( xl, yb-0.15*(yt-yb)*px/py, runname )
    END IF

  ELSE

    IF(nopic == nypic*nxpic -(nypic-1)) THEN
      CALL xcharl( xl, yb-0.15*(yt-yb)*px/py, runname )
    END IF

  END IF

  RETURN
END SUBROUTINE runlab
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE VPROFIL                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE vprofil(nx,ny,nz,nzprofc,var,xc,yc,zpc,plwr,pupr,            &
           xpnt,ypnt,npoints,zlwr,zupr,xcaptn,ycaptn,npicprof,          &
           profil,height)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  This subroutine will plot the vertical profiles of a given
!  variable through points (xpnt(i),ypnt(i),i=1,npoints).
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!    Adwait Sathye (2/28/94)
!
!  Modification history:
!
!    4/18/94, (Ming Xue)
!    Major overhaul. Many temporary arrays removed. New frame option
!    added.
!
!    9/18/1995 (Ming Xue)
!    Fixed a problem in the code that determines kbgn and kend.
!
!    10/8/1996 (Y. Richardson)
!    Corrected a bug in the interpolation weights.
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx, ny, nz    Array dimensions
!    nzprofc       the maximum vertical index in height (zpc/zpsoilc)
!                  and variables to be profiled when calling vprofil
!                  subroutine. 06/10/2002, Zuwen He
!                  In the atmosphere model, the vertical index is
!                  typically nz-1, while in the soil model, it's nzsoil.
!    var           Variable data array
!    xc,yc,zpc     The coordinate of input data var.
!    plwr,pupr     Lower and upper bounds for the horiz. axis of profile
!    xpnt, ypnt    Arrays containing the X and Y locations of the
!                  mulitple profiles to be plotted
!    npoints       Number of profile points to be plotted
!    zlwr, zupr    Bounds in the vertical direction
!    xcaptn        Caption for the X axis
!    ycaptn        Caption for the Y axis
!
!  Work arrays:
!
!    profil,height Temporary arrays
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Variables passed in
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx, ny, nz, nzprofc

  REAL :: var(nx,ny,nz)
  REAL :: xc(nx,ny,nz), yc(nx,ny,nz), zpc(nx,ny,nz)

  REAL :: plwr,pupr
  REAL :: zlwr, zupr


  INTEGER :: npoints
  REAL :: xpnt(npoints), ypnt(npoints)

  CHARACTER (LEN=*) :: xcaptn
  CHARACTER (LEN=*) :: ycaptn

  INTEGER :: npicprof

  REAL :: profil(nz,npoints)
  REAL :: height(nz,npoints)

  LOGICAL :: multiprof
!
!-----------------------------------------------------------------------
!
!  Temporary local variables
!
!-----------------------------------------------------------------------
!
  REAL :: lower, upper, zmin, zmax
  REAL :: x1, x2, y1, y2
  REAL :: a(2,2)
  INTEGER :: i, j, k, ix, jy,kbgn,kend,ip,lchar
  REAL :: dx,dy,temp,hmaxk,hmink
  CHARACTER (LEN=90) :: ch

  REAL :: lblmag, ctrlbsiz, axlbsiz
  COMMON /labmag/ lblmag, ctrlbsiz, axlbsiz

  INCLUDE 'mp.inc'

  INTEGER :: nxlg, nylg
  INTEGER :: source, itags, itagr
  INTEGER, PARAMETER :: destination = 0
  INTEGER :: indxx,indxy,xp(2),yp(2)
  INTEGER :: ii,jj, ierr

  REAL :: xtem, ytem
  REAL, ALLOCATABLE :: vartem(:), zptem(:)
  REAL, ALLOCATABLE :: varctem(:,:,:), zpctem(:,:,:)
  REAL :: pl, pr, pb, pt, px, py

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  CALL xqpspc( pl, pr, pb, pt)
  px = pr - pl
  py = pt - pb

  nxlg = (nx-3)*nproc_x + 3
  nylg = (ny-3)*nproc_y + 3

  multiprof = .false.
  IF( npicprof == 0 .AND. npoints > 1 ) multiprof = .true.
!
!-----------------------------------------------------------------------
!
!  If lower boundary is bigger, then swap boundaries.
!
!-----------------------------------------------------------------------
!
  IF (plwr > pupr) THEN
    lower = pupr
    upper = plwr
  ELSE
    lower = plwr
    upper = pupr
  END IF
!
!-----------------------------------------------------------------------
!
!  Find corresponding coordinates for the boundaries in the Z
!  dimension. IF both have been set to 0, use the boundary values.
!  Else, loop through the zpc array and find the location of the
!  point.
!
!-----------------------------------------------------------------------
!
  dx = xc(2,1,1)-xc(1,1,1)
  dy = yc(1,2,1)-yc(1,1,1)

  zmin = zpc(1,1,1)
  zmax = zpc(1,1,nzprofc)
  DO j=1,ny-1
    DO i=1,nx-1
      zmin=MIN( zmin, zpc(i,j,1))
      zmax=MAX( zmax, zpc(i,j,nzprofc))
      zmin=MIN( zmin, zpc(i,j,nzprofc))  ! because of soil model, Zuwen
      zmax=MAX( zmax, zpc(i,j,1))        ! because of soil model, Zuwen
    END DO
  END DO

  CALL mpmax0(zmax,zmin)

  IF( zlwr /= zupr ) THEN
    zmin = MAX(zmin, zlwr)
    zmax = MIN(zmax, zupr)
  END IF

  IF( zmax < zmin) WRITE(6,'(a, f10.3, a, f10.3)')                      &
      'Warning: zmax is less then zmin. Check input data zprofbgn',     &
      zlwr , 'and zprofend',zupr

  ALLOCATE(vartem(nz), STAT = ierr)
  ALLOCATE(zptem(nz),  STAT = ierr)
  ALLOCATE(varctem(2,2,nz), STAT = ierr)
  ALLOCATE(zpctem(2,2,nz),  STAT = ierr)
  CALL check_alloc_status(ierr, "vprofil:zpctem")

  varctem = 0.0
  zpctem  = 0.0

  DO ip = 1, npoints

    ix = INT ( (xpnt(ip) - xc(1,1,1))/dx ) + 1
    jy = INT ( (ypnt(ip) - yc(1,1,1))/dy ) + 1
    ix = MIN(MAX(1,ix), nxlg-2)
    jy = MIN(MAX(1,jy), nylg-2)

    IF(myproc == 0) WRITE(6,'(1x,2a,2(a,f10.3),2(a,i4) )')              &
        'Plotting ',xcaptn,' profile through (',                        &
        xpnt(ip),',',ypnt(ip),') km, at i=',ix,' j=',jy

    CALL mpupdatei(ix,1)
    CALL mpupdatei(jy,1)

    xp(1) = (ix-2)/(nx-3) + 1
    xp(2) = (ix-1)/(nx-3) + 1
    yp(1) = (jy-2)/(ny-3) + 1
    yp(2) = (jy-1)/(ny-3) + 1
    IF(xp(2) > nproc_x) xp(2) = nproc_x
    IF(yp(2) > nproc_y) yp(2) = nproc_y

    DO jj = 1,2
      DO ii = 1,2
        indxx = MOD((ix-2+ii-1),(nx-3)) + 2
        indxy = MOD((jy-2+jj-1),(ny-3)) + 2
        IF(ix+ii-1 > nxlg-2) indxx = nx-1
        IF(jy+jj-1 > nylg-2) indxy = ny-1

        source = xp(ii) + (yp(jj)-1)*nproc_x -1

        vartem = 0.0
        zptem  = 0.0

        CALL inctag
        IF (myproc == source) THEN
          xtem = xc(indxx,indxy,1)
          ytem = yc(indxx,indxy,1)
          vartem(:) = var(indxx,indxy,1:nz)
          zptem(:) = zpc(indxx,indxy,1:nz)

          itags = gentag
          CALL mpsendr(xtem,1,destination,itags,ierr)
          itags = gentag + 1
          CALL mpsendr(ytem,1,destination,itags,ierr)
          itags = gentag + 2
          CALL mpsendr(vartem,nz,destination,itags,ierr)
          itags = gentag + 3
          CALL mpsendr(zptem,nz,destination,itags,ierr)
        END IF

        IF(myproc == 0) THEN
          itagr = gentag
          CALL mprecvr(xtem,1,source,itagr,ierr)
          itagr = gentag + 1
          CALL mprecvr(ytem,1,source,itagr,ierr)
          itagr = gentag + 2
          CALL mprecvr(vartem,nz,source,itagr,ierr)
          itagr = gentag + 3
          CALL mprecvr(zptem,nz,source,itagr,ierr)

          a(ii,jj) = ABS( (xtem-xpnt(ip))*(ytem-ypnt(ip)) )
          varctem(ii,jj,:) = vartem(:)
          zpctem(ii,jj,:)  = zptem(:)
        END IF

      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Interpolate the data value and its height to the specified point.
!
!-----------------------------------------------------------------------
!
    IF( myproc == 0) THEN

      DO k = 1,nzprofc

        profil(k,ip)= (a(2,2)*varctem(1,1,k) + a(2,1)*varctem(1,2,k)+    &
                       a(1,2)*varctem(2,1,k) + a(1,1)*varctem(2,2,k))    &
                      /(a(1,1) + a(1,2) + a(2,1) + a(2,2))

        height(k,ip)= (a(2,2)*zpctem(1,1,k) + a(2,1)*zpctem(1,2,k)+      &
                       a(1,2)*zpctem(2,1,k) + a(1,1)*zpctem(2,2,k))      &
                      /(a(1,1) + a(1,2) + a(2,1) + a(2,2))

      END DO

    END IF  ! myproc == 0

  END DO

  DEALLOCATE(vartem,  zptem)
  DEALLOCATE(varctem, zpctem)

  IF(myproc == 0) THEN

    kbgn = nzprofc
    DO k=nzprofc,1,-1

      hmaxk = height(k,1)
      DO ip=1,npoints
        hmaxk = MAX(hmaxk,height(k,ip))
      END DO

      IF( hmaxk >= zmin) kbgn = k

    END DO

    kend = 1
    DO k=1,nzprofc

      hmink = height(k,1)
      DO ip=1,npoints
        hmink = MIN(hmink,height(k,ip))
      END DO

      IF( hmink <= zmax) kend=k

    END DO

!
!-----------------------------------------------------------------------
!
!  If input bounds for the profile are zero, use the min. and max.
!  in the profile as the lower and upper bounds for the horizontal
!  axis.
!
!-----------------------------------------------------------------------
!
    IF( plwr == 0.0 .AND. pupr == 0.0 ) THEN

      lower = profil(kbgn,1)
      upper = profil(kend,1)
      DO ip=1,npoints
        DO k = kbgn,kend
          lower = MIN(lower, profil(k,ip))
          upper = MAX(upper, profil(k,ip))
        END DO
      END DO

    ELSE

      lower = plwr
      upper = pupr

    END IF

!
!-----------------------------------------------------------------------
!
!    If the lower and upper bounds are equal, set the horizontal
!    axis scale to 1.0.
!
!-----------------------------------------------------------------------
!

   IF ((lower == 0.0 .AND. upper == 0.0).OR.upper == lower) upper = lower+1.0

!
!-----------------------------------------------------------------------
!
!  Start to plot the profile...
!
!-----------------------------------------------------------------------
!
    DO ip=1,npoints

      IF( (.NOT.multiprof) .OR. (multiprof.AND.ip == 1) ) THEN

        CALL xnwpic
        CALL xaxtik(1, 1)
        CALL xaxant(-1, -1)
        CALL xmap (lower, upper, zmin, zmax)

        CALL xaxnsz ( axlbsiz*(zmax-zmin)*lblmag )

        CALL xqmap(x1,x2,y1,y2)
        CALL xchsiz(0.03*(y2-y1)*lblmag)
        CALL xchori(0.0)

        IF( .NOT.multiprof ) THEN
          lchar = LEN( xcaptn)
          ch = xcaptn
          WRITE(ch(lchar+1:lchar+33), '(a,f13.3,a,f13.3,a)')              &
              ' at (',xpnt(ip),',',ypnt(ip),')'
          lchar = lchar+33
          CALL strmin(ch(1:lchar), lchar)
          CALL xcharc((x1+x2)*0.5, y1-(y2-y1)*0.10*px/py, ch(1:lchar))
        ELSE
          CALL xcharc((x1+x2)*0.5, y1-(y2-y1)*0.10*px/py, xcaptn )
        END IF
!
!-----------------------------------------------------------------------
!
!  Check if the points lie on one side of the axis. If the points are
!  all positive, draw the y-axis on the left border, if all points are
!  negative, draw the y-axis on the right border. If points lie on both
!  sides, draw the y-axis through x=0.0.
!
!-----------------------------------------------------------------------
!
        temp = lower * upper

        IF (temp > 0.0) THEN

          IF (lower > 0.0) THEN
            CALL xaxes(lower,0.0,zmin,0.0)
            CALL xchori(90.0)
            CALL xcharc(x1-0.12*(x2-x1), (y1+y2)*0.5, ycaptn)
          ELSE
            CALL xaxant(-1, 1)
            CALL xaxtik(1, -1)
            CALL xaxes(upper,0.0,zmin,0.0)
            CALL xchori(90.0)
            CALL xcharc(x1-0.05*(x2-x1), (y1+y2)*0.5, ycaptn)
          END IF

        ELSE
          CALL xaxes(0.0,0.0,zmin,0.0)
          CALL xchori(90.0)
          CALL xcharc(x1-0.10*(x2-x1), (y1+y2)*0.5, ycaptn)
        END IF

      END IF

      CALL xchori(0.0)
      CALL xbordr
      CALL xfull
!
!-----------------------------------------------------------------------
!
!  The first plot is labeled `A'. The subsequent plots will be `B'...
!
!-----------------------------------------------------------------------
!
      IF( multiprof ) THEN
        CALL xlbon
        CALL xlabel(CHAR(64+ip))

        ch(1:1) =  CHAR(64+ip)
        ch(2:2) =  ' '
        lchar = 2
        WRITE(ch(lchar+1:lchar+33), '(a,f13.2,a,f13.2,a)')                &
            ' at (',xpnt(ip),',',ypnt(ip),')'
        lchar = lchar+33
        CALL strmin(ch(1:lchar), lchar)

        CALL xqmap(x1,x2,y1,y2)
        CALL xchsiz(0.025*(y2-y1)*lblmag)
        CALL xchori(0.0)
        CALL xcharl(x1+(x2-x1)*0.03, y2-(y2-y1)*px/py*(0.03+0.035*ip),          &
            ch(1:lchar))

        CALL xlbsiz( ctrlbsiz*(y2-y1)*lblmag )
      ELSE
        CALL xlboff
      END IF

      CALL xwindw(lower, upper, zmin, zmax)

      CALL xqmap(x1,x2,y1,y2)

      CALL xcurve(profil(kbgn,ip),height(kbgn,ip),kend-kbgn+1,0)
      CALL xwdwof

    END DO  ! ip

  END IF   ! myproc == 0

  RETURN
END SUBROUTINE vprofil


!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SPLTPARA                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE spltpara(inc,MIN,MAX,ovr,hlf,zro,col1,col2,pltvar)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Set some parameters for one plot.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!    Min Zou (3/2/98)
!
!  Modification history:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!    inc      interval of the contour
!    min      the minimum value for the contour
!    max      the maximum valur foj the contour
!    ovr      overlay option
!    hlf      the contour highlight frequency
!    zro      define the attributes of zero contours
!    col1     the start color index for contour
!    col2     the end color index for contour
!    pltvar   the plot name
!    len      the length of pltvar
!
!-----------------------------------------------------------------------
!
  INTEGER :: ovr,hlf,zro,col1,col2
  REAL :: inc, MIN, MAX
  CHARACTER (LEN=12) :: pltvar

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  CALL ctrinc( inc, MIN, MAX )
  CALL overlay(ovr)
  CALL xhlfrq (hlf)
  CALL xczero (zro)
  CALL ctrcol( col1,col2)
  CALL varplt( pltvar )

  RETURN
END SUBROUTINE spltpara

SUBROUTINE fillmissval (m, n, xl, xr, yb,yt )

  REAL :: x1,x2,y1,y2
  REAL :: xra(4), yra(4)
  INTEGER :: missval_colind, missfill_opt    ! miss value color index
  COMMON /multi_value/ missfill_opt,missval_colind

  x1 = xl + (xr-xl)/REAL(m)*0.5
  x2 = xr - (xr-xl)/REAL(m)*0.5
  y1 = yb + (yt-yb)/REAL(n)*0.5
  y2 = yt - (yt-yb)/REAL(n)*0.5

  xra(1) = x1
  xra(2) = x2
  xra(3) = x2
  xra(4) = x1
  yra(1) = y1
  yra(2) = y1
  yra(3) = y2
  yra(4) = y2

  CALL xcolor(missval_colind)
  CALL xfilarea(xra, yra, 4)

  RETURN
END SUBROUTINE fillmissval
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE HINTRP                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE hintrp(nx,ny,nz,a3din,z3d,zlevel, a2dout)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Interpolate a 3-D array to horizontal level z=zlevel.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  Based on original SECTHRZ.
!  12/10/98.
!
!  MODIFICATION HISTORY:
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
!    a3din    3-d input array
!    z3d      z-coordinate of data in a3din
!    zlevel   Level to which data is interpolated.
!
!  OUTPUT:
!    a2dout   2-d output array interpolated to zlevel
!
!-----------------------------------------------------------------------
!
!  Parameters of output
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nx,ny,nz

  REAL :: a3din(nx,ny,nz) ! 3-d input array
  REAL :: z3d  (nx,ny,nz) ! z-coordinate of data in a3din
  REAL :: zlevel          ! Level to which data is interpolated.

  REAL :: a2dout(nx,ny)   ! 2-d output array interpolated to zlevel

  INTEGER :: i,j,k
!
!-----------------------------------------------------------------------
!
!  Find index for interpolation
!
!-----------------------------------------------------------------------
!
  DO i=1,nx-1
    DO j=1,ny-1
      IF(zlevel <= z3d(i,j,1)) GO TO 11
      IF(zlevel >= z3d(i,j,nz-1)) GO TO 12
      DO k=2,nz-2
        IF(zlevel >= z3d(i,j,k).AND.zlevel < z3d(i,j,k+1)) GO TO 15
      END DO

      11    k=2
      GO TO 15
      12    k=nz-1
      GO TO 15

      15    a2dout(i,j)=a3din(i,j,k)+(a3din(i,j,k+1)-a3din(i,j,k))*     &
                        (zlevel-z3d(i,j,k))/(z3d(i,j,k+1)-z3d(i,j,k))

!-----------------------------------------------------------------------
!
!  If the data point is below the ground level, set the
!  data value to the missing value.
!
!-----------------------------------------------------------------------

      IF( zlevel < z3d(i,j,2)   ) a2dout(i,j) = -9999.0
      IF( zlevel > z3d(i,j,nz-1)) a2dout(i,j) = -9999.0

    END DO
  END DO

  RETURN
END SUBROUTINE hintrp
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE HINTRP1                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE hintrp1(nx,ny,nz, kbgn,kend,a3din,z3d,zlevel, a2dout)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Interpolate a 3-D array to horizontal level z=zlevel.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  Based on original SECTHRZ.
!  12/10/98.
!
!  MODIFICATION HISTORY:
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
!    kbgn
!    kend
!
!    a3din    3-d input array
!    z3d      z-coordinate of data in a3din
!    zlevel   Level to which data is interpolated.
!
!  OUTPUT:
!    a2dout   2-d output array interpolated to zlevel
!
!-----------------------------------------------------------------------
!
!  Parameters of output
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  INTEGER :: kbgn, kend

  REAL :: a3din(nx,ny,nz) ! 3-d input array
  REAL :: z3d  (nx,ny,nz) ! z-coordinate of data in a3din
  REAL :: zlevel          ! Level to which data is interpolated.

  REAL :: a2dout(nx,ny)   ! 2-d output array interpolated to zlevel

  INTEGER :: i,j,k
!
!-----------------------------------------------------------------------
!
!  Find index for interpolation
!
!-----------------------------------------------------------------------
!
  DO i=1,nx-1
    DO j=1,ny-1
      IF(zlevel <= z3d(i,j,kbgn)) GO TO 11
      IF(zlevel >= z3d(i,j,kend)) GO TO 12
      DO k=kbgn,kend-1
        IF(zlevel >= z3d(i,j,k) .AND. zlevel < z3d(i,j,k+1)) GO TO 15
      END DO

      11    k=kbgn
      GO TO 15
      12    k=kend-1
      GO TO 15

      15    a2dout(i,j)=a3din(i,j,k)+(a3din(i,j,k+1)-a3din(i,j,k))*     &
                        (zlevel-z3d(i,j,k))/(z3d(i,j,k+1)-z3d(i,j,k))

!-----------------------------------------------------------------------
!
!  If the data point is below the ground level, set the
!  data value to the missing value.
!
!-----------------------------------------------------------------------

      IF( zlevel < z3d(i,j,kbgn) ) a2dout(i,j) = -9999.0
      IF( zlevel > z3d(i,j,kend) ) a2dout(i,j) = -9999.0

    END DO
  END DO

  RETURN
END SUBROUTINE hintrp1
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE HINTRP2                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE hintrp2(nx,ny,nz, kbgn,kend,a3din,z3d,zlevel,hterain, a2dout)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Interpolate a 3-D array to horizontal level z=zlevel.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  Based on original SECTHRZ.
!  12/10/98.
!
!  MODIFICATION HISTORY:
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
!    kbgn
!    kend
!
!    a3din    3-d input array
!    z3d      z-coordinate of data in a3din
!    zlevel   Level to which data is interpolated.
!
!  OUTPUT:
!    a2dout   2-d output array interpolated to zlevel
!
!-----------------------------------------------------------------------
!
!  Parameters of output
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  INTEGER :: kbgn, kend

  REAL :: a3din(nx,ny,nz) ! 3-d input array
  REAL :: z3d  (nx,ny,nz) ! z-coordinate of data in a3din
  REAL :: zlevel          ! Level to which data is interpolated.
  REAL :: hterain(nx,ny)   ! terrain height in meters

  REAL :: a2dout(nx,ny)   ! 2-d output array interpolated to zlevel

  INTEGER :: i,j,k
  REAL :: zlevel1
!
!-----------------------------------------------------------------------
!
!  Find index for interpolation
!
!-----------------------------------------------------------------------
!
  DO i=1,nx-1
    DO j=1,ny-1
      zlevel1 = zlevel + hterain(i,j)*0.001

!     if( i ==2 .and. j==ny-2 ) print*,'zlevel1 at 2, ny-2 =', zlevel1, hterain(i,j)

      IF(zlevel1 <= z3d(i,j,kbgn)) GO TO 11
      IF(zlevel1 >= z3d(i,j,kend)) GO TO 12


      DO k=kbgn,kend-1
        IF(zlevel1 >= z3d(i,j,k).AND.zlevel1 < z3d(i,j,k+1)) GO TO 15
      END DO


      11    k=kbgn
      GO TO 15
      12    k=kend-1
      GO TO 15

      15    a2dout(i,j)=a3din(i,j,k)+(a3din(i,j,k+1)-a3din(i,j,k))*     &
                        (zlevel1-z3d(i,j,k))/(z3d(i,j,k+1)-z3d(i,j,k))

!-----------------------------------------------------------------------
!
!  If the data point is below the ground level, set the
!  data value to the missing value.
!
!-----------------------------------------------------------------------

      IF( zlevel1 < z3d(i,j,kbgn) ) a2dout(i,j) = -9999.0
      IF( zlevel1 > z3d(i,j,kend) ) a2dout(i,j) = -9999.0

    END DO
  END DO

  RETURN
END SUBROUTINE hintrp2

SUBROUTINE indxbnds(xc,yc,zpc,zpsoilc,nx,ny,nz,nzsoil,                  &
           xbgn,xend,ybgn,yend,zbgn,zend,zsoilbgn,zsoilend,             &
           ibgn,iend,jbgn,jend,kbgn,kend,ksoilbgn,ksoilend)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Return index bounds of the domain to be plotted
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  INTEGER :: nzsoil
  REAL :: xc   (nx,ny,nz)   ! x-coor of sacalar point (km)
  REAL :: yc   (nx,ny,nz)   ! y-coor of sacalar point (km)
  REAL :: zpc  (nx,ny,nz)   ! z-coor of sacalar point in physical
                            ! space (km)
  REAL :: zpsoilc(nx,ny,nzsoil)   ! z-coor of sacalar point in physical
                            ! space (m) for soil model
  REAL :: xbgn,xend,ybgn,yend,zbgn,zend,zsoilbgn,zsoilend
  INTEGER :: ibgn,iend,jbgn,jend,kbgn,kend,ksoilbgn,ksoilend
  INTEGER :: i,j,k

!----------------------------------------------------------------------
!
! Include files
!
!----------------------------------------------------------------------
  INCLUDE 'mp.inc'
!
!----------------------------------------------------------------------
  INTEGER :: nxlg,nylg
  INTEGER :: istatus


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begining of executable code ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  nxlg = (nx-3)*nproc_x+3
  nylg = (ny-3)*nproc_y+3

  IF(xbgn /= xend) THEN

    IF (xbgn >= xc(1,2,2) .AND. xbgn < xc(nx-2,2,2)) THEN
      DO i = 1,nx-2
        IF (xc(i,2,2) >= xbgn) EXIT  ! find local ibgn
      END DO
      ibgn = (nx-3)*(loc_x-1) + i    ! find global ibgn
    ELSE
      ibgn = -1
    END IF

    IF (xend > xc(2,2,2) .AND. xend <= xc(nx-1,2,2)) THEN
      DO i = nx-1,2,-1
        IF (xc(i,2,2) <= xend) EXIT
      END DO
      iend = (nx-3)*(loc_x-1) + i
    ELSE
      iend = -1
    END IF

    CALL mpmaxi(ibgn)
    CALL mpmaxi(iend)

  ELSE

    ibgn = 2
    iend = nxlg- 2
    IF (myproc == 0)          xbgn = xc(   2,2,2)
    IF (myproc == nproc_x-1)  xend = xc(nx-2,2,2)
    CALL mpbcastr(xbgn,0)
    CALL mpbcastr(xend,nproc_x-1)

  END IF

  IF(ybgn /= yend) THEN

    IF (ybgn >= yc(2,1,2) .AND. ybgn < yc(2,ny-2,2) ) THEN
      DO j = 1,ny-2
        IF (yc(2,j,2) >= ybgn) EXIT
      END DO
      jbgn = (ny-3)*(loc_y-1) + j
    ELSE
      jbgn = -1
    END IF

    IF (yend > yc(2,2,2) .AND. yend <= yc(2,ny-1,2) ) THEN
      DO j = ny-1,2,-1
        IF (yc(2,j,2) <= yend) EXIT
      END DO
      jend = (ny-3)*(loc_y-1) + j
    ELSE
      jend = -1
    END IF

    CALL mpmaxi(jbgn)
    CALL mpmaxi(jend)

  ELSE

    jbgn = 2
    jend = nylg-2
    IF (loc_x == 1 .AND. loc_y == 1) THEN  ! processor 0
      ybgn = yc(2,   2,2)
    END IF
    IF (loc_x == 1 .AND. loc_y == nproc_y) THEN  ! processor (nporc_y-1)*nproc_x
      yend = yc(2,ny-2,2)
    END IF
    CALL mpbcastr(ybgn,0)
    CALL mpbcastr(yend,(nproc_y-1)*nproc_x)
  END IF

  IF(zbgn /= zend) THEN
    kend = 2
    loop_kend: DO k = 2,nz-1
      loop_jend: DO j = 2,ny-2
        loop_iend: DO i = 2,nx-2
          IF(zpc(i,j,k) < zend) THEN
            kend=k
            CYCLE loop_kend
          END IF
        END DO loop_iend
      END DO loop_jend
      EXIT loop_kend
    END DO loop_kend
    kend = MIN(kend+1, nz-1)

    kbgn= nz-1
    loop_kbgn: DO k = nz-1,2,-1
      loop_jbgn: DO j = 2,ny-2
        loop_ibgn: DO i = 2,nx-2
          IF(zpc(i,j,k) > zbgn) THEN
            kbgn=k
            CYCLE loop_kbgn
          END IF
        END DO loop_ibgn
      END DO loop_jbgn
      EXIT loop_kbgn
    END DO loop_kbgn
    kbgn = MAX(kbgn,2)

    CALL mpmax0i(kend,kbgn)

  ELSE

    kbgn = 2
    kend = nz-2

  END IF

  IF(zsoilbgn /= zsoilend) THEN
!
! 05/31/2002 Zuwen He
!
! Note: k is 1 at the surface in the soil model,
!       and k increase when zpsoilc decrease.
!       zpsoilc=zpsoil(k)-zpsoil(1) < 0.
!
      ksoilend = 1
      DO k = 1,nzsoil
        DO j = 2,ny-2
          DO i = 2,nx-2
            IF(zpsoilc(i,j,k) > zsoilend) THEN
              ksoilend=k
              GO TO 325
            END IF
          END DO
        END DO
        GO TO 335
  325   CONTINUE
      END DO
  335 ksoilend = MIN(ksoilend+1, nzsoil)

      ksoilbgn= nzsoil
      DO k = nzsoil,1,-1
        DO j = 2,ny-2
          DO i = 2,nx-2
            IF(zpsoilc(i,j,k) < zsoilbgn) THEN
              ksoilbgn=k
              GO TO 350
            END IF
          END DO
        END DO
        GO TO 345
  350   CONTINUE
      END DO
  345 ksoilbgn = MAX(ksoilbgn-1,1)


    CALL mpmax0i(ksoilbgn,ksoilend)

  ELSE

    ksoilbgn = 1
    ksoilend = nzsoil

  END IF

  IF(myproc == 0) WRITE(6,'(/1x,a,i3,a,i5)') 'ibgn =',ibgn,', iend =',iend
  IF(iend < ibgn) THEN
    IF(myproc == 0) WRITE(6,'(1x,a,/1x,a)')                             &
        'iend was found smaller than ibgn. Check the input',            &
        'domain bounds in x direction. Program stopped.'
    CALL arpsstop('ibgn & iend error inside indxbnds.',1)
  END IF

  IF(myproc == 0) WRITE(6,'(1x,a,i3,a,i5)') 'jbgn =',jbgn,', jend =',jend
  IF(jend < jbgn) THEN
    IF(myproc == 0) WRITE(6,'(1x,a,/1x,a)')                             &
        'jend was found smaller than jbgn. Check the input',            &
        'domain bounds in y direction. Program stopped.'
    CALL arpsstop('jbgn & jend error inside indxbnds.',1)
  END IF

  IF(myproc == 0)  WRITE(6,'(1x,a,i3,a,i5)') 'kbgn =',kbgn,', kend =',kend
  IF(kend < kbgn) THEN
    IF(myproc == 0) WRITE(6,'(1x,a,/1x,a)')                             &
        'kend was found smaller than kbgn. Check the input',            &
        'domain bounds in z direction. Program stopped.'
    CALL arpsstop('kbgn & kend error inside indxbnds.',1)
  END IF

  IF(myproc == 0) WRITE(6,'(1x,a,i2,a,i2)') 'ksoilbgn =',  ksoilbgn,    &
                                            ', ksoilend =',ksoilend
  IF(ksoilend < ksoilbgn) THEN
    IF(myproc == 0)  WRITE(6,'(1x,a,/1x,a)')                            &
        'ksoilend was found smaller than ksoilbgn. Check the input',    &
        'domain bounds in zpsoil direction. Program stopped.'
    CALL arpsstop('ksoilbgn & ksoilend error inside indxbnds.',1)
  END IF

  RETURN
END SUBROUTINE indxbnds

SUBROUTINE ctrsetup(zinc,zminc,zmaxc,                                   &
                    zovr,zhlf,zzro,zcol1,zcol2,zlabel)

  IMPLICIT NONE
  REAL    :: zinc,zminc,zmaxc
  INTEGER :: zovr,zhlf,zzro,zcol1,zcol2
  CHARACTER (LEN=*) :: zlabel

  IF(zhlf <= 0.0) THEN
    WRITE(6,'(/4a/a)') 'ERROR: ZHLF must be a positive value for "',    &
               TRIM(zlabel),'".',                                       &
               'Please check your input file. Program Stopping...'
    STOP
  END IF

  CALL ctrinc (zinc,zminc,zmaxc )
  CALL overlay(zovr)
  CALL xhlfrq (zhlf)
  CALL xczero (zzro)
  CALL ctrcol (zcol1,zcol2 )
  CALL varplt (zlabel)

  RETURN
END SUBROUTINE ctrsetup
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE PLTTRN                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE plttrn(hterain,x,y,m,n,slicopt,iwrk,xwk,ywk)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Generate terrain contours
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    hterain  2-D terrain data to contour
!    x        x coordinate of grid points in plot space (over on page)
!    y        y coordinate of grid points in plot space (up on page)
!    m        first dimension
!    n        second dimension
!
!    slicopt  slice orientation indicator
!       slicopt = 1, x-y slice of u,v at z index kslice is plotted.
!       slicopt = 2, x-z slice of u,w at y index jslice is plotted.
!       slicopt = 3, y-z slice of v,w at x index islice is plotted.
!       slicopt = 4, x-y slice of u,v at z index islice is plotted.
!       slicopt = 5, xy-z cross section of wind islice is plotted.
!       slicopt = 6, data field on constant p-level is plotted.
!       slicopt = 0, all of the three slices above are plotted.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: m,n
  REAL    :: hterain(m,n)
  REAL    :: x(m,n)
  REAL    :: y(m,n)
  INTEGER :: slicopt
  INTEGER, INTENT(INOUT) :: iwrk(m,n)
  REAL   , INTENT(INOUT) :: xwk(m,n),ywk(m,n)
!
!-----------------------------------------------------------------------
!
!  Plotting control common blocks
!
!-----------------------------------------------------------------------
!
  REAL :: ctinc,ctmin,ctmax,vtunt  ! contour interval and vector unit
  COMMON /incunt/ ctinc,ctmin,ctmax,vtunt

  INTEGER :: icolor,icolor1,lbcolor,trcolor                ! required color
  COMMON /recolor/icolor,icolor1,lbcolor,trcolor

  REAL :: ztmin,ztmax
  INTEGER :: ovrtrn ,trnplt       ! overlay terrain option (0/1)
  REAL :: trninc,trnmin, trnmax   ! terrain interval minimum, maximum
  COMMON /trnpar/ trnplt,ovrtrn,trninc,trnmin, trnmax,ztmin,ztmax
  REAL :: zlevel
  COMMON /sliceh/zlevel
  INTEGER :: col_table,pcolbar
  COMMON /coltable/col_table,pcolbar

  REAL :: pl, pr, pb, pt
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------

  REAL :: cl(900)
  INTEGER :: ncl

  REAL :: z02,xl,xr,yt,yb,xfinc
  INTEGER :: mode1

  INTEGER :: istatus
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF( ovrtrn == 0 .OR. ztmax-ztmin < 1.0E-20  ) RETURN
!
!-----------------------------------------------------------------------
!
!  Overlay terrain contour if required in x-y level
!  or Plot terrain outline in slice zlevel
!
!-----------------------------------------------------------------------
!
  CALL xqmap(xl,xr,yb,yt)
  CALL xqpspc( pl, pr, pb, pt)

  cl(1)=0.0
  IF(slicopt == 1 .OR. slicopt == 8 .OR. slicopt == 9)  THEN

    CALL ctrinc( trninc,trnmin, trnmax )
    IF( trninc == 0.0) THEN
      cl(2)=cl(1)+ xfinc(ztmax-ztmin)/2
      IF(cl(2)-cl(1) == 0.0) cl(2)=cl(1)+1.0
      mode1=1
      CALL xnctrs(6,18)
    ELSE
      cl(2)=cl(1)+trninc
      CALL xnctrs(1,900)
      IF(ztmin == 0.0 .AND. ztmax == 0.0) THEN
        mode1=1
      ELSE
        mode1=3
      END IF
    END IF

    CALL xctrlim(ctmin,ctmax)
    IF (trnplt == 1) THEN
      CALL xthick(2)
      CALL xctrclr(trcolor, trcolor)
      IF(mode1 == 3) THEN
        ncl=FLOOR( (ztmax-ztmin)/trninc ) + 1
        cl(1)=ztmin
        cl(2)=cl(1)+trninc
      END IF
      CALL xconta(hterain,x,y,iwrk,m,m,n,cl,ncl,mode1)
    ELSE IF (trnplt == 2) THEN
      CALL xctrclr(icolor, icolor1)
      IF(mode1 == 3) THEN
        ncl=FLOOR( (ztmax-ztmin)/trninc ) + 1
        cl(1)=ztmin
        cl(2)=cl(1)+trninc
      END IF
      CALL xcolfil(hterain,x,y,iwrk,xwk,ywk,m,m,n,cl,ncl,mode1)
      CALL xchmag(0.025*sqrt((pr-pl)*(pt-pb)))
      CALL xcpalet(pcolbar)
    ELSE IF (trnplt == 4) THEN
      CALL xctrclr(icolor, icolor1)
      CALL xconta(hterain,x,y,iwrk,m,m,n,cl,ncl,mode1)
    END IF
  ELSE IF(slicopt == 4.OR.slicopt == 6.OR.slicopt == 7) THEN
    CALL xcolor(trcolor)
    z02=zlevel*1000.
    CALL xthick(2)
    CALL xcontr(hterain,x,y,iwrk,m,m,n,z02)
    CALL xthick(1)
  END IF

  RETURN
END SUBROUTINE plttrn

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE PLTAXES                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE pltaxes(slicopt,dx,dy)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!-----------------------------------------------------------------------
!
!  AUTHOR:  M. Xue
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!
!  INPUT:
!    slicopt slice orientation indicator
!         = 1, x-y slice of at k=kslice is plotted.
!         = 2, x-z slice of at j=jslice is plotted.
!         = 3, y-z slice of at i=islice is plotted.
!         = 4, horizontal slice at z index islice is plotted.
!         = 5, xy-z cross section of wind islice is plotted.
!         = 6, data field on constant p-level is plotted.
!         = 0, all of the three slices above are plotted.
!    dx   Spacing between the x-axis tick marks
!    dy   Spacing between the y-axis tick marks
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  REAL,    INTENT(IN) :: dx,dy
  INTEGER, INTENT(IN) :: slicopt
!
!-----------------------------------------------------------------------
!
!  Plotting control common blocks
!
!-----------------------------------------------------------------------
!
  INTEGER :: layover
  COMMON /laypar/ layover

  INTEGER :: presaxis_no
  REAL :: pres_val(20), pres_z(20)
  COMMON /pressbar_par/presaxis_no,pres_val,pres_z

  REAL :: lblmag, ctrlbsiz, axlbsiz
  COMMON /labmag/ lblmag, ctrlbsiz, axlbsiz

  INTEGER :: timeovr
  COMMON /timover/ timeovr

  REAL :: x101, y101, x102,y102
  COMMON /slicev1/x101, y101, x102,y102

  INTEGER :: xfont   ! the font of character
  INTEGER :: haxisu, vaxisu
  INTEGER :: lbaxis
  INTEGER :: tickopt
  INTEGER :: axlbfmt
  REAL :: hmintick,vmajtick,vmintick,hmajtick
  COMMON /var_par/ xfont,haxisu,vaxisu,lbaxis,tickopt,hmintick,         &
         vmajtick, vmintick,hmajtick,axlbfmt

  INTEGER :: icolor,icolor1,lbcolor,trcolor       ! required color
  COMMON /recolor/icolor,icolor1,lbcolor,trcolor

  INTEGER :: mapgrid,mapgridcol
  REAL    :: latgrid,longrid
  COMMON /mappar2/mapgrid,mapgridcol,latgrid,longrid
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  REAL    :: axnmag
  REAL    :: xl,xr,yb,yt,pl,pr,pb,pt

  REAL    :: xtem1, xtem2               !local temporary variable
  REAL    :: x1,x2, y1,y2, xstep, ystep, xmstep, ymstep
  INTEGER :: LEN
  CHARACTER (LEN=16) :: ylabel
  CHARACTER (LEN=16) :: xlabel

  REAL :: px, py

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!
!-----------------------------------------------------------------------
!
!  Set-up variables for tick marks and draw axes
!
!-----------------------------------------------------------------------
!
  CALL xqmap(xl,xr, yb,yt)

  CALL xqpspc( pl, pr, pb, pt)
  px = pr - pl
  py = pt - pb

  CALL setcords(xl,xr,yb,yt,dx,dy, slicopt,                             &
       x1,x2,y1,y2,xlabel,ylabel,xstep,ystep,xmstep,ymstep)

  CALL xqpspc( pl, pr, pb, pt)
  axnmag = axlbsiz*MIN(pt-pb, pr-pl)*lblmag

  CALL xaxnmg( axnmag )

  IF(slicopt == 5) THEN
    IF( ABS(y101-y102) <= 1.0E-3 ) THEN
      xtem1 = x101
      xtem2 = x102
    ELSE IF(ABS(x101-x102) <= 1.0E-3 ) THEN
      xtem1 = y101
      xtem2 = y102
    ELSE
!
!     xtem1 = SQRT(x101*x101 + y101*y101)

! Recommend changing the origin of the horizontal axis to zero in future
! version as follows.
       xtem1 = 0.0

      xtem2 = xtem1 + SQRT( (y102-y101)*(y102-y101) +                   &
                            (x102-x101)*(x102-x101) )
    END IF
  ELSE
    xtem1 = x1
    xtem2 = x2
  END IF

  CALL xaxant(-1,-1)    ! explicitly set default axis annotation
  CALL xaxtik( 1, 1)    ! set default axis ticks
  CALL xmap(xtem1,xtem2,y1,y2)
  IF( layover == 0) THEN

    IF ( (latgrid <   0.0  .OR.  longrid <    0.0) .AND.                &
         (latgrid > -900.0 .AND. longrid > -900.0) .AND.                &
         (slicopt == 1 .OR. slicopt == 4 .OR. slicopt == 6 .OR.         &
          slicopt == 7 .OR. slicopt == 8 .OR. slicopt == 9 ) ) THEN
      CALL xaxant(0,0)  ! suppress axis annotation if lat/lon grid is desired
      CALL xaxtik(0,0)  ! suppress axis ticks
    END IF

    CALL xaxsor(0.0, 0.0)
!    call xthick(2)
    CALL xaxsca1( xtem1,xtem2,xstep,xmstep, y1,y2,ystep,ymstep )
    CALL xthick(1)
  END IF

!
!  Plot pressure axis
!
  CALL xqpspc(pl, pr, pb, pt)
  IF(presaxis_no > 0 .AND. timeovr == 0 .AND.                           &
     (slicopt == 2 .OR. slicopt == 3 .OR. slicopt == 5 .OR.             &
      slicopt == 10 .OR. slicopt ==11) ) THEN
    x1 = pl - (pr-pl)*0.25
    x2 = pl
    y1 = pb
    y2 = pt
    CALL xpspac(x1,x2,y1,y2)
    y1 = yb
    y2 = yt
    CALL xmap(x1,x2,y1,y2)
    CALL xaxfmt('(I4)')
    CALL xyaxis(x1+0.40*(x2-x1),pres_z,pres_val,presaxis_no)
    CALL xchori(90.)
    CALL xcharc(x1-0.10*(xr-xl),(yt+yb)*0.5 ,'Pressure(mb)')
    CALL xchori(0.)
  END IF
!
!  Restore the original plotting scape
!
  CALL xpspac( pl, pr, pb, pt)
  CALL xmap(xl,xr, yb,yt)

  IF(layover > 1) THEN
    CALL xchmag( 0.018*sqrt((pr-pl)*(pt-pb)) * lblmag )
  ELSE
    CALL xchmag( 0.020*sqrt((pr-pl)*(pt-pb)) * lblmag )
  END IF

  IF(lbaxis == 1 .AND. timeovr == 0) THEN
    CALL xcolor(lbcolor)
    LEN=LEN_TRIM(xlabel)
    CALL strmin(xlabel,LEN)
    CALL xcharc( xl+(xr-xl)*0.5,yb-0.08*(yt-yb)*px/py,xlabel(1:LEN))
    LEN=LEN_TRIM(ylabel)
    CALL strmin(ylabel,LEN)
    CALL xchori(90.)
    CALL xcharc(xl-0.10*(xr-xl),(yt+yb)*0.5,ylabel(1:LEN))
    CALL xchori(0.)
  END IF

  RETURN
END SUBROUTINE pltaxes
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE PLTEXTRA                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE pltextra(slicopt, pltopt)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Plot extra things such as map, boxes, polygons and stations
!  in a 2D plot
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!
!  MODIFICATION HISTORY:
!
!  6/08/92 (K. Brewster)
!  Added full documentation.
!
!  8/28/94 (M. Zou)
!  Added color routing , overlay terrain.
!
!  1/24/96 (J. Zong and M. Xue)
!  Fixed a problem related to finding the minimum and maximum of the
!  2D array, a, when there exist missing data. Initial min. and max.
!  should be set to values other than the missing value, -9999.0.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!    slicopt     slice orientation indicator
!             = 1, x-y slice of at k=kslice is plotted.
!             = 2, x-z slice of at j=jslice is plotted.
!             = 3, y-z slice of at i=islice is plotted.
!             = 4, horizontal slice at z index islice is plotted.
!             = 5, xy-z cross section of wind islice is plotted.
!             = 6, data field on constant p-level is plotted.
!             = 0, all of the three slices above are plotted.
!    plot
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: slicopt, pltopt

  INCLUDE 'arpsplt.inc'
!
!-----------------------------------------------------------------------
!
!  Plotting control common blocks
!
!-----------------------------------------------------------------------
!
  INTEGER :: ovrstaopt
  INTEGER :: ovrstam,staset,ovrstan,ovrstav,stacol,markprio,wrtstax
  INTEGER :: nsta_typ,sta_typ(30),sta_marktyp(30),sta_markcol(30)
  REAL :: sta_marksz(30),wrtstad
  CHARACTER (LEN=256) :: stalofl
  COMMON /sta_par/ ovrstaopt,ovrstam,staset,ovrstan,ovrstav,stacol,     &
         markprio,nsta_typ,sta_typ,sta_marktyp,                         &
         sta_markcol,sta_marksz,stalofl,wrtstax,wrtstad

  INTEGER :: icolor,icolor1,lbcolor,trcolor                ! required color
  COMMON /recolor/icolor,icolor1,lbcolor,trcolor
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: timeovr
  COMMON /timover/ timeovr

  INTEGER :: ovrmap
  COMMON /mappar / ovrmap

  INTEGER :: number_of_boxes,boxcol
  REAL :: bx1(10),bx2(10),by1(10),by2(10)
  COMMON /boxesopt/number_of_boxes,boxcol,bx1,bx2,by1,by2

  INTEGER :: num_of_verts
  INTEGER :: number_of_polys,polycol
  REAL :: vertx(max_verts,max_polys),verty(max_verts,max_polys)
  COMMON /polysopt/number_of_polys,polycol,vertx,verty

  REAL :: lblmag, ctrlbsiz, axlbsiz
  COMMON /labmag/ lblmag, ctrlbsiz, axlbsiz

  INTEGER :: nunit
  INTEGER :: i,j
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
!  Plot map boundaries.
!
!-----------------------------------------------------------------------
!
  CALL xcolor(lbcolor)

  IF((slicopt == 1 .OR. slicopt == 4 .OR. slicopt == 6 .OR.             &
      slicopt == 7 .OR. slicopt == 8 .OR. slicopt == 9)                 &
        .AND.ovrmap == 1                                                &
        .AND.(timeovr == 0 .OR. (timeovr == 1 .AND. pltopt == 2) ))THEN

    CALL getunit(nunit)

    CALL drawmap(nunit)

    CLOSE (UNIT=nunit)
    CALL retunit(nunit)

    CALL xthick(1)
  END IF
!
!-----------------------------------------------------------------------
!
!  Draw boxes
!
!-----------------------------------------------------------------------
!
  IF(number_of_boxes /= 0 .AND.  &
     (slicopt == 1 .OR. slicopt == 4 .OR. slicopt == 6 .OR.  &
      slicopt == 7 .OR. slicopt == 8 .OR. slicopt == 9)                &
        .AND. timeovr == 0 ) THEN
    CALL xthick(1)
    CALL xcolor(boxcol)
    CALL xbrokn(6,3,6,3)
    DO i=1,number_of_boxes
      CALL xbox(bx1(i),bx2(i),by1(i),by2(i))
    END DO
    CALL xthick(1)
    CALL xfull
  END IF
!
!-----------------------------------------------------------------------
!
!  Draw polylines
!
!-----------------------------------------------------------------------
!
  IF(number_of_polys /= 0 .AND.  &
     (slicopt == 1 .OR. slicopt == 4 .OR. slicopt == 6 .OR.  &
      slicopt == 7 .OR. slicopt == 8 .OR. slicopt == 9)                &
        .AND. timeovr == 0 ) THEN
    CALL xthick(2)
    CALL xcolor(polycol)
!    CALL xbrokn(6,3,6,3)
    DO j=1,number_of_polys
      num_of_verts=0
      DO i=1,max_verts
        IF(vertx(i,j) /= -9999. .AND. verty(i,j) /= -9999.)             &
                  num_of_verts = num_of_verts +1
      END DO
      IF(num_of_verts /= 0 ) CALL xcurve(vertx(1,j),verty(1,j),num_of_verts, 0)
    END DO
    CALL xthick(1)
    CALL xfull
  END IF

  RETURN
END SUBROUTINE pltextra
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SMOOTH9PMV                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!

SUBROUTINE smooth9pmv( arr, nx,ny,ibgn,iend,jbgn,jend, tem1 )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!                                        1 2 1
!  Smooth a 2-D array by the filter of { 2 4 2 }
!                                        1 2 1
!
!-----------------------------------------------------------------------
!
!  AUTHOR:       Yuhe Liu
!
!  5/3/94
!
!  Modification History
!  8/20/1995 (M. Xue)
!  Fixed errors in the index bound of loops 100 and 200.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  nx       Number of grid points in the x-direction
!  ny       Number of grid points in the y-direction
!  ibgn     First index in x-direction in the soomthing region.
!  iend     Last  index in x-direction in the soomthing region.
!  jbgn     First index in j-direction in the soomthing region.
!  jend     Last  index in j-direction in the soomthing region.
!
!  arr    2-D array
!
!  OUTPUT:
!
!  arr    2-D array
!
!  TEMPORARY:
!
!  tem1     Temporary 2-D array
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: nx         ! Number of grid points in the x-direction
  INTEGER :: ny         ! Number of grid points in the y-direction
  INTEGER :: ibgn
  INTEGER :: iend
  INTEGER :: jbgn
  INTEGER :: jend
!
  REAL :: arr (nx,ny)   ! 2-D array
!
  REAL :: tem1(nx,ny)   ! Temporary array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,ip,im,jp,jm
  REAL :: wtf,mv
!
!-----------------------------------------------------------------------
!
!  Include files:
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
  wtf    = 1.0/16.0

  mv = -9999.0 ! missing value flag

  DO i=1,nx
    DO j=1,ny
      IF( ABS(arr(i,j)-mv) <= 0.1) arr(i,j)=mv
    END DO
  END DO

  DO j = jbgn,jend
    DO i = ibgn,iend
      ip=MIN(nx,i+1)
      im=MAX( 1,i-1)
      jp=MIN(ny,j+1)
      jm=MAX( 1,j-1)

      tem1(i,j) = wtf                                                   &
          * (    arr(im,jm) + 2.*arr(i,jm) +    arr(ip,jm)              &
          + 2.*arr(im,j ) + 4.*arr(i,j ) + 2.*arr(ip,j )                &
          +    arr(im,jp) + 2.*arr(i,jp) +    arr(ip,jp))

      IF(arr(im,jm) == mv.OR.arr(i,jm) == mv.OR.arr(ip,jm) == mv.OR.    &
            arr(im,j ) == mv.OR.arr(i,j ) == mv.OR.arr(ip,j ) == mv.OR. &
            arr(im,jp) == mv.OR.arr(i,jp) == mv.OR.arr(ip,jp) == mv)THEN
        tem1(i,j)=mv
      END IF

    END DO
  END DO

  DO j = jbgn,jend
    DO i = ibgn,iend
      arr(i,j) = tem1(i,j)
    END DO
  END DO

  RETURN
END SUBROUTINE smooth9pmv

SUBROUTINE buoycy_plt(nx,ny,nz,ptprt,pprt,qv,qscalar,            &
           ptbar,pbar,rhobar,qvbar, wbuoy, tem1)
!
!-----------------------------------------------------------------------
!
!     PURPOSE:
!
!     Calculate the total buoyancy including liquid and solid water
!     loading.
!
!-----------------------------------------------------------------------
!
!     AUTHOR: Ming Xue
!     10/10/91.
!
!     MODIFICATION HISTORY:
!
!     5/05/92 (M. Xue)
!     Added full documentation.
!
!     3/10/93 (M. Xue)
!     The buoyancy term is reformulated. The previous formula was
!     in error. The water loading was calculated wrong, resulting in
!     a value of the water loading that is typically an order of
!     magnitude too small.
!
!     3/25/94 (G. Bassett & M. Xue)
!     The buoyancy terms are reformulated for better numerical accuracy.
!     Instead of storing numbers which had the form (1+eps)*(1+eps1)
!     (eps << 1 and eps1 <<1), terms were expanded out, and most of the
!     high order terms neglected, except for the second order terms
!     in ptprt, pprt and qvbar.
!
!     9/10/94 (D. Weber & Y. Lu)
!     Cleaned up documentation.
!
!     6/21/95 (Alan Shapiro)
!     Fixed bug involving missing qvpert term in buoyancy formulation.
!
!     10/15/97 (Donghai Wang)
!     Added a new option for including the second order terms.
!
!     11/05/97 (D. Weber)
!     Changed lower loop bounds in DO LOOP 400 for computing the
!     buoyancy term from k=3,nz-2 to k=2,nz-1.  Level k=2 data will be
!     used in the hydrostatic pprt lower boundary condition (removed
!     DO LOOP 410 used to set wbuoy = 0.0 for k= 2 and nz-1).
!
!-----------------------------------------------------------------------
!
!     INPUT :
!
!       nx       Number of grid points in the x-direction (east/west)
!       ny       Number of grid points in the y-direction (north/south)
!       nz       Number of grid points in the vertical direction.
!
!       ptprt    Perturbation potential temperature at a time level (K)
!       pprt     Perturbation pressure at a given time level (Pascal)
!       qv       Water vapor specific humidity at a given time level
!                (kg/kg)
!       qc       Cloud water mixing ratio at a given time level (kg/kg)
!       qr       Rainwater mixing ratio at a given time level (kg/kg)
!       qi       Cloud ice mixing ratio at a given time level (kg/kg)
!       qs       Snow mixing ratio at a given time level (kg/kg)
!       qh       Hail mixing ratio at a given time level (kg/kg)
!
!       ptbar    Base state potential temperature (K)
!       pbar     Base state pressure (Pascal)
!       rhobar   Base state density rhobar
!       qvbar    Base state water vapor specific humidity (kg/kg)
!
!     OUTPUT:
!
!       wbuoy    The total buoyancy force (kg/(m*s)**2)
!
!     WORK ARRAYS:
!
!       tem1     Temporary work array.
!
!-----------------------------------------------------------------------
!
!     Variable Declarations
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!     Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'      ! Physical constants
  INCLUDE 'globcst.inc'

!-----------------------------------------------------------------------

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions

  REAL :: ptprt (nx,ny,nz)     ! Perturbation potential temperature
                               ! at a given time level (K)
  REAL :: pprt  (nx,ny,nz)     ! Perturbation pressure at a given time
                               ! level (Pascal)
  REAL :: qv    (nx,ny,nz)     ! Water vapor specific humidity (kg/kg)
  REAL :: qscalar(nx,ny,nz,nscalar)

  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal).
  REAL :: rhobar(nx,ny,nz)     ! Base state density rhobar
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific
                               ! humidity(kg/kg)

  REAL :: wbuoy(nx,ny,nz)      ! Total buoyancy in w-eq. (kg/(m*s)**2)

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array.

!
!-----------------------------------------------------------------------
!
!     Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  REAL    :: g5
  REAL    :: pttem,tema
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!     Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
!-----------------------------------------------------------------------
!
!     The total buoyancy
!
!       wbuoy = rhobar*g ( ptprt/ptbar-pprt/(rhobar*csndsq)+
!       qvprt/(rddrv+qvbar)-(qvprt+qc+qr+qs+qi+qh)/(1+qvbar)
!       -(ptprt*ptprt)/(ptbar*ptbar)                        !2nd-order
!       +0.5*(ptprt*pprt)/(cpdcv*ptbar*pbar))               !2nd-order
!
!     and rddrv=rd/rv, cp, cv, rd and rv are defined in phycst.inc.
!
!     Here, the contribution from pprt (i.e., term pprt/(rhobar*csndsq))
!     is evaluated inside the small time steps, therefore wbuoy
!     does not include this part.
!
!     The contribution from ptprt is calculated inside the small time
!     steps if the potential temperature equation is solved inside
!     small time steps, i.e., if ptsmlstp=1.
!
!-----------------------------------------------------------------------
!
  tema = 1.0/cpdcv
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        pttem = ptprt(i,j,k)/ptbar(i,j,k)
        tem1(i,j,k) = pttem*                                            &
            (1.0-pttem+0.5*pprt(i,j,k)*(tema/pbar(i,j,k)))
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!     Add on the contributions to the buoyancy from the water vapor
!     content and the liquid and ice water loading.
!
!-----------------------------------------------------------------------
!

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem1(i,j,k) = tem1(i,j,k)                                       &
            + (qv(i,j,k) - qvbar(i,j,k))/(rddrv + qvbar(i,j,k))         &
            - (qv(i,j,k) - qvbar(i,j,k))/(1 + qvbar(i,j,k))
!            - (qv(i,j,k) - qvbar(i,j,k) + qc(i,j,k) + qr(i,j,k) +       &
!            qs(i,j,k) + qi(i,j,k) + qh(i,j,k))/(1 + qvbar(i,j,k))

        tema = 1 + qvbar(i,j,k)
        IF (P_QC > 0) tem1(i,j,k) = tem1(i,j,k) - qscalar(i,j,k,P_QC)/tema
        IF (P_QR > 0) tem1(i,j,k) = tem1(i,j,k) - qscalar(i,j,k,P_QR)/tema
        IF (P_QI > 0) tem1(i,j,k) = tem1(i,j,k) - qscalar(i,j,k,P_QI)/tema
        IF (P_QS > 0) tem1(i,j,k) = tem1(i,j,k) - qscalar(i,j,k,P_QS)/tema
        IF (P_QG > 0) tem1(i,j,k) = tem1(i,j,k) - qscalar(i,j,k,P_QG)/tema
        IF (P_QH > 0) tem1(i,j,k) = tem1(i,j,k) - qscalar(i,j,k,P_QH)/tema

      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!     Then the total buoyancy:
!
!       wbuoy = tem1 * rhobar * g
!
!     averged to the w-point on the staggered grid.
!
!-----------------------------------------------------------------------
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        wbuoy(i,j,k)= tem1(i,j, k )*rhobar(i,j, k )*g
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE buoycy_plt

!#######################################################################
!
!   TO DETERMINE THE CONTOUR INCRMENT AND CONTOUR VALUES FOR a
!   max and min from a data set.
!
!   It is suppose that only ROOT processor call this subroutine.
!
!#######################################################################

SUBROUTINE setcontr(zmin,zmax,nminctr,nmaxctr,cl,ncl)

  IMPLICIT NONE

  REAL,    INTENT(IN)    :: zmin, zmax       ! field min, and max
  INTEGER, INTENT(IN)    :: nminctr, nmaxctr ! contour number limits
  REAL,    INTENT(INOUT) :: cl(*)            ! contour levels
  INTEGER, INTENT(OUT)   :: ncl              ! actual contour number

  INTEGER :: nch
  COMMON /XOUTCH/ nch
  INTEGER :: LCPTN,LABTYP,ICLF,LHILIT,IHLF,KCT01
  REAL    :: clref
  COMMON /XCRF17/clref,LCPTN,LABTYP,ICLF,LHILIT,IHLF,KCT01

  REAL    :: zminc, zmaxc     ! first and last contours

  REAL    :: diff, zinc, clv
  REAL    :: eps, zref
  INTEGER :: kcount, kzinc

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  zinc = cl(2) - cl(1)
  diff = zmax-zmin
  IF(diff <= ABS(zinc)*1.0E-6) THEN
     WRITE(nch,'(1x,a,/1x,a)')                                      &
       'Bad first guess of contour increment or field is constant', &
       'number of contours is one.'
     ncl = 1
     cl(1) = zmin
     cl(2) = zmin + 1.0
     RETURN
  ENDIF

  zref = clref
  IF (abs(cl(1)-clref) < abs(zinc)) zref = cl(1)

  kcount = 0
  1 CONTINUE
  eps    = 0.001*zinc
  kcount = kcount + 1
  IF (kcount > 20) GOTO 998    ! too many loops, abort the program
  kzinc = (zmin-zref)/zinc
  zminc = kzinc*zinc + zref
  kzinc = (zmax-zref)/zinc
  zmaxc = kzinc*zinc + zref
!   IF(zmin-clref > 0.0) ZMINC=ZMINC+ZINC
!   IF(zmax-clref < 0.0) ZMAXC=ZMAXC-ZINC

  CLV = ZMINC-ZINC
  NCL = 0
  6 CLV = CLV + ZINC
  IF(CLV-ZMAXC-EPS > 0.0) GOTO 8  ! Reach zmax, check the contour levels
  NCL = NCL + 1

  IF(NCL > nmaxctr) THEN
    ZINC=ZINC*2
    WRITE(nch,'(a,I3,a,E10.3)')           &
      ' Number of contours > ',nmaxctr,' ,Zinc is doubled. Zinc=',zinc
    GO TO 1
  ENDIF
  IF( ABS( CLV-zref ) < EPS ) CLV=zref
  CL(NCL) = CLV
!  WRITE(6,*) 'ncl = ',ncl,', cl(ncl) = ',clv
  GOTO 6

  8 CONTINUE

  IF( NCL < nminctr) THEN
    ZINC=ZINC/2
    WRITE(nch,'(a,I3,a,E10.3)')                  &
       ' Number of contours < ',nminctr,' ,Zinc is halved. Zinc=',zinc
    GO TO 1
  ENDIF

  WRITE(nch,'(a,I5,2(a,E12.4),a,E12.5)')                         &
      ' * Number of contours= ',ncl,                             &
      '  Min=',zminc, ' Max=', zmaxc,' Inc=',zinc

  RETURN

  998 CONTINUE
  WRITE(NCH,*)' Contour levels can not be selected by XCNTLV.'
  WRITE(NCH,*)                                                &
    ' Plz alter input contour interval or limits of contour number'

   RETURN
END SUBROUTINE setcontr

SUBROUTINE plt_trajc(slicopt, curtim, x01,x02,y01,y02,x1, lblmag, iorig,&
                     xor_current,yor_current,xgrdorg,ygrdorg,xorig,yorig)

  IMPLICIT NONE
  INTEGER :: slicopt
  REAL :: curtim
  REAL :: x01,x02,y01,y02,x1,lblmag
  INTEGER, INTENT(IN) :: iorig
  REAL,    INTENT(IN) :: xor_current,yor_current
  REAL,    INTENT(IN) :: xgrdorg, ygrdorg
  REAL,    INTENT(IN) :: xorig, yorig

  INCLUDE 'arpstrajc.inc'

  CHARACTER(LEN=25) :: timestring

  REAL :: xpnt, ypnt

  INTEGER :: npoints_cur,itrajc1,jtrajc1
  INTEGER :: istat
  INTEGER :: ireturn

! REAL :: tstart,tzero, tend
! REAL :: xlow, xhigh, ylow, yhigh, zlow, zhigh, pttem
! REAL :: dx,dy,dz

  INTEGER :: pen_status
  REAL :: pl1, pr1, pb1, pt1, trajc_lbl_mag, trajc_mkr_mag
  REAL :: dist
  REAL :: xl, xr, yb, yt, vshift, hshift

  INTEGER :: i,j,k, itrajc_start, itrajc_end
  CHARACTER(LEN=100) :: trajcfn_header

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF( trajc_plt_opt == 0 ) RETURN

  WRITE(6,'(/1x,a,i2,a,f8.2,a/)')                                         &
    'Plotting trajectories for slicopt= ', slicopt,', at time ',curtim,' seconds.'

  CALL xqmap( xl, xr, yb, yt)
  CALL xwindw(xl, xr, yb, yt)

  DO k=1,ntimes

!-----------------------------------------------------------------------
! Start plotting trajectories
!-----------------------------------------------------------------------

    IF( trajc_plt_bgn_time == 0.0 .AND. trajc_plt_end_time == 0.0 ) THEN
      trajc_plotting_startime = ttrajc(1)
      trajc_plotting_stoptime = ttrajc(npoints_in(k))
    ELSE

      IF( ABS( trajc_plt_bgn_time + 9999.0 ) < 0.001 ) THEN
        trajc_plotting_startime = curtim
      ELSE
        trajc_plotting_startime = trajc_plt_bgn_time
      ENDIF

      IF( ABS( trajc_plt_end_time + 9999.0 )< 0.001 ) THEN
        trajc_plotting_stoptime = curtim
      ELSE
        trajc_plotting_stoptime = trajc_plt_end_time
      ENDIF

    ENDIF

    CALL xthick(2)

    IF(trajc_lbl_opt /= 0) THEN
      CALL xqpspc( pl1, pr1, pb1, pt1)
      trajc_lbl_mag = trajc_lbl_siz*MIN(pt1-pb1, pr1-pl1)*lblmag
      CALL xchmag(trajc_lbl_mag)
    END IF
    hshift = (xr-xl)*trajc_lbl_mag*0.5
    vshift = (yt-yb)*trajc_lbl_mag*0.4

    IF(trajc_mkr_typ /= 0) THEN
      CALL xqpspc( pl1, pr1, pb1, pt1)
      trajc_mkr_mag = trajc_mkr_siz*MIN(pt1-pb1, pr1-pl1)
    END IF

    itrajc_start = min(ntrajc_start, ntrajcs(k))
    itrajc_end   = min(ntrajc_end  , ntrajcs(k))
    IF( ntrajc_end < 0 ) itrajc_end = ntrajcs(k)

    DO i=itrajc_start, itrajc_end, ntrajc_stride

      CALL xcolor(traj_col(i))

      IF(slicopt==1 .OR. slicopt==4 .OR. slicopt==6 .OR. slicopt==7) THEN ! horizontal slice

! plot only the trajectory corresponding to the current slice when constant height level
! slices follow trajectory
!
        IF((.not.slicopt==4) .OR. (h_follow_trajc==0) .OR.                &
           (slicopt==4 .AND. i==itrajc_index) ) THEN

          pen_status = 0

          DO j=1,npoints_in(k)

            IF( ttrajc(j) < trajc_plotting_startime-0.001 ) CYCLE
            IF( ttrajc(j) > trajc_plotting_stoptime+0.001 ) EXIT

            IF( ( abs( xtrajc(j,i,k) + 99999.0 ) < 0.001 .OR.               &
                  abs( ytrajc(j,i,k) + 99999.0 ) < 0.001  ) ) THEN

              pen_status = 0

            ELSE

              IF( pen_status == 0 ) then
                CALL XPENUP(xtrajc(j,i,k)*0.001,ytrajc(j,i,k)*0.001)
              ELSE
                CALL XPENDN(xtrajc(j,i,k)*0.001,ytrajc(j,i,k)*0.001)
              ENDIF
              pen_status = 1

            ENDIF

          END DO

          DO j=1,npoints_in(k)

            IF( ttrajc(j) < trajc_plotting_startime-0.001 ) CYCLE
            IF( ttrajc(j) > trajc_plotting_stoptime+0.001 ) EXIT

            IF(  ( abs( xtrajc(j,i,k) + 99999.0 ) .lt. 0.001 .or.  &
                   abs( ytrajc(j,i,k) + 99999.0 ) .lt. 0.001  ) ) CYCLE

            IF(trajc_lbl_opt /= 0 .AND. MOD(j-1,trajc_lbl_frq) == 0) THEN

              IF(trajc_lbl_opt == 3) then
                CALL get_time_string ( ttrajc(j), timestring,'Z  ',0 )
                CALL XCHARL(xtrajc(j,i,k)*0.001+hshift,ytrajc(j,i,k)*0.001-vshift,timestring(1:6))
              ELSE

                IF(trajc_lbl_opt == 1) then
                  trajc_lbl_number = ztrajc(j,i,k)*0.001 ! Z MSL in km
                ELSEIF(trajc_lbl_opt == 2) then
                  trajc_lbl_number = ttrajc(j)/60.0 ! time label in min
                ENDIF

                CALL label_trajc(xtrajc(j,i,k)*0.001+hshift,ytrajc(j,i,k)*0.001-vshift, &
                                 trajc_lbl_number,trajc_lbl_fmt)
              ENDIF

            END IF

            IF((trajc_mkr_typ /= 0) .and. (MOD(j-1,trajc_mkr_frq) == 0 &
                .or.abs( ttrajc(j)-curtim ) < 0.01)) THEN
              IF( abs( ttrajc(j)-curtim ) < 0.01 ) then
                CALL xmrksz(trajc_mkr_mag*3)
              ELSE
                CALL xmrksz(trajc_mkr_mag)
              ENDIF
              CALL xmarker(xtrajc(j,i,k)*0.001,ytrajc(j,i,k)*0.001,trajc_mkr_typ)
            ENDIF

          END DO

        END IF !

      ELSE IF(slicopt == 2) THEN ! y-z slice

        pen_status = 0

        DO j=1,npoints_in(k)

          IF( ttrajc(j) < trajc_plotting_startime-0.001 ) CYCLE
          IF( ttrajc(j) > trajc_plotting_stoptime+0.001 ) EXIT

          IF( ( abs( xtrajc(j,i,k) + 99999.0 ) .lt. 0.001 .or.  &
                abs( ztrajc(j,i,k) + 99999.0 ) .lt. 0.001  ) ) then
            pen_status = 0
          ELSE

            IF( pen_status == 0 ) then
              CALL XPENUP(xtrajc(j,i,k)*0.001,ztrajc(j,i,k)*0.001)
            ELSE
              CALL XPENDN(xtrajc(j,i,k)*0.001,ztrajc(j,i,k)*0.001)
            ENDIF
            pen_status = 1

          ENDIF

        END DO

        DO j=1,npoints_in(k)

          IF( ttrajc(j) < trajc_plotting_startime-0.001 ) CYCLE
          IF( ttrajc(j) > trajc_plotting_stoptime+0.001 ) EXIT

          IF( ( abs( xtrajc(j,i,k) + 99999.0 ) .lt. 0.001 .or.  &
                abs( ztrajc(j,i,k) + 99999.0 ) .lt. 0.001  ) ) CYCLE

          IF(trajc_lbl_opt /= 0 .AND. MOD(j-1,trajc_lbl_frq) == 0) THEN

            IF(trajc_lbl_opt == 3) then
              CALL get_time_string ( ttrajc(j), timestring,'Z  ',0 )
              CALL XCHARL(xtrajc(j,i,k)*0.001+hshift,ztrajc(j,i,k)*0.001-vshift, timestring(1:6))
            ELSE

              IF(trajc_lbl_opt == 1) THEN
                trajc_lbl_number = ytrajc(j,i,k)*0.001 ! Z MSL in km
              ELSE IF(trajc_lbl_opt == 2) THEN
                trajc_lbl_number = ttrajc(j)/60.0 ! time label in min
              END IF

              CALL label_trajc(xtrajc(j,i,k)*0.001+hshift,ztrajc(j,i,k)*0.001-vshift, &
                               trajc_lbl_number,trajc_lbl_fmt)
            END IF

          END IF

          IF((trajc_mkr_typ /= 0) .and. (MOD(j-1,trajc_mkr_frq) == 0 &
              .or.abs( ttrajc(j)-curtim ) < 0.01)) THEN
            IF( abs( ttrajc(j)-curtim ) < 0.01 ) then
              CALL xmrksz(trajc_mkr_mag*3)
            ELSE
              CALL xmrksz(trajc_mkr_mag)
            ENDIF
            CALL xmarker(xtrajc(j,i,k)*0.001,ztrajc(j,i,k)*0.001,trajc_mkr_typ)
          ENDIF

        END DO

      ELSE IF(slicopt == 3) THEN

        pen_status = 0

        DO j=1,npoints_in(k)

          IF( ttrajc(j) < trajc_plotting_startime-0.001 ) CYCLE
          IF( ttrajc(j) > trajc_plotting_stoptime+0.001 ) EXIT

          IF( ( abs( ytrajc(j,i,k) + 99999.0 ) .lt. 0.001 .or.  &
                abs( ztrajc(j,i,k) + 99999.0 ) .lt. 0.001  ) ) then
            pen_status = 0
          ELSE

            IF( pen_status == 0 ) then
              CALL XPENUP(ytrajc(j,i,k)*0.001,ztrajc(j,i,k)*0.001)
            ELSE
              CALL XPENDN(ytrajc(j,i,k)*0.001,ztrajc(j,i,k)*0.001)
            ENDIF
            pen_status = 1

          ENDIF
        END DO

        DO j=1,npoints_in(k)

          IF( ttrajc(j) < trajc_plotting_startime-0.001 ) CYCLE
          IF( ttrajc(j) > trajc_plotting_stoptime+0.001 ) EXIT

          IF( ( abs( ytrajc(j,i,k) + 99999.0 ) .lt. 0.001 .or.  &
                abs( ztrajc(j,i,k) + 99999.0 ) .lt. 0.001  ) ) CYCLE

          IF(trajc_lbl_opt /= 0 .AND. MOD(j-1,trajc_lbl_frq) == 0) THEN

            IF(trajc_lbl_opt == 3) then
              CALL get_time_string ( ttrajc(j), timestring,'Z  ',0 )
              CALL XCHARL(ytrajc(j,i,k)*0.001+hshift,ztrajc(j,i,k)*0.001-vshift, timestring(1:6))
            ELSE

              IF(trajc_lbl_opt == 1) then
                trajc_lbl_number = xtrajc(j,i,k)*0.001 ! Z MSL in km
              ELSEIF(trajc_lbl_opt == 2) then
                trajc_lbl_number = ttrajc(j)/60.0 ! time label in min
              ENDIF

              CALL label_trajc(ytrajc(j,i,k)*0.001+hshift,ztrajc(j,i,k)*0.001-vshift, &
                               trajc_lbl_number,trajc_lbl_fmt)
            ENDIF

          END IF

          IF((trajc_mkr_typ /= 0) .and. (MOD(j-1,trajc_mkr_frq) == 0 &
              .or.abs( ttrajc(j)-curtim ) < 0.01)) THEN
            IF( abs( ttrajc(j)-curtim ) < 0.01 ) then
              CALL xmrksz(trajc_mkr_mag*3)
            ELSE
              CALL xmrksz(trajc_mkr_mag)
            ENDIF
            CALL xmarker(ytrajc(j,i,k)*0.001,ztrajc(j,i,k)*0.001,trajc_mkr_typ)
          ENDIF

        END DO

      ELSE IF(slicopt == 5) THEN  ! arbitary vertical cross section

        dist    = SQRT((x01-x02)**2+(y01-y02)**2)

        pen_status = 0

        DO j=1,npoints_in(k)

          IF( ttrajc(j) < trajc_plotting_startime-0.001 ) CYCLE
          IF( ttrajc(j) > trajc_plotting_stoptime+0.001 ) EXIT

          IF( ( abs( xtrajc(j,i,k) + 99999.0 ) .lt. 0.001 .or.  &
                abs( ytrajc(j,i,k) + 99999.0 ) .lt. 0.001 .or.  &
                abs( ztrajc(j,i,k) + 99999.0 ) .lt. 0.001   ) ) then
            pen_status = 0
          ELSE

            xpnt =  x1+( (xtrajc(j,i,k)*0.001-x01)*(x02-x01)+   &
                         (ytrajc(j,i,k)*0.001-y01)*(y02-y01))/dist

            IF( pen_status == 0 ) then
              CALL XPENUP(xpnt,ztrajc(j,i,k)*0.001)
            ELSE
              CALL XPENDN(xpnt,ztrajc(j,i,k)*0.001)
            ENDIF
            pen_status = 1

!           print*,'j,xpnt,ztrajc(j,i,k)*0.001=',j,xpnt,ztrajc(j,i,k)*0.001

          ENDIF
        END DO

        DO j=1,npoints_in(k)

          IF( ttrajc(j) < trajc_plotting_startime-0.001 ) CYCLE
          IF( ttrajc(j) > trajc_plotting_stoptime+0.001 ) EXIT

          IF( ( abs( xtrajc(j,i,k) + 99999.0 ) .lt. 0.001 .or.  &
                abs( ytrajc(j,i,k) + 99999.0 ) .lt. 0.001 .or.  &
                abs( ztrajc(j,i,k) + 99999.0 ) .lt. 0.001   ) ) CYCLE

          xpnt =  x1+( (xtrajc(j,i,k)*0.001-x01)*(x02-x01)+   &
                       (ytrajc(j,i,k)*0.001-y01)*(y02-y01))/dist

          IF((trajc_lbl_opt /= 0) .and. MOD(j-1,trajc_lbl_frq) == 0) THEN

            IF(trajc_lbl_opt == 3) then
              CALL get_time_string ( ttrajc(j), timestring,'Z  ',0 )
              CALL XCHARL(xpnt+hshift,ztrajc(j,i,k)*0.001-vshift, timestring(1:6))
            ELSE

              IF(trajc_lbl_opt == 1) then
                ypnt = (-(xtrajc(j,i,k)*0.001-x01)*(y02-y01)+   &
                         (ytrajc(j,i,k)*0.001-y01)*(x02-x01))/dist
                trajc_lbl_number = ypnt
              ELSEIF(trajc_lbl_opt == 2) then
                trajc_lbl_number = ttrajc(j)/60.0 ! time label in min
              ENDIF

              CALL label_trajc(xpnt+hshift,ztrajc(j,i,k)*0.001-vshift, &
                               trajc_lbl_number,trajc_lbl_fmt)
            ENDIF

          END IF

          IF((trajc_mkr_typ /= 0) .and. (MOD(j-1,trajc_mkr_frq) == 0 &
              .or.abs( ttrajc(j)-curtim ) < 0.01)) THEN
            IF( abs( ttrajc(j)-curtim ) < 0.01 ) then
              CALL xmrksz(trajc_mkr_mag*3)
            ELSE
              CALL xmrksz(trajc_mkr_mag)
            ENDIF
            CALL xmarker(xpnt,ztrajc(j,i,k)*0.001,trajc_mkr_typ)
          ENDIF

        END DO

      ELSE IF(slicopt == 12) THEN  ! vertical cross section along trajectory

        !print*,'npoints_bgn(i,k),npoints_end(i,k)=', npoints_bgn(i,k),npoints_end(i,k)

        IF( i == itrajc_index ) THEN ! plot only the trajectory corresponding to the current slice

          DO j=npoints_bgn(i,k),npoints_end(i,k)

            IF( j == npoints_bgn(i,k) ) then
              xpnt = 0.0
              CALL XPENUP(xpnt,ztrajc(j,i,k)*0.001)
            ELSE
              xpnt=xpnt+sqrt((xtrajc(j,i,k)-xtrajc(j-1,i,k))**2+(ytrajc(j,i,k)-ytrajc(j-1,i,k))**2)*0.001
              CALL XPENDN(xpnt,ztrajc(j,i,k)*0.001)
            ENDIF


          END DO

          DO j=npoints_bgn(i,k),npoints_end(i,k)

            IF( j == npoints_bgn(i,k) ) then
              xpnt = 0.0
            ELSE
              xpnt=xpnt+sqrt((xtrajc(j,i,k)-xtrajc(j-1,i,k))**2+(ytrajc(j,i,k)-ytrajc(j-1,i,k))**2)*0.001
            ENDIF

            IF((trajc_lbl_opt /= 0) .and. MOD(j-1,trajc_lbl_frq) == 0) THEN

              IF(trajc_lbl_opt == 3) then
                CALL get_time_string ( ttrajc(j), timestring,'Z  ',0 )
                CALL XCHARL(xpnt+hshift,ztrajc(j,i,k)*0.001-vshift, timestring(1:6))
              ELSEIF(trajc_lbl_opt == 2) then
                trajc_lbl_number = ttrajc(j)/60.0 ! time label in min
                CALL label_trajc(xpnt+hshift,ztrajc(j,i,k)*0.001-vshift, &
                                 trajc_lbl_number,trajc_lbl_fmt)
              ENDIF

            END IF

            IF((trajc_mkr_typ /= 0) .and. (MOD(j-1,trajc_mkr_frq) == 0  &
                .OR. abs( ttrajc(j)-curtim ) < 0.01)) THEN
              IF( abs( ttrajc(j)-curtim ) < 0.01 ) then
                CALL xmrksz(trajc_mkr_mag*3)
              ELSE
                CALL xmrksz(trajc_mkr_mag)
              ENDIF
              CALL xmarker(xpnt,ztrajc(j,i,k)*0.001,trajc_mkr_typ)
            ENDIF

          END DO

        ENDIF

      END IF ! slicopt

    END DO ! i=1,ntrajcs

  END DO ! k=1,ntimes

  CALL xthick(1)
  CALL xfull
  CALL xwdwof

  RETURN
END SUBROUTINE plt_trajc

SUBROUTINE label_trajc(x,y, trajc_lbl_number, trajc_lbl_fmt )
!
! Label points along a trajectory
!
  REAL :: x,y,trajc_lbl_number
  INTEGER :: trajc_lbl_fmt, lch
  CHARACTER (LEN=6) :: trajc_lblfmt
  CHARACTER (LEN=1) :: stem1
  CHARACTER (LEN=10):: trajc_lbl_string

  IF(trajc_lbl_fmt < 0 ) then
    CALL XRNUMB(x,y,trajc_lbl_number,'*')
  ELSE IF(trajc_lbl_fmt == 0 ) then
    CALL XINUMB(x,y,NINT(trajc_lbl_number),'*')
  ELSE
    WRITE(stem1,'(i1)') trajc_lbl_fmt
    WRITE(trajc_lblfmt,'(a,a1,a)') '(f8.',stem1,')'
    write(trajc_lbl_string,trajc_lblfmt) trajc_lbl_number
    LCH=20
    CALL XCHLJ( trajc_lbl_string, LCH)
    CALL XCHARL(x,y, trim(trajc_lbl_string(1:lch)) )
  ENDIF

  RETURN
END SUBROUTINE label_trajc

SUBROUTINE read_trajc(ireturn)

  IMPLICIT NONE

  INTEGER :: k
  INTEGER :: ireturn

  INCLUDE 'arpstrajc.inc'

  INTEGER :: ntrajcs_in,npoints, npoints_cur,itrajc1,jtrajc1
  INTEGER :: istat, j, i, nunit

  REAL :: tstart,tzero, tend
  REAL :: dx, dy, dz
  REAL :: xlow, xhigh, ylow, yhigh, zlow, zhigh, pttem
  CHARACTER(LEN=100) :: trajcfn_header

  ireturn = 0

  DO k=1,ntimes

!-----------------------------------------------------------------------
! Read kth group of trajectories
!-----------------------------------------------------------------------

    CALL getunit (nunit)

    WRITE(6,'(1x,a,a)') 'Opening ',trim(trajc_fn_in(k))
    OPEN(UNIT=nunit,FILE=trim(trajc_fn_in(k)),STATUS='old',             &
         FORM='formatted',IOSTAT= istat )

    IF(istat == 0) THEN
!      print*,'Trajectory file ',trim(trajc_fn_in(k)),' successfully openned.'

      READ(nunit,'(a)') trajcfn_header
      READ(nunit,'(6e17.6)') xlow, xhigh, ylow, yhigh, zlow, zhigh

!     write(6,'(6e17.6)') xlow, xhigh, ylow, yhigh, zlow, zhigh
      READ(nunit,'(3e17.6)') dx, dy, dz
!     write(6,'(3e17.6)') dx, dy, dz

      READ(nunit,'(3e17.6)') tstart, tzero, tend
      READ(nunit,'(i10)') npoints
      READ(nunit,'(i10)') ntrajcs(k)

      DO j=1,npoints
        READ(nunit,'(4e17.6)',err=115,end=115) ttrajc(j)
        READ(nunit,'(i10)',err=115,end=115) ntrajcs_in
        IF( ntrajcs_in /= ntrajcs(k) ) then
          PRINT*,'ntrajcs read in .ne. ntrajcs in program.'
          PRINT*,'Job stopped'
          CALL arpsstop('Reading error for trajc data.',1)
        END IF
        READ(nunit,'(6e17.6)',err=115,end=115)                          &
           (xtrajc(j,i,k),ytrajc(j,i,k),ztrajc(j,i,k),i=1,ntrajcs(k))
      END DO

      CLOSE(UNIT=nunit)
      CALL retunit(nunit)

      WRITE(6,'(/1x,a,a,a/)') 'Trajectory file ',trim(trajc_fn_in(k)),' successfully read.'

    ELSE
      CALL retunit(nunit)
      135 WRITE(6,'(1x,a)') 'Failed to open trajectory data file.'
      ireturn = 1
      RETURN
    END IF

    npoints_in(k) = npoints

    GOTO 125
    115 continue

    npoints_in(k) = max(1,j-1)

    125 continue

  ENDDO ! k

  RETURN

END SUBROUTINE read_trajc

SUBROUTINE select_trajc_points(x1,x2,y1,y2,time)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Determining the starting and ending points of each trajectory to be
!  plotted.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  REAL :: x1,x2,y1,y2 ! The x and y domain bounds of plotting space.
  REAL :: time ! time of the data being plotted

  REAL :: xptr1,xptr2,dxptr,xtr,ytr
  INTEGER :: itrajc_start, itrajc_end, iptr
  INTEGER :: ireturn,j_bgn,j_end
  INTEGER :: j_bgn_last,j_end_last,istatus,kz
  INTEGER :: i,j,k

!----------------------------------------------------------------------
!
! Include files
!
!---------------------------------------------------------------------
  INCLUDE 'mp.inc'
  INCLUDE 'arpstrajc.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

    DO k=1,ntimes

      IF( trajc_plt_bgn_time == 0.0 .and. trajc_plt_end_time == 0.0 ) THEN
        trajc_plotting_startime = ttrajc(1)
        trajc_plotting_stoptime = ttrajc(npoints_in(k))
      ELSE

        IF( ABS( trajc_plt_bgn_time + 9999.0 )< 0.001 ) THEN
          trajc_plotting_startime = time
        ELSE
          trajc_plotting_startime = trajc_plt_bgn_time
        ENDIF

        IF( ABS( trajc_plt_end_time + 9999.0 )< 0.001 ) THEN
          trajc_plotting_stoptime = time
        ELSE
          trajc_plotting_stoptime = trajc_plt_end_time
        ENDIF

      ENDIF

      itrajc_start = min(ntrajc_start, ntrajcs(k))
      itrajc_end   = min(ntrajc_end  , ntrajcs(k))
      IF( ntrajc_end < 0 ) itrajc_end = ntrajcs(k)

      DO i=itrajc_start, itrajc_end, ntrajc_stride

        j_end  = npoints_in(k)
        j_end_last  = npoints_in(k)

        DO j=npoints_in(k),1,-1

          xtr = xtrajc(j,i,k)
          ytr = ytrajc(j,i,k)
          IF(abs(xtr+99999.0)<0.0001.or.abs(ytr+99999.0)<0.0001.or. &
             abs(ztrajc(j,i,k)+99999.0)<0.0001) CYCLE

          xtr = xtr*0.001
          ytr = ytr*0.001

          IF(.not.( xtr<x1.or.xtr>x2.or.ytr<y1.or.ytr>y2 ) .and. &
              ttrajc(j) < trajc_plotting_stoptime+0.0010 ) then
            j_end=j
            exit
          ELSE
            j_end_last = j
          ENDIF

        ENDDO

!       IF( j_end /= j_end_last ) then
!         CALL clipwd(x01,y01,x02,y02, idisplay )
!       ENDIF

        j_bgn = 1
        j_bgn_last = 1
        DO j=1,npoints_in(k)

          xtr = xtrajc(j,i,k)
          ytr = ytrajc(j,i,k)
          IF(abs(xtr+99999.0)<0.0001.or.abs(ytr+99999.0)<0.0001.or. &
             abs(ztrajc(j,i,k)+99999.0)<0.0001) CYCLE

          xtr = xtr*0.001
          ytr = ytr*0.001

          IF(.not.( xtr<x1.or.xtr>x2.or.ytr<y1.or.ytr>y2 ) .and. &
              ttrajc(j) > trajc_plotting_startime+0.0010 ) then
            j_bgn=j
            exit
          ELSE
            j_bgn_last = j
          ENDIF

        ENDDO

!       IF( j_bgn /= j_bgn_last ) then
!       CALL clipwd(x01,y01,x02,y02, idisplay )
!       ENDIF

        npoints_bgn(i,k)=j_bgn
        npoints_end(i,k)=j_end

!       print*,'i,k,npoints_bgn(i,k),npoints_end(i,k)=',i,k,npoints_bgn(i,k),npoints_end(i,k)

      ENDDO ! i

    ENDDO ! k

RETURN

END SUBROUTINE select_trajc_points
