!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE INITPARA                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE initpara(nx,ny,nz,nzsoil,nstyps,namelist_filename)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Initialize the model control parameters. Most of them are read in
!  from an input file.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!
!  3/17/1991.
!
!  MODIFICATION HISTORY:
!
!  5/02/92 (M. Xue)
!    Added full documentation.
!
!  5/25/92 (M. Xue)
!    Reworked to provide a friendly user interface for control
!    parameter input, and write out a log file.
!
!  6/04/92 (M. Xue)
!    Further facelift.
!
!  8/03/92 (M. Xue)
!    The grid scale in the formula of divergecne damping coefficient
!    is changed from (dx*dy*dz)**(2/3) to min(dx,dy,dz)**2.
!    This will affect the results as compared with previous runs.
!    Added output control parameter input etc.
!
!  4/16/93 (M. Xue and H. Jin)
!    Added parameter inputs that are related to terrain.
!
!  9/20/93 (A. Sathye and M. Xue)
!    Changed to the NAMELIST input format.
!
!  9/27/93 (M. Xue)
!    For non-stretched case, dzmin is set to dz. For 2-d mode,
!    appropriate LBC's are automatially set to periodic B.C.
!
!  12/3/93 (M. Xue)
!    Added parameters for automatical grid translation.
!
!  2/12/94 (Yuhe Liu)
!    Added parameters for surface energy budget model.
!
!  10/26/94 (Y. Liu)
!    Add lbcopt to namelist &boundary_condition_options for the
!    lateral boundary condition option.
!
!  12/12/94 (Y. Liu)
!    Changed default values of variables in namelists to the same as
!    in User's Guide and corrected the log file output.
!
!  12/22/94 (Y. Liu)
!    Added more parameters into the namelist blocks, including
!    ubar0 and vbar0 in initialization, and strhtune in grid.
!
!  01/28/95 (G. Bassett)
!    Added new parameter, buoyopt, to input namelist &initialization.
!
!  08/24/95 (K. Brewster)
!    Changed informative opening message to tell user about namelist.
!
!  2/2/96 (Donghai Wang & yuhe Liu)
!    Added parameters for map projection factor.
!
!  3/26/96 (Yuhe Liu)
!    Added a namelist, &radiation, and parameters for radiation.
!
!  4/2/96  (Donghai Wang, X. Song and M. Xue)
!    Added parameters for implicit treatment of vertical mixing.
!
!  5/7/96  (Donghai Wang and M. Xue)
!    Added a parameter for Rayleigh damping.
!
!  7/31/96 (Ming Xue and Yuhe Liu)
!    Added the isotropic option for divergence damping. Parameter
!    divdmpnd changed to divdmpndh for horizontal and divdmpndv for
!    vertical.
!
!  3/23/97 (Ming Xue)
!    Parameter scmixfctr added to namelist block computational_mixing.
!
!  3/23/97 (Ming Xue)
!    Modifications made so that the program will complete reading in
!    input parameters and check their validity even when error is
!    encountered before it stops at the end of this subroutine.
!
!  7/27/97 (Dan Weber)
!    Added fftopt to the list of specified parameters.
!
!  10/21/97 (Donghai Wang)
!    Added two parameters,buoy2nd and rhofctopt.
!
!  04/15/98 (Donghai Wang)
!    Added a new fraction factor for Kain-Fritsch scheme.
!
!  08/31/98 (K. Brewster)
!    Added nudging NAMELIST to version 4.4.
!
!  1999/10/21 (Gene Bassett)
!    Separated the reading in from the computation of derived variables
!    (since some values were not being set when reading was aborted due
!    to namelist read errors).
!
!  2000/04/13 (Gene Bassett)
!    Added dumping options for HDF formats.
!
!  2000/04/24 (Gene Bassett)
!    Update message passing version and added grid_dims &
!    message_passing namelist blocks.
!
!  2002/03/20 (Dan Weber, M. Xue, X. Jin)
!   Added one option for vertically implicit fall velocity and
!   weight coefficients for two time levels
!
!  05/14/2002 (J. Brotzge)
!    Modified code for multiple soil schemes
!
!  08/27/02 (Dan Weber)
!     Added option for using Ferrier (1994) fall velocity coefficients.
!     They replace the Lin constants when fallopt=2, otherwise fallopt
!     =1 for lin scheme values.  Added fallopt variable.
!
!  01/30/07 (Mingjing Tong)
!  Added option to add perturbation to soil model forcing term shflx
!  and lhflx
!
!  12/07/10 (Bryan Putnam)
!  Added option to turn off hail.
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction
!    ny       Number of grid points in the y-direction
!    nz       Number of grid points in the z-direction
!    nzsoil   Number of grid points in the soil profile
!
!  OUTPUT:
!
!    Control parameters defined in include files.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations. (Local Variables)
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: nx,ny,nz      ! The number of grid points in 3 dimensions.
  INTEGER, INTENT(OUT) :: nzsoil        ! The number of grid points in the soil
  INTEGER, INTENT(OUT) :: nstyps        ! Maximum number of soil types per grid point.

  CHARACTER(LEN=*), INTENT(INOUT) :: namelist_filename
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  REAL    :: wrmax            ! Maximun value of canopy moisture
  INTEGER :: i

  INTEGER :: lenstr      ! Length of a string
  LOGICAL :: iexist      ! Flag set by inquire statement for file
                         ! existence
  REAL    :: temr
  REAL    :: dtsml0,dtsfc0        ! Temporary variable

  CHARACTER (LEN=19)  :: initime  ! Real time in form of 'year-mo-dy:hr:mn:ss'

!  INTEGER, PARAMETER  :: unum=5   ! unit number for reading in namelist
  INTEGER :: unum         ! unit number for reading in namelist
  CHARACTER(LEN=256)   :: nlfile, outdirname

  INTEGER :: inisplited, dmp_out_joined
  INTEGER :: idummy

!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Global constants and parameters, most of them specify the
!  model run options.
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'soilcst.inc'
  INCLUDE 'nudging.inc'
  INCLUDE 'dfilter.inc'
  INCLUDE 'radcst.inc'
!
!-----------------------------------------------------------------------
!
!  Grid and map parameters.
!
!-----------------------------------------------------------------------
!
  INCLUDE 'grid.inc'
!
!-----------------------------------------------------------------------
!
!  Control parameters defining the boundary condition types.
!
!-----------------------------------------------------------------------
!
  INCLUDE 'bndry.inc'
!
!-----------------------------------------------------------------------
!
!  Universal physical constants such as gas constants.
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
!
!-----------------------------------------------------------------------
!
!  External boundary parameters and variables.
!
!-----------------------------------------------------------------------
!
  INCLUDE 'exbc.inc'
!
!-----------------------------------------------------------------------
!
!  Message passing parameters.
!
!-----------------------------------------------------------------------
!
  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
! Flags for fields readin from history format initial condition file
! e.g., tkein is used to determine whether tke exists in history IC file.
!-----------------------------------------------------------------------

  INCLUDE 'indtflg.inc'
!
!-----------------------------------------------------------------------
!
!  namelist Declarations:
!
!-----------------------------------------------------------------------
!
  NAMELIST /grid_dims/ nx,ny,nz

  NAMELIST /message_passing/ nproc_x,nproc_y,max_fopen,nproc_x_out,nproc_y_out
!
!-----------------------------------------------------------------------
!
!  Define the namelist, &arpsagr, for ARPS AGR. Not used by other programs
!  except for ARPSagr. Included here to avoid problem on e.g., Cray's
!  when reading arps.input.
!
!-----------------------------------------------------------------------
!
!  INCLUDE 'agricst.inc'
!  INCLUDE 'nodal.inc'
!  NAMELIST /arpsagr/  runold, rstime,                                  &
!                      levfix, intrat, intratt,                         &
!                      intrpodr, kcheck,                                &
!                      verbose1, verbose2, verbose3,                    &
!                      verbose4, verbose5, verbose6,                    &
!                      rstart,rstdump,grdsrt,                           &
!                      nfinelv, ngrdnew,                                &
!                      ixc,jyc,ixln,jyln,gangle

  NAMELIST /comment_lines/ nocmnt, cmnt

  NAMELIST /jobname/ runname

  NAMELIST /model_configuration/ runmod, hxopt, memid,hx_interval

  NAMELIST /initialization/ initime,initopt,inibasopt,viniopt,ubar0,    &
            vbar0,pt0opt,ptpert0,pt0radx,pt0rady,pt0radz,pt0ctrx,       &
            pt0ctry,pt0ctrz,rstinf,inifmt,inisplited,inifile,inigbf,    &
            sndfile,soilinitopt,soiltintv,timeopt,tsfcopt,              &
            ptground,pttrop,ttrop,htrop,qvmixed,rhmixed,mixtop,zshear

  NAMELIST /nudging/ nudgopt,ndstart,ndstop,ndintvl,ndgain,incrfnam,    &
                     nudgu,nudgv,nudgw,nudgp,nudgpt,nudgqv,             &
                     nudgqc,nudgqr,nudgqi,nudgqs,nudgqh,incrfmt,        &
                     dfilter_opt,df_tstart,df_tinv,df_nstps,df_wght

  NAMELIST /equation_formulation/ buoyopt,buoy2nd,rhofctopt,bsnesq,     &
            peqopt

  NAMELIST /terrain/ ternopt,mntopt,hmount,mntwidx,mntwidy,             &
            mntctrx,mntctry,terndta,ternfmt

  NAMELIST /grid/ dx,dy,dz,strhopt,dzmin,zrefsfc,dlayer1,dlayer2,       &
            strhtune,zflat,ctrlat,ctrlon, crdorgnopt

  NAMELIST /projection/ mapproj, trulat1,trulat2,trulon, sclfct,        &
            mpfctopt,mptrmopt,maptest

  NAMELIST /timestep/ dtbig,tstart,tstop

  NAMELIST /acoustic_wave/csopt,csfactr,csound,dtsml,vimplct,           &
            ptsmlstp,tacoef

  NAMELIST /numerics/ tintegopt, madvopt, sadvopt,fctorderopt,fctadvptprt

  NAMELIST /boundary_condition_options/ lbcopt, wbc,ebc,sbc,nbc,        &
            tbc,fftopt,bbc, rbcopt,rbc_plbc,c_phase,rlxlbc,pdetrnd

  NAMELIST /exbcpara/exbcname,tinitebd,tintvebd,                        &
            ngbrz,brlxhw,cbcdmp,exbcfmt

  NAMELIST /coriolis_force/ coriopt,earth_curvature,coriotrm

  NAMELIST /turbulence/ tmixopt,trbisotp,tkeopt,trbvimp,tmixvert,       &
            tmixcst, prantl, kmlimit

  NAMELIST /computational_mixing/                                       &
            cmix2nd,cfcm2h,cfcm2v,cmix4th,cfcm4h,cfcm4v,scmixfctr,      &
            cmix_opt

  NAMELIST /divergence_damping/ divdmp,divdmpndh, divdmpndv

  NAMELIST /concentration/ ccin,cpoint,icc,jcc,kcc,ccemit,ccstart,ccend !michi

  NAMELIST /rayleigh_damping/ raydmp,cfrdmp,zbrdmp

  NAMELIST /asselin_time_filter/ flteps

  NAMELIST /microphysics/ mphyopt,moist,nmphystp,cnvctopt,              &
            kffbfct,kfsubsattrig,wcldbs,confrq,qpfgfrq,idownd,          &
            subsatopt,rhsat,rhsatmin,dx_rhsatmin,dx_rhsat100,           &
            impfallopt,fallopt,                                         &
            dsdpref,ntcloud,n0rain,n0snow,n0grpl,n0hail,                &
            rhoice,rhosnow,rhogrpl,rhohail,alpharain,alphaice,          &
            alphasnow,alphagrpl,alphahail,graupel_ON,hail_ON,mpthermdiag

  NAMELIST /radiation/ radopt, radstgr, rlwopt, radshade, dtrad, raddiag

  NAMELIST /surface_physics/ sfcphy,landwtr,cdhwtropt,                  &
            cdmlnd,cdmwtr,cdhlnd,cdhwtr,cdqlnd,cdqwtr,                  &
            pbldopt,pbldpth0,lsclpbl0,tqflxdis,dtqflxdis,               &
            smthflx,numsmth,sfcdiag

  NAMELIST /soil_ebm/ sfcdat,soilmodel_forced,soilmodel_option,         &
            soilstrhopt,soilinit,dtsfc,styp,vtyp,lai0,roufns0,veg0,     &
            sitemeso,siteflux,siternet,sitesoil,siteveg,                &
            nzsoil,dzsoil,zrefsoil,tsoilint,qsoilint,                   &
            soildzmin,soildlayer1,soildlayer2,soilstrhtune,             &
            ptslnd0,ptswtr0,wetcanp0,snowdpth0,                         &
            ttprt,tbprt,wgrat,w2rat,                                    &
            sfcdtfl,soilinfl,sfcfmt,soilfmt,nstyp,                      &
            tsoil_offset, tsoil_offset_amplitude, prtsoilflx

  NAMELIST /grdtrans/cltkopt,grdtrns,umove,vmove,chkdpth,               &
            twindow,tceltrk,tcrestr

  NAMELIST /history_dump/ hdmpopt,dmp_out_joined,hdmpfmt,grbpkbit,      &
            thisdmp,tstrtdmp,numhdmp,hdmptim,istager,hdfcompr

  NAMELIST /output/ dirname,tfmtprt,exbcdmp,exbchdfcompr,extdadmp,      &
            grdout,basout,varout,mstout,rainout,prcout,                 &
            iceout,tkeout, trbout,sfcout,landout,totout,                &
            radout,flxout,                                              &
            qcexout,qrexout,qiexout,qsexout,qhexout,                    &
            qgexout,nqexout,zqexout,                                    &
            trstout,tmaxmin,tenergy,imgopt,                             &
            timgdmp,pltopt,tplots,filcmprs,readyfl,                     &
            sfcdmp,soildmp,terndmp

  NAMELIST /debug/ lvldbg

  REAL    :: dh
  INTEGER :: err_no
  DATA err_no /0/

  INTEGER :: ip
  INTEGER :: ebcsv,wbcsv,nbcsv,sbcsv

  INTEGER :: istatus

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

   istatus = 0
!-----------------------------------------------------------------------
!
!  Set the ARPS version number, which will be printed in the log
!  file in the comment line. The string can be up to 20 character long.
!
!-----------------------------------------------------------------------

  arpsversion = 'ARPS 5.3'
!
!-----------------------------------------------------------------------
!
!  Now we begin to read in the values of parameters:
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Set up the default values for all the variables to be read in
!  using the namelist method. In case the user does not specify a
!  particular value, this value will be used.
!
!-----------------------------------------------------------------------
!
  nx = 67
  ny = 67
  nz = 35

  nproc_x = 1
  nproc_y = 1
  max_fopen = 1

  nocmnt = 10

  cmnt(1) = ' '
  cmnt(2) = 'A zero perturbation run                    '
  cmnt(3) = ' '
  cmnt(4) = ' '
  cmnt(5) = ' '
  cmnt(6) = ' '
  cmnt(7) = ' '
  cmnt(8) = ' '
  cmnt(9) = ' '
  cmnt(10) =' '

  runname = 'may20'

  runmod = 1
  hxopt = 0
  memid = 1
  hx_interval = 60

  initime = '1977-05-20.21:00:00'
  timeopt = 0
  initopt = 1
  inibasopt  = 1
  inisplited = 0
  viniopt = 1
  pt0opt  = 1

  ubar0   = 0.0
  vbar0   = 0.0

  zshear = 3000.0
  pttrop = 343.0           ! Tropopause pot.temp.
  ttrop  = 213.0           ! Tropopause temp.
  ptground = 300.0         ! Groud urface pot.temp.
  htrop  = 12000.0         ! Tropopause height
  qvmixed= 0.015           ! Mixed layer mixing ratio
  rhmixed= 0.95
  mixtop = 1200.0          ! Mixed layer height

  ptpert0 = 0.0
  pt0radx = 0.0
  pt0rady = 0.0
  pt0radz = 0.0
  pt0ctrx = 0.0
  pt0ctry = 0.0
  pt0ctrz = 0.0

  ptpert0(1) = 4.0
  pt0radx(1) = 10000.0
  pt0rady(1) = 10000.0
  pt0radz(1) =  1500.0
  pt0ctrx(1) = 32000.0
  pt0ctry(1) = 32000.0
  pt0ctrz(1) =  1500.0

  rstinf  = 'may20.rst003600'
  inifmt  = 1
  inifile = 'may20.bin003600'
  inigbf  = 'may20.bingrdbas'
  sndfile = 'may20.snd'

  soilinitopt = 0
  soiltintv   = 0.0
  tsfcopt = 0

  nudgopt  = 0
  ndstart  = 0.
  ndstop   = 0.
  ndintvl  = 600.
  ndgain   = 1.9
  incrfnam = 'nudge.spam'
  incrfmt  = 1
  nudgu    = 1
  nudgv    = 1
  nudgw    = 1
  nudgp    = 1
  nudgpt   = 1
  nudgqv   = 1
  nudgqc   = 0
  nudgqr   = 0
  nudgqi   = 0
  nudgqs   = 0
  nudgqh   = 0

  dfilter_opt = 0
  df_tstart   = 0
  df_tinv  = 10
  df_nstps = 10
  df_wght  = 1.0/df_nstps

  ternopt = 0
  mntopt =  1
  hmount =  0.0
  mntwidx = 1.0E4
  mntwidy = 1.0E4
  mntctrx = 1.0E4
  mntctry = 1.0E4
  terndta ='arpstern.data'
  ternfmt = 1

  dx = 1000.0
  dy = 1000.0
  dz =  500.0

  strhopt  = 0
  dzmin    = 500.0
  zrefsfc  =   0.0
  dlayer1  =   0.0
  dlayer2  =   1.0E5
  strhtune =   1.0
  zflat    =   1.0E5

  ctrlat  =  35.0
  ctrlon  = -100.0
  crdorgnopt = 0

  mapproj  = 0
  trulat1  =   30.0
  trulat2  =   60.0
  trulon   = -100.0
  sclfct   =    1.0
  mpfctopt = 1
  mptrmopt = 1
  maptest  = 0

  dtbig = 6.0
  tstart= 0.0
  tstop = 3600.0

  vimplct  = 1
  ptsmlstp = 0
  csopt    = 1
  csfactr  =   0.5
  csound   = 150.0
  tacoef   =   0.6
  dtsml    =   1.0

  buoyopt   = 1
  buoy2nd   = 1
  rhofctopt = 1
  bsnesq    = 0
  peqopt    = 1

  tintegopt = 1
  madvopt  = 1
  sadvopt  = 1
  fctorderopt=1
  fctadvptprt=1

  lbcopt = 1
  wbc = 4
  ebc = 4
  sbc = 4
  nbc = 4
  tbc = 1
  fftopt = 2
  bbc = 1
  rbcopt  = 1
  rbc_plbc = 1
  c_phase = 30.0
  rlxlbc  =  0.0
  pdetrnd = 0

  radopt  = 0
  radstgr = 1
  rlwopt  = 1
  radshade = 0
  dtrad   = 600.0
  raddiag = 1

  moist    = 0
  mphyopt  = 0
  nscalar  = 0
  nscalarq = 0
  P_QC = -1; P_QR = -1; P_QI = -1; P_QS = -1; P_QH = -1; P_QG = -1
  P_NC = -1; P_NR = -1; P_NI = -1; P_NS = -1; P_NG = -1; P_NH = -1
             P_ZR = -1; P_ZI = -1; P_ZS = -1; P_ZG = -1; P_ZH = -1;
  P_CC = -1
  qnames(:)= '    '; qdescp(:)= ' '
  nmphystp = 1
  cnvctopt = 0
  subsatopt = 0
  kffbfct  = 0.0
  kfsubsattrig = 0
  ice      = 0
  wcldbs   = 0.005
  confrq   = 600.0
  qpfgfrq  = 120.0
  idownd   = 1
  dsdpref  = 0
  ntcloud  = 1.0e8 ! Default for MY scheme (not used for LIN scheme)
  n0rain   = 8.0e6 ! Default for LIN scheme
  n0snow   = 3.0e6 ! Default for LIN scheme
  n0grpl   = 4.0e5 ! Default for MY scheme (not used for LIN scheme)
  n0hail   = 4.0e4 ! Default for LIN scheme
  rhoice   = 500.0 ! Default for MY scheme (not used for LIN scheme)
  rhosnow  = 100.0 ! Default for LIN and MY scheme
  rhogrpl  = 400.0 ! Default for MY scheme (not used for LIN scheme)
  rhohail  = 913.0 ! Default for LIN scheme

  alpharain = 0.0  ! Default for MY scheme (not used for LIN scheme)
  alphaice = 0.0   ! Default for MY scheme (not used for LIN scheme)
  alphasnow = 0.0  ! Default for MY scheme (not used for LIN scheme)
  alphagrpl = 0.0  ! Default for MY scheme (not used for LIN scheme)
  alphahail = 0.0  ! Default for MY scheme (not used for LIN scheme)

  graupel_ON = 1
  hail_ON = 1      !BJP OCT 2010 hail switch

  mpthermdiag = 0

  impfallopt = 0
  fallopt = 1

  rhsat    = 0.80
  rhsatmin = 0.80
  dx_rhsatmin = 50000.
  dx_rhsat100 = 5000.

  sfcphy   = 0
  landwtr  = 1
  cdhwtropt= 0
  cdmlnd   = 3.0E-3
  cdmwtr   = 1.0E-3
  cdhlnd   = 3.0E-3
  cdhwtr   = 1.0E-3
  cdqlnd   = 2.1E-3
  cdqwtr   = 0.7E-3
  pbldopt  = 0
  pbldpth0 = 1400.0
  lsclpbl0 = 0.15
  sflxdis  = 0
  tqflxdis = 0
  dtqflxdis= 200.0
  smthflx  = 0
  numsmth  = 1
  sfcdiag  = 0

  sfcdat  = 1
  nstyp  = 4
  styp    = 3
  vtyp    = 10
  lai0    = 0.31
  roufns0 = 0.1
  veg0    = 0.0
  sfcdtfl = 'arpssfc.data'
  sfcfmt  = 1

  soilmodel_forced = 0
  sitemeso = '../../arpsdata.dir/mtsnorm.dir/Meso'
  siteflux = '../../arpsdata.dir/mtsnorm.dir/Flux'
  siternet = '../../arpsdata.dir/mtsnorm.dir/Radd'
  sitesoil = '../../arpsdata.dir/mtsnorm.dir/Soil'
  siteveg =  '../../arpsdata.dir/mtsnorm.dir/Veg'

  soilmodel_option = 1
  nzsoil     = 2
  dzsoil     = 1.00
  zrefsoil   = 0.0
  DO i=1,nzsoil
    tsoilint(i) = 273.15
    qsoilint(i) = 0.50
  END DO

  soilstrhopt = 0
  soildzmin   = 0.025
  soildlayer1 = 0.0
  soildlayer2 = 1.0
  soilstrhtune = 1.0

  soilinit = 1
  ptslnd0  = 300.16
  ptswtr0  = 288.16
  wetcanp0 = 0.0
  snowdpth0 = 0
  soilinfl = 'may20.soilinit'
  soilfmt  = 1

  dtsfc = 10.0

  prtsoilflx = 0

  coriopt = 0
  earth_curvature = 0
  coriotrm= 0

  tmixopt  = 2
  trbisotp = 1
  tkeopt   = 1
  trbvimp  = 0
  tmixvert = 1
  prantl   = 1.0
  tmixcst  = 0.0
  kmlimit  = 0.5

  cmix2nd = 1
  cfcm2h  = 0.0
  cfcm2v  = 1.0E-3
  cmix4th = 1
  cfcm4h  = 1.0E-3
  cfcm4v  = 0.0
  scmixfctr = 1.0
  cmix_opt = 0

  divdmp    = 1
  divdmpndh = 0.05
  divdmpndv = 0.05

  tmaxmin  = 60.0
  tenergy  = 360000.0
  imgopt   = 0
  timgdmp  = 60.0
  pltopt   = 0
  tplots   = 1800.0
  filcmprs = 1
  readyfl = 0

  raydmp = 0
  cfrdmp = 1./300.
  zbrdmp = 10000.0

  flteps = 0.10

  cltkopt = 0
  grdtrns = 0
  umove = 0.0
  vmove = 0.0
  chkdpth = 2500.0
  twindow = 33333
  tceltrk  = 120.0
  tcrestr  = 1800.0

  lvldbg = 0

  hdmpopt  = 1
  dmp_out_joined = 0
  hdmpfmt  = 3
  grbpkbit = 16
  hdfcompr = 0
  thisdmp  = 3600.0
  tstrtdmp = 0.0
  numhdmp  = 1
  DO i=1,numhdmp
    hdmptim(i) = 0.
  END DO
  istager = 0

  dirname  = ' '
  tfmtprt  = 3600.0
  exbcdmp  = 0
  exbchdfcompr = 0
  extdadmp = 0
  grdout   = 0
  basout   = 0
  varout   = 1
  mstout   = 1
  rainout  = 0
  prcout   = 0
  iceout   = 0
  totout   = 1
  tkeout   = 0
  trbout   = 0
  sfcout   = 0
  snowout  = 0
  landout  = 0
  radout   = 0
  flxout   = 0

  grdin = 0
  basin = 0
  varin = 0
  mstin = 0
  rainin= 0
  prcin = 0
  icein = 0
  tkein = 0
  trbin = 0
  sfcin = 0
  landin= 0
  totin = 0
  radin = 0
  flxin = 0
  snowcin=0
  snowin= 0

  qcexout = 0
  qrexout = 0
  qiexout = 0
  qsexout = 0
  qhexout = 0
  qgexout = 0
  nqexout = 0
  zqexout = 0

  sfcdmp   = 0
  soildmp  = 0
  terndmp  = 0

  trstout  = 3600.0

  exbcname = 'arpsexbc'
  tinitebd = '1977-05-20.15:00:00'
  tintvebd = 10800
  ngbrz    = 5
  brlxhw   = 2.3
  cbcdmp   = 0.0033333333
  cbcmix   = 0.0                 ! dropped from namelist since ARPS5.2.4
  exbcfmt  = 1

  mgrid = 1
!
!-----------------------------------------------------------------------
!
!  Initialize message passing processors.
!
!-----------------------------------------------------------------------
!
  ! Non-MPI defaults:
  mp_opt = 0
  myproc = 0
  loc_x  = 1
  loc_y  = 1

  nproc_x_in  = 1
  nproc_y_in  = 1
  nproc_x_out = 1
  nproc_y_out = 1

  readsplit(:) = 0
  joindmp(:)   = 0
!
!-----------------------------------------------------------------------
!
!      Initialize the processors for an MPI job.
!
!-----------------------------------------------------------------------
!

  IF (nstyps == -884889488) THEN    ! 88d2arps (C code) specific
    CALL mpinit_proc(1)
  ELSE
    CALL mpinit_proc(0)
  END IF

  IF (myproc == 0)THEN

    WRITE(6,'(/ 16(/5x,a)//)')                                            &
       '###############################################################', &
       '###############################################################', &
       '#####                                                     #####', &
       '#####                      Welcome to                     #####', &
       '#####                                                     #####', &
       '#####   The Advanced Regional Prediction System  (ARPS)   #####', &
       '#####                                                     #####', &
       '#####                     Version 5.3                     #####', &
       '#####                                                     #####', &
       '#####                     Developed by                    #####', &
       '#####     Center for Analysis and Prediction of Storms    #####', &
       '#####                University of Oklahoma               #####', &
       '#####                                                     #####', &
       '###############################################################', &
       '###############################################################'

    WRITE(6,'(5(/5x,a),/6(/5x,a)/)')                                      &
        'The model begins by reading a number of control parameters,',    &
        'which are specified in namelist format through the standard',    &
        'input stream (unit 5).  See the ARPS Users Guide and the',       &
        'sample input file, arps.input, for guidance on specifying',      &
        'these parameters',                                               &
        'At the end of all parameter input, a log file is produced',      &
        'which can be directly used as the input file when you want',     &
        'to replicate the same job. This file is named runnam.log.nn',    &
        'where runnam is a standard prefix for all output files that',    &
        'are produced by this job and nn is a number appended to the',    &
        'file name when file runnam.log.nn-1 already exists.'

    unum = 0
    IF (LEN_TRIM(namelist_filename) <= 0 .OR. namelist_filename == ' ') THEN
      !CALL getarg(1,nlfile)
      CALL GET_COMMAND_ARGUMENT(1, nlfile, lenstr, istatus )
      IF ( nlfile(1:1) == ' ' .OR. istatus /= 0 ) THEN  ! Use standard input to be backward-compatible
        unum = 5
      ELSE
        INQUIRE(FILE=nlfile,EXIST=iexist)
        IF (.NOT. iexist) THEN
          WRITE(6,'(1x,3a)') 'WARNING: namelist file - ',               &
                TRIM(nlfile),' does not exist. Falling back to standard input.'
          unum = 5
        END IF
      END IF
    ELSE
      nlfile = namelist_filename
    END IF

    IF (unum /= 5) THEN
      namelist_filename = nlfile
      CALL getunit( unum )
      OPEN(unum,FILE=TRIM(nlfile),STATUS='OLD',FORM='FORMATTED')
      WRITE(6,'(1x,3a,/,1x,a,/)') 'Reading ARPS namelist from file - ',   &
              TRIM(nlfile),' ... ','========================================'
    ELSE
      WRITE(6,'(2(1x,a,/))') 'Waiting namelist from standard input ... ', &
                             '========================================'
    END IF

  END IF

!
!-----------------------------------------------------------------------
!
!  Read in grid dimensions.
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) THEN
    READ (unum,grid_dims, END=100)
    WRITE(6,'(a)')'Namelist block grid_dims sucessfully read.'

    WRITE(6,'(5x,a,i5)') "nx =",nx
    WRITE(6,'(5x,a,i5)') "ny =",ny
    WRITE(6,'(5x,a,i5)') "nz =",nz
  END IF

  CALL mpupdatei(nx,1)
  CALL mpupdatei(ny,1)
  CALL mpupdatei(nz,1)

!
!-----------------------------------------------------------------------
!
!  Read in message passing options.
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) READ (unum,message_passing, END=100)

  IF (myproc == 0) &
    WRITE(6,'(a)')'Namelist block message_passing sucessfully read.'

  CALL mpupdatei(nproc_x,1)
  CALL mpupdatei(nproc_y,1)
  CALL mpupdatei(max_fopen,1)
  CALL mpupdatei(nproc_x_out,1)
  CALL mpupdatei(nproc_y_out,1)

  nproc_x_in = nproc_x        ! Should be run-time parameters, it is hard-coded here
  nproc_y_in = nproc_y        ! for backward-compatibility

  IF (myproc == 0)THEN
    WRITE(6,'(5x,a,i4)')                                                &
         "Number of processors in the x-direction is:",nproc_x
    WRITE(6,'(5x,a,i4)')                                                &
        "Number of processors in the y-direction is:",nproc_y
    WRITE(6,'(5x,a,i4)')                                                &
        "Maximum number of files open:",max_fopen
  END IF

  ! Note that for MP version nx & ny here are global values.  They will
  ! be reassigned to their per-processor value below.

  IF (mp_opt > 0) THEN

    IF (nx /= nproc_x*int((nx-3)/nproc_x)+3) THEN
      nx = nproc_x*int((nx-3)/nproc_x+0.9999999999999) + 3
      IF (myproc == 0) THEN
        WRITE (6,*) "WARNING: adjusting nx to fit on ",nproc_x," processors:"
        WRITE(6,'(5x,a,i5)') "   new nx =",nx
      ENDIF
    ENDIF
    IF (ny /= nproc_y*int((ny-3)/nproc_y)+3) THEN
      ny = nproc_y*int((ny-3)/nproc_y+0.9999999999999) + 3
      IF (myproc == 0) THEN
        WRITE (6,*) "WARNING: adjusting ny to fit on ",nproc_y," processors:"
        WRITE(6,'(5x,a,i5)') "   new ny =",ny
      ENDIF
    ENDIF

  ELSE

    nproc_x = 1
    nproc_y = 1
    nprocs = 1
    max_fopen = 1

  ENDIF

!
!-----------------------------------------------------------------------
!
!  Initialize message passing variables.
!
!-----------------------------------------------------------------------
!
  CALL mpinit_var

!-----------------------------------------------------------------------
!  Read in agri name list block - not used by single grid arps run
!-----------------------------------------------------------------------

!  IF (myproc == 0) READ (unum,arpsagr, END=100)
!  IF (myproc == 0) &
!    WRITE(6,'(a)')'Namelist block arpsagr sucessfully read.'
!
!-----------------------------------------------------------------------
!
!  Read in some comment lines on this job and the name of
!  this run designated by a string at least 6 character long.
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) READ (unum,comment_lines, END=100)
  IF (myproc == 0) &
    WRITE(6,'(a)')'Namelist block comment_lines sucessfully read.'
  CALL mpupdatei(nocmnt,1)
  CALL mpupdatec(cmnt,80*nocmnt)

  IF (myproc == 0) READ (unum,jobname,END=100)
  IF (myproc == 0) &
    WRITE(6,'(a)')'Namelist block jobname sucessfully read.'
  CALL mpupdatec(runname,80)

  IF (myproc == 0) &
    WRITE(6,'(/5x,a,a)') 'The name of this run is: ', runname
!
!-----------------------------------------------------------------------
!
!  Find out the number of characters to be used to construct file
!  names.
!
!-----------------------------------------------------------------------
!
  CALL gtlfnkey( runname, lfnkey )
!
!-----------------------------------------------------------------------
!
!  Read in the parameter that controls the model run mode:
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0)THEN
    READ (unum,model_configuration,END=100)
    WRITE(6,'(a)')                                                      &
        'Namelist block model_configuration sucessfully read.'
  END IF
  CALL mpupdatei(runmod,1)
  CALL mpupdatei(hxopt,1)
  CALL mpupdatei(memid,1)
  CALL mpupdatei(hx_interval,1)

  IF (myproc == 0) WRITE (6,'(/5x,a,i4)') 'The run mode is: ', runmod
!
!-----------------------------------------------------------------------
!
!  Read in control parameter INITOPT for model initialization
!
!  INITOPT = 1, Self initialization (e.g. specify perturbation using
!                 analytical functions),
!            = 2, Restart run, initialize the model using previous
!                 model output,
!            = 3, Initialize the model using external input data file.
!            = 4, Read in restart file then overwrite current time level
!                 with those in history file.
!
!  For options 2-4, the names of input files need to be provided.
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) THEN
    READ (unum,initialization,END=100)
    WRITE(6,'(a)')'Namelist block initialization sucessfully read.'
  END IF

  DO i = FINDX_NUM, 1, -1
    idummy = inisplited / 10**(i-1)     ! ith digit is odd (1,3,5..), no readsplit
    readsplit(i) = MOD(idummy+1, 2)     ! ith digit is even(0,2,4..), readsplit
  END DO

  CALL mpupdatec(initime,19)
  CALL mpupdatei(initopt,1)
  CALL mpupdatei(inibasopt,1)
  CALL mpupdatei(viniopt,1)
  CALL mpupdater(ubar0,1)
  CALL mpupdater(vbar0,1)

  CALL mpupdater(zshear,1)
  CALL mpupdater(pttrop,1)
  CALL mpupdater(ttrop,1)
  CALL mpupdater(ptground,1)
  CALL mpupdater(htrop,1)
  CALL mpupdater(qvmixed,1)
  CALL mpupdater(rhmixed,1)
  CALL mpupdater(mixtop,1)

  CALL mpupdatei(pt0opt,1)
  CALL mpupdater(ptpert0,maxpert)
  CALL mpupdater(pt0radx,maxpert)
  CALL mpupdater(pt0rady,maxpert)
  CALL mpupdater(pt0radz,maxpert)
  CALL mpupdater(pt0ctrx,maxpert)
  CALL mpupdater(pt0ctry,maxpert)
  CALL mpupdater(pt0ctrz,maxpert)
  CALL mpupdatec(rstinf,256)
  CALL mpupdatei(inifmt,1)
  CALL mpupdatei(readsplit,FINDX_NUM)
  CALL mpupdatec(inifile,256)
  CALL mpupdatec(inigbf, 256)
  CALL mpupdatec(sndfile,256)
  CALL mpupdatei(soilinitopt,1)
  CALL mpupdater(tsfcopt,1)
  CALL mpupdater(soiltintv,1)
  CALL mpupdatei(timeopt,1)

  READ (initime, '(i4.4,1x,i2.2,1x,i2.2,1x,i2.2,1x,i2.2,1x,i2.2)' )     &
       year,month,day,hour,minute,second

  IF (myproc == 0)THEN
    WRITE(6,'(a,a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)')            &
        '     The initial local time for this run is ',                 &
        '     year-mo-dy:hr:mn:ss = ',                                  &
        year,'-',month,'-',day,'.',hour,':',minute,':',second

    WRITE(6,'(5x,a,i4/)') 'Perturbation option was ', pt0opt

    ip = 1

    DO WHILE (ptpert0(ip) /= 0.0)
      WRITE(6,'(5x,a,i2,a,f10.3,a)')                                           &
          'The magnitude of the initial perturbation ',ip,' is ', ptpert0(ip),' K.'

      WRITE(6,'(5x,a,3e10.3,a,/5x,a,a,/5x,a,e10.3,a,e10.3,a,e10.3,a)')    &
          'The input radii of the thermal bubble are ',                   &
          pt0radx(ip), pt0rady(ip), pt0radz(ip),' (m)',                               &
          'in x, y and z direction respectively, and the center is ',     &
          'located at','x=',pt0ctrx(ip),' y=',pt0ctry(ip),' z=',pt0ctrz(ip),          &
          ' (m).'

      ip = ip + 1
    END DO
    WRITE(6,'(5x,a,i4/)') 'tsfc option was ', tsfcopt
  END IF
!
!-----------------------------------------------------------------------
!
!  Input data files for initialization:
!
!-----------------------------------------------------------------------
!
  lenstr = 256
  CALL strlnth( rstinf, lenstr)
  IF (myproc == 0)THEN
    WRITE(6,'(5x,a,a)')                                                 &
        'The two time level restart data to be read in is ',            &
        rstinf(1:lenstr)

    WRITE(6,'(5x,a,i4)')                                                &
        'The history dump type restart data format was ',inifmt
  END IF

  lenstr = 256
  CALL strlnth( inifile, lenstr)
  IF (myproc == 0)THEN
    WRITE(6,'(5x,a,a)')                                                 &
      'The t-dependent history dump format restart data to be read is ',&
      inifile(1:lenstr)
  END IF

  lenstr = 256
  CALL strlnth( inigbf, lenstr)
  IF (myproc == 0)THEN
    WRITE(6,'(5x,a,a,a)')                                               &
        'The base state/grid history dump ',                            &
        'format restart data to be read is ', inigbf(1:lenstr)
  END IF
!
!-----------------------------------------------------------------------
!
!  Input the environmental sounding.
!
!-----------------------------------------------------------------------
!
  lenstr = 256
  CALL strlnth( sndfile, lenstr )
  IF (myproc == 0)THEN
    WRITE(6,'(5x,a,a)')                                                 &
         'Sounding file to be used is ', sndfile(1:lenstr)
  END IF
!
!-----------------------------------------------------------------------
!
!  Read the nudging options.
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) THEN
    READ (unum,nudging,END=100)
    WRITE(6,'(a)') 'Namelist block nudging sucessfully read.'
  END IF
  CALL mpupdatei(nudgopt,1)
  CALL mpupdater(ndstart,1)
  CALL mpupdater(ndstop,1)
  CALL mpupdater(ndintvl,1)
  CALL mpupdater(ndgain,1)
  CALL mpupdatec(incrfnam,256)
  CALL mpupdatei(nudgu,1)
  CALL mpupdatei(nudgv,1)
  CALL mpupdatei(nudgw,1)
  CALL mpupdatei(nudgp,1)
  CALL mpupdatei(nudgpt,1)
  CALL mpupdatei(nudgqv,1)
  CALL mpupdatei(nudgqc,1)
  CALL mpupdatei(nudgqr,1)
  CALL mpupdatei(nudgqi,1)
  CALL mpupdatei(nudgqs,1)
  CALL mpupdatei(nudgqh,1)
  CALL mpupdatei(incrfmt,1)

  CALL mpupdatei(dfilter_opt,1)
  CALL mpupdater(df_tstart,1)
  CALL mpupdater(df_tinv,1)
  CALL mpupdatei(df_nstps,1)
  CALL mpupdater(df_wght,df_maxnstps)

  IF (myproc == 0)THEN
    WRITE(6,'(5x,a,i4)')                                                &
          'The nudging assimilation option was ',nudgopt

    WRITE(6,'(5x,a,f9.2)') 'Nudging assimilation start: ',ndstart

    WRITE(6,'(5x,a,f9.2)') 'Nudging assimilation stop: ',ndstop

    WRITE(6,'(5x,a,f9.2)') 'Nudging assimilation interval: ',ndintvl

    WRITE(6,'(5x,a,f9.2)') 'Nudging assimilation gain: ',ndgain

  END IF

  IF (myproc == 0) WRITE(6,'(5x,a,a)')                                  &
        'The nudging increment file is ',TRIM(incrfnam)
!
!-----------------------------------------------------------------------
!
!  Specify the types of terrain option:
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0)THEN
    READ (unum,terrain,END=100)
    WRITE(6,'(a)')'Namelist block terrain sucessfully read.'
  END IF
  CALL mpupdatei(ternopt,1)
  CALL mpupdatei(mntopt, 1)
  CALL mpupdater(hmount, 1)
  CALL mpupdater(mntwidx,1)
  CALL mpupdater(mntwidy,1)
  CALL mpupdater(mntctrx,1)
  CALL mpupdater(mntctry,1)
  CALL mpupdatec(terndta,256)
  CALL mpupdatei(ternfmt,1)

  IF (myproc == 0)THEN

    WRITE(6,'(5x,a,i4)') 'The mountain type option was ', mntopt

    WRITE(6,'(5x,a,f10.3,a)')                                           &
        'The height of mountain is ', hmount,' (m).'

    WRITE(6,'(5x,a,2e10.3,a,/5x,a,a,/5x,a,e10.3,a,e10.3,a)')            &
        'The input half-width of bell-shaped mountain are ',            &
        mntwidx, mntwidy, ' (m)',                                       &
        'in x and y direction respectively, and the center is ',        &
        'located at','x=',mntctrx,' y=',mntctry,' (m).'
    WRITE(6,'(5x,a,i4)') 'The terrain option was ', ternopt
  END IF

  lenstr = 256
  CALL strlnth( terndta, lenstr)

  IF (myproc == 0)THEN
    WRITE(6,'(5x,a,a)')                                              &
      'The terrain data file is ',terndta(1:lenstr)
    WRITE(6,'(5x,a,i4)') 'The terrain data file format is ', ternfmt
  END IF
!
!-----------------------------------------------------------------------
!
!  Input horizontal grid size
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) THEN
    READ (unum,grid,END=100)
    WRITE(6,'(a)')'Namelist block grid sucessfully read.'
  END IF

  CALL mpupdater(dx,1)
  CALL mpupdater(dy,1)
  CALL mpupdater(dz,1)
  CALL mpupdatei(strhopt,1)
  CALL mpupdater(dzmin,1)
  CALL mpupdater(zrefsfc,1)
  CALL mpupdater(dlayer1,1)
  CALL mpupdater(dlayer2,1)
  CALL mpupdater(strhtune,1)
  CALL mpupdater(zflat,1)
  CALL mpupdater(ctrlat,1)
  CALL mpupdater(ctrlon,1)
  CALL mpupdatei(crdorgnopt,1)

  IF( strhopt == 0.AND.dzmin /= dz ) THEN
    dzmin = dz
    IF (myproc == 0)  WRITE(6,'(5x,a)')                               &
         'For non-stretched case, dzmin was reset to dz.'
  END IF

  IF (myproc == 0)THEN
    WRITE(6,'(5x,a,f10.3,a)')                                         &
       'Input dx was',dx,' meters'

    WRITE(6,'(5x,a,f10.3,a)')                                         &
       'Input dy was ',dy,' meters'

    WRITE(6,'(5x,a,i4)') 'The stretch option was ', strhopt

    WRITE(6,'(5x,a,f10.3,a)')                                         &
       'Input dz was ',dz,' meters'

    WRITE(6,'(5x,a,f10.3,a)')                                         &
       'Input ctrlat was ',ctrlat,' degrees North'

    WRITE(6,'(5x,a,f10.3,a)')                                         &
       'Input ctrlon was ',ctrlon,' degrees East'

    WRITE(6,'(5x,a,i4)') 'domain origin option was', crdorgnopt

    WRITE(6,'(5x,a,f10.3,a)')                                         &
       'dzmin is ',dzmin,' meters'

    WRITE(6,'(5x,a,f10.3,a)')                                         &
       'zrefsfc is ',zrefsfc ,' meters'

    WRITE(6,'(5x,a,f10.3,a)')                                         &
       'dlayer1 is ',dlayer1 ,' meters'

    WRITE(6,'(5x,a,f10.3,a)')                                         &
       'dlayer2 is ',dlayer2 ,' meters'

    WRITE(6,'(5x,a,f10.3,a)')                                         &
       'zflat   is ',zflat   ,' meters'
  END IF
!
!-----------------------------------------------------------------------
!
!  Input map projection parameters
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) THEN
    READ (unum,projection,END=100)
    WRITE(6,'(a)')'Namelist block projection sucessfully read.'
  END IF

  CALL mpupdatei(mapproj,1)
  CALL mpupdater(trulat1,1)
  CALL mpupdater(trulat2,1)
  CALL mpupdater(trulon,1)
  CALL mpupdater(sclfct,1)
  CALL mpupdatei(mpfctopt,1)
  CALL mpupdatei(mptrmopt,1)
  CALL mpupdatei(maptest,1)

  IF (myproc == 0)THEN
    WRITE(6,'(5x,a,i4)')                                            &
         'Input mapproj was ',mapproj

    WRITE(6,'(5x,a,f10.3,a)')                                       &
         'Input trulat1 was ',trulat1,' degrees North'

    WRITE(6,'(5x,a,f10.3,a)')                                       &
         'Input trulat2 was ',trulat2,' degrees North'

    WRITE(6,'(5x,a,f10.3)')                                         &
        'The latitude of the center of the model domain was ',ctrlat

    WRITE(6,'(5x,a,f10.3)')                                         &
        'The longitude of the center of the model domain was ',ctrlon

    WRITE(6,'(5x,a,f10.3,a)')                                       &
         'Input trulon was ',trulon,' degrees East'

    WRITE(6,'(5x,a,e15.5)')                                         &
         'Input sclfct was ',sclfct

    WRITE(6,'(5x,a,i5)')                                            &
         'The option for map factor was ', mpfctopt

    WRITE(6,'(5x,a,i5)')                                            &
       'The option for map factor term in u and v advection was ',  &
       mptrmopt
  END IF

  IF (myproc == 0) THEN
    READ (unum,timestep,END=100)
    WRITE(6,'(a)')'Namelist block timestep sucessfully read.'
  END IF

  CALL mpupdater(dtbig,1)
  CALL mpupdater(tstart,1)
  CALL mpupdater(tstop,1)

  IF (myproc == 0)THEN
  WRITE(6,'(5x,a,f10.3,a)')                                             &
       'The big timestep was ',dtbig,' seconds.'

  WRITE(6,'(5x,a,f10.3,a)')                                             &
       'The model startup time was ',tstart, ' seconds.'

  WRITE(6,'(5x,a,f10.3,a)')                                             &
       'The termination time was ',tstop, ' seconds.'
  END IF

  IF (myproc == 0) READ (unum,acoustic_wave,END=100)
  IF (myproc == 0)  WRITE(6,'(a)')'Namelist block acoustic_wave sucessfully read.'
  CALL mpupdatei(csopt,1)
  CALL mpupdater(csfactr,1)
  CALL mpupdater(csound,1)
  CALL mpupdater(dtsml,1)
  CALL mpupdatei(vimplct,1)
  CALL mpupdatei(ptsmlstp,1)
  CALL mpupdater(tacoef,1)

  IF (myproc == 0)THEN
  WRITE(6,'(5x,a,i4)') 'The sound speed option was ',csopt

  WRITE(6,'(5x,a,f10.3)')                                               &
       'The reduction factor for sound speed was ', csfactr

  WRITE(6,'(5x,a,f10.3,a)')                                             &
       'The constant sound speed was ', csound,' m/s.'

  WRITE(6,'(5x,a,a,i5)')                                                &
       'The vertical implicit integration option for ',                 &
       'w and p equations was ', vimplct

  WRITE(6,'(5x,a,a,i5)')                                                &
       'The option for potential temperature equation integration',     &
       'was ', ptsmlstp

  WRITE(6,'(5x,a,a,f10.3)')                                             &
       'The time average coefficient for vertically ',                  &
       'implicit solver was ', tacoef

  WRITE(6,'(5x,a,f10.3,a)')                                             &
       'The input small timestep was ',dtsml,' seconds.'

  WRITE(6,'(5x,a,f10.3,a)')                                             &
       'The actual small time step size to be used is ',                &
       dtsml,' seconds.'
  END IF
!
!-----------------------------------------------------------------------
!
!  Read in parameters related to equation formaulation
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) THEN
    READ (unum,equation_formulation,END=100)
    WRITE(6,'(a)')                                                 &
      'Namelist block equation_formulation sucessfully read.'
  END IF

  CALL mpupdatei(buoyopt,1)
  CALL mpupdatei(buoy2nd,1)
  CALL mpupdatei(rhofctopt,1)
  CALL mpupdatei(bsnesq,1)
  CALL mpupdatei(peqopt,1)

  IF ( buoyopt == 0 ) THEN
    IF(myproc == 0) WRITE(6,*) 'WARNING: buoyancy terms turned off',  &
               ' by selecting buoyopt=0.'
  END IF

  IF (myproc == 0) WRITE(6,'(5x,a,i5)')                               &
       'The option for pressure equation formulation was ',peqopt
!
!-----------------------------------------------------------------------
!
!  Read in parameters related to numerical schemes
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) THEN
    READ (unum,numerics,END=100)
    WRITE(6,'(a)')'Namelist block numerics sucessfully read.'
  END IF
  CALL mpupdatei(tintegopt,1)
  CALL mpupdatei(madvopt,1)
  CALL mpupdatei(sadvopt,1)
  CALL mpupdatei(fctorderopt,1)
  CALL mpupdatei(fctadvptprt,1)

  IF (myproc == 0)THEN
  WRITE(6,'(5x,a,i5)')                                                  &
       'The option for time integration was ', tintegopt

  IF( tintegopt < 1 .OR. tintegopt > 3 ) THEN
    WRITE(6,'(5x,a,i3,a,2(/5x,a))')                                     &
        'Input value of tintegopt= ', tintegopt,' was invalid.',            &
        'Program will try to complete reading in input parameters, ',   &
        'but will stop at the end of subroutine INITPARA.'
    err_no = err_no + 1
  END IF

  WRITE(6,'(5x,a,i5)')                                                  &
       'The option for momentum advection was ', madvopt

  IF( madvopt < 1 .OR. madvopt > 3 ) THEN
    WRITE(6,'(5x,a,i3,a,2(/5x,a))')                                     &
        'Input value of madvopt= ', madvopt,' was invalid.',            &
        'Program will try to complete reading in input parameters, ',   &
        'but will stop at the end of subroutine INITPARA.'
    err_no = err_no + 1
  END IF

  WRITE(6,'(5x,a,i5)')                                                  &
       'The option for scalar   advection was ', sadvopt

  IF( sadvopt < 1 .OR. sadvopt > 5 ) THEN
    WRITE(6,'(5x,a,i3,a,2(/5x,a))')                                     &
        'Input value of sadvopt= ', sadvopt,' was invalid.',            &
        'Program will try to complete reading in input parameters, ',   &
        'but will stop at the end of subroutine INITPARA.'
    err_no = err_no + 1
  END IF

  IF( fctorderopt /= 1 .AND. fctorderopt /= 2 ) THEN
    WRITE(6,'(5x,a,i3,a,2(/5x,a))')                                     &
        'Input value of fctorderopt= ', fctorderopt,' was invalid.',    &
        'Program will try to complete reading in input parameters, ',   &
        'but will stop at the end of subroutine INITPARA.'
    err_no = err_no + 1
  END IF

  IF( fctadvptprt < 0 .OR. fctadvptprt > 2 ) THEN
    WRITE(6,'(5x,a,i3,a,2(/5x,a))')                                     &
        'Input value of fctadvptprt= ', fctadvptprt,' was invalid.',    &
        'Program will try to complete reading in input parameters, ',   &
        'but will stop at the end of subroutine INITPARA.'
    err_no = err_no + 1
  END IF

  IF(sadvopt == 4.AND.ptsmlstp == 1.AND.fctadvptprt /= 1) THEN
    WRITE(6,'(2(/5x,a))')                                               &
        'When sadvopt=4, and ptsmlstp=1, fctadvptprt has to be 1',      &
        'fctadvptprt reset to 1'

    fctadvptprt = 1

  END IF
  END IF
!
!-----------------------------------------------------------------------
!
!  Input boundary condition control parameters:
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) READ (unum,boundary_condition_options,END=100)
  IF (myproc == 0)THEN
  WRITE(6,'(a)')                                                        &
      'Namelist block boundary_condition_options sucessfully read.'
  END IF
  CALL mpupdatei(lbcopt,1)
  CALL mpupdatei(wbc,1)
  CALL mpupdatei(ebc,1)
  CALL mpupdatei(sbc,1)
  CALL mpupdatei(nbc,1)
  CALL mpupdatei(tbc,1)
  CALL mpupdatei(fftopt,1)
  CALL mpupdatei(bbc,1)
  CALL mpupdatei(rbcopt,1)
  CALL mpupdatei(rbc_plbc,1)
  CALL mpupdater(c_phase,1)
  CALL mpupdater(rlxlbc,1)
  CALL mpupdatei(pdetrnd,1)
!
!-----------------------------------------------------------------------
!
!  For 2-D or 1-D runs, appropriate boundary conditions are
!  automatically set to periodic.
!
!-----------------------------------------------------------------------
!
  IF( runmod == 2 .OR. runmod == 4 ) THEN
    IF( nbc /= 2 ) THEN
  IF (myproc == 0)THEN
      WRITE(6,'(5x,a,i2)') 'nbc reset to 2 for runmod=',runmod
  END IF
      nbc = 2
    END IF

    IF( sbc /= 2 ) THEN
  IF (myproc == 0)THEN
      WRITE(6,'(5x,a,i2)') 'sbc reset to 2 for runmod=',runmod
  END IF
      sbc = 2
    END IF
  END IF

  IF( runmod == 3 .OR. runmod == 4 ) THEN
    IF( wbc /= 2 ) THEN
  IF (myproc == 0)THEN
      WRITE(6,'(5x,a,i2)') 'wbc reset to 2 for runmod=',runmod
  END IF
      wbc = 2
    END IF
    IF( ebc /= 2 ) THEN
  IF (myproc == 0)THEN
      WRITE(6,'(5x,a,i2)') 'ebc reset to 2 for runmod=',runmod
  END IF
      ebc = 2
    END IF
  END IF

  IF ( lbcopt == 1 .AND.                                                &
         (wbc == 5 .OR. ebc == 5 .OR. sbc == 5 .OR. nbc == 5) ) THEN
  IF (myproc == 0)THEN
    WRITE(6,'(5x,a/5x,a,4(/5x,a,i2),2(/5x,a))')                         &
        'The lateral boundary conditions were set to internal ',        &
        'determined, but one of them was set to external forced.',      &
        'wbc = ',wbc,'  ebc = ',ebc,'  sbc = ',sbc,'  nbc = ',nbc,      &
        'Program will try to complete reading in input parameters, ',   &
        'but will stop at the end of subroutine INITPARA.'
  END IF
    err_no = err_no + 1
  ELSE IF ( lbcopt == 2 ) THEN
    IF (myproc == 0)THEN
      WRITE(6,'(5x,a/5x,a)')                                            &
        'The lateral boundary conditions were set to external forced.', &
        'All lateral boundary conditions will be reset to 5 accordingly.'
    END IF
    wbc = 5
    ebc = 5
    sbc = 5
    nbc = 5
  ELSE IF( lbcopt /= 1 ) THEN
    IF (myproc == 0)THEN
      WRITE(6,'(5x,a,i3,a,2(/5x,a))')                                   &
        'Input value of lbcopt = ', lbcopt,' was invalid.',             &
        'Program will try to complete reading in input parameters, ',   &
        'but will stop at the end of subroutine INITPARA.'
    END IF
    err_no = err_no + 1
  END IF

  IF ( ebc == 3 .OR. wbc == 3 ) THEN
    ebc = 3
    wbc = 3
  END IF

  IF ( nbc == 3 .OR. sbc == 3 ) THEN
    nbc = 3
    sbc = 3
  END IF

  IF (myproc == 0)THEN
  WRITE(6,'(5x,a, 4(/5x,i3,a))')                                        &
       'The boundary options are:',                                     &
       wbc,' for  west boundary,',ebc,' for  east boundary,',           &
       sbc,' for south boundary,',nbc,' for north boundary.'


  IF ( vimplct == 1 .AND. (tbc == 2 .OR. bbc == 2) ) THEN
    WRITE(6,'(5x,a/5x,a,2(/5x,a,i2),2(/5x,a))')                         &
        'The small time step integration scheme was set to implicit ',  &
        'which is not valid for periodic vertical boundary conditions.',&
        'tbc = ',tbc,'  bbc = ',bbc,                                    &
        'Program will try to complete reading in input parameters, ',   &
        'but will stop at the end of subroutine INITPARA.'
    err_no = err_no + 1
  END IF

  WRITE(6,'(5x,a, /2(5x,i3,a))')                                        &
       'The boundary options are:',                                     &
       tbc,' for top boundary,',bbc,' for bottom boundary.'

  WRITE(6,'(5x,a, /2(5x,i3,a))')                                        &
       'The upper boundary fft transform option is: ', fftopt

  IF( vimplct == 0.AND.tbc == 4)THEN
    WRITE(6,'(5x,a,/5x,a,/5x,a,2(/5x,a))')                              &
        'The upper radiation condition boundary option was chosen and', &
        'is not compatible with the vertically explicit option',        &
        'reset tbc or vimplct.',                                        &
        'Program will try to complete reading in input parameters, ',   &
        'but will stop at the end of subroutine INITPARA.'
    err_no = err_no + 1
  END IF

  IF( tbc == 4.AND.inibasopt == 2)THEN
    WRITE(6,'(5x,a,/5x,a,4(/5x,a))')                                    &
        'The upper radiation boundary condition option was chosen ',    &
        'with a neutral environment.  This upper radiation ',           &
        'is not compatible with a neutral base state. ',                &
        'Program will try to complete reading in input parameters, ',   &
        'but will stop at the end of subroutine INITPARA.'
    err_no = err_no + 1
  END IF
!
!-----------------------------------------------------------------------
!
!  Input radiation lateral boundary condition options:
!
!-----------------------------------------------------------------------
!
  IF( rbcopt < 1.OR.rbcopt > 5 ) THEN
    WRITE(6,'(1x,a,/1x,a,i3,a,2(5x,a))')                                &
        'Only option rbcopt=1,2,3,4,5 is available in current version.',  &
        'The input was ',rbcopt,', Please reset rbcopt and rerun ARPS.', &
        'Program will try to complete reading in input parameters, ',   &
        'but will stop at the end of subroutine INITPARA.'
    err_no = err_no + 1
  END IF

  WRITE(6,'(5x,a,i4)')                                                  &
       'The radiation boundary condition option was',rbcopt

  WRITE(6,'(5x,a,/5x,a,f13.6)')                                         &
       'The constant gravity phase speed used by radiation ',           &
       'lateral boundary condition option 2 was ',c_phase

  WRITE(6,'(5x,a,f13.6)')                                               &
       'The relaxation coefficient used at the inflow boundaries is'    &
        ,rlxlbc

  WRITE(6,'(5x,a,I4)') 'rbc_plbc = ', rbc_plbc

!  IF ( initopt.ne.1 .and. lbcopt.ne.1 ) THEN
!    pdetrnd = 0
!  ENDIF

  WRITE(6,'(5x,a,i4)')                                                  &
       'Option for pressure detrending was', pdetrnd
  END IF
!
!-----------------------------------------------------------------------
!
!  Input external boundary condition control parameters:
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) READ (unum,exbcpara)
  IF (myproc == 0) WRITE(6,'(a)')'Namelist block exbcpara sucessfully read.'

  CALL mpupdatec(exbcname,80)
  CALL mpupdatec(tinitebd,19)
  CALL mpupdatei(tintvebd,1)
  CALL mpupdatei(ngbrz,1)
  CALL mpupdatei(brlxhw,1)
  CALL mpupdatei(cbcdmp,1)
  CALL mpupdatei(exbcfmt,1)

  IF (myproc == 0)THEN
  WRITE(6,'(5x,a,a9)')                                                  &
      'The initial external boundary time string was ',tinitebd

  WRITE(6,'(5x,a,i10)')                                                 &
      'The time interval to update external boundary conditions was ',  &
      tintvebd

  WRITE(6,'(5x,a,i10)')                                                 &
      'The number of boundary relaxation zone grids was ', ngbrz

  WRITE(6,'(5x,a,e15.8)')                                               &
      'The real grid number where BC relaxation is half weighted was ', &
      brlxhw

  WRITE(6,'(5x,a,e15.8)')                                               &
      'The magnitude of the boundary relaxation damping was ',cbcdmp

  WRITE(6,'(5x,a,i4)')                                                  &
      'The external boundary file format was ', exbcfmt
  END IF
!
!-----------------------------------------------------------------------
!
!  Coriolis parameters:
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) THEN
    READ (unum,coriolis_force,END=100)
    WRITE(6,'(a)')'Namelist block coriolis_force sucessfully read.'

    WRITE (6,'(5x,a,i4)') 'The Coriolis term option was ', coriopt
    WRITE (6,'(5x,a,i4)') 'The earth curvature option was ',            &
                          earth_curvature
    WRITE (6,'(5x,a,i4)') 'The flag for Coriolis formulation was',      &
                          coriotrm
  END IF

  CALL mpupdatei(coriopt,1)
  CALL mpupdatei(earth_curvature,1)
  CALL mpupdatei(coriotrm,1)

!
!-----------------------------------------------------------------------
!
!  Input parameters for turbulent mixing.
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) READ (unum,turbulence,END=100)
  IF (myproc == 0)THEN
  WRITE(6,'(a)')'Namelist block turbulence sucessfully read.'
  END IF
  CALL mpupdatei(tmixopt,1)
  CALL mpupdatei(trbisotp,1)
  CALL mpupdatei(tkeopt,1)
  CALL mpupdatei(trbvimp,1)
  CALL mpupdatei(tmixvert,1)
  CALL mpupdater(tmixcst,1)
  CALL mpupdater(prantl,1)
  CALL mpupdater(kmlimit,1)

  IF (myproc == 0)THEN
    WRITE (6,'(5x,a,i4)') 'The turbulence option was ', tmixopt

    WRITE (6,'(5x,a,i4)') 'The isotropic turbulence option was ', trbisotp

    WRITE (6,'(5x,a,i4)') 'The 1.5 order TKE option was ', tkeopt

    WRITE (6,'(5x,a,i4)')                                               &
          'The implicit treatment of vertical mixing option was',trbvimp

    WRITE (6,'(5x,a,f10.3)')                                            &
         'The nondimensional turbulent prandtl number was ', prantl

    WRITE (6,'(5x,a,f10.3)')                                            &
         'The constant mixing coeff was ', tmixcst

    WRITE (6,'(5x,a,f10.3)')                                            &
         'The parameter used to limit km was ', kmlimit

    WRITE (6,'(5x,a,i4)') 'The PBL parameterization option was ', tmixopt

  END IF

  IF (tmixopt == 4 .AND. (tkeopt <= 0 .OR. tkeopt >= 4) ) THEN
    IF (myproc == 0)THEN
      WRITE (6,'(5x,a,i3,2(/5x,a))')                                      &
          'tkeopt should be 1, 2, or 3 for tmixopt=4, input was=',tkeopt, &
          'Program will try to complete reading in input parameters, ',   &
          'but will stop at the end of subroutine INITPARA.'
    END IF
    err_no = err_no + 1
  END IF
!
!-----------------------------------------------------------------------
!
!  Input control parameters for computational mixing
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) READ (unum,computational_mixing,END=100)
  IF (myproc == 0)THEN
    WRITE(6,'(a)') 'Namelist block computational_mixing sucessfully read.'
  END IF
  CALL mpupdatei(cmix2nd,1)
  CALL mpupdater(cfcm2h,1)
  CALL mpupdater(cfcm2v,1)
  CALL mpupdatei(cmix4th,1)
  CALL mpupdater(cfcm4h,1)
  CALL mpupdater(cfcm4v,1)
  CALL mpupdater(scmixfctr,1)
  CALL mpupdatei(cmix_opt,1)

  IF (myproc == 0)THEN
    WRITE(6,'(5x,a,i4)')                                                &
        'The second order computational mixing option was ',cmix2nd

    WRITE(6,'(5x,a,e15.5)')                                             &
        'The coeff for second order horizontal mixing was ',cfcm2h

    WRITE(6,'(5x,a,e15.5)')                                             &
        'The coeff for second order vertical mixing was ',cfcm2v

    WRITE(6,'(5x,a,i4)')                                                &
        'The fourth order computational mixing option was ',cmix4th

    WRITE(6,'(5x,a,e15.5)')                                             &
        'The coeff for fourth order horizontal mixing was ',cfcm4h

    WRITE(6,'(5x,a,e15.5)')                                             &
        'The coeff for fourth order vertical mixing was ',cfcm4v

    WRITE(6,'(5x,a,e15.5)')                                             &
        'The c-mixing reduction factor for scalars was ',scmixfctr

    WRITE(6,'(5x,a,i4)')                                                &
        'The c-mix monotonic option was ',cmix_opt

    IF( cmix2nd /= 0 .AND. cfcm2h > 0.125/dtbig ) THEN
      WRITE(6,'(5x,a,a,2(/5x,a))')                                      &
          'Value of cfcm2h was too large. ',                            &
          'It has to be less than 1/(8*dtbig).',                        &
          'Program will try to complete reading in input parameters, ', &
          'but will stop at the end of subroutine INITPARA.'
      err_no = err_no + 1
      CALL arpsstop('arpsstop called from INITPARA with cdcm2h selection',1)
    END IF

    IF( cmix2nd /= 0 .AND. cfcm2v > 0.125/dtbig ) THEN
      WRITE(6,'(5x,a,a,2(/5x,a))')                                      &
          'Value of cfcm2v was too large. ',                            &
          'It has to be less than 1/(8*dtbig).',                        &
          'Program will try to complete reading in input parameters, ', &
          'but will stop at the end of subroutine INITPARA.'
      err_no = err_no + 1
    END IF

    IF( cmix4th /= 0 .AND. cfcm4h > 0.125/dtbig ) THEN
      WRITE(6,'(5x,a,a,2(/5x,a))')                                      &
          'Value of cfcm4h was too large. ',                            &
          'It has to be less than 1/(8*dtbig).',                        &
          'Program will try to complete reading in input parameters, ', &
          'but will stop at the end of subroutine INITPARA.'
      err_no = err_no + 1
    END IF

    IF( cmix4th /= 0 .AND. cfcm4v > 0.125/dtbig ) THEN
      WRITE(6,'(5x,a,a,2(/5x,a))')                                      &
          'Value of cfcm4v was too large. ',                            &
          'It has to be less than 1/(8*dtbig).',                        &
          'Program will try to complete reading in input parameters, ', &
          'but will stop at the end of subroutine INITPARA.'
      err_no = err_no + 1
    END IF
  END IF
!
!-----------------------------------------------------------------------
!
!  Input controls for divergence damping on acoustic waves
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) READ (unum,divergence_damping,END=100)
  IF (myproc == 0)THEN
    WRITE(6,'(a)')                                                      &
        'Namelist block divergence_damping sucessfully read.'
  END IF
  CALL mpupdatei(divdmp,1)
  CALL mpupdater(divdmpndh,1)
  CALL mpupdater(divdmpndv,1)

  IF (myproc == 0)THEN
    WRITE (6,'(5x,a,i4)')                                               &
        'The acoustic wave damping option was ', divdmp

    WRITE (6,'(5x,a,f10.3,a,f10.3,a)')                                  &
        'The non-dimensional divergence damping coeff was ',            &
        divdmpndh, ' for horizontal and ',                              &
        divdmpndv, ' for vertical'
  END IF
!
!-----------------------------------------------------------------------
!
!  Rayleigh damping parameters:
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) READ (unum,rayleigh_damping,END=100)
  IF (myproc == 0)THEN
  WRITE(6,'(a)')'Namelist block rayleigh_damping sucessfully read.'
  END IF
  CALL mpupdatei(raydmp,1)
  CALL mpupdater(cfrdmp,1)
  CALL mpupdater(zbrdmp,1)

  IF (myproc == 0)THEN
    WRITE(6,'(5x,a,i4)') 'The rayleigh damping option was ',raydmp

    WRITE(6,'(5x,a,e15.5)') 'The rayleigh damping coeff was ',cfrdmp

    WRITE(6,'(5x,a,e15.5)')                                             &
        'The altitude of base of rayleigh damping was ',zbrdmp

    IF ( raydmp == 2 .AND. lbcopt /= 2 ) THEN
      WRITE (6,'(5x,a,i3,/5x,a,a,2(/5x,a))')                            &
          'You can use raydmp=2 only when lbcopt=2. lbcopt=', lbcopt,   &
          'had been chosen. Please reset raydmp or lbcopt in the',      &
          'input file.',                                                &
          'Program will try to complete reading in input parameters, ', &
          'but will stop at the end of subroutine INITPARA.'
      err_no = err_no + 1
      CALL arpsstop('arpsstop called from INITPARA with raydmp/lbcopt   &
           & selection',1)
    END IF
  END IF
!
!-----------------------------------------------------------------------
!
!  Robert-Asselin time filter coefficient:
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) THEN
    READ (unum,asselin_time_filter,END=100)

    WRITE(6,'(a)')                                                      &
      'Namelist block asselin_time_filter sucessfully read.'
  END IF
  CALL mpupdater(flteps,1)

  IF (myproc == 0)THEN
    WRITE(6,'(5x,a,e15.5)')                                             &
        'The non-dimensional coeff of asselin time filter was ',flteps
  END IF
!
!-----------------------------------------------------------------------
!
!  Input the control parameters for microphysics parameterizations
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) THEN
    READ (unum,microphysics,END=100)
    WRITE(6,'(a)')'Namelist block microphysics sucessfully read.'
  END IF

  CALL mpupdatei(mphyopt,1)
  CALL mpupdatei(moist,1)
  CALL mpupdatei(nmphystp,1)
  CALL mpupdatei(cnvctopt,1)
  CALL mpupdatei(subsatopt,1)
  CALL mpupdatei(kfsubsattrig,1)
  CALL mpupdater(kffbfct,1)
  CALL mpupdater(wcldbs,1)
  CALL mpupdater(confrq,1)
  CALL mpupdater(qpfgfrq,1)
  CALL mpupdatei(idownd,1)

  CALL mpupdatei(dsdpref,1)
  CALL mpupdater(ntcloud,1)
  CALL mpupdater(n0rain,1)
  CALL mpupdater(n0snow,1)
  CALL mpupdater(n0grpl,1)
  CALL mpupdater(n0hail,1)
  CALL mpupdater(rhoice,1)
  CALL mpupdater(rhosnow,1)
  CALL mpupdater(rhogrpl,1)
  CALL mpupdater(rhohail,1)
  CALL mpupdater(alpharain,1)
  CALL mpupdater(alphaice,1)
  CALL mpupdater(alphasnow,1)
  CALL mpupdater(alphagrpl,1)
  CALL mpupdater(alphahail,1)
  CALL mpupdatei(graupel_ON,1)
  CALL mpupdatei(hail_ON,1)
  CALL mpupdatei(mpthermdiag,1)

  CALL mpupdatei(impfallopt,1)
  CALL mpupdatei(fallopt,1)

  CALL mpupdater(rhsat,1)
  CALL mpupdater(rhsatmin,1)
  CALL mpupdater(dx_rhsatmin,1)
  CALL mpupdater(dx_rhsat100,1)

  IF (myproc == 0) THEN
    WRITE (6,'(5x,a,i4)') 'The microphysics option was ', mphyopt
    WRITE (6,'(5x,a,i4,a)') 'And the microphysics will be called every ', &
                          nmphystp, ' time steps.'

    WRITE (6,'(5x,a,I2)'   ) 'DSD preference was ',dsdpref
    WRITE (6,'(5x,a,e15.5)') 'Cloud number concentration was ',ntcloud
    WRITE (6,'(5x,a,e15.5)') 'Rain DSD intercept parameters was ',n0rain
    WRITE (6,'(5x,a,e15.5)') 'Snow DSD intercept parameters was ',n0snow
    WRITE (6,'(5x,a,e15.5)') 'Graupel DSD intercept parameters was ',n0grpl
    WRITE (6,'(5x,a,e15.5)') 'Hail DSD intercept parameters was ',n0hail
    WRITE (6,'(5x,a,e15.5)') 'Ice density was ', rhoice
    WRITE (6,'(5x,a,e15.5)') 'Snow density was ', rhosnow
    WRITE (6,'(5x,a,e15.5)') 'Graupel density was ', rhogrpl
    WRITE (6,'(5x,a,e15.5)') 'Hail density was ', rhohail

    WRITE (6,'(5x,a,e15.5)') 'Rain shape parameter was ', alpharain
    WRITE (6,'(5x,a,e15.5)') 'Ice shape parameter was ', alphaice
    WRITE (6,'(5x,a,e15.5)') 'Snow shape parameter was ', alphasnow
    WRITE (6,'(5x,a,e15.5)') 'Graupel shape parameter was ', alphagrpl
    WRITE (6,'(5x,a,e15.5)') 'Hail shape parameter was ', alphahail

    WRITE (6,'(5x,a,I2)'   ) 'Graupel initiation option was ', graupel_ON
    WRITE (6,'(5x,a,I2)'   ) 'Hail initiation option was ', hail_ON !BJP OCT 2010 hail switch

    WRITE (6,'(5x,a,I2)'   ) 'Microphysics thermal diagnostic option was ',mpthermdiag

    WRITE (6,'(5x,a,i4)') 'The moist phyics option was ', moist


    WRITE (6,'(5x,a,i4)')                                                 &
        'The convective cumulus option was ',cnvctopt

    WRITE (6,'(5x,a,i4)')                                                 &
        'The sub-saturation option was ',subsatopt

    WRITE (6,'(5x,a,f15.5)')                                              &
        'The K-F rainwater feedback option was ',kffbfct

    WRITE (6,'(5x,a,i4)')                                                 &
        'The K-F sub-saturation trigger was kfsubsattrig=', kfsubsattrig

    WRITE (6,'(5x,a,f10.5)') 'The vertical motion was ', wcldbs

    WRITE (6,'(5x,a,f10.5)')                                              &
        'The frequency of conv. para. updated in seconds was',confrq

    WRITE (6,'(5x,a,f10.5)')                                              &
        'The frequency of grid  para. updated in seconds was',qpfgfrq

    WRITE (6,'(5x,a,i4)') 'The downdraft flag was ', idownd

    WRITE (6,'(5x,a,i4)')                                                 &
        'The vertically implicit fall velocity option was ', impfallopt

    WRITE (6,'(5x,a,i4)') 'The fall velocity option was ', fallopt

    WRITE (6,'(5x,a,f10.5)')                                              &
        'The threshold of RH for condensation to occur: rhsat = ',        &
        rhsat

    WRITE (6,'(5x,a,f10.5)')                                              &
        'The threshold of RH for a grid size of dx_rhsatmin: rhsatmin = ',&
        rhsatmin

    WRITE (6,'(5x,a,f15.5)')                                              &
        'The grid size for condensation to occur (RH=rhsatmin): dx_rhsatmin = ', &
        dx_rhsatmin

    WRITE (6,'(5x,a,f15.5)')                                              &
        'The grid size for condensation to occur (RH=100%): dx_rhsat100 = ', &
        dx_rhsat100

    IF (cnvctopt == 1) THEN
      WRITE(6,'(6(5x,a/))')                                             &
           'Option cnvctopt=1 for Kuo convective parameterization ',    &
           'scheme is little tested in the ARPS, and is disabled.',     &
           'It is recommended that you choose a different option. ',    &
           'If you know what you are doing, you can continue to use ',  &
           'this option by commenting out the following ',              &
           'STOP statement in subroutine INITPARA.'
      CALL arpsstop('arpsstop called because cnvctopt == 1.',1)
    END IF

    IF ( moist == 0 .AND. cnvctopt == 1 ) THEN
      WRITE (6,'(5x,a/5x,a,a,2(/5x,a))')                                  &
          'Since cnvctopt = 1, ',                                         &
          'moist has to be set to 1 in order to use cumulus physics',     &
          'Program will try to complete reading in input parameters, ',   &
          'but will stop at the end of subroutine INITPARA.'
      err_no = err_no + 1
      CALL arpsstop('arpsstop called from INITPARA with moist/cnvctopt    &
                    & selection',1)
    END IF

    IF ( cnvctopt == 1 .AND. mphyopt /= 0 ) THEN
      WRITE (6,'(5x,a/5x,a,a,2(/5x,a))')                                  &
          'Use cnvctopt=2 if you wish to use',                            &
          'both cumulus parameterization and microphysics.',              &
          'Program will try to complete reading in input parameters, ',   &
          'but will stop at the end of subroutine INITPARA.'
      err_no = err_no + 1
    END IF

  END IF

  IF ( mphyopt < 0 .OR. mphyopt > 12 ) THEN
    IF (myproc == 0)THEN
      WRITE (6,'(5x,a/5x,a/5x,a)')                                      &
        'No option for mphyopt > 11.',                                  &
        'Program will try to complete reading in input parameters, ',   &
        'but will stop at the end of subroutine INITPARA.'
    END IF
    err_no = err_no + 1
  ELSE IF ( mphyopt == 2 .OR. mphyopt == 3 .OR. mphyopt == 4 ) THEN
    nscalar = 5; nscalarq =  nscalar
    P_QC = 1; P_QR = 2; P_QI = 3; P_QS = 4; P_QH = 5
    qnames(P_QC) = 'qc'; qdescp(P_QC) = 'Cloud water mixing ratio (kg/kg)'
    qnames(P_QR) = 'qr'; qdescp(P_QR) = 'Rain  water mixing ratio (kg/kg)'
    qnames(P_QI) = 'qi'; qdescp(P_QI) = 'Cloud ice   mixing ratio (kg/kg)'
    qnames(P_QS) = 'qs'; qdescp(P_QS) = 'Snow mixing ratio (kg/kg)'
    qnames(P_QH) = 'qh'; qdescp(P_QH) = 'Hail mixing ratio (kg/kg)'
    ice = 1
  ELSE IF ( mphyopt == 5 .OR. mphyopt == 6 .OR. mphyopt == 7 ) THEN

    nscalar = 5; nscalarq =  nscalar
    P_QC = 1; P_QR = 2; P_QI = 3; P_QS = 4; P_QG = 5
    qnames(P_QC) = 'qc'; qdescp(P_QC) = 'Cloud water mixing ratio (kg/kg)'
    qnames(P_QR) = 'qr'; qdescp(P_QR) = 'Rain  water mixing ratio (kg/kg)'
    qnames(P_QI) = 'qi'; qdescp(P_QI) = 'Cloud ice   mixing ratio (kg/kg)'
    qnames(P_QS) = 'qs'; qdescp(P_QS) = 'Snow mixing ratio (kg/kg)'
    qnames(P_QG) = 'qg'; qdescp(P_QG) = 'Graupel mixing ratio (kg/kg)'
    ice = 1

  ELSE IF ( mphyopt == 8 ) THEN
    nscalar = 6; nscalarq = 6
    P_QC = 1; P_QR = 2; P_QI = 3; P_QS = 4; P_QG = 5; P_QH = 6;
    qnames(P_QC) = 'qc'; qdescp(P_QC) = 'Cloud water mixing ratio (kg/kg)'
    qnames(P_QR) = 'qr'; qdescp(P_QR) = 'Rain  water mixing ratio (kg/kg)'
    qnames(P_QI) = 'qi'; qdescp(P_QI) = 'Cloud ice   mixing ratio (kg/kg)'
    qnames(P_QS) = 'qs'; qdescp(P_QS) = 'Snow mixing ratio (kg/kg)'
    qnames(P_QG) = 'qg'; qdescp(P_QG) = 'Graupel mixing ratio (kg/kg)'
    qnames(P_QH) = 'qh'; qdescp(P_QH) = 'Hail mixing ratio (kg/kg)'
    ice = 1
  ELSE IF ( mphyopt == 9 .OR. mphyopt == 10 .OR. mphyopt == 12) THEN
    nscalar = 12; nscalarq = 6
    P_QC = 1; P_QR = 2; P_QI = 3; P_QS =  4; P_QG =  5; P_QH =  6;
    P_NC = 7; P_NR = 8; P_NI = 9; P_NS = 10; P_NG = 11; P_NH = 12;
    qnames(P_QC) = 'qc'; qdescp(P_QC) = 'Cloud water mixing ratio (kg/kg)'
    qnames(P_QR) = 'qr'; qdescp(P_QR) = 'Rain  water mixing ratio (kg/kg)'
    qnames(P_QI) = 'qi'; qdescp(P_QI) = 'Cloud ice   mixing ratio (kg/kg)'
    qnames(P_QS) = 'qs'; qdescp(P_QS) = 'Snow mixing ratio (kg/kg)'
    qnames(P_QG) = 'qg'; qdescp(P_QG) = 'Graupel mixing ratio (kg/kg)'
    qnames(P_QH) = 'qh'; qdescp(P_QH) = 'Hail mixing ratio (kg/kg)'
    qnames(P_NC) = 'nc'; qdescp(P_NC) = 'Cloud water concentrations (#/m3)'
    qnames(P_NR) = 'nr'; qdescp(P_NR) = 'Rain water concentrations (#/m3)'
    qnames(P_NI) = 'ni'; qdescp(P_NI) = 'Cloud ice concentrations (#/m3)'
    qnames(P_NS) = 'ns'; qdescp(P_NS) = 'Snow concentrations (#/m3)'
    qnames(P_NG) = 'ng'; qdescp(P_NG) = 'Graupel concentrations (#/m3)'
    qnames(P_NH) = 'nh'; qdescp(P_NH) = 'Hail concentrations (#/m3)'
    ice = 1
  ELSE IF ( mphyopt == 11) THEN
    nscalar = 17; nscalarq = 6
    P_QC = 1; P_QR =  2; P_QI =  3; P_QS =  4; P_QG =  5; P_QH =  6;
    P_NC = 7; P_NR =  8; P_NI =  9; P_NS = 10; P_NG = 11; P_NH = 12;
              P_ZR = 13; P_ZI = 14; P_ZS = 15; P_ZG = 16; P_ZH = 17;
    qnames(P_QC) = 'qc  '; qdescp(P_QC) = 'Cloud water mixing ratio (kg/kg)'
    qnames(P_QR) = 'qr  '; qdescp(P_QR) = 'Rain  water mixing ratio (kg/kg)'
    qnames(P_QI) = 'qi  '; qdescp(P_QI) = 'Cloud ice   mixing ratio (kg/kg)'
    qnames(P_QS) = 'qs  '; qdescp(P_QS) = 'Snow mixing ratio (kg/kg)'
    qnames(P_QG) = 'qg  '; qdescp(P_QG) = 'Graupel mixing ratio (kg/kg)'
    qnames(P_QH) = 'qh  '; qdescp(P_QH) = 'Hail mixing ratio (kg/kg)'
    qnames(P_NC) = 'nc  '; qdescp(P_NC) = 'Cloud water concentrations (#/kg)'
    qnames(P_NR) = 'nr  '; qdescp(P_NR) = 'Reain water concentrations (#/kg)'
    qnames(P_NI) = 'ni  '; qdescp(P_NI) = 'Cloud ice concentrations (#/kg)'
    qnames(P_NS) = 'ns  '; qdescp(P_NS) = 'Snow concentrations (#/kg)'
    qnames(P_NG) = 'ng  '; qdescp(P_NG) = 'Graupel concentrations (#/kg)'
    qnames(P_NH) = 'nh  '; qdescp(P_NH) = 'Hail concentrations (#/kg)'
    qnames(P_ZR) = 'zr  '; qdescp(P_ZR) = 'Rain reflectivity (m6/kg)'
    qnames(P_ZI) = 'zi  '; qdescp(P_ZI) = 'Ice reflectivity (m6/kg)'
    qnames(P_ZS) = 'zs  '; qdescp(P_ZS) = 'Snow reflectivity (m6/kg)'
    qnames(P_ZG) = 'zg  '; qdescp(P_ZG) = 'Graupel reflectivity (m6/kg)'
    qnames(P_ZH) = 'zh  '; qdescp(P_ZH) = 'Hail reflectivity (m6/kg)'
    ice = 1
  ELSE
    nscalar = 2; nscalarq = 2
    P_QC = 1; P_QR = 2
    qnames(P_QC) = 'qc'; qdescp(P_QC) = 'Cloud water mixing ratio (kg/kg)'
    qnames(P_QR) = 'qr'; qdescp(P_QR) = 'Rain  water mixing ratio (kg/kg)'
    ice = 0
  END IF

  IF (kfsubsattrig < 0 .OR. kfsubsattrig > 1) THEN
    WRITE (6,'(5x,a/5x,a/5x,a)')                                        &
      'ERROR: No option for kfsubsattrig < 0 or > 1. ',                 &
      'Program will try to complete reading in input parameters, ',     &
      'but will stop at the end of subroutine INITPARA.'
    err_no = err_no + 1
  END IF

  IF ( subsatopt < 0 .OR. subsatopt > 2 ) THEN
    WRITE (6,'(5x,a/5x,a/5x,a)')                                        &
      'ERROR: No option for subsatopt < 0 or > 2. ',                    &
      'Program will try to complete reading in input parameters, ',     &
      'but will stop at the end of subroutine INITPARA.'
    err_no = err_no + 1
  ELSE IF (subsatopt == 0) THEN
    rhsat = 1.0
  ELSE IF (subsatopt == 1) THEN
    IF (rhsat < 0.) THEN
      WRITE (6,'(5x,a/5x,a/5x,a/5x,a)')                                 &
        'ERROR: The threshold of RH can not be less than 0. ',          &
        'But you have chosen rhsat < 0. ',                              &
        'Program will try to complete reading in input parameters, ',   &
        'but will stop at the end of subroutine INITPARA.'
      err_no = err_no + 1
    ELSE IF (rhsat > 1.) THEN
      WRITE (6,'(5x,a/5x,a)')                                           &
        'WARNING: The subsatopt is designed for RH less than 1. ',      &
        'But you have chosen rhsat > 1. It is now set to 1. '
    END IF
    rhsat = max(0., min(1.0, rhsat))
  ELSE IF (subsatopt == 2) THEN
    IF (rhsatmin < 0.) THEN
      WRITE (6,'(5x,a/5x,a/5x,a)')                                      &
        'ERROR: The threshold of RH can not be less than 0. ',          &
        'But you have chosen rhsatmin < 0. ',                           &
        'Program will try to complete reading in input parameters, ',   &
        'but will stop at the end of subroutine INITPARA.'
      err_no = err_no + 1
    ELSE IF (rhsatmin > 1.) THEN
      WRITE (6,'(5x,a/5x,a/5x,a)')                                      &
        'WARNING: the subsatopt is designed for RH less than 1. ',      &
        'But you have chosen rhsatmin greater than 1. ',                &
        'It is now re-set to 1. '
      rhsatmin = 1.
    END IF
    IF (dx_rhsatmin < 0. .OR. dx_rhsat100 < 0.) THEN
      WRITE (6,'(5x,a/5x,a/5x,a)')                                      &
        'ERROR: dx_rhsatmin or dx_rhsat100 can not be less than 0. ',   &
        'Program will try to complete reading in input parameters, ',   &
        'but will stop at the end of subroutine INITPARA.'
      err_no = err_no + 1
    ELSE IF (dx_rhsatmin < dx_rhsat100) THEN
      WRITE (6,'(5x,a/5x,a/5x,a)')                                      &
        'ERROR: dx_rhsatmin can not be greater than dx_rhsat100. ',     &
        'Program will try to complete reading in input parameters, ',   &
        'but will stop at the end of subroutine INITPARA.'
      err_no = err_no + 1
    END IF
    rhsat = max(rhsatmin,                                               &
              min(1.0,1.0+(rhsatmin-1.0)*(dx-dx_rhsat100)               &
                              /max(0.1,dx_rhsatmin-dx_rhsat100)))
  ENDIF

  IF (myproc == 0) WRITE (6,'(5x,a,f10.5)')                             &
      'rhsat for model integration is re-adjusted to ', rhsat
!
!michi
!-----------------------------------------------------------------------
!
!  Input the control parameters for concentration parameterizations
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) THEN
    READ (unum,concentration,END=100)
    WRITE(6,'(a)')'Namelist block concentration sucessfully read.'
  END IF

  CALL mpupdatei(ccin,1)
  CALL mpupdatei(cpoint,1)

  IF (cpoint > 0) THEN   ! TINA - to accommodate bubble cases
    CALL mpupdatei(icc,cpoint)
    CALL mpupdatei(jcc,cpoint)
    CALL mpupdatei(kcc,cpoint)
    CALL mpupdater(ccemit,cpoint)
    CALL mpupdater(ccstart,cpoint)
    CALL mpupdater(ccend,cpoint)
  ELSE                   ! for bubble case
    CALL mpupdater(ccemit,1)
  END IF

  IF (myproc == 0) THEN
    WRITE (6,'(5x,a,i4)') 'The concentration option was ', ccin
    WRITE (6,'(5x,a,i4)') 'Number of release points was ', cpoint
    IF (cpoint < 0) THEN
      WRITE (6,'(5x,a,i4)') 'Bubble shaped initial conc. field'
    ELSE
      WRITE (6,'(5x,a,i4)') 'Number of emission points ',cpoint
      WRITE(6,'(5x,a)')     '     Release point  Emitted     Emitted time'
      WRITE(6,'(5x,a)')     ' No.  icc  jcc  kcc Concentrat. sta. end. '
      WRITE(6,'(5x,a)')     '---- ---- ---- ---- ----------- ---- ---- '
      DO i = 1, cpoint      ! Y. Wang reformatted
        WRITE(6,'(5x,4i5,f10.5,2i5)')    i,icc(i),jcc(i),kcc(i),ccemit(i),ccstart(i),ccend(i)
        ! TINA or the initial concentration for a bubble
      END DO
    END IF
  ENDIF
!
!michi
!
  IF (ccin > 0) THEN
    nscalar = nscalar + 1
    P_CC = nscalar
    qnames(P_CC) = 'cc'; qdescp(P_CC) = 'Concentration (-)'
  END IF
!
!-----------------------------------------------------------------------
!
!  Input the control parameters for radiation parameterizations
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) THEN
    READ (unum,radiation,END=100)
    WRITE(6,'(a)')'Namelist block radiation sucessfully read.'
  END IF

  CALL mpupdatei(radopt,1)
  CALL mpupdatei(radstgr,1)
  CALL mpupdatei(rlwopt,1)
  CALL mpupdatei(radshade,1)    ! augustin add radshade
  CALL mpupdater(dtrad,1)
  CALL mpupdatei(raddiag,1)

  IF (myproc == 0)THEN
    WRITE (6,'(5x,a,i4)')                                               &
        'The radiation phyics option was ', radopt

    WRITE (6,'(5x,a,i4)')                                               &
        'The radiation staggering option was ', radstgr

    WRITE (6,'(5x,a,i4)')                                               &
        'The option for longwave schemes was ', rlwopt
! augustin add radshade
  WRITE (6,'(5x,a,i4)')                                                 &
      'The option to compute the shade was ', radshade

    WRITE (6,'(5x,a,f10.5)')                                            &
        'The Time interval to update the radiation forcing was ', dtrad

    WRITE (6,'(5x,a,i4)')                                               &
        'The radiation diagnostic output option was ', raddiag
  END IF

  nradstp = nint( dtrad/dtbig )
  IF( nradstp /= 0 ) THEN
    IF (myproc == 0) WRITE(6,'(5x,a,i6,a)')                             &
        'Radiation physics will be calculated every ', nradstp,         &
        ' time steps'
    dtrad = dtbig*nradstp
  ELSE
    nradstp = -1
    radopt  = 0
    IF (myproc == 0)  WRITE(6,'(5x,a)')      &
                                'Radiation physics is switched off.'
  END IF
!
!-----------------------------------------------------------------------
!
!  Input surface physics options
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) THEN
    READ (unum,surface_physics,END=100)
    WRITE(6,'(a)') 'Namelist block surface_physics sucessfully read.'
  END IF
  CALL mpupdatei(sfcphy,1)
  CALL mpupdatei(landwtr,1)
  CALL mpupdatei(cdhwtropt,1)
  CALL mpupdater(cdmlnd,1)
  CALL mpupdater(cdmwtr,1)
  CALL mpupdater(cdhlnd,1)
  CALL mpupdater(cdhwtr,1)
  CALL mpupdater(cdqlnd,1)
  CALL mpupdater(cdqwtr,1)
  CALL mpupdatei(pbldopt,1)
  CALL mpupdater(pbldpth0,1)
  CALL mpupdater(lsclpbl0,1)
  CALL mpupdatei(tqflxdis,1)
  CALL mpupdater(dtqflxdis,1)
  CALL mpupdatei(smthflx,1)
  CALL mpupdatei(numsmth,1)
  CALL mpupdatei(sfcdiag,1)

  IF (myproc == 0)THEN
    WRITE (6,'(5x,a,i4)') 'The surface physics option was ', sfcphy
    WRITE (6,'(5x,a,i4)') 'The land/water option was ', landwtr
    WRITE (6,'(5x,a,i4)') 'The constant water cdh option was ', cdhwtropt
    WRITE (6,'(5x,a,f10.3)')                                            &
      'The user specified drag coeff for momentun over land was ',cdmlnd
    WRITE (6,'(5x,a,f10.3)')                                            &
      'The user specified drag coeff for momentun over water was ',cdmwtr
    WRITE (6,'(5x,a,f10.3)')                                              &
        'The user specified drag coeff for heat over land was ',          &
        cdhlnd
    WRITE (6,'(5x,a,f10.3)')                                              &
        'The user specified drag coeff for heat over water was ',         &
        cdhwtr
    WRITE (6,'(5x,a,f10.3)')                                              &
        'The user specified drag coeff for moisture over land was ',      &
        cdqlnd
    WRITE (6,'(5x,a,f10.3)')                                              &
        'The user specified drag coeff for moisture over water was ',     &
        cdqwtr

    WRITE (6,'(5x,a,i3)')                                                 &
        'The option for determining PBL depth was ', pbldopt

    WRITE (6,'(5x,a,f10.3)')                                              &
        'The user specified PBL depth was ', pbldpth0

    WRITE (6,'(5x,a,f10.3)')                                              &
        'The PBL length scale ', lsclpbl0

    WRITE (6,'(5x,a,i4)') 'The flux distribution option was ',            &
                         sflxdis
  END IF

  IF ( sfcphy == 0 ) THEN
    sflxdis = 0
    IF (myproc == 0) WRITE (6,'(5x,a/5x,a)')                            &
        'When sfcphy=0, there is no surface flux to be distributed.',   &
        ' Set sflxdis=0.'
  ELSE IF ( sflxdis < 0 .OR. sflxdis > 3 ) THEN
    IF (myproc == 0) WRITE (6,'(5x,a,/5x,a,/5x,a)')                     &
        'The options for sflxdis must be between 0 and 3.',             &
        'Program will try to complete reading in input parameters, ',   &
        'but will stop at the end of subroutine INITPARA.'
    err_no = err_no + 1
  END IF

  IF (myproc == 0) WRITE (6,'(5x,a,i4)')                                &
      'The heat and moisture distribution option was ',tqflxdis

  IF ( sfcphy == 0 ) THEN
    tqflxdis = 0
    IF (myproc == 0) WRITE (6,'(5x,a/5x,a)')                            &
        'When sfcphy=0, there is no surface flux to be distributed.',   &
        ' Set tqflxdis=0.'
  ELSE IF ( tqflxdis < 0 .OR. tqflxdis > 2 ) THEN
  IF (myproc == 0)THEN
    WRITE (6,'(5x,a,/5x,a,/5x,a)')                                      &
        'The options for tqflxdis must be 0, 1, or 2.',                 &
        'Program will try to complete reading in input parameters, ',   &
        'but will stop at the end of subroutine INITPARA.'
  END IF
    err_no = err_no + 1
  END IF

  IF( tqflxdis /= 0 .AND. sflxdis /= 0 ) THEN
  IF (myproc == 0)THEN
    WRITE (6,'(5x,a,/5x,a,/5x,a,/5x,a)')                                &
        'Options tqflxdis and sflxdis should not be turned on at the',  &
        'same time. Please turn one of them off',                       &
        'Program will try to complete reading in input parameters, ',   &
        'but will stop at the end of subroutine INITPARA.'
  END IF
    err_no = err_no + 1
  END IF

  IF ( smthflx >= 0 ) THEN
    numsmth = MAX( 1, numsmth )
  END IF

  IF ( radopt == 0 .AND. (sfcphy == 3 .OR. sfcphy == 4) ) THEN
    radopt = 1
    IF (myproc == 0) WRITE (6,'(5x,a/5x,a/5x,a,i2)')                    &
        'Since soil-vegetation process was switched on, we must',       &
        'compute the surface radiation flux for energy balance.',       &
        'radopt was reset to 1 in INITPARA.'
  END IF
!
!-----------------------------------------------------------------------
!
!  Input soil and vegetation parameters
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) THEN
    READ (unum,soil_ebm,END=100)
    WRITE(6,'(a)')'Namelist block soil_ebm sucessfully read.'
  END IF

  CALL mpupdatei(sfcdat,1)
  CALL mpupdatei(soilinit,1)
  CALL mpupdater(dtsfc,1)
  CALL mpupdatei(styp,1)
  CALL mpupdatei(vtyp,1)
  CALL mpupdater(lai0,1)
  CALL mpupdater(roufns0,1)
  CALL mpupdater(veg0,1)
  CALL mpupdatei(nzsoil,1)
  CALL mpupdatei(soilmodel_option,1)
  CALL mpupdater(dzsoil,1)
  CALL mpupdater(zrefsoil,1)
  CALL mpupdater(tsoilint,nzsoil)
  CALL mpupdater(qsoilint,nzsoil)
  CALL mpupdater(soildzmin,1)
  CALL mpupdater(soildlayer1,1)
  CALL mpupdater(soildlayer2,1)
  CALL mpupdater(soilstrhtune,1)
  CALL mpupdater(ptslnd0,1)
  CALL mpupdater(ptswtr0,1)
  CALL mpupdater(wetcanp0,1)
  CALL mpupdater(snowdpth0,1)
  CALL mpupdater(ttprt,1)
  CALL mpupdater(tbprt,1)
  CALL mpupdater(wgrat,1)
  CALL mpupdater(w2rat,1)
  CALL mpupdatec(sfcdtfl, 256)
  CALL mpupdatec(soilinfl,256)
  CALL mpupdatei(sfcfmt,1)
  CALL mpupdatei(soilfmt,1)
  CALL mpupdatei(nstyp,1)
  CALL mpupdatei(tsoil_offset, 1)
  CALL mpupdater(tsoil_offset_amplitude, 1)
  CALL mpupdatei(prtsoilflx, 1)

!EMK 15 July 2002
  CALL mpupdatei(soilmodel_forced,1)
  CALL mpupdatec(sitemeso,256)
  CALL mpupdatec(siteflux,256)
  CALL mpupdatec(siternet,256)
  CALL mpupdatec(sitesoil,256)
  CALL mpupdatec(siteveg, 256)
  CALL mpupdatei(soilstrhopt,1)
!EMK END 15 July 2002

  nstyp = MAX(1,nstyp)
  nstyps = nstyp

  IF( soilstrhopt == 0 .AND. soildzmin /= dzsoil ) THEN
    IF (myproc == 0)THEN
      WRITE(6,'(5x,a)')                                                   &
         'For non-stretched case, dzmin was reset to dz.'
    END IF
    soildzmin = dzsoil
  END IF

  IF (myproc == 0)THEN

  WRITE(6,'(5x,a,i4)') 'The stretch option was ', soilstrhopt

  WRITE(6,'(5x,a,i5,a)')                                             &
       'Input nzsoil was ',nzsoil,' levels'

  WRITE(6,'(5x,a,f10.3,a)')                                             &
       'Input dzsoil was ',dzsoil,' meters'

  WRITE(6,'(5x,a,f10.3,a)')                                             &
       'Input zrefsoil was ',zrefsoil,' meters'

  WRITE(6,'(5x,a,f10.3,a)')                                             &
       'soildzmin is ',soildzmin,' meters'

  WRITE(6,'(5x,a,f10.3,a)')                                             &
       'soildlayer1 is ',soildlayer1 ,' meters'

  WRITE(6,'(5x,a,f10.3,a)')                                             &
       'soildlayer2 is ',soildlayer2 ,' meters'
  END IF


  IF (myproc == 0)THEN
  WRITE (6,'(5x,a,i4)')                                                 &
       'The surface data input option was ',sfcdat

  WRITE (6,'(5x,a,i4)')                                                 &
       'The soil model forced option was ',soilmodel_forced

  WRITE (6,'(5x,a,i4)')                                                 &
       'The soil scheme input option was ',soilmodel_option

  WRITE (6,'(5x,a,i4)')                                                 &
       'The surface initial data input option was ',soilinit
  END IF

  IF ( sfcphy == 0 ) THEN
    sfcdat = 1
    soilinit = 1
    landwtr = 0
    sfcdiag = 0
    IF (myproc == 0)       &
    WRITE (6,'(5x,a/5x,a/5x,a)')                                        &
        'Since sfcphy = 0, sfcdat and soilinit are set to 1 and',       &
        'landwtr to 0 to avoid reading the surface data.',              &
        'Diagnostics printing is turned off.'
  END IF

  IF (soilmodel_option == 1) THEN
    nzsoil = 2
    dzsoil = 1.0
    IF (myproc == 0)   &
      WRITE(6,'(5x,a)')                                                &
        'Since soilmodel_option = 1, nzsoil is set to 2 and dzsoil is set to 1.0'
  END IF

  IF (myproc == 0)THEN
  WRITE (6,'(5x,a,f10.3)')                                              &
       'The time step for surface energy budget model was ',dtsfc

  WRITE (6,'(5x,a,i4)') 'The surface soil type is ',styp

  WRITE (6,'(5x,a,i4)') 'The surface vegtation type is ',vtyp

  WRITE (6,'(5x,a,f10.3)') 'The leaf area index is ', lai0

  WRITE (6,'(5x,a,f10.3)')                                              &
       'The user specified land surface roughness was ', roufns0

  WRITE (6,'(5x,a,f10.3)')                                              &
       'The user specified vegetation fraction was ', veg0



  WRITE (6,'(5x,a,i4)') 'The option for the soil scheme being forced is ', &
        soilmodel_forced

  WRITE (6,'(5x,a,i4)') 'The soil scheme is ',soilmodel_option

  WRITE (6,'(5x,a,i4)') 'The number of soil layers are ',nzsoil

  WRITE (6,'(5x,a,f10.3)') 'The average soil layer depth is ',dzsoil

  WRITE (6,'(5x,a,f10.3)') 'The reference height of the soil depth is ',zrefsoil

  IF(soilinit == 1)THEN
    DO i=1,nzsoil
      WRITE (6,'(5x,a,f10.3)') 'The profile soil temperature is ', &
               tsoilint(i)

      WRITE (6,'(5x,a,f10.3)') 'The profile soil moisture is ',    &
               qsoilint(i)
    END DO
  END IF

  WRITE (6,'(5x,a,a,f10.3)')                                            &
       'The initial ground level soil potential temperature ',          &
       'over land is ',ptslnd0

  WRITE (6,'(5x,a,a,f10.3)')                                            &
       'The initial ground level soil potential temperature ',          &
       'over water is ',ptswtr0

  WRITE (6,'(5x,a,f10.3)') 'The canopy moisture is ',wetcanp0

  WRITE (6,'(5x,a,f10.3)') 'The snow depth is ',snowdpth0

  WRITE (6,'(5x,a,f10.3)')                                              &
      'The offset of top soil layer from surface air temperature is ',ttprt

  WRITE (6,'(5x,a,f10.3)')                                              &
      'The offset of bottom soil layer from surface air temperature is ',tbprt

  WRITE (6,'(5x,a,f10.3)')                                              &
      'The saturation ratio of surface soil moisture is ',wgrat

  WRITE (6,'(5x,a,f10.3)')                                              &
      'The saturation ratio of deep soil moisture is ',w2rat
  END IF

  IF(soilmodel_forced == 1)THEN

    lenstr = 256
    CALL strlnth( sitemeso, lenstr )
    IF (myproc == 0) WRITE(6,'(5x,a,a)')                                &
         'Surface data file to be used is ', sitemeso(1:lenstr)

    lenstr = 256
    CALL strlnth( siteflux, lenstr )
    IF (myproc == 0) WRITE(6,'(5x,a,a)')                                &
         'Surface data file to be used is ', siteflux(1:lenstr)

    lenstr = 256
    CALL strlnth( siternet, lenstr )
    IF (myproc == 0) WRITE(6,'(5x,a,a)')                                &
         'Surface data file to be used is ', siternet(1:lenstr)

    lenstr = 256
    CALL strlnth( sitesoil, lenstr )
    IF (myproc == 0) WRITE(6,'(5x,a,a)')                                &
         'Surface data file to be used is ', sitesoil(1:lenstr)

    lenstr = 256
    CALL strlnth( siteveg, lenstr )
    IF (myproc == 0) WRITE(6,'(5x,a,a)')                                &
         'Surface data file to be used is ', siteveg(1:lenstr)

    lenstr = 256
    CALL strlnth( sfcdtfl, lenstr )
    IF (myproc == 0) THEN
      WRITE(6,'(5x,a,a)')                                               &
         'Surface data file to be used is ', sfcdtfl(1:lenstr)

      WRITE (6,'(5x,a,i4)')                                             &
         'The surface data file format was ',sfcfmt
    END IF

    lenstr = 256
    CALL strlnth( soilinfl, lenstr )
    IF (myproc == 0)THEN
      WRITE(6,'(5x,a,a)')                                               &
         'Soil data file to be used is ', soilinfl(1:lenstr)

      WRITE (6,'(5x,a,i4)')                                             &
         'The soil data file format was ',soilfmt
    END IF

  END IF ! end of soilmodel_forced if block....

  IF ( sfcphy /= 0 ) THEN

    dtsfc = MIN( dtsfc, dtbig )
    dtsfc0 = dtsfc
    nsfcst = nint(dtbig/dtsfc)
    dtsfc  = dtbig/nsfcst

    IF ( dtsfc > dtsfc0 ) THEN
      nsfcst = nsfcst + 1
      dtsfc = dtbig/nsfcst
  IF (myproc == 0)THEN
      WRITE(6,'(/a,1x,f10.4,a)')                                        &
          '     The actual surface model time step to be used is ',     &
           dtsfc, ' seconds'
  END IF
    END IF

    IF (styp <= 0 .OR. styp >= 14) THEN
  IF (myproc == 0)THEN
      WRITE(6,'(5x,a,i3,a,2(/5x,a))')                                   &
          'The input styp =',styp, 'not acceptable.',                   &
          'Program will try to complete reading in input parameters, ', &
          'but will stop at the end of subroutine INITPARA.'
  END IF
      err_no = err_no + 1
    END IF

    IF (vtyp <= 0 .OR. vtyp >= 15) THEN
  IF (myproc == 0)THEN
      WRITE(6,'(5x,a,i3,a,2(/5x,a))')                                   &
          'The input vtyp =',vtyp, 'not acceptable.',                   &
          'Program will try to complete reading in input parameters, ', &
          'but will stop at the end of subroutine INITPARA.'
  END IF
      err_no = err_no + 1
    END IF

    IF ( lai0 < 0.0 ) THEN
  IF (myproc == 0)THEN
      WRITE(6,'(5x,a,f6.2,a,2(/5x,a))')                                 &
          'The input lai0=',lai0, 'not acceptable.',                    &
          'Program will try to complete reading in input parameters, ', &
          'but will stop at the end of subroutine INITPARA.'
  END IF
      err_no = err_no + 1
    END IF

    IF ( roufns0 < 0.0 ) THEN
  IF (myproc == 0)THEN
      WRITE(6,'(5x,a,f6.2,a,2(/5x,a))')                                 &
          'The input roufns0=',roufns0, 'not acceptable.',              &
          'Program will try to complete reading in input parameters, ', &
          'but will stop at the end of subroutine INITPARA.'
  END IF
      err_no = err_no + 1
    END IF

    IF ( veg0 < 0.0 .OR. veg0 > 1.0 ) THEN
  IF (myproc == 0)THEN
      WRITE(6,'(5x,a,f6.2,a,2(/5x,a))')                                 &
          'The input veg0=',veg0, 'not acceptable.',                    &
          'Program will try to complete reading in input parameters, ', &
          'but will stop at the end of subroutine INITPARA.'
  END IF
      err_no = err_no + 1
    END IF

  IF ( nzsoil < 1 ) THEN
   IF (myproc == 0)THEN
      WRITE(6,'(5x,a,f6.2,a,2(/5x,a))')                                 &
          'The input nzsoil=',nzsoil, 'not acceptable.',                &
          'Program will try to complete reading in input parameters, ', &
          'but will stop at the end of subroutine INITPARA.'
   END IF
   err_no = err_no + 1
  END IF

  IF(soilinit == 1)THEN
    DO i=1,nzsoil

    IF ( tsoilint(i) < 173.13 .OR. tsoilint(i) > 373.16 ) THEN
      IF (myproc == 0)THEN
      WRITE(6,'(5x,a,f10.5,a,2(/5x,a))')                                &
          'The input tsoilint=',tsoilint(i), 'not acceptable.',         &
          'Program will try to complete reading in input parameters, ', &
          'but will stop at the end of subroutine INITPARA.'
      END IF
      err_no = err_no + 1
    END IF

    END DO
  END IF


    IF ( ptslnd0 < 173.13 .OR. ptslnd0 > 373.16 ) THEN
  IF (myproc == 0)THEN
      WRITE(6,'(5x,a,f10.5,a,2(/5x,a))')                                &
          'The input ptslnd0=',ptslnd0, 'not acceptable.',              &
          'Program will try to complete reading in input parameters, ', &
          'but will stop at the end of subroutine INITPARA.'
  END IF
      err_no = err_no + 1
    END IF

    IF ( ptswtr0 < 173.13 .OR. ptswtr0 > 373.16 ) THEN
  IF (myproc == 0)THEN
      WRITE(6,'(5x,a,f10.5,a,2(/5x,a))')                                &
          'The input ptswtr0=',ptswtr0, 'not acceptable.',              &
          'Program will try to complete reading in input parameters, ', &
          'but will stop at the end of subroutine INITPARA.'
  END IF
      err_no = err_no + 1
    END IF

    IF ( styp == 12 ) THEN
      ptslnd0  = MIN( ptslnd0, 273.16 )
      DO i=1,nzsoil
        tsoilint(i) = MIN( tsoilint(i), 273.16 )
        qsoilint(i) = 0.
      END DO

      wetcanp0 = 0.
  IF (myproc == 0)THEN
      WRITE(6,'(/a/a/a)')                                               &
          '     The soil type is ice.',                                 &
          '     The soil temperatures are set to ice point 273.16 K.',  &
          '     And the moisture variables are set to 0.'
  END IF
    ELSE IF ( styp == 13 ) THEN
      DO i=1,nzsoil
      qsoilint(i) = 1.
      wetcanp0 = 1.
  IF (myproc == 0)THEN
      WRITE(6,'(/a/a)')                                                 &
          '     The soil type is water.',                               &
          '     The moisture variables are set to 1.'
  END IF
      END DO

    ELSE
      DO i=1,nzsoil
      IF (qsoilint(i) < 0.0 .OR. qsoilint(i) > 1.0) THEN
  IF (myproc == 0)THEN
        WRITE(6,'(5x,a,f10.5,a,2(/5x,a))')                              &
            'The input qsoil=',qsoilint(i), 'not acceptable.',            &
            'Program will try to complete reading in input parameters, ', &
            'but will stop at the end of subroutine INITPARA.'
  END IF
        err_no = err_no + 1
      ELSE IF ( qsoilint(i) > wsat(styp) ) THEN
        qsoilint(i) = wsat(styp)
  IF (myproc == 0)THEN
        WRITE(6,'(/a/a,f10.4)')                                         &
            '     The input qsoil is greater than the saturated value', &
            '     and actually adjusted to the saturated value: ',      &
             qsoilint(i)
  END IF
      END IF
      END DO



      wgrat = MAX( 0.0, MIN(1.0,wgrat) )
      w2rat = MAX( 0.0, MIN(1.0,w2rat) )

      wrmax = .2*veg0*lai0
      IF (wetcanp0 < 0.0 .OR. wetcanp0 > 1.0) THEN
  IF (myproc == 0)THEN
        WRITE(6,'(5x,a,f10.5,a,2(/5x,a))')                              &
            'The input wetcanp0=',wetcanp0, 'not acceptable.',          &
            'Program will try to complete reading in input parameters, ', &
            'but will stop at the end of subroutine INITPARA.'
  END IF
        err_no = err_no + 1
      ELSE IF ( wetcanp0 > wrmax ) THEN
        wetcanp0 = wrmax
  IF (myproc == 0)THEN
        WRITE(6,'(/a/a,f10.4)')                                         &
            '     The input wetcanp0 is greater than the maximun value ', &
            '     and actually adjusted to the maximun value: ',wetcanp0
  END IF
      END IF

    END IF

  END IF

  nstyps = MAX(1,nstyps)
  IF (myproc == 0)THEN
  WRITE(6,'(5x,a,i5)') "nstyps =",nstyps
  END IF

!
! soil_offset amplitude
!
  SELECT CASE (tsoil_offset)
    CASE (0:2)
        IF (myproc == 0) &
          WRITE(6, '(a/, a, I2)')                                           &
                '     Option for including seasonal deep and skin'//        &
                ' layer temperature offset in the two-layer soil model',    &
                '     Your choice is ', tsoil_offset
    CASE DEFAULT
       err_no = err_no + 1
       IF (myproc == 0) &
         WRITE(6, '(a, I2, 3(a/))')                                          &
            '     The input tsoil_offset =',tsoil_offset, 'not acceptable.', &
            '     Program will try to complete reading in input parameters,',&
            '     but will stop at the end of subroutine INITPARA.'
  END SELECT


  IF (myproc == 0) &
     WRITE(6, '(a/, a, F5.3)')                                                &
          '    The amplitude of the annual cycle of the difference (offset)', &
          '    in deep and skin layer soil seasonal-mean temperatures is ',   &
          tsoil_offset_amplitude

!
!-----------------------------------------------------------------------
!
!  Read in parameters for automatic grid translation.
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) THEN
    READ (unum,grdtrans,END=100)
    WRITE(6,'(a)')'Namelist block grdtrans sucessfully read.'
  END IF

  CALL mpupdatei(cltkopt,1)
  CALL mpupdatei(grdtrns,1)
  CALL mpupdater(umove,1)
  CALL mpupdater(vmove,1)
  CALL mpupdater(chkdpth,1)
  CALL mpupdater(twindow,1)
  CALL mpupdater(tceltrk,1)
  CALL mpupdater(tcrestr,1)

  IF (myproc == 0)THEN
    WRITE(6,'(5x,a,i10)') 'Cell tracking option was ',cltkopt
    WRITE(6,'(5x,a,i10)') 'The grid translation option was ',grdtrns
  END IF

  IF( grdtrns == 2 .AND. cltkopt == 0 ) THEN
    cltkopt = 1
    IF (myproc == 0) WRITE(6,'((5x,a)/)')                               &
                           'Since grdtrns =2, cltkopt was reset to 1.'
  END IF
!
!-----------------------------------------------------------------------
!
!  Ground-relative domain translation speed:
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0)THEN
  WRITE(6,'(5x,a,f10.3)')                                               &
       'The domain translation speed in x was ', umove

  WRITE(6,'(5x,a,f10.3)')                                               &
       'The domain translation speed in y was ', vmove


  WRITE(6,'(5x,a,f10.3)')                                               &
      'The depth of domain to check for grid translation was',chkdpth
  WRITE(6,'(5x,a,f10.3)')                                               &
      'The time window for updating umove and vmove was ',twindow

  WRITE(6,'(5x,a,f10.3)') 'Cell tracking interval was ', tceltrk

  END IF
  IF( tceltrk > 0.0 .AND. tceltrk < dtbig ) THEN
    nceltrk = 1
  ELSE
    nceltrk = nint(tceltrk/dtbig)
  END IF

  IF (myproc == 0) THEN
    WRITE(6,'(5x,a,i6,a)')                                              &
       'Cell-tracking routine will be called every', nceltrk,           &
       ' time steps.'

    WRITE(6,'(5x,a,f10.3)') 'Cell restore time was ', tcrestr
  END IF
!
!-----------------------------------------------------------------------
!
!  Read in namelist &history_dump
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) THEN
    READ (unum,history_dump,END=100)
    WRITE(6,'(a)')'Namelist block history_dump sucessfully read.'
  END IF

  DO i = FINDX_NUM,1,-1
    idummy = dmp_out_joined/10**(i-1)
    joindmp(i) = MOD(idummy,2)
  END DO

  CALL mpupdatei(hdmpopt,1)
  CALL mpupdatei(hdmpfmt,1)
  CALL mpupdatei(joindmp,FINDX_NUM)
  CALL mpupdatei(grbpkbit,1)
  CALL mpupdater(thisdmp,1)
  CALL mpupdater(tstrtdmp,1)
  CALL mpupdatei(numhdmp,1)
  CALL mpupdater(hdmptim,numhdmp)
  CALL mpupdatei(istager,1)
  CALL mpupdatei(hdfcompr,1)

  splitdmp = hdmpfmt/100
  hdmpfmt  = MOD(hdmpfmt,100)
  IF (hdmpfmt /= 3 .AND. hdmpfmt /= 7) splitdmp = 0

  IF (myproc == 0)THEN
    WRITE(6,'(5x,a,i4)') 'The history dump option was ',hdmpopt

    WRITE(6,'(5x,a,i4)') 'The history dump format was ',hdmpfmt

    !IF (mp_opt > 0) THEN
    !  IF (joindmp == 1) THEN
    !    WRITE(6, '(5x,a)') 'Joined dump (Do not need to run "joinfiles" any more.)'
    !  ELSE
    !    WRITE(6, '(5x,a)') 'Need to run "joinfiles" to get joined files'
    !  END IF
    !END IF
  END IF

  IF( hdmpfmt < 0 .OR. hdmpfmt > 11) THEN
    IF (myproc == 0) WRITE(6,'(5x,a,i4,a,2(/5x,a))')                  &
        'The option hdmpfmt=', hdmpfmt, ' not valid.',                &
        'Program will try to complete reading in input parameters, ', &
        'but will stop at the end of subroutine INITPARA.'
    err_no = err_no + 1
  END IF

!  IF( hdmpfmt.eq.10 .and. nz.ge.256 ) THEN
!    write(6,'(5x,a/5x,a)')
!    :  'The GRIB format can only handle number of vertical levels',
!    :  'less than 256 (8-bits). Reset hdmpfmt to 1 (binary format)'
!    hdmpfmt = 1
!  ENDIF

  IF (myproc == 0)THEN
    WRITE(6,'(5x,a,i6)')                                             &
         'Number of bits in packing GRIB dump data was ',grbpkbit

    WRITE(6,'(5x,a,i4)')                                             &
         'HDF4 compression option was ',hdfcompr

    WRITE(6,'(5x,a,f10.3,a)')                                        &
        'The history dump time interval was ',thisdmp,' seconds '

    WRITE(6,'(5x,a,f10.3,a)')                                        &
        'The history dump starting time was ',tstrtdmp,' seconds '
  END IF

  IF ( hdmpopt == 2 ) THEN

    IF(numhdmp > 0) THEN
      DO i=1,numhdmp
        hdmpstp(i) = nint(hdmptim(i)/dtbig)
      END DO

      nhisdmp = 1
      IF (myproc == 0)THEN
        WRITE(6,'(5x,i3,a,a)') numhdmp,                               &
          ' history data dumps will be produced at',                  &
          ' the following time steps:'
        WRITE(6,'(5x,10i6)') (hdmpstp(i),i=1,numhdmp)
      END IF
    ELSE
      nhisdmp = -1
      IF (myproc == 0)  WRITE(6,'(5x,a)')                             &
          'History data dump is switched off.'
    END IF

  ELSE

    hdmpopt  = 1
    nhisdmp  = nint(thisdmp/dtbig)
    nstrtdmp = nint(tstrtdmp/dtbig)

    thisdmp  = nhisdmp*dtbig
    tstrtdmp = nstrtdmp*dtbig

    IF(nhisdmp > 0) THEN
      IF (myproc == 0) WRITE(6,'(5x,a,i6,a,i6,a/5x,a,a)')             &
          'History data dumps will be produced every ', nhisdmp,      &
          ' time steps after the first ',nstrtdmp,' time steps ',     &
          'where tttttt (if any) stands for the time of the data ',   &
          'in seconds.'
    ELSE
      nhisdmp = -1
      IF (myproc == 0) WRITE(6,'(5x,a)')                              &
          'History data dump is switched off.'
    END IF

  END IF
!
!-----------------------------------------------------------------------
!
!  Read in namelist &output
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) THEN
    READ (unum,output,END=100)
    WRITE(6,'(a)')'Namelist block output sucessfully read.'
  END IF
  CALL mpupdatec(dirname,256)
  CALL mpupdater(tfmtprt,1)
  CALL mpupdatei(exbcdmp,1)
  CALL mpupdatei(exbchdfcompr,1)
  CALL mpupdatei(extdadmp,1)
  CALL mpupdatei(grdout,1)
  CALL mpupdatei(basout,1)
  CALL mpupdatei(varout,1)
  CALL mpupdatei(mstout,1)
  CALL mpupdatei(rainout,1)
  CALL mpupdatei(prcout,1)
  CALL mpupdatei(iceout,1)
  CALL mpupdatei(tkeout,1)
  CALL mpupdatei(trbout,1)
  CALL mpupdatei(sfcout,1)
  CALL mpupdatei(landout,1)
  CALL mpupdatei(totout,1)
  CALL mpupdatei(radout,1)
  CALL mpupdatei(flxout,1)
  CALL mpupdatei(qcexout,1)
  CALL mpupdatei(qrexout,1)
  CALL mpupdatei(qiexout,1)
  CALL mpupdatei(qsexout,1)
  CALL mpupdatei(qhexout,1)
  CALL mpupdatei(qgexout,1)
  CALL mpupdatei(nqexout,1)
  CALL mpupdatei(zqexout,1)
  CALL mpupdater(trstout,1)
  CALL mpupdater(tmaxmin,1)
  CALL mpupdater(tenergy,1)
  CALL mpupdatei(imgopt,1)
  CALL mpupdater(timgdmp,1)
  CALL mpupdatei(pltopt,1)
  CALL mpupdater(tplots,1)
  CALL mpupdatei(filcmprs,1)
  CALL mpupdatei(readyfl,1)
  CALL mpupdatei(sfcdmp,1)
  CALL mpupdatei(soildmp,1)
  CALL mpupdatei(terndmp,1)

  splitexbc = exbcdmp/100
  exbcdmp = MOD(exbcdmp,100)
  IF (exbcdmp /= 3 .AND. exbcdmp /= 7) splitexbc = 0

  splitsoil = soildmp/100
  soildmp   = MOD(soildmp,100)
  IF (soildmp /= 3 .AND. soildmp /= 7) splitsoil = 0

  grdbasfout = MOD ( basout/100+1, 2 )
  basout = MOD(basout,100)
  IF (myproc == 0) WRITE(6,'(5x,a,I10,a)') 'Grid and Base file output is ',grdbasfout,'.'
  IF (myproc == 0) WRITE(6,'(5x,a,f10.3,a)')                       &
      'Formatted printout time interval was ',tfmtprt,' seconds.'

  nfmtprt =  nint(tfmtprt/dtbig)

  IF( nfmtprt /= 0) THEN
    IF (myproc == 0) WRITE(6,'(5x,a,i6,a)')                         &
        'Formatted printing is done every ', nfmtprt,' time steps.'
  ELSE
    nfmtprt = -1
    IF (myproc == 0) WRITE(6,'(5x,a)')                              &
        'Formatted printing is switched off.'
  END IF
!
!-----------------------------------------------------------------------
!
!  Input model output parameters:
!
!  First, give the name of the directory into which output files
!  will be written:
!
!-----------------------------------------------------------------------
!
  ldirnam = LEN_TRIM(dirname)

  IF( ldirnam == 0 ) THEN
    dirname = './'
    ldirnam=2

    IF (myproc == 0) WRITE(6,'(5x,a)')                                  &
        'Output files will be in the current work directory.'
  !ELSE
  !  CALL get_output_dirname(1,dirname,tstart,1,outdirname,istatus)
  END IF

  IF (myproc == 0) THEN
    WRITE(6,'(5x,a,i4)')                                                &
      'The flag to dump out ARPS array into EXBC fields was ',exbcdmp

    WRITE(6,'(5x,a,i4)')                                                &
      'The flag for HDF4 compression was ',exbchdfcompr
  END IF

  IF ( lbcopt /= 2 ) THEN
    extdadmp = 0
    IF (myproc == 0) WRITE(6,'(5x,a,i4)')                            &
     'The flag to dump out EXBC array into ARPS history file was ',  &
     extdadmp
  END IF

  IF (myproc == 0)THEN
    WRITE(6,'(5x,a,i4)')                                             &
      'The flag to dump out ARPS surface data files was ',sfcdmp
    WRITE(6,'(5x,a,i4)')                                             &
      'The flag to dump out ARPS soil data files was ',soildmp
    WRITE(6,'(5x,a,i4)')                                             &
      'The flag to dump out an ARPS terrain data file was ',terndmp
  END IF
!
!-----------------------------------------------------------------------
!
!  Set the control parameters for the output of selected fields.
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0)THEN
    WRITE(6,'(5x,a,i4)')                                             &
         'The input grid coordinate dump option was ', grdout

    WRITE(6,'(5x,a,i4)')                                             &
        'The input base state array dump option was ', basout

    WRITE(6,'(5x,a,i4)')                                             &
        'The input mass-velocity array dump option was ', varout

    WRITE(6,'(5x,a,i4)')                                             &
        'The input non-ice water array dump option was ',mstout

    WRITE(6,'(5x,a,i4)')                                             &
        'The input rain array dump option was ', rainout
  END IF
  rainout = rainout * mstout

  IF (myproc == 0)  WRITE(6,'(5x,a,i4)')                             &
      'The input precipitation rates array dump option was ',prcout
  prcout = prcout * mstout

  IF (ice == 1) THEN
    IF (iceout == 0) THEN
      IF (myproc == 0) WRITE(6,'(1x,a,I4,a,/10x,a)')                 &
      'WARNING: Since option mphyopt = ',mphyopt,                    &
      ', ice array dump should be truned on.',                       &
      'iceout is reset to 1.'
    END IF
    iceout = 1
  END IF

  IF (myproc == 0)THEN
    WRITE(6,'(5x,a,i4)')                                             &
        'The input ice array dump option was ', iceout

    WRITE(6,'(5x,a,i4)')                                             &
        'The input TKE dump option was ', tkeout

    WRITE(6,'(5x,a,i4)')                                             &
        'The input eddy mixing coeff dump option was ', trbout

    WRITE(6,'(5x,a,i4)')                                             &
        'The soil variable dump option was ', sfcout

    WRITE(6,'(5x,a,i4)')                                             &
        'The surface property array dump option was ', landout

    WRITE(6,'(5x,a,i4)')                                             &
        'The radiation arrays dump option was ', radout

    WRITE(6,'(5x,a,i4)')                                             &
        'The surface fluxes dump option was ', flxout

    WRITE(6,'(5x,a,i4)')                                             &
        'The qc EXBC dump option was ', qcexout

    WRITE(6,'(5x,a,i4)')                                             &
        'The qr EXBC dump option was ', qrexout

    WRITE(6,'(5x,a,i4)')                                             &
        'The qi EXBC dump option was ', qiexout

    WRITE(6,'(5x,a,i4)')                                             &
        'The qs EXBC dump option was ', qsexout

    WRITE(6,'(5x,a,i4)')                                             &
        'The qh EXBC dump option was ', qhexout

    WRITE(6,'(5x,a,i4)')                                             &
        'The qg EXBC dump option was ', qgexout

    WRITE(6,'(5x,a,i4)')                                             &
        'Concentration number EXBC dump option was ', nqexout

    WRITE(6,'(5x,a,i4)')                                             &
        'Reflectivity EXBC dump option was ', zqexout
  END IF
!
!-----------------------------------------------------------------------
!
!  Input restart data dump time:
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0)THEN
    WRITE(6,'(5x,a)')                                                 &
       'Specify the time interval between restart data dumps (s):'

    WRITE(6,'(5x,a,f10.3,a)')                                         &
      'Time interval between restart dumps was ',trstout,' seconds '
  END IF

  nrstout =  nint(trstout/dtbig)
  IF( nrstout > 0) THEN
    IF (myproc == 0)  WRITE(6,'(5x,a,/5x,a,i6,a,/5x,a)')              &
        'Restart data files '//runname(1:lfnkey)//'.rsttttttt',       &
        'will be produced every ', nrstout,' time steps',             &
        'where tttttt stands for the time of the data in seconds.'
  ELSE
    nrstout = -1
    IF (myproc == 0) WRITE(6,'(5x,a)')                                &
        'Restart data dump is switched off.'
  END IF
!
!-----------------------------------------------------------------------
!
!  Input parameters for maximum and minimum statistics calculations:
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0)  WRITE(6,'(5x,a,f10.3,a)')                         &
     'Interval between max/min calculations was ',tmaxmin,' seconds '

  nmaxmin = nint(tmaxmin/dtbig)

  IF( tmaxmin > 0.0 .AND. tmaxmin < dtbig ) THEN
    nmaxmin = 1
  ELSE
    nmaxmin = nint(tmaxmin/dtbig)
  END IF

  IF(nmaxmin /= 0)THEN
    IF (myproc == 0)  WRITE(6,'(5x,a,i6,a,/5x,a/)')                   &
        'Max. min. statistics are calculated every ', nmaxmin,        &
        ' time steps.',                                               &
        'and the results are written into file '//runname(1:lfnkey)   &
        //'.maxmin '
  ELSE
    nmaxmin = -1
  IF (myproc == 0)  WRITE(6,'(5x,a)')                                 &
        'Max. min. statistics calculations are switched off.'
  END IF
!
!-----------------------------------------------------------------------
!
!  Input parameter for energy/ptprt variance statistics calculations:
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) WRITE(6,'(5x,a,f10.3)')                            &
      'Interval between energy stats calculations was ',tenergy

  nenergy = nint(tenergy/dtbig)

  IF(nenergy /= 0)THEN
    IF (myproc == 0) WRITE(6,'(5x,a,i6,a,/5x,a)')                     &
    'Energy statistics are calculated every',nenergy,' time steps.',  &
    'and the results are written into file '//runname(1:lfnkey)//'.eng '
  ELSE
    nenergy = -1
    IF (myproc == 0)  WRITE(6,'(5x,a)')                               &
        'Energy statistics calculations are switched off.'
  END IF
!
!-----------------------------------------------------------------------
!
!  More output control parameters, for HDF imgage generation, cell-
!  tracking calls and graphic plotting.
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0)THEN
    WRITE(6,'(5x,a,i4)') 'Image dump option was ',imgopt
    WRITE(6,'(5x,a,f10.3)') 'Image dump interval was ',timgdmp
  END IF

  IF( timgdmp == 0.0) THEN
    imgopt = 0
    nimgdmp = 1
  ELSE IF( timgdmp > 0.0 .AND. timgdmp < dtbig ) THEN
    nimgdmp = 1
  ELSE
    nimgdmp = nint(timgdmp/dtbig)
  END IF

  IF (myproc == 0)THEN
    WRITE(6,'(5x,a,i6,a)')                                          &
       'Data files for images will be written every ', nimgdmp,' time steps.'

    WRITE(6,'(5x,a,i4)') 'Plotting option was ', pltopt

    WRITE(6,'(5x,a,f10.3)') 'Plotting interval was ', tplots
  END IF

  nplots  = nint(tplots /dtbig)
  IF(nplots > 0)THEN
    IF (myproc == 0) WRITE(6,'(5x,a,i6,a)')                          &
        'Plotting routine will be called every',nplots,' time steps.'
  ELSE
    nplots = -1
  END IF

  IF (myproc == 0) WRITE(6,'(5x,a,i6)')                              &
    'File compression option was ',filcmprs
!
!-----------------------------------------------------------------------
!
!  Input debug information print controls:
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) THEN
    READ (unum,debug)
    WRITE(6,'(a)')'Namelist block debug sucessfully read.'
  END IF
  CALL mpupdatei(lvldbg,1)

  IF (myproc == 0) WRITE(6,'(5x,a,i4)') 'The debug printing level was ', lvldbg

  GO TO 102
!
!-----------------------------------------------------------------------
!
!  Print out the input parameters.
!  Write out a log file of model parameters which can be used as
!  the input file to re-run the model.
!
!-----------------------------------------------------------------------
!

  100  CONTINUE
  IF (myproc == 0) CALL wrtcomment('Error reading NAMELIST file. Default values used',1)
  CALL arpsstop('ERROR: please check namelist file.',1)

  102  CONTINUE

  IF (unum /= 5 .AND. myproc == 0) THEN
    CLOSE( unum )
    CALL retunit( unum )
  END IF

  IF(myproc ==0) WRITE(6,'(1x,a)') '======================================='

!-----------------------------------------------------------------------
!
! Hardcode some parameters
! (they has been removed from the namelist)
!
!----------------------------------------------------------------------

  fallvalpha = 0.5
  alfcoef    = 0.5
  IF ( trbvimp == 0 ) alfcoef = 1.0
!
!-----------------------------------------------------------------------
!
!  Adjust max_fopen if readsplit > 0 and/or joindmp > 0
!
!-----------------------------------------------------------------------

  IF (mp_opt > 0 ) THEN

    readstride = max_fopen
    IF (ANY(readsplit > 0) ) readstride = nprocs

    dumpstride = max_fopen
    IF (ANY(joindmp > 0) )   dumpstride = nprocs

    IF (nproc_x_out > nproc_x .OR. nproc_y_out > nproc_y) THEN
      IF (MOD((nx-3),nproc_x_out) /= 0) THEN
        WRITE(6,'(/,1x,2a,I4,a,/,8x,a,I5,a,I4,a,/)') 'ERROR: ',         &
        'wrong size of nproc_x_out = ',nproc_x_out,'.',                 &
        'The grid size is ',nx,                                         &
        ' and it''s physical size is not dividable by ',nproc_x_out,'.'
        CALL arpsstop('Wrong value of nproc_x_out.',1)
      END IF

      IF (MOD((ny-3),nproc_y_out) /= 0) THEN
        WRITE(6,'(/,1x,2a,I4,a,/,8x,a,I5,a,I4,a,/)') 'ERROR: ',         &
        'wrong size of nproc_y_out = ',nproc_y_out,'.',                 &
        'The grid size is ',ny,                                         &
        ' and it''s physical size is not dividable by ',nproc_y_out,'.'
        CALL arpsstop('Wrong value of nproc_y_out.',1)
      END IF

      IF (splitdmp  > 0) joindmp(FINDX_H) = 0    ! We should not join history file
      IF (splitexbc > 0) joindmp(FINDX_B) = 0    ! boundary file
      IF (splitsoil > 0) joindmp(FINDX_S) = 0    ! soil file even
    ELSE
      !WRITE(6,'(/,1x,3a,/)') 'WARNING: ',                               &
      !        'Output multiple patches is only supported for HDF4 outputs', &
      !        'nproc_x_out/nproc_y_out was reset.'
      !IF (joindmp > 0) THEN
      !  nproc_x_out = 1
      !  nproc_y_out = 1
      !ELSE
      !  nproc_x_out = nproc_x
      !  nproc_y_out = nproc_y
      !END IF
    END IF

    !IF (readsplit > 0) THEN
    !  nproc_x_in = 1
    !  nproc_y_in = 1
    !ELSE
    !  nproc_x_in = nproc_x
    !  nproc_y_in = nproc_y
    !END IF

  ELSE
    readsplit(:)  = 0
    joindmp(:)    = 0
    readstride = max_fopen
    dumpstride = max_fopen
    nproc_x_out = 1
    nproc_y_out = 1
    !nproc_x_in  = 1
    !nproc_y_in  = 1
  END IF
!
!-----------------------------------------------------------------------
!
!  Compute derived variables.
!
!-----------------------------------------------------------------------
!
  ebc_global = ebc
  wbc_global = wbc
  nbc_global = nbc
  sbc_global = sbc

  IF (mp_opt > 0) THEN  ! Convert from global to processor specific values.
    nx = (nx - 3)/nproc_x + 3
    ny = (ny - 3)/nproc_y + 3
    IF (myproc == 0) THEN
      WRITE(6,'(5x,a,i5)') "Processor nx =",nx
      WRITE(6,'(5x,a,i5)') "Processor ny =",ny
    END IF

    IF (loc_x /= 1) wbc = 0
    IF (loc_x /= nproc_x) ebc = 0
    IF (loc_y /= 1) sbc = 0
    IF (loc_y /= nproc_y) nbc = 0
  END IF

  IF( initopt == 2 .or.initopt == 4 ) THEN
    restrt = 1
  ELSE
    restrt = 0
  END IF

  CALL julday( year, month, day, jday )         ! Get the Julian day

  nudgstp=1
  temr = ndintvl/dtbig
  nudgstp=MAX(nint(temr),1)
  ndintvl=dtbig*nudgstp
  temr = (ndstop-ndstart)/ndintvl
  ndscale=ndgain/MAX(nint(temr),1)

  IF ( mapproj == 0 ) mpfctopt = 0
  mptrmopt = mptrmopt * mpfctopt
  latitud = ctrlat
  longitud= ctrlon
  IF ( mapproj == 0 ) THEN
    trulat1 = ctrlat
    trulat2 = ctrlat
    trulon  = ctrlon
  END IF

  IF(tintegopt == 1) THEN
    temr = 2.0
  ELSE IF(tintegopt == 2 .or. tintegopt == 3) THEN
    temr = 1.0
  END IF

  dtsml0 = dtsml
  nsmstp = MAX( 1, NINT(temr*dtbig/dtsml) )
  IF (tintegopt == 1) THEN
    dtsml  = temr*dtbig/nsmstp
    IF (dtsml > dtsml0) THEN
      nsmstp = nsmstp + 1
      dtsml  = temr*dtbig/nsmstp
    END IF
  ELSE IF(tintegopt == 2 .or. tintegopt == 3) THEN
    IF(nsmstp /= 1) THEN
      DO WHILE(MOD(nsmstp,2) /= 0)
        nsmstp = nsmstp + 1
        dtsml  = temr*dtbig/nsmstp
      END DO
    END IF
    IF (myproc == 0) THEN
      print*,'Number of small steps per big step must be 1 or even for RK3'
      print*,'dtsml may have been adjusted from the requested value (check below)'
    END IF
  END IF
  IF (myproc == 0) THEN
    WRITE(6,'(1x,2(a,F12.5))') 'Input dtsml = ',dtsml0, ', Adjusted dtsml = ',dtsml
  END IF

  dxinv = 1.0/dx
  dyinv = 1.0/dy
  dzinv = 1.0/dz
  dzsoilinv = 1.0/dzsoil
  dzsoilinv2 = dzsoilinv * dzsoilinv

  xl = (nx-3)*dx
  yl = (ny-3)*dy
  zh = (nz-3)*dz

  IF( ternopt /= 0 .OR. strhopt /= 0 ) THEN
    crdtrns = 1
  ELSE
    crdtrns = 0
  END IF

  cbcmixh = cbcmix * dx*dy

  IF( runmod == 2 ) THEN
    dh = dx
  ELSE IF( runmod == 3 ) THEN
    dh = dy
  ELSE
    dh = SQRT(dx*dy)
  END IF

  cfcmh2 = cfcm2h * dh*dh
  cfcmh4 = cfcm4h * dh**4
  cfcmv2 = cfcm2v * dz*dz
  cfcmv4 = cfcm4v * dz**4

  IF ( divdmp == 1 ) THEN    ! isotropic, cdvdmph=cdvdmpv

    IF ( runmod == 1 ) THEN
      temr = MIN(dx,dy,dzmin)
    ELSE IF( runmod == 2 ) THEN
      temr = MIN(dx,dzmin)
    ELSE IF( runmod == 3 ) THEN
      temr = MIN(dy,dzmin)
    ELSE IF( runmod == 4 ) THEN
      temr = dzmin
    END IF

    cdvdmph = divdmpndh * temr **2 / dtsml
    cdvdmpv = cdvdmph

  ELSE IF ( divdmp == 2 ) THEN

    IF ( runmod == 1 ) THEN
      temr = MIN( SQRT(dx*dy), 5000.0 )
    ELSE IF( runmod == 2 ) THEN
      temr = MIN( dx, 5000.0 )
    ELSE IF( runmod == 3 ) THEN
      temr = MIN( dy, 5000.0 )
    ELSE IF( runmod == 4 ) THEN
      temr = dzmin
    END IF

    cdvdmph = divdmpndh * temr**2 / dtsml
    cdvdmpv = divdmpndv * dzmin **2 / dtsml

  END IF

  IF( err_no /= 0 ) THEN
    IF (myproc == 0) WRITE(6,'(5x,i4,a,/5x,a,/5x,a,a)')              &
        err_no, ' fatal errors found with the input parameters.',    &
        'Please check the ARPS input parameters carefully.',         &
        'The values of parameters you have used can be found',       &
        ' in the log file.'
    CALL arpsstop('arpsstop called from INITPARA with an option',1)
  END IF

  RETURN
END SUBROUTINE initpara
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE PRTLOG                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE prtlog(nx,ny,nz,nzsoil,nunit)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Print a log file compatible in the namelist format
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Adwait Sathye
!  9/15/93
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

  INTEGER :: nx,ny,nz      ! The number of grid points in 3 directions
  INTEGER :: nzsoil        ! The number of grid points in the soil
  INTEGER :: nunit         ! The I/O unit to be used for the log file output
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: logfn, outdirname
  INTEGER :: llogfn
  INTEGER :: logfunt
  INTEGER :: istat
  INTEGER :: lenstr,i,j,ncmnt,ip,k
!
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid parameters
  INCLUDE 'bndry.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'exbc.inc'
  INCLUDE 'nudging.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.
  INCLUDE 'agricst.inc'
  INCLUDE 'nodal.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  Write out a log file in namelist format which can be used as
!  the input file for replicating this run.
!
!  First get a name for the log file:
!
!-----------------------------------------------------------------------
!
  IF( nunit == 6 ) THEN

    logfunt = 6
    WRITE(logfunt,'(///2x,a,i3//)')                                     &
        'PRINT OUT OF MODEL PARAMETERS FOR GRID ',mgrid

  ELSE

    CALL get_output_dirname(0,dirname,tstart,1,outdirname,istat)

    CALL gtlogfn(TRIM(outdirname)//runname(1:lfnkey),mgrid, nestgrd, logfn, llogfn)

    CALL getunit( logfunt )

    OPEN (UNIT=logfunt, FILE=trim(logfn(1:llogfn)),      &
          STATUS='new', IOSTAT=istat)

    IF(istat /= 0) THEN

      WRITE(6,'(/3x,a)')'Error in opening log file ',logfn(1:llogfn)
      WRITE(6,'(3x,a/)')'Job stopped in subroutine INITPARA.'
      CALL arpsstop('arpsstop called from PRTLOG with opening log file' &
                     ,1)

    END IF

  END IF
!
!-----------------------------------------------------------------------
!
!  Construct a namelist input file in the namelog file. default
!  format used by the write function for namelist is to write all the
!  data in a single line, as opposed to a single value per line.
!
!-----------------------------------------------------------------------
!
  WRITE (logfunt, '(1x,a)') '&grid_dims'
  WRITE (logfunt, '(3x,a,i4,a)')    'nx = ', nx,  ','
  WRITE (logfunt, '(3x,a,i4,a)')    'ny = ', ny,  ','
  WRITE (logfunt, '(3x,a,i4,a)')    'nz = ', nz,  ','
  WRITE (logfunt, '(1x,a)') '/'

  WRITE (logfunt, '(1x,a)') '&message_passing'
  WRITE (logfunt, '(3x,a,i4,a)') 'nproc_x   = ', nproc_x,   ','
  WRITE (logfunt, '(3x,a,i4,a)') 'nproc_y   = ', nproc_y,   ','
  WRITE (logfunt, '(3x,a,i4,a)') 'max_fopen = ', max_fopen, ','
  WRITE (logfunt, '(1x,a)') '/'

!  nxc = nx ! Base grid dimensions in ARPS AGR
!  nyc = ny ! Base grid dimensions in ARPS AGR
!  nzc = nz ! Base grid dimensions in ARPS AGR
!
!  WRITE (logfunt, '(1x,a)') '&arpsagr'
!  WRITE (logfunt, '(3x,a,i4,a)')    'levfix = ', levfix,  ','
!  WRITE (logfunt, '(3x,a,i4,a)')    '  intrat = ', intrat,  ','
!  WRITE (logfunt, '(3x,a,i4,a)')    '  intratt= ', intratt,  ','
!  WRITE (logfunt, '(3x,a,i4,a)')    'intrpodr= ', intrpodr,  ','
!  WRITE (logfunt, '(3x,a,i10,a)')   'kcheck  = ', kcheck,  ','
!  WRITE (logfunt, '(3x,a,l7,a)')     'verbose1= ', verbose1,  ','
!  WRITE (logfunt, '(3x,a,l7,a)')     'verbose2= ', verbose2,  ','
!  WRITE (logfunt, '(3x,a,l7,a)')     'verbose3= ', verbose3,  ','
!  WRITE (logfunt, '(3x,a,l7,a)')     'verbose4= ', verbose4,  ','
!  WRITE (logfunt, '(3x,a,l7,a)')     'verbose5= ', verbose5,  ','
!  WRITE (logfunt, '(3x,a,l7,a)')     'verbose6= ', verbose6,  ','
!  WRITE (logfunt, '(3x,a,l7,a)')     'rstart  = ', rstart  ,  ','
!  WRITE (logfunt,'(3x,a,a,a)')      'runold = ''', trim(runold), ''','
!  WRITE (logfunt, '(3x,a,f16.4,a)') 'rstime    = ',rstime, ','
!  WRITE (logfunt, '(3x,a,l7,a)')     'rstdump = ', rstdump ,  ','
!  WRITE (logfunt, '(3x,a,l7,a)')     'grdsrt  = ', grdsrt  ,  ','
!  WRITE (logfunt, '(3x,a,i4,a)')    'nfinelv = ', nfinelv, ','
!
!  DO i=1,nfinelv
!    WRITE (logfunt, '(5x,a,i3.3,a,i3,a)')'ngrdnew(',i,') =',ngrdnew(nfinelv),','
!    DO j=1,ngrdnew(nfinelv)
!      WRITE (logfunt,'(7x,a,i3.3,a,i3.3,a,f10.3,a)')'ixc(',j,',',i,') =',ixc(j,i),','
!      WRITE (logfunt,'(7x,a,i3.3,a,i3.3,a,f10.3,a)')'jyc(',j,',',i,') =',jyc(j,i),','
!      WRITE (logfunt,'(7x,a,i3.3,a,i3.3,a,f10.3,a)')'ixln(',j,',',i,') =',ixln(j,i),','
!      WRITE (logfunt,'(7x,a,i3.3,a,i3.3,a,f10.3,a)')'jyln(',j,',',i,') =',jyln(j,i),','
!      WRITE (logfunt,'(7x,a,i3.3,a,i3.3,a,f10.3,a)')'gangle(',j,',',i,') =',gangle(j,i),','
!    ENDDO
!  ENDDO
!  WRITE (logfunt, '(1x,a)') '/'

  ncmnt = MAX( 1, MIN(9,nocmnt) )

  WRITE (logfunt, '(1x,a)') '&comment_lines'
  WRITE (logfunt, '(3x,a,i4,a)') 'nocmnt    = ', ncmnt , ','

  WRITE (cmnt(ncmnt),'(a,i4,a,i4,a,i4,a,a)')                            &
                    ' nx =',nx,', ny =',ny,', nz =',nz,                 &
                    ' ',arpsversion

  DO i=1,ncmnt
    WRITE (logfunt,'(3x,a,i2.2,a,a,a)')                                 &
         'cmnt(',i,')  = ''', trim(cmnt(i)), ''','
  END DO

  WRITE (logfunt, '(1x,a)') '/'
!
!-----------------------------------------------------------------------
!
!  Write the jobname namelist values into the log file.
!
!-----------------------------------------------------------------------
!
  WRITE (logfunt, '(/1x,a)') '&jobname'
  WRITE (logfunt,'(3x,a,a,a)') 'runname   = ''', trim(runname), ''','
  WRITE (logfunt, '(1x,a)') '/'
!
!-----------------------------------------------------------------------
!
!  Write the model config namelist values into the log file.
!
!-----------------------------------------------------------------------
!
  WRITE (logfunt, '(/1x,a)') '&model_configuration'
  WRITE (logfunt, '(3x,a,i4,a)') 'runmod    = ', runmod, ','
  WRITE (logfunt, '(1x,a)') '/'
!
!-----------------------------------------------------------------------
!
!  Write the initialization namelist values into the namelist logfile
!
!-----------------------------------------------------------------------
!
  WRITE (logfunt, '(/1x,a)')    '&initialization'
  WRITE (logfunt, '(3x,a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a)')   &
                   'initime   = ''', year,'-',month,'-',day,'.',        &
                               hour,':',minute,':',second, ''','
  WRITE (logfunt, '(3x,a,i4,a)')    'timeopt   = ', timeopt,  ','
  WRITE (logfunt, '(3x,a,i4,a)')    'initopt   = ', initopt,  ','
  WRITE (logfunt, '(3x,a,i4,a)')    'inibasopt = ', inibasopt,','
  WRITE (logfunt, '(3x,a,i4,a)')    'viniopt   = ', viniopt,  ','
  WRITE (logfunt, '(3x,a,i4,a)')    'soilinitopt = ',soilinitopt,','
  WRITE (logfunt, '(3x,a,i4,a)')    'tsfcopt = ',tsfcopt,','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'soiltintv   = ',soiltintv,  ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'ubar0     = ', ubar0,    ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'vbar0     = ', vbar0,    ','
  WRITE (logfunt, '(3x,a,i4,a)')    'pt0opt    = ', pt0opt,   ','

  ip = 1
  DO WHILE (ptpert0(ip) /= 0.0)
    WRITE (logfunt, '(3x,a,i2,a,f16.4,a)') 'ptpert0(',ip,')   = ', ptpert0(ip),  ','
    WRITE (logfunt, '(3x,a,i2,a,f16.4,a)') 'pt0radx(',ip,')   = ', pt0radx(ip),  ','
    WRITE (logfunt, '(3x,a,i2,a,f16.4,a)') 'pt0rady(',ip,')   = ', pt0rady(ip),  ','
    WRITE (logfunt, '(3x,a,i2,a,f16.4,a)') 'pt0radz(',ip,')   = ', pt0radz(ip),  ','
    WRITE (logfunt, '(3x,a,i2,a,f16.4,a)') 'pt0ctrx(',ip,')   = ', pt0ctrx(ip),  ','
    WRITE (logfunt, '(3x,a,i2,a,f16.4,a)') 'pt0ctry(',ip,')   = ', pt0ctry(ip),  ','
    WRITE (logfunt, '(3x,a,i2,a,f16.4,a)') 'pt0ctrz(',ip,')   = ', pt0ctrz(ip),  ','
    ip = ip + 1
  END DO

  WRITE (logfunt, '(3x,a,a,a)') 'sndfile   = ''', trim(sndfile), ''','
  WRITE (logfunt,'(3x,a,a,a)')  'rstinf    = ''', trim(rstinf), ''','
  WRITE (logfunt, '(3x,a,i4,a)')'inifmt    = ', inifmt, ','
  WRITE (logfunt, '(3x,a,7(i1),a)')'inisplited= ', (MOD(readsplit(i)+1,2),i=FINDX_NUM,1,-1), ','
  WRITE (logfunt, '(3x,a,a,a)') 'inifile   = ''', trim(inifile), ''','
  WRITE (logfunt, '(3x,a,a,a)') 'inigbf    = ''', trim(inigbf),  ''','
  WRITE (logfunt, '(1x,a)')     '/'
!
!-----------------------------------------------------------------------
!
!  Write the nudging namelist values into the namelist logfile
!
!-----------------------------------------------------------------------
!
  WRITE (logfunt, '(/1x,a)')    '&nudging'
  WRITE (logfunt, '(3x,a,i4,a)')    'nudgopt   = ', nudgopt,  ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'ndstart   = ', ndstart,  ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'ndstop    = ', ndstop,   ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'ndintvl   = ', ndintvl,  ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'ndgain    = ', ndgain,   ','
  WRITE (logfunt, '(3x,a,a,a)')     'incrfnam = ''', trim(incrfnam), ''','
  WRITE (logfunt, '(3x,a,i4,a)')    'incrfmt   = ', incrfmt,  ','
  WRITE (logfunt, '(3x,a,i4,a)')    'nudgu     = ', nudgu,    ','
  WRITE (logfunt, '(3x,a,i4,a)')    'nudgv     = ', nudgv,    ','
  WRITE (logfunt, '(3x,a,i4,a)')    'nudgw     = ', nudgw,    ','
  WRITE (logfunt, '(3x,a,i4,a)')    'nudgp     = ', nudgp,    ','
  WRITE (logfunt, '(3x,a,i4,a)')    'nudgpt    = ', nudgpt,   ','
  WRITE (logfunt, '(3x,a,i4,a)')    'nudgqv    = ', nudgqv,   ','
  WRITE (logfunt, '(3x,a,i4,a)')    'nudgqc    = ', nudgqc,   ','
  WRITE (logfunt, '(3x,a,i4,a)')    'nudgqr    = ', nudgqr,   ','
  WRITE (logfunt, '(3x,a,i4,a)')    'nudgqi    = ', nudgqi,   ','
  WRITE (logfunt, '(3x,a,i4,a)')    'nudgqs    = ', nudgqs,   ','
  WRITE (logfunt, '(3x,a,i4,a)')    'nudgqh    = ', nudgqh,   ','
  WRITE (logfunt, '(1x,a)')      '/'
!
!-----------------------------------------------------------------------
!
!  Write the terrain namelist values into the namelist logfile
!
!-----------------------------------------------------------------------
!
  WRITE (logfunt, '(/1x,a)')      '&terrain'
  WRITE (logfunt, '(3x,a,i4,a)')    'ternopt   = ', ternopt, ','
  WRITE (logfunt, '(3x,a,i4,a)')    'mntopt    = ', mntopt,  ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'hmount    = ', hmount,  ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'mntwidx   = ', mntwidx, ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'mntwidy   = ', mntwidy, ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'mntctrx   = ', mntctrx, ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'mntctry   = ', mntctry, ','
  WRITE (logfunt, '(3x,a,i4,a)')    'ternfmt   = ', ternfmt, ','
  WRITE (logfunt, '(3x,a,a,a)')     'terndta  = ''', trim(terndta), ''','
  WRITE (logfunt, '(1x,a)')       '/'
!
!-----------------------------------------------------------------------
!
!  Write the grid namelist values into the namelist logfile.
!
!-----------------------------------------------------------------------
!
  WRITE (logfunt, '(/1x,a)')      '&grid'
  WRITE (logfunt, '(3x,a,f16.4,a)') 'dx        = ', dx,       ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'dy        = ', dy,       ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'dz        = ', dz,       ','
  WRITE (logfunt, '(3x,a,i4,a)')    'strhopt   = ', strhopt,  ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'dzmin     = ', dzmin,    ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'zrefsfc   = ', zrefsfc,  ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'dlayer1   = ', dlayer1,  ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'dlayer2   = ', dlayer2,  ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'strhtune  = ', strhtune, ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'zflat     = ', zflat,    ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'ctrlat    = ', ctrlat,   ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'ctrlon    = ', ctrlon,   ','
  WRITE (logfunt, '(3x,a,i4,a)')    'crdorgnopt= ', crdorgnopt,  ','
  WRITE (logfunt, '(1x,a)')       '/'
!
!-----------------------------------------------------------------------
!
!  Write the map projection namelist values into the namelist logfile.
!
!-----------------------------------------------------------------------
!
  WRITE (logfunt, '(/1x,a)')      '&projection'
  WRITE (logfunt, '(3x,a,i4,a)')    'mapproj   = ',mapproj, ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'trulat1   = ',trulat1, ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'trulat2   = ',trulat2, ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'trulon    = ',trulon,  ','
  WRITE (logfunt, '(3x,a,e15.5,a)') 'sclfct    = ',sclfct,  ','
  WRITE (logfunt, '(3x,a,i4,a)')    'mpfctopt  = ',mpfctopt,','
  WRITE (logfunt, '(3x,a,i4,a)')    'mptrmopt  = ',mptrmopt,','
  WRITE (logfunt, '(3x,a,i4,a)')    'maptest   = ',maptest, ','
  WRITE (logfunt, '(1x,a)')       '/'
!
!-----------------------------------------------------------------------
!
!  Write the timestep namelist values into the namelist logfile.
!
!-----------------------------------------------------------------------
!
  WRITE (logfunt, '(/1x,a)')      '&timestep'
  WRITE (logfunt, '(3x,a,f16.4,a)') 'dtbig     = ', dtbig, ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'tstart    = ', tstart,','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'tstop     = ', tstop, ','
  WRITE (logfunt, '(1x,a)')       '/'
!
!-----------------------------------------------------------------------
!
!  Write the acoustic_wave namelist values into the namelist logfile.
!
!-----------------------------------------------------------------------
!
  WRITE (logfunt, '(/1x,a)')      '&acoustic_wave'
  WRITE (logfunt, '(3x,a,i4,a)')    'vimplct   = ', vimplct, ','
  WRITE (logfunt, '(3x,a,i4,a)')    'ptsmlstp  = ', ptsmlstp,','
  WRITE (logfunt, '(3x,a,i4,a)')    'csopt     = ', csopt,   ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'csfactr   = ', csfactr, ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'csound    = ', csound,  ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'tacoef    = ', tacoef,  ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'dtsml     = ', dtsml,   ','
  WRITE (logfunt, '(1x,a)')       '/'
!
!-----------------------------------------------------------------------
!
!  Write out equation formulation related parameters
!
!-----------------------------------------------------------------------
!
  WRITE (logfunt, '(/1x,a)')      '&equation_formulation'
  WRITE (logfunt, '(3x,a,i4,a)')    'buoyopt   = ', buoyopt,  ','
  WRITE (logfunt, '(3x,a,i4,a)')    'buoy2nd   = ', buoy2nd,  ','
  WRITE (logfunt, '(3x,a,i4,a)')    'rhofctopt = ', rhofctopt,','
  WRITE (logfunt, '(3x,a,i4,a)')    'bsnesq    = ', bsnesq,   ','
  WRITE (logfunt, '(3x,a,i4,a)')    'peqopt    = ', peqopt,   ','
  WRITE (logfunt, '(1x,a)')       '/'
!
!-----------------------------------------------------------------------
!
!  Write out numerics related parameters
!
!-----------------------------------------------------------------------
!
  WRITE (logfunt, '(/1x,a)')      '&numerics'
  WRITE (logfunt, '(3x,a,i4,a)')    'tintegopt = ', tintegopt, ','
  WRITE (logfunt, '(3x,a,i4,a)')    'madvopt   = ', madvopt, ','
  WRITE (logfunt, '(3x,a,i4,a)')    'sadvopt   = ', sadvopt, ','
  WRITE (logfunt, '(3x,a,i3,a)')    'fctorderopt= ',fctorderopt, ','
  WRITE (logfunt, '(3x,a,i3,a)')    'fctadvptprt= ',fctadvptprt, ','
  WRITE (logfunt, '(1x,a)')       '/'
!
!-----------------------------------------------------------------------
!
!  Write all the boundary conditions in the namelist format into
!  the namelist logfile.
!
!-----------------------------------------------------------------------
!
  WRITE (logfunt, '(/1x,a)')      '&boundary_condition_options'
  WRITE (logfunt, '(3x,a,i4,a)')    'lbcopt    = ', lbcopt, ','
  WRITE (logfunt, '(3x,a,i4,a)')    'wbc       = ', wbc_global,    ','
  WRITE (logfunt, '(3x,a,i4,a)')    'ebc       = ', ebc_global,    ','
  WRITE (logfunt, '(3x,a,i4,a)')    'sbc       = ', sbc_global,    ','
  WRITE (logfunt, '(3x,a,i4,a)')    'nbc       = ', nbc_global,    ','
  WRITE (logfunt, '(3x,a,i4,a)')    'tbc       = ', tbc,    ','
  WRITE (logfunt, '(3x,a,i4,a)')    'fftopt    = ', fftopt, ','
  WRITE (logfunt, '(3x,a,i4,a)')    'bbc       = ', bbc,    ','
  WRITE (logfunt, '(3x,a,i4,a)')    'rbcopt    = ', rbcopt, ','
  WRITE (logfunt, '(3x,a,i4,a)')    'rbc_plcb  = ', rbc_plbc,','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'c_phase   = ', c_phase,','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'rlxlbc    = ', rlxlbc ,','
  WRITE (logfunt, '(3x,a,i4,a)')    'pdetrnd   = ', pdetrnd,','
  WRITE (logfunt, '(1x,a)')       '/'
!
!-----------------------------------------------------------------------
!
!  Write the exbcpara namelist value into the log file.
!
!-----------------------------------------------------------------------
!
  lenstr = 80
  CALL strlnth( exbcname, lenstr)

  WRITE (logfunt, '(/1x,a)')      '&exbcpara'
  WRITE (logfunt, '(3x,a)')                                             &
                      'exbcname  = '''//exbcname(1:lenstr)//''','
  WRITE (logfunt, '(3x,a)')     'tinitebd  = '''//tinitebd//''','
  WRITE (logfunt, '(3x,a,i10,a)')   'tintvebd  = ', tintvebd, ','
  WRITE (logfunt, '(3x,a,i10,a)')   'ngbrz     = ', ngbrz,    ','
  WRITE (logfunt, '(3x,a,e15.5,a)') 'brlxhw    = ', brlxhw,   ','
  WRITE (logfunt, '(3x,a,e15.5,a)') 'cbcdmp    = ', cbcdmp,   ','
  WRITE (logfunt, '(3x,a,i4,a)')    'exbcfmt   = ', exbcfmt,  ','
  WRITE (logfunt, '(1x,a)')       '/'
!
!-----------------------------------------------------------------------
!
!  Write the coriolis force namelist values into the logfile.
!
!-----------------------------------------------------------------------
!
  WRITE (logfunt, '(/1x,a)')   '&coriolis_force'
  WRITE (logfunt, '(3x,a,i4,a)') 'coriopt   = ', coriopt, ','
  WRITE (logfunt, '(3x,a,i4,a)') 'earth_curvature   = ', earth_curvature, ','
  WRITE (logfunt, '(3x,a,i4,a)') 'coriotrm  = ', coriotrm,','
  WRITE (logfunt, '(1x,a)')    '/'
!
!-----------------------------------------------------------------------
!
!  Write the turbulence namelist values into the namelist logfile.
!
!-----------------------------------------------------------------------
!
  WRITE (logfunt, '(/1x,a)')      '&turbulence'
  WRITE (logfunt, '(3x,a,i4,a)')    'tmixopt   = ', tmixopt, ','
  WRITE (logfunt, '(3x,a,i4,a)')    'trbisotp  = ', trbisotp,','
  WRITE (logfunt, '(3x,a,i4,a)')    'tkeopt    = ', tkeopt,  ','
  WRITE (logfunt, '(3x,a,i4,a)')    'trbvimp   = ', trbvimp, ','
  WRITE (logfunt, '(3x,a,i4,a)')    'tmixvert  = ', tmixvert,','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'prantl    = ', prantl,  ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'tmixcst   = ', tmixcst, ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'kmlimit   = ', kmlimit, ','
  WRITE (logfunt, '(1x,a)')       '/'
!
!-----------------------------------------------------------------------
!
!  Write a computational_mixing namelist values into the logfile.
!
!-----------------------------------------------------------------------
!
  WRITE (logfunt, '(/1x,a)')      '&computational_mixing'
  WRITE (logfunt, '(3x,a,i4,a)')    'cmix2nd   = ', cmix2nd,','
  WRITE (logfunt, '(3x,a,e15.5,a)') 'cfcm2h    = ', cfcm2h, ','
  WRITE (logfunt, '(3x,a,e15.5,a)') 'cfcm2v    = ', cfcm2v, ','
  WRITE (logfunt, '(3x,a,i4,a)')    'cmix4th   = ', cmix4th,','
  WRITE (logfunt, '(3x,a,e15.5,a)') 'cfcm4h    = ', cfcm4h, ','
  WRITE (logfunt, '(3x,a,e15.5,a)') 'cfcm4v    = ', cfcm4v, ','
  WRITE (logfunt, '(3x,a,e15.5,a)') 'scmixfctr = ', scmixfctr, ','
  WRITE (logfunt, '(3x,a,i4,a)') 'cmix_opt = ', cmix_opt, ','
  WRITE (logfunt, '(1x,a)')       '/'
!
!-----------------------------------------------------------------------
!
!  Calculate divdmpnd and write the divergence namelist data into
!  the namelist log file.
!
!-----------------------------------------------------------------------
!
  WRITE (logfunt, '(/1x,a)')      '&divergence_damping'
  WRITE (logfunt, '(3x,a,i4,a)')    'divdmp    = ', divdmp,   ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'divdmpndh = ', divdmpndh,','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'divdmpndv = ', divdmpndv,','
  WRITE (logfunt, '(1x,a)')       '/'
!
!-----------------------------------------------------------------------
!
!  Write the rayleigh_damping namelist values into the logfile.
!
!-----------------------------------------------------------------------
!
  WRITE (logfunt, '(/1x,a)')      '&rayleigh_damping'
  WRITE (logfunt, '(3x,a,i4,a)')    'raydmp    = ', raydmp, ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'cfrdmp    = ', cfrdmp, ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'zbrdmp    = ', zbrdmp, ','
  WRITE (logfunt, '(1x,a)')       '/'
!
!-----------------------------------------------------------------------
!
!  Write the asselin_time_filter namelist data into the logfile.
!
!-----------------------------------------------------------------------
!
  WRITE (logfunt, '(/1x,a)')      '&asselin_time_filter'
  WRITE (logfunt, '(3x,a,f16.4,a)') 'flteps    = ', flteps, ','
  WRITE (logfunt, '(1x,a)')       '/'

!
!-----------------------------------------------------------------------
!
!  Write the Concentration namelist data into the logfile.
!
!-----------------------------------------------------------------------
!
  WRITE (logfunt, '(/1x,a)')      '&concentration'
  WRITE (logfunt,'(3x,a,i4,a)') ' ccin   = ',ccin, ','
  WRITE (logfunt,'(3x,a,i4,a)') ' cpoint = ',cpoint, ','
  WRITE (logfunt,'(3x,a,50(i4,a))') '   icc     = ',(icc(k),     ', ', k=1,cpoint)
  WRITE (logfunt,'(3x,a,50(i4,a))') '   jcc     = ',(jcc(k),     ', ', k=1,cpoint)
  WRITE (logfunt,'(3x,a,50(i4,a))') '   kcc     = ',(kcc(k),     ', ', k=1,cpoint)
  WRITE (logfunt,'(3x,a,50(i8,a))') '   ccstart = ',(ccstart(k), ', ', k=1,cpoint)
  WRITE (logfunt,'(3x,a,50(i8,a))') '   ccend   = ',(ccend(k),   ', ', k=1,cpoint)
  WRITE (logfunt,'(3x,a,50(f16.4,a))') ' ccemit = ',(ccemit(k),  ', ', k=1,cpoint)
  WRITE (logfunt, '(1x,a)')       '/'
!
!-----------------------------------------------------------------------
!
!  Write the microphysics namelist values into the log file.
!
!-----------------------------------------------------------------------
!
  WRITE (logfunt, '(/1x,a)')      '&microphysics'
  WRITE (logfunt, '(3x,a,i4,a)')    'moist     = ', moist,   ','
  WRITE (logfunt, '(3x,a,i4,a)')    'mphyopt   = ', mphyopt, ','
  WRITE (logfunt, '(3x,a,i4,a)')    'nmphystp   = ', nmphystp, ','
  WRITE (logfunt, '(3x,a,i4,a)')    'dsdpref    = ', dsdpref , ','
  WRITE (logfunt, '(3x,a,e15.5,a)') '  ntcloud    = ', ntcloud , ','
  WRITE (logfunt, '(3x,a,e15.5,a)') '  n0rain    = ', n0rain , ','
  WRITE (logfunt, '(3x,a,e15.5,a)') '  n0snow    = ', n0snow , ','
  WRITE (logfunt, '(3x,a,e15.5,a)') '  n0grpl    = ', n0grpl , ','
  WRITE (logfunt, '(3x,a,e15.5,a)') '  n0hail    = ', n0hail , ','
  WRITE (logfunt, '(3x,a,e15.5,a)') '  rhoice   = ', rhoice, ','
  WRITE (logfunt, '(3x,a,e15.5,a)') '  rhosnow   = ', rhosnow, ','
  WRITE (logfunt, '(3x,a,e15.5,a)') '  rhogrpl   = ', rhogrpl, ','
  WRITE (logfunt, '(3x,a,e15.5,a)') '  rhohail   = ', rhohail, ','
  WRITE (logfunt, '(3x,a,e15.5,a)') '  alpharain   = ', alpharain, ','
  WRITE (logfunt, '(3x,a,e15.5,a)') '  alphaice   = ', alphaice, ','
  WRITE (logfunt, '(3x,a,e15.5,a)') '  alphasnow   = ', alphasnow, ','
  WRITE (logfunt, '(3x,a,e15.5,a)') '  alphagrpl   = ', alphagrpl, ','
  WRITE (logfunt, '(3x,a,e15.5,a)') '  alphahail   = ', alphahail, ','
  WRITE (logfunt, '(3x,a,i4,a)')    'graupel_ON    = ', graupel_ON , ','
  WRITE (logfunt, '(3x,a,i4,a)')    'hail_ON       = ', hail_ON,  ','
  WRITE (logfunt, '(3x,a,i4,a)')    'mpthermdiag    = ', mpthermdiag , ','
  WRITE (logfunt, '(3x,a,i4,a)')    'cnvctopt  = ', cnvctopt,','
  WRITE (logfunt, '(3x,a,e15.5,a)') 'kffbfct   = ', kffbfct, ','
  WRITE (logfunt, '(3x,a,i4,a)')    'kfsubsattrig = ', kfsubsattrig, ','
  WRITE (logfunt, '(3x,a,e15.5,a)') 'wcldbs    = ', wcldbs,  ','
  WRITE (logfunt, '(3x,a,e15.4,a)') 'confrq    = ', confrq,  ','
  WRITE (logfunt, '(3x,a,e15.4,a)') 'qpfgfrq   = ', qpfgfrq, ','
  WRITE (logfunt, '(3x,a,i4,a)')    'idownd    = ', idownd,  ','
  WRITE (logfunt, '(3x,a,i4,a)')    'impfallopt   = ', impfallopt,  ','
  WRITE (logfunt, '(3x,a,i4,a)')    'fallopt   = ', fallopt,  ','
  WRITE (logfunt, '(3x,a,i4,a)')    'subsatopt = ', subsatopt,  ','
  WRITE (logfunt, '(3x,a,e15.5,a)') 'rhsat     = ', rhsat, ','
  WRITE (logfunt, '(3x,a,e15.5,a)') 'rhsatmin  = ', rhsatmin, ','
  WRITE (logfunt, '(3x,a,e15.5,a)') 'dx_rhsatmin  = ', dx_rhsatmin, ','
  WRITE (logfunt, '(3x,a,e15.5,a)') 'dx_rhsat100  = ', dx_rhsat100, ','
  WRITE (logfunt, '(1x,a)')       '/'
!
!-----------------------------------------------------------------------
!
!  Write the radiation namelist values into the log file.
!
!-----------------------------------------------------------------------
!
  WRITE (logfunt, '(/1x,a)')      '&radiation'
  WRITE (logfunt, '(3x,a,i4,a)')    'radopt    = ', radopt,  ','
  WRITE (logfunt, '(3x,a,i4,a)')    'radstgr   = ', radstgr, ','
  WRITE (logfunt, '(3x,a,i4,a)')    'rlwopt    = ', rlwopt,  ','
  WRITE (logfunt, '(3x,a,i4,a)')    'radshade   = ', radshade,  ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'dtrad     = ', dtrad,   ','
  WRITE (logfunt, '(3x,a,i4,a)')    'raddiag   = ', raddiag, ','
  WRITE (logfunt, '(1x,a)')       '/'
!
!-----------------------------------------------------------------------
!
!  Write the surface physics namelist values into the namelist
!  log file.
!
!-----------------------------------------------------------------------
!
  WRITE (logfunt, '(/1x,a)')      '&surface_physics'
  WRITE (logfunt, '(3x,a,i4,a)')    'sfcphy    = ', sfcphy, ','
  WRITE (logfunt, '(3x,a,i4,a)')    'landwtr   = ', landwtr,','
  WRITE (logfunt, '(3x,a,i4,a)')    'cdhwtropt = ', cdhwtropt,','
  WRITE (logfunt, '(3x,a,e15.5,a)') 'cdmlnd    = ', cdmlnd, ','
  WRITE (logfunt, '(3x,a,e15.5,a)') 'cdmwtr    = ', cdmwtr, ','
  WRITE (logfunt, '(3x,a,e15.5,a)') 'cdhlnd    = ', cdhlnd, ','
  WRITE (logfunt, '(3x,a,e15.5,a)') 'cdhwtr    = ', cdhwtr, ','
  WRITE (logfunt, '(3x,a,e15.5,a)') 'cdqlnd    = ', cdqlnd, ','
  WRITE (logfunt, '(3x,a,e15.5,a)') 'cdqwtr    = ', cdqwtr, ','
  WRITE (logfunt, '(3x,a,i4,a)')    'pbldopt   = ', pbldopt,','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'pbldpth0  = ', pbldpth0,','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'lsclpbl0  = ', lsclpbl0,','
!  write (logfunt, '(3x,a,i4,a)')    'sflxdis   = ', sflxdis,','
  WRITE (logfunt, '(3x,a,i4,a)')    'tqflxdis  = ', tqflxdis,','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'dtqflxdis = ', dtqflxdis,','
  WRITE (logfunt, '(3x,a,i4,a)')    'sfcdiag   = ', sfcdiag,','
  WRITE (logfunt, '(1x,a)')       '/'
!
!-----------------------------------------------------------------------
!
!  Write the surface energy budget model (EBM) namelist values into
!  the namelist log file.
!
!-----------------------------------------------------------------------
!
  WRITE (logfunt, '(/1x,a)')      '&soil_ebm'
  WRITE (logfunt, '(3x,a,i4,a)')    'sfcdat    = ', sfcdat,   ','
  WRITE (logfunt, '(3x,a,i4,a)')    'styp      = ', styp,     ','
  WRITE (logfunt, '(3x,a,i4,a)')    'vtyp      = ', vtyp,     ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'lai0      = ', lai0,     ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'roufns0   = ', roufns0,  ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'veg0      = ', veg0,     ','

  lenstr = 256
  CALL strlnth( sfcdtfl,lenstr )
  WRITE (logfunt, '(3x,a,a,a)')                                         &
                       'sfcdtfl   = ''', sfcdtfl(1:lenstr), ''','
!  lenstr = 4
!  CALL strlnth( sitedir,lenstr )
!  WRITE (logfunt, '(3x,a,a,a)')                                         &
!                       'sitedir   = ''', sitedir(1:lenstr), ''','
!
!  lenstr = 4
!  CALL strlnth( sitefile,lenstr )
!  WRITE (logfunt, '(3x,a,a,a)')                                         &
!                       'sitefile   = ''', sitefile(1:lenstr), ''','

  lenstr = 256
  CALL strlnth( sitemeso,lenstr )
  WRITE (logfunt, '(3x,a,a,a)')                                         &
                       'sitemeso   = ''', sitemeso(1:lenstr), ''','

  lenstr = 256
  CALL strlnth( siteflux,lenstr )
  WRITE (logfunt, '(3x,a,a,a)')                                         &
                       'siteflux   = ''', siteflux(1:lenstr), ''','

  lenstr = 256
  CALL strlnth( siternet,lenstr )
  WRITE (logfunt, '(3x,a,a,a)')                                         &
                       'siternet   = ''', siternet(1:lenstr), ''','

  lenstr = 256
  CALL strlnth( sitesoil,lenstr )
  WRITE (logfunt, '(3x,a,a,a)')                                         &
                       'sitesoil   = ''', sitesoil(1:lenstr), ''','

  lenstr = 256
  CALL strlnth( siteveg,lenstr )
  WRITE (logfunt, '(3x,a,a,a)')                                         &
                       'siteveg   = ''', siteveg(1:lenstr), ''','


  WRITE (logfunt, '(3x,a,i4,a)')    'sfcfmt    = ', sfcfmt,   ','

  WRITE (logfunt, '(3x,a,i4,a)') 'soilmodel_forced = ',soilmodel_forced, ','

  WRITE (logfunt, '(3x,a,i4,a)') 'soilmodel_option = ',soilmodel_option, ','
  WRITE (logfunt, '(3x,a,i4,a)') 'nzsoil       = ', nzsoil,  ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'dzsoil    = ', dzsoil,  ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'zrefsoil    = ', zrefsoil,  ','

  IF ( nzsoil > 0 ) THEN
    DO i=1,nzsoil
      WRITE (logfunt, '(3x,a,i3.3,a,f16.4,a)')     &
                      'tsoilint(',i,') = ', tsoilint(i),','
      WRITE (logfunt, '(3x,a,i3.3,a,f16.4,a)')     &
                      'qsoilint(',i,') = ', qsoilint(i),','
    END DO
  END IF

  WRITE (logfunt, '(3x,a,i4.4,a)')  'soilstrhopt    = ', soilstrhopt,  ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'soildzmin      = ', soildzmin,    ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'soildlayer1    = ', soildlayer1,  ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'soildlayer2    = ', soildlayer2,  ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'soilstrhtune   = ', soilstrhtune, ','


  lenstr = 256
  CALL strlnth( sfcdtfl,lenstr )
  WRITE (logfunt, '(3x,a,a,a)')                                         &
                       'sfcdtfl   = ''', sfcdtfl(1:lenstr), ''','
  WRITE (logfunt, '(3x,a,i4,a)')    'sfcfmt    = ', sfcfmt,   ','

  WRITE (logfunt, '(3x,a,i4,a)')    'soilinit  = ', soilinit, ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'ptslnd0   = ', ptslnd0,  ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'ptswtr0   = ', ptswtr0,  ','
  WRITE (logfunt, '(3x,a,e15.5,a)') 'wetcanp0  = ', wetcanp0, ','
  WRITE (logfunt, '(3x,a,e15.5,a)') 'snowdpth0 = ', snowdpth0,','
  WRITE (logfunt, '(3x,a,e15.5,a)') 'ttprt     = ', ttprt,    ','
  WRITE (logfunt, '(3x,a,e15.5,a)') 'tbprt     = ', tbprt,    ','
  WRITE (logfunt, '(3x,a,e15.5,a)') 'wgrat     = ', wgrat,    ','
  WRITE (logfunt, '(3x,a,e15.5,a)') 'w2rat     = ', w2rat,    ','

  lenstr = 256
  CALL strlnth( soilinfl,lenstr )
  WRITE (logfunt, '(3x,a,a,a)')                                         &
                      'soilinfl  = ''', soilinfl(1:lenstr), ''','
  WRITE (logfunt, '(3x,a,i4,a)')    'soilfmt   = ',soilfmt,   ','

  WRITE (logfunt, '(3x,a,i4,a)')    'tsoil_offset   = ',tsoil_offset,   ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'tsoil_offset_amplitude   = ',tsoil_offset_amplitude,   ','

  WRITE (logfunt, '(3x,a,f16.4,a)') 'dtsfc     = ', dtsfc,    ','
  WRITE (logfunt, '(3x,a,i4,a)')    'prtsoilflx  = ', prtsoilflx, ','

  WRITE (logfunt, '(1x,a)')       '/'
!
!-----------------------------------------------------------------------
!
!  Parameters for automatic grid translation.
!
!-----------------------------------------------------------------------
!
  WRITE (logfunt, '(/1x,a)')      '&grdtrans'
  WRITE (logfunt, '(3x,a,i4,a)')    'cltkopt   =',cltkopt,','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'tceltrk   =',tceltrk,','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'tcrestr   =',tcrestr,','
  WRITE (logfunt, '(3x,a,i4,a)')    'grdtrns   =',grdtrns,','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'umove     =',umove,  ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'vmove     =',vmove,  ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'twindow   =',twindow,','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'chkdpth   =',chkdpth,','
  WRITE (logfunt, '(1x,a)')       '/'
!
!-----------------------------------------------------------------------
!
!  Write the history_dump namelist data into the log file.
!
!-----------------------------------------------------------------------
!
  WRITE (logfunt, '(/1x,a)')      '&history_dump'

  WRITE (logfunt, '(3x,a,i4,a)')    'hdmpopt   = ', hdmpopt, ','
  WRITE (logfunt, '(3x,a,i4,a)')    'hdmpfmt   = ', hdmpfmt, ','
  WRITE (logfunt, '(3x,a,7(i1),a)') 'dmp_out_joined   = ', (joindmp(i),i=FINDX_NUM,1,-1), ','
  WRITE (logfunt, '(3x,a,i4,a)')    'grbpkbit  = ', grbpkbit,','
  WRITE (logfunt, '(3x,a,i4,a)')    'hdfcompr  = ', hdfcompr,','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'thisdmp   = ', thisdmp, ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'tstrtdmp  = ', tstrtdmp,','
  WRITE (logfunt, '(3x,a,i4,a)')    'numhdmp   = ', numhdmp, ','
  IF ( numhdmp > 0 ) THEN
    DO i=1,numhdmp
      WRITE (logfunt, '(3x,a,i3.3,a,f16.4,a)')                          &
                              'hdmptim(',i,') = ', hdmptim(i),','
    END DO
  END IF
  WRITE (logfunt, '(1x,a)')       '/'
!
!-----------------------------------------------------------------------
!
!  Write the output namelist data into the log file.
!
!-----------------------------------------------------------------------
!
  WRITE (logfunt, '(/1x,a)')      '&output'

  WRITE (logfunt, '(3x,a)')                                             &
                 'dirname   = '''//dirname(1:ldirnam)//''','

  WRITE (logfunt, '(3x,a,i4,a)')    'exbcdmp   = ', exbcdmp, ','
  WRITE (logfunt, '(3x,a,i4,a)')    'exbchdfcompr = ', exbchdfcompr, ','
  WRITE (logfunt, '(3x,a,i4,a)')    'extdadmp  = ', extdadmp,','
  WRITE (logfunt, '(3x,a,i4,a)')    'filcmprs  = ', filcmprs,','
  WRITE (logfunt, '(3x,a,i4,a)')    'basout    = ', basout,  ','
  WRITE (logfunt, '(3x,a,i4,a)')    'grdout    = ', grdout,  ','
  WRITE (logfunt, '(3x,a,i4,a)')    'varout    = ', varout,  ','
  WRITE (logfunt, '(3x,a,i4,a)')    'mstout    = ', mstout,  ','
  WRITE (logfunt, '(3x,a,i4,a)')    'iceout    = ', iceout,  ','
  WRITE (logfunt, '(3x,a,i4,a)')    'tkeout    = ', tkeout,  ','
  WRITE (logfunt, '(3x,a,i4,a)')    'trbout    = ', trbout,  ','
  WRITE (logfunt, '(3x,a,i4,a)')    'sfcout    = ', sfcout,  ','
  WRITE (logfunt, '(3x,a,i4,a)')    'rainout   = ', rainout, ','
  WRITE (logfunt, '(3x,a,i4,a)')    'prcout    = ', prcout,  ','
  WRITE (logfunt, '(3x,a,i4,a)')    'landout   = ', landout, ','
  WRITE (logfunt, '(3x,a,i4,a)')    'radout    = ', radout,  ','
  WRITE (logfunt, '(3x,a,i4,a)')    'flxout    = ', flxout,  ','
  WRITE (logfunt, '(3x,a,i4,a)')    'qcexout   = ', qcexout, ','
  WRITE (logfunt, '(3x,a,i4,a)')    'qrexout   = ', qrexout, ','
  WRITE (logfunt, '(3x,a,i4,a)')    'qiexout   = ', qiexout, ','
  WRITE (logfunt, '(3x,a,i4,a)')    'qsexout   = ', qsexout, ','
  WRITE (logfunt, '(3x,a,i4,a)')    'qhexout   = ', qhexout, ','
  WRITE (logfunt, '(3x,a,i4,a)')    'qgexout   = ', qgexout, ','
  WRITE (logfunt, '(3x,a,i4,a)')    'nqexout   = ', nqexout, ','
  WRITE (logfunt, '(3x,a,i4,a)')    'zqexout   = ', zqexout, ','
  WRITE (logfunt, '(3x,a,i4,a)')    'sfcdmp    = ',  sfcdmp, ','
  WRITE (logfunt, '(3x,a,i4,a)')    'soildmp   = ', soildmp, ','
  WRITE (logfunt, '(3x,a,i4,a)')    'terndmp   = ', terndmp, ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'tfmtprt   = ', tfmtprt, ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'trstout   = ', trstout, ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'tmaxmin   = ', tmaxmin, ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'tenergy   = ', tenergy, ','
  WRITE (logfunt, '(3x,a,i4,a)')    'imgopt    = ', imgopt,  ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'timgdmp   = ', timgdmp, ','
  WRITE (logfunt, '(3x,a,i4,a)')    'pltopt    = ', pltopt,  ','
  WRITE (logfunt, '(3x,a,f16.4,a)') 'tplots    = ', tplots,  ','
  WRITE (logfunt, '(1x,a)')       '/'
!
!-----------------------------------------------------------------------
!
!  Write the debug namelist value into the log file.
!
!-----------------------------------------------------------------------
!
  WRITE (logfunt, '(/1x,a)')      '&debug'
  WRITE (logfunt, '(3x,a,i4,a)')    'lvldbg    = ', lvldbg, ','
  WRITE (logfunt, '(1x,a)')       '/'

  IF( logfunt /= 6 ) THEN

    CLOSE(UNIT=logfunt)
    CALL retunit( logfunt )

    WRITE(6,'(/3x,a,a,a/)')                                             &
         'Log file ',logfn(1:llogfn),' was produced for this job.'

  END IF

  RETURN
END SUBROUTINE prtlog
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SETGRD                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE setgrd( nx,ny, x, y )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set up the ARPS model grid.
!
!  The structure of this program is as follows:
!
!  1. Get the map projection information.
!     (call subroutine setmapr)
!
!  2. Get the absolute coordinates of the model grid origin on map
!     grid with the origin at north pole.
!     (call subroutine lltoxy)
!
!  3. Set up the model origin.
!     (call subroutine setorig)
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  1/26/94
!
!  MODIFICATIONS:
!
!  7/15/94
!  Change the model grid reference point from the southwest corner to
!  the center of model domain.
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
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
!  OUTPUT:
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
  REAL :: x   (nx)             ! 2-D model grid points east-west
                               ! direction (model grid units)
  REAL :: y   (ny)             ! 2-D model grid points north-south
                               ! direction (model grid units)
!
!-----------------------------------------------------------------------
!
!  Include files: globcst.inc phycst.inc
!
!  dx       Model grid spacing in the x-direction east-west
!           (meters)
!  dy       Analysis grid spacing in the y-direction north-south
!           (meters)
!
!  ctrlat   Latitude of the center of the model grid (deg. N)
!  ctrlon   Longitude of the center of the model grid (deg. E)
!
!wdt update
!  mapproj  type of map projection used to setup the model grid.
!           mapproj = 1  Polar Stereographic projection
!                   = 2  Lambert Contformal
!                   = 3  Mercator projection
!                   = 4  Lat, Lon Projection
!                   = 5  User defined
!  trulat1  The 1st real true latitude of map projection.
!  trulat2  The 2nd real true latitude of map projection.
!  trulon   Real true longitude of map projection.
!  sclfct   Map scale factor (eg. sclfct=1/1000000)
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid parameters
  INCLUDE 'phycst.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j
  REAL :: alatpro(2)
  REAL :: sclf
  REAL :: dxscl             ! Model x-direction grid spacing
                            ! normalized by the map scale
                            ! dxscl=dx/sclf
  REAL :: dyscl             ! Model y-direction grid spacing
                            ! normalized by the map scale
                            ! dyscl=dy/sclf
  REAL :: ctrx, ctry, swx, swy

  REAL :: xsub0, ysub0

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  alatpro(1) = trulat1
  alatpro(2) = trulat2

  IF( sclfct /= 1.0) THEN
    sclf  = 1.0/sclfct
    dxscl = dx*sclf
    dyscl = dy*sclf
  ELSE
    sclf  = 1.0
    dxscl = dx
    dyscl = dy
  END IF

  xsub0 = dx * (nx-3) * (loc_x-1)
  ysub0 = dy * (ny-3) * (loc_y-1)
!
!-----------------------------------------------------------------------
!
!  Note IMPORTANT!!!!: dx and dy are in meters...and the grid is
!  oriented so that the y-axis line through the true longitude of
!  a map projection runs along a longitude line towards
!  the northpole and the x-axis is perpendicular to the y-axis.
!  Create the x,y grid in grid meters (multiplied by sclf), the
!  origin is the southwest corner of the model physical domain as
!  translated from the center point specified by user
!  (ctrlat/ctrlon).
!
!-----------------------------------------------------------------------
!
  CALL setmapr( mapproj,sclf,alatpro,trulon )
                               ! set up parameters for map projection

!
!-----------------------------------------------------------------------
!
!  Find the absolute coordinate (ctrx,ctry) of point (ctrlat,ctrlon)
!  in the latitude-longitude space.
!
!-----------------------------------------------------------------------
!
  CALL lltoxy( 1,1, ctrlat,ctrlon, ctrx, ctry )
!
!-----------------------------------------------------------------------
!
!  Translate the center point to the first physical point, i.e.,
!  the origin of the model grid.
!
!-----------------------------------------------------------------------
!
!     swx = ctrx - (float(nx-3)/2.) * dxscl
!     swy = ctry - (float(ny-3)/2.) * dyscl
  swx = ctrx - (REAL(nproc_x*(nx-3))/2.) * dxscl
  swy = ctry - (REAL(nproc_y*(ny-3))/2.) * dyscl


  IF(crdorgnopt == 0)THEN
    CALL setorig( 1, swx, swy)   ! set up the model origin to the coord.
  ELSE IF(crdorgnopt == 1)THEN   ! Used by arpsEnKF only so far
    CALL setorig( 1, ctrx, ctry)
  ELSE
    WRITE(6,'(1x,a)') 'crdorgnopt is not set correctly'
    CALL arpsstop('ERROR: wrong option of crdorgnopt.',1)
  END IF

  xgrdorg = 0.0
  ygrdorg = 0.0
!
!-----------------------------------------------------------------------
!
!  Calculate the rest of the model grid points in earth meters*sclf
!
!-----------------------------------------------------------------------
!
  DO i=1,nx

!       x(i) = dxscl * (i-2)
    x(i) = sclf*xsub0 + dxscl * (i-2)

  END DO

  DO j=1,ny

!       y(j) = dyscl * (j-2)
    y(j) = sclf*ysub0 + dyscl * (j-2)

  END DO

  CALL setcornerll(nx,ny,x,y)

  RETURN
END SUBROUTINE setgrd

SUBROUTINE getnscalar(nscalarout)
  IMPLICIT NONE
!-----------------------------------------------------------------------
!
! PURPOSE: get the value of nscalar. The variable is defined in the
!          include file globcst.inc. Fortran procedure can get the value
!          directly by including that file, but C procedure cannot
!          include Fortran include file.
!          So it is provided for using with C procedure, 88d2arps only
!          at present
!
! NOTE: it must be called after the call of initpara or get_dims_from_data
!       because it is these two Fortran procedures that initialize this
!       variable.
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(OUT) :: nscalarout

  INCLUDE 'globcst.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  nscalarout = nscalar

  RETURN
END SUBROUTINE getnscalar
