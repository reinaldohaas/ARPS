SUBROUTINE arpsout( mptr , initdata , nestgrd0 )
!*********************************************************************
! call the model output routine to produce normal data dump
! and printing for grid mptr.
!
!  Implementation for ARPS 4.0
!  Author: Ming Xue, 10/27/1992
!  Updates: E.J. Adlerman
!             August 1995 for Arps 4.0.22
!*********************************************************************
  INTEGER :: mptr   ! pointer to a grid
  INTEGER :: initdata  ! flag indicating if output is for data at initial time

  INCLUDE 'nodal.inc'
  INCLUDE 'agrialloc.inc'
  INCLUDE 'agrigrid.inc'
  INCLUDE 'agricpu.inc'
  INCLUDE 'agricst.inc'
!
  INCLUDE 'phycst.inc'      ! Unchanging physical constants
  INCLUDE 'globcst.inc'     ! Global constants that control model
  INCLUDE 'bndry.inc'       ! Boundary condition control parameters
  INCLUDE 'sfcphycst.inc'   ! Unchanging physical constants
  INCLUDE 'cumucst.inc'     ! Work arrays for cumulus parameterization

  INTEGER :: nx,ny,nz
  INTEGER :: nxy,nxyz

  INTEGER :: nstyps               ! Number of soil type
  PARAMETER (nstyps=4)

  INTEGER :: exbcbufsz            ! EXBC buffer size

  INTEGER :: nxp,nyp,nzp
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  mgrid = mptr
  nestgrd = nestgrd0

  mptrbase = 1

!*********************************************************************
!  here we get the pointers to the data
!*********************************************************************

  CALL resett

!
!-----------------------------------------------------------------------
!
!  Get those parameters from coarse grid which may be updated by
!  input file. The update should be made to all grids.
!
!-----------------------------------------------------------------------
!
  IF ( initdata == 1 ) THEN
    nxp = node(5,mptrbase)
    nyp = node(6,mptrbase)
    nzp = node(14,mptrbase)

    iip = igtint(parent,1)
    irp = igtrel(mptrbase,1)

    CALL getcnts(nxp,nyp,nzp, a(iip),nsint, a(irp),nsreal)

    tstartp = tstart
    tstrtdmpp = tstrtdmp
    thisdmpp = thisdmp
    tfmtprtp = tfmtprt
    tmaxminp = tmaxmin
    trstoutp = trstout
  END IF
!
!-----------------------------------------------------------------------
!
!  Get parameters for this grid.
!
!-----------------------------------------------------------------------
!
  nx = node(5,mptr)
  ny = node(6,mptr)
  nz = node(14,mptr)

  nxy  = nx*ny
  nxyz = nxy*nz

  ii = igtint(mptr,1)
  ir = igtrel(mptr,1)

!*********************************************************************
! retrieve constant values from the constant arrays
!*********************************************************************

  CALL getcnts(nx,ny,nz, a(ii),nsint, a(ir),nsreal)

  IF ( initdata == 1 ) THEN
    tstart = tstartp
    curtim = tstart
    tstrtdmp = tstrtdmpp
    thisdmp = thisdmpp
    tfmtprt = tfmtprtp
    tmaxmin = tmaxminp
    trstout = trstoutp
  END IF

  nfmtprt = nint(tfmtprt/dtbig )
  IF( nfmtprt == 0 ) nfmtprt = -1

  nhisdmp = nint(thisdmp/dtbig)
  IF( nhisdmp == 0 ) nhisdmp = -1

  nstrtdmp = nint(tstrtdmp/dtbig)

  nrstout = nint(trstout/dtbig)
  IF( nrstout == 0 ) nrstout = -1

  nmaxmin = nint(tmaxmin/dtbig)
  IF( nmaxmin == 0 ) nmaxmin = -1

  nimgdmp = nint(timgdmp/dtbig)
  IF( nimgdmp == 0 ) nimgdmp = -1

  nceltrk = nint(tceltrk/dtbig)
  IF( nceltrk == 0 ) nceltrk = -1

  nplots  = nint(tplots /dtbig)
  IF( nplots == 0 ) nplots = -1

  nstep = node(13,mptr)

  WRITE(6,'(1x,a,i2,a,f8.1)')                                           &
      'Calling ARPSOUT for field output for grid ',mptr,                &
      ' at time =', curtim

!*********************************************************************
! 1-d x-Arrays
!*********************************************************************

  ix = igtnx1(mptr,id_x)

!*********************************************************************
! 1-d y-Arrays
!*********************************************************************

  iy = igtny1(mptr,id_y)

!*********************************************************************
! 1-d z-Arrays
!*********************************************************************

  iz = igtnz1(mptr,id_z)

!*********************************************************************
! 2-d xy arrays
!*********************************************************************

  ihterain  = igtnxy(mptr,id_hterain, 1)
  isoiltyp  = igtnxy(mptr,id_soiltyp, 4)
  istypfrct = igtnxy(mptr,id_stypfrct,4)
  ivegtyp   = igtnxy(mptr,id_vegtyp,  1)
  ilai      = igtnxy(mptr,id_lai,     1)
  iroufns   = igtnxy(mptr,id_roufns,  1)
  iveg      = igtnxy(mptr,id_veg,     1)
  itsfc     = igtnxy(mptr,id_tsfc,    5)
  itsoil    = igtnxy(mptr,id_tsoil,   5)
  iwetsfc   = igtnxy(mptr,id_wetsfc,  5)
  iwetdp    = igtnxy(mptr,id_wetdp,   5)
  iwetcanp  = igtnxy(mptr,id_wetcanp, 5)
  iqvsfc    = igtnxy(mptr,id_qvsfc,   5)
  isnowdpth = igtnxy(mptr,id_snowdpth, 1)
  iraing    = igtnxy(mptr,id_raing,   1)
  irainc    = igtnxy(mptr,id_rainc,   1)
  imapfct   = igtnxy(mptr,id_mapfct,  8)
  iradsw    = igtnxy(mptr,id_radsw,   1)
  irnflx    = igtnxy(mptr,id_rnflx,   1)
  iprcrate  = igtnxy(mptr,id_prcrate, 4)
  iraincv   = igtnxy(mptr,id_raincv,  1)
  inca      = igtnxy(mptr,id_nca,     1)
  iusflx    = igtnxy(mptr,id_usflx,   1)
  ivsflx    = igtnxy(mptr,id_vsflx,   1)
  iptsflx   = igtnxy(mptr,id_ptsflx,  1)
  iqvsflx   = igtnxy(mptr,id_qvsflx,  1)

!*********************************************************************
! 2-d xz arrays
!*********************************************************************

  ivdtnb  = igtnxz(mptr,id_vdtnb,1)
  ivdtsb  = igtnxz(mptr,id_vdtsb,1)
  ipdtnb  = igtnxz(mptr,id_pdtnb,1)
  ipdtsb  = igtnxz(mptr,id_pdtsb,1)

!*********************************************************************
! 2-d yz arrays
!*********************************************************************

  iudteb  = igtnyz(mptr,id_udteb,1)
  iudtwb  = igtnyz(mptr,id_udtwb,1)
  ipdteb  = igtnyz(mptr,id_pdteb,1)
  ipdtwb  = igtnyz(mptr,id_pdtwb,1)

!*********************************************************************
! 3-d arrays
!*********************************************************************

  iu     = igtxyz(mptr,id_u,    3)
  iv     = igtxyz(mptr,id_v,    3)
  iw     = igtxyz(mptr,id_w,    3)
  iwcont = igtxyz(mptr,id_wcont,3)
  iptprt = igtxyz(mptr,id_ptprt,3)
  ipprt  = igtxyz(mptr,id_pprt, 3)
  iqv    = igtxyz(mptr,id_qv,   3)
  iqc    = igtxyz(mptr,id_qc,   3)
  iqr    = igtxyz(mptr,id_qr,   3)
  iqi    = igtxyz(mptr,id_qi,   3)
  iqs    = igtxyz(mptr,id_qs,   3)
  iqh    = igtxyz(mptr,id_qh,   3)
  itke   = igtxyz(mptr,id_tke,  3)

  iubar  = igtxyz(mptr,id_ubar,  1)
  ivbar  = igtxyz(mptr,id_vbar,  1)
  iwbar  = igtxyz(mptr,id_wbar,  1)
  iptbar = igtxyz(mptr,id_ptbar, 1)
  ipbar  = igtxyz(mptr,id_pbar,  1)
  irhostr= igtxyz(mptr,id_rhostr,1)
  iqvbar = igtxyz(mptr,id_qvbar, 1)

  izp    = igtxyz(mptr,id_zp,   1)
  ij1    = igtxyz(mptr,id_j1,   1)
  ij2    = igtxyz(mptr,id_j2,   1)
  ij3    = igtxyz(mptr,id_j3,   1)
  ij3inv = igtxyz(mptr,id_j3inv,1)

  ikmh      = igtxyz(mptr,id_kmh,     1)
  ikmv      = igtxyz(mptr,id_kmv,     1)
  iptcumsrc = igtxyz(mptr,id_ptcumsrc,1)
  iqcumsrc  = igtxyz(mptr,id_qcumsrc, 5)
  iradfrc   = igtxyz(mptr,id_radfrc,  1)
  iw0avg    = igtxyz(mptr,id_w0avg,   1)

  IF ( lexbc == 1 ) THEN
    IF ( mptr == 1 ) THEN
!
!-----------------------------------------------------------------------
!
!  3-d EXBC arrays for base grid
!
!-----------------------------------------------------------------------
!
      exbcbufsz = nxyz*nexbc3d
      iexbcbuf  = igtexbc(mptr,1,nexbc3d)
    ELSE

!
!-----------------------------------------------------------------------
!
!  3-d EXBC temporary arrays for fine grids
!
!-----------------------------------------------------------------------
!
      exbcbufsz = 1
      iexbcbuf  = igetsp( exbcbufsz )
    END IF
  ELSE
    exbcbufsz = 1
    iexbcbuf  = igetsp( exbcbufsz )
  END IF
!
!-----------------------------------------------------------------------
!
!  3-d temporary arrays for all grids
!
!-----------------------------------------------------------------------
!
  irhobar = igetsp( nxyz )
  item1   = igetsp( nxyz )
  item2   = igetsp( nxyz )
  item3   = igetsp( nxyz )
  item4   = igetsp( nxyz )
  item5   = igetsp( nxyz )
  item6   = igetsp( nxyz )
  item7   = igetsp( nxyz )
  item8   = igetsp( nxyz )
  item9   = igetsp( nxyz )
  item10  = igetsp( nxyz )
  item11 =  igetsp( nxyz )

  cpu0 = f_cputime()

  IF( initdata == 1) THEN
!
!-----------------------------------------------------------------------
!
!  Reset runname to runnew and store it back to constant arrays
!  before dumping, no matter whether this run is a restart or not.
!
!  This approach should be somewhere for initialization instead of
!  here. But I couldn't find a better place to reset runname in every
!  grid for restart run.
!
!-----------------------------------------------------------------------
!
    IF ( lfnkey /= nmlnthn .OR.  runname(1:lfnkey) /= runnew(1:nmlnthn) ) THEN
      lfnkey = nmlnthn
      runname(1:lfnkey) = runnew(1:nmlnthn)
    END IF
!
!-----------------------------------------------------------------------
!
!  Output the initial fields
!
!-----------------------------------------------------------------------
!
    CALL initout(mptr,nx,ny,nz,nstyps,exbcbufsz,                        &
         a(iu),a(iv),a(iw),a(iwcont),a(iptprt),a(ipprt),                &
         a(iqv),a(iqc),a(iqr),a(iqi),a(iqs),a(iqh),a(itke),             &
         a(iubar),a(ivbar),a(iptbar),a(ipbar),a(irhostr),a(iqvbar),     &
         a(ikmh),a(ikmv),a(ix),a(iy),a(iz),a(izp),a(ihterain),          &
         a(imapfct),a(ij1),a(ij2),a(ij3),a(ij3inv),                     &
         a(isoiltyp),a(istypfrct),                                      &
         a(ivegtyp),a(ilai),a(iroufns),a(iveg),                         &
         a(itsfc),a(itsoil),a(iwetsfc),a(iwetdp),a(iwetcanp),           &
         a(isnowdpth),a(iraing),a(irainc),a(iprcrate),a(iexbcbuf),      &
         a(iradfrc),a(iradsw),a(irnflx),                                &
         a(iusflx),a(ivsflx),a(iptsflx),a(iqvsflx),                     &
         a(item1),a(item2), a(item3),a(item4),a(item5),                 &
         a(item6),a(item7),a(item8),a(item9),a(item10),                 &
         a(item11) )

    IF ( rstart ) THEN
      CALL prtlog(nx,ny,nz,0)
    END IF

    CALL strcnts(nx,ny,nz, a(ii),nsint, a(ir),nsreal)

  ELSE
!
!-----------------------------------------------------------------------
!
!  Output the data at any model time
!
!-----------------------------------------------------------------------
!
    CALL output(mptr,nx,ny,nz,nstyps,exbcbufsz,                         &
         a(iu),a(iv),a(iw),a(iwcont),a(iptprt),a(ipprt),                &
         a(iqv),a(iqc),a(iqr),a(iqi),a(iqs),a(iqh),a(itke),             &
         a(iudteb),a(iudtwb),a(ivdtnb),a(ivdtsb),                       &
         a(ipdteb),a(ipdtwb),a(ipdtnb),a(ipdtsb),                       &
         a(iubar),a(ivbar),a(iptbar),a(ipbar),a(irhostr),a(iqvbar),     &
         a(ikmh),a(ikmv),a(ix),a(iy),a(iz),a(izp),a(ihterain),          &
         a(imapfct), a(ij1),a(ij2),a(ij3),a(ij3inv),                    &
         a(isoiltyp),a(istypfrct),                                      &
         a(ivegtyp),a(ilai),a(iroufns),a(iveg),                         &
         a(itsfc),a(itsoil),a(iwetsfc),a(iwetdp),a(iwetcanp),           &
         a(isnowdpth),a(iqvsfc),a(iptcumsrc),a(iqcumsrc),               &
         a(iw0avg),a(inca),a(iraincv),                                  &
         a(iraing),a(irainc),a(iprcrate),a(iexbcbuf),                   &
         a(iradfrc),a(iradsw),a(irnflx),                                &
         a(iusflx),a(ivsflx),a(iptsflx),a(iqvsflx),                     &
         a(item1),a(item2),a(item3),a(item4),a(item5),                  &
         a(item6),a(item7),a(item8),a(item9),a(item10),                 &
         a(item11) )

  END IF

  cpu_usrout= cpu_usrout + f_cputime() - cpu0

!*********************************************************************
!  return 2-d and 3-d arrays to their permanent storage when packed.
!*********************************************************************

!*********************************************************************
! 2-d xy arrays
!*********************************************************************

  CALL retnxy(mptr,id_hterain, 1,ihterain, .false.)
  CALL retnxy(mptr,id_soiltyp, 4,isoiltyp, .false.)
  CALL retnxy(mptr,id_stypfrct,4,istypfrct,.false.)
  CALL retnxy(mptr,id_vegtyp,  1,ivegtyp,  .false.)
  CALL retnxy(mptr,id_lai,     1,ilai,     .false.)
  CALL retnxy(mptr,id_roufns,  1,iroufns,  .false.)
  CALL retnxy(mptr,id_veg,     1,iveg,     .false.)
  CALL retnxy(mptr,id_tsfc,    5,itsfc,    .false.)
  CALL retnxy(mptr,id_tsoil,   5,itsoil,   .false.)
  CALL retnxy(mptr,id_wetsfc,  5,iwetsfc,  .false.)
  CALL retnxy(mptr,id_wetdp,   5,iwetdp,   .false.)
  CALL retnxy(mptr,id_wetcanp, 5,iwetcanp, .false.)
  CALL retnxy(mptr,id_qvsfc,   5,iqvsfc,   .false.)
  CALL retnxy(mptr,id_snowdpth, 1,isnowdpth, .false.)
  CALL retnxy(mptr,id_raing,   1,iraing,   .false.)
  CALL retnxy(mptr,id_rainc,   1,irainc,   .false.)
  CALL retnxy(mptr,id_mapfct,  8,imapfct,  .false.)
  CALL retnxy(mptr,id_radsw,   1,iradsw,   .false.)
  CALL retnxy(mptr,id_rnflx,   1,irnflx,   .false.)
  CALL retnxy(mptr,id_prcrate, 4,iprcrate, .false.)
  CALL retnxy(mptr,id_raincv,  1,iraincv,  .false.)
  CALL retnxy(mptr,id_nca,     1,inca,     .false.)
  CALL retnxy(mptr,id_usflx,   1,iusflx,   .false.)
  CALL retnxy(mptr,id_vsflx,   1,ivsflx,   .false.)
  CALL retnxy(mptr,id_ptsflx,  1,iptsflx,  .false.)
  CALL retnxy(mptr,id_qvsflx,  1,iqvsflx,  .false.)

!*********************************************************************
! 2-d xz arrays
!*********************************************************************

  CALL retnxz(mptr,id_vdtnb,1,ivdtnb,.false.)
  CALL retnxz(mptr,id_vdtsb,1,ivdtsb,.false.)
  CALL retnxz(mptr,id_pdtnb,1,ipdtnb,.false.)
  CALL retnxz(mptr,id_pdtsb,1,ipdtsb,.false.)

!*********************************************************************
! 2-d yz arrays
!*********************************************************************

  CALL retnyz(mptr,id_udteb,1,iudteb,.false.)
  CALL retnyz(mptr,id_udtwb,1,iudtwb,.false.)
  CALL retnyz(mptr,id_pdteb,1,ipdteb,.false.)
  CALL retnyz(mptr,id_pdtwb,1,ipdtwb,.false.)

!*********************************************************************
! 3-d arrays
!*********************************************************************

  CALL retxyz(mptr,id_u,    3,iu    ,.false.)
  CALL retxyz(mptr,id_v,    3,iv    ,.false.)
  CALL retxyz(mptr,id_w,    3,iw    ,.false.)
  CALL retxyz(mptr,id_ptprt,3,iptprt,.false.)
  CALL retxyz(mptr,id_pprt, 3,ipprt ,.false.)
  CALL retxyz(mptr,id_qv,   3,iqv   ,.false.)
  CALL retxyz(mptr,id_qc,   3,iqc   ,.false.)
  CALL retxyz(mptr,id_qr,   3,iqr   ,.false.)
  CALL retxyz(mptr,id_qi,   3,iqi   ,.false.)
  CALL retxyz(mptr,id_qs,   3,iqs   ,.false.)
  CALL retxyz(mptr,id_qh,   3,iqh   ,.false.)
  CALL retxyz(mptr,id_tke,  3,itke  ,.false.)

  CALL retxyz(mptr,id_ubar,  1,iubar  ,.false.)
  CALL retxyz(mptr,id_vbar,  1,ivbar  ,.false.)
  CALL retxyz(mptr,id_wbar,  1,iwbar  ,.false.)
  CALL retxyz(mptr,id_ptbar, 1,iptbar ,.false.)
  CALL retxyz(mptr,id_pbar,  1,ipbar  ,.false.)
  CALL retxyz(mptr,id_rhostr,1,irhostr,.false.)
  CALL retxyz(mptr,id_qvbar, 1,iqvbar ,.false.)

  CALL retxyz(mptr,id_zp,   1,izp,   .false.)
  CALL retxyz(mptr,id_j1,   1,ij1,   .false.)
  CALL retxyz(mptr,id_j2,   1,ij2,   .false.)
  CALL retxyz(mptr,id_j3,   1,ij3,   .false.)
  CALL retxyz(mptr,id_j3inv,1,ij3inv,.false.)

  CALL retxyz(mptr,id_kmh,     1,ikmh,     .false.)
  CALL retxyz(mptr,id_kmv,     1,ikmv,     .false.)
  CALL retxyz(mptr,id_ptcumsrc,1,iptcumsrc,.false.)
  CALL retxyz(mptr,id_qcumsrc, 5,iqcumsrc, .false.)
  CALL retxyz(mptr,id_radfrc,  1,iradfrc,  .false.)
  CALL retxyz(mptr,id_w0avg,   1,iw0avg,   .false.)

  IF ( lbcopt == 2 .AND. mptr == 1 ) THEN
    CALL retexbc(mptr,1,nexbc3d,iexbcbuf,.true.)
  END IF

!*********************************************************************
! re-set all tem. space
!*********************************************************************

  CALL resett

!*********************************************************************
!  we're done here
!*********************************************************************

  RETURN
END SUBROUTINE arpsout
