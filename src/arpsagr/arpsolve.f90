SUBROUTINE arpsolve( mptr )

!*********************************************************************
! call the model solver to step grid mptr forward for one time step.
!
!  Implementation for ARPS 4.0
!  Author: Ming Xue, 10/27/1992
!
!  Updates: E.J.Adlerman
!              August 1995 for Arps4.0.22
!
!           Yuhe Liu
!              Nov. 14 1996 for ARPS 4.2.2
!
!*********************************************************************

  INCLUDE 'nodal.inc'
  INCLUDE 'agrialloc.inc'
  INCLUDE 'agrigrid.inc'
  INCLUDE 'agricpu.inc'

  INCLUDE 'globcst.inc'
  INCLUDE 'bndry.inc'
  INCLUDE 'exbc.inc'
  INCLUDE 'radcst.inc'

  INTEGER :: ncalls
  DATA ncalls /0/
  SAVE ncalls

  INTEGER :: frstep            ! Flag for the initial step

  INTEGER :: nx,ny,nz
  INTEGER :: nxy,nxyz,nxyz0

  INTEGER :: rbufsz            ! buffer size for working
                               ! arrays in radiation package

  INTEGER :: nstyps            ! Number of soil type
  PARAMETER (nstyps=4)

  INTEGER :: exbcbufsz         ! EXBC buffer size

  INTEGER :: ncall1, ncall3, ncall4
  DATA ncall1, ncall3, ncall4/ 0,0,0/
  SAVE ncall1, ncall3, ncall4

!*********************************************************************
!  here we get the pointers to the data and pass the
!  stuff to step to actually take the timestep
!*********************************************************************
  CALL resett

  nx = node(5,mptr)
  ny = node(6,mptr)
  nz = node(14,mptr)

  nxy   = nx*ny
  nxyz  = nxy*nz
  nxyz0 = (nx+1)*(ny+1)*(nz+1)

!*****************************************************************
! retrieve constant values from the constant arrays
!*****************************************************************

  ii = igtint(mptr,1)
  ir = igtrel(mptr,1)

  CALL getcnts(nx,ny,nz, a(ii),nsint, a(ir),nsreal)

  mgrid = mptr

  IF ( radopt == 2 ) THEN
    rbufsz = n2d_radiat*nx*ny+n3d_radiat*nx*ny*nz
  ELSE
    rbufsz = 1
  END IF

!*********************************************************************
! 1-d constant-Arrays
!*********************************************************************

  iifax1 = igtns1(mptr,id_ifax1,13)
  iifax2 = igtns1(mptr,id_ifax2,13)

!*********************************************************************
! 1-d x-Arrays
!*********************************************************************

  ix = igtnx1(mptr,id_x)
  itrigs1 = igtnx1(mptr,id_trigs1)
  iwsave2 = igtnx1(mptr,id_wsave2)

!*********************************************************************
! 1-d y-Arrays
!*********************************************************************

  iy = igtny1(mptr,id_y)
  itrigs2 = igtny1(mptr,id_trigs2)
  iwsave1 = igtny1(mptr,id_wsave1)

!*********************************************************************
! 1-d z-Arrays
!*********************************************************************

  iz = igtnz1(mptr,id_z)

!*********************************************************************
! 2-d xy arrays
!*********************************************************************

  ihterain  = igtnxy(mptr,id_hterain, 1)
  isinlat   = igtnxy(mptr,id_sinlat,  1)
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
  iptsfc    = igtnxy(mptr,id_ptsfc,   1)
  ipbldpth  = igtnxy(mptr,id_pbldpth, 3)
  imapfct   = igtnxy(mptr,id_mapfct,  8)
  irad2d    = igtnxy(mptr,id_rad2d,  10)
  iradsw    = igtnxy(mptr,id_radsw,   1)
  irnflx    = igtnxy(mptr,id_rnflx,   1)
  iprcrate  = igtnxy(mptr,id_prcrate, 4)
  ibcrlx    = igtnxy(mptr,id_bcrlx,   1)
  ivwork1   = igtnxy(mptr,id_vwork1,  2)
  ivwork2   = igtnxy(mptr,id_vwork2,  2)
  iraincv   = igtnxy(mptr,id_raincv,  1)
  inca      = igtnxy(mptr,id_nca,     1)
  iusflx    = igtnxy(mptr,id_usflx,   1)
  ivsflx    = igtnxy(mptr,id_vsflx,   1)
  iptsflx   = igtnxy(mptr,id_ptsflx,  1)
  iqvsflx   = igtnxy(mptr,id_qvsflx,  1)

  itemxy1 = igetsp(nxy)

!*********************************************************************
! 2-d xz arrays
!*********************************************************************\

  iudtnb  = igtnxz(mptr,id_udtnb,1)
  iudtsb  = igtnxz(mptr,id_udtsb,1)
  ivdtnb  = igtnxz(mptr,id_vdtnb,1)
  ivdtsb  = igtnxz(mptr,id_vdtsb,1)
  isdtnb  = igtnxz(mptr,id_sdtnb,1)
  isdtsb  = igtnxz(mptr,id_sdtsb,1)
  ipdtnb  = igtnxz(mptr,id_pdtnb,1)
  ipdtsb  = igtnxz(mptr,id_pdtsb,1)
  iwdtnb  = igtnxz(mptr,id_wdtnb,1)
  iwdtsb  = igtnxz(mptr,id_wdtsb,1)

!*********************************************************************
! 2-d yz arrays
!*********************************************************************

  iudteb  = igtnyz(mptr,id_udteb,1)
  iudtwb  = igtnyz(mptr,id_udtwb,1)
  ivdteb  = igtnyz(mptr,id_vdteb,1)
  ivdtwb  = igtnyz(mptr,id_vdtwb,1)
  isdteb  = igtnyz(mptr,id_sdteb,1)
  isdtwb  = igtnyz(mptr,id_sdtwb,1)
  ipdteb  = igtnyz(mptr,id_pdteb,1)
  ipdtwb  = igtnyz(mptr,id_pdtwb,1)
  iwdteb  = igtnyz(mptr,id_wdteb,1)
  iwdtwb  = igtnyz(mptr,id_wdtwb,1)

!*********************************************************************
! 3-d arrays
!*********************************************************************

  iu     = igtxyz(mptr,id_u,    3)
  iv     = igtxyz(mptr,id_v,    3)
  iw     = igtxyz(mptr,id_w,    3)
  iptprt = igtxyz(mptr,id_ptprt,3)
  ipprt  = igtxyz(mptr,id_pprt, 3)
  iqv    = igtxyz(mptr,id_qv,   3)
  iqc    = igtxyz(mptr,id_qc,   3)
  iqr    = igtxyz(mptr,id_qr,   3)
  iqi    = igtxyz(mptr,id_qi,   3)
  iqs    = igtxyz(mptr,id_qs,   3)
  iqh    = igtxyz(mptr,id_qh,   3)
  itke   = igtxyz(mptr,id_tke,  3)

  iubar    = igtxyz(mptr,id_ubar,   1)
  ivbar    = igtxyz(mptr,id_vbar,   1)
  iwbar    = igtxyz(mptr,id_wbar,   1)
  iptbar   = igtxyz(mptr,id_ptbar,  1)
  ipbar    = igtxyz(mptr,id_pbar,   1)
  iptbari  = igtxyz(mptr,id_ptbari, 1)
  ipbari   = igtxyz(mptr,id_pbari,  1)
  irhostr  = igtxyz(mptr,id_rhostr, 1)
  irhostri = igtxyz(mptr,id_rhostri,1)
  iqvbar   = igtxyz(mptr,id_qvbar,  1)

  izp    = igtxyz(mptr,id_zp,   1)
  ij1    = igtxyz(mptr,id_j1,   1)
  ij2    = igtxyz(mptr,id_j2,   1)
  ij3    = igtxyz(mptr,id_j3,   1)
  ij3inv = igtxyz(mptr,id_j3inv,1)
  iaj3x  = igtxyz(mptr,id_aj3x, 1)
  iaj3y  = igtxyz(mptr,id_aj3y, 1)
  iaj3z  = igtxyz(mptr,id_aj3z, 1)

  iwcont    = igtxyz(mptr,id_wcont,   1)
  ikmh      = igtxyz(mptr,id_kmh,     1)
  ikmv      = igtxyz(mptr,id_kmv,     1)
  irprntl   = igtxyz(mptr,id_rprntl,  1)
  ippi      = igtxyz(mptr,id_ppi,     1)
  icsndsq   = igtxyz(mptr,id_csndsq,  1)
  iptcumsrc = igtxyz(mptr,id_ptcumsrc,1)
  iqcumsrc  = igtxyz(mptr,id_qcumsrc, 5)
  iradfrc   = igtxyz(mptr,id_radfrc,  1)
  iw0avg    = igtxyz(mptr,id_w0avg,   1)
!
!-----------------------------------------------------------------------
!
!  3-d temporary arrays for all grids
!
!-----------------------------------------------------------------------
!
  item1  = igetsp( nxyz )
  item2  = igetsp( nxyz )
  item3  = igetsp( nxyz )
  item4  = igetsp( nxyz )
  item5  = igetsp( nxyz )
  item6  = igetsp( nxyz )
  item7  = igetsp( nxyz )
  item8  = igetsp( nxyz )
  item9  = igetsp( nxyz )
  item10 = igetsp( nxyz )
  item11 = igetsp( nxyz )
  item12 = igetsp( nxyz )
  item13 = igetsp( nxyz )
  item14 = igetsp( nxyz )
  item15 = igetsp( nxyz )
  item16 = igetsp( nxyz )
  item17 = igetsp( nxyz )
  item18 = igetsp( nxyz )
  item19 = igetsp( nxyz )
  item20 = igetsp( nxyz )
  item21 = igetsp( nxyz )
  item22 = igetsp( nxyz )
  item23 = igetsp( nxyz )
  item24 = igetsp( nxyz )
  item25 = igetsp( nxyz )
  item26 = igetsp( nxyz )

  item1_0= igetsp( nxyz0 )
  item2_0= igetsp( nxyz0 )
  item3_0= igetsp( nxyz0 )
!
!-----------------------------------------------------------------------
!
!  Get radiation buffer for working arrays
!
!-----------------------------------------------------------------------
!
  iradbuf  = igetsp( rbufsz )
!
!-----------------------------------------------------------------------
!
!  3-d EXBC arrays if lbcopt=2
!
!-----------------------------------------------------------------------
!
  IF ( lbcopt == 2 ) THEN
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

    ELSE IF ( raydmp == 2 ) THEN
!
!-----------------------------------------------------------------------
!
!  3-d EXBC temporary arrays for fine grids
!
!-----------------------------------------------------------------------
!
      exbcbufsz = 4*nxyz
      iexbcbuf  = igetsp( exbcbufsz )

      CALL updexbc(mptr,iexbcbuf)

    ELSE
      exbcbufsz = 1
      iexbcbuf  = igetsp( exbcbufsz )
    END IF
  ELSE
    exbcbufsz = 1
    iexbcbuf  = igetsp( exbcbufsz )
  END IF

!*********************************************************************
!  now call the solver
!*********************************************************************

  nstep = node(13,mptr)
  nstep = nstep + 1
  node(13,mptr) = nstep

  IF( (nstep == 1) .AND. (restrt /= 1) ) THEN
    frstep=1           ! Indicate that this is the initial step of
                       ! model integration. For this step forward
                       ! time integration scheme will be used.
  ELSE                 ! For non-first step or restart run
    frstep=0           ! When frstep=0, leapfrog scheme is used.
  END IF

!  print*,' To call cordintg for grid ',mptr,' at nstep=',nstep
!  print*,'Frstep = ', frstep

  cpu0 = f_cputime()

  CALL cordintg(mptr,frstep,nx,ny,nz,rbufsz,nstyps,exbcbufsz,           &
       a(iu),a(iv),a(iw),a(iwcont),a(iptprt),a(ipprt),                  &
       a(iqv),a(iqc),a(iqr),a(iqi),a(iqs),a(iqh),                       &
       a(itke),a(ipbldpth),                                             &
       a(iudteb),a(iudtwb),a(iudtnb),a(iudtsb),                         &
       a(ivdteb),a(ivdtwb),a(ivdtnb),a(ivdtsb),                         &
       a(iwdteb),a(iwdtwb),a(iwdtnb),a(iwdtsb),                         &
       a(ipdteb),a(ipdtwb),a(ipdtnb),a(ipdtsb),                         &
       a(isdteb),a(isdtwb),a(isdtnb),a(isdtsb),                         &
       a(iubar),a(ivbar),a(iptbar),a(ipbar),a(iptbari),a(ipbari),       &
       a(irhostr),a(irhostri),a(iqvbar),a(ippi),a(icsndsq),             &
       a(ix),a(iy),a(iz),a(izp),a(imapfct),                             &
       a(ij1),a(ij2),a(ij3),a(iaj3x),a(iaj3y),a(iaj3z),a(ij3inv),       &
       a(itrigs1),a(itrigs2),a(iifax1),a(iifax2),                       &
       a(iwsave1),a(iwsave2),a(ivwork1),a(ivwork2),                     &
       a(isinlat),a(ikmh),a(ikmv),a(irprntl),                           &
       a(isoiltyp),a(istypfrct),                                        &
       a(ivegtyp),a(ilai),a(iroufns),a(iveg),                           &
       a(itsfc),a(itsoil),a(iwetsfc),a(iwetdp),a(iwetcanp),             &
       a(isnowdpth),a(iptsfc),a(iqvsfc),                                &
       a(iptcumsrc),a(iqcumsrc),a(iraing),a(irainc),a(iprcrate),        &
       a(iw0avg),a(inca),a(iraincv),                                    &
       a(iradfrc),a(iradsw),a(irnflx),a(irad2d),a(iradbuf),             &
       a(iexbcbuf),a(ibcrlx),                                           &
       a(iusflx),a(ivsflx),a(iptsflx),a(iqvsflx),a(itemxy1),            &
       a(item1),a(item2),a(item3),a(item4),a(item5),                    &
       a(item6),a(item7),a(item8),a(item9),a(item10),                   &
       a(item11),a(item12),a(item13),a(item14),a(item15),               &
       a(item16),a(item17),a(item18),a(item19),a(item20),               &
       a(item21),a(item22),a(item23),a(item24),a(item25),               &
       a(item26),a(item1_0),a(item2_0),a(item3_0) )

  cpu_main = cpu_main + f_cputime() - cpu0

  curtim = curtim + dtbig
  rnode(20,mptr) = curtim

  WRITE(6,'(1x,a,i6,a,f10.2,a,i2)')                                     &
       'nstep=',nstep,', curtim =',curtim,                              &
       ' after calling arpsolve for grid ',mptr

!*********************************************************************
! put back the values constants into the integer and real arrays
!*********************************************************************

  CALL strcnts(nx,ny,nz, a(ii),nsint, a(ir),nsreal)

!*********************************************************************
!  return 2-d and 3-d arrays to their permanent storage when unpacked.
!*********************************************************************

!*********************************************************************
! 2-d xy arrays
!*********************************************************************

  CALL retnxy(mptr,id_hterain, 1,ihterain, .false.)
  CALL retnxy(mptr,id_sinlat,  1,isinlat,  .false.)
  CALL retnxy(mptr,id_soiltyp, 4,isoiltyp, .false.)
  CALL retnxy(mptr,id_stypfrct,4,istypfrct,.false.)
  CALL retnxy(mptr,id_vegtyp,  1,ivegtyp,  .false.)
  CALL retnxy(mptr,id_lai,     1,ilai,     .false.)
  CALL retnxy(mptr,id_roufns,  1,iroufns,  .false.)
  CALL retnxy(mptr,id_veg,     1,iveg,     .false.)
  CALL retnxy(mptr,id_tsfc,    5,itsfc,    .true.)
  CALL retnxy(mptr,id_tsoil,   5,itsoil,   .true.)
  CALL retnxy(mptr,id_wetsfc,  5,iwetsfc,  .true.)
  CALL retnxy(mptr,id_wetdp,   5,iwetdp,   .true.)
  CALL retnxy(mptr,id_wetcanp, 5,iwetcanp, .true.)
  CALL retnxy(mptr,id_qvsfc,   5,iqvsfc,   .true.)
  CALL retnxy(mptr,id_snowdpth, 1,isnowdpth, .true.)
  CALL retnxy(mptr,id_raing,   1,iraing,   .true.)
  CALL retnxy(mptr,id_rainc,   1,irainc,   .true.)
  CALL retnxy(mptr,id_ptsfc,   1,iptsfc,   .true.)
  CALL retnxy(mptr,id_pbldpth, 3,ipbldpth, .true.)
  CALL retnxy(mptr,id_mapfct,  8,imapfct,  .true.)
  CALL retnxy(mptr,id_rad2d,  10,irad2d,   .true.)
  CALL retnxy(mptr,id_radsw,   1,iradsw,   .true.)
  CALL retnxy(mptr,id_rnflx,   1,irnflx,   .true.)
  CALL retnxy(mptr,id_prcrate, 4,iprcrate, .true.)
  CALL retnxy(mptr,id_bcrlx,   1,ibcrlx,   .true.)
  CALL retnxy(mptr,id_vwork1,  2,ivwork1,  .true.)
  CALL retnxy(mptr,id_vwork2,  2,ivwork2,  .true.)
  CALL retnxy(mptr,id_raincv,  1,iraincv,  .true.)
  CALL retnxy(mptr,id_nca,     1,inca,     .true.)
  CALL retnxy(mptr,id_usflx,   1,iusflx,   .true.)
  CALL retnxy(mptr,id_vsflx,   1,ivsflx,   .true.)
  CALL retnxy(mptr,id_ptsflx,  1,iptsflx,  .true.)
  CALL retnxy(mptr,id_qvsflx,  1,iqvsflx,  .true.)

!*********************************************************************
! 2-d xz arrays
!*********************************************************************

  CALL retnxz(mptr,id_udtnb,1,iudtnb,.false.)
  CALL retnxz(mptr,id_udtsb,1,iudtsb,.false.)
  CALL retnxz(mptr,id_vdtnb,1,ivdtnb,.true.)
  CALL retnxz(mptr,id_vdtsb,1,ivdtsb,.true.)
  CALL retnxz(mptr,id_sdtnb,1,isdtnb,.false.)
  CALL retnxz(mptr,id_sdtsb,1,isdtsb,.false.)
  CALL retnxz(mptr,id_pdtnb,1,ipdtnb,.true.)
  CALL retnxz(mptr,id_pdtsb,1,ipdtsb,.true.)
  CALL retnxz(mptr,id_wdtnb,1,iwdtnb,.false.)
  CALL retnxz(mptr,id_wdtsb,1,iwdtsb,.false.)

!*********************************************************************
! 2-d yz arrays
!*********************************************************************

  CALL retnyz(mptr,id_udteb,1,iudteb,.true.)
  CALL retnyz(mptr,id_udtwb,1,iudtwb,.true.)
  CALL retnyz(mptr,id_vdteb,1,ivdteb,.false.)
  CALL retnyz(mptr,id_vdtwb,1,ivdtwb,.false.)
  CALL retnyz(mptr,id_sdteb,1,isdteb,.false.)
  CALL retnyz(mptr,id_sdtwb,1,isdtwb,.false.)
  CALL retnyz(mptr,id_pdteb,1,ipdteb,.true.)
  CALL retnyz(mptr,id_pdtwb,1,ipdtwb,.true.)
  CALL retnyz(mptr,id_wdteb,1,iwdteb,.false.)
  CALL retnyz(mptr,id_wdtwb,1,iwdtwb,.false.)


!*********************************************************************
! 3-d arrays
!*********************************************************************

  CALL retxyz(mptr,id_u,    3,iu,    .true.)
  CALL retxyz(mptr,id_v,    3,iv,    .true.)
  CALL retxyz(mptr,id_w,    3,iw,    .true.)
  CALL retxyz(mptr,id_ptprt,3,iptprt,.true.)
  CALL retxyz(mptr,id_pprt, 3,ipprt, .true.)
  CALL retxyz(mptr,id_qv,   3,iqv,   .true.)
  CALL retxyz(mptr,id_qc,   3,iqc,   .true.)
  CALL retxyz(mptr,id_qr,   3,iqr,   .true.)
  CALL retxyz(mptr,id_qi,   3,iqi,   .true.)
  CALL retxyz(mptr,id_qs,   3,iqs,   .true.)
  CALL retxyz(mptr,id_qh,   3,iqh,   .true.)
  CALL retxyz(mptr,id_tke,  3,itke,  .true.)

  CALL retxyz(mptr,id_ubar,   1,iubar   ,.false.)
  CALL retxyz(mptr,id_vbar,   1,ivbar   ,.false.)
  CALL retxyz(mptr,id_wbar,   1,iwbar   ,.false.)
  CALL retxyz(mptr,id_ptbar,  1,iptbar  ,.false.)
  CALL retxyz(mptr,id_pbar,   1,ipbar   ,.false.)
  CALL retxyz(mptr,id_ptbari, 1,iptbari ,.false.)
  CALL retxyz(mptr,id_pbari,  1,ipbari  ,.false.)
  CALL retxyz(mptr,id_rhostr, 1,irhostr ,.false.)
  CALL retxyz(mptr,id_rhostri,1,irhostri,.false.)
  CALL retxyz(mptr,id_qvbar,  1,iqvbar  ,.false.)

  CALL retxyz(mptr,id_zp,   1,izp,   .false.)
  CALL retxyz(mptr,id_j1,   1,ij1,   .false.)
  CALL retxyz(mptr,id_j2,   1,ij2,   .false.)
  CALL retxyz(mptr,id_j3,   1,ij3,   .false.)
  CALL retxyz(mptr,id_j3inv,1,ij3inv,.false.)
  CALL retxyz(mptr,id_aj3x, 1,iaj3x, .false.)
  CALL retxyz(mptr,id_aj3y, 1,iaj3y, .false.)
  CALL retxyz(mptr,id_aj3z, 1,iaj3z, .false.)

  CALL retxyz(mptr,id_wcont,   1,iwcont,   .true.)
  CALL retxyz(mptr,id_kmh,     1,ikmh,     .true.)
  CALL retxyz(mptr,id_kmv,     1,ikmv,     .true.)
  CALL retxyz(mptr,id_rprntl,  1,irprntl,  .true.)
  CALL retxyz(mptr,id_ppi,     1,ippi,     .false.)
  CALL retxyz(mptr,id_csndsq,  1,icsndsq,  .false.)
  CALL retxyz(mptr,id_ptcumsrc,1,iptcumsrc,.true.)
  CALL retxyz(mptr,id_qcumsrc, 5,iqcumsrc, .true.)
  CALL retxyz(mptr,id_radfrc,  1,iradfrc,  .true.)
  CALL retxyz(mptr,id_w0avg,   1,iw0avg,   .false.)

  IF ( mptr == 1 .AND. lbcopt == 2 ) THEN
    CALL retexbc(1,1,nexbc3d,iexbcbuf, .true.)
  END IF

!*********************************************************************
! re-set all tem. space
!*********************************************************************

  CALL resett

!*********************************************************************
!  we're done here
!*********************************************************************

  IF(.true.)WRITE(6,'(''ARPSOLVE CALLED FOR GRID '',I3)') mptr
  WRITE(6,'('' GRID '',I3,'' STEPPED FORWARD ONE STEP.'')')mptr
  RETURN
END SUBROUTINE arpsolve


