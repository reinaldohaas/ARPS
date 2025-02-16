SUBROUTINE arpsinit( mptr , nestgrd0 )
!
!-----------------------------------------------------------------------
!
!  Call the model initialization routine to initialize the model
!  variables for grid mptr.
!
!  Author: Ming Xue, 10/27/1992
!
!  Updates: E.J. Adlerman
!               August 1995 for Arps 4.0.22
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: mptr, nestgrd0

  INCLUDE 'nodal.inc'
  INCLUDE 'agrialloc.inc'
  INCLUDE 'agrigrid.inc'
  INCLUDE 'grddsc.inc'
  INCLUDE 'agricpu.inc'
  INCLUDE 'agricst.inc'

  INCLUDE 'bndry.inc'       ! Boundary condition control parameters
  INCLUDE 'cumucst.inc'     ! cumulus parametization
  INCLUDE 'globcst.inc'     ! Global constants that control model
  INCLUDE 'phycst.inc'      ! Unchanging physical constants
  INCLUDE 'radcst.inc'      ! radiation parametization
  INCLUDE 'sfcphycst.inc'   ! Unchanging physical constants
  INCLUDE 'soilcst.inc'     ! soil-vegetation parametization
  INCLUDE 'grid.inc'

  INTEGER :: nx,ny,nz
  INTEGER :: nxy,nxyz,nxyz0

  INTEGER :: i,j,k,ijk

  INTEGER :: nstyps            ! Number of soil type
  PARAMETER (nstyps=4)

  INTEGER :: exbcbufsz

  INTEGER :: iexbcbuf

  INTEGER :: ii,ir
  INTEGER :: ix,iy,iz
  INTEGER :: iifax1,iifax2,itrigs1,itrigs2
  INTEGER :: ivwork1,ivwork2,iwsave1,iwsave2
  INTEGER :: ihterain,isinlat
  INTEGER :: isoiltyp,istypfrct,ivegtyp,ilai,iroufns,iveg
  INTEGER :: itsfc,itsoil,iwetsfc,iwetdp,iwetcanp,iqvsfc,isnowdpth
  INTEGER :: iptsfc,iraing,irainc,ipbldpth,imapfct
  INTEGER :: irad2d,iradsw,irnflx

  INTEGER :: iprcrate,ibcrlx
  INTEGER :: iraincv,inca
  INTEGER :: iusflx,ivsflx,iptsflx,iqvsflx

  INTEGER :: itemxy1

  INTEGER :: iudtnb,iudtsb,ivdtnb,ivdtsb,isdtnb
  INTEGER :: isdtsb,ipdtnb,ipdtsb,iwdtnb,iwdtsb

  INTEGER :: iudteb,iudtwb,ivdteb,ivdtwb,isdteb
  INTEGER :: isdtwb,ipdteb,ipdtwb,iwdteb,iwdtwb

  INTEGER :: iu,iv,iw
  INTEGER :: iptprt,ipprt,iqv,iqc,iqr,iqi,iqs,iqh,itke

  INTEGER :: iubar, ivbar, iwbar
  INTEGER :: iptbar,ipbar,iptbari,ipbari,ippi,iqvbar
  INTEGER :: irhostr,irhostri
  INTEGER :: izp, ij1,ij2,ij3,ij3inv,iaj3x,iaj3y,iaj3z
  INTEGER :: iwcont,ikmh,ikmv,irprntl,icsndsq
  INTEGER :: iptcumsrc,iqcumsrc,iradfrc,iw0avg

  INTEGER :: item1,item2,item3,item4,item5
  INTEGER :: item6,item7,item8,item9,item10
  INTEGER :: item11,item12,item13,item14,item15
  INTEGER :: item16,item17,item18,item19,item20
  INTEGER :: item21,item22,item23,item24,item25,item26
  INTEGER :: item1_0,item2_0,item3_0

  REAL :: pmax,pmin,qvmax,qvmin,qrmax,qrmin

  REAL :: f_cputime
  INTEGER :: igtint,igtrel,igtns1,igtnx1,igtny1,igtnz1
  INTEGER :: igtnxy,igetsp,igtnxz,igtnyz,igtxyz,igtexbc
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  mgrid = mptr
  nestgrd = nestgrd0

  WRITE(6,'(/1x,a,i3/)')                                                &
      'Calling ARPSINIT to initialized grid ',mptr

!*********************************************************************
!  here we get the pointers to the data
!*********************************************************************

  CALL resett

  nx = node(5,mptr)
  ny = node(6,mptr)
  nz = node(14,mptr)

  nxy   = nx*ny
  nxyz  = nxy*nz
  nxyz0 = (nx+1)*(ny+1)*(nz+1)

  ii = igtint(mptr,1)
  ir = igtrel(mptr,1)

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
  iptsfc    = igtnxy(mptr,id_ptsfc,   1)
  isnowdpth = igtnxy(mptr,id_snowdpth, 1)
  iraing    = igtnxy(mptr,id_raing,   1)
  irainc    = igtnxy(mptr,id_rainc,   1)
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
!*********************************************************************

  iudtnb = igtnxz(mptr,id_udtnb,1)
  iudtsb = igtnxz(mptr,id_udtsb,1)
  ivdtnb = igtnxz(mptr,id_vdtnb,1)
  ivdtsb = igtnxz(mptr,id_vdtsb,1)
  isdtnb = igtnxz(mptr,id_sdtnb,1)
  isdtsb = igtnxz(mptr,id_sdtsb,1)
  ipdtnb = igtnxz(mptr,id_pdtnb,1)
  ipdtsb = igtnxz(mptr,id_pdtsb,1)
  iwdtnb = igtnxz(mptr,id_wdtnb,1)
  iwdtsb = igtnxz(mptr,id_wdtsb,1)

!*********************************************************************
! 2-d yz arrays
!*********************************************************************

  iudteb = igtnyz(mptr,id_udteb,1)
  iudtwb = igtnyz(mptr,id_udtwb,1)
  ivdteb = igtnyz(mptr,id_vdteb,1)
  ivdtwb = igtnyz(mptr,id_vdtwb,1)
  isdteb = igtnyz(mptr,id_sdteb,1)
  isdtwb = igtnyz(mptr,id_sdtwb,1)
  ipdteb = igtnyz(mptr,id_pdteb,1)
  ipdtwb = igtnyz(mptr,id_pdtwb,1)
  iwdteb = igtnyz(mptr,id_wdteb,1)
  iwdtwb = igtnyz(mptr,id_wdtwb,1)

!*********************************************************************
! 3-d arrays
!*********************************************************************

  iu = igtxyz(mptr,id_u,    3)
  iv = igtxyz(mptr,id_v,    3)
  iw = igtxyz(mptr,id_w,    3)

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
  ippi     = igtxyz(mptr,id_ppi,    1)

  izp    = igtxyz(mptr,id_zp,   1)
  ij1    = igtxyz(mptr,id_j1,   1)
  ij2    = igtxyz(mptr,id_j2,   1)
  ij3    = igtxyz(mptr,id_j3,   1)
  ij3inv = igtxyz(mptr,id_j3inv,1)
  iaj3x  = igtxyz(mptr,id_aj3x ,1)
  iaj3y  = igtxyz(mptr,id_aj3y ,1)
  iaj3z  = igtxyz(mptr,id_aj3z ,1)

  iwcont    = igtxyz(mptr,id_wcont,   1)
  ikmh      = igtxyz(mptr,id_kmh,     1)
  ikmv      = igtxyz(mptr,id_kmv,     1)
  irprntl   = igtxyz(mptr,id_rprntl,  1)
  icsndsq   = igtxyz(mptr,id_csndsq,  1)
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
  item11  = igetsp( nxyz )
  item12  = igetsp( nxyz )
  item13  = igetsp( nxyz )
  item14  = igetsp( nxyz )
  item15  = igetsp( nxyz )
  item16  = igetsp( nxyz )
  item17  = igetsp( nxyz )
  item18  = igetsp( nxyz )
  item19  = igetsp( nxyz )
  item20  = igetsp( nxyz )
  item21  = igetsp( nxyz )
  item22  = igetsp( nxyz )
  item23  = igetsp( nxyz )
  item24  = igetsp( nxyz )
  item25  = igetsp( nxyz )
  item26  = igetsp( nxyz )

  item1_0 = igetsp( nxyz0 )
  item2_0 = igetsp( nxyz0 )
  item3_0 = igetsp( nxyz0 )

  cpu0 = f_cputime()

  CALL initial(mptr,nx,ny,nz,nstyps,exbcbufsz,                          &
       a(iu),a(iv),a(iw),a(iwcont),a(iptprt),a(ipprt),                  &
       a(iqv),a(iqc),a(iqr),a(iqi),a(iqs),a(iqh),a(itke),               &
       a(iudteb),a(iudtwb),a(iudtnb),a(iudtsb),                         &
       a(ivdteb),a(ivdtwb),a(ivdtnb),a(ivdtsb),                         &
       a(iwdteb),a(iwdtwb),a(iwdtnb),a(iwdtsb),                         &
       a(ipdteb),a(ipdtwb),a(ipdtnb),a(ipdtsb),                         &
       a(isdteb),a(isdtwb),a(isdtnb),a(isdtsb),                         &
       a(iubar),a(ivbar),a(iptbar),a(ipbar),                            &
       a(iptbari),a(ipbari),a(irhostr),a(irhostri),                     &
       a(iqvbar),a(ippi),a(icsndsq),                                    &
       a(ix),a(iy),a(iz),a(izp),a(ihterain),a(imapfct),                 &
       a(ij1),a(ij2),a(ij3),a(iaj3x),a(iaj3y),a(iaj3z),a(ij3inv),       &
       a(itrigs1),a(itrigs2),a(iifax1),a(iifax2),                       &
       a(iwsave1),a(iwsave2),a(ivwork1),a(ivwork2),                     &
       a(isinlat),a(ikmh),a(ikmv),a(irprntl),                           &
       a(isoiltyp),a(istypfrct),                                        &
       a(ivegtyp),a(ilai),a(iroufns),a(iveg),                           &
       a(itsfc),a(itsoil),a(iwetsfc),a(iwetdp),a(iwetcanp),             &
       a(isnowdpth),a(iptsfc),a(iqvsfc),                                &
       a(iptcumsrc),a(iqcumsrc),a(iw0avg),a(inca),a(iraincv),           &
       a(iraing),a(irainc),a(iprcrate),a(iexbcbuf),a(ibcrlx),           &
       a(iradfrc),a(iradsw),a(irnflx),                                  &
       a(iusflx),a(ivsflx),a(iptsflx),a(iqvsflx),a(itemxy1),            &
       a(item1),a(item2),a(item3),a(item4),a(item5),                    &
       a(item6),a(item7),a(item8),a(item9),a(item10),                   &
       a(item11),a(item12),a(item13),a(item14),a(item15),               &
       a(item16),a(item17),a(item18),a(item19),a(item20),               &
       a(item21),a(item22),a(item23),a(item24),a(item25),a(item26),     &
       a(item1_0),a(item2_0),a(item3_0) )

  WRITE(6,'(/1x,a,i3,a/)')                                              &
      'Grid ',mptr,' was initialized by calling ARPS routine INITIAL.'

  cpu_init0 = cpu_init0 + f_cputime() - cpu0

  hxposs(1) = dx
  hyposs(1) = dy
  possk(1)  = dtbig
  rnode(31,mstart) = dz
!
!-----------------------------------------------------------------------
!
!  The current model time (initial time):
!
!-----------------------------------------------------------------------
!
  curtim = tstart
  node(13,mptr) = nint( curtim/dtbig )
  nstep = node(13,mptr)

  WRITE(6,'(1x,a,f13.3,a)')                                             &
       'The initial model time is at ', curtim,' seconds.'
!
!-----------------------------------------------------------------------
!
!  Initialize the external boundary data array
!
!-----------------------------------------------------------------------
!
  PRINT *, 'In arpsinit'
  PRINT *, 'lbcopt = ', lbcopt

!
!-----------------------------------------------------------------------
!
!  Total number of time steps to be taken:
!  This will be not be by the interface to control the time
!  integration. The time period of integration is instead
!  set in the input file for the interface.
!
!-----------------------------------------------------------------------
!

  nsteps = nint( (tstop-tstart)/dtbig )
!
!-----------------------------------------------------------------------
!
!  nsteps is now unused in the model either. It's calculated
!  here for safty reason only.
!
!-----------------------------------------------------------------------
!
!  return 2-d and 3-d arrays to their permanent storage when unpacked.
!
!-----------------------------------------------------------------------
!
  PRINT*,' calling retnxy'

!*********************************************************************
! 2-d xy arrays
!*********************************************************************

  CALL retnxy(mptr,id_hterain, 1,ihterain, .true.)
  CALL retnxy(mptr,id_sinlat,  1,isinlat,  .true.)
  CALL retnxy(mptr,id_soiltyp, 4,isoiltyp, .true.)
  CALL retnxy(mptr,id_stypfrct,4,istypfrct,.true.)
  CALL retnxy(mptr,id_vegtyp,  1,ivegtyp,  .true.)
  CALL retnxy(mptr,id_lai,     1,ilai,     .true.)
  CALL retnxy(mptr,id_roufns,  1,iroufns,  .true.)
  CALL retnxy(mptr,id_veg,     1,iveg,     .true.)
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

  CALL retnxz(mptr,id_udtnb,1,iudtnb,.true.)
  CALL retnxz(mptr,id_udtsb,1,iudtsb,.true.)
  CALL retnxz(mptr,id_vdtnb,1,ivdtnb,.true.)
  CALL retnxz(mptr,id_vdtsb,1,ivdtsb,.true.)
  CALL retnxz(mptr,id_sdtnb,1,isdtnb,.true.)
  CALL retnxz(mptr,id_sdtsb,1,isdtsb,.true.)
  CALL retnxz(mptr,id_pdtnb,1,ipdtnb,.true.)
  CALL retnxz(mptr,id_pdtsb,1,ipdtsb,.true.)
  CALL retnxz(mptr,id_wdtnb,1,iwdtnb,.true.)
  CALL retnxz(mptr,id_wdtsb,1,iwdtsb,.true.)

!*********************************************************************
! 2-d yz arrays
!*********************************************************************

  CALL retnyz(mptr,id_udteb,1,iudteb,.true.)
  CALL retnyz(mptr,id_udtwb,1,iudtwb,.true.)
  CALL retnyz(mptr,id_vdteb,1,ivdteb,.true.)
  CALL retnyz(mptr,id_vdtwb,1,ivdtwb,.true.)
  CALL retnyz(mptr,id_sdteb,1,isdteb,.true.)
  CALL retnyz(mptr,id_sdtwb,1,isdtwb,.true.)
  CALL retnyz(mptr,id_pdteb,1,ipdteb,.true.)
  CALL retnyz(mptr,id_pdtwb,1,ipdtwb,.true.)
  CALL retnyz(mptr,id_wdteb,1,iwdteb,.true.)
  CALL retnyz(mptr,id_wdtwb,1,iwdtwb,.true.)

!*********************************************************************
! 3-d arrays
!*********************************************************************

  PRINT*,' calling retxyz for u'
  CALL retxyz(mptr,id_u ,3,iu     , .true.)

  PRINT*,' calling retxyz for v'
  CALL retxyz(mptr,id_v ,3,iv     , .true.)

  PRINT*,' calling retxyz for w'
  CALL retxyz(mptr,id_w ,3,iw     , .true.)

  PRINT*,' calling retxyz for ptprt'
  CALL retxyz(mptr,id_ptprt,3,iptprt , .true.)

  PRINT*,' calling retxyz for pprt'
  pmax = 0.0
  pmin = 0.0
  PRINT*,'ipprt =', ipprt
  DO ijk = ipprt, ipprt+nx*ny*nz*3-1
!    print*,'ijk,ip ,nx*ny*nz =', ijk,ipprt,nx*ny*nz, a(ijk)
    pmax = MAX( pmax, a(ijk))
    pmin = MIN( pmax, a(ijk))
  END DO
  PRINT*,' pmin,  pmax = ',  pmin,  pmax, ipprt
  CALL retxyz(mptr,id_pprt,3,ipprt  , .true.)

  PRINT*,' calling retxyz for qv  '
  CALL retxyz(mptr,id_qv,3,iqv    , .true.)
  qvmax = 0.0
  qvmin = 0.0
  DO ijk = iqv, iqv + nx*ny*nz*3 - 1
    qvmax = MAX(qvmax,a(ijk))
    qvmin = MIN(qvmin,a(ijk))
  END DO
  PRINT*, 'qvmax, qvmin, i#', qvmax,qvmin, iqv

  PRINT*,' calling retxyz for qc  '
  CALL retxyz(mptr,id_qc,3,iqc    , .true.)

  qrmax = 0.0
  qrmin = 0.0
  DO ijk = iqr, iqr+nx*ny*nz*3-1
!    print*,'ijk,iqr,nx*ny*nz =', ijk,iqr,nx*ny*nz, a(ijk)
    qrmax = MAX(qrmax, a(ijk))
    qrmin = MIN(qrmax, a(ijk))
  END DO
  PRINT*,'qrmin, qrmax = ', qrmin, qrmax, iqr

  PRINT*,' calling retxyz for qr  '
  CALL retxyz(mptr,id_qr,3,iqr    , .true.)

  PRINT*,' calling retxyz for qi  '
  CALL retxyz(mptr,id_qi,3,iqi    , .true.)

  PRINT*,' calling retxyz for qs  '
  CALL retxyz(mptr,id_qs,3,iqs    , .true.)

  PRINT*,' calling retxyz for qh  '
  CALL retxyz(mptr,id_qh,3,iqh    , .true.)

  PRINT*,' calling retxyz for tke  '
  CALL retxyz(mptr,id_tke,3,itke    , .true.)

  PRINT*,' calling retxyz for *bar'
  CALL retxyz(mptr,id_ubar,   1,iubar,   .true.)
  CALL retxyz(mptr,id_vbar,   1,ivbar,   .true.)
  CALL retxyz(mptr,id_wbar,   1,iwbar,   .true.)
  CALL retxyz(mptr,id_ptbar,  1,iptbar,  .true.)
  CALL retxyz(mptr,id_pbar,   1,ipbar,   .true.)
  CALL retxyz(mptr,id_ptbari, 1,iptbari, .true.)
  CALL retxyz(mptr,id_pbari,  1,ipbari,  .true.)
  CALL retxyz(mptr,id_rhostr, 1,irhostr, .true.)
  CALL retxyz(mptr,id_rhostri,1,irhostri,.true.)
  CALL retxyz(mptr,id_qvbar,  1,iqvbar,  .true.)

  PRINT*,' calling retxyz for zp  '
  CALL retxyz(mptr,id_zp,   1,izp,   .true.)
  CALL retxyz(mptr,id_j1,   1,ij1,   .true.)
  CALL retxyz(mptr,id_j2,   1,ij2,   .true.)
  CALL retxyz(mptr,id_j3,   1,ij3,   .true.)
  CALL retxyz(mptr,id_j3inv,1,ij3inv,.true.)
  CALL retxyz(mptr,id_aj3x, 1,iaj3x, .true.)
  CALL retxyz(mptr,id_aj3y, 1,iaj3y, .true.)
  CALL retxyz(mptr,id_aj3z, 1,iaj3z, .true.)

  PRINT*,' calling retxyz for wcont'
  CALL retxyz(mptr,id_wcont,    1,iwcont,   .true.)
  CALL retxyz(mptr,id_kmh,     1,ikmh,     .true.)
  CALL retxyz(mptr,id_kmv,     1,ikmv,     .true.)
  CALL retxyz(mptr,id_rprntl,  1,irprntl,  .true.)
  CALL retxyz(mptr,id_ppi,     1,ippi,     .true.)
  CALL retxyz(mptr,id_csndsq,  1,icsndsq,  .true.)
  CALL retxyz(mptr,id_ptcumsrc,1,iptcumsrc,.true.)
  CALL retxyz(mptr,id_qcumsrc, 5,iqcumsrc, .true.)
  CALL retxyz(mptr,id_radfrc,  1,iradfrc,  .true.)
  CALL retxyz(mptr,id_w0avg,   1,iw0avg,   .true.)

  IF ( lbcopt == 2 .AND. mptr == 1 ) THEN
    CALL retexbc(mptr,1,nexbc3d,iexbcbuf,.true.)
  END IF

  PRINT*,' calling resett'
!
! re-set all tem. space
!
  CALL resett
!
!-----------------------------------------------------------------------
!
! Store constant variables of the model into constant arrays
!
!-----------------------------------------------------------------------
!
  PRINT*,' calling strcnts'

  CALL strcnts(nx,ny,nz, a(ii),nsint, a(ir),nsreal)

  PRINT*,'ii, nsint, ir, nsreal ', ii, nsint, ir, nsreal
  DO i = izp, izp+1859,169
    PRINT*,'izp values..',i, a(i)
  END DO
  DO i = iu, iu+1859,169
    PRINT*,'iu values..',i, a(i)
  END DO

!
  RETURN
END SUBROUTINE arpsinit

!*********************************************************************
!*********************************************************************
!*********************************************************************


SUBROUTINE initngrd( mptr, mparent )

!*********************************************************************
! To initialize constants and uninitialized fields for a
! new grid mptr.
!
! Author: Ming Xue, 11/15/1992
!  updated: E.J.Adlerman, May 1995, arps4.0.18
!*********************************************************************

  INCLUDE 'nodal.inc'
  INCLUDE 'agrialloc.inc'
  INCLUDE 'agrigrid.inc'
  INCLUDE 'grddsc.inc'
  INCLUDE 'agricst.inc'
!
  INCLUDE 'bndry.inc'       ! Boundary condition control parameters
  INCLUDE 'cumucst.inc'     ! Cumulus paramet.
  INCLUDE 'globcst.inc'     ! Global constants that control model
  INCLUDE 'phycst.inc'      ! Unchanging physical constants
  INCLUDE 'sfcphycst.inc'   ! Unchanging physical constants
  INCLUDE 'grid.inc'

  REAL :: ctrx, ctry
!
!-----------------------------------------------------------------------
!
!  Set lat/lon for each grid
!
!-----------------------------------------------------------------------
!
  IF(verbose6) WRITE(6,'(1x,a,i3)') 'Calling INITNGRD to initialized grid ',mptr

  mgrid = mptr

!  mparent = node(1,mptr)
!*********************************************************************
! it seems that node(1,mptr) has not seen set.
! For now we assume the first grid at one level lower to be the
! parent grid.
!*********************************************************************

  level = node(4,mptr)

  IF( level <= 1  ) THEN
    WRITE(6,'(1x,a,i3,/1x,a)')                                          &
        'The new grid to be initialized is on an improper level ',      &
        level,'Job stopped in INITNGRD.'
    STOP
  END IF

!  mparent = lstart( level-1 )

  IF( mparent == 0 ) THEN
    WRITE(6,'(1x,a,i3,/1x,a)') 'No parent grid found for grid ',        &
         mptr,'Job stopped in INITNGRD.'
    STOP
  END IF

!*********************************************************************
!  Initialize the constants and constant arrays for new grid mptr.
!*********************************************************************

  ii = igtint(mptr,1)
  ir = igtrel(mptr,1)

  iip = igtint(mparent,1)
  irp = igtrel(mparent,1)

  PRINT*,' iip, irp for grid ',mparent,iip,irp

  IF(.true.) THEN

    PRINT*,' integers for grid ', mparent
    DO iii = 1,50
      PRINT*,' iii, a(iip),a(irp)',iii, a(iip+iii-1),a(irp+iii-1)
!      print*,'iii,a(iip+irp)',iip+irp+iii,a(iip+irp+iii-1)
    END DO


  END IF

  CALL nwgrdcst(mptr, mparent,                                          &
       a(ii),a(iip),nsint, a(ir),a(irp),nsreal )


!*********************************************************************
!  Then initialize not yet initialized arrays.
!  Others have be initialized by interpolation from its paraent grid.
!*********************************************************************

  CALL resett

  nx = node(5,mptr)
  ny = node(6,mptr)
  nz = node(14,mptr)

  nxy  = nx*ny
  nxyz = nxy*nz

!*********************************************************************
! 1-d constant arrays
!*********************************************************************

  iifax1 = igtns1(mptr,id_ifax1,13)
  iifax2 = igtns1(mptr,id_ifax2,13)

!*********************************************************************
! 1-d Arrays
!*********************************************************************

  ix = igtnx1(mptr,id_x)

  itrigs1 = igtnx1(mptr,id_trigs1)
  iwsave2 = igtnx1(mptr,id_wsave2)

  iy = igtny1(mptr,id_y)

  itrigs2 = igtny1(mptr,id_trigs2)
  iwsave1 = igtny1(mptr,id_wsave1)

  iz = igtnz1(mptr,id_z)

!*********************************************************************
! 2-d Arrays
!*********************************************************************

  imapfct1 = igtnxy(mptr,id_mapfct,1)
  imapfct2 = igtnxy(mptr,id_mapfct+1,1)
  imapfct3 = igtnxy(mptr,id_mapfct+2,1)
  imapfct4 = igtnxy(mptr,id_mapfct+3,1)
  imapfct5 = igtnxy(mptr,id_mapfct+4,1)
  imapfct6 = igtnxy(mptr,id_mapfct+5,1)
  imapfct7 = igtnxy(mptr,id_mapfct+6,1)
  imapfct8 = igtnxy(mptr,id_mapfct+7,1)

  ivwork1  = igtnxy(mptr,id_vwork1,2)
  ivwork2  = igtnxy(mptr,id_vwork2,2)

!*********************************************************************
! 3-d arrays
!*********************************************************************

  izp    = igtxyz(mptr,id_zp,   1)
  ij1    = igtxyz(mptr,id_j1,   1)
  ij2    = igtxyz(mptr,id_j2,   1)
  ij3    = igtxyz(mptr,id_j3,   1)
  ij3inv = igtxyz(mptr,id_j3inv,1)
  iaj3x  = igtxyz(mptr,id_aj3x ,1)
  iaj3y  = igtxyz(mptr,id_aj3y ,1)
  iaj3z  = igtxyz(mptr,id_aj3z ,1)

  iwcont = igtxyz(mptr,id_wcont, 1)
  ikmh   = igtxyz(mptr,id_kmh,   1)
  ikmv   = igtxyz(mptr,id_kmv,   1)
  irprntl= igtxyz(mptr,id_rprntl,1)
!
!-----------------------------------------------------------------------
!
!  Initialize the computational coordinate arrays x, y and z.
!
!-----------------------------------------------------------------------
!
  DO i=1,nx
    a(ix+i-1) = xorig + (i-2) * dx
  END DO

  DO j=1,ny
    a(iy+j-1) = yorig + (j-2) * dy
  END DO

  DO k=1,nz
    a(iz+k-1) = zorig + (k-2) * dz
  END DO

!
!-----------------------------------------------------------------------
!
!  Set lat/lon for each grid
!
!-----------------------------------------------------------------------
!
  ctrx = 0.5*(a(ix)+a(ix+nx-1))
  ctry = 0.5*(a(iy)+a(iy+ny-1))
  CALL xytoll(1,1, ctrx,ctry, ctrlat,ctrlon)
  CALL setcornerll( nx,ny, a(ix),a(iy) )       ! set corner lat/lon

  WRITE(6,'(a,i3,a,2f15.2)')                                            &
      'The x and y at the center of grid ',mptr,' is ',ctrx,ctry
  WRITE(6,'(a,i3,a,2f15.2)')                                            &
      'The lat/lon at the center of grid ',mptr,' is ',ctrlat,ctrlon
!
!-----------------------------------------------------------------------
!
!  Calculate transformation Jacobians J1, J2 and J3.
!  zp and hterrain are interpolated from the parent grid.
!
!-----------------------------------------------------------------------
!
  CALL jacob(nx,ny,nz,a(ix),a(iy),a(iz),a(izp),a(ij1),a(ij2),a(ij3))

  DO k=1,nz
    DO j=1,ny-1
      DO i=1,nx-1
        ijk = nxy*(k-1) + nx*(j-1) + i
        a(ij3inv+ijk-1)=1.0/a(ij3+ijk-1)
      END DO
    END DO
  END DO

  DO k=1,nz
    DO j=1,ny-1
      DO i=2,nx-1
        ijk  = nxy*(k-1) + nx*(j-1) + i
        ijk1 = nxy*(k-1) + nx*(j-1) + i-1
        a(iaj3x+ijk-1)=0.5*(a(ij3+ijk-1)+a(ij3+ijk1-1))
      END DO
    END DO
  END DO

  CALL bcsu(nx,ny,nz,1,ny-1,1,nz,ebc,wbc,a(iaj3x))

  DO k=1,nz
    DO j=2,ny-1
      DO i=1,nx-1
        ijk  = nxy*(k-1) + nx*(j-1) + i
        ijk1 = nxy*(k-1) + nx*(j-2) + i
        a(iaj3y+ijk-1)=0.5*(a(ij3+ijk-1)+a(ij3+ijk1-1))
      END DO
    END DO
  END DO

  CALL bcsv(nx,ny,nz,1,nx-1,1,nz,nbc,sbc,a(iaj3y))

  DO k=2,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        ijk  = nxy*(k-1) + nx*(j-1) + i
        ijk1 = nxy*(k-2) + nx*(j-1) + i
        a(iaj3z+ijk-1)=0.5*(a(ij3+ijk-1)+a(ij3+ijk1-1))
      END DO
    END DO
  END DO

  CALL bcsw(nx,ny,nz,1,nx-1,1,ny-1,tbc,bbc,a(iaj3z))
!
!-----------------------------------------------------------------------
!
!  Set map factors for this grid
!
!-----------------------------------------------------------------------
!
  IF ( mpfctopt /= 0 ) THEN

    DO j=1,ny-1
      DO i=1,nx-1
        xs = 0.5*(a(ix+i-1)+a(ix+i))
        ys = 0.5*(a(iy+j-1)+a(iy+j))
        ij = nx*(j-1) + i
        CALL xytomf( 1,1,xs,ys,a(imapfct1+ij-1) )
        a(imapfct4+ij-1) = 1.0/a(imapfct1+ij-1)
        a(imapfct7+ij-1) = a(imapfct1+ij-1)*a(imapfct1+ij-1)
        a(imapfct8+ij-1) = 0.25*a(imapfct1+ij-1)  ! for use in sovlwpim
                                                  ! and wcontra...
      END DO
    END DO

    DO j=1,ny-1
      DO i=1,nx
        ys = 0.5*(a(iy+j-1)+a(iy+j))
        ij = nx*(j-1) + i
        CALL xytomf( 1,1,a(ix+i-1),ys,a(imapfct2+ij-1) )
        a(imapfct5+ij-1) = 1.0/a(imapfct2+ij-1)
      END DO
    END DO

    DO j=1,ny
      DO i=1,nx-1
        xs = 0.5*(a(ix+i-1)+a(ix+i))
        ij = nx*(j-1) + i
        CALL xytomf( 1,1,xs,a(iy+j-1),a(imapfct3+ij-1) )
        a(imapfct6+ij-1) = 1.0/a(imapfct3+ij-1)
      END DO
    END DO

  ELSE

    DO j=1,ny
      DO i=1,nx
        ij = nx*(j-1) + i
        a(imapfct1+ij-1) = 1.0
        a(imapfct2+ij-1) = 1.0
        a(imapfct3+ij-1) = 1.0
        a(imapfct4+ij-1) = 1.0
        a(imapfct5+ij-1) = 1.0
        a(imapfct6+ij-1) = 1.0
        a(imapfct7+ij-1) = 1.0
        a(imapfct8+ij-1) = 0.25
      END DO
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Set the arrays for top boundary option 4
!
!-----------------------------------------------------------------------
!
  DO i=1,13
    a(iifax1+i-1) = 0
    a(iifax2+i-1) = 0
  END DO

  DO i=1,3*(nx-1)/2+1
    a(itrigs1+i-1) = 0
  END DO

  DO j=1,3*(ny-1)/2+1
    a(itrigs2+j-1) = 0
  END DO

  DO j=1,ny+1
    DO i=1,nx+1
      ij = (nx+1)*(j-1) + i
      a(ivwork1+ij-1) = 0.0
    END DO
  END DO

  DO j=1,ny
    DO i=1,nx+1
      ij = (nx+1)*(j-1) + i
      a(ivwork2+ij-1) = 0.0
    END DO
  END DO

  DO i=1,3*(nx-1)+15
    a(iwsave2+i-1) = 0.0
  END DO

  DO j=1,3*(ny-1)+15
    a(iwsave1+j-1) = 0.0
  END DO

  IF ( tbc == 4 ) THEN  ! set up the fft work arrays for use in the
                        ! upper radiation boundary condition.
    IF ( fftopt == 1 ) THEN     ! set up periodic work arrays...
      IF ( ny == 4 ) THEN       ! set up trigs in x direction only
        CALL set99(a(itrigs1),a(iifax1),nx-1)    ! NOTE: nx must be ODD!!!!
                                                 ! and of special character...
                                                 ! see fft99f.f for details....
      ELSE IF ( nx == 4 ) THEN     ! set up trigs in y direction only

        CALL set99(a(itrigs2),a(iifax2),ny-1)    ! NOTE: ny must be ODD!!!!
                                                 ! and of special character...
                                                 ! see fft99f.f for details....
      ELSE    ! set up for 2-d transform...

        CALL set99(a(itrigs1),a(iifax1),nx-1)    ! NOTE: nx must be ODD!!!!
                                                 ! and of special character...
                                                 ! see fft99f.f for details....
        CALL set99(a(itrigs2),a(iifax2),ny-1)    ! NOTE: ny must be ODD!!!!
                                                 ! and of special character...
                                                 ! see fft99f.f for details....
      END IF

    ELSE IF ( fftopt == 2 ) THEN   ! set up the cos fft arrays...

      IF(ny == 4)THEN   ! set up function in x direction only...

        CALL vcosti(nx-1,a(iwsave2))         ! nx should be even.

      ELSE IF(nx == 4)THEN ! set up function in y direction only...

        CALL vcosti(ny-1,a(iwsave1))         ! ny should be even.

      ELSE   ! set up functions for 2-d transform...

        CALL vcosti(ny-1,a(iwsave1))         ! ny should be even.
        CALL vcosti(nx-1,a(iwsave2))         ! nx should be even.

      END IF  ! end of run type if block...

    END IF   ! end of fftopt if block.....

  END IF
!
!-----------------------------------------------------------------------
!
!  Fill soiltyp and vegtyp which are integers and can not be
!  interpolated.
!
!-----------------------------------------------------------------------
!
  IF ( sfcphy > 0 ) THEN
    CALL initindex(mptr, mparent)
  END IF
!
!-----------------------------------------------------------------------
!
!  Fill wcont and km with zeros.
!  They are work arrays between time steps, therefore their values
!  here do not matter.
!
!-----------------------------------------------------------------------
!
  CALL filzero(a(iwcont), nxyz)
  CALL filzero(a(ikmh   ), nxyz)
  CALL filzero(a(ikmv   ), nxyz)
  CALL filzero(a(irprntl), nxyz)

!*********************************************************************
! Specify the initial bubble on the new grid using the analytical
! function. This supercedes the value interpolated from base grid.
!
!*********************************************************************
!  if(.true.) then
  IF(.false.) THEN

    IF( curtim == 0.0 ) THEN
!
      PRINT*,' Reinitializing the thermal bubble at time ',curtim,      &
             'for grid ', mptr
!
      iptprt = igtxyz(mptr,id_ptprt,3)

      CALL setbubble(mptr, nx,ny,nz,a(ix),a(iy),a(izp), a(iptprt) )

      CALL retxyz(mptr,id_ptprt,3,iptprt , .true.)

    END IF

  END IF
!
!-----------------------------------------------------------------------
!
!  return the reset arrays to their permanent storage when unpacked.
!
!-----------------------------------------------------------------------
!
  CALL retnxy(mptr,id_mapfct,8,imapfct1,.true.)
  CALL retnxy(mptr,id_vwork1,2,ivwork1, .true.)
  CALL retnxy(mptr,id_vwork2,2,ivwork2, .true.)

  CALL retxyz(mptr,id_zp,   1,izp,   .false.)
  CALL retxyz(mptr,id_j1,   1,ij1,   .true.)
  CALL retxyz(mptr,id_j2,   1,ij2,   .true.)
  CALL retxyz(mptr,id_j3,   1,ij3,   .true.)
  CALL retxyz(mptr,id_j3inv,1,ij3inv,.true.)
  CALL retxyz(mptr,id_aj3x, 1,iaj3x, .true.)
  CALL retxyz(mptr,id_aj3y, 1,iaj3y, .true.)
  CALL retxyz(mptr,id_aj3z, 1,iaj3z, .true.)

  CALL retxyz(mptr,id_wcont, 1,iwcont ,.true.)
  CALL retxyz(mptr,id_kmh,   1,ikmh   ,.true.)
  CALL retxyz(mptr,id_kmv,   1,ikmv   ,.true.)
  CALL retxyz(mptr,id_rprntl,1,irprntl,.true.)

!*********************************************************************
! re-set all tem. space
!*********************************************************************

  CALL resett
!
!-----------------------------------------------------------------------
!
!  Print out the parameters.
!
!  Write out a log file of model parameters which can be used as
!  the input file to re-run the model for this grid.
!
!-----------------------------------------------------------------------
!
  CALL prtlog(nx,ny,nz,0)

  CALL strcnts(nx,ny,nz, a(ii),nsint, a(ir),nsreal)

  PRINT*, '...done with INITNGRD on grid ', mptr

  RETURN
END SUBROUTINE initngrd

!*********************************************************************
!*********************************************************************
!*********************************************************************


SUBROUTINE nwgrdcst(mptr, mparent,                                      &
           iconst,iconstp,nsint, rconst,rconstp,nsreal)
!
!-----------------------------------------------------------------------
!
!  Initialize the constants and constant arrays for new grid mptr.
!
!  AUTHOR: Ming Xue
!
!  11/15/1992
!   Updated: E.J.Adlerman, May 1995, arps4.0.18
!
!  09/26/1996 (Yuhe Liu)
!   Added the isotropic option for divergence damping. Parameter
!   divdmpnd changed to divdmpndh for horizontal and divdmpndv for
!   vertical.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nsint,nsreal
  INTEGER :: iconst(nsint ),iconstp(nsint)
  REAL :: rconst(nsreal),rconstp(nsreal)
!
  INCLUDE 'nodal.inc'
  INCLUDE 'agrialloc.inc'

  INCLUDE 'bndry.inc'       ! Boundary condition control parameters
  INCLUDE 'cumucst.inc'     ! Cumulus param.
  INCLUDE 'globcst.inc'     ! Global constants that control model
  INCLUDE 'phycst.inc'      ! Unchanging physical constants
  INCLUDE 'sfcphycst.inc'   ! Unchanging physical constants
  INCLUDE 'grid.inc'

  INTEGER :: i
  REAL :: dxbase,dybase
  REAL :: temscl
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  PRINT*, '...starting NWGRDCST on grid',  mptr
!
!-----------------------------------------------------------------------
!
!  First, copy the values of constants from the parent grid.
!  Then reset the inappropriate values.
!
!-----------------------------------------------------------------------
!
  DO i=1,nsint
    iconst(i) = iconstp(i)
  END DO

  DO i=1,nsreal
    rconst(i) = rconstp(i)
  END DO

!
!-----------------------------------------------------------------------
!
!  Retrieve constant values from the constant arrays
!
!-----------------------------------------------------------------------
!
  CALL getcnts(nx,ny,nz, iconst,nsint, rconst,nsreal)

  mgrid = mptr
!
!-----------------------------------------------------------------------
!
!  Let's save some of the constants of the parent grid.
!
!-----------------------------------------------------------------------
!

!
! for run mx03a only
!
  nocmnt = 2

  dxp = dx
  dyp = dy
  dzp = dz

  xorigp  = xorig
  yorigp  = yorig
  cfcmh2p = cfcmh2
  cfcmv2p = cfcmv2
  cfcmh4p = cfcmh4
  cfcmv4p = cfcmv4

  cdvdmphp = cdvdmph
  cdvdmpvp = cdvdmpv
  dtsmlp  = dtsml
  dtbigp  = dtbig
!
!-----------------------------------------------------------------------
!
!  Reset some parameters for the new grid.
!
!-----------------------------------------------------------------------
!
  nx = node(5,mptr)
  ny = node(6,mptr)
  nz = node(14,mptr)

  dx = rnode(9 ,mptr)
  dy = rnode(10,mptr)
  dz = dzp

  dxinv = 1.0/dx
  dyinv = 1.0/dy
  dzinv = 1.0/dz

  xl = (nx-3)*dx
  yl = (ny-3)*dy
  zh = (nz-3)*dz

!*********************************************************************
! here we have assumed that for the base grid xorig=yorig=0.0
! which is not necessarily true.
!*********************************************************************

  dxbase = rnode(9 ,mstart)
  dybase = rnode(10,mstart)

  xorig = rnode(1,mptr)-dxbase+dx
  yorig = rnode(2,mptr)-dybase+dy
  zorig = zorig

  dtbig = rnode(11,mptr)

  IF( vimplct == 1 ) THEN
!*********************************************************************
!  For implicit case, dtsml does not depend on dz.
!*********************************************************************
    print*,'dtsmlp,dxp,dyp=', dtsmlp, dxp, dyp

    dtsml = dtsmlp/(1.0/SQRT(1.0/dxp**2+1.0/dyp**2))                    &
                  *(1.0/SQRT(1.0/dx **2+1.0/dy **2))

    print*,'dtsml ,dx ,dy =', dtsml , dx , dy 
  ELSE
    dtsml = dtsmlp/(1.0/SQRT(1.0/dxp**2+1.0/dyp**2+1.0/dzmin**2))       &
                  *(1.0/SQRT(1.0/dx **2+1.0/dy **2+1.0/dzmin**2))
  END IF

  nsmstp = nint( 2*dtbig/dtsml )
  IF( MOD(2*dtbig,dtsml) /= 0.0 ) nsmstp =nsmstp +1
  dtsml = (2*dtbig)/nsmstp

  cfcmh2 = cfcmh2p*(dx*dy)/(dxp*dyp)
!    :         *dtbigp/dtbig
!    :         *dtbigp/dtbig
  cfcmv2 = cfcmv2p*(dz/dzp)**2
!    :         *dtbigp/dtbig
  cfcmh4 = cfcmh4p*((dx*dy)/(dxp*dyp))**2
!    :         *dtbigp/dtbig
!    :         *dtbigp/dtbig
  cfcmv4 = cfcmv4p*(dz/dzp)**4
!    :         *dtbigp/dtbig

  IF ( divdmp == 1 ) THEN    ! isotropic, cdvdmph=cdvdmpv
    IF ( runmod == 1 ) THEN
      temsclp = MIN(dxp,dyp,dzmin)
      temscl  = MIN(dx,dy,dzmin)
    ELSE IF( runmod == 2 ) THEN
      temsclp = MIN(dxp,dzmin)
      temscl  = MIN(dx,dzmin)
    ELSE IF( runmod == 3 ) THEN
      temsclp = MIN(dyp,dzmin)
      temscl  = MIN(dy,dzmin)
    ELSE IF( runmod == 4 ) THEN
      temsclp = dzmin
      temscl  = dzmin
    END IF

    divdmpndh = cdvdmphp / ( temsclp**2 / dtsmlp)
    divdmpndv = divdmpndh

    cdvdmph = divdmpndh * temscl**2 / dtsml
    cdvdmpv = cdvdmph
  ELSE IF ( divdmp == 2 ) THEN
    IF ( runmod == 1 ) THEN
      temsclp = MIN( SQRT(dxp*dyp), 5000.0 )
      temscl  = MIN( SQRT(dx*dy), 5000.0 )
    ELSE IF( runmod == 2 ) THEN
      temsclp = MIN( dxp, 5000.0 )
      temscl  = MIN( dx, 5000.0 )
    ELSE IF( runmod == 3 ) THEN
      temsclp = MIN( dyp, 5000.0 )
      temscl  = MIN( dy, 5000.0 )
    ELSE IF( runmod == 4 ) THEN
      temsclp = dzmin
      temscl  = dzmin
    END IF

    divdmpndh = cdvdmphp / ( temsclp**2 / dtsmlp )
    divdmpndv = cdvdmpvp / ( dzmin**2 / dtsmlp )

    cdvdmph = divdmpndh * temscl**2 / dtsml
    cdvdmpv = divdmpndv * dzmin**2 / dtsml
  END IF

  tstart = curtim
  tstop  = tstop

  nsteps = nint( (tstop-tstart)/dtbig )

!*********************************************************************
!  nsteps is now unused in the model either. It's calculated
!  here for safty reason only.
!*********************************************************************

  node(13,mptr) = 0
  nstep = node(13,mptr)

  restrt = 0

  nfmtprt  = nint(tfmtprt/dtbig)
  IF ( nfmtprt == 0 ) nfmtprt = -1

  nhisdmp  = nint(thisdmp/dtbig)
  IF ( nhisdmp == 0 ) nhisdmp = -1

  nstrtdmp = nint(tstrtdmp/dtbig)

  nrstout  = nint(trstout/dtbig)
  IF ( nrstout == 0 ) nrstout = -1

  nmaxmin  = nint(tmaxmin/dtbig)
  IF ( nmaxmin == 0 ) nmaxmin = -1

  nenergy  = nint(tenergy/dtbig)
  IF ( nenergy == 0 ) nenergy = -1
!
!-----------------------------------------------------------------------
!
!  If this is not a base grid, set the lateral boundary
!  condition to option 6.
!
!-----------------------------------------------------------------------
!
  IF( mptr /= 1 ) THEN  ! if not the base grid
    ebc = 6
    wbc = 6
    sbc = 6
    nbc = 6
  END IF
!
!-----------------------------------------------------------------------
!
!  Store constant variables of the model into constant arrays
!
!-----------------------------------------------------------------------
!
  CALL strcnts(nx,ny,nz, iconst,nsint, rconst,nsreal)

  PRINT *, '...done with NWGRDCST on grid..', mptr



  RETURN
END SUBROUTINE nwgrdcst

!*********************************************************************
!*********************************************************************
!*********************************************************************



SUBROUTINE setbubble(mptr,nx,ny,nz,x,y,zp,ptprt)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Initialize ptprt with a thermal bubble for grid mptr.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  11/30/1992.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    mptr     Pointer to the grid.
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    x        x-coordinate of grid points in computational space (m)
!    y        y-coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space (m)
!
!  OUTPUT:
!
!    ptprt    Perturbation potential temperature at all time levels (K)
!
!-----------------------------------------------------------------------
!
  INTEGER :: mptr
  INTEGER :: nx,ny,nz          ! The number of grid points in 3 directions

  REAL :: x     (nx)
  REAL :: y     (ny)
  REAL :: zp    (nx,ny,nz)
  REAL :: ptprt (nx,ny,nz,3)

  INCLUDE 'globcst.inc'
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  REAL :: xs, ys, zs
  REAL :: beta, pi2
  INTEGER :: i,j,k

  IF( pt0opt == 1 ) THEN  ! Bubble shaped perturbation
!
!-----------------------------------------------------------------------
!
!  Set the potential temperature perturbation inside a bubble volume.
!
!-----------------------------------------------------------------------
!

    pi2 = 1.5707963267949

    IF( ptpert0 == 0.0) THEN

      DO k= 1,nz-1
        DO j= 1,ny-1
          DO i= 1,nx-1
            ptprt(i,j,k,1) = 0.0
            ptprt(i,j,k,2) = 0.0
            ptprt(i,j,k,3) = 0.0
          END DO
        END DO
      END DO

    ELSE

      DO k= 1,nz-1
        DO j= 1,ny-1
          DO i= 1,nx-1

            xs = (x(i)+x(i+1))*0.5
            ys = (y(j)+y(j+1))*0.5
            zs = (zp(i,j,k)+zp(i,j,k+1))*0.5

            IF( pt0rady < 0.0 ) THEN   ! 2-d bubble in x-z plane.

              beta = SQRT( ((xs-pt0ctrx)/pt0radx)**2                    &
                         + ((zs-pt0ctrz)/pt0radz)**2 )

            ELSE IF( pt0radx < 0.0 ) THEN ! 2-d bubble in y-z plane.

              beta = SQRT( ((ys-pt0ctry)/pt0rady)**2                    &
                         + ((zs-pt0ctrz)/pt0radz)**2 )

            ELSE                        ! 3-d bubble

              beta = SQRT( ((xs-pt0ctrx)/pt0radx)**2 +                  &
                           ((ys-pt0ctry)/pt0rady)**2 +                  &
                           ((zs-pt0ctrz)/pt0radz)**2 )
            END IF

            IF(beta >= 1.0) THEN
              ptprt(i,j,k,1) = 0.0
            ELSE
              ptprt(i,j,k,1) = ptpert0*(COS(pi2*beta)**2)
            END IF

            ptprt(i,j,k,2) = ptprt(i,j,k,1)
            ptprt(i,j,k,3) = ptprt(i,j,k,1)

          END DO
        END DO
      END DO

    END IF

  END IF

  RETURN
END SUBROUTINE setbubble

SUBROUTINE getdzmin( zp, nx,ny,nz, dzpmin)

  INTEGER :: nx,ny,nz
  REAL :: zp(nx,ny,nz)

  dzpmin = zp(1,1,3)-zp(1,1,2)
  PRINT *, 'zp1 and zp2', zp(1,1,3), zp(1,1,2)

  DO k=2,nz-2
    DO j=1,ny-1
      DO i=1,nx-1
        dzpmin=MIN(dzpmin, zp(i,j,k+1)-zp(i,j,k))
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE getdzmin
