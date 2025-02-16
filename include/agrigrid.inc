!***************************************************************
!
! *************** data for each grid *******************
!
!  for ARPS 4.0.18 (Ming Xue, 10/25/1992)
!               (updated 5/1995 by Edwin Adlerman)
!  nsint  - number of integer constants
!  nsreal - number of real constants
!  nslogc - number of logical constants
!  nx1d   - number of 1-d arrays in x
!  ny1d   - number of 1-d arrays in y
!  nz1d   - number of 1-d arrays in z
!  nxy2d  - number of 2-d xy arrays
!  nxz2d  - number of 2-d xz arrays
!  nyz2d  - number of 2-d yz arrays
!  nxyz3d - number of 3-d xyz arrays
!  nexbc3d- number of 3-d EXBC arrays which are only applied to base grid
!
!***************************************************************

  INTEGER :: nsint,nsreal,nslogc,ns1d,                                  &
          nx1d,ny1d,nz1d,nxy2d,nxz2d,nyz2d,nxyz3d,nexbc3d

  PARAMETER( nsint  = 10000,                                            &
             nsreal = 10000,                                            &
             nslogc = 0,                                                &
             ns1d   = 1000,                                             &
             nx1d   = 10,                                               &
             ny1d   = 10,                                               &
             nz1d   = 1,                                                &
             nxy2d  = 86,                                               &
             nxz2d  = 10,                                               &
             nyz2d  = 10,                                               &
             nxyz3d = 68,                                               &
             nexbc3d= 22 )
!
!  so, this is the story concerning the variables for this
!  solver
!
!----------------------------------------------------------------------
!  interger constants (see list in subroutine SETCNTS).
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!  real constants (see list in subroutine SETCNTS).
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!  logical constants (see list in subroutine SETCNTS).
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!  ns1d 1-d fixed length arrays:
!----------------------------------------------------------------------

  INTEGER :: id_ifax1,id_ifax2

  PARAMETER(id_ifax1 = 1)              ! 1-13  ifax1(13)
  PARAMETER(id_ifax2 = id_ifax1 + 13)  ! 14-27 ifax2(13)

!----------------------------------------------------------------------
!  nx1d 1-d arrays in x
!----------------------------------------------------------------------

  INTEGER :: id_x,id_trigs1,id_wsave2

  PARAMETER(id_x = 1)                  ! 1     x(nx)
  PARAMETER(id_trigs1 = id_x + 1)      ! trigs1(3*(nx-1)/2+1), 2-4
                                       ! Use 3 1-D arrays to store 1
                                       ! array with size 3*(nx-1)/2+1
  PARAMETER(id_wsave2 = id_trigs1 + 3) ! wsave2(3*(nx-1)+15), 5-10
                                       ! Use 6 1-D arrays to store 1
                                       ! array with size 3*(nx-1)+15

!----------------------------------------------------------------------
!  ny1d 1-d arrays in y:
!----------------------------------------------------------------------

  INTEGER :: id_y,id_trigs2,id_wsave1

  PARAMETER(id_y = 1)                  ! 1     y(ny)

  PARAMETER(id_trigs2 = id_y + 1)      ! 2-4   trigs2(3*(nx-1)/2+1)
                                       ! Use 3 1-D arrays to store
                                       ! 1 array with size 3*(ny-1)/2+1
  PARAMETER(id_wsave1 = id_trigs2 + 3) ! 5-10  wsave1(3*(nx-1)+15)
                                       ! Use 6 1-D arrays to store
                                       ! 1 array with size 3*(ny-1)+15

!----------------------------------------------------------------------
!  nz1d 1-D arrays in z:
!----------------------------------------------------------------------

  INTEGER :: id_z

  PARAMETER(id_z = 1)                  ! 1     z(nz)

!----------------------------------------------------------------------
!  2-d xy arrays:
!----------------------------------------------------------------------

  INTEGER :: id_hterain,id_sinlat,id_soiltyp,id_stypfrct,               &
          id_vegtyp,id_lai,id_roufns,id_veg,id_tsfc,                    &
          id_tsoil,id_wetsfc,id_wetdp,id_wetcanp,id_snowdpth,           &
          id_qvsfc,id_raing,id_rainc,id_ptsfc,                          &
          id_pbldpth,id_mapfct,id_rad2d,id_radsw,                       &
          id_rnflx,id_prcrate,id_bcrlx,id_vwork1,                       &
          id_vwork2,id_raincv,id_nca,id_usflx,id_vsflx,                 &
          id_ptsflx,id_qvsflx

  PARAMETER(id_hterain  = 1)               ! 1     hterain (nx,ny)
  PARAMETER(id_sinlat   = id_hterain  + 1) ! 2     sinlat  (nx,ny)
  PARAMETER(id_soiltyp  = id_sinlat   + 1) ! 3-6   soiltyp (nx,ny,4)
  PARAMETER(id_stypfrct = id_soiltyp  + 4) ! 7-10  stypfrct(nx,ny,4)
  PARAMETER(id_vegtyp   = id_stypfrct + 4) ! 11    vegtyp  (nx,ny)
  PARAMETER(id_lai      = id_vegtyp   + 1) ! 12    lai     (nx,ny)
  PARAMETER(id_roufns   = id_lai      + 1) ! 13    roufns  (nx,ny)
  PARAMETER(id_veg      = id_roufns   + 1) ! 14    veg     (nx,ny)
  PARAMETER(id_tsfc     = id_veg      + 1) ! 15-19 tsfc    (nx,ny,0:4)
  PARAMETER(id_tsoil    = id_tsfc     + 5) ! 20-24 tsoil   (nx,ny,0:4)
  PARAMETER(id_wetsfc   = id_tsoil    + 5) ! 25-29 wetsfc  (nx,ny,0:4)
  PARAMETER(id_wetdp    = id_wetsfc   + 5) ! 30-34 wetdp   (nx,ny,0:4)
  PARAMETER(id_wetcanp  = id_wetdp    + 5) ! 35-39 wetcanp (nx,ny,0:4)
  PARAMETER(id_qvsfc    = id_wetcanp  + 5) ! 40-44 qvsfc   (nx,ny,0:4)
  PARAMETER(id_raing    = id_qvsfc    + 5) ! 45    raing   (nx,ny)
  PARAMETER(id_rainc    = id_raing    + 1) ! 46    rainc   (nx,ny)
  PARAMETER(id_ptsfc    = id_rainc    + 1) ! 47    ptsfc   (nx,ny)
  PARAMETER(id_pbldpth  = id_ptsfc    + 1) ! 48-50 pbldpth (nx,ny,3)
  PARAMETER(id_mapfct   = id_pbldpth  + 3) ! 51-58 mapfct  (nx,ny,8)
  PARAMETER(id_rad2d    = id_mapfct   + 8) ! 59-68 rad2d   (nx,ny,10)
  PARAMETER(id_radsw    = id_rad2d    +10) ! 69    radsw   (nx,ny)
  PARAMETER(id_rnflx    = id_radsw    + 1) ! 70    rnflx   (nx,ny)
  PARAMETER(id_prcrate  = id_rnflx    + 1) ! 71-74 prcrate (nx,ny,4)
  PARAMETER(id_bcrlx    = id_prcrate  + 4) ! 75    prcrate (nx,ny,4)
  PARAMETER(id_vwork1   = id_bcrlx    + 1) ! 76-77 vwork1  (nx+1,ny+1)
  PARAMETER(id_vwork2   = id_vwork1   + 2) ! 78-79 vwork2  (nx+1,ny+1)
  PARAMETER(id_raincv   = id_vwork2   + 2) ! 80    raincv  (nx,ny)
  PARAMETER(id_nca      = id_raincv   + 1) ! 81    nca     (nx,ny)
  PARAMETER(id_usflx    = id_nca      + 1) ! 82    usflx   (nx,ny)
  PARAMETER(id_vsflx    = id_usflx    + 1) ! 83    vsflx   (nx,ny)
  PARAMETER(id_ptsflx   = id_vsflx    + 1) ! 84    ptsflx  (nx,ny)
  PARAMETER(id_qvsflx   = id_ptsflx   + 1) ! 85    qvsflx  (nx,ny)
  PARAMETER(id_snowdpth = id_qvsflx   + 1) ! 86    snowdpth(nx,ny)

!----------------------------------------------------------------------
!  nxz2d 2-d xz arrays
!----------------------------------------------------------------------

  INTEGER :: id_udtnb,id_udtsb,                                         &
          id_vdtnb,id_vdtsb,                                            &
          id_sdtnb,id_sdtsb,                                            &
          id_pdtnb,id_pdtsb,                                            &
          id_wdtnb,id_wdtsb

  PARAMETER(id_udtnb = 1)               ! 1     udtnb(nx,nz)
  PARAMETER(id_udtsb = id_udtnb + 1)    ! 2     udtsb(nx,nz)
  PARAMETER(id_vdtnb = id_udtsb + 1)    ! 3     vdtnb(nx,nz)
  PARAMETER(id_vdtsb = id_vdtnb + 1)    ! 4     vdtsb(nx,nz)
  PARAMETER(id_sdtnb = id_vdtsb + 1)    ! 5     sdtnb(nx,nz)
  PARAMETER(id_sdtsb = id_sdtnb + 1)    ! 6     sdtsb(nx,nz)
  PARAMETER(id_pdtnb = id_sdtsb + 1)    ! 7     pdtnb(nx,nz)
  PARAMETER(id_pdtsb = id_pdtnb + 1)    ! 8     pdtsb(nx,nz)
  PARAMETER(id_wdtnb = id_pdtsb + 1)    ! 9     wdtnb(nx,nz)
  PARAMETER(id_wdtsb = id_wdtnb + 1)    ! 10    wdtsb(nx,nz)

!----------------------------------------------------------------------
!  nyz2d 2-d yz arrays:
!----------------------------------------------------------------------

  INTEGER :: id_udteb,id_udtwb,                                         &
          id_vdteb,id_vdtwb,                                            &
          id_sdteb,id_sdtwb,                                            &
          id_pdteb,id_pdtwb,                                            &
          id_wdteb,id_wdtwb

  PARAMETER(id_udteb = 1)               ! 1     udteb(nx,nz)
  PARAMETER(id_udtwb = id_udteb + 1)    ! 2     udtwb(nx,nz)
  PARAMETER(id_vdteb = id_udtwb + 1)    ! 3     vdteb(nx,nz)
  PARAMETER(id_vdtwb = id_vdteb + 1)    ! 4     vdtwb(nx,nz)
  PARAMETER(id_sdteb = id_vdtwb + 1)    ! 5     sdteb(nx,nz)
  PARAMETER(id_sdtwb = id_sdteb + 1)    ! 6     sdtwb(nx,nz)
  PARAMETER(id_pdteb = id_sdtwb + 1)    ! 7     pdteb(nx,nz)
  PARAMETER(id_pdtwb = id_pdteb + 1)    ! 8     pdtwb(nx,nz)
  PARAMETER(id_wdteb = id_pdtwb + 1)    ! 9     wdteb(nx,nz)
  PARAMETER(id_wdtwb = id_wdteb + 1)    ! 10    wdtwb(nx,nz)

!----------------------------------------------------------------------
!  nxyz3d 3-d xyz arrays
!----------------------------------------------------------------------

  INTEGER :: id_u,id_v,id_w

  PARAMETER(id_u = 1)                ! 1-3   u    (nx,ny,nz,3)
  PARAMETER(id_v = id_u + 3)         ! 4-6   v    (nx,ny,nz,3)
  PARAMETER(id_w = id_v + 3)         ! 7-9   w    (nx,ny,nz,3)

!***********************************************************************
!  The following variables, ptprt,pprt, qv,qc,qr,qi,qs,qh, tke, have
!  to stick together because they will be configured by do-loops.
!***********************************************************************

  INTEGER :: id_ptprt,id_pprt,id_qv,id_qc,id_qr,                        &
          id_qi,id_qs,id_qh,id_tke

  PARAMETER(id_ptprt    = id_w        + 3)  ! 10-12 ptprt(nx,ny,nz,3)
  PARAMETER(id_pprt     = id_ptprt    + 3)  ! 13-15 pprt (nx,ny,nz,3)
  PARAMETER(id_qv       = id_pprt     + 3)  ! 16-18 qv   (nx,ny,nz,3)
  PARAMETER(id_qc       = id_qv       + 3)  ! 19-21 qc   (nx,ny,nz,3)
  PARAMETER(id_qr       = id_qc       + 3)  ! 22-24 qr   (nx,ny,nz,3)
  PARAMETER(id_qi       = id_qr       + 3)  ! 25-27 qi   (nx,ny,nz,3)
  PARAMETER(id_qs       = id_qi       + 3)  ! 28-30 qs   (nx,ny,nz,3)
  PARAMETER(id_qh       = id_qs       + 3)  ! 31-33 qh   (nx,ny,nz,3)
  PARAMETER(id_tke      = id_qh       + 3)  ! 34-36 tke  (nx,ny,nz,3)

  INTEGER :: id_ubar,id_vbar,id_wbar,id_ptbar,id_pbar,                  &
          id_rhostr,id_qvbar,id_zp,id_j1,id_j2,id_j3,                   &
          id_wcont,id_kmh,id_kmv,id_rprntl,id_ppi,                      &
          id_csndsq,id_ptcumsrc,id_qcumsrc,id_radfrc,                   &
          id_j3inv,id_w0avg,id_aj3x,id_aj3y,id_aj3z,                    &
          id_ptbari,id_pbari,id_rhostri

  PARAMETER(id_ubar     = id_tke      + 3)  ! 37    ubar    (nx,ny,nz)
  PARAMETER(id_vbar     = id_ubar     + 1)  ! 38    vbar    (nx,ny,nz)
  PARAMETER(id_wbar     = id_vbar     + 1)  ! 39    wbar    (nx,ny,nz)
  PARAMETER(id_ptbar    = id_wbar     + 1)  ! 40    ptbar   (nx,ny,nz)
  PARAMETER(id_pbar     = id_ptbar    + 1)  ! 41    pbar    (nx,ny,nz)
  PARAMETER(id_rhostr   = id_pbar     + 1)  ! 42    rhostr  (nx,ny,nz)
  PARAMETER(id_qvbar    = id_rhostr   + 1)  ! 43    qvbar   (nx,ny,nz)
  PARAMETER(id_zp       = id_qvbar    + 1)  ! 44    zp      (nx,ny,nz)
  PARAMETER(id_j1       = id_zp       + 1)  ! 45    j1      (nx,ny,nz)
  PARAMETER(id_j2       = id_j1       + 1)  ! 46    j2      (nx,ny,nz)
  PARAMETER(id_j3       = id_j2       + 1)  ! 47    j3      (nx,ny,nz)
  PARAMETER(id_wcont    = id_j3       + 1)  ! 48    wcont   (nx,ny,nz)
  PARAMETER(id_kmh      = id_wcont    + 1)  ! 49    kmh     (nx,ny,nz)
  PARAMETER(id_kmv      = id_kmh      + 1)  ! 50    kmv     (nx,ny,nz)
  PARAMETER(id_rprntl   = id_kmv      + 1)  ! 51    rprntl  (nx,ny,nz)
  PARAMETER(id_ppi      = id_rprntl   + 1)  ! 52    ppi     (nx,ny,nz)
  PARAMETER(id_csndsq   = id_ppi      + 1)  ! 53    csndsq  (nx,ny,nz)
  PARAMETER(id_ptcumsrc = id_csndsq   + 1)  ! 54    ptcumsrc(nx,ny,nz)
  PARAMETER(id_qcumsrc  = id_ptcumsrc + 1)  ! 55    qcumsrc (nx,ny,nz,5)
  PARAMETER(id_radfrc   = id_qcumsrc  + 5)  ! 60    radfrc  (nx,ny,nz)
  PARAMETER(id_j3inv    = id_radfrc   + 1)  ! 61    j3inv   (nx,ny,nz)
  PARAMETER(id_w0avg    = id_j3inv    + 1)  ! 62    w0avg   (nx,ny,nz)
  PARAMETER(id_aj3x     = id_w0avg    + 1)  ! 63    aj3x    (nx,ny,nz)
  PARAMETER(id_aj3y     = id_aj3x     + 1)  ! 64    aj3y    (nx,ny,nz)
  PARAMETER(id_aj3z     = id_aj3y     + 1)  ! 65    aj3z    (nx,ny,nz)
  PARAMETER(id_ptbari   = id_aj3z     + 1)  ! 66    ptbari  (nx,ny,nz)
  PARAMETER(id_pbari    = id_ptbari   + 1)  ! 67    pbari   (nx,ny,nz)
  PARAMETER(id_rhostri  = id_pbari    + 1)  ! 68    rhostri (nx,ny,nz)

!----------------------------------------------------------------------
!  nexbc3d 3-d xyz arrays for EXBC
!----------------------------------------------------------------------

  INTEGER :: id_u0exb,id_v0exb,id_w0exb,id_pt0exb,id_pr0exb,            &
          id_qv0exb,id_qc0exb,id_qr0exb,                                &
          id_qi0exb,id_qs0exb,id_qh0exb

  PARAMETER(id_u0exb  = 1)               ! 1     u0exb (nx,ny,nz)
  PARAMETER(id_v0exb  = id_u0exb  + 1)   ! 2     v0exb (nx,ny,nz)
  PARAMETER(id_w0exb  = id_v0exb  + 1)   ! 3     w0exb (nx,ny,nz)
  PARAMETER(id_pt0exb = id_w0exb  + 1)   ! 4     pt0exb(nx,ny,nz)
  PARAMETER(id_pr0exb = id_pt0exb + 1)   ! 5     pr0exb(nx,ny,nz)
  PARAMETER(id_qv0exb = id_pr0exb + 1)   ! 6     qv0exb(nx,ny,nz)
  PARAMETER(id_qc0exb = id_qv0exb + 1)   ! 7     qc0exb(nx,ny,nz)
  PARAMETER(id_qr0exb = id_qc0exb + 1)   ! 8     qr0exb(nx,ny,nz)
  PARAMETER(id_qi0exb = id_qr0exb + 1)   ! 9     qi0exb(nx,ny,nz)
  PARAMETER(id_qs0exb = id_qi0exb + 1)   ! 10    qs0exb(nx,ny,nz)
  PARAMETER(id_qh0exb = id_qs0exb + 1)   ! 11    qh0exb(nx,ny,nz)

  INTEGER :: id_udtexb,id_vdtexb,id_wdtexb,id_ptdtexb,id_prdtexb,       &
          id_qvdtexb,id_qcdtexb,id_qrdtexb,                             &
          id_qidtexb,id_qsdtexb,id_qhdtexb

  PARAMETER(id_udtexb  = id_qh0exb  + 1) ! 12    udtexb (nx,ny,nz)
  PARAMETER(id_vdtexb  = id_udtexb  + 1) ! 13    vdtexb (nx,ny,nz)
  PARAMETER(id_wdtexb  = id_vdtexb  + 1) ! 14    wdtexb (nx,ny,nz)
  PARAMETER(id_ptdtexb = id_wdtexb  + 1) ! 15    ptdtexb(nx,ny,nz)
  PARAMETER(id_prdtexb = id_ptdtexb + 1) ! 16    prdtexb(nx,ny,nz)
  PARAMETER(id_qvdtexb = id_prdtexb + 1) ! 17    qvdtexb(nx,ny,nz)
  PARAMETER(id_qcdtexb = id_qvdtexb + 1) ! 18    qcdtexb(nx,ny,nz)
  PARAMETER(id_qrdtexb = id_qcdtexb + 1) ! 19    qrdtexb(nx,ny,nz)
  PARAMETER(id_qidtexb = id_qrdtexb + 1) ! 20    qidtexb(nx,ny,nz)
  PARAMETER(id_qsdtexb = id_qidtexb + 1) ! 21    qsdtexb(nx,ny,nz)
  PARAMETER(id_qhdtexb = id_qsdtexb + 1) ! 22    qhdtexb(nx,ny,nz)
