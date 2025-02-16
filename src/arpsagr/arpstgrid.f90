SUBROUTINE stgrid

  INCLUDE 'agrigrid.inc'
  INCLUDE 'grddsc.inc'

  INTEGER :: iv

!***********************************************************************
!  This routine sets the generic grid description data
!  in the common block grddsc.  The user must alter this
!  when putting in a new solver or new equations.
!
!    UPDATED: August 1995, by E.J.Adlerman for arps4.0.22
!
!  They have the conventions (e.g. for 3-d arrays):
!
!  ivar = 1
!  stgxyz(1,ivar) = 0.    grid staggerring relative to grid origin in x.
!  stgxyz(2,ivar) = 0.5   grid staggerring relative to grid origin in y.
!  stgxyz(3,ivar) = 0.5   grid staggerring relative to grid origin in z.
!  idmxyz(1,ivar) = 0     dimension off set in x dir. i=1,nx-idmxyz(1,ivar).
!  idmxyz(2,ivar) = 1     dimension off set in y dir. j=1,ny-idmyyz(2,ivar).
!  idmxyz(3,ivar) = 1     dimension off set in z dir. k=1,nz-idmzyz(3,ivar).
!  ipkxyz(  ivar) = 1     word packing factor of this array
!  inixyz(1,ivar) = 1     0 or 1 to determine if this array needs init. on
!                      new grid.
!  inixyz(2,ivar) = 5     when init., 0 for scalar, or the other comp. of
!                      the vector, positive if pointing to the comp. in
!                      y-dir. of the grid, negative if pointing back to
!                      comp. in x-dir.
!  inixyz(3,ivar) = 1     The corresponding array at time dt earlier.
!                      if 0, no earlier array available.
!  iupxyz(1,ivar) = 1     0 or 1 to determine if this array needs updating
!                      from fine grid(s).
!  iupxyz(2,ivar) = 5     When updating, 0 for scalar, or the variable number
!                      of the other comp. of a vector. sign convention
!                      same as for inixyz(2,ivar).
!  ibdxyz(1,ivar) = 1     0 or 1 to determine if this array needs B.C.
!                      by interpolation.
!  ibdxyz(2,ivar) = 5     When interpolating for BC, 0 for scalar, or the
!                      variable number of the other comp. of a vector.
!                       sign convention same as for inixyz(2,ivar).
!  ibdxyz(3,ivar) = 1     The var. number of IVAR at the previous time level.
!***********************************************************************

!***********************************************************************
!  set 2d xy grid info
!***********************************************************************

  DO iv = 1, nxy2d   ! hterain, sinlat, soiltyp, vegtyp

    stgxy(1,iv) = 0.5    ! lai, roufns, veg, tsfc, qvsfc, tsoil
    stgxy(2,iv) = 0.5    ! wetsfc, wetdp, wetcanp, raing, rainc
    idmxy(1,iv) = 1      ! stabchk, ptsfc, snowcvr
    idmxy(2,iv) = 1
    ipkxy(  iv) = 1
    inixy(1,iv) = 1
    inixy(2,iv) = 0
    inixy(3,iv) = 0
    iupxy(1,iv) = 0
    iupxy(2,iv) = 0
    ibdxy(1,iv) = 0
    ibdxy(2,iv) = 0
    ibdxy(3,iv) = 0

  END DO

!
!-----------------------------------------------------------------------
!
!  Don't initialize soiltyp and vegtyp in normal procedure because
!  these two are integers and can not be interpolated. Special care
!  has to be given to the two variables in their initialization.
!
!  The map factors also need to be taken care of because they are
!  for three different grid points, i.e., u, v and scalar points.
!
!  See subroutine nwgrdcst in file arpsinit.f
!
!-----------------------------------------------------------------------
!
  inixy(1,id_soiltyp  )  = 0     ! soiltyp(:,:,1)
  inixy(1,id_soiltyp+1)  = 0     ! soiltyp(:,:,2)
  inixy(1,id_soiltyp+2)  = 0     ! soiltyp(:,:,3)
  inixy(1,id_soiltyp+3)  = 0     ! soiltyp(:,:,4)
  inixy(1,id_stypfrct  ) = 0     ! stypfrct(:,:,1)
  inixy(1,id_stypfrct+1) = 0     ! stypfrct(:,:,2)
  inixy(1,id_stypfrct+2) = 0     ! stypfrct(:,:,3)
  inixy(1,id_stypfrct+3) = 0     ! stypfrct(:,:,4)
  inixy(1,id_vegtyp   )  = 0     ! vegtyp

  inixy(1,id_mapfct  ) = 0      ! mapfct    at scalar points
  inixy(1,id_mapfct+6) = 0      ! mapfct**2 at scalar points
  inixy(1,id_mapfct+7) = 0      ! mapfct/2  at scalar points

  stgxy(1,id_mapfct+1) = 0.0    ! mapfct at u-points
  idmxy(1,id_mapfct+1) = 0      !
  inixy(1,id_mapfct+1) = 0

  stgxy(2,id_mapfct+2) = 0.0    ! mapfct at v-points
  idmxy(2,id_mapfct+2) = 0      !
  inixy(1,id_mapfct+2) = 0

  inixy(1,id_mapfct+3) = 0      ! 1/mapfct at scalar points

  stgxy(1,id_mapfct+4) = 0.0    ! 1/mapfct at u-points
  idmxy(1,id_mapfct+4) = 0      !
  inixy(1,id_mapfct+4) = 0

  stgxy(2,id_mapfct+5) = 0.0    ! 1/mapfct at v-points
  idmxy(2,id_mapfct+5) = 0      !
  inixy(1,id_mapfct+5) = 0

  IF(nxz2d > 0) THEN

!***********************************************************************
!  set 2d xz grid info
!***********************************************************************

    DO iv = 1,nxz2d
      stgxz(1,iv) = 0.  ! not correct for all, but does not matter.
      stgxz(2,iv) = 0.5 ! udtnb, udtsb, vdtnb, vdtsb, sdtnb, sdtsb
      idmxz(1,iv) = 0   ! pdtnb, pdtsb, wdtnb, wdtsb
      idmxz(2,iv) = 0
      ipkxz(  iv) = 1
      inixz(1,iv) = 0
      inixz(2,iv) = 0
      inixz(3,iv) = 0
      iupxz(1,iv) = 0
      iupxz(2,iv) = 0
      ibdxz(1,iv) = 0
      ibdxz(2,iv) = 0
      ibdxz(3,iv) = 0
    END DO
  END IF

  IF(nyz2d > 0) THEN

!***********************************************************************
!  set 2d yz grid info
!***********************************************************************

    DO iv = 1,nyz2d
      stgyz(1,iv) = 0.  ! not correct for all, but does not matter.
      stgyz(2,iv) = 0.5 ! udteb, udtwb, vdteb, vdtwb, sdteb, sdtwb
      idmyz(1,iv) = 0   ! pdteb, pdtwb, wdteb, wdtwb
      idmyz(2,iv) = 0
      ipkyz(  iv) = 1
      iniyz(1,iv) = 0
      iniyz(2,iv) = 0
      iniyz(3,iv) = 0
      iupyz(1,iv) = 0
      iupyz(2,iv) = 0
      ibdyz(1,iv) = 0
      ibdyz(2,iv) = 0
      ibdyz(3,iv) = 0
    END DO
  END IF



  IF( nxyz3d <= 0) THEN
    RETURN
  ELSE

!***********************************************************************
!  set the 3-d grid information
!
!  variable 1   u(*,*,*,1)
!***********************************************************************

    iv = id_u
    stgxyz(1,iv) = 0.0
    stgxyz(2,iv) = 0.5
    stgxyz(3,iv) = 0.5
    idmxyz(1,iv) = 0
    idmxyz(2,iv) = 1
    idmxyz(3,iv) = 1
    ipkxyz(iv)   = 1
    inixyz(1,iv) = 1
    inixyz(2,iv) = id_v
    inixyz(3,iv) = 0
    iupxyz(1,iv) = 0
    iupxyz(2,iv) = 0
    ibdxyz(1,iv) = 0
    ibdxyz(2,iv) = 0
    ibdxyz(3,iv) = 0

!***********************************************************************
!  variable 2   u(*,*,*,2)
!***********************************************************************

    iv = id_u + 1
    stgxyz(1,iv) = 0.0
    stgxyz(2,iv) = 0.5
    stgxyz(3,iv) = 0.5
    idmxyz(1,iv) = 0
    idmxyz(2,iv) = 1
    idmxyz(3,iv) = 1
    ipkxyz(iv)   = 1
    inixyz(1,iv) = 1
    inixyz(2,iv) = id_v+1
    inixyz(3,iv) = id_u
    iupxyz(1,iv) = 1
    iupxyz(2,iv) = id_v+1
    ibdxyz(1,iv) = 1
    ibdxyz(2,iv) = id_v+1
    ibdxyz(3,iv) = id_u

!***********************************************************************
!  variable 3   u(*,*,*,3)
!***********************************************************************

    iv = id_u + 2
    stgxyz(1,iv) = 0.0
    stgxyz(2,iv) = 0.5
    stgxyz(3,iv) = 0.5
    idmxyz(1,iv) = 0
    idmxyz(2,iv) = 1
    idmxyz(3,iv) = 1
    ipkxyz(iv)   = 1
    inixyz(1,iv) = 0
    inixyz(2,iv) = 0
    inixyz(3,iv) = 0
    iupxyz(1,iv) = 0
    iupxyz(2,iv) = 0
    ibdxyz(1,iv) = 0
    ibdxyz(2,iv) = 0
    ibdxyz(3,iv) = 0

!***********************************************************************
!  variable 4   v(*,*,*,1)
!***********************************************************************

    iv = id_v
    stgxyz(1,iv) = 0.5
    stgxyz(2,iv) = 0.0
    stgxyz(3,iv) = 0.5
    idmxyz(1,iv) = 1
    idmxyz(2,iv) = 0
    idmxyz(3,iv) = 1
    ipkxyz(iv)   = 1
    inixyz(1,iv) = 1
    inixyz(2,iv) = -id_u
    inixyz(3,iv) = 0
    iupxyz(1,iv) = 0
    iupxyz(2,iv) = 0
    ibdxyz(1,iv) = 0
    ibdxyz(2,iv) = 0
    ibdxyz(3,iv) = 0

!***********************************************************************
!  variable 5   v(*,*,*,2)
!***********************************************************************

    iv = id_v + 1
    stgxyz(1,iv) = 0.5
    stgxyz(2,iv) = 0.0
    stgxyz(3,iv) = 0.5
    idmxyz(1,iv) = 1
    idmxyz(2,iv) = 0
    idmxyz(3,iv) = 1
    ipkxyz(iv)   = 1
    inixyz(1,iv) = 1
    inixyz(2,iv) = -(id_u+1)
    inixyz(3,iv) = id_v
    iupxyz(1,iv) = 1
    iupxyz(2,iv) = -(id_u+1)
    ibdxyz(1,iv) = 1
    ibdxyz(2,iv) = -(id_u+1)
    ibdxyz(3,iv) = id_v

!***********************************************************************
!  variable 6   v(*,*,*,3)
!***********************************************************************

    iv = id_v + 2
    stgxyz(1,iv) = 0.5
    stgxyz(2,iv) = 0.0
    stgxyz(3,iv) = 0.5
    idmxyz(1,iv) = 1
    idmxyz(2,iv) = 0
    idmxyz(3,iv) = 1
    ipkxyz(iv)   = 1
    inixyz(1,iv) = 0
    inixyz(2,iv) = 0
    inixyz(3,iv) = 0
    iupxyz(1,iv) = 0
    iupxyz(2,iv) = 0
    ibdxyz(1,iv) = 0
    ibdxyz(2,iv) = 0
    ibdxyz(3,iv) = 0

!***********************************************************************
!  variable 7   w(*,*,*,1)
!***********************************************************************

    iv = id_w
    stgxyz(1,iv) = 0.5
    stgxyz(2,iv) = 0.5
    stgxyz(3,iv) = 0.0
    idmxyz(1,iv) = 1
    idmxyz(2,iv) = 1
    idmxyz(3,iv) = 0
    ipkxyz(iv)   = 1
    inixyz(1,iv) = 1
    inixyz(2,iv) = 0
    inixyz(3,iv) = 0
    iupxyz(1,iv) = 0
    iupxyz(2,iv) = 0
    ibdxyz(1,iv) = 0
    ibdxyz(2,iv) = 0
    ibdxyz(3,iv) = 0

!***********************************************************************
!  variable 8   w(*,*,*,2)
!***********************************************************************

    iv = id_w + 1
    stgxyz(1,iv) = 0.5
    stgxyz(2,iv) = 0.5
    stgxyz(3,iv) = 0.0
    idmxyz(1,iv) = 1
    idmxyz(2,iv) = 1
    idmxyz(3,iv) = 0
    ipkxyz(iv)   = 1
    inixyz(1,iv) = 1
    inixyz(2,iv) = 0
    inixyz(3,iv) = id_w
    iupxyz(1,iv) = 1
    iupxyz(2,iv) = 0
    ibdxyz(1,iv) = 1
    ibdxyz(2,iv) = 0
    ibdxyz(3,iv) = id_w

!***********************************************************************
!  variable 9   w(*,*,*,3)
!***********************************************************************

    iv = id_w + 2
    stgxyz(1,iv) = 0.5
    stgxyz(2,iv) = 0.5
    stgxyz(3,iv) = 0.0
    idmxyz(1,iv) = 1
    idmxyz(2,iv) = 1
    idmxyz(3,iv) = 0
    ipkxyz(iv)   = 1
    inixyz(1,iv) = 0
    inixyz(2,iv) = 0
    inixyz(3,iv) = 0
    iupxyz(1,iv) = 0
    iupxyz(2,iv) = 0
    ibdxyz(1,iv) = 0
    ibdxyz(2,iv) = 0
    ibdxyz(3,iv) = 0

!***********************************************************************
!  variables ptprt, pprt, qv,qc,qr,qi,qs,qh, tke at t-dt
!  (assume these variables' id pointer numbers are ranged together)
!***********************************************************************

    DO iv = id_ptprt,id_tke,3

      stgxyz(1,iv) = 0.5
      stgxyz(2,iv) = 0.5
      stgxyz(3,iv) = 0.5
      idmxyz(1,iv) = 1
      idmxyz(2,iv) = 1
      idmxyz(3,iv) = 1
      ipkxyz(iv)   = 1
      inixyz(1,iv) = 1
      inixyz(2,iv) = 0
      inixyz(3,iv) = 0
      iupxyz(1,iv) = 0
      iupxyz(2,iv) = 0
      ibdxyz(1,iv) = 0
      ibdxyz(2,iv) = 0
      ibdxyz(3,iv) = 0

    END DO

!  inixyz(1,id_tke) = 0  ! adjust for TKE scheme

!***********************************************************************
!  variables ptprt, pprt, qv,qc,qr,qi,qs,qh, tke at t
!***********************************************************************

    DO iv =id_ptprt+1,id_tke+1,3

      stgxyz(1,iv) = 0.5
      stgxyz(2,iv) = 0.5
      stgxyz(3,iv) = 0.5
      idmxyz(1,iv) = 1
      idmxyz(2,iv) = 1
      idmxyz(3,iv) = 1
      ipkxyz(iv)   = 1
      inixyz(1,iv) = 1
      inixyz(2,iv) = 0
      inixyz(3,iv) = iv-1
      iupxyz(1,iv) = 1
      iupxyz(2,iv) = 0
      ibdxyz(1,iv) = 1
      ibdxyz(2,iv) = 0
      ibdxyz(3,iv) = iv-1

    END DO

!  inixyz(1,id_tke+1) = 0   ! adjust for TKE scheme
!  iupxyz(1,id_tke+1) = 0
!  ibdxyz(1,id_tke+1) = 0

!***********************************************************************
!  variables ptprt, pprt, qv,qc,qr,qi,qs,qh, tke at t+dt
!***********************************************************************

    DO iv =id_ptprt+2,id_tke+2,3

      stgxyz(1,iv) = 0.5
      stgxyz(2,iv) = 0.5
      stgxyz(3,iv) = 0.5
      idmxyz(1,iv) = 1
      idmxyz(2,iv) = 1
      idmxyz(3,iv) = 1
      ipkxyz(iv)   = 1
      inixyz(1,iv) = 0
      inixyz(2,iv) = 0
      inixyz(3,iv) = 0
      iupxyz(1,iv) = 0
      iupxyz(2,iv) = 0
      ibdxyz(1,iv) = 0
      ibdxyz(2,iv) = 0
      ibdxyz(3,iv) = 0

    END DO

!***********************************************************************
!  All non-time-dependent 3d variables, ubar,vbar,wbar,ptbar,pbar,
!  rhostr,qvbar,zp, j1,j2,j3,wcont,kmh,kmv, rprntl, ppi, csndsq,
!  ptcumsrc, qcumsrc, radfrc, j3inv,aj3x,aj3y,aj3z,
!  ptbari,pbari,rhostri
!***********************************************************************

    DO iv = id_ubar,nxyz3d
      stgxyz(1,iv) = 0.5
      stgxyz(2,iv) = 0.5
      stgxyz(3,iv) = 0.5
      idmxyz(1,iv) = 1
      idmxyz(2,iv) = 1
      idmxyz(3,iv) = 1
      ipkxyz(  iv) = 1
      inixyz(1,iv) = 1
      inixyz(2,iv) = 0
      inixyz(3,iv) = 0
      iupxyz(1,iv) = 0
      iupxyz(2,iv) = 0
      ibdxyz(1,iv) = 0
      ibdxyz(2,iv) = 0
      ibdxyz(3,iv) = 0
    END DO

!***********************************************************************

    iv = id_ubar ! ubar

    stgxyz(1,iv) = 0.0
    idmxyz(1,iv) = 0
    inixyz(2,iv) = id_vbar

!***********************************************************************

    iv = id_vbar ! vbar

    stgxyz(2,iv) = 0.0
    idmxyz(2,iv) = 0
    inixyz(2,iv) = -id_ubar

!***********************************************************************

    iv = id_wbar ! wbar

    stgxyz(3,iv) = 0.0
    idmxyz(3,iv) = 0

!***********************************************************************

    iv = id_zp ! zp

    stgxyz(3,iv) = 0.0
    idmxyz(3,iv) = 0
    inixyz(1,iv) = 1

!***********************************************************************

    iv = id_j1 ! j1

    stgxyz(1,iv) = 0.0
    idmxyz(1,iv) = 0
    inixyz(1,iv) = 0

!***********************************************************************

    iv = id_j2 ! j2

    stgxyz(2,iv) = 0.0
    idmxyz(2,iv) = 0
    inixyz(1,iv) = 0

!***********************************************************************

    iv = id_j3 ! j3

    stgxyz(3,iv) = 0.0
    idmxyz(3,iv) = 0
    inixyz(1,iv) = 0

!***********************************************************************

    iv = id_wcont ! wcont

    stgxyz(3,iv) = 0.0
    idmxyz(3,iv) = 0
    inixyz(1,iv) = 0

!***********************************************************************

    iv = id_kmh ! kmh

    inixyz(1,iv) = 0

    iv = id_kmv ! kmv

    inixyz(1,iv) = 0

!***********************************************************************

    iv = id_rprntl ! rprntl

    inixyz(1,iv) = 0

!***********************************************************************

    iv = id_ptcumsrc ! ptcumsrc

    inixyz(1,iv) = 0

!***********************************************************************

    DO iv = id_qcumsrc,id_qcumsrc+4  ! qcumsrc(nx,ny,nz,5)
      inixyz(1,iv  ) = 0
    END DO

!***********************************************************************

    iv = id_radfrc ! radfrc

    inixyz(1,iv) = 0

!***********************************************************************

    iv = id_j3inv ! j3inv

    stgxyz(3,iv) = 0.0
    idmxyz(3,iv) = 0
    inixyz(1,iv) = 0

    iv = id_aj3x  ! aj3x

    stgxyz(1,iv) = 0.0
    stgxyz(3,iv) = 0.0
    idmxyz(1,iv) = 0
    idmxyz(3,iv) = 0
    inixyz(1,iv) = 0

    iv = id_aj3y  ! aj3y

    stgxyz(2,iv) = 0.0
    stgxyz(3,iv) = 0.0
    idmxyz(2,iv) = 0
    idmxyz(3,iv) = 0
    inixyz(1,iv) = 0

    iv = id_aj3y  ! aj3y

    inixyz(1,iv) = 0

!***********************************************************************

    iv = id_w0avg ! w0avg

    stgxyz(3,iv) = 0.0
    idmxyz(3,iv) = 0
    inixyz(1,iv) = 0

!***********************************************************************

    DO iv=1,nexbc3d
      stgexbc(1,iv) = 0.5
      stgexbc(2,iv) = 0.5
      stgexbc(3,iv) = 0.5
      idmexbc(1,iv) = 1
      idmexbc(2,iv) = 1
      idmexbc(3,iv) = 1
      ipkexbc(  iv) = 1
      iniexbc(1,iv) = 0
      iniexbc(2,iv) = 0
      iniexbc(3,iv) = 0
      iupexbc(1,iv) = 0
      iupexbc(2,iv) = 0
      ibdexbc(1,iv) = 0
      ibdexbc(2,iv) = 0
      ibdexbc(3,iv) = 0
    END DO

    iv = id_u0exb  ! u0exb

    stgexbc(1,iv) = 0.0
    idmexbc(1,iv) = 0

    iv = id_v0exb  ! v0exb

    stgexbc(2,iv) = 0.0
    idmexbc(2,iv) = 0

    iv = id_w0exb  ! w0exb

    stgexbc(3,iv) = 0.0
    idmexbc(3,iv) = 0

    iv = id_udtexb ! udtexb

    stgexbc(1,iv) = 0.0
    idmexbc(1,iv) = 0

    iv = id_vdtexb ! vdtexb

    stgexbc(2,iv) = 0.0
    idmexbc(2,iv) = 0

    iv = id_wdtexb ! wdtexb

    stgexbc(3,iv) = 0.0
    idmexbc(3,iv) = 0

  END IF

  RETURN
END SUBROUTINE stgrid
