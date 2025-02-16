!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE STRCNTS                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE strcnts(nx,ny,nz, iconst,nicnst, rconst,nrcnst )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Store the integer, real and logical constants into 1-D vector
!  arrays.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/27/1992
!
!  MODIFICATIONS:
!
!  UPDATES:   E.J. Adlerman
!             August 1995 for ARPS 4.0.22
!
!  09/04/1996 (Yuhe Liu)
!    Fixed bugs in assigning ARPS arrays of constants into the 1-D
!    vector array.
!
!  04/09/1997 (Yuhe Liu)
!    Rewrote this subroutine so that users may be easier to update
!    when there are new constants added to the model
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz,nicnst,nrcnst
  INTEGER :: iconst(nicnst)
  REAL :: rconst(nrcnst)

  INCLUDE 'agricst.inc'
  INCLUDE 'globcst.inc'     ! Global constants that control model
  INCLUDE 'bndry.inc'       ! Boundary condition control parameters
  INCLUDE 'exbc.inc'        ! EXBC control parameters
  INCLUDE 'grid.inc'

  INTEGER :: i,nc, niptr,nrptr,ncptr
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
!  Store the values of real constants into array rconst.
!  These are from the *.inc files, ordered as they appear in
!  an alphabetical order of those files.
!
!-----------------------------------------------------------------------
!
  nrptr = 0
!
!-----------------------------------------------------------------------
!
!  Real constants (everything common blocks) from globcst.inc and grid.inc 
!
!-----------------------------------------------------------------------
!
  rconst(nrptr+1) = ptpert0
  rconst(nrptr+2) = pt0radx
  rconst(nrptr+3) = pt0rady
  rconst(nrptr+4) = pt0radz
  rconst(nrptr+5) = pt0ctrx
  rconst(nrptr+6) = pt0ctry
  rconst(nrptr+7) = pt0ctrz
  rconst(nrptr+8) = ubar0
  rconst(nrptr+9) = vbar0

  nrptr = nrptr + 9

  rconst(nrptr+1) = hmount
  rconst(nrptr+2) = mntwidx
  rconst(nrptr+3) = mntwidy
  rconst(nrptr+4) = mntctrx
  rconst(nrptr+5) = mntctry

  nrptr = nrptr + 5

  rconst(nrptr+1) = dx
  rconst(nrptr+2) = dy
  rconst(nrptr+3) = dz
  rconst(nrptr+4) = dzmin
  rconst(nrptr+5) = zrefsfc
  rconst(nrptr+6) = dlayer1
  rconst(nrptr+7) = dlayer2
  rconst(nrptr+8) = strhtune
  rconst(nrptr+9) = zflat
  rconst(nrptr+10) = dxinv
  rconst(nrptr+11) = dyinv
  rconst(nrptr+12) = dzinv
  rconst(nrptr+13) = xl
  rconst(nrptr+14) = yl
  rconst(nrptr+15) = zh
  rconst(nrptr+16) = ctrlat
  rconst(nrptr+17) = ctrlon
  rconst(nrptr+18) = xorig
  rconst(nrptr+19) = yorig
  rconst(nrptr+20) = zorig
  rconst(nrptr+21) = xgrdorg
  rconst(nrptr+22) = ygrdorg

  nrptr = nrptr + 22

  rconst(nrptr+1) = csfactr
  rconst(nrptr+2) = csound

  nrptr = nrptr + 2

  rconst(nrptr+1) = dtbig
  rconst(nrptr+2) = dtsml
  rconst(nrptr+3) = tstart
  rconst(nrptr+4) = tstop
  rconst(nrptr+5) = curtim
  rconst(nrptr+6) = tacoef

  nrptr = nrptr + 6

  rconst(nrptr+1) = flteps

  nrptr = nrptr + 1

  rconst(nrptr+1) = alfcoef
  rconst(nrptr+2) = prantl
  rconst(nrptr+3) = tmixcst
  rconst(nrptr+4) = kmlimit

  nrptr = nrptr + 4

  rconst(nrptr+1) = cfcm2h
  rconst(nrptr+2) = cfcm2v
  rconst(nrptr+3) = cfcmh2
  rconst(nrptr+4) = cfcmv2
  rconst(nrptr+5) = cfcm4h
  rconst(nrptr+6) = cfcm4v
  rconst(nrptr+7) = cfcmh4
  rconst(nrptr+8) = cfcmv4
  rconst(nrptr+9) = scmixfctr

  nrptr = nrptr + 9

  rconst(nrptr+1) = zbrdmp
  rconst(nrptr+2) = cfrdmp
  rconst(nrptr+3) = cdvdmph
  rconst(nrptr+4) = cdvdmpv
  rconst(nrptr+5) = divdmpndh
  rconst(nrptr+6) = divdmpndv

  nrptr = nrptr + 6

  rconst(nrptr+1) = wcldbs
  rconst(nrptr+2) = confrq
  rconst(nrptr+3) = qpfgfrq
  rconst(nrptr+4) = kffbfct

  nrptr = nrptr + 4

  rconst(nrptr+1) = umove
  rconst(nrptr+2) = vmove

  nrptr = nrptr + 2

  rconst(nrptr+1) = cdmlnd
  rconst(nrptr+2) = cdmwtr
  rconst(nrptr+3) = cdhlnd
  rconst(nrptr+4) = cdhwtr
  rconst(nrptr+5) = cdqlnd
  rconst(nrptr+6) = cdqwtr
  rconst(nrptr+7) = pbldpth0
  rconst(nrptr+8) = lsclpbl0
  rconst(nrptr+9) = dtqflxdis

  nrptr = nrptr + 9

  rconst(nrptr+1) = dtsfc
  rconst(nrptr+2) = lai0
  rconst(nrptr+3) = roufns0
  rconst(nrptr+4) = veg0
  rconst(nrptr+5) = ptslnd0
  rconst(nrptr+6) = ptswtr0
  rconst(nrptr+7) = tsoil0
  rconst(nrptr+8) = wetsfc0
  rconst(nrptr+9) = wetdp0
  rconst(nrptr+10) = wetcanp0
  rconst(nrptr+11) = snowdpth0

  nrptr = nrptr + 11

  rconst(nrptr+1) = latitud
  rconst(nrptr+2) = longitud

  nrptr = nrptr + 2

  rconst(nrptr+1) = tfmtprt
  rconst(nrptr+2) = tmaxmin
  rconst(nrptr+3) = tenergy
  rconst(nrptr+4) = trstout
  rconst(nrptr+5) = timgdmp
  rconst(nrptr+6) = tplots

  nrptr = nrptr + 6

  rconst(nrptr+1) = chkdpth
  rconst(nrptr+2) = twindow
  rconst(nrptr+3) = tceltrk
  rconst(nrptr+4) = tcrestr

  nrptr = nrptr + 4

  rconst(nrptr+1) = trulat1
  rconst(nrptr+2) = trulat2
  rconst(nrptr+3) = trulon
  rconst(nrptr+4) = sclfct
  rconst(nrptr+5) = swlatu
  rconst(nrptr+6) = swlonu
  rconst(nrptr+7) = nelatu
  rconst(nrptr+8) = nelonu
  rconst(nrptr+9)  = swlatv
  rconst(nrptr+10) = swlonv
  rconst(nrptr+11) = nelatv
  rconst(nrptr+12) = nelonv
  rconst(nrptr+13) = swlats
  rconst(nrptr+14) = swlons
  rconst(nrptr+15) = nelats
  rconst(nrptr+16) = nelons

  nrptr = nrptr + 16

  rconst(nrptr+1) = tstrtdmp
  rconst(nrptr+2) = thisdmp

  nrptr = nrptr + 2

  DO i=1,hdmpmax
    rconst(nrptr+i) = hdmptim(i)
  END DO

  nrptr = nrptr + hdmpmax

  rconst(nrptr+1) = dtrad

  nrptr = nrptr + 1

  DO i=1,max_acct
    rconst(nrptr+i)  = acct_cpu_time(i)
  ENDDO
  nrptr = nrptr + max_acct 

  DO i=1,max_acct
    rconst(nrptr+i)  = acct_wall_time(i)
  ENDDO
  nrptr = nrptr + max_acct 

  rconst(nrptr+1)  = wall_init
  rconst(nrptr+2)  = cpu_init
  rconst(nrptr+3)  = total_cpu
  rconst(nrptr+4)  = total_wall
  rconst(nrptr+5)  = wall_acct
  rconst(nrptr+6)  = cpu_acct
  rconst(nrptr+7)  = wall_inter
  rconst(nrptr+8) = cpu_inter

  nrptr = nrptr + 8 
!
!-----------------------------------------------------------------------
!
!  Real constants from bndry.inc
!
!-----------------------------------------------------------------------
!
  rconst(nrptr+1) = c_phase
  rconst(nrptr+2) = rlxlbc

  nrptr = nrptr + 2
!
!-----------------------------------------------------------------------
!
!  Write out the number of real constants stored in rconst
!
!-----------------------------------------------------------------------
!
  IF ( verbose6 ) THEN
    WRITE(6,'(a,i6)')                                                   &
        ' The number of real constants stored in rconst is ', nrptr
    DO i=1,nrptr
      WRITE(6,'(a,i5,a,e20.10)') 'rconst(',i,')=',rconst(i)
    END DO
  END IF
!
!-----------------------------------------------------------------------
!
!  Set rest of rconst array to 0.0
!
!-----------------------------------------------------------------------
!
  DO i=nrptr+1,nrcnst
    rconst(i) = 0.0
  END DO

!
!************************************************************************
!
!  Store the values of integer constants into array iconst.
!
!************************************************************************
!
  niptr = 0
!
!-----------------------------------------------------------------------
!
!  Integer constants from (everything in common block) globcst.inc and 
!  grid.inc
!
!-----------------------------------------------------------------------
!
  iconst(niptr+1) = mgrid
  iconst(niptr+2) = nestgrd

  niptr = niptr + 2

  iconst(niptr+1) = lfnkey
  iconst(niptr+2) = nocmnt
  iconst(niptr+3) = runmod

  niptr = niptr + 3

  iconst(niptr+1) = initopt
  iconst(niptr+2) = inibasopt
  iconst(niptr+3) = viniopt
  iconst(niptr+4) = pt0opt
  iconst(niptr+5) = inifmt
  iconst(niptr+6) = inifunt
  iconst(niptr+7) = buoyopt
  iconst(niptr+8) = bsnesq
  iconst(niptr+9) = buoy2nd
  iconst(niptr+10) = rhofctopt

  niptr = niptr + 10

  iconst(niptr+1) = crdtrns
  iconst(niptr+2) = ternopt
  iconst(niptr+3) = mntopt
  iconst(niptr+4) = strhopt

  niptr = niptr + 4

  iconst(niptr+1) = sfunit

  niptr = niptr + 1

  iconst(niptr+1) = vimplct
  iconst(niptr+2) = ptsmlstp
  iconst(niptr+3) = peqopt
  iconst(niptr+4) = csopt

  niptr = niptr + 4

  iconst(niptr+1) = nstep
  iconst(niptr+2) = nsteps
  iconst(niptr+3) = nsmstp

  niptr = niptr + 3

  iconst(niptr+1) = year
  iconst(niptr+2) = month
  iconst(niptr+3) = day
  iconst(niptr+4) = hour
  iconst(niptr+5) = minute
  iconst(niptr+6) = second
  iconst(niptr+7) = jday

  niptr = niptr + 7

  iconst(niptr+1) = madvopt
  iconst(niptr+2) = sadvopt
  iconst(niptr+3) = fctorderopt
  iconst(niptr+4) = fctadvptprt

  niptr = niptr + 4

  iconst(niptr+1) = tmixopt
  iconst(niptr+2) = trbisotp
  iconst(niptr+3) = tkeopt
  iconst(niptr+4) = trbvimp
  iconst(niptr+5) = tmixvert

  niptr = niptr + 5

  iconst(niptr+1) = cmix2nd
  iconst(niptr+2) = cmix4th

  niptr = niptr + 2

  iconst(niptr+1) = raydmp
  iconst(niptr+2) = divdmp

  niptr = niptr + 2

  iconst(niptr+1) = mphyopt
  iconst(niptr+2) = moist
  iconst(niptr+3) = cnvctopt
  iconst(niptr+4) = ice
  iconst(niptr+5) = idownd

  niptr = niptr + 5

  iconst(niptr+1) = sfcphy
  iconst(niptr+2) = landwtr
  iconst(niptr+3) = cdhwtropt
  iconst(niptr+4) = wtrexist
  iconst(niptr+5) = pbldopt
  iconst(niptr+6) = sflxdis
  iconst(niptr+7) = sfcdiag
  iconst(niptr+8) = tqflxdis

  niptr = niptr + 8

  iconst(niptr+1) = sfcdat
  iconst(niptr+2) = soilinit
  iconst(niptr+3) = sfcunit
  iconst(niptr+4) = nsfcst
  iconst(niptr+5) = nstyp
  iconst(niptr+6) = styp
  iconst(niptr+7) = vtyp

  niptr = niptr + 7

  iconst(niptr+1 ) = coriopt
  iconst(niptr+2 ) = coriotrm

  niptr = niptr + 2

  iconst(niptr+1) = ldirnam
  iconst(niptr+2) = nfmtprt
  iconst(niptr+3) = nmaxmin
  iconst(niptr+4) = nenergy
  iconst(niptr+5) = nrstout
  iconst(niptr+6) = totout
  iconst(niptr+7) = grdout
  iconst(niptr+8) = basout
  iconst(niptr+9) = varout
  iconst(niptr+10) = mstout
  iconst(niptr+11) = rainout
  iconst(niptr+12) = iceout
  iconst(niptr+13) = trbout
  iconst(niptr+14) = tkeout
  iconst(niptr+15) = sfcout
  iconst(niptr+16) = exbcdmp
  iconst(niptr+17) = extdadmp
  iconst(niptr+18) = sfcdmp
  iconst(niptr+19) = soildmp
  iconst(niptr+20) = terndmp
  iconst(niptr+21) = exbcfmt
  iconst(niptr+22) = sfcfmt
  iconst(niptr+23) = soilfmt
  iconst(niptr+24) = ternfmt
  iconst(niptr+25) = incrfmt

  niptr = niptr + 25

  iconst(niptr+1) = prcout
  iconst(niptr+2) = radout
  iconst(niptr+3) = flxout
  niptr = niptr + 3

  iconst(niptr+1) = hdmpopt
  iconst(niptr+2) = hdmpfmt
  iconst(niptr+3) = grbpkbit
  iconst(niptr+4) = hdfcompr
  iconst(niptr+5) = nhisdmp
  iconst(niptr+6) = nstrtdmp
  iconst(niptr+7) = numhdmp

  niptr = niptr + 7

  DO i=1,numhdmp
    iconst(niptr+i) = hdmpstp(i)
  END DO
  niptr = niptr + numhdmp

  iconst(niptr+1) = nchdmp
  iconst(niptr+2) = ldmpf
  iconst(niptr+3) = imgopt
  iconst(niptr+4) = nimgdmp
  iconst(niptr+5) = pltopt
  iconst(niptr+6) = nplots
  iconst(niptr+7) = filcmprs
  iconst(niptr+8) = readyfl

  niptr = niptr + 8

  iconst(niptr+1) = cltkopt
  iconst(niptr+2) = grdtrns
  iconst(niptr+3) = nceltrk

  niptr = niptr + 3

  iconst(niptr+1) = lvldbg

  niptr = niptr + 1

  iconst(niptr+1) = restrt
  iconst(niptr+2) = rstiunt
  iconst(niptr+3) = rstount

  niptr = niptr + 3

  iconst(niptr+1) = nchmax
  iconst(niptr+2) = ncheng
  iconst(niptr+3) = grafh
  iconst(niptr+4) = gridh
  iconst(niptr+5) = dsindex
  iconst(niptr+6) = gridid

  niptr = niptr + 6

  iconst(niptr+1) = mapproj
  iconst(niptr+2) = mpfctopt
  iconst(niptr+3) = mptrmopt
  iconst(niptr+4) = maptest

  niptr = niptr + 4

  iconst(niptr+1) = radopt
  iconst(niptr+2) = radstgr
  iconst(niptr+3) = rlwopt
  iconst(niptr+4) = raddiag
  iconst(niptr+5) = nradstp
  iconst(niptr+6) = ict
  iconst(niptr+7) = icb

  niptr = niptr + 7
!
!-----------------------------------------------------------------------
!
!  Integer constants from bndry.inc
!
!-----------------------------------------------------------------------
!
  iconst(niptr+1 ) = lbcopt       ! from bndry.inc
  iconst(niptr+2 ) = ebc
  iconst(niptr+3 ) = wbc
  iconst(niptr+4 ) = nbc
  iconst(niptr+5 ) = sbc
  iconst(niptr+6 ) = tbc
  iconst(niptr+7 ) = bbc
  iconst(niptr+8 ) = rbcopt
  iconst(niptr+9 ) = pdetrnd

  niptr = niptr + 9
!
!-----------------------------------------------------------------------
!
!  Integer constants from exbc.inc
!
!-----------------------------------------------------------------------
!
  iconst(niptr+1 ) = abstinit       ! from exbc.inc
  iconst(niptr+2 ) = abststop
  iconst(niptr+3 ) = abstfcst0
  iconst(niptr+4 ) = abstfcst

  niptr = niptr + 4

  iconst(niptr+1 ) = ubcrd
  iconst(niptr+2 ) = vbcrd
  iconst(niptr+3 ) = wbcrd
  iconst(niptr+4 ) = ptbcrd
  iconst(niptr+5 ) = prbcrd
  iconst(niptr+6 ) = qvbcrd
  iconst(niptr+7 ) = qcbcrd
  iconst(niptr+8 ) = qrbcrd
  iconst(niptr+9 ) = qibcrd
  iconst(niptr+10) = qsbcrd
  iconst(niptr+11) = qhbcrd

  niptr = niptr + 11

  iconst(niptr+1)  = use_acct
  iconst(niptr+2)  = current_acct
  iconst(niptr+3)  = interrupt_acct
  niptr = niptr + 3 
!
!-----------------------------------------------------------------------
!
!  Write out the number of integer constants stored in iconst
!
!-----------------------------------------------------------------------
!
  IF ( verbose6 ) THEN
    WRITE (6,'(a,i6)')                                                  &
        'The number of integer constants stored in iconst is ',niptr
    DO i=1,niptr
      WRITE(6,'(a,i5,a,i10)') 'iconst(',i,')=',iconst(i)
    END DO
  END IF
!
!-----------------------------------------------------------------------
!
!  Encode the character strings defined in globcst.inc into integers
!  and then store them into iconst
!
!-----------------------------------------------------------------------
!
  ncptr = niptr

  DO i=1,20
    iconst(ncptr+i)  = ICHAR( arpsversion(i:i) )
  END DO

  ncptr = ncptr + 20

  DO i=1,80
    iconst(ncptr+i+0*80)  = ICHAR( runname (i:i) )
    iconst(ncptr+i+1*80)  = ICHAR( inifile (i:i) )
    iconst(ncptr+i+2*80)  = ICHAR( inigbf  (i:i) )
    iconst(ncptr+i+3*80)  = ICHAR( terndta (i:i) )
    iconst(ncptr+i+4*80)  = ICHAR( sndfile (i:i) )
    iconst(ncptr+i+5*80)  = ICHAR( sfcdtfl (i:i) )
    iconst(ncptr+i+6*80)  = ICHAR( soilinfl(i:i) )
    iconst(ncptr+i+7*80)  = ICHAR( dirname (i:i) )
    iconst(ncptr+i+8*80)  = ICHAR( hdmpfn  (i:i) )
    iconst(ncptr+i+9*80)  = ICHAR( rstinf  (i:i) )
    iconst(ncptr+i+10*80) = ICHAR( rstoutf (i:i) )
  END DO

  ncptr = ncptr + 11*80

  DO nc =1,MIN(nocmnt,50)
    DO i=1,80
      iconst(ncptr+i+(nc-1)*80) = ICHAR( cmnt(nc)(i:i) )
    END DO
  END DO

  DO nc = MIN(nocmnt,50)+1, 50
    DO i=1,80
      iconst(ncptr+i+(nc-1)*80) = ICHAR( ' ' )
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Write out the number of integers in encoding character constants
!  stored in iconst
!
!-----------------------------------------------------------------------
!
  IF ( verbose6 ) THEN
    WRITE (6,'(a,a,i6)')                                                &
        ' The number of integers in encoding cahracter constants',      &
        ' stored in iconst is ',ncptr-niptr
  END IF

  RETURN
END SUBROUTINE strcnts

!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE GETCNTS                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE getcnts(nx,ny,nz, iconst,nicnst, rconst,nrcnst)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Retrieve the integer, real and logical constants stored in
!  1-D vector arrays. These constants will be passed to the model
!  solver through common blocks.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!
!  MODIFICATIONS:
!
!  10/27/1992
!    UPDATED: E.J.Adlerman, May 1995, arps4.0.18
!
!  09/04/1996 (Yuhe Liu)
!    Fixed bugs in assigning ARPS arrays of constants into the 1-D
!    vector array.
!
!  04/09/1997 (Yuhe Liu)
!    Rewrote this subroutine so that users may be easier to update
!    when there are new constants added to the model
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz,nicnst,nrcnst

  INTEGER :: iconst(nicnst)
  REAL :: rconst(nrcnst)

  INCLUDE 'agricst.inc'

  INCLUDE 'phycst.inc'      ! Unchanging physical constants
  INCLUDE 'globcst.inc'     ! Global constants that control model
  INCLUDE 'bndry.inc'       ! Boundary condition control parameters
  INCLUDE 'exbc.inc'        ! EXBC control parameters
  INCLUDE 'sfcphycst.inc'   ! Unchanging physical constant
  INCLUDE 'cumucst.inc'     ! cumulus parameterization
  INCLUDE 'soilcst.inc'     ! soil-vegetation parameterization
  INCLUDE 'grid.inc'

  INTEGER :: i,nc, niptr,nrptr, ncptr
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
!  Retrieve the values of real constants stored in array rconst.
!
!-----------------------------------------------------------------------
!
  nrptr = 0
!
!-----------------------------------------------------------------------
!
!  Integer constants from globcst.inc
!
!-----------------------------------------------------------------------
!
  ptpert0   = rconst(nrptr+1)
  pt0radx   = rconst(nrptr+2)
  pt0rady   = rconst(nrptr+3)
  pt0radz   = rconst(nrptr+4)
  pt0ctrx   = rconst(nrptr+5)
  pt0ctry   = rconst(nrptr+6)
  pt0ctrz   = rconst(nrptr+7)
  ubar0     = rconst(nrptr+8)
  vbar0     = rconst(nrptr+9)

  nrptr = nrptr + 9

  hmount    = rconst(nrptr+1)
  mntwidx   = rconst(nrptr+2)
  mntwidy   = rconst(nrptr+3)
  mntctrx   = rconst(nrptr+4)
  mntctry   = rconst(nrptr+5)

  nrptr = nrptr + 5

  dx        = rconst(nrptr+1)
  dy        = rconst(nrptr+2)
  dz        = rconst(nrptr+3)
  dzmin     = rconst(nrptr+4)
  zrefsfc   = rconst(nrptr+5)
  dlayer1   = rconst(nrptr+6)
  dlayer2   = rconst(nrptr+7)
  strhtune  = rconst(nrptr+8)
  zflat     = rconst(nrptr+9)
  dxinv     = rconst(nrptr+10)
  dyinv     = rconst(nrptr+11)
  dzinv     = rconst(nrptr+12)
  xl        = rconst(nrptr+13)
  yl        = rconst(nrptr+14)
  zh        = rconst(nrptr+15)
  ctrlat    = rconst(nrptr+16)
  ctrlon    = rconst(nrptr+17)
  xorig     = rconst(nrptr+18)
  yorig     = rconst(nrptr+19)
  zorig     = rconst(nrptr+20)
  xgrdorg   = rconst(nrptr+21)
  ygrdorg   = rconst(nrptr+22)

  nrptr = nrptr + 22

  csfactr   = rconst(nrptr+1)
  csound    = rconst(nrptr+2)

  nrptr = nrptr + 2

  dtbig     = rconst(nrptr+1)
  dtsml     = rconst(nrptr+2)
  tstart    = rconst(nrptr+3)
  tstop     = rconst(nrptr+4)
  curtim    = rconst(nrptr+5)
  tacoef    = rconst(nrptr+6)

  nrptr = nrptr + 6

  flteps    = rconst(nrptr+1)

  nrptr = nrptr + 1

  alfcoef   = rconst(nrptr+1)
  prantl    = rconst(nrptr+2)
  tmixcst   = rconst(nrptr+3)
  kmlimit   = rconst(nrptr+4)

  nrptr = nrptr + 4

  cfcm2h    = rconst(nrptr+1)
  cfcm2v    = rconst(nrptr+2)
  cfcmh2    = rconst(nrptr+3)
  cfcmv2    = rconst(nrptr+4)
  cfcm4h    = rconst(nrptr+5)
  cfcm4v    = rconst(nrptr+6)
  cfcmh4    = rconst(nrptr+7)
  cfcmv4    = rconst(nrptr+8)
  scmixfctr = rconst(nrptr+9)

  nrptr = nrptr + 9

  zbrdmp    = rconst(nrptr+1)
  cfrdmp    = rconst(nrptr+2)
  cdvdmph   = rconst(nrptr+3)
  cdvdmpv   = rconst(nrptr+4)
  divdmpndh = rconst(nrptr+5)
  divdmpndv = rconst(nrptr+6)

  nrptr = nrptr + 6

  wcldbs    = rconst(nrptr+1)
  confrq    = rconst(nrptr+2)
  qpfgfrq   = rconst(nrptr+3)
  kffbfct   = rconst(nrptr+4)

  nrptr = nrptr + 4

  umove     = rconst(nrptr+1)
  vmove     = rconst(nrptr+2)

  nrptr = nrptr + 2

  cdmlnd    = rconst(nrptr+1)
  cdmwtr    = rconst(nrptr+2)
  cdhlnd    = rconst(nrptr+3)
  cdhwtr    = rconst(nrptr+4)
  cdqlnd    = rconst(nrptr+5)
  cdqwtr    = rconst(nrptr+6)
  pbldpth0  = rconst(nrptr+7)
  lsclpbl0  = rconst(nrptr+8)
  dtqflxdis = rconst(nrptr+9)

  nrptr = nrptr + 9

  dtsfc     = rconst(nrptr+1)
  lai0      = rconst(nrptr+2)
  roufns0   = rconst(nrptr+3)
  veg0      = rconst(nrptr+4)
  ptslnd0   = rconst(nrptr+5)
  ptswtr0   = rconst(nrptr+6)
  tsoil0    = rconst(nrptr+7)
  wetsfc0   = rconst(nrptr+8)
  wetdp0    = rconst(nrptr+9)
  wetcanp0  = rconst(nrptr+10)
  snowdpth0  = rconst(nrptr+11)

  nrptr = nrptr + 11

  latitud   = rconst(nrptr+1)
  longitud  = rconst(nrptr+2)

  nrptr = nrptr + 2

  tfmtprt   = rconst(nrptr+1)
  tmaxmin   = rconst(nrptr+2)
  tenergy   = rconst(nrptr+3)
  trstout   = rconst(nrptr+4)
  timgdmp   = rconst(nrptr+5)
  tplots    = rconst(nrptr+6)

  nrptr = nrptr + 6

  chkdpth   = rconst(nrptr+1)
  twindow   = rconst(nrptr+2)
  tceltrk   = rconst(nrptr+3)
  tcrestr   = rconst(nrptr+4)

  nrptr = nrptr + 4

  trulat1   = rconst(nrptr+1)
  trulat2   = rconst(nrptr+2)
  trulon    = rconst(nrptr+3)
  sclfct    = rconst(nrptr+4)
  swlatu    = rconst(nrptr+5)
  swlonu    = rconst(nrptr+6)
  nelatu    = rconst(nrptr+7)
  nelonu    = rconst(nrptr+8)
  swlatv    = rconst(nrptr+9)
  swlonv    = rconst(nrptr+10)
  nelatv    = rconst(nrptr+11)
  nelonv    = rconst(nrptr+12)
  swlats    = rconst(nrptr+13)
  swlons    = rconst(nrptr+14)
  nelats    = rconst(nrptr+15)
  nelons    = rconst(nrptr+16)

  nrptr = nrptr + 16

  tstrtdmp  = rconst(nrptr+1)
  thisdmp   = rconst(nrptr+2)

  nrptr = nrptr + 2

  DO i=1,hdmpmax
    hdmptim(i) = rconst(nrptr+i)
  END DO
  nrptr = nrptr + hdmpmax

  dtrad     = rconst(nrptr+1)

  nrptr = nrptr + 1

  DO i=1,max_acct
    acct_cpu_time(i) = rconst(nrptr+i)
  ENDDO
  nrptr = nrptr + max_acct

  DO i=1,max_acct
    acct_wall_time(i) = rconst(nrptr+i)
  ENDDO
  nrptr = nrptr + max_acct

  wall_init        = rconst(nrptr+1)
  cpu_init         = rconst(nrptr+2)
  total_cpu        = rconst(nrptr+3)
  total_wall       = rconst(nrptr+4)
  wall_acct        = rconst(nrptr+5)
  cpu_acct         = rconst(nrptr+6)
  wall_inter       = rconst(nrptr+7)
  cpu_inter        = rconst(nrptr+8)

  nrptr = nrptr + 8 
!
!-----------------------------------------------------------------------
!
!  Real constants from bndry.inc
!
!-----------------------------------------------------------------------
!
  c_phase   = rconst(nrptr+1)
  rlxlbc    = rconst(nrptr+2)

  nrptr = nrptr + 2
!
!-----------------------------------------------------------------------
!
!  Write out the number of real constants retrieved from rconst
!
!-----------------------------------------------------------------------
!
  IF ( verbose6 ) THEN
    WRITE(6,'(a,i3,a,i6)')                                              &
        ' The number of real constants in getcnst for grid ',mgrid,     &
        ' is ', nrptr
    DO i=1,nrptr
      WRITE(6,'(a,i5,a,e20.10)') 'rconst(',i,')=',rconst(i)
    END DO
  END IF
!
!***********************************************************************
!
!  Retrieve the values of integer constants stored in array iconst.
!
!***********************************************************************
!
  niptr = 0
!
!-----------------------------------------------------------------------
!
!  Integer constants from globcst.inc
!
!-----------------------------------------------------------------------
!
  mgrid     = iconst(niptr+1)
  nestgrd   = iconst(niptr+2)

  niptr = niptr + 2

  lfnkey    = iconst(niptr+1)
  nocmnt    = iconst(niptr+2)
  runmod    = iconst(niptr+3)

  niptr = niptr + 3

  initopt   = iconst(niptr+1)
  inibasopt = iconst(niptr+2)
  viniopt   = iconst(niptr+3)
  pt0opt    = iconst(niptr+4)
  inifmt    = iconst(niptr+5)
  inifunt   = iconst(niptr+6)
  buoyopt   = iconst(niptr+7)
  bsnesq    = iconst(niptr+8)
  buoyopt   = iconst(niptr+9)
  rhofctopt = iconst(niptr+10)

  niptr = niptr + 10

  crdtrns   = iconst(niptr+1)
  ternopt   = iconst(niptr+2)
  mntopt    = iconst(niptr+3)
  strhopt   = iconst(niptr+4)

  niptr = niptr + 4

  sfunit    = iconst(niptr+1)

  niptr = niptr + 1

  vimplct   = iconst(niptr+1)
  ptsmlstp  = iconst(niptr+2)
  peqopt    = iconst(niptr+3)
  csopt     = iconst(niptr+4)

  niptr = niptr + 4

  nstep     = iconst(niptr+1)
  nsteps    = iconst(niptr+2)
  nsmstp    = iconst(niptr+3)

  niptr = niptr + 3

  year      = iconst(niptr+1)
  month     = iconst(niptr+2)
  day       = iconst(niptr+3)
  hour      = iconst(niptr+4)
  minute    = iconst(niptr+5)
  second    = iconst(niptr+6)
  jday      = iconst(niptr+7)

  niptr = niptr + 7

  madvopt   = iconst(niptr+1)
  sadvopt   = iconst(niptr+2)
  fctorderopt = iconst(niptr+3)
  fctadvptprt = iconst(niptr+4)

  niptr = niptr + 4

  tmixopt   = iconst(niptr+1)
  trbisotp  = iconst(niptr+2)
  tkeopt    = iconst(niptr+3)
  trbvimp   = iconst(niptr+4)
  tmixvert  = iconst(niptr+5)

  niptr = niptr + 5

  cmix2nd   = iconst(niptr+1)
  cmix4th   = iconst(niptr+2)

  niptr = niptr + 2

  raydmp    = iconst(niptr+1)
  divdmp    = iconst(niptr+2)

  niptr = niptr + 2

  mphyopt   = iconst(niptr+1)
  moist     = iconst(niptr+2)
  cnvctopt  = iconst(niptr+3)
  ice       = iconst(niptr+4)
  idownd    = iconst(niptr+5)

  niptr = niptr + 5

  sfcphy    = iconst(niptr+1)
  landwtr   = iconst(niptr+2)
  cdhwtropt = iconst(niptr+3)
  wtrexist  = iconst(niptr+4)
  pbldopt   = iconst(niptr+5)
  sflxdis   = iconst(niptr+6)
  sfcdiag   = iconst(niptr+7)
  tqflxdis  = iconst(niptr+8)

  niptr = niptr + 8

  sfcdat    = iconst(niptr+1)
  soilinit  = iconst(niptr+2)
  sfcunit   = iconst(niptr+3)
  nsfcst    = iconst(niptr+4)
  nstyp     = iconst(niptr+5)
  styp      = iconst(niptr+6)
  vtyp      = iconst(niptr+7)

  niptr = niptr + 7

  coriopt   = iconst(niptr+1)
  coriotrm  = iconst(niptr+2)

  niptr = niptr + 2

  ldirnam   = iconst(niptr+1)
  nfmtprt   = iconst(niptr+2)
  nmaxmin   = iconst(niptr+3)
  nenergy   = iconst(niptr+4)
  nrstout   = iconst(niptr+5)
  totout    = iconst(niptr+6)
  grdout    = iconst(niptr+7)
  basout    = iconst(niptr+8)
  varout    = iconst(niptr+9)
  mstout    = iconst(niptr+10)
  rainout   = iconst(niptr+11)
  iceout    = iconst(niptr+12)
  trbout    = iconst(niptr+13)
  tkeout    = iconst(niptr+14)
  sfcout    = iconst(niptr+15)
  exbcdmp   = iconst(niptr+16)
  extdadmp  = iconst(niptr+17)
  sfcdmp    = iconst(niptr+18)
  soildmp   = iconst(niptr+19)
  terndmp   = iconst(niptr+20)
  exbcfmt   = iconst(niptr+21)
  sfcfmt    = iconst(niptr+22)
  soilfmt   = iconst(niptr+23)
  ternfmt   = iconst(niptr+24)
  incrfmt   = iconst(niptr+25)

  niptr = niptr + 25

  prcout  = iconst(niptr+1)
  radout  = iconst(niptr+2)
  flxout  = iconst(niptr+3)
  niptr = niptr + 3

  hdmpopt   = iconst(niptr+1)
  hdmpfmt   = iconst(niptr+2)
  grbpkbit  = iconst(niptr+3)
  hdfcompr  = iconst(niptr+4)
  nhisdmp   = iconst(niptr+5)
  nstrtdmp  = iconst(niptr+6)
  numhdmp   = iconst(niptr+7)

  niptr = niptr + 7

  DO i=1,numhdmp
    hdmpstp(i) = iconst(niptr+i)
  END DO
  niptr = niptr + numhdmp

  nchdmp    = iconst(niptr+1)
  ldmpf     = iconst(niptr+2)
  imgopt    = iconst(niptr+3)
  nimgdmp   = iconst(niptr+4)
  pltopt    = iconst(niptr+5)
  nplots    = iconst(niptr+6)
  filcmprs  = iconst(niptr+7)
  readyfl   = iconst(niptr+8)

  niptr = niptr + 8

  cltkopt   = iconst(niptr+1)
  grdtrns   = iconst(niptr+2)
  nceltrk   = iconst(niptr+3)

  niptr = niptr + 3

  lvldbg    = iconst(niptr+1)

  niptr = niptr + 1

  restrt    = iconst(niptr+1)
  rstiunt   = iconst(niptr+2)
  rstount   = iconst(niptr+3)

  niptr = niptr + 3

  nchmax    = iconst(niptr+1)
  ncheng    = iconst(niptr+2)
  grafh     = iconst(niptr+3)
  gridh     = iconst(niptr+4)
  dsindex   = iconst(niptr+5)
  gridid    = iconst(niptr+6)

  niptr = niptr + 6

  mapproj   = iconst(niptr+1)
  mpfctopt  = iconst(niptr+2)
  mptrmopt  = iconst(niptr+3)
  maptest   = iconst(niptr+4)

  niptr = niptr + 4

  radopt    = iconst(niptr+1)
  radstgr   = iconst(niptr+2)
  rlwopt    = iconst(niptr+3)
  raddiag   = iconst(niptr+4)
  nradstp   = iconst(niptr+5)
  ict       = iconst(niptr+6)
  icb       = iconst(niptr+7)

  niptr = niptr + 7
!
!-----------------------------------------------------------------------
!
!  Integer constants from bndry.inc
!
!-----------------------------------------------------------------------
!
  lbcopt    = iconst(niptr+1)
  ebc       = iconst(niptr+2)
  wbc       = iconst(niptr+3)
  nbc       = iconst(niptr+4)
  sbc       = iconst(niptr+5)
  tbc       = iconst(niptr+6)
  bbc       = iconst(niptr+7)
  rbcopt    = iconst(niptr+8)
  pdetrnd   = iconst(niptr+9)

  niptr = niptr + 9
!
!-----------------------------------------------------------------------
!
!  Integer constants from exbc.inc
!
!-----------------------------------------------------------------------
!
  abstinit  = iconst(niptr+1)
  abststop  = iconst(niptr+2)
  abstfcst0 = iconst(niptr+3)
  abstfcst  = iconst(niptr+4)

  niptr = niptr + 4

  ubcrd  = iconst(niptr+1)
  vbcrd  = iconst(niptr+2)
  wbcrd  = iconst(niptr+3)
  ptbcrd = iconst(niptr+4)
  prbcrd = iconst(niptr+5)
  qvbcrd = iconst(niptr+6)
  qcbcrd = iconst(niptr+7)
  qrbcrd = iconst(niptr+8)
  qibcrd = iconst(niptr+9)
  qsbcrd = iconst(niptr+10)
  qhbcrd = iconst(niptr+11)

  niptr = niptr + 11

  use_acct       = iconst(niptr+1)
  current_acct   = iconst(niptr+2)
  interrupt_acct = iconst(niptr+3)
  niptr = niptr + 3
!
!-----------------------------------------------------------------------
!
!  Write out the number of integer constants stored in iconst
!
!-----------------------------------------------------------------------
!
  IF ( verbose6 ) THEN
    WRITE(6,'(a,i3,a,i6)')                                              &
        ' The number of integer constants retrieved for grid ',mgrid,   &
        ' is ', niptr
    DO i=1,niptr
      WRITE(6,'(a,i5,a,i10)') 'iconst(',i,')=',iconst(i)
    END DO
  END IF
!
!-----------------------------------------------------------------------
!
!  Decode the character strings stored as integers
!
!-----------------------------------------------------------------------
!
  ncptr = niptr

  DO i=1,20
    arpsversion(i:i) = CHAR( iconst(ncptr+i) )
  END DO
  ncptr = ncptr + 20

  DO i=1,80
    runname (i:i) = CHAR( iconst(ncptr+i+0*80) )
    inifile (i:i) = CHAR( iconst(ncptr+i+1*80) )
    inigbf  (i:i) = CHAR( iconst(ncptr+i+2*80) )
    terndta (i:i) = CHAR( iconst(ncptr+i+3*80) )
    sndfile (i:i) = CHAR( iconst(ncptr+i+4*80) )
    sfcdtfl (i:i) = CHAR( iconst(ncptr+i+5*80) )
    soilinfl(i:i) = CHAR( iconst(ncptr+i+6*80) )
    dirname (i:i) = CHAR( iconst(ncptr+i+7*80) )
    hdmpfn  (i:i) = CHAR( iconst(ncptr+i+8*80) )
    rstinf  (i:i) = CHAR( iconst(ncptr+i+9*80) )
    rstoutf (i:i) = CHAR( iconst(ncptr+i+10*80) )
  END DO

  ncptr = ncptr + 11*80

  DO nc =1,50
    DO i=1,80
      cmnt(nc)(i:i) = CHAR( iconst(ncptr+i+(nc-1)*80) )
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Write out the total number of integers encoding character
!  constants stored in iconst
!
!-----------------------------------------------------------------------
!
  IF ( verbose6 ) THEN
    WRITE(6,'(a,a,i3,a,i6)')                                            &
        ' The number of integers in encoding character constants',      &
        ' retrieved for grid ',mgrid, ' is ', ncptr-niptr
  END IF

  RETURN
END SUBROUTINE getcnts
