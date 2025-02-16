!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE RADTRNS                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE radtrns(nx,ny,nz,rbufsz,                                     &
           ptprt,pprt,qv,qscalar,                                       &
           ptbar,pbar,ppi, rhostr, tsfc,                                &
           x,y,z,zp, j3inv,                                             &
           radfrc, radsw,rnflx,radswnet,radlwin, cosss,                 &
           rsirbm,rsirdf,rsuvbm,rsuvdf, cosz,                           &
           fdirir,fdifir,fdirpar,fdifpar, st4,                          &
           plinv,tinv,qvinv,o3a,ccld,                                   &
           flxir,flcir,flxuv,flcuv,dfdts,                               &
           tauir,taual, tauswi,tauswl,reffi,reffl,                      &
           radbuf,sh, tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  This subroutine is to compute the atmospheric radiation forcing
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  03/11/1996
!
!  MODIFICATION HISTORY:
!
!  04/09/1997 (Yuhe Liu)
!  Removed dimension size check statement. For normal ARPS run, the
!  dimension size check is now done in the main driver, arps##. For
!  nested runs, no need to check the dimension size because the
!  working array is allocated automatically from existing space.
!
!  Added call of subroutine SETRADWRK to set the indeces of working
!  arrays in the radiation buffer. Those indeces are passed by
!  common blocks in radcst.inc.
!
!  10/11/1998 (Keith Brewster)
!  Added option for using RH in cloud optical depth calculation.
!
!  10/30/1998 (Keith Brewster)
!  Added calculation of aerosol density.
!
!  11/18/98 (Keith Brewster)
!  Changed pibar to ppi (full pi).
!
!  07/16/03  (J. Brotzge)
!  Fixed new radswnet bug (detected by T. Katopodes)
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'radcst.inc'

  INCLUDE 'bndry.inc'                !  add by Yunheng
  INCLUDE 'mp.inc'

  INTEGER :: nx,ny,nz
  INTEGER :: rbufsz
!
!-----------------------------------------------------------------------
!
!  Define ARPS variables
!
!-----------------------------------------------------------------------
!
  REAL :: x(nx)
  REAL :: y(ny)
  REAL :: z(nz)
  REAL :: zp(nx,ny,nz)

  REAL :: ptprt (nx,ny,nz)
  REAL :: pprt  (nx,ny,nz)
  REAL :: qv    (nx,ny,nz)

  REAL :: qscalar (nx,ny,nz,nscalar)

  REAL :: ptbar (nx,ny,nz)
  REAL :: pbar  (nx,ny,nz)
  REAL :: ppi   (nx,ny,nz)

  REAL :: rhostr(nx,ny,nz)
  REAL :: j3inv (nx,ny,nz)  ! Inverse of j3

  REAL :: tsfc  (nx,ny)     ! Surface temperature (K)

  REAL :: radfrc(nx,ny,nz)  ! Radiation forcing (K/s)

  REAL :: radsw (nx,ny)     ! Solar radiation reaching the surface
  REAL :: rnflx (nx,ny)     ! Net upward radiation flux
  REAL :: radswnet(nx,ny)   ! Net solar radiation at surface
  REAL :: radlwin(nx,ny)    ! Incoming longwave radiation at surface
  REAL :: flxd(nx,ny)

  REAL :: cosz  (nx,ny)     ! Cosine of zenith
  REAL :: cosss (nx,ny)     ! Cosine of angle between sun light and
                            ! surface terrain slope

  REAL :: sh(nx,ny)         ! augustin add sh
!
!-----------------------------------------------------------------------
!
!  Define 2-D variables for radiation calculation.
!
!-----------------------------------------------------------------------
!
  REAL :: rsirbm(nx,ny)     ! Solar IR surface albedo for beam radiation
  REAL :: rsirdf(nx,ny)     ! Solar IR surface albedo for diffuse radiation
  REAL :: rsuvbm(nx,ny)     ! Solar UV surface albedo for beam radiation
  REAL :: rsuvdf(nx,ny)     ! Solar UV surface albedo for diffuse radiation

  REAL :: fdirir (nx,ny)    ! all-sky direct downward IR flux
                            ! (0.7-10 micron) at the surface
  REAL :: fdifir (nx,ny)    ! all-sky diffuse downward IR flux
                            ! at the surface
  REAL :: fdirpar(nx,ny)    ! all-sky direct downward par flux
                            ! (0.4-0.7 micron) at the surface
  REAL :: fdifpar(nx,ny)    ! all-sky diffuse downward par flux
                            ! at the surface

  REAL :: st4(nx,ny)        ! Emission by the surface
!
!-----------------------------------------------------------------------
!
!  Arrays which have the vertical coordinate inversed, that
!  is, k=1 is for top while k=nz is at the surface.
!
!-----------------------------------------------------------------------
!
  REAL :: plinv (nx,ny,nz)  ! Pressure in mb at scalar points
  REAL :: tinv  (nx,ny,nz)  ! Temperature
  REAL :: qvinv (nx,ny,nz)  ! Water vapor mixing ratio (g/g)

  REAL :: o3a   (nx,ny,nz)  ! Ozone (o3) mixing ratio (g/g)
  REAL :: ccld  (nx,ny,nz)  ! Cloud coverage (fraction)

  REAL :: flxir (nx,ny,nz)  ! all-sky net downward flux
  REAL :: flcir (nx,ny,nz)  ! clear-sky net downward flux

  REAL :: flxuv (nx,ny,nz)  ! all-sky solar flux (downward minus upward)
  REAL :: flcuv (nx,ny,nz)  ! clear-sky solar flux (downward minus upward)

  REAL :: dfdts (nx,ny,nz)  ! Sensitivity of net downward flux to surface
                            ! temperature

  REAL :: tauir (nx,ny,nz)  ! Cloud optical depth for LW IR
  REAL :: taual (nx,ny,nz)  ! Aerosol optical thickness

  REAL :: tauswi(nx,ny,nz)  ! Cloud optical depth for solar IR for
                            ! ice particles
  REAL :: tauswl(nx,ny,nz)  ! Cloud optical depth for solar IR for
                            ! liquid particles

  REAL :: reffi (nx,ny,nz)  ! Effective cloud-particle size for
                            ! ice particles
  REAL :: reffl (nx,ny,nz)  ! Effective cloud-particle size for
                            ! liquid particles

  REAL :: tem1 (nx,ny,nz)   ! Work array for message passing
!
!-----------------------------------------------------------------------
!
!  Include file radcst.inc which contains the definition of dimension
!  sizes for 2-d and 3-d temporary arrays represented by a buffer,
!  radbuf. When radopt is NOT set to 2, the dimensions and buffer
!  sizes can be 1. Otherwise, the dimension sizes should be the same
!  as nx, ny, and nz, and the buffer size should be larger than the
!  total size of 27 2-d arrays and 44 3-d arrays.
!
!  integer n2d_rad  ! number of 2-d arrays in the buffer
!  integer n3d_rad  ! number of 3-d arrays in the buffer
!
!  integer rbufsz   ! nx*ny*(n2d_rad+n3d_rad*nz)
!  real radbuf( rbufsz )
!
!  The 2-d arrays should be always at the beginning of radbuf and
!  the 3-d arrays then following.
!
!-----------------------------------------------------------------------
!
  REAL :: radbuf( rbufsz )
!
!-----------------------------------------------------------------------
!
!  Local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: nxy, nxyz
  INTEGER :: i,j,k
  INTEGER :: nq
  INTEGER :: im,jm,km, ij, ijm

  INTEGER :: istgr, jstgr
  INTEGER :: nxodd, nyodd, jeven, jodd

  INTEGER :: m, n

  INTEGER :: night           ! Flag for night time

  REAL :: tqe, dp, coolrate,heatrate
  REAL :: pk,psfc,tk,qvsat,rh,hgtagl,ccld1,ccld2,aden,psqc,aersum

  LOGICAL :: high

  INTEGER, PARAMETER :: rh2cldopt = 0
  REAL,    PARAMETER :: rhcldwgt  = 0.667
  REAL    :: qcwgt
!
!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------
!
  REAL :: f_qvsat,rh_to_cldcv,aeroden

!fpp$ expand (f_qvsat)
!!dir$ inline always f_qvsat
!*$*  inline routine (f_qvsat)

!
!-----------------------------------------------------------------------
!
!  Define additional layer to the top of model atmosphere
!
!-----------------------------------------------------------------------
!
  REAL, PARAMETER :: padd =  1.0  ! Additional pressure level (mb) to TOA
  REAL, PARAMETER :: tadd =223.0  ! Temperature for additional layer to TOA
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF ( rlwopt == 0 ) THEN
    high = .false.
  ELSE
    high = .true.
  END IF

  CALL setradwrk(nx,ny,nz)

  nxy  = nx*ny
  nxyz = nxy*nz
  qcwgt=1.-rhcldwgt

  IF ( (pbar(1,1,nz-1)*0.01) <= padd ) THEN
    WRITE (6,'(a/a)')                                                   &
        'The pressure at the top of atmosphere was too low to add ',    &
        'additional levels. Check the sounding profile.',               &
        'Program stopped in RADTRNS.'
    CALL arpsstop('arpsstop stopped RADTRNS pressure to low at the top ',1)
  END IF

  IF (P_QC > 0) THEN
    DO j=1,ny-1
      DO i=1,nx-1
        DO k=1,nz-1
          o3a(i,j,k)=qscalar(i,j,k,P_QC)
        END DO
      END DO
    END DO
  ELSE
    o3a(:,:,:) = 0.0
  END IF
!
!-----------------------------------------------------------------------
!
!  In the staggering case, only those points of (i,j) = (even,even)
!  and (odd,odd) will be calculated. Do defragment to move the points
!  to the left half domain, i=1,nx/2 and j=1,ny-1
!
!-----------------------------------------------------------------------
!
  nyodd = MOD(ny,2)
  nxodd = MOD(nx,2)

  IF ( radstgr == 0 ) THEN
    m = nx-1
    n = ny-1
!    WARNING: this loop is repeated for radstgr=1, changes made here must
!    also be made there as well.
    DO km=3,nz-1
      k = nz+1-km
      DO j=1,n
        DO i=1,m
          pk = pbar(i,j,k)+pprt(i,j,k)
          plinv(i,j,km) = 0.005*( ( pbar(i,j,k  ) + pprt(i,j,k  ) )     &
                                 +( pbar(i,j,k+1) + pprt(i,j,k+1) ) )
          tk = (ptbar(i,j,k)+ptprt(i,j,k))*ppi(i,j,k)
          tinv (i,j,km) = MAX(tk, 190.)
          qvinv(i,j,km) = MAX(qv(i,j,k), 1.0E-6)

          tqe = 0.0
          DO nq = 1,nscalarq
            tqe = tqe + qscalar(i,j,k,nq)
          END DO
          tqe = MAX(1.0E-20, tqe)

          IF( rh2cldopt > 0 ) THEN
            qvsat=MAX(f_qvsat(pk,tk),1.0E-20)
            rh = MIN(1.0,(qv(i,j,k)/qvsat))
            hgtagl = (0.5*(zp(i,j,k)+zp(i,j,k+1)))-zp(i,j,2)
            ccld1 = rh_to_cldcv(rh,hgtagl)
            ccld2 = MIN(1.0, MAX(0.0, 0.25*ALOG10(tqe)+1.5))
            ccld1 = MIN(1.0, MAX(ccld2, ccld1))
            ccld (i,j,km) = rhcldwgt*ccld1 + qcwgt*ccld2
!  borrow o3a to store psuedo-cloud
            psqc= 0.10*ccld1*qvsat
            IF (P_QC > 0) THEN
              o3a(i,j,k) = MAX(qscalar(i,j,k,P_QC),psqc)
            ELSE
              o3a(i,j,k) = MAX(0.0,psqc)
            END IF
          ELSE
            ccld(i,j,km) = MIN(1.0, MAX(0.0, 0.25*ALOG10(tqe)+1.5))
!  borrow o3a to store psuedo-cloud
            IF (P_QC > 0) THEN
              o3a(i,j,k) = qscalar(i,j,k,P_QC)
            ELSE
              o3a(i,j,k) = 0.0
            END IF
          END IF

          psfc = 0.5* (pbar(i,j,1)+pprt(i,j,1) +                        &
                       pbar(i,j,2)+pprt(i,j,2) )
          aden=2.0*aeroden(pk,psfc)
          taual(i,j,km)=aden*0.5*((pbar(i,j,k-1)+pprt(i,j,k-1))         &
                                 -(pbar(i,j,k+1)+pprt(i,j,k+1)))
        END DO
      END DO
    END DO

!    WARNING: this loop is repeated for radstgr=1, changes made here must
!    also be made there as well.
    DO j=1,n
      DO i=1,m
        ij = nx*(j-1) + i
        radbuf(ij) = tsfc(i,j)

        plinv(i,j,1) = padd
        tinv (i,j,1) = tadd
        qvinv(i,j,1) = 1.0E-6
        ccld (i,j,1) = 0.0
        taual(i,j,1) = 0.0

        plinv(i,j,2) = 0.5*(plinv(i,j,1)+plinv(i,j,3))
        tinv (i,j,2) = MAX(0.5*(tinv(i,j,1)+tinv(i,j,3)), 190.0)
        qvinv(i,j,2) = MAX(qv(i,j,nz-1), 1.0E-6)
        ccld (i,j,2) = 0.0
        taual(i,j,2) = 0.0

        plinv(i,j,nz) = 0.005*((pbar(i,j,1) + pprt(i,j,1))              &
                              +(pbar(i,j,2) + pprt(i,j,2)))
      END DO
    END DO

  ELSE IF ( radstgr == 1 ) THEN

    istgr = 1
    jstgr = 2
    m = nx-1
    n = ny/2
!    WARNING: this loop is repeated above for radstgr=0, changes made here
!    must also be made there as well.
    DO km=3,nz-1
      k = nz+1-km
      DO j=1,ny-1
        jodd = istgr*MOD(j,2)
        DO i=jstgr-jodd, nx-1, jstgr
          im = i
          jm = j/jstgr + istgr*MOD(i,2)

          pk = pbar(i,j,k)+pprt(i,j,k)
          plinv(im,jm,km) = 0.005*((pbar(i,j,k  ) + pprt(i,j,k  ))      &
                                  +(pbar(i,j,k+1) + pprt(i,j,k+1)))
          tk = (ptbar(i,j,k)+ptprt(i,j,k))*ppi(i,j,k)
          tinv (im,jm,km) = MAX(tk, 190.)
          qvinv(im,jm,km) = MAX(qv(i,j,k), 1.0E-6)

          tqe = 0.0
          DO nq = 1,nscalarq
            tqe = tqe + qscalar(i,j,k,nq)
          END DO
          tqe = MAX(1.0E-20, tqe)

          IF( rh2cldopt > 0 ) THEN
            qvsat=MAX(f_qvsat(pk,tk),1.0E-20)
            rh = MIN(1.0,(qv(i,j,k)/qvsat))
            hgtagl = (0.5*(zp(i,j,k)+zp(i,j,k+1)))-zp(i,j,2)
            ccld1 = rh_to_cldcv(rh,hgtagl)
            ccld2 = MIN(1.0, MAX(0.0, 0.25*ALOG10(tqe)+1.5))
            ccld1 = MIN(1.0, MAX(ccld2, ccld1))
            ccld (im,jm,km) = rhcldwgt*ccld1 + qcwgt*ccld2
!  borrow o3a to store psuedo-cloud
            psqc= 0.10*ccld1*qvsat
            IF (P_QC > 0) THEN
              o3a(i,j,k) = MAX(qscalar(i,j,k,P_QC),psqc)
            ELSE
              o3a(i,j,k) = MAX(0.0,psqc)
            END IF
          ELSE
            ccld (im,jm,km) = MIN(1.0, MAX(0.0, 0.25*ALOG10(tqe)+1.5))
!  borrow o3a to store psuedo-cloud
            IF (P_QC > 0) THEN
              o3a(i,j,k) = qscalar(i,j,k,P_QC)
            ELSE
              o3a(i,j,k) = 0.0
            END IF
          END IF

          psfc = 0.5* (pbar(i,j,1)+pprt(i,j,1) +                        &
                       pbar(i,j,2)+pprt(i,j,2) )
          aden=2.0*aeroden(pk,psfc)
          taual(im,jm,km)=aden*0.5*((pbar(i,j,k-1)+pprt(i,j,k-1))       &
                                   -(pbar(i,j,k+1)+pprt(i,j,k+1)))

        END DO
      END DO

      IF ( nyodd == 0 ) THEN
        DO im=2,m,2
          plinv(im,n,km) = plinv(im,n-1,km)
          tinv (im,n,km) = tinv (im,n-1,km)
          qvinv(im,n,km) = qvinv(im,n-1,km)
          ccld (im,n,km) = ccld (im,n-1,km)
          taual(im,n,km) = taual(im,n-1,km)
        END DO
      END IF

    END DO

!    WARNING: loops 180 & 190 are repeated above for radstgr=0, changes made
!    here must also be made there as well.
    DO j=1,ny-1
      jodd = istgr*MOD(j,2)
      DO i=jstgr-jodd, nx-1, jstgr
        im = i
        jm = j/jstgr + istgr*MOD(i,2)

        ijm = nx*(jm-1) + im
        radbuf(ijm) = tsfc(i,j)

        cosz (im,jm) = cosz(i,j)

        rsirbm(im,jm) = rsirbm(i,j)
        rsuvbm(im,jm) = rsuvbm(i,j)
        rsirdf(im,jm) = rsirdf(i,j)
        rsuvdf(im,jm) = rsuvdf(i,j)

        plinv(im,jm,1) = padd
        tinv (im,jm,1) = tadd
        qvinv(im,jm,1) = 1.0E-6
        ccld (im,jm,1) = 0.0
        taual(im,jm,1) = 0.0

        plinv(im,jm,2) = 0.5*(plinv(im,jm,1)+plinv(im,jm,3))
        tinv (im,jm,2) = MAX(0.5*(tinv(im,jm,1)+tinv(im,jm,3)),         &
                             190.0)
        qvinv(im,jm,2) = MAX(qv(i,j,nz-1), 1.0E-6)
        ccld (im,jm,2) = 0.0
        taual(im,jm,2) = 0.0

        plinv(im,jm,nz) = 0.005*((pbar(i,j,1) + pprt(i,j,1))            &
                                +(pbar(i,j,2) + pprt(i,j,2)))
      END DO
    END DO

    IF ( nyodd == 0 ) THEN
      DO im=2,m,2
        ijm = nx*(n-1) + im
        radbuf(ijm) = tsfc(im,n-1)

        cosz (im,n) = cosz(im,n-1)

        rsirbm(im,n) = rsirbm(im,n-1)
        rsuvbm(im,n) = rsuvbm(im,n-1)
        rsirdf(im,n) = rsirdf(im,n-1)
        rsuvdf(im,n) = rsuvdf(im,n-1)

        plinv(im,n,1) = padd
        tinv (im,n,1) = tadd
        qvinv(im,n,1) = 1.0E-6
        ccld (im,n,1) = 0.0
        taual(im,n,1) = 0.0

        plinv(im,n,2) = plinv(im,n-1,2)
        tinv (im,n,2) = tinv (im,n-1,2)
        qvinv(im,n,2) = qvinv(im,n-1,2)
        ccld (im,n,2) = 0.0
        taual(im,n,2) = 0.0

        plinv(im,n,nz) = plinv(im,n-1,nz)
      END DO
    END IF

  END IF

  DO km=1,nz-1
    DO j=1,n
      DO i=1,m
        dp = plinv(i,j,km+1) - plinv(i,j,km)
        IF ( dp <= 0.0 ) THEN
          WRITE (6,'(a,i3,a,i3,a,i3,a,i3,a/a,a/a)')                     &
              'ERROR: The pressure gradient between level k = ',nz+1-km, &
              ' and ',nz-km, ' at i = ',i,' and j = ',j,' was <=0.',    &
              'Please check the sounding file, ',                       &
              runname(1:lfnkey)//'.sound, or the data sets.',           &
              'Program stopped in RADTRNS.'
          CALL arpsstop('arpsstop stopped RADTRNS problem with sounding',1)
        END IF
      END DO
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Notes: The arguments in subroutine CLDOPTD have different vertical
!         coordinates orders:
!
!         plinv  -- 1 for top
!         tinv   -- 1 for top
!
!         qc     -- 1 for bottom
!         qr     -- 1 for bottom
!         qi     -- 1 for bottom
!         qs     -- 1 for bottom
!         qh     -- 1 for bottom
!
!         tauir  -- 1 for top
!         tauswi -- 1 for top
!         tauswl -- 1 for top
!         reffi  -- 1 for top
!         reffl  -- 1 for top
!
!-----------------------------------------------------------------------
!
!  note borrowed o3a to store psuedo-cloud in place of qc.
  CALL cldoptd(nx,ny,m,n,nz, radstgr,                                   &
               plinv,tinv,o3a,qscalar, rhostr,j3inv,                    &
               tauir,tauswi,tauswl,reffi,reffl)

  IF ( radstgr == 1 ) THEN
    IF ( nyodd == 0 ) THEN
      DO km=1,nz-1
        DO im=2,m,2
          tauir (im,n,km) = tauir (im,n-1,km)
          tauswi(im,n,km) = tauswi(im,n-1,km)
          tauswl(im,n,km) = tauswl(im,n-1,km)
          reffi (im,n,km) = reffi (im,n-1,km)
          reffl (im,n,km) = reffl (im,n-1,km)
        END DO
      END DO
    END IF
  END IF
!
!-----------------------------------------------------------------------
!
!  Fit the ozone concentration by interpolating a standard o3 profile
!  to ARPS pressure levels.
!
!-----------------------------------------------------------------------
!
  CALL fito3(nx,ny,m,n, nz-1, plinv,o3a)
!
!-----------------------------------------------------------------------
!
!  Calculate the downward longwave IR radiation.
!
!  Positions of 2-d arrays in the buffer:
!
!    fclr (m,n)     -- radbuf(1+ 1*nxy)
!    dbs  (m,n)     -- radbuf(1+ 2*nxy)
!    trant(m,n)     -- radbuf(1+ 3*nxy)
!
!    th2o (m,n,6)   -- radbuf(1+ 4*nxy)
!    tcon (m,n,3)   -- radbuf(1+10*nxy)
!    tco2 (m,n,6,2) -- radbuf(1+13*nxy)
!
!  Positions of 3-d arrays in the buffer:
!
!    pa    (m,n,np)     -- radbuf(1+25*nxy)
!    dt    (m,n,np)     -- radbuf(1+25*nxy+ 1*nxyz)
!    sh2o  (m,n,np+1)   -- radbuf(1+25*nxy+ 2*nxyz)
!    swpre (m,n,np+1)   -- radbuf(1+25*nxy+ 3*nxyz)
!    swtem (m,n,np+1)   -- radbuf(1+25*nxy+ 4*nxyz)
!    sco3  (m,n,np+1)   -- radbuf(1+25*nxy+ 5*nxyz)
!    scopre(m,n,np+1)   -- radbuf(1+25*nxy+ 6*nxyz)
!    scotem(m,n,np+1)   -- radbuf(1+25*nxy+ 7*nxyz)
!    dh2o  (m,n,np)     -- radbuf(1+25*nxy+ 8*nxyz)
!    dcont (m,n,np)     -- radbuf(1+25*nxy+ 9*nxyz)
!    dco2  (m,n,np)     -- radbuf(1+25*nxy+10*nxyz)
!    do3   (m,n,np)     -- radbuf(1+25*nxy+11*nxyz)
!    flxu  (m,n,np+1)   -- radbuf(1+25*nxy+12*nxyz)
!    flxd  (m,n,np+1)   -- radbuf(1+25*nxy+13*nxyz)
!    clr   (m,n,0:np+1) -- radbuf(1+25*nxy+14*nxyz) ! nz+1
!    blayer(m,n,0:np+1) -- radbuf(1+26*nxy+15*nxyz) ! nz+1
!
!    h2oexp(m,n,np,6)   -- radbuf(1+27*nxy+16*nxyz)
!    conexp(m,n,np,3)   -- radbuf(1+27*nxy+22*nxyz)
!
!    co2exp(m,n,np,6,2) -- radbuf(1+27*nxy+25*nxyz)
!
!-----------------------------------------------------------------------
!
  CALL irrad(nx,ny,m,n,nz-1,                                            &
             tauir,ccld, plinv,tinv,qvinv,o3a, co2,radbuf(1),           &
             high,flxd,flxir,flcir,dfdts,st4,                        &
             radbuf(ir2d1),radbuf(ir2d2),radbuf(ir2d3),                 &
             radbuf(ir2d4),radbuf(ir2d5),radbuf(ir2d6),                 &
             radbuf(ir3d1),radbuf(ir3d2),radbuf(ir3d3),                 &
             radbuf(ir3d4),radbuf(ir3d5),radbuf(ir3d6),                 &
             radbuf(ir3d7),radbuf(ir3d8),radbuf(ir3d9),                 &
             radbuf(ir3d10),radbuf(ir3d11),radbuf(ir3d12),              &
             radbuf(ir3d13),radbuf(ir3d14),radbuf(ir3d15),              &
             radbuf(ir3d16),radbuf(ir4d1),radbuf(ir4d2),                &
             radbuf(ir5d1))
!
!-----------------------------------------------------------------------
!
!  Calculate solar radiation fluxes.
!
!  Output flxuv and flcuv are the fractions to incoming solar flux
!  at the top of atmosphere
!
!  Positions of 2-d arrays used in subroutine SORAD
!
!    sdf  (m,n)     -- radbuf(1+ 1*nxy    )
!    sclr (m,n)     -- radbuf(1+ 2*nxy)
!    csm  (m,n)     -- radbuf(1+ 3*nxy)
!    cc   (m,n,3)   -- radbuf(1+ 4*nxy)
!
!  Positions of 3-d arrays used in subroutine SORAD
!
!    tauclb(m,n,np)   -- radbuf(1+ 7*nxy            )
!    tauclf(m,n,np)   -- radbuf(1+ 7*nxy+ 1*nxyz)
!    dp    (m,n,np)   -- radbuf(1+ 7*nxy+ 2*nxyz)
!    wh    (m,n,np)   -- radbuf(1+ 7*nxy+ 3*nxyz)
!    oh    (m,n,np)   -- radbuf(1+ 7*nxy+ 4*nxyz)
!    scal  (m,n,np)   -- radbuf(1+ 7*nxy+ 5*nxyz)
!    swh   (m,n,np+1) -- radbuf(1+ 7*nxy+ 6*nxyz)
!    so2   (m,n,np+1) -- radbuf(1+ 7*nxy+ 7*nxyz)
!    df    (m,n,np+1) -- radbuf(1+ 7*nxy+ 8*nxyz)
!
!  Positions of temporary arrays used in subroutine SOLIR, SOLUV, and
!  CLDFLX:
!
!    tem2d1 (m,n)         -- radbuf(1+ 7*nxy+ 9*nxyz)
!    tem2d2 (m,n)         -- radbuf(1+ 8*nxy+ 9*nxyz)
!    tem2d3 (m,n)         -- radbuf(1+ 9*nxy+ 9*nxyz)
!    tem2d4 (m,n)         -- radbuf(1+10*nxy+ 9*nxyz)
!    tem2d5 (m,n)         -- radbuf(1+11*nxy+ 9*nxyz)
!    tem2d6 (m,n)         -- radbuf(1+12*nxy+ 9*nxyz)
!    tem2d7 (m,n)         -- radbuf(1+13*nxy+ 9*nxyz)
!    tem2d8 (m,n)         -- radbuf(1+14*nxy+ 9*nxyz)
!    tem2d9 (m,n)         -- radbuf(1+15*nxy+ 9*nxyz)
!    tem2d10(m,n)         -- radbuf(1+16*nxy+ 9*nxyz)
!    tem2d11(m,n)         -- radbuf(1+17*nxy+ 9*nxyz)
!    tem2d12(m,n)         -- radbuf(1+18*nxy+ 9*nxyz)
!    tem2d13(m,n)         -- radbuf(1+19*nxy+ 9*nxyz)
!    tem2d14(m,n)         -- radbuf(1+20*nxy+ 9*nxyz)
!    tem2d15(m,n)         -- radbuf(1+21*nxy+ 9*nxyz)
!    tem2d16(m,n)         -- radbuf(1+22*nxy+ 9*nxyz)
!    tem2d17(m,n)         -- radbuf(1+23*nxy+ 9*nxyz)
!    tem2d18(m,n)         -- radbuf(1+24*nxy+ 9*nxyz)
!    tem2d19(m,n)         -- radbuf(1+25*nxy+ 9*nxyz)
!
!    tem3d1 (m,n,np+1)    -- radbuf(1+26*nxy+ 9*nxyz)
!    tem3d2 (m,n,np+1)    -- radbuf(1+26*nxy+10*nxyz)
!    tem3d3 (m,n,np+1)    -- radbuf(1+26*nxy+11*nxyz)
!    tem3d4 (m,n,np+1)    -- radbuf(1+26*nxy+12*nxyz)
!    tem3d5 (m,n,np+1)    -- radbuf(1+26*nxy+13*nxyz)
!
!    tem4d1 (m,n,np+1,2)  -- radbuf(1+26*nxy+14*nxyz)
!    tem4d2 (m,n,np+1,2)  -- radbuf(1+26*nxy+16*nxyz)
!    tem4d3 (m,n,np+1,2)  -- radbuf(1+26*nxy+18*nxyz)
!    tem4d4 (m,n,np+1,2)  -- radbuf(1+26*nxy+20*nxyz)
!    tem4d5 (m,n,np+1,2)  -- radbuf(1+26*nxy+22*nxyz)
!
!    tem5d1(m,n,np+1,2,2) -- radbuf(1+26*nxy+24*nxyz)
!    tem5d2(m,n,np+1,2,2) -- radbuf(1+26*nxy+28*nxyz)
!    tem5d3(m,n,np+1,2,2) -- radbuf(1+26*nxy+32*nxyz)
!    tem5d4(m,n,np+1,2,2) -- radbuf(1+26*nxy+36*nxyz)
!    tem5d5(m,n,np+1,2,2) -- radbuf(1+26*nxy+40*nxyz)
!
!-----------------------------------------------------------------------
!
  night = 1
  DO j=1,n
    DO i=1,m
      IF ( cosz(i,j) > 0.0 ) THEN
        night = 0
        GO TO 500
      END IF
    END DO
  END DO

  500   CONTINUE

!  aersum=0.
!  DO 505 k=1,nz-1
!    aersum=aersum+taual(2,2,k)
! 505 CONTINUE
!  write(6,'(a,f10.4))') ' Total aerosol optical depth: ',aersum

  IF ( night == 0 ) THEN

    CALL sorad(nx,ny,m,n,nz-1, plinv,tinv,qvinv,o3a,co2,                &
         tauswi,tauswl,reffi,reffl,ccld,ict,icb,                        &
         taual,rsirbm,rsirdf,rsuvbm,rsuvdf,cosz,                        &
         flxuv,flcuv,fdirir,fdifir,fdirpar,fdifpar,                     &
         radbuf(so2d1),radbuf(so2d2),radbuf(so2d3),radbuf(so2d4),       &
         radbuf(so3d1),radbuf(so3d2),radbuf(so3d3),radbuf(so3d4),       &
         radbuf(so3d5),radbuf(so3d6),radbuf(so3d7),radbuf(so3d8),       &
         radbuf(so3d9),radbuf(so2d5),radbuf(so2d6),radbuf(so2d7),       &
         radbuf(so2d8),radbuf(so2d9),radbuf(so2d10),radbuf(so2d11),     &
         radbuf(so2d12),radbuf(so2d13),radbuf(so2d14),                  &
         radbuf(so2d15),radbuf(so2d16),radbuf(so2d17),                  &
         radbuf(so2d18),radbuf(so2d19),radbuf(so2d20),                  &
         radbuf(so2d21),radbuf(so2d22),radbuf(so2d23),                  &
         radbuf(so3d10),radbuf(so3d11),radbuf(so3d12),                  &
         radbuf(so3d13),radbuf(so3d14),radbuf(so4d1),                   &
         radbuf(so4d2),radbuf(so4d3),radbuf(so4d4),radbuf(so4d5),       &
         radbuf(so5d1), radbuf(so5d2),                                  &
         radbuf(so5d3), radbuf(so5d4), radbuf(so5d5))

  ELSE

    DO k=1,nz-1
      DO j=1,n
        DO i=1,m
          flxuv(i,j,k) = 0.0
        END DO
      END DO
    END DO

    DO j=1,n
      DO i=1,m
        fdirir (i,j) = 0.0
        fdifir (i,j) = 0.0
        fdirpar(i,j) = 0.0
        fdifpar(i,j) = 0.0
      END DO
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Added the heating rate of solar radiation to total radiation
!  forcing (K/s)
!
!  Constant 9.770687e-05 is equal to g/cp in cgs unit, where g = 980
!  and cp = 1.003e7.
!
!  Outputs from SORAD such as flxuv, flcuv, etc., are fractions of
!  solar flux at the top of atmosphere. Therefore we need to multipy
!  solar constant and cosine of zenith angle to obtain the solar
!  radiation flux.
!
!-----------------------------------------------------------------------
!
  IF ( radstgr == 0 ) THEN
    DO k=2,nz-2
      km=nz+1-k          ! inverse vertical coordinates to ARPS grid
      DO j=1,ny-1
        DO i=1,nx-1
          coolrate = 9.770687E-05               & ! = g/cp in cgs unit
               * ( flxir(i,j,km+1) - flxir(i,j,km) )                    &
                   / ( plinv(i,j,km) - plinv(i,j,km+1) )

!        coolr (i,j,k) = -coolrate*86400.      ! cooling rate, C/day

          heatrate = solarc * cosz(i,j) * 9.770687E-05                  &
                   * ( flxuv(i,j,km+1) - flxuv(i,j,km) )                &
                   / ( plinv(i,j,km) - plinv(i,j,km+1) )

!        heatr (i,j,k) = heatrate*86400.0      ! heating rate, C/day

          radfrc(i,j,k) = (coolrate + heatrate)/ppi(i,j,k)
        END DO
      END DO
    END DO

    DO j=1,ny-1
      DO i=1,nx-1
! augustin add shade
! SHADE ONLY FOR DIRECT BEAM

        rnflx(i,j) = solarc * radsw(i,j) * cosss(i,j)       & ! radsw = a2dr2
                 * ( sh(i,j)*(1.0-rsirbm(i,j)) * fdirir(i,j)            &
                       + sh(i,j)*(1.0-rsuvbm(i,j)) * fdirpar(i,j)       &
                       + (1.0-rsirdf(i,j)) * fdifir(i,j)                &
                       + (1.0-rsuvdf(i,j)) * fdifpar(i,j) )             &
                   + flxir(i,j,nz)       ! net downward LW flux at sfc

        radswnet(i,j) = solarc * radsw(i,j) * cosss(i,j)       & ! radsw = a2dr2
                 * (sh(i,j)*(1.0-rsirbm(i,j)) * fdirir(i,j)             &
                 + sh(i,j)*(1.0-rsuvbm(i,j)) * fdirpar(i,j)             &
                          + (1.0-rsirdf(i,j)) * fdifir(i,j)             &
                          + (1.0-rsuvdf(i,j)) * fdifpar(i,j) )


! augustin
! BUG January 03
! The next line has been modified since the version of june 02
        radsw(i,j) = solarc * radsw(i,j) * cosss(i,j)                   &
                     * ( sh(i,j)*(fdirir(i,j) + fdirpar(i,j))           &
                       + fdifir(i,j) + fdifpar(i,j) )

        radlwin(i,j) = flxd(i,j)

      END DO
    END DO

  ELSE IF ( radstgr == 1 ) THEN

    DO k=2,nz-2
      km=nz+1-k          ! inverse vertical coordinates to ARPS grid
      DO j=1,ny-1
        jodd = istgr*MOD(j,2)
        DO i=jstgr-jodd, nx-1, jstgr
          im = i
          jm = j/jstgr + istgr*MOD(i,2)

          coolrate = 9.770687E-05               & ! = g/cp in cgs unit
               * ( flxir(im,jm,km+1) - flxir(im,jm,km) )                &
                   / ( plinv(im,jm,km) - plinv(im,jm,km+1) )

!        coolr(i,j,k) = -coolrate*86400.     ! cooling rate, C/day

          heatrate = solarc * cosz(im,jm) * 9.770687E-05                &
                   * ( flxuv(im,jm,km+1) - flxuv(im,jm,km) )            &
                   / ( plinv(im,jm,km) - plinv(im,jm,km+1) )

!        heatr(i,j,k) = heatrate*86400.0      ! heating rate, C/day

          radfrc(i,j,k) = (coolrate + heatrate)/ppi(i,j,k)
                          ! total radiation forcing (K/s)
        END DO
      END DO

      IF ( nx-1 > 3 ) THEN
        DO j=2,ny-2
          jeven = MOD(j+1,2)

          DO i=2+jeven,nx-2,2

!            coolr (i,j,k) = 0.25
!    :                        * (coolr (i+1,j,k)+coolr (i-1,j,k)
!    :                          +coolr (i,j+1,k)+coolr (i,j-1,k))

!            heatr (i,j,k) = 0.25
!    :                        * (heatr (i+1,j,k)+heatr (i-1,j,k)
!    :                          +heatr (i,j+1,k)+heatr (i,j-1,k))

            radfrc(i,j,k) = 0.25                                        &
                          * (radfrc(i+1,j,k)+radfrc(i-1,j,k)            &
                            +radfrc(i,j+1,k)+radfrc(i,j-1,k))
          END DO
        END DO
      END IF

      DO j=2,ny-2,2
!        coolr (1,j,k) = 0.25
!    :             *(coolr(2,j,k)+coolr(2,j,k)
!    :              +coolr(1,j+1,k)+coolr(1,j-1,k))

!        coolr (nx-1,j+nxodd,k) = 0.25
!    :             *(coolr(nx-2,j+nxodd,k)+coolr(nx-2,j+nxodd,k)
!    :              +coolr(nx-1,j+nxodd+1,k)+coolr(nx-1,j+nxodd-1,k))

!        heatr (1,j,k) = 0.25
!    :             *(heatr(2,j,k)+heatr(2,j,k)
!    :              +heatr(1,j+1,k)+heatr(1,j-1,k))

!        heatr (nx-1,j+nxodd,k) = 0.25
!    :             *(heatr(nx-2,j+nxodd,k)+heatr(nx-2,j+nxodd,k)
!    :              +heatr(nx-1,j+nxodd+1,k)+heatr(nx-1,j+nxodd-1,k))

        radfrc(1,j,k) = 0.25                                            &
               *(radfrc(2,j,k)+radfrc(2,j,k)                            &
                +radfrc(1,j+1,k)+radfrc(1,j-1,k))

        radfrc(nx-1,j+nxodd,k) = 0.25                                   &
               *(radfrc(nx-2,j+nxodd,k)+radfrc(nx-2,j+nxodd,k)          &
                +radfrc(nx-1,j+nxodd+1,k)+radfrc(nx-1,j+nxodd-1,k))
      END DO

      DO i=2,nx-2,2
!        coolr (i,1,k) = 0.25
!    :             *(coolr(i,2,k)+coolr(i,2,k)
!    :              +coolr(i+1,1,k)+coolr(i-1,1,k))

!        coolr (i+nyodd,ny-1,k) = 0.25
!    :             *(coolr(i+nyodd,ny-2,k)+coolr(i+nyodd,ny-2,k)
!    :              +coolr(i+nyodd+1,ny-1,k)+coolr(i+nyodd-1,ny-1,k))

!        heatr (i,1,k) = 0.25
!    :             *(heatr(i,2,k)+heatr(i,2,k)
!    :              +heatr(i+1,1,k)+heatr(i-1,1,k))

!        heatr (i+nyodd,ny-1,k) = 0.25
!    :             *(heatr(i+nyodd,ny-2,k)+heatr(i+nyodd,ny-2,k)
!    :              +heatr(i+nyodd+1,ny-1,k)+heatr(i+nyodd-1,ny-1,k))

        radfrc(i,1,k) = 0.25                                            &
               *(radfrc(i,2,k)+radfrc(i,2,k)                            &
                +radfrc(i+1,1,k)+radfrc(i-1,1,k))

        radfrc(i+nyodd,ny-1,k) = 0.25                                   &
               *(radfrc(i+nyodd,ny-2,k)+radfrc(i+nyodd,ny-2,k)          &
                +radfrc(i+nyodd+1,ny-1,k)+radfrc(i+nyodd-1,ny-1,k))
      END DO

      IF ( nxodd == 1 ) THEN
!        coolr (nx-1,1,k) = 0.5*(coolr (nx-2,1,k)+coolr (nx-1,2,k))

!        heatr (nx-1,1,k) = 0.5*(heatr (nx-2,1,k)+heatr (nx-1,2,k))

        radfrc(nx-1,1,k) = 0.5*(radfrc(nx-2,1,k)+radfrc(nx-1,2,k))
      END IF

      IF ( nyodd == 1 ) THEN
!        coolr (1,ny-1,k) = 0.5*(coolr (1,ny-2,k)+coolr (2,ny-1,k))

!        heatr (1,ny-1,k) = 0.5*(heatr (1,ny-2,k)+heatr (2,ny-1,k))

        radfrc(1,ny-1,k) = 0.5*(radfrc(1,ny-2,k)+radfrc(2,ny-1,k))
      END IF

      IF ( MOD(nyodd+nyodd,2) == 1 ) THEN
!        coolr (nx-1,ny-1,k) = 0.5*(coolr (nx-1,ny-2,k)
!    :                                +coolr (nx-2,ny-1,k))

!        heatr (nx-1,ny-1,k) = 0.5*(heatr (nx-1,ny-2,k)
!    :                                +heatr (nx-2,ny-1,k))

        radfrc(nx-1,ny-1,k) = 0.5*(radfrc(nx-1,ny-2,k)                  &
                                  +radfrc(nx-2,ny-1,k))
      END IF

    END DO

    DO j=1,ny-1
      jodd = istgr*MOD(j,2)
      DO i=jstgr-jodd, nx-1, jstgr
        im = i
        jm = j/jstgr + istgr*MOD(i,2)
!augustin
        rnflx(i,j) = solarc * radsw(i,j) * cosss(i,j)  & ! radsw = a2dr2
               * ( sh(i,j)*(1.0-rsirbm(im,jm)) * fdirir (im,jm)         &
                     + sh(i,j)*(1.0-rsuvbm(im,jm)) * fdirpar(im,jm)     &
                     + (1.0-rsirdf(im,jm)) * fdifir (im,jm)             &
                     + (1.0-rsuvdf(im,jm)) * fdifpar(im,jm) )           &
                   + flxir(im,jm,nz)   ! net downward LW flux at sfc



        radswnet(i,j) = solarc * radsw(i,j) * cosss(i,j)       & ! radsw = a2dr2
                 * (sh(i,j)*(1.0-rsirbm(im,jm)) * fdirir(im,jm)       &
                 + sh(i,j)*(1.0-rsuvbm(im,jm)) * fdirpar(im,jm)       &
                          + (1.0-rsirdf(im,jm)) * fdifir(im,jm)       &
                          + (1.0-rsuvdf(im,jm)) * fdifpar(im,jm) )

        radsw(i,j) = solarc * radsw(i,j) * cosss(i,j)                   &
                   * ( sh(i,j)*(fdirir(im,jm) + fdirpar(im,jm))         &
                     + fdifir(im,jm) + fdifpar(im,jm) )

        radlwin(i,j) = flxd(im,jm)

      END DO
    END DO

    IF ( nx-1 > 3 ) THEN
      DO j=2,ny-2
        jeven = MOD(j+1,2)
        DO i=2+jeven,nx-2,2
          radsw(i,j) = 0.25*(radsw(i+1,j)+radsw(i-1,j)                  &
                            +radsw(i,j+1)+radsw(i,j-1))

          radswnet(i,j) = 0.25*(radswnet(i+1,j)+radswnet(i-1,j)         &
                               +radswnet(i,j+1)+radswnet(i,j-1))

          rnflx(i,j) = 0.25*(rnflx(i+1,j)+rnflx(i-1,j)                  &
                            +rnflx(i,j+1)+rnflx(i,j-1))

          radlwin(i,j) = 0.25*(radlwin(i+1,j)+radlwin(i-1,j)            &
                            +radlwin(i,j+1)+radlwin(i,j-1))

        END DO
      END DO
    END IF

    DO j=2,ny-2,2
      radsw(1,j) = 0.25                                                 &
          * (radsw(2,j)+radsw(2,j)                                      &
           +radsw(1,j+1)+radsw(1,j-1))
      radsw(nx-1,j+nxodd) = 0.25                                        &
          * (radsw(nx-2,j+nxodd)+radsw(nx-2,j+nxodd)                    &
           +radsw(nx-1,j+nxodd+1)+radsw(nx-1,j+nxodd-1))

      radswnet(1,j) = 0.25                                              &
          * (radswnet(2,j)+radswnet(2,j)                                &
           +radswnet(1,j+1)+radswnet(1,j-1))
      radswnet(nx-1,j+nxodd) = 0.25                                     &
          * (radswnet(nx-2,j+nxodd)+radswnet(nx-2,j+nxodd)              &
           +radswnet(nx-1,j+nxodd+1)+radswnet(nx-1,j+nxodd-1))

      rnflx(1,j) = 0.25                                                 &
          * (rnflx(2,j)+rnflx(2,j)                                      &
           +rnflx(1,j+1)+rnflx(1,j-1))
      rnflx(nx-1,j+nxodd) = 0.25                                        &
          * (rnflx(nx-2,j+nxodd)+rnflx(nx-2,j+nxodd)                    &
           +rnflx(nx-1,j+nxodd+1)+rnflx(nx-1,j+nxodd-1))

      radlwin(1,j) = 0.25                                               &
          * (radlwin(2,j)+radlwin(2,j)                                  &
           +radlwin(1,j+1)+radlwin(1,j-1))
      radlwin(nx-1,j+nxodd) = 0.25                                      &
          * (radlwin(nx-2,j+nxodd)+radlwin(nx-2,j+nxodd)                &
           +radlwin(nx-1,j+nxodd+1)+radlwin(nx-1,j+nxodd-1))

    END DO

    DO i=2,nx-2,2
      radsw(i,1) = 0.25                                                 &
          * (radsw(i,2)  + radsw(i,2)                                      &
           + radsw(i+1,1)+ radsw(i-1,1))
      radsw(i+nyodd,ny-1) = 0.25                                        &
          * (radsw(i+nyodd,ny-2)+radsw(i+nyodd,ny-2)                    &
           +radsw(i+nyodd+1,ny-1)+radsw(i+nyodd-1,ny-1))

      radswnet(i,1) = 0.25                                              &
          * (radswnet(i,2)  + radswnet(i,2)                                &
           + radswnet(i+1,1)+ radswnet(i-1,1))
      radswnet(i+nyodd,ny-1) = 0.25                                     &
          * (radswnet(i+nyodd,ny-2)+radswnet(i+nyodd,ny-2)              &
           +radswnet(i+nyodd+1,ny-1)+radswnet(i+nyodd-1,ny-1))

      rnflx(i,1) = 0.25                                                 &
          * (rnflx(i,2)+rnflx(i,2)                                      &
           +rnflx(i+1,1)+rnflx(i-1,1))
      rnflx(i+nyodd,ny-1) = 0.25                                        &
          * (rnflx(i+nyodd,ny-2)+rnflx(i+nyodd,ny-2)                    &
           +rnflx(i+nyodd+1,ny-1)+rnflx(i+nyodd-1,ny-1))

      radlwin(i,1) = 0.25                                               &
          * (radlwin(i,2)+radlwin(i,2)                                  &
           +radlwin(i+1,1)+radlwin(i-1,1))
      radlwin(i+nyodd,ny-1) = 0.25                                      &
          * (radlwin(i+nyodd,ny-2)+radlwin(i+nyodd,ny-2)                &
           +radlwin(i+nyodd+1,ny-1)+radlwin(i+nyodd-1,ny-1))

    END DO

    IF ( nxodd == 1 ) THEN
      radsw(nx-1,1) = 0.5*(radsw(nx-2,1)+radsw(nx-1,2))
      radswnet(nx-1,1) = 0.5*(radswnet(nx-2,1)+radswnet(nx-1,2))
      rnflx(nx-1,1) = 0.5*(rnflx(nx-2,1)+rnflx(nx-1,2))
      radlwin(nx-1,1) = 0.5*(radlwin(nx-2,1)+radlwin(nx-1,2))
    END IF

    IF ( nyodd == 1 ) THEN
      radsw(1,ny-1) = 0.5*(radsw(1,ny-2)+radsw(2,ny-1))
      radswnet(1,ny-1) = 0.5*(radswnet(1,ny-2)+radswnet(2,ny-1))
      rnflx(1,ny-1) = 0.5*(rnflx(1,ny-2)+rnflx(2,ny-1))
      radlwin(1,ny-1) = 0.5*(radlwin(1,ny-2)+radlwin(2,ny-1))
    END IF

    IF ( MOD(nyodd+nyodd,2) == 1 ) THEN
      radsw(nx-1,ny-1) = 0.5*(radsw(nx-1,ny-2)+radsw(nx-2,ny-1))
      radswnet(nx-1,ny-1) = 0.5*(radswnet(nx-1,ny-2)+radswnet(nx-2,ny-1))
      rnflx(nx-1,ny-1) = 0.5*(rnflx(nx-1,ny-2)+rnflx(nx-2,ny-1))
      radlwin(nx-1,ny-1) = 0.5*(radlwin(nx-1,ny-2)+radlwin(nx-2,ny-1))
    END IF
  END IF

!---------------------------------------------------------------------
!
! Added by Yunheng to update the fake zone for radfrc, radsw, rnflx.
! two more variables (radswnet and radlwin) since IHOP_3
!
!--------------------------------------------------------------------

  IF (mp_opt > 0) THEN
    CALL acct_interrupt(mp_acct)
    CALL mpsendrecv2dew(radfrc, nx, ny, nz, ebc, wbc, 0, tem1)
    CALL mpsendrecv2dns(radfrc, nx, ny, nz, nbc, sbc, 0, tem1)

    CALL mpsendrecv1dew(radsw,  nx, ny, ebc, wbc, 0, tem1)
    CALL mpsendrecv1dns(radsw,  nx, ny, nbc, sbc, 0, tem1)

    CALL mpsendrecv1dew(rnflx,  nx, ny, ebc, wbc, 0, tem1)
    CALL mpsendrecv1dns(rnflx,  nx, ny, nbc, sbc, 0, tem1)

    CALL mpsendrecv1dew(radswnet, nx, ny, ebc, wbc, 0, tem1)
    CALL mpsendrecv1dns(radswnet, nx, ny, nbc, sbc, 0, tem1)

    CALL mpsendrecv1dew(radlwin,  nx, ny, ebc, wbc, 0, tem1)
    CALL mpsendrecv1dns(radlwin,  nx, ny, nbc, sbc, 0, tem1)

    CALL acct_stop_inter
  END IF

  IF ( raddiag == 1 ) THEN
    WRITE(6,'(a,i8,a,f10.2,a)')                                         &
        ' Dump radiation variables at time step,', nstep,               &
        ', model time=',curtim,' (s)'
!
!-----------------------------------------------------------------------
!
!  Write out results to GrADS file for display
!
!-----------------------------------------------------------------------
!
    CALL wrtrad(nx,ny,nz,m,n,x,y,z,                                     &
                plinv,tinv,qvinv,qscalar, o3a,radbuf(1),                &
                ccld, tauir,taual,tauswi,tauswl,reffi,reffl,            &
                rsirbm,rsirdf,rsuvbm,rsuvdf,                            &
                fdirir,fdifir,fdirpar,fdifpar,                          &
                dfdts, radsw,rnflx, cosz,                               &
                flxir,flcir, flxuv,flcuv,                               &
                radfrc)
!    :              radfrc, coolr,heatr)

  END IF

  RETURN

END SUBROUTINE radtrns
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SETRADWRK                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE setradwrk( nx,ny,nz )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Set the indeces for radiation working arrays
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  04/09/1997
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the z-direction (vertical)
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

  INTEGER :: nx,ny,nz       ! The number grid points in 3 directions
  INTEGER :: nxy, nxyz
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'radcst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  nxy  = nx*ny
  nxyz = nxy*nz
!
!-----------------------------------------------------------------------
!
!  Define indices which determine the positions of temporary arrays
!  used in subroutine IRRAD and the subroutine IRRAD calls.
!
!-----------------------------------------------------------------------
!
  ir2d1  = 1+ 1*nxy
  ir2d2  = 1+ 2*nxy
  ir2d3  = 1+ 3*nxy
  ir2d4  = 1+ 4*nxy
  ir2d5  = 1+10*nxy
  ir2d6  = 1+13*nxy

  ir3d1  = 1+25*nxy+ 0*nxyz
  ir3d2  = 1+25*nxy+ 1*nxyz
  ir3d3  = 1+25*nxy+ 2*nxyz
  ir3d4  = 1+25*nxy+ 3*nxyz
  ir3d5  = 1+25*nxy+ 4*nxyz
  ir3d6  = 1+25*nxy+ 5*nxyz
  ir3d7  = 1+25*nxy+ 6*nxyz
  ir3d8  = 1+25*nxy+ 7*nxyz
  ir3d9  = 1+25*nxy+ 8*nxyz
  ir3d10 = 1+25*nxy+ 9*nxyz
  ir3d11 = 1+25*nxy+10*nxyz
  ir3d12 = 1+25*nxy+11*nxyz
  ir3d13 = 1+25*nxy+12*nxyz
  ir3d14 = 1+25*nxy+13*nxyz
  ir3d15 = 1+25*nxy+14*nxyz
  ir3d16 = 1+26*nxy+15*nxyz

  ir4d1  = 1+27*nxy+16*nxyz
  ir4d2  = 1+27*nxy+22*nxyz

  ir5d1  = 1+27*nxy+25*nxyz
!
!-----------------------------------------------------------------------
!
!  Define indices which determine the positions of temporary arrays
!  used in subroutine SOLIR, SOLUV, and CLDFLX.
!
!-----------------------------------------------------------------------
!
  so2d1  = 1+ 1*nxy
  so2d2  = 1+ 2*nxy
  so2d3  = 1+ 3*nxy
  so2d4  = 1+ 4*nxy
  so2d5  = 1+ 7*nxy
  so2d6  = 1+ 8*nxy
  so2d7  = 1+ 9*nxy
  so2d8  = 1+10*nxy
  so2d9  = 1+11*nxy
  so2d10 = 1+12*nxy
  so2d11 = 1+13*nxy
  so2d12 = 1+14*nxy
  so2d13 = 1+15*nxy
  so2d14 = 1+16*nxy
  so2d15 = 1+17*nxy
  so2d16 = 1+18*nxy
  so2d17 = 1+19*nxy
  so2d18 = 1+20*nxy
  so2d19 = 1+21*nxy
  so2d20 = 1+22*nxy
  so2d21 = 1+23*nxy
  so2d22 = 1+24*nxy
  so2d23 = 1+25*nxy

  so3d1  = 1+26*nxy+ 0*nxyz
  so3d2  = 1+26*nxy+ 1*nxyz
  so3d3  = 1+26*nxy+ 2*nxyz
  so3d4  = 1+26*nxy+ 3*nxyz
  so3d5  = 1+26*nxy+ 4*nxyz
  so3d6  = 1+26*nxy+ 5*nxyz
  so3d7  = 1+26*nxy+ 6*nxyz
  so3d8  = 1+26*nxy+ 7*nxyz
  so3d9  = 1+26*nxy+ 8*nxyz
  so3d10 = 1+26*nxy+ 9*nxyz
  so3d11 = 1+26*nxy+10*nxyz
  so3d12 = 1+26*nxy+11*nxyz
  so3d13 = 1+26*nxy+12*nxyz
  so3d14 = 1+26*nxy+13*nxyz

  so4d1  = 1+26*nxy+14*nxyz
  so4d2  = 1+26*nxy+16*nxyz
  so4d3  = 1+26*nxy+18*nxyz
  so4d4  = 1+26*nxy+20*nxyz
  so4d5  = 1+26*nxy+22*nxyz

  so5d1  = 1+26*nxy+24*nxyz
  so5d2  = 1+26*nxy+28*nxyz
  so5d3  = 1+26*nxy+32*nxyz
  so5d4  = 1+26*nxy+36*nxyz
  so5d5  = 1+26*nxy+40*nxyz

  RETURN
END SUBROUTINE setradwrk
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRTRAD                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE wrtrad(nx,ny,nz,m,n,x,y,z,                                   &
           plinv,tinv,qvinv,qscalar,o3a,tsfc,                           &
           ccld, tauir,taual,tauswi,tauswl,reffi,reffl,                 &
           rsirbm,rsirdf,rsuvbm,rsuvdf,                                 &
           fdirir,fdifir,fdirpar,fdifpar,                               &
           dfdts, radsw,rnflx, cosz,                                    &
           flxir,flcir,flxuv,flcuv,                                     &
           radfrc)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Write surface fields in GrADS format for diagnostic purpose.
!
!  In the staggering case, only radfrc, coolr, and heatr have values
!  at entire domain. Other variables have their values on every
!  (i,j) = (even,even) and (odd,odd) point and they are defragmented
!  with index i in the range between 1 and nx/2.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  03/13/1996
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the z-direction (vertical)
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz       ! The number grid points in 3 directions
  INTEGER :: m, n

  REAL :: x(nx)
  REAL :: y(ny)
  REAL :: z(nz)

  REAL :: plinv (nx,ny,nz)
  REAL :: tinv  (nx,ny,nz)
  REAL :: qvinv (nx,ny,nz)

  REAL :: qscalar(nx,ny,nz,nscalar)

  REAL :: o3a   (nx,ny,nz)
  REAL :: tsfc  (nx,ny)

  REAL :: ccld  (nx,ny,nz)

  REAL :: tauir (nx,ny,nz)
  REAL :: taual (nx,ny,nz)
  REAL :: tauswi(nx,ny,nz)
  REAL :: tauswl(nx,ny,nz)
  REAL :: reffi (nx,ny,nz)
  REAL :: reffl (nx,ny,nz)

  REAL :: rsirbm(nx,ny)
  REAL :: rsirdf(nx,ny)
  REAL :: rsuvbm(nx,ny)
  REAL :: rsuvdf(nx,ny)

  REAL :: fdirir(nx,ny)
  REAL :: fdifir(nx,ny)
  REAL :: fdirpar(nx,ny)
  REAL :: fdifpar(nx,ny)

  REAL :: dfdts(nx,ny,nz)
  REAL :: cosz  (nx,ny)

  REAL :: flxir(nx,ny,nz)
  REAL :: flcir(nx,ny,nz)
  REAL :: flxuv(nx,ny,nz)
  REAL :: flcuv(nx,ny,nz)

  REAL :: radfrc(nx,ny,nz)

  REAL :: radsw(nx,ny)
  REAL :: rnflx(nx,ny)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k
  INTEGER :: nq

  CHARACTER (LEN=5)   :: qcloud(5) = (/'Cloud','Rain ','Ice  ','Snow ','Hail '/)

  INTEGER :: tnum,tint
  INTEGER :: year1,month1,day1, hour1,minute1,second1
  INTEGER :: jday1, loopdy

  CHARACTER (LEN=2) :: dtunit

  INTEGER :: mndys(12)                 ! days for each months

  CHARACTER (LEN=3)   :: monnam(12)
  CHARACTER (LEN=256) :: flnctl, flnrad
  INTEGER :: radunit, flnctlen, radlen
  INTEGER :: ierr
  LOGICAL :: firstcall
!
!-----------------------------------------------------------------------
!
!  Save and initialize variables.
!
!-----------------------------------------------------------------------
!
  SAVE firstcall, radunit
  DATA firstcall/.true./
  DATA mndys/0,31,59,90,120,151,181,212,243,273,304,334/
  DATA monnam/'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',                 &
              'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'/
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF ( firstcall ) THEN

    IF ( dtrad <= 0.0 ) THEN
      WRITE (6, '(/a,a)')                                               &
          'Since dtrad <= 0, only data at the first time step ',        &
          'will be dumped.'
      tnum = 1
      tint = 1
      dtunit = 'MN'
    ELSE IF ( dtrad < 60.0 ) THEN
      WRITE (6, '(/a/a)')                                               &
          'GrADS reqiures the smallest uint minute for time interval.', &
          'Here we use uint MN to represent the second.'
      tnum = nint(tstop/dtrad) + 1
      tint = nint(dtrad)
      dtunit = 'MN'
    ELSE IF ( dtrad < 3600.0 ) THEN
      tnum = nint(tstop/dtrad) + 1
      tint = nint(dtrad/60.)
      dtunit = 'MN'
    ELSE IF ( dtrad < 86400.0 ) THEN
      tnum = nint(tstop/dtrad) + 1
      tint = nint(dtrad/3600.)
      dtunit = 'HR'
    ELSE
      tnum = nint(tstop/dtrad) + 1
      tint = nint(dtrad/86400.)
      dtunit = 'DY'
    END IF

    IF ( initopt /= 2 ) THEN
      second1 = second
      minute1 = minute
      hour1   = hour
      day1    = day
      month1  = month
      year1   = year
    ELSE
      second1 = MOD( second + nint(tstart), 60 )
      minute1 = ( second + nint(tstart) ) / 60
      minute1 = MOD( minute + minute1, 60 )
      hour1   = ( minute + ( second + nint(tstart) ) / 60 ) /60
      hour1   = MOD( hour + hour1, 24 )
      day1    = ( hour + ( minute                                       &
            + ( second + nint(tstart) ) / 60 ) /60 ) / 24
      jday1   = jday + day1

      loopdy  = 0
      IF ( MOD( year, 4 ) == 0 ) loopdy = 1
      year1 = year + jday1 / ( 365 + loopdy )
      jday1 = MOD( jday1, 365 + loopdy )

      month1 = 1

      DO m = 2, 11
        IF ( jday1 > mndys(m) .AND. jday1 <= mndys(m+1) + loopdy ) month1 = m
      END DO
      day1 = jday1 - mndys(month1)

    END IF

    flnctlen = lfnkey + 7
    flnctl(1:flnctlen) = runname(1:lfnkey)//'.radctl'
    CALL fnversn( flnctl, flnctlen )

    flnrad(1:ldirnam) = dirname(1:ldirnam)

    radlen = ldirnam + lfnkey + 8
    flnrad(1:radlen) = flnrad(1:ldirnam)//'/'//runname(1:lfnkey)        &
                     //'.radout'
    CALL fnversn( flnrad, radlen )
!
!-----------------------------------------------------------------------
!
!  Open GrADS data control file for surface variables.
!
!-----------------------------------------------------------------------
!
    CALL getunit (radunit)
    OPEN (UNIT = radunit, FILE = flnctl(1:flnctlen),                    &
          FORM = 'formatted', STATUS = 'new')

    WRITE (6,'(a,a)')                                                   &
        'Radiation output control file ', flnctl(1:flnctlen)

    WRITE (radunit,'(a)')                                               &
        'TITLE   Radiation output for '//runname(1:lfnkey)
    WRITE (radunit,'(a)')                                               &
        '*'
    WRITE (radunit,'(a,a)')                                             &
        'DSET    ', flnrad(1:radlen)
    WRITE (radunit,'(a)')                                               &
        'OPTIONS sequential cray_32bit_ieee'
    WRITE (radunit,'(a)')                                               &
        'UNDEF   -9.e+33'
    WRITE (radunit,'(a,i8,a,2f10.4)')                                   &
        'XDEF    ', nx, '  LINEAR   ',(x(1)+x(2))/2000.,dx/1000.
    WRITE (radunit,'(a,i8,a,2f10.4)')                                   &
        'YDEF    ', ny, '  LINEAR   ',(y(1)+y(2))/2000.,dy/1000.
    WRITE (radunit,'(a,1x,i8,a)')                                       &
        'ZDEF',nz,'  LEVELS  '
    WRITE (radunit,'(8f10.4)')                                          &
        ((z(k)+z(k+1))/2000.,k=1,nz-1),z(nz)/1000.
    WRITE (radunit,'(a,i8,a,i2.2,a,i2.2,a,i2.2,a3,i4.4,3X,i2.2,a)')     &
        'TDEF    ', tnum, '  LINEAR   ',                                &
        hour1,':',minute1,'Z',day1,monnam(month1),year1,                &
        tint,dtunit

    WRITE (radunit,'(a)')                                               &
        '*'
    WRITE (radunit,'(a,I3)')                                            &
        'VARS  ',29+nscalar

    WRITE (radunit,'(a,2x,i4,2x,a)')                                    &
        'p    ', nz,'99  Pressure (mb)'
    WRITE (radunit,'(a,2x,i4,2x,a)')                                    &
        't    ', nz,'99  Temperature (K)'
    WRITE (radunit,'(a,2x,i4,2x,a)')                                    &
        'qv   ', nz,'99  Specific humidity (g/g)'

    DO nq = 1,nscalar
      WRITE (radunit,'(a,2x,i4,2x,a,a)')                                  &
          qnames(nq), nz,'99  ',TRIM(qdescp(nq))
    END DO

    WRITE (radunit,'(a,2x,i4,2x,a)')                                    &
        'o3   ', nz,'99  Ozone concentration (g/g)'
    WRITE (radunit,'(a,2x,i4,2x,a)')                                    &
        'cld  ', nz,'99  Cloud coverage (%)'
    WRITE (radunit,'(a,2x,i4,2x,a)')                                    &
        'taual', nz,'99  Cloud optical depth for aeresol'
    WRITE (radunit,'(a,2x,i4,2x,a)')                                    &
        'tauir', nz,'99  Cloud optical depth for IR'
    WRITE (radunit,'(a,2x,i4,2x,a)')                                    &
        'tausi', nz,'99  Ice cloud optical depth for SW'
    WRITE (radunit,'(a,2x,i4,2x,a)')                                    &
        'tausl', nz,'99  Liquid cloud optical depth for SW'
    WRITE (radunit,'(a,2x,i4,2x,a)')                                    &
        'reffi', nz,'99  Eff. ice-cloud-particle size (mm)'
    WRITE (radunit,'(a,2x,i4,2x,a)')                                    &
        'reffl', nz,'99  Eff. liquid-cloud-particle size (mm)'

    WRITE (radunit,'(a,2x,i4,2x,a)')                                    &
        'flxir', nz,'99  all-sky IR flux'
    WRITE (radunit,'(a,2x,i4,2x,a)')                                    &
        'flcir', nz,'99  Cloud IR flux'
    WRITE (radunit,'(a,2x,i4,2x,a)')                                    &
        'flxuv', nz,'99  all-sky UV flux'
    WRITE (radunit,'(a,2x,i4,2x,a)')                                    &
        'flcuv', nz,'99  Clear-sky UV flux'
    WRITE (radunit,'(a,2x,i4,2x,a)')                                    &
        'dfdts', nz,'99  sensitivity to surface temperature'

    WRITE (radunit,'(a,2x,i4,2x,a)')                                    &
        'radfc', nz,'99  radiation forcing (K/s)'

!      write (radunit,'(a,2x,i4,2x,a)')
!    :      'coolr', nz,'99  radiation cooling rate (K/day)'
!      write (radunit,'(a,2x,i4,2x,a)')
!    :      'heatr', nz,'99  radiation heating rate (K/day)'

    WRITE (radunit,'(a)')                                               &
        'ts        0  99  Soil temperature (K)'

    WRITE (radunit,'(a)')                                               &
        'albd1     0  99  solar ir surface albedo for beam'
    WRITE (radunit,'(a)')                                               &
        'albd2     0  99  solar ir surface albedo for diffuse'
    WRITE (radunit,'(a)')                                               &
        'albd3     0  99  solar uv surface albedo for beam'
    WRITE (radunit,'(a)')                                               &
        'albd4     0  99  solar uv surface albedo for diffuse'

    WRITE (radunit,'(a)')                                               &
        'cosz      0  99  Cosine of zenith'
    WRITE (radunit,'(a)')                                               &
        'radsw     0  99  solar flux reaching the surface'
    WRITE (radunit,'(a)')                                               &
        'rnflx     0  99  net radiation flux absorbed by surface'

    WRITE (radunit,'(a)')                                               &
        'dirir     0  99  all-sky direct downward ir at sfc'
    WRITE (radunit,'(a)')                                               &
        'difir     0  99  all-sky diffuse downward ir at sfc'
    WRITE (radunit,'(a)')                                               &
        'dirpa     0  99  all-sky direct downward par at sfc'
    WRITE (radunit,'(a)')                                               &
        'difpa     0  99  all-sky diffuse downward par at sfc'

    WRITE (radunit,'(a)')                                               &
        'ENDVARS'

    CLOSE (radunit)
    CALL retunit (radunit)
!
!-----------------------------------------------------------------------
!
!  Open GrADS data file for surface variables.
!
!-----------------------------------------------------------------------
!
    CALL getunit (radunit)

    CALL asnctl ('NEWLOCAL', 1, ierr)
    CALL asnfile(flnrad(1:radlen), '-F f77 -N ieee', ierr)

    OPEN (UNIT = radunit, FILE = flnrad(1:radlen),                      &
          FORM = 'unformatted', STATUS = 'new',                         &
        ACCESS = 'sequential')

    firstcall = .false.

  END IF

  CALL edgfill(plinv, 1,nx,1,m, 1,ny,1,n, 1,nz,1,nz)
  DO k=nz,1,-1
    WRITE(radunit) ((plinv(i,j,k),i=1,nx),j=1,ny)
  END DO
  CALL edgfill(tinv, 1,nx,1,m, 1,ny,1,n, 1,nz,1,nz-1)
  DO k=nz,1,-1
    WRITE(radunit) ((tinv(i,j,k),i=1,nx),j=1,ny)
  END DO

  CALL edgfill(qvinv, 1,nx,1,m, 1,ny,1,n, 1,nz,1,nz-1)
  DO k=nz,1,-1
    WRITE(radunit) ((qvinv(i,j,k),i=1,nx),j=1,ny)
  END DO

  DO nq = 1,nscalar
    CALL edgfill(qscalar(:,:,:,nq),1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    DO k=1,nz
      WRITE(radunit) ((qscalar(i,j,k,nq),i=1,nx),j=1,ny)
    END DO
  END DO

  CALL edgfill(o3a, 1,nx,1,m, 1,ny,1,n, 1,nz,1,nz-1)
  DO k=nz,1,-1
    WRITE(radunit) ((o3a(i,j,k),i=1,nx),j=1,ny)
  END DO
  CALL edgfill(ccld, 1,nx,1,m, 1,ny,1,n, 1,nz,1,nz-1)
  DO k=nz,1,-1
    WRITE(radunit) ((ccld(i,j,k),i=1,nx),j=1,ny)
  END DO
  CALL edgfill(taual, 1,nx,1,m, 1,ny,1,n, 1,nz,1,nz-1)
  DO k=nz,1,-1
    WRITE(radunit) ((taual(i,j,k),i=1,nx),j=1,ny)
  END DO
  CALL edgfill(tauir, 1,nx,1,m, 1,ny,1,n, 1,nz,1,nz-1)
  DO k=nz,1,-1
    WRITE(radunit) ((tauir(i,j,k),i=1,nx),j=1,ny)
  END DO
  CALL edgfill(tauswi, 1,nx,1,m, 1,ny,1,n, 1,nz,1,nz-1)
  DO k=nz,1,-1
    WRITE(radunit) ((tauswi(i,j,k),i=1,nx),j=1,ny)
  END DO
  CALL edgfill(tauswl, 1,nx,1,m, 1,ny,1,n, 1,nz,1,nz-1)
  DO k=nz,1,-1
    WRITE(radunit) ((tauswl(i,j,k),i=1,nx),j=1,ny)
  END DO
  CALL edgfill(reffi, 1,nx,1,m, 1,ny,1,n, 1,nz,1,nz-1)
  DO k=nz,1,-1
    WRITE(radunit) ((reffi(i,j,k),i=1,nx),j=1,ny)
  END DO
  CALL edgfill(reffl, 1,nx,1,m, 1,ny,1,n, 1,nz,1,nz-1)
  DO k=nz,1,-1
    WRITE(radunit) ((reffl(i,j,k),i=1,nx),j=1,ny)
  END DO

  CALL edgfill(flxir, 1,nx,1,m, 1,ny,1,n, 1,nz,1,nz)
  DO k=nz,1,-1
    WRITE(radunit) ((flxir(i,j,k),i=1,nx),j=1,ny)
  END DO
  CALL edgfill(flcir, 1,nx,1,m, 1,ny,1,n, 1,nz,1,nz)
  DO k=nz,1,-1
    WRITE(radunit) ((flcir(i,j,k),i=1,nx),j=1,ny)
  END DO
  CALL edgfill(flxuv, 1,nx,1,m, 1,ny,1,n, 1,nz,1,nz)
  DO k=nz,1,-1
    WRITE(radunit) ((flxuv(i,j,k),i=1,nx),j=1,ny)
  END DO
  CALL edgfill(flcuv, 1,nx,1,m, 1,ny,1,n, 1,nz,1,nz)
  DO k=nz,1,-1
    WRITE(radunit) ((flcuv(i,j,k),i=1,nx),j=1,ny)
  END DO
  CALL edgfill(dfdts, 1,nx,1,m, 1,ny,1,n, 1,nz,1,nz)
  DO k=nz,1,-1
    WRITE(radunit) ((dfdts(i,j,k),i=1,nx),j=1,ny)
  END DO

  CALL edgfill(radfrc, 1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,2,nz-2)
  DO k=1,nz
    WRITE(radunit) ((radfrc(i,j,k),i=1,nx),j=1,ny)
  END DO

!     CALL edgfill(coolr,  1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,2,nz-2)
!     DO 300 k=1,nz
!       write(radunit) ((coolr(i,j,k),i=1,nx),j=1,ny)
!300     CONTINUE
!     CALL edgfill(heatr,  1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,2,nz-2)
!     DO 310 k=1,nz
!       write(radunit) ((heatr(i,j,k),i=1,nx),j=1,ny)
!310     CONTINUE

  CALL edgfill(tsfc,   1,nx,1,m, 1,ny,1,n, 1,1,1,1)
  WRITE(radunit) tsfc

  CALL edgfill(rsirbm, 1,nx,1,m, 1,ny,1,n, 1,1,1,1)
  WRITE(radunit) rsirbm
  CALL edgfill(rsirdf, 1,nx,1,m, 1,ny,1,n, 1,1,1,1)
  WRITE(radunit) rsirdf
  CALL edgfill(rsuvbm, 1,nx,1,m, 1,ny,1,n, 1,1,1,1)
  WRITE(radunit) rsuvbm
  CALL edgfill(rsuvdf, 1,nx,1,m, 1,ny,1,n, 1,1,1,1)
  WRITE(radunit) rsuvdf

  CALL edgfill(cosz, 1,nx,1,m, 1,ny,1,n, 1,1,1,1)
  WRITE(radunit) cosz
  CALL edgfill(radsw, 1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
  WRITE(radunit) radsw
  CALL edgfill(rnflx, 1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
  WRITE(radunit) rnflx

  CALL edgfill(fdirir, 1,nx,1,m, 1,ny,1,n, 1,1,1,1)
  WRITE(radunit) fdirir
  CALL edgfill(fdifir, 1,nx,1,m, 1,ny,1,n, 1,1,1,1)
  WRITE(radunit) fdifir
  CALL edgfill(fdirpar, 1,nx,1,m, 1,ny,1,n, 1,1,1,1)
  WRITE(radunit) fdirpar
  CALL edgfill(fdifpar, 1,nx,1,m, 1,ny,1,n, 1,1,1,1)
  WRITE(radunit) fdifpar

  RETURN
END SUBROUTINE wrtrad
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CLDOPTD                    ######
!######                                                      ######
!######                     Developed by                     ######
!######                                                      ######
!######    Goddard Cumulus Ensemble Modeling Group, NASA     ######
!######                                                      ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE cldoptd(nx,ny,m,n,nz, radstgrl,                              &
           pres,temp,qc,qscalar, rhostr, j3inv,                         &
           tauir,tauswi,tauswl,reffi,reffl)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the optical depth
!
!-----------------------------------------------------------------------
!
!  AUTHOR: NASA Goddard Center
!
!  MODIFICATION:
!
!  03/11/1996 (Yuhe Liu)
!  Modified the original code from 1-D to 3-D
!
!  07/25/2011 (Yunheng Wang)
!  Upgraded to the latest multimoment module.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the z-direction (vertical)
!
!                                               Vertical index order
!    pres     Pressure (mb)                       1 for top
!    temp     Temperature (K),                    1 for top
!    qc       Cloud water mixing ratio (g/g),     1 for bottom
!    qr       Rain water mixing ratio (g/g),      1 for bottom
!    qi       cloud ice mixing ratio (g/g),       1 for bottom
!    qs       Snow mixing ratio (g/g),            1 for bottom
!    qh       Hail mixing ratio (g/g),            1 for bottom
!    rhostr   Density multiply by j3,             1 for bottom
!    j3inv    1/j3,                               1 for bottom
!
!  OUTPUT:
!
!    tauir    Cloud optical depth for longwave,   1 for bottom
!    tauswi   Cloud optical depth for ice cloud
!             for shortwave,                      1 for bottom
!    tauswl   Cloud optical depth for liquid cloud
!             for shortwave,                      1 for bottom
!    reffi    Effective cloud-particle size (mm)
!             for ice cloud for shortwave,        1 for bottom
!    reffl    Effective cloud-particle size (mm)
!             for liquid cloud for shortwave,     1 for bottom
!
!  WORK ARRAY:
!
!-----------------------------------------------------------------------
!

! DTD: Added to make available functions related to multi-moment code

  !USE MM_FUNCTIONS
  USE my3mom_fncs_mod, ONLY: gammaDP, diagAlpha_v33, solveAlpha_v33

  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
  INCLUDE 'radcst.inc'
  INCLUDE 'globcst.inc'

!-----------------------------------------------------------------------
!
! Variable declaration
!
!-----------------------------------------------------------------------

  INTEGER :: nx,ny,nz
  INTEGER :: m,n
  INTEGER :: radstgrl           ! local value to avoid conflict with global
                                ! value declared in globcst.inc
  REAL :: pres  (nx,ny,nz)      ! in mb
  REAL :: temp  (nx,ny,nz)      ! in K
  REAL :: qc    (nx,ny,nz)      ! in g/g  Note that qc is passed in separately

  REAL :: qscalar (nx,ny,nz,nscalar)      ! in g/g

  REAL :: rhostr(nx,ny,nz)      ! in kg/m**3
  REAL :: j3inv (nx,ny,nz)

  REAL :: tauir (nx,ny,nz)
  REAL :: tauswi(nx,ny,nz)
  REAL :: tauswl(nx,ny,nz)
  REAL :: reffi (nx,ny,nz)
  REAL :: reffl (nx,ny,nz)

  REAL :: tauqc
  REAL :: tauqr
  REAL :: tauqi
  REAL :: tauqs
  REAL :: tauqg
  REAL :: tauqh
  REAL :: reff1
  REAL :: reff2

  REAL :: w1, effrad, w2

  REAL :: cpi,twco, dpg, rho

  INTEGER :: i,j,k, im,jm,km
  INTEGER :: istgr, jstgr, jodd

  REAL :: cx,rhoMKS,Dm,alpha
  REAL*8, PARAMETER :: muc = 3.0, alphac = 1.0  ! Assumed fixed shape parameters in
                                              ! MY cloud distribution

!  REAL*8 :: gamma
!  REAL*8 :: diagAlpha
!  REAL*8 :: solveAlpha

  ! DTD: Added below parameter to test the treatment of ice and snow
  ! in the radiation calculations.  Originally (option 1), ice and snow are combined
  ! and the optical depths for each are derived according to a temperature-based
  ! formula for effective radius.
  !
  ! The new option (option 2) is to treat the ice and snow separately, and similarly to how
  ! hail and rain are treated (i.e. with effrad calculated from the ratio of the 3rd
  ! and 2nd moments), but only for the MY MM schemes, which have gamma
  ! distributions for ice and snow.  The MY1 scheme, and the other single-moment
  ! schemes will still use the original treatment.  It is unclear how accurate this
  ! modification will be, so it is added with an in-code testing switch for now.
  ! EDIT: 08/08/08  It appears the new option has problems so it is recommended to use the old
  ! formulation (option 1) for now.

  INTEGER, PARAMETER :: radiceopt = 1
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  cpi  = 4.*ATAN(1.)
  twco = 1.e-6

  IF ( radstgrl == 0 ) THEN
    istgr = 0
    jstgr = 1
  ELSE
    istgr = 1
    jstgr = 2
  END IF

  DO km=3,nz-1
    k = nz+1-km

    DO j=1,ny-1
      jodd = istgr*MOD(j,2)
      DO i=jstgr-jodd, nx-1, jstgr

        im = i
        jm = j/jstgr + istgr*MOD(im,2)

        rho = rhostr(i,j,k)*j3inv(i,j,k)*1.0E-3        ! g/cm**3
        dpg = 10.0*(pres(im,jm,km+1)-pres(im,jm,km))/g ! g/cm**2

        IF ( qc(i,j,k) >= twco ) THEN
          w1     = dpg*qc(i,j,k)
          ! DTD: Calculate effrad of cloud for MY scheme to be consistent
          ! with the assumed generalized gammaDP distribution.  For other schemes
          ! leave it the way it was before (i.e. constant effrad = 0.0015)
          ! Note that the computation of effrad is slightly more complicated in
          ! the case of cloud than of rain, due to the extra parameter in the
          ! gammaDP distribution.
          IF (mphyopt == 8) THEN ! MY1 scheme
            ! The following calculation is in MKS units
            effrad = sngl(gammaDP(1.d0+alphac+3.d0/muc)/gammaDP(1.d0+alphac+2.d0/muc)) * &
                     ((sngl(gammaDP(1.d0+alphac))*rho*1.0e3*qscalar(i,j,k,P_QC))/ &
                     (sngl(gammaDP(1.d0+alphac+3.d0/muc))*(cpi/6.)*roqr*1000.*ntcloud))** &
                     (1./3.)
            ! Convert to cm
            effrad = effrad*100.
          ELSE IF (mphyopt >=9 .and. mphyopt <= 11) THEN ! MY2, 2.5, or 3 scheme
            ! The following calculation is in MKS units
            IF(qscalar(i,j,k,P_QC) > 0.0 .and. qscalar(i,j,k,P_NC) > 0.0) THEN
              effrad = sngl(gammaDP(1.d0+alphac+3.d0/muc)/gammaDP(1.d0+alphac+2.d0/muc)) * &
                       ((sngl(gammaDP(1.d0+alphac))*rho*1.0e3*qscalar(i,j,k,P_QC))/ &
                       (sngl(gammaDP(1.d0+alphac+3.d0/muc))*(cpi/6.)*roqr*1000.*qscalar(i,j,k,P_NC)))** &
                       (1./3.)
            ELSE
              !print*,'QC,NC',qscalar(i,j,k,P_QC),qscalar(i,j,k,P_NC)
              effrad = 0.0
            END IF
            ! Convert to cm
            effrad = effrad*100.
          ELSE ! Other schemes (assumed constant effrad as originally for now)
            effrad = 0.0015
          END IF
          !effrad = 0.0015
          IF(effrad == 0) THEN
            !print*,'effradc = 0!'
          END IF
          IF(effrad /= 0.0) THEN
            tauqc  = w1/effrad
          ELSE
            tauqc = 0.0
          END IF
          reff2  = effrad*1.0E4
        ELSE
          tauqc = 0.0
          reff2 = 0.0
        END IF

        tauqr=0.0
        IF (P_QR > 0) THEN
          IF ( qscalar(i,j,k,P_QR) >= twco ) THEN
            w1     = dpg*qscalar(i,j,k,P_Qr)
            !DTD: Changed the following calculation for effrad to be consistent
            !with the appropriate microphysics package
            IF (mphyopt == 8) THEN ! MY1 scheme
              ! The following calculation is in MKS units
              effrad = sngl(gammaDP(4.d0+dble(alpharain))/gammaDP(3.d0+dble(alpharain))) * &
                       ((rho*1.0e3*qscalar(i,j,k,P_QR))/(sngl(gammaDP(4.d0+alpharain)) * &
                       (cpi/6.)*roqr*1000.*n0rain))**(1./(4.+alpharain))
              ! Convert to cm
              effrad = effrad*100.
            ELSE IF (mphyopt == 9) THEN ! MY2 scheme
              ! The following calculation is in MKS units
              IF(qscalar(i,j,k,P_QR) > 0.0 .and. qscalar(i,j,k,P_NR) > 0.0) THEN
                effrad = sngl(gammaDP(4.d0+dble(alpharain))/gammaDP(3.d0+dble(alpharain))) * &
                         ((sngl(gammaDP(1.d0+dble(alpharain)))*rho*1.0e3*qscalar(i,j,k,P_QR))/ &
                         (sngl(gammaDP(4.d0+alpharain))*(cpi/6.)*roqr*1000.*qscalar(i,j,k,P_NR)))** &
                         (1./3.)
              ELSE
                effrad = 0.0
              END IF
              ! Convert to cm
              effrad = effrad*100.
            ELSE IF (mphyopt == 10) THEN ! MY2.5 scheme
              ! First calculate mean-mass diameter
              Dm = (rho*1.0e3*qscalar(i,j,k,P_QR)/((cpi/6.)*roqr*1000.*qscalar(i,j,k,P_NR)))**(1./3.)
              ! Diagnose alpha
              alpha = sngl(diagAlpha_v33(dble(Dm),1))
              ! The following calculation is in MKS units
              IF(qscalar(i,j,k,P_QR) > 0.0 .and. qscalar(i,j,k,P_NR) > 0.0) THEN
                effrad = sngl(gammaDP(4.d0+dble(alpha))/gammaDP(3.d0+dble(alpha))) * &
                         ((sngl(gammaDP(1.d0+dble(alpha)))*rho*1.0e3*qscalar(i,j,k,P_QR))/ &
                         (sngl(gammaDP(4.d0+alpha))*(cpi/6.)*roqr*1000.*qscalar(i,j,k,P_NR)))** &
                         (1./3.)
              ELSE
                effrad = 0.0
              END IF
              ! Convert to cm
              effrad = effrad*100.
            ELSE IF (mphyopt == 11) THEN ! MY3 scheme
              ! First calculate alpha
              cx = (cpi/6.)*roqr*1000.
              rhoMKS = rho*1.0e3
              IF(qscalar(i,j,k,P_QR) > 0.0 .and. qscalar(i,j,k,P_NR) > 0.0 .and. qscalar(i,j,k,P_ZR) > 0.0 ) THEN
                alpha = max(0., sngl(solveAlpha_v33(qscalar(i,j,k,P_QR),qscalar(i,j,k,P_NR), &
                        qscalar(i,j,k,P_ZR),cx,rhoMKS)))
              ELSE
                alpha = 0.0
              END IF
              IF(alpha /= 0.0) THEN
                !print*,'alphar = ',alpha
              END IF
              ! The following calculation is in MKS units
              IF(qscalar(i,j,k,P_QR) > 0.0 .and. qscalar(i,j,k,P_NR) > 0.0 .and. qscalar(i,j,k,P_ZR) > 0.0 ) THEN
                effrad = sngl(gammaDP(4.d0+dble(alpha))/gammaDP(3.d0+dble(alpha))) * &
                         ((sngl(gammaDP(1.d0+dble(alpha)))*rho*1.0e3*qscalar(i,j,k,P_QR))/ &
                         (sngl(gammaDP(4.d0+alpha))*(cpi/6.)*roqr*1000.*qscalar(i,j,k,P_NR)))** &
                         (1./3.)
              ELSE
                effrad = 0.0
              END IF
              ! Convert to cm
              effrad = effrad*100.
            ELSE ! Other schemes (assumed exponential distribution for now)
              effrad = 3./((cpi*tnw*roqr/(rho*qscalar(i,j,k,P_Qr)))**.25)  ! this appears to be CGS units
            END IF                                                             ! it would be nice if it was
                                                                         ! explicitly stated so.
            !effrad = 3./((cpi*tnw*roqr/(rho*qscalar(i,j,k,P_Qr)))**.25)
            IF(effrad /= 0.0) THEN
              tauqr  = w1/effrad
            ELSE
              tauqr = 0.0
            END IF
          END IF
        END IF

        tauqi = 0.0
        tauqs = 0.0
        reff1 = 0.0
        IF (((P_QI > 0 .OR. P_QS > 0) .and. radiceopt == 1) .or. mphyopt <= 8) THEN  ! DTD: Original formulation
          w2 = 0.0               ! newly added for sum of QI+QS
          IF (P_QI > 0) w2 = w2+qscalar(i,j,k,P_QI)
          IF (P_QS > 0) w2 = w2+qscalar(i,j,k,P_QS)

          IF ( w2 >= twco ) THEN
            w1 = 1.e4*dpg*w2

            IF ( temp(im,jm,km) > 243.16 ) THEN
              effrad = 0.0125
            ELSE IF ( temp(im,jm,km) < 223.16 ) THEN
              effrad = 0.0025
            ELSE
              effrad = 0.0125+(temp(im,jm,km)-243.16)*0.00050
            END IF

            tauqs = w1*(-0.006656 +  3.686E-4/effrad)
            tauqi = w1*(-0.011500 +  4.110E-4/effrad                      &
                                  + 17.300E-8/(effrad*effrad))
            reff1 = effrad*1.0E4
          END IF
        END IF

        IF (P_QI > 0 .and. radiceopt == 2 .and. mphyopt >= 9) THEN ! New option for ice (treated same as cloud,rain,graupel,hail)
          IF ( qscalar(i,j,k,P_QI) >= twco ) THEN
            w1     = dpg*qscalar(i,j,k,P_QI)
            !DTD: Changed the following calculation for effrad to be consistent
            !with the appropriate microphysics package
            IF (mphyopt == 9) THEN ! MY2 scheme
              ! The following calculation is in MKS units
              effrad = sngl(gammaDP(4.d0+dble(alphaice))/gammaDP(3.d0+dble(alphaice))) * &
                       ((sngl(gammaDP(1.d0+dble(alphaice)))*rho*1.0e3*qscalar(i,j,k,P_QI))/ &
                       (sngl(gammaDP(4.d0+alphaice))*440.0*qscalar(i,j,k,P_NI)))** &
                       (1./3.)
              ! Convert to cm
              effrad = effrad*100.
            ELSE IF (mphyopt == 10) THEN ! MY2.5 scheme
              ! First calculate mean-mass diameter
              Dm = (rho*1.0e3*qscalar(i,j,k,P_QI)/(440.0*qscalar(i,j,k,P_NI)))**(1./3.)
              ! Diagnose alpha
              alpha = sngl(diagAlpha_v33(dble(Dm),2))
              ! The following calculation is in MKS units
              effrad = sngl(gammaDP(4.d0+dble(alpha))/gammaDP(3.d0+dble(alpha))) * &
                       ((sngl(gammaDP(1.d0+dble(alpha)))*rho*1.0e3*qscalar(i,j,k,P_QI))/ &
                       (sngl(gammaDP(4.d0+alpha))*440.0*qscalar(i,j,k,P_NI)))** &
                       (1./3.)
              ! Convert to cm
              effrad = effrad*100.
            ELSE IF (mphyopt == 11) THEN ! MY3 scheme
              ! First calculate alpha
              cx = 440.0
              rhoMKS = rho*1.0e3
              IF(qscalar(i,j,k,P_NI) > 0.0 .and. qscalar(i,j,k,P_ZI) > 0.0 ) THEN
                alpha = max(0., sngl(solveAlpha_v33(qscalar(i,j,k,P_QI),qscalar(i,j,k,P_NI), &
                        qscalar(i,j,k,P_ZI),cx,rhoMKS)))
              ELSE
                alpha = 0.0
              END IF
              ! The following calculation is in MKS units
              effrad = sngl(gammaDP(4.d0+dble(alpha))/gammaDP(3.d0+dble(alpha))) * &
                       ((sngl(gammaDP(1.d0+dble(alpha)))*rho*1.0e3*qscalar(i,j,k,P_QI))/ &
                       (sngl(gammaDP(4.d0+alpha))*cx*qscalar(i,j,k,P_NI)))** &
                       (1./3.)
              ! Convert to cm
              effrad = effrad*100.
!            ELSE ! Other schemes (assumed exponential distribution for now)
!              effrad = 3./((cpi*tng*roqg/(rho*qscalar(i,j,k,P_QH)))**.25)
            END IF
            tauqi  = w1/effrad
            IF(tauqi /= 0.0) THEN
              !print*,'tauqi = ',tauqi
            END IF
            IF(effrad /= 0.0) THEN
              !print*,'effrad = ',effrad
            END IF
            reff1 = effrad*1.0e4
          END IF
        END IF

        IF (P_QS > 0 .and. radiceopt == 2 .and. mphyopt >= 9) THEN ! New option for snow (treated same as cloud,rain,graupel,hail)
          IF ( qscalar(i,j,k,P_QS) >= twco ) THEN
            w1     = dpg*qscalar(i,j,k,P_QS)
            !DTD: Changed the following calculation for effrad to be consistent
            !with the appropriate microphysics package
            IF (mphyopt == 9) THEN ! MY2 scheme
              ! The following calculation is in MKS units
              effrad = sngl(gammaDP(4.d0+dble(alphasnow))/gammaDP(3.d0+dble(alphasnow))) * &
                       ((sngl(gammaDP(1.d0+dble(alphasnow)))*rho*1.0e3*qscalar(i,j,k,P_QS))/ &
                       (sngl(gammaDP(4.d0+alphasnow))*(cpi/6.)*rhosnow*qscalar(i,j,k,P_NS)))** &
                       (1./3.)
              ! Convert to cm
              effrad = effrad*100.
            ELSE IF (mphyopt == 10) THEN ! MY2.5 scheme
              ! First calculate mean-mass diameter
              Dm = (rho*1.0e3*qscalar(i,j,k,P_QS)/((cpi/6.)*rhosnow*qscalar(i,j,k,P_NS)))**(1./3.)
              ! Diagnose alpha
              alpha = sngl(diagAlpha_v33(dble(Dm),3))
              ! The following calculation is in MKS units
              effrad = sngl(gammaDP(4.d0+dble(alpha))/gammaDP(3.d0+dble(alpha))) * &
                       ((sngl(gammaDP(1.d0+dble(alpha)))*rho*1.0e3*qscalar(i,j,k,P_QS))/ &
                       (sngl(gammaDP(4.d0+alpha))*(cpi/6.)*rhosnow*qscalar(i,j,k,P_NS)))** &
                       (1./3.)
              ! Convert to cm
              effrad = effrad*100.
            ELSE IF (mphyopt == 11) THEN ! MY3 scheme
              ! First calculate alpha
              cx = (cpi/6.)*rhosnow
              rhoMKS = rho*1.0e3
              IF(qscalar(i,j,k,P_NS) > 0.0 .and. qscalar(i,j,k,P_ZS) > 0.0 ) THEN
                alpha = max(0., sngl(solveAlpha_v33(qscalar(i,j,k,P_QS),qscalar(i,j,k,P_NS), &
                        qscalar(i,j,k,P_ZS),cx,rhoMKS)))
              ELSE
                alpha = 0.0
              END IF
              ! The following calculation is in MKS units
              IF(qscalar(i,j,k,P_NS) /= 0.0) THEN
                effrad = sngl(gammaDP(4.d0+dble(alpha))/gammaDP(3.d0+dble(alpha))) * &
                         ((sngl(gammaDP(1.d0+dble(alpha)))*rho*1.0e3*qscalar(i,j,k,P_QS))/ &
                         (sngl(gammaDP(4.d0+alpha))*(cpi/6.)*rhosnow*qscalar(i,j,k,P_NS)))** &
                         (1./3.)
              ELSE
                effrad = 0.0
              END IF
              ! Convert to cm
              effrad = effrad*100.
!            ELSE ! Other schemes (assumed exponential distribution for now)
!              effrad = 3./((cpi*tng*roqg/(rho*qscalar(i,j,k,P_QH)))**.25)
            END IF
            tauqs  = w1/effrad
            IF(tauqs /= 0.0) THEN
              !print*,'tauqs = ',tauqs
            END IF
            reff1 = (reff1 + effrad*1.0e4)/2.  ! DTD: For now use simple average of effective radius for ice and snow
          END IF
        END IF
        ! DTD: added block for graupel
        tauqg  = 0.0
        IF (P_QG > 0) THEN
          IF ( qscalar(i,j,k,P_QG) >= twco ) THEN
            w1     = dpg*qscalar(i,j,k,P_QG)
            !DTD: Changed the following calculation for effrad to be consistent
            !with the appropriate microphysics package
            IF (mphyopt == 8) THEN ! MY1 scheme
              ! The following calculation is in MKS units
              effrad = sngl(gammaDP(4.d0+dble(alphagrpl))/gammaDP(3.d0+dble(alphagrpl))) * &
                       ((rho*1.0e3*qscalar(i,j,k,P_QG))/(sngl(gammaDP(4.d0+alphagrpl)) * &
                       (cpi/6.)*rhogrpl*n0grpl))**(1./(4.+alphagrpl))
              ! Convert to cm
              effrad = effrad*100.
            ELSE IF (mphyopt == 9) THEN ! MY2 scheme
              IF(qscalar(i,j,k,P_QG) > 0.0 .and. qscalar(i,j,k,P_NG) > 0.0) THEN
                ! The following calculation is in MKS units
                effrad = sngl(gammaDP(4.d0+dble(alphagrpl))/gammaDP(3.d0+dble(alphagrpl))) * &
                       ((sngl(gammaDP(1.d0+dble(alphagrpl)))*rho*1.0e3*qscalar(i,j,k,P_QG))/ &
                       (sngl(gammaDP(4.d0+alphagrpl))*(cpi/6.)*rhogrpl*qscalar(i,j,k,P_NG)))** &
                       (1./3.)
                ! Convert to cm
                effrad = effrad*100.
              ELSE
                effrad = 0.0
              END IF
            ELSE IF (mphyopt == 10) THEN ! MY2.5 scheme
              IF(qscalar(i,j,k,P_QG) > 0.0 .and. qscalar(i,j,k,P_NG) > 0.0) THEN
                ! First calculate mean-mass diameter
                Dm = (rho*1.0e3*qscalar(i,j,k,P_QG)/((cpi/6.)*rhogrpl*qscalar(i,j,k,P_NG)))**(1./3.)
                ! Diagnose alpha
                alpha = sngl(diagAlpha_v33(dble(Dm),4))
                ! The following calculation is in MKS units
                effrad = sngl(gammaDP(4.d0+dble(alpha))/gammaDP(3.d0+dble(alpha))) * &
                       ((sngl(gammaDP(1.d0+dble(alpha)))*rho*1.0e3*qscalar(i,j,k,P_QG))/ &
                       (sngl(gammaDP(4.d0+alpha))*(cpi/6.)*rhogrpl*qscalar(i,j,k,P_NG)))** &
                       (1./3.)
                ! Convert to cm
                effrad = effrad*100.
              ELSE
                effrad = 0.0
              END IF
            ELSE IF (mphyopt == 11) THEN ! MY3 scheme
              ! First calculate alpha
              cx = (cpi/6.)*rhogrpl
              rhoMKS = rho*1.0e3
              IF(qscalar(i,j,k,P_QG) > 0.0 .and. qscalar(i,j,k,P_NG) > 0.0 .and. qscalar(i,j,k,P_ZG) > 0.0 ) THEN
                alpha = max(0., sngl(solveAlpha_v33(qscalar(i,j,k,P_QG),qscalar(i,j,k,P_NG), &
                        qscalar(i,j,k,P_ZG),cx,rhoMKS)))
              ELSE
                alpha = 0.0
              END IF
              ! The following calculation is in MKS units
              IF(qscalar(i,j,k,P_QG) > 0.0 .and. qscalar(i,j,k,P_NG) > 0.0 .and. qscalar(i,j,k,P_ZG) > 0.0) THEN
                effrad = sngl(gammaDP(4.d0+dble(alpha))/gammaDP(3.d0+dble(alpha))) * &
                         ((sngl(gammaDP(1.d0+dble(alpha)))*rho*1.0e3*qscalar(i,j,k,P_QG))/ &
                         (sngl(gammaDP(4.d0+alpha))*(cpi/6.)*rhogrpl*qscalar(i,j,k,P_NG)))** &
                         (1./3.)
              ELSE
                effrad = 0.0
              END IF
              ! Convert to cm
              effrad = effrad*100.
!            ELSE ! Other schemes (assumed exponential distribution for now)
!              effrad = 3./((cpi*tng*roqg/(rho*qscalar(i,j,k,P_QH)))**.25)
            END IF
            IF(effrad /= 0.0) THEN
              tauqg  = w1/effrad
            ELSE
              tauqg = 0.0
            END IF
          END IF
        END IF

        tauqh  = 0.0
        IF (P_QH > 0) THEN
          IF ( qscalar(i,j,k,P_QH) >= twco ) THEN
            w1     = dpg*qscalar(i,j,k,P_QH)
            !DTD: Changed the following calculation for effrad to be consistent
            !with the appropriate microphysics package
            IF (mphyopt == 8) THEN ! MY1 scheme
              ! The following calculation is in MKS units
              effrad = sngl(gammaDP(4.d0+dble(alphahail))/gammaDP(3.d0+dble(alphahail))) * &
                       ((rho*1.0e3*qscalar(i,j,k,P_QH))/(sngl(gammaDP(4.d0+alphahail)) * &
                       (cpi/6.)*rhohail*n0hail))**(1./(4.+alphahail))
              ! Convert to cm
              effrad = effrad*100.
            ELSE IF (mphyopt == 9) THEN ! MY2 scheme
              IF(qscalar(i,j,k,P_QH) > 0.0 .and. qscalar(i,j,k,P_NH) > 0.0) THEN
                ! The following calculation is in MKS units
                effrad = sngl(gammaDP(4.d0+dble(alphahail))/gammaDP(3.d0+dble(alphahail))) * &
                       ((sngl(gammaDP(1.d0+dble(alphahail)))*rho*1.0e3*qscalar(i,j,k,P_QH))/ &
                       (sngl(gammaDP(4.d0+alphahail))*(cpi/6.)*rhohail*qscalar(i,j,k,P_NH)))** &
                       (1./3.)
                ! Convert to cm
                effrad = effrad*100.
              ELSE
                effrad = 0.0
              END IF
            ELSE IF (mphyopt == 10) THEN ! MY2.5 scheme
              IF(qscalar(i,j,k,P_QH) > 0.0 .and. qscalar(i,j,k,P_NH) > 0.0) THEN
                ! First calculate mean-mass diameter
                Dm = (rho*1.0e3*qscalar(i,j,k,P_QH)/((cpi/6.)*rhohail*qscalar(i,j,k,P_NH)))**(1./3.)
                ! Diagnose alpha
                alpha = sngl(diagAlpha_v33(dble(Dm),5))
                ! The following calculation is in MKS units
                effrad = sngl(gammaDP(4.d0+dble(alpha))/gammaDP(3.d0+dble(alpha))) * &
                       ((sngl(gammaDP(1.d0+dble(alpha)))*rho*1.0e3*qscalar(i,j,k,P_QH))/ &
                       (sngl(gammaDP(4.d0+alpha))*(cpi/6.)*rhohail*qscalar(i,j,k,P_NH)))** &
                       (1./3.)
                ! Convert to cm
                effrad = effrad*100.
              ELSE
                effrad = 0.0
              END IF
            ELSE IF (mphyopt == 11) THEN ! MY3 scheme
              ! First calculate alpha
              cx = (cpi/6.)*rhohail
              rhoMKS = rho*1.0e3
              IF(qscalar(i,j,k,P_QH) > 0.0 .and. qscalar(i,j,k,P_NH) > 0.0 .and. qscalar(i,j,k,P_ZH) > 0.0 ) THEN
                alpha = max(0., sngl(solveAlpha_v33(qscalar(i,j,k,P_QH),qscalar(i,j,k,P_NH), &
                        qscalar(i,j,k,P_ZH),cx,rhoMKS)))
              ELSE
                alpha = 0.0
              END IF
              ! The following calculation is in MKS units
              IF(qscalar(i,j,k,P_QH) > 0.0 .and. qscalar(i,j,k,P_NH) > 0.0 .and. qscalar(i,j,k,P_ZH) > 0.0) THEN
                effrad = sngl(gammaDP(4.d0+dble(alpha))/gammaDP(3.d0+dble(alpha))) * &
                         ((sngl(gammaDP(1.d0+dble(alpha)))*rho*1.0e3*qscalar(i,j,k,P_QH))/ &
                         (sngl(gammaDP(4.d0+alpha))*(cpi/6.)*rhohail*qscalar(i,j,k,P_NH)))** &
                         (1./3.)
              ELSE
                effrad = 0.0
              END IF
              ! Convert to cm
              effrad = effrad*100.
            ELSE ! Other schemes (assumed exponential distribution for now)
              effrad = 3./((cpi*tng*roqg/(rho*qscalar(i,j,k,P_QH)))**.25)
            END IF
!            effrad = 3./((cpi*tng*roqg/(rho*qscalar(i,j,k,P_QH)))**.25)
            IF(effrad /= 0.0) THEN
              tauqh  = w1/effrad
            ELSE
              tauqh = 0.0
            END IF
          END IF
        END IF

        tauswi(im,jm,km) = tauqs + tauqh + tauqg
        tauswl(im,jm,km) = 1.5 * ( tauqc + tauqr )
        reffi (im,jm,km) = reff1
        reffl (im,jm,km) = reff2
        tauir (im,jm,km)  = 0.5 * tauswl(im,jm,km) + tauqi + tauqg + tauqh
      END DO
    END DO
  END DO

  IF ( radstgrl /= 0 .AND. MOD(ny,2) == 0 ) THEN
    DO km=3,nz-1
      DO im=2,m,2
        tauswi(im,n,km) = tauswi(im,n-1,km)
        tauswl(im,n,km) = tauswl(im,n-1,km)
        reffi (im,n,km) = reffi (im,n-1,km)
        reffl (im,n,km) = reffl (im,n-1,km)
        tauir (im,n,km) = tauir (im,n-1,km)
      END DO
    END DO
  END IF

  DO jm=1,n
    DO im=1,m
      tauswi(im,jm,1) = 0.0
      tauswl(im,jm,1) = 0.0
      reffi (im,jm,1) = 0.0
      reffl (im,jm,1) = 0.0
      tauir (im,jm,1) = 0.0

      tauswi(im,jm,2) = 0.0
      tauswl(im,jm,2) = 0.0
      reffi (im,jm,2) = 0.0
      reffl (im,jm,2) = 0.0
      tauir (im,jm,2) = 0.0
    END DO
  END DO

  RETURN
END SUBROUTINE cldoptd
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE FITO3                      ######
!######                                                      ######
!######                     Developed by                     ######
!######                                                      ######
!######    Goddard Cumulus Ensemble Modeling Group, NASA     ######
!######                                                      ######
!######     Center for Analysis and Prediction of Storms     ######
!######               University of Oklahoma                 ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE fito3(nx,ny,m,n,np,pl,ao)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  This subroutine is to fit o3 to the model grid
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (a) Radiative Transfer Model: M.-D. Chou and M. Suarez
!          (b) Cloud Optics:Tao, Lang, Simpson, Sui, Ferrier and
!              Chou (1996)
!
!  MODIFICATION HISTORY:
!
!  03/15/1996 (Yuhe Liu)
!  Modified the subroutine from 1-D to 3-D
!
!-----------------------------------------------------------------------
!
!fpp$ expand (terp1)
!!dir$ inline always terp1
!*$*  inline routine (terp1)
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,np
  INTEGER :: m,n

  REAL :: pl(nx,ny,np+1)      ! Model pressure (mb)
  REAL :: ao(nx,ny,np)        ! Model o3 mixing ratio (g/g)
!
!-----------------------------------------------------------------------
!
!  Local definitions
!
!-----------------------------------------------------------------------
!
  INTEGER :: lay
  PARAMETER (lay=75)

!  integer iop(2),itab(3),iflag(5)
  INTEGER :: iop(2),itab(3)

  REAL :: pa(lay)       ! Local layer pressure (mb)
  REAL :: ta(lay)       ! Local layer temperature (K)
  REAL :: wa(lay)       ! Local layer water vapor mixing ratio (g/g)
  REAL :: oa(lay)       ! Local layer o3 mixing ratio (g/g)

!  real tsi(5)
!  real w(lay)
  REAL :: tab(3)
  REAL :: wk(lay,4)

  INTEGER :: i,j,k
  INTEGER :: IINT, ix
!  integer lun

!  real y
  REAL :: p
!
!----- ix=        1:trp; 2:mls; 3:mlw; 4:sas; 5:saw
!  data iflag/  1   ,  0   ,  0   ,  0   ,  0   /
!  data   tsi/ 300.0, 294.0, 272.2, 287.0, 257.1/
!
!-----------------------------------------------------------------------
!
!  Local data bank for pressure, temperature, moisture, and o3
!
!-----------------------------------------------------------------------
!
  DATA (pa(k),k=1,lay)/                                                 &
           .0003,    .0008,    .0011,    .0015,     .0021,              &
           .0029,    .0041,    .0058,    .0081,     .0113,              &
           .0158,    .0221,    .0310,    .0435,     .0609,              &
           .0855,    .1200,    .1700,    .2400,     .3350,              &
           .4650,    .6500,    .9150,   1.2850,    1.8000,              &
          2.5250,   3.5450,   4.9700,   6.9700,    9.7800,              &
         13.7150,  19.2350,  26.9850,  37.8550,   53.1000,              &
         73.8900,  97.6650, 121.4350, 145.2100,  168.9900,              &
        192.7650, 216.5400, 240.3150, 264.0900,  287.8650,              &
        311.6350, 335.4100, 359.1900, 382.9650,  406.7400,              &
        430.5150, 454.2850, 478.0600, 501.8350,  525.6100,              &
        549.3900, 573.1650, 596.9400, 620.7150,  644.4900,              &
        668.2650, 692.0350, 715.8100, 739.5850,  763.3600,              &
        787.1400, 810.9150, 834.6900, 858.4650,  882.2400,              &
        906.0150, 929.7850, 953.5600, 977.3350, 1001.1100/

  DATA (ta(k),k=1,lay)/                                                 &
        209.86, 210.20, 210.73, 211.27, 211.81,                         &
        212.35, 212.89, 213.44, 213.98, 214.53,                         &
        215.08, 215.62, 216.17, 216.74, 218.11,                         &
        223.20, 230.04, 237.14, 244.46, 252.00,                         &
        259.76, 267.70, 274.93, 274.60, 269.38,                         &
        262.94, 256.45, 250.12, 244.31, 238.96,                         &
        233.74, 228.69, 224.59, 221.75, 219.10,                         &
        216.64, 215.76, 215.75, 215.78, 216.22,                         &
        219.15, 223.79, 228.29, 232.45, 236.33,                         &
        239.92, 243.32, 246.53, 249.56, 252.43,                         &
        255.14, 257.69, 260.11, 262.39, 264.57,                         &
        266.66, 268.67, 270.60, 272.48, 274.29,                         &
        276.05, 277.75, 279.41, 281.02, 282.59,                         &
        284.09, 285.53, 286.86, 288.06, 289.13,                         &
        290.11, 291.03, 291.91, 292.76, 293.59/


  DATA (wa(k),k=1,lay)/                                                 &
        0.400E-05, 0.400E-05, 0.400E-05, 0.400E-05, 0.400E-05,          &
        0.400E-05, 0.400E-05, 0.400E-05, 0.400E-05, 0.400E-05,          &
        0.400E-05, 0.400E-05, 0.400E-05, 0.400E-05, 0.400E-05,          &
        0.400E-05, 0.400E-05, 0.400E-05, 0.400E-05, 0.400E-05,          &
        0.400E-05, 0.400E-05, 0.400E-05, 0.400E-05, 0.400E-05,          &
        0.400E-05, 0.400E-05, 0.400E-05, 0.400E-05, 0.400E-05,          &
        0.400E-05, 0.400E-05, 0.400E-05, 0.400E-05, 0.400E-05,          &
        0.400E-05, 0.400E-05, 0.400E-05, 0.406E-05, 0.520E-05,          &
        0.115E-04, 0.275E-04, 0.572E-04, 0.107E-03, 0.166E-03,          &
        0.223E-03, 0.285E-03, 0.360E-03, 0.446E-03, 0.547E-03,          &
        0.655E-03, 0.767E-03, 0.890E-03, 0.103E-02, 0.118E-02,          &
        0.136E-02, 0.159E-02, 0.190E-02, 0.225E-02, 0.264E-02,          &
        0.306E-02, 0.351E-02, 0.399E-02, 0.450E-02, 0.504E-02,          &
        0.560E-02, 0.619E-02, 0.680E-02, 0.742E-02, 0.805E-02,          &
        0.869E-02, 0.935E-02, 0.100E-01, 0.107E-01, 0.113E-01/

  DATA (oa(k),k=1,lay)/                                                 &
       .643E-07, .202E-06, .246E-06, .290E-06, .334E-06,                &
       .378E-06, .422E-06, .467E-06, .512E-06, .557E-06,                &
       .603E-06, .648E-06, .694E-06, .740E-06, .793E-06,                &
       .101E-05, .131E-05, .164E-05, .198E-05, .234E-05,                &
       .272E-05, .312E-05, .359E-05, .465E-05, .590E-05,                &
       .765E-05, .910E-05, .960E-05, .994E-05, .101E-04,                &
       .990E-05, .853E-05, .710E-05, .576E-05, .423E-05,                &
       .260E-05, .152E-05, .102E-05, .786E-06, .598E-06,                &
       .448E-06, .352E-06, .302E-06, .252E-06, .212E-06,                &
       .193E-06, .176E-06, .160E-06, .147E-06, .137E-06,                &
       .127E-06, .118E-06, .109E-06, .103E-06, .975E-07,                &
       .924E-07, .883E-07, .846E-07, .810E-07, .778E-07,                &
       .749E-07, .721E-07, .694E-07, .671E-07, .648E-07,                &
       .626E-07, .607E-07, .593E-07, .579E-07, .565E-07,                &
       .552E-07, .540E-07, .528E-07, .517E-07, .505E-07/
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ix=1

  iop(1)=4
  iop(2)=4
  IINT=1
  itab(1)=1
  itab(2)=0
  itab(3)=0

  CALL coeff(lay,pa,oa,wa,iop,IINT,wk)

  DO k=1,np
    DO j=1,n
      DO i=1,m
        p = 0.5 * ( pl(i,j,k+1) + pl(i,j,k) )
        CALL terp1(lay,pa,oa,wa,p,IINT,tab,itab)
        ao(i,j,k)=tab(1)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE fito3
!
!##################################################################
!##################################################################
!######                                                      ######
!######                   FUNCTION AERODEN                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

  FUNCTION aeroden(p,psfc)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Assign aerosol optical density (per Pa)
!  For now we use an estimate based on the aerosol density used
!  in simple climate models, where 0.05 of optical depth is
!  assumed uniformly distributed over mass sfc to 800 mb and
!  0.025 is uniformly distributed over mass 800 mb to 225 mb.
!
!  Here we normalize that distribution according to the
!  surface pressure so that the mountains don't lose aerosols.
!
!  Aerosol optical depth is computed by calling routine as
!  layer depth (Pa) times this aerosol density.
!
!  AUTHOR:   Keith Brewster
!
!  MODIFICATION HISTORY
!
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  REAL :: p,psfc
  REAL :: aeroden
!
  IF (p < (.225*psfc)) THEN
    aeroden=0.
  ELSE IF (p < (.80*psfc)) THEN
    aeroden=0.025/(0.575*psfc)
  ELSE
    aeroden=0.050/(0.2*psfc)
  END IF
  RETURN
  END FUNCTION aeroden
!
!##################################################################
!##################################################################
!######                                                      ######
!######            FUNCTION RH_TO_CLDCV                      ######
!######                                                      ######
!##################################################################
!##################################################################
!

  FUNCTION rh_to_cldcv(rh,hgt)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Obtain first guess cloud cover field from relative humidity.
!
!
!  AUTHOR:   Jian Zhang
!  07/95
!
!  MODIFICATION HISTORY
!
!  04/08/97  J. Zhang
!            Added the empirical relationship between RH and
!            cloud cover used by Koch et al. (1997).
!  Reference:
!  Reference:
!  Koch, S.E., A. Aksakal, and J.T. McQueen, 1997:
!      The influence of mesoscale humidity and evapotranspiration
!      fields on a model forecast of a cold-frontal squall line.
!      Mon. Wea. Rev.,  Vol.125,  384-409
!  09/10/97  J. Zhang
!            Modified the empirical relationship between cloud
!            fraction and relative humidity from quadratic
!            to one-fourth-power.
!
!
!-----------------------------------------------------------------------
!
!  INPUT:
!    rh               ! relative humidity
!    hgt              ! height (AGL)
!
!  OUTPUT:
!    rh_to_cld_cv     ! cloud fractional cover value
!
!  LOCAL:
!    rh0              ! the critical RH value that seperate clear
                      ! air condition and cloudy condition
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: rh2cform
  PARAMETER (rh2cform=2)

  REAL :: rh,hgt,rh_to_cldcv
  REAL :: rh0

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF(rh2cform == 1) THEN
!
!-----------------------------------------------------------------------
!
!   A quadratic relationship between relative humidity and cloud
!   fractional cover.
!
!-----------------------------------------------------------------------
!
    IF (hgt < 600.0) THEN
      rh0 = 0.9
    ELSE IF (hgt < 1500.0) THEN
      rh0 = 0.8
    ELSE IF (hgt < 2500.0) THEN
      rh0 = 0.6
    ELSE
      rh0 = 0.5
    END IF

    IF (rh < rh0) THEN
      rh_to_cldcv = 0.0
    ELSE
      rh_to_cldcv = (rh - rh0)/(1.0 - rh0)
      rh_to_cldcv = rh_to_cldcv*rh_to_cldcv
    END IF

  ELSE IF(rh2cform == 2) THEN
!
!-----------------------------------------------------------------------
!
!   A quadratic relationship between relative humidity and cloud
!   fractional cover with fixed rh0=0.75
!
!-----------------------------------------------------------------------
!
    IF (rh < 0.75) THEN
      rh_to_cldcv = 0.0
    ELSE
      rh_to_cldcv = 16.*(rh - 0.75)*(rh - 0.75)
    END IF

  ELSE
!
!-----------------------------------------------------------------------
!
!   A modified version of the sqrt relationship between
!   relative humidity and cloud fractional cover used in Eta model.
!
!-----------------------------------------------------------------------
!
    IF (hgt < 600.) THEN
      rh0 = 0.8
    ELSE
      rh0 = 0.75
    END IF

    IF (rh < rh0) THEN
      rh_to_cldcv = 0.0
    ELSE
      rh_to_cldcv = 1.0 - SQRT((1.0 - rh)/(1.0 - rh0))
    END IF

  END IF

  RETURN
  END FUNCTION rh_to_cldcv
