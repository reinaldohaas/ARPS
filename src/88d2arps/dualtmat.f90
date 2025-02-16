MODULE TMATRIX
!
! Fortran module with specifications for the 
! T-matrix scattering coefficient arrays.  The coefficients
! are complex variables for forward (0 deg) and back (180 deg)
! scattering.  For convenience in some calculations the square
! of the magnitude of foward and backward components are also stored.
!
! For snow and hail the coefficients are dependent on the 
! relative water content which is binned into 21 (nwfrac)
! water fraction ranging from 0 to 100%.
!
! Keith Brewster, CAPS  July, 2008
!
  IMPLICIT NONE
  INTEGER, PARAMETER :: nfile = 11
  INTEGER, PARAMETER :: nsize = 100
  INTEGER, PARAMETER :: nwfrac = 21
  REAL, PARAMETER :: fw_step = 0.05
  LOGICAL :: InitTmat = .FALSE.
  REAL :: wfrac(nwfrac)
  REAL :: diamrain(nsize)
  REAL :: diamsnow(nsize)
  REAL :: diamhail(nsize)
  COMPLEX :: fa_rain_fw(nsize)
  COMPLEX :: fa_rain_bk(nsize)
  COMPLEX :: fb_rain_fw(nsize)
  COMPLEX :: fb_rain_bk(nsize)
  COMPLEX :: fa_snow_fw(nsize,nwfrac)
  COMPLEX :: fa_snow_bk(nsize,nwfrac)
  COMPLEX :: fb_snow_fw(nsize,nwfrac)
  COMPLEX :: fb_snow_bk(nsize,nwfrac)
  COMPLEX :: fa_hail_fw(nsize,nwfrac)
  COMPLEX :: fa_hail_bk(nsize,nwfrac)
  COMPLEX :: fb_hail_fw(nsize,nwfrac)
  COMPLEX :: fb_hail_bk(nsize,nwfrac)
  REAL :: fa_rain_fw2(nsize)
  REAL :: fa_rain_bk2(nsize)
  REAL :: fb_rain_fw2(nsize)
  REAL :: fb_rain_bk2(nsize)
  REAL :: fa_snow_fw2(nsize,nwfrac)
  REAL :: fa_snow_bk2(nsize,nwfrac)
  REAL :: fb_snow_fw2(nsize,nwfrac)
  REAL :: fb_snow_bk2(nsize,nwfrac)
  REAL :: fa_hail_fw2(nsize,nwfrac)
  REAL :: fa_hail_bk2(nsize,nwfrac)
  REAL :: fb_hail_fw2(nsize,nwfrac)
  REAL :: fb_hail_bk2(nsize,nwfrac)
END MODULE TMATRIX

SUBROUTINE tmdualref(nx,ny,nz,tmatrix_dir,dbz_out,rhoair,qr,qs,qh, &
                     zh,zv,zhv,zdr,rhohv,kdp)
!
! Keith Brewster, CAPS  July, 2008
!
  USE TMATRIX
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nx,ny,nz
  CHARACTER(LEN=256), INTENT(IN) :: tmatrix_dir
  INTEGER, INTENT(IN) :: dbz_out
  REAL, INTENT(IN) :: rhoair(nx,ny,nz)     ! kg/m3
  REAL, INTENT(IN) :: qr(nx,ny,nz)         ! kg/kg
  REAL, INTENT(IN) :: qs(nx,ny,nz)         ! kg/kg
  REAL, INTENT(IN) :: qh(nx,ny,nz)         ! kg/kg
  REAL, INTENT(INOUT) :: zh(nx,ny,nz)        ! dBZ
  REAL, INTENT(INOUT) :: zv(nx,ny,nz)        ! dBZ
  REAL, INTENT(INOUT) :: zhv(nx,ny,nz)       ! dBZ
  REAL, INTENT(INOUT) :: zdr(nx,ny,nz)       ! dBZ
  REAL, INTENT(INOUT) :: rhohv(nx,ny,nz)     ! dimensionless
  REAL, INTENT(INOUT) :: kdp(nx,ny,nz)       ! deg/km

  REAL, PARAMETER :: wavel = 0.03          ! m
  REAL, PARAMETER :: kw2 = 0.91
  REAL, PARAMETER :: phibar = 0.0          ! degrees
  REAL, PARAMETER :: n0rain = 8.0E06       ! m**-4
  REAL, PARAMETER :: n0snow = 3.0E06       ! m**-4 
  REAL, PARAMETER :: n0hail = 4.0E04       ! m**-4
  REAL, PARAMETER :: rhorain = 800.0       ! kg/m3
  REAL, PARAMETER :: rhosnow = 100.0       ! kg/m3
  REAL, PARAMETER :: rhohail = 913.0       ! kg/m3
  REAL, PARAMETER :: epsilon = 1.0E-10     ! kg/kg

  REAL, SAVE :: pi,deg2rad
  REAL, SAVE :: fwcst
  REAL, SAVE :: zhcst,zvcst,zhvcst,kdpcst,zhvmag2,zhvmag
  REAL, SAVE :: sgrain,sgsnow
  REAL, SAVE :: arain,brain,crain,ckrain,ddrain
  REAL, SAVE :: erain,frain,grain,lcrain
  REAL, SAVE :: asnow,bsnow,csnow,cksnow,ddsnow
  REAL, SAVE :: esnow,fsnow,gsnow,lcsnow
  REAL, SAVE :: ahail,bhail,chail,ckhail,ddhail
  REAL, SAVE :: ehail,fhail,ghail,lchail

  REAL :: arg2,arg3,qwf,qtot
  REAL :: zhsum,zvsum
  REAL :: zhsumr,zvsumr,kdpsumr
  REAL :: zhsums,zvsums,kdpsums
  REAL :: zhsumh,zvsumh,kdpsumh
  REAL :: sghail,cfsig
  REAL :: lamrain,nrain
  REAL :: lamsnow,nsnow
  REAL :: lamhail,nhail
  COMPLEX :: zhvsum,zhvsumr,zhvsums,zhvsumh
  INTEGER :: i,j,k
  INTEGER :: isize,iwf

  IF(.NOT. InitTmat) THEN
    pi=acos(-1.0)
    deg2rad=pi/180.
    fwcst=1.0/fw_step
    zhcst=4.0E11*((wavel/pi)**4)/kw2
    zvcst=zhcst
    zhvcst=zhcst
    kdpcst=180.*wavel/pi
    lcrain=(pi*rhorain*n0rain)
    lcsnow=(pi*rhosnow*n0snow)
    lchail=(pi*rhohail*n0hail)
    sgrain=0.
    sgsnow=20.*deg2rad
    IF(phibar == 0.) THEN
      arg2=exp(-2.0*sgrain*sgrain)
      arg3=exp(-8.0*sgrain*sgrain)
    ELSE
      arg2=exp(-2.0*sgrain*sgrain)*cos(2.0*phibar*deg2rad)
      arg3=exp(-8.0*sgrain*sgrain)*cos(4.0*phibar*deg2rad)
    END IF
    arain=0.125*(3.0+4.0*arg2+arg3)
    brain=0.125*(3.0-4.0*arg2+arg3)
    crain=0.125*(1.0-arg3)
    ckrain=arg2
    erain=0.125*(1.0-arg3)
    frain=0.125*(3.0+4.0*arg2+arg3)
    grain=0.125*(3.0-4.0*arg2+arg3)
    print *, '     rain A,B,C,Ck=',arain,brain,crain,ckrain
    print *, ' erain,frain,grain=',erain,frain,grain

    arg2=exp(-2.0*sgsnow*sgsnow)
    arg3=exp(-8.0*sgsnow*sgsnow)
    asnow=0.125*(3.0+4.0*arg2+arg3)
    bsnow=0.125*(3.0-4.0*arg2+arg3)
    csnow=0.125*(1.0-arg3)
    cksnow=arg2
    esnow=0.125*(1.0-arg3)
    fsnow=0.125*(3.0+4.0*arg2+arg3)
    gsnow=0.125*(3.0-4.0*arg2+arg3)
    print *, '     snow A,B,C,Ck=',asnow,bsnow,csnow,cksnow
    print *, ' esnow,fsnow,gsnow=',esnow,fsnow,gsnow

    print *, ' Calling fill tmatrix'
    CALL FILLTM(tmatrix_dir)
    print *, ' Back from fill tmatrix'

    ddrain=diamrain(2)-diamrain(1)
    ddsnow=diamsnow(2)-diamsnow(1)
    ddhail=diamhail(2)-diamhail(1)
    print *, ' ddrain= ',ddrain
    print *, ' ddsnow= ',ddsnow
    print *, ' ddhail= ',ddhail

    print *, ' Done initializing T-matrix constants'
    InitTmat=.TRUE.
  END IF

! WRITE(57,'(a)') '  qwf   iqwf    qr     qs     qh     zh     zv     zdr'

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1

        qwf=1.0
        qtot=qr(i,j,k)+qs(i,j,k)+qh(i,j,k)
        IF(qtot > epsilon) THEN
          qwf=qr(i,j,k)/qtot
        END IF
        iwf=1+NINT(qwf*fwcst)
!       print *, ' qwf = ',qwf,' iwf = ',iwf,' wfrac(iwf)=',wfrac(iwf)
!       print *, ' HERE 1 qr,qs,qh:',qr(i,j,k),qs(i,j,k),qh(i,j,k)

        zhsumr=0.
        zvsumr=0.
        zhvsumr=(0.,0.)
        kdpsumr=0.

        IF(qr(i,j,k) > epsilon) THEN
          lamrain=(lcrain/(rhoair(i,j,k)*qr(i,j,k)))**0.25
!         print *, ' lamrain= ',lamrain
          DO isize=1,nsize
            nrain=n0rain*exp(-lamrain*diamrain(isize))
!           print *, ' isize,diam,nrain: ',isize,(1000.*diamrain(isize)),nrain
            zhsumr=zhsumr+nrain*(arain*fa_rain_bk2(isize) + &
                                 brain*fb_rain_bk2(isize) + &
                                 crain*fa_rain_bk2(isize))
            zvsumr=zvsumr+nrain*(brain*fa_rain_bk2(isize) + &
                                 arain*fb_rain_bk2(isize) + &
                                 crain*fa_rain_bk2(isize))
            zhvsumr=zhvsumr+nrain*(erain*(fa_rain_bk2(isize)+fb_rain_bk2(isize)) + &
                                   frain*(fa_rain_bk(isize)*CONJG(fb_rain_bk(isize))) + &
                                   grain*(fb_rain_bk(isize)*CONJG(fa_rain_bk(isize))))
            kdpsumr=kdpsumr+nrain*(real(fa_rain_fw(isize)) -            &
                                   real(fb_rain_fw(isize)))
          END DO
          zhsumr=zhsumr*ddrain
          zvsumr=zvsumr*ddrain
          zhvsumr=zhvsumr*ddrain
          kdpsumr=kdpsumr*ckrain*ddrain
        END IF

        zhsums=0.
        zvsums=0.
        zhvsums=(0.,0.)
        kdpsums=0.

        IF(qs(i,j,k) > epsilon) THEN
          lamsnow=(lcsnow/(rhoair(i,j,k)*qs(i,j,k)))**0.25
!         print *, ' lamsnow= ',lamsnow
          DO isize=1,nsize
            nsnow=n0snow*exp(-lamsnow*diamsnow(isize))
            zhsums=zhsums+nsnow*                                      &
                 (asnow*fa_snow_bk2(isize,iwf) +                      &
                  bsnow*fb_snow_bk2(isize,iwf) +                      &
                  csnow*fa_snow_bk2(isize,iwf))
            zvsums=zvsums+nsnow*                                      &
                 (bsnow*fa_snow_bk2(isize,iwf) +                      &
                  asnow*fb_snow_bk2(isize,iwf) +                      &
                  csnow*fa_snow_bk2(isize,iwf))
            zhvsums=zhvsums+nsnow*(esnow*(fa_snow_bk2(isize,iwf)+fb_snow_bk2(isize,iwf)) + &
                                   fsnow*(fa_snow_bk(isize,iwf)*CONJG(fb_snow_bk(isize,iwf))) + &
                                   gsnow*(fb_snow_bk(isize,iwf)*CONJG(fa_snow_bk(isize,iwf))))
            kdpsums=kdpsums+nsnow*(real(fa_snow_fw(isize,iwf)) -   &
                                   real(fb_snow_fw(isize,iwf)))

          END DO
          zhsums=zhsums*ddsnow
          zvsums=zvsums*ddsnow
          zhvsums=zhvsums*ddsnow
          kdpsums=kdpsums*cksnow*ddsnow
        END IF

        zhsumh=0.
        zvsumh=0.
        zhvsumh=(0.,0.)
        kdpsumh=0.

        IF(qh(i,j,k) > epsilon) THEN
          lamhail=(lchail/(rhoair(i,j,k)*qh(i,j,k)))**0.25
!         print *, ' lamhail = ',lamhail
          IF( qh(i,j,k) < 0.0002 ) THEN ! small hail
            cfsig=0.8*qwf
          ELSE
            cfsig=4000.0*qh(i,j,k)*qwf
          END IF
          cfsig=min(1.0,cfsig)
          sghail=deg2rad*(60.0*(1.0-cfsig))
          IF(phibar == 0.0) THEN
            arg2=exp(-2.0*sghail*sghail)
            arg3=exp(-8.0*sghail*sghail)
          ELSE
            arg2=exp(-2.0*sghail*sghail)*cos(2.0*phibar*deg2rad)
            arg3=exp(-8.0*sghail*sghail)*cos(4.0*phibar*deg2rad)
          END IF
          ahail=0.125*(3.0+4.0*arg2+arg3)
          bhail=0.125*(3.0-4.0*arg2+arg3)
          chail=0.125*(1.0-arg3)
          ckhail=arg2
          ehail=0.125*(1.0-arg3)
          fhail=0.125*(3.0+4.0*arg2+arg3)
          ghail=0.125*(3.0-4.0*arg2+arg3)

          DO isize=1,nsize
            nhail=n0hail*exp(-lamhail*diamhail(isize))
!           print *, ' isize,diamhail,nhail=',isize,(1000.*diamhail(isize)),nhail
            zhsumh=zhsumh+nhail*                                 &
                 (ahail*fa_hail_bk2(isize,iwf) +                      &
                  bhail*fb_hail_bk2(isize,iwf) +                      &
                  chail*fa_hail_bk2(isize,iwf))
            zvsumh=zvsumh+nhail*                                 &
                 (bhail*fa_hail_bk2(isize,iwf) +                      &
                  ahail*fb_hail_bk2(isize,iwf) +                      &
                  chail*fa_hail_bk2(isize,iwf))
            zhvsumh=zhvsumh+nhail*(ehail*(fa_hail_bk2(isize,iwf)+fb_hail_bk2(isize,iwf)) + &
                                   fhail*(fa_hail_bk(isize,iwf)*CONJG(fb_hail_bk(isize,iwf))) + &
                                   ghail*(fb_hail_bk(isize,iwf)*CONJG(fa_hail_bk(isize,iwf))))
            kdpsumh=kdpsumh+nhail*(real(fa_hail_fw(isize,iwf)) - &
                                   real(fb_hail_fw(isize,iwf)))
          END DO
          zhsumh=zhsumh*ddhail
          zvsumh=zvsumh*ddhail
          zhvsumh=zhvsumh*ddhail
          kdpsumh=kdpsumh*ckhail*ddhail
        END IF

        zhsum=zhcst*(zhsumr+zhsums+zhsumh)
        zvsum=zvcst*(zvsumr+zvsums+zvsumh)
        zhvsum=zhvcst*(zhvsumr+zhvsums+zhvsumh)

        kdp(i,j,k)=kdpcst*(kdpsumr+kdpsums+kdpsumh)

        IF(dbz_out > 0) THEN
          IF(zhsum > epsilon) THEN
            zh(i,j,k)=10.0*LOG10(zhsum)
          ELSE
            zh(i,j,k)=-30.
          END IF

          IF(zvsum > epsilon) THEN
            zv(i,j,k)=10.0*LOG10(zvsum)
          ELSE
            zv(i,j,k)=-30.
          END IF

          IF(zhsum > epsilon .AND. zvsum > epsilon) THEN
            zdr(i,j,k)=10.0*LOG10(zhsum/zvsum)
          ELSE
            zdr(i,j,k)=0.
          END IF

          zhvmag2=(real(zhvsum)*real(zhvsum))+(aimag(zhvsum)*aimag(zhvsum))
          IF(zhvmag2 > epsilon ) THEN
            zhvmag=sqrt(zhvmag2)
            zhv(i,j,k)=10.0*LOG10(zhvmag)
          ELSE
            zhvmag=0.
            zhv(i,j,k)=0.
          END IF

          IF(zhsum > epsilon .AND. zvsum > epsilon) THEN
            rhohv(i,j,k)=zhvsum/sqrt(zhsum*zhvsum)
          ELSE IF(zhsum == 0. .AND. zvsum == 0.) THEN
            rhohv(i,j,k)=1.
          ELSE
            rhohv(i,j,k)=0.
          END IF

!         IF(zh(i,j,k) > 30. .OR. zv(i,j,k) > 30.) &
!           WRITE(57,'(f9.2,i4,6f9.2)') qwf,iwf,(1000.*qr(i,j,k)), &
!             (1000.*qs(i,j,k)),(1000.*qh(i,j,k)), &
!             zh(i,j,k),zv(i,j,k),zdr(i,j,k)

        ELSE   ! dbz_out=0

          IF(zhsum > epsilon) THEN
            zh(i,j,k)=zhsum
          ELSE
            zh(i,j,k)=0.
          END IF

          IF(zvsum > epsilon) THEN
            zv(i,j,k)=zvsum
          ELSE
            zv(i,j,k)=0.
          END IF

          IF(zhsum > epsilon .AND. zvsum > epsilon) THEN
            zdr(i,j,k)=zhsum/zvsum
          ELSE
            zdr(i,j,k)=0.
          END IF

          zhvmag2=(real(zhvsum)*real(zhvsum))+(aimag(zhvsum)*aimag(zhvsum))
          IF(zhvmag2 > epsilon ) THEN
            zhv(i,j,k)=sqrt(zhvmag2)
          ELSE
            zhv(i,j,k)=0.
          END IF

          IF(zhsum > epsilon .AND. zvsum > epsilon) THEN
            rhohv(i,j,k)=zhvsum/sqrt(zhsum*zhvsum)
          ELSE IF(zhsum == 0. .AND. zvsum == 0.) THEN
            rhohv(i,j,k)=1.
          ELSE
            rhohv(i,j,k)=0.
          END IF

        END IF  ! dbz_out

      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE tmdualref

SUBROUTINE FILLTM (tmatrix_dir)
!
! Keith Brewster, CAPS  July, 2008
!
  USE TMATRIX
  IMPLICIT NONE
  CHARACTER(LEN=256), INTENT(IN) :: tmatrix_dir
!
!  Misc local variables
!
  INTEGER, PARAMETER :: iunit=51
  INTEGER :: iwf,nwf,isize
  REAL :: fa_real,fa_imag,fb_real,fb_imag
  CHARACTER(LEN=1) :: dummy
  CHARACTER(LEN=64) :: fwname
  CHARACTER(LEN=320) :: fwfile
!
!  Read rain scatter, 1 file
!
  nwf=100
  WRITE(fwname,'(a,i3.3,a)') 'SCTT_RAIN_fw',nwf,'.dat'
  WRITE(fwfile,'(a,a,a)') TRIM(tmatrix_dir),'/',TRIM(fwname)
  print *, ' Opening :',TRIM(fwfile)
  OPEN(iunit,file=fwfile,status='old',form='formatted')
  READ(iunit,'(a1)') dummy
  fa_rain_bk=(0.,0.)
  fb_rain_bk=(0.,0.)
  fa_rain_fw=(0.,0.)
  fb_rain_fw=(0.,0.)
  DO isize=1,nsize
    READ(iunit,'(f5.2,8g13.5)')  diamrain(isize),                       &
                 fa_rain_bk(isize),fb_rain_bk(isize),                   &
                 fa_rain_fw(isize),fb_rain_fw(isize)

    fa_real=real(fa_rain_bk(isize))
    fa_imag=aimag(fa_rain_bk(isize))
    fa_rain_bk2(isize)=(fa_real*fa_real + fa_imag*fa_imag)

    fb_real=real(fb_rain_bk(isize))
    fb_imag=aimag(fb_rain_bk(isize))
    fb_rain_bk2(isize)=(fb_real*fb_real + fb_imag*fb_imag)

    fa_real=real(fa_rain_fw(isize))
    fa_imag=aimag(fa_rain_fw(isize))
    fa_rain_fw2(isize)=(fa_real*fa_real + fa_imag*fa_imag)

    fb_real=real(fb_rain_fw(isize))
    fb_imag=aimag(fb_rain_fw(isize))
    fb_rain_fw2(isize)=(fb_real*fb_real + fb_imag*fb_imag)
  END DO
  CLOSE(iunit)
!
!  Read snow scatter, nwfrac files
!
  fa_snow_bk=(0.,0.)
  fb_snow_bk=(0.,0.)
  fa_snow_fw=(0.,0.)
  fb_snow_fw=(0.,0.)
  DO iwf=1,nwfrac
    wfrac(iwf)=(iwf-1)*fw_step
    nwf=NINT(wfrac(iwf)*100.0)
    WRITE(fwname,'(a,i3.3,a)') 'SCTT_SNOW_fw',nwf,'.dat'
    WRITE(fwfile,'(a,a,a)') TRIM(tmatrix_dir),'/',TRIM(fwname)
    print *, ' Opening :',TRIM(fwfile)
    OPEN(iunit,file=fwfile,status='old',form='formatted')
    READ(iunit,'(a1)') dummy
    DO isize=1,nsize
      READ(iunit,'(f5.2,8g13.5)')  diamsnow(isize),                       &
                 fa_snow_bk(isize,iwf),fb_snow_bk(isize,iwf),           &
                 fa_snow_fw(isize,iwf),fb_snow_fw(isize,iwf)

      fa_real=real(fa_snow_bk(isize,iwf))
      fa_imag=aimag(fa_snow_bk(isize,iwf))
      fa_snow_bk2(isize,iwf)=(fa_real*fa_real + fa_imag*fa_imag)

      fb_real=real(fb_snow_bk(isize,iwf))
      fb_imag=aimag(fb_snow_bk(isize,iwf))
      fb_snow_bk2(isize,iwf)=(fb_real*fb_real + fb_imag*fb_imag)

      fa_real=real(fa_snow_fw(isize,iwf))
      fa_imag=aimag(fa_snow_fw(isize,iwf))
      fa_snow_fw2(isize,iwf)=(fa_real*fa_real + fa_imag*fa_imag)

      fb_real=real(fb_snow_fw(isize,iwf))
      fb_imag=aimag(fb_snow_fw(isize,iwf))
      fb_snow_fw2(isize,iwf)=(fb_real*fb_real + fb_imag*fb_imag)
    END DO
    CLOSE(iunit)
  END DO
!
!  Read hail scatter, nwfrac files
!
  DO iwf=1,nwfrac
    nwf=NINT(wfrac(iwf)*100.0)
    WRITE(fwname,'(a,i3.3,a)') 'SCTT_HAIL_fw',nwf,'.dat'
    WRITE(fwfile,'(a,a,a)') TRIM(tmatrix_dir),'/',TRIM(fwname)
    print *, ' Opening :',TRIM(fwfile)
    OPEN(iunit,file=fwfile,status='old',form='formatted')
    READ(iunit,'(a1)') dummy
    DO isize=1,nsize
      READ(iunit,'(f5.2,8g13.5)')  diamhail(isize),                       &
                 fa_hail_bk(isize,iwf),fb_hail_bk(isize,iwf),           &
                 fa_hail_fw(isize,iwf),fb_hail_fw(isize,iwf)

      fa_real=real(fa_hail_bk(isize,iwf))
      fa_imag=aimag(fa_hail_bk(isize,iwf))
      fa_hail_bk2(isize,iwf)=(fa_real*fa_real + fa_imag*fa_imag)

      fb_real=real(fb_hail_bk(isize,iwf))
      fb_imag=aimag(fb_hail_bk(isize,iwf))
      fb_hail_bk2(isize,iwf)=(fb_real*fb_real + fb_imag*fb_imag)

      fa_real=real(fa_hail_fw(isize,iwf))
      fa_imag=aimag(fa_hail_fw(isize,iwf))
      fa_hail_fw2(isize,iwf)=(fa_real*fa_real + fa_imag*fa_imag)

      fb_real=real(fb_hail_fw(isize,iwf))
      fb_imag=aimag(fb_hail_fw(isize,iwf))
      fb_hail_fw2(isize,iwf)=(fb_real*fb_real + fb_imag*fb_imag)
    END DO
    CLOSE(iunit)
  END DO
!
!  Convert size to m from mm
!
  DO isize=1,nsize
    diamrain(isize)=0.001*diamrain(isize)
    diamsnow(isize)=0.001*diamsnow(isize)
    diamhail(isize)=0.001*diamhail(isize)
  END DO

END SUBROUTINE FILLTM
