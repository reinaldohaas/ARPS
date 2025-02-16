MODULE module_interpolation
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! This module is for intermediate interpolation variables
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE module_wrfgrid_constants
  USE wrf_llxy_module

  REAL, TARGET :: psnd(lvlprof),   zsnd(lvlprof),   tsnd(lvlprof),      &
          ptsnd(lvlprof),  rhssnd(lvlprof), qvsnd(lvlprof),             &
          rhosnd(lvlprof), usnd(lvlprof),   vsnd(lvlprof),              &
          dumsnd(lvlprof), plsnd(lvlprof)

  REAL, ALLOCATABLE, TARGET :: iur(:,:), jur(:,:)    ! relative to the Full E-grid
  REAL, ALLOCATABLE, TARGET :: ivr(:,:), jvr(:,:)
  REAL, ALLOCATABLE, TARGET :: isr(:,:), jsr(:,:)

  REAL, ALLOCATABLE, TARGET :: nmmi(:,:), nmmj(:,:)  ! Relative to the NMM E-grid

  REAL, ALLOCATABLE, TARGET :: zp_tmp(:,:,:), zps_tmp(:,:,:)
  REAL, POINTER             :: zpv_tmp(:,:,:)
  REAL, ALLOCATABLE, TARGET :: zpsoil_tmp(:)

  REAL, ALLOCATABLE :: workarr(:,:,:)

  INTEGER :: iEms, iEme         ! Full E-grid size in X direction

  CONTAINS

  SUBROUTINE set_interpolation_arrays( ims,ime,jms,jme, kms,kme,        &
             nx, ny, arpsxlat, arpsxlon, arpsylat, arpsylon,            &
             arpsslat, arpsslon,istatus )

!#######################################################################

    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: ims, ime, jms, jme, kms, kme
    INTEGER, INTENT(IN)  :: nx, ny
    REAL,    INTENT(IN)  :: arpsxlat(nx,ny), arpsxlon(nx,ny)
    REAL,    INTENT(IN)  :: arpsylat(nx,ny), arpsylon(nx,ny)
    REAL,    INTENT(IN)  :: arpsslat(nx,ny), arpsslon(nx,ny)

    INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

    INCLUDE 'mp.inc'

    INTEGER :: i,j
    !, ii, jj
!    REAL    :: xmin, ymin, xmax, ymax

!    REAL, ALLOCATABLE :: glatH(:,:), glonH(:,:)
!    REAL, ALLOCATABLE :: glatV(:,:), glonV(:,:)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    iEms = 2*ims -1
    iEme = 2*ime

!-----------------------------------------------------------------------
!
! This block is because wps_lltoxy does not return REAL number.
! If WPS fixes the problem later, it will be not necessary any more.
!
!-----------------------------------------------------------------------

!    ALLOCATE(glatH(ims:ime,jms:jme), STAT = istatus)
!    ALLOCATE(glonH(ims:ime,jms:jme), STAT = istatus)
!    ALLOCATE(glatV(ims:ime,jms:jme), STAT = istatus)
!    ALLOCATE(glonV(ims:ime,jms:jme), STAT = istatus)
!
!    DO j = jms,jme
!      DO i = ims,ime
!        CALL wps_xytoll(REAL(i),REAL(j),glatH(i,j),glonH(i,j),HH)
!        CALL wps_xytoll(REAL(i),REAL(j),glatV(i,j),glonV(i,j),VV)
!      END DO
!    END DO

!------------------ U Stagger ------------------------------------------

    ALLOCATE(iur(nx,ny), STAT = istatus)
    ALLOCATE(jur(nx,ny), STAT = istatus)

    DO j = 1,ny
      DO i = 1,nx
        CALL wps_lltoxy(arpsxlat(i,j),arpsxlon(i,j),iur(i,j),jur(i,j),VV)

        !write(0,*) 'U ',i,j,iur(i,j),jur(i,j)

        IF (iur(i,j) < iEms+1 .OR. iur(i,j) > iEme-1) THEN
          WRITE(6,'(1x,a,/,3(a,I6),a,F10.3)')       &
            'ARPS grid exceeds NMM grid. Please check the configuration or ' &
          //'increase number of NMM grid fake zone in namelist.',            &
            'Processor: ',myproc,' ims = ',iEms, ', ime = ',iEme,            &
            '; ARPS iur = ',iur(i,j)
          istatus = 11
          EXIT
        END IF

        IF (jur(i,j) < jms .OR. jur(i,j) > jme) THEN
          WRITE(6,'(1x,a,/,3(a,I6),a,F10.3)')        &
            'ARPS grid exceeds NMM grid. Please check the configuration or ' &
          //'increase number of NMM grid fake zone in namelist.',            &
            'Processor: ',myproc,' jms = ',jms, ', jme = ',jme,              &
            '; ARPS jur = ',jur(i,j)
          istatus = 12
          EXIT
        END IF

      END DO
    END DO

    CALL mpmaxi(istatus)
    IF (istatus > 0) THEN
      WRITE(6,'(1x,a,I2)') 'ERROR in set_interpolation_arrays, istatus = ',istatus
      CALL arpsstop('ARPS domain extends outside available external data.',1)
    END IF

!------------------ V Stagger ------------------------------------------

    ALLOCATE(ivr(nx,ny), STAT = istatus)
    ALLOCATE(jvr(nx,ny), STAT = istatus)

    DO j = 1,ny
       DO i = 1,nx
         CALL wps_lltoxy(arpsylat(i,j),arpsylon(i,j),ivr(i,j),jvr(i,j),VV)

         !write(0,*) 'V ',i,j,ivr(i,j),jvr(i,j)

         IF (ivr(i,j) < iEms+1 .OR. ivr(i,j) > iEme-1) THEN
           WRITE(6,'(1x,a,/,3(a,I6),a,F10.3)')        &
             'ARPS grid exceeds NMM grid. Please check the configuration or ' &
           //'increase number of NMM grid fake zone in namelist.',            &
             'Processor: ',myproc,' ims = ',iEms, ', ime = ',iEme,            &
             '; ARPS ivr = ',ivr(i,j)
           istatus = 21
           EXIT
         END IF

         IF (jvr(i,j) < jms .OR. jvr(i,j) > jme) THEN
           WRITE(6,'(1x,a,/,3(a,I6),a,F10.3)')        &
             'ARPS grid exceeds NMM grid. Please check the configuration or ' &
           //'increase number of NMM grid fake zone in namelist.',            &
             'Processor: ',myproc,' jms = ',jms, ', jme = ',jme,              &
             '; ARPS jvr = ',jvr(i,j)
           istatus = 22
           EXIT
         END IF

       END DO
    END DO

    CALL mpmaxi(istatus)
    IF (istatus > 0) THEN
      WRITE(6,'(1x,a,I2)') 'ERROR in set_interpolation_arrays, istatus = ',istatus
      CALL arpsstop('ARPS domain extends outside available external data.',1)
    END IF

!------------------ Scalar    ------------------------------------------

    ALLOCATE(isr(nx,ny), STAT = istatus)
    ALLOCATE(jsr(nx,ny), STAT = istatus)

    DO j = 1,ny
       DO i = 1,nx
         CALL wps_lltoxy(arpsslat(i,j),arpsslon(i,j),isr(i,j),jsr(i,j),HH)
         !write(0,*) 'H ',i,j,isr(i,j),jsr(i,j)

         IF (isr(i,j) < iEms+1 .OR. isr(i,j) > iEme-1) THEN
           WRITE(6,'(1x,a,/,3(a,I6),a,F10.3)')        &
             'ARPS grid exceeds NMM grid. Please check the configuration or ' &
           //'increase number of NMM grid fake zone in namelist.',            &
             'Processor: ',myproc,' ims = ',iEms, ', ime = ',iEme,            &
             '; ARPS isr = ',isr(i,j)
           istatus = 01
           EXIT
         END IF

         IF (jsr(i,j) < jms .OR. jsr(i,j) > jme) THEN
           WRITE(6,'(1x,a,/,3(a,I6),a,F10.3)')        &
             'ARPS grid exceeds NMM grid. Please check the configuration or ' &
           //'increase number of NMM grid fake zone in namelist.',            &
             'Processor: ',myproc,' jms = ',jms, ', jme = ',jme,              &
             '; ARPS jsr = ',jsr(i,j)
           istatus = 02
           EXIT
         END IF

      END DO
    END DO

    CALL mpmaxi(istatus)
    IF (istatus > 0) THEN
      WRITE(6,'(1x,a,I2)') 'ERROR in set_interpolation_arrays, istatus = ',istatus
      CALL arpsstop('ARPS domain extends outside available external data.',1)
    END IF

!------------------ NMM grid for Nearest Neighbor interpolation  -------

    ALLOCATE(nmmi(nx,ny), STAT = istatus)
    ALLOCATE(nmmj(nx,ny), STAT = istatus)

    DO j = 1,ny
      DO i = 1,nx
        CALL wps_lltoxy_nmm(arpsslat(i,j),arpsslon(i,j),nmmi(i,j),nmmj(i,j),HH)

        IF (nmmi(i,j) < ims .OR. nmmi(i,j) > ime) THEN
          WRITE(6,'(1x,a,/,3(a,I6),a,F10.3)')        &
            'ARPS grid exceeds NMM grid. Please check the configuration or ' &
          //'increase number of NMM grid fake zone in namelist.',            &
            'Processor: ',myproc,' ims = ',ims, ', ime = ',ime,              &
            '; ARPS nmmi = ',nmmi(i,j)
          istatus = 41
          EXIT
        END IF

        IF (nmmj(i,j) < jms .OR. nmmj(i,j) > jme) THEN
          WRITE(6,'(1x,a,/,3(a,I6),a,F10.3)')        &
            'ARPS grid exceeds NMM grid. Please check the configuration or ' &
          //'increase number of NMM grid fake zone in namelist.',            &
            'Processor: ',myproc,' jms = ',jms, ', jme = ',jme,              &
            '; ARPS nmmj = ',nmmj(i,j)
          istatus = 42
          EXIT
        END IF
      END DO
    END DO

    CALL mpmaxi(istatus)
    IF (istatus > 0) THEN
      WRITE(6,'(1x,a,I2)') 'ERROR in set_interpolation_arrays, istatus = ',istatus
      CALL arpsstop('ARPS domain extends outside available external data.',1)
    END IF

!------------------ Return  ------------------------------------------


    IF (.NOT. ALLOCATED(workarr))   &
      ALLOCATE(workarr(iEms:iEme,jms:jme,kms:kme), STAT = istatus)

!    DEALLOCATE(glatH, glatV, glonH, glonV )

    RETURN
  END SUBROUTINE set_interpolation_arrays

  SUBROUTINE set_interpolation_heights( nx,ny,intrpopt,ims,ime,jms,jme, &
                     kps,kpe,kpse,kzs,kzse,zp_ext,hgt_ext,zpsoil_ext, &
                     dbglvl,istatus )

!#######################################################################

    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: nx, ny
    INTEGER, INTENT(IN)  :: intrpopt
    INTEGER, INTENT(IN)  :: ims, ime, jms, jme
    INTEGER, INTENT(IN)  :: kps, kpe, kpse, kzs, kzse
    REAL,    INTENT(IN)  :: zp_ext(ims:ime,jms:jme,kps:kpse)
    REAL,    INTENT(IN)  :: hgt_ext(ims:ime,jms:jme,kps:kpe)
    REAL,    INTENT(IN)  :: zpsoil_ext(kzs:kzse)
    INTEGER, INTENT(IN)  :: dbglvl
    INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

    INTEGER :: i,j,k

    REAL    :: intrp_pnt_from_2d

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    ALLOCATE( zp_tmp    (nx,ny,kps:kpse),   STAT = istatus )
    ALLOCATE( zps_tmp   (nx,ny,kps:kpe),    STAT = istatus )
    ALLOCATE( zpsoil_tmp(kzs:kzse),         STAT = istatus )

!-----------------------------------------------------------------------
!
!  Process height data
!
!  Interpolate the external heights horizontally to the
!  ARPS x,y.  This is for scalar pts and W points.
!
!-----------------------------------------------------------------------
!
    CALl egridfill(zp_ext,ims,ime,jms,jme,kps,kpse,HH,istatus)

    DO k=kps,kpse         ! W points
      DO j = 1,ny
        DO i = 1,nx
          zp_tmp(i,j,k) = intrp_pnt_from_2d(intrpopt,isr(i,j),jsr(i,j), &
                            iEms,iEme,jms,jme,kps,kpse,k,workarr,       &
                            dbglvl,istatus)
        END DO
      END DO
    END DO

    CALL egridfill(hgt_ext,ims,ime,jms,jme,kps,kpe,HH,istatus)

    DO k=kps,kpe         ! scalar points
      DO j = 1,ny
        DO i = 1,nx
          zps_tmp(i,j,k) = intrp_pnt_from_2d(intrpopt,isr(i,j),jsr(i,j),&
                            iEms,iEme,jms,jme,kps,kpse,k,workarr,       &
                            dbglvl,istatus)
        END DO
      END DO
    END DO

    zpv_tmp => zps_tmp    ! Assume the same as scalar grid
    ! Acutally it should be interpolated according to E-grid structure

    ! NOTE:  zpsoil_ext is defined as soil depth!!!
    DO k=kzs,kzse
      zpsoil_tmp(k) = zpsoil_ext(k)
    END DO

    RETURN
  END SUBROUTINE set_interpolation_heights

  SUBROUTINE set_interpolation_basesnd(nx,ny,ims,ime,jms,jme,kps,kpe,   &
                      xscl,yscl,lat_ext,lon_ext,                        &
                      hgt_ext,pln_ext,t_ext,rhs_ext,u_ext,v_ext,        &
                      dbglvl, istatus )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Finds a mean sounding of all external data within the ARPS
!  domain.  This is used to define the model base state.
!  The data are interpolated to constant height levels and
!  averaged in the horizontal. The final sounding is hydrostatic.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  August, 1995  Replaces MNSOUND
!
!  MODIFICATION HISTORY:
!  Modified from extmnsnd in file src/arps/arpsplib3d.f90 for nmm2arps.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
    IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Sizing variables
!
!-----------------------------------------------------------------------
!
    INTEGER, INTENT(IN) :: nx,ny
    INTEGER, INTENT(IN) :: ims,ime,jms,jme
    INTEGER, INTENT(IN) :: kps,kpe

    INTEGER, INTENT(IN) :: dbglvl
    INTEGER, INTENT(OUT) :: istatus
!
!-----------------------------------------------------------------------
!
!  ARPS grid variables
!
!-----------------------------------------------------------------------
!
    REAL, INTENT(IN) :: xscl(nx)
    REAL, INTENT(IN) :: yscl(ny)
!
!-----------------------------------------------------------------------
!
!  External grid variables
!
!-----------------------------------------------------------------------
!
    REAL, INTENT(IN) :: lat_ext(ims:ime,jms:jme)
    REAL, INTENT(IN) :: lon_ext(ims:ime,jms:jme)

    REAL, INTENT(IN) :: hgt_ext(ims:ime,jms:jme,kps:kpe)   ! Height MSL
    REAL, INTENT(IN) :: t_ext  (ims:ime,jms:jme,kps:kpe)   ! Temperature (K)
    REAL, INTENT(IN) :: pln_ext(ims:ime,jms:jme,kps:kpe)   ! ln of total pressure
    REAL, INTENT(IN) :: rhs_ext(ims:ime,jms:jme,kps:kpe)   ! RHstar=SQRT(1.-RH)
    REAL, INTENT(IN) :: u_ext  (ims:ime,jms:jme,kps:kpe)   ! u wind component
    REAL, INTENT(IN) :: v_ext  (ims:ime,jms:jme,kps:kpe)   ! v wind component
  !
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
    INCLUDE 'phycst.inc'
    INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  Lapse rate and a constant for hydrostatic integration
!
!-----------------------------------------------------------------------
!
    REAL, PARAMETER :: gamma  = 0.0068         ! 6.8 degrees per km
    REAL, PARAMETER :: pconst = (g/rd)

    REAL, PARAMETER :: rhmin=0.05
    REAL, PARAMETER :: rhmax=1.0
!
!-----------------------------------------------------------------------
!
!  Output sounding variables, double precision version
!
!  Running a sum on lots of small numbers ends up with a precision
!  problem.  MPI and non-MPI runs will NEVER have the same answers.
!  Non-MPI runs will not even have the same answers in the numbers are
!  computed in a different order!
!
!-----------------------------------------------------------------------
!
    DOUBLE PRECISION :: ptsnd_dbl(lvlprof)
    DOUBLE PRECISION :: plsnd_dbl(lvlprof)
    DOUBLE PRECISION :: tsnd_dbl(lvlprof)
    DOUBLE PRECISION :: rhssnd_dbl(lvlprof)
    DOUBLE PRECISION :: usnd_dbl(lvlprof)
    DOUBLE PRECISION :: vsnd_dbl(lvlprof)

!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------

    INTEGER :: ibgn, jbgn, iend, jend

    REAL    :: deltaz

    INTEGER :: i,j,k,kext
    INTEGER :: knt,knttot,ktop,kbot
    REAL    :: wlow,whigh,tvbar,accept
    REAL    :: c1,c2,pres,qvsat,qvbot,qvtop,rh,repsilon
    REAL    :: temp, qvval

    REAL    :: xx, yy

!    REAL    :: imid,jmid
!    REAL    :: xmid,ymid,dist2,distmn

!-----------------------------------------------------------------------
!
!  Function f_qvsat and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
    REAL :: f_qvsat

!fpp$ expand (f_qvsat)
!!dir$ inline always f_qvsat
!*$*  inline routine (f_qvsat)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
    istatus = 0
!
!-----------------------------------------------------------------------
!
!  Set up the height levels on which to define the base-state.
!
!-----------------------------------------------------------------------
!
    deltaz = depthp/(lvlprof-1)

    IF ( myproc == 0 ) WRITE(6,'(/2x,2a,f7.2,a/)')                      &
      'The base state is formed from a mean sounding created',          &
      ' with a ',deltaz,' meter interval.'

    DO k=1,lvlprof
      zsnd(k) = (k-1)*deltaz
    END DO

    repsilon=0.01*(zsnd(2)-zsnd(1))
!
!-----------------------------------------------------------------------
!
!  Zero-out all sounding arrays
!  ptsnd array is used temporarily to store a count.
!
!-----------------------------------------------------------------------
!
    DO k=1,lvlprof
      ptsnd(k) = 0.
      plsnd(k) = 0.
      tsnd(k)  = 0.
      rhssnd(k)= 0.
      qvsnd(k) = 0.
      usnd(k)  = 0.
      vsnd(k)  = 0.

      ptsnd_dbl(k)  = 0.
      plsnd_dbl(k)  = 0.
      tsnd_dbl(k)   = 0.
      rhssnd_dbl(k) = 0.
      usnd_dbl(k)   = 0.
      vsnd_dbl(k)   = 0.
    END DO

!
!-----------------------------------------------------------------------
!
!  Find the external data grid location that is closest
!  to the middle of the grid.
!  For the case where no external grid point lies within the
!  grid this grid column will be used as the mean sounding.
!
!-----------------------------------------------------------------------
!
!    CALL wps_lltoxy(ctrlat,ctrlon,xmid,ymid,HH)
!
!    imid = NINT(xmid)
!    jmid = NINT(ymid)
!
!    iextmn=MIN(iextmn,imid)
!    iextmx=MAX(iextmx,imid)
!    jextmn=MIN(jextmn,jmid)
!    jextmx=MAX(jextmx,jmid)
!
!
!-----------------------------------------------------------------------
!
!  Look at all external points possibly within domain
!
!-----------------------------------------------------------------------
!
    knt=0
    knttot=0
    !
    !  Make sure we do the computations only once!
    !
    ibgn = 1
    jbgn = 1
    iend = nx
    jend = ny

    IF (loc_x > 1)         ibgn = 2
    IF (loc_x < nproc_x)   iend = nx-1
    IF (loc_y > 1)         jbgn = 2
    IF (loc_y < nproc_y)   jend = ny-1

    DO j = jms, jme
      DO i = ims,ime

        CALL lltoxy(1,1,lat_ext(i,j),lon_ext(i,j),xx,yy)
!
!-----------------------------------------------------------------------
!
!  Is this point within the ARPS domain?
!  Since the ARPS grid is Cartesian, need only compare
!  the external grid coordinates to the x and y limits
!  of the ARPS grid.
!
!-----------------------------------------------------------------------
!
        knttot=knttot+1
        IF( (xx >= xscl(ibgn) .AND. xx < xscl(iend)) .AND.              &
            (yy >= yscl(jbgn) .AND. yy < yscl(jend)) ) THEN

          knt=knt+1
!
!-----------------------------------------------------------------------
!
!  Interpolate external data in vertical onto profile
!  arrays.
!
!-----------------------------------------------------------------------
!
          DO k=1,lvlprof
            IF(zsnd(k) >= hgt_ext(i,j,kpe)) EXIT
            IF((hgt_ext(i,j,kps)-zsnd(k)) < repsilon) THEN
              DO kext= kps+1,kpe
                IF(hgt_ext(i,j,kext) >= zsnd(k)) EXIT
              END DO
              whigh=(zsnd(k)-hgt_ext(i,j,kext-1))/                        &
                    (hgt_ext(i,j,kext)-hgt_ext(i,j,kext-1))
              wlow=1.-whigh
              ptsnd_dbl(k)=ptsnd_dbl(k)+1.
              plsnd_dbl(k)=plsnd_dbl(k)+                                  &
                    whigh*pln_ext(i,j,kext)+wlow*pln_ext(i,j,kext-1)
              tsnd_dbl(k)=tsnd_dbl(k)+                                    &
                    whigh*t_ext(i,j,kext) + wlow*t_ext(i,j,kext-1)
              rhssnd_dbl(k)=rhssnd_dbl(k)+                                &
                    whigh*rhs_ext(i,j,kext)+wlow*rhs_ext(i,j,kext-1)
              usnd_dbl(k)=usnd_dbl(k)+                                    &
                    whigh*u_ext(i,j,kext) + wlow*u_ext(i,j,kext-1)
              vsnd_dbl(k)=vsnd_dbl(k)+                                    &
                    whigh*v_ext(i,j,kext) + wlow*v_ext(i,j,kext-1)
            END IF
          END DO
        END IF
      END DO
    END DO

    !
    ! mpi sum  ptsnd(k), knt, knttot
    !          plsnd(k),tsnd(k),rhssnd(k),usnd(k),vsnd(k)
    !
    IF(mp_opt > 0) THEN

      CALL mptotali(knt)
      CALL mptotali(knttot)

      CALL mpsumdp(ptsnd_dbl, lvlprof)
      CALL mpsumdp(plsnd_dbl, lvlprof)
      CALL mpsumdp(tsnd_dbl,  lvlprof)
      CALL mpsumdp(rhssnd_dbl,lvlprof)
      CALL mpsumdp(usnd_dbl,  lvlprof)
      CALL mpsumdp(vsnd_dbl,  lvlprof)
    END IF

    !
    !  Return to single precision.
    !
    DO k=1,lvlprof
      ptsnd(k)  = ptsnd_dbl(k)
      plsnd(k)  = plsnd_dbl(k)
      tsnd(k)   = tsnd_dbl(k)
      rhssnd(k) = rhssnd_dbl(k)
      usnd(k)   = usnd_dbl(k)
      vsnd(k)   = vsnd_dbl(k)
    END DO
!
! NOTE: plsnd, tsnd, etc should be exactly the same as
!       non-mpi mode because of the problem variables were made
!       double precision, which get around the round-off error.
!       The variables summed had sufficient precision.
!
!       "knttot" will be different in mpi and non-mpi modes, however, this
!       is totally harmless, as the last use is the WRITE below.  The numbers
!       can be made to be identical, however, doing it "the right way" makes
!       the run time longer for something that isn't of any value.
!


    IF(myproc == 0) WRITE(6,'(/a,i10,a,/12x,a,i10,a/)')                 &
         '  extmnsnd found ',knt,' points within ARPS domain',          &
         '  of ',knttot,' points in external domain checked.'

!
!-----------------------------------------------------------------------
!
!  Find lowest height with data
!
!-----------------------------------------------------------------------
!
    IF(myproc == 0) WRITE(6,'(a)') '  Finding range of mean sounding data ...'
    accept=0.3*knt
    DO k=1,lvlprof-1
      IF(myproc == 0) WRITE(6,'(3x,a,f10.2,a,f10.0,a,f10.0)')           &
                ' z = ',zsnd(k),' knt = ',ptsnd(k),' accept = ',accept
      IF(ptsnd(k) > accept) EXIT
    END DO
    kbot=k
!
!-----------------------------------------------------------------------
!
!  Find highest height with data
!
!-----------------------------------------------------------------------
!
    DO k=lvlprof,2,-1
      IF(myproc == 0) WRITE(6,'(3x,a,f10.2,a,f10.0,a,f10.0)')           &
               ' z = ',zsnd(k),' knt = ',ptsnd(k),' accept = ',accept
      IF(ptsnd(k) > accept) EXIT
    END DO
    ktop=k

    IF(myproc == 0) WRITE(6,'(3x,a,f10.2,a,f10.2,a)')                   &
                 ' Height of external data for mean spans from ',       &
                   zsnd(kbot),' to ',zsnd(ktop),' meters.'
!
!-----------------------------------------------------------------------
!
!  Divide through to find average.  We return to single precision here.
!
!-----------------------------------------------------------------------
!
    DO k=kbot,ktop
      plsnd(k)  = plsnd(k)/ptsnd(k)
      tsnd(k)   = tsnd(k)/ptsnd(k)
      rhssnd(k) = rhssnd(k)/ptsnd(k)
      usnd(k)   = usnd(k)/ptsnd(k)
      vsnd(k)   = vsnd(k)/ptsnd(k)
    END DO
!
!-----------------------------------------------------------------------
!
!  Set variables "below-ground"
!  Use a constant lapse rate, gamma.
!  plsnd is a sort of first-guess log(pressure), needed for qv
!  calculation from rhstar.
!
!-----------------------------------------------------------------------
!

    pres = EXP(plsnd(kbot))
    rh=MAX(rhmin,rhmax-(rhssnd(kbot)*rhssnd(kbot)))
    qvsat = f_qvsat( pres, tsnd(kbot) )
    qvbot=rh*qvsat
    c1=g/(rd*gamma)
    DO k=kbot-1,1,-1
      tsnd(k)=tsnd(kbot)-gamma*(zsnd(k)-zsnd(kbot))
      plsnd(k)=plsnd(kbot)+                                             &
          c1*ALOG((tsnd(kbot)-gamma*(zsnd(k)-zsnd(kbot)))/tsnd(kbot))
      psnd(k) = EXP(plsnd(k))
      qvsat = f_qvsat( psnd(k), tsnd(k) )
      rh=qvbot/qvsat
      rhssnd(k)=SQRT(MAX(0.,(rhmax-rh)))
      usnd(k)=usnd(kbot)
      vsnd(k)=vsnd(kbot)
    END DO
!
!-----------------------------------------------------------------------
!
!  Set variables "above-top"
!  Use a constant temperature. We're assuming stratosphere here.
!
!-----------------------------------------------------------------------
!

    pres = EXP(plsnd(ktop))
    rh=MAX(rhmin,rhmax-(rhssnd(ktop)*rhssnd(ktop)))
    qvsat = f_qvsat( pres, tsnd(ktop) )
    qvtop=rh*qvsat
    c2=g/rd
    DO k=ktop+1,lvlprof
      tsnd(k)=tsnd(ktop)
      plsnd(k)=plsnd(ktop)-c2*(zsnd(k)-zsnd(ktop))/tsnd(ktop)
      psnd(k) = EXP(plsnd(k))
      qvsat = f_qvsat( psnd(k), tsnd(k) )
      rh=qvtop/qvsat
      rhssnd(k)=SQRT(MAX(0.,(rhmax-rh)))
      usnd(k)=usnd(ktop)
      vsnd(k)=vsnd(ktop)
    END DO
!
!-----------------------------------------------------------------------
!
!  Calculate qv profile from RH-star.
!  Temporarily use rhosnd to store virtual temperature
!
!-----------------------------------------------------------------------
!

    DO k=1,lvlprof
      pres = EXP(plsnd(k))
      rh=MAX(rhmin,rhmax-(rhssnd(k)*rhssnd(k)))
      qvsat = f_qvsat( pres, tsnd(k) )
      qvsnd(k)=rh*qvsat
      rhosnd(k) = tsnd(k)*(1.0+rvdrd*qvsnd(k))/(1.0+qvsnd(k))
    END DO
!
!-----------------------------------------------------------------------
!
!  Make sure that the sounding is hydrostatic by integrating
!  from kbot. rhosnd is really virtual temperature here.
!
!-----------------------------------------------------------------------
!

    psnd(kbot)=EXP(plsnd(kbot))
    DO k=kbot-1,1,-1
      tvbar=0.5*(rhosnd(k+1)+rhosnd(k))
      psnd(k)=psnd(k+1)*EXP(pconst*(zsnd(k+1)-zsnd(k))/tvbar)
    END DO

    DO k=kbot+1,lvlprof
      tvbar=0.5*(rhosnd(k-1)+rhosnd(k))
      psnd(k)=psnd(k-1)*EXP(pconst*(zsnd(k-1)-zsnd(k))/tvbar)
    END DO
!
!-----------------------------------------------------------------------
!
!  Derived variable calculations
!  Compute density from virtual temperature.
!  Compute potential temperature.
!  Compute log of pressure.
!
!-----------------------------------------------------------------------

    DO k=1,lvlprof
      rhosnd(k)=psnd(k)/(rd*rhosnd(k))
      ptsnd(k)=tsnd(k)*((p0/psnd(k))**rddcp)
  !    pres=exp(plsnd(k))
  !    print *,'psnd old, psnd new: ',pres,psnd(k),(psnd(k)-pres)
      plsnd(k)=ALOG(psnd(k))
    END DO

!
!-----------------------------------------------------------------------
!
!  Print the mean sounding that was used in setting the
!  mean ARPS variables.
!
!-----------------------------------------------------------------------
!
    IF(dbglvl > 0 .AND. myproc == 0)  THEN

      WRITE(6,'(3x,a/,3x,a,a)')                                         &
          ' Mean sounding for ARPS base state variables',               &
          '  k     p(mb)     z(m)    pt(mb)    T(C)   qv(kg/kg) ',      &
          'RH %  u(m/s)    v(m/s)'

      DO k=lvlprof,1,-50
        pres = psnd(k)
        temp = ptsnd(k)*((pres/p0)**rddcp)
        rh=AMAX1(0.01,rhmax-(rhssnd(k)*rhssnd(k)))
        qvsat=f_qvsat( pres, temp )
        qvval=rh*qvsat
        WRITE(6,'(3x,i4,f9.1,f9.1,f9.1,f9.1,f9.5,f9.1,f9.1,f9.1)')      &
            k,(0.01*psnd(k)),zsnd(k),ptsnd(k),(temp-273.15),            &
            qvval,(100.*rh),usnd(k),vsnd(k)
      END DO

    END IF

    IF (dbglvl > 3) WRITE(6,'(3x,a)') '--- Done with interpolation base sounding.'

    RETURN
  END SUBROUTINE set_interpolation_basesnd

  SUBROUTINE deallocate_interpolation_arrays( istatus )

!#######################################################################

    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: istatus

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    DEALLOCATE( iur, jur )
    DEALLOCATE( ivr, jvr )
    DEALLOCATE( isr, jsr )

    DEALLOCATE( nmmi, nmmj )

    DEALLOCATE( zp_tmp )
    DEALLOCATE( zps_tmp )
    !DEALLOCATE( zpv_tmp )
    DEALLOCATE( zpsoil_tmp )

    DEALLOCATE( workarr )

    RETURN
  END SUBROUTINE deallocate_interpolation_arrays

  SUBROUTINE egridfill(arr,ims,ime,jms,jme,kms,kme,stagger, istatus)

    IMPLICIT NONE
!
!   This routine takes an array that contains data only on E-grid
!   staggered points (H or V), and spreads the data out to all
!   E-grid points (H and V).
!

    INTEGER, INTENT(IN)  :: ims,ime,jms,jme,kms,kme
    REAL,    INTENT(IN)  :: arr(ims:ime,jms:jme,kms:kme)
    INTEGER, INTENT(IN)  :: stagger
    INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

    INTEGER :: i,j,k, iec, jj, jje

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    jj  = MOD(jms,2)            ! start with even or odd rows?
    jje = MOD(jme,2)

!write(6,*) ims,ime,jms,jme,iEms,iEme
!call flush(6)

    IF (stagger == VV) THEN      ! E-grid velocity data

      DO k = kms, kme
        !
        ! First spread out the data to temporary array
        !
        DO j = jms+1-jj,jme, 2                    ! Odd rows
          DO i = iEme,iEms+1, -2                  ! iEms always odd, iEme always Even
            iec = i/2
            workarr(i,j,k)=arr(iec,j,k)
          END DO
        END DO

        DO j = jms+jj, jme, 2                     ! Even rows
          DO i = iEme-1, iEms,-2
             iec = (i+1)/2
             workarr(i,j,k) = arr(iec,j,k)
          END DO
        END DO
        !
        ! Next, fill in the corner points
        !
        IF (jj == 1) THEN
          workarr(iEms,jms,k)=.25*(3.*(arr(ims,jms,k)+arr(ims,jms+1,k))-  &
                                      (arr(ims+1,jms+1,k)+arr(ims,jms+2,k)))
        ELSE
          workarr(iEme,jms,k)=.25*(3.*(arr(ime,jms,k)+arr(ime,jms+1,k))-  &
                                      (arr(ime-1,jms+1,k)+arr(ime,jms+2,k)))

        END IF

        IF (jje == 0) THEN
          workarr(iEme,jme,k)=.25*(3.*(arr(ime,jme,k)+arr(ime,jme-1,k))-  &
                                      (arr(ime-1,jme-1,k)+arr(ime,jme-2,k)))
        ELSE
          workarr(iEms,jme,k)=.25*(3.*(arr(ims,jme,k)+arr(ims,jme-1,k))-  &
                                      (arr(ims+1,jme-1,k)+arr(ims,jme-2,k)))

        END IF

        DO j = jms+jj+1,jme-1,2  ! Fill in left and right side
           workarr(iEms,j,k)=.5*(workarr(iEms,j+1,k)+workarr(iEms,j-1,k))
           workarr(iEme,j,k)=.5*(workarr(iEme,j+1,k)+workarr(iEme,j-1,k))
        END DO

        DO i=iEms+jj+1,iEme-1,2  ! Fill in bototm sides.
           workarr(i,jms,k)=.5*(workarr(i+1,jms,k)+workarr(i-1,jms,k))
        END DO

        DO i=iEms+jje+1,iEme-1,2  ! Fill in top sides.
           workarr(i,jme,k)=.5*(workarr(i+1,jme,k)+workarr(i-1,jme,k))
        END DO

        DO j = jms+2-jj, jme-1, 2      ! Fill in interior even rows
          DO i = iEms+1, iEme-1, 2
            workarr(i,j,k) = 0.25*(workarr(i+1,j,k)+workarr(i-1,j,k)+   &
                                   workarr(i,j+1,k)+workarr(i,j-1,k) )
          END DO
        END DO

        DO j=jms+1+jj,jme-1,2  ! Fill in interior odd rows
          DO i= iEms+2, iEme-1, 2
            workarr(i,j,k) = 0.25*(workarr(i,j+1,k)+workarr(i,j-1,k)+   &
                                   workarr(i+1,j,k)+workarr(i-1,j,k))
          END DO
        END DO
      END DO ! end of k-loop

    ELSE IF ( stagger == HH ) THEN  ! E-grid mass data

      DO k = kms, kme
        !
        !  First spread out the data
        !
        DO j=jms+jj,jme,2           ! even rows
           DO i=iEme,iEms,-2
              iec=i/2
              workarr(i,j,k)=arr(iec,j,k)
           END DO
        END DO

        DO j = jms+1-jj,jme,2     ! odd rows
          DO i = iEme-1,iEms,-2
             iec=(i+1)/2
             workarr(i,j,k)=arr(iec,j,k)
          END DO
        END DO

        !
        ! Fill corners
        !
        IF (jj == 1) THEN
          workarr(iEme,jms,k) = 0.25*(3.*(arr(ime,jms,k)+arr(ime,jms+1,k)) -   &
                                         (arr(ime-1,jms+1,k)+arr(ime,jms+2,k)) )
        ELSE
          workarr(iEms,jms,k) = 0.25*(3.*(arr(ims,jms,k)+arr(ims,jms+1,k)) -   &
                                         (arr(ims,jms+2,k)+arr(ims+1,jms+1,k)) )
        END IF

        IF (jje == 1) THEN
          workarr(iEme,jme,k) = 0.25*(3.*(arr(ime,jme,k)+arr(ime,jme-1,k)) -   &
                                         (arr(ime-1,jme-1,k)+arr(ime,jme-2,k)) )
        ELSE
          workarr(iEms,jme,k) = 0.25*(3.*(arr(ims,jme,k)+arr(ims,jme-1,k)) -   &
                                         (arr(ims,jme-2,k)+arr(ims+1,jme-1,k)) )
        END IF

        DO j=jms+2-jj,jme-1,2    ! Fill in left and right side
           workarr(iEms,j,k)=.5*(workarr(iEms,j+1,k)+workarr(iEms,j-1,k))
           workarr(iEme,j,k)=.5*(workarr(iEme,j+1,k)+workarr(iEme,j-1,k))
        END DO

        DO i=iEms+2-jj,iEme-1,2  ! Fill in bototm sides.
           workarr(i,jms,k)=.5*(workarr(i+1,jms,k)+workarr(i-1,jms,k))
        END DO

        DO i=iEms+2-jje,iEme-1,2  ! Fill in top sides.
           workarr(i,jme,k)=.5*(workarr(i+1,jme,k)+workarr(i-1,jme,k))
        END DO

        DO j=jms+2-jj,jme-1,2  ! Fill in interior even rows
          DO i=iEms+2,iEme-1,2
             workarr(i,j,k)=.25*(workarr(i,j+1,k)+workarr(i,j-1,k)+     &
                                 workarr(i+1,j,k)+workarr(i-1,j,k))
          END DO
        END DO

        DO j=jms+1+jj,jme-1,2  ! Fill in interior odd rows
          DO i=iEms+1,iEme-1,2
             workarr(i,j,k)=.25*(workarr(i,j+1,k)+workarr(i,j-1,k)+     &
                                 workarr(i+1,j,k)+workarr(i-1,j,k))
          END DO
        END DO
      END DO ! end of k-loop

    END IF

    RETURN
  END SUBROUTINE egridfill

END MODULE module_interpolation
