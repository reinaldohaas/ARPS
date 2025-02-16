MODULE module_input_arpsgrid
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! This module contains data for the ARPS grid and it also handles ARPS
! related IO (reads ARPS history files).
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  USE module_arpsgrid_constants

  TYPE type_arpsgrid
    INTEGER :: nx, ny, nz, nzsoil, nstyps
    INTEGER :: nxlg, nylg

    INTEGER :: nzsoilh   ! nzsoil in history file
                         ! to handle extra soil data sets
    REAL    :: dx, dy

    INTEGER :: mapproj
    REAL    :: ctrlat,  ctrlon, ctrx, ctry
    REAL    :: trulat1, trulat2, trulon

    INTEGER :: nscalar, nscalarq
    INTEGER :: P_QC, P_QR, P_QI, P_QS, P_QG, P_QH
    INTEGER :: P_NC, P_NR, P_NI, P_NS, P_NG, P_NH
    INTEGER ::       P_ZR, P_ZI, P_ZS, P_ZG, P_ZH

!    REAL, ALLOCATABLE :: xlat(:,:),ylat(:,:),slat(:,:)
!    REAL, ALLOCATABLE :: xlon(:,:),ylon(:,:),slon(:,:)

    REAL, ALLOCATABLE :: x(:), xs(:)
    REAL, ALLOCATABLE :: y(:), ys(:)
    REAL, ALLOCATABLE :: z(:)
    REAL, ALLOCATABLE :: zp(:,:,:), zpsoil(:,:,:), zps(:,:,:)

    REAL, ALLOCATABLE :: p (:,:,:)
    REAL, ALLOCATABLE :: pt(:,:,:)
    REAL, ALLOCATABLE :: u (:,:,:)
    REAL, ALLOCATABLE :: v (:,:,:)
    REAL, ALLOCATABLE :: w (:,:,:)
    REAL, ALLOCATABLE :: qv(:,:,:)
    REAL, POINTER     :: qscalar(:,:,:,:)
    ! POINTER to be used as a target for handling the situation when QH is missing.
    REAL, ALLOCATABLE :: tke(:,:,:)

    INTEGER, ALLOCATABLE :: soiltyp(:,:,:)
    REAL,    ALLOCATABLE :: stypfrct(:,:,:)
    INTEGER, ALLOCATABLE :: vegtyp(:,:)
    REAL,    ALLOCATABLE :: veg(:,:)
    REAL,    ALLOCATABLE :: tsoil(:,:,:,:)
    REAL,    ALLOCATABLE :: qsoil(:,:,:,:)
    REAL,    ALLOCATABLE :: wetcanp(:,:,:)
    REAL,    ALLOCATABLE :: snowdpth(:,:)
    REAL,    ALLOCATABLE :: psfc(:,:)

  END TYPE type_arpsgrid

  TYPE (type_arpsgrid) :: arpsgrid

  LOGICAL :: arps_allocated = .FALSE.
  LOGICAL :: arps_readed    = .FALSE.

  INTEGER, PRIVATE :: ntile_x, ntile_y    ! used within this module only

  PRIVATE :: integ_moist  ! private subroutine

  CONTAINS

!##################### module level subroutines ########################

  SUBROUTINE alloc_init_arpsgrid(arpsfile,hinfmt,iamroot,               &
                                 extsoil,soilfile,                      &
                                 nprocx_in,nprocy_in,istatus)

!#######################################################################
!
!  PURPOSE:
!
!    Allocate arrays in "arpsgrid" and initialize them to some values
!    (zero).
!
!#######################################################################

    IMPLICIT NONE
    CHARACTER(LEN=MAXFILELEN), INTENT(IN)  :: arpsfile
    INTEGER,                   INTENT(IN)  :: hinfmt
    LOGICAL,                   INTENT(IN)  :: IAMROOT
    LOGICAL,                   INTENT(IN)  :: extsoil
    CHARACTER(LEN=MAXFILELEN), INTENT(IN)  :: soilfile
    INTEGER,                   INTENT(IN)  :: nprocx_in, nprocy_in
    INTEGER,                   INTENT(OUT) :: istatus

    INCLUDE 'globcst.inc'
    INCLUDE 'mp.inc'

!-----------------------------------------------------------------------

    CHARACTER(LEN=MAXFILELEN) :: filename
    INTEGER                   :: lenf

    INTEGER                   :: nx, ny, nz, nzsoil, nstyps, nzsoilh

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    lenf     = LEN_TRIM(arpsfile)
    filename = arpsfile
    ntile_x = 1
    ntile_y = 1
    IF ( nprocx_in > 1 .OR. nprocy_in > 1 ) THEN
      CALL gtsplitfn(arpsfile,1,1,loc_x,loc_y,1,1,         &
                     0,0,1,2,filename,istatus)
      lenf = LEN_TRIM(filename)
      ntile_x = nprocx_in / nproc_x
      ntile_y = nprocy_in / nproc_y
    END IF


    IF (IAMROOT) THEN
      CALL get_dims_from_data(hinfmt,filename,                  &
                              nx,ny,nz,nzsoil,nstyps, istatus)

      IF( istatus /= 0 ) THEN
        PRINT*,'Problem occured when trying to get dimensions from data.'
        PRINT*,'Program stopped.'
        istatus = -1
      END IF

      IF (mp_opt > 0) THEN
        IF (ntile_x > 1 .OR. ntile_y > 1) THEN     ! read and join tiles
          readsplit(:) = 0
          nx = (nx-3)*ntile_x + 3
          ny = (ny-3)*ntile_y + 3
          readstride = max_fopen
        ELSE IF (nprocx_in > 1 .OR. nprocy_in > 1) THEN ! read patches
          readsplit(:) = 0
          readstride = max_fopen
        ELSE                                       ! Read domain data and split
          readsplit(:) = 1
          IF( MOD(nx-3,nproc_x) /= 0 .OR. MOD(ny-3,nproc_y) /= 0) THEN
            WRITE(6,'(a/,a/,4(a,i5))')      &
             'The specification of nproc_x or nproc_y is not matched with nx or ny.',&
             'nx-3 and ny-3 must be multiples of nproc_x and nproc_y respectively.', &
             'nx = ', nx, ' ny = ', ny, ' nproc_x = ',nproc_x, ' nproc_y = ',nproc_y
             istatus = -2
          ELSE
            nx = (nx-3) / nproc_x + 3
            ny = (ny-3) / nproc_y + 3
          END IF
        END IF
        readstride = nprocs
      ELSE IF (ntile_x > 1 .OR. ntile_y > 1) THEN     ! nompi read and join
        nx = (nx-3)*ntile_x + 3
        ny = (ny-3)*ntile_y + 3
      END IF

      nzsoilh = nzsoil
      IF (extsoil) THEN
        filename = soilfile
        IF ( nprocx_in > 1 .OR. nprocy_in > 1 ) THEN
          CALL gtsplitfn(soilfile,1,1,loc_x,loc_y,1,1,         &
                         0,0,1,2,filename,istatus)
          lenf = LEN_TRIM(filename)
        END IF

        CALL get_ext_nzsoil(hinfmt,filename,nzsoil,istatus)
      END IF

    END IF
    CALL mpupdatei(istatus,1)
    IF (istatus /= 0) THEN
      CALL arpsstop('ERROR: in alloc_init_arpsgrid',1)
    END IF
    CALL mpupdatei(nx,1)
    CALL mpupdatei(ny,1)
    CALL mpupdatei(nz,1)
    CALL mpupdatei(nzsoil,1)
    CALL mpupdatei(nzsoilh,1)
    CALL mpupdatei(nstyps,1)
    CALL mpupdatei(readstride,1)
    CALL mpupdatei(readsplit,FINDX_NUM)

    CALL mpupdatei(nscalar, 1)
    CALL mpupdatei(nscalarq,1)
    CALL mpupdatei(P_QC,1);
    CALL mpupdatei(P_QR,1)
    CALL mpupdatei(P_QI,1)
    CALL mpupdatei(P_QS,1)
    CALL mpupdatei(P_QG,1)
    CALL mpupdatei(P_QH,1)
    CALL mpupdatei(P_NC,1)
    CALL mpupdatei(P_NR,1)
    CALL mpupdatei(P_NI,1)
    CALL mpupdatei(P_NS,1)
    CALL mpupdatei(P_NG,1)
    CALL mpupdatei(P_NH,1)
    CALL mpupdatei(P_ZR,1)
    CALL mpupdatei(P_ZI,1)
    CALL mpupdatei(P_ZS,1)
    CALL mpupdatei(P_ZG,1)
    CALL mpupdatei(P_ZH,1)

    CALL mpupdatec(qnames,nscalar*40)
    CALL mpupdatec(qdescp,nscalar*40)

    arpsgrid%nx = nx
    arpsgrid%ny = ny
    arpsgrid%nz = nz
    arpsgrid%nzsoil = nzsoil
    arpsgrid%nzsoilh = nzsoilh
    arpsgrid%nstyps = nstyps
    arpsgrid%nxlg = (nx-3)*nproc_x + 3
    arpsgrid%nylg = (ny-3)*nproc_y + 3

    arpsgrid%nscalar  = nscalar
    arpsgrid%nscalarq = nscalarq
    arpsgrid%P_QC     = P_QC
    arpsgrid%P_QR     = P_QR
    arpsgrid%P_QI     = P_QI
    arpsgrid%P_QS     = P_QS
    arpsgrid%P_QG     = P_QG
    arpsgrid%P_QH     = P_QH
    arpsgrid%P_NC     = P_NC
    arpsgrid%P_NR     = P_NR
    arpsgrid%P_NI     = P_NI
    arpsgrid%P_NS     = P_NS
    arpsgrid%P_NG     = P_NG
    arpsgrid%P_NH     = P_NH
    arpsgrid%P_ZR     = P_ZR
    arpsgrid%P_ZI     = P_ZI
    arpsgrid%P_ZS     = P_ZS
    arpsgrid%P_ZG     = P_ZG
    arpsgrid%P_ZH     = P_ZH

!-----------------------------------------------------------------------
!
! Allocate arrays in the structure "arpsgrid"
!
!-----------------------------------------------------------------------

!    ALLOCATE( arpsgrid%xlat(nx,ny), STAT = istatus )
!    ALLOCATE( arpsgrid%ylat(nx,ny), STAT = istatus )
!    ALLOCATE( arpsgrid%slat(nx,ny), STAT = istatus )
!    ALLOCATE( arpsgrid%xlon(nx,ny), STAT = istatus )
!    ALLOCATE( arpsgrid%ylon(nx,ny), STAT = istatus )
!    ALLOCATE( arpsgrid%slon(nx,ny), STAT = istatus )

    ALLOCATE( arpsgrid%x(nx),  STAT = istatus )
    ALLOCATE( arpsgrid%y(ny),  STAT = istatus )
    ALLOCATE( arpsgrid%xs(nx), STAT = istatus )
    ALLOCATE( arpsgrid%ys(ny), STAT = istatus )
    ALLOCATE( arpsgrid%z(nz),  STAT = istatus )
    ALLOCATE( arpsgrid%zp (nx,ny,nz),        STAT = istatus)
    ALLOCATE( arpsgrid%zps(nx,ny,nz),        STAT = istatus)
    ALLOCATE( arpsgrid%zpsoil(nx,ny,nzsoil), STAT = istatus )

    ALLOCATE( arpsgrid%p (nx,ny,nz), STAT = istatus )
    ALLOCATE( arpsgrid%pt(nx,ny,nz), STAT = istatus )
    ALLOCATE( arpsgrid%u (nx,ny,nz), STAT = istatus )
    ALLOCATE( arpsgrid%v (nx,ny,nz), STAT = istatus )
    ALLOCATE( arpsgrid%w (nx,ny,nz), STAT = istatus )
    ALLOCATE( arpsgrid%qv(nx,ny,nz), STAT = istatus )
    ALLOCATE( arpsgrid%qscalar(nx,ny,nz,nscalar), STAT = istatus )
    arpsgrid%qscalar(:,:,:,:) = 0.0
    ALLOCATE( arpsgrid%tke(nx,ny,nz), STAT = istatus )
    arpsgrid%tke(:,:,:) = 0.0

    ALLOCATE( arpsgrid%soiltyp (nx,ny,nstyps), STAT = istatus )
    ALLOCATE( arpsgrid%stypfrct(nx,ny,nstyps), STAT = istatus )
    arpsgrid%soiltyp (:,:,:) = 0
    arpsgrid%stypfrct(:,:,:) = 0.0
    ALLOCATE( arpsgrid%vegtyp(nx,ny),   STAT = istatus )
    ALLOCATE( arpsgrid%veg   (nx,ny),   STAT = istatus )
    ALLOCATE( arpsgrid%tsoil(nx,ny,nzsoil,0:nstyps), STAT = istatus )
    ALLOCATE( arpsgrid%qsoil(nx,ny,nzsoil,0:nstyps), STAT = istatus )
    ALLOCATE( arpsgrid%wetcanp     (nx,ny,0:nstyps), STAT = istatus )
    ALLOCATE( arpsgrid%snowdpth(nx,ny), STAT = istatus )
    ALLOCATE( arpsgrid%psfc(nx,ny),     STAT = istatus )

    arps_allocated = .TRUE.

    RETURN
  END SUBROUTINE alloc_init_arpsgrid

  SUBROUTINE read_arpsgrid(grdbasfn,hisfile,hinfmt,                     &
                           extsoil,soilfile,                            &
                           iamroot,datestr,istatus)

!#######################################################################
!
!  PURPOSE:
!
!    READ in "arpsgrid" at one time level from ARPS history file
!
!#######################################################################

    IMPLICIT NONE
    CHARACTER(LEN=MAXFILELEN), INTENT(IN) :: grdbasfn
    CHARACTER(LEN=MAXFILELEN), INTENT(IN) :: hisfile
    INTEGER,                   INTENT(IN) :: hinfmt
    LOGICAL,                   INTENT(IN) :: extsoil
    CHARACTER(LEN=MAXFILELEN), INTENT(IN) :: soilfile
    LOGICAL,                   INTENT(IN) :: IAMROOT
    CHARACTER(LEN=19), INTENT(OUT) :: datestr
    INTEGER,           INTENT(OUT) :: istatus

    INCLUDE 'mp.inc'
    INCLUDE 'globcst.inc'
    INCLUDE 'grid.inc'

!-----------------------------------------------------------------------

    REAL    :: time
    INTEGER :: abstime
    INTEGER :: lyear, lmonth, lday, lhour, lminute, lsecond

    INTEGER :: nx, ny, nz, nzsoil, nstyps

    REAL, ALLOCATABLE :: kmh(:,:,:),kmv(:,:,:)
    REAL, ALLOCATABLE :: ubar(:,:,:),vbar(:,:,:),wbar(:,:,:),pbar(:,:,:), &
                         ptbar(:,:,:),rhobar(:,:,:),qvbar(:,:,:)
    REAL, ALLOCATABLE :: lai(:,:),roufns(:,:)
    REAL, ALLOCATABLE :: raing(:,:),rainc(:,:),prcrate(:,:,:)
    REAL, ALLOCATABLE :: radfrc(:,:,:), radsw(:,:),radswnet(:,:), radlwin(:,:)
    REAL, ALLOCATABLE :: rnflx(:,:), usflx(:,:),vsflx(:,:),ptsflx(:,:),qvsflx(:,:)
    REAL, ALLOCATABLE :: tem1(:,:,:),tem2(:,:,:),tem3(:,:,:)

    REAL :: swx, swy
    REAL :: lattru(2), lontru

    REAL, ALLOCATABLE :: xs(:), ys(:)

    INTEGER :: i, j, k

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    nx = arpsgrid%nx
    ny = arpsgrid%ny
    nz = arpsgrid%nz
    nzsoil = arpsgrid%nzsoil
    nstyps = arpsgrid%nstyps

!-----------------------------------------------------------------------

    ALLOCATE( kmh(nx,ny,nz), STAT = istatus )
    ALLOCATE( kmv(nx,ny,nz), STAT = istatus )
    ALLOCATE( ubar(nx,ny,nz), STAT = istatus )
    ALLOCATE( vbar(nx,ny,nz), STAT = istatus )
    ALLOCATE( wbar(nx,ny,nz), STAT = istatus )
    ALLOCATE( pbar(nx,ny,nz), STAT = istatus )
    ALLOCATE( ptbar(nx,ny,nz),  STAT = istatus )
    ALLOCATE( rhobar(nx,ny,nz), STAT = istatus )
    ALLOCATE( qvbar(nx,ny,nz),  STAT = istatus )
    ALLOCATE( lai(nx,ny),   STAT = istatus )
    ALLOCATE( roufns(nx,ny),STAT = istatus )
    ALLOCATE( raing(nx,ny), STAT = istatus )
    ALLOCATE( rainc(nx,ny), STAT = istatus )
    ALLOCATE( prcrate(nx,ny,8), STAT = istatus )
    ALLOCATE( radfrc(nx,ny,nz), STAT = istatus )
    ALLOCATE( radsw(nx,ny),     STAT = istatus )
    ALLOCATE( radswnet(nx,ny),  STAT = istatus )
    ALLOCATE( radlwin(nx,ny),   STAT = istatus )
    ALLOCATE( rnflx(nx,ny),  STAT = istatus )
    ALLOCATE( usflx(nx,ny),  STAT = istatus )
    ALLOCATE( vsflx(nx,ny),  STAT = istatus )
    ALLOCATE( ptsflx(nx,ny), STAT = istatus )
    ALLOCATE( qvsflx(nx,ny), STAT = istatus )
    ALLOCATE( tem1(nx,ny,nz), STAT = istatus )
    ALLOCATE( tem2(nx,ny,nz), STAT = istatus )
    ALLOCATE( tem3(nx,ny,nz), STAT = istatus )

!-----------------------------------------------------------------------

    CALL readarpsmp(ntile_x, ntile_y,1,1,1,hinfmt,                          &
                    grdbasfn,LEN_TRIM(grdbasfn),hisfile,LEN_TRIM(hisfile),  &
                    nx,ny,nz,arpsgrid%nzsoilh,nstyps,time,                  &
       arpsgrid%x,arpsgrid%y,arpsgrid%z,arpsgrid%zp,arpsgrid%zpsoil,        &
       arpsgrid%u,arpsgrid%v,arpsgrid%w,arpsgrid%pt,arpsgrid%p,arpsgrid%qv, &
       arpsgrid%qscalar,                                                    &
       arpsgrid%tke,kmh,kmv,ubar,vbar,wbar,ptbar,pbar,rhobar,qvbar,         &
                    arpsgrid%soiltyp,arpsgrid%stypfrct,arpsgrid%vegtyp,     &
                    lai,roufns,arpsgrid%veg,arpsgrid%tsoil,arpsgrid%qsoil,  &
                    arpsgrid%wetcanp,arpsgrid%snowdpth,raing,rainc,prcrate, &
                    radfrc, radsw,rnflx,radswnet, radlwin,                  &
                    usflx,vsflx,ptsflx,qvsflx,istatus,tem1,tem2,tem3)

    arpsgrid%mapproj = mapproj
    arpsgrid%ctrlat  = ctrlat
    arpsgrid%ctrlon  = ctrlon
    arpsgrid%trulat1 = trulat1
    arpsgrid%trulat2 = trulat2
    arpsgrid%trulon  = trulon
    arpsgrid%dx      = dx
    arpsgrid%dy      = dy

    arpsgrid%u  = arpsgrid%u  + ubar
    arpsgrid%v  = arpsgrid%v  + vbar
    arpsgrid%w  = arpsgrid%w  + wbar
    arpsgrid%pt = arpsgrid%pt + ptbar
    arpsgrid%p  = arpsgrid%p  + pbar
    arpsgrid%qv = arpsgrid%qv + qvbar

    CALL ctim2abss( year,month,day,hour,minute,second, abstime )
    abstime = abstime + INT(time)
    CALL abss2ctim( abstime, lyear, lmonth, lday, lhour, lminute, lsecond )

    WRITE(datestr,'(I4.4,5(a,I2.2))') lyear,'-',lmonth,'-',lday,'_',    &
                                      lhour,':',lminute,':',lsecond

    arps_readed = .TRUE.

!-----------------------------------------------------------------------

    lattru(1) = trulat1
    lattru(2) = trulat2
    lontru    = trulon
    CALL setmapr(mapproj,sclfct,lattru,lontru)

    CALL lltoxy( 1,1, ctrlat,ctrlon, arpsgrid%ctrx, arpsgrid%ctry )
    swx = arpsgrid%ctrx - 0.5*(arpsgrid%nxlg-3) * dx
    swy = arpsgrid%ctry - 0.5*(arpsgrid%nylg-3) * dy
    CALL setorig( 1, swx, swy)

    DO i = 1,nx-1
      arpsgrid%xs(i) = 0.5*( arpsgrid%x(i) + arpsgrid%x(i+1) )
    END DO
    arpsgrid%xs(nx) = arpsgrid%xs(nx-1) + dx

    DO j = 1,ny-1
      arpsgrid%ys(j) = 0.5*( arpsgrid%y(j) + arpsgrid%y(j+1) )
    END DO
    arpsgrid%ys(ny) = arpsgrid%ys(ny-1) + dy

    DO k = 1,nz-1
      DO j = 1,ny
        DO i = 1,nx
          arpsgrid%zps(i,j,k)= (arpsgrid%zp(i,j,k)+arpsgrid%zp(i,j,k+1))*0.5
        END DO
      END DO
    END DO
    DO j = 1,ny
      DO i = 1,nx
        arpsgrid%zps(i,j,nz) = arpsgrid%zps(i,j,nz-1)+                  &
                           (arpsgrid%zp(i,j,nz)-arpsgrid%zp(i,j,nz-1))
      END DO
    END DO

!    CALL xytoll( nx,ny,arpsgrid%x,        ys, arpsgrid%xlat, arpsgrid%xlon )
!    CALL xytoll( nx,ny,        xs,arpsgrid%y, arpsgrid%ylat, arpsgrid%ylon )
!    CALL xytoll( nx,ny,        xs,        ys, arpsgrid%slat, arpsgrid%slon )
!
!    DEALLOCATE( xs, ys )

!-----------------------------------------------------------------------
! Calculate surface pressure
! Note the physical ground and the vertical staggering
!-----------------------------------------------------------------------

    DO j = 1, ny
      DO i = 1, nx
        arpsgrid%psfc(i,j) = EXP( 0.5*( LOG(arpsgrid%p(i,j,1)) + LOG(arpsgrid%p(i,j,2)) ) )
      END DO
    END DO

!-----------------------------------------------------------------------
!
! Read soil variable from another set of ARPS history files
!
!-----------------------------------------------------------------------

    IF ( extsoil) THEN
      CALL get_ext_soil(ntile_x, ntile_y,1,hinfmt,soilfile,LEN_TRIM(soilfile),  &
                    nx,ny,nzsoil,nstyps,arpsgrid%zpsoil,                &
                    arpsgrid%tsoil,arpsgrid%qsoil,istatus)
    END IF

!-----------------------------------------------------------------------

    DEALLOCATE( kmh,kmv )
    DEALLOCATE( ubar,vbar,wbar,pbar,ptbar,rhobar,qvbar )
    DEALLOCATE( lai,roufns )
    DEALLOCATE( raing,rainc,prcrate )
    DEALLOCATE( radfrc, radsw,radswnet, radlwin )
    DEALLOCATE( rnflx, usflx,vsflx,ptsflx,qvsflx )
    DEALLOCATE( tem1,tem2,tem3 )

    RETURN
  END SUBROUTINE read_arpsgrid

  SUBROUTINE dealloc_arpsgrid(istatus)
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: istatus

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!    DEALLOCATE( arpsgrid%xlat, arpsgrid%ylat, arpsgrid%slat )
!    DEALLOCATE( arpsgrid%xlon, arpsgrid%ylon, arpsgrid%slon )

    DEALLOCATE( arpsgrid%x )
    DEALLOCATE( arpsgrid%y )
    DEALLOCATE( arpsgrid%xs )
    DEALLOCATE( arpsgrid%ys )
    DEALLOCATE( arpsgrid%z )
    DEALLOCATE( arpsgrid%zp  )
    DEALLOCATE( arpsgrid%zps )
    DEALLOCATE( arpsgrid%zpsoil )

    DEALLOCATE( arpsgrid%p )
    DEALLOCATE( arpsgrid%pt )
    DEALLOCATE( arpsgrid%u )
    DEALLOCATE( arpsgrid%v )
    DEALLOCATE( arpsgrid%w )
    DEALLOCATE( arpsgrid%qv )
    DEALLOCATE( arpsgrid%qscalar )
    DEALLOCATE( arpsgrid%tke )

    DEALLOCATE( arpsgrid%soiltyp )
    DEALLOCATE( arpsgrid%stypfrct )
    DEALLOCATE( arpsgrid%vegtyp )
    DEALLOCATE( arpsgrid%veg )
    DEALLOCATE( arpsgrid%tsoil )
    DEALLOCATE( arpsgrid%qsoil )
    DEALLOCATE( arpsgrid%wetcanp )
    DEALLOCATE( arpsgrid%snowdpth )
    DEALLOCATE( arpsgrid%psfc )

    arps_allocated = .FALSE.

    RETURN
  END SUBROUTINE dealloc_arpsgrid

!#################### arpsgrid subroutines  ############################

  SUBROUTINE arpsgrid_getpsl(psl,istatus)

    IMPLICIT NONE

    REAL,    INTENT(OUT) :: psl(:,:)    ! sea level pressure
    INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

    REAL, PARAMETER :: gamma=6.5,              & ! 6.5 K/km
                       ex1=0.1903643,          & ! R*gamma/g
                       ex2=5.2558774             ! g/R/gamma

    REAL, ALLOCATABLE :: t700(:,:), tz(:,:,:)
    REAL, ALLOCATABLE :: algpzc(:,:,:)

    INTEGER :: nx, ny, nz
    INTEGER :: i, j, k

    REAL    :: p00, t00
    REAL    :: zlevel

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    nx = arpsgrid%nx
    ny = arpsgrid%ny
    nz = arpsgrid%nz
!-----------------------------------------------------------------------
!
!  Calculate temperature (K) at ARPS grid points and 700mb
!  pressure level
!
!-----------------------------------------------------------------------
!
    ALLOCATE(tz  (nx,ny,nz), STAT = istatus)
    ALLOCATE(t700(nx,ny),    STAT = istatus)

!-----------------------------------------------------------------------
!
!  Calculate the temperature using Poisson's formula.
!
!-----------------------------------------------------------------------
!
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tz(i,j,k) = arpsgrid%pt(i,j,k)*((arpsgrid%p(i,j,k) / p0) ** rddcp)
        END DO
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Compute temperature at 700 hPa
!
!-----------------------------------------------------------------------
!
    ALLOCATE(algpzc(nx,ny,nz), STAT = istatus)

    DO k = 1, nz
      DO j = 1, ny
        DO i = 1, nx
          algpzc(i,j,k) = -ALOG(arpsgrid%p(i,j,k))
        END DO
      END DO
    END DO

    zlevel=-ALOG(70000.0)

    DO i=1,nx-1
      DO j=1,ny-1
        IF(zlevel <= algpzc(i,j,1)) THEN
          k = 1
        ELSE IF(zlevel >= algpzc(i,j,nz-1)) THEN
          k = nz-1
        ELSE
          DO k=2,nz-2
            IF(zlevel >= algpzc(i,j,k) .AND. zlevel < algpzc(i,j,k+1)) EXIT
          END DO
        END IF

        t700(i,j)=tz(i,j,k)+(tz(i,j,k+1)-tz(i,j,k))*                    &
                  (zlevel-algpzc(i,j,k))/(algpzc(i,j,k+1)-algpzc(i,j,k))

!-----------------------------------------------------------------------
!
!  If the data point is below the ground level, set the
!  data value to the missing value.
!
!-----------------------------------------------------------------------

        IF( zlevel < algpzc(i,j,2)   ) t700(i,j) = -9999.0
        IF( zlevel > algpzc(i,j,nz-1)) t700(i,j) = -9999.0

      END DO
    END DO

!----------------------------------------------------------------------
!
!  Calculate sea level pressure (Pa)
!  Reduction method: Benjamin and Miller: 1990, MWR, vol.118, No.10,
!                   Page: 2100-2101
!
!-----------------------------------------------------------------------
!

    DO i=1,nx-1
      DO j=1,ny-1
        p00 = arpsgrid%p(i,j,2)
        IF(p00 <= 70000.0) THEN
          t00=tz(i,j,2)
        ELSE
          t00 = t700(i,j)*(p00/70000.0)**ex1
        END IF
        psl(i,j) = p00*((t00+gamma*arpsgrid%zps(i,j,2)*0.001)/t00)**ex2
      END DO
    END DO

    DEALLOCATE( tz, t700 )
    DEALLOCATE( algpzc )

    RETURN
  END SUBROUTINE arpsgrid_getpsl

  SUBROUTINE arpsgrid_gettmp( airtmp, istatus)
    IMPLICIT NONE
    REAL,    INTENT(OUT) :: airtmp(arpsgrid%nx,arpsgrid%ny,arpsgrid%nz)
    INTEGER, INTENT(OUT) :: istatus

  !---------------------------------------------------------------------
    INTEGER :: i, j, k

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    istatus = 0

    DO k = 1, arpsgrid%nz
      DO j = 1, arpsgrid%ny
        DO i = 1, arpsgrid%nx
          airtmp(i,j,k) = arpsgrid%pt(i,j,k) * ((arpsgrid%p(i,j,k) / p0) ** rddcp)
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE arpsgrid_gettmp

  SUBROUTINE arpsgrid_getgpht( gpht, istatus)
    IMPLICIT NONE
    REAL,    INTENT(OUT) :: gpht(arpsgrid%nx,arpsgrid%ny,arpsgrid%nz)
    INTEGER, INTENT(OUT) :: istatus

  !---------------------------------------------------------------------
    INTEGER :: i, j, k

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    istatus = 0

    DO k = 1, arpsgrid%nz
      DO j = 1, arpsgrid%ny
        DO i = 1, arpsgrid%nx
          gpht(i,j,k) = arpsgrid%zps(i,j,k)  ! * g
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE arpsgrid_getgpht

  SUBROUTINE arpsgrid_getdwpt( t, dwpt, istatus)
    IMPLICIT NONE
    REAL,    INTENT(IN)  :: t(arpsgrid%nx,arpsgrid%ny,arpsgrid%nz)
    REAL,    INTENT(OUT) :: dwpt(arpsgrid%nx,arpsgrid%ny,arpsgrid%nz)
    INTEGER, INTENT(OUT) :: istatus

  !---------------------------------------------------------------------

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    istatus = 0

    CALL getdew(arpsgrid%nx,arpsgrid%ny,arpsgrid%nz,                    &
                1,arpsgrid%nx,1,arpsgrid%ny,1,arpsgrid%nz,              &
                arpsgrid%p, t, arpsgrid%qv, dwpt)

    RETURN
  END SUBROUTINE arpsgrid_getdwpt

  SUBROUTINE arpsgrid_getvpress( vpres, istatus)
  !#####################################################################
!   this function returns the saturation vapor pressure over
!   water (mb) given the temperature (celsius).
!   the algorithm is due to nordquist, w.s.,1973: "numerical approxima-
!   tions of selected meteorlolgical parameters for cloud physics prob-
!   lems," ecom-5475, atmospheric sciences laboratory, u.s. army
!   electronics command, white sands missile range, new mexico 88002.

    IMPLICIT NONE
    REAL,    INTENT(OUT) :: vpres(arpsgrid%nx,arpsgrid%ny,arpsgrid%nz)
    INTEGER, INTENT(OUT) :: istatus

  !---------------------------------------------------------------------

    INTEGER :: i,j,k
    REAL    :: tk, p1, p2, c1

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    DO k = 1, arpsgrid%nz
      DO j = 1, arpsgrid%ny
        DO i = 1, arpsgrid%nx

          vpres(i,j,k) = arpsgrid%p(i,j,k)*arpsgrid%qv(i,j,k)/(arpsgrid%qv(i,j,k)+0.62197)
          vpres(i,j,k) = MAX( 1.0e-6, vpres(i,j,k) )

        END DO
      END DO
    END DO

    !DO k = 1, arpsgrid%nz
    !  DO j = 1, arpsgrid%ny
    !    DO i = 1, arpsgrid%nx
    !      tk = td(i,j,k)
    !      p1 = 11.344-0.0303998*tk
    !      p2 = 3.49149-1302.8844/tk
    !      c1 = 23.832241-5.02808*ALOG10(tk)
    !      vpres(i,j,k) = 10.**(c1-1.3816E-7*10.**p1+8.1328E-3*10.**p2-2949.076/tk)
    !      !vpres(i,j,k) = 6.112*exp(17.67*td(i,j,k)/(td(i,j,k)+243.5))
    !    END DO
    !  END DO
    !END DO

    RETURN
  END SUBROUTINE arpsgrid_getvpress

  SUBROUTINE arpsgrid_getrh( t, rh, tem1, istatus)
    IMPLICIT NONE
    REAL,    INTENT(IN)  :: t   (arpsgrid%nx,arpsgrid%ny,arpsgrid%nz)
    REAL,    INTENT(OUT) :: rh  (arpsgrid%nx,arpsgrid%ny,arpsgrid%nz)
    REAL,    INTENT(OUT) :: tem1(arpsgrid%nx,arpsgrid%ny,arpsgrid%nz)
    INTEGER, INTENT(OUT) :: istatus

  !---------------------------------------------------------------------

    INTEGER :: i,j,k
    INTEGER :: nx, ny, nz

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    nx = arpsgrid%nx
    ny = arpsgrid%ny
    nz = arpsgrid%nz

    CALL getqvs(nx,ny,nz,1,nx-1,1,ny-1,1,nz-1,arpsgrid%p,t,tem1)

    DO k = 1, nz
      DO j = 1, ny
        DO i = 1, nx
          rh (i,j,k) = MIN( MAX( arpsgrid%qv(i,j,k)/tem1(i,j,k), 0.0), 1.0)
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE arpsgrid_getrh

!############## COAMPS related diagnostic subroutines ##################

  SUBROUTINE compute_arpsgrid_sigma( istatus )
  !
  !-----------------------------------------------------------------------
  !
  !  PURPOSE:
  !
  !  Get corresponding COAMPS sigma levels on ARPS grid
  !
  !-----------------------------------------------------------------------
  !
  !  Author: Yunheng Wang
  !  Date: 05/10/11
  !
  !  MODIFICATION HISTORY:
  !
  !
  !-----------------------------------------------------------------------
  !
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: istatus
  !
  !-----------------------------------------------------------------------
  !
  !  Variable Declarations.
  !
  !-----------------------------------------------------------------------
  !
    INCLUDE 'mp.inc'

    INTEGER :: i,j,k, ktop, kbot,kka,kk
    INTEGER :: ibgn,iend,jbgn,jend
    INTEGER :: icount

    REAL, ALLOCATABLE :: sigmma(:), dsigma(:), sigmwa(:)
    REAL :: hmean, hmin, hmax, hdep
    REAL, ALLOCATABLE :: pmax(:), pmin(:), pmean(:), pgrd(:)
  !
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !
    istatus = 0

    ibgn = 2
    iend = arpsgrid%nx-2
    jbgn = 2
    jend = arpsgrid%ny-2

    IF (loc_x == 1) ibgn = 1
    IF (loc_y == 1) jbgn = 1
    IF (loc_x == nproc_x) iend = arpsgrid%nx-1
    IF (loc_y == nproc_y) jend = arpsgrid%ny-1

  !---------------------------------------------------------------------

    hmean = 0.0
    hmin  = 1E+7
    hmax  = 1E-7

    ALLOCATE(pgrd(arpsgrid%nz),  STAT = istatus)
    ALLOCATE(pmean(arpsgrid%nz), STAT = istatus)
    ALLOCATE(pmin(arpsgrid%nz),  STAT = istatus)
    ALLOCATE(pmax(arpsgrid%nz),  STAT = istatus)

    pmean(:) = 0.0
    pmax(:)  = 1E-7
    pmin(:)  = 1E+7

    ktop = arpsgrid%nz-1
    kbot = 2
    icount = 0
    DO j = jbgn,jend
      DO i = ibgn,iend
        hdep = arpsgrid%zp(i,j,ktop)-arpsgrid%zp(i,j,kbot)
        IF (hdep > hmax) hmax = hdep
        IF (hdep < hmin) hmin = hdep
        hmean = hmean+hdep
        icount = icount+1

        DO k = 1,arpsgrid%nz
          pgrd(k) = arpsgrid%p(i,j,k)
          IF (pgrd(k) > pmax(k)) pmax(k) = pgrd(k)
          IF (pgrd(k) < pmin(k)) pmin(k) = pgrd(k)
          pmean(k) = pmean(k) + pgrd(k)
        END DO
      END DO
    END DO

    hmean = hmean/icount
    DO k = 1,arpsgrid%nz
      pmean(k) = pmean(k)/icount
    END DO

    CALL mpmaxar(pmax,arpsgrid%nz)
    CALL mpminar(pmin,arpsgrid%nz)
    CALL mpsumr(pmean,arpsgrid%nz)
    DO k = 1,arpsgrid%nz
      pmean(k) = pmean(k)/(100.*nprocs)
    END DO
    pmax(:) = pmax(:)/100.
    pmin(:) = pmin(:)/100.

    IF (myproc == 0) THEN
      WRITE(*,'(1x,a)') 'ARPS mean/min/max pressure (mb) at each vertical levels:'
      WRITE(*,'(1x,4a)') '    k','      mean','       min','       max'
      WRITE(*,'(1x,4a)') '-----',' ---------',' ---------',' ---------'
      DO k = arpsgrid%nz,1,-1
        WRITE(*,'(1x,I5,3F10.2)') k,pmean(k),pmin(k),pmax(k)
      END DO
      WRITE(*,*)
    END IF

    CALL mpmaxr(hmax)
    CALL mpminr(hmin)
    CALL mptotal(hmean)
    hmean = hmean/nprocs

    IF (myproc == 0) THEN
      WRITE(*,'(1x,a)') 'ARPS mean/min/max atmospheric depth (meter) are:'
      WRITE(*,'(1x,4a)') '     ','      mean','       min','       max'
      WRITE(*,'(1x,4a)') '     ',' ---------',' ---------',' ---------'
      WRITE(*,'(6x,3F10.3)') hmean,hmin,hmax
      WRITE(*,*)
    END IF

    kka = ktop-kbot+1 - 1

    ALLOCATE(sigmma(kka+1), STAT = istatus)
    ALLOCATE(dsigma(kka+1), STAT = istatus)
    ALLOCATE(sigmwa(kka+1), STAT = istatus)

    DO k = 1, kka+1
      kk = kka-k+3
      sigmwa(k) = hmean*(arpsgrid%zp (2,2,kk)-arpsgrid%zp(2,2,2))/(hmean-arpsgrid%zp(2,2,2))
      sigmma(k) = hmean*(arpsgrid%zps(2,2,kk-1)-arpsgrid%zp(2,2,2))/(hmean-arpsgrid%zp(2,2,2))
    END DO
    DO k = 1, kka
      dsigma(k) = sigmwa(k)-sigmwa(k+1)
    END DO

    IF (myproc == 0) THEN
      WRITE(*,'(1x,a)') 'ARPS sigma levels (COAMPS definition, meter) at each vertical levels:'
      WRITE(*,'(1x,5a)') 'k(arps)','    k','    sigmma','    sigmwa','    dsigma'
      WRITE(*,'(1x,5a)') '-------',' ----',' ---------',' ---------',' ---------'
      DO k = 1,kka
        kk = kka-k+3
        WRITE(*,'(1x,I7,I5,3F10.2)') kk,k,sigmma(k),sigmwa(k),dsigma(k)
      END DO
      k = kka+1; kk = kka-k+3
      WRITE(*,'(1x,I7,2a,F10.3,a)') kk,' kka+1','         ',sigmwa(k),'          '
    END IF

  !---------------------------------------------------------------------

    DEALLOCATE(sigmma, dsigma, sigmwa)

    RETURN
  END SUBROUTINE compute_arpsgrid_sigma

!############## WRF related diagnostic subroutines #####################

  SUBROUTINE compute_arpsgrid_eta( iamroot, nprocs,istatus )
  !#####################################################################
  !
  ! PURPOSE:
  !   This subroutine computes the eta values that specify
  !   WRF ARW vertical coordinates. The purpose is to provide guide
  !   for a close WRF vertcial layers as the orginal ARPS vertical layers.
  !
  !   It is based on the corresponding algorithm of REAL.EXE in WRFV3.3.1.
  !
  !#####################################################################

    IMPLICIT NONE

    LOGICAL, INTENT(IN)  :: iamroot
    INTEGER, INTENT(IN)  :: nprocs
    INTEGER, INTENT(OUT) :: istatus

  !---------------------------------------------------------------------
    INTEGER :: nx, ny, nz
    INTEGER :: i,j,k
    INTEGER :: ia,ja

    REAL    :: ptop
    REAL, ALLOCATABLE :: wrfeta(:), zovrsfc(:)

    REAL, ALLOCATABLE :: arps_t(:,:,:), arps_pd(:,:,:), arps_eta(:,:,:)

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    nx = arpsgrid%nx
    ny = arpsgrid%ny
    nz = arpsgrid%nz

    ia = 2
    ja = 2
  !---------------------------------------------------------------------

    ALLOCATE(arps_t  (nx,ny,nz), STAT = istatus)
    ALLOCATE(arps_pd (nx,ny,nz), STAT = istatus)
    ALLOCATE(arps_eta(nx,ny,nz), STAT = istatus)
    ALLOCATE(wrfeta  (nz),       STAT = istatus)
    ALLOCATE(zovrsfc (nz),       STAT = istatus)

    DO k = 1, nz
      DO j = 1, ny
        DO i = 1, nx
          arps_t(i,j,k) = arpsgrid%pt(i,j,k)
        END DO
      END DO
    END DO

    CALL theta_to_t(arps_t,1,nx,1,ny,1,nz,arpsgrid%p,istatus)

    CALL integ_moist( nx,ny,nz,arpsgrid%qv,arpsgrid%p,arps_t,arpsgrid%zp,  &
                      arps_pd,istatus)

    ptop = MINVAL(arps_pd(:,:,nz-2))
    CALL mpminr(ptop)

    wrfeta(:) = 0.0
    DO k = 2, nz-1
      DO j = 2, ny-1
        DO i = 2, nx-1
          arps_eta(i,j,k) = (arps_pd(i,j,k)-ptop)/(arps_pd(i,j,2)-ptop)
          !write(*,*) i,j,k,arps_eta(i,j,k)
          wrfeta(k) = wrfeta(k) + arps_eta(i,j,k)
        END DO
      END DO
      wrfeta(k) = wrfeta(k) / ((nx-2)*(ny-2))
    END DO

    CALL mpsumr(wrfeta,nz)
    wrfeta(:) = wrfeta(:)/nprocs

    DO k = 1,nz
      zovrsfc(k) = arpsgrid%zp(ia,ja,k)-arpsgrid%zp(ia,ja,2)
    END DO

    IF (iamroot) THEN
      WRITE(*,'(1x,a,/,1x,a,F7.2,a,/,1x,a)') '===================',     &
      'NOTE: The ARPS model top has pressure P_TOP is = ', ptop/100.0,' hPa.',   &
      '==================='
      WRITE(*,'(1x,a)') 'Average WRF eta values for each ARPS vertical layers'
      WRITE(*,'(1x,a,2(I0,a))') '  No.   ETA (avg)     ARPSZ(m)  dz(m) at (',ia,',',ja,')'
      WRITE(*,'(1x,a)')         '-----   ---------    ---------  ---------'
      WRITE(*,'(1x,I5,3F12.3)') (k,wrfeta(k),zovrsfc(k),                &
                                 zovrsfc(k+1)-zovrsfc(k), k=2,nz-1)
      WRITE(*,'(1x,a,200(F6.3,a))') 'eta_levels = ',(wrfeta(k),',', k=2,nz-1)
    END IF

    DEALLOCATE (wrfeta, zovrsfc)
    DEALLOCATE (arps_t, arps_pd, arps_eta)

    RETURN
  END SUBROUTINE compute_arpsgrid_eta

!#################### Private subroutines ##############################

  SUBROUTINE integ_moist ( nx,ny,nz,q_in, p_in, t_in, ght_in,           &
                           pd_out, istatus )

  !#####################################################################
  !
  !  Integrate the moisture field vertically.  Mostly used to get the total
  !  vapor pressure, which can be subtracted from the total pressure to get
  !  the dry pressure.
  !
  !#####################################################################

     IMPLICIT NONE

     INTEGER, INTENT(IN) :: nx,ny,nz

     REAL, DIMENSION(nx,ny,nz), INTENT(IN)  :: q_in, p_in, t_in ! mass grid
     REAL, DIMENSION(nx,ny,nz), INTENT(IN)  :: ght_in           ! W grid
     REAL, DIMENSION(nx,ny,nz), INTENT(OUT) :: pd_out           ! W grid

     INTEGER, INTENT(OUT) :: istatus

  !-----------------------------------------------------------------------
     !  Local vars

     INTEGER :: i , j , k

     REAL :: rhobar , qbar , dz
     REAL :: intq

     REAL , PARAMETER :: Rd = 287.0  ! Gas constant for dry air  (m**2/(s**2*K))

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

     istatus = 0

     DO j = 1,ny
        DO i = 1,nx

           intq = 0. !  Initialize the integrated quantity of moisture to zero.

           !  Account for the moisture above the ground.

           pd_out(i,j,nz)   = p_in(i,j,nz)
           pd_out(i,j,nz-1) = EXP(0.5*(LOG(p_in(i,j,nz-1))+LOG(p_in(i,j,nz-2))))   ! physical model top

           DO k = nz-2,2,-1
              rhobar =  p_in(i,j,k) / ( Rd * t_in(i,j,k) )  ! W grid
              qbar   =  q_in(i,j,k)                         ! W grid
              dz     =  ght_in(i,j,k+1) - ght_in(i,j,k)
              intq   =  intq + g * qbar * rhobar  * dz
              pd_out(i,j,k) = EXP(0.5*(LOG(p_in(i,j,k))+LOG(p_in(i,j,k-1)))) - intq
              !write(*,*) i,j,k,intq,p_in(i,j,k),pd_out(i,j,k)
           END DO

           pd_out(i,j,1) =  EXP(2.*LOG(p_in(i,j,1))-LOG(p_in(i,j,2))) - intq  ! fake level underground

        END DO
     END DO

  END SUBROUTINE integ_moist

  !#######################################################################

  SUBROUTINE theta_to_t(t,ibgn,iend,jbgn,jend,kbgn,kend,p,istatus)

  !#######################################################################

    IMPLICIT NONE

    INTEGER, INTENT(IN)     :: ibgn, iend, jbgn,jend,kbgn,kend

    REAL,    INTENT(INOUT)  :: t(ibgn:iend,jbgn:jend,kbgn:kend)
    REAL,    INTENT(IN)     :: p(ibgn:iend,jbgn:jend,kbgn:kend)
    INTEGER, INTENT(OUT)    :: istatus

  !-----------------------------------------------------------------------

    INTEGER :: i,j,k

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    DO k = kbgn,kend
      DO j = jbgn, jend
        DO i = ibgn, iend
          t(i,j,k) = t(i,j,k)/( (p0/p(i,j,k))**rddcp )
        END DO
      END DO
    END DO

    RETURN
  END SUBROUTINE theta_to_t

END MODULE module_input_arpsgrid
