!
!##################################################################
!##################################################################
!######                                                      ######
!######            SUBROUTINE MULTIMOMENT_DRIVE              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE multimoment_driver(mscheme,nx,ny,nz,dtbig1,zp,w,ptprt,ptbar, &
              pprt,pbar,ppi,qv,qvbar,qscalar,raing,prcrate,             &
              tk,p,q,tkm,qm,pm,ww,zp2d,lr,sr,qnz,qnzm,tem)
              !,mpteqnterms,N0x)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  The driver for the multimoment microphysics scheme developed by
!  Jason Milbrandt (McGill University / Meteorological Services of Canada).
!
!  REFERENCE:
!
!  Milbrandt, J. A. and M. K. Yau, 2005: A multi-moment bulk microphysics
!  parameterization. Part I: Analysis of the role of the spectral shape
!  parameter. J. Atmos. Sci., 62, 3051-3064.
!
!  Milbrandt, J. A. and M. K. Yau, 2005: A multi-moment bulk microphysics
!  parameterization. Part II: A proposed three-moment closure and scheme
!  description. J. Atmos. Sci., 62, 3065-3081.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  02/21/2006.
!
!  MODIFICATION HISTORY:
!     Youngsun Jung   11/18/2009
!       Added the grpl_ON switch to suppress graupel production.
!
!     Bryan Putnam    12/7/2010
!       Added the hail_ON switch to suppress hail production.
!
!     Dan Dawson      07/19/2011
!       Rewrote/organized interface for latest version of full MY (with switches)
!       scheme (version 2.20):
!       Renamed driver file to micro_MY.f90.  Now calls the latest version
!       of the MY scheme contained in the files my3mom_main_mod.f90,
!       my3mom_fncs_mod.f90, and my3mom_sedi_mod.f90, replacing the
!       original interface.  Incorporated all changes from previous
!       interface to new one as appropriate.  The new interface should
!       facilitate easier updating of the MY scheme to newer versions as
!       needed.
!
!-----------------------------------------------------------------------

  USE my_tmom_mod

  IMPLICIT NONE

  INCLUDE 'phycst.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'timelvls.inc'

  INTEGER, INTENT(IN)    :: mscheme              ! multimoment schemes, 1-4
  INTEGER, INTENT(IN)    :: nx, ny, nz
  REAL,    INTENT(IN)    :: dtbig1
  REAL,    INTENT(IN)    :: zp(nx,ny,nz)
  REAL,    INTENT(IN)    :: w (nx,ny,nz)

  REAL,    INTENT(INOUT) :: ptprt(nx,ny,nz,nt)
  REAL,    INTENT(INOUT) :: pprt (nx,ny,nz,nt)
  REAL,    INTENT(INOUT) :: qv   (nx,ny,nz,nt)
  REAL,    INTENT(IN)    :: ptbar(nx,ny,nz)
  REAL,    INTENT(IN)    :: pbar (nx,ny,nz)
  REAL,    INTENT(IN)    :: qvbar(nx,ny,nz)

  REAL,    INTENT(IN)    :: ppi  (nx,ny,nz)

  REAL,    INTENT(INOUT) :: qscalar(nx,ny,nz,nt,nscalar)
  REAL,    INTENT(INOUT) :: raing  (nx,ny)
  REAL,    INTENT(INOUT) :: prcrate(nx,ny)

  !
  ! local arrays, Memory order in i-k-j
  !
  REAL     :: tk (nx,nz,ny)   ! temperature at time tfuture
  REAL     :: tkm(nx,nz,ny)   ! temperature at time tpast

  REAL     :: p (nx,nz,ny)    ! total pressure at time tfuture
  REAL     :: pm(nx,nz,ny)    ! total pressure at time tpast

  REAL     :: q (nx,nz,ny)    ! total mixing ratio at time tfuture
  REAL     :: qm(nx,nz,ny)    ! total mixing ratio at time tpast

  REAL     :: ww  (nx,nz)
  REAL     :: zp2d(nx,nz)

  REAL     :: lr(nx,ny)       ! liquid precipitation rate
  REAL     :: sr(nx,ny)       ! soild precipitation rate

  REAL     :: qnz (nx,nz,nscalar) ! hydrometeor arrays at time tfuture
  REAL     :: qnzm(nx,nz,nscalar) ! at time tpresent, qc, qr, qi, qs, qg, qh
                                  !                   nc, nr, ni, ns, ng, nh
                                  !                       zr, zi, zs, zg, zh
  REAL     :: tem(nx,nz,12)       ! used only when scheme == 1

! Added by DTD, 01/07: rain evaporation and intercept parameters added
! Cleaned up 01/08.

  !REAL     :: mpteqnterms(nx,ny,nz,28)
  !REAL     :: mpteqnterms_ikj(nx,nz,28,ny)
  !REAL(KIND=8)     :: N0x(nx,ny,nz,6)
  !REAL(KIND=8)     :: N0x_ikj(nx,nz,6,ny)

!
!-----------------------------------------------------------------------
!
! Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER :: i,j,k,nq
  REAL    :: deltat

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  !mpteqnterms = 0.0
  !mpteqnterms_ikj = 0.0
  !N0x = 0.0
  !N0x_ikj = 0.0

!-----------------------------------------------------------------------
!
! Compute total temperature and total pressure
!
!-----------------------------------------------------------------------

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1

        p  (i,k,j) =  pprt (i,j,k,tfuture) + pbar(i,j,k)
        tk (i,k,j) = (ptprt(i,j,k,tfuture) + ptbar(i,j,k)) * ppi(i,j,k)
        q  (i,k,j) =  qv   (i,j,k,tfuture)

        pm (i,k,j) =  pprt (i,j,k,tpresent)+ pbar(i,j,k)
        tkm(i,k,j) = (ptprt(i,j,k,tpresent)+ ptbar(i,j,k)) * ppi(i,j,k)
        qm (i,k,j) =  qv   (i,j,k,tpresent)

      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
! Call multimoment microphysics scheme
!
!-----------------------------------------------------------------------

  !IF( sadvopt /= 4 .and. tintegopt == 1) THEN         ! Leapfrog scheme
  !  deltat = 2.*dtbig1
  !ELSE                                                ! Forward scheme
    deltat = dtbig1
  !END IF

!$OMP PARALLEL DO   &
!$OMP PRIVATE ( j )

  DO j = 1,ny-1

    DO k = 1,nz-1
      DO i = 1,nx-1
        ww  (i,k) = w(i,j,k)
        zp2d(i,k) = zp(i,j,k)
        qnz (i,k,:) = qscalar(i,j,k,tfuture, :)
        qnzm(i,k,:) = qscalar(i,j,k,tpresent,:)
      END DO
    END DO

!   CALL multimoment(mscheme,nx,nz,deltat,ww,zp2d,           &
!            p (:,:,j),tk (:,:,j),q (:,:,j),qnz,             &
!            pm(:,:,j),tkm(:,:,j),qm(:,:,j),qnzm,            &
!            lr(:,j),  sr(:,j), tem, j, mpteqnterms_ikj(:,:,:,j),N0x_ikj(:,:,:,j), &
!            ntcloud,n0rain,n0snow,n0grpl,n0hail,rhoice,rhosnow,rhogrpl,rhohail,   &
!            alpharain,alphaice,alphasnow,alphagrpl,alphahail)

    CALL MYTMOM_MAIN(ww,tk(:,:,j),q(:,:,j),qnz,p(:,:,j),                &
                     tkm(:,:,j),qm(:,:,j),qnzm,pm(:,:,j),lr(:,j),       &
                     sr(:,j),tem,zp2d,deltat,nx,nz,j,mscheme,ntcloud,n0rain,n0snow,  &
                     n0grpl,n0hail,rhoice,rhosnow,rhogrpl,rhohail,alpharain,         &
                     alphaice,alphasnow,alphagrpl,alphahail)

    DO nq = 1,nscalar
      DO k = 1,nz-1
        DO i = 1,nx-1
          qscalar(i,j,k,tfuture, nq) = qnz (i,k,nq)
          qscalar(i,j,k,tpresent,nq) = qnzm(i,k,nq)
        END DO
      END DO
    END DO

  END DO

!$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
! Convert back to potential temperature
!
!-----------------------------------------------------------------------

  DO k = 1,nz-1
    DO j = 1,ny-1
      DO i = 1,nx-1
        ptprt(i,j,k,tfuture) = tk(i,k,j)/ppi(i,j,k) - ptbar(i,j,k)
        qv   (i,j,k,tfuture) = q (i,k,j)

        ptprt(i,j,k,tpresent) = tkm(i,k,j)/ppi(i,j,k) - ptbar(i,j,k)
        qv   (i,j,k,tpresent) = qm (i,k,j)

        ! DTD: convert ikj memory-order to regular ijk for mp temp eqn terms
        ! and intercept parameter arrays

        !DO nq=1,28
        !  mpteqnterms(i,j,k,nq) = mpteqnterms_ikj(i,k,nq,j)
        !END DO
        !
        !DO nq=1,6
        !  N0x(i,j,k,nq) = N0x_ikj(i,k,nq,j)
        !END DO

      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
!
! Compute rainfall rate and grid scale precipitation
!
!-----------------------------------------------------------------------

  DO j = 1,ny-1
    DO i = 1,nx-1
      prcrate(i,j) = lr(i,j) + sr(i,j)
      raing  (i,j) = raing(i,j) + prcrate(i,j)*dtbig
    END DO
  END DO

  RETURN
END SUBROUTINE multimoment_driver
