  module Vmodule_3dvar
!
!
! 1.1 background state vector
!     -----------------------
!
  IMPLICIT NONE
  SAVE
!
  INTEGER :: mgra, nxm, nwork
  INTEGER :: numctr
!
  REAL, DIMENSION(:),      ALLOCATABLE :: ctrv,grad,xgus
  REAL, DIMENSION(:),      ALLOCATABLE :: swork,ywork,diag,work
  REAL, DIMENSION(:,:,:),  ALLOCATABLE :: gdscal
  REAL, DIMENSION(:,:,:),  ALLOCATABLE :: gdu_err,gdv_err,               &
                                          gdp_err,gdt_err,gdq_err,gdw_err
  REAL, DIMENSION (:,:,:), ALLOCATABLE :: u_ctr,v_ctr,p_ctr,t_ctr,q_ctr,w_ctr
  REAL, DIMENSION (:,:,:), ALLOCATABLE :: psi,phi
  REAL, DIMENSION (:,:,:), ALLOCATABLE :: rhostru, rhostrv,rhostrw,div3
!

  CONTAINS

    SUBROUTINE allocateArray4Var(nx,ny,nz)
      INTEGER, INTENT(IN) :: nx,ny,nz
      INTEGER :: ierr
      INTEGER :: ctrsize

      INTEGER :: ibgn, iend, jbgn, jend
      INTEGER :: iuend,jvend
      INCLUDE 'mp.inc'

      mgra    = 2
      ctrsize = 6*nx*ny*nz
      nxm     = ctrsize*mgra
      nwork   = ctrsize+2*mgra

      ibgn = 1
      iend = nx-1
      jbgn = 1
      jend = ny-1

      !iuend = nx
      !jvend = ny

      IF (loc_x > 1) ibgn = 2
      IF (loc_y > 1) jbgn = 2

      IF (loc_x < nproc_x) THEN
        iend  = nx-2
        !iuend = iend
      END IF

      IF (loc_y < nproc_y) THEN
        jend  = ny-2
        !jvend = jend
      END IF

      numctr = 3*(nz-1)*(iend-ibgn+1) *(jend-jbgn+1)     & ! 3 scalar arrays
             +   (nz-1)*(iend-ibgn+1) *(jend-jbgn+1)     & ! U
             +   (nz-1)*(iend-ibgn+1) *(jend-jbgn+1)     & ! V
             +       nz*(iend-ibgn+1) *(jend-jbgn+1)       ! W

      ALLOCATE ( swork(nxm)   )
      ALLOCATE ( ywork(nxm)   )
      ALLOCATE ( diag(ctrsize)  )
      ALLOCATE ( work(nwork) )

      IF ( ALLOCATED(xgus) ) DEALLOCATE (xgus)
      ALLOCATE ( xgus(ctrsize), stat=ierr )
      CALL check_alloc_status(ierr, "3dvar:xgus")

      IF ( ALLOCATED(ctrv) ) DEALLOCATE (ctrv)
      ALLOCATE ( ctrv(ctrsize), stat=ierr )
      CALL check_alloc_status(ierr, "3dvar:ctrv")

      IF ( ALLOCATED(grad) ) DEALLOCATE (grad)
      ALLOCATE ( grad(ctrsize), stat=ierr)
      CALL check_alloc_status(ierr, "3dvar:grad")

      IF ( ALLOCATED(gdu_err) ) DEALLOCATE (gdu_err)
      ALLOCATE ( gdu_err(nx,ny,nz), stat=ierr)
      CALL check_alloc_status(ierr, "3dvar:gdu_err")
      gdu_err(:,:,:) = 0.0

      IF ( ALLOCATED(gdv_err) ) DEALLOCATE (gdv_err)
      ALLOCATE ( gdv_err(nx,ny,nz), stat=ierr )
      CALL check_alloc_status(ierr, "3dvar:gdv_err")
      gdv_err(:,:,:) = 0.0

      IF ( ALLOCATED(gdp_err) ) DEALLOCATE (gdp_err)
      ALLOCATE ( gdp_err(nx,ny,nz), stat=ierr )
      CALL check_alloc_status(ierr, "3dvar:gdp_err")
      gdp_err(:,:,:) = 0.0

      IF ( ALLOCATED(gdt_err) ) DEALLOCATE (gdt_err)
      ALLOCATE ( gdt_err(nx,ny,nz), stat=ierr )
      CALL check_alloc_status(ierr, "3dvar:gdt_err")
      gdt_err(:,:,:) = 0.0

      IF ( ALLOCATED(gdq_err) ) DEALLOCATE (gdq_err)
      ALLOCATE ( gdq_err(nx,ny,nz) , stat=ierr)
      CALL check_alloc_status(ierr, "3dvar:gdq_err")
      gdq_err(:,:,:) = 0.0

      IF ( ALLOCATED(gdw_err) ) DEALLOCATE (gdw_err)
      ALLOCATE ( gdw_err(nx,ny,nz) , stat=ierr)
      CALL check_alloc_status(ierr, "3dvar:gdw_err")
      gdw_err(:,:,:) = 0.0

      IF ( ALLOCATED(  u_ctr) ) DEALLOCATE (  u_ctr)
      ALLOCATE (   u_ctr(nx,ny,nz) , stat=ierr)
      CALL check_alloc_status(ierr, "3dvar:u_ctr")
      u_ctr(:,:,:) = 0.0

      IF ( ALLOCATED(  v_ctr) ) DEALLOCATE (  v_ctr)
      ALLOCATE (   v_ctr(nx,ny,nz) , stat=ierr)
      CALL check_alloc_status(ierr, "3dvar:v_ctr")
      v_ctr(:,:,:) = 0.0

      IF ( ALLOCATED(    psi) ) DEALLOCATE (    psi)
      ALLOCATE (     psi(nx,ny,nz) , stat=ierr)
      CALL check_alloc_status(ierr, "3dvar:psi")
      psi(:,:,:) = 0.0

      IF ( ALLOCATED(    phi) ) DEALLOCATE (    phi)
      ALLOCATE (     phi(nx,ny,nz) , stat=ierr)
      CALL check_alloc_status(ierr, "3dvar:phi")
      phi(:,:,:) = 0.0

      IF ( ALLOCATED(  p_ctr) ) DEALLOCATE (  p_ctr)
      ALLOCATE (   p_ctr(nx,ny,nz) , stat=ierr)
      CALL check_alloc_status(ierr, "3dvar:p_ctr")
      p_ctr(:,:,:) = 0.0

      IF ( ALLOCATED(  t_ctr) ) DEALLOCATE (  t_ctr)
      ALLOCATE (   t_ctr(nx,ny,nz) , stat=ierr)
      CALL check_alloc_status(ierr, "3dvar:t_ctr")
      t_ctr(:,:,:) = 0.0

      IF ( ALLOCATED(  q_ctr) ) DEALLOCATE (  q_ctr)
      ALLOCATE (   q_ctr(nx,ny,nz) , stat=ierr)
      CALL check_alloc_status(ierr, "3dvar:q_ctr")
      q_ctr(:,:,:) = 0.0

      IF ( ALLOCATED(  w_ctr) ) DEALLOCATE (  w_ctr)
      ALLOCATE (   w_ctr(nx,ny,nz) , stat=ierr)
      CALL check_alloc_status(ierr, "3dvar:w_ctr")
      w_ctr(:,:,:) = 0.0
    !
    !-----------------------------------------------------------------------
    !
      IF ( ALLOCATED(gdscal) ) DEALLOCATE (gdscal)
      ALLOCATE ( gdscal(nx,ny,nz) , stat=ierr)
      CALL check_alloc_status(ierr, "3dvar:gdscal")
      gdscal(:,:,:) = 0.0

      IF ( ALLOCATED(rhostru) ) DEALLOCATE (rhostru)
      ALLOCATE ( rhostru(nx,ny,nz) , stat=ierr)
      CALL check_alloc_status(ierr, "3dvar:rhostru")
      rhostru(:,:,:) = 0.0

      IF ( ALLOCATED(rhostrv) ) DEALLOCATE (rhostrv)
      ALLOCATE ( rhostrv(nx,ny,nz) , stat=ierr)
      CALL check_alloc_status(ierr, "3dvar:rhostrv")
      rhostrv(:,:,:) = 0.0

      IF ( ALLOCATED(rhostrw) ) DEALLOCATE (rhostrw)
      ALLOCATE ( rhostrw(nx,ny,nz) , stat=ierr)
      CALL check_alloc_status(ierr, "3dvar:rhostrw")
      rhostrw(:,:,:) = 0.0

      IF ( ALLOCATED(div3   ) ) DEALLOCATE (div3   )
      ALLOCATE ( div3   (nx,ny,nz) , stat=ierr)
      CALL check_alloc_status(ierr, "3dvar:div3")
      div3 (:,:,:) = 0.0

    END  SUBROUTINE allocatearray4Var

!#######################################################################

    SUBROUTINE deallocateArray4Var

      DEALLOCATE ( swork  )
      DEALLOCATE ( ywork  )
      DEALLOCATE (  diag  )
      DEALLOCATE (  work  )

      DEALLOCATE (u_ctr)
      DEALLOCATE (v_ctr)
      DEALLOCATE (p_ctr)
      DEALLOCATE (t_ctr)
      DEALLOCATE (q_ctr)
      DEALLOCATE (w_ctr)

      DEALLOCATE ( xgus )
      DEALLOCATE ( ctrv )
      DEALLOCATE ( grad )

      DEALLOCATE ( gdu_err )
      DEALLOCATE ( gdv_err )
      DEALLOCATE ( gdp_err )
      DEALLOCATE ( gdt_err )
      DEALLOCATE ( gdq_err )
      DEALLOCATE ( gdw_err )
      DEALLOCATE ( gdscal )

      DEALLOCATE ( rhostru )
      DEALLOCATE ( rhostrv )
      DEALLOCATE ( rhostrw )
      DEALLOCATE ( div3    )

    END SUBROUTINE deallocateArray4Var

END module Vmodule_3dvar

