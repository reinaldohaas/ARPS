  module module_soilArray
!
!
! 1.1 background state vector
!     -----------------------
!
  IMPLICIT NONE
  SAVE
!
!
!-----------------------------------------------------------------------
!
!  ARPS Surface variables:
!
!-----------------------------------------------------------------------
!
  INTEGER, ALLOCATABLE :: soiltyp (:,:,:) ! Soil type
  REAL,    ALLOCATABLE ::    stypfrct(:,:,:) ! Soil type fraction
  INTEGER, ALLOCATABLE :: vegtyp (:,:)    ! Vegetation type
  REAL,    ALLOCATABLE ::    lai    (:,:)    ! Leaf Area Index
  REAL,    ALLOCATABLE ::    roufns (:,:)    ! Surface roughness
  REAL,    ALLOCATABLE ::    veg    (:,:)    ! Vegetation fraction

  REAL,    ALLOCATABLE :: tsoil(:,:,:,:) ! Soil temperature (K)
  REAL,    ALLOCATABLE :: qsoil(:,:,:,:) ! Soil moisture
  REAL,    ALLOCATABLE :: qvsfc(:,:,:)   ! Effective qv at sfc.

  REAL,    ALLOCATABLE :: wetcanp(:,:,:) ! Canopy water amount
  REAL,    ALLOCATABLE :: snowdpth(:,:)  ! Snow depth (:)

  REAL,    ALLOCATABLE :: ptcumsrc(:,:,:)! Source term in pt-equation due
                                         ! to cumulus parameterization
  REAL,    ALLOCATABLE :: qcumsrc(:,:,:,:) ! Source term in Qv equation due
                                         ! to cumulus parameterization

  REAL,    ALLOCATABLE :: w0avg(:,:,:)   ! a closing running average vertical
                                         ! velocity in 10min for K-F scheme
  REAL,    ALLOCATABLE :: kfraincv(:,:)  ! K-F convective rainfall (:)
  INTEGER, ALLOCATABLE :: nca(:,:)    ! K-F counter for CAPE release
  REAL,    ALLOCATABLE :: cldefi(:,:)    ! BMJ cloud efficiency
  REAL,    ALLOCATABLE :: xland(:,:)     ! BMJ land mask
                                      !   (1.0 = land, 2.0 = sea)
  REAL,    ALLOCATABLE :: bmjraincv(:,:) ! BMJ convective rainfall (cm)


  REAL,    ALLOCATABLE :: raing(:,:)     ! Grid supersaturation rain
  REAL,    ALLOCATABLE :: rainc(:,:)     ! Cumulus convective rain
  REAL,    ALLOCATABLE :: prcrate(:,:,:) ! precipitation rate (kg/(m**2*s))
                                         ! prcrate(:,:,:) = total precipitation rate
                                         ! prcrate(:,:,:) = grid scale precip. rate
                                         ! prcrate(:,:,:) = cumulus precip. rate
                                         ! prcrate(:,:,:) = microphysics precip. rate

  REAL,    ALLOCATABLE :: usflx (:,:)    ! Surface flux of u-momentum (kg/(m*s**2))
  REAL,    ALLOCATABLE :: vsflx (:,:)    ! Surface flux of v-momentum (kg/(m*s**2))
  REAL,    ALLOCATABLE :: ptsflx(:,:)    ! Surface heat flux (K*kg/(m*s**2))
  REAL,    ALLOCATABLE :: qvsflx(:,:)    ! Surface moisture flux (kg/(m**2*s))


  CONTAINS

  SUBROUTINE allocatesoilArray(nx,ny,nz,nstyps,nzsoil)

    INTEGER:: nx,ny,nz,nstyps,nzsoil,istatus
!
!----------------------------------------------------------------------
!
!  Allocate adas arrays
!
!-----------------------------------------------------------------------
!
    ALLOCATE(soiltyp(nx,ny,nstyps),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:soiltyp")
    ALLOCATE(stypfrct(nx,ny,nstyps),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:stypfrct")
    ALLOCATE(vegtyp (nx,ny),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:vegtyp")
    ALLOCATE(lai    (nx,ny),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:lai")
    ALLOCATE(roufns (nx,ny),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:roufns")
    ALLOCATE(veg    (nx,ny),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:veg")

    ALLOCATE(qvsfc  (nx,ny,0:nstyps),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:qvsfc")
    ALLOCATE(tsoil  (nx,ny,nzsoil,0:nstyps),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:tsoil")
    ALLOCATE(qsoil  (nx,ny,nzsoil,0:nstyps),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:qsoil")
    ALLOCATE(wetcanp(nx,ny,0:nstyps),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:wetcanp")
    ALLOCATE(snowdpth(nx,ny),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:snowdpth")

    ALLOCATE(ptcumsrc(nx,ny,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:ptcumsrc")
    ALLOCATE(qcumsrc (nx,ny,nz,5),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:qcumsrc")
    ALLOCATE(w0avg   (nx,ny,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:w0avg")
    ALLOCATE(nca    (nx,ny),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:nca")
    ALLOCATE(kfraincv (nx,ny),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:kfraincv")
    ALLOCATE(cldefi (nx,ny),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:cldefi")
    ALLOCATE(xland  (nx,ny),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:xland")
    ALLOCATE(bmjraincv (nx,ny),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:bmjraincv")
    ALLOCATE(raing  (nx,ny),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:raing")
    ALLOCATE(rainc  (nx,ny),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:rainc")
    ALLOCATE(prcrate(nx,ny,4),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:prcrate")

    ALLOCATE(usflx (nx,ny),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:usflx")
    ALLOCATE(vsflx (nx,ny),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:vsflx")
    ALLOCATE(ptsflx(nx,ny),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:ptsflx")
    ALLOCATE(qvsflx(nx,ny),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:qvsflx")

  END  SUBROUTINE allocatesoilArray

  SUBROUTINE deallocatesoilArray

    DEALLOCATE(soiltyp )
    DEALLOCATE(stypfrct )
    DEALLOCATE(vegtyp  )
    DEALLOCATE(lai     )
    DEALLOCATE(roufns  )
    DEALLOCATE(veg     )

    DEALLOCATE(qvsfc   )
    DEALLOCATE(tsoil   )
    DEALLOCATE(qsoil   )
    DEALLOCATE(wetcanp )
    DEALLOCATE(snowdpth )

    DEALLOCATE(ptcumsrc )
    DEALLOCATE(qcumsrc  )
    DEALLOCATE(w0avg    )
    DEALLOCATE(nca     )
    DEALLOCATE(kfraincv  )
    DEALLOCATE(cldefi  )
    DEALLOCATE(xland  )
    DEALLOCATE(bmjraincv )
    DEALLOCATE(raing  )
    DEALLOCATE(rainc  )
    DEALLOCATE(prcrate )

    DEALLOCATE(usflx )
    DEALLOCATE(vsflx )
    DEALLOCATE(ptsflx )
    DEALLOCATE(qvsflx )

  END SUBROUTINE deallocatesoilArray

END module module_soilArray
