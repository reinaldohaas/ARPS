MODULE module_analysisArrays
!-----------------------------------------------------------------------
!
! MODULE contains variables for
!
!   Surface Station   / Surface observation    (mx_sng)
!   Upper-ari Station / Upper-air observation  (mx_ua)
!   Radar site        / Radar observation      (mx_rad, mx_colrad)
!
! HISTORY:
!
!   11/06/2007 (Yunheng Wang)
!   Initial version.
!
!-----------------------------------------------------------------------
!
!  Surface Station variables
!
!-----------------------------------------------------------------------
!
  INTEGER, ALLOCATABLE :: isrcsng(:)
  INTEGER, ALLOCATABLE :: icatsng(:)
  REAL,    ALLOCATABLE :: latsng(:,:)
  REAL,    ALLOCATABLE :: lonsng(:,:)
  REAL,    ALLOCATABLE :: hgtsng(:,:)
  REAL,    ALLOCATABLE :: xsng(:)
  REAL,    ALLOCATABLE :: ysng(:)
  REAL,    ALLOCATABLE :: trnsng(:)
  INTEGER, ALLOCATABLE :: timesng(:,:)
  CHARACTER(LEN=5), ALLOCATABLE :: stnsng(:,:)
!
!-----------------------------------------------------------------------
!
!  Surface (single-level) read-in observation variables
!
!-----------------------------------------------------------------------
!
  CHARACTER(LEN=8), ALLOCATABLE :: wx(:,:)
  CHARACTER(LEN=8), ALLOCATABLE :: csrcsng(:,:)
  CHARACTER(LEN=1), ALLOCATABLE :: store_emv(:,:,:)
  CHARACTER(LEN=4), ALLOCATABLE :: store_amt(:,:,:)
  INTEGER, ALLOCATABLE :: kloud(:,:),idp3(:,:)
  REAL,    ALLOCATABLE :: store_hgt(:,:,:)
  REAL,    ALLOCATABLE :: obrdsng(:,:,:)

  REAL,    ALLOCATABLE :: obsng  (:,:)
  REAL,    ALLOCATABLE :: odifsng(:,:)
  REAL,    ALLOCATABLE :: oanxsng(:,:)
  REAL,    ALLOCATABLE :: corsng(:,:)
  REAL,    ALLOCATABLE :: thesng(:)
  INTEGER, ALLOCATABLE :: ival(:)

!-----------------------------------------------------------------------
!
!  Upper Air Station variables
!
!-----------------------------------------------------------------------
!
  INTEGER, ALLOCATABLE :: isrcua(:)
  REAL,    ALLOCATABLE :: xua(:)
  REAL,    ALLOCATABLE :: yua(:)
  REAL,    ALLOCATABLE :: trnua(:)
  REAL,    ALLOCATABLE :: elevua(:)
  REAL,    ALLOCATABLE :: hgtua(:,:)
  INTEGER, ALLOCATABLE :: nlevsua(:)
  CHARACTER (LEN=5), ALLOCATABLE :: stnua(:)
!
!-----------------------------------------------------------------------
!
!  Upper-air observation variables
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: obsua (:,:,:)
  REAL, ALLOCATABLE :: odifua(:,:,:)
  REAL, ALLOCATABLE :: oanxua(:,:,:)
  REAL, ALLOCATABLE :: corua(:,:,:)
  REAL, ALLOCATABLE :: theua(:,:)
!
!-----------------------------------------------------------------------
!
!  Radar site variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=5), ALLOCATABLE :: stnrad(:)
  INTEGER, ALLOCATABLE :: isrcrad(:)
  REAL,    ALLOCATABLE :: latrad(:),lonrad(:)
  REAL,    ALLOCATABLE :: elvrad(:)

!
!-----------------------------------------------------------------------
!
!  Radar observation variables
!
!-----------------------------------------------------------------------
  INTEGER, ALLOCATABLE :: irad(:)
  INTEGER, ALLOCATABLE :: nlevrad(:)
  REAL,    ALLOCATABLE :: latradc(:)
  REAL,    ALLOCATABLE :: lonradc(:)
  REAL,    ALLOCATABLE :: xradc(:)
  REAL,    ALLOCATABLE :: yradc(:)
  REAL,    ALLOCATABLE :: trnradc(:)
  REAL,    ALLOCATABLE :: distrad(:)
  REAL,    ALLOCATABLE :: uazmrad(:)
  REAL,    ALLOCATABLE :: vazmrad(:)
  REAL,    ALLOCATABLE :: hgtradc(:,:)
  REAL,    ALLOCATABLE :: theradc(:,:)
  REAL,    ALLOCATABLE :: obsrad(:,:,:)

  REAL,    ALLOCATABLE :: odifrad(:,:,:)
  REAL,    ALLOCATABLE :: oanxrad(:,:,:)
  REAL,    ALLOCATABLE :: corrad(:,:,:)
  REAL,    ALLOCATABLE :: dsdr(:,:)
  REAL,    ALLOCATABLE :: dhdr(:,:)

!
!-----------------------------------------------------------------------
!
!  Quality Control Variables
!
!-----------------------------------------------------------------------
!
  INTEGER, ALLOCATABLE :: qualrdsng(:,:,:)
  INTEGER, ALLOCATABLE :: qualsng(:,:)
  REAL   , ALLOCATABLE :: qobsng(:,:)

  INTEGER, ALLOCATABLE :: qualua(:,:,:)
  REAL   , ALLOCATABLE :: qobsua(:,:,:)

  INTEGER, ALLOCATABLE :: qualrad(:,:,:)
  REAL,    ALLOCATABLE :: qobsrad(:,:,:)

  CONTAINS

  !
  ! Allocate them
  !
  SUBROUTINE allocate_adas_3dvar_arrays(progname,mx_sng,nvar_sng,nvar_anx,ntime, &
                mx_ua, nz_ua, mx_rad,                                            &
                mx_colrad,nvar_rad,nvar_radin,nz_rdr,raduvobs,radrhobs,istatus)

    IMPLICIT NONE
    CHARACTER(LEN=4), INTENT(IN) :: progname
    INTEGER,          INTENT(IN) :: mx_ua, nz_ua, mx_rad
    INTEGER,          INTENT(IN) :: mx_sng, nvar_sng, nvar_anx, ntime
    INTEGER,          INTENT(IN) :: mx_colrad,nz_rdr,nvar_rad,nvar_radin
    INTEGER,          INTENT(IN) :: raduvobs,radrhobs
    INTEGER,          INTENT(OUT) :: istatus
    INTEGER       :: n_mx_rad, n_mx_colrad, n_nz_rdr,n_nvar_radin, n_nvar_rad

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
!  Surface Station variables and Surface (single-level) read-in observation variables
!
!-----------------------------------------------------------------------

    ALLOCATE( isrcsng(mx_sng), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:isrcsng")
    ALLOCATE( icatsng(mx_sng), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:icatsng")
    ALLOCATE( latsng(mx_sng,ntime), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:latsng")
    ALLOCATE( lonsng(mx_sng,ntime), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:longsng")
    ALLOCATE( hgtsng(mx_sng,ntime), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:hgtsng")
    ALLOCATE( xsng(mx_sng), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:xsng")
    ALLOCATE( ysng(mx_sng), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:ysng")
    ALLOCATE( trnsng(mx_sng), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:trnsng")
    ALLOCATE( timesng(mx_sng,ntime), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:timesng")
    ALLOCATE( stnsng(mx_sng,ntime), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:stnsng")

    ALLOCATE( wx(mx_sng,ntime), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:wx")
    ALLOCATE( csrcsng(mx_sng,ntime), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:csrcsng")
    ALLOCATE( store_emv(mx_sng,5,ntime), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:store_emv")
    ALLOCATE( store_amt(mx_sng,5,ntime), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:store_amt")
    ALLOCATE( kloud(mx_sng,ntime), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:kloud")
    ALLOCATE( idp3(mx_sng,ntime), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:kdp3")
    ALLOCATE( store_hgt(mx_sng,5,ntime), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:store_hgt")
    ALLOCATE( obrdsng(mx_sng,nvar_sng,ntime), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:obrdsng")

    ALLOCATE( obsng  (nvar_anx,mx_sng), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:obsng")
    ALLOCATE( odifsng(nvar_anx,mx_sng), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:odifsng")
    ALLOCATE( oanxsng(nvar_anx,mx_sng), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:oanxsng")

    ALLOCATE( thesng(mx_sng), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:thesng")
    ALLOCATE( ival(mx_sng), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:ival")

!-----------------------------------------------------------------------
!
!  Upper Air Station variables
!
!-----------------------------------------------------------------------
!
    ALLOCATE( isrcua(mx_ua), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:isrcua")
    ALLOCATE( xua(mx_ua), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:xua")
    ALLOCATE( yua(mx_ua), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:yua")
    ALLOCATE( trnua(mx_ua), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:trnua")
    ALLOCATE( elevua(mx_ua), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:elevua")
    ALLOCATE( hgtua(nz_ua,mx_ua), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:hgtua")
    ALLOCATE( nlevsua(mx_ua), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:nlevsua")
    ALLOCATE( stnua(mx_ua), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:stnua")

!
!-----------------------------------------------------------------------
!
!  Upper-air observation variables
!
!-----------------------------------------------------------------------
!
    ALLOCATE( obsua (nvar_anx,nz_ua,mx_ua), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:obsua")
    ALLOCATE( odifua(nvar_anx,nz_ua,mx_ua), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:odifua")
    ALLOCATE( oanxua(nvar_anx,nz_ua,mx_ua), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:oanxua")
    ALLOCATE( theua(nz_ua,mx_ua), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:theua")
!
!-----------------------------------------------------------------------
!
!  Radar site variables
!
!  Radial velocity processing is less common then other data types.  If
!  we're not doing it, don't allocated unnecessary variables, especially
!  the really big ones.  The critical ones are stil left.
!
!-----------------------------------------------------------------------
!
    ALLOCATE( isrcrad(0:mx_rad), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:isrcrad")

    n_mx_rad     = mx_rad
    IF (raduvobs > 0 .OR. radrhobs > 0 ) THEN
!      n_mx_rad     = mx_rad
      n_mx_colrad  = mx_colrad
      n_nz_rdr     = nz_rdr
      n_nvar_radin = nvar_radin
      n_nvar_rad   = nvar_rad
    ELSE
!      n_mx_rad     = 1
      n_mx_colrad  = 1
      n_nz_rdr     = 1
      n_nvar_radin = 1
      n_nvar_rad   = 1
    END IF

    ALLOCATE( stnrad(n_mx_rad), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:stnrad")
    ALLOCATE( latrad(n_mx_rad), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:latrad")
    ALLOCATE(lonrad(n_mx_rad),  STAT = istatus )
    CALL check_alloc_status(istatus, "adas:longrad")
    ALLOCATE( elvrad(n_mx_rad), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:elevrad")
!
!-----------------------------------------------------------------------
!
! Radar observations
!
!-----------------------------------------------------------------------

    ALLOCATE( irad(n_mx_colrad) , STAT = istatus )
    CALL check_alloc_status(istatus, "adas:irad")
    ALLOCATE( nlevrad(n_mx_colrad) , STAT = istatus )
    CALL check_alloc_status(istatus, "adas:nlevrad")
    ALLOCATE( latradc(n_mx_colrad) , STAT = istatus )
    CALL check_alloc_status(istatus, "adas:latradc")
    ALLOCATE( lonradc(n_mx_colrad) , STAT = istatus )
    CALL check_alloc_status(istatus, "adas:lonradc")
    ALLOCATE( xradc(n_mx_colrad) , STAT = istatus )
    CALL check_alloc_status(istatus, "adas:xradc")
    ALLOCATE( yradc(n_mx_colrad) , STAT = istatus )
    CALL check_alloc_status(istatus, "adas:yradc")
    ALLOCATE( trnradc(n_mx_colrad) , STAT = istatus )
    CALL check_alloc_status(istatus, "adas:trnradc")
    ALLOCATE( distrad(n_mx_colrad) , STAT = istatus )
    CALL check_alloc_status(istatus, "adas:distrad")
    ALLOCATE( uazmrad(n_mx_colrad) , STAT = istatus )
    CALL check_alloc_status(istatus, "adas:uazmrad")
    ALLOCATE( vazmrad(n_mx_colrad) , STAT = istatus )
    CALL check_alloc_status(istatus, "adas:vazmrad")
    ALLOCATE( hgtradc(n_nz_rdr,n_mx_colrad) , STAT = istatus )
    CALL check_alloc_status(istatus, "adas:hgtradc")
    ALLOCATE( theradc(n_nz_rdr,n_mx_colrad) , STAT = istatus )
    CALL check_alloc_status(istatus, "adas:theradc")
    ALLOCATE( obsrad(n_nvar_radin,n_nz_rdr,n_mx_colrad) , STAT = istatus )
    CALL check_alloc_status(istatus, "adas:obsrad")

    ALLOCATE( odifrad(n_nvar_rad,n_nz_rdr,n_mx_colrad) , STAT = istatus )
    CALL check_alloc_status(istatus, "adas:odifrad")
    ALLOCATE( oanxrad(n_nvar_rad,n_nz_rdr,n_mx_colrad) , STAT = istatus )
    CALL check_alloc_status(istatus, "adas:oanxrad")

    ALLOCATE(    dsdr(n_nz_rdr,n_mx_colrad) , STAT = istatus )
    CALL check_alloc_status(istatus, "adas:dsdr")
    ALLOCATE(    dhdr(n_nz_rdr,n_mx_colrad) , STAT = istatus )
    CALL check_alloc_status(istatus, "adas:dhdr")

!-----------------------------------------------------------------------
!
! Different handling of ADAS and ARPS 3dvar
!
!-----------------------------------------------------------------------

    IF (progname == 'ADAS') THEN
      ALLOCATE( corsng(nvar_anx,mx_sng), STAT = istatus )
      CALL check_alloc_status(istatus, "adas:corsng")
      ALLOCATE(  corrad(nvar_rad,nz_rdr,mx_colrad) , STAT = istatus )
      CALL check_alloc_status(istatus, "adas:corrad")
      ALLOCATE( corua(nvar_anx,nz_ua,mx_ua), STAT = istatus )
      CALL check_alloc_status(istatus, "adas:corua")
    ELSE
      ALLOCATE( corsng(mx_sng,nvar_anx), STAT = istatus )
      CALL check_alloc_status(istatus, "adas:corsng")
      ALLOCATE( corrad(nz_rdr,mx_colrad,nvar_rad) , STAT = istatus )
      CALL check_alloc_status(istatus, "adas:corrad")
      ALLOCATE( corua(nz_ua,mx_ua,nvar_anx), STAT = istatus )
      CALL check_alloc_status(istatus, "adas:corua")
    END IF

!-----------------------------------------------------------------------
!
! Quality control variables
!
!-----------------------------------------------------------------------

    ALLOCATE( qualrdsng(mx_sng,nvar_sng,ntime), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:qualrdsng")
    ALLOCATE( qualsng(nvar_anx,mx_sng),    STAT = istatus )
    CALL check_alloc_status(istatus, "adas:qualsng")
    ALLOCATE( qobsng(nvar_anx,mx_sng),     STAT = istatus )
    CALL check_alloc_status(istatus, "adas:qobsng")

    ALLOCATE( qualua(nvar_anx,nz_ua,mx_ua),    STAT = istatus )
    CALL check_alloc_status(istatus, "adas:qualua")
    ALLOCATE(qobsua (nvar_anx,nz_ua,mx_ua),    STAT = istatus )
    CALL check_alloc_status(istatus, "adas:qobsua")

    ALLOCATE(  qualrad(n_nvar_rad,n_nz_rdr,n_mx_colrad), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:qualrad")
    ALLOCATE(  qobsrad(n_nvar_rad,n_nz_rdr,n_mx_colrad), STAT = istatus )
    CALL check_alloc_status(istatus, "adas:qobsrad")

    RETURN
  END SUBROUTINE allocate_adas_3dvar_arrays

  !
  ! DEALLOCATE
  !
  SUBROUTINE deallocate_adas_3dvar_arrays( istatus )

    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: istatus

!
!-----------------------------------------------------------------------
!
!  Surface Station variables
!
!-----------------------------------------------------------------------
!
    DEALLOCATE( isrcsng, STAT = istatus )
    DEALLOCATE( icatsng, STAT = istatus )
    DEALLOCATE( latsng,  STAT = istatus )
    DEALLOCATE( lonsng,  STAT = istatus )
    DEALLOCATE( hgtsng,  STAT = istatus )
    DEALLOCATE( xsng,    STAT = istatus )
    DEALLOCATE( ysng,    STAT = istatus )
    DEALLOCATE( trnsng,  STAT = istatus )
    DEALLOCATE( timesng, STAT = istatus )
    DEALLOCATE( stnsng,  STAT = istatus )
!
!-----------------------------------------------------------------------
!
!  Surface (single-level) read-in observation variables
!
!-----------------------------------------------------------------------
!
    DEALLOCATE( wx,        STAT = istatus )
    DEALLOCATE( csrcsng,   STAT = istatus )
    DEALLOCATE( store_emv, STAT = istatus )
    DEALLOCATE( store_amt, STAT = istatus )
    DEALLOCATE( kloud,     STAT = istatus )
    DEALLOCATE( store_hgt, STAT = istatus )
    DEALLOCATE( obrdsng,   STAT = istatus )

    DEALLOCATE( obsng  ,   STAT = istatus )
    DEALLOCATE( odifsng,   STAT = istatus )
    DEALLOCATE( oanxsng,   STAT = istatus )
    DEALLOCATE( corsng,    STAT = istatus )
    DEALLOCATE( thesng,    STAT = istatus )
    DEALLOCATE( ival,      STAT = istatus )

!-----------------------------------------------------------------------
!
!  Upper Air Station variables
!
!-----------------------------------------------------------------------
!
    DEALLOCATE( isrcua,  STAT = istatus )
    DEALLOCATE( xua,     STAT = istatus )
    DEALLOCATE( yua,     STAT = istatus )
    DEALLOCATE( trnua,   STAT = istatus )
    DEALLOCATE( elevua,  STAT = istatus )
    DEALLOCATE( hgtua,   STAT = istatus )
    DEALLOCATE( nlevsua, STAT = istatus )
    DEALLOCATE( stnua,   STAT = istatus )

!
!-----------------------------------------------------------------------
!
!  Upper-air observation variables
!
!-----------------------------------------------------------------------
!
    DEALLOCATE( obsua , STAT = istatus )
    DEALLOCATE( odifua, STAT = istatus )
    DEALLOCATE( oanxua, STAT = istatus )
    DEALLOCATE( corua,  STAT = istatus )
    DEALLOCATE( theua,  STAT = istatus )
!
!-----------------------------------------------------------------------
!
!  Radar site variables
!
!-----------------------------------------------------------------------
!
    DEALLOCATE( isrcrad, STAT = istatus )
    DEALLOCATE( stnrad,  STAT = istatus )
    DEALLOCATE( latrad,  STAT = istatus )
    DEALLOCATE( elvrad,  STAT = istatus )

!
!-----------------------------------------------------------------------
!
!  Radar observation variables
!
!-----------------------------------------------------------------------
    DEALLOCATE( irad,    STAT = istatus )
    DEALLOCATE( nlevrad, STAT = istatus )
    DEALLOCATE( latradc, STAT = istatus )
    DEALLOCATE( lonradc, STAT = istatus )
    DEALLOCATE( xradc,   STAT = istatus )
    DEALLOCATE( yradc,   STAT = istatus )
    DEALLOCATE( trnradc, STAT = istatus )
    DEALLOCATE( distrad, STAT = istatus )
    DEALLOCATE( uazmrad, STAT = istatus )
    DEALLOCATE( vazmrad, STAT = istatus )
    DEALLOCATE( hgtradc, STAT = istatus )
    DEALLOCATE( theradc, STAT = istatus )
    DEALLOCATE( obsrad,  STAT = istatus )

    DEALLOCATE( odifrad, STAT = istatus )
    DEALLOCATE( oanxrad, STAT = istatus )
    DEALLOCATE( corrad,  STAT = istatus )
    DEALLOCATE( dsdr,    STAT = istatus )
    DEALLOCATE( dhdr,    STAT = istatus )

!
!-----------------------------------------------------------------------
!
!  Quality Control Variables
!
!-----------------------------------------------------------------------
!
    DEALLOCATE( qualrdsng, STAT = istatus )
    DEALLOCATE( qualsng,   STAT = istatus )
    DEALLOCATE( qobsng,    STAT = istatus )

    DEALLOCATE( qualua,    STAT = istatus )
    DEALLOCATE( qobsua,    STAT = istatus )

    DEALLOCATE( qualrad,   STAT = istatus )
    DEALLOCATE( qobsrad,   STAT = istatus )

    RETURN
  END SUBROUTINE deallocate_adas_3dvar_arrays

END MODULE module_analysisArrays
