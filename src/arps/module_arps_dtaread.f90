MODULE arps_dtaread
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! This is a module designed to replace subroutine dtaread in arpsio3d.f90.
!
! It's main purpose is to provide convenience and flexibility to the
! pre-processing and post-process utilities in the ARPS package. For
! example, arpsenkf, arpsplt, arps4wrf, arps2coamps, arpscvt etc.
!
! This module provides the following flexibility:
!
!   o. Read from either split files or joined file in both MPI and
!      Non-MPI mode;
!   o. Separate base fields, total fields, as well other auxiliary fields
!      for performing minumum I/O;
!   o. Implicit opening and closing of files.
!
! Any program can use this module and follow the follow procedure:
!
!   1. Add USE statement as
!      USE arps_dtaread
!   2. Initialize this module for file conventions and format
!      CALL dtaread_init( ... )
!   3. Read field sets in any order, but dtaread_get_dim must always be
!      called first.
!
!      CALL dtaread_get_dim ( ... )     ! Get dimensions in the data file
!                                       ! patch dimensions
!
!      CALL dtaread_get_grid( ... )     ! Grid arrays
!      CALL dtaread_get_base( ... )     ! Base arrays
!      CALL dtaread_get_static( ... )   ! surface static fields
!      CALL dtaread_get_optional( ... ) ! radiation and flux arrays
!      CALL dtaread_get_3d( ... )       ! model forecast states
!
!   4. Terminate the module nicely by calling
!      CALL dtaread_final( ... )
!
! Note:
!
!   1. User should always initialize the module with two files even not
!      both of them may be read. One is the time-independent grid and base
!      file. The other is time-dependent history file.
!
!   2. This module has divided the ARPS arrays in the history file into
!      5 field sets and users can retrieve any number of sets in any order.
!
!      grid (x,y,z,zp,zpsoil)    retrieve from  grid and base file
!      static                    ...            grid and base file
!      base                      ...            grid and base file
!      3d state                  ...            time-dependent file
!      optional                  ...            time-dependent file
!
!   3. This module only initialize hydrometeor dimensions variables in
!      "globcst.inc" with "dtaread_get_dim". All other variables should
!      user's responsibility to initialize them explicitly.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: Yunheng Wang (05/24/2011)
!
! MODIFICATION HISOTRY:
!
!-----------------------------------------------------------------------

  USE arps_precision

  IMPLICIT NONE
  PRIVATE

  INTEGER :: dbglvl
  LOGICAL :: IAMROOT, IREAD, ISPLIT, IJOIN, IONE
  INTEGER :: ncompressx, ncompressy

  CHARACTER(LEN=256) :: basfile, hisfile
  INTEGER            :: filefmt
  INTEGER, ALLOCATABLE :: fileid(:,:)
  LOGICAL :: basopened, hisopened

  INTEGER :: nxd,nyd,nzd,nzsoild,nstypsd,nscalard
  INTEGER :: nxin, nyin
  CHARACTER(LEN=4) :: q_names(20)

  INCLUDE 'mp.inc'
!-----------------------------------------------------------------------
!
! Public interface
!
!-----------------------------------------------------------------------

  PUBLIC :: dtaread_init, dtaread_final
  PUBLIC :: dtaread_get_dim, dtaread_get_head, dtaread_get_dsd
  PUBLIC :: dtaread_get_grid, dtaread_get_static, dtaread_get_base,     &
            dtaread_get_3d, dtaread_get_optional, dtaread_get_var

  CONTAINS

  !##################################################################

  SUBROUTINE dtaread_init( grdbas,hfile,hinfmt,lvldbg,istatus )
  !
  !-----------------------------------------------------------------------
  !
  !  PURPOSE:
  !
  !  General initialize for this module. It determines how the file will
  !  read. There are the following variations:
  !                                                              IREAD   ISPLIT   IJOIN  IONE
  !  1. Non-mpi read one file;                                   .T.     .F.      .F.    .T.
  !  2. Non-mpi read multiple patches;                           .T.     .F.      .T.    .F.
  !  3. Root process read one file and then split on-the-fly.    .X.     .T.      .F.    .T.
  !  4. All processes read its own patches;                      .T.     .F.      .F.    .F.
  !  5. All processes read and join several small patches        .T.     .F.      .T.    .F.
  !
  !  NOTE: module variable basfile & hisfile are just the file name base
  !  and the exact file name is determined later while reading.
  !
  !-----------------------------------------------------------------------
  !
  !  Author: Yunheng Wang
  !  Date: 05/31/11
  !
  !  MODIFICATION HISTORY:
  !
  !
  !-----------------------------------------------------------------------
  !
    IMPLICIT NONE
    CHARACTER(LEN=256), INTENT(IN)  :: grdbas, hfile
    INTEGER,            INTENT(IN)  :: hinfmt
    INTEGER, INTENT(IN)  :: lvldbg
    INTEGER, INTENT(OUT) :: istatus
  !
  !-----------------------------------------------------------------------
  !
  !  Variable Declarations.
  !
  !-----------------------------------------------------------------------
  !


  !
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !
    istatus = 0

    dbglvl = lvldbg

    IAMROOT = .FALSE.
    IF (myproc == 0) IAMROOT = .TRUE.

  !---------------------------------------------------------------------

    IREAD  = .FALSE.
    ISPLIT = .FALSE.
    IJOIN  = .FALSE.
    IONE   = .TRUE.
    IF (mp_opt == 0) THEN     ! no-mpi specific
      ncompressx = nproc_x_in
      ncompressy = nproc_y_in
      IREAD = .TRUE.
    ELSE                      ! mpi specific
      readstride = max_fopen

      IF (nproc_x_in == 1 .AND. nproc_y_in == 1) THEN
        readsplit = 1
        readstride = nprocs

        ncompressx = 1
        ncompressy = 1
        IF (myproc == 0) IREAD = .TRUE.
        ISPLIT = .TRUE.
      ELSE
        readsplit = 0

        ncompressx = nproc_x_in/nproc_x
        ncompressy = nproc_y_in/nproc_y

        IF ( (MOD(nproc_x_in,nproc_x) /= 0 .OR. MOD(nproc_y_in,nproc_y) /= 0)   &
             .OR. (ncompressx < 1 .OR. ncompressy < 1)  ) THEN
          IF (myproc == 0) WRITE(6,'(3x,a/,2(3x,2(a,I2)/))')                    &
            'nprocx_in (nprocy_in) must be multiples of nproc_x(nproc_y)',      &
            'nprocx_in = ',nproc_x_in, ',nprocy_in = ',nproc_y_in,              &
            '  nproc_x = ',nproc_x,  ',  nproc_y = ', nproc_y
          CALL mpexit(1);
          STOP
        END IF

        IREAD = .TRUE.
        IONE  = .FALSE.
      END IF
    END IF

    IF (ncompressx > 1 .OR. ncompressy > 1) THEN
      IJOIN = .TRUE.
      IONE  = .FALSE.
    END IF

  !---------------------------------------------------------------------
  !
  ! Set file name base
  !
  !---------------------------------------------------------------------

    filefmt = hinfmt

    IF (filefmt /= 3) THEN
      WRITE(*,'(1x,a,I2)') 'ERROR: unsupported file format at present, filefmt = ',filefmt
      CALL arpsstop('ERROR: unsupported file format.',1)
    END IF

    ALLOCATE(fileid(ncompressx,ncompressy), STAT = istatus)

    basfile = grdbas
    hisfile = hfile

    basopened = .FALSE.
    hisopened = .FALSE.

    RETURN
  END SUBROUTINE dtaread_init

  !##################################################################

  SUBROUTINE dtaread_final( istatus )
  !
  !-----------------------------------------------------------------------
  !
  !  PURPOSE:
  !
  !  Close opened files and release this module
  !
  !-----------------------------------------------------------------------
  !
  !  Author: Yunheng Wang
  !  Date: 06/01/11
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


  !
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !
    istatus = 0

    IF (basopened .OR. hisopened) THEN
      CALL close_files( istatus )
    END IF

    DEALLOCATE(fileid)

    RETURN
  END SUBROUTINE dtaread_final

  !##################################################################

  SUBROUTINE dtaread_get_dim( nx,ny,nz,nzsoil,nstyps,istatus )
  !
  !---------------------------------------------------------------------
  !
  !  PURPOSE:
  !
  !  Retrieve data dimenions after file name is set.
  !
  !---------------------------------------------------------------------
  !
  !  Author: Yunheng Wang
  !  Date: 05/31/11
  !
  !  MODIFICATION HISTORY:
  !
  !
  !---------------------------------------------------------------------
  !
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: nx,ny,nz,nzsoil,nstyps
    INTEGER, INTENT(OUT) :: istatus
  !
  !---------------------------------------------------------------------
  !
  !  Variable Declarations.
  !
  !---------------------------------------------------------------------
  !
    CHARACTER(LEN=256) :: filename
    INTEGER :: nq

    INCLUDE 'globcst.inc'
  !
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !
    istatus = 0

    IF (IAMROOT) THEN   ! Suppose all files contains the same dimensions
      IF (IONE) THEN
        filename = hisfile
      ELSE
        CALL gtsplitfn(hisfile,ncompressx,ncompressy,loc_x,loc_y,1,1,   &
                     0,0,1,dbglvl,filename,istatus)
      END IF

      SELECT CASE (filefmt)
      CASE (3)
        CALL hdf_getdims(filename,nxin,nyin,nz,nzsoil,nstyps,istatus)
      END SELECT
    END IF
    CALL mpupdatei(istatus,1)

    CALL mpupdatei(nxin,1)
    CALL mpupdatei(nyin,1)
    CALL mpupdatei(nz,1)
    CALL mpupdatei(nzsoil,1)
    CALL mpupdatei(nstyps,1)
    CALL mpupdatei(nscalar,1)
    CALL mpupdatei(nscalarq,1)
    CALL mpupdatec(qnames,40*20)
    CALL mpupdatec(qdescp,40*20)

    CALL mpupdatei(P_QC, 1)
    CALL mpupdatei(P_QR, 1)
    CALL mpupdatei(P_QI, 1)
    CALL mpupdatei(P_QS, 1)
    CALL mpupdatei(P_QG, 1)
    CALL mpupdatei(P_QH, 1)

    CALL mpupdatei(P_NC, 1)
    CALL mpupdatei(P_NR, 1)
    CALL mpupdatei(P_NI, 1)
    CALL mpupdatei(P_NS, 1)
    CALL mpupdatei(P_NG, 1)
    CALL mpupdatei(P_NH, 1)

    CALL mpupdatei(P_ZR, 1)
    CALL mpupdatei(P_ZI, 1)
    CALL mpupdatei(P_ZS, 1)
    CALL mpupdatei(P_ZG, 1)
    CALL mpupdatei(P_ZH, 1)

    IF (ISPLIT) THEN
      nxd = (nxin-3)/nproc_x+3
      nyd = (nyin-3)/nproc_y+3
    ELSE
      nxd = (nxin-3)*ncompressx+3
      nyd = (nyin-3)*ncompressy+3
    END IF
    nzd = nz
    nzsoild = nzsoil
    nstypsd = nstyps

    nscalard = nscalar

    DO nq = 1, nscalard
      q_names(nq) = qnames(nq)
    END DO

    nx = nxd       ! Return patch size
    ny = nyd
    RETURN
  END SUBROUTINE dtaread_get_dim

  !##################################################################

  SUBROUTINE dtaread_get_grid(x,y,z,zp,zpsoil,dx,dy,dz,istatus)

  !
  !---------------------------------------------------------------------
  !
  !  PURPOSE:
  !
  !  Reads in grid arrays from the grid and base file.
  !
  !---------------------------------------------------------------------
  !
  !  AUTHOR: Yunheng Wang
  !    06/02/2011.
  !
  !  MODIFICATION HISTORY:
  !
  !  DATA ARRAYS READ IN:
  !
  !    x        x coordinate of grid points in physical/comp. space (m)
  !    y        y coordinate of grid points in physical/comp. space (m)
  !    z        z coordinate of grid points in computational space (m)
  !    zp       z coordinate of grid points in physical space (m)
  !    zpsoil   z coordinate of grid points in the soil (m)
  !
  !  OUTPUT:
  !
  !    dx       grid resolution in X direction
  !    dy       grid resolution in Y direction
  !    dz
  !    ireturn  Return status indicator
  !
  !---------------------------------------------------------------------

    USE arps_precision

    IMPLICIT NONE

  !---------------------------------------------------------------------
  !
  !  Variable Declarations.
  !
  !---------------------------------------------------------------------
  !
    REAL(P), INTENT(OUT) :: x (nxd)           ! x-coord. of the physical and compu
                                              ! -tational grid. Defined at u-point(m).
    REAL(P), INTENT(OUT) :: y (nyd)           ! y-coord. of the physical and compu
                                              ! -tational grid. Defined at v-point(m).
    REAL(P), INTENT(OUT) :: z (nzd)           ! z-coord. of the computational grid.
                                              ! Defined at w-point on the staggered
                                              ! grid(m).
    REAL(P), INTENT(OUT) :: zp(nxd,nyd,nzd)   ! Physical height coordinate defined at
                                              ! w-point of the staggered grid(m).
    REAL(P), INTENT(OUT) :: zpsoil(nxd,nyd,nzsoild)   ! Physical height coordinate defined at
                                                      ! w-point of the soil (m).
    REAL(P), INTENT(OUT) :: dx, dy, dz

    INTEGER, INTENT(OUT) :: istatus

  !---------------------------------------------------------------------
  !
  ! Misc. Local variables
  !
  !---------------------------------------------------------------------

    REAL(P) :: amin, amax
    INTEGER :: i,j,k
  !
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !
  !  Beginning of executable code...
  !
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !
    istatus = 0

    IF (.NOT. basopened) CALL open_files( 0,istatus )

    SELECT CASE (filefmt)
    CASE ( 1 )

      ! To be implemented

    CASE( 3 )  ! HDF 4 Format

      CALL hdf_read_grid(IAMROOT, IREAD, ISPLIT, ncompressx, ncompressy,&
                         fileid,nxd,nyd,nzd,nzsoild,nxin,nyin,          &
                         x,y,z,zp,zpsoil,istatus)

    CASE (7)    ! NetCDF format

      ! To be implemented

    CASE DEFAULT
      WRITE(*,'(1x,a,I2)') 'ERROR: unsupported file format at present, filefmt = ',filefmt
      CALL arpsstop('ERROR: unsupported file format.',1)
    END SELECT

    IF(IAMROOT .AND. istatus == 0)  WRITE(6,'(/3a/)')                   &
           ' Grid arrays from ',TRIM(basfile),' were successfully read.'
  !
  !-----------------------------------------------------------------------
  !
  !  Set dx, dy and dz from read-in data
  !
  !-----------------------------------------------------------------------
  !
    dx = x(2) - x(1)
    dy = y(2) - y(1)
    dz = z(2) - z(1)

    IF (dbglvl > 90) THEN
      IF (IAMROOT)  &
        WRITE(6,'(/1x,a/)') 'Min. and max. of the grid arrays read in:'

      CALL edgfill(x,1,nxd,1,nxd,1,1,1,1,1,1,1,1)
      CALL a3dmax0(x,1,nxd,1,nxd,1,1,1,1,1,1,1,1,amax,amin)
      IF (IAMROOT)  &
         WRITE(6,'(/1x,2(a,e13.6))') 'xmin    = ', amin,',  xmax    =',amax

      CALL edgfill(y,1,nyd,1,nyd,1,1,1,1,1,1,1,1)
      CALL a3dmax0(y,1,nyd,1,nyd,1,1,1,1,1,1,1,1,amax,amin)
      IF (IAMROOT)  &
         WRITE(6,'(1x,2(a,e13.6))') 'ymin    = ', amin,',  ymax    =',amax

      CALL edgfill(z,1,nzd,1,nzd,1,1,1,1,1,1,1,1)
      CALL a3dmax0(z,1,nzd,1,nzd,1,1,1,1,1,1,1,1,amax,amin)
      IF (IAMROOT)  &
         WRITE(6,'(1x,2(a,e13.6))') 'zmin    = ', amin,',  zmax    =',amax

      CALL edgfill(zp,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,nzd,1,nzd)
      CALL a3dmax0(zp,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,nzd,1,nzd,amax,amin)
      IF (IAMROOT)  &
         WRITE(6,'(1x,2(a,e13.6))') 'zpmin   = ', amin,',  zpmax   =',amax

      CALL edgfill(zpsoil,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,nzsoild,1,nzsoild)
      CALL a3dmax0(zpsoil,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,nzsoild,1,nzsoild,amax,amin)
      IF (IAMROOT)  &
         WRITE(6,'(1x,2(a,e13.6))') 'zpsoilmin   = ', amin,',  zpsoilmax   =',amax

    END IF

    RETURN
  END SUBROUTINE dtaread_get_grid

  !##################################################################

  SUBROUTINE dtaread_get_base(ubar,vbar,wbar,ptbar,pbar,rhobar,qvbar,   &
                              istatus)

  !
  !---------------------------------------------------------------------
  !
  !  PURPOSE:
  !
  !   Read base state from time independent grid and base file.
  !
  !---------------------------------------------------------------------
  !
  !  AUTHOR: Yunheng Wang (06/01/2011).
  !
  !  MODIFICATION HISTORY:
  !
  !---------------------------------------------------------------------
  !  DATA ARRAYS READ IN:
  !
  !    ubar     Base state x velocity component (m/s)
  !    vbar     Base state y velocity component (m/s)
  !    wbar     Base state z velocity component (m/s)
  !    ptbar    Base state potential temperature (K)
  !    pbar     Base state pressure (Pascal)
  !    rhobar   Base state air density (kg/m**3)
  !    qvbar    Base state water vapor mixing ratio (kg/kg)
  !
  !---------------------------------------------------------------------

    USE arps_precision

    IMPLICIT NONE

  !---------------------------------------------------------------------
  !
  !  Variable Declarations.
  !
  !---------------------------------------------------------------------
  !
    REAL(P), INTENT(OUT) :: ubar  (nxd,nyd,nzd)     ! Base state x velocity component (m/s)
    REAL(P), INTENT(OUT) :: vbar  (nxd,nyd,nzd)     ! Base state y velocity component (m/s)
    REAL(P), INTENT(OUT) :: wbar  (nxd,nyd,nzd)     ! Base state z velocity component (m/s)
    REAL(P), INTENT(OUT) :: ptbar (nxd,nyd,nzd)     ! Base state potential temperature (K)
    REAL(P), INTENT(OUT) :: rhobar(nxd,nyd,nzd)     ! Base state air density (kg/m**3)
    REAL(P), INTENT(OUT) :: pbar  (nxd,nyd,nzd)     ! Base state pressure (Pascal)
    REAL(P), INTENT(OUT) :: qvbar (nxd,nyd,nzd)     ! Base state water vapor mixing ratio
                                       ! (kg/kg)
    INTEGER, INTENT(OUT) :: istatus

  !---------------------------------------------------------------------
  !
  ! Misc. Local variables
  !
  !---------------------------------------------------------------------

    REAL(P) :: p0inv, amin, amax
    INTEGER :: i,j,k

    INCLUDE 'phycst.inc'
  !
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !
  !  Beginning of executable code...
  !
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !
    istatus = 0

    p0inv=1./p0
  !
  !
  !---------------------------------------------------------------------
  !
  !  Read grid and base state fields.
  !
  !---------------------------------------------------------------------
  !
    IF (.NOT. basopened) CALL open_files( 0,istatus )

    SELECT CASE (filefmt)
    CASE ( 1 )

      ! To be implemented

    CASE( 3 )  ! HDF 4 Format

      CALL hdf_read_base(IAMROOT, IREAD, ISPLIT, ncompressx, ncompressy,&
           fileid,nxd,nyd,nzd,nxin,nyin,                                &
           ubar, vbar, wbar, ptbar, pbar, qvbar,                        &
           istatus)

    CASE (7)    ! NetCDF format

      ! To be implemented

    CASE DEFAULT
      WRITE(*,'(1x,a,I2)') 'ERROR: unsupported file format at present, filefmt = ',filefmt
      CALL arpsstop('ERROR: unsupported file format.',1)
    END SELECT

    IF (istatus /= 0) THEN
      WRITE(*,'(1x,a,I2)') 'ERROR: base file reading status = ',istatus
      CALL arpsstop('ERROR: base file problem.',1)
    ELSE
      IF(IAMROOT)  WRITE(6,'(/3a/)')                                    &
           ' Base states from ',TRIM(basfile),' were successfully read.'
    END IF

    DO k=1,nzd-1
      DO j=1,nyd-1
        DO i=1,nxd-1
          rhobar(i,j,k)= pbar(i,j,k)/                                   &
                 ( rd * ptbar(i,j,k)*(pbar(i,j,k)*p0inv)**rddcp )
        END DO
      END DO
    END DO

    IF (dbglvl > 90) THEN

      IF (IAMROOT)  &
        WRITE(6,'(/1x,a/)') 'Min. and max. of the base arrays read in:'

      CALL edgfill(ubar,1,nxd,1,nxd,1,nyd,1,nyd-1,1,nzd,1,nzd-1)
      CALL a3dmax0(ubar,1,nxd,1,nxd,1,nyd,1,nyd-1,1,nzd,1,nzd-1,amax,amin)
      IF (IAMROOT)  &
         WRITE(6,'(1x,2(a,e13.6))') 'ubarmin = ', amin,',  ubarmax =',amax

      CALL edgfill(vbar,1,nxd,1,nxd-1,1,nyd,1,nyd,1,nzd,1,nzd-1)
      CALL a3dmax0(vbar,1,nxd,1,nxd-1,1,nyd,1,nyd,1,nzd,1,nzd-1,amax,amin)
      IF (IAMROOT)  &
         WRITE(6,'(1x,2(a,e13.6))') 'vbarmin = ', amin,',  vbarmax =',amax

      CALL edgfill(wbar,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,nzd,1,nzd)
      CALL a3dmax0(wbar,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,nzd,1,nzd,amax,amin)
      IF (IAMROOT)  &
         WRITE(6,'(1x,2(a,e13.6))') 'wbarmin = ', amin,',  wbarmax =',amax

      CALL edgfill(ptbar,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,nzd,1,nzd-1)
      CALL a3dmax0(ptbar,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,nzd,1,nzd-1,amax,amin)
      IF (IAMROOT)  &
         WRITE(6,'(1x,2(a,e13.6))') 'ptbarmin= ', amin,',  ptbarmax=',amax

      CALL edgfill(pbar,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,nzd,1,nzd-1)
      CALL a3dmax0(pbar,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,nzd,1,nzd-1,amax,amin)
      IF (IAMROOT)  &
         WRITE(6,'(1x,2(a,e13.6))') 'pbarmin = ', amin,',  pbarmax =',amax

      CALL edgfill(qvbar,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,nzd,1,nzd-1)
      CALL a3dmax0(qvbar,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,nzd,1,nzd-1,amax,amin)
      IF (IAMROOT)  &
         WRITE(6,'(1x,2(a,e13.6))') 'qvbarmin= ', amin, ',  qvbarmax=',amax

    END IF

    RETURN
  END SUBROUTINE dtaread_get_base

  !##################################################################

  SUBROUTINE dtaread_get_static(landin,soiltyp,stypfrct,vegtyp,lai,roufns,veg, &
                                istatus)

  !
  !---------------------------------------------------------------------
  !
  !  PURPOSE:
  !
  !  Reads in surface static arrays from the grid and base file.
  !
  !---------------------------------------------------------------------
  !
  !  AUTHOR: Yunheng Wang
  !    06/02/2011.
  !
  !  MODIFICATION HISTORY:
  !
  !---------------------------------------------------------------------
  !
  !  DATA ARRAYS READ IN:
  !
  !    soiltyp  Soil type
  !    stypfrct  Soil type fraction
  !    vegtyp   Vegetation type
  !    lai      Leaf Area Index
  !    roufns   Surface roughness
  !    veg      Vegetation fraction
  !
  !---------------------------------------------------------------------

    USE arps_precision

    IMPLICIT NONE

  !---------------------------------------------------------------------

    INTEGER, INTENT(OUT) :: landin
    INTEGER, INTENT(OUT) :: soiltyp (nxd,nyd,nstypsd)    ! Soil type
    REAL(P), INTENT(OUT) :: stypfrct(nxd,nyd,nstypsd)    ! Soil type
    INTEGER, INTENT(OUT) :: vegtyp  (nxd,nyd)            ! Vegetation type
    REAL(P), INTENT(OUT) :: lai     (nxd,nyd)            ! Leaf Area Index
    REAL(P), INTENT(OUT) :: roufns  (nxd,nyd)            ! Surface roughness
    REAL(P), INTENT(OUT) :: veg     (nxd,nyd)            ! Vegetation fraction

    INTEGER, INTENT(OUT) :: istatus

  !---------------------------------------------------------------------
  !
  ! Misc. Local variables
  !
  !---------------------------------------------------------------------

    REAL(P) :: amin, amax
    INTEGER :: imin, imax
    INTEGER :: is
  !
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !
  !  Beginning of executable code...
  !
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !
    istatus = 0

    IF (.NOT. basopened) CALL open_files( 0,istatus )

    SELECT CASE (filefmt)
    CASE ( 1 )

      ! To be implemented

    CASE( 3 )  ! HDF 4 Format

      CALL hdf_read_static(IAMROOT, IREAD, ISPLIT, ncompressx, ncompressy,&
                 fileid,nxd,nyd,nstypsd,nxin,nyin,                        &
                 landin,soiltyp,stypfrct,vegtyp,lai,roufns,veg,istatus)

    CASE (7)    ! NetCDF format

      ! To be implemented

    CASE DEFAULT
      WRITE(*,'(1x,a,I2)') 'ERROR: unsupported file format at present, filefmt = ',filefmt
      CALL arpsstop('ERROR: unsupported file format.',1)
    END SELECT

    IF(IAMROOT .AND. istatus == 0)  WRITE(6,'(/3a/)')                 &
           ' Static arrays from ',TRIM(basfile),' were successfully read.'

    IF (dbglvl > 90) THEN
      IF (IAMROOT)  &
        WRITE(6,'(/1x,a/)') 'Min. and max. of the static arrays read in:'

      CALL iedgfill(soiltyp,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,nstypsd,1,nstypsd)
      imax = MAXVAL(soiltyp)
      imin = MINVAL(soiltyp)
      CALL mpmaxi(imax)
      CALL mpmini(imin)
      IF (IAMROOT)  WRITE(6,'(1x,2(a,tr1,I13))')                        &
          'soiltypmin = ', imin,',  soiltypmax = ',imax

      CALL edgfill(stypfrct,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,nstypsd,1,nstypsd)
      CALL a3dmax0(stypfrct,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,nstypsd,1,nstypsd,amax,amin)
      IF (IAMROOT)  WRITE(6,'(1x,2(a,tr1,e13.4))')                      &
          'stypfrctmin = ', amin,',  stypfrctmax = ',amax

      CALL iedgfill(vegtyp,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,1,1,1)
      imax = MAXVAL(vegtyp)
      imin = MINVAL(vegtyp)
      CALL mpmaxi(imax)
      CALL mpmini(imin)
      IF (IAMROOT)  WRITE(6,'(1x,2(a,tr1,I13))')                        &
          'vegtypmin = ', imin,',  vegtypmax = ',imax

      CALL edgfill(lai,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,1,1,1)
      CALL a3dmax0(lai,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,1,1,1,amax,amin)
      IF (IAMROOT)  WRITE(6,'(1x,2(a,tr1,e13.4))')                      &
          'laimin = ', amin,',  laimax = ',amax

      CALL edgfill(roufns,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,1,1,1)
      CALL a3dmax0(roufns,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,1,1,1,amax,amin)
      IF (IAMROOT)  WRITE(6,'(1x,2(a,tr1,e13.4))')                      &
          'roufnsmin = ', amin,',  roufnsmax = ',amax

      CALL edgfill(veg,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,1,1,1)
      CALL a3dmax0(veg,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,1,1,1,amax,amin)
      IF (IAMROOT)  WRITE(6,'(1x,2(a,tr1,e13.4))')                      &
          'vegmin = ', amin,',  vegmax = ',amax

    END IF

    RETURN
  END SUBROUTINE dtaread_get_static

  !##################################################################

  SUBROUTINE dtaread_get_3d( time,u,v,w,pt, pp ,                        &
                             qv, qscalar, tke, kmh,kmv,                 &
                             tsoil,qsoil,wetcanp,snowdpth,              &
                             raing,rainc,prcrate,                       &
                             istatus)
  !
  !---------------------------------------------------------------------
  !
  !  PURPOSE:
  !
  !  Coordinate the reading of history data of various formats.
  !
  !---------------------------------------------------------------------
  !
  !  AUTHOR: Yunheng Wang  (5/31/2011).
  !
  !  MODIFICATION HISTORY:
  !
  !
  !---------------------------------------------------------------------
  !
  !  OUTPUT:
  !
  !    time     The time of the input data (s)
  !
  !    u     x component of total velocity (m/s)
  !    v     y component of total velocity (m/s)
  !    w     vertical component of total velocity
  !             in Cartesian coordinates (m/s).
  !
  !    pt    total potential temperature (K)
  !    pp    total pressure (Pascal)
  !    qv    total water vapor mixing ratio (kg/kg)
  !
  !    qscalar  hydrometeors
  !    tke      Turbulent Kinetic Energy ((m/s)**2)
  !
  !    kmh      Horizontal turb. mixing coef. for momentum ( m**2/s )
  !    kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
  !
  !    tsoil    Soil temperature (K)
  !    qsoil    Soil moisture (m**3/m**3)
  !    wetcanp  Canopy water amount
  !
  !    raing    Grid supersaturation rain
  !    rainc    Cumulus convective rain
  !    prcrate  Precipitation rates
  !
  !---------------------------------------------------------------------

    IMPLICIT NONE

    REAL(P), INTENT(OUT) :: time                 ! The time of the input data (s)

    REAL(P), INTENT(OUT) :: u (nxd,nyd,nzd)      ! x component of perturbation velocity
                                                 ! (m/s)
    REAL(P), INTENT(OUT) :: v (nxd,nyd,nzd)      ! y component of perturbation velocity
                                                 ! (m/s)
    REAL(P), INTENT(OUT) :: w (nxd,nyd,nzd)      ! vertical component of perturbation
                                                 ! velocity in Cartesian coordinates
                                                 ! (m/s).
    REAL(P), INTENT(OUT) :: pt(nxd,nyd,nzd)      ! perturbation potential temperature (K)
    REAL(P), INTENT(OUT) :: pp(nxd,nyd,nzd)      ! perturbation pressure (Pascal)
    REAL(P), INTENT(OUT) :: qv(nxd,nyd,nzd)      ! perturbation water vapor mixing ratio
                                                 ! (kg/kg)
    REAL(P), INTENT(OUT) :: qscalar(nxd,nyd,nzd,nscalard)

    REAL(P), INTENT(OUT) :: tke  (nxd,nyd,nzd)      ! Turbulent Kinetic Energy ((m/s)**2)

    REAL(P), INTENT(OUT) :: kmh   (nxd,nyd,nzd)     ! Horizontal turb. mixing coef. for
                                                    ! momentum. ( m**2/s )
    REAL(P), INTENT(OUT) :: kmv   (nxd,nyd,nzd)     ! Vertical turb. mixing coef. for
                                                    ! momentum. ( m**2/s )
                                                    ! (kg/kg)

    REAL(P), INTENT(OUT) :: tsoil  (nxd,nyd,nzsoild,0:nstypsd) ! Soil temperature (K)
    REAL(P), INTENT(OUT) :: qsoil  (nxd,nyd,nzsoild,0:nstypsd) ! Soil moisture (m**3/m**3)
    REAL(P), INTENT(OUT) :: wetcanp(nxd,nyd,0:nstypsd)         ! Canopy water amount
    REAL(P), INTENT(OUT) :: snowdpth(nxd,nyd)                  ! Snow depth (m)

    REAL(P), INTENT(OUT) :: raing(nxd,nyd)         ! Grid supersaturation rain
    REAL(P), INTENT(OUT) :: rainc(nxd,nyd)         ! Cumulus convective rain
    REAL(P), INTENT(OUT) :: prcrate(nxd,nyd,4)     ! precipitation rate (kg/(m**2*s))
                                   ! prcrate(1,1,1) = total precip. rate
                                   ! prcrate(1,1,2) = grid scale precip. rate
                                   ! prcrate(1,1,3) = cumulus precip. rate
                                   ! prcrate(1,1,4) = microphysics precip. rate

    INTEGER, INTENT(OUT) :: istatus

  !---------------------------------------------------------------------
  !
  !  Variable Declarations.
  !
  !---------------------------------------------------------------------
  !
    INTEGER :: nq,is
    REAL(P) :: amin, amax
  !
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !
  !  Beginning of executable code...
  !
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !

    IF (.NOT. hisopened) CALL open_files( 1,istatus )

    SELECT CASE (filefmt)
    CASE ( 1 )

      ! To be implemented

    CASE( 3 )  ! HDF 4 Format

      CALL hdf_read_3d(IAMROOT, IREAD, ISPLIT, ncompressx, ncompressy,  &
           fileid,nxd,nyd,nzd,nzsoild,nstypsd,nscalard,q_names,nxin,nyin,&
           time,u, v, w, pt, pp,                                        &
           qv, qscalar, tke, kmh,kmv,                                   &
           tsoil,qsoil,wetcanp,snowdpth,                                &
           raing,rainc,prcrate,dbglvl,                                  &
           istatus)

    CASE (7)    ! NetCDF format

      ! To be implemented

    CASE DEFAULT
      WRITE(*,'(1x,a,I2)') 'ERROR: unsupported file format at present, filefmt = ',filefmt
      CALL arpsstop('ERROR: unsupported file format.',1)
    END SELECT

    IF(IAMROOT .AND. istatus == 0)  WRITE(6,'(/a,F8.1,a/)')             &
      ' Data at time=', time/60,' (min) were successfully read.'

    IF (dbglvl > 90) THEN
      IF (IAMROOT)  &
        WRITE(6,'(/1x,a/)') 'Min. and max. of the data arrays read in:'

      CALL edgfill(u,1,nxd,1,nxd,1,nyd,1,nyd-1,1,nzd,1,nzd-1)
      CALL a3dmax0(u,1,nxd,1,nxd,1,nyd,1,nyd-1,1,nzd,1,nzd-1,amax,amin)
      IF (IAMROOT)  &
         WRITE(6,'(1x,2(a,e13.6))') 'umin = ', amin,',  umax =',amax

      CALL edgfill(v,1,nxd,1,nxd-1,1,nyd,1,nyd,1,nzd,1,nzd-1)
      CALL a3dmax0(v,1,nxd,1,nxd-1,1,nyd,1,nyd,1,nzd,1,nzd-1,amax,amin)
      IF (IAMROOT)  &
         WRITE(6,'(1x,2(a,e13.6))') 'vmin = ', amin,',  vmax =',amax

      CALL edgfill(w,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,nzd,1,nzd)
      CALL a3dmax0(w,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,nzd,1,nzd,amax,amin)
      IF (IAMROOT)  &
         WRITE(6,'(1x,2(a,e13.6))') 'wmin = ', amin,',  wmax =',amax

      CALL edgfill(pt,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,nzd,1,nzd-1)
      CALL a3dmax0(pt,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,nzd,1,nzd-1,amax,amin)
      IF (IAMROOT)  &
         WRITE(6,'(1x,2(a,e13.6))') 'ptmin= ', amin,',  ptmax=',amax

      CALL edgfill(pp,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,nzd,1,nzd-1)
      CALL a3dmax0(pp,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,nzd,1,nzd-1,amax,amin)
      IF (IAMROOT)  &
         WRITE(6,'(1x,2(a,e13.6))') 'pmin = ', amin,',  pmax =',amax

      CALL edgfill(qv,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,nzd,1,nzd-1)
      CALL a3dmax0lcl(qv,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,nzd,1,nzd-1,amax,amin)
      IF (IAMROOT) WRITE(6,'(1x,2(a,e13.6))')                           &
                               'qvmin= ', amin, ',  qvmax= ', amax

      DO nq = 1,nscalard
        CALL edgfill   (qscalar(:,:,:,nq),1,nxd,1,nxd-1,1,nyd,1,nyd-1,  &
                        1,nzd,1,nzd-1)
        CALL a3dmax0lcl(qscalar(:,:,:,nq),1,nxd,1,nxd-1,1,nyd,1,nyd-1,  &
                        1,nzd,1,nzd-1,amax,amin)

        IF (IAMROOT)  &
           WRITE(6,'(1x,2(a,e13.6))') TRIM(q_names(nq))//'min   = ', amin, &
                               ',  '//TRIM(q_names(nq))//'max   = ', amax
      END DO

      CALL edgfill(raing,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,1,1,1)
      CALL a3dmax0(raing,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,1,1,1,amax,amin)
      IF (IAMROOT)  &
         WRITE(6,'(1x,2(a,e13.6))')                                     &
          'raingmin= ', amin,',  raingmax=',amax

      CALL edgfill(rainc,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,1,1,1)
      CALL a3dmax0(rainc,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,1,1,1,amax,amin)
      IF (IAMROOT)  &
         WRITE(6,'(1x,2(a,e13.6))')                                     &
          'raincmin= ', amin,',  raincmax=',amax

      CALL edgfill(prcrate(1,1,1),1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,1,1,1)
      CALL a3dmax0(prcrate(1,1,1),1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,1,1,1,amax,amin)
      IF (IAMROOT)  &
         WRITE(6,'(1x,2(a,e13.6))')                                     &
          'prcr1min= ', amin,',  prcr1max=',amax

      CALL edgfill(prcrate(1,1,2),1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,1,1,1)
      CALL a3dmax0(prcrate(1,1,2),1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,1,1,1,amax,amin)
      IF (IAMROOT)  &
         WRITE(6,'(1x,2(a,e13.6))')                                     &
          'prcr2min= ', amin,',  prcr2max=',amax

      CALL edgfill(prcrate(1,1,3),1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,1,1,1)
      CALL a3dmax0(prcrate(1,1,3),1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,1,1,1,amax,amin)
      IF (IAMROOT)  &
         WRITE(6,'(1x,2(a,e13.6))')                                     &
          'prcr3min= ', amin,',  prcr3max=',amax

      CALL edgfill(prcrate(1,1,4),1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,1,1,1)
      CALL a3dmax0(prcrate(1,1,4),1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,1,1,1,amax,amin)
      IF (IAMROOT)  &
         WRITE(6,'(1x,2(a,e13.6))')                                     &
          'prcr4min= ', amin,',  prcr4max=',amax

      CALL edgfill(tke,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,nzd,1,nzd-1)
      CALL a3dmax0(tke,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,nzd,1,nzd-1,amax,amin)
      IF (IAMROOT)  &
         WRITE(6,'(1x,2(a,e13.4))')                                     &
          'tkemin  = ', amin,',  tkemax  =',amax

      CALL edgfill(kmh,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,nzd,1,nzd-1)
      CALL a3dmax0(kmh,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,nzd,1,nzd-1,amax,amin)
      IF (IAMROOT)  &
        WRITE(6,'(1x,2(a,e13.4))')                                          &
          'kmhmin  = ', amin,',  kmhmax  =',amax

      CALL edgfill(kmv,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,nzd,1,nzd-1)
      CALL a3dmax0(kmv,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,nzd,1,nzd-1,amax,amin)
      IF (IAMROOT)  &
        WRITE(6,'(1x,2(a,e13.4))')                                          &
          'kmvmin  = ', amin,',  kmvmax  =',amax

      DO is = 0, nstypsd

        IF(IAMROOT)  &
           WRITE(6,'(1x,a,I2.2)') ,'Max/min of soil model variables for type index ', is

        CALL edgfill(tsoil(1,1,1,is),1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,nzsoild,1,     &
              nzsoild)
        CALL a3dmax0(tsoil(1,1,1,is),1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,nzsoild,1,  &
              nzsoild,amax,amin)
        IF(IAMROOT)  &
           WRITE(6,'(1x,2(a,e13.6))')                                        &
            'tsoilmin = ', amin,',  tsoilmax =',amax

        CALL edgfill(qsoil(1,1,1,is),1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,nzsoild, &
               1,nzsoild)
        CALL a3dmax0(qsoil(1,1,1,is),1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,nzsoild, &
               1,nzsoild,amax,amin)
        IF(IAMROOT)  &
           WRITE(6,'(1x,2(a,e13.6))')                                        &
            'qsoilmin= ', amin,',  qsoilmax=',amax

        CALL edgfill(wetcanp(1,1,is),1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,1,1,1)
        CALL a3dmax0(wetcanp(1,1,is),1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,1,1,1,amax,amin)
        IF(IAMROOT)  &
           WRITE(6,'(1x,2(a,e13.6))')                                        &
            'wetcmin = ', amin,',  wetcmax =',amax

      END DO
    END IF

    RETURN
  END SUBROUTINE dtaread_get_3d

  !##################################################################

  SUBROUTINE dtaread_get_optional(radfrc,radsw,rnflx,radswnet,radlwin,  &
                                  usflx,vsflx,ptsflx,qvsflx,istatus)

  !
  !---------------------------------------------------------------------
  !
  !  PURPOSE:
  !
  !  Reads in flux and radiation variables from time dependent history
  !  file.
  !
  !---------------------------------------------------------------------
  !
  !  AUTHOR: Yunheng Wang
  !    06/02/2011.
  !
  !  MODIFICATION HISTORY:
  !
  !
  !---------------------------------------------------------------------
  !
  !    radfrc   Radiation forcing (K/s)
  !    radsw    Solar radiation reaching the surface
  !    rnflx    Net radiation flux absorbed by surface
  !    radswnet Net shortwave radiation, SWin - SWout
  !    radlwin  Incoming longwave radiation
  !
  !    usflx    Surface flux of u-momentum (kg/(m*s**2))
  !    vsflx    Surface flux of v-momentum (kg/(m*s**2))
  !    ptsflx   Surface heat flux (K*kg/(m**2 * s ))
  !    qvsflx   Surface moisture flux of (kg/(m**2 * s))
  !
  !---------------------------------------------------------------------

    USE arps_precision

    IMPLICIT NONE

    REAL(P), INTENT(OUT) :: radfrc  (nxd,nyd,nzd)    ! Radiation forcing (K/s)
    REAL(P), INTENT(OUT) :: radsw   (nxd,nyd)        ! Solar radiation reaching the surface
    REAL(P), INTENT(OUT) :: rnflx   (nxd,nyd)        ! Net radiation flux absorbed by surface
    REAL(P), INTENT(OUT) :: radswnet(nxd,nyd)        ! Net shortwave radiation
    REAL(P), INTENT(OUT) :: radlwin (nxd,nyd)        ! Incoming longwave radiation

    REAL(P), INTENT(OUT) :: usflx (nxd,nyd)          ! Surface flux of u-momentum (kg/(m*s**2))
    REAL(P), INTENT(OUT) :: vsflx (nxd,nyd)          ! Surface flux of v-momentum (kg/(m*s**2))
    REAL(P), INTENT(OUT) :: ptsflx(nxd,nyd)          ! Surface heat flux (K*kg/(m*s**2))
    REAL(P), INTENT(OUT) :: qvsflx(nxd,nyd)          ! Surface moisture flux (kg/(m**2*s))

    INTEGER :: istatus           ! Return status indicator

  !---------------------------------------------------------------------
  !
  ! Misc. Local variables
  !
  !---------------------------------------------------------------------

    REAL(P) :: amin, amax
  !
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !
  !  Beginning of executable code...
  !
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !
    istatus = 0

    IF (.NOT. hisopened) CALL open_files( 1,istatus )

    SELECT CASE (filefmt)
    CASE ( 1 )

      ! To be implemented

    CASE( 3 )  ! HDF 4 Format

      CALL hdf_read_optional(IAMROOT, IREAD, ISPLIT, ncompressx, ncompressy,&
                 fileid,nxd,nyd,nzd,nxin,nyin,                          &
                 radfrc,radsw,rnflx,radswnet,radlwin,                   &
                 usflx,vsflx,ptsflx,qvsflx,istatus)

    CASE (7)    ! NetCDF format

      ! To be implemented

    CASE DEFAULT
      WRITE(*,'(1x,a,I2)') 'ERROR: unsupported file format at present, filefmt = ',filefmt
      CALL arpsstop('ERROR: unsupported file format.',1)
    END SELECT

    IF(IAMROOT .AND. istatus == 0)  WRITE(6,'(/3a/)')                   &
           ' Radiation and Flux arrays from ',TRIM(hisfile),' were successfully read.'

    IF (dbglvl > 90) THEN
      IF (IAMROOT)  &
        WRITE(6,'(/1x,a/)') 'Min. and max. of the radiation and flux data arrays read in:'

      CALL edgfill(radfrc,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,nzd,1,nzd-1)
      CALL a3dmax0(radfrc,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,nzd,1,nzd-1,amax,amin)
      IF (IAMROOT)  WRITE(6,'(1x,2(a,e13.4))')                          &
          'radfmin = ', amin,',  radfmax =',amax

      CALL edgfill(radsw,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,1,1,1)
      CALL a3dmax0(radsw,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,1,1,1,amax,amin)
      IF (IAMROOT) WRITE(6,'(1x,2(a,e13.4))')                           &
          'radswmin= ', amin,',  radswmax=',amax

      CALL edgfill(rnflx,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,1,1,1)
      CALL a3dmax0(rnflx,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,1,1,1,amax,amin)
      IF (IAMROOT) WRITE(6,'(1x,2(a,e13.4))')                           &
          'rnflxmin= ', amin,',  rnflxmax=',amax

      CALL edgfill(radswnet,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,1,1,1)
      CALL a3dmax0(radswnet,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,1,1,1,amax,amin)
      IF (IAMROOT) WRITE(6,'(1x,2(a,e13.4))')                           &
          'radswnetmin= ', amin,',  radswnetmax=',amax

      CALL edgfill(radlwin,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,1,1,1)
      CALL a3dmax0(radlwin,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,1,1,1,amax,amin)
      IF (IAMROOT) WRITE(6,'(1x,2(a,e13.4))')                           &
          'radlwinmin= ', amin,',  radlwinmax=',amax

      CALL edgfill(usflx,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,1,1,1)
      CALL a3dmax0(usflx,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,1,1,1,amax,amin)
      IF (IAMROOT) WRITE(6,'(1x,2(a,e13.4))')                           &
          'usflxmin= ', amin,',  usflxmax=',amax

      CALL edgfill(vsflx,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,1,1,1)
      CALL a3dmax0(vsflx,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,1,1,1,amax,amin)
      IF (IAMROOT) WRITE(6,'(1x,2(a,e13.4))')                           &
          'vsflxmin= ', amin,',  vsflxmax=',amax

      CALL edgfill(ptsflx,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,1,1,1)
      CALL a3dmax0(ptsflx,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,1,1,1,amax,amin)
      IF (IAMROOT) WRITE(6,'(1x,2(a,e13.4))')                           &
          'ptflxmin= ', amin,',  ptflxmax=',amax

      CALL edgfill(qvsflx,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,1,1,1)
      CALL a3dmax0(qvsflx,1,nxd,1,nxd-1,1,nyd,1,nyd-1,1,1,1,1,amax,amin)
      IF (IAMROOT) WRITE(6,'(1x,2(a,e13.4))')                           &
          'qvflxmin= ', amin,',  qvflxmax=',amax

    END IF

    RETURN
  END SUBROUTINE dtaread_get_optional

  !##################################################################

  SUBROUTINE dtaread_get_dsd( ntcloud,n0rain, n0snow,  n0grpl,   n0hail,&
                                   rhoice,  rhosnow,  rhogrpl,  rhohail,&
                       alpharain,alphaice,alphasnow,alphagrpl,alphahail,&
                       istatus )
  !
  !-----------------------------------------------------------------------
  !
  !  PURPOSE:
  !
  !  Read DSD parameters from history file
  !
  !-----------------------------------------------------------------------
  !
  !  Author: Yunheng Wang
  !  Date: 06/03/11
  !
  !  MODIFICATION HISTORY:
  !
  !
  !-----------------------------------------------------------------------
  !
    IMPLICIT NONE
    REAL(P), INTENT(OUT) :: ntcloud
    REAL(P), INTENT(OUT) :: n0rain
    REAL(P), INTENT(OUT) :: n0snow
    REAL(P), INTENT(OUT) :: n0grpl
    REAL(P), INTENT(OUT) :: n0hail
    REAL(P), INTENT(OUT) :: rhoice
    REAL(P), INTENT(OUT) :: rhosnow
    REAL(P), INTENT(OUT) :: rhogrpl
    REAL(P), INTENT(OUT) :: rhohail
    REAL(P), INTENT(OUT) :: alpharain
    REAL(P), INTENT(OUT) :: alphaice
    REAL(P), INTENT(OUT) :: alphasnow
    REAL(P), INTENT(OUT) :: alphagrpl
    REAL(P), INTENT(OUT) :: alphahail
    INTEGER, INTENT(OUT) :: istatus
  !
  !-----------------------------------------------------------------------
  !
  !  Variable Declarations.
  !
  !-----------------------------------------------------------------------
  !


  !
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !
    istatus = 0

    IF (.NOT. hisopened) CALL open_files( 1,istatus )

    SELECT CASE (filefmt)
    CASE ( 1 )

      ! To be implemented

    CASE( 3 )  ! HDF 4 Format

      CALL hdf_read_dsd(IAMROOT, IREAD, ISPLIT, ncompressx, ncompressy, &
                      fileid,                                           &
                      ntcloud,n0rain, n0snow,  n0grpl,   n0hail,        &
                      rhoice,  rhosnow,  rhogrpl,  rhohail,             &
                      alpharain,alphaice,alphasnow,alphagrpl,alphahail, &
                      istatus)

    CASE (7)    ! NetCDF format

      ! To be implemented

    CASE DEFAULT
      WRITE(*,'(1x,a,I2)') 'ERROR: unsupported file format at present, filefmt = ',filefmt
      CALL arpsstop('ERROR: unsupported file format.',1)
    END SELECT

    RETURN
  END SUBROUTINE dtaread_get_dsd

  !##################################################################

  SUBROUTINE dtaread_get_head(month,day,year, hour,minute,second,       &
                        umove,vmove,xgrdorg,ygrdorg,                    &
                        mapproj,trulat1,trulat2,trulon, sclfct,         &
                        tstop,thisdmp,latitud,ctrlat,ctrlon,            &
                        grdin,basin, varin,mstin,icein,trbin,totin,     &
                        sfcin, rainin,tkein, prcin, radin,flxin,snowin, &
                        istatus )
  !
  !-----------------------------------------------------------------------
  !
  !  PURPOSE:
  !
  !  Read history file header.
  !
  !-----------------------------------------------------------------------
  !
  !  Author: Yunheng Wang
  !  Date: 06/03/11
  !
  !  MODIFICATION HISTORY:
  !
  !
  !-----------------------------------------------------------------------
  !
    IMPLICIT NONE

    INTEGER, INTENT(OUT) :: month,day,year
    INTEGER, INTENT(OUT) :: hour,minute,second

    REAL(P), INTENT(OUT) :: umove,vmove
    REAL(P), INTENT(OUT) :: xgrdorg,ygrdorg

    INTEGER, INTENT(OUT) :: mapproj
    REAL(P), INTENT(OUT) :: trulat1,trulat2,trulon
    REAL(P), INTENT(OUT) :: sclfct
    REAL(P), INTENT(OUT) :: tstop,thisdmp
    REAL(P), INTENT(OUT) :: latitud
    REAL(P), INTENT(OUT) :: ctrlat,ctrlon

    INTEGER, INTENT(OUT) :: grdin,basin, varin,mstin,icein,trbin,totin
    INTEGER, INTENT(OUT) :: sfcin, rainin,tkein, prcin, radin,flxin,snowin

    INTEGER, INTENT(OUT) :: istatus
  !
  !-----------------------------------------------------------------------
  !
  !  Variable Declarations.
  !
  !-----------------------------------------------------------------------
  !


  !
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !
    istatus = 0

    IF (.NOT. hisopened) CALL open_files( 1,istatus )

    SELECT CASE (filefmt)
    CASE ( 1 )

      ! To be implemented

    CASE( 3 )  ! HDF 4 Format

      CALL hdf_read_head(IAMROOT, IREAD, ISPLIT, ncompressx, ncompressy,&
                        fileid,                                         &
                        month,day,year, hour,minute,second,             &
                        umove,vmove,xgrdorg,ygrdorg,                    &
                        mapproj,trulat1,trulat2,trulon, sclfct,         &
                        tstop,thisdmp,latitud,ctrlat,ctrlon,            &
                        grdin,basin, varin,mstin,icein,trbin,totin,     &
                        sfcin, rainin,tkein, prcin, radin,flxin,snowin, &
                        istatus)

    CASE (7)    ! NetCDF format

      ! To be implemented

    CASE DEFAULT
      WRITE(*,'(1x,a,I2)') 'ERROR: unsupported file format at present, filefmt = ',filefmt
      CALL arpsstop('ERROR: unsupported file format.',1)
    END SELECT

    RETURN
  END SUBROUTINE dtaread_get_head

  !##################################################################

  SUBROUTINE dtaread_get_var(varname,varrank,vardim,varjoin,varout,istatus)

  !
  !---------------------------------------------------------------------
  !
  !  PURPOSE:
  !
  !  Reads in a variable given the variable name and rank.
  !
  !---------------------------------------------------------------------
  !
  !  AUTHOR: Yunheng Wang
  !    09/21/2012.
  !
  !  MODIFICATION HISTORY:
  !
  !---------------------------------------------------------------------

    USE arps_precision

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: varname
    INTEGER,          INTENT(IN) :: varrank
    INTEGER,          INTENT(IN) :: varjoin  ! 0 : no join
                                             ! 1 : Join in X direction
                                             ! 2 : Join in Y direction
                                             ! 3 : Join in Bother X and Y direction for the first 2 dimensions
    INTEGER,          INTENT(IN) :: vardim(varrank)

    REAL(P), INTENT(OUT) :: varout  (:,:,:,:)

    INTEGER, INTENT(OUT) :: istatus           ! Return status indicator

  !---------------------------------------------------------------------
  !
  ! Misc. Local variables
  !
  !---------------------------------------------------------------------

    INTEGER :: k, arrsize

    REAL(P) :: amin, amax
  !
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !
  !  Beginning of executable code...
  !
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !
    istatus = 0

    IF (.NOT. hisopened) CALL open_files( 1,istatus )

    arrsize = 1
    DO k = 1,varrank
      arrsize = arrsize*vardim(k)
    END DO

    SELECT CASE (filefmt)
    CASE ( 1 )

      ! To be implemented

    CASE( 3 )  ! HDF 4 Format

      CALL hdf_read_var(IAMROOT, IREAD, ISPLIT, ncompressx, ncompressy, &
                 fileid,nxd,nyd,nxin,nyin,                              &
                 varname,varrank,vardim,varjoin,varout,arrsize,istatus)

    CASE (7)    ! NetCDF format

      ! To be implemented

    CASE DEFAULT
      WRITE(*,'(1x,a,I2)') 'ERROR: unsupported file format at present, filefmt = ',filefmt
      CALL arpsstop('ERROR: unsupported file format.',1)
    END SELECT

    IF(IAMROOT .AND. istatus == 0)  WRITE(6,'(/4a/)')                   &
           ' Variable ',TRIM(varname),' read from ',TRIM(hisfile),'.'

    IF (dbglvl > 90) THEN
      IF (IAMROOT)  &
        WRITE(6,'(/1x,a/)') 'Min. and max. of the reading data arrays are:'

      SELECT CASE (varrank)
      CASE (1)
        CALL a3dmax0(varout,1,vardim(1),1,vardim(1),                    &
                            1,1,1,1, 1,1,1,1, amax,amin)
        IF (IAMROOT)  WRITE(6,'(1x,2(a,e13.4))')                        &
            'minmum = ', amin,',  maximum = ', amax
      CASE (2)
        CALL a3dmax0(varout,1,vardim(1),1,vardim(1),                    &
                            1,vardim(2),1,vardim(2),                    &
                            1,1,1,1, amax,amin)
        IF (IAMROOT)  WRITE(6,'(1x,2(a,e13.4))')                        &
            'minmum = ', amin,',  maximum = ', amax
      CASE (3)
        CALL a3dmax0(varout,1,vardim(1),1,vardim(1),                    &
                            1,vardim(2),1,vardim(2),                    &
                            1,vardim(3),1,vardim(3),                    &
                            amax,amin)
        IF (IAMROOT)  WRITE(6,'(1x,2(a,e13.4))')                        &
            'minmum = ', amin,',  maximum = ', amax
      CASE (4)
        DO k = 1,vardim(4)
          CALL a3dmax0(varout(:,:,:,k),1,vardim(1),1,vardim(1),         &
                                       1,vardim(2),1,vardim(2),         &
                                       1,vardim(3),1,vardim(3),         &
                                       amax,amin)
          IF (IAMROOT)  WRITE(6,'(1x,a,i2,2(a,e13.4))')                 &
              'Level ',k,' : minmum = ', amin,',  maximum = ', amax
        END DO
      END SELECT

    END IF

    RETURN
  END SUBROUTINE dtaread_get_var

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Private subroutines
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  SUBROUTINE open_files( fileflag,istatus )
  !
  !---------------------------------------------------------------------
  !
  !  PURPOSE:
  !
  !  Open file handlers
  !
  !---------------------------------------------------------------------
  !
  !  Author: Yunheng Wang
  !  Date: 06/01/11
  !
  !  MODIFICATION HISTORY:
  !
  !
  !---------------------------------------------------------------------
  !
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: fileflag
    INTEGER, INTENT(OUT) :: istatus
  !
  !---------------------------------------------------------------------
  !
  !  Variable Declarations.
  !
  !---------------------------------------------------------------------
  !
    INTEGER :: ii, jj
    CHARACTER(LEN=256) :: filebase, filename

  !
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !
    istatus = 0

    IF (fileflag > 0) THEN   ! History file
      filebase = hisfile
      IF (basopened ) CALL close_files ( istatus )
    ELSE
      filebase = basfile
      IF (hisopened) CALL close_files( istatus)
    END IF

    SELECT CASE (filefmt)
    CASE ( 3 )
      DO jj = 1, ncompressy
        DO ii = 1, ncompressx
          IF (IREAD) THEN
            IF (IONE) THEN
              filename = filebase
            ELSE
              CALL gtsplitfn(filebase,ncompressx,ncompressy,loc_x,loc_y,ii,jj, &
                     0,0,1,dbglvl,filename,istatus)
            END IF

            IF (dbglvl > 2) WRITE(*,'(1x,a,I4,2a)') 'Process: ',myproc,' opening HDF file: ', trim(filename)

            CALL hdfopen(filename,1,fileid(ii,jj))
            IF (fileid(ii,jj) < 0) THEN
              WRITE (*,'(1x,3a)') "module_dtaread: ERROR opening ",     &
                             trim(filename)," for reading."
              istatus = -1
              RETURN
            END IF
          END IF
        END DO
      END DO
    END SELECT

    IF (fileflag > 0) THEN
      hisopened = .TRUE.
    ELSE
      basopened = .TRUE.
    END IF

    RETURN
  END SUBROUTINE open_files

  !##################################################################

  SUBROUTINE close_files( istatus )
  !
  !---------------------------------------------------------------------
  !
  !  PURPOSE:
  !
  !  Close opened files
  !
  !---------------------------------------------------------------------
  !
  !  Author: Yunheng Wang
  !  Date: 06/01/11
  !
  !  MODIFICATION HISTORY:
  !
  !
  !---------------------------------------------------------------------
  !
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: istatus
  !
  !---------------------------------------------------------------------
  !
  !  Variable Declarations.
  !
  !---------------------------------------------------------------------
  !
    INTEGER :: ii, jj
  !
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !
    istatus = 0

    SELECT CASE (filefmt)
    CASE ( 3 )
      DO jj = 1, ncompressy
        DO ii = 1, ncompressx
          IF (IREAD) THEN
            CALL hdfclose(fileid(ii,jj),istatus)
            IF (istatus == 0 ) THEN
              IF (dbglvl > 0) WRITE(*,'(1x,a,I4,2a)') "Process: ",myproc," closed file ",TRIM(hisfile)
            ELSE
              WRITE(*,'(1x,a,I4,a,I3,2a)') "Process: ",myproc,          &
                " ERROR (status = ", istatus, ") on closing ",TRIM(hisfile)
            END IF
          END IF
        END DO
      END DO
    END SELECT

    hisopened = .FALSE.
    basopened = .FALSE.

    RETURN
  END SUBROUTINE close_files

END MODULE arps_dtaread
