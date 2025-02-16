!########################################################################
!########################################################################
!######                                                            ######
!######                    PROGRAM QPFMASK                         ######
!######                                                            ######
!######                      Developed by                          ######
!######         Center for Analysis and Prediction of Storms       ######
!######                   University of Oklahoma                   ######
!######                                                            ######
!########################################################################
!########################################################################

PROGRAM QPFMASK

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Reads in HDF file created by INTQPF, and creates a binary bitmap
! file listing which grid points have missing data.  The bitmap
! file is then read in by QPFSTATS and used to exclude certain
! grid points from the statistical calculations.
!
! AUTHOR:  Eric Kemp, March 2000.
!
! MODIFICATION HISTORY:
! Eric Kemp, 31 March 2000.
! Added lat/lon coordinates of four corners for NCL plotting.
!     
!-----------------------------------------------------------------------
!     
! Variable Declarations:
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: date(6)
  INTEGER,ALLOCATABLE :: timesec(:)  

  INTEGER :: mapproj
  REAL :: scale
  REAL :: trulat1, trulat2, trulon
  REAL :: dx,dy

  CHARACTER*8 :: model
  CHARACTER*4 :: grid
  CHARACTER*4,ALLOCATABLE,DIMENSION(:) :: id_p,id_sfc,idmask_sfc
  CHARACTER*32,ALLOCATABLE,DIMENSION(:) :: name_p,unit_p, &
                                           name_sfc,unit_sfc, &
                                           namemask_sfc,unitmask_sfc
  
  REAL,ALLOCATABLE, DIMENSION(:,:,:,:,:) :: var_p
  REAL,ALLOCATABLE, DIMENSION(:,:,:,:) :: var_sfc  
  REAL,ALLOCATABLE, DIMENSION(:) :: pressure
  REAL,ALLOCATABLE, DIMENSION(:,:) :: lat2d,lon2d
  REAL :: scswlat,scswlon
  REAL :: mapswlat,mapswlon,mapnwlat,mapnwlon, &
          mapselat,mapselon,mapnelat,mapnelon

  INTEGER :: npreslevel,nvar_p,nvar_sfc

  INTEGER :: flag_p,flag_sfc
  INTEGER :: status

  REAL,ALLOCATABLE :: mask(:,:,:,:)

  INTEGER :: i,j,k,l

  REAL,PARAMETER :: missing = -9999.

!-----------------------------------------------------------------------
!
! Namelists
!
!-----------------------------------------------------------------------

  INTEGER :: nx,ny,ntime
  CHARACTER*256 :: infilename
  NAMELIST /hdf_input/ nx,ny,ntime,infilename 

  CHARACTER*256 :: outfilename
  NAMELIST /output/ outfilename 

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  WRITE(6,'(//5x,a)')                                                   &
 &'###################################################################'
  WRITE(6,'(5x,a,/5x,a)')                                               &
 &'#                                                                 #',&
 &'# Welcome to QPFMASK, a program that reads in an interpolated     #'  
  WRITE(6,'(5x,a)')                                                     &
 &'# precipitation HDF file from INTQPF and creates a binary bitmap  #'
  WRITE(6,'(5x,a)')                                                     &
 &'# file for use by QPFSTATS.                                       #', &
 &'#                                                                 #',&
 &'###################################################################'
  WRITE(6,*)

!-----------------------------------------------------------------------
!
! Read in namelist
!
!-----------------------------------------------------------------------

  PRINT *, 'Reading NAMELIST hdf_input'
  READ(5,hdf_input)
  PRINT *, 'Reading NAMELIST output'
  READ(5,output)

!-----------------------------------------------------------------------
!
! Read precipitation data.
!
!-----------------------------------------------------------------------

  CALL rd_verif_dims(nx,ny,npreslevel,ntime,nvar_p,nvar_sfc, &
                     infilename)

  ALLOCATE (timesec(ntime), pressure(npreslevel),        &
            var_p(nx,ny,npreslevel,ntime,nvar_p), &
            var_sfc(nx,ny,ntime,nvar_sfc), &
            lat2d(nx,ny),lon2d(nx,ny), &
            id_p(nvar_p),name_p(nvar_p),unit_p(nvar_p), &
            id_sfc(nvar_sfc),name_sfc(nvar_sfc), &
            unit_sfc(nvar_sfc), &
            STAT=status)

  IF (status /= 0) CALL alloc_fail (status, 'f1')

  CALL rdverif ( nx,ny,npreslevel,ntime,nvar_p,nvar_sfc,         &
                 infilename,model,grid,date,timesec,pressure,        &
                 mapproj, scale, trulat1, trulat2,trulon, dx,dy,   &
                 scswlat,scswlon, &
                 mapswlat,mapswlon,mapnwlat,mapnwlon, &
                 mapselat,mapselon,mapnelat,mapnelon, &
                 flag_p,var_p, id_p, name_p, unit_p, &
                 flag_sfc,var_sfc, id_sfc, name_sfc, unit_sfc )

  IF (flag_sfc == 0) THEN
    WRITE(6,*)'ERROR:  Could not find surface data in HDF file.'
    WRITE(6,*)'Aborting...'
    STOP
  ENDIF

!-----------------------------------------------------------------------
!
! Make bit map
!
!-----------------------------------------------------------------------

  ALLOCATE(mask(nx,ny,ntime,nvar_sfc),idmask_sfc(nvar_sfc), &
           namemask_sfc(nvar_sfc),unitmask_sfc(nvar_sfc),&
           STAT=status)

  DO l = 1,nvar_sfc
    DO k = 1,ntime
      DO j = 1,ny
        DO i = 1,nx
          IF (var_sfc(i,j,k,l).ne.missing) THEN
            mask(i,j,k,l) = 1
          ELSE
            mask(i,j,k,l) = 0
          END IF
        END DO
      END DO
    END DO

    idmask_sfc(l) = 'MASK'
    namemask_sfc(l) = 'Mask for verification'
    unitmask_sfc(l) = 'none'
  END DO

!-----------------------------------------------------------------------
!
! Output bit map.
!
!-----------------------------------------------------------------------
  
  flag_p = 0   ! Don't output pressure level data
  flag_sfc = 1

  CALL wrtverif(nx,ny,npreslevel,ntime,nvar_p,nvar_sfc, missing,      &
                outfilename,model,grid,date,timesec,pressure,         &
                mapproj, scale, trulat1,trulat2,trulon,               &
                dx,dy,scswlat,scswlon,                                &
                mapswlat,mapswlon,mapnwlat,mapnwlon,                  &
                mapselat,mapselon,mapnelat,mapnelon,                  &
                flag_p,id_p, name_p, unit_p, var_p,                   &
                flag_sfc,idmask_sfc, namemask_sfc, unitmask_sfc,      &
                mask)

!-----------------------------------------------------------------------
!
! The end.
!
!-----------------------------------------------------------------------

  WRITE(6,*)'Program QPFMASK successfully completed.'

END PROGRAM qpfmask
