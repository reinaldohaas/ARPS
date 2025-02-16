!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE RADREAD                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE radread(filnam,lendtf,                                       &
           isrc,tim,ireturn,                                            &
           nx,ny,nz,dx,dy,dz,                                           &
           xor,yor,zor,                                                 &
           tem1,tem2,tem3,tem4,tem5,                                    &
           ubar,vbar,tem6,tem7,tem8)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Write history data into channel nchanl as unformatted binary data.
!
!  All output data are located at the grid cell centers.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Steven Lazarus
!  8/10/95.
!
!  MODIFICATION HISTORY:
!
!  02/22/96 (Limin Zhao)
!  Added ubar and vbar for mean velocity. Changes are also made to
!  incorporate the retrieval data.
!
!  08/08/96 (Limin Zhao)
!  Clean up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    filnam  Name of the input file
!    nx      Number of grid points in the x-direction (east/west)
!    ny      Number of grid points in the y-direction (north/south)
!    nz      Number of grid points in the vertical
!
!    dx      x grid spacing
!    dy      y grid spacing
!    dz      z grid spacing
!
!    xor     the location of radar station in x-direction
!    yor     the location of radar station in y-direction
!    zor     the location of radar station in z-direction
!
!    tim     time of retrieval data
!    isrc    Location of calling routine
!
!  OUTPUT:
!
!    tem1   Observed radial velocity
!    tem6   Retrieved u-component
!    tem7   Retrieved v-component
!    tem8   Retrieved w-component
!    tem5   Observed reflectivity
!    ubar   Retrieved mean velocity of u-component
!    vbar   Retrieved mean velocity of v-component
!
!  WORK ARRAYS:
!
!    tem2      Temporary work array.
!    tem3      Temporary work array.
!    tem4      Temporary work array.

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

  INTEGER :: nx,ny,nz       ! Number of model grid points in 3 directions
  INTEGER :: nxr,nyr,nzr    ! Number of real data in 3 directions

  INTEGER :: lendtf,i,j,k,ireturn

  CHARACTER (LEN=256) :: filnam
  CHARACTER (LEN=10)  :: radarsp
  INTEGER :: iyrrr,imorr,idarr,iharr,imarr,isarr

  REAL :: dxr,dyr,dzr
  REAL :: dx,dy,dz
  REAL :: xor, yor, zor
  REAL :: bad
  REAL :: tim
  REAL :: umean, vmean

  REAL :: tem1  (nx,ny,nz)  ! Temporary work array
  REAL :: tem2  (nx,ny,nz)  ! Temporary work array
  REAL :: tem3  (nx,ny,nz)  ! Temporary work array
  REAL :: tem4  (nx,ny,nz)  ! Temporary work array
  REAL :: tem5  (nx,ny,nz)  ! Temporary work array
  REAL :: tem6  (nx,ny,nz)  ! Temporary work array
  REAL :: tem7  (nx,ny,nz)  ! Temporary work array
  REAL :: tem8  (nx,ny,nz)  ! Temporary work array
  REAL :: ubar  (nx,ny,nz)  ! Temporary work array
  REAL :: vbar  (nx,ny,nz)  ! Temporary work array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  REAL :: tref

  INTEGER :: l
  INTEGER :: nchanl            ! FORTRAN I/O channel number for output
  INTEGER :: istat
  INTEGER :: ierr
  INTEGER :: isrc

  LOGICAL :: fexist
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'assim.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
!-----------------------------------------------------------------------
!
!  Assign a proper filename, and check to see if the file exist.
!
!  Note: Cray routines are used to force binary data file to be
!        in the IEEE format or compatable with IBM.
!
!-----------------------------------------------------------------------
!
  CALL asnctl ('NEWLOCAL', 1, ierr)
  CALL asnfile(filnam(1:lendtf), '-F f77 -N ieee', ierr)

  CALL getunit(nchanl)

  WRITE(6,*) 'nchanl,filnam: ',nchanl,filnam(1:lendtf)

  OPEN(UNIT=nchanl,FILE=filnam(1:lendtf),                               &
       STATUS='unknown', FORM='unformatted',IOSTAT= istat )

  IF( istat /= 0 ) GO TO 998

  100   CONTINUE

!
!-----------------------------------------------------------------------
!
!    Read data in.
!
!-----------------------------------------------------------------------
!
!
  READ(nchanl,ERR=110,END=120) radarsp,iyrrr,imorr,idarr,               &
                                       iharr,imarr,isarr
  READ(nchanl,ERR=110,END=120) tim

  WRITE(6,'(3x,a,f10.2,a)') 'Data read at ',tim,'seconds'

  IF(isrc == 1.OR.isrc == 2) THEN
    CALL retunit( nchanl )
    CLOSE(nchanl)
    ireturn=0
    RETURN
  END IF

  READ(nchanl,ERR=110,END=120) xor, yor, zor

  READ(nchanl,ERR=110,END=120) nxr,nyr,nzr,dxr,dyr,dzr

!
!-----------------------------------------------------------------------
!
!  Check to make sure the dimensions of the input data match that of
!  the current model stream.
!
!-----------------------------------------------------------------------
!
  IF((nx /= nxr).OR.(ny /= nyr).OR.(nz /= nzr)) THEN

    WRITE(6,'(a,/a,i5,a,i5,a,i5,/a,i5,a,i5,a,i5)')                      &
        ' Array dimension(s) of the input file inconsistent with ',     &
        ' model definitions, dimensions in input data were nx=',nxr,    &
        ', ny=',nyr,', nz=',nzr,' the model definitions were nx=',      &
        nx,' ny= ', ny, ' nz= ',nz
    WRITE(6,'(a)') ' Job stopped in subroutine RADREAD.'
    STOP

  END IF

  IF(ABS(dx -  dxr) > 0.1.OR.ABS(dy - dyr) > 0.1 .OR.ABS(dz -  dzr) > 0.1) THEN

    WRITE(6,'(a,/a,f10.2,a,f10.2,a,f10.2,/a,f10.2,a,f10.2,a,f10.2)')    &
        'Grid interval in the input data inconsisent with',             &
        'model definitions. In the input data dx=',dxr,                 &
        ', dy=',dyr,', dz=',dzr,' the model definitions were dx=',      &
        dx,' dy= ', dy, ' dz= ',dz
    WRITE(6,'(a)') ' Job stopped in subroutine RADREAD.'
    STOP

  END IF
!
!-----------------------------------------------------------------------
!
!  Continue reading
!
!-----------------------------------------------------------------------
!
  READ(nchanl,ERR=110,END=120) tem1      ! input Vr

  READ(nchanl,ERR=110,END=120) tem2      ! input retrieved u component

  READ(nchanl,ERR=110,END=120) tem3      ! input retrieved v component

  READ(nchanl,ERR=110,END=120) tem4      ! input retrieved w component

  READ(nchanl,ERR=110,END=120) tem5      ! input Zr

  CALL retunit(nchanl)

  CLOSE(nchanl)

  WRITE (*,*) "XXX RADREAD have in nx,ny,nz,dx,dy,dz",                  &
              nx,ny,nz,dx,dy,dz
  WRITE (*,*) "XXX RADREAD read in nx,ny,nz,dx,dy,dz,xor,yor,zor",      &
              nxr,nyr,nzr,dxr,dyr,dzr

!   write(31) tem2
!   write(31) tem3
!   write(31) tem4
!
!-----------------------------------------------------------------------
!
!  Calculate the mean profiles.
!  Note: the retrieval winds are total velocity. No overwritten of
!        background wind for current version.
!
!-----------------------------------------------------------------------
!
!       DO 160 k = 1,nz-1
!       DO 160 i = 1,nx
!       DO 160 j = 1,ny-1

  DO k = 1,nz
    DO i = 1,nx
      DO j = 1,ny
        ubar(i,j,k)   =  0.0
        tem6(i,j,k)  = tem2(i,j,k) - ubar(i,j,k)
      END DO
    END DO
  END DO

!       DO 170 k = 1,nz-1
!       DO 170 i = 1,nx-1
!       DO 170 j = 1,ny

  DO k = 1,nz
    DO i = 1,nx
      DO j = 1,ny
        vbar(i,j,k)  = 0.0
        tem7(i,j,k) = tem3(i,j,k) - vbar(i,j,k)
      END DO
    END DO
  END DO


!       DO 180 k = 1,nz-1
!       DO 180 i = 1,nx-1
!       DO 180 j = 1,ny

  DO k = 1,nz
    DO i = 1,nx
      DO j = 1,ny
        tem8(i,j,k) = tem4(i,j,k)
      END DO
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Friendly exit message
!
!----------------------------------------------------------------------
!
  930   CONTINUE

  WRITE(6,'(/a/)') 'Reading was successfully in RADREAD'
  ireturn = 0
  RETURN
!
!-----------------------------------------------------------------------
!
!  Error during read
!
!----------------------------------------------------------------------
!
  110   CONTINUE
  WRITE(6,'(/a/)') ' Error reading data in RADREAD'
  ireturn=1
  RETURN
!
!-----------------------------------------------------------------------
!
!  End-of-file during read
!
!----------------------------------------------------------------------
!

  120   CONTINUE
  WRITE(6,'(/a/)') ' End of file reached in RADREAD'
  ireturn=2
  RETURN
  998   CONTINUE

  WRITE(6,'(/1x,a,/1x,a/)')                                             &
       'File '//filnam(1:lendtf)                                        &
       //' not found.',                                                 &
       'Program returned from RADREAD.'
  ireturn=1

  RETURN
END SUBROUTINE radread
