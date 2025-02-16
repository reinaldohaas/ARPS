!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE IMG3D                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE img3d(b,x,y,z,dx,dy,                                         &
           nx,ibgn,iend,isk, ny,jbgn,jend,jsk, nz,kbgn,kend,ksk,        &
           bmax,bmin, label,fnkey,time,mode,kslice,jslice,islice,       &
           n,xp,yp,b1,b2,zs2, mgrid, nestgrd, tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Produce 2-D HDF images for specified slices of 3-d array b.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!
!  MODIFICATION HISTORY:
!
!  9/1/94 (Y. Lu)
!  Cleaned up documentation
!
!  2/1/1995 (J. Zhang)
!    Added constant p-level plots (mode.eq.6)
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    b        3-dimensional array of data
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!
!    nx       first dimension of b
!    ibgn     index of first i grid point to be used.
!    iend     index of last  i grid point to be used.
!    isk      skip or interpolation parameter in i-dir.
!             If isk >0, every isk'th data is used,
!             If isk <0, each grid zone is sub-divided into isk zones.
!
!    ny       second dimension of b
!    jbgn     index of first j grid point to be used.
!    jend     index of last  j grid point to be used.
!    jsk      skip or interpolation parameter in j-dir.
!             If jsk >0, every jsk'th data is used,
!             If jsk <0, each grid zone is sub-divided into jsk zones.
!
!    nz       third dimension of b
!    kbgn     index of first k grid point to be used.
!    kend     index of last  k grid point to be used.
!    ksk      skip or interpolation parameter in k-dir.
!             If ksk >0, every ksk'th data is used,
!             If ksk <0, each grid zone is sub-divided into ksk zones.
!
!    label    character string describing the contents of b
!    fnkey    file name key
!    time     time of data in seconds
!
!    mode     slice orientation indicator
!             mode = 1, x-y slice at z index kslice is plotted.
!             mode = 2, x-z slice at y index jslice is plotted.
!             mode = 3, y-z slice at x index islice is plotted.
!             mode = 0, all of the three slices above are plotted.
!
!    kslice   k index of plane for mode=1 x-y slice
!    jslice   j index of plane for mode=2 x-z slice
!    islice   i index of plane for mode=1 y-z slice
!
!  OUTPUT :
!
!    n        horizontal dimension of arbitary vertical cross section
!    xp       x-coordinate of grid points on arbitary
!             vertical cross-section
!    yp       y-coordinate of grid points on arbitary
!             vertical cross-section
!    b1       2-D field interpolated to a horizontal cross-section
!    b2       2-D field interpolated to a vertical corss-section
!    zs2      z-coordinate of grid points on arbitary vertical
!             cross-section
!    mgrid    the grid number
!    nestgrd  flag for nested grid run.
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz       ! dimemsions of 3-D input data arrray b
  INTEGER :: n              ! horizontal dimension of arbitary
                            ! vertical cross section.

  REAL :: b(nx,ny,nz)       ! 3-dimensional array of data
  REAL :: x(nx,ny,nz)       ! x coordinate of grid points in physical
                            ! /comp. space (m)
  REAL :: y(ny,ny,nz)       ! y coordinate of grid points in physical
                            ! /comp. space (m)
  REAL :: z(nz,ny,nz)       ! z coordinate of grid points in physical
                            ! /comp. space (m)
  REAL :: dx, dy

  REAL :: bmax,bmin         ! upper and lower bounds of b used in IMG2D

  INTEGER :: kslice         ! k index of plane for mode=1 x-y slice
  INTEGER :: jslice         ! j index of plane for mode=1 x-z slice
  INTEGER :: islice         ! i index of plane for mode=1 y-z slice

  REAL :: b1(nx,ny)         ! 2-D field interpolated to a horizontal
                            ! cross-section
  REAL :: b2(n,nz)          ! 2-D field interpolated to a vertical
                            ! corss-section.
  REAL :: zs2(n,nz)         ! z-coordinate of grid points on arbitary
                            ! vertical cross-section
  REAL :: xp(n)             ! x-coordinate of grid points on arbitary
                            ! vertical cross-section
  REAL :: yp(n)             ! y-coordinate of grid points on arbitary
                            ! vertical cross-section

  INTEGER :: ibgn,iend,isk  ! index of first and last i grid point;
                            ! skip/interpolation parameter in i-dir
  INTEGER :: jbgn,jend,jsk  ! index of first and last j grid point;
                            ! skip/interpolation parameter in j-dir
  INTEGER :: kbgn,kend,ksk  ! index of first and last k grid point;
                            ! skip/interpolation parameter in k-dir
  INTEGER :: mode           ! slice orientation indicator
  INTEGER :: mgrid          ! the grid number
  INTEGER :: nestgrd        ! flag for nested grid run.
  REAL :: tem1(*)           ! temporary work array

  INTEGER :: length         ! the length of string

  CHARACTER (LEN=*   ) :: label ! character string describing the contents
                                ! of b
  CHARACTER (LEN=80  ) :: timsnd  ! character string describing the time of b
  CHARACTER (LEN=*   ) :: fnkey ! file name key
  INTEGER :: tmstrln
  REAL :: time              ! time of data (s)
!
!-----------------------------------------------------------------------
!
!  Common blocks for plotting control parameters
!
!-----------------------------------------------------------------------
!
  REAL :: x01,y01            ! the first  point of interpolation
  REAL :: x02,y02            ! the second point of interpolation
  REAL :: z01                ! the given height of the slice
  REAL :: sinaf,cosaf,dist,sqrtdxy
  COMMON /slicev/x01,y01,x02,y02,sinaf,cosaf,dist,sqrtdxy
  COMMON /sliceh/z01
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,ij,ik,jk,idata,jdata,kdata
  CHARACTER (LEN=256) :: imgfn
  CHARACTER (LEN=35)  :: gridnum
  DATA gridnum /'123456789abcdefghijklmnopqrstuvwxyz'/
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
!  mode=1   Plot x-y cross-section
!
!-----------------------------------------------------------------------
!
  CALL cvttsnd( time, timsnd, tmstrln )

  IF(mode == 1 .OR. mode == 0 ) THEN

    k = kslice
    DO j=jbgn,jend
      DO i=ibgn,iend
        ij = i-ibgn+1 + (j-jbgn)*(iend-ibgn+1)
        tem1(ij)=b(i,j,k)
      END DO
    END DO

    imgfn = fnkey
    length = LEN( fnkey )
    WRITE(imgfn(length+1:length+5),'(a,i3.3)') '.k',k
    imgfn(length+6:256) = label//timsnd(1:tmstrln)
    length = 256
    CALL strlnth( imgfn, length)

    IF( nestgrd == 1 ) THEN
      WRITE(imgfn(length+1:length+3),'(''.G'',A)')                      &
            gridnum(mgrid:mgrid)
      length = length + 3
    END IF

    idata = iend-ibgn+1
    jdata = jend-jbgn+1

    CALL img2d(tem1,idata,1,idata,isk,jdata,1,jdata,isk,                &
               bmax,bmin,                                               &
               imgfn(1:length))

!
!-----------------------------------------------------------------------
!
!  mode=2   Plot x-z cross-section
!
!-----------------------------------------------------------------------
!
  ELSE IF (mode == 2 .OR. mode == 0 ) THEN

    j = jslice
    DO k=kbgn,kend
      DO i=ibgn,iend
        ik = i-ibgn+1 + (k-kbgn)*(iend-ibgn+1)
        tem1(ik)=b(i,j,k)
      END DO
    END DO

    imgfn = fnkey
    length = LEN( fnkey )
    WRITE(imgfn(length+1:length+5),'(a,i3.3)') '.j',j
    imgfn(length+6:256) = label//timsnd(1:tmstrln)
    length = 256
    CALL strlnth( imgfn, length)

    IF( nestgrd == 1 ) THEN
      WRITE(imgfn(length+1:length+3),'(''.G'',A)')                      &
            gridnum(mgrid:mgrid)
      length = length + 3
    END IF

    idata = iend-ibgn+1
    kdata = kend-kbgn+1

    CALL img2d(tem1,idata,1,idata,isk,kdata,1,kdata,ksk,                &
               bmax,bmin,                                               &
               imgfn(1:length))

!
!-----------------------------------------------------------------------
!
!  mode=3   Plot y-z cross-section
!
!-----------------------------------------------------------------------
!
  ELSE IF ( mode == 3 .OR. mode == 0) THEN

    i = islice
    DO k=kbgn,kend
      DO j=jbgn,jend
        jk = j-jbgn+1 + (k-kbgn)*(jend-jbgn+1)
        tem1(jk)=b(i,j,k)
      END DO
    END DO

    imgfn = fnkey
    length = LEN( fnkey )
    WRITE(imgfn(length+1:length+5),'(a,i3.3)') '.i',i
    imgfn(length+6:256) = label//timsnd(1:tmstrln)
    length = 256
    CALL strlnth( imgfn, length)

    IF( nestgrd == 1 ) THEN
      WRITE(imgfn(length+1:length+3),'(''.G'',A)')                      &
            gridnum(mgrid:mgrid)

      length = length + 3
    END IF

    jdata = jend-jbgn+1
    kdata = kend-kbgn+1

    CALL img2d(tem1,jdata,1,jdata,jsk,kdata,1,kdata,ksk,                &
               bmax,bmin,                                               &
               imgfn(1:length))
!
!-----------------------------------------------------------------------
!
!  mode=4   Plot horizontal slice at given height
!
!-----------------------------------------------------------------------
!
  ELSE IF( mode == 4 ) THEN

    CALL secthrz(nx,ny,nz,b,z,b1)

    DO j=jbgn,jend
      DO i=ibgn,iend
        ij = i-ibgn+1 + (j-jbgn)*(iend-ibgn+1)
        tem1(ij)=b1(i,j)
      END DO
    END DO

    imgfn = fnkey
    length = LEN( fnkey )
    WRITE(imgfn(length+1:length+5),'(a,i3.3)') '.z',nint(z01)
    imgfn(length+6:256) = label//timsnd(1:tmstrln)
    length = 256
    CALL strlnth( imgfn, length)

    IF( nestgrd == 1 ) THEN
      WRITE(imgfn(length+1:length+3),'(''.G'',A)')                      &
            gridnum(mgrid:mgrid)

      length = length + 3
    END IF

    idata = iend-ibgn+1
    jdata = jend-jbgn+1

    CALL img2d(tem1,idata,1,idata,isk,jdata,1,jdata,jsk,                &
               bmax,bmin,                                               &
               imgfn(1:length))
!
!-----------------------------------------------------------------------
!
!  mode=5   Plot vectical slice through two given points
!
!-----------------------------------------------------------------------
!
  ELSE IF( mode == 5 ) THEN

    CALL sectvrt(nx,ny,nz,b,x,y,z,dx,dy,b2,zs2,n,xp,yp)

    DO k=kbgn,kend
      DO i=ibgn,iend
        ik = i-ibgn+1 + (k-kbgn)*(iend-ibgn+1)
        tem1(ik)=b2(i,k)
      END DO
    END DO

    imgfn = fnkey
    length = LEN( fnkey )
    WRITE(imgfn(length+1:length+5),'(a,i3.3)') '.v',nint(x01)
    imgfn(length+6:256) = label//timsnd(1:tmstrln)
    length = 256
    CALL strlnth( imgfn, length)

    IF( nestgrd == 1 ) THEN
      WRITE(imgfn(length+1:length+3),'(''.G'',A)')                      &
            gridnum(mgrid:mgrid)

      length = length + 3
    END IF

    idata = iend-ibgn+1
    kdata = kend-kbgn+1

    CALL img2d(tem1,idata,1,idata,isk,kdata,1,kdata,ksk,                &
               bmax,bmin,                                               &
               imgfn(1:length+4))

!
!-----------------------------------------------------------------------
!
!  mode=6   Plot on constant pressure-level
!
!-----------------------------------------------------------------------
!
  ELSE IF( mode == 6 ) THEN

    CALL secthrz(nx,ny,nz,b,z,b1)

    DO j=jbgn,jend
      DO i=ibgn,iend
        ij = i-ibgn+1 + (j-jbgn)*(iend-ibgn+1)
        tem1(ij)=b1(i,j)
      END DO
    END DO

    imgfn = fnkey
    length = LEN( fnkey )
    WRITE(imgfn(length+1:length+5),'(a,i3.3)') '.z',nint(z01)
    imgfn(length+6:256) = label//timsnd(1:tmstrln)
    length = 256
    CALL strlnth( imgfn, length)

    IF( nestgrd == 1 ) THEN
      WRITE(imgfn(length+1:length+3),'(''.G'',A)')                      &
            gridnum(mgrid:mgrid)

      length = length + 3
    END IF

    idata = iend-ibgn+1
    jdata = jend-jbgn+1

    CALL img2d(tem1,idata,1,idata,isk,jdata,1,jdata,jsk,                &
               bmax,bmin,                                               &
               imgfn(1:length))
!

  END IF

  RETURN
END SUBROUTINE img3d
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE IMG3D0                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE img3d0(b,nx,ibgn,iend,isk, ny,jbgn,jend,jsk,                 &
           nz,kbgn,kend,ksk,                                            &
           bmax,bmin, label,fnkey,time,mode,kslice,jslice,islice,       &
           mgrid,nestgrd, tem1)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Produce 2-D HDF images for specified slices of 3-d array b.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!
!  MODIFICATION HISTORY:
!
!  9/1/94 (Y. Lu)
!  Cleaned up documentation
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    b        3-dimensional array of data
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!
!    nx       first dimension of b
!    ibgn     index of first i grid point to plot
!    iend     index of last  i grid point to plot
!             if isk >0, every isk'th data is used,
!             if isk <0, each grid zone is sub-divided into isk zone.
!
!    ny       second dimension of b
!    jbgn     index of first j grid point to plot
!    jend     index of last  j grid point to plot
!             if jsk >0, every jsk'th data is used,
!             if jsk <0, each grid zone is sub-divided into jsk zone.
!
!    nz       third dimension of b
!    kbgn     index  f first k grid point to plot
!    kend     index of last  k grid point to plot
!             if ksk >0, every ksk'th data is used,
!             if ksk <0, each grid zone is sub-divided into ksk zone.
!
!    label    character string describing the contents of b
!    fnkey    File name key
!    time     time of data in seconds
!
!    mode     slice orientation indicator
!             mode = 1, x-y slice at z index kslice is plotted.
!             mode = 2, x-z slice at y index jslice is plotted.
!             mode = 3, y-z slice at x index islice is plotted.
!             mode = 0, all of the three slices above are plotted.
!
!    kslice   k index of plane for mode=1 x-y slice
!    jslice   j index of plane for mode=2 x-z slice
!    islice   i index of plane for mode=1 y-z slice
!
!    mgrid    The grid number
!    nestgrd  Flag for nested grid run.
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz       ! 3 dimemsions of b
!  integer n              ! dimension on horizontal dirction

  REAL :: b(nx,ny,nz)       ! 3-dimensional array of data

  REAL :: bmax,bmin         ! the maximum and minimum value of b

  INTEGER :: ibgn,iend,isk  ! index of first and last i grid point;
                            ! skip/interpolation parameter in i-dir
  INTEGER :: jbgn,jend,jsk  ! index of first and last j grid point;
                            ! skip/interpolation parameter in j-dir
  INTEGER :: kbgn,kend,ksk  ! index of first and last k grid point;
                            ! skip/interpolation parameter in k-dir
  INTEGER :: length         ! the length of string

  INTEGER :: idata          ! total number of grid points in i-dir
  INTEGER :: jdata          ! total number of grid points in j-dir
  INTEGER :: kdata          ! total number of grid points in k-dir

  CHARACTER (LEN=*   ) :: label ! character string describing the contents
                                ! of b
  CHARACTER (LEN=80  ) :: timsnd  ! character string describing thetime of b
  CHARACTER (LEN=*   ) :: fnkey ! file name key
  INTEGER :: tmstrln        ! slice orientation indicator
  REAL :: time              ! time of data (s)

  INTEGER :: mode           ! slice orientation indicator
  INTEGER :: kslice         ! k index of plane for mode=1 x-y slice
  INTEGER :: jslice         ! j index of plane for mode=1 x-z slice
  INTEGER :: islice         ! i index of plane for mode=1 y-z slice

  INTEGER :: mgrid          ! the grid number
  INTEGER :: nestgrd        ! flag for nested grid run.

  REAL :: tem1(*)           ! temporary work array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,ij,ik,jk
  CHARACTER (LEN=256) :: imgfn
  CHARACTER (LEN=35)  :: gridnum
  DATA gridnum /'123456789abcdefghijklmnopqrstuvwxyz'/
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
!  mode=1   Plot x-y cross-section
!
!-----------------------------------------------------------------------
!
  CALL cvttsnd( time, timsnd, tmstrln )

  IF(mode == 1 .OR. mode == 0 ) THEN

    k = kslice
    DO j=jbgn,jend
      DO i=ibgn,iend
        ij = i-ibgn+1 + (j-jbgn)*(iend-ibgn+1)
        tem1(ij)=b(i,j,k)
      END DO
    END DO

    imgfn = fnkey
    length = LEN( fnkey )
    WRITE(imgfn(length+1:length+5),'(a,i3.3)') '.k',k
    imgfn(length+6:256) = label//timsnd(1:tmstrln)
    length = 256
    CALL strlnth( imgfn, length)

    IF( nestgrd == 1 ) THEN
      WRITE(imgfn(length+1:length+3),'(''.G'',A)')                      &
            gridnum(mgrid:mgrid)

      length = length + 3
    END IF

    idata = iend-ibgn+1
    jdata = jend-jbgn+1

    CALL img2d(tem1,idata,1,idata,isk,jdata,1,jdata,isk ,bmax,bmin,     &
               imgfn(1:length+4))
!
!-----------------------------------------------------------------------
!
!  mode=2   Plot x-z cross-section
!
!-----------------------------------------------------------------------
!
  ELSE IF (mode == 2 .OR. mode == 0 ) THEN

    j = jslice
    DO k=kbgn,kend
      DO i=ibgn,iend
        ik = i-ibgn+1 + (k-kbgn)*(iend-ibgn+1)
        tem1(ik)=b(i,j,k)
      END DO
    END DO

    imgfn = fnkey
    length = LEN( fnkey )
    WRITE(imgfn(length+1:length+5),'(a,i3.3)') '.j',j
    imgfn(length+6:256) = label//timsnd(1:tmstrln)
    length = 256
    CALL strlnth( imgfn, length)

    IF( nestgrd == 1 ) THEN
      WRITE(imgfn(length+1:length+3),'(''.G'',A)')                      &
            gridnum(mgrid:mgrid)

      length = length + 3
    END IF


    idata = iend-ibgn+1
    kdata = kend-kbgn+1

    CALL img2d(tem1,idata,1,idata,isk,kdata,1,kdata,ksk,                &
               bmax,bmin,                                               &
               imgfn(1:length+4))
!
!-----------------------------------------------------------------------
!
!  mode=3   Plot y-z cross-section
!
!-----------------------------------------------------------------------
!
  ELSE IF ( mode == 3 .OR. mode == 0) THEN

    i = islice
    DO k=kbgn,kend
      DO j=jbgn,jend
        jk = j-jbgn+1 + (k-kbgn)*(jend-jbgn+1)
        tem1(jk)=b(i,j,k)
      END DO
    END DO

    imgfn = fnkey
    length = LEN( fnkey )
    WRITE(imgfn(length+1:length+5),'(a,i3.3)') '.i',i
    imgfn(length+6:256) = label//timsnd(1:tmstrln)
    length = 256
    CALL strlnth( imgfn, length)

    IF( nestgrd == 1 ) THEN
      WRITE(imgfn(length+1:length+3),'(''.G'',A)')                      &
            gridnum(mgrid:mgrid)

      length = length + 3
    END IF

    jdata = jend-jbgn+1
    kdata = kend-kbgn+1

    CALL img2d(tem1,jdata,1,jdata,jsk,kdata,1,kdata,ksk,                &
               bmax,bmin,                                               &
               imgfn(1:length+4))

  ELSE IF( mode == 4 .OR. mode == 5 .OR. mode == 6 ) THEN

    WRITE(6,'(a)')                                                      &
        'Mode 4 or 5 is not available with IMG3D0, call IGM3D instead.'

  END IF

  RETURN
END SUBROUTINE img3d0
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE IMG2D                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE img2d(z,nx,ibgn,iend,isk,ny,jbgn,jend,jsk,                   &
           zmax,zmin,imgfn)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Create an HDF image file 'imgfn' for a 2-D field stored in z.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  4/14/93
!
!  MODIFICATION HISTORY:
!
!  9/1/94 (Y. Lu)
!  Cleaned up documentation
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    z        a two array containing the field for which a image file
!             is to be created.
!    nx, ny   the dimensions of array z.
!    ibgn     index of first i grid point to be used.
!    iend     index of last  i grid point to be used.
!    isk      skip or interpolation parameter in i-dir.
!             If isk >0, every isk'th data is used,
!             If isk <0, each grid zone is sub-divided into isk zones.
!
!    jbgn     index of first j grid point to be used.
!    jend     index of last  j grid point to be used.
!    jsk      skip or interpolation parameter in j-dir.
!             If jsk >0, every jsk'th data is used,
!             If jsk <0, each grid zone is sub-divided into jsk zones.
!
!    zmax,zmin
!             prespecified upper and lower bounds of the array value
!             to be used to scale the image.
!
!    imgfn    the name of the output image file.
!
!  OUTPUT:
!
!    None.
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

  INTEGER :: nx, ny         ! dimensions of array z
  REAL :: z(nx,ny)       ! A two array containing the field for
                         ! which a image file is to be created.

  INTEGER :: ibgn,iend,isk  ! index of first and last i grid point;
                            ! skip/interpolation parameter in i-dir
  INTEGER :: jbgn,jend,jsk  ! index of first and last j grid point;
                            ! skip/interpolation parameter in j-dir

  REAL :: zmax,zmin      ! prespecified upper and lower bounds
                         ! of the array value to be used to scale
                         ! the image.

  CHARACTER (LEN=*) :: imgfn    ! name of the output image file


  INTEGER :: maxdim         ! maximum number of gird points in work
                            ! array image
  INTEGER :: icompres
  INTEGER :: idata          ! total number of grid points in i-dir
  INTEGER :: jdata          ! total number of grid points in j-dir

  INTEGER :: d8pimg, iret, imax,imin, nimage

  PARAMETER (maxdim= 512 , icompres=11)
  CHARACTER (LEN=1) :: image(maxdim*maxdim)

  REAL :: zx, zn
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  IF( isk > 0 )THEN
    idata = (iend-ibgn)/isk + 1
  ELSE
    idata = (iend-ibgn)*ABS(isk) + 1
  END IF

  IF( jsk > 0 )THEN
    jdata = (jend-jbgn)/jsk + 1
  ELSE
    jdata = (jend-jbgn)*ABS(jsk) + 1
  END IF


  IF( idata*jdata > maxdim*maxdim ) THEN
    WRITE(6,'(1x,a,a/)') 'Character work array image ',                 &
        'defined in IMG2D not big enough. Job stopped in IMG2D.'
    CALL arpsstop('arpsstop called from IMG2D due to character length',1)
  END IF
!
!-----------------------------------------------------------------------
!
!  RASTERIZE and store as HDF image
!
!-----------------------------------------------------------------------
!
  CALL rasteriz(z, nx,ibgn,iend,isk, ny,jbgn,jend,jsk,zmax,zmin,        &
                image,nimage,imax,imin)

!mx
  print*,'D8pimg is not called. No HDF image produced.'

  RETURN  ! D8pimg not found in HDF library for some reason

! iret = D8pimg (imgfn , image, idata,jdata, icompres)
!
!  IF (iret /= 0) THEN
!    WRITE(6,'(/1x,a,a,a,i3)')                                           &
!        'Error writing HDF file ',imgfn,'. Error flag was ',iret
!  END IF

!
!-----------------------------------------------------------------------
!
!  Find Max & Min (for output to screen only)
!
!-----------------------------------------------------------------------
!
  CALL a3dmax0(z, 1,nx,ibgn,iend, 1,ny,jbgn,jend,1,1,1,1, zx, zn)

  WRITE(6,'(/1x,3a)') 'Image file ',imgfn,' created.'
  WRITE(6,'(/1x,a,e15.5,a,e15.5,a,i4,a,i4/)')                           &
        'zmax=',zx,' zmin=',zn,' imin=',imin,' imax=',imax

!  100  FORMAT (1X,(a,2X),1P,2G16.7,0P, 2X,2I4,1X,2F5.2)

END SUBROUTINE img2d
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE RASTERIZ                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE rasteriz (z, nx,ibgn,iend,isk0, ny,jbgn,jend,jsk0,           &
           zmax,zmin, image, nimage,imax,imin)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Convert a 2-D real array into a raster form.
!
!  This routine should be equivalent to the C routine  floattor8.c,
!  p. 2.10 of NCSA HDF manual.
!
!  Values are scaled to be in the range 1 to 254 (Icharmin to Icharmax).
!  Image is a 1-D array and must be stored upside-down.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  4/14/93
!
!  MODIFICATION HISTORY:
!
!  9/1/94 (Y. Lu)
!  Cleaned up documentation
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    z          a two array containing the field for which a image
!               file is to be created.
!    nx, ny     the dimensions of array z.
!    zmax,zmin  the maximum or minimum value of z
!
!  OUTPUT:
!
!    image      1-D array
!    imax,imin  the maximum and minimum value of image
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

  INTEGER :: nx, ny         ! dimensions of array z
  REAL :: z(nx,ny)       ! A two array containing the field for
                         ! which a image file is to be created.

  INTEGER :: ibgn,iend,isk  ! index of first and last i grid point;
                            ! skip/interpolation parameter in i-dir
  INTEGER :: jbgn,jend,jsk  ! index of first and last j grid point;
                            ! skip/interpolation parameter in j-dir

  INTEGER :: isk0,jsk0
  REAL :: zmax, zmin     ! the maximum or minimum value of z

  INTEGER :: imax, imin, nimage
                             ! the maximum and minimum value and
! dimension of imag
  CHARACTER (LEN=1) :: image(nimage)
!
!-----------------------------------------------------------------------
!
!  Misc. variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: icharmax,icharmin,icp1,iblack
  PARAMETER (icharmax=254, icharmin=2, icp1=icharmax+1, iblack=255)

  INTEGER :: i,j,k, ivalue, rtoi, isub, jsub
  REAL :: step, rvalue ,arg, zj1,zj2
!
!-----------------------------------------------------------------------
!
!  Inline funtion
!
!-----------------------------------------------------------------------
!
  rtoi( arg )  = nint ((arg-zmin)/step) + icharmin
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  step      = (zmax - zmin) / icp1

  IF (step <= 0.) THEN
    PRINT *, '.....Error in RASTERIZ:  Zmax <= Zmin'
    RETURN
  END IF


  imax = -32000
  imin =  32000

  isk = isk0
  jsk = jsk0

  IF( isk*jsk < 0 ) THEN
    IF( isk < 0 .AND. jsk == 1 ) THEN
      jsk = -1
    ELSE IF( jsk < 0 .AND. isk == 1 ) THEN
      isk = -1
    ELSE
      WRITE(6,'(/1x,a,/1x,a,i2,a,i2,/1x,a/)')                           &
          'Error in RASTERIZ: isk and jsk must have the same sign.',    &
          'isk = ', isk, ' jsk= ', jsk,                                 &
          'Job stopped in the subroutine.'
      CALL arpsstop('arpsstop called from RASTERIZ due to sign error',1)
    END IF
  END IF


  IF( isk > 0 .AND. jsk > 0 ) THEN

    k = 0
    DO j=jend,jbgn,-jsk
      DO i=ibgn,iend,isk
        rvalue = z(i,j)
        ivalue = MIN( MAX( rtoi(rvalue),icharmin ), icharmax )

        imax   = MAX (imax, ivalue)
        imin   = MIN (imin, ivalue)

        k      = k + 1
        image(k)= CHAR (ivalue)
      END DO
    END DO

    nimage = k

  ELSE IF( isk < 0 .AND. jsk < 0 ) THEN
                                  ! Case with binlinear interpolation.

    k = 0
    DO j=jend,jbgn+1,-1
      DO jsub=1,ABS(jsk)

        DO i=ibgn,iend-1
          DO isub=1,ABS(isk)
            zj1=z(i,j  )+(isub-1.0)/ABS(isk)*(z(i+1,j  )-z(i,j  ))
            zj2=z(i,j-1)+(isub-1.0)/ABS(isk)*(z(i+1,j-1)-z(i,j-1))
            rvalue = zj1+(jsub-1.0)/ABS(jsk)*(zj2-zj1)

            ivalue = MIN( MAX( rtoi(rvalue),icharmin ), icharmax )
            imax   = MAX (imax, ivalue)
            imin   = MIN (imin, ivalue)
            k      = k + 1
            image(k)= CHAR (ivalue)
          END DO
        END DO
        rvalue = z(iend,j)+(jsub-1.0)/ABS(jsk) *                        &
                (z(iend,j-1)-z(iend,j))
        ivalue = MIN( MAX( rtoi(rvalue),icharmin ), icharmax )
        imax   = MAX (imax, ivalue)
        imin   = MIN (imin, ivalue)
        k      = k + 1
        image(k)= CHAR (ivalue)
      END DO
    END DO

    DO i=ibgn,iend-1
      DO isub=1,ABS(isk)
        zj1=z(i,jbgn)+(isub-1.0)/ABS(isk)*(z(i+1,jbgn)-z(i,jbgn))
        rvalue = zj1

        ivalue = MIN( MAX( rtoi(rvalue),icharmin ), icharmax )
        imax   = MAX (imax, ivalue)
        imin   = MIN (imin, ivalue)
        k      = k + 1
        image(k)= CHAR (ivalue)
      END DO
    END DO

    rvalue = z(iend,jbgn)
    ivalue = MIN( MAX( rtoi(rvalue),icharmin ), icharmax )
    imax   = MAX (imax, ivalue)
    imin   = MIN (imin, ivalue)
    k      = k + 1
    image(k)= CHAR (ivalue)

    nimage = k

  END IF

  RETURN
END SUBROUTINE rasteriz
