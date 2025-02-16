!
!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE build_wrf_grid               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE build_wrf_grid(domid,parent_id,use_arps_grid,                &
                    dx_wrf_parent, dy_wrf_parent,                       &
                    xsub0_parent, ysub0_parent,                         &
                    i_parent_start, j_parent_start,                     &
                    nx_wrf,ny_wrf,nxlg_wrf,nylg_wrf,dx_wrf,dy_wrf,      &
                    mapproj_wrf,sclfct_wrf,lattru_wrf,lontru_wrf,       &
                    ctrlat_wrf,ctrlon_wrf,swx_wrf,swy_wrf,              &
                    xsub0, ysub0,                                       &
                    lat_wrf,lon_wrf,lat_ll,lat_ul,lat_ur,lat_lr,        &
                    lon_ll,lon_ul,lon_ur,lon_lr,istatus)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
!   Set up WRF mapprojection and get lat/lon of the WRF grid points.
!
!   Get lat/lon at the four corners of the WRF domain. After the call
!   Only processor 0 has valid lat/lon at corners.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: nx_wrf,   ny_wrf
  INTEGER, INTENT(IN)  :: nxlg_wrf, nylg_wrf
  REAL,    INTENT(IN)  :: dx_wrf,   dy_wrf
  INTEGER, INTENT(IN)  :: mapproj_wrf
  REAL,    INTENT(IN)  :: lattru_wrf(2),lontru_wrf
  REAL,    INTENT(INOUT) :: sclfct_wrf
  REAL,    INTENT(INOUT) :: ctrlat_wrf,ctrlon_wrf
  INTEGER, INTENT(INOUT) :: i_parent_start, j_parent_start
  REAL,    INTENT(INOUT) :: swx_wrf,swy_wrf
  INTEGER, INTENT(IN)  :: use_arps_grid
  INTEGER, INTENT(IN)  :: domid, parent_id
  REAL,    INTENT(IN)  :: dx_wrf_parent,dy_wrf_parent
  REAL,    INTENT(IN)  :: xsub0_parent, ysub0_parent
  REAL,    INTENT(OUT) :: xsub0,        ysub0
  REAL,    INTENT(OUT) :: lat_wrf(nx_wrf,ny_wrf,3)
  REAL,    INTENT(OUT) :: lon_wrf(nx_wrf,ny_wrf,3)

  REAL,    INTENT(OUT) :: lat_ll(4), lon_ll(4)
  REAL,    INTENT(OUT) :: lat_ul(4), lon_ul(4)
  REAL,    INTENT(OUT) :: lat_ur(4), lon_ur(4)
  REAL,    INTENT(OUT) :: lat_lr(4), lon_lr(4)

  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER :: i, j
  REAL :: sclf_wrf
  REAL :: dx_wrfscl, dy_wrfscl
  REAL :: ctrx, ctry, swx, swy
  REAL :: xsub, ysub

  INTEGER   :: mptag, source
  INTEGER, PARAMETER :: fzone_wrf = 1
  REAL,    PARAMETER :: a_small_number = 1.0E-5
    
  REAL, ALLOCATABLE ::  x_wrf(:),  y_wrf(:)
  REAL, ALLOCATABLE :: xs_wrf(:), ys_wrf(:)

  INCLUDE 'mp.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code ... ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF( sclfct_wrf /= 1.0) THEN
    sclf_wrf  = 1.0/sclfct_wrf
    dx_wrfscl = dx_wrf*sclf_wrf
    dy_wrfscl = dy_wrf*sclf_wrf
  ELSE
    sclf_wrf  = 1.0
    dx_wrfscl = dx_wrf
    dy_wrfscl = dy_wrf
  END IF
  sclfct_wrf = sclf_wrf
 
  CALL setmapr(mapproj_wrf,sclf_wrf,lattru_wrf,lontru_wrf)
  IF (domid == 1) THEN
    CALL lltoxy( 1,1, ctrlat_wrf,ctrlon_wrf, ctrx, ctry )
    swx_wrf = ctrx - 0.5*(nxlg_wrf-1) * dx_wrfscl
    swy_wrf = ctry - 0.5*(nylg_wrf-1) * dy_wrfscl
    xsub0 = 0.0
    ysub0 = 0.0
    CALL setorig( 1, swx_wrf, swy_wrf)
    WRITE(6,*)
  ELSE
    CALL setorig( 1, swx_wrf, swy_wrf)

    IF (use_arps_grid == 1) THEN   ! use ctrlat_wrf, ctrlon_wrf
      CALL lltoxy( 1,1, ctrlat_wrf,ctrlon_wrf, ctrx, ctry )
      swx = ctrx - 0.5*(nxlg_wrf-1) * dx_wrfscl
      swy = ctry - 0.5*(nylg_wrf-1) * dy_wrfscl
      xsub0 = swx 
      ysub0 = swy

      ctrx = ( xsub0-xsub0_parent ) / dx_wrf_parent + 1
      ctry = ( ysub0-ysub0_parent ) / dy_wrf_parent + 1
      i_parent_start = NINT(ctrx)
      j_parent_start = NINT(ctry)

      IF ( ABS(ctrx - i_parent_start) > a_small_number*i_parent_start .OR.   &
           ABS(ctry - j_parent_start) > a_small_number*j_parent_start ) THEN
        WRITE(6,'(1x,a,I2,2(a,F12.2),/,1x,2(a,F12.2),a,I2,/,2(1x,2(a,F12.2),a,/))')    &
        'Cannot align subdomain ',domid,' with central lat/lon = ',ctrlat_wrf,', ',ctrlon_wrf, &
        'and low-left corner = ',swx,', ',swy,', inside the parent grid ',parent_id,   &
        'with grid spacing = ',dx_wrf_parent,', ',dy_wrf_parent,'.',    &
        'The parent grid is starting from ',xsub0_parent,', ',ysub0_parent,'.'
        CALL arpsstop('Inconsistent nesting grid',1)
      END IF

      WRITE(6,'(14x,2(a,i6),a,/)')                                      &
      'Starting index are reset to: ',i_parent_start,                   &
                                 ', ',j_parent_start,'.'


    ELSE                           ! use i_parent_start, j_parent_start
      xsub0 = xsub0_parent + (i_parent_start-1)*dx_wrf_parent 
      ysub0 = ysub0_parent + (j_parent_start-1)*dy_wrf_parent 

      ctrx = xsub0 + 0.5*(nxlg_wrf-1)*dx_wrfscl
      ctry = ysub0 + 0.5*(nylg_wrf-1)*dy_wrfscl
      CALL xytoll( 1,1, ctrx, ctry, ctrlat_wrf, ctrlon_wrf)

      WRITE(6,'(14x,2(a,F12.2),a,/)')                                   &
              'Center lat/lon is reset to : ',ctrlat_wrf,               &
                                         ', ',ctrlon_wrf,'.'
    END IF
  END IF

  xsub = xsub0 + dx_wrf * (nx_wrf-fzone_wrf) * (loc_x-1)
  ysub = ysub0 + dy_wrf * (ny_wrf-fzone_wrf) * (loc_y-1)

!-----------------------------------------------------------------------
!
! Assign WRF grid arrays
!
!-----------------------------------------------------------------------

  ALLOCATE(x_wrf (nx_wrf), STAT = istatus)
  ALLOCATE(y_wrf (ny_wrf), STAT = istatus)
  ALLOCATE(xs_wrf(nx_wrf), STAT = istatus)
  ALLOCATE(ys_wrf(ny_wrf), STAT = istatus)

  DO i=1,nx_wrf
    x_wrf(i)= sclf_wrf*xsub + (i-1)*dx_wrfscl
  END DO
  DO j=1,ny_wrf
    y_wrf(j)= sclf_wrf*ysub + (j-1)*dy_wrfscl
  END DO
 
  DO i=1,nx_wrf-1
    xs_wrf(i)=0.5*(x_wrf(i)+x_wrf(i+1))
  END DO
  xs_wrf(nx_wrf)=2.*xs_wrf(nx_wrf-1)-xs_wrf(nx_wrf-2)

  DO j=1,ny_wrf-1
    ys_wrf(j)=0.5*(y_wrf(j)+y_wrf(j+1))
  END DO
  ys_wrf(ny_wrf)=2.*ys_wrf(ny_wrf-1)-ys_wrf(ny_wrf-2)

!-----------------------------------------------------------------------
!
!  Find latitude and longitude of WRF grid.
!
!   1 -- T point, 2 -- U point, 3 -- V point, 4 -- massless point
!
!-----------------------------------------------------------------------
 
  CALL xytoll(nx_wrf,ny_wrf,xs_wrf,ys_wrf,lat_wrf(:,:,1),lon_wrf(:,:,1))
  CALL xytoll(nx_wrf,ny_wrf, x_wrf,ys_wrf,lat_wrf(:,:,2),lon_wrf(:,:,2))
  CALL xytoll(nx_wrf,ny_wrf,xs_wrf, y_wrf,lat_wrf(:,:,3),lon_wrf(:,:,3))

!-----------------------------------------------------------------------
!
! Get lat/lon at each corners
!
!    ul - upper left corner           ur - upper right corner
!    ll - lower left corner           lr - lower right corner
!
!-----------------------------------------------------------------------

  lat_ll(1) = lat_wrf(       1,       1,1)
  lat_ul(1) = lat_wrf(       1,ny_wrf-1,1)
  lat_ur(1) = lat_wrf(nx_wrf-1,ny_wrf-1,1)
  lat_lr(1) = lat_wrf(nx_wrf-1,       1,1)

  lon_ll(1) = lon_wrf(       1,       1,1)
  lon_ul(1) = lon_wrf(       1,ny_wrf-1,1)
  lon_ur(1) = lon_wrf(nx_wrf-1,ny_wrf-1,1)
  lon_lr(1) = lon_wrf(nx_wrf-1,       1,1)

  lat_ll(2) = lat_wrf(       1,       1,2)
  lat_ul(2) = lat_wrf(       1,ny_wrf-1,2)
  lat_ur(2) = lat_wrf(nx_wrf,  ny_wrf-1,2)
  lat_lr(2) = lat_wrf(nx_wrf,         1,2)

  lon_ll(2) = lon_wrf(       1,       1,2)
  lon_ul(2) = lon_wrf(       1,ny_wrf-1,2)
  lon_ur(2) = lon_wrf(nx_wrf,  ny_wrf-1,2)
  lon_lr(2) = lon_wrf(nx_wrf,         1,2)

  lat_ll(3) = lat_wrf(       1,       1,3)
  lat_ul(3) = lat_wrf(       1,  ny_wrf,3)
  lat_ur(3) = lat_wrf(nx_wrf-1,  ny_wrf,3)
  lat_lr(3) = lat_wrf(nx_wrf-1,       1,3)

  lon_ll(3) = lon_wrf(       1,       1,3)
  lon_ul(3) = lon_wrf(       1,  ny_wrf,3)
  lon_ur(3) = lon_wrf(nx_wrf-1,  ny_wrf,3)
  lon_lr(3) = lon_wrf(nx_wrf-1,       1,3)

  CALL xytoll(1,1,x_wrf(1),     y_wrf(1),     lat_ll(4),lon_ll(4))
  CALL xytoll(1,1,x_wrf(1),     y_wrf(ny_wrf),lat_ul(4),lon_ul(4))
  CALL xytoll(1,1,x_wrf(nx_wrf),y_wrf(ny_wrf),lat_ur(4),lon_ur(4))
  CALL xytoll(1,1,x_wrf(nx_wrf),y_wrf(1),     lat_lr(4),lon_lr(4))

  IF (mp_opt > 0 ) THEN

    CALL inctag

    mptag  = gentag + 2                      ! passing UL to 0
    source = (nproc_y-1)*nproc_x
    IF (myproc == source)  THEN
      CALL mpsendr(lat_ul,4,0,mptag,  istatus)
      CALL mpsendr(lon_ul,4,0,mptag+1,istatus)
    END IF

    IF (myproc == 0) THEN
      CALL mprecvr(lat_ul,4,source,mptag,  istatus)
      CALL mprecvr(lon_ul,4,source,mptag+1,istatus)
    END IF

    mptag  = gentag + 4                     ! Pasing UR to 0
    source = (nproc_y-1)*nproc_x + nproc_x-1
    IF (myproc == source) THEN
      CALL mpsendr(lat_ur,4,0,mptag,  istatus)
      CALL mpsendr(lon_ur,4,0,mptag+1,istatus)
    END IF

    IF (myproc == 0) THEN
      CALL mprecvr(lat_ur,4,source,mptag,  istatus)
      CALL mprecvr(lon_ur,4,source,mptag+1,istatus)
    END IF

    mptag  = gentag + 6                     ! Passing LR to 0
    source = nproc_x-1
    IF (myproc == source) THEN
      CALL mpsendr(lat_lr,4,0,mptag,istatus)
      CALL mpsendr(lon_lr,4,0,mptag+1,istatus)
    END IF

    IF (myproc == 0) THEN
      CALL mprecvr(lat_lr,4,source,mptag,  istatus)
      CALL mprecvr(lon_lr,4,source,mptag+1,istatus)
    END IF

    CALL mpbarrier

  END IF

!-----------------------------------------------------------------------
!
! Returning
!
!-----------------------------------------------------------------------

  DEALLOCATE( x_wrf, y_wrf)
  DEALLOCATE(xs_wrf,ys_wrf)

  istatus = 0
  RETURN
END SUBROUTINE build_wrf_grid
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE prepinterp                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE prepinterp(nx,ny,nxlg,nylg,nx_wrf,ny_wrf,dx,dy,x,y,xs,ys,    &
                mapproj,scalef,latnot,trulon,ctrlat,ctrlon,swx,swy,     &
                lat_wrf,lon_wrf,x2d,y2d,iloc,jloc,                      &
                dxfld,rdxfld,dyfld,rdyfld,istatus)
!
!-----------------------------------------------------------------------
!
! PURPOSE:
!
!   Prepare horizontal interpolation arrays.
!
!   Do NOT support MPI.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: nx,     ny
  INTEGER, INTENT(IN)  :: nxlg,   nylg
  INTEGER, INTENT(IN)  :: nx_wrf, ny_wrf
  REAL,    INTENT(IN)  :: dx,     dy
  REAL,    INTENT(IN)  :: x(nx),  y(ny)
  REAL,    INTENT(IN)  :: xs(nx), ys(ny)
  INTEGER, INTENT(IN)  :: mapproj
  REAL,    INTENT(IN)  :: scalef
  REAL,    INTENT(IN)  :: latnot(2)
  REAL,    INTENT(IN)  :: trulon
  REAL,    INTENT(IN)  :: ctrlat,ctrlon
  REAL,    INTENT(OUT) :: swx,  swy
  REAL,    INTENT(IN)  :: lat_wrf(nx_wrf,ny_wrf,3)
  REAL,    INTENT(IN)  :: lon_wrf(nx_wrf,ny_wrf,3)
  REAL,    INTENT(OUT) :: x2d(nx_wrf,ny_wrf,3)
  REAL,    INTENT(OUT) :: y2d(nx_wrf,ny_wrf,3)
  INTEGER, INTENT(OUT) :: iloc(nx_wrf,ny_wrf,3)
  INTEGER, INTENT(OUT) :: jloc(nx_wrf,ny_wrf,3)
  REAL,    INTENT(OUT) :: dxfld(nx,3)
  REAL,    INTENT(OUT) :: rdxfld(nx,3)
  REAL,    INTENT(OUT) :: dyfld(ny,3)
  REAL,    INTENT(OUT) :: rdyfld(ny,3)

  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER  :: n

  INTEGER  :: iproj
  REAL     :: scl, trlat(2), trlon,x0,y0

  REAL     :: ctrx, ctry

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code ... ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  Save previous map projection values
!
!-----------------------------------------------------------------------
!
  CALL getmapr(iproj,scl,trlat,trlon,x0,y0)

!-----------------------------------------------------------------------
!
!  Set ARPS map projection
!
!-----------------------------------------------------------------------

  CALL setmapr(mapproj,scalef,latnot,trulon)
  CALL lltoxy( 1,1, ctrlat,ctrlon, ctrx, ctry )
  swx = ctrx - (FLOAT((nxlg-3))/2.) * dx*scalef
  swy = ctry - (FLOAT((nylg-3))/2.) * dy*scalef
  CALL setorig(1,swx,swy)

  DO n = 1,3   ! For T, U, V three staggers
    CALL lltoxy(nx_wrf,ny_wrf,lat_wrf(:,:,n),lon_wrf(:,:,n),        &
                x2d(:,:,n),y2d(:,:,n))
  END DO

  IF(MAXVAL(x2d) > x(nx) .OR. MINVAL(x2d) < x(1) .OR.          &
     MAXVAL(y2d) > y(ny) .OR. MINVAL(y2d) < y(1) ) THEN

    WRITE(6,'(/,a,/,2(4(a,F15.2),a,/))')                                &
     'It seems that WRF domain is outside of the ARPS physical domain.',&
     'ARPS domain [',x(1),' -- ',x(nx),' ; ',y(1), ' -- ', y(ny), '  ]',&
     'WRF  domain [',MINVAL(x2d),' -- ', MAXVAL(x2d), ' ; ',            &
                     MINVAL(y2d),' -- ', MAXVAL(y2d), '  ]'

    CALL arpsstop('Domain error.',1)
  END IF

  CALL setijloc(nx_wrf,ny_wrf,nx,ny,x2d(:,:,1),y2d(:,:,1),          &
                xs,ys,iloc(:,:,1),jloc(:,:,1))  ! T points
  CALL setijloc(nx_wrf,ny_wrf,nx,ny,x2d(:,:,2),y2d(:,:,2),          &
                x,ys,iloc(:,:,2),jloc(:,:,2))   ! U points
  CALL setijloc(nx_wrf,ny_wrf,nx,ny,x2d(:,:,3),y2d(:,:,3),          &
                xs,y,iloc(:,:,3),jloc(:,:,3))   ! V points

  CALL setdxdy(nx,ny,1,nx,1,ny,xs,ys,                               &
               dxfld(:,1),dyfld(:,1),rdxfld(:,1),rdyfld(:,1))
  CALL setdxdy(nx,ny,1,nx,1,ny,x,ys,                                &
               dxfld(:,2),dyfld(:,2),rdxfld(:,2),rdyfld(:,2)) 
  CALL setdxdy(nx,ny,1,nx,1,ny,xs,y,                                &
               dxfld(:,3),dyfld(:,3),rdxfld(:,3),rdyfld(:,3))
!
!-----------------------------------------------------------------------
!
!  Reset map projection to previous values
!
!-----------------------------------------------------------------------
!
  CALL setmapr(iproj,scl,trlat,trlon)
  CALL setorig(1,x0,y0)

  istatus = 0

  RETURN
END SUBROUTINE prepinterp
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE hinterp                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE hinterp(nx,ny,nz,nx_ext,ny_ext,iorder,                    &
                   iloc,jloc,x,y,x2d,y2d,                            &
                   varin,dxfld,dyfld,rdxfld,rdyfld,                  &
                   varout,slopey,alphay,betay)

!-----------------------------------------------------------------------
!  PURPOSE:
! 
!    Interpolates ARPS 3D array varin horizontally to WRF or other
!    grids. The output array varout is a 3D interplation array with
!    the same vertical coordinate as the input array varin. 
! 
!  NOTE:
!    Refer to function pntint2d in src/adas/intfield.f90.
!
!-----------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nx,ny,nz            ! Grid size for input array
  INTEGER, INTENT(IN) :: nx_ext,ny_ext       ! Horizontal grid size of the 
                                             ! interpolated array
  INTEGER, INTENT(IN) :: iorder              ! polynomial order of the interpolation
                                             ! = 1 bi-linear
                                             ! > 2 quadratic
  INTEGER, INTENT(IN) :: iloc(nx_ext,ny_ext) ! ARPS index of external grid points
  INTEGER, INTENT(IN) :: jloc(nx_ext,ny_ext) ! ARPS index of external grid points
  REAL,    INTENT(IN) :: x(nx), y(ny)        ! Grid coordiate of the input array
  REAL,    INTENT(IN) :: x2d(nx_ext,ny_ext), y2d(nx_ext,ny_ext)
             ! x and y coordinate of external grid point in ARPS coordinate
  REAL,    INTENT(IN) :: varin(nx,ny,nz)
  REAL,    INTENT(IN) :: dxfld(nx),rdxfld(nx)
  REAL,    INTENT(IN) :: dyfld(ny),rdyfld(ny)

  REAL,    INTENT(OUT) :: varout(nx_ext,ny_ext,nz)
  REAL,    INTENT(OUT) :: slopey(nx,ny,nz)
  REAL,    INTENT(OUT) :: alphay(nx,ny,nz)
  REAL,    INTENT(OUT) :: betay(nx,ny,nz)

!-----------------------------------------------------------------------
!
! Misc. local variables
!
!-----------------------------------------------------------------------
  INTEGER :: i,j,k 

  INTEGER :: ii,jj
  REAL    :: delx,dely
  REAL    :: varm1,var00,varp1
  REAL    :: rtwodx, alpha, beta

!----------------------------------------------------------------------
!
! External interpolation function
! Interpolate a 2-d field for a single point on that plane.
!
!----------------------------------------------------------------------
!  REAL :: pntint2d             ! the function doing horizontal 
                                ! interpolation for one grid point
!
! On some platform, it is hard to make pntint2d inline for efficiency.
! So we hardcoded the function pntint2d inside this subroutine.
!
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
!  Compute derivative terms
!
!-----------------------------------------------------------------------
!
  CALL setdrvy(nx,ny,nz,1,nx,1,ny,1,nz,                                 &
               dyfld,rdyfld,varin,slopey,alphay,betay)

!
!-----------------------------------------------------------------------
!
!    Horizontal interpolation
!
!-----------------------------------------------------------------------
!
  IF(iorder == 1) THEN
!
!-----------------------------------------------------------------------
!
!  Loop through all WRF grid points
!
!-----------------------------------------------------------------------
!
    DO k = 1,nz
      DO j = 1,ny_ext
        DO i = 1,nx_ext
!
!-----------------------------------------------------------------------
!
!  Compute bilinear interpolated value
!
!-----------------------------------------------------------------------
!
           ii=MIN(MAX(iloc(i,j),2),(nx-2))
           jj=MIN(MAX(jloc(i,j),2),(ny-2))
           delx=(x2d(i,j)-x(ii))
           dely=(y2d(i,j)-y(jj))
           varout(i,j,k)= (1.-delx*rdxfld(ii))*                         &
                          (varin(ii,  jj,k)+slopey(ii,  jj,k)*dely)     &
                        + (delx*rdxfld(ii))*                            &
                          (varin(ii+1,jj,k)+slopey(ii+1,jj,k)*dely)
        END DO
      END DO
    END DO

  ELSE
!
!-----------------------------------------------------------------------
!
!  Loop through all WRF grid points
!
!-----------------------------------------------------------------------
!
    DO k = 1,nz
      DO j = 1,ny_ext
        DO i = 1,nx_ext
!
!-----------------------------------------------------------------------
!
!  Compute biquadratic
!
!-----------------------------------------------------------------------
!
           ii=MIN(MAX(iloc(i,j),(2+1)),(nx-2))
           jj=MIN(MAX(jloc(i,j),(2+1)),(ny-2))
           delx=(x2d(i,j)-x(ii))
           dely=(y2d(i,j)-y(jj))
!
!-----------------------------------------------------------------------
!
!    Stencil is ii-1 to ii+1 and jj-1 to jj + 1
!
!    Interpolate in y.
!
!-----------------------------------------------------------------------
!
          varm1=(alphay(ii-1,jj,k)*dely+betay(ii-1,jj,k))*dely          &
               + varin(ii-1,jj,k)
          var00=(alphay(ii  ,jj,k)*dely+betay(ii  ,jj,k))*dely          &
               + varin(ii  ,jj,k)
          varp1=(alphay(ii+1,jj,k)*dely+betay(ii+1,jj,k))*dely          &
               + varin(ii+1,jj,k)
!
!-----------------------------------------------------------------------
!
!    Interpolate intermediate results in x.
!
!-----------------------------------------------------------------------
!
          rtwodx=1./(dxfld(ii-1)+dxfld(ii))
          alpha=((varp1-var00)*rdxfld(ii  )                            &
               + (varm1-var00)*rdxfld(ii-1))*rtwodx
          beta= (varp1-var00)*rdxfld(ii)                               &
              - dxfld(ii)*alpha
          varout(i,j,k)=(alpha*delx+beta)*delx+var00
        END DO
      END DO
    END DO

  END IF

!         varout(i,j,k) = pntint2d(nx,ny,2,nx-1,2,ny-1,                  &
!               iorder,x,y,x2d(i,j),y2d(i,j),iloc(i,j),jloc(i,j),        &
!               varin(:,:,k),dxfld,dyfld,rdxfld,rdyfld,                  &
!               slopey(:,:,k),alphay(:,:,k),betay(:,:,k))

  RETURN
END SUBROUTINE hinterp
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE compute_eta_3d              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE compute_eta_3d(nx,ny,nz,p,t,mix,zp,ztopo,ptop,eta,mu)

!-------------------------------------------------------------------
!
! PURPOSE
!
! Subroutine to compute eta on a 3D grid.  Also computes mu as a 
! function of the terrain passed in.
!
! The Eta-P coordinate.
!
!                Pdry - Ptop       Pdry = Dry pressure (Pa)
!         EtaP = ----------        Ptop = Pressure at model top (Pa)
!                    mu            mu   = Psfc - Pvapor_bottom - Ptop
!
!  General Procedure:
!
!     1.  Compute downward integrated vapor-pressure for each gridpoint
!          ( Pvapor)
!     2.  Subtract the Pvapor value from the pressure at each point to
!         get the dry pressure (Pdry) at each point
!     3.  Compute mu (dry surface pressure minus top pressure)
!     4.  Compute EtaP value using newly compute dry pressure for each
!         point
!
! NOTE:
!   refer to subroutine compute_eta_3d in WRFSI src/mod/module_vinterp_utils.F
!   ARPS physical domain is from 2 to nz-1.
!   Water mixing ratio "mix" is passed in instead of RH.
!
!-------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx,ny        ! Horizontal dimensions of array 
                                      ! in WRF grid
  INTEGER, INTENT(IN) :: nz           ! Horizontal dimensions of array 
                                      ! in ARPS grid
  REAL,    INTENT(IN) :: p(nx,ny,nz)  ! 3D Pressure (Pa)
  REAL,    INTENT(IN) :: t(nx,ny,nz)  ! 3D Temperature (K)
  REAL,    INTENT(IN) :: mix(nx,ny,nz)! 3D water vapor mixing ratio
  REAL,    INTENT(IN) :: zp(nx,ny,nz) ! 3D GPH (gpm)
  REAL,    INTENT(IN) :: ztopo(nx,ny) ! WRF Topography
  REAL,    INTENT(IN) :: ptop         ! Top pressure (Pa)

  REAL,    INTENT(OUT) :: eta(nx,ny,nz) ! Eta values
  REAL,    INTENT(OUT) :: mu(nx,ny)     ! Mu

  INCLUDE 'mp.inc'  ! require variable myproc

!------------------------------------------------------------------
!
! Misc. Local variables
!
!------------------------------------------------------------------
  REAL, ALLOCATABLE :: rho(:)
  REAL, ALLOCATABLE :: pdry(:)
  REAL, ALLOCATABLE :: iqvp(:)

  REAL    :: tsfc, psfc, qvsfc, dz, qvbar,iqvp_sfc
  REAL    :: wgt0, wgt1
  LOGICAL :: found_level

  INTEGER :: i,j,k

  REAL, PARAMETER :: g = 9.81

  INTEGER :: imid, jmid, imidproc, jmidproc

!------------------------------------------------------------------
!
! External functions
!
!------------------------------------------------------------------
  REAL :: compute_density  ! function to compute grid point density
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  imid = ( (nx-3)*nproc_x+3 )/2     ! Added by Y. Wang for proper diagnostic outputs
  jmid = ( (ny-3)*nproc_y+3 )/2
  imidproc = (imid-2) / (nx-3) + 1
  jmidproc = (jmid-2) / (ny-3) + 1

  IF (loc_x == imidproc) THEN
    imid = MOD((imid-2),(nx-3)) + 2   ! Local index for global middle point
  ELSE
    jmid = -999
  END IF

  IF (loc_y == jmidproc) THEN
    jmid = MOD((jmid-2),(ny-3)) + 2   ! Local index for global middle point
  ELSE
    jmid = -999
  END IF

  ALLOCATE (rho (nz))      ! Density for each level in a column
  ALLOCATE (pdry(nz))      ! Dry air pressure (iqvp removed)
  ALLOCATE (iqvp(nz))      ! Integrated (downward) vapor pressure

!-----------------------------------------------------------------------
!
! Loops over horizontal grid points
!
!-----------------------------------------------------------------------

  DO j = 1, ny
    DO i = 1, nx

      ! Compute density and qv in the column
      DO k = 1,nz
        rho(k)= compute_density(t(i,j,k),p(i,j,k))
      ENDDO

      ! Compute the dry pressure values and iqvp for the column
      ! (pdry+iqvp = p)
      CALL compute_column_pdry(nz,p(i,j,:),rho,mix(i,j,:),zp(i,j,:),    &
                               pdry,iqvp)

      ! Compute mu at this point using the terrain height. Mu is just
      ! the dry surface pressure, which would be the model-adjusted
      ! surface pressure minus the integrated water vapor.

      IF (ztopo(i,j) .LE. zp(i,j,2) ) THEN  ! ARPS physical domain is 
                                            ! 2 -> (nz-1)

        ! WRF terrain below lowest BG model level.  Make crude adjustment
        ! using 10 Pa per meter to get surface pressure and assume a
        ! constant qv for moisture.  Assume a standard 6.5 K/km lapse rate
        ! for temperature
        k = 2
        dz    = zp(i,j,2)-ztopo(i,j)
        psfc  = p(i,j,2) + dz*10.
        qvbar = mix(i,j,2)
        tsfc  = t(i,j,2) + dz*.0065
      ELSE
        found_level = .false.
        find_level: DO k = 3, nz-1
          IF ( (ztopo(i,j) .LE. zp(i,j,k)) .AND.                        &
               (ztopo(i,j) .GT. zp(i,j,k-1)) ) THEN
            found_level = .true.
            EXIT find_level
          END IF
        END DO find_level
        IF (.NOT. found_level) THEN
          PRINT *, 'COMPUTE_ETA_3D: This should not happen.'
          PRINT *, 'Could not find where to insert WRF surface.'
          PRINT *, 'I/J = ',i,j
          PRINT *, 'Topo = ', ztopo(i,j)
          PRINT *, 'Z in column = ', zp(i,j,:)
          STOP 'level_problem'
        END IF

        ! Determine pressure, temperature, and qv
        psfc = EXP ( ( ztopo(i,j)*ALOG(p(i,j,k-1)/p(i,j,k)) -           &
                       zp(i,j,k)*ALOG(p(i,j,k-1))           +           &
                       zp(i,j,k-1)*ALOG(p(i,j,k))  )          /         &
                     ( zp(i,j,k-1) - zp(i,j,k)     )    )
        IF ( (zp(i,j,k)-zp(i,j,k-1)) .GE. 1.) THEN
          wgt0 = (zp(i,j,k)-ztopo(i,j)) / (zp(i,j,k)-zp(i,j,k-1))
          wgt1 = 1 - wgt0
        ELSE
          wgt0 = 0.0
          wgt1 = 1.0
        ENDIF

        tsfc = t(i,j,k-1)*wgt0 + t(i,j,k)*wgt1
        qvsfc = mix(i,j,k-1)*wgt0 + mix(i,j,k)*wgt1
        qvbar = (qvsfc + mix(i,j,k))*0.5
        dz = zp(i,j,k) - ztopo(i,j)
      ENDIF

      ! Now compute the integrated vapor pressure between
      ! the surface and the next lowest level 
      iqvp_sfc = iqvp(k)+g*qvbar*rho(k)*dz/(1.+qvbar)
      mu(i,j) = psfc - iqvp_sfc - ptop

      ! We now have everything we need to compute eta for the column
      DO k = 1,nz
        eta(i,j,k) = (pdry(k) - ptop) / mu(i,j)
      END DO

      IF ((i == imid) .AND. (j == jmid)) THEN
        ! Some diagnostic prints from the center of the column
        PRINT *, ' '
        PRINT *, 'Eta values on background grid center column:'
        PRINT *, 'TOPO/PSFC/MU = ',ztopo(i,j), psfc, mu(i,j)
        PRINT *, 'Z          PRESS      PRESSDRY   IQVAPORP   ETA'
        PRINT *, '---------- ---------- ---------- ---------- ----------'
        DO k=1,nz-1
          PRINT '(F10.1,1x,F10.1,1x,F10.1,1x,F10.6,1x,F10.6)',zp(i,j,k), &
            p(i,j,k), pdry(k),iqvp(k),eta(i,j,k)
        END DO
        PRINT *, ' '
      END IF

    END DO
  END DO

!-----------------------------------------------------------------------
!
! Just before returning
!
!-----------------------------------------------------------------------

  DEALLOCATE (rho)
  DEALLOCATE (pdry)
  DEALLOCATE (iqvp)

  RETURN
END SUBROUTINE compute_eta_3d

!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE vinterp                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE vinterp(nx,ny,nz_in,nz_out,korder,eta_in, eta_out,           &
                   data_in, data_out, extrapolate)

!-------------------------------------------------------------------
!
! PURPOSE:
!
!   Interpolates each column of data from ARPS height surface 
!   to WRF mass vertical coordinate.
!
! NOTE:
!   Refer to subroutine interp_eta2eta_lin in WRFSI src/mod/module_vinterp_utils.F
!    o ARPS vertical physical domain is from 2 to nz-1
!    o We added two supported interpolations method, linear and 
!      quadratic while WRFSI only support linear interoplation so far.
!
!-------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx      ! X dimension of WRF grid
  INTEGER, INTENT(IN) :: ny      ! Y dimension of WRF grid
  INTEGER, INTENT(IN) :: nz_in   ! Vertical dimension of ARPS grid
  INTEGER, INTENT(IN) :: nz_out  ! Vertical dimension of WRF grid (mass)
  INTEGER, INTENT(IN) :: korder  ! polyn. interpolation order
  REAL,    INTENT(IN) :: eta_in(nx,ny,nz_in)
                                 ! mass value at ARPS grid
  REAL,    INTENT(IN) :: eta_out(nz_out)
                                 ! Desired mass levels
  REAL,    INTENT(IN) :: data_in(nx,ny,nz_in)
  REAL,    INTENT(OUT):: data_out(nx,ny,nz_out)
  LOGICAL, INTENT(IN) :: extrapolate
                                 ! Whether to do extrapolation

!------------------------------------------------------------------
!
! External function
!
!------------------------------------------------------------------
  REAL :: intpnt_lin             ! do vertical interpolation linearly
  REAL :: intpnt_quad            ! do vertical interpolation quadratically
!
!-----------------------------------------------------------------
!
! Misc. Local variables
!
!-----------------------------------------------------------------

  INTEGER :: i,j,k,kk, ki
  REAL    :: desired_eta
  REAL    :: dvaldeta
  INTEGER :: nz_top

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  data_out(:,:,:) = -99999.9
  nz_top = nz_in - 1

  DO j = 1, ny
    DO i = 1, nx
      output_loop: DO k = 1, nz_out

        desired_eta = eta_out(k)
        IF (desired_eta .GT. eta_in(i,j,2)) THEN

          IF ((desired_eta - eta_in(i,j,2)).LT. 0.0001) THEN
             data_out(i,j,k) = data_in(i,j,2)
          ELSE
            IF (extrapolate) THEN
              ! Extrapolate downward because desired eta level is below
              ! the lowest level in our input data.  Extrapolate using simple
              ! 1st derivative of value with respect to eta for the bottom 2
              ! input layers.

              ! Add a check to make sure we are not using the gradient of
              ! a very thin layer

              IF ( (eta_in(i,j,2)-eta_in(i,j,3)) .GT. 0.001) THEN
                dvaldeta = (data_in(i,j,2) - data_in(i,j,3)) /          &
                           (eta_in(i,j,2) - eta_in(i,j,3) )
              ELSE
                dvaldeta = (data_in(i,j,2) - data_in(i,j,4)) /          &
                            (eta_in(i,j,2) - eta_in(i,j,4) )
              ENDIF
              data_out(i,j,k) = data_in(i,j,2) +                        &
                              dvaldeta * (desired_eta-eta_in(i,j,2))
            ELSE
              data_out(i,j,k) = data_in(i,j,2)

            ENDIF
          ENDIF
        ELSE IF (desired_eta .LE. eta_in(i,j,nz_top)) THEN
          IF ( abs(eta_in(i,j,nz_top) - desired_eta) .LT. 0.0001) THEN
             data_out(i,j,k) = data_in(i,j,nz_top)
          ELSE
            IF (extrapolate) THEN
              ! Extrapolate upward
              IF ( (eta_in(i,j,nz_top-1)-eta_in(i,j,nz_top)) .GT. 0.0005) THEN
                dvaldeta = (data_in(i,j,nz_top) - data_in(i,j,nz_top-1)) / &
                           (eta_in(i,j,nz_top)  - eta_in(i,j,nz_top-1))
              ELSE
                dvaldeta = (data_in(i,j,nz_top) - data_in(i,j,nz_top-2)) / &
                           (eta_in(i,j,nz_top)  - eta_in(i,j,nz_top-2))
              ENDIF
              data_out(i,j,k) =  data_in(i,j,nz_top) +                     &
                                 dvaldeta * (desired_eta-eta_in(i,j,nz_top))
            ELSE
              data_out(i,j,k) = data_in(i,j,nz_top)
            ENDIF
          ENDIF
        ELSE
          ! We can trap between two levels and linearly interpolate

          input_loop:  DO kk = 2, nz_top-1
            IF (desired_eta .EQ. eta_in(i,j,kk) )THEN
              data_out(i,j,k) = data_in(i,j,kk)
              EXIT input_loop
            ELSE IF ( (desired_eta .LT. eta_in(i,j,kk)) .AND.           &
                      (desired_eta .GT. eta_in(i,j,kk+1)) ) THEN
              IF(korder == 1) THEN      
                data_out(i,j,k) = intpnt_lin( eta_in(i,j,kk),           &
                                    eta_in(i,j,kk+1), desired_eta,      &
                                    data_in(i,j,kk), data_in(i,j,kk+1) )
              ELSE IF(korder == 2)  THEN
                IF ( ABS(desired_eta - eta_in(i,j,kk)) <=               &
                     ABS(desired_eta - eta_in(i,j,kk+1)) ) THEN
                  ki = kk
                ELSE IF ( kk == nz_top-1 ) THEN
                  ki = kk
                ELSE
                  ki = kk+1
                END IF
                data_out(i,j,k) = intpnt_quad( eta_in(i,j,ki-1),        &
                                    eta_in(i,j,ki),eta_in(i,j,ki+1),    &
                                    desired_eta, data_in(i,j,ki-1),     &
                                    data_in(i,j,ki),data_in(i,j,ki+1) )
              ELSE
                WRITE(6,*) 'Unsupported polynomial interpolation order:',korder
                STOP 'Unsupported korder'
              END IF

              EXIT input_loop
            ENDIF
          ENDDO input_loop
        ENDIF

      ENDDO output_loop

    ENDDO   ! i
  ENDDO   ! j

  RETURN
END SUBROUTINE vinterp

!
!##################################################################
!##################################################################
!######                                                      ######
!######                  FUNCTION intpnt_lin                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
FUNCTION intpnt_lin(x1,x2,x,d1,d2) RESULT(d)

!------------------------------------------------------------------
!
! PURPOSE:
!   Fuction to do vertical interpolation linearly
!
!-----------------------------------------------------------------------
!
! Author: Yunheng Wang.
!  
!-----------------------------------------------------------------

  IMPLICIT NONE

  REAL, INTENT(IN)  :: x1    ! coordinate of first point
  REAL, INTENT(IN)  :: x2    ! coordinate of second point
  REAL, INTENT(IN)  :: x     ! coordinate of desired point
  REAL, INTENT(IN)  :: d1    ! data at first point
  REAL, INTENT(IN)  :: d2    ! data at second point
  REAL              :: d     ! Interpolated value
  
!------------------------------------------------------------------
!
!  Misc. local variables
!
!------------------------------------------------------------------

  REAL    :: wgt

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  wgt = (x - x2) / (x1 - x2)
  d   = wgt*d1 + (1.-wgt)*d2

  RETURN
END FUNCTION intpnt_lin

!##################################################################
!##################################################################
!######                                                      ######
!######                  FUNCTION intpnt_quad                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
FUNCTION intpnt_quad(x0,x1,x2,x,d0,d1,d2) RESULT(d)

!-------------------------------------------------------------------
!
! PURPOSE:
!
!   Quadratically interpolates from the given three points 
!
!-----------------------------------------------------------------------
!
! Author: Yunheng Wang.
!
!-------------------------------------------------------------------

  IMPLICIT NONE

  REAL, INTENT(IN) :: x0     ! coordinate of the first point
  REAL, INTENT(IN) :: x1     ! coordinate of the second point
  REAL, INTENT(IN) :: x2     ! coordinate of the third point
  REAL, INTENT(IN) :: x      ! coordinate of the desired point
  REAL, INTENT(IN) :: d0     ! value at the first point
  REAL, INTENT(IN) :: d1     ! value at the second point
  REAL, INTENT(IN) :: d2     ! value at the third point

  REAL             :: d      ! value to be calculated
!
!-----------------------------------------------------------------
!
! Misc. Local variables
!
!-----------------------------------------------------------------

  REAL    :: b1,b2, tmp

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
!  Quadratic interpolation
!
!             f(x) = f(x0) + b1*(x-x0) + b2*(x-x0)(x-x1)
!
!                    f(x1)-f(x0)
!             b1   = ----------
!                     x1 - x0
!
!                    f(x2)-f(x1)   f(x1)-f(x0)
!                    ----------- - ------------
!                     x2 - x1        x1 - x0
!             b2   = -----------------------------
!                             x2 - x0
!
!------------------------------------------------------------------------
!
  b1  = (d1  - d0) / (x1 - x0)
  tmp = (d2  - d1) / (x2 - x1)
  b2  = (tmp - b1) / (x2 - x0)
  
  d   = d0 + b1*(x-x0) + b2*(x-x0)*(x-x1)

  RETURN
END FUNCTION intpnt_quad

!
!##################################################################
!##################################################################
!######                                                      ######
!######                  FUNCTION compute_density            ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
FUNCTION compute_density(t_k, p_pa) RESULT(rho)

!-----------------------------------------------------------------
!
! PURPOSE:
! Computes density (kg/m{-3}) given temperature (K) and pressure (Pa)
!
! NOTE:
!   Refer to the same function in WRFSI src/mode/module_diagnostic_vars.F.
!
!-----------------------------------------------------------------

  IMPLICIT NONE
  REAL, INTENT(IN) :: t_k
  REAL, INTENT(IN) :: p_pa
  REAL             :: rho

  REAL, PARAMETER  :: Rd =  287.04  ! J deg{-1} kg{-1}
                                    ! Gas constant for dry air

  rho = p_pa / (Rd * t_k)

  RETURN
END FUNCTION compute_density

!
!##################################################################
!##################################################################
!######                                                      ######
!######               FUNCTION compute_column_pdry           ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE compute_column_pdry(nz,p,rho,qv,z,pdry,pvapor)

!-----------------------------------------------------------------
!
! PURPOSE:
!
! Given full pressure in Pascals (p), dry air density in kg/m^3 (rho),
! water vapor mixing ratio in kg/kg (qv), compute the dry pressure value
! at each grid point in the column of nz elements.
! This routine assumes the pressure array decreases monotonically from
! 1->nz and that the rho and qv arrays match the pressure array.
!
! NOTE:
!   Refer to subroutine compute_column_pdry in WRFSI src/mod/module_vinterp_utils.F.
!
!-----------------------------------------------------------------

  IMPLICIT NONE
  INTEGER, INTENT(IN)                    :: nz
  REAL, INTENT(IN)                       :: p  (nz)
  REAL, INTENT(IN)                       :: rho(nz)
  REAL, INTENT(IN)                       :: qv (nz)
  REAL, INTENT(IN)                       :: z  (nz)
  REAL, INTENT(OUT)                      :: pdry  (nz)
  REAL, INTENT(OUT)                      :: pvapor(nz)

  INTEGER                                :: kbot, k
  REAL                                   :: qv_mean,dz

  REAL, PARAMETER :: g = 9.8

! Set top vapor pressure to zero and top dry pressure equal to
! top total pressure

    pvapor(nz) = 0.
    pdry(nz) = p(nz)

    ! Integrate moisture downward

    main_loop:  DO kbot = nz-1, 1, -1

      ! Initialize for upcoming sums.
      pvapor(kbot) = 0.

      ! Integrate downward
      down_to_here:  DO k = nz-1,kbot,-1
        ! Compute delta-Z and mean Qv for this layer
        dz = z(k+1) - z(k)
        qv_mean = (qv(k) + qv(k+1)) * 0.5

        ! Compute pvapor for this layer and sum with previous layer
        pvapor(kbot) = pvapor(kbot)+g*qv_mean*rho(k)*dz/(1.+qv_mean)
      END DO down_to_here

      ! Subtract the integrated vapor pressure from the total pressure
      pdry(kbot) = p(kbot) - pvapor(kbot)
    END DO main_loop

END SUBROUTINE compute_column_pdry
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE compute_ph                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE compute_ph(nx,ny,nz,ht,znw,ptop,mu_input,                    &
                      pt_input,qv,mu,mub,ph,phb,pt,p,pb,                &
                      t_init,alb,alt,al)
!------------------------------------------------------------------
!
! PURPOSE:
!
!   Compute WRF variables, such as geopotential, air mass and potential
!   temperature etc.
!
! NOTE:
!   Refered to subroutine init_domain_rk in WRF dyn_em/module_initialize_real.f.
!
!------------------------------------------------------------------

  USE wrf_metadata

  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: nx,ny,nz     ! WRF dimension size
  REAL,    INTENT(IN)  :: ht(nx,ny)    ! Topographic height
  REAL,    INTENT(IN)  :: znw(nz)      ! full vertical levels
  REAL,    INTENT(IN)  :: ptop         ! model top pressure (Pa)
  REAL,    INTENT(IN)  :: mu_input(nx,ny)    ! Total air mass in column
  REAL,    INTENT(IN)  :: pt_input(nx,ny,nz) ! potential temperature
  REAL,    INTENT(IN)  :: qv(nx,ny,nz) ! water vapor mixing ratio
  
  REAL,    INTENT(OUT) :: ph(nx,ny,nz) ! base-state geopotential
  REAL,    INTENT(OUT) :: phb(nx,ny,nz)! perturbation geopotential
  REAL,    INTENT(OUT) :: mu(nx,ny)    ! perturbation dry air mass in column
  REAL,    INTENT(OUT) :: mub(nx,ny)   ! base state dry air mass in column
  REAL,    INTENT(OUT) :: pt(nx,ny,nz) ! perturbation potential temperature
  REAL,    INTENT(OUT) :: p(nx,ny,nz)
  REAL,    INTENT(OUT) :: pb(nx,ny,nz) ! base state pressure
                        ! work arrays below
  REAL,    INTENT(OUT) :: t_init(nx,ny,nz) ! 
  REAL,    INTENT(OUT) :: alb(nx,ny,nz)    ! inverse of reference density
  REAL,    INTENT(OUT) :: alt(nx,ny,nz)    ! inverse of perturb. density
  REAL,    INTENT(OUT) :: al(nx,ny,nz)     ! inverse of density
!------------------------------------------------------------------
!
!  Misc. local variables
!
!------------------------------------------------------------------

  REAL    :: p00,t00,a
  REAL    :: p_surf                ! surface pressure
  REAL    :: qvf1,qvf2,qvf

  REAL    :: znu(nz-1)             ! half level locations
  REAL    :: dnw(nz-1)             ! layer thickness
  REAL    :: rdn(nz-1)

  INTEGER :: i,j,k

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  DO k = 1, nz-1
    znu(k) = 0.5*( znw(k+1) + znw(k) )
    dnw(k) = znw(k+1) - znw(k)
  END DO
  DO k = 2, nz-1
    rdn(k) = 1.0/(0.5*(dnw(k)+dnw(k-1)))
  END DO

  !  To define the base state, we call a USER MODIFIED routine to set the three
  !  necessary constants:  p00 (sea level pressure, Pa), t00 (sea level temperature, K),
  !  and A (temperature difference, from 1000 mb to 300 mb, K).

  CALL const_module_initialize ( p00 , t00 , a )

  ! bast state pressure is a function of eta level and terrain only

  DO j = 1, ny-1
    DO i = 1,nx-1
      p_surf = p00*EXP(-t00/a + ((t00/a)**2 - 2.*g_wrf*ht(i,j)/a/r_d)**0.5)

      DO k = 1,nz-1
        pb(i,j,k) = znu(k)*(p_surf-ptop) + ptop
        t_init(i,j,k) = (t00 + a*LOG(pb(i,j,k)/p00))*(p00/pb(i,j,k))**(r_d/cp_wrf) - t0
        alb(i,j,k) = (r_d/p1000mb)*(t_init(i,j,k)+t0)*(pb(i,j,k)/p1000mb)**cvpm
      END DO
      
      mub(i,j) = p_surf - ptop            ! base state dry air mass
      mu(i,j)  = mu_input(i,j) - mub(i,j) ! perturbation dry air mass

      phb(i,j,1) = g_wrf*ht(i,j)  ! base geopotential at the lowest level 
                                  ! is terrain elevation * gravity
      ! Integerate base geopotential. this assures that the base state
      ! is in exact hydrostatic balance with respect to the model equations.
      ! this field is on full levels.
      DO k = 2, nz
        phb(i,j,k) = phb(i,j,k-1) - dnw(k-1)*mub(i,j)*alb(i,j,k-1)
      END DO
    END DO
  END DO

  ! Filling edges
  mub(nx,:)           = mub(nx-1,:)
  mu (nx,:)           = mu(nx-1,:)
  pb (nx,:,1:nz-1)    = pb(nx-1,:,1:nz-1)
  alb(nx,:,1:nz-1)    = alb(nx-1,:,1:nz-1)
  t_init(nx,:,1:nz-1) = t_init(nx-1,:,1:nz-1)
  phb(nx,:,:)         = pb(nx-1,:,:)

  mub(:,ny)           = mub(:,ny-1)
  mu (:,ny)           = mu(:,ny-1)
  pb (:,ny,1:nz-1)    = pb(:,ny-1,1:nz-1)
  alb(:,ny,1:nz-1)    = alb(:,ny-1,1:nz-1)
  t_init(:,ny,1:nz-1) = t_init(:,ny-1,1:nz-1)
  phb (:,ny,:)        = pb(:,ny-1,:)

  pt = pt_input - t0
  DO j = 1, ny-1
    DO i = 1, nx-1

      ! first get the pressure perturbation, moisture and inverse density
      ! at the top-most level

      qvf1 = qv(i,j,nz-1)
      qvf2 = 1./(1.+qvf1)
      qvf1 = qvf1*qvf2
      p(i,j,nz-1) = - 0.5*(mu(i,j) + qvf1*mub(i,j))*dnw(nz-1)/qvf2
      qvf  = 1.0 + rvovrd*qv(i,j,nz-1)
      alt(i,j,nz-1) = (r_d/p1000mb)*(pt(i,j,nz-1)+t0)*qvf*              &
                      (((p(i,j,nz-1)+pb(i,j,nz-1))/p1000mb)**cvpm)
      al(i,j,nz-1)  = alt(i,j,nz-1) - alb(i,j,nz-1)
    
     
      ! now, integrate down the column.
      DO k = nz-2,1,-1
        qvf1 = 0.5*(qv(i,j,k)+ qv(i,j,k+1))
        qvf2 = 1./(1.+qvf1)
        qvf1 = qvf1*qvf2
        p(i,j,k) = p(i,j,k+1) - (mu(i,j) + qvf1*mub(i,j))/qvf2/rdn(k+1)
        qvf  = 1.0 + rvovrd*qv(i,j,k)
        alt(i,j,k) = (r_d/p1000mb)*(pt(i,j,k)+t0)*qvf*                  &
                        (((p(i,j,k)+pb(i,j,k))/p1000mb)**cvpm)
        al(i,j,k)  = alt(i,j,k) - alb(i,j,k)
      END DO

      ph(i,j,1) = 0.0
      DO k = 2, nz
        ph(i,j,k) = ph(i,j,k-1) - dnw(k-1) *                            &
                    ( (mub(i,j)+mu(i,j))*al(i,j,k-1)+mu(i,j)*alb(i,j,k-1))
      END DO

    END DO
  END DO

  RETURN
END SUBROUTINE compute_ph
!
!##################################################################
!##################################################################
!######                                                      ######
!######            SUBROUTINE init_soil_depth                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE init_soil_depth(sfcphys, zs,dzs,nzsoil)

!------------------------------------------------------------------
!
!  Define layers (top layer = 0.01 m).  Double the thicknesses at 
!  each step (dzs values). The distance from the ground level to 
!  the midpoint of the layer is given by zs.
!
!    -------   Ground Level   ----------      ||      ||   ||  ||
!                                             ||      ||   ||  || zs(1) = 0.005 m
!    --  --  --  --  --  --  --  --  --       ||      ||   ||  \/
!                                             ||      ||   ||
!    -----------------------------------      ||  ||  ||   \/   dzs(1) = 0.01 m
!                                             ||  ||  ||
!                                             ||  ||  || zs(2) = 0.02
!    --  --  --  --  --  --  --  --  --       ||  ||  \/
!                                             ||  ||
!                                             ||  ||
!    -----------------------------------  ||  ||  \/   dzs(2) = 0.02 m
!                                         ||  ||
!                                         ||  ||
!                                         ||  ||
!                                         ||  || zs(3) = 0.05
!    --  --  --  --  --  --  --  --  --   ||  \/
!                                         ||
!                                         ||
!                                         ||
!                                         ||
!    -----------------------------------  \/   dzs(3) = 0.04 m
!
! See file share/module_soil_pre.F and subroutine init_soil_depth_?
!
!----------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: sfcphys
  REAL,    INTENT(OUT) :: zs(9), dzs(9)
  INTEGER, INTENT(OUT) :: nzsoil
  
!------------------------------------------------------------------
!
!  Misc. local variables
!
!------ -----------------------------------------------------------
   INTEGER :: k

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF(sfcphys == 1) THEN
    nzsoil = 5
    dzs(1) = 0.01
    zs(1)  = 0.5*dzs(1)
    DO k=2,nzsoil
       dzs(k)=2*dzs(k-1)
       zs(k)=zs(k-1)+.5*dzs(k-1)+.5*dzs(k)
    END DO

  ELSE IF(sfcphys == 2) THEN
    nzsoil = 4
    dzs(1:4) = (/ 0.1 , 0.3 , 0.6 , 1.0 /)
    zs(1)=.5*dzs(1)
    DO k=2,nzsoil
       zs(k)=zs(k-1)+.5*dzs(k-1)+.5*dzs(k)
    END DO

  ELSE IF(sfcphys == 3) THEN
    nzsoil = 6
    zs(1:6)  = (/ 0.00 , 0.05 , 0.20 , 0.40 , 1.60 , 3.00 /)
    dzs(1:6) = (/ 0.00 , 0.125, 0.175 , 0.70 , 1.30 , 1.40 /)

!    nzsoil = 9
!    zs  = (/ 0.00 , 0.05 , 0.20 , 0.40 , 0.60, 1.00, 1.60 , 2.20, 3.00 /)
!!   dzs = (/ 0.00 , 0.125, 0.175 , 0.70 , 1.30 , 1.40 /)

  ELSE
    WRITE(6,*) ' Wrong bl_surface_physics option', sfcphys
    WRITE(6,*) ' Must be 1, 2 or 3'
    STOP
  END IF

  RETURN
END SUBROUTINE init_soil_depth
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE vinterp_soil               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE vinterp_soil(nx,ny,nzsoil_in,zs_in,tsoil_in,qsoil_in,        &
                      nzsoil_out,zs_out,tsoil_out,qsoil_out)

!------------------------------------------------------------------
!
! PURPOSE:
!
!  Vertically interpolate ARPS soil variables to WRF soil layers.
!
!------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: nx,ny                     ! WRF horizontal grid size
  INTEGER, INTENT(IN)  :: nzsoil_in                 ! ARPS soil layer size
  REAL,    INTENT(IN)  :: zs_in(nx,ny,nzsoil_in)    ! ARPS soil layer depth
  REAL,    INTENT(IN)  :: tsoil_in(nx,ny,nzsoil_in)
  REAL,    INTENT(IN)  :: qsoil_in(nx,ny,nzsoil_in)
 
  INTEGER, INTENT(IN)  :: nzsoil_out                ! WRF soil layer size
  REAL,    INTENT(IN)  :: zs_out(nzsoil_out)        ! WRF soil layer depth
  REAL,    INTENT(OUT) :: tsoil_out(nx,ny,nzsoil_out)
  REAL,    INTENT(OUT) :: qsoil_out(nx,ny,nzsoil_out)

!------------------------------------------------------------------
!
!  Misc. local variables
!
!------------------------------------------------------------------

  INTEGER :: i,j,k,m
  REAL    :: w1,w2
  REAL    :: dist,wt,wq

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF(nzsoil_in == 2) THEN    ! suppose ARPS old soil model
                             ! contains two layers, 
                             ! soil layer -- 0-10 cm,  (0.05m)and 
                             ! deep soil  -- 10-100 cm (0.55m)
    DO k = 1, nzsoil_out
      IF(zs_out(k) < 0.05) THEN            ! extrapolation
        dist = 0.05 - zs_out(k)
        DO j = 1,ny
          DO i = 1,nx
            wt = (tsoil_in(i,j,1) - tsoil_in(i,j,2))/0.5
            tsoil_out(i,j,k) = tsoil_in(i,j,1) + wt*dist
            wq = (qsoil_in(i,j,1) - qsoil_in(i,j,2))/0.5
            qsoil_out(i,j,k) = qsoil_in(i,j,1) + wq*dist
          END DO
        END DO
      ELSE IF(zs_out(k) >= 0.55) THEN      ! constant layers
        tsoil_out(:,:,k) = tsoil_in(:,:,2)
        qsoil_out(:,:,k) = qsoil_in(:,:,2)
      ELSE
        w1 = ( 0.55      - zs_out(k) )/0.5        
        w2 = ( zs_out(k) - 0.05      )/0.5        
        DO j = 1,ny
          DO i = 1,nx
             tsoil_out(i,j,k) = tsoil_in(i,j,1)*w1 + tsoil_in(i,j,2)*w2
             qsoil_out(i,j,k) = qsoil_in(i,j,1)*w1 + qsoil_in(i,j,2)*w2
          END DO
        END DO
      END IF
    END DO

  ELSE IF(nzsoil_in > 2) THEN  ! ARPS new soil model

    DO j = 1,ny
      DO i = 1,nx
        DO k = 1, nzsoil_out
          IF(zs_out(k) < zs_in(i,j,1)) THEN            ! extrapolation
            dist = zs_in(i,j,1) - zs_out(k)
            wt = (tsoil_in(i,j,1) - tsoil_in(i,j,2))/(zs_in(i,j,2)-zs_in(i,j,1))
            wq = (qsoil_in(i,j,1) - qsoil_in(i,j,2))/(zs_in(i,j,2)-zs_in(i,j,1))
            tsoil_out(i,j,k) = tsoil_in(i,j,1) + wt*dist
            qsoil_out(i,j,k) = qsoil_in(i,j,1) + wq*dist
          ELSE IF(zs_out(k) >= zs_in(i,j,nzsoil_in)) THEN ! constant layers
            tsoil_out(i,j,k) = tsoil_in(i,j,nzsoil_in)
            qsoil_out(i,j,k) = qsoil_in(i,j,nzsoil_in)
          ELSE                         ! linear interpolation
            DO m = 1, nzsoil_in
              IF(zs_in(i,j,m) > zs_out(k)) EXIT
            END DO
            w1 = (zs_in(i,j,m) - zs_out(k))/(zs_in(i,j,m)-zs_in(i,j,m-1))
            w2 = (zs_out(k) - zs_in(i,j,m-1) )/(zs_in(i,j,m)-zs_in(i,j,m-1))
            tsoil_out(:,:,k) = tsoil_in(:,:,m-1)*w1 + tsoil_in(:,:,m)*w2
            qsoil_out(:,:,k) = qsoil_in(:,:,m-1)*w1 + qsoil_in(:,:,m)*w2
          END IF

        END DO   ! k
      END DO   ! i
    END DO  ! j
  ELSE
    WRITE(6,*) ' Wrong soil levels, nzsoil_arps = ',nzsoil_in
    STOP
  END IF

  RETURN
END SUBROUTINE vinterp_soil
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE get_sfcdt                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE get_sfcdt(sfcinitopt,sfcdtfl,ncompressx,ncompressy,nxsm,nysm,&
            nx,ny,nstyps,dx,dy,mapproj,trulat1,trulat2,trulon,sclfct,  &
            ctrlat,ctrlon,soiltyp,stypfrct,vegtyp,veg,tem1,            &
            hinterp_needed,nx_wrf,ny_wrf,dx_wrf,dy_wrf,                &
            mapproj_wrf,trulat1_wrf,trulat2_wrf,trulon_wrf,jday,       &
            x2d,y2d,xs,ys,iloc,jloc,                                   &
            hgt_wrf,soiltyp_wrf,vegtyp_wrf,vegfrct_wrf,xice,xland,     &
            tmn_wrf,shdmin,shdmax,albbck,snoalb,istatus)

!-----------------------------------------------------------------------
!
! PURPOSE:
!   Read in surface characteristic data to intialize
!
!     ARPS         WRF              ARPS          WRF
!  ===========   ============    ============   ============== 
!   soiltyp      soiltyp_wrf      stypfrct       stypfrct_wrf
!   vegtyp       vegtyp_wrf       veg            veg_wrf
!                vegfrct_wrf
!                xice                            xland
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: sfcinitopt
  CHARACTER(LEN=*), INTENT(IN) :: sfcdtfl
  INTEGER,          INTENT(IN) :: ncompressx,ncompressy
  INTEGER,          INTENT(IN) :: nxsm, nysm
  INTEGER,          INTENT(IN) :: nx,   ny,  nstyps, nx_wrf, ny_wrf
  REAL,             INTENT(IN) :: dx,   dy,          dx_wrf, dy_wrf
  INTEGER,          INTENT(IN) :: mapproj,           mapproj_wrf
  REAL,             INTENT(IN) :: trulat1,trulat2, trulat1_wrf,trulat2_wrf
  REAL,             INTENT(IN) :: trulon,          trulon_wrf
  REAL,             INTENT(IN) :: ctrlat,ctrlon
  INTEGER,          INTENT(IN) :: sclfct, jday
  LOGICAL,          INTENT(IN) :: hinterp_needed
  REAL,             INTENT(IN) :: x2d(nx_wrf,ny_wrf),  y2d(nx_wrf,ny_wrf)
  INTEGER,          INTENT(IN) :: iloc(nx_wrf,ny_wrf), jloc(nx_wrf,ny_wrf)
  REAL,             INTENT(IN) :: xs(nx), ys(ny)
  INTEGER                      :: soiltyp (nx,ny,nstyps)  ! ARPS arrays
  REAL                         :: stypfrct(nx,ny,nstyps)
  INTEGER                      :: vegtyp(nx,ny)
  REAL                         :: veg(nx,ny)
  REAL                         :: tem1(nx,ny,3)           ! ARPS work arrays
  REAL,            INTENT(OUT) :: hgt_wrf(nx_wrf,ny_wrf)  ! WRF arrays
  INTEGER,         INTENT(OUT) :: soiltyp_wrf(nx_wrf,ny_wrf)
  INTEGER,         INTENT(OUT) :: vegtyp_wrf(nx_wrf,ny_wrf)
  REAL,            INTENT(OUT) :: vegfrct_wrf(nx_wrf,ny_wrf)
  REAL,            INTENT(OUT) :: xice(nx_wrf,ny_wrf)
  REAL,            INTENT(OUT) :: xland(nx_wrf,ny_wrf)
  REAL,            INTENT(OUT) :: tmn_wrf(nx_wrf,ny_wrf)
  REAL,            INTENT(OUT) :: shdmax(nx_wrf,ny_wrf)
  REAL,            INTENT(OUT) :: shdmin(nx_wrf,ny_wrf)
  REAL,            INTENT(OUT) :: albbck(nx_wrf,ny_wrf)
  REAL,            INTENT(OUT) :: snoalb(nx_wrf,ny_wrf)
  INTEGER,         INTENT(OUT) :: istatus

  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
! Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER, ALLOCATABLE :: soiltypsm (:,:,:)
  REAL,    ALLOCATABLE :: stypfrctsm(:,:,:)
  INTEGER, ALLOCATABLE :: vegtypsm  (:,:)
  REAL,    ALLOCATABLE :: vegsm     (:,:)

  INTEGER :: i, j, ii, jj, ia, ja
  INTEGER :: lenfil

  CHARACTER(LEN=256) :: tmpstr

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  lenfil = LEN_TRIM(sfcdtfl)
  WRITE(tmpstr,'(a)') TRIM(sfcdtfl)

  IF(sfcinitopt(1:4) == 'ARPS') THEN
    
    DO jj = 1, ncompressy
      DO ii = 1, ncompressx
        IF (mp_opt > 0 .AND. readsplit(FINDX_T) <= 0) THEN
          CALL gtsplitfn(sfcdtfl,ncompressx,ncompressy,loc_x,loc_y,ii,jj, &
                         0,0,1,2,tmpstr,istatus)
          lenfil = LEN_TRIM(tmpstr)
        END IF

        IF (mp_opt > 0 .AND. readsplit(FINDX_T) > 0) THEN
          CALL readsplitsfcdt(nx,ny,nstyps,TRIM(tmpstr),dx,dy,          &
                   mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon, &
                   soiltyp,stypfrct,vegtyp,tem1(:,:,1),tem1(:,:,2),     &
                   veg, tem1(:,:,3))
        ELSE

          IF (.NOT. ALLOCATED(soiltypsm)) THEN
            ALLOCATE(soiltypsm (nxsm,nysm,nstyps), STAT = istatus)
            ALLOCATE(stypfrctsm(nxsm,nysm,nstyps), STAT = istatus)
            ALLOCATE(vegtypsm  (nxsm,nysm),        STAT = istatus)
            ALLOCATE(vegsm     (nxsm,nysm),        STAT = istatus)
          END IF

          CALL readsfcdt(nxsm,nysm,nstyps,TRIM(tmpstr),dx,dy,           &
                 mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,   &
                 soiltypsm,stypfrctsm,vegtypsm,tem1(:,:,1),tem1(:,:,2), &
                 vegsm, tem1(:,:,3))

          DO j = 1, nysm
            ja = (jj-1)*(nysm-3)+j
            DO i = 1, nxsm
              ia = (ii-1)*(nxsm-3)+i
              soiltyp(ia,ja,:)  = soiltypsm(i,j,:)
              stypfrct(ia,ja,:) = stypfrctsm(i,j,:)
              vegtyp(ia,ja)     = vegtypsm(i,j)
              veg(ia,ja)        = vegsm(i,j)
            END DO
          END DO

        END IF
      END DO
    END DO

    IF (ALLOCATED(soiltypsm))  DEALLOCATE(soiltypsm)
    IF (ALLOCATED(stypfrctsm)) DEALLOCATE(stypfrctsm)
    IF (ALLOCATED(vegtypsm))   DEALLOCATE(vegtypsm)
    IF (ALLOCATED(vegsm))      DEALLOCATE(vegsm)

    CALL sfcinterp(nx_wrf,ny_wrf,nx,ny,nstyps,hinterp_needed,           &
                   x2d,y2d,xs,ys,iloc,jloc,                             &
                   soiltyp,stypfrct,vegtyp,veg,                         &
                   soiltyp_wrf,vegtyp_wrf,vegfrct_wrf,xice,xland)

    vegfrct_wrf(:,:) = vegfrct_wrf(:,:)*100.

    IF (myproc == 0) THEN
      WRITE(6,'(/,2(a,/),3(10x,a/),a,/)')                               &
        ' ************************** WARNING ************************', &             
  ' WARNING: Noah land-surface model, sf_surface_physics = 2, requires',&
  'annual MAX/MIN veg fraction, MAX snow albedo, background albedo.',   &
  'However, ARPSSFC does not provide such variables. You had better ',  &
  'use sfcinitopt = "GEOGRID", "WRFSI" or "WRF". see arps2wrf.input.',  &
        ' ************************** WARNING ************************'
    END IF
    hgt_wrf(:,:) = 0.0
    tmn_wrf(:,:) = -999.0 ! indicates it is not available
    shdmax(:,:)  = 0.0
    shdmin(:,:)  = 0.0
    snoalb(:,:)  = 0.0
    albbck(:,:)  = 0.0

    istatus = 1           ! change to 0 if you want to read ARPSSFC data

  ELSE IF(sfcinitopt(1:5) == 'WRFSI') THEN

    CALL get_sfcdt_wrfsi(TRIM(tmpstr),nx_wrf,ny_wrf,dx_wrf,dy_wrf,      &
                   mapproj_wrf,trulat1_wrf,trulat2_wrf,trulon_wrf,jday, &
                   hgt_wrf,soiltyp_wrf,vegtyp_wrf,vegfrct_wrf,xice,     &
                   xland,tmn_wrf,shdmin,shdmax,albbck,snoalb,istatus)

  ELSE IF (sfcinitopt(1:3) == 'WRF') THEN

    CALL get_sfcdt_wrf(TRIM(tmpstr),nx_wrf,ny_wrf,dx_wrf,dy_wrf,        &
                   mapproj_wrf,trulat1_wrf,trulat2_wrf,trulon_wrf,jday, &
                   hgt_wrf,soiltyp_wrf,vegtyp_wrf,vegfrct_wrf,xice,     &
                   xland,tmn_wrf,shdmin,shdmax,albbck,snoalb,istatus)

  ELSE IF (sfcinitopt(1:7) == 'GEOGRID') THEN

    CALL get_sfcdt_geogrid(TRIM(tmpstr),nx_wrf,ny_wrf,dx_wrf,dy_wrf,    &
                   mapproj_wrf,trulat1_wrf,trulat2_wrf,trulon_wrf,jday, &
                   hgt_wrf,soiltyp_wrf,vegtyp_wrf,vegfrct_wrf,xice,     &
                   xland,tmn_wrf,shdmin,shdmax,albbck,snoalb,istatus)

  ELSE
    WRITE(6,*) 'Unsupported surface data format: ',sfcinitopt
    istatus = -1
  END IF

  RETURN
END SUBROUTINE get_sfcdt
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE sfcinterp                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE sfcinterp(nx,ny,nx_ext,ny_ext,nstyps_ext,interp,             &
                    x2d,y2d,x_ext,y_ext,iloc,jloc,                      &
                    soiltyp_ext,stypfrct_ext,vegtyp_ext,vegfrct_ext,    &
                    soiltyp,vegtyp,vegfrct,xice,xland)

!------------------------------------------------------------------
!
! Purpose:
!   Interpolate ARPS surface variables to WRF grid. Those variables
!   include soil types, vegtation types, vegatation fraction
!   ice flage and land flage etc.
!
!------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: nx,ny          ! WRF grid size
  INTEGER, INTENT(IN)  :: nx_ext,ny_ext  ! ARPS grid size
  INTEGER, INTENT(IN)  :: nstyps_ext     ! ARPS soil type size

  LOGICAL, INTENT(IN)  :: interp         ! whether interpolation is needed

  REAL,    INTENT(IN)  :: x2d(nx,ny)
  REAL,    INTENT(IN)  :: y2d(nx,ny)
  INTEGER, INTENT(IN)  :: iloc(nx,ny)
  INTEGER, INTENT(IN)  :: jloc(nx,ny)
  REAL,    INTENT(IN)  :: x_ext(nx_ext)  ! ARPS grid coordinate
  REAL,    INTENT(IN)  :: y_ext(ny_ext)

  INTEGER, INTENT(IN)  :: soiltyp_ext(nx_ext,ny_ext,nstyps_ext)
  INTEGER, INTENT(IN)  :: vegtyp_ext(nx_ext,ny_ext)
  REAL,    INTENT(IN)  :: stypfrct_ext(nx_ext,ny_ext,nstyps_ext)
  REAL,    INTENT(IN)  :: vegfrct_ext(nx_ext,ny_ext)

  INTEGER, INTENT(OUT) :: soiltyp(nx,ny)
  INTEGER, INTENT(OUT) :: vegtyp(nx,ny)
  REAL,    INTENT(OUT) :: vegfrct(nx,ny)

  REAL,    INTENT(OUT) :: xice(nx,ny)
  REAL,    INTENT(OUT) :: xland(nx,ny)

!------------------------------------------------------------------
!
!  Misc. local variables
!
!------------------------------------------------------------------
  INTEGER :: i,j
  INTEGER :: ii,jj
  INTEGER :: k
  REAL    :: dx,dy
  REAL    :: dxmid,dymid
  INTEGER :: itmp(1)

  INTEGER, PARAMETER :: veg_table(14)  = (/19,22, 7, 9, 6,11,14,        &
                                           13,24, 3,17, 8,19,16  /)
!
!  ARPS                                     WRF
!  ====                                    ====
!     1  ! Desert                            19      ! Bar. Sparse Veg.
!     2  ! Tundra                            22      ! Mixed Tundra
!     3  ! Grassland                          7      ! Grassland
!     4  ! Grassland with shrub cover         9      ! Mix Shrb./Grs.
!     5  ! Grassland with tree cover          6      ! Crop./Wood Mosc
!     6  ! Deciduous forest                  11      ! Decids. Broadlf.
!     7  ! Evergreen forest                  14      ! Evergrn. Needlf.
!     8  ! Rain forest                       13      ! Evergn Broadlf Forest
!     9  ! Ice                               24      ! Snow or Ice
!    10  ! Cultivation                        3      ! Irrg. Crop. Past.
!    11  ! Bog or marsh                      17      ! Herb. Wetland
!    12  ! Dwarf shrub                        8      ! Shrubland
!    13  ! Semidesert                        19      ! Bar. Sparse Veg.
!    14  ! Water                             16      ! Water bodies
!

  INTEGER, PARAMETER :: soil_table(13) = (/ 1, 2, 3, 4, 6, 7, 8,        &
                                            9,10,11,12,14,14  /)
!
!  ARPS    WRF                       ||   ARPS     WRF
! =====   ====                       ||  =====    ====
!     1      1  ! Sand               ||      8       9  ! Clay loam
!     2      2  ! Loamy sand         ||      9      10  ! Sandy clay
!     3      3  ! Sandy loam         ||     10      11  ! Silty clay
!     4      4  ! Silt loam          ||     11      12  ! Clay
!     5      6  ! Loam               ||     12      16  ! Ice     
!     6      7  ! Sandy clay loam    ||                 ! Changed to 14 by Jerry, 10/20/2003
!     7      8  ! Silty clay loam    ||     13      14  ! Water
!

  INTEGER, PARAMETER :: ice   = 24            ! vegetation type
  INTEGER, PARAMETER :: water = 16            ! vegetation type

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  dxmid = (x_ext(3)-x_ext(2))/2.0
  dymid = (y_ext(3)-y_ext(2))/2.0

!
!(ii,jj+1)                (ii+1,jj+1)
!       _____________________
!       |         |          |
!       |         |          |
!       |         |          |
!       |---------|----------|
!       |         |          |
! dymid{|  dxmid  |          |
!       |_________|__________|
!      (ii,jj)             (ii+1,jj)
!

! nearest neighbor interpolation

  DO j = 1,ny
    DO i = 1,nx
      IF(interp) THEN
        ii = iloc(i,j)
        jj = jloc(i,j)
        dx = x2d(i,j) - x_ext(ii)
        dy = y2d(i,j) - y_ext(jj)
        IF(dx > dxmid) ii = ii + 1
        IF(dy > dymid) jj = jj + 1
      ELSE
        ii = i+1
        jj = j+1
      END IF

      vegtyp(i,j)  = veg_table(vegtyp_ext(ii,jj))
      vegfrct(i,j) = vegfrct_ext(ii,jj)
      itmp = MAXLOC(stypfrct_ext(ii,jj,:))
      k = itmp(1)
      soiltyp(i,j) = soil_table(soiltyp_ext(ii,jj,k))

      xice(i,j) = 0.0
      IF(vegtyp(i,j) == ice) xice(i,j) = 1.

      xland(i,j) = 1.0
      IF(vegtyp(i,j) == water) THEN
        xland(i,j)   = 2.0
        vegfrct(i,j) = 0
      END IF

      IF(soiltyp(i,j) == 14 .AND. xland(i,j) < 1.5) THEN
         soiltyp(i,j) = 7
      ELSE IF(xland(i,j) > 1.5) THEN
         soiltyp(i,j) = 14
      END IF

!      IF(vegtyp(i,j) == 0)   &
!        WRITE(6,*) 'Vegtyp is zero at i = ',i,' j = ',j
!      IF(soiltyp(i,j) == 0)  &
!        WRITE(6,*) 'Soiltyp is zero at i = ',i,' j = ',j

    END DO
  END DO
  
  RETURN
END SUBROUTINE sfcinterp
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE get_sfcdt_wrfsi            ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE get_sfcdt_wrfsi(filnam,nx,ny,dx,dy,                          &
                     mapproj,trulat1,trulat2,trulon,jday,               &
                     hgt,soiltyp,vegtyp,vegfrct,xice,                   &
                     xland,tmn,shdmin,shdmax,albbck,snoalb,istatus)

!------------------------------------------------------------------
!
! PURPOSE:
!   Read surface variables from NetCDF static file created by
!   gridgen_model of WRFSI package
!   
! NOTE:
!   Those data must be at the same WRF grid as that specified
!   by the arguments.
!
!------------------------------------------------------------------

  USE wrf_metadata

  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN) :: filnam           ! static file name
  INTEGER, INTENT(IN)  :: nx,ny                    ! WRF grid size
  INTEGER, INTENT(IN)  :: mapproj                  ! Map projection
  REAL,    INTENT(IN)  :: dx,dy                    ! WRF grid length
  REAL,    INTENT(IN)  :: trulat1,trulat2,trulon   ! Map projection specification
  INTEGER, INTENT(IN)  :: jday                     ! julian day of the year

  REAL,    INTENT(OUT) :: hgt    (nx,ny)
  INTEGER, INTENT(OUT) :: soiltyp(nx,ny)
  INTEGER, INTENT(OUT) :: vegtyp (nx,ny)
  REAL,    INTENT(OUT) :: vegfrct(nx,ny)
  REAL,    INTENT(OUT) :: xice   (nx,ny)
  REAL,    INTENT(OUT) :: xland  (nx,ny)
  REAL,    INTENT(OUT) :: tmn    (nx,ny)
  REAL,    INTENT(OUT) :: shdmin (nx,ny)
  REAL,    INTENT(OUT) :: shdmax (nx,ny)
  REAL,    INTENT(OUT) :: albbck (nx,ny)
  REAL,    INTENT(OUT) :: snoalb (nx,ny)
  INTEGER, INTENT(OUT) :: istatus

!------------------------------------------------------------------
!
!  Misc. local variables
!
!------------------------------------------------------------------

  INTEGER :: nscid

  REAL, ALLOCATABLE :: landusef(:,:,:)
  REAL, ALLOCATABLE :: soiltyp_top(:,:,:),soiltyp_bot(:,:,:)

  INTEGER :: i,j,k

  INTEGER :: nxlg, nylg

  REAL, ALLOCATABLE :: var3d1(:,:,:), var3d2(:,:,:), var3d3(:,:,:)
  REAL, ALLOCATABLE :: var2d(:,:)

  REAL              :: rdummy(1,1,12)

  INCLUDE   'mp.inc'

!-----------------------------------------------------------------

  INTERFACE

    SUBROUTINE get_static_landusef(ncid,var3d)
      INTEGER,          INTENT(IN)  :: ncid
      REAL,             INTENT(OUT) :: var3d(:,:,:)
    END SUBROUTINE get_static_landusef

    SUBROUTINE get_static_soil(ncid,var3d1,var3d2)
      INTEGER,          INTENT(IN)  :: ncid
      REAL,             INTENT(OUT) :: var3d1(:,:,:)
      REAL,             INTENT(OUT) :: var3d2(:,:,:)
    END SUBROUTINE get_static_soil

    SUBROUTINE get_static_monthly(staticopt,ncid, vartype, jday,        &
                                  var2d, varin3d, istatus)
      INTEGER,          INTENT(IN)  :: staticopt
      INTEGER,          INTENT(IN)  :: ncid
      CHARACTER(LEN=1), INTENT(IN)  :: vartype
      INTEGER,          INTENT(IN)  :: jday
      REAL,             INTENT(OUT) :: var2d(:,:)
      REAL,             INTENT(IN)  :: varin3d(:,:,:)
      INTEGER,          INTENT(OUT) :: istatus
    END SUBROUTINE get_static_monthly

    SUBROUTINE get_static_2d(ncid,varname,var2d)
      INTEGER,          INTENT(IN)  :: ncid
      CHARACTER(LEN=3), INTENT(IN)  :: varname
      REAL,             INTENT(OUT) :: var2d(:,:)
    END SUBROUTINE get_static_2d

  END INTERFACE

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  nxlg = (nx - 1)*nproc_x + 1
  nylg = (ny - 1)*nproc_y + 1

!  write(0,*) 'nx = ',nx, 'ny = ',ny, 'nxlg = ',nxlg, 'nylg = ',nylg

  ALLOCATE(landusef(nx,ny,LanduseCategories), STAT = istatus)
  ALLOCATE(soiltyp_top(nx,ny,SoilCategories), STAT = istatus)
  ALLOCATE(soiltyp_bot(nx,ny,SoilCategories), STAT = istatus)

  ALLOCATE(var3d1(nxlg,nylg,LanduseCategories), STAT = istatus)
  ALLOCATE(var3d2(nxlg,nylg,SoilCategories),    STAT = istatus)
  ALLOCATE(var3d3(nxlg,nylg,SoilCategories),    STAT = istatus)
  ALLOCATE(var2d (nxlg,nylg),                   STAT = istatus)
!
! Open the static NetCDF file which is created by gridgen_model.
!

  IF(myproc == 0) CALL open_static_file(filnam,nscid)

! we need to check whether the static file is
! on the same grid as WRF grid

  IF(myproc == 0) THEN 
    CALL check_static_grid(0,nscid,nx,ny,dx,dy,mapproj,                   &
                         trulat1,trulat2,trulon,istatus)
    IF(istatus /= 0) RETURN
  END IF

  IF (myproc == 0) THEN
    CALL get_static_landusef(nscid,var3d1)
    CALL get_static_soil(nscid,var3d2,var3d3)
  END IF
  CALL wrf_split3d(var3d1,nx,ny,LanduseCategories,1,landusef)
  CALL wrf_split3d(var3d2,nx,ny,SoilCategories,1,soiltyp_top)
  CALL wrf_split3d(var3d3,nx,ny,SoilCategories,1,soiltyp_bot)

  IF (myproc == 0) CALL get_static_monthly (0,nscid, 'g', jday, var2d, rdummy,istatus)
  IF (istatus /= 0) CALL mpexit(1)
  CALL wrf_split2d(var2d,nx,ny,1,vegfrct)

  IF (myproc == 0) CALL get_static_monthly (0,nscid, 'a', jday, var2d, rdummy,istatus)
  IF (istatus /= 0) CALL mpexit(1)
  CALL wrf_split2d(var2d,nx,ny,1,albbck)
  
  IF (myproc == 0) CALL get_static_2d(nscid,'gnn',var2d)
  CALL wrf_split2d(var2d,nx,ny,1,shdmin)

  IF (myproc == 0) CALL get_static_2d(nscid,'gnx',var2d)
  CALL wrf_split2d(var2d,nx,ny,1,shdmax)

  IF (myproc == 0) CALL get_static_2d(nscid,'alb',var2d)
  CALL wrf_split2d(var2d,nx,ny,1,snoalb)

  IF (myproc == 0) CALL get_static_2d(nscid,'avc',var2d)
  CALL wrf_split2d(var2d,nx,ny,1,hgt)

  IF (myproc == 0) CALL get_static_2d(nscid,'tmp',var2d)
  CALL wrf_split2d(var2d,nx,ny,1,tmn)

  IF(myproc == 0) CALL close_static_file(nscid)

  DEALLOCATE(var3d1,var3d2,var3d3)
  DEALLOCATE(var2d)

!--------------------------------------------------------------------
!
! Vegetation index and soil texture indexSoilCategories
!
!--------------------------------------------------------------------

  CALL land_soil_cat(nx,ny,LanduseCategories,SoilCategories,            &
                     landusef,soiltyp_top,soiltyp_bot,                  &
                     vegtyp,soiltyp,xland,xice,istatus)

  DEALLOCATE(landusef, soiltyp_top, soiltyp_bot)

  RETURN
END SUBROUTINE get_sfcdt_wrfsi
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE get_sfcdt_wrf              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE get_sfcdt_wrf(filnam,nx,ny,dx,dy,                            &
                     mapproj,trulat1,trulat2,trulon,jday,               &
                     hgt,soiltyp,vegtyp,vegfrct,xice,                   &
                     xland,tmn,shdmin,shdmax,albbck,snoalb,istatus)

!------------------------------------------------------------------
!
! PURPOSE:
!   Read surface variables from NetCDF static file created by
!   staticpost of WRFSI package or wrfstatic from ARPS package
!   
! NOTE:
!   Those data must be at the same WRF grid as that specified
!   by the arguments.
!
!------------------------------------------------------------------

  USE wrf_metadata

  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN) :: filnam           ! static file name
  INTEGER, INTENT(IN)  :: nx,ny                    ! WRF grid size
  INTEGER, INTENT(IN)  :: mapproj                  ! Map projection
  REAL,    INTENT(IN)  :: dx,dy                    ! WRF grid length
  REAL,    INTENT(IN)  :: trulat1,trulat2,trulon   ! Map projection specification
  INTEGER, INTENT(IN)  :: jday                     ! julian day of the year

  REAL,    INTENT(OUT) :: hgt    (nx,ny)
  INTEGER, INTENT(OUT) :: soiltyp(nx,ny)
  INTEGER, INTENT(OUT) :: vegtyp (nx,ny)
  REAL,    INTENT(OUT) :: vegfrct(nx,ny)
  REAL,    INTENT(OUT) :: xice   (nx,ny)
  REAL,    INTENT(OUT) :: xland  (nx,ny)
  REAL,    INTENT(OUT) :: tmn    (nx,ny)
  REAL,    INTENT(OUT) :: shdmin (nx,ny)
  REAL,    INTENT(OUT) :: shdmax (nx,ny)
  REAL,    INTENT(OUT) :: albbck (nx,ny)
  REAL,    INTENT(OUT) :: snoalb (nx,ny)
  INTEGER, INTENT(OUT) :: istatus

!------------------------------------------------------------------
!
!  Misc. local variables
!
!------------------------------------------------------------------

  INTEGER :: nscid

  REAL, ALLOCATABLE :: landusef(:,:,:)
  REAL, ALLOCATABLE :: soiltyp_top(:,:,:),soiltyp_bot(:,:,:)

  INTEGER :: i,j,k

  INTEGER :: nxlg, nylg

  REAL, ALLOCATABLE :: temdm(:,:,:)
  REAL, ALLOCATABLE :: var3d1(:,:,:), var3d2(:,:,:), var3d3(:,:,:)
  REAL, ALLOCATABLE :: var2d(:,:)

  INCLUDE   'mp.inc'

!-----------------------------------------------------------------

  INTERFACE 

    SUBROUTINE get_static_monthly(staticopt,ncid, vartype, jday,        &
                                  var2d, varin3d, istatus)
      INTEGER,          INTENT(IN)  :: staticopt
      INTEGER,          INTENT(IN)  :: ncid
      CHARACTER(LEN=1), INTENT(IN)  :: vartype
      INTEGER,          INTENT(IN)  :: jday
      REAL,             INTENT(OUT) :: var2d(:,:)
      REAL,             INTENT(IN)  :: varin3d(:,:,:)
      INTEGER,          INTENT(OUT) :: istatus
    END SUBROUTINE get_static_monthly

  END INTERFACE

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  nxlg = (nx - 1)*nproc_x + 1
  nylg = (ny - 1)*nproc_y + 1

!  write(0,*) 'nx = ',nx, 'ny = ',ny, 'nxlg = ',nxlg, 'nylg = ',nylg

  ALLOCATE(landusef(nx,ny,LanduseCategories), STAT = istatus)
  ALLOCATE(soiltyp_top(nx,ny,SoilCategories), STAT = istatus)
  ALLOCATE(soiltyp_bot(nx,ny,SoilCategories), STAT = istatus)

  ALLOCATE(var3d1(nxlg,nylg,LanduseCategories), STAT = istatus)
  ALLOCATE(var3d2(nxlg,nylg,SoilCategories),    STAT = istatus)
  ALLOCATE(var3d3(nxlg,nylg,SoilCategories),    STAT = istatus)
  ALLOCATE(var2d (nxlg,nylg),                   STAT = istatus)

  ALLOCATE(temdm(nxlg-1,nylg-1,MAX(LanduseCategories,SoilCategories)), STAT = istatus)
!
! Open the static NetCDF file which is created by wrfstatic.
!

  IF(myproc == 0) CALL open_static_file(filnam,nscid)

! we need to check whether the static file is
! on the same grid as WRF grid

  IF(myproc == 0) THEN 
    CALL check_static_grid(1,nscid,nx,ny,dx,dy,mapproj,                 &
                         trulat1,trulat2,trulon,istatus)
    IF(istatus /= 0) RETURN
  END IF

  IF (myproc == 0) THEN
    CALL get_ncd_3d(nscid,1,'LANDUSEF',nxlg-1,nylg-1,LanduseCategories,temdm,istatus)
    DO k = 1, LanduseCategories
      DO j = 1,nylg-1
        DO i = 1,nxlg-1
          var3d1(i,j,k) = temdm(i,j,k)
        END DO
      END DO
    END DO

    CALL get_ncd_3d(nscid,1,'SOILCTOP',nxlg-1,nylg-1,SoilCategories,temdm,istatus)
    DO k = 1, SoilCategories
      DO j = 1,nylg-1
        DO i = 1,nxlg-1
          var3d2(i,j,k) = temdm(i,j,k)
        END DO
      END DO
    END DO

    CALL get_ncd_3d(nscid,1,'SOILCBOT',nxlg-1,nylg-1,SoilCategories,temdm,istatus)
    DO k = 1, SoilCategories
      DO j = 1,nylg-1
        DO i = 1,nxlg-1
          var3d3(i,j,k) = temdm(i,j,k)
        END DO
      END DO
    END DO

    CALL edgfill(var3d1,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,LanduseCategories,1,LanduseCategories)
    CALL edgfill(var3d2,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,SoilCategories,1,SoilCategories)
    CALL edgfill(var3d3,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,SoilCategories,1,SoilCategories)
  END IF
  CALL wrf_split3d(var3d1,nx,ny,LanduseCategories,1,landusef)
  CALL wrf_split3d(var3d2,nx,ny,SoilCategories,1,soiltyp_top)
  CALL wrf_split3d(var3d3,nx,ny,SoilCategories,1,soiltyp_bot)

  IF (myproc == 0) THEN
    CALL get_ncd_3d(nscid,1,'GREEN12M',nxlg-1,nylg-1,12,temdm,istatus)
    DO k = 1, 12
      DO j = 1,nylg-1
        DO i = 1,nxlg-1
          var3d1(i,j,k) = temdm(i,j,k)
        END DO
      END DO
    END DO
    CALL edgfill(var3d1,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,12,1,12)

    CALL get_static_monthly (1,nscid, 'g', jday, var2d, var3d1,istatus)
  END IF
  IF (istatus /= 0) CALL mpexit(1)
  CALL wrf_split2d(var2d,nx,ny,1,vegfrct)

  IF (myproc == 0) THEN
    CALL get_ncd_3d(nscid,1,'ALBDO12M',nxlg-1,nylg-1,12,temdm,istatus)
    DO k = 1, 12
      DO j = 1,nylg-1
        DO i = 1,nxlg-1
          var3d2(i,j,k) = temdm(i,j,k)
        END DO
      END DO
    END DO
    CALL edgfill(var3d2,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,12,1,12)

    CALL get_static_monthly (1,nscid, 'a', jday, var2d, var3d2,istatus)
  END IF
  IF (istatus /= 0) CALL mpexit(1)
  CALL wrf_split2d(var2d,nx,ny,1,albbck)
  
  IF (myproc == 0) THEN
    CALL get_ncd_2d(nscid,1,'SHDMIN',nxlg-1,nylg-1,temdm,istatus)
    DO j = 1,nylg-1
      DO i = 1,nxlg-1
        var2d(i,j) = temdm(i,j,1)
      END DO
    END DO
    CALL edgfill(var2d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,1,1,1)
  END IF
  CALL wrf_split2d(var2d,nx,ny,1,shdmin)

  IF (myproc == 0) THEN
    CALL get_ncd_2d(nscid,1,'SHDMAX',nxlg-1,nylg-1,temdm,istatus)
    DO j = 1,nylg-1
      DO i = 1,nxlg-1
        var2d(i,j) = temdm(i,j,1)
      END DO
    END DO
    CALL edgfill(var2d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,1,1,1)
  END IF
  CALL wrf_split2d(var2d,nx,ny,1,shdmax)

  IF (myproc == 0) THEN
    CALL get_ncd_2d(nscid,1,'SNOALB',nxlg-1,nylg-1,temdm,istatus)
    DO j = 1,nylg-1
      DO i = 1,nxlg-1
        var2d(i,j) = temdm(i,j,1)
      END DO
    END DO
    CALL edgfill(var2d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,1,1,1)
  END IF
  CALL wrf_split2d(var2d,nx,ny,1,snoalb)

  IF (myproc == 0) THEN
    CALL get_ncd_2d(nscid,1,'HGT',nxlg-1,nylg-1,temdm,istatus)
    DO j = 1,nylg-1
      DO i = 1,nxlg-1
        var2d(i,j) = temdm(i,j,1)
      END DO
    END DO
    CALL edgfill(var2d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,1,1,1)
  END IF
  CALL wrf_split2d(var2d,nx,ny,1,hgt)

  IF (myproc == 0) THEN
    CALL get_ncd_2d(nscid,1,'TMN',nxlg-1,nylg-1,temdm,istatus)
    DO j = 1,nylg-1
      DO i = 1,nxlg-1
        var2d(i,j) = temdm(i,j,1)
      END DO
    END DO
    CALL edgfill(var2d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,1,1,1)
  END IF
  CALL wrf_split2d(var2d,nx,ny,1,tmn)

  IF(myproc == 0) CALL close_static_file(nscid)

  DEALLOCATE(temdm)
  DEALLOCATE(var3d1,var3d2,var3d3)
  DEALLOCATE(var2d)

  CALL land_soil_cat(nx,ny,LanduseCategories,SoilCategories,            &
                     landusef,soiltyp_top,soiltyp_bot,                  &
                     vegtyp,soiltyp,xland,xice,istatus)

  DEALLOCATE(landusef, soiltyp_top, soiltyp_bot)

  RETURN
END SUBROUTINE get_sfcdt_wrf
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE get_sfcdt_geogrid          ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE get_sfcdt_geogrid(filnam,nx,ny,dx,dy,                        &
                     mapproj,trulat1,trulat2,trulon,jday,               &
                     hgt,soiltyp,vegtyp,vegfrct,xice,                   &
                     xland,tmn,shdmin,shdmax,albbck,snoalb,istatus)

!------------------------------------------------------------------
!
! PURPOSE:
!   Read surface variables from NetCDF static file created by
!   program geogrid from the WPS package.
!   
! NOTE:
!   Those data must be at the same WRF grid as that specified
!   by the arguments.
!
!------------------------------------------------------------------

  USE wrf_metadata

  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN) :: filnam           ! static file name
  INTEGER, INTENT(IN)  :: nx,ny                    ! WRF grid size
  INTEGER, INTENT(IN)  :: mapproj                  ! Map projection
  REAL,    INTENT(IN)  :: dx,dy                    ! WRF grid length
  REAL,    INTENT(IN)  :: trulat1,trulat2,trulon   ! Map projection specification
  INTEGER, INTENT(IN)  :: jday                     ! julian day of the year

  REAL,    INTENT(OUT) :: hgt    (nx,ny)
  INTEGER, INTENT(OUT) :: soiltyp(nx,ny)
  INTEGER, INTENT(OUT) :: vegtyp (nx,ny)
  REAL,    INTENT(OUT) :: vegfrct(nx,ny)
  REAL,    INTENT(OUT) :: xice   (nx,ny)
  REAL,    INTENT(OUT) :: xland  (nx,ny)
  REAL,    INTENT(OUT) :: tmn    (nx,ny)
  REAL,    INTENT(OUT) :: shdmin (nx,ny)
  REAL,    INTENT(OUT) :: shdmax (nx,ny)
  REAL,    INTENT(OUT) :: albbck (nx,ny)
  REAL,    INTENT(OUT) :: snoalb (nx,ny)
  INTEGER, INTENT(OUT) :: istatus

!------------------------------------------------------------------
!
!  Misc. local variables
!
!------------------------------------------------------------------

  INTEGER :: nscid

  REAL, ALLOCATABLE :: landusef(:,:,:)
  REAL, ALLOCATABLE :: soiltyp_top(:,:,:),soiltyp_bot(:,:,:)

  INTEGER :: i,j,k

  INTEGER :: nxlg, nylg

  REAL, ALLOCATABLE :: temdm(:,:,:)
  REAL, ALLOCATABLE :: var3d1(:,:,:), var3d2(:,:,:), var3d3(:,:,:)
  REAL, ALLOCATABLE :: var2d(:,:)

  INCLUDE   'mp.inc'

!-----------------------------------------------------------------

  INTERFACE 

    SUBROUTINE get_static_monthly(staticopt,ncid, vartype, jday,        &
                                  var2d, varin3d, istatus)
      INTEGER,          INTENT(IN)  :: staticopt
      INTEGER,          INTENT(IN)  :: ncid
      CHARACTER(LEN=1), INTENT(IN)  :: vartype
      INTEGER,          INTENT(IN)  :: jday
      REAL,             INTENT(OUT) :: var2d(:,:)
      REAL,             INTENT(IN)  :: varin3d(:,:,:)
      INTEGER,          INTENT(OUT) :: istatus
    END SUBROUTINE get_static_monthly

  END INTERFACE

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  nxlg = (nx - 1)*nproc_x + 1
  nylg = (ny - 1)*nproc_y + 1

!  write(0,*) 'nx = ',nx, 'ny = ',ny, 'nxlg = ',nxlg, 'nylg = ',nylg

  ALLOCATE(landusef(nx,ny,LanduseCategories), STAT = istatus)
  ALLOCATE(soiltyp_top(nx,ny,SoilCategories), STAT = istatus)
  ALLOCATE(soiltyp_bot(nx,ny,SoilCategories), STAT = istatus)

  ALLOCATE(var3d1(nxlg,nylg,LanduseCategories), STAT = istatus)
  ALLOCATE(var3d2(nxlg,nylg,SoilCategories),    STAT = istatus)
  ALLOCATE(var3d3(nxlg,nylg,SoilCategories),    STAT = istatus)
  ALLOCATE(var2d (nxlg,nylg),                   STAT = istatus)

  ALLOCATE(temdm(nxlg-1,nylg-1,MAX(LanduseCategories,SoilCategories)), STAT = istatus)
!
! Open the static NetCDF file which is created by geogrid.
!

  IF(myproc == 0) CALL open_static_file(filnam,nscid)

! we need to check whether the static file is
! on the same grid as WRF grid

  IF(myproc == 0) THEN 
    CALL check_static_grid(1,nscid,nx,ny,dx,dy,mapproj,                 &
                         trulat1,trulat2,trulon,istatus)
    IF(istatus /= 0) RETURN
  END IF

  IF (myproc == 0) THEN
    CALL get_ncd_3d(nscid,1,'LANDUSEF',nxlg-1,nylg-1,LanduseCategories,temdm,istatus)
    DO k = 1, LanduseCategories
      DO j = 1,nylg-1
        DO i = 1,nxlg-1
          var3d1(i,j,k) = temdm(i,j,k)
        END DO
      END DO
    END DO

    CALL get_ncd_3d(nscid,1,'SOILCTOP',nxlg-1,nylg-1,SoilCategories,temdm,istatus)
    DO k = 1, SoilCategories
      DO j = 1,nylg-1
        DO i = 1,nxlg-1
          var3d2(i,j,k) = temdm(i,j,k)
        END DO
      END DO
    END DO

    CALL get_ncd_3d(nscid,1,'SOILCBOT',nxlg-1,nylg-1,SoilCategories,temdm,istatus)
    DO k = 1, SoilCategories
      DO j = 1,nylg-1
        DO i = 1,nxlg-1
          var3d3(i,j,k) = temdm(i,j,k)
        END DO
      END DO
    END DO

    CALL edgfill(var3d1,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,LanduseCategories,1,LanduseCategories)
    CALL edgfill(var3d2,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,SoilCategories,1,SoilCategories)
    CALL edgfill(var3d3,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,SoilCategories,1,SoilCategories)
  END IF
  CALL wrf_split3d(var3d1,nx,ny,LanduseCategories,1,landusef)
  CALL wrf_split3d(var3d2,nx,ny,SoilCategories,1,soiltyp_top)
  CALL wrf_split3d(var3d3,nx,ny,SoilCategories,1,soiltyp_bot)

  IF (myproc == 0) THEN
    CALL get_ncd_3d(nscid,1,'GREENFRAC',nxlg-1,nylg-1,12,temdm,istatus)
    DO k = 1, 12
      DO j = 1,nylg-1
        DO i = 1,nxlg-1
          var3d1(i,j,k) = temdm(i,j,k)
        END DO
      END DO
    END DO
    CALL edgfill(var3d1,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,12,1,12)

    CALL get_static_monthly (1,nscid, 'g', jday, var2d, var3d1,istatus)

    DO j = 1,nylg-1
      DO i = 1,nxlg-1
        var3d2(i,j,1) = temdm(i,j,1)   ! find monthly minimum
        var3d2(i,j,2) = temdm(i,j,1)   ! find monthly maximum
        DO k = 2,12
          IF (var3d2(i,j,1) > temdm(i,j,k) ) THEN
            var3d2(i,j,1) = temdm(i,j,k)
          END IF
          IF (var3d2(i,j,2) < temdm(i,j,k) ) THEN
            var3d2(i,j,2) = temdm(i,j,k)
          END IF
        END DO
      END DO
    END DO
    CALL edgfill(var3d2,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,2,1,2)

  END IF
  IF (istatus /= 0) CALL mpexit(1)
  CALL wrf_split2d(var2d,nx,ny,1,vegfrct)
  CALL wrf_split2d(var3d2(:,:,1),nx,ny,1,shdmin)
  CALL wrf_split2d(var3d2(:,:,2),nx,ny,1,shdmax)
!-----------------------------------------------------------------------
!
! Green-ness values are expected in percent, but it is in fraction from geo_grid
!
!-----------------------------------------------------------------------

  vegfrct(:,:) = vegfrct(:,:) * 100.
  shdmin(:,:)  = shdmin (:,:) * 100.
  shdmax(:,:)  = shdmax (:,:) * 100.

  IF (myproc == 0) THEN
    CALL get_ncd_3d(nscid,1,'ALBEDO12M',nxlg-1,nylg-1,12,temdm,istatus)
    DO k = 1, 12
      DO j = 1,nylg-1
        DO i = 1,nxlg-1
          var3d2(i,j,k) = temdm(i,j,k)
        END DO
      END DO
    END DO
    CALL edgfill(var3d2,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,12,1,12)

    CALL get_static_monthly (1,nscid, 'a', jday, var2d, var3d2,istatus)
  END IF
  IF (istatus /= 0) CALL mpexit(1)
  CALL wrf_split2d(var2d,nx,ny,1,albbck)
  
  IF (myproc == 0) THEN
    CALL get_ncd_2d(nscid,1,'SNOALB',nxlg-1,nylg-1,temdm,istatus)
    DO j = 1,nylg-1
      DO i = 1,nxlg-1
        var2d(i,j) = temdm(i,j,1)
      END DO
    END DO
    CALL edgfill(var2d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,1,1,1)
  END IF
  CALL wrf_split2d(var2d,nx,ny,1,snoalb)
!-----------------------------------------------------------------------
!
! Albedo is expected in fraction, but it is percent from geo_grid
!
!-----------------------------------------------------------------------

  albbck(:,:) = albbck(:,:) / 100.
  snoalb(:,:) = snoalb(:,:) / 100.

  IF (myproc == 0) THEN
    CALL get_ncd_2d(nscid,1,'HGT_M',nxlg-1,nylg-1,temdm,istatus)
    DO j = 1,nylg-1
      DO i = 1,nxlg-1
        var2d(i,j) = temdm(i,j,1)
      END DO
    END DO
    CALL edgfill(var2d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,1,1,1)
  END IF
  CALL wrf_split2d(var2d,nx,ny,1,hgt)

  IF (myproc == 0) THEN
    CALL get_ncd_2d(nscid,1,'SOILTEMP',nxlg-1,nylg-1,temdm,istatus)
    DO j = 1,nylg-1
      DO i = 1,nxlg-1
        var2d(i,j) = temdm(i,j,1)
      END DO
    END DO
    CALL edgfill(var2d,1,nxlg,1,nxlg-1,1,nylg,1,nylg-1,1,1,1,1)
  END IF
  CALL wrf_split2d(var2d,nx,ny,1,tmn)

  IF(myproc == 0) CALL close_static_file(nscid)

  DEALLOCATE(temdm)
  DEALLOCATE(var3d1,var3d2,var3d3)
  DEALLOCATE(var2d)

  CALL land_soil_cat(nx,ny,LanduseCategories,SoilCategories,            &
                     landusef,soiltyp_top,soiltyp_bot,                  &
                     vegtyp,soiltyp,xland,xice,istatus)

  DEALLOCATE(landusef, soiltyp_top, soiltyp_bot)

  RETURN
END SUBROUTINE get_sfcdt_geogrid
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE land_soil_cat              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE land_soil_cat(nx,ny,land_cat,soil_cat,                       &
                     landusef,soiltyp_top,soiltyp_bot,                  &
                     vegtyp,soiltyp,xland,xice,istatus)

!-----------------------------------------------------------------------
!
! PURPOSE
!   Extract dominant landuse and soil texture index
!
!   See subroutine process_percent_cat_new in share/module_soil_pre.F
!
!-----------------------------------------------------------------------

  use wrf_metadata

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx, ny
  INTEGER, INTENT(IN) :: land_cat, soil_cat
  REAL,    INTENT(IN) :: landusef(nx,ny,land_cat)
  REAL,    INTENT(IN) :: soiltyp_top(nx,ny,soil_cat)
  REAL,    INTENT(IN) :: soiltyp_bot(nx,ny,soil_cat)
  INTEGER, INTENT(OUT) :: soiltyp(nx,ny)
  INTEGER, INTENT(OUT) :: vegtyp (nx,ny)
  REAL,    INTENT(OUT) :: xice   (nx,ny)
  REAL,    INTENT(OUT) :: xland  (nx,ny)
  INTEGER, INTENT(OUT) :: istatus

  INTEGER  :: i,j,k
  REAL     :: dominant_value
  INTEGER  :: dominant_index
  
!--------------------------------------------------------------------
!
! Vegetation index and soil texture index SoilCategories
!
!--------------------------------------------------------------------

  DO j = 1, ny
    DO i = 1,nx
      dominant_value = landusef(i,j,1)
      dominant_index = 1
      DO k = 2,LanduseCategories
        IF(k == ISWATER .AND. landusef(i,j,k) >= 0.5) THEN
          dominant_value = landusef(i,j,ISWATER)
          dominant_index = k
        ELSE IF(k /= ISWATER .AND. landusef(i,j,k) > dominant_value) THEN
          dominant_value = landusef(i,j,k)
          dominant_index = k
        END IF
      END DO
      IF(dominant_index == ISWATER) THEN
        xland(i,j) = 2.
      ELSE
        xland(i,j) = 1.
      END IF
      vegtyp(i,j) = dominant_index
    END DO
  END DO

  !  Compute the dominant SOIL TEXTURE INDEX, TOP.

  DO j = 1, ny
    DO i = 1,nx
      dominant_value = soiltyp_top(i,j,1)
      dominant_index = 1
      IF(xland(i,j) < 1.5) THEN
        DO k = 2,SoilCategories 
          IF(k /= ISWATER_SOIL .AND. soiltyp_top(i,j,k) > dominant_value) THEN
            dominant_value = soiltyp_top(i,j,k)
            dominant_index = k
          END IF
        END DO
        IF(dominant_value < 0.01) THEN
          dominant_index = 8
        END IF
      ELSE
        dominant_index = ISWATER_SOIL
      END IF
      soiltyp(i,j) = dominant_index
    END DO
  END DO

  xice = 0.
  WHERE(vegtyp == ISICE) xice = 1.

  istatus = 0

  RETURN
END SUBROUTINE land_soil_cat
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE couple                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE couple(nx,ny,nz,vname,mu,mub,var3d,msf,out3d,tem)

!-----------------------------------------------------------------
!
! PURPOSE
!   Couple the input variable with the air mass and the map scale
!   factor.
!
! NOTE: Refer to the same subroutine in WRF dyn_em/module_big_step_utilities_em.f.
!
!----------------------------------------------------------------
  IMPLICIT NONE
  INTEGER,      INTENT(IN)  :: nx,ny,nz
  CHARACTER(1), INTENT(IN)  :: vname
  REAL,         INTENT(IN)  :: mu(nx,ny), mub(nx,ny)
  REAL,         INTENT(IN)  :: var3d(nx,ny,nz)
  REAL,         INTENT(IN)  :: msf(nx,ny)
  REAL,         INTENT(OUT) :: out3d(nx,ny,nz)
  REAL,         INTENT(OUT) :: tem(nx,ny)

!------------------------------------------------------------------
!
!  Misc. local variables
!
!------------------------------------------------------------------
!
  INTEGER :: i,j,k

  INCLUDE   'mp.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  tem = 0.0

  IF(vname == 'U') THEN
    CALL wrf_mpsendrecv1de(mu, nx,ny,tem)  ! to make sure the fake zone, mu(nx,:)
    CALL wrf_mpsendrecv1de(mub,nx,ny,tem)  ! & mub(nx,:) contains valid data
    DO j = 1,ny
      DO i = 2,nx
        tem(i,j) = 0.5*(mu(i-1,j)+mub(i-1,j)+mu(i,j)+mub(i,j))
      END DO
    END DO
    tem(1,:) = mu(1,:) + mub(1,:)      ! The first row may be invalid
    !tem(nx,:) = mu(nx-1,:) + mub(nx-1,:)

    DO k = 1, nz
      DO j = 1, ny
        DO i = 1, nx
          out3d(i,j,k) = var3d(i,j,k)*tem(i,j)/msf(i,j)
        END DO
      END DO
    END DO

  ELSE IF(vname == 'V') THEN

    CALL wrf_mpsendrecv1dn(mu, nx,ny,tem)  ! to make sure the fake zone, mu(:,ny)
    CALL wrf_mpsendrecv1dn(mub,nx,ny,tem)  ! & mub(:,ny) contains valid data
    DO j = 2,ny
      DO i = 1,nx
        tem(i,j) = 0.5*(mu(i,j-1)+mub(i,j-1)+mu(i,j)+mub(i,j))
      END DO
    END DO
    tem(:,1)  = mu(:,1) + mub(:,1)     ! The firs column may be invalid
    !tem(:,ny) = mu(:,ny-1) + mub(:,ny-1)

    DO k = 1, nz
      DO j = 1, ny
        DO i = 1, nx
          out3d(i,j,k) = var3d(i,j,k)*tem(i,j)/msf(i,j)
        END DO
      END DO
    END DO

  ELSE IF(vname == 'W') THEN          ! No variable uses this "vname" so far

    DO k = 1, nz
      DO j = 1, ny
        DO i = 1, nx
          out3d(i,j,k) = var3d(i,j,k)*(mu(i,j)+mub(i,j))/msf(i,j)
        END DO
      END DO
    END DO

  ELSE      ! vname == 'H', or 'T', vertical does not matter here

    DO k = 1, nz
      DO j = 1, ny
        DO i = 1, nx
          out3d(i,j,k) = var3d(i,j,k)*(mu(i,j)+mub(i,j))
        END DO
      END DO
    END DO
  END IF
  
  RETURN
END SUBROUTINE couple
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  FUNCTION stuff_bdy                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE stuff_bdy(nx,ny,nz,bdy_width,stag,var3d,bs,bn,bw,be)

!-----------------------------------------------------------------
!
! Assign boundary arrays
!
! NOTE: Refer to the same subroutine in WRF share/module_bc.f.
!
!-----------------------------------------------------------------

  IMPLICIT NONE
 
  ! nx,ny and nz must be the non-staggered values

  INTEGER,          INTENT(IN) :: nx,ny,nz,bdy_width
  CHARACTER(LEN=1), INTENT(IN) :: stag
  REAL,             INTENT(IN) :: var3d(nx,ny,nz)       ! coupled variable

  REAL,    INTENT(OUT) :: bs(nx,nz,bdy_width)   ! south boundary
  REAL,    INTENT(OUT) :: bn(nx,nz,bdy_width)   ! north boundary
  REAL,    INTENT(OUT) :: bw(ny,nz,bdy_width)   ! west boundary
  REAL,    INTENT(OUT) :: be(ny,nz,bdy_width)   ! east boundary

!------------------------------------------------------------------
!
!  Misc. local variables
!
!------------------------------------------------------------------
!
  INTEGER :: i,j,k
  INTEGER :: ii,jj

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ! X start boundary  WEST
  DO k = 1,nz
    DO j = 1,ny
      DO i = 1, bdy_width
        bw(j,k,i) = var3d(i,j,k)
      END DO
    END DO
  END DO

  ! X end boundary    EAST
  IF(stag == 'U') THEN
    DO k = 1,nz
      DO j = 1,ny
        DO i = nx, nx-bdy_width+1,-1
          ii = nx - i + 1
          be(j,k,ii) = var3d(i,j,k)
        END DO
      END DO
    END DO
  ELSE
    DO k = 1,nz
      DO j = 1,ny
        DO i = nx-1, nx-bdy_width,-1
          ii = nx - i
          be(j,k,ii) = var3d(i,j,k)
        END DO
      END DO
    END DO
  END IF

  ! Y start boundary  SOUTH
  DO k = 1, nz
    DO j = 1, bdy_width
      DO i = 1, nx
        bs(i,k,j) = var3d(i,j,k)
      END DO
    END DO
  END DO

  ! Y end boundary    NORTH
  IF(stag == 'V') THEN
    DO k = 1, nz
      DO j = ny, ny-bdy_width+1,-1
        DO i = 1, nx
          jj = ny - j + 1
          bn(i,k,jj) = var3d(i,j,k)
        END DO
      END DO
    END DO
  ELSE
    DO k = 1, nz
      DO j = ny-1, ny-bdy_width,-1
        DO i = 1, nx
          jj = ny - j
          bn(i,k,jj) = var3d(i,j,k)
        END DO
      END DO
    END DO
  END IF

  RETURN
END SUBROUTINE stuff_bdy
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROTINE stuff_bdytend              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE stuff_bdytend(nx,ny,nz,bdy_width,stag,var3dold,var3dnew,     &
                         time_diff,bs,bn,bw,be)

!-----------------------------------------------------------------
!
! Compute and assign boundary tendency.
!
! NOTE: Refer to the same subroutine in WRF share/module_bc.f.
!
!-----------------------------------------------------------------

  IMPLICIT NONE
 
  ! nx,ny and nz must be the non-staggered values

  INTEGER, INTENT(IN)  :: nx,ny,nz,bdy_width
  CHARACTER(LEN=1), INTENT(IN) :: stag
  REAL,    INTENT(IN)  :: var3dold(nx,ny,nz)
  REAL,    INTENT(IN)  :: var3dnew(nx,ny,nz)
  REAL,    INTENT(IN)  :: time_diff
  REAL,    INTENT(OUT) :: bs(nx,nz,bdy_width)  ! south boundary tendency
  REAL,    INTENT(OUT) :: bn(nx,nz,bdy_width)  ! north
  REAL,    INTENT(OUT) :: bw(ny,nz,bdy_width)  ! west
  REAL,    INTENT(OUT) :: be(ny,nz,bdy_width)  ! east

!------------------------------------------------------------------
!
!  Misc. local variables
!
!------------------------------------------------------------------
!
  INTEGER :: i,j,k
  INTEGER :: ii,jj

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ! X start boundary  WEST
  DO k = 1,nz
    DO j = 1,ny
      DO i = 1, bdy_width
        bw(j,k,i) = (var3dnew(i,j,k)-var3dold(i,j,k))/time_diff
      END DO
    END DO
  END DO

  ! X end boundary    EAST
  IF(stag == 'U') THEN
    DO k = 1,nz
      DO j = 1,ny
        DO i = nx, nx-bdy_width+1,-1
          ii = nx - i + 1
          be(j,k,ii) = (var3dnew(i,j,k)-var3dold(i,j,k))/time_diff
        END DO
      END DO
    END DO
  ELSE
    DO k = 1,nz
      DO j = 1,ny
        DO i = nx-1, nx-bdy_width,-1
          ii = nx - i
          be(j,k,ii) = (var3dnew(i,j,k)-var3dold(i,j,k))/time_diff
        END DO
      END DO
    END DO
  END IF

  ! Y start boundary  SOUTH
  DO k = 1, nz
    DO j = 1, bdy_width
      DO i = 1, nx
        bs(i,k,j) = (var3dnew(i,j,k)-var3dold(i,j,k))/time_diff
      END DO
    END DO
  END DO

  ! Y end boundary    NORTH
  IF(stag == 'V') THEN
    DO k = 1, nz
      DO j = ny, ny-bdy_width+1,-1
        DO i = 1, nx
          jj = ny - j + 1
          bn(i,k,jj) = (var3dnew(i,j,k)-var3dold(i,j,k))/time_diff
        END DO
      END DO
    END DO
  ELSE
    DO k = 1, nz
      DO j = ny-1, ny-bdy_width,-1
        DO i = 1, nx
          jj = ny - j
          bn(i,k,jj) = (var3dnew(i,j,k)-var3dold(i,j,k))/time_diff
        END DO
      END DO
    END DO
  END IF

  RETURN
END SUBROUTINE stuff_bdytend
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE stuff_bdy2d                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE stuff_bdy2d(nx,ny,nz,bdy_width,var2d,bs,bn,bw,be)

!------------------------------------------------------------------
!
!  Assign boundary arrays for 2D array.
!
!-----------------------------------------------------------------

  IMPLICIT NONE
 
  ! nx,ny and nz must be the non-staggered values

  INTEGER, INTENT(IN)  :: nx,ny,nz,bdy_width
  REAL,    INTENT(IN)  :: var2d(nx,ny)
  REAL,    INTENT(OUT) :: bs(nx,nz,bdy_width)  ! south boundary array
  REAL,    INTENT(OUT) :: bn(nx,nz,bdy_width)  ! north
  REAL,    INTENT(OUT) :: bw(ny,nz,bdy_width)  ! west
  REAL,    INTENT(OUT) :: be(ny,nz,bdy_width)  ! east

!------------------------------------------------------------------
!
!  Misc. local variables
!
!------------------------------------------------------------------
!
  INTEGER :: i,j,k
  INTEGER :: ii,jj

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ! X start boundary  WEST
  DO k = 1,nz
    DO j = 1,ny
      DO i = 1, bdy_width
        bw(j,k,i) = var2d(i,j)
      END DO
    END DO
  END DO

  ! X end boundary    EAST
  DO k = 1,nz
    DO j = 1,ny
      DO i = nx-1, nx-bdy_width,-1
        ii = nx - i
        be(j,k,ii) = var2d(i,j)
      END DO
    END DO
  END DO

  ! Y start boundary  SOUTH
  DO k = 1, nz
    DO j = 1, bdy_width
      DO i = 1, nx
        bs(i,k,j) = var2d(i,j)
      END DO
    END DO
  END DO

  ! Y end boundary    NORTH
  DO k = 1, nz
    DO j = ny-1, ny-bdy_width,-1
      DO i = 1, nx
        jj = ny - j
        bn(i,k,jj) = var2d(i,j)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE stuff_bdy2d
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE stuff_bdytend2d            ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE stuff_bdytend2d(nx,ny,nz,bdy_width,var2dold,var2dnew,   &
                         time_diff,bs,bn,bw,be)

  IMPLICIT NONE
 
  ! nx,ny and nz must be the non-staggered values

  INTEGER, INTENT(IN)  :: nx,ny,nz,bdy_width
  REAL,    INTENT(IN)  :: var2dold(nx,ny)
  REAL,    INTENT(IN)  :: var2dnew(nx,ny)
  REAL,    INTENT(IN)  :: time_diff
  REAL,    INTENT(OUT) :: bs(nx,nz,bdy_width)
  REAL,    INTENT(OUT) :: bn(nx,nz,bdy_width)
  REAL,    INTENT(OUT) :: bw(ny,nz,bdy_width)
  REAL,    INTENT(OUT) :: be(ny,nz,bdy_width)

!------------------------------------------------------------------
!
!  Misc. local variables
!
!------------------------------------------------------------------
!
  INTEGER :: i,j,k
  INTEGER :: ii,jj

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ! X start boundary  WEST
  DO k = 1,nz
    DO j = 1,ny
      DO i = 1, bdy_width
        bw(j,k,i) = (var2dnew(i,j)-var2dold(i,j))/time_diff
      END DO
    END DO
  END DO

  ! X end boundary    EAST
  DO k = 1,nz
    DO j = 1,ny
      DO i = nx-1, nx-bdy_width,-1
        ii = nx - i
        be(j,k,ii) = (var2dnew(i,j)-var2dold(i,j))/time_diff
      END DO
    END DO
  END DO

  ! Y start boundary  SOUTH
  DO k = 1, nz
    DO j = 1, bdy_width
      DO i = 1, nx
        bs(i,k,j) = (var2dnew(i,j)-var2dold(i,j))/time_diff
      END DO
    END DO
  END DO

  ! Y end boundary    NORTH
  DO k = 1, nz
    DO j = ny-1, ny-bdy_width,-1
      DO i = 1, nx
        jj = ny - j
        bn(i,k,jj) = (var2dnew(i,j)-var2dold(i,j))/time_diff
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE stuff_bdytend2d
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE rotate_UV                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE rotate_UV(nx_wrf,ny_wrf,nz_wrf,mapproj,sclfct,latnot,trulon, &
                     swx_arps,swy_arps,mapproj_wrf,                     &
                     sclfct_wrf,lattru_wrf,lontru_wrf,swx_wrf,swy_wrf,  &
                     lonu_wrf,lonv_wrf,u_wrf,v_wrf,uatv_wrf,vatu_wrf,   &
                     tem1,tem2,tem3,tem4,istatus)
!
!-----------------------------------------------------------------------
!
! PURPOSE:
!
!   Get uatv_wrf & vatu_wrf for doing horizontal interpolation.
!   Do not support mpi.
!
!-----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: nx_wrf, ny_wrf, nz_wrf
  INTEGER, INTENT(IN)  :: mapproj, mapproj_wrf
  REAL,    INTENT(IN)  :: sclfct,  sclfct_wrf
  REAL,    INTENT(IN)  :: latnot(2),lattru_wrf(2)
  REAL,    INTENT(IN)  :: trulon,  lontru_wrf
  REAL,    INTENT(IN)  :: swx_arps,  swx_wrf
  REAL,    INTENT(IN)  :: swy_arps,  swy_wrf
  REAL,    INTENT(IN)  :: lonu_wrf(nx_wrf,ny_wrf)
  REAL,    INTENT(IN)  :: lonv_wrf(nx_wrf,ny_wrf)
  REAL,    INTENT(INOUT) :: u_wrf(nx_wrf,ny_wrf,nz_wrf)
  REAL,    INTENT(INOUT) :: v_wrf(nx_wrf,ny_wrf,nz_wrf)
  REAL,    INTENT(OUT)   :: uatv_wrf(nx_wrf,ny_wrf,nz_wrf)
  REAL,    INTENT(OUT)   :: vatu_wrf(nx_wrf,ny_wrf,nz_wrf)
  REAL,    INTENT(OUT)   :: tem1(nx_wrf,ny_wrf,nz_wrf)
  REAL,    INTENT(OUT)   :: tem2(nx_wrf,ny_wrf,nz_wrf)
  REAL,    INTENT(OUT)   :: tem3(nx_wrf,ny_wrf,nz_wrf)
  REAL,    INTENT(OUT)   :: tem4(nx_wrf,ny_wrf,nz_wrf)
  INTEGER, INTENT(OUT)   :: istatus

!-----------------------------------------------------------------------
!
! Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER  :: i, j, k

  INTEGER  :: iproj
  REAL     :: scl, trlat(2), trlon,x0,y0

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code ... ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  Save previous map projection values
!
!-----------------------------------------------------------------------
!
  CALL getmapr(iproj,scl,trlat,trlon,x0,y0)
!
!-----------------------------------------------------------------------
!
!  Get uatv and vatu
!
!-----------------------------------------------------------------------
!
  DO k=1,nz_wrf

    ! get u at v grid point locations
    DO j=2,ny_wrf
      DO i=1,nx_wrf-1
        uatv_wrf(i,j,k) = 0.25*(u_wrf(i,j-1,k)+u_wrf(i+1,j-1,k)  &
                                + u_wrf(i,j,k)+u_wrf(i+1,j,k))
      END DO
    END DO
    DO i=1,nx_wrf-1
      uatv_wrf(i,1,k) = 0.5*(u_wrf(i,1,k)+u_wrf(i+1,1,k))
    END DO
    DO j=2,ny_wrf
      uatv_wrf(nx_wrf,j,k) = 0.5*(u_wrf(nx_wrf,j-1,k)+u_wrf(nx_wrf,j,k))
    END DO
    uatv_wrf(nx_wrf,1,k) = u_wrf(nx_wrf,1,k)

    ! get v at u grid point locations
    DO j=1,ny_wrf-1
      DO i=2,nx_wrf
        vatu_wrf(i,j,k) = 0.25*(v_wrf(i-1,j,k)+v_wrf(i,j,k)  &
                                +v_wrf(i-1,j+1,k)+v_wrf(i,j+1,k))
      END DO
    END DO
    DO j=1,ny_wrf-1
      vatu_wrf(1,j,k) = 0.5*(v_wrf(1,j,k)+v_wrf(1,j+1,k))
    END DO
    DO i=2,nx_wrf
      vatu_wrf(i,ny_wrf,k) = 0.5*(v_wrf(i-1,ny_wrf,k)+v_wrf(i,ny_wrf,k))
    END DO
    vatu_wrf(1,ny_wrf,k) = v_wrf(1,ny_wrf,k)

  END DO

!
!-----------------------------------------------------------------------
!
!  Establish ARPS map projection
!
!-----------------------------------------------------------------------
!
  CALL setmapr(mapproj,sclfct,latnot,trulon)
  CALL setorig(1,swx_arps,swy_arps)
!
!-----------------------------------------------------------------------
!
!  Rotate from ARPS grid to earth-relative
!
!-----------------------------------------------------------------------
!
  DO k = 1,nz_wrf
    CALL uvmptoe(nx_wrf,ny_wrf,u_wrf(:,:,k),vatu_wrf(:,:,k),        &
                 lonu_wrf,tem1(:,:,k),tem2(:,:,k))

    CALL uvmptoe(nx_wrf,ny_wrf,uatv_wrf(:,:,k),v_wrf(:,:,k),        &
                 lonv_wrf,tem3(:,:,k),tem4(:,:,k))
  END DO
  ! tem1 = u,    tem2 = vatu
  ! tem3 = uatv, tem4 = v
!
!-----------------------------------------------------------------------
!
!  Establish WRF map projection
!
!-----------------------------------------------------------------------
!
  CALL setmapr(mapproj_wrf,sclfct_wrf,lattru_wrf,lontru_wrf)
  CALL setorig( 1, swx_wrf, swy_wrf)

!
!-----------------------------------------------------------------------
!
!  Rotate from earth to WRF grid-relative
!
!-----------------------------------------------------------------------
!
  DO k = 1,nz_wrf
    CALL uvetomp(nx_wrf,ny_wrf,tem1(:,:,k),tem2(:,:,k),        &
                 lonu_wrf,u_wrf(:,:,k),vatu_wrf(:,:,k))

    CALL uvetomp(nx_wrf,ny_wrf,tem3(:,:,k),tem4(:,:,k),        &
                 lonv_wrf,uatv_wrf(:,:,k),v_wrf(:,:,k))
  END DO
!
!-----------------------------------------------------------------------
!
!  Reset map projection to previous values
!
!-----------------------------------------------------------------------
!
  CALL setmapr(iproj,scl,trlat,trlon)
  CALL setorig(1,x0,y0)

  istatus = 0

  RETURN
END SUBROUTINE rotate_UV


