PROGRAM mrgtrn
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Program to merge terrain files to insure continuity for
!  external grids near boundaries while accepting smaller-scale
!  terrain features in the interior.
!
!  Steps:  Read in 2 ARPS terrain files
!             1) large-scale external grid
!             2) fine-scale terrain on final grid
!          Interpolate large scale terrain to fine scale grid
!          Correct for any biases in the large-scale terrain
!          Blend the two grids, using the large-scale data near
!                the domain edge.
!  Compilation/linking:
!
!   f77 -o mergetrn mergetrn.o maproj3d.o \
!          extlib3d.o outlib3d.o iolib3d.o ibmlib3d.o
!
!  AUTHOR: Keith Brewster
!  01/03/95
!
!  MODIFICATION HISTORY:
!  10/09/96 (K. Brewster)
!  Corrected input of map parameters by adding routine rdtrnall.
!
!  06/10/97 (K. Brewster)
!  Added option to leave a border of external terrain data before
!  beginning the merging.  The border is specific in the input file
!  using iborder and jborder.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'

  INTEGER :: nx, ny
  INTEGER :: nx_ext, ny_ext
!
  REAL, allocatable :: xs(:),ys(:)
  REAL, allocatable :: xgr(:,:),ygr(:,:)
  REAL, allocatable :: lat(:,:),lon(:,:)
  REAL, allocatable :: hterain(:,:)
  REAL, allocatable :: hterain1(:,:)
  REAL, allocatable :: hterain2(:,:)
  REAL, allocatable :: hterout(:,:)
  REAL, allocatable :: work(:,:)
  INTEGER, allocatable :: iloc(:,:),jloc(:,:)
!
  REAL, allocatable :: xs_ext(:),ys_ext(:)
  REAL, allocatable :: terext(:,:)

  REAL, allocatable :: dxfld(:)
  REAL, allocatable :: dyfld(:)
  REAL, allocatable :: rdxfld(:)
  REAL, allocatable :: rdyfld(:)
  REAL, allocatable :: slopey(:,:)
  REAL, allocatable :: alphay(:,:)
  REAL, allocatable :: betay(:,:)
  INTEGER, allocatable :: iw(:,:)
!
!-----------------------------------------------------------------------
!
!  Terrain namelist variables
!
!-----------------------------------------------------------------------
!
  REAL :: dx_ext,dy_ext

  CHARACTER (LEN=256) :: terndta_ext
  NAMELIST /terrain/ nx,ny,dx,dy,nx_ext,ny_ext,dx_ext,dy_ext,           &
                     terndta,terndta_ext

  INTEGER :: iborder,jborder
  REAL :: rlen
  INTEGER :: rmvbias,nsmth
  NAMELIST /ternmrg/ iborder,jborder,rlen,rmvbias,nsmth

  INTEGER :: ovrmap,mapcol,mapgrid
  REAL :: latgrid,longrid
  CHARACTER (LEN=256) :: mapfile
  NAMELIST /map_plot/ ovrmap,mapcol,mapgrid,latgrid,longrid,            &
            mapfile

!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: ternfn
  REAL :: latnot(2)
  INTEGER :: i,j,lenstr,ksmth,iproj_ext,korder
  INTEGER :: istart,iend,jstart,jend
  REAL :: scale_ext,trulat1_ext,trulat2_ext
  REAL :: trulon_ext,ctrlat_ext,ctrlon_ext
  REAL :: xctr,yctr,x0_ext,y0_ext
  REAL :: x0,y0,sum,sum_ext,bias,wt
  REAL :: w1,w2,w3,w4,scalex,scaley

  INTEGER :: nxpic,nypic
  PARAMETER (nxpic=1,nypic=1)
  REAL :: sm
  PARAMETER (sm=0.5)
  INTEGER :: iorder
  PARAMETER (iorder=1)

  INTEGER :: ncl
  REAL :: cl(100)

  INTEGER :: iplt,mode
  REAL :: angl,xlimit,ylimit,xymax
  INTEGER :: istatus

  CHARACTER(LEN=256) :: outfilename
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  iborder=0
  jborder=0
  rlen=10000.
  rmvbias=0
  nsmth=0

  ovrmap=1
  mapcol=1
  mapgrid=0
  latgrid=10.
  longrid=10.

  WRITE(6,'(/a/)')                                                      &
      ' Please input control parameters for terrain merging.'

  READ(5,terrain,END=600)
  READ(5,ternmrg,END=600)
  READ(5,map_plot,END=600)

  allocate(xs(nx),ys(ny),stat=istatus)
  allocate(xgr(nx,ny),ygr(nx,ny),stat=istatus)
  allocate(lat(nx,ny),lon(nx,ny),stat=istatus)
  allocate(hterain(nx,ny),stat=istatus)
  allocate(hterain1(nx,ny),stat=istatus)
  allocate(hterain2(nx,ny),stat=istatus)
  allocate(hterout(nx,ny),stat=istatus)
  allocate(work(nx,ny),stat=istatus)
  allocate(iloc(nx,ny),jloc(nx,ny),stat=istatus)
!
  allocate(xs_ext(nx_ext),ys_ext(ny_ext),stat=istatus)
  allocate(terext(nx_ext,ny_ext),stat=istatus)

  allocate(dxfld(nx_ext),stat=istatus)
  allocate(dyfld(ny_ext),stat=istatus)
  allocate(rdxfld(nx_ext),stat=istatus)
  allocate(rdyfld(ny_ext),stat=istatus)
  allocate(slopey(nx_ext,ny_ext),stat=istatus)
  allocate(alphay(nx_ext,ny_ext),stat=istatus)
  allocate(betay(nx_ext,ny_ext),stat=istatus)
  allocate(iw(nx,ny),stat=istatus)
!
!-----------------------------------------------------------------------
!
!  Read in external terrain data
!
!-----------------------------------------------------------------------
!
  lenstr =  LEN(terndta_ext)
  CALL strlnth( terndta_ext, lenstr)
  mapproj=-1
!
  CALL rdtrnall( nx_ext,ny_ext,dx_ext,dy_ext,                           &
                 iproj_ext,trulat1_ext,trulat2_ext,                     &
                 trulon_ext,scale_ext,                                  &
                 ctrlat_ext,ctrlon_ext,                                 &
                 terndta_ext(1:lenstr), terext)

  IF( scale_ext == 0.0 ) scale_ext = 1.0
!
  IF(ctrlon_ext > 180.) ctrlon_ext=ctrlon_ext-360.
  IF(trulon_ext > 180.) trulon_ext=trulon_ext-360.
!
!-----------------------------------------------------------------------
!
!  Read in fine-scale terrain data
!
!-----------------------------------------------------------------------
!
  lenstr =  LEN(terndta)
  CALL strlnth( terndta, lenstr)
  mapproj=-1
!
  CALL rdtrnall( nx,ny,dx,dy,                                           &
                 mapproj,trulat1,trulat2,trulon,sclfct,                 &
                 ctrlat,ctrlon,                                         &
                 terndta(1:lenstr), hterain )

  IF( sclfct == 0.0 ) sclfct = 1.0
!
  IF(ctrlon > 180.) ctrlon=ctrlon-360.
  IF(trulon > 180.) trulon=trulon-360.
  latnot(1)=trulat1
  latnot(2)=trulat2
  CALL setmapr(mapproj,sclfct,latnot,trulon)
!
  CALL lltoxy(1,1,ctrlat,ctrlon,xctr,yctr)
  x0=xctr - 0.5*(nx-3)*dx
  y0=yctr - 0.5*(ny-3)*dy
!
  DO i=1,nx
    xs(i)=x0+(i-1)*dx
  END DO
  DO j=1,ny
    ys(j)=y0+(j-1)*dy
  END DO
!
  CALL xytoll(nx,ny,xs,ys,lat,lon)
!
!-----------------------------------------------------------------------
!
!  Interpolate the external data to the fine-scale grid
!
!-----------------------------------------------------------------------
!
  latnot(1)=trulat1_ext
  latnot(2)=trulat2_ext
  CALL setmapr(iproj_ext,scale_ext,latnot,trulon_ext)
  CALL lltoxy(1,1,ctrlat_ext,ctrlon_ext,xctr,yctr)

  x0_ext=xctr - 0.5*(nx_ext-3)*dx_ext
  y0_ext=yctr - 0.5*(ny_ext-3)*dy_ext
!
!-----------------------------------------------------------------------
!
!  Set up external grid
!
!-----------------------------------------------------------------------
!
  DO i=1,nx_ext
    xs_ext(i)=x0_ext+(i-1)*dx_ext
  END DO
  DO j=1,ny_ext
    ys_ext(j)=y0_ext+(j-1)*dy_ext
  END DO
!
!-----------------------------------------------------------------------
!
!  Find x,y locations of ARPS scalar grid in terms of the
!  external grid.
!
!-----------------------------------------------------------------------
!
  CALL lltoxy(nx,ny,lat,lon,xgr,ygr)
!
!-----------------------------------------------------------------------
!
!  Find i,j indices in the external data of the fine grid points
!
!-----------------------------------------------------------------------
!
  CALL setijloc(nx,ny,nx_ext,ny_ext,xgr,ygr,xs_ext,ys_ext,              &
                iloc,jloc)
  CALL setdxdy(nx_ext,ny_ext,                                           &
               1,nx_ext,1,ny_ext,                                       &
               xs_ext,ys_ext,dxfld,dyfld,rdxfld,rdyfld)
  korder=MIN(iorder,2)
  CALL fldint2d(nx,ny,nx_ext,ny_ext,                                    &
                1,nx,1,ny,                                              &
                1,nx_ext,1,ny_ext,                                      &
                korder,xgr,ygr,terext,xs_ext,ys_ext,iloc,jloc,          &
                dxfld,dyfld,rdxfld,rdyfld,                              &
                slopey,alphay,betay,                                    &
                hterain1)
!
!-----------------------------------------------------------------------
!
!  Correct for any bias in the external data
!  This is to avoid artifical hills in the terrain due to
!  any bias in the external data
!
!-----------------------------------------------------------------------
!
  IF(rmvbias > 0) THEN
    sum=0.
    sum_ext=0.
    DO j=1,ny-1
      DO i=1,nx-1
        sum=sum+hterain(i,j)
        sum_ext=sum_ext+hterain1(i,j)
      END DO
    END DO
    bias=(sum_ext-sum)/FLOAT((nx-1)*(ny-1))
    WRITE(6,'(/a,f6.2,a/)') ' Removing bias of ',bias,' meters'
!
    DO j=1,ny-1
      DO i=1,nx-1
        hterain1(i,j)=hterain1(i,j)-bias
      END DO
    END DO
  ELSE
    WRITE(6,'(/a,i4,a/)')                                               &
        ' rmvbias= ',rmvbias,' --skipping bias removal.'
  END IF
!
!-----------------------------------------------------------------------
!
!  Blend the two terrain arrays
!
!-----------------------------------------------------------------------
!
!  scalex=-dx*dx/(rlen*rlen)
!  scaley=-dy*dy/(rlen*rlen)
  scalex=dx/rlen
  scaley=dy/rlen
  WRITE(6,'(/a,f10.2)')                                                 &
      ' Merging terrain fields using length scale',rlen
  WRITE(6,'(a,f10.2,a)')                                                &
      ' Which is equal to ',(1./scalex),' grid lengths in x'
  WRITE(6,'(a,f10.2,a/)')                                               &
      '   and is equal to ',(1./scaley),' grid lengths in y'
  istart=iborder+1
  iend=(nx-1)-iborder
  jstart=jborder+1
  jend=(ny-1)-jborder
!
!  Start with all points set to external terrain height
!
  DO j=1,ny-1
    DO i=1,nx-1
      hterain2(i,j) = hterain1(i,j)
    END DO
  END DO
!
  DO j=jstart,jend
    DO i=istart,iend
!    w1=exp((i-1)*(i-1)*scalex)
!    w2=exp((nx-1-i)*(nx-1-i)*scalex)
!    w3=exp((j-1)*(j-1)*scaley)
!    w4=exp((ny-1-j)*(ny-1-j)*scaley)
!    wt=amin1(1.0,amax1(0.0,w1,w2,w3,w4))
!    hterain2(i,j) = (1.0-wt)*hterain(i,j) + wt*hterain1(i,j)
      w1=(i-istart)*scalex
      w2=(iend-i)*scalex
      w3=(j-jstart)*scaley
      w4=(jend-j)*scaley
      wt=AMAX1(0.0,AMIN1(1.0,w1,w2,w3,w4))
      hterain2(i,j) = wt*hterain(i,j) + (1.-wt)*hterain1(i,j)
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Call edgehammer to filter noise near the edges caused
!  by merging
!
!-----------------------------------------------------------------------
!
  IF(nsmth > 0) THEN
    ksmth=MAX((nsmth/4),4)
    WRITE(6,'(/a)') ' Calling edge hammer'
    WRITE(6,'(a,i4,a,i4,a//)')                                          &
        ' nsmth= ',nsmth,' affecting ',ksmth,' points from boundary'
    CALL edgeham(nx,ny,hterain2,work,hterout,nsmth,sm)
  ELSE
    WRITE(6,'(a)') ' nsmth= ',nsmth,' --skipping edge hammer.'
    DO j=1,ny-1
      DO i=1,nx-1
        hterout(i,j) = hterain2(i,j)
      END DO
    END DO
  END IF
!
!-----------------------------------------------------------------------
!
!  Write the final terrain field
!
!-----------------------------------------------------------------------
!
  ternfn = terndta(1:lenstr)//'.new'
  lenstr = lenstr+4

  terndmp = 1

  CALL fnversn(ternfn, lenstr)
  CALL writtrn(nx,ny,ternfn(1:lenstr),dx,dy,                            &
               mapproj,trulat1,trulat2,trulon,sclfct,                   &
               ctrlat,ctrlon,hterout)
!
!-----------------------------------------------------------------------
!
!  Whip up some plots
!
!-----------------------------------------------------------------------
!
  CALL xdevic_new(1,outfilename,0,0)
  CALL xdspac(0.9)
  CALL xafstyl(1)
  CALL xartyp(2)

  angl = 0.0
  CALL xspace(nxpic, nypic, angl , xlimit,ylimit)

  CALL xaxfmt( '(i10)' )
  CALL xpmagn(0.10*xlimit, 0.10*ylimit)
!
!-----------------------------------------------------------------------
!
!  Set map
!
!-----------------------------------------------------------------------
!

  xl = (nx-3)*dx * 0.001
  yl = (ny-3)*dy * 0.001
  xorig = 0.0
  yorig = 0.0

  CALL xstpjgrd(mapproj,trulat1,trulat2,trulon,                         &
                ctrlat,ctrlon,xl,yl,xorig,yorig)

  DO j=1,ny
    DO i=1,nx
      xgr(i,j)=(i-1)*dx * 0.001
      ygr(i,j)=(j-1)*dy * 0.001
    END DO
  END DO

  DO iplt=1,3

    CALL xnwpic

    xl = (nx-3)*dx * 0.001
    yl = (ny-3)*dy * 0.001

    xymax=AMAX1(xl,yl)

    CALL xmap(0.0,xymax, 0.0, xymax)
    CALL xaxsca(0.0,xl, (dx*0.001), 0.0, yl, (dy*0.001) )

    cl(1)=0.0
    cl(2)=50.
    CALL xnctrs(10,30)
    mode=1

    IF(iplt == 1) THEN
      CALL xconta(hterain(2,2),xgr(2,2),ygr(2,2),iw,nx,nx-3,ny-3,       &
          cl, ncl,mode )
    ELSE IF (iplt == 2) THEN
      CALL xconta(hterain2(2,2),xgr(2,2),ygr(2,2),iw,nx,nx-3,ny-3,      &
          cl, ncl,mode )
    ELSE IF (iplt == 3) THEN
      CALL xconta(hterout(2,2),xgr(2,2),ygr(2,2),iw,nx,nx-3,ny-3,       &
          cl, ncl,mode )
    END IF
    IF(ovrmap == 1) THEN
      CALL xdrawmap(10,mapfile,mapgrid,latgrid,longrid)
    END IF

  END DO

  CALL xgrend

  STOP
!
  600 CONTINUE
  PRINT *, ' Error reading input data, stopping'
  STOP
!
END PROGRAM mrgtrn
!

SUBROUTINE edgeham(nx,ny,hterain,work,newtern,nsmth,sm)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  01/03/95
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: nx,ny
  REAL :: hterain(nx,ny)
  REAL :: work(nx,ny)
  REAL :: newtern(nx,ny)
  INTEGER :: nsmth
  REAL :: sm
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,ksm,kwid
  REAL :: wcen,wsid
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO j=1,ny-1
    DO i=1,nx-1
      newtern(i,j)=hterain(i,j)
    END DO
  END DO
!
  wcen=1.-sm
  wsid=0.5*sm
  DO ksm=nsmth,1,-1
    kwid=ksm/4
    kwid=MAX(kwid,4)
    DO j=1,ny-1
      DO i=1,nx-1
        work(i,j)=newtern(i,j)
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Smooth in the x direction near west boundary
!
!-----------------------------------------------------------------------
!
    DO j=1,ny-1
      DO i=2,kwid
        work(i,j)=wsid*(newtern(i-1,j)+newtern(i+1,j)) +                &
                  wcen* newtern(i,j)
      END DO
    END DO
    DO j=1,ny-1
      work(1,j)=wsid*(newtern(1,j)+newtern(2,j)) +                      &
                wcen* newtern(1,j)
    END DO
!
!-----------------------------------------------------------------------
!
!  Smooth in the x direction near east boundary
!
!-----------------------------------------------------------------------
!
    DO j=1,ny-1
      DO i=(nx-kwid),nx-2
        work(i,j)=wsid*(newtern(i-1,j)+newtern(i+1,j)) +                &
                  wcen* newtern(i,j)
      END DO
    END DO
    DO j=1,ny-1
      work(nx-1,j)=wsid*(newtern(nx-1,j)+newtern(nx-2,j)) +             &
                   wcen* newtern(nx-1,j)
    END DO
!
!-----------------------------------------------------------------------
!
!  Smooth in the y direction near south boundary
!
!-----------------------------------------------------------------------
!
    DO j=1,ny-1
      DO i=1,nx-1
        newtern(i,j)=work(i,j)
      END DO
    END DO
    DO j=2,kwid
      DO i=1,nx-1
        newtern(i,j)=wsid*(work(i,j-1)+work(i,j+1)) +                   &
                     wcen* work(i,j)
      END DO
    END DO
    DO i=1,nx-1
      newtern(i,1)=wsid*(work(i,1)+work(i,2)) +                         &
                   wcen* work(i,1)
    END DO
!
!-----------------------------------------------------------------------
!
!  Smooth in the y direction near north boundary
!
!-----------------------------------------------------------------------
!
    DO j=(ny-kwid),ny-2
      DO i=1,nx-1
        newtern(i,j)=wsid*(work(i,j-1)+work(i,j+1)) +                   &
                     wcen* work(i,j)
      END DO
    END DO
    DO i=1,nx-1
      newtern(i,ny-1)=wsid*(work(i,ny-1)+work(i,ny-2)) +                &
                      wcen* work(i,ny-1)
    END DO
  END DO
  RETURN
END SUBROUTINE edgeham
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE RDTRNALL                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE rdtrnall(nx,ny, dx,dy,                                       &
           mapproj,trulat1,trulat2,trulon,sclfct,                       &
           ctrlat,ctrlon,                                               &
           terndta, hterain )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read the terrain data into model array hterain from a specified
!  terrain data file.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  2/27/1994.
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!
!    dx       Grid interval in x-direction
!    dy       Grid interval in y-direction
!
!    terndta     Terrain data file name
!
!  OUTPUT:
!
!    hterain  Terrain height (m)
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

  INTEGER :: nx            ! Number of grid points in the x-direction
  INTEGER :: ny            ! Number of grid points in the y-direction
  REAL :: dx               ! Grid interval in x-direction
  REAL :: dy               ! Grid interval in y-direction
  INTEGER :: mapproj
  REAL :: trulat1
  REAL :: trulat2
  REAL :: trulon
  REAL :: sclfct
  REAL :: ctrlat
  REAL :: ctrlon
  CHARACTER (LEN=*) :: terndta ! Terrain data file name

  REAL :: hterain(nx,ny)   ! Terrain height.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: inunit,istat
  INTEGER :: nxin,nyin,idummy,ierr
  REAL :: dxin,dyin,rdummy,amin,amax

!MP insert:      character*80 savename            !Message passing code.
!MP insert:      integer lenfl                    !Message passing code.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!
!-----------------------------------------------------------------------
!
!  Define a uniform model grid.
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Read in the terrain data.
!
!-----------------------------------------------------------------------
!
  CALL getunit( inunit )

!MP insert:        lenfl = 80                     !Message passing code.
!MP insert:        CALL strlnth( terndta, lenfl ) !Message passing code.
!MP insert:        savename(1:80) = terndta(1:80) !Message passing code.
!MP insert:        write(terndta, '(a,a,2i2.2)')  !Message passing code.
!MP insert:     :          savename(1:lenfl),'_',loc_x,loc_y

  CALL asnctl ('NEWLOCAL', 1, ierr)
  CALL asnfile(terndta, '-F f77 -N ieee', ierr)

  OPEN(UNIT=inunit,FILE=terndta,FORM='unformatted',STATUS='old',        &
       IOSTAT=istat)

!MP insert:        terndta(1:80) = savename(1:80) !Message passing code.

  IF( istat /= 0) THEN

    WRITE(6,'(/1x,a,a,/1x,a/)')                                         &
        'Error occured when opening terrain data file ',                &
        terndta,' Job stopped in READTRN.'
    STOP

  END IF

  READ(inunit,ERR=999) nxin,nyin


  IF((nx /= nxin).OR.(ny /= nyin))THEN

    WRITE(6,'(a,/a,i5,a,i5,/a,i5,a,i5)')                                &
        ' Array size in the terrain data does not match that of the',   &
        ' model grid. Dimensions in the data were nx=',nxin,            &
        ', ny=',nyin,' the model grid size were nx=',nx,' ny= ', ny
    WRITE(6,'(a)') ' Job stopped in subroutine READTRN.'
    STOP

  END IF

  READ(inunit,ERR=999) idummy,mapproj,idummy,idummy,idummy,             &
               idummy,idummy,idummy,idummy,idummy,                      &
               idummy,idummy,idummy,idummy,idummy,                      &
               idummy,idummy,idummy,idummy,idummy

  READ(inunit,ERR=999) dxin  ,dyin  ,ctrlat,ctrlon,rdummy,              &
               rdummy,trulat1,trulat2,trulon,sclfct,                    &
               rdummy,rdummy,rdummy,rdummy,rdummy,                      &
               rdummy,rdummy,rdummy,rdummy,rdummy

  IF(ABS((dx-dxin)/dx) > 0.01.OR.ABS((dy-dyin)/dy) > 0.01)THEN

    WRITE(6,'(a,a,/2(a,f13.4),/2(a,f13.4))')                            &
        'Grid intervals in the terrain data do not match those ',       &
        'in the model.','In the data  dx=',dxin,', dy=',dyin,           &
        'In the model dx=',dx,' dy= ', dy
    WRITE(6,'(a)') ' Job stopped in subroutine READTRN.'


    STOP

  END IF

  READ(inunit,ERR=999) hterain

  WRITE(6,'(1x,a/)') 'Minimum and maximum terrain height:'

  CALL a3dmax0(hterain,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1,amax,amin)
  WRITE(6,'(1x,2(a,e13.6))') 'htrnmin = ', amin,', htrnmax=',amax


  RETURN

  999   WRITE(6,'(1x,a)')                                               &
        'Error in reading terrain data. Job stopped in READTRN.'
  STOP

END SUBROUTINE rdtrnall
