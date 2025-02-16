!
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE ASSIMOUT                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE assimout(nx,ny,nz,retfname,                                  &
           x,y,z,zp,                                                    &
           u,v,                                                         &
           ptprt,pprt,qvprt,qr,                                         &
           ubar,vbar,ptbar,pbar,qvbar,                                  &
           tem1,tem2,tem3,tem4,tem5)
!
!--------------------------------------------------------------------------
!
!  PURPOSE:
!
!  This subroutine outputs the adjusted wind, thermo and moisture
!  fields to ARPS data analysis system (ADAS) as "pseudo-soundings".
!
!---------------------------------------------------------------------------
!
!  AUTHOR: Limin Zhao
!  5/06/96.
!
!  MODIFICATION HISTORY:
!
!  06/10/96 (Keith Brewster and Limin Zhao)
!  Found a bug in calculating the location of "pseudo-soundings".
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE     ! Force explicit declarations
!
!-----------------------------------------------------------------------
!
!   Arrays to be read in:
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz

  REAL :: x(nx)
  REAL :: y(ny)
  REAL :: z(nz)
  REAL :: zp(nx,ny,nz)

  REAL :: u(nx,ny,nz)
  REAL :: v(nx,ny,nz)

  REAL :: ptprt(nx,ny,nz)
  REAL :: pprt(nx,ny,nz)
  REAL :: qvprt(nx,ny,nz)
  REAL :: qr(nx,ny,nz)

  REAL :: ubar(nx,ny,nz)
  REAL :: vbar(nx,ny,nz)
  REAL :: ptbar(nx,ny,nz)
  REAL :: pbar(nx,ny,nz)
  REAL :: qvbar(nx,ny,nz)

  REAL :: tem1(nx,ny,nz)
  REAL :: tem2(nx,ny,nz)
  REAL :: tem3(nx,ny,nz)
  REAL :: tem4(nx,ny,nz)
  REAL :: tem5(nx,ny,nz)
!
!-----------------------------------------------------------------------
!

  INTEGER :: ncl
  PARAMETER (ncl=500)
! NOTE: ncl must be greater than or equal to nx,ny, & nz
!
!-----------------------------------------------------------------------
!
  REAL :: tem1d1(ncl),tem1d2(ncl),tem1d3(ncl),tem1d4(ncl)
  REAL :: tem1d5(ncl),tem1d6(ncl),tem1d7(ncl),tem1d8(ncl)
  REAL :: tem1d9(ncl)
!

  REAL :: xsc(ncl),ysc(ncl)

!
!-----------------------------------------------------------------------
!
!    "Fake" radar id stuff
!
!-----------------------------------------------------------------------
!
  INTEGER :: iretfmt
  CHARACTER (LEN=256) :: retfname
  PARAMETER(iretfmt=1)
!  character*4 radid
!  real latrad,lonrad,elvrad
!
!
!  parameter(radid='KTLX',          !need change for different radar
!    :          latrad= 35.3331,
!    :          lonrad=-97.2778,
!cc  :          elvrad=389.4)
!    :          elvrad=384.0)
!
!   parameter(radid='KDDC',
!  :          latrad= 37.7597,
!  :          lonrad=-99.9678,
!  :          elvrad=814.0)
!
!   parameter(radid='KICT',
!  :          latrad= 37.654444,
!  :          lonrad=-97.442778,
!  :          elvrad=427.0)
!
!
!-----------------------------------------------------------------------
!
!    Misc internal variables.
!
!-----------------------------------------------------------------------
!
  REAL :: latnot(2)

  INTEGER :: i,j,k,ireturn

  CHARACTER (LEN=256) :: filename
  CHARACTER (LEN=256) :: grdbasfn
  INTEGER :: ngchan,nchan1,lenfil,lengbf
  INTEGER :: nchin,iyr

  REAL :: xctr,yctr,x0,y0
  REAL :: flag
  REAL :: temv,rhobar,refdbz,svnfrth,arg
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'assim.inc'
  INCLUDE 'grid.inc'
!
!-----------------------------------------------------------------------
!
!  Routines called:
!
!-----------------------------------------------------------------------
!
  EXTERNAL wtretcol
!
!
!=======================================================================
!
!    Special mods for AMS99 TD retr from SDVR
!
!  open(unit=59,file='radloc.input',form='formatted',status='old')
!  read (59,83) radid
!  read (59,*) latrad
!  read (59,*) lonrad
!  read (59,*) elvrad
!  close(unit=59)
!83   format(a4)
!
!=======================================================================
!
!
!-----------------------------------------------------------------------
!
!    Check the dimension for local 1-D arrays. It must be
!    defined larger than the model's dimension.
!
!-----------------------------------------------------------------------
!
  spval = 999.0
  WRITE (*,*) "XXX IASSIM nx,ny,nz,ncl",nx,ny,nz,ncl

  IF(ncl < nx .OR. ncl < ny .OR. ncl < nz) THEN
    WRITE(6,'(3(/1x,a))')                                               &
            'Dimension of ncl in ASSIMOUT is less than ARPS model,',    &
            'please redefine it. it should be larger than OR',          &
            'equal TO arps model dimensions.'
    WRITE(6,'(1x,a)')                                                   &
             'ncl,nx,ny,nz:', ncl,nx,ny,nz
    WRITE(6,'(1x,a)')                                                   &
             'Program aborted in ASSIMOUT'
    STOP
  END IF
!
!-----------------------------------------------------------------------
!
!    Put the (u,v) at the scalar points.
!
!-----------------------------------------------------------------------
!

  DO i=1,nx-1
    xsc(i) = 0.5*(x(i)+x(i+1))
  END DO
  xsc(nx)=2.0*xsc(nx-1)-xsc(nx-2)

  DO j=1,ny-1
    ysc(j) = 0.5*(y(j)+y(j+1))
  END DO
  ysc(ny) = 2.0*ysc(ny-1)-ysc(ny-2)

  DO i=1,nx
    x(i) = xsc(i)
  END DO

  DO j=1,ny
    y(j) = ysc(j)
  END DO

  svnfrth = 7./4.
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem4(i,j,k) = 0.5*(zp(i,j,k) + zp(i,j,k+1))
        tem1(i,j,k)=0.5*(u(  i,j,k)+u(  i+1,j,k))
        tem2(i,j,k)=0.5*(v(i,  j,k)+v(i,  j+1,k))
!
!
        temv = ptbar(i,j,k) * (pbar(i,j,k)/p0)**rddcp
        temv = temv * (rddrv+qvbar(i,j,k)) /                            &
                        ( rddrv*(1.0+qvbar(i,j,k)) )
                             ! Base state virtual temperature
        rhobar = pbar(i,j,k) / (rd*temv)
        arg = 17300.0*( 1000.0*rhobar                                   &
                 * MAX(0.0,qr(i,j,k)) )**svnfrth  !! needs a pcntg
        refdbz = 10.0 * ALOG10( MAX(arg,1.0) )

!
!
        flag=tem5(i,j,k)
        IF ( flag == spval ) THEN
          ptprt(i,j,k) = -999.0
          pprt(i,j,k)  = -999.0
          qr(i,j,k)    = -999.0
          qvprt(i,j,k) = -999.0
          ptbar(i,j,k) = -999.0
          pbar(i,j,k)  = -999.0
          qvbar(i,j,k) = -999.0
          tem1(i,j,k)  = -999.0
          tem2(i,j,k)  = -999.0
          tem3(i,j,k)  = -999.0
        ELSE
          tem3(i,j,k) = 1.0
        END IF
!
! test #1: no u & v output for re-analysis; i.e., the original ones are used.
!       Others are output every other point.
!
!      tem1(i,j,k)  = -999.0  ! no output for these variables
!      tem2(i,j,k)  = -999.0
!      qr(i,j,k)    = -999.0
!      qvprt(i,j,k) = -999.0
!      qvbar(i,j,k) = -999.0
!
!      if( refdbz .lt. 20.0 ) then
!        ptprt(i,j,k) = -999.0
!        pprt (i,j,k) = -999.0
!        ptbar(i,j,k) = -999.0
!        pbar (i,j,k) = -999.0
!
!        tem3 (i,j,k) = -999.0
!      else
!        tem3 (i,j,k) = 1.0
!      endif
!
! test #2: u & v output for re-analysis; output every other point.
!
!
!      if( i.eq.(i/2)*2 .or. j.eq.(j/2)*2 ) then
!        tem1(i,j,k)  = -999.0
!        tem2(i,j,k)  = -999.0
!        qr(i,j,k)    = -999.0
!        qvprt(i,j,k) = -999.0
!        ptprt(i,j,k) = -999.0
!        pprt(i,j,k)  = -999.0
!        ptbar(i,j,k) = -999.0
!        pbar(i,j,k)  = -999.0
!        qvbar(i,j,k) = -999.0
!
!        tem3(i,j,k)  = -999.0
!      else
!        tem3(i,j,k) = 1.0
!      endif
!
!
      END DO
    END DO
  END DO

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        zp(i,j,k) = tem4(i,j,k)
      END DO
    END DO
  END DO


  latnot(1) = trulat1
  latnot(2) = trulat2

  CALL setmapr(mapproj,1.0,latnot,trulon)
  CALL lltoxy(1,1,ctrlat,ctrlon,xctr,yctr)
  x0=xctr-(((nx-3)/2)*(x(2)-x(1)))
  y0=yctr-(((ny-3)/2)*(y(2)-y(1)))
  CALL setorig(1,x0,y0)


  iyr=MOD(year,100)
  WRITE(retfname,'(a,a,i4,2i2.2,a,2i2.2)')                              &
       radid,'ret.',year,month,day,'.',hour,minute
  PRINT *, ' Writing data into ',retfname


  CALL wtretcol(nx,ny,nz,                                               &
              2,nx-2,2,ny-2,2,nz-2,                                     &
              iyr,month,day,hour,minute,second,                         &
              iretfmt,retfname,radid,latrad,lonrad,elvrad,              &
              x,y,zp,                                                   &
              tem1,tem2,ptprt,pprt,qvprt,qr,                            &
              ptbar,pbar,qvbar,tem3,                                    &
              tem1d1,tem1d2,tem1d3,tem1d4,                              &
              tem1d5,tem1d6,tem1d7,tem1d8,tem1d9)


  RETURN
END SUBROUTINE assimout

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WTRETCOL                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######    university of Oklahoma.  All rights reserved.     ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE wtretcol(nx,ny,nz,                                           &
           ibeg,iend,jbeg,jend,kbeg,kend,                               &
           iyr,imon,iday,ihr,imin,isec,                                 &
           iretfmt,retfname,radid,latrad,lonrad,elvrad,                 &
           xsc,ysc,zpsc,                                                &
           us,vs,ptprt,pprt,qvprt,qr,                                   &
           ptbar,pbar,qvbar,retrflg,                                    &
           outk,outhgt,outu,outv,                                       &
           outpr,outpt,outqv,outqr,outret)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Test writing of retrieval columns.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  Writes gridded radar data to a file
!
!  INPUT
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (km)
!    zp       z coordinate of grid points in computational space (m)
!
!    uprt     x component of perturbation velocity (m/s)
!    vprt     y component of perturbation velocity (m/s)
!    wprt     Vertical component of perturbation velocity in Cartesian
!             coordinates (m/s).
!
!    ptprt    Perturbation potential temperature (K)
!    pprt     Perturbation pressure (Pascal)
!
!    qvprt    Perturbation water vapor mixing ratio (kg/kg)
!    qc       Cloud water mixing ratio (kg/kg)
!    qr       Rainwater mixing ratio (kg/kg)
!    qi       Cloud ice mixing ratio (kg/kg)
!    qs       Snow mixing ratio (kg/kg)
!    qh       Hail mixing ratio (kg/kg)
!
!    ubar     Base state x velocity component (m/s)
!    vbar     Base state y velocity component (m/s)
!    wbar     Base state z velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    rhobar   Base state air density (kg/m**3)
!    qvbar    Base state water vapor mixing ratio (kg/kg)
!    qr       Rainwater mixing ratio (kg/kg)
!
!-----------------------------------------------------------------------
!
!
  IMPLICIT NONE
!
  INTEGER :: nx,ny,nz
  INTEGER :: ibeg,iend,jbeg,jend,kbeg,kend
!
  REAL :: xsc(nx)
  REAL :: ysc(ny)
  REAL :: zpsc(nx,ny,nz)
  REAL :: us(nx,ny,nz)         ! total u velocity component at scalar points
  REAL :: vs(nx,ny,nz)         ! total v velocity component at scalar points
  REAL :: ptprt (nx,ny,nz)     ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz)     ! Perturbation pressure (Pascal)
  REAL :: qvprt (nx,ny,nz)     ! Perturbation water vapor specific
                               ! humidity (kg/kg)
  REAL :: qr    (nx,ny,nz)     ! Rain water mixing ratio (kg/kg)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal)
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific humidity
                               ! (kg/kg)
  REAL :: retrflg(nx,ny,nz)
!
  INTEGER :: iretfmt
  CHARACTER (LEN=256) :: retfname
  CHARACTER (LEN=4)   :: radid
  REAL :: latrad
  REAL :: lonrad
  REAL :: elvrad
  INTEGER :: iyr,imon,iday,ihr,imin,isec
!
!-----------------------------------------------------------------------
!
!  Retrieval output variables
!
!-----------------------------------------------------------------------
!
  REAL :: outk(nz)
  REAL :: outhgt(nz)
  REAL :: outu(nz)
  REAL :: outv(nz)
  REAL :: outpr(nz)
  REAL :: outpt(nz)
  REAL :: outqv(nz)
  REAL :: outqr(nz)
  REAL :: outret(nz)
!
!-----------------------------------------------------------------------
!
!  Retrieved data threshold
!
!-----------------------------------------------------------------------
!
  REAL :: retrthr
  PARAMETER(retrthr=0.0)
!
!-----------------------------------------------------------------------
!
!  Include file
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iunit,myr,itime
  INTEGER :: i,j,k,klev,kk,kntcol
  INTEGER :: ireftim,idummy
  INTEGER :: ierr
  REAL :: gridlat,gridlon,elev,rdummy
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  WRITE(6,*) '===> into WTRETCOL'
  ireftim = 0
  myr=1900+iyr
  IF(myr < 1960) myr=myr+100
  CALL ctim2abss(myr,imon,iday,ihr,imin,isec,itime)

  CALL asnctl ('NEWLOCAL', 1, ierr)
  CALL asnfile(retfname, '-F f77 -N ieee', ierr)
  CALL getunit(iunit)
!
!  Open file for output
!
  OPEN(iunit,FILE=retfname,STATUS='unknown',                            &
       FORM='unformatted')
!
!  Write retrieval description variables
!
  idummy = 0
  WRITE(iunit) radid
  WRITE(iunit) ireftim,itime,idummy,idummy,idummy,                      &
               idummy,idummy,idummy,idummy,idummy
!
!  Write grid description variables
!  This should provide enough info to verify that the
!  proper grid has been chosen.  To recreate the grid,
!  icluding elevation information,
!  the reading program should get a grid-base-file
!  named runname.grdbasfil
!

  idummy=0
  rdummy=0.
  WRITE(iunit) runname
  WRITE(iunit) hdmpfmt,strhopt,mapproj,nx,ny,                           &
               nz,idummy,idummy,idummy,idummy
  WRITE(iunit) dx,dy,dz,dzmin,ctrlat,                                   &
               ctrlon,trulat1,trulat2,trulon,sclfct,                    &
               latrad,lonrad,elvrad,rdummy,rdummy

!cc
!cc

  WRITE(6,*) runname
  WRITE(6,*) hdmpfmt,strhopt,mapproj,idummy,idummy,                     &
              idummy,idummy,idummy,idummy,idummy
  WRITE(6,*) dx,dy,dz,dzmin,ctrlat,                                     &
      ctrlon,trulat1,trulat2,trulon,sclfct,                             &
      latrad,lonrad,elvrad,rdummy,rdummy


!
!  For each horizontal grid point form a column of remapped
!  data containing the non-missing grid points
!
  kntcol=0
  DO j=jbeg,jend
    DO i=ibeg,iend
      klev=0
      DO k=kbeg,kend
        IF(retrflg(i,j,k) > retrthr) THEN
          klev=klev+1
          outk(klev)=FLOAT(k)
          outhgt(klev)=zpsc(i,j,k)
          outu(klev)=us(i,j,k)
          outv(klev)=vs(i,j,k)
          outpr(klev)=pprt(i,j,k)+pbar(i,j,k)
          outpt(klev)=ptprt(i,j,k)+ptbar(i,j,k)
          outqv(klev)=qvprt(i,j,k)+qvbar(i,j,k)
          outqr(klev)=qr(i,j,k)
          outret(klev)=retrflg(i,j,k)
!------------------------------------------------------
!  print some diagnostics
!
          WRITE(6,*) i,j,k,outu(klev),outpt(klev),outqr(klev),outqv(klev), &
              outret(klev)
!
!------------------------------------------------------
        END IF
      END DO
!
!  If there are data in this column, write them to the file.
!
      IF(klev > 0) THEN
        kntcol=kntcol+1
        CALL xytoll(1,1,xsc(i),ysc(j),gridlat,gridlon)
        elev=0.5*(zpsc(i,j,1)+zpsc(i,j,2))
        WRITE(iunit) i,j,xsc(i),ysc(j),                                 &
                     gridlat,gridlon,elev,klev
        WRITE(iunit) (outk(kk),kk=1,klev)
        WRITE(iunit) (outhgt(kk),kk=1,klev)
        WRITE(iunit) (outu(kk),kk=1,klev)
        WRITE(iunit) (outv(kk),kk=1,klev)
        WRITE(iunit) (outpr(kk),kk=1,klev)
        WRITE(iunit) (outpt(kk),kk=1,klev)
        WRITE(iunit) (outqv(kk),kk=1,klev)
        WRITE(iunit) (outqr(kk),kk=1,klev)
        WRITE(iunit) (outret(kk),kk=1,klev)

!      write(99,*) i,j,xsc(i),ysc(j),gridlat,gridlon,elev,klev

!       write(93,*)kntcol,xsc(i),ysc(j),gridlat,gridlon
!       write(93,*)(kk,outu(kk),outv(kk),outqv(kk),kk=1,klev)
      END IF
      IF(i == 92 .AND. j == 92) THEN
        CALL xytoll(1,1,xsc(i),ysc(j),gridlat,gridlon)
      END IF
    END DO
  END DO

  CLOSE(iunit)
  CALL retunit(iunit)
!
!  Report on what data were written
!
  WRITE(6,'(//a,i2.2,i2.2,i2.2,a1,i2.2,a1,i2.2)')                       &
                    ' Output statistics for time ',                     &
                      iyr,imon,iday,' ',ihr,':',imin
  WRITE(6,'(a,i6,a,/a,i6,a//)')                                         &
           ' There were ',kntcol,' columns written ',                   &
           ' of a total ',(nx*ny),' possible.'

  RETURN
END SUBROUTINE wtretcol
